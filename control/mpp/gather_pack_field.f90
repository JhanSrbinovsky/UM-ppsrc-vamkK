! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Gathers a field from many processors to one processor and packs

! Subroutine Interface:
  SUBROUTINE gather_pack_field(                                   &
                          local_field,global_field,               &
                          local_row_len,local_rows,               &
                          global_row_len,global_rows,             &
                          grid_type,halo_type,                    &
                          global_gather_pe,proc_group,            &
                          packing, im_ident, packing_type,  &
                          num_out,                                &
                          comp_accrcy, rmdi)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE UM_ParVars
USE unite_output_files_mod, ONLY : &
    unite_coex_files
USE gcom_mod
USE Submodel_Mod
IMPLICIT NONE


! Description:
!  Takes a model field that has been decomposed over a group of
!  processors, and gathers the data together so that one processor
!  contains the entire global field. Optionally WGDOS packs the
!  gathered data as it is sent to get parallelism of the packing
!  process.

! Method:
!  In the most simple situation,
!  a send and receive map is constructed which instructs the GCOM
!  permute operation to do a gather from all processors in the
!  group to the GATHER_PE.

!  In the more ususal version, there is a 2 stage gather. The
!  first gathers a "block_factor" number of rows together. These
!  are packed (optionally) and then a final gather gets the full
!  field into place. This will require some adjustment to ensure
!  the packed data is correct.


! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Subroutine Arguments:

  INTEGER, INTENT(IN) ::                                          &
  local_row_len                                                   &
                     ! IN length of rows in local part of field
, local_rows                                                      &
                     ! IN number of rows in local part of field
, global_row_len                                                  &
                     ! IN length of rows in global field
, global_rows                                                     &
                     ! IN number of rows in global field
, grid_type                                                       &
                     ! IN type (P,U or V) of grid
, halo_type                                                       &
                     ! IN halo type (hence width) of grid
, global_gather_pe                                                &
                     ! IN processor to gather global field to
, proc_group         ! IN group ID of processors involved here


! Optional Arguments to handle the COEX packing if necessary

LOGICAL, INTENT(IN) ::                                            &
   packing           ! IN: Set .true. if packing of the input
                     !     field is to be packed!

INTEGER, INTENT(IN) ::                                            &
  im_ident           ! IN: Internal model identifier

INTEGER, INTENT(INOUT) ::                                         &
  packing_type       ! IN/OUT: This flag is zero on input,
                     !         then stash packing is selected,
                     !         and the routine returns the
                     !         packing flag.
                     !
                     !         If the variable is set to 2 on
                     !         input then 32-bit packing for
                     !         dumpfiles is selected

INTEGER, INTENT(OUT) ::                                           &
  num_out            ! OUT: Number of 32-bit IBM words in the
                     !      Packed field for WDGOS packing

INTEGER, INTENT(IN) ::                                            &
  comp_accrcy        ! IN: Packing Accuracy in Power of 2

REAL, INTENT(IN) ::                                               &
  rmdi               ! IN: Missing data indicator

! Remaining Non-Optional Arguments

REAL, INTENT(IN) ::                                               &
  local_field(local_row_len*local_rows)
                     ! IN local part of field

REAL, INTENT(OUT) ::                                              &
  global_field(global_row_len*global_rows)
                     ! OUT (on PE GATHER_PE) global field

! Local variables

INTEGER                                                           &
   send_map(7,1)                                                  &
,  n_mess_to_recv                                                 &
,  n_mess_to_send                                                 &
,  send_map_2(7,1)                                                &
,  n_mess_to_recv_2                                               &
,  n_mess_to_send_2
INTEGER, POINTER :: receive_map(:,:)=>NULL()
INTEGER, POINTER :: receive_map_2(:,:)=>NULL()

INTEGER                                                           &
  old_global_row_len                                              &
                        ! value on last call
, old_global_rows                                                 &
                        ! value on last call
, old_proc_group                                                  &
                        ! value on last call
, old_gather_pe                                                   &
                        ! value on last call
, old_decomp                                                      &
                        ! value on last call
, old_grid_type                                                   &
                        ! value on last call
, old_halo_type         ! value on last call

  SAVE send_map,n_mess_to_send,receive_map,n_mess_to_recv,        &
     old_global_row_len,old_global_rows,old_proc_group,           &
     old_gather_pe,old_decomp,                                    &
     old_grid_type,old_halo_type
DATA old_global_row_len,old_global_rows,old_proc_group,           &
     old_gather_pe,old_decomp,                                    &
     old_grid_type,old_halo_type                                  &
   / -1234, -1234, -1234, -1234, -1234, -1234, -1234/

INTEGER ::                                                        &
  iproc                                                           &
, jproc                                                           &
, kproc                                                           &
, lproc                                                           &
, info                                                            &
, flag                                                            &
, gather_pe                                                       &
                       ! Local gather PE for a group/block of
                       ! rows
, row_start_pe                                                    &
                       ! First PE in a block - may or may not
                       ! the gather PE
, data_address                                                    &
, data_size                                                       &
, block_pe                                                        &
                       ! The PE holding the current block of
                       ! rows on the GLOBAL_GATHER_PE
, pes_per_block                                                   &
                       ! Number of PE's in a block (nproc_x*
                       ! block_factor)
, n_local_rows                                                    &
                       ! Number of rows per block
, packed_buffer_size                                              &
                       ! Size of the buffer to hold packed data
, length_fullwrd                                                  &
                       ! Size of full word on this machine 64-bit
                       ! for Cray
, local_packing_type   ! Local copy of the packing_type, which
                       ! maybe not present

INTEGER :: icode       ! error code

INTEGER, ALLOCATABLE :: icomp(:)
                       ! Local array to hold compressed data

! The block_factor is the number of rows that are gathered together
! for packing purposes before the final gather. This is best set to
! 1 where maximum parallelism can be exploited (eg on the NEC SX-6).
! Massively parallel systems may benefit from a higher value.
INTEGER, PARAMETER   :: block_factor = 1

REAL  :: my_global_buffer(global_row_len*global_rows+2)

CHARACTER (LEN=*), PARAMETER :: routinename='GATHER_PACK_FIELD'
CHARACTER (LEN=80)           :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-------------------------------------------------------

! Compute the local gather PE - head of row PE
IF (lhook) CALL dr_hook('GATHER_PACK_FIELD',zhook_in,zhook_handle)

IF (.NOT.ASSOCIATED(receive_map)) &
    ALLOCATE(receive_map(7,nproc_max))
IF (.NOT.ASSOCIATED(receive_map_2)) &
    ALLOCATE(receive_map_2(7,nproc_max))


pes_per_block=block_factor*nproc_x
gather_pe=(mype/pes_per_block)*pes_per_block
row_start_pe=gather_pe

! Check if the global_gather_pe is in the block - if so
! make sure it is the gather_pe
IF(global_gather_pe >= gather_pe .AND.                            &
   global_gather_pe <= MIN(nproc, gather_pe+pes_per_block)-1) THEN
  gather_pe=global_gather_pe
END IF

! 0.0 Can we use the same send/receive map that we calculated
!     last time round?

IF ((global_row_len  /=  old_global_row_len) .OR.                 &
    (global_rows     /=  old_global_rows   ) .OR.                 &
    (proc_group      /=  old_proc_group    ) .OR.                 &
    (gather_pe       /=  old_gather_pe     ) .OR.                 &
    (grid_type       /=  old_grid_type     ) .OR.                 &
    (halo_type       /=  old_halo_type     ) .OR.                 &
    (current_decomp_type  /=  old_decomp  )) THEN

!       Different arguments from the last call so we need
!       to calculate a new send/receive map
  IF (grid_type  ==  fld_type_unknown) THEN
    WRITE(6,*) 'GATHER_PACK_FIELD : Bad field type'
    WRITE(6,*) 'Field will not be scattered.'
    cmessage='GATHER_PACK_FIELD : Bad field type'
    icode=1


    CALL ereport(routinename,icode,cmessage)
  END IF


! 2.0 Set up send map

  send_map(s_destination_pe,1) = gather_pe

  send_map(s_base_address_in_send_array,1) =                      &
    halosize(2,halo_type)*local_row_len+                          &
    1+halosize(1,halo_type)

  send_map(s_number_of_elements_in_item,1)=blsize(2,grid_type)

  send_map(s_stride_in_send_array,1) = local_row_len

  send_map(s_element_length,1) =                                  &
    local_row_len-2*halosize(1,halo_type)

  send_map(s_base_address_in_recv_array,1) =                      &
    datastart_f(1,grid_type) +                                    &
   (datastart_f(2,grid_type)-1)*global_row_len

  send_map(s_stride_in_recv_array,1) = global_row_len

  n_mess_to_send=1

! 3.0 Set up the receive map (for PE GATHER_PE only)

  n_mess_to_recv=0

  IF (mype  ==  gather_pe) THEN

! Loop over PE's in this block
    DO jproc=row_start_pe,                                        &
     MIN(nproc, row_start_pe+pes_per_block)-1
      iproc=jproc-row_start_pe

      receive_map(r_source_pe,iproc+1) = jproc

      receive_map(r_base_address_in_recv_array,iproc+1) =         &
          g_datastart_f(1,grid_type,jproc)+                       &
          (g_datastart_f(2,grid_type,jproc)-1)                    &
           *global_row_len

      receive_map(r_number_of_elements_in_item,iproc+1) =         &
          g_blsize(2,grid_type,jproc)

      receive_map(r_stride_in_recv_array,iproc+1) =               &
        global_row_len

      receive_map(r_element_length,iproc+1) =                     &
          g_blsize(1,grid_type,jproc)

      receive_map(r_base_address_in_send_array,iproc+1) =         &
        halosize(2,halo_type)*                                    &
          g_lasize(1,grid_type,halo_type,jproc)+                  &
        halosize(1,halo_type)+1

      receive_map(r_stride_in_send_array,iproc+1) =               &
          g_lasize(1,grid_type,halo_type,jproc)

      n_mess_to_recv=n_mess_to_recv+1

    END DO
  END IF

  old_global_row_len=global_row_len
  old_global_rows=global_rows
  old_proc_group=proc_group
  old_gather_pe=gather_pe
  old_decomp=current_decomp_type
  old_grid_type=grid_type
  old_halo_type=halo_type

END IF  ! we need to recalculate send/receive maps.

! 4.0 Do the exchange of data

flag=0  ! This is currently ignored at GCG v1.1
info=gc_none

local_packing_type = 0

! Only the gather PE's need to do anything now

!      Run length encoding applies only to unpacked ocean fields, but
!      most ocean fields remain unpacked even when packing profiles
!      5 or 6 are set. Hence selecting, for example, both packing
!      profile 5 and run length encoding makes sense.

IF(packing .AND. global_rows  >=  2) THEN
  ! Climate wgdos packing has been selected for current
  !    file stream via UMUI
  IF(comp_accrcy  >   -99) THEN
     ! 2. STASH packing profile for the field is set.
     packing_type = 1
  END IF
END IF

num_out=0
local_packing_type = packing_type

! If there is no packing, then use the output buffer
! if I am the gather PE, or if it has been reduced to a single
! stage gather.  Otherwise, use a temporary buffer.
IF(local_packing_type /= 1 .AND.                                  &
 ((gather_pe == global_gather_pe) .OR.                            &
  (block_factor >= nproc_y))) THEN

  CALL gcg_ralltoalle(local_field,send_map,n_mess_to_send,        &
                      local_row_len*local_rows,                   &
                      global_field,receive_map,n_mess_to_recv,    &
                      global_row_len*global_rows,                 &
                      proc_group,flag,icode)
ELSE

  CALL gcg_ralltoalle(local_field,send_map,                       &
                      n_mess_to_send,                             &
                      local_row_len*local_rows,                   &
                      my_global_buffer,receive_map,               &
                      n_mess_to_recv,                             &
                      global_row_len*global_rows,                 &
                      proc_group,flag,icode)
END IF

IF (icode /= 0) THEN
  cmessage = ' Error in GCG_RALLTOALLE'


  CALL ereport(routinename,icode,cmessage)
END IF


! Now we must check if the packing_type is 1
IF(local_packing_type == 1) THEN

! Only gather_pe's need to do anything
  IF(mype == gather_pe) THEN

! Work out how much local data we have
    n_local_rows=0
    DO kproc=1, block_factor
      IF(row_start_pe+(kproc-1)*nproc_x <  nproc) THEN
        n_local_rows=n_local_rows+                                &
         g_blsize(2,grid_type,row_start_pe+(kproc-1)*nproc_x)
      END IF
    END DO

! Setup a buffer for the packed data
    packed_buffer_size=n_local_rows*global_row_len+2
    ALLOCATE (icomp(packed_buffer_size))

! Pack the data
    icode=0
    length_fullwrd=64

! Check if the2 stage gather has been eliminated - if so
! put the data straight into the output array
    IF(block_factor >= nproc_y) THEN
! DEPENDS ON: coex
      CALL coex(                                                  &
         my_global_buffer(                                        &
                   g_datastart_f(1,grid_type,row_start_pe)+       &
                   (g_datastart_f(2,grid_type,row_start_pe)-1)*   &
                          global_row_len),                        &
                          global_row_len*global_rows,             &
                          global_field, packed_buffer_size,       &
                          global_row_len, n_local_rows,           &
                          num_out,                                &
                          comp_accrcy, .TRUE., rmdi,              &
                          length_fullwrd,                         &
                          icode, cmessage)

! Still doing 2 stage gather
    ELSE

! DEPENDS ON: coex
      CALL coex(                                                  &
         my_global_buffer(g_datastart_f(1,grid_type,row_start_pe)+&
                   (g_datastart_f(2,grid_type,row_start_pe)-1)*   &
                          global_row_len),                        &
                          global_row_len*global_rows,             &
                          icomp, packed_buffer_size,              &
                          global_row_len, n_local_rows,           &
                          num_out,                                &
                          comp_accrcy, .TRUE., rmdi,              &
                          length_fullwrd,                         &
                          icode, cmessage)
    END IF

    IF (icode /= 0) THEN
      cmessage='Error in COEX'


      CALL ereport(routinename, icode, cmessage)
    END IF


  END IF ! mype == gather_pe

END IF !packing_type is 1

IF(local_packing_type == 1) THEN

  IF(mype == gather_pe) THEN

    IF(block_factor <  nproc_y) THEN

      CALL gc_rsend(1001+mype, (num_out+1)/2, global_gather_pe,   &
       info,                                                      &
       my_global_buffer(                                          &
       g_datastart_f(1,grid_type,row_start_pe)+                   &
      (g_datastart_f(2,grid_type,row_start_pe)-1)*global_row_len),&
       icomp)

    END IF ! are we doing 2 stage gather?

  END IF ! mype == gather_pe

END IF !  packing_type is 1

IF(local_packing_type == 1) THEN

! Now pack the into the output buffer on the global_gather_pe
  IF(mype == global_gather_pe) THEN

! Check if we are doing 2 stage gather
    IF(block_factor <  nproc_y) THEN

! Loop over each processors Contribution
      DO jproc=0, nproc_y-1, block_factor
        block_pe=jproc*nproc_x

! Preserve the row_start_pe for this block of data
        iproc=block_pe

! If the global_gather_pe is in the current block, then the block_pe
! needs to be set to global_gather_pe, not the row/block leader

        IF(row_start_pe == block_pe) THEN
          block_pe=global_gather_pe
        END IF

        n_local_rows=0
        DO kproc=1, block_factor
          IF(iproc+(kproc-1)*nproc_x <  nproc) THEN
            n_local_rows=n_local_rows+                            &
             g_blsize(2,grid_type,iproc+(kproc-1)*nproc_x)
          END IF
        END DO

        CALL gc_rrecv(1001+block_pe,                              &
         n_local_rows*global_row_len+2,                           &
         block_pe, info,                                          &
         my_global_buffer(                                        &
         g_datastart_f(1,grid_type,iproc)+                        &
         (g_datastart_f(2,grid_type,iproc)-1)*global_row_len),    &
         icomp)

        CALL unite_coex_files(                                    &
          my_global_buffer(g_datastart_f(1,grid_type,iproc)+      &
         (g_datastart_f(2,grid_type,iproc)-1)*global_row_len),    &
          global_field, num_out, iproc)
      END DO

    END IF ! are we doing 2 stage gather?

  END IF ! mype == global_gather_pe

ELSE

! Normal Gather Operation, without Packing
! 5.0 Set up the second send map

  IF(block_factor <  nproc_y) THEN

    IF(mype == gather_pe .AND. mype /= global_gather_pe) THEN

      send_map_2(s_destination_pe,1) = global_gather_pe

      send_map_2(s_base_address_in_send_array,1) =                &
            g_datastart_f(1,grid_type,mype)+                      &
            (g_datastart_f(2,grid_type,mype)-1)*global_row_len

      send_map_2(s_number_of_elements_in_item,1)=0
      DO kproc=1, block_factor
        IF(row_start_pe+(kproc-1)*nproc_x <  nproc) THEN
          send_map_2(s_number_of_elements_in_item,1)=             &
            send_map_2(s_number_of_elements_in_item,1)+           &
            g_blsize(2,grid_type,mype+(kproc-1)*nproc_x)
        END IF
      END DO

      send_map_2(s_stride_in_send_array,1) =                      &
            global_row_len

      send_map_2(s_element_length,1) =                            &
            global_row_len

      send_map_2(s_base_address_in_recv_array,1) =                &
            g_datastart_f(1,grid_type,mype)+                      &
            (g_datastart_f(2,grid_type,mype)-1)*global_row_len

      send_map_2(s_stride_in_recv_array,1) =                      &
            global_row_len

      n_mess_to_send_2=1

    ELSE

      n_mess_to_send_2=0

    END IF

! 6.0 Set up the second receive map (for PE GLOBAL_GATHER_PE only)

    n_mess_to_recv_2=0

    IF(mype == global_gather_pe) THEN

      iproc=0
      DO jproc=0, nproc_y-1, block_factor

! Compute the block PE for this group of rows, and check it is not
! not the global_gather_pe
        block_pe=jproc*nproc_x
        lproc=block_pe

! Check if the global_gather_pe is in the block - if so
! make sure it is the block_pe
        IF(global_gather_pe >= block_pe .AND.                     &
          global_gather_pe <=                                     &
          MIN(nproc, block_pe+pes_per_block)-1) THEN
           block_pe=global_gather_pe
        END IF

        IF (block_pe  /=  global_gather_pe) THEN

          receive_map_2(r_source_pe,iproc+1) = block_pe

          receive_map_2(r_base_address_in_recv_array,iproc+1) =   &
            g_datastart_f(1,grid_type,block_pe)+                  &
            (g_datastart_f(2,grid_type,block_pe)-1)*global_row_len

          receive_map_2(r_number_of_elements_in_item,iproc+1) = 0
          DO kproc=1, block_factor
            IF(lproc+(kproc-1)*nproc_x <  nproc) THEN
              receive_map_2(r_number_of_elements_in_item,         &
                            iproc+1) =                            &
               receive_map_2(r_number_of_elements_in_item,        &
                            iproc+1) +                            &
               g_blsize(2,grid_type,lproc+(kproc-1)*              &
                        nproc_x)
            END IF
          END DO

          receive_map_2(r_stride_in_recv_array,iproc+1) =         &
            global_row_len

          receive_map_2(r_element_length,iproc+1) =               &
            global_row_len

          receive_map_2(r_base_address_in_send_array,iproc+1) =   &
            g_datastart_f(1,grid_type,block_pe)+                  &
            (g_datastart_f(2,grid_type,block_pe)-1)*global_row_len

          receive_map_2(r_stride_in_send_array,iproc+1) =         &
            global_row_len

          iproc=iproc+1
          n_mess_to_recv_2=n_mess_to_recv_2+1

        END IF ! block_pe  /=  global_gather_pe

      END DO ! gather_pe's

    END IF ! mype == global_gather_pe

! 7.0 Do the exchange of data

    flag=0  ! This is currently ignored at GCG v1.1
    info=gc_none

    CALL gcg_ralltoalle(my_global_buffer,send_map_2,              &
                        n_mess_to_send_2,                         &
                        global_row_len*global_rows,               &
                        global_field,receive_map_2,               &
                        n_mess_to_recv_2,                         &
                        global_row_len*global_rows,               &
                        proc_group,flag,icode)

  END IF ! packing_type is 1

END IF ! are we doing 2 stage gather

! Deallocate the temporary buffer for coex processed data
IF(ALLOCATED(icomp)) DEALLOCATE(icomp)


IF (lhook) CALL dr_hook('GATHER_PACK_FIELD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE gather_pack_field
