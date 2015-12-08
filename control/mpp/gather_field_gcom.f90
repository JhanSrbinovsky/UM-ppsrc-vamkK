! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Gathers a field from many processors to one processor

! Subroutine Interface:
SUBROUTINE gather_field_gcom(local_field,global_field,       &
                             local_row_len,local_rows,       &
                             global_row_len,global_rows,     &
                             grid_type,halo_type,            &
                             gather_pe,proc_group,           &
                             icode,cmessage)


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE UM_ParVars
  USE gcom_mod
  IMPLICIT NONE


! Description:
!  Takes a model field that has been decomposed over a group of
!  processors, and gathers the data together so that one processor
!  contains the entire global field.

! Method:
!  A send and receive map is constructed which instructs the GCOM
!  permute operation to do a gather from all processors in the
!  group to the GATHER_PE

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Subroutine Arguments:

INTEGER, INTENT(IN) ::                                            &
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
, gather_pe                                                       &
                   ! IN processor to gather global field to
, proc_group                                                      
                   ! IN group ID of processors involved here

INTEGER, INTENT(OUT) :: icode            ! OUT return code

REAL, INTENT(IN) ::                                               &
  local_field(local_row_len*local_rows)                            
!                        ! IN local part of field
REAL, INTENT(OUT) ::                                              &
  global_field(global_row_len*global_rows)
!                        ! OUT (on PE GATHER_PE) global field

CHARACTER(LEN=80) ::  cmessage         ! OUT error message

! Parameters and Common blocks

! Local variables

INTEGER                                                           &
  send_map(7,1)                                                   &
, n_mess_to_rec

INTEGER, POINTER :: receive_map(:,:)=>NULL()

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

SAVE send_map,receive_map,n_mess_to_rec,                          &
     old_global_row_len,old_global_rows,old_proc_group,           &
     old_gather_pe,old_decomp,                                    &
     old_grid_type,old_halo_type
DATA old_global_row_len,old_global_rows,old_proc_group,           &
     old_gather_pe,old_decomp,                                    &
     old_grid_type,old_halo_type                                  &
   / -1234, -1234, -1234, -1234, -1234, -1234, -1234/

INTEGER  ::                                                         &
  iproc                                                           &
, info                                                            &
, flag

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-------------------------------------------------------

! 0.0 Can we use the same send/receive map that we calculated
!     last time round?

IF (lhook) CALL dr_hook('GATHER_FIELD_GCOM',zhook_in,zhook_handle)

IF (.NOT.ASSOCIATED(receive_map)) &
    ALLOCATE(receive_map(7,nproc_max))

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
    WRITE(6,*) 'GATHER_FIELD_GCOM : Bad field type'
    WRITE(6,*) 'Field will not be gathered.'
    cmessage='GATHER_FIELD_GCOM : Bad field type'
    icode=1
    GO TO 9999
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

! 3.0 Set up the receive map (for PE GATHER_PE only)

  n_mess_to_rec=0

  IF (mype  ==  gather_pe) THEN
    DO iproc=0,nproc-1
      receive_map(r_source_pe,iproc+1) = iproc

      receive_map(r_base_address_in_recv_array,iproc+1) =         &
        g_datastart_f(1,grid_type,iproc)+                         &
        (g_datastart_f(2,grid_type,iproc)-1)*global_row_len

      receive_map(r_number_of_elements_in_item,iproc+1) =         &
        g_blsize(2,grid_type,iproc)

      receive_map(r_stride_in_recv_array,iproc+1) =               &
        global_row_len

      receive_map(r_element_length,iproc+1) =                     &
        g_blsize(1,grid_type,iproc)

      receive_map(r_base_address_in_send_array,iproc+1) =         &
        halosize(2,halo_type)*                                    &
        g_lasize(1,grid_type,halo_type,iproc)+                    &
        halosize(1,halo_type)+1

      receive_map(r_stride_in_send_array,iproc+1) =               &
        g_lasize(1,grid_type,halo_type,iproc)

    END DO
    n_mess_to_rec=nproc
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

CALL gcg_ralltoalle(local_field,send_map,1,                       &
                    local_row_len*local_rows,                     &
                    global_field,receive_map,n_mess_to_rec,       &
                    global_row_len*global_rows,                   &
                    proc_group,flag,info)

 9999 CONTINUE

IF (lhook) CALL dr_hook('GATHER_FIELD_GCOM',zhook_out,zhook_handle)
RETURN
END SUBROUTINE gather_field_gcom
