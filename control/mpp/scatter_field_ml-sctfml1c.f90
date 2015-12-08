! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
! This routine is a direct inverse of GATHER_field_ML. It takes
! full fields (selected levels - chosen by way of a map) and
! scatters them to decomposed data.

! Method:
!  This routine copies all global data (selected levels) into a single
!  array which is then sent to all CPUs where it is unpacked into
!  local data. MPL used for speed.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Subroutine Interface:
SUBROUTINE scatter_field_ml(                                      &
    local_field,    global_field,                                 &
    local_row_len,  local_rows,    local_levs,                    &
    global_row_len, global_rows,   global_levs,                   &
    pe_for_level,                                                 &
    fld_type,       halo_type)

USE mpl, ONLY :                                                   &
         mpl_real,                                                &
         mpl_status_size

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
IMPLICIT NONE

! Subroutine Arguments:

INTEGER, INTENT(IN) :: local_row_len   ! local field row length
INTEGER, INTENT(IN) :: local_rows      ! local field rows
INTEGER, INTENT(IN) :: local_levs      ! local field levels
INTEGER, INTENT(IN) :: global_row_len  ! global field row length
INTEGER, INTENT(IN) :: global_rows     ! global field rows
INTEGER, INTENT(IN) :: global_levs     ! global field levels
INTEGER, INTENT(IN) :: fld_type        ! field type of grid
INTEGER, INTENT(IN) :: halo_type       ! halo type of field
INTEGER, INTENT(IN) :: pe_for_level(local_levs)  ! PE to scatter
                                                 ! level from

! Scattered data
REAL, INTENT(OUT)   :: local_field( local_row_len,                &
                                    local_rows, local_levs )
! Original data
REAL, INTENT(IN)    :: global_field( global_row_len,              &
                                     global_rows, global_levs )


! Parameters and Common blocks


! Local variables

INTEGER :: i                 ! loop index  - cols
INTEGER :: j                 ! loop index  - rows
INTEGER :: k                 ! loop index  - levels
INTEGER :: iproc             ! loop index  - processors
INTEGER :: halo_x            ! halo size - x
INTEGER :: halo_y            ! halo size - y
INTEGER :: local_row_len_nh  ! local row length without halos
INTEGER :: local_rows_nh     ! local rows without halos
INTEGER :: pos               ! buffer position

INTEGER :: levs_to_send(0 : nproc-1) ! num of levs to send
INTEGER :: kpos(0 : nproc-1)         ! buffer position
INTEGER :: send_size(0 : nproc-1)    ! size to send
INTEGER :: recv_size(0 : nproc-1)    ! size to receive

INTEGER :: ierr                      ! error flag
INTEGER :: statu(mpl_status_size)   ! MPL status
INTEGER :: my_comm                   ! MPL Communicator


! Array to hold all local data - note contains space for halos
! that won't be used
REAL :: local_buffer(local_row_len * local_rows * global_levs,    &
                     0 : nproc-1)

! Array to hold received data - much too big but can't think of
! how to make it better sized at the moment.
REAL :: send_buff(global_row_len * global_rows * global_levs,     &
                  0 : nproc -1)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!-------------------------------------------------------
! 0) Calculate a few sizes I'll need later
!-------------------------------------------------------
IF (lhook) CALL dr_hook('SCATTER_FIELD_ML',zhook_in,zhook_handle)

halo_x = halosize(1, halo_type)
halo_y = halosize(2, halo_type)

! Non halo sizes
local_row_len_nh = local_row_len - (2 * halo_x)
local_rows_nh = local_rows - (2 * halo_y)

! Find sizes to send and receive
levs_to_send(:) = 0
DO k = 1, local_levs
  levs_to_send(pe_for_level(k)) = levs_to_send(pe_for_level(k))+1
END DO

!-------------------------------------------------------
! 0) Setup  - get communicator from GCOM
!-------------------------------------------------------
CALL gc_get_communicator(my_comm, ierr)

!-------------------------------------------------------
! 1) Copy data from global fields into send buffer
!-------------------------------------------------------
DO iproc = 0, nproc - 1
  DO k = 1, levs_to_send(mype)
    DO j = 1, g_blsize(2,fld_type,iproc)
      DO i = 1, g_blsize(1,fld_type,iproc)
        send_buff(i + (j - 1) * g_blsize(1,fld_type,iproc) +      &
                    (k - 1) * g_blsize(1,fld_type,iproc) *        &
                    g_blsize(2,fld_type,iproc), iproc)            &
        = global_field( g_datastart_f(1,fld_type,iproc) + i - 1,  &
                      g_datastart_f(2,fld_type,iproc) + j - 1 ,k)
      END DO
    END DO
  END DO
END DO

!-------------------------------------------------------
! 2) Find sizes for send/recv and do the communications
!    Use MPL_Sendrecv to pair up comms.
!-------------------------------------------------------
DO iproc = 0, nproc - 1
  recv_size(iproc) = local_row_len_nh * local_rows_nh *           &
                     levs_to_send(iproc)
  send_size(iproc) = g_blsize(1,fld_type,iproc) *                 &
                     g_blsize(2,fld_type,iproc) *                 &
                     levs_to_send(mype)
END DO

! Do communications using MPL directly
DO iproc = 0, nproc - 1
  CALL mpl_sendrecv( send_buff(1,iproc), send_size(iproc),        &
                     mpl_real, iproc, 999, local_buffer(1,iproc), &
                     recv_size(iproc), mpl_real, iproc, 999,      &
                     my_comm, statu, ierr)

END DO

!-------------------------------------------------------
! 3) Copy data from received buffer into proper
!    decomposed data locations
!-------------------------------------------------------
DO iproc = 0, nproc - 1
  kpos(iproc) = 0
END DO


! Copy local_buffer (no halos) into local_field (with halos)
! Need to get levels right too.
DO k = 1, local_levs
  DO j = 1+halo_y, local_rows - halo_y
    DO i = 1+halo_x, local_row_len - halo_x
      pos = i - halo_x +                                          &
            (j - halo_y - 1) * local_row_len_nh +                 &
            kpos(pe_for_level(k)) * local_rows_nh *               &
                                    local_row_len_nh

      local_field(i,j,k)                                          &
      = local_buffer(pos,pe_for_level(k))
    END DO
  END DO
  kpos(pe_for_level(k)) = kpos(pe_for_level(k)) + 1
END DO


IF (lhook) CALL dr_hook('SCATTER_FIELD_ML',zhook_out,zhook_handle)
RETURN
END SUBROUTINE scatter_field_ml
