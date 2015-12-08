! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!  Takes 1 or more levels of a model field that have been decomposed
!  over a group of processors, and gathers the data together so that
!  one complete global level is contained on one processor.  

! Method:
!  This routine copies all local data (all levels) into a single
!  array which is then sent to the gathering CPU where it is unpacked.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP


! Subroutine Interface:
SUBROUTINE gather_field_ml(                                       &
    local_field,    global_field,                                 &
    local_row_len,  local_rows,   local_levs,                     &
    global_row_len, global_rows,  global_levs,                    &
    pe_for_level,                                                 &
    fld_type,       halo_type)

USE mpl, ONLY :           &
         mpl_real,        &
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
INTEGER, INTENT(IN) :: pe_for_level(local_levs)  ! PE to gather
                                                 ! each level to


! Decomposed data
REAL, INTENT(IN)    :: local_field( local_row_len,                &
                                    local_rows, local_levs )
! Gathered data
REAL, INTENT(OUT)   :: global_field( global_row_len,              &
                                     global_rows, global_levs )

! Parameters and Common blocks


! Local variables

INTEGER :: i       ! loop index  - cols
INTEGER :: j       ! loop index  - rows
INTEGER :: k       ! loop index  - levels
INTEGER :: iproc   ! loop index  - processors
INTEGER :: halo_x  ! halo size - x
INTEGER :: halo_y  ! halo size - y
INTEGER :: local_row_len_nh  ! row length without halo
INTEGER :: local_rows_nh     ! row count without halo
INTEGER :: pos               ! position in array

INTEGER :: levs_to_send(0 : nproc-1) ! num of levels to each PE
INTEGER :: kpos(0 : nproc-1)         ! array position
INTEGER :: send_size(0 : nproc-1)    ! size of sent data
INTEGER :: recv_size(0 : nproc-1)    ! size of received data

INTEGER :: ierr                      ! error code
INTEGER :: statu(mpl_status_size)   ! MPL status
INTEGER :: my_comm                   ! MPL communicator


! Array to hold all local data - note contains space for halos
! that won't be used
REAL :: local_buffer(local_row_len * local_rows * global_levs,    &
                     0 : nproc-1)

! Array to hold received data - much too big but best I can
! think of
REAL :: recv_buff(global_row_len * global_rows * global_levs,     &
                  0 : nproc -1)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!-------------------------------------------------------
! 0) Setup  - get communicator from GCOM
!-------------------------------------------------------
IF (lhook) CALL dr_hook('GATHER_FIELD_ML',zhook_in,zhook_handle)
CALL gc_get_communicator(my_comm, ierr)

!-------------------------------------------------------
! 1) Copy data from multiple levels into 1 buffer
!    for sending - this reduces the number of messages
!    sent.
!-------------------------------------------------------
halo_x = halosize(1, halo_type)
halo_y = halosize(2, halo_type)
! Non halo sizes
local_row_len_nh = local_row_len - (2 * halo_x)
local_rows_nh = local_rows - (2 * halo_y)

DO iproc = 0, nproc - 1
  kpos(iproc) = 0
END DO


! Copy local_field into local_buffer with halos stripped off
! and with contiguous levels (ie 1 message per PE)
DO k = 1, local_levs
  DO j = 1+halo_y, local_rows - halo_y
    DO i = 1+halo_x, local_row_len - halo_x
      pos = i - halo_x +                                          &
            (j - halo_y - 1) * local_row_len_nh +                 &
            kpos(pe_for_level(k)) * local_rows_nh *               &
                                    local_row_len_nh
      local_buffer(pos,pe_for_level(k)) =                         &
      local_field(i,j,k)
    END DO
  END DO
  kpos(pe_for_level(k)) = kpos(pe_for_level(k)) + 1
END DO


!-------------------------------------------------------
! 2) Find sizes for send/recv and do the communications
!    Use MPL_Sendrecv to pair up comms.
!-------------------------------------------------------
! Find sizes to send and receive
levs_to_send(:) = 0
DO k = 1, local_levs
  levs_to_send(pe_for_level(k)) = levs_to_send(pe_for_level(k))+1
END DO

DO iproc = 0, nproc - 1
  send_size(iproc) = local_row_len_nh * local_rows_nh *           &
                     levs_to_send(iproc)
  recv_size(iproc) = g_blsize(1,fld_type,iproc) *                 &
                     g_blsize(2,fld_type,iproc) *                 &
                     levs_to_send(mype)
END DO

! Do communications using MPL
DO iproc = 0, nproc - 1
  CALL mpl_sendrecv( local_buffer(1,iproc), send_size(iproc),     &
                     mpl_real, iproc, 999, recv_buff(1,iproc),    &
                     recv_size(iproc), mpl_real, iproc, 999,      &
                     my_comm, statu, ierr)

END DO

!-------------------------------------------------------
! 3) Copy data from recv buffers into final data field
!-------------------------------------------------------
DO iproc = 0, nproc - 1
  DO k = 1, levs_to_send(mype)
    DO j = 1, g_blsize(2,fld_type,iproc)
      DO i = 1, g_blsize(1,fld_type,iproc)
        global_field( g_datastart_f(1,fld_type,iproc) + i - 1,    &
                      g_datastart_f(2,fld_type,iproc) + j - 1 ,k) &
        = recv_buff(i + (j - 1) * g_blsize(1,fld_type,iproc) +    &
                    (k - 1) * g_blsize(1,fld_type,iproc) *        &
                    g_blsize(2,fld_type,iproc), iproc)
      END DO
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook('GATHER_FIELD_ML',zhook_out,zhook_handle)
RETURN
END SUBROUTINE gather_field_ml
