! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine solve_forward_y2
! Purpose : Portable solver for dynamics
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver

! ---------------------------------------------------------

SUBROUTINE solve_forward_y2(                                      &
                             soln,factor_y,                       &
                             row_length, rows, model_levels,      &
                             off_x,off_y)

USE mpl, ONLY : &
    mpl_real,   &
    mpl_status_size

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
IMPLICIT NONE

! Arguments:
INTEGER, INTENT(IN)    ::  row_length
INTEGER, INTENT(IN)    ::  rows
INTEGER, INTENT(IN)    ::  model_levels
INTEGER, INTENT(IN)    ::  off_x
INTEGER, INTENT(IN)    ::  off_y

REAL,    INTENT(IN)    ::  factor_y(row_length,rows,model_levels)
REAL,    INTENT(INOUT) ::  soln(1-off_x:row_length+off_x,         &
                                1-off_y:rows+off_y,               &
                                model_levels)

! Include files

! Locals:
INTEGER :: i         ! looper
INTEGER :: j         ! looper
INTEGER :: k         ! looper
INTEGER :: j_start   ! loop lower bound
INTEGER :: j_end     ! loop upper bound

INTEGER :: info                                      ! mpi return code
INTEGER :: my_comm                                   ! communicator
INTEGER :: status(mpl_status_size)                   ! mpi status for recv
INTEGER :: send_status(mpl_status_size,model_levels) ! mpi status for sends
INTEGER :: request(model_levels)                     ! mpi requests


REAL :: receive_data(row_length,model_levels)  ! Received data

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------

IF (lhook) CALL dr_hook('SOLVE_FORWARD_Y2',zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm,info)

j_start=2
j_end=rows
IF (at_extremity(psouth)) j_start=3
IF (at_extremity(pnorth)) j_end=rows-1

DO k=1,model_levels

  IF (.NOT. at_extremity(psouth)) THEN
    ! Everyone except Southernmost processor waits to receive
    ! data from their South
    CALL mpl_recv(receive_data(1,k), row_length, mpl_real, &
                  neighbour(psouth), k, my_comm, status, info)
  END IF

  IF (.NOT. at_extremity(psouth)) THEN
    ! First row, use the data I've received from
    ! neighbouring processor

    DO i=1,row_length
      soln(i,1,k)=soln(i,1,k)-factor_y(i,1,k)*receive_data(i,k)
    END DO
  END IF

  DO j=j_start,j_end
    DO i=1,row_length
      soln(i,j,k)=soln(i,j,k)-factor_y(i,j,k)*soln(i,j-1,k)
    END DO
  END DO

  IF (.NOT. at_extremity(pnorth)) THEN
    ! Everybody except Northernmost processor sends
    ! data to the processor to their North
    CALL mpl_isend(soln(1,rows,k), row_length, mpl_real, &
                   neighbour(pnorth), k, my_comm, request(k), info)
  END IF ! IF (.NOT. at_extremity(PNorth))

END DO ! k

! Ensure all comms have finished
IF (.NOT. at_extremity(pnorth)) THEN
  CALL mpl_waitall(model_levels, request, send_status, info)
END IF

      IF (lhook) CALL dr_hook('SOLVE_FORWARD_Y2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE solve_forward_y2
