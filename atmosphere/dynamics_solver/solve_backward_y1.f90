! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine solve_backward_y1
! Purpose : Portable solver for dynamics
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver

! ---------------------------------------------------------
SUBROUTINE solve_backward_y1(                                     &
                              soln, a0_y, a1_y,                   &
                              row_length, rows, model_levels,     &
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

REAL,    INTENT(IN)    :: a0_y(row_length,rows,model_levels) 
REAL,    INTENT(IN)    :: a1_y(row_length,rows,model_levels) 
REAL,    INTENT(INOUT) :: soln(1-off_x:row_length+off_x,          &
                               1-off_y:rows+off_y,                &
                               model_levels) 


! Locals:
INTEGER :: i        ! Looper
INTEGER :: j        ! Looper
INTEGER :: k        ! Looper
INTEGER :: j_start  ! Loop lower bound
INTEGER :: j_end    ! Loop upper bound

INTEGER :: info                                      ! mpi return code
INTEGER :: my_comm                                   ! communicator
INTEGER :: status(mpl_status_size)                   ! mpi status for recv
INTEGER :: send_status(mpl_status_size,model_levels) ! mpi status for sends
INTEGER :: request(model_levels)                     ! mpi requests

REAL    :: receive_data(row_length,model_levels)     ! Received data

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SOLVE_BACKWARD_Y1',zhook_in,zhook_handle)

!------------------------------------------------------------
CALL gc_get_communicator(my_comm,info)

j_end=1
j_start=rows-1
IF (at_extremity(psouth)) j_end=2
IF (at_extremity(pnorth)) j_start=rows-2

DO k=1,model_levels

  IF (.NOT. at_extremity(pnorth)) THEN
    ! Everyone except Northernmost processor waits to receive
    ! data from their North
    CALL mpl_recv(receive_data(1,k), row_length, mpl_real, &
                  neighbour(pnorth), k, my_comm, status, info)
  END IF

  IF (.NOT. at_extremity(pnorth)) THEN
    ! First row, use the data I've received from
    ! neighbouring processor

    j=rows
    DO i=1,row_length
      soln(i,j,k) = a0_y(i,j,k)*(soln(i,j,k)-                     &
                                 a1_y(i,j,k)*                     &
                                 receive_data(i,k))
    END DO
  END IF

  DO j=j_start,j_end,-1
    DO i=1,row_length
      soln(i,j,k)=a0_y(i,j,k)*(soln(i,j,k)-                       &
                               a1_y(i,j,k)*soln(i,j+1,k))
    END DO
  END DO

  IF (.NOT. at_extremity(psouth)) THEN
     ! Everybody except Southernmost processor sends
     ! data to the processor to their South
    CALL mpl_isend(soln(1,1,k), row_length, mpl_real, &
                   neighbour(psouth), k, my_comm, request(k), info)
  END IF ! IF (.NOT. at_extremity(PSouth))

END DO ! k

! Ensure all comms finished
IF (.NOT. at_extremity(psouth)) THEN
  CALL mpl_waitall(model_levels, request, send_status, info)
END IF

      IF (lhook) CALL dr_hook('SOLVE_BACKWARD_Y1',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE solve_backward_y1
