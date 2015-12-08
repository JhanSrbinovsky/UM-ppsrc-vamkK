! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine solve_forward_x1
! Purpose : Portable solver for dynamics
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver

! ---------------------------------------------------------
SUBROUTINE solve_forward_x1(                                      &
                             a0_x, a1_x, a2_x,                    &
                             factor_x, f_vector_x,                &
                             row_length, rows, model_levels,      &
                             j_start, j_end)

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
INTEGER, INTENT(IN)    ::  j_start
INTEGER, INTENT(IN)    ::  j_end


REAL,    INTENT(IN)    ::  a1_x(row_length+1,rows,model_levels)
REAL,    INTENT(IN)    ::  a2_x(row_length  ,rows,model_levels)
REAL,    INTENT(  OUT) ::  factor_x(row_length+1,rows,model_levels)
REAL,    INTENT(INOUT) ::  a0_x(row_length+1,rows,model_levels)
REAL,    INTENT(INOUT) ::  f_vector_x(row_length+1,rows,model_levels)

! Include files

! Locals:
INTEGER :: i         ! looper
INTEGER :: j         ! looper
INTEGER :: k         ! looper
INTEGER :: i_start   ! loop lower bound
INTEGER :: i_end     ! loop upper bound

INTEGER :: info                                      ! mpi return code
INTEGER :: my_comm                                   ! communicator
INTEGER :: status(mpl_status_size)                   ! mpi status for recv
INTEGER :: send_status(mpl_status_size,model_levels) ! mpi status for sends
INTEGER :: request(model_levels)                     ! mpi requests

REAL    :: send_data(3,rows,model_levels)     ! Data to send
REAL    :: receive_data(3,rows,model_levels)  ! Received data
                                              ! Index 1: a0_x
                                              ! Index 2: a1_x
                                              ! Index 3: F_vector_x

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------

IF (lhook) CALL dr_hook('SOLVE_FORWARD_X1',zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm,info)

i_start=2
IF (at_extremity(peast)) THEN
  i_end=row_length-1
ELSE
  i_end=row_length
END IF

DO k=1,model_levels

  IF (.NOT. at_extremity(pwest)) THEN
    ! Everyone except Westernmost processor waits to receive
    ! data from their West
    CALL mpl_recv(receive_data(1,j_start,k), 3*(j_end-j_start+1), mpl_real, &
                  neighbour(pwest), k, my_comm, status, info)
  END IF

  DO j=j_start,j_end

    IF (.NOT. at_extremity(pwest)) THEN
      ! First point, use the data I've received from
      ! neighbouring processor
      factor_x(1,j,k)   = a2_x(1,j,k)*receive_data(1,j,k)
      a0_x(1,j,k)       = 1.0/(a0_x(1,j,k)-factor_x(1,j,k)*       &
                                           receive_data(2,j,k))
      f_vector_x(1,j,k) = f_vector_x(1,j,k)-                      &
                          factor_x(1,j,k)*receive_data(3,j,k)
    END IF

    DO i=i_start,i_end
      factor_x(i,j,k)   = a2_x(i,j,k) * a0_x(i-1,j,k)
      a0_x(i,j,k)       = 1.0/(a0_x(i,j,k) - factor_x(i,j,k)*     &
                                             a1_x(i-1,j,k))
      f_vector_x(i,j,k) = f_vector_x(i,j,k) -                     &
                          factor_x(i,j,k)*f_vector_x(i-1,j,k)
    END DO

  END DO

  IF (.NOT. at_extremity(peast)) THEN
    ! Everybody except Easternmost processor sends
    ! data to the processor to their East
    DO j=j_start,j_end
      send_data(1,j,k)=a0_x(row_length,j,k)
      send_data(2,j,k)=a1_x(row_length,j,k)
      send_data(3,j,k)=f_vector_x(row_length,j,k)
    END DO

    CALL mpl_isend(send_data(1,j_start,k),3*(j_end-j_start+1), mpl_real, &
                   neighbour(peast), k, my_comm, request(k), info)
  END IF ! IF (.NOT. at_extremity(PEast))

END DO ! k

! Ensure all comms are finished
IF (.NOT. at_extremity(peast) ) THEN
  CALL mpl_waitall(model_levels, request, send_status, info)
END IF

IF (lhook) CALL dr_hook('SOLVE_FORWARD_X1',zhook_out,zhook_handle)
RETURN
      END SUBROUTINE solve_forward_x1
