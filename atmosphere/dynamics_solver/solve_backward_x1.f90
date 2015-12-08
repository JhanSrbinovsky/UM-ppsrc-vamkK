! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine solve_backward_x1
! Purpose : Portable solver for dynamics
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver


SUBROUTINE solve_backward_x1(                                     &
                              soln, a0_x, a1_x,                   &
                              row_length, rows, model_levels,     &
                              j_start, j_end, off_x, off_y)

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
INTEGER, INTENT(IN)    ::  off_x
INTEGER, INTENT(IN)    ::  off_y

REAL,    INTENT(IN)    ::  a0_x(row_length+1,rows,model_levels)
REAL,    INTENT(IN)    ::  a1_x(row_length+1,rows,model_levels)
REAL,    INTENT(INOUT) ::  soln(1-off_x:row_length+off_x+1,       &
                                1-off_y:rows+off_y,               &
                                model_levels)  

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

INTEGER, SAVE :: my_type         = -1234    ! mpi datatype 
INTEGER, SAVE :: j_start_old     = -1234    ! previous j_start
INTEGER, SAVE :: j_end_old       = -1234    ! previous j_end
INTEGER, SAVE :: off_x_old       = -1234    ! previous off_x
INTEGER, SAVE :: row_length_old  = -1234    ! previous row_length

REAL :: receive_data(rows,model_levels)     ! Received data

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------

IF (lhook) CALL dr_hook('SOLVE_BACKWARD_X1',zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm,info)

! Test if we need a new type
IF (j_start    /= j_start_old .OR. &
    j_end      /= j_end_old   .OR. &
    off_x      /= off_x_old   .OR. &
    row_length /= row_length_old) THEN

  ! If not first call free up old type
  IF (my_type /= -1234) THEN
    CALL mpl_type_free(my_type, info)
  END IF

  CALL mpl_type_vector(j_end-j_start+1, 1, row_length+(2*off_x)+1, &
                       mpl_real, my_type, info)
  CALL mpl_type_commit(my_type, info)

  j_start_old = j_start
  j_end_old   = j_end
  row_length_old = row_length
  off_x_old = off_x
END IF


i_end=1
IF (at_extremity(peast)) THEN
  i_start=row_length-2
ELSE
  i_start=row_length-1
END IF

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

DO k=1,model_levels


!$OMP MASTER

  IF (.NOT. at_extremity(peast)) THEN
    ! Everyone except Easternmost processor waits to receive
    ! data from their East
    CALL mpl_recv(receive_data(j_start,k), j_end-j_start+1, mpl_real, &
                  neighbour(peast), k, my_comm, status, info)
  END IF

!$OMP END MASTER 

!$OMP BARRIER

  IF (.NOT. at_extremity(peast)) THEN

!$OMP DO SCHEDULE(STATIC)
    DO j = j_start, j_end
      ! First point, use the data I've received from
      ! neighbouring processor
      i=i_start+1
      soln(i,j,k)=a0_x(i,j,k)*(soln(i,j,k)-                       &
                  a1_x(i,j,k)*                                    &
                  receive_data(j,k))
    END DO
!$OMP END DO NOWAIT
  END IF


!$OMP DO SCHEDULE(STATIC)
  DO j = j_start, j_end
    DO i=i_start,i_end,-1
      soln(i,j,k)=a0_x(i,j,k)*(soln(i,j,k)-                       &
                               a1_x(i,j,k)*                       &
                               soln(i+1,j,k))
    END DO
  END DO
!$OMP END DO  

!$OMP MASTER

  IF (.NOT. at_extremity(pwest)) THEN
     ! Everybody except Westernmost processor sends
     ! data to the processor to their West
    CALL mpl_isend(soln(1,j_start,k), 1, my_type, neighbour(pwest), &
                   k, my_comm, request(k), info)

  END IF ! IF (.NOT. at_extremity(PWest))

!$OMP END MASTER 

END DO ! k

!$OMP END PARALLEL

! Ensure all comms are finished.
IF (.NOT. at_extremity(pwest) ) THEN
  CALL mpl_waitall(model_levels, request, send_status, info)
END IF

IF (lhook) CALL dr_hook('SOLVE_BACKWARD_X1',zhook_out,zhook_handle)
RETURN
END SUBROUTINE solve_backward_x1
