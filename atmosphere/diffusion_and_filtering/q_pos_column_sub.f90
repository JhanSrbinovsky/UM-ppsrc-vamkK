! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!----------------------------------------------------------------------
! Subroutine q_pos_column_sub - this version does a proportional rebalance
!                               within a column
!----------------------------------------------------------------------
SUBROUTINE q_pos_column_sub(q, row_length, rows, off_x, off_y,          &
                            q_levels, qlimit_in)

!
!   PURPOSE:   Removes values of Q below qlimit_in. This version
!              rebalances within a column when possible. This is done
!              by taking from levels proportionally from the amount of
!              Q they have. This may not always be possible, so there
!              may be conservation issues. Taking moisture from a column
!              may introduce biases so needs to be considered carefully.
!              There is little communication used so this should scale 
!              well.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE PrintStatus_mod
USE UM_ParVars
IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)    :: off_x          ! x halo size
INTEGER, INTENT(IN)    :: off_y          ! y halo size
INTEGER, INTENT(IN)    :: q_levels       ! number of levels
INTEGER, INTENT(IN)    :: rows           ! number of rows
INTEGER, INTENT(IN)    :: row_length     ! length of rows

REAL, INTENT(IN)       :: qlimit_in      ! minimum value for q
REAL,    INTENT(INOUT) :: q( 1-off_x : row_length+off_x,                &
                             1-off_y : rows+off_y, q_levels)

! Local variables
INTEGER :: i                 ! looper
INTEGER :: j                 ! looper
INTEGER :: k                 ! looper
INTEGER :: l                 ! looper
INTEGER :: istat             ! GCOM error code
INTEGER :: non_conserve      ! count of non-conserved columns

REAL :: negative             ! column total of negative q
REAL :: positive             ! column total of positive q

LOGICAL :: q_fix(row_length, rows) ! flag of colums to fix
 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('Q_POS_COLUMN_SUB',zhook_in,zhook_handle)


! Set logical flag if a column has points that need correcting

q_fix(:,:) = .FALSE.
non_conserve = 0

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, positive, negative)

!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i = 1, row_length
    DO k = 1, q_levels
      IF (q(i,j,k) < qlimit_in) THEN ! Need to fix this column
        q_fix(i,j) = .TRUE.
        q(i,j,:) = q(i,j,:) - qlimit_in  ! normalise 
        EXIT
      END IF
    END DO
  END DO
END DO
!$OMP END DO


!$OMP DO SCHEDULE(STATIC) REDUCTION(+:non_conserve) 
DO j = 1, rows
  DO i = 1, row_length
    positive = 0.0
    negative = 0.0
    IF ( q_fix(i,j) ) THEN
      DO k = 1, q_levels
        IF (q(i,j,k) > 0.0) positive = positive + q(i,j,k)
        IF (q(i,j,k) < 0.0) negative = negative - q(i,j,k)
      END DO

      ! can we be conservative? If so, fix up.
      IF (positive > negative) THEN
        DO k = 1, q_levels
          IF (q(i,j,k) < 0.0) THEN
            q(i,j,k) = 0.0
          ELSE IF (q(i,j,k) > 0.0) THEN
            ! Take proportionally from those that have
            q(i,j,k) = q(i,j,k) - (q(i,j,k) * negative/positive)
          END IF
        END DO
      ELSE
        ! all column is set to 0
        q(i,j,:) = 0.0
        non_conserve = non_conserve + 1
      END IF

      ! go back to q_limit rather than 0
      q(i,j,:) = q(i,j,:) + qlimit_in

    END IF
  END DO
END DO
!$OMP END DO

!$OMP END PARALLEL

! Find total number of non-conserved columns
IF (printstatus >= prstatus_diag) THEN
  CALL gc_isum(1, nproc, istat, non_conserve)
  IF (non_conserve > 0  .AND. mype == 0) THEN
    WRITE (6,'(A,I8,A)') 'Q_POS: unable to conserve in ', non_conserve, &
                         ' columns'
  END IF
END IF

IF (lhook) CALL dr_hook('Q_POS_COLUMN_SUB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE q_pos_column_sub
