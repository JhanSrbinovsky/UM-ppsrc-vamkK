! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!----------------------------------------------------------------------
! Subroutine q_pos_level_sub- this version does a proportional rebalance
!                             within a level
!----------------------------------------------------------------------
SUBROUTINE q_pos_level_sub(q, row_length, rows, off_x, off_y,          &
                           q_levels, qlimit_in)

!
!   PURPOSE:   Removes values of Q below qlimit_in. This version
!              rebalances within a level when possible. This is done
!              by taking from points proportionally the amount of
!              Q they have. This may not always be possible, so there
!              may be conservation issues. Taking moisture from a column
!              is effectively "Method 1" in the old Q_POS.
!              In this version data remains in place so scalability shouldn't
!              be too bad.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90

USE global_2d_sums_mod, ONLY : global_2d_sums
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
INTEGER :: istat             ! GCOM error code
INTEGER :: non_conserve      ! count of non-conserved columns

REAL :: negative             ! column total of negative q
REAL :: positive             ! column total of positive q
REAL :: posneg(q_levels,2)   ! temp array for summing
REAL :: posneg_full(row_length, rows, q_levels, 2)  ! temp array for pos and 
                                                    ! neg values
!REAL :: posneg_rows(rows, q_levels, 2)              ! temp array for rowsums



INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('Q_POS_LEVEL_SUB',zhook_in,zhook_handle)


!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, j, k)     &
!$OMP& SHARED(q_levels, rows, row_length, q, qlimit_in, posneg_full)
DO k = 1, q_levels
  DO j = 1, rows
    DO i = 1, row_length

      ! Convert to a test against 0.0
      q(i,j,k) = q(i,j,k) - qlimit_in

      ! posneg(:,:,:,1) is the positive parts of q
      IF (q(i,j,k) > 0.0) THEN 
        posneg_full(i,j,k,1) = q(i,j,k)
      ELSE
        posneg_full(i,j,k,1) = 0.0
      END IF

      ! posneg(:,:,:,2) is the negative parts of q
      IF (q(i,j,k) < 0.0) THEN 
        posneg_full(i,j,k,2) = -q(i,j,k)
      ELSE
        posneg_full(i,j,k,2) = 0.0
      END IF
    END DO
  END DO
END DO
!$OMP END PARALLEL DO 

CALL global_2d_sums(posneg_full, row_length, rows, 0, 0, 2*q_levels,    &
                    posneg)

non_conserve = 0

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(positive,     &
!$OMP& negative, i, j, k) SHARED(q_levels, posneg, q, qlimit_in, rows,  &
!$OMP& row_length) REDUCTION(+:non_conserve)
DO k = 1, q_levels

  positive = posneg(k,1)
  negative = posneg(k,2)
  
  ! can we be conservative? If so, fix up.
  IF (positive > negative) THEN
    DO j = 1, rows
      DO i = 1, row_length
        IF (q(i,j,k) < 0.0) THEN
          q(i,j,k) = 0.0
        ELSE IF (q(i,j,k) > 0.0) THEN
          ! Take proportionally from those that have
          q(i,j,k) = q(i,j,k) - (q(i,j,k) * negative/positive)
        END IF
      END DO
    END DO

  ELSE    ! non-conserve set whole level to 0.0
    q(:,:,k) = 0.0
    non_conserve = non_conserve + 1
  END IF

  DO j = 1, rows
    DO i = 1, row_length      
      ! Convert back to qlimit_in
      q(i,j,k) = q(i,j,k) + qlimit_in
    END DO
  END DO

END DO   ! k
!$OMP END PARALLEL DO 
    
! Report number of non_conserved levels
IF (non_conserve > 0 .AND. printstatus >= prstatus_diag) THEN
  WRITE (6,'(A,I4,A)') 'Q_POS: unable to conserve on', non_conserve, ' levels'
END IF

IF (lhook) CALL dr_hook('Q_POS_LEVEL_SUB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE q_pos_level_sub
