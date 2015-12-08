! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!----------------------------------------------------------------------
! Subroutine q_pos_reset_sub - this version does a simple reset
!----------------------------------------------------------------------
SUBROUTINE q_pos_reset_sub (q, row_length, rows, off_x, off_y,          &
                            q_levels, qlimit_in)

!
!   PURPOSE:   Removes values of Q below qlimit_in. This version
!              is a simple one that just reinstates values below the
!              minimum to that minimum. There are conservation issues
!              with this, but it may be of use for a forecast as it is
!              cheap and scales well
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)    :: off_x          ! x halo size
INTEGER, INTENT(IN)    :: off_y          ! y halo size
INTEGER, INTENT(IN)    :: q_levels       ! number of levels
INTEGER, INTENT(IN)    :: rows           ! number of rows
INTEGER, INTENT(IN)    :: row_length     ! length of rows


REAL,    INTENT(IN)    :: qlimit_in      ! minimum value for q
REAL,    INTENT(INOUT) :: q( 1-off_x : row_length+off_x,                &
                             1-off_y : rows+off_y, q_levels)

! Local variables
INTEGER :: i                 ! looper
INTEGER :: j                 ! looper
INTEGER :: k                 ! looper

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('Q_POS_RESET_SUB',zhook_in,zhook_handle)

! Simple 3d loop. Reset any values that are below minimum
DO k = 1, q_levels
  DO j = 1, rows
    DO i = 1, row_length
      IF (q(i,j,k) < qlimit_in) THEN
        q(i,j,k) = qlimit_in
      END IF
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook('Q_POS_RESET_SUB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE q_pos_reset_sub
