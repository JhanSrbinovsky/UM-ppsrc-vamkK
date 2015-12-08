! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!----------------------------------------------------------------------
! Subroutine q_pos_hybrid_sub- this version does a hybrid of other methods.
!                              column scavenging for bottom bl_levels and
!                              level scavenging above this
!----------------------------------------------------------------------
SUBROUTINE q_pos_hybrid_sub(q, row_length, rows, off_x, off_y,          &
                            q_levels, variables, bl_levels, qlimit_in)

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

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)    :: off_x          ! x halo size
INTEGER, INTENT(IN)    :: off_y          ! y halo size
INTEGER, INTENT(IN)    :: q_levels       ! number of levels
INTEGER, INTENT(IN)    :: variables      ! number of variables 
INTEGER, INTENT(IN)    :: bl_levels      ! number of boundary layer levels
INTEGER, INTENT(IN)    :: rows           ! number of rows
INTEGER, INTENT(IN)    :: row_length     ! length of rows

REAL, INTENT(IN)       :: qlimit_in      ! minimum value for q
REAL,    INTENT(INOUT) :: q( 1-off_x : row_length+off_x,                &
                             1-off_y : rows+off_y, q_levels/variables,  &
                             variables)

! Local variables
INTEGER                       :: actual_levels
INTEGER                       :: i
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('Q_POS_HYBRID_SUB',zhook_in,zhook_handle)
actual_levels = q_levels/variables

DO i = 1, variables
! DEPENDS ON : q_pos_column_sub
  CALL q_pos_column_sub(q(:,:,1:bl_levels,i), row_length, rows, off_x,    &
                        off_y, bl_levels, qlimit_in)

! DEPENDS ON : q_pos_level_sub
  CALL q_pos_level_sub(q(:,:,bl_levels+1:,i), row_length, rows, off_x,    &
                       off_y, actual_levels-bl_levels, qlimit_in)
END DO

IF (lhook) CALL dr_hook('Q_POS_HYBRID_SUB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE q_pos_hybrid_sub
