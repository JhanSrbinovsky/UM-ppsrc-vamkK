! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!----------------------------------------------------------------------
! Subroutine q_pos_diag_pr - diagnostic prints for q_pos
!----------------------------------------------------------------------
SUBROUTINE q_pos_diag_pr(q, row_length, rows, levels, off_x, off_y,     &
                         qpos_diag_lim, ident)

!
!   PURPOSE:  A diagnostic print routine. Checks values of q to see
!             if any lie below the given limit. If so prints some 
!             helpful diagnostic information
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

  USE ereport_mod, ONLY : ereport
  USE UM_ParVars
  IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)    :: off_x          ! x halo size
INTEGER, INTENT(IN)    :: off_y          ! y halo size
INTEGER, INTENT(IN)    :: levels         ! number of levels
INTEGER, INTENT(IN)    :: rows           ! number of rows
INTEGER, INTENT(IN)    :: row_length     ! length of rows


REAL,    INTENT(IN)    :: qpos_diag_lim ! minimum value for q
REAL,    INTENT(IN)    :: q( 1-off_x : row_length+off_x,                &
                             1-off_y : rows+off_y, levels)

CHARACTER(LEN=*), INTENT(IN) :: ident    ! identifier for qpos call

! Local variables
INTEGER :: i                 ! looper
INTEGER :: j                 ! looper
INTEGER :: k                 ! looper
INTEGER :: small_q_count     ! count of number of q points below qpos_diag_lim
INTEGER :: istat             ! GCOM error code
INTEGER :: icode             ! Warning code

REAL    :: q_min             ! smallest q below qpos_diag_lim

CHARACTER (LEN=256) :: c_message      ! error message
CHARACTER (LEN=*), PARAMETER :: routinename='q_pos_diag_pr'

! Include files

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('Q_POS_DIAG_PR',zhook_in,zhook_handle)

q_min         = qpos_diag_lim
small_q_count = 0

DO k = 1, levels
  DO j = 1, rows
    DO i = 1, row_length
      IF (q(i,j,k) < qpos_diag_lim) THEN
        ! This is new minimum value, and 1 to add to count of small q values
        q_min = MIN(q(i,j,k), q_min)
        small_q_count = small_q_count + 1
      END IF
    END DO
  END DO
END DO

! Get the minimum Q across all processors
CALL gc_isum(1, nproc, istat, small_q_count)
CALL gc_rmin(1, nproc, istat, q_min)

! is min q below the limit?
IF (q_min < qpos_diag_lim) THEN
  WRITE(c_message,*) 'QPOS: In ',ident,' call there are ',              &
                    small_q_count, ' values less than ', qpos_diag_lim, &
                   '. The smallest is ', q_min
  icode = -100

  CALL ereport(routinename, icode, c_message)
END IF


IF (lhook) CALL dr_hook('Q_POS_DIAG_PR',zhook_out,zhook_handle)
RETURN
END SUBROUTINE q_pos_diag_pr
