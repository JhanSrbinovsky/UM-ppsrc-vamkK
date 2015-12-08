! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_destroy_vert_damp_mod

! Subroutine: eg_destroy_vert_damp
!
! Description: deallocated mu_w at the end of last time-step (called
!              from atm_step)
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: DYNAMICS ADVECTION
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

CONTAINS

SUBROUTINE eg_destroy_vert_damp()

USE eg_vert_damp_mod, ONLY : mu_w

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

! local variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('EG_DESTROY_VERT_DAMP',zhook_in,zhook_handle)

DEALLOCATE (mu_w)

IF (lhook) CALL dr_hook('EG_DESTROY_VERT_DAMP',zhook_out,zhook_handle)
RETURN

END SUBROUTINE eg_destroy_vert_damp

END MODULE eg_destroy_vert_damp_mod
