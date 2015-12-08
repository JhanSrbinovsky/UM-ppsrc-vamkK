! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Returns an ensemble vertical velocity profile for shallow cumulus
! 

REAL FUNCTION F_W(X)   ! ENSEMBLE VERTICAL VELOCITY PROFILE

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE
!
! Description:
!   Calculates non-dimensional vertical velocity profile for shallow cumulus
!
! Method:
!   Derived from large-eddy simulations.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!-----------------------------------------------------------------------------

REAL, INTENT(IN) :: x             ! non-dimesional height in cumulus layer.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------------
IF (lhook) CALL dr_hook('F_W',zhook_in,zhook_handle) 

IF(x <= 1.0) THEN
  f_w=6.0*x
ELSE
  f_w=6.0
END IF

IF (lhook) CALL dr_hook('F_W',zhook_out,zhook_handle)
RETURN
END FUNCTION F_W
