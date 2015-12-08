! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Calculate zonal wind for baroclinic wave test case.

REAL FUNCTION eg_baro_u_p(xi1, xi2, eta)

USE parkind1,        ONLY: jpim, jprb       !DrHook
USE yomhook,         ONLY: lhook, dr_hook   !DrHook
USE conversions_mod, ONLY: pi
USE eg_idl_baro_mod, ONLY: u0, eta0, L_channel, Ly, b

IMPLICIT NONE
!
! Description:
!  
! Calculates zonal wind for steady-state solution.
!
! Method:
!
! Baroclinic instability test case (QJRMS 132, 2943--2975).
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! Declarations:

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Function arguments
REAL, INTENT(IN) :: xi1 ! Longitude (radians)
REAL, INTENT(IN) :: xi2 ! Latitude (radians)
REAL, INTENT(IN) :: eta ! p/p_surface

! Local variables
REAL :: eta_v

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_BARO_U_P',zhook_in,zhook_handle)

IF ( L_channel ) THEN
  eg_baro_u_p = -u0 * SIN(pi*xi2/Ly)**2 &
                *LOG(eta)*EXP(-(LOG(eta)/b)**2)
ELSE

  eta_v = (eta - eta0) * pi / 2.0
  eg_baro_u_p = u0 * COS(eta_v)**(3.0/2.0) * SIN(2.0*xi2)**2
END IF

IF (lhook) CALL dr_hook('EG_BARO_U_P',zhook_out,zhook_handle)
RETURN
END FUNCTION eg_baro_u_p
