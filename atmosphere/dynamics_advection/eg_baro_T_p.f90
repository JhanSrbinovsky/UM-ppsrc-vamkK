! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates temperature for baroclinic wave test

REAL FUNCTION eg_baro_t_p(xi1, xi2, eta)

USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook
USE atmos_constants_mod, ONLY: R
USE earth_constants_mod, ONLY: Earth_radius, g, omega
USE conversions_mod,     ONLY: pi
USE eg_idl_baro_mod

IMPLICIT NONE
!
! Description:
!                Used in initialisation of baroclinic test case.  
!
! Method:
!  
!   Coded to closely follow formulae as given in
!   Jablonowski & Williamson (2006), QJRMS 132, pp.2943--2975.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Function arguments
REAL, INTENT(IN) :: xi1 ! longitude (radians)
REAL, INTENT(IN) :: xi2 ! latitude  (radians)
REAL, INTENT(IN) :: eta ! p/p_surface

! Local
REAL :: eta_v
REAL :: av_temp
REAL :: const
REAL :: wind_term
REAL :: rotation_term
REAL :: phi_prime

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_BARO_T_P',zhook_in,zhook_handle)

IF ( L_channel ) THEN
  av_temp = T0 * eta ** (R * eg_lapse / g)
  
  const = 2.0*pi*xi2/Ly
  phi_prime = 0.5*u0*((f0-beta0*Ly/2.0)*(xi2-Ly/2.0-Ly/(2.0*pi)*SIN(const)) &
            + 0.5*beta0*(xi2**2-Ly*xi2/pi*SIN(const)                        &
            - Ly**2/(2.0*pi**2)*COS(const) - Ly**2/3.0 - Ly**2/(2.0*pi**2)))
  
  eg_baro_t_p = av_temp                                                     &
              + phi_prime/R*(2.0/b**2*LOG(eta)**2-1.0)*EXP(-(LOG(eta)/b)**2)
  
ELSE

  eta_v = (eta - eta0) * pi / 2.0

  IF (eta >= eta_t) THEN
    av_temp = T0 * eta ** (R * eg_lapse / g)
  ELSE IF (eta < eta_t) THEN
    av_temp = T0 * eta ** (R * eg_lapse / g) + DT * (eta_t - eta) ** 5
  END IF

  const=(3.0/4.0)*(pi*u0/R)*SIN(eta_v)*COS(eta_v)**(1.0/2.0)

  wind_term=(-2.0*SIN(xi2)**6 * (COS(xi2)**2+1.0/3.0) + 10.0/63.0) *      &
            2.0*u0*COS(eta_v)**(3.0/2.0)

  rotation_term=((8.0/5.0)*COS(xi2)**3 * (SIN(xi2)**2+2.0/3.0) - pi/4.0) *&
                Earth_radius*omega
                
  eg_baro_t_p = av_temp + const*eta*(wind_term + rotation_term)                
END IF


IF (lhook) CALL dr_hook('EG_BARO_T_P',zhook_out,zhook_handle)
RETURN
END FUNCTION eg_baro_t_p
