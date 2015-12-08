! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
FUNCTION eg_deep_baro_pert(x,y,option, dir)

USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook
USE conversions_mod,     ONLY: pi
USE earth_constants_mod, ONLY: a=>earth_radius

IMPLICIT NONE
!
! Description:
!
! Applies initial perturbation to zonal wind field.
!
! Method:
!  
! Baroclinic instability test case (QJRMS 132, 2943--2975).
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

REAL            :: r, xc,yc,R0,q
REAL            :: x,y
INTEGER         :: option, dir
REAL            :: denom, drdx, drdy, psi0
REAL            :: eg_deep_baro_pert


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_DEEP_BARO_PERT',zhook_in,zhook_handle)

R0 = a / 10.0
xc = pi / 9.0
yc = 2 * pi / 9.0

r = a * ACOS(SIN(yc)*SIN(y) + COS(yc) * COS(y) * COS(x - xc))

SELECT CASE ( option ) 

  CASE ( 1 ) ! JW06 pertubation

    IF ( dir == 1) THEN
      ! u
      eg_deep_baro_pert = EXP(-(r / R0)**2)
    ELSE
      ! v
      eg_deep_baro_pert = 0.0
    END IF

  CASE ( 2 )
    R0 = a/6.0
    psi0 = -R0/2.0
    ! Stream function
    q = 0.0
    IF ( r <= R0 ) THEN
      q  = psi0*EXP( -2.0*r**2/(r - R0)**2 )
      q = -4.0*r/(r-R0)**2*(1.0 - r/(r - R0))*q
    END IF
    denom = 1.0/SQRT(1.0-COS(r/a)**2);

    IF (dir == 1) THEN
    ! u
      drdy = -a*(SIN(yc)*COS(y)-COS(yc)*SIN(y)*COS(x-xc))*denom
      eg_deep_baro_pert = -1.0/a*drdy*q
    ELSE
    ! v
      drdx = a*(COS(yc)*COS(y)*SIN(x-xc))*denom;
      eg_deep_baro_pert = 1.0/(a*COS(y))*drdx*q
    END IF

  CASE ( 3 )
    R0 = a/6.0
!     psi0 = -R0/2.0

    psi0 = -8.0*R0/(3.0*SQRT(3.0)*pi)
    ! Stream function
    q = 0.0
    IF ( r <= R0 ) q  = 0.5*pi*r/R0
! Cos^4 Stream function
    eg_deep_baro_pert = -4.0*psi0*q/r*COS(q)**3*SIN(q)
! Cos^6 Stream function
!     eg_deep_baro_pert = -6.0*psi0*q/r*COS(q)**5*SIN(q)

    denom = 1.0/SQRT(1.0-COS(r/a)**2);

    IF (dir == 1) THEN
    ! u
      drdy = -a*(SIN(yc)*COS(y)-COS(yc)*SIN(y)*COS(x-xc))*denom
      eg_deep_baro_pert =  -1.0/a*eg_deep_baro_pert*drdy 
    ELSEIF (dir == 2) THEN
    ! v
      drdx = a*(COS(yc)*COS(y)*SIN(x-xc))*denom;
      eg_deep_baro_pert =  1.0/(a*COS(y))*eg_deep_baro_pert*drdx
    ELSE
    ! p
      q = 0.5*pi
      IF ( r <= R0 ) q  = 0.5*pi*r/R0
      eg_deep_baro_pert = -psi0*COS(q)**4
    END IF
   
  CASE DEFAULT
    WRITE(0,*) 'Invalid perturbation used in eg_deep_baro_pert'
END SELECT


IF (lhook) CALL dr_hook('EG_DEEP_BARO_PERT',zhook_out,zhook_handle)
RETURN
END FUNCTION eg_deep_baro_pert
