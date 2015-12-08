! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
FUNCTION eg_baro_pert(x,y)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE conversions_mod

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

REAL, PARAMETER :: a = 6.371229e6
REAL, PARAMETER :: up = 1.0
REAL            :: r, xc,yc,rg
REAL            :: x,y
REAL            :: eg_baro_pert


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_BARO_PERT',zhook_in,zhook_handle)

r = a / 10.0
xc = pi / 9.0
yc = 2 * pi / 9.0
rg = a * ACOS(SIN(yc)*SIN(y) + COS(yc) * COS(y) * COS(x - xc))

eg_baro_pert = up * EXP(-(rg / r)**2)

IF (lhook) CALL dr_hook('EG_BARO_PERT',zhook_out,zhook_handle)
RETURN
END FUNCTION eg_baro_pert
