! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!

REAL FUNCTION eg_baro_eta_conv(xi1,xi2,z)

USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook
USE atmos_constants_mod, ONLY: R
USE earth_constants_mod, ONLY: g

USE ereport_mod, ONLY : ereport
IMPLICIT NONE
!
! Description:
!
! Converts height coordinate (z) to normalised pressure coordinate (eta)
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
!
! Declarations:

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Function arguments
REAL, INTENT(IN) :: xi1 ! Longitude (radians)
REAL, INTENT(IN) :: xi2 ! Latitude (radians)
REAL, INTENT(IN) :: z   ! Distance above Earth's surface

! Local constants
INTEGER, PARAMETER :: maxits = 1000

! Local variables
LOGICAL :: conv
INTEGER :: itrn, total_its
REAL    :: error
REAL    :: eta_n0, eta_n1, checker
REAL    :: f, f_der,dn

! External functions
REAL, EXTERNAL     :: eg_baro_t_p, eg_baro_geo_p

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_BARO_ETA_CONV',zhook_in,zhook_handle)

conv = .FALSE.
error = 1.0e-14
eta_n0 = 1.0e-7

DO itrn=0,maxits

! DEPENDS ON: eg_baro_geo_p
  f= -g * z + eg_baro_geo_p(xi1,xi2,eta_n0)
! DEPENDS ON: eg_baro_T_p
  f_der = -R / eta_n0 * eg_baro_t_p(xi1,xi2,eta_n0)
  eta_n1 = eta_n0 - f / f_der
  dn = ABS(eta_n1 - eta_n0)
  IF (dn < error) THEN
    conv = .TRUE.
    total_its = itrn
    EXIT
  END IF
  eta_n0 = eta_n1
END DO

IF (conv) THEN
  eg_baro_eta_conv = eta_n1
ELSE

     Call ereport("eg_baro_eta_conv", 1,                              &
                  "did not converge" )
END IF

IF (lhook) CALL dr_hook('EG_BARO_ETA_CONV',zhook_out,zhook_handle)
RETURN
END FUNCTION eg_baro_eta_conv
