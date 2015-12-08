! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Subroutine LLTOEQ-------------------------------------------------
!
!    Purpose:  Calculates latitude and longitude on equatorial
!              latitude-longitude (eq) grid used in regional
!              models from input arrays of latitude and
!              longitude on standard grid. Both input and output
!              latitudes and longitudes are in degrees.
!
!
!   Programming standard :
!
!   Logical components covered : S132
!
!   Project task :
!
!    Documentation: The transformation formulae are described in
!                   unified model on-line documentation paper S1.
!
!                   Code Owner: See Unified Model Code Owners HTML page
!                   This file belongs in section: Grids
MODULE lltoeq_mod

IMPLICIT NONE

CONTAINS


SUBROUTINE lltoeq                                                 &
(phi,lambda,phi_eq,lambda_eq,phi_pole,lambda_pole,points)


USE conversions_mod, ONLY: pi_over_180, recip_pi_over_180
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN) ::  points   !IN  Number of points to be processed

REAL, INTENT(IN) ::     phi(points)      !IN  Latitude
REAL, INTENT(IN) ::     lambda(points)   !IN  Longitude

REAL, INTENT(IN) ::     phi_pole !IN  Latitude of equatorial lat-lon pole

REAL, INTENT(OUT) ::     lambda_eq(points)
                                 !Longitude in equatorial lat-lon coords
REAL, INTENT(OUT) ::     phi_eq(points)
                                 !Latitude in equatorial lat-lon coords
                                 
REAL, INTENT(IN) ::   lambda_pole   !INOUT  Longitude of equatorial lat-lon pole                              



! Define local varables
REAL  ::  a_lambda,a_phi,e_lambda,arg,e_phi,sin_phi_pole,cos_phi_pole
REAL  ::  term1,term2,lambda_zero,lambda_pole_local
REAL, PARAMETER :: small=1.0e-6
INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle



!  1. Initialise local constants
! Scale lambda pole to range -180 to 180 degs
IF (lhook) CALL dr_hook('LLTOEQ',zhook_in,zhook_handle)

IF (lambda_pole >   180.0) THEN
  lambda_pole_local=lambda_pole-360.0
ELSE
  lambda_pole_local=lambda_pole
END IF

! Latitude of zeroth meridian
lambda_zero=lambda_pole_local+180.
! Sine and cosine of latitude of eq pole
IF (phi_pole >= 0.0) THEN
  sin_phi_pole =  SIN(pi_over_180*phi_pole)
  cos_phi_pole =  COS(pi_over_180*phi_pole)
ELSE
  sin_phi_pole = -SIN(pi_over_180*phi_pole)
  cos_phi_pole = -COS(pi_over_180*phi_pole)
END IF

!  2. Transform from standard to equatorial latitude-longitude

DO i= 1,points

! Scale longitude to range -180 to +180 degs

  a_lambda=lambda(i)-lambda_zero
  IF(a_lambda >   180.0)a_lambda=a_lambda-360.
  IF(a_lambda <= -180.0)a_lambda=a_lambda+360.

! Convert latitude & longitude to radians

  a_lambda=pi_over_180*a_lambda
  a_phi=pi_over_180*phi(i)

! Compute eq latitude using equation (4.4)

  arg=-cos_phi_pole*COS(a_lambda)*COS(a_phi)                      &
                     +SIN(a_phi)*sin_phi_pole
  arg=MIN(arg, 1.0)
  arg=MAX(arg,-1.0)
  e_phi=ASIN(arg)
  phi_eq(i)=recip_pi_over_180*e_phi

! Compute eq longitude using equation (4.6)

  term1 =(COS(a_phi)*COS(a_lambda)*sin_phi_pole                   &
         +SIN(a_phi)*cos_phi_pole)
  term2=COS(e_phi)
  IF(term2 <  small) THEN
    e_lambda=0.0
  ELSE
    arg=term1/term2
    arg=MIN(arg, 1.0)
    arg=MAX(arg,-1.0)
    e_lambda=recip_pi_over_180*ACOS(arg)
    e_lambda=SIGN(e_lambda,a_lambda)
  END IF

! Scale longitude to range 0 to 360 degs

  IF(e_lambda >= 360.0) e_lambda=e_lambda-360.0
  IF(e_lambda <  0.0) e_lambda=e_lambda+360.0
  lambda_eq(i)=e_lambda

END DO

IF (lhook) CALL dr_hook('LLTOEQ',zhook_out,zhook_handle)
RETURN
END SUBROUTINE lltoeq
END MODULE lltoeq_mod
