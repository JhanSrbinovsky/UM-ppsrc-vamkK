! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Subroutine W_COEFF------------------------------------------------
!
!    Purpose:
!            Calculates coefficients used to translate u and v compo-
!            nents of wind between equatorial (eq) latitude-longitude
!            grid and standard latitude-longitude grid (or visa versa).
!            Input latitudes and longitudes are in degrees.
!
!
!    Documentation: The transformation formulae are described in
!                   unified model on-line documentation paper S1.
!
!    ------------------------------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Grids
MODULE w_coeff_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE w_coeff(                                               &
 coeff1,coeff2,lambda,lambda_eq,phi_pole,lambda_pole,points)

USE conversions_mod, ONLY: pi_over_180
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN) ::  points  !IN  Number of points to be processed

REAL, INTENT(OUT)   ::  coeff1(points) !OUT Coefficient of rotation no 1
REAL, INTENT(OUT)   ::  coeff2(points) !OUT Coefficient of rotation no 2

REAL, INTENT(IN)    ::  lambda(points) !IN  Longitude
REAL, INTENT(IN)    ::  lambda_eq(points)
                            !IN  Longitude in equatorial lat-lon coords
REAL, INTENT(IN)    ::   phi_pole  !IN  Latitude of equatorial lat-lon pole
REAL, INTENT(IN)    ::   lambda_pole  !IN  Longitude of equatorial lat-lon pole

! Define local varables:-----------------------------------------------
REAL  ::  a_lambda,e_lambda,sin_e_lambda,sin_phi_pole 
REAL  ::  cos_phi_pole,c1,c2,lambda_zero
INTEGER  ::  i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------

!  1. Initialise local constants

! Longitude of zeroth meridian
IF (lhook) CALL dr_hook('W_COEFF',zhook_in,zhook_handle)
lambda_zero=lambda_pole+180.

! Sine and cosine of latitude of eq pole

sin_phi_pole=SIN(pi_over_180*phi_pole)
cos_phi_pole=COS(pi_over_180*phi_pole)

!  2. Evaluate translation coefficients
DO i=1,points

! Actual longitude converted to radians

  a_lambda=pi_over_180*(lambda(i)-lambda_zero)

! Convert eq longitude to radians and take sine

  e_lambda=lambda_eq(i)*pi_over_180
  sin_e_lambda=SIN(e_lambda)

! Formulae used are from eqs (4.19) and (4.21)

  c1=SIN(a_lambda)*sin_e_lambda*sin_phi_pole                      &
             +COS(a_lambda)*COS(e_lambda)
  coeff1(i)=c1




! avoid rounding error problems
  IF (ABS(c1) >= 1.0) THEN
    c1 = 1.0
    c2 = 0.0
  ELSE
    c2=SQRT(1.0-c1*c1)
  END IF


! Set the sign of C2. See ticket #796 for an explanation.
  coeff2(i)=SIGN(c2,sin_e_lambda*cos_phi_pole)

END DO

IF (lhook) CALL dr_hook('W_COEFF',zhook_out,zhook_handle)
RETURN
END SUBROUTINE w_coeff
END MODULE w_coeff_mod
