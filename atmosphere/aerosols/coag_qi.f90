! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

REAL FUNCTION coag_qi(a, b, rmedait, rmedacc, saitsq, saccsq)

!---------------------------------------------------------------------
! Purpose: To calculate coagulation coefficients for sulphate aerosol.
!          Called by SULPHR

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols

! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP3 v8 programming standards

!---------------------------------------------------------------------

USE conversions_mod, ONLY: pi
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Arguments with intent IN:
REAL :: a
REAL :: b
REAL :: rmedait
REAL :: rmedacc
REAL :: saitsq
REAL :: saccsq

! Local variables:
REAL :: power1
REAL :: power3
REAL :: power4
REAL :: q1
REAL :: q3
REAL :: q4
REAL :: q5
REAL :: spread

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('COAG_QI',zhook_in,zhook_handle)
power1 = a + b - 2.0
power3 = 3.0*(a - 1.0)
power4 = 3.0*(b - 1.0)

q1 = (1.333333*pi)**power1
q3 = rmedait**power3
q4 = rmedacc**power4

spread = saitsq*(a*a - 1.0) + saccsq*(b*b - 1.0)
spread = spread*4.5
q5 = EXP(spread)

coag_qi = q1*q3*q4*q5

IF (lhook) CALL dr_hook('COAG_QI',zhook_out,zhook_handle)
RETURN
END FUNCTION coag_qi
