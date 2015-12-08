! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

SUBROUTINE hygro_fact(humidity, alpha)

!---------------------------------------------------------------------
! Purpose: To calculate hygroscopic growth factor for sulphate aerosol.
!          Called by SULPHR

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Aerosols

! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP3 v8 programming standards
!---------------------------------------------------------------------


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Arguments with intent IN:
REAL  ::   humidity              ! Relative humidity (0-1)

! Arguments with intent OUT:
REAL  ::   alpha                 ! Fitzgerald's alpha parameter

! Local variables:
REAL  ::   alpha81               ! value of alpha at humidity 0.81

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! This routine calculates growth factors for AMMONIUM SULPHATE
! (based on routine GROW_PARTICLES by D. L. Roberts (1995)).
! Note no growth takes place below a humidity of 0.3.
! Above the deliquescence point, taken as 0.81, the scheme is the
! one due to Fitzgerald (1975).
! This is a simplified version for use with RH < 0.97, assuming that
! BETA=1, and only ALPHA needs to be calculated.

IF (lhook) CALL dr_hook('HYGRO_FACT',zhook_in,zhook_handle)

alpha = 1.0e+00

IF ( humidity  >=  0.3e+00 ) THEN

  IF ( humidity  <   0.81e+00 ) THEN

    ! Calculate ALPHA
    ! We have to be careful here to use the actual value of alpha
    ! at 0.81 otherwise the function would become discontinuous at 0.81.

    alpha81 = afunc1(0.81e+00)
    alpha = 1.0e+00 + (humidity-0.3e+00)                                       &
      *(alpha81-1.0e+00)/0.51e+00

  ELSE IF ( humidity  <=  0.97e+00 ) THEN

    alpha = afunc1(humidity)

  END IF

END IF           ! End RH > 0.3 condition

IF (lhook) CALL dr_hook('HYGRO_FACT',zhook_out,zhook_handle)
RETURN

CONTAINS

REAL FUNCTION afunc1(h)
IMPLICIT NONE
REAL, INTENT(IN) :: h
afunc1 = 1.2e+00*EXP( (0.066e+00*h)/(1.058e+00 - h) )
RETURN
END FUNCTION afunc1

END SUBROUTINE hygro_fact
