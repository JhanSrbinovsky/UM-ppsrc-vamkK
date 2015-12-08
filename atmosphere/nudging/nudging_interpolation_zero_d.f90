! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
! One dimensional interpolation (i.e 0 extra dimensions)

!  Part of the Nudged model (see nudging_main.F90)

!  Called from various

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_interpolation_zero_d(  &
  input_limit1,                           & ! Input variable at the start
  input_limit2,                           & ! Input variable at the end
  step,                                   & ! Fraction between start and end
  variable,                               & ! Variable at this step
  interp_type,                            & ! Type of interpolation (eg linear =0)
  debug)                                    ! Debug flag

USE nudging_control                        ! Use standard nudging switches

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

REAL, INTENT(IN)      :: input_limit1      ! Variable at lower limit
REAL, INTENT(IN)      :: input_limit2      ! Variable at upper limit
REAL, INTENT(OUT)     :: variable          ! Interpolated variable
REAL, INTENT(IN)      :: step              ! Fraction between limits
INTEGER, INTENT(IN)   :: interp_type       ! Type of interpolation (eg linear)

! Miscellanous variables
INTEGER, INTENT(IN)   :: debug             ! debug flag

REAL                  :: delta_variable    ! Slope for linear interp.

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_INTERPOLATION_ZERO_D',zhook_in,zhook_handle)

!************************************************************************
! End of Header

!************************************
! Linear Interpolation  (At present only have linear interpolation)

delta_variable = input_limit2 - input_limit1
variable       = input_limit1 + delta_variable*step

IF (lhook) CALL dr_hook('NUDGING_INTERPOLATION_ZERO_D',zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_interpolation_zero_d

