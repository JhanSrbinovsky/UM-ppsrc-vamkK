! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
! One dimensional bi-linear interpolation (i.e 0 extra dimensions)

!  Part of the Nudged model (see nudging_main.F90)

!  Called from various

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_bi_interpolation_zero_d(  &
 input_limit11,                              &  ! Variable at x_0, y_0
 input_limit12,                              &  ! Variable at x_0, y_1
 input_limit21,                              &  ! Variable at x_1, y_0
 input_limit22,                              &  ! Variable at x_1, y_1
 step1,                                      &  ! fraction between x_0 and x_1
 step2,                                      &  ! fraction between y_0 and y_1
 variable,                                   &  ! intrepolated variable
 interp_type,                                &  ! Type of interpolation
 debug)                                         ! Debug flag

USE nudging_control

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook


IMPLICIT NONE

!****************************************************************************

REAL, INTENT(IN)    :: input_limit11
REAL, INTENT(IN)    :: input_limit12
REAL, INTENT(IN)    :: input_limit21
REAL, INTENT(IN)    :: input_limit22
REAL, INTENT(OUT)   :: variable
REAL, INTENT(IN)    :: step1
REAL, INTENT(IN)    :: step2
INTEGER, INTENT(IN) :: interp_type
INTEGER, INTENT(IN) :: debug


!*****************************************
! Variables for Linear Interpolation
REAL :: variable_11
REAL :: variable_21
REAL :: variable_12
REAL :: variable_22

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_BI_INTERPOLATION_ZERO_D',zhook_in,zhook_handle)

!***************************************************************************
! End of Header

! Standard Subroutine Entry Comment
IF(debug > 5) THEN
  WRITE(OUT,*)                                                        &
  'NUDGING_BI_INTERPOLATION_ZERO_D: Entered Routine'
END IF

!************************************
! Linear Interpolation  (At present only have linear interpolation)

! Drop a line to let us know what we are doing
IF(debug > 5) THEN
  WRITE(OUT, *)                                                       &
   ' NUDGING_BI_INTERPOLATION_ZERO_D:',                               &
   ' Using Linear Interpolation'
END IF

! Performing the interpolation
variable_11 =  input_limit11 *(1-step1)*(1-step2)
variable_12 =  input_limit12 *(1-step1)*(step2)
variable_21 =  input_limit21 *(step1)*(1-step2)
variable_22 =  input_limit22 *(step1)*(step2)

variable       = variable_11 + variable_12                            &
               + variable_21 + variable_22

! Drop a line to let us know what we are doing
IF(debug > 15) THEN
  WRITE(OUT, *)                                                       &
    ' NUDGING_BI_INTERPOLATION_ZERO_D: Adding',                       &
    step1, step2,                                                     &
    variable_11, variable_12, variable_21, variable_22,               &
    input_limit11, input_limit12,input_limit21,input_limit22
END IF

!**************************************
! Standard Subroutine Exit Comment
IF(debug > 5) THEN
  WRITE(OUT,*)                                                        &
    ' NUDGING_BI_INTERPOLATION_ZERO_D: Leaving routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_BI_INTERPOLATION_ZERO_D',zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_bi_interpolation_zero_d

