! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!   Inetrpolates a 4D variable (ie 1D + 3D)

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_MAIN.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
     SUBROUTINE nudging_interpolation_three_d(    &
 input_dim1,                                & ! extra dimension 1
 input_dim2,                                & ! extra dimension 2
 input_dim3,                                & ! extra dimension 3
 input_limit1,                              & ! starting variable
 input_limit2,                              & ! finishing variable
 step,                                      & ! Fraction between start and end
 variable,                                  & ! Interpolated variable
 interp_type,                               & ! Interp type
 debug)                                       ! Debug flag

USE nudging_control

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

!*********************************************************************************

INTEGER, INTENT(IN)   :: input_dim1  ! extra dimension 1
INTEGER, INTENT(IN)   :: input_dim2  ! extra dimension 2
INTEGER, INTENT(IN)   :: input_dim3  ! extra dimension 3

! Variable at start at end of interpolation
REAL, INTENT(IN)      :: input_limit1                           &
                 (1: input_dim1, 1: input_dim2, 1: input_dim3)
REAL, INTENT(IN)      :: input_limit2                           &
                 (1: input_dim1, 1: input_dim2, 1: input_dim3)
! Interpolated variable
REAL, INTENT(OUT)     :: variable                               &
                 (1: input_dim1, 1: input_dim2, 1: input_dim3)

REAL, INTENT(IN)      :: step        ! Fraction between steps
INTEGER, INTENT(IN)   :: interp_type ! Interpolation type
INTEGER, INTENT(IN)   :: debug       ! Debug Flag


! Slope for linear interpolation
REAL                  :: delta_variable                         &
                  (1: input_dim1, 1: input_dim2, 1: input_dim3)
INTEGER               :: i, j, k     ! Looping variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_INTERPOLATION_THREE_D',zhook_in,zhook_handle)

!****************************************************************************
! End of Header

! Standard Subroutine Entry Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
 ': NUDGING_INTERPOLATION_THREE_D: Entered Routine'
END IF

!************************************
! Linear Interpolation  (At present only have linear interpolation)

! Drop a line to let us know what we are doing
IF(debug > 10) THEN
  WRITE(OUT, *)                                                        &
   ': NUDGING_INTERPOLATION_3D: Using Linear Interpolation',           &
   ' Dimensions equal ',                                               &
   input_dim1, input_dim2, input_dim3, step
END IF

! Loop over the arrays performing the interpolation
DO i=1,input_dim1
  DO j=1,input_dim2
    DO k=1, input_dim3

      delta_variable(i,j,k) = input_limit2(i,j,k)                      &
                              - input_limit1(i,j,k)
      variable(i,j,k)       = input_limit1(i,j,k)                      &
                              + delta_variable(i,j,k)                  &
                              * step
    END DO
  END DO
END DO ! 3D Linear Interpolation

!**************************************
! Standard Subroutine Exit Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                        &
   ' NUDGING_INTERPOLATION_THREE_D: Leaving routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_INTERPOLATION_THREE_D',zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_interpolation_three_d

