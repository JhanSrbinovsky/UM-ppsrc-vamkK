! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Given a year, is it a leap year?

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_GETDATE.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_getleap(            &
i_year,                                & ! Year
leap,                                  & ! output leap year or not
debug)                                   ! Miscellanea

USE nudging_control

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN)  :: i_year           ! Year
LOGICAL, INTENT(OUT) :: leap             ! string to return

INTEGER, INTENT(IN)  :: debug            ! Debug flag

! Remainders when we divide by 4, 100 & 400
INTEGER :: rem1
INTEGER :: rem2
INTEGER :: rem3

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_GETLEAP',zhook_in,zhook_handle)

!**********************************************************

! Standard Subroutine Entry Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                        &
 ' NUDGING_GETLEAP: Entered Routine'
END IF

! Calculate remainders when we divide by 4, 100 & 400
rem1 = MOD(i_year, 4)
rem2 = MOD(i_year, 100)
rem3 = MOD(i_year, 400)

! Leap year if divisible by 4
! Unless divisible by 100 when not
! Unless divisible by 400 when is
IF(rem1 == 0)  THEN
  IF(rem2 == 0) THEN
    IF(rem3 == 0) THEN
      leap = .TRUE.
    ELSE
      leap = .FALSE.
    END IF
  ELSE
    leap = .TRUE.
  END IF
ELSE
  leap = .FALSE.
END IF

! Standard Subroutine Exit Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
 ' NUDGING_GETLEAP: Leaving Routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_GETLEAP',zhook_out,zhook_handle)

RETURN

END SUBROUTINE nudging_getleap
