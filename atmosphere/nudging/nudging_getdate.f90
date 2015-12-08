! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Given date and time produce a string in format
!  YYYYMMDDHH for files required at a given timestep

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_GETFILENAME.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------

SUBROUTINE nudging_getdate(          &
  i_year,                            & ! Model year
  i_month,                           & ! Model month
  i_day,                             & ! Model day
  i_hour,                            & ! Model hour
  date1,                             & ! Return first data date string
  date2,                             & ! Return second data date string
  debug)                               ! Miscellanea

USE nudging_control

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: i_year     ! Model year
INTEGER, INTENT(IN) :: i_month    ! Model month
INTEGER, INTENT(IN) :: i_day      ! Model day
INTEGER, INTENT(IN) :: i_hour     ! Model hour

CHARACTER(LEN=10), INTENT(OUT) :: date1  ! First date string
CHARACTER(LEN=10), INTENT(OUT) :: date2  ! Second date string

INTEGER, INTENT(IN) :: debug       ! Debug flag

! Variables used to assemble the date strings
CHARACTER(LEN=4) :: year1          ! First year string
CHARACTER(LEN=4) :: year2          ! Second year string
CHARACTER(LEN=2) :: month1         ! First month string
CHARACTER(LEN=2) :: month2         ! Second month string
CHARACTER(LEN=2) :: day1           ! First day string
CHARACTER(LEN=2) :: day2           ! Second day string
CHARACTER(LEN=2) :: hour1          ! First hour string
CHARACTER(LEN=2) :: hour2          ! Second hour string

INTEGER :: month_length(12)        ! Length of months (in days)
LOGICAL :: leap                    ! Are we in a leap year?

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_GETDATE',zhook_in,zhook_handle)

!****************************************************************************
!  End of Header

! Standard Entry Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
 ' NUDGING_GETDATE: Entering Routine'
END IF

! Are we in a leap year?
! DEPENDS ON: nudging_getleap
CALL nudging_getleap(         &
     i_year,                  &     ! Year
     leap,                    &     ! Return whether a leap one
     debug)                         ! Miscellanea

! Define arrays of month length
IF (leap) THEN
  month_length=(/                                                      &
   31, 29, 31, 30, 31, 30,                                             &
   31, 31, 30, 31, 30, 31 /)

ELSE
  month_length=(/                                                      &
   31, 28, 31, 30, 31, 30,                                             &
   31, 31, 30, 31, 30, 31 /)

END IF

!********************************************************
! Load hours, days, months and years

IF(i_hour < 18) THEN

! If not in last ECMWF step of the day then always return
! what you put in out
  WRITE(day1,  "(i2.2)") i_day
  WRITE(day2,  "(i2.2)") i_day
  WRITE(month1,"(i2.2)") i_month
  WRITE(month2,"(i2.2)") i_month
  WRITE(year1, "(i4.4)") i_year
  WRITE(year2, "(i4.4)") i_year

!****************************************************
! Evaluate quartetr day divsions
  IF(i_hour < 6) THEN
    WRITE(hour1,  "(i2.2)") 0
    WRITE(hour2,  "(i2.2)") 6
  ELSE IF(i_hour < 12) THEN
    WRITE(hour1,  "(i2.2)") 6
    WRITE(hour2,  "(i2.2)") 12
  ELSE
    WRITE(hour1,  "(i2.2)") 12
    WRITE(hour2,  "(i2.2)") 18
  END IF
!*************************************************

! If in ,last quarter of the day can overun
ELSE
  WRITE(hour1,  "(i2.2)") 18
  WRITE(hour2,  "(i2.2)") 0

! If not last day of the month only day overruns
  IF(i_day < month_length(i_month)) THEN
    WRITE(day1,  "(i2.2)") i_day
    WRITE(day2,  "(i2.2)") (i_day+1)
    WRITE(month1,"(i2.2)") i_month
    WRITE(month2,"(i2.2)") i_month
    WRITE(year1, "(i4.4)") i_year
    WRITE(year2, "(i4.4)") i_year
  ELSE

!  Day overuns to first day of the next month
    WRITE(day1,  "(i2.2)") i_day
    WRITE(day2,  "(i2.2)") 1

! If not last month of teh year only month overruns
    IF(i_month < 12) THEN
      WRITE(month1,"(i2.2)") i_month
      WRITE(month2,"(i2.2)") i_month+1
      WRITE(year1, "(i4.4)") i_year
      WRITE(year2, "(i4.4)") i_year

! Else year overruns
    ELSE
      WRITE(month1,"(i2.2)") i_month
      WRITE(month2,"(i2.2)") 1
      WRITE(year1, "(i4.4)") i_year
      WRITE(year2, "(i4.4)") i_year+1

    END IF ! Last month
  END IF   ! Last day
END IF     ! Last hour

! concatenate together to produce the ECMWF date
date1 = year1 // month1 // day1 // hour1
date2 = year2 // month2 // day2 // hour2

! Standard Exit Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
 ' NUDGING_GETDATE: Leaving Routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_GETDATE',zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_getdate

