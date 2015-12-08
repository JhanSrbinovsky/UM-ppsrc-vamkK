! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Given the number of model steps per hour, the hour and minute
!  and the time between data convert into a fraction between
!  data files for interpolation purposes.
!  NB we assume that the data has one file always at midnight.

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_MAIN.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_getfrac(       &
model_perhour,                    & ! Number of model steps per hour
data_hours,                       & ! Number of hours between data
i_hour,                           & ! Model hour
i_minute,                         & ! Model minute
frac,                             & ! Fraction between two data timesteps
debug)                              ! Miscellanea

USE nudging_control

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

!     Inpts/ outs
INTEGER, INTENT(IN) :: model_perhour      ! Number of model steps per hour
INTEGER, INTENT(IN) :: data_hours         ! Hours between data
INTEGER, INTENT(IN) :: i_hour             ! Model hour
INTEGER, INTENT(IN) :: i_minute           ! Model minute
REAL, INTENT(OUT)   :: frac               ! Return Fraction of ECMWF timestep
                                          ! for interpolation purposes
INTEGER, INTENT(IN) :: debug              ! Debug flag

INTEGER     :: minflag                    ! Flag for timestep in hour
INTEGER     :: hour_offset                ! Accounting tool
INTEGER     :: offset_cut                 ! Accounting tool

INTEGER    :: i                           ! Loop variable
INTEGER    :: total                       ! No. data timesteps per day
INTEGER    :: errstatus                   ! error status

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_GETFRAC',zhook_in,zhook_handle)

!***********************************************************************

! Standard Subroutine Entry Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
 ' NUDGING_GETFRAC: Entered Routine'
END IF

! Flag which timestep we are in in the hour
minflag = INT((i_minute*model_perhour)/60)

!**********************************************************
! Given the hour of the day calculate which strings we need

! Calculate the number of data timesteps in the day
! Subtracting one excludes couting both 24 and 0
! Note we assume that there are an exact integer number
! of steps in a day
total = INT(24/data_hours)  - 1

! Loop over data timesteps in the day
DO i=0, total

! Calculate the time that each of these periods ends
! eg if 4 periods in a day, the first would end at 0600
  offset_cut = data_hours*(i+1)

! If the model hour is less than this time then calculate
! the offset and exit. If not then do nothing
  IF(i_hour < offset_cut) THEN
    hour_offset = data_hours*i
    EXIT
  END IF

END DO !data timesteps

!********************************************************************

! Calculate fraction from components derived above
frac = REAL(i_hour) - REAL(hour_offset)
frac = frac*(REAL(model_perhour)) + REAL(minflag)
frac = frac / (REAL(model_perhour)*REAL(data_hours))

! Drop a line to out the variables used in the above calculation
IF(debug > 10) THEN
  WRITE(OUT, *)                                                        &
   ' NUDGING_GETFRAC: Components of frac equal',                       &
   i_hour, i_minute, hour_offset, model_perhour,                       &
   minflag, data_hours, frac
END IF

! Standard Subroutine Entry Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
 ' NUDGING_GETFRAC: Entered Routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_GETFRAC',zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_getfrac

