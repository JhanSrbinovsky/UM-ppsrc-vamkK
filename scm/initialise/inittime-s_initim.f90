! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    INITTIME_SCMA
!
!    Purpose: Initialises the model time relative to the calender zero
!             time. The model basis time (MBT) is converted to a time
!             in seconds since T=0 with respect to the calender. The
!             model initialisation time in days and seconds is then
!             calculated relative to this.
!
!    INPUTS - Year
!           - Month
!           - Day
!           - Hour
!           - Minute
!           - Second
!           - Flag for 360 day year calendar
!
!    OUTPUTS - Daynumber in year
!            - Time within daynumber
!
!
!   Phil Hopwood <- programmer of some or all of previous code changes
!
!
!
!   Documentation: SSFM Project Documentation - Implementation of UM
!                  date structure within the SCM
!
!  -------------------------------------------------------------------

!   Arguments --------------------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Single Column Model

SUBROUTINE inittime_scm                                                      &
  ( dayno, daytime, lcal360 )

  USE scm_utils, ONLY:                                                        &
    zhook_in, zhook_out, jprb, lhook, dr_hook

  USE s_main_force, ONLY:                                                     &
    year_init, month_init, day_init, hour_init, min_init, sec_init

  IMPLICIT NONE

  EXTERNAL time2sec

  LOGICAL ::         &
    lcal360           ! IN Flag for 360 year calender

  INTEGER ::         &
    dayno            &! OUT Daynumber in year
  , daytime           ! OUT Time in daynumber

  INTEGER ::         &
    basis_time_days  &! LOCAL whole days to basis time from start of calender
  , basis_time_secs   ! LOCAL secs in day at basis time

  ! Dr Hook
  !==============================
  REAL(KIND=jprb) :: zhook_handle

  IF (lhook) CALL dr_hook('INITTIME_SCM',zhook_in,zhook_handle)

! Calculate Model Basis Time (i.e. beginning of current year)
! relative to calender zero

! DEPENDS ON: time2sec
  CALL time2sec                                                               &
    ( year_init-1, 12, 31, 0, 0, 0, 0, 0, basis_time_days, basis_time_secs    &
    , lcal360 )

! Calculate daynumber and initial time relative to Model Basis Time

! DEPENDS ON: time2sec
  CALL time2sec                                                               &
    ( year_init, month_init, day_init, hour_init, min_init, sec_init          &
    , basis_time_days, basis_time_secs, dayno, daytime, lcal360 )

  IF (lhook) CALL dr_hook('INITTIME_SCM',zhook_out,zhook_handle)

  RETURN
END SUBROUTINE inittime_scm

