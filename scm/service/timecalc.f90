! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!     Routine to calculate the year in run and actual time and
!     daynumber in year.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
!---------------------------------------------------------------------
SUBROUTINE timecalc                                                           &
  ( dayno_init, time_init, time_string, lcal360, year, day, time_sec          &
  , previous_time, timehr, timemin )

  USE scm_utils, ONLY:                                                        &
    zhook_in, zhook_out, jprb, lhook, dr_hook

  USE s_main_force, ONLY:                                                     &
    year_init, timestep

  IMPLICIT NONE

  INTEGER ::               &
    dayno_init              ! In daynumber in year of start of of the model

  REAL ::                  &
    time_init               ! In start time in day of the model (seconds).


  LOGICAL ::               &
    lcal360                 ! In true in 360 idealized year

  INTEGER ::               &
    time_sec                ! Out actual time of day in secs. timesteps.

  INTEGER ::               &
    year                   &! Out year
  , day                    &! Out day in year
  , previous_time(7)        ! Out

  CHARACTER(LEN=8) ::           &
    time_string             ! Out actual time in day XX..XX..XX
                            !                        hr..mn..se

! Local Variables
  INTEGER ::               &
    time                   &! Seconds elapsed since model ref time.
  , timesec, timemin       &
  , timehr                 &
  , elapsed_days_prev      &! Days elapsed, end of previous step
  , elapsed_secs_prev      &! Seconds elapsed, end of previous step
  , month, day_m           &
  , elapsed_days           &! Days elapsed since model reference time
  , elapsed_secs           &! Seconds in day
  , basis_time_days        &
  , basis_time_secs        &! Basis time of the model.
  , timestep_count          ! Counter for timesteps

  ! Dr Hook
  !==============================
  REAL(KIND=jprb) :: zhook_handle

  DATA timestep_count /0/   ! Initialise to 1

  IF (lhook) CALL dr_hook('TIMECALC',zhook_in,zhook_handle)

  time_string = '00.00.00'  ! total secs. up to start of this timestep

  time = 24*3600*(dayno_init) + time_init + (timestep_count * timestep)
  elapsed_days_prev = INT(time/(24*3600))
  elapsed_secs_prev = time - (elapsed_days_prev * (24*3600))

  timestep_count = timestep_count + 1

  time = 24*3600*(dayno_init) + time_init + (timestep_count * timestep)
  elapsed_days = INT(time/(24*3600))
  elapsed_secs = time - (elapsed_days * (24*3600))

! DEPENDS ON: time2sec
  CALL time2sec                                                               &
    ( year_init-1, 12, 31, 0, 0, 0, 0, 0, basis_time_days, basis_time_secs    &
    , lcal360 )

! DEPENDS ON: sec2time
  CALL sec2time                                                               &
    ( elapsed_days_prev, elapsed_secs_prev, basis_time_days, basis_time_secs  &
    , previous_time(1), previous_time(2), previous_time(3), previous_time(4)  &
    , previous_time(5), previous_time(6), previous_time(7), lcal360 )

! DEPENDS ON: sec2time
  CALL sec2time                                                               &
    ( elapsed_days, elapsed_secs, basis_time_days, basis_time_secs, year      &
    , month, day_m, timehr, timemin, timesec, day, lcal360 )



! Set up time string for O/P with diagnostics
  WRITE (time_string(1:2), '(i2)') timehr
  WRITE (time_string(4:5), '(i2)') timemin
  WRITE (time_string(7:8), '(i2)') timesec
  time_sec = timehr*3600 + timemin*60 + timesec

  IF (lhook) CALL dr_hook('TIMECALC',zhook_out,zhook_handle)

  RETURN
END SUBROUTINE timecalc

