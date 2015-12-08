! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Subroutine GAS_CALC ----------------------------------------------

!   Purpose :
!     Calculates the trace gas mixing ratio (or weighting factor for
!   aerosol forcing fields.  Rates of increase (yearly compound factors)
!   can be supplied, or spot values (which will be linearly
!   interpolated) or a mixture of these.  It is designed so it can be
!   called each time step, but when rates of increase are being used,
!   values are in fact only updated at New Year.
!   The rules are:
!     If rates exist (i.e. are positive) for the first & current years
!   then all concentrations are ignored, except for the initial value.
!     If there is a positive rate for the current year but not for the
!   start, the current rate & most recent concentration are used.
!     If rates do not exist for the current year then the concentration
!   is calculated by linear interpolation between the concentrations at
!   the given years.
!     The mixing ratios calculated after the last given year use the
!   rate for the final given year.
!     The 360-day year is assumed.
!   CARE should be taken if solitary rates are specified, as this can
!   result in discontinuities in the concentration time profile at
!   the next given year without a corresponding given rate.

!   ------------------------------------------------------------------

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

      SUBROUTINE gas_calc(gas_now                                       &
                         ,gas_index_max                                 &
                         ,gas_year                                      &
                         ,gas_conc                                      &
                         ,gas_rate                                      &
                         ,max_scenario_pts                              &
                         ,icode)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Control_Max_Sizes
      USE Submodel_Mod
      IMPLICIT NONE

! Maximum size of common block arrays
! Submodel parameters for array sizes

      REAL          gas_now      !OUT Gas concentration at time step
      INTEGER       gas_index_max                                       &
                                 !IN
                  , max_scenario_pts ! IN
      INTEGER       gas_year(max_scenario_pts)   !IN
      REAL          gas_conc(max_scenario_pts)   !IN
      REAL          gas_rate(max_scenario_pts)   !IN
      INTEGER       icode        !OUT Return code: successful=0


! Common blocks

! Model time
! CTIME ----------------------------------------------------
!
!  Purpose: Derived model time/step information including start/end
!           step numbers and frequencies (in steps) of interface field
!           generation, boundary field updating, ancillary field
!           updating; and assimilation start/end times.
!           NB: Last three are set by IN_BOUND, INANCCTL, IN_ACCTL.
!           Also contains current time/date information, current
!           step number (echoed in history file) and steps-per-group.
!
!END -----------------------------------------------------------------

      INTEGER :: I_YEAR               ! Current model time (years)
      INTEGER :: I_MONTH              ! Current model time (months)
      INTEGER :: I_DAY                ! Current model time (days)
      INTEGER :: I_HOUR               ! Current model time (hours)
      INTEGER :: I_MINUTE             ! Current model time (minutes)
      INTEGER :: I_SECOND             ! Current model time (seconds)
      INTEGER :: I_DAY_NUMBER         ! Current model time (day no)
      INTEGER :: PREVIOUS_TIME(7)     ! Model time at previous step
      INTEGER :: IAU_DTResetStep      ! Data time reset step for IAU run

      INTEGER :: BASIS_TIME_DAYS  ! Integral no of days to basis time
      INTEGER :: BASIS_TIME_SECS  ! No of seconds-in-day at basis time

      LOGICAL :: L_C360DY

! UM6.5MODEL_ANALYSIS_HRS changed to REAL - 
!   requires FORECAST_HRS and DATA_MINUS_BASIS_HRS to REAL also 
      REAL    :: FORECAST_HRS     ! Hours since Data Time (ie T+nn)
      REAL    :: DATA_MINUS_BASIS_HRS ! Data time - basis time (hours)

      COMMON /CTIMED/ I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND,    &
        I_DAY_NUMBER,PREVIOUS_TIME,                                     &
        BASIS_TIME_DAYS,BASIS_TIME_SECS,                                &
        FORECAST_HRS,DATA_MINUS_BASIS_HRS,                              &
        IAU_DTResetStep, L_C360DY

      INTEGER :: STEPim(INTERNAL_ID_MAX)  ! Step no since basis time
      INTEGER :: GROUPim(INTERNAL_ID_MAX) ! Number of steps per group

      ! Finish step number this run
      INTEGER :: TARGET_END_STEPim(INTERNAL_ID_MAX)

      REAL :: SECS_PER_STEPim(INTERNAL_ID_MAX) ! Timestep length in secs

      ! Frequency of interface field generation in steps
      INTEGER :: INTERFACE_STEPSim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Start steps for interface field generation
      INTEGER :: INTERFACE_FSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! End steps for interface field generation
      INTEGER :: INTERFACE_LSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Frequency of  updating boundary fields in steps
      INTEGER :: BOUNDARY_STEPSim(INTERNAL_ID_MAX)

      ! No of steps from boundary data prior to basis time to model
      ! basis time
      INTEGER :: BNDARY_OFFSETim(INTERNAL_ID_MAX)

      ! Lowest frequency for updating of ancillary fields in steps
      INTEGER :: ANCILLARY_STEPSim(INTERNAL_ID_MAX)

      ! Start steps for assimilation
      INTEGER :: ASSIM_FIRSTSTEPim(INTERNAL_ID_MAX)

      ! Number of assimilation steps to analysis
      INTEGER :: ASSIM_STEPSim(INTERNAL_ID_MAX)

      ! Number of assimilation steps after analysis
      INTEGER :: ASSIM_EXTRASTEPSim(INTERNAL_ID_MAX)

      COMMON/CTIMEE/                                                    &
     &  STEPim,GROUPim,TARGET_END_STEPim,INTERFACE_STEPSim,             &
     &  INTERFACE_FSTEPim,INTERFACE_LSTEPim,BOUNDARY_STEPSim,           &
     &  BNDARY_OFFSETim,ANCILLARY_STEPSim,ASSIM_FIRSTSTEPim,            &
     &  ASSIM_STEPSim,ASSIM_EXTRASTEPSim,SECS_PER_STEPim

! CTIME end

!     Local variables

      INTEGER       index       ! to subscript gas concs for NOW_TIME
      INTEGER       i           ! Loop over indices
      INTEGER       year_in_secs! Year length in seconds
      INTEGER       now_time_day, now_time_sec
!                               ! Time now in days/secs from time zero
      INTEGER       gas_yr_day1,  gas_yr_sec1
!                               ! Time in days/secs of current GAS_YEAR
      INTEGER       time1
!                               ! The same converted to seconds
      INTEGER       gas_yr_day2,  gas_yr_sec2
!                               ! Time in days/secs of next GAS_YEAR

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!     Check that GASCNST namelist is defined for this year
      IF (lhook) CALL dr_hook('GAS_CALC',zhook_in,zhook_handle)
      IF ( i_year  <   gas_year(1) ) THEN
        icode = 8325
        WRITE (6, *) 'GAS_CALC: no gas data for this year'
        IF (lhook) CALL dr_hook('GAS_CALC',zhook_out,zhook_handle)
        RETURN
      END IF

!     Loop over I to find correct index for current NOW_TIME
      index = 0
      DO i=1, gas_index_max
        IF ( i_year  >=  gas_year(i) ) index = index+1
      END DO

!     Calculate time now in seconds
! DEPENDS ON: time2sec
      CALL time2sec (i_year, i_month, i_day, i_hour, i_minute, i_second,&
                    0, 0, now_time_day, now_time_sec, .TRUE.)

!     If gas rate at current year is non zero calculate new GAS_NOW
!     by considering compound increases of GAS_RATE(1:INDEX)
      IF ( gas_rate(index)  >   0. ) THEN
        year_in_secs = 360 * 86400
! DEPENDS ON: time2sec
        CALL time2sec (gas_year(index), 1, 1, 0, 0, 0,                  &
                      0, 0, gas_yr_day1, gas_yr_sec1, .TRUE.)
        gas_now = gas_conc(1)
        DO i=1, index-1
          IF ( gas_rate(i)  <   0. ) THEN
             gas_now = gas_conc(i+1)
           ELSE
             gas_now = gas_now *                                        &
                  ( gas_rate(i) ** REAL(gas_year(i+1)-gas_year(i)) )
          END IF
        END DO
!       GAS_NOW now holds the concentration in year INDEX - need only
!       update it to the current year.
        gas_now=gas_now*(gas_rate(index)**                              &
          REAL(((now_time_day-gas_yr_day1)*86400+                       &
                now_time_sec-gas_yr_sec1)/year_in_secs))

!     Otherwise calculate by linear interpolation between respective
!     GAS concentrations of given years.
      ELSE
! DEPENDS ON: time2sec
        CALL time2sec (gas_year(index), 1, 1, 0, 0, 0,                  &
                      0, 0, gas_yr_day1, gas_yr_sec1, .TRUE.)
! DEPENDS ON: time2sec
        CALL time2sec (gas_year(index+1), 1, 1, 0, 0, 0,                &
                      0, 0, gas_yr_day2, gas_yr_sec2, .TRUE.)
        time1   = gas_yr_day1*86400 + gas_yr_sec1
        gas_now = gas_conc(index) +                                     &
                ( gas_conc(index+1) - gas_conc(index) )                 &
       * REAL ( now_time_day*86400 + now_time_sec - time1 )             &
            / REAL ( gas_yr_day2*86400 + gas_yr_sec2 - time1 )
      END IF

      IF (lhook) CALL dr_hook('GAS_CALC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE gas_calc
