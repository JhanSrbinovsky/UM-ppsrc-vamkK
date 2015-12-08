! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Routine: GET_NAME -------------------------------------------------
!
!    Purpose: Generates an output file name of up to 20 characters using
!             the defined file naming convention, taking account of
!             file type, validity time, etc., then adds the full path
!             of the Rose cycling task directory to it.
!             Obeys new filenaming convention introduced at version 2.7.
!
!    Programming standard: UM Doc Paper 3, version 8.2 (25/3/2009)
!
!    Project task: S51
!
!    External documentation: UM documentation paper 7 - Filenaming
!                            conventions for the Unified Model
!
!    -------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Top Level

SUBROUTINE get_name(expt_id,job_id,isubmodel,meanlev,toggle,      &
     reinit_steps,filetype,letter_3, model_status,                &
     time_convention,analysis_hrs,filename,icode,cmessage,        &
     lcal360)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE filenamelength_mod, ONLY : filenamelength
USE Control_Max_Sizes
USE nlstgen_mod, ONLY : steps_per_periodim, secs_per_periodim,    &
                        meanfreqim, dumpfreqim
USE Submodel_Mod
IMPLICIT NONE
LOGICAL lcal360

CHARACTER(LEN=4)   expt_id     ! IN  - Experiment ident or alias
CHARACTER(LEN=1)   job_id      ! IN  - Job ident within experiment
INTEGER       isubmodel   ! IN  - Submodel indicator
INTEGER       meanlev     ! IN  - Mean level indicator
INTEGER       toggle      ! IN  - Alternately 1/2 for partial sums
REAL          analysis_hrs! IN  - Hrs from basis time to analysis
                          !       Allow for fractional ANALYSIS_HRS
                          !       Changed from INTEGER to REAL
INTEGER       reinit_steps! IN  - timesteps between file reinitial-
!                                       isation for non-mean pp files
!                                       or -ve for Gregorian reinit.

CHARACTER(LEN=1)   filetype    ! IN  - Code for file type
CHARACTER(LEN=1)   letter_3    ! IN  - character for use in position 9
!                                       of non-mean pp files.
CHARACTER(LEN=14)  model_status! IN  - Operational/NonOperational

CHARACTER(LEN=17)  time_convention ! IN  - Relative/Timestep/
!                                     Absolute_standard/Absolute_long/
!                                     Absolute_short/Sub-hourly
CHARACTER(LEN=filenamelength) :: filename  ! OUT - Generated file name
INTEGER icode             ! OUT - Error return code
CHARACTER(LEN=80) cmessage

!*----------------------------------------------------------------------
!  Common blocks
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

!  Local variables

INTEGER                                                           &
 yyyy,mm,dd,hh,imin,isec   ! Current time values for filename
INTEGER     count,                                                &
                           ! Counter for steps or hours
            dayno,                                                &
                           ! day number
            days,                                                 &
                           ! Number of days for period
            hours,                                                &
                           ! Number of hours for period
            i,                                                    &
                           ! loop counter
            steps          ! number of steps
INTEGER     end_days       ! number of whole days from run start
INTEGER     end_seconds    ! number of extra secs from run start
INTEGER     mon            ! month for mean period
INTEGER     a_steps_per_hr ! steps per hour for atmos sub-model
INTEGER     whole_day_flag ! in period for reinitialisation

CHARACTER(LEN=2) qw             ! Operational file prefix

CHARACTER(LEN=1) submodel_id    ! character for model a or o
CHARACTER(LEN=1) filetype_2     ! letter after FILETYPE in name
CHARACTER(LEN=1) mean_period(4) ! default letter for mean period

CHARACTER(LEN=1) y_hunds        ! Character year identifier (hundreds)
CHARACTER(LEN=1) y_tens         ! Character year identifier (tens)
CHARACTER(LEN=1) y_units        ! Character year identifier (units)
CHARACTER(LEN=1) m              ! Character month identifier
CHARACTER(LEN=1) d              ! Character day-of-month identifier
CHARACTER(LEN=1) h              ! Character hour identifier
CHARACTER(LEN=1) hundreds       ! Character for hundreds counter
CHARACTER(LEN=1) tens           ! Character for tens counter
CHARACTER(LEN=1) units          ! Character for units counter
CHARACTER(LEN=1) deci           ! Character for tenths (=mins)
CHARACTER(LEN=1) char_id(36)    ! Valid characters for above (lookup)
CHARACTER(LEN=1) separator      ! character used as separator in name
CHARACTER(LEN=1) style          ! style of date in filename
CHARACTER(LEN=3) cdayno         ! character day number
CHARACTER(LEN=3) month_3char(12)! 3 character month identifier
CHARACTER(LEN=2) month_2char(12)! 2 character month identifier
CHARACTER(LEN=3) season_3char(12)! 3 character season identifier
CHARACTER(LEN=2) season_2char(12)! 2 character season identifier

CHARACTER(LEN=filenamelength) :: datam ! Path to dir containing filename

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

DATA qw / 'qw'/
DATA mean_period / '1', '2', '3', '4' /
DATA char_id/ '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',   &
              'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',   &
              'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',   &
              'u', 'v', 'w', 'x', 'y', 'z' /
DATA month_3char/ 'jan', 'feb', 'mar', 'apr', 'may', 'jun',       &
                  'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/
DATA month_2char/ 'ja', 'fb', 'mr', 'ar', 'my', 'jn',             &
                  'jl', 'ag', 'sp', 'ot', 'nv', 'dc'/
DATA season_3char/ 'ndj', 'djf', 'jfm', 'fma', 'mam', 'amj',      &
                   'mjj', 'jja', 'jas', 'aso', 'son', 'ond'/
DATA season_2char/ 'nj', 'df', 'jm', 'fa', 'mm', 'aj',            &
                   'mj', 'ja', 'js', 'ao', 'sn', 'od'/
!
! ----------------------------------------------------------------------
!  1. Determine submodel id - (used in filechar 6 or 7 if operational)
!
IF (lhook) CALL dr_hook('GET_NAME',zhook_in,zhook_handle)

IF (isubmodel == atmos_sm) THEN
  submodel_id= 'a'
ELSE
  icode=2
  cmessage='GET_NAME: Illegal sub-model specified'
  GO TO 9999
END IF

! 1.1 Compute steps per hour for atmosphere sub_model
a_steps_per_hr = 3600*steps_per_periodim(a_im)/                   &
                       secs_per_periodim(a_im)

! ----------------------------------------------------------------------
!  2. Determine style filename and separator
!
IF (filetype /= 's') THEN

!
!  2.1 Relative time convention
!
  IF (time_convention == 'Relative        ') THEN
    separator='_'
    style='A'
    count = stepim(a_im) / a_steps_per_hr - analysis_hrs
!
!  2.2 Step time convention
!
  ELSE IF (time_convention == 'Timestep         ') THEN
    separator='_'
    style='A'
    count = stepim(a_im)
!
!  2.3 Absolute time convention -standard version
!
  ELSE IF (time_convention == 'Absolute_standard') THEN
    separator='.'
    style='B'
!
!  2.4 Absolute time convention - short
!
  ELSE IF (time_convention == 'Absolute_short   ') THEN
    separator='-'
    style='B'
!
!  2.5 Absolute time convention - long
!
  ELSE IF (time_convention == 'Absolute_long    ') THEN
    separator='@'
    style='B'
!
!  2.5.1 Absolute time convention - Date Stamp
!    
  ELSE IF (time_convention == 'Absolute_Dstamp  ') THEN
    separator='.'
    style='D'
!
!  2.6 Sub-hourly filenaming time convention
!
  ELSE IF (time_convention == 'Sub-hourly      ') THEN
    separator='_'
    style='A'
    count = stepim(a_im) * secs_per_stepim(a_im)
    WRITE(6,*)'GET_NAM - count ',a_im,count,stepim(a_im),         &
    secs_per_stepim(a_im)

!   ANALYSIS_HRS as REAL needs careful coding of sub-hourly file-naming.
    count = 100*(REAL(count)/3600.0 - analysis_hrs)
    count = (count/100)*100 +                                     &
             NINT((REAL(MOD(count,100))/100.0)*60)

    WRITE(6,*)'GET_NAM - hhmm ',count,analysis_hrs
!
  ELSE
    icode=1
    cmessage='GET_NAME: Illegal TIME_CONVENTION specified'
    GO TO 9999
  END IF
! ----------------------------------------------------------------------
!
!  3.0 work out encoding of date time and filetype_2
!
  IF (style == 'A') THEN
    IF (count <  0) THEN
      count = -count
      filetype_2='z'
      IF (time_convention == 'Sub-hourly      ') THEN
!!!            FILETYPE_2='h'
        IF (count  >=  36000) THEN
          icode = count
          cmessage='GET_NAME: count too big for '//               &
                   'sub-hourly filenaming convention'
          GO TO 9999
        END IF
        hundreds = char_id(mod(count/1000,36)+1)
        tens =     char_id(mod(count/100, 10)+1)
        units =    char_id(mod(count/10,  10)+1)
        deci =     char_id(mod(count,     10)+1)
      ELSE
        IF (count  >=  3600) THEN
          icode = count
          cmessage='GET_NAME: count too big for '//               &
                   'hourly or timestep filenaming convention'
          GO TO 9999
        END IF
        hundreds = char_id(mod(count/100,36)+1)
        tens =     char_id(mod(count/10 ,10)+1)
        units =    char_id(mod(count,    10)+1)
      END IF
    ELSE
      IF (filetype == 'p') THEN
        filetype_2=letter_3
      ELSE IF (filetype == 'b') THEN   ! Boundary File
        filetype_2=letter_3
      ELSE IF (filetype == 'c') THEN   ! Macro File
        filetype_2=letter_3
      ELSE
        filetype_2='a'
      END IF
      IF (time_convention == 'Sub-hourly      ') THEN
!!!            FILETYPE_2='h'
        IF (count  >=  36000) THEN
          icode = count
          cmessage='GET_NAME: count too big for '//               &
                   'sub-hourly filenaming convention'
          GO TO 9999
        END IF
        hundreds = char_id(mod(count/1000,36)+1)
        tens =     char_id(mod(count/100, 10)+1)
        units =    char_id(mod(count/10,  10)+1)
        deci =     char_id(mod(count,     10)+1)
      ELSE
        IF (count  >=  3600) THEN
          icode = count
          cmessage='GET_NAME: count too big for '//                &
                   'hourly or timestep filenaming convention'
          GO TO 9999
        END IF
        hundreds = char_id(mod(count/100,36)+1)
        tens =     char_id(mod(count/10 ,10)+1)
        units =    char_id(mod(count,    10)+1)
      END IF
    END IF
  ELSE IF (style == 'B' .OR. style == 'D') THEN   
! Some sort of absolute time, C transformed from style B
! D is a new, E transformed from style D 
    whole_day_flag = 0

! Current date time is

    yyyy  = i_year
    mm    = i_month
    dd    = i_day
    hh    = i_hour
    dayno = i_day_number

! Instantaneous files

    IF (meanlev == 0) THEN
    
      IF (filetype == 'd') THEN
        filetype_2 = 'a'      
      END IF

! Work out reintialisation period for pp and boundary files.
! Note: assumes reinitialisation period is whole number of hours.
! This is not strictly true but is probably ok for this purpose.

      IF (filetype == 'p'                                          &
                           !  PP File
     .OR. filetype == 'b'                                          &
                           !  Boundary File
     .OR. filetype == 'c'                                          &
                           !  Macro File
    ) THEN

        filetype_2=letter_3      
    
        hours = reinit_steps/a_steps_per_hr
        IF (reinit_steps <  0) THEN ! Gregorian reinitialisation
          hours=720 ! dummy: could be anything divisible by 24
        END IF
!   do further checks if multiple of 1 day
        IF (MOD(hours,24) == 0) THEN  ! whole days in period, or
          days=hours/24               ! Gregorian reinit.
          whole_day_flag = 1

          IF (dd == 1 .AND. time_convention /= 'Absolute_short   '&
            .AND. style == 'B'                                    &
            .AND. (days == 30 .OR. reinit_steps <  0)) THEN
!  Original code didn't allow style=C for 3-month reinit files but new
!  code (Gregorian reinit) does, at least in section 3.0.
            style='C'                ! month in characters
          END IF

          IF (dd == 1 .AND. time_convention == 'Absolute_Dstamp  '&
         .AND. (days == 30 .OR. reinit_steps <  0)) THEN
!  Style=E is identical to style C but derived from the new style  
!  D and, has yyyy stamp and month or season in 3 characters
            style='E'                ! month in characters
          END IF          
        END IF

!  For instantaneous pp file which does not use new naming style
!  with DateStamp we need to work out the end of period as call to
!  this routine occurrs on the first output timestep.

        IF (filetype == 'p' .AND.                                 &
            time_convention /= 'Absolute_Dstamp  ' ) THEN
          filetype_2 = letter_3
          IF (style /= 'C') THEN
            IF (stepim(a_im) == 0) THEN
              steps = reinit_steps
            ELSE
              steps = stepim(a_im) + reinit_steps - 1
            END IF
! DEPENDS ON: stp2time
            CALL stp2time(steps,                                  &
                          a_steps_per_hr*24,86400,                &
                          end_days,end_seconds)
! DEPENDS ON: sec2time
            CALL sec2time(end_days,end_seconds,                   &
                          basis_time_days,basis_time_secs,        &
                          yyyy,mm,dd,hh,imin,isec,dayno,lcal360)
          END IF
        END IF        
      END IF  ! End of block for b,c d,pp files   
  
    ELSE         !  MEANS

!  determine if special mean period
      hours=dumpfreqim(a_im)/a_steps_per_hr
      DO i=1,meanlev
        hours=hours*meanfreqim(i,a_im) !hours per meaning period
      END DO
      IF (mod(hours,24) == 0) THEN
        days=hours/24
        whole_day_flag = 1  
! DEPENDS ON: day2char
        CALL day2char(days,filetype_2)
        IF (filetype_2 == '0') THEN
          filetype_2=mean_period(meanlev)
        ELSE IF (filetype_2 == 'm'.AND.dd == 1                    &
       .AND.time_convention /= 'Absolute_short    ') THEN
          IF (style == 'B') THEN
            style='C'      ! period starts at beginning of a month
          ELSE IF (style == 'D') THEN 
            style='E'      ! period starts at beginning of a month  
          END IF          
          IF (mm == 1) THEN ! correct year if month december
            yyyy=yyyy-1
          END IF
        ELSE IF (filetype_2 == 's'.AND.dd == 1                    &
       .AND.time_convention /= 'Absolute_short    ') THEN
          IF (style == 'B') THEN       
            style='C'      ! period starts at beginning of a season
          ELSE IF (style == 'D')  THEN 
            style='E'      ! period starts at beginning of a season  
          END IF                      
        END IF
      ELSE
        filetype_2=mean_period(meanlev)
      END IF
      
    END IF

    y_units = char_id(mod(yyyy,10)+1)
    m = char_id(mm+1)
    d = char_id(dd+1)
    h = char_id(hh+1)

  END IF

ELSE
! partial sum files - no date time information required
  separator='_'
  filetype_2 = mean_period(meanlev)
END IF

! ----------------------------------------------------------------------
!  3.1 Construct filename from the various components
!
! Generated filenames are up to 20 characters in length. 
! The remainder of the string is used to add a directory path later.
filename="                    "
IF (model_status == 'Operational'.OR.model_status == 'SCS') THEN
  filename(1:2)  =qw
  filename(3:6)  =expt_id
  filename(7:7)  =submodel_id
  filename(8:8)  =separator
  filename(9:9)  =filetype
  filename(10:10)=filetype_2
  filename(11:11)=hundreds
  filename(12:12)=tens
  filename(13:13)=units
  IF (time_convention == 'Sub-hourly      ') THEN
    filename(14:14)=deci
  END IF
ELSE
  filename(1:4)  =expt_id
  filename(5:5)  =job_id
  filename(6:6)  =submodel_id
  filename(7:7)  =separator
  filename(8:8)  =filetype
  filename(9:9)  =filetype_2
  IF (filetype == 's') THEN
    IF (toggle == 1) THEN
      filename(10:10)='a'
    ELSE
      filename(10:10)='b'
    END IF
  ELSE IF (style == 'A') THEN
    filename(10:10)=hundreds
    filename(11:11)=tens
    filename(12:12)=units
    IF (time_convention == 'Sub-hourly      ') THEN
      filename(13:13)=deci
    END IF
  ELSE IF (style == 'B') THEN
    IF (time_convention == 'Absolute_standard') THEN

! decades measured relative to 1800
      y_tens  = char_id(mod(yyyy/10,36)+1)
      filename(10:10)=y_tens
      filename(11:11)=y_units
      filename(12:12)=m
      filename(13:13)=d
      filename(14:14)=h
      
    ELSE IF (time_convention == 'Absolute_long    ') THEN

! centuries measured from 0 i.e. 1992 as j92, with wraparound at 3600
      y_hunds = char_id(mod(yyyy/100,36)+1)
      y_tens  = char_id((yyyy-(yyyy/100)*100)/10+1)
      filename(10:10)=y_hunds
      filename(11:11)=y_tens
      filename(12:12)=y_units
      filename(13:13)=m
      filename(14:14)=d

    ELSE IF (time_convention == 'Absolute_short   ') THEN
      filename(10:10)=y_units
    1       FORMAT (i3)
    2       FORMAT ('0',i2)
    3       FORMAT ('00',i1)
      IF (dayno <  100) THEN
        IF (dayno <  10) THEN
          WRITE(cdayno,3) dayno
        ELSE
          WRITE(cdayno,2) dayno
        END IF
      ELSE
        WRITE(cdayno,1) dayno
      END IF
      filename(11:13)= cdayno
      filename(14:14)=h
    END IF
  ELSE IF (style =='C') THEN     ! style C - Character date
    IF (time_convention == 'Absolute_standard') THEN

! decades meansured relative to 1800
      y_tens  = char_id(mod(yyyy/10,36)+1)
      filename(10:10)=y_tens
      filename(11:11)=y_units
            
      IF (meanlev == 0) THEN
        filename(12:14) = month_3char(mm)
              
      ELSE ! means date routine called is at beginning of next period
        mon=mm-1
        IF (mon == 0) THEN
          mon = 12
        END IF
        IF (filetype_2 == 'm') THEN
          filename(12:14) = month_3char(mon)        
        ELSE
          filename(12:14) = season_3char(mon)
        END IF
      END IF
    ELSE IF (time_convention == 'Absolute_long    ') THEN

! centuries measured from 0 i.e. 1992 as j92, with wraparound at 3600
      y_hunds = char_id(mod(yyyy/100,36)+1)
      y_tens  = char_id((yyyy-(yyyy/100)*100)/10+1)
      filename(10:10)=y_hunds
      filename(11:11)=y_tens
      filename(12:12)=y_units
      
      IF (meanlev == 0) THEN
        filename(13:14) = month_2char(mm)      
        
      ELSE ! means date routine called is at beginning of next period
        mon=mm-1
        IF (mon == 0) THEN
          mon = 12
        END IF
        IF (filetype_2 == 'm') THEN
          filename(13:14) = season_2char(mon)        
        ELSE
          filename(13:14) = season_2char(mon)        
        END IF
      END IF
    END IF

  ELSE IF (style =='D') THEN
! style D - for time convention Absolute Date stamp 
! Replace with data stamp
      Write( filename(10:13), '(i4.4)' )  yyyy  
      Write( filename(14:15), '(i2.2)' )  mm  
      Write( filename(16:17), '(i2.2)' )  dd  
      IF (whole_day_flag == 0) THEN
        filename(18:18)='_'
        Write( filename(19:20), '(i2.2)' )  hh    
      END IF 
  
  ELSE IF  (style =='E') THEN
! Style E: 4 digit year and 3 char month or season
    Write( filename(10:13), '(i4.4)' )  yyyy       
    IF (meanlev == 0) THEN
      filename(14:16)=month_3char(mm)        
      
    ELSE ! means date routine called is at beginning of next period
      mon=mm-1
      IF (mon == 0) THEN
        mon = 12
      END IF
      IF (filetype_2 == 'm') THEN
        filename(14:16) = month_3char(mon)
      ELSE
        filename(14:16) = season_3char(mon)
      END IF
    END IF  
  END IF
END IF

! Append the path to DATAM to the generated filename
CALL fort_get_env('DATAM', 5, datam, filenamelength, icode)
IF (icode /= 0) THEN
  icode =  10  ! Force an abort
  cmessage = 'GET_NAME: Failed to get value of $DATAM from environment'
  GO TO 9999
END IF

! Check we won't overrun when we merge the two
IF (LEN_TRIM(datam) + LEN_TRIM(filename) > filenamelength) THEN
  icode = 20
  cmessage = 'GET_NAME: Full filename (including path) too long to store.'
  GO TO 9999
END IF

filename = TRIM(datam) // '/' // TRIM(filename)

 9999 CONTINUE

IF (lhook) CALL dr_hook('GET_NAME',zhook_out,zhook_handle)
RETURN
! ----------------------------------------------------------------------
END SUBROUTINE get_name
