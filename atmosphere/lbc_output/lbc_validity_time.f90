! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate Validity Time for LBCs.
!
! Subroutine Interface:

      SUBROUTINE LBC_Validity_Time ( LCal360 )

      USE Submodel_Mod
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE Control_Max_Sizes
      USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim
      IMPLICIT NONE

!
! Description:
!   Calculates the validity time for the LBCs.
!
! Method:
!   Takes the model basis time and derives the validity time
!   through a call to SEC2TIME.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:
!
! Global variables

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
! Include file : parlbcs.h
!
! Must be called after parvars.h
!
! Description:
!   Contains variables in connection with generating LBCs.
!
! -----------------------------------------------------------
! Stash Codes for LBCs in Section 32 (and section 31 except tracers) 
!
      Integer, Parameter :: lbc_stashcode_orog    = 1 
      Integer, Parameter :: lbc_stashcode_u       = 2 
      Integer, Parameter :: lbc_stashcode_v       = 3 
      Integer, Parameter :: lbc_stashcode_w       = 4 
      Integer, Parameter :: lbc_stashcode_density = 5 
      Integer, Parameter :: lbc_stashcode_theta   = 6 
      Integer, Parameter :: lbc_stashcode_q       = 7 
      Integer, Parameter :: lbc_stashcode_qcl     = 8 
      Integer, Parameter :: lbc_stashcode_qcf     = 9 
      Integer, Parameter :: lbc_stashcode_exner   = 10 
      Integer, Parameter :: lbc_stashcode_u_adv   = 11 
      Integer, Parameter :: lbc_stashcode_v_adv   = 12 
      Integer, Parameter :: lbc_stashcode_w_adv   = 13 
      Integer, Parameter :: lbc_stashcode_qcf2    = 14 
      Integer, Parameter :: lbc_stashcode_qrain   = 15 
      Integer, Parameter :: lbc_stashcode_qgraup  = 16 
      Integer, Parameter :: lbc_stashcode_cf_bulk = 17 
      Integer, Parameter :: lbc_stashcode_cf_liquid = 18 
      Integer, Parameter :: lbc_stashcode_cf_frozen = 19 
      Integer, Parameter :: lbc_stashcode_murk      = 20 
      Integer, Parameter :: lbc_stashcode_free_tracer = 21 
      Integer, Parameter :: lbc_stashcode_ukca_tracer = 22 
      Integer, Parameter :: lbc_stashcode_dust_div1 = 23
      Integer, Parameter :: lbc_stashcode_dust_div2 = 24
      Integer, Parameter :: lbc_stashcode_dust_div3 = 25
      Integer, Parameter :: lbc_stashcode_dust_div4 = 26
      Integer, Parameter :: lbc_stashcode_dust_div5 = 27
      Integer, Parameter :: lbc_stashcode_dust_div6 = 28
      Integer, Parameter :: lbc_stashcode_so2      = 29
      Integer, Parameter :: lbc_stashcode_dms      = 30
      Integer, Parameter :: lbc_stashcode_so4_aitken = 31
      Integer, Parameter :: lbc_stashcode_so4_accu = 32
      Integer, Parameter :: lbc_stashcode_so4_diss = 33
      Integer, Parameter :: lbc_stashcode_nh3      = 35
      Integer, Parameter :: lbc_stashcode_soot_new = 36
      Integer, Parameter :: lbc_stashcode_soot_agd = 37
      Integer, Parameter :: lbc_stashcode_soot_cld = 38
      Integer, Parameter :: lbc_stashcode_bmass_new = 39
      Integer, Parameter :: lbc_stashcode_bmass_agd = 40
      Integer, Parameter :: lbc_stashcode_bmass_cld = 41
      Integer, Parameter :: lbc_stashcode_ocff_new = 42
      Integer, Parameter :: lbc_stashcode_ocff_agd = 43
      Integer, Parameter :: lbc_stashcode_ocff_cld = 44
      Integer, Parameter :: lbc_stashcode_nitr_acc = 45
      Integer, Parameter :: lbc_stashcode_nitr_diss = 46

! -----------------------------------------------------------
!     Data Time for LBC data
      Integer :: LBC_DT_Year
      Integer :: LBC_DT_Month
      Integer :: LBC_DT_Day
      Integer :: LBC_DT_Hour
      Integer :: LBC_DT_Min
      Integer :: LBC_DT_Sec
      Integer :: LBC_DT_DayNo

      COMMON /LBC_DT/ LBC_DT_Year, LBC_DT_Month, LBC_DT_Day,            &
         LBC_DT_Hour, LBC_DT_Min,  LBC_DT_Sec,   LBC_DT_DayNo

! -----------------------------------------------------------

!     Validity Time for LBC data
      Integer :: LBC_VT_Year
      Integer :: LBC_VT_Month
      Integer :: LBC_VT_Day
      Integer :: LBC_VT_Hour
      Integer :: LBC_VT_Min
      Integer :: LBC_VT_Sec
      Integer :: LBC_VT_DayNo

      COMMON /LBC_VT/ LBC_VT_Year, LBC_VT_Month, LBC_VT_Day,            &
         LBC_VT_Hour, LBC_VT_Min,  LBC_VT_Sec,   LBC_VT_DayNo

! -----------------------------------------------------------

      Integer, Parameter :: P_Src_Grid = 2
      Integer, Parameter :: P_LBC_Grid = 4

!     1 : Start Latitude
!     2 : Start Longitude
!     3 : Row Length
!     4 : Rows

      Real :: Src_Grid (Nfld_max, P_Src_Grid)
      Real :: LBC_Grid (Nfld_max, P_LBC_Grid)

      COMMON /LBC_Grids/ Src_Grid, LBC_Grid

! -------------------------------------------------------------

      Integer :: LBC_Global_LenRimA (Nfld_max, Nhalo_max)
      Integer :: LBC_Interp_LenRimA (Nfld_max, Nhalo_max)

      COMMON /LBC_Sizes/ LBC_Global_LenRimA, LBC_Interp_LenRimA

! -------------------------------------------------------------

! Subroutine arguments

      Logical :: LCal360     !  360/365 day calander indicator

! Local scalars:

      Integer  :: yy,mm,dd,hr,mn,ss,day_no,sec

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Function & Subroutine calls:

      External SEC2TIME

!- End of header

! ---------------------------------------------
! Determine Validity Time from Model Basis Time
! ---------------------------------------------

      IF (lhook) CALL dr_hook('LBC_VALIDITY_TIME',zhook_in,zhook_handle)
      sec = STEPim(a_im) * secs_per_periodim(a_im) /                    &
                steps_per_periodim(a_im)

! DEPENDS ON: sec2time
      Call SEC2TIME(0,SEC,BASIS_TIME_DAYS,BASIS_TIME_SECS,              &
                        YY,MM,DD,HR,MN,SS,DAY_NO,LCAL360)

! ------------------------------
! Store validity time in parlbcs
! ------------------------------
      LBC_VT_Year  = YY
      LBC_VT_Month = MM
      LBC_VT_Day   = DD
      LBC_VT_Hour  = HR
      LBC_VT_Min   = MN
      LBC_VT_DayNo = Day_No

      IF (lhook) CALL dr_hook('LBC_VALIDITY_TIME',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE LBC_Validity_Time
