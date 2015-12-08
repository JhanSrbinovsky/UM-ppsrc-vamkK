! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DIAGO --------------------------------------------------
!LL
!LL  Purpose : Provide Statistics on 'Observation-Model' Increments
!LL
!LL            Calculate and print out the following :
!LL            Mean Increment
!LL            Mean Square Increment
!LL            Extreme Increment and position (Lat/Long) of obs
!LL
!LL            Can be used on Individual or group of AC Obs Types
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE diago_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE DIAGO (CSUBR,KTYPE,KCALL,                              &
     &           OBS_INCR,NORMF,OBS_LAT,OBS_LONG,                       &
     &           LMISSD,LENOB,LENOBT,NPTOBT,                            &
     &           NO_ANAL_LEVS,NO_ANAL_VAR)

      USE um_input_control_mod, ONLY : model_domain
      USE domain_params, ONLY: mt_global

      USE conversions_mod, ONLY: recip_pi_over_180
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx
      USE eqtoll_mod, ONLY: eqtoll
      USE w_coeff_mod, ONLY: w_coeff
      USE w_eqtoll_mod, ONLY: w_eqtoll
      IMPLICIT NONE

      CHARACTER(LEN=*) CSUBR   !  Subroutine where DIAGO is called from
      INTEGER KTYPE,KCALL,LENOB,LENOBT,NPTOBT,NO_ANAL_LEVS
      INTEGER NO_ANAL_VAR
      REAL                                                              &
     & OBS_INCR(LENOB+1,NO_ANAL_LEVS,NO_ANAL_VAR),                      &
     & NORMF (LENOB+1,NO_ANAL_LEVS), OBS_LAT (LENOBT+1),                &
     & OBS_LONG (LENOBT+1)
      LOGICAL LMISSD(LENOB+1,NO_ANAL_LEVS)

!     INTENT=IN---------------------------------------------------------
!     KTYPE    : AC Observation Type
!     KCALL    : Indicator to where DIAGO called from
!              : 1 = FROM AC
!                2 = FROM VAN## ; IPASS=1 ; BEFORE VERTICAL FILTERING
!                3 = FROM VAN## ; IPASS=1 ; AFTER  VERTICAL FILTERING
!                4 = FROM VAN## ; IPASS=2 ; BEFORE VERTICAL FILTERING
!                5 = FROM VAN## ; IPASS=2 ; AFTER  VERTICAL FILTERING
!                6 = from DIAGOCC ; IPASS=1 ; for multi-level cloud
!                7 = from DIAGOCC ; IPASS=1 ; for low cloud
!                8 = from DIAGOCC ; IPASS=1 ; for medium cloud
!                9 = from DIAGOCC ; IPASS=1 ; for high cloud
!               10 = from DIAGOCC ; IPASS=1 ; for total cloud
!                VAN## is Vertical Analysis routine calling DIAGO.
!     INTENT=INOUT------------------------------------------------------
!     INC      : Increments - P*, T, U, RH
!     VINC     : Increments - V
!     NORMF    : Normalisation factor
!     OBS_LAT  : Observation co-latitudes
!     OBS_LONG : Observation longitudes
!     LMISSD   : Array to control which obs are used on each level
!     LENOB    : No of obs in group of observation types
!     LENOBT   : No of obs for Observation Type KTYPE
!     NPTOBT   : Offset to first observation for type KTYPE
!     NO_ANAL_LEVS : No of analysis levels
!     INTENT=OUT--------------------------------------------------------
!*   ------------------------------------------------------------------

!-----AC common blocks
! MPPAC start
      ! dimensions for obs allocated in subroutine AC
      ! dynamically allocated in AC and RDOBS

      ! dimension inxdim allocated in subroutine HORINF
      INTEGER,PARAMETER:: inxdim    = 15000

      ! common for Statistics Calcs in DIAGO ; Prints in RDOBS,GETOBS
      REAL :: R_STAT(MODEL_LEVELS_MAX,0:8)
      REAL :: S_STAT(MODEL_LEVELS_MAX,0:8)
      INTEGER :: COUNTA(NOBTYPMX)
      INTEGER :: COUNTB(NOBTYPMX)
      INTEGER :: COUNTC(NOBTYPMX)

      COMMON /mpp_ac/ R_STAT,S_STAT,COUNTA,COUNTB,COUNTC

      ! common to pass longitudes and latitudes for edges of local
      ! box from setcona to RDOBS and HINTCF

      REAL :: LONG_E
      REAL :: LONG_W
      REAL :: LAT_N,LAT_S
      REAL :: LONG_W_MODEL
      REAL :: LONG_E_MODEL
      COMMON/latlonmax/                                                 &
     &  LONG_E,LONG_W,LAT_N,LAT_S,LONG_W_MODEL,LONG_E_MODEL
! MPPAC end
      real a_maxinc(MODEL_LEVELS_MAX)
      integer iproc,istat
!-----------------------------------------------------------------------
!*LCOMDECK COMACP
! Description:
!   Declares variables for the common block COMACP. This controls the
!  execution of the assimilation code.
!-----------------------------------------------------------------------

! Declare parameters:
      INTEGER      MODEACP
        PARAMETER (MODEACP = 36)
      INTEGER      NANALTYP
        PARAMETER (NANALTYP = 30)
      INTEGER      NRADARS
        PARAMETER (NRADARS  = 15)

! Declare variables:
      INTEGER NACT,                        NPROG
      INTEGER AC_OBS_TYPES(NOBTYPMX),      LACT(NOBTYPMX)
      INTEGER GROUP_INDEX(NOBTYPMX),       TYPE_INDEX(NOBTYPMX)
      INTEGER GROUP_FIRST(NOBTYPMX),       GROUP_LAST(NOBTYPMX)
      INTEGER OBS_UNITNO,                  OBS_FORMAT
      INTEGER NO_OBS_FILES,                DIAG_RDOBS
      INTEGER IUNITNO,                     MGEOWT
      INTEGER N_GROUPS,                    GROUP_NO(NOBTYPMX)
      INTEGER MHCORFN,                     MACDIAG(MODEACP)
      INTEGER MWTFN,                       MDATADFN
      INTEGER NPASS_RF,                    NSLABS_SCFACT(MODEL_LEVELS_MAX)
      INTEGER NO_SCFACT(NOBTYPMX),         IOMITOBS(NANALTYP)
      INTEGER MASTER_AC_TYPES(NOBTYPMX),   DEF_AC_ORDER(NOBTYPMX)
      INTEGER DEF_NO_ITERATIONS(NOBTYPMX), DEF_INTERVAL_ITER(NOBTYPMX)
      INTEGER DEF_NO_ANAL_LEVS(NOBTYPMX),  DEF_NO_WT_LEVS(NOBTYPMX)
      INTEGER DEF_MODE_HANAL(NOBTYPMX),    LENACT(NOBTYPMX)
      INTEGER DEF_OBTHIN(NOBTYPMX),        MVINT205
      INTEGER MRAMPFN,                     MGLOSSFN
      INTEGER      LHN_RANGE
      INTEGER      NPASS_RF_LHN

      INTEGER WB_LonOffset,                WB_LonPts
      INTEGER WB_LatOffset,                WB_LatPts

      REAL    OBTIME_NOM
      REAL    VERT_FILT
      REAL    GEOWT_H(ROWS_MAX -1)
      REAL    TROPLAT
      REAL    GEOWT_V(MODEL_LEVELS_MAX)
      REAL    VERT_COR_SCALE(MODEL_LEVELS_MAX, 4)
      REAL    VERT_CUTOFF_SL
      REAL    VERT_CUTOFF_BW
      REAL    VERT_CUTOFF_BH
      REAL    NON_DIV_COR
      REAL NON_DIV_COR_10M
      REAL    SPEED_LIMIT305
      REAL    TROPINT
      REAL    TIMEF_START
      REAL    TIMEF_OBTIME
      REAL    TIMEF_END
      REAL    CSCFACT_H(ROWS_MAX)
      REAL    CSCFACT_V(MODEL_LEVELS_MAX)
      REAL    DEF_TIMEB(NOBTYPMX)
      REAL    DEF_TIMEA(NOBTYPMX)
      REAL    DEF_TGETOBB(NOBTYPMX)
      REAL    DEF_TGETOBA(NOBTYPMX)
      REAL    DEF_CSCALE_START(NOBTYPMX)
      REAL    DEF_CSCALE_OBTIME(NOBTYPMX)
      REAL    DEF_CSCALE_END(NOBTYPMX)
      REAL    DEF_RADINF(NOBTYPMX)
      REAL    WB_LAT_CC(ROWS_MAX)
      REAL    WB_VERT_V(MODEL_LEVELS_MAX)
      REAL    WB_LAND_FACTOR
      REAL         RADAR_LAT(NRADARS)
      REAL         RADAR_LON(NRADARS)
      REAL         RADAR_RANGE_MAX
      REAL         EPSILON_LHN
      REAL         RELAX_CF_LHN
      REAL         F1_506 , F2_506 , F3_506
      REAL         ALPHA_LHN
      REAL         LHN_LIMIT
      REAL         FI_SCALE_LHN

      REAL    DEF_NUDGE_NH(NOBTYPMX)
      REAL    DEF_NUDGE_TR(NOBTYPMX)
      REAL    DEF_NUDGE_SH(NOBTYPMX)
      REAL    DEF_NUDGE_LAM(NOBTYPMX)

      REAL    DEF_FI_VAR_FACTOR(NOBTYPMX)
      REAL    FI_SCALE
      REAL    FI_SCALE_FACTOR(MODEL_LEVELS_MAX)
      REAL    DF_SCALE
      REAL    DF_SCALE_LEV(MODEL_LEVELS_MAX)
      REAL    DF_COEFF(MODEL_LEVELS_MAX)
      REAL    THRESH_DL
      REAL    THRESH_LM
      REAL    THRESH_MH
      REAL    THRESH_RMSF
      REAL    RADAR_RANGE
      REAL    NORTHLAT, SOUTHLAT, WESTLON, EASTLON
      REAL    VERT_COR_AERO

      LOGICAL LGEO
      LOGICAL LHYDR
      LOGICAL LHYDROL
      LOGICAL LSYN
      LOGICAL LTIMER_AC
      LOGICAL LAC_UARS
      LOGICAL LAC_MES
      LOGICAL LWBAL_SF,     LWBAL_UA
      LOGICAL WB_THETA_UA, WB_LAND_SCALE, WB_THETA_SF
      LOGICAL LRADAR (NRADARS)
      LOGICAL L_LATLON_PRVER
      LOGICAL L_MOPS_EQUALS_RH
      LOGICAL LCHECK_GRID
      LOGICAL      L_506_OBERR
      LOGICAL      L_LHN , L_LHN_SCALE
      LOGICAL      L_LHN_SEARCH , LHN_DIAG
      LOGICAL      L_VERIF_RANGE
      LOGICAL      L_LHN_LIMIT
      LOGICAL      L_LHN_FACT
      LOGICAL      L_LHN_FILT
      LOGICAL L_OBS_CHECK
      LOGICAL  REMOVE_NEG_LH
      LOGICAL USE_CONV_IN_MOPS

      COMMON /COMACP/ NACT,N_GROUPS,NPROG,                              &
     &  AC_OBS_TYPES,     LACT,              GROUP_NO,                  &
     &  LENACT,           LWBAL_SF,          LWBAL_UA,                  &
     &  LTIMER_AC,        LGEO,              LHYDR,                     &
     &  MGEOWT,           LSYN,              LAC_UARS,                  &
     &  OBS_UNITNO,       OBS_FORMAT,        NO_OBS_FILES,              &
     &  L_OBS_CHECK,                                                    &
     &  DIAG_RDOBS,       IUNITNO,           MVINT205,                  &
     &  MHCORFN,          MACDIAG,                                      &
     &  DEF_AC_ORDER,     DEF_NO_ITERATIONS, DEF_INTERVAL_ITER,         &
     &  MWTFN,            MDATADFN,          NSLABS_SCFACT,             &
     &  NO_SCFACT,        NPASS_RF,          MRAMPFN,                   &
     &  IOMITOBS,         TROPINT,           SPEED_LIMIT305,            &
     &  GEOWT_H,          GEOWT_V,           MGLOSSFN,                  &
     &  NON_DIV_COR,      TROPLAT,           VERT_FILT,                 &
     &  NON_DIV_COR_10M,                                                &
     &  VERT_COR_SCALE,                                                 &
     &  VERT_CUTOFF_SL,   VERT_CUTOFF_BW,    VERT_CUTOFF_BH,            &
     &  TIMEF_START,      TIMEF_OBTIME,      TIMEF_END,                 &
     &  CSCFACT_H,        CSCFACT_V,                                    &
     &  MASTER_AC_TYPES,                                                &
     &  DEF_NO_ANAL_LEVS, DEF_NO_WT_LEVS,    DEF_MODE_HANAL,            &
     &  DEF_TIMEB,        DEF_TIMEA,         DEF_TGETOBB,               &
     &  DEF_TGETOBA,      OBTIME_NOM,        DEF_OBTHIN,                &
     &  DEF_RADINF,       DEF_CSCALE_START,  DEF_CSCALE_OBTIME,         &
     &  DEF_CSCALE_END,                                                 &
     &  DEF_NUDGE_NH,     DEF_NUDGE_TR,      DEF_NUDGE_SH,              &
     &  DEF_NUDGE_LAM,                                                  &
     &  WB_LonOffset,     WB_LonPts,         WB_LatOffset,              &
     &  WB_LatPts,                                                      &
     &  GROUP_INDEX,      GROUP_FIRST,       GROUP_LAST,                &
     &  TYPE_INDEX,                                                     &
     &  FI_SCALE,         FI_SCALE_FACTOR,   DEF_FI_VAR_FACTOR,         &
     &  DF_SCALE,         DF_SCALE_LEV,      DF_COEFF,                  &
     &  LAC_MES,                                                        &
     &  THRESH_DL,        THRESH_LM,         THRESH_MH,                 &
     &  THRESH_RMSF,                                                    &
     &  RADAR_RANGE,      LRADAR,            LHYDROL,                   &
     &  L_LATLON_PRVER,   NORTHLAT,          SOUTHLAT,                  &
     &  WESTLON,          EASTLON,           L_MOPS_EQUALS_RH,          &
     &  LHN_RANGE ,  L_LHN , L_LHN_SCALE ,                              &
     &  L_LHN_SEARCH , LHN_DIAG , REMOVE_NEG_LH,                        &
     &  RADAR_LAT , RADAR_LON , RADAR_RANGE_MAX ,                       &
     &  EPSILON_LHN , ALPHA_LHN , RELAX_CF_LHN , LHN_LIMIT ,            &
     &  F1_506 , F2_506 , F3_506 ,                                      &
     &  L_506_OBERR , L_VERIF_RANGE , L_LHN_LIMIT , L_LHN_FACT ,        &
     &  L_LHN_FILT , FI_SCALE_LHN , NPASS_RF_LHN ,                      &
     &  VERT_COR_AERO,    LCHECK_GRID,                                  &
     &  WB_LAT_CC,        WB_VERT_V,         WB_LAND_FACTOR,            &
     & WB_THETA_UA,      WB_LAND_SCALE,   WB_THETA_SF, USE_CONV_IN_MOPS

! COMACP end
! COMACDG

      INTEGER,PARAMETER :: NLDACP=5
      INTEGER,PARAMETER :: NDACOP=200
      INTEGER,PARAMETER :: NDACP=100
      INTEGER,PARAMETER :: NDACVP=11

      LOGICAL :: LDIAGAC
      LOGICAL :: LLDAC(NLDACP)
      LOGICAL :: LLDAG0
      LOGICAL :: LLDAG(NLDACP)
      LOGICAL :: LLBAND
      LOGICAL :: LTITLE
      LOGICAL :: LRMS
      LOGICAL :: LNORMF
      LOGICAL :: LTEMP
      LOGICAL :: LVERIF

      INTEGER :: MDIAGTPS(NOBTYPMX)
      INTEGER :: MDIAG
      INTEGER :: NDGPRT
      INTEGER :: NDACO
      INTEGER :: MDACO(NDACOP)
      INTEGER :: NDAC
      INTEGER :: NDACPRT
      INTEGER :: MDAC(NDACP)
      INTEGER :: NDACV
      INTEGER :: MODACO

      REAL :: DLATN
      REAL :: DLATS
      REAL :: DLONGW
      REAL :: DLONGE
      REAL :: DAGLAT
      REAL :: DAGLON
      REAL :: DACV(NDACP,NDACVP)

      CHARACTER(LEN=10) :: CACV(NDACVP)
      CHARACTER(LEN=12) :: CACT1
      CHARACTER(LEN=32) :: CACT2
      CHARACTER(LEN=30) :: CACT3
      CHARACTER(LEN=55) :: CACT4

      COMMON /COMACDG/ LDIAGAC,LLDAC,LLDAG0,LLDAG,                      &
     &  LRMS,LNORMF,LTEMP,LVERIF,                                       &
     &  NDACO,MDACO,NDAC,NDACPRT,MDAC,NDACV,MODACO,                     &
     &  MDIAGTPS,MDIAG,LLBAND,LTITLE,NDGPRT,                            &
     &  DLATN,DLATS,DLONGW,DLONGE,DAGLAT,DAGLON,                        &
     &  DACV,CACV,CACT1,CACT2,CACT3,CACT4

! COMACDG end
!-----------------------------------------------------------------------
!L  COMDECK DOCACAG
!L  ---------------
!L    COMACDG CONTAINS PARAMETERS CONTROLLING DIAGNOSTIC PRINTOUT
!L
!L    LDIAGAC =.FALSE. SWITCHES OFF ALL DIAGNOSTICS
!L    LLDAC(#)=.FALSE. SWITCHES OFF DIAGNOSTICS OF TYPE #
!L
!L    DIAGNOSTIC TYPES & THEIR ASSOCIATED PARAMETERS ARE :-
!L    --------------------------------------------------
!L
!L    LLDAC(1) = .TRUE.
!L    -----------------
!L
!L    DETAILED DIAGNOSTICS ON A SHORT LIST OF OBSERVATIONS
!L    TO BE DONE FROM AC ETC.
!L    THIS TYPE OF DIAGNOSTIC IS SET-UP & USED BY TYPE 5 IF LLDAG0=T.
!L    AFTER EACH ITERATION THE LIST OF OBS INFLUENCING THE SELECTED PNT
!L    IS ADDED TO MDACO. THUS FOR SUBSEQUENT ITERATIONS DETAILED
!L    DIAGNOSTICS ON THESE OBS ARE OBTAINED. TO GET THIS MODE SET
!L    LLDAC(5)=T LLDAG0=T.
!L
!L    THIS TYPE CAN ALSO BE SWITCHED ON INDEPENDENTLY OF TYPE 5 BY :-
!L    MODACO  MODE FOR SETTING UP LIST MDAC:-
!L            1 TAKE FIRST NDACP OBS IN CURRENT VECTOR
!L            2 SEARCH FOR THOSE OBS IN MDACO DEFINED IN NAMELIST ADIAG
!L            3 SEARCH FOR THOSE OBS IN LTD AREA FROM    NAMELIST ADIAG
!L
!L    NDACOP  MAX NO OF OBS ON WHICH DIAGNOSTICS ARE REQUIRED
!L    NDACO       NO OF OBS ON WHICH DIAGNOSTICS ARE REQUIRED
!L    MDACO     LIST OF OBS ON WHICH DIAGNOSTICS ARE REQUIRED
!L              POINTS TO POSITION OF OB IN COMOBS.
!L
!L    NDACP  MAX NO OF OBS OF CURRENT TYPE ON WHICH DIAGS REQD.
!L    NDAC       NO OF OBS OF CURRENT TYPE ON WHICH DIAGS REQD.
!L    MDAC     LIST OF OBS OF CURRENT TYPE ON WHICH DIAGS REQD.
!L              POINTS TO POSITION OF OB IN CURRENT VECTORS FROM GETOBS.
!L
!L    NDACVP  MAX NO OF PARAMETERS WHICH CAN BE STORED FOR EACH OB
!L    NDACV       NO OF PARAMETERS WHICH ARE    STORED FOR EACH OB
!L    DACV        STORED PARAMETERS FOR EACH OB
!L    CACV        DESCRIPTION OF EACH PARAMETER
!L    CACT1/2/3/4 TITLES SET UP BY DACOIN
!L    LTITLE      CONTROLS TITLING OF OUTPUT FROM DACOUT
!L
!L    LLDAC(2) = .TRUE.
!L    -----------------
!L    STATISTICS OF OBSERVATION-MODEL INCREMENTS
!L    PRINTED OUT IN DIAGO
!L    MDIAGTPS: LIST OF TYPES TO BE PROCESSED BY DIAGO.
!L              SET MDIAGTPS(1)=0 FOR 'ALL'(DEFAULT).
!L    MDIAG:    1 = ONLY CALLS FROM AC PROCESSED
!L              2 = ONLY CALLS FROM Van### BEFORE VRTF PROCESSED
!L              4 = ONLY CALLS FROM Van### AFTER VRTF PROCESSED
!L                  (Van### is group of vertical analysis routines)
!L              0 = ALL CALLS PROCESSED (DEFAULT)
!L              BINARY COMBINATIONS SUPPORTED.
!L    LLBAND:   F = GLOBAL STATISTICS (DEFAULT).
!L              T = SEPARATE STATISTICS FOR BANDS 90N-22N,
!L                                                22N-22S,
!L                                                22S-90S.
!L    LRMS   :   T/F = Print RMS/Mean Square Values. Default = T.
!L    LNORMF :   T : Use Normalisation Factors (NF) as weights (Default)
!L           :   F : Set NF to 1.0 (ie. no weights used).
!L    LTEMP  :   T : Print Temperature statistics. Default. Theta
!L                   Increments are converted to temperature increments
!L                   and p* is assumed to be 1000mb for all obs.
!L           :   F : Print Theta statistics.
!L    LVERIF :   T : Sets parameters to get statistics for verification
!L                   purposes. LRMS=T,LTEMP=T,LNORMF=F,LGEO=F,LHYDR=F
!L           :   F : Default.
!L
!L    LLDAC(3) = .TRUE.
!L    -----------------
!L    STATISTICS OF INCREMENTS ON ANALYSIS GRID
!L
!L    LLDAC(4) = .TRUE.
!L    -----------------
!L    STATISTICS OF INCREMENTS ON MODEL GRID
!L    PRINTED OUT IN MMSPT
!L    MEAN AND MEAN-SQUARE OF INCREMENTS PRINTED FOR :-
!L        1. THE WHOLE MODEL(NOT WITH GEOSTROPHIC INCREMENTS)
!L    AND 2. EACH LEVEL OF THE MODEL
!L
!L    LLDAC(5) = .TRUE.
!L    -----------------
!L    DETAILS OF OBS INFLUENCING A SPECIFIED POINT.
!L    PARAMETERS IN NAMELIST ADIAG :-
!L    DAGLAT  -  LATITUDE  OF SPECIFIED POINT. DEFAULT =  57.0
!L    DAGLON  -  LONGITUDE OF SPECIFIED POINT. DEFAULT = 340.0
!L    LLDAG0      = SWITCH ON OPTION 1 FOR THOSE OBS FOUND, SO THAT
!L                   THEIR DETAILS ARE PRINTED NEXT TIME THROUGH AC.
!L    LLDAG(NAVT) = DO DIAGNOSTIC  FOR THIS VARIABLE TYPE.
!L
!-----------------------------------------------------------------------
!--------------------------------------------------------------------
!LCOMDECK COMMG
!L-------------
      REAL          DLAT,DLONG,XLATN,XLONGW
      REAL          ELFPLAT,ELFPLON
      COMMON/COMMG/ DLAT,DLONG,XLATN,XLONGW,ELFPLAT,ELFPLON
!--------------------------------------------------------------------

!*L   Workspace Usage:--------------------------------------------------
!     -----------------------------------------------------------
!     Local work arrays
!     -----------------------------------------------------------
!**   Dynamic allocation
      LOGICAL OBSUSE(LENOBT+1)
      REAL WLT(LENOBT+1), WLN(LENOBT+1), WLN2(LENOBT+1)
      REAL COEFF1(LENOBT+1), COEFF2(LENOBT+1)
      REAL UWK(LENOBT+1,NO_ANAL_LEVS), VWK(LENOBT+1,NO_ANAL_LEVS)
!     WLT/WLN  : Real lat/long of wind obs
!     WLN2     : ELF  longitude of wind obs
!     COEFF1/2 : Coefficients to rotate wind components
!     UWK/VWK  : Rotated u/v components on real lat/long grid
!*    ------------------------------------------------------------------

!     Local arrays and variables

      INTEGER JVAR,JBAND,JLEV,JOB,JOBT,JNANL
      INTEGER NBAND,INOBS,INOBSLV,IOBT,ANAL_VAR
      REAL STAT(0:NO_ANAL_LEVS,8),LAT,LONG,UI,VI,NF,SQINC
      REAL BANDLAT(4),MAXINC,SUMUI,SUMVI,SUMNF,SUMSQI
      CHARACTER NS*1,WE*1,TITLE1*9,TITLE2*25,LAB1*10,LAB2*10
      CHARACTER(LEN=3) BANDC(4)
      LOGICAL LWIND

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     ANAL_VAR : Analysis varibale (=KTYPE/100)
!     BANDLAT  : Boundaries of latitude bands
!     BANDC    : Labels to print out boundaries of latitude bands
!     INOBS    : Number of observations in latitude band
!     INOBSLV  : Number of observations on level
!     IOBT     : Index to list of obs types in MDIAGTPS
!     JBAND    : Loop variable for latitude bands
!     JLEV     : Loop variable for analysis levels
!     JNANL    : Number of analysis levels to be processed
!     JOB      : Loop variable for observations
!     JOBT     : Loop variable for observation types
!     JVAR     : Loop variable for statistic variable
!     LAB1/2   : Labels for print out
!     LAT      : Absolute value of latitude of maximum increment
!     LONG     : Absolute value of longitude of maximum increment
!     LWIND    : Indicator for working with wind obs type
!     MAXINC   : Maximum increment for level
!     NBAND    : Number of latitude bands
!     NF       : Normalisation Factors
!     NS/WE    : Labels for printing out lat/long of maximum increment
!     SQINC    : Square of increments
!     STAT     : Array storing statistics data
!                STAT(0,#) Overall
!                STAT(JLEV,#) For level JLEV
!                STAT(#,1) No of observations on  level JLEV
!                STAT(#,2) Mean Norm factor
!                STAT(#,3) Mean Increment - p*, t, u, rh
!                STAT(#,4) Mean Increment - v
!                STAT(#,5) Mean square Increment
!                STAT(#,6) Extreme increment
!                STAT(#,7) Lat of extreme increment
!                STAT(#,8) Long of extreme increment
!     SUMNF    : Sum of norm. factors
!     SUMUI    : Sum of increments - p*, t, u, rh
!     SUMVI    : Sum of increments - v
!     SUMSQI   : Sum of squared increments
!     TITLE1/2 : Titles for print out
!     UI       : Increments - p*, t, u, rh
!     VI       : Increments - v

!     ------------------------------------------------------------------
      DATA BANDLAT / 0.0, 1.0472, 2.0944, 3.1416 /
      DATA BANDC/'90N','30N','30S','90S'/
!     -----------------------------------------------------------

      IF (lhook) CALL dr_hook('DIAGO',zhook_in,zhook_handle)
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('DIAGO   ',3)

!     Check that 'Obs-Model' statistics required
      IF (.NOT.LLDAC(2)) GO TO 9999

!     For verification purposes, only get statistics on first pass.
      IF (LVERIF.AND.(KCALL /= 2.AND.KCALL /= 3)) GO TO 9999

!     Check that statistics required for calling routine
      IF (MDIAG >  0 .AND. KCALL /= MDIAG) GO TO 9999

!     Check that statistics required for observation type KTYPE
      IF (MDIAGTPS(1) >  0) THEN
        IOBT=0
        DO JOBT=1,NOBTYPMX
          IF (KTYPE == MDIAGTPS(JOBT)) IOBT=JOBT
        END DO
        IF (IOBT == 0) GO TO 9999
      ENDIF

!---- Set up titles

      LAB1='    '
      LAB2='    '
      ANAL_VAR = KTYPE/100
      IF (ANAL_VAR == 1) THEN
        LAB1 = 'P* OBS INC'
      ELSEIF (ANAL_VAR == 2) THEN
        IF (LTEMP) THEN
          LAB1 = ' T OBS INC'
        ELSE
          LAB1 = 'TH OBS INC'
        ENDIF
      ELSEIF (ANAL_VAR == 3) THEN
        LAB1 = ' U OBS INC'
        LAB2 = ' V OBS INC'
      ELSEIF (ANAL_VAR == 4) THEN
        IF ( KCALL >  5 ) THEN
          LAB1 = 'CC OBS INC'
        ELSE
          LAB1 = 'RH OBS INC'
        ENDIF
      ELSEIF (ANAL_VAR == 9) THEN
        LAB1 = 'LOG(VIS) OBS INC'
      ELSE
      GOTO 9999
      ENDIF

      IF (KCALL == 2 .OR. KCALL == 3) THEN
        TITLE1 = 'IPASS=1'
      ELSEIF (KCALL == 4 .OR. KCALL == 5) THEN
        TITLE1 = 'IPASS=2'
      ELSE
        TITLE1 = '  '
      ENDIF

      TITLE2 = ' '
      IF (KTYPE == 201 .OR. KTYPE == 205 .OR. KTYPE == 206 .OR.         &
     &    KTYPE == 207 .OR. KTYPE == 208 .OR. KTYPE == 209 .OR.         &
     &    KTYPE == 301 .OR. KTYPE == 401 ) THEN
        IF (KCALL == 2 .OR. KCALL == 4) THEN
          TITLE2 = 'BEFORE VERTICAL FILTERING'
        ELSEIF (KCALL == 3 .OR. KCALL == 5) THEN
          TITLE2 = 'AFTER VERTICAL FILTERING'
        ELSE
          TITLE2 = ' '
        ENDIF
      ENDIF

      LWIND = ANAL_VAR == 3

      JNANL = NO_ANAL_LEVS

      IF (KTYPE == 202 .OR. KTYPE == 204 .OR.                           &
     &    KTYPE == 302 .OR. KTYPE == 304 .OR.                           &
     &    KTYPE == 305 .OR. KTYPE == 306 .OR.                           &
     &    KTYPE == 402 .OR. KTYPE == 404                                &
     &    .OR. KTYPE  ==  901) THEN
        JNANL=1
      ENDIF

!     ONLY ONE BAND OF STATISTICS ALLOWED FOR ELF GRID
      NBAND=1

!     Loop over latitude bands
      DO JBAND=1,NBAND

!     Initialise statistics data
      DO JVAR=1,8
        DO JLEV=0,NO_ANAL_LEVS
          STAT(JLEV,JVAR)=0.0
        ENDDO
      ENDDO

!L    1.0 Determine which observations to be used
      IF (LENOBT <= 0) THEN
       do jlev=1,jnanl
        do jvar=0,8
         s_stat(jlev,jvar)=0.0
        enddo
       enddo
      ELSE

!     1.1 Initialise so that all observations to be used
      DO JOB=1,LENOBT
        OBSUSE(JOB)=.TRUE.
      ENDDO

      IF (model_domain == mt_global) THEN
!     1.2 Skip observations outside latitude band
      IF (LLBAND) THEN
        DO JOB=1,LENOBT
          OBSUSE(JOB) = OBS_LAT(JOB) >= BANDLAT(JBAND)   .AND.          &
     &                  OBS_LAT(JOB) <  BANDLAT(JBAND+1)
        ENDDO
      ENDIF
      END IF

!     1.3 Count number of observations to be used
      INOBS = 0
      DO JOB=1,LENOBT
        IF (OBSUSE(JOB)) INOBS = INOBS+1
      ENDDO

      IF (model_domain /= mt_global) THEN

!       2. Co-ordinate tranformation  (ELF version only)
!       ------------------------------------------------

!       Get real Lat/Long of ELF observations

!       Lat/Lon values in OBS_LAT/OBS_LONG are for ELF grid.
!       First convert from radians to degrees.
        DO JOB=1,LENOBT
          WLT (JOB) = 90.0-OBS_LAT(JOB)*RECIP_PI_OVER_180
          WLN (JOB) = OBS_LONG(JOB)*RECIP_PI_OVER_180
          WLN2(JOB) = WLN(JOB)
        ENDDO
!       Get real Lat/Lon of wind obs in WLT/WLN, keep ELF lon in WLN2
        CALL EQTOLL (WLT,WLN,WLT,WLN,ELFPLAT,ELFPLON,LENOBT)

        IF (LWIND) THEN

!         Rotate wind components on ELF grid to real lat/lon system.

!         Calculate coefficients of rotation.
          CALL W_COEFF(COEFF1,COEFF2,WLN,WLN2,ELFPLAT,ELFPLON,LENOBT)
!         Rotation of wind components in loop over levels.
          DO JLEV=1,JNANL
            CALL W_EQTOLL (COEFF1,COEFF2,                               &
     &                     OBS_INCR(NPTOBT+1,JLEV,1),                   &
     &                     OBS_INCR(NPTOBT+1,JLEV,2),                   &
     &                     UWK(1,JLEV),VWK(1,JLEV),LENOBT,.TRUE.)
          ENDDO

        ENDIF
      END IF


!       3.0 Calculate statistics required
!       ---------------------------------

!-----  Loop over levels
        DO JLEV=1,JNANL

        INOBSLV=0
        SUMNF  = 0.0
        SUMUI  = 0.0
        SUMVI  = 0.0
        SUMSQI = 0.0
        MAXINC = 0.0

        IF (model_domain == mt_global) THEN
!-----  Loop over observations
        DO JOB=1,LENOBT

        NF    = 0.0
        UI    = 0.0
        VI    = 0.0
        SQINC = 0.0

        IF (OBSUSE(JOB) .AND. .NOT.LMISSD(NPTOBT+JOB,JLEV)) THEN

          INOBSLV = INOBSLV+1

          IF (LNORMF) THEN
!           Use NF as weighting
            NF = NORMF (NPTOBT+JOB,JLEV)
          ELSE
!           Do not use NF as weighting (set to 1.0)
            NF = 1.0
          ENDIF
          UI = OBS_INCR(NPTOBT+JOB,JLEV,1)
          IF (LWIND) VI = OBS_INCR(NPTOBT+JOB,JLEV,2)


!-----    3.1 Sum of Normalistion factors used as weights

!         For IPASS=1, the preliminary norm factor (Q1M) calculated in
!                      VAN### is used to weight the obs increments.
!         For IPASS=2, the final norm factor (Q2M) calculated in
!                      VAN### is used to weight the obs increments.

          SUMNF = SUMNF + NF

!-----    3.2 Weighted sum of increments - P*, T, U, RH
          SUMUI = SUMUI + NF*UI

!-----    3.3 Weighted sum of increments - V
          IF (LWIND) SUMVI = SUMVI + NF*VI

!-----    3.4 Weighted sum of squared increments
          SQINC = UI*UI
          IF (LWIND) SQINC = SQINC + VI*VI
          SUMSQI = SUMSQI + NF*SQINC

!-----    3.5 Determine maximum (squared) increment
          IF (SQINC  >   MAXINC) THEN
            MAXINC = SQINC
            IF (LWIND) THEN
!             For wind, extreme is maximum speed of vector increment
              STAT(JLEV,6) = SQRT(SQINC)
            ELSE
!             Extreme is increment
              STAT(JLEV,6) = UI
            ENDIF
!           Lat/Long of obs
            STAT(JLEV,7) = OBS_LAT(JOB)
            STAT(JLEV,8) = OBS_LONG(JOB)
!-----      Convert to degrees
            STAT(JLEV,7) = 90.0 - STAT(JLEV,7) * RECIP_PI_OVER_180
            STAT(JLEV,8) = STAT(JLEV,8) * RECIP_PI_OVER_180
            IF (STAT(JLEV,8) >  180.0) STAT(JLEV,8)=STAT(JLEV,8)-360.0
          END IF

        END IF

      END DO !JOB
      ELSE
!-----  Loop over observations
        DO JOB=1,LENOBT

        NF    = 0.0
        UI    = 0.0
        VI    = 0.0
        SQINC = 0.0

        IF (OBSUSE(JOB) .AND. .NOT.LMISSD(NPTOBT+JOB,JLEV)) THEN

          INOBSLV = INOBSLV+1

          IF (LNORMF) THEN
!           Use NF as weighting
            NF = NORMF (NPTOBT+JOB,JLEV)
          ELSE
!           Do not use NF as weighting (set to 1.0)
            NF = 1.0
          END IF
          IF (LWIND) THEN
            UI = UWK(JOB,JLEV)
            VI = VWK(JOB,JLEV)
          ELSE
            UI = OBS_INCR(NPTOBT+JOB,JLEV,1)
          END IF

!-----    3.1 Sum of Normalistion factors used as weights

!         For IPASS=1, the preliminary norm factor (Q1M) calculated in
!                      VAN### is used to weight the obs increments.
!         For IPASS=2, the final norm factor (Q2M) calculated in
!                      VAN### is used to weight the obs increments.

          SUMNF = SUMNF + NF

!-----    3.2 Weighted sum of increments - P*, T, U, RH
          SUMUI = SUMUI + NF*UI

!-----    3.3 Weighted sum of increments - V
          IF (LWIND) SUMVI = SUMVI + NF*VI

!-----    3.4 Weighted sum of squared increments
          SQINC = UI*UI
          IF (LWIND) SQINC = SQINC + VI*VI
          SUMSQI = SUMSQI + NF*SQINC

!-----    3.5 Determine maximum (squared) increment
          IF (SQINC  >   MAXINC) THEN
            MAXINC = SQINC
            IF (LWIND) THEN
!             For wind, extreme is maximum speed of vector increment
              STAT(JLEV,6) = SQRT(SQINC)
            ELSE
!             Extreme is increment
              STAT(JLEV,6) = UI
            END IF
!           Lat/Long of obs
            STAT(JLEV,7) = WLT(JOB)
            STAT(JLEV,8) = WLN(JOB)
!-----      Already in degrees
          END IF

        END IF

      END DO !JOB
      END IF  !  if GLOBAL
        s_stat(jlev,0)=maxinc
        s_stat(jlev,1)=real(inobslv)
        s_stat(jlev,2)=sumnf
        s_stat(jlev,3)=sumui
        if(lwind)s_stat(jlev,4)=sumvi
        s_stat(jlev,5)=sumsqi
        s_stat(jlev,6)=STAT(JLEV,6)
        s_stat(jlev,7)=STAT(JLEV,7)
        s_stat(jlev,8)=STAT(JLEV,8)
      END DO ! JLEV
      END IF   ! LENOBT <= 0
      if(mype == 0)then
!
!   set stat to zero
!
        do jlev=1,jnanl
          stat(jlev,1)=0.0
          stat(jlev,2)=0.0
          stat(jlev,3)=0.0
          stat(jlev,4)=0.0
          stat(jlev,5)=0.0
          a_maxinc(jlev)=0.0
        enddo
      endif  ! mype=0

! Gather stats on PE0
      if(mype == 0) then ! stats off PE0
        do jlev=1,jnanl
          stat(jlev,1)=s_stat(jlev,1)
          stat(jlev,2)=s_stat(jlev,2)
          stat(jlev,3)=s_stat(jlev,3)
          if(lwind)stat(jlev,4)=s_stat(jlev,4)
          stat(jlev,5)=s_stat(jlev,5)
          if(s_stat(jlev,0) >  a_maxinc(jlev))then
            a_maxinc(jlev)=s_stat(jlev,0)
            stat(jlev,6)=s_stat(jlev,6)
            stat(jlev,7)=s_stat(jlev,7)
            stat(jlev,8)=s_stat(jlev,8)
          endif
        enddo
      endif


      DO IPROC=1,NPROC-1
        IF(mype == 0) THEN  ! PE0 receives
          CALL GC_RRECV(IPROC,9*model_levels_max,IPROC,           &
                       ISTAT,r_stat,s_stat)
          do jlev=1,jnanl
            stat(jlev,1)=stat(jlev,1)+r_stat(jlev,1)
            stat(jlev,2)=stat(jlev,2)+r_stat(jlev,2)
            stat(jlev,3)=stat(jlev,3)+r_stat(jlev,3)
            if(lwind)stat(jlev,4)=stat(jlev,4)+r_stat(jlev,4)
            stat(jlev,5)=stat(jlev,5)+r_stat(jlev,5)
            if(r_stat(jlev,0) >  a_maxinc(jlev))then
              a_maxinc(jlev)=r_stat(jlev,0)
              stat(jlev,6)=r_stat(jlev,6)
              stat(jlev,7)=r_stat(jlev,7)
              stat(jlev,8)=r_stat(jlev,8)
            endif
          enddo
        ELSEIF(mype == IPROC) THEN  ! other PEs send
          CALL GC_RSEND(IPROC,9*model_levels_max,0,               &
                       ISTAT,r_stat,s_stat)
        ENDIF
      ENDDO



      if(mype == 0) then

!-----  3.6 Get means for level JLEV and accumulate for overall means.
      do jlev=1,jnanl
      inobslv=int(stat(jlev,1))
      sumnf=stat(jlev,2)
      sumui=stat(jlev,3)
      if(lwind)sumvi=stat(jlev,4)
      sumsqi=stat(jlev,5)

        IF (INOBSLV >  0) THEN

          STAT(JLEV,1) = REAL( INOBSLV )
          STAT(0,1)    = STAT(0,1)+STAT(JLEV,1)
          IF (SUMNF >  0.0) THEN
!           Mean weight
            STAT(JLEV,2) = SUMNF/STAT(JLEV,1)
            STAT(0,2)    = STAT(0,2)+SUMNF
!           Mean observation increment - P*, T, U, RH
            STAT(JLEV,3) = SUMUI/SUMNF
            STAT(0,3)    = STAT(0,3)+SUMUI
            IF (LWIND) THEN
!             Mean observation increment - V
              STAT(JLEV,4) = SUMVI/SUMNF
              STAT(0,4)    = STAT(0,4)+SUMVI
            ENDIF
!           Mean square observation increment
            STAT(JLEV,5) = SUMSQI/SUMNF
            STAT(0,5)    = STAT(0,5)+SUMSQI
!           Maximum Increment
            MAXINC = ABS(STAT(JLEV,6))
            IF (MAXINC  >=  ABS(STAT(0,6))) THEN
              STAT(0,6) = STAT(JLEV,6)
              STAT(0,7) = STAT(JLEV,7)
              STAT(0,8) = STAT(JLEV,8)
            ENDIF
          ENDIF

        ENDIF

      enddo

!       4. Get means for all levels
!       ---------------------------
        IF (STAT(0,1) >  0.0 .AND. STAT(0,2) >  0.0) THEN

          SUMNF = STAT(0,2)
!         Mean weight
          STAT(0,2) = STAT(0,2)/STAT(0,1)
!         Mean observation increment - P*, T, U, RH
          STAT(0,3) = STAT(0,3)/SUMNF
!         Mean observation increment - V
          IF (LWIND) STAT(0,4) = STAT(0,4)/SUMNF
!         Mean square observation increment
          STAT(0,5) = STAT(0,5)/SUMNF

        ENDIF

!       If p*, convert from pascals to mb
        IF (ANAL_VAR == 1) THEN
          STAT(0,3) = STAT(0,3) * 0.01    ! Mean Increment
          STAT(0,5) = STAT(0,5) * 0.0001  ! Mean Square Increment
          STAT(0,6) = STAT(0,6) * 0.01    ! Maximum Increment
        ENDIF

!       Get RMS values if required
        IF (LRMS) THEN
          DO JLEV=0,JNANL
            IF (STAT(JLEV,5) >  0.0) THEN
              STAT(JLEV,5) = SQRT (STAT(JLEV,5))
            ENDIF
          ENDDO
        ENDIF

      endif  !mype=0

!     5. Print out results
!     --------------------

      if(mype == 0)then
!---- 5.1 Print out titles

      IF (JBAND == 1) THEN
        IF (KCALL >  1 .AND. KCALL  <  6) THEN
          PRINT '(/,'' DIAGO called from '',A,T30,''Obs Type '',I5)',   &
     &    CSUBR,KTYPE
          PRINT '('' '',A9,2X,A25)', TITLE1,TITLE2
        ELSEIF (KCALL >= 6) THEN
          PRINT '(/,'' Fit to '',A,T22,''cloud (oktas)'')', CSUBR
        ELSE
          PRINT '(/,'' DIAGO called from '',A)', CSUBR
        ENDIF
      ENDIF
      IF (model_domain == mt_global) THEN
      IF (LLBAND) THEN
        PRINT '(/,'' LATITUDES '',A3,'' to '',A3,5X,                    &
     &  ''No of obs '',I6/)', BANDC(JBAND),BANDC(JBAND+1),INOBS
      ELSE
        PRINT '(/,'' LATITUDES '',A3,'' to '',A3/)',BANDC(1),BANDC(4)
      ENDIF
      ELSE
      PRINT '(/,'' WHOLE LIMITED AREA '',/)'
      END IF  ! if GLOBAL

      IF (LRMS) THEN
        PRINT '(A,3X,A10,3X,A10,2X,A)',                                 &
     &  ' LEVEL     N  NORM FACTOR',LAB1,LAB2,                          &
     &  '    RMS       EXTREME INC & POSITION'
      ELSE
        PRINT '(A,3X,A10,3X,A10,2X,A)',                                 &
     &  ' LEVEL     N  NORM FACTOR',LAB1,LAB2,                          &
     &  'MEAN SQUARE   EXTREME INC & POSITION'
      ENDIF

!---- 5.2 Print overall means

      NS = ' '
      IF (STAT(0,7) >  0.0) THEN
        NS = 'N'
      ELSE
        NS = 'S'
      ENDIF
      WE = ' '
      IF (STAT(0,8) >  0.0) THEN
        WE = 'E'
      ELSE
        WE = 'W'
      ENDIF
      INOBSLV = INT(STAT(0,1))
      LAT     = ABS(STAT(0,7))
      LONG    = ABS(STAT(0,8))
      PRINT '(A,I6,5E13.5,F5.1,A1,F6.1,A1)','   ALL',INOBSLV,           &
     & (STAT(0,JVAR),JVAR=2,6),LAT,NS,LONG,WE

!---- 5.3 Print means for individual levels

       IF (JNANL >  1) THEN
        DO JLEV=JNANL,1,-1
          INOBSLV = INT( STAT(JLEV,1) )
          IF (INOBSLV >  0) THEN
            NS = ' '
            IF (STAT(JLEV,7) >  0.0) THEN
              NS = 'N'
            ELSE
              NS = 'S'
            ENDIF
            WE = ' '
            IF (STAT(JLEV,8) >  0.0) THEN
              WE = 'E'
            ELSE
              WE = 'W'
            ENDIF
            LAT  = ABS(STAT(JLEV,7))
            LONG = ABS(STAT(JLEV,8))
            PRINT '(2I6,5E13.5,F5.1,A1,F6.1,A1)',                       &
     &      JLEV,INOBSLV,(STAT(JLEV,JVAR),JVAR=2,6),LAT,NS,LONG,WE
          ENDIF
        ENDDO
      ENDIF

      endif
      END DO ! jband

!   if jumping to 9999 a problem has occurred but treat as non fatal
!   since this is a diagnostic routine (Fixme: needs an ereport?)
9999   CONTINUE

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('DIAGO   ',4)
      IF (lhook) CALL dr_hook('DIAGO',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE DIAGO
END MODULE diago_mod
