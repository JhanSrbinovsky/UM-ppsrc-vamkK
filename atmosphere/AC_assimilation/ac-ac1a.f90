! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: AC
!
!    Purpose : Main control routine for AC Scheme
!     This makes the analysis correction to model fields at the start
!   of each timestep of model integration.  The corrections are
!   calculated at analysis grid points,then interpolated to the model
!   grid before being added to the model fields. The corrections are
!   weighted sums of increments (ob-model) obtained at analysis grid
!   points. The weights depend on the distance between observation and
!   grid point, Observational error and the timeliness of the
!   observation.  Observations are analysed sequentially by observation
!   type in the order given by the array lact in comdeck comacp. A
!   vertical profile of observational increments and errors at
!   observation points are obtained before the weighted sums are
!   calculated. Horizontal filtering of increments on the model grid may
!   be done,and hydrostatic and geostrophic increments are calculated on
!   the model grid if appropriate.mean and rms diagnostics may be
!   calculated.
!
!    This section of code is now only used for the assimilation of
!    rainfall data from MOPS.
!
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!*
!*L  Arguments and declarations:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE ac_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE AC (                                                   &
     &  P_LEVELS, Q_LEVELS, ROW_LENGTH, P_ROWS, BL_LEVELS,              &
     &  TNDVMAX, NOBSMAX, P_FIELD,                                      &
     &  TIMESTEP_NO, TIMESTEP,                                          &
     &  EXNER, PSTAR, PRESSURE,                                         &
     &  THETA, RH, QCL, QCF,                                            &
     &  CONV_CLD, LS_RAIN, LS_SNOW, CONV_RAIN, CONV_SNOW,               &
     &  LAYER_CLOUD, D_THETA_DT_CONV,D_THETA_DT_LS, RHCRIT,             &
     & OBS_FLAG,OBS,                                                    &
     &  STINDEX,                                                        &
     &  STLIST, LEN_STLIST, SI, SF,                                     &
     &  STASHWORK, STASH_LEVELS,                                        &
     &  NUM_STASH_LEVELS, STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,        &
     &  lambda_p,phi_p,L_regular,                                       &
     &  ICODE, CMESSAGE)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE PrintStatus_mod
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx, obs_info
      USE domain_params
      USE ac2_mod, ONLY: ac2
      USE getobs_mod, ONLY: getobs
      USE hmrtorh_mod, ONLY: hmrtorh
      USE rdobs_mod, ONLY: rdobs

      IMPLICIT NONE


! Imported arguments (INTENT=IN):
      INTEGER        P_LEVELS                ! Total number of levels
      INTEGER        Q_LEVELS                ! Number of wet levels
      INTEGER        BL_LEVELS               ! Bdy layer levels
      INTEGER        ROW_LENGTH              ! Number of points on row
      INTEGER        P_ROWS                  ! Number of rows (pstar)
      INTEGER        TNDVMAX                 ! Total no of data values
                                             ! in all obs files
      INTEGER        NOBSMAX                 ! Total no of obs
                                             ! in all obs files
      INTEGER        P_FIELD                 ! Number of points in
                                             ! mass field
      INTEGER        TIMESTEP_NO             ! Current model timestep
      REAL           CONV_CLD(P_FIELD, Q_LEVELS)  ! conv cld amount
      REAL           LS_RAIN(P_FIELD)        ! large scale rain rate
      REAL           LS_SNOW(P_FIELD)        ! large scale snow rate
      REAL           CONV_RAIN(P_FIELD)      ! convective rain rate
      REAL           CONV_SNOW(P_FIELD)      ! convective snow rate
                                             ! above rates diagnostic
      REAL           TIMESTEP                ! Timestep in seconds
      REAL           RHCRIT(Q_LEVELS)        ! Critical rh array
                                             ! for cloud
      REAL           D_THETA_DT_CONV(P_FIELD,Q_LEVELS)
                                       ! convective latent heating rate
      REAL           D_THETA_DT_LS(P_FIELD,Q_LEVELS)
                                      ! large scale latent heating rate
! Stash variables:
      REAL           STASHWORK(*)

      INTEGER        LEN_STLIST              ! Length of STLIST
      INTEGER        NUM_STASH_LEVELS        ! No of Stash levels lists
      INTEGER        NUM_STASH_PSEUDO        ! No of Stash pseudo lists
      INTEGER        STINDEX(2, *)
      INTEGER        STLIST(LEN_STLIST, *)
      INTEGER        SI(*)                   ! Stash Index
      INTEGER        STASH_LEVELS(NUM_STASH_LEVELS +1, *)
                                             ! Levels lists
      INTEGER        STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO +1, *)
                                             ! Pseudo lists

      LOGICAL        SF(*)                   ! Stash Flags

! Import/export arguments (INTENT=INOUT):
      REAL           EXNER(P_FIELD, P_LEVELS) ! exner on theta levels
      REAL           PRESSURE(P_FIELD, P_LEVELS) ! p  on theta levels
      REAL           LAYER_CLOUD (P_FIELD, Q_LEVELS) ! as name says
      REAL           PSTAR(P_FIELD)          ! Prognostic variable pstar
      REAL           THETA(P_FIELD, P_LEVELS)! Prognostic variable
                                             ! theta
      REAL           RH(P_FIELD, Q_LEVELS)   ! Prognostic variable hmr
      REAL           QCL(P_FIELD, Q_LEVELS)   ! Prognostic variable qcl
      REAL           QCF(P_FIELD, Q_LEVELS)   ! Prognostic variable qcf

! Exported arguments (INTENT=OUT):
      INTEGER        ICODE                   ! Non zero for failure

      CHARACTER(LEN=256)  CMESSAGE                ! Reason for failure

! Local (dynamic) arrays:
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

! Global variables:
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
!L Description:
!L----------------------------------------------------------------------
!L  COMDECK DOCACP
!L  --------------
!L  DOCACP is list of variables in comdeck COMACP
!L
!L  COMACP contains parameters controlling the assimilation.
!L
!L  AC_OBS_TYPES      - List of AC Obs Types to be processed.
!L                    - Order of processing is stored in LACT.
!L  AC_ORDER          - Order which AC Obs Types MUST be processed.
!L                    - (Coded as group*1000+Obstype)
!L  CSCFACT_H         - Horizontal Correlation Scale Factor.
!L                    - Varies with latitude.
!L  CSCFACT_V         - Vertical Correlation Scale Factor.
!L                    - Varies with level
!L                    - NSLABS_SCFACT(J-1)+1 TO NSLABS_SCFACT(J)
!L  DEF_AC_ORDER      - Default Order and groupings of Obs Types.
!L  DEF_CSCALE_START  - Default Correlation Scale at start of
!L  DEF_CSCALE_OBTIME - insertion period/observation time/end of
!L  DEF_CSCALE_END    - insertion period for each group. Use
!L                    - CSCALE_START/OBTIME/END to change defaults.
!L  DEF_MODE_HANAL    - Default Mode of Horizontal Analysis.
!L  DEF_FI_VAR_FACTOR - Default group dep scaling in FI
!L  DEF_NO_ANAL_LEVS  - Default No of Analysis Levels for each group.
!L  DEF_NO_WT_LEVS    - Default No of Weight Levels for each group.
!                    - Use N_ANAL_LEVS/N_WT_LEVS to change defaults.
!L  DEF_NO_ITERATIONS - Default No of Iterations for groups. Use
!L                    - NO_ITERATIONS in namelist to change defaults.
!L  DEF_INTERVAL_ITER - Default No of Iterations for groups. Use
!L                    - INTERVAL_ITER in namelist to change defaults.
!L  DEF_OBTHIN        - Default ob thinning (use OBTHIN in namelist)
!L                    - (values of 1 imply no thinning, N implies
!L                    - 1/N reports assimilated every Nth step)
!L  DEF_NUDGE_NH      - Default Nudging Coeffs for NH for groups.
!L  DEF_NUDGE_TR      - Default Nudging Coeffs for TR for groups.
!L  DEF_NUDGE_SH      - Default Nudging Coeffs for SH for groups.
!                    - Use NUDGE_NH/TR/SH in namelist to change
!                    - defaults.
!L  DEF_NUDGE_LAM     - Default Nudging Coeffs for LAM for groups.
!                    - Use NUDGE_LAM in namelist to change defaults.
!L  DEF_RADINF        - Default Max Normalised Influence Radius.
!L                    - Use RADINF in namelist to change defaults.
!L  DEF_TGETOBB )     - Default Time Window before/after obs time to
!L  DEF_TGETOBA )     - fetch obs for groups. Use TGETOBB/TGETOBA in
!L                    - namelist to change deafults.
!L  DEF_TIMEB )       - Default Insertion Period before/after obs time
!L  DEF_TIMEA )       - for groups. Use TIMEB/TIMEA in namelist
!L                    - to change defaults.
!L  DF_COEFF          - Coefficient for DIVFILT
!L  DF_SCALE          - DIVFILT scale (metres)
!L  DF_SCALE_LEV      - DIVFILT scale for each level
!L  DIAG_RDOBS        - Diagnostic Control for subroutine RDOBS.
!L  EPSILON_LHN       - Epsilon value for use in LHN
!L  F1_506           \                                                 .
!L  F2_506            } Parameters for 506 ob weight evaluation
!L  F3_506           /
!L  ALPHA_LHN         - Alpha value for use in LHN
!L  LHN_LIMIT         - Limit on + or - Theta incr rate in LHN (K/day)
!L  FI_SCALE_LHN      - Recursive filter scale in m
!L  NPASS_RF_LHN      - Number of passes through filter
!L  FI_SCALE          - FI (Filtered Increments) scale (metres)
!L  FI_SCALE_FACTOR   - FI Scale Factor
!L  GEOWT_H           - Latitude weights for Geostrophic Increments.
!L  GEOWT_V           - Vertical weights for Geostrophic Increments.
!L  GROUP_NO          - Group No of each obs type in LACT.
!L  GROUP_FIRST       - Position in LACT of first type of each group.
!L  GROUP_LAST        - Position in LACT of last  type of each group.
!L  GROUP_INDEX       - Corresponding group in DEF_AC_ORDER for
!L                    - groups in GROUP_NO.
!L  IOMITOBS          - List of Observations not to be assimilated.
!L                    - Use Model Observation Numbers to omit obs.
!L  IUNITNO           - Unit No of Cache file to store obs
!L                    - between timesteps.
!L  L_506_OBERR       - Logical switch to control 506 ob weight calc'n
!L  L_LHN             - Logical switch to perform latent heat nudging
!L  L_LHN_SCALE       - Logical switch to control scaling within LHN
!L  L_LHN_SEARCH      - Logical switch to control use of LHN_SEARCH
!L  L_VERIF_RANGE     - Logical switch to control verification range
!L  L_LHN_LIMIT       - Logical switch to control limiting of increments
!L  L_LHN_FACT        - Logical switch to control limiting by 1/alpha
!L  L_LHN_FILT        - Logical switch to control filtering of incrs
!L  LACT              - List of Obs Types to be processed in order
!L                    - of processing.
!L  LAC_UARS          - Logical switch for UARS assimilation.
!L  LAC_MES           - Logical switch for Mesoscale assimilation.
!L  LCHECK_GRID       - Logical switch to control CHECK_OBS call
!L  LENACT            - No of obs for each type in group fetched
!L                    - this timestep.
!L  LGEO  )           - Logical switches to calculate
!L  LHN_DIAG          - Logical switch for detailed LHN diagnostics
!L  LHN_RANGE         - Max search radius used in LHN_SEARCH. Do not set
!L                      to zero, set L_LHN_SEARCH=.FALSE. instead.
!L  LHYDR )           - Geostrophic/Hydrstatic Increments.
!L  LHYDROL           - Logical switch to calc hydrology incrs.
!L  LRADAR            - Logical array to determine which radars to use
!L  L_LATLON_PRVER    - Logical switch to verify precip in lat/lon area
!L    NORTHLAT        - Co-ords in
!L    SOUTHLAT        -            real lat/lon
!L    WESTLON         -                         for rain
!L    EASTLON         -                                  verif area.
!L  L_MOPS_EQUALS_RH  - If .TRUE. then MOPS cloud obs are
!L                    - rh values (%), else cloud fractions
!L  L_OBS_CHECK       - If .FALSE. then skip check to see if there
!L                    - are any obs to assimilate (non-oper run only)
!L  LSYN              - Logical switch for Synoptic Insertion.
!L  LWBAL_SF          - Controls use of WINDBAL routine for surface wind
!L  LWBAL_UA          - Controls use of WINDBAL routine for uair wind
!L  MASTER_AC_TYPES   - Master list of AC Observation Types
!L                    - known to AC Scheme. See DEF_TYPE.
!L  MACDIAG           - Diagnostics control for each timestep.
!L                    - See AC on use of this array.
!L  MDATADFN          - Mode for Data Density Formula.
!L  MGEOWT            - Mode of Latitude Weighting for
!L                    - Geostrophic Increments.
!L  MGLOSSFN          - GLOSS processing Function Option. (see VANLASS)
!L  MHCORFN           - Correlation Function Option.
!L  MODEACP           - No of timesteps allowed in dimensioning of
!L                    - MACDIAG. Code loops back to MACDIAG(1) if
!L                    - TIMESTEP_NO.GT.MODEACP
!L  MVINT205          - Options for vertical interp    (LASS/GLOSS  )
!L  MRAMPFN           - Mode for Time ramp in HORINF
!L  MWTFN             - Mode for Weights Formula.
!L  NACT              - No of obs types in LACT.
!L  N_GROUPS          - No of groups being processed.
!L  NO_OBS_FILES      - No of observation files to be used.
!L  NO_SCFACT         - List of obs types on which correlation
!L                    - scale factor is not to be applied.
!L  NON_DIV_COR       - Factor for Non-Divergent Correction.
!L  NON_DIV_COR_10M - As NON_DIV_COR but for 10m wind data
!L  NPASS_RF          - Number of passes in RFILT
!L  NPROG             - Number set by individual programmers
!L                    - doing test work. Numbers in use : 1001 (DR)
!L  NRADARS           - No of radars in network
!L  NSLABS_SCFACT     - Slab for each level
!L                    - (Slab is group of levels with same CSCFACT_V)
!L  OBS_FORMAT        - Format of AC Obs file (=2, only one format)
!L  OBS_UNITNO        - Unit No of first AC Obs file (=70)
!L  OBTIME_NOM        - Nominal Observation Time for Synoptic Insertion
!L                    - Mode. Relative time from start of assimilation.
!L  RADAR_LAT         - Coordinates of radars
!L  RADAR_LON         -     "        "    "
!L  RADAR_RANGE       - Max range (km) of reliable radar rain rates
!L  RADAR_RANGE_MAX   - Max. range of radar data used for LHN (km)
!L  RELAX_CF_LHN      - Relaxation coef used for theta incrs in LHN
!L  REMOVE_NEG_LH - Logical switch for removing -ve LH values
!L                                      Mark Dixon 23/02/06
!L  SPEED_LIMIT305    - Min speed of scatwinds for which observed
!L                    - direction is used. (below limit speed only
!L                    - is assimilated.
!L  THRESH_DL         - threshold (mm/hr) between dry/light rain
!L  THRESH_LM         - threshold (mm/hr) between light/moderate rain
!L  THRESH_DL         - threshold (mm/hr) between moderate/heavy rain
!L  THRESH_RMSF       - threshold (mm/hr) for calcn of rms factor score
!L  TIMEF_START       - Start    ) Values for
!L  TIMEF_OBTIME      - Obs time ) Time Factor
!L  TIMEF_END         - End      ) Ramp Function
!L  TROPLAT           - Latitude at which parameters start to change
!L                    - from their mid-latitude to tropical values.
!L  TROPINT           - Interval over which parameters start to change
!L                    - from their mid-latitude to tropical values.
!L  TYPE_INDEX        - Position of obs types in LACT in MASTER_AC_TYPES
!L  USE_CONV_IN_MOPS  - Logical switch for using convection in MOPS
!L  VERT_CUTOFF_SL    - No of scale hts over which single level obs used
!L  VERT_CUTOFF_BW    - as VERT_CUTOFF_SL but for bogus winds
!L  VERT_CUTOFF_BH    - as VERT_CUTOFF_SL but for bogus humidity
!L  VERT_COR_AERO     - Vertical correlation scale for aerosol incrs
!L  VERT_COR_SCALE    - Vertical Correlation Scale Coefficient.
!L                    - calculated from ACP namelist array CSCALE_VERT
!L                      (n,1)=extra tropical temps,(n,2)=tropical temps
!L                      (n,3)=extra tropical winds,(n,4)=tropical winds
!L  VERT_FILT         - Vertical Filtering of Increments
!L                    - from soundings.
!L  WB_THETA_UA       - If T WINDBAL will calc theta incs from upper air
!L                    - winds.
!L  WB_THETA_SF       - If T WINDBAL will calc theta incs from surface
!L                    - winds.
!L  WB_LAT_CC         - horizontal correlation coeff for WINDBAL
!L  WB_VERT_V         - vertical variation in WINDBAL correlations
!L  WB_LAND_FACTOR    - extra scaling factor of WINDBAL inc over land
!L  WB_LAND_SCALE     - apply WB_LAND_FACTOR scaling if true
!L  WB_LonOffset      -) These define a subset of a limited area model
!L  WB_LonPts         -) within which WINDBAL will work. They exist to
!L  WB_LatOffset      -) select a region on which a multigrid Poisson
!L  WB_LatPts         -) solver can be used efficiently. Offsets are
!L                    -) from start of full LAM grid (so offsets of
!L                    -) zero mean no offset). WB_LonPts & WB_LatPts
!L                    -) define the length of the subset in their
!L                    -) respective directions.
!L
!-----------------------------------------------------------------------
! COMAG start

      INTEGER :: DEF_AGRES_ROWS(NOBTYPMX)
      INTEGER :: DEF_AGRES_PTS(NOBTYPMX)
      INTEGER :: NROWSAG
      INTEGER :: NPTSAGMX
      INTEGER :: NPTSAGMN
      INTEGER :: NPTSAG(ROWS_MAX)
      INTEGER :: NPTS0AG(ROWS_MAX+1)
      INTEGER :: MIN_AGPTS

      REAL :: STAGROW1
      REAL :: STAGPT1
      REAL :: ROW1MG
      REAL :: ROW1MGTH
      REAL :: ROW1MGUV
      REAL :: ROW1AG
      REAL :: DLATAG
      REAL :: DLATMG
      REAL :: PT1MGTH
      REAL :: PT1MGUV
      REAL :: PT1AG
      REAL :: DLONGMG
      REAL :: AGLATDEC
      REAL :: AGROWLEN
      REAL :: DLONGAG(ROWS_MAX)
      REAL :: COSROWAG(ROWS_MAX)

      LOGICAL :: LAGNP
      LOGICAL :: LAGSP

      COMMON /COMAG/ DEF_AGRES_ROWS, DEF_AGRES_PTS,                     &
      ! Analysis grid variables
     &  NROWSAG,                                                        &
     &  NPTSAG,   NPTS0AG,  NPTSAGMX,  NPTSAGMN,                        &
     &  MIN_AGPTS,                                                      &
     &  LAGNP,    LAGSP,                                                &
     &  STAGROW1, STAGPT1,                                              &
     &  ROW1AG,   PT1AG,                                                &
     &  DLATAG,   DLONGAG,                                              &
     &  AGLATDEC, AGROWLEN, COSROWAG,                                   &
      ! Model grid variables
     &  ROW1MG,   ROW1MGTH, ROW1MGUV,                                   &
     &  PT1MGTH,  PT1MGUV,                                              &
     &  DLATMG,   DLONGMG

! COMAG end
!-----------------------------------------------------------------------
!L  COMDECK DOCAG
!L  -------------
!L  COMAG contains variables to define the Analysis Grid.
!L
!L  AG = Analysis Grid
!L  MG = Model    Grid
!L
!L  DEF_AGRES_ROWS - Default Ratio of MG Rows to AG Rows for each group.
!L                 - Use AGRES_ROWS in namelist to change defaults.
!L  DEF_AGRES_PTS  - Default Ratio of MG Points to AG Points for each
!L                 - group. (Equatorwards from AGLATDEC).
!L                 - Use AGRES_ROWS and AGRES_PTS in namelist to
!L                 - change defaults.
!L
!L  COLATITUDES & LONGITUDES are in RADIANS, not DEGREES.
!L
!L  Variables defining Analysis Grid.
!L  ---------------------------------
!L
!L  NROWSAG        - No of rows in AG.
!L  NPTSAG         - No of points in each AG row.
!L  NPTS0AG        - Offset to first point in each row in AG.
!L  NPTSAGMX       - Maximum no of points in AG rows.
!L  NPTSAGMN       - Minimum no of points in AG rows.
!L  MIN_AGPTS      - Minimum no of points on any AG ROW.
!L  LAGNP          - Logical set to T if first row of AG = North Pole
!L  LAGSP          - Logical set to T if last  row of AG = South Pole
!L  STAGROW1       - Stagger of first row of MG to first row of AG,
!L                 - (expressed as fraction of AG row spacing).
!L  STAGPT1        - Stagger of first point of MG to first point of AG,
!L                 - (expressed as fraction of AG row spacing).
!L  ROW1AG         - Co-latitude of first row   of AG.
!L  PT1AG          - Longitude   of first point of AG.
!L  DLATAG         - Co-latitude of row spacing of AG.
!L  DLONGAG        - Longitude spacing for each AG row.
!L  AGLATDEC       - Co-latitude at which no of pts in AG rows start
!L                 - to decrease.
!L  AGROWLEN       - Length of row in radians of latitude.
!L  COSROWAG       - Array storing 1/(DLONGAG(x)*cos(lat of row x))
!L                 - where x is the AG row no.
!L
!L  Variables defining Model Grid.
!L  ------------------------------
!L
!L  ROW1MG         - Co-latitude of first row of MG in use.
!L  ROW1MGTH       - Co-latitude of first row of MG for Theta.
!L  ROW1MGUV       - Co-latitude of first row of MG for Winds.
!L  PT1MGTH        - Longitude of first point on MG for Theta.
!L  PT1MGUV        - Longitude of first point on MG for Winds.
!L  DLATMG         - Co-latitude of row spacing of MG (+ve for N->S)
!L  DLONGMG        - Longitude spacing of MG.
!-----------------------------------------------------------------------
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


! These are used in variable resolution runs
      REAL lambda_p(1-halo_i:row_length+halo_i)
      REAL phi_p(1-halo_i:row_length+halo_i, 1-halo_j:p_rows+halo_j)
      LOGICAL L_regular

!  For load balancing
      integer :: lenob_total

      INTEGER :: OBS_FLAG(NOBSMAX)       ! Observation flags
      INTEGER        OBS_NO(NOBSMAX)       ! Observation numbers
      REAL  ::   OBS(TNDVMAX)            ! Observation data
! Local variables:
      INTEGER        TNDV                    ! Total no of data values
                                             ! for obs in assimilation
                                             ! time window
      INTEGER        TNOBS                   ! Total no of obs in
                                             ! assimilation time window
      INTEGER        NAVT                    ! Analysis variable type
                                             ! 1/2/3/4/5/6
                                             ! p*/theta/winds/rh/precip
                                             ! /tracer
      INTEGER        WKLEN                   ! Length of vertical
                                             ! dimension of array for
                                             ! derived increments made
                                             ! by HYDRST, GEOSTR &
                                             ! WINDBAL.
      INTEGER        NPTOBT                  ! Offset to first obs of
                                             ! each type in group
      INTEGER        NO_WT_LEVS              ! No of weight levels
                                             ! for group of obs types
      INTEGER        NO_ANAL_LEVS            ! No of analysis levels
                                             ! for group of obs types
      INTEGER        NO_ANAL_VAR             ! No of variables being
                                             ! analysed (2 for winds)
      INTEGER        LENMG                   ! Length of model grid
      INTEGER        LENAG                   ! Length of analysis grid
      INTEGER        LENOB                   ! No of obs in group
      INTEGER        LENOBT                  ! No of obs for obs type
      INTEGER        ITOTAL                  ! Total no of iterations so
                                             ! far in whole assimilation
      INTEGER        IDIAG                   ! Diagnostic control for
                                             ! this timestep
      INTEGER        TOTAL_NO_ITERS          ! Total no of iterations
                                             ! to be done this timestep
                                             ! (including diagnostic
                                             ! only iterations)
      INTEGER        ITER_NO                 ! Iteration number
                                             ! excluding diagnostic
                                             ! only iterations
      INTEGER        KITER                   ! Iteration no used in
                                             ! diagnostic output.
      INTEGER        IACTF                   ! First obs type in group
      INTEGER        IACTL                   ! Last  obs type in group
      INTEGER        INOBS                   ! No of obs for this type
      INTEGER        INOBSDIM                !  "  " (for dimensioning)
      INTEGER        INC_TYPE                ! Pointer into MODEL_INCR
                                             ! if FI used
      INTEGER        ITNOBS(NOBTYPMX)        ! No of obs in each group
                                             ! assimilated on timestep
      INTEGER        I             !? tracer
      INTEGER        IPT_TRACER    !? tracer ! pointer to one tracer
                                             ! within full array TRACERS
      INTEGER        MODE_HANAL              ! FI or HORINF for this
                                             ! group
      INTEGER        JITER                   ! Loop conter for iteration
      INTEGER        J
      INTEGER        JGROUP,JJGROUP          ! Loop counter over groups
      INTEGER        JACT,JJACT              ! Loop counter over obs
                                             ! types in group
      INTEGER        IGROUP                  ! Group index in group
                                             ! dependent arrays
      INTEGER        HMRMODE                 ! Mode for action in
                                             ! HMRTORH
                                             ! 1 = Convert RH to HMR
                                             ! 2 = Convert HMR to RH
      INTEGER        VMODE                  ! Mode for action in
                                            ! VISTOLV and AEROTOLA
                                            ! 1 = Convert VIS to LOGVIS
                                            ! 2 = Convert LOGVIS to VIS
      INTEGER        ISTAT                  ! status for GCOM

      REAL           ASSM_TIME               ! Assimilation time
                                             ! Relative to start

      LOGICAL        LWIND                   ! Indicator if group
                                             ! has wind obs types
      LOGICAL        DG_THIS                 ! ) Switches to control
      LOGICAL        DG_BETWEEN              ! ) whether diagnostics
      LOGICAL        DG_END                  ! ) required this timestep
                                             ! ) , between iterations
                                             ! ) and at end of timestep
      LOGICAL        DG_ONLY                 ! Indicator if diagnostic
                                             ! only iteration
      LOGICAL        OLD_LDIAGAC             ! Original value of LDIAGAC
                                             ! stored during timestep
      LOGICAL        L_AC_TIMESTEP           ! Switch to control if
                                             ! iteration done this
                                             ! timestep

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Subroutine calls:

! Extra bits:
      DATA ITOTAL/0/
      SAVE ITOTAL

!- End of header

! If no valid ACOB files available then jump straight out
      IF(NO_OBS_FILES.EQ.0) THEN
        WRITE(6,*) '***WARNING: NO VALID ACOB FILES TO ASSIMILATE'
        WRITE(6,*) 'TIMESTEP_NO TIMESTEP ',TIMESTEP_NO,TIMESTEP
        RETURN
      ENDIF

!*
!L 1.0 Call RDOBS to read in observations either from cache file or
!L     from observation files depending on timestep number.
      IF (lhook) CALL dr_hook('AC',zhook_in,zhook_handle)
      CALL RDOBS (NO_OBS_FILES, TIMESTEP_NO, TIMESTEP, OBS, OBS_FLAG,   &
     &               TNDV, TNOBS, P_LEVELS, Q_LEVELS, TNDVMAX, NOBSMAX, &
     &               lambda_p,phi_p,L_regular,P_ROWS,ROW_LENGTH,        &
     &                                                 ICODE, CMESSAGE)

      IF (MYPE == 0 .AND. PRINTSTATUS >= PRSTATUS_NORMAL) THEN
        WRITE(6,*)'RDOBS done, TIMESTEP_NO = ',TIMESTEP_NO
      END IF

! Check for errors (drop out if one has occured):
      IF (ICODE  >   0) GOTO 999

!L 2.0 Set up Diagnostic Iteration Control from MACDIAG(TIMESTEP_NO)
!     MACDIAG = 0   no diagnostics this timestep
!     MACDIAG = 32  to get diagnostics on each iteration this timestep
!     MACDIAG = +64 for extra diagnostic only iteration between
!                   iterations
!     MACDIAG = + 8 for extra diagnostic only iteration at end
!                   of timestep

      IDIAG = MACDIAG(1 + MOD(TIMESTEP_NO -1, MODEACP))

! Remember LDIAGAC for re-use at end of timestep
      OLD_LDIAGAC = LDIAGAC

! Diagnostics required on each iteration this timestep ?
      DG_THIS    = MOD(IDIAG/32, 2)  ==  1

! Extra diagnostic only iteration required between iterations ?
      DG_BETWEEN = DG_THIS .AND. MOD(IDIAG/64, 2)  ==  1

! Extra diagnostic only iteration required at end of timestep ?
      DG_END     = DG_THIS .AND. MOD(IDIAG/8, 2)  ==  1

! Any diagnostics required this timestep ?
      LDIAGAC = LDIAGAC .AND. (DG_THIS .OR. DG_BETWEEN .OR. DG_END)

!L 3.0 Set up iteration control for ac
!     NO_ITERATIONS : no of iterations for each group.
!     INTERVAL_ITER : interval (in timesteps) between iterations
!                   : for each group.

! Determine total no of iterations.
      TOTAL_NO_ITERS = 0

      DO JGROUP = 1, N_GROUPS
        IGROUP = GROUP_INDEX(JGROUP)

        IF (MOD(TIMESTEP_NO -1, DEF_INTERVAL_ITER(IGROUP))  ==  0) THEN
         ! Do iteration(s) from group JGROUP this timestep.
          TOTAL_NO_ITERS = MAX(TOTAL_NO_ITERS,DEF_NO_ITERATIONS(IGROUP))

        END IF
      END DO

      IF (DG_BETWEEN) TOTAL_NO_ITERS = TOTAL_NO_ITERS*2
      IF (DG_END)     TOTAL_NO_ITERS = TOTAL_NO_ITERS+1

!L 4.0 Start loop iterating ac.
      ITER_NO = 0   !  Iteration No excluding Diagnostic only iters.

! Loop over total no of iterations.
      DO JITER =1, TOTAL_NO_ITERS
       ! Diagnostics only this iteration
        DG_ONLY = (DG_BETWEEN .AND. MOD(JITER,2)  ==  1) .OR.           &
     &            (DG_END     .AND. JITER  ==  TOTAL_NO_ITERS)

        IF (.NOT. DG_ONLY) THEN
          ITOTAL = ITOTAL +1
          ITER_NO = ITER_NO +1

        END IF

       ! KITER is used for printed output.
        KITER = ITER_NO

        IF (DG_ONLY) KITER = 0

        DO J = 1, NOBTYPMX
         ! ITNOBS holds a count of obs used this step of each type
          ITNOBS(J) = 0

        END DO
        NAVT =0

!L      4.1 Loop over groups of ob types
        DO JGROUP = 1, N_GROUPS
        JJGROUP=JGROUP
          IGROUP     = GROUP_INDEX(JGROUP)
          L_AC_TIMESTEP = ITER_NO  <=  DEF_NO_ITERATIONS(IGROUP) .AND.  &
     &                MOD(TIMESTEP_NO-1,DEF_INTERVAL_ITER(IGROUP)) == 0

          IF (DG_ONLY .OR. (.NOT.DG_ONLY .AND. L_AC_TIMESTEP) ) THEN

           ! First obs type in group
            IACTF = GROUP_FIRST(JGROUP)

           ! Last obs type in group
            IACTL = GROUP_LAST (JGROUP)

!L          4.1.1 Define analysis variable indiciator for this
!L                  group (NAVT)
           ! NAVT=4 for humidity (cloud), 5 for precip rate
           ! defined by (1st OB type in group)/100
            NAVT  = LACT(IACTF) / 100
            LWIND = NAVT  ==  3

            IF (LDIAGAC.AND.mype == 0) THEN
              PRINT '(/,A,I3,A,I3,A,I2,A,(20I4),/)',                    &
     &              ' AC STEP =', TIMESTEP_NO, ',', KITER,              &
     &              ' STARTING GROUP NO ', JGROUP,                      &
     &              ' , OBS TYPES ', (LACT(J), J = IACTF, IACTL)

            END IF

           ! Convert Model MIXING RATIO to RELATIVE HUMIDITY
  ! DO for all variables (ensures RH preserved when theta changed)
! (except for humidity with 2A cloud microphysics in use and
!  doing a group of MOPS cloud data
!  -assumes obs type 406 will always be in group on its own)


            IF(NAVT /= 4 ) THEN
              HMRMODE = 1
              CALL HMRTORH (HMRMODE, EXNER, PRESSURE, THETA, RH,        &
     &               P_FIELD, P_LEVELS, Q_LEVELS, ICODE, CMESSAGE)
            ENDIF

             ! Check for error - drop out if bad.
              IF (ICODE >  0) GO TO 999

           ! Get Horizontal analysis mode : FI or HORINF?
            MODE_HANAL = DEF_MODE_HANAL(IGROUP)

           ! Get No of Weight/Analysis Levels for this group.
            NO_ANAL_LEVS = DEF_NO_ANAL_LEVS(IGROUP)
            NO_WT_LEVS   = DEF_NO_WT_LEVS  (IGROUP)

           ! Set up work area for derived theta incrs from lhn

            WKLEN = 1

            IF(NAVT == 5 .AND. L_LHN) WKLEN = Q_LEVELS

              NO_ANAL_VAR = 1

           ! Set size of model grid to that used in lower
           ! level routines
              LENMG = P_FIELD

!L          4.1.2  Decide on analysis grid for this group of types.
! parameters for FI method
              LENAG       = 1   ! ANAL_INCR not used in FI method
              INC_TYPE = NO_ANAL_VAR +1


!L          4.1.3 Get a list of obs relevant to this time (getobs)
!L                for each type in the current group of types
           ! ASSM_TIME is time since start in minutes (for GETOBS)
            ASSM_TIME = REAL(TIMESTEP_NO -1) * TIMESTEP / 60.0

            NPTOBT = 0

           ! Loop over observation types in group
            DO JACT = IACTF, IACTL
            JJACT=JACT
             ! LENOBT is no of obs. for type to be used on
             ! this timestep and is determined in GETOBS
              LENOBT = 0

             ! INOBS is total no of observations for type
              INOBS = OBS_INFO % NOBS(JACT)
              INOBSDIM=MAX(INOBS,1)

               ! Check on available workspace in OBS_NO
                IF (NPTOBT+INOBS  <=  NOBSMAX) THEN
                 CALL GETOBS (JJACT, ASSM_TIME, INOBS, OBS_NO(NPTOBT+1),&
     &                   LENOBT, OBS(OBS_INFO % MDISPOBT(JACT)+1),      &
     &                   OBS_FLAG(OBS_INFO % OBS_NO_ST(JACT)+1),        &
     &                   TIMESTEP_NO, DG_ONLY,INOBSDIM,ICODE, CMESSAGE)

      IF (LDIAGAC.AND.mype == 0) THEN
      WRITE(6,*)'AC af GETOBS - LENOBT ',LENOBT
      END IF
                 ! Check for error - drop out if bad
                  IF (ICODE >  0) GO TO 999

                ELSE
                  ICODE = 2
      WRITE(6,*)'TEST FAILED - NPTOBT,INOBS,NOBSMAX ',                  &
     & NPTOBT,INOBS,NOBSMAX

                 ! Drop out of assimilation
                  GOTO 990

                END IF

              LENACT(JACT) = LENOBT
              NPTOBT       = NPTOBT + LENOBT

            END DO   ! End of loop over obs types in group (JACT)

           ! LENOB is now the no of observations in this group
           ! to be used on this timestep
            LENOB = NPTOBT

!
!   HOW MANY OBS (lenob_total) for all pe's?
!
      lenob_total = lenob
      CALL gc_isum(1, nproc, istat, lenob_total)

      IF (LDIAGAC.AND.mype == 0) THEN
        WRITE(6,*)'b4 AC2 - lenob_total,TNDV,LENOB ',    lenob_total,   &
                                                         tndv, lenob
      END IF
            if(lenob_total >  0)then !  Skip if no obs for this group
             ! increment ob counter array
              ITNOBS(JGROUP) = ITNOBS(JGROUP) + lenob_total


             ! Lower level routine of ac to begin here now that
             ! all the array dimensions are known - lenob, lenag,
             ! no_anal_levs
              IF(LWIND) THEN
                WRITE(6,*) 'FATAL ERROR LWIND=T NO LONGER SUPPORTED'
                CMESSAGE='FATAL ERROR LWIND=T NO LONGER SUPPORTED'
                ICODE=101
                RETURN
              ENDIF
              CALL AC2(P_LEVELS, Q_LEVELS, BL_LEVELS,                   &
     &            ROW_LENGTH, P_ROWS,                                   &
     &            P_FIELD, TIMESTEP_NO, ITER_NO,                        &
     &            TIMESTEP, OBS, TNDV, EXNER, PSTAR,                    &
     &            THETA, RH, QCL, QCF,                                  &
     &            CONV_CLD, LS_RAIN, LS_SNOW, CONV_RAIN, CONV_SNOW,     &
     &               D_THETA_DT_CONV,D_THETA_DT_LS,                     &
     &               LAYER_CLOUD,PRESSURE,                              &
     &            RHCRIT,                                               &
     &            OBS_NO, LENOB, NO_ANAL_LEVS, NO_WT_LEVS, NO_ANAL_VAR, &
     &            LENAG, LENMG, WKLEN, INC_TYPE, NAVT, JJGROUP, LWIND,  &
     &            IACTF, IACTL, DG_ONLY, STINDEX, STLIST, LEN_STLIST,   &
     &            SI, SF, STASHWORK, STASH_LEVELS, NUM_STASH_LEVELS,    &
     &            STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,                &
     &            lambda_p,phi_p,L_regular,                             &
     &            ICODE, CMESSAGE)

             ! Check for error - drop out if bad.
              IF (ICODE >  0) GO TO 999

            END IF

           ! Convert RELATIVE HUMIDITY back to MIXING RATIO
          ! (except for MOPS with 2A cld microphys)

            IF(NAVT /= 4 ) THEN
              HMRMODE = 2
              CALL HMRTORH (HMRMODE, EXNER, PRESSURE, THETA, RH,        &
     &               P_FIELD, P_LEVELS, Q_LEVELS, ICODE, CMESSAGE)
            ENDIF

             ! Check for error - drop out if bad
              IF (ICODE  >   0) GO TO 999


            IF (LDIAGAC.AND.mype == 0) THEN
              PRINT '(A,I3,A,I3,A,I4)', ' AC STEP',                     &
     &             TIMESTEP_NO, ',', KITER, ' END OF GROUP NO', JGROUP

            END IF
          END IF
        END DO   ! End of loop over groups (JGROUP)
      END DO   ! End of loop over total no of iterations (JITER)

       IF (mype == 0) THEN
        PRINT *, ' '
        PRINT '(A,I3)',  ' End of AC for time step:',TIMESTEP_NO
        PRINT '(A,(10I8))', ' Group No   ',(J,J=1,N_GROUPS)
        PRINT '(A,(10I8))', ' No of obs  ',(ITNOBS(J),J=1,N_GROUPS)
        PRINT *, ' '

      END IF

      LDIAGAC = OLD_LDIAGAC

! Fix to allow vectorization of JACT loop by removing character
! operations.
 990  CONTINUE
      IF (ICODE  ==  2) THEN
        ICODE    = 1
        CMESSAGE = 'AC : Insufficient space in array OBS_NO.'//         &
     &                                               'Increase NOBSMAX'

      END IF

 999  CONTINUE
      IF (lhook) CALL dr_hook('AC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE AC

!*L  Arguments & declarations:
END MODULE ac_mod
