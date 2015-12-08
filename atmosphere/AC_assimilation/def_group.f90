! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DEF_GROUP ----------------------------------------------
!LL
!LL  Purpose : Initialise Group Dependent arrays
!LL            in comdeck COMACP and COMAG.
!LL
!LL  For Global runs : Enable defs GLOBAL
!LL
!LL  Project Task : P3
!LL
!LL  Documentation:
!LL
!
!*L  ARGUMENTS :--------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE def_group_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE DEF_GROUP (P_LEVELS,Q_LEVELS,BL_LEVELS,TR_LEVELS,      &
     &                      ICODE,CMESSAGE)

      USE um_input_control_mod, ONLY : model_domain
      USE domain_params, ONLY: mt_global

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Atmos_Max_Sizes
      USE UM_ParParams
      USE comobs_mod, ONLY: nobtypmx
      IMPLICIT NONE

      INTEGER                                                           &
     &   P_LEVELS                                                       &
                       ! IN - No of Model Levels.
     & , Q_LEVELS                                                       &
                       ! IN - No of Model Wet Levels.
     & , BL_LEVELS                                                      &
                       ! IN - No of levels in boundary layer.
     & , TR_LEVELS                                                      &
                       ! IN - No of Tracer levels.
     & , ICODE         ! OUT - Return Code
      CHARACTER(LEN=256) CMESSAGE  ! OUT - Error Message

      
!     AC Comdecks
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
!     UM Comdecks
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

!*L  Workspace Usage ---------------------------------------------------
!     Local array.
      INTEGER AC_GROUPS (NOBTYPMX)  ! Stores default groups of AC Types.
!     Local variables.
      INTEGER JOBT  !  Loop counter.
      INTEGER I     !  group identifier

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     Order of processing and grouping of Observation Types
!     =====================================================

!     Numbers = (Group No*1000) + Obs Type

!L Default groupings are:
!L 1:pstar
!L 2:upper level temperatures
!L 3:surface temperatures
!L 4:upper level winds
!L 5:surface winds
!L 6:upper level RH
!L 7:surface RH
!L 8:MOPS RH
!L 9:Cloud histograms
!L 10:MOPS precip rate/phase
!L 11:Tracers (NB. when the assimilation is actually run, each tracer
!L    must be in a separate group; they are only put together here for
!L    convenience, to define common assimilation parameters.
!L    Normally only a few would be used at once.)

      DATA AC_GROUPS /                                                  &
     &  1101,                                                           &
     &  2201,2203,2205,2206,2207,2208,2209,2211,                        &
     &  3202,3204,                                                      &
     &  4301,4303,4311,                                                 &
     &  5302,5304,5305,5306,                                            &
     &  6401,6403,6405,                                                 &
     &  7402,7404,                                                      &
     &  8406,                                                           &
     &  9407,                                                           &
     &  10506,                                                          &
     &  11601, 11602, 11603, 11604, 11605,                              &
     &  11606, 11607, 11608, 11609, 11610,                              &
     &  11611, 11612, 11613, 11614, 11615,                              &
     &  11616, 11617, 11618, 11619, 11620,                              &
     &  11621, 11622, 11623, 11624, 11625,                              &
     &  11626, 11627, 11628, 11629,                                     &
     &  12901, 126*0/
!     NB. the number of tracers (group 11) should correspond to
!     A_MAX_TRVARS (in CTRACERA), currently 29

      IF (lhook) CALL dr_hook('DEF_GROUP',zhook_in,zhook_handle)

!     Copy above list into COMACP array DEF_AC_ORDER

      DO JOBT=1,NOBTYPMX
        DEF_AC_ORDER(JOBT) = AC_GROUPS(JOBT)
      ENDDO

!     Group Dependent Variables
!     =========================

!     DEF_NO_ANAL_LEVS  No of analysis levels
!     DEF_NO_WT_LEVS    No of weight levels
!     DEF_NO_ITERATIONS No of iterations
!     DEF_INTERVAL_ITER Interval in timesteps between iterations
!     DEF_AGRES_ROWS Ratio of No of rows in Model Grid to Analysis Grid
!     DEF_AGRES_PTS  Ratio of No of pts in Model Grid to Analysis Grid
!     DEF_MODE_HANAL Mode of Horizontal Analysis
!     DEF_FI_VAR_FACTOR group dep scaling factor in FI
!     DEF_NUDGE_NH   Nudging Coefficients for NH
!     DEF_NUDGE_TR   Nudging Coefficients for TR
!     DEF_NUDGE_SH   Nudging Coefficients for SH
!     DEF_NUDGE_LAM  Nudging Coefficients for LAM

      IF (model_domain == mt_global) THEN

      I = 1                             ! pstar
      DEF_NO_ANAL_LEVS(I)  = 1
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.8E-4
      ENDIF

      I = 2                             ! upper level temps
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = P_LEVELS
      DEF_NO_WT_LEVS(I)    = P_LEVELS
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.8E-4
      ENDIF

      I = 3                             ! surf temps
      DEF_NO_ANAL_LEVS(I)  = MAX(BL_LEVELS-2,3)
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.8E-4
      ENDIF

      I = 4                             !upper level winds
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35 * 3.
      DEF_NO_ANAL_LEVS(I)  = P_LEVELS
      DEF_NO_WT_LEVS(I)    = P_LEVELS
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 6.0E-4
        DEF_NUDGE_TR(I) = 6.0E-4
        DEF_NUDGE_SH(I) = 6.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 6.6E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.3E-4
      ENDIF

      I = 5                             ! surf winds (inc scat)
      DEF_NO_ANAL_LEVS(I)  = BL_LEVELS
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35 * 3.
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 6.0E-4
        DEF_NUDGE_TR(I) = 6.0E-4
        DEF_NUDGE_SH(I) = 6.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 6.3E-4
        DEF_NUDGE_TR(I) = 3.8E-4

        DEF_NUDGE_SH(I) = 4.0E-4
      ENDIF

      I = 6                             ! upper level RH
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_WT_LEVS(I)    = MIN(Q_LEVELS,P_LEVELS)
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.5E-4
      ENDIF

      I = 7                             ! surface RH
      DEF_NO_ANAL_LEVS(I)  = MAX(BL_LEVELS-2,3)
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.5E-4
      ENDIF

      I = 8                             ! MOPS RH
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 2
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_WT_LEVS(I)    = MIN(Q_LEVELS,P_LEVELS)
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.5E-4
      ENDIF

      I = 9                             ! Cloud histograms
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 2
      DEF_FI_VAR_FACTOR(I) = 2.0
      DEF_NO_ANAL_LEVS(I)  = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_WT_LEVS(I)    = MIN(Q_LEVELS,P_LEVELS)
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.5E-4
      ENDIF

      I = 10                  ! MOPS precip rate/phase
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 2
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = 1
      DEF_NO_WT_LEVS(I)    = 1
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 1.0E6
        DEF_NUDGE_TR(I) = 1.0E6
        DEF_NUDGE_SH(I) = 1.0E6
      ENDIF

      I = 11                            ! tracers
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = TR_LEVELS
      DEF_NO_WT_LEVS(I)    = TR_LEVELS
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 4.0E-4
        DEF_NUDGE_TR(I) = 4.0E-4
        DEF_NUDGE_SH(I) = 4.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.8E-4
      ENDIF

      I = 12                            ! surface LOG Visibility
      DEF_NO_ANAL_LEVS(I)  = BL_LEVELS
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      IF (LAC_UARS) THEN
        DEF_NUDGE_NH(I) = 3.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.0E-4
      ELSE
        DEF_NUDGE_NH(I) = 5.0E-4
        DEF_NUDGE_TR(I) = 3.0E-4
        DEF_NUDGE_SH(I) = 3.5E-4
      ENDIF

      ELSE

      I = 1                             !pstar
      DEF_NO_ANAL_LEVS(I)  = 1
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 5.0E-4

      I = 2                             !upper level temps
      DEF_NO_ANAL_LEVS(I)  = P_LEVELS
      DEF_NO_WT_LEVS(I)    = P_LEVELS
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 5.0E-4

      I = 3                             !surf temps
      DEF_NO_ANAL_LEVS(I)  = MAX(BL_LEVELS-2,3)
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 5.0E-4
      IF (LAC_MES) THEN
        DEF_NO_ANAL_LEVS(I)  = 6
      ENDIF

      I = 4                             !upper level winds
      DEF_NO_ANAL_LEVS(I)  = P_LEVELS
      DEF_NO_WT_LEVS(I)    = P_LEVELS
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35 * 3.
      DEF_NUDGE_LAM(I)   = 6.6E-4

      I = 5                            !surf winds (inc scat)
      DEF_NO_ANAL_LEVS(I)  = BL_LEVELS
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35 * 3.
      DEF_NUDGE_LAM(I)   = 6.3E-4
      IF (LAC_MES) THEN
        DEF_NO_ANAL_LEVS(I)  = 6
      ENDIF

      I = 6                            !upper level RH
      DEF_NO_ANAL_LEVS(I)  = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_WT_LEVS(I)    = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 5.0E-4

      I = 7                            !surface RH
      DEF_NO_ANAL_LEVS(I)  = MAX(BL_LEVELS-2,3)
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 5.0E-4
      IF (LAC_MES) THEN
        DEF_NO_ANAL_LEVS(I)  = 6
      ENDIF

      I = 8                            !MOPS RH
      DEF_NO_ANAL_LEVS(I)  = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_WT_LEVS(I)    = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 2
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 1.0E-3

      I = 9                             ! Cloud histograms
      DEF_NO_ANAL_LEVS(I)  = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_WT_LEVS(I)    = MIN(Q_LEVELS,P_LEVELS)
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 2
      DEF_FI_VAR_FACTOR(I) = 2.0
      DEF_NUDGE_LAM(I)   = 5.0E-4

      I = 10                  ! MOPS precip rate/phase
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 2
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = 1
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NUDGE_LAM(I)     = 1.0E6

      I = 11                            ! tracers
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NO_ANAL_LEVS(I)  = TR_LEVELS
      DEF_NO_WT_LEVS(I)    = TR_LEVELS
      DEF_NUDGE_LAM(I)  = 5.0E-4

      I = 12                           !surface LOG Visibility
      DEF_NO_ANAL_LEVS(I)  = BL_LEVELS
      DEF_NO_WT_LEVS(I)    = 1
      DEF_NO_ITERATIONS(I) = 1
      DEF_INTERVAL_ITER(I) = 1
      DEF_AGRES_ROWS(I)    = 1
      DEF_AGRES_PTS(I)     = 1
      DEF_MODE_HANAL(I)    = 1
      DEF_FI_VAR_FACTOR(I) = 1.35
      DEF_NUDGE_LAM(I)   = 5.0E-4
      IF (LAC_MES) THEN
        DEF_NO_ANAL_LEVS(I)  = 6
      ENDIF

      END IF  ! if GLOBAL

      DO JOBT=I+1,NOBTYPMX
        DEF_NO_ANAL_LEVS(JOBT)  = IMDI
        DEF_NO_WT_LEVS(JOBT)    = IMDI
        DEF_NO_ITERATIONS(JOBT) = IMDI
        DEF_INTERVAL_ITER(JOBT) = IMDI
        DEF_AGRES_ROWS(JOBT)    = IMDI
        DEF_AGRES_PTS(JOBT)     = IMDI
        DEF_MODE_HANAL(JOBT)    = IMDI
        DEF_FI_VAR_FACTOR(JOBT) = RMDI
        IF (model_domain == mt_global) THEN
        DEF_NUDGE_NH(JOBT) = RMDI
        DEF_NUDGE_TR(JOBT) = RMDI
        DEF_NUDGE_SH(JOBT) = RMDI

        ELSE
        DEF_NUDGE_LAM(JOBT) = RMDI
        END IF
      ENDDO


      IF (lhook) CALL dr_hook('DEF_GROUP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE DEF_GROUP
END MODULE def_group_mod
