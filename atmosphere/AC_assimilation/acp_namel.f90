! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINES ACP_NAMEL,ACDIAG_NAMEL  -------------------------------
!LL
!LL  Purpose : Read in AC Namelist (&ACP) and process.
!LL
!LL ACDIAG_NAMEL:
!LL Set defaults for ACDIAG namelist variables.
!LL                  Read in namelist and process.
!LL
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Logical system components covered:
!LL
!LL  Project Task : P3
!LL
!LL  External documentation:
!LL
!LLEND -------------------------------------------------------------
!LLEND
!*L  Arguments:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE acp_namel_mod

USE missing_data_mod, ONLY: IMDI

IMPLICIT NONE

! Switches for assm mode
LOGICAL :: l_ac = .TRUE.

CONTAINS

      SUBROUTINE ACP_NAMEL (P_LEVELS, Q_LEVELS, BL_LEVELS, TR_LEVELS,   &
     & P_ROWS, U_ROWS,ROW_LEN,                                          &
     & TIMESTEP, ICODE, CMESSAGE)

      USE um_input_control_mod, ONLY : model_domain
      USE domain_params, ONLY: mt_global
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx
      USE def_group_mod, ONLY: def_group
      USE def_type_mod, ONLY: def_type
      USE group_dep_var_mod, ONLY: group_dep_var
      USE type_dep_var_mod, ONLY: type_dep_var
      IMPLICIT NONE

! Imported variables (INTENT=IN):L
      INTEGER     P_LEVELS            ! Total number of levels
      INTEGER     Q_LEVELS            ! Total number of wet levels
      INTEGER     BL_LEVELS           ! total number of boundary layer
                                      ! levels
      INTEGER     TR_LEVELS           ! total number of tracer levels
      INTEGER     P_ROWS              ! Number of rows (for pstar)
      INTEGER     U_ROWS              ! Number of rows (for wind)
      INTEGER     ROW_LEN             ! Length of each row.

      REAL        TIMESTEP            ! timestep in seconds

! Exported variables (INTENT=OUT)
      INTEGER ICODE ! Non zero for failure

      CHARACTER(LEN=256) CMESSAGE          ! Reason for failure
!*
!L---------------------------------------------------------------------
!L UM Comdecks
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
!L---------------------------------------------------------------------
!L AC Comdecks
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
      INTEGER P_ROWS_GLOBAL,U_ROWS_GLOBAL
!L---------------------------------------------------------------------
      REAL    GEOWT_NH
      REAL    GEOWT_SH
      REAL    NUDGE_NH(NOBTYPMX)
      REAL    NUDGE_TR(NOBTYPMX)
      REAL    NUDGE_SH(NOBTYPMX)
      REAL    SCFACT_NH
      REAL    SCFACT_TR
      REAL    SCFACT_SH
      REAL    WB_CC_N                     ! Northern hemisphere
                                          ! corln. coeff. for WINDBAL.
      REAL    WB_CC_S                     ! Southern hemisphere
                                          ! corln. coeff. for WINDBAL.
      REAL    GEOWT_LAM
      REAL    NUDGE_LAM(NOBTYPMX)
      REAL    WB_CC_LAM                   ! WINDBAL correlation

      REAL    GEOWT_VERT(21)
      REAL    CSCALE_VERT(21,4)
      REAL    CSCALE_VERT_G(21,4)
      REAL    DEF_CSCALE_VERT_MES(21,4)
      REAL    SCFACT_VERT(MODEL_LEVELS_MAX)
      REAL    FI_VAR_FACTOR(NOBTYPMX)
      REAL    TEMP_COORD
      REAL    CSCALE_VERT_AERO

      INTEGER AC_ORDER     (NOBTYPMX)
      INTEGER NO_ITERATIONS(NOBTYPMX)
      INTEGER INTERVAL_ITER(NOBTYPMX)
      INTEGER N_ANAL_LEVS  (NOBTYPMX)
      INTEGER N_WT_LEVS    (NOBTYPMX)
      INTEGER AGRES_ROWS   (NOBTYPMX)
      INTEGER AGRES_PTS    (NOBTYPMX)
      INTEGER MODE_HANAL   (NOBTYPMX)
      INTEGER TIMEB        (NOBTYPMX)
      INTEGER TIMEA        (NOBTYPMX)
      INTEGER TGETOBB      (NOBTYPMX)
      INTEGER TGETOBA      (NOBTYPMX)
      INTEGER RADINF       (NOBTYPMX)
      INTEGER OBTHIN       (NOBTYPMX)
      INTEGER CSCALE_START (NOBTYPMX)
      INTEGER CSCALE_OBTIME(NOBTYPMX)
      INTEGER CSCALE_END   (NOBTYPMX)
      INTEGER J,JJ,JLEV,JROW,JTYP,JOBT
      INTEGER NCOUNT
      INTEGER WB_START_LEV                 ! First level for non zero
                                           ! vertical correlation coeff.
      INTEGER WB_END_LEV                   ! Last  level for non zero
                                           ! vertical correlation coeff.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Vertical correlation scale defaults; 21 values at 50hPa intervals
! for extratropics T,tropics T,for extratropics V,tropics V
      DATA CSCALE_VERT_G/                                                 &
     &       12.0, 8.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
     &        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
     &       12.0, 8.0, 8.0, 8.0, 7.0, 6.0, 5.5, 5.0, 4.5, 4.0,         &
     &        3.5, 3.0, 3.0, 3.0, 3.0, 4.0, 3.5, 3.0, 2.5, 1.5, 1.5,    &
     &        6.0, 6.0, 5.0, 5.0, 5.0, 5.0, 4.0, 3.5, 3.5, 3.0,         &
     &        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5,    &
     &        6.5, 6.5, 6.0, 6.0, 5.5, 5.0, 5.0, 4.0, 3.5, 3.5,         &
     &        3.5, 3.5, 3.0, 3.0, 3.0, 3.0, 3.0, 2.5, 2.0, 1.5, 1.5/
! for LAM Tropics and Extra-tropics set the same
      DATA CSCALE_VERT/                                                 &
     &       12.0, 8.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
     &        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
     &       12.0, 8.0, 8.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
     &        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
     &        9.0, 6.0, 5.0, 5.0, 5.0, 5.0, 4.0, 3.5, 3.5, 3.0,         &
     &        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5,    &
     &        9.0, 6.0, 5.0, 5.0, 5.0, 5.0, 4.0, 3.5, 3.5, 3.0,         &
     &        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5/

      DATA DEF_CSCALE_VERT_MES/                                         &
     &       18.0,18.0,18.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
     &        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
     &       18.0,18.0,18.0, 8.0, 7.0, 6.0, 5.0, 4.5, 4.0, 4.0,         &
     &        4.0, 4.5, 5.0, 5.0, 5.5, 6.0, 5.5, 4.0, 3.5, 2.0, 2.0,    &
     &       27.0,27.0,27.0,12.0, 8.0, 6.0, 4.0, 3.5, 3.5, 3.0,         &
     &        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5,    &
     &       27.0,27.0,27.0,12.0, 8.0, 6.0, 4.0, 3.5, 3.5, 3.0,         &
     &        3.0, 2.5, 2.5, 2.5, 2.5, 3.0, 2.5, 2.0, 2.0, 1.5, 1.5/
!     NAMELIST ACP
!     ------------
      NAMELIST /ACP/                                                    &
      AC_ORDER, AC_OBS_TYPES,                                          &
      AGRES_ROWS, AGRES_PTS, MIN_AGPTS, AGLATDEC,                      &
      TIMEA,TIMEB,TGETOBA,TGETOBB,                                     &
      N_ANAL_LEVS,N_WT_LEVS,MODE_HANAL,FI_VAR_FACTOR,                  &
      OBTIME_NOM, LTIMER_AC, LGEO, LHYDR, LSYN, LAC_UARS,              &
      VERT_FILT,NO_ITERATIONS,INTERVAL_ITER,                           &
      IOMITOBS,OBTHIN,SPEED_LIMIT305,MGLOSSFN,                         &
      MHCORFN,MACDIAG,MWTFN,MDATADFN,MRAMPFN,                          &
      FI_SCALE,FI_SCALE_FACTOR,DF_COEFF,DF_SCALE,DF_SCALE_LEV,         &
      NPASS_RF,TIMEF_START,TIMEF_OBTIME,TIMEF_END,CSCALE_VERT_G,       &
      CSCALE_VERT,VERT_CUTOFF_SL,VERT_CUTOFF_BW,VERT_CUTOFF_BH,        &
      CSCALE_START,CSCALE_OBTIME,CSCALE_END,RADINF,                    &
      SCFACT_VERT,NO_SCFACT,                                           &
      MGEOWT,GEOWT_VERT,MVINT205,                                      &
      NON_DIV_COR, LAC_MES,                                            &
      NON_DIV_COR_10M,                                                 &
      OBS_FORMAT,NO_OBS_FILES,DIAG_RDOBS,NPROG,L_OBS_CHECK,            &
      LWBAL_SF, LWBAL_UA, WB_THETA_UA, WB_THETA_SF,                    &
      WB_START_LEV, WB_END_LEV, WB_LAND_SCALE, WB_LAND_FACTOR,         &
      THRESH_DL, THRESH_LM, THRESH_MH, THRESH_RMSF,                    &
      RADAR_RANGE  ,  LRADAR  ,  LHYDROL ,   L_LATLON_PRVER ,          &
      NORTHLAT ,  SOUTHLAT ,  WESTLON ,  EASTLON ,                     &
      LHN_RANGE ,  L_LHN , L_LHN_SCALE ,                               &
      L_LHN_SEARCH , LHN_DIAG ,                                        &
      RADAR_LAT , RADAR_LON , RADAR_RANGE_MAX ,                        &
      EPSILON_LHN , ALPHA_LHN , RELAX_CF_LHN , LHN_LIMIT ,             &
      F1_506 , F2_506 , F3_506 ,                                       &
      L_506_OBERR , L_VERIF_RANGE , L_LHN_LIMIT , L_LHN_FACT ,         &
      L_LHN_FILT , FI_SCALE_LHN , NPASS_RF_LHN ,                       &
      L_MOPS_EQUALS_RH , CSCALE_VERT_AERO, LCHECK_GRID,                &
! Global model terms
      GEOWT_NH,GEOWT_SH,                                               &
      NUDGE_NH,NUDGE_TR,NUDGE_SH,                                      &
      SCFACT_NH,SCFACT_TR,SCFACT_SH,                                   &
      TROPLAT,TROPINT,                                                 &
      WB_CC_N, WB_CC_S,                                                &
! Non-global model terms
      WB_LonOffset, WB_LonPts, WB_LatOffset, WB_LatPts,                &
      GEOWT_LAM, NUDGE_LAM, WB_CC_LAM,                                 &
! common
      REMOVE_NEG_LH, USE_CONV_IN_MOPS

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('ACP_NAMEL',zhook_in,zhook_handle)

      P_ROWS_GLOBAL=glsize(2,fld_type_p)
      U_ROWS_GLOBAL=glsize(2,fld_type_p)-1
!L    1. ACP NAMELIST
!L    ---------------
!L 1.1 Set namelist defaults for variables defined by UI

      DO JOBT = 1, NOBTYPMX
! Group dependent arrays in namelist.
        AC_ORDER(JOBT)       = IMDI

! Obs Type dependent arrays in namelist.
        AC_OBS_TYPES(JOBT)   = 0
      END DO

! Switch  for UARS assimilation.
      LAC_UARS = .FALSE.

! Switch  for Mesoscale UM assimilation.
      LAC_MES = .FALSE.
!  Latent Heat Nudging
      L_LHN = .FALSE.

! Format of ac observation file
! 2 - Current ; 1 - Old Format now redundant
      OBS_FORMAT = 2

! Number of observation files to be used
      NO_OBS_FILES = 2

! Switch for 'non-oper no obs' test in AC1
      L_OBS_CHECK = .TRUE.

!L 1.2 Read in Namelist set by UI
      REWIND 5
      READ (5, ACP, END=120, ERR=120)

! QA Fortran recommends use of IOSTAT but it doesn't work with namelist
! Jump over error trap
      GOTO 121

! NAMELIST read error trap:
120   CONTINUE
        ICODE    = 1
        CMESSAGE = ' ACP_NAMEL : Error reading 1st ACP namelist'
        if(mype == 0) PRINT *, CMESSAGE


       ! Jump to end of routine:
        GOTO 999

! No NAMELIST read error so continue...
121   CONTINUE

!L 1.3 Set namelist defaults for variables defined by user
      DO JOBT = 1, NOBTYPMX
       ! Group dependent arrays in namelist.

        NO_ITERATIONS(JOBT)  = IMDI
        INTERVAL_ITER(JOBT)  = IMDI
        N_ANAL_LEVS(JOBT)    = IMDI
        N_WT_LEVS(JOBT)      = IMDI
        AGRES_ROWS(JOBT)     = IMDI
        AGRES_PTS (JOBT)     = IMDI
        MODE_HANAL(JOBT)     = IMDI
        FI_VAR_FACTOR(JOBT)  = RMDI

        IF (model_domain == mt_global) THEN
        NUDGE_NH(JOBT)       = RMDI
        NUDGE_TR(JOBT)       = RMDI
        NUDGE_SH(JOBT)       = RMDI

        ELSE
        NUDGE_LAM(JOBT)      = RMDI

        END IF

! Obs Type dependent arrays in namelist.
        TIMEB(JOBT)          = 0
        TIMEA(JOBT)          = 0
        TGETOBB(JOBT)        = 0
        TGETOBA(JOBT)        = 0
        RADINF(JOBT)         = 0
        OBTHIN(JOBT)         = 0
        CSCALE_START(JOBT)   = 0
        CSCALE_OBTIME(JOBT)  = 0
        CSCALE_END(JOBT)     = 0
        NO_SCFACT(JOBT)      = 0

      END DO
      NO_SCFACT(1) = 305

      IF (model_domain == mt_global) THEN

! Global version:
! Nominal ob time relative to start of assm for synoptic insertion
        OBTIME_NOM = 300.0

! Correlation Scale Factors (Horizontal)
      IF(LAC_UARS)THEN
        SCFACT_NH = 1.0
        SCFACT_TR = 1.0
        SCFACT_SH = 1.0

      ELSE
        SCFACT_NH = 1.0
        SCFACT_TR = 1.3
        SCFACT_SH = 1.5

      ENDIF

! Scaling for geostrophic increments
      IF(LAC_UARS)THEN
        GEOWT_NH  = 0.7
        GEOWT_SH  = 0.7

      ELSE
        GEOWT_NH  = 0.5
        GEOWT_SH  = 0.7

      END IF

! Latitude at which mid-latitude Analysis Parameters
! begin to change to tropical values.
        TROPLAT=30.0

! Interval (in degrees Lat) over which Relaxation Coefficients
! and Correlation Scale Factor change.
        TROPINT = 10.0

! Latitude at which model and Analysis Grid become identical.
! analysis grid=model grid in MPP mode
        AGLATDEC = 90.0

      ELSE
! Limited Area Version

! Nominal ob time relative to start of assm for synoptic insertion
        OBTIME_NOM = 150.0

! Scaling for geostrophic increments
        GEOWT_LAM = 0.5

! Force Model and Analysis Grid to be identical for Limited Area
        AGLATDEC = 90.0
      END IF

! Parameters for time factor
      IF (LAC_UARS) THEN
        TIMEF_START  = 0.01
        TIMEF_END    = 0.01

      ELSE
        TIMEF_START  = 0.05
        TIMEF_END    = 0.05

      END IF

      TIMEF_OBTIME = 1.0

! Min no. of points on analysis grid row
      MIN_AGPTS = 14

! Geostrophic & hydrostatic incs.
      LGEO  = .TRUE.
      LHYDR = .TRUE.

! Hydrology incs.
      LHYDROL = .TRUE.

! Variable (rh or cloud fraction) used in MOPS obs
      L_MOPS_EQUALS_RH = .FALSE.

! Lat/lon area selection for precip verif
      L_LATLON_PRVER = .FALSE.
! Default lat/lon area when above logical set to .TRUE.
! (this area is larger than current FRONTIERS area)
      NORTHLAT = 60.0
      SOUTHLAT = 45.0
      WESTLON  = -20.0
      EASTLON  = 10.0

! Balanced surface pressure & pot. temp. increments for wind ?
      IF (model_domain == mt_global) THEN
      LWBAL_SF = .TRUE.
      LWBAL_UA = .TRUE.

      ELSE
      LWBAL_SF = .FALSE.
      LWBAL_UA = .FALSE.

      END IF

! Enable calculation of balanced pot. temp. incs. for wind ?
! if true WINDBAL theta incrs are calculated from upper air winds.
      WB_THETA_UA = .TRUE.

! if true WINDBAL theta incrs are calculated from surface winds
! (using HYDRST).
      WB_THETA_SF = .TRUE.

! Initialize variables to define the correlation coefficients for
! WINDBAL.
      IF (model_domain == mt_global) THEN
      WB_CC_N       = 0.70
      WB_CC_S       = 0.80
      ELSE
      WB_CC_LAM     = 0.70

! Initialize variables to define a subset of limited area model on
! which a multigrid can be run efficiently (for use in WINDBAL).
! The LAM values are given (with mesoscale values after the comment).
! N.B. the LAM is 229 E/W by 132 N/S, and the mesoscale model is 92 by
! 92.

      END IF
      WB_START_LEV   = BL_LEVELS - 1   !first non zero level
      WB_END_LEV     = P_LEVELS - 2    !last  non zero level
      WB_LAND_SCALE  = .TRUE.  ! rescale WINDBAL incs over land ?
      WB_LAND_FACTOR = 0.5     ! factor for rescaling   "   "

! Synoptic insertion mode
      LSYN = .FALSE.

! Vertical filtering parameter (used in vrtf)
      VERT_FILT = 0.0

! vertical correlation function coeff for aerosol
      CSCALE_VERT_AERO = 18.0

! Vertical correlation function coeff
! EXP(-(CSCALE_VERT*(LN(P1)-LN(P2)))**2)
! VERT_COR_SCALE is interpolated from  CSCALE_VERT
! which is defined for a regular set of 21 eta levels
! (1 through 0 step 0.05)
! 4 copies for extratropical non-wind (,1);tropical non-wind (,2)
!          for extratropical     wind (,3);tropical     wind (,4)
      IF (LAC_UARS) THEN
       ! Altered to decay more sharply than before
       ! (b in equation 5.1.3 is 9 rather than 3 (old default)
       ! or 16 in tropics - same at all levels)
        DO J = 1, 21
          CSCALE_VERT(J,1) = 3.0
          CSCALE_VERT(J,2) = 4.0
          CSCALE_VERT(J,3) = 3.0
          CSCALE_VERT(J,4) = 4.0
          CSCALE_VERT_G(J,1) = 3.0
          CSCALE_VERT_G(J,2) = 4.0
          CSCALE_VERT_G(J,3) = 3.0
          CSCALE_VERT_G(J,4) = 4.0
        END DO
      ELSE IF (LAC_MES) THEN
       ! use revised mes values
       DO JJ = 1, 4
          DO J = 1, 21
            CSCALE_VERT(J,JJ) = DEF_CSCALE_VERT_MES(J,JJ)
            CSCALE_VERT_G(J,JJ) = DEF_CSCALE_VERT_MES(J,JJ)

          END DO
        END DO
!     ELSE for LAM and GL (oper) use DATA table
      END IF

! Cut off in vertical correlation function
! single lvl obs(203,303,403) infl only vert_cutoff_sl scale hts
! same cut-off correl. applies to extrap'n of incomplete soundings.
      IF (LAC_UARS) THEN
        VERT_CUTOFF_SL = 1.0

      ELSE
        VERT_CUTOFF_SL = 2.0

      END IF

! Cut off in vertical correlation function for bogus wind data
      IF (LAC_UARS) THEN
        VERT_CUTOFF_BW = 1.0

      ELSE
        VERT_CUTOFF_BW = 1.0

      END IF

! Cut off in vertical correlation function for bogus humidity data
      VERT_CUTOFF_BH = 0.10

! Correlation function option
!     MHCORFN = 3 For auto-regressive function
!             = 4 with mod which tends to zero at edge on infl
!               +16 For time factor option in INITAC & AC
      IF (LAC_MES) THEN
        MHCORFN = 3

      ELSE
        MHCORFN = 4

      END IF

! Set ac diagnostic mode (each timestep)
! +8 FOR DIAG ONLY ITERATION AT END OF ITERATIONS
! +32 FOR DIAGS ON (OVERRULING LDIAGAC)
! +64 FOR DIAG ONLY ITERATION BETWEEN ITERATIONS
      DO J = 1, MODEACP
        MACDIAG(J) = 0

      END DO

      IF (LAC_MES) THEN
       ! Set MACDIAG to get diagnostics on timesteps 1, 12, 24, 36.
        MACDIAG(1)   = 32
        MACDIAG(12)  = 32
        MACDIAG(24)  = 32
        MACDIAG(36)  = 32

      ELSE  !  global and lam
       ! Set MACDIAG to get diagnostics on timesteps 1, 18 and 36.
        MACDIAG(1)  = 32
        MACDIAG(18) = 32
        MACDIAG(36) = 32

      END IF

! set mode for weights formula
      MWTFN = 2

! set mode for time ramp (1=linear,2=quadratic)
      IF (LAC_UARS) THEN
        MRAMPFN = 1

      ELSE
        MRAMPFN = 2

      END IF

! set mode for data density formula
      MDATADFN = 2

! set mode for GLOSS (type 208) processing option
! 1 as 205,2 use GLOSS error ratio,
! 3 do constrained retrieval, 4 as 3 without a background correction
      MGLOSSFN = 1

!  Min speed limit for ERS-1 winds above which direction is believed
      SPEED_LIMIT305 = 4.0

! set parameters to define LASS/GLOSS vertical interpolation
      MVINT205 = 2     ! 1 for V_INT_T, 2 for V_INT_TP

! set mode for latiude weighting of geostrophic increments
      MGEOWT = 3

!     VERTICAL SCALING FACTOR FOR GEOSTROPHIC INCREMENTS
!     this is defined for a regular set of 21 eta levels
!       (1 through 0 step 0.05)
      IF (LAC_UARS) THEN
! UARS defaults give full weights, except near jet level
        DO J=1,21
          GEOWT_VERT(J) = 1.0
        END DO

        GEOWT_VERT(13)=0.8
        GEOWT_VERT(14)=0.6
        GEOWT_VERT(15)=0.4
        GEOWT_VERT(16)=0.3
        GEOWT_VERT(17)=0.4
        GEOWT_VERT(18)=0.6
        GEOWT_VERT(19)=0.8
      ELSE
! the defaults here give full wt below 500mb and zero wt above 250mb
        DO J=1,10
          GEOWT_VERT(J) = 1.0
        END DO

        GEOWT_VERT(11) = 0.820
        GEOWT_VERT(12) = 0.720
        GEOWT_VERT(13) = 0.620
        GEOWT_VERT(14) = 0.500
        GEOWT_VERT(15) = 0.367
        GEOWT_VERT(16) = 0.200

        DO J = 17, 21
          GEOWT_VERT(J) = 0.0
        END DO
      END IF

! Vertical scaling for horizontal correlation scale
      IF (LAC_MES) THEN
        DO J = 1, BL_LEVELS
          SCFACT_VERT(J)   = 0.875
          NSLABS_SCFACT(J) = 0
          CSCFACT_V(J)     = 0.0
        END DO
        DO J = BL_LEVELS+1, MODEL_LEVELS_MAX
          SCFACT_VERT(J)   = 1.0
          NSLABS_SCFACT(J) = 0
          CSCFACT_V(J)     = 0.0
        END DO

      ELSE
        DO J = 1, MODEL_LEVELS_MAX
          SCFACT_VERT(J)   = 1.0
          NSLABS_SCFACT(J) = 0
          CSCFACT_V(J)     = 0.0
        END DO

      ENDIF

! Defaults for Filtered Increment Mode of Horizontal Analysis
! For testing, FI_SCALE_FACTOR=1.0 will match current HORINF
! Eventually, higher levels can be given a larger scale.
! DF_COEFF=1.0 is initial estimate. DF_COEFF=0.0 turns off DIVFILT
! NPASS_RF=2 matches SOAR and is default.

      NPASS_RF = 2
      IF (LAC_MES) THEN
!       For MOPS data, use FI method without filtering on model grid.
!       (Assumes that no other mes data groups require standard FI)
        NPASS_RF = 0
      END IF

! Only MOPS data are assimilated now so that LAC_MES is T and 
! NPASS_RF is 0. If this is not the case abort
      IF(NPASS_RF.NE.0) THEN

         WRITE(6,*) 'NPASS_RF NOT ZERO'
         WRITE(6,*) 'Check that MOPS data assimilation is on in umui'
         CMESSAGE='MOPS data not checked in umui'
         ICODE=1
         RETURN
      ENDIF

      FI_SCALE = 400000.0
      DF_SCALE = 400000.0

      DO JLEV = 1, MODEL_LEVELS_MAX
        FI_SCALE_FACTOR(JLEV) = 1.0
        DF_COEFF(JLEV)        = 1.0
        DF_SCALE_LEV(JLEV)    = DF_SCALE * FI_SCALE_FACTOR(JLEV)

      END DO

! Factor for non-divergence correction term in correlations
      NON_DIV_COR = 0.8
!  same factor for 10m wind data
      IF(LAC_MES) THEN
        NON_DIV_COR_10M = 0.0
      ELSE
        NON_DIV_COR_10M = 0.8
      ENDIF

! Control of diagnostic output from RDOBS/RDOBS2
! 0 - No listing ; 1 - Standard Listing ; 2 - More detailed listing.
      DIAG_RDOBS = 1

! Thresholds (mm/hr) for rainfall verification
      THRESH_DL   = 0.03
      THRESH_LM   = 0.125
      THRESH_MH   = 0.5
      THRESH_RMSF = 0.125

! Radar data reliability ranges.
! RADAR_RANGE is limit of high reliability (used in LHN & HCS)
! RADAR_RANGE_MAX is limit of usability (used in LHN)
      RADAR_RANGE = 100.0
!  Max Radar range
      RADAR_RANGE_MAX = 200.0
!  --- SPECIFY RADAR LOCATIONS ------
      RADAR_LAT(1)   =  58.21
      RADAR_LON(1)   =  353.81      !  BEACON HILL

      RADAR_LAT(2)   =  57.43
      RADAR_LON(2)   =  357.97      !  HILL OF DUDWICK

      RADAR_LAT(3)   =  55.69
      RADAR_LON(3)   =  355.77      !  CORSE HILL

      RADAR_LAT(4)   =  54.50
      RADAR_LON(4)   =  353.65      !  CASTOR BAY

      RADAR_LAT(5)   =  53.43
      RADAR_LON(5)   =  353.75      !  DUBLIN

      RADAR_LAT(6)   =  52.69
      RADAR_LON(6)   =  351.07      !  SHANNON

      RADAR_LAT(7)   =  53.75
      RADAR_LON(7)   =  357.72      !  HAMELDON

      RADAR_LAT(8)   =  53.33
      RADAR_LON(8)   =  359.45      !  INGHAM

      RADAR_LAT(9)   =  52.40
      RADAR_LON(9)   =  357.40      !  CLEE HILL

      RADAR_LAT(10)  =  51.69
      RADAR_LON(10)  =  359.47      !  CHENIES

      RADAR_LAT(11)  =  51.98
      RADAR_LON(11)  =  355.55      !  CRUGYGORLLWYN

      RADAR_LAT(12)  =  50.96
      RADAR_LON(12)  =  356.55      !  COBBACOMBE

      RADAR_LAT(13)  =  50.82
      RADAR_LON(13)  =  357.44      !  WARDON HILL

      RADAR_LAT(14)  =  50.00
      RADAR_LON(14)  =  354.77      !  PREDANNACK

      RADAR_LAT(15)  =  49.18
      RADAR_LON(15)  =  357.78      !  JERSEY

! List of which radars to use
! (.TRUE. implies MOPS precip data close to that radar will be used
! in rainfall verification and hydrology correction)
      LRADAR(1)   = .TRUE.     !  Beacon Hill
      LRADAR(2)   = .TRUE.     !  Hill of Dudwick
      LRADAR(3)   = .TRUE.     !  Corse Hill
      LRADAR(4)   = .TRUE.     !  Castor Bay
      LRADAR(5)   = .TRUE.     !  Dublin
      LRADAR(6)   = .TRUE.     !  Shannon
      LRADAR(7)   = .TRUE.     !  Hameldon
      LRADAR(8)   = .TRUE.     !  Ingham
      LRADAR(9)   = .TRUE.     !  Clee
      LRADAR(10)  = .TRUE.     !  Chenies
      LRADAR(11)  = .TRUE.     !  Crugygorllwyn
      LRADAR(12)  = .TRUE.     !  Cobbacombe
      LRADAR(13)  = .TRUE.     !  Wardon Hill
      LRADAR(14)  = .TRUE.     !  Predannack
      LRADAR(15)  = .TRUE.     !  Jersey


!  506 observation error calculations (WP171, eqn10)
      L_506_OBERR = .FALSE.
      F1_506 = 10.0
      F2_506 = 1.0
      F3_506 = 0.0
!  Scale points by 1/EPSILON, in LHN, if no near rain found
      L_LHN_SCALE = .TRUE.
!  Use search routine for nearby rain in LHN
      L_LHN_SEARCH = .TRUE.
!  Display detailed diagnostics for LHN routine
      LHN_DIAG = .TRUE.
!  Max range (in grid points) to search to in Latent Heat Nudging code
      LHN_RANGE = 4
!  Verify up to reliable range, or max range
      L_VERIF_RANGE = .TRUE.
!  Minimum acceptable ratio of model ppn to observed, for LHN
      EPSILON_LHN = 0.333
!  Relaxation coefficient for theta increments from the LHN scheme
      RELAX_CF_LHN = 1.0
!  Minimum factor by which model ppn can be scaled to fit obs ppn
      ALPHA_LHN = 0.333
!  Set default for switch to limit size of increments
      L_LHN_LIMIT = .FALSE.
!  Maximum size of increment to theta within LHN (K/day)
      LHN_LIMIT = 1.0
!  Set default for switch to use limit in ALPHA_LHN
      L_LHN_FACT = .TRUE.
!  Set default for switch to filter LHN theta incrs
      L_LHN_FILT = .TRUE.
!  Recursive filter scale in metres
      FI_SCALE_LHN = 17000.0
!  Number of passes through recursive filter, for LHN increments
      NPASS_RF_LHN = 2
! List of observations to be omitted from assimilation
! ** Specify model observation type to omit obs **
      DO JTYP = 1, NANALTYP
        IOMITOBS(JTYP) = 0

      END DO

! Diagnostics options for programmers
      NPROG = 0

! Set timing switch
      LTIMER_AC = .FALSE.

! Set control for calling CHECK_OBS
      LCHECK_GRID = .FALSE.

! Set REMOVE_NEG_LH to false
      REMOVE_NEG_LH =  .FALSE.

! Set USE_CONV_IN_MOPS to true
      USE_CONV_IN_MOPS = .TRUE.


!L 1.4 Read in Namelist set by user
      READ (5, ACP, END=140, ERR=140)

! QA Fortran recommends use of IOSTAT but it doesn't work with namelist
! Jump over NAMELIST read error trap:
      GOTO 141

! Trap NAMELIST read error
140   CONTINUE
        ICODE    = 1
        CMESSAGE = ' ACP_NAMEL : Error reading 2nd ACP namelist'
        if(mype == 0) PRINT *, CMESSAGE

       ! Jump to end of routine
        GOTO 999

! No NAMELIST read error so continue...
141   CONTINUE

!L 1.5 Print out Namelist for this run
      if(mype == 0)then
       PRINT *, ' '
       PRINT *, ' Parameters in Namelist ACP for this run'
       WRITE(6,ACP)
      endif

!L 2   Check Validity of Namelist Parameters
! Check AGLATDEC
      IF(AGLATDEC <  89.)THEN
          AGLATDEC = 90.
          if(mype == 0)                                                 &
     &     PRINT *, ' AGLATDEC reset to 90. when running MPP '
      END IF
! Check that AC_OBS_TYPES has any obs types
      NCOUNT = 0

      DO JOBT = 1, NOBTYPMX
        IF (AC_OBS_TYPES(JOBT)  >   0) THEN
          NCOUNT = NCOUNT + 1

        END IF
      END DO

      IF (NCOUNT  ==  0) THEN
        ICODE    = 1
        CMESSAGE = ' ACP_NAMEL : No obs types in AC_OBS_TYPES.'
        GOTO 999

      END IF

! Convert any old type numbers to new type numbers.
      DO JOBT = 1, NOBTYPMX
        IF (AC_OBS_TYPES(JOBT)  ==  501) THEN
          AC_OBS_TYPES(JOBT) = 302
          if(mype == 0)                                                 &
     &     PRINT *, 'Type 501 in AC_OBS_TYPES changed to Type 302'

        END IF

        IF (AC_OBS_TYPES(JOBT) == 502) THEN
          AC_OBS_TYPES(JOBT) = 305
          if(mype == 0)                                                 &
     &     PRINT *, 'Type 502 in AC_OBS_TYPES changed to Type 305'

        END IF
      END DO

      IF (model_domain == mt_global) THEN
! Check that troplat and tropint are correctly specified
      IF (TROPLAT + TROPINT  >   90.0) THEN
        ICODE    = 1
        CMESSAGE =                                                      &
     &          ' ACP_NL : TROPLAT/TROPINT bug, is TROPLAT-TROPINT >= 0'
        GOTO 999

      END IF

! Check geowt_nh non zero if geostrophic increments are to be used
      IF (GEOWT_NH  ==  0.0 .AND. (LGEO)) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : GEOWT_NH = 0 not allowed with LGEO = T'
        GOTO 999

      END IF

      ELSE
! Check geowt_lam non zero if geostrophic increments are to be used
      IF (GEOWT_LAM  ==  0.0 .AND. (LGEO)) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : GEOWT_LAM=0 not allowed with LGEO = T'
        GOTO 999

      END IF
      END IF  ! if GLOBAL

! Check windbal variables
      IF (model_domain == mt_global) THEN
      IF (WB_CC_N  <   0.0 .OR. WB_CC_N  >   1.0) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid WB_CC_N'
        GOTO 999

      ENDIF
      IF (WB_CC_S  <   0.0 .OR. WB_CC_S  >   1.0) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid WB_CC_S'
        GOTO 999

      ENDIF
      ELSE
      IF (WB_CC_LAM  <   0.0 .OR. WB_CC_LAM  >   1.0) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid WB_CC_LAM'
        GOTO 999

      ENDIF
      END IF  ! if GLOBAL

! Check MRAMPFN
      IF (MRAMPFN  /=  1 .AND. MRAMPFN  /=  2) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid MRAMPFN'
        GOTO 999

      ENDIF

! Check MDATADFN
      IF (MDATADFN  /=  1 .AND. MDATADFN  /=  2) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid MDATADFN'
        GOTO 999

      ENDIF

! Check MGLOSSFN
      IF (MGLOSSFN  <   1 .OR. MGLOSSFN  >   4) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid MGLOSSFN'
        GOTO 999

      ENDIF

! Check MWTFN
      IF (MWTFN /= 1 .AND. MWTFN /= 2) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid MWTFN'
        GOTO 999

      ENDIF

! Check MHCORFN
      IF (MOD(MHCORFN,8) <  1 .OR. MOD(MHCORFN,8) >  4) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid MHCORFN'
        GOTO 999

      END IF

! Check MVINT205
      IF (MVINT205  /=  1 .AND. MVINT205  /=  2) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid MVINT205'
        GOTO 999

      ENDIF

! Check NO_OBS_FILES
      IF (NO_OBS_FILES >  10) THEN
        ICODE    = 1
        CMESSAGE = 'ACP_NL : invalid NO_OBS_FILES (max allowed is 10)'
        GOTO 999
      END IF

! Check OBS_FORMAT - only allow portable UM format
      IF (OBS_FORMAT /= 2 .AND. OBS_FORMAT /= 3) THEN
       ICODE    = 1
       CMESSAGE = 'ACP_NL : OBS_FORMAT not equal to 2 or 3'
       GOTO 999
      END IF

!L 3. Set Defaults for Group Dependent arrays
      CALL DEF_GROUP (P_LEVELS, Q_LEVELS, BL_LEVELS, TR_LEVELS,         &
     & ICODE, CMESSAGE)

      IF (ICODE >  0) GO TO 999

!L 4. Set Defaults for Obs Type Dependent arrays
      CALL DEF_TYPE (ICODE,CMESSAGE)

      IF (ICODE >  0) GO TO 999

!L 5. Check/Process Group dependent variables
      CALL GROUP_DEP_VAR (AC_ORDER, NO_ITERATIONS, INTERVAL_ITER,       &
     &     N_ANAL_LEVS, N_WT_LEVS,                                      &
     &     NUDGE_NH, NUDGE_TR, NUDGE_SH,                                &
     &     NUDGE_LAM,                                                   &
     &     AGRES_ROWS, AGRES_PTS, MODE_HANAL, FI_VAR_FACTOR,            &
     &     ICODE,CMESSAGE)

      IF (ICODE >  0) GO TO 999

!L 6. Check/Process Obs Type dependent variables
      CALL TYPE_DEP_VAR (TIMEB,TIMEA,TGETOBB,TGETOBA,RADINF,OBTHIN,     &
     &                   CSCALE_START,CSCALE_OBTIME,CSCALE_END,         &
     &                   ICODE,CMESSAGE)

      IF (ICODE >  0) GO TO 999


!L 11. Convert rainfall thresholds from mm/hr to kgm-2s-1
       THRESH_DL   = THRESH_DL   / 3600
       THRESH_LM   = THRESH_LM   / 3600
       THRESH_MH   = THRESH_MH   / 3600
       THRESH_RMSF = THRESH_RMSF / 3600

!L  12. Check Lat/Lon verification coordinates

! Need longitudes in range :  -180 < LON <= 180
      IF (EASTLON  >   180.0) THEN
        EASTLON = EASTLON - 360.0
        if(mype == 0) PRINT *,                                          &
     &    "EASTLON changed from ",EASTLON+360.0," to ",EASTLON
      ENDIF

      IF (WESTLON  >   180.0) THEN
        WESTLON = WESTLON - 360.0
        if(mype == 0) PRINT *,                                          &
     &    "WESTLON changed from ",WESTLON+360.0," to ",WESTLON
      ENDIF

      IF (NORTHLAT  <   SOUTHLAT) THEN
        if(mype == 0) PRINT *,                                          &
     &    "Swapping NORTHLAT and SOUTHLAT so that North>South"
        TEMP_COORD = NORTHLAT
        NORTHLAT   = SOUTHLAT
        SOUTHLAT   = TEMP_COORD
      ENDIF

      IF (EASTLON  <   WESTLON) THEN
        if(mype == 0) PRINT *,                                          &
     &   "Swapping EASTLON and WESTLON so that East > West"
        TEMP_COORD = EASTLON
        EASTLON    = WESTLON
        WESTLON    = TEMP_COORD
      ENDIF

 999  CONTINUE
      IF (lhook) CALL dr_hook('ACP_NAMEL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ACP_NAMEL
!*L  Arguments:---------------------------------------------------
END MODULE acp_namel_mod
