! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE DEF_TYPE ----------------------------------------------
!LL
!LL  Purpose : Initialise Observation Type Dependent Arrays
!LL            in comdeck COMACP
!LL
!LL  Programming Standard : UM Doc Paper No 3 ; Version 4 ; 5/2/92
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!
!*L  ARGUMENTS:---------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE def_type_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE DEF_TYPE (ICODE,CMESSAGE)

      USE um_input_control_mod, ONLY : model_domain
      USE domain_params, ONLY: mt_global

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Atmos_Max_Sizes
      USE UM_ParParams
      USE comobs_mod, ONLY: nobtypmx
      IMPLICIT NONE

      INTEGER ICODE           !  Return Code
      CHARACTER(LEN=256) CMESSAGE  !  Error Message


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
! CTRACERA start
!  Vn    Date    Modification History
! 6.1  23/06/04  Prognostic tracers now in section 33, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.
! 6.2  13/07/05  Also increase A_MAX_TRVARS to 150. R Barnes.
! 6.2  10/11/05  UKCA tracers put into section 34, but limited
!                to 150 to allow space there for emissions and
!                diagnostics too.  R Barnes.

      ! First atmospheric tracer (STASH No)
      INTEGER,PARAMETER:: A_TRACER_FIRST = 1
      !First UKCA tracer (STASH No)
      INTEGER,PARAMETER:: A_UKCA_FIRST = 1

      ! Last atmospheric tracer  (STASH No)
      INTEGER,PARAMETER:: A_TRACER_LAST = 150
      !Last UKCA tracer  (STASH No)
      INTEGER,PARAMETER:: A_UKCA_LAST = 150

      ! Maximum number of atmospheric tracers
      INTEGER,PARAMETER:: A_MAX_TRVARS  = 150
      !Maximum number of UKCA tracers
      INTEGER,PARAMETER:: A_MAX_UKCAVARS  = 150

      ! Index to relative position.
      ! A_TR_INDEX(N) gives position in JTRACER for tracer number N.
      ! Set in SET_ATM_POINTERS.
      ! A_TR_INDEX(N) is the position, in the list of tracers
      ! actually present in D1, that tracer number N (in the list
      ! of all tracers selectable from the user interface) occupies,
      ! if it is present.
      ! If tracer number N is absent then A_TR_INDEX(N) is -1.
      ! Similarly for A_UKCA_INDEX.

      INTEGER :: A_TR_INDEX(A_MAX_TRVARS)
      ! A_TR_StashItem is set up in SET_ATM_POINTERS 
      INTEGER :: A_TR_StashItem(A_MAX_TRVARS)

      INTEGER :: A_UKCA_INDEX(A_MAX_UKCAVARS)
      ! UKCA_tr_StashItem is set up in SET_ATM_POINTERS 
      INTEGER :: UKCA_tr_StashItem(A_MAX_UKCAVARS) 

      ! A_TR_LBC_StashItem is set up in INBOUNDA and is only 
      ! referenced if LBC code is active. 
      INTEGER :: A_TR_LBC_StashItem(A_MAX_TRVARS) 
      INTEGER :: A_TR_active_lbc_index(A_MAX_TRVARS) 

      ! UKCA_tr_LBC_StashItem is set up in INBOUNDA and is only 
      ! referenced if LBC code is active. 
      INTEGER :: UKCA_tr_LBC_StashItem(A_MAX_UKCAVARS) 
      INTEGER :: UKCA_tr_active_lbc_index(A_MAX_UKCAVARS)

      COMMON/ATRACER/A_TR_INDEX, A_TR_StashItem,                        &
     &               A_TR_LBC_StashItem, A_TR_active_lbc_index,         &
     &               A_UKCA_INDEX, UKCA_tr_StashItem,                   &
     &               UKCA_tr_LBC_StashItem, UKCA_tr_active_lbc_index

! CTRACERA end

!*L  WORKSPACE USAGE:---------------------------------------------------
!     Local variables
      INTEGER I,J,JOBT  !  Array index, loop counters.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     Arrays initialised for each observation type.
!     MASTER_AC_TYPES   - AC Observation Types known in AC Scheme
!     DEF_TIMEB         - Insertion period before observation time
!     DEF_TIMEA         - Insertion period after  observation time
!     DEF_TGETOBB       - Time Window before obs time to fetch obs
!     DEF_TGETOBA       - Time Window after  obs time to fetch obs
!     DEF_RADINF        - Maximum Normalised Influence Radius.
!     DEF_OBTHIN        - Observation thinning factor
!     DEF_CSCALE_START  - Correlation Scale at start of insertion period
!     DEF_CSCALE_OBTIME -                      observation time
!     DEF_CSCALE_END    -                      end of insertion period.

      IF (lhook) CALL dr_hook('DEF_TYPE',zhook_in,zhook_handle)
      I=0

      IF (model_domain == mt_global) THEN

      I = I + 1              !  Defaults for Type 101 pstar
      MASTER_AC_TYPES(I)   = 101
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 201 sonde temp
      MASTER_AC_TYPES(I)   = 201
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 202 ship surf temp
      MASTER_AC_TYPES(I)   = 202
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 203 airep temp
      MASTER_AC_TYPES(I)   = 203
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 204 synop surf temp
      MASTER_AC_TYPES(I)   = 204
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 200.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 205 LASS temp
      MASTER_AC_TYPES(I)   = 205
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 3
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 206 SATEM temp
      MASTER_AC_TYPES(I)   = 206
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 207 SAT120 temp
      MASTER_AC_TYPES(I)   = 207
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 3
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 208 constrained LASS
      MASTER_AC_TYPES(I)   = 208
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 3
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 209 BOGUS 1000-500
      MASTER_AC_TYPES(I)   = 209
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 211 UARS Temp
      MASTER_AC_TYPES(I)   = 211
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 301 Sonde winds
      MASTER_AC_TYPES(I)   = 301
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 302 Ship surf wind
      MASTER_AC_TYPES(I)   = 302
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 303 Airep wind
      MASTER_AC_TYPES(I)   = 303
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 304 Synop surf wind
      MASTER_AC_TYPES(I)   = 304
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 240.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 305 Scatwind
      MASTER_AC_TYPES(I)   = 305
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 5
      DEF_CSCALE_START(I)  = 240.0
      DEF_CSCALE_OBTIME(I) = 120.0
      DEF_CSCALE_END(I)    = 120.0

      I = I + 1              !  Defaults for Type 306 Drifter winds
      MASTER_AC_TYPES(I)   = 306
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 311 UARS wind
      MASTER_AC_TYPES(I)   = 311
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 200.0
      DEF_CSCALE_END(I)    = 200.0

      I = I + 1              !  Defaults for Type 401 Sonde RH
      MASTER_AC_TYPES(I)   = 401
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 402 Ship surf RH
      MASTER_AC_TYPES(I)   = 402
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 403 Bogus RH
      MASTER_AC_TYPES(I)   = 403
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 404 Synop surf RH
      MASTER_AC_TYPES(I)   = 404
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 200.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 405 LASS RH
      MASTER_AC_TYPES(I)   = 405
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 3
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 406 MOPS RH
      MASTER_AC_TYPES(I)   = 406
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 0.01
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 999.0
      DEF_CSCALE_OBTIME(I) = 999.0
      DEF_CSCALE_END(I)    = 999.0

      I = I + 1              !  Defaults for Type 407 Cloud histograms
      MASTER_AC_TYPES(I)   = 407
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 100.0
      DEF_CSCALE_OBTIME(I) = 100.0
      DEF_CSCALE_END(I)    = 100.0

      I = I + 1              !  Defaults for type 506 MOPS precip
      MASTER_AC_TYPES(I)   = 506
      DEF_TIMEB(I)         = 180.
      DEF_TIMEA(I)         = 180.
      DEF_TGETOBB(I)       = 180.
      DEF_TGETOBA(I)       = 180.
      DEF_RADINF(I)        = 0.01
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 999.0
      DEF_CSCALE_OBTIME(I) = 999.0
      DEF_CSCALE_END(I)    = 999.0

      !  Defaults for Type 6nn (tracers)
      !  Allow for A_MAX_TRVARS (currently 29) possible tracers
      DO J=1,A_MAX_TRVARS
        I = I + 1
        MASTER_AC_TYPES(I) = 600+J
        DEF_TIMEB(I)       = 240.0
        DEF_TIMEA(I)       = 60.0
        DEF_TGETOBB(I)     = 240.0
        DEF_TGETOBA(I)     = 60.0
        DEF_RADINF(I)      = 3.5
        DEF_OBTHIN(I)      = 1
        DEF_CSCALE_START(I) = 360.0
        DEF_CSCALE_OBTIME(I) = 200.0
        DEF_CSCALE_END(I)  = 200.0
      END DO

      I = I + 1              !  Defaults for Type 901 LOG Visibility
      MASTER_AC_TYPES(I)   = 901
      DEF_TIMEB(I)         = 240.0
      DEF_TIMEA(I)         = 60.0
      DEF_TGETOBB(I)       = 240.0
      DEF_TGETOBA(I)       = 60.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 360.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

!  Defaults for Limited Area.
      ELSE

      I = I + 1              !  Defaults for Type 101 pstar
      MASTER_AC_TYPES(I)   = 101
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 225.0
        DEF_CSCALE_OBTIME(I) = 150.0
        DEF_CSCALE_END(I)    = 165.0
        DEF_TIMEB(I)         = 120.0
        DEF_TIMEA(I)         = 24.0
        DEF_TGETOBB(I)       = 120.0
        DEF_TGETOBA(I)       = 24.0
        DEF_RADINF(I)        = 1.75
      ENDIF

      I = I + 1              !  Defaults for Type 201 sonde temp
      MASTER_AC_TYPES(I)   = 201
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 202 ship surf temp
      MASTER_AC_TYPES(I)   = 202
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 225.0
        DEF_CSCALE_OBTIME(I) = 150.0
        DEF_CSCALE_END(I)    = 165.0
      ENDIF

      I = I + 1              !  Defaults for Type 203 airep temp
      MASTER_AC_TYPES(I)   = 203
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 204 synop surf temp
      MASTER_AC_TYPES(I)   = 204
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 150.0
      DEF_CSCALE_OBTIME(I) = 120.0
      DEF_CSCALE_END(I)    = 120.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 105.0
        DEF_CSCALE_OBTIME(I) = 85.0
        DEF_CSCALE_END(I)    = 85.0
        DEF_TIMEB(I)         = 120.0
        DEF_TIMEA(I)         = 24.0
        DEF_TGETOBB(I)       = 120.0
        DEF_TGETOBA(I)       = 24.0
        DEF_RADINF(I)        = 1.75
      ENDIF

      I = I + 1              !  Defaults for Type 205 LASS temp
      MASTER_AC_TYPES(I)   = 205
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 206 SATEM temp
      MASTER_AC_TYPES(I)   = 206
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 207 SAT120 temp
      MASTER_AC_TYPES(I)   = 207
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 208 constrained LASS
      MASTER_AC_TYPES(I)   = 208
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0

      I = I + 1              !  Defaults for Type 209 BOGUS 1000-500
      MASTER_AC_TYPES(I)   = 209
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 211 UARS Temp
      MASTER_AC_TYPES(I)   = 211
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0

      I = I + 1              !  Defaults for Type 301 Sonde winds
      MASTER_AC_TYPES(I)   = 301
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 302 Ship surf wind
      MASTER_AC_TYPES(I)   = 302
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  =  95.0
        DEF_CSCALE_OBTIME(I) =  75.0
        DEF_CSCALE_END(I)    =  75.0
      ENDIF

      I = I + 1              !  Defaults for Type 303 Airep wind
      MASTER_AC_TYPES(I)   = 303
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 120.0
        DEF_CSCALE_OBTIME(I) = 100.0
        DEF_CSCALE_END(I)    = 100.0
      ENDIF

      I = I + 1              !  Defaults for Type 304 Synop surf wind
      MASTER_AC_TYPES(I)   = 304
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 240.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  =  95.0
        DEF_CSCALE_OBTIME(I) =  75.0
        DEF_CSCALE_END(I)    =  75.0
        DEF_TIMEB(I)         = 120.0
        DEF_TIMEA(I)         = 24.0
        DEF_TGETOBB(I)       = 120.0
        DEF_TGETOBA(I)       = 24.0
        DEF_RADINF(I)        = 1.75
      ENDIF

      I = I + 1              !  Defaults for Type 305 Scatwind
      MASTER_AC_TYPES(I)   = 305
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 240.0
      DEF_CSCALE_OBTIME(I) = 120.0
      DEF_CSCALE_END(I)    = 120.0

      I = I + 1              !  Defaults for Type 306 Ship surf wind
      MASTER_AC_TYPES(I)   = 306
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  =  95.0
        DEF_CSCALE_OBTIME(I) =  75.0
        DEF_CSCALE_END(I)    =  75.0
      ENDIF

      I = I + 1              !  Defaults for Type 311 UARS wind
      MASTER_AC_TYPES(I)   = 311
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 3.5
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 180.0
      DEF_CSCALE_END(I)    = 180.0

      I = I + 1              !  Defaults for Type 401 Sonde RH
      MASTER_AC_TYPES(I)   = 401
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 115.0
        DEF_CSCALE_OBTIME(I) = 85.0
        DEF_CSCALE_END(I)    = 85.0
      ENDIF

      I = I + 1              !  Defaults for Type 402 Ship surf RH
      MASTER_AC_TYPES(I)   = 402
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 240.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 225.0
        DEF_CSCALE_OBTIME(I) = 150.0
        DEF_CSCALE_END(I)    = 165.0
      ENDIF

      I = I + 1              !  Defaults for Type 403 Bogus RH
      MASTER_AC_TYPES(I)   = 403
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 115.0
        DEF_CSCALE_OBTIME(I) =  85.0
        DEF_CSCALE_END(I)    =  85.0
      ENDIF

      I = I + 1              !  Defaults for Type 404 Synop surf RH
      MASTER_AC_TYPES(I)   = 404
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 150.0
      DEF_CSCALE_OBTIME(I) = 120.0
      DEF_CSCALE_END(I)    = 120.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 100.0
        DEF_CSCALE_OBTIME(I) =  75.0
        DEF_CSCALE_END(I)    =  75.0
        DEF_TIMEB(I)         = 120.0
        DEF_TIMEA(I)         = 24.0
        DEF_TGETOBB(I)       = 120.0
        DEF_TGETOBA(I)       = 24.0
        DEF_RADINF(I)        = 1.75
      ENDIF

      I = I + 1              !  Defaults for Type 405 LASS RH
      MASTER_AC_TYPES(I)   = 405
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0

      I = I + 1              !  Defaults for Type 406 MOPS RH
      MASTER_AC_TYPES(I)   = 406
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 0.01
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 999.0
      DEF_CSCALE_OBTIME(I) = 999.0
      DEF_CSCALE_END(I)    = 999.0

      I = I + 1              !  Defaults for Type 407 Cloud histograms
      MASTER_AC_TYPES(I)   = 407
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 50.0
      DEF_CSCALE_OBTIME(I) = 50.0
      DEF_CSCALE_END(I)    = 50.0

      I = I + 1              !  Defaults for type 506 MOPS precip
      MASTER_AC_TYPES(I)   = 506
      DEF_TIMEB(I)         = 180.
      DEF_TIMEA(I)         = 180.
      DEF_TGETOBB(I)       = 180.
      DEF_TGETOBA(I)       = 180.
      DEF_RADINF(I)        = 0.01
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 999.0
      DEF_CSCALE_OBTIME(I) = 999.0
      DEF_CSCALE_END(I)    = 999.0

      I = I + 1              !  Defaults for Type 901 LOG Visibility
      MASTER_AC_TYPES(I)   = 901
      DEF_TIMEB(I)         = 150.0
      DEF_TIMEA(I)         = 30.0
      DEF_TGETOBB(I)       = 150.0
      DEF_TGETOBA(I)       = 30.0
      DEF_RADINF(I)        = 2.25
      DEF_OBTHIN(I)        = 1
      DEF_CSCALE_START(I)  = 300.0
      DEF_CSCALE_OBTIME(I) = 150.0
      DEF_CSCALE_END(I)    = 150.0
      IF (LAC_MES) THEN
        DEF_CSCALE_START(I)  = 100.0
        DEF_CSCALE_OBTIME(I) =  75.0
        DEF_CSCALE_END(I)    =  75.0
        DEF_TIMEB(I)         = 120.0
        DEF_TIMEA(I)         = 24.0
        DEF_TGETOBB(I)       = 120.0
        DEF_TGETOBA(I)       = 24.0
        DEF_RADINF(I)        = 1.75
      ENDIF

      !  Defaults for Type 6nn (tracers)
      !  Allow for A_MAX_TRVARS (currently 29) possible tracers
      DO J=1,A_MAX_TRVARS
        I = I + 1
        MASTER_AC_TYPES(I) = 600+J
        DEF_TIMEB(I)       = 150.0
        DEF_TIMEA(I)       = 30.0
        DEF_TGETOBB(I)     = 150.0
        DEF_TGETOBA(I)     = 30.0
        DEF_RADINF(I)      = 3.5
        DEF_OBTHIN(I)      = 1
        DEF_CSCALE_START(I) = 300.0
        DEF_CSCALE_OBTIME(I) = 180.0
        DEF_CSCALE_END(I)  = 180.0
      END DO

      END IF  !  if GLOBAL

!     Modify defaults for UARS assimilation.
      IF (LAC_UARS) THEN
        DO JOBT=1,I
          DEF_CSCALE_START(JOBT)  = 600.0
          DEF_CSCALE_OBTIME(JOBT) = 400.0
          DEF_CSCALE_END(JOBT)    = 400.0
          IF (MASTER_AC_TYPES(JOBT) == 207) THEN
            DEF_OBTHIN(JOBT)      = 2
          ELSE
            DEF_OBTHIN(JOBT)      = 1
          END IF
        ENDDO
      ENDIF

!     Initialise rest of arrays.
      IF (I <  NOBTYPMX) THEN
        DO JOBT = I+1,NOBTYPMX
          MASTER_AC_TYPES(JOBT)   = 0
          DEF_TIMEB(JOBT)         = 0.0
          DEF_TIMEA(JOBT)         = 0.0
          DEF_TGETOBB(JOBT)       = 0.0
          DEF_TGETOBA(JOBT)       = 0.0
          DEF_RADINF(JOBT)        = 0.0
          DEF_OBTHIN(JOBT)        = 0
          DEF_CSCALE_START(JOBT)  = 0.0
          DEF_CSCALE_OBTIME(JOBT) = 0.0
          DEF_CSCALE_END(JOBT)    = 0.0
        ENDDO
      ENDIF

      IF (lhook) CALL dr_hook('DEF_TYPE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE DEF_TYPE
END MODULE def_type_mod
