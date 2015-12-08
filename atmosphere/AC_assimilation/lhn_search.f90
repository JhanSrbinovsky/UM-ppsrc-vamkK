! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE LHN_SEARCH --------------------------------------------
!LL
!LL  Purpose : Search for suitable nearby latent heating profiles for
!LL            use in the LHN sceme.
!LL
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL  Documentation :  FR  WP 171
!LL
!LLEND------------------------------------------------------------------
!L*  ARGUMENTS
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE lhn_search_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE LHN_SEARCH(POINT,NEAR,RANGE                            &
     &                           ,SEARCH,IROWS,ICOLS,L_FIRST_SEARCH     &
     &                           ,TOT_PR,ANAL_PR,P_FIELD,RADIUS,jpts)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Atmos_Max_Sizes
      USE UM_ParParams
      USE comobs_mod, ONLY: nobtypmx
      IMPLICIT NONE
!-----AC common blocks
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

!-----DECLARE VARIABLES
      INTEGER       POINT,NEAR,P_FIELD,RANGE
      INTEGER       SEARCH(4*RANGE*(RANGE+1),2)
      INTEGER       IROWS, ICOLS,jpts
      INTEGER       RADIUS(5)
      REAL          TOT_PR(P_FIELD),ANAL_PR
      LOGICAL       L_FIRST_SEARCH

!-INTENT=IN---------------------------------------------------
!     P_FIELD                    - No. of points
!     POINT                      - Current point
!     RADIUS(5)                  - Diagnostic for breakdown of searches
!     RANGE                      - Search range in grid points
!     SEARCH                     - Template holding relative positions
!                                    of search points from POINT
!     TOT_PR(P_FIELD)            - Model total precip rates
!     ANAL_PR                    - Analysed ppn rate at POINT
!     IROWS                      - Number of rows in the model grid
!     ICOLS                      -    "    " columns  "    "     "
!-INTENT=INOUT-----------------------------------------------
!     L_FIRST_SEARCH             - True if this is the first search
!                                    this timestep.
!-INTENT=OUT-------------------------------------------------
!     NEAR                       - Nearest point with suitable profile
!*-----------------------------------------------------------

! Local arrays and variables
      INTEGER       XPOINT       ! X coord of POINT
      INTEGER       YPOINT       ! Y   "    "   "
      INTEGER       X            ! X coord of point to check
      INTEGER       Y            ! Y   "    "   "    "   "
      INTEGER       P            ! Number of point at (X,Y)
      INTEGER       JR , JN      ! Loop counters
      REAL          RATIO        ! ratio of model to obs ppn
      REAL          BEST         ! closest RATIO to 1
      LOGICAL       L_FOUND      ! false until a suitable point is found

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!
!*----------------------------------------------------------------------
!
!
!C 1.0     Set up search template if this is the first time
!C         round this timestep.
!
      IF (lhook) CALL dr_hook('LHN_SEARCH',zhook_in,zhook_handle)
      IF (L_FIRST_SEARCH) THEN
        DO JR = 1 , RANGE
          DO JN = 1 , 2*JR
            SEARCH(4*(JR-1)*JR+JN,1) = -JR-1+JN
            SEARCH(4*(JR-1)*JR+JN,2) = -JR
            SEARCH(4*(JR-1)*JR+JN+2*JR,1) = JR
            SEARCH(4*(JR-1)*JR+JN+2*JR,2) = -JR-1+JN
            SEARCH(4*(JR-1)*JR+JN+4*JR,1) = JR+1-JN
            SEARCH(4*(JR-1)*JR+JN+4*JR,2) = JR
            SEARCH(4*(JR-1)*JR+JN+6*JR,1) = -JR
            SEARCH(4*(JR-1)*JR+JN+6*JR,2) = JR+1-JN
          ENDDO     ! JN
        ENDDO       ! JR
        L_FIRST_SEARCH=.FALSE.
      ENDIF

!
!C 2.0     Loop round nearest points in ever increasing radius,
!C         looking for best fit to observed rain at POINT
!
!C 2.1     Initialise variables
!
      L_FOUND = .FALSE.
      BEST = 0.0
      RATIO = 0.0
      YPOINT = INT ((POINT-1)/ICOLS) + 1
      XPOINT = POINT - ( (YPOINT-1) * ICOLS )
!
!C 2.2     Loop over possible ranges
!
      DO JR = 1 , RANGE
        IF (.NOT. L_FOUND) THEN
!
!C  If not found one yet then loop over points at this range
!
          DO JN = 1 , 8*JR
!
!C  Make sure search point is within the grid.
!C    for global code this means P from 1 to number of pts,
!C    for limited area code, also need to make sure X within limits.
!
            X = XPOINT + SEARCH( 4*JR*(JR-1) +JN , 1)
            Y = YPOINT + SEARCH( 4*JR*(JR-1) +JN , 2)
            IF ( X  >=  1 .AND. X  <=  ICOLS .AND.                      &
     &           Y  >=  1 .AND. Y  <=  IROWS )      THEN
              P = (Y-1) * ICOLS + X
! count points
              RADIUS(5)=RADIUS(5)+1
!
!C  Test model ppn at point, P
!
              IF ( TOT_PR(P)  >=  (EPSILON_LHN * ANAL_PR) ) THEN
                RATIO = TOT_PR(P) / ANAL_PR
                IF (RATIO  >   1.0) RATIO = 1.0/RATIO
!
!C  Keep record of best match at this range
!
                IF (RATIO  >   BEST) THEN
                  BEST    = RATIO
                  NEAR    = P
                  IF (.NOT.L_FOUND .AND. LHN_DIAG) THEN
                    IF (JR == 1) RADIUS(2)=RADIUS(2)+1
                    IF (JR == 2) RADIUS(3)=RADIUS(3)+1
                    IF (JR >= 3) RADIUS(4)=RADIUS(4)+1
                  ENDIF
                  L_FOUND = .TRUE.
                ENDIF       ! Ratio test
! Diagnostics of where pt is found
              ENDIF         ! Ppn rate test
            ENDIF           ! Within domain test
          ENDDO      ! JN
        ENDIF             ! Found one test
      ENDDO          ! JR
      IF (lhook) CALL dr_hook('LHN_SEARCH',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE LHN_SEARCH
END MODULE lhn_search_mod
