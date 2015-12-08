! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE VERTANL ------------------------------------------------
!LL
!LL  Purpose : High level vertical analysis routine.
!LL            Calls lower level VAN--- routines which are specific
!LL            to particular AC Observation types.
!LL
!LL   Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL   Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!LL
!LL   Called by: AC2
!*L
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE vertanl_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE VERTANL (PSTGS,                                        &
     &                    KACT,IPASS,LENOBT,NDV,OBS_NO,                 &
     &                    MODEL_OBS_TYPE,OBS_LAT,OBS_LONG,OBS,INOBS,    &
     &                    EXNER,PSTAR,THETA,RH,QCL,QCF,                 &
     &                    CONV_CLD,LAYER_CLOUD,PRESSURE,                &
     &                    LS_RAIN,LS_SNOW,CONV_RAIN,CONV_SNOW,          &
     &                    RHCRIT,                                       &
     &                    P_FIELD,WTSFLD,LENWTS,NLEVWT,                 &
     &                    CF1PT,CF2PT,CF3PT,CF4PT,                      &
     &                    NP1PT,NP2PT,NP3PT,NP4PT,                      &
     &                    OBS_INCR,NORMF,LMISSD,                        &
     &                    P_LEVELS,Q_LEVELS,BL_LEVELS,                  &
     &                    ROW_LENGTH,P_ROWS,                            &
     &                    LENOB,NPTOBT,NO_ANAL_LEVS,NO_ANAL_VAR,        &
     &                    ICODE,CMESSAGE)
!
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Atmos_Max_Sizes
      USE UM_ParParams
      USE comobs_mod, ONLY: nobtypmx, obs_info
      USE getob3_mod, ONLY: getob3
      USE vanmops_mixed_phase_mod, ONLY: vanmops_mixed_phase
      USE vanrain_mod, ONLY: vanrain
      IMPLICIT NONE

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
!
      EXTERNAL  TIMER  
!
      INTEGER P_LEVELS,Q_LEVELS,BL_LEVELS,ROW_LENGTH,P_ROWS
      INTEGER P_FIELD,LENWTS,NLEVWT
      INTEGER LENOBT,NDV,INOBS
      INTEGER LENOB,NO_ANAL_LEVS,NO_ANAL_VAR

      INTEGER NP1PT(LENOBT),NP2PT(LENOBT),NP3PT(LENOBT)
      INTEGER NP4PT(LENOBT)
      INTEGER OBS_NO(LENOBT)
      INTEGER MODEL_OBS_TYPE(LENOBT)

      REAL                                                              &
     &               EXNER    (P_FIELD,P_LEVELS),                       &
     &               PRESSURE (P_FIELD,P_LEVELS),                       &
     &               LAYER_CLOUD(P_FIELD,Q_LEVELS),                     &
     &               PSTGS    (P_FIELD),                                &
     &               PSTAR    (P_FIELD),                                &
     &               THETA    (P_FIELD,P_LEVELS),                       &
     &               RH       (P_FIELD,Q_LEVELS),                       &
     &               QCL      (P_FIELD,Q_LEVELS),                       &
     &               QCF      (P_FIELD,Q_LEVELS),                       &
     &               CONV_CLD (P_FIELD,Q_LEVELS),                       &
     &               LS_RAIN(P_FIELD),                                  &
     &               LS_SNOW(P_FIELD),                                  &
     &               CONV_RAIN(P_FIELD),                                &
     &               CONV_SNOW(P_FIELD),                                &
     &               WTSFLD   (LENWTS,NLEVWT),                          &
     &               RHCRIT(Q_LEVELS),                                  &
     &               CF1PT(LENOBT),CF2PT(LENOBT),CF3PT(LENOBT),         &
     &               CF4PT(LENOBT),                                     &
     &               OBS_INCR(LENOB+1,NO_ANAL_LEVS,NO_ANAL_VAR),        &
     &               NORMF(LENOB+1,NO_ANAL_LEVS),                       &
     &               OBS_LAT(LENOBT),OBS_LONG(LENOBT),                  &
     &               OBS(INOBS,*)

      LOGICAL LMISSD (LENOB+1,NO_ANAL_LEVS)
!
      INTEGER KACT, IPASS, KTYPE, NPTOBT, NLEVOB
      INTEGER ICODE
      CHARACTER(LEN=256) CMESSAGE
!
!     DYNAMIC ALLOCATION
!
      REAL           OBDATA(LENOBT,NDV)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!*
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('VERTANL',zhook_in,zhook_handle)
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER ('VERTANL ',3)

!L*** 1       Get obs and error values using GETOB3
      KTYPE  = LACT(KACT)
      NLEVOB = OBS_INFO % NOBLEV(KACT)

      CALL GETOB3 (KACT,OBS,OBDATA,                                     &
     &             OBS_LAT,OBS_LONG,                                    &
     &             LENOBT,INOBS,NDV,OBS_NO,ICODE,CMESSAGE)
      IF (ICODE >  0) GO TO 999
!
!L*** 2      Call appropriate vertical analysis routine VAN---
!     NB. ICODE is not checked here, because only one VAN routine is
!     called followed by a jump to the end of the routine.
!
      IF (KTYPE == 406) THEN
! for MOPS cloud data

         IF (IPASS == 1) THEN

             CALL VANMOPS_MIXED_PHASE (                                 &
     &                    KACT,IPASS,THETA,EXNER,CONV_CLD,              &
     &                    LAYER_CLOUD,PRESSURE,                         &
     &                    RH,QCL,QCF,P_FIELD,NO_ANAL_LEVS,RHCRIT,       &
     &                    OBDATA,CF1PT,CF2PT,CF3PT,CF4PT,               &
     &                    NP1PT,NP2PT,NP3PT,NP4PT,                      &
     &                    OBS_INCR,NORMF,OBS_LAT,OBS_LONG,              &
     &                    OBS_NO,LMISSD,NPTOBT,LENOBT,NDV,LENOB,        &
     &                    NO_ANAL_LEVS,NO_ANAL_VAR,                     &
     &                    P_LEVELS,BL_LEVELS,ICODE,CMESSAGE)
         ELSEIF(IPASS == 2) THEN

             CALL VANMOPS_MIXED_PHASE (                                 &
     &                    KACT,IPASS,THETA,EXNER,CONV_CLD,              &
     &                    LAYER_CLOUD,PRESSURE,                         &
     &                    WTSFLD,QCL,QCF,LENWTS,NLEVWT,RHCRIT,          &
     &                    OBDATA,CF1PT,CF2PT,CF3PT,CF4PT,               &
     &                    NP1PT,NP2PT,NP3PT,NP4PT,                      &
     &                    OBS_INCR,NORMF,OBS_LAT,OBS_LONG,              &
     &                    OBS_NO,LMISSD,NPTOBT,LENOBT,NDV,LENOB,        &
     &                    NO_ANAL_LEVS,NO_ANAL_VAR,                     &
     &                    P_LEVELS,BL_LEVELS,ICODE,CMESSAGE)
         ENDIF


      ELSEIF (KTYPE == 506) THEN
!
        IF (IPASS == 1) THEN

          CALL VANRAIN (KACT,IPASS,LS_RAIN,                             &
     &      LS_SNOW,CONV_RAIN,CONV_SNOW,P_FIELD,                        &
     &      OBDATA,CF1PT,CF2PT,CF3PT,CF4PT,                             &
     &      NP1PT,NP2PT,NP3PT,NP4PT,                                    &
     &      OBS_INCR,NORMF,OBS_LAT,OBS_LONG,                            &
     &      OBS_NO,LMISSD,NPTOBT,                                       &
     &      P_LEVELS,LENOBT,NDV,LENOB,                                  &
     &      NO_ANAL_LEVS,NO_ANAL_VAR,                                   &
     &      ICODE,CMESSAGE)

        ELSEIF(IPASS == 2) THEN

          CALL VANRAIN (KACT,IPASS,WTSFLD,                              &
     &      WTSFLD,WTSFLD,WTSFLD,LENWTS,                                &
     &      OBDATA,CF1PT,CF2PT,CF3PT,CF4PT,                             &
     &      NP1PT,NP2PT,NP3PT,NP4PT,                                    &
     &      OBS_INCR,NORMF,OBS_LAT,OBS_LONG,                            &
     &      OBS_NO,LMISSD,NPTOBT,                                       &
     &      P_LEVELS,LENOBT,NDV,LENOB,                                  &
     &      NO_ANAL_LEVS,NO_ANAL_VAR,                                   &
     &      ICODE,CMESSAGE)

        ENDIF


      ELSE
        ICODE=1
        CMESSAGE = 'VERTANL : Obs Type not known'
        WRITE(6,*) 'VERTANL KTYPE not processed ',KTYPE
      ENDIF
!
 999  CONTINUE

      IF(ICODE >  0) THEN
        WRITE(6,*) 'VERTANL Error code',ICODE
        WRITE(6,*) CMESSAGE
      ENDIF

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER ('VERTANL ',4)
      IF (lhook) CALL dr_hook('VERTANL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE VERTANL
END MODULE vertanl_mod
