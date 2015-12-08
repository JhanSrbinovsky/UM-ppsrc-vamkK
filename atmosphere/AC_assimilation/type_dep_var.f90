! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE TYPE_DEP_VAR -------------------------------------------
!LL
!LL  Purpose : Process &ACP Namelist arrays which are
!LL            Observation Type Dependent.
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: AC Assimilation
MODULE type_dep_var_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE TYPE_DEP_VAR(TIMEB,TIMEA,TGETOBB,TGETOBA,RADINF,OBTHIN,&
     &                         CSCALE_START,CSCALE_OBTIME,CSCALE_END,   &
     &                         ICODE,CMESSAGE)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx
      IMPLICIT NONE


      INTEGER TIMEA   (NOBTYPMX)         !IN time window after ob time
      INTEGER TIMEB   (NOBTYPMX)         !IN time window before ob time
      INTEGER TGETOBA (NOBTYPMX)         !IN window for getting obs (a)
      INTEGER TGETOBB (NOBTYPMX)         !IN window for getting obs (b)
      INTEGER RADINF  (NOBTYPMX)         !IN horiz infl radius
      INTEGER OBTHIN  (NOBTYPMX)         !IN obs thinning factors
      INTEGER CSCALE_START  (NOBTYPMX)   !IN horiz cor scale at start
      INTEGER CSCALE_OBTIME (NOBTYPMX)   !IN   "    "    "   "  ob time
      INTEGER CSCALE_END    (NOBTYPMX)   !IN   "    "    "   "   end
      INTEGER ICODE                      !OUT error code and message
      CHARACTER(LEN=256) CMESSAGE

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


!     Local arrays/variables.

      INTEGER ISCFACT(NOBTYPMX)
      INTEGER JOBT,JOBT2,AC_TYPE,J
      INTEGER FIRST_TYPE,LAST_TYPE,N_TYPES
      LOGICAL LFOUND

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


  30  FORMAT (A,F5.2,A,I4)

!L    Process TIMEB
!     -------------
      IF (lhook) CALL dr_hook('TYPE_DEP_VAR',zhook_in,zhook_handle)
      DO JOBT=1,NOBTYPMX
        IF (TIMEB(JOBT) >  0) THEN
          AC_TYPE = TIMEB(JOBT)/1000
          TIMEB(JOBT) = MOD(TIMEB(JOBT),1000)
          IF (TIMEB(JOBT) <  1 .OR. TIMEB(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE = 'TYPE_DEP_VAR : Invalid value given for TIMEB'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_TIMEB(JOBT2)=TIMEB(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, 'TIMEB : Defaults changed to ',                    &
     &                  TIMEB(JOBT),' mins for all obs types.'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_TIMEB(JOBT2)=TIMEB(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'TIMEB : Default changed to ',                 &
     &                TIMEB(JOBT),' mins for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in TIMEB'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process TIMEA
!     -------------
      DO JOBT=1,NOBTYPMX
        IF (TIMEA(JOBT) >  0) THEN
          AC_TYPE = TIMEA(JOBT)/1000
          TIMEA(JOBT) = MOD(TIMEA(JOBT),1000)
          IF (TIMEA(JOBT) <  1 .OR. TIMEA(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE = 'TYPE_DEP_VAR : Invalid value given for TIMEA'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_TIMEA(JOBT2)=TIMEA(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, 'TIMEA : Defaults changed to ',                    &
     &                  TIMEA(JOBT),' mins for all obs types.'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_TIMEA(JOBT2)=TIMEA(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'TIMEA : Default changed to ',                 &
     &                TIMEA(JOBT),' mins for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in TIMEA'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process TGETOBB
!     ---------------
      DO JOBT=1,NOBTYPMX
        IF (TGETOBB(JOBT) >  0) THEN
          AC_TYPE = TGETOBB(JOBT)/1000
          TGETOBB(JOBT) = MOD(TGETOBB(JOBT),1000)
          IF (TGETOBB(JOBT) <  1 .OR. TGETOBB(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE = 'TYPE_DEP_VAR : Invalid value given in TGETOBB'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_TGETOBB(JOBT2)=TGETOBB(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, 'TGETOBB : Default changed to ',                   &
     &                   TGETOBB(JOBT),' mins for all obs types'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_TGETOBB(JOBT2)=TGETOBB(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'TGETOBB : Default changed to ',               &
     &                TGETOBB(JOBT),' mins for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in TGETOBB'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process TGETOBA
!     ---------------
      DO JOBT=1,NOBTYPMX
        IF (TGETOBA(JOBT) >  0) THEN
          AC_TYPE = TGETOBA(JOBT)/1000
          TGETOBA(JOBT) = MOD(TGETOBA(JOBT),1000)
          IF (TGETOBA(JOBT) <  1 .OR. TGETOBA(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE = 'TYPE_DEP_VAR : Invalid value given in TGETOBA'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_TGETOBA(JOBT2)=TGETOBA(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, 'TGETOBA : Default changed to ',                   &
     &                   TGETOBA(JOBT),' mins for all obs types'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_TGETOBA(JOBT2)=TGETOBA(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'TGETOBA : Default changed to ',               &
     &                       TGETOBA(JOBT),' mins for obs type',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in TGETOBA'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process RADINF
!     --------------
      DO JOBT=1,NOBTYPMX
        IF (RADINF(JOBT) >  0) THEN
          AC_TYPE = RADINF(JOBT)/1000
          RADINF(JOBT) = MOD(RADINF(JOBT),1000)
          IF (RADINF(JOBT) <  1 .OR. RADINF(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE = 'TYPE_DEP_VAR : Invalid value given in RADINF'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_RADINF(JOBT2)=RADINF(JOBT)/100.0
            ENDDO
            if(mype == 0)                                               &
     &      PRINT 30, ' RADINF  : Defaults changed to ',                &
     &                  DEF_RADINF(1),' for all obs types'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_RADINF(JOBT2)=RADINF(JOBT)/100.0
                if(mype == 0)                                           &
     &          PRINT 30, ' RADINF  : Default changed to ',             &
     &                      DEF_RADINF(JOBT2),' for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in RADINF'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO


!L    Process OBTHIN
!     --------------
      DO JOBT=1,NOBTYPMX
        IF (OBTHIN(JOBT) >  0) THEN
          AC_TYPE = OBTHIN(JOBT)/1000
          OBTHIN(JOBT) = MOD(OBTHIN(JOBT),1000)
          IF (OBTHIN(JOBT) <  1 .OR. OBTHIN(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE = 'TYPE_DEP_VAR : Invalid value given in OBTHIN'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_OBTHIN(JOBT2)=OBTHIN(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, ' OBTHIN  : Defaults changed to ',                 &
     &                  DEF_OBTHIN(1),' for all obs types'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_OBTHIN(JOBT2)=OBTHIN(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, ' OBTHIN  : Default changed to ',              &
     &                      DEF_OBTHIN(JOBT2),' for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in OBTHIN'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process CSCALE_START
!     --------------------
      DO JOBT=1,NOBTYPMX
        IF (CSCALE_START(JOBT) >  0) THEN
          AC_TYPE = CSCALE_START(JOBT)/1000
          CSCALE_START(JOBT) = MOD(CSCALE_START(JOBT),1000)
          IF (CSCALE_START(JOBT) <  1 .OR.                              &
     &        CSCALE_START(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE =                                                  &
     &      'TYPE_DEP_VAR : Invalid value given in CSCALE_START'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_CSCALE_START(JOBT2)=CSCALE_START(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &               PRINT *, 'CSCALE_START  : Defaults changed to ',   &
     &                   CSCALE_START(JOBT),' km for all obs types.'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_CSCALE_START(JOBT2)=CSCALE_START(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'CSCALE_START  : Default changed to ',         &
     &          CSCALE_START(JOBT),' km for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in CSCALE_START'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process CSCALE_OBTIME
!     ---------------------
      DO JOBT=1,NOBTYPMX
        IF (CSCALE_OBTIME(JOBT) >  0) THEN
          AC_TYPE = CSCALE_OBTIME(JOBT)/1000
          CSCALE_OBTIME(JOBT) = MOD(CSCALE_OBTIME(JOBT),1000)
          IF (CSCALE_OBTIME(JOBT) <  1 .OR.                             &
     &        CSCALE_OBTIME(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE =                                                  &
     &      'TYPE_DEP_VAR : Invalid value given in CSCALE_OBTIME'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_CSCALE_OBTIME(JOBT2)=CSCALE_OBTIME(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, 'CSCALE_OBTIME : Defaults changed to ',            &
     &                  CSCALE_OBTIME(JOBT),' km for all obs types.'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_CSCALE_OBTIME(JOBT2)=CSCALE_OBTIME(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'CSCALE_OBTIME : Default changed to ',         &
     &          CSCALE_OBTIME(JOBT),' km for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid Obs type used in CSCALE_OBTIME'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process CSCALE_END
!     ------------------
      DO JOBT=1,NOBTYPMX
        IF (CSCALE_END(JOBT) >  0) THEN
          AC_TYPE = CSCALE_END(JOBT)/1000
          CSCALE_END(JOBT) = MOD(CSCALE_END(JOBT),1000)
          IF (CSCALE_END(JOBT) <  1 .OR.                                &
     &        CSCALE_END(JOBT) >  999) THEN
            ICODE=1
            CMESSAGE =                                                  &
     &      'TYPE_DEP_VAR : Invalid value given in CSCALE_END'
            GO TO 999
          ENDIF
          IF (AC_TYPE == 999) THEN
!           Use new value for ALL Obs types
            DO JOBT2=1,NOBTYPMX
              DEF_CSCALE_END(JOBT2)=CSCALE_END(JOBT)
            ENDDO
            if(mype == 0)                                               &
     &      PRINT *, 'CSCALE_END    : Defaults changed to ',            &
     &                   CSCALE_END(JOBT),' km for all obs types.'
          ELSE
            LFOUND = .FALSE.
            DO JOBT2=1,NOBTYPMX
              IF (AC_TYPE == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                DEF_CSCALE_END(JOBT2)=CSCALE_END(JOBT)
                if(mype == 0)                                           &
     &          PRINT *, 'CSCALE_END    : Default changed to ',         &
     &          CSCALE_END(JOBT),' km for obs type ',AC_TYPE
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE =                                                &
     &        'TYPE_DEP_VAR : Invalid AC Obs type used in CSCALE_END'
              GO TO 999
            ENDIF
          ENDIF
        ENDIF
      ENDDO

!L    Process NO_SCFACT
!     -----------------
      DO JOBT=1,NOBTYPMX
        ISCFACT(JOBT)   = NO_SCFACT(JOBT)
        NO_SCFACT(JOBT) = 1   !  Apply scale factor to all types.
      ENDDO
      DO JOBT=1,NOBTYPMX
        IF (ISCFACT(JOBT) >  0) THEN
          LFOUND = .FALSE.
          DO JOBT2=1,NOBTYPMX
            IF (ISCFACT(JOBT) == MASTER_AC_TYPES(JOBT2)) THEN
                LFOUND = .TRUE.
                NO_SCFACT(JOBT2) = 0   !  No scale factor for this type
            ENDIF
          ENDDO
          IF (.NOT.LFOUND) THEN
            ICODE=1
            CMESSAGE =                                                  &
     &      'TYPE_DEP_VAR : Invalid Obs type used in NO_SCFACT'
            GO TO 999
          ENDIF
        ENDIF
      ENDDO

      N_TYPES = 0
      DO JOBT=1,NOBTYPMX
        IF (MASTER_AC_TYPES(JOBT) >  0) THEN
          N_TYPES = N_TYPES+1
        ENDIF
      ENDDO

      DO J=1,(N_TYPES+6)/7
        FIRST_TYPE=(J-1)*7+1
        LAST_TYPE =MIN(J*7,N_TYPES)
          if(mype == 0)then
        PRINT *, ' '
        PRINT '(A,7I8)',   ' AC Obs Types      ',                       &
     &   (MASTER_AC_TYPES(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' TIMEB (mins)      ',                       &
     &   (DEF_TIMEB(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' TIMEA (mins)      ',                       &
     &   (DEF_TIMEA(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' TGETOBB (mins)    ',                       &
     &   (DEF_TGETOBB(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' TGETOBA (mins)    ',                       &
     &   (DEF_TGETOBA(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.2)', ' RADINF            ',                       &
     &   (DEF_RADINF(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7I8)', ' OBTHIN            ',                         &
     &   (DEF_OBTHIN(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' CSCALE_START (km) ',                       &
     &   (DEF_CSCALE_START(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' CSCALE_OBTIME (km)',                       &
     &   (DEF_CSCALE_OBTIME(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7F8.1)', ' CSCALE_END (km)   ',                       &
     &   (DEF_CSCALE_END(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
        PRINT '(A,7I8)',   ' NO_SCFACT         ',                       &
     &   (NO_SCFACT(JOBT),JOBT=FIRST_TYPE,LAST_TYPE)
          endif
      ENDDO

!L    Convert Horizontal Correlation Scales from Km to metres.
      DO JOBT=1,NOBTYPMX
        DEF_CSCALE_START(JOBT)  = DEF_CSCALE_START(JOBT)  * 1000.0
        DEF_CSCALE_OBTIME(JOBT) = DEF_CSCALE_OBTIME(JOBT) * 1000.0
        DEF_CSCALE_END(JOBT)    = DEF_CSCALE_END(JOBT)    * 1000.0
      ENDDO

 999  CONTINUE
      IF (lhook) CALL dr_hook('TYPE_DEP_VAR',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE TYPE_DEP_VAR
END MODULE type_dep_var_mod
