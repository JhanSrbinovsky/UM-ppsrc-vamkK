! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE GROUP_DEP_VAR ------------------------------------------
!LL
!LL  Purpose : Process &ACP Namelist arrays which are Group Dependent.
!LL
!LL
!LL
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!*L  Arguments:---------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE group_dep_var_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE GROUP_DEP_VAR (AC_ORDER,NO_ITERATIONS,INTERVAL_ITER,   &
     &                          N_ANAL_LEVS,N_WT_LEVS,                  &
     &                          NUDGE_NH,NUDGE_TR,NUDGE_SH,             &
     &                          NUDGE_LAM,                              &
     &                          AGRES_ROWS,AGRES_PTS,                   &
     &                          MODE_HANAL,FI_VAR_FACTOR,               &
     &                          ICODE,CMESSAGE)

      USE um_input_control_mod, ONLY : model_domain
      USE domain_params, ONLY: mt_global

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx
      IMPLICIT NONE

      INTEGER AC_ORDER     (NOBTYPMX)  !IN groupings
      INTEGER NO_ITERATIONS(NOBTYPMX)  !IN iterations per step
      INTEGER INTERVAL_ITER(NOBTYPMX)  !IN steps between AC
      INTEGER N_ANAL_LEVS  (NOBTYPMX)  !IN analysis levels
      INTEGER N_WT_LEVS    (NOBTYPMX)  !IN weights levels
      INTEGER AGRES_ROWS   (NOBTYPMX)  !IN ratio modelrows/ anl rows
      INTEGER AGRES_PTS    (NOBTYPMX)  !IN ratio model pts/anl pts
      INTEGER MODE_HANAL   (NOBTYPMX)  !IN mode of horizonatl analysis
      REAL    FI_VAR_FACTOR(NOBTYPMX)  !IN group dep scaling in FI
      REAL    NUDGE_NH(NOBTYPMX)       !IN Nudging coeff NH
      REAL    NUDGE_TR(NOBTYPMX)       !IN Nudging coeff TR
      REAL    NUDGE_SH(NOBTYPMX)       !IN Nudging coeff SH
      REAL    NUDGE_LAM(NOBTYPMX)      !IN Nudging coeff LAM
      INTEGER ICODE                    !OUT error code and message
      CHARACTER(LEN=256) CMESSAGE
!*---------------------------------------------------------------------
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
!*---------------------------------------------------------------------
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


!L---------------------------------------------------------------------
!L Local work space and variables
      INTEGER JOBT,JOBT2,JG,N_OBTYP,LAST_GROUP,FIRST_TYPE,LAST_TYPE
      INTEGER THIS_TYPE,THIS_TYPE_DEF
      INTEGER THIS_GROUP,NO_GROUPS,NCOUNT
      INTEGER NEW_NO_ITERS(NOBTYPMX),NEW_INT_ITER(NOBTYPMX)
      INTEGER NEW_AGRES_ROWS(NOBTYPMX),NEW_AGRES_PTS(NOBTYPMX)
      INTEGER NEW_ANAL_LEVS(NOBTYPMX),NEW_WT_LEVS(NOBTYPMX)
      INTEGER NEW_MODE_HANAL(NOBTYPMX)
      LOGICAL L_NEW_GROUPS,L_NO_ITERS,L_INT_ITER
      LOGICAL L_AGRES_ROWS,L_AGRES_PTS,L_ANAL_LEVS,L_WT_LEVS
      LOGICAL L_MODE_HANAL,L_FI_VAR_FACTOR
      REAL    NEW_FI_VAR_FACTOR(NOBTYPMX)
      REAL    NEW_NUDGE_NH(NOBTYPMX),NEW_NUDGE_TR(NOBTYPMX)
      REAL    NEW_NUDGE_SH(NOBTYPMX)
      LOGICAL L_NUDGE_NH,L_NUDGE_TR,L_NUDGE_SH
      REAL    NEW_NUDGE_LAM(NOBTYPMX)
      LOGICAL L_NUDGE_LAM
      LOGICAL LKNOWN,LFOUND
      LOGICAL L_DUMMY  ! dummy logical

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!L---------------------------------------------------------------------

!     Check and count values in namelist AC_ORDER
      IF (lhook) CALL dr_hook('GROUP_DEP_VAR',zhook_in,zhook_handle)
      N_OBTYP=0
      DO JOBT=1,NOBTYPMX
        IF (AC_ORDER(JOBT) >  0) THEN
          N_OBTYP = N_OBTYP+1
        ELSEIF (AC_ORDER(JOBT) /= IMDI) THEN
          ICODE = 1
          CMESSAGE = 'GROUPDEP : Negative Value in AC_ORDER ?'
          GO TO 999
        ELSE
        ENDIF
      ENDDO

      L_NEW_GROUPS = N_OBTYP >  0
      IF (L_NEW_GROUPS) THEN

!       Check validity of obs types in namelist AC_ORDER.
        DO JOBT=1,N_OBTYP
          IF (AC_ORDER(JOBT) >  0) THEN
            THIS_TYPE = MOD(AC_ORDER(JOBT),1000)
            IF (THIS_TYPE == 501) THEN
              THIS_TYPE = 302
              THIS_GROUP = AC_ORDER(JOBT)/1000
              AC_ORDER(JOBT) = THIS_GROUP*1000 + THIS_TYPE
              if(mype == 0)                                             &
     &        print *, 'Type 501 in AC_ORDER changes to Type 302'
            ENDIF
            IF (THIS_TYPE == 502) THEN
              THIS_TYPE = 305
              THIS_GROUP = AC_ORDER(JOBT)/1000
              AC_ORDER(JOBT) = THIS_GROUP*1000 + THIS_TYPE
              if(mype == 0)                                             &
     &        print *, 'Type 502 in AC_ORDER changes to Type 305'
            ENDIF
            LKNOWN  = .FALSE.
            DO JOBT2 = 1,NOBTYPMX
              THIS_TYPE_DEF = MOD(DEF_AC_ORDER(JOBT2),1000)
              IF (THIS_TYPE == THIS_TYPE_DEF) THEN
                LKNOWN = .TRUE.
              ENDIF
            ENDDO
            IF (.NOT.LKNOWN) THEN
              ICODE = 1
              CMESSAGE =                                                &
     &        'GROUPDEP : Obs Type in AC_ORDER not known.'
              GO TO 999
            ENDIF
          ENDIF
        ENDDO

!       Check that all obs types in AC_OBS_TYPES are in AC_ORDER.
        DO JOBT=1,NOBTYPMX
          IF (AC_OBS_TYPES(JOBT) >  0) THEN
            LFOUND = .FALSE.
            DO JOBT2=1,N_OBTYP
              THIS_TYPE = MOD(AC_ORDER(JOBT2),1000)
              IF (THIS_TYPE == AC_OBS_TYPES(JOBT)) THEN
                LFOUND = .TRUE.
              ENDIF
            ENDDO
            IF (.NOT.LFOUND) THEN
              ICODE=1
              CMESSAGE ='GROUPDEP : AC_OBS_TYPES / AC_ORDER mismatch'
              GO TO 999
            ENDIF
          ENDIF
        ENDDO

!       Check validity and order of group numbers in AC_ORDER.
        LAST_GROUP = 0
        DO JOBT =1,NOBTYPMX
          IF (AC_ORDER(JOBT) >  0) THEN
            THIS_GROUP = AC_ORDER(JOBT)/1000
            IF (THIS_GROUP == 0) THEN
              ICODE = 1
              CMESSAGE =                                                &
     &        'GROUPDEP : Obs Type in AC_ORDER with no group number.'
              GO TO 999
            ENDIF
            IF (THIS_GROUP /= LAST_GROUP .AND.                          &
     &          THIS_GROUP /= LAST_GROUP+1 )THEN
              ICODE = 1
              CMESSAGE =                                                &
     &        'GROUPDEP : Order of groups in AC_ORDER incorrect.'
              GO TO 999
            ENDIF
            LAST_GROUP = THIS_GROUP
          ENDIF
        ENDDO

      ELSE   !  AC_ORDER not used ; Get no of obs types.

        N_OBTYP = 0
        DO JOBT=1,NOBTYPMX
          AC_ORDER(JOBT) = DEF_AC_ORDER(JOBT)
          IF (AC_ORDER(JOBT) >  0) THEN
            N_OBTYP = N_OBTYP+1
          ENDIF
        ENDDO

      ENDIF

      LAST_GROUP = 0
      NO_GROUPS  = 0
      DO JOBT=1, NOBTYPMX
        IF (AC_ORDER(JOBT) >  0) THEN
          THIS_GROUP = AC_ORDER(JOBT)/1000
          IF (THIS_GROUP /= LAST_GROUP) THEN
            NO_GROUPS = NO_GROUPS+1
          ENDIF
          LAST_GROUP = THIS_GROUP
        ENDIF
      ENDDO

      IF (NPROG == 1001.AND.mype == 0) THEN

        print '(10I7)', (AC_ORDER(JOBT),JOBT=1,NOBTYPMX)
        print *, 'N_OBTYP = ',N_OBTYP
        print *, 'NO_ITERATIONS'
        print '(10I7)', (NO_ITERATIONS(JOBT),JOBT=1,NOBTYPMX)
        print *, 'INTERVAL_ITER'
        print '(10I7)', (INTERVAL_ITER(JOBT),JOBT=1,NOBTYPMX)
      ENDIF

!     Process NO_ITERATIONS
!     ---------------------
      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (NO_ITERATIONS(JG) >  0) THEN
          IF (JG >  NO_GROUPS) THEN
            ICODE = 1
            CMESSAGE =                                                  &
     &      'GROUPDEP : NO_ITERATIONS used incorrectly.'
            GO TO 999
          ENDIF
          NCOUNT = NCOUNT+1
        ELSEIF (NO_ITERATIONS(JG) /= IMDI) THEN
          ICODE = 1
          CMESSAGE =                                                    &
     &    'GROUPDEP : Invalid value given in NO_ITERATIONS'
          GO TO 999
        ELSE
        ENDIF
      ENDDO

      L_NO_ITERS = NCOUNT >  0
      IF (L_NEW_GROUPS .AND. L_NO_ITERS .AND. NCOUNT /= NO_GROUPS) THEN
        ICODE = 1
        CMESSAGE =                                                      &
     &  'GROUPDEP : Wrong no of values given in NO_ITERATIONS'
        GO TO 999
      ENDIF

!     Process INTERVAL_ITER
!     ---------------------
      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (INTERVAL_ITER(JG) >  0) THEN
          IF (JG >  NO_GROUPS) THEN
            ICODE = 1
            CMESSAGE =                                                  &
     &      'GROUPDEP : INTERVAL_ITER used incorrectly.'
            GO TO 999
          ENDIF
          NCOUNT = NCOUNT+1
        ELSEIF (INTERVAL_ITER(JG) /= IMDI) THEN
          ICODE = 1
          CMESSAGE = 'GROUPDEP : Invalid value in INTERVAL_ITER'
          GO TO 999
        ELSE
        ENDIF
      ENDDO

      L_INT_ITER = NCOUNT >  0
      IF (L_NEW_GROUPS .AND. L_INT_ITER .AND. NCOUNT /= NO_GROUPS) THEN
        ICODE = 1
        CMESSAGE = 'GROUPDEP : Wrong no of values in INTERVAL_ITER'
        GO TO 999
      ENDIF

!     Process N_ANAL_LEVS
!     --------------------
      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (N_ANAL_LEVS(JG) >  0) THEN
          IF (JG >  NO_GROUPS) THEN
            ICODE = 1
            CMESSAGE =                                                  &
     &      'GROUPDEP : N_ANAL_LEVS used incorrectly.'
            GO TO 999
          ENDIF
          NCOUNT = NCOUNT+1
        ELSEIF (N_ANAL_LEVS(JG) /= IMDI) THEN
          ICODE = 1
          CMESSAGE = 'GROUPDEP : Invalid value in N_ANAL_LEVS'
          GO TO 999
        ELSE
        ENDIF
      ENDDO

      L_ANAL_LEVS = NCOUNT >  0
      IF (L_NEW_GROUPS .AND. L_ANAL_LEVS .AND. NCOUNT /= NO_GROUPS)THEN
        ICODE = 1
        CMESSAGE = 'GROUPDEP : Wrong no of values in N_ANAL_LEVS'
        GO TO 999
      ENDIF

!     Process N_WT_LEVS
!     --------------------
      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (N_WT_LEVS(JG) >  0) THEN
          IF (JG >  NO_GROUPS) THEN
            ICODE = 1
            CMESSAGE =                                                  &
     &      'GROUPDEP : N_WT_LEVS used incorrectly.'
            GO TO 999
          ENDIF
          NCOUNT = NCOUNT+1
        ELSEIF (N_WT_LEVS(JG) /= IMDI) THEN
          ICODE = 1
          CMESSAGE = 'GROUPDEP : Invalid value in N_WT_LEVS'
          GO TO 999
        ELSE
        ENDIF
      ENDDO

      L_WT_LEVS = NCOUNT >  0
      IF (L_NEW_GROUPS .AND. L_WT_LEVS .AND. NCOUNT /= NO_GROUPS)THEN
        ICODE = 1
        CMESSAGE = 'GROUPDEP : Wrong no of values in N_WT_LEVS'
        GO TO 999
      ENDIF

      IF (model_domain == mt_global) THEN
!     Process NUDGE_NH
!     ----------------
      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (NUDGE_NH(JG) >  0.0) THEN
          IF (JG >  NO_GROUPS) THEN
            ICODE = 1
            CMESSAGE =                                                  &
     &      'GROUPDEP : NUDGE_NH used incorrectly.'
            GO TO 999
          ENDIF
          NCOUNT = NCOUNT+1
        ELSEIF (NUDGE_NH(JG) /= RMDI) THEN
          ICODE = 1
          CMESSAGE = 'GROUPDEP : Invalid value in NUDGE_NH'
          GO TO 999
        ELSE
        ENDIF
      ENDDO

      L_NUDGE_NH = NCOUNT >  0

      IF (L_NEW_GROUPS .AND. L_NUDGE_NH .AND. NCOUNT /= NO_GROUPS) THEN
        ICODE = 1
        CMESSAGE = 'GROUPDEP : Wrong no of values in NUDGE_NH'
        GO TO 999
      ENDIF

!     Process NUDGE_TR
!     ----------------
      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (NUDGE_TR(JG) >  0.0) THEN
          IF (JG >  NO_GROUPS) THEN
            ICODE = 1
            CMESSAGE = 'GROUPDEP : NUDGE_TR used incorrectly.'
            GO TO 999
          ENDIF
          NCOUNT = NCOUNT+1
        ELSEIF (NUDGE_TR(JG) /= RMDI) THEN
          ICODE = 1
          CMESSAGE = 'GROUPDEP : Invalid value in NUDGE_TR'
          GO TO 999
        ELSE
        ENDIF
      ENDDO

      L_NUDGE_TR = NCOUNT >  0

      IF (L_NEW_GROUPS .AND. L_NUDGE_TR .AND. NCOUNT /= NO_GROUPS) THEN
        ICODE = 1
        CMESSAGE = 'GROUPDEP : Wrong no of values in NUDGE_TR'
        GO TO 999
      ENDIF

!     Process NUDGE_SH
!     ----------------
      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (NUDGE_SH(JG) >  0.0) THEN
          IF (JG >  NO_GROUPS) THEN
            ICODE = 1
            CMESSAGE = 'GROUPDEP : NUDGE_SH used incorrectly.'
            GO TO 999
          ENDIF
          NCOUNT = NCOUNT+1
        ELSEIF (NUDGE_SH(JG) /= RMDI) THEN
          ICODE = 1
          CMESSAGE = 'GROUPDEP : Invalid value in NUDGE_SH'
          GO TO 999
        ELSE
        ENDIF
      ENDDO

      L_NUDGE_SH = NCOUNT >  0

      IF (L_NEW_GROUPS .AND. L_NUDGE_SH .AND. NCOUNT /= NO_GROUPS) THEN
        ICODE = 1
        CMESSAGE = 'GROUPDEP : Wrong no of values in NUDGE_SH'
        GO TO 999
      ENDIF
      ELSE
!     Process NUDGE_LAM
!     -----------------
      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (NUDGE_LAM(JG) >  0.0) THEN
          IF (JG >  NO_GROUPS) THEN
            ICODE = 1
            CMESSAGE = 'GROUPDEP : NUDGE_LAM used incorrectly.'
            GO TO 999
          ENDIF
          NCOUNT = NCOUNT+1
        ELSEIF (NUDGE_LAM(JG) /= RMDI) THEN
          ICODE = 1
          CMESSAGE = 'GROUPDEP : Invalid value in NUDGE_LAM'
          GO TO 999
        ELSE
        ENDIF
      ENDDO

      L_NUDGE_LAM = NCOUNT >  0

      IF (L_NEW_GROUPS .AND. L_NUDGE_LAM .AND. NCOUNT /= NO_GROUPS)THEN
        ICODE = 1
        CMESSAGE = 'GROUPDEP : Wrong no of values in NUDGE_LAM'
        GO TO 999
      ENDIF
      END IF  ! if GLOBAL

!     Process AGRES_ROWS
!     ------------------
      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (AGRES_ROWS(JG) >  0) THEN
          IF (JG >  NO_GROUPS) THEN
            ICODE = 1
            CMESSAGE = 'GROUPDEP : AGRES_ROWS used incorrectly.'
            GO TO 999
          ENDIF
          NCOUNT = NCOUNT+1
        ELSEIF (AGRES_ROWS(JG) /= IMDI) THEN
          ICODE = 1
          CMESSAGE = 'GROUPDEP : Invalid value in AGRES_ROWS'
          GO TO 999
        ELSE
        ENDIF
      ENDDO

      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (AGRES_ROWS(JG) >  1) THEN
         if(mype == 0)then
          print *, ' AGRES_ROWS /= 1 disallowed for MPP'
         endif
        ENDIF
      ENDDO
      L_AGRES_ROWS = NCOUNT >  0
      IF (L_NEW_GROUPS .AND. L_AGRES_ROWS .AND. NCOUNT /= NO_GROUPS)THEN
        ICODE = 1
        CMESSAGE = 'GROUPDEP : Wrong no of values in AGRES_ROWS'
        GO TO 999
      ENDIF

!     Process AGRES_PTS
!     -----------------
      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (AGRES_PTS(JG) >  0) THEN
          IF (JG >  NO_GROUPS) THEN
            ICODE = 1
            CMESSAGE = 'GROUPDEP : AGRES_PTS used incorrectly.'
            GO TO 999
          ENDIF
          NCOUNT = NCOUNT+1
        ELSEIF (AGRES_PTS(JG) /= IMDI) THEN
          ICODE = 1
          CMESSAGE = 'GROUPDEP : Invalid value in AGRES_PTS'
          GO TO 999
        ELSE
        ENDIF
      ENDDO

      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (AGRES_PTS(JG) >  1) THEN
         if(mype == 0)then
          print *, ' AGRES_PTS /= 1 disallowed for MPP'
         endif
        ENDIF
      ENDDO

      L_AGRES_PTS = NCOUNT >  0

      IF (L_NEW_GROUPS .AND. L_AGRES_PTS .AND. NCOUNT /= NO_GROUPS)THEN
        ICODE = 1
        CMESSAGE = 'GROUPDEP : Wrong no of values in AGRES_PTS'
        GO TO 999
      ENDIF

!     Process MODE_HANAL
!     ------------------
      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (MODE_HANAL(JG) == 1 .OR. MODE_HANAL(JG) == 2) THEN
          IF (JG >  NO_GROUPS) THEN
            ICODE = 1
            CMESSAGE = 'GROUPDEP : MODE_HANAL used incorrectly.'
            GO TO 999
          ENDIF
          NCOUNT = NCOUNT+1
        ELSEIF (MODE_HANAL(JG) /= IMDI) THEN
          ICODE = 1
          CMESSAGE = 'GROUPDEP : Invalid value in MODE_HANAL'
          GO TO 999
        ELSE
        ENDIF
      ENDDO

      L_MODE_HANAL = NCOUNT >  0

      IF(L_NEW_GROUPS .AND. L_MODE_HANAL .AND. NCOUNT /= NO_GROUPS)THEN
        ICODE = 1
        CMESSAGE = 'GROUPDEP : Wrong no of values in MODE_HANAL'
        GO TO 999
      ENDIF

!     Process FI_VAR_FACTOR
!     ---------------------
      NCOUNT = 0
      DO JG=1,NOBTYPMX
        IF (FI_VAR_FACTOR(JG) >= 0.0)THEN
          IF (JG >  NO_GROUPS) THEN
            ICODE = 1
            CMESSAGE = 'GROUPDEP : FI_VAR_FACTOR used incorrectly.'
            GO TO 999
          ENDIF
          NCOUNT = NCOUNT+1
        ELSEIF (FI_VAR_FACTOR(JG) /= RMDI) THEN
          ICODE = 1
          CMESSAGE = 'GROUPDEP : Invalid value in FI_VAR_FACTOR'
          GO TO 999
        ELSE
        ENDIF
      ENDDO

      L_FI_VAR_FACTOR = NCOUNT >  0

      IF(L_NEW_GROUPS .AND. L_FI_VAR_FACTOR                             &
     &    .AND. NCOUNT /= NO_GROUPS)THEN
        ICODE = 1
        CMESSAGE = 'GROUPDEP : Wrong no of values in FI_VAR_FACTOR'
        GO TO 999
      ENDIF

      if(mype == 0)then
      print *
      IF (L_NEW_GROUPS) print *, 'AC_ORDER      used in namelist'
      IF (L_NO_ITERS)   print *, 'NO_ITERATIONS used in namelist'
      IF (L_INT_ITER)   print *, 'INTERVAL_ITER used in namelist'
      IF (L_ANAL_LEVS)  print *, 'N_ANAL_LEVS   used in namelist'
      IF (L_WT_LEVS)    print *, 'N_WT_LEVS     used in namelist'
      IF (L_AGRES_ROWS) print *, 'AGRES_ROWS    used in namelist'
      IF (L_AGRES_PTS)  print *, 'AGRES_PTS     used in namelist'
      IF (L_FI_VAR_FACTOR)print *,'FI_VAR_FACTOR used in namelist'
      IF (model_domain == mt_global) THEN
      IF (L_NUDGE_NH)   print *, 'NUDGE_NH      used in namelist'
      IF (L_NUDGE_TR)   print *, 'NUDGE_TR      used in namelist'
      IF (L_NUDGE_SH)   print *, 'NUDGE_SH      used in namelist'
      ELSE
      IF (L_NUDGE_LAM)  print *, 'NUDGE_LAM     used in namelist'
      END IF
      endif

!     --------------------------------------------------------------
!     If AC_ORDER has been used to change the order or groups
!     but none of the group dependent arrays used, then this routine
!     attempts to derive new defaults for the new order or groups
!     from the existing defaults.
!     --------------------------------------------------------------
      IF (L_NEW_GROUPS) THEN  !  New order or groups.
        IF (model_domain == mt_global) THEN
           L_DUMMY = .NOT.L_NUDGE_NH  .OR. .NOT.L_NUDGE_TR  .OR.        &
                     .NOT.L_NUDGE_SH
        ELSE
           L_DUMMY = .NOT.L_NUDGE_LAM
        END IF

        IF (.NOT.L_NO_ITERS  .OR. .NOT.L_INT_ITER  .OR.                 &
     &      .NOT.L_ANAL_LEVS .OR. .NOT.L_WT_LEVS   .OR.                 &
            L_DUMMY .OR.                                                &
     &      .NOT.L_AGRES_ROWS .OR. .NOT.L_AGRES_PTS .OR.                &
     &      .NOT.L_MODE_HANAL .OR. .NOT.L_FI_VAR_FACTOR)THEN

!         Set up new defaults in NEW_ arrays from existing
!         defaults. New defaults are set up for all observation
!         types in AC_ORDER from the corresponding default array.
!         eg. NEW_NO_ITERS is set up from DEF_NO_ITERATIONS.
!         New defaults are set up for all group dependent arrays
!         here. Those not required are ignored later in this routine.

          IF (model_domain == mt_global) THEN
          DO JOBT=1,N_OBTYP
            THIS_TYPE = MOD(AC_ORDER(JOBT),1000)
            DO JOBT2=1,NOBTYPMX
              THIS_TYPE_DEF = MOD(DEF_AC_ORDER(JOBT2),1000)
              IF (THIS_TYPE == THIS_TYPE_DEF) THEN
                THIS_GROUP = DEF_AC_ORDER(JOBT2)/1000
                NEW_NO_ITERS(JOBT)   = DEF_NO_ITERATIONS(THIS_GROUP)
                NEW_INT_ITER(JOBT)   = DEF_INTERVAL_ITER(THIS_GROUP)
                NEW_ANAL_LEVS(JOBT)  = DEF_NO_ANAL_LEVS(THIS_GROUP)
                NEW_WT_LEVS(JOBT)    = DEF_NO_WT_LEVS(THIS_GROUP)
                NEW_NUDGE_NH(JOBT)   = DEF_NUDGE_NH(THIS_GROUP)
                NEW_NUDGE_TR(JOBT)   = DEF_NUDGE_TR(THIS_GROUP)
                NEW_NUDGE_SH(JOBT)   = DEF_NUDGE_SH(THIS_GROUP)
                NEW_AGRES_ROWS(JOBT) = DEF_AGRES_ROWS(THIS_GROUP)
                NEW_AGRES_PTS(JOBT)  = DEF_AGRES_PTS(THIS_GROUP)
                NEW_MODE_HANAL(JOBT) = DEF_MODE_HANAL(THIS_GROUP)
                NEW_FI_VAR_FACTOR(JOBT) = DEF_FI_VAR_FACTOR(THIS_GROUP)
              END IF
            END DO
          END DO
          ELSE
          DO JOBT=1,N_OBTYP
            THIS_TYPE = MOD(AC_ORDER(JOBT),1000)
            DO JOBT2=1,NOBTYPMX
              THIS_TYPE_DEF = MOD(DEF_AC_ORDER(JOBT2),1000)
              IF (THIS_TYPE == THIS_TYPE_DEF) THEN
                THIS_GROUP = DEF_AC_ORDER(JOBT2)/1000
                NEW_NO_ITERS(JOBT)   = DEF_NO_ITERATIONS(THIS_GROUP)
                NEW_INT_ITER(JOBT)   = DEF_INTERVAL_ITER(THIS_GROUP)
                NEW_ANAL_LEVS(JOBT)  = DEF_NO_ANAL_LEVS(THIS_GROUP)
                NEW_WT_LEVS(JOBT)    = DEF_NO_WT_LEVS(THIS_GROUP)
                NEW_NUDGE_LAM(JOBT)  = DEF_NUDGE_LAM(THIS_GROUP)
                NEW_AGRES_ROWS(JOBT) = DEF_AGRES_ROWS(THIS_GROUP)
                NEW_AGRES_PTS(JOBT)  = DEF_AGRES_PTS(THIS_GROUP)
                NEW_MODE_HANAL(JOBT) = DEF_MODE_HANAL(THIS_GROUP)
                NEW_FI_VAR_FACTOR(JOBT) = DEF_FI_VAR_FACTOR(THIS_GROUP)
              END IF
            END DO
          END DO
          END IF  ! if GLOBAL

        IF (NPROG == 1001.AND.mype == 0) THEN
          print *, 'NEW_NO_ITERS'
          print '(10I7)', (NEW_NO_ITERS(JOBT),JOBT=1,N_OBTYP)
          print *, 'NEW_INT_ITER'
          print '(10I7)', (NEW_INT_ITER(JOBT),JOBT=1,N_OBTYP)
          print *, 'NEW_ANAL_LEVS'
          print '(10I7)', (NEW_ANAL_LEVS(JOBT),JOBT=1,N_OBTYP)
          print *, 'NEW_WT_LEVS'
          print '(10I7)', (NEW_WT_LEVS(JOBT),JOBT=1,N_OBTYP)
          IF (model_domain == mt_global) THEN
          print *, 'NEW_NUDGE_NH'
          print '(10E10.3)', (NEW_NUDGE_NH(JOBT),JOBT=1,N_OBTYP)
          print *, 'NEW_NUDGE_TR'
          print '(10E10.3)', (NEW_NUDGE_TR(JOBT),JOBT=1,N_OBTYP)
          print *, 'NEW_NUDGE_SH'
          print '(10E10.3)', (NEW_NUDGE_SH(JOBT),JOBT=1,N_OBTYP)
          ELSE
          print *, 'NEW_NUDGE_LAM'
          print '(10E10.3)', (NEW_NUDGE_LAM(JOBT),JOBT=1,N_OBTYP)
          END IF
          print *, 'NEW_AGRES_ROWS'
          print '(10I7)', (NEW_AGRES_ROWS(JOBT),JOBT=1,N_OBTYP)
          print *, 'NEW_AGRES_PTS'
          print '(10I7)', (NEW_AGRES_PTS(JOBT),JOBT=1,N_OBTYP)
          print *, 'NEW_MODE_HANAL'
          print '(10I7)', (NEW_MODE_HANAL(JOBT),JOBT=1,N_OBTYP)
          print *, 'NEW_FI_VAR_FACTOR'
          print '(10E10.3)',(NEW_FI_VAR_FACTOR(JOBT),JOBT=1,N_OBTYP)
        ENDIF

!         Go through groups and check that the new default for each
!         obs type in the group matches. If there is a mismatch, then
!         the routine will abort and prompt the user to use the
!         appropriate namelist array to specify new values for the
!         new groups. If the new defaults for all the obs types in the
!         group match, then this becomes the default for the new group.

          DO JG=1,NO_GROUPS

            FIRST_TYPE = 0
            LAST_TYPE  = 0
            DO JOBT=1,N_OBTYP
              THIS_GROUP = AC_ORDER(JOBT)/1000
              IF (THIS_GROUP == JG) THEN
                IF (FIRST_TYPE == 0) THEN
                  FIRST_TYPE = JOBT
                ENDIF
                LAST_TYPE = JOBT
              ENDIF
            ENDDO

            IF (FIRST_TYPE == 0 .OR. LAST_TYPE == 0) THEN
              ICODE = 1
              CMESSAGE ='GROUPDEP : FIRST_TYPE=0 or LAST_TYPE=0 ?'
              GO TO 999
            ENDIF

            IF (.NOT.L_NO_ITERS) THEN
              IF (FIRST_TYPE <  LAST_TYPE) THEN
                DO JOBT=FIRST_TYPE+1,LAST_TYPE
                IF (NEW_NO_ITERS(FIRST_TYPE) /= NEW_NO_ITERS(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =                                            &
     &            'GROUPDEP : NO_ITERATIONS must be used in namelist'
                  GO TO 999
                ENDIF
                ENDDO
              ENDIF
              DEF_NO_ITERATIONS(JG) = NEW_NO_ITERS(FIRST_TYPE)
            ENDIF

            IF (.NOT.L_INT_ITER) THEN
              IF (FIRST_TYPE <  LAST_TYPE) THEN
                DO JOBT=FIRST_TYPE+1,LAST_TYPE
                IF (NEW_INT_ITER(FIRST_TYPE) /= NEW_INT_ITER(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =                                            &
     &            'GROUPDEP : INTERVAL_ITER must be used in namelist'
                  GO TO 999
                ENDIF
                ENDDO
              ENDIF
              DEF_INTERVAL_ITER(JG) = NEW_INT_ITER(FIRST_TYPE)
            ENDIF

            IF (.NOT.L_ANAL_LEVS) THEN
              IF (FIRST_TYPE <  LAST_TYPE) THEN
                DO JOBT=FIRST_TYPE+1,LAST_TYPE
                IF(NEW_ANAL_LEVS(FIRST_TYPE) /= NEW_ANAL_LEVS(JOBT))THEN
                  ICODE = 1
                  CMESSAGE =                                            &
     &            'GROUPDEP : N_ANAL_LEVS must be used in namelist'
                  GO TO 999
                ENDIF
                ENDDO
              ENDIF
              DEF_NO_ANAL_LEVS(JG) = NEW_ANAL_LEVS(FIRST_TYPE)
            ENDIF

            IF (.NOT.L_WT_LEVS) THEN
              IF (FIRST_TYPE <  LAST_TYPE) THEN
                DO JOBT=FIRST_TYPE+1,LAST_TYPE
                IF (NEW_WT_LEVS(FIRST_TYPE) /= NEW_WT_LEVS(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =                                            &
     &            'GROUPDEP : N_WT_LEVS must be used in namelist'
                  GO TO 999
                ENDIF
                ENDDO
              ENDIF
              DEF_NO_WT_LEVS(JG) = NEW_WT_LEVS(FIRST_TYPE)
            ENDIF

            IF (.NOT.L_AGRES_ROWS) THEN
              IF (FIRST_TYPE <  LAST_TYPE) THEN
                DO JOBT=FIRST_TYPE+1,LAST_TYPE
                IF (NEW_AGRES_ROWS(FIRST_TYPE)  /=                      &
     &              NEW_AGRES_ROWS(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =                                            &
     &            'GROUPDEP : AGRES_ROWS must be used in namelist'
                  GO TO 999
                ENDIF
                ENDDO
              ENDIF
              DEF_AGRES_ROWS(JG) = NEW_AGRES_ROWS(FIRST_TYPE)
            ENDIF

            IF (.NOT.L_AGRES_PTS) THEN
              IF (FIRST_TYPE <  LAST_TYPE) THEN
                DO JOBT=FIRST_TYPE+1,LAST_TYPE
                IF (NEW_AGRES_PTS(FIRST_TYPE)  /=                       &
     &              NEW_AGRES_PTS(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =                                            &
     &            'GROUPDEP : AGRES_PTS must be used in namelist'
                  GO TO 999
                ENDIF
                ENDDO
              ENDIF
              DEF_AGRES_PTS(JG) = NEW_AGRES_PTS(FIRST_TYPE)
            ENDIF

            IF (.NOT.L_MODE_HANAL) THEN
              IF (FIRST_TYPE <  LAST_TYPE) THEN
                DO JOBT=FIRST_TYPE+1,LAST_TYPE
                IF (NEW_MODE_HANAL(FIRST_TYPE)  /=                      &
     &              NEW_MODE_HANAL(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =                                            &
     &            'GROUPDEP : MODE_HANAL must be used in namelist'
                  GO TO 999
                ENDIF
                ENDDO
              ENDIF
              DEF_MODE_HANAL(JG) = NEW_MODE_HANAL(FIRST_TYPE)
            ENDIF

            IF (.NOT.L_FI_VAR_FACTOR) THEN
              IF (FIRST_TYPE <  LAST_TYPE) THEN
                DO JOBT=FIRST_TYPE+1,LAST_TYPE
                IF (NEW_FI_VAR_FACTOR(FIRST_TYPE)  /=                   &
     &              NEW_FI_VAR_FACTOR(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =                                            &
     &            'GROUPDEP : FI_VAR_FACTOR must be used in namelist'
                  GO TO 999
                ENDIF
                ENDDO
              ENDIF
              DEF_FI_VAR_FACTOR(JG) = NEW_FI_VAR_FACTOR(FIRST_TYPE)
            ENDIF

            IF (model_domain == mt_global) THEN
            IF (.NOT.L_NUDGE_NH) THEN
              IF (FIRST_TYPE <  LAST_TYPE) THEN
                DO JOBT=FIRST_TYPE+1,LAST_TYPE
                IF (NEW_NUDGE_NH(FIRST_TYPE) /= NEW_NUDGE_NH(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =                                            &
     &            'GROUPDEP : NUDGE_NH must be used in namelist'
                  GO TO 999
                ENDIF
                ENDDO
              ENDIF
              DEF_NUDGE_NH(JG) = NEW_NUDGE_NH(FIRST_TYPE)
            ENDIF

            IF (.NOT.L_NUDGE_TR) THEN
              IF (FIRST_TYPE <  LAST_TYPE) THEN
                DO JOBT=FIRST_TYPE+1,LAST_TYPE
                IF (NEW_NUDGE_TR(FIRST_TYPE) /= NEW_NUDGE_TR(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =                                            &
     &            'GROUPDEP : NUDGE_TR must be used in namelist'
                  GO TO 999
                ENDIF
                ENDDO
              ENDIF
              DEF_NUDGE_TR(JG) = NEW_NUDGE_TR(FIRST_TYPE)
            ENDIF

            IF (.NOT.L_NUDGE_SH) THEN
              IF (FIRST_TYPE <  LAST_TYPE) THEN
                DO JOBT=FIRST_TYPE+1,LAST_TYPE
                IF (NEW_NUDGE_SH(FIRST_TYPE) /= NEW_NUDGE_SH(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =                                            &
     &            'GROUPDEP : NUDGE_SH must be used in namelist'
                  GO TO 999
                ENDIF
                ENDDO
              ENDIF
              DEF_NUDGE_SH(JG) = NEW_NUDGE_SH(FIRST_TYPE)
            ENDIF

            ELSE
            IF (.NOT.L_NUDGE_LAM) THEN
              IF (FIRST_TYPE <  LAST_TYPE) THEN
                DO JOBT=FIRST_TYPE+1,LAST_TYPE
                IF (NEW_NUDGE_LAM(FIRST_TYPE)  /=                       &
     &              NEW_NUDGE_LAM(JOBT)) THEN
                  ICODE = 1
                  CMESSAGE =                                            &
     &            'GROUPDEP : NUDGE_LAM must be used in namelist'
                  GO TO 999
                ENDIF
                ENDDO
              ENDIF
              DEF_NUDGE_LAM(JG) = NEW_NUDGE_LAM(FIRST_TYPE)
            ENDIF

            END IF  ! if GLOBAL
          ENDDO

        ENDIF
      ENDIF

!     Overwrite existing defaults with new defaults.

      IF (L_NEW_GROUPS) THEN   !  New values for DEF_AC_ORDER
        DO JOBT=1,NOBTYPMX
          DEF_AC_ORDER(JOBT) = AC_ORDER(JOBT)
        ENDDO
      ENDIF

      IF (L_NO_ITERS) THEN   !  New values for DEF_NO_ITERATIONS
        DO JG=1,NO_GROUPS
          IF (NO_ITERATIONS(JG) >  0) THEN
            DEF_NO_ITERATIONS(JG) = NO_ITERATIONS(JG)
          ENDIF
        ENDDO
        IF (NO_GROUPS <  NOBTYPMX) THEN
          DO JG=NO_GROUPS+1,NOBTYPMX
            DEF_NO_ITERATIONS(JG) = IMDI
          ENDDO
        ENDIF
      ENDIF

      IF (L_INT_ITER) THEN   !   New values for DEF_INTERVAL_ITER
        DO JG=1,NO_GROUPS
          IF (INTERVAL_ITER(JG) >  0) THEN
            DEF_INTERVAL_ITER(JG) = INTERVAL_ITER(JG)
          ENDIF
        ENDDO
        IF (NO_GROUPS <  NOBTYPMX) THEN
          DO JG=NO_GROUPS+1,NOBTYPMX
            DEF_INTERVAL_ITER(JG) = IMDI
          ENDDO
        ENDIF
      ENDIF

      IF (L_ANAL_LEVS) THEN   !   New values for DEF_NO_ANAL_LEVS
        DO JG=1,NO_GROUPS
          IF (N_ANAL_LEVS(JG) >  0) THEN
            DEF_NO_ANAL_LEVS(JG) = N_ANAL_LEVS(JG)
          ENDIF
        ENDDO
        IF (NO_GROUPS <  NOBTYPMX) THEN
          DO JG=NO_GROUPS+1,NOBTYPMX
            DEF_NO_ANAL_LEVS(JG) = IMDI
          ENDDO
        ENDIF
      ENDIF

      IF (L_WT_LEVS) THEN   !   New values for DEF_NO_WT_LEVS
        DO JG=1,NO_GROUPS
          IF (N_WT_LEVS(JG) >  0) THEN
            DEF_NO_WT_LEVS(JG) = N_WT_LEVS(JG)
          ENDIF
        ENDDO
        IF (NO_GROUPS <  NOBTYPMX) THEN
          DO JG=NO_GROUPS+1,NOBTYPMX
            DEF_NO_WT_LEVS(JG) = IMDI
          ENDDO
        ENDIF
      ENDIF

      IF (L_AGRES_ROWS) THEN  !   New values for DEF_AGRES_ROWS
        DO JG=1,NO_GROUPS
          IF (AGRES_ROWS(JG) >  0) THEN
            DEF_AGRES_ROWS(JG) = AGRES_ROWS(JG)
          ENDIF
        ENDDO
        IF (NO_GROUPS <  NOBTYPMX) THEN
          DO JG=NO_GROUPS+1,NOBTYPMX
            DEF_AGRES_ROWS(JG) = IMDI
          ENDDO
        ENDIF
      ENDIF

      IF (L_AGRES_PTS) THEN  !   New values for DEF_AGRES_PTS
        DO JG=1,NO_GROUPS
          IF (AGRES_PTS(JG) >  0) THEN
            DEF_AGRES_PTS(JG) = AGRES_PTS(JG)
          ENDIF
        ENDDO
        IF (NO_GROUPS <  NOBTYPMX) THEN
          DO JG=NO_GROUPS+1,NOBTYPMX
            DEF_AGRES_PTS(JG) = IMDI
          ENDDO
        ENDIF
      ENDIF

      IF (L_MODE_HANAL) THEN  !   New values for DEF_MODE_HANAL
        DO JG=1,NO_GROUPS
          IF (MODE_HANAL(JG) >  0) THEN
            DEF_MODE_HANAL(JG) = MODE_HANAL(JG)
          ENDIF
        ENDDO
        IF (NO_GROUPS <  NOBTYPMX) THEN
          DO JG=NO_GROUPS+1,NOBTYPMX
            DEF_MODE_HANAL(JG) = IMDI
          ENDDO
        ENDIF
      ENDIF

      IF (L_FI_VAR_FACTOR) THEN ! New values for DEF_FI_VAR_FACTOR
        DO JG=1,NO_GROUPS
          IF (FI_VAR_FACTOR(JG) >  0.0) THEN
            DEF_FI_VAR_FACTOR(JG) = FI_VAR_FACTOR(JG)
          ENDIF
        ENDDO
        IF (NO_GROUPS <  NOBTYPMX) THEN
          DO JG=NO_GROUPS+1,NOBTYPMX
            DEF_FI_VAR_FACTOR(JG) = RMDI
          ENDDO
        ENDIF
      ENDIF

      IF (model_domain == mt_global) THEN
      IF (L_NUDGE_NH) THEN   !   New values for DEF_NUDGE_NH
        DO JG=1,NO_GROUPS
          IF (NUDGE_NH(JG) >  0.0) THEN
            DEF_NUDGE_NH(JG) = NUDGE_NH(JG)
          ENDIF
        ENDDO
        IF (NO_GROUPS <  NOBTYPMX) THEN
          DO JG=NO_GROUPS+1,NOBTYPMX
            DEF_NUDGE_NH(JG) = RMDI
          ENDDO
        ENDIF
      ENDIF

      IF (L_NUDGE_TR) THEN   !   New values for DEF_NUDGE_TR
        DO JG=1,NO_GROUPS
          IF (NUDGE_TR(JG) >  0.0) THEN
            DEF_NUDGE_TR(JG) = NUDGE_TR(JG)
          ENDIF
        ENDDO
        IF (NO_GROUPS <  NOBTYPMX) THEN
          DO JG=NO_GROUPS+1,NOBTYPMX
            DEF_NUDGE_TR(JG) = RMDI
          ENDDO
        ENDIF
      ENDIF

      IF (L_NUDGE_SH) THEN   !   New values for DEF_NUDGE_SH
        DO JG=1,NO_GROUPS
          IF (NUDGE_SH(JG) >  0.0) THEN
            DEF_NUDGE_SH(JG) = NUDGE_SH(JG)
          ENDIF
        ENDDO
        IF (NO_GROUPS <  NOBTYPMX) THEN
          DO JG=NO_GROUPS+1,NOBTYPMX
            DEF_NUDGE_SH(JG) = RMDI
          ENDDO
        ENDIF
      ENDIF
      ELSE

      IF (L_NUDGE_LAM) THEN   !   New values for DEF_NUDGE_LAM
        DO JG=1,NO_GROUPS
          IF (NUDGE_LAM(JG) >  0.0) THEN
            DEF_NUDGE_LAM(JG) = NUDGE_LAM(JG)
          ENDIF
        ENDDO
        IF (NO_GROUPS <  NOBTYPMX) THEN
          DO JG=NO_GROUPS+1,NOBTYPMX
            DEF_NUDGE_LAM(JG) = RMDI
          ENDDO
        ENDIF
      ENDIF
      END IF  ! if GLOBAL

!     Check that no of wt levs is < or = of no of anal levs.
      IF (L_ANAL_LEVS .OR. L_WT_LEVS) THEN
        DO JG = 1,NO_GROUPS
          IF (DEF_NO_ANAL_LEVS(JG) <  DEF_NO_WT_LEVS(JG)) THEN
            ICODE = 1
            CMESSAGE = 'GROUPDEP: No of Anal Levs < No of Wt Levs ?'
            if(mype == 0)then
            print *, 'Group No ',JG,' No of anal/wt levs = ',           &
     &      DEF_NO_ANAL_LEVS(JG),DEF_NO_WT_LEVS(JG)
            endif
            GO TO 999
          ENDIF
        ENDDO
      ENDIF

      IF (NPROG == 1001.AND.mype == 0) THEN

      print *, 'DEF_AC_ORDER'
      print '(10I7)', (DEF_AC_ORDER(JOBT),JOBT=1,NOBTYPMX)

      print *, 'GROUPNO'
      print '(10I7)', (JG,JG=1,NO_GROUPS)
      print *, 'NO_ITERATIONS'
      print '(10I7)', (NO_ITERATIONS(JG),JG=1,NO_GROUPS)
      print *, 'INTERVAL_ITER'
      print '(10I7)', (INTERVAL_ITER(JG),JG=1,NO_GROUPS)
      print *, 'DEF_NO_ITERATIONS'
      print '(10I7)', (DEF_NO_ITERATIONS(JG),JG=1,NO_GROUPS)
      print *, 'DEF_INTERVAL_ITER'
      print '(10I7)', (DEF_INTERVAL_ITER(JG),JG=1,NO_GROUPS)
      print *, 'DEF_AGRES_ROWS'
      print '(10I7)', (DEF_AGRES_ROWS(JG),JG=1,NO_GROUPS)
      print *, 'DEF_AGRES_PTS'
      print '(10I7)', (DEF_AGRES_PTS(JG),JG=1,NO_GROUPS)
      print *, 'DEF_MODE_HANAL'
      print '(10I7)', (DEF_MODE_HANAL(JG),JG=1,NO_GROUPS)
      print *, 'DEF_FI_VAR_FACTOR'
      print '(10E10.3)', (DEF_FI_VAR_FACTOR(JG),JG=1,NO_GROUPS)
      IF (model_domain == mt_global) THEN
      print *, 'DEF_NUDGE_NH'
      print '(10E10.3)', (DEF_NUDGE_NH(JG),JG=1,NO_GROUPS)
      print *, 'DEF_NUDGE_TR'
      print '(10E10.3)', (DEF_NUDGE_TR(JG),JG=1,NO_GROUPS)
      print *, 'DEF_NUDGE_SH'
      print '(10E10.3)', (DEF_NUDGE_SH(JG),JG=1,NO_GROUPS)
      ELSE
      print *, 'DEF_NUDGE_LAM'
      print '(10E10.3)', (DEF_NUDGE_LAM(JG),JG=1,NO_GROUPS)
      END IF
      ENDIF

 999  CONTINUE
      IF (lhook) CALL dr_hook('GROUP_DEP_VAR',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GROUP_DEP_VAR
END MODULE group_dep_var_mod
