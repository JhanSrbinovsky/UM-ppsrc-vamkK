! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE SETTPS--------------------------------------------------
!LL
!LL  Purpose : Sets up the list of AC Observation Types in the
!LL            order they are to be processed in the assimilation.
!LL            This routine is called each time the AC Observation
!LL            files are read in.
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!*L  Arguments----------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE settps_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE SETTPS (ICODE,CMESSAGE)

      USE um_input_control_mod, ONLY : model_domain
      USE domain_params, ONLY: mt_global

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx
      IMPLICIT NONE
      INTEGER       ICODE
      CHARACTER(LEN=256) CMESSAGE

!    INTENT=OUT--------------------------------------------------------
!     ICODE        : Return Code
!     CMESSAGE     : Reason for failure
!*   ------------------------------------------------------------------

!     The variables/arrays (all in comdeck COMACP) set up are :-

!     NACT  : Number of AC Observation Types to be processed.
!     LACT  : List of the Obs Types in the order to be processed.
!     N_GROUPS : Number of groups for processing in AC.
!     GROUP_NO : Group in which each type is to be processed.

!     The order of processing is controlled by the array AC_ORDER.

!     The array AC_OBS_TYPES in the ACP namelist is used to control
!     which observation types are to be processed.

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


!     Local variables

      INTEGER THIS_TYPE,J,JOBT,JTYPE
      INTEGER LAST_GROUP,THIS_GROUP
      LOGICAL LUSE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     THIS_TYPE  : Obs type in AC_ORDER
!     J          : Loop counter for obs types in LACT
!     JOBT       : Loop counter for obs types in AC_ORDER
!     JTYPE      : Loop counter for obs types in AC_ORDER
!     THIS_GROUP : Indicator from AC_ORDER of grouping (current type)
!     LAST_GROUP : Indicator from AC_ORDER of grouping (previous type)

!     ------------------------------------------------------------------
!L    1. Initialise arrays and variables set up by SETTPS
!     ------------------------------------------------------------------
      IF (lhook) CALL dr_hook('SETTPS',zhook_in,zhook_handle)
      NACT  = 0
      N_GROUPS = 0
      LAST_GROUP =0
      DO JOBT=1,NOBTYPMX
        LACT (JOBT) = 0
        GROUP_NO(JOBT) = 0
        GROUP_FIRST(JOBT) = 0
        GROUP_LAST(JOBT) = 0
      ENDDO
!     ----------------------------------------------------------------
!L    2. Set up order of processing in LACT
!     -----------------------------------------------------------------
!     Loop over all AC Obs types known to AC Scheme

      DO JTYPE=1,NOBTYPMX
      THIS_TYPE  = MOD(DEF_AC_ORDER(JTYPE),1000)
      THIS_GROUP = (DEF_AC_ORDER(JTYPE)-THIS_TYPE)/1000

!     This loop determines whether the observation type - THIS_TYPE -
!     is to be used or not from the AC_OBS_TYPES array.

      IF (THIS_TYPE >  0) THEN

!       Use observation type if in namelist array AC_OBS_TYPES

        LUSE = .FALSE.
        DO JOBT=1,NOBTYPMX
          IF (THIS_TYPE == AC_OBS_TYPES(JOBT)) THEN
            LUSE = .TRUE.
          ENDIF
        ENDDO

        IF (LUSE) THEN

!         Set up to process this observation type
          NACT = NACT+1
          LACT(NACT) = THIS_TYPE

!         Group observation types ; Set up GROUP_NO and N_GROUPS

          IF (NACT == 1 .OR. THIS_GROUP /= LAST_GROUP) THEN

!           Start a new group.
            N_GROUPS = N_GROUPS+1
            GROUP_INDEX(N_GROUPS) = THIS_GROUP
            GROUP_FIRST(N_GROUPS) = NACT

          ENDIF
          LAST_GROUP = THIS_GROUP
          GROUP_NO(NACT) = N_GROUPS
          GROUP_LAST(N_GROUPS) = NACT

!         Find this type in MASTER_AC_TYPES ; Abort if not found.
          TYPE_INDEX(NACT)=0
          DO JOBT=1,NOBTYPMX
            IF (THIS_TYPE == MASTER_AC_TYPES(JOBT)) THEN
              TYPE_INDEX(NACT)=JOBT
            ENDIF
          ENDDO
          IF (TYPE_INDEX(NACT) == 0) THEN
            ICODE = 1
            CMESSAGE = 'SETTPS : Observation Type not in Master List ?'
            if(mype == 0)                                               &
     &      PRINT *, ' Observation Type ',THIS_TYPE,                    &
     &                  ' not in Master List ?'
            GO TO 999
          ENDIF

        ENDIF

      ENDIF
      ENDDO   !   End of JTYPE loop.

!     -----------------------------------------------------------------
!     3. Print out list of AC Obs types to be processed
!     -----------------------------------------------------------------

      IF (NACT >  0.AND.mype == 0) THEN

        PRINT '(/,A,/)', ' AC Obs Types to be processed this run'
        PRINT '(A,(T12,13I5))', ' Type  No ',                           &
     &        (LACT(J),J=1,NACT)
        PRINT '(A,(T12,13I5))', ' Group No ',                           &
     &        (GROUP_NO(J),J=1,NACT)
!       PRINT '(A,15I5)', ' Position in Obs Type List    ',
!    +  (TYPE_INDEX(J),J=1,NACT)

        PRINT '(/,A,15I5)', ' Group Number                ',            &
     &  (J,J=1,N_GROUPS)
        PRINT '(A,15I5)', ' Group Index                 ',              &
     &  (GROUP_INDEX(J),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' First Type in Group         ',              &
     &  (GROUP_FIRST(J),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' Last Type in Group          ',              &
     &  (GROUP_LAST (J),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' No of iterations            ',              &
     &  (DEF_NO_ITERATIONS(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' Interval between Iterations ',              &
     &  (DEF_INTERVAL_ITER(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' Ratio of MG Rows to AG Rows ',              &
     &  (DEF_AGRES_ROWS(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' Ratio of MG Pts  to AG Pts  ',              &
     &  (DEF_AGRES_PTS(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' No of analysis levels       ',              &
     &  (DEF_NO_ANAL_LEVS(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' No of weight levels         ',              &
     &  (DEF_NO_WT_LEVS(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15I5)', ' Horizontal Analysis Mode    ',              &
     &  (DEF_MODE_HANAL(GROUP_INDEX(J)),J=1,N_GROUPS)

        PRINT '(/,A,15I5)', ' Group Dep scaling FACTORS in FI  '
        PRINT '(A,15I11)',' Group No ',(J,J=1,N_GROUPS)
        PRINT '(A,15E11.4)', '          ',                              &
     &  (DEF_FI_VAR_FACTOR(GROUP_INDEX(J)),J=1,N_GROUPS)

        PRINT '(/,A,15I5)', ' Nudging Coefficients '
        PRINT '(A,15I11)',' Group No ',(J,J=1,N_GROUPS)
        IF (model_domain == mt_global) THEN
        PRINT '(A,15E11.4)', ' NH       ',                              &
     &  (DEF_NUDGE_NH(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15E11.4)', ' TR       ',                              &
     &  (DEF_NUDGE_TR(GROUP_INDEX(J)),J=1,N_GROUPS)
        PRINT '(A,15E11.4)', ' SH       ',                              &
     &  (DEF_NUDGE_SH(GROUP_INDEX(J)),J=1,N_GROUPS)
        ELSE
        PRINT '(A,15E11.4)', '          ',                              &
     &  (DEF_NUDGE_LAM(GROUP_INDEX(J)),J=1,N_GROUPS)
        END IF

      ELSEIF (NACT == 0.AND.mype == 0) THEN

        PRINT *,' SETTPS : No observation types to process ?'
        ICODE = 1
        CMESSAGE = 'SETTPS : No obs types to process ?'
        GO TO 999

      ENDIF

 999  CONTINUE
      IF (lhook) CALL dr_hook('SETTPS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SETTPS
END MODULE settps_mod
