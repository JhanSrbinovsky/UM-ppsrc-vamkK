! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE NUM_OBS2 -----------------------------------------------
!
!  Programming Standard : UM Doc Paper No 3 ; Version 4 ; 5/2/92
!
!  Project Task : P3
!
!  Purpose : Read in header section of observation file and get
!            number of observations and data values in this obs file.
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: AC Assimilation
MODULE num_obs2_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE NUM_OBS2 (UNIT_NO,P_LEVELS,Q_LEVELS,P_ROWS,ROW_LENGTH, &
     &                     AK,BK,REALHD1,REALHD2,REALHD3,REALHD4,       &
     &                     REALHD5,REALHD6,                             &
! DUMP_AR2 array dimensions
     &  LEN_FIXHD, LEN_INTHD, LEN_REALHD,                               &
     &  LEN1_LEVDEPC, LEN2_LEVDEPC, LEN1_ROWDEPC, LEN2_ROWDEPC,         &
     &  LEN1_COLDEPC, LEN2_COLDEPC, LEN1_FLDDEPC,LEN2_FLDDEPC,          &
     &  LEN_EXTCNST,  LEN_DUMPHIST,                                     &
     &  LEN_CFI1, LEN_CFI2, LEN_CFI3,                                   &
     &  LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS,                               &
! DUMP_AR2 end
     &                     LEN_DATA,NOBTYP,TNOBS,TNDV,                  &
     &                     ICODE,CMESSAGE)
!L----------------------------------------------------------------------
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE IO
      USE Atmos_Max_Sizes
      USE UM_ParParams
      USE comobs_mod, ONLY: nobtypmx
      IMPLICIT NONE

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
!-----------------------------------------------------------------------
!     ARGUMENTS
!-----------------------------------------------------------------------
      INTEGER UNIT_NO      ! IN  : Unit no of Observation file
      INTEGER P_LEVELS     ! IN  : No of model levels
      INTEGER Q_LEVELS     ! IN  : No of model wet levels
      INTEGER P_ROWS       ! IN  : No of model rows
      INTEGER ROW_LENGTH   ! IN  : No of points on row
      REAL    AK(P_LEVELS) ! IN  : Vertical grid
      REAL    BK(P_LEVELS)
      REAL REALHD1,REALHD2 ! IN  : Horizontal grid
      REAL REALHD3,REALHD4
      REAL REALHD5,REALHD6
      INTEGER LEN_DATA     ! IN  : Dimension of data section
      INTEGER NOBTYP       ! OUT : No of observation types
      INTEGER TNOBS        ! OUT : Total no of observations
      INTEGER TNDV         ! OUT : Total no of data values
      INTEGER ICODE        ! OUT : Return code
      CHARACTER(LEN=256) CMESSAGE  !  OUT : Error message if ICODE > 0
!-----------------------------------------------------------------------
!     LEVEL/GRID VARIABLES
!-----------------------------------------------------------------------
      INTEGER OBS_ROW_LENGTH,OBS_P_ROWS,OBS_P_LEVELS,OBS_Q_LEVELS
      REAL OBS_AK(P_LEVELS),OBS_BK(P_LEVELS)
      REAL OBS_LONG_RES,OBS_LAT_RES,OBS_START_LAT
      REAL OBS_START_LONG,OBS_LAT_PSEUDO_POLE,OBS_LONG_PSEUDO_POLE
!-----------------------------------------------------------------------
!     LOCAL VARIABLES
!-----------------------------------------------------------------------
      INTEGER JLEV
      INTEGER START_BLOCK

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
!*L----------------- COMDECK DUMP_LEN --------------------------------

      INTEGER LEN_FIXHD
      INTEGER LEN_INTHD
      INTEGER LEN_REALHD
      INTEGER LEN1_LEVDEPC, LEN2_LEVDEPC
      INTEGER LEN1_ROWDEPC, LEN2_ROWDEPC
      INTEGER LEN1_COLDEPC, LEN2_COLDEPC
      INTEGER LEN1_FLDDEPC, LEN2_FLDDEPC
      INTEGER LEN_EXTCNST
      INTEGER LEN_DUMPHIST
      INTEGER LEN_CFI1, LEN_CFI2, LEN_CFI3
      INTEGER LEN1_LOOKUP_OBS, LEN2_LOOKUP_OBS
!*----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Dynamic allocated arrays
!*L----------------- COMDECK DUMP_DIM ----------------------------------

      INTEGER FIXHD(LEN_FIXHD)
      INTEGER INTHD(LEN_INTHD)
      REAL    REALHD(LEN_REALHD)
      REAL    LEVDEPC(LEN1_LEVDEPC,LEN2_LEVDEPC)
      REAL    ROWDEPC(LEN1_ROWDEPC,LEN2_ROWDEPC)
      REAL    COLDEPC(LEN1_COLDEPC,LEN2_COLDEPC)
      REAL    FLDDEPC(LEN1_FLDDEPC,LEN2_FLDDEPC)
      REAL    EXTCNST(LEN_EXTCNST)
      REAL    DUMPHIST(LEN_DUMPHIST)
      INTEGER CFI1(LEN_CFI1), CFI2(LEN_CFI2), CFI3(LEN_CFI3)
      INTEGER LOOKUP(LEN1_LOOKUP_OBS,LEN2_LOOKUP_OBS)
!*----------------------------------------------------------------------
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('NUM_OBS2',zhook_in,zhook_handle)

!     Go to start of obs file

      CALL SETPOS (UNIT_NO,0,ICODE)

!     Read in headers from obs file
! DEPENDS ON: readhead
      CALL READHEAD (UNIT_NO,                                           &
!*L----------------- COMDECK DUMP_AR1--------------------------------

!      Array addresses and dimensions
     &  FIXHD,    LEN_FIXHD,                                            &
     &  INTHD,    LEN_INTHD,                                            &
     &  REALHD,   LEN_REALHD,                                           &
     &  LEVDEPC,  LEN1_LEVDEPC,LEN2_LEVDEPC,                            &
     &  ROWDEPC,  LEN1_ROWDEPC,LEN2_ROWDEPC,                            &
     &  COLDEPC,  LEN1_COLDEPC,LEN2_COLDEPC,                            &
     &  FLDDEPC,  LEN1_FLDDEPC,LEN2_FLDDEPC,                            &
     &  EXTCNST,  LEN_EXTCNST,                                          &
     &  DUMPHIST, LEN_DUMPHIST,                                         &
     &  CFI1,     LEN_CFI1,                                             &
     &  CFI2,     LEN_CFI2,                                             &
     &  CFI3,     LEN_CFI3,                                             &
     &  LOOKUP,   LEN1_LOOKUP_OBS,LEN2_LOOKUP_OBS,                      &
!*----------------------------------------------------------------------
     &               LEN_DATA,                                          &
     &               START_BLOCK,ICODE,CMESSAGE)
      IF (ICODE >  0) THEN
        IF (lhook) CALL dr_hook('NUM_OBS2',zhook_out,zhook_handle)
        RETURN
      END IF

      OBS_ROW_LENGTH = INTHD(6)         !  No of points on row
      OBS_P_ROWS     = INTHD(7)         !  No of rows
      OBS_P_LEVELS   = INTHD(8)         !  No of model levels
      OBS_Q_LEVELS   = INTHD(9)         !  No of model wet levels

      TNOBS  = INTHD(28)                !  Total no of observations
      TNDV   = INTHD(29)                !  Total no of data values
      NOBTYP = INTHD(32)                !  No of observation types

      DO JLEV=1,P_LEVELS
        OBS_AK(JLEV) = LEVDEPC(JLEV+2,NOBTYP+1) !  Vertical grid
        OBS_BK(JLEV) = LEVDEPC(JLEV+2,NOBTYP+2) !
      ENDDO

      OBS_LONG_RES         = REALHD(1)  !  Horizontal grid
      OBS_LAT_RES          = REALHD(2)  !
      OBS_START_LAT        = REALHD(3)  !
      OBS_START_LONG       = REALHD(4)  !
      OBS_LAT_PSEUDO_POLE  = REALHD(5)  !
      OBS_LONG_PSEUDO_POLE = REALHD(6)  !

!     Check model and acobs file are consistent
! *** NB CHECK_OBS needs rewriting for ND, so commented out
!      IF(LCHECK_GRID)THEN
!      CALL CHECK_OBS (OBS_ROW_LENGTH,OBS_P_ROWS,OBS_P_LEVELS,
!     +                OBS_Q_LEVELS,OBS_AK,OBS_BK,OBS_LONG_RES,
!     +                OBS_LAT_RES,OBS_START_LAT,OBS_START_LONG,
!     +                OBS_LAT_PSEUDO_POLE,OBS_LONG_PSEUDO_POLE,
!     +                P_LEVELS,Q_LEVELS,P_ROWS,ROW_LENGTH,AK,BK,
!     +                REALHD1,REALHD2,REALHD3,REALHD4,REALHD5,REALHD6,
!     +                ICODE,CMESSAGE)
!      IF (ICODE >  0) RETURN
!      END IF

      IF (lhook) CALL dr_hook('NUM_OBS2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE NUM_OBS2
!-----------------------------------------------------------------------
END MODULE num_obs2_mod
