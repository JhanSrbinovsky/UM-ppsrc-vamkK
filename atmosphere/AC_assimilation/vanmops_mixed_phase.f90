! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINES VANMOPS_MIXED_PHASE and CC_TO_RHTOT-------------------
!
!  Purpose : Performs vertical analysis of MOPS cloud data
!            for (2A) mixed phase cloud microphysics scheme
!     When IPASS=1,
!     model field is interpolated to ob locations and increments
!     calculated. Preliminary weight normalisation factors are also
!     calculated.
!
!     When IPASS=2,
!     data density is interpolated to ob locations and
!     final weight normalisation factors are calculated.
!
!  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!  Project Task : P3
!
!  Arguments:---------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE vanmops_mixed_phase_mod

IMPLICIT NONE

CONTAINS

       SUBROUTINE VANMOPS_MIXED_PHASE(                                  &
     &                    KACT,IPASS,THETA,EXNER,CONV_CLD,              &
     &                    CF,PRESSURE,                                  &
     &                    Q,QCL,QCF,LENFLD,NLVFLD,RHCRIT,               &
     &                    OBDATA,CF1PT,CF2PT,CF3PT,CF4PT,               &
     &                    NP1PT,NP2PT,NP3PT,NP4PT,                      &
     &                    QINC,NORMF,OBS_LAT,OBS_LONG,                  &
     &                    OBS_NO,LMISSD,NPTOBT,LENOBT,NDV,LENOB,        &
     &                    NO_ANAL_LEVS,NO_ANAL_VAR,                     &
     &                    P_LEVELS,BL_LEVELS,ICODE,CMESSAGE)

      USE atmos_constants_mod, ONLY: cp

      USE water_constants_mod, ONLY: lc
      USE conversions_mod, ONLY: zerodegc
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx, obs_info
      USE cc_to_rhtot_mod, ONLY: cc_to_rhtot
      USE diago_mod, ONLY: diago
      USE hintmo_mod, ONLY: hintmo
      IMPLICIT NONE
!-----------------------------------------------------------------
!  UM comdecks and functions
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!  Analysis Correction comdecks
!-----------------------------------------------------------------
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
!-----------------------------------------------------------------
      INTEGER LENFLD,NLVFLD,LENOBT,LENOB,NDV,NO_ANAL_LEVS,              &
     &        P_LEVELS,BL_LEVELS,NO_ANAL_VAR
      REAL                                                              &
     &               THETA(LENFLD,P_LEVELS),                            &
     &               EXNER(LENFLD,P_LEVELS),                            &
     &               PRESSURE (LENFLD,P_LEVELS),                        &
     &               CONV_CLD(LENFLD,NLVFLD),                           &
     &               CF(LENFLD,NLVFLD),                                 &
     &               Q   (LENFLD,NLVFLD),                               &
     &               QCL (LENFLD,NLVFLD),                               &
     &               QCF (LENFLD,NLVFLD),                               &
     &               RHCRIT  (NLVFLD),                                  &
     &               OBDATA  (LENOBT,NDV),                              &
     &               QINC(LENOB+1,NO_ANAL_LEVS),                        &
     &               NORMF(LENOB+1,NO_ANAL_LEVS),                       &
     &               OBS_LAT(LENOBT),OBS_LONG(LENOBT),                  &
     &               CF1PT(LENOBT), CF2PT(LENOBT),                      &
     &               CF3PT(LENOBT), CF4PT(LENOBT)
      LOGICAL LMISSD(LENOB+1,NO_ANAL_LEVS)
      INTEGER OBS_NO(LENOBT)
      INTEGER NP1PT(LENOBT),NP2PT(LENOBT),NP3PT(LENOBT)
      INTEGER NP4PT(LENOBT)

      INTEGER KACT,IPASS,NPTOBT
      INTEGER ICODE
      CHARACTER(LEN=256) CMESSAGE

!-INTENT=IN--------------------------------------------------------
!     KACT      - Observation type
!     IPASS     - Calculate incrs and preliminary norm factors for
!                 IPASS=1 and final norm factors for IPASS=2.
!     LENOBT    - no of obs in type
!     LENOB     - no of obs in group
!     NPTOBT    - pointer to present ob type within group
!     NDV       - no of data values in observation
!     NO_ANAL_LEVS - no of analysis levels
!     NO_ANAL_VAR  - no of analysis variables
!     P_LEVELS     - no of model levels
!     EXNER     - diagnostic variable exner pressure
!     CONV_CLD  - convective cloud amount on model levels
!     CF         - layer cloud amount (liq+ice) on model levels
!     PRESSURE   - p on theta levels
!     Q         - model specific humidity      (IPASS=1)
!               - data density on model grid   (IPASS=2)
!     QCF       - model cloud ice
!     QCL       - model cld liq water
!     THETA     - model potential temperature
!     LENFLD    - length of model field
!     NLVFLD    - no of levels in Q field
!     RHCRIT    - critical rh values for cloud formation (fraction)
!     OBDATA    - observed values
!     OBS_LAT   - ob co-lats
!     OBS_LONG  - ob longitudes
!     OBS_NO    - observation numbers
!     NP1PT,NP2PT,NP3PT,NP4PT,
!               - pointers to nearest model grid points for each ob
!     CF1PT,CF2PT,CF3PT,CF4PT,
!               - weighting of model points in bilinear interpolation
!
!-INTENT=INOUT-----------------------------------------------------
!     QINC        -  ob-model humidity increments
!     NORMF       -  normalisation factors
!-INTENT=OUT-----------------------------------------------------
!     LMISSD          -  logical to indicate missing data
!     ICODE,CMESSAGE  - error code and message
!*
!----------------------------------------------------------------------
!*L   Workspace usage
!-----------------------------------------------------------------------
!     DYNAMIC ALLOCATION ON CRAY
      REAL           TL(LENFLD,NLVFLD)

      REAL           CLD_INC (LENOB+1,NO_ANAL_LEVS)
      REAL           CLD_NORM(LENOB+1,NO_ANAL_LEVS)

      REAL           CF_LYROB(LENOBT), CF_CONVB(LENOBT)
      REAL           QB(LENOBT), RHTOT_OB(LENOBT)
      REAL           CFB(LENOBT), TLB(LENOBT), QCLB(LENOBT)
      REAL           QCFB(LENOBT), PB(LENOBT)
      REAL           QSAT_ICEB(LENOBT), QSAT_WATB(LENOBT)

!     TL         - model liquid water temperature
!     CLD_INC    - obs - model cloud fraction (layer+convective)
!     CLD_NORM   - normalisation array needed for input to DIAGO call
!     CF_LYROB   - target layer cloud fraction derived from cloud ob
!     CF_CONVB   - model convective cloud fraction interp'd to obs pts
!     QB         - model humidity (or obs density) at obs point
!     RHTOT_OB   - target rhtot derived from CF_LYROB
!     CFB        - CF at obs points
!     TLB        - TL at obs points
!     QCLB,QCFB  - QCL,QCF at obs points
!     PB         - model level pressure at obs pt
!     QSAT_ICEB  - QSAT wrt ice at obs pt
!     QSAT_WATB  - QSAT wrt water at obs pt
!*
!----------------------------------------------------------------------
!*L   External subroutine calls
!-----------------------------------------------------------------------
      EXTERNAL TIMER
      EXTERNAL QSAT,QSAT_WAT
!*
!----------------------------------------------------------------------
!     Define local variables
!----------------------------------------------------------------------
      REAL    TQSATICE,TQSATWAT,DTQSAT
      INTEGER KTYPE,NOBIPT,NOBLPT,NEROPT,JK,JOB,J,ERROR

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!     KTYPE    -   Observation type
!     NOBIPT   -   Pointer to data values within OBDATA array
!     NOBLPT   -   Pointer to obs levels within OBDATA array
!     NEROPT   -   Pointer to obs errors within OBDATA array
!     JK       -   Loop counter in loops over levels
!     JOB      -         "      in loops over obs
!     J        -         "      in loops over points
!     ERROR    -   error code
!     TQSATICE -   new cloud likely to be ice below this temp
!                  (see THOMO in C_LSPMIC)
!     TQSATWAT -   new cloud likely to be liquid above this temp
!                  (see TNUC in C_LSPMIC)
!     DTQSAT   -   difference between two previous temperatures
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('VANMOPS_MIXED_PHASE',zhook_in,zhook_handle)
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('VANMMP  ',3)

!L
!L*** 1. PRELIMINARIES
!L       -------------

      KTYPE = LACT (KACT)

      IF (KTYPE == 406) THEN
        NOBIPT = 1
      ENDIF
      NEROPT = OBS_INFO % NERLEV1(KACT) - OBS_INFO % NDVHDR

!     Define temps used in choosing qsat value when creating cloud
!     (These temps are the same as the homogeneous and heterogeneous
!      nucleation thresholds in the mixed phase precip scheme,
!      but need not be, hence separate definition)
      TQSATICE = ZERODEGC - 40.0
      TQSATWAT = ZERODEGC - 10.0
      DTQSAT   = TQSATWAT - TQSATICE

      IF (KTYPE == 406) THEN
!     =================
      IF(IPASS == 1) THEN
!L*** 2.1 Current model cloud fraction is consistent with
!L        latest humidity, since ls_cld called just before AC.
!L        But need to get current TL for use in qsat calcns.


!     2.1.2 calculate TL
         DO JK=1,NLVFLD
           DO J=1,LENFLD
            TL(J,JK) = THETA(J,JK)*EXNER(J,JK)- LC * QCL(J,JK) /CP
           END DO
         END DO

!L*** 2.2  HORIZONTAL INTERP OF FIELDS ON MODEL GRID


        ENDIF ! IPASS == 1

! begin loop over levels
      DO JK=1,NO_ANAL_LEVS

!       PUT Q OR OBS DENSITY INTERPOLATED TO OB POSITIONS IN QB
        CALL HINTMO( Q(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,                   &
     &               NP1PT,NP2PT,NP3PT,NP4PT,                           &
     &               LENFLD,1,LENOBT,QB,ICODE,CMESSAGE)
        IF (ICODE >  0) GO TO 9999

!L*** 2.3  CALCULATION OF INCRS AND WEIGHT FACTORS FOR OBS INC

!     2.3.1  FIRST SET UP BIT ARRAY FOR NO OBSERVATIONAL DATA
      DO JOB=1,LENOBT
        LMISSD(NPTOBT+JOB,JK) = OBDATA(JOB,NOBIPT+JK-1) == OBS_INFO % MISSD &
                           .OR. OBDATA(JOB,NEROPT+JK-1) == OBS_INFO % MISSD
      ENDDO   !JOB

      IF (IPASS == 1) THEN

!      2.3.2 put conv_cld interp'd to ob locations in CF_CONVB
       CALL HINTMO( CONV_CLD(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,             &
     &              NP1PT,NP2PT,NP3PT,NP4PT,                            &
     &              LENFLD,1,LENOBT,CF_CONVB,ICODE,CMESSAGE)
       IF(ICODE >  0) GO TO 9999

!      2.3.3 calculate target layer cloud cover in CF_LYROB
!            based on total cloud = conv_cld + ls_cld *(1-conv_cld)
!            as assumed in radiation scheme.
        DO JOB=1,LENOBT
          IF(CF_CONVB(JOB) /=  1.) THEN
            CF_LYROB(JOB)=(OBDATA(JOB,NOBIPT+JK-1)-CF_CONVB(JOB))/      &
     &                                           (1.-CF_CONVB(JOB))
          ELSE
!         set target layer cloud of zero where conv cover = 1
            CF_LYROB(JOB)=0.
          ENDIF
!         set target layer cloud of zero where conv cover > MOPS
!         (will apply also where MOPS data is missing data value<0,
!          though obs cloud not used then)
          CF_LYROB(JOB)=MAX(0.,CF_LYROB(JOB) )
        ENDDO

!      2.3.4 convert target layer cloud to rhtot
        CALL CC_TO_RHTOT (CF_LYROB,LENOBT,RHCRIT(JK),RHTOT_OB,JK,       &
     &                                             BL_LEVELS)

!      2.3.5 interp model b'ground fields to ob pts for incr calcns
!            calculate derived quantities

!        2.3.5.1 put TL    interp'd to ob pts in TLB
          CALL HINTMO( TL(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,                &
     &                 NP1PT,NP2PT,NP3PT,NP4PT,                         &
     &                 LENFLD,1,LENOBT,TLB,ICODE,CMESSAGE)
          IF (ICODE >  0) GO TO 9999

!        2.3.5.2 put QCL   interp'd to ob pts in QCLB
          CALL HINTMO( QCL(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,               &
     &                 NP1PT,NP2PT,NP3PT,NP4PT,                         &
     &                 LENFLD,1,LENOBT,QCLB,ICODE,CMESSAGE)
          IF (ICODE >  0) GO TO 9999

!        2.3.5.3 put QCF   interp'd to ob pts in QCFB
          CALL HINTMO( QCF(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,               &
     &                 NP1PT,NP2PT,NP3PT,NP4PT,                         &
     &                 LENFLD,1,LENOBT,QCFB,ICODE,CMESSAGE)
          IF (ICODE >  0) GO TO 9999

!        2.3.5.4 put TOTAL LAYER CLOUD FRACTION interp'd to
!                                                   ob pts in CFB
          CALL HINTMO( CF(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,                &
     &                 NP1PT,NP2PT,NP3PT,NP4PT,                         &
     &                 LENFLD,1,LENOBT,CFB,ICODE,CMESSAGE)
          IF (ICODE >  0) GO TO 9999

!        2.3.5.5 put MODEL LEVEL PRESSURE interp'd to ob pts in PB
      CALL HINTMO( PRESSURE(1,JK),CF1PT,CF2PT,CF3PT,CF4PT,              &
     &                 NP1PT,NP2PT,NP3PT,NP4PT,                         &
     &                 LENFLD,1,LENOBT,PB,ICODE,CMESSAGE)
          IF (ICODE >  0) GO TO 9999

!        2.3.5.6 get QSAT(TL,P) wrt ice in QSAT_ICEB
! DEPENDS ON: qsat
           CALL QSAT (QSAT_ICEB, TLB, PB, LENOBT)

!        2.3.5.7 get QSAT(TL,P) wrt water in QSAT_WATB
! DEPENDS ON: qsat_wat
           CALL QSAT_WAT (QSAT_WATB, TLB, PB, LENOBT)

!        2.3.5.8 update RHTOT_OB where ice present
!                covers eqns 12 & 13 in working paper
           DO JOB=1,LENOBT
             IF (QCFB(JOB)  >   0.0) THEN
               RHTOT_OB(JOB) = RHTOT_OB(JOB) *                          &
     &          (QCLB(JOB) + QCFB(JOB)*                                 &
     &                          QSAT_ICEB(JOB)/QSAT_WATB(JOB) )         &
     &         /( QCLB(JOB) + QCFB(JOB) )
             ENDIF
           END DO

!      2.3.6 calculate increments to Q and normalisation factors
          DO JOB=1,LENOBT
            IF(.NOT.LMISSD(NPTOBT+JOB,JK)) THEN
!             Data is valid
!             Calculate incrs to Q
              IF( CFB(JOB) == CF_LYROB(JOB) ) THEN
!               model and ob cloud equal, so no humidity incr
                QINC(NPTOBT+JOB,JK) = 0.0
              ELSEIF(CFB(JOB) >  0.0) THEN
!               eqns 7 & 8 in working paper
                QINC(NPTOBT+JOB,JK) = RHTOT_OB(JOB)*QSAT_WATB(JOB) -    &
     &                                ( QB(JOB)+QCLB(JOB)+QCFB(JOB) )
              ELSEIF(CFB(JOB) == 0.0 .AND. CF_LYROB(JOB) >  0.0) THEN
!               need to create cloud in model,
!               refer to section 2.5 of working paper
                IF    (TLB(JOB) >  TQSATWAT) THEN
                  QINC(NPTOBT+JOB,JK) = RHTOT_OB(JOB)                   &
     &                                 *QSAT_WATB(JOB) - QB(JOB)
                ELSEIF(TLB(JOB) <  TQSATICE) THEN
                  QINC(NPTOBT+JOB,JK) = RHTOT_OB(JOB)                   &
     &                                 *QSAT_ICEB(JOB) - QB(JOB)
                ELSE
                  QINC(NPTOBT+JOB,JK)  = RHTOT_OB(JOB) *                &
     &              (QSAT_ICEB(JOB)*(1-(TLB(JOB)-TQSATICE)/DTQSAT) +    &
     &               QSAT_WATB(JOB)*(TLB(JOB)-TQSATICE)/DTQSAT )        &
     &              - QB(JOB)
                ENDIF
              ENDIF
!             check humidity incr reasonable, or reset to zero
!             (humidity and cloud incrs must have same sign)
              IF( QINC(NPTOBT+JOB,JK)*                                  &
     &                      (CF_LYROB(JOB)-CFB(JOB)) <  0.0 ) THEN
                QINC(NPTOBT+JOB,JK) = 0.0
              ENDIF
!             Calculate normalisation factor
              NORMF(NPTOBT+JOB,JK) = 1.0 /                              &
     &            ( OBDATA(JOB,NEROPT+JK-1)*OBDATA(JOB,NEROPT+JK-1) )
            ELSE
!             Missing data
              QINC(NPTOBT+JOB,JK)  = 0.0
              NORMF(NPTOBT+JOB,JK) = 0.0
            ENDIF ! test on missing data
          ENDDO  !  JOB

!     2.3.7  store cloud increments (oktas) for diagnostics
      IF(LDIAGAC) THEN
        DO JOB=1,LENOBT
          IF(.NOT.LMISSD(NPTOBT+JOB,JK)) THEN
            CLD_INC (JOB,JK) = 8*(OBDATA(JOB,NOBIPT+JK-1)-CF_CONVB(JOB) &
     &                           -CFB(JOB)*(1. - CF_CONVB(JOB) )  )
            CLD_NORM(JOB,JK) = 1.0
          ELSE
            CLD_INC (JOB,JK) = 0.0
            CLD_NORM(JOB,JK) = 0.0
          ENDIF
        END DO

      ENDIF


      ELSEIF (IPASS == 2) THEN
          DO JOB=1,LENOBT
            IF (.NOT.LMISSD(NPTOBT+JOB,JK)) THEN
              NORMF(NPTOBT+JOB,JK) = NORMF(NPTOBT+JOB,JK) /             &
     &                              ( QB(JOB) + 1.0 )
            ELSE
              NORMF(NPTOBT+JOB,JK) = 0.0
            ENDIF
          ENDDO !  JOB
      ENDIF ! IPASS


      ENDDO   ! JK

!L*** 3.  Diagnostics
!         ===========
!     measure fit of model to MOPS cloud cover
      IF(LDIAGAC .AND. IPASS == 1 ) THEN
          if(mype == 0)                                                 &
              PRINT '(/,'' DIAGO called from VANMOPS_MIXED_PHASE'//&
              ' - Obs Type'',I5,T50,''No of obs '',I6)', KTYPE,LENOBT

          CALL DIAGO ('MULTI-LEVEL',KTYPE,6,                            &
     &                CLD_INC,CLD_NORM,OBS_LAT,OBS_LONG,LMISSD,         &
     &                LENOBT,LENOBT,0,NO_ANAL_LEVS,NO_ANAL_VAR)


      ENDIF

      ENDIF  ! KTYPE
!     =====

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('VANMMP ',4)

 9999 CONTINUE
      IF (lhook) CALL dr_hook('VANMOPS_MIXED_PHASE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE VANMOPS_MIXED_PHASE



END MODULE vanmops_mixed_phase_mod
