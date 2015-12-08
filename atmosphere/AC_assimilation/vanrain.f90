! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE VANRAIN-------------------------------------------------
!LL
!LL  Purpose : Performs 'vertical' analysis for precip rate/phase data.
!LL
!LL     When IPASS=1,
!LL     model field is interpolated to ob locations and increments
!LL     calculated. Preliminary weight normalisation factors are also
!LL     calculated.
!LL
!LL     When IPASS=2,
!LL     data density is interpolated to ob locations and
!LL     final weight normalisation factors are calculated.
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LLEND------------------------------------------------------------------
!
!*L  Arguments:---------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE vanrain_mod

IMPLICIT NONE

CONTAINS

       SUBROUTINE VANRAIN (KACT,IPASS,LS_RAIN,LS_SNOW,                  &
     &                     CONV_RAIN,CONV_SNOW,LENFLD,                  &
     &                     OBDATA,CF1PT,CF2PT,CF3PT,CF4PT,              &
     &                     NP1PT,NP2PT,NP3PT,NP4PT,                     &
     &                     PRINC,NORMF,OBS_LAT,OBS_LONG,                &
     &                     OBS_NO,LMISSD,NPTOBT,                        &
     &                     P_LEVELS,LENOBT,NDV,LENOB,                   &
     &                     NO_ANAL_LEVS,NO_ANAL_VAR,                    &
     &                     ICODE,CMESSAGE)
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Atmos_Max_Sizes
      USE UM_ParParams
      USE comobs_mod, ONLY: nobtypmx, obs_info
      USE diagopr_mod, ONLY: diagopr
      USE hintmo_mod, ONLY: hintmo
      IMPLICIT NONE

!-----------------------------------------------------------------
!  Analysis Correction comdecks
!
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
!--------------------------------------------------------------------
!LCOMDECK COMMG
!L-------------
      REAL          DLAT,DLONG,XLATN,XLONGW
      REAL          ELFPLAT,ELFPLON
      COMMON/COMMG/ DLAT,DLONG,XLATN,XLONGW,ELFPLAT,ELFPLON
!--------------------------------------------------------------------
!
! ----------------------------------------------------------------------
      INTEGER LENFLD,LENOBT,LENOB,NDV
      INTEGER P_LEVELS
      INTEGER KACT,IPASS,NPTOBT,NO_ANAL_LEVS,NO_ANAL_VAR
      REAL                                                              &
     &               LS_RAIN(LENFLD),LS_SNOW(LENFLD),                   &
     &               CONV_RAIN(LENFLD),CONV_SNOW(LENFLD),               &
     &               OBDATA(LENOBT,NDV),                                &
     &               PRINC(LENOB),NORMF(LENOB),                         &
     &               OBS_LAT(LENOBT),OBS_LONG(LENOBT),                  &
     &               CF1PT(LENOBT),CF2PT(LENOBT),CF3PT(LENOBT),         &
     &               CF4PT(LENOBT)
!
      LOGICAL LMISSD(LENOB+1,NO_ANAL_LEVS)
!
      INTEGER OBS_NO(LENOBT)
      INTEGER NP1PT(LENOBT),NP2PT(LENOBT),NP3PT(LENOBT)
      INTEGER NP4PT(LENOBT)
!
      INTEGER ICODE
      CHARACTER(LEN=256) CMESSAGE
!
!-INTENT=IN--------------------------------------------------------
!     KACT      - Observation type
!     IPASS     - Calculate incrs and preliminary norm factors for
!                 IPASS=1 and final norm factors for IPASS=2.
!     LENOBT    - no of obs in type
!     LENOB     - no of obs in group
!     NPTOBT    - pointer to present ob type within group
!     NDV       - no of data values in observation
!     NO_ANAL_LEVS - no of analysis levels
!     P_LEVELS  - no of model level
!     LS_RAIN   - }  model precipitation rates (IPASS=1) for
!     LS_SNOW   - }  large-scale and convective rain and snow
!     CONV_RAIN - }  or data density on model grid (IPASS=2)
!     CONV_SNOW - }
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
!     PRINC        -  ob-model increments
!     NORMF        -  normalisation factors
!-INTENT=OUT-----------------------------------------------------
!     LMISSD          - logical to indicate missing data
!     ICODE,CMESSAGE  - error code and message
!*
!---------------------------------------------------------------------
!*L   Workspace usage
!-----------------------------------------------------------------------
!     DYNAMIC ALLOCATION ON CRAY
      REAL           INTPR(LENOBT),TOTPR(LENFLD)
      REAL           PHSINT(LENOBT),PHSFLD(LENFLD)
      REAL           WLAT(LENOBT),WLON(LENOBT)
!     INTPR   - model rate field interpolated to ob locations
!     TOTPR   - total rate field, =sum of LS/CONV_RAIN/SNOW
!     WLAT    - Equatorial latitude of obs in degrees
!     WLON    - Equatorial longitude of obs in degrees
!*
!----------------------------------------------------------------------
!*L   External subroutine calls
!-----------------------------------------------------------------------
      EXTERNAL TIMER
      EXTERNAL EQTOLL
!
!----------------------------------------------------------------------
!     Define local variables
!-----------------------------------------------------------------------
      INTEGER KTYPE,NLEV,NOBIPT,NEROPT,JOB,JPT
      INTEGER JRAD
      REAL  OB_TO_RADAR, OB_TO_RADAR_MIN(LENOBT)
      REAL  F2

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!     KTYPE    -   Observation type
!     NLEV     -   No of observed levels
!     NOBIPT   -   Pointer to data values within OBDATA array
!     NEROPT   -   Pointer to obs errors within OBDATA array
!     JOB      -   Loop counter in loops over obs
!     JPT      -   Loop counter in loops over model points
!     JRAD     -   Loop counter over radars
!     OB_TO_RADAR      - Distance (km) calculated from ob to radar
!     OB_TO_RADAR_MIN  - Min value of OB_TO_RADAR
!----------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook('VANRAIN',zhook_in,zhook_handle)
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('VANRAIN ',3)
!
      KTYPE = LACT(KACT)

      NLEV  = OBS_INFO % NOBLEV(KACT)
      IF (KTYPE == 506) THEN
        NOBIPT = 1
      ENDIF
      NEROPT = OBS_INFO % NERLEV1(KACT) - OBS_INFO % NDVHDR
!
!   COMBINE model precipitation rates to give TOTPR
!

      IF (IPASS  ==  1) THEN
        DO JPT=1,LENFLD
          TOTPR(JPT)=  LS_RAIN(JPT)   + LS_SNOW(JPT)                  &
     &                 + CONV_RAIN(JPT) + CONV_SNOW(JPT)
        END DO
      ELSEIF (IPASS  ==  2) THEN   ! TOTPR set to obs density
        DO JPT=1,LENFLD
          TOTPR(JPT) = LS_RAIN(JPT)
        END DO
      ENDIF
!
!
!L*** 1.   OBSERVATIONS GIVING DATA FOR PRECIP RATE/PHASE
!L         ----------------------------------------------
!
!L*** 1.1  HORIZONTAL INTERP OF RATE FIELD ON MODEL GRID
!
!
        CALL HINTMO(TOTPR,CF1PT,CF2PT,CF3PT,CF4PT,                      &
     &              NP1PT,NP2PT,NP3PT,NP4PT,                            &
     &              LENFLD,1,LENOBT,INTPR,                              &
     &              ICODE,CMESSAGE)
        IF (ICODE >  0) GO TO 9999

!
!L*** 1.2  VERTICAL INTERP/PROCESSING OF MODEL
!          NONE.
!L*** 1.3  SUBTRACTION, CALCULATION OF PRELIM NORM FACTOR (IPASS=1)
!L                      CALCULATION OF FINAL NORM FACTOR (IPASS=2)
!
      IF (IPASS == 1) THEN
!       PR Increment = PR Observation - Background PR
        DO JOB=1,LENOBT
          PRINC(NPTOBT+JOB) = OBDATA(JOB,NOBIPT) - INTPR(JOB)
        END DO
      ENDIF

!
!     SET WEIGHT FACTOR TO ZERO WHERE OB. VALUE OR ERROR IS FLAGGED
      DO JOB=1,LENOBT
        LMISSD(NPTOBT+JOB,1) = OBDATA(JOB,NOBIPT) == OBS_INFO % MISSD .OR.  &
     &                         OBDATA(JOB,NEROPT) == OBS_INFO % MISSD
      END DO
!
      IF (IPASS == 1) THEN
        IF (MDATADFN == 1) THEN
          DO JOB=1,LENOBT
            IF (LMISSD(NPTOBT+JOB,1)) THEN
              NORMF(NPTOBT+JOB) = 0.0
            ELSE
              NORMF(NPTOBT+JOB) = 1.0 /                                 &
     &                  SQRT(OBDATA(JOB,NEROPT)*OBDATA(JOB,NEROPT)+1.0)
            ENDIF
          END DO !job
        ELSEIF (MDATADFN == 2) THEN
          DO JOB=1,LENOBT
            IF (LMISSD(NPTOBT+JOB,1)) THEN
              NORMF(NPTOBT+JOB) = 0.0
            ELSE
              NORMF(NPTOBT+JOB) = 1.0 /                                 &
     &                        ( OBDATA(JOB,NEROPT)*OBDATA(JOB,NEROPT) )
            ENDIF
          END DO ! job
        ENDIF
!
!
      ELSEIF (IPASS == 2) THEN
        IF (MDATADFN == 1) THEN
          DO JOB=1,LENOBT
            IF (LMISSD(NPTOBT+JOB,1)) THEN
              NORMF(NPTOBT+JOB) = 0.0
            ELSE
              NORMF(NPTOBT*JOB) = NORMF(NPTOBT+JOB)*NORMF(NPTOBT+JOB) / &
     &         ( INTPR(JOB)*NORMF(NPTOBT+JOB) +1.0 -                    &
     &           NORMF(NPTOBT+JOB)*NORMF(NPTOBT+JOB) )
            ENDIF
          END DO ! job
        ELSEIF (MDATADFN == 2) THEN
          DO JOB=1,LENOBT
            IF (LMISSD(NPTOBT+JOB,1)) THEN
              NORMF(NPTOBT+JOB) = 0.0
            ELSE
              NORMF(NPTOBT+JOB) = NORMF(NPTOBT+JOB)/( INTPR(JOB)+1.0 )
            ENDIF
          END DO ! job
        ENDIF
!
      ENDIF
!
!L*** 1.4  VERTICAL INTERP/PROCESSING OF INCS. CALC'N OF WEIGHT FACTORS.
!          NONE.
!
!     DIAGNOSTICS ON FIT OF MODEL PRECIP RATE/PHASE TO OBS
      IF (LDIAGAC .AND. (IPASS == 1) ) THEN

        CALL DIAGOPR (lenobt,INTPR,OBDATA,LENOBT,LENOB,NDV,             &
     &                LMISSD,NOBIPT,NPTOBT,NO_ANAL_LEVS)
      ENDIF

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('VANRAIN ',4)

 9999 CONTINUE
      IF (lhook) CALL dr_hook('VANRAIN',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE VANRAIN
END MODULE vanrain_mod
