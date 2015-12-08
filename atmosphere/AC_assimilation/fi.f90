! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE FI -------------------------------------------------
!LL
!LL  Purpose :
!LL
!LL     Perform the horizontal spreading of observational increments.
!LL     There are three stages:
!LL     1. calculate the horizontal scale to be used
!LL     2. spread the observation weights to get normalization factors
!LL     3. spread the normalized increments.
!LL     Spreading is performed using the adjoint of the model->ob interp
!LL     followed by a filter on the model grid using RF.
!LL
!LL     This is an alternative to the method using HORINF.
!LL
!LL NB: This routine is now only used for MOPS rain and since 
!LL NPASS_RF=0 FILT=F so recursive filtering routines are no longer
!LL called. For this reason the code has not been upgraded to
!LL handle variable resolution grids in this mode.
!LL
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Logical components covered:
!LL
!LL  Project Task : P3
!LL
!LL  Documentation: UM Doc Paper 30 section 6.4
!LL
!LLEND------------------------------------------------------------------

!*L  Arguments:---------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE fi_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE FI(IPASS,TIMESTEP_NO,TIMESTEP,NAVT,JGROUP,             &
     &  NO_ANAL_VAR,LENOB,NO_ANAL_LEVS,NO_WT_LEVS,FIRST_LEVEL,          &
     &  LENMG,ROW_LENGTH,P_ROWS,OBS_NO,                                 &
     &  OBS_LAT,OBS_LONG,OBS_TIME,NORMF,OBS_INCR,CFPT,NPPT,             &
     &  IP_SCALE,MODEL_INCR,PSTAR,                                      &
     &  P_FIELD,P_LEVELS,                                               &
     &  ICODE,CMESSAGE)

       USE um_input_control_mod, ONLY : model_domain
       USE domain_params, ONLY: mt_global

       USE earth_constants_mod, ONLY: earth_radius

       USE yomhook, ONLY: lhook, dr_hook
       USE parkind1, ONLY: jprb, jpim
       USE Atmos_Max_Sizes
       USE UM_ParParams
       USE comobs_mod, ONLY: nobtypmx
       USE rfcsl_mod, ONLY: rfcsl
       USE rfvsg_mod, ONLY: rfvsg
       USE rfvsl_mod, ONLY: rfvsl
       USE rfvvg_mod, ONLY: rfvvg
       USE rfvvl_mod, ONLY: rfvvl
       IMPLICIT NONE

      INTEGER   IPASS,TIMESTEP_NO
      INTEGER   JGROUP,FIRST_LEVEL
      INTEGER   NAVT,NO_ANAL_VAR,NO_ANAL_LEVS,NO_WT_LEVS,IP_SCALE
      INTEGER   LENOB,LENMG,ROW_LENGTH,P_ROWS
      REAL      TIMESTEP
      INTEGER   ICODE  ! 0 if OK
      CHARACTER(LEN=256) CMESSAGE

      INTEGER OBS_NO(LENOB),NPPT(LENOB+1,4)
      REAL                                                              &
     &     OBS_LAT(LENOB),OBS_TIME(LENOB),CFPT(LENOB+1,4),              &
     &     OBS_LONG(LENOB),NORMF(LENOB+1,NO_ANAL_LEVS),                 &
     &     OBS_INCR(LENOB+1,NO_ANAL_LEVS,NO_ANAL_VAR),                  &
     &     MODEL_INCR(LENMG,NO_ANAL_LEVS,IP_SCALE),                     &
     &     PSTAR(ROW_LENGTH*P_ROWS)

      INTEGER                                                           &
     & P_FIELD,P_LEVELS !IN PTS ON U GRID;P GRID ,NO OF LEVELS


!-INTENT=IN--------------------------------------------------------
!     IPASS : mode
!            1- do stages 1. calculate scale
!                       & 2. spread NORMF
!            0- do stage  3. spread INCR
!     TIMESTEP_NO: current model timestep (TIMESTEP_NO=1,.....N)
!     NAVT: observation variable type
!     JGROUP : ob group code
!     NO_ANAL_LEVS: number of levels analysed ( >=  NO_WT_LEVS)
!     FIRST_LEVEL : model level corresponding to analysis level 1
! N.B. current AC code assumes this is 1; this could be relaxed later
!     NO_WT_LEVS: number of levels for which NORMF is needed
!     LENOB: number of observations in group
!     LENMG : NUMBER OF POINTS IN FIELD ON model GRID
!     ROW_LENGTH : length of model grid row
!     P_ROWS : no of row on model p grid (assumed > U_ROWS)
!     TIMESTEP   : is model timestep in secs
!     OBS_NO : list of ob numbers for this group within ACOBS file
!     OBS_LAT  :  co-lat of obs
!     OBS_LONG : longitude of obs
!     OBS_TIME : time of obs rel to assm start
!     NORMF    : normalization factors for obs (from routine VERTANL)
!     OBS_INCR : increment (ob-model) at ob pt
!     CFPT : model->ob interpolation coeffs
!     NPPT : model->ob interpolation pointers
!     PSTAR : surface pressure (used in DIVFILT for mass weighting)
!     IP_SCALE : 3rd index to MODEL_INCR, pointing to space for scale
!     IP_SCALE=2 when doing a scalar, =3 when doing a vector analysis.
!-INTENT=OUT---------------------------------------------------------
!    stage 1:
!     MODEL_INCR(,,IP_SCALE)     : horizontal scale for stages 2 & 3
!     ( MODEL_INCR(,,1:2) )      : is used as work space in stage 1
!    stage 2:
!     MODEL_INCR(,,1)            : data density field D for norm.factor
!                                : (only first NO_WT_LEVS levels done)
!    stage 3: (second call to FI)
!     MODEL_INCR(,,1:NO_ANAL_VAR): increments
!-INTENT=IN   (second call to FI)
!     MODEL_INCR(,,IP_SCALE)     : horizontal scale
!*---------------------------------------------------------------------

!*L  Workspace usage:-------------------------------------------------
      REAL COS_LAT(P_ROWS) ! sin of co-latitude of model grid rows
      REAL S_INCR(LENOB,NO_ANAL_LEVS,2) ! incs at obs to be scattered
      REAL TF(LENOB) ! Time Factor for each ob  as in HORINF
      REAL S(LENOB)  ! Scale for each ob        as in HORINF
      REAL WRK1(LENOB)
      INTEGER IWK1(LENOB)
!
!*---------------------------------------------------------------------


!----analysis correction comdecks-------------------------------------
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
!*---------------------------------------------------------------------

!-- constants comdecks-----------------------------------------------
! ---------------------------------------------------------------------

!----------------------------------------------------------------------
!    Define local variables
      INTEGER JROW,JOB,JOBT,J,JLEV,JC ! loop indices
      INTEGER I,ITYPE
      INTEGER IROWS ! number of rows on model grid for current variable
      LOGICAL ILUV  ! .true. if winds
      LOGICAL FILT  ! .true. if filtering on model grid is reqd.
      INTEGER M_GRID ! 0=Ltd Area, 1=pts at pole, 2=staggd from pole
      REAL    TIMEB,TIMEA,RADINF
      REAL    CSCALE_START,CSCALE_OBTIME,CSCALE_END
      INTEGER FIRST_OBS,LAST_OBS,FIRST_TYPE,LAST_TYPE
      REAL Z,ZZ,ZINTEGRAL,ZIOA
      REAL FI_VAR_FACTOR  ! to change SCALE to match tuning to HORINF
      REAL, ALLOCATABLE :: DELTA_LAT(:) 
      REAL, ALLOCATABLE :: DELTA_LONG(:)

      CHARACTER(LEN=10) LABEL

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------

!L***     PRELIMINARIES
!         -------------

      IF (lhook) CALL dr_hook('FI',zhook_in,zhook_handle)
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('FI      ',3)
      
      ICODE=0
      FILT = NPASS_RF >  0
      IF (FILT) THEN
      ZIOA=(25.1327-18.2655/NPASS_RF)/(DLATMG*DLONGMG)
!     ZIOA is a constant in the formula for correlation area/box area
!     ZIOA=v (6.5.15) / (dlat*dlong)
!     For NPASS_RF =2, v=16.; NPASS_RF=infinity, v=8*PI.
!     v for other NPASS_RF is got from an approximate interpolation.
!     the COS_LAT(JROW) term in the grid box area is done in the loops
!
!     HORINF truncated the influence area at RADINF corr scales
!     and uses an EXP(-r/SCALE) correlation function for wind -
!     FI_VAR_FACTOR compensates for this so FI gives similar
!     results to HORINF with the default values for scale etc..
      FI_VAR_FACTOR=DEF_FI_VAR_FACTOR(GROUP_INDEX(JGROUP))
      FI_VAR_FACTOR=SQRT(FI_VAR_FACTOR)
      ELSE
        ZIOA=1.
        FI_VAR_FACTOR=1.
      ENDIF  !  FILT

      IF(FILT) THEN ! these only required if FILT=T so using recursive filters
      ILUV=NAVT == 3
      IF(ILUV)THEN  ! model U grid
        Z=ROW1MGUV
        IROWS=P_ROWS-1
        M_GRID=2
      ELSE          ! model P grid
        Z=ROW1MGTH
        IROWS=P_ROWS
        M_GRID=1
      ENDIF
!L    INITIALIZE COS_LAT
      DO      J=1,IROWS
        COS_LAT(J)=SIN(Z)
        Z=Z+DLATMG
      ENDDO ! J
      ENDIF
      M_GRID=0              ! no pole in ltd area version

!**** define time-factors and scales for each observation
!     The following code for section 2.4 was copied from HORINF (2.6)
!     (preliminary version 1/6/92) and updated 16/2/93
!     so that FI behaves as similarly as possible to HORINF for testing
!
!L2.4 TF which is time factor R(dT) in 3.18 and 3.21
!L##  AND S(DT) BEING THE CORRELATION SCALE AS IN EQUATION 3.19

!     for synoptic insertion mode treat all obs times as though at
!     the nearest synoptic hour I.E. T(OB) = T+0
      IF(LSYN)THEN
       DO 2400 JOB=1,LENOB
       OBS_TIME(JOB) = OBTIME_NOM
!      value of OBTIME_NOM is relative to start of assimilation
2400   ENDDO ! JOB
      END IF

      FIRST_OBS = 1

!     Get first and last types for group being processed.
      FIRST_TYPE = GROUP_FIRST(JGROUP)
      LAST_TYPE  = GROUP_LAST (JGROUP)

      DO JOBT = FIRST_TYPE,LAST_TYPE   !  Loop over types in group

        IF (LENACT(JOBT) >  0) THEN   !  Any obs for this type

          LAST_OBS = FIRST_OBS + LENACT(JOBT) - 1

!         Get analysis parameters for this type
          ITYPE = TYPE_INDEX(JOBT)

          TIMEB         = DEF_TIMEB(ITYPE)
          TIMEA         = DEF_TIMEA(ITYPE)
          RADINF        = DEF_RADINF(ITYPE)
          CSCALE_START  = DEF_CSCALE_START(ITYPE)
          CSCALE_OBTIME = DEF_CSCALE_OBTIME(ITYPE)
          CSCALE_END    = DEF_CSCALE_END(ITYPE)

      DO JOB = FIRST_OBS,LAST_OBS
!       WRK1= time relative obs time (+ve before and -ve after)
        WRK1(JOB) = OBS_TIME(JOB) - (TIMESTEP_NO-1)*TIMESTEP/60.0
        TF(JOB)=0.0
        S(JOB)=1.0
!       note that these initialised values will not change for an ob
!       with relative time outside the insertion period (as required)
      ENDDO ! JOB

      DO JOB = FIRST_OBS,LAST_OBS
!     calculate ramp functions between start of insertion & ob time
!     data on edge of insertion period not fetched by GETOBS, so may
!     be discounted here too.
      IF (WRK1(JOB) >  0.0 .AND. WRK1(JOB) <= TIMEB) THEN
       IF(MRAMPFN == 2)THEN
! quadratic function
       TF(JOB) =                                                        &
     & WRK1(JOB)*WRK1(JOB)*(TIMEF_START-TIMEF_OBTIME)/(TIMEB*TIMEB)     &
     &  + TIMEF_OBTIME
       ELSEIF(MRAMPFN == 1)THEN
! linear function
       TF(JOB) =                                                        &
     & WRK1(JOB)*(TIMEF_START-TIMEF_OBTIME)/TIMEB + TIMEF_OBTIME
       ENDIF

       S(JOB) =                                                         &
     & WRK1(JOB) * ((CSCALE_START-CSCALE_OBTIME)/TIMEB) + CSCALE_OBTIME
      ENDIF
      ENDDO ! JOB

      DO JOB = FIRST_OBS,LAST_OBS
!     calculate ramp functions between ob time and end of insertion
      IF (WRK1(JOB) <= 0.0 .AND. WRK1(JOB) >= -TIMEA) THEN
       IF(MRAMPFN == 2)THEN
! quadratic function
       TF(JOB) =                                                        &
     & WRK1(JOB)*WRK1(JOB)*(TIMEF_END-TIMEF_OBTIME)/(TIMEA*TIMEA)       &
     &  + TIMEF_OBTIME
       ELSEIF(MRAMPFN == 1)THEN
! linear function
       TF(JOB) =                                                        &
     & WRK1(JOB) * (TIMEF_OBTIME-TIMEF_END)/TIMEA + TIMEF_OBTIME
       ENDIF

       S(JOB) =                                                         &
     & WRK1(JOB) * ((CSCALE_OBTIME-CSCALE_END)/TIMEA) + CSCALE_OBTIME
      ENDIF
      ENDDO ! JOB

!     WRK1 can be reused.

!     IWK1,WRK1 CAN BE REUSED.

      IF (MOD(MHCORFN/16,2) == 1) THEN
        DO JOB = FIRST_OBS,LAST_OBS
!         reduce time factor where scale is large so that the ramp
!         function specified holds for the area integrated weight
          WRK1(JOB) = CSCALE_OBTIME/S(JOB)
          WRK1(JOB) = WRK1(JOB)*WRK1(JOB)
          TF(JOB)   = TF(JOB)*WRK1(JOB)
        ENDDO ! JOB
      ENDIF

      ENDIF   !  Jump here if no obs for this type.

      FIRST_OBS = FIRST_OBS + LENACT(JOBT)
      ENDDO    !  End of JOBT loop over obs types in group.

      IF (LDIAGAC) THEN
       IF (LLDAC(1) .AND. NDAC >  0 .AND. IPASS == 1)THEN
       ENDIF
      ENDIF


      IF(IPASS == 1)THEN ! do stages 1 & 2
!L***                  <<<<    STAGE 1: SCALE    >>>>
!L*** 1.  calculate scale as smoothed average of that of local obs
!         using UM Doc 30 (6.4.17)
!     1.1 Initialize MODEL_INCR
      DO      JC=1,2
      DO      JLEV=1,NO_ANAL_LEVS
      DO      J=1,LENMG
        MODEL_INCR(J,JLEV,JC)=0.
      ENDDO ! J
      ENDDO ! JLEV
      ENDDO ! JC

!     1.2 put values for each ob,level into S_INCR
      DO      JLEV=1,NO_ANAL_LEVS
      DO      JOB=1,LENOB
!       see remark in UM Doc P30 5.1 about the multiple uses of NORMF
!       NORMF here is Q1m of (5.1.2) i.e. 1/E**2
!       TF is used as in (6.1.4)
        S_INCR(JOB,JLEV,1)=NORMF(JOB,JLEV)*TF(JOB)
        S_INCR(JOB,JLEV,2)=NORMF(JOB,JLEV)*TF(JOB)*                     &
     &    S(JOB)*FI_SCALE_FACTOR(FIRST_LEVEL-1+JLEV) ! in metres
!       FI allows a different scale for each level
!       To get results comparable with HORINF use FI_SCALE_FACTOR=1.
      ENDDO ! JOB
      ENDDO ! JLEV

      IF (FILT) THEN
!     1.3 add S_INCR into MODEL_INCR using transpose of model->obs intrp

      DO     JOB=1,LENOB
!       the 4 passes through the loop over JC can be concurrent as
!       long as NPPT(J,1) NPPT(J,2) NPPT(J,3) & NPPT(J,4) are different.
!       The version of HINTCF introduced with FI ensures that they are.
        DO      JC=1,4
!       DO      J=1,2               ! the loops on J & JLEV should
!       DO      JLEV=1,NO_ANAL_LEVS !collapse into a single vector loop
                J=1                     ! collapse loops by hand

        DO      JLEV=1,NO_ANAL_LEVS*2   ! collapse loops by hand
         MODEL_INCR(NPPT(JOB,JC),JLEV,J)=                               &
     &   MODEL_INCR(NPPT(JOB,JC),JLEV,J)+S_INCR(JOB,JLEV,J)*CFPT(JOB,JC)
        ENDDO ! JLEV
!       ENDDO ! J                       ! collapse loops by hand
        ENDDO ! JC
      END DO

!     1.4 divide by grid-box areas, to convert from transpose to adjoint
!         (The grid box area is Px in (6.4.15) ).
!         multiply by area integral of correlation function
      ZINTEGRAL=ZIOA*(FI_SCALE/Earth_Radius)**2 ! convert to radians
      I=0
      DO JROW=1,IROWS
        ZZ=ZINTEGRAL/COS_LAT(JROW)
        DO      JC=1,2
        DO      JLEV=1,NO_ANAL_LEVS
        Z=ZZ*FI_SCALE_FACTOR(FIRST_LEVEL-1+JLEV)**2
        DO      J=1,ROW_LENGTH
          MODEL_INCR(I+J,JLEV,JC)=MODEL_INCR(I+J,JLEV,JC)*Z
        ENDDO ! J
        ENDDO ! JLEV
        ENDDO ! JC
        I=I+ROW_LENGTH
      ENDDO ! JROW

!     1.5 filter weighted scales, and weights
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFCSloop',3)
! Allocate arrays to store row/col deltas.  This is needed since rfcsl
! routine was updated for variable resolution while this routine is not going
! to be (see comment in header)

      ALLOCATE(delta_lat(irows))
      ALLOCATE(delta_long(row_length))
      delta_lat(:) = dlatmg
      delta_long(:)= dlongmg

!     collapse JC & JLEV loops so that multitasking can be over both
!     DO      JC=1,2
!     DO      JLEV=1,NO_ANAL_LEVS
      DO      J=   0,NO_ANAL_LEVS*2-1 ! J loop replaces JC & JLEV loops
      JC=MOD(J,2)+1                   ! J loop replaces JC & JLEV loops
      JLEV=J/2+1                      ! J loop replaces JC & JLEV loops
      Z=FI_SCALE*FI_SCALE_FACTOR(FIRST_LEVEL-1+JLEV)/Earth_Radius
      
!       These calls can be done concurrently
        CALL RFCSL(MODEL_INCR(1,JLEV,JC),IROWS,ROW_LENGTH,M_GRID,       &
     &   Z,COS_LAT,DELTA_LAT,DELTA_LONG,Z,NPASS_RF)
!     ENDDO ! JLEV
!     ENDDO ! JC
      ENDDO ! J                       ! J loop replaces JC & JLEV loops
      DEALLOCATE(delta_lat)
      DEALLOCATE(delta_long)

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFCSloop',4)

!     1.6 divide to give SCALE, using (6.4.15). metre=>radians using A.
      DO      JLEV=1,NO_ANAL_LEVS
      DO      J=1,LENMG
        MODEL_INCR(J,JLEV,IP_SCALE)= (MODEL_INCR(J,JLEV,2) +            &
     &  0.1*FI_SCALE*FI_SCALE_FACTOR(FIRST_LEVEL-1+JLEV))/              &
     &     ((MODEL_INCR(J,JLEV,1)+0.1)*Earth_Radius*FI_VAR_FACTOR)
      ENDDO ! J
      ENDDO ! JLEV


      ENDIF  !  FILT

!L***                  <<<<    STAGE 2: NORMF    >>>>
!L*** 2.  calculate data-density and hence Q to be used in VERTANL
!         using UM Doc 30 (6.4.15)
!     2.1 Initialize MODEL_INCR
      DO      JLEV=1,NO_WT_LEVS
      DO      J=1,LENMG
        MODEL_INCR(J,JLEV,1)=0.
      ENDDO ! J
      ENDDO ! JLEV

!     2.2 put values for each ob,level into S_INCR
!     This need not be done; the values are already there from step 1.2
!     DO JLEV=1,NO_WT_LEVS
!     DO JOB=1,LENOB
!       S_INCR(JOB,JLEV,1)=NORMF(JOB,JLEV)*TF(JOB)
!     ENDDO ! JOB
!     ENDDO ! JLEV

!     2.3 add S_INCR into MODEL_INCR using transpose of model->obs intrp

      DO  JOB=1,LENOB
!       the 4 passes through the loop over JC could be concurrent as
!       long as NPPT(J,1) NPPT(J,2) NPPT(J,3) & NPPT(J,4) are different.
!       The HINTCF introducted with FI makes sure that they differ.
        DO      JC=1,4
        DO      JLEV=1,NO_WT_LEVS
         MODEL_INCR(NPPT(JOB,JC),JLEV,1)=                               &
     &   MODEL_INCR(NPPT(JOB,JC),JLEV,1)+S_INCR(JOB,JLEV,1)*CFPT(JOB,JC)
        ENDDO ! JLEV
        ENDDO ! JC
      END DO

      IF (FILT) THEN
!     2.4 divide by grid-box areas, to convert from transpose to adjoint
!         (The grid box area is Px in (6.4.15) ).
!         multiply by area integral of correlation function
!         Do not add weight given to background; this is done in VERTANL
      I=0
      DO JROW=1,IROWS
        Z=ZIOA/COS_LAT(JROW)
        DO      JLEV=1,NO_WT_LEVS
        DO      J=1,ROW_LENGTH
          MODEL_INCR(I+J,JLEV,1)=MODEL_INCR(I+J,JLEV,1)*Z*              &
     &                           MODEL_INCR(I+J,JLEV,IP_SCALE)**2
        ENDDO ! J
        ENDDO ! JLEV
        I=I+ROW_LENGTH
      ENDDO ! JROW

!     2.5 filter weights to get effective data density
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFVSloop',3)
      DO      JLEV=1,NO_WT_LEVS
!       These calls can be done concurrently
        CALL RFVSL(MODEL_INCR(1,JLEV,1),IROWS,ROW_LENGTH,M_GRID,0.,     &
     &   COS_LAT,DLATMG,DLONGMG,MODEL_INCR(1,JLEV,IP_SCALE),NPASS_RF)
      ENDDO ! JLEV
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFVSloop',4)

      ENDIF  !  FILT

      ELSE ! (IPASS == 2)THEN ! do stage  3

!L***                  <<<<    STAGE 3: INCR     >>>>
!L*** 3.  calculate INCR on model grid
!         using UM Doc 30 (6.4.16)
!     3.1 Initialize MODEL_INCR
      DO      JC=1,NO_ANAL_VAR
      DO      JLEV=1,NO_ANAL_LEVS
      DO      J=1,LENMG
        MODEL_INCR(J,JLEV,JC)=0.
      ENDDO ! J
      ENDDO ! JLEV
      ENDDO ! JC

!     3.2 put values for each ob,level into S_INCR

      DO   JC=1,NO_ANAL_VAR
      DO      JLEV=1,NO_ANAL_LEVS
      DO      JOB=1,LENOB
!       NORMF (from VERTANL) is now Q2m of (6.1.3) i.e. (1/E**2) * Q
!       TF is used as in (6.1.1)
        S_INCR(JOB,JLEV,JC)=OBS_INCR(JOB,JLEV,JC)*TF(JOB)**2*           &
     &                         NORMF(JOB,JLEV)
      ENDDO ! JLEV
      ENDDO ! JOB
      END DO ! JOB


!     3.3 add S_INCR into MODEL_INCR using transpose of model->obs intrp
      DO JOB=1,LENOB
!       the 4 passes through the loop over JC could be concurrent as
!       long as NPPT(J,1) NPPT(J,2) NPPT(J,3) & NPPT(J,4) are different.
        DO      JC=1,4
!       DO      J=1,NO_ANAL_VAR     ! the loops on J & JLEV should
!       DO      JLEV=1,NO_ANAL_LEVS !collapse into a single vector loop
!       the JLEV*J loop is too short to multitask
                J=1                              ! collapse explicitly
        DO      JLEV=1,NO_ANAL_LEVS*NO_ANAL_VAR  ! collapse explicitly
         MODEL_INCR(NPPT(JOB,JC),JLEV,J)=                               &
     &   MODEL_INCR(NPPT(JOB,JC),JLEV,J)+S_INCR(JOB,JLEV,J)*CFPT(JOB,JC)
        ENDDO ! JLEV  & collapsed J
!       ENDDO ! J                                ! collapse explicitly
        ENDDO ! JC
      END DO


      IF (FILT) THEN
!     3.4 divide by grid-box areas, to convert from transpose to adjoint
!         (The grid box area is Px in (6.4.15) ).
!         multiply by area integral of correlation function
!             which is about 16*SCALE**2 for a 2-pass filter.
      I=0
      DO JROW=1,IROWS
        Z=ZIOA/COS_LAT(JROW)
        DO      JC=1,NO_ANAL_VAR ! JC is never equal to IP_SCALE
        DO      JLEV=1,NO_ANAL_LEVS
        DO      J=1,ROW_LENGTH
          MODEL_INCR(I+J,JLEV,JC)=MODEL_INCR(I+J,JLEV,JC)*Z*            &
     &                            MODEL_INCR(I+J,JLEV,IP_SCALE)**2
        ENDDO ! J
        ENDDO ! JLEV
        ENDDO ! JC
        I=I+ROW_LENGTH
      ENDDO ! JROW
      IF (LDIAGAC) THEN
       IF (LLDAC(1).AND.NDAC >  0)THEN
        DO JLEV=1,NO_ANAL_LEVS
         DO JOB=1,LENOB ! interpolate incr at ob position
          WRK1(JOB)=MODEL_INCR(NPPT(JOB,1 ),JLEV,1)*CFPT(JOB, 1)        &
     &             +MODEL_INCR(NPPT(JOB,2 ),JLEV,1)*CFPT(JOB, 2)        &
     &             +MODEL_INCR(NPPT(JOB,3 ),JLEV,1)*CFPT(JOB, 3)        &
     &             +MODEL_INCR(NPPT(JOB,4 ),JLEV,1)*CFPT(JOB, 4)
         ENDDO ! JOB
         IF(ILUV)THEN
         DO JOB=1,LENOB ! interpolate V incr at ob position
          WRK1(JOB)=MODEL_INCR(NPPT(JOB,1 ),JLEV,2)*CFPT(JOB, 1)        &
     &             +MODEL_INCR(NPPT(JOB,2 ),JLEV,2)*CFPT(JOB, 2)        &
     &             +MODEL_INCR(NPPT(JOB,3 ),JLEV,2)*CFPT(JOB, 3)        &
     &             +MODEL_INCR(NPPT(JOB,4 ),JLEV,2)*CFPT(JOB, 4)
         ENDDO ! JOB
         ENDIF ! ILUV
        ENDDO ! JLEV
       ENDIF
      ENDIF

!     3.5 filter increments
      IF(ILUV)THEN  ! filter winds as a vector field
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFVVloop',3)
       DO JLEV=1,NO_ANAL_LEVS
!       These calls can be done concurrently
        IF (model_domain == mt_global) THEN
        CALL RFVVG(MODEL_INCR(1,JLEV,1),MODEL_INCR(1,JLEV,2),           &
     &   IROWS,ROW_LENGTH,M_GRID,COS_LAT,ROW1MGUV,                      &
     &   DLATMG,DLONGMG,MODEL_INCR(1,JLEV,IP_SCALE),NPASS_RF)
        ELSE
        CALL RFVVL(MODEL_INCR(1,JLEV,1),MODEL_INCR(1,JLEV,2),           &
     &   IROWS,ROW_LENGTH,M_GRID,COS_LAT,ROW1MGUV,                      &
     &   DLATMG,DLONGMG,MODEL_INCR(1,JLEV,IP_SCALE),NPASS_RF)
        END IF
       ENDDO ! JLEV
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFVVloop',4)
      ELSE          ! filter scalar increment field
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFVSloop',3)
       DO JLEV=1,NO_ANAL_LEVS
!       These calls can be done concurrently
        IF (model_domain == mt_global) THEN
        CALL RFVSG(MODEL_INCR(1,JLEV,1),IROWS,ROW_LENGTH,M_GRID,        &
     &   COS_LAT,DLATMG,DLONGMG,MODEL_INCR(1,JLEV,IP_SCALE),NPASS_RF)
        ELSE
        CALL RFVSL(MODEL_INCR(1,JLEV,1),IROWS,ROW_LENGTH,M_GRID,0.,     &
     &   COS_LAT,DLATMG,DLONGMG,MODEL_INCR(1,JLEV,IP_SCALE),NPASS_RF)
        END IF
       ENDDO ! JLEV
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('RFVSloop',4)
      ENDIF ! ILUV

      ENDIF ! FILT

      IF (LDIAGAC) THEN
       IF (LLDAC(1).AND.NDAC >  0)THEN
        DO JLEV=1,NO_ANAL_LEVS
         DO JOB=1,LENOB ! interpolate incr at ob position
          WRK1(JOB)=MODEL_INCR(NPPT(JOB,1 ),JLEV,1)*CFPT(JOB, 1)        &
     &             +MODEL_INCR(NPPT(JOB,2 ),JLEV,1)*CFPT(JOB, 2)        &
     &             +MODEL_INCR(NPPT(JOB,3 ),JLEV,1)*CFPT(JOB, 3)        &
     &             +MODEL_INCR(NPPT(JOB,4 ),JLEV,1)*CFPT(JOB, 4)
         ENDDO ! JOB
         IF(ILUV)THEN
         DO JOB=1,LENOB ! interpolate V incr at ob position
          WRK1(JOB)=MODEL_INCR(NPPT(JOB,1 ),JLEV,2)*CFPT(JOB, 1)        &
     &             +MODEL_INCR(NPPT(JOB,2 ),JLEV,2)*CFPT(JOB, 2)        &
     &             +MODEL_INCR(NPPT(JOB,3 ),JLEV,2)*CFPT(JOB, 3)        &
     &             +MODEL_INCR(NPPT(JOB,4 ),JLEV,2)*CFPT(JOB, 4)
         ENDDO ! JOB
         ENDIF ! ILUV
        ENDDO ! JLEV
       ENDIF
      ENDIF


      ENDIF ! IPASS

! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('FI      ',4)
      IF (lhook) CALL dr_hook('FI',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FI
END MODULE fi_mod
