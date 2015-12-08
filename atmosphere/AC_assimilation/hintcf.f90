! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE HINTCF -------------------------------------------------
!LL
!LL  Purpose :
!LL
!LL THIS CALCULATES INTERPOLATION COEFFICIENTS FOR BI-LINEAR
!LL INTERPOLATION OF VALUES FROM MODEL GRID POINTS TO OBSERVATION
!LL POINTS.IT STORES THE COEFFICIENTS IN ANALYSIS WORK ARRAY ANWORK.
!LL IN CALCULATING THE COEFFICIENTS,DISTANCES ARE MEASURED RELATIVE
!LL TO THE NEAREST MODEL POINT TO THE OBSERVATION.INCREMENT VECTORS
!LL (OVER OBSERVATIONS) ARE CALCULATED TO DEFINE THE POSITION
!LL RELATIVE TO THE NEAREST POINT,OF THE FOUR POINTS SURROUNDING
!LL THE OBSERVATION.
!LL
!LL  Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!LL
!LL  Project Task : P3
!LL
!LL WORKS FOR ARAKAWA 'B' GRID ONLY
!LL
!LLEND------------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: AC Assimilation
MODULE hintcf_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE HINTCF (LWIND,LENOB,OBS_LAT,OBS_LONG,OBS_NO,           &
     &                   ROW_LENGTH,P_ROWS,                             &
     &                   CF1PT,CF2PT,CF3PT,CF4PT,                       &
     &                   NP1PT,NP2PT,NP3PT,NP4PT,                       &
     &                   ICODE,lambda_p,phi_p,L_regular,CMESSAGE)
!L    --------------------------------------------------
!
!
      USE conversions_mod, ONLY: pi_over_180, pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx
      IMPLICIT NONE
!--------------------------------------------------------------------
!LCOMDECK COMMG
!L-------------
      REAL          DLAT,DLONG,XLATN,XLONGW
      REAL          ELFPLAT,ELFPLON
      COMMON/COMMG/ DLAT,DLONG,XLATN,XLONGW,ELFPLAT,ELFPLON
!--------------------------------------------------------------------
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
! MPPAC start
      ! dimensions for obs allocated in subroutine AC
      ! dynamically allocated in AC and RDOBS

      ! dimension inxdim allocated in subroutine HORINF
      INTEGER,PARAMETER:: inxdim    = 15000

      ! common for Statistics Calcs in DIAGO ; Prints in RDOBS,GETOBS
      REAL :: R_STAT(MODEL_LEVELS_MAX,0:8)
      REAL :: S_STAT(MODEL_LEVELS_MAX,0:8)
      INTEGER :: COUNTA(NOBTYPMX)
      INTEGER :: COUNTB(NOBTYPMX)
      INTEGER :: COUNTC(NOBTYPMX)

      COMMON /mpp_ac/ R_STAT,S_STAT,COUNTA,COUNTB,COUNTC

      ! common to pass longitudes and latitudes for edges of local
      ! box from setcona to RDOBS and HINTCF

      REAL :: LONG_E
      REAL :: LONG_W
      REAL :: LAT_N,LAT_S
      REAL :: LONG_W_MODEL
      REAL :: LONG_E_MODEL
      COMMON/latlonmax/                                                 &
     &  LONG_E,LONG_W,LAT_N,LAT_S,LONG_W_MODEL,LONG_E_MODEL
! MPPAC end
      EXTERNAL TIMER
!     -----------------------------------------------------------
      LOGICAL LWIND              !IN switch to identify wind grid
      INTEGER LENOB              !IN number of obs
      INTEGER ROW_LENGTH,P_ROWS  !IN size of model grid
!     -----------------------------------------------------------
      REAL     OBS_LAT(LENOB)    !IN ob co-latitudes between 0 & PI
      REAL     OBS_LONG(LENOB)   !IN ob longitudes between 0 & 2*PI
!     NP#PT are pointers to the 4 gridpoints surrounding each ob.
      INTEGER  NP1PT(LENOB)      !OUT pointer to nearest model point
      INTEGER  NP2PT(LENOB)      !OUT pointer to same row, +or-1 pt.
      INTEGER  NP3PT(LENOB)      !OUT pointer to +or-1 row, same pt.
      INTEGER  NP4PT(LENOB)      !OUT pointer to +or-1 row, +or-1 pt.
!     CF#PT are the weights given to each pt pointed to by NP#PT.
      REAL     CF1PT(LENOB)      !OUT interpolation coeffs
      REAL     CF2PT(LENOB)      !OUT interpolation coeffs
      REAL     CF3PT(LENOB)      !OUT interpolation coeffs
      REAL     CF4PT(LENOB)      !OUT interpolation coeffs
      INTEGER  OBS_NO(LENOB)     !IN pointers to obs

      Real lambda_p(1-halo_i:row_length+halo_i)
      Real phi_p(1-halo_i:row_length+halo_i, 1-halo_j:p_rows+halo_j)
      LOGICAL L_regular

!     -----------------------------------------------------------
      INTEGER ICODE              !OUT error code and message
      CHARACTER(LEN=256) CMESSAGE

!     -----------------------------------------------------------
! DECLARATIONS FOR VARIABLE GRID
      REAL GRID_LATS_THISPE(0:P_ROWS)
      REAL GRID_LONS_THISPE(0:ROW_LENGTH)
      REAL R_DELTA_LATS_THISPE(0:P_ROWS)
      REAL R_DELTA_LONS_THISPE(0:ROW_LENGTH)
      REAL DELTA_LATS_THISPE(0:P_ROWS)
      REAL DELTA_LONS_THISPE(0:ROW_LENGTH)
      LOGICAL VARIABLE
      INTEGER I,K
      INTEGER LOCPOS(1)
      INTEGER N_ROW_TEST(LENOB) ! nearest row to ob     minus 1
      INTEGER N_PNT_TEST(LENOB) ! nearest point to ob   minus 1
      INTEGER IOBS_TEMP(LENOB)
      REAL OBS_LONG_TEST,OBS_LAT_TEST
!     -----------------------------------------------------------
!**   DYNAMIC ALLOCATION WITH LENOB
      REAL           WKLON(LENOB)
      REAL           WORK5(LENOB)
!     CF#PT and NP1PT are used for workspace during the calculation
      INTEGER N_ROW(LENOB) ! nearest row to ob     minus 1
      INTEGER N_PNT(LENOB) ! nearest point to ob   minus 1
      INTEGER I_ROW(LENOB) ! increment to row other side of ob   (+or-1)
      INTEGER I_PNT(LENOB) ! increment to point other side of ob (+or-1)
!     -----------------------------------------------------------
      REAL                                                              &
     &               ZDLAT,ZDLONG,ZLATN,ZLONGW,R_ZDLAT,R_ZDLONG,ZLONGE
      INTEGER JOB
!     -----------------------------------------------------------
      REAL PI2P
      PARAMETER (PI2P = 2.0*PI)
!     -----------------------------------------------------------
      REAL Tiny

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      PARAMETER (Tiny = 1.6E-7)
      REAL Tiny1
      PARAMETER (Tiny1 = 9.9E-13)
!     -----------------------------------------------------------

      IF (lhook) CALL dr_hook('HINTCF',zhook_in,zhook_handle)
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('HINTCF  ',3)

!
! PART 1
! ------
! SET UP LATITUDE OF NORTHERN ROW AND LONGITUDE OF WESTERN POINT
! FOR EACH MODEL GRID
!
! PARAMETERS FOR MODEL GRID
!
!     GET MODEL GRID INFORMATION FROM COMMG AND CONVERT TO RADIANS
!     CONVERT XLATN TO CO-LATITUDE FIRST
!
      ZDLAT  = DLAT         * PI_OVER_180
      R_ZDLAT = 1.0/ZDLAT
      ZDLONG = DLONG        * PI_OVER_180
      R_ZDLONG = 1.0/ZDLONG
      ZLATN  = (90.0-LAT_N) * PI_OVER_180
      ZLONGW = LONG_W_MODEL * PI_OVER_180
      ZLONGE = LONG_E_MODEL * PI_OVER_180

!     MAKE SURE ZLONGW,ZLONGE ARE IN RANGE 0-->2*PI
      IF (ZLONGW <  0.0)  THEN
        ZLONGW=ZLONGW+PI2P
      ELSEIF (ZLONGW >  PI2P) THEN
        ZLONGW=ZLONGW-PI2P
      ENDIF
      IF (ZLONGE <  0.0)  THEN
        ZLONGE=ZLONGE+PI2P
      ELSEIF (ZLONGE >  PI2P) THEN
        ZLONGE=ZLONGE-PI2P
      ENDIF

! Set up cooridinates of grid (in radians)
! If a regular grid these are calculated from input parameters otherwise
! copy from arrays phi_p and lambda_p
! Also note that variable grids have monotonically increasing longitudes

      IF(L_regular) THEN
        DO I=1,P_ROWS
          GRID_LATS_THISPE(I)=(90.0-LAT_N)+REAL(I-1)*DLAT
          GRID_LATS_THISPE(I)=GRID_LATS_THISPE(I) * PI_OVER_180

! spacing and reciprocal between row
          DELTA_LATS_THISPE(I)= DLAT * PI_OVER_180
          R_DELTA_LATS_THISPE(I)=1.0/DELTA_LATS_THISPE(I)
        ENDDO

        DO I=1,ROW_LENGTH
          GRID_LONS_THISPE(I)=LONG_W_MODEL+(I-1)*DLONG
! allow for Meridian
          IF(GRID_LONS_THISPE(I).GT.360.0) THEN
            GRID_LONS_THISPE(I)=GRID_LONS_THISPE(I)-360.0
          ENDIF
          GRID_LONS_THISPE(I)=GRID_LONS_THISPE(I) * PI_OVER_180

! spacing and reciprocal between cols
          DELTA_LONS_THISPE(I)=DLONG * PI_OVER_180
          R_DELTA_LONS_THISPE(I)=1.0/DELTA_LONS_THISPE(I)
        ENDDO

      ELSE ! variable grid
         
        DO I=1,ROW_LENGTH
          GRID_LONS_THISPE(I)=lambda_p(I)
        ENDDO
        DO I=1,P_ROWS
! colatitudes in radians are required
          GRID_LATS_THISPE(I)=phi_p(1,p_rows-I+1)/PI_OVER_180
          GRID_LATS_THISPE(I)=90.0-GRID_LATS_THISPE(I)
          GRID_LATS_THISPE(I)=GRID_LATS_THISPE(I)*PI_OVER_180
        ENDDO

! Calculate reciprocal latitude spacings between rows
        DO I=1,P_ROWS-1
          DELTA_LATS_THISPE(I)=GRID_LATS_THISPE(I+1)-GRID_LATS_THISPE(I)
          R_DELTA_LATS_THISPE(I)=1.0/DELTA_LATS_THISPE(I)
        ENDDO
        DELTA_LATS_THISPE(P_ROWS)=DELTA_LATS_THISPE(P_ROWS-1)
        R_DELTA_LATS_THISPE(P_ROWS)=R_DELTA_LATS_THISPE(P_ROWS-1)

! Calculate reciprocal of longitude spacings between columns
        DO I=1,ROW_LENGTH-1
          DELTA_LONS_THISPE(I)=GRID_LONS_THISPE(I+1)-GRID_LONS_THISPE(I)
          R_DELTA_LONS_THISPE(I)=1.0/DELTA_LONS_THISPE(I)
        ENDDO
        DELTA_LONS_THISPE(ROW_LENGTH)=DELTA_LONS_THISPE(ROW_LENGTH-1)
        R_DELTA_LONS_THISPE(ROW_LENGTH)=R_DELTA_LONS_THISPE(ROW_LENGTH-1)

      ENDIF ! L_regular

! Calculate zeroth row/col details
      DELTA_LATS_THISPE(0)=DELTA_LATS_THISPE(1)
      DELTA_LONS_THISPE(0)=DELTA_LONS_THISPE(1)
      R_DELTA_LATS_THISPE(0)=R_DELTA_LATS_THISPE(1)
      R_DELTA_LONS_THISPE(0)=R_DELTA_LONS_THISPE(1)
      GRID_LATS_THISPE(0)=GRID_LATS_THISPE(1)-DELTA_LATS_THISPE(0)
      GRID_LONS_THISPE(0)=GRID_LONS_THISPE(1)-DELTA_LONS_THISPE(0)
!
!
! PART 2
! ------
! CALCULATE (NEAREST-1) MODEL GRID PT TO OBSERVATION
!
!L    INPUT  VECTOR IN OBS_LAT  - COLATITUDE (RADIANS) SET UP IN AC
!L                  IN OBS_LONG - LONGITUDE  (RADIANS) SET UP IN AC
!
!L    OUTPUT VECTORS
!
!L    OUTPUT VECTORS OF INTERPOLATION COEFFICIENTS ARE IN
!L    CF1PT CF2PT CF3PT CF4PT NP1PT NP2PT NP3PT NP4PT

! The TINY1 offset is required to prevent invalid 
! NP1PT NP2PT NP3PT NP4PT on variable grids (-ve values)
!
! 2.1
! ---
!
      DO JOB=1,LENOB
  
!       FIND NEAREST MODEL ROW minus 1.  in N_ROW
        IF(L_regular) THEN
          N_ROW(JOB) = NINT((OBS_LAT(JOB)-ZLATN)*R_ZDLAT)
        ELSE
          OBS_LAT_TEST=OBS_LAT(JOB)-TINY1

          N_ROW(JOB)=P_ROWS-1
          IF(OBS_LAT_TEST.LT.GRID_LATS_THISPE(1)) THEN
            N_ROW(JOB)=0
          ELSE
           
            DO K=1,P_ROWS-1
              IF(OBS_LAT_TEST .GE. GRID_LATS_THISPE(K).AND.                 &
     &           OBS_LAT_TEST .LT. GRID_LATS_THISPE(K+1)) THEN
                OBS_LAT_TEST=GRID_LATS_THISPE(K+1)-OBS_LAT_TEST
! To mimic regular grid check which half of grid box point lies in
                IF(OBS_LAT_TEST.LT.DELTA_LATS_THISPE(K)/2.0) THEN
                  N_ROW(JOB)=K
                ELSE
                  N_ROW(JOB)=K-1
                ENDIF
                EXIT
             ENDIF
           ENDDO
          ENDIF
        ENDIF
!
! 2.2
! ---
! FIND NEAREST POINT minus 1 ON ROW.   in N_PNT
!
!
!       SPECIAL TREATMENT WHEN LIMITED AREA SPANS MERIDIAN.
!       OBS EAST OF THE MERIDIAN HAVE LONGITUDES < WESTERN BOUNDARY.
!       CALCULATE THEIR LONGITUDE DISPLACEMENT WITH EXTRA 2*PI SHIFT.
!
!       INITIALISE WORK ARRAY OF LONGITUDE SHIFTS
        WKLON(JOB)=0.0
        IF(L_regular) THEN
          IF (ZLONGW  >   ZLONGE) THEN
!         ASSUME AREA STRADDLES MERIDIAN
!         SET UP LONGITUDE SHIFT FOR OBS BETWEEN MERIDIAN
!                                      AND EASTERN BOUNDARY.
            IF (OBS_LONG(JOB) >= 0.0.AND.OBS_LONG(JOB)                  &
     &                    <= (ZLONGE+ZDLONG/2.0))   WKLON(JOB)= PI2P
        ENDIF
          N_PNT(JOB) = NINT((OBS_LONG(JOB)+WKLON(JOB)-ZLONGW)*R_ZDLONG)
        ELSE
          OBS_LONG_TEST=OBS_LONG(JOB)-TINY1

! Check to see if any part of this PE crosses Meridan
! (Variable grids work in monotonically increasing longitude)
          IF(MAXVAL(GRID_LONS_THISPE) > PI2P) THEN  
! Observations east of the Meridian need 2PI added to longitude
            IF((OBS_LONG_TEST-GRID_LONS_THISPE(0)) <= 0.0) THEN
              WKLON(JOB)= PI2P
              OBS_LONG_TEST=OBS_LONG_TEST+WKLON(JOB)
            ENDIF
          ENDIF
    
          N_PNT(JOB)=ROW_LENGTH-1
          IF(OBS_LONG_TEST.LT.GRID_LONS_THISPE(1)) THEN
            N_PNT(JOB)=0
          ELSE
           
            DO K=1,ROW_LENGTH-1
              IF(OBS_LONG_TEST .GE. GRID_LONS_THISPE(K).AND.                 &
     &           OBS_LONG_TEST .LT. GRID_LONS_THISPE(K+1)) THEN
                OBS_LONG_TEST=GRID_LONS_THISPE(K+1)-OBS_LONG_TEST
! To mimic regular grid check which half of grid box point lies in
                IF(OBS_LONG_TEST .LT. DELTA_LONS_THISPE(K)/2.0) THEN
                  N_PNT(JOB)=K
                ELSE
                  N_PNT(JOB)=K-1
                ENDIF
                EXIT
             ENDIF
           ENDDO
          ENDIF
        ENDIF 
!
! PART 3
! ------
! CALCULATE INTERPOLATION COEFFS. In lat in CF1PT, in long in CF2PT.
!
        IF(L_Regular) THEN
          CF1PT(JOB) = OBS_LAT(JOB) -(ZLATN +N_ROW(JOB)*ZDLAT )
          CF2PT(JOB) = OBS_LONG(JOB)-(ZLONGW+N_PNT(JOB)*ZDLONG)         &
     &                 + WKLON(JOB)
       ELSE
         IF(N_ROW(JOB).LT.P_ROWS) THEN
            CF1PT(JOB)=OBS_LAT(JOB)-GRID_LATS_THISPE(N_ROW(JOB)+1)
         ELSE
            CF1PT(JOB)=1.0  ! can this happen?
         ENDIF
      
 
          IF(N_PNT(JOB).LT.ROW_LENGTH) THEN
            CF2PT(JOB)=OBS_LONG(JOB)-GRID_LONS_THISPE(N_PNT(JOB)+1)     &
     &                + WKLON(JOB)
          ELSE
            CF2PT(JOB)=1.0   ! can this happen?
          ENDIF
        ENDIF
!
! GET INCREMENT VECTORS TO OBTAIN FOUR SURROUNDING GRID POINTS
! INCLUDE POSSIBILITY OF OBS BEING ON SOUTHERN OR EASTERN HALO
! BOUNDARIES. (Cannot be on northern or western because of the way
! obs are distributed in RDOBS)
!
!   N.B Value of Tiny is chosen so that observations within
!   approx 1 metre of grid edge are affected.

!
        IF (CF1PT(JOB) >= 0.0) THEN
          I_ROW(JOB)=1
!         Set I_ROW to 0 if CF1PT is very small.
          IF (CF1PT(JOB)  <   Tiny) THEN
            I_ROW(JOB) = 0
          END IF
        ELSE
          I_ROW(JOB)=-1
        ENDIF
!
        IF (CF2PT(JOB) >= 0.0) THEN
          I_PNT(JOB) = 1
!         Set I_PNT to 0 if CF2PT is very small.
          IF (CF2PT(JOB)  <   Tiny) THEN
            I_PNT(JOB) = 0
          END IF
        ELSE
          I_PNT(JOB) = -1
        ENDIF
!
        CF1PT(JOB) = ABS(CF1PT(JOB))
        CF2PT(JOB) = ABS(CF2PT(JOB))

!
!       Make sure that all pointers are within grid, and different
        NP1PT(JOB)=N_ROW(JOB)+I_ROW(JOB) ! use NP1PT as workspace
        IF(NP1PT(JOB) <  0)THEN          ! row-1 is outside grid
          I_ROW(JOB)=1                   ! point instead to row+1
          CF1PT(JOB)=0.                  ! but do not use row+1.
        ENDIF

        IF(NP1PT(JOB) >= P_ROWS)THEN     ! row+1 is outside grid
          I_ROW(JOB)=-1                  ! point instead to row-1
          CF1PT(JOB)=0.                  ! but do not use row-1.
        ENDIF
        IF ((N_PNT(JOB) + I_PNT(JOB)) <  1) THEN
          I_PNT(JOB) = 1
          CF2PT(JOB) = 0.
        ENDIF

        IF ((N_PNT(JOB) + I_PNT(JOB)) >=  ROW_LENGTH) THEN
          I_PNT(JOB) = -1
          CF2PT(JOB) = 0.
        ENDIF

        N_ROW(JOB)=MAX(N_ROW(JOB),0)
        N_PNT(JOB)=MAX(N_PNT(JOB),0)
!
!       Convert to 2-D coeffs for 4 surrounding gridpoints
        IF(L_regular) THEN
           CF4PT(JOB)  = CF1PT(JOB)*R_ZDLAT   ! lat coeff for row +or-1
           WORK5 (JOB) = CF2PT(JOB)*R_ZDLONG  ! long coeff for pt +or-1
        ELSE
          CF4PT(JOB) = CF1PT(JOB)*R_DELTA_LATS_THISPE(N_ROW(JOB))
          WORK5(JOB) = CF2PT(JOB)*R_DELTA_LONS_THISPE(N_ROW(JOB))
        ENDIF
        CF2PT(JOB)  = 1.0-CF4PT(JOB)       ! lat  coeff for row
        CF3PT(JOB)  = 1.0-WORK5(JOB)       ! long coeff for pt
!
        CF1PT(JOB) = CF2PT(JOB)*CF3PT(JOB) ! row, pt
        CF2PT(JOB) = CF2PT(JOB)*WORK5(JOB) ! row, pt +or-1
        CF3PT(JOB) = CF4PT(JOB)*CF3PT(JOB) ! row +or-1, pt
        CF4PT(JOB) = CF4PT(JOB)*WORK5(JOB) ! row +or-1, pt +or-1
!
!  up to now have assumed north most row is first
!  but ND order is reversed, so adjust use of N_ROW
        NP1PT(JOB) =N_PNT(JOB) + 1 +                                    &
     &                       (P_ROWS -1 - N_ROW(JOB) )*ROW_LENGTH
        NP2PT(JOB) =NP1PT(JOB) + I_PNT(JOB)

!  respecting S->N order of ND fields,
!  use -I_ROW to point to other row adjacent to ob
        NP3PT(JOB) =NP1PT(JOB) - I_ROW(JOB)*ROW_LENGTH
        NP4PT(JOB) =NP3PT(JOB) + I_PNT(JOB)

      if(NP1PT(JOB) >  ROW_LENGTH*P_ROWS)then
          CMESSAGE = 'NP1PT(JOB) too big!!'
        write(0,*)'NP1PT(JOB) too big!!',NP1PT(JOB),job,mype
          ICODE=1
        endif
        if(NP1PT(JOB) <  1)then
          CMESSAGE = 'NP1PT(JOB) too small!!'
        write(0,*)'NP1PT(JOB) too small!!',NP1PT(JOB),job,mype
          ICODE=1
        endif

      if(NP2PT(JOB) >  ROW_LENGTH*P_ROWS)then
          CMESSAGE = 'NP2PT(JOB) too big!!'
          write(0,*)'NP2PT(JOB) too big!!',NP2PT(JOB),job,mype
          ICODE=1
        endif
        if(NP2PT(JOB) <  1)then
          CMESSAGE = 'NP2PT(JOB) too small!!'
          write(0,*)'NP2PT(JOB) too small!!',NP2PT(JOB),job,mype
          ICODE=1
        endif

      if(NP3PT(JOB) >  ROW_LENGTH*P_ROWS)then
          CMESSAGE = 'NP3PT(JOB) too big!!'
          write(0,*)'NP3PT(JOB) too big!!',NP3PT(JOB),job,mype
          ICODE=1
        endif
        if(NP3PT(JOB) <  1)then
          CMESSAGE = 'NP3PT(JOB) too small!!'
          write(0,*)'NP3PT(JOB) too small!!',NP3PT(JOB),job,mype
          ICODE=1
        endif

      if(NP4PT(JOB) >  ROW_LENGTH*P_ROWS)then
          CMESSAGE = 'NP4PT(JOB) too big!!'
          write(0,*)'NP4PT(JOB) too big!!',NP4PT(JOB),job,mype
          ICODE=1
        endif
        if(NP4PT(JOB) <  1)then
          CMESSAGE = 'NP4PT(JOB) too small!!'
          write(0,*)'NP4PT(JOB) too small!!',NP4PT(JOB),job,mype
          ICODE=1
        endif

      ENDDO ! JOB
!
!
! DEPENDS ON: timer
      IF (LTIMER_AC) CALL TIMER('HINTCF  ',4)
      IF (lhook) CALL dr_hook('HINTCF',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE HINTCF
END MODULE hintcf_mod
