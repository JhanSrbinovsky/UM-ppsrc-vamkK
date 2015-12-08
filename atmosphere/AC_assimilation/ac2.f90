! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: AC
!
!    Purpose : Main control routine for AC Scheme
!     This makes the analysis correction to model fields at the start
!   of each timestep of model integration.  The corrections are
!   calculated at analysis grid points,then interpolated to the model
!   grid before being added to the model fields. The corrections are
!   weighted sums of increments (ob-model) obtained at analysis grid
!   points. The weights depend on the distance between observation and
!   grid point, Observational error and the timeliness of the
!   observation.  Observations are analysed sequentially by observation
!   type in the order given by the array lact in comdeck comacp. A
!   vertical profile of observational increments and errors at
!   observation points are obtained before the weighted sums are
!   calculated. Horizontal filtering of increments on the model grid may
!   be done,and hydrostatic and geostrophic increments are calculated on
!   the model grid if appropriate.mean and rms diagnostics may be
!   calculated.
!
!
!    Programming Standard : UM Doc Paper No 4 ; Version 3 ; 7/9/90
!
!    Project Task : P3
!*
!*L  Arguments and declarations:

!*L  Arguments & declarations:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
MODULE ac2_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE AC2(P_LEVELS, Q_LEVELS, BL_LEVELS,                     &
     &               ROW_LENGTH, P_ROWS,                                &
     &               P_FIELD,                                           &
     &               TIMESTEP_NO, ITER_NO, TIMESTEP,                    &
     &               OBS, TNDV,                                         &
     &               EXNER, PSTAR, THETA, RH, QCL, QCF,                 &
     &               CONV_CLD, LS_RAIN, LS_SNOW, CONV_RAIN, CONV_SNOW,  &
     &               D_THETA_DT_CONV,D_THETA_DT_LS,                     &
     &               LAYER_CLOUD,PRESSURE,                              &
     &               RHCRIT,                                            &
     &               OBS_NO, LENOB, NO_ANAL_LEVS, NO_WT_LEVS,           &
     &               NO_ANAL_VAR,                                       &
     &               LENAG, LENMG, WKLEN, INC_TYPE,                     &
     &               NAVT, JGROUP, LWIND, IACTF, IACTL, DG_ONLY,        &
     &               STINDEX, STLIST, LEN_STLIST, SI, SF,               &
     &               STASHWORK, STASH_LEVELS, NUM_STASH_LEVELS,         &
     &               STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,             &
     &               lambda_p,phi_p,L_regular,                          &
     &               ICODE, CMESSAGE)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE comobs_mod, ONLY: nobtypmx, obs_info
      USE domain_params
      USE ac_stash_mod, ONLY: ac_stash
      USE addinc_mod, ONLY: addinc
      USE diago_mod, ONLY: diago
      USE diagopr_mod, ONLY: diagopr
      USE fi_mod, ONLY: fi
      USE getob2_mod, ONLY: getob2
      USE hintcf_mod, ONLY: hintcf
      USE lhn_inc_mod, ONLY: lhn_inc
      USE mmspt_mod, ONLY: mmspt
      USE relaxc_mod, ONLY: relaxc
      USE set_relax_cf_mod, ONLY: set_relax_cf
      USE vertanl_mod, ONLY: vertanl
      IMPLICIT NONE

! Global parameters:

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

! Import arguments (INTENT=IN):
      INTEGER        P_LEVELS                ! Total number of levels
      INTEGER        Q_LEVELS                ! Number of wet levels
      INTEGER        BL_LEVELS               ! Bdy layer levels
      INTEGER        ROW_LENGTH              ! Number of points on row
      INTEGER        P_ROWS                  ! Number of rows (pstar)
      INTEGER        P_FIELD                 ! Number of points in
                                             ! mass field
      INTEGER        TNDV
      INTEGER        TIMESTEP_NO             ! Current model timestep
      INTEGER        ITER_NO                 ! Iteration number
                                             ! excluding diagnostic
      INTEGER        LENMG                   ! Length of model grid
      INTEGER        LENAG                   ! Length of analysis grid
      INTEGER        WKLEN                   ! Length of 2nd dimension
                                             ! of array for derived incs
                                             ! produced by HTDRST,
                                             ! GEOSTR & WINDBAL.
      INTEGER        INC_TYPE                ! Index for data types in
                                             ! MODEL_INCR
      INTEGER        LENOB                   ! No of obs in group
      INTEGER I,IPROC
      INTEGER        NO_ANAL_LEVS            ! No of analysis levels
                                             ! for group of obs types
      INTEGER        NO_ANAL_VAR             ! No of variables being
                                             ! analysed (2 for winds)
      INTEGER        NO_WT_LEVS              ! No of weight levels
                                             ! for group of obs types
      INTEGER        NAVT                    ! Analysis variable type
                                             ! 1/2/3/4/5 =
                                             ! p*/theta/winds/rh/precip
      INTEGER        JGROUP                  ! Loop counter over groups
      INTEGER        IACTF                   ! First obs type in group
      INTEGER        IACTL                   ! Last  obs type in group
      INTEGER        OBS_NO(LENOB+1)

      LOGICAL        LWIND                   ! Indicator if group
                                             ! has wind obs types
      LOGICAL        DG_ONLY                 ! Indicator if diagnostic
                                             ! only iteration
      REAL           TIMESTEP                ! Timestep in seconds
      REAL           OBS(TNDV+1)             ! Observation data (set up
                                             ! in RDOBS)
      REAL           CONV_CLD(P_FIELD, Q_LEVELS)  ! conv cld amount
      REAL           LS_RAIN(P_FIELD)        ! large scale rain rate
      REAL           LS_SNOW(P_FIELD)        ! large scale snow rate
      REAL           CONV_RAIN(P_FIELD)      ! convective rain rate
      REAL           CONV_SNOW(P_FIELD)      ! convective snow rate
                                             ! above rates diagnostic
      REAL           D_THETA_DT_CONV(P_FIELD,Q_LEVELS)
                                       ! convective latent heating rate
      REAL           D_THETA_DT_LS(P_FIELD,Q_LEVELS)
                                      ! large scale latent heating rate

      REAL           RHCRIT(Q_LEVELS)        ! Critical RH array
                                             ! for cloud
! These are used in variable resolution runs
      REAL lambda_p(1-halo_i:row_length+halo_i)
      REAL phi_p(1-halo_i:row_length+halo_i, 1-halo_j:p_rows+halo_j)
      LOGICAL L_regular

! Import / export arguments (INTENT=INOUT):
      REAL           EXNER(P_FIELD, P_LEVELS) ! exner on theta levels
      REAL           LAYER_CLOUD(P_FIELD, Q_LEVELS) ! as name says
      REAL           PRESSURE (P_FIELD, P_LEVELS) ! p on theta levels
      REAL           PSTAR(P_FIELD)          ! Prognostic variable P*
      REAL           THETA(P_FIELD, P_LEVELS)! Prognostic variable
                                             ! theta
      REAL           RH(P_FIELD, Q_LEVELS)   ! Prognostic variable HMR
! for 2A cloud microphysics, but otherwise contains rh at this stage
      REAL           QCL(P_FIELD, Q_LEVELS)   ! Prognostic variable QCL
      REAL           QCF(P_FIELD, Q_LEVELS)   ! Prognostic variable QCF

! Export arguments (INTENT=OUT):
      INTEGER        ICODE                   ! Error code (0 = OK)

      CHARACTER(LEN=256)  CMESSAGE                ! Error message

! Stash variables
      REAL           STASHWORK(*)

      INTEGER        LEN_STLIST
      INTEGER        NUM_STASH_LEVELS
      INTEGER        NUM_STASH_PSEUDO
      INTEGER        STINDEX(2, *)
      INTEGER        STLIST(LEN_STLIST, *)
      INTEGER        SI(*)
      INTEGER        STASH_LEVELS(NUM_STASH_LEVELS +1, *)
      INTEGER        STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO +1, *)

      LOGICAL        SF(*)

!*
!L Single iteration of analysis correction data assimilation
!L for observations or group of observations of variable type
!L defined by NAVT

! More global variables:
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

! Local (dynamic) arrays:
      INTEGER        MODEL_OBS_TYPE(LENOB+1)
                                         ! Model observation type No
      INTEGER        NE_AG_PT(LENOB+1)     ! Nearest analysis grid
                                            ! point to obs.
      INTEGER        NPPT(LENOB+1, 4)       ! horizontal model ->
                                            ! obs interp pointers

      REAL           OBS_LAT (LENOB+1)      ! Observation latitude
      REAL           OBS_LONG(LENOB+1)      ! "longitude
      REAL           OBS_TIME(LENOB+1)      ! "time
      REAL           OBS_INCR(LENOB+1, NO_ANAL_LEVS, NO_ANAL_VAR)!
                                            ! Obs - Model increments
                                            ! at obs points.
      REAL           NORMF(LENOB+1, NO_ANAL_LEVS)
                                               ! Normalization factor
      REAL           ANAL_INCR_LOCAL(LENMG*NO_ANAL_LEVS*NO_ANAL_VAR) !
                                             ! Accumulated increments
                                             ! on analysis grid
      REAL           MODEL_INCR(LENMG, NO_ANAL_LEVS, INC_TYPE) !
                                             ! Accumulated increments
                                             ! on model grid
      REAL           DRINCS(P_FIELD, WKLEN)  ! Array to hold derived
                                             ! increment fields
      REAL           CFPT(LENOB+1, 4)
                                             ! obs interp coeffs
      REAL           RELAX_CF(P_ROWS)        ! Relaxation coeffs for
                                             ! current group

      LOGICAL        LMISSD(LENOB+1, NO_ANAL_LEVS)
                                             ! missing data in obs.
      REAL           pstgs(P_FIELD)          ! PSTAR at U points
      INTEGER dum3
      REAL    dum2(1,1), dum1(1)
      INTEGER bc,bcount
      INTEGER lev_ave,lev_start,lev_end,num_levs,ntimes,lev_times
      INTEGER lev_rem
      INTEGER num_incr_levs
      INTEGER ic
      INTEGER pe_group,iii,jjj
      INTEGER num_groups,                                               &
     &    g_pe_start(0:3),g_pe_end(0:3),                                &
     &    g_lev_start(0:3),g_lev_end(0:3)
      COMMON/split_pe/num_groups,                                       &
     &     g_pe_start,g_pe_end,g_lev_start,g_lev_end

! Local scalars:
      INTEGER        J                       ! Miscelaneous loop cntr
      INTEGER        JACT,JJACT              ! Loop counter over obs
                                             ! types in group
      INTEGER        JANL                    ! Loop counter over
                                             ! analysis levels
      INTEGER        JLEV,JJLEV              ! Loop counter over model
                                             ! levels
      INTEGER        JL,JJL
      INTEGER        JK,JJK
      INTEGER        JPOINT                  ! Loop counter over points

      INTEGER        IPASS                   ! 1 if weights pass
                                             ! 2 if increment pass
      INTEGER        IAG                     ! Offset for level in
                                             ! ANAL_INCR array.
      INTEGER        IJK                     ! Level number in HYDRST
                                             ! & GEOSTR loops
      INTEGER        INOBS                   ! No of obs for this type
      INTEGER        NPTOBT                  ! Offset to first obs of
                                             ! each type in group
      INTEGER        LENOBT                  ! No of obs for obs type
      INTEGER        NDV                     ! No of data values
      INTEGER        N_ROWS                  ! No of rows
      INTEGER        NO_ANAL_INCR_LEVS       ! No of analysis increment
                                             ! levels.
      INTEGER        MODE_HANAL              ! FI or HORINF for this
                                             ! group
      INTEGER        STASH_ITEM_NO           ! Stash Item Number
      LOGICAL        SURFACE_WIND            ! True if surface wind data
      LOGICAL ILUV

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! Subroutine calls:
      EXTERNAL  TIMER

!- End of header section.


!L Subroutine structure:
!L 1.0 Initialize:
! Mode of Horizontal Analysis for this group
      IF (lhook) CALL dr_hook('AC2',zhook_in,zhook_handle)
      MODE_HANAL = DEF_MODE_HANAL(GROUP_INDEX(JGROUP))

!L 2.0 Weights calculation starts here

      IF (LDIAGAC.AND.mype == 0) THEN
        PRINT '(/,A,(10I4))',                                           &
     &                       ' START OF WTS ANALYSIS PHASE FOR TYPES ', &
     &                                      (LACT(J), J = IACTF, IACTL)

      END IF

!L 2.1 Setup for interpolation model -> obs & vertical analysis
!L     Get common format data for all obs in group       (getob2)
! Vectorisation over all observation types in group
! Loop over obs types in group
      if(lenob /= 0)then
      NPTOBT = 0

      DO JACT = IACTF, IACTL
      JJACT=JACT
        INOBS  = OBS_INFO % NOBS(JACT)
        LENOBT = LENACT(JACT)

        IF (INOBS >  0 .AND. LENOBT >  0) THEN
          CALL GETOB2 (JJACT, OBS(OBS_INFO % MDISPOBT(JACT) +1), INOBS, &
     &                OBS_LAT(NPTOBT +1), OBS_LONG(NPTOBT +1),          &
     &                OBS_TIME(NPTOBT +1), MODEL_OBS_TYPE(NPTOBT +1),   &
     &                OBS_NO(NPTOBT +1), LENOBT, ICODE, CMESSAGE)

         ! Check for error - exit routine if occured
          IF (ICODE  >   0) GOTO 999

        END IF

       ! Move pointer to next obs type in group
        NPTOBT = NPTOBT + LENOBT

      END DO

!L 2.2  Horizontal interpolation coefficients model -> obs  (hintcf)
      CALL HINTCF(LWIND, LENOB, OBS_LAT, OBS_LONG, OBS_NO, ROW_LENGTH,  &
     &           P_ROWS, CFPT(1,1), CFPT(1,2), CFPT(1,3), CFPT(1,4),    &
     &           NPPT(1,1), NPPT(1,2), NPPT(1,3), NPPT(1,4), ICODE,     &
     &           lambda_p,phi_p,L_regular,CMESSAGE)

! Check for error - exit routine if occured
      IF (ICODE  >   0) GOTO 999

!L 2.3 Do vertical analysis
!L     Make increment & error factor vectors for each type in this
!L     group. The details of processing method depend on ob type
!L     so vectorization is over LENOBT obs in one type
!L     rather than LENOB obs in the group. NPTOBT gives the
!L     increment to point to each section in turn of those
!L     vectors which go over all obs in group.

! Loop over obs types in group
      endif ! lenob /= 0
      NPTOBT = 0

      DO JACT = IACTF, IACTL
      JJACT=JACT
        INOBS  = OBS_INFO % NOBS(JACT)
        LENOBT = LENACT(JACT)
!
!
!
      if(inobs == 0.or.lenobt == 0)then
      if(LDIAGAC)then
      if(LACT(JACT) == 506)then
        CALL DIAGOPR (1,dum1,dum2,LENOBT,LENOB,NDV,                     &

     &                LMISSD,dum3,NPTOBT,NO_ANAL_LEVS)
      else
      BCOUNT=2
      if(LACT(JACT) == 101.or.                                          &
     &   (LACT(JACT) >= 202.and.LACT(JACT) <= 204).or.                  &
     &   (LACT(JACT) >= 302.and.LACT(JACT) <= 306).or.                  &
     &   (LACT(JACT) >= 402.and.LACT(JACT) <= 404).or.                  &
     &   LACT(JACT) == 406.or.                                          &
     &   LACT(JACT) == 901)BCOUNT=1
!
! BARRIERS FOR DIAGO
!
      do bc=1,bcount
      do jlev=1,glsize(3,fld_type_p)
      do i=0,8
      s_stat(jlev,i)=0.0
      enddo
      enddo
      if(LACT(JACT) == 406)then
! dummy diago call
       CALL DIAGO ('MULTI-LEVEL', LACT(JACT),6,OBS_INCR,NORMF,          &
     &              OBS_LAT, OBS_LONG, LMISSD, LENOB, LENOBT, NPTOBT,   &
     &              NO_ANAL_LEVS, NO_ANAL_VAR)
      else
! dummy diago call
       CALL DIAGO ('VAN?', LACT(JACT), 3+bc-bcount,                     &
     &              OBS_INCR, NORMF,                                    &
     &              OBS_LAT, OBS_LONG, LMISSD, LENOB, LENOBT, NPTOBT,   &
     &              NO_ANAL_LEVS, NO_ANAL_VAR)
      endif ! LACT(JACT) == 406
      enddo
      endif ! LACT(JACT) == 506
      endif ! LDIAGAC
      else !inobs == 0.or.lenobt == 0

        IF (INOBS  >   0 .AND. LENOBT  >   0) THEN
          NDV   = OBS_INFO % NDATAV(JACT) - OBS_INFO % NDVHDR
          IPASS = 1

          CALL VERTANL (PSTGS,                                          &
     &          JJACT, IPASS, LENOBT, NDV, OBS_NO(NPTOBT +1),           &
     &          MODEL_OBS_TYPE(NPTOBT +1), OBS_LAT(NPTOBT +1),          &
     &          OBS_LONG(NPTOBT +1), OBS(OBS_INFO % MDISPOBT(JACT) +1), &
     &          INOBS, EXNER, PSTAR, THETA, RH, QCL, QCF, CONV_CLD,     &
     &          LAYER_CLOUD,PRESSURE,                                   &
     &          LS_RAIN, LS_SNOW, CONV_RAIN, CONV_SNOW,                 &
     &          RHCRIT,                                                 &
     &          P_FIELD, MODEL_INCR, LENMG, NO_WT_LEVS,                 &
     &          CFPT(NPTOBT+1, 1), CFPT(NPTOBT+1, 2),                   &
     &          CFPT(NPTOBT+1, 3), CFPT(NPTOBT+1, 4),                   &
     &          NPPT(NPTOBT+1, 1), NPPT(NPTOBT+1, 2),                   &
     &          NPPT(NPTOBT+1, 3), NPPT(NPTOBT+1, 4),                   &
     &          OBS_INCR, NORMF, LMISSD,                                &
     &          P_LEVELS, Q_LEVELS, BL_LEVELS, ROW_LENGTH, P_ROWS,      &
     &          LENOB, NPTOBT, NO_ANAL_LEVS, NO_ANAL_VAR,               &
     &          ICODE, CMESSAGE)

         ! Check for error - exit routine if occured
          IF (ICODE  >   0) GOTO 999

        END IF

       ! Move pointer to next obs type in group
        NPTOBT = NPTOBT + LENOBT

      endif !inobs == 0.or.lenobt == 0
      END DO

! Vectorisation elsewhere over all obs in group
      NPTOBT = 0
      LENOBT = LENOB

      IF (LDIAGAC) THEN
        IF (IACTF  /=  IACTL) THEN
         ! The group has >1 type in it, so group stats are worthwhile.
         ! Print mean & rms stats for all obs in group

         IF (LLDAC(2)) THEN
           CALL DIAGO ('AC', LACT(IACTF), 1, OBS_INCR, NORMF,           &
     &                OBS_LAT, OBS_LONG, LMISSD, LENOB, LENOBT, NPTOBT, &
     &                NO_ANAL_LEVS, NO_ANAL_VAR)

          END IF
        END IF
      END IF


!L 2.4 Do horizontal analysis (for weights)
          IPASS = 1

      IF (MODE_HANAL  ==  2) THEN   !  Use FI method
      if(lenob /= 0)then

        CALL FI (IPASS, TIMESTEP_NO, TIMESTEP, NAVT, JGROUP,            &
     &          NO_ANAL_VAR, LENOB, NO_ANAL_LEVS, NO_WT_LEVS, 1, LENMG, &
     &          ROW_LENGTH, P_ROWS, OBS_NO, OBS_LAT, OBS_LONG, OBS_TIME,&
     &          NORMF, OBS_INCR, CFPT, NPPT, INC_TYPE, MODEL_INCR,      &
     &          PSTAR,                                                  &
     &          P_FIELD, P_LEVELS, ICODE,                               &
     &          CMESSAGE)
      else
      do jlev=1,NO_ANAL_LEVS
      do j=1,lenmg
      MODEL_INCR(j,jlev,1)=0.0
      enddo
      enddo
      endif


       ! Check for error - exit from routine if occured
        IF (ICODE  >   0) GOTO 999

      END IF ! End of FI horizontal analysis step

      IF (LDIAGAC) THEN
        IF(mype == 0)PRINT *, 'END OF WTS ANALYSIS PHASE'

      END IF

! Save all levels of wts on model grid

      DO JLEV = 1, NO_WT_LEVS   !  Loop over no of weight levels.
      JJLEV=JLEV

       ! Save weights field

        STASH_ITEM_NO = 200+NAVT
        IF (SF(STASH_ITEM_NO)) THEN
          CALL AC_STASH (STASH_ITEM_NO,JJLEV, P_LEVELS, JGROUP,         &
     &      N_GROUPS, TIMESTEP_NO, STINDEX, STLIST, LEN_STLIST, SI, SF, &
     &      STASHWORK, STASH_LEVELS, NUM_STASH_LEVELS,                  &
     &      STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO, MODEL_INCR(1,JLEV,1),&
     &      LENMG, 'Weights Field', ICODE, CMESSAGE)

        END IF

       ! Get statistics on sum of weights
        IF (LDIAGAC .AND. LLDAC(4)) THEN
         IF (NAVT  ==  4) THEN
            CALL MMSPT (MODEL_INCR(1,JLEV,1), JJLEV, 0,                 &
     &                          'SUM OF WTS - RH ', ROW_LENGTH, P_ROWS, &
     &                           phi_p,L_regular)

          ELSE IF (NAVT  ==  5) THEN
            CALL MMSPT (MODEL_INCR(1,JLEV,1), JJLEV, 0,                 &
     &                          'SUM OF WTS - PR ', ROW_LENGTH, P_ROWS, &
     &                           phi_p,L_regular)
          ELSE      ! NAVT values other than 4 or 5 excluded
             ICODE = 1
             CMESSAGE = 'AC : Should not reach here (3)'
             GOTO 999   ! RETURN

          END IF
        END IF
      END DO   ! End of loop over no of weight levels.

!L 2.6 Calculate ob normalisation factors
! Vectorization here is over types within groups
      NPTOBT = 0

      DO JACT = IACTF, IACTL
      JJACT=JACT
        INOBS  = OBS_INFO % NOBS(JACT)
        LENOBT = LENACT(JACT)
!
!
!
      if(inobs == 0.or.lenobt == 0)then
      if(LDIAGAC)then
      if(LACT(JACT) /= 406.and.LACT(JACT) /= 506)then
      BCOUNT=2
      if(LACT(JACT) == 101.or.                                          &
     &   (LACT(JACT) >= 202.and.LACT(JACT) <= 204).or.                  &
     &   (LACT(JACT) >= 302.and.LACT(JACT) <= 306).or.                  &
     &   (LACT(JACT) >= 402.and.LACT(JACT) <= 404).or.                  &
     &   LACT(JACT) == 901)BCOUNT=1
!
! BARRIERS FOR DIAGO
!
      do bc=1,bcount
      do jlev=1,glsize(3,fld_type_p)
      do i=0,8
      s_stat(jlev,i)=0.0
      enddo
      enddo
! dummy diago call
       CALL DIAGO ('VAN?', LACT(JACT), 5+bc-bcount,                     &
     &              OBS_INCR, NORMF,                                    &
     &              OBS_LAT, OBS_LONG, LMISSD, LENOB, LENOBT, NPTOBT,   &
     &              NO_ANAL_LEVS, NO_ANAL_VAR)
      enddo
      endif ! LACT(JACT) /= 406 (and 506)
      endif ! LDIAGAC
      else ! lenob == 0

        IF (INOBS >  0 .AND. LENOBT >  0) THEN
          NDV   = OBS_INFO % NDATAV(JACT) - OBS_INFO % NDVHDR
          IPASS = 2

          CALL VERTANL (PSTGS,                                          &
     &           JJACT, IPASS, LENOBT, NDV, OBS_NO(NPTOBT +1),          &
     &           MODEL_OBS_TYPE(NPTOBT +1), OBS_LAT(NPTOBT +1),         &
     &           OBS_LONG(NPTOBT+1), OBS(OBS_INFO % MDISPOBT(JACT) +1), &
     &           INOBS, EXNER, PSTAR, THETA, RH, QCL, QCF, CONV_CLD,    &
     &           LAYER_CLOUD,PRESSURE,                                  &
     &           LS_RAIN, LS_SNOW, CONV_RAIN, CONV_SNOW,                &
     &           RHCRIT,                                                &
     &           P_FIELD, MODEL_INCR, LENMG, NO_ANAL_LEVS,              &
     &           CFPT(NPTOBT+1, 1), CFPT(NPTOBT+1, 2),                  &
     &           CFPT(NPTOBT+1, 3), CFPT(NPTOBT+1, 4),                  &
     &           NPPT(NPTOBT+1, 1), NPPT(NPTOBT+1, 2),                  &
     &           NPPT(NPTOBT+1, 3), NPPT(NPTOBT+1, 4),                  &
     &           OBS_INCR, NORMF, LMISSD,                               &
     &           P_LEVELS, Q_LEVELS, BL_LEVELS, ROW_LENGTH, P_ROWS,     &
     &           LENOB, NPTOBT, NO_ANAL_LEVS, NO_ANAL_VAR,              &
     &           ICODE, CMESSAGE)

         ! Check for errors:
          IF (ICODE  >   0) GOTO 999

        END IF

        NPTOBT = NPTOBT + LENOBT
      endif ! lenob == 0

      END DO

! Vectorisation elsewhere over all obs in group
      NPTOBT = 0
      LENOBT = LENOB

!L 3.0 Increments calculation starts here

      IF (LDIAGAC.AND.mype == 0) THEN
        PRINT '(/,A,(10I4))',                                           &
     &           ' START OF INCS ANALYSIS PHASE FOR TYPES  ',           &
     &           (LACT(J),J=IACTF,IACTL)

      END IF

!L 3.1 Horizontal analysis step (increments)
          IPASS = 2

      IF (MODE_HANAL  ==  2) THEN   ! Use FI method

      if(lenob /= 0)then
        CALL FI (IPASS, TIMESTEP_NO, TIMESTEP, NAVT, JGROUP,            &
     &         NO_ANAL_VAR, LENOB, NO_ANAL_LEVS, NO_WT_LEVS, 1, LENMG,  &
     &         ROW_LENGTH, P_ROWS, OBS_NO, OBS_LAT, OBS_LONG, OBS_TIME, &
     &         NORMF, OBS_INCR, CFPT, NPPT, INC_TYPE, MODEL_INCR,       &
     &         PSTAR,                                                   &
     &         P_FIELD, P_LEVELS, ICODE,                                &
     &         CMESSAGE)
      else
      do jlev=1,NO_ANAL_LEVS
      do j=1,lenmg
      MODEL_INCR(j,jlev,1)=0.0
      enddo
      enddo
      endif

       ! Check for errors:
        IF (ICODE  >   0) GOTO 999

      END IF ! End of FI horizontal analysis step

      IF (LDIAGAC) THEN
        IF(mype == 0)PRINT *, 'END OF INCS ANALYSIS PHASE'

      END IF

!L 3.3 Set up Relaxation Coefficients for this group
        N_ROWS = P_ROWS

      CALL SET_RELAX_CF (JGROUP, N_ROWS, RELAX_CF, LWIND, TIMESTEP,     &
     &     TIMESTEP_NO, ITER_NO, ICODE, CMESSAGE)

! Check for errors:
      IF (ICODE  >   0) GOTO 999

!L 3.4 Interpolation from analysis grid to model grid
!L     Loop over levels to interpolate incrs from analysis grid to
!L     model grid. Add incrs to model.
      DO JLEV = 1, NO_ANAL_LEVS
      JJLEV=JLEV

!L     3.4.1 Scale increment by relaxation coefficent

          CALL RELAXC (MODEL_INCR(1,JLEV,1), LENMG, ROW_LENGTH,         &
     &      P_ROWS, RELAX_CF, ICODE, CMESSAGE)

         ! Check for errors:
          IF (ICODE  >   0) GOTO 999


!L     3.4.2 Save increment fields into output file

          STASH_ITEM_NO = 210+NAVT
          IF (SF(STASH_ITEM_NO)) THEN
            CALL AC_STASH (STASH_ITEM_NO, JJLEV, P_LEVELS, JGROUP,      &
     &        N_GROUPS, TIMESTEP_NO, STINDEX, STLIST, LEN_STLIST, SI,   &
     &        SF, STASHWORK, STASH_LEVELS, NUM_STASH_LEVELS,            &
     &        STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,                    &
     &        MODEL_INCR(1,JLEV,1), LENMG, 'Increment Field', ICODE,    &
     &        CMESSAGE)

          END IF


!L     3.4.10 Humidity increments
        IF (NAVT  ==  4) THEN

         ! Add RH increments to RH fields
          CALL ADDINC (RH(1,JLEV), MODEL_INCR(1,JLEV,1), LENMG,         &
     &      ROW_LENGTH, P_ROWS, NAVT, ICODE, CMESSAGE)

         ! Check for errors:
          IF (ICODE  >   0) GOTO 999

         ! Statistics of increments on model grid
          IF (LDIAGAC .AND. LLDAC(4)) THEN
            CALL MMSPT (MODEL_INCR(1,JLEV,1), JJLEV, 0,                 &
     &        'RH INCREMENTS   ', ROW_LENGTH, P_ROWS,                   &
     &        phi_p,L_regular)

          END IF
        END IF   ! NAVT = 4

!L     3.4.11 latent heat nudging
        IF (NAVT  ==  5) THEN
!  Stash LS latent heating theta increments if required (rates in K/s)
            STASH_ITEM_NO = 271

            IF (SF(STASH_ITEM_NO)) THEN
              DO JL = 1 , Q_LEVELS
                JJL = JL
                CALL AC_STASH (STASH_ITEM_NO, JJL, Q_LEVELS,            &
     &            JGROUP, N_GROUPS, TIMESTEP_NO,                        &
     &            STINDEX, STLIST, LEN_STLIST, SI, SF, STASHWORK,       &
     &            STASH_LEVELS, NUM_STASH_LEVELS,                       &
     &            STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,                &
     &            D_THETA_DT_LS(1,JL),P_FIELD,                          &
     &            'LS latent heating rates',                            &
     &            ICODE, CMESSAGE)
              ENDDO    ! JL
            ENDIF

          IF (L_LHN) THEN
!C Call to new subroutine LHN_INC
      CALL LHN_INC(D_THETA_DT_CONV,D_THETA_DT_LS,                       &
     &             LS_RAIN,LS_SNOW,CONV_RAIN,                           &
     &             CONV_SNOW,P_FIELD,Q_LEVELS,MODEL_INCR,TIMESTEP,      &
     &             P_ROWS,ROW_LENGTH,                                   &
     &             LHN_RANGE, phi_p,lambda_p,L_regular,                 & 
     &             DRINCS,LENMG,ICODE,CMESSAGE)
!C Add on increments
            DO J=1,Q_LEVELS
              CALL ADDINC(THETA(1,J), DRINCS(1,J), P_FIELD, ROW_LENGTH, &
     &                    P_ROWS, NAVT,ICODE,CMESSAGE)
            ENDDO     !J

!  Stash LHN Theta Increments if required (amounts in K)
            STASH_ITEM_NO = 272

            IF (SF(STASH_ITEM_NO)) THEN
              DO JL = 1 , Q_LEVELS
                JJL = JL
                CALL AC_STASH (STASH_ITEM_NO, JJL, Q_LEVELS,            &
     &            JGROUP, N_GROUPS, TIMESTEP_NO,                        &
     &            STINDEX, STLIST, LEN_STLIST, SI, SF, STASHWORK,       &
     &            STASH_LEVELS, NUM_STASH_LEVELS,                       &
     &            STASH_PSEUDO_LEVELS, NUM_STASH_PSEUDO,                &
     &            DRINCS(1,JL),P_FIELD,'LHN theta increments',          &
     &            ICODE, CMESSAGE)
              ENDDO    ! JL
            ENDIF

            IF (LDIAGAC .AND. LLDAC(4) ) THEN
            DO JL = 1 , Q_LEVELS
            JJL = JL
              CALL MMSPT(DRINCS(1,JL), JJL, 0,                          &
     &                    'LHN Theta incrs ', ROW_LENGTH, P_ROWS,       &
     &                     phi_p,L_regular)
            ENDDO     ! JL
            ENDIF
          ENDIF       ! L_LHN

        ENDIF    ! NAVT = 5

      END DO   ! JLEV


!L 5.0 Exit from AC2
 999  CONTINUE
      IF (lhook) CALL dr_hook('AC2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE AC2
END MODULE ac2_mod
