! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE AC_INIT ------------------------------------------------
!LL
!LL  Purpose : Initialise and set up for assimilation.
!LL          : via calls to ACP_NAMEL,ACDIAG_NAMEL,NUM_OBS,SETTPS
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
MODULE ac_init_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE AC_INIT (P_LEVELS, Q_LEVELS, BL_LEVELS, TR_LEVELS,     &
     &                    P_ROWS, U_ROWS, ROW_LENGTH,                   &
     &                    TNDVMAX, NOBSMAX, TIMESTEP,                   &
     &                    BASIS_TIME_YY, BASIS_TIME_MM, BASIS_TIME_DD,  &
     &                    BASIS_TIME_HH, BASIS_TIME_MIN,                &
     &                    REALHD1, REALHD2, REALHD3,                    &
     &                    REALHD4, REALHD5, REALHD6,                    &
     &                    AK, BK, lambda_p,phi_p, no_lambda_p,          &
     &                    ICODE, CMESSAGE)

USE um_input_control_mod, ONLY : model_domain
USE domain_params, ONLY: mt_global
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE conversions_mod, ONLY: pi_over_180

USE filenamelength_mod, ONLY :                                          & 
    filenamelength
    
USE dynamics_input_mod, ONLY: L_regular
    
      USE UM_ParVars
      USE Control_Max_Sizes
      USE comobs_mod, ONLY: nobtypmx, obs_info
      USE acdiag_namel_mod, ONLY: acdiag_namel
      USE acp_namel_mod, ONLY: acp_namel
      USE num_obs_mod, ONLY: num_obs
      USE settps_mod, ONLY: settps
      IMPLICIT NONE
!*L --------------------- Comdeck: CENVIR   ----------------------------
!
!   Purpose: COMDECK defining Character enviroment variables used
!            by portable IO to open and close files
!
!
! Type declarations
!
      CHARACTER(LEN=8) FT_ENVIRON(199)  ! Array holding enviroment variables
!                                  for filenames
      INTEGER     LEN_FT_ENVIR(199) ! character length of each variable
!
!
!Common Blocks for character and integer arrays
!
      COMMON/CENVIR/FT_ENVIRON
      COMMON/CLENVIR/LEN_FT_ENVIR
!

      INTEGER P_LEVELS
      INTEGER Q_LEVELS
      INTEGER BL_LEVELS
      INTEGER TR_LEVELS
      INTEGER P_ROWS
      INTEGER U_ROWS
      INTEGER ROW_LENGTH
      INTEGER TNDVMAX
      INTEGER TNDVMAX_TOTAL
      INTEGER NOBSMAX
      INTEGER NOBSMAX_TOTAL
      REAL TIMESTEP
      INTEGER BASIS_TIME_YY
      INTEGER BASIS_TIME_MM
      INTEGER BASIS_TIME_DD
      INTEGER BASIS_TIME_HH
      INTEGER BASIS_TIME_MIN
      REAL REALHD1,REALHD2,REALHD3,REALHD4,REALHD5,REALHD6
      REAL AK(P_LEVELS),BK(P_LEVELS)
      CHARACTER(LEN=256) CMESSAGE
      INTEGER ICODE

      integer no_lambda_p
      Real lambda_p(no_lambda_p)
      Real phi_p(1-halo_i:row_length+halo_i, 1-halo_j:p_rows+halo_j)

!-INTENT=IN--------------------------------------------------------
!     P_LEVELS      : TOTAL NUMBER OF LEVELS
!     Q_LEVELS      : TOTAL NUMBER OF wet LEVELS
!     BL_LEVELS     : TOTAL NUMBER OF boundary layer LEVELS
!     TR_LEVELS     : TOTAL NUMBER OF tracer LEVELS
!     ROW_LENGTH    : NUMBER OF POINTS ON ROW
!     P_ROWS        : NUMBER OF ROWS (FOR PSTAR)
!     U_ROWS        : NUMBER OF ROWS (FOR wind)
!     TIMESTEP      : TIMESTEP IN SECONDS
!     BASIS_TIME_## : DEFINES DATA TIME
!     REALHD1-6     : DEFINES HORIZONTAL GRID
!     AK,BK         : DEFINES VERTICAL GRID
!-INTENT=OUT-----------------------------------------------------
!     TNDVMAX       : MAX SIZE OF OBS ARRAY
!     NOBSMAX       : MAX NUMBER OF OBS
!     ICODE         : NON ZERO FOR FAILURE
!     CMESSAGE      : REASON FOR FAILURE
!*---------------------------------------------------------------------
      EXTERNAL GET_FILE
!L---------------------------------------------------------------------

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
!--------------------------------------------------------------------
!LCOMDECK COMMG
!L-------------
      REAL          DLAT,DLONG,XLATN,XLONGW
      REAL          ELFPLAT,ELFPLON
      COMMON/COMMG/ DLAT,DLONG,XLATN,XLONGW,ELFPLAT,ELFPLON
!--------------------------------------------------------------------

! Other comdecks
! ----------------------- Header file CRUNTIMC  -----------------------
! Description: Run-time constants for the Atmosphere model (read only).
!              Contains variables that define parametrization values
!              chosen for atmosphere physics and dynamics schemes.
!              [Note that CNTLATM holds accompanying control switches
!              needed for addressing.]
!
! This file belongs in section: Top Level

!
!------------------   Physics:   --------------------------------------
! Generalised physics switches:

!------------------   End of Physics   ---------------------------------

      INTEGER :: rpemax ! array size needed for diagnostic printing
      INTEGER :: rpemin ! array size needed for diagnostic printing
      INTEGER :: rpesum ! array size needed for diagnostic printing
      INTEGER :: ipesum ! array size needed for diagnostic printing
      INTEGER :: time_theta1_min      ! Timestep of min level 1 theta
      INTEGER :: time_w_max(model_levels_max) ! Timestep of max w
      INTEGER :: time_div_max(model_levels_max) ! Timestep of max div
      INTEGER :: time_div_min(model_levels_max) ! Timestep of min div
      INTEGER :: time_lapse_min(model_levels_max) ! Timestep of min
      INTEGER :: time_max_shear(model_levels_max) !Timestep max shear
      INTEGER :: time_max_wind(model_levels_max) ! Timestep of max wind
      INTEGER :: time_KE_max(model_levels_max) ! Timestep of max KE
      INTEGER :: time_KE_min(model_levels_max) ! Timestep of min KE
      INTEGER :: time_noise_max(model_levels_max) ! Timestep of max

      REAL:: frictional_timescale(model_levels_max) ! For idealised case
      REAL :: tropics_deg  ! define latitude for tropics
      REAL :: min_theta1_run                  ! Min theta level 1
      REAL :: dtheta1_run   ! Largest -ve delta theta at min theta1
      REAL :: max_w_run(0:model_levels_max) ! Max w at a level
      REAL :: max_div_run(model_levels_max) ! Max divergence at a level
      REAL :: min_div_run(model_levels_max) ! Min divergence at a level
      REAL :: min_lapse_run(model_levels_max) ! Min dtheta/dz at a level
      REAL :: max_shear_run(model_levels_max) ! Max shear at a level
      REAL :: max_wind_run(model_levels_max) ! Max wind at a level
      REAL :: max_KE_run(model_levels_max)   ! Max KE at a level
      REAL :: min_KE_run(model_levels_max)   ! Min KE at a level
      REAL :: max_noise_run(model_levels_max) ! Max noise at a level

!     Problem_number not set here  ! Now controlled by namelist input
!     Instability_diagnostics      ! Now controlled by namelist input
!     frictional_timescale         ! Now intitialised in SETCONA
!------------------   Diagnostics:   --------------------------------

      COMMON  /RUN_Diagnostics/                                         &
        rpemax, rpemin, ipesum, rpesum,                                 &
        max_w_run, min_theta1_run, dtheta1_run,                         &
        max_div_run, min_div_run, min_lapse_run,                        &
        max_shear_run, max_wind_run,                                    &
        max_noise_run, max_KE_run, min_KE_run,                          &
        time_KE_max, time_KE_min,                                       &
        time_w_max, time_div_max, time_div_min, time_lapse_min,         &
        time_max_shear, time_max_wind,                                  &
        time_theta1_min, time_noise_max

!------------------   Dynamics:   --------------------------------------
! Suarez-Held variables
      REAL :: SuHe_newtonian_timescale_ka
      REAL :: SuHe_newtonian_timescale_ks
      REAL :: SuHe_pole_equ_deltaT
      REAL :: SuHe_static_stab
      REAL :: base_frictional_timescale
      REAL :: SuHe_sigma_cutoff
      REAL :: SuHe_level_weight(model_levels_max)
      REAL :: friction_level(model_levels_max)

      INTEGER :: SuHe_relax
      INTEGER :: SuHe_fric

      LOGICAL :: L_SH_Williamson

      COMMON/Run_Dyncore/                                              &
       SuHe_newtonian_timescale_ka, SuHe_newtonian_timescale_ks,       &
       SuHe_pole_equ_deltaT, SuHe_static_stab,                         &
       base_frictional_timescale, SuHe_sigma_cutoff,                   &
       L_SH_Williamson, SuHe_relax, SuHe_fric,                         &
       SuHe_level_weight, frictional_timescale, friction_level

!------------------  Idealised model   ----------------------------

      INTEGER,PARAMETER:: max_num_profile_data = 100
      INTEGER,PARAMETER:: max_num_force_times = 100
      INTEGER,PARAMETER:: idl_max_num_bubbles = 3

! Idealised  variables
      REAL :: h_o
      REAL :: h_o_actual  ! height of growing mountain
      REAL :: h_o_per_step ! height change per step of growing mountain
      REAL :: lambda_fraction
      REAL :: phi_fraction
      REAL :: half_width_x
      REAL :: half_width_y
      REAL :: Witch_power
      REAL :: plat_size_x
      REAL :: plat_size_y
      REAL :: height_domain
      REAL :: delta_x
      REAL :: delta_y
      REAL :: big_factor
      REAL :: mag
      REAL :: vert_grid_ratio
      REAL :: first_theta_height
      REAL :: thin_theta_height
      REAL :: p_surface
      REAL :: theta_surface
      REAL :: dtheta_dz1(3)
      REAL :: height_dz1(3)
      REAL :: Brunt_Vaisala
      REAL :: u_in(4)
      REAL :: v_in(4)
      REAL :: height_u_in(3)
      REAL :: u_ramp_start
      REAL :: u_ramp_end
      REAL :: ujet_lat
      REAL :: ujet_width
      REAL :: t_horizfn_data(10)
      REAL :: q1
      REAL :: theta_ref(model_levels_max)
      REAL :: rho_ref(model_levels_max)
      REAL :: exner_ref(model_levels_max + 1)
      REAL :: q_ref(model_levels_max)
      REAL :: u_ref(model_levels_max)
      REAL :: v_ref(model_levels_max)
      REAL :: z_orog_print(0:model_levels_max)
      REAL :: f_plane
      REAL :: ff_plane
      REAL :: r_plane
      REAL :: zprofile_data(max_num_profile_data)
      REAL :: tprofile_data(max_num_profile_data)
      REAL :: qprofile_data(max_num_profile_data)
      REAL :: z_uvprofile_data(max_num_profile_data)
      REAL :: uprofile_data(max_num_profile_data)
      REAL :: vprofile_data(max_num_profile_data)
      REAL :: tforce_time_interval
      REAL :: qforce_time_interval
      REAL :: uvforce_time_interval
      REAL :: newtonian_timescale
      REAL :: z_tforce_data(max_num_profile_data)
      REAL :: z_qforce_data(max_num_profile_data)
      REAL :: z_uvforce_data(max_num_profile_data)
      REAL :: tforce_data(max_num_profile_data, max_num_force_times)
      REAL :: qforce_data(max_num_profile_data, max_num_force_times)
      REAL :: uforce_data(max_num_profile_data, max_num_force_times)
      REAL :: vforce_data(max_num_profile_data, max_num_force_times)
      REAL :: tforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: qforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: uforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: vforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: pforce_time_interval
      REAL :: p_surface_data(max_num_force_times)
      REAL :: perturb_factor
      REAL :: perturb_magnitude_t
      REAL :: perturb_magnitude_q
      REAL :: perturb_height(2)
      REAL :: orog_hgt_lbc
      REAL :: zprofile_orog
      REAL :: hf
      REAL :: cool_rate
      REAL :: IdlSurfFluxSeaParams(10) ! Idealised surface flux params
      REAL :: roughlen_z0m   
      REAL :: roughlen_z0h
      ! Idealised bubbles
      REAL :: idl_bubble_max(idl_max_num_bubbles) ! Bubble max amplitude
      REAL :: idl_bubble_height(idl_max_num_bubbles)  ! Bubble height
      REAL :: idl_bubble_width(idl_max_num_bubbles)   ! Bubble width
      REAL :: idl_bubble_depth(idl_max_num_bubbles)   ! Bubble depth
      ! Bubble x-offset, y-offset in normalised units (0:1)
      ! (0.5=domain centre)
      REAL :: idl_bubble_xoffset(idl_max_num_bubbles)
      REAL :: idl_bubble_yoffset(idl_max_num_bubbles)
      REAL :: DMPTIM, HDMP, ZDMP   ! Damping layer values
      REAL :: u_geo, v_geo         ! Geostrophic wind

! ENDGAME
      REAL :: T_surface
      REAL :: Eccentricity
      ! Following two variables used only if L_rotate_grid=.true.
      REAL :: grid_NP_lon ! Longitude (degrees) of grid's north pole
      REAL :: grid_NP_lat ! Latitude (degrees) of grid's north pole
      REAL :: AA_jet_u0   ! See QJRMS 133,1605--1614
      REAL :: AA_jet_A    !
      REAL :: theta_pert
      REAL :: ring_height
      REAL :: angular_velocity ! Planet's angular velocity
      REAL :: T0_P, T0_E ! deep atmosphere baroclinic wave surface temperatures
      INTEGER :: Trefer_number
      INTEGER :: tstep_plot_frequency
      INTEGER :: tstep_plot_start
      INTEGER :: AA_jet_m  ! See QJRMS 133,1605--1614
      INTEGER :: AA_jet_n  !
      INTEGER :: chain_number ! Run continuation number
      LOGICAL :: L_rotate_grid    ! .true. for rotating North pole of grid
      LOGICAL :: L_baro_Perturbed ! Used for baroclinic test to specify
                                  ! pert or steady
      LOGICAL :: L_shallow, L_const_grav, L_HeldSuarez,L_HeldSuarez1_drag,  &
                 L_HeldSuarez2_drag,                                        &
                 L_baro_inst, L_isothermal, L_exact_profile, L_balanced,    &
                 L_solid_body
      LOGICAL :: L_deep_baro_inst ! deep atmosphere baroclinic wave switch          


      INTEGER :: surface_type
      INTEGER :: grow_steps
      INTEGER :: grid_number
      INTEGER :: grid_flat
      INTEGER :: tprofile_number
      INTEGER :: qprofile_number
      INTEGER :: uvprofile_number
      INTEGER :: num_profile_data
      INTEGER :: num_uvprofile_data
      INTEGER :: t_horizfn_number
      INTEGER :: uv_horizfn_number
      INTEGER :: pforce_option
      INTEGER :: num_pforce_times
      INTEGER :: tforce_option
      INTEGER :: qforce_option
      INTEGER :: uvforce_option
      INTEGER :: num_tforce_levels
      INTEGER :: num_tforce_times
      INTEGER :: num_qforce_levels
      INTEGER :: num_qforce_times
      INTEGER :: num_uvforce_levels
      INTEGER :: num_uvforce_times
      INTEGER :: IdlSurfFluxSeaOption  ! Idealised surface flux option
      INTEGER :: first_constant_r_rho_level_new
      INTEGER :: big_layers
      INTEGER :: transit_layers
      INTEGER :: mod_layers
      INTEGER :: idl_bubble_option(idl_max_num_bubbles) ! Bubble option
      INTEGER :: idl_interp_option  ! Profile interpolation option
      INTEGER :: perturb_type
      INTEGER :: b_const, k_const ! deep atmosphere baroclinic wave parameters

      LOGICAL :: L_initialise_data
      LOGICAL :: L_constant_dz
      LOGICAL :: L_trivial_trigs !.false. for Cartesian coords (lat=0.0)
      LOGICAL :: L_idl_bubble_saturate(idl_max_num_bubbles)
      LOGICAL :: L_fixed_lbcs
      LOGICAL :: L_fix_orog_hgt_lbc
      LOGICAL :: L_pressure_balance
      LOGICAL :: L_wind_balance
      LOGICAL :: L_rotate_winds
      LOGICAL :: L_polar_wind_zero
      LOGICAL :: L_vert_Coriolis
      LOGICAL :: L_rotating     ! .true. for Earth's rotation
      LOGICAL :: L_perturb      ! add random perturb. to surface theta
      LOGICAL :: L_code_test    ! User switch for testing code
      LOGICAL :: L_pforce
      LOGICAL :: L_baroclinic
      LOGICAL :: L_cyclone
      LOGICAL :: L_force
      LOGICAL :: L_force_lbc
      LOGICAL :: L_perturb_t
      LOGICAL :: L_perturb_q
      LOGICAL :: L_perturb_correlate_tq
      LOGICAL :: L_perturb_correlate_vert
      LOGICAL :: L_perturb_correlate_time
      LOGICAL :: L_damp      ! Logical for damping layer
      LOGICAL :: L_geo_for ! Logical for geostrophic wind forcing
      LOGICAL :: L_bomex     ! Logical for BOMEX set up
      LOGICAL :: L_spec_z0   ! specification of roughness length    

      COMMON  /RUN_Ideal/                                              &
       h_o, h_o_actual, h_o_per_step,                                  &
       lambda_fraction, phi_fraction, half_width_x, half_width_y,      &
       Witch_power, plat_size_x, plat_size_y,                          &
       height_domain, delta_x, delta_y, big_factor, mag, vert_grid_ratio, &
       first_theta_height, thin_theta_height, p_surface,               &
       theta_surface, dtheta_dz1, height_dz1, Brunt_Vaisala,           &
       u_in, v_in, height_u_in, u_ramp_start, u_ramp_end, q1,          &
       ujet_lat, ujet_width,                                           &
       t_horizfn_number, t_horizfn_data, uv_horizfn_number,            &
       u_ref, v_ref, theta_ref, exner_ref, rho_ref, q_ref,             &
       z_orog_print, grow_steps,                                       &
       surface_type, grid_number, grid_flat,                           &
       tprofile_number, qprofile_number, uvprofile_number,             &
       num_profile_data, num_uvprofile_data,                           &
       tforce_option, qforce_option, uvforce_option,                   &
       num_tforce_levels, num_tforce_times,                            &
       num_qforce_levels, num_qforce_times,                            &
       num_uvforce_levels, num_uvforce_times,                          &
       L_pforce, pforce_option, num_pforce_times,                      &
       first_constant_r_rho_level_new,                                 &
       big_layers, transit_layers, mod_layers,                         &
       zprofile_data, tprofile_data, qprofile_data,                    &
       z_uvprofile_data, uprofile_data, vprofile_data,                 &
       tforce_time_interval, qforce_time_interval,                     &
       uvforce_time_interval, newtonian_timescale,                     &
       z_tforce_data, z_qforce_data, z_uvforce_data,                   &
       tforce_data, qforce_data, uforce_data, vforce_data,             &
       tforce_data_modlev, qforce_data_modlev,                         &
       uforce_data_modlev, vforce_data_modlev,                         &
       pforce_time_interval, p_surface_data,                           &
       L_initialise_data,                                              &
       L_perturb_t, perturb_magnitude_t,                               &
       L_perturb_q, perturb_magnitude_q,                               &
       L_perturb_correlate_tq,                                         &
       L_perturb_correlate_vert,                                       &
       L_perturb_correlate_time,                                       &
       perturb_type, perturb_height,                                   &
       L_constant_dz, L_polar_wind_zero,                               &
       L_wind_balance, L_rotate_winds,                                 &
       IdlSurfFluxSeaOption, IdlSurfFluxSeaParams,                     &
       L_spec_z0, roughlen_z0m, roughlen_z0h,                          &
       L_pressure_balance, L_vert_Coriolis,                            &
       cool_rate, L_force, L_force_lbc,                                &
       zprofile_orog, idl_interp_option, hf,                           &
       L_fix_orog_hgt_lbc, orog_hgt_lbc,                               &
       L_trivial_trigs, f_plane, ff_plane, r_plane,                    &
       idl_bubble_option, idl_bubble_max                               &
      , idl_bubble_height, idl_bubble_width, idl_bubble_depth          &
      , idl_bubble_xoffset,idl_bubble_yoffset                          &
      , L_idl_bubble_saturate,                                         &
       L_rotating, L_fixed_lbcs, L_code_test,                          &
       L_baroclinic, L_cyclone,                                        &
       L_damp, L_geo_for, L_bomex,                                     &
       DMPTIM, HDMP, ZDMP,                                             &
       u_geo, v_geo,                                                   &
!ENDGAME
       T_surface, chain_number,                                        &
       Trefer_number,                                                  &
       tstep_plot_frequency, tstep_plot_start, Eccentricity,           &
       L_rotate_grid, grid_NP_lon, grid_NP_lat,                        &
       AA_jet_m, AA_jet_n, AA_jet_u0, AA_jet_A, L_baro_Perturbed,      &
       L_shallow, L_const_grav, L_HeldSuarez,L_HeldSuarez1_drag,       &
       L_HeldSuarez2_drag,                                             &
       L_baro_inst, L_deep_baro_inst, T0_P, T0_E, b_const, k_const,    &      
       ring_height, theta_pert, L_isothermal,                          &
       L_exact_profile, L_balanced, L_solid_body, angular_velocity
! CRUNTIMC end

!-----------------------------------------------------------------------
      INTEGER JFILE
      INTEGER ICCODE
      INTEGER IUNIT7
      CHARACTER(LEN=filenamelength) :: filename
      REAL, ALLOCATABLE :: PHI_P_GLOBAL(:,:)
      INTEGER ISTAT,i,j
      REAL PHI_P_LOCAL(row_length,p_rows)
      REAL XLATS,XLONGE ! used in variable grid runs
      REAL DLAT_N,DLAT_S,DLON_W,DLON_E
      INTEGER ROW_LENGTH_GLOBAL,P_ROWS_GLOBAL

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('AC_INIT',zhook_in,zhook_handle)

      ROW_LENGTH_GLOBAL=glsize(1,fld_type_p)
      P_ROWS_GLOBAL=glsize(2,fld_type_p)


!-----------------------------------------------------------------------

      if(mype == 0) PRINT *, ' IN AC_INIT'

!L 1. Initialise Unit Numbers used in AC
! Unit No for first AC Observation File
      OBS_UNITNO = 70

! Unit No for CACHE FILE to store observations between timesteps
! Stores arrays OBS and OBS_FLAG
      IUNITNO = 15

!     Unit No for printed output used to inform operators of
!     data problems
      IUNIT7 = 7

!L 2. Set up COMMG Variables and check dimensions
! For regular grids read values from header otherwise from coordinate arrays
      IF(L_regular) THEN

        DLONG  = REALHD1
        DLAT   = REALHD2
        XLATN  = REALHD3 + (P_ROWS_GLOBAL-1)*DLAT
        XLONGW = REALHD4

       ELSE

! Gather phi and lambda coordinates for entire domain
         DO I=1,P_ROWS
           DO J=1,ROW_LENGTH
             PHI_P_LOCAL(J,I)=PHI_P(J,I)
           ENDDO
         ENDDO

        ALLOCATE (PHI_P_GLOBAL(row_length_global,p_rows_global))

! DEPENDS ON: gather_field
        CALL GATHER_FIELD(PHI_P_LOCAL,   PHI_P_GLOBAL,                  &
     &                    ROW_LENGTH,    P_ROWS,                        &
     &                    row_length_global,                            &
     &                    p_rows_global,                                &
     &                    fld_type_p,    halo_type_no_halo,             &
     &                    0,             gc_all_proc_group,             &
     &                    ICODE,         CMESSAGE)

! In a variable grid set boundaries and DLAT and DLONG to be 
! outer values
! Calculate on PE0 then broadcast
! DLAT and DLON are the outermost grid lengths
        IF(mype == 0) THEN
          DLONG = (LAMBDA_P(2)-LAMBDA_P(1))/PI_OVER_180
          DLAT = (PHI_P_GLOBAL(1,2)-PHI_P_GLOBAL(1,1))/PI_OVER_180
          XLATN = PHI_P_GLOBAL(1,P_ROWS_GLOBAL)/PI_OVER_180
          XLONGW = LAMBDA_P(1)/PI_OVER_180
          XLATS = PHI_P_GLOBAL(1,1)/PI_OVER_180
          XLONGE = LAMBDA_P(ROW_LENGTH_GLOBAL)/PI_OVER_180
        ENDIF

! Broadcast to all PEs
        CALL GC_RBCAST(1000,1,0,NPROC,ISTAT,DLONG)
        CALL GC_RBCAST(1001,1,0,NPROC,ISTAT,DLAT)
        CALL GC_RBCAST(1002,1,0,NPROC,ISTAT,XLATN)
        CALL GC_RBCAST(1003,1,0,NPROC,ISTAT,XLATS)
        CALL GC_RBCAST(1004,1,0,NPROC,ISTAT,XLONGW)
        CALL GC_RBCAST(1005,1,0,NPROC,ISTAT,XLONGE)

      ENDIF  

      IF(.NOT.L_regular) THEN
        DLAT_N=(PHI_P(1,P_ROWS)-PHI_P(1,P_ROWS-1))/PI_OVER_180
        DLAT_S=(PHI_P(1,2)-PHI_P(1,1))/PI_OVER_180
        DLON_W=(LAMBDA_P(2)-LAMBDA_P(1))/PI_OVER_180
        DLON_E=(LAMBDA_P(ROW_LENGTH)-LAMBDA_P(ROW_LENGTH-1))/PI_OVER_180
      ENDIF

      IF(P_LEVELS >   MODEL_LEVELS_MAX)THEN
       ICODE=1
       CMESSAGE = ' ACINIT: MODEL_LEVELS_MAX too small'
       if(mype == 0) PRINT *,  CMESSAGE,                                &
     &      ' P_LEVELS=',P_LEVELS,' MODEL_LEVELS_MAX=',MODEL_LEVELS_MAX
       GOTO 999
      ENDIF

      IF(Q_LEVELS >   WET_LEVELS_MAX)THEN
       ICODE=1
       CMESSAGE = ' ACINIT: WET_LEVELS_MAX too small'
       if(mype == 0) PRINT *,  CMESSAGE,                                &
     &      ' Q_LEVELS=',Q_LEVELS,' WET_LEVELS_MAX=',WET_LEVELS_MAX
       GOTO 999
      ENDIF

      IF(P_ROWS >   ROWS_MAX)THEN
       ICODE=1
       CMESSAGE = ' ACINIT: ROWS_MAX too small'
       if(mype == 0) PRINT *,  CMESSAGE,                                &
     &      ' P_ROWS=',P_ROWS,' ROWS_MAX=',ROWS_MAX
       GOTO 999
      ENDIF

      IF(ROW_LENGTH >   ROW_LENGTH_MAX)THEN
       ICODE=1
       CMESSAGE = ' ACINIT: ROW_LENGTH_MAX too small'
       if(mype == 0) PRINT *,  CMESSAGE,                                &
     &      ' ROW_LENGTH=',ROW_LENGTH,' ROW_LENGTH_MAX=',ROW_LENGTH_MAX
       GOTO 999
      ENDIF

! Make sure western boundary longitude in range 0 - 360 degrees
      IF (XLONGW <  0.0) XLONGW = XLONGW + 360.0

! Real lat/lon of pseudo N. pole in degrees

      IF (model_domain /= mt_global) THEN
      ELFPLAT = REALHD5
      ELFPLON = REALHD6
      END IF  ! .NOT. GLOBAL

      WRITE (6,'(/,A,(T25,10F10.6))') ' DLAT,DLONG,XLATN,XLONGW',       &
     &                                  DLAT,DLONG,XLATN,XLONGW

!L 3. ACP namelist. Set defaults, read in and process.
      CALL ACP_NAMEL (P_LEVELS, Q_LEVELS, BL_LEVELS, TR_LEVELS,         &
     &  P_ROWS, U_ROWS, ROW_LENGTH,                                     &
     &  TIMESTEP, ICODE, CMESSAGE)

      IF (ICODE >  0) GO TO 999

!L 4. ADIAG namelist. Set defaults, read in and process.
      CALL ACDIAG_NAMEL (ICODE,CMESSAGE)
      IF (ICODE >  0) GO TO 999

!L 5. Open Cache file (unit 15)

      if(mype == 0)then
      CALL GET_FILE(IUNITNO,FILENAME,filenamelength,ICODE)
        OPEN(IUNITNO,FILE=FILENAME,FORM='UNFORMATTED')
      CALL GET_FILE(IUNIT7,FILENAME,filenamelength,ICODE)
        OPEN(IUNIT7,FILE=FILENAME)
      endif

!L 6. Read in AC Obs Files and compute NOBSMAX and TNDVMAX

      CALL NUM_OBS (NO_OBS_FILES,NOBSMAX,TNDVMAX,P_LEVELS,Q_LEVELS,     &
     &              P_ROWS,ROW_LENGTH,AK,BK,REALHD1,REALHD2,            &
     &              REALHD3,REALHD4,REALHD5,REALHD6,                    &
     &              ICODE,CMESSAGE)
      IF (ICODE >  0) GO TO 999

!L 7. Set up list of obs types and groups to be processed this run
      CALL SETTPS (ICODE,CMESSAGE)
      IF (ICODE >  0) GO TO 999

!L 8. Set up Analysis Grid in COMAG
! Initialise COMAG variables which don't change during a run.

! Model Grid Spacings
      DLATMG  = DLAT *PI_OVER_180
      DLONGMG = DLONG*PI_OVER_180

! First row/point on p*/theta grid
      ROW1MGTH = ( 90.0 -XLATN )*PI_OVER_180
      PT1MGTH  = XLONGW*PI_OVER_180

! Convert Lat at which pts/row starts decreasing to co-lat/radians.
      AGLATDEC = (90.0-AGLATDEC)*PI_OVER_180

! Width of area according to model grid dimensions.
      AGROWLEN = DLONGMG*(ROW_LENGTH_GLOBAL-1)


! The COMAG variables DEF_AGRES_ROWS and DEF_AGRES_PTS are
! initialised in DEF_GROUP. Use &ACP namelist arrays AGRES_ROWS
! and AGRES_PTS to change initialised values.

! The remaining COMAG variables are set in SETAG.

!L 9. Set up COMOBS Variables
! Time Interval (mins) between reading AC Obs files
      OBS_INFO % TIMEINT = 180.0

! Reference Time/Date which is start of assimilation
      OBS_INFO % OBS_REF_YY  = BASIS_TIME_YY
      OBS_INFO % OBS_REF_MM  = BASIS_TIME_MM
      OBS_INFO % OBS_REF_DD  = BASIS_TIME_DD
      OBS_INFO % OBS_REF_HH  = BASIS_TIME_HH
      OBS_INFO % OBS_REF_MIN = BASIS_TIME_MIN

! Set up Latitudes and Longitudes of grid boundaries
! for observations to be used in assimilation.

! For Limited area assimilations, reject observations within
! one grid length of boundary.
      IF(L_regular) THEN
        OBS_INFO % OBS_LAT_N  = XLATN + 0.01 - DLAT
        OBS_INFO % OBS_LAT_S  = XLATN - 0.01 - (P_ROWS_GLOBAL-2)*DLAT
      ELSE
        OBS_INFO % OBS_LAT_N = XLATN + 0.01 - DLAT_N
        OBS_INFO % OBS_LAT_S = XLATS - 0.01 + DLAT_S
      ENDIF
      IF(L_regular) THEN
        OBS_INFO % OBS_LONG_W = XLONGW - 0.01 + DLONG
        OBS_INFO % OBS_LONG_E = XLONGW + 0.01 + (ROW_LENGTH_GLOBAL-2)*DLONG
      ELSE
        OBS_INFO % OBS_LONG_W = XLONGW - 0.01 + DLON_W
        OBS_INFO % OBS_LONG_E = XLONGE + 0.01 - DLON_E
      ENDIF

! After rotation of the Lat/Long of Obs to ELF co-ords,
! longitude values will be in range 0-360 degrees, so
! make sure that boundary values are consistent. Note that
! if this leaves ZLONMN > ZLONMX, the test for obs in the area
! will assume that the Limited Area grid straddles the Meridian
! between ZLONMN and ZLONMX.
      IF (OBS_INFO % OBS_LONG_W  <   0.0) THEN
        OBS_INFO % OBS_LONG_W = OBS_INFO % OBS_LONG_W + 360.0

      ELSEIF (OBS_INFO % OBS_LONG_W  >   360.0) THEN
        OBS_INFO % OBS_LONG_W = OBS_INFO % OBS_LONG_W - 360.0

      ENDIF

      IF (OBS_INFO % OBS_LONG_E  <   0.0) THEN
        OBS_INFO % OBS_LONG_E = OBS_INFO % OBS_LONG_E + 360.0

      ELSEIF (OBS_INFO % OBS_LONG_E  >   360.0) THEN
        OBS_INFO % OBS_LONG_E = OBS_INFO % OBS_LONG_E - 360.0

      ENDIF

! Initialise time for next read in of observation files.
! Forces read on first timestep.
      OBS_INFO % TIMENEXT = -1440.0

 999  CONTINUE
      IF (lhook) CALL dr_hook('AC_INIT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE AC_INIT

END MODULE ac_init_mod
