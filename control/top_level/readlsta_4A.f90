! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Read run-time control information from namelists for atmos model
!  ENDGAME VERSION
!
! Subroutine Interface:
      SUBROUTINE Readlsta_4A()

      USE umsections_mod, ONLY: atmos_sr
      USE switches, ONLY: l_flake_model 
      USE rad_input_mod
      USE sw_rad_input_mod
      USE lw_rad_input_mod
      USE MAX_CALLS, ONLY : npd_swcall, npd_lwcall
      USE clmchfcg_scenario_mod, ONLY: L_clmchfcg, clim_fcg_nyears,     &
                   clim_fcg_rates, nclmfcgs, clmchfcg
      USE CORADOCA

      USE um_input_control_mod, ONLY:                                   &
           l_oasis
      USE switches, ONLY: l_ssice_albedo
      USE cv_run_mod                         ! Access to all variables
      USE cv_param_mod                       ! Access to all variables
      USE check_iostat_mod
      
      USE missing_data_mod, ONLY: rmdi, imdi
      
      USE sl_input_mod
      USE Submodel_Mod

! module for RUN_PRECIP namelist
      USE mphys_inputs_mod, ONLY:                                       &
      ! Logicals block 1
        l_cry_agg_dep, l_it_melting, l_autoc_3b, l_warm_new,            &
        l_autolim_3b,  l_psd, l_psd_global, l_autoconv_murk,            &
      ! cloud rain correlation coefficient
        c_r_correl,                                                      &
      ! Drop-size distribution parameters
        x1r, x2r,                                                       &
        ai,  bi,  aic, bic,                                             &
        lsp_eic, lsp_fic,                                               &
      ! Number of iterations of microphysics
        lsiter,                                                         &
      ! Droplet taper parameters
        z_peak_nd, ndrop_surf,                                          &
      ! Logicals block 2
        l_droplet_tpr, l_rainfall_as,                                   &
        l_mcr_iter,      l_mcr_qrain,      l_mcr_qgraup,                &
        l_mcr_qrain_lbc, l_mcr_qgraup_lbc, l_rain,                      &
        RUN_PRECIP

! module for RUN_CLOUD namelist
      USE cloud_inputs_mod, ONLY: pc2_falliceshear_method,              &
       cloud_fraction_method,  i_fixbug_pc2_checks,                     &
       i_pc2_conv_coupling, i_pc2_erosion_method, L_eacf,               &
       l_ensure_min_in_cloud_qcf, l_fixbug_pc2_qcl_incr,                &
       l_fixbug_pc2_mixph, l_micro_eros,                                &
       dbsdtbs_turb_0,                                                  &
       starticeTKelvin, alliceTdegC, cff_spread_rate, ice_width, rhcrit,&
       l_pc2, l_rhcpt, RUN_CLOUD, check_run_cloud,                      &
       l_filter_cloud, tau_thresh

! module for RUN_RIVERS namelist
      USE river_inputs_mod, ONLY: RUN_Rivers

! module for RUN_Eng_Corr namelist
      USE eng_corr_inputs_mod, ONLY: RUN_Eng_Corr

! module for RUN_Murk namelist
      USE murk_inputs_mod, ONLY: RUN_Murk

! module for RUN_LAND, RUN_BLVEG and RUN_PFT namelists
      USE land_surf_mod

! module for mineral dust scheme options
      USE dust_parameters_mod, ONLY:                                    &
         RUN_Dust, dust_parameters_check, dust_parameters_load,         &
         dust_size_dist_initialise

! module for aerosol emission options
      USE run_aerosol_mod, ONLY: RUN_Aerosol
 
! module for Boundary Layer options
      USE bl_option_mod, ONLY: run_bl, cor_mo_iter, ishear_bl,          &
           alpha_cd_batches, alpha_cd_items, alpha_cd_vals, alpha_cd,   &
           l_quick_ap2, l_lambdam2, l_full_lambdas, nl_bl_levels,       &
           seasalinityfactor, iseaz0t, puns, pstb, sbl_op,              &
           variable_ric, idyndiag, isrfexcnvgust, fric_heating,         &
           subs_couple_fix, l_us_blsol, bl_option_off=>off, check_run_bl

! module for UKCA options
      USE ukca_option_mod, ONLY: run_ukca,                        &
         l_ukca, l_ukca_aie1, l_ukca_aie2,                        &
         i_ukca_chem, L_ukca_achem,                               & 
         L_ukca_mode, L_ukca_dust, L_ukca_ageair,                 &
         L_ukca_qch4inter,                                        &
         L_ukca_useumuivals,                                      &
         L_ukca_het_psc, L_ukca_sa_clim,                          &
         L_ukca_h2o_feedback,                                     &
         L_ukca_rado3, L_ukca_radch4, L_ukca_radn2o,              &
         L_ukca_radf11, L_ukca_radf12, L_ukca_radf113,            &
         L_ukca_radf22,                                           &
         L_ukca_intdd, L_ukca_trophet, L_ukca_prescribech4,       &
         L_ukca_set_trace_gases, L_ukca_use_background_aerosol,   &
         L_ukca_primsu, L_ukca_primss,                            &
         L_ukca_primbcoc, L_ukca_primdu, L_ukca_use_2dtop,        &
         L_bcoc_ff, L_bcoc_bf, L_bcoc_bm, L_mode_bhn_on,          &
         L_mode_bln_on, L_ukca_arg_act,                           &
         L_ukca_sfix, i_mode_setup, i_mode_nzts,                  &
         i_mode_bln_param_method, mode_parfrac, dts0, nit,        &
         jvspec_dir, jvspec_file, jvscat_file, phot2d_dir,        &
         strat2d_dir, fastjx_numwl, fastjx_mode,                  &
         fastjx_prescutoff, dir_strat_aer, file_strat_aer,        &
         ukca_MeBrmmr, ukca_MeClmmr, ukca_CH2Br2mmr, ukca_H2mmr,  &
         ukca_N2mmr, ukca_CFC115mmr, ukca_CCl4mmr,                &
         ukca_MeCCl3mmr, ukca_HCFC141bmmr, ukca_HCFC142bmmr,      &
         ukca_H1211mmr, ukca_H1202mmr, ukca_H1301mmr,             &
         ukca_H2402mmr, ukca_COSmmr

      USE ukca_init_mod, ONLY: ukca_init

! module for TKE schemes options
      USE mym_option_mod, ONLY:                                         &
       bdy_tke, tke_dlen, l_tke_dlen_blackadar, tke_cm_mx, tke_cm_fa,   &
       l_my3_improved_closure, my_lowest_pd_surf,                       &
       l_my_condense, l_shcu_buoy, l_adv_turb_field,                    &
       l_my_extra_level, my_z_extra_fact, l_my_prod_adj,                &
       my_prod_adj_fact, tke_levels,                                    &
       l_my_initialize, l_my_ini_zero, l_local_above_tkelvs,            &
       none, my_length, no_pd_surf,                                     &
       high_order_scheme_adv_turb, monotone_scheme_adv_turb,            &
       l_high_adv_turb, l_mono_adv_turb, l_conserv_adv_turb,            &
       l_print_max_tke, l_my_lowest_pd_surf_tqc, my_ini_dbdz_min,       &
       my_z_limit_elb, wb_ng_max, shcu_levels
         

! Module for stochastic physics UM section 35 : SKEB2
      USE stochastic_physics_run_mod, ONLY:                             &
       l_skeb2, l_rp2, run_stochastic,                                  &
      ! Stochastic parameters for ls precip
        m_ci, m_ci_max, m_ci_min, rhcrit_max,  rhcrit_min,              &
      ! Logic control subroutine
        check_run_stochastic


! Module for IAU scheme:
      USE IAU_mod, ONLY : &
          L_IAU,          &
          IAU_nl

! Module for Nudging scheme:
      USE nudging_input_mod, ONLY : Run_Nudging

! Module for GWD scheme:
      USE g_wave_input_mod         

      USE conversions_mod, ONLY: pi
      
      USE q_pos_method_mod, ONLY : &
          q_pos_local

! Modules for reading JULES namelists:
      USE read_jules_namelists_mod, ONLY:                                     &
          read_jules_switches,   read_jules_nvegparm,   read_jules_pftparm,   &
          read_jules_triffid,    read_jules_snow_param, read_jules_soil_param,&
          read_jules_surf_param, read_jules_elevate,    read_jules_rad_param, &
          read_jules_csmin,      read_jules_seed,       read_jules_sigm,      &
          read_urban_switches,   read_urban2t_param

      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook

      USE turb_diff_mod, ONLY: run_diffusion,                                 &
          l_print_shear, l_print_max_wind, l_print_theta1, l_diag_l2norms,    &
          l_diag_l2helm, l_print_lapse, l_diag_print, l_flush6, l_print_pe,   &
          l_diag_print_ops, l_print_div, l_print_wmax, l_print_w,             &
          diag_interval, print_step, first_norm_print, up_diff_scale,         &
          polar_filter_north_lat_limit, polar_filter_step_per_sweep,          &
          polar_filter_lat_limit, polar_filter_south_lat_limit, norm_lev_end, &
          norm_lev_start, w_print_limit, l_diag_wind, dom_n_in, dom_s_in,     &
          dom_e_in, dom_w_in, blocky_in, blockx_in, l_diag_noise, vdiffuv_end,&
          vdiffuv_start, l_vdiff_uv, l_adjust_theta, adjust_lapse_min,        &
          adjust_theta_end, adjust_theta_start, vdiffuv_test, l_filter_incs,  &
          l_filter, ramp_lat_radians, l_diff_auto, diff_coeff_ref,            &
          ref_lat_deg, max_sweeps, vdiffuv_timescale, tar_horizontal,         &
          sponge_power, sponge_ns, q_pos_method, qpos_diag_limit,             &
          l_qpos_diag_pr, q_pos_tracer_method, sponge_ew, top_filt_end,       &
          top_filt_start, l_upper_ramp, top_diff, l_sponge, l_pofil_new   

      USE eg_alpha_mod
      USE eg_alpha_ramp_mod
      USE ref_pro_mod
      USE eg_parameters_mod

      USE ereport_mod, ONLY : ereport
      USE PrintStatus_mod
      USE UM_ParVars
      USE Control_Max_Sizes
      USE science_fixes_mod
      USE blopt8a, ONLY : Limit_ObukhovL
! Module for COSP
      USE cosp_input_mod, ONLY: run_cosp

      USE c_gwave_mod, ONLY: nsigma, amplitude_saturation, stress_saturation,&
                             beta_fix, frac_wl, lambdaz_min, lambdaz_max,    &
                             nsq_neutral, zav_converge, zav_iterate
      
      USE dynamics_input_mod      
      USE dynamics_testing_mod
      USE dynamics_grid_mod, ONLY: l_vatpoles
      
      USE nlstcall_mod, ONLY : run_assim_mode

      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE
!
! Description: Read run-time control information passed as namelists
!   from the UMUI, as required to set up parametrization constants and
!   logical switches needed by physics and dynamics schemes for the
!   Atmosphere model.
! Method:  Sequential read of namelists. Note that defaults are set to
!   missing data indicators. 
!   Namelist variables/declarations which are NOT in modules listed above 
!   are held in cruntimc.h include file.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component: Control Atmos
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
! TYPSIZE start
!   Description:
!     This file contains sizes needed for dynamic allocation of
!   main data arrays within the model. Sizes read in from the user
!   interface via NAMELISTs are passed by /COMMON/. Other control
!   sizes that are fundamental in the definition of data structures
!   are assigned by PARAMETER statements.
!
!   Declarations for the NLSIZES namelist are also held in the module
!   nlsizes_namelist_mod. That module is currently only used by the
!   reconfiguration, while the UM uses this include file.
!
! All sizes
! Not dependent on sub-model
! DATA IN NAMLST#x MEMBER OF THE JOB LIBRARY
! ATMOS START
! Main sizes of fields for each submodel
! Grid-related sizes for ATMOSPHERE submodel.
INTEGER :: ROW_LENGTH           ! No of points per local row
INTEGER :: global_ROW_LENGTH    ! Points per global row
INTEGER :: ROWS                 ! No of local (theta) rows
INTEGER :: global_ROWS          ! No of global (theta) rows
INTEGER :: MODEL_LEVELS         ! No of model levels
INTEGER :: LAND_FIELD           ! No of land points in field
INTEGER :: NTILES               ! No of land surface tiles
INTEGER :: NICE                 ! No. of sea ice thickness categories
INTEGER :: NICE_USE             ! No. of sea ice categories used fully
                                !  in surface exchange and radiation
                                !  (If nice>1 & nice_use=1, categories only 
                                !  partially used in surface exchange)

! Physics-related sizes for ATMOSPHERE submodel
INTEGER :: WET_LEVELS          ! No of moist-levels
INTEGER :: CLOUD_LEVELS        ! No of cloud-levels
INTEGER :: ST_LEVELS           ! No of soil temperature levels
INTEGER :: SM_LEVELS           ! No of soil moisture levels
INTEGER :: BL_LEVELS           ! No of boundary-layer-levels
INTEGER :: OZONE_LEVELS        ! No of ozone-levels
INTEGER :: TPPS_OZONE_LEVELS   ! No of tropopause-ozone-levels
INTEGER :: RIVER_ROWS          ! No of rows for river routing
INTEGER :: RIVER_ROW_LENGTH    ! Row length for river routing
! Dynamics-related sizes for ATMOSPHERE submodel

INTEGER :: TR_LEVELS            ! No of tracer-levels
INTEGER :: TR_VARS              ! No of passive tracers
INTEGER :: TR_LBC_VARS          ! No of tracers in lbcs 
INTEGER :: TR_UKCA              ! No of UKCA tracers
INTEGER :: TR_LBC_UKCA          ! No of UKCA tracer lbcs 

! For Small executables

! Grid related sizes for data structure
! Data structure sizes for ATMOSPHERE submodel
INTEGER :: A_PROG_LOOKUP     ! No of prognostic fields
INTEGER :: A_PROG_LEN        ! Total length of prog fields
INTEGER :: A_LEN_INTHD       ! Length of INTEGER header
INTEGER :: A_LEN_REALHD      ! Length of REAL header
INTEGER :: A_LEN2_LEVDEPC    ! No of LEVEL-dependent arrays
INTEGER :: A_LEN2_ROWDEPC    ! No of ROW-dependent arrays
INTEGER :: A_LEN2_COLDEPC    ! No of COLUMN-dependent arrays
INTEGER :: A_LEN2_FLDDEPC    ! No of FIELD arrays
INTEGER :: A_LEN_EXTCNST     ! No of EXTRA scalar constants
INTEGER :: A_LEN_CFI1        ! Length of compressed fld index 1
INTEGER :: A_LEN_CFI2        ! Length of compressed fld index 2
INTEGER :: A_LEN_CFI3        ! Length of compressed fld index 3
! atmos end

! Data structure sizes for ATMOSPHERE ANCILLARY file control
! routines
INTEGER :: NANCIL_LOOKUPSA  ! Max no of fields to be read

! Data structure sizes for ATMOSPHERE INTERFACE file control
! routines
INTEGER :: N_INTF_A          ! No of atmosphere interface areas
INTEGER :: MAX_INTF_MODEL_LEVELS ! Max no of model levs in all areas
INTEGER :: MAX_LBCROW_LENGTH ! Max no of lbc row length in all areas
INTEGER :: MAX_LBCROWS ! Max no of lbc rows in all areas

!  Data structure sizes for ATMOSPHERE BOUNDARY file control
! routines

! Sizes applicable to all configurations (DUMPS/FIELDSFILE)

INTEGER :: PP_LEN_INTHD   ! Length of PP file integer header
INTEGER :: PP_LEN_REALHD  ! Length of PP file real    header


      ! Grid related sizes for COUPLING between ATMOS and OCEAN
      ! submodels [For MPP, sizes are 'global' values over all
      ! PEs.Also needed for river routing]
      INTEGER:: AOCPL_IMT                ! Ocean rowlength
      INTEGER:: AOCPL_JMT                ! Ocean no. of rows
      INTEGER:: AOCPL_ROW_LENGTH         ! Atmos rowlength
      INTEGER:: AOCPL_P_ROWS             ! Atmos no. of p rows

      COMMON/SIZE_AOCPL/                                                &
        AOCPL_IMT, AOCPL_JMT, AOCPL_ROW_LENGTH, AOCPL_P_ROWS

! Other sizes passed from namelist into common blocks
! Any additions to this common block must be mirrored in nlsizes_namelist_mod.
COMMON/NLSIZES/                                                     &
    ROW_LENGTH,global_ROW_LENGTH,ROWS,global_ROWS,                  &
    LAND_FIELD,MODEL_LEVELS,WET_LEVELS,                             &
    NTILES, NICE, NICE_USE,                                         &
    CLOUD_LEVELS,TR_LEVELS,ST_LEVELS,SM_LEVELS,BL_LEVELS,           &
    OZONE_LEVELS,TPPS_OZONE_LEVELS,TR_VARS,TR_LBC_VARS,             &
    TR_UKCA,TR_LBC_UKCA,RIVER_ROWS,RIVER_ROW_LENGTH,                &
    A_PROG_LOOKUP,A_PROG_LEN,                                       &
    A_LEN_INTHD,A_LEN_REALHD,                                       &
    A_LEN2_LEVDEPC,A_LEN2_ROWDEPC,A_LEN2_COLDEPC,                   &
    A_LEN2_FLDDEPC,A_LEN_EXTCNST,                                   &
    A_LEN_CFI1,A_LEN_CFI2,A_LEN_CFI3,                               &    
    NANCIL_LOOKUPSA,                                                &    
    N_INTF_A, MAX_INTF_MODEL_LEVELS, MAX_LBCROW_LENGTH,             &
    MAX_LBCROWS, PP_LEN_INTHD,PP_LEN_REALHD

!-----------------------------------------------------------------
! data in STASHC#x member of the job library

! Data structure sizes for ATMOSPHERE submodel (config dependent)
INTEGER :: A_LEN2_LOOKUP   ! Total no of fields (incl diags)
INTEGER :: A_LEN_DATA      ! Total no of words of data
INTEGER :: A_LEN_D1        ! Total no of words in atmos D1

! Size of main data array for this configuration

INTEGER :: LEN_TOT             ! Length of D1 array
INTEGER :: N_OBJ_D1_MAX         ! No of objects in D1 array

COMMON/STSIZES/                                                     &
    A_LEN2_LOOKUP,A_LEN_DATA,A_LEN_D1,                              &
    LEN_TOT,N_OBJ_D1_MAX
! global (ie. dump version) of *_LEN_DATA
INTEGER :: global_A_LEN_DATA

COMMON /MPP_STSIZES_extra/ global_A_LEN_DATA
! Sizes of Stash Auxillary Arrays and associated index arrays
! Initialised in UMINDEX and UMINDEX_A/O/W
INTEGER :: LEN_A_IXSTS
INTEGER :: LEN_A_SPSTS

COMMON /DSIZE_STS/                                                  &
    LEN_A_IXSTS, LEN_A_SPSTS
!     The number of land points is computed for each PE
!     before the addressing section. All prognostics on land
!     points in the D1 space are now dimensioned by the local
!     no of land points rather than the global no of land points.

      INTEGER:: global_land_field    !  Global no of land points
      INTEGER:: local_land_field     !  Local no of land points
      COMMON /mpp_landpts/ global_land_field,local_land_field
      ! ----------------------------------------------------------------
      ! extra variables not passed through user interface

      ! fundamental data sizes :
      ! Fundamental parameter  sizes of data structure
      ! Sizes applicable to all configurations (HISTORY FILE)

      ! Length of history file in dump
      INTEGER, PARAMETER :: LEN_DUMPHIST = 0

      ! Sizes applicable to all configurations (DUMPS/FIELDSFILE)
      ! Length of dump fixed header
      INTEGER, PARAMETER :: LEN_FIXHD = 256

      ! Size of a single LOOKUP header
      INTEGER, PARAMETER :: LEN1_LOOKUP  = 64
      INTEGER, PARAMETER :: MPP_LEN1_LOOKUP= 2

      ! Size of compressed LBC LOOKUP (only used internally and
      ! contains just the items which change between each set of LBCs
      INTEGER, PARAMETER :: LEN1_LBC_COMP_LOOKUP = 8

      ! Sizes applicable to all configurations (STASH)
      ! Moved to typstsz.h

      INTEGER:: INTF_LEN2_LEVDEPC !1st dim of interface out lev dep cons
      INTEGER:: INTF_LEN2_ROWDEPC !2nd dim of interface out Row dep cons
      INTEGER:: INTF_LEN2_COLDEPC !2nd dim of interface out Col dep cons
      
      COMMON /DSIZE/                                                    &
        INTF_LEN2_LEVDEPC,INTF_LEN2_ROWDEPC,INTF_LEN2_COLDEPC
      ! sub-model atmosphere   :
      ! Data structure sizes derived from grid size
      INTEGER:: A_LEN1_LEVDEPC ! IN: 1st dim of level  dep const
      INTEGER:: A_LEN1_ROWDEPC ! IN: 1st dim of row    dep const
      INTEGER:: A_LEN1_COLDEPC ! IN: 1st dim of column dep const
      INTEGER:: A_LEN1_FLDDEPC ! IN: 1st dim of field  dep const

      ! Data structure sizes for ATMOSPHERE INTERFACE file control
      ! routines
      INTEGER :: INTF_LOOKUPSA        ! No of interface lookups.
      COMMON /DSIZE_A/                                                  &
        A_LEN1_LEVDEPC,A_LEN1_FLDDEPC,A_LEN1_ROWDEPC,A_LEN1_COLDEPC,    &
        INTF_LOOKUPSA

      ! sub-model atmosphere   : derived sizes
      ! derived from model grid/levels. Arakawa B-grid

                                  ! Size of fields on THETA grid:
      INTEGER :: THETA_FIELD_SIZE     ! IN: with no halos
      INTEGER :: THETA_OFF_SIZE       ! IN: with simple halos
      INTEGER :: THETA_HALO_SIZE      ! IN: with extended halos

                                  ! Size of fields on U grid:
      INTEGER :: U_FIELD_SIZE         ! IN: with no halos
      INTEGER :: U_OFF_SIZE           ! IN: with simple halos
      INTEGER :: U_HALO_SIZE          ! IN: with extended halos

                                  ! Size of fields on V grid
      INTEGER :: V_FIELD_SIZE         ! IN: with no halos
      INTEGER :: V_OFF_SIZE           ! IN: with simple halos
      INTEGER :: V_HALO_SIZE          ! IN: with extended halos

      INTEGER :: N_ROWS               ! IN: No of V-rows
      INTEGER :: N_CCA_LEV            ! IN: No of CCA Levels
      COMMON/DRSIZE_A/                                                  &
        N_ROWS,N_CCA_LEV,THETA_FIELD_SIZE,U_FIELD_SIZE,V_FIELD_SIZE,    &
        THETA_OFF_SIZE,U_OFF_SIZE,V_OFF_SIZE,                           &
        THETA_HALO_SIZE,U_HALO_SIZE,V_HALO_SIZE
      ! boundary updating      : derived values
      ! Variables describing the Atmosphere Lateral Boundary Conditions
      ! Local (per processor) information


! TYPSIZE end
!-----------------Start of TYPCONA------------------------------------

! Constants for routines independent of resolution.
! Requires Use Control_Max_Sizes for MAX_REQ_THPV_LEVS in cconsts.h
! CCONSTS start
! Description:
!    This file contains declarations for derived constants within
!   the atmospheric model. Where necessary PARAMETERS are defined to
!   dimension these constants. All constants are placed in the common
!   block CDERIVED, except hardwired constants, e.g. ETA_SPLIT and LENs.
!   file CMAXSIZE must be included first
!
!   The derived constants are calculated in the routine SETCONA.
!
      ! No of cloud types ie low/med/high
      INTEGER, PARAMETER :: NUM_CLOUD_TYPES = 3

      ! derived constants:
      INTEGER :: LOW_BOT_LEVEL      ! Bottom level of lowest cloud type
      INTEGER :: LOW_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: MED_BOT_LEVEL      ! Bottom   "    "  med      "    "
      INTEGER :: MED_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: HIGH_BOT_LEVEL     ! Bottom   "    "  top      "    "
      INTEGER :: HIGH_TOP_LEVEL     ! Top      "    "   "       "    "

      ! height values to split model levels into l/m/h cloud
      REAL ::    h_split(NUM_CLOUD_TYPES+1)

      LOGICAL :: ELF                ! T if atmosphere model on LAM grid

      ! Constants for dynamics output independent of resolution but
      ! dependent on choice of levels for output.
      REAL :: REQ_THETA_PV_LEVS(MAX_REQ_THPV_LEVS)

      COMMON /CDERIVED/                                                 &
        h_split,LOW_BOT_LEVEL,LOW_TOP_LEVEL,MED_BOT_LEVEL,MED_TOP_LEVEL,&
        HIGH_BOT_LEVEL, HIGH_TOP_LEVEL,ELF,REQ_THETA_PV_LEVS
! CCONSTS end

! typcona.h originally contained constants for the atmosphere.
! Almost all of these constants have moved to a set of modules:
! LEVEL_HEIGHTS_MOD, TRIGNOMETRIC_MOD, DYN_CORIOLIS_MOD, DYN_VAR_RES_MOD,
! DIFF_COEFF_MOD, RAD_MASK_TROP_MOD, ROT_COEFF_MOD

! The following common block did not correspond to constants specified in 
! argcona.h, so it remains here, even though argcona.h has been deleted.
! cderv_trig and CDERIVED in cconsts.h should be moved to modules too.

      ! Trigonometric co-ordinates in radians
      REAL:: Delta_lambda       ! EW (x) grid spacing in radians
      REAL:: Delta_phi          ! NS (y) grid spacing in radians
      REAL:: Base_phi           ! Lat of first theta point in radians
      REAL:: Base_lambda        ! Long of first theta point in radians
      REAL:: lat_rot_NP         ! Real lat of 'pseudo' N pole in radians
      REAL:: long_rot_NP        ! Real long of 'pseudo' N pole in radians

      COMMON/cderv_trig/                                                &
     &  Delta_lambda,Delta_phi,Base_phi,Base_lambda,                    &
     &  lat_rot_NP,long_rot_NP
!-------------End of TYPCONA---------------------------------------
! TYPLNDM
! Formerly integral part of TYPCONA, the variables below have been
! separated from the rest of TYPCONA as they are required by some
! of the Ocean routines in the Ocean-Atmosphere configuration of
! the UM whilest TYPCONA is not.

      ! Primary Arrays
      INTEGER::land_points     ! No. of land points  (can be 0)
      INTEGER::land_ice_points ! Number of land ice points
      INTEGER::soil_points     ! Number of soil points

      ! Do not allow these arrays to have zero size
      INTEGER::land_index    (max(1,land_field)) ! set from land_sea_mask
      INTEGER::land_ice_index(max(1,land_field)) ! Array of land ice points.
      INTEGER::soil_index    (max(1,land_field)) ! Array of soil points.

      ! Gets some sizes transported around the model :
      COMMON /land_soil_dimensions/                                     &
     &  land_points , land_ice_points , soil_points

! TYPLNDM end
!LL  Comdeck: CCONTROL -------------------------------------------------
!LL
!LL  Purpose: COMMON block for top level switches and 2nd level switches
!LL           needed by the top level (C0) and 2nd level routines, but
!LL           not held in the history COMMON block.
!LL
!LLEND ---------------------------------------------------------------

!#include "cntlall.h"
! cntlgen.h was replaced by control/top_level/nlstgen_mod.F90
! #include "cntlgen.h"

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

!
! Subroutine arguments
!   Scalar arguments with intent(in):
!   Array  arguments with intent(in):
!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):
!   Array  arguments with intent(out):

! Local parameters:
      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='Readlsta_4A')

! Local scalars:
      INTEGER                                                           &
     & ErrorStatus      ! Return code : 0 Normal Exit : >0 Error
      CHARACTER(LEN=256)                                                     &
     & CMessage         ! Error message if Errorstatus >0
      INTEGER :: level,j,jj,i,k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('READLSTA_4A',zhook_in,zhook_handle)
      ErrorStatus = 0

! Read atmosphere run consts. control data into COMMON


!     **********    Dynamics defaults   *********


! Dynamics MDI defaults
      ramp_lat_radians = RMDI
      
! Skip all calculations in atmos_physics2 which are identical between
! the first and subsequent iterations
      l_quick_ap2 = .TRUE.

! End of ENDGAME defaults

! Diffusion defaults
      L_filter = .false.
      L_filter_incs = .false.
      L_diff_auto = .false.
      max_sweeps = 8
      ref_lat_deg = 0.0
      diff_coeff_ref = 0.25
      vdiffuv_test = 100.0
      L_vdiff_uv = .false.
      vdiffuv_start = imdi
      vdiffuv_end = imdi
      L_adjust_theta = .false.
      adjust_theta_start = imdi
      adjust_theta_end = imdi
      adjust_lapse_min = 0.0
      vdiffuv_timescale = 1
      L_upper_ramp = .false.
      top_filt_start = imdi
      top_filt_end = imdi
      top_diff = 0.1
      up_diff_scale = 0.5
      L_pofil_new = .false.
      L_sponge = .false.
      sponge_ew = 0
      sponge_ns = 0
      sponge_power = 1

      tar_horizontal = 0


!    QPOS defaults
      q_pos_method           = q_pos_local   ! original with local
      q_pos_tracer_method    = q_pos_local   ! original with local
      l_qpos_diag_pr         = .FALSE.       ! QPOS diagnostics off
      qpos_diag_limit        = -1.0e-8


!     **********    Dynamics defaults   *********

!   Diagnostic printing defaults
      L_print_pe = .false.
      L_flush6 = .false.
      L_diag_print = .false.
      L_diag_print_ops = .false.
      L_print_w = .false.
      L_print_wmax = .false.
      L_print_div = .false.
      L_print_lapse = .false.
      L_print_theta1 = .false.
      L_print_max_wind = .false.
      L_print_shear = .false.
      L_diag_L2norms = .false.
      L_diag_L2helm = .false.
      print_step = 1
      diag_interval = 1
      first_norm_print = 1
      w_print_limit = 1000.0
      norm_lev_start = 1
      norm_lev_end = imdi
      L_print_pe = .false.
      L_diag_wind = .false.
      L_diag_noise = .false.
      dom_w_in = 0
      dom_e_in = 0
      dom_s_in = 0
      dom_n_in = 0
      blockx_in = 0
      blocky_in = 0

!     In common block CDERIVED (in comdeck CCONSTS called in TYPCONA)

! Default settings for h_split (can be overridden by UMUI-supplied
!  values if this is required at a future release, when rmdi defaults
!  should be reinstated):
!  low:middle:high cloud model levels =(1)->(2):(2)->(3):(3)->(4)
      h_split(1) =   111.        ! ICAO 1000mb height (m)
      h_split(2) =  1949.        ! ICAO  800mb height (m)
      h_split(3) =  5574.        ! ICAO  500mb height (m)
      h_split(4) = 13608.        ! ICAO  150mb height (m)

      DO J = 1,MAX_REQ_THPV_LEVS
        REQ_THETA_PV_LEVS(J) = RMDI
      END DO
 
      REWIND(UNIT=5)
! -----------------------------------------
! read in scientific fixes namelist
! controls what fixes are applied 
! each of these fixes are anticipated 
! to become the default code in the future
! -----------------------------------------

      READ (UNIT=5, NML=temp_fixes, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist temp_fixes")
      
      CALL warn_temp_fixes()

! UKCA Sub-model
      READ (UNIT=5, NML=run_ukca, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_UKCA")

! Read in Physics/Dynamics namelists

! Gravity wave drag physics
      READ (UNIT=5, NML=RUN_GWD, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_GWD")
      
! Murk aerosol physics
      READ (UNIT=5, NML=RUN_Murk, IOSTAT=ErrorStatus) 
      CALL check_iostat(errorstatus, "namelist RUN_Murk")

! Convection physics
      READ (UNIT=5, NML=RUN_Convection, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_Convection")
      
! Set other dependent convective switches valid for whole run
! DEPENDS ON: cv_set_dependent_switches
      CALL cv_set_dependent_switches

! Boundary layer physics
      READ (UNIT=5, NML=RUN_BL, IOSTAT=ErrorStatus) 
      CALL check_iostat(errorstatus, "namelist RUN_BL")
      CALL check_run_bl()
      
! River routing
      READ (UNIT=5, NML=RUN_RIVERS, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_RIVERS")

! Large scale precipitation physics 
      READ (UNIT=5, NML=RUN_Precip, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_Precip")

! Radiation physics
      READ (UNIT=5, NML=RUN_Radiation, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_Radiation")
      CALL check_run_radiation()

! Mineral dust modelling
! Initialise sizes in case not set in namelist
      CALL dust_size_dist_initialise
      READ (UNIT=5, NML=RUN_Dust, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_Dust")
      CALL dust_parameters_load
      CALL dust_parameters_check
      
! Large scale cloud physics
      READ (UNIT=5, NML=RUN_Cloud, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_Cloud")
      CALL check_run_cloud()

! Surface type characteristics
      READ (UNIT=5, NML=RUN_BLVEG, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_BLVEG")
      
! Energy correction physics
      READ (UNIT=5, NML=RUN_Eng_Corr, IOSTAT=ErrorStatus) 
      CALL check_iostat(errorstatus, "namelist RUN_Eng_Corr")

! Stochastic physics
      READ (UNIT=5, NML=RUN_Stochastic, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_Stochastic")
      CALL check_run_stochastic()
     
! Aerosol Modelling
      READ (UNIT=5, NML=RUN_Aerosol, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_Aerosol")
      
! Check UKCA logicals are consistent and 
! set internal UKCA values from UKCA namelists
! Needs to be after CLASSIC aerosol namelist is read as 
! well as after UKCA
      CALL ukca_init()
      
! UM nudging
      READ (UNIT=5, NML=RUN_Nudging, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_Nudging")
      
! Generalised integration and GCR dynamics
      READ (UNIT=5, NML=RUN_Dyn, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_Dyn")
      CALL check_run_dyn()

! Generalised integration and GCR dynamics
      READ (UNIT=5, NML=RUN_Dyntest, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_Dyntest")
      CALL check_run_dyntest()

! Semi-Lagrangian advection dynamics
      READ (UNIT=5, NML=RUN_SL, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_SL")
      CALL check_run_sl()
 
! Diffusion, divergence damping and filtering dynamics
      READ (UNIT=5, NML=RUN_Diffusion, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_Diffusion")
      
! Call to COSP
      READ (UNIT=5, NML=RUN_COSP, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RUN_COSP")

! Diagnostic double call to radiation
      READ (UNIT=5, NML=RADFCDIA, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist RADFCDIA")

!------------------------------------------------------- 
! Check that if L_CCRAD is selected certain other
! switches are not set so that they conflict.
!------------------------------------------------------- 
! Error capture 
!------------------------------------------------------- 
                 
      If (L_ccrad) Then 

        If (.NOT. L_3D_CCA) Then 
          ErrorStatus = 100 
          CMessage    = '**ERROR**: CCRad is not yet available without'// &
                                  ' the anvil scheme (L_3D_CCA = .True.)'
 
          CALL ereport(RoutineName, ErrorStatus, CMessage) 
        End If

        If (l_fix_udfactor) Then 
          ErrorStatus = 101 
          CMessage    = '**ERROR**: L_CCRad and l_fix_udfactor'//         &
                                  ' should not be both set to true.'
 
          CALL ereport(RoutineName, ErrorStatus, CMessage) 
        End If

        If (l_pc2_diag_sh) Then 
          ErrorStatus = 102 
          CMessage    = '**ERROR**: L_CCRad and l_pc2_diag_sh'//          &
                                  ' should not be both set to true.'
 
          CALL ereport(RoutineName, ErrorStatus, CMessage) 
        End If

      End If       ! l_ccrad 

      IF (l_pc2) THEN 
        ! So that PC2 does NOT affect Section 5 diagnostics
        l_dcpl_cld4pc2=.TRUE.
      END IF

!------------------------------------------------------- 
! End error capture 
!-------------------------------------------------------

      IF (L_Backwards) L_Physics = .false.
      
      IF (l_endgame) THEN
          l_vatpoles=.TRUE.  ! must be the case for the 4A scheme
          NumCycles = 2       ! Only option for EG so set it here
      END IF    

      IF(PrintStatus >= PrStatus_Normal) THEN

        WRITE(6,'(/,A,A,A)') '******************** ',RoutineName,         &
                    ': Atmosphere run-time constants *******************'
                    
        WRITE(6,temp_fixes)          
        WRITE(6,RUN_Cloud)
        WRITE(6,RUN_BLVEG)
        WRITE(6,RUN_BL)
        WRITE(6,RUN_Eng_Corr)
        WRITE(6,RUN_Murk)
        WRITE(6,RUN_Precip)
        WRITE(6,RUN_Convection)
        WRITE(6,RUN_Stochastic)
        WRITE(6,RUN_Radiation)
        WRITE(6,RUN_GWD)
        WRITE(6,RUN_Aerosol)
        WRITE(6,RUN_Dust)
        WRITE(6,RUN_UKCA)
        WRITE(6,RUN_Nudging)
        WRITE(6,RUN_Dyn)
        WRITE(6,RUN_Dyntest)
        WRITE(6,RUN_SL)
        WRITE(6,RUN_Diffusion)
        WRITE(6,RUN_RIVERS)
        WRITE(6,RUN_COSP)
        WRITE(6,RADFCDIA)

      END IF ! PrintStatus



! Convert from degrees to radians
      polar_filter_north_lat_limit= polar_filter_north_lat_limit*pi/180.
      polar_filter_south_lat_limit= polar_filter_south_lat_limit*pi/180.
      polar_filter_lat_limit      = polar_filter_lat_limit      *pi/180.
      polar_filter_step_per_sweep = polar_filter_step_per_sweep *pi/180.
      ramp_lat_radians = ramp_lat_radians *pi/180.



!  Multiply input req_theta_pv_levs by 100. to convert mb to pascals.
      DO LEVEL = 1,MAX_REQ_THPV_LEVS
        REQ_THETA_PV_LEVS(LEVEL) = REQ_THETA_PV_LEVS(LEVEL)*100.
      END DO

! Set radiation aerosol switches
      IF (cusack_aero==2 .OR.  cusack_aero==3    ) l_climat_aerosol    = .TRUE.
      IF (cusack_aero==3 .AND. cusack_aero_hgt==1) l_clim_aero_hgt     = .TRUE.
      IF (cusack_aero==2)                          l_HadGEM1_clim_aero = .TRUE.

!
!     Read controlling information for version 3C/Z of the SW or LW
!     radiation schemes.
!
      READ (UNIT=5, NML=R2SWNCAL, IOSTAT=ErrorStatus)   
      CALL check_iostat(errorstatus, "namelist R2SWNCAL")
      
      IF (n_swcall > npd_swcall) THEN
        ErrorStatus = 100 
        CMessage    = '**ERROR**: Too many calls to SW radiation:'//  &
                      ' n_swcall > npd_swcall. Increase npd_swcall'// &
                      ' in max_calls.F90 and recompile.'
 
        CALL ereport(RoutineName, ErrorStatus, CMessage) 
      END IF

! Options for the shortwave
      CALL sw_input
        
      READ (UNIT=5, NML=R2LWNCAL, IOSTAT=ErrorStatus)   
      CALL check_iostat(errorstatus, "namelist R2LWNCAL")
      IF (n_lwcall > npd_lwcall) THEN
        ErrorStatus = 100 
        CMessage    = '**ERROR**: Too many calls to LW radiation:'//  &
                      ' n_lwcall > npd_lwcall. Increase npd_lwcall'// &
                      ' in max_calls.F90 and recompile.'
 
        CALL ereport(RoutineName, ErrorStatus, CMessage) 
      END IF

!       Options for the longwave 
      CALL lw_input
          
      READ (UNIT=5, NML=clmchfcg, IOSTAT=ErrorStatus)   
      CALL check_iostat(errorstatus, "namelist CLMCHFCG")
      
      IF ( L_clmchfcg ) THEN
         !  Convert rates from percent to multiplicative factors:
         DO j = 1, nclmfcgs
!          !  This is a null loop, as it should be, if clim_fcg_nyears=0
           DO jj = 1, clim_fcg_nyears(j)
             IF ( clim_fcg_rates(jj,j)  >   -100. ) THEN
                clim_fcg_rates(jj,j) = 1. + 0.01 * clim_fcg_rates(jj,j)
             END IF
           END DO
         END DO
      ELSE
!        ! If the namelist is not to be read, set number of designated
!        !   years to zero for all possible forcings, as this may be
!        !   used to test if this system is being used for each forcing.
         DO J = 1, nclmfcgs
           clim_fcg_nyears(J) = 0
         END DO
      END IF !  L_clmchfcg

      CALL coradoca_defaults

! Read the JULES namelists
      REWIND( 5 )

! Unit number is passed as argument
      CALL read_jules_switches( 5 )

      CALL read_urban_switches( 5 )

      CALL read_jules_snow_param( 5 )

      CALL read_jules_soil_param( 5 )

      CALL read_jules_nvegparm( 5 )

      CALL read_jules_pftparm( 5 )

      CALL read_jules_triffid( 5 )

      CALL read_jules_surf_param( 5 )

      CALL read_jules_elevate( 5 )

      CALL read_jules_rad_param( 5 )

      CALL read_jules_csmin( 5 )

      CALL read_jules_seed( 5 )

      CALL read_jules_sigm( 5 )

      CALL read_urban2t_param( 5 )
      
! In the coupled model case, set alpham/c to ssalpham/c and dtice to ssdtice 
! if l_ssice_albedo == T. Note, ssalpham/c, ssdtice accessed via rad_input_mod
! Do this after JULES namelist reads, as l_ssice_albedo is in JULES_SWITCHES
      IF (l_oasis .AND. l_ssice_albedo ) THEN
        alpham = ssalpham
        alphac = ssalphac
        dtice  = ssdtice
      END IF

! ROSE changes for BL_option_mod - logic control moved to after JULES namelist
! reads as l_flake_model now in jules_switches
      IF (l_flake_model) THEN
        cor_mo_iter = Limit_ObukhovL
        CMessage =&
       ': WARNING, cor_mo_iter set to Limit_ObukhovL since l_flake_model on'
        IF(PrintStatus >= PrStatus_Normal) THEN
          ErrorStatus = -100
          CALL ereport(RoutineName, ErrorStatus, CMessage)
        END IF
      END IF

      IF (atmos_sr(3) == '1A') THEN
        ishear_bl = bl_option_off
        CMessage =&
       ': WARNING, ishear_bl set to off since boundary layer version 1A chosen'
        IF(PrintStatus >= PrStatus_Normal) THEN
          ErrorStatus = -100
          CALL ereport(RoutineName, ErrorStatus, CMessage)
        END IF
      END IF

      k = 0
      DO i = 1, alpha_cd_batches
        DO j = 1, alpha_cd_items(i)
          k = k + 1
          alpha_Cd(k) = alpha_cd_vals(i)
        END DO
      END DO
! ROSE changes for BL end

      REWIND (UNIT=5)

      ! Read IAU namelist:
      READ (UNIT=5, NML=IAU_nl, IOSTAT=ErrorStatus)
      IF (ErrorStatus > 0) THEN
        CMessage(1:)   = 'Error reading IAU namelist IAU_nl.'
        CMessage(81:)  = 'Note that the namelist was heavily revised at vn7.6.'
        CMessage(161:) = 'See the UMUI help panel, or section A.2 in '// &
                         'UMDP 31 for details.'
        CALL EReport (RoutineName, ErrorStatus, CMessage)
      END IF

      ! Turn off IAU if in FASTRUN mode:
      IF (run_assim_mode == "FastRun   ") THEN
        IF (L_IAU) THEN
          run_assim_mode = "NoIAU     "
          L_IAU = .FALSE.
        ELSE
          run_assim_mode = "None      "
        END IF
      END IF

      REWIND (UNIT=5)

! Check error condition
      IF(ErrorStatus >  0) THEN

         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      END IF

      IF (lhook) CALL dr_hook('READLSTA_4A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Readlsta_4A
