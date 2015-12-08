! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! scm_shell is the main calling program for the Single Column Model.
! It sets up the vertical level information read in from the UMUI and
! passes the information down to Scm_Main which then performs the run.
!
! Program scm_shell
!=====================================================================
!                     SCM
!           Single Column Unified Model
!                  Master Deck
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
!=====================================================================
PROGRAM scm_shell

  USE um_input_control_mod, ONLY:                                          &
    nlstcatm, l_mr_physics2,  lcal360

  ! Dependencies for interfacing with C routines
  USE io_dependencies
  
  USE dynamics_input_mod
  USE dynamics_testing_mod

  ! SCM Modules
  !---------------------------------------------------------------------------
  USE scm_utils
  USE scm_cntl_mod
  USE s_main_force, ONLY: netcdf_file, obs_t0_prf, source
  USE s_nc_obs, ONLY: read_nc_obs
  USE global_SCMop, ONLY: incdf
  USE netcdf

  ! Physics modules
  !---------------------------------------------------------------------------
  ! Convection modules
  !---------------------------------------------------------------------------
  USE cv_run_mod                         ! Access to all variables
  USE cv_param_mod                       ! Access to all variables

  ! Module for RUN_Precip namelist
  USE mphys_inputs_mod, ONLY: RUN_precip

  ! Module for RUN_Cloud namelist
  USE cloud_inputs_mod, ONLY: rhcrit, pc2_falliceshear_method,      &
   cloud_fraction_method, ice_fraction_method, i_fixbug_pc2_checks, &
   i_pc2_conv_coupling, i_pc2_erosion_method, L_eacf,               &
   l_ensure_min_in_cloud_qcf, l_fixbug_pc2_qcl_incr,                &
   l_fixbug_pc2_mixph, l_micro_eros, overlap_ice_liquid,            &
   ctt_weight, t_weight, qsat_fixed, sub_cld, dbsdtbs_turb_0,       &
   starticeTKelvin, alliceTdegC, cff_spread_rate, ice_width,        &
   l_pc2, l_rhcpt, RUN_cloud, check_run_cloud

  ! Module for RUN_LAND, RUN_BLVEG and RUN_PFT namelists
  USE land_surf_mod
  USE switches, ONLY: can_rad_mod, l_flake_model

  ! Module for RUN_Murk namelist
  USE murk_inputs_mod, ONLY: RUN_murk

  ! Module for RUN_Rivers namelist
  USE river_inputs_mod, ONLY: RUN_Rivers

  ! Module for RUN_Eng_Corr namelist
  USE eng_corr_inputs_mod, ONLY: RUN_Eng_Corr

  ! Boundary layer modules
  USE bl_option_mod, ONLY: run_bl, cor_mo_iter, ishear_bl,          &
        alpha_cd_batches, alpha_cd_items, alpha_cd_vals, alpha_cd,  &
        tke_levels, shcu_levels, buddy_sea, off, check_run_bl
  USE blopt8a, ONLY : Limit_ObukhovL

  ! module for UKCA options
  USE ukca_option_mod, ONLY: run_ukca,                            &
         l_ukca, l_ukca_aie1, l_ukca_aie2,                        &
         i_ukca_chem,                                             & 
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

  ! GWD modules
  USE g_wave_input_mod

  ! JULES
  USE switches, ONLY:                                                         &
     l_360, l_aggregate

! Modules for reading JULES namelists:
  USE read_jules_namelists_mod, ONLY:                                         &
      read_jules_switches,   read_jules_nvegparm,   read_jules_pftparm,       &
      read_jules_triffid,    read_jules_snow_param, read_jules_soil_param,    &
      read_jules_surf_param, read_jules_elevate,    read_jules_rad_param,     &
      read_jules_csmin,      read_jules_seed,       read_jules_sigm,          &
      read_urban_switches,   read_urban2t_param
  USE init_from_jules_namelists_mod, ONLY: init_from_jules_namelists
  USE nstypes,                       ONLY: npft, nnvg
  USE calc_ntiles_mod,               ONLY: calc_ntiles

  ! Radiation modules
  USE rad_input_mod
  USE sw_rad_input_mod
  USE lw_rad_input_mod
  
  USE missing_data_mod, ONLY: imdi

  USE nlstcall_mod

  ! UM Modules
  !---------------------------------------------------------------------------
  USE atm_fields_bounds_mod, ONLY: atm_fields_bounds_init
  USE run_aerosol_mod
  USE turb_diff_mod, ONLY: run_diffusion
  USE ereport_mod, ONLY : ereport
  USE Control_Max_Sizes
  USE um_parvars
  USE filenamelength_mod, ONLY: filenamelength
  USE UM_Config, ONLY : &
       appInit,         &
       exe_scm
  USE PrintStatus_mod
  USE Submodel_Mod, ONLY: n_internal_model

  USE chsunits_mod, ONLY : nunits

  IMPLICIT NONE

  INTEGER :: me_gc
  INTEGER :: nproc_gc

!include file: s_dims.h
! Description:
!   This include file defines variables from the UM NLSIZES / STSHCOMP
!   namelists which are used in the SCM.
!
!
! Declarations:

  INTEGER ::       &
    rows           &! No of rows
  , row_length     &! Row length  
  , land_field     &! No of land points in field
  , ntiles         &! No of land surface tiles
  , nice           &! No of sea-ice thickness categories
  , nice_use        ! No of sea-ice thickness categories used fully

  INTEGER ::       &
    model_levels   &! No of model levels
  , wet_levels     &! No of moist-levels
  , cloud_levels   &! No of cloud-levels
  , st_levels      &! No of soil temperature levels
  , sm_levels      &! No of soil moisture levels
  , bl_levels      &! No of boundary-layer-levels
  , ozone_levels   &! No of ozone-levels
  , tr_levels      &! No of tracer-levels
  , tr_vars        &! No of passive tracers
  , tr_ukca         ! No of UKCA tracers

! Max. no. of STASH sections  per internal model (44 in practice)
! NOTE:  Not used for STASH purposes in the SCM.
! Copied from version.h
  INTEGER,PARAMETER :: nsectp = 99

! String array to hold version choices for various physics schemes.
! Only radiation choices (indices 1 and 2) are used in SCM.
! Copied from model.h
  CHARACTER(LEN=2) :: atmos_sr(0:nsectp)

!---------------------------------------------------------------------
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

! Local Variables

  INTEGER, PARAMETER :: &
    ntrop       = 20    &! Max number of levels in the troposphere
                         ! STATS forcing
  , sec_day     = 86400 &
  , nsprog      = 10    &! No. of single level prognostics
  , ntab        = 32     ! Dimension of array used in random
                         ! generator (Do not change this value as it
                         ! is hard coded into the S_RANDOM deck)

  INTEGER, PARAMETER :: &
    co2_dim_len = 1     &! Length of a CO2 field row.
  , co2_dim_row = 1      ! Number of CO2 field rows.

  REAL ::               &
    dummy

  INTEGER ::            &
    Istatus             &
  , Icode               &
  , first_blank         &
  , level               &
  , j,i,k

  LOGICAL :: l_ts_log ! Option for timestep information

  CHARACTER (LEN=filenamelength) ::      &
    dname               &! Directory where the sizes file is kept
  , filename            &! Sizes filename to read in basic model dimensions
  , vert_lev             ! Vertical level file

  CHARACTER(LEN=100) ::  dummy_env

  CHARACTER(LEN=256) :: Cmessage  ! Error message if ErrorStatus > 0
  CHARACTER (Len=*), PARAMETER :: RoutineName = 'scm_shell'

!---------------------------------------------------------------------
! Variables read in which refer to the SCM forcing namelist
!---------------------------------------------------------------------
  INTEGER :: land_points      = 0    ! Default number of land points
  INTEGER :: nfor             = IMDI ! Number terms for observational
                                     ! forcing, requires setting in
                                     ! namelist
  INTEGER :: model_levels_nml = IMDI ! Number of model levels
                                     ! specified in supplied
                                     ! namelist. Must be set in namelist.

  ! Variables for NC_obs namelist and to read NetCDF_file
  LOGICAL        :: l_netcdf_obs = .FALSE.

  INTEGER(incdf) :: status
  INTEGER(incdf) :: ncid
  INTEGER(incdf) :: time_dimid, id_lsm
  INTEGER(incdf) :: nt_in

  REAL, ALLOCATABLE :: rdum1d(:)

  ! Dr Hook
  !==============================
  REAL(KIND=jprb) :: zhook_handle

  NAMELIST/VERT_RES/ vert_lev
  NAMELIST/CNTLSCM/                                                    &
          nfor, model_levels_nml, l_ts_log, land_points, l_netcdf_obs

! Initialise GCOM
  CALL gc_init(dummy_env, me_gc, nproc_gc)

  IF (lhook) CALL dr_hook('SCM_SHELL',zhook_in,zhook_handle)
  CALL appInit(exe_scm)


!=====================================================================
! Assign variables normally for MPP configs
!=====================================================================
  halo_i = 0   ! Halo in i
  halo_j = 0   ! Halo in j
  offx   = 0   ! Small halo in i.
  offy   = 0   ! Small halo in j.

!=====================================================================
!     First read in directory where UMUI files
!=====================================================================

  CALL fort_get_env('JOBDIR', 6, dname, filenamelength, istatus)

  IF (istatus /= 0) THEN
    icode    = -506
    cmessage = ' Environment variable $JOBDIR '//                         &
               'not set, use current directory'
    CALL ereport(routinename, icode, cmessage)
    dname = '.'
  END IF

  first_blank = LEN_TRIM(dname)+1
  dname       = TRIM(dname)


!=====================================================================
! Read in model dimensions and sections from SIZES file
!=====================================================================

  filename = dname(1:first_blank-1)//'/SIZES'

  ! DEPENDS ON: read_um_nml
  CALL read_um_nml ( filename, rows, row_length                           &
    , model_levels, wet_levels, cloud_levels, bl_levels, ozone_levels     &
    , st_levels, sm_levels, tr_levels, tr_vars, tr_ukca, nice             &
    , atmos_sr )

   nice = 1              ! No. of sea ice categories
   nice_use  = 1         ! No. of sea ice categories used fully in
                         !  surface exchange

!=====================================================================
! Grid definitions for NewDynamics/EndGame
!=====================================================================
  CALL atm_fields_bounds_init                                             &
     ( offx, offy, halo_i, halo_j, row_length, rows, rows                 &
     , model_levels, wet_levels, tr_levels, bl_levels, ozone_levels )


!=====================================================================
!     Now read in the name of the file containing the vertical level
!     information
!=====================================================================
  filename = dname(1:first_blank-1)//'/SCM_SET'

  OPEN(10, File=Filename, Iostat=IstatuS, Status='old')
  IF (Istatus  /=  0) THEN
    Icode = 500
    WRITE(Cmessage,*) " Error opening " //TRIM(ADJUSTL(Filename))//       &
                      " on unit 10"
    
    CALL ereport (RoutineName, Icode, Cmessage)
  END IF

  READ(10,VERT_RES)
  READ(10,SCM_CNTL)

  CLOSE(10)

!=====================================================================
!     Read SCM runtime namelist CNTLSCM
!=====================================================================

  OPEN(10, File=TRIM(ADJUSTL(scm_nml)), iostat=istatus, status='old')

  IF (istatus /= 0) THEN

    WRITE(6,*)  " Error opening " //TRIM(ADJUSTL(scm_nml))
    WRITE(6,*)  " Checking for namelist.scm in current directory"

    scm_nml = './namelist.scm'
    OPEN(10, File=TRIM(ADJUSTL(scm_nml)), iostat=istatus, status='old')

    IF (istatus /= 0) THEN
      icode = 500
      WRITE(cmessage,*)  " Error opening namelist.scm"
      CALL ereport (routinename, icode, cmessage)
    END IF

  END IF

  READ(10,CNTLSCM)

  CLOSE(10)

  ! If using external netcdf file, read the nc_obs namelist
  ! in scm_nml for information the netcdf file.
  IF (l_netcdf_obs) THEN
    CALL read_nc_obs
  END IF

!=====================================================================
!     Read in the namelists from the CNTLALL file
!=====================================================================
  filename = dname(1:first_blank-1)//'/CNTLALL'

  OPEN(10, File=Filename, Iostat=Istatus, Status='old')
  IF (Istatus  /=  0) THEN
    Icode = 500
    WRITE(Cmessage,*) " Error opening " //TRIM(ADJUSTL(Filename))//       &
                      " on unit 10"
    CALL ereport (RoutineName, Icode, Cmessage)
  END IF

  READ(10,NLSTCALL)

  CLOSE(10)

!=====================================================================
!     Read in namelists from the SHARED file
!=====================================================================
  filename = dname(1:first_blank-1)//'/SHARED'

  OPEN(10, File=Filename, Iostat=Istatus, Status='old')
  IF (Istatus  /=  0) THEN
    Icode = 500
    WRITE(Cmessage,*) " Error opening " //TRIM(ADJUSTL(Filename))//       &
                      " on unit 10"
    CALL ereport (RoutineName, Icode, Cmessage)
  END IF

! Read a shared JULES namelist
! Unit number is passed as argument
  CALL read_jules_switches( 10 )
  CALL read_urban_switches( 10 )

  READ(10,NLSTCATM)

  l_360      = lcal360        ! Overwrite the JULES module values

  READ(10,RUN_UKCA)
  READ(10,RUN_GWD)
! Murk not yet used by SCM - read in for consistency with UM
  READ(10,RUN_Murk)
  READ(10,RUN_Convection)
! Set other dependent convective switches valid for whole run
! DEPENDS ON: cv_set_dependent_switches
  CALL cv_set_dependent_switches

  READ(10,RUN_BL)
  CALL check_run_bl()

! ROSE changes for BL
  IF (l_flake_model) THEN
    cor_mo_iter = Limit_ObukhovL
    cmessage =&
   ': WARNING, cor_mo_iter set to Limit_ObukhovL since l_flake_model on'
    IF(PrintStatus >= PrStatus_Normal) THEN
      icode = -100
      CALL ereport(RoutineName, icode, CMessage)
    END IF
  END IF
  IF (atmos_sr(3) == '1A') THEN
    ishear_bl = off
    cmessage =&
   ': WARNING, ishear_bl set to off since boundary layer version 1A chosen'
    IF(PrintStatus >= PrStatus_Normal) THEN
      icode = -100
      CALL ereport(RoutineName, icode, CMessage)
    END IF
  END IF

  k = 0
  DO i = 1, alpha_cd_batches
    DO j = 1, alpha_Cd_items(i)
      k = k + 1
      alpha_Cd(k) = alpha_Cd_vals(i)
    END DO
  END DO
! ROSE changes for BL end

! River routing not yet used by SCM - read in for consistency with UM
  READ(10,RUN_Rivers)
  READ(10,RUN_Precip)
  READ(10,RUN_Radiation)
  CALL check_run_radiation()
  READ(10,RUN_Cloud)
  CALL check_run_cloud()

! Read shared JULES namelists
! Unit number is passed as argument
  CALL read_jules_snow_param( 10 )
  CALL read_jules_soil_param( 10 )

  CLOSE(10)

!=====================================================================
!     Calculate the value to use for NTILES
!=====================================================================

  ntiles=9
  CALL calc_ntiles(l_aggregate,npft,nnvg,ntiles)

!=====================================================================
!     Read in the namelists from the CNTLATM file
!=====================================================================
  filename = dname(1:first_blank-1)//'/CNTLATM'

  OPEN(10, File=Filename, Iostat=Istatus, Status='old')
  IF (Istatus  /=  0) THEN
    Icode = 500
    WRITE(Cmessage,*) " Error opening " //TRIM(ADJUSTL(Filename))//       &
                      " on unit 10"
    CALL ereport (RoutineName, Icode, Cmessage)
  END IF

! Default value of CAN_RAD_MOD which is defined in RUN_BLVEG
  can_rad_mod = 1

  READ(10,RUN_BLVEG)

! A negative "tke_levels" should be equal to bl_levels.
  IF (tke_levels < 0) THEN
    tke_levels = bl_levels
  END IF
! A negative "shcu_levels" should be equal to tke_levels.
  IF (shcu_levels < 0) THEN
    shcu_levels = tke_levels
  END IF

! Eng_Corr not yet used by SCM - read in for consistency with UM
  READ(10,RUN_Eng_Corr)
  READ(10,RUN_Aerosol)
  READ(10,RUN_Dyn)
  READ(10,RUN_DynTest)
  READ(10,RUN_Diffusion)

! Read in Namelists R2SWCLNL and R2LWCLNL and transfer data to
! the structure SW_CONTROL and LW_CONTROL. This part of the
! code (under 3C and 3Z) uses modules to pass arguments around.

  READ(10, R2SWNCAL)

! Set radiation aerosol switches
      IF (cusack_aero==2 .OR.  cusack_aero==3    ) l_climat_aerosol    = .TRUE.
      IF (cusack_aero==3 .AND. cusack_aero_hgt==1) l_clim_aero_hgt     = .TRUE.
      IF (cusack_aero==2)                          l_HadGEM1_clim_aero = .TRUE.
      IF (cusack_aero_hgt==2) aero_bl_levels = bl_levels

! Options for the shortwave
  CALL sw_input


  READ(10, R2LWNCAL)
! Options for the longwave 
  CALL lw_input
  
! Read the JULES namelists
! Unit number is passed as argument
  CALL read_jules_nvegparm( 10 )

  CALL read_jules_pftparm( 10 )

  CALL read_jules_triffid( 10 )

  CALL read_jules_surf_param( 10 )

  CALL read_jules_elevate( 10 )

  CALL read_jules_rad_param( 10 )

  CALL read_jules_csmin( 10 )

  CALL read_jules_seed( 10 )

  CALL read_jules_sigm( 10 )

  CALL read_urban2t_param( 10 )

! Initialise after everything has been read
  CALL init_from_jules_namelists( land_points, ntiles, sm_levels,   &
                                  nice ,nice_use )

  CLOSE(10)

  n_internal_model = 1


!---------------------------------------------------------------------
! Set dimensions held in scm module

  nobs    = nfor
  nmod_lv = model_levels
  nbl_lv  = bl_levels
  nwet_lv = wet_levels
  ntile   = ntiles
  o3_lv   = ozone_levels
  ntr_lv  = tr_levels
  ntr_var = tr_vars
  nlnd    = land_points
  rw      = rows
  rw_lng  = row_length
  obs_t0  = obs_t0_prf
  st_lv   = st_levels
  sm_lv   = sm_levels

  nml_nmod_lv = model_levels_nml

  IF (nfor == IMDI) nc_obsfor = .TRUE.

!-------------------------------------------------------
! Error capture
!-------------------------------------------------------

!-----------------------------------------------------------------------------
! Check settings
!-----------------------------------------------------------------------------

  IF (   l_use_orog_corr  .OR. l_use_grad_corr                   &
    .OR. l_use_skyview    .OR. l_orog_unfilt    ) THEN

    WRITE (*,'(3(TR1,A52/),TR1,A52)')                            &
      '====================================================',    &
      '| Invalid control options for SCM. Altering the    |',    &
      '| following NLSTCATM namelist variables:           |',    &
      '|                                                  |'

    IF (l_use_orog_corr) THEN
      WRITE (*,'(TR1,A52)')                                      &
      '|  l_use_orog_corr = .FALSE.                       |'
      l_use_orog_corr=.FALSE.
    END IF
    IF (l_use_grad_corr) THEN
      WRITE (*,'(TR1,A52)')                                      &
      '|  l_use_grad_corr = .FALSE.                       |'
      l_use_grad_corr=.FALSE.
    END IF
    IF (l_use_skyview) THEN
      WRITE (*,'(TR1,A52)')                                      &
      '|  l_use_skyview = .FALSE.                         |'
      l_use_skyview=.FALSE.
    END IF
    IF (l_orog_unfilt) THEN
      WRITE (*,'(TR1,A52)')                                      &
      '|  l_orog_unfilt = .FALSE.                         |'
      l_orog_unfilt=.FALSE.
    END IF

    WRITE (*,'(2(TR1,A52/))')                                    &
      '|                                                  |',    &
      '===================================================='

  END IF ! Test for invalid NLSTCATM options

  IF (l_netcdf_obs) THEN

    status = 0
    IF (source == 2) THEN
      status = nf90_open(TRIM(netcdf_file), nf90_noWrite, ncid)

      IF (status /= nf90_noerr) THEN
        PRINT*,'Error opening netcdf file'
      END IF
      status = nf90_inq_dimid (ncid, 'time', time_dimid)
      status = nf90_inq_dimid (ncid, 'time', time_dimid)
      status = nf90_inquire_dimension (ncid, time_dimid, len=nt_in)
      nfor   = nt_in - obs_t0 + 1
      nobs   = nfor

      ALLOCATE( rdum1d (nobs))

      status = nf90_inq_varid (ncid, 'lsm',  id_lsm)
      status = nf90_get_var   (ncid, id_lsm, rdum1d)

      IF (rdum1d(obs_t0) > 0.5) THEN
        land_points = 1
      ELSE
        land_points = 0
      END IF

      nlnd = land_points

      DEALLOCATE(rdum1d)

      status = nf90_close(ncid)
    END IF

  ELSE

    IF (model_levels_nml == IMDI) THEN
      icode=501
      WRITE(6,*) " "
      WRITE(6,*) "============================================="
      WRITE(6,*) " Number of model levels (model_levels_nml) in"
      WRITE(6,*) " SCM namelist (&CNTLSCM) has not been set.   "
      WRITE(6,*) "============================================="
      WRITE(6,*) " "
      WRITE(6,*) "Run ABORTED"
      WRITE(6,*) " "

      WRITE(cmessage,*) ' Variable MODEL_LEVELS_NML has not been set'
      CALL ereport(routinename, icode, cmessage)
    END IF

    IF (land_points > row_length*rows) THEN
      icode=502
      WRITE(6,*) " "
      WRITE(6,*) "============================================="
      WRITE(6,*) " Specified number of land points greater than"
      WRITE(6,*) " row_length*rows.                            "
      WRITE(6,*) "============================================="
      WRITE(6,*) " "
      WRITE(6,*) "Run ABORTED"
      WRITE(6,*) " "

      WRITE(cmessage,*) ' Too many land points specified'
      CALL ereport(routinename, icode, cmessage)
    END IF

  END IF ! l_netcdf_obs

  IF (buddy_sea /= off) THEN

    WRITE (*,'(8(TR1,A52/))')                                   &
      '===================================================='    &
    , '| Invalid Boundary Layer options for SCM           |'    &
    , '| Altering the following RUN_BL namelist variable  |'    &
    , '|                                                  |'    &
    , '|  Buddy_sea                                       |'    &
    , '|                                                  |'    &
    , '===================================================='

    buddy_sea = off

  END IF ! Test for invalid RUN_BL options

  !---------------------------------------------------------------
  ! Check that number of Wet and Cloud levels do not exceed the
  ! number of model levels
  !---------------------------------------------------------------
  IF (wet_levels > model_levels) THEN

    WRITE (*,'(4(TR1,A52/))')                                   &
      '===================================================='    &
    , '| Warning : Wet_levels > Model_levels              |'    &
    , '| Setting : Wet_levels = Model_levels              |'    &
    , '===================================================='

    wet_levels = model_levels

  END IF

  IF (cloud_levels > wet_levels) THEN
    WRITE (*,'(4(TR1,A52/))')                                   &
      '===================================================='    &
    , '| Warning : Cloud_levels > Wet_levels              |'    &
    , '| Setting : Cloud_levels = Wet_levels              |'    &
    , '===================================================='

    cloud_levels = wet_levels
  END IF

  If ((tr_levels /= model_levels) .AND. &
      (tr_levels /= 0)            .AND. &
      (tr_levels /= 1)) THEN 

    ! This test is enforced because at present, tracers
    ! from different areas of the UM are treated inconsistently.
    ! This causes issues for tr_levels set to anything other 
    ! that model_levels
    icode    = -99
    cmessage = ' Setting tracer levels = model levels'
    CALL ereport(routinename, icode, cmessage)

    tr_levels = model_levels
    ntr_lv    = tr_levels

  END IF


  !-------------------------------------------------------
  ! Check that if l_CCRAD is selected certain other
  ! switches are not set so that they conflict.
  !-------------------------------------------------------
  IF (l_ccrad) THEN

    IF (.NOT. l_3d_cca) THEN
      icode    = 100
      cmessage = 'Error: CCRad only available with 3D'//        &
                       ' cloud field (l_3d_cca = .TRUE.)'
      CALL ereport(routinename, icode, cmessage)
    END IF

    IF (l_fix_udfactor) THEN
      icode    = 101
      cmessage = 'Error: l_ccrad and l_fix_udfactor'//          &
                       ' should not be both set to .TRUE.'
      CALL ereport(routinename, icode, cmessage)
    END IF

    IF (l_pc2_diag_sh) THEN
      icode    = 102
      cmessage = 'Error: l_ccrad and l_pc2_diag_sh'//           &
                       ' should not be both set to .TRUE.'
      CALL ereport(routinename, icode, cmessage)
    END IF

  END IF      ! l_ccrad

  IF (l_anvil) THEN
    IF (.NOT. l_3d_cca) THEN
      icode    = 104
      cmessage = 'Error: Anvil scheme requires 3D CCA'//        &
                       ' field (l_3d_cca = .TRUE.)'
      CALL ereport(routinename, icode, cmessage)
    END IF
  END IF

  IF (l_pc2) THEN
    ! So that PC2 does NOT affect Section 5 diagnostics
    l_dcpl_cld4pc2 = .TRUE.
  END IF

!-------------------------------------------------------
! End error capture
!-------------------------------------------------------


!=====================================================================
!     Call the main part of the model
!=====================================================================
! DEPENDS ON: scm_main
  CALL scm_main                                                               &
    ( vert_lev, nfor, l_ts_log, ntrop, sec_day                                &
    , land_points, nsprog, ntab, co2_dim_len, co2_dim_row                     &
    , cloud_levels, model_levels, wet_levels, tr_levels                       &
    , tr_vars, tr_ukca, st_levels, sm_levels, bl_levels                       &
    , ozone_levels, ntiles, nice, nice_use                                    &
    , atmos_sr, nsectp, l_netcdf_obs )

  IF (lhook) CALL dr_hook('SCM_SHELL',zhook_out,zhook_handle)

! Shut down parallel comms
  CALL gc_exit()

END PROGRAM scm_shell

