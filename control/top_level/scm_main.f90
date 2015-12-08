
! *********************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *********************************COPYRIGHT*********************************
!
!+ Routine to run the SCM.
!
! Subroutine Interface:
SUBROUTINE scm_main                                                           &
  ( vert_lev, nfor, l_ts_log, ntrop, sec_day, land_points, nsprog, ntab       &
  , co2_dim_len, co2_dim_row, cloud_levels, model_levels, wet_levels          &
  , tr_levels, tr_vars, tr_ukca, st_levels, sm_levels, bl_levels              &
  , ozone_levels, ntiles, nice, nice_use, atmos_sr, nsectp, l_netcdf_obs )

  ! Physics modules
  !---------------------------------------------------------------------------
  USE ancil_info,  ONLY: nsmax
  USE rad_input_mod
  USE sw_control_struct
  USE lw_control_struct
  USE ukca_radaer_struct_mod
  USE mphys_inputs_mod, ONLY: l_psd, l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup
  USE stochastic_physics_run_mod, ONLY: rhcrit_max, rhcrit_min                &
    , m_ci, m_ci_max, m_ci_min, g0_rp, par_mezcla, lambda_min_rp              &
    , ricrit_rp, a_ent_1_rp, g1_rp
  USE mphys_constants_mod, ONLY: ntot_land, ntot_sea
  USE cloud_inputs_mod, ONLY: rhcrit,        l_rhcpt,      l_cld_area,        &
                              l_acf_cusack,  l_acf_brooks, l_pc2
  USE river_inputs_mod, ONLY: river_vel, river_mcoef, l_inland, river_step,  &
                              i_river_vn
  USE mphys_bypass_mod
  USE g_wave_input_mod, ONLY: Gsharp, fbcd,l_smooth,l_nonhydro,l_dynbeta      &
    , l_gw_heating
  USE bl_option_mod
  USE turb_diff_mod, ONLY: qlimit, l_qpos
  USE turb_diff_ctl_mod, ONLY:                                                &
      visc_m, visc_h, rneutml_sq, max_diff, delta_smag, shear
  USE run_aerosol_mod, ONLY: so2_high_level
  USE arcl_mod,  ONLY: npd_arcl_compnts, npd_arcl_species, ip_arcl_sulp,      &
                       ip_arcl_dust, ip_arcl_sulp, ip_arcl_dust, ip_arcl_sslt,&
                       ip_arcl_blck, ip_arcl_biom, ip_arcl_ocff, ip_arcl_dlta
  USE ukca_option_mod, ONLY: L_ukca_chem, l_ukca_useumuivals,                 &
                             l_ukca_set_trace_gases, l_ukca_strat,            &
                             l_ukca_strattrop, l_ukca_prescribech4, l_ukca,   &
                             l_ukca_radaer
  USE nstypes, ONLY: npft
  ! SCM modules
  !---------------------------------------------------------------------------
  USE netcdf
  USE netcdf_obs, ONLY: get_netcdf_obs
  USE mcc_data
  USE global_scmop
  USE scm_utils
  USE s_main_force, l_spec_z0_nml=>l_spec_z0, timestep_nml=>timestep
  USE scm_cntl_mod

  ! Constants modules
  !---------------------------------------------------------------------------
  USE visbty_constants_mod, ONLY: n_vis_thresh, vis_thresh
  USE atmos_constants_mod,  ONLY: cp, r
  USE water_constants_mod,  ONLY: lc
  USE conversions_mod,      ONLY: pi_over_180, pi
  USE earth_constants_mod,  ONLY: g, earth_radius, omega, two_omega

  ! UM Modules
  !---------------------------------------------------------------------------
  USE Submodel_Mod, ONLY: n_internal_model
  USE atm_fields_bounds_mod
  USE vertnamelist_mod, ONLY:                                                 &
      first_constant_r_rho_level, z_top_of_model, eta_theta, eta_rho          &
    , vertlevs, max_number_alpha_cds, max_121_rows          &
    , max_updiff_levels, max_bl_levels, max_req_thpv_levs, max_sponge_width   &
    , max_look

  ! Trig arrays (requires dynamics A12_2A code)
  USE trignometric_mod, ONLY:                                                 &
      cos_theta_latitude, sec_theta_latitude, FV_cos_theta_latitude           &
    , sin_theta_latitude, cos_theta_longitude, sin_theta_longitude            &
    , true_longitude, true_latitude

  ! Model level heights from centre of Earth  (requires dynamics A12_2A code)
  USE level_heights_mod, ONLY:                                                &
      eta_theta_levels, eta_rho_levels, r_theta_levels, r_rho_levels,         &
      z_top_theta

  USE timestep_mod, ONLY :                                                    &
      timestep_number, timestep, radiation_timestep, radiation_tstep_diag     &
    , radiation_tstep_prog

  USE yomhook,  ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim


  USE ereport_mod, ONLY : ereport
  USE UM_ParVars
  USE domain_params
  USE trophgt1_mod, ONLY: z_min_trop, z_max_trop

  USE dynamics_input_mod, ONLY:                                &
      NumCycles,n_rims_to_do,l_regular,l_lbc_old

  USE dynamics_testing_mod, ONLY:                              &
      L_dry    
  USE dust_parameters_mod, ONLY: l_dust, l_dust_diag
  USE um_input_control_mod,  ONLY:                                            &
       model_domain,                                                          &
                            l_soot,             l_biomass,                    &
                            l_triffid,          l_phenol,                     &
                            l_ocff,             l_use_ocff_autoconv,          &
       l_nitrate,                               l_use_nitrate_autoconv,       &
       l_use_soot_autoconv, l_use_bmass_autoconv,                             &
       l_mr_physics1,       l_mr_physics2,      h_sect,                       &
                            lcal360

  USE eng_corr_inputs_mod, ONLY: lflux_reset
  USE murk_inputs_mod, ONLY: l_murk, l_murk_bdry, l_murk_rad
  USE cv_run_mod, ONLY: l_ccrad, l_3d_cca
  USE cosp_input_mod, ONLY: l_cosp
  USE switches, ONLY: l_snow_albedo, can_model

  USE nlstcall_mod, ONLY : ltimer
  USE chsunits_mod, ONLY : nunits
  IMPLICIT NONE

!
! Description:
!   Subroutine that runs the single column model.  This replaces the previous
!   version of scm_main that was the top level of the model.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code description:
!   FORTRAN 90 This code is written to UM programming standards version 8.2

! Arguments with Intent In. ie: Input variables.

  INTEGER, INTENT(In) :: nfor        ! Number terms for observational forcing
  INTEGER, INTENT(In) :: ntrop       ! Max number of levels in the troposphere
  INTEGER, INTENT(In) :: co2_dim_len
  INTEGER, INTENT(In) :: co2_dim_row
  INTEGER, INTENT(In) :: cloud_levels
  INTEGER, INTENT(In) :: model_levels
  INTEGER, INTENT(In) :: wet_levels
  INTEGER, INTENT(In) :: tr_levels
  INTEGER, INTENT(In) :: tr_vars
  INTEGER, INTENT(In) :: tr_ukca
  INTEGER, INTENT(In) :: st_levels
  INTEGER, INTENT(In) :: sm_levels
  INTEGER, INTENT(In) :: bl_levels
  INTEGER, INTENT(In) :: ozone_levels
  INTEGER, INTENT(In) :: ntiles
  INTEGER, INTENT(In) :: nice
  INTEGER, INTENT(In) :: nice_use
  INTEGER, INTENT(In) :: nsectp
  INTEGER, INTENT(In) :: sec_day
  INTEGER, INTENT(In) :: nsprog      ! No. of single level prognostics
  INTEGER, INTENT(In) :: ntab        ! Dimension of array used in random
                                     ! generator (Do not change this value as
                                     ! it is hard coded in the S_RANDOM deck)

  CHARACTER(LEN=200), INTENT(In) :: vert_lev
  CHARACTER(LEN=200)             :: sw_spec_file
  CHARACTER(LEN=200)             :: lw_spec_file
  CHARACTER(LEN=200)             :: file

  CHARACTER(LEN=2), INTENT(In) :: atmos_sr(0:nsectp)
  LOGICAL,          INTENT(In) :: &
    l_ts_log                      &! Output timestep information
  , l_netcdf_obs                   ! Using netcdf driver file for forcing

!=====================================================================

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

! degrees to radians & vice versa
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

! Required for tropopause min/max tropopause level calculation

! for N_INTERNAL_MODEL in typsts.h
! TYPSTS starts
! submodel_mod must be included before this file
!Applicable to all configurations
!STASH related variables for describing output requests and space
!management.
! Include sizes for dimensioning arrays in this deck
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS


! TYPSTSZ end

! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! Whether a calculation is needed for SF above
      LOGICAL :: SF_CALC(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
! TYPSTS end



! Local variables

  INTEGER  :: istatus

  INTEGER, PARAMETER :: height_gen_method = 2

  INTEGER  :: icode
  INTEGER  :: n_cca_lev
  INTEGER  :: errorstatus
  REAL     :: a2out(row_length,rows,model_levels)

  ! Arrays for PC2 initiation control arg list.
  REAL ::                                     &
    t_work   (row_length, rows, model_levels) &
  , q_work   (row_length, rows, wet_levels)   &
  , qcl_work (row_length, rows, wet_levels)   &
  , qcf_work (row_length, rows, wet_levels)   &
  , cf_work  (row_length, rows, wet_levels)   &
  , cfl_work (row_length, rows, wet_levels)   &
  , cff_work (row_length, rows, wet_levels)

!-----------------------------------------------------------------------------

  INTEGER, PARAMETER :: nconvars  = 20*mx_mod_lv        &
                                  + 12*mx_wet_lv        &
                                  + mx_tr_lv*mx_tr_vars &
                                  + 3

! Minimum no. of variables required to restart from a dump and is equal to :
  INTEGER, PARAMETER :: nprimvars = 4*mx_mod_lv         &
                                  + 4*mx_wet_lv         &
                                  + mx_st_lv            &
                                  + mx_sm_lv            &
                                  + mx_nsprog

!-----------------------------------------------------------------------------
! Dummy/Fixed variables
!-----------------------------------------------------------------------------
! These either have not impact on the SCM, or must be a set value
! for use with the SCM.

  INTEGER :: CycleNo = 1  ! Define CycleNo default

  INTEGER, PARAMETER ::  &
    n_proc   = 1         &! Total number of processors
  , n_procx  = 1         &! Number of processors in longitude
  , n_procy  = 1          ! Number of processors in latitude

!------------------------------------------------------------


  REAL, PARAMETER ::         &
    lat_rot_NP        = 0.0  &! Diagnostic variable
  , long_rot_NP       = 0.0   ! Diagnostic variable

  REAL, PARAMETER ::         &! STASH space for
    stashwork1        = 0.0  &! Section 1  (sw)
  , stashwork2        = 0.0  &! Section 2  (lw)
  , stashwork3        = 0.0  &! Section 3  (b layer)
  , stashwork4        = 0.0  &! Section 4  (ls precipitation)
  , stashwork5        = 0.0  &! Section 5  (convection)
  , stashwork6        = 0.0  &! Section 6  (gravity wave drag)
  , stashwork8        = 0.0  &! Section 8  (hydrology)
  , stashwork9        = 0.0  &! Section 9  (LS Cloud)
  , stashwork14       = 0.0  &! Section 14 (energy correction)
  , stashwork19       = 0.0  &! Section 19 (Veg)
  , stashwork26       = 0.0   ! Section 26 (River Routing)


! Duplication of variables needed to pass UKCA gases to radiation scheme.
! This option will not work in the SCM.
  INTEGER, PARAMETER         :: ngrgas     = 8
  INTEGER, DIMENSION(ngrgas) :: grgas_addr = -1


  INTEGER ::            &
    land_pts_trif  = 1  &! To dimension land fields
  , npft_trif      = 1  &! To dimension pft  fields
  , co2_dim_lev    = 1   ! 3-D CO2 field passed to NI_rad_ctl


! Declare pointers for dummy arrays
  REAL, POINTER ::           &
    lambda_p                 &! VarRes hor. co-ordinate information
  , lambda_u                 &! VarRes hor. co-ordinate information
  , dlambda_p                &! VarRes hor. co-ordinate information
  , wt_lambda_p              &! VarRes hor. co-ordinate information
  , wt_lambda_u               ! VarRes hor. co-ordinate information


  INTEGER, POINTER ::        &
    g_p_field                &! Size of global atmos field
  , g_r_field                &! Size of global river field
  , a_steps_since_riv        &! Number of physics timsteps since
                              ! last call to river routing
  , river_row_length         &! Local  river row length
  , river_rows               &! Local  river rows
  , global_river_row_length  &! Global river row length
  , global_river_rows         ! Global river rows


  REAL, POINTER, DIMENSION(:) :: &
    inlandout_atm                &! Inland basin flow (kg/m2/s)
  , tot_surf_runoff              &! Accumulated runoff over river
  , tot_sub_runoff               &! routing timestep (Kg/m2/s)
  , xpa                          &! Atmosphere TP long coordinates
  , xua                          &! Atmosphere U  long coordinates
  , xva                          &! Atmosphere V  long coordinates
  , ypa                          &! Atmosphere TP lat  coordinates
  , yua                          &! Atmosphere U  lat  coordinates
  , yva                           ! Atmosphere V  lat  coordinates


  REAL, POINTER, DIMENSION(:,:) :: &
    acc_lake_evap &! Accumulated lake evap over river routing timestep
  , r_area        &! Accumulated areas file
  , slope         &
  , flowobs1      &! Initialisation for flows
  , r_inext       &! x-coordinate of downstream grid point
  , r_jnext       &! y-coordinate of downstream grid point
  , r_land        &! Land/river/sea
  , substore      &! Routing sub_surface store (mm)
  , surfstore     &! Routing surface store (mm)
  , flowin        &! Rurface lateral inflow (mm)
  , bflowin       &! Rub-surface lateral inflow (mm)
  , trivdir       &! River direction file
  , trivseq       &! River sequence file
  , twatstor      &! Water storage
  , soil_clay     &! Fields for mineral dust
  , soil_silt     &! source flux calculations
  , soil_sand     &
  , dust_mrel1    &
  , dust_mrel2    &
  , dust_mrel3    &
  , dust_mrel4    &
  , dust_mrel5    &
  , dust_mrel6    &
  , phi_p         &! VarRes horizontal co-ordinate information
  , phi_v         &! VarRes horizontal co-ordinate information
  , dphi_p        &! VarRes horizontal co-ordinate information
  , wt_phi_p      &! VarRes horizontal co-ordinate information
  , wt_phi_v       ! VarRes horizontal co-ordinate information


  REAL, POINTER, DIMENSION(:,:,:) :: &
    dust_div1     &! Tracer variables
  , dust_div2     &
  , dust_div3     &
  , dust_div4     &
  , dust_div5     &
  , dust_div6

  ! UKCA_RADAER structure, for consistent call to Atmos_Physics1
  TYPE (ukca_radaer_struct) :: ukca_radaer

  ! Declare dummy target arrays

  REAL, TARGET ::                          &
    rdum0                                  &
  , rdum1(row_length)                      &
  , rdum2(row_length+1)                    &
  , rdum3(row_length, rows)                &
  , rdum4(row_length, rows, model_levels)  &
  , rdum5(row_length, rows, bl_levels)     &
  , rdum6(row_length, rows, bl_levels-1)   &
  , dum_land(land_points)

  INTEGER, TARGET :: idum0

  REAL ::                                  &
    frac_control    (land_points,ntype)    &
  , es_space_interp (4,row_length, rows)   &
  , o3_trop_level   (row_length,rows)      &
  , o3_trop_height  (row_length,rows)      &
  , t_trop_level    (row_length,rows)      &
  , t_trop_height   (row_length,rows)


  ! CDNC from UKCA-MODE, for consistent call to Atmos_Physics1
  ! This is currently a dummy array but could be used
  ! in future if UKCA is to be used with the SCM

  ! Declare dummy array and dimensions

  INTEGER, PARAMETER :: cdnc_dim1 = 1
  INTEGER, PARAMETER :: cdnc_dim2 = 1
  INTEGER, PARAMETER :: cdnc_dim3 = 1

  REAL :: cdnc_ukca_dummy( cdnc_dim1, cdnc_dim2, cdnc_dim3 )

!-----------------------------------------------------------------------------
! &INOBSFOR information
!-----------------------------------------------------------------------------
! We need to keep the variables that are read in from &INOBSFOR and reset
! them to their correct values at the start of each timestep so use dummy
! equivalent variables.

  REAL ::                                    &
    tls     (row_length,rows,model_levels)   &
  , uls     (row_length,rows,model_levels)   &
  , vls     (row_length,rows,model_levels)   &
  , wls     (row_length,rows,0:model_levels) &
  , qls     (row_length,rows,wet_levels)     &
  , qcl_inc (row_length,rows,model_levels)   &
  , qcf_inc (row_length,rows,model_levels)



!-----------------------------------------------------------------------------
! Time information
!-----------------------------------------------------------------------------
  INTEGER ::     &
    ihour        &! Current hour
  , imin         &! Current min
  , isec = 0      ! Current sec, not implemented

!-----------------------------------------------------------------------------
! Initial day of the year and initial time (in seconds) in that day
!-----------------------------------------------------------------------------
  INTEGER ::      &
    dayno_init    &! Initial day
  , time_initi     ! Initial time in seconds (integer)

  REAL ::         &
    time_init      ! Initial time in seconds (real)


!-----------------------------------------------------------------------------
! Primary Model Variables plus P_star (UMDP No1)
!-----------------------------------------------------------------------------

  INTEGER ::                &
    iccb(row_length, rows)  &! Convective cloud base
  , icct(row_length, rows)   ! Convective cloud top

  REAL ::                                 &
    canopy_gb(land_points)                &! Canopy water content(kg/m2)
  , cca(row_length,rows,model_levels)     &! Convective cloud amount
  , q(row_length,rows,wet_levels)         &! Specific humidity (kg/kg)
  , qcf(row_length,rows,wet_levels)       &! Cloud ice content (kg/kg)
  , qcl(row_length,rows,wet_levels)       &! Cloud water content(kg/kg)
  , qcf2(row_length,rows,wet_levels)      &! 2nd ice water content (kg/kg)
  , qrain(row_length,rows,wet_levels)     &! Rain water content (kg/kg)
  , qgraup(row_length,rows,wet_levels)    &! Graupel water content(kg/kg)
  , mix_v(row_length,rows,wet_levels)     &! Vapour mixing ratio (kg/kg)
  , mix_cf(row_length,rows,wet_levels)    &! Cloud ice content (kg/kg)
  , mix_cl(row_length,rows,wet_levels)    &! Cloud water content(kg/kg)
  , mix_cf2(row_length,rows,wet_levels)   &! 2nd ice water content (kg/kg)
  , mix_rain(row_length,rows,wet_levels)  &! Rain water content (kg/kg)
  , mix_graup(row_length,rows,wet_levels)  ! Graupel water content(kg/kg)

  REAL ::                                            &
    exner_rho_levels(row_length,rows,model_levels+1) &
  , exner_theta_levels(row_length,rows,model_levels) &
  , exner_prime(row_length,rows,model_levels)         ! Increment to exner

  REAL ::                                         &
    dtheta_dr_term(row_length,rows, model_levels) &
  , rho(row_length,rows,model_levels)             &
  , rho_n(row_length,rows,model_levels)           &
  , p(row_length,rows,model_levels+1)             &! Pressure on rho levels
  , p_theta_levels(row_length,rows,model_levels)  &! Pressure on theta levels
  , rp(row_length,rows,model_levels+1)            &! 1/p on rho levels
  , rp_theta(row_length,rows,model_levels)         ! 1/p on theta levels

  REAL ::                          &
    p_star(row_length,rows)        &! Surface pressure
  , smc(land_points)               &! Soil moisture content (Kg/m^2)
  , smcl(land_points,sm_levels)    &! Soil moisture content in layers (Kg/m^2)
  , sthf(land_points,sm_levels)    &! Frozen soil moisture content of each
                                    ! layer as a fraction of saturation.
  , sthu(land_points,sm_levels)     ! Unfrozen soil moisture content of each
                                    ! layer as a fraction of saturation.

!
!---------------------------------------------------------------------
!     Water
!---------------------------------------------------------------------
  REAL ::                                 &
    snodep(row_length,rows)               &! Snow depth (Kg/m^2)
  , t(row_length,rows,model_levels)       &! Temperature (K)
  , t_deep_soil(land_points,st_levels)    &! Deep soil temperatures (K)
                                           ! top level not included, =surface
  , tsi(row_length,rows)                  &! Temperature of sea-ice
  , tstar(row_length,rows)                &! Surface temperature (K)
  , u(row_length,rows,model_levels)       &! Zonal wind (m/s)
  , v(row_length,rows,model_levels)       &! Meridional wind (m/s)
  , w(row_length,rows,0:model_levels)     &! Vertical velocity (m/s)
  , z0msea(row_length,rows)               &! Sea surface roughness length
  , w_adv(row_length,rows,0:model_levels)  ! Advective w component of wind

  REAL :: &
    ozone_tracer(row_length,rows,model_levels)  ! Ozone tracer

  REAL ::                         &
    deep_flag(row_length,rows)    &! Indicator of how long since last
                                   ! deep convective event.
  , past_precip(row_length,rows)  &! Decayed memory of past total
                                   ! convective precipitation
  , past_conv_ht(row_length,rows)  ! Memory of past height of convection


  REAL ::                                            &
    theta_star     (row_length, rows, model_levels)  &
  , qcl_star       (row_length, rows, wet_levels)    &
  , qcf_star       (row_length, rows, wet_levels)    &
  , qcf2_star      (row_length, rows, wet_levels)    &
  , qrain_star     (row_length, rows, wet_levels)    &
  , qgraup_star    (row_length, rows, wet_levels)    &
  , mix_v_star     (row_length, rows, wet_levels)    &
  , mix_cl_star    (row_length, rows, wet_levels)    &
  , mix_cf_star    (row_length, rows, wet_levels)    &
  , mix_cf2_star   (row_length, rows, wet_levels)    &
  , mix_rain_star  (row_length, rows, wet_levels)    &
  , mix_graup_star (row_length, rows, wet_levels)    &
  , cf_star        (row_length, rows, wet_levels)    &
  , cfl_star       (row_length, rows, wet_levels)    &
  , cff_star       (row_length, rows, wet_levels)    &
  , uinc_geo       (row_length, rows, model_levels)  &
  , vinc_geo       (row_length, rows, model_levels)

!-----------------------------------------------------------------------------
! CCRad Prognostics
!-----------------------------------------------------------------------------
  INTEGER :: &
    lcbase (row_length, rows) ! Model level of lowest convective cloud base
                              ! in profile(theta level) passed to radiation

  REAL :: &
    ccw (row_length, rows, wet_levels)
                              ! Convective Cloud Water profile passed to
                              ! radiation scheme (theta levels). (kg/kg)
                              ! NOTE: May be different to ccw produced by
                              !       Convection Scheme it Radiative
                              !       convective cloud decay is used.


!-----------------------------------------------------------------------------
! PC2 cloud and condensation scheme
!-----------------------------------------------------------------------------
  REAL ::                                      &
    q_forcing   (row_length, rows, wet_levels) &
  , qcl_forcing (row_length, rows, wet_levels) &
  , t_forcing   (row_length, rows, wet_levels) &
  , p_forcing   (row_length, rows, wet_levels) &
  , q_earliest  (row_length, rows, wet_levels) &
  , qcl_earliest(row_length, rows, wet_levels) &
  , t_earliest  (row_length, rows, wet_levels)


!-----------------------------------------------------------------------------
! Water
!-----------------------------------------------------------------------------
  REAL ::                     &
    ls_rain(row_length, rows) &! Large scale rainfall rate (Kg/m^2)
  , ls_snow(row_length, rows)  ! Large scale snowfall rate (Kg/m^2/s)

!---------------------------------------------------------------------
! JULES2 snow scheme prognostics
!---------------------------------------------------------------------
  REAL ::                             &
    snowdepth(    land_points,ntiles) &! Snow depth on ground on tiles (m)
  , rho_snow_grnd(land_points,ntiles) &! Snowpack bulk density (kg/m3)
! , nsnow(        land_points,ntiles) &! No. of snow layers on ground on tiles
                                       ! NOTE: this is converted to an integer

  , ds(     land_points,ntiles,nsmax) &! Snow layer thickness (m)
  , sice(   land_points,ntiles,nsmax) &! Snow layer ice  mass on tiles (Kg/m2)
  , sliq(   land_points,ntiles,nsmax)  ! Snow layer liq. mass on tiles (Kg/m2)

  REAL ::                                &
    rho_snow  (land_points,ntiles,nsmax) &! Snow layer densities (kg/m3)
  , tsnowlayer(land_points,ntiles,nsmax) &! Snow layer temperature (K)
  , rgrainl   (land_points,ntiles,nsmax)  ! Snow layer grain size on tiles
                                          ! (microns)

!-----------------------------------------------------------------------------
! Large scale statistical forcing
!-----------------------------------------------------------------------------
! Random generator variables

  INTEGER ::            &
    iv(ntab), iy, idum  &! Contains info on generator
  , iseed                ! Seed for random number generator

  INTEGER ::            &
    dayno_wint           ! Day number relative to winter solstice

  REAL ::                                &
    ad(row_length, rows,wet_levels-1)    &! Term a of equation 2.22
                                          !  for dewpoint depression
  , alfada(row_length,rows)              &! Amplitude and mean of seasonal
  , alfadb(row_length,rows)              &!  variation of tuning factor
  , at(row_length,rows,model_levels-1)   &! Term a of equation 2.22
                                          !  for dewpoint depression
  , atime, btime                         &! Constants for calculating annual
                                          ! cycle
  , avn(row_length,rows,model_levels-1)  &! Term a of equation 2.22
  , aw(row_length,rows,ntrop-1)          &!  for horiz. and vert. vel.
  , cdbar(row_length,rows,wet_levels)    &! Mean and SD of random variable
  , cdsd(row_length,rows,wet_levels)     &!  for dew pt. depression
  , ctbar(row_length,rows,model_levels)  &! Mean and SD of random variable
  , ctsd(row_length,rows,model_levels)   &!  for temp.
  , cvnbar(row_length,rows,model_levels) &! Mean and SD of random variable
  , cvnsd(row_length,rows,model_levels)  &!  for velocity VN
  , cwbar(row_length,rows,ntrop)         &! Mean and SD of random variable
  , cwsd(row_length,rows,ntrop)          &!  for vertical velocity
  , dbar(row_length,rows,wet_levels)     &! Mean dewpoint depression at
                                          !  daycount days from winter
                                          !  solstice (K)
  , dbara(row_length,rows,wet_levels)    &! Amplitude and mean of seasonal
  , dbarb(row_length,rows,wet_levels)    &!  variation of mean dew pt.
                                          ! depression (K)
  , ddash(row_length,rows,wet_levels)    &! Dew pt depression correction (K)
  , deltan(row_length,rows)              &! Radius of area (m)
  , dgrada(row_length,rows,wet_levels)   &! Amplitude and mean of seasonal
  , dgradb(row_length,rows,wet_levels)   &!  variation of dew pt. depression
                                          !  gradient (K/km)
  , dsd(row_length,rows,wet_levels)      &! SD dewpoint depression at daycount
                                          !  days from winter solstice (K)
  , pa(row_length,rows,model_levels+1)   &! Amplitude and mean of seasonal
  , pb(row_length,rows,model_levels+1)   &!  variation of pressure (Pa)
  , px(row_length,rows,ntrop)            &! Reciprocal log functions for
  , py(row_length,rows,ntrop-1)          &!  calc. of vert. advection
                                          !  used in eqns 2.12 and 2.13
  , qr(row_length,rows,wet_levels,2)     &! Randomly sampled specific
                                          !  humidity (kg/kg)
  , tbar(row_length,rows,model_levels)   &! Mean temperature at daycount days
                                          !  from winter solstice (K)
  , tbara(row_length,rows,model_levels)  &! Amplitude and mean of seasonal
  , tbarb(row_length,rows,model_levels)  &!  variation of temp. (K)
  , tdash(row_length,rows,model_levels)  &! Temp correction (K)
  , tgrada(row_length,rows,model_levels) &! Amplitude and mean of seasonal
  , tgradb(row_length,rows,model_levels) &!  variation of temp. gradient
                                          !  (k/Km)
  , tr(row_length,rows,model_levels,2)   &! InOut Randomly sampled temp. (K)
  , tsd(row_length,rows,model_levels)    &! SD temp. at daycount days
                                          ! from winter solstice (K)
  , tsda(row_length,rows,model_levels)   &! Amplitude and mean of seasonal
  , tsdb(row_length,rows,model_levels)   &!  variation of SD of temp. (K)
  , vnbar(row_length,rows,model_levels)  &! Mean velocity VN at daycount days
                                          !  from winter solstice
  , vnbara(row_length,rows,model_levels) &! Amplitude and mean of seasonal
  , vnbarb(row_length,rows,model_levels) &!  variation of velocity VN (m/s)
  , vnr(row_length,rows,model_levels,2)  &! InOut Randomly sampled horizontal
                                          !       velocity (m/s)
  , vnsd(row_length,rows,model_levels)   &! SD velocity VN at daycount days
                                          !  from winter solstice (m/s)
  , vnsda(row_length,rows,model_levels)  &! Amplitude and mean of seasonal
  , vnsdb(row_length,rows,model_levels)  &!  variation of SD of velocity VN
                                          !  (m/s)
  , vpbar(row_length,rows,model_levels)  &! Mean velocity VP at daycount days
                                          !  from winter solstice
  , vpbara(row_length,rows,model_levels) &! Amplitude and mean of seasonal
  , vpbarb(row_length,rows,model_levels) &!  variation of velocity VP (m/s)
  , vpr(row_length,rows,model_levels,2)  &! InOut Randomly sampled horizontal
                                          !  velocity (m/s)
  , wbar(row_length,rows,ntrop)          &! Mean vertical velocity at daycount
                                          !  days from winter solstice (m/s)
  , wbara(row_length,rows,ntrop)         &! Amplitude and mean of seasonal
  , wbarb(row_length,rows,ntrop)         &!  variation of SD of vert. vel.
                                          !  (mb or HPa/s)
  , wr(row_length,rows,ntrop,2)          &! InOut Randomly sampled vertical
                                          !  velocity (mb/s)
  , wsd(row_length,rows,ntrop)           &! SD vertical velocity at daycount
                                          !  days from winter solstice (m/s)
  , wsda(row_length,rows,ntrop)          &! Amplitude and mean of seasonal
  , wsdb(row_length,rows,ntrop)           !  variation of SD of vert. vel.
                                          !  (mb/s) roughness length (m)

!-----------------------------------------------------------------------------
! Large scale observational forcing
!-----------------------------------------------------------------------------
! Variables for diagnostic output for observational forcing

  REAL ::                                        &
    dap1(row_length,rows,36,model_levels)        &! Instantaneous profiles
  , dap2(row_length,rows,36,model_levels)        &! Mean profiles
  , dap3(row_length,rows,36,nfor-1,model_levels) &! Mean profiles - timeseries
  , dab1(row_length,rows,44)                     &! Instantaneous budgets
  , dab2(row_length,rows,44)                     &! Mean budgets
  , dab3(row_length,rows,44,nfor-1)              &! Mean budgets - timeseries
  , deltap(row_length,rows,model_levels)         &! Layer thickness (Pa)
  , factor_rhokh(row_length,rows)                 ! Used to specify surface
                                                  ! flux from observation

  REAL ::                              &
    resdump(row_length,rows,nprimvars)  ! Dump array of restart variables


!-----------------------------------------------------------------------------
! Variables enabling diagnostic output
!-----------------------------------------------------------------------------

  INTEGER, PARAMETER ::   &
    nSCMdpkgs = 13         ! No of diags packages

! Note that if nSCMDpkgs is changed it should also be changed
! in PC2_PRES and ATMSTEP2

  LOGICAL ::              &
    l_scmdiags(nscmdpkgs)  ! Logicals for diagnostics packages

  LOGICAL ::              &
    kill_scmdiags(nscmdpkgs)
                           ! Disable diagnostics for geostrophic
                           ! initialisation call

! Parameters used in calls to scmoutput
! Start of include file: s_scmop.h
! Description:
!  Declares and defines some parameters necessary for calling SCMoutput
!
!

! Integers to represent the different time profiles. All must
! be non-negative and less than "only_radsteps".

  INTEGER, PARAMETER :: &
    t_inst        = 1   &! Give the instantaneous value
  , t_avg         = 2   &! Construct the average value
  , t_max         = 3   &! " maximum value
  , t_min         = 4   &! " minimum value
  , t_acc         = 5   &! " accumulated value
  , t_div         = 7   &! " average value divided
                         !   by another diagnostic
  , t_mult        = 8   &! " average value multiplied
                         !   by another diagnostic
  , t_acc_div     = 9   &! " accumulated value divided
                         !   by another diagnostic
  , t_acc_mult    = 10  &! " accumulated value multiplied
                         !   by another diagnostic
  , t_const       = 11  &! The value is constant.
  , only_radsteps = 100  ! When added to one of the above parameters,
                         ! flags that the diagnostic is only available
                         ! on radiation timesteps

! Integers to represent the different domain profiles
  INTEGER, PARAMETER :: &
    d_sl      = 1       &
  , d_soilt   = 2       &
  , d_bl      = 3       &
  , d_wet     = 4       &
  , d_all     = 5       &
  , d_soilm   = 6       &
  , d_tile    = 7       &
  , d_vis     = 9       &
  , d_point   = 13      &
  , d_allxtra = 14      &
  , d_land    = 15      &
  , d_cloud   = 16

! Statement function to encode a stream number into an integer
  INTEGER :: &
    Stream   &
  , strm

  Stream(strm) = 2**(strm-1)

! The default streams for diagnostics to go to will be 1,2,3,4,5 and 6.
! The following should thus be equal to:
!
! Stream(1) [2^0=1] + Stream(2) [2^1=2]  + Stream(3) [2^2=4]
! Stream(4) [2^3=8] + Stream(5) [2^4=16] + Stream(6) [2^5=32]
! Total = 63
!
! where Stream() is the statement function defined above.
! default is 63 (all)

  INTEGER, PARAMETER :: &
    default_streams = 63

! Integers to represent the different diagnostics packages
  INTEGER, PARAMETER :: &
    SCMDiag_gen   = 1   & ! General diagnostics
  , SCMDiag_rad   = 2   & ! Radiation
  , SCMDiag_bl    = 3   & ! Boundary layer
  , SCMDiag_surf  = 4   & ! Surface
  , SCMDiag_land  = 5   & ! Land points only
  , SCMDiag_sea   = 6   & ! Sea points only
  , SCMDiag_lsp   = 7   & ! Large scale precip
  , SCMDiag_conv  = 8   & ! Convection
  , SCMDiag_lscld = 9   & ! Large scale cloud
  , SCMDiag_pc2   = 10  & ! PC2
  , SCMDiag_forc  = 11  & ! Forcing
  , SCMDiag_incs  = 12  & ! Increments
  , SCMDiag_gwd   = 13    ! Gravity Wave Drag

! End of include file: s_scmop.h

  INTEGER ::              &
    site(row_length*rows)  ! SSFM site WMO number

! Arrays to store initial fields for calculating total increments and
! total increment fields
  REAL ::                                      &
    t_start(row_length,rows,model_levels)      &
  , u_start(row_length,rows,model_levels)      &
  , v_start(row_length,rows,model_levels)      &
  , w_start(row_length,rows,0:model_levels)    &
  , q_start(row_length,rows,wet_levels)        &
  , qcl_start(row_length,rows,wet_levels)      &
  , qcf_start(row_length,rows,wet_levels)      &
  , cf_start(row_length,rows,wet_levels)       &
  , cfl_start(row_length,rows,wet_levels)      &
  , cff_start(row_length,rows,wet_levels)

  REAL ::                                      &
    t_totalinc(row_length,rows,model_levels)   &
  , u_totalinc(row_length,rows,model_levels)   &
  , v_totalinc(row_length,rows,model_levels)   &
  , w_totalinc(row_length,rows,0:model_levels) &
  , q_totalinc(row_length,rows,wet_levels)     &
  , qcl_totalinc(row_length,rows,wet_levels)   &
  , qcf_totalinc(row_length,rows,wet_levels)   &
  , cf_totalinc(row_length,rows,wet_levels)    &
  , cfl_totalinc(row_length,rows,wet_levels)   &
  , cff_totalinc(row_length,rows,wet_levels)


! Qpos arrays
!---------------------------------------------------------------------
! Hold incs to q to prevent -ve values.
  REAL ::                                      &
    dq_qpos_for (row_length,rows,wet_levels)   &
  , dq_qpos     (row_length,rows,wet_levels)   &
  , dqcl_qpos   (row_length,rows,wet_levels)   &
  , dqcf_qpos   (row_length,rows,wet_levels)


!-----------------------------------------------------------------------------
! Extra diagnostics boundary layer code
!-----------------------------------------------------------------------------

  REAL ::                             &
    resp_w_ft(row_length, rows,npft)   ! Out Wood maintenance respiration
                                       !     (Kg C m^-2 s^-1).

!-----------------------------------------------------------------------------
! Clouds
!-----------------------------------------------------------------------------

  REAL ::                                    &
    cf(row_length,rows,wet_levels)           &! Layer cloud amount
                                              ! (decimal fraction)
  , cfl(row_length,rows,wet_levels)          &! liquid layer cloud amount
  , cff(row_length,rows,wet_levels)          &! frozen layer cloud amount
  , area_cloud_fraction(row_length,rows,wet_levels)  &
                                              ! area cloud fraction
  , ccwpin(row_length,rows)                  &! Condensed water path (Kg/m^2)
  , rhts(row_length,rows,wet_levels)         &! Initial relative humidity
                                              ! wrt TL
  , qtts(row_length,rows,wet_levels)         &! Initial q_T
  , tlts(row_length,rows,wet_levels)         &! Initial TL
  , ptts(row_length,rows,wet_levels)         &! Initial p on theta levs
  , tl(row_length,rows)                      &! Liquid water temperature (K)
  , qsl_tl(row_length,rows)                   ! Qsat wrt liquid water
                                              ! at temp TL

!-----------------------------------------------------------------------------
! Radiation
!-----------------------------------------------------------------------------
  INTEGER ::         &
    daynumber        &! Day in the year (default=1)
  , year              ! Year


!-----------------------------------------------------------------------------
! Used in calculation of min/max tropopause
!-----------------------------------------------------------------------------
  REAL :: r_ref_rho(model_levels)  ! height of model rho levels above surface.


!-----------------------------------------------------------------------------
! Loop Counters and limits
!-----------------------------------------------------------------------------
  INTEGER ::       &
    daycount       &! Counts through days
  , daysteps       &! Number of timestep in a day
  , full_daysteps  &! Number of timesteps in a full day
  , nstepsin       &! Number of steps in final day
  , total_nsteps   &! Total number of steps in run
  , i, j, j2, k    &! General loop counters, array dumpmean
  , itype          &
  , m1              ! No. of dumps

!-----------------------------------------------------------------------------
! Miscellaneous
!-----------------------------------------------------------------------------
  CHARACTER(LEN=80) :: &
    cmessage         ! Error message if Icode >0

  CHARACTER(LEN=*), PARAMETER :: routinename = 'scm_main'

  CHARACTER(LEN=8) ::    &
    time_string      ! String containing actual time

  INTEGER ::        &
    error_code      &
  , day             &
  , yearno          &
  , time_sec         ! Actual time of day in secs.

  REAL ::                        &
    modug(row_length,rows)       &! Magnitude of geostrophic wind (m/s)
  , f_coriolis(row_length,rows)  &! 2*omega*sin(latitude)
  , f3_at_u(row_length,rows)     &
  , rccb(row_length,rows)        &! Real val. cloud base geoint only
  , rcct(row_length,rows)        &! Real cloud top geoint only
  , maxinc                        ! Maximum wind increment from geoinit

  REAL ::                                     &
    sw_incs(row_length,rows,0:model_levels+1) &
  , lw_incs(row_length,rows,0:model_levels)   &
  , dirpar_incs(row_length,rows)              &! PAR flux variable
  , t1_sd(row_length,rows)                    &! Set to zero initially
  , q1_sd(row_length,rows)                    &! Set to zero initially
  , sw_tile_rts(row_length*rows,ntiles)       &! Surface net SW on land tiles
  , rhokh(row_length,rows,bl_levels)
    ! Exchange coeffs for moisture.
    !       surface:out of SF_EXCH
    !               contains only RHOKH.
    ! above surface:out of KMKH
    !               contains GAMMA(1)*RHOKH(,1)*RDZ(,1)

  REAL ::                           &
    gridbox_area_m(row_length,rows)  ! Gridbox area in m^2

  INTEGER ::                        &
    land_points                     &! No of land points - can be 0
  , land_ice_points                 &! No of ice land points
  , soil_points                     &! No of soil points
  , land_index(row_length*rows)     &! Should be defined by land_points
  , land_ice_index(row_length*rows) &! but land_points has not yet been
  , soil_index(row_length*rows)      ! given a value

  ! Land point specific surface soil parameters
  REAL ::                &
    b_exp  (land_points) &! Clapp-Hornberger exponent
  , hcap   (land_points) &! Soil heat capacity
  , hcon   (land_points) &! Soil thermal conductivity
  , satcon (land_points) &! Saturated hydrological conductivity
  , sathh  (land_points) &! Saturated soil water suction
  , v_sat  (land_points) &! Volumetric SMC at saturation
  , v_wilt (land_points) &! Volumetric SMC at wilting point
  , v_crit (land_points)  ! Volumetric SMC at critical point


  REAL ::                                               &
    energy_correction                                   &
  , ukca_tracers(row_length,rows,tr_levels,tr_ukca)

  REAL ::                                               &
    biogenic (row_length,rows,model_levels)

! Local arrays for using the aerosol climatology for NWP

! Internal model switches
  LOGICAL :: &
    l_use_arcl(npd_arcl_species)

! Internal array of mass-mixing ratios
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: arcl


  INTEGER ::        &
    n_arcl_species  &! Number of requested species within the climatology
  , n_arcl_compnts   ! Corresponding number of requested components

  INTEGER :: &
    i_arcl_compnts(npd_arcl_compnts) ! Array index of each component

  REAL :: &
    unscl_dry_rho(row_length,rows,model_levels) ! unscaled dry density

  REAL ::                                       &
    rhokm   (row_length,rows,0:bl_levels-1)     &
  , ch_term (row_length,rows,model_levels-1)

  REAL ::                                       &
    photosynth_act_rad(row_length, rows)        &! Net downward shortwave
                                                 ! radiation in band 1 (w/m2).
  , rad_hr(row_length, rows, 2, bl_levels)      &! BL radiative heating rates
  , micro_tends(row_length, rows, 2, bl_levels) &! BL Microphys tendencies
  , dOLR(row_length, rows)                      &! TOA - surface upward LW
  , sw_tile(row_length*rows, ntiles)            &! Surface net SW on land tiles
  , cos_zenith_angle(row_length, rows)

  !---------------------------------------------------------------------
  ! Rate of change in forcing increments (1/s)
  !---------------------------------------------------------------------
  REAL ::                                            &
    ch_t_inc (row_length,rows,nfor-1,model_levels)   &! Temperature
  , ch_u_inc (row_length,rows,nfor-1,model_levels)   &! Zonal wind
  , ch_v_inc (row_length,rows,nfor-1,model_levels)   &! Meridional wind
  , ch_w_inc (row_length,rows,nfor-1,0:model_levels) &! Vertical wind
  , ch_q_star(row_length,rows,nfor-1,wet_levels)     &! Specific Humidity
  , ch_flux_e(row_length,rows,nfor-1)                &! Latent heat flux
  , ch_flux_h(row_length,rows,nfor-1)                &! Sensible heat flux
  , ch_tstar_forcing(row_length,rows,nfor-1)          ! Srf. Temperature

  REAL ::                                            &
    ch_ug(row_length,rows,nfor-1,model_levels)       &! Geostrophic u-wind
  , ch_vg(row_length,rows,nfor-1,model_levels)        ! Geostrophic v-wind

  REAL ::                                            &
    ch_q_bg (row_length,rows,nfor-1,wet_levels)      &! Observed q
  , ch_t_bg (row_length,rows,nfor-1,model_levels)    &! Observed T
  , ch_u_bg (row_length,rows,nfor-1,model_levels)    &! Observed u-wind
  , ch_v_bg (row_length,rows,nfor-1,model_levels)    &! Observed v-wind
  , ch_w_bg (row_length,rows,nfor-1,0:model_levels)   ! Observed w-velocity

  !---------------------------------------------------------------------
  ! Forcings passed to physics routines
  !---------------------------------------------------------------------

  REAL ::                          &
    flux_e_scm (row_length, rows)  &! Srf. Lat heat fluxes
  , flux_h_scm (row_length, rows)   ! Srf. Sen heat fluxes

  REAL ::                                       &
    ug_scm (row_length,rows,model_levels)       &! Geostrophic u-wind (m/s)
  , vg_scm (row_length,rows,model_levels)        ! Geostrophic v-wind (m/s)

  ! Increments/timestep due to large-scale horizontal and vertical advection
  REAL ::                                       &
    t_inc_scm  (row_length,rows,model_levels)   &! Temperature inc (K)
  , u_inc_scm  (row_length,rows,model_levels)   &! Zonal wind inc (m/s)
  , v_inc_scm  (row_length,rows,model_levels)   &! Merid wind inc (m/s)
  , w_inc_scm  (row_length,rows,0:model_levels) &! Vert  wind inc (m/s)
  , q_star_scm (row_length,rows,wet_levels)      ! Spec. humid. inc (kg/kg)


  ! Background field profiles (for relaxation)
  REAL ::                                    &
    t_bg_scm(row_length,rows,model_levels)   &! Temperature (K)
  , u_bg_scm(row_length,rows,model_levels)   &! Zonal wind (m/s)
  , v_bg_scm(row_length,rows,model_levels)   &! Merid wind (m/s)
  , w_bg_scm(row_length,rows,0:model_levels) &! Vert  wind (m/s)
  , q_bg_scm(row_length,rows,wet_levels)      ! Spec. humid (kg/kg)

  LOGICAL ::         &
    l_rad_step       &! True if radiation  timestep
  , l_rad_step_diag  &! True if diagnostic timestep
  , l_rad_step_prog  &! True if prognostic timestep
  , l_mr_pc2         &! True if using mixing ratios in PC2 section
  , l_calc_exner     &
  , l_calc_rho

  LOGICAL ::         &
    rad_mask(row_length,rows)

! Land based variables
  REAL ::                    &
    fexp       (land_points) &! Decay factor in Sat. conductivity in water
                              ! table layer.
  , gamtot     (land_points) &! Integrated complete Gamma function.
  , ti_mean    (land_points) &! Mean topographic index.
  , ti_sig     (land_points) &! St dev. of topographic index. in water
                              ! table layer
  , fsat       (land_points) &! Surface saturation fraction.
  , fwetl      (land_points) &! Wetland fraction.
  , zw         (land_points) &! Water table depth (m).
  , sthzw      (land_points) &! Soil moist fract. in deep-zw layer.
  , a_fsat     (land_points) &! Fitting parameter for Fsat in LSH model
  , c_fsat     (land_points) &! Fitting parameter for Fsat in LSH model
  , a_fwet     (land_points) &! Fitting parameter for Fwet in LSH model
  , c_fwet     (land_points) &! Fitting parameter for Fwet in LSH model
  , catch_snow (land_points, ntiles) &! Snow capacity for NLT tile (kg/m2)
  , snow_grnd  (land_points, ntiles)  ! Snow below canopy (kg/m2).


  CHARACTER(LEN=300) :: sdum0, sdum1, sdum2

  INTEGER :: asteps_since_triffid
  INTEGER :: previous_time(7)      ! Time information for current timestep

  REAL ::         &
    delta_lambda  &
  , delta_phi

! Variables for COSP, for consistency in calls to atmos_physics1/2
! Convective rainfall
  REAL :: cosp_crain_3d(row_length,rows,model_levels)
! Convective snowfall
  REAL :: cosp_csnow_3d(row_length,rows,model_levels)

! Ancillary fields and fields needed to be kept from timestep to timestep
! Assuming here that nice=nice_use
  REAL ::                                &
    ice_fract_n(row_length, rows, nice)  &! Ice categories
  , ice_thick_n(row_length, rows, nice)  &
  , ti_n       (row_length, rows)        & ! Ice temperature as grid box mean
  , ice_k_n(row_length, rows, nice)      & ! Effective cond of sea ice
  , tstar_sice_n(row_length, rows, nice_use) & ! Surface ice T
  , snodep_sice_n(row_length, rows, nice_use)  ! Snow depth on ice

! Add new variables from extra args in calls to physics routines
  REAL :: &
    rhcpt(row_length, rows, wet_levels)

! Work array for scmoutput for ug and vg
  REAL :: &
    geo_diag(row_length, rows, model_levels)

! TKE based turbulence scheme
  REAL ::                                            &
    e_trb(row_length, rows, model_levels)            &
!                   ! TKE defined on theta levels K-1
  , tsq_trb(row_length, rows, model_levels)          &
!                   ! Self covariance of liquid potential temperature
!                   ! (thetal'**2) defined on theta levels K-1
  , qsq_trb(row_length, rows, model_levels)          &
!                   ! Self covariance of total water
!                   ! (qw'**2) defined on theta levels K-1
  , cov_trb(row_length, rows, model_levels)          &
!                   ! Correlation between thetal and qw
!                   ! (thetal'qw') defined on theta levels K-1
  , zhpar_shcu(row_length, rows)
!                   ! Height of mixed layer used to evaluate
!                   ! the non-gradient buoyancy flux


  INTEGER :: &
    ichgf    &! No. of timesteps between change in observational forcing
  , ilscnt    ! Count for observational forcing

  INTEGER :: &
    nout(19)  ! Output units to receive printout of initial
              ! data by PRINT_INITDATA

  REAL ::    &
    ti (row_length,rows,model_levels) ! Initial temperature profile, (K)


  ! Dr Hook
  !==============================
  REAL(KIND=jprb) :: zhook_handle

  IF (lhook) CALL dr_hook('SCM_MAIN',zhook_in,zhook_handle)

  ! Assign variables normally for MPP configs


  CALL decompose(1,1,0,0,model_levels,noSmallHalos=.TRUE.)
  CALL Change_Decomposition(decomp_smexe)
  at_extremity      = .FALSE.

  ! Assign Pointers to dummy target arrays
  g_p_field               => idum0
  g_r_field               => idum0
  a_steps_since_riv       => idum0
  river_row_length        => idum0
  river_rows              => idum0
  global_river_row_length => idum0
  global_river_rows       => idum0
  lambda_p                => rdum0
  phi_p                   => rdum3
  lambda_u                => rdum0
  phi_v                   => rdum3
  dlambda_p               => rdum0
  dphi_p                  => rdum3
  wt_lambda_p             => rdum0
  wt_lambda_u             => rdum0
  wt_phi_p                => rdum3
  wt_phi_v                => rdum3
  dust_div1               => rdum4
  dust_div2               => rdum4
  dust_div3               => rdum4
  dust_div4               => rdum4
  dust_div5               => rdum4
  dust_div6               => rdum4
  acc_lake_evap           => rdum3
  r_area                  => rdum3
  slope                   => rdum3
  flowobs1                => rdum3
  r_inext                 => rdum3
  r_jnext                 => rdum3
  r_land                  => rdum3
  substore                => rdum3
  surfstore               => rdum3
  flowin                  => rdum3
  bflowin                 => rdum3
  trivdir                 => rdum3
  trivseq                 => rdum3
  twatstor                => rdum3
  soil_clay               => rdum3
  soil_silt               => rdum3
  soil_sand               => rdum3
  dust_mrel1              => rdum3
  dust_mrel2              => rdum3
  dust_mrel3              => rdum3
  dust_mrel4              => rdum3
  dust_mrel5              => rdum3
  dust_mrel6              => rdum3
  inlandout_atm           => dum_land
  tot_surf_runoff         => dum_land
  tot_sub_runoff          => dum_land
  xpa                     => rdum2
  xua                     => rdum2
  xva                     => rdum2
  ypa                     => rdum1
  yua                     => rdum1
  yva                     => rdum2

! Initialise Dummy/Hardwired values

  ice_fract_n  = 0.0    ! Was uninitialised before
  ice_thick_n  = 0.0    ! Was uninitialised before
  ti_n         = 273.0  ! Ice temperature as grid box mean
  ice_k_n      = 41.8   ! 2*kappai/de
  tstar_sice_n = 273.0
  snodep_sice_n = 0.0
  model_domain = 5      ! Specifies model type
  n_rims_to_do = 1      ! Size of LAM rim zone
  l_regular    = .TRUE. ! True if NOT variable resolution
  l_lbc_old    = .TRUE. ! True if old lbcs (not active scm)
  numcycles    = 1      ! cycling dynamics-physics to be compatible
                        ! with UM

  i_river_vn         = 1
  fexp               = 1.0
  gamtot             = 1.0
  ti_mean            = 1.0
  ti_sig             = 1.0
  fsat               = 1.0
  fwetl              = 1.0
  zw                 = 1.0
  sthzw              = 1.0
  a_fsat             = 1.0
  c_fsat             = 1.0
  a_fwet             = 1.0
  c_fwet             = 1.0
  catch_snow(:,:)    = 1.0
  snow_grnd(:,:)     = 1.0

  cdnc_ukca_dummy(:,:,:) = 0.0

! Allocate Smagorinsky variables
  ALLOCATE ( visc_h(1,1,1) )
  ALLOCATE ( visc_m(1,1,1) )
  ALLOCATE ( rneutml_sq(1,1,1) )
  ALLOCATE ( shear(1,1,1) )
  ALLOCATE ( max_diff(1,1) )
  ALLOCATE ( delta_smag(1,1) )

! Hardwired model switches for SCM

  l_mr_pc2     = .FALSE. ! Use mixing ratios in PC2 section
  l_murk_bdry  = .FALSE. ! Atmos_physics1 call
  lflux_reset  = .FALSE. ! Atmos_physics1 call
  l_rad_deg    = .FALSE. ! Atmos_physics1 call
  l_rhcpt      = .FALSE. ! Atmos_physics2 call
  l_inland     = .FALSE. ! Re-route inland basin flow to soil moisture

  l_phenol  = .FALSE.
  rad_mask  = .TRUE.


! Set dimensions of aerosol climatologies - aerosols are not used
! in the SCM so this are hardwired to be unused.
!
! If aerosols are enabled in the future then these should be correctly
! dimensioned using set_arcl_dimensions and set_arcl_clim as in
! atm_step

  n_arcl_species = 0
  n_arcl_compnts = 1  ! Set to one to prevent out of bounds error in radiation

  l_use_arcl(ip_arcl_biom) = .FALSE.
  l_use_arcl(ip_arcl_blck) = .FALSE.
  l_use_arcl(ip_arcl_sslt) = .FALSE.
  l_use_arcl(ip_arcl_sulp) = .FALSE.
  l_use_arcl(ip_arcl_dust) = .FALSE.
  l_use_arcl(ip_arcl_ocff) = .FALSE.
  l_use_arcl(ip_arcl_dlta) = .FALSE.

  ALLOCATE ( arcl(1,1,1,1) )

! Mid Latitude values, (ref J-C-Thelen)
  o3_trop_level  = 22
  o3_trop_height = 11000
  t_trop_level   = 22
  t_trop_height  = 11000

! Initialise dummy arguments
  rdum0    = 1.0
  rdum1    = 1.0
  rdum2    = 1.0
  rdum3    = 1.0

  rdum4    = 1.0
  rdum5    = 1.0
  rdum6    = 1.0

  idum0    = 1
  dum_land = 1.0

  WRITE(sdum0,*) nfor
  WRITE(sdum1,*) mx_nobs

  CALL alloc_common_scm
  CALL alloc_forcing
  CALL set_forcing_defaults

  ti = 0.0

!=============================================================================
! Set up vertical level information
!=============================================================================

  OPEN(10, FILE=vert_lev, IOSTAT=istatus, STATUS='old')

  IF (istatus /= 0) THEN

    icode=500
    WRITE(Cmessage,*)  " Error opening " //TRIM(ADJUSTL(vert_lev))//     &
                       " file on unit 10"
    CALL ereport (routinename, icode, cmessage)

  END IF

  READ(10,VERTLEVS)

  CLOSE(10)

! Initialise error code
  error_code = 0

!-----------------------------------------------------------------------------
! Calculate values of grid radius/heights
!-----------------------------------------------------------------------------
  ! Note full model does this in subroutine setcona
  ALLOCATE( r_rho_levels   (row_length,rows,  model_levels))
  ALLOCATE( r_theta_levels (row_length,rows,0:model_levels))

  IF (.NOT.allocated(eta_theta_levels)) THEN
        ALLOCATE (eta_theta_levels(0:model_levels))
  END IF
  IF (.NOT.allocated(eta_rho_levels)) THEN
        ALLOCATE (eta_rho_levels(model_levels))
  END IF

  eta_theta_levels(0:model_levels) = eta_theta (1:model_levels+1)
  eta_rho_levels(:) = eta_rho (:)
  z_top_theta = z_top_of_model

  ! DEPENDS ON: calc_levels
  CALL calc_levels                                                       &
     ! In
     ( orog, height_gen_method                                           &
     , bl_levels, model_levels                                           &
     , rows, row_length )


  z_th(1:nmod_lv) = r_theta_levels(1,1,1:nmod_lv) - r_theta_levels(1,1,0)
  z_rh(1:nmod_lv) = r_rho_levels(1,1,1:nmod_lv)   - r_theta_levels(1,1,0)
  z_rh(nmod_lv+1) = 2*z_th(nmod_lv) - z_rh(nmod_lv)


! Set radiation switches in rad_input_mod

  lrad_ccrad         = l_ccrad

! Set microphysics switches used in mphys_bypass_mod
l_crystals      = l_mcr_qcf2 .OR. (.NOT. l_psd)
mphys_mod_top   = z_top_of_model


!---------------------------------------------------------------------
!     Initialise the array giving the unit nos. for output of the
!     initial data and anything else that needs setting
!---------------------------------------------------------------------

  asteps_since_triffid = 0
  omega                = 7.292116E-5
  two_omega            = 2 * 7.292116E-5
  nout(:)              = 0

! Initialise
  rad_hr       = 0.0
  micro_tends  = 0.0

  ccwpin   = 1.0
  ls_rain  = 0.0
  ls_snow  = 0.0

  t_bg_scm = 0.0
  q_bg_scm = 0.0
  u_bg_scm = 0.0
  v_bg_scm = 0.0
  w_bg_scm = 0.0

  ug_scm(:,:,:) = 0.0
  vg_scm(:,:,:) = 0.0

  q_star_scm = 0.0
  t_inc_scm  = 0.0
  u_inc_scm  = 0.0
  v_inc_scm  = 0.0
  w_inc_scm  = 0.0

  flux_h_scm = 0.0
  flux_e_scm = 0.0
  qcl_inc    = 0.0
  qcf_inc    = 0.0

! Initialise past history for convection

  DO j=1, rows
    DO i=1, row_length
      deep_flag(i,j)    = 0.0      ! No previous deep convection
      past_precip(i,j)  = 0.0      ! No previous convective precip.
      past_conv_ht(i,j) = 0.0      ! No previous convection.
    END DO
  END DO

! Initialise the prognostic variables in the TKE schemes
! The missing values are automatically initialised in the first call of
! the TKE schemes.
  e_trb   = rmdi
  tsq_trb = rmdi
  qsq_trb = rmdi
  cov_trb = rmdi

! ilscnt hardwired initialisation to zero until further
! development to correct functionality. Also removed from
! &INOBSFOR namelist to prevent incorrect usage.

  ilscnt = 0

  ! No of levels for Convective Cloud Amount.
  IF (l_3d_cca) THEN
    n_cca_lev = wet_levels
  ELSE
    n_cca_lev = 1
  END IF

!-----------------------------------------------------------------------------
! Write out user notification of ancillary files
!-----------------------------------------------------------------------------

  sw_spec_file = sw_control(1)%spectral_file
  lw_spec_file = lw_control(1)%spectral_file

  WRITE(6,*) ' '
  WRITE(6,*) ' Using Ancillary files:'
  WRITE(6,'(A)')                                                  &
             ' ----------------------------------------'//        &
             '-----------------------------------------'
  WRITE(6,'(A)') ' Vertical levels  : '//TRIM(ADJUSTL(vert_lev))
  WRITE(6,'(A)') ' Forcing namelist : '//TRIM(ADJUSTL(scm_nml))
  WRITE(6,'(A)')                                                  &
             ' ----------------------------------------'//        &
             '-----------------------------------------'
  WRITE(6,*) ' '


  IF (l_netcdf_obs) THEN
    CALL get_netcdf_obs
  END IF

  ! DEPENDS ON: read_scm_nml
  CALL read_scm_nml

  IF (timestep_nml /= RMDI) THEN
    timestep = timestep_nml
  END IF

  DO i=1, row_length
    DO j=1, rows
      gridbox_area_m(i,j) = gridbox_area(i,j) * 1e+6
      deltan(i,j)         = SQRT(gridbox_area_m(i,j)/pi)
    END DO
  END DO

  !---------------------------------------------------------------------------
  ! Use McClatchey atmospheric profiles where there is misssing data
  !---------------------------------------------------------------------------
  CALL alloc_mcc(nmod_lv)
  CALL get_mcc(lat, month_init)

  ! UM levels ABOVE available data/obs_top
  !---------------------------------------
  DO k=1, nmod_lv
    IF (z_th(k) > obs_top) THEN
      ozone   (:,:,k) = mcc_o3(k)
      theta   (:,:,k) = mcc_th(k)
      t_bg  (:,:,:,k) = mcc_th_t(k)
      t_inc (:,:,:,k) = 0.0
      qi      (:,:,k) = mcc_q(k)
      q_bg  (:,:,:,k) = mcc_q(k)
      q_star(:,:,:,k) = 0.0
    ELSE
      WHERE (ozone  (:,:,k) == RMDI) ozone  (:,:,k) = mcc_o3(k)
      WHERE (theta  (:,:,k) == RMDI) theta  (:,:,k) = mcc_th(k)
      WHERE (t_bg (:,:,:,k) == RMDI) t_bg (:,:,:,k) = mcc_th_t(k)
      WHERE (qi     (:,:,k) == RMDI) qi     (:,:,k) = mcc_q(k)
      WHERE (q_bg (:,:,:,k) == RMDI) q_bg (:,:,:,k) = mcc_q(k)
    END IF

    ! Check variables on rho levels
    IF (z_rh(k) > obs_top) THEN
      ui(:,:,k) = ui(:,:,k-1)
      vi(:,:,k) = vi(:,:,k-1)
      wi(:,:,k) = 0.0

      u_inc(:,:,:,k) = u_inc(:,:,:,k-1)
      v_inc(:,:,:,k) = v_inc(:,:,:,k-1)
      w_inc(:,:,:,k) = 0.0
    END IF
  END DO

  DO k=1, nmod_lv+1
    IF (z_rh(k) > obs_top) THEN
      p_in(:,:,k) = mcc_rh_p(k)
    ELSE
      WHERE (p_in(:,:,k) == RMDI) p_in(:,:,k) = mcc_rh_p(k)
    END IF
  END DO


  ! UM levels BELOW available data/obs_top
  !---------------------------------------
  ! NOTE: Wind profiles and forcing data on UM levels below the obs_bot
  !       (default lowest level in namelist.scm) are set to values at
  !       obs_bot.  Users wishing to prescribe values should add
  !       additional levels to their namelist.scm
  !---------------------------------------------------------------------------
  DO k=nmod_lv, 1, -1
    IF (z_th(k) < obs_bot) THEN
      t_inc (:,:,:,k) = t_inc (:,:,:,k+1)
      q_star(:,:,:,k) = q_star(:,:,:,k+1)
      t_bg  (:,:,:,k) = t_bg  (:,:,:,k+1)
      q_bg  (:,:,:,k) = q_bg  (:,:,:,k+1)
      w_inc (:,:,:,k) = w_inc (:,:,:,k+1)
      wi    (:,:,k)   = wi    (:,:,k+1)
    END IF

    IF (z_rh(k) < obs_bot) THEN
      ui(:,:,k) = ui(:,:,k+1)
      vi(:,:,k) = vi(:,:,k+1)

      u_inc(:,:,:,k) = u_inc(:,:,:,k+1)
      v_inc(:,:,:,k) = v_inc(:,:,:,k+1)
    END IF
  END DO

  ! Impose zero w wind at top of model
  wi(:,:,nmod_lv)      = 0.0
  w_inc(:,:,:,nmod_lv) = 0.0


  CALL dealloc_mcc

  !---------------------------------------------------------------------------
  ! Now all input has been read in from combination of netcdf/scm
  ! forcing namelist, check settings
  !---------------------------------------------------------------------------

  IF    ((stats    .AND. obs)     .OR.(stats .AND. noforce)         &
    .OR. (stats    .AND. geoforce).OR.(obs   .AND. geoforce)        &
    .OR. (geoforce .AND. noforce) .OR.(obs   .AND. noforce)) THEN

    WRITE (*,'(4(TR1,A52/),4(TR1,A12,TR1,L1,T53,A1/),4(TR1,A52/))') &
         '====================================================',    &
         '| Warning: More than one forcing logical is set to |',    &
         '| true. This may produce unexpected results.       |',    &
         '|                                                  |',    &
         '|  noforce: ', noforce,                           '|',    &
         '| geoforce: ', geoforce,                          '|',    &
         '|    stats: ', stats,                             '|',    &
         '|      obs: ', obs,                               '|',    &
         '|                                                  |',    &
         '===================================================='

  END IF ! Test for multiple forcing options


  IF     ((.NOT. stats)   .AND. (.NOT. obs)                         &
    .AND. (.NOT. noforce) .AND. (.NOT. geoforce)) THEN

    WRITE (*,'(4(TR1,A52/))')                                       &
         '====================================================',    &
         '| No forcing logical set.                          |',    &
         '| Setting noforce to true as default               |',    &
         '===================================================='
    noforce = .TRUE.

  END IF ! Test for at least one forcing option.

  IF (l_use_dust .OR. l_dust .OR. l_dust_diag) THEN
    WRITE (*,'(8(TR1,A52/))')                                       &
          '====================================================',   &
          '| Mineral dust scheme not yet implemented in SCM.  |',   &
          '| Setting following logicals to false:             |',   &
          '|                                                  |',   &
          '|  l_use_dust                                      |',   &
          '|  l_dust                                          |',   &
          '|  l_dust_diag                                     |',   &
          '|                                                  |',   &
          '===================================================='
    l_use_dust  = .FALSE.
    l_dust      = .FALSE.
    l_dust_diag = .FALSE.

  END IF ! Test for dust scheme options


  IF    (l_use_soot_direct   .OR. l_use_soot_indirect               &
    .OR. l_use_soot_autoconv .OR. l_soot) THEN
    WRITE (*,'(8(TR1,A52/))')                                       &
          '====================================================',   &
          '| Soot chemistry not yet implemented in the SCM.   |',   &
          '| Setting following logicals to false:             |',   &
          '|                                                  |',   &
          '|  l_use_soot                                      |',   &
          '|  l_soot                                          |',   &
          '|                                                  |',   &
          '===================================================='
    l_soot              = .FALSE.
    l_use_soot_direct   = .FALSE.
    l_use_soot_indirect = .FALSE.
    l_use_soot_autoconv = .FALSE.

  END IF ! Test for soot chemistry options


  IF    (l_biomass            .OR. l_use_bmass_direct               &
    .OR. l_use_bmass_indirect .OR. l_use_bmass_autoconv) THEN

    WRITE (*,'(10(TR1,A52/))')                                      &
          '====================================================',   &
          '| Biomass scheme not available in the SCM.         |',   &
          '| Setting following logicals to false:             |',   &
          '|                                                  |',   &
          '|  l_biomass                                       |',   &
          '|  l_use_bmass_direct                              |',   &
          '|  l_use_bmass_indirect                            |',   &
          '|  l_use_bmass_autoconv                            |',   &
          '|                                                  |',   &
          '===================================================='
    l_biomass             = .FALSE.
    l_use_bmass_direct    = .FALSE.
    l_use_bmass_indirect  = .FALSE.
    l_use_bmass_autoconv  = .FALSE.

  END IF

  IF    (l_ocff               .OR. l_use_ocff_direct                &
    .OR. l_use_ocff_indirect  .OR. l_use_ocff_autoconv) THEN

    WRITE (*,'(10(TR1,A52/))')                                      &
         '====================================================',    &
         '| Fossil-fuel OC scheme not available in the SCM.  |',    &
         '| Setting following logicals to false:             |',    &
         '|                                                  |',    &
         '|  l_ocff                                          |',    &
         '|  l_use_ocff_direct                               |',    &
         '|  l_use_ocff_indirect                             |',    &
         '|  l_use_ocff_autoconv                             |',    &
         '|                                                  |',    &
         '===================================================='
    l_ocff                = .FALSE.
    l_use_ocff_direct     = .FALSE.
    l_use_ocff_indirect   = .FALSE.
    l_use_ocff_autoconv   = .FALSE.

  END IF

  IF    (l_nitrate              .OR. l_use_nitrate_direct           &
    .OR. l_use_nitrate_indirect .OR. l_use_nitrate_autoconv) THEN

    WRITE (*,'(10(TR1,A52/))')                                      &
          '====================================================',   &
          '| Nitrate scheme not available in the SCM.         |',   &
          '| Setting following logicals to false:             |',   &
          '|                                                  |',   &
          '|  l_nitrate                                       |',   &
          '|  l_use_nitrate_direct                            |',   &
          '|  l_use_nitrate_indirect                          |',   &
          '|  l_use_nitrate_autoconv                          |',   &
          '|                                                  |',   &
          '===================================================='
    l_nitrate               = .FALSE.
    l_use_nitrate_direct    = .FALSE.
    l_use_nitrate_indirect  = .FALSE.
    l_use_nitrate_autoconv  = .FALSE.

  END IF

  IF (l_ukca_radaer) THEN

    WRITE (*,'(7(TR1,A52/))')                                       &
          '====================================================',   &
          '| UKCA_RADAER scheme not available in the SCM.     |',   &
          '| Setting following logical to false:              |',   &
          '|                                                  |',   &
          '|  l_ukca_radaer                                   |',   &
          '|                                                  |',   &
          '===================================================='
    l_ukca_radaer           = .FALSE.

  END IF

!-----------------------------------------------------------------------------
! Set l_cosp to false because it is not currently available
! for use in the SCM
!-----------------------------------------------------------------------------
  IF (l_cosp) THEN
    l_cosp = .FALSE.
    icode = -1
    cmessage=' COSP is not currently available in the SCM.'//       &
             ' Setting l_cosp to false.'
    CALL ereport (routinename, icode, cmessage)
  END IF

  ! Check consistency on land masks/land_points
  !============================================
  land_ice_points = 0
  soil_points     = 0
  k = 0

  DO j=1, rows
    DO i=1, row_length
      IF (land_sea_mask(i,j)) THEN
        k = k + 1
        land_index(k) = (j-1)*row_length + i

        IF (land_ice_mask(i,j)) THEN
          land_ice_points = land_ice_points + 1
          land_ice_index(land_ice_points) = (j-1)*row_length+i
        END IF

        IF (soil_mask(i,j)) THEN
          soil_points = soil_points + 1
          soil_index(soil_points) = (j-1)*row_length+i
        END IF

      ELSE
        IF (land_ice_mask(i,j) .OR. soil_mask(i,j)) THEN
          WRITE (*,'(9(TR1,A52/))')                                       &
            '===================================================='        &
          , '| LAND_ICE_MASK/SOIL_MASK in inconsistent with sea |'        &
          , '| point (land_sea_mask=.FALSE.)                    |'        &
          , '| Setting following logicals to false:             |'        &
          , '|                                                  |'        &
          , '|  land_ice_mask                                   |'        &
          , '|  soil_mask                                       |'        &
          , '|                                                  |'        &
          , '===================================================='

          land_ice_mask(i,j) = .FALSE.
          soil_mask(i,j)     = .FALSE.

        END IF
      END IF
    END DO
  END DO


  IF (land_points /= k) THEN
    Icode = 507
    WRITE (6,'(4(TR1,A52/))')                                             &
      '===================================================='              &
    , '| Specified total number of land points and        |'              &
    , '| land_sea_mask are inconsisent                    |'              &
    , '===================================================='
    WRITE(6,*) " "
    WRITE(6,*) "Run ABORTED"
    WRITE(6,*) " "

    WRITE(Cmessage,*) ' Check that land_points and land_sea_mask in'//    &
                      ' namelist are consistent.'
    CALL ereport (routinename, icode, cmessage)

  END IF

  IF (land_points > 0) THEN

    DO j=1, ntype
      DO i=1, land_points
        IF (frac_typ(i,j) == RMDI) THEN
          WRITE (6,'(2(TR1,A52/), (TR1,A13,I2,A30,T53,A1/),'//         &
                   ' 2(TR1,A52/))')                                    &
            '====================================================',    &
            '| Land tile fractions (frac_typ) has not been      |',    &
            '| fully set. ',ntype,' fractions must be specified','|',  &
            '| for each land point in namelist.                 |',    &
            '===================================================='
          cmessage = 'Land tile fractions improperly set'
          icode = 1
          CALL ereport (routinename, icode, cmessage)
        END IF
      END DO
    END DO

    DO i=1, land_points
      IF (SUM(frac_typ(i,:)) /= 1.0) THEN
        WRITE (6,'(4(TR1,A52/))')                                     &
          '====================================================',     &
          '| Total tile_fractions (frac_typ) must sum to 1.0  |',     &
          '| for each land point.                             |',     &
          '===================================================='
        cmessage = 'Land tile fractions do not total 1.0'
        icode = 1
        CALL ereport (routinename, icode, cmessage)
      END IF

      IF (soil_type(i) == IMDI) THEN
        WRITE (6,'(4(TR1,A52/))')                                     &
          '====================================================',     &
          '| Soil types (soil_type) must be specified by user |',     &
          '| for each land point.                             |',     &
          '===================================================='
        cmessage = 'land point has unspecified soil type'
        icode = 1
        CALL ereport (routinename, icode, cmessage)
      END IF
    END DO

  END IF ! test on land_points > 0

  IF (obs .AND. .NOT. l_netcdf_obs) THEN
    IF (nfor == IMDI) THEN
      icode=503
      WRITE(6,*) " "
      WRITE(6,*) "============================================="
      WRITE(6,*) " Number of observational forcings (nfor) in  "
      WRITE(6,*) " SCM namelist (&CNTLSCM) has not been set.   "
      WRITE(6,*) "============================================="
      WRITE(6,*) " "

      CLOSE(10)

      WRITE(cmessage,*) ' Variable NFOR has not been set'
      CALL ereport(routinename, icode, cmessage)

    ELSE IF (nfor > mx_nobs) THEN
      icode=504
      WRITE(6,*) " "
      WRITE(6,*) "============================================="
      WRITE(6,*) " Specified nfor(" // TRIM(ADJUSTL(sdum0)) //              &
                 ") > mx_nobs("     // TRIM(ADJUSTL(sdum1)) // ")."
      WRITE(6,*) " This WILL produce incorrect forcings."
      WRITE(6,*) " Please check your forcing file: "
      WRITE(6,*) "   " // TRIM(ADJUSTL(scm_nml))
      WRITE(6,*) "============================================="
      WRITE(6,*) " "

      CLOSE(10)

      WRITE(cmessage,*)                                                     &
                 " Specified nfor(" // TRIM(ADJUSTL(sdum0)) //              &
                 ") > mx_nobs("     // TRIM(ADJUSTL(sdum1)) // ")."
      CALL ereport(routinename, icode, cmessage)

    END IF
  END IF !  obs .and. .not. l_netcdf_obs


  IF ((ndayin == IMDI) .OR. (nminin == IMDI) .OR. (nsecin == IMDI)) THEN

    ! Run full length of forcing provided using
    ! observed profiles and observation period.
    ndayin =  INT((nfor-1)*obs_pd/86400.0)
    nminin = (INT((nfor-1)*obs_pd - ndayin*86400.0))/60.0
    nsecin = (nfor-1)*obs_pd - ndayin*86400.0 - nminin*60.0

  END IF

  IF (INT(sec_day/timestep) /= REAL(sec_day)/timestep) THEN

    ! Print a warning if we're going to run for more than a day.
    IF (ndayin + REAL(nminin)/(60.0*24.0)                           &
               + REAL(nsecin)/(3600.0*24.0) >  1.0) THEN
      PRINT*,'************************************************'
      PRINT*,'* Warning: Your timestep is such that you have *'
      PRINT*,'*          have requested a non-integer number *'
      PRINT*,'*          of steps per day.                   *'
      PRINT*,'************************************************'
    END IF
  END IF

  full_daysteps = INT(sec_day/timestep)
  nstepsin      = INT((nminin*60 + nsecin)/timestep)
  total_nsteps  = ndayin * full_daysteps + nstepsin
  scm_timestep  = timestep

  IF (MOD(obs_pd,timestep) == 0.0) THEN
    ichgf = INT(obs_pd/timestep)
  ELSE
    PRINT*,'Specified run timestep is not'
    PRINT*,'a factor of the obs period.  '
    RETURN
  END IF

  WRITE(6,'(A80)')                                                  &
    ' ======================================='//                    &
    '========================================'
  WRITE(6,'(A11,T22,F7.1)')                                         &
    ' Timestep: ', timestep
  WRITE(6,'(A21,T22,F7.1)')                                         &
    ' Observation period: ', obs_pd
  WRITE(6,'(A13,T22,I2,A1,I2,A1,I4,A3,I2,A1,I2)')                   &
    ' Start date: ', day_init,'/', month_init,'/', year_init, ',  ' &
                   , hour_init, ':', min_init
  WRITE(6,'(A13,T22,I2,A6,TR1,I2,A5,I2,A5)')                        &
    ' Run length: ', ndayin, ' days,', INT(nminin/60.0), ' hrs,'    &
                   , nminin - INT(nminin/60.0)*60, ' mins'
  WRITE(6,'(A80)')                                                  &
    ' ======================================='//                    &
    '========================================'


!-----------------------------------------------------------------------------
! Code taken from SETCONA
! Calculation of min_trop_level and max_trop_level
! NOTE: min and max_trop_level are used if climatological aerosols
! are chosen.
!-----------------------------------------------------------------------------
!   The tropopause diagnosed for radiative purposes divides theta-levels
!   considered to lie in the stratosphere from those considered to lie
!   in the troposphere: the tropopause is therefore taken as a
!   rho-level. This level is constrained to lie between heights of
!   z_min_trop and z_max_trop. The level is used in the setting of the
!   background aerosol climatology and of ozone profiles, subject to the
!   choice of appropriate options; additionally it is used in the
!   calculation of diagnostics defined at the tropopause.
!
!   Start at the second rho-level because the first is omitted from the
!   grid seen in the physics.
!-----------------------------------------------------------------------------
!   This code is only used if min_trop_level and max_trop_level in
!   &RUNDATA namelist are set both set to 0.  Otherwise values in scm
!   namelist are used.
!-----------------------------------------------------------------------------

  IF ((min_trop_level == 0 .AND. max_trop_level == 0) .OR.                &
       min_trop_level <  0  .OR. max_trop_level <  0) THEN

    min_trop_level=2

    DO k=1, model_levels
      r_ref_rho(k) = eta_rho(k)*z_top_of_model
    END DO

    DO ; IF ((r_ref_rho(min_trop_level) >= z_min_trop) .OR.               &
            (min_trop_level == model_levels)) EXIT
         min_trop_level = min_trop_level+1
    END DO

    max_trop_level=min_trop_level
    DO ; IF ( (r_ref_rho(max_trop_level) > z_max_trop) .OR.               &
            (max_trop_level == model_levels) ) EXIT
         max_trop_level = max_trop_level+1
    END DO

    max_trop_level = max_trop_level-1
  END IF

!-----------------------------------------------------------------------------

! Derive the initial daynumber in the year and the initial time
! from the UM data structure supplied

! DEPENDS ON: inittime_scm
  CALL inittime_scm( dayno_init, time_initi, lcal360 )

  time_init = time_initi    ! type real expected elsewhere.

! Initial tape daynumber in the year is the same as the initial
! daynumber

  tapeday_init = dayno_init

  previous_time(1) = year_init
  previous_time(2) = month_init
  previous_time(3) = day_init
  previous_time(4) = hour_init
  previous_time(5) = min_init
  previous_time(6) = sec_init
  previous_time(7) = dayno_init

  ! Set time in scm_utils
  time_info = previous_time


!-----------------------------------------------------------------------------
! Set longitude to zero so that diagnostics apply to local time
! rather than GMT and convert to radians
!-----------------------------------------------------------------------------

 ALLOCATE( true_latitude(row_length,rows))
 ALLOCATE( true_longitude(row_length,rows))

  DO j=1, rows
    DO i=1, row_length
      true_latitude(i,j) = pi_over_180 * lat(i,j)
      IF (local_time) THEN
        true_longitude(i,j) = 0.0
      ELSE
        true_longitude(i,j) = pi_over_180 * long(i,j)
      END IF
    END DO
  END DO


! Allocate and Calculate trig arrays - stored in module
! Full model does this in routine setcona
  ALLOCATE( cos_theta_latitude (row_length,rows) )
  ALLOCATE( sec_theta_latitude (row_length,rows) )
  ALLOCATE( sin_theta_latitude (row_length,rows) )

  ALLOCATE( cos_theta_longitude (row_length,rows) )
  ALLOCATE( sin_theta_longitude (row_length,rows) )
  ALLOCATE( fv_cos_theta_latitude (row_length,rows) )

  DO j=1, rows
    DO i=1, row_length
      cos_theta_latitude    (i,j) = COS(true_latitude(i,j))
      sec_theta_latitude    (i,j) = 1.0/cos_theta_latitude(i,j)
      cos_theta_longitude   (i,j) = COS(true_longitude(i,j))
      sin_theta_latitude    (i,j) = SIN(true_latitude(i,j))
      sin_theta_longitude   (i,j) = SIN(true_longitude(i,j))
      FV_cos_theta_latitude (i,j) = cos_theta_latitude(i,j)
    END DO
  END DO

! Initialise various arrays
  SW_incs     = 0.0
  LW_incs     = 0.0
  dirpar_incs = 0.0
  t1_sd       = 0.0
  q1_sd       = 0.0

! Set h_sect - for use in choosing radiation scheme.
! Code copied from set_h_sect (readlsa2.dk)
  DO i=0, nsectp
    h_sect(i)(1:1) = '0'
    h_sect(i)(2:3) = atmos_sr(i)
  END DO

! Determine number of advection timesteps between the radiation 
! diagnostics/prognostics steps
  a_sw_radstep_diag = ( secs_per_day / i_sw_radstep_perday_diag ) / timestep
  a_sw_radstep_prog = ( secs_per_day / i_sw_radstep_perday_prog ) / timestep
  a_lw_radstep_diag = ( secs_per_day / i_lw_radstep_perday_diag ) / timestep
  a_lw_radstep_prog = ( secs_per_day / i_lw_radstep_perday_prog ) / timestep
  IF (i_rad_extra_call == 0) THEN 
    a_sw_radstep_diag = a_sw_radstep_prog
    a_lw_radstep_diag = a_lw_radstep_prog
  END IF
! Determine radiation diagnostic and prognostic timesteps
  radiation_tstep_diag = timestep * a_sw_radstep_diag 
  radiation_tstep_prog = timestep * a_sw_radstep_prog

! In SCM A_LW_SEGMENTS should be set to 1 and segment sizes for
! LW and SW should be unset (i.e. -99)
  a_lw_segments = 1
  a_lw_seg_size = -99
  a_sw_seg_size = -99

!-----------------------------------------------------------------------------
! Set l_triffid to false because it is not currently available
! for use in the SCM (UM vn6.2). See Chris D Jones of the
! Terrestrial Carbon Cycle Group, Hadley Centre if you wish to
! use dynamical vegetation.
!-----------------------------------------------------------------------------
  IF (l_triffid) THEN
    WRITE(6,*) '*************************************************'
    WRITE(6,*) '* Warning: you have requested 19_2A interactive *'
    WRITE(6,*) '* vegetation distribution, which is not         *'
    WRITE(6,*) '* currently available in the SCM.               *'
    WRITE(6,*) '* Setting l_triffid to false                    *'
    WRITE(6,*) '*************************************************'
    l_triffid = .FALSE.
  END IF

!-----------------------------------------------------------------------------
! Calculate values of delta_lambda and delta_phi
!-----------------------------------------------------------------------------
! The SCM sets delta_lambda and delta_phi to 0 as a default.
! Calculate these, instead, from gridbox_area_m

   delta_lambda = 0
   delta_phi    = 0

   DO j=1, rows
     DO i=1, row_length
       delta_lambda = SQRT ( gridbox_area_m(i,j) /                            &
                    ( r_theta_levels(i,j,0) * r_theta_levels(i,j,0)           &
                     * fv_cos_theta_latitude(i,j) ) )
       delta_phi = delta_lambda
     END DO
   END DO

! Set values of mp_dell and mp_delp for microphysics
mp_dell = delta_lambda
mp_delp = delta_phi

!-----------------------------------------------------------------------------
! Calculate values of exner
!-----------------------------------------------------------------------------
  l_calc_exner = .TRUE.
  l_calc_rho   = .TRUE.

! DEPENDS ON: calc_press
  CALL calc_press                                                             &
    ! (In)
    ( model_levels, wet_levels, rows, row_length, p_in, theta, qi             &
    , l_calc_exner, l_calc_rho                                                &
    ! (InOut)
    , rho                                                                     &
    ! (Out)
    , exner_theta_levels, exner_rho_levels, p_theta_levels, rp, rp_theta      &
    , p_star )


!-----------------------------------------------------------------------------
! Convert thetai to ti
!-----------------------------------------------------------------------------
  ti(:,:,:) = theta(:,:,:)*exner_theta_levels(:,:,:)


!-----------------------------------------------------------------------------
! Set up the unit nos. for output.
!-----------------------------------------------------------------------------

  nout(1)=6
  IF (test)          nout(2) = 22
  IF (prindump_step) nout(3) = 30
  IF (prindump_day)  nout(4) = 31

  IF (prindump_days) THEN
    IF (dump_days(1) > 1) nout(5) = 32
    IF (dump_days(2) > 1) nout(6) = 33
    IF (dump_days(3) > 1) nout(7) = 34
    IF (dump_days(4) > 1) nout(8) = 35
  END IF

  IF (grafdump_step) nout(9)  = 37
  IF (grafdump_day)  nout(10) = 38

  IF (grafdump_days) THEN
    IF (dump_days(1) > 1) nout(11) = 39
    IF (dump_days(2) > 1) nout(12) = 40
    IF (dump_days(3) > 1) nout(13) = 41
    IF (dump_days(4) > 1) nout(14) = 42
  END IF

  IF (obs .AND. prindump_obs) THEN
    DO i=1, 5
      nout(i+14) = 42 + i
    END DO
  END IF

!-----------------------------------------------------------------------------
! Write out initial data for run to standard output and to all
! the units to which diagnostics will be written.
!-----------------------------------------------------------------------------

! DEPENDS ON: print_initdata
  CALL print_initdata                                                          &
    ( row_length, rows, land_points, model_levels, wet_levels                  &
    , ozone_levels, nfor, bl_levels, st_levels, sm_levels                      &
    , ntrop, dayno_init, a_sw_radstep_diag, a_sw_radstep_prog, t_inc           &
    , q_star, u_inc, v_inc, w_inc, ichgf, ilscnt, ti, time_init, nout, 19 )

!-----------------------------------------------------------------------------
! Large scale cloud
!-----------------------------------------------------------------------------
! qcl/qcf values initialised in Run_Init, except for when
! geostropic forcing is used.

  cf                    = 0.0
  area_cloud_fraction   = 0.0
  cfl                   = 0.0
  cff                   = 0.0
  qcl                   = 0.0
  qcf                   = 0.0

! Initialise extra microphysics variables to zero as they are not
! set in Run_Init. Forcing options are currently not available for
! these variables.

  qcf2        = 0.0
  qrain       = 0.0
  qgraup      = 0.0
  qcf2_star   = 0.0
  qrain_star  = 0.0
  qgraup_star = 0.0


  daynumber = dayno_init
  year = 1

!=============================================================================
! Options to set initial profiles
!=============================================================================
!     (i)   Observational large scale forcing (OBS=TRUE of
!           Namelist LOGIC)
!           Initial data is then from namelist INPROF
!     (ii)  Statistical large scale forcing (STATS=TRUE of
!           Namelist LOGIC)
!           Initial data can either be derived from climate datasets
!           using subroutine INITSTAT or set from namelist
!           INPROF (set ALTDAT=TRUE in namelist LOGIC)
!     (iii) No large-scale forcing initial data is set from namelist
!           INPROF
!     (iv)  Continuation from previous run stored on tape
!     (Set TAPEIN=TRUE in namelist LOGIC).  All other initial data
!     is overwritten
!=============================================================================

! DEPENDS ON: run_init
  CALL run_init                                                               &
    ! (In)
    ( row_length, rows, model_levels, wet_levels, land_points                 &
    , nfor, bl_levels, st_levels, sm_levels                                   &
    , ntiles, nice, nice_use, ntrop, n_cca_lev, ntab                          &
    , nprimvars, dayno_init, can_model,ichgf                                  &
    , sec_day, rhcrit(1:wet_levels), l_mr_physics2, lcal360, l_snow_albedo    &
    , l_use_sulpc_direct, l_use_seasalt_direct, l_use_soot_direct             &
    , l_use_biogenic, l_use_dust, l_use_bmass_direct, l_use_ocff_direct       &
    , l_use_nitrate_direct, l_use_arclbiom, l_use_arclblck, l_use_arclsslt    &
    , l_use_arclsulp, l_use_arcldust, l_use_arclocff, l_use_arcldlta          &
    , l_murk_rad, l_climat_aerosol, l_use_aod                                 &

    ! (InOut)
    , rho, ti, smcl                                                           &

    ! (Out)
    , iv, iy, idum, iseed, dayno_wint, iccb, icct, cca, cf, cfl, cff, t, q    &
    , qcl, qcf, theta_star, p_star, tstar, u, v, w, w_adv, z0msea, zh         &
    , t_deep_soil, canopy_gb, tsi, smc, sthf, sthu, snodep, catch, infil_tile &
    , z0_tile, z0h_tile, catch_snow                                           &
    , exner_rho_levels, exner_theta_levels, deltap, p                         &
    , p_theta_levels, rp, rp_theta, ch_t_inc, ch_q_star, ch_u_inc, ch_v_inc   &
    , ch_w_inc, ch_flux_e, ch_flux_h, ch_tstar_forcing                        &
    , ch_t_bg, ch_q_bg, ch_u_bg, ch_v_bg, ch_w_bg, ch_ug, ch_vg               &
    , flux_e_scm, flux_h_scm, t_inc_scm                                       &
    , u_inc_scm, v_inc_scm, w_inc_scm, q_star_scm, ug_scm, vg_scm             &
    , t_bg_scm, q_bg_scm, u_bg_scm, v_bg_scm, w_bg_scm                        &
    , tls, qls, uls, vls, wls                                                 &
    , resdump, dap1, dap2, dap3, dab1, dab2, dab3                             &
    , b_exp, hcap, hcon, satcon, sathh, v_sat, v_wilt                         &
    , v_crit, atime, btime, alfada, dbara, dgrada, pa, tbara, tgrada, tsda    &
    , vnbara, vnsda, vpbara, wbara, wsda, alfadb, dbarb, dgradb, pb, tbarb    &
    , tgradb, tsdb, vnbarb, vnsdb, vpbarb, wbarb, wsdb )


  ! Should be set in JULES module
  rho_snow_grnd = i_rho_snow_grnd

!==========================================================
! cf is initialised in run_init,
! area_cloud_fraction initialised to the same value as cf (bulk fraction)
!==========================================================
  DO i=1, row_length
    DO j=1, rows
      DO  k=1, wet_levels
        area_cloud_fraction(i,j,k) = cf(i,j,k)
      END DO  ! k
    END DO  ! j
  END DO  ! i


! Initialise all surf temperature variables to TSTARI from
! namelist &INPROF
  DO i=1, row_length
    DO j=1, rows
      IF (tstar_sea(i,j) == RMDI) THEN
        tstar_sea(i,j)  = tstar(i,j)
      END IF

      IF (tstar_sice(i,j) == RMDI) THEN
        tstar_sice(i,j) = tstar(i,j)
      END IF

      IF (tstar_land(i,j) == RMDI) THEN
        tstar_land(i,j) = tstar(i,j)
      END IF
    END DO
  END DO

  DO i=1, nlnd
    DO itype=1, ntype
      IF (tstar_tile(i,itype) == RMDI) THEN
        tstar_tile(i,itype) = tstar(i,j)
      END IF
    END DO
  END DO

!-----------------------------------------------------------------------------
! For geostrophic forcing
!-----------------------------------------------------------------------------

  DO j=1, rows
    DO i=1, row_length
      f_coriolis(i,j) = 2.0 * omega * SIN(lat(i,j) * pi_over_180)
      f3_at_u(i,j) = two_omega * sin_theta_latitude(i,j)
    END DO
  END DO

  IF (geoinit .AND. geoforce) THEN

    kill_scmdiags(:) = .FALSE.

    DO i=1, row_length
      DO j=1, rows

        modug(i,j) = 0.0
        DO k=1, model_levels
          modug(i,j) = MAX( modug(i,j)                                   &
                          , SQRT(  ug_scm(i,j,k)*ug_scm(i,j,k)           &
                                 + vg_scm(i,j,k)*vg_scm(i,j,k)) )
          u(i,j,k) = ug_scm(i,j,k)
          v(i,j,k) = vg_scm(i,j,k)
        END DO

      END DO
    END DO

!-----------------------------------------------------------------------------
! Form restart dump
!-----------------------------------------------------------------------------

    resdump(:,:,:) = 0.0

! DEPENDS ON: restart_dump
    CALL restart_dump                                                         &
      ( row_length, rows, model_levels, wet_levels, nprimvars                 &
      , land_points, bl_levels, st_levels, sm_levels                          &
      , n_cca_lev, land_sea_mask, resdump, u, v, w, t, theta                  &
      , q, qcl, qcf, cf, p, rho, t_deep_soil, smc, canopy_gb                  &
      , snodep, tstar, zh, z0msea, cca, iccb, icct, smcl )

    DO i=1, row_length
      DO j=1, rows
        rccb(i,j) = REAL(iccb(i,j))
        rcct(i,j) = REAL(icct(i,j))
      END DO
    END DO

    timestep_number = 0
    maxinc = 9999.0

    DO WHILE (maxinc  >   0.1                                                 &
      .AND.  timestep_number  <   full_daysteps)

      timestep_number = timestep_number + 1
      daycount  = 1

      ! DEPENDS ON: timecalc
      CALL timecalc                                                           &
        ( dayno_init, time_init, time_string, lcal360, yearno, day, time_sec  &
        , previous_time, ihour, imin )

      time_info = previous_time

      ! Call the pre-physics routine to set up variables
      ! DEPENDS ON: pre_physics
      CALL pre_physics                                                        &
        ! (In)
        ( row_length, rows, model_levels, wet_levels, nfor, ichgf, qcl, qcf   &
        , ch_ug, ch_vg, ilscnt, f_coriolis, lcal360, daycount                 &
        , timestep_number, a_sw_radstep_diag, a_sw_radstep_prog               &
        , l_triffid, npft                                                     &
        ! (InOut)
        , u, v, ug_scm, vg_scm, npft_trif                                     &
        ! (Out)
        , co2_mmr, l_rad_step, L_rad_step_prog, l_rad_step_diag )

      l_rad_step      = .FALSE.
      l_rad_step_prog = .FALSE.
      l_rad_step_diag = .FALSE.

!         Now call the physics routines with flags set to
!         only enable the boundary layer call.

! DEPENDS ON: atmos_physics2
      CALL atmos_physics2                                                     &
! Parallel variables
        ( row_length, rows, n_proc, n_procx, n_procy, rows, row_length        &
        , numcycles, cycleno                                                  &

! Model dimensions.
        , row_length, rows, rows, land_points, model_levels, nice, nice_use   &
        , wet_levels                                                          &
        , bl_levels, st_levels, sm_levels, cloud_levels, land_ice_points      &
        , soil_points, n_cca_lev, ntiles, tr_levels                           &
        , first_constant_r_rho_level, dim_cs1, dim_cs2                        &

! Model switches
        , l_regular, l_mr_physics2, l_dry, l_lbc_old                          &
        , lcal360, ltimer                                                     &
        , l_cld_area                                                          &

! Model Parameters
        , rhcrit(1:wet_levels), co2_mmr                                       &
        , tr_vars, tr_ukca                                                    &

! In: coordinate information
        , unscl_dry_rho                                                       &
        , delta_lambda, delta_phi, dlambda_p, dphi_p                          &
        , wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v                        &
        , lat_rot_np, long_rot_np, f3_at_u                                    &

! In: time stepping information.
        , yearno, day, ihour, imin, time_sec                                  &

! River routing
        , row_length, row_length, xpa, xua, xva, ypa, yua, yva                &
        , g_p_field, g_r_field, a_steps_since_riv, river_row_length           &
        , river_rows, global_river_row_length, global_river_rows              &
        , river_vel, river_mcoef, i_river_vn                                  &
        , trivdir, trivseq, twatstor, inlandout_atm                           &

! Lake evaporation:
        , acc_lake_evap                                                       &

! Grid-to-grid river routing
        , r_area, slope, flowobs1, r_inext, r_jnext, r_land                   &
        , substore, surfstore, flowin, bflowin,                               &

! diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
          stashwork3, stashwork5, stashwork8, stashwork9, stashwork19         &
        , stashwork26                                                         &

! SCM diagnostics
        , nscmdpkgs, kill_scmdiags                                            &

! In: data fields.
        , theta, q, qcl, qcf, qrain, qgraup, qcf2                             &
        , rho, u, v, w, w_adv, p, p_star                                      &
        , exner_rho_levels, exner_theta_levels, land_sea_mask, p_theta_levels &

! Ancillary fields and fields needed to be kept from timestep to
! timestep
        , land_index, land_ice_index, soil_index, canopy_gb, snodep, hcon     &
        , hcap, v_crit, v_wilt, v_sat, sthf, sthu, sil_orog_land              &
        , ho2r2_orog, sd_orog_land, di, ice_fract, u_0, v_0, u_0, v_0         &
        , cca, iccb, icct, cclwp, ccw, lcbase                                 &
        , t_deep_soil, tsi, ti_n, ice_k_n, tstar                              &
        , z0msea, ice_fract_n, ice_thick_n, satcon, sathh, b_exp, smcl        &
        , t1_sd, q1_sd, zh, ddmfx, area_cloud_fraction, cf                    &
        , cfl, cff, ls_rain, ls_snow                                          &
        , micro_tends, photosynth_act_rad, rad_hr                             &
        , soil_clay, soil_silt, soil_sand, dust_mrel1, dust_mrel2, dust_mrel3 &
        , dust_mrel4, dust_mrel5, dust_mrel6                                  &
        , so2_high_level, so2_em, nh3_em, dms_em, soot_hilem                  &
        , soot_em, ocff_hilem, ocff_em, co2_emits, co2flux                    &
        , deep_flag, past_precip, past_conv_ht                                &

! In/Out
        , theta_star                                                          &
        , q_star_scm, qcl_star, qcf_star, qrain_star, qgraup_star, qcf2_star  &
        , cf_star, cfl_star, cff_star                                         &
        , u_inc_scm, v_inc_scm, w_inc_scm(:,:,1:model_levels)                 &
        , sum_eng_fluxes, sum_moist_flux                                      &

! Tracer fields(InOut)
       , aerosol, free_tracers, ukca_tracers                                  &
       , dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6     &
       , so2, dms, so4_aitken, so4_accu, so4_diss, nh3, soot_new              &
       , soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld                &
       , ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss, co2              &

! RIVERS (InOut)
       , tot_surf_runoff, tot_sub_runoff                                      &

! Out: fields
       , rhokm, cH_term, ntml, cumulus, nbdsc, ntdsc                          &
       , rhcpt,row_length, rows                                               &

! Additional variables for Land Surface Scheme
       , frac_typ, frac_disturb, canht, lai                                   &
       , canopy(1:land_points,1:ntiles)                                       &
       , catch(1:land_points,1:ntiles), catch_snow, snow_grnd, snow_tile      &
       , z0_tile, z0h_tile, tstar_tile                                        &
       , infil_tile(1:land_points,1:ntiles)                                   &
       , rgrain(1:land_points,1:ntiles)                                       &
       , cs, gs, co2_dim_row, co2_dim_len                                     &
       , asteps_since_triffid, timestep_number                                &
       , g_leaf_acc                                                           &
       , g_leaf_phen_acc, npp_ft_acc, resp_w_ft, resp_s_acc                   &
       , land_pts_trif, npft_trif, dolr, lw_down, sw_tile, fland_ctile        &
       , tstar_land, tstar_sea, tstar_sice_n, tstar_sice                      &
       , albsoil, cos_zenith_angle                                            &

! INOUT variables for TKE based turbulence schemes
       , e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu                         &

! Additional variables required for large-scale hydrology:
       , fexp, gamtot, ti_mean, ti_sig, fsat, fwetl, zw                       &
       , sthzw, a_fsat, c_fsat, a_fwet, c_fwet                                &

! JULES 2 prognostics (InOut)
       , snowdepth, rho_snow_grnd                                             &
       , nsnow, ds, sice, sliq, tsnowlayer, rho_snow, rgrainl                 &

! FLake lake scheme prognostics (InOut)
       , lake_depth, lake_fetch, lake_t_mean, lake_t_mxl                      &
       , lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape                      &
       , lake_g_dt                                                            &

! Cariolle ozone
       , ozone_tracer                                                         &

! Additional screen-level variables
       , tscrndcl_ssi, tscrndcl_tile, tstbtrans                               &

! Variables required for COSP (out)
        , cosp_crain_3d, cosp_csnow_3d                                        &

! Prescribed surface forcing and roughness lengths
       , l_flux_bc, flux_e_scm, flux_h_scm, l_spec_z0, z0m_scm, z0h_scm       &

! Additional variables for SCM diagnostics
       , nfor, ntrop, time_string, daycount, cf                               &
       , nconvars, ntrad1, conv_mode, .FALSE.                                 &
! Error information
       , error_code )


! Calculate max increment and copy winds for safe-keeping
      maxinc = 0.0
      DO i=1, row_length
        DO j=1, rows
          DO k=1, model_levels
            maxinc = MAX(maxinc,(u_inc_scm(i,j,k)**2 +                        &
                                 v_inc_scm(i,j,k)**2))
            ui(i,j,k) = u(i,j,k)
            vi(i,j,k) = v(i,j,k)
          END DO
          maxinc = SQRT(maxinc)/(f_coriolis(i,j) * timestep*modug(i,j))
        END DO
      END DO


!-----------------------------------------------------------------------------
! Copy initial data back from DUMP.
!-----------------------------------------------------------------------------

! DEPENDS ON: dumpinit
      CALL dumpinit                                                           &
        ( row_length, rows, nprimvars, land_points, model_levels              &
        , wet_levels, bl_levels, st_levels, sm_levels, ntrop, n_cca_lev       &
        , land_sea_mask, resdump, u, v, w, t, theta, q, qcl, qcf              &
        , cf, p, rho, t_deep_soil, smc, canopy_gb, snodep, tstar, zh          &
        , z0msea, cca, rccb, rcct, smcl )

      DO i=1, row_length
        DO j=1, rows
          iccb(i,j) = INT(rccb(i,j))
          icct(i,j) = INT(rcct(i,j))

! Copy saved U,V back

          DO k=1, model_levels
            u(i,j,k) = ui(i,j,k)
            v(i,j,k) = vi(i,j,k)
          END DO
        END DO
      END DO
    END DO                   ! maxinc and timestep_number < daysteps

    DO i=1, row_length
      DO j=1, rows
        DO k=1, nprimvars
          resdump(i,j,k) = 0.0
        END DO
      END DO
    END DO

    WRITE (6,*) "Geostrophic wind initialised."
    WRITE (6,*) "max relative wind change at end=",maxinc,       &
      " after ",timestep_number," steps"
    daynumber = dayno_init
    year = 1

  END IF                     ! Geostrophic forcing initialising.

!-----------------------------------------------------------------------------
! Initialise the output diagnostic system
!-----------------------------------------------------------------------------
! DEPENDS ON: setup_diags
  CALL setup_diags                                               &
    ( row_length, rows, model_levels                             &! In
    , wet_levels, bl_levels, sm_levels                           &! In
    , st_levels, land_points, ntiles, n_vis_thresh, cloud_levels &! In
    , total_nsteps, timestep, full_daysteps                      &! In
    , a_sw_radstep_prog, a_sw_radstep_diag, ntrad1               &! In
    , daycount, timestep_number, nscmdpkgs                       &! In
    , l_scmdiags, scmop)                ! Out

! If PC2 is off then l_SCMDiags(SCMDiag_pc2) must be false
  IF (.NOT. l_pc2) THEN
    IF (l_scmdiags(scmdiag_pc2)) THEN
      WRITE(6,*) '***********************************'
      WRITE(6,*) '* Warning: you have requested PC2 *'
      WRITE(6,*) '* diagnostics but PC2 is off,     *'
      WRITE(6,*) '* resetting l_SCMDiag_PC2 logical *'
      WRITE(6,*) '***********************************'
    END IF
    l_scmdiags(scmdiag_pc2) = .FALSE.
  END IF

! The availability of Surface based diagnostics packages is determined
! by the surface type.
!
! Surface Package     (SCMDiags_surf) - always available
! Land Points Package (SCMDiags_land) - Only if land_sea_mask TRUE
! Sea Points Package  (SCMDiags_sea)  - Only if land_sea_mask FALSE
!
! This only works in the SCM because rows, row_length are size 1 here

  DO j=1, rows
    DO i=1, row_length
      IF (land_sea_mask(i,j)) THEN

        ! Land point
        IF (l_scmdiags(scmdiag_sea)) THEN
          WRITE(6,*) '********************************************'
          WRITE(6,*) '* Warning: you have requested sea          *'
          WRITE(6,*) '* diagnostics but this is a land point,    *'
          WRITE(6,*) '* resetting l_SCMDiag_sea logical to false *'
          WRITE(6,*) '********************************************'
        END IF

        l_scmdiags(scmdiag_sea) = .FALSE.

      ELSE

        ! Sea point
        IF (l_scmdiags(scmdiag_land)) THEN
          WRITE(6,*)'*********************************************'
          WRITE(6,*)'* Warning: you have requested land          *'
          WRITE(6,*)'* diagnostics but this is a sea point,      *'
          WRITE(6,*)'* resetting l_SCMDiag_land logical to false *'
          WRITE(6,*)'*********************************************'
        END IF

        l_scmdiags(scmdiag_land) = .FALSE.

      END IF ! land_sea_mask
    END DO
  END DO

!-----------------------------------------------------------------------------
! LOOP OVER DAYS
!-----------------------------------------------------------------------------

! Timestepping proper is about to begin, switch the
! diagnostic system on (as long as at least one stream is open)

  IF (main_diag_switch /= 0 .AND. ANY(scmop%strm%switch /= 0)) THEN
        ! any() is a Fortran90 intrinsic funtion
    scmop%on = .TRUE.
  END IF


  scm_timestep_count = 0

  DO daycount=1, ndayin+1

!-----------------------------------------------------------------------------
! Reset sinusoidal distribution every 10 days if climate stats required
!-----------------------------------------------------------------------------

    IF (stats .AND. (daycount  ==  1                                          &
      .OR. (ancyc .AND. MOD(daycount, change_clim)  ==  0))) THEN
! DEPENDS ON: statday
      CALL StatDay                                                            &
        ! (In)
        ( row_length, rows, model_levels, wet_levels, ntrop                   &
        , atime, btime, dayno_wint, deltan, daycount                          &
        , tbara, tbarb, tsda, tsdb, dbara, dbarb, vnbara, vnbarb              &
        , vnsda, vnsdb, vpbara, vpbarb, wbara, wbarb, wsda, wsdb              &
        , alfada, alfadb, pa, pb, p, tgrada                                   &
        , tgradb, dgrada, dgradb, cort, cord, corvn, corw                     &
        , tdash, ddash, ctbar, ctsd, at, cdbar, cdsd, ad                      &
        , cvnbar, cvnsd, avn, cwbar, cwsd, aw                                 &
        , tbar, tsd, dbar, dsd                                                &
        , vnbar, vnsd, vpbar, wbar, wsd )

!-----------------------------------------------------------------------------
!     Calculate values of exner
!-----------------------------------------------------------------------------

      l_calc_exner = .TRUE.
      l_calc_rho   = .FALSE.

! DEPENDS ON: calc_press
      CALL calc_press                                                         &
        ! (In)
        ( model_levels, wet_levels, rows, row_length, p                       &
        , theta, q, l_calc_exner, l_calc_rho                                  &
        ! (InOut)
        , rho                                                                 &
        ! (Out)
        , exner_theta_levels, exner_rho_levels, p_theta_levels                &
        , rp, rp_theta, p_star)

!-----------------------------------------------------------------------------
! Initialise PX and PY arrays for calculation of vertical fluxes later
!-----------------------------------------------------------------------------

      DO i=1, row_length
        DO j=1, rows
          DO k=1, ntrop
            px(i,j,k) = 1.0 / alog(p(i,j,k+1)/ p(i,j,k))
          END DO

          DO k=1, ntrop-1
            py(i,j,k) = 1.0 / alog(p(i,j,k+2)/ p(i,j,k))
          END DO
        END DO
      END DO
    END IF                   ! stats

!=============================================================================
! Options for forcing
!-----------------------------------------------------------------------------
!
!       Observational forcing : use observed values of T,Q,U,V
!       to calculate their increments due to large scale effects.
! OR
!       Statistical forcing : take random sample from Normal
!       (Gaussian) distribution with mean and sd climlogical
!       average to calculate increments to T and Q due to large scale
!       effects.
!=============================================================================
!
!       Loop over timesteps
!
!
!       If it is the last day in the run and a full day is not
!       required, loop over the number of timesteps required
!       otherwise Do the full number of timesteps in a day.

    IF (daycount  ==  ndayin+1 .AND. nstepsin  /=  full_daysteps) THEN
      daysteps = nstepsin
    ELSE
      daysteps = full_daysteps
    END IF

    DO timestep_number=1, daysteps

      ! Update local_timestep_count
      scm_timestep_count = scm_timestep_count + 1

!-----------------------------------------------------------------------------
! VARIABLE VALUES: q     is vapour         at timelevel n
!                  qcl   is liquid         at timelevel n
!                  t     is temperature    at timelevel n
!                  theta is potential temp at timelevel n
!-----------------------------------------------------------------------------

      IF (main_diag_switch /= 0) THEN
!
!-----------------------------------------------------------------------------
! SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------------
        IF (l_scmdiags(scmdiag_pc2)) THEN

! DEPENDS ON: scmoutput
          CALL scmoutput(q, 'q_timen'                                         &
            , 'Vapour at start of timestep', 'kg/kg'                          &
            , t_inst, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(qcl, 'qcl_timen'                                     &
            , 'Liquid at start of timestep', 'kg/kg'                          &
            , t_inst, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(t, 't_timen'                                         &
            , 'Temperature at start of timestep', 'K'                         &
            , t_inst, d_all, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(theta, 'th_timen'                                    &
            , 'Potential temperature at start of timestep', 'K'               &
            , t_inst, d_all, default_streams, '', routinename)

        END IF ! l_SCMDiags(SCMDiag_PC2)

      END IF ! main_diag_switch /= 0

! Reset the wind increments to zero -
! The other increments are reset in physics1

      u_inc_scm(:,:,:) = 0.0
      v_inc_scm(:,:,:) = 0.0
      w_inc_scm(:,:,:) = 0.0

!-----------------------------------------------------------------------------
! Calculate the year (in run) and actual time and day number
! for labelling  of the diagnostics only.
!-----------------------------------------------------------------------------

! DEPENDS ON: timecalc
      CALL timecalc                                                           &
        ( dayno_init, time_init, time_string, lcal360, yearno, day, time_sec  &
        , previous_time, ihour, imin)

      time_info = previous_time

! If there is no annual cycle, the year and day numbers
! will be the init ones.
      IF (.NOT. ancyc) THEN
        year = 1
        day = dayno_init
      END IF

      !-------------------------------------------------------------
      ! Convert temperature to potential temperature
      !         t_inc_scm   to theta_star
      !-------------------------------------------------------------
      theta(:,:,:)      = t(:,:,:)         / exner_theta_levels(:,:,:)
      theta_star(:,:,:) = t_inc_scm(:,:,:) / exner_theta_levels(:,:,:)


      IF (l_pc2) THEN

!-----------------------------------------------------------------------------
! Calculate initial relative humidity w.r.t. Liquid temperature TL
! Is used by PC2 routines at end of timestep. Jeremy Price Feb 2005.
!-----------------------------------------------------------------------------
        DO k=1, wet_levels

          DO j=1, rows
            DO i=1, row_length
              tl(i,j) = t(i,j,k) - (lc/cp)*qcl(i,j,k)
            END DO ! i
          END DO ! j

! DEPENDS ON: qsat_wat
          CALL qsat_wat (qsl_tl, tl, p_theta_levels(1,1,k), row_length*rows )

          DO j=1, rows
            DO i=1, row_length
              rhts(i,j,k) = (q(i,j,k) + qcl(i,j,k))/qsl_tl(i,j)
              tlts(i,j,k) = tl(i,j)
              ptts(i,j,k) = p_theta_levels(i,j,k)
              qtts(i,j,k) = q(i,j,k) + qcl(i,j,k)
            END DO ! i
          END DO ! j

        END DO ! k

      END IF  ! l_pc2


!-----------------------------------------------------------------------------
! Save initial fields for later to calculate total increments
!-----------------------------------------------------------------------------
      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length
            t_start(i,j,k) = t(i,j,k)
            u_start(i,j,k) = u(i,j,k)
            v_start(i,j,k) = v(i,j,k)
          END DO
        END DO
      END DO

      DO k=0, model_levels
        DO j=1, rows
          DO i=1, row_length
            w_start(i,j,k) = w(i,j,k)
          END DO
        END DO
      END DO

      DO k=1, wet_levels
        DO j=1, rows
          DO i=1, row_length
            q_start(i,j,k)   = q(i,j,k)
            qcl_start(i,j,k) = qcl(i,j,k)
            qcf_start(i,j,k) = qcf(i,j,k)
            cf_start(i,j,k)  = cf(i,j,k)
            cfl_start(i,j,k) = cfl(i,j,k)
            cff_start(i,j,k) = cff(i,j,k)
          END DO
        END DO
      END DO

!-----------------------------------------------------------------------------
! Call physics one before advection step
!-----------------------------------------------------------------------------

! DEPENDS ON: pre_physics
      CALL pre_physics                                                        &
        ! (In)
        ( row_length, rows, model_levels, wet_levels, nfor, ichgf, qcl, qcf   &
        , ch_ug, ch_vg, ilscnt, f_coriolis, lcal360, daycount                 &
        , timestep_number, a_sw_radstep_diag, a_sw_radstep_prog               &
        , l_triffid, npft                                                     &
        ! (InOut)
        , u, v, ug_scm, vg_scm, npft_trif                                     &
        ! (Out)
        , co2_mmr, l_rad_step, l_rad_step_prog, l_rad_step_diag )

! Set flags for forcing to false

      l_forcing = .FALSE.

!------------------------------------------------------------------|
! UM5.x timestepping stores all wind increments in u_inc_scm and   |
! v_inc_scm. These are then added to u,v in atmos_physics2.        |
!------------------------------------------------------------------|

      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length

            !---------------------------------------------------|
            ! Calculate increments from geostrophic forcing     |
            ! These are then added to u_inc_scm and v_inc_scm   |
            ! before atmos_physics2 as atmos_physics1 sets      |
            ! u_inc_scm and v_inc_scm to zero.                  |
            !---------------------------------------------------|

            uinc_geo(i,j,k) = u(i,j,k)-u_start(i,j,k)
            vinc_geo(i,j,k) = v(i,j,k)-v_start(i,j,k)

            !--------------------------------|
            ! Reset u,v back to time-level n |
            !--------------------------------|

            u(i,j,k) = u_start(i,j,k)
            v(i,j,k) = v_start(i,j,k)

          END DO  ! i
        END DO  ! j
      END DO  ! k

      !-----------------------------------
      ! Convert to mixing ratios if needed
      !-----------------------------------
      IF (l_mr_physics1) THEN

        WRITE(6,*) 'Convert to mixing ratios'

        ! Convert to mixing ratios
! DEPENDS ON: q_to_mix
        CALL q_to_mix                                                         &
          ( row_length, rows, wet_levels, halo_i, halo_j                      &
          , q, qcl, qcf, qcf2, qrain, qgraup                                  &
          , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                             &
          , mix_v, mix_cl, mix_cf, mix_cf2, mix_rain, mix_graup )

        ! Now place mixing ratio values back into d1
        DO k=1, wet_levels
          DO j=1, rows
            DO i=1, row_length

              q(i,j,k)     = mix_v(i,j,k)     ! Vapour
              qcl(i,j,k)   = mix_cl(i,j,k)    ! Liquid
              qcf(i,j,k)   = mix_cf(i,j,k)    ! Ice
              qcf2(i,j,k)  = mix_cf2(i,j,k)   ! Ice 2
              qrain(i,j,k) = mix_rain(i,j,k)  ! Rain
              qgraup(i,j,k)= mix_graup(i,j,k) ! Graupel

            END DO
          END DO
        END DO

      END IF  ! l_mr_physics1

!=============================================================================

! DEPENDS ON: atmos_physics1
      CALL atmos_physics1                                                     &
! Parallel variables
        ( row_length, rows, n_proc, n_procx, n_procy, rows, row_length        &

! model dimensions.
        , row_length, rows, rows, land_points, model_levels                   &
        , wet_levels, bl_levels, st_levels, sm_levels                         &
        , ozone_levels, cloud_levels, land_ice_points, soil_points            &
        , n_cca_lev, ntiles, nice_use, salt_dim1, salt_dim2, salt_dim3        &
        , tr_levels, tr_ukca, cdnc_dim1, cdnc_dim2, cdnc_dim3                 &
        , co2_dim_len, co2_dim_row, co2_dim_lev                               &
        , n_arcl_species, n_arcl_compnts, i_arcl_compnts                      &

! Model switches
        , l_regular, l_lbc_old                                                &
        , l_rad_step, l_rad_step_diag, l_rad_step_prog                        &
        , lcal360, ltimer                                                     &

! Setting lflux_reset to false but, if using energy correction, will
! need to be worked out every timestep
        , l_mr_physics1                                                       &
        , L_ukca_chem, l_ukca_useumuivals, l_ukca_set_trace_gases             &
        , l_ukca_strat, l_ukca_strattrop                                      &
        , l_ukca_prescribech4, l_use_arcl                                     &

! Model Parameters
        , rhcrit(1:wet_levels)                                                &
        , min_trop_level, max_trop_level                                      &
        , ngrgas, grgas_addr                                                  &
        , Ntot_land, Ntot_sea                                                 &

! Parameter for stochastic physics random parameters2
        , m_ci                                                                &

! In: coordinate information
        , delta_lambda, delta_phi, lat_rot_NP, long_rot_NP                    &

! In: time stepping information.
        , yearno, day, ihour, imin, isec, previous_time,                      &

! Diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
          stashwork1, stashwork2, stashwork4, stashwork6, stashwork14         &

! SCM diagnostics
        , nscmdpkgs, l_scmdiags                                               &

! In data fields.
        , theta, q, qcl, qcf, qcf2, qrain, qgraup, rho,  u, v, p              &
        , p_star, exner_rho_levels, exner_theta_levels, land_sea_mask         &
        , p_theta_levels, fland_ctile, frac_control, cdnc_ukca_dummy          &

! ancillary fields and fields needed to be kept from timestep to
! timestep
        , land_index, rgrain(1:land_points,1:ntiles), soot, ntml, cumulus     &
        , ice_fract_n, ice_fract_n, ice_thick_n                               &
        , cca, iccb, icct, cclwp, ccw, lcbase                                 &
        , tstar, tstar_land, tstar_sea, tstar_sice_n                          &
        , sice_alb, land_alb, snodep, snodep_sice_n                           &
        , ozone, sw_incs, lw_incs                                             &
        , dirpar_incs, o3_trop_level, o3_trop_height, t_trop_level            &
        , t_trop_height, zh, sd_orog_land, orog_grad_xx_land                  &
        , orog_grad_xy_land, orog_grad_yy_land, area_cloud_fraction           &
        , cf, cfl, cff, aerosol_em, arcl, albsoil, albobs_sw, albobs_vis      &
        , albobs_nir, lai, snow_tile, frac_typ, tstar_tile, z0_tile           &
        , dolr_rts, lw_down, sw_tile_rts, es_space_interp, rad_mask           &
        , cos_zenith_angle                                                    &

! Variables for COSP
        , cosp_crain_3d,cosp_csnow_3d                                         &

! JULES 2 prognostics (In)
        , snowdepth,lake_h_ice                                                &

! InOut
        , theta_star                                                          &
        , q_star_scm, qcl_star, qcf_star, qcf2_star, qrain_star, qgraup_star  &
        , cf_star, cfl_star, cff_star, u_inc_scm, v_inc_scm                   &
        , energy_correction, sum_eng_fluxes, sum_moist_flux, aerosol          &
        , dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6    &
        , so2, so4_aitken, so4_accu, so4_diss, nh3                            &
        , soot_new, soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld     &
        , ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss                  &
        , co2, ukca_tracers, biogenic, asteps_since_triffid, ukca_radaer      &

! Out
        , ls_rain, ls_snow, micro_tends, unscl_dry_rho                        &
        , photosynth_act_rad, rad_hr, dolr, sw_tile                           &

! Error information
        , error_code )

      !-------------------------------
      ! Convert to specific humidities
      !-------------------------------

      IF (l_mr_physics1) THEN
        WRITE(6,*) 'In l_mr_physics1 branch of reconversion'

        ! Copy q and qstar variables to mix variables
        mix_v_star(:,:,:)     = q_star_scm(:,:,:)  ! Vapour
        mix_cl_star(:,:,:)    = qcl_star(:,:,:)    ! Liquid
        mix_cf_star(:,:,:)    = qcf_star(:,:,:)    ! Ice
        mix_cf2_star(:,:,:)   = qcf2_star(:,:,:)   ! Ice2
        mix_rain_star(:,:,:)  = qrain_star(:,:,:)  ! Rain
        mix_graup_star(:,:,:) = qgraup_star(:,:,:) ! Graupel


        ! Convert mixing ratios back to specific humidities
! DEPENDS ON: mix_to_q
        CALL mix_to_q                                                         &
          ( row_length, rows, wet_levels                                      &
          , halo_i, halo_j, mix_v, mix_cl, mix_cf                             &
          , mix_cf2, mix_rain, mix_graup                                      &
          , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                             &
          , q, qcl, qcf, qcf2, qrain, qgraup )

        ! Convert mixing ratio increments (mix_star) back
        ! to specific humidity increments (q_star_scm)
! DEPENDS ON: calc_q_star
        CALL calc_q_star                                                      &
          ( row_length, rows, wet_levels                                      &
          , halo_i, halo_j, offx, offy                                        &
          , mix_v, mix_cl, mix_cf                                             &
          , mix_cf2, mix_rain, mix_graup                                      &
          , mix_v_star, mix_cl_star, mix_cf_star                              &
          , mix_cf2_star, mix_rain_star, mix_graup_star                       &
          , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                             &
          , q, qcl, qcf, qcf2, qrain, qgraup                                  &
          , q_star_scm, qcl_star, qcf_star                                    &
          , qcf2_star, qrain_star, qgraup_star )

      END IF  ! l_mr_physics1


!------------------------------------------------------------------|
! Add on increments from geostrophic wind forcing.                 |
!------------------------------------------------------------------|

      IF (geoforce) THEN
        u_inc_scm(:,:,:) = u_inc_scm(:,:,:) + uinc_geo(:,:,:)
        v_inc_scm(:,:,:) = v_inc_scm(:,:,:) + vinc_geo(:,:,:)
      END IF

!-----------------------------------------------------------------------------
! Convert theta_star to t_inc_scm for call to forcing
!-----------------------------------------------------------------------------
      t_inc_scm(:,:,:) = theta_star(:,:,:) * exner_theta_levels(:,:,:)

      ! set qcl_inc and qcf_inc to increments
      qcl_inc(:,:,:) = qcl_star(:,:,:)
      qcf_inc(:,:,:) = qcf_star(:,:,:)


!-----------------------------------------------------------------------------
! VARIABLE VALUES:
!   qcl_inc and qcl_star = liq. water incs      (atmos_physics1)
!   q_star_scm           = vapour incs          (atmos_physics1)
!   t_inc_scm            = temp. incs           (atmos_physics1)
!   theta_star           = pot. temp.incs       (atmos_physics1)
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! PC2 section.
! Store increments from atmos_physics1. We can then work out the forcing
! by subtraction of these _earliest values from the values returned after
! forcing.
!-----------------------------------------------------------------------------

      IF (l_pc2) THEN
        q_earliest   (:,:,:) = q_star_scm (:,:,:)
        qcl_earliest (:,:,:) = qcl_inc    (:,:,:)
        t_earliest   (:,:,:) = t_inc_scm  (:,:,:)
      END IF  ! l_pc2

      IF (main_diag_switch /= 0) THEN

!-----------------------------------------------------------------------------
! SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------------
        IF (l_scmdiags(scmdiag_pc2)) THEN

! DEPENDS ON: scmoutput
          CALL scmoutput(q_earliest, 'dq_earliest'                            &
            , 'q_star vapour incs from atmos_phys1', 'kg/kg'                  &
            , t_inst, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(qcl_earliest, 'dqcl_earliest'                        &
            , 'qcl_inc liq water incs atmos_phys1', 'kg/kg'                   &
            , t_inst, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(t_earliest, 'dt_earliest'                            &
            , 't_inc temp incs from atmos_phys1', 'K'                         &
            , t_inst, d_all, default_streams, '', routinename)

        END IF ! l_SCMDiags(SCMDiag_PC2)

      END IF ! main_diag_switch /= 0

!-----------------------------------------------------------------------------
! Include any forcing required
!-----------------------------------------------------------------------------

      IF (stats .OR. obs) THEN
! DEPENDS ON: forcing
        CALL forcing                                                          &
          ! (In)
          ( row_length, rows, model_levels, wet_levels, nfor, bl_levels       &
          , st_levels, sm_levels, ntrop, sec_day, timestep_number, daycount   &
          , dayno_wint, daysteps, nscmdpkgs, ntab, ichgf, t, q, qcl, qcf      &
          , u, v, w, l_scmdiags, p, exner_theta_levels, rp, r_theta_levels    &
          , ch_tstar_forcing, ch_flux_h, ch_flux_e                            &
          , ch_t_inc, ch_q_star, ch_u_inc, ch_v_inc, ch_w_inc                 &
          , ch_t_bg,  ch_q_bg,   ch_u_bg,  ch_v_bg,  ch_w_bg                  &
          , ad, at, avn, aw, cdbar, ctbar, cvnbar, cwbar, dbar                &
          , tbar, vnbar, wbar, vpbar, cdsd, ctsd, cvnsd, cwsd, dsd, tsd, vnsd &
          , wsd, tdash, ddash, deltan, px, py                                 &

          ! (InOut)
          , ilscnt, flux_h_scm, flux_e_scm, ti, qi, t_inc_scm                 &
          , q_star_scm, qcl_inc, qcf_inc, u_inc_scm, v_inc_scm, w_inc_scm     &
          , t_bg_scm, q_bg_scm, u_bg_scm, v_bg_scm, w_bg_scm                  &
          , tls, qls, uls, vls, wls, tr, qr, vnr, vpr, wr                     &

          ! (Out)
          , iv, iy, idum, tstar, factor_rhokh, rhokh, dab1, dap1 )
      END IF

!-----------------------------------------------------------------------------
! VARIABLE VALUES:
!   qcl_inc    = liquid water increments        (atmos_physics1)
!              + forcing increments.
!   t_inc_scm  = temperature increments         (atmos_physics1)
!              + forcing increments
!   q_star_scm = vapour increments              (atmos_physics1)
!              + forcing increments
!   wls        = The current vertical velocity forcing tendency (m/s)/day
!----------------------------------------------------------------------------
! Update model vertical velocities to be consistent with
! subsidence forcing
!----------------------------------------------------------------------------

      ! Want latest w whether vertical advection on or off.
      ! w = Current w from start of timestep
      w(:,:,:) = w(:,:,:) + w_inc_scm(:,:,:)
      w_adv(:,:,:) = w(:,:,:)


!-----------------------------------------------------------------------------
! Update SST to be consistent with tstar_forcing (Coastal tiling)
!-----------------------------------------------------------------------------
      DO j=1, rows
        DO i=1, row_length
          tstar_sea(i,j) = tstar(i,j)
        END DO ! i
      END DO ! j

!-----------------------------------------------------------------------------
! This is a PC2 section.
! Calculate the forcing of vapour, liquid, temperature and pressure.
!-----------------------------------------------------------------------------

      IF (l_pc2) THEN
        q_forcing(:,:,:)   = q_star_scm(:,:,:) - q_earliest(:,:,:)
        qcl_forcing(:,:,:) = qcl_inc(:,:,:)    - qcl_earliest(:,:,:)
        t_forcing(:,:,:)   = t_inc_scm(:,:,:)  - t_earliest(:,:,:)
        p_forcing(:,:,:)   = 0.0
      END IF  ! l_pc2

!-----------------------------------------------------------------------------
! Increments need to be converted to increment+value for t, q, qcl and qcf
!-----------------------------------------------------------------------------

!  Before going into physics2, t_inc_scm, q_inc, qcl_inc and qcf_inc are
!  converted to increment plus value whilst the winds stay as increments

      t_inc_scm(:,:,:) = t_inc_scm(:,:,:) + t(:,:,:)

      DO k=1, wet_levels
        DO j=1, rows
          DO i=1, row_length
            q_star_scm(i,j,k) = q_star_scm(i,j,k) + q(i,j,k)
            qcl_star(i,j,k)   = qcl_inc(i,j,k)    + qcl(i,j,k)
            qcf_star(i,j,k)   = qcf_inc(i,j,k)    + qcf(i,j,k)

            ! At present qcf2_inc etc are not used as there is
            ! no forcing option.

            qcf2_star(i,j,k)   = qcf2_star(i,j,k)   + qcf2(i,j,k)
            qrain_star(i,j,k)  = qrain_star(i,j,k)  + qrain(i,j,k)
            qgraup_star(i,j,k) = qgraup_star(i,j,k) + qgraup(i,j,k)

            cf_star (i,j,k) = cf_star (i,j,k) + cf(i,j,k)
            cfl_star(i,j,k) = cfl_star(i,j,k) + cfl(i,j,k)
            cff_star(i,j,k) = cff_star(i,j,k) + cff(i,j,k)
          END DO
        END DO
      END DO


      IF (l_qpos_for) THEN
        ! Check that forcing has not caused q to < qlimit

        ! NOTE: This has no equivalent in the full UM. It is justified
        !       as the forcing is decoupled from the model and so in the
        !       SCM there is no other mechanism to prevent a q < qlimit

        WRITE(sdum1,'(ES8.2)') qlimit
        DO k=1, nwet_lv
          DO j=1, rows
            DO i=1, row_length

              IF (q_star_scm(i,j,k) < qlimit) THEN
                ! The increment to q from the forcing routine will cause
                ! q to below qlimit. Don't allow this and output the
                ! required q to stop it going below qlimit
                dq_qpos_for(i,j,k) = qlimit - q_star_scm(i,j,k)
                q_star_scm(i,j,k)  = qlimit

                WRITE(sdum0,*) k

                CALL scm_message                                        &
                   ( 'q < '//TRIM(ADJUSTL(sdum1))//                     &
                     ': q reset on level '//TRIM(ADJUSTL(sdum0)) )
              ELSE
                dq_qpos_for(i,j,k) = 0.0
              END IF
            END DO      ! i
          END DO      ! j
        END DO      ! k

        IF (L_SCMDiags(SCMDiag_forc) .OR.                               &
            L_SCMDiags(SCMDiag_incs)) THEN

          ! DEPENDS ON: scmoutput
          CALL scmoutput                                                &
             ( dq_qpos_for, 'dq_qpos_for'                               &
             , 'Specific humidity adjustment to LS forcing to '//       &
               'maintain qlimit','kg/kg'                                &
             , t_avg, d_wet, default_streams, '', routinename )

        END IF ! l_scmdiags

      END IF ! l_qpos_for


      ! Convert t_inc_scm back to theta_star for call to Atmos_Physics2.
      theta_star(:,:,:) = t_inc_scm(:,:,:) / exner_theta_levels(:,:,:)

      !-----------------------------------------------------------------------
      ! theta_star, q_star_scm, qcl_star and qcf_star now all store
      ! Increment + Whole Values
      !-----------------------------------------------------------------------



        ! Convert X_star to mixing ratios
! DEPENDS ON: q_to_mix
      CALL q_to_mix                                                           &
        ( row_length, rows, wet_levels                                        &
        , halo_i, halo_j, q_star, qcl_star, qcf_star                          &
        , qcf2_star, qrain_star, qgraup_star                                  &
        , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                               &
        , mix_v_star, mix_cl_star, mix_cf_star                                &
        , mix_cf2_star, mix_rain_star, mix_graup_star )

!-----------------------------------------------------------------------------
! VARIABLE VALUES:
!   qcl_star   = liq. water       (at timelevel n)
!              + liq. water incs  (atmos_physics1 & forcing)
!   t_inc_scm  = temp.            (at timelevel n)
!              + temp. incs       (atmos_physics1 & forcing)
!   q_star_scm = vapour           (at timelevel n)
!              + vapour incs      (atmos_physics1 & forcing)
!   theta_star = pot. temp.       (at timelevel n)
!              + pot. temp. incs  (atmos_physics1 & forcing)
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Call physics2
!-----------------------------------------------------------------------------

! DEPENDS ON: atmos_physics2
      CALL atmos_physics2                                                     &
! Parallel variables
        ( row_length, rows, n_proc, n_procx, n_procy,rows, row_length         &
        , numcycles, cycleno                                                  &

! Model dimensions.
        , row_length, rows, rows, land_points, model_levels                   &
        , nice, nice_use, wet_levels, bl_levels                               &
        , st_levels, sm_levels, cloud_levels, land_ice_points                 &
        , soil_points, n_cca_lev, ntiles, tr_levels                           &
        , first_constant_r_rho_level, dim_cs1, dim_cs2                        &

! Model switches
        , l_regular, l_mr_physics2, l_dry, l_lbc_old                          &
        , lcal360, ltimer                                                     &
        , l_cld_area                                                          &

! Model parameters
        , rhcrit(1:wet_levels), co2_mmr                                       &
        , tr_vars, tr_ukca                                                    &

! In: coordinate information
        , unscl_dry_rho                                                       &
        , delta_lambda, delta_phi, dlambda_p, dphi_p                          &
        , wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v                        &
        , lat_rot_np, long_rot_np, f3_at_u                                    &

! In: time stepping information.
        , yearno, day, ihour, imin, time_sec                                  &

! River routing
        , row_length, row_length, xpa, xua, xva, ypa, yua, yva                &
        , g_p_field, g_r_field, a_steps_since_riv, river_row_length           &
        , river_rows, global_river_row_length, global_river_rows              &
        , river_vel, river_mcoef, i_river_vn                                  &
        , trivdir, trivseq, twatstor, inlandout_atm                           &

! Lake evaporation:
        , acc_lake_evap                                                       &

! Grid-to-grid river routing
        , r_area, slope, flowobs1, r_inext, r_jnext, r_land                   &
        , substore, surfstore, flowin, bflowin,                               &

! Diagnostics info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
          stashwork3, stashwork5, stashwork8, stashwork9, stashwork19         &
        , stashwork26                                                         &

! SCM diagnostics
        , nSCMdpkgs, l_SCMdiags                                               &

! In: data fields.
        , theta, q, qcl, qcf, qrain, qgraup, qcf2                             &
        , rho, u, v, w, w_adv, p, p_star                                      &
        , exner_rho_levels, exner_theta_levels                                &
        , land_sea_mask, p_theta_levels                                       &

! Ancillary fields and fields needed to be kept from timestep to
! timestep
        , land_index, land_ice_index, soil_index, canopy_gb                   &
        , snodep, hcon, hcap, v_crit, v_wilt, v_sat , sthf                    &
        , sthu, sil_orog_land                                                 &
        , ho2r2_orog, sd_orog_land, di, ice_fract, u_0, v_0, u_0, v_0         &
        , cca, iccb, icct, cclwp, ccw, lcbase                                 &
        , t_deep_soil, tsi, ti_n, ice_k_n, tstar                              &
        , z0msea, ice_fract_n, ice_thick_n, satcon, sathh, b_exp, smcl        &
        , t1_sd, q1_sd, zh, ddmfx, area_cloud_fraction, cf                    &
        , cfl, cff, ls_rain, ls_snow                                          &
        , micro_tends, photosynth_act_rad, rad_hr                             &
        , soil_clay, soil_silt, soil_sand, dust_mrel1 ,dust_mrel2             &
        , dust_mrel3, dust_mrel4, dust_mrel5, dust_mrel6                      &
        , so2_high_level, so2_em, nh3_em, dms_em, soot_hilem                  &
        , soot_em, ocff_hilem, ocff_em, co2_emits, co2flux                    &
        , deep_flag, past_precip, past_conv_ht                                &

! InOut
        , theta_star, q_star_scm                                              &
        , qcl_star, qcf_star, qrain_star, qgraup_star, qcf2_star              &
        , cf_star, cfl_star, cff_star                                         &
        , u_inc_scm, v_inc_scm, w_inc_scm(:,:,1:model_levels)                 &
        , sum_eng_fluxes, sum_moist_flux                                      &

! InOut: tracer fields
        , aerosol, free_tracers, ukca_tracers                                 &
        , dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6    &
        , so2, dms, so4_aitken, so4_accu, so4_diss, nh3, soot_new             &
        , soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld               &
        , ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss, co2             &

! InOut: RIVERS
        , tot_surf_runoff, tot_sub_runoff                                     &

! Out: fields
        , rhokm, ch_term, ntml, cumulus, nbdsc, ntdsc                         &
        , rhcpt,row_length, rows                                              &

! Additional variables for Land Surface Scheme
        , frac_typ, frac_disturb, canht, lai                                  &
        , canopy(1:land_points,1:ntiles)                                      &
        , catch(1:land_points,1:ntiles), catch_snow, snow_grnd, snow_tile     &
        , z0_tile, z0h_tile, tstar_tile                                       &
        , infil_tile(1:land_points,1:ntiles)                                  &
        , rgrain(1:land_points,1:ntiles)                                      &
        , cs, gs, co2_dim_row, co2_dim_len                                    &
        , asteps_since_triffid, timestep_number                               &
        , g_leaf_acc                                                          &
        , g_leaf_phen_acc, npp_ft_acc, resp_w_ft, resp_s_acc                  &
        , land_pts_trif, npft_trif, dolr, lw_down, sw_tile, fland_ctile       &
        , tstar_land, tstar_sea, tstar_sice_n, tstar_sice                     &
        , albsoil, cos_zenith_angle                                           &

! INOUT variables for TKE based turbulence schemes
        , e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu                        &

! Additional variables required for large-scale hydrology:
        , fexp, gamtot, ti_mean, ti_sig, fsat, fwetl, zw                      &
        , sthzw, a_fsat, c_fsat, a_fwet, c_fwet                               &

! JULES 2 prognostics (InOut)
        , snowdepth, rho_snow_grnd                                            &
        , nsnow, ds, sice, sliq, tsnowlayer, rho_snow, rgrainl                &

! FLake lake scheme prognostics (InOut)
       , lake_depth, lake_fetch, lake_t_mean, lake_t_mxl                      &
       , lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape                      &
       , lake_g_dt                                                            &

! Cariolle ozone
        , ozone_tracer                                                        &

! Additional screen-level variables
        , tscrndcl_ssi, tscrndcl_tile, tstbtrans                              &

! Variables required for COSP (out)
        , cosp_crain_3d, cosp_csnow_3d                                        &

! Prescribed surface forcing and roughness lengths
        , l_flux_bc, flux_e_scm, flux_h_scm, l_spec_z0, z0m_scm, z0h_scm      &

! Additional variables for SCM diagnostics
        , nfor, ntrop, time_string, daycount, cf                              &
        , nconvars, ntrad1, conv_mode, .TRUE.                                 &

! Error information
        , error_code                                                          &
        )

!-----------------------------------------------------------------------------
! VARIABLE VALUES
!   qcl_star   = liq. water        (at timelevel n)
!              + liq. water incs   (atmos_physics1, forcing & atmos_physics2)
!   q_star_scm = vapour            (at timelevel n)
!              + vapour incs       (atmos_physics1, forcing & atmos_physics2)
!   theta_star = pot. temp.        (at timelevel n)
!              + pot. temp. incs   (atmos_physics1, forcing & atmos_physics2)
!-----------------------------------------------------------------------------


      IF (stats .OR. obs .OR. noforce) THEN

! rho and pressure out of sync with T, but can't update with IdealGL
! settings i.e. P=rho.R.T as this will cause a crash in the physics
! routines. I assume that Pressure and density need to be updated
! using exner_prime


        ! Update theta, q and winds
        DO k=1, model_levels
          DO j=1, rows
            DO i=1, row_length

              theta(i,j,k) = theta_star(i,j,k)

              !  Convert theta back to t
              t(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)

              ! add wind increments to winds
              u(i,j,k) = u(i,j,k) + u_inc_scm(i,j,k)
              v(i,j,k) = v(i,j,k) + v_inc_scm(i,j,k)

              ! Vertical wind either prescribed for vertical advection
              ! forcing and updated after s_forcing or left at initial
              ! value.

            END DO
          END DO
        END DO

        DO k=1, wet_levels
          DO j=1, rows
            DO i=1, row_length
              q      (i,j,k) = q_star_scm  (i,j,k)
              qcl    (i,j,k) = qcl_star    (i,j,k)
              qcf    (i,j,k) = qcf_star    (i,j,k)
              qcf2   (i,j,k) = qcf2_star   (i,j,k)
              qrain  (i,j,k) = qrain_star  (i,j,k)
              qgraup (i,j,k) = qgraup_star (i,j,k)

              IF (l_pc2) THEN
                cf (i,j,k) = cf_star (i,j,k)
                cfl(i,j,k) = cfl_star(i,j,k)
                cff(i,j,k) = cff_star(i,j,k)
              END IF  ! l_pc2

            END DO
          END DO
        END DO

      END IF   ! stats .OR. obs .OR. noforce

!-----------------------------------------------------------------------------
! VARIABLE VALUES FOR THE NON-PC2 SITUATION.
!   (FOR PC2 SEE NEXT SECTION OF CODE FOR ADDITIONAL TERMS)
!
!   qcl   = liq. water  (at timelevel n+1)
!   q     = vapour      (at timelevel n+1)
!   theta = pot. temp.  (at timelevel n+1)
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! PC2 section
!-----------------------------------------------------------------------------

      IF (l_pc2) THEN

! The PC2 section needs to
! 1) increment condensation and cloud fractions due to the forcing
! 2) initiate cloud.
! 3) set area_cloud_fraction to bulk_cloud_fraction, cf, and correct
!    theta.
!
! 1. Calculate condensation increments which result from the forcing
! by using the homogeneous forcing approach

!-----------------------------------------------------------------------------
! SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------------
      IF (main_diag_switch /= 0) THEN
        IF (l_scmdiags(scmdiag_pc2)) THEN

!         ! save fields before incrementing with PC2 response to forcing
          q_earliest   (:,:,:) = q  (:,:,:)
          qcl_earliest (:,:,:) = qcl(:,:,:)
          t_earliest   (:,:,:) = t  (:,:,:)

        END IF ! l_scmdiags(scmdiag_pc2)
      END IF ! main_diag_switch /= 0

! DEPENDS ON: pc2_homog_plus_turb
        CALL pc2_homog_plus_turb                                              &
          ( p_theta_levels(1,1,1:wet_levels), wet_levels                      &
          , timestep, t, cf, cfl, cff                                         &
          , q, qcl, t_forcing, q_forcing, qcl_forcing, p_forcing              &
          , 0.0, 0.0, l_mr_pc2)

! 1a. We have already applied the forcing to q, qcl and t in the
!     forcing section. PC2_homog_plus_turb has now done this again
!     so we need to subtract off these increments. Just the
!     condensation will therefore remain.

        DO k=1, wet_levels
          DO j=1, rows
            DO i=1, row_length
              q(i,j,k)   = q(i,j,k)   - q_forcing(i,j,k)
              qcl(i,j,k) = qcl(i,j,k) - qcl_forcing(i,j,k)
              t(i,j,k)   = t(i,j,k)   - t_forcing(i,j,k)
            END DO ! i
          END DO ! j
        END DO ! k

      IF (main_diag_switch /= 0) THEN
        IF (l_scmdiags(scmdiag_pc2)) THEN
          qcl_inc (:,:,:)   = qcl(:,:,:) - qcl_earliest (:,:,:)
          q_star_scm(:,:,:) = q(:,:,:)   - q_earliest   (:,:,:)
          t_inc_scm(:,:,:)  = t(:,:,:)   - t_earliest   (:,:,:)
          cf_work (:,:,:)   = cf (:,:,:) - cf_star (:,:,:)
          cfl_work(:,:,:)   = cfl(:,:,:) - cfl_star(:,:,:)
          cff_work(:,:,:)   = cff(:,:,:) - cff_star(:,:,:)

! DEPENDS ON: scmoutput
          CALL scmoutput(qcl_inc, 'dqcl_pc2forc'                              &
            , 'PC2 qcl increment response to forcing', 'kg/kg'                &
            , t_inst, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(q_star_scm, 'dq_pc2forc'                             &
            , 'PC2 q increment response to forcing', 'kg/kg'                  &
            , t_inst, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(t_inc_scm, 'dt_pc2forc'                              &
            , 'PC2 T increment response to forcing', 'K'                      &
            , t_inst, d_all, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(cf_work, 'dbcf_pc2forc'                              &
            , 'PC2 bulk cf increment response to forcing', ''                 &
            , t_inst, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(cfl_work, 'dcfl_pc2forc'                             &
            , 'PC2 cfl increment response to forcing', ''                     &
            , t_inst, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(cff_work, 'dcff_pc2forc'                             &
            , 'PC2 cff increment response to forcing', ''                     &
            , t_inst, d_wet, default_streams, '', routinename)

        END IF ! l_scmdiags(scmdiag_pc2)
      END IF ! main_diag_switch /= 0

! 2. Initiate cloud if necessary and check for sensible values.
!    Use a dummy argument to receive the increment information.
!    NOTE: SCM increments are output within pc2_initiation_ctl

! DEPENDS ON: pc2_initiation_ctl
        CALL pc2_initiation_ctl                                               &
          ( row_length, rows                                                  &
          , .FALSE.                                                           &
          , l_mr_pc2, L_ACF_Cusack, L_cld_area                                &
          , timestep                                                          &
          , nscmdpkgs, l_scmdiags                                             &
          , t, q, qcl, qcf, cf, cfl, cff                                      &
          , rhts, tlts, qtts, ptts                                            &
          , area_cloud_fraction, p, p_star, p_theta_levels(1,1,1)             &
          , iccb, cumulus                                                     &
          , rhcpt                                                             &
          , t_work, q_work, qcl_work, qcf_work, cf_work, cfl_work, cff_work)


        IF (.NOT. l_cld_area) THEN

! 3. Set area_cloud_fraction to bulk_cloud_fraction, cf, and update
!    theta from new temperature

          DO k=1, wet_levels
            DO j=1, rows
              DO i=1, row_length
                area_cloud_fraction(i,j,k) = cf(i,j,k)
                theta(i,j,k) = t(i,j,k)/exner_theta_levels(i,j,k)
              END DO ! i
            END DO ! j
          END DO ! k

        ELSE IF (l_cld_area) THEN

          IF (l_acf_brooks) THEN

! DEPENDS ON: ls_acf_brooks
            CALL ls_acf_brooks(                                               &
                delta_lambda, delta_phi                                       &
              , fv_cos_theta_latitude, cf, cfl, cff                           &
              , cumulus, area_cloud_fraction )

          END IF ! l_ACF_brooks

        END IF ! l_cld_area


      END IF  ! l_pc2

!-----------------------------------------------------------------------------
! End of PC2 section.
!
! Variables q, qcl, t, theta and cloud fractions
! are now all at timelevel n+1
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------
! Check for negative q, qcl and qcf
! This would be done by q_pos in the full UM, l_qpos is specified
! in UMUI dyn_run namelist
!-----------------------------------------------------------------------
      IF (l_qpos) THEN
        ! Check q, qcl, and qcf. are not -ve
        dq_qpos   (:,:,:) = 0.0
        dqcl_qpos (:,:,:) = 0.0
        dqcf_qpos (:,:,:) = 0.0

        DO k=1, wet_levels
          DO j=1, rows
            DO i=1, row_length

              IF (q(i,j,k) < 0.0) THEN
                WRITE(sdum0,*) k
                CALL scm_message                                       &
                   ( 'q < qlimit: q added on level '//                 &
                     TRIM(ADJUSTL(sdum0)))
                dq_qpos(i,j,k) = qlimit - q(i,j,k)
                q(i,j,k) = qlimit
              END IF

              IF (qcl(i,j,k) < 0.0) THEN
                WRITE(sdum0,*) k
                CALL scm_message                                       &
                   ( 'qcl < 0.0: qcl added on level '//                &
                     TRIM(ADJUSTL(sdum0)))
                dqcl_qpos(i,j,k) = - qcl(i,j,k)
                qcl(i,j,k) = 0.0
              END IF

              IF (qcf(i,j,k) < 0.0) THEN
                WRITE(sdum0,*) k
                CALL scm_message                                       &
                   ( 'qcf < 0.0: qcl added on level '//                &
                     TRIM(ADJUSTL(sdum0)))
                dqcf_qpos(i,j,k) = - qcf(i,j,k)
                qcf(i,j,k) = 0.0
              END IF

            END DO      ! i
          END DO      ! j
        END DO      ! k

        IF (main_diag_switch /= 0) THEN
          IF (l_scmdiags(scmdiag_forc) .OR.                            &
              l_scmdiags(scmdiag_incs)) THEN

            ! DEPENDS ON: scmoutput
            CALL scmoutput(dq_qpos,'dq_qpos'                           &
               , 'Q inc to prevent q < qlimit','kg/kg'                 &
               ,  t_avg,d_wet,default_streams,'',routinename)

            ! DEPENDS ON: scmoutput
            CALL scmoutput(dqcl_qpos,'dqcl_qpos'                       &
               , 'Qcl inc to prevent -ve qcl','kg/kg'                  &
               ,  t_avg,d_wet,default_streams,'',routinename)

            ! DEPENDS ON: scmoutput
            CALL scmoutput(dqcf_qpos,'dqcf_qpos'                       &
               , 'Qcf inc to prevent -ve qcf','kg/kg'                  &
               ,  t_avg,d_wet,default_streams,'',routinename)

          END IF        ! l_scmdiags
        END IF        ! main_diag_switch /= 0

      END IF        ! l_qpos



!-----------------------------------------------------------------------------
! Calculate some final diagnostics, and write the lot out if needs be.
!-----------------------------------------------------------------------------

      IF (main_diag_switch /= 0) THEN
!
!-----------------------------------------------------------------------------
! SCM PC2 Diagnostics Package
!-----------------------------------------------------------------------------
        IF (l_scmdiags(scmdiag_pc2)) THEN

! DEPENDS ON: scmoutput
          CALL scmoutput(qcl, 'qcl_n1_afterpc2'                               &
            , 'qcl at end of timestep after pc2', 'kg/kg'                     &
            , t_inst, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(q, 'q_n1_afterpc2'                                   &
            , 'q at end of timestep after pc2', 'kg/kg'                       &
            , t_inst, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(theta, 'th_n1_afterpc2'                              &
            , 'theta at end of timestep after pc2', 'K'                       &
            , t_inst, d_all, default_streams, '', routinename)

        END IF ! l_scmdiags(scmdiag_pc2)

!-----------------------------------------------------------------------------
! SCM Increments Diagnostics Package
!-----------------------------------------------------------------------------
        IF (l_scmdiags(scmdiag_incs)) THEN

          t_totalinc(:,:,:) = t(:,:,:) - t_start(:,:,:)
          u_totalinc(:,:,:) = u(:,:,:) - u_start(:,:,:)
          v_totalinc(:,:,:) = v(:,:,:) - v_start(:,:,:)
          w_totalinc(:,:,:) = w(:,:,:) - w_start(:,:,:)

          q_totalinc   (:,:,:) = q  (:,:,:) - q_start  (:,:,:)
          qcl_totalinc (:,:,:) = qcl(:,:,:) - qcl_start(:,:,:)
          qcf_totalinc (:,:,:) = qcf(:,:,:) - qcf_start(:,:,:)
          cf_totalinc  (:,:,:) = cf (:,:,:) - cf_start (:,:,:)
          cfl_totalinc (:,:,:) = cfl(:,:,:) - cfl_start(:,:,:)
          cff_totalinc (:,:,:) = cff(:,:,:) - cff_start(:,:,:)

! DEPENDS ON: scmoutput
          CALL scmoutput(t_totalinc, 'dt_total'                               &
            , 'Total increment to T', 'K'                                     &
            , t_avg, d_all, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(u_totalinc, 'du_total'                               &
            , 'Total increment to u', 'm/s'                                   &
            , t_avg, d_all, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(v_totalinc, 'dv_total'                               &
            , 'Total increment to v', 'm/s'                                   &
            , t_avg, d_all, default_streams, '', routinename)

          DO k=1, model_levels
            a2out(:,:,k) = w_totalinc(:,:,k)
          END DO

! DEPENDS ON: scmoutput
          CALL scmoutput(a2out, 'dw_total'                                    &
            , 'Total increment to w', 'm/s'                                   &
            , t_avg, d_all, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(q_totalinc, 'dq_total'                               &
            , 'Total increment to q', 'kg/kg'                                 &
            , t_avg, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(qcl_totalinc, 'dqcl_total'                           &
            , 'Total increment to qcl', 'kg/kg'                               &
            , t_avg, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(qcf_totalinc, 'dqcf_total'                           &
            , 'Total increment to qcf', 'kg/kg'                               &
            , t_avg, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(cf_totalinc, 'dbcf_total'                            &
            , 'Total increment to bulk cf', ' '                               &
            , t_avg, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(cfl_totalinc, 'dcfl_total'                           &
            , 'Total increment to cfl', ' '                                   &
            , t_avg, d_wet, default_streams, '', routinename)

! DEPENDS ON: scmoutput
          CALL scmoutput(cff_totalinc, 'dcff_total'                           &
            , 'Total increment to cff', ' '                                   &
            , t_avg, d_wet, default_streams, '', routinename)

          IF (geoforce) THEN
! DEPENDS ON: scmoutput
            CALL scmoutput(uinc_geo, 'du_geo'                                 &
              , 'Geostrophic forcing increment to u', 'm/s'                   &
              , t_avg, d_all, default_streams, '', routinename)

! DEPENDS ON: scmoutput
            CALL scmoutput(vinc_geo, 'dv_geo'                                 &
              , 'Geostrophic forcing increment to v', 'm/s'                   &
              , t_avg, d_all, default_streams, '', routinename)

            IF (ug_opt >=1) THEN
              geo_diag(:,:,:) = ug_scm(:,:,:)

              ! DEPENDS ON: scmoutput
              CALL scmoutput(geo_diag, 'u_g'                                  &
                 , 'Zonal geostrophic wind', 'm/s'                            &
                 , t_avg, d_all, default_streams, '', routinename)
            END IF

            IF (vg_opt >= 1) THEN
              geo_diag(:,:,:) = vg_scm(:,:,:)

              ! DEPENDS ON: scmoutput
              CALL scmoutput(geo_diag, 'v_g'                                  &
                 , 'Meridional geostrophic wind', 'm/s'                       &
                 , t_avg, d_all, default_streams, '', routinename)
            END IF

          END IF

        END IF ! l_SCMDiags(SCMDiag_incs)


         ! Store some diagnostic
! DEPENDS ON: dgnstcs_scm_main
        CALL dgnstcs_scm_main                                                 &
          ( row_length, rows, land_points, model_levels                       &
          , wet_levels, tr_levels, tr_vars, sm_levels, st_levels, ntype       &
          , rho, timestep, u, v, t                                            &
          , theta, q, qcl, qcf, cf, cca, ccw                                  &
          , t_deep_soil, p_star, tstar, smc, canopy_gb, snodep, zh            &
          , z0msea, smcl, sthu, sthf, gs, lw_incs                             &
          , photosynth_act_rad, tstar_tile, aerosol, free_tracers             &
          , p_theta_levels, p, iccb, icct, w, w_adv, area_cloud_fraction      &
          , cf, cfl, cff, cclwp, nscmdpkgs, l_scmdiags)

        ! Initialise the diagnostic output files if this is the
        ! end of the first timestep for which the system was on
        IF (scmop%first_pass) THEN
          scmop%first_pass = .FALSE.

          ! The list of diagnostics should be finalised now, so
          ! dump_streams_init can be called.
! DEPENDS ON: dump_streams_init
          CALL dump_streams_init                                              &
            ! (InOut)
            ( scmop                                                           &
            ! (In)
            , row_length, rows, model_levels, wet_levels                      &
            , bl_levels, cloud_levels                                         &
            , ozone_levels, st_levels, sm_levels, ntiles                      &
            , year_init, month_init, day_init                                 &
            , hour_init, min_init, sec_init                                   &
            , timestep, ndayin,nminin, nsecin, sec_day                        &
            , ndayin*full_daysteps+nstepsin                                   &
            , a_sw_radstep_prog, a_sw_radstep_diag                            &
            , z_top_of_model, first_constant_r_rho_level                      &
            , eta_theta, eta_rho, orog, r_theta_levels                        &
            , r_rho_levels, netcdf_chunksize )

          !=======================================================
          ! ADD REMINDER TO CHECK NFOR IS CORRECT
          ! Added here it doesn't get push off screen
          ! and missed by user
          !=======================================================

          WRITE(6,*) "========================================"
          WRITE(6,*) " Namelist NFOR = " // TRIM(ADJUSTL(sdum0))
          WRITE(6,*) " NOTE: INCORRECT SPECIFICATION OF NFOR  "
          WRITE(6,*) "       WILL PRODUCE UNINTENDED RESULTS. "
          WRITE(6,*) "========================================"
          WRITE(6,*) " "


        END IF

        ! Write output diagnostics to file(s)
! DEPENDS ON: dump_streams
        CALL dump_streams                                                     &
          ( scmop, day, time_sec, row_length, rows, model_levels, dayno_init  &
          , INT(timestep_number*timestep), site, lat, long                    &
          , time_initi, year, lcal360 )

      END IF ! main_diag_switch /= 0

      IF (l_ts_log) THEN
        CALL scm_message('Complete')
      END IF

    END DO   ! timestep_number

  END DO   ! daycount

! Close the output files
! DEPENDS ON: dump_streams_end
  CALL dump_streams_end(scmop)


  DEALLOCATE( arcl )                  ! Aerosol climatology array for NWP
  DEALLOCATE( true_latitude )
  DEALLOCATE( true_longitude )
  DEALLOCATE( cos_theta_latitude )
  DEALLOCATE( sec_theta_latitude )
  DEALLOCATE( sin_theta_latitude )
  DEALLOCATE( cos_theta_longitude )
  DEALLOCATE( sin_theta_longitude )
  DEALLOCATE( FV_cos_theta_latitude )
  DEALLOCATE( r_rho_levels )
  DEALLOCATE( r_theta_levels )
  DEALLOCATE( eta_rho_levels )
  DEALLOCATE( eta_theta_levels )


  CALL dealloc_common_scm
  CALL dealloc_forcing

  WRITE(6,*) '------------------------'
  WRITE(6,*) 'Run completed'
  WRITE(6,*) '------------------------'

9999 CONTINUE

  IF (lhook) CALL dr_hook('SCM_MAIN',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE scm_main
