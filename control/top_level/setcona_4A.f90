! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine SETCONA (ENDGAME VERSION) ------------------------------
!  
!   Purpose : to calculate additional constants derived from
!        information in the input dump and pass out into the model via
!        modules and the argument lists in argcona and arglndm.
!  
!   Method:
!      0. Initialise vertical coordinates from eta values and surface
!         orography.
!      0.1 Initialise data arrays held in secondary dump space.
!      1. Trigonometric functions:
!          Convert from degrees to radians.
!          Calculate trig functions for this grid.
!          Update halos for trig functions.
!      2. Set up semi-lagrangian advection options.
!         Initialise land/soil points in index arrays.
!      3. Determine model level boundaries for L/M/H cloud diagnostics.
!      4. Add locally derived parameters.
!  
!   Language: FORTRAN 77 + common extensions also in Fortran 90.
!   Programming standard; Unified Model Documentation Paper No. 3
!  
!   Documentation : Unified Model Documentation Paper No P0
!  
!  -------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

      SUBROUTINE Setcona_4A                                             &
! Input data (including some logical fields)
     &     (eta_theta_in,eta_rho_in,Smvcst,Land_sea_mask,Orog,          &
     &      grad_x, grad_y, orog_unfilt,                                &
     &      rho,exner_rho_levels,                                       &
     &      orog_lbc,exner_lbc,                                         &
! Input size and control variables
     &      LENRIMA,LBC_SIZEA,LBC_STARTA,                               &
     &      RIMWIDTHA, RIMWEIGHTSA,                                     &
     &      global_row_length, global_rows,                             &
     &      MODEL_LEVELS,ROWS,N_ROWS,ROW_LENGTH,                        &
     &      LAND_FIELD,WET_LEVELS,boundary_layer_levels,                &
     &      first_constant_r_rho_level,cloud_levels,z_top_of_model,     &
     &      ozone_levels,bl_levels,tr_levels, tr_vars, tr_ukca,         &
     &      height_gen_method,                                          &
! other constants
     &      Delta_lambda_in,Delta_phi_in,Base_phi_in,Base_lambda_in,    &
     &      lat_rot_NP_in,long_rot_NP_in,                               &
! VarRes grid info in degrees
     &      Lambda_p_in, Lambda_u_in, Phi_p_in, Phi_v_in,               &
     &      RowDepCStart,                                               &
! Initialise variables in secondary D1:
     &      exner_theta_levels,p_theta_levels,p,p_star,                 &
! Output data
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
     &     A_LEN2_COLDEPC, A_LEN2_ROWDEPC,                              &
           icode,Cmessage)

USE conversions_mod, ONLY: pi_over_180, pi
USE rimtypes
USE solinc_data, ONLY: L_orog, slope_angle, horiz_limit, l_skyview
USE level_heights_Mod
USE trignometric_Mod
USE dyn_coriolis_Mod
USE dyn_var_res_Mod
USE diff_coeff_Mod
USE turb_diff_mod
USE rad_mask_trop_Mod
USE rad_input_mod, ONLY:                                                &
    L_use_orog_corr, L_use_grad_corr,                                   &
    l_use_skyview,L_rad_deg,l_orog_unfilt, l_use_cariolle
USE switches, ONLY: l_ctile    
USE rot_coeff_Mod
USE swapable_field_mod, ONLY : swapable_field_pointer_type
USE atm_fields_bounds_Mod
USE earth_constants_mod, ONLY: g, earth_radius, omega, two_omega
USE atmos_constants_mod, ONLY: recip_kappa, p_zero
USE dust_parameters_mod, ONLY: l_twobin_dust 

USE global_2d_sums_mod, ONLY: global_sum_method, global_sum_reprod_dd,  &
                              global_sum_reprod_orig

USE bl_option_mod, ONLY: l_us_blsol, puns,pstb

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE Control_Max_Sizes
USE UM_ParVars
USE ereport_mod, ONLY : ereport
USE PrintStatus_mod

USE integrity_mod
USE set_horiz_grid_mod
USE set_var_horiz_grid_mod
USE horiz_grid_mod
USE ref_pro_mod ! needed due to the inclusion of cruntimc
USE eg_alpha_mod ! needed due to the inclusion of cruntimc
USE eg_alpha_ramp_mod ! needed due to the inclusion of cruntimc
USE comobs_mod, ONLY: nobtypmx
USE eg_parameters_mod  ! needed due to the inclusion of cruntimc
USE alloc_grid_mod
USE calc_global_grid_spacing_mod
 
USE domain_params
USE eg_swap_bounds_mod

USE lam_lbc_weights_mod, ONLY: init_lbc_wghts

USE problem_mod, ONLY: standard, monsoon, dynamical_core,               &
                       idealised_problem, standard_namelist
USE trophgt1_mod, ONLY: z_min_trop, z_max_trop
USE mcica_mod, ONLY: gridbox_size
USE sl_input_mod, ONLY:                                                 &
          thmono_levels,halo_lam,halo_phi,look_lam,look_phi,            &
          L_conserv,L_mono,L_high,high_order_scheme,thmono_height,      &
          recip_dphi,recip_dlam



USE dynamics_input_mod, ONLY:   L_mix_ratio,                            &
                                L_LBC_balance,L_lbc_new,L_transparent,  &
                                L_regular,n_rims_to_do,                 &
                                GCR_diagnostics,GCR_its_avg_step,       &
                                L_lbc_old

USE dynamics_testing_mod, ONLY:                                         &
                L_idealised_data,                                       &
                lambda_p_end,phi_p_end,lambda_u_end,phi_v_end,          &
                dlambda_p_end,dphi_p_end,dlambda_u_end,dphi_v_end,      &      
                lam_var, phi_var, var_ratio,lam_ratio,phi_ratio,        &
                phi_frac,lam_frac

USE gcr_input_mod, ONLY:                                                &
          GCR_its_switch, GCR_max_its,GCR_min_its,GCR_max_time,         &
          GCR_min_time,GCR_sum_its

USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain

USE cloud_inputs_mod, ONLY: l_pc2

USE murk_inputs_mod, ONLY: l_murk, l_murk_advect

USE dust_parameters_mod, ONLY: l_dust

USE um_input_control_mod,  ONLY:                                        &
     model_domain,         problem_number,                              &
     l_int_uvw_lbc,                                                     &
                           l_sulpc_so2,       l_sulpc_dms,              &
     l_sulpc_nh3,          l_soot,            l_biomass,                &
     l_co2_interactive,    l_ocff,            l_nitrate,                &
     super_array_size,     moisture_array_size

USE missing_data_mod, ONLY: rmdi  
USE highos_mod,       ONLY: quinticLagrange

IMPLICIT NONE

!     INCLUDED COMDECKS
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

      INTEGER                                                           &
     &       LENRIMA(Nfld_max,NHalo_max,Nrima_max),                     &
                                                     ! IN LBC data len
     &       LBC_SIZEA(4,Nfld_max,NHalo_max,Nrima_max),                 &
                                                       ! IN LBC size
     &       LBC_STARTA(4,Nfld_max,NHalo_max,Nrima_max),                &
                                                         ! IN LBC start
     &       RIMWIDTHA(Nrima_max),                                      &
                                    ! IN RIM width
     &       global_row_length,                                         &
                                 ! IN total number of point in a row
     &       global_rows,                                               &
                           ! IN total number of rows in model
     &       MODEL_LEVELS,ROWS,N_ROWS,ROW_LENGTH,ICODE,                 &
     &       LAND_FIELD,                                                &
                           ! IN Number of land points in model from umui
     &       WET_LEVELS,                                                &
                           ! IN Number of wet levels in model, from umui
     &       first_constant_r_rho_level,                                &
                                         !(IN) 1st constant height level
     &       boundary_layer_levels,                                     &
                                    ! (IN) Num. of boundary layer levels
     &       height_gen_method,                                         &
                                  ! (IN) Method for generating heights
     &       CLOUD_LEVELS                                               &
                           ! IN No of cloudy levels in the model
     &      ,tr_levels, tr_vars, tr_ukca

      INTEGER RowDepCStart
                        ! IN Start of Row dependent constant
      CHARACTER(LEN=80)                                                    &
     &       CMESSAGE              ! Error message if ICODE >0
      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='SETCONA_4A')

! The following are used in the decalaration of Lambda_u_in, Phi_v_in
!  to avoid out-of-bounds errors
      INTEGER :: A_LEN2_COLDEPC, A_LEN2_ROWDEPC

      REAL                                                              &
            ! Input arguments (grid constants)
     &       Delta_lambda_in                                            &
                             ! EW (x) grid spacing in degrees
     &      ,Delta_phi_in                                               &
                             ! NS (y) grid spacing in degrees
     &      ,Base_phi_in                                                &
                             ! Latitude of first theta point in degrees
     &      ,Base_lambda_in                                             &
                             ! Longitude of first theta point in degs
     &      ,lat_rot_NP_in                                              & 
                           ! Real latitude of 'pseudo' N pole in degs
     &      ,long_rot_NP_in                                             & 
                           ! Real longitude of 'pseudo' N pole in degs
     &      ,z_top_of_model                                             &
                            ! (IN) Height of top of model in metres
     &,      RIMWEIGHTSA(RIMWIDTHA(rima_type_norm))  ! IN RIM weights

! In the regular grid case A_LEN2_COLDEPC, A_LEN2_ROWDEPC are zero, but
! Lambda_u_in, Phi_v_in must have size 1 to avoid out-of-bounds errors
      REAL                                                              &
            ! Input VarRes grid info in degrees
     &   Lambda_p_in(global_row_length)                                 &
                                               !IN EW and NS VarRes grid
     &  ,Lambda_u_in(MIN(global_row_length, &
                         global_row_length *A_LEN2_COLDEPC +1))         &
                                               !IN EW and NS u,v grid
     &  ,Phi_p_in(global_rows)                                          &
                                               !location in degrees
     &  ,Phi_v_in(MIN(global_rows+1,                                    &
                      global_rows *A_LEN2_ROWDEPC +1)) !location in degrees

      REAL                                                              &
            ! Array arguments with intent(IN):
     &     Smvcst(land_field)                                           &
                                      ! IN Volumetric saturation point
     &    ,eta_theta_in(0:model_levels)                                 &
                                        ! IN Eta values for theta levs
     &    ,eta_rho_in(model_levels)                                     &
                                        ! IN Eta values for rho levels
     &    ,Orog(row_length, rows)                                       &
                                        ! IN Orography (on all points)
     &    ,orog_unfilt(row_length, rows)                                &
                                        ! IN Unfiltered orography
     &    ,grad_x(row_length*rows)                                      &
                                        ! IN Orographic X-gradient
     &    ,grad_y(row_length*rows)                                      &
                                        ! IN Orographic Y-gradient
          ,rho(pdims_s%i_start:pdims_s%i_end,                           &
               pdims_s%j_start:pdims_s%j_end,                           &
               pdims_s%k_start:pdims_s%k_end)                           &
                           ! density*(r**2): used in call to Calc_P_star
          ,exner_rho_levels(pdims_s%i_start:pdims_s%i_end,              &
                            pdims_s%j_start:pdims_s%j_end,              &
                            pdims_s%k_start:pdims_s%k_end+1)            &
                                                     ! Exner at p levels
     &,      orog_lbc(LENRIMA(fld_type_p,halo_type_extended,            &
     &                        rima_type_orog))                          &
                                                 ! Orography LBC
     &,      exner_lbc(LENRIMA(fld_type_p,halo_type_extended,           &
                                                              !Exner LBC
     &                         rima_type_norm),MODEL_LEVELS+1)

      REAL                                                              &
            ! Output args (secondary arrays calculated from dump)
        p_theta_levels(tdims_s%i_start:tdims_s%i_end,                   &
                       tdims_s%j_start:tdims_s%j_end,                   &
                       tdims_s%k_start:tdims_s%k_end)                   &
                                               ! press on theta levs Pa
      , p(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end,                                &
          pdims_s%k_start:pdims_s%k_end+1)                              &
                                                              ! in Pa
     &, p_star(row_length, rows)                                        &
                                 ! surface pressure (Pa)
      , exner_theta_levels(tdims_s%i_start:tdims_s%i_end,               &
                           tdims_s%j_start:tdims_s%j_end,               &
                           tdims_s%k_start:tdims_s%k_end)
                                                   ! Exner at theta levs

      LOGICAL                                                           &
     &       Land_sea_mask(row_length,rows)

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
! C_ETA_PMSL start
! ETA_PMSL is the ETA value used to determine the model level
! used in PMSL reduction calculations.
      REAL,PARAMETER:: ETA_PMSL=0.795
      INTEGER LEVNO_PMSL_CALC   !  Model level for PMSL reduction calc.
      COMMON /PMSLCALC/ LEVNO_PMSL_CALC
! C_ETA_PMSL end
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

! Local variables

      REAL :: orog_halos(1-halo_i:row_length+halo_i,                    &
                         1-halo_j:rows+halo_j)

      REAL :: orog_horiz(1-horiz_limit:row_length+horiz_limit,          &
                         1-horiz_limit:rows+horiz_limit)

      REAL :: orog_global(global_row_length,global_rows)

      REAL :: phi_horiz(1-horiz_limit:rows+horiz_limit)
      REAL :: dphi_horiz(1-horiz_limit:rows+horiz_limit)
      REAL :: dlambda_horiz(1-horiz_limit:row_length+horiz_limit)

      REAL                                                              &
     &   Delta_lambda_wk                                                &
                             ! EW (x) grid spacing in degrees
     &  ,Delta_phi_wk                                                   &
                             ! NS (y) grid spacing in degrees
     &  ,Base_phi_wk                                                    &
                             ! Latitude of first theta point in degrees
     &  ,Base_lambda_wk                                                 &
                             ! Longitude of first theta point in degrees
     &  ,latitude_rot(row_length, rows)                                 &
                                         ! rot. latit. in degrees (LAM)
     &  ,longitude_rot(row_length, rows)                                &
                                         ! rot. longit. in degrees (LAM)
     &  ,true_latitude_b(row_length, rows)                              &
                                           ! true latitude on B-grid
     &  ,true_longitude_b(row_length,rows)                              &
                                           ! true longitude on B-grid
     &  ,bear_rot_NP(row_length, rows)                                  & 
                                         ! Bearing of 'pseudo' N pole
     &  ,temp1                                                          &
                   ! Used in calculation of Coriolis terms (LAM only)
     &  ,temp2                                                          &
                   ! Used in calculation of Coriolis terms (LAM only)
     &   , f_plane_rad                                                  &
                        ! f_plane latitude in radians
     &   , ff_plane_rad                                                 &
                         ! f_plane latitude in radians
     &   , f1_temp                                                      &
                    ! f1 term for f-plane or ff-plane
     &   , f2_temp                                                      &
                    ! f2 term for f-plane or ff-plane
     &   , f3_temp                                                      &
                    ! f3 term for f-plane or ff-plane
     &  , scale                                                         &
     &  ,CLOUD_BOUND(NUM_CLOUD_TYPES+1)                                 &
                                        ! boundaries of cloud types
     &  ,r_ref_theta(model_levels)                                      &
                                   ! Local dynamic array
        ,r_ref_rho(model_levels)                                        &
                                   ! Local dynamic array
        ,xscale(row_length,rows)                                        &
!           Width of gridbox (x) in m.
        ,yscale
!           Width of gridbox (y) in m.

      LOGICAL                                                           &
     &   landmask(1-offx:row_length+offx, 1-offy:rows+offy)             &
     &  ,L_error                                                        &
     &  ,L_varres

      INTEGER                                                           &
              ! Mostly loop counters, but j also used for interp points
     &   I,KK,                                                          &
     &   J,GI,GJ,                                                       &
     &   j0,j1,k,                                                       &
     &   LEVEL                                                          &
                         ! Used to set up cloud type boundaries
     &  ,first_row                                                      &
     &  ,last_row                                                       &
     &  ,info                                                           &
     &  , active_levels

      INTEGER, PARAMETER :: height_gen_original = 1 ! methods for height
      INTEGER, PARAMETER :: height_gen_smooth   = 2 ! generation

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      INTEGER :: i_field

      TYPE(swapable_field_pointer_type) :: fields_to_swap(2)

      INTEGER :: ozone_levels,bl_levels
      REAL    :: z_offset


!----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SETCONA_4A',zhook_in,zhook_handle)

!  Set processor properties in module for diagnostic  printing
!  to avoid having to include in argument list from calling routine
!               calls for printing need to be inserted by user
! DEPENDS ON: init_proc_info
        CALL init_proc_info(                                            &
     &                      global_row_length, global_rows,             &
     &                      model_domain,                               &
     &                      mype, nproc, nproc_x, nproc_y,              &
     &                      gc_proc_col_group, gc_proc_row_group,       &
     &                      datastart, at_extremity) 


! Set flag for orography correction in module solinc_data
      L_orog = L_use_orog_corr .or. L_use_grad_corr
      l_skyview = l_orog .AND. l_use_skyview
! Grid type in dump: L_varres = T / F for variable / unif grid
      L_varres = .FALSE.
      If (RowDepCStart >= 1) L_varres = .TRUE.
! ----------------------------------------------------------------------
! 0. Set-up vertical co-ordinate arrays
! ----------------------------------------------------------------------

      z_offset = Earth_Radius

      DO j=1,ROWS
        DO i=1,ROW_LENGTH
          orog_halos(i,j)=orog(i,j)
        END DO
      END DO

! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS(orog_halos,row_length,rows,1,                    &
                       halo_i,halo_j,fld_type_p,.FALSE.)

! ----------------------------------------------------------------------
! 0.1 Write out some useful information at start of run
! ----------------------------------------------------------------------

      WRITE(6,'(A)') '  '
      IF (global_sum_method == global_sum_reprod_orig .OR.              &
          global_sum_method == global_sum_reprod_dd) THEN
        WRITE(6,'(A)') '**********************  This run uses bit-'//   &
                       'reproducible code ********************'
      ELSE
        WRITE(6,'(A)') '***************** This run uses fast code, so'//&
                       ' different processor **************'
        WRITE(6,'(A)') '*****************     configurations may not '//&
                       'bit-reproduce        **************'
      END IF

      WRITE(6,'(A)') '  '
      WRITE(6,'(A)') '************************  PE configuration for '// &
                     'this run  ***********************'
      IF ( nproc_x > 1 .and. nproc_y > 1 ) THEN
        WRITE(6,'(A,I5,A,I5,A)') '****  PE configuration: ', nproc_x,    &
            ' processors East-West ', nproc_y, ' processors North-South *'
      ELSE IF ( nproc_x > 1 .and. nproc_y == 1 ) THEN
        WRITE(6,'(A,I5,A,I5,A)') '****  PE configuration: ', nproc_x,    &
            ' processors East-West ', nproc_y, ' processor North-South *'
      ELSE IF ( nproc_x == 1 .and. nproc_y > 1 ) THEN
        WRITE(6,'(A,I5,A,I5,A)') '****  PE configuration:', nproc_x,     &
            ' processor East-West ', nproc_y, ' processors North-South *'
      ELSE ! nproc_x = 1 .and. nproc_y =1 ) then
        WRITE(6,'(A)') '****  PE configuration:  Single processor 1x1 * '
      END IF ! nproc_x > 1 .and. nproc_y > 1 )

      if ( offx >= halo_i .or. offy >= halo_j ) then
          WRITE(6,'(A)') '  '
          WRITE(6,'(A)') '  ****  DANGER  ****   '
          WRITE(6,'(A)') '  *** SMALL halo >= LARGE halo  *****   '
          WRITE(6,'(A)') '  This could result in overwriting  '
          ICODE    = 20
          CMESSAGE ='SETCONA: Large halo size is too small'
          CALL Ereport(RoutineName,ICODE,Cmessage)
      end if ! offx >= halo_i .or. offy >= halo_j

      WRITE(6,'(A)') '  '
      IF (PrintStatus >= PrStatus_Oper) THEN
        IF (L_mix_ratio) THEN
          WRITE(6,'(A)') '*****  This run uses mixing ratios for the'//  &
                         ' moist variables in the dynamics *****'
        ELSE
          WRITE(6,'(A)') '**  This run uses specific humidities for '//  &
                         'the moist variables in the dynamics **'
        END IF !  L_mix_ratio
      END IF ! PrintStatus >= PrStatus_Oper

      If (model_domain  ==  mt_LAM) Then
      
        L_lbc_old = .true.
        If ( L_lbc_new ) L_lbc_old = .false.

        write(6,*) '  '
        Write ( Unit = 6, fmt=*) ' ***   LAM set-up ***' 
        
        Write ( Unit = 6, fmt=*) 'LBC frame size for solver, '          &
     &,                          ' n_rims_to_do = ', n_rims_to_do
        If( L_LBC_balance ) then
          Write ( Unit = 6, fmt=*) 'L_LBC_balance = ', L_LBC_balance    &
     &,          ' Impose vertically balanced Exner pressures'          &
     &,          ' and rho and set w=0 in lbcs '
        Else !  L_LBC_balance = .false.
          Write ( Unit = 6, fmt=*) 'L_LBC_balance = ', L_LBC_balance    &
     &,                            ' No balancing of lbcs '
        EndIf ! L_LBC_balance      
        If( L_lbc_new ) then
          Write ( Unit = 6, fmt=*) 'L_lbc_new = ', L_lbc_new            &
     &,                              'Use new lbc algorithm '
          If ( L_transparent ) Then
            Write ( Unit = 6, fmt=*) 'L_transparent = ', L_transparent  &
     &,         'Does not apply lbcs to horizontal winds in the solver'
          Else !  L_transparent = .false.
            Write ( Unit = 6, fmt=*) 'L_transparent = ', L_transparent  &
     &,         'Applies lbcs to horizontal winds in the solver'
          End If ! L_transparent
        Else !  L_lbc_new = .false.
          Write ( Unit = 6, fmt=*) 'L_lbc_new = ', L_lbc_new            &
     &,                              'Use old lbc algorithm '
        EndIf ! L_lbc_new
        If( L_int_uvw_lbc) then
          Write ( Unit = 6, fmt=*) 'L_int_uvw_lbc =', L_int_uvw_lbc     &
     &,          ' Advecting winds interpolated in lateral boundaries '
        Else !  L_int_uvw_lbc = .false.
          Write ( Unit = 6, fmt=*) 'L_int_uvw_lbc =', L_int_uvw_lbc     &
     &,                ' Extrapolated advecting winds from lbc file '
        EndIf ! L_int_uvw_lbc
        write(6,*) '  '

! DEPENDS ON: set_lateral_boundaries
        CALL set_lateral_boundaries(                                    &
     &    ROW_LENGTH,ROWS,halo_i,halo_j,1,fld_type_p,orog_halos,        &
     &    LENRIMA(fld_type_p,halo_type_extended,rima_type_orog),        &
     &    LBC_SIZEA(1,fld_type_p,halo_type_extended,rima_type_orog),    &
     &    LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_orog),   &
     &    halo_i,halo_j,orog_lbc,RIMWIDTHA(rima_type_orog),             &
     &    RIMWIDTHA(rima_type_orog),RIMWEIGHTSA,                        &
     &    at_extremity,                                                 &
     &    .FALSE.,.TRUE.)

! DEPENDS ON: set_lateral_boundaries
        CALL set_lateral_boundaries(                                    &
     &    ROW_LENGTH,ROWS,Offx,Offy,MODEL_LEVELS+1,fld_type_p,          &
     &    exner_rho_levels,                                             &
     &    LENRIMA(fld_type_p,halo_type_extended,rima_type_norm),        &
     &    LBC_SIZEA(1,fld_type_p,halo_type_extended,rima_type_norm),    &
     &    LBC_STARTA(1,fld_type_p,halo_type_extended,rima_type_norm),   &
     &    halo_i,halo_j,exner_lbc,RIMWIDTHA(rima_type_norm),            &
     &    RIMWIDTHA(rima_type_norm),RIMWEIGHTSA,                        &
     &    at_extremity,                                                 &
     &    .FALSE.,.TRUE.)

        CALL init_lbc_wghts(row_length,rows,RIMWIDTHA(rima_type_norm),  &
                         halo_i, halo_j, at_extremity, RIMWEIGHTSA )

      END IF ! IF (model_domain  ==  mt_LAM)

! Set reference height profile

! Allocate arrays for level_heights_mod module
      If (.not.allocated(eta_theta_levels)) Then
        Allocate (eta_theta_levels(0:model_levels))
      End If
      If (.not.allocated(eta_rho_levels)) Then
        Allocate (eta_rho_levels(model_levels))
      End If
      If (.not.allocated(r_theta_levels)) Then
        Allocate (r_theta_levels(1-halo_i:row_length+halo_i,           &
     &                           1-halo_j:rows+halo_j, 0:model_levels))
      End If
      If (.not.allocated(r_rho_levels)) Then
        Allocate (r_rho_levels(1-halo_i:row_length+halo_i,             &
     &                         1-halo_j:rows+halo_j, model_levels))
      End If

      eta_theta_levels(0) = eta_theta_in(0)
      z_top_theta = z_top_of_model
      Do k = 1, model_levels
        eta_theta_levels(k) = eta_theta_in(k)
        eta_rho_levels(k) = eta_rho_in(k)
        r_ref_theta(k) = eta_theta_levels(k) * z_top_of_model
      End Do
      Do k = 1, model_levels
        r_ref_rho(k) = eta_rho_levels(k) * z_top_of_model
      End Do
! set bottom level, ie: orography
      Do j = 1-halo_j, rows+halo_j
        Do i= 1-halo_i, row_length+halo_i
          r_theta_levels(i,j,0) = Orog_halos(i,j) + z_offset
        End Do
      End Do
! For constant levels set r to be a constant on the level
      Do k = first_constant_r_rho_level, model_levels
        Do j = 1-halo_j, rows+halo_j
          Do i= 1-halo_i, row_length+halo_i
            r_theta_levels(i,j,k) = z_offset + r_ref_theta(k)
            r_rho_levels(i,j,k)   = z_offset + r_ref_rho(k)
          End Do
        End Do
      End Do

      Select Case( height_gen_method )
        Case( height_gen_original )
! The original version of height generation used in the SI dynamics
!
! For boundary layer levels set depth to be constant.
      Do k = 1, boundary_layer_levels
        Do j = 1-halo_j, rows+halo_j
          Do i= 1-halo_i, row_length+halo_i
            r_theta_levels(i,j,k) = r_theta_levels(i,j,0)               &
     &                                 + r_ref_theta(k)
            r_rho_levels(i,j,k) = r_theta_levels(i,j,0)                 &
     &                               + r_ref_rho(k)
          End Do
        End Do
      End Do
! For intermediate levels use linear relaxation to constant value.
! set orographic heights.
      Do k = boundary_layer_levels+1,                                   &
     &       first_constant_r_rho_level-1
        Do j = 1-halo_j, rows+halo_j
          Do i= 1-halo_i, row_length+halo_i
              r_rho_levels(i,j,k) =                                     &
     &          ( r_rho_levels(i,j,first_constant_r_rho_level) -        &
     &            r_theta_levels(i,j,boundary_layer_levels) ) *         &
     &          ( eta_rho_levels(k) -                                   &
     &            eta_theta_levels(boundary_layer_levels) ) /           &
     &          ( eta_rho_levels(first_constant_r_rho_level) -          &
     &            eta_theta_levels(boundary_layer_levels) )             &
     &           +  r_theta_levels(i,j,boundary_layer_levels)
            r_theta_levels(i,j,k) =                                     &
     &          ( r_rho_levels(i,j,first_constant_r_rho_level) -        &
     &            r_theta_levels(i,j,boundary_layer_levels) ) *         &
     &          ( eta_theta_levels(k) -                                 &
     &            eta_theta_levels(boundary_layer_levels) ) /           &
     &          ( eta_rho_levels(first_constant_r_rho_level) -          &
     &            eta_theta_levels(boundary_layer_levels) )             &
     &           +  r_theta_levels(i,j,boundary_layer_levels)
          End Do
        End Do
      End Do

        Case( height_gen_smooth )
! A smooth quadratic height generation
          Do k = 1, first_constant_r_rho_level-1
            Do j = 1-halo_j, rows+halo_j
              Do i= 1-halo_i, row_length+halo_i

              r_rho_levels(i,j,k) = eta_rho_levels(k) * z_top_of_model +&
     &         z_offset + Orog_halos(i,j) * (1.0 - eta_rho_levels(k)    &
     &              /eta_rho_levels(first_constant_r_rho_level))**2

              r_theta_levels(i,j,k) = eta_theta_levels(k) *             &
     &             z_top_of_model + z_offset + Orog_halos(i,j) *        &
     &             (1.0 - eta_theta_levels(k) /                         &
     &              eta_rho_levels(first_constant_r_rho_level))**2

              End Do
            End Do
          End Do

        Case Default
          icode = 10
          Write (Cmessage,*) 'Unrecognised height generation method - ',&
     &                       'Dump needs to be reconfigured'
          Call Ereport( RoutineName, icode, Cmessage )
      End Select

! ENDGAME-only

CALL eg_swap_bounds(r_theta_levels,tdims_l,fld_type_p,.FALSE.)
CALL eg_swap_bounds(r_rho_levels,pdims_l,fld_type_p,.FALSE.)

! calculate r_at_u points on rho levels and
! r_at_v points on rho levels.
!

! Allocate r_at_u/v arrays for level_heights_mod module
      IF (.NOT. ALLOCATED(r_at_u)) Then
        ALLOCATE (r_at_u (udims_l%i_start:udims_l%i_end,                &
                          udims_l%j_start:udims_l%j_end,                &
                          udims_l%k_start:udims_l%k_end))
      END IF
      IF (.NOT. ALLOCATED(r_at_v)) Then
        ALLOCATE (r_at_v (vdims_l%i_start:vdims_l%i_end,                &
                          vdims_l%j_start:vdims_l%j_end,                &
                          vdims_l%k_start:vdims_l%k_end))
      END IF

      IF (.NOT. ALLOCATED(r_at_u_w)) THEN
        ALLOCATE (r_at_u_w(udims_l%i_start:udims_l%i_end,               &
                           udims_l%j_start:udims_l%j_end,               &
                                 0:model_levels))
      END IF

      IF (.NOT. ALLOCATED(r_at_v_w)) THEN
        ALLOCATE (r_at_v_w(vdims_l%i_start:vdims_l%i_end,               &
                           vdims_l%j_start:vdims_l%j_end,               &
                                  0:model_levels))
      END IF

! NOTE: Simple regular grid is assumed
!
        DO k=1, model_levels
          DO j = udims%j_start, udims%j_end
            DO i = udims%i_start, udims%i_end
              r_at_u(i,j,k) = .5 * (r_rho_levels(i,j,k) +               &
                                    r_rho_levels(i+1,j,k) )
            END DO
          END DO
        END DO
        DO k = 0, model_levels
          DO j = udims%j_start, udims%j_end
            DO i = udims%i_start, udims%i_end
              r_at_u_w(i,j,k) = .5 * (r_theta_levels(i,j,k) +            &
                                      r_theta_levels(i+1,j,k) )
            END DO
          END DO
        END DO

        DO k = 1, model_levels
          DO j = vdims%j_start,vdims%j_end
            DO i = vdims%i_start, vdims%i_end
              r_at_v(i,j,k) = .5 * (r_rho_levels(i,j,k) +                &
                                    r_rho_levels(i,j+1,k) )
            END DO
          END DO
        END DO
        DO k = 0, model_levels
          DO j = vdims%j_start, vdims%j_end
            DO i = vdims%i_start, vdims%i_end
              r_at_v_w(i,j,k) = .5 * (r_theta_levels(i,j,k) +            &
                                      r_theta_levels(i,j+1,k) )
            END DO
          END DO
        END DO

CALL eg_swap_bounds(r_at_u,udims_l,fld_type_u,.FALSE.)
CALL eg_swap_bounds(r_at_v,vdims_l,fld_type_v,.FALSE.)
! DEPENDS ON: swap_bounds
CALL swap_bounds(r_at_u_w,row_length,rows,                            &
                 model_levels+1,halo_i,halo_j,                        &
                 fld_type_u,.FALSE.)
! DEPENDS ON: swap_bounds
CALL swap_bounds(r_at_v_w,row_length,n_rows,                          &
                 model_levels+1,halo_i,halo_j,                        &
                 fld_type_v,.FALSE.)
CALL swap_bounds(exner_rho_levels,row_length,rows,                    &
                 model_levels+1,offx,offy,                            &
                 fld_type_p,.FALSE.)

! End of ENDGAME-only

IF ( PrintStatus == PrStatus_Diag .AND. mype == 0 ) THEN
  write(6,fmt='(A)') '================================================'
  write(6,fmt='(A)') 'level  | r_theta_levls,r_rho_levls,etc @ i=1,j=1'
  write(6,fmt='(A)') '------------------------------------------------'

  write(6,fmt='(I4,A,E15.5)') 0,'   |', r_theta_levels(1,1,0)-earth_radius

  DO k=1,model_levels

    write(6,fmt='(I4,A,6E15.5)') k,'   |',                            &
                              r_theta_levels(1,1,k)-earth_radius,     &
                              r_rho_levels  (1,1,k)-earth_radius,     &
                              r_at_u        (1,1,k)-earth_radius,     &
                              r_at_v        (1,1,k)-earth_radius,     &
                              r_at_u_w        (1,1,k)-earth_radius,   &
                              r_at_v_w        (1,1,k)-earth_radius



  END DO
  write(6,fmt='(A)') '================================================'
END IF

! 0.1 Initialise secondary arrays.
! Exner at p (=rho) levels is obtained from dump.
! calculate p from exner_rho_levels
! [halos required for diagnostic calculations at T+0 in INITDIAG.]

      DO k = 1, model_levels+1
        DO j = 1-offy, rows+offy
          DO i = 1-offx, row_length+offx
            p(i,j,k)= (exner_rho_levels(i,j,k) ** recip_kappa)         &
     &                      * p_zero
          END DO
        END DO
      END DO

! calculate exner at theta_levels which is then used to get
! p at theta_levs
! DEPENDS ON: calc_exner_at_theta
      CALL Calc_Exner_at_theta( r_theta_levels, r_rho_levels,           &
                       exner_rho_levels,                                &
                       row_length, rows, model_levels,                  &
                       offx, offy, halo_i, halo_j,                      &
                       exner_theta_levels,.TRUE.)

! Calculate pressure from Exner at theta levels.

! DEPENDS ON: calc_p_from_exner
      CALL Calc_P_from_Exner(                                           &
                            p_theta_levels,                             &
                            row_length, rows,                           &
                            tdims_s%k_end-tdims_s%k_start+1,            &
                            offx, offy,                                 &
                            exner_theta_levels,.TRUE.)

! calculate p_star using rho (from dump) and p on model levels
! DEPENDS ON: calc_p_star
      CALL Calc_P_star (r_theta_levels, r_rho_levels, p, rho,           &
                        row_length, rows, model_levels,                 &
                        offx, offy, halo_i, halo_j,                     &
                        p_star)

! ----------------------------------------------------------------------
! Initialise Constants and trigonometric functions.
! ----------------------------------------------------------------------

      IF( L_rotating) THEN
        omega = 7.292116E-5 ! Angular speed of Earth's rotation
                            ! = 2*pi/siderial day (23h56m04s)
      ELSE
        omega = 0.0         ! Zero angular speed of rotation for planet
      END IF    ! L_rotating
      two_omega = 2. * omega

! Convert grid-spacing and base values from degrees to radians
      IF( L_regular ) THEN 
        base_lambda_wk  = base_lambda_in
        base_phi_wk     = base_phi_in
        delta_lambda_wk = delta_lambda_in
        delta_phi_wk    = delta_phi_in
      ELSE
        base_lambda_wk  = lambda_p_in(1)
        base_phi_wk     = phi_p_in(1)
        delta_lambda_wk = (lambda_p_in(global_row_length) -             &
                          lambda_p_in(1))/(global_row_length - 1) 
        delta_phi_wk    = (phi_p_in(global_rows) -                      &
                          phi_p_in(1))/(global_rows - 1) 
      END IF 
      delta_lambda = delta_lambda_wk * Pi_over_180
      delta_phi    = delta_phi_wk    * Pi_over_180
      base_phi     = base_phi_wk     * Pi_over_180
      base_lambda  = base_lambda_wk  * Pi_over_180
      f_plane_rad  = f_plane         * Pi_over_180
      ff_plane_rad = ff_plane        * Pi_over_180
      lat_rot_NP   = lat_rot_NP_in   * Pi_over_180
      long_rot_NP  = long_rot_NP_in  * Pi_over_180


! ENDGAME-only

      IF (.NOT.L_idealised_data) THEN
        base_xi1  = base_lambda_in
        base_xi2  = base_phi_in
        delta_xi1 = delta_lambda_in
        delta_xi2 = delta_phi_in
      END IF

      CALL alloc_grid()
      IF ( l_regular ) THEN
        CALL set_horiz_grid( L_Cartesian )
      ELSE 
        CALL set_var_horiz_grid( L_Cartesian,                                 &
                              A_LEN2_COLDEPC,A_LEN2_ROWDEPC,                  &
                              Lambda_p_in, Lambda_u_in, Phi_p_in, Phi_v_in )
      END IF
      CALL calc_global_grid_spacing()

! End of ENDGAME-only

!  Set trig fields and Coriolis components

! DEPENDS ON: set_trigs_4A
      CALL set_trigs_4A(                                                &
                     row_length, rows, n_rows,model_levels,             &
                     base_lambda, base_phi,                             &
                     delta_lambda, delta_phi,                           &
                     base_lambda_wk, base_phi_wk,                       &
                     delta_lambda_wk, delta_phi_wk,                     &
                     f_plane_rad, ff_plane_rad,                         &
                     two_Omega,                                         &
                     lat_rot_NP, long_rot_NP,                           &
                     lat_rot_NP_in, long_rot_NP_in,                     &
                     lat_n, lat_s, long_e, long_w,                      &
                     long_e_model, long_w_model,                        &
                     model_domain, rmdi)

! ----------------------------------------------------------------------
! 1. Set polar filtering and diffusion
! ----------------------------------------------------------------------

! Allocate arrays for diff_coeff_mod module
      IF (.NOT.ALLOCATED(diff_coeff_u)) THEN
        ALLOCATE (diff_coeff_u ( 1-Offx : row_length+Offx,              &
                                 1-Offy : rows+Offy ))
      END IF
      IF (.NOT.ALLOCATED(diff_coeff_v)) THEN
        ALLOCATE (diff_coeff_v ( 1-Offx : row_length+Offx,              &
                                 1-Offy : n_rows+Offy ))
      END IF
!
!
!  THIS NEEDS OUT_OF_BOUNDS ERROR FIXED FOR ENDGAME
!
!
     IF ( L_pofil_new ) THEN
  WRITE(6,'(A)') ' Newest combi-polar filter used '
! DEPENDS ON: eg_setdiff
        Call eg_Setdiff(                                                &
! Input size and control variables
            global_rows, row_length, rows, n_rows, model_levels,        &
            model_domain, at_extremity, datastart,                      &
            offx, offy, mype, nproc_y,                                  &
            model_levels_max, max_121_rows, max_sweeps,                 &
! other constants
            delta_lambda, delta_phi,                                    &
            polar_cap, scale_ratio,                                     &
            ref_lat_deg, diff_coeff_ref,                                &
            cos_theta_latitude, sin_theta_latitude,                     &
            cos_v_latitude, sin_v_latitude,                             &
! Output data
             global_u_filter, global_v_filter,                          &
             u_sweeps, v_sweeps,                                        &
             u_begin, u_end, v_begin, v_end,                            &
             diff_coeff_u, diff_coeff_v, diff_coeff_phi,                &
             diffusion_coefficient_thermo, diffusion_order_thermo,      &
             diffusion_coefficient_wind, diffusion_order_wind,          &
             diff_order_thermo, diff_order_wind,                        &
             diff_timescale_thermo, diff_timescale_wind,                &
             L_sponge, sponge_power, sponge_ew, sponge_ns,              &
             sponge_wts_ew, sponge_wts_ns, L_subfilter_horiz,           &
             L_diffusion, L_cdiffusion, L_filter, L_diff_auto,          &
             L_pfcomb, L_pfexner, L_pftheta, L_pfuv, L_pfw, L_pfincs,   &
             L_diff_thermo, L_diff_wind, L_diff_w, L_diff_incs,         &
             halo_i, halo_j,xi2_p, xi2_v                              )
     ELSE
       WRITE(6,'(A)') 'WARNING - '
       WRITE(6,'(A)') '   only New combi-filter is compatible with ENDGAME'
     END IF

! ----------------------------------------------------------------------
! 1.5 Set filtering control variable
! ----------------------------------------------------------------------
       L_diff_active = .FALSE.

      If ( L_polar_filter .or. L_pfcomb ) Then
        L_diff_active = .TRUE.
        L_filter = .true.
        if ( model_domain  /=  mt_global) then
          write(6,*)' ******   Polar filter NOT ALLOWED in LAMs ****** '
          write(6,*)' ******       RESET UMUI or NAMELIST       ****** '
          ICODE    = 1
          CMESSAGE ='SETCONA:Polar filter NOT ALLOWED in LAMs'
          CALL Ereport(RoutineName,ICODE,Cmessage)
        else if ( L_pofil_new ) then
          write(6,*)' ******  Newest filtering code being used ****** '
        end if ! model_domain  /=  mt_global
      Else
        write(6,*) ' ******   Polar filter is OFF   ****** '
      End If ! L_polar_filter .or. L_pfcomb
      If ( L_diff_thermo .or. L_diff_wind .or. L_diff_w ) Then
        L_diff_active = .TRUE.
        if ( L_pofil_new ) then
          write(6,*) ' ******  Newest diffusion code being used ****** '
        endif !  L_pofil_new
        L_filter = .true.
      End If ! L_diff_thermo .or. L_diff_wind .or. L_diff_w
      If ( L_diff_incs .or. L_pfincs) Then
        L_filter_incs = .true.
      End If ! L_diff_incs .or. L_pfincs
      If ( L_pfexner  ) Then
        L_diff_active = .TRUE.
        write(6,*) ' ** Polar filtering of Exner pressure is active ** '
      End If ! L_pfexner
      If ( L_diff_exner  ) Then
        L_diff_active = .TRUE.
        write(6,*) ' ** Diffusion of Exner pressure is active ** '
      End If ! L_diff_exner

      If( L_tardiff_q ) Then
        write(6,*) ' '
        write(6,*) '  ***   Targeted diffusion of q is ACTIVE  *** '
        write(6, 913) w_conv_limit, tardiffq_factor
        write(6, 914) tardiffq_test
        write(6, 915) tardiffq_start, tardiffq_end
        write(6, 916) tar_horizontal
      Else
        write(6,*) ' '
        write(6,*) '  ***   Targeted diffusion of q is NOT ACTIVE  *** '
      End If !   L_tardiff_q

! ----------------------------------------------------------------------
! Section 1.6 Set up turb_diff control
! ----------------------------------------------------------------------
      IF ( L_subfilter_horiz .OR. L_subfilter_vert ) THEN

        L_diff_active = .TRUE.

! DEPENDS ON: init_turb_diff
        CALL init_turb_diff(                                            &
                            model_levels, bl_levels,                    &
                            turb_startlev_horiz, turb_endlev_horiz,     &
                            turb_startlev_vert, turb_endlev_vert,       &
                            global_u_filter, global_v_filter,           &
                            diff_factor, mix_factor,                    &
                            L_subfilter_horiz, L_subfilter_vert,        &
                            L_subfilter_blend )

      END IF ! L_subfilter_horiz .OR. L_subfilter_vert

!L ----------------------------------------------------------
!  Print information to tell user what is active
!L ----------------------------------------------------------
!  vertical diffusion coefficients
      If ( L_vdiff_uv ) Then
        vdiffuv_factor = 0.25 *                                         &
     &                     (1.0 - EXP( -1.0 / vdiffuv_timescale) )
       Write ( Unit = 6, fmt=*) 'A factor delta_z**2/timestep is '      &
     & ,' allowed for in the application of the vertical diffusion.'
       Write ( Unit = 6, fmt=*) 'Vertical diffusion of order 1 with'    &
     &        , ' vertical diffusion coefficient = ', vdiffuv_factor
       Write ( Unit = 6, fmt=*) 'Vertical diffusion timescale is '      &
     &                           , vdiffuv_timescale,' steps '
      End If ! L_vdiff_uv

      if ( top_filt_start > model_levels ) then
        Write ( Unit = 6, fmt=*) 'No additional upper-level diffusion'
      else  ! top_filt_start <= model_levels
        active_levels = top_filt_end - top_filt_start + 1
        if ( active_levels > max_updiff_levels) then
          top_filt_start  = top_filt_end - max_updiff_levels + 1
          write(6,*) ' Max of uppermost ', max_updiff_levels            &
     &          , ' levels allowed for upper-level'                     &
     &          , ' diffusion - start level reset to ', top_filt_start
          active_levels = max_updiff_levels
        end if ! active_levels > max_updiff_levels
        Write ( Unit = 6, fmt=*) 'Extra upper-level diffusion applied'
        scale = 1.0
        if ( L_upper_ramp ) then
          scale = up_diff_scale
          write(6, 911) top_diff, top_filt_end
          write(6, 912) scale, top_filt_start
        else  !
          write(6, 910) top_diff, top_filt_start, top_filt_end
        end if  !  L_upper_ramp
        kk = active_levels
        up_diff(active_levels) = top_diff
        do k = top_filt_end - 1, top_filt_start, - 1
          kk = kk - 1
          up_diff(kk) = scale * up_diff(kk + 1)
        enddo !  k = top_filt_end - 1, top_filt_start, - 1
        Do k = 1, active_levels
          kk = k + top_filt_start - 1
          write(6,*)'Level', kk,' Diffusion factor = ',up_diff(k)
        End do !  k = 1, active_levels
      end if ! top_filt_start > model_levels

      if (L_adjust_theta) then
        if( adjust_theta_start < 2 ) then
          adjust_theta_start = 2
          Write ( Unit = 6, fmt=*) 'Start level for convective '        &
     & , ' adjustment reset to ', adjust_theta_start
        end if ! adjust_theta_start < 2
        Write ( Unit = 6, fmt=*) 'Convective adjustment applied when'   &
     &    ,' lapse rate is negative above level ',  adjust_theta_start  &
     &    ,' up to level ', adjust_theta_end
      end if  ! L_adjust_theta

! ----------------------------------------------------------------------
! Section 1.5a . Print FORMATTING
! ----------------------------------------------------------------------

 910  FORMAT(' Diffusion coefficient = ',F5.2                           &
     &      ,' from level ', I3,' to level ', I3)
 911  FORMAT(' Diffusion coefficient = ',F5.2,' at level ', I3)
 912  FORMAT(' decreasing by ', F6.3, ' at each level down to ', I3)
 913  FORMAT('   Vertical velocity test value is  ',F5.2,' m/s and ',   &
     &                                 'diffusion factor = ',F6.3)
 914  FORMAT('   Start testing at level  ',I3)
 915  FORMAT('   Apply from level  ',I3,' to level ',I3)
 916  FORMAT('   Slope test applied up to level ',I4)

! ----------------------------------------------------------------------
! 2. Semi-Lagrangian scheme options.
! ----------------------------------------------------------------------
! j holds number of points required for interpolation scheme
      j = 2
      Do i = 1, 3
        If (high_order_scheme(i)  ==  quinticLagrange) Then
          j = 3
          If (halo_j  <   5) Then
            write(6,*)' error halo_j too small ',halo_j
          End If
        Else
          If (halo_j  <   4) Then
            write(6,*)' error halo_j too small ',halo_j
          End If
        End If
      End Do

! Calculate max cfl possible in LAM
! Y-direction also applies in global model though variable is not
! used in code.

      LAM_max_cfl(1) = halo_i - j
      LAM_max_cfl(2) = halo_j - j

! ----------------------------------------------------------------------
! Set up data.
! ----------------------------------------------------------------------

! Calculate land_index
          land_points = 0
          Do j =1, rows
            Do i = 1, row_length
              If (land_sea_mask(i,j)) then
                land_points = land_points + 1
                land_index(land_points) = (j-1)*row_length + i
              End If
            End Do
          End Do

! set-up code for hydrology
          soil_points = 0
          land_ice_points = 0
          Do i= 1, land_points
! Test on soil moisture concentration at saturation
            If (smvcst(i) >   0.0) Then       ! Soil points
              soil_points = soil_points+1
              soil_index(soil_points)=i
            Else If (smvcst(i)  ==  0.0) Then   ! Land-ice points
              land_ice_points = land_ice_points+1
              land_ice_index(land_ice_points)=i
            End If
          End Do

! call swap_bounds to set extra points
      i_field = 0
      i_field = i_field + 1
      fields_to_swap(i_field) % field      => r_at_u(:,:,:)
      fields_to_swap(i_field) % field_type = fld_type_u
      fields_to_swap(i_field) % levels     = model_levels
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field      => r_at_v(:,:,:)
      fields_to_swap(i_field) % field_type = fld_type_v
      fields_to_swap(i_field) % levels     = model_levels
      fields_to_swap(i_field) % rows       = n_rows
      fields_to_swap(i_field) % vector     = .FALSE.

! DEPENDS ON: swap_bounds_mv
       CALL swap_bounds_mv(fields_to_swap, i_field, row_length,         &
                              halo_i, halo_j)

      If ( L_diag_L2norms ) then
        Write ( Unit = 6, fmt=*) 'Printing of norms during timestep'    &
     &, ' activated'
      endIf ! L_diag_L2norms
      If ( L_diag_L2helm ) then
        Write ( Unit = 6, fmt=*) 'Printing of coefficient norms  '      &
     &, ' in solver'
      endIf ! L_diag_L2helm
!  Initialise diagnostic printing items
! Sampling interval must not be larger than print interval
      if ( diag_interval > print_step) then
        diag_interval = print_step
      end if !  diag_interval > print_step
! Set array sizes needed for chosen diagnostic printing
! These are needed to do sums and max/mins across processors
      rpemax = 0
      rpemin = 0
      ipesum = 0
      rpesum = 0
      if (L_print_theta1 ) then
        rpemin = rpemin + 1
        ipesum = ipesum + 1
        rpesum = rpesum + 3
      end if !  L_print_theta1
      if (L_print_lapse ) then
        rpemin = rpemin + (model_levels - 1)
        ipesum = ipesum + 2 * (model_levels - 1)
        rpesum = rpesum + 2 * (model_levels - 1)
      end if !  L_print_lapse
      if (L_print_div) then
        rpemax = rpemax + model_levels
        rpemin = rpemin + model_levels
        ipesum = ipesum + 2 * model_levels
        rpesum = rpesum + 4 * model_levels
      end if !  L_print_div
      if (L_print_w .or. L_print_wmax) then
        rpemax = rpemax + model_levels
        ipesum = ipesum + 2 * (model_levels)
        rpesum = rpesum + 2 * (model_levels)
      end if !  L_print_w .or. L_print_wmax
      if (L_print_shear) then
        rpemax = rpemax + model_levels - 1
        ipesum = ipesum + model_levels - 1
        rpesum = rpesum + 2 * (model_levels - 1)
      end if !  L_print_shear
      if (L_print_max_wind) then
        rpemax = rpemax + model_levels
        ipesum = ipesum + model_levels
        rpesum = rpesum + 2 * model_levels
      end if !  L_print_max_wind
      if ( L_diag_wind ) then
        rpesum = rpesum + 2 * model_levels
      end if  ! L_diag_wind
        min_theta1_run = 1000.0
        time_theta1_min = 0
      max_w_run(0) =  0.0
      Do k = 1, model_levels
         max_w_run(k) =  0.0
         time_w_max(k) = 0
         max_div_run(k) =  0.0
         time_div_max(k) = 0
         min_div_run(k) =  0.0
         time_div_min(k) = 0
        min_lapse_run(k) =  1.0e6
         time_lapse_min(k) = 0
        max_shear_run(k) = 0.0
        time_max_shear(k) = 0
        max_wind_run(k) = 0.0
        time_max_wind(k) = 0
        max_KE_run(k) = 0.0
        min_KE_run(k) = 1.0e30
        time_KE_max(k) = 0
        time_KE_min(k) = 0
        time_noise_max(k) = 0
        max_noise_run(k) = 0.0
      End Do  ! k = 1, model_levels
      max_KE_run(model_levels + 1) = 0.0
      min_KE_run(model_levels + 1) = 1.0e30
      time_KE_max(model_levels + 1) = 0
      time_KE_min(model_levels + 1) = 0

!  Set values for block printing:
!               calls for printing need to be inserted by user
! DEPENDS ON: init_block4_pr
        CALL init_block4_pr(                                            &
     &                      global_row_length, global_rows,             &
     &                      blockx_in, blocky_in,                       &
     &                      dom_w_in, dom_e_in, dom_s_in, dom_n_in)
 

!L ----------------------------------------------------------
!L 3. Set up cloud type boundaries for low/medium/high cloud.
!L ----------------------------------------------------------
       IF (NUM_CLOUD_TYPES  >   3) THEN
         ICODE    = 1
         CMESSAGE = 'SETCONA: Parameter NUM_CLOUD_TYPES exceeds 3'
         WRITE(6,*) 'NUM_CLOUD_TYPES=',NUM_CLOUD_TYPES
         CALL Ereport(RoutineName,ICODE,Cmessage)

       END IF

!  Diagnostics of cloud amount for low/medium/high cloud are calculated
!  by finding the max cloud cover over a range of model levels. The
!  ranges are determined by heights [ 1->2:2->3:3->4 = L:M:H ] held in
!  array h_split by comparison with heights of model levels (above
!  surface) in r_ref_theta.

      DO KK = 1, NUM_CLOUD_TYPES + 1
        LEVEL = 1
!
!       r_ref_theta turns out to be A_theta(k) - Earth_radius at every
!       model level because z_top_of_model = z_rho_top is chosen here.
!
        DO WHILE ((r_ref_theta(LEVEL)  <=  h_split(KK)) .AND.           &
     &                                       (LEVEL  <=  model_levels))
          LEVEL = LEVEL + 1

        END DO

        IF (LEVEL  >   model_levels) THEN
          ICODE    = -1
          CMESSAGE ='SETCONA:Error in locating levels for cloud layers'
          CALL Ereport(RoutineName,ICODE,Cmessage)

        END IF
        CLOUD_BOUND(KK) = LEVEL

      END DO

      LOW_BOT_LEVEL  = CLOUD_BOUND(1)
      LOW_TOP_LEVEL  = CLOUD_BOUND(2) - 1
      MED_BOT_LEVEL  = CLOUD_BOUND(2)
      MED_TOP_LEVEL  = CLOUD_BOUND(3) - 1
      HIGH_BOT_LEVEL = CLOUD_BOUND(3)
      HIGH_TOP_LEVEL = CLOUD_BOUND(4) - 1

      IF (LOW_TOP_LEVEL  >   CLOUD_LEVELS) THEN
        ICODE    = 1
        CMESSAGE = 'SETCONA: No of cloud levels less than Top of Low'
        CALL Ereport(RoutineName,ICODE,Cmessage)

      END IF

      IF (MED_TOP_LEVEL >  CLOUD_LEVELS) THEN
        ICODE    = 1
        CMESSAGE = 'SETCONA:  No of cloud levels less than Top of Med'
        CALL Ereport(RoutineName,ICODE,Cmessage)

      END IF

      IF (HIGH_TOP_LEVEL  >   CLOUD_LEVELS) THEN
        ICODE    = 1
        CMESSAGE = 'SETCONA: No of cloud levels less than Top of High'
        CALL Ereport(RoutineName,ICODE,Cmessage)

      END IF

!L ----------------------------------------------------------
!L 4. Set up locally-derived parameters.
!L ----------------------------------------------------------

!     The tropopause diagnosed for radiative purposes
!     divides theta-levels considered to lie in the stratosphere
!     from those considered to lie in the troposphere: the
!     tropopause is therefore taken as a rho-level. This level
!     is constrained to lie between heights of z_min_trop
!     and z_max_trop. The level is used in the setting of the
!     background aerosol climatology and of ozone profiles,
!     subject to the choice of appropriate options; additionally
!     it is used in the calculation of diagnostics defined at
!     the tropopause.
!
!     Start at the second rho-level because the first is omitted
!     from the grid seen in the physics.
      min_trop_level=2
      Do ; If ( (r_ref_rho(min_trop_level) >= z_min_trop) .OR.          &
     &          (min_trop_level == model_levels) ) Exit
        min_trop_level = min_trop_level+1
      End Do
!
      max_trop_level=min_trop_level
      Do ; If ( (r_ref_rho(max_trop_level) > z_max_trop) .OR.           &
     &          (max_trop_level == model_levels) ) Exit
        max_trop_level = max_trop_level+1
      End Do
      max_trop_level = max_trop_level-1

!L ----------------------------------------------------------
!L 5.0 Set up tracer advection parameters
!L ----------------------------------------------------------

! set up super tracer array size

          super_array_size=0

          If (l_CO2_interactive) super_array_size = super_array_size+1
          If (l_Soot) super_array_size = super_array_size+3
          If (l_Biomass) super_array_size = super_array_size+3
          If (l_Sulpc_so2)super_array_size = super_array_size+4
          If (L_sulpc_nh3)super_array_size = super_array_size+1
          If (L_sulpc_dms)super_array_size = super_array_size+1
          If (L_DUST) THEN 
            IF (l_twobin_dust) THEN 
              super_array_size = super_array_size + 2 
            ELSE 
              super_array_size = super_array_size + 6 
            END IF 
          END IF 
          If (L_ocff) super_array_size = super_array_size+3
          IF (L_nitrate) super_array_size=super_array_size+2
          If (L_Murk_advect) super_array_size = super_array_size + 1
          If (L_USE_CARIOLLE) super_array_size = super_array_size + 1

! ----------------------------------------------------------------------
! Section 5.1   limit for tracer_level required
! ----------------------------------------------------------------------

       If ( tr_levels == model_levels ) Then

         If ( tr_vars > 0 ) super_array_size = super_array_size + tr_vars

         If ( tr_ukca > 0 ) super_array_size = super_array_size + tr_ukca

       Else   !  tr_levels /= model_levels
       
         If (tr_ukca > 0 .or. tr_vars > 0 ) Then
       
           WRITE(6,*) 'Tracer Levels: tr_levels =',tr_levels 
           WRITE(6,*) 'Model Levels: model_levels =', model_levels
           ICODE = 1
           CMESSAGE='Tracer levels must equal model levels, when '  //  &
     &              'running with Free or UKCA tracers'
 
           Call Ereport(RoutineName, ICODE, CMESSAGE)

         End If   ! tr_ukca > 0 .or. tr_vars > 0

       End If ! tr_levels == model_levels

       If (mype == 0) write(6,*)'tracer super_array_size =  ',          &
     &                          super_array_size

!L ----------------------------------------------------------
!L 6. Set up moisture array sizes    always do q,qcl,qcf
!L ----------------------------------------------------------

          moisture_array_size=3
          if (l_pc2) moisture_array_size=moisture_array_size+4
          if (l_mcr_qrain)moisture_array_size=moisture_array_size+1
          if (l_mcr_qcf2 )moisture_array_size=moisture_array_size+1
          if (l_mcr_qgraup)moisture_array_size=moisture_array_size+1

          if(mype == 0)write(6,*)'moisture_array_size =  ',             &
     &                          moisture_array_size

! ----------------------------------------------------------------
! 7. Calculate the coeffs for the interpolation of the radiation
!    quantities if spatial degradation of the radiation code is on.
!    Calculate the mask which ensures a chequer-board pattern over
!    the whole domain (and not just for a single PE).
! -----------------------------------------------------------------

! Create a land mask with a halo.
      IF (L_rad_deg) THEN

        DO j=1,ROWS
          DO i=1,ROW_LENGTH
            landmask(i,j)=Land_sea_mask(i,j)
          END DO
        END DO

! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(landmask,row_length,rows,1,                    &
     &                   offx,offy,fld_type_p,.FALSE.)

        first_row=1
        last_row=rows
        IF (model_domain  ==  mt_global) THEN
          IF (at_extremity(PNorth)) THEN
            last_row=rows-1
          END IF
          IF (at_extremity(PSouth)) THEN
            first_row=2
          END IF
        END IF

! Allocate arrays for rad_mask_trop_mod module
      If (.not.allocated(es_space_interp)) Then
        Allocate (es_space_interp(4, row_length, rows))
      End If
      If (.not.allocated(rad_mask)) Then
        Allocate (rad_mask(row_length, rows))
      End If

! DEPENDS ON: coeffs_degrade
        CALL COEFFS_DEGRADE(ES_SPACE_INTERP, landmask,                  &
     &                      first_row, last_row, model_domain,          &
     &                      mt_lam, at_extremity(PNorth),               &
     &                      at_extremity(PSouth), at_extremity(PWest),  &
     &                      at_extremity(PEast),                        &
     &                      row_length*rows,row_length, rows,offx, offy)
! DEPENDS ON: rad_degrade_mask
        CALL RAD_DEGRADE_MASK(RAD_MASK, datastart,                      &
     &                      Ndim_max, first_row, last_row, offx,        &
     &                      row_length, rows)

      ELSE
        If (.not.allocated(es_space_interp)) Then
          Allocate (es_space_interp(1,1,1))
        End If
        If (.not.allocated(rad_mask)) Then
          Allocate (rad_mask(1,1))
        End If
      END IF  !  If L_rad_deg loop

        If ( thmono_height >= 0.5 * r_ref_theta(1) ) then
          If ( thmono_height < r_ref_theta(1) ) then
            k = 1
          else if ( thmono_height > r_ref_theta(model_levels) ) then
            k = model_levels
          else
            k = 1
            do  !  cycle to find nearest level to  thmono_height
              k = k + 1
              if ( thmono_height < r_ref_theta(k) ) then
                if ( r_ref_theta(k) - thmono_height >                   &
     &               thmono_height - r_ref_theta(k-1) ) k = k - 1
                EXIT
              end if ! thmono_height < r_ref_theta(k)
              CYCLE
            end do !  cycle to find nearest level to  thmono_height
          End If ! thmono_height < r_ref_theta(1)
          thmono_levels = k
          Write ( Unit = 6, fmt=*) '3D Monotone limiter will be applied'&
     &,                            ' to advection of theta'
          Write ( Unit = 6, fmt=*) 'thmono_height was set to '          &
     &,                             thmono_height,' in the UMUI'
          Write ( Unit = 6, fmt=*) 'Limiter will be applied up to level'&
     &,   thmono_levels,', the nearest level to',thmono_height,'metres'
        Else   ! thmono_height < 0.5 * r_ref_theta(1)
          thmono_levels = 0
        End If ! thmono_height >= 0.5 * r_ref_theta(1)
!L ----------------------------------------------------------
!L 8. Set up dynamical core parameters
!L ----------------------------------------------------------
!    Set frictional_timescale = 0 everywhere since this value
!    is passed into PE_HELMHOLTZ from ATMSTEP
!    The value used in simple friction is calculated internally
!           ( stored in friction_level)
!     as is SuHe_level_weight (temp1) in the temperature relaxation

      Do k = 1, model_levels
        frictional_timescale(k) = 0.0
      End Do

      If( problem_number  ==  dynamical_core) Then

! Note: p_theta_levels is passed into IDL_Set_SuHe_params with a levels
!       mismatch - but it is not used, so it doesn't matter
! DEPENDS ON: idl_set_suhe_params
        Call IDL_Set_SuHe_params                                        &
     &                      (row_length, rows, model_levels             &
     &,                     offx, offy, halo_i, halo_j                  &
     &,                     mype, nproc, at_extremity                   &
     &,                     datastart, gc_all_proc_group                &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     SuHe_pole_equ_deltaT, SuHe_static_stab      &
     &,                     p, p_theta_levels                           &
     &,                     base_frictional_timescale                   &
     &,                     friction_level, SuHe_sigma_cutoff           &
     &,                     SuHe_level_weight )

      else if( problem_number  ==  idealised_problem) Then
!     Clear arrays
        Do k = 1, model_levels
          SuHe_level_weight(k) = 0.0
          friction_level(k) = 0.0
        End Do    ! k = 1, model_levels

      end if    ! problem_number  ==  dynamical_core
! ----------------------------------------------------------------------
! 9. GCR diagnostics initialisation.
! ----------------------------------------------------------------------
      if (GCR_diagnostics == 3 )then
        L_error = .false.
        Do i = 1, 3
          if(GCR_its_avg_step(i) == 0 )L_error = .true.
        EndDo  !  i = 1, 3
        if ( L_error ) then
          Write ( Unit = 6, fmt=*) 'WARNING GCR iteration counting at'  &
     &,         ' timestep 0 or interval of 0 NOT PERMITTED '
! Following values will output iteration info at 6, 12 and
!  1440 timesteps (30 days for a 30 minute timestep)
          GCR_its_avg_step(1) = 6
          GCR_its_avg_step(2) = 12
          GCR_its_avg_step(3) = 1440
          write(6,*) ' Iteration count diagnostics reset for timesteps '&
     &          , GCR_its_avg_step(1), GCR_its_avg_step(2)              &
     &                               , GCR_its_avg_step(3)
          write(6,*)  '  and at intervals of ', GCR_its_avg_step(3),    &
     &            ' timesteps thereafter.'
          write(6,*)  ' Change in UMUI if you want different values '
        else ! L_error .false.
          write(6,*)  ' '
          write(6,*)  'Iteration count diagnostics at timesteps '       &
     &          , GCR_its_avg_step(1), GCR_its_avg_step(2)              &
     &          , GCR_its_avg_step(3)
          write(6,*)  ' and at intervals of ', GCR_its_avg_step(3),     &
     &            ' timesteps thereafter.'
        endif ! L_error
        GCR_its_switch =  1
        GCR_sum_its = 0
        GCR_min_its = 1000
        GCR_max_its = 0
        GCR_max_time = 0
        GCR_min_time = 0
      endif ! GCR_diagnostics == 3

! ----------------------------------------------------------------------
! 10. Error trapping for advection choices
! ----------------------------------------------------------------------

      do i=1,3
        if (.not. L_mono(i) .and. .not. L_high(i)) then
          if(i == 1)write ( Unit = 6, fmt=*)                            &
     &      'WARNING Advection choices incompatible for theta'
          if(i == 2)write ( Unit = 6, fmt=*)                            &
     &      'WARNING Advection choices incompatible for moisture'
          if(i == 3)write ( Unit = 6, fmt=*)                            &
     &      'WARNING Advection choices incompatible for winds'
          Write ( Unit = 6, fmt=*)                                      &
     &      'Both L_high and L_mono switches set to false '

        endif
      end do

! ----------------------------------------------------------------------
! 11. Orographic correction parameters for the radiation code.
! ----------------------------------------------------------------------

      IF (l_orog) THEN

!       Calculate (rotated) grid coordinates with a halo out to the
!       horizon limit, to be used in calculation of gridpoint separation.
        IF (l_regular) THEN
          DO j = 1-horiz_limit, rows+horiz_limit
            gj = datastart(2) + j - 1
            phi_horiz(j) = base_phi+(gj-1)*delta_phi
          END DO
          dphi_horiz=delta_phi
          dlambda_horiz=delta_lambda
        ELSE
          DO i = 1-horiz_limit, row_length+horiz_limit
            gi = MAX(MIN(datastart(1)+i-1,global_row_length-1),1)
            dlambda_horiz(i) = (lambda_p_in(gi+1)-lambda_p_in(gi))      &
              * pi_over_180
          END DO
          DO j = 1-horiz_limit, rows+horiz_limit
            gj = MAX(MIN(datastart(2)+j-1,global_rows-1),1)
            dphi_horiz(j) = (phi_p_in(gj+1)-phi_p_in(gj))*pi_over_180
          END DO
          DO j = 1-horiz_limit, rows+horiz_limit
            gj = datastart(2) + j - 1
            IF (gj >= 1 .AND. gj <= global_rows) THEN
              phi_horiz(j) = phi_p_in(gj)*pi_over_180
            ELSE IF (gj < 1) THEN
              phi_horiz(j) = phi_p_in(1)*pi_over_180 +                  &
                (gj-1)*dphi_horiz(j)
            ELSE
              phi_horiz(j) = phi_p_in(global_rows)*pi_over_180 +        &
                (gj-global_rows)*dphi_horiz(j)                 
            END IF
          END DO
        END IF


!       Calculate local bearing of the rotated pole.
        IF (model_domain == mt_global) THEN
          bear_rot_NP = 0.0
        ELSE
! DEPENDS ON: pole_bearing
          CALL pole_bearing(row_length, rows,                           &
            lat_rot_NP, long_rot_NP, bear_rot_NP)
        END IF


        IF (l_skyview) THEN
!         Set orographic heights for calculation of horizon angles
!         (using unsmoothed height data if available).
          orog_horiz=0.0
          IF (l_orog_unfilt) THEN
! DEPENDS ON: all_gather_field
            CALL all_gather_field(orog_unfilt, orog_global,             &
              row_length, rows, global_row_length, global_rows,         &
              fld_type_p, halo_type_no_halo, gc_all_proc_group,         &
              icode, cmessage)
          ELSE
            CALL all_gather_field(orog, orog_global,                    &
              row_length, rows, global_row_length, global_rows,         &
              fld_type_p, halo_type_no_halo, gc_all_proc_group,         &
              icode, cmessage)
          END IF
          DO j = 1-horiz_limit, rows+horiz_limit
            gj = datastart(2) + j - 1
            IF (gj >= 1 .AND. gj <= global_rows) THEN
              DO i = 1-horiz_limit, row_length+horiz_limit
                gi = datastart(1) + i - 1
                IF (gi >= 1 .AND. gi <= global_row_length) THEN
                  orog_horiz(i,j) = orog_global(gi,gj)
                ELSE IF (model_domain == mt_global) THEN
                  IF (gi < 1) THEN
                    orog_horiz(i,j)=orog_global(gi+global_row_length,gj)
                  ELSE
                    orog_horiz(i,j)=orog_global(gi-global_row_length,gj)
                  END IF
                END IF
              END DO
            END IF
          END DO
        END IF


!       Calculate slope aspect and angle.
        IF (l_use_grad_corr) THEN
!         Use ancillary slopes calculated from the raw orography.
! DEPENDS ON: aspang_ancil
          CALL aspang_ancil(row_length, rows, land_points,              &
            land_sea_mask, grad_x, grad_y, bear_rot_NP)
        ELSE IF (l_skyview) THEN
!         Use the same orographic heights as used for the
!         calculation of the horizon angles.
! DEPENDS ON: aspang
          CALL aspang(row_length, rows, bear_rot_NP,                    &
            phi_horiz, dphi_horiz, dlambda_horiz,                       &
            orog_horiz(0:row_length+1,0:rows+1))
        ELSE
!         Use the model orography (generally smoothed).
          CALL aspang(row_length, rows, bear_rot_NP,                    &
            phi_horiz, dphi_horiz, dlambda_horiz,                       &
            orog_halos(0:row_length+1,0:rows+1))
        END IF

        IF (model_domain == mt_global) THEN
           IF (at_extremity(PNorth)) slope_angle(:,rows) = 0.0
           IF (at_extremity(PSouth)) slope_angle(:,1) = 0.0
        END IF


        IF (l_skyview) THEN
!         Calculate horizon angles and sky-view correction factor.
! DEPENDS ON: skyview
          CALL skyview(row_length, rows, bear_rot_NP,                   &
            phi_horiz, dphi_horiz, dlambda_horiz, orog_horiz)
        END IF

      END IF

!     Calculate the area within each gridbox. Note that this calculation
!     assumes that the length of each gridbox in the North-South
!     direction is fixed, so it should not be used in variable
!     resolution models.
      yscale=delta_phi*earth_radius
      DO k=1,rows
        DO j=1,row_length 
          xscale(j,k)=delta_lambda*earth_radius*fv_cos_theta_latitude(j,k)
        END DO
      END DO
      WHERE (xscale == 0.0) xscale=EPSILON(xscale)

      ALLOCATE(gridbox_size(row_length,rows))
      DO k=1,rows
        DO j=1,row_length 
          gridbox_size(j,k)=0.001*SQRT(xscale(j,k)*yscale)
        END DO
      END DO

! ----------------------------------------------------------------------
! 12. Boundary layer solver options
! ----------------------------------------------------------------------
      If ( L_us_blsol ) Then
        If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
          write(6,*)
          write(6,*)  '*** New stable and non-oscillatory'//&
              ' boundary-layer solver is ACTIVE ***'
          write(6,*)  '     It will run with Puns=',Puns,'Pstb=',Pstb
        End If 
      Else
        If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
          write(6,*)
          write(6,*)'*** Standard boundary-layer solver is ACTIVE ***'
        End If 
      End If

IF (integrity_test)                                                   &
  CALL update_hash_m(                                                 &
                  eta_theta_levels,SIZE(eta_theta_levels),'etatl',  &
                  eta_rho_levels,  SIZE(eta_rho_levels),  'etarl')


      IF (lhook) CALL dr_hook('SETCONA_4A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Setcona_4A
