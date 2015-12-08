! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      SUBROUTINE CHECK_IDEALISE_4A                                      &
     &  ( LEN_INTHD, LEN_REALHD, LEN_LEVDEPC1, LEN_LEVDEPC2,            &
     &  INTHD, REALHD, LEVDEPC )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE filenamelength_mod, ONLY:                                     &
          filenamelength

      USE Control_Max_Sizes
      USE UM_ParVars
      USE Ereport_mod, ONLY : Ereport      
      USE PrintStatus_mod

      USE eg_alpha_mod
      USE eg_alpha_ramp_mod

      USE ref_pro_mod ! needed due to the inclusion of cruntimc
      USE eg_parameters_mod  ! needed due to the inclusion of cruntimc

      USE horiz_grid_mod, ONLY :  delta_xi1, delta_xi2, base_xi1,       &
                                  base_xi2,                             &
                                  Nxi1L, Nxi1V, Nxi2L, Nxi2V,           &
                                  delta_xi1_H, delta_xi1_L,             &
                                  delta_xi2_H, delta_xi2_L

      USE domain_params
      USE problem_mod, ONLY: standard, monsoon, dynamical_core,         &
                             idealised_problem, standard_namelist
      USE qprofile_mod, ONLY: qp_dry, qp_qsat, qp_namelist_rh,          &
                              qp_namelist, qp_dump
      USE tprofile_mod, ONLY: tp_dthetadz, tp_isothermal, tp_bruntv,    &
                              tp_bv_isoth, tp_dyn_core, tp_dyn_core_lam,&
                              tp_namelist, tp_dump
      USE uvhoriz_mod, ONLY: uv_horiz_const, uv_horiz_ramp,             &
                             uv_horiz_balance, uv_horiz_deform,         &
                             uv_vert_const, uv_vert_interp,             &
                             uv_vert_namelist, uv_vert_dump
      
      USE sl_input_mod, ONLY:  Instability_diagnostics
                           
      USE dynamics_testing_mod, ONLY: L_idealised_data 
      USE check_iostat_mod
                        
      USE eg_idl_baro_mod, ONLY: baro_phi0   
      USE um_input_control_mod,  ONLY: problem_number
      
      IMPLICIT NONE

!  Subroutine  - check the idealised namelist.
!
! Description:
!   Read the namelist, assigning and freeing Fortran units as
!   required.
!
! Method:
!   The namelists is provided in one file
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

! Arguments
      Integer :: Len_IntHd
      Integer :: Len_RealHd
      Integer :: Len_LevDepC1
      Integer :: Len_LevDepC2

      Integer :: IntHd   (Len_IntHd)
      Real    :: RealHd  (Len_RealHd)
      Real    :: LevDepC (Len_LevDepC1, Len_LevDepC2)

! Local Variables/Paramters
      Character (Len=*), Parameter :: RoutineName =                     &
                                      'Check_Idealise_4A'
      Integer, Parameter           :: max_levels = 160
      Integer, Parameter           :: jtheta     = 1
      Integer, Parameter           :: jrho       = 2

      Character (Len=80)  :: Cmessage
      CHARACTER (LEN=80)  :: cmessage_invalid_flag
      CHARACTER (LEN=filenamelength)  :: filename
      Integer             :: ErrorStatus
      Integer             :: status
      Integer             :: model_levels
      Integer             :: i             ! looper
      Integer             :: j             ! looper
      Integer             :: nft
      Logical             :: l_exist

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Comdecks
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
! Description: COMDECK containing surface types
!  for use in idealised problems
!

      INTEGER, PARAMETER :: surface_zero=0
      INTEGER, PARAMETER :: surface_ellipse=1
      INTEGER, PARAMETER :: surface_ridge=2
      INTEGER, PARAMETER :: surface_plateau=3
      INTEGER, PARAMETER :: surface_massif=4
      INTEGER, PARAMETER :: surface_mask=5
      INTEGER, PARAMETER :: surface_gauss=6
      INTEGER, PARAMETER :: surface_ridge_series=7
! ENDGAME-only parameters
      INTEGER, PARAMETER :: surface_schar_ridge=8
      INTEGER, PARAMETER :: surface_baroclinic=9
! End of ENDGAME-only parameters
      INTEGER, PARAMETER :: surface_dump=10
! Description: COMDECK containing vertical grid types
!  for use in idealised problems
!
      INTEGER, PARAMETER :: vert_regular=1
      INTEGER, PARAMETER :: vert_quadratic_theta=21
      INTEGER, PARAMETER :: vert_bi_quadratic=22
      INTEGER, PARAMETER :: vert_quadratic_uvtheta=23
      INTEGER, PARAMETER :: vert_schar=3
      INTEGER, PARAMETER :: vert_dwd=4
      INTEGER, PARAMETER :: vert_stretch_plus_regular=5
      INTEGER, PARAMETER :: vert_quad_stretch_thin=6
      INTEGER, PARAMETER :: vert_regular_thin=7
      INTEGER, PARAMETER :: vert_geometric_theta=8
      INTEGER, PARAMETER :: vert_dump=10
      INTEGER, PARAMETER :: vert_idl_um_grid=11

      Namelist /Idealise/ L_idealised_data,                             &
             L_initialise_data, surface_type, grid_number, grid_flat,   &
             first_theta_height, thin_theta_height, height_domain,      &
             tprofile_number, qprofile_number, uvprofile_number,        &
             theta_surface, p_surface, Brunt_Vaisala,                   &
             dtheta_dz1, height_dz1,                                    &
             t_horizfn_number, t_horizfn_data,                          &
             uv_horizfn_number, u_in, v_in, height_u_in, q1,            &
             u_ramp_start, u_ramp_end, ujet_lat, ujet_width, r_plane,   &
             h_o, grow_steps, Witch_power,                              &
             lambda_fraction, phi_fraction,                             &
             half_width_x, half_width_y, plat_size_x, plat_size_y,      &
             first_constant_r_rho_level_new, L_constant_dz,             &
             L_polar_wind_zero, L_rotating,                             &
             L_trivial_trigs, f_plane, ff_plane, L_vert_Coriolis,       &
             L_fixed_lbcs, L_pressure_balance,                          &
             L_wind_balance, L_rotate_winds,                            &
             IdlSurfFluxSeaOption, IdlSurfFluxSeaParams,                &
             L_spec_z0, roughlen_z0m, roughlen_z0h,                     &
             big_layers, transit_layers, big_factor, vert_grid_ratio,   &
             SuHe_newtonian_timescale_ka, SuHe_newtonian_timescale_ks,  &
             SuHe_pole_equ_deltaT, SuHe_static_stab,                    &
             base_frictional_timescale, SuHe_sigma_cutoff,              &
             SuHe_relax, SuHe_fric,                                     &
             L_SH_williamson, Instability_diagnostics,                  &
             tropics_deg,                                               &
             num_profile_data,                                          &
             zprofile_data, tprofile_data, qprofile_data,               &
             num_uvprofile_data,                                        &
             z_uvprofile_data, uprofile_data, vprofile_data,            &
             tforce_option, qforce_option, uvforce_option,              &
             num_tforce_levels, num_tforce_times,                       &
             tforce_time_interval,                                      &
             num_qforce_levels, num_qforce_times,                       &
             qforce_time_interval,                                      &
             num_uvforce_levels, num_uvforce_times,                     &
             uvforce_time_interval,newtonian_timescale,                 &
             z_tforce_data, tforce_data,                                &
             z_qforce_data, qforce_data,                                &
             z_uvforce_data, uforce_data, vforce_data,                  &
             tforce_data_modlev, qforce_data_modlev,                    &
             uforce_data_modlev, vforce_data_modlev,                    &
             pforce_option,                                             &
             num_pforce_times, pforce_time_interval,                    &
             p_surface_data, L_pforce,                                  &
             L_perturb_t, perturb_magnitude_t,                          &
             L_perturb_q, perturb_magnitude_q,                          &
             L_perturb_correlate_tq,                                    &
             L_perturb_correlate_vert,                                  &
             L_perturb_correlate_time,                                  &
             perturb_type, perturb_height,                              &
             L_code_test, L_force, cool_rate,                           &
             L_cyclone, L_baroclinic,                                   &
             L_fix_orog_hgt_lbc, orog_hgt_lbc, L_force_lbc,             &
             idl_interp_option, zprofile_orog, hf,                      &
             idl_bubble_option, idl_bubble_max,                         &
             idl_bubble_width,  idl_bubble_depth, idl_bubble_height,    &
             idl_bubble_xoffset,idl_bubble_yoffset,                     &
             L_idl_bubble_saturate,                                     &
             L_damp, L_geo_for,  L_bomex,                               &
             DMPTIM, HDMP, ZDMP,                                        &
             u_geo, v_geo,                                              &
             delta_xi1, delta_xi2, base_xi1, base_xi2,                  &
             Trefer_number, T_surface,                                  &
             L_cartesian, eccentricity,                                 &
             tstep_plot_frequency, tstep_plot_start, L_rotate_grid,     &
             grid_NP_lon, grid_NP_lat, AA_jet_m, AA_jet_n, AA_jet_u0,   &
             AA_jet_A, L_baro_Perturbed, L_shallow, L_const_grav,       &
             L_HeldSuarez,L_HeldSuarez1_drag, L_HeldSuarez2_drag,       &
             L_baro_inst, L_deep_baro_inst, T0_P, T0_E, b_const,        &
             k_const, ring_height,                                      &
             theta_pert, L_isothermal, L_exact_profile, L_balanced,     &
             L_solid_body, angular_velocity,                            &
             Nxi1L, Nxi1V, Nxi2L, Nxi2V,                                &
             delta_xi1_H, delta_xi1_L, delta_xi2_H, delta_xi2_L         
! End of ENDGAME-only code

      IF (lhook) CALL dr_hook('CHECK_IDEALISE_4A',zhook_in,zhook_handle)
      model_levels = IntHd(8)

      nft=106      ! set unit number explicitly

! Header values used in subroutine
!     IntHd(13)    : Number of boundary layer levels
!     IntHd(24)    : First rho level at which height is constant


! ENDGAME-only code
      L_cartesian = .FALSE.
      Eccentricity = 0.0
      Trefer_number = 2
      T_surface = 280.0
      base_xi1 = 0.0
      base_xi2 = 0.0
      delta_xi1 = 0.0
      delta_xi2 = 0.0
! EG Variable resolution params      
       Nxi1L = 0
       Nxi1V = 0
       Nxi2L = 0
       Nxi2V = 0
       delta_xi1_H = 0.0
       delta_xi1_L = 0.0
       delta_xi2_H = 0.0
       delta_xi2_L = 0.0
! Parameters for solid-body rotation/jets
      L_rotate_grid=.FALSE.
      grid_NP_lon=0.0
      grid_NP_lat=90.0
      AA_jet_m=1
      AA_jet_n=1
      AA_jet_u0=0.0
      AA_jet_A=0.0
      L_baro_Perturbed = .FALSE.
      L_shallow = .FALSE.
      L_balanced = .FALSE.
      L_exact_profile = .FALSE.
      L_HeldSuarez = .FALSE.
      L_HeldSuarez1_drag = .TRUE.
      L_HeldSuarez2_drag = .FALSE.
      L_solid_body = .FALSE.
      L_baro_inst = .FALSE.
      L_const_grav = .TRUE.
      L_isothermal = .TRUE.
      ring_height = 5000.0
      theta_pert = 1.0
      angular_velocity = 7.292116E-5
      tstep_plot_frequency = 1
      tstep_plot_start = -1
      ! deep atmosphere baroclinic wave      
      L_deep_baro_inst = .FALSE.
      T0_P = 240.0
      T0_E = 310.0
      b_const = 2
      k_const = 3
! End of ENDGAME-only code

! Set model defaults
      L_initialise_data = .FALSE.
       tropics_deg = 30.0
      L_code_test = .FALSE.

! Set IDEALISE defaults
      surface_type = 10
      grid_number = 10
      first_theta_height = 10.0
      thin_theta_height = 1.0
      grid_flat = 3
      first_constant_r_rho_level_new = -1
      tprofile_number = 10
      qprofile_number = 10
      uvprofile_number = 10
      lambda_fraction = 0.5
      phi_fraction = 0.5
      h_o = 0.0
      grow_steps = 0
      plat_size_x = 0.0
      plat_size_y = 0.0
      Witch_power = 1.5
! Half-widths for Witch of Agnesi.  Also use to mask real orography
      half_width_x = 2500000.0
      half_width_y = 2500000.0
      theta_surface = 280.0
      p_surface = 100000.0
      height_domain = 10000.
      height_domain = 80000.
      Brunt_Vaisala =  0.01
      L_constant_dz = .TRUE.
      L_trivial_trigs = .FALSE.
      L_vert_Coriolis = .FALSE.
      r_plane = -90.0
      f_plane = -90.0
      ff_plane = -90.0
      L_rotating = .TRUE.
      L_fixed_lbcs = .FALSE.
      L_force_lbc = .FALSE.
      L_fix_orog_hgt_lbc = .FALSE.
      orog_hgt_lbc = 0.0
      idl_interp_option = 1  ! Interpolation -> const on height levels
      zprofile_orog = 0.0    ! Sea level
      hf = 0.0
      L_wind_balance = .FALSE.
      L_rotate_winds = .FALSE.
      L_pressure_balance = .FALSE.
      L_polar_wind_zero= .FALSE.
      L_perturb = .FALSE.
      perturb_factor = 1.0
      L_perturb_t         = .False.
      perturb_magnitude_t = 0.5      ! Kelvin
      L_perturb_q         = .False.
      perturb_magnitude_q = 0.5E-3   ! kg/kg
      L_perturb_correlate_tq   = .True.
      L_perturb_correlate_vert = .True.
      L_perturb_correlate_time = .True.
      perturb_type        = 1        ! random
      perturb_height      = 0.0
      L_force = .FALSE.
      L_cyclone = .FALSE.
      L_baroclinic = .FALSE.
      L_damp = .FALSE.
      L_geo_for = .FALSE.
      L_bomex = .FALSE.
      DMPTIM = 0.0
      HDMP = 0.0
      ZDMP = 0.0 
      u_geo = 0.0 
      v_geo = 0.0 

      cool_rate = 0.0
      q1 = 70.0  ! 70% relative humidity
      t_horizfn_number = 0
      Do i = 1, 10
        t_horizfn_data(i) = 0.0
      End Do
      uv_horizfn_number = 0
      Do i = 1, 4
        u_in(i) = 0.0
        v_in(i) = 0.0
      End Do
      u_ramp_start = 0.0
      u_ramp_end = -90.0
      ujet_lat = -90.0
      ujet_width = 0.0
      Do i = 1, 3
        height_u_in(i) = 0.0  ! dummy value - not used for constant u
! ENDGAME-only code
        dtheta_dz1(i) = 0.0
! End of ENDGAME-only code
        height_dz1(i) = 0.0
      End Do
      big_layers = 0
      transit_layers = 0
      big_factor = 1.0
      vert_grid_ratio = 1.0 ! Used by grid_number = vert_geometric_theta
      num_profile_data   = 0
      Do i = 1, max_num_profile_data
        zprofile_data(i) = 0.0
        tprofile_data(i) = 0.0
        qprofile_data(i) = 0.0
      End Do
      num_uvprofile_data = 0
      Do i = 1, max_num_profile_data
        z_uvprofile_data(i) = 0.0
        uprofile_data(i)    = 0.0
        vprofile_data(i)    = 0.0
      End Do
      pforce_option = 0
      num_pforce_times = 1
      pforce_time_interval = 600.0
      L_pforce = .FALSE.
      tforce_option = 0
      qforce_option = 0
      uvforce_option = 0
      num_tforce_levels = 1
      num_tforce_times = 1
      tforce_time_interval = 600.0
      num_qforce_levels = 1
      num_qforce_times = 1
      qforce_time_interval = 600.0
      num_uvforce_levels = 1
      num_uvforce_times = 1
      uvforce_time_interval = 600.0
      newtonian_timescale = 3600.0
      Do i = 1, max_num_profile_data
        z_tforce_data(i) = 0.0
        z_qforce_data(i) = 0.0
        z_uvforce_data(i) = 0.0
        Do j = 1, max_num_force_times
          tforce_data(i,j)=0.0
          qforce_data(i,j)=0.0
          uforce_data(i,j)=0.0
          vforce_data(i,j)=0.0
        End Do
      End Do
      Do j = 1, max_num_force_times
        Do i = 1, model_levels_max
          tforce_data_modlev(i,j)= 0.0
          qforce_data_modlev(i,j)= 0.0
          uforce_data_modlev(i,j)= 0.0
          vforce_data_modlev(i,j)= 0.0
        End Do
      End Do
      Do j = 1, max_num_force_times
        p_surface_data(j)=0.0
      End Do

      Do i=1,idl_max_num_bubbles
        idl_bubble_option(i) = 0      ! Default no bubbles
        idl_bubble_max(i)    = 1.0    ! 1 K
        idl_bubble_height(i) = 1000.  ! 1 km
        idl_bubble_xoffset(i)= 0.5    ! Centre of domain
        idl_bubble_yoffset(i)= 0.5    ! Centre of domain
        idl_bubble_width(i)  = 1000.  ! 1 km
        idl_bubble_depth(i)  = 1000.  ! 1 km
        L_idl_bubble_saturate(i) = .False.
      End Do
! Set dynamical core defaults
      L_SH_williamson = .FALSE.
      SuHe_pole_equ_deltaT = 60.
      SuHe_static_stab = 10.
      SuHe_sigma_cutoff = 0.7
      base_frictional_timescale = 1.1574e-5
      SuHe_newtonian_timescale_ka = 2.893e-7
      SuHe_newtonian_timescale_ks = 2.893e-6
      SuHe_relax = 2
      SuHe_fric = 2

      IdlSurfFluxSeaOption = 0
      IdlSurfFluxSeaParams(:) = 0.0

! defaults for specification of roughness length
      L_spec_z0 = .FALSE.
      roughlen_z0m = 0.0
      roughlen_z0h = 0.0

      if(Problem_number  /=  0) then

! Now find the namelist filename from Env Vars

      Call Fort_Get_Env( 'IDEALISE', 8, FileName, filenamelength,       &
           ErrorStatus )

      If ( ErrorStatus /= 0 ) Then
        ErrorStatus = 10
        Cmessage =                                                      &
     &  'Unable to Obtain Idealise Filename from Environment'
        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      End If

      FileName = Trim( FileName )

! Check file exists - Do all this by hand (not using file_open)
! as wish to do this on *all* PEs together

      Inquire( file=FileName, exist=l_exist )

      If ( .Not. l_exist ) Then
        Write (cmessage,fmt='(3A)') 'Idealise file: ',FileName,' does not exist!'
        ErrorStatus = 20
        Call Ereport( routineName, errorStatus, cmessage )
      End If

! Open the file containing Idealised model settings
      Open( Unit=nft, File=FileName, IOstat=status )

      If ( PrintStatus >= PrStatus_Oper ) Then
        Write (6,*) '-Idealise settings file: ',FileName
      End If

! Quick error check to make sure parameter max_levels isn't too small
      If ( max_levels < model_levels ) Then
        ErrorStatus = 10
        Cmessage = 'Internal paramter max_levels is too small!'
        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      End If

! Read Idealise Namelist
      READ( Unit = nft, Nml=Idealise, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist Idealise")
      

prst1: IF ( PrintStatus >= PrStatus_Normal ) THEN

        cmessage_invalid_flag='invalid idealised flag set'

!     Write out namelist variables for diagnostic
        WRITE ( UNIT = 6, FMT='(A)') 'Values in IDEALISE Namelist.'
        WRITE ( UNIT = 6, FMT='(A,L1)') 'L_initialise_data ',L_initialise_data

        IF ( L_code_test ) THEN
            errorstatus = 106
            CALL ereport( routinename, errorstatus, cmessage )
        END IF
       END IF prst1
        If( L_initialise_data )then
        
! ENDGAME-only code
          IF (L_baro_inst) THEN
           IF ( PrintStatus >= PrStatus_Normal ) THEN
             WRITE(6,fmt='(A)')
             WRITE(6,fmt='(A)') '*************************************'
             WRITE(6,fmt='(A)') '      Baroclinic Wave Test           '
             WRITE(6,fmt='(A)') '   Jablonowski & Williamson (2006)   '
             WRITE(6,fmt='(A)') '       QJRMS 132, 2943--2975         '
             WRITE(6,fmt='(A)') '*************************************'
             WRITE(6,fmt='(A)')
             IF (.NOT. L_shallow) THEN
               WRITE(6,fmt='(A)') ' WARNING: L_shallow <> .TRUE.'
               WRITE(6,fmt='(A)') ' ======='
             END IF
             IF (surface_type /= surface_baroclinic) THEN
               WRITE(6,fmt='(A,I1)') ' WARNING: surface_type <> ',     &
                                      surface_baroclinic
               WRITE(6,fmt='(A)') ' ======='
             END IF
             WRITE(6,fmt='(A)') '*************************************'
           END IF
         ! Channel baroclinic
           baro_phi0 = f_plane  
         ELSE IF ( L_deep_baro_inst ) THEN
! deep atmosphere baroclinic wave   
           IF ( PrintStatus >= PrStatus_Normal ) THEN
             WRITE(6,fmt='(A)')
             WRITE(6,fmt='(A)') '*************************************'
             WRITE(6,fmt='(A)') ' Deep Atmosphere Baroclinic Wave Test' 
             WRITE(6,fmt='(A)') '   Staniforth (2011)                 '
             WRITE(6,fmt='(A)') '*************************************'
             WRITE(6,fmt='(A,E16.8,A,E16.8)')                          &
                         ' T0_Pole = ',T0_P,', T0_Equator = ',T0_E   
             WRITE(6,fmt='(A,I5,A,I5)') ' k = ',k_const,', b = ',b_const  
             WRITE(6,fmt='(A)')
             IF (surface_type /= surface_zero) THEN
               WRITE(6,fmt='(A,I1)') ' WARNING: surface_type <> ',     &
                                      surface_zero
               WRITE(6,fmt='(A)') ' ======='
             END IF
             WRITE(6,fmt='(A)') '*************************************'
           END IF
         END IF
! End of ENDGAME-only code

prst2: IF ( PrintStatus >= PrStatus_Normal ) THEN
          WRITE ( UNIT = 6, fmt='(A,I2)') 'surface_type ',surface_type
          WRITE ( UNIT = 6, fmt='(A,I3)') 'grid_number ',grid_number
          SELECT CASE(grid_number)
          CASE(vert_stretch_plus_regular)
            WRITE ( UNIT = 6, fmt='(A,I5)') '..big_layers ',big_layers
            WRITE ( UNIT = 6, fmt='(A,I5)') '..transit_layers ',transit_layers
            WRITE ( UNIT = 6, fmt='(A,E16.8)') '..big_factor ',big_factor
          CASE(vert_quad_stretch_thin)
            WRITE ( UNIT = 6, fmt='(A,E16.8)') '..first_theta_height ',       &
                                        first_theta_height
            WRITE ( UNIT = 6, fmt='(A,E16.8)') '..thin_theta_height ',        &
                                        thin_theta_height
          CASE(vert_geometric_theta)
            WRITE ( UNIT = 6, fmt='(A,E16.8)') '..vert_grid_ratio ',          &
                                        vert_grid_ratio
          ENDSELECT
          WRITE ( UNIT = 6, fmt='(A,I1)') 'grid_flat ',grid_flat
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'height_domain ',height_domain
          WRITE ( UNIT = 6, fmt='(A,I5)') 'first_constant_r_rho_level_new ',  &
                                    first_constant_r_rho_level_new
          WRITE ( UNIT = 6, fmt='(A,I2)') 'tprofile_number ',tprofile_number
          WRITE ( UNIT = 6, fmt='(A,I2)') 'qprofile_number ',qprofile_number
          WRITE ( UNIT = 6, fmt='(A,I2)') 'uvprofile_number ',uvprofile_number
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'theta_surface ',theta_surface
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'p_surface ',p_surface
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'Brunt_Vaisala ',Brunt_Vaisala
          WRITE ( UNIT = 6, fmt='(A  )') 'dtheta_dz1'
          WRITE ( UNIT = 6, Fmt='(3F10.7)' )( dtheta_dz1(i),i=1,3 )
          WRITE ( UNIT = 6, fmt='(A  )') 'height_dz1'
          WRITE ( UNIT = 6, Fmt='(3F10.3)' )( height_dz1(i),i=1,3 )
          WRITE ( UNIT = 6, fmt='(A,I5)') 't_horizfn_number ',t_horizfn_number
          If (t_horizfn_number  /=  0) Then
            WRITE ( UNIT = 6, fmt='(A  )') 't_horizfn_data'
            WRITE ( UNIT = 6, Fmt='(10F10.7)' )(t_horizfn_data(i),i=1,10)
          End If
          WRITE ( UNIT = 6, fmt='(A,I5)') 'uv_horizfn_number ',uv_horizfn_number
          WRITE ( UNIT = 6, fmt='(A  )') 'height_u_in'
          WRITE ( UNIT = 6, Fmt='(3F10.3)' )( height_u_in(i),i=1,3 )
          WRITE ( UNIT = 6, fmt='(A  )') 'u_in'
          WRITE ( UNIT = 6, Fmt='(4F10.3)' )( u_in(i),i=1,4 )
          WRITE ( UNIT = 6, fmt='(A  )') 'v_in'
          WRITE ( UNIT = 6, Fmt='(4F10.3)' )( v_in(i),i=1,4 )
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'q1 ',q1
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'orog_height h_o ',h_o
          WRITE ( UNIT = 6, fmt='(A,I5)') 'grow_steps ',grow_steps
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'Witch_power ',Witch_power
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'lambda_fraction ',lambda_fraction
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'phi_fraction ',phi_fraction
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'half_width_x ',half_width_x
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'half_width_y ',half_width_y
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'plat_size_x ',plat_size_x
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'plat_size_y ',plat_size_y
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_constant_dz ',L_constant_dz
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'u_ramp_start ',u_ramp_start
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'u_ramp_end ',u_ramp_end
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'ujet_lat ',ujet_lat
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'ujet_width ',ujet_width
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'r_plane ',r_plane
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'f_plane ',f_plane
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'ff_plane ',ff_plane
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_trivial_trigs ',L_trivial_trigs
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_rotating ',L_rotating
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_vert_Coriolis ',L_vert_Coriolis
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_fixed_lbcs ',L_fixed_lbcs
          WRITE (UNIT = 6, fmt='(A,L1)') 'L_force_lbc ', L_force_lbc
          WRITE (UNIT=6, fmt='(A,L1)') 'L_fix_orog_hgt_lbc ',L_fix_orog_hgt_lbc
          WRITE (UNIT = 6, fmt='(A,E16.8)') 'orog_hgt_lbc ', orog_hgt_lbc
          WRITE (UNIT = 6, fmt='(A,I3)') 'idl_interp_option ',idl_interp_option
          WRITE (UNIT = 6, fmt='(A,E16.8)') 'zprofile_orog ', zprofile_orog
          WRITE (UNIT = 6, fmt='(A,E16.8)') 'hf ', hf
          WRITE (UNIT = 6, fmt='(A,E16.8)') 'p_surface_data ', p_surface_data
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_wind_balance ',L_wind_balance
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_rotate_winds ',L_rotate_winds
          WRITE ( UNIT = 6, fmt='(A,L1)')'L_pressure_balance ',   &
                                          L_pressure_balance
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_polar_wind_zero ',L_polar_wind_zero
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_perturb ',L_perturb
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'tropics_deg ',tropics_deg

! ENDGAME-only code
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'base_xi1 ',base_xi1
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'base_xi2 ',base_xi2
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'delta_xi1 ',delta_xi1
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'delta_xi2 ',delta_xi2
          WRITE ( UNIT = 6, fmt='(A,I5)') 'Trefer_number ',Trefer_number
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 't_surface ',t_surface
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_shallow ', L_shallow
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_balanced ', L_balanced
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_const_grav ', L_const_grav
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_HeldSuarez ', L_HeldSuarez
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_baro_inst ', L_baro_inst
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_deep_baro_inst ', L_deep_baro_inst
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'T0_P = ',T0_P
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'T0_E = ',T0_E
          WRITE ( UNIT = 6, fmt='(A,I5)') 'k_const = ',k_const
          WRITE ( UNIT = 6, fmt='(A,I5)') 'b_const = ',b_const
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_solid_body ', L_solid_body
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'Ring height = ', ring_height     
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'theta_pert = ', theta_pert
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_isothermal ', L_isothermal
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'angular_velocity ',           &
                                              angular_velocity
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_exact_profile ', L_exact_profile
          If (L_rotate_grid) Then
            WRITE ( UNIT = 6, fmt='(A,E16.8)') 'grid_NP_lon ', grid_NP_lon
            WRITE ( UNIT = 6, fmt='(A,E16.8)') 'grid_NP_lat ', grid_NP_lat
          Endif
! Following should be under a L_AA_jet switch
          WRITE ( UNIT = 6, fmt='(A,I5)') 'AA_jet_m ', AA_jet_m
          WRITE ( UNIT = 6, fmt='(A,I5)') 'AA_jet_n ', AA_jet_n
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'AA_jet_A ', AA_jet_A
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'AA_jet_u0 ', AA_jet_u0
! End of ENDGAME-only code
          If (L_perturb) Then
            WRITE ( UNIT = 6, fmt='(A,E16.8)') 'perturb_factor ',perturb_factor
            WRITE ( UNIT = 6, fmt='(A,E16.8)') 'perturb_height ',perturb_height
          EndIf
END IF prst2

        If (tprofile_number == tp_namelist                              &
           .or. qprofile_number == qp_namelist                         &
           .or. qprofile_number == qp_namelist_rh) Then

           If ( PrintStatus >= PrStatus_Normal ) &
             WRITE (UNIT=6, Fmt='(A,E16.8)') 'num_profile_data ',      &
                                              num_profile_data

          ! Check to make sure data points in profile is less than max
          If ((num_profile_data == 0) .or.                              &
             (num_profile_data > max_num_profile_data)) Then
            WRITE (6,'(A,I5)') 'max_num_profile_data ',max_num_profile_data
            WRITE(Cmessage,'(A)')                                           &
             'Idealised namelist vertical profile data:'               &
             //'Zero or too many points. '                             &
             //'num_profile_data must be 0 <= max_num_profile_data'
            ErrorStatus = 1
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If

          IF ( PrintStatus >= PrStatus_Normal ) THEN
            WRITE (UNIT=6, Fmt='(A)') 'zprofile_data'
            WRITE (UNIT=6, Fmt='(F10.3)')                               &
                   (zprofile_data(i),i=1,num_profile_data)
          END IF

        End If



        If (tprofile_number == tp_namelist .AND.                        &
            PrintStatus >= PrStatus_Normal ) Then
          WRITE (UNIT=6, Fmt='(A)') 'tprofile_data'
          WRITE (UNIT=6, Fmt='(F10.3)')                                 &
                   (tprofile_data(i),i=1,num_profile_data)
        End If
        If ((qprofile_number == qp_namelist                             &
           .or. qprofile_number == qp_namelist_rh).AND.                &
            PrintStatus >= PrStatus_Normal ) Then
          WRITE (UNIT=6, Fmt='(A)') 'qprofile_data'
          WRITE (UNIT=6, Fmt='(F10.7)')                                 &
                   (qprofile_data(i),i=1,num_profile_data)
        End If


        If (uvprofile_number == uv_vert_namelist) Then

          ! If num_uvprofile_data not set in namelist then assume
          ! data is on same levels as t,q data
          If (num_uvprofile_data == 0 .AND. num_profile_data /= 0)      &
           THEN
           IF ( PrintStatus >= PrStatus_Normal )                        &
            WRITE (UNIT=6, Fmt='(A)') 'Assuming uv data is on the same'     &
                                //' height levels as tq data.'
            num_uvprofile_data = num_profile_data
            Do i = 1, max_num_profile_data
              z_uvprofile_data(i) = zprofile_data(i)
            End Do
          End If

          IF ( PrintStatus >= PrStatus_Normal )                        &
            WRITE (UNIT=6, Fmt='(A,E16.8)') 'num_uvprofile_data ',      &
                                         num_uvprofile_data

          ! Check to make sure no. points in profile less than max
          IF ((num_uvprofile_data == 0) .OR.                            &
              (num_uvprofile_data > max_num_profile_data)) THEN
            WRITE (6,fmt='(A,E16.8)') 'max_num_profile_data ',          &
                                       max_num_profile_data
            WRITE(Cmessage,fmt='(A)')                                           &
              'Idealised namelist vertical profile data:'               &
              //'Zero or too many points. '                             &
              //'num_uv_profile_data must be 0 <= max_num_profile_data'
            ErrorStatus = 1
            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If

          ! WRITE out arrays
          IF ( PrintStatus >= PrStatus_Normal ) THEN
            WRITE (UNIT=6, Fmt='(A)') 'z_uvprofile_data'
            WRITE (UNIT=6, Fmt='(F10.3)')                               &
                    (z_uvprofile_data(i),i=1,num_uvprofile_data)
            WRITE (UNIT=6, Fmt='(A)') 'uprofile_data'
            WRITE (UNIT=6, Fmt='(F10.4)')                               &
                    (uprofile_data(i),i=1,num_uvprofile_data)
            WRITE (UNIT=6, Fmt='(A)') 'vprofile_data'
            WRITE (UNIT=6, Fmt='(F10.4)')                               &
                    (vprofile_data(i),i=1,num_uvprofile_data)
          END IF

        End If
        IF ( PrintStatus >= PrStatus_Normal ) THEN
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_force =',L_force
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_cyclone =',L_cyclone
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_baroclinic =',L_baroclinic
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_damp =',L_damp
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_geo_for =',L_geo_for
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_bomex =',L_bomex
          WRITE ( UNIT = 6, fmt='(A)') 'Damping layer settings'
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'DMPTIM =',DMPTIM
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'HDMP =',HDMP
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'ZDMP =',ZDMP
          WRITE ( UNIT = 6, fmt='(A)') 'Geostrophic forcings'
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'u_geo =',u_geo
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'v_geo =',v_geo
        END IF
        ! -------------------------------------------------
        ! Theta forcing
        ! -------------------------------------------------

        If (tforce_option == 1.AND.PrintStatus >= PrStatus_Normal) Then
          WRITE (UNIT = 6, fmt='(A)')                                       &
                'Forcing increments added to theta field'
        ELSE IF (tforce_option == 2.AND.PrintStatus >= PrStatus_Normal) &
          Then
          WRITE ( UNIT = 6, fmt='(A)')                                      &
                'Relaxation forcing for theta field'
! Option not yet implemented
!          Elseif (tforce_option  ==  3) Then
!            WRITE (UNIT = 6, fmt='(A,E16.8)')
!     &            'Theta reset to specified forcing data'
        Else
           IF ( PrintStatus >= PrStatus_Normal ) WRITE(UNIT=6,fmt='(A)')    &
           'No forcing of theta'
        Endif

        ! -------------------------------------------------
        ! Humidity forcing
        ! -------------------------------------------------

        IF ( PrintStatus >= PrStatus_Normal ) THEN
          If (qforce_option  ==  1) Then
            WRITE ( UNIT = 6, fmt='(A)')                                    &
                'Forcing increments added to q field'
          ELSE IF (qforce_option  ==  2) Then
            WRITE ( UNIT = 6, fmt='(A)')                                    &
                'Relaxation forcing for q field'
! Option not yet implemented
!          Elseif (qforce_option  ==  3) Then
!            WRITE ( UNIT = 6, fmt='(A,E16.8)')
!     &            'q reset to specified forcing data'
          Else
            WRITE ( UNIT = 6, fmt='(A)') 'No forcing of q'
          Endif
        END IF
        ! -------------------------------------------------
        ! Horizontal wind forcing
        ! -------------------------------------------------

        IF ( PrintStatus >= PrStatus_Normal ) THEN
          If (uvforce_option  ==  1) Then
            WRITE ( UNIT = 6, fmt='(A)')                                    &
                'Forcing increments added to u and v fields'
          ELSE IF (uvforce_option  ==  2) Then
            WRITE ( UNIT = 6, fmt='(A)')                                    &
                'Relaxation forcing for u and v fields'
! Option not yet implemented
!          Elseif (uvforce_option  ==  3) Then
!            WRITE ( UNIT = 6, fmt='(A,E16.8)')
!     &            'u and v reset to specified forcing data'
          Else
            WRITE ( UNIT = 6, fmt='(A)')                                    &
               'No forcing of winds'
          Endif
        END IF

      ELSE IF( Problem_number  ==  dynamical_core)then

        IF ( PrintStatus >= PrStatus_Normal ) THEN
          WRITE ( UNIT = 6, fmt='(A,L1)') 'L_SH_williamson ',L_SH_williamson
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'SuHe_newtonian_timescale_ka ',  &
                               SuHe_newtonian_timescale_ka
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'SuHe_newtonian_timescale_ks ',  &
                               SuHe_newtonian_timescale_ks
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'SuHe_pole_equ_deltaT ',         &
                               SuHe_pole_equ_deltaT
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'SuHe_static_stab ',             &
                               SuHe_static_stab
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'base_frictional_timescale ',    &
                               base_frictional_timescale
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'SuHe_sigma_cutoff ',            &
                               SuHe_sigma_cutoff
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'SuHe_relax ',SuHe_relax
          WRITE ( UNIT = 6, fmt='(A,E16.8)') 'SuHe_fric ',SuHe_fric
        END IF

      EndIf ! L_initialise_data

        IF ( PrintStatus >= PrStatus_Normal ) THEN
          WRITE ( UNIT = 6,fmt='(A,I5)') 'Instability_diagnostics ',          &
                                Instability_diagnostics
        END IF


    Close( UNIT=nft )

    endif

    IF (lhook) CALL dr_hook('CHECK_IDEALISE_4A',zhook_out,zhook_handle)
    RETURN
    END SUBROUTINE CHECK_IDEALISE_4A
