! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Initialises trigonometric and Coriolis dynamical arrays (ENDGAME version)
!
! Subroutine Interface:
      SUBROUTINE set_trigs_4A(                                          &
                           row_length, rows, n_rows,model_levels,       &
                           base_lambda, base_phi,                       &
                           delta_lambda, delta_phi,                     &
                           base_lambda_wk, base_phi_wk,                 &
                           delta_lambda_wk, delta_phi_wk,               &
                           f_plane_rad, ff_plane_rad,                   &
                           two_Omega,                                   &
                           lat_rot_NP, long_rot_NP,                     &
                           lat_rot_NP_deg, long_rot_NP_deg,             &
                           lat_n, lat_s, long_e, long_w,                &
                           long_e_model, long_w_model,                  &
                           model_domain, rmdi) 

      USE trignometric_mod
      USE rimtypes
      USE dyn_coriolis_mod
      USE dyn_var_res_Mod,    ONLY: glambda_p, phi_p
      USE rot_coeff_mod,      ONLY: rot_coeff1, rot_coeff2 
      USE conversions_mod,    ONLY: pi_over_180

      USE swapable_field_mod, ONLY : swapable_field_pointer_type

      USE parkind1, ONLY: jpim, jprb       !DrHook
      USE yomhook,  ONLY: lhook, dr_hook   !DrHook
      USE eg_alpha_mod
      USE eg_alpha_ramp_mod
      USE horiz_grid_mod
      USE ref_pro_mod ! needed due to the inclusion of cruntimc
 
      USE Control_Max_Sizes
      USE UM_ParVars
      USE domain_params
      USE ereport_mod, ONLY : ereport
      USE integrity_mod
      USE eg_alpha_mod
      USE eg_alpha_ramp_mod
      USE horiz_grid_mod

      USE ref_pro_mod ! needed due to the inclusion of cruntimc
      USE eg_parameters_mod  ! needed due to the inclusion of cruntimc

      USE atm_fields_bounds_mod
      USE coriolis_mod
      
      USE dynamics_input_mod, ONLY:L_regular

      USE eqtoll_mod, ONLY: eqtoll
      USE w_coeff_mod, ONLY: w_coeff
      IMPLICIT NONE

!
! Description:
!   Initialises trigonometric and Coriolis dynamic arrays
!   
! Method:
!   Sets values from grid information
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

! Input variables.
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      INTEGER  ::  row_length
      INTEGER  ::  rows
      INTEGER  ::  n_rows
      INTEGER  ::  model_levels
      INTEGER  ::  model_domain


      REAL  :: base_lambda
      REAL  :: base_phi
      REAL  :: delta_lambda
      REAL  :: delta_phi
      REAL  :: f_plane_rad
      REAL  :: ff_plane_rad
      REAL  :: lat_rot_NP_deg
      REAL  :: long_rot_NP_deg
      REAL  :: lat_rot_NP
      REAL  :: long_rot_NP
      REAL  :: base_lambda_wk ! Longitude of first theta point in degs
      REAL  :: base_phi_wk    ! Latitude of first theta point in degs
      REAL  :: delta_lambda_wk ! EW (x) grid spacing in degrees
      REAL  :: delta_phi_wk    ! NS (y) grid spacing in degrees
      REAL  :: rmdi            ! For unset (missing) data
      REAL  :: two_Omega

      REAL :: lat_n
      REAL :: lat_s
      REAL :: long_e
      REAL :: long_w
      REAL :: long_w_model
      REAL :: long_e_model

! Local variables.
      INTEGER  ::  i, j
      INTEGER  ::  gi, gj
      INTEGER  ::  j0, j1

      REAL  ::  f1_temp        ! f1 term for f-plane or ff-plane
      REAL  ::  f2_temp        ! f2 term for f-plane or ff-plane
      REAL  ::  f3_temp        ! f3 term for f-plane or ff-plane
      REAL  ::  temp1 ! Used in calculation of Coriolis terms (LAM only)
      REAL  ::  temp2 ! Used in calculation of Coriolis terms (LAM only)

      REAL  ::  latitude_rot(row_length, n_rows)
      REAL  ::  longitude_rot(row_length, n_rows)
      REAL  ::  true_latitude_b(row_length, n_rows)
                                             ! true latitude on B-grid
      REAL  ::  true_longitude_b(row_length,n_rows)
                                             ! true longitude on B-grid

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


      INTEGER :: i_field

      TYPE(swapable_field_pointer_type) :: fields_to_swap(10)

! -----------------------------------------
! 1. set trig fields
! -----------------------------------------

    IF (lhook) CALL dr_hook('SET_TRIGS_4A',zhook_in,zhook_handle)

! Quick fix to get cartesian LAMS working
    IF ( L_cartesian ) L_trivial_trigs = .TRUE.

!!! Allocate latitude arrays for trignometric_mod module
      IF (.NOT. ALLOCATED(cos_theta_latitude)) THEN
        ALLOCATE (cos_theta_latitude (tdims_s%i_start:tdims_s%i_end,      &
                                      tdims_s%j_start:tdims_s%j_end))
      END IF
      IF (.NOT. ALLOCATED(sec_theta_latitude)) THEN
        ALLOCATE (sec_theta_latitude (tdims_s%i_start:tdims_s%i_end,      &
                                      tdims_s%j_start:tdims_s%j_end))
      END IF
      IF (.NOT. ALLOCATED(FV_cos_theta_latitude)) THEN
        ALLOCATE (FV_cos_theta_latitude (tdims_s%i_start:tdims_s%i_end,   &
                                         tdims_s%j_start:tdims_s%j_end))
      END IF
      IF (.NOT. ALLOCATED(FV_sec_theta_latitude)) THEN
        ALLOCATE (FV_sec_theta_latitude (tdims_s%i_start:tdims_s%i_end,   &
                                         tdims_s%j_start:tdims_s%j_end))
      END IF
      IF (.NOT. ALLOCATED(sin_theta_latitude)) THEN
        ALLOCATE (sin_theta_latitude (row_length, rows))
      END IF
      IF (.NOT. ALLOCATED(tan_theta_latitude)) THEN
        ALLOCATE (tan_theta_latitude (row_length, rows))
      END IF
      IF (.NOT. ALLOCATED(sin_v_latitude)) THEN
        ALLOCATE (sin_v_latitude (vdims%i_start:vdims%i_end,              &
                                  vdims%j_start:vdims%j_end))
      END IF
      IF (.NOT. ALLOCATED(tan_v_latitude)) THEN
        ALLOCATE (tan_v_latitude (vdims%i_start:vdims%i_end,              &
                                  vdims%j_start:vdims%j_end))
      END IF
      IF (.NOT. ALLOCATED(cos_v_latitude)) THEN
        ALLOCATE (cos_v_latitude (vdims_s%i_start:vdims_s%i_end,          &
                                  vdims_s%j_start:vdims_s%j_end))
      END IF
      IF (.NOT. ALLOCATED(sec_v_latitude)) THEN
        ALLOCATE (sec_v_latitude (vdims_s%i_start:vdims_s%i_end,          &
                                  vdims_s%j_start:vdims_s%j_end))
      END IF

!!! Allocate longitude and 'true' arrays for trignometric_mod module
      IF (.NOT. ALLOCATED(cos_theta_longitude)) THEN
        ALLOCATE (cos_theta_longitude (row_length, rows))
      END IF
      IF (.NOT. ALLOCATED(sin_theta_longitude)) THEN
        ALLOCATE (sin_theta_longitude (row_length, rows))
      END IF
      IF (.NOT. ALLOCATED(cos_u_longitude)) THEN
        ALLOCATE (cos_u_longitude (udims%i_start:udims%i_end,             &
                                   udims%j_start:udims%j_end))
      END IF
      IF (.NOT. ALLOCATED(sin_u_longitude)) THEN
        ALLOCATE (sin_u_longitude (udims%i_start:udims%i_end,             &
                                   udims%j_start:udims%j_end))
      END IF


      IF (.NOT.l_cartesian) THEN
        DO j = tdims%j_start,tdims%j_end
          DO i = tdims%i_start,tdims%i_end
            cos_theta_latitude(i,j) = COS(xi2_p(j))
            sin_theta_latitude(i,j) = SIN(xi2_p(j))
            FV_cos_theta_latitude(i,j) = cos_theta_latitude(i,j)
          END DO
        END DO
      ELSE
        cos_theta_latitude    = 1.0
        sin_theta_latitude    = 0.0
        FV_cos_theta_latitude = 1.0
      END IF

 

      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
          sec_theta_latitude(i, j) = 1. / cos_theta_latitude(i, j)
          tan_theta_latitude(i, j) = sin_theta_latitude(i,j) /            &
                                       cos_theta_latitude(i,j)
        END DO
      END DO
 

      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
          FV_sec_theta_latitude(i, j) = 1./FV_cos_theta_latitude(i,j)
        END DO
      END DO
  
      IF (.NOT.l_cartesian) THEN

        DO j = vdims%j_start,vdims%j_end
          DO i = vdims%i_start,vdims%i_end
            sin_v_latitude(i,j) = SIN(xi2_v(j))
            cos_v_latitude(i,j) = COS(xi2_v(j))
            sec_v_latitude(i,j) = 1./cos_v_latitude(i,j)
            tan_v_latitude(i,j) = sin_v_latitude(i,j) /                   &
                                    cos_v_latitude(i,j)
          END DO
        END DO

        DO j = tdims%j_start,tdims%j_end
          DO i = tdims%i_start,tdims%i_end
            sin_theta_longitude(i,j) = SIN(xi1_p(i))
            cos_theta_longitude(i,j) = COS(xi1_p(i))
           END DO
        END DO   
                
        DO j = udims%j_start,udims%j_end
          DO i = udims%i_start,udims%i_end           
            sin_u_longitude(i,j)     = SIN(xi1_u(i))
            cos_u_longitude(i,j)     = COS(xi1_u(i))
          END DO
        END DO

      ELSE

        DO j = vdims%j_start,vdims%j_end
          DO i = vdims%i_start,vdims%i_end
            sin_v_latitude(i,j) = 0.0
            cos_v_latitude(i,j) = 1.0
            sec_v_latitude(i,j) = 1.0
            tan_v_latitude(i,j) = 0.0
          END DO
        END DO

        DO j = tdims%j_start,tdims%j_end
          DO i = tdims%i_start,tdims%i_end
            sin_theta_longitude(i,j) = 0.0
            cos_theta_longitude(i,j) = 1.0
           END DO
        END DO   
                
        DO j = udims%j_start,udims%j_end
          DO i = udims%i_start,udims%i_end                               
            sin_u_longitude(i,j)     = 0.0
            cos_u_longitude(i,j)     = 1.0
          END DO
        END DO

      END IF

! -----------------------------------------
! 2. set Coriolis components
! -----------------------------------------

! ENDGAME/VATPOLES-only
!----------------------------------------------------------------------
!  Initialize fi_star coriolis terms: fi_star(k) = fi where fi is the
!  actual Coriolis term. Will be converted to fi_star after the
!  computing the metric terms.
!----------------------------------------------------------------------

      IF ( model_domain  ==  mt_global ) THEN

        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end

            IF ( L_shallow ) THEN
              f1_comp(i,j) = 0.0
              f2_comp(i,j) = 0.0 !cos_theta_latitude(i,j)
              f3_comp(i,j) = two_omega*SIN(xi2_p(j)) !sin_theta_latitude(i,j)
            ELSE
              f1_comp(i,j) = 0.0
              f2_comp(i,j) = two_omega*COS(xi2_p(j)) !cos_theta_latitude(i,j)
              f3_comp(i,j) = two_omega*SIN(xi2_p(j)) !sin_theta_latitude(i,j)
            END IF

          END DO
        END DO

      ELSE

! limited area model
        temp1 = cos( lat_rot_NP)
        temp2 = sin( lat_rot_NP)

        IF ( L_Cartesian .or. f_plane  >   -89.0) THEN

          f_plane_rad = f_plane * Pi_over_180
          ff_plane_rad = ff_plane * Pi_over_180
          f1_temp = - two_omega * temp1 * SIN(ff_plane_rad)
          f2_temp = two_omega * (COS(f_plane_rad) * temp2                 &
                     - SIN(f_plane_rad) * temp1 * COS(ff_plane_rad))
          f3_temp = two_omega * ( SIN(f_plane_rad) * temp2                &
                     + COS(f_plane_rad) * temp1 * COS(ff_plane_rad))
          IF ( L_vert_Coriolis) THEN
            f1_temp = 0.0
            f2_temp = 0.0
            f3_star = two_omega * ( SIN(f_plane_rad) * temp2              &
                    + COS(f_plane_rad) * temp1 * COS(ff_plane_rad))
          END IF  !  L_vert_Coriolis

          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              f1_comp(i,j) = f1_temp
              f2_comp(i,j) = f2_temp
              f3_comp(i,j) = f3_temp
            END DO
          END DO

        ELSE      !  L_Cartesian = .false.

          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              f1_comp(i,j) = - two_omega * temp1 *                        &
                               sin_theta_longitude(i,j)
              f2_comp(i,j) = two_omega*(cos_theta_latitude(i,j)*          &
                                 temp2 -sin_theta_latitude(i,j)*temp1     &
                                        *cos_theta_longitude(i,j) )
              f3_comp(i,j) = two_omega*(sin_theta_latitude(i,j)*          &
                                 temp2 +cos_theta_latitude(i,j)*temp1     &
                                        *cos_theta_longitude(i,j) )
            END DO
          END DO

        END IF !  L_Cartesian

      END IF ! If ( model_domain  ==  mt_global )
! End of ENDGAME/VATPOLES-only

!!! Allocate arrays for dyn_coriolis_mod module
      IF (.NOT. ALLOCATED(f1_at_v)) THEN
        ALLOCATE (f1_at_v (vdims_s%i_start:vdims_s%i_end,                 &
                           vdims_s%j_start:vdims_s%j_end ))
      END IF
      IF (.NOT. ALLOCATED(f2_at_u)) THEN
        ALLOCATE (f2_at_u (udims_s%i_start:udims_s%i_end,                 &
                           udims_s%j_start:udims_s%j_end ))
      END IF
      IF (.NOT. ALLOCATED(f3_at_u)) THEN
        ALLOCATE (f3_at_u (udims_s%i_start:udims_s%i_end,                 &
                           udims_s%j_start:udims_s%j_end))
      END IF
      IF (.NOT. ALLOCATED(f3_at_v)) THEN
        ALLOCATE (f3_at_v (vdims_s%i_start:vdims_s%i_end,                 &
                           vdims_s%j_start:vdims_s%j_end))
      END IF
      IF (.NOT. ALLOCATED(true_latitude)) THEN
        ALLOCATE (true_latitude (row_length, rows))
      END IF
      IF (.NOT. ALLOCATED(true_longitude)) THEN
        ALLOCATE (true_longitude (row_length, rows))
      END IF

      IF (model_domain ==  mt_global) THEN

        DO j = vdims%j_start,vdims%j_end
          DO i = vdims%i_start,vdims%i_end
            f1_at_v(i,j) = 0.
            f3_at_v(i,j) = two_omega * sin_v_latitude(i,j)
          END DO
        END DO
        DO j = udims%j_start,udims%j_end
          DO i = udims%i_start,udims%i_end
            f2_at_u(i,j) = two_omega * cos_theta_latitude(i+1,j)             
            f3_at_u(i,j) = two_omega * sin_theta_latitude(i+1,j)
          END DO
        END DO
        IF (L_regular) THEN
          DO j = 1, rows
            DO i = 1, row_length
              gi = datastart(1) + i - 1
              true_longitude(i,j) = (gi-1)*delta_lambda
            END DO
          END DO
        ELSE !  variable resolution
          DO j = 1, rows
            DO i = 1, row_length
              gi = datastart(1) + i - 1
              true_longitude(i,j) = glambda_p(gi) - base_lambda
            END DO
          END DO
        END IF !  L_regular

        true_latitude = ASIN(sin_theta_latitude)

! set polar values of true longitude to be all the same

        IF (at_extremity(PNorth)) THEN
          DO i = 1, row_length
            true_longitude(i,rows) = 0.
          END DO
        END IF

        IF (at_extremity(PSouth)) THEN
          DO i = 1, row_length
            true_longitude(i,1) = 0.
          END DO
        END IF

! get parameters for common latlonmax within mppac
! these need to be in degrees
      IF (L_regular) THEN
        lat_n = Base_phi_wk + (datastart(2)+rows-2) * Delta_phi_wk
        long_e = Base_lambda_wk +                                         &
                    (datastart(1)+row_length-2) * Delta_lambda_wk
      ELSE !  variable resolution
        lat_n =  phi_p(1,rows)  / Pi_over_180  
        gi = datastart(1) + row_length - 1
        long_e = glambda_p(gi) / Pi_over_180 
      END IF !  L_regular

      lat_s = Base_phi_wk
      long_w       = Base_lambda_wk
      long_w_model = long_w
      long_e_model = long_e

      IF(long_w >  180.0) long_w = long_w-360.0
      IF(long_e >  180.0) long_e = long_e-360.0

      ELSE     ! limited area model

        temp1 = COS(lat_rot_NP)
        temp2 = SIN(lat_rot_NP)

        IF( L_trivial_trigs .OR. f_plane  >   -89.0) THEN
          f1_temp = - two_omega * temp1 * SIN(ff_plane_rad)
          f2_temp = two_omega * (COS(f_plane_rad) * temp2                 &
                       - SIN(f_plane_rad) * temp1 * COS(ff_plane_rad))
          f3_temp = two_omega * ( SIN(f_plane_rad) * temp2                &
                     + COS(f_plane_rad) * temp1 * COS(ff_plane_rad))
          IF (L_vert_Coriolis) THEN
            f1_temp = 0.0
            f2_temp = 0.0
            f3_temp = two_omega * ( SIN(f_plane_rad) * temp2              &
                       + COS(f_plane_rad) * temp1 * COS(ff_plane_rad))
          END IF  !  L_vert_Coriolis

          DO j = vdims%j_start,vdims%j_end
            DO i = vdims%i_start,vdims%i_end
              f1_at_v(i,j) = f1_temp
              f3_at_v(i,j) = f3_temp
            END DO
          END DO
          DO j = udims%j_start,udims%j_end
            DO i = udims%i_start,udims%i_end
              f2_at_u(i,j) = f2_temp
              f3_at_u(i,j) = f3_temp
            END DO
          END DO

        ELSE    !  L_trivial_trigs = .false.

          DO j = vdims%j_start,vdims%j_end
            DO i = vdims%i_start,vdims%i_end
              f1_at_v(i,j) = - two_omega * temp1 *                        &
                               sin_theta_longitude(i,MIN(j+1,rows))             
              f3_at_v(i,j) = two_omega * ( sin_v_latitude(i,j) * temp2    &
                                        +cos_v_latitude(i,j) * temp1      &
                                        *cos_theta_longitude(i,MIN(j+1,rows))) 
            END DO
          END DO

          DO j = udims%j_start,udims%j_end
            DO i = udims%i_start,udims%i_end
              f2_at_u(i,j) = two_omega * (                                &
                                   cos_theta_latitude(i+1,j) * temp2 -    &       
                                   sin_theta_latitude(i+1,j) * temp1 *    &
                                      cos_u_longitude(i,j) )
              f3_at_u(i,j) = two_omega * (                                &
                                   sin_theta_latitude(i+1,j) * temp2 +    &       
                                   cos_theta_latitude(i+1,j) * temp1 *    &
                                      cos_u_longitude(i,j) )
            END DO
          END DO

        END IF   !  L_trivial_trigs

! calculate true longitude in radians

          DO j = 1, rows
            gj = datastart(2) + j - 1
            DO i = 1, row_length
              gi = datastart(1) + i - 1
              longitude_rot(i,j) = (base_lambda + (gi-1) * delta_lambda)  &
                                   / Pi_over_180
              latitude_rot(i,j) = (base_phi + (gj-1) * delta_phi)         &
                                   / Pi_over_180
            END DO
          END DO

        CALL eqtoll(latitude_rot, longitude_rot,                          &
                    true_latitude, true_longitude,                        &
                    lat_rot_NP_deg, long_rot_NP_deg, rows*row_length)

        DO j = 1, rows
          DO i = 1, row_length
            true_longitude(i,j) = true_longitude(i,j) * Pi_over_180
            true_latitude(i,j) = true_latitude(i,j) * Pi_over_180
          END DO
        END DO

! get parameters for common latlonmax within mppac
! these need to be in degrees
        lat_n = latitude_rot(1,rows)
        lat_s = latitude_rot(1,1)

        long_w       = longitude_rot(1,1)
        long_w_model = long_w
        long_e       = longitude_rot(row_length,1)
        long_e_model = long_e

        IF(long_w >  180.0) long_w = long_w-360.0
        IF(long_e >  180.0) long_e = long_e-360.0
        IF(long_w_model >= 360.0)long_w_model = long_w_model-360.
        IF(long_e_model >= 360.0)long_e_model = long_e_model-360.

! calculate lat/longitude for points on equatorial grid for B grid
          DO j = 1, n_rows
            gj = datastart(2) + j - 1
            DO i = 1, row_length
              gi = datastart(1) + i - 1
              longitude_rot(i,j) = (base_lambda + (gi-.5) *               &
                                         delta_lambda) / Pi_over_180
              latitude_rot(i,j) = (base_phi + (gj-.5) * delta_phi)        &
                                                      / Pi_over_180
            END DO
          END DO

        CALL eqtoll(                                                      &
                    latitude_rot, longitude_rot,                          &
                    true_latitude_b, true_longitude_b,                    &
                    lat_rot_NP_deg, long_rot_NP_deg,                      &
                    n_rows*row_length)

! Calculate rotation coefficients for wind

!!! Allocate arrays for rot_coeff_mod module
        IF (.NOT. ALLOCATED(rot_coeff1)) THEN
          ALLOCATE (rot_coeff1 ( row_length, n_rows ))
        END IF
        IF (.NOT. ALLOCATED(rot_coeff2)) THEN
          ALLOCATE (rot_coeff2 ( row_length, n_rows ))
        END IF

        CALL w_coeff(rot_coeff1, rot_coeff2,                              &
                     true_longitude_b, longitude_rot,                     &
                     lat_rot_NP_deg, long_rot_NP_deg,                     &
                     n_rows*row_length )

      END IF ! model_domain


     IF (l_cartesian) THEN
       true_latitude    = 0.0
       true_longitude   = 0.0
       true_latitude_b  = 0.0
       true_longitude_b = 0.0
       longitude_rot    = 0.0
       latitude_rot     = 0.0
     END IF


! ------------------------------------------------
! 3. swap_bounds for fields which havee halos
! ------------------------------------------------

      i_field = 0

      i_field = i_field + 1
      fields_to_swap(i_field) % field_2d   => fv_cos_theta_latitude(:,:)
      fields_to_swap(i_field) % field_type = fld_type_p
      fields_to_swap(i_field) % levels     = 1
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field_2d   => cos_theta_latitude(:,:)
      fields_to_swap(i_field) % field_type = fld_type_p
      fields_to_swap(i_field) % levels     = 1
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field_2d   => sec_theta_latitude(:,:)
      fields_to_swap(i_field) % field_type = fld_type_p
      fields_to_swap(i_field) % levels     = 1
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field_2d   => fv_sec_theta_latitude(:,:)
      fields_to_swap(i_field) % field_type = fld_type_p
      fields_to_swap(i_field) % levels     = 1
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field_2d   => cos_v_latitude(:,:)
      fields_to_swap(i_field) % field_type = fld_type_v
      fields_to_swap(i_field) % levels     = 1
      fields_to_swap(i_field) % rows       = n_rows
      fields_to_swap(i_field) % vector     = .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field_2d   => sec_v_latitude(:,:)
      fields_to_swap(i_field) % field_type = fld_type_v
      fields_to_swap(i_field) % levels     = 1
      fields_to_swap(i_field) % rows       = n_rows
      fields_to_swap(i_field) % vector     = .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field_2d   => f3_at_u(:,:)
      fields_to_swap(i_field) % field_type = fld_type_u
      fields_to_swap(i_field) % levels     = 1
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field_2d   => f3_at_v(:,:)
      fields_to_swap(i_field) % field_type = fld_type_v
      fields_to_swap(i_field) % levels     = 1
      fields_to_swap(i_field) % rows       = n_rows
      fields_to_swap(i_field) % vector     = .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field_2d   => f1_at_v(:,:)
      fields_to_swap(i_field) % field_type = fld_type_v
      fields_to_swap(i_field) % levels     = 1
      fields_to_swap(i_field) % rows       = n_rows
      fields_to_swap(i_field) % vector     = .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field_2d   => f2_at_u(:,:)
      fields_to_swap(i_field) % field_type = fld_type_u
      fields_to_swap(i_field) % levels     = 1
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.

! DEPENDS ON: swap_bounds_2d_mv
      CALL swap_bounds_2d_mv(fields_to_swap, i_field, row_length,         &
                              offx, offy)

      IF (integrity_test)                                                 &
         CALL update_hash_m(                                              &
          cos_theta_latitude, SIZE(cos_theta_latitude), 'cthla')

      IF (lhook) CALL dr_hook('SET_TRIGS_4A',zhook_out,zhook_handle)
 
      END SUBROUTINE set_trigs_4A
