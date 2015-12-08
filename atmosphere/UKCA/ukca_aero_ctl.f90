! *****************************COPYRIGHT*******************************
! 
! (c) [University of Leeds] [2008]. All rights reserved. 
! This routine has been licensed to the Met Office for use and 
! distribution under the UKCA collaboration agreement, subject  
! to the terms and conditions set out therein. 
! [Met Office Ref SC138]  
! 
! *****************************COPYRIGHT*******************************
!
!  Description:
!     UKCA-MODE aerosol code: interface routine called from
!     UKCA_MAIN1 to perform a 1-timestep integration.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
! Subroutine Interface:
      SUBROUTINE UKCA_AERO_CTL(i_month, i_day_number,                   &
                      i_hour,i_minute, DTC,                             &
                      model_levels, rows, row_length,                   &
                      wet_levels,                                       &
                      global_row_length,global_rows,                    &
                      n_chemistry_tracers,                              &
                      n_mode_tracers,                                   &
                      het_dimn, nhet_std_trop,                          &
                      area,                                             &
                      sinlat,                                           &
                      coslat,                                           &
                      true_longitude,                                   &
                      pres,                                             &
                      temp,                                             &
                      q,                                                &
                      rh3d,                                             &
                      p_bdrs,                                           &
                      chemistry_tracers,                                &
                      mode_tracers,                                     &
                      bl_levels,                                        &
                      t_surf,                                           &
                      sea_ice_frac,                                     &
                      z0m,                                              &
                      u_s,                                              &
                      u_10m,                                            &
                      drain, crain,                                     &
                      land_fraction,                                    &
                      nbox,                                             &
                      delso2_wet_h2o2,                                  &
                      delso2_wet_o3,                                    &
                      delso2_dry_oh,                                    &
                      mode_diags,                                       &
                      het_rates,                                        &
                      emissions,                                        &
                      em_index,                                         &
                      SO2_volc_3D,                                      &
                      BC_biom_3D,                                       &
                      OC_biom_3D,                                       &
                      cloud_frac,                                       &
                      cloud_liq_frac,                                   &
                      cloud_liq_wat,                                    &
                      offx, offy,                                       &
                      r_rho_levels, r_theta_levels,                     &
                      z_half, z_half_alllevs, ml_depth, delta_r,        &
                      rhokh_mix,                                        &
                      dtrdz_charney_grid, kent, we_lim,                 &
                      t_frac, zrzi, kent_dsc,                           &
                      we_lim_dsc, t_frac_dsc,                           &
                      zrzi_dsc, zhsc,                                   &
                      volume,mass,zbl,                                  &
                      DRYOX_IN_AER,                                     &
                      WETOX_IN_AER,                                     &
                      chem_diag_cdnc,                                   &
                      aerosol_surface_area                              &
                      )

      USE UKCA_D1_DEFS,     ONLY: n_mode_diags, nmax_mode_diags,        &
                                  ukca_sect, Nukca_D1items,             &
                                  ukcaD1codes, mode_diag_sect,          &
                                  item1_mode_diags, L_ukca_mode_diags,  &
                                  n_chem_emissions, nm_spec, tr_index,  &
                                  n_chem_tracers, n_aero_tracers
      USE UKCA_CONSTANTS,   ONLY: PPI, AVC, ZBOLTZ, VKARMN, RA, RR,     &
                                  GG, RAD_E, MM_DA, MMSUL, NMOL,        &
                                  EMS_EPS, CONC_EPS, DN_EPS,            &
                                  RHOSUL
      USE UKCA_CSPECIES,    ONLY: n_h2so4,n_h2o2,n_so2,n_o3,n_sec_org
      USE ukca_option_mod,  ONLY: i_mode_setup, i_mode_ss_scheme,       &
                                  i_mode_nucscav,                       &
                                  i_mode_nzts, i_mode_bln_param_method, &
                                  L_ukca_primsu,                        &
                                  L_ukca_primbcoc, L_bcoc_ff, L_bcoc_bf,&
                                  L_bcoc_bm,                            &
                                  L_mode_bhn_on, L_mode_bln_on,         &
                                  L_ukca_arg_act, L_ukca_trophet,       &
                                  mode_parfrac, l_ukca, l_ukca_aie1,    &
                                  l_ukca_aie2
      USE UKCA_MODE_SETUP
      USE UKCA_SETUP_INDICES
      USE ASAD_MOD,  ONLY: advt
      USE run_aerosol_mod,  ONLY: SO2_high_level
!
      USE bl_option_mod,    ONLY: alpha_cd 
      USE parkind1,         ONLY: jprb, jpim 
      USE yomhook,          ONLY: lhook, dr_hook 
      USE ereport_mod,      ONLY: ereport 
      USE PrintStatus_mod 
      USE Control_Max_Sizes 
      USE tr_mix_mod,            ONLY: tr_mix

      IMPLICIT NONE

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
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! Inputs
      INTEGER, INTENT(IN) :: i_month           ! month
      INTEGER, INTENT(IN) :: i_day_number      ! day
      INTEGER, INTENT(IN) :: i_hour            ! hour
      INTEGER, INTENT(IN) :: i_minute          ! minute
      INTEGER, INTENT(IN) :: model_levels      ! # of model levels
      INTEGER, INTENT(IN) :: wet_levels        ! # of wet levels
      INTEGER, INTENT(IN) :: rows              ! # of rows in patch
      INTEGER, INTENT(IN) :: row_length        ! # of pts in a patch row
      INTEGER, INTENT(IN) :: global_rows       ! # of rows (global)
      INTEGER, INTENT(IN) :: global_row_length ! # of pts in a row (global)
      INTEGER, INTENT(IN) :: n_chemistry_tracers ! # of aerosol precursor gas tracers
      INTEGER, INTENT(IN) :: n_mode_tracers    ! # of mode tracers
      INTEGER, INTENT(IN) :: bl_levels         ! # of levels in BL
      INTEGER, INTENT(IN) :: nbox              ! dimension of slice
      INTEGER, INTENT(IN) :: het_dimn          ! 1st dimension of het_rates array
      INTEGER, INTENT(IN) :: nhet_std_trop     ! No. of het. equations on aerosol
      INTEGER, INTENT(IN) :: offx              ! standard halo size
      INTEGER, INTENT(IN) :: offy              ! standard halo size
      INTEGER, INTENT(IN) :: kent(1:row_length, 1:rows)
      INTEGER, INTENT(IN) :: kent_dsc(1:row_length, 1:rows)
      INTEGER, INTENT(IN) :: em_index(n_chem_emissions)

! switch for updating of condensables in MODE or in UKCA-CHEMISTRY
      INTEGER, INTENT(IN) :: DRYOX_IN_AER      ! 0 external, 1 internal
! switch for doing aqueous SO4 production in MODE or in UKCA-CHEMISTRY
      INTEGER, INTENT(IN) :: WETOX_IN_AER      ! 0 external, 1 internal

      REAL, INTENT(IN) :: dtc                                ! timestep(s)
      REAL, INTENT(IN) :: area(row_length,rows,model_levels) ! area (m^2)
      REAL, INTENT(IN) :: sinlat(row_length, rows)           ! sin(latitude)
      REAL, INTENT(IN) :: coslat(row_length, rows)           ! cos(latitude)
      REAL, INTENT(IN) :: true_longitude(row_length,rows)    ! longitude
      REAL, INTENT(IN) :: pres(row_length,rows,model_levels) ! pressure
      REAL, INTENT(IN) :: temp(row_length,rows,model_levels) ! temperature
      REAL, INTENT(IN) :: q(row_length,rows,wet_levels)      ! sp humidity
      REAL, INTENT(IN) :: rh3d(row_length,rows,wet_levels)   ! RH (frac)
      REAL, INTENT(IN) :: p_bdrs(row_length,rows,0:model_levels) ! pressure on interfaces
      REAL, INTENT(IN) :: emissions(row_length,rows,n_chem_emissions) ! 2D emissions fields
! 3-D volcanic SO2 emissions, Biomass burning emissions for BC and OC
      REAL, INTENT(IN) :: SO2_volc_3D(row_length,rows,model_levels)
      REAL, INTENT(IN) :: BC_biom_3D(row_length,rows,model_levels)
      REAL, INTENT(IN) :: OC_biom_3D(row_length,rows,model_levels)
      REAL, INTENT(IN) :: t_surf(row_length, rows)           ! surface temp
      REAL, INTENT(IN) :: sea_ice_frac(row_length, rows)     ! sea ice
      REAL, INTENT(IN) :: u_s(row_length, rows)              ! friction velocity
      REAL, INTENT(IN) :: u_10m(row_length, rows)            ! wind at 10m
      REAL, INTENT(IN) :: z0m(row_length, rows)              ! roughness length
      REAL, INTENT(IN) :: drain(row_length,rows, model_levels) ! 3-D LS rain rate
      REAL, INTENT(IN) :: crain(row_length,rows, model_levels) ! 3-D conv rain
      REAL, INTENT(IN) :: land_fraction(row_length,rows)     ! land_fraction
! in-cloud oxidation rates (molecules/cc/DTC) from h2o2 & o3 (UKCA):
      REAL, INTENT(IN) :: delso2_wet_h2o2(row_length,rows,model_levels)
      REAL, INTENT(IN) :: delso2_wet_o3  (row_length,rows,model_levels)
! in-air   oxidation rate  (molecules/cc/DTC) from oh        (UKCA):
      REAL, INTENT(IN) :: delso2_dry_oh   (row_length,rows,model_levels)
! in-cloud oxidation rates (kgS/kgair/s     ) from h2o2 & o3 (CLASSIC):
!      REAL, INTENT(IN) :: delso2_wet_h2o2C(row_length,rows,model_levels)
!      REAL, INTENT(IN) :: delso2_wet_o3C  (row_length,rows,model_levels)
! in-air   oxidation rate  (kgS/kgair/s     ) from oh        (CLASSIC):
!      REAL, INTENT(IN) :: delso2_dry_ohC  (row_length,rows,model_levels)
! cloud fraction
      REAL, INTENT(IN) :: cloud_frac(row_length, rows, wet_levels)
      REAL, INTENT(IN) :: cloud_liq_frac(row_length, rows, wet_levels)
      REAL, INTENT(IN) :: cloud_liq_wat(row_length, rows, wet_levels)
      REAL, INTENT(IN) :: volume(row_length,rows, model_levels)
      REAL, INTENT(IN) :: mass(row_length,rows, model_levels)
      REAL, INTENT(IN) :: zbl(row_length, rows)  ! BL height

! gas phase aerosol precursor tracer mass mixing ratios
      REAL, INTENT(INOUT) :: chemistry_tracers(row_length,rows,           &
                                  model_levels,n_chemistry_tracers)
! aerosol tracers mass mixing ratio
      REAL, INTENT(INOUT) :: mode_tracers(row_length,rows,              &
                                  model_levels,n_mode_tracers)
! 3-D diagnostic array
      REAL, INTENT(INOUT) :: mode_diags(row_length,rows,                &
                                  model_levels,n_mode_diags)

! 1-D diagnostic to hold heterogenous rates fro tropospheric chemistry
      REAL, INTENT(OUT)   :: het_rates(het_dimn,nhet_std_trop)

! Items for tr_mix call:
      REAL, INTENT(IN) :: r_rho_levels(1:row_length, 1:rows,            &
                             1:model_levels)        ! ht of rho levs
      REAL, INTENT(IN) :: r_theta_levels(1:row_length, 1:rows,          &
                             0:model_levels)        ! ht of theta levs
      REAL, INTENT(IN) :: z_half_alllevs(1:row_length,1:rows,           &
                                         1:model_levels)
      REAL, INTENT(IN) :: z_half(1:row_length,1:rows,1:bl_levels)
      REAL, INTENT(IN) :: ml_depth(1:row_length,1:rows)
      REAL, INTENT(IN) :: delta_r(row_length,rows,model_levels)
      REAL, INTENT(IN) :: rhokh_mix(1:row_length,1:rows,1:bl_levels)
      REAL, INTENT(IN) :: dtrdz_charney_grid(1:row_length,1:rows,       &
                                               1:bl_levels)
      REAL, INTENT(IN) :: we_lim(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: t_frac(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: zrzi(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: we_lim_dsc(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: t_frac_dsc(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: zrzi_dsc(1:row_length,1:rows,1:3)
      REAL, INTENT(IN) :: zhsc(1:row_length,1:rows)
! Add variable to chem_diags array to store CDNC in dump for indirect effects
      REAL, INTENT(INOUT):: chem_diag_cdnc(row_length,rows,model_levels)
      REAL, INTENT(OUT) :: aerosol_surface_area(1:row_length,1:rows,    &
                                  1:model_levels)

! Local variables
      INTEGER, PARAMETER :: NMTS=1
! No. of microphysical sub-steps per DTC
      INTEGER :: NZTS
! No. of condensation-nucleation competition sub-steps per DTM
      INTEGER :: PRIMSU_ON
! Switch for whether primary sulfate particle emissions are on/off
      INTEGER :: PRIMBCOC_ON
! Switch for whether primary carbonaceous particle emissions are on/off
      INTEGER, PARAMETER :: RAINOUT_ON=1
! Switch for whether rainout (nucl. scav.) is on/off
      INTEGER, PARAMETER :: IMSCAV_ON=1
! Switch for whether impaction scavenging is on/off
      INTEGER, PARAMETER :: WETOX_ON=1
! Switch for whether wet oxidation (cloud processing) is on/off
      INTEGER, PARAMETER :: DDEPAER_ON=1
! Switch for whether aerosol dry deposition is on/off
      INTEGER, PARAMETER :: SEDI_ON=1
! Switch for whether aerosol sedimentation is on/off
      INTEGER, PARAMETER :: ISO2WETOXBYO3=1
! Switch for whether SO2 wet oxidation by ozone is on/off
! Note that this switch is only used if WETOX_IN_AER=1
! When code used in UM    , WETOX_IN_AER is always set to 0
! When code used in TOMCAT, WETOX_IN_AER is always set to 1
      INTEGER, PARAMETER :: COND_ON=1
! Switch for whether vapour condensation is  on/off
      INTEGER :: NUCL_ON
! Switch for whether binary nucleation is on/off
      INTEGER :: BLN_ON
! Switch for whether binary BL nucleation is on/off
      INTEGER, PARAMETER :: COAG_ON=1
! Switch for whether coagulation is on/off
      INTEGER, PARAMETER :: ICOAG=1
! Switch for KIJ method (1:GLOMAP, 2: M7, 3: UMorig, 4:UMorig MFPP)
!   =3 Cunnigham scheme as in UM, =4 as in UM but computing values)
      INTEGER, PARAMETER :: IMERGE=2
! Switch to use mid-pts (=1), edges (2) or dynamic (=3) in remode
      INTEGER, PARAMETER :: IFUCHS=2
! Switch for Fuchs(1964) (=1) or Fuchs-Sutugin(1971) for CC (=2)
      INTEGER, PARAMETER :: IDCMFP=2
! Switch for vapour-diffusion-method (1=as bin v1, 2=as bin v1.1)
      INTEGER, PARAMETER :: ICONDIAM=2
! Switch for what diameter to use for CONDEN (1=g.m.diam, 2=conden-diam)
      INTEGER, PARAMETER :: I_NUC_METHOD=2
!  I_NUC_METHOD: Switch for nucleation (how to combine BHN and BLN)
! (1=initial Pandis94 approach (no BLN even if switched on) -- Do not use!!
! (2=binary homogeneous nucleation applying BLN to BL only if switched on)
!   note there is an additional switch i_bhn_method (local to CALCNUCRATE)
!   to switch between using Kulmala98 or Vehkamakki02 for BHN rate
! (3=use same rate at all levels either activation(IBLN=1), kinetic(IBLN=2),
!  PNAS(IBLN=3),eucaari-kinetic(IBLN=4),eucaari-org1(IBLN=5),eucaari-org2(IBLN=6)
!  note that if I_NUC_METHOD=3 and IBLN=3 then also add on BHN rate as in PNAS.
      INTEGER            :: IBLN  
! Switch for BLN parametrisation rate calc (1=activation,2=kinetic,3=PNAS)
      INTEGER :: IACTMETHOD
! Switch for activation method (0=off,1=fixed ract,2=NSO3 scheme)
      INTEGER, PARAMETER :: IWVOLMETHOD=2
! Switch for wet volume method (1=behave-as-H2SO4,2=multi-cpt)
      INTEGER :: INUCSCAV
! Switch for nucl scav method (1=as GLOMAP Spr05, 2=use scav coeffs)
      INTEGER, PARAMETER :: IDDEPAER=2
! Switch for dry dep method (1=as GLOMAP Spr05, 2=incl. sedi)
!      INTEGER, PARAMETER :: IDDEPAER=1
! Switch for dry dep method (1=as GLOMAP Spr05, 2=incl. sedi)
      INTEGER :: VERBOSE
! Switch to determine level of debug o/p (0=none, 1, 2)
! Now set via printstatus from UMUI Output Choices panel, see below
      INTEGER, PARAMETER :: CHECKMD_ND=0
! Switch for whether to check for bad values of MD and ND
      INTEGER, PARAMETER :: INTRAOFF=0
! Switch to turn off intra-modal coagulation
      INTEGER, PARAMETER :: INTEROFF=0
! Switch to turn off inter-modal coagulation
      INTEGER, PARAMETER :: IDUSTEMS=0
! Switch for using Pringle scheme (=1) or AEROCOMdaily (=2)

      INTEGER  :: LDAY(NBOX)
! Switch for day/night (1/0)
      INTEGER  :: JLABOVE(NBOX)
! Index of box directly above this grid box
      INTEGER  :: ILSCAT(NBOX)
! Land surface category (based on 9 UM landsurf types)
      INTEGER  :: IARR(NBOX)
! Index of LAT for grid box
      INTEGER  :: KARR(NBOX)
! Index of LON for grid box
      INTEGER  :: LARR(NBOX)
! Index of vertical level for grid box

      REAL, PARAMETER  :: emfactor=1.4
! To convert emissions of OC and BC from kg[OM] to kg[C]

      REAL     :: DLAT3D(row_length,rows,model_levels) ! 3D lat array
      REAL     :: DLON3D(row_length,rows,model_levels) ! 3D lon array

      REAL     :: EMANSO2(NBOX,6)
! Anthrop. SO2 ems rates, low sources (kgSO2/box/s)
      REAL     :: EMVOLCONSO2(NBOX)
! Volcanic SO2 ems rates (cont. src) (kgSO2/box/s)
      REAL     :: EMVOLEXPSO2(NBOX)
! Volcanic SO2 ems rates (expl. src) (kgSO2/box/s)
      REAL     :: EMBIOMSO2(NBOX)
! Biomass SO2 ems rates (kgSO2/box/s)
      REAL     :: EMC(NBOX,4)
! BC/OC emission rates from bio- & fossil-fuels (kgC/box/s)
      REAL     :: EMCBM(NBOX,2)
! BC/OC emission rates from biomass burning (kgC/box/s)
      REAL     :: ND(NBOX,NMODES)
! Aerosol ptcl number density for mode (cm^-3)
      REAL     :: MDT(NBOX,NMODES)
! Avg tot mass of aerosol ptcl in mode (particle^-1)
      REAL     :: MDWAT(NBOX,NMODES)
! Molecular concentration of water (molecules per particle)
      REAL     :: DRYDP(NBOX,NMODES)
! Geometric mean dry diameter of particles in each mode (m)
      REAL     :: WETDP(NBOX,NMODES)
! Geometric mean wet diameter of particles in each mode (m)
      REAL     :: RHOPAR(NBOX,NMODES)
! Total particle density [incl. H2O & insoluble cpts] (kgm^-3)
      REAL     :: DVOL(NBOX,NMODES)
! Geometric mean dry volume of particles in each mode (m^3)
      REAL     :: WVOL(NBOX,NMODES)
! Geometric mean wet volume of particles in each mode (m^3)
      REAL     :: MD(NBOX,NMODES,NCP)
! Avg cpt mass of aerosol particle in mode (particle^-1)
      REAL     :: S0(NBOX,NADVG)
! Partial masses of gas phase species (kg per gridbox)
      REAL     :: S0_DOT_CONDENSABLE(NBOX,NCHEMG)
! ASAD tendencies for condensable gas phase species (vmr per s)
      REAL     :: SM(NBOX)     ! Mass of air in gridbox (kg)
! Grid box mass of air (kg)
      REAL     :: AIRD(NBOX)
! Number density of air (per cm3)
      REAL     :: AIRDM3(NBOX)
! Number density of air (per m3)
      REAL     :: RHOA(NBOX)
! Air density (kg/m3)
      REAL     :: VBA(NBOX)
! Mean free speed of air molecules (m/s)
      REAL     :: TSQRT(NBOX)
! Square-root of centre level temperature (K)
      REAL     :: DVISC(NBOX)
! Dynamic viscosity of air (kg m^-1 s^-1)
      REAL     :: MFPA(NBOX)
! Mean free path of air (m)
      REAL     :: T(NBOX)
! Air temperature at mid-point (K)
      REAL     :: RH(NBOX)
! Relative humidity (fraction)
      REAL     :: S(NBOX)
! Specific humidity (kg/kg)
      REAL     :: PMID(NBOX)
! Air pressure at mid-point (Pa)
      REAL     :: PUPPER(NBOX)
! Air pressure at upper interface (Pa)
      REAL     :: PLOWER(NBOX)
! Air pressure at lower interface (Pa)
      REAL     :: ZO3(NBOX)
! Background vmr of O3 (dimensionless)
      REAL     :: ZHO2(NBOX)
! Background conc. of HO2 (molecules per cc)
      REAL     :: ZH2O2(NBOX)
! Background conc. of H2O2 (molecules per cc)
      REAL     :: USTR(NBOX)
! Surface friction velocity (m/s)
      REAL     :: US10M(NBOX)
! Scalar wind at 10m (ms-1)
      REAL     :: ZNOTG(NBOX)
! Roughness length (m)
      REAL     :: DELTA_Z(NBOX)
! Grid-box thickness (m)
      REAL     :: SURTP(NBOX)
! Surface type: 0=seasurf,1=landsurf,2=above-seasurf,3=above-landsurf
      REAL     :: SURF(NBOX)
! Surface area of box (horizontal) (m^2)
      REAL     :: LAND_FRAC(NBOX)
! Fraction of horizontal gridbox area covered by land
      REAL     :: SEAICE(NBOX)
! Fraction of horizontal gridbox area containing seaice
      REAL     :: CRAING(NBOX)
! Rain rate for conv precip. in box (kgm^-2s^-1)
      REAL     :: DRAING(NBOX)
! Rain rate for dyn. precip. in box (kgm^-2s^-1)
      REAL     :: CRAING_UP(NBOX)
! Rain rate for conv precip. in box above (kgm^-2s^-1)
      REAL     :: DRAING_UP(NBOX)
! Rain rate for dyn. precip. in box above (kgm^-2s^-1)
      REAL     :: FCONV_CONV(NBOX)
! Fraction of box condensate --> rain in 6 hours (conv)
      REAL     :: LOWCLOUD(NBOX)
! Horizontal low cloud fraction
      REAL     :: VFAC(NBOX)
! Vertical low cloud fraction
      REAL     :: LWC(NBOX)
! Cloud liquid water content [kg/m3]
      REAL     :: CLWC(NBOX)
! Cloud liquid water content [kg/kg]
      REAL     :: CLF(NBOX)
! Liquid Cloud fraction
      REAL :: PVOL(NBOX,NMODES,NCP)
! Aerosol partial volume of each cpt in each mode
      REAL :: PVOL_WAT(NBOX,NMODES)
! Aerosol partial volume of water in each mode
      REAL :: tr_rs(NBOX)
! Local variable to hold re-shaped aerosol tracers
!
      REAL     :: DTM
! Microphysics time step (s)
      REAL     :: DTZ
! Competition (cond/nucl) time step (s)
      REAL     :: RMOIS
! Month of year
      REAL     :: RJOUR
! Day of month
      REAL     :: DLON
! Delta longitude
      REAL     :: DLAT
! Delta latitude
      REAL     :: FAC(NBOX)
! Conversion factor for budget diagnostics
      REAL     :: FAC_MMRCONV(NBOX)
! Conversion factor to convert MMRSO2/s to molecules/cc/DTC

      REAL     :: SILT(NBOX)
! Silt fraction for gridbox
      REAL     :: SAND(NBOX)
! Sand fraction for gridbox
      REAL     :: CLAY(NBOX)
! Clay fraction for gridbox
      REAL     :: SWC(NBOX)
! Soil water content for gridbox
      REAL     :: SNOWICE(NBOX)
! Snow/ice fraction for gridbox
      REAL     :: PSOURCE(NBOX)
! Preferential source grid
      REAL     :: LAI(NBOX)
! Leaf area index
      REAL     :: TEX(NBOX)
! Soil texture
      REAL     :: MODE2_RADIUS(NBOX,10)
! Mode 2 g.m. radius (for AEROCOM prescribed daily dust fluxes)
      REAL     :: MODE3_RADIUS(NBOX,10)
! Mode 3 g.m. radius (for AEROCOM prescribed daily dust fluxes)
      REAL     :: MODE2_NUMBER(NBOX,10)
! Mode 2 number conc (for AEROCOM prescribed daily dust fluxes)
      REAL     :: MODE3_NUMBER(NBOX,10)
! Mode 3 number conc (for AEROCOM prescribed daily dust fluxes)
      INTEGER  :: NDUSTEMINBOX(NBOX)
! No. of 1x1 emissions in box (for AEROCOM daily dust fluxes)
      REAL     :: DLONARR(NBOX)
! Longitude (radians) of gridbox
      REAL     :: DLATARR(NBOX)
! Latitude  (radians) of gridbox
      REAL     :: HEIGHT(NBOX)
! Mid-level height of gridbox
      REAL     :: HTPBLG(NBOX)
! Height of boundary-layer in gridbox vertical-column
      REAL     :: MDTFIXFLAG(NBOX,NMODES)
! Array to find where MDT is too low (and ND->0 is being applied).
! Gets to set = 100.0 when fix applied so that gives percentage of
! timesteps where fix is applied in each gridbox when take mean
! Store this in STASH by temporarily over-writing a budget diag
      INTEGER  :: myproc
! Index of PE element (only used for debugging -- dummy as 1 for now)

      REAL     :: Y(NMODES)
! EXP(2*LN(SIGMA)*LN(SIGMA)) for each mode (for surf area conc calc)
      REAL     :: SAREA(NBOX,NMODES)
! Surface area concentration for each mode (cm^2 / cm^3)
      REAL     :: VCONC(NBOX,NMODES)
! Volume concentration for each mode
      REAL     :: CN_3NM(NBOX)
! CN concentration (dry diameter > 3nm)
      REAL     :: CCN_1(NBOX)
! CCN concentration (acc-sol + cor-sol)
      REAL     :: CCN_2(NBOX)
! CCN concentration (acc-sol + cor-sol + Aitsol>25nm drydp)
      REAL     :: CCN_3(NBOX)
! CCN concentration (acc-sol + cor-sol + Aitsol>35nm drydp)
      REAL     :: CDN(NBOX)
! CDN concentration

! Outputs for budget calculations
      REAL :: BUD_AER_MAS(NBOX,0:NBUDAER)

! Items for tr_mix:
      REAL :: res_factor(row_length,rows)                       ! Ra/(Ra+Rb+Rc)
      REAL :: rho_aresist(1:row_length,1:rows)
      REAL :: tr_flux(row_length,rows,bl_levels,n_mode_tracers) ! tracer fluxes
      REAL :: surf_dep_flux(row_length,rows,n_mode_tracers)     ! surf depn. flux
      REAL :: em_field(row_length,rows,n_mode_tracers)          ! surface emissions

! Items for CN,CCN,CDN calculation
      REAL :: ERFTERM(NBOX)
      REAL :: ERF_ARG(NBOX)                                     ! Error Fn argument
      REAL :: DP0                                               ! Diam (nm)
      REAL, PARAMETER :: CDNMIN = 5.0                           ! Min CDN (no/cm^-3)
!               5.0 for marine & ice-continental
!               35.0 for ice-free continental

! Local declaration for molar mass of sulphur
      REAL :: MM_SULPHUR
      
! Items for MDT too low/hi check
      REAL :: MDTMIN
      REAL :: SUMMDMODE(NBOX)

      INTEGER :: field_size              ! size of 2D field
      INTEGER :: field_size3d            ! size of 3D field
      REAL, ALLOCATABLE :: z(:,:,:)   ! temporary array

      INTEGER :: I,J,K,L,N,JL,IMODE,ICP
      INTEGER :: n_reqd_tracers           ! No of tracers required
      INTEGER :: ITRA

      INTEGER :: ISO2EMS ! setting for emissions (now from namelist)

! Now use ISO2EMS=1 when chosen only sulphate and sea-salt in 4 modes
!                   (no BC/OC) -- i_mode_setup=1 -- small sizes
!                   make up for lack of primary BC/OC emissions
! 
! For all other options use ISO2EMS=3 as standard as in GLOMAP.
!
! ISO2EMS == 1 : 3.0% of SO2 ems --> primary SO4 ems --> 15%/85% to
!                10/70 nm g.m.diam. modes as in Spracklen (2005)
!                and Binkowski & Shankar (1995). (UMUI value not used)
!
! ISO2EMS == 2 : 2.5% of SO2 ems --> primary SO4 ems (Now set in UMUI - mode_parfrac)
!                road/off-road/domestic      all -->   30nm gm.diam mode
!                industrial/power-plant/ship all --> 1000nm gm.diam mode
!                as for original AEROCOM size recommendations.
!
! ISO2EMS >= 3 : 2.5% of SO2 ems --> primary SO4 ems (Now set in UMUI - mode_parfrac)
!                --> 50%/50% to
!                150/1500nm  g.m.diam. modes as for Stier et al (2005)
!                modified AEROCOM sizdis recommendations.
!
! Note also that currently, ISO2EMS also controls size assumptions for
! primary carbonaceous aerosol (this needs to be changed):
!
! ISO2EMS /= 3 : biofuel & biomass BC/OC emissions -->  80nm g.m.diam.
!                fossil-fuel BC/OC emissions --> 30nm g.m.diam.
!
! ISO2EMS == 3 : biofuel & biomass BC/OC emissions --> 150nm g.m.diam.
!                fossil-fuel BC/OC emissions --> 60nm g.m.diam.
!
      REAL    :: PARFRAC
! PARFRAC is fraction of mass of SO2 emitted as primary SO4
      REAL    :: DELSO2(NBOX)
!  S(IV) --> S(VI) by H2O2 (molecules per cc) [input if WETOX_IN_AER=0]
      REAL    :: DELSO2_2(NBOX)
!  S(IV) --> S(VI) by O3   (molecules per cc) [input if WETOX_IN_AER=0]

      REAL, PARAMETER :: MA=4.78E-26 ! mass of air molecule (kg)
      REAL :: ACT                    ! radius for activation (m)

! used for debug output
      LOGICAL,           ALLOCATABLE, SAVE :: mode_tracer_debug(:)
      CHARACTER(LEN=10), ALLOCATABLE, SAVE :: chemistry_tracer_names(:)
      CHARACTER(LEN=10), ALLOCATABLE, SAVE :: mode_tracer_names(:)

      CHARACTER(LEN=72) :: cmessage     ! Error message

      LOGICAL :: MASK1(NBOX)
      LOGICAL :: MASK_NOTLC(NBOX)
      LOGICAL :: MASK_NDGTE(NBOX)
      LOGICAL :: MASK_TRLT0(NBOX)
      LOGICAL :: MASK_NDLOW(NBOX)
      LOGICAL :: LOGIC  ! for use in storing of aerosol budget terms
      LOGICAL :: LOGIC1 ! for use in storing of aerosol budget terms
      LOGICAL :: LOGIC2 ! for use in storing of aerosol budget terms
      LOGICAL, SAVE :: firstcall=.TRUE.
      INTEGER :: N_MERGE_1D(NBOX,NMODES) ! counter: mode-merges applied
      INTEGER :: N_MERGE_3D(row_length,rows,model_levels,NMODES)
      INTEGER :: NBADMDT(NBOX,NMODES)
      INTEGER :: IDX(NBOX)
      INTEGER :: NFIX
      INTEGER :: JV
      INTEGER :: IA
      INTEGER, SAVE :: lBCff    ! Index for BC fossil fuel emissions
      INTEGER, SAVE :: lBCbf    ! Index for BC biofuel emissions
      INTEGER, SAVE :: lOCff    ! Index for OC fossil fuel emissions
      INTEGER, SAVE :: lOCbf    ! Index for OC biofuel emissions
      INTEGER, SAVE :: lso2emlo ! Index for SO2-low emissions
      INTEGER, SAVE :: lso2emhi ! Index for SO2-low emissions

! These taken out of run_ukca as set here for now
      INTEGER :: i_mode_act_method
      REAL    :: mode_actdryr

      INTEGER :: errorstatus
      INTEGER :: errcode
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0 
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1 
      REAL(KIND=jprb)               :: zhook_handle 

      IF (lhook) CALL dr_hook('UKCA_AERO_CTL',zhook_in,zhook_handle) 

!-------------------------------------
! Set VERBOSE from printstatus as defined from UMUI panel "Output choices"
      VERBOSE = MIN(PrintStatus-1,2)  ! sets it to 0-2
!      VERBOSE=0 ! set to 0 for now

! Set BL nucleation parametrisation rate calculation method
      IBLN = I_MODE_BLN_PARAM_METHOD
! Set local variables to values in run_ukca namelist

! below are parameters set in the UKCA panel for MODE

      IF (firstcall .AND. PrintStatus >= PrStatus_Normal) THEN
!
        IF (VERBOSE > 0) THEN
          WRITE(6,'(A25,I6)') 'i_mode_setup =      ', i_mode_setup
          WRITE(6,'(A25,L7)') 'L_UKCA_PRIMSU =     ', l_ukca_primsu
          WRITE(6,'(A25,F6.2)') 'mode_parfrac =      ', mode_parfrac
          WRITE(6,'(A25,I6)') 'i_mode_ss_scheme =  ', i_mode_ss_scheme
          WRITE(6,'(A25,L7)') 'L_UKCA_PRIMBCOC =   ', l_ukca_primbcoc
          WRITE(6,'(A25,L7)') 'L_BCOC_ff =         ', l_bcoc_ff
          WRITE(6,'(A25,L7)') 'L_BCOC_bf =         ', l_bcoc_bf
          WRITE(6,'(A25,L7)') 'L_BCOC_bm =         ', l_bcoc_bm
          WRITE(6,'(A25,I6)') 'i_mode_nucscav =    ', i_mode_nucscav
          WRITE(6,'(A25,I6)') 'i_mode_nzts =       ', i_mode_nzts
          WRITE(6,'(A25,L7)') 'L_mode_bhn_on =     ', L_mode_bhn_on
          WRITE(6,'(A25,L7)') 'L_mode_bln_on =     ', L_mode_bln_on
          WRITE(6,'(A25,I6)') 'i_mode_bln_param_method',                &
                               i_mode_bln_param_method
        END IF   ! VERBOSE > 0

      END IF

! set ISO2EMS according to i_mode_setup in UMUI
      IF(i_mode_setup == 1) ISO2EMS=1 ! SUSS_4mode
      IF(i_mode_setup == 2) ISO2EMS=3 ! SUSSBCOC_5mode
      IF(i_mode_setup == 3) ISO2EMS=3 ! SUSSBCOC_4mode
      IF(i_mode_setup == 4) ISO2EMS=3 ! SUSSBCOCSO_5mode
      IF(i_mode_setup == 5) ISO2EMS=3 ! SUSSBCOCSO_4mode
      IF(i_mode_setup == 6) ISO2EMS=0 ! DUonly_2mode
      IF(i_mode_setup == 7) ISO2EMS=0 ! DUonly_3mode
      IF(i_mode_setup == 8) ISO2EMS=3 ! SUSSBCOCDU_7mode
      IF(i_mode_setup == 9) ISO2EMS=3 ! SUSSBCOCDU_4mode
!
      IF(firstcall .AND. VERBOSE > 0)                                   &
        write(6,*) 'i_mode_setup,ISO2EMS=',i_mode_setup,ISO2EMS
!
! set PRIMSU_ON according to L_UKCA_PRIMSU in UMUI
      IF (L_UKCA_PRIMSU) THEN
        PRIMSU_ON = 1
      ELSE
        PRIMSU_ON = 0
      ENDIF
      IF(firstcall .AND. VERBOSE > 0) THEN
        write(6,*) 'L_UKCA_PRIMSU,PRIMSU_ON=',L_UKCA_PRIMSU,PRIMSU_ON
      ENDIF
!
! set PARFRAC according to mode_parfrac in UMUI
      PARFRAC=mode_parfrac/100.0
      IF (firstcall .AND. VERBOSE > 0)                                  &
        write(6,*) 'mode_parfrac,PARFRAC=',mode_parfrac,PARFRAC
!
      IF(firstcall .AND. VERBOSE > 0) THEN
       write(6,*) 'i_mode_ss_scheme in UMUI (as passed)=',              &
                   i_mode_ss_scheme
       IF(i_mode_ss_scheme == 1)                                        &
        write(6,*) 'Gong-Monahan scheme selected'
      ENDIF

      IF(firstcall) THEN
       IF((i_mode_ss_scheme >= 2).AND.(i_mode_ss_scheme <= 4)) THEN
        write(6,*) 'You have selected I_MODE_SS_SCHEME=',               &
                    I_MODE_SS_SCHEME
        write(6,*) 'This has not been implemented !!!!!!!'
        cmessage='Bad value of I_MODE_SS_SCHEME'
        CALL EREPORT('UKCA_AERO_CTL',I_MODE_SS_SCHEME,cmessage)
       ENDIF
      ENDIF
! 
! set PRIMBCOC_ON according to L_UKCA_PRIMBCOC in UMUI
      IF (L_UKCA_PRIMBCOC) THEN
        PRIMBCOC_ON = 1
      ELSE
        PRIMBCOC_ON = 0
      ENDIF
      IF(firstcall .AND. VERBOSE > 0) THEN
       write(6,*) 'L_UKCA_PRIMBCOC,PRIMBCOC_ON=',                       &
                   L_UKCA_PRIMBCOC,PRIMBCOC_ON
       write(6,*) 'L_BCOC_ff defined in UMUI =',                        &
                   l_bcoc_ff
       write(6,*) 'L_BCOC_bf defined in UMUI =',                        &
                   l_bcoc_bf
       write(6,*) 'L_BCOC_bm defined in UMUI =',                        &
                   l_bcoc_bm
      ENDIF
 
! set INUCSCAV according to i_mode_nucscav in UMUI
      INUCSCAV = i_mode_nucscav
      IF(firstcall .AND. VERBOSE > 0) THEN
       write(6,*) 'I_MODE_NUCSCAV,INUCSCAV=',I_MODE_NUCSCAV,INUCSCAV
      ENDIF
!
! set NZTS according to i_mode_nzts in UMUI
      NZTS = i_mode_nzts
      IF(firstcall .AND. VERBOSE > 0) THEN
        write(6,*) 'I_MODE_NZTS,NZTS=',I_MODE_NZTS,NZTS
      END IF
!
! currently (18/11/08) problem in UMUI at 7.1 which results in 
! I_MODE_ACT_METHOD not being set to the value specified in the UMUI
! so for now, set I_MODE_ACT_METHOD=1 for constant activation
! Also MODE_ACTDRYR is not being set to the value specified in the UMUI
! so for now, set MODE_ACTDRY to be 37.5nm to match that used in GLOMAP'
      IF(firstcall .AND. VERBOSE > 0) THEN
       write(6,*) 'There is currently a problem with the UMUI in that'
       write(6,*) 'UMUI-set values of MODE_ACTDRYR,I_ACT_METHOD are'//  &
                  ' not set'
       write(6,*) 'setting MODE_ACTDRYR,I_MODE_ACT_METHOD to default'// &
                  ' values'
      ENDIF
      I_MODE_ACT_METHOD=1
      MODE_ACTDRYR=37.5 ! MODE_ACTDRYR is in nm
!
! set IACTMETHOD according to i_mode_act_method in UMUI
      IACTMETHOD=i_mode_act_method
      IF(firstcall .AND. VERBOSE > 0) THEN
        write(6,*) 'I_MODE_ACT_METHOD,IACTMETHOD=',                     &
                            i_mode_act_method,iactmethod
      END IF
!
! set ACT according to mode_actdryr in UMUI
      ACT = mode_actdryr*1.0e-9 ! convert nm to m
      IF(firstcall .AND. VERBOSE > 0) THEN
       write(6,*) 'MODE_ACTDRYR,ACT=',MODE_ACTDRYR,ACT
       write(6,*) 'L_mode_bhn_on=',L_mode_bhn_on
       write(6,*) 'L_mode_bln_on=',L_mode_bln_on
       write(6,*) 'i_mode_bln_param_method',I_MODE_BLN_PARAM_METHOD
      ENDIF
!
! set NUCL_ON according to L_MODE_BHN_ON & L_MODE_BLN_ON in UMUI
      IF (L_mode_BHN_ON) THEN
        NUCL_ON = 1 ! BHN
      ELSE
        NUCL_ON = 0 ! BHN
      ENDIF

      IF (L_mode_BLN_ON) THEN
        BLN_ON = 1 ! BLN
      ELSE
        BLN_ON = 0 ! BLN
      ENDIF
!
      IF(firstcall .AND. VERBOSE > 0) THEN
       write(6,*) 'NUCL_ON,BLN_ON=',NUCL_ON,BLN_ON
       write(6,*) 'IBLN=',IBLN
      ENDIF
!
      DTM=DTC/REAL(NMTS)
      DTZ=DTM/REAL(NZTS)
      RMOIS=REAL(i_month)
      RJOUR=REAL(i_day_number)
      DLAT=PPI/REAL(global_rows)
      DLON=2.0*PPI/REAL(global_row_length)
      field_size=row_length*rows
      field_size3d=field_size*model_levels

      IF (firstcall) THEN
         ALLOCATE(chemistry_tracer_names(n_chemistry_tracers))
         ALLOCATE(mode_tracer_names(n_mode_tracers))
         ALLOCATE(mode_tracer_debug(n_mode_tracers))
         DO j=1,n_chemistry_tracers
            chemistry_tracer_names(j)=                                  &
!                 nm_spec(tr_index(j))
                 advt(j)
         ENDDO
         DO j=1,n_mode_tracers
            mode_tracer_names(j)=nm_spec(tr_index(n_chem_tracers+       &
                 n_aero_tracers+j))
         ENDDO
         mode_tracer_debug(:)=.true.     ! all tracers with debug o/p
      END IF

      IF (VERBOSE > 1) THEN

      write(6,*) 'nbox,field_size,field_size3d=',                       &
                  nbox,field_size,field_size3d

      write(6,*) 'UKCA_MODE INPUTS from UMUI : '
      write(6,*) 'i_mode_setup=',i_mode_setup
      write(6,*) 'L_UKCA_PRIMSU=',l_ukca_primsu
      write(6,*) 'mode_parfrac=',mode_parfrac
      write(6,*) 'i_mode_ss_scheme=',i_mode_ss_scheme
      write(6,*) 'L_UKCA_PRIMBCOC=',l_ukca_primbcoc
      write(6,*) 'L_BCOC_ff=',l_bcoc_ff
      write(6,*) 'L_BCOC_bf=',l_bcoc_bf
      write(6,*) 'L_BCOC_bm=',l_bcoc_bm
      write(6,*) 'i_mode_nucscav=',i_mode_nucscav
      write(6,*) 'i_mode_nzts=',i_mode_nzts
      write(6,*) 'L_mode_bhn_on=',L_mode_bhn_on
      write(6,*) 'L_mode_bln_on=',L_mode_bln_on
      write(6,*) 'i_mode_bln_param_method',I_MODE_BLN_PARAM_METHOD

      write(6,*) 'i_month: ',i_month
      write(6,*) 'i_day_number: ',i_day_number
      write(6,*) 'i_hour: ',i_hour
      write(6,*) 'i_minute: ',i_minute
      write(6,*) 'DTC: ',DTC
      write(6,*) 'model_levels: ',model_levels
      write(6,*) 'rows: ',rows
      write(6,*) 'row_length: ',row_length
      write(6,*) 'global_rows: ',global_rows
      write(6,*) 'global_row_length: ',global_row_length
      write(6,*) 'n_chemistry_tracers: ',n_chemistry_tracers
      write(6,*) 'n_mode_tracers: ',n_mode_tracers

      write(6,*) 'Array:     MIN        MAX         MEAN'
      write(6,*) 'area: ',minval(area),maxval(area),                    &
                 sum(area)/real(size(area))
      l=0
      write(6,*) 'Level: ',l
      write(6,*) 'p_bdrs: ',minval(p_bdrs(:,:,l)),maxval(p_bdrs(:,:,l)) &
                           ,sum(p_bdrs(:,:,l))/real(size(p_bdrs(:,:,l)))
      do l=1,2            ! model_levels
      write(6,*) 'Level: ',l
      write(6,*) 'pres: ',minval(pres(:,:,l)),maxval(pres(:,:,l)),      &
                 sum(pres(:,:,l))/real(size(pres(:,:,l)))
      write(6,*) 'temp: ',minval(temp(:,:,l)),maxval(temp(:,:,l)),      &
                 sum(temp(:,:,l))/real(size(temp(:,:,l)))
      write(6,*) 'q: ',minval(q(:,:,l)),maxval(q(:,:,l)),               &
                 sum(q(:,:,l))/real(size(q(:,:,l)))
      write(6,*) 'rh3d: ',minval(rh3d(:,:,l)),maxval(rh3d(:,:,l)),      &
                 sum(rh3d(:,:,l))/real(size(rh3d(:,:,l)))
      write(6,*) 'p_bdrs: ',minval(p_bdrs(:,:,l)),maxval(p_bdrs(:,:,l)) &
                           ,sum(p_bdrs(:,:,l))/real(size(p_bdrs(:,:,l)))
      write(6,*) 'delso2_wet_h2o2: ',                                   &
                 minval(delso2_wet_h2o2(:,:,l)),                        &
                 sum(delso2_wet_h2o2(:,:,l))/                           &
                 real(size(delso2_wet_h2o2(:,:,l)))
      write(6,*) 'delso2_wet_o3  : ',                                   &
                 minval(delso2_wet_o3  (:,:,l)),                        &
                 sum(delso2_wet_o3  (:,:,l))/                           &
                 real(size(delso2_wet_o3  (:,:,l)))
!      write(6,*) 'delso2_wet_h2o2C: ',                                  &
!                 minval(delso2_wet_h2o2C(:,:,l)),                       &
!                 sum(delso2_wet_h2o2C(:,:,l))/                          &
!                 real(size(delso2_wet_h2o2C(:,:,l)))
!      write(6,*) 'delso2_wet_o3C  : ',                                  &
!                 minval(delso2_wet_o3C  (:,:,l)),                       &
!                 sum(delso2_wet_o3C  (:,:,l))/                          &
!                 real(size(delso2_wet_o3C  (:,:,l)))
      enddo
      l = 8
      write(6,*) 'delso2_wet_h2o2: ',l,                                 &
                 minval(delso2_wet_h2o2(:,:,l)),                        &
                 maxval(delso2_wet_h2o2(:,:,l)),                        &
                 sum(delso2_wet_h2o2(:,:,l))/                           &
                 real(size(delso2_wet_h2o2(:,:,l)))
      write(6,*) 'delso2_wet_o3  : ',l,                                 &
                 minval(delso2_wet_o3  (:,:,l)),                        &
                 maxval(delso2_wet_o3  (:,:,l)),                        &
                 sum(delso2_wet_o3  (:,:,l))/                           &
                 real(size(delso2_wet_o3  (:,:,l)))
      write(6,*) 'delso2_dry_oh  : ',l,                                 &
                 minval(delso2_dry_oh(:,:,l)),                          &
                 maxval(delso2_dry_oh(:,:,l)),                          &
                 sum(delso2_dry_oh(:,:,l)), size(delso2_dry_oh(:,:,l))


      do l=1,2       ! model_levels
      do j=1,n_chemistry_tracers
      write(6,*) 'Level: ',l,' Tracer: ',j,chemistry_tracer_names(j)
      write(6,*) 'chemistry_tracers: ',                                 &
           minval(chemistry_tracers(:,:,l,j)),                          &
           maxval(chemistry_tracers(:,:,l,j)),                          &
           sum(chemistry_tracers(:,:,l,j))/                             &
           real(size(chemistry_tracers(:,:,l,j)))
      enddo
      do j=1,n_mode_tracers
        if (mode_tracer_debug(j)) then
         write(6,*) 'Level: ',l,' Tracer: ',j,nm_spec(100+j)
         write(6,*) 'mode_tracers: ',minval(mode_tracers(:,:,l,j)),     &
                   maxval(mode_tracers(:,:,l,j)),                       &
                   sum(mode_tracers(:,:,l,j))/                          &
                   real(size(mode_tracers(:,:,l,j)))
        endif
      enddo
      enddo     ! model_levels

      ENDIF ! IF (verbose > 1)

! Calculate number of aerosol tracers required
      n_reqd_tracers = 0 
      DO imode=1,nmodes 
        DO icp=1,ncp 
          IF (component(imode,icp)) n_reqd_tracers = n_reqd_tracers + 1 
        ENDDO 
      ENDDO 
      n_reqd_tracers = n_reqd_tracers + sum(mode_choice) 

      IF (firstcall) THEN

! .. Set indices for emissions array
        DO ITRA=1,n_chem_emissions
          IF (em_index(ITRA) == 58)  lso2emlo=ITRA
          IF (em_index(ITRA) == 126) lso2emhi=ITRA
          IF (em_index(ITRA) == 310) lBCff   =ITRA
          IF (em_index(ITRA) == 311) lBCbf   =ITRA
          IF (em_index(ITRA) == 312) lOCff   =ITRA
          IF (em_index(ITRA) == 313) lOCbf   =ITRA
        ENDDO

! .. Check the number of tracers, warn if too many, stop if too few
        IF (n_mode_tracers > n_reqd_tracers) THEN
          errcode=-1
          cmessage=' Too many tracers input'
          write(6,*) cmessage,n_mode_tracers,n_reqd_tracers
          CALL EREPORT('UKCA_AERO_CTL',errcode,cmessage)
        ENDIF
        IF (n_mode_tracers < n_reqd_tracers) THEN
          errcode=1
          cmessage=' Too few advected aerosol tracers input'
          write(6,*) cmessage,n_mode_tracers,n_reqd_tracers
          CALL EREPORT('UKCA_AERO_CTL',errcode,cmessage)
        ENDIF

      END IF   ! firstcall

      JLABOVE(:)=-1
      DO l=1,(model_levels-1) ! let top level have JLABOVE=-1
       DO i=1,rows ! this is latitude
        DO k=1,row_length ! this is longitude
         JL         =k+(i-1)*row_length+(l-1)*field_size
         JLABOVE(JL)=k+(i-1)*row_length+    l*field_size
         KARR(JL)=k ! this is longitude
         IARR(JL)=i ! this is latitude
         LARR(JL)=l
        ENDDO
       ENDDO
      ENDDO    ! l, model_levels

      N_MERGE_1D(:,:)=0
      N_MERGE_3D(:,:,:,:)=0

! Reshape input quantities
! ========================
      T     (:)=RESHAPE(temp  (:,:,:),(/field_size3d/))
      PMID  (:)=RESHAPE(pres  (:,:,:),(/field_size3d/))
      PUPPER(:)=RESHAPE(p_bdrs(:,:,1:model_levels),(/field_size3d/))
      PLOWER(:)=RESHAPE(p_bdrs(:,:,0:model_levels-1),(/field_size3d/))
      SURF  (:)=RESHAPE(area  (:,:,:),(/field_size3d/))
      CRAING(:)=RESHAPE(crain (:,:,:),(/field_size3d/))
      DRAING(:)=RESHAPE(drain (:,:,:),(/field_size3d/))
!
! .. set CRAING_UP using JLABOVE as calculated above
      DO JL=1,NBOX
       IF(JLABOVE(JL) > 0) THEN 
        CRAING_UP(JL)=CRAING(JLABOVE(JL))
        DRAING_UP(JL)=DRAING(JLABOVE(JL))
       ELSE
        CRAING_UP(JL)=0.0
        DRAING_UP(JL)=0.0
       ENDIF
      ENDDO
!
! .. currently set FCONV_CONV=0.99 -- need to change to take as input
      FCONV_CONV(:)=0.99 ! fraction of condensate-->rain in 6 hrs
! .. weakened rainout -- FCONV_CONV=0.5
!!      FCONV_CONV(:)=0.50 ! fraction of condensate-->rain in 6 hrs
!
! calculate molecular concentration of air
      AIRD(:)=PMID(:)/(T(:)*ZBOLTZ*1.0E6)  ! no conc of air (/cm3)

! copy from delso2_wet_xxx arrays as output from UKCA_CHEMISTRY_CTL
       DELSO2  (:) = RESHAPE(delso2_wet_h2o2(:,:,:),(/field_size3d/))
       DELSO2_2(:) = RESHAPE(delso2_wet_o3  (:,:,:),(/field_size3d/))
!
! set these to zero as will not be used in UKCA_AERO_STEP
      ZO3(:)=0.0     ! currently do wet ox separately in UM
      ZHO2(:)=0.0    ! currently do wet ox separately in UM
      ZH2O2(:)=0.0   ! currently do wet ox separately in UM
      LDAY(:)=0      ! currently do wet ox separately in UM
!
      IF (wet_levels == model_levels) THEN
        RH(:)=RESHAPE(rh3d(:,:,:),(/field_size3d/))
        S(:)=RESHAPE(q(:,:,:),(/field_size3d/))
        LWC(:)=RESHAPE(cloud_liq_wat(:,:,:),(/field_size3d/)) ! *CHECK UNITS IF USED
        CLWC(:)=RESHAPE(cloud_liq_wat(:,:,:),(/field_size3d/))
      ELSE
        RH(:)=0.0
        RH(1:row_length*rows*wet_levels)=reshape(rh3d(:,:,:),           &
                     (/row_length*rows*wet_levels/))
        S(:)=0.0
        S(1:row_length*rows*wet_levels)=reshape(q(:,:,:),               &
                     (/row_length*rows*wet_levels/))
        LWC(:)=0.0
        LWC(1:row_length*rows*wet_levels)=reshape(cloud_liq_wat(:,:,:), &
                     (/row_length*rows*wet_levels/))
        CLWC(1:row_length*rows*wet_levels)=reshape(cloud_liq_wat(:,:,:),&
                     (/row_length*rows*wet_levels/))
      ENDIF

      LOWCLOUD(:)=0.0
      VFAC(:)=0.0
      CLF(:) = 0.0
      CLF(1:row_length*rows*wet_levels)=RESHAPE(cloud_liq_frac(:,:,:), &
              (/row_length*rows*wet_levels/))
      LOWCLOUD(1:row_length*rows*wet_levels)=RESHAPE(cloud_frac(:,:,:), &
              (/row_length*rows*wet_levels/))
      WHERE (LOWCLOUD > 0.0)
        VFAC = 1.0 ! set to 1 so that VFAC*LOWCLOUD=cloud_frac
      ENDWHERE

      DO L=1,model_levels
       DLAT3D(:,:,L)=ASIN(sinlat(:,:))
       DLON3D(:,:,L)=true_longitude(:,:)
      ENDDO
      DLATARR(:)=RESHAPE(DLAT3D(:,:,:),(/field_size3d/))
      DLONARR(:)=RESHAPE(DLON3D(:,:,:),(/field_size3d/))

      HEIGHT (:)=RESHAPE(z_half_alllevs(:,:,:),(/field_size3d/))

      DELTA_Z(:)=RESHAPE(delta_r(:,:,:),(/field_size3d/))

      IF (VERBOSE >= 2) THEN
        WRITE(6,'(A35,2E12.3)') 'UKCA_AERO_CTL:MAX,MIN of delta_r= ',   &
                   maxval(delta_r(:,:,:)),minval(delta_r(:,:,:))
        WRITE(6,'(A35,2E12.3)') 'UKCA_AERO_CTL:MAX,MIN of delta_z= ',   &
                   maxval(delta_z(:)),minval(delta_z(:))
      END IF

      myproc=1 ! just set to dummy value for now

! GM added code here to only set cloud fraction > 0 if in low cloud
!    here low cloud is defined as being cloud with p>=680hPa
      MASK_NOTLC(:)=(PMID(:) < 680.0e2) ! PMID is in Pa
      WHERE(MASK_NOTLC(:))
       LOWCLOUD(:)=0.0
       LWC(:)=0.0
      ENDWHERE
!
      ALLOCATE(z(1:row_length,1:rows,1:model_levels))
      z(:,:,:)=SPREAD(u_s(:,:),DIM=3,NCOPIES=model_levels)
      USTR(:)=RESHAPE(z(:,:,:),(/field_size3d/))
      z(:,:,:)=SPREAD(u_10m(:,:),DIM=3,NCOPIES=model_levels)
      US10M(:)=RESHAPE(z(:,:,:),(/field_size3d/))
      z(:,:,:)=SPREAD(z0m(:,:),DIM=3,NCOPIES=model_levels)
      ZNOTG(:)=RESHAPE(z(:,:,:),(/field_size3d/))
      z(:,:,:)=SPREAD(sea_ice_frac(:,:),DIM=3,NCOPIES=model_levels)
      SEAICE(:)=RESHAPE(z(:,:,:),(/field_size3d/))
! fraction of land at surface
      z(:,:,:)=SPREAD(land_fraction(:,:),DIM=3,NCOPIES=model_levels)
      LAND_FRAC(:)=RESHAPE(z(:,:,:),(/field_size3d/))
      z(:,:,:)=SPREAD(zbl(:,:),DIM=3,NCOPIES=model_levels)
      HTPBLG(:)=RESHAPE(z(:,:,:),(/field_size3d/))
      DEALLOCATE(z)

      WHERE ( LAND_FRAC(1:row_length*rows) < 0.5 )   ! >50% sea at surface
        SURTP(1:field_size)=0.0
      ELSEWHERE
        SURTP(1:field_size)=1.0
      ENDWHERE

      WHERE ( LAND_FRAC(field_size+1:) < 0.5 )      ! >50% sea, above surface
        SURTP(field_size+1:)=2.0
      ELSEWHERE
        SURTP(field_size+1:)=3.0
      ENDWHERE

!---------------------------------------------------------------
! Put in section here to set land-surface category ILSCAT based on ZNOTG
! ILSCAT=1-9 based on 9 UM landsurf types but use existing approach
! to set 4 possible types according to roughness length

! Find out what category (water,forest,grass,desert)
! based on roughness length ZNOTG (desert not used at present)
! This should be updated in later version to read land type
! category directly. Desert not represented here.
!
! Now set to match 9 UM Land-use types (ILSCAT)      YR  ALPHA  CR
! 1=BL tree  (Zhang cat=avg of 2,4) [Evrgrn,Dec BL] 0.57  0.70 0.005
! 2=NL tree  (Zhang cat=avg of 1,3) [Evrgrn,Dec NL] 0.56  1.05 0.002
! 3=C3 grass (Zhang cat=6)          [grass        ] 0.54  1.20 0.002
! 4=C4 grass (Zhang cat=6)          [grass        ] 0.54  1.20 0.002
! 5=Shrub    (Zhang cat=10)         [shrub i. wood] 0.54  1.30 0.010
! 6=Urban    (Zhang cat=15)         [urban        ] 0.56  1.50 1.500
! 7=Water    (Zhang cat=13/14)      [inl wat/ocean] 0.50 100.0 0.000
! 8=Soil     (Zhang cat=8)          [desert       ] 0.54 50.00 0.000
! 9=Ice      (Zhang cat=12)         [ice cap/glac.] 0.54 50.00 0.000
! water/sea - z0<0.001m
!
! Previously used code below to set ILSCAT based on roughness length
! as either 1 = water/sea, 2= forest, 3 = all other land
!           4 = desert (not used), 5 = sea ice
! and used values from Zhang et al (2001) for YR (gamma in paper),
!                                             ALPHA
!                                             CR (A in paper)
! BELOW IS OLD TREATMENT                            YR   ALPHA    CR
! 1=water  (Zhang cat=13/14)[inland water/ocean  ] 0.50  100.0   0.000
! 2=forest (Zhang cat=1    )[evergreen needleleaf] 0.56    1.0   0.005
! 3=o land (Zhang cat=6/7  )[grass/crops         ] 0.54    1.2   0.002
! 4=desert (Zhang cat=8    )[desert              ] 0.54   50.0   0.000
! 5=seaice (Zhang cat=12   )[ice cap & glacier   ] 0.54   50.0   0.000

! Find out what category (water,forest,grass,desert)
! based on roughness length ZNOTG (desert not used at present)
! This should be updated in later version to read land type
! category directly. Desert not represented here.

      MASK1(:)=(ZNOTG(:) < 1.0E-3) ! water/sea
      WHERE(MASK1(:)) ILSCAT(:)=7

! forests - z0>0.1m
      MASK1(:)=(ZNOTG(:) > 1.0E-1) ! forest
      WHERE(MASK1(:)) ILSCAT(:)=1

! all other lands, grass 0.001<z0<0.1m
      MASK1(:)=((ZNOTG(:) >= 1.0E-3).AND.(ZNOTG(:) <= 1.0E-1)) ! grass
      WHERE(MASK1(:)) ILSCAT(:)=3

! If sea ice covers > 50% of sea surface, treat as sea ice
      MASK1(:)=(SEAICE(:) > 0.5) ! seaice
      WHERE(MASK1(:)) ILSCAT(:)=9

!---------------------------------------------------------------

! Derived quantities
! ==================
      AIRDM3(:)=AIRD(:)*1.0e6              ! no conc of air (/m3)
      RHOA(:)=PMID(:)/(T(:)*RA)
      VBA(:)=SQRT(8.0*ZBOLTZ*T(:)/(PPI*MA))
      TSQRT(:)=SQRT(T(:))
      DVISC(:)=1.83E-5*(416.16/(T(:)+120.0))*(SQRT(T(:)/296.16)**3)
      MFPA(:)=2.0*DVISC(:)/(RHOA(:)*VBA(:))
      SM(:)   =RESHAPE(mass(:,:,:),(/field_size3d/)) ! mass air (kg/box)
!
!-----------------------

! DO BL MIXING HERE with (for now) emissions, and resistance factor set to 
!   zero, TR_MIX is called with no halos and Ltimer set to false.

      em_field(:,:,:) = 0.0
      res_factor(:,:) = 0.0
      rho_aresist(:,:) = 0.0
      DO k=1,n_mode_tracers
        CALL TR_MIX(                                                    &
             bl_levels,alpha_cd(1:bl_levels),                           &
             rhokh_mix(:,:,2:), rho_aresist,                            &
             dtrdz_charney_grid,em_field(:,:,k), res_factor(:,:),       &
             kent, we_lim, t_frac, zrzi,                                &
             kent_dsc, we_lim_dsc, t_frac_dsc,                          &
             zrzi_dsc, ml_depth, zhsc, z_half,                          &
             mode_tracers(:,:,1:bl_levels,k),tr_flux(:,:,1:bl_levels,k),&
             surf_dep_flux(:,:,k)                                       &
             )

      ENDDO          ! end of loop over n_mode_tracers

! Gas-phase tracers required in aerosol code
! ==========================================
!  S0 array that is passed in to aerosol code uses MH2SO4, MSEC_ORG, etc
!     to index tracers (set in UKCA_SETUP_INDICES module procedures)

      S0(:,:)=0.0                 ! set gas phase tracer masses to zero

      IF(DRYOX_IN_AER == 0) THEN
       S0_DOT_CONDENSABLE(:,:)=0.0 ! condensable tracer tendencies->0
! DRYOX_IN_AER=0 -> update of condensables done in UKCA_CHEMISTRY_CTL 
      ENDIF ! if DRYOX_IN_AER=0

      IF(DRYOX_IN_AER == 1) THEN
! Initialise S0_DOT_CONDENSABLE to 0
       S0_DOT_CONDENSABLE(:,:)=0.0 ! condensable tracer tendencies->0
! DRYOX_IN_AER=1 -> update of condensables done in UKCA_AERO_STEP
! condensable tracer tendencies need to be set here to pass in as input

        S0_DOT_CONDENSABLE(:,MH2SO4)=                                   &
           RESHAPE(delso2_dry_oh (:,:,:),(/field_size3d/))              &
          /(AIRD(:)*DTC)
! .. delso2_dry_oh  is in units of molecules/cc/DTC
! .. need S0_DOT_CONDENSABLE to be in units of vmr/s
! .. so need to divide by (AIRD(:)*DTC)

       IF(MSEC_ORG > 0) THEN
        S0_DOT_CONDENSABLE(:,MSEC_ORG)=0.0
! for L_classSO2_inAer=T or F set SEC_ORG prodn--> 0 (done in UKCA)
       ENDIF
!
      ENDIF      ! dryox_in_aer
! Set H2O2, O3 and SO2 for aqueous oxidation
      IF (WETOX_IN_AER == 1) THEN
! .. set h2o2 when required
       IF (MH2O2F  > 0 .AND. MM_GAS(MH2O2F) > 1e-3 ) THEN
         S0(:,MH2O2F  )=SM(:)*(MM_DA/MM_GAS(MH2O2F ))*                  &
           RESHAPE(chemistry_tracers(:,:,:,n_h2o2),(/field_size3d/))
       ELSE
         cmessage=' H2O2 needs updating, but MH2O2F'//                  &
                    'or MM_GAS(MH2O2F) is wrong'
         write(6,*) cmessage,MH2O2F,MM_GAS(MH2O2F)
         errcode = 1
         CALL EREPORT('UKCA_AERO_CTL',errcode,cmessage)
       ENDIF

! .. set O3 
       IF (MOX > 0 .AND. MM_GAS(MOX) > 1e-3 ) THEN
         ZO3(:) = SM(:)*(MM_DA/MM_GAS(MOX))*                            &
           RESHAPE(chemistry_tracers(:,:,:,n_o3),(/field_size3d/))
       ELSE
         cmessage=' O3 needs updating, but MOX or MM_GAS(MOX) is wrong'
         write(6,*) cmessage,' MOX = ',MOX,MM_GAS(MOX)
         errcode = 1
         CALL EREPORT('UKCA_AERO_CTL',errcode,cmessage)
       ENDIF

! .. SO2 mmr in kg[SO2]/kg[dryair]
       IF (MSOTWO > 0 .AND. MM_GAS(MSOTWO) > 1e-3) THEN
         S0(:,MSOTWO)=SM(:)*(MM_DA/MM_GAS(MSOTWO))*                     &
           RESHAPE(chemistry_tracers(:,:,:,n_so2),(/field_size3d/))
       ELSE
         cmessage=' SO2 needs updating, but MSOTWO or MM_GAS(MSOTWO)'// &
                  'is wrong'
         write(6,*) cmessage,' MSOTWO = ',MSOTWO,MM_GAS(MSOTWO)
         errcode = 1
         CALL EREPORT('UKCA_AERO_CTL',errcode,cmessage)
       ENDIF

     ENDIF      ! wetox_in_aer=1

      IF (MH2SO4   > 0) THEN
        S0(:,MH2SO4  )=SM(:)*(MM_DA/MM_GAS(MH2SO4  ))*                  &
          RESHAPE(chemistry_tracers(:,:,:,n_h2so4),(/field_size3d/))
      ELSE
        errcode=1
        cmessage='MH2SO4 <= 0'
        write(6,*) cmessage,' MH2SO4 = ',MH2SO4
        CALL EREPORT('UKCA_AERO_CTL',errcode,cmessage)
      ENDIF

! .. Secondary Organic tracer mmr in kg[SEC_ORG]/kg[dryair]
      IF (MSEC_ORG > 0) THEN
        S0(:,MSEC_ORG)=SM(:)*(MM_DA/MM_GAS(MSEC_ORG))*                  &
          RESHAPE(chemistry_tracers(:,:,:,n_sec_org),(/field_size3d/))
      ELSE
         ! may not have Sec_Org tracer - e.g. StratAer
        errcode=-1
        cmessage='MSEC_ORG <= 0'
        write(6,*) cmessage,' MSEC_ORG = ',MSEC_ORG
        CALL EREPORT('UKCA_AERO_CTL',errcode,cmessage)
      ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dust:
! Section below just initialises dust emission inputs to constant
! values everywhere for now --- currently dust is switched off.

      SILT(:)=0.6
      SAND(:)=0.2
      CLAY(:)=0.2
      SWC(:)=0.005
      SNOWICE(:)=1.0  ! SNOWICE is the fraction which is *not* ice
      PSOURCE(:)=0.15 !        PSOURCE(:)=0.0
      LAI(:)=0.0
      TEX(:)=6        ! Silt sized aggregates

      NDUSTEMINBOX(:)=1
      MODE2_NUMBER(:,1)=2.42594e20 ! in number/1x1gridbox/day
      MODE3_NUMBER(:,1)=5.34868e20 ! in number/1x1gridbox/day
      MODE2_RADIUS(:,1)=0.215908   ! in microns
      MODE3_RADIUS(:,1)=0.588690   ! in microns
      MODE2_NUMBER(:,2:10)=0.0
      MODE3_NUMBER(:,2:10)=0.0
      MODE2_RADIUS(:,2:10)=0.215908 ! in microns
      MODE3_RADIUS(:,2:10)=0.588690 ! in microns

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Surface Emissions
! =================
      EMANSO2(:,:)=0.0
      EMVOLCONSO2(:)=0.0
      EMVOLEXPSO2(:)=0.0
      EMBIOMSO2(:)=0.0
      IF(PRIMSU_ON == 1) THEN      ! do primary S emissions
! .. multiplying by SURF below converts from kgSO2/m2/s to kgSO2/box/s
        EMANSO2(1:field_size,1)=RESHAPE(emissions(:,:,lso2emlo),        &
                            (/field_size/))*SURF(1:field_size)

! .. Add high-level primary SO4 emissions into level specified in UMUI
! .. multiplying by SURF below converts from kgSO2/m2/s to kgSO2/box/s
        i=row_length*rows*(SO2_high_level-1)+1
        j=row_length*rows*SO2_high_level
        EMANSO2(i:j,2)=RESHAPE(emissions(:,:,lso2emhi),                 &
                              (/field_size/))*SURF(i:j)

! .. Add volcanic SO4 emissions,
! ..  multiplying by SURF below converts from kgSO2/m2/s to kgSO2/box/s
        EMVOLCONSO2(:)=RESHAPE(SO2_volc_3D(:,:,:),                      &
                              (/field_size3d/))*SURF(:)

      ENDIF ! if PRIMSU_ON = 1

      EMC(:,:)=0.0
      EMCBM(:,:)=0.0
      IF (PRIMBCOC_ON == 1) THEN
! Set carbonaceous aerosol emissions:
! ..  emissions(:,lBCbf/lBCff/lOCbf/lOCff/BC_biom_3D/OC_biom_3D) are in kgC/m2/s
! ..  multiplying by SURF below converts from kgC/m2/s to kgC/box/s

        IF(L_bcoc_bf) THEN
         IF(PrintStatus >= PrStatus_Oper)                               &
           WRITE(6,'(A35)') 'Setting biofuel BC/OC emissions'
         EMC(1:field_size,1)=RESHAPE(emissions(:,:,lBCbf),              &
                    (/field_size/))*SURF(1:field_size)  ! 2D BC bf ems
         EMC(1:field_size,3)=RESHAPE(emissions(:,:,lOCbf),              &
                    (/field_size/))*SURF(1:field_size)  ! 2D OC bf ems
        ELSE
         IF(PrintStatus >= PrStatus_Oper)                               &
           WRITE(6,'(A35)') 'Not set biofuel BC/OC emissions'
        ENDIF
        IF(L_bcoc_ff) THEN
         IF(PrintStatus >= PrStatus_Oper)                               &
           WRITE(6,'(A35)') 'Setting fossil-fuel BC/OC emissions'
         EMC(1:field_size,2)=RESHAPE(emissions(:,:,lBCff),              &
                    (/field_size/))*SURF(1:field_size)  ! 2D BC ff ems
         EMC(1:field_size,4)=RESHAPE(emissions(:,:,lOCff),              &
                    (/field_size/))*SURF(1:field_size)  ! 2D OC ff ems
        ELSE
         IF(PrintStatus >= PrStatus_Oper)                               &
           WRITE(6,'(A35)') 'Not set fossil-fuel BC/OC emissions'
        ENDIF
        IF(L_bcoc_bm) THEN
         IF(PrintStatus >= PrStatus_Oper)                               &
           WRITE(6,'(A35)') 'Setting biomass burning BC/OC ems'
         EMCBM(:,1)=RESHAPE(BC_biom_3D(:,:,:),(/field_size3d/))*SURF(:)
         EMCBM(:,2)=RESHAPE(OC_biom_3D(:,:,:),(/field_size3d/))*        &
                    SURF(:)/emfactor
        ELSE
         IF(PrintStatus >= PrStatus_Oper)                               &
           WRITE(6,'(A35)') 'Not set biomass burning BC/OC ems'
        ENDIF
!
      ENDIF ! PRIMBCOC_ON=1

! Aerosol Tracers
! ===============
      MDTFIXFLAG(:,:) = 0.0
      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
        ITRA=II_ND(IMODE)
        tr_rs(:)=RESHAPE(mode_tracers(:,:,:,ITRA),(/field_size3d/))
        MASK_TRLT0(:)=(tr_rs(:) < 0.0)
        WHERE(MASK_TRLT0(:)) 
         tr_rs(:)=0.0
        ENDWHERE
! .. above sets tr_rs to zero if negative
        ND(:,IMODE)=tr_rs(:)*AIRD(:)
! .. above sets ND (particles per cc) from advected number-mixing-ratio
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          ITRA=II_MD(IMODE,ICP)
          tr_rs(:)=RESHAPE(mode_tracers(:,:,:,ITRA),(/field_size3d/))
          MASK_TRLT0(:)=(tr_rs(:) < 0.0)
          WHERE(MASK_TRLT0(:)) 
           tr_rs(:)=0.0
          ENDWHERE
! .. above sets tr_rs to zero if negative
          MASK_NDGTE(:)=(ND(:,IMODE) > NUM_EPS(IMODE))
          WHERE(MASK_NDGTE(:)) 
           MD(:,IMODE,ICP)=(MM_DA/MM(ICP))*AIRD(:)*tr_rs(:)/ND(:,IMODE)
          ELSEWHERE
           MD(:,IMODE,ICP)=MMID(IMODE)*MFRAC_0(IMODE,ICP)
          ENDWHERE
! .. above sets MD (molecules per particle) from advected mass-mixing-ratio
! .. note that only "trusts" advected values where ND>NUM_EPS
         ELSE
          MD(:,IMODE,ICP)=0.0
         ENDIF
        ENDDO ! loop over cpts
!
! Set total mass array MDT from sum over individual component MDs
        MDT(:,IMODE)=0.0
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
           MDT(:,IMODE)=MDT(:,IMODE)+MD(:,IMODE,ICP)
         ENDIF
        ENDDO
!
! below checks if MDT is coming out too low after advection (af BLMIX)
!
        MDTMIN=MLO(IMODE)*0.001 ! set equiv. to DPLIM0*0.1
!
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
! where MDT too low after advection set ND to zero and set default MD
          WHERE(MDT(:,IMODE) < MDTMIN)
           MD(:,IMODE,ICP)=MMID(IMODE)*MFRAC_0(IMODE,ICP)
          ENDWHERE
         ENDIF
        ENDDO

! below count occurrences (NBADMDT) and store % where
        NBADMDT(:,IMODE)=0
        MDTFIXFLAG(:,IMODE)=0.0
        WHERE(MDT(:,IMODE) < MDTMIN)
         NBADMDT(:,IMODE)=1
         MDTFIXFLAG(:,IMODE)=100.0
        ENDWHERE
!
! below set ND->0 where MDT toolo (& set MDT->MMID) & count occurrences
        WHERE(MDT(:,IMODE) < MDTMIN)
         ND (:,IMODE)=0.0
         MDT(:,IMODE)=MMID(IMODE)
! where MDT too low after advection (but +ve) set to MMID
        ENDWHERE
!
        IF(sum(NBADMDT(:,IMODE)) > 0 .AND. VERBOSE > 0) THEN
! below print out total occurrences if > 0
          WRITE(6,'(A55)') 'MDT<MDTMIN, ND=0: IMODE,MDTMIN,NBADMDT'//   &
                           ' (after bl mix)'
          WRITE(6,'(I6,E12.3,I12)') imode,mdtmin,SUM(nbadmdt(:,imode))
        ENDIF
!
       ELSE
        NBADMDT(:,IMODE)=0
        ND(:,IMODE)=0.0
        MDT(:,IMODE)=MMID(IMODE)
        DO ICP=1,NCP
         IF(COMPONENT(IMODE,ICP)) THEN
          MD(:,IMODE,ICP)=MMID(IMODE)*MFRAC_0(IMODE,ICP)
         ELSE
          MD(:,IMODE,ICP)=0.0
         ENDIF
        ENDDO
       ENDIF
      ENDDO ! loop over modes
!
      IF(VERBOSE >= 2) THEN
       do itra=1,nadvg
        write(6,*) 'S0 : ',itra,minval(s0(:,itra)),maxval(s0(:,itra))
       enddo
       do imode=1,nmodes
        if(mode(imode)) then
         write(6,*) 'ND : ',IMODE,minval(ND (:,IMODE)),                 &
                                  maxval(ND (:,IMODE))
         write(6,*) 'MDT: ',IMODE,minval(MDT(:,IMODE)),                 &
                                  maxval(MDT(:,IMODE))
         do icp=1,ncp
          if(component(imode,icp)) then
           write(6,*) 'MD : ',IMODE,minval(MD(:,IMODE,ICP))             &
                                   ,maxval(MD(:,IMODE,ICP))
          endif
         enddo
        endif
       enddo
      END IF  ! verbose >
!
! .. zero aerosol budget terms before calling UKCA_AERO_STEP
      BUD_AER_MAS(:,:)=0.0

! .. below is call to UKCA_AERO_STEP as at ukca_mode_v1_gm1.f90

      IF(firstcall .AND. VERBOSE > 0) THEN
       write(6,*) 'Values of input variables passed to mode:'
       write(6,*) 'PRIMSU_ON=',PRIMSU_ON
       write(6,*) 'PRIMBCOC_ON=',PRIMBCOC_ON
       write(6,*) 'RAINOUT_ON=',RAINOUT_ON
       write(6,*) 'IMSCAV_ON=',IMSCAV_ON
       write(6,*) 'WETOX_ON=',WETOX_ON
       write(6,*) 'DDEPAER_ON=',DDEPAER_ON
       write(6,*) 'SEDI_ON=',SEDI_ON
       write(6,*) 'I_MODE_SS_SCHEME=',I_MODE_SS_SCHEME
       write(6,*) 'DRYOX_IN AER=',DRYOX_IN_AER
       write(6,*) 'WETOX_IN AER=',WETOX_IN_AER
       write(6,*) 'COND_ON=',COND_ON
       write(6,*) 'NUCL_ON=',NUCL_ON
       write(6,*) 'BLN_ON=',BLN_ON
       write(6,*) 'COAG_ON=',COAG_ON
       write(6,*) 'ICOAG=',ICOAG
       write(6,*) 'IMERGE=',IMERGE
       write(6,*) 'IFUCHS=',IFUCHS
       write(6,*) 'IWVOLMETHOD=',IWVOLMETHOD
       write(6,*) 'IDCMFP=',IDCMFP
       write(6,*) 'ICONDIAM=',ICONDIAM
       write(6,*) 'IBLN=',IBLN
       write(6,*) 'IACTMETHOD=',IACTMETHOD
       write(6,*) 'I_NUC_METHOD=',I_NUC_METHOD
       write(6,*) 'IDDEPAER=',IDDEPAER
       write(6,*) 'INUCSCAV=',INUCSCAV
       write(6,*) 'VERBOSE=',VERBOSE
       write(6,*) 'CHECKMD_ND=',CHECKMD_ND
       write(6,*) 'INTRAOFF=',INTRAOFF
       write(6,*) 'INTEROFF=',INTEROFF
       write(6,*) 'IDUSTEMS=',IDUSTEMS
      ENDIF
!
! DEPENDS ON: ukca_aero_step
      CALL UKCA_AERO_STEP(NBOX,                                         &
       ND,MDT,MD,MDWAT,S0,DRYDP,WETDP,RHOPAR,DVOL,WVOL,SM,              &
       AIRD,AIRDM3,RHOA,MFPA,DVISC,T,TSQRT,RH,S,PMID,PUPPER,PLOWER,     &
       EMC,EMCBM,ZO3,ZHO2,ZH2O2,USTR,US10M,ZNOTG,DELTA_Z,               &
       SURTP,LAND_FRAC,SURF,SEAICE,                                     &
       CRAING,DRAING,CRAING_UP,DRAING_UP,FCONV_CONV,LOWCLOUD,VFAC,CLF,  &
       EMANSO2,EMVOLCONSO2,EMVOLEXPSO2,EMBIOMSO2,ISO2EMS,               &
       SILT,CLAY,SAND,SWC,SNOWICE,                                      &
       PSOURCE,LAI,TEX,MODE2_RADIUS,MODE3_RADIUS,                       &
       MODE2_NUMBER,MODE3_NUMBER,NDUSTEMINBOX,                          &
       DTC,DTM,DTZ,NMTS,NZTS,LDAY,RMOIS,RJOUR,ACT,BUD_AER_MAS,          &
       RAINOUT_ON,                                                      &
       IMSCAV_ON,WETOX_ON,DDEPAER_ON,SEDI_ON,ISO2WETOXBYO3,             &
       DRYOX_IN_AER,WETOX_IN_AER,DELSO2,DELSO2_2,                       &
       COND_ON,NUCL_ON,COAG_ON,BLN_ON,ICOAG,IMERGE,IFUCHS,IWVOLMETHOD,  &
       IDCMFP,ICONDIAM,IBLN,I_NUC_METHOD,                               &
       IACTMETHOD,IDDEPAER,INUCSCAV,VERBOSE,CHECKMD_ND,INTRAOFF,        &
       INTEROFF,IDUSTEMS,S0_DOT_CONDENSABLE,LWC,CLWC,PVOL,PVOL_WAT,     &
       JLABOVE,ILSCAT,N_MERGE_1D,DLONARR,DLATARR,HEIGHT,HTPBLG,myproc)

!
! Update tracers
! ==============
      IF (MH2O2F > 0 .AND. WETOX_IN_AER > 0) THEN
! .. update gas phase H2O2 mmr following SO2 aqueous phase oxidation
        chemistry_tracers(:,:,:,n_h2o2)=RESHAPE((S0(:,MH2O2F  )/SM(:)), &
        (/row_length,rows,model_levels/))*(MM_GAS(MH2O2F)/MM_DA)
      ENDIF
      IF (MSOTWO > 0 .AND. WETOX_IN_AER > 0) THEN
! .. update gas phase SO2 mmr following aqueous phase oxidation
        chemistry_tracers(:,:,:,n_so2)=RESHAPE((S0(:,MSOTWO)/SM(:)),    &
        (/row_length,rows,model_levels/))*(MM_GAS(MSOTWO)/MM_DA)
      ENDIF
      IF (MH2SO4   > 0) THEN
! .. update gas phase H2SO4 mmr following H2SO4 condensation/nucleation
        chemistry_tracers(:,:,:,n_h2so4)=RESHAPE((S0(:,MH2SO4  )/SM(:)),&
        (/row_length,rows,model_levels/))*(MM_GAS(MH2SO4  )/MM_DA)
      ENDIF
      IF (MSEC_ORG > 0) THEN
! .. update gas phase SEC_ORG mmr following condensation/nucleation
        chemistry_tracers(:,:,:,n_sec_org)=                             &
             RESHAPE((S0(:,MSEC_ORG)/SM(:)),                            &
             (/row_length,rows,model_levels/))*(MM_GAS(MSEC_ORG)/MM_DA)
      ENDIF

! Set H2O2, O3 and SO2 for aqueous oxidation
      IF (WETOX_IN_AER == 1) THEN
! .. set h2o2 from h2o2_tracer when required
       IF (MH2O2F  > 0 .AND. MM_GAS(MH2O2F) > 1e-3 ) THEN
         S0(:,MH2O2F  )=SM(:)*(MM_DA/MM_GAS(MH2O2F ))*                  &
           RESHAPE(chemistry_tracers(:,:,:,n_h2o2),(/field_size3d/))
       ELSE
         cmessage=' H2O2 needs updating, but MH2O2F'//                  &
                    'or MM_GAS(MH2O2F) is wrong'
         write(6,*) cmessage,MH2O2F,MM_GAS(MH2O2F)
         errcode = 1
         CALL EREPORT('UKCA_AERO_CTL',errcode,cmessage)
       END IF

! .. set O3 from O3_tracer
       IF (MOX > 0 .AND. MM_GAS(MOX) > 1e-3 ) THEN
         ZO3(:) = SM(:)*(MM_DA/MM_GAS(MOX))*                            &
           RESHAPE(chemistry_tracers(:,:,:,n_o3),(/field_size3d/))
       ELSE
         cmessage=' O3 needs updating, but MOX or MM_GAS(MOX) is wrong'
         write(6,*) cmessage,' MOX = ',MOX,MM_GAS(MOX)
         errcode = 1
         CALL EREPORT('UKCA_AERO_CTL',errcode,cmessage)
       END IF

       IF (MSOTWO > 0 .AND. MM_GAS(MSOTWO) > 1e-3) THEN
         S0(:,MSOTWO)=SM(:)*(MM_DA/MM_GAS(MSOTWO))*                     &
           RESHAPE(chemistry_tracers(:,:,:,n_so2),(/field_size3d/))
       ELSE
         cmessage=' SO2 needs updating, but MSOTWO or MM_GAS(MSOTWO)'// &
                  'is wrong'
         write(6,*) cmessage,' MSOTWO = ',MSOTWO,MM_GAS(MSOTWO)
         errcode = 1
         CALL EREPORT('UKCA_AERO_CTL',errcode,cmessage)
       ENDIF
      END IF      ! wetox_in_aer

      DO IMODE=1,NMODES
        IF(MODE(IMODE)) THEN
         ITRA=II_ND(IMODE)
! .. update aerosol no. conc. following aerosol microphysics
         mode_tracers(:,:,:,ITRA)=RESHAPE((ND(:,IMODE)/AIRD(:)),        &
                      (/row_length,rows,model_levels/))
         DO ICP=1,NCP
          IF(COMPONENT(IMODE,ICP)) THEN
           ITRA=II_MD(IMODE,ICP)
           WHERE (ND(:,IMODE) <= NUM_EPS(IMODE))
             MD(:,IMODE,ICP)=MMID(IMODE)*MFRAC_0(IMODE,ICP)
           ENDWHERE
! .. update aerosol mmr following aerosol microphysics
           mode_tracers(:,:,:,ITRA)=(MM(ICP)/MM_DA)*                    &
              RESHAPE((MD(:,IMODE,ICP)*ND(:,IMODE)/AIRD(:)),            &
                      (/row_length,rows,model_levels/))
          ENDIF
         ENDDO ! loop over cpts

         N_MERGE_3D(:,:,:,IMODE)=RESHAPE(N_MERGE_1D(:,IMODE),           &
                 (/row_length,rows,model_levels/))
         VCONC(:,IMODE)=WVOL(:,IMODE)*ND(:,IMODE)
         Y(IMODE)=EXP(2.0*LOG(SIGMAG(IMODE))*LOG(SIGMAG(IMODE)))
! .. multiply by 1e4 to give units of cm^2 / cm^3
         SAREA(:,IMODE)=PPI*1E4*ND(:,IMODE)*(WETDP(:,IMODE)**2)*Y(IMODE)

        ENDIF     ! mode(imode)
      ENDDO       ! imode=1,nmodes

!--------------------------------------------------------
! below sets CN,CCN,CDN diagnostics

! intialise to zero
      CN_3NM(:)=0.0
      CCN_1 (:)=0.0
      CCN_2 (:)=0.0
      CCN_3 (:)=0.0

      DO IMODE=1,NMODES
       IF(MODE(IMODE)) THEN
!
        DP0=3.0e-9
        ERF_ARG(:)=LOG(DP0/DRYDP(:,IMODE))/SQRT(2.0)/LOG(SIGMAG(IMODE))
        ERFTERM(:)=0.5*ND(:,IMODE)*(1.0-ERF(ERF_ARG(:)))
        CN_3NM(:)=CN_3NM(:)+ERFTERM(:)
!
        IF(MODESOL(IMODE).EQ.1) THEN
!
! for CCN_1 take CCN for particles > ACT dry radius
         DP0=2.0*ACT
         ERF_ARG(:)=LOG(DP0/DRYDP(:,IMODE))/SQRT(2.0)/LOG(SIGMAG(IMODE))
         ERFTERM(:)=0.5*ND(:,IMODE)*(1.0-ERF(ERF_ARG(:)))
         CCN_1(:)=CCN_1(:)+ERFTERM(:)
!
! for CCN_2 take CCN for particles > 25nm dry radius
         DP0=50.0e-9
         ERF_ARG(:)=LOG(DP0/DRYDP(:,IMODE))/SQRT(2.0)/LOG(SIGMAG(IMODE))
         ERFTERM(:)=0.5*ND(:,IMODE)*(1.0-ERF(ERF_ARG(:)))
         CCN_2(:)=CCN_2(:)+ERFTERM(:)
!
! for CCN_3 take CCN for particles > 35nm dry radius
         DP0=70.0e-9
         ERF_ARG(:)=LOG(DP0/DRYDP(:,IMODE))/SQRT(2.0)/LOG(SIGMAG(IMODE))
         ERFTERM(:)=0.5*ND(:,IMODE)*(1.0-ERF(ERF_ARG(:)))
         CCN_3(:)=CCN_3(:)+ERFTERM(:)
!
        ENDIF
       ENDIF
      ENDDO
!
      CDN(:)=375.0*(1.0-EXP(-0.0025*CCN_1(:)))
!
      WHERE(CDN(:) < CDNMIN)
       CDN(:)=CDNMIN
      ENDWHERE
!
!----------------------------------------------------
!
      IF (VERBOSE > 1) THEN

        write(6,*) 'AFTER CALL TO UKCA_AERO_STEP'
        do itra=1,nadvg
         write(6,*) 'S0 : ',itra,minval(s0(:,itra)),maxval(s0(:,itra))
        enddo
        do imode=1,nmodes
         if(mode(imode)) then
          write(6,*) 'ND : ',IMODE,minval(ND (:,IMODE)),                &
                                   maxval(ND (:,IMODE))
          write(6,*) 'MDT: ',IMODE,minval(MDT(:,IMODE)),                &
                                   maxval(MDT(:,IMODE))
          do icp=1,ncp
           if(component(imode,icp)) then
            write(6,*) 'MD : ',IMODE,minval(MD(:,IMODE,ICP))            &
                                    ,maxval(MD(:,IMODE,ICP))
           endif
          enddo
         endif
        enddo

! DEPENDS ON: ukca_calc_drydiam
        CALL UKCA_CALC_DRYDIAM(NBOX,ND,MD,MDT,DRYDP,DVOL,VERBOSE)

        do imode=1,nmodes
         if(mode(imode)) then
          write(6,*) 'DRYDP',IMODE,minval(DRYDP(:,IMODE)),              &
                                   maxval(DRYDP(:,IMODE))
         endif
        enddo

      ENDIF ! if PrintStatus > PrStatus_oper

! Calculate heterogenous rate coeffs for tropospheric chemistry 
      IF (L_UKCA_trophet) THEN 
! DEPENDS ON: ukca_trop_hetchem 
        CALL UKCA_TROP_HETCHEM(nbox, nhet_std_trop, t, rh, aird, pvol,  &
                               wetdp, sarea, het_rates) 
      ENDIF

! Copy MODE diagnostics from BUD_AER_MAS to BUD_AER_MAS1D
! and convert all budget variables to be in kg/DTC

      FAC(:)=SM(:)/AIRD(:)
! FAC converts aerosol mass fluxes from kg(dryair)/box/tstep to moles/gridbox/s
      DO JV=1,NBUDAER
! NMASCLPR variables are in kg(dryair)/DTC, others in molecules/cc/DTC
        LOGIC1=(JV.LT.NMASCLPRSUAITSOL1)
        LOGIC2=(JV.GT.NMASCLPRSUCORSOL2)
        LOGIC=LOGIC1.OR.LOGIC2
        IF(LOGIC) THEN
! when passed out from AERO_STEP, BUD_AER_MAS is in molecules/cc/DTC
         BUD_AER_MAS(:,JV)=BUD_AER_MAS(:,JV)*FAC(:)/MM_DA/DTC ! moles/s
        ELSE
! when passed out from AERO_STEP, BUD_AER_MAS (NMASCLPR) is in kg(dryair)/box/DTC
         BUD_AER_MAS(:,JV)=BUD_AER_MAS(:,JV)/MM_DA/DTC  ! moles/s
        ENDIF
      ENDDO
!
! Write 3_D diagnostics to mode_diags array
!  N.B. L_ukca_mode_diags is set whenever a STASH request for a relevant
!  item is found, and the ukcaD1codes(N)%item is then set, otherwise it is IMDI

      IF (L_ukca_mode_diags) THEN  ! fill 3D array
       k=0
       DO N=1,nukca_D1items
        IF (ukcaD1codes(N)%section == MODE_diag_sect .AND.              &
            ukcaD1codes(N)%item >= item1_mode_diags .AND.               &
            ukcaD1codes(N)%item <= item1_mode_diags+                    &
                                   nmax_mode_diags-1 .AND.              &
            ukcaD1codes(N)%item /= IMDI) THEN
         k=k+1

! number for user psm  mode_fluxdiagsv6.6_gm3_SUSSBCOC_reduced_gm1
! item1_mode_diags=201
! prim SU -- 201-203        ! Note that primary emission diagnostics
! prim SS -- 204-205        ! are now in routine UKCA_MODE_EMS_UM.
! prim BC -- 206-207        !                "
! prim OC -- 208-209        !                "
! prim DU -- 210-213        !                "
! ddep SU -- 214-217
! ddep SS -- 218-219
! ddep BC -- 220-223
! ddep OC -- 224-228
! ddep SO -- 229-232
! ddep DU -- 233-236
! nusc SU -- 234-240
! nusc SS -- 241-242
! nusc BC -- 243-246
! nusc OC -- 247-251
! nusc SO -- 252-256
! nusc DU -- 257-260
! imsc SU -- 261-264
! imsc SS -- 265-266
! imsc BC -- 267-270
! imsc OC -- 271-275
! imsc SO -- 276-279
! imsc DU -- 280-283
! clpr SU -- 284-289 (this is wet oxidation of SO2)
! proc SU,BC,OC,SO -- 290-293 (this is processing Aitsol-->accsol)
! cond SU -- 294-300
! cond OC -- 301-307
! cond SO -- 308-314
! htox SU -- 315-318
! nucl SU -- 319
! coag SU,SS,BC,OC,SO,DU -- 320-370
! aged SU,SS,BC,OC,SO,DU -- 371-374
! merg SU,SS,BC,OC,SO,DU -- 375-387
! drydp1-7 - 401-407
! wetdp1-4 - 408-411
! mdwat1-4 - 412-415
! sarea1-7 - 416-422
! vconc1-7 - 423-429
! rhop 1-7 - 430-436
! cnccncdn - 437-441
! pvol,wat - 442-468
!
         IF(firstcall) THEN
          write(6,*) 'About to set mode diagnostics for N,k,item=',     &
                  N,k,ukcaD1codes(N)%item
         ENDIF
!
         SELECT CASE(UkcaD1codes(N)%item)
          CASE(item1_mode_diags:item1_mode_diags+12)
!           Do nothing, primary emissions not handled here now
          CASE(item1_mode_diags+13)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPSUNUCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+14)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPSUAITSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+15)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPSUACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+16)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPSUCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+17)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPSSACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+18)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPSSCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+19)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPBCAITSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+20)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPBCACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+21)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPBCCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+22)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPBCAITINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+23)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPOCNUCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+24)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPOCAITSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+25)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPOCACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+26)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPOCCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+27)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPOCAITINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+28)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPSONUCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+29)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPSOAITSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+30)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPSOACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+31)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPSOCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+32)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPDUACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+33)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPDUCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+34)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPDUACCINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+35)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASDDEPDUCORINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+36)
           mode_diags(:,:,:,k)=RESHAPE(MDTFIXFLAG(:,1),                 &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+37)
           mode_diags(:,:,:,k)=RESHAPE(MDTFIXFLAG(:,2),                 &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+38)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCSUACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+39)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCSUCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+40)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCSSACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+41)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCSSCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+42)
           mode_diags(:,:,:,k)=RESHAPE(MDTFIXFLAG(:,3),                 &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+43)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCBCACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+44)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCBCCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+45)
           mode_diags(:,:,:,k)=RESHAPE(MDTFIXFLAG(:,4),                 &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+46)
           mode_diags(:,:,:,k)=RESHAPE(MDTFIXFLAG(:,5),                 &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+47)
           mode_diags(:,:,:,k)=RESHAPE(MDTFIXFLAG(:,6),                 &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+48)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCOCACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+49)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCOCCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+50)
           mode_diags(:,:,:,k)=RESHAPE(MDTFIXFLAG(:,7),                 &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+51)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCSONUCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+52)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCSOAITSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+53)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCSOACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+54)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCSOCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+55)
           mode_diags(:,:,:,k)=0.0
!! .. there is no MASNUSCSOAITINS --- erroneously included in 
!! .. UKCA_mode stash section so set it to zero here
          CASE(item1_mode_diags+56)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCDUACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+57)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCDUCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+58)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCDUACCINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+59)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUSCDUCORINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+60)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCSUNUCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+61)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCSUAITSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+62)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCSUACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+63)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCSUCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+64)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCSSACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+65)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCSSCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+66)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCBCAITSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+67)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCBCACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+68)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCBCCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+69)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCBCAITINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+70)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCOCNUCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+71)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCOCAITSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+72)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCOCACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+73)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCOCCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+74)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCOCAITINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+75)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCSONUCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+76)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCSOAITSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+77)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCSOACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+78)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCSOCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+79)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCDUACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+80)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCDUCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+81)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCDUACCINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+82)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASIMSCDUCORINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+83)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCLPRSUAITSOL1),&
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+84)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCLPRSUACCSOL1),&
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+85)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCLPRSUCORSOL1),&
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+86)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCLPRSUAITSOL2),&
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+87)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCLPRSUACCSOL2),&
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+88)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCLPRSUCORSOL2),&
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+89)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASPROCSUINTR23), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+90)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASPROCBCINTR23), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+91)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASPROCOCINTR23), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+92)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASPROCSOINTR23), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+93)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSUNUCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+94)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSUAITSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+95)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSUACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+96)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSUCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+97)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSUAITINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+98)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSUACCINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+99)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSUCORINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+100)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDOCNUCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+101)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDOCAITSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+102)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDOCACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+103)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDOCCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+104)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDOCAITINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+105)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDOCACCINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+106)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDOCCORINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+107)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSONUCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+108)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSOAITSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+109)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSOACCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+110)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSOCORSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+111)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSOAITINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+112)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSOACCINS), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+113)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCONDSOCORINS), &
                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+114)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASHTOXSUACCSOL), &
!!                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+115)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASHTOXSUCORSOL), &
!!                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+116)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASHTOXSUACCINS), &
!!                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+117)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASHTOXSUCORINS), &
!!                                 (/row_length,rows,model_levels/))
!!
!! Code for heterogeneous oxidation of SO2 --> SO4 on dust not yet in
!!
          CASE(item1_mode_diags+114)
           mode_diags(:,:,:,k)=0.0
          CASE(item1_mode_diags+115)
           mode_diags(:,:,:,k)=0.0
          CASE(item1_mode_diags+116)
           mode_diags(:,:,:,k)=0.0
          CASE(item1_mode_diags+117)
           mode_diags(:,:,:,k)=0.0

          CASE(item1_mode_diags+118)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASNUCLSUNUCSOL), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+119)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSUINTR12), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+120)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSUINTR13), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+121)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSUINTR14), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+122)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSUINTR15), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+123)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSUINTR16), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+124)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSUINTR17), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+125)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR12), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+126)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR13), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+127)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR14), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+128)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR15), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+129)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR16), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+130)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR17), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+131)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSOINTR12), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+132)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSOINTR13), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+133)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSOINTR14), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+134)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSOINTR15), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+135)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSOINTR16), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+136)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSOINTR17), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+137)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSUINTR23), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+138)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSUINTR24), &
                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+139)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSUINTR26), &
!!                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+140)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSUINTR27), &
!!                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+139)
           mode_diags(:,:,:,k)=0.0
          CASE(item1_mode_diags+140)
           mode_diags(:,:,:,k)=0.0
!! .. above fluxes not included (only have mode 1 coag to ins modes)
          CASE(item1_mode_diags+141)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGBCINTR23), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+142)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGBCINTR24), &
                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+143)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGBCINTR26), &
!!                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+144)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGBCINTR27), &
!!                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+143)
           mode_diags(:,:,:,k)=0.0
          CASE(item1_mode_diags+144)
           mode_diags(:,:,:,k)=0.0
!! .. above fluxes not included (only have mode 1 coag to ins modes)
          CASE(item1_mode_diags+145)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR23), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+146)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR24), &
                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+147)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR26), &
!!                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+148)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR27), &
!!                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+147)
           mode_diags(:,:,:,k)=0.0
          CASE(item1_mode_diags+148)
           mode_diags(:,:,:,k)=0.0
!! .. above fluxes not included (only have mode 1 coag to ins modes)
          CASE(item1_mode_diags+149)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSOINTR23), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+150)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSOINTR24), &
                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+151)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSOINTR26), &
!!                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+152)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSOINTR27), &
!!                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+151)
           mode_diags(:,:,:,k)=0.0
          CASE(item1_mode_diags+152)
           mode_diags(:,:,:,k)=0.0
!! .. above fluxes not included (only have mode 1 coag to ins modes)
          CASE(item1_mode_diags+153)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSUINTR34), &
                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+154)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSUINTR37), &
!!                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+154)
           mode_diags(:,:,:,k)=0.0
!! .. above fluxes not included (only have mode 1 coag to ins modes)
          CASE(item1_mode_diags+155)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGBCINTR34), &
                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+156)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGBCINTR37), &
!!                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+156)
           mode_diags(:,:,:,k)=0.0
!! .. above fluxes not included (only have mode 1 coag to ins modes)
          CASE(item1_mode_diags+157)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR34), &
                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+158)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR37), &
!!                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+158)
           mode_diags(:,:,:,k)=0.0
!! .. above fluxes not included (only have mode 1 coag to ins modes)
          CASE(item1_mode_diags+159)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSSINTR34), &
                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+160)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSSINTR37), &
!!                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+160)
           mode_diags(:,:,:,k)=0.0
!! .. above fluxes not included (only have mode 1 coag to ins modes)
          CASE(item1_mode_diags+161)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSOINTR34), &
                                 (/row_length,rows,model_levels/))
!!          CASE(item1_mode_diags+162)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGSOINTR37), &
!!                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+162)
           mode_diags(:,:,:,k)=0.0
!! .. above fluxes not included (only have mode 1 coag to ins modes)
          CASE(item1_mode_diags+163)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGDUINTR34), &
                                 (/row_length,rows,model_levels/))
!!  CASE(item1_mode_diags+164)
!!           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGDUINTR37), &
!!                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+164)
           mode_diags(:,:,:,k)=0.0
!! .. above fluxes not included (only have mode 1 coag to ins modes)
          CASE(item1_mode_diags+165)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGBCINTR53), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+166)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR53), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+167)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGBCINTR54), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+168)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGOCINTR54), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+169)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASCOAGDUINTR64), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+170)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASAGEDSUINTR52), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+171)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASAGEDBCINTR52), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+172)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASAGEDOCINTR52), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+173)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASAGEDSOINTR52), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+174)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASMERGSUINTR12), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+175)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASMERGOCINTR12), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+176)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASMERGSOINTR12), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+177)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASMERGSUINTR23), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+178)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASMERGBCINTR23), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+179)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASMERGOCINTR23), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+180)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASMERGSOINTR23), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+181)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASMERGSUINTR34), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+182)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASMERGBCINTR34), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+183)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASMERGOCINTR34), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+184)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASMERGSSINTR34), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+185)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASMERGDUINTR34), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+186)
           mode_diags(:,:,:,k)=RESHAPE(BUD_AER_MAS(:,NMASMERGSOINTR34), &
                                 (/row_length,rows,model_levels/))
          CASE(item1_mode_diags+200)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(DRYDP(:,1),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+201)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(DRYDP(:,2),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+202)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(DRYDP(:,3),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+203)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(DRYDP(:,4),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+204)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(DRYDP(:,5),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+205)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(DRYDP(:,6),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+206)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(DRYDP(:,7),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+207)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(WETDP(:,1),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+208)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(WETDP(:,2),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+209)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(WETDP(:,3),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+210)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(WETDP(:,4),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+211)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(MDWAT(:,1)/AVC,(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+212)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(MDWAT(:,2)/AVC,(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+213)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(MDWAT(:,3)/AVC,(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+214)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(MDWAT(:,4)/AVC,(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+215)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(SAREA(:,1),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+216)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(SAREA(:,2),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+217)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(SAREA(:,3),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+218)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(SAREA(:,4),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+219)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(SAREA(:,5),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+220)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(SAREA(:,6),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+221)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(SAREA(:,7),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+222)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(VCONC(:,1),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+223)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(VCONC(:,2),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+224)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(VCONC(:,3),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+225)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(VCONC(:,4),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+226)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(VCONC(:,5),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+227)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(VCONC(:,6),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+228)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(VCONC(:,7),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+229)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(RHOPAR(:,1),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+230)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(RHOPAR(:,2),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+231)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(RHOPAR(:,3),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+232)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(RHOPAR(:,4),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+233)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(RHOPAR(:,5),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+234)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(RHOPAR(:,6),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+235)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(RHOPAR(:,7),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+236)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(CN_3NM(:),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+237)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(CCN_1(:),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+238)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(CCN_2(:),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+239)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(CCN_3(:),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+240)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(CDN(:),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+241)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,1,1),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+242)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,1,3),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+243)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,1,6),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+244)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL_WAT(:,1),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+245)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,2,1),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+246)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,2,2),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+247)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,2,3),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+248)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,2,6),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+249)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL_WAT(:,2),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+250)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,3,1),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+251)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,3,2),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+252)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,3,3),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+253)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,3,4),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+254)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,3,5),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+255)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,3,6),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+256)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL_WAT(:,3),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+257)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,4,1),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+258)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,4,2),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+259)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,4,3),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+260)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,4,4),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+261)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,4,5),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+262)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,4,6),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+263)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL_WAT(:,4),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+264)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,5,2),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+265)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,5,3),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+266)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,6,5),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+267)
           mode_diags(:,:,:,k)=                                         &
              RESHAPE(PVOL(:,7,5),(/row_length,rows,model_levels/))
          CASE(item1_mode_diags+268:item1_mode_diags+284)
! Do nothing - these are used by ukca_activate
          CASE DEFAULT
           cmessage=' Item not found in CASE statement'
          CALL EREPORT('UKCA_AERO_CTL',UkcaD1codes(N)%item,cmessage)
         END SELECT
         IF (PrintStatus >= PrStatus_Diag .AND. l==1) THEN
          WRITE(6,'(A16,3i4,3e14.4)') 'UKCA_MODE diag: ',               &
                        k,N,ukcaD1codes(N)%item,                        &
                        sum(mode_diags(:,:,:,k)),                       &
                     maxval(mode_diags(:,:,:,k)),                       &
                     minval(mode_diags(:,:,:,k))
         ENDIF
        ENDIF    ! section == MODE_diag_sec etc
       ENDDO      ! nukca_D1items
      ENDIF        ! L_UKCA_MODE_diags

      IF (.NOT. L_UKCA_ARG_ACT) THEN
        IF (L_UKCA_AIE1 .OR. L_UKCA_AIE2) THEN
! set ukca_cdnc prognostic without activation scheme
! change units from cm^-3 to m^-3
           chem_diag_cdnc(:,:,:) = 1.0E+6*                                &
                RESHAPE(CDN(:),(/row_length,rows,model_levels/))
        END IF
      END IF
! Calculate aerosol surface area for use in UKCA chemistry for heterogeneous reactions
      aerosol_surface_area(:,:,:) = 0.0
      DO imode=1,4 ! only take the first 4 (soluble) modes
         aerosol_surface_area(:,:,:) = aerosol_surface_area(:,:,:) + &
              1.0E4*RESHAPE(SAREA(:,imode),&
              (/row_length,rows,model_levels/))
         ! convert to cm^2/cm^3 from m^2/cm^3
      END DO

      IF (PrintStatus >= PrStatus_Diag) THEN

       WRITE(6,*) ' Tracers at end of UKCA_MODE:'
       DO i=1,2           !model_levels
        DO j=1,n_chemistry_tracers
        WRITE(6,*) 'Level: ',i,' Tracer: ',j
        WRITE(6,*) 'chemistry_tracers:',                                & 
                minval(chemistry_tracers(:,:,i,j)),                     &
                maxval(chemistry_tracers(:,:,i,j)),                     &
                sum(chemistry_tracers(:,:,i,j))/                        &
                real(size(chemistry_tracers(:,:,i,j)))
        ENDDO
        DO j=1,n_mode_tracers
         IF (mode_tracer_debug(j)) then
          WRITE(6,*) 'Level: ',i,' Tracer: ',j,mode_tracer_names(j)
          WRITE(6,*) 'mode_tracers: ',minval(mode_tracers(:,:,i,j)),    &
                 maxval(mode_tracers(:,:,i,j)),                         &
                 sum(mode_tracers(:,:,i,j))/                            &
                 real(size(mode_tracers(:,:,i,j)))
         ENDIF
        ENDDO
        WRITE(6,*) 'Number of merges for Level: ',i
        DO j=1,nmodes
         WRITE(6,*) j,sum(N_MERGE_3D(:,:,i,j))
        ENDDO
       ENDDO      ! i

      write(6,'(1a24,1i9,1e12.3)')                                      &
           'Total Number of merges=:',                                  &
            sum(N_MERGE_3D),real(sum(N_MERGE_3D))/real(size(N_MERGE_3D))
      DO j=1,nmodes
       write(6,'(2i9,1e12.3)') j,                                       &
            sum(N_MERGE_3D(:,:,:,j)),                                   &
       real(sum(N_MERGE_3D(:,:,:,j)))/real(size(N_MERGE_3D(:,:,:,j)))
      ENDDO
      ENDIF      ! PrintStatus

      IF(firstcall) firstcall=.FALSE.

      IF (lhook) CALL dr_hook('UKCA_AERO_CTL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_AERO_CTL
