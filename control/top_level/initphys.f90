! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine INITPHYS
!
! Purpose:
!   Calls the routines required to read in spectral files for the 
!   radiation scheme.
!
! Documentation:
!   UM documentation paper no 23
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
!- ---------------------------------------------------------------------
      SUBROUTINE INITPHYS(ICODE,CMESSAGE)

      USE rad_input_mod, ONLY:                                          &
          l_climat_aerosol, l_use_dust, l_use_arcldust,                 &
          l_use_sulpc_direct, l_use_arclsulp,                           &
          l_use_soot_direct, l_use_arclblck,                            &
          l_use_bmass_direct, l_use_arclbiom,                           &
          l_use_seasalt_direct, l_use_arclsslt,                         &
          l_use_ocff_direct, l_use_arclocff,                            &
          l_use_biogenic, l_use_arcldlta,                               &
          l_use_nitrate_direct, l_use_aod
          
      USE SW_CONTROL_STRUCT
      USE LW_CONTROL_STRUCT
      USE SPEC_SW_LW
      USE MCICA_MOD
      USE rad_pcf, ONLY: i_err_io, ip_cloud_mcica
      USE filenamelength_mod, ONLY :                                    & 
          filenamelength

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Control_Max_Sizes
      USE domain_params
      USE dust_parameters_mod, ONLY: l_dust
      USE um_input_control_mod,  ONLY:                                  &
           l_sulpc_so2,    l_soot,           l_biomass,                 &
           l_ocff,         l_nitrate
      USE ukca_option_mod, ONLY: l_ukca, l_ukca_radaer
      USE murk_inputs_mod, ONLY: l_murk, l_murk_rad

      IMPLICIT NONE

      INTEGER    ICODE         ! Return code : 0 Normal exit
!                              !             : >0 Error
      CHARACTER(LEN=80) CMESSAGE    ! Error message if ICODE > 0

!     Local variables

      INTEGER :: J, PATH_END, FILE_START
      INTEGER :: IERR_GET_FILE
!             Error flag returned by GET_FILE (not necessarily
!             consistent with the flags in rad_pcf).
      CHARACTER (LEN=filenamelength) :: spectral_file
!             Full path to spectral files (returned from GET_FILE)
      CHARACTER (LEN=200) :: MCICA_DATA
!             Path to McICA data file

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

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

      IF (lhook) CALL dr_hook('INITPHYS',zhook_in,zhook_handle)

!
!     ------------- Shortwave Radiation -------------------------
!
!        Get the full path to the main spectral file after expansion
!        of environment variables.
         CALL GET_FILE(57, SPECTRAL_FILE, filenamelength, IERR_GET_FILE)
         IF (IERR_GET_FILE /= 0) THEN
!           Convert the error flag from get_file to a flag recognised
!           by the radiation code.
            ICODE=I_ERR_IO
            CMESSAGE='Error reading name of shortwave spectral file.'
            IF (lhook) CALL dr_hook('INITPHYS',zhook_out,zhook_handle)
            RETURN
         END IF

         IF (SPECTRAL_FILE /= SW_CONTROL(1)%SPECTRAL_FILE) THEN
!           We need to expand the environment variables for all of
!           the spectral file paths.
            PATH_END=SCAN(SPECTRAL_FILE,'/',.TRUE.)
            DO J=1, N_SWCALL
               FILE_START=SCAN(SW_CONTROL(J)%SPECTRAL_FILE,'/',.TRUE.)
               SW_CONTROL(J)%SPECTRAL_FILE(PATH_END+1:) =               &
                 SW_CONTROL(J)%SPECTRAL_FILE(FILE_START+1:)
               SW_CONTROL(J)%SPECTRAL_FILE(:PATH_END) =                 &
                 SPECTRAL_FILE(:PATH_END)
            END DO
         END IF

         WRITE (6,'(/,A)') '********************** Reading SW '//       &
                         'spectral files **********************'
         DO J=1, N_SWCALL
            FILE_START=SCAN(SW_CONTROL(J)%SPECTRAL_FILE,'/',.TRUE.)
            SELECT CASE (                                               &
              SW_CONTROL(J)%SPECTRAL_FILE(FILE_START+1:FILE_START+3) )

            CASE ('ses','SES')
! DEPENDS ON: ses_inisw
              CALL ses_inisw( ICODE, CMESSAGE                           &
              , SW_CONTROL(J)%SPECTRAL_FILE                             &
              , SW_CONTROL(J)%L_O2                                      &
              , SW_CONTROL(J)%L_CH4                                     &
              , SW_CONTROL(J)%L_N2O                                     &
              , L_CLIMAT_AEROSOL                                        &
              , L_DUST, L_USE_ARCLDUST                                  &
              , L_SULPC_SO2, L_USE_ARCLSULP                             &
              , L_SOOT, L_USE_ARCLBLCK                                  &
              , L_BIOMASS, L_USE_ARCLBIOM                               &
              , L_USE_SEASALT_DIRECT, L_USE_ARCLSSLT                    &
              , L_OCFF, L_USE_ARCLOCFF                                  &
              , L_USE_BIOGENIC, L_USE_ARCLDLTA                          &
              , L_NITRATE                                               &
              , L_MURK_RAD                                              &
              , SW_CONTROL(J)%L_RAYLEIGH, SW_CONTROL(J)%L_GAS           &
              , SW_CONTROL(J)%L_CONTINUUM, SW_CONTROL(J)%L_DROP         &
              , SW_CONTROL(J)%L_AEROSOL, SW_CONTROL(J)%L_ICE            &
              , SW_CONTROL(J)%i_solar_src, SW_SPECTRUM(J)               &
                )

            CASE DEFAULT
! DEPENDS ON: r2_sw_specin
              CALL R2_SW_SPECIN(ICODE, CMESSAGE                         &
              , SW_CONTROL(J)%SPECTRAL_FILE                             &
              , SW_CONTROL(J)%L_O2                                      &
              , L_CLIMAT_AEROSOL                                        &
              , L_DUST, L_USE_ARCLDUST                                  &
              , L_SULPC_SO2, L_USE_ARCLSULP                             &
              , L_SOOT, L_USE_ARCLBLCK                                  &
              , L_BIOMASS, L_USE_ARCLBIOM                               &
              , L_USE_SEASALT_DIRECT, L_USE_ARCLSSLT                    &
              , L_OCFF, L_USE_ARCLOCFF                                  &
              , L_USE_BIOGENIC, L_USE_ARCLDLTA                          &
              , L_NITRATE                                               &
              , L_MURK_RAD                                              &
              , SW_CONTROL(J)%L_RAYLEIGH, SW_CONTROL(J)%L_GAS           &
              , SW_CONTROL(J)%L_CONTINUUM, SW_CONTROL(J)%L_DROP         &
              , SW_CONTROL(J)%L_AEROSOL, SW_CONTROL(J)%L_ICE            &
              , SW_SPECTRUM(J)                                          &
              )
            END SELECT
         END DO

         IF (ICODE /= 0) THEN 
           IF (lhook) CALL dr_hook('INITPHYS',zhook_out,zhook_handle)
           RETURN
         END IF
!
!     ------------- Longwave Radiation -------------------------
!
!        Get the full path to the main spectral file after expansion
!        of environment variables.
         CALL GET_FILE(80, SPECTRAL_FILE, filenamelength, IERR_GET_FILE)
         IF (IERR_GET_FILE /= 0) THEN
!           Convert the error flag from get_file to a flag recognised
!           by the radiation code.
            ICODE=I_ERR_IO
            CMESSAGE='Error reading name of longwave spectral file.'
            IF (lhook) CALL dr_hook('INITPHYS',zhook_out,zhook_handle)
            RETURN
         END IF

         IF (SPECTRAL_FILE /= LW_CONTROL(1)%SPECTRAL_FILE) THEN
!           We need to expand the environment variables for all of
!           the spectral file paths.
            PATH_END=SCAN(SPECTRAL_FILE,'/',.TRUE.)
            DO J=1, N_LWCALL
               FILE_START=SCAN(LW_CONTROL(J)%SPECTRAL_FILE,'/',.TRUE.)
               LW_CONTROL(J)%SPECTRAL_FILE(PATH_END+1:) =               &
                 LW_CONTROL(J)%SPECTRAL_FILE(FILE_START+1:)
               LW_CONTROL(J)%SPECTRAL_FILE(:PATH_END) =                 &
                 SPECTRAL_FILE(:PATH_END)
            END DO
         END IF

         WRITE (6,'(/,A)') '********************** Reading LW '//       &
                         'spectral files **********************'
         DO J=1,N_LWCALL
            FILE_START=SCAN(LW_CONTROL(J)%SPECTRAL_FILE,'/',.TRUE.)
            SELECT CASE (                                               &
              LW_CONTROL(J)%SPECTRAL_FILE(FILE_START+1:FILE_START+3) )

            CASE ('ses','SES')
! DEPENDS ON: ses_inilw
              CALL ses_inilw( ICODE, CMESSAGE                           &
              , LW_CONTROL(J)%SPECTRAL_FILE                             &
              , LW_CONTROL(J)%L_CH4, LW_CONTROL(J)%L_N2O                &
              , LW_CONTROL(J)%L_CFC11, LW_CONTROL(J)%L_CFC12            &
              , LW_CONTROL(J)%L_CFC113, LW_CONTROL(J)%L_CFC114          &
              , LW_CONTROL(J)%L_HCFC22, LW_CONTROL(J)%L_HFC125          &
              , LW_CONTROL(J)%L_HFC134A                                 &
              , L_CLIMAT_AEROSOL                                        &
              , L_DUST, L_USE_ARCLDUST                                  &
              , L_SULPC_SO2, L_USE_ARCLSULP                             &
              , L_SOOT, L_USE_ARCLBLCK                                  &
              , L_BIOMASS, L_USE_ARCLBIOM                               &
              , L_USE_SEASALT_DIRECT, L_USE_ARCLSSLT                    &
              , L_OCFF, L_USE_ARCLOCFF                                  &
              , L_USE_BIOGENIC, L_USE_ARCLDLTA                          &
              , L_NITRATE                                               &
              , L_MURK_RAD, L_USE_AOD                                   &
              , LW_CONTROL(J)%L_GAS, LW_CONTROL(J)%L_CONTINUUM          &
              , LW_CONTROL(J)%L_DROP, LW_CONTROL(J)%L_AEROSOL           &
              , LW_CONTROL(J)%L_ICE                                     &
              , LW_CONTROL(J)%i_solar_src, LW_SPECTRUM(J)               &
              )

            CASE DEFAULT
! DEPENDS ON: r2_lw_specin
              CALL R2_LW_SPECIN(ICODE, CMESSAGE                         &
              , LW_CONTROL(J)%SPECTRAL_FILE                             &
              , LW_CONTROL(J)%L_CH4, LW_CONTROL(J)%L_N2O                &
              , LW_CONTROL(J)%L_CFC11, LW_CONTROL(J)%L_CFC12            &
              , LW_CONTROL(J)%L_CFC113                                  &
              , LW_CONTROL(J)%L_HCFC22, LW_CONTROL(J)%L_HFC125          &
              , LW_CONTROL(J)%L_HFC134A                                 &
              , L_CLIMAT_AEROSOL                                        &
              , L_DUST, L_USE_ARCLDUST                                  &
              , L_SULPC_SO2, L_USE_ARCLSULP                             &
              , L_SOOT, L_USE_ARCLBLCK                                  &
              , L_BIOMASS, L_USE_ARCLBIOM                               &
              , L_USE_SEASALT_DIRECT, L_USE_ARCLSSLT                    &
              , L_OCFF, L_USE_ARCLOCFF                                  &
              , L_USE_BIOGENIC, L_USE_ARCLDLTA                          &
              , L_NITRATE                                               &
              , L_MURK_RAD                                              &
              , L_USE_AOD                                               &
              , LW_CONTROL(J)%L_GAS, LW_CONTROL(J)%L_CONTINUUM          &
              , LW_CONTROL(J)%L_DROP, LW_CONTROL(J)%L_AEROSOL           &
              , LW_CONTROL(J)%L_ICE                                     &
              , LW_SPECTRUM(J)                                          &
              )
            END SELECT
         END DO

         IF (ICODE /= 0) THEN 
           IF (lhook) CALL dr_hook('INITPHYS',zhook_out,zhook_handle)
           RETURN
         END IF

      IF ((sw_control(1)%i_cloud == ip_cloud_mcica) .OR.                &
          (lw_control(1)%i_cloud == ip_cloud_mcica)) THEN

        path_end=SCAN(sw_control(1)%spectral_file,'/',.TRUE.)
        mcica_data(:PATH_END)=sw_control(1)%spectral_file(:PATH_END)
        mcica_data(PATH_END+1:)='mcica_data'

        WRITE (6,'(/,A,A)') 'Reading McICA data file: ',mcica_data
        CALL read_mcica_data(mcica_data)
      END IF

!
!
!     ------------- UKCA_RADAER  -------------------------------
!
      IF (l_ukca_radaer) THEN
      
! DEPENDS ON: ukca_radaer_read_luts
        CALL ukca_radaer_read_luts(icode, cmessage)
        IF (icode /= 0) THEN
          IF (lhook) CALL dr_hook('INITPHYS',zhook_out,zhook_handle)
          RETURN
        END IF
        
! DEPENDS ON: ukca_radaer_read_precalc
        CALL ukca_radaer_read_precalc(icode, cmessage)
        IF (icode /= 0) THEN
          IF (lhook) CALL dr_hook('INITPHYS',zhook_out,zhook_handle)
          RETURN
        END IF
      
      END IF
!
!
      WRITE (6,'(A,/)') '**********************************'//       &
                        '*************************************'
      IF (lhook) CALL dr_hook('INITPHYS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE INITPHYS
