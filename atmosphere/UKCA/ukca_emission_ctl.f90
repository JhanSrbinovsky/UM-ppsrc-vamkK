! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
! Purpose: Subroutine to do boundary layer mixing of UKCA tracers,
!          simultaneously add surface emissions fields onto tracer
!          fields. Also, adds aircraft and lightning emissions.
!          Adapted from version supplied by Olaf Morgenstern.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_MAIN1.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
      SUBROUTINE UKCA_EMISSION_CTL(                                    &
          n_chem_tracers, n_mode_tracers, n_use_emissions,             &
          n_chem_emissions, timestep, em_chem_spec,                    &
          f3_at_u, r_rho_levels, r_theta_levels, sin_theta_latitude,   &
          FV_cos_theta_latitude,                                       &
          tan_theta_latitude, cos_zenith_angle, int_zenith_angle,      &
          true_longitude, delta_lambda, delta_phi,                     &
          iyear, imonth, iday,                                         &
          tropopause_height,                                           &
          ls_mask, conv_cloud_base, conv_cloud_top,                    &
          theta, q, qcl, qcf,                                          &
          exner_rho_levels, rho_r2,                                    &
          p_layer_boundaries, p_theta_levels, t_theta_levels,          &
          all_emissions, aircraftems,                                  &
          em_index,                                                    &
          BC_biom_3D,                                                  &
          OC_biom_3D,                                                  &
          SO2_volc_3D,                                                 &
          dust_flux, u_scalar_10m, rough_length,                       &
          land_fraction, seaice_frac, area,                            &
          z_half, alpha_cdx, ml_depth,                                 &
          rhokh_mix,                                                   &
          dtrdz_charney_grid, kent, we_lim,                            &
          t_frac, zrzi, kent_dsc,                                      &
          we_lim_dsc, t_frac_dsc,                                      &
          zrzi_dsc, zhsc,                                              &
          ch4_wetl_emiss, tracers, mode_tracers,                       &
          n_boundary_vals, lbc_spec, lbc_mmr,                          &
          mass, totnodens, volume,                                     &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
          len_stashwork38,STASHwork38)

      USE bl_option_mod,           ONLY: alpha_cd

      USE asad_mod,                ONLY: advt, method
      USE ukca_chem_schemes_mod,   ONLY: int_method_nr
      USE ASAD_CHEM_FLUX_DIAGS
      USE UKCA_CONSTANTS
      USE UKCA_TRACE_GAS_MIXRATIO, ONLY: um_ch4_for_ukca
      USE ukca_mode_ems_um_mod,    ONLY: ukca_mode_ems_um
      USE dust_parameters_mod,     ONLY: ndiv
      USE ukca_option_mod,         ONLY: L_ukca_stratcfc, L_ukca_strat,&
                                         L_ukca_prescribech4,          &
                                         L_ukca_strattrop,             &
                                         L_ukca_achem, L_ukca_aerchem, &
                                         L_ukca_chem, L_ukca_mode,     &
                                         jpctr, mode_parfrac
      USE run_aerosol_mod,         ONLY: SO2_high_level
      USE parkind1,                ONLY: jprb, jpim
      USE yomhook,                 ONLY: lhook, dr_hook
      USE ereport_mod,             ONLY: ereport
      USE tr_mix_mod,              ONLY: tr_mix
      USE trsrce_mod,              ONLY: trsrce
      USE UM_ParVars
      USE Control_Max_Sizes
! version_mod items required by cstash.h
      USE version_mod,             ONLY: nproftp, nprofdp, nprofup,    &
                                         ndiagpm, ntimep, NTimSerP,    &
                                         nlevp, npslevp, npslistp,     &
                                         outfile_s, outfile_e 
      USE Submodel_Mod
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
! COMDECK CSTASH
! Description:
!  Contains variables and arrays involved in STASH request specification
!  and STASH processing in the UM. Includes namelist STASH.
!
!
! Declarations:
! Imported global variables:
!    None, but see note above.

! Global parameters:

! Global scalars:
      INTEGER      NDIAG   ! No. of diagnostics
      INTEGER      NTPROF  ! No. of time profiles
      INTEGER      NSERIES ! No. of stash time series
      INTEGER      NDPROF  ! No. of domain profiles
      INTEGER      NUPROF  ! No. of useage profiles

! Global dynamic arrays:

!   STASH specification table (JSTASH file):
!   NDIAGPM set in VERSION_MOD module
      INTEGER      MODL_B(NDIAGPM)  ! Internal model no.
      INTEGER      ISEC_B(NDIAGPM)  ! Section
      INTEGER      ITEM_B(NDIAGPM)  ! Item
      INTEGER      ITIM_B(NDIAGPM)  ! Time profile number
      INTEGER      IDOM_B(NDIAGPM)  ! Domain profile number
      INTEGER      IUSE_B(NDIAGPM)  ! Useage profile number

!   Time profile information:

      CHARACTER(LEN=8)  TIMPRO(NPROFTP)         ! Name of profile
      INTEGER      ITYP_T(NPROFTP)         ! Type of profile
      INTEGER      INTV_T(NPROFTP)         ! Time Interval
      CHARACTER(LEN=2)  UNT1_T(NPROFTP)         ! Units for time interval
      INTEGER      ISAM_T(NPROFTP)         ! Sampling period
      CHARACTER(LEN=2)  UNT2_T(NPROFTP)         ! Units for sampling period
      INTEGER      IOPT_T(NPROFTP)         ! Output option
      INTEGER      ISTR_T(NPROFTP)         ! Output Start time
      INTEGER      IEND_T(NPROFTP)         ! Output End time
      INTEGER      ISDT_T(6, NPROFTP)      ! Output Start date
      INTEGER      IEDT_T(6, NPROFTP)      ! Output End date
      INTEGER      IFRE_T(NPROFTP)         ! Output frequency
      INTEGER      IOFF_T(NPROFTP)         ! Offset for sampling
      CHARACTER(LEN=2)  UNT3_T(NPROFTP)         ! Units for output times
      INTEGER      ITIM_T(NPROFTP)         ! No. of times in times table
      INTEGER      ISER_T(NTIMEP ,NPROFTP) ! Times table (with units)
      INTEGER      MODL_T(NPROFTP)         ! Indicates internal model
                                           !  for each times table

!   Domain profile information:

      CHARACTER(LEN=8) DOMPRO  (NPROFDP)           ! Name of domain profile
      INTEGER     IOPL_D  (NPROFDP)           ! Levels option
      INTEGER     LEVB_D  (NPROFDP)           ! Bottom level
      INTEGER     LEVT_D  (NPROFDP)           ! Top level
      INTEGER     IOPA_D  (NPROFDP)           ! Area option
      INTEGER     INTH_D  (NPROFDP)           ! North boundary
      INTEGER     ISTH_D  (NPROFDP)           ! South boundary
      INTEGER     IEST_D  (NPROFDP)           ! East boundary
      INTEGER     IWST_D  (NPROFDP)           ! West boundary
      INTEGER     IMSK_D  (NPROFDP)           ! Mask type
      INTEGER     IMN_D   (NPROFDP)           ! Meaning option
      INTEGER     IWT_D   (NPROFDP)           ! Weighting option
      CHARACTER(LEN=1) TS_D    (NPROFDP)           ! Time series profile
      INTEGER     IG_TS
      INTEGER     I1_TS
      INTEGER     I51_TS
      INTEGER     BLIM_TS (NTimSerP)
      INTEGER     TLIM_TS (NTimSerP)
      REAL        BLIMR_TS(NTimSerP)
      REAL        TLIMR_TS(NTimSerP)
      INTEGER     NLIM_TS (NTimSerP)
      INTEGER     SLIM_TS (NTimSerP)
      INTEGER     ELIM_TS (NTimSerP)
      INTEGER     WLIM_TS (NTimSerP)
      INTEGER     ILEV_D  (NPROFDP)           ! Output levels code
      INTEGER     LEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      REAL       RLEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      INTEGER     PLT_D   (NPROFDP)
      INTEGER     PLLEN_D (NPROFDP)
      INTEGER     PLPOS_D (NPROFDP)
      INTEGER     PSLIST_D(NPSLEVP ,NPSLISTP)
      INTEGER     NPSLISTS
      EQUIVALENCE        (RLEVLST_D,LEVLST_D)

! Useage information:

      CHARACTER(LEN=8) USEPRO(NPROFUP)   ! Name of useage profile
      INTEGER     LOCN_U(NPROFUP)   ! Storage location of profile
      INTEGER     IUNT_U(NPROFUP)   ! Unit no.

! Information from ppxref file:

      INTEGER      MODEL_ST       ! Internal model number
      INTEGER      ISPACE         ! Space code
      INTEGER      ITIMA          ! Time availability code
      INTEGER      IGP            ! Grid of data code
      INTEGER      ILEV           ! Level type code
      INTEGER      IBOT           ! First level code
      INTEGER      ITOP           ! Last level code
      INTEGER      IFLAG          ! Level compression flag
      INTEGER      IOPN(6)        ! Sectional option code
      INTEGER      VMSK           ! Integer equiv of bin vers mask
      INTEGER      IPSEUDO        ! Pseudo dimension type
      INTEGER      IPFIRST        ! First pseudo dim code
      INTEGER      IPLAST         ! Last pseudo dim code
      INTEGER      PTR_PROG       ! Section zero point back
      INTEGER      HALO_TYPE      ! Type of halo the field has

! PP output file units
      INTEGER      PPlen2LkUp(OUTFILE_S:OUTFILE_E)
      CHARACTER(LEN=1)  FTOutUnit (OUTFILE_S:OUTFILE_E)

! COMMON blocks:
      COMMON/STCHA/ TIMPRO,UNT1_T,UNT2_T,UNT3_T,DOMPRO,TS_D,            &
     &  USEPRO,FTOutUnit

      COMMON/STSH/                                                      &
     &  NDIAG   ,MODL_B  ,ISEC_B ,ITEM_B  ,ITIM_B  ,IDOM_B  ,IUSE_B,    &
     &  NTPROF  ,ITYP_T  ,INTV_T ,ISAM_T  ,ITIM_T  ,                    &
     &  IOPT_T  ,ISTR_T  ,IEND_T ,IFRE_T  ,IOFF_T, ISER_T  ,MODL_T  ,   &
     &  NDPROF  ,IOPL_D  ,LEVB_D ,ISDT_T  ,IEDT_T  ,                    &
     &  IOPA_D  ,INTH_D  ,ISTH_D ,IEST_D  ,IWST_D  ,                    &
     &  IMSK_D  ,IMN_D   ,IWT_D  ,                                      &
     &  LEVT_D  ,LEVLST_D,                                              &
     &  PLT_D   ,PLLEN_D ,PLPOS_D,PSLIST_D,NPSLISTS,                    &
     &  BLIM_TS ,TLIM_TS ,BLIMR_TS,TLIMR_TS,IG_TS   ,I1_TS   ,          &
     &  NLIM_TS ,SLIM_TS ,ELIM_TS ,WLIM_TS ,I51_TS  ,NSERIES ,          &
     &  NUPROF  ,LOCN_U  ,IUNT_U ,                                      &
     &  MODEL_ST,ISPACE  ,ITIMA  ,IGP     ,                             &
     &  ILEV    ,IBOT    ,ITOP   ,IFLAG   ,IOPN    ,VMSK    ,           &
     &  IPSEUDO ,IPFIRST ,IPLAST ,PTR_PROG, HALO_TYPE,                  &
     & PPlen2LkUp

! CSTASH end
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

      INTEGER, PARAMETER  :: surface_level = 1

      INTEGER, INTENT(IN) :: n_chem_tracers          ! No of chemical tracers
      INTEGER, INTENT(IN) :: n_mode_tracers          ! No of aerosol tracers
      INTEGER, INTENT(IN) :: n_use_emissions         ! Total No of emissions
      INTEGER, INTENT(IN) :: n_chem_emissions
! No of emissions in all_emissions and em_index arrays (chemical and aerosol)
      INTEGER, INTENT(IN) :: n_boundary_vals         ! No species with  b.c's
      INTEGER, INTENT(IN) :: conv_cloud_base(1:row_length, 1:rows)
      INTEGER, INTENT(IN) :: conv_cloud_top(1:row_length, 1:rows)
      INTEGER, INTENT(IN) :: kent(1:row_length, 1:rows)
      INTEGER, INTENT(IN) :: kent_dsc(1:row_length, 1:rows)
      INTEGER, INTENT(IN) :: em_index(n_chem_emissions) ! Index to item numbers
      INTEGER, INTENT(IN) :: iyear
      INTEGER, INTENT(IN) :: imonth
      INTEGER, INTENT(IN) :: iday
      INTEGER, INTENT(IN) :: len_stashwork38      ! length of diagnostic array

      LOGICAL, INTENT(IN) :: ls_mask(1:row_length, 1:rows)

      REAL, INTENT(IN) :: timestep
      REAL, INTENT(IN) :: all_emissions(row_length,rows,n_chem_emissions)
! All emissions from 2-D ancillaries
      REAL, INTENT(IN) :: aircraftems(row_length,rows,model_levels)  
! Aircraft NOx ems
      REAL, INTENT(IN) :: BC_biom_3D(row_length,rows,model_levels)   
! BC biomass ems
      REAL, INTENT(IN) :: OC_biom_3D(row_length,rows,model_levels)   
! OC biomass ems
      REAL, INTENT(IN) :: SO2_volc_3D(row_length,rows,model_levels)  
! SO2 volcanic ems
      REAL, INTENT(IN) :: area(row_length, rows, model_levels)       
! area of grid cell
      REAL, INTENT(IN) :: dust_flux(row_length,rows,ndiv)            
! dust emissions (kg/m2/s)
      REAL, INTENT(IN) :: u_scalar_10m(row_length,rows)              
! scalar wind at 10m (m/s)
      REAL, INTENT(IN) :: rough_length(row_length, rows)             
! roughness length (m)
      REAL, INTENT(IN) :: land_fraction(row_length,rows)             
! land_fraction
      REAL, INTENT(IN) :: seaice_frac(row_length, rows)              
! sea ice fraction
      REAL, INTENT(IN) :: lbc_mmr(n_boundary_vals)
      REAL, INTENT(IN) :: f3_at_u(1:row_length,1:rows)
      REAL, INTENT(IN) :: r_rho_levels(1:row_length, 1:rows,            &
                             1:model_levels)        ! ht of rho levs
      REAL, INTENT(IN) :: r_theta_levels(1:row_length, 1:rows,          &
                             0:model_levels)        ! ht of theta levs
      REAL, INTENT(IN) :: sin_theta_latitude(1:row_length,1:rows)
      REAL, INTENT(IN) :: FV_cos_theta_latitude(1:row_length,1:rows)
      REAL, INTENT(IN) :: tan_theta_latitude(1:row_length,1:rows)
      REAL, INTENT(IN) :: cos_zenith_angle(1:row_length,1:rows)
      REAL, INTENT(IN) :: int_zenith_angle(1:row_length,1:rows)
      REAL, INTENT(IN) :: true_longitude(row_length, rows)
      REAL, INTENT(IN) :: delta_lambda
      REAL, INTENT(IN) :: delta_phi
      REAL, INTENT(IN) :: theta(1:row_length,1:rows,model_levels)
      REAL, INTENT(IN) :: q(1:row_length,1:rows,1:model_levels)
      REAL, INTENT(IN) :: qcl(1:row_length,1:rows,1:model_levels)
      REAL, INTENT(IN) :: qcf(1:row_length,1:rows,1:model_levels)
      REAL, INTENT(IN) :: exner_rho_levels(1:row_length,1:rows,         &
                                             1:model_levels+1)
      REAL, INTENT(IN) :: rho_r2(1:row_length,1:rows,model_levels)
      REAL, INTENT(IN) :: p_layer_boundaries(1:row_length,1:rows,       &
                                               0:model_levels)
      REAL, INTENT(IN) :: p_theta_levels(1:row_length,1:rows,           &
                                           1:model_levels)
      REAL, INTENT(IN) :: t_theta_levels(1:row_length,1:rows,           &
                                           1:model_levels)
      REAL, INTENT(IN) :: z_half(1:row_length,1:rows,1:bl_levels)
      REAL, INTENT(IN) :: alpha_cdx(1:bl_levels)
      REAL, INTENT(IN) :: ml_depth(1:row_length,1:rows)
      REAL, INTENT(IN) :: tropopause_height(row_length, rows)
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
      REAL, INTENT(IN) :: mass(row_length, rows, model_levels)
      REAL, INTENT(IN) :: totnodens(row_length, rows, model_levels)
      REAL, INTENT(IN) :: volume(row_length, rows, model_levels)

      CHARACTER(LEN=10), INTENT(IN) :: em_chem_spec(n_use_emissions)
      CHARACTER(LEN=10), INTENT(IN) :: lbc_spec(n_boundary_vals)

! Diagnostic array
      REAL,INTENT(INOUT) :: STASHwork38(len_stashwork38)
! Methane wetland emissions
      REAL,INTENT(INOUT) :: ch4_wetl_emiss(1:row_length,1:rows)
! Tracer MMRs
      REAL,INTENT(INOUT) :: tracers(row_length,rows,model_levels,       &
                                      n_chem_tracers)
      REAL,INTENT(INOUT) :: mode_tracers(row_length,rows,model_levels,  &
                                      n_mode_tracers)

! Local variables

      INTEGER, SAVE :: inox           ! Index for NO/NOx tracer
      INTEGER, SAVE :: iso2           ! Index for SO2 tracer
      INTEGER, SAVE :: jso2_high      ! HL SO2 emissn index
      INTEGER :: j                    ! Loop variable
      INTEGER :: k                    ! Loop variable
      INTEGER :: kaer                 ! Loop variable for UKCA aerosol tracers
      INTEGER :: l                    ! Loop variable
      INTEGER :: n                    ! counter
      INTEGER :: errcode              ! Error code for ereport
      INTEGER :: icode                ! Variable passed to asad diags routines


      INTEGER :: ils_mask(row_length,rows)  ! Land/sea mask (1/0)
      INTEGER :: em_count                   ! counter

      REAL, PARAMETER :: meoh_factor=0.215
! Factor to convert
!  GEIA NVOC emissions of methanol (260.95 Tg as C) - assuming Biogenic
!  emission of methanol is 151 Tg (56 Tg as C)
      REAL :: res_factor(row_length,rows)                      
! dry
      REAL :: tr_flux(row_length,rows,bl_levels,n_chem_tracers)
! tracer flux (chem)
      REAL :: tr_flux_mode(row_length,rows,bl_levels)          
! tracer flux (mode)
      REAL :: surf_dep_flux(row_length,rows,n_chem_tracers)    
! surface deposition flux
      REAL :: surf_dep_flux_mode(row_length,rows)              
! surface deposition flux (mode)
      REAL :: em_field(row_length,rows,n_chem_tracers)         
! surf
      REAL :: em_field_mode(row_length,rows,model_levels,n_mode_tracers)
! all levels emissions (Mode)
      REAL :: conv_aircraftems(row_length,rows,model_levels)   
! aircraft emissions
      REAL :: lightningems(row_length,rows,model_levels)       
! nox lightning emissions
      REAL :: conv_SO2emiss_3D(row_length,rows,model_levels)   
! 3D SO2 emissions
      REAL :: tmp_in_em_field(row_length,rows)
      REAL :: tmp_out_em_field(row_length,rows)
      REAL :: rho_aresist(1:row_length,1:rows)

      REAL, SAVE, ALLOCATABLE :: surf_area(:,:)                ! gridbox area
      REAL, SAVE, ALLOCATABLE :: theta_latitude(:,:)
      REAL, SAVE, ALLOCATABLE :: molmass(:)
      REAL, SAVE, ALLOCATABLE :: lbc_molmass(:)

! required for diagnostics
      REAL, DIMENSION(:,:,:), ALLOCATABLE  :: tmp3dems   ! volcanic ems
       
      INTEGER :: ierr                                    ! error indicator

      LOGICAL, SAVE :: firstcall = .TRUE.
      LOGICAL, SAVE :: testdcycl = .FALSE.

      CHARACTER(LEN=72)  :: cmessage                          ! Error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!     Initialise variables

      IF (lhook) CALL dr_hook('UKCA_EMISSION_CTL',zhook_in,zhook_handle)
      rho_aresist(:,:)  = 0.0     ! Set dry deposition to zero
      res_factor(:,:)  = 0.0      ! Set dry deposition to zero
      em_field (:,:,:) = 0.0      ! Initial ems field for each tracer

      IF (firstcall) THEN
        inox             = -99      ! Initial index for nox tracer
        iso2             = -99      ! Initial index for so2 tracer
        jso2_high        = -99      ! Initial index for SO2 HL emissions

! Find index for SO2 emissions 
        IF(L_ukca_aerchem .OR. L_ukca_achem) THEN
          DO k=1,n_use_emissions
            IF (em_chem_spec(k)(1:8) == 'SO2_high') THEN
              jso2_high=k
              EXIT
            END IF
          END DO
        END IF

        ALLOCATE(theta_latitude(row_length,rows))

        theta_latitude = asin(sin_theta_latitude)

!       Calculate the gridbox surface area

        ALLOCATE(surf_area(row_length,rows))
        DO k = 1, rows
          DO j = 1, row_length
            surf_area(j,k) = r_theta_levels(j,k,0)                   &
                           * r_theta_levels(j,k,0)                   &
                           * delta_lambda * delta_phi                &
                           * FV_cos_theta_latitude(j,k)
         END DO
        END DO


!       Find index for SO2 and NOx tracer

        DO k = 1,jpctr
          SELECT CASE (advt(k))
            CASE('NOx       ')
              inox = k
           CASE('NO        ')
              inox = k
           CASE('SO2       ')
              iso2 = k
          END SELECT
        ENDDO

        IF (inox == -99) then
          cmessage = 'Did not find NO or NOx tracer'
          errcode=1
          CALL EREPORT('UKCA_EMISSION_CTL',errcode,cmessage)
        ENDIF

        IF ((L_ukca_aerchem .OR. L_ukca_achem)   &
          .AND. iso2 == -99) THEN
          cmessage = 'Did not find SO2 tracer'
          errcode=1
          CALL EREPORT('UKCA_EMISSION_CTL',errcode,cmessage)
        ENDIF

        IF (L_ukca_chem .AND. method == int_method_NR) THEN    ! N-R solver

          IF (ANY(em_chem_spec == 'NO_aircrft')) THEN
            ALLOCATE(molmass(n_use_emissions-1))
          ELSE
            ALLOCATE(molmass(n_use_emissions))
          END IF

            molmass = 0.
            WHERE (em_chem_spec == 'NOx       ') molmass = m_no2
            WHERE (em_chem_spec == 'NO2       ') molmass = m_no2
            WHERE (em_chem_spec == 'NO        ') molmass = m_no
            WHERE (em_chem_spec == 'NO_aircrft') molmass = m_no2
            WHERE (em_chem_spec == 'CH4       ') molmass = m_ch4
            WHERE (em_chem_spec == 'CO        ') molmass = m_co
            WHERE (em_chem_spec == 'HCHO      ') molmass = m_hcho
            WHERE (em_chem_spec == 'C2H6      ') molmass = m_c2h6
            WHERE (em_chem_spec == 'C3H8      ') molmass = m_c3h8
            WHERE (em_chem_spec == 'Me2CO     ') molmass = m_me2co
            WHERE (em_chem_spec == 'MeCHO     ') molmass = m_mecho
            WHERE (em_chem_spec == 'DMS       ') molmass = m_dms
            WHERE (em_chem_spec == 'SO2       ') molmass = m_so2
            WHERE (em_chem_spec == 'SO2_low   ') molmass = m_so2
            WHERE (em_chem_spec == 'SO2_high  ') molmass = m_so2
            WHERE (em_chem_spec == 'SO2_nat   ') molmass = m_so2
            WHERE (em_chem_spec == 'Me2S      ') molmass = m_me2s
            WHERE (em_chem_spec == 'COS       ') molmass = m_ocs
            WHERE (em_chem_spec == 'H2S       ') molmass = m_h2s
            WHERE (em_chem_spec == 'CS2       ') molmass = m_cs2
            WHERE (em_chem_spec == 'NH3       ') molmass = m_nh3
            WHERE (em_chem_spec == 'N2O       ') molmass = m_n2o
            WHERE (em_chem_spec == 'CF2Cl2    ') molmass = m_cf2cl2
            WHERE (em_chem_spec == 'CFCl3     ') molmass = m_cfcl3
            WHERE (em_chem_spec == 'MeBr      ') molmass = m_mebr
            WHERE (em_chem_spec == 'CF2ClCFCl2') molmass = m_cf2clcfcl2
            WHERE (em_chem_spec == 'MeCl      ') molmass = m_mecl
            WHERE (em_chem_spec == 'MeCCl3    ') molmass = m_meccl3
            WHERE (em_chem_spec == 'CHF2Cl    ') molmass = m_chf2cl
            WHERE (em_chem_spec == 'CFCl3     ') molmass = m_cfcl3
            WHERE (em_chem_spec == 'MeBr      ') molmass = m_mebr
            WHERE (em_chem_spec == 'CF2ClCFCl2') molmass = m_cf2clcfcl2
            WHERE (em_chem_spec == 'MeCl      ') molmass = m_mecl
            WHERE (em_chem_spec == 'MeCCl3    ') molmass = m_meccl3
            WHERE (em_chem_spec == 'CHF2Cl    ') molmass = m_chf2cl
            WHERE (em_chem_spec == 'CCl4      ') molmass = m_ccl4
            WHERE (em_chem_spec == 'CF2ClBr   ') molmass = m_cf2clbr
            WHERE (em_chem_spec == 'CF3Br     ') molmass = m_cf3br
            WHERE (em_chem_spec == 'C5H8      ') molmass = m_isop
            WHERE (em_chem_spec == 'SO4       ') molmass = m_so4
            WHERE (em_chem_spec == 'H2        ') molmass = m_h2
            WHERE (em_chem_spec == 'Monoterp  ') molmass = m_monoterp
            WHERE (em_chem_spec == 'NVOC      ') molmass = m_c
            WHERE (em_chem_spec == 'C4H10     ') molmass = m_c4h10
!          WHERE (em_chem_spec == 'C10H16    ') molmass = m_c10h16
!          WHERE (em_chem_spec == 'N2O       ') molmass = m_n2o
            WHERE (em_chem_spec == 'MeOH      ') molmass = m_meoh
            WHERE (em_chem_spec == 'C2H4      ') molmass = m_c2h4
            WHERE (em_chem_spec == 'C3H6      ') molmass = m_c3h6
            WHERE (em_chem_spec == 'TOLUENE   ') molmass = m_toluene
            WHERE (em_chem_spec == 'oXYLENE   ') molmass = m_oxylene
            WHERE (em_chem_spec == 'BC_fossil ') molmass = m_c
            WHERE (em_chem_spec == 'BC_biofuel') molmass = m_c
            WHERE (em_chem_spec == 'BC_biomass') molmass = m_c
            WHERE (em_chem_spec == 'OC_fossil ') molmass = m_c
            WHERE (em_chem_spec == 'OC_biofuel') molmass = m_c
            WHERE (em_chem_spec == 'OC_biomass') molmass = m_c

! Check if all the emitted species have a valid molecular weight
            IF (ANY(molmass(:) < 0.00001)) THEN
              n = 0
              DO j=1,SIZE(molmass)
                IF (molmass(j) < 0.00001) THEN
                  cmessage = ' Species: '//em_chem_spec(j)//              &
                       ' missing from molmass list.'
                  errcode = -j
                  CALL EREPORT('UKCA_EMISSION_CTL',errcode,cmessage)
                  n = n+1
                END IF
              END DO
              IF (n > 0) THEN
                cmessage = ' Species missing from molmass list'
                errcode = n
                CALL EREPORT('UKCA_EMISSION_CTL',errcode,cmessage)
              END IF
            END IF

          IF (L_ukca_strat .OR. L_ukca_strattrop .OR. L_ukca_stratcfc)  &
                 THEN

! Set MMRs of N2O, CFCs and halons to specified global constants,
! which may follow a time dependence. Set values of inorganic
! Cl and Br compounds to 0 at lower boundary.
! Adjust this if new stratospheric species are introduced.
! The diagnostic (surface emission) contains the additional
! amount of tracer in mols added globally, which is negative
! in case of sink gases.

           ALLOCATE(lbc_molmass(n_boundary_vals))
           lbc_molmass = 0.
           WHERE (lbc_spec == 'no_N2O    ') lbc_molmass = m_n2o
           WHERE (lbc_spec == 'no_H2     ') lbc_molmass = m_h2
           WHERE (lbc_spec == 'N2O       ') lbc_molmass = m_n2o
           WHERE (lbc_spec == 'CH4       ') lbc_molmass = m_ch4
           WHERE (lbc_spec == 'CFCl3     ') lbc_molmass = m_cfcl3
           WHERE (lbc_spec == 'CF2Cl2    ') lbc_molmass = m_cf2cl2
           WHERE (lbc_spec == 'HCl       ') lbc_molmass = m_hcl
           WHERE (lbc_spec == 'HOCl      ') lbc_molmass = m_hocl
           WHERE (lbc_spec == 'ClONO2    ') lbc_molmass = m_clono2
           WHERE (lbc_spec == 'Clx       ') lbc_molmass = m_clo
           WHERE (lbc_spec == 'OClO      ') lbc_molmass = m_oclo
           WHERE (lbc_spec == 'MeBr      ') lbc_molmass = m_mebr
           WHERE (lbc_spec == 'HBr       ') lbc_molmass = m_hbr
           WHERE (lbc_spec == 'HOBr      ') lbc_molmass = m_hobr
           WHERE (lbc_spec == 'BrONO2    ') lbc_molmass = m_brono2
           WHERE (lbc_spec == 'Brx       ') lbc_molmass = m_bro
           WHERE (lbc_spec == 'BrCl      ') lbc_molmass = m_brcl
           WHERE (lbc_spec == 'CF2ClCFCl2') lbc_molmass = m_cf2clcfcl2
           WHERE (lbc_spec == 'MeCl      ') lbc_molmass = m_mecl
           WHERE (lbc_spec == 'MeCCl3    ') lbc_molmass = m_meccl3
           WHERE (lbc_spec == 'CHF2Cl    ') lbc_molmass = m_chf2cl
           WHERE (lbc_spec == 'CCl4      ') lbc_molmass = m_ccl4
           WHERE (lbc_spec == 'CF2ClBr   ') lbc_molmass = m_cf2clbr
           WHERE (lbc_spec == 'CF3Br     ') lbc_molmass = m_cf3br
           WHERE (lbc_spec == 'CH2Br2    ') lbc_molmass = m_ch2br2
           WHERE (lbc_spec == 'H2        ') lbc_molmass = m_h2
           WHERE (lbc_spec == 'COS       ') lbc_molmass = m_ocs
           WHERE (lbc_spec == 'Cl        ') lbc_molmass = m_cl
           WHERE (lbc_spec == 'ClO       ') lbc_molmass = m_clo
           WHERE (lbc_spec == 'Cl2O2     ') lbc_molmass = m_cl2o2
           WHERE (lbc_spec == 'Br        ') lbc_molmass = m_br
           WHERE (lbc_spec == 'TOT_Cl    ') lbc_molmass = m_cl
           WHERE (lbc_spec == 'TOT_Br    ') lbc_molmass = m_br
           WHERE (lbc_spec == 'BrO       ') lbc_molmass = m_bro
           WHERE (lbc_spec == 'AGE       ') lbc_molmass = 1.
           WHERE (lbc_spec == 'AGE OF AIR') lbc_molmass = 1.
           WHERE (lbc_spec == 'PASSIVE O3') lbc_molmass = 1.
           WHERE (lbc_spec == 'XXX       ') lbc_molmass = 1.   ! unused species

           DO j=1,SIZE(lbc_molmass)
             IF (lbc_molmass(j) < 0.00001) THEN
               cmessage=' Species '//lbc_spec(j)//                       &
                        ' missing from LBC molmass list.'
               CALL EREPORT('UKCA_EMISSION_CTL',j,cmessage)
             END IF
           END DO
         END IF      ! L_ukca_strat etc
       END IF      ! L_ukca_chem etc

        firstcall = .FALSE.
      ENDIF      ! firstcall etc

      IF (L_ukca_mode) THEN       ! primary emissions for UKCA-MODE

      CALL ukca_mode_ems_um(imonth, iday,                               &
                      row_length, rows, model_levels,                   &
                      n_mode_tracers, ndiv,                             &
                      n_use_emissions, n_chem_emissions,                &
                      area,                                             &
                      mass,                                             &
                      p_theta_levels,                                   &
                      t_theta_levels,                                   &
                      seaice_frac,                                      &
                      rough_length,                                     &
                      u_scalar_10m,                                     &
                      land_fraction,                                    &
                      all_emissions,                                    &
                      em_index,                                         &
                      SO2_volc_3D,                                      &
                      BC_biom_3D,                                       &
                      OC_biom_3D,                                       &
                      dust_flux,                                        &
                      em_field_mode,                                    &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                      len_stashwork38,STASHwork38)

      END IF     ! l_ukca_mode

!       Set methane emissions to zero over non-land surfaces

      WHERE(ch4_wetl_emiss < 0.0) ch4_wetl_emiss = 0.0

      em_count = 0   ! set emission counter to 0


!   Check if tracer has surface emissions and set emission.
!  otherwise emission field is zero from initialisation.

      DO k = 1,jpctr              ! loop over tracers

        DO l=1,n_chem_emissions
          IF (advt(k) == em_chem_spec(l) .AND.                          &
            em_chem_spec(l) == 'NO      ' ) THEN
!          Convert from kg NO2/m2/s to kg NO/m2/s
            em_field(:,:,k) = all_emissions(:,:,l)*m_no/m_no2
          ELSE IF (advt(k) == em_chem_spec(l)(1:3) .AND.                &
            em_chem_spec(l) == 'SO2_low ' ) THEN
!          Convert from kg S/m2/s to kg SO2/m2/s and take off sulphate fraction
            em_field(:,:,k) = all_emissions(:,:,l)*                     &
                               (1.0 - mode_parfrac/100.0)*m_so2/m_s
          ELSE IF (advt(k) == em_chem_spec(l) .AND.                     &
            em_chem_spec(l) == 'DMS     ' ) THEN
!          Convert from kg S/m2/s to kg DMS/m2/s
            em_field(:,:,k) = all_emissions(:,:,l)*m_dms/m_s
          ELSE IF (advt(k) == 'MeOH      ' .AND.                        &
            em_chem_spec(l) == 'NVOC      ' ) THEN
!          Convert from kg C/m2/s to kg CH3OH/m2/s
            em_field(:,:,k) = all_emissions(:,:,l)*meoh_factor*         &
                m_meoh/(m_c*3.0)
          ELSE IF (advt(k) == em_chem_spec(l) .AND.                     &
            em_chem_spec(l) == 'Monoterp  ' ) THEN
!          Convert from kg C/m2/s to kg C10H16/m2/s
            em_field(:,:,k) = all_emissions(:,:,l)*m_monoterp/(m_c*10.0)
!           === biogenic emissions ===
          ELSE IF (advt(k) == em_chem_spec(l) .AND.                     &
                  em_chem_spec(l) == 'C5H8      ') THEN
            IF (L_ukca_diurnal_isopems) THEN
              tmp_in_em_field(:,:) = all_emissions(:,:,l)*(m_isop/      &
                                                            (5.0*m_c))
! DEPENDS ON: ukca_diurnal_isop_ems
!                testdcycl = .TRUE.
                CALL UKCA_DIURNAL_ISOP_EMS(row_length, rows,             &
                          tmp_in_em_field,  cos_zenith_angle,            &
                          int_zenith_angle,                              &
                          sin_theta_latitude, FV_cos_theta_latitude,     &
                          tan_theta_latitude, timestep, tmp_out_em_field,&
                          testdcycl)
                em_field(:,:,k) = tmp_out_em_field(:,:)
            ELSE
                em_field(:,:,k) = all_emissions(:,:,l)*(m_isop/(5.0*m_c))
            END IF
          ELSE IF (advt(k) == em_chem_spec(l) ) THEN
            em_field(:,:,k) = all_emissions(:,:,l)
          ENDIF             ! end advt(k)
 
        END DO       ! l=1,n_use_emissions

!       Add on wetland methane emissions, converting from kgC/m2/s 
!       to kgCH4/m2/s, when L_ukca_prescribech4 is false. Otherwise, 
!       fix the surface ch4 mass mixing ratio using um_ch4_for_ukca 
!       from UKCA_TRACE_GAS_MIXRATIO module 

        IF (advt(k) == 'CH4       ' .AND. (.NOT. L_ukca_prescribech4) ) THEN 
          em_field(:,:,k) = em_field(:,:,k) +                        & 
                            ch4_wetl_emiss*m_ch4/m_c 
        ELSE IF (advt(k) == 'CH4       ' .AND. L_ukca_prescribech4 ) THEN 
          tracers(:,:,surface_level,k) = um_ch4_for_ukca 
          em_field(:,:,k) = 0.0
        ENDIF

          IF (L_ukca_strat .OR. L_ukca_stratcfc .OR.                   &
              L_ukca_strattrop)  THEN
            DO l=1,n_boundary_vals
              IF (advt(k) == lbc_spec(l)) THEN
! Set emissions equal to difference of tracer at surface and intented value,
! scaled with mass in gridbox, area and timestep, to turn it into a
! surface emission rate.
                em_field(:,:,k) = (lbc_mmr(l) - tracers(:,:,1,k))      &
                  * mass(:,:,1) / surf_area / timestep
              END IF
            END DO
          END IF   ! L_ukca_strat etc

!      Call boundary layer mixing and add surface emissions.
!      Exclude H2O tracer here if an advected tracer. Note 
!      that in climate model q is not mixed either.

          IF (advt(k) /= 'H2O       ') THEN
            CALL TR_MIX(                                               &
            bl_levels, alpha_cd,                                       &
            rhokh_mix(:,:,2:), rho_aresist,                            &
            dtrdz_charney_grid, em_field(:,:,k), res_factor(:,:),      &
            kent, we_lim, t_frac, zrzi,                                &
            kent_dsc, we_lim_dsc, t_frac_dsc,                          &
            zrzi_dsc, ml_depth, zhsc, z_half,                          &
            tracers(:,:,1:bl_levels,k),tr_flux(:,:,1:bl_levels,k),     &
            surf_dep_flux(:,:,k)                                       &
            )

          ENDIF

        ENDDO                       ! end of loop over tracers



! Calculate emissions diagnostics
        IF (L_asad_use_chem_diags .AND. L_asad_use_surf_ems)            &
             CALL asad_emissions_diagnostics(                           &
             row_length, rows, model_levels,                            &
             n_chem_tracers, em_field, surf_area,                       &
             timestep, n_use_emissions, n_boundary_vals,                &
             em_chem_spec, lbc_spec, molmass,                           &
             lbc_molmass, ierr) 

        IF (L_ukca_mode) THEN
!  Call boundary layer mixing and add surface emissions for GLOMAP-mode
!   aerosol tracers

         DO kaer = 1,n_mode_tracers              ! loop over tracers

           tr_flux_mode(:,:,:) = 0.0
           surf_dep_flux_mode(:,:) = 0.0

          CALL TR_MIX(                                                  &
          bl_levels, alpha_cd,                                          &
          rhokh_mix(:,:,2:), rho_aresist,                               &
          dtrdz_charney_grid,em_field_mode(:,:,1,kaer),res_factor(:,:), &
          kent, we_lim, t_frac, zrzi,                                   &
          kent_dsc, we_lim_dsc, t_frac_dsc,                             &
          zrzi_dsc, ml_depth, zhsc, z_half,                             &
          mode_tracers(:,:,1:bl_levels,kaer),                           &
          tr_flux_mode(:,:,1:bl_levels), surf_dep_flux_mode(:,:)        &
          )

! Add in emission fluxes not at surface level
          DO k = 2,model_levels
            IF (SUM(em_field_mode(:,:,k,kaer)) > 1e-30) THEN
              CALL TRSRCE( rows, row_length, 0, 0, 0, 0,                &
                model_levels, wet_levels, 0, 0,                         &
                theta, q , qcl , qcf , exner_rho_levels, rho_r2,        &
                mode_tracers(:,:,k,kaer), em_field_mode(:,:,k,kaer),    &
                k, timestep, 1, 1, 0.0)
            END IF
          END DO ! loop over model levels 2 to top (k)

         END DO                       ! end of loop over tracers (kaer)

        END IF ! if L_UKCA_MODE

!       Set up integer land/sea mask 
         
        DO l = 1,rows 
          DO k = 1,row_length 
            IF (ls_mask(k,l)) THEN 
              ils_mask(k,l) = 1.0 
            ELSE 
              ils_mask(k,l) = 0.0 
            ENDIF 
          ENDDO 
        ENDDO 

!     Diagnose NO2 lightning emissions 

! DEPENDS ON: ukca_light_ctl 
      CALL UKCA_LIGHT_CTL(                                             &
          rows,row_length,delta_lambda,delta_phi,model_levels,         & 
          conv_cloud_base(1:row_length,1:rows),                        &
          conv_cloud_top(1:row_length,1:rows),                         &
          ils_mask(1:row_length,1:rows),                               &
          Recip_Pi_Over_180                                            &
          *asin(f3_at_u(1:row_length,1:rows)/two_omega),               &
          surf_area(1:row_length,1:rows),                              &
          r_theta_levels(1:row_length,1:rows,0:model_levels),          &
          r_rho_levels(1:row_length, 1:rows,1:model_levels),           &
          p_theta_levels(1:row_length,1:rows,1:model_levels),          &
          p_layer_boundaries(1:row_length,1:rows,0:model_levels),      &
          lightningems(1:row_length,1:rows,1:model_levels) ) 

          ! diagnostics - lightning emissions in kg(NO2)/kg(air)/s
          ! Note: lightningems        ==> L-NOx emissions
          !       lightning_emissions ==> type flag used in asad diag-routines
          !                               (c.f., asad_chem_flugs_diags.F90)
          IF (L_asad_use_chem_diags .AND. L_asad_use_light_ems) &
               CALL asad_3D_emissions_diagnostics( &
               row_length, &
               rows, &
               model_levels, &
               inox, &
               lightningems, &
               surf_area, &
               totnodens, &
               volume, &
               mass, &
               m_no2, &
               timestep, &
               lightning_emissions, &
               ierr)
         

! Convert aircraft emissions from kg NO2/gridbox/s to kg NO/m2/s
        DO l=1,model_levels
          DO k=1,rows
            DO j=1,row_length
              conv_aircraftems(j,k,l)  = aircraftems(j,k,l)             &
                                       * m_no/(surf_area(j,k)*m_no2) 
            END DO
          END DO
        END DO

! Setup SO2 emissions
        IF(L_ukca_aerchem .OR. L_ukca_achem) THEN

! Check if tracer and emissions are present
          DO k = 1,jpctr
            IF(advt(k) == 'SO2       ') iso2 = k
          ENDDO
          IF (iso2 == -99) THEN
            cmessage = 'Did not find SO2 tracer'
            errcode=1
            CALL EREPORT('UKCA_EMISSION_CTL',errcode,cmessage)
          ENDIF
          IF (jso2_high == -99) THEN
            cmessage = 'Did not find High Level SO2 emissions'
            errcode=1
            CALL EREPORT('UKCA_EMISSION_CTL',errcode,cmessage)
          ENDIF

        ENDIF    ! L_ukca_aerchem/achem

!       Add aircraft emissions to NO or NOx tracer

        DO k = 1,model_levels   ! Loop over levels

          CALL TRSRCE(                                                 &
          rows, row_length, 0, 0, 0, 0,                                &
          model_levels, wet_levels,                                    &
          0, 0,                                                        &
          theta, q , qcl , qcf , exner_rho_levels, rho_r2,             &
          tracers(:,:,k,inox), conv_aircraftems(:,:,k), k,             &
          timestep, 1, 1, 0.0)

        ENDDO                  ! End of looping over levels

          ! diagnostics - aircraft emissions (conv_aircraftems) in kg(NO)/m^2/s
          IF (L_asad_use_chem_diags .AND. L_asad_use_air_ems)           &
               CALL asad_3D_emissions_diagnostics(                      &
               row_length, rows, model_levels,                          &
               inox, conv_aircraftems,  surf_area,                      &
               totnodens, volume, mass,                                 &
               m_no, timestep, aircraft_emissions,                      &
               ierr)



! Update tracer fields with NO/NOx lightning emissions 
        tracers(1:row_length,1:rows,:,inox) =                           &
        tracers(1:row_length,1:rows,:,inox) +                           &
        timestep*lightningems*m_no/m_no2 

! Add 3-D volcanic and high-level anthropogenic emissions to SO2 tracer
!  remove direct sulphate fraction of emissions and convert from [S] to [SO2]
        IF(L_ukca_aerchem .OR. L_ukca_achem) THEN
          DO k = 1,model_levels
            IF (k == SO2_high_level) THEN
              conv_SO2emiss_3D(:,:,k) = (SO2_volc_3D(:,:,k) +           &
              all_emissions(:,:,jso2_high))*(1.0 - mode_parfrac/100.0)* &
                  m_so2/m_s 
            ELSE
              conv_SO2emiss_3D(:,:,k) = SO2_volc_3D(:,:,k)*             &
                               (1.0 - mode_parfrac/100.0)*m_so2/m_s
            ENDIF
            CALL TRSRCE( rows, row_length, 0, 0, 0, 0,                 &
             model_levels, wet_levels, 0, 0,                           &
             theta, q , qcl , qcf , exner_rho_levels, rho_r2,          &
             tracers(:,:,k,iso2), conv_SO2emiss_3D(:,:,k),             &
             k, timestep, 1, 1, 0.0)
          ENDDO
! Asad emission diagnostics
         IF (L_asad_use_chem_diags .AND. L_asad_use_sulp_ems) THEN
           CALL ASAD_3D_EMISSIONS_DIAGNOSTICS(row_length, rows,         &
                        model_levels, iso2,                             &
                        conv_SO2emiss_3D, surf_area, totnodens,         &
                        volume, mass,                                   &
                        m_so2, timestep,                                &
                        so2_emchar, icode)
         END IF
        ENDIF

      IF ((L_ukca_strat .OR. L_ukca_stratcfc .OR.                       &
           L_ukca_strattrop) .AND. iso2 > 0)  THEN
! Perform emission of volcanic SO2 from explosive volcanic eruptions into
!  Stratosphere

         ! diagnostics
        IF (L_asad_use_chem_diags .AND. L_asad_use_volc_ems) THEN
          ALLOCATE(tmp3dems(row_length, rows, model_levels))
          tmp3dems(:,:,:)=                                              &
                 tracers(1:row_length, 1:rows, 1:model_levels, iso2)
        END IF


! DEPENDS ON: ukca_volcanic_so2
        CALL UKCA_VOLCANIC_SO2(                                         &
                 tracers(1:row_length, 1:rows, :, iso2),                &
                 mass, theta_latitude, true_longitude,                  &
                 delta_phi, delta_lambda,                               &
                 row_length, rows, model_levels,                        &
                 iyear, imonth, iday, timestep,                         &
                 tropopause_height,                                     &
                 r_theta_levels(1:row_length, 1:rows, 1:model_levels))

! Diagnostics - volcanic emissions in kg(SO2)/kg(air)/gridcell/timestep
        IF (L_asad_use_chem_diags .AND. L_asad_use_volc_ems) THEN
          tmp3dems(:,:,:) =                                             &
                tracers(1:row_length, 1:rows, :, iso2) - tmp3dems(:,:,:)

          CALL asad_3D_emissions_diagnostics(                           &
                row_length, rows, model_levels,                         &
                iso2, tmp3dems, surf_area,                              &
                totnodens, volume, mass,                                &
                m_so2, timestep, volcanic_emissions,                    &
                ierr)
          DEALLOCATE(tmp3dems)
        END IF

      END IF    ! L_ukca_strat etc.



      IF (lhook) CALL dr_hook('UKCA_EMISSION_CTL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_EMISSION_CTL
