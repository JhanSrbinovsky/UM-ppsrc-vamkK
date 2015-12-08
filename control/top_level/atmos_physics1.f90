! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interface to Atmos Physics parametrizations before S-L advection.
!
! Subroutine Interface:
      SUBROUTINE Atmos_Physics1(                                        &

! Parallel variables
     &  global_row_length, global_rows, n_proc, n_procx, n_procy        &
     &, g_rows, g_row_length                                            &

! model dimensions
     &, row_length, rows, n_rows, land_points, model_levels, wet_levels &
     &, bl_levels, dst_levels, dsm_levels, Ozone_levels, cloud_levels   &
     &, land_ice_points, soil_points, n_cca_levels, ntiles, nice_use    &
     &, salt_dim1, salt_dim2, salt_dim3, tr_levels, tr_ukca             &
      , cdnc_dim1, cdnc_dim2, cdnc_dim3                                 &
     &, co2_dim_len, co2_dim_row, co2_dim_lev                           &
     &, n_arcl_species, n_arcl_compnts, i_arcl_compnts                  &

! model switches
     &, L_regular, L_lbc_old                                            &
     &, L_Rad_Step, L_Rad_Step_diag, L_Rad_Step_prog                    &
     &, L_CAL360, Ltimer                                                &
     &, L_mixing_ratio                                                  &
     &, L_ukca_chem, L_ukca_useumuivals, L_ukca_set_trace_gases         &
     &, L_ukca_strat, L_ukca_strattrop                                  &
     &, L_ukca_prescribech4, L_USE_ARCL                                 &

! model Parameters
     &, rhcrit                                                          &
     &, min_trop_level, max_trop_level                                  &

! Position of greenhouse gases in tracer_ukca array 
     &, ngrgas, grgas_addr, Ntot_land, Ntot_sea                         &
! parameter for stochastic physics random parameters
     &, M_CI                                                            &

! in coordinate information
     &, delta_lambda, delta_phi, lat_rot_NP, long_rot_NP                &

! in time stepping information.
     &, val_year, val_day_number, val_hour, val_minute                  &
     &, val_second, PREVIOUS_TIME,                                      &

! diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork1,STASHwork2,STASHwork4,STASHwork6,STASHwork14          &
!
! Additional variables for SCM diagnostics
     &, nSCMDpkgs, L_SCMDiags                                           &
!
! in data fields.
     &, theta, q, qcl, qcf, qcf2, qrain, qgraup, rho, u, v, p, p_star   &
     &, exner_rho_levels, exner_theta_levels                            &
     &, land_sea_mask, p_theta_levels,fland_ctile                       &

     &, frac_control                                                    &
      , ukca_cdnc                                                       &
! ancillary fields and fields needed to be kept from timestep to
! timestep

     &, land_index,rgrain,soot,ntml,cumulus                             &
     &, ice_fract,ice_fract_cat,ice_thick_cat                           &
     &, cca,ccb,cct,cclwp                                               &
     &, ccw,lcbase                                                      &
     &, t_surf, tstar_land_ctile, tstar_sea_ctile, tstar_sice_ctile     &
     &, sice_alb_ctile,land_alb_ctile,snow_depth,snow_depth_sea_cat     &
     &, ozone, SW_incs, LW_incs, dirpar_inc                             &
     &, O3_trop_level, O3_trop_height, T_trop_level, T_trop_height      &
     &, zh, sd_orog_land , orog_grad_xx_land, orog_grad_xy_land         &
     &, orog_grad_yy_land, area_cloud_fraction, cf, cfl, cff            &
     &, aerosol_em, arcl                                                &
      , albsoil, albobs_sw, albobs_vis, albobs_nir, lai, snow_tile      &
      , tile_frac, tstar_tile, z0_tile                                  &
     &, dOLR_rts, LW_down, SW_tile_rts, ES_SPACE_INTERP, RAD_MASK       &
     &, cos_zenith_angle                                                & 

! Variables for COSP
      , cosp_crain_3d,cosp_csnow_3d                                     &

! IN/OUT JULES 2 prognostics
      , snowdepth_p,lake_h_ice_p                                        &

! in/out
     &, theta_star, q_star, qcl_star, qcf_star, qcf2_star, qrain_star   &
     &, qgraup_star, cf_star, cfl_star, cff_star, u_inc, v_inc          &
     &, energy_correction, sum_eng_fluxes, sum_moist_flux, AEROSOL      &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, SO2, SO4_AITKEN, SO4_ACCU, SO4_DISS, NH3                        &
     &, soot_new, soot_agd, soot_cld, bmass_new, bmass_agd, bmass_cld   &
     &, ocff_new, ocff_agd, ocff_cld, nitr_acc, nitr_diss               &
     &, co2, tracer_ukca, BIOGENIC, ASTEPS_SINCE_TRIFFID, ukca_radaer   &
! out fields
     &, ls_rain, ls_snow, micro_tends, unscaled_dry_rho                 &
      , photosynth_act_rad, rad_hr, dOLR, SW_tile                       &

! error information
     &, Error_code  )

      USE dynamics_input_mod, ONLY: l_endgame
      USE dynamics_grid_mod,  ONLY: l_vatpoles

      USE atmos_constants_mod, ONLY: kappa, p_zero, cp
     
      USE water_constants_mod, ONLY: lc, lf, tfs

      USE atm_fields_bounds_mod
      
      USE cv_run_mod, ONLY: l_pc2_diag_sh

      USE cloud_inputs_mod,  ONLY: l_micro_eros,  l_rhcpt,              &
                                   pc2_falliceshear_method, l_pc2
      USE pc2_constants_mod, ONLY: real_shear

      USE ukca_trace_gas_mixratio, ONLY: ukca_set_trace_gas_mixratio
      
      USE timestep_mod, ONLY:                                           &
          timestep_number, timestep, radiation_timestep,                &
          radiation_tstep_diag, radiation_tstep_prog

      USE rad_input_mod
      USE clmchfcg_scenario_mod, ONLY: clim_fcg_nyears, clim_fcg_years, &
                                       clim_fcg_levls, clim_fcg_rates,  &
           s_co2, s_n2o, s_ch4, s_cfc11, s_cfc12, s_cfc113, s_cfc114,   &
           s_hcfc22, s_hfc125, s_hfc134a, lenscen
!        histories/scenarios of climate change forcings

      USE arcl_mod, ONLY: npd_arcl_species, npd_arcl_compnts

      USE NI_gwd_ctl_mod, ONLY : NI_gwd_ctl

      USE earth_constants_mod, ONLY: g, earth_radius 
      
      USE trignometric_mod, ONLY:                                       &
          cos_theta_latitude, sin_theta_latitude,                       &
          cos_theta_longitude, sin_theta_longitude,                     &
          FV_cos_theta_latitude, true_latitude, true_longitude

      USE global_2d_sums_mod, ONLY: global_2d_sums

      USE g_wave_input_mod, ONLY: l_gwd, l_use_ussp
      
!     JULES
      USE switches, ONLY :                                              &
        l_flake_model, can_rad_mod, l_snow_albedo, l_ssice_albedo,      &
        l_sice_meltponds, l_sice_scattering, l_sice_hadgem1a,           &
        l_cice_alb, l_ctile
      USE prognostics, ONLY :                                           &
        snowdepth
      USE lake_mod, ONLY :                                              &
        lake_h_ice

      USE ancil_info, ONLY:                                             &
          ssi_pts, sea_pts, sice_pts                                    &
        , sice_pts_ncat, ssi_index, sea_index, sice_index               &
        , sice_index_ncat, fssi, sea_frac, sice_frac, sice_frac_ncat

      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE ukca_radaer_struct_mod
      USE swapable_field_mod,    ONLY: swapable_field_pointer_type

      USE level_heights_mod, ONLY:                                      &
          eta_theta_levels, eta_rho_levels,                             &
          r_theta_levels, r_rho_levels
      
      USE nstypes
      USE Submodel_Mod

      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE domain_params
      USE cosp_types_mod, ONLY: cosp_gridbox, cosp_config
      USE cosp_init_mod
      USE cosp_main_mod
      USE ukca_feedback_mod, ONLY: p_o3, p_ch4, p_n2o, p_f11, p_f12,    &
                                   p_f113, p_f22, p_h2os

      USE uv_p_pnts_mod,   ONLY: uv_p_pnts
      USE set_seasalt_mod, ONLY: set_seasalt_4A


      USE add_eng_corr_mod, ONLY: add_eng_corr
      USE dust_parameters_mod, ONLY: l_dust
      USE um_input_control_mod,  ONLY:                                  &
           model_domain,                                                &
                               l_sulpc_so2,      l_sulpc_nh3,           &
           l_soot,             l_biomass,                               &
           l_co2_interactive,  l_triffid,        l_use_seasalt_autoconv,&
           l_ocff,             l_nitrate,                               &
                                                 l_use_methox,          &
           l_mr_physics1,      l_mr_physics2,    h_sect,                &
                               l_consistent_cdnc
      USE ukca_option_mod, ONLY: l_ukca, l_ukca_aie1, l_ukca_radaer
      USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain, l_rain
      USE eng_corr_inputs_mod, ONLY: l_emcorr, lflux_reset
      USE murk_inputs_mod, ONLY: l_murk, l_murk_source, l_murk_bdry, l_murk_rad

      USE cosp_input_mod, ONLY: l_cosp

      USE diagnostics_pc2checks_mod, ONLY: diagnostics_pc2checks
      USE u_to_p_mod, ONLY: u_to_p
      USE v_to_p_mod, ONLY: v_to_p
      USE trsrce_mod, ONLY: trsrce
      IMPLICIT NONE
!
! Description: This version interfaces to physics schemes in sequence:
!    energy correction              (optional)
!    microphysics (cloud and large scale precipitation schemes)
!    radiation
!    gravity wave drag
!
!          CALLed before Semi-Lagrangian in atmosphere timestep.
! Method:
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! System component: Control Atmos
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
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

! Subroutine arguments ===============================================

! arguments with intent in. ie: input variables.

! Parallel setup variables
      Integer                                                           &
     &  global_row_length                                               &
                           ! number of points on a row
     &, global_rows                                                     &
                           ! NUMBER OF global rows
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       


! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, salt_dim1                                                       &
                      !
     &, salt_dim2                                                       &
                      ! Dimensions of sea-salt aerosol arrays
     &, salt_dim3                                                       &
                      !
     &, land_points                                                     &
                    ! IN No.of land points being processed, can be 0.
     &, model_levels                                                    &
     &, wet_levels                                                      &
     &, bl_levels                                                       &
     &, dst_levels                                                      &
                    ! number of deep soil temperature levels
     &, dsm_levels                                                      &
                    ! number of deep soil moisture levels
     &, Ozone_levels                                                    &
     &, cloud_levels                                                    &
     &, land_ice_points                                                 &
                        ! number of land ice points
     &, soil_points                                                     &
                        ! number of soil points
     &, n_cca_levels                                                    &
                      ! Number of levels for conv cloud
     &, tr_levels                                                       &
                      ! number of tracer levels
     &, tr_ukca                                                         &
                      ! number of UKCA tracers
                      ! amount: 1 for 2D, nlevs for 3D.
     &, ntiles                                                          &
     &, nice_use                                                        &
                      ! No. of sea ice catagories used in radiation
     &,CO2_DIM_LEN                                                      &
                     !\ For dimension 3-D CO2 field to be passed
     &,CO2_DIM_ROW                                                      &
                     !/ to NI_rad_ctl
     &,CO2_DIM_LEV   !/

      Logical                                                           &
     &  Ltimer   ! true then output some timing information

      Logical                                                           &
     &  L_Rad_Step                                                      &
                        ! true if a radiation timestep
     &, L_Rad_Step_diag                                                 &
                         ! true if fast radiation timestep    (3C)
     &, L_Rad_Step_prog                                                 &
                         ! true if slow radiation timestep    (3C)
     &, L_CAL360                                                        &
                        ! true if using 360 day calender
     &, L_regular                                                       &
                     ! True if NOT variable resolution
     &, L_lbc_old                                                       &
                        !  false for new lbc treatment
     &, L_mixing_ratio                                                  &
                       ! true if mixing ratios used
     &, L_ukca_chem                                                     &
                ! T for ukca chemistry
     &, L_ukca_useumuivals                                              &
                ! Switch for prescribing CFC/CH4/N2O concs in UKCA
     &, L_ukca_set_trace_gases                                          &
                ! Switch for prescribing atmospheric gases in UKCA
     &, L_ukca_strat                                                    &
                ! True for stratospheric chemistry scheme
     &, L_ukca_strattrop                                                &
                ! True for stratospheric+tropospheric chemistry scheme
     &, L_ukca_prescribech4                                             
                ! Switch for prescribing surface ch4 in UKCA


      Real                                                              &
                        !, Intent(IN)
     &    Ntot_land                                                     &
                               ! Number of droplets over land / m-3
     &   ,Ntot_sea             ! Number of droplets over sea / m-3

      Real                                                              &
                        !, Intent(IN)
     &    M_CI      ! variable to modify ice fall speed for LSPCON3C
                    !  for stochastic physics random parameters     

      Integer                                                           &
     &  min_trop_level                                                  &
                        ! Lowest permitted level for the tropopause
!                       ! used for radiative purposes.
     &, max_trop_level  ! Highest permitted level for the tropopause
!                       ! used for radiative purposes.


      Real                                                              &
     &  rhcrit(wet_levels)   ! IN Critical relative humidity.
                             ! the values need to be tuned
                             ! for the given set of levels.

! Diagnostics info
      REAL                                                              &
     & STASHWORK1(*)                                                    &
                         ! STASH workspace for section 1 (SW rad)
     &,STASHWORK2(*)                                                    &
                         ! STASH workspace for section 2 (LW rad)
     &,STASHWORK4(*)                                                    &
                         ! STASH workspace for section 4 (LS precip)
     &,STASHWORK6(*)                                                    &
                         ! STASH workspace for section 6 (gw drag)
     &,STASHWORK14(*)    ! STASH workspace for section 14 (Energy cor)


! Data arrays

! Small halo
      REAL, INTENT (INOUT) :: u  (udims_s%i_start:udims_s%i_end,        &
                                  udims_s%j_start:udims_s%j_end,        &
                                  udims_s%k_start:udims_s%k_end)
      REAL, INTENT (INOUT) :: v  (vdims_s%i_start:vdims_s%i_end,        &
                                  vdims_s%j_start:vdims_s%j_end,        &
                                  vdims_s%k_start:vdims_s%k_end)
      REAL, INTENT (INOUT) :: rho(pdims_s%i_start:pdims_s%i_end,        &
                                  pdims_s%j_start:pdims_s%j_end,        &
                                  pdims_s%k_start:pdims_s%k_end)  
      REAL, INTENT (INOUT) :: p  (pdims_s%i_start:pdims_s%i_end,        &
                                  pdims_s%j_start:pdims_s%j_end,        &
                                  pdims_s%k_start:pdims_s%k_end) 
      REAL, INTENT (INOUT) :: p_theta_levels                            &  
                                 (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT (INOUT) :: theta                                     &
                                 (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT (INOUT) :: exner_rho_levels                          &
                                 (pdims_s%i_start:pdims_s%i_end,        &
                                  pdims_s%j_start:pdims_s%j_end,        &
                                  pdims_s%k_start:pdims_s%k_end + 1)      
      REAL, INTENT (INOUT) :: exner_theta_levels                        &
                                 (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end) 

!No halo
      REAL, INTENT (INOUT) :: p_star (pdims%i_start:pdims%i_end,        &
                                      pdims%j_start:pdims%j_end)

!Large halo
      REAL, INTENT (INOUT) :: q  (qdims_l%i_start:qdims_l%i_end,        &
                                  qdims_l%j_start:qdims_l%j_end,        &
                                  qdims_l%k_start:qdims_l%k_end)
      REAL, INTENT (INOUT) :: qcl(qdims_l%i_start:qdims_l%i_end,        &
                                  qdims_l%j_start:qdims_l%j_end,        &
                                  qdims_l%k_start:qdims_l%k_end)       
      REAL, INTENT (INOUT) :: qcf(qdims_l%i_start:qdims_l%i_end,        &
                                  qdims_l%j_start:qdims_l%j_end,        &
                                  qdims_l%k_start:qdims_l%k_end)
      REAL, INTENT (INOUT) :: qcf2                                      &
                                 (qdims_l%i_start:qdims_l%i_end,        &
                                  qdims_l%j_start:qdims_l%j_end,        &
                                  qdims_l%k_start:qdims_l%k_end)
      REAL, INTENT (INOUT) :: qrain                                     &
                                 (qdims_l%i_start:qdims_l%i_end,        &
                                  qdims_l%j_start:qdims_l%j_end,        &
                                  qdims_l%k_start:qdims_l%k_end)
      REAL, INTENT (INOUT) :: qgraup                                    &
                                 (qdims_l%i_start:qdims_l%i_end,        &
                                  qdims_l%j_start:qdims_l%j_end,        &
                                  qdims_l%k_start:qdims_l%k_end)
      REAL, INTENT (INOUT) :: cf (qdims_l%i_start:qdims_l%i_end,        &
                                  qdims_l%j_start:qdims_l%j_end,        &
                                  qdims_l%k_start:qdims_l%k_end)       
      REAL, INTENT (INOUT) :: cfl(qdims_l%i_start:qdims_l%i_end,        &
                                  qdims_l%j_start:qdims_l%j_end,        &
                                  qdims_l%k_start:qdims_l%k_end)
      REAL, INTENT (INOUT) :: cff(qdims_l%i_start:qdims_l%i_end,        &
                                  qdims_l%j_start:qdims_l%j_end,        &
                                  qdims_l%k_start:qdims_l%k_end)

      Real                                                              &
     &  energy_correction

! ancillary arrays and fields required to be saved from timestep to
! timestep.

      Real                                                              &
     &  T_surf(row_length, rows)
! The following are only used if coastal tiling is switched on:
      REAL                                                              &
     & FLAND_CTILE(LAND_POINTS)                                         &
                                   ! IN Land fraction on land points.
     &,TSTAR_LAND_CTILE(ROW_LENGTH,ROWS)                                &
!                                  ! IN Land mean sfc temperature (K)
     &,TSTAR_SEA_CTILE(ROW_LENGTH,ROWS)                                 &
!                                  ! IN Open sea sfc temperature (K).
     &,TSTAR_SICE_CTILE(ROW_LENGTH,ROWS,NICE_USE)                       &
!                                  ! IN Sea-ice sfc temperature (K).
     &,LAND_ALB_CTILE(ROW_LENGTH,ROWS)                                  &
!                                  ! INOUT Mean land albedo.
     &,SICE_ALB_CTILE(ROW_LENGTH,ROWS)
!                                  ! INOUT Sea-ice albedo.

      logical                                                           &
     &  land_sea_mask(row_length, rows)                                 &
     &, RAD_MASK(row_length, rows)
!  A mask which ensures a chequerboard pattern of radiation calculations
!  over the whole domain (not just one PE)

      Integer                                                           &
     &  land_index (land_points)                                        &
                                      ! set from land_sea_mask
     &, ntml(row_length, rows)

      Logical                                                           &
     &  cumulus (row_length, rows) ! bl convection flag

      Real                                                              &
     &  snow_depth (row_length, rows)                                   &
                                      ! snow/qrclim.snow.(month)
     &, snow_depth_sea_cat (row_length, rows, nice_use)                 &
                                      ! snow depth on sea ice
     &, sd_orog_land (land_points)                                      &
                                   ! orog/qrparm.orog.stdev
     &, orog_grad_xx_land(land_points)                                  &
                                       ! orog/qrparm.orog.sigmaxx
     &, orog_grad_xy_land(land_points)                                  &
                                       ! orog/qrparm.orog.sigmaxy
     &, orog_grad_yy_land(land_points)                                  &
                                       ! orog/qrparm.orog.sigmayy
     &, zh(row_length, rows)                                            &
                             ! boundary layer height
     &, aerosol_em( 1 : row_length, 1 : rows, 1 :model_levels )

      REAL ::                                                           &
        snowdepth_p(land_points,ntiles)                                 &
                                  ! Snow depth on ground on tiles (m)
      , lake_h_ice_p(land_points)
                                  ! FLake lake ice thickness (m)

      Real                                                              &
     &  ozone(row_length, rows, ozone_levels)                           &
     &, O3_trop_level(row_length,rows)                                  &
     &, O3_trop_height(row_length,rows)                                 &
     &, cos_zenith_angle(row_length, rows)                              &
     &, T_trop_level(row_length,rows)                                   &
     &, T_trop_height(row_length,rows)                                  &
     &, SW_incs(row_length, rows, 0:model_levels+1)                     &
     &, LW_incs(row_length, rows, 0:model_levels)                       &
     &, dirpar_inc(row_length, rows)                                    &
     &, soot(row_length, rows)                                          &
                                       ! Snow soot
     &, ice_fract (row_length, rows)                                    &
                                       ! ice/qrclim.ice.(month)
     &, ice_fract_cat (row_length, rows, nice_use)                      &
                                       ! ice fraction per cat
     &, ice_thick_cat (row_length, rows, nice_use)
                                       ! effective ice thickness per cat

      Real                                                              &
     &  albsoil(land_points)                                            &
      , albobs_sw(land_points)                                          & 
      , albobs_vis(land_points)                                         & 
      , albobs_nir(land_points)                                         &
     &, lai(land_points, npft)                                          &
     &, rgrain(land_points, ntiles)                                     &
     &, snow_tile(land_points, ntiles)                                  &
     &, tile_frac(land_points, ntype)                                   &
     &, tstar_tile(land_points, ntiles)                                 &
     &, z0_tile(land_points, ntiles)                                    &
     &, dOLR_rts(row_length, rows)                                      &
                                         ! TOA - surface upward LW
     &, LW_down(row_length, rows)                                       &
                                         ! Surface downward LW
     &, SW_tile_rts(land_points, ntiles)                                &
                                         ! Surface net SW on land tiles
     &, ES_SPACE_INTERP(4, row_length, rows)
!              Coeffs for spatial interpolation of radiation quantities
!

! Aerosol climatology for NWP
      
      ! Number of requested species within the climatology
      Integer n_arcl_species
      
      ! Corresponding number of requested components
      Integer n_arcl_compnts
      
      ! Model switch for each species
      Logical L_USE_ARCL(NPD_ARCL_SPECIES)
      
      ! Mass-mixing ratios
      Real                                                              &
     &   arcl(row_length, rows, model_levels, n_arcl_compnts)
     
      ! Array index of each component
       Integer i_arcl_compnts(NPD_ARCL_COMPNTS)

! Convection Scheme

      Real    :: ccw(row_length, rows, wet_levels)
      Integer :: lcbase(row_length, rows)
 
      Real                                                              &
     &  cca (row_length, rows, n_cca_levels)                            &
     &, cclwp(row_length, rows) ! condensed water path (KG/M**2)

      Integer                                                           &
     &  ccb (row_length, rows)                                          &
     &, cct (row_length, rows)

! Co-ordinate arrays

      REAL ::                                                           &
     &  delta_lambda                                                    &
     &, delta_phi

! time information for current timestep
      Integer                                                           &
     &  val_year                                                        &
     &, val_day_number                                                  &
     &, val_hour                                                        &
     &, val_minute                                                      &
     &, val_second                                                      &
     &,         PREVIOUS_TIME(7)                                        &
     &,     ASTEPS_SINCE_TRIFFID   ! INOUT  Number of atmospheric
!                                  !        timesteps since last call
!                                  !        to TRIFFID.
! UKCA_RADAER structure
      TYPE (ukca_radaer_struct) :: ukca_radaer

! Diagnostic variables
      Real                                                              &
     &  lat_rot_NP                                                      &
     &, long_rot_NP

      Real                                                              &
     &  area_cloud_fraction(row_length, rows, wet_levels)

! Variables for COSP
! 3D convective rainfall flux COSP
      REAL,INTENT(IN) :: cosp_crain_3d(row_length,rows,model_levels)
! 3D convective snowfall flux COSP
      REAL,INTENT(IN) :: cosp_csnow_3d(row_length,rows,model_levels)

! arguments with intent in/out. ie: input variables changed on output.

      Real, Intent (InOut) ::                                           &
     &  sum_eng_fluxes(row_length,rows)                                 &
                                         ! sum atmosphere fluxes
     &, sum_moist_flux(row_length,rows)  ! sum moist fluxes

      REAL, INTENT (INOUT) :: theta_star                                &
                                 (tdims_s%i_start:tdims_s%i_end,        &
                                  tdims_s%j_start:tdims_s%j_end,        &
                                  tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT (INOUT) :: q_star                                    &
                                 (qdims_s%i_start:qdims_s%i_end,        &
                                  qdims_s%j_start:qdims_s%j_end,        &
                                  qdims_s%k_start:qdims_s%k_end)      
      REAL, INTENT (INOUT) :: qcl_star                                  &
                                 (qdims_s%i_start:qdims_s%i_end,        &
                                  qdims_s%j_start:qdims_s%j_end,        &
                                  qdims_s%k_start:qdims_s%k_end)      
      REAL, INTENT (INOUT) :: qcf_star                                  &
                                 (qdims_s%i_start:qdims_s%i_end,        &
                                  qdims_s%j_start:qdims_s%j_end,        &
                                  qdims_s%k_start:qdims_s%k_end)
      REAL, INTENT (INOUT) :: qcf2_star                                 &
                                 (qdims_s%i_start:qdims_s%i_end,        &
                                  qdims_s%j_start:qdims_s%j_end,        &
                                  qdims_s%k_start:qdims_s%k_end) 
      REAL, INTENT (INOUT) :: qrain_star                                &
                                 (qdims_s%i_start:qdims_s%i_end,        &
                                  qdims_s%j_start:qdims_s%j_end,        &
                                  qdims_s%k_start:qdims_s%k_end) 
      REAL, INTENT (INOUT) :: qgraup_star                               &
                                 (qdims_s%i_start:qdims_s%i_end,        &
                                  qdims_s%j_start:qdims_s%j_end,        &
                                  qdims_s%k_start:qdims_s%k_end) 
      REAL, INTENT (INOUT) :: cf_star                                   &
                                 (qdims_s%i_start:qdims_s%i_end,        &
                                  qdims_s%j_start:qdims_s%j_end,        &
                                  qdims_s%k_start:qdims_s%k_end)       
      REAL, INTENT (INOUT) :: cfl_star                                  &
                                 (qdims_s%i_start:qdims_s%i_end,        &
                                  qdims_s%j_start:qdims_s%j_end,        &
                                  qdims_s%k_start:qdims_s%k_end)
      REAL, INTENT (INOUT) :: cff_star                                  &
                                 (qdims_s%i_start:qdims_s%i_end,        &
                                  qdims_s%j_start:qdims_s%j_end,        &
                                  qdims_s%k_start:qdims_s%k_end)
      REAL, INTENT (INOUT) :: u_inc                                     &
                                 (udims_s%i_start:udims_s%i_end,        &
                                  udims_s%j_start:udims_s%j_end,        &
                                  udims_s%k_start:udims_s%k_end)
      REAL, INTENT (INOUT) :: v_inc                                     &
                                 (vdims_s%i_start:vdims_s%i_end,        &
                                  vdims_s%j_start:vdims_s%j_end,        &
                                  vdims_s%k_start:vdims_s%k_end)   

      REAL, INTENT(INOUT) ::  aerosol  (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 

      REAL, INTENT(INOUT) ::  DUST_DIV1(tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 

      REAL, INTENT(INOUT) ::  DUST_DIV2(tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 

      REAL, INTENT(INOUT) ::  DUST_DIV3(tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 

      REAL, INTENT(INOUT) ::  DUST_DIV4(tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 

      REAL, INTENT(INOUT) ::  DUST_DIV5(tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 

      REAL, INTENT(INOUT) ::  DUST_DIV6(tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 

      REAL, INTENT(INOUT) ::  so2      (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  so4_aitken                                &
                                       (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  so4_accu (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  so4_diss (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  nh3      (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  soot_new (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  soot_agd (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  soot_cld (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  bmass_new(tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  bmass_agd(tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  bmass_cld(tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  ocff_new (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  ocff_agd (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  ocff_cld (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  nitr_acc (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  nitr_diss(tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  co2      (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end) 
      REAL, INTENT(INOUT) ::  tracer_ukca                               &
                                       (tdims_s%i_start:tdims_s%i_end,  &
                                        tdims_s%j_start:tdims_s%j_end,  &
                                        tdims_s%k_start:tdims_s%k_end,  &
                                        tr_ukca) 
     
      Real, Intent(In) ::                                               &
     &  biogenic(row_length, rows, model_levels)

! arguments with intent out. ie: output variables.

! arrays holding information to be passed between physics
! routines.

      Real                                                              &
     &  ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, micro_tends(row_length, rows, 2, bl_levels)
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)

! Radiation fields 1. SW & common with LW.
      Real                                                              &
     &  photosynth_act_rad(row_length, rows)                            &
                                             ! Net downward
!                                 shortwave radiation in band 1 (w/m2).
     &, rad_hr(row_length, rows, 2, bl_levels)                          &
                                               !
!                               ! BL radiative (LW,SW) heating rates
     &, dOLR(row_length, rows)                                          &
                                    ! TOA - surface upward LW
     &, SW_tile(land_points, ntiles)! Surface net SW on land tiles

! Position of greenhouse gases in tracer_ukca array
      INTEGER, INTENT(IN) :: ngrgas
      INTEGER, INTENT(IN) :: grgas_addr(ngrgas)

! UKCA cloud drop number concentration dimensions
      INTEGER, INTENT(IN) :: cdnc_dim1
      INTEGER, INTENT(IN) :: cdnc_dim2
      INTEGER, INTENT(IN) :: cdnc_dim3

!  Additional variables for SCM diagnostics (dummy in full UM)
      Integer                                                           &
     &  nSCMDpkgs             ! No of SCM diagnostics packages

      Logical                                                           &
     &  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

      CHARACTER(LEN=72)                                                      &
     &   cmessage
      Integer                                                           &
     &  Error_code

! === End of arguments ==============================================

! Local parameters:
      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='Atm_Physics1')
      Real                                                              &
     &  amp                                                             &
                        ! amplitude of diurnal variation in tracers
     &, tau_decay                                                       &
                        ! time constant for decay of tracer to clim
     &, clim_murk_land                                                  &
                        ! climatological murk value over land points
     &, clim_murk_sea   ! climatological murk value over sea points
      Parameter (                                                       &
     &            amp=0.7                                               &
     &,           tau_decay=1.728E5                                     &
     &,           clim_murk_land=25.0                                   &
     &,           clim_murk_sea=12.5                                    &
     &           )

! Local scalars:
      REAL :: z0_sea    ! Roughness length over sea (m) 
      REAL :: P1 
      REAL :: windspeed_1(salt_dim1, salt_dim2)
      REAL :: windspeed_10m(salt_dim1, salt_dim2)

      REAL ::                                                           &
        height_theta(salt_dim1, salt_dim2, salt_dim3),                  &
!          theta level centre height above surface
        height_rho_1(salt_dim1, salt_dim2)
!          theta level 1 rho centre height above surface

! loop counters
      Integer                                                           &
     &  i, j, k                                                         &
     &, l, n

! local variables
      Integer                                                           &
     &  rhc_row_length                                                  &
                        ! Row length for RHcrit array
     &, rhc_rows                                                        &
                        ! Row number for RHcrit array
     &, lspice_dim1,lspice_dim2,lspice_dim3                             &
!                       ! Required for array dimensions MCR_CTL2
     &, sulp_dim1                                                       &
                                ! dimensions for sulphate arrays in
     &, DUST_DIM1                                                       &
                                ! dimensions for mineral dust arrays i
     &, DUST_DIM2                                                       &
                                ! in rad_ctl
     &, BIOGENIC_DIM1                                                   &
                                ! dimensions of biogenic arrays in
     &, BIOGENIC_DIM2                                                   &
                                !
     &, sulp_dim2                                                       &
                                !   RAD_CTL
     &, soot_dim1, soot_dim2                                            &
                                ! dimensions of soot arrays in RAD_CTL
     &, bmass_dim1, bmass_dim2                                          &
                                ! dimensions of biomass arrays in radn
     &, ocff_dim1, ocff_dim2                                            &
                                ! dimensions of OCFF arrays in radiation
     &, nitrate_dim1, nitrate_dim2                                      &
                                ! dimensions of nitrate arrays in radn
     &, arcl_dim1, arcl_dim2                                            &
                     ! dimensions of aerosol clim for NWP arrays in radn
     &, ukca_dim1, ukca_dim2                                            &
                     ! dimensions of UKCA_RADAER arrays in RAD_CTL
     &, i_start                                                         &
                        ! Row start point for polar row tidying
     &, istat           ! Status (error code) indicator

! Local data arrays

      Real                                                              &
     &  dirpar_local(row_length, rows)                                  &
     &, u_1(salt_dim1, salt_dim2)                                       &
     &, v_1(salt_dim1, salt_dim2)                                       &
     &, co2_3D(co2_dim_len, co2_dim_row, co2_dim_lev)                   &
     &, u_1_mean                                                        &
     &, v_1_mean                                                        &
     &, frac_control(land_points,ntype)

      REAL :: u_1_arr(1), v_1_arr(1)

      REAL :: T_n   (tdims%i_start:tdims%i_end,        &
                     tdims%j_start:tdims%j_end,        &
                                 1:tdims%k_end) 
      REAL :: q_n   (qdims%i_start:qdims%i_end,        &
                     qdims%j_start:qdims%j_end,        &
                                 1:qdims%k_end)  
      REAL :: qcl_n (qdims%i_start:qdims%i_end,        &
                     qdims%j_start:qdims%j_end,        &
                                 1:qdims%k_end)       
      REAL :: qcf_n (qdims%i_start:qdims%i_end,        &
                     qdims%j_start:qdims%j_end,        &
                                 1:qdims%k_end)       
      REAL :: cf_n  (qdims%i_start:qdims%i_end,        &
                     qdims%j_start:qdims%j_end,        &
                                 1:qdims%k_end)
      REAL :: cfl_n (qdims%i_start:qdims%i_end,        &
                     qdims%j_start:qdims%j_end,        &
                                 1:qdims%k_end)
      REAL :: cff_n (qdims%i_start:qdims%i_end,        &
                     qdims%j_start:qdims%j_end,        &
                                 1:qdims%k_end)
      REAL :: u_on_p(pdims%i_start:pdims%i_end,            &
                     pdims%j_start:pdims%j_end,            & 
                     pdims%k_start:pdims%k_end) 
      REAL :: v_on_p(pdims%i_start:pdims%i_end,            & 
                     pdims%j_start:pdims%j_end,            & 
                     pdims%k_start:pdims%k_end)
      REAL :: theta_inc (tdims%i_start:tdims%i_end,        &
                         tdims%j_start:tdims%j_end,        &
                                     1:tdims%k_end) 
      REAL :: q_inc     (qdims%i_start:qdims%i_end,        &
                         qdims%j_start:qdims%j_end,        &
                                     1:qdims%k_end)  
      REAL :: qcl_inc   (qdims%i_start:qdims%i_end,        &
                         qdims%j_start:qdims%j_end,        &
                                     1:qdims%k_end)       
      REAL :: qcf_inc   (qdims%i_start:qdims%i_end,        &
                         qdims%j_start:qdims%j_end,        &
                                     1:qdims%k_end)       
      REAL :: cf_inc    (qdims%i_start:qdims%i_end,        &
                         qdims%j_start:qdims%j_end,        &
                                     1:qdims%k_end)
      REAL :: cfl_inc   (qdims%i_start:qdims%i_end,        &
                         qdims%j_start:qdims%j_end,        &
                                     1:qdims%k_end)
      REAL :: cff_inc   (qdims%i_start:qdims%i_end,        &
                         qdims%j_start:qdims%j_end,        &
                                     1:qdims%k_end)

! Potential droplet number (includes values where cloud not present)
      REAL :: n_drop_pot(qdims%i_start : qdims%i_end,                   &
                         qdims%j_start : qdims%j_end,                   &
                                     1 : qdims%k_end )

! COSP variables
      TYPE(cosp_config)     :: cosp_cfg     ! Configuration options
      TYPE(cosp_gridbox)    :: cosp_gbx     ! Gridbox-mean inputs
      INTEGER :: cosp_npoints
      LOGICAL :: L_cosp_call

! Local Arrays to store additional microphysics fields if in use
      Real, Dimension (:,:,:), Allocatable ::                           &
     &  qcf2_n,   qrain_n,   qgraup_n                                   &
     &, qcf2_inc, qrain_inc, qgraup_inc                                 &
     &, T_inc_diag,  q_inc_diag,   qcl_inc_diag, qcf_inc_diag           &
     &, cf_inc_diag, cfl_inc_diag, cff_inc_diag

      Real                                                              &
     &  p_layer_boundaries(row_length, rows, 0:model_levels)            &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, p_layer_centres(row_length, rows, 0:model_levels)               &
              ! pressure at layer centres. Same as p_theta_levels
              !except bottom level = pstar, and at top = 0.
     &, exner_layer_boundaries(row_length, rows, 0:model_levels)        &
     &, exner_layer_centres(row_length, rows, 0:model_levels)

      Real :: moist                                                     &
                                 (qdims_l%i_start:qdims_l%i_end,        &
                                  qdims_l%j_start:qdims_l%j_end,        &
                                  qdims_l%k_start:qdims_l%k_end)

      ! holds total moisture for conversion from wet to dry density
      REAL :: unscaled_dry_rho(pdims_s%i_start:pdims_s%i_end,  &
                               pdims_s%j_start:pdims_s%j_end,  &
                               pdims_s%k_start:pdims_s%k_end) 

      ! unscaled dry density
      Real :: weight1
      Real :: weight2
      Real :: weight3
      Real :: temp

      REAL, TARGET, ALLOCATABLE ::                                      &
        ext_p_layer_centres(:,:,:),                                     &
        ext_tl(:,:,:),                                                  &
        ext_ql(:,:,:),                                                  &
        ext_qcf(:,:,:),                                                 &
        ext_ice_frac(:,:),                                              &
        ext_land_frac(:,:)

      INTEGER :: i_field   ! field counter for multivar swapbounds
      TYPE(swapable_field_pointer_type) :: fields_to_swap(4) 
                           ! mv swapbounds

      REAL                                                              &
     & FLAND(LAND_POINTS)                                               &
                                   ! Land fraction on land points.
     &,TSTAR_LAND(ROW_LENGTH,ROWS)                                      &
                                   ! Land mean sfc temperature (K)
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                   ! Open sea sfc temperature (K).
     &,TSTAR_SICE_CAT(ROW_LENGTH,ROWS,NICE_USE)                         &
                                   ! Sea-ice sfc temperature (K).
     &,LAND_ALB(ROW_LENGTH,ROWS)                                        &
                                   ! Mean land albedo.
     &,SICE_ALB(ROW_LENGTH,ROWS)   ! Sea-ice albedo.

      LOGICAL                                                           &
     & land0p5(row_length, rows)

! Diagnostics controlled by Diagnostic switches

! fields output by ls_cld not using stash flag

! Energy correction work variables

      Real                                                              &
     &  tot_precip_scaled(row_length, rows)                             &
     &, lclf

! Temporary logicals
      Logical                                                           &
     &  L_zero_boundaries                                               &
     &, L_zero_halos                                                    &
     &, L_use_dirpar

      INTEGER                                                           &
     &  level

      REAL, ALLOCATABLE :: grgas_field(:,:,:,:)

!     Global topmost cloudy level
      INTEGER :: global_cloud_top

!     Level of tropopause
      INTEGER :: trindx(tdims%i_start:tdims%i_end,                      &
                        tdims%j_start:tdims%j_end)

!     SCM Dummy variables to keep call to tropin consistent.
      REAL :: scm_dummy_1d(1,1), scm_dummy_2d(1,1,0:model_levels)

      REAL ::                                                           &
        sea_salt_film(salt_dim1, salt_dim2, salt_dim3),                 &
!          Film-mode sea-salt aerosol number concentration
        sea_salt_jet(salt_dim1, salt_dim2, salt_dim3),                  &
!          Jet-mode sea-salt aerosol number concentration
        height(salt_dim1, salt_dim2, salt_dim3)
!          Layer-centre height above surface

      REAL :: ukca_cdnc(cdnc_dim1, cdnc_dim2, cdnc_dim3) 
!          Cloud droplet number from UKCA-MODE

      REAL :: land_fract(tdims%i_start:tdims%i_end,                     &
                         tdims%j_start:tdims%j_end)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! needed for vatpoles for fv_cos_theta_latitude vs cos_theta_latitude
REAL, POINTER :: xx_cos_theta_latitude (:,:)

      IF (lhook) CALL dr_hook('ATMOS_PHYSICS1',zhook_in,zhook_handle)

IF ( l_vatpoles ) THEN
   xx_cos_theta_latitude => cos_theta_latitude
ELSE
   xx_cos_theta_latitude => fv_cos_theta_latitude
END IF ! vatpoles

! ----------------------------------------------------------------------
! Section INI. Initialisation of variables.
! ----------------------------------------------------------------------
      IF ( l_endgame ) THEN
      z0_sea = 2.5E-04    ! Roughness length over sea (m) 
      P1     = ALOG(10.0/Z0_SEA)
      END IF

      If (L_mcr_qcf2) Then  ! Second cloud ice variable in use
        Allocate ( qcf2_inc(row_length, rows, wet_levels) )
        qcf2_inc(:,:,:) = 0.0
      Else
        Allocate ( qcf2_inc(1,1,1) )
      End If

      If (L_mcr_qrain) Then  ! Prognostic rain in use
        Allocate ( qrain_inc(row_length, rows, wet_levels) )
        qrain_inc(:,:,:) = 0.0
      Else
        Allocate ( qrain_inc(1,1,1) )
      End If

      If (L_mcr_qgraup) Then  ! Prognostic graupel in use
        Allocate ( qgraup_inc(row_length, rows, wet_levels) )
        qgraup_inc(:,:,:) = 0.0
      Else
        Allocate ( qgraup_inc(1,1,1) )
      End If

      IF (L_TRIFFID) THEN
! Increment counter for number of atmosphere timesteps since last
! call to TRIFFID vegetation model
        ASTEPS_SINCE_TRIFFID = ASTEPS_SINCE_TRIFFID + 1
      ENDIF

! map JULES prognostics to module prognostics before use
        snowdepth = snowdepth_p
      IF (l_flake_model) THEN
        lake_h_ice = lake_h_ice_p
      END IF

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, l, weight2, weight1,   & 
!$OMP& weight3, temp)

! set p at layer boundaries.
! NB: some arrays have haloes but are unset, if never used they will
!     be removed.

!$OMP DO SCHEDULE(STATIC)
      Do j = 1, rows
        Do i = 1, row_length
          p_layer_boundaries(i,j,0) = p_star(i,j)
          p_layer_centres(i,j,0) = p_star(i,j)
          exner_layer_boundaries(i,j,0) = (p_layer_boundaries(i,j,0)/   &
     &                                     p_zero)**kappa
          exner_layer_centres(i,j,0) = exner_layer_boundaries(i,j,0)

        End Do
      End Do
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
      Do k = 1, model_levels - 1
        Do j = 1, rows
          Do i = 1, row_length
            p_layer_boundaries(i,j,k) = p(i,j,k+1)
            p_layer_centres(i,j,k) = p_theta_levels(i,j,k)
            exner_layer_boundaries(i,j,k) = exner_rho_levels(i,j,k+1)
            exner_layer_centres(i,j,k) = exner_theta_levels(i,j,k)
          End Do
        End Do
      End Do
!$OMP END DO NOWAIT

      k = model_levels

!$OMP DO SCHEDULE(STATIC)
      Do j = 1, rows
        Do i = 1, row_length
          p_layer_boundaries(i,j,model_levels) = 0.0
          p_layer_centres(i,j,model_levels) =                           &
     &                           p_theta_levels(i,j,model_levels)
          exner_layer_boundaries(i,j,k) = 0.
          exner_layer_centres(i,j,k) = exner_theta_levels(i,j,k)
        End Do
      End Do
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
      DO k = 1, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            theta_inc(i,j,k) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
      DO k = 1, qdims%k_end
        DO j = qdims%j_start, qdims%j_end
          DO i = qdims%i_start, qdims%i_end
            q_inc(i,j,k) = 0.0
            qcl_inc(i,j,k) = 0.0
            qcf_inc(i,j,k) = 0.0
            cf_inc(i,j,k) = 0.0
            cfl_inc(i,j,k) = 0.0
            cff_inc(i,j,k) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
      DO k = udims_s%k_start, udims_s%k_end
        DO j = udims_s%j_start, udims_s%j_end
          DO i = udims_s%i_start, udims_s%i_end
            u_inc(i,j,k) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
      
!$OMP DO SCHEDULE(STATIC)
      DO k = vdims_s%k_start, vdims_s%k_end
        DO j = vdims_s%j_start, vdims_s%j_end
          DO i = vdims_s%i_start, vdims_s%i_end
            v_inc(i,j,k) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          land_fract(i,j) = 0.0
        END DO
      END DO
!$OMP END DO

! ----------------------------------------------------------------------
! Set Coastal tiling dependent prognostics:
! ----------------------------------------------------------------------
      IF(L_CTILE)THEN

!$OMP DO SCHEDULE(STATIC)
        DO J = 1, ROWS
          DO I = 1, ROW_LENGTH
            TSTAR_LAND(I,J)=TSTAR_LAND_CTILE(I,J)
            TSTAR_SEA(I,J)=TSTAR_SEA_CTILE(I,J)
            LAND_ALB(I,J)=LAND_ALB_CTILE(I,J)
            SICE_ALB(I,J)=SICE_ALB_CTILE(I,J)
            LAND0P5(I,J)=.FALSE.
          ENDDO
        ENDDO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
        DO K = 1, NICE_USE
          DO J = 1, ROWS
            DO I = 1, ROW_LENGTH
              TSTAR_SICE_CAT(I,J,K)=TSTAR_SICE_CTILE(I,J,K)
            ENDDO
          ENDDO
        ENDDO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
        DO L=1,LAND_POINTS
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          FLAND(L)=FLAND_CTILE(L)
          IF(FLAND(L) >= 0.5)LAND0P5(I,J)=.TRUE.
          land_fract(i,j) = fland(l)
        ENDDO
!$OMP END DO NOWAIT

      ELSE
! No coastal tiling so land fraction is just 1 for all land points.
!$OMP DO SCHEDULE(STATIC)
        DO L=1,LAND_POINTS
          J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
          I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
          fland(l) = 1.0
          land_fract(i,j) = 1.0
        ENDDO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC) 
!(Note, nice_use must = 1 if not using coastal tiling)
        DO J = 1, ROWS
          DO I = 1, ROW_LENGTH
            TSTAR_LAND(I,J)=T_SURF(I,J)
            IF(.NOT.LAND_SEA_MASK(I,J))THEN
              IF(ICE_FRACT(I,J) <= 0.0)THEN
                TSTAR_SEA(I,J)=T_SURF(I,J)
                TSTAR_SICE_CAT(I,J,1)=T_SURF(I,J)
              ELSE
                TSTAR_SEA(I,J)=TFS
                TSTAR_SICE_CAT(I,J,1)=(T_SURF(I,J)             &
     &            -(1.-ICE_FRACT(I,J))*TSTAR_SEA(I,J))         &
                      /ICE_FRACT(I,J)
              ENDIF
            ELSE
              TSTAR_SEA(I,J)=T_SURF(I,J)
              TSTAR_SICE_CAT(I,J,1)=T_SURF(I,J)
            ENDIF
            LAND0P5(I,J)=LAND_SEA_MASK(I,J)

            LAND_ALB(I,J)=RMDI
            SICE_ALB(I,J)=RMDI
          ENDDO
        ENDDO
!$OMP END DO  NOWAIT
      ENDIF

      IF (.NOT. l_endgame ) THEN
! ----------------------------------------------------------------------
! Create an unscaled dry density variable
! ----------------------------------------------------------------------
! NOTE: Although this density variable is calculated here, it is not
!       actually used in this subroutine. It is, though, passed to
!       atmos_physics2 where it is used if L_mr_physics2 is true.
!       Therefore, if L_mr_physics2 is true, it is also necessary that
!       L_mr_physics1 also be true, so that this variable is calculated
!       here.
      IF (L_mixing_ratio) THEN

!$OMP DO SCHEDULE(STATIC)
      DO k = 1, wet_levels
        DO j = 1-halo_j, rows+halo_j
          DO i = 1-halo_i, row_length+halo_i
            moist(i,j,k) = q(i,j,k) + qcl(i,j,k) + qcf(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO 

      IF(L_mcr_qcf2)THEN

!$OMP DO SCHEDULE(STATIC)
         DO k = 1, wet_levels
           DO j = 1-halo_j, rows+halo_j
             DO i = 1-halo_i, row_length+halo_i
               moist(i,j,k)      = moist(i,j,k) + qcf2(i,j,k)
             END DO
           END DO
        END DO
!$OMP END DO 

      ENDIF

      IF(L_mcr_qrain)THEN

!$OMP DO SCHEDULE(STATIC)
        DO k = 1, wet_levels
          DO j = 1-halo_j, rows+halo_j
            DO i = 1-halo_i, row_length+halo_i
              moist(i,j,k)      = moist(i,j,k) + qrain(i,j,k)
            END DO
          END DO
        END DO
!$OMP END DO 

      END IF
      
      IF(L_mcr_qgraup)THEN

!$OMP DO SCHEDULE(STATIC)
        DO k = 1, wet_levels
          DO j = 1-halo_j, rows+halo_j
            DO i = 1-halo_i, row_length+halo_i
              moist(i,j,k)      = moist(i,j,k) + qgraup(i,j,k)
            END DO
          END DO
        END DO
!$OMP END DO 
      END IF


      k = 1

!$OMP DO SCHEDULE(STATIC)
        DO j = pdims_s%j_start, pdims_s%j_end
          DO i = pdims_s%i_start, pdims_s%i_end
            unscaled_dry_rho(i,j,k) = rho(i,j,k) /                      &
               ((1.0 + moist(i,j,k)) * r_rho_levels(i,j,k)**2)
          END DO
        END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      Do k = 2, qdims_s%k_end
        Do j = pdims_s%j_start, pdims_s%j_end
          Do i = pdims_s%i_start, pdims_s%i_end
            weight2 = r_rho_levels(i,j,k)                               &
                      - r_theta_levels(i,j,k-1)
            weight1 = r_theta_levels(i,j,k)                             &
                      - r_rho_levels(i,j,k)
            weight3 = r_theta_levels(i,j,k)                             &
                      - r_theta_levels(i,j,k-1)
            temp = (weight2 * moist(i,j,k) + weight1 * moist(i,j,k-1))  &
                   / weight3
            unscaled_dry_rho(i,j,k) = rho(i,j,k) /                      &
               ((1.0 + temp) * r_rho_levels(i,j,k)**2)
          End Do
        End Do
      End Do
!$OMP END DO

      k = qdims_s%k_end + 1
      If ( k <= pdims_s%k_end ) Then

!$OMP DO SCHEDULE(STATIC)
        Do j = pdims_s%j_start, pdims_s%j_end
          Do i = pdims_s%i_start, pdims_s%i_end
            weight2 = r_rho_levels(i,j,k)                               &
                      - r_theta_levels(i,j,k-1)
            weight1 = r_theta_levels(i,j,k)                             &
                      - r_rho_levels(i,j,k)
            weight3 = r_theta_levels(i,j,k)                             &
                      - r_theta_levels(i,j,k-1)
            temp = weight1 *  moist(i,j,k-1) / weight3
            unscaled_dry_rho(i,j,k) = rho(i,j,k) /                      &
               ((1.0 + temp) * r_rho_levels(i,j,k)**2)
          End Do
        End Do
!$OMP END DO

      End If !  k <= model_levels

!$OMP DO SCHEDULE(STATIC)
      Do k = qdims_s%k_end+2, pdims_s%k_end
        Do j = pdims_s%j_start, pdims_s%j_end
          Do i = pdims_s%i_start, pdims_s%i_end
            unscaled_dry_rho(i,j,k) = rho(i,j,k)                        &
                                      /(r_rho_levels(i,j,k)**2)
          End Do
        End Do
      End Do
!$OMP END DO NOWAIT

      END IF
      END IF  ! .NOT. l_endgame

!$OMP END PARALLEL

!-----------------------------------------------------------------------
! Initialise sea and sea-ice variables for JULES that require a land
! mask to be available  (also used by MOSES in the radiation scheme)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Set up index for sea and sea-ice
!-----------------------------------------------------------------------
      ssi_pts = 0
      ssi_index(:)=0
      DO j=1,rows
        DO i=1,row_length
          ssi_pts=ssi_pts + 1
          IF ( land_fract(i,j) < 1.0 ) THEN
            ssi_index(ssi_pts) = (j-1)*row_length + i
          END IF
          fssi(i,j)=1.0 - land_fract(i,j)
        END DO
      END DO

!-----------------------------------------------------------------------
! Allocate space for and initialise sea and sea-ice indices.
!-----------------------------------------------------------------------      

      sea_pts = 0
      sice_pts = 0
      sea_index(:)=0
      sice_index(:)=0
      sice_frac(:)=0.0
      sea_frac(:)=0.0
      DO l=1,ssi_pts
        j=(ssi_index(l)-1)/row_length + 1
        i = ssi_index(l) - (j-1)*row_length
        IF (ssi_index(l) > 0) THEN
          IF (ice_fract(i,j) > 0.0) THEN
            sice_pts=sice_pts+1
            sice_index(sice_pts)=l
            sice_frac(l)=ice_fract(i,j)
          END IF
          IF (ice_fract(i,j) < 1.0) THEN
            sea_pts=sea_pts+1
            sea_index(sea_pts)=l
            sea_frac(l)=1.0 - sice_frac(l)
          END IF
        END IF
      END DO

! NOTE these settings use NICE_USE not NICE so are suitable for use here and
! in the explicit part of the surface exchange code.  They are then updated
! using NICE before the call to the implicit part of the surface exchange 
! code
      sice_pts_ncat(:)=0
      sice_index_ncat(:,:)=0
      sice_frac_ncat(:,:)=0.0
      DO n=1,nice_use
        DO l=1,ssi_pts
          j=(ssi_index(l)-1)/row_length + 1
          i = ssi_index(l) - (j-1)*row_length
          IF (ssi_index(l) > 0) THEN
            IF (ice_fract_cat(i,j,n) > 0.0) THEN
              sice_pts_ncat(n)=sice_pts_ncat(n)+1
              sice_index_ncat(sice_pts_ncat(n),n)=l
              sice_frac_ncat(l,n)=ice_fract_cat(i,j,n)
            END IF
          END IF
        END DO
      END DO

! ----------------------------------------------------------------------
! Section ENG.1  Add energy correction increments to temperature
! ----------------------------------------------------------------------
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Energy Correct.',5)

      If (L_emcorr) then

        If (Lflux_reset) then
! reinitialise net flux field at beginning of energy correction period
            Do j = 1, rows
              Do i = 1, row_length
              sum_eng_fluxes(i,j)=0.0
              sum_moist_flux(i,j)=0.0
            Enddo
          Enddo
        Endif

! Add energy correction increments every timestep.
! This is a temperature increment

        Call add_eng_corr (energy_correction,Theta_inc,                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &                         STASHwork14)

      End If    ! (L_emcorr)

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Energy Correct.',6)

! ----------------------------------------------------------------------
! Section Set-up time-level n
! ----------------------------------------------------------------------

      Do k =             1, tdims%k_end
        Do j = tdims%j_start, tdims%j_end
          Do i = tdims%i_start, tdims%i_end
            T_n(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
          End Do
        End Do
      End Do

      Do k =             1, qdims%k_end
        Do j = qdims%j_start, qdims%j_end
          Do i = qdims%i_start, qdims%i_end
            q_n(i,j,k) = q(i,j,k)
            qcl_n(i,j,k) = qcl(i,j,k)
            qcf_n(i,j,k) = qcf(i,j,k)
            cf_n(i,j,k)  = cf(i,j,k)
            cfl_n(i,j,k) = cfl(i,j,k)
            cff_n(i,j,k) = cff(i,j,k)
          End Do
        End Do
      End Do

      IF (pc2_falliceshear_method == real_shear) THEN 

        CALL u_to_p (u, udims_s%i_start,udims_s%i_end,                  &  
                     udims_s%j_start,udims_s%j_end,                     & 
                     pdims%i_start,pdims%i_end,                         &  
                     pdims%j_start,pdims%j_end,                         & 
                     model_levels,                                      & 
                     model_domain,at_extremity, u_on_p) 

        CALL v_to_p (v, vdims_s%i_start,vdims_s%i_end,                  &  
                     vdims_s%j_start,vdims_s%j_end,                     & 
                     pdims%i_start,pdims%i_end,                         &  
                     pdims%j_start,pdims%j_end,                         & 
                     model_levels,                                      & 
                     model_domain,at_extremity, v_on_p) 

IF (.NOT. l_vatpoles) THEN 
    ! set polar values of u and v to zero. 
        IF (model_domain  ==  mt_global) THEN 
          IF (at_extremity(psouth)) THEN 
            DO k = 1, model_levels 
              DO i = 1, row_length 
                v_on_p(i,pdims%j_start,k) = 0.0 
                u_on_p(i,pdims%j_start,k) = 0.0 
              END DO 
            END DO 
          END IF 
          IF (at_extremity(pnorth)) THEN 
            DO k = 1, model_levels 
              DO i = 1, row_length 
                v_on_p(i,pdims%j_end,k) = 0.0 
                u_on_p(i,pdims%j_end,k) = 0.0 
              END DO 
            END DO 
          END IF 
        END IF  ! model_domain  ==  mt_global 
END IF ! vatpoles

      END IF ! pc2_falliceshear_method 

      If (L_CO2_INTERACTIVE) then
       Do k = 1, co2_dim_lev
         Do j = 1, co2_dim_row
           Do i = 1, co2_dim_len
             co2_3D(i,j,k) = co2(i,j,k)
           End Do
         End Do
       End Do
      Else
       co2_3D(1,1,1) = 0.0
      End If

      If (L_mcr_qcf2) Then  ! Second cloud ice variable in use
        Allocate ( qcf2_n(row_length, rows, wet_levels) )
        qcf2_n(:,:,:) = qcf2(1:row_length, 1:rows, 1:wet_levels)
      Else
        Allocate ( qcf2_n(1,1,1) )
      End If

      If (L_mcr_qrain) Then  ! Prognostic rain in use
        Allocate ( qrain_n(row_length, rows, wet_levels) )
        qrain_n(:,:,:) = qrain(1:row_length, 1:rows, 1:wet_levels)
      Else
        Allocate ( qrain_n(1,1,1) )
      End If

      If (L_mcr_qgraup) Then  ! Prognostic graupel in use
        Allocate ( qgraup_n(row_length, rows, wet_levels) )
        qgraup_n(:,:,:) = qgraup(1:row_length, 1:rows, 1:wet_levels)
      Else
        Allocate ( qgraup_n(1,1,1) )
      End If


      If ( model_domain == mt_LAM .and. L_lbc_old ) Then

! In the LAM set qcl_n, qcf_n and cloud fractions to zero on model
! boundaries. This avoids failures due to inconsistences in these fields.
! As the increments on the boundary are purely from the boundary
! conditions we are free to do anything sensible to these values.
! Note that this only applies to the variables in the physics.

        L_zero_boundaries=.TRUE.
        L_zero_halos=.FALSE.

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,QCL_N,               &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,QCF_N,               &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

        If (L_mcr_qcf2)                                                 &
                        ! prognostic second cloud ice in use
! DEPENDS ON: zero_lateral_boundaries
     &   CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &    row_length,rows,0,0,wet_levels,fld_type_p,qcf2_n,             &
     &    1, AT_EXTREMITY,                                              &
     &    L_zero_boundaries,L_zero_halos)

        If (L_mcr_qrain)                                                &
                         ! prognostic rain in use
! DEPENDS ON: zero_lateral_boundaries
     &   CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &    row_length,rows,0,0,wet_levels,fld_type_p,qrain_n,            &
     &    1, AT_EXTREMITY,                                              &
     &    L_zero_boundaries,L_zero_halos)

        If (L_mcr_qgraup)                                               &
                          ! prognostic graupel in use
! DEPENDS ON: zero_lateral_boundaries
     &   CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &    row_length,rows,0,0,wet_levels,fld_type_p,qgraup_n,           &
     &    1, AT_EXTREMITY,                                              &
     &    L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,                     &
     &   AREA_CLOUD_FRACTION,                                           &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,                     &
     &   cf_n,                                                          &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,                     &
     &   cfl_n,                                                         &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,                     &
     &   cff_n,                                                         &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)


      End If !  model_domain == mt_LAM .and. L_lbc_old


! ----------------------------------------------------------------------
! Section Communications.
! ----------------------------------------------------------------------


!     Find the tropopause level
!     -------------------------
      IF ( (l_rad_step_diag.OR.l_rad_step_prog).AND.                    &
!       If climatological aerosols are present or any diagnostics
!       requested for fluxes at the tropopause then it is necessary
!       to find the tropopause.
        ( l_climat_aerosol.OR.                                          &
          sf_calc(237,1).OR.sf_calc(238,1).OR.                          &
          sf_calc(237,2).OR.sf_calc(238,2) ) ) THEN

!       The variables min_trop_level and max_trop_level
!       index rho-levels in the control routines,
!       which start from a zeroth level at the surface.
!       Layer boundaries in physical routines are indexed with
!       the convention that the lowest is indexed by 1.
!       Since the first rho-level is omitted from the physics
!       grid, the two methods of indexing refer to the same
!       horizontal level from the second rho-level upwards,
!       so there is actually no need to adjust these variables
!       for the change of indexing convention.

!       Set SCM dummy values to zero
        scm_dummy_1d(:,:)   = 0.0
        scm_dummy_2d(:,:,:) = 0.0

! DEPENDS ON: tropin
        CALL tropin (t_n, exner_rho_levels,                             &
          exner_theta_levels(:,:,1:tdims%k_end),                        &
          row_length, rows, model_levels, offx, offy,                   &
          at_extremity,scm_dummy_1d,scm_dummy_2d,                       &
          min_trop_level, max_trop_level, trindx )

      END IF


!     Generate sub-grid cloud field
!     -----------------------------
! DEPENDS ON: open_cloud_gen
      CALL open_cloud_gen (                                             &

!       Parallel variables
        global_row_length, global_rows,                                 &
        mype,n_proc, at_extremity,                                      &

!       Model dimensions.
        row_length, rows, model_levels, wet_levels,                     &
        row_length*rows, p_layer_centres,                               &

!       Properties of clouds
        area_cloud_fraction, dp_corr_strat, cct, global_cloud_top,      &

!       Model switches
        L_Rad_Step_diag, L_Rad_Step_prog, model_domain,                 &

!       Time stepping information.
        val_year, val_day_number, val_hour, val_minute, val_second,     &

!       Error information
        Error_code  )


!     Diagnostic RHcrit
!     -----------------
      IF (l_rhcpt) THEN

!       Dimension diagnostic 3D RHcrit array
        rhc_row_length = qdims%i_end - qdims%i_start + 1
        rhc_rows       = qdims%j_end - qdims%j_start + 1

        ALLOCATE(ext_p_layer_centres(0:rhc_row_length+1,0:rhc_rows+1,     &
                                               0:qdims%k_end))
        ALLOCATE(ext_tl(0:rhc_row_length+1, 0:rhc_rows+1, 1: qdims%k_end))
        ALLOCATE(ext_ql(0:rhc_row_length+1, 0:rhc_rows+1, 1: qdims%k_end))
        ALLOCATE(ext_qcf(0:rhc_row_length+1,0:rhc_rows+1, 1: qdims%k_end))
        ALLOCATE(ext_ice_frac(0:rhc_row_length+1,0:rhc_rows+1))
        ALLOCATE(ext_land_frac(0:rhc_row_length+1,0:rhc_rows+1))

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
        DO k = 1, qdims%k_end
          DO j = 1, rhc_rows
            DO i = 1, rhc_row_length

              ext_p_layer_centres(i,j,k) = p_layer_centres(i,j,k)

              IF (l_pc2) THEN
                ext_tl(i,j,k) = t_n(i,j,k)                               &
                                - (lc * qcl_n(i,j,k) ) / cp
                ext_ql(i,j,k) = q_n(i,j,k) + qcl_n(i,j,k)
              ELSE  ! l_pc2
                ext_tl(i,j,k) = t_n(i,j,k)
                ext_ql(i,j,k) = q_n(i,j,k)
              END IF  ! l_pc2
              ext_qcf(i,j,k) = qcf_n(i,j,k)

            END DO ! i
          END DO ! j
        END DO ! k
!$OMP END DO

!       Rhc_rows_do2:
!$OMP DO SCHEDULE(STATIC) 
        DO j = 1, rhc_rows

!         Rhc_rowlen_do2:
          DO i = 1, rhc_row_length

            ext_p_layer_centres(i,j,0) = p_layer_centres(i,j,0)
            ext_ice_frac(i,j) = ice_fract(i,j)
            ext_land_frac(i,j) = 0.0

          END DO ! Rhc_rowlen_do2

        END DO ! Rhc_rows_do2
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
        DO k = 1, land_points

          j = (land_index(k)-1)/( tdims%i_end - tdims%i_start + 1) + 1
          i = land_index(k) - (j-1)*( tdims%i_end - tdims%i_start + 1)
          ext_land_frac(i,j) = fland(k)

        END DO
!$OMP END DO

!$OMP END PARALLEL 
!       Synchronize haloes.

        i_field = 0

        i_field = i_field + 1
        fields_to_swap(i_field) % field    => ext_p_layer_centres(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  qdims%k_end+1
        fields_to_swap(i_field) % rows        =  rhc_rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => ext_tl(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  qdims%k_end
        fields_to_swap(i_field) % rows        =  rhc_rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => ext_ql(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  qdims%k_end
        fields_to_swap(i_field) % rows        =  rhc_rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => ext_qcf(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  qdims%k_end
        fields_to_swap(i_field) % rows        =  rhc_rows
        fields_to_swap(i_field) % vector      =  .FALSE.

! DEPENDS ON: swap_bounds_mv
        CALL swap_bounds_mv( fields_to_swap, i_field,                   &
                         rhc_row_length, 1, 1)

        i_field = 0
        i_field = i_field + 1
        fields_to_swap(i_field) % field_2d    => ext_land_frac(:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  1
        fields_to_swap(i_field) % rows        =  rhc_rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field_2d    => ext_ice_frac(:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  1
        fields_to_swap(i_field) % rows        =  rhc_rows
        fields_to_swap(i_field) % vector      =  .FALSE.

! DEPENDS ON: swap_bounds_2d_mv
        CALL swap_bounds_2d_mv( fields_to_swap, i_field,                &
                         rhc_row_length, 1, 1)

      ELSE ! l_rhcpt

!       RHcrit will be a 1D parametrized array input from user interface
        rhc_row_length = 1
        rhc_rows = 1

        ALLOCATE(ext_p_layer_centres(1,1,1))
        ALLOCATE(ext_tl(1,1,1))
        ALLOCATE(ext_ql(1,1,1))
        ALLOCATE(ext_qcf(1,1,1))
        ALLOCATE(ext_ice_frac(1,1))
        ALLOCATE(ext_land_frac(1,1))

      END IF ! l_rhcpt

!     Set sea salt arrays
!     -------------------
sea_salt: IF (((l_use_seasalt_direct .OR. l_use_seasalt_indirect).AND.  &
           (l_rad_step_diag.OR.l_rad_step_prog)) .OR.                   &
            L_use_seasalt_autoconv) THEN

       IF ( l_endgame ) THEN
! calculate level 1 winds on p points using subroutine
       CALL uv_p_pnts(u(0:row_length-1,1:rows,1),v(1:row_length,0:n_rows-1,1), &
         cos_theta_longitude,sin_theta_longitude,                              &
         model_domain,global_row_length,gc_proc_row_group,u_1,v_1)

! calculate level 1 windspeed (NB on rho levels)
        DO j = 1, rows
          DO i = 1, row_length
            windspeed_1(i,j)=SQRT(u_1(i,j)**2+v_1(i,j)**2)
          END DO
        END DO

! Q - should the height calculations all be in Section INI. Initialisation of variables.
! Height of rho level 1 above the surface 
        DO j = 1, rows
          DO i = 1, row_length
              height_rho_1(i,j) = r_rho_levels(i,j,1) - r_theta_levels(i,j,0)
          END DO
        END DO

! Height of theta levels above the surface (used in set_seasalt_4A)
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED(height_theta, r_theta_levels, salt_dim1, salt_dim2,       &
!$OMP& salt_dim3) PRIVATE(i, j, k)
        DO k = 1, salt_dim3
          DO j = 1, salt_dim2
            DO i = 1, salt_dim1

              height_theta(i,j,k) = r_theta_levels(i,j,k)               &
                            - r_theta_levels(i,j,0)

            END DO
          END DO
        END DO
!$OMP END PARALLEL DO

! Adjust level 1 windspeed to be 10m windspeed 
        DO j = 1, rows
          DO i = 1, row_length
            windspeed_10m(i,j)=windspeed_1(i,j)*p1/(ALOG(height_theta(i,j,1)/z0_sea)) 
          END DO
        END DO

        CALL set_seasalt_4A(windspeed_10m, height_theta, land_fract,    &
                        ice_fract,row_length, rows, model_levels,       &
                        salt_dim1, salt_dim2, salt_dim3,                &
                        bl_levels, sea_salt_film, sea_salt_jet)


        ELSE

         DO j = udims%j_start, udims%j_end
           DO i = udims%i_start, udims%i_end

             u_1(i-udims%i_start+1, j) = u(i,j,1)

           END DO
         END DO

         DO j = vdims%j_start, vdims%j_end
           DO i = vdims%i_start, vdims%i_end

             v_1(i, j-vdims%j_start+1) = v(i,j,1)

           END DO
         END DO
       
!       Tidy up at North Pole; not done at South Pole as it's land.
!       Polar values are mean of 1st interior row:

         IF (at_extremity(PNorth)) THEN

!          Start point of first interior (i.e. non-polar) row:
           i_start  = rows - 1

!          Sum over points on PEs in order along first interior row:
           CALL global_2d_sums(u_1(:,i_start:i_start), row_length,      &
                              1, 0, 0, 1, u_1_arr, gc_proc_row_group)

           CALL global_2d_sums(v_1(:,i_start:i_start), row_length,      &
                              1, 0, 0, 1, v_1_arr, gc_proc_row_group)

           u_1_mean = u_1_arr(1) / global_row_length
           v_1_mean = v_1_arr(1) / global_row_length

          DO i = 1, row_length

            u_1(i, rows) = u_1_mean
            v_1(i, rows) = v_1_mean

          END DO

        END IF ! at_extremity(PNorth)

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED(height, r_theta_levels, salt_dim1, salt_dim2,             &
!$OMP& salt_dim3) PRIVATE(i, j, k)
        DO k = 1, salt_dim3
          DO j = 1, salt_dim2
            DO i = 1, salt_dim1

              height(i,j,k) = r_theta_levels(i,j,k)                     &
                            - r_theta_levels(i,j,0)

            END DO
          END DO
        END DO
!$OMP END PARALLEL DO

! DEPENDS ON: set_seasalt
        CALL set_seasalt(u_1, v_1, height, land_fract, ice_fract,       &
                        row_length, rows, model_levels,                 &
                        salt_dim1, salt_dim2, salt_dim3,                &
                        bl_levels, sea_salt_film, sea_salt_jet)

      END IF

      ELSE ! l_use_seasalt_direct etc

        sea_salt_film(1, 1, 1) = 0.0
        sea_salt_jet( 1, 1, 1) = 0.0

      END IF sea_salt ! l_use_seasalt_direct etc


!  ^^ Keep this last in the communications section
!     as call to set_seasalt has a load imbalance.

!-----------------------------------------------------------------------
! COSP initialization before calling mycrophysics and radiation
!-----------------------------------------------------------------------
      L_cosp_call = .FALSE.
      IF (L_cosp) CALL cosp_init(L_cosp,L_Rad_Step_prog,L_radiation,       &
         row_length,rows,model_levels,cosp_crain_3d,cosp_csnow_3d,p,       &
         q_n,T_n,t_surf,p_star,p_theta_levels,r_theta_levels,r_rho_levels, &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
         L_cosp_call,cosp_npoints,cosp_cfg,cosp_gbx)

! ----------------------------------------------------------------------
! Section Microphysics. Call microphys_ctl routine
! ----------------------------------------------------------------------

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Microphys (AP1M)',5)
!
      IF(L_rain)then
!
! Only allow dynamic allocation of full space for arrays for 3D
! precipitation diagnostics if they are being used. Otherwise save space
! and give them a minimum size of 1 by 1.
      IF ( SF(222,4) .OR. SF(223,4) .OR. SF(224,4) .OR. SF(225,4)       &
     &               .OR. L_DUST                                        &
     &               .OR. SF(227,4) .OR. L_SULPC_SO2                    &
     &               .OR. L_SOOT .OR. L_BIOMASS .OR. L_OCFF             &
     &               .OR. L_cosp_call) THEN
        LSPICE_DIM1 = row_length
        LSPICE_DIM2 = rows
        LSPICE_DIM3 = wet_levels
      ELSE
        LSPICE_DIM1 = 1
        LSPICE_DIM2 = 1
        LSPICE_DIM3 = 1
      END IF

! DEPENDS ON: microphys_ctl
      Call microphys_ctl (                                              &

! Parallel variables
     &  halo_i, halo_j, offx, offy, global_row_length                   &
     &, at_extremity                                                    &

! model dimensions.
     &, row_length, rows, rhc_row_length, rhc_rows, n_rows, land_points &
     &, model_levels, wet_levels, bl_levels                             &
     &, lspice_dim1,lspice_dim2,lspice_dim3                             &
     &, salt_dim1, salt_dim2, salt_dim3                                 &
      , cdnc_dim1, cdnc_dim2, cdnc_dim3                                 &

! Model switches
     &, Ltimer, L_rhcpt, L_dust                                         &
     &, L_sulpc_so2, L_sulpc_nh3, L_soot, L_biomass, L_ocff, L_nitrate  &
     &, L_cosp_call                                                     &

! Model parameters
     &, rhcrit                                                          &

! Primary fields passed in
     &, T_n, q_n, qcl_n, qcf_n, qcf2_n, qrain_n, qgraup_n               &
     &, cf_n, cfl_n, cff_n                                              &
     &, u_on_p, v_on_p                                                  &
     &, snow_depth                                                      &
     &, land_sea_mask, ice_fract                                        &
     &, p_layer_centres, p_layer_boundaries                             &
     &, rho                                                             &
     &, aerosol                                                         &
      , ukca_cdnc                                                       &
     &, dust_div1, dust_div2, dust_div3,                                &
        dust_div4, dust_div5, dust_div6                                 &
     &, so2, nh3, so4_aitken, so4_accu, so4_diss                        &
     &, soot_agd, soot_cld, bmass_agd, bmass_cld                        &
     &, ocff_agd, ocff_cld, nitr_acc, nitr_diss, biogenic               &
     &, sea_salt_film, sea_salt_jet, arcl                               &

! Other fields passed in
     &, ntml, cumulus                                                   &
     &, fland, land_index                                               &
! Variables for stochastic physics random parameters
     &, m_ci                                                            &
 

! diagnostic info
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &  stashwork4                                                      &
!
! SCM diagnostics switches (dummy in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &

! Increment fields passed in/out
     &, Theta_inc, q_inc, qcl_inc, qcf_inc                              &
     &, qcf2_inc, qrain_inc, qgraup_inc                                 &
     &, cf_inc, cfl_inc, cff_inc                                        &


! Fields required elsewhere
      , ls_rain, ls_snow, micro_tends, cosp_gbx, n_drop_pot             &
! Field for Rh crit parametrization
      , ext_p_layer_centres                                             &
      , ext_tl                                                          &
      , ext_ql                                                          &
      , ext_qcf                                                         &
      , ext_ice_frac                                                    &
      , ext_land_frac                                                   &

! error information
     &, Error_code  )

      else
        ls_rain(:,:) = 0.0
        ls_snow(:,:) = 0.0
        micro_tends(:,:,:,:) = 0.0
      endif

! Deallocate arrays used for RH crit parametrization
      DEALLOCATE(ext_land_frac)
      DEALLOCATE(ext_ice_frac)
      DEALLOCATE(ext_qcf)
      DEALLOCATE(ext_ql)
      DEALLOCATE(ext_tl)
      DEALLOCATE(ext_p_layer_centres)

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Microphys (AP1M)',6)

! ----------------------------------------------------------------------
! Section Turbulence. Call pc2_turbulence_ctl routine
! ----------------------------------------------------------------------
! Earlier versions of PC2 had the erosion term included here,
! parallel with the rest of the slow physics. We have since decided to
! move the erosion to be parallel with the PC2 response to convection
! since this results in improved numerical balances in shallow
! convection. 
! In the case of PC2 using diagnostic shallow cloud it is better
! to do the erosion term here so the call is now under a switch for this.

      IF (l_pc2 .AND. l_pc2_diag_sh .AND. .NOT. l_micro_eros ) THEN

! DEPENDS ON: pc2_turbulence_ctl
        Call pc2_turbulence_ctl (                                       &

! Primary fields passed in, unchanged on exit
     &    T_n, q_n, qcl_n, cf_n, cfl_n, cff_n                           &
     &,   p_layer_centres(1,1,1),                                       &

! diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &    STASHwork4                                                    &
!
! SCM diagnostics switches (dummy in full UM)
     &,   nSCMDpkgs, L_SCMDiags                                         &

! Increment fields passed in/out, updated on exit
     &,   Theta_inc, q_inc, qcl_inc, cf_inc, cfl_inc)

      END IF  ! L_pc2 and l_pc2_diag_sh
!
! ----------------------------------------------------------------------
! Section RAD   Radiation scheme.
!               This incorporates radiation code for non-radiation
!               timesteps which is in CLDCTL1.dk in the UM.
!-----------------------------------------------------------------------

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Radiation (AP1R)',5)
!  Set dimensions of mineral dust arrays for passing to rad_ctl
      IF (L_DUST) THEN
        DUST_DIM1 = ROW_LENGTH*ROWS
        DUST_DIM2 = MODEL_LEVELS
      ELSE
        DUST_DIM1 = 1
        DUST_DIM2 = 1
      END IF
!  Set dimensions of _SULPHATE arrays for passing to RAD_CTL
      IF (L_SULPC_SO2) THEN
        SULP_DIM1 = rows*row_length
        SULP_DIM2 = model_levels
      ELSE
        SULP_DIM1 = 1
        SULP_DIM2 = 1
      END IF
!  Set dimensions of soot arrays for passing to RAD_CTL
      If (l_soot) then
        soot_dim1 = rows*row_length
        soot_dim2 = model_levels
      Else
        soot_dim1 = 1
        soot_dim2 = 1
      End If
!  Set dimensions of array for aerosol climatology for NWP
      If (n_arcl_species > 0) then
        arcl_dim1 = rows*row_length
        arcl_dim2 = model_levels
      Else
        arcl_dim1 = 1
        arcl_dim2 = 1
      End If
!  Set dimensions of biomass arrays for passing to RAD_CTL
      If (l_biomass) then
        bmass_dim1 = rows*row_length
        bmass_dim2 = model_levels
      Else
        bmass_dim1 = 1
        bmass_dim2 = 1
      End If
!  Set dimensions of biogenic array for passing to RAD_CTL
      IF (L_USE_BIOGENIC) THEN
        biogenic_dim1 = rows*row_length
        biogenic_dim2 = model_levels
      ELSE
        biogenic_dim1 = 1
        biogenic_dim2 = 1
      ENDIF
!  Set dimensions of OCFF array for passing to RAD_CTL
      If (l_ocff) then
        ocff_dim1 = rows*row_length
        ocff_dim2 = model_levels
      Else
        ocff_dim1 = 1
        ocff_dim2 = 1
      End If
!  Set dimensions of nitrate aerosol array for passing to RAD_CTL
      If (L_nitrate) then
        nitrate_dim1 = rows*row_length
        nitrate_dim2 = model_levels
      Else
        nitrate_dim1 = 1
        nitrate_dim2 = 1
      End If
!  Set dimensions of UKCA_RADAER arrays for passing to RAD_CTL
      IF (L_ukca_radaer) THEN
        ukca_dim1 = rows*row_length
        ukca_dim2 = model_levels
      ELSE
        ukca_dim1 = 1
        ukca_dim2 = 1
      END IF

!     Code to calculate mixing ratio of well-mixed greenhouse gases
!       if scenarios for their time variation have been prescribed.

        If ( clim_fcg_nyears(s_co2) > 0 ) Then  !  CO2 level calculated
! depends on: gas_calc
          Call gas_calc ( co2_mmr,                                      &
            clim_fcg_nyears(s_co2),  clim_fcg_years(1,s_co2),           &
            clim_fcg_levls(1,s_co2), clim_fcg_rates(1,s_co2),           &
            lenscen, error_code)
          If ( error_code /= 0 ) Then 
            IF (lhook) CALL dr_hook('ATMOS_PHYSICS1',zhook_out,zhook_handle)
            RETURN
          End if
        End if

        If ( clim_fcg_nyears(s_n2o) > 0 ) Then  !  same for N2O
! depends on: gas_calc
          Call gas_calc ( n2ommr,                                       &
             clim_fcg_nyears(s_n2o),  clim_fcg_years(1,s_n2o),          &
             clim_fcg_levls(1,s_n2o), clim_fcg_rates(1,s_n2o),          &
             lenscen, error_code)
          If ( error_code /= 0 ) Then 
            IF (lhook) CALL dr_hook('ATMOS_PHYSICS1',zhook_out,zhook_handle)
            RETURN
          End if
        End if

        If ( clim_fcg_nyears(s_ch4) > 0 ) Then  !  CH4
! depends on: gas_calc
          Call gas_calc ( ch4mmr,                                       &
             clim_fcg_nyears(s_ch4),  clim_fcg_years(1,s_ch4),          &
             clim_fcg_levls(1,s_ch4), clim_fcg_rates(1,s_ch4),          &
             lenscen, error_code)
          If ( error_code /= 0 ) Then 
            IF (lhook) CALL dr_hook('ATMOS_PHYSICS1',zhook_out,zhook_handle)
            RETURN
          End if
        End if

        If ( clim_fcg_nyears(s_cfc11) > 0 ) Then   !  "CFC11"
! depends on: gas_calc
          Call gas_calc ( c11mmr,                                       &
             clim_fcg_nyears(s_cfc11),  clim_fcg_years(1,s_cfc11),      &
             clim_fcg_levls(1,s_cfc11), clim_fcg_rates(1,s_cfc11),      &
             lenscen, error_code)
          If ( error_code /= 0 ) Then 
            IF (lhook) CALL dr_hook('ATMOS_PHYSICS1',zhook_out,zhook_handle)
            RETURN
          End if
        End if

        If ( clim_fcg_nyears(s_cfc12) > 0 ) Then   !  "CFC12"
! depends on: gas_calc
          Call gas_calc ( c12mmr,                                       &
             clim_fcg_nyears(s_cfc12),  clim_fcg_years(1,s_cfc12),      &
             clim_fcg_levls(1,s_cfc12), clim_fcg_rates(1,s_cfc12),      &
             lenscen, error_code)
          If ( error_code /= 0 ) Then 
            IF (lhook) CALL dr_hook('ATMOS_PHYSICS1',zhook_out,zhook_handle)
            RETURN
          End if
        End if

        If ( clim_fcg_nyears(s_cfc113) > 0 ) Then  !  "CFC113"
! depends on: gas_calc
          Call gas_calc ( c113mmr,                                      &
             clim_fcg_nyears(s_cfc113),  clim_fcg_years(1,s_cfc113),    &
             clim_fcg_levls(1,s_cfc113), clim_fcg_rates(1,s_cfc113),    &
             lenscen, error_code)
          If ( error_code /= 0 ) Then 
            IF (lhook) CALL dr_hook('ATMOS_PHYSICS1',zhook_out,zhook_handle)
            RETURN
          End if
        End if

        If ( clim_fcg_nyears(s_cfc114) > 0 ) Then  !  "CFC114"
! depends on: gas_calc
          Call gas_calc ( c114mmr,                                      &
             clim_fcg_nyears(s_cfc114),  clim_fcg_years(1,s_cfc114),    &
             clim_fcg_levls(1,s_cfc114), clim_fcg_rates(1,s_cfc114),    &
             lenscen, error_code)
          If ( error_code /= 0 ) Then
            IF (lhook) CALL dr_hook('ATMOS_PHYSICS1',zhook_out,zhook_handle)
            RETURN
          End if
        End if

        If ( clim_fcg_nyears(s_hcfc22) > 0 ) Then  !  "HCFC22"
! depends on: gas_calc
          Call gas_calc ( hcfc22mmr,                                    &
             clim_fcg_nyears(s_hcfc22),  clim_fcg_years(1,s_hcfc22),    &
             clim_fcg_levls(1,s_hcfc22), clim_fcg_rates(1,s_hcfc22),    &
             lenscen, error_code)
          If ( error_code /= 0 ) Then 
            IF (lhook) CALL dr_hook('ATMOS_PHYSICS1',zhook_out,zhook_handle)
            RETURN
          End if
        End if

        If ( clim_fcg_nyears(s_hfc125) > 0 ) Then  !  "HFC125"
! depends on: gas_calc
          Call gas_calc ( hfc125mmr,                                    &
             clim_fcg_nyears(s_hfc125),  clim_fcg_years(1,s_hfc125),    &
             clim_fcg_levls(1,s_hfc125), clim_fcg_rates(1,s_hfc125),    &
             lenscen, error_code)
          If ( error_code /= 0 ) Then 
            IF (lhook) CALL dr_hook('ATMOS_PHYSICS1',zhook_out,zhook_handle)
            RETURN
          End if
        Endif

        If ( clim_fcg_nyears(s_hfc134a) > 0 ) Then  ! "HFC134a"
! depends on: gas_calc
          Call gas_calc ( hfc134ammr,                                   &
             clim_fcg_nyears(s_hfc134a),  clim_fcg_years(1,s_hfc134a),  &
             clim_fcg_levls(1,s_hfc134a), clim_fcg_rates(1,s_hfc134a),  &
             lenscen, error_code)
          If ( error_code /= 0 ) Then
            IF (lhook) CALL dr_hook('ATMOS_PHYSICS1',zhook_out,zhook_handle)
            RETURN
          End if
        End if

! Copy greenhouse gases into greenhouse gas array
      ALLOCATE(grgas_field(row_length, rows, model_levels, ngrgas))

! check to see if UKCA is specifying trace gases, and fill fields
      DO i=1,ngrgas
        IF (grgas_addr(i) > 0) grgas_field(:,:,:,i) =                   &
         tracer_ukca(1:row_length,1:rows,1:model_levels,grgas_addr(i))  

! For some configurations CFC-11 is held as a lumped species, so
!  rescale field by the ratio of CFC-11/field at the surface
        IF (i == p_f11 .AND. grgas_addr(i) > 0 .AND.                    &
           (L_ukca_strat .OR. L_ukca_strattrop)) THEN
          DO k=model_levels,1,-1
            grgas_field(:,:,k,p_f11) = grgas_field(:,:,k,p_f11) *       &
              c11mmr / grgas_field(:,:,1,p_f11)
          END DO
        END IF

! For some configurations CFC-12 is held as a lumped species, so
!  rescale field by the ratio of CFC-12/field at the surface
        IF (i == p_f12 .AND. grgas_addr(i) > 0 .AND.                    &
           (L_ukca_strat .OR. L_ukca_strattrop)) THEN
          DO k=model_levels,1,-1
            grgas_field(:,:,k,p_f12) = grgas_field(:,:,k,p_f12) *       &
              c12mmr / grgas_field(:,:,1,p_f12)
          END DO
        END IF

        ! Make sure the field if used is all positive.
        IF (grgas_addr(i) > 0) THEN
          grgas_field(:,:,:,i) = MAX(grgas_field(:,:,:,i),0.0)
        ELSE
          ! Set to RMDI if not used.
          grgas_field(:,:,:,i) = rmdi
        END IF
      END DO    ! Loop over grgas species

      IF (grgas_addr(p_o3) <= 0)                                        &
        grgas_field(:,:,1:ozone_levels,p_o3) = ozone

      IF (grgas_addr(p_ch4) <= 0)                                       &
        grgas_field(:,:,:,p_ch4) = ch4mmr          ! global constant

      IF (grgas_addr(p_n2o) <= 0)                                       &
        grgas_field(:,:,:,p_n2o) = n2ommr          ! global constant

      IF (grgas_addr(p_f11) <= 0)                                       &
        grgas_field(:,:,:,p_f11) = c11mmr          ! global constant

      IF (grgas_addr(p_f12) <= 0)                                       &
        grgas_field(:,:,:,p_f12) = c12mmr          ! global constant

      IF (grgas_addr(p_f113) <= 0)                                      &
        grgas_field(:,:,:,p_f113) = c113mmr        ! global constant

      IF (grgas_addr(p_f22) <= 0)                                       &
        grgas_field(:,:,:,p_f22) = hcfc22mmr       ! global constant

      grgas_field(:,:,1:wet_levels,p_h2os) =                            & 
        q_n(1:row_length, 1:rows, :)

!  Write trace gas mixing ratios to module for use in UKCA. The value can be 
!  either that which was passed to this routine or the value calculated from 
!  the time interpolation, depending on whether L_CLMCHFCG is true or false.  
      IF (L_ukca .AND.  ((L_ukca_useumuivals .OR.                       &
          L_ukca_set_trace_gases) .OR. L_ukca_prescribech4)) THEN
        IF (L_ukca_prescribech4 .AND. (ch4mmr == rmdi)) THEN
          cmessage = 'Missing value for ch4_mix_ratio'
          Error_code = 1
          WRITE(6,'(A)')'Set value in UMUI panel: atmos_Science_Section_LW_Meth'
          CALL EREPORT('Atmos_Physics1',Error_code,cmessage)
        END IF
        CALL ukca_set_trace_gas_mixratio(                               &
              ch4mmr, co2_mmr, n2ommr, o2mmr,                           &
              c11mmr, c12mmr, c113mmr, c114mmr,                         &
              hcfc22mmr, hfc125mmr, hfc134ammr)
      END IF

!     Set number of cloud droplets to be used for the 1st indirect effect
      IF (l_consistent_cdnc) THEN
!       Use the n_drop_pot array from microphysics in radiation
        l_use_ndrop = .TRUE.
      ELSE IF (l_ukca_aie1) THEN
!       Use the UKCA cdnc array in radiation
        n_drop_pot = ukca_cdnc
        l_use_ndrop = .TRUE.
      END IF
 
      If (L_radiation) then

! DEPENDS ON: ni_rad_ctl
       Call NI_rad_ctl (                                                &

! Parallel variables
     &  halo_i, halo_j, offx, offy, global_row_length, global_rows      &
     &, gc_proc_row_group, gc_proc_col_group, at_extremity              &
     &, n_proc, n_procx, n_procy                                        &
     &, neighbour, g_rows, g_row_length, mype                           &
     &, global_cloud_top                                                &

! model dimensions.
     &, row_length, rows, n_rows                                        &
     &, model_levels, wet_levels, bl_levels                             &
     &, Ozone_levels, cloud_levels, N_cca_levels                        &
     &, NTILES, LAND_POINTS, nice_use, DUST_DIM1, DUST_DIM2             &
     &, biogenic_dim1                                                   &
     &, biogenic_dim2, sulp_dim1, sulp_dim2, soot_dim1, soot_dim2       &
     &, bmass_dim1, bmass_dim2, salt_dim1, salt_dim2, salt_dim3         &
     &, co2_dim_len, co2_dim_row, co2_dim_lev, arcl_dim1, arcl_dim2     &
     &, n_arcl_species, n_arcl_compnts, i_arcl_compnts                  &
     &, ocff_dim1, ocff_dim2, nitrate_dim1, nitrate_dim2                &
     &, ukca_dim1, ukca_dim2, ukca_radaer%n_mode, ukca_radaer%n_cpnt    &

! Model switches
     &, model_domain                                                    &
     &, L_Rad_Step,L_Rad_Step_diag, L_Rad_Step_prog                     &
     &, L_CAL360                                                        &
     &, L_emcorr                                                        &
     &, Ltimer, L_ssice_albedo, L_snow_albedo                           &
     &, L_sice_meltponds,L_sice_scattering,L_sice_hadgem1a,l_cice_alb   &
     &, l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                           &
     &, L_pc2, l_mixing_ratio                                           &
     &, L_MURK_RAD                                                      &
     &, L_co2_interactive                                               &
     &, L_ukca                                                          &
     &, L_USE_ARCL                                                      &
     &, L_DUST, l_sulpc_so2, l_soot, l_biomass, l_ocff, l_nitrate       &
      , L_cosp_call                                                     &
     
! model Parameters
     &, min_trop_level, max_trop_level                                  &
     &, Ntot_land, Ntot_sea                                             &

! in coordinate information
     &, rho                                                             &
     &, delta_lambda, delta_phi                                         &
     &, lat_rot_NP, long_rot_NP                                         &
     &, trindx                                                          &

! in time stepping information.
     &, timestep, radiation_timestep                                    &
     &, radiation_tstep_diag, radiation_tstep_prog                      &
     &, val_year, val_day_number, val_hour, val_minute                  &
     &, val_second, timestep_number                                     &
     &,         PREVIOUS_TIME                                           &

! diagnostic info
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork1,                                                      &
     & STASHwork2                                                       &
!
! SCM diagnostics switches (dummy in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &

! in data fields.
     &, p_star                                                          &
     &, p_layer_boundaries, p_layer_centres                             &
     &, p, p_theta_levels(tdims_s%i_start,tdims_s%j_start,1)            &
     &, exner_rho_levels                                                &
     &, exner_theta_levels(tdims_s%i_start,tdims_s%j_start,1)           &
     &, land_sea_mask, fland,land0p5                                    &
     &, T_SURF,TSTAR_LAND,TSTAR_SEA,TSTAR_SICE_CAT, AREA_CLOUD_FRACTION &
     &, DUST_DIV1(tdims_s%i_start,tdims_s%j_start,1)                    &
     &, DUST_DIV2(tdims_s%i_start,tdims_s%j_start,1)                    &
     &, DUST_DIV3(tdims_s%i_start,tdims_s%j_start,1)                    &
     &, DUST_DIV4(tdims_s%i_start,tdims_s%j_start,1)                    &
     &, DUST_DIV5(tdims_s%i_start,tdims_s%j_start,1)                    &
     &, DUST_DIV6(tdims_s%i_start,tdims_s%j_start,1)                    &
     &, biogenic, so4_aitken(tdims_s%i_start,tdims_s%j_start,1)         &
     &, so4_accu(tdims_s%i_start,tdims_s%j_start,1)                     &
     &, so4_diss(tdims_s%i_start,tdims_s%j_start,1)                     &
     &, soot_new(tdims_s%i_start,tdims_s%j_start,1)                     &
     &, soot_agd(tdims_s%i_start,tdims_s%j_start,1)                     &
     &, soot_cld(tdims_s%i_start,tdims_s%j_start,1)                     &
     &, bmass_new(tdims_s%i_start,tdims_s%j_start,1)                    &
     &, bmass_agd(tdims_s%i_start,tdims_s%j_start,1)                    &
     &, bmass_cld(tdims_s%i_start,tdims_s%j_start,1)                    &
     &, ocff_new(tdims_s%i_start,tdims_s%j_start,1)                     &
     &, ocff_agd(tdims_s%i_start,tdims_s%j_start,1)                     &
     &, ocff_cld(tdims_s%i_start,tdims_s%j_start,1)                     &
     &, nitr_acc(tdims_s%i_start,tdims_s%j_start,1)                     &
     &, nitr_diss(tdims_s%i_start,tdims_s%j_start,1)                    &
     &, aerosol(1:row_length, 1:rows, 1:model_levels), arcl, ukca_radaer&
     &, sea_salt_film, sea_salt_jet, co2_3D                             &
      , cos_zenith_angle, can_rad_mod, frac_control, n_drop_pot         &
! chemical greenhouse gas fields
     &, ngrgas, grgas_field                                             &

! ancillary fields and fields needed to be kept from timestep to
! timestep
     &, snow_depth, snow_depth_sea_cat, ice_fract, ice_fract_cat        &
     &, ice_thick_cat, rgrain, soot                                     &
     &, cca, ccb, cct, cclwp, ccw, lcbase                               &
! chemical ozone; replaces climatological ozone
     &, grgas_field(:,:,1:ozone_levels,p_o3)                            &
     &, SW_incs, LW_incs, dirpar_local                                  &
     &, O3_trop_level, O3_trop_height                                   &
      , T_trop_level, T_trop_height, zh, land_index, albsoil, albobs_sw & 
      , albobs_vis, albobs_nir, lai, snow_tile, tile_frac, tstar_tile   & 
      , z0_tile, dOLR_rts, LW_down, SW_tile_rts                         & 
     &, u_1, v_1                                                        &
     &, land_alb,sice_alb                                               &
     &, ES_SPACE_INTERP, RAD_MASK                                       &

! in/out
! chemical water replaces hydrological water
     &, T_n, grgas_field(:,:,1:wet_levels,p_h2os)                       &
     &, qcl_n, qcf_n, cf_n, cfl_n, cff_n                                &
     &, qcf2_n, qrain_n, qgraup_n                                       &
     &, Theta_inc, q_inc, qcl_inc, cf_inc, cfl_inc                      &
     &, sum_eng_fluxes                                                  &

! out.
      , photosynth_act_rad, rad_hr, dOLR, SW_tile                       &
! COSP arguments
      , cosp_gbx                                                        &
! error information
     &, Error_code  )

! Check error condition
        IF (Error_code > 0) THEN

          CALL ereport(RoutineName, Error_code,                         &
            "Error on return from radiation code.")
        END IF

      else
        photosynth_act_rad = 0.0
        rad_hr = 0.0
        dolr = 0.0
        sw_tile = 0.0
      endif

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Radiation (AP1R)',6)

! DEPENDS ON: close_cloud_gen
      CALL close_cloud_gen

      DEALLOCATE(grgas_field)

!-----------------------------------------------------------------------
!     Call to COSP
!-----------------------------------------------------------------------
      IF (L_cosp_call) THEN
       CALL cosp_main(Ltimer,model_domain,at_extremity,row_length,rows, &
                   n_rows,model_levels,cosp_cfg,cosp_gbx,               &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                   STASHwork2)
        CALL free_cosp_gridbox(cosp_gbx)
      END IF

!-----------------------------------------------------------------------
! Set dirpar_inc dependent prognostic for output to D1 array:
!-----------------------------------------------------------------------

! Never use prognostic dirpar 

!      IF ( H_SECT(34) == "01A" ) THEN
!         L_USE_DIRPAR=.TRUE.
!      ELSE
         L_USE_DIRPAR=.FALSE.
!      ENDIF              
         
      IF (L_USE_DIRPAR) THEN     
         DO J=1, ROWS
           DO I=1, ROW_LENGTH
              DIRPAR_INC(I,J)=DIRPAR_LOCAL(I,J)
           ENDDO
         ENDDO  
      ENDIF      

! ----------------------------------------------------------------------
! Set Coastal tiling dependent prognostics for output to D1 array:
! ----------------------------------------------------------------------
      IF(L_CTILE)THEN
        DO J = 1, ROWS
          DO I = 1, ROW_LENGTH
            LAND_ALB_CTILE(I,J)=LAND_ALB(I,J)
            SICE_ALB_CTILE(I,J)=SICE_ALB(I,J)
          ENDDO
        ENDDO
      ENDIF

      If (Error_code  ==  0) Then

! Convert temperature held in t_n to potential temperature.
! Convert increment after call to gravity wave drag, to allow for
! inclusion of dissipative heating term.

        Do k =             1, tdims%k_end
          Do j = tdims%j_start, tdims%j_end
            Do i = tdims%i_start, tdims%i_end
              T_n(i,j,k) = Theta(i,j,k)
            End Do
          End Do
        End Do

      End If ! on error code equal to zero

! ----------------------------------------------------------------------
! Section CNV.2 Energy correction code
! ----------------------------------------------------------------------

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Conv Eng Corr',5)
      If (L_emcorr .and. Error_code  ==  0) Then

! Add convective + large scale, rain and snow, at the surface to the
! diabatic heating for use in the energy correction
! procedure.
! Scale variables by conversion factor so that only one call is required

        lclf = lc + lf
        Do j = 1, rows
          Do i = 1, row_length
            tot_precip_scaled(i,j) = ls_rain(i,j) * lc +                &
     &                               ls_snow(i,j) * lclf
          End Do
        End Do

! DEPENDS ON: flux_diag
        CALL flux_diag(tot_precip_scaled, xx_cos_theta_latitude,        &
                       row_length, rows ,offx, offy, 1.0,               &
                       sum_eng_fluxes,timestep)

! moist fluxes
        Do j = 1, rows
          Do i = 1, row_length
            tot_precip_scaled(i,j) = -ls_rain(i,j) -ls_snow(i,j)
          End Do
        End Do

! DEPENDS ON: flux_diag
        CALL flux_diag(tot_precip_scaled, xx_cos_theta_latitude,        &
                       row_length, rows ,offx, offy, 1.0,               &
                       sum_moist_flux,timestep)

      End If   ! L_emcorr

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP1 Conv Eng Corr',6)

! ----------------------------------------------------------------------
! Section GWD.1 Call gravity wave drag
! l_gwd      = orographic GWD
! l_use_ussp = middle atmosphere non-orographic GWD scheme
! ----------------------------------------------------------------------

      If ( l_gwd .or. l_use_ussp ) Then

        If (error_code  ==  0 ) Then
! DEPENDS ON: timer
          If (Ltimer) Call timer ('AP1 G-wave drag',5)

! reset p at layer boundaries for spectral GWD code
          Do k = 1, model_levels - 1
            Do j = 1, rows
              Do i = 1, row_length
                p_layer_boundaries(i,j,k) = p_theta_levels(i,j,k)
              End Do
            End Do
          End Do
          Do j = 1, rows
            Do i = 1, row_length
              p_layer_boundaries(i,j,model_levels) = 0.0
            End Do
          End Do

        Call NI_gwd_ctl(                                                &
     &                   halo_i, halo_j, offx, offy                     &
     &,  global_row_length,n_proc, n_procy, gc_proc_row_group           &
     &,                  at_extremity, neighbour                        &

! model dimensions.
     &,                  row_length, rows, n_rows, land_points          &
     &,                  model_levels                                   &

! Model switches
     &,                  model_domain                                   &
! trig arrays
     &,                  sin_theta_longitude, sin_theta_latitude        &

! in coordinate information
     &,                  delta_lambda,delta_phi,true_latitude           &
     &,                  exner_theta_levels                             &

! in time stepping information.
     &,                  timestep, timestep_number                      &

! diagnostic info
     &     ,                                                            &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork6                                                       &
! SCM diagnostics (dummy in full UM)
     &,                  nSCMDpkgs, L_SCMDiags                          &
! in data fields.
     &,                  u, v, land_sea_mask                            &
     &,                  p_layer_boundaries                             &
     &,                  rho, t_n, sd_orog_land                         &
     &,                  orog_grad_xx_land, orog_grad_xy_land           &
     &,                  orog_grad_yy_land, land_index                  &

! in/out
     &,                  u_inc, v_inc,Theta_inc                         &

! error information
     &, Error_code  )

! DEPENDS ON: timer
          If (Ltimer) Call timer ('AP1 G-wave drag',6)
        End If

      End If

      IF (Error_code  ==  0) THEN

! Now convert temperature increment to theta increment
        DO k =             1, tdims%k_end
          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              Theta_inc(i,j,k) = Theta_inc(i,j,k) /                     &
                                 exner_theta_levels(i,j,k)
            END DO
          END DO
        END DO

      END IF! on error code equal to zero


!-----------------------------------------------------------------------
! Tracer Source and Boundary updating where applicable
!-----------------------------------------------------------------------
! Aerosol
      If ( L_murk_source ) Then

        Do level = 1, model_levels
          Call trsrce(                                                  &
     &       rows, row_length, offx, offy                             &
     &,      halo_i, halo_j, model_levels, wet_levels                   &
     &, 0, 0                                                            &
     &,      theta, q_n, qcl_n, qcf_n, exner_rho_levels, rho            &
     &,      aerosol(:,:,level), aerosol_em (:,:,level )                &
     &,      level, timestep, val_hour, val_minute, amp                 &
     &)
        End Do
      End If

      If ( L_murk_bdry ) Then
! DEPENDS ON: trbdry
        Call trbdry(                                                    &
     &       row_length, rows, n_rows, model_levels                     &
     &,      offx, offy, at_extremity                                 &
     &,      p, u, v                                                    &
     &,      aerosol, timestep                                          &
     &)
      End If

      If ( L_murk ) Then
        If (.NOT. L_murk_bdry) Then
          ! Set the external halos if not specified by the
          ! formula for UK mes.

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &         aerosol, row_length, rows,                               &
     &         model_levels, offx, offy, fld_type_p, .false.)

        End If

      End If

!---------------------------------------------------------------------
! Section METHOX Call methane oxidation
!---------------------------------------------------------------------

      If (L_USE_METHOX) Then

        If (error_code == 0) Then
! DEPENDS ON: timer
          If (Ltimer) Call timer ('AP1 NI_methox',5)
!
!    call methane oxidation directly (no interface routine needed)
!
! DEPENDS ON: ni_methox
          call NI_methox(                                               &
! Parallel variables
     &  halo_i, halo_j                                                  &

! model dimensions.
     &, row_length, rows                                                &
     &, model_levels, wet_levels                                        &

! model levels
     &, eta_theta_levels                                                &

! in time stepping information.
     &, timestep                                                        &

! in/out
     &,  q_n,q_inc                                                      &

! error information
     &, Error_code  )

! DEPENDS ON: timer
          If (Ltimer) Call timer ('AP1 NI_methox',6)
        End If
      End If

! ----------------------------------------------------------------------
! Copy increment arrays into star locations
! ----------------------------------------------------------------------

! In the LAM set physics increments on boundary to zero.
      If (model_domain == mt_LAM .and. L_lbc_old) Then

        L_zero_boundaries=.TRUE.
        L_zero_halos=.FALSE.

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,MODEL_LEVELS,fld_type_p,THETA_INC,         &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,Q_INC,               &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,QCL_INC,             &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,QCF_INC,             &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

        If (L_mcr_qcf2)                                                 &
                        ! prognostic second cloud ice in use
! DEPENDS ON: zero_lateral_boundaries
     &   CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &    row_length,rows,0,0,wet_levels,fld_type_p,qcf2_inc,           &
     &    1, AT_EXTREMITY,                                              &
     &    L_zero_boundaries,L_zero_halos)

        If (L_mcr_qrain)                                                &
                         ! prognostic rain in use
! DEPENDS ON: zero_lateral_boundaries
     &   CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &    row_length,rows,0,0,wet_levels,fld_type_p,qrain_inc,          &
     &    1, AT_EXTREMITY,                                              &
     &    L_zero_boundaries,L_zero_halos)

        If (L_mcr_qgraup)                                               &
                          ! prognostic graupel in use
! DEPENDS ON: zero_lateral_boundaries
     &   CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &    row_length,rows,0,0,wet_levels,fld_type_p,qgraup_inc,         &
     &    1, AT_EXTREMITY,                                              &
     &    L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,CF_INC,              &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,CFL_INC,             &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,CFF_INC,             &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,n_ROWS,offx,offy,MODEL_LEVELS,fld_type_v,V_INC,   &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,offx,offy,MODEL_LEVELS,fld_type_u,U_INC,     &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

      End If  !  model_domain == mt_LAM .and. L_lbc_old

      If (L_pc2) then
!
      Allocate( T_inc_diag  (row_length,rows,model_levels) )
      Allocate( q_inc_diag  (row_length,rows,wet_levels) )
      Allocate( qcl_inc_diag(row_length,rows,wet_levels) )
      Allocate( qcf_inc_diag(row_length,rows,wet_levels) )
      Allocate( cf_inc_diag (row_length,rows,wet_levels) )
      Allocate( cfl_inc_diag(row_length,rows,wet_levels) )
      Allocate( cff_inc_diag(row_length,rows,wet_levels) )

!
! Now call a consistency check for the moisture and cloud fields.
! The catch is that _inc variables hold the increments and not the
! full variables yet, so temporarily form them as _inc, check them,
! then rewrite them as increments.
!
! 1. Create full variables (not increments). Only need to do this
!    on wet_levels for theta_inc (theta_inc and T_n hold potential
!    temperature) since the checking is only over wet levels.
!
      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            theta_inc(i,j,k) = (t_n(i,j,k)   + theta_inc(i,j,k))        &
     &                         * exner_theta_levels(i,j,k)
! theta_inc now holds temperature, t_n still holds potential temp
            q_inc(i,j,k)     = q_n(i,j,k)   + q_inc(i,j,k)
            qcl_inc(i,j,k)   = qcl_n(i,j,k) + qcl_inc(i,j,k)
            qcf_inc(i,j,k)   = qcf_n(i,j,k) + qcf_inc(i,j,k)
            cf_inc(i,j,k)    = cf_n(i,j,k)  + cf_inc(i,j,k)
            cfl_inc(i,j,k)   = cfl_n(i,j,k) + cfl_inc(i,j,k)
            cff_inc(i,j,k)   = cff_n(i,j,k) + cff_inc(i,j,k)

            t_inc_diag(i,j,k)   = theta_inc(i,j,k)
            q_inc_diag(i,j,k)   = q_inc(i,j,k)
            qcl_inc_diag(i,j,k) = qcl_inc(i,j,k)
            qcf_inc_diag(i,j,k) = qcf_inc(i,j,k)
            cf_inc_diag(i,j,k)  = cf_inc(i,j,k)
            cfl_inc_diag(i,j,k) = cfl_inc(i,j,k)
            cff_inc_diag(i,j,k) = cff_inc(i,j,k)
          End Do
        End Do
      End Do
!
! 2. Call the checking routine (needs to use temperature, not
!    potential temperature)
!
! DEPENDS ON: pc2_checks
      call pc2_checks(p_layer_centres(1,1,1)                            &
                     ,theta_inc, cf_inc, cfl_inc                        &
                     ,cff_inc, q_inc, qcl_inc, qcf_inc                  &
                     ,l_mixing_ratio)
!
! Now form the increment diagnostics in t_inc_diag etc.
!
      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            t_inc_diag(i,j,k)   = theta_inc(i,j,k) - t_inc_diag(i,j,k)
            q_inc_diag(i,j,k)   = q_inc(i,j,k)     - q_inc_diag(i,j,k)
            qcl_inc_diag(i,j,k) = qcl_inc(i,j,k)   - qcl_inc_diag(i,j,k)
            qcf_inc_diag(i,j,k) = qcf_inc(i,j,k)   - qcf_inc_diag(i,j,k)
            cf_inc_diag(i,j,k)  = cf_inc(i,j,k)    - cf_inc_diag(i,j,k)
            cfl_inc_diag(i,j,k) = cfl_inc(i,j,k)   - cfl_inc_diag(i,j,k)
            cff_inc_diag(i,j,k) = cff_inc(i,j,k)   - cff_inc_diag(i,j,k)
          End Do ! i
        End Do ! j
      End Do ! k
!
! 3. Update fields to produce the net increments which are written to
!    the _star variables. Note that theta_inc
!    currently holds the full temperature value up to wet_levels, and
!    the increment in potential temp from there to model_levels
!
      Do k = wet_levels+1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            theta_star(i,j,k) = theta_inc(i,j,k)
            t_inc_diag(i,j,k) = 0.0 
          End Do
        End Do
      End Do
!
!                                                                           
! 4. Call diagnostic writing subroutine.                                    
! Not called for SCM. 
!                                                                           
! Check that checking diagnostics requested this timestep                   
      If (sf(0,4)) Then                                                     
          Call diagnostics_pc2checks(                                   &
               &                row_length, rows, model_levels          &    
     &,                      wet_levels                                 &    
     &,                      mype,timestep, at_extremity                &    
     &,                      T_inc_diag,q_inc_diag,qcl_inc_diag         &    
     &,                      qcf_inc_diag,cf_inc_diag,cfl_inc_diag      &    
     &,                      cff_inc_diag                               &    
     &,                                                                 &    
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork4                                                       &    
     &       )                                                              
                                                                            
      End If   ! on error_code and sf(0,4)                                  
!                                                                           
! 5. Deallocate fields                                                      
!                                                                           
      Deallocate( T_inc_diag )                                              
      Deallocate( q_inc_diag )                                              
      Deallocate( qcl_inc_diag )                                            
      Deallocate( qcf_inc_diag )                                            
      Deallocate(  cf_inc_diag )                                            
      Deallocate( cfl_inc_diag )                                            
      Deallocate( cff_inc_diag )    
!
! 6. Update fields to produce the net increments which are written to
!    the _star variables.

      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            theta_star(i,j,k) = theta_inc(i,j,k)                        &
     &                         /exner_theta_levels(i,j,k) - t_n(i,j,k)
! theta_star now holds the potential temperature increment
            q_star(i,j,k)     = q_inc(i,j,k)     - q_n(i,j,k)
            qcl_star(i,j,k)   = qcl_inc(i,j,k)   - qcl_n(i,j,k)
            qcf_star(i,j,k)   = qcf_inc(i,j,k)   - qcf_n(i,j,k)
            cf_star(i,j,k)    = cf_inc(i,j,k)    - cf_n(i,j,k)
            cfl_star(i,j,k)   = cfl_inc(i,j,k)   - cfl_n(i,j,k)
            cff_star(i,j,k)   = cff_inc(i,j,k)   - cff_n(i,j,k)
          End Do
        End Do
      End Do

      Else  ! L_pc2

      Do k = 1, model_levels
        Do j = 1, rows
          Do i = 1, row_length
            theta_star(i,j,k) = theta_inc(i,j,k)
          End Do
        End Do
      End Do

      Do k = 1, wet_levels
        Do j = 1, rows
          Do i = 1, row_length
            q_star  (i,j,k) = q_inc  (i,j,k)
            qcl_star(i,j,k) = qcl_inc(i,j,k)
            qcf_star(i,j,k) = qcf_inc(i,j,k)
            cf_star (i,j,k) = cf_inc (i,j,k)
            cfl_star(i,j,k) = cfl_inc(i,j,k)
            cff_star(i,j,k) = cff_inc(i,j,k)
          End Do
        End Do
      End Do

      End If  ! L_pc2

! prognostic second cloud ice in use
      IF (L_mcr_qcf2) THEN
        qcf2_star(1:row_length, 1:rows, 1:wet_levels) =                      &
                                                   qcf2_inc(:,:,1:wet_levels)
      END IF

! prognostic rain in use
      IF (L_mcr_qrain) THEN
        qrain_star(1:row_length, 1:rows, 1:wet_levels) =                     &
                                                   qrain_inc(:,:,1:wet_levels)
      END IF

! prognostic graupel in use
      IF (L_mcr_qgraup)THEN
        qgraup_star(1:row_length, 1:rows, 1:wet_levels) =                    &
                                                  qgraup_inc(:,:,1:wet_levels)
      END IF

      ! Deallocate additional microphysics variables
      Deallocate ( qcf2_inc )
      Deallocate ( qrain_inc )
      Deallocate ( qgraup_inc )
      Deallocate ( qcf2_n )
      Deallocate ( qrain_n )
      Deallocate ( qgraup_n )

! end of routine Atmos_physics1
      IF (lhook) CALL dr_hook('ATMOS_PHYSICS1',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Atmos_Physics1

