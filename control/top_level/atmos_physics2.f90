! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Interface to Atmos Physics parametrizations after S-L advection.
!
! Subroutine Interface:
      SUBROUTINE Atmos_Physics2(                                        &
! Parallel variables
        global_row_length, global_rows, n_proc, n_procx, n_procy        &
      , g_rows, g_row_length, NumCycles, CycleNo                        &

! model dimensions.
      , row_length, rows, n_rows, land_points, model_levels, nice       &
      , nice_use                                                        &
      , wet_levels, bl_levels, dst_levels, dsm_levels, cloud_levels     &
      , land_ice_points, soil_points, n_cca_levels, ntiles, tr_levels   &
      , first_constant_r_rho_level, DIM_CS1, DIM_CS2                    &

! Model switches
      , L_regular, l_mixing_ratio, L_dry, L_lbc_old                     &
      , L_CAL360,  Ltimer                                               &
      , L_area_cloud                                                    &

! Model Parameters
      , rhcrit, CO2_MMR, tr_vars, tr_ukca                               &

! in coordinate information
      , unscaled_dry_rho                                                &
      , delta_lambda, delta_phi                                         &
      , dlambda_p, dphi_p, wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v &
      , lat_rot_NP, long_rot_NP , f3_at_u                               &

! in time stepping information.
      , val_year, val_day_number, val_hour, val_minute                  &
      , val_second,                                                     &

! River routing
       AOCPL_ROW_LENGTH,AOCPL_P_ROWS,XPA,XUA,XVA,YPA,YUA,YVA,           &
       G_P_FIELD, G_R_FIELD, A_STEPS_SINCE_RIV, RIVER_ROW_LENGTH,       &
       RIVER_ROWS, global_river_row_length, global_river_rows,          &
       RIVER_VEL, RIVER_MCOEF, I_RIVER_VN,                              &
! Add inland basin outflow to arguments
       TRIVDIR, TRIVSEQ, TWATSTOR,INLANDOUT_ATM,                        &

!  Add lake evaporation: 
        ACC_LAKE_EVAP,                                                  & 

! Grid-to-grid river routing
        R_AREA, SLOPE, FLOWOBS1, R_INEXT, R_JNEXT, R_LAND,              &
        SUBSTORE, SURFSTORE, FLOWIN, BFLOWIN,                           &
!
! diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
       STASHwork3,STASHwork5,STASHwork8,STASHwork9,STASHwork19,         &
       STASHwork26                                                      &
!
! SCM Diagnostics (dummy values in full UM)
      , nSCMDpkgs, L_SCMDiags                                           &
!
! in data fields.
      , theta, q, qcl, qcf, qrain, qgraup, qcf2                         &
      , rho_rsq, u, v, w, w_adv, p, p_star                              &
     &, exner_rho_levels, exner_theta_levels                            &
     &, land_sea_mask, p_theta_levels                                   &
! ancillary fields and fields needed to be kept from timestep to
! timestep

      , land_index, land_ice_index, soil_index, canopy_gb, snow_depth   &
      , hcon, hcap, smvccl, smvcwt, smvcst, sthf, sthu                  &
      , sil_orog_land, ho2r2_orog, sd_orog, di, ice_fract               &
      , u_0, v_0, u_0_p, v_0_p                                          &
      , cca0, ccb0, cct0, cclwp0, ccw_out, lcbase_out, t_soil           &
      , ti, ti_gb, k_sice_ml                                            &
      , t_surf, z0msea, ice_fract_ncat,di_ncat,satcon,sathh,clapp       &
      , soil_layer_moisture, t1_sd, q1_sd, zh, ddmfx                    &
      , area_cloud_fraction, bulk_cloud_fraction_halos                  &
      , cloud_fraction_liquid_halos, cloud_fraction_frozen_halos        &
      , ls_rain, ls_snow, micro_tends                                   &
      , photosynth_act_rad, rad_hr                                      &
      , SOIL_CLAY,SOIL_SILT,SOIL_SAND,DUST_MREL1,DUST_MREL2,DUST_MREL3  &
      , DUST_MREL4,DUST_MREL5,DUST_MREL6                                &
      , so2_hilem, so2_em, nh3_em, dms_em, soot_hilem, soot_em          &
      , ocff_hilem, ocff_em, co2_emits, co2flux                         &
      , deep_flag, past_precip, past_conv_ht                            & 

! in/out
      , theta_star, q_star, qcl_star, qcf_star, qrain_star, qgraup_star &
      , qcf2_star, cf_star, cfl_star                                    &
      , cff_star, R_u, R_v, R_w, sum_eng_fluxes, sum_moist_flux         &

! In/Out tracer fields
      , aerosol, free_tracers, ukca_tracers                             &
      , DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &

      , so2, dms, so4_aitken, so4_accu, so4_diss, nh3, soot_new         &
      , soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld           &
      , ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss, co2         &

! IN/OUT RIVERS
      , TOT_SURF_RUNOFF, TOT_SUB_RUNOFF                                 &
!
! out fields
      , rhokm, cH_term, ntml, cumulus, nbdsc, ntdsc ,rhcpt              &
      , rhc_row_length, rhc_rows                                        &

! Additional variables for MOSES II
      , frac, frac_disturb, canht_ft, lai_ft, canopy, catch, catch_snow &
      , snow_grnd, snow_tile, z0_tile, z0h_tile_bare                    &
      , t_surf_tile, infil_tile, rgrain                                 &
      , cs, gs, co2_dim_row, co2_dim_len                                &
      , asteps_since_triffid, a_step                                    &
      , g_leaf_acc, g_leaf_phen_acc, npp_ft_acc, resp_w_ft_acc          &
      , resp_s_acc, land_pts_trif, npft_trif, olr, lw_down, sw_tile     &
      , FLAND_CTILE,TSTAR_LAND_CTILE,TSTAR_SEA_CTILE                    &
      , TSTAR_SICE_CAT_CTILE,TSTAR_SICE_CTILE                           &
      , ALBSOIL,COS_ZENITH_ANGLE                                        &

! INOUT variables for TKE based turbulence schemes
      , e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu                    &

! Additional variables required for large-scale hydrology:
      , FEXP,GAMTOT,TI_MEAN,TI_SIG,FSAT,FWETL,ZW                        &
      , STHZW,A_FSAT,C_FSAT,A_FWET,C_FWET                               &
! JULES 2 Prognostics
      , SNOWDEPTH_P, RHO_SNOW_GRND_P                                    &
      , NSNOW_P                                                         &
      , DS_P, SICE_P, SLIQ_P, TSNOWLAYER_P, RHO_SNOW_P, RGRAINL_P       &
! FLake lake scheme prognostics
      , lake_depth_p, lake_fetch_p, lake_t_mean_p, lake_t_mxl_p         &
      , lake_t_ice_p, lake_h_mxl_p, lake_h_ice_p,  lake_shape_p         &
      , lake_g_dt_p                                                     &
! Cariolle ozone and associated ancillaries
      , OZONE_TRACER                                                    &
!
! Additional screen-level variables
      , TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans                          &

! Variables required for COSP (out)
      , cosp_crain_3d, cosp_csnow_3d                                    &

! Variables for SCM and idealised UM
      , L_flux_bc,flux_e,flux_h,L_spec_z0,z0m_scm,z0h_scm,              &
!

! error information
        Error_code  )

 
      USE dynamics_input_mod, ONLY: l_endgame
      USE dynamics_grid_mod,  ONLY: l_vatpoles
      
      USE atmos_physics2_alloc_mod

      USE bl_option_mod, ONLY: l_bl, l_us_blsol, buddy_sea,             &
                               formdrag, ng_stress, explicit_stress,on, &
                               l_quick_ap2

      USE atm_fields_bounds_mod

      USE vertnamelist_mod, ONLY: z_top_of_model

      USE earth_constants_mod, ONLY: g 

      USE atmos_constants_mod, ONLY: kappa, p_zero

      USE water_constants_mod, ONLY: lc, lf, tfs, rho_water

      USE timestep_mod, ONLY: timestep_number, timestep
      
      USE level_heights_mod, ONLY:                                      &
                      r_theta_levels, r_rho_levels,                     &
                      eta_theta_levels, eta_rho_levels
      
      USE rad_input_mod, ONLY: a_sw_radstep_diag, a_lw_radstep_diag,    &
                               a_sw_radstep_prog, a_lw_radstep_prog,    &
                               l_use_cariolle
                                                    

      USE trignometric_mod, ONLY : sec_theta_latitude,                  &
                            sin_theta_longitude, cos_theta_longitude,   &
                            cos_theta_latitude,                         &
                            FV_cos_theta_latitude

      USE cv_run_mod,  ONLY:                                            &
          i_convection_vn,                                              &
          i_convection_vn_0a,                                           &
          i_convection_vn_4a,                                           &
          i_convection_vn_5a,                                           &
          i_convection_vn_6a,                                           &
          cape_opt, cape_bottom, cape_top, l_ccrad, l_3d_cca,           &
          rad_cloud_decay_opt, l_mom, l_rediagnosis, l_pc2_diag_sh,     &
          l_param_conv


      Use cv_param_mod,  Only:                                          &
          rad_decay_off, rad_decay_full_timestep, rad_decay_conv_substep

      USE cv_stash_flg_mod,  ONLY:                                      &
        ! subroutine
         set_convection_output_flags,                                   &
        ! variables
          flg_up_flx, flg_dwn_flx, l_u_incr_conv, l_v_incr_conv,        &
          l_apply_diag 

      Use swapable_field_mod, Only:                                     &
          swapable_field_pointer_type

      USE switches, ONLY:                                               &
          l_aggregate                                                   &
         ,l_flake_model                                                 &
         ,IScrnTDiag, l_snow_albedo, l_sice_multilayers,                &
         l_sice_heatflux, l_top, l_pdm, l_soil_sat_down, can_model,     &
         l_anthrop_heat_src, l_ctile, l_cable

      USE jules_mod, ONLY :  clapp_levs                                 &
                           , sathh_levs                                 &
                           ,  hcap_levs                                 &
                           ,  hcon_levs                                 &
                           ,satcon_levs                                 &
                           ,smvccl_levs                                 &
                           ,smvcwt_levs                                 &
                           ,smvcst_levs
! FLake model
   USE lake_mod, ONLY:  h_snow_min_flk                                  &
                      , u_s_lake                                        &
                      , surf_ht_flux_lake                               &
                      , surf_ht_flux_lk                                 &
                      , sw_down                                         &
                      , lake_depth                                      &
                      , lake_fetch                                      &
                      , coriolis_param                                  &
                      , lake_albedo                                     &
                      , lake_t_snow                                     &
                      , lake_t_ice                                      &
                      , lake_t_mean                                     &
                      , lake_t_mxl                                      &
                      , lake_shape_factor                               &
                      , lake_h_snow                                     &
                      , lake_h_ice                                      &
                      , lake_h_mxl                                      &
                      , lake_t_sfc                                      &
                      , ts1_lake                                        &
                      , g_dt                                            &
                      , trap_frozen                                     &
                      , trap_unfrozen

   USE ancil_info, ONLY:                                                &
    ssi_pts                                                             &
   ,sea_pts                                                             &
! arrays which JULES needs to allocate
   ,ssi_index                                                           &
   ,sea_index                                                           &
   ,sice_pts_ncat                                                       &
   ,sice_index_ncat                                                     &
   ,fssi                                                                &
   ,sea_frac                                                            &
   ,sice_frac_ncat                                                      &
   ,nsmax

  USE nstypes

! arrays which JULES needs to allocate
      USE prognostics, ONLY :                                           &
        SNOWDEPTH                                                       &
       ,RHO_SNOW_GRND                                                   &
       ,NSNOW                                                           &
       ,SICE                                                            &
       ,SLIQ                                                            &
       ,TSNOW                                                           &
       ,RHO_SNOW                                                        &
       ,RGRAINL

      USE snow_param, ONLY :                                            &
        ds                                                              &
       ,rho_snow_const

! Convective diagnostic output arrays
      USE cv_diagnostic_array_mod, ONLY:                                  &
        cape_out, up_flux, dwn_flux, u_incr_diag_conv, v_incr_diag_conv   &
       ,dubydt_pout ,dvbydt_pout, t_incr_diag_conv                        &
       ,conv_rain_3d ,conv_snow_3d

  USE nvegparm, ONLY :                                                  &
    albsnc_nvg                                                          &
   ,albsnf_nvg

  USE c_0_dg_c, Only :                                                  &
    tm

  USE c_kappai, Only: kappai,de

      USE conv_diag_0a_mod, ONLY: conv_diag_0a
      USE conv_diag_4a_mod, ONLY: conv_diag_4a
      USE conv_diag_5a_mod, ONLY: conv_diag_5a
      USE conv_diag_6a_mod, ONLY: conv_diag_6a

! Copy of arrays from Convection required by SKEB2
! Dev Note: request owners of APP code to store orig arrays in modules 
      USE stochastic_physics_run_mod,  ONLY:                            &
          skeb2_up_flux, skeb2_dwn_flux, skeb2_cape

      USE dust_parameters_mod, ONLY: ndiv, ndivh, l_dust, l_dust_diag

      USE um_input_control_mod,  ONLY:                                   &
           model_domain,                                                 &
           l_bl_tracer_mix,                                              &
           l_aero_classic,  l_sulpc_so2,       l_sulpc_dms,              &
           l_sulpc_nh3,     l_soot,            l_biomass,                &
                            l_co2_interactive, l_co2_emits,              &
           l_q10,                              l_triffid,                &
           l_phenol,        l_trif_eq,                                   &
           phenol_period,   triffid_period,    l_ocff,                   &
           l_nitrate,                                                    &
           l_hydrology

           
      USE ukca_option_mod, ONLY: l_ukca
      USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
      USE cloud_inputs_mod, ONLY: l_acf_cusack, l_pc2, l_pc2_reset,      &
                                  l_acf_brooks,  l_rhcpt
      USE river_inputs_mod, ONLY: l_rivers, l_inland, river_step
      USE eng_corr_inputs_mod, ONLY: l_emcorr
      USE murk_inputs_mod, ONLY: l_murk, l_murk_advect, l_murk_source
      USE cosp_input_mod, ONLY: l_cosp 
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE Field_Types
      USE ereport_mod, ONLY : ereport
      USE PrintStatus_mod
      USE horiz_grid_mod, ONLY : intw_rho2w,intw_u2p,intw_v2p
      USE metric_terms_mod
      USE domain_params
      USE riv_intctl_mod_1A,      ONLY: riv_intctl_1A
      USE riv_intctl_mod_2A,      ONLY: riv_intctl_2A

      USE diagnostics_riv_mod, ONLY: diagnostics_riv

      USE p_to_t_mod,      ONLY: p_to_t
      USE p_to_u_mod,      ONLY: p_to_u
      USE p_to_u_land_mod, ONLY: p_to_u_land
      USE p_to_u_sea_mod,  ONLY: p_to_u_sea
      USE p_to_v_mod,      ONLY: p_to_v
      USE p_to_v_land_mod, ONLY: p_to_v_land
      USE p_to_v_sea_mod,  ONLY: p_to_v_sea
      USE u_to_p_mod,      ONLY: u_to_p
      USE v_to_p_mod,      ONLY: v_to_p

      USE science_fixes_mod, ONLY: l_riverinland_fix

      USE Submodel_Mod
      
      use cable_data_mod, only : cable_control2      
  
      IMPLICIT NONE
!
! Description: This version interfaces to physics schemes in sequence:
!    convection                     (optional)
!    boundary layer
!    convection                     (optional)
!    hydrology
!    river routing                  (optional)
!
!          CALLed after Semi-Lagrangian in atmosphere timestep.
! Method:
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90 
!   This code is written to UMDP3 v8.3 programming standards.
!
! System component: Control Atmos
!
! Declarations:
!
! Global variables:

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
! Need INVERT_OCEAN (hence N->S configuration of ATMOS grid)
! for river routing
      LOGICAL, PARAMETER :: INVERT_OCEAN=.FALSE.

!     File to set heights for screen diagnostics and options for
!     diagnosis.
!
      REAL Z_OBS_TQ,Z_OBS_WIND
      PARAMETER (                                                       &
     & Z_OBS_TQ = 1.5                                                   &
                         ! Height of screen observations of temperature
!                        ! and humidity.
     &,Z_OBS_WIND = 10.0                                                &
                         ! Height of surface wind observations.
     &)
!
      INTEGER, PARAMETER :: IP_ScrnSurfSim = 0
!                           ! Diagnose the screen temperature using
!                           ! pure surface similarity theory
      INTEGER, PARAMETER :: IP_ScrnDecpl1 = 1
!                           ! Diagnose the screen temperature using 
!                           ! surface similarity theory, but allow 
!                           ! decoupling in very stable conditions
!                           ! based on the quasi-equilibrium radiative
!                           ! solution.
      INTEGER, PARAMETER :: IP_ScrnDecpl2 = 2
!                           ! Diagnose the screen temperature using 
!                           ! including transient effects and radiative
!                           ! cooling

! arguments with intent in. ie: input variables.

! Parallel setup variables
      Integer                                                           &
                   ! Size of small halo in j.
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
     &, g_row_length (0:n_proc-1)                                       &
     &, NumCycles                                                       &
                   ! Number of cycles (iterations) for iterative SISL
     &, CycleNo
                   ! Sweep number

! Parameters

! Model dimensions
      Integer                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, n_rows                                                          &
     &, land_points                                                     &
                    ! IN No.of land points being processed, can be 0.
     &, model_levels                                                    &
     &, wet_levels                                                      &
     &, bl_levels                                                       &
     &, dst_levels                                                      &
                    ! number of deep soil temperature levels
     &, dsm_levels                                                      &
                    ! number of deep soil moisture levels
     &, cloud_levels                                                    &
     &, land_ice_points                                                 &
                        ! number of land ice points
     &, soil_points                                                     &
                        ! number of soil points
     &, n_cca_levels                                                    &
                      ! Number of levels for conv cloud
                      ! amount: 1 for 2D, nlevs for 3D.
     &, ntiles                                                          &
                      ! No. of land-surface tiles ( MOSES II )
     &, tr_levels                                                       &
                      ! No. of free tracer levels
     &, first_constant_r_rho_level                                      &
                                   ! 1st rho level on which r constant
     &, nice                                                            &
                      ! No. of sea ice categories
     &, nice_use                                                        &
                      ! No. of sea ice categories used fully in sfc exch  
     &, DIM_CS1, DIM_CS2  ! soil carbon dimensions

! Model switches

      Logical                                                           &
     &  L_regular                                                       &
                     ! True if NOT variable resolution
     &, L_lbc_old                                                       &
                     !  false for new lbc treatment
     &, Ltimer                                                          &
                 ! true then output some timing information
     &, L_dry       ! true if model to be run with no moisture


      Logical                                                           &
     &  L_CAL360                                                        &
                        ! true if using 360 day calender
     &, L_area_cloud                                                    &
                          ! true if using area cloud fraction (ACF)
     &, L_mixing_ratio                                                  &
                          ! Use mixing ratios (if available)
     &, l_combcld_cca0
                          ! Combined cld diag uses Sec 0 convective cld
!
                        
! model parameters
      Real                                                              &
     &  rhcrit(wet_levels)                                              &
                            ! IN Critical relative humidity.
                                  ! the values need to be tuned
                                  ! for the given set of levels.
     &, CO2_MMR
                        ! set equal to co2_start

      Integer                                                           &
     &  tr_vars                                                         &
                          ! number of free tracer variables
     &, tr_ukca                                             
                          ! number of ukca tracer variables

! RIVER routing

      INTEGER                                                           &
     & AOCPL_P_ROWS, AOCPL_ROW_LENGTH                                   &
     &, G_P_FIELD                                                       &
                                  ! IN size of global ATMOS field
     &, G_R_FIELD                                                       &
                                  ! IN Size of global river field
     &, RIVER_ROW_LENGTH                                                &
                                  ! IN local river row length
     &, RIVER_ROWS                                                      &
                                  ! IN local river rows
     &, GLOBAL_RIVER_ROW_LENGTH                                         &
                                  ! IN global river row length
     &, GLOBAL_RIVER_ROWS                                               &
                                        ! IN GLOBAL river rows
     &, A_STEPS_SINCE_RIV         ! IN No. Physics timsteps since last
!                                 ! call to river routing

! Local parameters:
      INTEGER                                                           &
     &       swap_levels                                                &
                                  ! no. levels for SWAPBOUNDS
     &, gather_pe_trip                                                  &
                                  ! pe River routing to be run on
     &, info                                                            &
                              ! Return code from MPP
     &, icode                 ! Return code : 0 Normal Exit : >0 Error
      PARAMETER(swap_levels=1)              ! by definition for A- T

! data to regrid runoff from ATMOS to river routing grid
      REAL                                                              &
     & XPA(AOCPL_ROW_LENGTH+1)                                          &
                               ! IN Atmosphere TP long coordinates
     &,XUA(0:AOCPL_ROW_LENGTH)                                          &
                               ! IN Atmosphere U long coordinates
     &,XVA(AOCPL_ROW_LENGTH+1)                                          &
                               ! IN Atmosphere V long coordinates
     &,YPA(AOCPL_P_ROWS)                                                &
                               ! IN Atmosphere TP lat coordinates
     &,YUA(AOCPL_P_ROWS)                                                &
                               ! IN Atmosphere U lat coordinates
     &,YVA(0:AOCPL_P_ROWS)     ! IN Atmosphere V lat coordinates
! Data to run river routing
      REAL                                                              &
     & TRIVDIR(RIVER_ROW_LENGTH, RIVER_ROWS)                            &
                                              ! IN River direction file
     &,TRIVSEQ(RIVER_ROW_LENGTH, RIVER_ROWS)                            &
                                              ! IN River sequence file
     &,TWATSTOR(RIVER_ROW_LENGTH, RIVER_ROWS)                           &
                                              ! IN/OUT Water storage
!                                             ! file (Kg)
     &,RIVER_VEL                                                        &
                                              ! IN river velocity (m/s)
     &,RIVER_MCOEF                                                      
                                              ! IN meander coefficient

      INTEGER                                                           &
       I_RIVER_VN
                                              ! IN river model type

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
!
       REAL                                                             &
     & r_area(ROW_LENGTH,ROWS),                                         &
!               accumulated areas file
     & r_inext(ROW_LENGTH,ROWS),                                        &
!               x-coordinate of downstream grid point
     & r_jnext(ROW_LENGTH,ROWS),                                        &
!               y-coordinate of downstream grid point
     &  slope(ROW_LENGTH,ROWS),                                         &
!             slopes (not used yet)
     &  flowobs1(ROW_LENGTH,ROWS),                                      &
!             initialisation for flows
     & r_land(ROW_LENGTH,ROWS),                                         &
!             land/river/sea
     & substore(ROW_LENGTH,ROWS),                                       &
!          routing sub_surface store (mm)
     & surfstore(ROW_LENGTH,ROWS),                                      &
!          routing surface store (mm)
     & flowin(ROW_LENGTH,ROWS),                                         &
!          surface lateral inflow (mm)
     & bflowin(ROW_LENGTH,ROWS)
!          sub-surface lateral inflow (mm)
! Diagnostics info
      REAL                                                              &
     & STASHWORK3(*)                                                    &
                         ! STASH workspace for section 3 (b layer)
     &,STASHWORK5(*)                                                    &
                         ! STASH workspace for section 5 (convection)
     &,STASHWORK8(*)                                                    &
                         ! STASH workspace for section 8 (hydrology)
     &,STASHwork9(*)                                                    &
                         ! STASH workspace for section 9 (LS Cloud)
     &,STASHwork19(*)                                                   &
                         ! STASH workspace for section 19 (Veg)
     &,STASHwork26(*)    ! STASH workspace for sect. 26 (River routing)
!
! Data arrays
      REAL :: u      (udims_s%i_start:udims_s%i_end,        &
                      udims_s%j_start:udims_s%j_end,        &
                      udims_s%k_start:udims_s%k_end)
      REAL :: v      (vdims_s%i_start:vdims_s%i_end,        &
                      vdims_s%j_start:vdims_s%j_end,        &
                      vdims_s%k_start:vdims_s%k_end)
      REAL :: w      (wdims_s%i_start:wdims_s%i_end,        &
                      wdims_s%j_start:wdims_s%j_end,        &
                      wdims_s%k_start:wdims_s%k_end)
      REAL :: w_adv  (wdims_l%i_start:wdims_l%i_end,        &
                      wdims_l%j_start:wdims_l%j_end,        &
                      wdims_l%k_start:wdims_l%k_end)
      REAL :: rho_rsq(pdims_s%i_start:pdims_s%i_end,        &
                      pdims_s%j_start:pdims_s%j_end,        &
                      pdims_s%k_start:pdims_s%k_end)
      REAL :: p      (pdims_s%i_start:pdims_s%i_end,        &
                      pdims_s%j_start:pdims_s%j_end,        &
                      pdims_s%k_start:pdims_s%k_end)
      REAL :: p_theta_levels                                &
                     (tdims_s%i_start:tdims_s%i_end,        &
                      tdims_s%j_start:tdims_s%j_end,        &
                      tdims_s%k_start:tdims_s%k_end)
      REAL :: theta  (tdims_s%i_start:tdims_s%i_end,        &
                      tdims_s%j_start:tdims_s%j_end,        &
                      tdims_s%k_start:tdims_s%k_end)
      REAL :: exner_rho_levels                              &
                     (pdims_s%i_start:pdims_s%i_end,        &
                      pdims_s%j_start:pdims_s%j_end,        &
                      pdims_s%k_start:pdims_s%k_end + 1)      
      REAL :: exner_theta_levels                            &
                     (tdims_s%i_start:tdims_s%i_end,        &
                      tdims_s%j_start:tdims_s%j_end,        &
                      tdims_s%k_start:tdims_s%k_end) 

      REAL :: p_star (pdims%i_start:pdims%i_end,            &
                      pdims%j_start:pdims%j_end)

      REAL :: q      (qdims_l%i_start:qdims_l%i_end,        &
                      qdims_l%j_start:qdims_l%j_end,        &
                      qdims_l%k_start:qdims_l%k_end)
      REAL :: qcl    (qdims_l%i_start:qdims_l%i_end,        &
                      qdims_l%j_start:qdims_l%j_end,        &
                      qdims_l%k_start:qdims_l%k_end)       
      REAL :: qcf    (qdims_l%i_start:qdims_l%i_end,        &
                      qdims_l%j_start:qdims_l%j_end,        &
                      qdims_l%k_start:qdims_l%k_end)

      REAL, INTENT (IN) ::                       &
      qrain  (qdims_l%i_start:qdims_l%i_end,     & ! prognostic rain if present
              qdims_l%j_start:qdims_l%j_end,     &
              qdims_l%k_start:qdims_l%k_end)     &
     ,qgraup (qdims_l%i_start:qdims_l%i_end,     & ! prognostic graupel
              qdims_l%j_start:qdims_l%j_end,     & ! if present
              qdims_l%k_start:qdims_l%k_end)     &  
     ,qcf2   (qdims_l%i_start:qdims_l%i_end,     & ! 2nd ice type if present
              qdims_l%j_start:qdims_l%j_end,     &
              qdims_l%k_start:qdims_l%k_end)

      REAL :: unscaled_dry_rho                                   &
                          (pdims_s%i_start:pdims_s%i_end,        &
                           pdims_s%j_start:pdims_s%j_end,        &
                           pdims_s%k_start:pdims_s%k_end)
             ! unscaled dry density

! ancillary arrays and fields required to be saved from timestep to
! timestep.

      Real                                                              &
     &  T_surf(row_length, rows)

      logical                                                           &
     &  land_sea_mask(row_length, rows)

      Integer                                                           &
     &  land_index (land_points)                                        &
                                      ! set from land_sea_mask
     &, land_ice_index (land_points)                                    &
                                      ! Array of land ice points.
     &, soil_index(land_points)       ! Array of soil points.

      REAL :: u_0  (udims%i_start:udims%i_end,                          &
                    udims%j_start:udims%j_end)
                                ! set to zero
      REAL :: v_0  (vdims%i_start:vdims%i_end,                          &
                    vdims%j_start:vdims%j_end)
                                ! set to zero
      REAL :: u_0_p(pdims%i_start:pdims%i_end,                          &
                    pdims%j_start:pdims%j_end)
                                  ! set to zero
      REAL :: v_0_p(pdims%i_start:pdims%i_end,                          &
                    pdims%j_start:pdims%j_end)                                
                                ! set to zero

      REAL ::                                                           &
     &  hcon (land_points)                                              &
                             ! soil/qrparm.soil.hcond
     &, hcap (land_points)                                              &
                             ! soil/qrparm.soil.hcap
     &, smvccl (land_points)                                            &
                             ! soil/qrparm.soil.crit
     &, smvcwt (land_points)                                            &
                             ! soil/qrparm.soil.wilt
     &, smvcst (land_points)                                            &
                             ! soil/qrparm.soil.satn
     &, sthf(land_points,dsm_levels)                                    &
                                      ! IN Frozen soil moisture content
!                of each layer as a fraction of saturation.
     &, sthu(land_points,dsm_levels)                                    &
                                      ! IN Unfrozen soil moisture
!                content of each layer as a fraction of saturation.
     &, canopy_gb (land_points)                                         &
                                ! set to zero.
     &, snow_depth (row_length, rows) ! snow/qrclim.snow.(month)

      Real                                                              &
     &  ice_fract (row_length, rows)                                    &
                                     ! ice/qrclim.ice.(month)
     &, di(row_length, rows)                                            &
                              ! ice/qrclim.ice_thick.(month)
     &, ice_fract_ncat(row_length, rows, nice)                          &
     &, di_ncat(row_length, rows, nice)                                 &
     &, z0msea(row_length, rows)                                        &
                                   ! Sea surface roughness
     &, z0m_scm(row_length, rows)                                       &
                                   ! Fixed sea surface roughness
                                   ! length(m) for MOMENTUM,
                                   ! used in SCM
     &, z0h_scm(row_length, rows)  ! Fixed sea surface roughness
                                   ! length(m) for HEAT,
                                   ! used in SCM

      Real                                                              &
     &  sil_orog_land (land_points)                                     &
                                   ! orog/qrparm.orog.as
     &, ho2r2_orog (land_points)                                        &
                                   ! orog/qrparm.orog.h2root2
     &, sd_orog (land_points)                                           &
                                   ! orog/qrparm.orog.stdev
     &, t_soil(land_points,dsm_levels)                                  &
                                   ! slt/qrclim.slt_pm(lev).(month)
     &, ti_gb(row_length, rows)                                         &
                                   ! sea ice sfc layer temp (ice mean)
     &, ti(row_length, rows, nice)                                      &
                                   ! sea ice sfc layer temp on categories
     &, k_sice_ml(row_length, rows, nice)                                &
                                   ! sea ice effective conductivity in
                                   !  sfc layer on categories (W/m2/K)
                                   !  (only set if l_sice_multilayers=T)
     &, t1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, q1_sd(row_length, rows)                                         &
                                ! set to zero initially
     &, zh(row_length, rows)                                            &
                             ! boundary layer height
      , zh_prev(row_length, rows)                                       & 
                                  ! boundary layer height from previous 
!                                 ! timestep 
      , ddmfx(row_length, rows) 
!                                 ! Convective downdraught mass-flux
!                                 ! at cloud-base 

      Real                                                              &
     &  clapp(land_points)                                              &
                             !  qrparm.soil.bwag ?
!                               Clapp-Hornberger exponent.
     &, satcon(land_points)                                             &
                                 !  qrparm.soil.satcon
     &, sathh(land_points)                                              &
                             !  soil water suction
     &, soil_layer_moisture(land_points, dsm_levels)
                             !  qrclim.smc_pm(lev).(month)

! CLoud fields

! local variables.
      Integer                                                           &
     &  rhc_row_length                                                  &
                        ! Row length for RHcrit array
     &, rhc_rows        ! Row number for RHcrit array

     REAL :: area_cloud_fraction (qdims%i_start:qdims%i_end,        &
                                  qdims%j_start:qdims%j_end,        &
                                              1:qdims%k_end)
     REAL :: bulk_cloud_fraction_halos                              &
                                 (qdims_l%i_start:qdims_l%i_end,    &
                                  qdims_l%j_start:qdims_l%j_end,    &
                                  qdims_l%k_start:qdims_l%k_end)
     REAL :: cloud_fraction_liquid_halos                            &
                                 (qdims_l%i_start:qdims_l%i_end,    &
                                  qdims_l%j_start:qdims_l%j_end,    &
                                  qdims_l%k_start:qdims_l%k_end)
     REAL :: cloud_fraction_frozen_halos                            &
                                 (qdims_l%i_start:qdims_l%i_end,    &
                                  qdims_l%j_start:qdims_l%j_end,    &
                                  qdims_l%k_start:qdims_l%k_end)
     REAL :: rhcpt (rhc_row_length, rhc_rows,1:qdims%k_end)

! Rain fields
      Real                                                              &
     &  ls_rain(row_length, rows)                                       &
     &, ls_snow(row_length, rows)                                       &
     &, micro_tends(row_length, rows, 2, bl_levels)                     &
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)
     &, conv_rain(row_length, rows)                                     &
     &, conv_snow(row_length, rows)

! Radiation fields
      Real                                                              &
     &  photosynth_act_rad(row_length, rows)                            &
                                             ! Net downward
!                                 shortwave radiation in band 1 (w/m2).
      , rad_hr(row_length, rows, 2, bl_levels)
                                               !
!                               ! BL radiative (LW,SW) heating rates

! Fields for mineral dust source flux calculations
      REAL, INTENT(IN) ::                                               &
     &  soil_clay ( row_length, rows )                                  &
     &, soil_silt ( row_length, rows )                                  &
     &, soil_sand ( row_length, rows )                                  &
     &, dust_mrel1 ( row_length, rows )                                 &
     &, dust_mrel2 ( row_length, rows )                                 &
     &, dust_mrel3 ( row_length, rows )                                 &
     &, dust_mrel4 ( row_length, rows )                                 &
     &, dust_mrel5 ( row_length, rows )                                 &
     &, dust_mrel6 ( row_length, rows )

! Tracer emissions
      REAL, INTENT(IN) ::                                               &
     &  so2_hilem ( row_length, rows )                                  &
     &, so2_em    ( row_length, rows )                                  &
     &, nh3_em    ( row_length, rows )                                  &
     &, dms_em    ( row_length, rows )                                  &
     &, soot_hilem( row_length, rows )                                  &
     &, soot_em   ( row_length, rows )                                  &
     &, ocff_hilem( row_length, rows )                                  &
     &, ocff_em   ( row_length, rows )

! CO2 fields
      Real                                                              &
     &  co2_emits ( row_length, rows )                                  &
     &, co2flux ( row_length, rows )

! Fields holding information on past history of convection. Note if 
! prognostics not requested these fields are unset and there is no
! space allocated for them at the top level of the UM.
! Fields are reduced over time or reset if convection occurs.
      Real, Intent(Inout) ::        &
       deep_flag(row_length,rows)   &  ! Value between 0.0 and 1.0
                                       ! 1 if deep convection last timestep
                                       ! Reduced according to a decay period
                                       ! if no deep convection
      ,past_precip(row_length,rows) &  ! Rate of convective precip (kg/m2/s)
                                       ! decayed over time.
      ,past_conv_ht(row_length,rows)   ! Height of convection (m)

! Convection
      Real                                                              &
     &  cca (row_length, rows, n_cca_levels)                            &
     &, cclwp(row_length, rows) ! condensed water path (KG/M**2)

      Integer                                                           &
     &  ccb (row_length, rows)                                          &
     &, cct (row_length, rows)


      Real                                                              &
           ! local vertical co-ordinate information
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
           !VarRes horizontal co-ordinate information
     &  dlambda_p(1-halo_i:row_length+halo_i)                           &
     &, dphi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)        &
     &, wt_lambda_p(1-halo_i:row_length+halo_i)                         &
     &, wt_lambda_u(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)      &
     &, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

      
      REAL,INTENT(In) ::                                                &
     & f3_at_u(1-offx:row_length+offx, 1-offy:rows+offy)

! time information for current timestep
      Integer                                                           &
     &  val_year                                                        &
     &, val_day_number                                                  &
     &, val_hour                                                        &
     &, val_minute                                                      &
     &, val_second

! Diagnostic variables
      Real                                                              &
     &  lat_rot_NP                                                      &
     &, long_rot_NP

! Variables to be used by COSP
!  Convective rainfall
      REAL,INTENT(OUT) :: cosp_crain_3d(row_length,rows,model_levels)
!  Convective snowfall
      REAL,INTENT(OUT) :: cosp_csnow_3d(row_length,rows,model_levels)

! arguments with intent in/out. ie: input variables changed on output.


      REAL :: theta_star(tdims_s%i_start:tdims_s%i_end,        &
                         tdims_s%j_start:tdims_s%j_end,        &
                         tdims_s%k_start:tdims_s%k_end)
      REAL :: q_star    (qdims_s%i_start:qdims_s%i_end,        &
                         qdims_s%j_start:qdims_s%j_end,        &
                         qdims_s%k_start:qdims_s%k_end)
      REAL :: qcl_star  (qdims_s%i_start:qdims_s%i_end,        &
                         qdims_s%j_start:qdims_s%j_end,        &
                         qdims_s%k_start:qdims_s%k_end)
      REAL :: qcf_star  (qdims_s%i_start:qdims_s%i_end,        &
                         qdims_s%j_start:qdims_s%j_end,        &
                         qdims_s%k_start:qdims_s%k_end)

! Extra prognostics - currently not updated by atmos_physics2 so intent in.
! Passed down so convection can do conversions to/from mixing to specific
! accurately if any are present in the run.
      REAL, INTENT (IN) ::                          &
      qrain_star(qdims_s%i_start:qdims_s%i_end,     & ! prognostic rain 
                 qdims_s%j_start:qdims_s%j_end,     & ! if present
                 qdims_s%k_start:qdims_s%k_end)     &
     ,qgraup_star(qdims_s%i_start:qdims_s%i_end,    & ! prognostic graupel
                 qdims_s%j_start:qdims_s%j_end,     & ! if present
                 qdims_s%k_start:qdims_s%k_end)     &  
     ,qcf2_star (qdims_s%i_start:qdims_s%i_end,     & ! 2nd ice type if present
                 qdims_s%j_start:qdims_s%j_end,     &
                 qdims_s%k_start:qdims_s%k_end)

      REAL :: cf_star   (qdims_s%i_start:qdims_s%i_end,        &
                         qdims_s%j_start:qdims_s%j_end,        &
                         qdims_s%k_start:qdims_s%k_end)
      REAL :: cfl_star  (qdims_s%i_start:qdims_s%i_end,        &
                         qdims_s%j_start:qdims_s%j_end,        &
                         qdims_s%k_start:qdims_s%k_end)
      REAL :: cff_star  (qdims_s%i_start:qdims_s%i_end,        &
                         qdims_s%j_start:qdims_s%j_end,        &
                         qdims_s%k_start:qdims_s%k_end)
      REAL, TARGET :: R_u(udims_s%i_start:udims_s%i_end,       &
                          udims_s%j_start:udims_s%j_end,       &
                          udims_s%k_start:udims_s%k_end)
      REAL, TARGET :: R_v(vdims_s%i_start:vdims_s%i_end,       &
                          vdims_s%j_start:vdims_s%j_end,       &
                          vdims_s%k_start:vdims_s%k_end) 
      REAL :: R_w       (wdims%i_start:wdims%i_end,            &
                         wdims%j_start:wdims%j_end,            &
                         wdims%k_start:wdims%k_end)

      REAL :: sum_eng_fluxes(row_length, rows)
      REAL :: sum_moist_flux(row_length, rows)

! Used in SCM for prescribed surface flux forcing
      Real                                                              &
     &  flux_e(row_length, rows)                                        &
                                 ! Surface latent heat flux (W/m^2)
     &, flux_h(row_length, rows) ! Surface sensible heat flux (W/m^2)

!    IN logicals for surface forcing
      LOGICAL                                                           &
     & L_flux_bc                                                        &
                    ! T if prescribed surface fluxes to be used
     &,L_spec_z0    ! T if roughness lengths have been specified

! Tracer variables
      REAL, INTENT(INOUT) ::                                            &
     &  aerosol     ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end)
      REAL, INTENT(INOUT) ::                                            &
     &  free_tracers( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
     &                tr_levels, tr_vars)
      REAL, INTENT(INOUT) ::                                            &
        ukca_tracers( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end,tr_ukca)

      REAL, INTENT(INOUT) ::                                            &
     &  dust_div1   ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  dust_div2   ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  dust_div3   ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  dust_div4   ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  dust_div5   ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  dust_div6   ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )


      REAL, INTENT(INOUT) ::                                            &
     &  so2         ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  dms         ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  so4_aitken  ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  so4_accu    ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  so4_diss    ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  nh3         ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )

      REAL, INTENT(INOUT) ::                                            &
     &  soot_new    ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  soot_aged   ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  soot_cld    ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  bmass_new   ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  bmass_aged  ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  bmass_cld   ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  ocff_new    ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  ocff_aged   ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  ocff_cld    ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  nitr_acc    ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  nitr_diss   ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )
      REAL, INTENT(INOUT) ::                                            &
     &  co2         ( tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end )

! Add definitions for the cariolle scheme.

      REAL, INTENT(INOUT) ::                                            &
     &   ozone_tracer(tdims_s%i_start:tdims_s%i_end,                    &
                      tdims_s%j_start:tdims_s%j_end,                    &
                      tdims_s%k_start:tdims_s%k_end  )         

! arguments with intent out. ie: output variables.

      REAL, TARGET  :: rhokm  (rkmdims%i_start:rkmdims%i_end,           &
                               rkmdims%j_start:rkmdims%j_end,           &
                               rkmdims%k_start:rkmdims%k_end)   
      REAL  :: cH_term(chdims%i_start:chdims%i_end,                     &
                       chdims%j_start:chdims%j_end,                     &
                       chdims%k_start:chdims%k_end)

      Integer                                                           &
     &  ntml (row_length, rows)                                         &
     &, nbdsc(row_length, rows)                                         &
     &, ntdsc(row_length, rows)

      Logical                                                           &
     &  cumulus (row_length, rows) ! bl convection flag

!     Variables for screen-level diagnostics
      REAL, INTENT(INOUT) :: TScrnDcl_SSI(row_length,rows)
!                           !    Decoupled screen-level temperature
!                           !    over sea or sea-ice
      REAL, INTENT(INOUT) :: TScrnDcl_TILE(LAND_POINTS,NTILES)
!                           !    Decoupled screen-level temperature
!                           !    over land tiles
      REAL, INTENT(INOUT) :: tStbTrans(row_length,rows)
!                           !    Time since the transition

      Integer                                                           &
     &  Error_code

      Integer                                                           &
     & CO2_DIM_LEN                                                      &
                                   ! IN Length of a CO2 field row.
     &,CO2_DIM_ROW                                                      &
                                   ! IN Number of CO2 field rows.
     &,LAND_PTS_TRIF                                                    &
                                   ! IN For dimensioning land fields
     &,NPFT_TRIF                                                        &
                                   ! IN For dimensioning PFT fields
!                                  !    available only with TRIFFID.
!                                  !    Set to NPFT when TRIFFID on,
!                                  !    set to 1 when TRIFFID off.
     &,ASTEPS_SINCE_TRIFFID                                             &
                                   ! IN Number of atmospheric
!                                  !    timesteps since last call
!                                  !    to TRIFFID.
     &,A_STEP                                                           
                                   ! IN Atmospheric timestep number.

      Real                                                              &
     & FRAC(LAND_POINTS,NTYPE)                                          &
                                      ! IN Fractions of surface types.
     &,FRAC_DISTURB(LAND_POINTS)                                        &
                                    ! IN Fraction of gridbox in which
!                                   !    vegetation is disturbed.
     &,CANHT_FT(LAND_POINTS,NPFT)                                       &
                                      ! IN Canopy height (m)
     &,LAI_FT(LAND_POINTS,NPFT)                                         &
                                      ! IN Leaf area index
     &,CANOPY(LAND_POINTS,NTILES)                                       &
                                      ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
     &,CATCH(LAND_POINTS,NTILES)                                        &
                                      ! IN Surface/canopy water capacity
!                                  !    of snow-free land tiles (kg/m2).
     &,CATCH_SNOW(LAND_POINTS,NTILES)                                   &
                                   ! IN Coniferous canopy snow capacity
!                                  !    (kg/m2).
     &,SNOW_TILE(LAND_POINTS,NTILES)                                    &
                                      ! IN Lying snow on tiles (kg/m2)
     &,Z0_TILE(LAND_POINTS,NTILES)                                      &
                                      ! IN Tile roughness lengths (m).
      ,z0h_tile_bare(land_points,ntiles)                                &
                                      ! IN Tile thermal roughness 
                                      ! lengths without snow (m).
                                      ! c.f. z0h_tile which is after
                                      ! adjustment for snow.
     &,T_SURF_TILE(LAND_POINTS,NTILES)                                  &
                                      ! IN Surface tile temperatures
     &,INFIL_TILE(LAND_POINTS,NTILES)                                   &
!                          ! IN Maximum surface infiltration
     &,RGRAIN(LAND_POINTS,NTILES)                                       &
                                 ! INOUT Snow grain size (microns).
     &,SNOW_GRND(LAND_POINTS,NTILES)                                    &
                                   ! INOUT Snow below canopy (kg/m2).
     &,GS(LAND_POINTS)                                                  &
                                      ! INOUT "Stomatal" conductance to
!                                  !        evaporation (m/s).
     &,CS(LAND_POINTS,DIM_CS1)                                          &
                                        ! IN Soil carbon (kg C/m2).
     &,G_LEAF_ACC(LAND_POINTS,NPFT)                                     &
                                      ! INOUT Accumulated G_LEAF
     &,G_LEAF_PHEN_ACC(LAND_POINTS,NPFT)                                &
!                                  ! INOUT Accumulated leaf turnover
!                                  !       rate including phenology.
     &,NPP_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                              &
!                                  ! INOUT Accumulated NPP_FT
     &,RESP_W_FT_ACC(LAND_PTS_TRIF,NPFT_TRIF)                           &
!                                  ! INOUT Accum RESP_W_FT
     &,RESP_S_ACC(LAND_PTS_TRIF,DIM_CS1)                                &
                                        ! INOUT Accumulated RESP_S
     !CABLE: kdcorbin, 10/10 
     &,RESP_S_TILE(LAND_POINTS,NTILES)                                  &
                                  ! Soil respiration on tiles (kg C/m2/s).
     &,OLR(ROW_LENGTH,ROWS)                                             &
                                ! IN    TOA - surface upward LW on
!                               !       last radiation timestep
!                               ! OUT   Corrected TOA outward LW
     &,LW_DOWN(ROW_LENGTH,ROWS)                                         &
                                   ! IN Surface downward LW radiation
!                                  !    (W/m2).
     &,SW_TILE(LAND_POINTS,NTILES)    ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
! The following are only used if coastal tiling is switched on:
      REAL                                                              &
     & FLAND_CTILE(LAND_POINTS)                                         &
                                   ! IN Land fraction on land tiles.
     &,TSTAR_LAND_CTILE(ROW_LENGTH,ROWS)                                &
!                                  ! INOUT Land mean sfc temperature (K)
     &,TSTAR_SEA_CTILE(ROW_LENGTH,ROWS)                                 &
!                                  ! IN    Open sea sfc temperature (K).
     &,TSTAR_SICE_CAT_CTILE(ROW_LENGTH,ROWS,NICE_USE)                   &
!                                  ! INOUT Category seaice sfc temp (K).
     &,TSTAR_SICE_CTILE(ROW_LENGTH,ROWS)                                &
!                                  ! INOUT Sea-ice sfc temperature (K).
     &,ALBSOIL(LAND_POINTS)                                             &
                                   ! Soil albedo.
     &,COS_ZENITH_ANGLE(ROW_LENGTH,ROWS)
!                                  ! Cosine of the zenith angle

! INOUT variables for TKE based turbulence schemes
      REAL                                                              &
     &  e_trb(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,          &
     &      bl_levels)                                                  &
!                   ! TKE defined on theta levels K-1
     &, tsq_trb(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,        &
     &      bl_levels)                                                  &
!                   ! Self covariance of liquid potential temperature
!                   ! (thetal'**2) defined on theta levels K-1
     &, qsq_trb(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,        &
     &      bl_levels)                                                  &
!                   ! Self covariance of total water
!                   ! (qw'**2) defined on theta levels K-1
     &, cov_trb(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j,        &
     &      bl_levels)                                                  &
!                   ! Correlation between thetal and qw
!                   ! (thetal'qw') defined on theta levels K-1
     &, zhpar_shcu(row_length, rows)
!                   ! Height of mixed layer used to evaluate
!                   ! the non-gradient buoyancy flux


! Additional variables required for large-scale hydrology:
      REAL                                                              &
     & FEXP(LAND_POINTS)                                                &
                            ! IN Decay factor in Sat. Conductivity
!                           !    in water table layer.
     &,GAMTOT(LAND_POINTS)                                              &
                            ! IN Integrated complete Gamma function.
     &,TI_MEAN(LAND_POINTS)                                             &
                            ! IN Mean topographic index.
     &,TI_SIG(LAND_POINTS)  ! IN Standard dev. of topographic index.
!                           !    in water table layer.
      REAL                                                              &
     & FSAT(LAND_POINTS)                                                &
                            ! INOUT Surface saturation fraction.
     &,FWETL(LAND_POINTS)                                               &
                            ! INOUT Wetland fraction.
     &,ZW(LAND_POINTS)                                                  &
                            ! INOUT Water table depth (m).
     &,STHZW(LAND_POINTS)                                               &
                            ! INOUT soil moist fract. in deep-zw layer.
     &,A_FSAT(LAND_POINTS)                                              &
                            ! IN Fitting parameter for Fsat in LSH model
     &,C_FSAT(LAND_POINTS)                                              &
                            ! IN Fitting parameter for Fsat in LSH model
     &,A_FWET(LAND_POINTS)                                              &
                            ! IN Fitting parameter for Fwet in LSH model
     &,C_FWET(LAND_POINTS)  ! IN Fitting parameter for Fwet in LSH model

      REAL                                                              &
     & DUN_ROFF(LAND_POINTS)                                            &
                            ! OUT Dunne part of sfc runoff (kg/m2/s).
     &,QBASE(LAND_POINTS)                                               &
                            ! OUT Base flow (kg/m2/s).
     &,QBASE_ZW(LAND_POINTS)                                            &
                            ! OUT Base flow from ZW layer (kg/m2/s).
     &,DRAIN(LAND_POINTS)                                               &
                         ! Drainage out of nshyd'th level (kg/m2/s).
     &,FCH4_WETL(LAND_POINTS)
!                           ! OUT Wetland methane flux. (kg C/m2/s).

! Additional variables for SCM diagnostics which are dummy in full UM
      Integer                                                           &
     &  nSCMDpkgs               ! No of SCM diagnostics packages

      Logical                                                           &
     &  L_SCMDiags(nSCMDpkgs)   ! Logicals for SCM diagnostics packages

! loop counters
      Integer                                                           &
     &  i, j, k, l, m, n                                                &
     &, j_begin, j_end

      Logical                                                           &
     &  L_poles   !  include poles in etadot calc (false for LAMs)

! Diagnostic switches
! a) hydrology
      Logical                                                           &
     &  stf_sub_surf_roff                                               &
     &, smlt

! Local parameters:
      CHARACTER (LEN=*), PARAMETER :: RoutineName = 'ATM_PHYSICS2'

! local variables

! Local data arrays


     REAL :: T_latest             (tdims%i_start:tdims%i_end,        &
                                   tdims%j_start:tdims%j_end,        &
                                               1:tdims%k_end) 
     REAL :: theta_conv           (tdims%i_start:tdims%i_end,        &
                                   tdims%j_start:tdims%j_end,        &
                                               1:tdims%k_end)
     REAL :: q_latest             (qdims%i_start:qdims%i_end,        &
                                   qdims%j_start:qdims%j_end,        &
                                               1:qdims%k_end)
     REAL :: qcl_latest           (qdims%i_start:qdims%i_end,        &
                                   qdims%j_start:qdims%j_end,        &
                                               1:qdims%k_end)
     REAL :: qcf_latest           (qdims%i_start:qdims%i_end,        &
                                   qdims%j_start:qdims%j_end,        &
                                               1:qdims%k_end)
     REAL :: cf_latest            (qdims%i_start:qdims%i_end,        &
                                   qdims%j_start:qdims%j_end,        &
                                               1:qdims%k_end)
     REAL :: cfl_latest           (qdims%i_start:qdims%i_end,        &
                                   qdims%j_start:qdims%j_end,        &
                                               1:qdims%k_end)
     REAL :: cff_latest           (qdims%i_start:qdims%i_end,        &
                                   qdims%j_start:qdims%j_end,        &
                                               1:qdims%k_end)
     REAL :: bulk_cloud_fraction  (qdims%i_start:qdims%i_end,        &
                                   qdims%j_start:qdims%j_end,        &
                                               1:qdims%k_end)
     REAL :: cloud_fraction_liquid(qdims%i_start:qdims%i_end,        &
                                   qdims%j_start:qdims%j_end,        &
                                               1:qdims%k_end)
     REAL :: cloud_fraction_frozen(qdims%i_start:qdims%i_end,        &
                                   qdims%j_start:qdims%j_end,        &
                                               1:qdims%k_end)

     REAL :: theta_inc     (tdims%i_start:tdims%i_end,        &
                            tdims%j_start:tdims%j_end,        &
                                        1:tdims%k_end)
     REAL :: q_conv        (qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: qcl_conv      (qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: qcf_conv      (qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: qrain_conv    (qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: qgraup_conv   (qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: qcf2_conv     (qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: cf_liquid_conv(qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: cf_frozen_conv(qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: bulk_cf_conv  (qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: q_inc         (qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: qcl_inc       (qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: qcf_inc       (qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: cf_liquid_inc (qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: cf_frozen_inc (qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)
     REAL :: bulk_cf_inc   (qdims%i_start:qdims%i_end,        &
                            qdims%j_start:qdims%j_end,        &
                                        1:qdims%k_end)

      REAL :: l_s_poles(row_length, first_constant_r_rho_level-1)
      REAL :: l_n_poles(row_length, first_constant_r_rho_level-1)

!-------------------------------------------------------------
      Real                                                              &
     &  p_layer_boundaries(row_length, rows, 0:model_levels)            &
              ! pressure at layer boundaries. Same as p except at
              ! bottom level = pstar, and at top = 0.
     &, p_layer_centres(row_length, rows, 0:model_levels)               &
              ! pressure at layer centres. Same as p_theta_levels
              !except bottom level = pstar, and at top = 0.
     &, exner_layer_boundaries(row_length, rows, 0:model_levels)        &
     &, exner_layer_centres(row_length, rows, 0:model_levels)

! Local arrays used by both conv_diag and convection. Calculations moved
! from conv_diag & convection to save CPU.

      REAL :: z_theta(row_length, rows, model_levels)
                        ! height of theta levels above surface (m)
      REAL :: z_rho  (row_length, rows, model_levels)
                        ! height of rho levels above surface (m)
      REAL :: rho_wet(row_length, rows, model_levels)
                        ! wet density on rho levels (kg/m3)
      REAL :: rho_wet_theta(row_length, rows, model_levels-1)
                        ! wet density on theta levels (not top) (kg/m3)
      REAL :: rho_dry(row_length, rows, model_levels)
                        ! dry density on rho levels (kg/m3)
      Real                                                              &
     &  etadot_copy(row_length,rows,0:model_levels)
! Local arrays holding information to be passed between physics
! routines.

      Real                                                              &
     &  ls_rain_land(land_points)                                       &
     &, ls_snow_land(land_points)                                       &
     &, conv_rain_land(land_points)                                     &
     &, conv_snow_land(land_points)
! Local Arrays for convection
      Integer                                                           &
     &  lcbase (row_length, rows)                                       &
                                     ! lowest convective cloud base leve
     &, lctop (row_length, rows)     ! lowest convective cloud top level

!                                        ! shallow convection
      Real                                                              &
     &  lcca(row_length, rows)                                          &
                                ! lowest convective cloud amount (%)
     &, ccw(row_length, rows, wet_levels)
                                  ! convective cloud liquid water
                                  ! (G/KG) on model levels

      !=================================================================
      ! DO NOT REMOVE OR WRITE TO if l_ccrad=.FALSE.
      !
      ! Intermediary versions of ccw lcbase. These are used when
      ! l_ccrad=.TRUE.. They exist because no space is allocated in the
      ! D1 array in atm_step when l_ccrad=.FALSE., without them
      ! overwriting will occur. They should not be removed unless an
      ! alternative in atm_step can be found.

      Real    :: ccw_out    (row_length, rows, wet_levels)
      Integer :: lcbase_out (row_length, rows)            

      ! DO NOT REMOVE OR WRITE TO if l_ccrad=.FALSE.
      !=================================================================

      REAL :: cca0    (row_length,rows,n_cca_levels)
      REAL :: ccw0    (row_length,rows,wet_levels)
      REAL :: cca0_2d (row_length,rows)
      REAL :: cclwp0  (row_length,rows)

      INTEGER :: ccb0    (row_length,rows)
      INTEGER :: cct0    (row_length,rows)
      INTEGER :: lcbase0 (row_length, rows)

      Integer                                                           &
     &  ntpar(row_length, rows)                                         &
                                     ! top of diagnostic parcel ascent
     &, nlcl(row_length, rows)       ! lifting condensation level

      ! Convective type array :: 
      Integer :: conv_type(row_length, rows) 
                      ! Integer index describing convective type: 
                      !    0=no convection 
                      !    1=non-precipitating shallow 
                      !    2=drizzling shallow 
                      !    3=warm congestus 
                      !    ... 
                      !    8=deep convection 

      Logical ::                      &
      l_shallow(row_length, rows)     & ! Logical indicator of shallow 
                                        ! shallow 
     , l_congestus(row_length, rows)  & ! Logical indicator of congestus 
     , l_congestus2(row_length, rows) & ! Logical indicator of congestus
     , l_pc2_diag_sh_pts(row_length, rows)  ! PC2 is using diagnostic
                                            ! shallow convection

      real                                                              &
     & CIN_undilute(row_length, rows)                                   &
                                       ! undilute CIN from parcel ascent
                                       ! (m2/s2)
     &,CAPE_undilute(row_length, rows) ! undilute CAPE from parcel
                                       ! ascent (m2/s2)

      Real ::                        &
       delthvu(row_length, rows)     & ! buoyancy integral
     , zhpar(row_length, rows)       & ! height of ntpar
     , dzh(row_length, rows)         & ! inversion thickness
     , zlcl(row_length, rows)        & ! height of nlcl
     , zlcl_uv(row_length,rows)      & ! height of nlcl for uv grid
     , ql_ad(row_length,rows)        & ! adiabatic liquid water content
                                       ! at inversion (kg/kg)
     , entrain_coef(row_length,rows) & ! entrainment factor
     , qsat_lcl(row_length,rows)       ! qsat at cloud base (kg/kg)

      Real                                                              &
     &  wstar(row_length, rows)                                         &
                                     ! surface-based mixed layer
!                                    ! velocity scale
     &, wthvs(row_length, rows)                                         &
                                     ! surface buoyancy flux
     &, w_max(row_length, rows) ! max w in column

       
      Real                                                              &
     &  sub_surf_roff(land_points)                                      &
                                    ! sub-surface runoff
     &, snomlt_sub_htf(land_points) ! subsurface snowmelt heat flux

      Real                                                              &
     &  infil(land_points)                                              &
                            ! max surface infiltration rate (kg/m2/s)
     &, snow_melt(land_points)                                          &
                                ! snowmelt (kg/m2/s).
     &, surf_roff(land_points)                                          &
                                ! surface runoff (kg/m2/s).
     &, tot_tfall(land_points)  ! total throughfall (kg/m2/s).

! Fields passed from BDY_LAYR to IMP_SOLVER
      Real                                                              &
      &  rhokh (row_length, rows, bl_levels)

! Additional fields
      REAL                                                              &
     & TSTAR_LAND(ROW_LENGTH,ROWS)                                      &
                                   ! Land mean sfc temperature (K)
     &,TSTAR_SEA(ROW_LENGTH,ROWS)                                       &
                                   ! Open sea sfc temperature (K).
     &,TSTAR_SICE_CAT(ROW_LENGTH,ROWS,NICE_USE)                         &
                                   ! Sea-ice category sfc temperature (K).
     &,TSTAR_SICE(ROW_LENGTH,ROWS)                                      &
                                   ! Sea-ice sfc temperature (K).
     &,TSTAR_SSI(ROW_LENGTH,ROWS)  ! Sea mean sfc temperature (K).
! diagnostics
      Real                                                              &
            ! output from bdy_layr.
     &  e_sea(row_length, rows)                                         &
     &, fqT(row_length, rows, bl_levels)                                &
     &, ftl(row_length, rows, bl_levels)                                &
     &, h_sea(row_length, rows)                                         &
     &, taux(udims%i_start:udims%i_end,udims%j_start:udims%j_end,       &
             bl_levels)                                                 &
     &, tauy(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,       &
             bl_levels)


! data required for tracer mixing :
      Real                                                              &
     &  DRYDEP2(ROW_LENGTH,ROWS,NDIV)
                                      !dry dep though grav. set.

! Fields passed out of boundary layer into hydrology
      Real                                                              &
             !
     &  ecan(row_length, rows)                                          &
                                        !output from sf_evap.
     &, ei(row_length, rows)                                            &
                                        !output from sf_evap.
     &, snowmelt(row_length, rows)                                      &
                                        !output from sf_evap.
     &, ext(land_points,dsm_levels) ! Extraction of water from each
!                                    soil layer (kg/m2/s).

      Real                                                              &
             ! ( needed for soil_hft )
     &  surf_ht_flux_land(row_length, rows)                             &
                                             !
     &, surf_ht_flux_ld(land_points)                                    &
                                             !
     &, snomlt_surf_htf(row_length, rows)

      Real                                                              &
             ! (STASH diagnostic)
     &  LYING_SNOW(LAND_POINTS)  ! Gridbox snowmass (kg/m2).
      Real                                                              &
     &  cca_2d(row_length, rows)                                        &
     &, lclf

      Real                                                              &
     &  tot_precip_scaled(row_length, rows)

! logicals
      Logical                                                           &
     &  L_zero_boundaries                                               &
     &, L_zero_halos
!
      Logical                                                           &
     & L_scrn                                                           &
                                 ! Logical to control output
                                 !    of screen level T,Q,QCL,QCF
     &,L_plsp                    ! Logical to control output
                                 !    of Probability of LS Precip
!
      Logical                                                           &
     &  L_calc_dxek                                                     &
                     ! Switch for calculation of condensate increment
     &, L_q_interact ! Switch allows overwriting of parcel variables
!                      when calculating condensate increments.
      Logical                                                           &
     &  L_cape_opt_345 !Logical to store whether cape_opt has value
                       ! 3, 4, 5 now also 6.

! Additional variables for MOSES II
      Real                                                              &
     & FTL_TILE(LAND_POINTS,NTILES)                                     &
                                      ! Surface FTL for land tiles
     &,LE_TILE(LAND_POINTS,NTILES)                                      &
                                      ! Surface latent heat flux for
!                                  !     land tiles
     &,RADNET_SICE(ROW_LENGTH,ROWS,NICE_USE)                            &
                                   ! Surface net radiation on
!                                  !     sea-ice (W/m2)
     &,RADNET_TILE(LAND_POINTS,NTILES)                                  &
                                      ! Surface net radiation on
!                                  !     land tiles (W/m2)
     &,FQT_TILE(LAND_POINTS,NTILES)                                     &
                                      ! Surface FQT for land tiles
     &,EPOT_TILE(LAND_POINTS,NTILES)                                    &
!                                  ! OUT Local EPOT for land tiles.
     &,FQT_ICE(ROW_LENGTH,ROWS,NICE_USE)                                &
                                   ! Surface FQT for sea-ice
     &,FTL_ICE(ROW_LENGTH,ROWS,NICE_USE)                                &
                                   ! Surface FTL for sea-ice
     &,TAUX_LAND(udims%i_start:udims%i_end,udims%j_start:udims%j_end)   &
                                   ! Taux over land part of grid box.
     &,TAUX_SSI(udims%i_start:udims%i_end,udims%j_start:udims%j_end)    &
                                   ! Taux over sea part of grid box.
     &,TAUY_LAND(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)   &
                                   ! Tauy over land part of grid box.
     &,TAUY_SSI(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)

      REAL                                                              &
     & ESOIL_TILE(LAND_POINTS,NTILES)                                   &
                                ! Evaporation from bare soil (kg/m2)
     &,ES(row_length,rows)                                              &
                                ! Surface evapotranspiration from
!                               !     soil moisture store (kg/m2/s).
     &,EI_TILE(LAND_POINTS,NTILES)                                      &
                                   ! EI for land tiles
     &,Q1P5M_TILE(LAND_POINTS,NTILES)                                   &
!                               ! Q1P5M over land tiles.
     &,T1P5M_TILE(LAND_POINTS,NTILES)                                   &
!                               ! T1P5M over land tiles.
     &,ECAN_TILE(LAND_POINTS,NTILES)                                    &
                                ! ECAN for land tiles
     &,MELT_TILE(LAND_POINTS,NTILES)                                    &
!                               ! Snowmelt on tiles (kg/m2/s).
     &,SURF_HTF_TILE(LAND_POINTS,NTILES)
!                               ! Net downward surface heat flux
!                               ! on tiles (W/m2)
      INTEGER                                                           &
     & PHENOL_CALL                                                      &
                    ! indicates whether phenology is to be called
     &,TRIFFID_CALL                                                     &
                    ! indicates whether TRIFFID is to be called
     &,NSTEP_TRIF   ! Number of atmospheric timesteps between calls to
!                   ! TRIFFID vegetation model.

! Additional variables for JULES

      REAL ::                                                           &
       DTSTAR(ROW_LENGTH,ROWS,NICE_USE)
      REAL ::                                                           &
     & DTSTAR_TILE(LAND_POINTS,NTILES)
                                   ! Change in TSTAR over timestep 
!                                  ! for land tiles
      REAL ::                                                           &
     & SNOWDEPTH_P(    land_points,ntiles)                              &
                             ! Snow depth on ground on tiles (m)
     &,RHO_SNOW_GRND_P(land_points,ntiles)                              &
                             ! Snowpack bulk density (kg/m3)
     &,NSNOW_P(        land_points,ntiles)                              &
                             ! Number of snow layers on ground on tiles
                             ! NOTE that this is converted to an integer.
     &,DS_P(        land_points,ntiles,nsmax)                           &
                             ! Snow layer thickness (m)
     &,SICE_P(      land_points,ntiles,nsmax)                           &
                             ! Snow layer ice mass on tiles (Kg/m2)
     &,SLIQ_P(      land_points,ntiles,nsmax)                           &
                             ! Snow layer liquid mass on tiles (Kg/m2)
     &,TSNOWLAYER_P(land_points,ntiles,nsmax)                           &
                             ! Snow layer temperature (K)
     &,RHO_SNOW_P(  land_points,ntiles,nsmax)                           &
                             ! Snow layer densities (kg/m3)
     &,RGRAINL_P(   land_points,ntiles,nsmax)
                             ! Snow layer grain size on tiles (microns)
! prognostic FLake fields
      REAL ::                                                           &
     &  lake_depth_p( land_points)                                      &
     &, lake_fetch_p( land_points)                                      &
     &, lake_t_mean_p(land_points)                                      &
     &, lake_t_mxl_p( land_points)                                      &
     &, lake_t_ice_p( land_points)                                      &
     &, lake_h_mxl_p( land_points)                                      &
     &, lake_h_ice_p( land_points)                                      &
      , lake_shape_p( land_points)                                      &
      , lake_g_dt_p(  land_points)

! Sea ice fields
      REAL ::                                                           &
       ice_fract_cat_use(row_length, rows, nice_use)                    &
                               ! Sea ice fraction passed to ni_bl_ctl
      ,k_sice(row_length, rows, nice) 
                               ! Sea ice effective conductivity

! River routing:
      INTEGER :: nstep_trip

      REAL :: TOT_SURF_RUNOFF(LAND_POINTS)! Accumulated runoff over
      REAL :: TOT_SUB_RUNOFF(LAND_POINTS) ! river routing timestep
                                         ! (Kg/m2/s)
      REAL :: A_BOXAREAS(ROW_LENGTH, ROWS) ! Atmos.gridbox areas (m2)

      REAL :: RIVEROUT(ROW_LENGTH, ROWS) ! river outflow at
                                         ! seapoints on the ATMOS grid
! Declare local inland basin variables
       REAL                                                             &
     & INLANDOUT_ATMOS(row_length,rows)                                 &
                                                ! INLAND
!                                    BASIN FLOW  ATMOS GRID  kg/m2/s
     &,INLANDOUT_ATM(land_points)         ! INLAND BASIN FLOW
!                                  land points only kg/m2/s

      REAL :: INLANDOUT_RIV(RIVER_ROW_LENGTH,RIVER_ROWS)
      ! inland basin OUTFLOW on the trip grid

      REAL :: BOX_OUTFLOW(RIVER_ROW_LENGTH, RIVER_ROWS)
                                         ! gridbox outflow on river
                                         ! routing grid (Kg/s)
      REAL :: BOX_INFLOW(RIVER_ROW_LENGTH, RIVER_ROWS)
                                         ! gridbox runoff on river
                                         ! routing grid (Kg/s)

      LOGICAL :: INVERT_ATMOS            ! marker whether ATMOS
!                                        ! grid is inverted N/S or not
      LOGICAL :: TRIP_CALL               !If River routing called
      LOGICAL :: FIRST_ROUTING      !.T. on first call river routing
! Water conservation 
! Remove global lake evaporation from soil moisture store  
      REAL :: ACC_LAKE_EVAP(ROW_LENGTH,ROWS) 
!                                       ! Accumulated lake evap over  
!                                       ! river routing timestep 
!                                       ! (Kg/m2) 
 
      real ::                                                           &
     & w_copy(row_length,rows,0:model_levels)  ! copy of w to pass to
                                               ! conv_diag, BL?

! Convection and BL INC diagnostics that
! need to be kept between substeps
! Moved from IMP_CTL2 to enable substepping.
      Real,Dimension(:,:,:),Allocatable ::                              &
       u_incr_diag_bl                                                   &
                                 ! u wind      increment for BL
     &,v_incr_diag_bl                                                   &
                                 ! v wind      increment for BL
     &,T_incr_diag_bl                                                   &
                                 ! temperature increment for BL
     &,q_incr_diag_bl                                                   &
                                 ! humidity    increment for BL
     &,qcl_incr_diag_bl                                                 &
                                 ! cl liq frac increment for BL
     &,qcf_incr_diag_bl                                                 &
                                 ! cl fro frac increment for BL
     &,bulk_cf_incr_diag_bl                                             &
                                 ! cf_l  increment for BL
     &,cf_liquid_incr_diag_bl                                           &
                                 ! cf_f  increment for BL
     &,cf_frozen_incr_diag_bl    ! bcf   increment for BL



! STASHflag switches for increment diagnostics:
      Logical                                                           &
     & l_u_incr_bl                                                      &
                             ! u wind
     &,l_v_incr_bl                                                      &
                             ! v wind
     &,L_T_incr_bl_lsc                                                  &
                             ! T across BL and LS CLD
     &,L_Tl_incr_bl_lsc                                                 &
                             ! Tl across BL (and LS CLD)
     &,L_q_incr_bl_lsc                                                  &
                             ! Q across BL and LS CLD
     &,L_qtl_incr_bl_lsc                                                &
                             ! QT (q+qCL) across BL (and LS CLD)
     &,L_qcl_incr_bl_lsc                                                &
                             ! qCL across BL and LS CLD
     &,L_qcl_incr_bl                                                    &
                             ! qcl across BL
     &,L_q_incr_bl                                                      &
                             ! q across BL
     &,L_T_incr_bl                                                      &
                             ! T across BL
     &,L_qcf_incr_bl_lsc                                                &
                             ! qCF across BL (and LS CLD)
     &,L_cf_incr_bl                                                     &
                             ! tot cl frac increment for BL
     &,L_cfl_incr_bl                                                    &
                             ! liq cl frac increment for BL
     &,L_cff_incr_bl 
                             ! ice cl frac increment for BL

!  Workspace :-

!  Local scalars :-

! Variables for multivariate swapbounds
      INTEGER :: i_field
      TYPE(swapable_field_pointer_type) :: fields_to_swap(7) 

      REAL                                                              &
     &  MAG_VECTOR_NP (model_levels)                                    &
     &, DIR_VECTOR_NP (model_levels)                                    &
     &, MAG_VECTOR_SP (model_levels)                                    &
     &, DIR_VECTOR_SP (model_levels)

      Real                                                              &
     & pptn_rate(row_length,rows)                                       &
                                    ! Total precipitation
                                    ! (convective+large scale) (kg/m2/s)
     &,accum_pptn(row_length,rows)                                      &
                                    ! Accumulated total precip (kg/m2)
     &,tot_rain(row_length,rows)                                        &
                                    ! Total rain (kg/m2/s)
     &,tot_snow(row_length,rows)    ! Total snow (kg/m2/s)
 
! Extra variable for JULES snow scheme
      LOGICAL, PARAMETER :: stf_hf_snow_melt = .TRUE.

  REAL ::                                &
    r_u_p(pdims%i_start:pdims%i_end,     &  ! r_u on P-grid
          pdims%j_start:pdims%j_end,     &
          pdims%k_start:pdims%k_end)     &
   ,r_v_p(pdims%i_start:pdims%i_end,     &  ! r_v on P-grid
          pdims%j_start:pdims%j_end,     &
          pdims%k_start:pdims%k_end)     &
   ,ustar_p(pdims%i_start:pdims%i_end,   &  ! u+r_u on P-grid
            pdims%j_start:pdims%j_end,   &
            pdims%k_start:pdims%k_end)   &
   ,vstar_p(pdims%i_start:pdims%i_end,   &  ! v+r_v on P-grid
            pdims%j_start:pdims%j_end,   &
            pdims%k_start:pdims%k_end)

! Work arrays for convection wind variables
REAL, TARGET, ALLOCATABLE ::             &
  work_dubydt_p(:,:,:)                   & ! work array 
 ,work_dvbydt_p(:,:,:) 
REAL, ALLOCATABLE ::                     &
  dubydt_u(:,:,:)                        & ! du/dt on u grid
 ,dvbydt_v(:,:,:)                        & ! dv/dt on v grid
 ,work_u_halo(:,:,:)                     & ! work array for u
 ,work_v_halo(:,:,:)                       ! work array for v

! New varibles added for message passing
  REAL, TARGET ::                                                       &
         flandfac(pdims_s%i_start:pdims_s%i_end,                        &
                  pdims_s%j_start:pdims_s%j_end),                       &
         fseafac(pdims_s%i_start:pdims_s%i_end,                         &
                pdims_s%j_start:pdims_s%j_end),                         &
         rhokm_land(pdims_s%i_start:pdims_s%i_end,                      &
                    pdims_s%j_start:pdims_s%j_end),                     &
         rhokm_ssi(pdims_s%i_start:pdims_s%i_end,                       &
                   pdims_s%j_start:pdims_s%j_end),                      &
         cdr10m(pdims_s%i_start:pdims_s%i_end,                          &
                pdims_s%j_start:pdims_s%j_end),                         &
         tau_fd_x(pdims_s%i_start:pdims_s%i_end,                        &
                  pdims_s%j_start:pdims_s%j_end,bl_levels),             &
         tau_fd_y(pdims_s%i_start:pdims_s%i_end,                        &
                  pdims_s%j_start:pdims_s%j_end,bl_levels)

  REAL ::                                                               &
    flandfac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),    &
    flandfac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),    &
    fseafac_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),     &
    fseafac_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),     &
    taux_fd_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,      &
              bl_levels),                                               &
    tauy_fd_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,      &
              bl_levels),                                               &
    rhokm_u_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end),  &
    rhokm_u_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end),   &
    rhokm_v_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),  &
    rhokm_v_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)

  LOGICAL :: su10, sv10

! 1A BL Scheme uses rho * gamma for non-gradient stress
  REAL, TARGET, ALLOCATABLE :: rhogamu(:,:,:),                          &
                               rhogamv(:,:,:)
  REAL, ALLOCATABLE :: rhogamu_u(:,:,:),                                &
                       rhogamv_v(:,:,:)
! Other BL schemes use dimensionless function for non-gradient stress
  REAL, TARGET, ALLOCATABLE :: f_ngstress(:,:,:)
  REAL, ALLOCATABLE :: f_ngstress_u(:,:,:),                             &
                       f_ngstress_v(:,:,:)

! Required for convective sub-stepping with rediagnosis

      Logical ::                        &
        cumulus_copy(row_length, rows)  & ! Copy of cumulus from conv_diag
       ,no_cumulus(row_length, rows)      ! Points changed from cumulus to 
                                          ! non-cumulus by BL routine

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

! Error reporting

  CHARACTER(LEN=256)       :: message
  INTEGER                  :: errorstatus

  INTEGER, PARAMETER :: i_river_vn_1A = 1
  INTEGER, PARAMETER :: i_river_vn_2A = 2


! needed for ENDGame etadot calc below
REAL u_at_w,v_at_w

! needed for vatpoles for fv_cos_theta_latitude vs cos_theta_latitude
REAL, POINTER :: xx_cos_theta_latitude (:,:)

      IF (lhook) CALL dr_hook('ATMOS_PHYSICS2',zhook_in,zhook_handle)

IF ( l_vatpoles ) THEN
   xx_cos_theta_latitude => cos_theta_latitude
ELSE 
   xx_cos_theta_latitude => fv_cos_theta_latitude
END IF ! vatpoles

! ----------------------------------------------------------------------
! Section INI. Initialisation of variables.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Communication Section
! ----------------------------------------------------------------------

! allocate arrays if this is the first time AP2 has been called
  IF ( .NOT. ALLOCATED(uhalo) ) THEN
     CALL atmos_physics2_alloc(land_points,ntiles,ndiv,ndivh,npft,ntype,&
               dim_cs1,dim_cs2,dsm_levels,bl_levels,nice_use)
  END IF

! only call on 1st cycle or if not fast running
  IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN

!     Include surface currents in coupling.
!     Was not included for old HadGEM1.
! SURFACE CURRENT INTERPOLATION START-
! interpolate U_0,V_0 onto U_0_P,V_0_P

        UHALO(udims%i_start:udims%i_end,                                &
              udims%j_start:udims%j_end) =                              &
          U_0(udims%i_start:udims%i_end,                                &
              udims%j_start:udims%j_end)

        VHALO(vdims%i_start:vdims%i_end,                                &
              vdims%j_start:vdims%j_end) =                              &
          V_0(vdims%i_start:vdims%i_end,                                &
              vdims%j_start:vdims%j_end)

        i_field = 0

        i_field = i_field + 1
        fields_to_swap(i_field) % field_2d   => uhalo(:,:)
        fields_to_swap(i_field) % field_type = fld_type_u
        fields_to_swap(i_field) % levels     = 1
        fields_to_swap(i_field) % rows       = rows
        fields_to_swap(i_field) % vector     = .TRUE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field_2d   => vhalo(:,:)
        fields_to_swap(i_field) % field_type = fld_type_v
        fields_to_swap(i_field) % levels     = 1
        fields_to_swap(i_field) % rows       = n_rows
        fields_to_swap(i_field) % vector     = .TRUE.

! Need to call swap bounds as halo points not set
! DEPENDS ON: swap_bounds_2d_mv
        CALL swap_bounds_2d_mv(fields_to_swap, i_field, row_length,     &
                                offx, offy)

  END IF !l_quick_ap2

        CALL v_to_p(vhalo,                                              &
                        vdims_s%i_start,vdims_s%i_end,                  & 
                        vdims_s%j_start,vdims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        1,                                              &
                        model_domain,at_extremity,v_0_p)

        CALL u_to_p(uhalo,                                              &
                        udims_s%i_start,udims_s%i_end,                  & 
                        udims_s%j_start,udims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        1,                                              &
                        model_domain,at_extremity,u_0_p)

! only call on 1st cycle or if not fast running
  IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN
!
! Calculate etadot, requires comms for ND only
!
      L_poles = .FALSE.
      j_begin = 1
      j_end = rows
      If (model_domain == mt_Global) Then
        If (at_extremity(PSouth)) j_begin = 2
        If (at_extremity(PNorth)) j_end = rows-1
        L_poles = .TRUE.
      End If ! model_domain == mt_Global

IF (.NOT. l_vatpoles) THEN
! Copied from PEHELMEU2A
! DEPENDS ON: etadot_calc
      Call Etadot_Calc (                                                &
     &                  r_theta_levels, r_rho_levels,                   &
     &                  eta_theta_levels, eta_rho_levels,               &
     &                  u, v, w,                                        &
     &                  sec_theta_latitude,                             &
     &                  row_length, rows, n_rows, model_levels,         &
     &                  delta_lambda, delta_phi,                        &
     &                  dlambda_p, dphi_p,                              &
     &                  wt_lambda_p, wt_lambda_u,                       &
     &                  wt_phi_p, wt_phi_v,                             &
     &                  model_domain, first_constant_r_rho_level,       &
     &               gc_proc_row_group, at_extremity, global_row_length,&
     &                  offx, offy, halo_i, halo_j,                     &
     &                  offx, offy, offx, offy, 0, 0,                   &
     &                  1, row_length, 1, rows, j_begin, j_end,         &
     &                  L_regular, L_poles,                             &
     &                  L_s_poles, l_n_poles, etadot_copy)
END IF ! vatpoles
IF (l_vatpoles) THEN
!----------------------------------------------------------------------
! Initialise etadot: Will be computed by the solver at timestep>1
!----------------------------------------------------------------------
 DO k = 1, model_levels-1
   DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end

       u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*u(i-1,j,k+1) +        &
                                  intw_u2p(i,2)*u(i,j,k+1) ) +        &
                intw_rho2w(k,2)*( intw_u2p(i,1)*u(i-1,j,k) +          &
                                  intw_u2p(i,2)*u(i,j,k) )

       v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*v(i,j-1,k+1) +        &
                                  intw_v2p(j,2)*v(i,j,k+1) ) +        &
                intw_rho2w(k,2)*( intw_v2p(j,1)*v(i,j-1,k) +          &
                                  intw_v2p(j,2)*v(i,j,k) )

       etadot_copy(i,j,k) = ( w(i,j,k)/h3_p_eta(i,j,k) -              &
                          u_at_w*dxi1_xi3(i,j,k)/                     &
                                        h1_p_eta(i,j,k) -             &
                          v_at_w*dxi2_xi3(i,j,k)/                     &
                                        h2_p_eta(i,j,k) ) /           &
                                          deta_xi3_theta(i,j,k)
     END DO
   END DO
 END DO

  etadot_copy(:,:,0) = 0.0
  etadot_copy(:,:,model_levels) = 0.0

END IF ! vatpoles

! u and v winds on all levels copied to p-grid, again only requires
! comms for ND runs
      CALL u_to_p(u,                                                    &
                        udims_s%i_start,udims_s%i_end,                  & 
                        udims_s%j_start,udims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        model_levels,                                   &
                        model_domain,at_extremity,u_p)                 

      CALL v_to_p(v,                                                    &
                        vdims_s%i_start,vdims_s%i_end,                  & 
                        vdims_s%j_start,vdims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        model_levels,                                   &
                        model_domain,at_extremity,v_p)    

IF (.NOT. l_vatpoles) THEN 

          If(model_domain  ==  mt_global) then

! DEPENDS ON: polar_vector_wind_n
            Call Polar_vector_wind_n(                                     &
                               v,                                         &
                               sin_theta_longitude,                       &
                               cos_theta_longitude, row_length,           &
                               n_rows, model_levels, mag_vector_np,       &
                               dir_vector_np, mag_vector_sp,              &
                               dir_vector_sp,                             &
                               offx, offy, global_row_length,             &
                               gc_proc_row_group, at_extremity)
            If (at_extremity(PSouth) ) Then
             Do k = 1, model_levels
              DO I=1,row_length
                V_P(i,1,k) = mag_vector_sp(k)
                U_P(i,1,k) = 0.0
              End Do
             End Do
            End If
            If (at_extremity(PNorth) ) Then
             Do k = 1, model_levels
              DO I=1,row_length
                V_P(i,rows,k) = mag_vector_np(k)
                U_P(i,rows,k) = 0.0
              End Do
             End Do
            End If

          End If  ! model_domain
END IF ! vatpoles

  END IF !l_quick_ap2

! Winds required for convection only if convective momentum being used
      IF (l_mom) THEN     

! Call swapbounds on r_u and r_v.
        i_field = 0
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => r_u(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_u
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .TRUE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => r_v(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_v
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  n_rows
        fields_to_swap(i_field) % vector      =  .TRUE.

! DEPENDS ON: swap_bounds_mv
        CALL swap_bounds_mv(fields_to_swap, i_field, row_length,        &
                             offx, offy)

        IF (l_rediagnosis) THEN
          CALL u_to_p (r_u, udims_s%i_start,udims_s%i_end,              & 
                        udims_s%j_start,udims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        model_levels,                                   &
                        model_domain,at_extremity, r_u_p)

          CALL v_to_p (r_v, vdims_s%i_start,vdims_s%i_end,              & 
                        vdims_s%j_start,vdims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        model_levels,                                   &
                        model_domain,at_extremity, r_v_p)

          IF (.NOT. l_vatpoles) THEN
          ! set polar values of u and v to zero.

          IF (model_domain  ==  mt_global) THEN
            IF (at_extremity(psouth)) THEN
              DO k = 1, tdims%k_end
                DO i = 1, row_length
                  r_u_p(i,1,k) = 0.0
                  r_v_p(i,1,k) = 0.0
                END DO
              END DO
            END IF
            IF (at_extremity(pnorth)) THEN
              DO k = 1, tdims%k_end
                DO i = 1, row_length
                  r_u_p(i,rows,k) = 0.0
                  r_v_p(i,rows,k) = 0.0
                END DO
              END DO
            END IF 
          END IF  ! domain
          END IF ! vatpoles

        ELSE     ! non rediagnosis option
          ALLOCATE(work_u_halo(udims_s%i_start:udims_s%i_end,      & 
                               udims_s%j_start:udims_s%j_end,      &
                               udims_s%k_start:udims_s%k_end) )
          ALLOCATE(work_v_halo(vdims_s%i_start:vdims_s%i_end,      & 
                               vdims_s%j_start:vdims_s%j_end,      &
                               vdims_s%k_start:vdims_s%k_end) )

          ! Work out ustar=u+du on u grid and then interpolate to p_grid

          DO k = udims_s%k_start,udims_s%k_end
            DO j = udims_s%j_start, udims_s%j_end
              DO i = udims_s%i_start,udims_s%i_end
                work_u_halo(i,j,k) = u(i,j,k) + r_u(i,j,k)
              END DO ! i
            END DO ! j
          END DO ! k
          CALL u_to_p (work_u_halo,                                     &   
                        udims_s%i_start,udims_s%i_end,                  & 
                        udims_s%j_start,udims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        model_levels,                                   &
                        model_domain,at_extremity, ustar_p)

          DO k = vdims_s%k_start,vdims_s%k_end
            DO j = vdims_s%j_start, vdims_s%j_end
              DO i = vdims_s%i_start,vdims_s%i_end
                work_v_halo(i,j,k) = v(i,j,k) + r_v(i,j,k)
              END DO ! i
            END DO ! j
          END DO ! k
          CALL v_to_p (work_v_halo,                                     &   
                        vdims_s%i_start,vdims_s%i_end,                  & 
                        vdims_s%j_start,vdims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        model_levels,                                   &
                        model_domain,at_extremity, vstar_p)

          DEALLOCATE(work_v_halo)
          DEALLOCATE(work_u_halo)

          IF (.NOT. l_vatpoles) THEN
          ! set polar values of u and v to zero.
          IF (model_domain  ==  mt_global) THEN
            IF (at_extremity(psouth)) THEN
              DO k = 1, tdims%k_end
                DO i = 1, row_length
                  ustar_p(i,1,k)        = 0.0
                  vstar_p(i,1,k)        = 0.0
                END DO
              END DO
            END IF
            IF (at_extremity(pnorth)) THEN
              DO k = 1, tdims%k_end
                DO i = 1, row_length
                  ustar_p(i,rows,k)     = 0.0
                  vstar_p(i,rows,k)     = 0.0
                END DO
              END DO
            END IF
          END IF
          END IF ! vatpoles
        END IF   ! l_rediagnosis



      END IF     ! l_mom

! Calculate extended halo land fraction
! Surface currents and level 1 winds needed for coastal tiling only

! only call on 1st cycle or if not fast running
  IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN

      IF(L_CTILE)THEN

        DO L=1,LAND_POINTS
          FLAND(L)=FLAND_CTILE(L)
        ENDDO

      ELSE

        DO L=1,LAND_POINTS
           FLAND(L)=1.0
        ENDDO

      ENDIF

!-----------------------------------------------------------------------
! Expand land fraction to global field:
!-----------------------------------------------------------------------
      DO J=1,ROWS
       DO I=1,ROW_LENGTH
           FLANDG(I,J)=0.0
       ENDDO
      ENDDO
      DO L=1,LAND_POINTS
        J=(LAND_INDEX(L)-1)/ROW_LENGTH + 1
        I = LAND_INDEX(L) - (J-1)*ROW_LENGTH
        FLANDG(I,J)=FLAND(L)
      ENDDO



      i_field = 0
      i_field = i_field + 1
      fields_to_swap(i_field) % field_2d   => flandg(:,:)
      fields_to_swap(i_field) % field_type = fld_type_p
      fields_to_swap(i_field) % levels     = 1
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.



      DO j = 1, rows
        DO i = 1, row_length
          u_1_px(i,j) = u_p(i,j,1)
          v_1_px(i,j) = v_p(i,j,1)
          u_0_px(i,j) = u_0_p(i,j)
          v_0_px(i,j) = v_0_p(i,j)
        END DO
      END DO



      IF (l_ctile .AND. buddy_sea == on) THEN

        i_field = i_field + 1
        fields_to_swap(i_field) % field_2d   => u_1_px(:,:)
        fields_to_swap(i_field) % field_type = fld_type_p
        fields_to_swap(i_field) % levels     = 1
        fields_to_swap(i_field) % rows       = rows
        fields_to_swap(i_field) % vector     = .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field_2d   => v_1_px(:,:)
        fields_to_swap(i_field) % field_type = fld_type_p
        fields_to_swap(i_field) % levels     = 1
        fields_to_swap(i_field) % rows       = rows
        fields_to_swap(i_field) % vector     = .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field_2d   => u_0_px(:,:)
        fields_to_swap(i_field) % field_type = fld_type_p
        fields_to_swap(i_field) % levels     = 1
        fields_to_swap(i_field) % rows       = rows
        fields_to_swap(i_field) % vector     = .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field_2d   => v_0_px(:,:)
        fields_to_swap(i_field) % field_type = fld_type_p
        fields_to_swap(i_field) % levels     = 1
        fields_to_swap(i_field) % rows       = rows
        fields_to_swap(i_field) % vector     = .FALSE.

      END IF  ! test on buddy_sea switch

! DEPENDS ON: swap_bounds_2d_mv
      CALL swap_bounds_2d_mv(fields_to_swap, i_field, row_length,       &
                           offx , offy )



  END IF !l_quick_ap2

! ----------------------------------------------------------------------
! End of communication Section
! ----------------------------------------------------------------------








       L_apply_diag = CycleNo == NumCycles


! Hydrology

! Sub-surface runoff must be calculated if running with river routeing
! or model will fail.
      stf_sub_surf_roff = sf(205,8) .or. sf(235,8) .or. l_rivers
      smlt = sf(202,8)





!-----------------------------------------------------------------------
! map JULES prognostics to module prognostics
!-----------------------------------------------------------------------
        SNOWDEPTH     = SNOWDEPTH_P
      IF(NSMAX > 0)THEN
        RHO_SNOW_GRND = RHO_SNOW_GRND_P
        NSNOW         = NSNOW_P
        DS            = DS_P
        SICE          = SICE_P
        SLIQ          = SLIQ_P
        TSNOW         = TSNOWLAYER_P
        RHO_SNOW      = RHO_SNOW_P
        RGRAINL       = RGRAINL_P
      ELSE
        rho_snow_grnd = rho_snow_const
        nsnow         = 0
      ENDIF

      IF ( l_flake_model                         &
           .AND. (.NOT.l_aggregate   )           &
           .AND. (land_points > 0    ) ) THEN
!
! map FLake prognostics to module prognostics
          lake_depth        = lake_depth_p
          lake_fetch        = lake_fetch_p
          lake_t_mean       = lake_t_mean_p
          lake_t_mxl        = lake_t_mxl_p
          lake_t_ice        = lake_t_ice_p
          lake_h_mxl        = lake_h_mxl_p
          lake_h_ice        = lake_h_ice_p
          lake_shape_factor = lake_shape_p
          g_dt              = lake_g_dt_p

! initialise FLake variables needed in surface exchange
        DO l = 1, land_points
! set surface T based on T* of the lake tile
          lake_t_sfc(l) = t_surf_tile(l, lake)

! set the snow depth from the mass and density of snow on the tile,
! BUT depending on the presence of ice
          lake_h_snow(l) = 0.0
          IF(lake_h_ice(l) > 0.0)THEN
            lake_h_snow(l) = snow_tile(l, lake) / rho_snow_const
          END IF

! set the snow surface T based on T* of the lake tile
          lake_t_snow(l) = t_surf_tile(l, lake)
        END DO

      END IF

! check BL solver compatibility with Flake
      IF (l_flake_model .AND. (.NOT.l_us_blsol)) THEN
        l_flake_model = .FALSE.
        IF ( printstatus > PrStatus_Normal ) THEN 
          errorstatus = 42
          message = "FLake is only available with the unconditionally stable BL solver."
          CALL ereport(RoutineName, errorstatus, message)
        END IF
      END IF

! Set convective outputs to zero in case convection scheme not called.
      conv_rain(:,:) = 0.0
      conv_snow(:,:) = 0.0


      ! Copy ccw/lcbase values for Radiation if l_ccrad in use.
      ! See comments at declaration of ccw_out/lcbase_out
      IF (l_ccrad) THEN

        DO k=1, wet_levels
          DO j=1, rows
            DO i=1, row_length
              ccw0(i,j,k) = ccw_out(i,j,k)
            END DO
          END DO
        END DO

        DO j=1, rows
          DO i=1, row_length
            lcbase0(i,j)  = lcbase_out(i,j)
          END DO
        END DO

      END IF ! l_ccrad

      IF (rad_cloud_decay_opt == rad_decay_off) THEN
        ! Zero section 0 convective cloud diagnostics going to radiation
        ccb0(:,:)    = 0
        cct0(:,:)    = 0
        cca0_2d(:,:) = 0.0
        cclwp0(:,:)  = 0.0
        lcbase0(:,:) = 0
        ccw0(:,:,:) = 0.0
        cca0(:,:,:) = 0.0
      END IF ! Rad decay disabled

      ccb(:,:)    = 0
      cct(:,:)    = 0
      cca_2d(:,:) = 0.0
      cclwp(:,:)  = 0.0
      lcbase(:,:) = 0
      lctop(:,:)  = 0
      lcca(:,:)   = 0.0
      ccw(:,:,:) = 0.0
      cca(:,:,:) = 0.0

! Set radiation outputs to zero in case convection scheme not called.
! This code can be removed later.
! set temporary logicals to disabled un-called physics.
      if(.Not. L_bl .and. l_param_conv) then
        write(*,*)' convection on and boundary layer off,'
        write(*,*)' not a sensible choice'
        write(*,*)' will try and run but results may be garbage'
      endif

! set p at layer boundaries.
! NB: some arrays have haloes but are unset, if never used they will
!     be removed.
      Do j = 1, rows
        Do i = 1, row_length
          zh_prev(i,j) = zh(i,j)  ! make a copy of zh, as passed in
          p_layer_boundaries(i,j,0) = p_star(i,j)
          p_layer_centres(i,j,0) = p_star(i,j)
          exner_layer_boundaries(i,j,0) = (p_layer_boundaries(i,j,0)/   &
     &                                     p_zero)**kappa
          exner_layer_centres(i,j,0) = exner_layer_boundaries(i,j,0)

        End Do
      End Do

! Start of OpenMP parallel region

!$OMP  PARALLEL DEFAULT(NONE)                                            &
!$OMP& SHARED(p, p_layer_boundaries, p_layer_centres, p_theta_levels,    &
!$OMP& exner_layer_boundaries, exner_layer_centres, exner_rho_levels,    &
!$OMP& exner_theta_levels, model_levels, rows, row_length,               &
!$OMP& z_theta, z_rho, r_theta_levels, r_rho_levels, rho_wet, rho_dry,   &
!$OMP& rho_rsq, unscaled_dry_rho, wet_levels, rho_wet_theta,             &
!$OMP& bulk_cloud_fraction, bulk_cloud_fraction_halos,                   &
!$OMP& cloud_fraction_liquid, cloud_fraction_liquid_halos,               &
!$OMP& cloud_fraction_frozen, cloud_fraction_frozen_halos, qdims, tdims) &
!$OMP& PRIVATE(i,j,k)

! k=model_levels never appears on the LHS, and is not modified by this loop.
! There are instances of p( , ,k+1) and exner_rho_levels( , ,k+1) in the first
! loop, but those arrays are not modified by the second. Hence nowait valid.

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
!$OMP END DO nowait

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
!$OMP END DO
!Implicit barrier

!
! Variables required by conv_diag and convection_control
! Heights of model levels, density; all without halos
!

!$OMP DO SCHEDULE(STATIC)
      Do k =               1, tdims%k_end
        Do j = tdims%j_start, tdims%j_end
          Do i = tdims%i_start, tdims%i_end
            z_theta(i,j,k) = r_theta_levels(i,j,k)                      &
     &                                   - r_theta_levels(i,j,0)
            z_rho(i,j,k)   = r_rho_levels(i,j,k)                        &
     &                                   - r_theta_levels(i,j,0)
            rho_wet(i,j,k) = rho_rsq(i,j,k) /(r_rho_levels(i,j,k)       &
     &                                    *r_rho_levels(i,j,k))
            rho_dry(i,j,k) = unscaled_dry_rho(i,j,k)

          End Do
        End Do
      End Do
!$OMP END DO 
!Implicit barrier

! density on theta levels (not top level).

! There are no private varibles passed to the following routine. Can therefore
! execute on a single thread to avoid waiting for the slowest thread, and
! any superfluous memory costs within the called subroutine.

!$OMP SINGLE
      call p_to_t(row_length,rows, 0,0,0,0, model_levels-1,             &
     &            r_theta_levels(1:row_length,1:rows,0:model_levels),   &
     &            r_rho_levels(1:row_length,1:rows,1:model_levels),     &
     &            rho_wet, rho_wet_theta)
!$OMP END SINGLE

!
!$OMP DO SCHEDULE(STATIC)
      Do k =             1, qdims%k_end
        Do j = qdims%j_start, qdims%j_end
          Do i = qdims%i_start, qdims%i_end
            bulk_cloud_fraction(i,j,k) =                                &
     &        bulk_cloud_fraction_halos(i,j,k)
            cloud_fraction_liquid(i,j,k) =                              &
     &        cloud_fraction_liquid_halos(i,j,k)
            cloud_fraction_frozen(i,j,k) =                              &
     &        cloud_fraction_frozen_halos(i,j,k)
          End Do
        End Do
      End Do
!$OMP END DO
!

! End of OpenMP parallel region
!$OMP END PARALLEL

      If (model_domain == mt_LAM .and. L_lbc_old ) Then

! zero qcl, qcf and bulk_cloud_fraction on LAM boundaries to avoid
! inconsistencies between cloud fields and cloud amounts.
! As results on boundaries do not affect model results this can be done.

        L_zero_boundaries=.TRUE.
        L_zero_halos=.FALSE.

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QCL,       &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
        CALL ZERO_LATERAL_BOUNDARIES(                                   &
     &   ROW_LENGTH,ROWS,HALO_I,HALO_J,WET_LEVELS,fld_type_p,QCF,       &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

! DEPENDS ON: zero_lateral_boundaries
         CALL ZERO_LATERAL_BOUNDARIES(                                  &
     &   ROW_LENGTH,ROWS,0,0,WET_LEVELS,fld_type_p,                     &
     &   BULK_CLOUD_FRACTION,                                           &
     &   1, AT_EXTREMITY,                                               &
     &   L_zero_boundaries,L_zero_halos)

      End If !  model_domain == mt_LAM .and. L_lbc_old

! copy of w or w_adv for use in conv_diag

    IF (i_convection_vn == i_convection_vn_4a) THEN
! to give bit comparison for HadGEM1 etc currently must continue
! to use w_adv

      IF ( .NOT. l_endgame ) THEN
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)             &
!$OMP& SHARED(model_levels, rows, row_length, w_copy, w_adv)  &
!$OMP& PRIVATE(i,j,k)
      Do k=0,model_levels
        Do j=1,rows
          Do i=1,row_length
            w_copy(i,j,k) = w_adv(i,j,k)
          End do
        End do
      End do
!$OMP END PARALLEL DO
      ELSE
! All new convection or no convection scheme switch to using w

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)              &
!$OMP& SHARED(model_levels, rows, row_length, w_copy, w)       &
!$OMP& PRIVATE(i,j,k)
      DO k=0,model_levels
        DO j=1,rows
          DO i=1,row_length
            w_copy(i,j,k) = w(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
      END IF

    ELSE
! All new convection or no convection scheme switch to using w

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)              &
!$OMP& SHARED(model_levels, rows, row_length, w_copy, w)       &
!$OMP& PRIVATE(i,j,k)
      Do k=0,model_levels
        Do j=1,rows
          Do i=1,row_length
            w_copy(i,j,k) = w(i,j,k)
          End do
        End do
      End do
!$OMP END PARALLEL DO
    END IF
! ----------------------------------------------------------------------
! Set Coastal tiling dependent prognostics:
! ----------------------------------------------------------------------
      IF (l_ctile) THEN

        DO j = 1, rows
          DO i = 1, row_length
            tstar_land(i,j)=tstar_land_ctile(i,j)
            tstar_sea(i,j)=tstar_sea_ctile(i,j)
            tstar_sice(i,j)=tstar_sice_ctile(i,j)
            IF (ice_fract(i,j) <= 0.0) THEN
              tstar_ssi(i,j) = tstar_sea(i,j)
            ELSE
              tstar_ssi(i,j)=ice_fract(i,j)*tstar_sice(i,j)             &
                      +(1.0-ice_fract(i,j))*tstar_sea(i,j)
            END IF
          END DO
        END DO

        DO k = 1, nice_use
          DO j = 1, rows
            DO i = 1, row_length
              tstar_sice_cat(i,j,k)=tstar_sice_cat_ctile(i,j,k)
            ENDDO
          ENDDO
        ENDDO

      ELSE

        DO j = 1, rows
          DO i = 1, row_length
            tstar_land(i,j)=t_surf(i,j)
            tstar_ssi(i,j)=t_surf(i,j)

            IF (.NOT.land_sea_mask(i,j)) THEN
              IF (ice_fract(i,j) <= 0.0) THEN
                tstar_sea(i,j)=t_surf(i,j)
                tstar_sice(i,j)=t_surf(i,j)
              ELSE
                tstar_sea(i,j)=tfs
                tstar_sice(i,j)=(t_surf(i,j)                            &
                  -(1.-ice_fract(i,j))*tstar_sea(i,j))/ice_fract(i,j)
              END IF
            ELSE
              tstar_sea(i,j)=t_surf(i,j)
              tstar_sice(i,j)=t_surf(i,j)
            END IF

          END DO
        END DO

! Note that nice_use=1 if not using coastal tiling
        DO k = 1, nice_use
          DO j = 1, rows
            DO i = 1, row_length
              tstar_sice_cat(i,j,k)=tstar_sice(i,j)
            ENDDO
          ENDDO
        ENDDO

      END IF

! Set up sea ice fraction and conductivity fields depending on science choices
     IF (l_sice_multilayers) THEN  !Note nice_use=nice in this case
       k_sice(:,:,:) = k_sice_ml(:,:,:)  ! Received from sea ice model
     ELSE    ! Note nice may not equal nice_use
       ! Set sea ice effective conductivity to constant value
       k_sice(:,:,:) = 2.0*kappai/de  
     ENDIF

     IF (nice_use > 1) THEN   ! Use categories fully
       ice_fract_cat_use(:,:,:) = ice_fract_ncat(:,:,:)
     ELSE
       ice_fract_cat_use(:,:,1) = ice_fract(:,:)
     ENDIF

! ---------------------------------------------
! Call CONV_DIAG to diagnose convection
! ---------------------------------------------

! latest values needed for substepping/fully sequential BL.

        ! Determine value of L_cape_opt_345
        L_cape_opt_345 = ((cape_opt == 3).or.(cape_opt == 4).or.             &
                          (cape_opt == 5).or.(cape_opt == 6))

        !If necessary, initialise the w_max array.
        If (L_cape_opt_345) Then
           w_max(:,:) = 0.0
        End If

! Start of OpenMP parallel region
!$OMP  PARALLEL DEFAULT(NONE)                                                 &
!$OMP& SHARED(rows, row_length,                                               &
!$OMP& zh_prev, zh, t_latest, theta_star,                                     &
!$OMP& exner_theta_levels, q_latest, q_star, qcl_latest, qcl_star,            &
!$OMP& qcf_latest, qcf_star, cf_latest, cf_star, cfl_latest, cfl_star,        &
!$OMP& cff_latest, cff_star, ntml, ntpar, nlcl, cumulus, l_shallow, delthvu,  &
!$OMP& ql_ad, zhpar, dzh, zlcl, zlcl_uv,                                      &
!$OMP& cape_bottom, L_cape_opt_345, cape_top, w_copy, conv_type,              &
!$OMP& no_cumulus, qdims, tdims)                                              &
!$OMP& PRIVATE(i,j,k)                                                         &
!$OMP& REDUCTION(MAX:w_max)

!$OMP DO SCHEDULE(STATIC)
        Do k =             1, tdims%k_end
          Do j = tdims%j_start, tdims%j_end
            Do i = tdims%i_start, tdims%i_end
              T_latest(i,j,k) = theta_star(i,j,k)                       &
     &                        * exner_theta_levels(i,j,k)
            End Do
          End Do
        End Do
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
        Do k =             1, qdims%k_end
          Do j = qdims%j_start, qdims%j_end
            Do i = qdims%i_start, qdims%i_end
              q_latest(i,j,k) = q_star(i,j,k)
              qcl_latest(i,j,k) = qcl_star(i,j,k)
              qcf_latest(i,j,k) = qcf_star(i,j,k)
              cf_latest(i,j,k)  = cf_star(i,j,k)
              cfl_latest(i,j,k) = cfl_star(i,j,k)
              cff_latest(i,j,k) = cff_star(i,j,k)
            End Do
          End Do
        End Do
!$OMP END DO nowait
!  Initialise  conv_diag output arrays
!$OMP DO SCHEDULE(STATIC)
          Do j = 1, rows
            Do i = 1, row_length
              NTML(i,j)=1
              NTPAR(i,j)=1
              NLCL(i,j)=1
              CUMULUS(i,j)=.FALSE.
              no_cumulus(i,j) = .false.
              L_SHALLOW(i,j)=.FALSE.
              conv_type(i,j)=0
              DELTHVU(i,j)=0.0
              ql_ad(i,j) = 0.0
              ZHPAR(i,j)=0.0
              dzh(i,j) =0.0
              ZLCL(i,j)=0.0
              ZLCL_UV(i,j)=0.0
            End do
          End do
!$OMP END DO nowait

          If(L_cape_opt_345) Then

! Find w_max for each column. The w_max array is initialised just
! before the start of the OpenMP parallel region. 
!
!$OMP DO SCHEDULE(STATIC)
            Do k =  cape_bottom, cape_top
              Do j = 1, rows
                Do i = 1, row_length
                  w_max(i,j) = MAX(w_max(i,j), w_copy(i,j,k))
                End Do
              End Do
            End Do !  k =  cape_bottom, cape_top
!$OMP END DO

          End If  !L_cape_opt_345

! End of OpenMP parallel region
!$OMP END PARALLEL
!
! only call on 1st cycle or if not fast running
          IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN
! Uses beginning of timestep values for Theta, q, qcl, qcf ,u and v etc

            SELECT CASE ( i_convection_vn )
            CASE ( i_convection_vn_0a )

            ! convection scheme is not called, while diagnosis is done here.
            IF (PrintStatus >= Prstatus_diag .AND. mype == 0 ) THEN
              WRITE(6,'(A)') ' No convection scheme called.'
            END IF

            CALL conv_diag_0a(                                          &

! IN Parallel variables
            row_length, rows                                            &

! IN model dimensions.
      ,     bl_levels, model_levels, wet_levels                         &
      ,     land_points                                                 &
      ,     p, p_layer_centres(1,1,1),exner_rho_levels                  &
      ,     rho_wet, rho_wet_theta, z_theta, z_rho                      &

! IN Model switches
      ,     l_mixing_ratio,l_ctile                                      & 
! IN  value of l_extra_call in conv_diag
      ,      .FALSE.                                                    &
      ,     no_cumulus                                                  &

! IN cloud data
      ,     qcf(1:row_length,1:rows,1:qdims%k_end)                      &
      ,     qcl(1:row_length,1:rows,1:qdims%k_end), bulk_cloud_fraction &

! IN everything not covered so far :

      ,     p_star, q(1:row_length,1:rows,1:qdims%k_end)                &
      ,     theta(tdims%i_start:tdims%i_end,                            &
              tdims%j_start:tdims%j_end,1:tdims%k_end)                  &
      ,     exner_theta_levels(tdims%i_start:tdims%i_end,               &
              tdims%j_start:tdims%j_end, 1:tdims%k_end)                 &
      ,     u_p, v_p, u_0_p, v_0_p                                      &
      ,     tstar_land, tstar_sea, tstar_sice, z0msea                   &
      ,     L_flux_bc, flux_e, flux_h, L_spec_z0, z0m_scm, z0h_scm      &
      ,     t_surf, land_sea_mask, flandg, ice_fract, timestep          &
      ,     w_copy, w_max, deep_flag, past_precip, past_conv_ht         &
!
! SCM Diagnostics (dummy values in full UM)
      ,     nSCMDpkgs,L_SCMDiags                                        &

! OUT data required elsewhere in UM system :
      ,     ZH,ZHPAR,dzh,ZLCL,ZLCL_UV,DELTHVU,ql_ad,NTML,NTPAR,NLCL     &
      ,     CUMULUS,L_SHALLOW,l_congestus,l_congestus2                  &
      ,     conv_type                                                   &
      ,     CIN_undilute,CAPE_undilute, wstar, wthvs                    &
      ,     entrain_coef, qsat_lcl                                      & 
      ,     Error_code                                                  &
           )

            CASE ( i_convection_vn_4a )

            CALL conv_diag_4a(                                          &

! IN Parallel variables
            row_length, rows                                            &

! IN model dimensions.
      ,     bl_levels, model_levels, wet_levels                         &
      ,     land_points                                                 &
      ,     p, p_layer_centres(1,1,1),exner_rho_levels                  &
      ,     rho_wet, rho_wet_theta, z_theta, z_rho                      &

! IN Model switches
      ,     l_mixing_ratio,l_ctile                                      & 
! IN  value of l_extra_call in conv_diag
      ,      .FALSE.                                                    &
      ,     no_cumulus                                                  &

! IN cloud data
      ,     qcf(1:row_length,1:rows,1:qdims%k_end)                      &
      ,     qcl(1:row_length,1:rows,1:qdims%k_end), bulk_cloud_fraction &

! IN everything not covered so far :

      ,     p_star, q(1:row_length,1:rows,1:qdims%k_end)                &
      ,     theta(tdims%i_start:tdims%i_end,                            &
              tdims%j_start:tdims%j_end,1:tdims%k_end)                  &
      ,     exner_theta_levels(tdims%i_start:tdims%i_end,               &
              tdims%j_start:tdims%j_end, 1:tdims%k_end)                 &
      ,     u_p, v_p, u_0_p, v_0_p                                      &
      ,     tstar_land, tstar_sea, tstar_sice, z0msea                   &
      ,     L_flux_bc, flux_e, flux_h, L_spec_z0, z0m_scm, z0h_scm      &
      ,     t_surf, land_sea_mask, flandg, ice_fract, timestep          &
      ,     w_copy, w_max, deep_flag, past_precip, past_conv_ht         &
!
! SCM Diagnostics (dummy values in full UM)
      ,     nSCMDpkgs,L_SCMDiags                                        &

! OUT data required elsewhere in UM system :
      ,     ZH,ZHPAR,dzh,ZLCL,ZLCL_UV,DELTHVU,ql_ad,NTML,NTPAR,NLCL     &
      ,     CUMULUS,L_SHALLOW,l_congestus,l_congestus2                  &
      ,     conv_type                                                   &
      ,     CIN_undilute,CAPE_undilute, wstar, wthvs                    &
      ,     entrain_coef, qsat_lcl                                      & 
      ,     Error_code                                                  &
           )


            CASE ( i_convection_vn_5a )

            CALL conv_diag_5a(                                          &

! IN Parallel variables
            row_length, rows                                            &

! IN model dimensions.
      ,     bl_levels, model_levels, wet_levels                         &
      ,     land_points                                                 &
      ,     p, p_layer_centres(1,1,1),exner_rho_levels                  &
      ,     rho_wet, rho_wet_theta, z_theta, z_rho                      &

! IN Model switches
      ,     l_mixing_ratio,l_ctile                                      & 
! IN  value of l_extra_call in conv_diag
      ,      .FALSE.                                                    &
      ,     no_cumulus                                                  &

! IN cloud data
      ,     qcf(1:row_length,1:rows,1:qdims%k_end)                      &
      ,     qcl(1:row_length,1:rows,1:qdims%k_end), bulk_cloud_fraction &

! IN everything not covered so far :

      ,     p_star, q(1:row_length,1:rows,1:qdims%k_end)                &
      ,     theta(tdims%i_start:tdims%i_end,                            &
              tdims%j_start:tdims%j_end,1:tdims%k_end)                  &
      ,     exner_theta_levels(tdims%i_start:tdims%i_end,               &
              tdims%j_start:tdims%j_end, 1:tdims%k_end)                 &
      ,     u_p, v_p, u_0_p, v_0_p                                      &
      ,     tstar_land, tstar_sea, tstar_sice, z0msea                   &
      ,     L_flux_bc, flux_e, flux_h, L_spec_z0, z0m_scm, z0h_scm      &
      ,     t_surf, land_sea_mask, flandg, ice_fract, timestep          &
      ,     w_copy, w_max, deep_flag, past_precip, past_conv_ht         &
!
! SCM Diagnostics (dummy values in full UM)
      ,     nSCMDpkgs,L_SCMDiags                                        &

! OUT data required elsewhere in UM system :
      ,     ZH,ZHPAR,dzh,ZLCL,ZLCL_UV,DELTHVU,ql_ad,NTML,NTPAR,NLCL     &
      ,     CUMULUS,L_SHALLOW,l_congestus,l_congestus2                  &
      ,     conv_type                                                   &
      ,     CIN_undilute,CAPE_undilute, wstar, wthvs                    &
      ,     entrain_coef, qsat_lcl                                      & 
      ,     Error_code                                                  &
           )


            CASE ( i_convection_vn_6a )

            CALL conv_diag_6a(                                          &

! IN Parallel variables
            row_length, rows                                            &

! IN model dimensions.
      ,     bl_levels, model_levels, wet_levels                         &
      ,     land_points                                                 &
      ,     p, p_layer_centres(1,1,1),exner_rho_levels                  &
      ,     rho_wet, rho_wet_theta, z_theta, z_rho                      &

! IN Model switches
      ,     l_mixing_ratio,l_ctile                                      & 
! IN  value of l_extra_call in conv_diag
      ,      .FALSE.                                                    &
      ,     no_cumulus                                                  &

! IN cloud data
      ,     qcf(1:row_length,1:rows,1:qdims%k_end)                      &
      ,     qcl(1:row_length,1:rows,1:qdims%k_end), bulk_cloud_fraction &

! IN everything not covered so far :

      ,     p_star, q(1:row_length,1:rows,1:qdims%k_end)                &
      ,     theta(tdims%i_start:tdims%i_end,                            &
              tdims%j_start:tdims%j_end,1:tdims%k_end)                  &
      ,     exner_theta_levels(tdims%i_start:tdims%i_end,               &
              tdims%j_start:tdims%j_end, 1:tdims%k_end)                 &
      ,     u_p, v_p, u_0_p, v_0_p                                      &
      ,     tstar_land, tstar_sea, tstar_sice, z0msea                   &
      ,     L_flux_bc, flux_e, flux_h, L_spec_z0, z0m_scm, z0h_scm      &
      ,     t_surf, land_sea_mask, flandg, ice_fract, timestep          &
      ,     w_copy, w_max, deep_flag, past_precip, past_conv_ht         &
!
! SCM Diagnostics (dummy values in full UM)
      ,     nSCMDpkgs,L_SCMDiags                                        &

! OUT data required elsewhere in UM system :
      ,     ZH,ZHPAR,dzh,ZLCL,ZLCL_UV,DELTHVU,ql_ad,NTML,NTPAR,NLCL     &
      ,     CUMULUS,L_SHALLOW,l_congestus,l_congestus2                  &
      ,     conv_type                                                   &
      ,     CIN_undilute,CAPE_undilute, wstar, wthvs                    &
      ,     entrain_coef, qsat_lcl                                      & 
      ,     Error_code                                                  &
           )

            CASE DEFAULT ! i_convection_vn

          errorstatus = 10
          WRITE (message,'(A)') 'Convection scheme version value not recognised'
          WRITE (message,'(A,I6)') '   i_convection_vn = ',i_convection_vn
          CALL Ereport ( RoutineName, errorstatus, message)

            END SELECT ! i_convection_vn


            IF (l_quick_ap2) THEN
! save outputs for second EG cycle
! N.B. any new OUT data added to conv_diag will need to be saved here
! and restored below. The array size will also need to be increased in
! atmos_physics2_alloc.
               conv_diag_reals(:,:,1)=zh
               conv_diag_reals(:,:,2)=zhpar
               conv_diag_reals(:,:,3)=zlcl
               conv_diag_reals(:,:,4)=zlcl_uv
               conv_diag_reals(:,:,5)=delthvu
               conv_diag_reals(:,:,6)=ql_ad
               conv_diag_reals(:,:,7)=cin_undilute
               conv_diag_reals(:,:,8)=cape_undilute
               conv_diag_reals(:,:,9)=wstar
               conv_diag_reals(:,:,10)=wthvs
               conv_diag_reals(:,:,11)=entrain_coef
               conv_diag_reals(:,:,12)=qsat_lcl
               conv_diag_reals(:,:,13)=dzh

               conv_diag_ints(:,:,1)=ntml
               conv_diag_ints(:,:,2)=ntpar
               conv_diag_ints(:,:,3)=nlcl
               conv_diag_ints(:,:,4)=conv_type

               conv_diag_logs(:,:,1)=cumulus
               conv_diag_logs(:,:,2)=l_shallow
               conv_diag_logs(:,:,3)=l_congestus
               conv_diag_logs(:,:,4)=l_congestus2
            ENDIF

         ELSE
! restore outputs on second EG cycle
            zh=conv_diag_reals(:,:,1)
            zhpar=conv_diag_reals(:,:,2)
            zlcl=conv_diag_reals(:,:,3)
            zlcl_uv=conv_diag_reals(:,:,4)
            delthvu=conv_diag_reals(:,:,5)
            ql_ad=conv_diag_reals(:,:,6)
            cin_undilute=conv_diag_reals(:,:,7)
            cape_undilute=conv_diag_reals(:,:,8)
            wstar=conv_diag_reals(:,:,9)
            wthvs=conv_diag_reals(:,:,10)
            entrain_coef=conv_diag_reals(:,:,11)
            qsat_lcl=conv_diag_reals(:,:,12)
            dzh=conv_diag_reals(:,:,13)

            ntml=conv_diag_ints(:,:,1)
            ntpar=conv_diag_ints(:,:,2)
            nlcl=conv_diag_ints(:,:,3)
            conv_type=conv_diag_ints(:,:,4)

            cumulus=conv_diag_logs(:,:,1)
            l_shallow=conv_diag_logs(:,:,2)
            l_congestus=conv_diag_logs(:,:,3)
            l_congestus2=conv_diag_logs(:,:,4)
         END IF !l_quick_ap2


! Make a copy of the cumulus array so that can tell which points are altered
! by the boundary layer scheme.

          Do j=1,rows
            Do i=1,row_length
              cumulus_copy(i,j) = cumulus(i,j)       
            End Do
          End Do

!----------------------------------------------------------------------
! Section BL Call Explicit part of Boundary Layer scheme.
! ---------------------------------------------------------------------

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Boundary Layer',5)
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Explicit BL',5)

      Do j = 1, rows
        Do i = 1, row_length
! Initialise arrays which will be passed from BL to convection
           wstar(i,j) = 0.0
           wthvs(i,j) = 0.0
        End do
      End do

! DEPENDS ON: Generate_Anthropogenic_Heat
      Call Generate_Anthropogenic_Heat(                                 &
     & val_year, val_day_number, val_hour, val_minute, val_second       &
     &, ntiles, land_points, frac, l_anthrop_heat_src                   &
     & )
!

      L_plsp=sf(281,3) .or. sf(282,3) .or. sf(283,3)
      L_scrn=L_plsp





   call cable_control2( &
                     npft, &
                     tile_frac, &
                     snow_tile, & 
                     vshr_land, & 
                     canopy, &
                     canht_ft, &
                     lai_ft, &
                     conv_rain, &
                     conv_snow, &
                     NPP,&
                     NPP_FT, &
                     GPP,&
                     GPP_FT,&
                     RESP_S,&
                     RESP_S_TOT,&
                     RESP_S_TILE,     &
                     RESP_P,&
                     RESP_P_FT, &
                     G_LEAF, &
                     Radnet_TILE,     &
                     Lying_snow,       &
                     surf_roff, &
                     sub_surf_roff, &
                     tot_tfall, t_soil &
          )

  IF (L_bl) THEN
!
! only call on 1st cycle or if not fast running
    IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN
! set local switch for whether model is using CLASSIC aerosol scheme
      l_aero_classic = (l_sulpc_so2 .OR. l_soot .or. l_biomass .OR.     &
                        l_ocff .OR. l_nitrate .OR. l_dust)

! allocate variables depending on which bl scheme is used







  ALLOCATE(rhogamu(1,1,1))
  ALLOCATE(rhogamv(1,1,1))
  ALLOCATE(f_ngstress(pdims_s%i_start:pdims_s%i_end,                    &
                      pdims_s%j_start:pdims_s%j_end,2:bl_levels))


! DEPENDS ON: ni_bl_ctl
      CALL NI_bl_ctl (                                                  &
! IN parameters for SISL scheme
 NumCycles,CycleNo,                                                     &
! IN model dimensions.
 land_points, ntiles, bl_levels,dst_levels, dsm_levels, nice_use,       &
! IN switches
 L_mixing_ratio, L_scrn, L_aero_classic, L_dust, L_dust_diag,           & 
! IN model Parameters
 CO2_MMR,                                                               &
! IN data fields.
 p, p_layer_centres, rho_rsq,rho_wet,rho_dry, u_p, v_p,                 &
 u_1_px, v_1_px, u_0_px, v_0_px,                                        &
 land_sea_mask, q, qcl, qcf, p_star, theta, EXNER_THETA_LEVELS, RAD_HR, &
 MICRO_TENDS, SOIL_LAYER_MOISTURE, RHO_WET_THETA, Z_RHO, Z_THETA,       &
! IN ancillary fields and fields needed to be kept from tstep to tstep
 hcon,smvccl_levs, smvcwt_levs, smvcst_levs, sthf,sthu,sil_orog_land,   &
!-----------------------------------------------------------------------
 ho2r2_orog, sd_orog, ice_fract_cat_use, k_sice(:,:,1:nice_use),        &
 land_index, photosynth_act_rad,                                        &
 SOIL_CLAY,SOIL_SAND,DUST_MREL1,DUST_MREL2,                             &
 DUST_MREL3,DUST_MREL4,DUST_MREL5,DUST_MREL6,                           &
! IN additional variables for JULES
 CANOPY,CATCH,CATCH_SNOW,SNOW_TILE,Z0_TILE,Z0H_TILE_BARE,               &
 LW_DOWN,SW_TILE,T_SURF_TILE,                                           &
 CO2(1:CO2_DIM_LEN,1:CO2_DIM_ROW,1),L_CO2_INTERACTIVE,                  &
 L_PHENOL,L_TRIFFID,ASTEPS_SINCE_TRIFFID,                               &
 CS,FRAC,CANHT_FT,LAI_FT,FLAND,FLANDG,ALBSOIL,COS_ZENITH_ANGLE,         &
! IN everything not covered so far
 t_soil,ti_gb,t_surf,zh_prev,ddmfx,bulk_cloud_fraction,nlcl,zhpar,zlcl, &
! IN SCM namelist data
 L_spec_z0, z0m_scm, z0h_scm, flux_e, flux_h, L_flux_bc,                &
! SCM diagnostics and STASH
 nSCMDpkgs, L_SCMDiags, BL_diag,                                        &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
! INOUT data
 gs,z0msea,w_copy,etadot_copy,TSTAR_SEA,TSTAR_SICE_CAT,zh,cumulus,      &
 ntpar,l_shallow,error_code,                                            &
! INOUT additional variables for JULES
 G_LEAF_ACC,NPP_FT_ACC,RESP_W_FT_ACC,RESP_S_ACC,l_q10,                  &
! INOUT variables for TKE based turbulence schemes
 e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                          &
! OUT variables for message passing
 flandfac, fseafac,rhokm_land, rhokm_ssi, cdr10m, tau_fd_x, tau_fd_y,   &
 rhogamu, rhogamv, f_ngstress,                                          &
! OUT variables required in IMP_SOLVER
 alpha1_sice, ashtf, BQ_GB, BT_GB, dtrdz_charney_grid, rdz_charney_grid,&
 dtrdz_u, dtrdz_v, rdz_u, rdz_v, z1_tq, uStarGBM,                       &
! OUT diagnostics (done after implicit solver)
 e_sea, fqT, ftl, h_sea, rib_gb, vshr, zht, shallowc, cu_over_orog,     &
 bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6,           &
 bl_type_7, z0m_eff_gb, z0h_eff_gb, fme,                                &
! OUT diagnostics required for soil moisture nudging scheme :
 WT_EXT,RA,                                                             &
! OUT data required for tracer mixing :
 rho_aresist,aresist,resist_b,                                          &
!OUT variables required for mineral dust scheme
 R_B_DUST,DUST_FLUX,DUST_EMISS_FRAC,                                    &
 U_S_T_TILE,U_S_T_DRY_TILE,U_S_STD_TILE, KENT, WE_LIM, T_FRAC, ZRZI,    &
 KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC, ZHSC,                      &
! OUT additional variables for JULES
 FTL_TILE,LE_TILE,RADNET_SICE,RADNET_TILE,RIB_TILE,RHO_ARESIST_TILE,    &
 ARESIST_TILE,RESIST_B_TILE,ALPHA1,ASHTF_TILE,FQT_TILE,EPOT_TILE,       &
 FQT_ICE,FTL_ICE,FRACA,RESFS,RESFT,RHOKH_TILE,RHOKH_SICE,               &
 Z0HSSI,Z0H_TILE,Z0M_GB,Z0MSSI,Z0M_TILE,CHR1P5M,CHR1P5M_SICE,SMC,       &
 GPP,NPP,RESP_P,G_LEAF,GPP_FT,NPP_FT,RESP_P_FT,RESP_S,RESP_S_TOT,       &
 RESP_W_FT,GC,CANHC_TILE,WT_EXT_TILE,FLAKE,TILE_INDEX,TILE_PTS,         &
 TILE_FRAC,FSMC,rib_ssi, vshr_land,vshr_ssi,tstar_land,tstar_ssi,       &
 RHOKH_MIX,DTSTAR_TILE,DTSTAR,HCONS,EMIS_TILE,EMIS_SOIL,                &
! OUT fields
 t1_sd,q1_sd,ntml,nbdsc,ntdsc,wstar,wthvs,uw0,vw0,rhokm,rhokh           &
 )

      IF (l_quick_ap2) THEN
! save outputs for second EG cycle
! N.B. any new OUT data added to ni_bl_ctl will need to be saved here
! and restored below. The array size will also need to be increased in
! atmos_physics2_alloc
         bl_ctl_2d(:,:,1)=zh
         bl_ctl_2d(:,:,2)=z0msea
         bl_ctl_2d(:,:,3)=e_sea
         bl_ctl_2d(:,:,4)=h_sea
         bl_ctl_2d(:,:,5)=wstar
         bl_ctl_2d(:,:,6)=wthvs

         bl_ctl_int2d(:,:,1)=ntml
     
         bl_ctl_log2d(:,:,1)=cumulus
         bl_ctl_log2d(:,:,2)=l_shallow

         bl_ctl_3d(:,:,:,1)=fqt
         bl_ctl_3d(:,:,:,2)=ftl
         bl_ctl_3d(:,:,:,3)=rhokh


         sice_save(:,:,:,1)=radnet_sice
         sice_save(:,:,:,2)=fqt_ice
         sice_save(:,:,:,3)=ftl_ice
         sice_save(:,:,:,4)=dtstar

         rhokm_save=rhokm

         tile_save(:,:,1)=ftl_tile
         tile_save(:,:,2)=le_tile
         tile_save(:,:,3)=radnet_tile
         tile_save(:,:,4)=fqt_tile
         tile_save(:,:,5)=epot_tile
         tile_save(:,:,6)=dtstar_tile

         land_save(:,1)=gs

      END IF

   ELSE
! restore outputs on second EG cycle
      zh     = bl_ctl_2d(:,:,1)
      z0msea = bl_ctl_2d(:,:,2)
      e_sea  = bl_ctl_2d(:,:,3)
      h_sea  = bl_ctl_2d(:,:,4)
      wstar  = bl_ctl_2d(:,:,5)
      wthvs  = bl_ctl_2d(:,:,6)

      ntml = bl_ctl_int2d(:,:,1)
     
      cumulus   = bl_ctl_log2d(:,:,1)
      l_shallow = bl_ctl_log2d(:,:,2)

      fqt         = bl_ctl_3d(:,:,:,1)
      ftl         = bl_ctl_3d(:,:,:,2)
      rhokh       = bl_ctl_3d(:,:,:,3)

      radnet_sice = sice_save(:,:,:,1)
      fqt_ice     = sice_save(:,:,:,2)
      ftl_ice     = sice_save(:,:,:,3)
      dtstar      = sice_save(:,:,:,4)

      rhokm = rhokm_save

      ftl_tile    = tile_save(:,:,1)
      le_tile     = tile_save(:,:,2)
      radnet_tile = tile_save(:,:,3)
      fqt_tile    = tile_save(:,:,4)
      epot_tile   = tile_save(:,:,5)
      dtstar_tile = tile_save(:,:,6)

      gs = land_save(:,1)

   END IF !l_quick_ap2

  END IF !L_bl
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Explicit BL',6)
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Boundary Layer',6)

!------------------------------------------------------------------------
! Which cumulus points has BL changed?
!  cumulus_copy - cumulus as diagnosed by conv_diag
!  cumulus      - cumulus array after BL 
!  no_cumulus   - .true. for those points which have been changed from
!                  .true. to .false.
! Note the BL scheme has good reasons for changing the diagnosis of
! cumulus to .false. for a small number of locations. The information
! used to make these changes is not available to the convection scheme.
!------------------------------------------------------------------------
      Do j= 1,rows
        Do i=1,row_length
          no_cumulus(i,j) = .not.cumulus(i,j) .and. cumulus_copy(i,j)  
        End Do
      End Do

! PC2 required mask
    IF (i_convection_vn == i_convection_vn_5a .OR.   &
        i_convection_vn == i_convection_vn_6a ) THEN
          ! Calculate for each point whether we want to carry
          ! convective cloud information for shallow conv in PC2
      IF ( l_pc2_diag_sh) THEN
        ! Carry convective cloud information
        Do j = 1, rows
          Do i = 1, row_length
            l_pc2_diag_sh_pts(i,j) = l_shallow(i,j)
          End Do
        End Do
      Else
        ! Do not carry convective cloud information
        Do j = 1, rows
          Do i = 1, row_length
            l_pc2_diag_sh_pts(i,j) = .false.
          End Do
        End Do
      Endif

    ELSE
          ! Calculate for each point whether we want to carry
          ! convective cloud information for shallow conv in PC2
      Do j = 1, rows
        Do i = 1, row_length
          IF (l_pc2_diag_sh .AND. l_shallow(i,j) .AND.                &
               bl_type_6(i,j) == 1.0) THEN
            ! Carry convective cloud information
            l_pc2_diag_sh_pts(i,j) = .true.
          Else
            ! Do not carry convective cloud information
            l_pc2_diag_sh_pts(i,j) = .false.
          End if
        End Do
      End Do
    END IF     
! ----------------------------------------------------------------------
! Section CNV.1 Call Convection scheme.
! ----------------------------------------------------------------------

      IF &
       (l_param_conv) &
       THEN
      If (Error_code  ==  0) Then
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Convection',5)

!
! Set up logicals for condensate increment calculation.
! The PC2 Convection code is implemented in three phases:
!  1. OFF-OFF No extra diagnostic space reserved, no calculations.
!  2.  ON-OFF NEW diagnostics calculated, no existing fields touched.
!  3.  ON-ON  Full prognostic interactions and accompanying diagnostics.
!
! Other parts of the PC2 code use L_PC2 to perform overwriting calcs,
! then L_PC2_RESET to write back the original fields for non-interacting
! mode if required. Hence use of different logicals here.
!
          L_calc_dxek  = ( L_pc2 )
          L_q_interact = ( L_pc2  .AND.  ( .NOT. L_pc2_reset ) )
!
!
! Store desired potential temperature in theta_conv array

! Start of OpenMP parallel region
!$OMP  PARALLEL DEFAULT(NONE)                                         &
!$OMP& SHARED(model_levels, rows, row_length, theta_conv, theta_star, &
!$OMP& theta_inc, q_conv, q_star, q_inc, qcl_conv, qcl_star, qcl_inc, &
!$OMP& qrain_star, qrain_conv, qgraup_star, qgraup_conv, qcf2_star,   &
!$OMP& qcf2_conv, l_mcr_qrain, l_mcr_qgraup, l_mcr_qcf2,              &
!$OMP& cf_liquid_conv, cf_liquid_inc, cf_frozen_conv, cf_frozen_inc,  &
!$OMP& cf_star, cfl_star, cff_star, bulk_cf_conv, bulk_cf_inc,        &
!$OMP& wet_levels, qcf_conv, qcf_inc, qcf_star, tdims, qdims)         &
!$OMP& PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
          Do k =             1, tdims%k_end
            Do j = tdims%j_start, tdims%j_end
              Do i = tdims%i_start, tdims%i_end
                theta_conv(i,j,k) = Theta_star(i,j,k)
                theta_inc(i,j,k) = 0.0
              End Do
            End Do
          End Do
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
          Do k =             1, qdims%k_end
            Do j = qdims%j_start, qdims%j_end
              Do i = qdims%i_start, qdims%i_end
                q_conv(i,j,k) = q_star(i,j,k)
                q_inc(i,j,k) = 0.0
                   qcl_conv(i,j,k) = qcl_star(i,j,k)
                qcl_inc(i,j,k)        = 0.0
                   qcf_conv(i,j,k) = qcf_star(i,j,k)
                qcf_inc(i,j,k)        = 0.0
                cf_liquid_conv(i,j,k) = cfl_star(i,j,k)
                cf_liquid_inc(i,j,k)  = 0.0
                cf_frozen_conv(i,j,k) = cff_star(i,j,k)
                cf_frozen_inc(i,j,k)  = 0.0
                bulk_cf_conv(i,j,k)   = cf_star(i,j,k)
                bulk_cf_inc(i,j,k)    = 0.0
              End Do
            End Do
          End Do
!$OMP END DO
          ! Only copy extra prognostics if present. 

          ! Version of array without halos - no inc array as convection
          ! does not update this variable
          IF (l_mcr_qrain) THEN
!$OMP DO SCHEDULE(STATIC)
            DO k =             1, qdims%k_end
              DO j = qdims%j_start, qdims%j_end
                DO i = qdims%i_start, qdims%i_end
                  qrain_conv(i,j,k) = qrain_star(i,j,k)
                END DO
              END DO
            END DO
!$OMP END DO
          END IF
         ! Version of array without halos - no inc array as convection
          ! does not update this variable
          IF (l_mcr_qgraup) THEN
!$OMP DO SCHEDULE(STATIC)
            DO k =             1, qdims%k_end
              DO j = qdims%j_start, qdims%j_end
                DO i = qdims%i_start, qdims%i_end
                  qgraup_conv(i,j,k) = qgraup_star(i,j,k)
                END DO
              END DO
            END DO
!$OMP END DO
          END IF
          ! Version of array without halos - no inc array as convection
          ! does not update this variable
          IF (l_mcr_qcf2) THEN
!$OMP DO SCHEDULE(STATIC)
            DO k =             1, qdims%k_end
              DO j = qdims%j_start, qdims%j_end
                DO i = qdims%i_start, qdims%i_end
                  qcf2_conv(i,j,k) = qcf2_star(i,j,k)
                END DO
              END DO
            END DO
!$OMP END DO
          END IF


! End of OpenMP parallel region
!$OMP END PARALLEL

!-----------------------------------------------------------------------
! Setup flags to control convection output diagnostics
!-----------------------------------------------------------------------

! Apply diags at last cycle only

  CALL set_convection_output_flags(nsects,nitems,l_calc_dxek,l_cosp,sf)


! Convection diagnostics needed for SKEB2
  flg_up_flx       = .TRUE.
  flg_dwn_flx      = .TRUE.

!---------------------------------------------------------------------
! Intercept values of physics increments before being updated by
! convection for optional output of convective increments
! Allocate arrays required for convective diagnostics
!---------------------------------------------------------------------

! DEPENDS ON: cv_alloc_diag_array 
        CALL cv_alloc_diag_array( row_length, rows,                          &
                               l_calc_dxek, l_cosp,                          &
                               l_dust, l_sulpc_so2, l_sulpc_nh3, l_soot,     &
                               l_biomass, l_ocff, l_nitrate,                 &
! Stash info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
                               ntml,ntpar,                                   &
                               theta_inc, q_inc, qcl_inc, qcf_inc ,          &
                               cf_liquid_inc, cf_frozen_inc, bulk_cf_inc,    &
                               conv_rain, conv_snow)

!-----------------------------------------------------------------------

! DEPENDS ON: ni_conv_ctl
        Call NI_conv_ctl (                                              &
! Parallel variables
        at_extremity, n_proc, mype                                      &
! parameters for cycling physics-dynamics
      , NumCycles, CycleNo                                              &
! model dimensions.
      , row_length, rows                                                &
      , rows*row_length                                                 &
      , model_levels, wet_levels, bl_levels, n_cca_levels               &
      , tr_levels, tr_vars, tr_ukca                                     &
      , land_points                                                     &

! Model switches
      , model_domain                                                    &
      , L_calc_dxek, L_q_interact, l_mixing_ratio                       &
      , l_mcr_qrain, l_mcr_qgraup, l_mcr_qcf2                           &
      , l_ctile, L_MURK_SOURCE, L_DUST, L_SULPC_SO2                     &
      , L_sulpc_dms, L_sulpc_nh3, L_soot, L_ocff, L_biomass, L_nitrate  &
      , L_co2_interactive, L_USE_CARIOLLE, l_flux_bc, l_spec_z0         &

! in coordinate information
      ,z_rho, z_theta                                                   &
! in time stepping information.
      , timestep, timestep_number                                       &
! SCM/Idealised UM 
      , flux_e,flux_h,z0m_scm,z0h_scm                                   &
      , t_surf, zh, u_0_p, v_0_p                                        &
      , ls_rain, ls_snow                                                &
!
! SCM diagnostics (dummy in full UM)
      , nSCMDpkgs,L_SCMDiags                                            &

! in data fields.
      , rho_rsq, rho_wet, rho_wet_theta                                 &
      , u_p,v_p, ustar_p, vstar_p,w_copy, p, p_star, exner_rho_levels   &
      , land_sea_mask, flandg, ice_fract                                &
      , tstar_land, tstar_sea, tstar_sice, z0msea                       & 
      , p_layer_boundaries, p_layer_centres                             &
      , exner_layer_boundaries, exner_layer_centres                     &
      , t1_sd, q1_sd, exner_theta_levels                                &
      , uw0, vw0, w_max, zlcl, zlcl_uv, zhpar, dzh, entrain_coef        &
      , conv_type                                                       &
      , cumulus, l_shallow, l_congestus, l_congestus2,l_pc2_diag_sh_pts &
      , no_cumulus                                                      &
      , ntml, ntpar                                                     &
      , wstar, wthvs, delthvu, ql_ad, qsat_lcl, ftl ,fqt                &
      , shallowc, cu_over_orog, cape_undilute, cin_undilute             &
      , deep_flag, past_precip, past_conv_ht                            &

! in start of time step prognostic values
      , theta, q, qcl, qcf, qrain, qgraup, qcf2                         &
      , cloud_fraction_liquid_halos                                     &
      , cloud_fraction_frozen_halos, bulk_cloud_fraction_halos          &
! in/out values with all increments added
      , theta_conv,q_conv,qcl_conv,qcf_conv, qrain_conv, qgraup_conv    &
      , qcf2_conv, cf_liquid_conv, cf_frozen_conv, bulk_cf_conv         &
! in/out Convective increments
      , theta_inc, q_inc, qcl_inc, qcf_inc, cf_liquid_inc, cf_frozen_inc&
      , bulk_cf_inc                                                     &
      , r_u_p, r_v_p, AEROSOL                                           &
      , DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
      , SO2, SO4_AITKEN, SO4_ACCU, SO4_DISS                             &
      , dms, nh3, soot_new, soot_aged, soot_cld, bmass_new, bmass_aged  &
      , bmass_cld, ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss   &
      , co2,  free_tracers, ukca_tracers                                &
      , OZONE_TRACER                                                    &

! out fields
      , cca0, ccw0, ccb0, cct0, cclwp0, lcbase0, cca0_2d, lctop, lcca   &
      , cca , ccw , ccb , cct , cclwp , lcbase,  cca_2d                 &
      , conv_rain, conv_snow, ddmfx                                     &

! error information
      , Error_code  )

! Arrays from Convection required by SKEB2
! Passed through module stochastic_physics_run_mod
! Deallocation occurs in stph_skeb2
      IF  (.NOT.ALLOCATED(skeb2_up_flux)) THEN
        ALLOCATE (skeb2_up_flux(row_length, rows, model_levels))
      END IF
      IF  (.NOT.ALLOCATED(skeb2_dwn_flux)) THEN
        ALLOCATE (skeb2_dwn_flux(row_length, rows, model_levels))
      END IF
      IF  (.NOT.ALLOCATED(skeb2_cape)) THEN
        ALLOCATE (skeb2_cape(row_length, rows))
      END IF
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            skeb2_up_flux(i,j,k) = up_flux(i,j,k)
            skeb2_dwn_flux(i,j,k) = dwn_flux(i,j,k)
          END DO
        END DO
      END DO
      DO j = 1, rows
        DO i = 1, row_length
          skeb2_cape(i,j) = cape_out(i,j)
        END DO
      END DO

! Copy variables required by COSP
      IF (L_cosp) THEN
        cosp_crain_3d(:,:,:) = 0.0
        cosp_csnow_3d(:,:,:) = 0.0
        cosp_crain_3d(:,:,1:wet_levels) = conv_rain_3d(:,:,:)
        cosp_csnow_3d(:,:,1:wet_levels) = conv_snow_3d(:,:,:)
        WHERE (cosp_crain_3d < 0.0) cosp_crain_3d = 0.0
        WHERE (cosp_csnow_3d < 0.0) cosp_csnow_3d = 0.0
      END IF



      IF (l_mom) THEN
        ALLOCATE(work_dubydt_p(pdims_s%i_start:pdims_s%i_end,      & 
                               pdims_s%j_start:pdims_s%j_end,      &
                               pdims_s%k_start:pdims_s%k_end) )
        ALLOCATE(work_dvbydt_p(pdims_s%i_start:pdims_s%i_end,      & 
                               pdims_s%j_start:pdims_s%j_end,      &
                               pdims_s%k_start:pdims_s%k_end) )
        ! U & V increments
        ! need to copy increments into arrays with one point haloes
        ! ready for swap bounds call (see much later in this routine)

        DO k = 1, tdims%k_end 
          DO j = 1, rows
            DO i = 1, row_length
              work_dubydt_p(i,j,k) = dubydt_pout(i,j,k)
              work_dvbydt_p(i,j,k) = dvbydt_pout(i,j,k)
            END DO
          END DO
        END DO
      END IF 

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Convection',6)

      End If ! on error code equal to zero

      End If ! on convection option

! The following needs to be executed regardless if the conv scheme
! is used.
! Restore thermo vars to their value before convection as they are
! needed in the BL for the calculation of NT terms and increments
! when L_phys2_substep=T

      IF &
      (l_param_conv) &
      THEN

      If ( Error_code  ==  0 ) Then

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Convection',5)

! update star variables by adding on increments

! Start of OpenMP parallel region
!$OMP  PARALLEL DEFAULT(NONE)                                        &
!$OMP& SHARED(model_levels, rows, row_length, wet_levels,            &
!$OMP& theta_star, theta_inc, q_star, q_inc, tdims, qdims)           &
!$OMP& PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
      Do k =             1, tdims%k_end
        Do j = tdims%j_start, tdims%j_end
          Do i = tdims%i_start, tdims%i_end
            theta_star(i,j,k) = theta_star(i,j,k)                       &
     &                        + theta_inc(i,j,k)
          End Do
        End Do
      End Do
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC) 
        Do k =             1, qdims%k_end
          Do j = qdims%j_start, qdims%j_end
            Do i = qdims%i_start, qdims%i_end
            q_star(i,j,k) = q_star(i,j,k) + q_inc(i,j,k)
          End Do
        End Do
      End Do
!$OMP END DO

! End of OpenMP parallel region
!$OMP END PARALLEL

! ----------------------------------------------------------------------
! Protected loop. Update increments only when PC2 scheme is fully ON-ON.
! ----------------------------------------------------------------------
!
! L_calc_dxek_if1:
      If (L_calc_dxek .AND. L_q_interact) Then
!
        Do k =             1, qdims%k_end
          Do j = qdims%j_start, qdims%j_end
            Do i = qdims%i_start, qdims%i_end
              qcl_star(i,j,k) = qcl_star(i,j,k) + qcl_inc(i,j,k)
              qcf_star(i,j,k) = qcf_star(i,j,k) + qcf_inc(i,j,k)
              cfl_star(i,j,k) = cfl_star(i,j,k) + cf_liquid_inc(i,j,k)
              cff_star(i,j,k) = cff_star(i,j,k) + cf_frozen_inc(i,j,k)
              cf_star(i,j,k)  = cf_star(i,j,k)  + bulk_cf_inc(i,j,k)
            End Do
          End Do
        End Do
!
      End If  ! L_calc_dxek_if1
!
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Convection',6)

      Endif ! Error_code

! ----------------------------------------------------------------------
! Other prognostics if present are not updated by convection so no
! copying back to qstar variables required.
! i.e. qrain, qgraup, qcf2
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Section CNV.2 Energy correction code
! ----------------------------------------------------------------------

      If ( CycleNo == NumCycles .and.                                   &
     &     L_emcorr .and. Error_code == 0 )  Then

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Conv Eng Corr',5)

! Add convective rain and snow, at the surface to the
! diabatic heating for use in the energy correction
! procedure.
! Scale variables by conversion factor so that only one call is required

        lclf = lc + lf
        Do j = 1, rows
          Do i = 1, row_length
            tot_precip_scaled(i,j) =  conv_rain(i,j)                    &
     &                               * lc +                             &
     &                                conv_snow(i,j)                    &
     &                               * lclf
          End Do
        End Do

! DEPENDS ON: flux_diag
        CALL flux_diag(tot_precip_scaled, xx_cos_theta_latitude,        &
     &                 row_length, rows ,offx,offy, 1.0,                &
     &                  sum_eng_fluxes,    timestep)

        Do j = 1, rows
          Do i = 1, row_length
            tot_precip_scaled(i,j) = -conv_rain(i,j)-conv_snow(i,j)
          End Do
        End Do

! DEPENDS ON: flux_diag
        CALL flux_diag(tot_precip_scaled, xx_cos_theta_latitude,        &
     &                 row_length, rows ,offx,offy, 1.0,                &
     &                  sum_moist_flux,    timestep)

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Conv Eng Corr',6)

      End If   ! L_emcorr

      End If ! on convection option

! ----------------------------------------------------------------------
! Section BL Call implicit solver
! ---------------------------------------------------------------------
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Boundary Layer',5)
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Implicit BL',5)

! Start of OpenMP parallel region
!$OMP  PARALLEL DEFAULT(NONE)                                           &
!$OMP& SHARED(model_levels, rows, row_length, T_latest, theta_star,     &
!$OMP& exner_theta_levels, wet_levels, q_latest, q_star, qcl_latest,    &
!$OMP& qcl_star, qcf_latest, qcf_star, cf_latest, cf_star, cfl_latest,  &
!$OMP& cfl_star, cff_latest, cff_star, tdims, qdims)                    &
!$OMP& PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
      Do k =             1, tdims%k_end
        Do j = tdims%j_start, tdims%j_end
          Do i = tdims%i_start, tdims%i_end
            T_latest(i,j,k) = theta_star(i,j,k)                         &
     &                      * exner_theta_levels(i,j,k)
          End Do
        End Do
      End Do
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
      Do k =             1, qdims%k_end
        Do j = qdims%j_start, qdims%j_end
          Do i = qdims%i_start, qdims%i_end
            q_latest(i,j,k) = q_star(i,j,k)
            qcl_latest(i,j,k) = qcl_star(i,j,k)
            qcf_latest(i,j,k) = qcf_star(i,j,k)
            cf_latest(i,j,k)  = cf_star(i,j,k)
            cfl_latest(i,j,k) = cfl_star(i,j,k)
            cff_latest(i,j,k) = cff_star(i,j,k)
          End Do
        End Do
      End Do
!$OMP END DO

! End of OpenMP parallel region
!$OMP END PARALLEL
 
      If ( Error_code  ==  0 ) Then

      IF (L_bl) Then

!---------------------------------------------------------------------
! Intercept values of physics increments before being updated by
! implicit solver for optional output of bl increments
!---------------------------------------------------------------------
! [ Note: u/v_incr_diag_bl has been redefined over bl_levels to
!   saved space, since upper level increments should be all zero. ]

! Set diagnostic flags required for boundary layer diagnostics from
! STASHflags.
        L_u_incr_bl = sf(185,3) .and. L_apply_diag
        L_v_incr_bl = sf(186,3) .and. L_apply_diag
        L_T_incr_bl_lsc = sf(181,9) .and. L_apply_diag
        L_Tl_incr_bl_lsc = sf(189,3) .and. L_apply_diag
        L_q_incr_bl_lsc = sf(182,9) .and. L_apply_diag
        L_qtl_incr_bl_lsc = sf(190,3) .and. L_apply_diag
        L_qcl_incr_bl_lsc = sf(183,9) .and. L_apply_diag
        L_qcf_incr_bl_lsc = (sf(184,3).OR.sf(172,3).OR.sf(173,3))       &
                            .and. L_apply_diag
        L_cf_incr_bl  = sf(192,3) .and. L_apply_diag
        L_cfl_incr_bl = (sf(193,3).OR.sf(176,3).OR.sf(177,3))           &
                            .and. L_apply_diag
        L_cff_incr_bl = (sf(194,3).OR.sf(178,3).OR.sf(179,3))           & 
                            .and. L_apply_diag
        L_qcl_incr_bl = (sf(183,3).OR.sf(170,3).OR.sf(171,3))           &
                            .and. L_apply_diag
        L_q_incr_bl   = sf(182,3) .and. L_apply_diag
        L_T_incr_bl   = sf(181,3) .and. L_apply_diag
!
! Allocate and initialize u,v,T,Q,Qcl,Qcf increments.
! Diagnostics for u,v,T,Q,Qcl,Qcf increments need to retain their values
! between calls, to allow their summation over the number of substeps,
! and thus had to be moved out of ni_imp_ctl().
!
        If ( L_u_incr_bl ) Then
          Allocate ( u_incr_diag_bl(udims%i_start:udims%i_end,        &
                                    udims%j_start:udims%j_end,        &
                                    udims%k_start:udims%k_end) )
          u_incr_diag_bl(:,:,:) = 0.0
        Else
          Allocate ( u_incr_diag_bl(1,1,1) )
        Endif                   ! on STASHflag

        If ( L_v_incr_bl ) Then
          Allocate ( v_incr_diag_bl(vdims%i_start:vdims%i_end,        &
                                    vdims%j_start:vdims%j_end,        &
                                    vdims%k_start:vdims%k_end) )
          v_incr_diag_bl(:,:,:) = 0.0
        Else
          Allocate ( v_incr_diag_bl(1,1,1) )
        Endif                   ! on STASHflag

        If ( L_T_incr_bl_lsc .OR. L_Tl_incr_bl_lsc                    &
     &      .or. L_T_incr_bl ) Then
          Allocate ( T_incr_diag_bl(tdims%i_start:tdims%i_end,        &
                                    tdims%j_start:tdims%j_end,        &
                                                1:tdims%k_end) )
          T_incr_diag_bl(:,:,:) = 0.0
        Else
          Allocate ( T_incr_diag_bl(1,1,1) )
        End if                   ! on STASHflags

        If ( L_q_incr_bl_lsc .OR. L_qtl_incr_bl_lsc                   &
     &      .or. L_q_incr_bl ) Then
          Allocate ( q_incr_diag_bl(qdims%i_start:qdims%i_end,        &
                                    qdims%j_start:qdims%j_end,        &
                                                1:qdims%k_end) )
          q_incr_diag_bl(:,:,:) = 0.0
        Else
          Allocate ( q_incr_diag_bl(1,1,1) )
        End if                  ! on STASHflags

        If ( L_qcl_incr_bl_lsc .OR. L_Tl_incr_bl_lsc .OR.             &
     &       L_qtl_incr_bl_lsc .or. L_qcl_incr_bl) Then
          Allocate ( qcl_incr_diag_bl(qdims%i_start:qdims%i_end,      &
                                      qdims%j_start:qdims%j_end,      &
                                                  1:qdims%k_end) )
          qcl_incr_diag_bl(:,:,:) = 0.0
        Else
          Allocate ( qcl_incr_diag_bl(1,1,1) )
        End if                  ! on STASHflag

        If ( L_qcf_incr_bl_lsc ) Then
          Allocate (qcf_incr_diag_bl(qdims%i_start:qdims%i_end,       &
                                     qdims%j_start:qdims%j_end,       &
                                                 1:qdims%k_end))
          qcf_incr_diag_bl(:,:,:) = 0.0
        Else
          Allocate ( qcf_incr_diag_bl(1,1,1) )
        End if                  ! on STASHflag!

        If ( L_cf_incr_bl ) Then
          Allocate (bulk_cf_incr_diag_bl(qdims%i_start:qdims%i_end,   &
                                         qdims%j_start:qdims%j_end,   &
                                                     1:qdims%k_end))
          bulk_cf_incr_diag_bl(:,:,:) = 0.0
        Else
          Allocate (bulk_cf_incr_diag_bl(1,1,1))
        Endif

        If ( L_cfl_incr_bl ) Then
          Allocate (cf_liquid_incr_diag_bl(qdims%i_start:qdims%i_end,  &
                                           qdims%j_start:qdims%j_end,  &
                                                       1:qdims%k_end))
          cf_liquid_incr_diag_bl(:,:,:) = 0.0
        Else
          Allocate (cf_liquid_incr_diag_bl(1,1,1))
        Endif

        If ( L_cff_incr_bl ) Then
          Allocate (cf_frozen_incr_diag_bl(qdims%i_start:qdims%i_end,  &
                                           qdims%j_start:qdims%j_end,  &
                                                       1:qdims%k_end))
         cf_frozen_incr_diag_bl(:,:,:) = 0.0
        Else
          Allocate (cf_frozen_incr_diag_bl(1,1,1))
        Endif

      End If ! L_bl
!--------------------------------------------------------------------------
! Second communication section before ni_imp_ctl
!--------------------------------------------------------------------------


! swap bounds for 3d fields first
    i_field = 0

  IF (L_bl) THEN
! only call on 1st cycle or if not fast running
     IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => rhokm(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  bl_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.

        IF(formdrag ==  explicit_stress)THEN

           i_field = i_field + 1
           fields_to_swap(i_field) % field       => tau_fd_x(:,:,:)
           fields_to_swap(i_field) % field_type  =  fld_type_p
           fields_to_swap(i_field) % levels      =  bl_levels
           fields_to_swap(i_field) % rows        =  rows
           fields_to_swap(i_field) % vector      =  .FALSE.

           i_field = i_field + 1
           fields_to_swap(i_field) % field       => tau_fd_y(:,:,:)
           fields_to_swap(i_field) % field_type  =  fld_type_p
           fields_to_swap(i_field) % levels      =  bl_levels
           fields_to_swap(i_field) % rows        =  rows
           fields_to_swap(i_field) % vector      =  .FALSE.

        END IF

        i_field = i_field + 1
        fields_to_swap(i_field) % field       => f_ngstress(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  bl_levels-1
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.

     END IF !l_quick_ap2

  END IF !L_bl

! Convection winds fields requiring swap bound if CMT being used

  IF (l_param_conv .AND. l_mom) THEN

    i_field = i_field + 1
    fields_to_swap(i_field) % field       => work_dubydt_p(:,:,:)
    fields_to_swap(i_field) % field_type  =  fld_type_p
    fields_to_swap(i_field) % levels      =  model_levels
    fields_to_swap(i_field) % rows        =  rows
    fields_to_swap(i_field) % vector      =  .FALSE.

    i_field = i_field + 1
    fields_to_swap(i_field) % field       => work_dvbydt_p(:,:,:)
    fields_to_swap(i_field) % field_type  =  fld_type_p
    fields_to_swap(i_field) % levels      =  model_levels
    fields_to_swap(i_field) % rows        =  rows
    fields_to_swap(i_field) % vector      =  .FALSE.

  END IF


! DEPENDS ON: swap_bounds_mv
  CALL swap_bounds_mv(fields_to_swap, i_field, row_length,              &
                      offx , offy )

  IF (L_bl) Then
! only call on 1st cycle or if not fast running
     IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN


! then 2d fields
        i_field = 0
        i_field = i_field + 1
        fields_to_swap(i_field) % field_2d   => rhokm_land(:,:)
        fields_to_swap(i_field) % field_type = fld_type_p
        fields_to_swap(i_field) % levels     = 1
        fields_to_swap(i_field) % rows       = rows
        fields_to_swap(i_field) % vector     = .FALSE.

        i_field = i_field + 1
        fields_to_swap(i_field) % field_2d   => rhokm_ssi(:,:)
        fields_to_swap(i_field) % field_type = fld_type_p
        fields_to_swap(i_field) % levels     = 1
        fields_to_swap(i_field) % rows       = rows
        fields_to_swap(i_field) % vector     = .FALSE.

        IF (l_ctile .AND. buddy_sea == on) THEN
        ! Interpolate wind speed factors to u and v columns

           i_field = i_field + 1
           fields_to_swap(i_field) % field_2d   => flandfac(:,:)
           fields_to_swap(i_field) % field_type = fld_type_p
           fields_to_swap(i_field) % levels     = 1
           fields_to_swap(i_field) % rows       = rows
           fields_to_swap(i_field) % vector     = .FALSE.

           i_field = i_field + 1
           fields_to_swap(i_field) % field_2d   => fseafac(:,:)
           fields_to_swap(i_field) % field_type = fld_type_p
           fields_to_swap(i_field) % levels     = 1
           fields_to_swap(i_field) % rows       = rows
           fields_to_swap(i_field) % vector     = .FALSE.

        END IF

        su10 = (sf(209,3) .OR. sf(225,3) .OR. sf(227,3) .OR.            &
                               sf(230,3) .OR. sf(463,3))                &
                .AND.  (L_apply_diag .OR. l_quick_ap2)
        sv10 = (sf(210,3) .OR. sf(226,3) .OR. sf(227,3) .OR.            &
                               sf(230,3) .OR. sf(463,3))                &
                .AND.  (L_apply_diag .OR. l_quick_ap2)

        IF (su10 .OR. sv10) THEN
           i_field = i_field + 1
           fields_to_swap(i_field) % field_2d   => cdr10m(:,:)
           fields_to_swap(i_field) % field_type = fld_type_p
           fields_to_swap(i_field) % levels     = 1
           fields_to_swap(i_field) % rows       = rows
           fields_to_swap(i_field) % vector     = .FALSE.
        END IF

! DEPENDS ON: swap_bounds_2d_mv
        CALL swap_bounds_2d_mv(fields_to_swap, i_field, row_length,     &
                               offx , offy )

        CALL p_to_u(rhokm,                                          &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     udims%i_start,udims%i_end,                     &
                     udims%j_start,udims%j_end,                     &
                     1,bl_levels, rhokm_u)
              
        CALL p_to_u_land(rhokm_land, flandg,                        &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     udims%i_start,udims%i_end,                     &
                     udims%j_start,udims%j_end,                     &
                     1,1, rhokm_u_land)    

        CALL p_to_u_sea(rhokm_ssi, flandg,                          &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     udims%i_start,udims%i_end,                     &
                     udims%j_start,udims%j_end,                     &
                     1,1, rhokm_u_ssi)    

        CALL p_to_u(flandg,                                         &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     udims%i_start,udims%i_end,                     &
                     udims%j_start,udims%j_end,                     &
                     1,1,flandg_u)

        CALL p_to_v(rhokm,                                          &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     vdims%i_start,vdims%i_end,                     &
                     vdims%j_start,vdims%j_end,                     &
                     1,bl_levels, rhokm_v)
                     
        CALL p_to_v_land(rhokm_land, flandg,                        &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     vdims%i_start,vdims%i_end,                     &
                     vdims%j_start,vdims%j_end,                     &
                     1,1, rhokm_v_land)    

        CALL p_to_v_sea(rhokm_ssi, flandg,                          &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     vdims%i_start,vdims%i_end,                     &
                     vdims%j_start,vdims%j_end,                     &
                     1,1, rhokm_v_ssi)    

        CALL p_to_v(flandg,                                         &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     vdims%i_start,vdims%i_end,                     &
                     vdims%j_start,vdims%j_end,                     &
                     1,1,flandg_v)

        IF (l_ctile .AND. buddy_sea == on) THEN

           CALL p_to_u(flandfac,                                    &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     udims%i_start,udims%i_end,                     &
                     udims%j_start,udims%j_end,                     &
                     1,1, flandfac_u)

           CALL p_to_u(fseafac,                                     &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     udims%i_start,udims%i_end,                     &
                     udims%j_start,udims%j_end,                     &
                     1,1, fseafac_u)

           CALL p_to_v(flandfac,                                    &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     vdims%i_start,vdims%i_end,                     &
                     vdims%j_start,vdims%j_end,                     &
                     1,1, flandfac_v)

           CALL p_to_v(fseafac,                                     &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     vdims%i_start,vdims%i_end,                     &
                     vdims%j_start,vdims%j_end,                     &
                     1,1, fseafac_v)

        END IF !c_tile and buddy_sea


        IF (su10)THEN
           CALL p_to_u(cdr10m,                                      &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     udims%i_start,udims%i_end,                     &
                     udims%j_start,udims%j_end,                     &
                     1,1, cdr10m_u)

        END IF

        IF (sv10)THEN
           CALL p_to_v(cdr10m,                                      &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     vdims%i_start,vdims%i_end,                     &
                     vdims%j_start,vdims%j_end,                     &
                     1,1, cdr10m_v)

        END IF


ALLOCATE(rhogamu_u(1,1,1))
ALLOCATE(rhogamv_v(1,1,1))
ALLOCATE(f_ngstress_u(udims%i_start:udims%i_end,                    &
                      udims%j_start:udims%j_end,2:bl_levels))
ALLOCATE(f_ngstress_v(vdims%i_start:vdims%i_end,                    &
                      vdims%j_start:vdims%j_end,2:bl_levels))


        IF(formdrag ==  explicit_stress)THEN

           CALL p_to_u (tau_fd_x,                                      &
                        pdims_s%i_start,pdims_s%i_end,                 &
                        pdims_s%j_start,pdims_s%j_end,                 &
                        udims%i_start,udims%i_end,                     &
                        udims%j_start,udims%j_end,                     &
                        1,bl_levels,taux_fd_u)

           CALL p_to_v (tau_fd_y,                                      &
                        pdims_s%i_start,pdims_s%i_end,                 &
                        pdims_s%j_start,pdims_s%j_end,                 &
                        vdims%i_start,vdims%i_end,                     &
                        vdims%j_start,vdims%j_end,                     &
                        1,bl_levels,tauy_fd_v)
                   
        END IF

        CALL p_to_u (f_ngstress,                                    &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     udims%i_start,udims%i_end,                     &
                     udims%j_start,udims%j_end,                     &
                     2,bl_levels,f_ngstress_u)

        CALL p_to_v (f_ngstress,                                    &
                     pdims_s%i_start,pdims_s%i_end,                 &
                     pdims_s%j_start,pdims_s%j_end,                 &
                     vdims%i_start,vdims%i_end,                     &
                     vdims%j_start,vdims%j_end,                     &
                     2,bl_levels,f_ngstress_v)   



! deallocate temp variables for message passing
  DEALLOCATE(rhogamu)
  DEALLOCATE(rhogamv)
  DEALLOCATE(f_ngstress)

     END IF !l_quick_ap2
  END IF !L_bl
!----------------------------------------------------------------------------
! Rest of convection section regridding fields after swap_bounds and 
! outputing convection diagnostics to stash
!----------------------------------------------------------------------------

  IF (l_param_conv) THEN
    IF (Ltimer) Call timer ('AP2 Convection',5)

    IF (l_mom) THEN

      ALLOCATE(dubydt_u(udims%i_start:udims%i_end,             & 
                        udims%j_start:udims%j_end,             &
                        udims%k_start:udims%k_end) )
      ALLOCATE(dvbydt_v(vdims%i_start:vdims%i_end,             & 
                        vdims%j_start:vdims%j_end,             &
                        vdims%k_start:vdims%k_end) )

      ! interpolate to u grid
      ! CONV have set work_dubydt_p to be on tdims grid.
      CALL p_to_u (work_dubydt_p,                                 &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   1,udims%k_end, dubydt_u)                                    
           
      ! add on to increment field
      IF ( l_u_incr_conv ) THEN
      ! update R_u and diagnostics
        DO k=1,udims%k_end
          DO j=udims%j_start,udims%j_end
            DO i=udims%i_start,udims%i_end
              r_u(i,j,k) = r_u(i,j,k)  + dubydt_u(i,j,k)*timestep

              u_incr_diag_conv(i,j,k) = u_incr_diag_conv(i,j,k)      &
                                       + dubydt_u(i,j,k)* timestep
            END DO ! i
          END DO ! j
        END DO ! k
      ELSE
        DO k = 1, udims%k_end
          DO j=udims%j_start,udims%j_end
            DO i=udims%i_start,udims%i_end
              r_u(i,j,k) = r_u(i,j,k) + dubydt_u(i,j,k)*timestep
            END DO
          END DO
        END DO
      END IF                   ! on STASHflag

      ! V increments
      ! first need to copy increments into arrays with one point haloes
      ! and swop boundaries

      ! interpolate to v grid
      ! CONV have set work_dubydt_p to be on tdims grid.
      CALL p_to_v (work_dvbydt_p,                                 &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   vdims%i_start,vdims%i_end,                     &
                   vdims%j_start,vdims%j_end,                     &
                   1,vdims%k_end, dvbydt_v)  
                   
      ! add on to increment field
      IF ( l_v_incr_conv ) THEN
        DO k=1, vdims%k_end
          DO j=vdims%j_start,vdims%j_end
            DO i=vdims%i_start,vdims%i_end
              r_v(i,j,k) = r_v(i,j,k)                            &
                         + dvbydt_v(i,j,k)*timestep
              v_incr_diag_conv(i,j,k) = v_incr_diag_conv(i,j,k)  &
                                      + dvbydt_v(i,j,k) * timestep
            END DO ! i
          END DO ! j
        END DO ! k
      ELSE
        DO k = 1, vdims%k_end
          DO j=vdims%j_start,vdims%j_end
            DO i=vdims%i_start,vdims%i_end
              r_v(i,j,k)=r_v(i,j,k)+dvbydt_v(i,j,k)*timestep
            END DO
          END DO
        END DO
      END IF                   ! on STASHflag

      DEALLOCATE(dvbydt_v)
      DEALLOCATE(dubydt_u)
      DEALLOCATE(work_dvbydt_p)
      DEALLOCATE(work_dubydt_p)
    END IF    !l_mom

    ! Output section 5 (convection) stash diagnostics. Moved to atmos_physics2
    ! so that swap_bounds calls for wind increments can be combined with 
    ! BL swap_bounds.
    ! Note SCM calls for output still in ni_conv_ctl


    IF ( sf(0,5) ) THEN

      IF ( l_apply_diag ) THEN
 
! DEPENDS ON: diagnostics_conv
          CALL diagnostics_conv(                                  &
                       row_length, rows, tdims%k_end, qdims%k_end &
,                      at_extremity                               &
,                      u, v, p, r_u, r_v                          &
,                      theta_inc, q_inc                           &
,                      qcl_inc, qcf_inc, cf_liquid_inc            &
,                      cf_frozen_inc, bulk_cf_inc                 &
,                      exner_theta_levels                         &
,                      ls_rain, ls_snow                           &
,                      ccw, conv_rain, conv_snow                  &
,                      cca_2d, cca, ccb, cct                      &
,                      cu_over_orog, cape_undilute, cin_undilute  &
,                      lcbase, lctop, lcca, n_cca_levels          &
,                      l_dust                                     &
,                      conv_type                                  &
,                      timestep                                   &
,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
   STASHwork5                                                     &
   )

      END IF       ! l_apply_diag
    END IF         ! test on sf(0,5)
   

    ! Deallocate convective diagnostics arrays allocated earlier

! DEPENDS ON: cv_dealloc_diag_array
      CALL cv_dealloc_diag_array( )

    IF (Ltimer) Call timer ('AP2 Convection',6)

  END IF ! l_param_conv  test

!----------------------------------------------------------------------
! End of second communication section
!----------------------------------------------------------------------

  IF (L_bl) THEN
! only call on 1st cycle or if not fast running
     IF ( cycleno == 1 .OR. .NOT. l_quick_ap2 ) THEN

! DEPENDS ON: bdy_expl3
        CALL bdy_expl3 (                                                &
! IN grid related variables 
             bl_levels,                                                 &
! IN SCM diags
             nSCMDpkgs,L_SCMDiags,                                      &
! IN variables used in flux calculations
             u, v, u_0, v_0, rhokm_u_land, rhokm_v_land, flandfac_u,    &
             flandfac_v, rhokm_u_ssi, rhokm_v_ssi, fseafac_u, fseafac_v,&
             flandg_u, flandg_v, zh, rdz_u, rdz_v, rhokm_u, rhokm_v,    &
             taux_fd_u, tauy_fd_v,                                      &
             rhogamu_u, rhogamv_v, f_ngstress_u, f_ngstress_v,          &
! OUT explicit momentum fluxes
             taux_land, tauy_land, taux_ssi, tauy_ssi, taux, tauy       &
             )

        IF (l_quick_ap2) THEN
! save outputs for second EG cycle
! N.B. any new OUT data added to bdy_expl3 will need to be saved here
! and restored below. The array size will also need to be increased in
! atmos_physics2_alloc
           bdy_expl3_u3d(:,:,:)=taux
           bdy_expl3_v3d(:,:,:)=tauy
           bdy_expl3_u2d(:,:,1)=taux_land
           bdy_expl3_v2d(:,:,1)=tauy_land
           bdy_expl3_u2d(:,:,2)=taux_ssi
           bdy_expl3_v2d(:,:,2)=tauy_ssi
        END IF

DEALLOCATE(rhogamu_u)
DEALLOCATE(rhogamv_v)
DEALLOCATE(f_ngstress_u)
DEALLOCATE(f_ngstress_v)

     ELSE
! restore outputs on second EG cycle
        taux=bdy_expl3_u3d(:,:,:)
        tauy=bdy_expl3_v3d(:,:,:)
        taux_land=bdy_expl3_u2d(:,:,1)
        tauy_land=bdy_expl3_v2d(:,:,1)
        taux_ssi=bdy_expl3_u2d(:,:,2)
        tauy_ssi=bdy_expl3_v2d(:,:,2)

     END IF !l_quick_ap2
!-----------------------------------------------------------------------
! Reset index for sea-ice categories
!-----------------------------------------------------------------------
! These indices are set up initially in atmos_physics1 using nice_use as
! there is the option of only partially using the categories in the 
! radation and explicit surface scheme.
! Here they are reset using nice as all the categories are always used in 
! the implicit surface scheme.

      sice_pts_ncat(:)=0
      sice_index_ncat(:,:)=0
      sice_frac_ncat(:,:)=0.0
      DO n=1,nice
        DO l=1,ssi_pts
          j=(ssi_index(l)-1)/row_length + 1
          i = ssi_index(l) - (j-1)*row_length
          IF (ssi_index(l) > 0) THEN
            IF (ice_fract_ncat(i,j,n) > 0.0) THEN
              sice_pts_ncat(n)=sice_pts_ncat(n)+1
              sice_index_ncat(sice_pts_ncat(n),n)=l
              sice_frac_ncat(l,n)=ice_fract_ncat(i,j,n)
            END IF
          END IF
        END DO
      END DO

! DEPENDS ON: ni_imp_ctl
      Call NI_imp_ctl (                                                 &

! Parallel variables
     &  global_row_length, global_rows, gc_proc_row_group               &
     &, gc_proc_col_group, at_extremity, n_proc, n_procx                &
     &, n_procy, neighbour, g_rows, g_row_length, mype                  &

! model dimensions.
     &, rhc_row_length, rhc_rows, land_points                           &
     &, ntiles, bl_levels, dst_levels                                   &
     &, dsm_levels, cloud_levels, n_cca_levels, nice, nice_use          &
     &, DIM_CS1, DIM_CS2                                                &
!
! Model switches
     &, model_domain, L_CAL360, L_area_cloud, L_ACF_Cusack              &
     &, L_ACF_Brooks, L_RHCPT, L_emcorr, Ltimer, L_DRY,L_MURK           &
     &, L_MURK_ADVECT,L_BL_TRACER_MIX                                   &
     &, L_DUST,L_DUST_DIAG,L_SULPC_SO2                                  &
     &, L_sulpc_nh3, L_sulpc_dms, L_soot, L_biomass, L_ocff, L_nitrate  &
     &, L_co2_interactive                                               &
     &, L_co2_emits, L_pc2                                              &
     &, NumCycles, CycleNo, L_mixing_ratio                              &
     &, L_ukca, L_sice_heatflux                                         &
     &, L_sice_multilayers                                              &     
     &, L_USE_CARIOLLE                                                  &

! model Parameters
     &, rhcrit, tr_vars, tr_ukca, co2_mmr                               &

! in coordinate information
     &, delta_lambda, delta_phi,lat_rot_NP,long_rot_NP                  &

! in time stepping information.
     &, timestep, val_year, val_day_number, val_hour, val_minute        &
      , val_second, timestep_number                                     &

! trig arrays
      , sin_theta_longitude, cos_theta_longitude, xx_cos_theta_latitude &
      , f3_at_u,                                                        &

! diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork3, STASHwork9                                           &
!
! SCM Diagnostics (dummy in full UM)
     &, nSCMDpkgs, L_SCMDiags                                           &

! in data fields.
     &, p, p_layer_centres, p_layer_boundaries, rho_rsq, rho_wet_theta  &
     &, u, v, w                                                         &
     &, land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels   &
     &, theta_conv,q_conv,qcl_conv,qcf_conv                             &
! ancillary fields and fields needed to be kept from timestep to
! timestep
     &, smvccl, smvcwt, smvcst, sthf, sthu,sil_orog_land, ho2r2_orog    &
     &, di, ice_fract, di_ncat, ice_fract_ncat, k_sice                  &
     &, u_0, v_0, land_index                                            &
     &, cca,  ccb,  cct,  ccw,  cca_2d, lcbase                          &
      , cca0, ccb0, cct0, ccw0, cca0_2d                                 &
     &, ls_rain, ls_snow, conv_rain, conv_snow, L_scrn, L_plsp          &

! in variables required from BDY_LAYR
     &, alpha1_sice, ashtf, BQ_GB, BT_GB, dtrdz_charney_grid            &
     &, rdz_charney_grid, dtrdz_u, dtrdz_v, rdz_u, rdz_v                &
     &, cdr10m_u, cdr10m_v, z1_tq                                       &
     &, uStarGBM                                                        &
!ajm   extra variable
     &, rhokm_u, rhokm_v                                                &

! in diagnostics (started or from) BDY_LAYR
! pass diagn arrays and reduce number of argument lines
      , T_incr_diag_bl, q_incr_diag_bl, qcl_incr_diag_bl                &
      , qcf_incr_diag_bl, u_incr_diag_bl, v_incr_diag_bl                &
      , cf_liquid_incr_diag_bl, cf_frozen_incr_diag_bl                  &
      , bulk_cf_incr_diag_bl                                            &
      , E_SEA,FQT,FTL,H_SEA,RIB_GB,TAUX,TAUY,VSHR,ZLCL,ZHT,dzh          &
      , bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6     &
      , bl_type_7, z0m_gb, z0m_eff_gb, z0h_eff_gb, fme, rhokh           &
      , TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans                          &

! in data required to calculate increment diagnostics
     &, theta_star,q_star                                               &

! IN logical for scm surface forcing
     &, L_flux_bc                                                       &

! in data required for tracer mixing :
     &, RHO_ARESIST,ARESIST,RESIST_B,R_B_DUST                           &
     &, KENT, WE_LIM, T_FRAC, ZRZI                                      &
     &, KENT_DSC, WE_LIM_DSC, T_FRAC_DSC, ZRZI_DSC                      &
     &, ZHSC,Z_RHO,DUST_FLUX,DUST_EMISS_FRAC                            &
     &, U_S_T_TILE,U_S_T_DRY_TILE,U_S_STD_TILE                          &
     &, so2_hilem, so2_em, nh3_em, dms_em, soot_hilem, soot_em          &
     &, ocff_hilem, ocff_em                                             &

! IN additional variables for MOSES II. Now includes lai_ft, canht_ft.
      , TILE_PTS,TILE_INDEX,TILE_FRAC,CANOPY                            &
     &, ALPHA1,FRACA,RHOKH_TILE,SMC,CHR1P5M,RESFS,Z0HSSI,Z0MSSI         &
     &, CANHC_TILE,FLAKE,WT_EXT_TILE,LW_DOWN,lai_ft,canht_ft            &
     &, SW_TILE,ASHTF_TILE,gc,aresist_tile,resist_b_tile                &
     &, FQT_ICE,FTL_ICE,RESFT,RHOKH_SICE,RHOKPM,RHOKPM_POT              &
     &, RHOKPM_SICE,Z0H_TILE,Z0M_TILE,CHR1P5M_SICE                      &
     &, FLAND, FLANDG                                                   &
     &, FLANDG_U,FLANDG_V,TSTAR_SEA,VSHR_LAND,VSHR_SSI                  &

! IN additional variables for JULES
     &, RHOKH_MIX,DTSTAR_TILE,DTSTAR,HCONS,EMIS_TILE,EMIS_SOIL          &

! IN MOSES II variables for STASH
     &, GS,GPP,NPP,RESP_P,GPP_FT,NPP_FT,RESP_P_FT,RESP_S                &
     &, RESP_S_TOT,CS                                                   &
     &, RIB_TILE,FSMC,CATCH,G_LEAF                                      &
     &, CO2_EMITS, CO2FLUX                                              &

!IN additional variables for soil moisture nudging scheme
     &, WT_EXT,RA                                                       &

! in/out
! (Note ti and ti_gb are IN only if l_sice_multilayers=T)
     &, t_soil, ti, t_surf, ti_gb                                       &
     &, area_cloud_fraction, bulk_cloud_fraction                        &
     &, T_latest, q_latest, qcl_latest, qcf_latest                      &
     &, cf_latest, cfl_latest, cff_latest                               &
     &, R_u, R_v, R_w, cloud_fraction_liquid, cloud_fraction_frozen     &
     &, zh, sum_eng_fluxes,sum_moist_flux, rhcpt                        &

! In/Out tracer fields
     &, aerosol, free_tracers                                           &
     &, DUST_DIV1,DUST_DIV2,DUST_DIV3,DUST_DIV4,DUST_DIV5,DUST_DIV6     &
     &, drydep2, so2, dms, so4_aitken, so4_accu, so4_diss, nh3          &
     &, soot_new, soot_aged, soot_cld, bmass_new, bmass_aged            &
     &, bmass_cld, ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss   &
     &, co2, OZONE_TRACER                                               &

! in/out fields
     &, ecan, ei, ext, snowmelt                                         &
     &, t1_sd, q1_sd, ntml, cumulus                                     &
     &, l_pc2_diag_sh_pts, nbdsc, ntdsc                                 &
     &, surf_ht_flux_land, snomlt_surf_htf, cH_term                     &

! INOUT additional variables for MOSES II
     &, T_SURF_TILE,FQT_TILE,EPOT_TILE,FTL_TILE                         &
     &, SNOW_TILE,LE_TILE,RADNET_SICE,RADNET_TILE,OLR                   &
     &, TSTAR_LAND,TSTAR_SICE,TSTAR_SICE_CAT,TSTAR_SSI                  &

! Additional diagnostics for MOSES II
     &, rib_ssi,taux_land,taux_ssi,tauy_land,tauy_ssi                   &

! OUT additional variables for MOSES II
     &, ESOIL_TILE,ES,EI_TILE                                           &
     &, Q1P5M_TILE,T1P5M_TILE,ECAN_TILE,MELT_TILE                       &
     &, SURF_HTF_TILE                                                   &

! error information
     &, Error_code,BL_diag  )
!
      Else

          Do k = 1, wet_levels
            RHCPT(1,1,k) = rhcrit(k)
          End Do
      End If

      End If ! error_code == 0

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Implicit BL',6)
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Boundary Layer',6)

! ----------------------------------------------------------------------
! Section RESTORE: Copy _latest fields into _star locations
! Area cloud fraction has been done in NI_IMP_CTL
! ----------------------------------------------------------------------
      If (error_code  ==  0) Then

! Start of OpenMP parallel region
!$OMP  PARALLEL DEFAULT(NONE)                                       &
!$OMP& SHARED(model_levels, rows, row_length, theta_star, T_latest, &
!$OMP& exner_theta_levels, wet_levels, q_star, q_latest, qcl_star,  &
!$OMP& qcl_latest, qcf_star, qcf_latest, cf_star, cf_latest,        &
!$OMP& cfl_star, cfl_latest, cff_star, cff_latest,                  &
!$OMP& bulk_cloud_fraction_halos, bulk_cloud_fraction,              &
!$OMP& cloud_fraction_liquid_halos, cloud_fraction_liquid,          &
!$OMP& cloud_fraction_frozen, cloud_fraction_frozen_halos,          &
!$OMP& l_pc2, tdims, qdims)                                         &
!$OMP& PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
        Do k =             1, tdims%k_end
          Do j = tdims%j_start, tdims%j_end
            Do i = tdims%i_start, tdims%i_end
              theta_star(i,j,k) = T_latest(i,j,k) /                     &
     &                            exner_theta_levels(i,j,k)
            End Do
          End Do
        End Do
!$OMP END DO nowait
!$OMP DO SCHEDULE(STATIC)
        Do k =             1, qdims%k_end
          Do j = qdims%j_start, qdims%j_end
            Do i = qdims%i_start, qdims%i_end
              q_star(i,j,k) = q_latest(i,j,k)
              qcl_star(i,j,k) = qcl_latest(i,j,k)
              qcf_star(i,j,k) = qcf_latest(i,j,k)
              cf_star(i,j,k)  = cf_latest(i,j,k)
              cfl_star(i,j,k) = cfl_latest(i,j,k)
              cff_star(i,j,k) = cff_latest(i,j,k)
            End Do
          End Do
        End Do
!$OMP END DO nowait
! Reset halo variables for writing back into D1 array on exit from
! this subroutine.
        If (.not. L_pc2 ) Then
!$OMP DO SCHEDULE(STATIC)
          Do k =             1, qdims%k_end
            Do j = qdims%j_start, qdims%j_end
              Do i = qdims%i_start, qdims%i_end
                  bulk_cloud_fraction_halos(i,j,k) =                    &
     &              bulk_cloud_fraction(i,j,k)
                  cloud_fraction_liquid_halos(i,j,k) =                  &
     &              cloud_fraction_liquid(i,j,k)
                  cloud_fraction_frozen_halos(i,j,k) =                  &
     &              cloud_fraction_frozen(i,j,k)
              End Do
            End Do
          End Do
!$OMP END DO
        End If  ! .not. L_pc2

! End of OpenMP parallel region
!$OMP END PARALLEL

      End If ! on error code equal to zero

      If ( L_bl ) Then
          Deallocate ( u_incr_diag_bl )
          Deallocate ( v_incr_diag_bl )
          Deallocate ( T_incr_diag_bl )
          Deallocate ( q_incr_diag_bl )
          Deallocate ( qcl_incr_diag_bl )
          Deallocate ( qcf_incr_diag_bl )
          Deallocate ( bulk_cf_incr_diag_bl )
          Deallocate ( cf_liquid_incr_diag_bl )
          Deallocate ( cf_frozen_incr_diag_bl )
      End If

! ----------------------------------------------------------------------
! Section HYD.1 Compress fields to land points, then call hydrology.
! ----------------------------------------------------------------------

! Call hydrology only at last cycle
      If ( CycleNo == NumCycles) Then
! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Hydrology',5)

      IF (error_code  ==  0 .AND. l_hydrology                           &
          .AND. land_points  /=  0) THEN

! Inland basin outflow is added to soil moisture at each timestep. 
! This flux changes its value only when the river routing scheme 
! has been called in the previous timestep 
        IF (l_rivers) THEN
          IF (l_riverinland_fix) THEN
!         Corrected code. 
!           Pass inland flow to soil moisture every timestep
            IF  (.NOT. l_inland) THEN
              inlandout_atmos=0.0
              inlandout_riv=0.0
            END IF
          ELSE
!         Bugged code.
!           This only passes inflow flux on TRIP timesteps.
!           On other timesteps, inland flow is zero.
            IF(l_inland)THEN 
              IF (a_steps_since_riv /= 0) THEN 
                ! If river routing called in previous timestep 
                inlandout_atmos=0.0 
                inlandout_atm=0.0 
                inlandout_riv=0.0   
              END IF
            ELSE 
! If inland river routing basin not required 
              inlandout_atmos=0.0
              inlandout_riv=0.0
            END IF
          END IF
        END IF  ! l_rivers

! Compress fields to land points

        Do l = 1, land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          ls_rain_land(l)=ls_rain(i,j)
          conv_rain_land(l)=conv_rain(i,j)
          conv_snow_land(l)=conv_snow(i,j)
          ls_snow_land(l)=ls_snow(i,j)
          surf_ht_flux_ld(l) = surf_ht_flux_land(i,j)
          snow_melt(l) = snowmelt(i,j)
        End Do

! call the JULES2.1 snow scheme
! DEPENDS ON: snow_intctl        
        CALL SNOW_INTCTL ( LAND_POINTS,TIMESTEP,SMLT,NTILES,TILE_PTS,   &
                           TILE_INDEX,CATCH_SNOW,CONV_SNOW_LAND,        &
                           TILE_FRAC,LS_SNOW_LAND,EI_TILE,              &
                           HCAP_LEVS(:,1),HCON_LEVS(:,0),               &
                           MELT_TILE,SOIL_LAYER_MOISTURE(:,1),          &
                           STHF(:,1),SURF_HTF_TILE,                     &
                           T_SOIL(:,1),T_SURF_TILE,SMVCST_LEVS(:,1),    &
                           RGRAIN,SNOW_GRND,SNOW_MELT,SNOW_TILE,        &
                           SNOMLT_SURF_HTF,LYING_SNOW,                  &
                           SNOMLT_SUB_HTF,SURF_HT_FLUX_LD )


!-----------------------------------------------------------------------
!     call to the FLake interface
!-----------------------------------------------------------------------
        IF (             l_flake_model                                  &
             .AND. (.NOT.l_aggregate   )                                &
             .AND. (land_points > 0    ) ) THEN

          DO k=1,tile_pts(lake)
            l = tile_index(k,lake)
            j=(land_index(l)-1)/row_length + 1
            i = land_index(l) - (j-1)*row_length

! U* of lake obtained by matching surface stress
            u_s_lake(l) =   u_s_std_tile(l, lake)                       &
                          * SQRT( rho_wet_theta(i,j,1) / rho_water )

! downwelling SW on lake tile
            sw_down(l) = sw_tile(l,lake) / (1. - lake_albedo(l))

! Take the net SW flux out of the surface heat flux
! since this is done separately within FLake.
            surf_ht_flux_lk(l) =   surf_ht_flux_lake(i,j)               &
                                 - sw_tile(l,lake)

          ENDDO

          trap_frozen   = 0
          trap_unfrozen = 0

! DEPENDS ON: flake_interface
      CALL flake_interface( land_points                                 &
                           ,tile_pts(lake)                              &
                           ,tile_index(:,lake)                          &
                           ,u_s_lake                                    &
                           ,surf_ht_flux_lk                             &
                           ,sw_down                                     &
                           ,lake_depth                                  &
                           ,lake_fetch                                  &
                           ,coriolis_param                              &
                           ,timestep                                    &
                           ,lake_albedo                                 &
                           ,lake_t_snow                                 &
                           ,lake_t_ice                                  &
                           ,lake_t_mean                                 &
                           ,lake_t_mxl                                  &
                           ,lake_shape_factor                           &
                           ,lake_h_snow                                 &
                           ,lake_h_ice                                  &
                           ,lake_h_mxl                                  &
                           ,lake_t_sfc                                  &
                           ,ts1_lake                                    &
                           ,g_dt                                        &
                           ,trap_frozen                                 &
                           ,trap_unfrozen )

          IF ( printstatus > PrStatus_Oper ) THEN 
            IF ( trap_frozen > 0 ) THEN
              WRITE(*,*)'AP2-FLake: # zero-divide (  frozen) avoided =',&
                        trap_frozen
            END IF
            IF ( trap_unfrozen > 0 ) THEN
              WRITE(*,*)'AP2-FLake: # zero-divide (unfrozen) avoided =',&
                        trap_unfrozen
            END IF
          END IF

! copy FROM FLake module variables TO Prognostics
            lake_depth_p  = lake_depth
            lake_fetch_p  = lake_fetch
            lake_t_mean_p = lake_t_mean
            lake_t_mxl_p  = lake_t_mxl
            lake_t_ice_p  = lake_t_ice
            lake_h_mxl_p  = lake_h_mxl
            lake_h_ice_p  = lake_h_ice
            lake_shape_p  = lake_shape_factor
            lake_g_dt_p   = g_dt

        ENDIF ! Flake

! DEPENDS ON: hyd_intctl
        Call HYD_INTCTL (                                               &
                     land_ice_points, land_ice_index,                   &
                     soil_points, soil_index,                           &
                     land_points, dsm_levels, clapp_levs,               &
                     conv_rain_land, conv_snow_land,                    &
                     ext, hcap_levs, hcon_levs,                         &
                     ls_rain_land, ls_snow_land,                        &
                     satcon_levs, sathh_levs,                           &
                     surf_ht_flux_ld, timestep,                         &
                     smvcst_levs, smvcwt_levs,                          &
                     canopy_gb, snomlt_surf_htf, smlt,                  &
                     soil_layer_moisture,                               &
                     sthf, sthu, t_soil, lying_snow,                    &
! output
! Beware snow_melt is actually in/out
                     infil, smc, snow_melt, snomlt_sub_htf,             &
                     stf_sub_surf_roff, sub_surf_roff, surf_roff,       &
                     tot_tfall,                                         &
! add inland basin outflow in call to hydctl
                     inlandout_atm,L_INLAND,                            &


! Additional variables for MOSES II
                     NTILES,TILE_PTS,TILE_INDEX,                        &
                     L_SNOW_ALBEDO,CAN_MODEL,                           &
                     CATCH,ECAN_TILE,INFIL_TILE,                        &
                     MELT_TILE,T_SURF_TILE,TILE_FRAC,CANOPY,            &
                     RGRAIN,SNOW_TILE,                                  &
                     CATCH_SNOW,SNOW_GRND,                              &

! Additional variables required for large-scale hydrology:
     &               L_TOP,L_PDM,FEXP,GAMTOT,TI_MEAN,TI_SIG,CS_CH4,     &
     &               DUN_ROFF,FSAT,FWETL,QBASE,QBASE_ZW,ZW,DRAIN,       &
     &               STHZW,A_FSAT,C_FSAT,A_FWET,C_FWET,FCH4_WETL,       &
     &               dim_cs1,L_SOIL_SAT_DOWN,l_triffid,                 &

! Timer diagnostics disabled
     &               Ltimer )

! Copy land points output back to full fields array.
        Do l = 1, land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          snow_depth(i,j) = LYING_SNOW(L)
          SNOWMELT(I,J) = SNOW_MELT(L)
        End Do

      End If

! map module prognostics back to JULES prognostics
        SNOWDEPTH_P     = SNOWDEPTH
      IF(NSMAX > 0)THEN
        RHO_SNOW_GRND_P = RHO_SNOW_GRND
        NSNOW_P         = NSNOW
        DS_P            = DS
        SICE_P          = SICE
        SLIQ_P          = SLIQ
        TSNOWLAYER_P    = TSNOW
        RHO_SNOW_P      = RHO_SNOW
        RGRAINL_P       = RGRAINL
      ENDIF

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 Hydrology',6)

!-------------------------------------------------------------------
! RIVER ROUTING


      IF ( L_RIVERS ) THEN

! Set the river routing to run on the 'last' PE as PE0 is very busy
        gather_pe_trip = n_proc-1

! Initialise diagnostics on non-river routing timesteps
          RIVEROUT = 0.0
          BOX_OUTFLOW = 0.0
          BOX_INFLOW = 0.0
! Initialise the gather fields at the beginning of TRIP timestep

        IF(A_STEPS_SINCE_RIV == 0)THEN
          TOT_SURF_RUNOFF=0.0
          TOT_SUB_RUNOFF=0.0
          ACC_LAKE_EVAP=0.0 
        ENDIF

        A_STEPS_SINCE_RIV = A_STEPS_SINCE_RIV + 1
        NSTEP_TRIP=INT(RIVER_STEP/TIMESTEP)
        IF (A_STEPS_SINCE_RIV == NSTEP_TRIP) THEN
          TRIP_CALL=.TRUE.
        ELSE
          TRIP_CALL=.FALSE.
        ENDIF

! Accumulate the runoff as Kg/m2/s over the TRIP period
        DO l = 1, land_points
          IF(surf_roff(L) <  0.0)THEN
            WRITE(6,*)'surf_roff(',L,')= ',surf_roff(L)
          ELSE
            TOT_SURF_RUNOFF(L) = TOT_SURF_RUNOFF(L) +                   &
     &                   (surf_roff(L)/ REAL(NSTEP_TRIP))
          ENDIF
        ENDDO
        DO l = 1, land_points
          IF(sub_surf_roff(L) <  0.0)THEN
            WRITE(6,*)'sub_surf_roff(',L,')= ',sub_surf_roff(L)
          ELSE
            TOT_SUB_RUNOFF(L) = TOT_SUB_RUNOFF(L) +                     &
     &                   (sub_surf_roff(L)/ REAL(NSTEP_TRIP))
          ENDIF
        ENDDO

        DO l = 1, land_points 
          j = (land_index(l)-1)/row_length +1 
          i = land_index(l) - (j-1)*row_length 
          ACC_LAKE_EVAP(I,J) = ACC_LAKE_EVAP(I,J) +                     & 
     &                   FRAC(L,7)*FQT_TILE(L,7)*TIMESTEP 
        ENDDO 

! detect first entry into river routing
      FIRST_ROUTING = .FALSE.
      IF (TIMESTEP_NUMBER == NSTEP_TRIP) FIRST_ROUTING = .TRUE.

        IF(TRIP_CALL)THEN

! DEPENDS ON: timer
        If (Ltimer) Call timer ('AP2 River Routing',5)

! If ATMOS fields are as Ocean (i.e. inverted NS) set INVERT_ATMOS
           INVERT_ATMOS = .FALSE.
           IF(.NOT.INVERT_OCEAN)INVERT_ATMOS = .TRUE.

! Calculate the Atmosphere gridbox areas
          DO J = 1, rows
            DO I = 1, row_length
              A_BOXAREAS(I,J) = r_theta_levels(i,j,0)                   &
                            * r_theta_levels(i,j,0)                     &
                            * delta_lambda * delta_phi                  &
                            * xx_cos_theta_latitude(i,j)
            ENDDO
          ENDDO

          SELECT CASE ( i_river_vn )

            CASE ( i_river_vn_1A )
          CALL RIV_INTCTL_1A(                                           &
       XPA, XUA, XVA, YPA, YUA, YVA,                                    &
       G_P_FIELD, G_R_FIELD, N_PROC, mype, RMDI,                        &
       GATHER_PE_TRIP,LAND_POINTS,LAND_INDEX,                           &
       INVERT_ATMOS, ROW_LENGTH, ROWS,                                  &
       GLOBAL_ROW_LENGTH, GLOBAL_ROWS,                                  &
       RIVER_ROW_LENGTH, RIVER_ROWS,                                    &
       GLOBAL_RIVER_ROW_LENGTH, GLOBAL_RIVER_ROWS,                      &
       FLANDG(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
       RIVER_STEP, RIVER_VEL, RIVER_MCOEF,                              &
       TRIVDIR, TRIVSEQ, TWATSTOR, A_BOXAREAS,                          &
        DELTA_PHI,FIRST_ROUTING,                                        &
        r_area, slope, flowobs1,r_inext,r_jnext,r_land,                 &
        substore,surfstore,flowin,bflowin,                              &
! IN/OUT accumulated runoff
       TOT_SURF_RUNOFF, TOT_SUB_RUNOFF,                                 &
! OUT
        BOX_OUTFLOW, BOX_INFLOW, RIVEROUT                               &
! Add inland basin arguments in call to rivctl
       ,INLANDOUT_ATMOS,INLANDOUT_RIV                                   & 
! Required for soil moisture correction for water conservation 
       ,DSM_LEVELS,ACC_LAKE_EVAP,SMVCST,SMVCWT                          & 
       ,SOIL_LAYER_MOISTURE(1:,DSM_LEVELS),STHU(1:,DSM_LEVELS)          & 
       )
            CASE ( i_river_vn_2A )
          CALL RIV_INTCTL_2A(                                           &
       XPA, XUA, XVA, YPA, YUA, YVA,                                    &
       G_P_FIELD, G_R_FIELD, N_PROC, mype, RMDI,                        &
       GATHER_PE_TRIP,LAND_POINTS,LAND_INDEX,                           &
       INVERT_ATMOS, ROW_LENGTH, ROWS,                                  &
       GLOBAL_ROW_LENGTH, GLOBAL_ROWS,                                  &
       RIVER_ROW_LENGTH, RIVER_ROWS,                                    &
       GLOBAL_RIVER_ROW_LENGTH, GLOBAL_RIVER_ROWS,                      &
       FLANDG(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
       RIVER_STEP, RIVER_VEL, RIVER_MCOEF,                              &
       TRIVDIR, TRIVSEQ, TWATSTOR, A_BOXAREAS,                          &
        DELTA_PHI,FIRST_ROUTING,                                        &
        r_area, slope, flowobs1,r_inext,r_jnext,r_land,                 &
        substore,surfstore,flowin,bflowin,                              &
! IN/OUT accumulated runoff
       TOT_SURF_RUNOFF, TOT_SUB_RUNOFF,                                 &
! OUT
        BOX_OUTFLOW, BOX_INFLOW, RIVEROUT                               &
! Add inland basin arguments in call to rivctl
       ,INLANDOUT_ATMOS,INLANDOUT_RIV                                   &
! Required for soil moisture correction for water conservation
       ,DSM_LEVELS,ACC_LAKE_EVAP,SMVCST,SMVCWT                          &
       ,SOIL_LAYER_MOISTURE(1:,DSM_LEVELS),STHU(1:,DSM_LEVELS)          &
       )

             CASE DEFAULT

          errorstatus = 10
          WRITE (message,'(A,I6,A)') 'River model type option ',       &
                             i_river_vn,' not recognised.'

          CALL Ereport ( RoutineName, errorstatus, message)

        END SELECT

! compress inland basin outputs to land points only
          if (L_INLAND) then
            do l = 1,land_points
              j = (land_index(l)-1)/row_length +1
              i = land_index(l) - (j-1)*row_length
              inlandout_atm(l) = inlandout_atmos(i,j)
            end do
          end if

! Mult RIVEROUT by the number of physics timesteps per River routing
! timestep as DAGHYD stores RIVEROUT every timestep. Non-routing
! timestep vals are passed in as 0.0

          A_STEPS_SINCE_RIV = 0

! DEPENDS ON: timer
      If (Ltimer) Call timer ('AP2 River Routing',6)


! -------------------------------------------------------------------
! Section RIVER Output diagnostics
! ------------------------------------------------------------------
! Check that river diagnostics requested this timestep
      If (error_code  ==  0  .and. sf(0,26) ) Then

        Call diagnostics_riv(                                           &
     &                       row_length, rows                           &
     &,                      river_row_length, river_rows               &
     &,                      at_extremity                               &
     &,                      at_extremity                               &
     &,                      RIVEROUT                                   &
     &,                      BOX_OUTFLOW, BOX_INFLOW                    &
! Put inland basin outflow in call to dagriv
     &,                      TWATSTOR,INLANDOUT_RIV                     &

     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork26                                                      &
     & )

      End If ! on error code .and. L_rivers .and. sf(0,26)
        ENDIF                        !TRIP_CALL

      ENDIF                        ! L_RIVERS (ATMOS)

! ------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section HYD.2 Output diagnostics
! ----------------------------------------------------------------------

! Check that hydrology diagnostics requested this timestep
      If (error_code  ==  0 .and. L_hydrology .and. sf(0,8) ) Then

! DEPENDS ON: diagnostics_hyd
        Call diagnostics_hyd(                                           &
     &                       row_length, rows, model_levels             &
     &,                      n_rows, global_row_length, global_rows     &
     &,                      halo_i, halo_j, offx, offy, mype           &
     &,                      n_proc, n_procx, n_procy                   &
     &,                      g_rows, g_row_length                       &
     &,                      at_extremity                               &
     &,                      land_points, dsm_levels                    &
 ! Put inland basin outflow in call to dagriv
     &,                      land_index,inlandout_atm                   &

     &,                      smc, surf_roff, sub_surf_roff              &
     &,                      lying_snow, snow_melt                      &
     &,                      canopy_gb,t_soil                           &
     &,                      soil_layer_moisture                        &
     &,                      ntiles, snomlt_surf_htf, sthu, sthf        &
     &,                      tot_tfall, snow_tile, melt_tile            &
     &,                      rgrain, land_sea_mask                      &
     &,                      dun_roff, drain, qbase, qbase_zw           &
     &,                      fch4_wetl                                  &
     &,                      fexp,gamtot,ti_mean,ti_sig                 &
     &,                      fsat,fwetl,zw,sthzw                        &
     &,                      timestep                                   &
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     & STASHwork8                                                       &
     & )

      End If ! on error code .and. L_hydrology .and. sf(0,8)

! ----------------------------------------------------------------------
! Section 19 -- VEGETATION DYNAMICS
! ----------------------------------------------------------------------

!-----------------------------------------------------------------------
! If leaf phenology is activated, check whether the atmosphere model
! has run an integer number of phenology calling periods.
!-----------------------------------------------------------------------
      PHENOL_CALL=1
      TRIFFID_CALL=1
      IF (L_PHENOL) THEN
        PHENOL_CALL = MOD ( FLOAT(A_STEP),(FLOAT(PHENOL_PERIOD)*        &
     &  (86400.0/TIMESTEP)) )
      ENDIF

      IF (L_TRIFFID) THEN


        NSTEP_TRIF=INT(86400.0*TRIFFID_PERIOD/TIMESTEP)
        IF (ASTEPS_SINCE_TRIFFID == NSTEP_TRIF) THEN
          TRIFFID_CALL=0
        ENDIF
      ENDIF

      IF ((PHENOL_CALL == 0).OR.(TRIFFID_CALL == 0)) THEN
! DEPENDS ON: veg_ctl
        CALL VEG_CTL(                                                   &
     &               row_length, rows, n_rows                           &
     &,              global_row_length, global_rows                     &
     &,              DIM_CS1, DIM_CS2                                   &
     &,              halo_i, halo_j, offx, offy, mype                   &
     &,              n_proc, n_procx, n_procy                           &
     &,              g_rows, g_row_length                               &
     &,              at_extremity                                       &
     &,              LAND_POINTS,LAND_INDEX,NTILES,CAN_MODEL            &
     &,              A_STEP,ASTEPS_SINCE_TRIFFID                        &
     &,              PHENOL_PERIOD,TRIFFID_PERIOD                       &
     &,              L_PHENOL,L_TRIFFID,L_TRIF_EQ                       &
     &,              TIMESTEP,FRAC_DISTURB,SATCON                       &
     &,              G_LEAF_ACC,G_LEAF_PHEN_ACC,NPP_FT_ACC              &
     &,              RESP_S_ACC,RESP_W_FT_ACC                           &
     &,              CS,FRAC,LAI_FT,SOIL_CLAY,CANHT_FT                  &
     &,              CATCH_SNOW,CATCH,INFIL_TILE,Z0_TILE                &
      ,              z0h_tile_bare                                      &
     &,                                                                 &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
     &               STASHwork19                                        &
     &                   )
      ENDIF

! ----------------------------------------------------------------------
! Set Coastal tiling dependent prognostics for output to D1 array:
! ----------------------------------------------------------------------
      IF(L_ctile)THEN
        DO J = 1, ROWS
          DO I = 1, ROW_LENGTH
            TSTAR_LAND_CTILE(I,J)=TSTAR_LAND(I,J)
            TSTAR_SICE_CTILE(I,J)=TSTAR_SICE(I,J)
          ENDDO
        ENDDO

        DO K = 1, NICE_USE
          DO J = 1, ROWS
            DO I = 1, ROW_LENGTH
              TSTAR_SICE_CAT_CTILE(I,J,K)=TSTAR_SICE_CAT(I,J,K)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      End If ! If ( CycleNo == NumCycles ) Then

      ! Copy ccw/lcbase values for Radiation if l_ccrad in use.
      ! See comments at declaration of ccw_out/lcbase_out
      IF (l_ccrad) THEN

        DO k=1, wet_levels
          DO j=1, rows
            DO i=1, row_length
              ccw_out(i,j,k) = ccw0(i,j,k)
            END DO
          END DO
        END DO

        DO j=1, rows
          DO i=1, row_length
            lcbase_out(i,j)  = lcbase0(i,j)
          END DO
        END DO

      END IF ! l_ccrad



      IF (lhook) CALL dr_hook('ATMOS_PHYSICS2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Atmos_Physics2
