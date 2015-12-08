! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Interface to Atmospheric Physics radiation code.

! Purpose:
!   This is the top-level radiation control routine. Major book-keeping
!   is carried out here. Separate calls are made to calculate SW and LW
!   radiative fluxes and heating rates on radiation time-steps.
!   Radiation increments are applied on all physics time-steps.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!-----------------------------------------------------------------------
SUBROUTINE glue_rad (                                                   &

! Parallel variables
  halo_i, halo_j, off_x, off_y, global_row_length, global_rows,         &
  proc_row_group, proc_col_group, at_extremity, n_proc,                 &
  n_procx, n_procy, neighbour, g_rows, g_row_length,                    &
  me, global_cloud_top,                                                 &

! model dimensions.
  row_length, rows, n_rows,                                             &
  model_levels, wet_model_levels, bl_levels,                            &
  ozone_levels, cloud_levels, n_cca_levels,                             &
  ntiles, land_field, nice_use, dust_dim1, dust_dim2, biogenic_dim1,    &
  biogenic_dim2, sulp_dim1, sulp_dim2, soot_dim1, soot_dim2,            &
  bmass_dim1, bmass_dim2, salt_dim1, salt_dim2, salt_dim3,              &
  co2_dim_len, co2_dim_row, co2_dim2, arcl_dim1, arcl_dim2,             &
  n_arcl_species, n_arcl_compnts, i_arcl_compnts,                       &
  ocff_dim1, ocff_dim2, nitrate_dim1, nitrate_dim2,                     &
  ukca_dim1, ukca_dim2, n_ukca_mode, n_ukca_cpnt,                       &

! Model switches
  model_domain, l_rad_step,                                             &
  l_rad_step_diag, l_rad_step_prog,                                     &
  l_forcing, l_timestep, l_radiance,                                    &
  l_cal360, l_sec_var, l_eqt,                                           &
  l_inhom_cloud, l_dust, l_use_dust, l_use_biogenic,                    &
  l_sulpc_so2, l_use_sulpc_direct,                                      &
  l_soot, l_use_soot_direct, l_use_soot_indirect,                       & 
  l_biomass, l_use_bmass_direct, l_use_bmass_indirect,                  &
  l_ocff, l_use_ocff_direct, l_use_ocff_indirect,                       &
  l_use_sulpc_indirect_sw, l_use_sulpc_indirect_lw,                     &
  l_nitrate, l_use_nitrate_direct, l_use_nitrate_indirect,              &
  l_emcorr,l_climat_aerosol,l_clim_aero_hgt,l_hadgem1_clim_aero,        &
  ltimer, l_ssice_albedo,                                               &
  l_snow_albedo,                                                        &
  l_sice_meltponds,l_sice_scattering,l_sice_hadgem1a,l_cice_alb,        &
  l_use_seasalt_indirect,                                               &
  l_use_seasalt_direct,                                                 &
  l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,                                &
  l_pc2, l_mixing_ratio,                                                &
  l_murk_rad,                                                           &
  l_rad_deg,                                                            &
  l_co2_interactive,                                                    &
  l_ukca, l_ukca_radaer,                                                &
  l_use_arcl, l_use_spec_sea, l_mod_barker_albedo, l_mod_k_flux,        &
  l_cosp_in,                                                            &
  can_rad_mod,                                                          &

! model Parameters
  a_sw_segments, a_sw_seg_size, a_lw_segments, a_lw_seg_size,           &
  inhom_cloud_sw, inhom_cloud_lw, dp_corr_strat, dp_corr_conv,          &
  co2_mmr, o2mmr, n2o_mix_ratio, ch4_mix_ratio,                         &
  cfc11_mix_ratio, cfc12_mix_ratio,                                     &
  c113mmr, c114mmr, hcfc22mmr, hfc125mmr, hfc134ammr,                   &
  sw_alpham, sw_alphac, sw_alphab, sw_dtice,                            &
  dt_bare,dalb_bare_wet,pen_rad_frac,sw_beta,                           &
  min_trop_level, max_trop_level,                                       &
  ntot_land, ntot_sea, aero_bl_levels,                                  &

! in coordinate information
  rho_r2,                                                               &
  delta_lambda, delta_phi,                                              &
  lat_rot_np, long_rot_np,                                              &
  trindx,                                                               &

! in time stepping information.
  timestep,radiation_timestep,                                          &
  radiation_tstep_diag,                                                 &
  radiation_tstep_prog,                                                 &
  val_year, val_day_number, val_hour, val_minute,                       &
  val_second, timestep_number,                                          &
  previous_time,                                                        &

! diagnostic info
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
  stashwork1,                                                           &
  stashwork2,                                                           &

! SCM diagnostics switches (dummy in full UM)
  nscmdpkgs,l_scmdiags,                                                 &

! in data fields.
  p_star,                                                               &
  p_layer_boundaries, p_layer_centres,                                  &
  p, p_theta_levels,                                                    &
  exner_rho_levels, exner_theta_levels,                                 &
  land_sea_mask, fland, land0p5, l_ctile,                               &
  t_surf,tstar_sea,tstar_sice_cat,area_cloud_fraction,                  &
  dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6,          &
  biogenic, so4_aitken, so4_accu, so4_diss,                             &
  soot_new, soot_agd, soot_cld, bmass_new, bmass_agd, bmass_cld,        &
  ocff_new, ocff_agd, ocff_cld, nitr_acc, nitr_diss, aerosol, arcl,     &
  ukca_radaer, sea_salt_film, sea_salt_jet,                             &
  co2_3d, frac_control, n_drop_pot,                                     &

! chemical greenhouse gas fields
  ngrgas, grgas_field,                                                  &

! ancillary fields and fields needed to be kept from timestep to
! timestep
  snow_depth, snow_depth_sea_cat, ice_fract, ice_fract_cat,             &
  ice_thick_cat, rgrain, soot,                                          &
  cca, ccb, cct, cclwp,ccw,lcbase,                                      &
  ozone, sw_incs, lw_incs, dirpar_inc,                                  &
  o3_trop_level, o3_trop_height,                                        &
  t_trop_level, t_trop_height, zh,                                      &
  land_index, albsoil, albobs_sw, albobs_vis, albobs_nir,               &
  lai, snow_tile, frac, tstar_tile, z0_tile,                            &
  dOLR_rts, lw_down, sw_tile_rts,                                       &
  land_alb,sice_alb,                                                    &
  es_space_interp, a_sw_radstep, a_lw_radstep,                          &
  a_sw_radstep_diag,a_sw_radstep_prog,                                  &
  a_lw_radstep_diag, a_lw_radstep_prog, rad_mask,                       &

! in/out
  t_n, q_n, qcl_n, qcf_n, cf_n, cfl_n, cff_n,                           &
  qcf2_n, qrain_n, qgraup_n,                                            &
  t_inc, q_inc, qcl_inc, cf_inc, cfl_inc,                               &
  sum_eng_fluxes,                                                       &
  cos_zenith_angle,                                                     &

! out.
  photosynth_act_rad, rad_hr, dOLR, sw_tile,                            &

! COSP arguments
  cosp_gbx,                                                             &
  
! error information
  error_code  )

  USE dynamics_grid_mod, ONLY: l_vatpoles

  USE nstypes
  USE earth_constants_mod, ONLY: two_omega
  USE level_heights_mod, ONLY:                                          &
    r_theta_levels, r_rho_levels
  USE trignometric_mod, ONLY:                                           &
    true_longitude, true_latitude,                                      &
    cos_theta_latitude,                                                 &
    fv_cos_theta_latitude

! Modules required for multiple calls to radiation
  USE max_calls
  USE diagoffset

! Modules required for reading in spectral files and control options
! as structures
  USE spec_sw_lw
  USE sw_control_struct
  USE lw_control_struct

! Modules required for radiative forcing
  USE mm_ratios
  USE coradoca, ONLY: c2c_o2, c2c_o3, c2c_co2, c2c_n2o, c2c_ch4,       &
       c2c_cfc11, c2c_cfc12, c2c_c113, c2c_hcfc22, c2c_hfc125,         &
       c2c_hfc134, c2c_aerosol, c2c_sulpc_d, c2c_seas_d, c2c_soot_d,   &
       c2c_bmb_d, c2c_ocff_d, c2c_land_s, c2c_all,                     &
       c2c_wmg, c2c_nitr_d, c2c_dust_d, c2c_biog_d, c2c_ukca_d

! Modules required for the orography scheme
  USE solinc_data, ONLY:                                                &
    sol_bearing, orog_corr, f_orog,                                     &
    slope_aspect, slope_angle, l_orog

! Modules for diagnostics
  USE sw_diag_mod
  USE lw_diag_mod

! Module for radiation switches
  USE rad_input_mod, ONLY: l_rad_perturb, l_rad_szacor,                 &
    l_rad_snow_emis, l_t_land_nosnow, l_quad_t_coast,                   &
    l_t_rad_solid, sc

! Modules for JULES
  USE switches, ONLY: l_aggregate, l_flake_model, l_dolr_land_black   
   
  USE nvegparm, ONLY: emis_nvg    
  USE pftparm,  ONLY: emis_pft    
  USE surf_param,  ONLY: emis_sea, emis_sice
  USE snow_param, ONLY : rho_snow_const   
  USE ancil_info, ONLY:                                                 &
    ssi_pts, sice_pts, sice_pts_ncat,                                   &
    ssi_index, sice_index, sice_index_ncat
  USE lake_mod, ONLY: lake_albedo

  USE fluxes, ONLY: sw_sice, sw_sice_rts, alb_sice

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

! Module for UKCA-MODE aerosol interaction with radiation
  USE ukca_radaer_struct_mod

! Module for two-bin or six-bin dust switch
  USE dust_parameters_mod, ONLY: l_twobin_dust

! Module for the aerosol climatologies
 USE arcl_mod,                ONLY: npd_arcl_compnts, npd_arcl_species, &
                                    ip_arcl_sulp_ac, ip_arcl_sulp_ak

  USE ereport_mod, ONLY : ereport
  USE um_parparams, ONLY: pnorth, psouth
  USE Field_Types
  USE domain_params
! Modules for COSP
  USE cosp_types_mod, ONLY: cosp_gridbox
  USE dimfix3a_mod, ONLY: npd_cloud_component, npd_cloud_type,          &
                          npd_cloud_representation, npd_overlap_coeff,  &
                          npd_source_coeff, npd_region

  USE segments_mod, ONLY:                                               &
     segment_type, meta_segment_type,                                   &
     segments_mod_seg_meta, segments_mod_segments

  USE science_fixes_mod, ONLY: l_rm_neg_par

  USE Submodel_Mod

!$ USE omp_lib
  use cable_data_mod, ONLY : cable_control5, cable_glue_rad_init
  IMPLICIT NONE

! Segmentation variables
  TYPE(segment_type),     ALLOCATABLE  :: segments(:)
  TYPE(meta_segment_type)              :: meta_segments
  INTEGER                              :: ipar
  INTEGER                              :: num_parallel_sw
  INTEGER                              :: num_parallel_isccp
  INTEGER                              :: num_parallel_lw

! Arguments with intent in. ie: input variables.

! Parallel setup variables
  INTEGER ::                                                            &
    halo_i,                                                             &
!     Size of halo in i direction.
    halo_j,                                                             &
!     Size of halo in j direction.
    off_x,                                                              &
!     Size of small halo in i
    off_y,                                                              &
!     Size of small halo in j.
    global_row_length,                                                  &
!     number of points on a row
    global_rows,                                                        &
!     NUMBER OF global rows
    proc_row_group,                                                     &
!     Group id for processors on the same row
    proc_col_group,                                                     &
!     Group id for processors on the same column
    n_proc,                                                             &
!     Total number of processors
    n_procx,                                                            &
!     Number of processors in longitude
    n_procy,                                                            &
!     Number of processors in latitude
    neighbour(4),                                                       &
!     Array with the Ids of the four neighbours
!     in the horizontal plane
    g_rows (0:n_proc-1),                                                &
    g_row_length (0:n_proc-1),                                          &
    me,                                                                 &
!     My processor number
    global_cloud_top
!     Global topmost cloudy level

  LOGICAL ::                                                            &
    at_extremity(4)
!     Indicates if this processor is at north, south
!     east or west of the processor grid

! Parameters
! STPARAM
!
!  Purpose: Meaningful PARAMETER names for STASH processing routines.
!           Both a long name and short name have been declared, to
!           reduce the existence of "magic" numbers in STASH.
!           Format is that first the address of the item is declare in
!           both long and short form. example is;
!             integer st_item_code,s_item  !Item number (declaration)
!             parameter(st_item_code=3,s_item=3)
!
!  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!  Logical components covered: D70
!
!  Project task: D7
!
!  External documentation:
!    Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                 system (STASH)
!--------------------------------------------------------------

      ! Internal model number address
      INTEGER,PARAMETER:: st_model_code = 28
      INTEGER,PARAMETER:: s_modl        = 28

      ! Section Number address
      INTEGER,PARAMETER:: st_sect_no_code = 2
      INTEGER,PARAMETER:: s_sect          = 2
      INTEGER,PARAMETER:: st_sect_code    = 2

      INTEGER,PARAMETER:: st_item_code=1,s_item=1 ! Item number address

      ! Processing Code address
      INTEGER,PARAMETER:: st_proc_no_code=3,s_proc=3

      ! subsidiary codes for st_proc_no_code now
      INTEGER,PARAMETER:: st_replace_code=1
      INTEGER,PARAMETER:: st_accum_code=2
      INTEGER,PARAMETER:: st_time_mean_code=3
      INTEGER,PARAMETER:: st_time_series_code=4
      INTEGER,PARAMETER:: st_max_code=5
      INTEGER,PARAMETER:: st_min_code=6
      INTEGER,PARAMETER:: st_append_traj_code=7
      INTEGER,PARAMETER:: st_time_series_mean=8
      INTEGER,PARAMETER:: st_variance_code=9

      ! Frequency (Input & output) addres
      INTEGER,PARAMETER:: st_freq_code=4,s_freq=4

      ! Offset for sampling
      INTEGER,PARAMETER:: st_offset_code=30,s_offs=30

      ! start timestep address
      INTEGER,PARAMETER:: st_start_time_code=5,s_times=5

      ! end timestep address
      INTEGER,PARAMETER:: st_end_time_code=6,s_timee=6

      ! period in timesteps address
      INTEGER,PARAMETER:: st_period_code=7,s_period=7

      ! infinite end/period value
      INTEGER,PARAMETER:: st_infinite_time=-1

      INTEGER,PARAMETER:: st_end_of_list=-1 !end-of-list marker in times

      ! grid point stuff
      ! gridpoint info address
      INTEGER,PARAMETER:: st_gridpoint_code=8,s_grid=8

      ! now subsid grid point stuff
      ! no masking done
      INTEGER,PARAMETER:: stash_null_mask_code=1,s_nomask=1

      ! land mask conds
      INTEGER,PARAMETER:: stash_land_mask_code=2,s_lndms=2

      ! sea mask code
      INTEGER,PARAMETER:: stash_sea_mask_code=3,s_seams =3

      ! processing options

      ! size of block for gridpoint code
      INTEGER,PARAMETER:: block_size=10

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: extract_base=block_size*0

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: extract_top=block_size*1

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_base=block_size*1

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_top=block_size*2

      ! max code for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_base=block_size*2

      ! base codes for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_top=block_size*3

      ! max code for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_base=block_size*3

      ! base codes for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_top=block_size*4

      ! max code for field mean subroutine
      INTEGER,PARAMETER:: field_mean_base=block_size*4

      ! base codes for field mean subroutine
      INTEGER,PARAMETER:: field_mean_top=block_size*5

      ! max code for global mean subroutine
      INTEGER,PARAMETER:: global_mean_base=block_size*5

      ! base codes for global mean subroutine
      INTEGER,PARAMETER:: global_mean_top=block_size*6

      ! Weighting

      ! weighting info address
      INTEGER,PARAMETER:: st_weight_code=9,s_weight=9

      INTEGER,PARAMETER:: stash_weight_null_code  =0,s_noweight  =0
      INTEGER,PARAMETER:: stash_weight_area_code  =1,s_areaweight=1
      INTEGER,PARAMETER:: stash_weight_volume_code=2,s_volweight =2
      INTEGER,PARAMETER:: stash_weight_mass_code  =3,s_massweight=3

      ! Domain definition

      ! row addresses
      INTEGER,PARAMETER:: st_north_code=12,s_north=12
      INTEGER,PARAMETER:: st_south_code=13,s_south=13
      INTEGER,PARAMETER:: st_west_code =14,s_west =14
      INTEGER,PARAMETER:: st_east_code =15,s_east =15

      ! Levels

      ! input bottom level address
      INTEGER,PARAMETER:: st_input_bottom=10,s_bottom =10

      ! special code
      INTEGER,PARAMETER:: st_special_code=100,s_special=100

      ! input top level address
      INTEGER,PARAMETER:: st_input_top=11,s_top=11

      ! output bottom level address
      INTEGER,PARAMETER:: st_output_bottom=21,s_outbot=21

      ! output top level address
      INTEGER,PARAMETER:: st_output_top=22,s_outtop=22

      INTEGER,PARAMETER:: st_model_level_code=1,s_model=1

      ! code for pressure leve
      INTEGER,PARAMETER:: st_pressure_level_code=2,s_press=2

      ! code for height levels
      INTEGER,PARAMETER:: st_height_level_code=3,s_height=3

      ! input code addres
      INTEGER,PARAMETER:: st_input_code=16,s_input=16

      ! input length of diagnostic address
      INTEGER,PARAMETER:: st_input_length=17,s_length=17

      ! output code address
      INTEGER,PARAMETER:: st_output_code=18,s_output=18

      ! Pointer to D1 addressing information
      ! Pos of item in D1 for relevant submodel
      INTEGER,PARAMETER:: st_position_in_d1=29,st_d1pos=29

      ! Output destination options

      INTEGER,PARAMETER:: st_dump=1
      INTEGER,PARAMETER:: st_secondary=2

      ! output length of diagnostic address
      INTEGER,PARAMETER:: st_output_length=19,s_outlen=19
         integer st_dump_output_length,s_doutlen ! output length on
         parameter(st_dump_output_length=32,s_doutlen=32)  ! dump
         integer st_dump_level_output_length,s_dlevoutlen
         parameter(st_dump_level_output_length=33,s_dlevoutlen=33)
! output length of a single level on dump

         integer st_output_addr,s_outadd ! start locn of diag after stas
         parameter(st_output_addr=20,s_outadd=20)       ! output address
         integer st_dump_output_addr,s_doutadd ! output address on
         parameter(st_dump_output_addr=31,s_doutadd=31)  ! dump

      ! ptr to dump lookup header address
      INTEGER,PARAMETER:: st_lookup_ptr=23

      ! ptr into stash_series where control data address
      INTEGER,PARAMETER:: st_series_ptr=24

      ! subsid stuff for time series
      INTEGER,PARAMETER:: series_grid_type=1
      INTEGER,PARAMETER:: series_grid_code=0
      INTEGER,PARAMETER:: series_long_code=1
      INTEGER,PARAMETER:: series_size=2
      INTEGER,PARAMETER:: series_proc_code=3
      INTEGER,PARAMETER:: series_north=4
      INTEGER,PARAMETER:: series_south=5
      INTEGER,PARAMETER:: series_west=6
      INTEGER,PARAMETER:: series_east=7
      INTEGER,PARAMETER:: series_list_start=8
      INTEGER,PARAMETER:: series_list_end=9
      INTEGER,PARAMETER:: record_size=9

      ! Miscellaneous parameters

      ! system/user tag field in stlist address
      INTEGER,PARAMETER:: st_macrotag=25

      ! Pseudo-level list pointers

      ! pseudo-levels input list address
      INTEGER,PARAMETER:: st_pseudo_in=26

      ! pseudo-levels output list address
      INTEGER,PARAMETER:: st_pseudo_out=27

      ! Internal horizontal gridtype codes common to all diagnostics

      INTEGER,PARAMETER:: st_tp_grid =1 ! T-p grid
      INTEGER,PARAMETER:: st_uv_grid =2 ! u-v grid
      INTEGER,PARAMETER:: st_cu_grid =3 ! C-grid u point
      INTEGER,PARAMETER:: st_cv_grid =4 ! C-grid v point
      INTEGER,PARAMETER:: st_zt_grid =5 ! Zonal T-grid
      INTEGER,PARAMETER:: st_zu_grid =6 ! Zonal u-grid
      INTEGER,PARAMETER:: st_mt_grid =7 ! Meridional T-grid
      INTEGER,PARAMETER:: st_mu_grid =8 ! Meridional u-grid
      INTEGER,PARAMETER:: st_riv_grid= 23    ! river_routing grid
      INTEGER,PARAMETER:: st_scalar  =9 ! Scalar (ie. single value)
      INTEGER,PARAMETER:: st_wam_all= 60    ! Wam Field on Full Grid
      INTEGER,PARAMETER:: st_wam_sea= 62    ! Wam Field on Sea Points

! STPARAM end

! Model dimensions
  INTEGER ::                                                            &
    row_length,                                                         &
    rows,                                                               &
    n_rows,                                                             &
    model_levels,                                                       &
    wet_model_levels,                                                   &
    bl_levels,                                                          &
    ozone_levels,                                                       &
    cloud_levels,                                                       &
    n_cca_levels,                                                       &
    ntiles,                                                             &
    land_field,                                                         &
    nice_use,                                                           &
    dust_dim1,                                                          &
!     dimensions of mineral dust arrays
    dust_dim2,                                                          &
    biogenic_dim1,                                                      &
!     dimensions of biogenic aerosol arrays
    biogenic_dim2,                                                      &
    sulp_dim1,                                                          &
!     dimensions of S Cyc arrays
    sulp_dim2,                                                          &
    soot_dim1,                                                          &
!     dimensions of soot arrays
    soot_dim2,                                                          &
    bmass_dim1,                                                         &
!     dimensions of biomass arrays
    bmass_dim2,                                                         &
    ocff_dim1,                                                          &
    ocff_dim2,                                                          &
!     dimensions of OCFF arrays
    salt_dim1,                                                          &
!     dimensions of sea-salt arrays
    salt_dim2,                                                          &
    salt_dim3,                                                          &
    co2_dim1,                                                           &
!     dimensions of CO2 array
    co2_dim2,                                                           &
    co2_dim_len,                                                        &
    co2_dim_row,                                                        &
!     dimensions of aerosol clim for NWP
    arcl_dim1,                                                          &
    arcl_dim2,                                                          &
!     dimensions of nitrate arrays
    nitrate_dim1,                                                       &
    nitrate_dim2,                                                       &
!     dimensions of UKCA_RADAER arrays
    ukca_dim1,                                                          &
    ukca_dim2,                                                          &
    n_ukca_mode,                                                        &
    n_ukca_cpnt

! Model switches
  INTEGER ::                                                            &
    can_rad_mod
!     Which canopy radiation model we're using

  INTEGER ::                                                            &
    model_domain

  LOGICAL ::                                                            &
    l_rad_step,                                                         &
!     dummy variable in version 3Z
    l_rad_step_diag,                                                    &
!     true if fast radiation timestep
    l_rad_step_prog,                                                    &
!     true if slow radiation timestep
    l_timestep,                                                         &
!     true if improved timestepping is used
    l_forcing,                                                          &
!     true if radiative frocing required
    l_radiance,                                                         &
!     true if radiances are to be calculated
    l_cal360,                                                           &
!     true if using 360 day calender
    l_sec_var,                                                          &
!     True if using time varying astronomy
    l_eqt,                                                              &
!     true if including the equation of time
    l_inhom_cloud,                                                      &
!     switch to simulate inhomogeneous cloud
    l_climat_aerosol,                                                   &
!     true
    l_clim_aero_hgt,                                                    &
!     True if using the prognostic BL depth to
!     specify the BL component of the aerosol
!     climatology.
    l_hadgem1_clim_aero,                                                &
!     True if using HadGEM1 setting for climatological aerosols
    l_dust,                                                             &
!     mineral dust available for use (for direct effect or diagnostics)
    l_use_dust,                                                         &
!     include mineral dust direct effect
    l_use_biogenic,                                                     &
!     include biogenic aerosol direct effect
    l_sulpc_so2,                                                        &
!     Sulphur C available for use (for direct/indirect or diagnostics)
    l_use_sulpc_direct,                                                 &
!     Include direct effect of sulphate
    l_soot,                                                             &
!     Soot available for use (for direct/indirect or diagnostics)
    l_use_soot_direct,                                                  &
!     Include direct effect of soot
    l_biomass,                                                          &
!     biomass smoke available for use (for direct/indirect or diagnostics)
    l_use_bmass_direct,                                                 &
!     Include direct effect of biomass smoke
    l_ocff,                                                             &
!     OCFF available for use (for direct/indirect or diagnostics)
    l_use_ocff_direct,                                                  &
!     Include direct effect of OCFF
    l_nitrate,                                                          &
!     nitrate available for use (for direct/indirect or diagnostics)
    l_use_nitrate_direct,                                               &
!     Include direct effect of nitrate
    l_use_sulpc_indirect_sw,                                            &
    l_use_sulpc_indirect_lw,                                            &
!     Use sulphate aerosol to determine
!     cloud droplet number in SW and LW
!     radiation respectively.
    l_use_soot_indirect,                                                &
!     Include indirect effect of soot
    l_use_bmass_indirect,                                               &
!     Include indirect effect of biomass
    l_use_ocff_indirect,                                                &
!     Include indirect effect of OCFF
    l_use_nitrate_indirect,                                             &
!     Include indirect effect of nitrate
    l_emcorr,                                                           &
!     true if energy correction scheme is to be used.
    ltimer,                                                             &
!     true then output some timing information
    l_ssice_albedo,                                                     &
!     Switch on the effect of snow on sea-ice albedo
    l_sice_meltponds,                                                   &
!     true if sea ice meltponds albedo scheme
    l_sice_scattering,                                                  &
!     true if sea ice internal scattering scheme
    l_sice_hadgem1a,                                                    &
!     true if HadGEM1 albedo bug corrected
    l_cice_alb,                                                         &
!     true if sea ice CICE albedo scheme
    l_snow_albedo,                                                      &
!     True if spectral albedo scheme selected
    l_use_seasalt_indirect,                                             &
!     Switch for indirect effect of sea-salt.
    l_use_seasalt_direct,                                               &
!     Switch for direct effect of sea-salt.
    l_mcr_qcf2,                                                         &
!     Use second ice category
    l_mcr_qrain,                                                        &
!     Use prognostic rain
    l_mcr_qgraup,                                                       &
!     Use graupel
    l_pc2,                                                              &
!     True if using the PC2 cloud scheme
    l_mixing_ratio,                                                     &
!     True if using mixing ratios
    l_murk_rad,                                                         &
!     True if using radiative effects of 'murk'.
    l_ctile,                                                            &
!     coastal tiling switch
    l_rad_deg,                                                          &
!     Controls the spatial degradation of E-S code
    l_co2_interactive,                                                  &
!     Controls the use of 3D CO2 field
    l_ukca,                                                             &
!     Switch for UKCA sub-model
    l_ukca_radaer,                                                      &
!     Switch for interaction UKCA aerosols-radn
    l_use_spec_sea,                                                     &
!     Switch for spectrally dependent sea albedos
    l_mod_barker_albedo
!     Use modified Barker albedo (open sea)

! Is COSP requested?
  LOGICAL, INTENT(IN) :: l_cosp_in

! Use modulus of fluxes to remove negative effective extinctions
  LOGICAL, INTENT(IN) :: l_mod_k_flux

! model parameters
  REAL ::                                                               &
    timestep,                                                           &
    radiation_timestep,                                                 &
!     Dummy variable in version 3Z
    radiation_tstep,                                                    &
    radiation_tstep_diag,                                               &
    radiation_tstep_prog,                                               &
    co2_mmr,                                                            &
!     set equal to co2_start
    o2mmr,                                                              &
    n2o_mix_ratio,                                                      &
    ch4_mix_ratio,                                                      &
    cfc11_mix_ratio,                                                    &
    cfc12_mix_ratio,                                                    &
    sw_alphab,                                                          &
    sw_alphac,                                                          &
    sw_alpham,                                                          &
    sw_dtice,                                                           &
    dt_bare,                                                            &
    dalb_bare_wet,                                                      &
    pen_rad_frac,                                                       &
    sw_beta

  REAL :: c113mmr     ! CFC113 mmr
  REAL :: c114mmr     ! CFC114 mmr
  REAL :: hcfc22mmr   ! HCFC22 mmr
  REAL :: hfc125mmr   ! HFC125 mmr
  REAL :: hfc134ammr  ! HFC134A mmr

  INTEGER ::                                                            &
    min_trop_level,                                                     &
!     Lowest permitted level for the tropopause
!     used for radiative purposes.
    max_trop_level
!     Highest permitted level for the tropopause
!     used for radiative purposes.

  INTEGER ::                                                            &
    a_sw_segments,                                                      &
    a_sw_seg_size,                                                      &
    a_lw_segments,                                                      &
    a_lw_seg_size,                                                      &
    a_sw_radstep,                                                       &
!     Dummy variable in version 3Z
    a_lw_radstep,                                                       &
!     Dummy variable in version 3Z
    a_sw_radstep_diag,                                                  &
    a_sw_radstep_prog,                                                  &
    a_lw_radstep_diag,                                                  &
    a_lw_radstep_prog,                                                  &
    aero_bl_levels
!     Common number of layers taken to be occupied by the
!     boundary-layer aerosol if the boundary layer
!     depth is not used to determine the number separately
!     at each grid-point


! Scaling factors to simulate inhomogeneous cloud.
  REAL :: inhom_cloud_sw(npd_cloud_component)
  REAL :: inhom_cloud_lw(npd_cloud_component)

! Decorrelation pressure scale for large scale cloud
  REAL :: dp_corr_strat
! Decorrelation pressure scale for convective cloud
  REAL :: dp_corr_conv

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

! Diagnostics info
  REAL ::                                                               &
    stashwork1(*),                                                      &
    stashwork2(*)
!     STASH workspace


! Data arrays
  REAL ::                                                               &
    p_layer_boundaries(row_length, rows, 0:model_levels),               &
!     pressure at layer boundaries. Same as p except at
!     bottom level = pstar, and at top = 0.
    p_layer_centres(row_length, rows, 0:model_levels),                  &
!     pressure at layer centres. Same as p_theta_levels
!     except bottom level = pstar, and at top = 0.
    p_star(row_length, rows),                                           &
    p_theta_levels(1-off_x:row_length+off_x,                            &
                   1-off_y:rows+off_y, model_levels),                   &
    p(1-off_x:row_length+off_x,                                         &
      1-off_y:rows+off_y, model_levels),                                &
    exner_theta_levels(1-off_x:row_length+off_x,                        &
                       1-off_y:rows+off_y, model_levels),               &
    exner_rho_levels(1-off_x:row_length+off_x,                          &
                       1-off_y:rows+off_y, model_levels)

  REAL    :: ccw   (row_length, rows, wet_model_levels)
  INTEGER :: lcbase(row_length, rows)

  REAL ::                                                               &
    cca (row_length, rows, n_cca_levels),                               &
    cclwp(row_length, rows),                                            &
!     condensed water path (KG/M**2)
    area_cloud_fraction(row_length, rows, wet_model_levels)

  INTEGER ::                                                            &
    ccb (row_length, rows),                                             &
    cct (row_length, rows)

! ancillary arrays and fields required to be saved from timestep to
! timestep.

  REAL ::                                                               &
    t_surf(row_length, rows),                                           &
    tstar_sea(row_length,rows),                                         &
!     IN Open sea sfc temperature (K).
    tstar_sice_cat(row_length,rows,nice_use)
!     IN Sea-ice sfc temperature (K).

  REAL ::                                                               &
    ice_fract_cat(row_length, rows, nice_use),                          &
!     Area fraction of sea ice categories
    ice_thick_cat(row_length, rows, nice_use)
!     Effective thickness of each sea ice categories


  LOGICAL ::                                                            &
    land_sea_mask(row_length, rows),                                    &
    land0p5(row_length, rows),                                          &
!     A mask set to .TRUE. if the fraction of land in the grid-box
!     exceeds 0.5.
    rad_mask(row_length, rows)
!     A mask which ensures a chequerboard pattern of radiation
!     calculations over the whole domain (not just one PE)


  REAL ::                                                               &
    fland(land_field)
!     Fractional amount of land at each land point

  REAL ::                                                               &
    ntot_land,                                                          &
!     Number of droplets over land / m-3
    ntot_sea
!     Number of droplets over sea / m-3

  REAL ::                                                               &
    snow_depth (row_length, rows),                                      &
!     snow/qrclim.snow.(month)
    snow_depth_sea_cat (row_length, rows, nice_use)
!     snow depth on sea ice categories

  REAL ::                                                               &
    ice_fract (row_length, rows),                                       &
!     ice/qrclim.ice.(month)
    soot (row_length, rows)

! Input ancillary data:
  REAL, INTENT(IN) ::                                                   &
    albsoil(land_field),                                                &
    albobs_sw(land_field),                                              &
    albobs_vis(land_field),                                             &
    albobs_nir(land_field),                                             &
    lai(land_field, npft)

  REAL ::                                                               &
    rgrain(land_field, ntiles),                                         &
    snow_tile(land_field, ntiles),                                      &
    frac(land_field, ntype),                                            &
    frac_control(land_field, ntype),                                    &
    tstar_tile(land_field, ntiles),                                     &
    z0_tile(land_field, ntiles),                                        &
    dOLR_rts(row_length, rows),                                         &
!     TOA - surface upward LW
    lw_down(row_length, rows),                                          &
!     Surface downward LW
    sw_tile_rts(land_field, ntiles),                                    &
!     Surface net SW on land tiles
    land_alb(row_length, rows),                                         &
!     Mean land albedo
    sice_alb(row_length, rows),                                         &
!     Mean sea-ice albedo
    es_space_interp(4, row_length, rows)
!     Coeffs for spatial interpolation of radiation quantities

! Number of requested species within the climatology
  INTEGER :: n_arcl_species

! Corresponding number of requested components
  INTEGER :: n_arcl_compnts

! Model switch for each species
  LOGICAL :: l_use_arcl(npd_arcl_species)

! Array indices of components
  INTEGER :: i_arcl_compnts(npd_arcl_compnts)

! Mass-mixing ratios
  REAL ::                                                               &
    arcl(row_length, rows, model_levels, n_arcl_compnts)

! UKCA_RADAER structure: Interaction between UKCA-MODE aerosols
!                        and radiation
  TYPE (ukca_radaer_struct), INTENT(IN) :: ukca_radaer


  INTEGER ::                                                            &
    land_index(land_field)

  REAL ::                                                               &
    ozone(row_length, rows, ozone_levels),                              &
    o3_trop_level(row_length,rows),                                     &
    o3_trop_height(row_length,rows),                                    &
    t_trop_level(row_length,rows),                                      &
    t_trop_height(row_length,rows),                                     &
    sw_incs(row_length, rows, 0:model_levels+1),                        &
    lw_incs(row_length, rows, 0:model_levels),                          &
    dirpar_inc(row_length, rows),                                       &
    zh(row_length, rows),                                               &
!     boundary layer height
    co2_3d(row_length, rows, model_levels)

! Potential droplet number (includes values where cloud not present)
  REAL, INTENT(IN) :: n_drop_pot(row_length, rows, model_levels)

! chemical greenhouse gas fields
  INTEGER, INTENT(IN) :: ngrgas
  REAL, INTENT(IN) :: grgas_field(row_length,rows,model_levels,ngrgas)

! Co-ordinate arrays
  REAL ::                                                               &
    rho_r2(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels), &
!     Density*radius^2
    delta_lambda,                                                       &
    delta_phi


! Allocatable arrays for diagnostic variables - required to save memory
! use during radiation routines
  REAL, ALLOCATABLE :: t_incr_diagnostic(:,:,:)
!     temperature increment for STASH


! time information for current timestep
  INTEGER ::                                                            &
    val_year,                                                           &
    val_day_number,                                                     &
    val_hour,                                                           &
    val_minute,                                                         &
    val_second,                                                         &
    timestep_number,                                                    &
    previous_time(7)

  REAL ::                                                               &
    recip_timestep

! Diagnostic variables
  REAL ::                                                               &
    lat_rot_np,                                                         &
    long_rot_np

! Level of tropopause
  INTEGER, INTENT(IN) :: trindx(row_length, rows)

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER ::                                                            &
    nscmdpkgs
!     No of diagnostics packages

  LOGICAL ::                                                            &
    l_scmdiags(nscmdpkgs)
!     Logicals for diagnostics packages

  INTEGER ::                                                            &
    error_code

! Variables with intent (in/out)

  REAL ::                                                               &
    t_n(row_length, rows, model_levels),                                &
    q_n(row_length, rows, wet_model_levels),                            &
    qcl_n(row_length, rows, wet_model_levels),                          &
    qcf_n(row_length, rows, wet_model_levels),                          &
    qcf2_n(row_length, rows, wet_model_levels),                         &
!     2nd ice prog
    qrain_n(row_length, rows, wet_model_levels),                        &
!     Rain prognostic
    qgraup_n(row_length, rows, wet_model_levels),                       &
!     Graupel
    cf_n(row_length, rows, wet_model_levels),                           &
    cfl_n(row_length, rows, wet_model_levels),                          &
    cff_n(row_length, rows, wet_model_levels),                          &
    t_inc(row_length, rows, model_levels),                              &
    q_inc(row_length, rows, wet_model_levels),                          &
    qcl_inc(row_length, rows, wet_model_levels),                        &
    cf_inc(row_length, rows, wet_model_levels),                         &
    cfl_inc(row_length, rows, wet_model_levels),                        &
    sea_salt_film(salt_dim1, salt_dim2, salt_dim3),                     &
    sea_salt_jet(salt_dim1, salt_dim2, salt_dim3),                      &
    sum_eng_fluxes(row_length, rows)

  REAL, INTENT(INOUT) ::                                                &
    so4_aitken(1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
                                                 model_levels),         &
    so4_accu(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
                                                 model_levels),         &
    so4_diss(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
                                                 model_levels),         &
    soot_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
                                                 model_levels),         &
    soot_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
                                                 model_levels),         &
    soot_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
                                                 model_levels),         &
    bmass_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                                 model_levels),         &
    bmass_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                                 model_levels),         &
    bmass_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                                 model_levels),         &
    ocff_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
                                                 model_levels),         &
    ocff_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
                                                 model_levels),         &
    ocff_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
                                                 model_levels),         &
    nitr_acc(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
                                                 model_levels),         &
    nitr_diss(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                                 model_levels),         &
    dust_div1(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                                model_levels),          &
    dust_div2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                                model_levels),          &
    dust_div3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                                model_levels),          &
    dust_div4(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                                model_levels),          &
    dust_div5(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                                model_levels),          &
    dust_div6(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
                                                model_levels),          &
    aerosol(row_length, rows, model_levels)

  REAL, INTENT(IN) ::                                                   &
    biogenic(row_length, rows, model_levels)

! Variables with intent out
  REAL ::                                                               &
    photosynth_act_rad(row_length, rows),                               &
!     Net downward shortwave radiation in band 1 (w/m2).
    rad_hr(row_length, rows, 2, bl_levels),                             &
!     BL radiative (LW,SW) heating rates
    dOLR(row_length, rows),                                             &
!     TOA - surface upward LW
    sw_tile(land_field, ntiles),                                        &
!     Surface net SW on land tiles
    cos_zenith_angle(row_length,rows)
! COSP variables 
    TYPE(cosp_gridbox),INTENT(OUT) :: cosp_gbx

! local variables.


  CHARACTER (LEN=*), PARAMETER :: routinename = 'glue_rad'
  CHARACTER (LEN=240) :: cmessage

! loop counters
  INTEGER ::                                                            &
    i, j, k, info,                                                      &
    l, n, point

  INTEGER :: ij

! local variables
  INTEGER ::                                                            &
    daylight_points,                                                    &
    lit_points,                                                         &
    start_point,                                                        &
    first_point,                                                        &
    ptr_local,                                                          &
!     Pointer for LW arrays used to reduce length
!     of lines, repeating current element of first_point_local
    list_lw_points(row_length*rows),                                    &
    lw_points,                                                          &
!     Variables to enable scatter/gather of LW (like SW): creates
!     a LIST of points to do a calculation, and also
!     facilitates segmenting.
    interp_index,                                                       &
!     Variable which alternates between 0 and 1 every other
!     radiation timestep
    first_data_interp,                                                  &
!     The first data point, in data co-ords, which needs to be
!     interpolated on a PE (for use in interpolation routine)
    first_data_interp_sw,                                               &
!     Saved value to use with ISCCP diagnostics
    first_row, last_row,                                                &
    tot_daylight_points
!     Total number of daylight points in whole domain

! LW segmentation allocatable arrays
  INTEGER, ALLOCATABLE :: first_point_local(:)
  INTEGER, ALLOCATABLE :: seg_points_local(:)

! SW segmentation allocatable arrays
  INTEGER, ALLOCATABLE :: first_point_temp(:)
  INTEGER, ALLOCATABLE :: seg_points_temp(:)

! Local arrays holding information to be passed between physics
! routines.

! The following arrays are required by the PC2 cloud scheme
  REAL ::                                                               &
    t_latest(row_length,rows,model_levels),                             &
    q_latest(row_length,rows,wet_model_levels),                         &
    qcl_latest(row_length,rows,wet_model_levels),                       &
    cf_latest(row_length,rows,wet_model_levels),                        &
    cfl_latest(row_length,rows,wet_model_levels),                       &
    cff_latest(row_length,rows,wet_model_levels),                       &
    delta_t(row_length,rows,model_levels),                              &
    zeros(row_length,rows,wet_model_levels),                            &
    work_3d(row_length,rows,model_levels)

! For calculation of masses of levels - these arrays will not have
! halo information
  REAL ::                                                               &
    r_rho_levels_nh(row_length*rows,model_levels),                      &
    r_theta_levels_nh(row_length*rows,0:model_levels),                  &
    rho_r2_nh(row_length*rows,model_levels)

! For albedo scaling diagnostic:
  REAL ::                                                               &
    albobs_sc(row_length, rows, ntiles, 2)

! Radiation fields 1. SW & common with LW.

  REAL, ALLOCATABLE, SAVE :: cos_zen_rts(:,:)
  REAL, ALLOCATABLE, SAVE :: day_frac_rts(:,:)
  REAL, ALLOCATABLE, SAVE :: surfdir_rts(:,:)
  REAL, ALLOCATABLE, SAVE :: surf_down_sw_rts(:,:)
  REAL :: eps
  REAL :: sza_cor(row_length,rows)

! Fields for zenith angle correction diagnostics:
  REAL, ALLOCATABLE, SAVE :: surfsw_rts(:,:)
  REAL, ALLOCATABLE, SAVE :: toasw_rts(:,:)
  LOGICAL :: l_surfsw_cor, l_toasw_cor
  LOGICAL :: l_surfdir_cor, l_surfdif_cor
  REAL :: surfsw_cor(row_length,rows)
  REAL :: toasw_cor(row_length,rows)
  REAL :: surfdir_cor(row_length,rows)
  REAL :: surfdif_cor(row_length,rows)

  LOGICAL :: l_surfsw_rts, l_toasw_rts, l_surfdir_rts

  REAL ::                                                               &
    day_fraction(row_length,rows),                                      &
    day_fraction_2(row_length,rows,n_swcall),                           &
    cos_zenith_angle_2(row_length,rows,n_swcall),                       &
    mean_cos_zenith_angle(row_length,rows),                             &
    flandg(row_length, rows),                                           &
    land_albedo(row_length, rows, 4, n_swcall),                         &
    sea_ice_albedo(row_length, rows, 4, n_swcall),                      &
    open_sea_albedo(row_length, rows, 2,n_swcall),                      &
    netsw(row_length, rows),                                            &
!     Net short-wave absorbed by planet
    swsea(row_length, rows),                                            &
!     Net short-wave absorbed by planet
    dust_1(dust_dim1,dust_dim2),                                        &
    dust_2(dust_dim1,dust_dim2),                                        &
    dust_3(dust_dim1,dust_dim2),                                        &
    dust_4(dust_dim1,dust_dim2),                                        &
    dust_5(dust_dim1,dust_dim2),                                        &
    dust_6(dust_dim1,dust_dim2),                                        &
    local_biogenic(biogenic_dim1,biogenic_dim2),                        &
    accum_sulphate(sulp_dim1,sulp_dim2),                                &
    aitken_sulphate(sulp_dim1,sulp_dim2),                               &
    diss_sulphate(sulp_dim1,sulp_dim2),                                 &
    fresh_soot(soot_dim1,soot_dim2),aged_soot(soot_dim1,soot_dim2),     &
    fresh_bmass(bmass_dim1,bmass_dim2),                                 &
    aged_bmass(bmass_dim1,bmass_dim2),                                  &
    cloud_bmass(bmass_dim1,bmass_dim2),                                 &
    fresh_ocff(ocff_dim1, ocff_dim2),                                   &
    aged_ocff(ocff_dim1, ocff_dim2),                                    &
    cloud_ocff(ocff_dim1, ocff_dim2),                                   &
    accum_nitrate(nitrate_dim1, nitrate_dim2),                          &
    diss_nitrate(nitrate_dim1, nitrate_dim2),                           &
    height_rho(row_length, rows, model_levels),                         &
    height_theta(row_length, rows, 0:model_levels),                     &
    sw_net_land(row_length, rows),                                      &
!     SW net local flux over land
    sw_net_sice(row_length, rows)
!     SW net local flux over sea-ice

! Aerosol climatology for NWP
  REAL ::                                                               &
    local_arcl(arcl_dim1, arcl_dim2, n_arcl_compnts)
  REAL ::                                                               &
    arcl_multf
!     Mass-mixing ratio conversion factor

! UKCA_RADAER: Interaction UKCA-MODE aerosols and radiation
  REAL ::                                                               &
    local_ukca_mmr(ukca_dim1, ukca_dim2, n_ukca_cpnt),                  &
    local_ukca_cvl(ukca_dim1, ukca_dim2, n_ukca_cpnt),                  &
    local_ukca_dry(ukca_dim1, ukca_dim2, n_ukca_mode),                  &
    local_ukca_wet(ukca_dim1, ukca_dim2, n_ukca_mode),                  &
    local_ukca_rho(ukca_dim1, ukca_dim2, n_ukca_mode),                  &
    local_ukca_vol(ukca_dim1, ukca_dim2, n_ukca_mode),                  &
    local_ukca_wtv(ukca_dim1, ukca_dim2, n_ukca_mode),                  &
    local_ukca_nbr(ukca_dim1, ukca_dim2, n_ukca_mode)

  INTEGER ::                                                            &
    list_daylight_points(row_length*rows,n_swcall),                     &
    list_daylight_points_start(row_length*rows),                        &
    diag_row_list(row_length*rows,MAX(n_swcall,n_lwcall)),              &
!     List of row indices of points where diagnostics
!     are calculated.
    diag_row_list_sw(row_length*rows),                                  &
!     Saved version for ISCCP
    diag_col_list(row_length*rows,MAX(n_swcall,n_lwcall)),              &
!     List of column indices of points where diagnostics
!     are calculated.
    diag_col_list_sw(row_length*rows)
!     Saved version for ISCCP

  LOGICAL ::                                                            &
    switch(row_length,rows,MAX(n_swcall,n_lwcall))

  REAL ::                                                               &
    sindec,                                                             &
!     sin(solar declination)
    hour_angle(row_length,rows),                                        &
!     Hour Angle
    cosz_beg(row_length,rows), cosz_end(row_length,rows),               &
!     cos(zenith_angle) at beginning and end of timestep
    scs,                                                                &
!     solar constant scaling factor
    seconds_since_midnight,                                             &
    eq_time
!     The Equation of time

  INTEGER ::                                                            &
    first_point_dust_a,                                                 &
    first_point_dust_b,                                                 &
    first_point_biogenic,                                               &
    first_point_sulpc,                                                  &
    first_point_soot,                                                   &
    first_point_biomass,                                                &
    first_point_ocff,                                                   &
    first_point_seasalt,                                                &
    first_point_arcl,                                                   &
    first_point_nitrate,                                                &
    first_point_ukca,                                                   &
    nd_rad_pts,                                                         &
    n_rad_layers

  LOGICAL :: l_complete_north
!     Flag to complete field on the northern polar row
  LOGICAL :: l_complete_south
!     Flag to complete field on the southern polar row
  LOGICAL :: l_complete_deg
!     Flag to complete field because of degradation

  LOGICAL :: l_co2_3d
!     Controls use of 3D co2 field

  LOGICAL :: l_t_incr_sw
!     Flag for SW temperature increments
    
  LOGICAL :: L_cosp
! Local flag for COSP 

! SW diagnostics not on stash flag
  REAL ::                                                               &
    itoasw(row_length,rows),                                            &
    surfsw(row_length,rows)

! SW diagnostics on always true stash flag
  LOGICAL, PARAMETER :: l_flux_below_690nm_surf = .TRUE.

  REAL ::                                                               &
    flux_below_690nm_surf(row_length,rows)

! Direct Photosynthetically Active Radiation diagnostic for all tsteps
  LOGICAL :: l_direct_par
  REAL :: flxdirparsurf(row_length, rows)


! Radiation fields 2. LW

  REAL ::                                                               &
    lwsea(row_length, rows),                                            &
    olr(row_length, rows),                                              &
    surflw(row_length,rows),                                            &
    top_absorption(row_length, rows)
!     needed by energy correction
  REAL ::                                                               &
    net_atm_flux (row_length, rows)

  LOGICAL :: l_t_incr_lw
!     Flag for LW temperature increments

! needed for 8A boundary layer
  REAL ::                                                               &
    frac_tile_alb(land_field, ntype),                                   &
    alb_tile(land_field,ntiles,4),                                      &
    surf_down_sw(row_length,rows,4),                                    &
    t_rad_surf(row_length, rows),                                       &
!     Effective radiative temperature over whole grid-box
    t_rad_land(row_length, rows),                                       &
!     Effective radiative temperature over land portion of grid-box
    t_rad_sice(row_length, rows),                                       &
!     Effective radiative temperature over sea ice
    t_rad_solid(row_length, rows),                                      &
!     Radiative temperature over solid surface (now deprecated)
    tstar_tile_tmp(land_field),                                         &
!     Temporary copy of tiled surface temperature to preserve existing
!     scientific behaviour within reorganized code
    surf_emission(row_length, rows),                                    &
!     LW emitted from the surface to calculate dOLR divided by sbcon
    tile_frac(land_field,ntiles),                                       &    
    emis_tiles(ntype),                                                  &    
                                     ! Emissivity for each surface type    
    emis_here,                                                          &
!     Emissivity of the current tile in the current grid-box
    emis_nosnow,                                                        &
!     Emissivity of the current tile in the current grid-box
!     ignoring snow
    emis_land(row_length, rows)
                                     ! JULES mean land emissivity (zero    
                                     ! for ocean only grid-boxes).
  INTEGER ::                                                            &
    land_index_i(land_field),                                           &
    land_index_j(land_field),                                           &
    ssi_index_i(ssi_pts),                                               &
    ssi_index_j(ssi_pts),                                               &
    type_pts(ntype),                                                    &
!     Number of points contining each functional category
    type_index(land_field,ntype),                                       &
!     Index over functional types for aggregating albedo and emissivity
    tile_pts(ntype),                                                    &
!     Number of points contining tiles of the current category
    tile_index(land_field,ntiles)
!     Index over tiles

! Local variables required for multiple time-stepping
  LOGICAL :: l_call_swrad
  LOGICAL :: l_call_lwrad

! Local variables required for advancing the model. These are required
! for the improved timestepping algorithm.
  REAL, ALLOCATABLE, SAVE :: sw_incs_local(:,:,:,:)
  REAL, ALLOCATABLE, SAVE :: lw_incs_local(:,:,:,:)
  REAL, ALLOCATABLE, SAVE :: netsw_local(:,:,:)
  REAL, ALLOCATABLE, SAVE :: swsea_local(:,:,:)
  REAL, ALLOCATABLE, SAVE :: lwsea_local(:,:,:)
  REAL, ALLOCATABLE, SAVE :: top_abs_sw(:,:,:)
  REAL, ALLOCATABLE, SAVE :: top_abs_lw(:,:,:)
  REAL, ALLOCATABLE, SAVE :: olr_local(:,:,:)
  REAL, ALLOCATABLE, SAVE :: lw_down_local(:,:,:)
  REAL, ALLOCATABLE, SAVE :: flux_b690nm_local(:,:,:)
  REAL, ALLOCATABLE, SAVE :: surf_down_sw_local(:,:,:,:)

! Dimension of SW_diag or LW_diag required for forcing,
! improved timestepping or radiances
  INTEGER :: n_swdiag
  INTEGER :: n_lwdiag

! Looping variables for extra calls to radiation
  INTEGER :: j_sw
  INTEGER :: j_lw

! Offset in STASH item number for different sets of diagnostics
  INTEGER :: i_off


! Local variables relating to the radiance code

! Dynamically defined dimensions used by the radiance code
  INTEGER :: nd_profile
  INTEGER :: nd_field_flux_diag
  INTEGER :: nd_field_rad_diag
  INTEGER :: nd_flux_profile
  INTEGER :: nd_radiance_profile
  INTEGER :: nd_viewing_level
  INTEGER :: nd_direction
  INTEGER :: nd_column
  INTEGER :: nd_cloud_component
  INTEGER :: nd_cloud_type
  INTEGER :: nd_brdf_basis_fnc
  INTEGER :: nd_brdf_trunc
  INTEGER :: nd_channel
  INTEGER :: nd_point_tile
  INTEGER :: nd_tile
  INTEGER :: id_ct

! Generic logicals
  LOGICAL :: l_scale_inc

! Mapping of spectral bands to satellite channels
  INTEGER :: n_channel
!     Number of viewing channels
  INTEGER, POINTER :: map_channel(:) => NULL()
!     Mapping of spectral bands to satellite channels

! Satellite viewing geometry
  LOGICAL :: l_viewed
!     Flag set to .TRUE. if point can be viewed by current satellite
  INTEGER :: n_viewing_direction
  INTEGER :: n_viewing_level
  REAL :: viewing_direction(row_length*rows, 1, 2)
  REAL :: viewing_level(model_levels+1)


! Local Variables required to calculate radiative forcings.

! Variables holding the mass mixing ratios of either the
! reference call or the time advancing call.
  REAL ::                                                               &
    j_co2_mmr,                                                          &
    j_o2_mmr,                                                           &
    j_n2o_mmr,                                                          &
    j_ch4_mmr,                                                          &
    j_cfc11_mmr,                                                        &
    j_cfc12_mmr,                                                        &
    j_c113_mmr,                                                         &
    j_c114_mmr,                                                         &
    j_hcfc22_mmr,                                                       &
    j_hfc125_mmr,                                                       &
    j_hfc134_mmr

! The variable j_ozone contains the reference ozone field for the
! diagnostic call and the 'normal' ozone field for the time-advancing
! call.
  REAL :: j_ozone(row_length, rows, ozone_levels)

! Same as above for logicals defining the use of aerosols
  LOGICAL ::                                                            &
    j_l_sulpc_so2,                                                      &
    j_l_use_sulpc_direct,                                               &
    j_l_use_seasalt_direct,                                             &
    j_l_soot,                                                           &
    j_l_use_soot_direct,                                                &
    j_l_biomass,                                                        &
    j_l_use_bmass_direct,                                               &
    j_l_ocff,                                                           &
    j_l_use_ocff_direct,                                                &
    j_l_nitrate,                                                        &
    j_l_use_nitrate_direct,                                             &
    j_l_murk_rad,                                                       &
    j_l_dust,                                                           &
    j_l_use_dust,                                                       &
    j_l_use_biogenic,                                                   &
    j_l_climat_aerosol,                                                 &
    j_l_use_arcl(npd_arcl_species),                                     &
    j_l_use_ukca_radaer

  INTEGER :: j_n_arcl_species

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

! CSIGMA start
      ! Stefan-Boltzmann constant (W/m**2/K**4).
      REAL, PARAMETER ::  SBCON=5.67E-8
! CSIGMA end

! needed for vatpoles for fv_cos_theta_latitude vs cos_theta_latitude
REAL, POINTER :: xx_cos_theta_latitude (:,:)
logical, save :: first_call=.TRUE.

  IF (lhook) CALL dr_hook('GLUE_RAD',zhook_in,zhook_handle)
  eps  = EPSILON(surfsw)

IF ( l_vatpoles ) THEN
   xx_cos_theta_latitude => cos_theta_latitude
ELSE
   xx_cos_theta_latitude => fv_cos_theta_latitude
END IF ! vatpoles

!---------------------------------------------------------------------
! Section: CHECKS. Check that flags and timesteps are in agreement
!---------------------------------------------------------------------
  IF (l_forcing) THEN
    IF ( (a_sw_radstep_prog > a_sw_radstep_diag) .OR.                   &
         (a_lw_radstep_prog > a_lw_radstep_diag) ) THEN
      cmessage =                                                        &
        "If radiative forcing is required the diagnostic " //           &
        "timestep has to be larger or equal to the "   //               &
        "prognostic timestep."
      error_code=20
      GO TO 9999
    END IF
    IF (MOD(a_sw_radstep_diag, a_sw_radstep_prog) /= 0) THEN
      cmessage =                                                        &
        "The diagnostic(fast) timestep needs to be " //                 &
        "a multiple of the prognostic (slow) timestep."
      error_code=20
      GO TO 9999
    END IF
  END IF

  IF (l_timestep) THEN
    IF ( (a_sw_radstep_prog < a_sw_radstep_diag) .OR.                   &
         (a_lw_radstep_prog < a_lw_radstep_diag) ) THEN
      cmessage =                                                        &
        "If improved timestepping is required the diagnostic " //       &
        "timestep has to be smaller than the prognostic "   //          &
        "timestep."
      error_code=20
      GO TO 9999
    END IF
    IF (MOD(a_sw_radstep_prog, a_sw_radstep_diag) /= 0) THEN
      cmessage =                                                        &
        "The prognostic(slow) timestep needs to be " //                 &
        "a multiple of the diagnostic (fast) timestep."
      error_code=20
      GO TO 9999
    END IF
  END IF

! ----------------------------------------------------------------------
! Section INI. Initialisation of stash output variables.
! ----------------------------------------------------------------------

  l_co2_3d = l_co2_interactive

! Map the stashflags to the internal controlling logicals for
! diagnostics that are available on all timesteps:
  l_t_incr_sw                 = (sf_calc(181,1) .OR. sf_calc(161,1))
  l_t_incr_lw                 = (sf_calc(181,2) .OR. sf_calc(161,2))
  l_direct_par                = sf_calc(291,1)
  l_surfsw_cor                = sf_calc(202,1)
  l_toasw_cor                 = sf_calc(205,1)
  l_surfdir_cor               = sf_calc(215,1)
  l_surfdif_cor               = sf_calc(216,1)

! Set STASH flags to false if orography scheme is off.
  IF (.NOT.l_orog) sf_calc(292:296,1) = .FALSE.

! If the forcing or the improved timestepping is not required then
! the diagnostic timestep is equal to the prognostic timestep.

  IF ((.NOT.l_forcing).AND.(.NOT.l_timestep)) THEN
    l_rad_step_diag=l_rad_step_prog
  END IF
  IF (l_forcing.OR.l_timestep.OR.l_radiance) THEN
    n_swdiag=n_swcall
    n_lwdiag=n_lwcall
  ELSE
    n_swdiag=1
    n_lwdiag=1
  END IF

! Set up the logicals for diagnostics that are only available on
! radiation timesteps, contained in the structures SW_diag and LW_diag.

! DEPENDS ON: init_swdiag_logic
  CALL init_swdiag_logic(1)
! DEPENDS ON: init_lwdiag_logic
  CALL init_lwdiag_logic(1)

  IF (l_rad_step_prog) THEN

    i_off=0

! The LW logicals need to be set before the SW ones since
! the latter ones depends on the LW ones.

! DEPENDS ON: set_lwdiag_logic
    CALL set_lwdiag_logic(sf_calc,nitems,nsects, l_radiance, 1, i_off)
! DEPENDS ON: set_swdiag_logic
    CALL set_swdiag_logic(sf_calc,nitems,nsects, l_radiance, 1, i_off)

  END IF

  IF (l_radiance) THEN

    DO j_lw = 2, n_lwcall

      i_off=90+(diagnostic_offset/20)*(j_lw-1)

      CALL init_lwdiag_logic(j_lw)
      CALL set_lwdiag_logic(sf_calc,nitems,nsects, l_radiance, j_lw, i_off)

    END DO
    DO j_sw = 2, n_swcall

      i_off=90+(diagnostic_offset/20)*(j_sw-1)

      CALL init_swdiag_logic(j_sw)
      CALL set_swdiag_logic(sf_calc,nitems,nsects, l_radiance, j_sw, i_off)

    END DO
  END IF

  IF (l_forcing.OR.l_timestep) THEN

    IF (l_timestep) i_off=0
    IF (l_forcing)  i_off=diagnostic_offset

    CALL init_swdiag_logic(2)
    CALL init_lwdiag_logic(2)

    IF (l_rad_step_diag) THEN

      CALL set_lwdiag_logic(sf_calc,nitems,nsects, l_radiance, 2, i_off)
      CALL set_swdiag_logic(sf_calc,nitems,nsects, l_radiance, 2, i_off)

    END IF
  END IF


! ----------------------------------------------------------------------
! Set tile_pts, type_index and tile_index.
! ----------------------------------------------------------------------
! Type index over functional types.
! DEPENDS ON: tilepts
  CALL tilepts(land_field,frac,type_pts,type_index)
  IF (l_aggregate) THEN
    tile_pts(1) = land_field
    DO l = 1, land_field
      tile_index(l,1) = l
      tile_frac(l,1) = 1.0
    END DO
  ELSE
    DO n = 1, ntype
      tile_pts(n)=type_pts(n)
      DO j = 1, type_pts(n)
!       Here tiles match types.
        tile_index(j,n)=type_index(j,n)
        l = tile_index(j,n)
        tile_frac(l,n) = frac(l,n)
      END DO
    END DO
  END IF

! Set land_index.
  DO l = 1, land_field
    j = (land_index(l)-1)/row_length + 1
    land_index_i(l) = land_index(l) - (j-1)*row_length
    land_index_j(l) = j
  END DO

! Set ssi_index.
  DO l = 1, ssi_pts
    j = (ssi_index(l)-1)/row_length + 1
    ssi_index_i(l) = ssi_index(l) - (j-1)*row_length
    ssi_index_j(l) = j
  END DO

! ----------------------------------------------------------------------
! Set global land fraction
! ----------------------------------------------------------------------
  DO j = 1, rows
    DO i = 1, row_length
      flandg(i,j)=0.0
    END DO
  END DO
  DO l = 1, land_field
    i = land_index_i(l)
    j = land_index_j(l)
    flandg(i,j)=fland(l)
  END DO

! set CO2 array dimensions for passing to SW & LW schemes
! If used:
! CO2_DIM1 = row_length*rows (calc'd from CO2_DIM_LEN & CO2_DIM_ROW)
! CO2_DIM2 = model_levels (already passed in from CO2_DIM_LEV)
  co2_dim1 = co2_dim_len * co2_dim_row

! ----------------------------------------------------------------------
! Section RAD.0 Set up required fields and constants.
!-----------------------------------------------------------------------

  IF (error_code  ==  0) THEN
! DEPENDS ON: timer
    IF (ltimer) CALL timer ('AP1R SW Rad  ',5)

! Calculate number of seconds since midnight to the beginning of
! the timetsep.
    seconds_since_midnight = REAL( previous_time(4) * 3600              &
             + previous_time(5) * 60  + previous_time(6))

! Take a temporary copy of the tiled surface temperature if using
! aggregated surface properties. Subsequently the actual tiled
! temperature will be over-written with by t_surf to enable the
! existing scientific behaviour to be replicated by the reorganized
! code that works with tstar_tile. (Currently the temperatures may
! differ after reconfiguration, though logically they should be
! equal, and will be equal once they have been through the
! surface scheme.) tstar_tile is restored from the temporary copy
! at the end of the routine. Once an acceptable method for making
! these temperatures consistent in reconfiguration has been 
! implemented this temporary copy can be removed.
  IF (l_aggregate) tstar_tile_tmp=tstar_tile(:,1)

! Set radiation fields to zero
    DO j = 1, rows
      DO i = 1, row_length
        flux_below_690nm_surf(i,j) = 0.0
        sw_net_land(i,j) = 0.0
        sw_net_sice(i,j) = 0.0
      END DO
    END DO

    IF (l_orog) THEN
      ALLOCATE(f_orog(row_length,rows))
      f_orog=0.0
    ELSE
      ALLOCATE(f_orog(1,1))
      ALLOCATE(slope_aspect(1,1))
      ALLOCATE(slope_angle(1,1))
    END IF

    DO n = 1, ntiles
      DO l = 1, land_field
        sw_tile(l,n) = 0.
      END DO
    END DO

    DO n = 1, nice_use
      DO l = 1, ssi_pts
        sw_sice(l,n) = 0.0
      END DO
    END DO

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          t_latest(i,j,k)= t_n(i,j,k)
        END DO
      END DO
    END DO
    DO k = 1, wet_model_levels
      DO j = 1, rows
        DO i = 1, row_length
          q_latest(i,j,k)= q_n(i,j,k)
          qcl_latest(i,j,k)=qcl_n(i,j,k)
          cf_latest(i,j,k)=cf_n(i,j,k)
          cfl_latest(i,j,k)=cfl_n(i,j,k)
          cff_latest(i,j,k)=cff_n(i,j,k)
          zeros(i,j,k)=0.0
        END DO
      END DO
    END DO

    IF (l_rad_step_diag.OR.l_rad_step_prog) THEN

! Code required to initialize and update snow soot content
      IF ( l_snow_albedo ) THEN
        DO j = 1, rows
          DO i = 1, row_length
            soot(i,j) = 0.0
          END DO
        END DO
      END IF

! Code for the mixing ratio calculation of air density.
! Copy coordinate information into arrays that do not contain
! halo information
      IF (l_mixing_ratio) THEN

!       Level 0 theta level information
        DO j = 1, rows
          DO i = 1, row_length
            ij=i+(j-1)*row_length
            r_theta_levels_nh(ij,0) = r_theta_levels(i,j,0)
          END DO
        END DO

        DO k = 1, model_levels
!         Other level information
          DO j = 1, rows
            DO i = 1, row_length
              ij=i+(j-1)*row_length
              rho_r2_nh(ij,k)   = rho_r2(i,j,k)
              r_theta_levels_nh(ij,k) = r_theta_levels(i,j,k)
              r_rho_levels_nh(ij,k)   = r_rho_levels(i,j,k)
            END DO
          END DO
        END DO

      END IF  ! l_mixing_ratio

!     Code for the mineral dust aerosol scheme. We copy
!     dust aerosol into local arrays if the direct effect is
!     switched on.
      IF (l_dust .AND. l_rad_step_prog) THEN
        IF (dust_dim1  ==  rows*row_length .AND.                        &
                           dust_dim2  ==  model_levels) THEN
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length
                ij=i+(j-1)*row_length
                dust_1(ij,k) = dust_div1(i,j,k)
                dust_2(ij,k) = dust_div2(i,j,k)
                IF (.NOT. l_twobin_dust) THEN
                  dust_3(ij,k) = dust_div3(i,j,k)
                  dust_4(ij,k) = dust_div4(i,j,k)
                  dust_5(ij,k) = dust_div5(i,j,k)
                  dust_6(ij,k) = dust_div6(i,j,k)
                END IF
              END DO
            END DO
          END DO
        ELSE
          cmessage = 'DUST_DIM INCONSISTENT WITH L_DUST'
          error_code=1
          GO TO 9999
        END IF
      END IF

!     Code for the biogenic aerosol. Copy into local arrays if
!     this aerosol was requested.
      IF (l_use_biogenic .AND. l_rad_step_prog) THEN
        IF (biogenic_dim1  ==  rows*row_length .AND.                    &
            biogenic_dim2  ==  model_levels) THEN
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length
                ij=i+(j-1)*row_length
                local_biogenic(ij,k) = biogenic(i,j,k)
              END DO
            END DO
          END DO
        ELSE
          cmessage = 'BIOGENIC_DIM INCONSISTENT WITH L_USE_BIOGENIC'
          error_code=1
          GO TO 9999
        END IF
      END IF

!     Code for the Sulphur Cycle. We multiply by 4.125 to convert from
!     mass mixing ratio of sulphur atoms to mass mixing ratio of
!     ammonium sulphate.
      IF ((l_sulpc_so2 .AND. l_rad_step_prog).OR.                       &
        l_use_sulpc_indirect_sw .OR. l_use_sulpc_indirect_lw) THEN
        IF (sulp_dim1  ==  rows*row_length .AND.                        &
                             sulp_dim2  ==  model_levels)  THEN
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length
                ij=i+(j-1)*row_length
                accum_sulphate(ij,k)=so4_accu(i,j,k)*4.125
                aitken_sulphate(ij,k)=so4_aitken(i,j,k)*4.125
                diss_sulphate(ij,k)=so4_diss(i,j,k)*4.125
              END DO
            END DO
          END DO
        ELSE
          cmessage = 'SULP_DIM INCONSISTENT WITH L_SULPC'
          error_code = 1
          GO TO 9999
        END IF
      END IF

!     Code for the Soot Scheme. As for the Sulphur Cycle (above), we
!     copy soot into local arrays if the direct effect is switched on,
!     but no multiplication is required.
      IF (l_soot .AND. l_rad_step_prog) THEN
        IF (soot_dim1  ==  rows*row_length .AND.                        &
                           soot_dim2  ==  model_levels) THEN
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length
                ij=i+(j-1)*row_length
                fresh_soot(ij,k) = soot_new(i,j,k)
                aged_soot(ij,k) = soot_agd(i,j,k)
              END DO
            END DO
          END DO
        ELSE
          cmessage = 'SOOT_DIM INCONSISTENT WITH L_SOOT'
          error_code=1
          GO TO 9999
        END IF
      END IF

!     Code for the biomass aerosol scheme. As for soot, we copy biomass
!     smoke aerosol into local arrays if the direct/indirect effect is
!     switched on.
      IF ((l_biomass .AND. l_rad_step_prog)                             &
                               .OR. l_use_bmass_indirect) THEN
        IF (bmass_dim1  ==  rows*row_length .AND.                       &
                           bmass_dim2  ==  model_levels) THEN
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length
                ij=i+(j-1)*row_length
                fresh_bmass(ij,k) = bmass_new(i,j,k)
                aged_bmass(ij,k)  = bmass_agd(i,j,k)
                cloud_bmass(ij,k) = bmass_cld(i,j,k)
              END DO
            END DO
          END DO
        ELSE
          cmessage = 'BMASS_DIM INCONSISTENT WITH L_BIOMASS'
          error_code=1
          GO TO 9999
        END IF
      END IF

!     Code for the fossil-fuel organic carbon aerosol scheme. As for
!     biomass, we copy fossil-fuel organic carbon aerosol into local
!     arrays if the direct/indirect effect is switched on.
      IF ((l_ocff .AND. l_rad_step_prog) .OR. l_use_ocff_indirect) THEN
        IF (ocff_dim1  ==  rows*row_length .AND.                        &
            ocff_dim2  ==  model_levels) THEN
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length
                ij=i+(j-1)*row_length
                fresh_ocff(ij,k) = ocff_new(i,j,k)
                aged_ocff(ij,k)  = ocff_agd(i,j,k)
                cloud_ocff(ij,k) = ocff_cld(i,j,k)
              END DO
            END DO
          END DO
        ELSE
          cmessage = 'OCFF_DIM INCONSISTENT WITH L_USE_OCFF'
          error_code=1
          GO TO 9999
        END IF
      END IF

!     Code for the NWP aerosol climatology. As above, we copy the input
!     into local arrays if the climatology was requested.
      IF (n_arcl_species > 0 .AND. l_rad_step_prog) THEN

        DO l = 1, npd_arcl_compnts

          IF (i_arcl_compnts(l) /= -1) THEN
            IF (arcl_dim1 == rows * row_length .AND.                    &
                arcl_dim2 == model_levels) THEN

!             Sulphur mass-mixing ratios have to be converted
!             into ammonium sulphate.
              IF ((l == ip_arcl_sulp_ac).OR.(l == ip_arcl_sulp_ak)) THEN
                arcl_multf = 4.125
              ELSE
                arcl_multf = 1.0
              END IF

              DO k = 1, model_levels
                DO j = 1, rows
                  DO i = 1, row_length
                    ij=i+(j-1)*row_length
                    local_arcl(ij,k,i_arcl_compnts(l)) =                &
                          arcl_multf *                                  &
                          arcl(i,j,k,i_arcl_compnts(l))
                  END DO ! i
                END DO ! j
              END DO ! k
            ELSE
              cmessage = 'ARCL_DIM IS INCONSISTENT WITH L_USE_ARCL'
              error_code=1
              GO TO 9999
            END IF
          END IF

        END DO
      END IF ! n_arcl_species

!     Code for the nitrate aerosol scheme. We copy
!     nitrate aerosol into local arrays if the direct/indirect effect is
!     switched on and we multiply by 5.714 to convert from
!     mass mixing ratio of nitrogen atoms to mass mixing ratio of
!     ammonium nitrate.
      IF ((l_nitrate .AND. l_rad_step_prog)                             &
                                 .OR. l_use_nitrate_indirect) THEN
        IF (nitrate_dim1  ==  rows*row_length .AND.                     &
            nitrate_dim2  ==  model_levels) THEN
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length
                ij=i+(j-1)*row_length
                accum_nitrate(ij,k) = nitr_acc(i,j,k)*5.714
                diss_nitrate(ij,k)  = nitr_diss(i,j,k)*5.714
              END DO
            END DO
          END DO
        ELSE
          cmessage = 'NITRATE_DIM INCONSISTENT WITH L_NITRATE'
          error_code=1
          GO TO 9999
        END IF
      END IF

!     Code for UKCA_RADAER. As for all other aerosols, we copy
!     UKCA aerosols into local arrays if they are switched on.
      IF (l_ukca_radaer) THEN

        IF (ukca_dim1 == rows*row_length .AND.                          &
            ukca_dim2 == model_levels) THEN

!         Component-related variables.
          DO l = 1, n_ukca_cpnt
            DO k = 1, model_levels
              DO j = 1, rows
                DO i = 1, row_length
                  ij = i + (j-1)*row_length
                  local_ukca_mmr(ij, k, l) =                            &
                        ukca_radaer%mix_ratio(i,j,k,l)
                  local_ukca_cvl(ij, k, l) =                            &
                        ukca_radaer%comp_vol(i,j,k,l)
                END DO ! i
              END DO ! j
            END DO ! k
          END DO ! l

!         Mode-related variables.
          DO l = 1, n_ukca_mode
            DO k = 1, model_levels
              DO j = 1, rows
                DO i = 1, row_length
                  ij = i + (j-1)*row_length
                  local_ukca_dry(ij, k, l) =                            &
                        ukca_radaer%dry_diam(i,j,k,l)
                  local_ukca_wet(ij, k, l) =                            &
                        ukca_radaer%wet_diam(i,j,k,l)
                  local_ukca_rho(ij, k, l) =                            &
                        ukca_radaer%modal_rho(i,j,k,l)
                  local_ukca_vol(ij, k, l) =                            &
                        ukca_radaer%modal_vol(i,j,k,l)
                  local_ukca_wtv(ij, k, l) =                            &
                        ukca_radaer%modal_wtv(i,j,k,l)
                  local_ukca_nbr(ij, k, l) =                            &
                        ukca_radaer%modal_nbr(i,j,k,l)
                END DO ! i
              END DO ! j
            END DO ! k
          END DO ! l

        ELSE
          cmessage = 'UKCA_DIM INCONSISTENT WITH L_UKCA_RADAER'
          error_code=1
          GO TO 9999
        END IF
      END IF


      height_rho(:,:,:)=r_rho_levels(1:row_length,1:rows,:)
      height_theta(:,:,:)=r_theta_levels(1:row_length,1:rows,:)

    END IF ! l_rad_step_diag .OR. l_rad_step_prog


!-----------------------------------------------------------------------
! Section Rad.0.1 Allocate SW variables and diagnostics
!-----------------------------------------------------------------------

    IF (l_rad_step_prog) THEN

! DEPENDS ON: allocate_sw_diag
      CALL allocate_sw_diag( row_length, rows, model_levels,            &
        cloud_levels, ntiles, 1 )

      IF (l_orog) ALLOCATE(orog_corr(row_length,rows))
      IF (l_rad_szacor.OR.l_rad_perturb) THEN
        ALLOCATE(cos_zen_rts(row_length, rows))
        ALLOCATE(day_frac_rts(row_length, rows))
      END IF

      l_surfdir_rts=l_rad_szacor.OR.l_surfdir_cor.OR.l_surfdif_cor
      l_surfsw_rts=l_surfsw_cor.OR.l_toasw_cor
      l_toasw_rts=l_toasw_cor
! Check to see if allocatable variables may be needed for the
! corrected flux diagnostics on non-radiation timesteps:
      DO i = 1, totitems
        IF (stlist(s_modl,i)==1.AND.stlist(s_sect,i)==1) THEN
          IF (stlist(s_item,i)==215.OR.stlist(s_item,i)==216) THEN
            l_surfdir_rts=.TRUE.
          ELSE IF (stlist(s_item,i)==202) THEN
            l_surfsw_rts=.TRUE.
          ELSE IF (stlist(s_item,i)==205) THEN
            l_surfsw_rts=.TRUE.
            l_toasw_rts=.TRUE.
          END IF
        END IF
      END DO
      ALLOCATE(surf_down_sw_rts(row_length, rows))
      IF (l_surfdir_rts) THEN
        ALLOCATE(surfdir_rts(row_length, rows))
      END IF
      IF (l_surfsw_rts) THEN
        ALLOCATE(surfsw_rts(row_length, rows))
      END IF
      IF (l_toasw_rts) THEN
        ALLOCATE(toasw_rts(row_length, rows))
      END IF

      IF (l_timestep) THEN
        ALLOCATE(sw_incs_local(row_length, rows,                        &
                                    0:model_levels+1,n_swcall))
        ALLOCATE(surf_down_sw_local(row_length,rows,4,n_swcall))
        ALLOCATE(netsw_local(row_length, rows, n_swcall))
        ALLOCATE(swsea_local(row_length, rows ,n_swcall))
        ALLOCATE(top_abs_sw(row_length, rows, n_swcall))
        ALLOCATE(flux_b690nm_local(row_length,rows,n_swcall))
      END IF

      IF (l_radiance) THEN
        DO j_sw = 2, n_swcall
          CALL allocate_sw_diag( row_length, rows,                      &
            model_levels, cloud_levels, ntiles, j_sw )
        END DO
      END IF

    END IF

    IF (l_rad_step_diag) THEN
      IF (l_timestep.OR.l_forcing) THEN
        CALL allocate_sw_diag( row_length, rows, model_levels,          &
          cloud_levels, ntiles, 2 )
      END IF
    END IF

! Calculates sine of the solar declination and the scaling
! factor for solar intensity from the day number and year.
! DEPENDS ON: solpos
    CALL solpos (previous_time(7), previous_time(1),                    &
      l_cal360, l_sec_var, l_eqt, eq_time, sindec, scs)

!------------------------------------------------------------
! Section Rad.0.2 Do preliminary timestep calculations
!------------------------------------------------------------
    IF (l_rad_step_prog.OR.l_rad_step_diag) THEN

      DO j_sw = 1, n_swcall

        l_call_swrad=.FALSE.

        IF (l_timestep) THEN

! Here we need to be careful with the timestepping. In the
! case of the improved timestepping the prognostic (slow)
! call is done with the 'astronomy' calculated on the
! prognostic timestep and the diagnostic (fast) call is done
! with the astronomy calculated on the diganostic timestep.
! This ensures that all the averaging over different timesteps
! is done accurately.

          IF ((j_sw==1).AND.(l_rad_step_prog)) THEN
            radiation_tstep = radiation_tstep_prog
            frac_tile_alb=frac
            l_call_swrad=.TRUE.
          ELSE IF ((j_sw==2).AND.(l_rad_step_diag)) THEN
            radiation_tstep = radiation_tstep_diag
            frac_tile_alb=frac
            l_call_swrad=.TRUE.
          END IF

        ELSE IF (l_forcing) THEN

! In the case of the forcing the astronomy is always calculated
! on the prognostic timescale, i.e. we always use the prognostic
! timestep to calculate the astronomy.

          IF ((j_sw==1).AND.(l_rad_step_prog)) THEN
            radiation_tstep = radiation_tstep_prog
            frac_tile_alb=frac
            l_call_swrad=.TRUE.
          END IF
          IF ((j_sw==2).AND.(l_rad_step_diag)) THEN
            radiation_tstep = radiation_tstep_prog
            IF (c2c_land_s(j_sw-1)) THEN
              frac_tile_alb=frac_control
            ELSE
              frac_tile_alb=frac
            END IF
            l_call_swrad=.TRUE.
          END IF

        ELSE
          radiation_tstep = radiation_tstep_prog
          frac_tile_alb=frac
          l_call_swrad=.TRUE.
        END IF

        IF (l_call_swrad) THEN

          IF (l_rad_perturb.AND.(j_sw==2)) THEN

            cos_zenith_angle_2(:,:,j_sw)=cos_zen_rts
            day_fraction_2(:,:,j_sw)=day_frac_rts

          ELSE

            IF (l_orog) ALLOCATE(sol_bearing(row_length,rows))

! DEPENDS ON: solang
            CALL solang(                                                &
! input constants
              sindec, seconds_since_midnight,                           &
              radiation_tstep, eq_time,                                 &
! row and column dependent constants
              true_latitude,                                            &
              true_longitude,                                           &
! size variables
              row_length*rows,                                          &
! output fields
              day_fraction_2(1,1,j_sw) ,                                &
              cos_zenith_angle_2(1,1,j_sw),                             &
              hour_angle, cosz_beg, cosz_end )

! Set rounding-error size values to zero - the criterion depends
! on the frequency of full SW calculations because on the physics
! timesteps which are not SW timesteps a test has to be done to
! avoid using the unset data for such points.
            DO j = 1, rows
              DO i = 1, row_length
                IF ( cos_zenith_angle_2(i,j,j_sw)*                      &
                                   day_fraction_2(i,j,j_sw) <           &
                     (1.e-10*timestep/radiation_tstep) ) THEN
                  cos_zenith_angle_2(i,j,j_sw) = 0.0
                  day_fraction_2(i,j,j_sw) = 0.0
                END IF
              END DO
            END DO

! Transfer solar bearing diagnostic before it is deallocated:
            IF (sw_diag(j_sw)%l_sol_bearing) THEN
              sw_diag(j_sw)%sol_bearing = sol_bearing
            END IF

! Calculate orography correction:
            IF (l_orog) THEN
! DEPENDS ON: solinc
              CALL solinc(row_length, rows,                             &
                cos_zenith_angle_2(:,:,j_sw),                           &
                cosz_beg, cosz_end)
              DEALLOCATE(sol_bearing)
              IF (l_rad_szacor) THEN
                WHERE (ABS(orog_corr-0.5) < eps)
                  orog_corr=0.5+eps
                END WHERE
              END IF
            END IF

            IF (l_rad_szacor.OR.l_rad_perturb) THEN
              cos_zen_rts=cos_zenith_angle_2(:,:,j_sw)
              day_frac_rts=day_fraction_2(:,:,j_sw)
            END IF

          END IF ! l_rad_perturb.and.(j_sw==2)

! DEPENDS ON: ftsa
          CALL ftsa (                                                   &
! input fields
            flandg, ice_fract, ice_fract_cat, t_surf, tstar_sice_cat,   &
            cos_zenith_angle_2(1,1,j_sw),                               &
            snow_depth, snow_depth_sea_cat, ice_thick_cat,              &
! max and min sea ice albedo specifications
            sw_alpham, sw_alphac, sw_alphab, sw_dtice,                  &
            l_ssice_albedo, l_mod_barker_albedo,                        &
            l_sice_meltponds, l_sice_scattering, l_sice_hadgem1a,       &
            l_cice_alb,                                                 &
            dt_bare, dalb_bare_wet, pen_rad_frac, sw_beta,              &
! size and control variables
            row_length*rows, row_length*rows, nice_use,                 &
! output arguments
            sea_ice_albedo(1,1,1,j_sw),                                 &
            open_sea_albedo(1,1,1,j_sw) )


!-----------------------------------------------------------------------
! Calculate tile albedos. ntype is taken from a module within the
! routine.
!-----------------------------------------------------------------------
! DEPENDS ON: tile_albedo
          CALL tile_albedo ( row_length*rows,                           &
            land_field,land_index,ntiles,type_pts,                      &
            type_index,                                                 &
            l_aggregate,                                                &
            l_snow_albedo,albsoil,albobs_sw,albobs_vis,albobs_nir,      &
            cos_zenith_angle_2(1,1,j_sw),frac_tile_alb,                 &
            lai,rgrain, snow_tile,soot,tstar_tile,                      &
            z0_tile, alb_tile,                                          &
            land_albedo(1,1,1,j_sw),albobs_sc,can_rad_mod )
! Write the albedo scalings to sw_diag, if requested:
          IF (SW_diag(j_sw)%l_vis_albedo_sc) THEN
            DO k =1, ntiles
              DO j = 1, rows
                DO i = 1, row_length
                  sw_diag(j_sw)%vis_albedo_sc(i,j,k) = albobs_sc(i,j,k,1)
                END DO
              END DO
            END DO
          END IF
          IF (SW_diag(j_sw)%l_nir_albedo_sc) THEN
            DO k =1, ntiles
              DO j = 1, rows
                DO i = 1, row_length
                  sw_diag(j_sw)%nir_albedo_sc(i,j,k) = albobs_sc(i,j,k,2)
                END DO
              END DO
            END DO
          END IF

        END IF ! if sw radiation call required
      END DO ! loop over radiation calls

      !CABLE
      !tile_frac = 0.

      DO j_sw = n_swcall, 1, -1
! Set COSP flag only on prognostic radiation steps
        IF (L_cosp_in.AND.(j_sw == 1).AND.L_rad_step_prog) THEN
          L_cosp=.TRUE.
        ELSE
          L_cosp=.FALSE.
        END IF

! For each SW call set the correct mass mixing ratios and switch the
! aerosols on/off.
! Note that in the case of the aerosols the reference state is
! 'no aerosols' and hence they are switched off in the diagnostic call.

! Since the "default" for each variable is to have it set in the
! diagnostic call to the same as in the prognostic, we set them all
! unconditionally to the prognostic & then re-set as needed.
        j_co2_mmr = co2_mmr
        j_o2_mmr  = o2mmr
        j_n2o_mmr = n2o_mix_ratio
        j_ch4_mmr = ch4_mix_ratio
        j_ozone = ozone  ! Note this is an array assignment
        j_l_sulpc_so2          = l_sulpc_so2
        j_l_use_sulpc_direct   = l_use_sulpc_direct
        j_l_use_seasalt_direct = l_use_seasalt_direct
        j_l_soot               = l_soot
        j_l_use_soot_direct    = l_use_soot_direct
        j_l_biomass            = l_biomass
        j_l_use_bmass_direct   = l_use_bmass_direct
        j_l_ocff               = l_ocff
        j_l_use_ocff_direct    = l_use_ocff_direct
        j_l_nitrate            = l_nitrate
        j_l_use_nitrate_direct = l_use_nitrate_direct
        j_l_dust               = l_dust
        j_l_use_dust           = l_use_dust
        j_l_use_biogenic       = l_use_biogenic
        j_l_climat_aerosol     = l_climat_aerosol
        j_l_murk_rad           = l_murk_rad
        j_l_use_arcl           = l_use_arcl
        j_n_arcl_species       = n_arcl_species
        j_l_use_ukca_radaer    = l_ukca_radaer

! Turn aerosols off for the "cloud only" call in order to minimise
! computational cost. Clear-sky and other diagnostics have already
! been turned off in set_swdiag_logic.
        IF (l_rad_perturb.AND.(j_sw==2)) THEN
          j_l_sulpc_so2          = .FALSE.
          j_l_use_sulpc_direct   = .FALSE.
          j_l_use_seasalt_direct = .FALSE.
          j_l_soot               = .FALSE.
          j_l_use_soot_direct    = .FALSE.
          j_l_biomass            = .FALSE.
          j_l_use_bmass_direct   = .FALSE.
          j_l_ocff               = .FALSE.
          j_l_use_ocff_direct    = .FALSE.
          j_l_nitrate            = .FALSE.
          j_l_use_nitrate_direct = .FALSE.
          j_l_dust               = .FALSE.
          j_l_use_dust           = .FALSE.
          j_l_use_biogenic       = .FALSE.
          j_l_climat_aerosol     = .FALSE.
          j_l_murk_rad           = .FALSE.
          j_l_use_arcl           = .FALSE.
          j_n_arcl_species       = 0
          j_l_use_ukca_radaer    = .FALSE.
        END IF

        l_call_swrad=.FALSE.

! if Improved timestepping required
        IF (l_timestep) THEN
          IF ((j_sw==1).AND.(l_rad_step_prog)) THEN
            a_sw_radstep=a_sw_radstep_prog
            l_call_swrad=.TRUE.
          END IF
          IF ((j_sw==2).AND.(l_rad_step_diag)) THEN
            IF (l_rad_perturb.AND.l_rad_step_prog) THEN
              a_sw_radstep=a_sw_radstep_prog
            ELSE
              a_sw_radstep=a_sw_radstep_diag
            END IF
            l_call_swrad=.TRUE.
          END IF

! if Forcing required
        ELSE IF (l_forcing) THEN

! Decide if it is necessary to call radiation
          IF ((j_sw==1).AND.(l_rad_step_prog)) THEN
            a_sw_radstep=a_sw_radstep_prog
            l_call_swrad=.TRUE.
          END IF
          IF ((j_sw==2).AND.(l_rad_step_diag)) THEN
            a_sw_radstep=a_sw_radstep_prog
            l_call_swrad=.TRUE.
          END IF

          IF ( j_sw > 1 ) THEN ! Diagnostic Call

! Gases
            IF ( c2c_co2(j_sw-1) ) j_co2_mmr = co2_mmr_d
            IF ( c2c_o2(j_sw-1) )  j_o2_mmr  = o2_mmr_d
            IF ( c2c_n2o(j_sw-1) ) j_n2o_mmr = n2o_mix_ratio_d
            IF ( c2c_ch4(j_sw-1) ) j_ch4_mmr = ch4_mix_ratio_d

! Ozone: Note that aerosol contains the MURK array.  It is in fact what
! is called d1(MURK) in Simon Tett's old vn4.5 mod.  Care needs to be
! taken so that MURK has the correct dimensions before it arrives in
! RAD_CTL2.
            IF ( c2c_o3(j_sw-1) ) j_ozone = aerosol

! Aerosols
            IF ( c2c_sulpc_d(j_sw-1) ) j_l_sulpc_so2          = .FALSE.
            IF ( c2c_sulpc_d(j_sw-1) ) j_l_use_sulpc_direct   = .FALSE.
            IF ( c2c_seas_d(j_sw-1) )  j_l_use_seasalt_direct = .FALSE.
            IF ( c2c_soot_d(j_sw-1) )  j_l_soot               = .FALSE.
            IF ( c2c_soot_d(j_sw-1) )  j_l_use_soot_direct    = .FALSE.
            IF ( c2c_bmb_d(j_sw-1) )   j_l_biomass            = .FALSE.
            IF ( c2c_bmb_d(j_sw-1) )   j_l_use_bmass_direct   = .FALSE.
            IF ( c2c_ocff_d(j_sw-1) )  j_l_ocff               = .FALSE.
            IF ( c2c_ocff_d(j_sw-1) )  j_l_use_ocff_direct    = .FALSE.
            IF ( c2c_nitr_d(j_sw-1) )  j_l_nitrate            = .FALSE.
            IF ( c2c_nitr_d(j_sw-1) )  j_l_use_nitrate_direct = .FALSE.
            IF ( c2c_dust_d(j_sw-1) )  j_l_dust               = .FALSE.
            IF ( c2c_dust_d(j_sw-1) )  j_l_use_dust           = .FALSE.
            IF ( c2c_biog_d(j_sw-1) )  j_l_use_biogenic       = .FALSE.
            IF ( c2c_ukca_d(j_sw-1) )  j_l_use_ukca_radaer    = .FALSE.

          END IF
        ELSE
          a_sw_radstep=a_sw_radstep_prog
          l_call_swrad=.TRUE.
        END IF

! If call to SW Radiation required
        IF (l_call_swrad) THEN

! DEPENDS ON: prelim_swrad
          CALL prelim_swrad(error_code,                                 &
            at_extremity, n_proc,                                       &
! Model Dimensions
            row_length,rows,model_levels,                               &
! Model Switches
            model_domain,mt_global,l_rad_deg,                           &
            sw_control(j_sw)%l_subsample,                               &
            sw_control(j_sw)%l_geostationary,                           &
! Time stepping Information
            timestep_number,a_sw_radstep,                               &
! Satellite Geometry
            sw_control(j_sw)%min_view_lon,                              &
            sw_control(j_sw)%max_view_lon,                              &
            sw_control(j_sw)%min_view_lat,                              &
            sw_control(j_sw)%max_view_lat,                              &
! Number of Call
            j_sw,                                                       &
! Other variables
            true_latitude,true_longitude,                               &
            seconds_since_midnight,                                     &
            tot_daylight_points,daylight_points,                        &
            day_fraction_2(1,1,j_sw),                                   &
            list_daylight_points(1,j_sw), rad_mask,                     &
            switch(1,1,j_sw), first_data_interp,                        &
            first_data_interp_sw,                                       &
            diag_row_list(1,j_sw), diag_row_list_sw,                    &
            diag_col_list(1,j_sw), diag_col_list_sw )

! Zero the calculated outputs: this will later simplify the
! treatment of the case when there are no lit points.
          sw_incs(:,:,:)=0.0
          netsw (:,:) = 0.0
          swsea (:,:) = 0.0
          IF (sw_control(1)%l_extra_top) THEN
            top_absorption(:,:) = 0.0
          END IF
          surf_down_sw (:,:,:) = 0.0

! Sunlit points for COSP
          IF (L_COSP) THEN
            DO j = 1, rows
              DO i = 1, row_length
                IF (day_fraction_2(i,j,j_sw) >  0.) &
                  cosp_gbx%sunlit((j-1)*row_length + i) = 1.0
              END DO
            END DO
          END IF

          IF ( daylight_points  >   0 ) THEN

!Parallelise over segments.
!$OMP  PARALLEL DEFAULT(SHARED)                                              &
!$OMP& PRIVATE(i, j, lit_points,start_point,first_point,first_point_dust_a,  &
!$OMP& first_point_dust_b,                                                   &
!$OMP& first_point_sulpc,first_point_soot,first_point_biomass,               &
!$OMP& first_point_seasalt,nd_rad_pts,n_rad_layers,error_code,               &
!$OMP& first_point_biogenic, first_point_ocff, first_point_arcl,             &
!$OMP& first_point_nitrate,l_scale_inc, map_channel,nd_cloud_type,           &
!$OMP& nd_cloud_component, nd_column, id_ct, nd_channel, nd_viewing_level,   &
!$OMP& nd_direction, nd_profile, nd_radiance_profile, nd_flux_profile,       &
!$OMP& nd_field_flux_diag, nd_field_rad_diag, nd_brdf_basis_fnc,             &
!$OMP& nd_brdf_trunc, nd_tile, nd_point_tile, n_viewing_level, viewing_level,&
!$OMP& n_viewing_direction, viewing_direction, n_channel, first_point_ukca,  &
!$OMP& meta_segments, segments, first_point_temp, seg_points_temp, ipar)

!$OMP SINGLE
!Determine the number of threads in this parallel region
num_parallel_sw=1
!$ num_parallel_sw = omp_get_num_threads()
!$OMP END SINGLE
!Implicit barrier

!Determine the team member number
!NB: this is numbered from 1 to follow the Fortran convention.
ipar=1
!$ ipar = omp_get_thread_num()+1

            !Set up meta-information for segments
            CALL segments_mod_seg_meta(meta_segments, ipar, num_parallel_sw, &
                    daylight_points, a_sw_seg_size, a_sw_segments)

            !Allocate storage
            ALLOCATE( segments        ( meta_segments%num_segments ) )
            ALLOCATE( first_point_temp( meta_segments%num_segments ) )
            ALLOCATE(  seg_points_temp( meta_segments%num_segments ) )

            !Find specific points in each segment
            CALL segments_mod_segments(segments,meta_segments,  &
              daylight_points, row_length, rows,                &
              list_points=list_daylight_points(:,j_sw))

            first_point_temp(:) = segments(:)%fp
            seg_points_temp(:)  = segments(:)%seg_points

            DO i = 1, meta_segments%num_segments
              lit_points  = segments(i)%use_points
              start_point = segments(i)%start_index
              first_point = segments(i)%fp

              !Additional work on list_points. Note that although work is done
              !on this array, each thread and each segment should only affect
              !one element.
              DO j = start_point, segments(i)%end_index
                list_daylight_points_start(j) = list_daylight_points(j,j_sw) &
                                               -segments(i)%start+1
              END DO

! Set the first point of the dust arrays to be used.
              IF (l_dust) THEN
                IF(l_twobin_dust) THEN
                  first_point_dust_a = first_point
                  first_point_dust_b = 1
                ELSE
                  first_point_dust_a = first_point
                  first_point_dust_b = first_point
                END IF
              ELSE
                first_point_dust_a = 1
                first_point_dust_b = 1
              END IF

! Set the first point of the biogenic array.
              IF (l_use_biogenic) THEN
                first_point_biogenic=first_point
              ELSE
                first_point_biogenic=1
              END IF

! Set the first point of the array of sulphate to be used.
! A separate assignment is necessary since this array will
! not be of the full size unless the sulphur cycle is on.
              IF (l_sulpc_so2 .OR. l_use_sulpc_indirect_sw) THEN
                first_point_sulpc=first_point
              ELSE
                first_point_sulpc=1
              END IF

              IF (l_soot) THEN
                first_point_soot=first_point
              ELSE
                first_point_soot=1
              END IF

              IF (l_biomass .OR. l_use_bmass_indirect) THEN
                first_point_biomass=first_point
              ELSE
                first_point_biomass=1
              END IF

              IF (l_ocff .OR. l_use_ocff_indirect) THEN
                first_point_ocff=first_point
              ELSE
                first_point_ocff=1
              END IF

              IF (l_use_seasalt_indirect .OR. l_use_seasalt_direct) THEN
                first_point_seasalt=first_point
              ELSE
                first_point_seasalt=1
              END IF

              IF (n_arcl_species > 0) THEN
                first_point_arcl = first_point
              ELSE
                first_point_arcl = 1
              END IF

              IF (l_nitrate .OR. l_use_nitrate_indirect) THEN
                first_point_nitrate=first_point
              ELSE
                first_point_nitrate=1
              END IF

              IF (l_ukca_radaer) THEN
                first_point_ukca=first_point
              ELSE
                first_point_ukca=1
              END IF

! Set the actual size of arrays in the radiation code:
! for some architectures (such as that of Cray vector
! machines) on odd size is preferred to avoid memory
! bank conflicts.
              nd_rad_pts=2*(lit_points/2)+1

! Set the number of layers seen in the radiation code.
! This may optionally be 1 greater than the number used
! in the rest of the model to avoid spurious effects
! resulting from the upper boundary (principally in
! stratospheric configurations).
              IF (sw_control(1)%l_extra_top) THEN
                n_rad_layers=model_levels+1
              ELSE
                n_rad_layers=model_levels
              END IF


! Set dynamically determined array sizes required for
! radiation calculations. There is redundancy here
! because the viewing direction is calculated for the
! whole PE within each segement.

              ALLOCATE(map_channel(sw_spectrum(j_sw)%n_band))

! DEPENDS ON: r2_set_rad_dim
              CALL r2_set_rad_dim(row_length, rows,                     &
                true_latitude, true_longitude,                          &
                .TRUE., sindec, seconds_since_midnight,                 &
                sw_control(j_sw)%l_geostationary,                       &
                sw_control(j_sw)%sat_hgt,                               &
                sw_control(j_sw)%sat_lat,                               &
                sw_control(j_sw)%sat_lon,                               &
                sw_control(j_sw)%i_cloud,                               &
                sw_control(j_sw)%i_angular_integration,                 &
                sw_control(j_sw)%i_sph_mode,                            &
                row_length*rows, seg_points_temp(i),                    &
                model_levels, global_cloud_top,                         &
                nd_cloud_type, nd_cloud_component,                      &
                sw_control(j_sw)%l_extra_top,                           &
                id_ct, n_rad_layers, nd_column,                         &
                sw_spectrum(j_sw)%npd_band,                             &
                sw_spectrum(j_sw)%n_band,                               &
                map_channel, n_channel, nd_channel,                     &
                nd_viewing_level, n_viewing_level,                      &
                viewing_level,                                          &
                nd_direction, n_viewing_direction,                      &
                viewing_direction,                                      &
                nd_brdf_basis_fnc, nd_brdf_trunc,                       &
                nd_profile, nd_flux_profile,                            &
                nd_radiance_profile,                                    &
                nd_field_flux_diag, nd_field_rad_diag,                  &
                l_ctile, nd_tile,                                       &
                nd_point_tile)

! To remain general heating rates are always
! converted to increments across a timestep:
              l_scale_inc = .TRUE.


! DEPENDS ON: r2_swrad3z
              CALL r2_swrad3z(error_code,                               &
! Mixing Ratios
  q_n(first_point,1,1), j_co2_mmr,j_ozone(first_point,1,1),             &
  j_o2_mmr, co2_dim1, co2_dim2, co2_3d(first_point,1,1),                &
  l_co2_3d,                                                             &
  j_n2o_mmr, j_ch4_mmr,                                                 &

! Pressures and Temperatures
  p_star(first_point,1), p_layer_boundaries(first_point,1,0),           &
  p_layer_centres(first_point,1,0),t_n(first_point,1,1),                &

! Options for COSP 
  L_cosp,                                                               &

! Options for treating clouds
  global_cloud_top,                                                     &
  l_inhom_cloud, inhom_cloud_sw, dp_corr_strat, dp_corr_conv,           &

! Stratiform Cloud Fields
  l_pc2, area_cloud_fraction(first_point,1,1), cf_n(first_point,1,1),   &
  qcl_n(first_point,1,1), qcf_n(first_point,1,1),                       &
  n_drop_pot(first_point,1,1),                                          &

! Convective Cloud Fields
  cca(first_point,1,1), cclwp(first_point,1),                           &
  ccw(first_point,1,1), lcbase(first_point,1),                          &
  ccb(first_point,1), cct(first_point,1),                               &

! Surface Fields
  land_albedo(first_point,1,1,j_sw), l_ctile, l_use_spec_sea,           &
  flandg(first_point,1), sea_ice_albedo(first_point,1,1,j_sw),          &
  open_sea_albedo(first_point,1,1,j_sw),                                &
  ice_fract(first_point,1), land_sea_mask(first_point,1),               &
  land0p5(first_point,1), snow_depth(first_point,1),                    &

! Solar Fields
  cos_zenith_angle_2(first_point,1,j_sw),                               &
  day_fraction_2(first_point,1,j_sw),                                   &
  list_daylight_points_start(start_point), scs,                         &

! Aerosol Fields
  j_l_climat_aerosol, l_clim_aero_hgt, l_hadgem1_clim_aero,             &
  zh(first_point,1), aero_bl_levels,                                    &
  j_l_dust, j_l_use_dust, dust_dim1, dust_dim2,                         &
  dust_1(first_point_dust_a,1), dust_2(first_point_dust_a,1),           &
  dust_3(first_point_dust_b,1), dust_4(first_point_dust_b,1),           &
  dust_5(first_point_dust_b,1), dust_6(first_point_dust_b,1),           &
  j_l_use_biogenic, biogenic_dim1, biogenic_dim2,                       &
  local_biogenic(first_point_biogenic, 1),                              &
  j_l_sulpc_so2, j_l_use_sulpc_direct, l_use_sulpc_indirect_sw,         &
  sulp_dim1, sulp_dim2,                                                 &
  accum_sulphate(first_point_sulpc, 1),                                 &
  aitken_sulphate(first_point_sulpc, 1),                                &
  diss_sulphate(first_point_sulpc, 1),                                  &
  sea_salt_film(first_point_seasalt,1,1),                               &
  sea_salt_jet(first_point_seasalt,1,1),                                &
  l_use_seasalt_indirect, j_l_use_seasalt_direct,                       &
  salt_dim1*salt_dim2, salt_dim3, j_l_soot, j_l_use_soot_direct,        &
  soot_dim1, soot_dim2, fresh_soot(first_point_soot, 1),                &
  aged_soot(first_point_soot, 1), j_l_biomass, j_l_use_bmass_direct,    &
  bmass_dim1, bmass_dim2, fresh_bmass(first_point_biomass, 1),          &
  aged_bmass(first_point_biomass, 1),                                   &
  cloud_bmass(first_point_biomass, 1), l_use_bmass_indirect,            &
  j_l_ocff, j_l_use_ocff_direct, ocff_dim1, ocff_dim2,                  &
  fresh_ocff(first_point_ocff, 1),aged_ocff(first_point_ocff, 1),       &
  cloud_ocff(first_point_ocff, 1), l_use_ocff_indirect,                 &
  j_l_nitrate, j_l_use_nitrate_direct, nitrate_dim1, nitrate_dim2,      &
  accum_nitrate(first_point_nitrate, 1),                                &
  diss_nitrate(first_point_nitrate, 1), l_use_nitrate_indirect,         &
  j_l_use_arcl, arcl_dim1, arcl_dim2, j_n_arcl_species,                 &
  n_arcl_compnts,i_arcl_compnts,local_arcl(first_point_arcl,1,1),       &
  aerosol(first_point,1,1), j_l_murk_rad, ntot_land, ntot_sea,          &
  j_l_use_ukca_radaer, ukca_radaer, ukca_dim1, ukca_dim2,               &
  local_ukca_mmr(first_point_ukca, 1, 1),                               &
  local_ukca_cvl(first_point_ukca, 1, 1),                               &
  local_ukca_dry(first_point_ukca, 1, 1),                               &
  local_ukca_wet(first_point_ukca, 1, 1),                               &
  local_ukca_rho(first_point_ukca, 1, 1),                               &
  local_ukca_vol(first_point_ukca, 1, 1),                               &
  local_ukca_wtv(first_point_ukca, 1, 1),                               &
  local_ukca_nbr(first_point_ukca, 1, 1),                               &

! time
  previous_time,                                                        &

! grid-dependent arrays
  true_latitude(first_point,1),                                         &

! Level of tropopause
  trindx(first_point,1),                                                &

! Spectrum
  sw_spectrum(j_sw),                                                    &

! Algorithmic options
  sw_control(j_sw), timestep, l_mod_k_flux, l_scale_inc,                &

! Satellite viewing geometry
  n_viewing_direction, viewing_direction(first_point,1, 1),             &
  viewing_direction(first_point,1,2),n_viewing_level,viewing_level,     &

! All diagnostics and associated arrays
  n_channel, map_channel, sw_diag(j_sw),                                &
  diag_row_list(start_point,j_sw),diag_col_list(start_point,j_sw),      &

! Dimensions
  lit_points,seg_points_temp(i),model_levels,n_rad_layers,              &
  global_cloud_top,wet_model_levels,ozone_levels,row_length,            &
  rows,rows*row_length,nd_field_flux_diag,                              &
  nd_field_rad_diag,nd_profile,n_rad_layers,nd_column,                  &
  n_cca_levels,nd_channel, nd_flux_profile,                             &
  nd_radiance_profile, nd_viewing_level, nd_direction,                  &
  nd_cloud_component, nd_cloud_type,                                    &
  nd_brdf_basis_fnc, nd_brdf_trunc,                                     &
  nd_point_tile, nd_tile, id_ct, n_ukca_mode, n_ukca_cpnt,              &

! Output data
  surf_down_sw(first_point,1,1),                                        &
  flux_below_690nm_surf(first_point,1),                                 &
  l_flux_below_690nm_surf,                                              &
  netsw(first_point,1),                                                 &
  top_absorption(first_point,1),                                        &
  swsea(first_point,1),                                                 &
  sw_incs(first_point,1,0),                                             &

! Variables needed to calculate layer masses
  rho_r2_nh(first_point,1),r_rho_levels_nh(first_point,1),              &
  r_theta_levels_nh(first_point,0),                                     &
  q_n(first_point,1,1), qcl_n(first_point,1,1),                         &
  qcf_n(first_point,1,1), qcf2_n(first_point,1,1),                      &
  qrain_n(first_point,1,1), qgraup_n(first_point,1,1),                  &
  l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio,                &

! COSP input arguments
  cosp_gbx)

              DEALLOCATE(map_channel)

            END DO              ! Loop over segments

   ! Deallocate the segmentation arrays
            DEALLOCATE( first_point_temp )
            DEALLOCATE(  seg_points_temp )
            DEALLOCATE(  segments )

!$OMP END PARALLEL

          END IF                 ! If daylight_points>0


! Radiative fluxes may not have been calculated at all
! points: we now fill in as required. In the case of
! spatial degradation calls to SWAPBOUNDS are made and
! these must be made on all PEs, which requires that
! this code shall executed even when the PE contains
! no daylit points.

        IF (.NOT. l_vatpoles) THEN
! At the North Pole in a global domain calculations
! are performed only at the first point if lit, so
! the rest of the row must be filled in if this point
! is lit.
          l_complete_north=(model_domain==mt_global).AND.               &
                       ( at_extremity(pnorth) ) .AND.                   &
                       ( switch(1,rows,j_sw) )

! At the South Pole in a global domain calculations
! are performed only at the first point if lit, so
! the rest of the row must be filled in if this point
! is lit.
          l_complete_south=(model_domain==mt_global).AND.               &
                       ( at_extremity(psouth) ) .AND.                   &
                       ( switch(1,1,j_sw) )
        ELSE
          l_complete_north=.FALSE.
          l_complete_south=.FALSE.
        END IF ! vatpoles

! When spatial degradation is performed fields must be
! filled in at alternate points.
          l_complete_deg = ( l_rad_deg ) .AND.                          &
                     ( tot_daylight_points > 0 )

! Set addressing limits for spatial degradation.
          IF ( l_complete_deg ) THEN
            first_row=1
            last_row=rows
            IF (.NOT. l_vatpoles) THEN
            IF (model_domain  ==  mt_global) THEN
              IF (at_extremity(pnorth)) THEN
                last_row=rows-1
              END IF
              IF (at_extremity(psouth)) THEN
                first_row=2
              END IF
            END IF
            END IF ! vatpoles
          END IF

! Call appropriate subroutines to fill in missing data
! as required.
          IF ( l_complete_north .OR. l_complete_south .OR.              &
               l_complete_deg ) THEN

! DEPENDS ON: fill_missing_data_sw
            CALL fill_missing_data_sw(                                  &
              off_x, off_y, row_length, rows, model_levels, ntiles,     &
              salt_dim1, salt_dim2, salt_dim3,                          &
              cloud_levels,first_row,last_row,                          &
              first_data_interp, es_space_interp,                       &
              l_complete_north, l_complete_south,                       &
              l_complete_deg,l_flux_below_690nm_surf,                   &
              n_channel, j_sw,                                          &
              sw_control(1)%l_extra_top,                                &
              sw_incs,                                                  &
              netsw,                                                    &
              swsea,                                                    &
              flux_below_690nm_surf,                                    &
              top_absorption,                                           &
              surf_down_sw,                                             &
              sea_salt_film, sea_salt_jet)
          END IF
          IF (l_rad_perturb.AND.(j_sw==1)) THEN
            sw_incs_local(:,:,:,1) = sw_incs                            &
                                   - sw_incs_local(:,:,:,2)
            netsw_local(:,:,1) = netsw                                  &
                               - netsw_local(:,:,2)
            swsea_local(:,:,1) = swsea                                  &
                               - swsea_local(:,:,2)
            surf_down_sw_local(:,:,:,1) = surf_down_sw                  &
                                     - surf_down_sw_local(:,:,:,2)
            flux_b690nm_local(:,:,1) = flux_below_690nm_surf            &
                                     - flux_b690nm_local(:,:,2)
            IF (sw_control(1)%l_extra_top) THEN
              top_abs_sw(:,:,1) = top_absorption                        &
                                - top_abs_sw(:,:,2)
            END IF
          ELSE IF (l_timestep) THEN
            sw_incs_local(:,:,:,j_sw)=sw_incs
            netsw_local(:,:,j_sw)=netsw
            swsea_local(:,:,j_sw)=swsea
            surf_down_sw_local(:,:,:,j_sw)=surf_down_sw
            flux_b690nm_local(:,:,j_sw)=flux_below_690nm_surf
            IF (sw_control(1)%l_extra_top) THEN
              top_abs_sw(:,:,j_sw)=top_absorption
            END IF
          END IF

        END IF ! l_call_swrad
      END DO ! radiation calls

      IF (l_rad_perturb.AND.l_rad_step_prog) THEN

!       For this case the full fluxes are already set.

        mean_cos_zenith_angle(:,:) =                                    &
           cos_zenith_angle_2(:,:,2)*day_fraction_2(:,:,2)

      ELSE IF (l_timestep.AND.l_rad_step_diag) THEN
        sw_incs(:,:,:) = sw_incs_local(:,:,:,1)                         &
                       + sw_incs_local(:,:,:,2)
        netsw(:,:)     = netsw_local(:,:,1)                             &
                       + netsw_local(:,:,2)
        swsea(:,:)     = swsea_local(:,:,1)                             &
                       + swsea_local(:,:,2)
        surf_down_sw(:,:,:) = surf_down_sw_local(:,:,:,1)               &
                            + surf_down_sw_local(:,:,:,2)
        flux_below_690nm_surf(:,:) = flux_b690nm_local(:,:,1)           &
                                   + flux_b690nm_local(:,:,2)
        IF (sw_control(1)%l_extra_top) THEN
          top_absorption(:,:) = top_abs_sw(:,:,1)                       &
                              + top_abs_sw(:,:,2)
        END IF
        mean_cos_zenith_angle(:,:) =                                    &
           cos_zenith_angle_2(:,:,2)*day_fraction_2(:,:,2)
      ELSE IF (l_rad_step_prog) THEN
        mean_cos_zenith_angle(:,:) =                                    &
           cos_zenith_angle_2(:,:,1)*day_fraction_2(:,:,1)
      END IF

      IF (l_rad_perturb.AND.(.NOT.l_rad_step_prog)) THEN
        WHERE (sw_incs(:,:,0) < 0.0)
          netsw(:,:) = netsw(:,:) -                                     &
            sw_incs(:,:,0)*mean_cos_zenith_angle(:,:)
          sw_incs(:,:,0) = 0.0
        END WHERE
        WHERE (swsea < 0.0)
          netsw = netsw - (1.-flandg)*swsea
          swsea = 0.0
        END WHERE
        IF (l_rm_neg_par) THEN
          WHERE (surf_down_sw < 0.0) surf_down_sw = 0.0
          WHERE (flux_below_690nm_surf < 0.0) flux_below_690nm_surf=0.0
        END IF
      END IF
      !print *, "jhan:_control5 tile_pts ",tile_pts
!if (.NOT. first_call) then
      call cable_glue_rad_init( surf_down_sw )
      call cable_control5( alb_tile, land_albedo,         &
                  TILE_PTS, TILE_INDEX, surf_down_sw )        
!endif
first_call=.FALSE.

!     Calculate net surface SW for diagnostic
      IF (l_ctile) THEN
        DO j = 1, rows
          DO i = 1, row_length
            surfsw(i,j) = sw_incs(i,j,0)                                &
                        * mean_cos_zenith_angle(i,j)                    &
                        + (1.-flandg(i,j))*swsea(i,j)
            IF (flandg(i,j) == 1.0) swsea(i,j)=rmdi
          END DO
        END DO
      ELSE
        DO j = 1, rows
          DO i = 1, row_length
            surfsw(i,j) = sw_incs(i,j,0)                                &
                        * mean_cos_zenith_angle(i,j)                    &
                        + swsea(i,j)
            IF (land_sea_mask(i,j)) swsea(i,j)=rmdi
          END DO
        END DO
      END IF

!     Set photosynthetically active radiation from flux below 690nm
      sw_incs(:,:,model_levels+1) = surf_down_sw(:,:,1)                 &
                                  + surf_down_sw(:,:,2)
      IF (l_direct_par) THEN
        dirpar_inc = surf_down_sw(:,:,1)
      END IF

      IF (l_rad_perturb.AND.(.NOT.l_rad_step_prog).AND.                 &
                            (.NOT.l_rm_neg_par) ) THEN
!       Zero the array after PAR values are set to reproduce old behaviour
        WHERE (surf_down_sw < 0.0) surf_down_sw = 0.0
      END IF

!     Calculate total downwards SW (sum of direct/diffuse vis/near-ir)
      surf_down_sw_rts=SUM(surf_down_sw,DIM=3)

      IF (ALLOCATED(surfdir_rts)) THEN
        surfdir_rts=surf_down_sw(:,:,1)+surf_down_sw(:,:,3)
      END IF
      IF (ALLOCATED(surfsw_rts)) THEN
        WHERE (mean_cos_zenith_angle > eps)
          surfsw_rts=surfsw/mean_cos_zenith_angle
        ELSEWHERE
          surfsw_rts=0.0
        END WHERE
      END IF
      IF (ALLOCATED(toasw_rts)) THEN
        WHERE (mean_cos_zenith_angle > eps)
          toasw_rts=scs * sc - netsw/mean_cos_zenith_angle
        ELSEWHERE
          toasw_rts=0.0
        END WHERE
      END IF


!     Calculate net surface SW on land tiles
      DO n = 1, ntiles
        DO l = 1, land_field
          sw_tile_rts(l,n) = 0.
        END DO
      END DO
      DO n = 1, ntiles
        DO point = 1, tile_pts(n)
          l = tile_index(point,n)
          i = land_index_i(l)
          j = land_index_j(l)
          DO k = 1, 4
            sw_tile_rts(l,n) = sw_tile_rts(l,n)                         &
               + (1. - alb_tile(l,n,k))*surf_down_sw(i,j,k)
          END DO
        END DO
      END DO

!     Calculate net surface SW on sea ice categories
      DO n = 1, nice_use
        DO l = 1, ssi_pts
          sw_sice_rts(l,n) = 0.0
        END DO
      END DO
      DO n = 1, nice_use
        DO point = 1, sice_pts_ncat(n)
          l = sice_index_ncat(point,n)
          i = ssi_index_i(l)
          j = ssi_index_j(l)
          DO k = 1, 4
            sw_sice_rts(l,n) = sw_sice_rts(l,n)                         &
              + (1.0 - alb_sice(l,n,k))*surf_down_sw(i,j,k)
          END DO
        END DO
      END DO

!     Set indexing for albedo arrays appropriate to timestep
      IF (l_rad_step_prog) THEN
        j_sw = 1
      ELSE
        j_sw = 2
      END IF

!     Calculate mean land albedo
      land_alb = 0.0
      DO l = 1, land_field
        i = land_index_i(l)
        j = land_index_j(l)
        IF (surf_down_sw_rts(i,j) > eps) THEN
          DO k = 1, 4
            land_alb(i,j) = land_alb(i,j)                               &
              + land_albedo(i,j,k,j_sw)*surf_down_sw(i,j,k)
          END DO
          land_alb(i,j) = land_alb(i,j)/surf_down_sw_rts(i,j)
        END IF
      END DO

!     Calculate mean albedo of lake tile
      IF (l_flake_model.AND.(.NOT.l_aggregate)) THEN
!     First zero the array
        lake_albedo(:) = 0.0
        n = lake
        DO point = 1, tile_pts(n)
          l = tile_index(point,n)
          i = land_index_i(l)
          j = land_index_j(l)
          IF (surf_down_sw_rts(i,j) > eps) THEN
            lake_albedo(l) = 1.0-(sw_tile_rts(l,n)/surf_down_sw_rts(i,j))
          END IF
        END DO
      END IF

!     Calculate mean sea ice albedo
      sice_alb = 0.0
      DO point = 1, sice_pts
        l = sice_index(point)
        i = ssi_index_i(l)
        j = ssi_index_j(l)
        IF (surf_down_sw_rts(i,j) > eps) THEN
          DO k = 1, 4
            sice_alb(i,j) = sice_alb(i,j)                               &
              + sea_ice_albedo(i,j,k,j_sw)*surf_down_sw(i,j,k)
          END DO
          sice_alb(i,j) = sice_alb(i,j)/surf_down_sw_rts(i,j)
        END IF
      END DO

    END IF ! on a radiation timestep


! Calculate day fraction and mean cos(solar zenith angle while
! the sun is up) for each grid point for this physics timestep:
! (if in fact full SW calculations are being done every timestep, this
! is of course unnecessary, as are various calculations later on)

    CALL solang(                                                        &
! input constants
      sindec, seconds_since_midnight,                                   &
      timestep, eq_time,                                                &
! row and column dependent constants
      true_latitude,                                                    &
      true_longitude,                                                   &
! size variables
      row_length*rows,                                                  &
! output fields
      day_fraction, cos_zenith_angle,                                   &
      hour_angle, cosz_beg, cosz_end)

    IF (l_rad_szacor) THEN
      WHERE (cos_zen_rts*day_frac_rts < eps)
        cos_zenith_angle = 0.0
        day_fraction = 0.0
      END WHERE
      IF (l_orog) THEN
        WHERE ((surfdir_rts > eps).AND.(surfdir_rts < orog_corr*scs*sc) &
          .AND.(cos_zen_rts > eps).AND.(cos_zenith_angle > eps)         &
          .AND.(orog_corr > SQRT(eps)))
          sza_cor = (orog_corr - 0.5) * (scs*sc                         &
                    * EXP(LOG(surfdir_rts/(orog_corr*scs*sc))           &
                    *(cos_zen_rts/cos_zenith_angle))                    &
                    - surfdir_rts/orog_corr)                            &
                    /surf_down_sw_rts
        ELSEWHERE
          sza_cor = 0.0
        END WHERE
      ELSE
        WHERE ((surfdir_rts > eps).AND.(surfdir_rts < scs*sc)           &
          .AND.(cos_zen_rts > eps).AND.(cos_zenith_angle > eps))
          sza_cor = 0.5*(scs*sc * EXP(LOG(surfdir_rts/(scs*sc))         &
                 *(cos_zen_rts/cos_zenith_angle)) - surfdir_rts)        &
                 /surf_down_sw_rts
        ELSEWHERE
          sza_cor = 0.0
        END WHERE
      END IF
    ELSE
      sza_cor = 0.0
    END IF

! Combine the two terms to give the mean cos zenith angle over the
! whole of the physics timestep.
! calculate incoming SW at top of atmosphere
    DO j = 1, rows
      DO i = 1, row_length
        cos_zenith_angle(i,j) = cos_zenith_angle(i,j) *                 &
                                day_fraction(i,j)
        itoasw(i,j) = cos_zenith_angle(i,j) * scs * sc
      END DO
    END DO

! calculate corrected net surface SW for diagnostic:
    IF ((l_surfsw_cor) .AND. ALLOCATED(surfsw_rts)) THEN
      surfsw_cor = (1.0+sza_cor)*cos_zenith_angle*surfsw_rts
    END IF
! calculate corrected TOA outgoing SW for diagnostic:
    IF ((l_toasw_cor) .AND. ALLOCATED(toasw_rts)                        &
                      .AND. ALLOCATED(surfsw_rts)) THEN
      toasw_cor=cos_zenith_angle*(toasw_rts-sza_cor*surfsw_rts)
    END IF

! calculate corrected direct surface SW for diagnostic:
    IF ((l_surfdir_cor.OR.l_surfdif_cor)                                &
         .AND. ALLOCATED(surfdir_rts)) THEN
      IF (l_orog.AND.l_rad_szacor) THEN
        surfdir_cor = cos_zenith_angle*(surfdir_rts +                   &
          orog_corr * sza_cor * surf_down_sw_rts/(orog_corr-0.5) )
      ELSE
        surfdir_cor = cos_zenith_angle*(surfdir_rts +                   &
                        2.0 * sza_cor * surf_down_sw_rts)
      END IF
      WHERE (surfdir_cor < 0.0) surfdir_cor = 0.0
    END IF
! calculate corrected diffuse surface SW for diagnostic:
    IF ((l_surfdif_cor) .AND. ALLOCATED(surfdir_rts)) THEN
      surfdif_cor = (1.0+sza_cor)*cos_zenith_angle*surf_down_sw_rts     &
                      - surfdir_cor
    END IF


! Is the PC2 cloud scheme being used?
    IF (l_pc2) THEN

!     Reset _latest values to _n values
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            delta_t(i,j,k) = sw_incs(i,j,k) * cos_zenith_angle(i,j)
            t_latest(i,j,k)   = t_n(i,j,k)
          END DO
        END DO
      END DO
      DO k = 1, wet_model_levels
        DO j = 1, rows
          DO i = 1, row_length
            q_latest(i,j,k)   = q_n(i,j,k)
            qcl_latest(i,j,k) = qcl_n(i,j,k)
            cf_latest(i,j,k)  = cf_n(i,j,k)
            cfl_latest(i,j,k) = cfl_n(i,j,k)
          END DO
        END DO
      END DO

! ----------------------------------------------------------------------
! Homogeneous forcing. Note the temperature increment from the shortwave
! is added in this routine
! ----------------------------------------------------------------------

! DEPENDS ON: pc2_homog_plus_turb
      CALL pc2_homog_plus_turb(p_layer_centres(1,1,1),                  &
        wet_model_levels,                                               &
        timestep, t_latest, cf_latest, cfl_latest,                      &
        cff_latest, q_latest, qcl_latest, delta_t(1,1,1),               &
        zeros, zeros, zeros, 0.0, 0.0,                                  &
        l_mixing_ratio)

! Add increments from the homogeneous forcing to the increment variables
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            t_inc(i,j,k) = t_inc(i,j,k) + t_latest(i,j,k)-t_n(i,j,k)
          END DO
        END DO
      END DO
      DO k = 1, wet_model_levels
        DO j = 1, rows
          DO i = 1, row_length
            q_inc(i,j,k) = q_inc(i,j,k) + q_latest(i,j,k)-q_n(i,j,k)
            qcl_inc(i,j,k) = qcl_inc(i,j,k)                             &
                             + qcl_latest(i,j,k)-qcl_n(i,j,k)
            cf_inc(i,j,k) = cf_inc(i,j,k)                               &
                             + cf_latest(i,j,k)-cf_n(i,j,k)
            cfl_inc(i,j,k) = cfl_inc(i,j,k)                             &
                             + cfl_latest(i,j,k)-cfl_n(i,j,k)
          END DO
        END DO
      END DO

    ELSE  ! L_pc2

! add SW radiative heating to temperatures
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            delta_t(i,j,k) = sw_incs(i,j,k) * cos_zenith_angle(i,j)
            t_inc(i,j,k) = t_inc(i,j,k) + delta_t(i,j,k)
            t_latest(i,j,k) = t_n(i,j,k) + delta_t(i,j,k)
          END DO
        END DO
      END DO

    END IF  ! L_pc2

! Get T_incr for output as STASH diagnostic
    IF ( ( l_t_incr_sw )                                                &
      .OR. ( sf_calc(232,1) )                                           &
      ) THEN
!     Increments will be calculated if to be diagnosed directly
!     or to be used to determine heating rates.
      ALLOCATE ( t_incr_diagnostic(row_length,rows,model_levels) )
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            t_incr_diagnostic(i,j,k) =  delta_t(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k

    ELSE
      ALLOCATE ( t_incr_diagnostic(1,1,1) )
    END IF                    ! on STASHflag


!   Set up photosynthetically active surface radiation, if calculated
    IF (l_flux_below_690nm_surf) THEN
      DO j = 1, rows
        DO i = 1, row_length
          photosynth_act_rad(i,j) = sw_incs(i,j,model_levels+1)         &
                                    * cos_zenith_angle(i,j)
          IF (l_direct_par) THEN
            flxdirparsurf(i,j) = dirpar_inc(i,j)                        &
                                   * cos_zenith_angle(i,j)
          END IF
        END DO
      END DO
      IF (l_rad_szacor) THEN
        IF (l_direct_par) THEN
          IF (l_orog) THEN
            flxdirparsurf = flxdirparsurf + orog_corr * sza_cor         &
                              * photosynth_act_rad/(orog_corr-0.5)
          ELSE
            flxdirparsurf = flxdirparsurf +                             &
                              2.0 * sza_cor * photosynth_act_rad
          END IF
        END IF
        photosynth_act_rad = photosynth_act_rad * (1.0+sza_cor)
      END IF
    END IF

!   Set up net surface SW on land tiles
    DO n = 1, ntiles
      DO point = 1, tile_pts(n)
        l = tile_index(point,n)
        i = land_index_i(l)
        j = land_index_j(l)
        sw_tile(l,n) = sw_tile_rts(l,n)                                 &
           * cos_zenith_angle(i,j) * (1.0+sza_cor(i,j))
        sw_net_land(i,j) = sw_net_land(i,j)                             &
           + sw_tile(l,n) * tile_frac(l,n)
      END DO
    END DO

!   Set up net surface SW on sea ice categories.

    DO n = 1, nice_use
      DO point = 1, sice_pts_ncat(n)
        l = sice_index_ncat(point,n)
        i = ssi_index_i(l)
        j = ssi_index_j(l)
        sw_sice(l,n) = sw_sice_rts(l,n)                                 &
           * cos_zenith_angle(i,j) * (1.0+sza_cor(i,j))
        sw_net_sice(i,j) = sw_net_sice(i,j)                             &
           + sw_sice(l,n) * ice_fract_cat(i,j,n)
      END DO
    END DO

! DEPENDS ON: timer
    IF (ltimer) CALL timer ('AP1R SW Rad  ',6)

! ----------------------------------------------------------------------
! Section RAD.1.1 Short Wave Radiation Energy correction code
!-----------------------------------------------------------------------

    IF ((l_rad_step_prog.OR.l_rad_step_diag) .AND. l_emcorr) THEN

! Sum short wave fluxes into the atmosphere and
! add into the net diabatic fluxes into the
! atmosphere for use in the energy correction
! procedure

      IF (sw_control(1)%l_extra_top) THEN
!       The energy absorbed above the top of the model in
!       the radiation scheme does not contribute to the
!       energy absorbed, but the diagnostics are calculated
!       at the top of the atmosphere, so the net atmospheric
!       flux must be adjusted.
        DO j = 1, rows
          DO i = 1, row_length
            net_atm_flux(i,j) = netsw(i,j) - surfsw(i,j)                &
                              - top_absorption(i,j)
          END DO
        END DO
      ELSE
        DO j = 1, rows
          DO i = 1, row_length
            net_atm_flux(i,j) = netsw(i,j) - surfsw(i,j)
          END DO
        END DO
      END IF

      IF (l_orog) net_atm_flux = net_atm_flux + f_orog


! DEPENDS ON: flux_diag
      CALL flux_diag(net_atm_flux, xx_cos_theta_latitude,               &
        row_length, rows ,off_x,off_y, 1.0,                             &
        sum_eng_fluxes,radiation_tstep)

    END IF

! ----------------------------------------------------------------------
! Section RAD.1.2 Short Wave Radiation diagnostics
!-----------------------------------------------------------------------

    DO j_sw = 1, n_swdiag

! Check that sw diagnostics requested this timestep

      IF (error_code  ==  0                                             &
        .AND. sf_calc(0,1)                                              &
        ) THEN

        IF (l_timestep) THEN

          i_off=0

          IF (j_sw == n_swdiag) THEN

            IF (l_rad_perturb.AND.l_rad_step_prog) THEN

              IF (sw_diag(2)%l_flux_up)                                 &
                  sw_diag(1)%flux_up =                                  &
                  sw_diag(1)%flux_up -                                  &
                  sw_diag(2)%flux_up

              IF (sw_diag(2)%l_flux_down)                               &
                  sw_diag(1)%flux_down =                                &
                  sw_diag(1)%flux_down -                                &
                  sw_diag(2)%flux_down

              IF (sw_diag(2)%l_solar_out_toa)                           &
                  sw_diag(1)%solar_out_toa =                            &
                  sw_diag(1)%solar_out_toa -                            &
                  sw_diag(2)%solar_out_toa

              IF (sw_diag(2)%l_surface_down_flux)                       &
                  sw_diag(1)%surface_down_flux =                        &
                  sw_diag(1)%surface_down_flux -                        &
                  sw_diag(2)%surface_down_flux

              IF (sw_diag(2)%l_net_flux_trop)                           &
                  sw_diag(1)%net_flux_trop =                            &
                  sw_diag(1)%net_flux_trop -                            &
                  sw_diag(2)%net_flux_trop

              IF (sw_diag(2)%l_up_flux_trop)                            &
                  sw_diag(1)%up_flux_trop =                             &
                  sw_diag(1)%up_flux_trop -                             &
                  sw_diag(2)%up_flux_trop

              IF (sw_diag(2)%l_flux_direct)                             &
                  sw_diag(1)%flux_direct =                              &
                  sw_diag(1)%flux_direct -                              &
                  sw_diag(2)%flux_direct

              IF (sw_diag(2)%l_flux_diffuse)                            &
                  sw_diag(1)%flux_diffuse =                             &
                  sw_diag(1)%flux_diffuse -                             &
                  sw_diag(2)%flux_diffuse

              IF (sw_diag(2)%l_FlxSolBelow690nmSurf)                    &
                  sw_diag(1)%FlxSolBelow690nmSurf =                     &
                  sw_diag(1)%FlxSolBelow690nmSurf -                     &
                  sw_diag(2)%FlxSolBelow690nmSurf

              IF (sw_diag(2)%l_FlxSeaBelow690nmSurf)                    &
                  sw_diag(1)%FlxSeaBelow690nmSurf =                     &
                  sw_diag(1)%FlxSeaBelow690nmSurf -                     &
                  sw_diag(2)%FlxSeaBelow690nmSurf

            END IF

            IF (l_rad_step_diag) THEN

              IF (sw_diag(2)%l_flux_up)                                 &
                  sw_diag(2)%flux_up =                                  &
                  sw_diag(2)%flux_up +                                  &
                  sw_diag(1)%flux_up

              IF (sw_diag(2)%l_flux_down)                               &
                  sw_diag(2)%flux_down =                                &
                  sw_diag(2)%flux_down +                                &
                  sw_diag(1)%flux_down

              IF (sw_diag(2)%l_solar_out_toa)                           &
                  sw_diag(2)%solar_out_toa =                            &
                  sw_diag(2)%solar_out_toa +                            &
                  sw_diag(1)%solar_out_toa

              IF (sw_diag(2)%l_surface_down_flux)                       &
                  sw_diag(2)%surface_down_flux =                        &
                  sw_diag(2)%surface_down_flux +                        &
                  sw_diag(1)%surface_down_flux

              IF (sw_diag(2)%l_net_flux_trop)                           &
                  sw_diag(2)%net_flux_trop =                            &
                  sw_diag(2)%net_flux_trop +                            &
                  sw_diag(1)%net_flux_trop

              IF (sw_diag(2)%l_up_flux_trop)                            &
                  sw_diag(2)%up_flux_trop =                             &
                  sw_diag(2)%up_flux_trop +                             &
                  sw_diag(1)%up_flux_trop

              IF (sw_diag(2)%l_flux_direct)                             &
                  sw_diag(2)%flux_direct =                              &
                  sw_diag(2)%flux_direct +                              &
                  sw_diag(1)%flux_direct

              IF (sw_diag(2)%l_flux_diffuse)                            &
                  sw_diag(2)%flux_diffuse =                             &
                  sw_diag(2)%flux_diffuse +                             &
                  sw_diag(1)%flux_diffuse

              IF (sw_diag(2)%l_FlxSolBelow690nmSurf)                    &
                  sw_diag(2)%FlxSolBelow690nmSurf =                     &
                  sw_diag(2)%FlxSolBelow690nmSurf +                     &
                  sw_diag(1)%FlxSolBelow690nmSurf

              IF (sw_diag(2)%l_FlxSeaBelow690nmSurf)                    &
                  sw_diag(2)%FlxSeaBelow690nmSurf =                     &
                  sw_diag(2)%FlxSeaBelow690nmSurf +                     &
                  sw_diag(1)%FlxSeaBelow690nmSurf

            END IF ! l_rad_step_diag
          END IF ! j_sw == n_swdiag


        ELSE IF (l_forcing.AND.l_rad_step_diag) THEN

! Calculate SW forcings. Note that forcings are only calculated for
! fluxes. Difference between j_sw=1 (call with forcing in place)
! and j_sw= 2 (call with reference values). N.B. Call to radiation
! schemes done first with j_sw=2 then with j_sw=1

          i_off=diagnostic_offset*(j_sw-1)

          IF (j_sw == n_swdiag) THEN

            IF (sw_diag(j_sw)%l_flux_up)                                &
                sw_diag(j_sw)%flux_up =                                 &
                sw_diag(1)%flux_up -                                    &
                sw_diag(j_sw)%flux_up

            IF (sw_diag(j_sw)%l_flux_down)                              &
                sw_diag(j_sw)%flux_down =                               &
                sw_diag(1)%flux_down -                                  &
                sw_diag(j_sw)%flux_down

            IF (sw_diag(j_sw)%l_flux_up_clear)                          &
                sw_diag(j_sw)%flux_up_clear =                           &
                sw_diag(1)%flux_up_clear -                              &
                sw_diag(j_sw)%flux_up_clear

            IF (sw_diag(j_sw)%l_flux_down_clear)                        &
                sw_diag(j_sw)%flux_down_clear =                         &
                sw_diag(1)%flux_down_clear -                            &
                sw_diag(j_sw)%flux_down_clear

            IF (sw_diag(j_sw)%l_solar_out_toa)                          &
                sw_diag(j_sw)%solar_out_toa =                           &
                sw_diag(1)%solar_out_toa -                              &
                sw_diag(j_sw)%solar_out_toa

            IF (sw_diag(j_sw)%l_solar_out_clear)                        &
                sw_diag(j_sw)%solar_out_clear =                         &
                sw_diag(1)%solar_out_clear -                            &
                sw_diag(j_sw)%solar_out_clear

            IF (sw_diag(j_sw)%l_surface_down_flux)                      &
                sw_diag(j_sw)%surface_down_flux =                       &
                sw_diag(1)%surface_down_flux -                          &
                sw_diag(j_sw)%surface_down_flux

            IF (sw_diag(j_sw)%l_surf_down_clr)                          &
                sw_diag(j_sw)%surf_down_clr =                           &
                sw_diag(1)%surf_down_clr -                              &
                sw_diag(j_sw)%surf_down_clr

            IF (sw_diag(j_sw)%l_surf_up_clr)                            &
                sw_diag(j_sw)%surf_up_clr =                             &
                sw_diag(1)%surf_up_clr -                                &
                sw_diag(j_sw)%surf_up_clr

            IF (sw_diag(j_sw)%l_clear_hr)                               &
                sw_diag(j_sw)%clear_hr =                                &
                sw_diag(1)%clear_hr -                                   &
                sw_diag(j_sw)%clear_hr

            IF (sw_diag(j_sw)%l_net_flux_trop)                          &
                sw_diag(j_sw)%net_flux_trop =                           &
                sw_diag(1)%net_flux_trop -                              &
                sw_diag(j_sw)%net_flux_trop

            IF (sw_diag(j_sw)%l_up_flux_trop)                           &
                sw_diag(j_sw)%up_flux_trop =                            &
                sw_diag(1)%up_flux_trop -                               &
                sw_diag(j_sw)%up_flux_trop

            IF (sw_diag(j_sw)%l_flux_direct)                            &
                sw_diag(j_sw)%flux_direct =                             &
                sw_diag(1)%flux_direct -                                &
                sw_diag(j_sw)%flux_direct

            IF (sw_diag(j_sw)%l_flux_diffuse)                           &
                sw_diag(j_sw)%flux_diffuse =                            &
                sw_diag(1)%flux_diffuse -                               &
                sw_diag(j_sw)%flux_diffuse

            IF (sw_diag(j_sw)%l_uvflux_direct)                          &
                sw_diag(j_sw)%uvflux_direct =                           &
                sw_diag(1)%uvflux_direct -                              &
                sw_diag(j_sw)%uvflux_direct

            IF (sw_diag(j_sw)%l_uvflux_up)                              &
                sw_diag(j_sw)%uvflux_up =                               &
                sw_diag(1)%uvflux_up -                                  &
                sw_diag(j_sw)%uvflux_up

            IF (sw_diag(j_sw)%l_uvflux_net)                             &
                sw_diag(j_sw)%uvflux_net =                              &
                sw_diag(1)%uvflux_net -                                 &
                sw_diag(j_sw)%uvflux_net

            IF (sw_diag(j_sw)%l_surf_uv)                                &
                sw_diag(j_sw)%surf_uv =                                 &
                sw_diag(1)%surf_uv -                                    &
                sw_diag(j_sw)%surf_uv

            IF (sw_diag(j_sw)%l_surf_uv_clr)                            &
                sw_diag(j_sw)%surf_uv_clr =                             &
                sw_diag(1)%surf_uv_clr -                                &
                sw_diag(j_sw)%surf_uv_clr

          END IF
        ELSE IF (l_radiance) THEN
          IF (j_sw == 1) THEN
            i_off = 0
          ELSE
            i_off=90+(diagnostic_offset/20)*(j_sw-1)
          END IF
        ELSE
          i_off=0
        END IF

! DEPENDS ON: diagnostics_sw
        CALL diagnostics_sw(                                            &
          row_length, rows, model_levels,                               &
          wet_model_levels, cloud_levels, ntiles,                       &
          n_rows, global_row_length, global_rows,                       &
          halo_i, halo_j, off_x, off_y, me,                             &
          n_proc, n_procx, n_procy,                                     &
          g_rows, g_row_length,                                         &
          at_extremity,                                                 &
          timestep,i_off,                                               &
          t_n, t_inc,                                                   &
          q_n, qcl_n, cf_n, cfl_n,                                      &
          t_latest, q_latest, qcl_latest,                               &
          cf_latest, cfl_latest,                                        &
          surfsw, itoasw, surfsw_cor, toasw_cor,                        &
          surfdir_cor, surfdif_cor,                                     &
          swsea, flux_below_690nm_surf,                                 &
          photosynth_act_rad, flxdirparsurf,                            &
          f_orog, slope_aspect, slope_angle,                            &
          sw_net_land,sw_net_sice,                                      &
          t_incr_diagnostic,                                            &
          sw_spectrum(j_sw)%n_band,                                     &
          sea_salt_film, sea_salt_jet,                                  &
          salt_dim1, salt_dim2, salt_dim3,                              &
          j_sw,                                                         &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
          stashwork1)

      END IF
    END DO


    DEALLOCATE ( t_incr_diagnostic )
    IF (l_orog) THEN
      DEALLOCATE(f_orog)
    ELSE
      DEALLOCATE(f_orog,slope_aspect,slope_angle)
    END IF

    IF (MOD(timestep_number,a_sw_radstep_prog) == 0) THEN
      IF (l_rad_szacor.OR.l_rad_perturb) THEN
        DEALLOCATE(cos_zen_rts)
        DEALLOCATE(day_frac_rts)
      END IF
      IF (ALLOCATED(surfdir_rts     )) DEALLOCATE(surfdir_rts)
      IF (ALLOCATED(surf_down_sw_rts)) DEALLOCATE(surf_down_sw_rts)
      IF (ALLOCATED(surfsw_rts      )) DEALLOCATE(surfsw_rts)
      IF (ALLOCATED(toasw_rts       )) DEALLOCATE(toasw_rts)
      IF (l_orog) THEN
        DEALLOCATE(orog_corr)
      END IF
    END IF

    IF (l_timestep) THEN
      IF (l_rad_step_diag) THEN
! DEPENDS ON: deallocate_swdiag
        CALL deallocate_swdiag(2)
      END IF
      IF (timestep_number >  1) THEN
        IF (MOD(timestep_number,a_sw_radstep_prog) == 0) THEN

          DEALLOCATE(sw_incs_local)
          DEALLOCATE(swsea_local)
          DEALLOCATE(netsw_local)
          DEALLOCATE(top_abs_sw)
          DEALLOCATE(flux_b690nm_local)
          DEALLOCATE(surf_down_sw_local)

          CALL deallocate_swdiag(1)

        END IF
      END IF
    ELSE IF (l_forcing) THEN
      IF (l_rad_step_prog) THEN
        CALL deallocate_swdiag(1)
      END IF
      IF (l_rad_step_diag) THEN
        CALL deallocate_swdiag(2)
      END IF
    ELSE IF (l_radiance) THEN
      IF (l_rad_step_prog) THEN
        DO j_sw = 1, n_swcall
          CALL deallocate_swdiag(j_sw)
        END DO
      END IF
    ELSE
      IF (l_rad_step_prog) THEN
        CALL deallocate_swdiag(1)
      END IF
    END IF

! ----------------------------------------------------------------------
! Section RAD.2 Long Wave Radiation Code.
!-----------------------------------------------------------------------

    IF (l_rad_step_prog) THEN

! DEPENDS ON: allocate_lw_diag
      CALL allocate_lw_diag(row_length, rows, model_levels,             &
        cloud_levels, 1 )

      IF (l_timestep) THEN
        ALLOCATE(lw_incs_local(row_length, rows, 0:model_levels,        &
                                                        n_lwcall))
        ALLOCATE(lwsea_local(row_length, rows, n_lwcall))
        ALLOCATE(olr_local(row_length, rows, n_lwcall))
        ALLOCATE(lw_down_local(row_length, rows, n_lwcall))
        ALLOCATE(top_abs_lw(row_length, rows, n_lwcall))
      END IF

      IF (l_radiance) THEN
        DO j_lw = 2, n_lwcall
          CALL allocate_lw_diag( row_length, rows,                      &
            model_levels, cloud_levels, j_lw )
        END DO
      END IF

    END IF

    IF (l_rad_step_diag) THEN
      IF (l_timestep.OR.l_forcing) THEN
        CALL allocate_lw_diag(row_length, rows, model_levels,           &
          cloud_levels, 2 )
      END IF
    END IF

! DEPENDS ON: timer
    IF (ltimer) CALL timer ('AP1R LW Rad  ',5)

    IF (l_rad_step_prog.OR.l_rad_step_diag) THEN

!     Calculate JULES mean land emissivity for use in t_rad_surf and
!     dOLR calculations and to pass down to R2_SET_SURFACE_FIELD    
!     Note that, unlike the albedo in the SW, the emissivity is not
!     recalculated from frac_control on a diagnsotic call. Given the
!     connection between surface temperatures and emissivity in the LW,
!     it is less clear how such surface forcing might be calculated. 
!     Since also the impact would be much smaller, no provision is 
!     made for calculation of forcing due to land-use changes in the LW.

!     For the following logic recall the definition of the mappings:
!     [tile,type]_index goes to the appropriate index over land points
!     and land_index then maps that to the grid-point.

      surf_emission(:,:) = 0.0
      t_rad_surf(:,:)    = 0.0
      t_rad_land(:,:)    = 0.0
      emis_land(:,:)     = 0.0

      IF (l_aggregate) THEN
!
!       First aggregate the emissivity under snow-free conditions.
        DO n=1, npft
          DO point = 1, type_pts(n)
            l = type_index(point, n)
            i = land_index_i(l)
            j = land_index_j(l)
            emis_land(i,j) = emis_land(i,j) +  &
              frac(l,n) * emis_pft(n)
          END DO
        END DO
        DO n=1, nnvg
          DO point = 1, type_pts(n+npft)
            l = type_index(point, n+npft)
            i = land_index_i(l)
            j = land_index_j(l)
            emis_land(i,j) = emis_land(i,j) +  &
              frac(l,n+npft) * emis_nvg(n)
          END DO
        END DO
!
!       Since JULES does not currently allow for the emissivity
!       of snow, the surface emission is calculated with the
!       snow-free emissivity.
!
        DO point = 1, tile_pts(1)
          l = tile_index(point, 1)
          i = land_index_i(l)
          j = land_index_j(l)
          IF (l_dolr_land_black) THEN
            surf_emission(i,j) = surf_emission(i,j) +  &    
              tstar_tile(l,1)**4    
          ELSE
            surf_emission(i,j) = surf_emission(i,j) +  &    
              emis_land(i,j) * tstar_tile(l,1)**4    
          END IF
          IF (l_t_land_nosnow) THEN
            emis_nosnow = emis_land(i,j)
          END IF
          IF (l_rad_snow_emis) THEN
            emis_land(i,j) = emis_land(i,j) +  &
              (emis_nvg(ice-npft) - emis_land(i,j)) *   &
              ( MAX(0.0, snow_tile(l,1)) /              &
              ( MAX(0.0, snow_tile(l,1)) +              &
              10.0 * z0_tile(l,1) * rho_snow_const ) )
          END IF
!         Now temporarily over-write tstar_tile to replicate
!         existing science.
          tstar_tile(l,1) = t_surf(i,j)
          IF (l_t_land_nosnow) THEN
            t_rad_land(i,j) = t_rad_land(i,j) +         &    
              emis_nosnow * tstar_tile(l,1)**4    
          ELSE
            t_rad_land(i,j) = t_rad_land(i,j) +         &    
              emis_land(i,j) * tstar_tile(l,1)**4    
          END IF
        END DO
!
      ELSE
!
!       Types match tiles here.
        emis_tiles(1:npft) = emis_pft
        emis_tiles(npft+1:ntype) = emis_nvg
        DO n = 1, ntype
          DO point = 1, tile_pts(n)
            l = type_index(point, n)
            i = land_index_i(l)
            j = land_index_j(l)
!           Calculate the surface emission without adjusting for
!           the emissivity of snow to match what JULES will add
!           back.
            emis_here = emis_tiles(n)
            IF (l_dolr_land_black) THEN
              surf_emission(i,j) = surf_emission(i,j) +   &
                tile_frac(l,n) * tstar_tile(l,n)**4
            ELSE
              surf_emission(i,j) = surf_emission(i,j) +   &
                tile_frac(l,n) * emis_here * tstar_tile(l,n)**4
            END IF
!           Adjust for snow and accumulate into calculation
!           of the radiative temperature.
            IF (l_t_land_nosnow) THEN
              emis_nosnow = emis_here
            END IF
            IF (l_rad_snow_emis) THEN
              emis_here = emis_here +                     &
                (emis_nvg(ice-npft) - emis_here) *        &
                ( MAX(0.0, snow_tile(l,n)) /              &
                ( MAX(0.0, snow_tile(l,n)) +              &
                10.0 * z0_tile(l,n) * rho_snow_const ) )
            END IF
            emis_land(i,j) = emis_land(i,j) +             &
              frac(l, n) * emis_here
            IF (l_t_land_nosnow) THEN
              t_rad_land(i,j) = t_rad_land(i,j) +         &
                tile_frac(l,n) * emis_nosnow * tstar_tile(l,n)**4    
            ELSE
              t_rad_land(i,j) = t_rad_land(i,j) +         &
                tile_frac(l,n) * emis_here * tstar_tile(l,n)**4    
            END IF
          END DO
        END DO
!
      END IF

      WHERE (flandg > 0.0)  
        surf_emission = flandg * surf_emission
        t_rad_surf =  t_rad_surf + flandg * t_rad_land
        t_rad_land = SQRT( SQRT( t_rad_land / emis_land ) )
      END WHERE


! Effective surface radiative temperature over sea ice
      t_rad_sice(:,:) = 0.0
      DO n = 1, nice_use
        DO point = 1, sice_pts_ncat(n)
          l = sice_index_ncat(point,n)
          i = ssi_index_i(l)
          j = ssi_index_j(l)
          t_rad_sice(i,j) = t_rad_sice(i,j) +                           & 
                 ice_fract_cat(i,j,n) * emis_sice *                     & 
                 tstar_sice_cat(i,j,n)**4
        END DO
      END DO
!     The test for the presence of sea ice in atmos_physics1 is
!     whether ice_fract exceeds 0.
      WHERE (ice_fract > 0.0)  
        t_rad_surf =  t_rad_surf + (1.0 - flandg) * t_rad_sice
        surf_emission =  surf_emission + (1.0 - flandg) * t_rad_sice
        t_rad_sice = SQRT( SQRT( t_rad_sice /                           &
          (emis_sice * ice_fract ) ) )
      END WHERE
!
!     Contributions from open sea
      WHERE (flandg < 1.0) 
        t_rad_surf =  t_rad_surf +                                      &
                      (1.0 - flandg) * (1.0 - ice_fract) *              &
                      emis_sea * tstar_sea**4
        surf_emission =  surf_emission +                                &
                      (1.0 - flandg) * (1.0 - ice_fract) *              &
                      emis_sea * tstar_sea**4
      END WHERE
!
      t_rad_surf = SQRT( SQRT( t_rad_surf / (flandg * emis_land +       &
                   (1.0 - flandg) * ( ice_fract * emis_sice +           &
                   (1.0 - ice_fract) * emis_sea) ) ) )
!
      IF ( l_ctile .AND. .NOT.l_quad_t_coast ) THEN
!       Use linear averaging at coastal points for consistency.
        WHERE ( (flandg > 0) .AND. (.NOT. land0p5) )
          t_rad_surf = (1.0 - flandg) * (1.0 - ice_fract) * tstar_sea
        END WHERE
        DO n = 1, nice_use
          DO point = 1, sice_pts_ncat(n)
            l = sice_index_ncat(point,n)
            i = ssi_index_i(l)
            j = ssi_index_j(l)
            IF ( (flandg(i,j) > 0) .AND. (.NOT. land0p5(i,j)) ) &
              t_rad_surf(i,j) = t_rad_surf(i,j) +  &
                (1.0 - flandg(i,j)) * &
                ice_fract_cat(i,j,n) * tstar_sice_cat(i,j,n)
          END DO
        END DO
        DO n = 1, ntype
          DO point = 1, tile_pts(n)
            l = type_index(point, n)
            i = land_index_i(l)
            j = land_index_j(l)
            IF ( (flandg(i,j) > 0) .AND. (.NOT. land0p5(i,j)) ) &
              t_rad_surf(i,j) = t_rad_surf(i,j) + &
                flandg(i,j) * tile_frac(l,n) * tstar_tile(l,n)
          END DO
        END DO
      END IF
!
      IF ( l_t_rad_solid ) THEN
!       Use deprecated common solid temperature for land and sea ice
!       Only coastal points require modification.
        t_rad_solid = 0.0
        DO n = 1, nice_use
          DO point = 1, sice_pts_ncat(n)
            l = sice_index_ncat(point,n)
            i = ssi_index_i(l)
            j = ssi_index_j(l)
            IF (flandg(i,j) > 0) &
              t_rad_solid(i,j) = t_rad_solid(i,j) +  &
                (1.0 - flandg(i,j)) * ice_fract_cat(i,j,n) * &
                emis_sice *tstar_sice_cat(i,j,n)**4
          END DO
        END DO
        DO j = 1, rows
          DO i = 1, row_length
            IF ( (flandg(i,j) > 0) .AND. (flandg(i,j) < 1) ) THEN
              t_rad_solid(i,j) = t_rad_solid(i,j) +  &
                flandg(i,j) * emis_land(i,j) * t_rad_land(i,j)**4
              t_rad_solid(i,j) = SQRT( SQRT( t_rad_solid(i,j) / &
                ( (1.0 - flandg(i,j)) * ice_fract(i,j) * emis_sice + &
                flandg(i,j) * emis_land(i,j) ) ) )
!             Over-write.
              t_rad_land(i,j) = t_rad_solid(i,j)
              t_rad_sice(i,j) = t_rad_solid(i,j)
            END IF
          END DO
        END DO
      END IF


      DO j_lw = n_lwcall, 1, -1
! Set COSP flag only on prognostic radiation steps
        IF (L_COSP_in.AND.(j_lw == 1).AND.L_rad_step_prog) THEN
          L_cosp=.TRUE.
        ELSE
          L_cosp=.FALSE.
        END IF

! For each LW call set the correct mass mixing ratios and switch the
! aerosols on/off, &c.  Note that as for the SW, in the case of the
! aerosols the reference  state is 'no aerosols' and hence they are
! switched off in the diagnostic call.

! As in the SW, since the "default" for each variable is to have it
! set in the diagnostic call to the same as in the prognostic, we set
! them all unconditionally to the prognostic & then re-set as needed.

! Set default values for greenhouse gases
        j_co2_mmr    = co2_mmr
        j_n2o_mmr    = n2o_mix_ratio
        j_ch4_mmr    = ch4_mix_ratio
        j_cfc11_mmr  = cfc11_mix_ratio
        j_cfc12_mmr  = cfc12_mix_ratio
        j_c113_mmr   = c113mmr
        j_c114_mmr   = c114mmr
        j_hcfc22_mmr = hcfc22mmr
        j_hfc125_mmr = hfc125mmr
        j_hfc134_mmr = hfc134ammr
        j_ozone      = ozone

! Set default flags for aerosols
        j_l_sulpc_so2          = l_sulpc_so2
        j_l_use_sulpc_direct   = l_use_sulpc_direct
        j_l_use_seasalt_direct = l_use_seasalt_direct
        j_l_soot               = l_soot
        j_l_use_soot_direct    = l_use_soot_direct
        j_l_biomass            = l_biomass
        j_l_use_bmass_direct   = l_use_bmass_direct
        j_l_ocff               = l_ocff
        j_l_use_ocff_direct    = l_use_ocff_direct
        j_l_nitrate            = l_nitrate
        j_l_use_nitrate_direct = l_use_nitrate_direct
        j_l_dust               = l_dust
        j_l_use_dust           = l_use_dust
        j_l_use_biogenic       = l_use_biogenic
        j_l_climat_aerosol     = l_climat_aerosol
        j_l_murk_rad           = l_murk_rad
        j_l_use_arcl           = l_use_arcl
        j_n_arcl_species       = n_arcl_species
        j_l_use_ukca_radaer    = l_ukca_radaer

! Turn aerosols off for the "cloud only" call in order to minimise
! computational cost. Clear-sky and other diagnostics have already
! been turned off in set_lwdiag_logic.
        IF (l_rad_perturb.AND.(j_lw==2)) THEN
          j_l_sulpc_so2          = .FALSE.
          j_l_use_sulpc_direct   = .FALSE.
          j_l_use_seasalt_direct = .FALSE.
          j_l_soot               = .FALSE.
          j_l_use_soot_direct    = .FALSE.
          j_l_biomass            = .FALSE.
          j_l_use_bmass_direct   = .FALSE.
          j_l_ocff               = .FALSE.
          j_l_use_ocff_direct    = .FALSE.
          j_l_nitrate            = .FALSE.
          j_l_use_nitrate_direct = .FALSE.
          j_l_dust               = .FALSE.
          j_l_use_dust           = .FALSE.
          j_l_use_biogenic       = .FALSE.
          j_l_climat_aerosol     = .FALSE.
          j_l_murk_rad           = .FALSE.
          j_l_use_arcl           = .FALSE.
          j_n_arcl_species       = 0
          j_l_use_ukca_radaer    = .FALSE.
        END IF

        l_call_lwrad=.FALSE.

! if timestepping required
        IF (l_timestep) THEN

          IF ((j_lw==1).AND.(l_rad_step_prog)) THEN
            a_lw_radstep=a_lw_radstep_prog
            l_call_lwrad=.TRUE.
          END IF
          IF ((j_lw==2).AND.(l_rad_step_diag)) THEN
            IF (l_rad_perturb.AND.l_rad_step_prog) THEN
              a_lw_radstep=a_lw_radstep_prog
            ELSE
              a_lw_radstep=a_lw_radstep_diag
            END IF
            l_call_lwrad=.TRUE.
          END IF

! If Forcing required
        ELSE IF (l_forcing) THEN

! Decide whether it is necessary to call radiation
          IF ((j_lw==1).AND.(l_rad_step_prog)) THEN
            a_lw_radstep=a_lw_radstep_prog
            l_call_lwrad=.TRUE.
          END IF
          IF ((j_lw==2).AND.(l_rad_step_diag)) THEN
            a_lw_radstep=a_lw_radstep_prog
            l_call_lwrad=.TRUE.
          END IF

          IF ( j_lw > 1 ) THEN ! Diagnostic Call

! Greenhouse Gases
            IF ( c2c_co2(j_lw-1) ) j_co2_mmr = co2_mmr_d
            IF ( c2c_n2o(j_lw-1) ) j_n2o_mmr = n2o_mix_ratio_d
            IF ( c2c_ch4(j_lw-1) ) j_ch4_mmr = ch4_mix_ratio_d
            IF ( c2c_cfc11(j_lw-1)) j_cfc11_mmr= cfc11_mix_ratio_d
            IF ( c2c_cfc12(j_lw-1)) j_cfc12_mmr= cfc12_mix_ratio_d
            IF ( c2c_c113(j_lw-1) ) j_c113_mmr = c113mmr_d
            IF ( c2c_hcfc22(j_lw-1) ) j_hcfc22_mmr = hcfc22mmr_d
            IF ( c2c_hfc125(j_lw-1) ) j_hfc125_mmr = hfc125mmr_d
            IF ( c2c_hfc134(j_lw-1) ) j_hfc134_mmr = hfc134ammr_d

! Ozone: Note that aerosol contains the MURK array.  It is in fact what
! is called d1(MURK) in Simon Tett's old vn4.5 mod.  Care needs to be
! taken so that MURK has the correct dimensions before it arrives in
! RAD_CTL2.
            IF ( c2c_o3(j_lw-1) ) j_ozone = aerosol

! Aerosols
            IF ( c2c_sulpc_d(j_lw-1) ) j_l_sulpc_so2          = .FALSE.
            IF ( c2c_sulpc_d(j_lw-1) ) j_l_use_sulpc_direct   = .FALSE.
            IF ( c2c_seas_d(j_lw-1) )  j_l_use_seasalt_direct = .FALSE.
            IF ( c2c_soot_d(j_lw-1) )  j_l_soot               = .FALSE.
            IF ( c2c_soot_d(j_lw-1) )  j_l_use_soot_direct    = .FALSE.
            IF ( c2c_bmb_d(j_lw-1) )   j_l_biomass            = .FALSE.
            IF ( c2c_bmb_d(j_lw-1) )   j_l_use_bmass_direct   = .FALSE.
            IF ( c2c_ocff_d(j_lw-1) )  j_l_ocff               = .FALSE.
            IF ( c2c_ocff_d(j_lw-1) )  j_l_use_ocff_direct    = .FALSE.
            IF ( c2c_nitr_d(j_lw-1) )  j_l_nitrate            = .FALSE.
            IF ( c2c_nitr_d(j_lw-1) )  j_l_use_nitrate_direct = .FALSE.
            IF ( c2c_dust_d(j_lw-1) )  j_l_dust               = .FALSE.
            IF ( c2c_dust_d(j_lw-1) )  j_l_use_dust           = .FALSE.
            IF ( c2c_biog_d(j_lw-1) )  j_l_use_biogenic       = .FALSE.
            IF ( c2c_ukca_d(j_lw-1) )  j_l_use_ukca_radaer    = .FALSE.

          END IF
        ELSE
          a_lw_radstep=a_lw_radstep_prog
          l_call_lwrad=.TRUE.
        END IF

        IF (l_call_lwrad) THEN

! DEPENDS ON: prelim_lwrad
          CALL prelim_lwrad(                                            &
! Parallel variables
            at_extremity, n_proc,                                       &
! Model Dimensions
            row_length,rows,model_levels,                               &
! Model Switches
            model_domain,mt_global,                                     &
            l_rad_deg, lw_control(1)%l_extra_top,                       &
            lw_control(j_lw)%l_subsample,                               &
            lw_control(j_lw)%l_geostationary,                           &
! Time stepping Information
            timestep_number,a_lw_radstep,                               &
! ancillary fields and fields needed to be kept from timestep to
! timestep
            lw_incs,                                                    &
! Satellite Geometry
            lw_control(j_lw)%min_view_lon,                              &
            lw_control(j_lw)%max_view_lon,                              &
            lw_control(j_lw)%min_view_lat,                              &
            lw_control(j_lw)%max_view_lat,                              &
! Number of Call
            j_lw,                                                       &
! Other variables
            true_latitude,true_longitude,                               &
            seconds_since_midnight,                                     &
            rad_mask, list_lw_points, first_data_interp,                &
            first_row,last_row,diag_row_list(1,j_lw),                   &
            diag_col_list(1,j_lw),                                      &
            lw_points,                                                  &
            olr, lw_down, lwsea, top_absorption )

!Parallelise over segments.
!$OMP  PARALLEL DEFAULT(SHARED)                                              &
!$OMP& PRIVATE(i,j, lit_points,start_point,first_point,first_point_dust_a,   &
!$OMP& first_point_dust_b,                                                   &
!$OMP& first_point_sulpc,first_point_soot,first_point_biomass,               &
!$OMP& first_point_seasalt,nd_rad_pts,n_rad_layers,error_code,               &
!$OMP& first_point_biogenic, first_point_ocff, first_point_arcl,             &
!$OMP& first_point_nitrate,l_scale_inc, n_channel, map_channel,nd_cloud_type,&
!$OMP& nd_cloud_component, nd_column, id_ct, nd_channel, nd_viewing_level,   &
!$OMP& nd_direction, nd_profile, nd_radiance_profile, nd_flux_profile,       &
!$OMP& nd_field_flux_diag, nd_field_rad_diag, nd_brdf_basis_fnc,             &
!$OMP& nd_brdf_trunc, nd_tile, nd_point_tile, n_viewing_level, viewing_level,&
!$OMP& n_viewing_direction, viewing_direction, ptr_local, first_point_ukca,  &
!$OMP& meta_segments, segments, first_point_local, seg_points_local, ipar)

!$OMP SINGLE
!Determine the number of threads in the parallel region
num_parallel_lw=1
!$ num_parallel_lw = omp_get_num_threads()
!$OMP END SINGLE
!Implicit barrier

!Determine the team member number.
!NB: this is numbered from 1 to follow the Fortran convention.
ipar=1
!$ ipar = omp_get_thread_num()+1

          !Set up segmentation information
          CALL segments_mod_seg_meta(meta_segments, ipar, num_parallel_lw, &
            lw_points, a_lw_seg_size, a_lw_segments)

! Allocate space for segmentation arrays
          ALLOCATE( first_point_local( meta_segments%num_segments ) )
          ALLOCATE(  seg_points_local( meta_segments%num_segments ) )
          ALLOCATE(          segments( meta_segments%num_segments ) )

          !Fill the segments array
          CALL segments_mod_segments( segments, meta_segments, &
                 lw_points, row_length, rows,                  &
                 list_points=list_lw_points)

          first_point_local(:) = segments(:)%fp
          seg_points_local(:)  = segments(:)%use_points

          DO i = 1, meta_segments%num_segments

            first_point = segments(i)%start_index 

            !Additional work on list_points
            DO j = first_point, segments(i)%end_index
              list_lw_points(j) = list_lw_points(j)-segments(i)%start+1
            END DO

! Set the pointer to the beginning of the current segment.
! This is done solely to reduce the number of continuation
! lines required by the call.
            ptr_local=first_point_local(i)

! Set the actual size of arrays in the radiation code:
! for some architectures (such as that of Cray vector
! machines) on odd size is preferred to avoid memory
! bank conflicts.
            nd_rad_pts=2*(seg_points_local(i)/2)+1

! Set the number of layers seen in the radiation code.
! This may optionally be 1 greater than the number used
! in the rest of the model to avoid spurious effects
! resulting from the upper boundary (principally in
! stratospheric configurations).
            IF (lw_control(1)%l_extra_top) THEN
              n_rad_layers=model_levels+1
            ELSE
              n_rad_layers=model_levels
            END IF

! Set the first point of the dust arrays to be used.
            IF (l_dust) THEN
              IF(l_twobin_dust) THEN
                first_point_dust_a = first_point_local(i)
                first_point_dust_b = 1
              ELSE
                first_point_dust_a = first_point_local(i)
                first_point_dust_b = first_point_local(i)
              END IF
            ELSE
              first_point_dust_a = 1
              first_point_dust_b = 1
            END IF

! Set the first point of the biogenic array
            IF (l_use_biogenic) THEN
              first_point_biogenic=first_point_local(i)
            ELSE
              first_point_biogenic=1
            END IF

! Set the first points of the arrays of sulphates to be used.

! A separate assignment is necessary since
! not be of the full size unless the sulphur cycle is on.
            IF (l_sulpc_so2 .OR. l_use_sulpc_indirect_lw) THEN
              first_point_sulpc=first_point_local(i)
            ELSE
              first_point_sulpc=1
            END IF

            IF (l_soot) THEN
              first_point_soot=first_point_local(i)
            ELSE
              first_point_soot=1
            END IF

            IF (l_biomass .OR. l_use_bmass_indirect) THEN
              first_point_biomass=first_point_local(i)
            ELSE
              first_point_biomass=1
            END IF

            IF (l_ocff .OR. l_use_ocff_indirect) THEN
              first_point_ocff=first_point_local(i)
            ELSE
              first_point_ocff=1
            END IF

            IF (l_use_seasalt_indirect .OR. l_use_seasalt_direct) THEN
              first_point_seasalt=first_point_local(i)
            ELSE
              first_point_seasalt=1
            END IF

            IF (n_arcl_species > 0) THEN
              first_point_arcl = first_point_local(i)
            ELSE
              first_point_arcl = 1
            END IF

            IF (l_nitrate .OR. l_use_nitrate_indirect) THEN
              first_point_nitrate=first_point_local(i)
            ELSE
              first_point_nitrate=1
            END IF

            IF (l_ukca_radaer) THEN
              first_point_ukca=first_point_local(i)
            ELSE
              first_point_ukca=1
            END IF
            
! Set dynamically determined array sizes required for
! radiation calculations. There is redundancy here in
! that the viewing directions are calculated for the
! whole PE in each segemnt: ideally, this will be
! tidied.
            ALLOCATE(map_channel(lw_spectrum(j_lw)%n_band))

! DEPENDS ON: r2_set_rad_dim
            CALL r2_set_rad_dim(row_length, rows,                       &
              true_latitude, true_longitude,                            &
              .FALSE., sindec, seconds_since_midnight,                  &
              lw_control(j_lw)%l_geostationary,                         &
              lw_control(j_lw)%sat_hgt,                                 &
              lw_control(j_lw)%sat_lat,                                 &
              lw_control(j_lw)%sat_lon,                                 &
              lw_control(j_lw)%i_cloud,                                 &
              lw_control(j_lw)%i_angular_integration,                   &
              lw_control(j_lw)%i_sph_mode,                              &
              row_length*rows, seg_points_local(i),                     &
              model_levels, global_cloud_top,                           &
              nd_cloud_type, nd_cloud_component,                        &
              lw_control(j_lw)%l_extra_top,                             &
              id_ct, n_rad_layers, nd_column,                           &
              lw_spectrum(j_lw)%npd_band,                               &
              lw_spectrum(j_lw)%n_band,                                 &
              map_channel, n_channel, nd_channel,                       &
              nd_viewing_level, n_viewing_level,                        &
              viewing_level,                                            &
              nd_direction, n_viewing_direction,                        &
              viewing_direction,                                        &
              nd_brdf_basis_fnc, nd_brdf_trunc,                         &
              nd_profile, nd_flux_profile,                              &
              nd_radiance_profile,                                      &
              nd_field_flux_diag, nd_field_rad_diag,                    &
              l_ctile, nd_tile,                                         &
              nd_point_tile)

! To remain general heating rates are always
! converted to increments across a timestep:
            l_scale_inc = .TRUE.


! DEPENDS ON: r2_lwrad3z
            CALL r2_lwrad3z(error_code,                                 &
! Input data
  q_n(ptr_local,1,1), j_co2_mmr, j_ozone(ptr_local,1,1),                &
  co2_dim1, co2_dim2, co2_3d(ptr_local,1,1), l_co2_3d,                  &

! chemical greenhouse gas fields
  ngrgas, grgas_field(ptr_local,1,1,1),                                 &
  j_n2o_mmr, j_ch4_mmr, j_cfc11_mmr, j_cfc12_mmr,                       &
  j_c113_mmr, j_c114_mmr, j_hcfc22_mmr, j_hfc125_mmr,                   &
  j_hfc134_mmr, t_n(ptr_local,1,1),t_rad_surf(ptr_local,1),             &
  t_rad_land(ptr_local,1), t_rad_sice(ptr_local,1),                     &
  tstar_sea(ptr_local,1),l_ctile,                                       &
  p_star(ptr_local,1),p_layer_boundaries(ptr_local,1,0),                &
  p_layer_centres(ptr_local,1,0),                                       &
  height_theta(ptr_local,1,0), height_rho(ptr_local,1,1),               &

! Options for COSP 
  L_cosp,                                                               & 

! Options for treating clouds
  global_cloud_top, l_inhom_cloud, inhom_cloud_lw,                      &
  dp_corr_strat, dp_corr_conv,                                          &

! Stratiform Cloud Fields
  l_pc2, area_cloud_fraction(ptr_local,1,1), cf_n(ptr_local,1,1),       &
  qcl_n(ptr_local,1,1),qcf_n(ptr_local,1,1),                            &
  n_drop_pot(ptr_local,1,1),                                            &

! Convective Cloud Fields
  cca(ptr_local,1,1),cclwp(ptr_local,1),                                &
  ccw(ptr_local,1,1), lcbase(ptr_local,1),                              &
  ccb(ptr_local,1), cct(ptr_local,1),                                   &

! Surface Fields. Only want the 0.5 threshold LAND mask
! and fractional land:
  land0p5(ptr_local,1),flandg(ptr_local,1),                             &
  ice_fract(ptr_local,1),snow_depth(ptr_local,1),                       &
  emis_land(ptr_local,1),                                               &

! Solar Fields
  cos_zenith_angle_2(ptr_local,1,j_lw),                                 &
  day_fraction_2(ptr_local,1,j_lw), scs,                                &

! Aerosol Fields
  j_l_climat_aerosol, l_clim_aero_hgt, l_hadgem1_clim_aero,             &
  zh(ptr_local,1), aero_bl_levels,                                      &
  j_l_dust, j_l_use_dust, dust_dim1, dust_dim2,                         &
  dust_1(first_point_dust_a,1), dust_2(first_point_dust_a,1),           &
  dust_3(first_point_dust_b,1), dust_4(first_point_dust_b,1),           &
  dust_5(first_point_dust_b,1), dust_6(first_point_dust_b,1),           &
  j_l_use_biogenic, biogenic_dim1, biogenic_dim2,                       &
  local_biogenic(first_point_biogenic, 1),                              &
  j_l_sulpc_so2, j_l_use_sulpc_direct, l_use_sulpc_indirect_lw,         &
  sulp_dim1, sulp_dim2,                                                 &
  accum_sulphate(first_point_sulpc, 1),                                 &
  aitken_sulphate(first_point_sulpc, 1),                                &
  diss_sulphate(first_point_sulpc, 1),                                  &
  sea_salt_film(first_point_seasalt,1,1),                               &
  sea_salt_jet(first_point_seasalt,1,1),                                &
  l_use_seasalt_indirect, j_l_use_seasalt_direct,                       &
  salt_dim1*salt_dim2, salt_dim3, j_l_soot, j_l_use_soot_direct,        &
  soot_dim1, soot_dim2, fresh_soot(first_point_soot, 1),                &
  aged_soot(first_point_soot, 1), j_l_biomass, j_l_use_bmass_direct,    &
  bmass_dim1, bmass_dim2, fresh_bmass(first_point_biomass, 1),          &
  aged_bmass(first_point_biomass, 1),                                   &
  cloud_bmass(first_point_biomass, 1), l_use_bmass_indirect,            &
  j_l_ocff, j_l_use_ocff_direct, ocff_dim1, ocff_dim2,                  &
  fresh_ocff(first_point_ocff, 1), aged_ocff(first_point_ocff, 1),      &
  cloud_ocff(first_point_ocff, 1), l_use_ocff_indirect,                 &
  j_l_nitrate, j_l_use_nitrate_direct, nitrate_dim1, nitrate_dim2,      &
  accum_nitrate(first_point_nitrate, 1),                                &
  diss_nitrate(first_point_nitrate, 1), l_use_nitrate_indirect,         &
  j_l_use_arcl, arcl_dim1, arcl_dim2, j_n_arcl_species,                 &
  n_arcl_compnts, i_arcl_compnts,local_arcl(first_point_arcl,1,1),      &
  aerosol(first_point_local(i),1,1), j_l_murk_rad,                      &
  ntot_land, ntot_sea,                                                  &
  j_l_use_ukca_radaer, ukca_radaer, ukca_dim1, ukca_dim2,               &
  local_ukca_mmr(first_point_ukca, 1, 1),                               &
  local_ukca_cvl(first_point_ukca, 1, 1),                               &
  local_ukca_dry(first_point_ukca, 1, 1),                               &
  local_ukca_wet(first_point_ukca, 1, 1),                               &
  local_ukca_rho(first_point_ukca, 1, 1),                               &
  local_ukca_vol(first_point_ukca, 1, 1),                               &
  local_ukca_wtv(first_point_ukca, 1, 1),                               &
  local_ukca_nbr(first_point_ukca, 1, 1),                               &

! time
  previous_time,                                                        &

! grid-dependent arrays
  true_latitude(ptr_local,1),                                           &

! Level of tropopause
  trindx(ptr_local, 1),                                                 &

! Spectral data
  lw_spectrum(j_lw),                                                    &

! Algorithmic options
  lw_control(j_lw), timestep, l_mod_k_flux, l_scale_inc,                &
  list_lw_points(first_point),                                          &

! Viewing geometry for satellite simulations
  n_viewing_direction,                                                  &
  viewing_direction(ptr_local, 1, 1),                                   &
  viewing_direction(ptr_local, 1, 2),                                   &
  n_viewing_level, viewing_level,                                       &

! All diagnostics
  n_channel, map_channel,                                               &
  lw_diag(j_lw), diag_row_list(first_point,j_lw),                       &
  diag_col_list(first_point,j_lw),                                      &

! Physical Dimensions
  seg_points_local(i),model_levels,n_rad_layers,                        &
  global_cloud_top,wet_model_levels,ozone_levels,row_length,            &
  rows,row_length*rows, nd_field_flux_diag,                             &
  nd_field_rad_diag, nd_profile, n_rad_layers, nd_column,               &
  n_cca_levels,                                                         &
  nd_channel,                                                           &
  nd_flux_profile, nd_radiance_profile,                                 &
  nd_viewing_level, nd_direction,                                       &
  nd_cloud_component, nd_cloud_type,                                    &
  nd_brdf_basis_fnc, nd_brdf_trunc,                                     &
  nd_point_tile, nd_tile, id_ct, n_ukca_mode, n_ukca_cpnt,              &

! Output data
  olr(ptr_local, 1), lw_down(ptr_local, 1),                             &
  top_absorption(ptr_local, 1), lwsea(ptr_local, 1),                    &
  lw_incs(ptr_local,1,0),                                               &

! Data for calculation of layer masses
  rho_r2_nh(ptr_local,1),                                               &
  r_rho_levels_nh(ptr_local,1), r_theta_levels_nh(ptr_local,0),         &
  q_n(ptr_local,1,1), qcl_n(ptr_local,1,1),qcf_n(ptr_local,1,1),        &
  qcf2_n(ptr_local,1,1),qrain_n(ptr_local,1,1),                         &
  qgraup_n(ptr_local,1,1),                                              &
  l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio,                &

! COSP input arguments
  cosp_gbx)

            DEALLOCATE(map_channel)

          END DO ! end loop over long-wave segments

! Deallocate the segmentation arrays
          DEALLOCATE( first_point_local )
          DEALLOCATE(  seg_points_local )
          DEALLOCATE(  segments         )

!$OMP END PARALLEL


! Radiative fluxes may not have been calculated at all
! points. We now fill in as required.

       IF (.NOT. l_vatpoles) THEN
! At the North Pole in a global domain calculations
! are performed only at the first point, so the rest of
! the row must be filled in.
          l_complete_north=(model_domain==mt_global) .AND.              &
                      ( at_extremity(pnorth) )

! At the South Pole in a global domain calculations
! are performed only at the first point, so the rest of
! the row must be filled in.
          l_complete_south=(model_domain==mt_global) .AND.              &
                      ( at_extremity(psouth) )
       ELSE
          l_complete_north=.FALSE.
          l_complete_south=.FALSE.
       END IF ! vatpoles

! When spatial degradation is performed fields must be
! filled in at alternate points.
          l_complete_deg = ( l_rad_deg )

!         Set addressing limits for spatial degradation.
          IF ( l_complete_deg ) THEN
            first_row=1
            last_row=rows
            IF (.NOT. l_vatpoles) THEN
            IF (model_domain  ==  mt_global) THEN
              IF (at_extremity(pnorth)) THEN
                last_row=rows-1
              END IF
              IF (at_extremity(psouth)) THEN
                first_row=2
              END IF
            END IF
            END IF ! vatpoles
          END IF

! DEPENDS ON: fill_missing_data_lw
          CALL fill_missing_data_lw(                                    &
            off_x, off_y, row_length, rows, model_levels,               &
            cloud_levels, lw_spectrum(j_lw)%n_aod_wavel,                &
            first_row,last_row,                                         &
            first_data_interp, es_space_interp,                         &
            l_complete_north, l_complete_south,                         &
            l_complete_deg, n_channel, j_lw,                            &
            lw_control(1)%l_extra_top,                                  &
            lw_incs, olr, lw_down, lwsea, top_absorption )

          IF (l_rad_perturb.AND.(j_lw==1)) THEN
            lw_incs_local(:,:,:,1) = lw_incs                            &
                                   - lw_incs_local(:,:,:,2)
            olr_local(:,:,1) = olr                                      &
                             - olr_local(:,:,2)
            lw_down_local(:,:,1) = lw_down                              &
                                 - lw_down_local(:,:,2)
            lwsea_local(:,:,1) = lwsea                                  &
                               - lwsea_local(:,:,2)
            IF (lw_control(1)%l_extra_top) THEN
              top_abs_lw(:,:,1) = top_absorption                        &
                                - top_abs_lw(:,:,2)
            END IF
          ELSE IF (l_timestep) THEN
            lw_incs_local(:,:,:,j_lw)=lw_incs
            olr_local(:,:,j_lw)=olr
            lw_down_local(:,:,j_lw)=lw_down
            lwsea_local(:,:,j_lw)=lwsea
            IF (lw_control(1)%l_extra_top) THEN
              top_abs_lw(:,:,j_lw)=top_absorption
            END IF
          END IF

        END IF ! if radiation call required
      END DO ! loop over radiation calls


      IF (l_rad_perturb.AND.l_rad_step_prog) THEN

!       For this case the full fluxes are already set.

      ELSE IF (l_timestep.AND.l_rad_step_diag) THEN
        lw_incs(:,:,:) = lw_incs_local(:,:,:,1)                         &
                       + lw_incs_local(:,:,:,2)
        olr(:,:) = olr_local(:,:,1)                                     &
                 + olr_local(:,:,2)
        lw_down(:,:) = lw_down_local(:,:,1)                             &
                     + lw_down_local(:,:,2)
        lwsea(:,:) = lwsea_local(:,:,1)                                 &
                   + lwsea_local(:,:,2)

        IF (lw_control(1)%l_extra_top) THEN
          top_absorption(:,:) = top_abs_lw(:,:,1)                       &
                              + top_abs_lw(:,:,2)
        END IF
      END IF


      IF ( l_ctile ) THEN
        DO j = 1, rows
          DO i = 1, row_length
            surflw(i,j) = lw_incs(i,j,0)                                &
                    + (1.-flandg(i,j))*lwsea(i,j)
            IF (flandg(i,j) == 1.0) lwsea(i,j)=rmdi
          END DO
        END DO
      ELSE
        DO j = 1, rows
          DO i = 1, row_length
            surflw(i,j) = lw_incs(i,j,0) + lwsea(i,j)
            IF (land_sea_mask(i,j)) lwsea(i,j)=rmdi
          END DO
        END DO
      END IF
!  
!     To allow for changes in surface temperature on non-radiation  
!     timesteps, dOLR is now defined as the grid-box mean OLR  
!     less the contributions to the upward LW flux at the surface  
!     from each portions of the surface, so conceptually  
!     emissivity * areal fraction * sigma T^4 for each portion.   
!     This is based on the effective temperatures seen by radiation,
!     contributions haveing been accumulated in surf_emission.
!  
      dOLR_rts = olr - sbcon * surf_emission
!  
    END IF ! end conditional on being a radiation timestep


! Is the PC2 cloud scheme being used?
    IF (l_pc2) THEN

! Reset _latest values to _n values
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            t_latest(i,j,k)   = t_n(i,j,k)
          END DO
        END DO
      END DO
      DO k = 1, wet_model_levels
        DO j = 1, rows
          DO i = 1, row_length
            q_latest(i,j,k)   = q_n(i,j,k)
            qcl_latest(i,j,k) = qcl_n(i,j,k)
            cf_latest(i,j,k)  = cf_n(i,j,k)
            cfl_latest(i,j,k) = cfl_n(i,j,k)
          END DO
        END DO
      END DO

! ----------------------------------------------------------------------
! Homogeneous forcing. Note the temperature increment from longwave
! is added in this routine
! ----------------------------------------------------------------------

      CALL pc2_homog_plus_turb(p_layer_centres(1,1,1),                  &
        wet_model_levels,                                               &
        timestep, t_latest, cf_latest, cfl_latest,                      &
        cff_latest, q_latest, qcl_latest, lw_incs(1,1,1),               &
        zeros, zeros, zeros, 0.0, 0.0,                                  &
        l_mixing_ratio)

! Add increments from the homogeneous forcing to the increment variables
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            t_inc(i,j,k) = t_inc(i,j,k) + t_latest(i,j,k)-t_n(i,j,k)
          END DO
        END DO
      END DO
      DO k = 1, wet_model_levels
        DO j = 1, rows
          DO i = 1, row_length
            q_inc(i,j,k) = q_inc(i,j,k) + q_latest(i,j,k)-q_n(i,j,k)
            qcl_inc(i,j,k) = qcl_inc(i,j,k)                             &
                             + qcl_latest(i,j,k)-qcl_n(i,j,k)
            cf_inc(i,j,k) = cf_inc(i,j,k)                               &
                             + cf_latest(i,j,k)-cf_n(i,j,k)
            cfl_inc(i,j,k) = cfl_inc(i,j,k)                             &
                             + cfl_latest(i,j,k)-cfl_n(i,j,k)
          END DO
        END DO
      END DO

    ELSE  ! L_pc2

      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            t_inc(i,j,k) = t_inc(i,j,k) + lw_incs(i,j,k)
            t_latest(i,j,k) = t_n(i,j,k) + lw_incs(i,j,k)

          END DO
        END DO
      END DO

    END IF  ! L_pc2

! Get T_incr for output as STASH diagnostic
    IF ( l_t_incr_lw ) THEN  ! STASHflag set
      ALLOCATE ( t_incr_diagnostic(row_length,rows,model_levels) )
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            t_incr_diagnostic(i,j,k) =  lw_incs(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
    ELSE
      ALLOCATE ( t_incr_diagnostic(1,1,1) )
    END IF ! on STASHflag

! Set up radiative heating rates for 6A boundary layer code
    recip_timestep = 1./ timestep
    DO k = 1, bl_levels
      DO j = 1, rows
        DO i = 1, row_length
          rad_hr(i,j,1,k) = lw_incs(i,j,k) * recip_timestep
          rad_hr(i,j,2,k) = sw_incs(i,j,k)                              &
                   * cos_zenith_angle(i,j) * recip_timestep
        END DO
      END DO
    END DO

! Copy dOLR from last LW timestep
    DO j = 1, rows
      DO i = 1, row_length
        dOLR(i,j) = dOLR_rts(i,j)
      END DO
    END DO

! DEPENDS ON: timer
    IF (ltimer) CALL timer ('AP1R LW Rad  ',6)

! ----------------------------------------------------------------------
! Section RAD.2.1 Isccp diagnostics
!-----------------------------------------------------------------------

    IF ( l_rad_step_prog .AND. (lw_diag(1)%l_isccp_cf .OR.              &
      lw_diag(1)%l_isccp_cf_tau_0_to_p3 .OR.                            &
      lw_diag(1)%l_isccp_cf_tau_p3_to_1p3 .OR.                          &
      lw_diag(1)%l_isccp_cf_tau_1p3_to_3p6 .OR.                         &
      lw_diag(1)%l_isccp_cf_tau_3p6_to_9p4 .OR.                         &
      lw_diag(1)%l_isccp_cf_tau_9p4_to_23 .OR.                          &
      lw_diag(1)%l_isccp_cf_tau_23_to_60 .OR.                           &
      lw_diag(1)%l_isccp_cf_tau_ge_60 .OR.                              &
      lw_diag(1)%l_meanalbedocld .OR.                                   &
      lw_diag(1)%l_meantaucld .OR.                                      &
      lw_diag(1)%l_meanptop .OR.                                        &
      lw_diag(1)%l_totalcldarea) .AND.                                  &
      error_code  ==  0) THEN


! This code is repeated from the earlier call to sw_rad (the
! diagnostics are calculated at lit points only)

      IF ( daylight_points  >   0 ) THEN

!$OMP  PARALLEL DEFAULT(SHARED)                             &
!$OMP& PRIVATE   (meta_segments, segments, ipar)            &
!$OMP& PRIVATE   (nd_rad_pts,n_rad_layers)                  &
!$OMP& PRIVATE   (i,j, lit_points,start_point,first_point)  &
!$OMP& PRIVATE   (error_code)

!$OMP SINGLE
!Determine the number of threads in this parallel region
num_parallel_isccp=1
!$ num_parallel_isccp = omp_get_num_threads()
!$OMP END SINGLE
!Implicit barrier

!Determine the team member number
!NB: this is numbered from 1 to follow the Fortran convention.
ipar=1
!$ ipar = omp_get_thread_num()+1

        !Set up meta-information for segments
        CALL segments_mod_seg_meta(meta_segments, ipar,   &
                num_parallel_isccp,                       &
                daylight_points, a_sw_seg_size, a_sw_segments)

        !Allocate storage
        ALLOCATE( segments( meta_segments%num_segments) )

        !Find specific points in each segment
        CALL segments_mod_segments(segments,meta_segments,  &
          daylight_points, row_length, rows,                &
          list_points=list_daylight_points(:,1))
 
        DO i = 1, meta_segments%num_segments

          lit_points  = segments(i)%use_points
          start_point = segments(i)%start_index
          first_point = segments(i)%fp

          !Additional work on list_points. Note that although work is done on
          !this array, each thread and each segment should only affect one
          !element.
          DO j = start_point, segments(i)%end_index
            list_daylight_points_start(j) = list_daylight_points(j,1) &
                                           -segments(i)%start+1
          END DO

!         Set the actual size of arrays in the radiation code:
!         for some architectures (such as that of Cray vector
!         machines) on odd size is preferred to avoid memory
!         bank conflicts.
          nd_rad_pts=2*(lit_points/2)+1

!         Set the number of layers seen in the radiation code.
!         This may optionally be 1 greater than the number used
!         in the rest of the model to avoid spurious effects
!         resulting from the upper boundary (principally in
!         stratospheric configurations).
          IF (sw_control(1)%l_extra_top) THEN
            n_rad_layers=model_levels+1
          ELSE
            n_rad_layers=model_levels
          END IF

! DEPENDS ON: isccp
          CALL isccp(error_code, i,                                     &
! Mixing Ratios
            q_n(first_point,1,1),                                       &
! Pressures and Temperatures
            t_rad_surf(first_point,1),p_star(first_point,1),            &
            p_layer_boundaries(first_point,1,0),                        &
            p_layer_centres(first_point,1,0),t_n(first_point,1,1),      &
! Stratiform Cloud Fields
            area_cloud_fraction(first_point,1,1),                       &
            cf_n(first_point,1,1),                                      &
! Convective Cloud Fields
            cca(first_point,1,1), cclwp(first_point,1),                 &
            ccb(first_point,1), cct(first_point,1),                     &
! Solar Fields 
            day_fraction_2(first_point,1,1),                            &
            list_daylight_points_start(start_point), scs,               &
! Level of tropopause
            trindx(first_point,1),                                      &
! All diagnostics and associated arrays
            lw_diag(1), sw_diag(1),                                     &
            diag_row_list_sw(start_point),                              &
            diag_col_list_sw(start_point),                              &
! Dimensions
            lit_points,segments(i)%seg_points,model_levels,             &
            n_rad_layers,cloud_levels,wet_model_levels,ozone_levels,    &
            rows*row_length,nd_rad_pts,n_rad_layers,1,                  &
            n_cca_levels,                                               &
! Variables needed to calculate layer masses
            rho_r2_nh(first_point,1),                                   &
            r_rho_levels_nh(first_point,1),                             &
            r_theta_levels_nh(first_point,0),                           &
            q_n(first_point,1,1), qcl_n(first_point,1,1),               &
            qcf_n(first_point,1,1), qcf2_n(first_point,1,1),            &
            qrain_n(first_point,1,1), qgraup_n(first_point,1,1),        &
            l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio )

       END DO !Segments

       !Deallocate storage
       DEALLOCATE( segments )

!$OMP END PARALLEL


! Radiative fluxes may not have been calculated at all
! points. We now fill in as required.

      IF (.NOT. l_vatpoles) THEN
! At the North Pole in a global domain calculations
! are performed only at the first point if lit, so
! the rest of the row must be filled in if this point
! is lit.
        l_complete_north = ( model_domain == mt_global ) .AND.          &
                           ( at_extremity(pnorth) ) .AND.               &
                           ( switch(1,rows,1) )

! At the South Pole in a global domain calculations
! are performed only at the first point if lit, so
! the rest of the row must be filled in if this point
! is lit.
        l_complete_south = ( model_domain == mt_global ) .AND.          &
                           ( at_extremity(psouth) ) .AND.               &
                           ( switch(1,1,1) )
     ELSE
        l_complete_north=.FALSE.
        l_complete_south=.FALSE.
     END IF ! vatpoles

! When spatial degradation is performed fields must be
! filled in at alternate points.
        l_complete_deg = ( l_rad_deg ) .AND.                            &
                         ( tot_daylight_points > 0 )

! Set addressing limits for spatial degradation.
        IF ( l_complete_deg ) THEN
          first_row=1
          last_row=rows
          IF (.NOT. l_vatpoles) THEN
          IF (model_domain  ==  mt_global) THEN
            IF (at_extremity(pnorth)) THEN
              last_row=rows-1
            ELSE IF (at_extremity(psouth)) THEN
              first_row=2
            END IF
          END IF
          END IF ! vatpoles
        END IF

        first_data_interp = first_data_interp_sw

! Call appropriate subroutines to fill in missing data
! as required.
        IF ( l_complete_north .OR. l_complete_south .OR.                &
             l_complete_deg ) THEN

! Isccp diagnostics
! DEPENDS ON: fill_missing_data_isccp
          CALL fill_missing_data_isccp(                                 &
            off_x, off_y, row_length, rows,                             &
            first_row,last_row,                                         &
            first_data_interp, 1, es_space_interp,                      &
            l_complete_north, l_complete_south, l_complete_deg)
        END IF

      END IF  ! On daylit points

    END IF  ! If a radiation timestep and want ISCCP diagnostics

! ----------------------------------------------------------------------
! Section RAD.2.2 Long Wave Radiation Energy correction code
!-----------------------------------------------------------------------

    IF ((l_rad_step_prog.OR.l_rad_step_diag) .AND. l_emcorr) THEN

! Sum long wave fluxes into the atmosphere and
! add into the net diabatic fluxes into the
! atmosphere for use in the energy correction
! procedure

      IF (lw_control(1)%l_extra_top) THEN
!       The energy absorbed above the top of the model in
!       the radiation scheme does not contribute to the
!       energy absorbed, but the diagnostics are calculated
!       at the top of the atmosphere, so the net atmospheric
!       flux must be adjusted.
        DO j = 1, rows
          DO i = 1, row_length
            net_atm_flux(i,j) = -olr(i,j) - surflw(i,j)                 &
                                -top_absorption(i,j)
          END DO
        END DO
      ELSE
        DO j = 1, rows
          DO i = 1, row_length
            net_atm_flux(i,j) = -olr(i,j) - surflw(i,j)
          END DO
        END DO
      END IF

! DEPENDS ON: flux_diag
      CALL flux_diag(net_atm_flux, xx_cos_theta_latitude,               &
        row_length, rows ,off_x,off_y,1.0,                              &
        sum_eng_fluxes,radiation_tstep)

    END IF

! ----------------------------------------------------------------------
! Section RAD.2.3 Long Wave Radiation diagnostics
!-----------------------------------------------------------------------

    DO j_lw = 1, n_lwdiag

! Check that lw diagnostics requested this timestep

      IF (error_code  ==  0                                             &
        .AND. sf_calc(0,2)                                                   &
        ) THEN

        IF (l_timestep) THEN

          i_off=0

          IF (j_lw == n_lwdiag) THEN

            IF (l_rad_perturb.AND.l_rad_step_prog) THEN

              IF (lw_diag(2)%l_flux_up)                                 &
                lw_diag(1)%flux_up =                                    &
                lw_diag(1)%flux_up -                                    &
                lw_diag(2)%flux_up

              IF (lw_diag(2)%l_flux_down)                               &
                lw_diag(1)%flux_down =                                  &
                lw_diag(1)%flux_down -                                  &
                lw_diag(2)%flux_down

              IF (lw_diag(2)%l_net_flux_trop)                           &
                lw_diag(1)%net_flux_trop =                              &
                lw_diag(1)%net_flux_trop -                              &
                lw_diag(2)%net_flux_trop

              IF (lw_diag(2)%l_down_flux_trop)                          &
                lw_diag(1)%down_flux_trop =                             &
                lw_diag(1)%down_flux_trop -                             &
                lw_diag(2)%down_flux_trop

            END IF

            IF (l_rad_step_diag) THEN

              IF (lw_diag(2)%l_flux_up)                                 &
                lw_diag(2)%flux_up =                                    &
                lw_diag(2)%flux_up +                                    &
                lw_diag(1)%flux_up

              IF (lw_diag(2)%l_flux_down)                               &
                lw_diag(2)%flux_down =                                  &
                lw_diag(2)%flux_down +                                  &
                lw_diag(1)%flux_down

              IF (lw_diag(2)%l_net_flux_trop)                           &
                lw_diag(2)%net_flux_trop =                              &
                lw_diag(2)%net_flux_trop +                              &
                lw_diag(1)%net_flux_trop

              IF (lw_diag(2)%l_down_flux_trop)                          &
                lw_diag(2)%down_flux_trop =                             &
                lw_diag(2)%down_flux_trop +                             &
                lw_diag(1)%down_flux_trop

            END IF ! l_rad_step_diag
          END IF ! j_lw == n_lwdiag

        ELSE IF (l_forcing.AND.l_rad_step_diag) THEN

! Calculate LW Forcings. Note that forcings are only calculated
! for fluxes. Difference between j_lw=1 (call with forcing in
! place) and j_lw= 2 (call with reference values). N.B. Call to
! radiation schemes done first with j_lw=2 then with j_lw =1

          i_off=diagnostic_offset*(j_lw-1)

          IF (j_lw >  1) THEN

            IF (lw_diag(j_lw)%l_flux_up)                                &
                lw_diag(j_lw)%flux_up =                                 &
                lw_diag(1)%flux_up -                                    &
                lw_diag(j_lw)%flux_up

            IF (lw_diag(j_lw)%l_flux_down)                              &
                lw_diag(j_lw)%flux_down =                               &
                lw_diag(1)%flux_down -                                  &
                lw_diag(j_lw)%flux_down

            IF (lw_diag(j_lw)%l_flux_up_clear)                          &
                lw_diag(j_lw)%flux_up_clear =                           &
                lw_diag(1)%flux_up_clear -                              &
                lw_diag(j_lw)%flux_up_clear

            IF (lw_diag(j_lw)%l_flux_down_clear)                        &
                lw_diag(j_lw)%flux_down_clear =                         &
                lw_diag(1)%flux_down_clear -                            &
                lw_diag(j_lw)%flux_down_clear

            IF (lw_diag(j_lw)%l_clear_olr)                              &
                lw_diag(j_lw)%clear_olr =                               &
                lw_diag(1)%clear_olr -                                  &
                lw_diag(j_lw)%clear_olr

            IF (lw_diag(j_lw)%l_surf_down_clr)                          &
                lw_diag(j_lw)%surf_down_clr =                           &
                lw_diag(1)%surf_down_clr -                              &
                lw_diag(j_lw)%surf_down_clr

            IF (lw_diag(j_lw)%l_clear_hr)                               &
                lw_diag(j_lw)%clear_hr =                                &
                lw_diag(1)%clear_hr -                                   &
                lw_diag(j_lw)%clear_hr

            IF (lw_diag(j_lw)%l_net_flux_trop)                          &
                lw_diag(j_lw)%net_flux_trop =                           &
                lw_diag(1)%net_flux_trop -                              &
                lw_diag(j_lw)%net_flux_trop

            IF (lw_diag(j_lw)%l_down_flux_trop)                         &
                lw_diag(j_lw)%down_flux_trop =                          &
                lw_diag(1)%down_flux_trop -                             &
                lw_diag(j_lw)%down_flux_trop

          END IF
        ELSE IF (l_radiance) THEN
          IF (j_lw == 1) THEN
            i_off=0
          ELSE
            i_off=90+(diagnostic_offset/20)*(j_lw-1)
          END IF
        ELSE
          i_off=0
        END IF

! DEPENDS ON: diagnostics_lw
        CALL diagnostics_lw(                                            &
          row_length, rows, model_levels,                               &
          wet_model_levels, ozone_levels,                               &
          cloud_levels,                                                 &
          n_rows, global_row_length, global_rows,                       &
          halo_i, halo_j, off_x, off_y, me,                             &
          n_proc, n_procx, n_procy,                                     &
          g_rows, g_row_length,                                         &
          at_extremity,                                                 &
          timestep,i_off,                                               &
          t_n, t_inc,                                                   &
          q_n, qcl_n, cf_n, cfl_n,                                      &
          t_latest, q_latest, qcl_latest,                               &
          cf_latest, cfl_latest,                                        &
          surflw, olr, lw_down,                                         &
          t_incr_diagnostic,                                            &
          lw_spectrum(j_lw)%n_band,                                     &
          ozone,                                                        &
          o3_trop_level,                                                &
          o3_trop_height,                                               &
          t_trop_level,                                                 &
          t_trop_height,                                                &
          lw_incs, lwsea,                                               &
          lw_spectrum(j_lw)%n_aod_wavel,                                &
          j_lw,                                                         &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
          stashwork2 )

      END IF
    END DO


    DEALLOCATE ( t_incr_diagnostic )

! Deallocate the diagnostic space that is no longer required.
    IF (l_timestep) THEN
      IF (l_rad_step_diag) THEN
! DEPENDS ON: deallocate_swdiag_ext
        CALL deallocate_swdiag_ext(2)
! DEPENDS ON: deallocate_lwdiag
        CALL deallocate_lwdiag(2)
      END IF
      IF (timestep_number >  1) THEN
        IF (MOD(timestep_number,a_lw_radstep_prog) == 0) THEN

          DEALLOCATE(lw_incs_local)
          DEALLOCATE(lwsea_local)
          DEALLOCATE(olr_local)
          DEALLOCATE(lw_down_local)
          DEALLOCATE(top_abs_lw)

          CALL deallocate_swdiag_ext(1)
          CALL deallocate_lwdiag(1)

        END IF
      END IF
    ELSE IF (l_forcing) THEN
      IF (l_rad_step_prog) THEN
        CALL deallocate_swdiag_ext(1)
        CALL deallocate_lwdiag(1)
      END IF
      IF (l_rad_step_diag) THEN
        CALL deallocate_swdiag_ext(2)
        CALL deallocate_lwdiag(2)
      END IF
    ELSE IF (l_radiance) THEN
      IF (l_rad_step_prog) THEN
        DO j_sw = 1, n_swcall
          CALL deallocate_swdiag_ext(j_sw)
        END DO
        DO j_lw = 1, n_lwcall
          CALL deallocate_lwdiag(j_lw)
        END DO
      END IF
    ELSE
      IF (l_rad_step_prog) THEN
        CALL deallocate_swdiag_ext(1)
        CALL deallocate_lwdiag(1)
      END IF
    END IF

  END IF ! on error_code

! Restore the true tiled temperature in aggregate runs.
  IF (l_aggregate) tstar_tile(:,1)=tstar_tile_tmp

  9999 CONTINUE
! Check error condition
  IF (error_code /= 0) THEN
    CALL ereport(routinename, error_code, cmessage)
  END IF

  IF (lhook) CALL dr_hook('GLUE_RAD',zhook_out,zhook_handle)
END SUBROUTINE glue_rad
