! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine NI_conv_ctl


SUBROUTINE ni_conv_ctl (                                          &

! Parallel variables
  at_extremity, n_proc, me                                        &

! Parameters for iterative SISL
, numcycles, cycleno                                              &

! Model dimensions.
, row_length, rows                                                &
, rowsrowlength                                                   &
, model_levels, wet_levels, bl_levels, n_cca_lev                  &
, tr_levels, tr_vars, tr_ukca                                     &
, land_points                                                     &

! Model switches
, model_domain, l_calc_dxek, l_q_interact                         &
, l_mixing_ratio, l_mcr_qrain, l_mcr_qgraup, l_mcr_qcf2, l_ctile  &
, l_murk_source, l_dust, l_sulpc_so2, l_sulpc_dms, l_sulpc_nh3    &
, l_soot, l_ocff, l_biomass, l_nitrate, l_co2_interactive         &
, l_use_cariolle, l_flux_bc, l_spec_z0                            &

! in coordinate information
, z_rho, z_theta                                                  &

! in time stepping information.
, timestep, timestep_number                                       &

! SCM/Idealised UM
, flux_e,flux_h,z0m_scm,z0h_scm                                   &
! require for conv_diag call
, t_surf, zh, u_0_p, v_0_p                                        &
, ls_rain, ls_snow                                                &

! SCM diagnostics (dummy in full UM)
, nscmdpkgs, l_scmdiags                                           &

! in data fields.
, rho,  rho_only, rho_theta                                       &
, u_p, v_p, ustar_p, vstar_p, w, p, p_star, exner_rho_levels      &
, land_sea_mask, flandg, ice_fract                                &
, tstar_land, tstar_sea, tstar_sice, z0msea                       &
, p_layer_boundaries, p_layer_centres                             &
, exner_layer_boundaries, exner_layer_centres                     &
, t1_sd, q1_sd, exner_theta_levels                                &
, uw0_p, vw0_p, w_max, zlcl, zlcl_uv, ztop, dzh, entrain_coef     &
, conv_type                                                       &
, cumulus, l_shallow, l_congestus,l_congestus2,l_pc2_diag_sh_pts  &
, no_cumulus                                                      &
, ntml, ntpar                                                     &
, wstar, wthvs, delthvu, ql_ad, qsat_lcl, ftl, fqt                &
, shallowc, cu_over_orog, cape_undilute, cin_undilute             &
, deep_flag, past_precip, past_conv_ht                            &

! in Prognostics at beginning of time step (winds passed in higher up)
, theta_n, q_n, qcl_n, qcf_n, qrain_n, qgraup_n, qcf2_n           & 
, cf_liquid_n, cf_frozen_n, bulk_cf_n                             &
! in/out values with all increments added
, theta_star, q_star, qcl_star, qcf_star, qrain_star, qgraup_star &
, qcf2_star, cf_liquid_star                                       &
, cf_frozen_star, bulk_cf_star                                    &
! in/out Convection increments
, theta_inc, q_inc, qcl_inc, qcf_inc, cf_liquid_inc, cf_frozen_inc&
, bulk_cf_inc                                                     &
!  In  - total wind increments before convection on p grid
, r_u_p, r_v_p                                                    &
! in/out
, aerosol                                                         &
, dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6     &
, so2, so4_aitken, so4_accu, so4_diss                             &
, dms, nh3, soot_new, soot_agd, soot_cld, bmass_new, bmass_agd    &
, bmass_cld, ocff_new, ocff_agd, ocff_cld, nitr_acc, nitr_diss    &
, co2,  free_tracers, ukca_tracers                                &
, ozone_tracer                                                    &

! out fields
, cca0, ccw0, ccb0, cct0, cclwp0, lcbase0, cca0_2d, lctop, lcca   &
, cca,  ccw,  ccb,  cct,  cclwp,  lcbase,  cca_2d                 &
, conv_rain, conv_snow,  ddmfx                                    &

! error information
, error_code  )

!---------------------------------------------------------------------------
! Description: Interface to Atmospheric Physics convection code.

! method:
! Note all working arrays other than prognostics are defined without surface 
! level values for the ENDGame case. This matches what currently happens for
! new dynamics, i.e. the convection scheme only operates on model levels
! above the surface. All loops over model levels should always work from 1
! to model levels not 0.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: convection

! code description:
!   language: fortran 90
!   This code is written to umdp3 programming standards version 8.3.
!---------------------------------------------------------------------------

!$ USE omp_lib

USE dynamics_grid_mod, ONLY: l_vatpoles

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                          &
   udims, vdims, wdims, tdims, pdims, qdims,                              &
   tdims_s, pdims_s,                                                      &
   qdims_l, trdims_xstl

USE segments_mod, ONLY:                                                   &
  segment_type, meta_segment_type, segment_dims_type,                     &
  segments_mod_seg_meta, segments_mod_segments, segments_mod_seg_dims

! Model level heights from centre of Earth
USE level_heights_mod, ONLY: &
  r_theta_levels             &  ! Radii on theta levels (m)
 ,r_rho_levels                  ! Radii on rho levels (m)

USE cv_run_mod,  ONLY:                                            &
    i_convection_vn,                                              &
    i_convection_vn_0a,                                           &
    i_convection_vn_4a,                                           &
    i_convection_vn_5a,                                           &
    i_convection_vn_6a,                                           &
    l_fix_udfactor,          l_mom,                               &
    l_safe_conv,             l_rediagnosis,         ud_factor,    &
    iconv_shallow,           iconv_congestus,       iconv_deep,   &
    a_convect_segments,      n_conv_calls,                        &
    rad_cloud_decay_opt,     cld_life_opt,          cca_min,      &
    fixed_cld_life,          a_convect_seg_size,    l_murk_conv,  &
    l_dcpl_cld4pc2,          qmin_conv,             l_ccrad,      &
    l_3d_cca,                l_pc2_diag_sh,       l_conv_hist

USE cv_hist_constants_mod, ONLY:                                  &
    decay_period

USE cv_param_mod,  ONLY:                                          &
    cld_life_constant,       cld_life_func_hgt,     rad_decay_off,&
    rad_decay_full_timestep, rad_decay_conv_substep

! Convective diagnostic output arrays - mainly required here for SCM output
USE cv_diagnostic_array_mod, ONLY:                                      &
        precip_deep, precip_shall ,precip_mid ,precip_cong ,cape_out    &
       ,deep_ind ,shallow_ind ,congestus_ind ,congestus_ind2 ,mid_ind   &
       ,ntml_diag ,ntpar_diag ,freeze_diag ,kterm_diag                  &
       ,wstar_up_diag ,wstar_dn_diag ,mb1_diag ,mb2_diag ,wqt_cb        &
       ,wthetal_cb ,wqt_inv ,wthetal_inv ,sh_top ,sh_base ,cg_top       &
       ,cg_base ,cg_term ,cape_ts_diag ,ind_cape_reduced_diag           &
       ,conscav_dust ,conscav_so4ait ,conscav_so4acc ,conscav_so4dis    &
       ,conscav_agedsoot ,conscav_agedbmass ,conscav_agedocff           &
       ,conscav_nitracc ,conscav_nitrdiss ,conwash_so2 ,conwash_nh3     &
       ,uw_dp, vw_dp ,uw_shall ,vw_shall ,uw_mid ,vw_mid ,up_flux       &
       ,dwn_flux ,entrain_up ,detrain_up ,entrain_dwn ,detrain_dwn      &
       ,up_flux_half ,T_incr_diag_conv ,q_incr_diag_conv                &
       ,qcl_incr_diag_conv ,qcf_incr_diag_conv                          &
       ,cf_liquid_incr_diag_conv ,cf_frozen_incr_diag_conv              &
       ,bulk_cf_incr_diag_conv,u_incr_diag_conv ,v_incr_diag_conv       &
       ,theta_diag, q_diag                                              &
       ,t_incr_conv_only, q_incr_conv_only                              &
       ,mf_deep ,mf_congest ,mf_shall ,mf_midlev                        &
       ,dt_deep ,dt_congest ,dt_shall ,dt_midlev                        &
       ,dq_deep ,dq_congest ,dq_shall ,dq_midlev                        &
       ,du_deep ,du_congest ,du_shall ,du_midlev                        &
       ,dv_deep ,dv_congest ,dv_shall ,dv_midlev                        &
       ,wqt_flux_sh ,wql_flux_sh ,wthetal_flux_sh ,wthetav_flux_sh      &
       ,dubydt_pout ,dvbydt_pout ,conv_rain_3d ,conv_snow_3d            &
       ,qcl_incr_inhom_diag ,qcf_incr_inhom_diag                        &
       ,bulk_cf_incr_inhom_diag ,cf_liquid_incr_inhom_diag              &
       ,cf_frozen_incr_inhom_diag, deep_cfl_limited, mid_cfl_limited 

! Subroutines
USE conv_diag_4a_mod, ONLY: conv_diag_4a
USE conv_diag_5a_mod, ONLY: conv_diag_5a
USE conv_diag_6a_mod, ONLY: conv_diag_6a

USE glue_conv_4a_mod, ONLY: glue_conv_4a
USE glue_conv_5a_mod, ONLY: glue_conv_5a
USE glue_conv_6a_mod, ONLY: glue_conv_6a

USE cv_stash_flg_mod, ONLY:                                       &
  ! variables
  l_qcl_incr_cinh, l_qcf_incr_cinh, l_cfl_incr_cinh,              &
  l_cff_incr_cinh, l_bcf_incr_cinh, l_apply_diag,                 &
  l_T_conv_only, l_q_conv_only

USE UM_ParVars, ONLY: mype

USE turb_diff_ctl_mod, ONLY: delta_smag

USE dust_parameters_mod, ONLY: krain_dust, ksnow_dust,            &
  l_twobin_dust

USE ereport_mod, ONLY : Ereport

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE bl_option_mod, ONLY:                                          &
  isrfexcnvgust
USE cloud_inputs_mod, ONLY: i_pc2_conv_coupling,                  &
                            i_pc2_erosion_method, l_micro_eros,   &
                            dbsdtbs_turb_0, l_fixbug_pc2_mixph
USE pc2_constants_mod,ONLY: pc2eros_hybrid_allfaces, dbsdtbs_turb_1
USE ereport_mod, ONLY : ereport
USE Field_Types
USE UM_ParParams
USE PrintStatus_mod
USE c_bm_con_mod, ONLY: krain_agedbmass, ksnow_agedbmass
USE c_ocff_con_mod, ONLY: krain_agedocff, ksnow_agedocff
USE c_st_con_mod, ONLY: krain_agedsoot, ksnow_agedsoot
USE ncnwsh2_mod, ONLY: ncnwsh2
USE scnscv2_mod, ONLY: scnscv2
USE scnwsh2_mod, ONLY: scnwsh2
USE Submodel_Mod
IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments
! arguments with intent in. ie: input variables.
! Parallel setup variables
  INTEGER, INTENT(IN)  ::      &
    n_proc                     & ! Total number of processors
  , me                         & ! My processor number
  , numcycles                  & ! Number of cycles
  , cycleno                      ! cycle number

  LOGICAL, INTENT(IN) ::  &
    at_extremity(4)         ! Indicates if this processor is at north, south
                            ! east or west of the processor grid

! Model dimensions
  INTEGER, INTENT(IN)  ::      &
    row_length                 & ! Row length
  , rows                       & ! Number of rows
  , rowsrowlength              & ! rows*row_length
  , model_levels               & ! Number of model levels
  , wet_levels                 & ! Number of wet levels
  , bl_levels                  & ! Number of boundary layer levels
  , n_cca_lev                  & ! Number of levels for conv cloud
                                 ! amount: 1 for 2D, nlevs for 3D.
  , tr_levels                  & ! free tracer levels
  , tr_vars                    & ! number of free tracers
  , tr_ukca                    & ! number of ukca tracers
  , land_points                  ! No.of land points being processed, can be 0.

! Model switches
  INTEGER, INTENT(IN)  ::  &
    model_domain
  LOGICAL, INTENT(IN)  ::  &
    l_calc_dxek            & ! Switch for calculation of condensate increment
  , l_q_interact           & ! Switch allows overwriting of parcel variables
                             ! when calculating condensate increments.
  , l_mixing_ratio         & ! Use mixing ratio formulation
  , l_mcr_qrain            & ! True if prognostic rain
  , l_mcr_qgraup           & ! True if prognostic graupel
  , l_mcr_qcf2             & ! True if prognostic qcf2 2nd ice type
  , l_ctile                & ! True if coastal tiling selected.
  , l_murk_source          & ! Switch for murk source/scavanging
  , l_dust                 & ! Switch for mineral dust
  , l_sulpc_so2            & ! Switch for sulphur cycle
  , l_sulpc_dms            & ! Switch for DMS in S-cycle
  , l_sulpc_nh3            & ! Switch for NH3 in S-cycle
  , l_soot                 & ! Switch for soot cycle
  , l_biomass              & ! Switch for biomass aerosol scheme
  , l_ocff                 & ! Switch for fossil-fuel organic carbon scheme
  , l_nitrate              & ! Switch for ammonium nitrate aerosol
  , l_co2_interactive      & ! Switch for CO2 cycle
  , l_use_cariolle         & ! true if Cariolle ozone tracer scheme being used
  , l_flux_bc              & ! true if SCM using specified surface fluxes
  , l_spec_z0                ! true if roughness length has been specified
 
! physical constants

! Co-ordinate arrays  ! local vertical co-ordinate information
  REAL, INTENT(IN) ::                               &
    z_rho(pdims%i_end,pdims%j_end,pdims%k_end)      & ! rho levels
  , z_theta(tdims%i_end,tdims%j_end,tdims%k_end)      ! theta levels

! model time stepping info
  REAL, INTENT(IN) ::      &
    timestep                 ! Full model time step

  INTEGER, INTENT(IN)  ::      &
    timestep_number              ! Model timestep number

! Used in SCM for prescribed surface flux forcing (required for call to
!  conv_diag)
  REAL, INTENT(IN) ::           &
    flux_e(row_length, rows)    & ! Surface latent heat flux (W/m^2)
   ,flux_h(row_length, rows)    & ! Surface sensible heat flux (W/m^2)
   ,z0m_scm(row_length, rows)   & ! Fixed sea surface roughness lengths
   ,z0h_scm(row_length, rows)     ! for momentum and heat, used in SCM
! Required for call to conv_diag
  REAL, INTENT(INOUT) ::        &
    t_surf(row_length, rows)      ! surface temperature (K)

  REAL, INTENT(IN) ::           &
    zh(row_length, rows)        & ! boundary layer height (m)
  , u_0_p(row_length, rows)     & ! set to zero (surface current m/s)
  , v_0_p(row_length, rows)       ! set to zero (surface current m/s)

  REAL, INTENT(INOUT) ::       &
    ls_rain(row_length, rows)  &  ! Large-scale rainfall from section 4
  , ls_snow(row_length, rows)     ! Large-scale snowfall from section 4

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER, INTENT(IN)  ::  &
    nscmdpkgs                ! No of SCM diagnostics packages

  LOGICAL, INTENT(IN)  ::  &
    l_scmdiags(nscmdpkgs)    ! Logicals for SCM diagnostics packages

! Data fields
  REAL, INTENT(IN) ::                   &
    rho(pdims_s%i_start:pdims_s%i_end,  & !  denisty *r*r (kg/m)
        pdims_s%j_start:pdims_s%j_end,  & !
        pdims_s%k_start:pdims_s%k_end)  & !
  , rho_only(pdims%i_end,pdims%j_end,pdims%k_end)    & ! density (kg/m3)
  , rho_theta(tdims%i_end,tdims%j_end,tdims%k_end-1) & ! on th lev (kg/m3)
  , u_p(pdims%i_end,pdims%j_end,pdims%k_end)         & ! U (m/s) P-grid
  , v_p(pdims%i_end,pdims%j_end,pdims%k_end)         & ! V (m/s) P-grid
  , ustar_p(pdims%i_end,pdims%j_end,pdims%k_end)     & ! U+du (m/s) P-grid
  , vstar_p(pdims%i_end,pdims%j_end,pdims%k_end)     & ! V+dv (m/s) P-grid
  , w(wdims%i_start:wdims%i_end,        & ! W (m/s) (copy of prognostic array)
      wdims%j_start:wdims%j_end,        &  
                  0:wdims%k_end)        & ! Not exact match to module values

  , p(pdims_s%i_start:pdims_s%i_end,    & ! pressure (Pa)
      pdims_s%j_start:pdims_s%j_end,    &
      pdims_s%k_start:pdims_s%k_end)    &
  , p_star(row_length, rows)                         & ! surface pressure
  , exner_rho_levels(pdims_s%i_start:pdims_s%i_end,  & ! Exner on rho level
                     pdims_s%j_start:pdims_s%j_end,  & !
                     pdims_s%k_start:pdims_s%k_end)  &
  , exner_theta_levels(tdims_s%i_start:tdims_s%i_end,  & ! Exner on theta level
                       tdims_s%j_start:tdims_s%j_end,  & !
                       tdims_s%k_start:tdims_s%k_end)

  LOGICAL, INTENT(IN) ::             &
    land_sea_mask(row_length, rows)    ! land sea land

  REAL, INTENT(IN) ::                &
    flandg(pdims_s%i_start:pdims_s%i_end, & ! Land fraction of gridbox
         pdims_s%j_start:pdims_s%j_end)   & ! on all points
   ,ice_fract(row_length,rows)              ! fraction of sea that has ice

  REAL, INTENT(IN) :: &
    tstar_land(row_length, rows) &
   ,tstar_sea(row_length, rows)  & ! Surface T on sea
   ,tstar_sice(row_length, rows) & ! Surface T on sea-ice 
   ,z0msea(row_length, rows)       ! Sea surface roughness

  REAL, INTENT(IN) ::                                                    &
    p_layer_boundaries(pdims%i_end,pdims%j_end,0:pdims%k_end)            &
        ! pressure at layer boundaries. Same as p except at
        ! bottom level = pstar, and at top = 0.
  , p_layer_centres(tdims%i_end,tdims%j_end,0:tdims%k_end)               &
        ! pressure at layer centres. Same as p_theta_levels
        !except bottom level = pstar, and at top = 0.
  , exner_layer_boundaries(pdims%i_end,pdims%j_end,0:pdims%k_end)        &
  , exner_layer_centres(tdims%i_end,tdims%j_end,0:tdims%k_end)

  REAL,INTENT(IN) ::           &
    t1_sd(row_length, rows)    & ! set to zero initially
  , q1_sd(row_length, rows)      ! set to zero initially

  REAL,INTENT(IN) ::               &
    uw0_p(row_length, rows)        & ! X-COMP SURFACE STRESS ON P GRID
  , vw0_p(row_length, rows)        & ! Y-COMP SURFACE STRESS ON P GRID
  , w_max(row_length, rows)          ! max w in column

! Data arrays coming in from conv_diag & BL

  REAL,INTENT(INOUT) ::            &
    zlcl(row_length, rows)         & ! accurate lifting condensation level(m)
  , zlcl_uv(row_length, rows)      & ! LCL at level for uv calculation (m)
  , ztop(row_length, rows)         & ! Top of cloud layer (m).
  , dzh(row_length, rows)          & ! Inversion thickness (m).
  , entrain_coef(row_length, rows)   ! Entrainment coefficient

! Convective type array ::
  INTEGER, INTENT(INOUT) :: conv_type(row_length, rows)
                ! Integer index describing convective type:
                !    0=no convection
                !    1=non-precipitating shallow
                !    2=drizzling shallow
                !    3=warm congestus
                !    ...
                !    8=deep convection

  LOGICAL, INTENT(INOUT) ::          &
    cumulus(row_length, rows)        & ! Logical switch from boundary
                                       !    layer for presence of Cu
  , l_shallow(row_length, rows)      & ! Logical switch for shallow Cu
  , l_congestus(row_length, rows)    & ! Logical switch for congestus
  , l_congestus2(row_length, rows)   & ! congestus in descending air
  , l_pc2_diag_sh_pts(row_length,rows) ! Carry diagnostic shallow convective
                                       ! information for PC2
  LOGICAL, INTENT(IN) ::             &
    no_cumulus(row_length, rows)       ! Points which BL says cannot be
                                       ! cumulus

  INTEGER, INTENT(INOUT) ::   &
    ntml(row_length, rows)    & ! Top level of surface mixed layer
  , ntpar(row_length, rows)     ! Top level of initial parcel ascent

  REAL,INTENT(INOUT) ::            &
    wstar(row_length, rows)        & ! Mixed layer convective velocity scale
  , wthvs(row_length, rows)        & ! Surface flux of THV
  , delthvu(row_length, rows)      & ! Integral of undilute parcel
                                     ! buoyancy over convective cloud layer
  , ql_ad(row_length, rows)        & ! adiabatic liquid water content
                                     ! at inversion (kg/kg)
  , qsat_lcl(row_length, rows)       ! qsat at LCL (kg/kg)

  REAL,INTENT(IN) ::               &
    ftl(row_length, rows)          & ! Surface sensible heat flux from BL
                                     ! (W/m2) i.e. cp*rho*w'tl'
  , fqt(row_length, rows)          & ! Total water flux from surface
                                     ! (kg/m2/s) i.e. rho*w'qT'
  , shallowc(row_length,rows)      & ! Indicator set to 1.0 if shallow,
                                     !  0.0 if not shallow or not cumulus
  , cu_over_orog(row_length,rows)    ! Indicator for cumulus over steep
                                     ! orography. Indicator set to 1.0 if
                                     ! true, 0.0 if false. Exclusive.

  REAL,INTENT(INOUT) ::            &
    cape_undilute(row_length,rows) & ! Undilute CAPE from parcel ascent m2/s2
  , cin_undilute(row_length,rows)    ! Undilute CIN from parcel ascent m2/s2

! Past history of convection - only used if prognostics requested otherwise
! these fields should remoain unset as no space is allocated at the top
! level of the UM.

  REAL,INTENT(INOUT) ::          &
    deep_flag(row_length,rows)   &  ! Value between 0.0 and 1.0
                                    ! 1 if deep convection last timestep
   ,past_precip(row_length,rows) &  ! Rate of convective precip (kg/m2/s)
                                    ! decayed over time.
   ,past_conv_ht(row_length,rows)   ! Height of previous convection (m)


! IN/OUT stash diagnostics held between sub-steps

! Prognostics at start of time step
  REAL,INTENT(INOUT) ::                        &
    theta_n(tdims_s%i_start:tdims_s%i_end,     & ! theta (K)
            tdims_s%j_start:tdims_s%j_end,     &
            tdims_s%k_start:tdims_s%k_end)     &
  , q_n(qdims_l%i_start:qdims_l%i_end,         & ! q (kg/kg)
        qdims_l%j_start:qdims_l%j_end,         &
        qdims_l%k_start:qdims_l%k_end)         &
  , qcl_n(qdims_l%i_start:qdims_l%i_end,       & ! qcl (kg/kg)
          qdims_l%j_start:qdims_l%j_end,       &
          qdims_l%k_start:qdims_l%k_end)       &
  , qcf_n(qdims_l%i_start:qdims_l%i_end,       & ! qcf (kg/kg)
          qdims_l%j_start:qdims_l%j_end,       &
          qdims_l%k_start:qdims_l%k_end)       &
  , qrain_n(qdims_l%i_start:qdims_l%i_end,     & ! qrain (kg/kg)
            qdims_l%j_start:qdims_l%j_end,     &
            qdims_l%k_start:qdims_l%k_end)     &
  , qgraup_n(qdims_l%i_start:qdims_l%i_end,    & ! qraup (kg/kg)
             qdims_l%j_start:qdims_l%j_end,    &
             qdims_l%k_start:qdims_l%k_end)    &
  , qcf2_n(qdims_l%i_start:qdims_l%i_end,      & ! qcf2 (kg/kg)
           qdims_l%j_start:qdims_l%j_end,      &
           qdims_l%k_start:qdims_l%k_end)      &
  , cf_liquid_n(qdims_l%i_start:qdims_l%i_end, & ! liquid cloud fraction
                qdims_l%j_start:qdims_l%j_end, &
                qdims_l%k_start:qdims_l%k_end) &
  , cf_frozen_n(qdims_l%i_start:qdims_l%i_end, & ! ice cloud fraction
                qdims_l%j_start:qdims_l%j_end, &
                qdims_l%k_start:qdims_l%k_end) &
  , bulk_cf_n(qdims_l%i_start:qdims_l%i_end,   & ! total cloud fraction
              qdims_l%j_start:qdims_l%j_end,   &
              qdims_l%k_start:qdims_l%k_end)

! arguments with intent in/out. ie: input variables changed on output.
! '_star' = '_n' + all increments up to now in time step
! _star variables held on same levels as prognostics in the vertical for 
! ENDGame. Note Convection will not update any variables stored on level 0.
  REAL,INTENT(INOUT) ::                                                    &
    theta_star(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
                           1:tdims%k_end)                                  &
  , q_star(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,            &
                       1:qdims%k_end)                                      &
  , qcl_star(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,          &
                         1:qdims%k_end)                                    &
  , qcf_star(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,          &
                         1:qdims%k_end)                                    &
  , cf_liquid_star(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,    &
                               1:qdims%k_end)                              &
  , cf_frozen_star(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,    &
                               1:qdims%k_end)                              &
  , bulk_cf_star(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,      &
                 1:qdims%k_end)
! Not altered by convection  
  REAL,INTENT(IN) ::                                                       &
    qrain_star(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,        &
                       1:qdims%k_end)                                      &
  , qgraup_star(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,       &
                         1:qdims%k_end)                                    &
  , qcf2_star(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,         &
                         1:qdims%k_end)

! '_inc' - all increments to a field 
  REAL,INTENT(INOUT) ::                                                    &
    theta_inc(tdims%i_end,tdims%j_end,tdims%k_end)                         &
  , q_inc(tdims%i_end,tdims%j_end,qdims%k_end)                             &
  , qcl_inc(tdims%i_end,tdims%j_end,qdims%k_end)                           &
  , qcf_inc(tdims%i_end,tdims%j_end,qdims%k_end)                           &
  , cf_liquid_inc(tdims%i_end,tdims%j_end,qdims%k_end)                     &
  , cf_frozen_inc(tdims%i_end,tdims%j_end,qdims%k_end)                     &
  , bulk_cf_inc(tdims%i_end,tdims%j_end,qdims%k_end)

  REAL,INTENT(IN)  ::                  &
  r_u_p(pdims%i_start:pdims%i_end,     &  ! r_u (m/s)on P-grid
        pdims%j_start:pdims%j_end,     &
        pdims%k_start:pdims%k_end)     &
 ,r_v_p(pdims%i_start:pdims%i_end,     &  ! r_v (m/s) on P-grid
        pdims%j_start:pdims%j_end,     &
        pdims%k_start:pdims%k_end)

! Tracers
  REAL,INTENT(INOUT) ::                         &
    aerosol(tdims_s%i_start:tdims_s%i_end,      & ! aerosol tracer
            tdims_s%j_start:tdims_s%j_end,      &
            tdims_s%k_start:tdims_s%k_end)      &
  , dust_div1(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div1
              tdims_s%j_start:tdims_s%j_end,    &
              tdims_s%k_start:tdims_s%k_end)    &
  , dust_div2(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div2
              tdims_s%j_start:tdims_s%j_end,    &
              tdims_s%k_start:tdims_s%k_end)    &
  , dust_div3(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div3
              tdims_s%j_start:tdims_s%j_end,    &
              tdims_s%k_start:tdims_s%k_end)    &
  , dust_div4(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div4
              tdims_s%j_start:tdims_s%j_end,    &
              tdims_s%k_start:tdims_s%k_end)    &
  , dust_div5(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div5
              tdims_s%j_start:tdims_s%j_end,    &
              tdims_s%k_start:tdims_s%k_end)    &
  , dust_div6(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div6
              tdims_s%j_start:tdims_s%j_end,    &
              tdims_s%k_start:tdims_s%k_end)
                                                  !DIV1,DIV2,...,DIV6 DUST
  REAL,INTENT(INOUT) ::                         &
    so2(tdims_s%i_start:tdims_s%i_end,          & ! SO2
        tdims_s%j_start:tdims_s%j_end,          &
        tdims_s%k_start:tdims_s%k_end)          &
  , so4_aitken(tdims_s%i_start:tdims_s%i_end,   & ! SO4 aitken mode
               tdims_s%j_start:tdims_s%j_end,   &
               tdims_s%k_start:tdims_s%k_end)   &
  , so4_accu(  tdims_s%i_start:tdims_s%i_end,   & ! SO4 accumulation mode
               tdims_s%j_start:tdims_s%j_end,   &
               tdims_s%k_start:tdims_s%k_end)   &
  , so4_diss(tdims_s%i_start:tdims_s%i_end,     & ! SO4 dissipation mode
             tdims_s%j_start:tdims_s%j_end,     &
             tdims_s%k_start:tdims_s%k_end)     &
  , dms(tdims_s%i_start:tdims_s%i_end,          & ! DMS
        tdims_s%j_start:tdims_s%j_end,          &
        tdims_s%k_start:tdims_s%k_end)          &
  , nh3(tdims_s%i_start:tdims_s%i_end,          & ! NH3
        tdims_s%j_start:tdims_s%j_end,          &
        tdims_s%k_start:tdims_s%k_end)          &
  , soot_new(tdims_s%i_start:tdims_s%i_end,     & ! New soot
             tdims_s%j_start:tdims_s%j_end,     &
             tdims_s%k_start:tdims_s%k_end)     &
  , soot_agd(tdims_s%i_start:tdims_s%i_end,     & ! Aged soot
             tdims_s%j_start:tdims_s%j_end,     &
             tdims_s%k_start:tdims_s%k_end)     &
  , soot_cld(tdims_s%i_start:tdims_s%i_end,     & ! soot in cloud
             tdims_s%j_start:tdims_s%j_end,     &
             tdims_s%k_start:tdims_s%k_end)     &
  , bmass_new(tdims_s%i_start:tdims_s%i_end,    & ! New biomass
              tdims_s%j_start:tdims_s%j_end,    &
              tdims_s%k_start:tdims_s%k_end)    &
  , bmass_agd(tdims_s%i_start:tdims_s%i_end,    & ! Aged biomass
              tdims_s%j_start:tdims_s%j_end,    &
              tdims_s%k_start:tdims_s%k_end)    &
  , bmass_cld(tdims_s%i_start:tdims_s%i_end,    & ! Biomass in cloud
              tdims_s%j_start:tdims_s%j_end,    &
              tdims_s%k_start:tdims_s%k_end)    &
  , ocff_new(tdims_s%i_start:tdims_s%i_end,     & ! New ocff
             tdims_s%j_start:tdims_s%j_end,     &
             tdims_s%k_start:tdims_s%k_end)     &
  , ocff_agd(tdims_s%i_start:tdims_s%i_end,     & ! Aged ocff
             tdims_s%j_start:tdims_s%j_end,     &
             tdims_s%k_start:tdims_s%k_end)     &
  , ocff_cld(tdims_s%i_start:tdims_s%i_end,     & ! Cloud ocff
             tdims_s%j_start:tdims_s%j_end,     &
             tdims_s%k_start:tdims_s%k_end)     &
  , nitr_acc(tdims_s%i_start:tdims_s%i_end,     & ! nitrate accumulation mode
             tdims_s%j_start:tdims_s%j_end,     &
             tdims_s%k_start:tdims_s%k_end)     &
  , nitr_diss(tdims_s%i_start:tdims_s%i_end,    & ! nitrate dissipation mode
              tdims_s%j_start:tdims_s%j_end,    &
              tdims_s%k_start:tdims_s%k_end)    &
  , co2(tdims_s%i_start:tdims_s%i_end,          & ! CO2
        tdims_s%j_start:tdims_s%j_end,          &
        tdims_s%k_start:tdims_s%k_end)

  REAL,INTENT(INOUT) ::                                                  &
    free_tracers(tdims_s%i_start:tdims_s%i_end,                          &
                 tdims_s%j_start:tdims_s%j_end,                          &
                 trdims_xstl%k_start:trdims_xstl%k_end, tr_vars)         &
  , ukca_tracers(tdims_s%i_start:tdims_s%i_end,                          &
                 tdims_s%j_start:tdims_s%j_end,                          &
                 trdims_xstl%k_start:trdims_xstl%k_end, tr_ukca)         &
  , ozone_tracer( tdims_s%i_start:tdims_s%i_end,          & ! Ozone tracer
                  tdims_s%j_start:tdims_s%j_end,          &
                  tdims_s%k_start:tdims_s%k_end)

! arguments with intent out. ie: output variables.
REAL, INTENT(INOUT) ::                                            &
  conv_rain(row_length, rows)                                     &
, conv_snow(row_length, rows)

REAL, INTENT(OUT) :: ddmfx(row_length, rows)
!           ! Convective downdraught mass-flux at cloud-base

! Section 0 cloud properties:
!   These are passed to radiation and are subject to PC2/CCRad
!   and/or Cloud Property decay.

REAL, INTENT(INOUT) ::                  &
  cca0    (row_length,rows,n_cca_lev)   &! Cnv.Cld Amount (0-1)
, ccw0    (tdims%i_end,tdims%j_end,qdims%k_end)  &! Cnv.Cld Water (kg/kg)
, cclwp0  (row_length,rows)              ! Cloud Condensed water path
                                         ! (kg/m^2)

REAL, INTENT(INOUT) ::                  &! Cnv.Cld Amount (2d) with
  cca0_2d (row_length,rows)              ! no anvil
          ! Only used if: l_3d_cca = .False. .and. l_pc2 .and.
          ! l_pc2_diag_sh

INTEGER, INTENT(INOUT) ::    &
  ccb0    (row_length,rows)  &! Cnv Cld Base of highest layer on gridpoint
, cct0    (row_length,rows)  &! Cnv Cld Top  of highest layer on gridpoint
, lcbase0 (row_length,rows)   ! Cnv Cld Top  of lowest  layer on gridpoint

! Section 5 cloud properties:
!   Refer to diagnostics from calls to convection only.
!   Intent InOut because they are passed upto atmos_physics2
!   so they can be held until the last convection call before
!   exiting atmos_physics2


REAL, INTENT(INOUT) ::                &
  cca    (row_length,rows,n_cca_lev)  &! Cnv.Cld Amount (0-1)
, ccw    (tdims%i_end,tdims%j_end,qdims%k_end) &! Cnv.Cld Water (kg/kg)
, cclwp  (row_length,rows)            &! Cloud Condensed water path
                                       ! (kg/m^2)
, cca_2d(row_length,rows)             &! Cnv.Cld Amount (2d)
                                       ! with no anvil
, lcca   (row_length,rows)             ! Conv. Cloud Amount of
                                       ! lowest Convective Cloud Layer
                                       ! on gridpoint. With no anvil
                                       ! scheme modification (0-1)

INTEGER, INTENT(INOUT) ::  &!
  ccb    (row_length,rows) &! Cnv Cld Base of highest layer on gridpoint
, cct    (row_length,rows) &! Cnv Cld Top  of highest layer on gridpoint
, lcbase (row_length,rows) &! Cnv Cld Base of lowest  layer on gridpoint
, lctop  (row_length,rows)  ! Cnv Cld Top  of lowest  layer on gridpoint

INTEGER, INTENT(INOUT) ::                                         &
  error_code

!------------------------------------------------------------------------------
! End of subroutine argument definitions
!------------------------------------------------------------------------------

! Parameters for soot cycle tracer scavenging
! Parameters for biomass aerosol scavenging
! Parameters for fossil-fuel organic carbon aerosol scavenging
! Parameters for ammonium nitrate scavenging
! Start c_nitrate_con
!
! Convective scavenging coefficients for ammonium nitrate
!
!
      REAL, PARAMETER ::  krain_nitracc = 0.3E-4
      REAL, PARAMETER ::  ksnow_nitracc = 0.3E-4
      REAL, PARAMETER ::  krain_nitrdiss = 0.3E-4
      REAL, PARAMETER ::  ksnow_nitrdiss = 0.3E-4
!
! End c_nitrate_con

! Start blopt8a

! Description:
!   Permissible settings for BL options.


      INTEGER, PARAMETER :: off = 0  ! Switch disabled
      INTEGER, PARAMETER :: on  = 1  ! Switch enabled

!     Options for non-gradient stress following
      INTEGER, PARAMETER :: BrownGrant97 = 1
      INTEGER, PARAMETER :: BrownGrant97_limited = 2
!       Brown and Grant (1997), version 2 including a limit on its size

!     Options for flux gradient formulation
      INTEGER, PARAMETER :: Locketal2000   = 0
!       Flux gradients as in Lock et al. (2000)
      INTEGER, PARAMETER :: HoltBov1993 = 1
!       Flux gradients as in Lock et al (2000) but using
!       coefficients from Holtslag and Boville (1993)
      INTEGER, PARAMETER :: LockWhelan2006 = 2
!       Flux gradients as in Lock and Whelan (2006)

!     Options for entrainment enhancement in Sc over Cu
      INTEGER, PARAMETER :: Buoyrev_feedback = 1

!     Options for form drag
      INTEGER, PARAMETER :: No_drag         = 0
      INTEGER, PARAMETER :: Effective_z0    = 1
      INTEGER, PARAMETER :: Explicit_stress = 2

!     Options for marine boundary layers
      INTEGER, PARAMETER :: Fixed_Z0T = 0
!       Stanard flixed value of thermal roughness length over sea
      INTEGER, PARAMETER :: SurfDivZ0T = 1
!       Thermal roughness length over sea defined from surface
!       divergence theory
      INTEGER, PARAMETER :: DynDiag_ZL = 1
      INTEGER, PARAMETER :: DynDiag_ZL_corrn = 2
!       The ratio of the height of the inversion to the surface
!       Obukhov length is used as a dynamic criterion in the
!       diagnosis of BL types: version 2 includes changes to
!       cope with BL_LEVELS >> 3km
      INTEGER, PARAMETER :: DynDiag_ZL_CuOnly = 3
!       As 2 but only applied to points diagnosed with Cumulus 
!       and strictly for sea points (fland<0.01, cf 0.5)
      INTEGER, PARAMETER :: DynDiag_Ribased = 4
!       As 3 but also overrides Cumulus diagnosis if 
!          ZH(Ri) > ZLCL+zhloc_depth_fac*(ZHPAR-ZLCL)
!       Note that here Ri accounts for gradient adjustment by the 
!       non-local scheme.

!     Options for surface exchange
      INTEGER, PARAMETER :: Use_Correct_Ustar = 2
!       Option under the COR_MO_ITER switch for the dust scheme
!       to use the correct ustar
      INTEGER, PARAMETER :: Limit_ObukhovL = 3
!       Option under the COR_MO_ITER switch for imposing
!       lower limit on L in very stable conditions.
      INTEGER, PARAMETER :: Limit_expl_ustar = 2
!       Option under the COR_UST switch to limit the magnitude of the
!       explicitly calculated ustar
      INTEGER, PARAMETER :: IP_SrfExWithCnv = 1
!       Option to include deep convective gustiness in the surface
!       transfer

!     Options for convective boundary layers
      INTEGER, PARAMETER ::  UM_std     = 0
      INTEGER, PARAMETER ::  neut_cbl   = 1
      INTEGER, PARAMETER ::  LEM_conven = 2
      INTEGER, PARAMETER ::  LEM_std    = 3

!     Options for stable boundary layers
      INTEGER, PARAMETER ::  Long_tails           = 0
      INTEGER, PARAMETER ::  Sharpest             = 1
      INTEGER, PARAMETER ::  Sharp_sea_long_land  = 2
      INTEGER, PARAMETER ::  Mes_tails            = 3
      INTEGER, PARAMETER ::  Louis_tails          = 4
      INTEGER, PARAMETER ::  Depth_based          = 5
      INTEGER, PARAMETER ::  Sharp_sea_mes_land   = 6
      INTEGER, PARAMETER ::  LEM_stability        = 7
      INTEGER, PARAMETER ::  Sharp_sea_Louis_land = 8

!     Options for Prandtl number (in local Ri scheme)
      INTEGER, PARAMETER ::  Constant_SBL = 0
      INTEGER, PARAMETER ::  LockMailhot2004 = 1

! End blopt8a

!---------------------------------------------------------------------------
! local variables.
!---------------------------------------------------------------------------

! loop counters
INTEGER                                                           &
  i, j, k, ii, jj                                                 &
, call_number

INTEGER                                                           &
  n_conv_levels                                                   &
, n_conv_points                                                   &
, this_point                                                      &
, info                                                            &
, step                                                            &
, seg_points                                                      &
, first_point                                                     &
, kp1, km1

! CCRad - Local variables and parameters
REAL    :: cld_life_3d(row_length, rows, n_cca_lev)
REAL    :: decay_time
INTEGER :: ndecay_steps

! To hold mean profiles over convection substeps until they
! are decayed on either each substep or at the end of a set of
! convection substeps(n_conv_calls)
INTEGER ::                         &
  ccb0_local    (row_length,rows)  &
, cct0_local    (row_length,rows)  &
, lcbase0_local (row_length,rows)

REAL ::                                              &
  cca0_local   (row_length,rows,n_cca_lev)           &
, ccw0_local   (tdims%i_end,tdims%j_end,qdims%k_end) &
, cclwp0_local (row_length,rows)


INTEGER ::                       &
  freeze_lev(row_length,rows)    & ! freezing level no
, it_kterm_deep(row_length,rows) & !lev no for terminating deep con
, it_kterm_shall(row_length,rows) & !lev no for terminating shallow con
, it_cg_term(row_length,rows)      !lev no for terminating congestus

LOGICAL                          &
  it_mid_level(row_length,rows)    ! true if mid level convection


LOGICAL ::                       &
  cumulus_1d(rowsrowlength)      & ! 1-d version of CUMULUS
 ,l_shallow_1d(rowsrowlength)    & ! 1-d version of SHALLOW_BL
 ,l_congestus_1d(rowsrowlength)  & ! 1-d version of congestus
 ,l_mid_1d(rowsrowlength)        & ! 1-d version of mid
 ,l_no_cumulus_1d(rowsrowlength)   ! 1-d version of no cumulus


! Total increments applied up to now across model timestep
! ie.  d(slow physics) + d(dyn),
REAL ::                                                                &
  u_inc_step(tdims%i_end,tdims%j_end,tdims%k_end)                      &
 ,v_inc_step(tdims%i_end,tdims%j_end,tdims%k_end)

LOGICAL ::       &
  l_tracer       & ! Switch for tracer variables used
, l_full_zero      ! True if a dummy zero full field is required


! Allocatable arrays for PC2 scheme increment calculation - used to save
! memory when increment calculations are not requested.
REAL,DIMENSION(:,:,:),ALLOCATABLE::                               &
  t_earliest                                                      &
                   !  temperature at start of convection
, t_inc_latest                                                    &
                   !  temperature increment on input to PC2 homog
, q_earliest                                                      &
                   !  humidity    at start of convection
, qcl_earliest                                                    &
                   !  qCL         at start of convection
, cfl_earliest                                                    &
                   !  cf_liquid   at start of convection
, cff_earliest                                                    &
                   !  cf_frozen   at start of convection
, bcf_earliest                                                    &
                   !  bulk cloud  at start of convection
, theta_inc_pc2                                                   &
                   !  pot temperature increment due to PC2 homog
, q_inc_pc2                                                       &
                   !  humidity        increment due to PC2 homog
, qcl_inc_pc2                                                     &
                   !  qCL             increment due to PC2 homog
, cfl_inc_pc2                                                     &
                   !  cf_liquid       increment due to PC2 homog
, bcf_inc_pc2                                                     &
                   !  bulk cloud      increment due to PC2 homog
, bcf_above                                                       &
                   !  Bulk cloud fraction in layer above
, bcf_below
                   !  Bulk cloud fraction in layer below

REAL,DIMENSION(:,:,:),ALLOCATABLE::                               &
  full_zero        !  a dummy array for a zero field

! Local arrays holding information to be passed between physics
! routines.

! Diagnostics controlled by Diagnostic switches

! Convection
REAL ::                                           &
  dthbydt(tdims%i_end,tdims%j_end,qdims%k_end)    & ! theta increment
, dqbydt(tdims%i_end,tdims%j_end,qdims%k_end)     & ! q increment
, dqclbydt(tdims%i_end,tdims%j_end,qdims%k_end)   & ! Q4 Increment qCL
, dqcfbydt(tdims%i_end,tdims%j_end,qdims%k_end)   & ! Q4 Increment qCF
, dcflbydt(tdims%i_end,tdims%j_end,qdims%k_end)   & ! Cloud Increment
, dcffbydt(tdims%i_end,tdims%j_end,qdims%k_end)   & ! Cloud Increment
, dbcfbydt(tdims%i_end,tdims%j_end,qdims%k_end)   & ! Cloud Increment
, dubydt_p(tdims%i_end,tdims%j_end,tdims%k_end)   & ! du/dt on p grid
, dvbydt_p(tdims%i_end,tdims%j_end,tdims%k_end)     ! dv/dt on p grid

REAL ::                                           &
  dq_add(tdims%i_end,tdims%j_end,qdims%k_end)      ! q added for safe profile

!  allocatable array holding copies of all tracers
REAL, DIMENSION(:,:,:,:), ALLOCATABLE ::                          &
  tot_tracer


!  sizes for temp. tracer calculations
INTEGER ::                                                        &
  ntra_fld                                                        &
, ntra_tmp                                                        &
, ntra_lev

! Parameters for S Cycle tracer scavenging

REAL, PARAMETER ::                                                &
    krain_so4ait=0.3e-4                                           &
,   ksnow_so4ait=0.3e-4                                           &
,   krain_so4acc=0.3e-4                                           &
,   ksnow_so4acc=0.3e-4                                           &
,   krain_so4dis=0.3e-4                                           &
,   ksnow_so4dis=0.3e-4


! Holding arrays for Section 0 cloud properties around convection substeps
REAL ::                                  &
  it_cca0   (row_length,rows,n_cca_lev)  &! Cnv.Cld Amount (0-1)
, it_ccw0   (tdims%i_end,tdims%j_end,qdims%k_end) &! Cnv.Cld Water (kg/kg)
, it_cclwp0 (row_length,rows)            &! Cloud Condensed water path
                                          ! (kg/m^2)
, it_cca0_2d(row_length,rows)             ! Cnv.Cld Amount (2d) with
                                          ! no anvil

INTEGER ::                     &
  it_ccb0   (row_length,rows)  &! Cnv Cld Base of highest layer on gridpoint
, it_cct0   (row_length,rows)  &! Cnv Cld Top  of highest layer on gridpoint
, it_lcbase0(row_length,rows)   ! Cnv Cld Base of lowest  layer on gridpoint


! Holding arrays for Section 5 cloud properties around convection substeps
REAL ::                                  &
  it_cca    (row_length,rows,n_cca_lev)  &! Cnv.Cld Amount (0-1)
, it_ccw    (tdims%i_end,tdims%j_end,qdims%k_end) &! Cnv.Cld Water (kg/kg)
, it_cclwp  (row_length,rows)            &! Cloud Condensed water path
                                          ! (kg/m^2)
, it_cca_2d (row_length,rows)            &! Cnv.Cld Amount (2d)
                                          ! with no anvil
, it_lcca   (row_length,rows)             ! Conv. Cloud Amount of
                                          ! lowest Convective Cloud Layer
                                          ! on gridpoint. With no anvil
                                          ! scheme modification (0-1)

INTEGER ::                     &
  it_ccb    (row_length,rows)  &! Cnv Cld Base of highest layer on gridpoint
, it_cct    (row_length,rows)  &! Cnv Cld Top  of highest layer on gridpoint
, it_lcbase (row_length,rows)  &! Cnv Cld Base of lowest  layer on gridpoint
, it_lctop  (row_length,rows)   ! Cnv Cld Top  of lowest  layer on gridpoin


REAL  ::                                                          &
  it_conv_rain    (row_length, rows)                              &
, it_conv_snow    (row_length, rows)                              &
, it_conv_rain_3d (tdims%i_end,tdims%j_end,qdims%k_end)           &
, it_conv_snow_3d (tdims%i_end,tdims%j_end,qdims%k_end)           &
, it_cape_out(row_length, rows)                                   &
                                ! CAPE
, it_precip_dp(row_length, rows)                                  &
                                 ! deep precip
, it_precip_sh(row_length, rows)                                  &
                                 ! shallow precip
, it_precip_md(row_length, rows)                                  &
                                 ! mid-level precip
, it_precip_cg(row_length, rows)                                  &
                                 ! congestus precip
, it_wstar_dn(row_length, rows)                                   &
, it_wstar_up(row_length, rows)                                   &
, it_mb1(row_length, rows)                                        &
, it_mb2(row_length, rows)                                        &
, it_up_flux_half(pdims%i_end,pdims%j_end,pdims%k_end)            &
                                !up flux on half levs.
, ind_cape_reduced(row_length, rows)                              &
                             ! indicator of reduced cape timescale
, cape_ts_used(row_length, rows)                                  &
                             ! Actual cape timescale used for deep (s)
, it_ind_deep(row_length, rows)        & ! indicator deep actually occurs
, it_ind_shall(row_length, rows)       & ! indicator shallow actually occurs
, it_dp_cfl_limited(row_length, rows)  & ! Indicator of CFL limited deep
, it_md_cfl_limited(row_length, rows)    ! Indicator of CFL limited mid


LOGICAL :: l_mid(row_length, rows)
                          ! true if mid level convection is
                          ! possible on a convection substep

! u'w' and v'w' fluxes are on theta levels
REAL ::                                             &
  it_uw_dp(tdims%i_end,tdims%j_end,tdims%k_end)     &
, it_vw_dp(tdims%i_end,tdims%j_end,tdims%k_end)     &
, it_uw_shall(tdims%i_end,tdims%j_end,tdims%k_end)  &
, it_vw_shall(tdims%i_end,tdims%j_end,tdims%k_end)  &
, it_uw_mid(tdims%i_end,tdims%j_end,tdims%k_end)    &
, it_vw_mid(tdims%i_end,tdims%j_end,tdims%k_end)

! Mass flux, entrain etc in our model is calculated on theta levels
REAL ::                                                                &
  it_up_flux(tdims%i_end,tdims%j_end,tdims%k_end)                      &
, it_dwn_flux(tdims%i_end,tdims%j_end,tdims%k_end)                     &
, it_entrain_up(tdims%i_end,tdims%j_end,tdims%k_end)                   &
, it_detrain_up(tdims%i_end,tdims%j_end,tdims%k_end)                   &
, it_entrain_dwn(tdims%i_end,tdims%j_end,tdims%k_end)                  &
, it_detrain_dwn(tdims%i_end,tdims%j_end,tdims%k_end)                  &
, r_rho_lev(pdims%i_end,pdims%j_end,pdims%k_end)                       &
, r_theta_lev(tdims%i_end,tdims%j_end,0:tdims%k_end)                   &
, zh_copy(row_length, rows)                                            &
, one_over_conv_calls                                                  &
, fraction_step                                                        &
, timestep_conv

INTEGER ::                     &
  nlcl(row_length, rows)         ! dummy array for conv_diag call

!-----------------------------------------------------------------------

REAL ::                                                            &
  it_mf_deep(tdims%i_end,tdims%j_end,qdims%k_end)                  &
, it_mf_congest (tdims%i_end,tdims%j_end,qdims%k_end)              &
, it_mf_shall(tdims%i_end,tdims%j_end,qdims%k_end)                 &
, it_mf_midlev(tdims%i_end,tdims%j_end,qdims%k_end)                &
, it_dt_deep(tdims%i_end,tdims%j_end,qdims%k_end)                  &
, it_dt_congest(tdims%i_end,tdims%j_end,qdims%k_end)               &
, it_dt_shall(tdims%i_end,tdims%j_end,qdims%k_end)                 &
, it_dt_midlev(tdims%i_end,tdims%j_end,qdims%k_end)                &
, it_dq_deep(tdims%i_end,tdims%j_end,qdims%k_end)                  &
, it_dq_congest(tdims%i_end,tdims%j_end,qdims%k_end)               &
, it_dq_shall(tdims%i_end,tdims%j_end,qdims%k_end)                 &
, it_dq_midlev (tdims%i_end,tdims%j_end,qdims%k_end)               &
, it_du_deep (tdims%i_end,tdims%j_end,qdims%k_end)                 &
, it_du_congest(tdims%i_end,tdims%j_end,qdims%k_end)               &
, it_du_shall(tdims%i_end,tdims%j_end,qdims%k_end)                 &
, it_du_midlev (tdims%i_end,tdims%j_end,qdims%k_end)               &
, it_dv_deep (tdims%i_end,tdims%j_end,qdims%k_end)                 &
, it_dv_congest(tdims%i_end,tdims%j_end,qdims%k_end)               &
, it_dv_shall (tdims%i_end,tdims%j_end,qdims%k_end)                &
, it_dv_midlev (tdims%i_end,tdims%j_end,qdims%k_end)

! These fluxes are on rho levels but only for wet levels so cannot
! put pdims as vertical dimension as would upset non ENDGame
REAL ::                                                            &
 it_wqt_flux(pdims%i_end,pdims%j_end,qdims%k_end)                  &
,it_wql_flux(pdims%i_end,pdims%j_end,qdims%k_end)                  &
,it_wthetal_flux(pdims%i_end,pdims%j_end,qdims%k_end)              &
,it_wthetav_flux(pdims%i_end,pdims%j_end,qdims%k_end)

! The '_conv' arrays hold the prognostic values input to the convection scheme.
! The convection scheme requires all the fields on the p_grid in the
! horizontal. Values held in _conv arrays depend on the choice of sub-stepping
! method and are never the start of time step prognostic values held in the
! '_n' arrays.

! If l_rediagnosis = T: '_conv' = '_n' + (1/n_conv_calls)*(slow phys+dyn incs)
! If l_rediagnosis = F: '_conv' = '_n' + (slow phys+dyn incs) = '_star'

REAL ::                                                                 &
  theta_conv(tdims%i_end,tdims%j_end,tdims%k_end)                       &
 ,q_conv(tdims%i_end,tdims%j_end,qdims%k_end)                           &
 ,qcl_conv(tdims%i_end,tdims%j_end,qdims%k_end)                         &
 ,qcf_conv(tdims%i_end,tdims%j_end,qdims%k_end)                         &
 ,qrain_conv(tdims%i_end,tdims%j_end,qdims%k_end)                       &
 ,qgraup_conv(tdims%i_end,tdims%j_end,qdims%k_end)                      &
 ,qcf2_conv(tdims%i_end,tdims%j_end,qdims%k_end)                        &
 ,cf_liquid_conv(tdims%i_end,tdims%j_end,qdims%k_end)                   &
 ,cf_frozen_conv(tdims%i_end,tdims%j_end,qdims%k_end)                   &
 ,bulk_cf_conv(tdims%i_end,tdims%j_end,qdims%k_end)                     &
 ,u_conv(pdims%i_end,pdims%j_end,pdims%k_end)                           &
 ,v_conv(pdims%i_end,pdims%j_end,pdims%k_end)


!----------------------------------------------------------------------
! Arrays to hold total increments to T, q (start of time step values) etc
! before calling convection. These increments are from the slow physics
! (atmos_physics1) and from semi-lagrangian dynamics.
! Note total wind increments are already held in r_u_p and r_v_p arrays.
! Want to know all these increments so they can be applied gradually
! over sub-steps in convection. Assuming l_rediagnosis = .true.
!----------------------------------------------------------------------
REAL, ALLOCATABLE ::         &
  theta_inc_step(:,:,:)      &
 ,q_inc_step(:,:,:)          &
 ,qcl_inc_step(:,:,:)        &
 ,qcf_inc_step(:,:,:)        &
 ,qrain_inc_step(:,:,:)      &
 ,qgraup_inc_step(:,:,:)     &
 ,qcf2_inc_step(:,:,:)       &
 ,cf_liquid_inc_step(:,:,:)  &
 ,cf_frozen_inc_step(:,:,:)  &
 ,bulk_cf_inc_step(:,:,:)

! Segmentation variables
TYPE(segment_type),      ALLOCATABLE  :: segments(:)
TYPE(meta_segment_type)               :: meta_segments
TYPE(segment_dims_type)               :: segment_dims

INTEGER, ALLOCATABLE  :: n_cumulus(:)   ! Number of CUMULUS points
INTEGER, ALLOCATABLE  :: n_deep(:)      ! Number of DEEP points
INTEGER, ALLOCATABLE  :: n_shallow(:)   ! Number of SHALLOW points
INTEGER, ALLOCATABLE  :: n_congestus(:) ! Number of congestus points
INTEGER, ALLOCATABLE  :: n_mid(:)       ! Number of mid-level points

REAL, PARAMETER :: a_land=0.3
REAL, PARAMETER :: a_sea=0.3
REAL, PARAMETER :: b_land=0.025
REAL, PARAMETER :: b_sea=0.025
REAL, PARAMETER :: reduction_factor=0.75

! convective history
REAL ::              &
  decay_amount       &  ! decay fraction
, tot_conv_precip       ! total convective precip rate

! Several options are available:
INTEGER, PARAMETER :: pc2_conv_original = 1        
! As used in Wilson et al (2008). Condensate and cloud fraction
! increments are calculated independently and there is no checking that
! the implied in-cloud condensate amounts are sensible.
!
! i_pc2_conv_coupling /= pc2_conv_original 
! Protect against generation of inconsistently low cloud 
! fraction implying very high in-cloud condensate amounts.
! If the in-cloud condensate amount is about to be made greater than
! 2.0e-3 then the cloud fraction is increased (up to a value of 1.0)

REAL    :: tmpb4 ! Temporary storing space for cloud fraction before
                 ! making consistency checks.

!-----------------------------------------------------------------------

! Error reporting variables
INTEGER ::                                                        &
  errorstatus

CHARACTER (LEN=80) ::                                             &
  cmessage

CHARACTER (LEN=*), PARAMETER ::                                   &
  routinename='NI_conv_ctl'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!Segmentation
INTEGER :: num_parallel
INTEGER :: ipar

!=============================================================================
! Start of convection control routine
!=============================================================================

IF (lhook) CALL dr_hook('NI_CONV_CTL',zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
IF (error_code  ==  0 ) THEN     ! continue with convection

!-----------------------------------------------------------------------------
! Section 1 - initialisation and setup
!-----------------------------------------------------------------------------

!-------------------------------------
! 1.1 Initialise local cloud arrays
!-------------------------------------
  DO j=1, rows
    DO i=1, row_length
      cclwp0_local(i,j) = 0.0
      ccb0_local(i,j)   = 0
      cct0_local(i,j)   = 0
    END DO
  END DO

  DO k=1, n_cca_lev
    DO j=1, rows
      DO i=1, row_length
        cca0_local(i,j,k) = 0.0
      END DO
    END DO
  END DO

  DO k=1, qdims%k_end
    DO j=1, rows
      DO i=1, row_length
        ccw0_local(i,j,k) = 0.0
      END DO
    END DO
  END DO


!-----------------------------------------------------------------------
! 1.3 Convection time step information and level number
!-----------------------------------------------------------------------

  timestep_conv=timestep/(n_conv_calls*1.0)

! Sub-timestep scheme

  fraction_step = 1.0/(FLOAT(n_conv_calls))

  one_over_conv_calls = 1.0/(n_conv_calls*1.0)

! Set number of levels to call convection for

  n_conv_levels = qdims%k_end

  IF (l_mom) THEN
! Limit convection calling levels to maximum of model_levels - 1
! This is because CMT increments to u and v exist on n_levels + 1
    IF (n_conv_levels  >   model_levels - 1 ) THEN
      n_conv_levels = model_levels - 1
    END IF
  END IF

!-----------------------------------------------------------------------
! 1.4  Initialise arrays
!-----------------------------------------------------------------------
  l_full_zero = l_calc_dxek


! L_full_zero_if1:
  IF (l_full_zero) THEN

!  Set up a field of zero values if dummy field is required

    ALLOCATE ( full_zero(row_length,rows,1) )

    DO j=1, rows
      DO i=1, row_length
        full_zero(i,j,1)  = 0.0
      END DO  ! I
    END DO  ! J

  END IF  ! L_full_zero_if1

! Initialise arrays for convection - without halos

  DO k = 1, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        r_rho_lev(i,j,k) = r_rho_levels(i,j,k)
      END DO
    END DO
  END DO

  DO k = 0, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        r_theta_lev(i,j,k) = r_theta_levels(i,j,k)
      END DO
    END DO
  END DO
!---------------------------------------------------------------------
! Initialise flag to indicate whether mid-level convection is possible
! On first convection substep it is possible everywhere.
!---------------------------------------------------------------------
  DO j=1, rows
    DO i=1, row_length
      l_mid(i,j) = .TRUE.
    END DO  !i
  END DO    !j

!---------------------------------------------------------------------
! 1.5 Set values of theta_conv, q_conv etc for input to convection
!-----------------------------------------------------------------------
  IF (l_rediagnosis .AND. n_conv_calls > 1) THEN

! Allocate arrays to hold increments
    ALLOCATE ( theta_inc_step(tdims%i_end,tdims%j_end,tdims%k_end) )
    ALLOCATE ( q_inc_step(tdims%i_end,tdims%j_end,qdims%k_end) )
    ALLOCATE ( qcl_inc_step(tdims%i_end,tdims%j_end,qdims%k_end) )
    ALLOCATE ( qcf_inc_step(tdims%i_end,tdims%j_end,qdims%k_end) )
    ALLOCATE ( cf_liquid_inc_step(tdims%i_end,tdims%j_end,qdims%k_end) )
    ALLOCATE ( cf_frozen_inc_step(tdims%i_end,tdims%j_end,qdims%k_end) )
    ALLOCATE ( bulk_cf_inc_step(tdims%i_end,tdims%j_end,qdims%k_end) )
    IF (l_mcr_qrain) THEN
      ALLOCATE ( qrain_inc_step(tdims%i_end,tdims%j_end,qdims%k_end) )
    END IF
    IF (l_mcr_qgraup) THEN
      ALLOCATE ( qgraup_inc_step(tdims%i_end,tdims%j_end,qdims%k_end) )
    END IF
    IF (l_mcr_qcf2) THEN
      ALLOCATE ( qcf2_inc_step(tdims%i_end,tdims%j_end,qdims%k_end) )
    END IF


!-----------------------------------------------------------------------
! Fill arrays
!-----------------------------------------------------------------------
! Work out timestep increments so far i.e. from slow physics and
! semi-lagrangian dynamics. Note may be possible in future to pass
! these in directly from atm_step rather than pass in theta_star etc.
!-----------------------------------------------------------------------
    DO k=1,tdims%k_end
      DO j=1,rows
        DO i=1,row_length
          theta_inc_step(i,j,k) = theta_star(i,j,k)-theta_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

    DO k=1,qdims%k_end
      DO j=1,rows
        DO i=1,row_length
          q_inc_step(i,j,k)   = q_star(i,j,k)   - q_n(i,j,k)
          qcl_inc_step(i,j,k) = qcl_star(i,j,k) - qcl_n(i,j,k)
          qcf_inc_step(i,j,k) = qcf_star(i,j,k) - qcf_n(i,j,k)
          cf_liquid_inc_step(i,j,k) = cf_liquid_star(i,j,k)                  &
                                                 -cf_liquid_n(i,j,k)
          cf_frozen_inc_step(i,j,k) = cf_frozen_star(i,j,k)                  &
                                                 -cf_frozen_n(i,j,k)
          bulk_cf_inc_step(i,j,k)   = bulk_cf_star(i,j,k)                    &
                                                 -bulk_cf_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

    ! Extra prognostics if present
    IF (l_mcr_qrain) THEN
      DO k=1,qdims%k_end
        DO j=1,rows
          DO i=1,row_length
            qrain_inc_step(i,j,k)   = qrain_star(i,j,k)   - qrain_n(i,j,k)
            qrain_conv(i,j,k)   = qrain_n(i,j,k)                            &
                                   + fraction_step*qrain_inc_step(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
    END IF
    IF (l_mcr_qgraup) THEN
      DO k=1,qdims%k_end
        DO j=1,rows
          DO i=1,row_length
            qgraup_inc_step(i,j,k)   = qgraup_star(i,j,k) - qgraup_n(i,j,k)
            qgraup_conv(i,j,k)   = qgraup_n(i,j,k)                            &
                                   + fraction_step*qgraup_inc_step(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
    END IF
    IF (l_mcr_qcf2) THEN
      DO k=1,qdims%k_end
        DO j=1,rows
          DO i=1,row_length
            qcf2_inc_step(i,j,k)   = qcf2_star(i,j,k) - qcf2_n(i,j,k)
            qcf2_conv(i,j,k)   = qcf2_n(i,j,k)                            &
                                   + fraction_step*qcf2_inc_step(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
    END IF
     

!-----------------------------------------------------------------------
! Set up input prognostic values for convection - sub-step 1
!-----------------------------------------------------------------------

! We want theta_conv = theta_n(start_step)+fraction_step*theta_inc_step(i,j,k)

    DO k=1,tdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          theta_conv(i,j,k) = theta_n(i,j,k)                                &
                            +fraction_step*theta_inc_step(i,j,k)
        END DO ! i
      END DO ! j
    END DO

    DO k=1,qdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          q_conv(i,j,k)   = q_n(i,j,k)  + fraction_step*q_inc_step(i,j,k)

          qcl_conv(i,j,k) = qcl_n(i,j,k)                                    &
                                        + fraction_step*qcl_inc_step(i,j,k)
          qcf_conv(i,j,k) = qcf_n(i,j,k)                                    &
                                        + fraction_step*qcf_inc_step(i,j,k)
          cf_liquid_conv(i,j,k) = cf_liquid_n(i,j,k)                        &
                                      + fraction_step*cf_liquid_inc_step(i,j,k)
          cf_frozen_conv(i,j,k) = cf_frozen_n(i,j,k)                        &
                                      + fraction_step*cf_frozen_inc_step(i,j,k)
          bulk_cf_conv(i,j,k) = bulk_cf_n(i,j,k)                            &
                                      + fraction_step*bulk_cf_inc_step(i,j,k)
        END DO ! i
      END DO ! j
    END DO

  ELSE     ! non rediagnosis option
!-----------------------------------------------------------------------
! Original code uses all increments added to start of time step for
! sweep one
!-----------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k, i, j)

!$OMP DO SCHEDULE(STATIC)
    DO k=1,tdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          theta_conv(i,j,k) = theta_star(i,j,k)
        END DO ! i
      END DO ! j
    END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
    DO k=1,qdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          q_conv(i,j,k)   = q_star(i,j,k)
          qcl_conv(i,j,k) = qcl_star(i,j,k)
          qcf_conv(i,j,k) = qcf_star(i,j,k)
          cf_liquid_conv(i,j,k) = cf_liquid_star(i,j,k)
          cf_frozen_conv(i,j,k) = cf_frozen_star(i,j,k)
          bulk_cf_conv(i,j,k)   = bulk_cf_star(i,j,k)
        END DO ! i
      END DO ! j
    END DO
!$OMP END DO

!$OMP END PARALLEL
    ! Extra prognostics if present
    IF (l_mcr_qrain) THEN
      DO k=1,qdims%k_end
        DO j=1,rows
          DO i=1,row_length
            qrain_conv(i,j,k)   = qrain_star(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
    END IF
    IF (l_mcr_qgraup) THEN
      DO k=1,qdims%k_end
        DO j=1,rows
          DO i=1,row_length
            qgraup_conv(i,j,k)   = qgraup_star(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
    END IF
    IF (l_mcr_qcf2) THEN
      DO k=1,qdims%k_end
        DO j=1,rows
          DO i=1,row_length
            qcf2_conv(i,j,k)   = qcf2_star(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
    END IF
  
  END IF   ! test on l_rediagnosis

! ----------------------------------------------------------------------
! 1.6 Check initial input for convection ok if PC2
! ----------------------------------------------------------------------
! Prevent negative condensate problems by condensing vapour if needed.
! Input fields are updated without updating increments so change will
! only affect PC2 increments and not the top-level QX_STAR fields.
! ----------------------------------------------------------------------
! L_calc_dxek_if0:
  IF (l_calc_dxek) THEN

! DEPENDS ON: conv_pc2_init
    CALL conv_pc2_init(rows, row_length, tdims%k_end, qdims%k_end,      &
                       exner_theta_levels,                              &
                       theta_conv, q_conv, qcl_conv, qcf_conv,          &
                       cf_liquid_conv, cf_frozen_conv, bulk_cf_conv)

  END IF  ! L_calc_dxek_if0

! ----------------------------------------------------------------------
! 1.7 Check for negative (less than a minimum) q being passed to convection
! ----------------------------------------------------------------------
! Note above PC2 code checking no negative qcl and qcf.
! Nothing currently preventing negative q going to convection
! Check for negative q if opt for safer convection or rediagnosis
! NOTE - NOT checking extra prognistics for negative values!
!---------------------------------------------------------------------

  IF (l_safe_conv .OR. l_rediagnosis) THEN
    DO k=1,qdims%k_end
      DO j=1,rows
        DO i=1,row_length
         IF (q_conv(i,j,k) < qmin_conv) THEN
           IF (printstatus >= prstatus_normal) THEN
             ! Print statement assumes i,j,k < 100000 ok at present
             ! Formatting for floating point double precision to get the same
             ! as unformatted print  
             WRITE(6,'(a42,3i5,a8,g26.18,a5,g26.18,a6,i10,i3)')                    &
                   ' Negative q at start of convection: i,j,k ',i,j,k,   &
                   ' q_conv ', q_conv(i,j,k),' q_n ',q_n(i,j,k),' step ', &
                     timestep_number,mype
           END IF
           ! store q added to ensure sensible profile
           dq_add(i,j,k) = qmin_conv - q_conv(i,j,k)
           ! Reset to qmin in non-conservative way
           q_conv(i,j,k) = qmin_conv
         ELSE
           dq_add(i,j,k) = 0.0    
         END IF
       END DO ! i
      END DO ! j
    END DO ! k
  END IF  ! l_rediagnosis

!-----------------------------------------------------------------------
! 1.9 Convective momemtum transport - l_mom=.true.
!-----------------------------------------------------------------------
! U and V required by convection but on the p-grid.

! Convection uses wind at beginning of time step plus increments
! For new sub-stepping option require only part of increment to be
! added each sub-step

! Are U and V already on p grid in atmos_physics2? If in future this is
! true just pass down and only need to interpolate R_u etc
!----------------------------------------------------------------------

  IF (l_mom) THEN

    IF (l_rediagnosis) THEN

      DO k = 1,tdims%k_end
        DO j = 1, rows
          DO i = 1, row_length
            u_conv(i,j,k) = u_p(i,j,k)+fraction_step*r_u_p(i,j,k)
            v_conv(i,j,k) = v_p(i,j,k)+fraction_step*r_v_p(i,j,k)
            ! Need to initialise u_inc_step and v_inc_step
            u_inc_step(i,j,k) = r_u_p(i,j,k)
            v_inc_step(i,j,k) = r_v_p(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k

    ELSE       ! original code

!$OMP  PARALLEL DO DEFAULT(NONE) SHARED(u_conv, v_conv, ustar_p, tdims, &
!$OMP& vstar_p, rows, row_length) PRIVATE(i, j, k)
      DO k = 1, tdims%k_end
        DO j = 1, rows
          DO i = 1, row_length
            u_conv(i,j,k) = ustar_p(i,j,k)
            v_conv(i,j,k) = vstar_p(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
!$OMP END PARALLEL DO

    END IF   ! test on l_rediagnosis



IF (.NOT. l_vatpoles) THEN
      ! Need to set cumulus diagnosis to false at poles to ensure that
      ! not running CMT on points without valid winds and stress

      IF (model_domain  ==  mt_global ) THEN
        IF ( at_extremity(psouth) ) THEN
          DO i = 1, row_length
            cumulus(i,1) = .FALSE.
            l_shallow(i,1) = .FALSE.
            l_congestus(i,1) = .FALSE.
            l_congestus2(i,1) = .FALSE.
          END DO
        END IF
        IF ( at_extremity(pnorth) ) THEN
          DO i = 1, row_length
            cumulus(i,rows) = .FALSE.
            l_shallow(i,rows) = .FALSE.
            l_congestus(i,rows) = .FALSE.
            l_congestus2(i,rows) = .FALSE.
          END DO
        END IF
      END IF
ELSE
  ! V-AT-POLES - p grid not on poles so no need to reset cumulus etc

END IF ! vatpoles


  END IF   !(l_mom)

!-----------------------------------------------------------------------
! Section 1.10  Setup total tracer variables - Aerosol scheme or UKCA
!-----------------------------------------------------------------------
  l_tracer = ( ( cycleno == numcycles ) .AND.                &
               ( l_soot       .OR.  l_co2_interactive  .OR.  &
                 l_dust       .OR.  l_biomass          .OR.  &
                 l_sulpc_so2  .OR.  l_use_cariolle     .OR.  &
                 l_ocff       .OR.  l_nitrate          .OR.  &
                 l_murk_conv  .OR.  tr_ukca > 0        .OR.  &
                 tr_vars > 0 ) )

  ! work with tracers only in final cycle

  IF ( cycleno == numcycles ) THEN

! Set up array dimensions for total tracer array (free + sulphur cycle
! tracers) so that convective transport of all tracers is done

    IF ( (l_sulpc_so2 .OR. l_soot .OR. l_co2_interactive .OR.       &
          l_dust .OR. l_biomass .OR. l_use_cariolle .OR. l_ocff .OR.&
          l_nitrate)                                                &
          .AND. (tr_vars > 0) .AND.                                 &
          (trdims_xstl%k_end /= tdims%k_end ) )  THEN ! exit
       WRITE(cmessage,*) 'Cannot call convect for tracer expts with'&
                      ,' tr_levels /= model_levels'
       errorstatus = -10
      CALL ereport(routinename, errorstatus, cmessage )
    END IF

    IF ( l_sulpc_so2        .OR.  l_soot          .OR.  &
         l_dust             .OR.  l_biomass       .OR.  &
         l_co2_interactive  .OR.  l_use_cariolle  .OR.  &
         l_ocff             .OR.  l_nitrate       .OR.  &
         l_murk_conv ) THEN

      ntra_lev = tdims%k_end
      ntra_fld = 0        !Initialise to zero

      IF (l_dust) THEN
        ntra_fld = ntra_fld + 6    !ADD 6 DUST SIZE CLASSES
      END IF

      IF (l_sulpc_so2) THEN
        ntra_fld = ntra_fld + 4    !Add SO2 + 3 SO4 modes
        IF (l_sulpc_nh3) THEN
          ntra_fld = ntra_fld + 1  !Add NH3 field
        END IF
        IF (l_sulpc_dms) THEN
          ntra_fld = ntra_fld + 1  !Add DMS field
        END IF
      END IF

      IF (l_soot) THEN
        ntra_fld = ntra_fld + 3    !Add 3 modes of soot
      END IF

      IF (l_biomass) THEN
        ntra_fld = ntra_fld + 3    !Add 3 modes of biomass aerosol
      END IF

      IF (l_co2_interactive) THEN
        ntra_fld = ntra_fld + 1    !Add CO2 field
      END IF

      IF (l_ocff) THEN
        ntra_fld = ntra_fld + 3    !Add 3 modes of fossil-fuel OC
      END IF

      IF (l_nitrate) THEN
        ntra_fld = ntra_fld + 2    !Add 2 modes of ammonium nitrate
      END IF

      IF (l_use_cariolle) THEN
        ntra_fld = ntra_fld + 1    !Add Cariolle Ozone tracer field
      END IF

      IF (l_murk_conv) THEN
        ntra_fld = ntra_fld + 1    !ADD MURK aerosol
      END IF

      IF (tr_vars > 0) THEN
        ntra_fld = ntra_fld + tr_vars
      END IF
  
      IF (tr_ukca > 0) THEN
        ntra_fld = ntra_fld + tr_ukca
      END IF

    ELSE
      IF (tr_vars == 0 .AND. tr_ukca == 0) THEN
        ntra_fld = 1             ! can't have 0 sized arrays
        ntra_lev = 1
      ELSE
        ntra_fld = tr_vars + tr_ukca
        ntra_lev = trdims_xstl%k_end
      END IF
    END IF

    ! Allocate the space in tot_tracer
    ALLOCATE( tot_tracer(row_length, rows, ntra_lev, ntra_fld) )

    ! copy arrays into tot_tracer

    IF ( l_sulpc_so2        .OR.  l_soot          .OR.  &
         l_dust             .OR.  l_biomass       .OR.  &
         l_co2_interactive  .OR.  l_use_cariolle  .OR.  &
         l_ocff             .OR.  l_nitrate       .OR.  &
         l_murk_conv ) THEN

      ntra_tmp = 0

      IF (l_dust) THEN
        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            dust_div1(1:row_length, 1:rows, 1:model_levels)

        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            dust_div2(1:row_length, 1:rows, 1:model_levels)

        IF (.NOT. l_twobin_dust) THEN
          ntra_tmp = ntra_tmp + 1
          tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
              dust_div3(1:row_length, 1:rows, 1:model_levels)

          ntra_tmp = ntra_tmp + 1
          tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
              dust_div4(1:row_length, 1:rows, 1:model_levels)

          ntra_tmp = ntra_tmp + 1
          tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
              dust_div5(1:row_length, 1:rows, 1:model_levels)

          ntra_tmp = ntra_tmp + 1
          tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
              dust_div6(1:row_length, 1:rows, 1:model_levels)
        END IF
      END IF

      IF (l_sulpc_so2) THEN
        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            so2(1:row_length, 1:rows, 1:model_levels)

        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            so4_aitken(1:row_length, 1:rows, 1:model_levels)

        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            so4_accu(1:row_length, 1:rows, 1:model_levels)

        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            so4_diss(1:row_length, 1:rows, 1:model_levels)

        IF (l_sulpc_dms) THEN
          ntra_tmp = ntra_tmp + 1
          tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)=&
              dms(1:row_length, 1:rows, 1:model_levels)
        END IF

        IF (l_sulpc_nh3) THEN
          ntra_tmp = ntra_tmp + 1
          tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp)=&
              nh3(1:row_length, 1:rows, 1:model_levels)
        END IF
      END IF   ! l_sulpc_so2

      IF (l_soot) THEN
        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            soot_new(1:row_length, 1:rows, 1:model_levels)

        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            soot_agd(1:row_length, 1:rows, 1:model_levels)

        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            soot_cld(1:row_length, 1:rows, 1:model_levels)
      END IF

      IF (l_biomass) THEN
        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            bmass_new(1:row_length, 1:rows, 1:model_levels)

        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            bmass_agd(1:row_length, 1:rows, 1:model_levels)

        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            bmass_cld(1:row_length, 1:rows, 1:model_levels)
      END IF

      IF (l_co2_interactive) THEN
        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            co2(1:row_length, 1:rows, 1:model_levels)
      END IF

      IF (l_ocff) THEN
        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            ocff_new(1:row_length, 1:rows, 1:model_levels)

        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            ocff_agd(1:row_length, 1:rows, 1:model_levels)

        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            ocff_cld(1:row_length, 1:rows, 1:model_levels)
      END IF

      IF (l_nitrate) THEN
        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
           nitr_acc(1:row_length, 1:rows, 1:model_levels)

        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
           nitr_diss(1:row_length, 1:rows, 1:model_levels)
      END IF

      IF (l_use_cariolle) THEN
        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            ozone_tracer(1:row_length, 1:rows, 1:model_levels)
      END IF

      IF (l_murk_conv) THEN
        ntra_tmp = ntra_tmp + 1
        tot_tracer(1:row_length, 1:rows, 1:model_levels, ntra_tmp) = &
            aerosol(1:row_length, 1:rows, 1:model_levels)
      END IF

      IF ( tr_vars > 0 ) THEN
        tot_tracer(1:row_length, 1:rows, 1:trdims_xstl%k_end,                &
                   ntra_tmp+1:ntra_tmp+tr_vars) =                            &
          free_tracers(1:row_length, 1:rows, 1:trdims_xstl%k_end, 1:tr_vars) 

        ntra_tmp = ntra_tmp + tr_vars
      END IF

      IF ( tr_ukca > 0 ) THEN
        tot_tracer(1:row_length, 1:rows, 1:trdims_xstl%k_end,                &
                   ntra_tmp+1:ntra_tmp+tr_ukca) =                            &
           ukca_tracers(1:row_length, 1:rows, 1:trdims_xstl%k_end, 1:tr_ukca) 

        ntra_tmp = ntra_tmp + tr_ukca
      END IF

    ELSE    ! no soot, co2, sulphur, dust, biomass, OCFF, nitrate
            ! or murk
      ntra_tmp = 0
      IF ( tr_vars > 0 ) THEN
        tot_tracer  (1:row_length, 1:rows, 1:trdims_xstl%k_end, 1:tr_vars) = &
           free_tracers(1:row_length, 1:rows, 1:trdims_xstl%k_end, 1:tr_vars)
        ntra_tmp = tr_vars
      END IF

      IF ( tr_ukca > 0 ) THEN
        tot_tracer(1:row_length, 1:rows, 1:trdims_xstl%k_end,                &
                   ntra_tmp+1:ntra_tmp+tr_ukca) =                            &
           ukca_tracers(1:row_length, 1:rows, 1:trdims_xstl%k_end,1:tr_ukca)
        ntra_tmp = ntra_tmp + tr_ukca
      END IF

      IF ( tr_vars + tr_ukca == 0) THEN ! make sure the value is set.
        tot_tracer(:,:,:,:) = 0.0
      END IF
    END IF !soot or co2, sulphur, dust, biomass, cariolle o3, ocff, 
           !nitrate or murk

  ELSE
    ntra_fld = 1             ! can't have 0 sized arrays
    ntra_lev = 1
    ! Allocate dummy space in tot_tracer
    ALLOCATE( tot_tracer(row_length, rows, ntra_fld, ntra_lev) )
  END IF       ! test on cycleno == numcycles


!--------------------------------------------------------------------------
! 1.11 Choices based on whether first call to ni_conv_ctl this model physics
! timestep
!--------------------------------------------------------------------------

! At the last substep keep theta, q values as they
! will be used by the diagnostics subroutine

    DO k = 1, tdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          theta_diag(i,j,k) = theta_conv(i,j,k)
        END DO
      END DO
    END DO
    DO k = 1, qdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          q_diag(i,j,k) = q_conv(i,j,k)
        END DO
      END DO
    END DO

  !---------------------------------------------------------------------
  ! Set up cloud decay lifetime
  !---------------------------------------------------------------------

  IF (rad_cloud_decay_opt /= rad_decay_off) THEN

    SELECT CASE(rad_cloud_decay_opt)

      CASE(rad_decay_conv_substep)
        decay_time   = timestep_conv
        ndecay_steps = 1

      CASE(rad_decay_full_timestep)
        decay_time   = timestep
        ndecay_steps = n_conv_calls

    END SELECT

    SELECT CASE(cld_life_opt)
      CASE(cld_life_func_hgt)
        !---------------------------------------------------------------
        ! Make lifetime a function of cloud size by making it a function
        ! of height: 30 minutes for shallow, 2 hours for high cloud
        !---------------------------------------------------------------
        DO k=1, n_cca_lev
          DO j=1, rows
            DO i=1, row_length
              cld_life_3d(i,j,k) = decay_time / (1800.0                 &
                                 + (fixed_cld_life*0.5 - 900.0)         &
                                 * (TANH((z_theta(i,j,k)/1000.0)-5.0)   &
                                 + 1.0))
            END DO
          END DO
        END DO

      CASE(cld_life_constant)
        !---------------------------------------------------------------
        ! Set Cloud_lifetime to a constant
        !---------------------------------------------------------------
        DO k=1, n_cca_lev
          DO j=1, rows
            DO i=1, row_length
              cld_life_3d(i,j,k) = decay_time / fixed_cld_life
            END DO
          END DO
        END DO
    END SELECT    ! cld_life_opt

  END IF        !  Rad_cloud_decay_opt

!==========================================================================
! Section 2 Call Convection scheme.
!==========================================================================
!============================================================================
!  Convection sub-stepping loop
!============================================================================
  DO call_number = 1, n_conv_calls

    !--------------------------------------------------------------------------
    ! 2.1 initialisation of work arrays for convection substeps
    !-------------------------------------------------------------------------- 
    DO k=1, qdims%k_end
      DO j=1, rows
        DO i=1, row_length
          it_ccw          (i,j,k) = 0.0
          it_ccw0         (i,j,k) = 0.0
          it_conv_rain_3d (i,j,k) = 0.0
          it_conv_snow_3d (i,j,k) = 0.0
        END DO
      END DO
    END DO

    DO k=1, n_cca_lev
      DO j=1, rows
        DO i=1, row_length
          it_cca  (i,j,k) = 0.0
          it_cca0 (i,j,k) = 0.0
        END DO
      END DO
    END DO

    DO j = 1, rows
      DO i = 1, row_length
        it_lcca(i,j)   = 0.0
        it_lcbase(i,j) = 0
        it_lctop(i,j)  = 0

        it_ccb(i,j)    = 0
        it_cct(i,j)    = 0
        it_cca_2d(i,j) = 0.0
        it_cclwp(i,j)  = 0.0

        it_ccb0(i,j)   = 0
        it_cct0(i,j)   = 0
        it_cca0_2d(i,j)= 0.0
        it_cclwp0(i,j) = 0.0

        it_conv_rain(i,j) = 0.0
        it_conv_snow(i,j) = 0.0
        it_precip_dp(i,j) = 0.0
        it_precip_sh(i,j) = 0.0
        it_precip_md(i,j) = 0.0
        it_cape_out(i,j) = 0.0
        it_kterm_deep(i,j) = 0
        it_kterm_shall(i,j) = 0
        it_mid_level(i,j) = .FALSE.
        it_dp_cfl_limited(i,j) = 0.0
        it_md_cfl_limited(i,j) = 0.0
      END DO
    END DO

    IF (i_convection_vn == i_convection_vn_5a .OR.  &
        i_convection_vn == i_convection_vn_6a ) THEN
! Only used by 5A and 6A schemes
      DO j = 1, rows
        DO i = 1, row_length
          it_precip_cg(i,j) = 0.0
          it_wstar_up(i,j) = 0.0
          it_mb1(i,j) = 0.0
          it_mb2(i,j) = 0.0
          it_cg_term(i,j) = 0
        END DO
      END DO
    END IF
!-----------------------------------------------------------------------------
! NB: increments to t and q are added on inside routine.
! Increments to qCL, qCF, CFl, CFf are calculated but only added at the
! control level (see below).
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! 2.2 segment the call to convection, to reduce the memory required in
! the subroutine convect.
!-----------------------------------------------------------------------------


   !--------------------------------------------------------------------------
   ! 2.3  Rediagnosis ?
   !--------------------------------------------------------------------------

    IF (l_rediagnosis) THEN
    !-------------------------------------------------------------------------
    ! Do we want to rediagnose where shallow and deep convection occurs on
    ! subsequent convection sub-steps?
    !    call_number = 1  - conv_diag called from atmos_physics2 to
    !                       determine cumulus_1d, l_shallow etc
    !    call_number > 1  - conv_diag call here to rediagnose cumulus etc
    !-------------------------------------------------------------------------

      IF (call_number > 1) THEN

       ! Take copy as do not want to overwrite otherwise use
       ! value from previous sweep already held in zh_copy
        IF (call_number == 2) THEN
          DO j=1,rows
            DO i=1,row_length
              zh_copy(i,j) = zh(i,j)
            END DO
          END DO
        END IF
        !  Initialise conv_diag output arrays
        DO j = 1, rows
          DO i = 1, row_length
          ! ntml(i,j)=1     May want previous value
            ntpar(i,j)=1
            nlcl(i,j)=1
            cumulus(i,j)=.FALSE.
            l_shallow(i,j)=.FALSE.
            delthvu(i,j)=0.0
            ql_ad(i,j) = 0.0
            ztop(i,j)=0.0
            dzh(i,j) =-9.9E9     ! initialise to large and negative
            zlcl(i,j)=0.0
            zlcl_uv(i,j)=0.0
            wstar(i,j)=0.0
            wthvs(i,j)=0.0
          END DO
        END DO


        ! --------------------------------------------------------------------
        ! Prevent negative condensate problems by condensing vapour if needed.
        ! Input fields are updated without updating increments so change will
        ! only affect PC2 increments and not the top-level QX_STAR fields.
        ! --------------------------------------------------------------------
        ! L_calc_dxek_if0:
        IF (l_calc_dxek) THEN

! DEPENDS ON: conv_pc2_init
          CALL conv_pc2_init(rows, row_length, tdims%k_end, qdims%k_end, &
                    exner_theta_levels,                                  &
                    theta_conv, q_conv, qcl_conv, qcf_conv,              &
                    cf_liquid_conv, cf_frozen_conv, bulk_cf_conv)

        END IF  ! L_calc_dxek_if0
        !---------------------------------------------------------------------
        ! Check for negative q being passed to convection
        ! Note above PC2 code checking no negative qcl and qcf
        ! Nothing currently preventing negative q going to convection
        !---------------------------------------------------------------------
        DO k=1,qdims%k_end
          DO j=1,rows
            DO i=1,row_length
              IF (q_conv(i,j,k) < qmin_conv) THEN
                IF (printstatus >= prstatus_normal) THEN
                  ! Print statement assumes i,j,k < 1000000 ok at present
                  WRITE(6,'(a42,3i6,a8,g26.18,a5,g26.18,a8,g26.18,a5,g26.18)') &
                   ' Negative q at start of convection: i,j,k ',i,j,k,         &
                   ' q_conv ', q_conv(i,j,k),' q_n ',q_n(i,j,k),               &
                   ' dqbydt ', dqbydt(i,j,k),                                  &
                   ' inc ',fraction_step*q_inc_step(i,j,k)
                END IF
                ! store added q 
                dq_add(i,j,k) = dq_add(i,j,k) +(qmin_conv - q_conv(i,j,k))   
                ! Reset to qmin_conv in a non-conservative way
                q_conv(i,j,k) = qmin_conv
              END IF
            END DO ! i
          END DO ! j
        END DO ! k
!---------------------------------------------------------------------


        SELECT CASE ( i_convection_vn )
           CASE ( i_convection_vn_0a )

              ! Convection not used. Diagnosis is done in atmos_physics2
              ! Should not be possible to get here as ni_conv_ctl is not called
              ! when i_convection_vn_0a as L_PARAM_CONV=.FALSE. in atmos_phys2
              
             errorstatus = 10
             WRITE (cmessage,'(A,A)') 'Convection scheme 0a requested',    &
                   ' when L_PARAM_CONV=.TRUE. This is an invalid option.'
             CALL Ereport ( RoutineName, errorstatus, cmessage)
              
           CASE ( i_convection_vn_4a )

              CALL conv_diag_4a(                                        &
             
              ! IN Parallel variables
               row_length, rows                                         &
             
              ! IN model dimensions.
              , bl_levels, tdims%k_end, qdims%k_end                     &
              , land_points                                             &
              ,p, p_layer_centres(1,1,1)                                &
              ,exner_rho_levels                                         &
              ,rho_only, rho_theta, z_theta, z_rho                      &
             
              ! IN Model switches
              ,l_mixing_ratio, l_ctile, .TRUE.                          &
              ,no_cumulus                                               &
              ! IN cloud data
              ,qcf_conv, qcl_conv, bulk_cf_conv                         &
             
              ! IN everything not covered so far :
              ,p_star, q_conv, theta_conv, exner_layer_centres(:,:,1:)  &
              ,u_conv, v_conv, u_0_p, v_0_p                             &
              ,tstar_land, tstar_sea, tstar_sice, z0msea                &
              ,l_flux_bc, flux_e, flux_h, l_spec_z0, z0m_scm, z0h_scm   &
              , t_surf, land_sea_mask, flandg, ice_fract, timestep      &
              , w, w_max, deep_flag, past_precip, past_conv_ht          &
             
              ! SCM Diagnostics (dummy values in full UM)
              , nscmdpkgs,l_scmdiags                                    &
             
              ! INOUT data required elsewhere in UM system :
              ,zh_copy,ztop,dzh,zlcl,zlcl_uv,delthvu,ql_ad              &
              ,ntml,ntpar,nlcl,cumulus                                  &
              ,l_shallow,l_congestus,l_congestus2, conv_type            &
              ,cin_undilute,cape_undilute, wstar, wthvs                 &
              ,entrain_coef, qsat_lcl                                   &
              ,error_code                                               &
                 )

           CASE ( i_convection_vn_5a )

              CALL conv_diag_5a(                                        &
             
              ! IN Parallel variables
               row_length, rows                                         &
             
              ! IN model dimensions.
              , bl_levels, tdims%k_end, qdims%k_end                     &
              , land_points                                             &
              ,p, p_layer_centres(1,1,1)                                &
              ,exner_rho_levels                                         &
              ,rho_only, rho_theta, z_theta, z_rho                      &
             
              ! IN Model switches
              ,l_mixing_ratio, l_ctile, .TRUE.                          &
              ,no_cumulus                                               &
              ! IN cloud data
              ,qcf_conv, qcl_conv, bulk_cf_conv                         &
             
              ! IN everything not covered so far :
              ,p_star, q_conv, theta_conv, exner_layer_centres(:,:,1:)  &
              ,u_conv, v_conv, u_0_p, v_0_p                             &
              ,tstar_land, tstar_sea, tstar_sice, z0msea                &
              ,l_flux_bc, flux_e, flux_h, l_spec_z0, z0m_scm, z0h_scm   &
              , t_surf, land_sea_mask, flandg, ice_fract, timestep      &
              , w, w_max, deep_flag, past_precip, past_conv_ht          &
             
              ! SCM Diagnostics (dummy values in full UM)
              , nscmdpkgs,l_scmdiags                                    &
             
              ! INOUT data required elsewhere in UM system :
              ,zh_copy,ztop,dzh,zlcl,zlcl_uv,delthvu,ql_ad              &
              ,ntml,ntpar,nlcl,cumulus                                  &
              ,l_shallow,l_congestus,l_congestus2, conv_type            &
              ,cin_undilute,cape_undilute, wstar, wthvs                 &
              ,entrain_coef, qsat_lcl                                   &
              ,error_code                                               &
                 )
           CASE ( i_convection_vn_6a )

              CALL conv_diag_6a(                                        &
             
              ! IN Parallel variables
               row_length, rows                                         &
             
              ! IN model dimensions.
              , bl_levels, tdims%k_end, qdims%k_end                     &
              , land_points                                             &
              ,p, p_layer_centres(1,1,1)                                &
              ,exner_rho_levels                                         &
              ,rho_only, rho_theta, z_theta, z_rho                      &
             
              ! IN Model switches
              ,l_mixing_ratio, l_ctile, .TRUE.                          &
              ,no_cumulus                                               &
              ! IN cloud data
              ,qcf_conv, qcl_conv, bulk_cf_conv                         &
             
              ! IN everything not covered so far :
              ,p_star, q_conv, theta_conv, exner_layer_centres(:,:,1:)  &
              ,u_conv, v_conv, u_0_p, v_0_p                             &
              ,tstar_land, tstar_sea, tstar_sice, z0msea                &
              ,l_flux_bc, flux_e, flux_h, l_spec_z0, z0m_scm, z0h_scm   &
              , t_surf, land_sea_mask, flandg, ice_fract, timestep      &
              , w, w_max, deep_flag, past_precip, past_conv_ht          &
             
              ! SCM Diagnostics (dummy values in full UM)
              , nscmdpkgs,l_scmdiags                                    &
             
              ! INOUT data required elsewhere in UM system :
              ,zh_copy,ztop,dzh,zlcl,zlcl_uv,delthvu,ql_ad              &
              ,ntml,ntpar,nlcl,cumulus                                  &
              ,l_shallow,l_congestus,l_congestus2, conv_type            &
              ,cin_undilute,cape_undilute, wstar, wthvs                 &
              ,entrain_coef, qsat_lcl                                   &
              ,error_code                                               &
                 )

         CASE DEFAULT ! i_convection_vn

         errorstatus = 10
         WRITE (cmessage,'(A,A,I6)') 'Convection scheme version value',  &
                                 'i_convection_vn = ',i_convection_vn
         CALL Ereport ( RoutineName, errorstatus, cmessage)

        END SELECT ! i_convection_vn

        ! need to set cumulus diagnosis to false at poles to ensure that
        ! not running cmt on points without valid winds and stress

IF (.NOT. l_vatpoles) THEN
          IF (model_domain  ==  mt_global) THEN
            IF ( at_extremity(psouth) ) THEN
              DO i = 1, row_length
                cumulus(i,1) = .FALSE.
                l_shallow(i,1) = .FALSE.
                l_congestus(i,1) = .FALSE.
                l_congestus2(i,1) = .FALSE.
              END DO
            END IF
            IF ( at_extremity(pnorth) ) THEN
              DO i = 1, row_length
                cumulus(i,rows) = .FALSE.
                l_shallow(i,rows) = .FALSE.
                l_congestus(i,rows) = .FALSE.
                l_congestus2(i,rows) = .FALSE.
              END DO
            END IF
          END IF
ELSE
        ! No pole points on p grid for V-AT-POLES
END IF ! vatpoles

        ! Reset l_mid(i,j) = .true. as mid can occur anywhere on
        ! rediagnosis
        DO j = 1, rows
          DO i = 1, row_length
            l_mid(i,j) = .TRUE.
          END DO
        END DO

      END IF  ! Call_number > 1
    END IF   ! test on l_rediagnosis

!---------------------------------------------------------------------------
! 2.4 WORK OUT NUMBER OF CUMULUS POINTS IN EACH CONVECTION SEGMENT
!---------------------------------------------------------------------------
! Also calculate number of shallow and deep points in each convection
! segment

    ii=1
    DO j = 1, rows
      DO i = 1, row_length
        cumulus_1d(ii) = cumulus(i,j)
        l_shallow_1d(ii) = l_shallow(i,j)
        l_congestus_1d(ii) = l_congestus(i,j)
        l_mid_1d(ii) = l_mid(i,j)
        l_no_cumulus_1d(ii) = no_cumulus(i,j)
        ii=ii+1
      END DO
    END DO

!===============================================================================
! Start of segmented region
!===============================================================================

!$OMP  PARALLEL  DEFAULT(SHARED)                           &
!$OMP& PRIVATE(ii,jj,i,n_cumulus,n_deep,n_shallow,         &
!$OMP& n_congestus,n_mid,segments,meta_segments, segment_dims, ipar) 

!$OMP SINGLE
!Determine the number of threads in the parallel region
num_parallel=1
!$ num_parallel = omp_get_num_threads()
!$OMP END SINGLE
!Implicit barrier

!Determine the team member number.
!NB: this is numbered from 1 to follow the Fortran convention.
ipar=1
!$ ipar = omp_get_thread_num()+1

!Set up segment meta-information
CALL segments_mod_seg_meta(meta_segments, ipar, num_parallel,      & 
  rows*row_length, a_convect_seg_size, a_convect_segments)

!Set up array dimension information
CALL segments_mod_seg_dims(segment_dims, row_length, rows,         &
  model_levels, wet_levels)

! Allocate space for segmentation arrays
ALLOCATE(    segments( meta_segments%num_segments ) )
ALLOCATE(   n_cumulus( meta_segments%num_segments ) )
ALLOCATE(      n_deep( meta_segments%num_segments ) )
ALLOCATE(   n_shallow( meta_segments%num_segments ) )
ALLOCATE( n_congestus( meta_segments%num_segments ) )
ALLOCATE(       n_mid( meta_segments%num_segments ) )

! Work out starting points and sizes of each segment individually
CALL segments_mod_segments(segments, meta_segments, row_length*rows, &
  row_length, rows)

    DO i=1, meta_segments%num_segments
      n_cumulus(i)   = 0
      n_deep(i)      = 0
      n_shallow(i)   = 0
      n_congestus(i) = 0
      n_mid(i)       = 0

      DO j=segments(i)%fp, segments(i)%fp+segments(i)%seg_points-1
        IF(cumulus_1d(j)) THEN
          n_cumulus(i)=n_cumulus(i)+1
          IF (iconv_deep >  0.AND..NOT. l_shallow_1d(j).AND.      &
                      .NOT. l_congestus_1d(j)) THEN
            n_deep(i) = n_deep(i)+1
          END IF
          IF (iconv_shallow >  0.AND.l_shallow_1d(j).AND.         &
               .NOT. l_congestus_1d(j)) THEN
            n_shallow(i) = n_shallow(i)+1
          END IF
          IF (iconv_congestus >  0) THEN
            IF (l_congestus_1d(j)) THEN
              n_congestus(i) = n_congestus(i)+1
            END IF
          ELSE
            n_congestus(i) = 1     ! may be required for dim
          END IF
        END IF
      END DO

      DO j=segments(i)%fp, segments(i)%fp+segments(i)%seg_points-1
        IF (l_mid_1d(j)) THEN
          n_mid(i)=n_mid(i)+1
        END IF
      END DO

    END DO   ! i (segments)

!---------------------------------------------------------------------------
!  2.5 Loop over convection segments calling convection
!---------------------------------------------------------------------------

    DO i=1, meta_segments%num_segments

      ii = segments(i)%first_x
      jj = segments(i)%first_y

      SELECT CASE ( i_convection_vn )
      
         CASE ( i_convection_vn_0a )
      
           ! Convection not used. 
           ! Should not be possible to get here as ni_conv_ctl is not called
           ! when i_convection_vn_0a as L_PARAM_CONV=.FALSE. in atmos_phys2
           
           errorstatus = 10
           WRITE (cmessage,'(A,A)') 'Convection scheme 0a requested',     &
                   ' when L_PARAM_CONV=.TRUE. This is an invalid option.'
           CALL Ereport ( RoutineName, errorstatus, cmessage)           
      
         CASE ( i_convection_vn_4a )

            CALL glue_conv_4a                                                 &
              ( rows*row_length, segments(i)%seg_points                       &
              , n_conv_levels, bl_levels                                      &
              , call_number, i                                                &
              , theta_conv(ii,jj,1), q_conv(ii,jj,1)                          &
              , qcl_conv(ii,jj,1), qcf_conv(ii,jj,1)                          &
              , qrain_conv(ii,jj,1), qgraup_conv(ii,jj,1), qcf2_conv(ii,jj,1) & 
              , cf_liquid_conv(ii,jj,1), cf_frozen_conv(ii,jj,1)              &
              , bulk_cf_conv(ii,jj,1)                                         &
              , p_star(ii,jj), land_sea_mask(ii,jj)                           &
              , u_conv(ii,jj,1), v_conv(ii,jj,1), w(ii,jj,1)                  &
              , tot_tracer(ii,jj,1,1), dthbydt(ii,jj,1)                       &
              , dqbydt(ii,jj,1),   dqclbydt(ii,jj,1), dqcfbydt(ii,jj,1)       &
              , dcflbydt(ii,jj,1), dcffbydt(ii,jj,1), dbcfbydt(ii,jj,1)       &
              , dubydt_p(ii,jj,1), dvbydt_p(ii,jj,1)                          &
              , it_conv_rain(ii,jj),      it_conv_snow(ii,jj)                 &
              , it_conv_rain_3d(ii,jj,1), it_conv_snow_3d(ii,jj,1)            &
              , it_cca0(ii,jj,1),  it_ccb0(ii,jj),   it_cct0(ii,jj)           &
              , it_cclwp0(ii,jj),  it_ccw0(ii,jj,1), it_lcbase0(ii,jj)        &
              , it_cca0_2d(ii,jj), it_lctop(ii,jj),  it_lcca(ii,jj)           &
              , it_cca(ii,jj,1),   it_ccb(ii,jj),    it_cct(ii,jj)            &
              , it_cclwp(ii,jj),   it_ccw(ii,jj,1),  it_lcbase(ii,jj)         &
              , it_cca_2d(ii,jj), freeze_lev(ii,jj), it_dp_cfl_limited(ii,jj) &
              , it_md_cfl_limited(ii,jj)                                      &
              , it_mid_level(ii,jj), it_kterm_deep(ii,jj)                     &
              , it_kterm_shall(ii,jj)                                         &
              , it_precip_dp(ii,jj), it_precip_sh(ii,jj)                      &
              , it_precip_md(ii,jj), it_precip_cg(ii,jj)                      &
              , it_wstar_dn(ii,jj),  it_wstar_up(ii,jj)                       &
              , it_mb1(ii,jj), it_mb2(ii,jj), it_cg_term(ii,jj), n_cumulus(i) &
              , uw0_p(ii,jj), vw0_p(ii,jj), w_max(ii,jj)                      &
              , zlcl(ii,jj), zlcl_uv(ii,jj), ztop(ii,jj), entrain_coef(ii,jj) &
              , deep_flag(ii,jj), past_precip(ii,jj), past_conv_ht(ii,jj)     &
              , it_cape_out(ii,jj), n_deep(i), n_congestus(i), n_shallow(i)   &
              , n_mid(i), r_rho_lev(ii,jj,1), r_theta_lev(ii,jj,0)            &
              , rho_only(ii,jj,1), rho_theta(ii,jj,1), delta_smag(ii,jj)      &
              , exner_layer_boundaries(ii,jj,0), exner_layer_centres(ii,jj,0) &
              , p_layer_boundaries(ii,jj,0), p_layer_centres(ii,jj,0)         &
              , z_theta(ii,jj,1), z_rho(ii,jj,1), timestep_conv               &
              , t1_sd(ii,jj), q1_sd(ii,jj), ntml(ii,jj), ntpar(ii,jj)         &
              , conv_type(ii,jj), l_shallow(ii,jj), l_pc2_diag_sh_pts(ii,jj)  &
              , l_congestus(ii,jj), l_mid(ii,jj), cumulus(ii,jj)              &
              , wstar(ii,jj), wthvs(ii,jj), delthvu(ii,jj), ql_ad(ii,jj)      &
              , qsat_lcl(ii,jj), ftl(ii,jj), fqt(ii,jj), l_tracer, ntra_fld   &
              , ntra_lev, n_cca_lev, l_mixing_ratio, l_mcr_qrain              &
              , l_mcr_qgraup, l_mcr_qcf2, l_calc_dxek                         &
              , l_q_interact, it_up_flux_half(ii,jj,1)                        &
              , it_up_flux(ii,jj,1),      it_dwn_flux(ii,jj,1)                &
              , it_entrain_up(ii,jj,1),   it_detrain_up(ii,jj,1)              &
              , it_entrain_dwn(ii,jj,1),  it_detrain_dwn(ii,jj,1)             &
              , it_uw_dp(ii,jj,1),        it_vw_dp(ii,jj,1)                   &
              , it_uw_shall(ii,jj,1),     it_vw_shall(ii,jj,1)                &
              , it_uw_mid(ii,jj,1),       it_vw_mid(ii,jj,1)                  &
              , it_wqt_flux(ii,jj,1),     it_wthetal_flux(ii,jj,1)            &
              , it_wthetav_flux(ii,jj,1), it_wql_flux(ii,jj,1)                &
              , it_mf_deep(ii,jj,1),      it_mf_congest(ii,jj,1)              &
              , it_mf_shall(ii,jj,1),     it_mf_midlev(ii,jj,1)               &
              , it_dt_deep(ii,jj,1),      it_dt_congest(ii,jj,1)              &
              , it_dt_shall(ii,jj,1),     it_dt_midlev(ii,jj,1)               &
              , it_dq_deep(ii,jj,1),      it_dq_congest(ii,jj,1)              &
              , it_dq_shall(ii,jj,1),     it_dq_midlev(ii,jj,1)               &
              , it_du_deep(ii,jj,1),      it_du_congest(ii,jj,1)              &
              , it_du_shall(ii,jj,1),     it_du_midlev(ii,jj,1)               &
              , it_dv_deep(ii,jj,1),      it_dv_congest(ii,jj,1)              &
              , it_dv_shall(ii,jj,1),     it_dv_midlev(ii,jj,1)               &
              , ind_cape_reduced(ii,jj),  cape_ts_used(ii,jj)                 &
              , it_ind_deep(ii,jj),       it_ind_shall(ii,jj)    )

         CASE ( i_convection_vn_5a )

            CALL glue_conv_5a                                                 &
              ( rows*row_length, segments(i)%seg_points                       &
              , n_conv_levels, bl_levels                                      &
              , call_number, i                                                &
              , theta_conv(ii,jj,1), q_conv(ii,jj,1)                          &
              , qcl_conv(ii,jj,1), qcf_conv(ii,jj,1)                          &
              , qrain_conv(ii,jj,1), qgraup_conv(ii,jj,1), qcf2_conv(ii,jj,1) & 
              , cf_liquid_conv(ii,jj,1), cf_frozen_conv(ii,jj,1)              &
              , bulk_cf_conv(ii,jj,1)                                         &
              , p_star(ii,jj), land_sea_mask(ii,jj)                           &
              , u_conv(ii,jj,1), v_conv(ii,jj,1), w(ii,jj,1)                  &
              , tot_tracer(ii,jj,1,1), dthbydt(ii,jj,1)                       &
              , dqbydt(ii,jj,1),   dqclbydt(ii,jj,1), dqcfbydt(ii,jj,1)       &
              , dcflbydt(ii,jj,1), dcffbydt(ii,jj,1), dbcfbydt(ii,jj,1)       &
              , dubydt_p(ii,jj,1), dvbydt_p(ii,jj,1)                          &
              , it_conv_rain(ii,jj),      it_conv_snow(ii,jj)                 &
              , it_conv_rain_3d(ii,jj,1), it_conv_snow_3d(ii,jj,1)            &
              , it_cca0(ii,jj,1),  it_ccb0(ii,jj),   it_cct0(ii,jj)           &
              , it_cclwp0(ii,jj),  it_ccw0(ii,jj,1), it_lcbase0(ii,jj)        &
              , it_cca0_2d(ii,jj), it_lctop(ii,jj),  it_lcca(ii,jj)           &
              , it_cca(ii,jj,1),   it_ccb(ii,jj),    it_cct(ii,jj)            &
              , it_cclwp(ii,jj),   it_ccw(ii,jj,1),  it_lcbase(ii,jj)         &
              , it_cca_2d(ii,jj), freeze_lev(ii,jj), it_dp_cfl_limited(ii,jj) &
              , it_md_cfl_limited(ii,jj)                                      &
              , it_mid_level(ii,jj), it_kterm_deep(ii,jj)                     &
              , it_kterm_shall(ii,jj)                                         &
              , it_precip_dp(ii,jj), it_precip_sh(ii,jj)                      &
              , it_precip_md(ii,jj), it_precip_cg(ii,jj)                      &
              , it_wstar_dn(ii,jj),  it_wstar_up(ii,jj)                       &
              , it_mb1(ii,jj), it_mb2(ii,jj), it_cg_term(ii,jj), n_cumulus(i) &
              , uw0_p(ii,jj), vw0_p(ii,jj), w_max(ii,jj)                      &
              , zlcl(ii,jj), zlcl_uv(ii,jj), ztop(ii,jj), entrain_coef(ii,jj) &
              , deep_flag(ii,jj), past_precip(ii,jj), past_conv_ht(ii,jj)     &
              , it_cape_out(ii,jj), n_deep(i), n_congestus(i), n_shallow(i)   &
              , n_mid(i), r_rho_lev(ii,jj,1), r_theta_lev(ii,jj,0)            &
              , rho_only(ii,jj,1), rho_theta(ii,jj,1), delta_smag(ii,jj)      &
              , exner_layer_boundaries(ii,jj,0), exner_layer_centres(ii,jj,0) &
              , p_layer_boundaries(ii,jj,0), p_layer_centres(ii,jj,0)         &
              , z_theta(ii,jj,1), z_rho(ii,jj,1), timestep_conv               &
              , t1_sd(ii,jj), q1_sd(ii,jj), ntml(ii,jj), ntpar(ii,jj)         &
              , conv_type(ii,jj), l_shallow(ii,jj), l_pc2_diag_sh_pts(ii,jj)  &
              , l_congestus(ii,jj), l_mid(ii,jj), cumulus(ii,jj)              &
              , wstar(ii,jj), wthvs(ii,jj), delthvu(ii,jj), ql_ad(ii,jj)      &
              , qsat_lcl(ii,jj), ftl(ii,jj), fqt(ii,jj), l_tracer, ntra_fld   &
              , ntra_lev, n_cca_lev, l_mixing_ratio, l_mcr_qrain              &
              , l_mcr_qgraup, l_mcr_qcf2, l_calc_dxek                         &
              , l_q_interact, it_up_flux_half(ii,jj,1)                        &
              , it_up_flux(ii,jj,1),      it_dwn_flux(ii,jj,1)                &
              , it_entrain_up(ii,jj,1),   it_detrain_up(ii,jj,1)              &
              , it_entrain_dwn(ii,jj,1),  it_detrain_dwn(ii,jj,1)             &
              , it_uw_dp(ii,jj,1),        it_vw_dp(ii,jj,1)                   &
              , it_uw_shall(ii,jj,1),     it_vw_shall(ii,jj,1)                &
              , it_uw_mid(ii,jj,1),       it_vw_mid(ii,jj,1)                  &
              , it_wqt_flux(ii,jj,1),     it_wthetal_flux(ii,jj,1)            &
              , it_wthetav_flux(ii,jj,1), it_wql_flux(ii,jj,1)                &
              , it_mf_deep(ii,jj,1),      it_mf_congest(ii,jj,1)              &
              , it_mf_shall(ii,jj,1),     it_mf_midlev(ii,jj,1)               &
              , it_dt_deep(ii,jj,1),      it_dt_congest(ii,jj,1)              &
              , it_dt_shall(ii,jj,1),     it_dt_midlev(ii,jj,1)               &
              , it_dq_deep(ii,jj,1),      it_dq_congest(ii,jj,1)              &
              , it_dq_shall(ii,jj,1),     it_dq_midlev(ii,jj,1)               &
              , it_du_deep(ii,jj,1),      it_du_congest(ii,jj,1)              &
              , it_du_shall(ii,jj,1),     it_du_midlev(ii,jj,1)               &
              , it_dv_deep(ii,jj,1),      it_dv_congest(ii,jj,1)              &
              , it_dv_shall(ii,jj,1),     it_dv_midlev(ii,jj,1)               &
              , ind_cape_reduced(ii,jj),  cape_ts_used(ii,jj)                 &
              , it_ind_deep(ii,jj),       it_ind_shall(ii,jj)    )

         CASE ( i_convection_vn_6a )

            CALL glue_conv_6a                                                 &
              ( rows*row_length, segments(i)%seg_points                       &
              , n_conv_levels, bl_levels                                      &
              , call_number, i                                                &
              , theta_conv(ii,jj,1), q_conv(ii,jj,1)                          &
              , qcl_conv(ii,jj,1), qcf_conv(ii,jj,1)                          &
              , qrain_conv(ii,jj,1), qgraup_conv(ii,jj,1), qcf2_conv(ii,jj,1) & 
              , cf_liquid_conv(ii,jj,1), cf_frozen_conv(ii,jj,1)              &
              , bulk_cf_conv(ii,jj,1)                                         &
              , p_star(ii,jj), land_sea_mask(ii,jj)                           &
              , u_conv(ii,jj,1), v_conv(ii,jj,1), w(ii,jj,1)                  &
              , tot_tracer(ii,jj,1,1), dthbydt(ii,jj,1)                       &
              , dqbydt(ii,jj,1),   dqclbydt(ii,jj,1), dqcfbydt(ii,jj,1)       &
              , dcflbydt(ii,jj,1), dcffbydt(ii,jj,1), dbcfbydt(ii,jj,1)       &
              , dubydt_p(ii,jj,1), dvbydt_p(ii,jj,1)                          &
              , it_conv_rain(ii,jj),      it_conv_snow(ii,jj)                 &
              , it_conv_rain_3d(ii,jj,1), it_conv_snow_3d(ii,jj,1)            &
              , it_cca0(ii,jj,1),  it_ccb0(ii,jj),   it_cct0(ii,jj)           &
              , it_cclwp0(ii,jj),  it_ccw0(ii,jj,1), it_lcbase0(ii,jj)        &
              , it_cca0_2d(ii,jj), it_lctop(ii,jj),  it_lcca(ii,jj)           &
              , it_cca(ii,jj,1),   it_ccb(ii,jj),    it_cct(ii,jj)            &
              , it_cclwp(ii,jj),   it_ccw(ii,jj,1),  it_lcbase(ii,jj)         &
              , it_cca_2d(ii,jj), freeze_lev(ii,jj), it_dp_cfl_limited(ii,jj) &
              , it_md_cfl_limited(ii,jj)                                      &
              , it_mid_level(ii,jj), it_kterm_deep(ii,jj)                     &
              , it_kterm_shall(ii,jj)                                         &
              , it_precip_dp(ii,jj), it_precip_sh(ii,jj)                      &
              , it_precip_md(ii,jj), it_precip_cg(ii,jj)                      &
              , it_wstar_dn(ii,jj),  it_wstar_up(ii,jj)                       &
              , it_mb1(ii,jj), it_mb2(ii,jj), it_cg_term(ii,jj), n_cumulus(i) &
              , uw0_p(ii,jj), vw0_p(ii,jj), w_max(ii,jj)                      &
              , zlcl(ii,jj), zlcl_uv(ii,jj), ztop(ii,jj), entrain_coef(ii,jj) &
              , deep_flag(ii,jj), past_precip(ii,jj), past_conv_ht(ii,jj)     &
              , it_cape_out(ii,jj), n_deep(i), n_congestus(i), n_shallow(i)   &
              , n_mid(i), r_rho_lev(ii,jj,1), r_theta_lev(ii,jj,0)            &
              , rho_only(ii,jj,1), rho_theta(ii,jj,1), delta_smag(ii,jj)      &
              , exner_layer_boundaries(ii,jj,0), exner_layer_centres(ii,jj,0) &
              , p_layer_boundaries(ii,jj,0), p_layer_centres(ii,jj,0)         &
              , z_theta(ii,jj,1), z_rho(ii,jj,1), timestep_conv               &
              , t1_sd(ii,jj), q1_sd(ii,jj), ntml(ii,jj), ntpar(ii,jj)         &
              , conv_type(ii,jj), l_shallow(ii,jj), l_pc2_diag_sh_pts(ii,jj)  &
              , l_congestus(ii,jj), l_mid(ii,jj), cumulus(ii,jj)              &
              , wstar(ii,jj), wthvs(ii,jj), delthvu(ii,jj), ql_ad(ii,jj)      &
              , qsat_lcl(ii,jj), ftl(ii,jj), fqt(ii,jj), l_tracer, ntra_fld   &
              , ntra_lev, n_cca_lev, l_mixing_ratio, l_mcr_qrain              &
              , l_mcr_qgraup, l_mcr_qcf2, l_calc_dxek                         &
              , l_q_interact, it_up_flux_half(ii,jj,1)                        &
              , it_up_flux(ii,jj,1),      it_dwn_flux(ii,jj,1)                &
              , it_entrain_up(ii,jj,1),   it_detrain_up(ii,jj,1)              &
              , it_entrain_dwn(ii,jj,1),  it_detrain_dwn(ii,jj,1)             &
              , it_uw_dp(ii,jj,1),        it_vw_dp(ii,jj,1)                   &
              , it_uw_shall(ii,jj,1),     it_vw_shall(ii,jj,1)                &
              , it_uw_mid(ii,jj,1),       it_vw_mid(ii,jj,1)                  &
              , it_wqt_flux(ii,jj,1),     it_wthetal_flux(ii,jj,1)            &
              , it_wthetav_flux(ii,jj,1), it_wql_flux(ii,jj,1)                &
              , it_mf_deep(ii,jj,1),      it_mf_congest(ii,jj,1)              &
              , it_mf_shall(ii,jj,1),     it_mf_midlev(ii,jj,1)               &
              , it_dt_deep(ii,jj,1),      it_dt_congest(ii,jj,1)              &
              , it_dt_shall(ii,jj,1),     it_dt_midlev(ii,jj,1)               &
              , it_dq_deep(ii,jj,1),      it_dq_congest(ii,jj,1)              &
              , it_dq_shall(ii,jj,1),     it_dq_midlev(ii,jj,1)               &
              , it_du_deep(ii,jj,1),      it_du_congest(ii,jj,1)              &
              , it_du_shall(ii,jj,1),     it_du_midlev(ii,jj,1)               &
              , it_dv_deep(ii,jj,1),      it_dv_congest(ii,jj,1)              &
              , it_dv_shall(ii,jj,1),     it_dv_midlev(ii,jj,1)               &
              , ind_cape_reduced(ii,jj),  cape_ts_used(ii,jj)                 &
              , it_ind_deep(ii,jj),       it_ind_shall(ii,jj)    )

         CASE DEFAULT ! i_convection_vn

         errorstatus = 10
         WRITE (cmessage,'(A,A,I6)') 'Convection scheme version value',     &
                   ' not recognised. i_convection_vn = ',i_convection_vn
         CALL Ereport ( RoutineName, errorstatus, cmessage)

      END SELECT ! i_convection_vn

    END DO  ! loop over number of segments

    ! Deallocate segmentation variables
    DEALLOCATE( segments  )
    DEALLOCATE( n_cumulus )
    DEALLOCATE( n_deep )
    DEALLOCATE( n_shallow )
    DEALLOCATE( n_congestus )
    DEALLOCATE( n_mid )

!$OMP END PARALLEL

!===============================================================================
! End of segmented region
!===============================================================================

    ! Update convective cloud diagnostics section 5
    ! and prognostics in section 0  after sub-step

! DEPENDS ON : update_conv_cloud
    CALL update_conv_cloud (rows, row_length, tdims%k_end, qdims%k_end, &
       n_cca_lev, ndecay_steps, n_conv_calls,                           &
       call_number,                                                     &
       it_ccb, it_cct, it_lcbase, it_lctop,it_ccb0, it_cct0, it_lcbase0,&
       one_over_conv_calls, timestep,decay_time,                        &
       p_layer_boundaries, cld_life_3d,                                 &
       it_cca_2d, it_lcca, it_cca, it_ccw, it_cclwp,                    &
       it_cca0_2d, it_cca0,it_ccw0,                                     &
       ! in/out
       ccb, cct, lcbase,lctop, ccb0, cct0, lcbase0,                     &
       ccb0_local, cct0_local,lcbase0_local,                            &
       lcca,  cca_2d,   cclwp,  cca,  ccw,                              &
       cca0_2d,  cclwp0, cca0, ccw0,                                    &
       it_cclwp0, cclwp0_local, cca0_local, ccw0_local)

    ! Update other diagnostics after sub-step

! DEPENDS ON : update_conv_diags
    CALL update_conv_diags(rows, row_length, tdims%k_end, qdims%k_end,    &
      n_conv_levels, call_number, n_conv_calls,                           &
      ntml,ntpar, freeze_lev, it_kterm_deep, it_kterm_shall, it_cg_term,  &
      cumulus, l_shallow, l_congestus, l_congestus2,it_dp_cfl_limited,    &
      it_md_cfl_limited, it_mid_level,                                    &
      one_over_conv_calls, timestep_conv , wstar,                         &
      exner_theta_levels, z_rho,                                          &
      ! Input (values output by latest convection call)
      it_cape_out, it_conv_rain, it_conv_snow,                            &
      it_precip_dp, it_precip_sh, it_precip_md, it_precip_cg,             &
      ind_cape_reduced, cape_ts_used, it_ind_deep, it_ind_shall,          &
      it_wstar_up, it_mb1, it_mb2,                                        &
      dthbydt, dqbydt, dqclbydt, dqcfbydt, dcflbydt, dcffbydt, dbcfbydt,  &
      dubydt_p, dvbydt_p,                                                 &
      it_up_flux, it_up_flux_half, it_dwn_flux, it_entrain_up,            &
      it_entrain_dwn, it_detrain_up, it_detrain_dwn,                      &
      it_conv_rain_3d, it_conv_snow_3d,                                   &
      it_wqt_flux,it_wthetal_flux,it_wthetav_flux,it_wql_flux,            &
      it_uw_dp, it_vw_dp, it_uw_shall, it_vw_shall, it_uw_mid, it_vw_mid, &
      it_mf_deep, it_mf_congest, it_mf_shall, it_mf_midlev,               &
      it_dt_deep, it_dt_congest, it_dt_shall, it_dt_midlev,               &
      it_dq_deep, it_dq_congest, it_dq_shall, it_dq_midlev,               &
      it_du_deep, it_du_congest, it_du_shall, it_du_midlev,               &
      it_dv_deep, it_dv_congest, it_dv_shall, it_dv_midlev,               &
      ! in/out diagnostics
      conv_rain, conv_snow)

    !-----------------------------------------------------------
    ! Mid-level conv is only possible on subsequent convective
    ! substeps if convection has occurred on the previous
    ! convective substep unless using rediagnosis option.
    !-----------------------------------------------------------
    IF (.NOT. l_rediagnosis) THEN
      DO j = 1, rows
        DO i= 1, row_length
          l_mid(i,j) = it_mid_level(i,j) .OR. cumulus(i,j)
        END DO
      END DO
    END IF


! ----------------------------------------------------------------------
! Section 2.8 Add on theta and q increments, qcl and qcf increments.
! ----------------------------------------------------------------------

    ! add on increments to theta and q for next convection call
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k)
    DO k = 1, n_conv_levels
      DO j = 1, rows
        DO i = 1, row_length
          theta_inc(i,j,k) = theta_inc(i,j,k)                     &
                           + dthbydt(i,j,k) * timestep_conv
          q_inc(i,j,k) = q_inc(i,j,k)                             &
                          + dqbydt(i,j,k) * timestep_conv
          qcl_inc(i,j,k) = qcl_inc(i,j,k) +                       &
                         (dqclbydt(i,j,k) * timestep_conv)
          qcf_inc(i,j,k) = qcf_inc(i,j,k) +                       &
                         (dqcfbydt(i,j,k) * timestep_conv)
          cf_liquid_inc(i,j,k) = cf_liquid_inc(i,j,k) +           &
                         (dcflbydt(i,j,k) * timestep_conv)
          cf_frozen_inc(i,j,k) = cf_frozen_inc(i,j,k) +           &
                         (dcffbydt(i,j,k) * timestep_conv)
          bulk_cf_inc(i,j,k)   = bulk_cf_inc(i,j,k) +             &
                         (dbcfbydt(i,j,k) * timestep_conv)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! Final sweep  reset theta_star to theta_conv

    IF (call_number == n_conv_calls) THEN

      IF (l_rediagnosis) THEN
      ! No longer need total increment arrays
        IF (l_mcr_qcf2) THEN
          DEALLOCATE( qcf2_inc_step )
        END IF
        IF (l_mcr_qgraup) THEN
          DEALLOCATE( qgraup_inc_step )
        END IF
        IF (l_mcr_qrain) THEN
          DEALLOCATE( qrain_inc_step )
        END IF
        DEALLOCATE( bulk_cf_inc_step )
        DEALLOCATE( cf_frozen_inc_step )
        DEALLOCATE( cf_liquid_inc_step )
        DEALLOCATE( qcf_inc_step )
        DEALLOCATE( qcl_inc_step )
        DEALLOCATE( q_inc_step )
        DEALLOCATE( theta_inc_step )
      END IF

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k, i, j, tmpB4)
      IF (l_safe_conv) THEN
        ! remove any q added for a safe qmin input profile
!$OMP DO SCHEDULE(STATIC)
        DO k = 1, n_conv_levels
          DO j = 1, rows
            DO i = 1, row_length
              q_star(i,j,k)   = q_conv(i,j,k) - dq_add(i,j,k)              &
                                  + dqbydt(i,j,k) * timestep_conv
!              q_star(i,j,k)   = q_conv(i,j,k)                               &
!                                  + dqbydt(i,j,k) * timestep_conv
            END DO
          END DO
        END DO
!$OMP END DO
      ELSE 
!$OMP DO SCHEDULE(STATIC)
        DO k = 1, n_conv_levels
          DO j = 1, rows
            DO i = 1, row_length
              q_star(i,j,k)   = q_conv(i,j,k)                              &
                                     + dqbydt(i,j,k) * timestep_conv
            END DO
          END DO
        END DO
!$OMP END DO
      END IF

!$OMP DO SCHEDULE(STATIC)
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length

            theta_star(i,j,k) = theta_conv(i,j,k)                        &
                          + dthbydt(i,j,k) * timestep_conv
            qcl_star(i,j,k) = qcl_conv(i,j,k)                            &
                            + (dqclbydt(i,j,k) * timestep_conv)
            qcf_star(i,j,k) = qcf_conv(i,j,k)                            &
                            + (dqcfbydt(i,j,k) * timestep_conv)
            cf_liquid_star(i,j,k) = cf_liquid_conv(i,j,k)                &
                                 + (dcflbydt(i,j,k) * timestep_conv)
            cf_frozen_star(i,j,k) = cf_frozen_conv(i,j,k)                &
                                 + (dcffbydt(i,j,k) * timestep_conv)
            bulk_cf_star(i,j,k)   = bulk_cf_conv(i,j,k)                  &
                                 + (dbcfbydt(i,j,k) * timestep_conv)

            IF (i_pc2_conv_coupling /= pc2_conv_original) THEN
              ! Protect against generation of inconsistently low cloud 
              ! fraction implying very high in-cloud condensate amounts.
              ! In-cloud condensate amounts above 2.0e-3 lead to
              ! cloud fraction being increased (up to a value of 1.0)

              ! Liquid cloud fraction
              IF (cf_liquid_star(i,j,k) > 0.0) THEN
                IF ((qcl_star(i,j,k)/cf_liquid_star(i,j,k))>2.0e-3) THEN
                  tmpB4 = cf_liquid_star(i,j,k) 
                  cf_liquid_star(i,j,k) = MIN(1.0,qcl_star(i,j,k)/2.0e-3) 
                  cf_liquid_inc(i,j,k)  = cf_liquid_inc(i,j,k)          & 
                    + cf_liquid_star(i,j,k) - tmpB4 
                  IF (l_fixbug_pc2_mixph) THEN
                    bulk_cf_star(i,j,k) = bulk_cf_star(i,j,k)           &
                    + cf_liquid_star(i,j,k) - tmpB4 
                    bulk_cf_inc(i,j,k) = bulk_cf_inc(i,j,k)             &
                    + cf_liquid_star(i,j,k) - tmpB4 
                    cf_liquid_incr_diag_conv(i,j,k) =                   &
                    cf_liquid_incr_diag_conv(i,j,k)                     &
                    + cf_liquid_star(i,j,k) - tmpB4 
                    bulk_cf_incr_diag_conv(i,j,k) =                     &
                    bulk_cf_incr_diag_conv(i,j,k)                       &
                    + cf_liquid_star(i,j,k) - tmpB4 
                  ELSE
                    dcflbydt(i,j,k)       = dcflbydt(i,j,k)             & 
                    + (cf_liquid_star(i,j,k) - tmpB4) / timestep_conv 
                  END IF
                END IF
              END IF

              ! Ice cloud fraction
              IF (cf_frozen_star(i,j,k) > 0.0) THEN
                IF ((qcf_star(i,j,k)/cf_frozen_star(i,j,k))>2.0e-3) THEN 
                  tmpB4 = cf_frozen_star(i,j,k) 
                  cf_frozen_star(i,j,k) = MIN(1.0,qcf_star(i,j,k)/2.0e-3) 
                  cf_frozen_inc(i,j,k)  = cf_frozen_inc(i,j,k)          & 
                    + cf_frozen_star(i,j,k) - tmpB4 
                  IF (l_fixbug_pc2_mixph) THEN
                    bulk_cf_star(i,j,k) = bulk_cf_star(i,j,k)           &
                    + cf_frozen_star(i,j,k) - tmpB4 
                    bulk_cf_inc(i,j,k) = bulk_cf_inc(i,j,k)             &
                    + cf_frozen_star(i,j,k) - tmpB4 
                    cf_frozen_incr_diag_conv(i,j,k) =                   &
                    cf_frozen_incr_diag_conv(i,j,k)                     &
                    + cf_frozen_star(i,j,k) - tmpB4 
                    bulk_cf_incr_diag_conv(i,j,k) =                     &
                    bulk_cf_incr_diag_conv(i,j,k)                       &
                    + cf_frozen_star(i,j,k) - tmpB4 
                  ELSE
                    dcffbydt(i,j,k)       = dcffbydt(i,j,k)             & 
                    + (cf_frozen_star(i,j,k) - tmpB4) / timestep_conv 
                  END IF
                END IF 
              END IF

            END IF ! i_pc2_conv_coupling

          END DO
        END DO
      END DO
!$OMP END DO

!$OMP END PARALLEL 

    ELSE

      IF (l_rediagnosis) THEN

      ! not final sub-step so add on part of increment from
      ! slow physics and dynamics

        DO k = 1, n_conv_levels
          DO j = 1, rows
            DO i = 1, row_length

              theta_conv(i,j,k) = theta_conv(i,j,k)                     &
                           + fraction_step*theta_inc_step(i,j,k)        &
                           + dthbydt(i,j,k) * timestep_conv
              q_conv(i,j,k) = q_conv(i,j,k)                             &
                            + fraction_step*q_inc_step(i,j,k)           &
                            + dqbydt(i,j,k) * timestep_conv
              qcl_conv(i,j,k) = qcl_conv(i,j,k)                         &
                           + fraction_step*qcl_inc_step(i,j,k)          &
                           +(dqclbydt(i,j,k) * timestep_conv)
              qcf_conv(i,j,k) = qcf_conv(i,j,k)                         &
                           + fraction_step*qcf_inc_step(i,j,k)          &
                           +(dqcfbydt(i,j,k) * timestep_conv)
              cf_liquid_conv(i,j,k) = cf_liquid_conv(i,j,k)             &
                           + fraction_step*cf_liquid_inc_step(i,j,k)    &
                           +(dcflbydt(i,j,k) * timestep_conv)
              cf_frozen_conv(i,j,k) = cf_frozen_conv(i,j,k)             &
                           + fraction_step*cf_frozen_inc_step(i,j,k)    &
                           + (dcffbydt(i,j,k) * timestep_conv)
              bulk_cf_conv(i,j,k)   = bulk_cf_conv(i,j,k)               &
                           + fraction_step*bulk_cf_inc_step(i,j,k)      &
                           +(dbcfbydt(i,j,k) * timestep_conv)

              IF (i_pc2_conv_coupling /= pc2_conv_original) THEN
                ! Protect against generation of inconsistently low cloud 
                ! fraction implying high in-cloud condensate amounts.
                ! In-cloud condensate amounts above 2.0e-3 lead to
                ! cloud fraction being increased (up to a value of 1.0)

                ! Liquid cloud fraction
                IF (cf_liquid_conv(i,j,k) > 0.0) THEN
                  IF ( ( qcl_conv(i,j,k) / cf_liquid_conv(i,j,k) )      &
                        > 2.0e-3) THEN 
                    tmpB4 = cf_liquid_conv(i,j,k) 
                    cf_liquid_conv(i,j,k) = MIN(1.0,                    &
                                                qcl_conv(i,j,k)/2.0e-3) 
                    cf_liquid_inc(i,j,k)  = cf_liquid_inc(i,j,k)        & 
                      + cf_liquid_conv(i,j,k) - tmpB4 
                    IF (l_fixbug_pc2_mixph) THEN
                      bulk_cf_conv(i,j,k) = bulk_cf_conv(i,j,k)         &
                     + cf_liquid_conv(i,j,k) - tmpB4 
                      bulk_cf_inc(i,j,k) = bulk_cf_inc(i,j,k)           &
                     + cf_liquid_conv(i,j,k) - tmpB4 
                      cf_liquid_incr_diag_conv(i,j,k) =                 &
                      cf_liquid_incr_diag_conv(i,j,k)                   &
                     + cf_liquid_conv(i,j,k) - tmpB4 
                      bulk_cf_incr_diag_conv(i,j,k) =                   &
                      bulk_cf_incr_diag_conv(i,j,k)                     &
                     + cf_liquid_conv(i,j,k) - tmpB4 
                    ELSE
                      dcflbydt(i,j,k)       = dcflbydt(i,j,k)           & 
                      + (cf_liquid_conv(i,j,k) - tmpB4) / timestep_conv 
                    END IF
                  END IF 
                END IF

                ! Ice cloud fraction
                IF (cf_frozen_conv(i,j,k) > 0.0) THEN
                  IF ( ( qcf_conv(i,j,k) / cf_frozen_conv(i,j,k) )      &
                        > 2.0e-3) THEN 
                    tmpB4 = cf_frozen_conv(i,j,k) 
                    cf_frozen_conv(i,j,k) = MIN(1.0,                    &
                                                qcf_conv(i,j,k)/2.0e-3) 
                    cf_frozen_inc(i,j,k) = cf_frozen_inc(i,j,k)         & 
                      + cf_frozen_conv(i,j,k) - tmpB4 
                    IF (l_fixbug_pc2_mixph) THEN
                      bulk_cf_conv(i,j,k) = bulk_cf_conv(i,j,k)         &
                      + cf_frozen_conv(i,j,k) - tmpB4 
                      bulk_cf_inc(i,j,k) = bulk_cf_inc(i,j,k)           &
                      + cf_frozen_conv(i,j,k) - tmpB4 
                      cf_frozen_incr_diag_conv(i,j,k) =                 &
                      cf_frozen_incr_diag_conv(i,j,k)                   &
                      + cf_frozen_conv(i,j,k) - tmpB4 
                      bulk_cf_incr_diag_conv(i,j,k) =                   &
                      bulk_cf_incr_diag_conv(i,j,k)                     &
                      + cf_frozen_conv(i,j,k) - tmpB4 
                    ELSE
                      dcffbydt(i,j,k) = dcffbydt(i,j,k)                 & 
                      + (cf_frozen_conv(i,j,k) - tmpB4) / timestep_conv 
                    END IF
                  END IF 
                END IF

              END IF ! i_pc2_conv_coupling

            END DO
          END DO
        END DO

        ! Extra prognostics if present
        IF (l_mcr_qrain) THEN
          DO k=1,qdims%k_end
            DO j=1,rows
              DO i=1,row_length
                qrain_conv(i,j,k)   = qrain_conv(i,j,k)                  &
                                     + fraction_step*qrain_inc_step(i,j,k)
              END DO ! i
            END DO ! j
          END DO ! k
        END IF
        IF (l_mcr_qgraup) THEN
          DO k=1,qdims%k_end
            DO j=1,rows
              DO i=1,row_length
                qgraup_conv(i,j,k)   = qgraup_conv(i,j,k)                &
                                      + fraction_step*qgraup_inc_step(i,j,k)
              END DO ! i
            END DO ! j
          END DO ! k
        END IF
        IF (l_mcr_qcf2) THEN
          DO k=1,qdims%k_end
            DO j=1,rows
              DO i=1,row_length
                qcf2_conv(i,j,k)   = qcf2_conv(i,j,k)                    &
                                    + fraction_step*qcf2_inc_step(i,j,k)
              END DO ! i
            END DO ! j
          END DO ! k
        END IF

        IF (l_mom) THEN
          DO k = 1, n_conv_levels
            DO j = 1, rows
              DO i = 1, row_length
                u_conv(i,j,k) = u_conv(i,j,k)                           &
                              + fraction_step*u_inc_step(i,j,k)         &
                              + dubydt_p(i,j,k) * timestep_conv
                v_conv(i,j,k) = v_conv(i,j,k)                           &
                              + fraction_step*v_inc_step(i,j,k)         &
                              + dvbydt_p(i,j,k) * timestep_conv

              END DO
            END DO
          END DO
        END IF   ! l_mom

      ELSE    ! original code
        ! Note no need to update extra prognostics if present 

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, tmpB4) 

!$OMP DO SCHEDULE(STATIC)
        DO k = 1, n_conv_levels
          DO j = 1, rows
            DO i = 1, row_length
              theta_conv(i,j,k) = theta_conv(i,j,k)                     &
                           + dthbydt(i,j,k) * timestep_conv
              q_conv(i,j,k) = q_conv(i,j,k)                             &
                           + dqbydt(i,j,k) * timestep_conv
              qcl_conv(i,j,k) = qcl_conv(i,j,k)                         &
                           +(dqclbydt(i,j,k) * timestep_conv)
              qcf_conv(i,j,k) = qcf_conv(i,j,k)                         &
                           +(dqcfbydt(i,j,k) * timestep_conv)
              cf_liquid_conv(i,j,k) = cf_liquid_conv(i,j,k)             &
                           +(dcflbydt(i,j,k) * timestep_conv)
              cf_frozen_conv(i,j,k) = cf_frozen_conv(i,j,k)             &
                           + (dcffbydt(i,j,k) * timestep_conv)
              bulk_cf_conv(i,j,k)   = bulk_cf_conv(i,j,k)               &
                           +(dbcfbydt(i,j,k) * timestep_conv)

              IF (i_pc2_conv_coupling /= pc2_conv_original) THEN
                ! Protect against generation of inconsistently low cloud 
                ! fraction implying high in-cloud condensate amounts.
                ! In-cloud condensate amounts above 2.0e-3 lead to
                ! cloud fraction being increased (up to a value of 1.0)

                ! Liquid cloud fraction
                IF ( cf_liquid_conv(i,j,k) > 0.0) THEN
                  IF ( ( qcl_conv(i,j,k) / cf_liquid_conv(i,j,k) )      &
                        > 2.0e-3) THEN 
                    tmpB4 = cf_liquid_conv(i,j,k) 
                    cf_liquid_conv(i,j,k) = MIN(1.0,                    &
                                                qcl_conv(i,j,k)/2.0e-3) 
                    cf_liquid_inc(i,j,k)  = cf_liquid_inc(i,j,k)        & 
                      + cf_liquid_conv(i,j,k) - tmpB4
                    IF (l_fixbug_pc2_mixph) THEN
                      bulk_cf_conv(i,j,k) = bulk_cf_conv(i,j,k)         &
                      + cf_liquid_conv(i,j,k) - tmpB4 
                      bulk_cf_inc(i,j,k) = bulk_cf_inc(i,j,k)           &
                      + cf_liquid_conv(i,j,k) - tmpB4 
                      cf_liquid_incr_diag_conv(i,j,k) =                 &
                      cf_liquid_incr_diag_conv(i,j,k)                   &
                      + cf_liquid_conv(i,j,k) - tmpB4 
                      bulk_cf_incr_diag_conv(i,j,k) =                   &
                      bulk_cf_incr_diag_conv(i,j,k)                     &
                      + cf_liquid_conv(i,j,k) - tmpB4 
                    ELSE
                      dcflbydt(i,j,k)       = dcflbydt(i,j,k)           & 
                      + (cf_liquid_conv(i,j,k) - tmpB4) / timestep_conv 
                    END IF
                  END IF 
                END IF

                ! Ice cloud fraction
                IF ( cf_frozen_conv(i,j,k) > 0.0) THEN
                  IF ( ( qcf_conv(i,j,k) / cf_frozen_conv(i,j,k) )       &
                       > 2.0e-3) THEN 
                    tmpB4 = cf_frozen_conv(i,j,k) 
                    cf_frozen_conv(i,j,k) = MIN(1.0,                    &
                                                qcf_conv(i,j,k)/2.0e-3) 
                    cf_frozen_inc(i,j,k)  = cf_frozen_inc(i,j,k)        &
                      + cf_frozen_conv(i,j,k) - tmpB4 
                    IF (l_fixbug_pc2_mixph) THEN
                      bulk_cf_conv(i,j,k) = bulk_cf_conv(i,j,k)         &
                      + cf_frozen_conv(i,j,k) - tmpB4 
                      bulk_cf_inc(i,j,k) = bulk_cf_inc(i,j,k)           &
                      + cf_frozen_conv(i,j,k) - tmpB4 
                      cf_frozen_incr_diag_conv(i,j,k) =                 &
                      cf_frozen_incr_diag_conv(i,j,k)                   &
                      + cf_frozen_conv(i,j,k) - tmpB4 
                      bulk_cf_incr_diag_conv(i,j,k) =                   &
                      bulk_cf_incr_diag_conv(i,j,k)                     &
                      + cf_frozen_conv(i,j,k) - tmpB4 
                    ELSE
                      dcffbydt(i,j,k)       = dcffbydt(i,j,k)           &
                      + (cf_frozen_conv(i,j,k) - tmpB4) / timestep_conv 
                    END IF
                  END IF
                END IF 

              END IF ! i_pc2_conv_coupling

            END DO
          END DO
        END DO
!$OMP END DO

        IF (l_mom) THEN

!$OMP DO SCHEDULE(STATIC)
          DO k = 1, n_conv_levels
            DO j = 1, rows
              DO i = 1, row_length
                u_conv(i,j,k) = u_conv(i,j,k)                           &
                             + dubydt_p(i,j,k) * timestep_conv
                v_conv(i,j,k) = v_conv(i,j,k)                           &
                             + dvbydt_p(i,j,k) * timestep_conv

              END DO
            END DO
          END DO
!$OMP END DO

        END IF    !l_mom

!$OMP END PARALLEL

      END IF    ! l_rediagnosis

    END IF    ! step number

    IF (isrfexcnvgust == ip_srfexwithcnv) THEN
!           Store convective mass fluxes at cloud base
!           if required for surface exchange.
      DO j = 1, rows
        DO i = 1, row_length
          IF (ccb(i,j) > 0) THEN
            ddmfx(i,j)=dwn_flux(i,j,ccb(i,j))
          ELSE
            ddmfx(i,j)=0.0
          END IF
        END DO
      END DO
    END IF

    ! diagnose number of convecting points
    IF (printstatus == prstatus_diag) THEN
      IF (n_conv_calls  >   1) THEN
        n_conv_points = 0
        DO j = 1, rows
          DO i = 1, row_length
            this_point = 0
            DO k = 1, n_conv_levels
              IF (ABS(dthbydt(i,j,k) * timestep)  >   0.0001) THEN
                this_point = 1
              END IF
            END DO
            n_conv_points = n_conv_points + this_point
          END DO
        END DO
        IF (n_proc  >   1) THEN
          CALL gc_isum(1, n_proc, info, n_conv_points )
        END IF
        IF (me  ==  0) THEN
          ! Note max points assumed to be in range i12 consistent with i & j
          ! dimensions being printed as i6
          WRITE(6,'(a11,i3,a5,i12,a19)') ' conv CALL ',call_number,' has ', &
                     n_conv_points,' convecting points '
        END IF

      END IF
    END IF

!--------------------------------------------------------------------------
! Wind increments - interpolation back to uv grid now in atmos_physics2
!                   so that swap-bound calls can be grouped with others
!                   to save CPU.
!--------------------------------------------------------------------------

  END DO ! loop over number of convection calls

!============================================================================
!  End of Convection sub-stepping loop
!============================================================================

!============================================================================
! Section 3. After convection substepping
!============================================================================

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)

! ----------------------------------------------------------------------
! 3.1 Check that CCA doesn't exceed 1.0
! ----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO k=1, n_cca_lev
    DO j=1, rows
      DO i=1, row_length
        cca(i,j,k) = MIN(cca(i,j,k), 1.0)
      END DO
    END DO
  END DO
!$OMP END DO


! ----------------------------------------------------------------------
! Section 3.2 Copy increment fields for PC2 inhomog diagnostics.
! ----------------------------------------------------------------------

  IF(l_qcl_incr_cinh) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, qdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          qcl_incr_inhom_diag(i,j,k) =                            &
                      qcl_inc(i,j,k) - qcl_incr_inhom_diag(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO

  END IF                   ! on STASHflag

  IF(l_qcf_incr_cinh) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, qdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          qcf_incr_inhom_diag(i,j,k) =                            &
                      qcf_inc(i,j,k) - qcf_incr_inhom_diag(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO
  END IF                   ! on STASHflag

  IF(l_bcf_incr_cinh) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, qdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          bulk_cf_incr_inhom_diag(i,j,k) =                        &
              bulk_cf_inc(i,j,k) - bulk_cf_incr_inhom_diag(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO
  END IF                   ! on STASHflag

  IF(l_cfl_incr_cinh .OR. l_calc_dxek) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, qdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          cf_liquid_incr_inhom_diag(i,j,k) =                      &
          cf_liquid_inc(i,j,k) - cf_liquid_incr_inhom_diag(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO

  END IF                   ! on STASHflag

  IF(l_cff_incr_cinh) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, qdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          cf_frozen_incr_inhom_diag(i,j,k) =                      &
          cf_frozen_inc(i,j,k) - cf_frozen_incr_inhom_diag(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO
  END IF                   ! on STASHflag

!==========================================================================
! Section 4 Treat tracers
!==========================================================================

!----------------------------------------------------------------------
! Section 4.1 Copy tracers back into variables from tot_tracers
!----------------------------------------------------------------------
  IF ( l_tracer ) THEN

!$OMP SINGLE
    ntra_tmp = 0
!$OMP END SINGLE

    IF (l_dust) THEN
!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        dust_div1(1:row_length, 1:rows, k) =            &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        dust_div2(1:row_length, 1:rows, k) =            &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO
      IF(.NOT.l_twobin_dust) THEN

!$OMP SINGLE
        ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
        DO k=1, model_levels
          dust_div3(1:row_length, 1:rows, k) =            &
               tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
        END DO
!$OMP END DO

!$OMP SINGLE
        ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
        DO k=1, model_levels
          dust_div4(1:row_length, 1:rows, k) =            &
               tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
        END DO
!$OMP END DO

!$OMP SINGLE
        ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
        DO k=1, model_levels
          dust_div5(1:row_length, 1:rows, k) =            &
               tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
        END DO
!$OMP END DO

!$OMP SINGLE
        ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
        DO k=1, model_levels
          dust_div6(1:row_length, 1:rows, k) =            &
               tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
        END DO
!$OMP END DO
      END IF
    END IF

    IF (l_sulpc_so2) THEN
!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        so2       (1:row_length, 1:rows, k) =           &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        so4_aitken(1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels  
        so4_accu  (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        so4_diss  (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO

      IF (l_sulpc_dms) THEN

!$OMP SINGLE
        ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        dms       (1:row_length, 1:rows, k) =         &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO
      END IF

      IF (l_sulpc_nh3) THEN

!$OMP SINGLE
        ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        nh3       (1:row_length, 1:rows, k) =         &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO
      END IF
    END IF

    IF (l_soot) THEN

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        soot_new  (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp) 
      END DO
!$OMP END DO

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        soot_agd  (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        soot_cld  (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO

    END IF

    IF (l_biomass) THEN
!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        bmass_new  (1:row_length, 1:rows, k) =          &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        bmass_agd  (1:row_length, 1:rows, k) =          &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        bmass_cld  (1:row_length, 1:rows, k) =          &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO
    END IF

    IF (l_co2_interactive) THEN

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        co2       (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO
    END IF

    IF (l_ocff) THEN

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        ocff_new  (1:row_length, 1:rows, k) =           &
           tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        ocff_agd  (1:row_length, 1:rows, k) =           &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        ocff_cld  (1:row_length, 1:rows, k) =           &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO
    END IF

    IF (l_nitrate) THEN

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        nitr_acc  (1:row_length, 1:rows, k) =           &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO

!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        nitr_diss  (1:row_length, 1:rows, k) =           &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO

    END IF

    IF (l_use_cariolle) THEN
!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        ozone_tracer (1:row_length, 1:rows, k) =        &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO
    END IF

    IF (l_murk_conv) THEN
!$OMP SINGLE
      ntra_tmp = ntra_tmp + 1
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
      DO k=1, model_levels
        aerosol   (1:row_length, 1:rows, k) =           &
             tot_tracer(1:row_length, 1:rows, k, ntra_tmp)
      END DO
!$OMP END DO
    END IF

    IF (tr_vars > 0 ) THEN

      i = 0

!$OMP DO SCHEDULE(STATIC)
      DO k= ntra_tmp+1 , ntra_tmp + tr_vars
        i = k - ntra_tmp
        free_tracers(1:row_length, 1:rows, 1:model_levels, i)=          &
             tot_tracer  (1:row_length, 1:rows, 1:model_levels,         &
             k)
      END DO
!$OMP END DO

!$OMP SINGLE
      ntra_tmp = ntra_tmp + tr_vars
!$OMP END SINGLE
    END IF

    IF (tr_ukca > 0 ) THEN
      
      i = 0

!$OMP DO SCHEDULE(STATIC)
      DO k= ntra_tmp+1, ntra_tmp + tr_ukca
        i = k - ntra_tmp
        ukca_tracers(1:row_length, 1:rows, 1:model_levels, i)=          &
             tot_tracer  (1:row_length, 1:rows, 1:model_levels,         &
             k)
      END DO
!$OMP END DO

!$OMP SINGLE
      ntra_tmp = ntra_tmp + tr_ukca
!$OMP END SINGLE

    END IF

  END IF

!$OMP END PARALLEL
  ! Free up the memory again
  DEALLOCATE( tot_tracer )

!----------------------------------------------------------------------
! Section 4.2 Scavenging for tracers
!----------------------------------------------------------------------

  ! work with tracers only in final cycle
  IF ( cycleno == numcycles ) THEN
    ! Scavenging for aerosol murk (Not part of section 17: aerosol scheme)
    IF (l_murk_source) THEN
! DEPENDS ON: con_scav
      CALL con_scav( timestep                                       & 
  ,    ccb, cct, conv_rain, conv_snow                               &
  ,    aerosol(tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,&
                           1:tdims%k_end) )
    END IF

    ! Section 17 : Aerosol scheme tracers
    ! Scavenging of other tracers should go here....
    ! Scavenge Mineral Dust tracers

    IF (l_dust) THEN

      CALL scnscv2( timestep,                                     &
           rho, q_star, qcl_star, qcf_star, dust_div1,            &
           ccb, cct, conv_rain, conv_snow,                        &
           .TRUE., krain_dust(1), ksnow_dust(1),                  &
           conscav_dust(1,1,1) )

      CALL scnscv2( timestep,                                     &
           rho, q_star, qcl_star, qcf_star, dust_div2,            &
           ccb, cct, conv_rain, conv_snow,                        &
           .TRUE., krain_dust(2), ksnow_dust(2),                  &
           conscav_dust(1,1,2) )
      
      ! Only do these when using six-bin dust     
      IF (.NOT.l_twobin_dust) THEN
        CALL scnscv2( timestep,                                   &
             rho, q_star, qcl_star, qcf_star, dust_div3,          &
             ccb, cct, conv_rain, conv_snow,                      &
             .TRUE., krain_dust(3), ksnow_dust(3),                &
             conscav_dust(1,1,3) )

        CALL scnscv2( timestep,                                   &
             rho, q_star, qcl_star, qcf_star, dust_div4,          &
             ccb, cct, conv_rain, conv_snow,                      &
             .TRUE., krain_dust(4), ksnow_dust(4),                &
             conscav_dust(1,1,4) )

        CALL scnscv2( timestep,                                   &
             rho, q_star, qcl_star, qcf_star, dust_div5,          &
             ccb, cct, conv_rain, conv_snow,                      &
             .TRUE., krain_dust(5), ksnow_dust(5),                &
             conscav_dust(1,1,5) )

        CALL scnscv2( timestep,                                   &
             rho, q_star, qcl_star, qcf_star, dust_div6,          &
             ccb, cct, conv_rain, conv_snow,                      &
             .TRUE., krain_dust(6), ksnow_dust(6),                &
             conscav_dust(1,1,6) )
      END IF
    END IF !L_DUST

    ! Scavenge Sulphur Cycle tracers

    IF (l_sulpc_so2) THEN

      ! Scavenge SO2
      CALL scnwsh2( timestep,                                     &
           rho, q_star, qcl_star, qcf_star, so2,                  &
           ccb, cct,                                              &
           conv_rain, conwash_so2 )

      ! Scavenge NH3 if present

      IF (l_sulpc_nh3) THEN

      ! Scavenge NH3
        CALL ncnwsh2( timestep,                                   &
           rho, q_star, qcl_star, qcf_star, nh3,                  &
           ccb, cct,                                              &
           conv_rain, conwash_nh3 )

      END IF              !End L_sulpc_nh3 condition

      ! Scavenge SO4_AIT
      CALL scnscv2( timestep,                                     &
           rho, q_star, qcl_star, qcf_star, so4_aitken,           &
           ccb, cct, conv_rain, conv_snow,                        &
           .FALSE., krain_so4ait, ksnow_so4ait,                   &
           conscav_so4ait )

      ! Scavenge SO4_ACC
      CALL scnscv2( timestep,                                     &
           rho, q_star, qcl_star, qcf_star, so4_accu,             &
           ccb, cct, conv_rain, conv_snow,                        &
           .FALSE., krain_so4acc, ksnow_so4acc,                   &
           conscav_so4acc )

      ! Scavenge SO4_DIS
      CALL scnscv2( timestep,                                     &
           rho, q_star, qcl_star, qcf_star, so4_diss,             &
           ccb, cct, conv_rain, conv_snow,                        &
           .FALSE., krain_so4dis, ksnow_so4dis,                   &
           conscav_so4dis )

    END IF              !End L_sulpc_so2 condition


    IF (l_soot) THEN  ! If soot modelling is included

    ! Scavenge aged soot
      CALL scnscv2( timestep,                                     &
           rho, q_star, qcl_star, qcf_star, soot_agd,             &
           ccb, cct, conv_rain, conv_snow,                        &
           .FALSE., krain_agedsoot, ksnow_agedsoot,               &
           conscav_agedsoot )

    END IF  ! l_soot

    IF (l_biomass) THEN  ! If biomass aerosol is included

! Scavenge aged biomass aerosol
      CALL scnscv2( timestep,                                     &
           rho, q_star, qcl_star, qcf_star, bmass_agd,            &
           ccb, cct, conv_rain, conv_snow,                        &
           .FALSE., krain_agedbmass, ksnow_agedbmass,             &
           conscav_agedbmass )

    END IF  ! l_biomass

    IF (l_ocff) THEN ! If fossil-fuel OC aerosol is included

! Scavenge aged fossil-fuel organic carbon aerosol
      CALL scnscv2( timestep,                                     &
           rho, q_star, qcl_star, qcf_star, ocff_agd,             &
           ccb, cct, conv_rain, conv_snow,                        &
           .FALSE., krain_agedocff, ksnow_agedocff,               &
           conscav_agedocff )

    END IF  ! l_ocff

    IF (l_nitrate) THEN ! If ammonium nitrate is included

! Scavenge accumulation-mode ammonium nitrate
      CALL scnscv2( timestep,                                     &
           rho, q_star, qcl_star, qcf_star, nitr_acc,             &
           ccb, cct, conv_rain, conv_snow,                        &
           .FALSE., krain_nitracc, ksnow_nitracc,                 &
           conscav_nitracc )
!
! Scavenge dissolved ammonium nitrate
      CALL scnscv2( timestep,                                     &
           rho, q_star, qcl_star, qcf_star, nitr_diss,            &
           ccb, cct, conv_rain, conv_snow,                        &
           .FALSE., krain_nitrdiss, ksnow_nitrdiss,               &
           conscav_nitrdiss )
!
    END IF  ! l_nitrate


  END IF ! CycleNo == NumCycles

! ----------------------------------------------------------------------
! Section 5 Calculate PC2 scheme increments to the increments.
! ----------------------------------------------------------------------
  ! ----------------------------------------------------------------------
  ! Copy increments to T and q before updating with PC2 inhomog
  ! ----------------------------------------------------------------------
  IF (l_T_conv_only) THEN 
    DO k=1, n_conv_levels
      DO j=1, rows
        DO i=1, row_length
          T_incr_conv_only(i,j,k) = T_incr_diag_conv(i,j,k)
        END DO ! i
      END DO ! j
    END DO  ! k
  END IF 
  IF (l_q_conv_only) THEN 
    DO k=1, n_conv_levels
      DO j=1, rows
        DO i=1, row_length
          q_incr_conv_only(i,j,k) = q_incr_diag_conv(i,j,k)
        END DO ! i
      END DO ! j
    END DO  ! k
  END IF 

  ! L_calc_dxek_if2:
  IF (l_calc_dxek) THEN

    ALLOCATE ( t_inc_latest(row_length,rows,1) )

    ALLOCATE ( theta_inc_pc2(row_length,rows,1) )
    ALLOCATE ( q_inc_pc2(row_length,rows,1) )
    ALLOCATE ( qcl_inc_pc2(row_length,rows,1) )
    ALLOCATE ( cfl_inc_pc2(row_length,rows,1) )
    ALLOCATE ( bcf_inc_pc2(row_length,rows,1) )

    ALLOCATE ( t_earliest(row_length,rows,1) )
    ALLOCATE ( q_earliest(row_length,rows,1) )
    ALLOCATE ( qcl_earliest(row_length,rows,1) )
    ALLOCATE ( cfl_earliest(row_length,rows,1) )
    ALLOCATE ( cff_earliest(row_length,rows,1) )
    ALLOCATE ( bcf_earliest(row_length,rows,1) )

    IF (i_pc2_erosion_method == pc2eros_hybrid_allfaces) THEN
      ALLOCATE ( bcf_above(row_length,rows,1) )
      ALLOCATE ( bcf_below(row_length,rows,1) )
    ELSE
      ALLOCATE ( bcf_above(1,1,1) )
      ALLOCATE ( bcf_below(1,1,1) )
    END IF

! Convert potential temperature increment to temperature increment

    DO k=1, n_conv_levels

      DO j=1, rows
        DO i=1, row_length
          t_inc_latest(i,j,1) = theta_inc(i,j,k) * exner_theta_levels(i,j,k)

          t_earliest(i,j,1)   = theta_star(i,j,k)*exner_theta_levels(i,j,k)  &
                                                    - t_inc_latest(i,j,1)
          q_earliest(i,j,1)   = q_star(i,j,k)       - q_inc(i,j,k)
          qcl_earliest(i,j,1) = qcl_star(i,j,k)     - qcl_inc(i,j,k)
          bcf_earliest(i,j,1) = bulk_cf_star(i,j,k) - bulk_cf_inc(i,j,k)
          cfl_earliest(i,j,1) = cf_liquid_star(i,j,k) - cf_liquid_inc(i,j,k)
          cff_earliest(i,j,1) = cf_frozen_star(i,j,k) - cf_frozen_inc(i,j,k)
        END DO ! i
      END DO ! j

      IF (i_pc2_erosion_method == pc2eros_hybrid_allfaces) THEN

        IF (k==1) THEN
          km1=k
        ELSE
          km1=k-1
        END IF

        IF (k==n_conv_levels) THEN
          kp1=k
        ELSE
          kp1=k+1
        END IF

        DO j=1, rows
          DO i=1, row_length
            bcf_above(i,j,1)    = bulk_cf_star(i,j,kp1) - bulk_cf_inc(i,j,kp1)      
            bcf_below(i,j,1)    = bulk_cf_star(i,j,km1) - bulk_cf_inc(i,j,km1) 
          END DO ! i
        END DO ! j

      END IF

      IF (l_pc2_diag_sh .OR. l_micro_eros) THEN

     ! call pc2_hom_conv without erosion term i.e. dbsdtbs_turb_0=0.0
     ! and dbsdtbs_turb_1=0.0

! DEPENDS ON: pc2_hom_conv
        CALL PC2_HOM_CONV(                                          &
  ! Input variables
         p_layer_centres(1,1,k), 1                                  &
   ,     timestep                                                   &
  ! INput variables
   ,     t_earliest(1,1,1), q_earliest(1,1,1), qcl_earliest(1,1,1)  &
   ,     bcf_earliest(1,1,1),cfl_earliest(1,1,1),cff_earliest(1,1,1)&
   ,     t_inc_latest,q_inc(1,1,k),full_zero(1,1,1),full_zero(1,1,1)&
   ,     cf_liquid_incr_inhom_diag(1,1,k)                           &
   ,     bcf_above(1,1,1),bcf_below(1,1,1)                          &
  ! OUTput variables
   ,     theta_inc_PC2,q_inc_PC2,qcl_inc_PC2,bcf_inc_PC2,cfl_inc_PC2&
  ! INput variables (other quantities)
   ,     0.0, 0.0                                                   &
  ! Model switches
   ,     l_mixing_ratio)
      ELSE
     ! call pc2_hom_conv with erosion term   

! DEPENDS ON: pc2_hom_conv
        CALL PC2_HOM_CONV(                                          &
   ! Input variables
         p_layer_centres(1,1,k), 1                                  &
   ,     timestep                                                   &
   ! INput variables
   ,     t_earliest(1,1,1), q_earliest(1,1,1), qcl_earliest(1,1,1)  &
   ,     bcf_earliest(1,1,1),cfl_earliest(1,1,1),cff_earliest(1,1,1)&
   ,     t_inc_latest,q_inc(1,1,k),full_zero(1,1,1),full_zero(1,1,1)&
   ,     cf_liquid_incr_inhom_diag(1,1,k)                           &
   ,     bcf_above(1,1,1),bcf_below(1,1,1)                          &
   ! OUTput variables
   ,     theta_inc_PC2,q_inc_PC2,qcl_inc_PC2,bcf_inc_PC2,cfl_inc_PC2&
   ! INput variables (other quantities)
   ,     dbsdtbs_turb_0, dbsdtbs_turb_1                             &
   ! Model switches
   ,     l_mixing_ratio)
          END IF 
!
! Calculate potential temperature increment (on convect levels only)
! from temperature increment output by PC2_Homog.

      DO j=1,rows
        DO i=1,row_length
          theta_inc_pc2(i,j,1) = theta_inc_pc2(i,j,1)/exner_theta_levels(i,j,k)
        END DO ! i
      END DO ! j

      ! L_q_interact_if1:
      IF (l_q_interact)  THEN

        DO j = 1, rows
          DO i = 1, row_length

        ! Update increments to theta, moisture and cloud fields with additional
        ! increments from the response to environment changes (homogenous),

            theta_inc(i,j,k) = theta_inc(i,j,k) + theta_inc_pc2(i,j,1)
            q_inc(i,j,k)   =   q_inc(i,j,k) +   q_inc_pc2(i,j,1)
            qcl_inc(i,j,k) = qcl_inc(i,j,k) + qcl_inc_pc2(i,j,1)
        ! Not updated   qcf_inc(i,j,k) = qcf_inc(i,j,k)
            cf_liquid_inc(i,j,k) = cf_liquid_inc(i,j,k) + cfl_inc_pc2(i,j,1)
        ! Not updated   cf_frozen_inc(i,j,k) = cf_frozen_inc(i,j,k)
            bulk_cf_inc(i,j,k)   = bulk_cf_inc(i,j,k) +  bcf_inc_pc2(i,j,1) 

        ! ... and update working version of theta, moisture and cloud fields.

            theta_star(i,j,k) = theta_star(i,j,k) + theta_inc_pc2(i,j,1)
            q_star(i,j,k)     =     q_star(i,j,k) +     q_inc_pc2(i,j,1)
            qcl_star(i,j,k)   =   qcl_star(i,j,k) +   qcl_inc_pc2(i,j,1)
        ! Not updated   qcf_star(i,j,k)   =   qcf_star(i,j,k)
            cf_liquid_star(i,j,k) = cf_liquid_star(i,j,k) + cfl_inc_pc2(i,j,1)
        ! Not updated   cf_frozen_n(i,j,k) = cf_frozen_n(i,j,k)
            bulk_cf_star(i,j,k)   = bulk_cf_star(i,j,k) +           &
                                  bcf_inc_pc2(i,j,1)
        ! ... and the diagnostics.

            t_incr_diag_conv(i,j,k) = t_incr_diag_conv(i,j,k)     &
                + exner_theta_levels(i,j,k) * theta_inc_pc2(i,j,1)
            q_incr_diag_conv(i,j,k) = q_incr_diag_conv(i,j,k)     &
                                    + q_inc_pc2(i,j,1)
            qcl_incr_diag_conv(i,j,k)= qcl_incr_diag_conv(i,j,k)  &
                                     + qcl_inc_pc2(i,j,1)
        ! Not updated   qcf_incr_diag_conv(i,j,k)= qcf_incr_diag_conv(i,j,k)
            cf_liquid_incr_diag_conv(i,j,k)                        &
              = cf_liquid_incr_diag_conv(i,j,k) + cfl_inc_pc2(i,j,1)
        ! Not updated   cf_frozen_incr_diag_conv(i,j,k)
        !    &        = cf_frozen_incr_diag_conv(i,j,k)
            bulk_cf_incr_diag_conv(i,j,k)                          &
              = bulk_cf_incr_diag_conv(i,j,k)   + bcf_inc_pc2(i,j,1)

          END DO  ! i loop
        END DO  ! j

      END IF  ! L_q_interact_if1
    END DO  ! k

  END IF  ! L_calc_dxek_if2

! ----------------------------------------------------------------------
! Section 6 Convective history - decay of convective precip
! ----------------------------------------------------------------------

  IF (l_conv_hist) THEN

    decay_amount = timestep/decay_period
    DO  j = 1, rows
      DO i = 1, row_length
        ! total convective precip rate for the timestep
        tot_conv_precip = precip_shall(i,j) + precip_cong(i,j) +   &
                          precip_deep(i,j)  + precip_mid(i,j)
        ! If convective precip > 0.0 reset history
        IF (tot_conv_precip > 0.0) THEN
          past_precip(i,j) = tot_conv_precip
        ELSE
          ! Decay value held  - note will never reset to zero
          past_precip(i,j) = (1.0 - decay_amount)* past_precip(i,j)
        END IF
      END DO
    END DO

  END IF


! ----------------------------------------------------------------------
! Section 7 Call Convection diagnostics, also returns total
!               precipitation diagnostics.
! ----------------------------------------------------------------------
! NOTE - this has now been moved to atmos_physics2 so that the interpolation
! of wind increments can be done with other BL fields to help save CPU.   
! Only the SCM diagnostic output remains here.


!----------------------------------------------------------------------------
! Clear up allocatable arrays
!----------------------------------------------------------------------------

  IF (l_full_zero) THEN
    DEALLOCATE ( full_zero )  
  END IF

  IF (l_calc_dxek) THEN
    DEALLOCATE ( bcf_below )
    DEALLOCATE ( bcf_above )

    DEALLOCATE ( t_earliest )
    DEALLOCATE ( q_earliest )
    DEALLOCATE ( qcl_earliest )
    DEALLOCATE ( cfl_earliest )
    DEALLOCATE ( cff_earliest )
    DEALLOCATE ( bcf_earliest )
    DEALLOCATE ( t_inc_latest )

    DEALLOCATE ( theta_inc_pc2 )
    DEALLOCATE ( q_inc_pc2 )
    DEALLOCATE ( qcl_inc_pc2 )
    DEALLOCATE ( cfl_inc_pc2 )
    DEALLOCATE ( bcf_inc_pc2 )
  END IF


!----------------------------------------------------------------------------
! end of routine NI_conv_ctl
!----------------------------------------------------------------------------
END IF ! on error code equal to zero

IF (lhook) CALL dr_hook('NI_CONV_CTL',zhook_out,zhook_handle)

RETURN
END SUBROUTINE ni_conv_ctl
