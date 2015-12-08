! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine NI_bl_ctl
! *********************************************************************
!
! purpose: Interface to boundary layer mixing coefficients calculation
!   language: fortran 90 + cray extensions
!   this code is written to umdp3 programming standards.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: boundary layer
!---------------------------------------------------------------------
SUBROUTINE ni_bl_ctl (                                                  &
! IN parameters for iterative SISL scheme
 numcycles,cycleno,                                                     &
! IN model dimensions.
 land_points, ntiles, bl_levels,dst_levels, dsm_levels, nice_use,       &
! IN switches
 l_mixing_ratio, l_scrn, l_aero_classic, l_dust, l_dust_diag,           &
! IN model Parameters
 co2_mmr,                                                               &
! IN data fields.
 p, p_layer_centres, rho_rsq, rho_wet, rho_dry, u_p, v_p,               &
 u_1_px, v_1_px, u_0_px, v_0_px,                                        &
 land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels, rad_hr, &
 micro_tends, soil_layer_moisture, rho_tq, z_uv, z_tq,                  &
! IN ancillary fields and fields needed to be kept from tstep to tstep
 hcon, smvccl, smvcwt, smvcst, sthf, sthu, sil_orog_land,               &
 ho2r2_orog, sd_orog, ice_fract_cat, k_sice,                            &
 land_index, photosynth_act_rad,                                        &
! IN variables required for mineral dust scheme
 soil_clay,soil_sand,dust_mrel1,dust_mrel2,dust_mrel3,                  &
 dust_mrel4,dust_mrel5,dust_mrel6,                                      &
! IN additional variables for JULES
 canopy,catch,catch_snow,snow_tile,z0_tile,z0h_tile_bare,               &
 lw_down,sw_tile,tstar_tile,                                            &
 co2_3d,l_co2_interactive,l_phenol,l_triffid,asteps_since_triffid,      &
 cs,frac,canht_ft,lai_ft,fland,flandg,albsoil,cos_zenith_angle,         &
! IN everything not covered so far
 t_soil, ti, t_surf,zh_prev,ddmfx,bulk_cloud_fraction,nlcl, zhpar, zlcl,&
! IN SCM namelist data
 l_spec_z0, z0m_scm, z0h_scm,flux_e, flux_h, l_flux_bc,                 &
! SCM Diagnostics (dummy values in full UM) and STASH
 nscmdpkgs,l_scmdiags, BL_diag,                                         &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
! INOUT data
 gs,z0msea,w,etadot,tstar_sea,tstar_sice_cat,zh,cumulus,ntpar,l_shallow,&
 error_code,                                                            &
! INOUT additional variables for JULES
 g_leaf_acc,npp_ft_acc,resp_w_ft_acc,resp_s_acc,l_q10,                  &
! INOUT variables for TKE based turbulence schemes
 e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                          &
! OUT variables for message passing
 flandfac, fseafac, rhokm_land, rhokm_ssi, cdr10m, tau_fd_x, tau_fd_y,  &
 rhogamu, rhogamv, f_ngstress,                                          &
! OUT variables required in IMP_SOLVER
 alpha1_sice, ashtf, bq_gb, bt_gb, dtrdz_charney_grid, rdz_charney_grid,&
 dtrdz_u, dtrdz_v, rdz_u, rdz_v, z1_tq, ustargbm,                       &
! OUT diagnostics (done after implicit solver)
 e_sea, fqt, ftl, h_sea, rib_gb, vshr, zht, shallowc, cu_over_orog,     &
 bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6,           &
 bl_type_7, z0m_eff_gb, z0h_eff_gb, fme,                                &
! OUT diagnostics required for soil moisture nudging scheme :
 wt_ext,ra,                                                             &
! OUT data required for tracer mixing :
 rho_aresist,aresist,resist_b,                                          &
!OUT variables required for mineral dust scheme
 r_b_dust,dust_flux,dust_emiss_frac,                                    &
 u_s_t_tile,u_s_t_dry_tile,u_s_std_tile,kent, we_lim, t_frac, zrzi,     &
 kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,                      &
! OUT additional variables for JULES
 ftl_tile,le_tile,radnet_sice,radnet_tile,rib_tile,rho_aresist_tile,    &
 aresist_tile,resist_b_tile,alpha1,ashtf_tile,fqt_tile,epot_tile,       &
 fqt_ice,ftl_ice,fraca,resfs,resft,rhokh_tile,rhokh_sice,               &
 z0hssi,z0h_tile,z0m_gb,z0mssi,z0m_tile,chr1p5m,chr1p5m_sice,smc,       &
 gpp,npp,resp_p,g_leaf,gpp_ft,npp_ft,resp_p_ft,resp_s,resp_s_tot,       &
 resp_w_ft,gc,canhc_tile,wt_ext_tile,flake,tile_index,tile_pts,         &
 tile_frac,fsmc,rib_ssi,vshr_land,vshr_ssi,tstar_land,tstar_ssi,        &
 rhokh_mix,dtstar_tile,dtstar,hcons,emis_tile,emis_soil,                &
! OUT fields
 t1_sd, q1_sd, ntml, nbdsc, ntdsc,wstar, wthvs,uw0,vw0, rhokm,rhokh     &
 )

USE nstypes, ONLY: ntype, npft
USE atm_fields_bounds_mod, ONLY:                                        &
   udims, vdims, tdims, tdims_s, qdims, qdims_l,                        &
   pdims, pdims_s, pdims_l, wdims
USE atm_step_local, ONLY: dim_cs1, dim_cs2, land_pts_trif, npft_trif,   &
       co2_dim_len,co2_dim_row
USE bl_diags_mod, ONLY : strnewbldiag
USE switches, ONLY: l_ctile
USE level_heights_mod, ONLY:   r_theta_levels
USE earth_constants_mod, ONLY: g, earth_radius
USE dust_parameters_mod, ONLY: ndiv, ndivh
USE bl_option_mod, ONLY: l_quick_ap2
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE Submodel_Mod

IMPLICIT NONE
!---------------------------------------------------------------------
! arguments with intent in. ie: input variables.
!---------------------------------------------------------------------
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
! Parallel setup variables
INTEGER, INTENT(IN) ::                                                  &
 numcycles,                                                             &
            ! Number of phys-dyn iterations per tstep
 cycleno    ! Iteration no

! Switch for calculating exchange coeffs from latest values.
LOGICAL, INTENT(IN) ::                                                  &
 l_mixing_ratio         ! TRUE if mixing ratios used in
!                             ! boundary layer code

LOGICAL, INTENT(IN) ::                                                  &
  l_aero_classic,                                                       &
         !switch for CLASSIC aerosol scheme
  l_dust,                                                               &
         !switch for prognostic mineral dust
  l_dust_diag
         !switch for diagnostic mineral dust lifting

! Model dimensions
INTEGER, INTENT(IN) ::                                                  &
  land_points,                                                          &
              ! IN No.of land points being processed, can be 0.
  ntiles,                                                               &
              ! IN No. of land-surface tiles ( MOSES II )
  bl_levels,                                                            &
  dst_levels,                                                           &
              ! number of deep soil temperature levels
  dsm_levels,                                                           &
              ! number of deep soil moisture levels
  nice_use    ! number of sea ice categories used fully in surface
              ! calculations

LOGICAL, INTENT(IN) ::                                                  &
  l_scrn
                           ! Logical to control output
                           !    of screen level T,Q,QCL,QCF
! model parameters
REAL, INTENT(IN) ::                                                     &
  co2_mmr         ! set equal to co2_start

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) ::                                                  &
 nscmdpkgs                ! No of SCM diagnostics packages

LOGICAL, INTENT(IN) ::                                                  &
 l_scmdiags(nscmdpkgs)    ! Logicals for SCM diagnostics packages

! Data arrays
REAL, INTENT(IN) ::                                                     &
  u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
      pdims%k_start:pdims%k_end),                                       &
  v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
      pdims%k_start:pdims%k_end),                                       &
  rho_rsq(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end,                                &
          pdims_s%k_start:pdims_s%k_end),                               &
!                       ! wet density times r^2 on rho levels (kg/m3)
  rho_wet(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          1:tdims%k_end),                                               &
!                       ! wet density on rho levels (kg/m3)
  rho_dry(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          pdims%k_start:pdims%k_end),                                   &
!                       ! dry density on rho levels (kg/m3)
  z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,             &
         bl_levels),                                                    &
                              ! Z_uv(*,K) is height of half
!                                   ! level k-1/2.
  z_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,             &
         bl_levels),                                                    &
  p(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,        &
    pdims_s%k_start:pdims_s%k_end),                                     &
  p_layer_centres(tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,0:tdims%k_end),             &
  p_star(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
  theta(tdims_s%i_start:tdims_s%i_end,                                  &
        tdims_s%j_start:tdims_s%j_end,                                  &
        tdims_s%k_start:tdims_s%k_end),                                 &
  exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end),                    &
  q(qdims_l%i_start:qdims_l%i_end,qdims_l%j_start:qdims_l%j_end,        &
    qdims_l%k_start:qdims_l%k_end),                                     &
  qcl(qdims_l%i_start:qdims_l%i_end,qdims_l%j_start:qdims_l%j_end,      &
      qdims_l%k_start:qdims_l%k_end),                                   &
  qcf(qdims_l%i_start:qdims_l%i_end,qdims_l%j_start:qdims_l%j_end,      &
      qdims_l%k_start:qdims_l%k_end),                                   &
  rad_hr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
                        2,bl_levels),                                   &
                        ! IN (LW,SW) radiative heating rate (K/s)
  micro_tends(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,      &
                     2, bl_levels),                                     &
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)
  u_1_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),  &
  v_1_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),  &
  u_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),  &
  v_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

!     Declaration of new BL diagnostics.
TYPE (strnewbldiag) :: bl_diag
REAL, INTENT(IN) ::                                                     &
 soil_layer_moisture(land_points,dsm_levels)!IN soil moisture
!                 ! per layer (kg m-2)
LOGICAL, INTENT(IN) ::                                                  &
  land_sea_mask(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! ancillary arrays and fields required to be saved from tstep to tstep
INTEGER, INTENT(IN) ::                                                  &
  land_index (land_points),                                             &
                             ! set from land_sea_mask
  nlcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                              ! IN lifting condensation level
 asteps_since_triffid

REAL, INTENT(IN) ::                                                     &
  k_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use), &
                       ! sea ice surface layer effective conductivity 
                       !  (W/m2/K)
  hcon (land_points),                                                   &
                       ! soil/qrparm.soil.hcond
  smvccl (land_points,dsm_levels),                                      &
                       ! soil/qrparm.soil.crit
  smvcwt (land_points,dsm_levels),                                      &
                       ! soil/qrparm.soil.wilt
  smvcst (land_points,dsm_levels),                                      &
                       ! soil/qrparm.soil.satn
  sthf(land_points,dsm_levels),                                         &
                          ! IN Frozen soil moisture content of
                          !    each layer as a fraction of
                          !    saturation.
  sthu(land_points,dsm_levels),                                         &
                          ! IN Unfrozen soil moisture content
                          !    of each layer as a fraction of
                          !    saturation.
  sil_orog_land (land_points),                                          &
                           ! orog/qrparm.orog.as
  ho2r2_orog (land_points),                                             &
                           ! orog/qrparm.orog.h2root2
  sd_orog (land_points),                                                &
                           ! orog/qrparm.orog.stdev
  ice_fract_cat (pdims%i_start:pdims%i_end,                             &
                 pdims%j_start:pdims%j_end,nice_use),                   &
                      ! ice/qrclim.ice.(month)
                      ! If nice_use=1, this is the sum of the categories
  photosynth_act_rad(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end)
                                       ! Net downward
!                                 shortwave radiation in band 1 (w/m2).

REAL, INTENT(IN) ::                                                     &
     ! mineral dust fields
  soil_clay (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
  soil_sand (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
  dust_mrel1 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel2 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel3 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel4 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel5 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel6 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

REAL, INTENT(IN) ::                                                     &
 rho_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
!                               ! RHO_TQ(*,K) is the density at half
!                               ! level k+1/2.
  t_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
  z0m_scm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                           ! Fixed Sea surface roughness
                           ! length(m) for momentum (SCM)
  z0h_scm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                           ! Fixed Sea surface roughness
                           ! length(m) for heat (SCM)
  t_soil(land_points,dsm_levels),                                       &
                                 ! slt/qrclim.slt_pm(lev).(month)
  ti(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
                        ! set equal to tstar
  bulk_cloud_fraction(qdims%i_start:qdims%i_end,                        &
                      qdims%j_start:qdims%j_end,qdims%k_end),           &
  zh_prev(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                            ! IN boundary layer height from
!                                 !    previous timestep
  ddmfx(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
!                                 ! IN Convective downdraught
!                                 !    mass-flux at cloud base
  flux_e(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
  flux_h(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
  zhpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                               ! height of ntpar
  zlcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! height of lcl accurate value not
                               ! a model level height (m)
LOGICAL, INTENT(IN) ::                                                  &
 l_flux_bc,                                                             &
              ! T if prescribed surface fluxes to be used
 l_spec_z0,                                                             &
              ! T is roughness lengths have been specified
 l_co2_interactive,                                                     &
                             ! IN Switch for 3D CO2 field
 l_phenol,                                                              &
                             ! IN Indicates whether phenology
!                                  !    in use
l_triffid
                             ! IN Indicates whether TRIFFID
!                                  !    in use.

REAL, INTENT(IN) ::                                                     &
 canopy(land_points,ntiles),                                            &
                                ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
 catch(land_points,ntiles),                                             &
                                ! IN Surface/canopy water capacity
!                                  !    of snow-free land tiles (kg/m2).
 catch_snow(land_points,ntiles),                                        &
                             ! IN Snow interception capacity of
!                                  !    tiles (kg/m2).
 snow_tile(land_points,ntiles),                                         &
                                ! IN Lying snow on tiles (kg/m2)
 z0_tile(land_points,ntiles),                                           &
                                ! IN Tile roughness lengths (m).
 z0h_tile_bare(land_points,ntiles),                                     &
                                ! IN Tile thermal roughness lengths (m)
                                ! without snow
 lw_down(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                             ! IN Surface downward LW radiation
!                                  !    (W/m2).
 sw_tile(land_points,ntiles),                                           &
                                ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
 tstar_tile(land_points,ntiles),                                        &
                                ! IN Surface tile temperatures
 co2_3d(co2_dim_len,co2_dim_row),                                       &
!                                  ! IN 3D CO2 field if required.
 cs(land_points,dim_cs1),                                               &
                           ! IN Soil carbon (kg C/m2).
 frac(land_points,ntype),                                               &
                                ! IN Fractions of surface types.
 canht_ft(land_points,npft),                                            &
                                ! IN Canopy height (m)
 lai_ft(land_points,npft),                                              &
                                ! IN Leaf area index
 fland(land_points),                                                    &
                             ! IN Land fraction on land points.
 flandg(pdims_s%i_start:pdims_s%i_end,                                  &
        pdims_s%j_start:pdims_s%j_end),                                 &
                             ! IN Land fraction on all points.
 albsoil(land_points),                                                  &
                             ! Soil albedo.
 cos_zenith_angle(pdims%i_start:pdims%i_end,                            &
                  pdims%j_start:pdims%j_end)
!                                  ! Cosine of the zenith angle
!---------------------------------------------------------------------
! arguments with intent INOUT. ie: input variables changed on output.
!---------------------------------------------------------------------
REAL, INTENT(INOUT) ::                                                  &
  z0msea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                           ! Sea surface roughness length(m)
                           ! for momentum
  zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
  w(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,                &
    0:wdims%k_end),                                                     &
  etadot(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,           &
         0:wdims%k_end),                                                &
 tstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                             ! INOUT Open sea sfc temperature (K).
 tstar_sice_cat(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,nice_use),                    &
                             ! INOUT Sea-ice sfc temperature (K).
 gs(land_points),                                                       &
                                ! INOUT "Stomatal" conductance to
!                                  !        evaporation (m/s).
 g_leaf_acc(land_points,npft),                                          &
                                ! INOUT Accumulated G_LEAF
 npp_ft_acc(land_pts_trif,npft_trif),                                   &
!                                  ! INOUT Accumulated NPP_FT
 resp_w_ft_acc(land_pts_trif,npft_trif),                                &
!                                  ! INOUT Accum RESP_W_FT
 resp_s_acc(land_pts_trif,dim_cs2)
                                   ! INOUT Accumulated RESP_S

! INOUT variables for TKE based turbulence schemes
REAL, INTENT(INOUT) ::                                                  &
  e_trb(pdims_l%i_start:pdims_l%i_end,                                  &
        pdims_l%j_start:pdims_l%j_end,bl_levels)                        &
!                   ! TKE defined on theta levels K-1
, tsq_trb(pdims_l%i_start:pdims_l%i_end,                                &
          pdims_l%j_start:pdims_l%j_end,bl_levels)                      &
!                   ! Self covariance of liquid potential temperature
!                   ! (thetal'**2) defined on theta levels K-1
, qsq_trb(pdims_l%i_start:pdims_l%i_end,                                &
          pdims_l%j_start:pdims_l%j_end,bl_levels)                      &
!                   ! Self covariance of total water
!                   ! (qw'**2) defined on theta levels K-1
, cov_trb(pdims_l%i_start:pdims_l%i_end,                                &
          pdims_l%j_start:pdims_l%j_end,bl_levels)                      &
!                   ! Correlation between thetal and qw
!                   ! (thetal'qw') defined on theta levels K-1
, zhpar_shcu(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
!                   ! Height of mixed layer used to evaluate
!                   ! the non-gradient buoyancy flux

LOGICAL, INTENT(INOUT) ::                                               &
  cumulus (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                             ! *APL bl convection flag
  l_shallow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
!                            ! Logical indicator of shallow convection
  l_q10                       ! INOUT Indicates Q10 for soil resp'n

INTEGER, INTENT(INOUT) ::                                               &
  ntpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! IN/OUT top of diagnostic parcel
!                                   !        ascent
INTEGER, INTENT(INOUT) ::                                               &
  error_code
!---------------------------------------------------------------------
! arguments with intent OUT. ie: output variables.
!---------------------------------------------------------------------
REAL, INTENT(OUT) ::                                                    &
  t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
                          ! set to zero initially
  q1_sd(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end),           &
                          ! set to zero initially
  z0m_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
  z0h_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                          ! Effective grid-box roughness
!                                 lengths for momentum and for
!                                 heat, moisture
  uw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
!                       ! U-component of surface wind stress (P-grid)
  vw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
!                       ! V-component of surface wind stress (P-grid)
  wstar(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end),           &
                               ! surface-based mixed layer
!                                    ! velocity scale
  wthvs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                               ! surface buoyancy flux

INTEGER, INTENT(OUT) ::                                                 &
  ntml (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                              ! OUT
  nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
  ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
  kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                              ! OUT grid-level of SML inversion
  kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! OUT grid-level of DSC inversion
REAL, INTENT(OUT) ::                                                    &
  rhokm(pdims_s%i_start:pdims_s%i_end,                                  &
        pdims_s%j_start:pdims_s%j_end,0:bl_levels-1)

! variables passed from BDY_LAYR to IMP_SOLVER
REAL, INTENT(OUT) ::                                                    &
  alpha1_sice(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,nice_use),                      &
  ashtf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use),  &
  dtrdz_charney_grid(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end,bl_levels),              &
  rdz_charney_grid(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,bl_levels),                &
  dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,          &
                          bl_levels),                                   &
  dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,          &
                            bl_levels),                                 &
  bq_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
  bt_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
  rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,            &
                        2:bl_levels),                                   &
  rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,            &
                          2:bl_levels),                                 &
  z1_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
  flandfac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),&
  fseafac(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end), &
  rhokm_land(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end),                            &
  rhokm_ssi(pdims_s%i_start:pdims_s%i_end,                              &
            pdims_s%j_start:pdims_s%j_end),                             &
  cdr10m(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),  &
  tau_fd_x(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end, &
           bl_levels),                                                  &
  tau_fd_y(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end, &
           bl_levels),                                                  &
  rhogamu(pdims_s%i_start:pdims_s%i_end,                                &
           pdims_s%j_start:pdims_s%j_end ,bl_levels),                   &
  rhogamv(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end ,bl_levels),                    &
  f_ngstress(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end,2:bl_levels)

REAL, INTENT(OUT) :: ustargbm(pdims%i_start:pdims%i_end,                &
                              pdims%j_start:pdims%j_end)
!       ! GBM surface friction velocity for diagnosis of decoupling
! Diagnostics needed in NI_imp_ctl
REAL, INTENT(OUT) ::                                                    &
  fme(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
  rho_aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                          ! OUT RHOSTAR*CD_STD*VSHR
                          !     for CLASSIC aerosol scheme
  aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                          ! OUT 1/(CD_STD*VSHR)
                          !     for CLASSIC aerosol scheme
  resist_b(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                          ! OUT (1/CH-1/(CD_STD)/VSHR
                          !     for CLASSIC aerosol scheme
  r_b_dust(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           ndiv),                                                       &
                                !OUT surface layer resist for dust
  dust_flux(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            ndiv),                                                      &
                                 !OUT dust emissions (kg m-2 s-1)
  dust_emiss_frac(land_points,ntiles),                                  &
                          ! OUT fraction of tile can emit dust
  u_s_t_tile(land_points,ntiles,ndivh),                                 &
                                     !OUT threshold frict. vel
  u_s_t_dry_tile(land_points,ntiles,ndivh),                             &
                                         !OUT dry soil value
  u_s_std_tile(land_points,ntiles),                                     &
                                  !OUT friction velocity
  we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                              ! OUT rho*entrainment rate implied b
!                                   !     placing of subsidence
  zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),          &
                              ! OUT (z-z_base)/(z_i-z_base)
  t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                              ! OUT a fraction of the timestep
  we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),    &
!                                   ! OUT rho*entrainment rate implied b
!                                   !     placing of subsidence
  zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                              ! OUT (z-z_base)/(z_i-z_base)
  t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),    &
!                                   ! OUT a fraction of the timestep
  zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                              ! OUT Top of decoupled layer
  zht(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                                ! Max height of turb mixing
  shallowc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                             !OUT Indicator set to 1.0 if shallow,
!                                  !   0.0 if not shallow or not cumulus
  cu_over_orog(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
!                                  ! OUT Indicator for cumulus
!                                  !     over steep orography
!                                  !   Indicator set to 1.0 if true,
!                                  !   0.0 if false. Exclusive.
  bl_type_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                            ! OUT Indicator set to 1.0 if stable
!                                  !     b.l. diagnosed, 0.0 otherwise.
  bl_type_2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                            ! OUT Indicator set to 1.0 if Sc over
!                                  !     stable surface layer diagnosed,
!                                  !     0.0 otherwise.
  bl_type_3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                            ! OUT Indicator set to 1.0 if well
!                                  !     mixed b.l. diagnosed,
!                                  !     0.0 otherwise.
  bl_type_4(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                            ! OUT Indicator set to 1.0 if
!                                  !     decoupled Sc layer (not over
!                                  !     cumulus) diagnosed,
!                                  !     0.0 otherwise.
  bl_type_5(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                            ! OUT Indicator set to 1.0 if
!                                  !     decoupled Sc layer over cumulus
!                                  !     diagnosed, 0.0 otherwise.
  bl_type_6(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                            ! OUT Indicator set to 1.0 if a
!                                  !     cumulus capped b.l. diagnosed,
!                                  !     0.0 otherwise.
  bl_type_7(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
!                                  !     shear-dominated  b.l.
!                                  !     diagnosed, 0.0 otherwise.
  vshr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
  ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
      bl_levels),                                                       &
                                       ! needed as diagnostic
  fqt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
      bl_levels),                                                       &
                                       ! needed as diagnostic ?
  h_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                       ! needed as diagnostic ?
  e_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                       ! needed as diagnostic ?
  rib_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                          ! Mean bulk Richardson number for
                               !  lowest layer.
  rhokh (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         bl_levels),                                                    &
  wt_ext(land_points,dsm_levels),                                       &
                                !OUT cumulative fract of transp'n
  ra(land_points),                                                      &
                                !OUT Aerodynamic resistance (s/m)
 tstar_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                             ! OUT Land mean sfc temperature (K)

 tstar_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                             ! OUT Sea mean sfc temperature (K).

INTEGER, INTENT(OUT) ::                                                 &
 tile_index(land_points,ntype),                                         &
                                ! OUT Index of tile points
 tile_pts(ntype)             ! OUT Number of tile points

REAL, INTENT(OUT) ::                                                    &
 ftl_tile(land_points,ntiles),                                          &
                                ! OUT Surface FTL for land tiles
 le_tile(land_points,ntiles),                                           &
                                ! OUT Surface latent heat flux for
!                                  !     land tiles
 radnet_sice(pdims%i_start:pdims%i_end,                                 &
             pdims%j_start:pdims%j_end,nice_use),                       &
                             ! OUT Surface net radiation on
!                                  !     sea-ice (W/m2)
 radnet_tile(land_points,ntiles),                                       &
                                ! OUT Surface net radiation on
!                                  !     land tiles (W/m2)
 rib_tile(land_points,ntiles),                                          &
                                ! OUT RIB for land tiles.
 rho_aresist_tile(land_points,ntiles),                                  &
!                                  ! OUT RHOSTAR*CD_STD*VSHR on land
!                                  !     tiles for CLASSIC aerosol scheme
 aresist_tile(land_points,ntiles),                                      &
!                                  ! OUT 1/(CD_STD*VSHR) on land tiles
                                   !     for CLASSIC aerosol scheme
 resist_b_tile(land_points,ntiles),                                     &
!                                  ! OUT (1/CH-1/CD_STD)/VSHR on land
!                                  !     tiles for CLASSIC aerosol scheme
 alpha1(land_points,ntiles),                                            &
                                ! OUT Mean gradient of saturated
!                                  !     specific humidity with respect
!                                  !     to temperature between the
!                                  !     bottom model layer and tile
!                                  !     surfaces
 ashtf_tile(land_points,ntiles),                                        &
                                !OUT Coefficient to calculate
!                                  !     surface heat flux into land
!                                  !     tiles.
 fqt_tile(land_points,ntiles),                                          &
                                ! OUT Surface FQT for land tiles
 epot_tile(land_points,ntiles),                                         &
                                ! OUT Local EPOT for land tiles.
 fqt_ice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use), &
                             ! OUT Surface FQT for sea-ice
 ftl_ice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use), &
                             ! OUT Surface FTL for sea-ice
 fraca(land_points,ntiles),                                             &
                                ! OUT Fraction of surface moisture
!                                  !     flux with only aerodynamic
!                                  !     resistance for snow-free land
!                                  !     tiles.
 resfs(land_points,ntiles),                                             &
                                ! OUT Combined soil, stomatal
!                                  !     and aerodynamic resistance
!                                  !     factor for fraction (1-FRACA)
!                                  !     of snow-free land tiles.
 resft(land_points,ntiles),                                             &
                                ! OUT Total resistance factor.
!                                  !     FRACA+(1-FRACA)*RESFS for
!                                  !     snow-free land, 1 for snow.
 rhokh_tile(land_points,ntiles),                                        &
                                ! OUT Surface exchange coefficient
!                                  !     for land tiles
 rhokh_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                             ! OUT Surface exchange coefficients
!                                  !     for sea and sea-ice
 z0hssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
 z0mssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                             ! OUT Roughness lengths over sea (m).
 z0h_tile(land_points,ntiles),                                          &
                             ! OUT Tile roughness lengths for h
!                                  !     and moisture (m).
 z0m_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                             ! OUT Grid-box mean roughness length
!                                  !      for momentum (m).
 z0m_tile(land_points,ntiles),                                          &
                             ! OUT Tile roughness lengths for
!                                  !     momentum.
 chr1p5m(land_points,ntiles),                                           &
                             ! OUT Ratio of coefffs for
!                                  !     calculation of 1.5m temp for
!                                  !     land tiles.
 chr1p5m_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
!                                  ! OUT CHR1P5M for sea and sea-ice
!                                  !     (leads ignored).
 smc(land_points),                                                      &
                                ! OUT Available moisture in the
!                                  !     soil profile (mm).
 gpp(land_points),                                                      &
                                ! OUT Gross primary productivity
!                                  !     (kg C/m2/s).
 npp(land_points),                                                      &
                                ! OUT Net primary productivity
!                                  !     (kg C/m2/s).
 resp_p(land_points),                                                   &
!CABLE: kdcorbin, 11/10 - changed from NPFT, currently setup for Leaf Resp
 g_leaf(land_points,ntiles),                                            &
                                ! OUT Leaf turnover rate (/360days
 gpp_ft(land_points,ntiles),                                            &
                                ! OUT Gross primary productivity
!                                  !     on PFTs (kg C/m2/s).
 npp_ft(land_points,ntiles),                                            &
                                ! OUT Net primary productivity
!                                  !     (kg C/m2/s).
 resp_p_ft(land_points,ntiles),                                           &
                                ! OUT Plant respiration on PFTs
!                                  !     (kg C/m2/s).
 resp_s(land_points,dim_cs1),                                           &
                            ! OUT Soil respiration (kg C/m2/s)
 resp_s_tot(dim_cs2),                                                   &
                           ! OUT total soil respiration
 resp_w_ft(land_points,npft),                                           &
                                ! OUT Wood maintenance respiration
!                                  !     (kg C/m2/s).
 gc(land_points,ntiles),                                                &
                                ! OUT "Stomatal" conductance to
!                                  !      evaporation for land tiles
!                                  !      (m/s).
 canhc_tile(land_points,ntiles),                                        &
                                ! OUT Areal heat capacity of canopy
!                                  !    for land tiles (J/K/m2).
 wt_ext_tile(land_points,dsm_levels,ntiles),                            &
!                                  ! OUT Fraction of evapotranspiration
!                                  !    which is extracted from each
!                                  !    soil layer by each tile.
 flake(land_points,ntiles),                                             &
                                ! OUT Lake fraction.
 tile_frac(land_points,ntiles),                                         &
                                ! OUT Tile fractions including
!                                  !     snow cover in the ice tile.
 fsmc(land_points,npft),                                                &
                                ! OUT Moisture availability
!                                  !     factor.
 rib_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                             ! OUT Sea mean bulk Richardson number
                             !     for lowest layer.
  vshr_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
  vshr_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! Additional variables for JULES
REAL, INTENT(OUT) ::                                                    &
 rhokh_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                             ! Exchange coeffs for moisture.
 dtstar_tile(land_points,ntiles),                                       &
                             ! Change in TSTAR over timestep
!                                  ! for land tiles
 dtstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use),  &
                             ! Change is TSTAR over timestep
!                                  ! for sea-ice
 hcons(land_points),                                                    &
                             ! Soil thermal conductivity
!                                  ! including water and ice
 emis_tile(land_points,ntiles),                                         &
                             ! Emissivity for land tiles
 emis_soil(land_points)
                             ! Emissivity of underlying soil
!---------------------------------------------------------------------
! local variables.
!---------------------------------------------------------------------
! loop counters
INTEGER                                                                 &
  i, j, k ,l

! Diagnostic switches
LOGICAL                                                                 &
  su10, sv10, slh, sq1p5, st1p5, sq_t1p5                                &
, sfme, sz0heff, l_apply_diag

REAL                                                                    &
  t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels)

REAL                                                                    &
 z_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
!                       !land height over fractional land points

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NI_BL_CTL',zhook_in,zhook_handle)

! error information
IF ( error_code  ==  0) THEN
! ----------------------------------------------------------------------
! Section BL.0 Initialisation of variables.
! ----------------------------------------------------------------------
! Apply diags at last cycle only or if l_quick_ap2 is .true.
  l_apply_diag = (cycleno == numcycles .OR. l_quick_ap2)
! Set diagnostic flags required for boundary layer diagnostics from
! STASHflags.
!        ! --------------------------------------------------------
!        ! Note that an equivalent block of code exists in routine
!        ! ni_imp_ctl, and needs to be kept consistent.
!        ! --------------------------------------------------------
!        ! Windspeed (227, 230) and u, v at 10m on 'B' or 'C' grid
   su10 = (sf(209,3) .OR. sf(225,3) .OR. sf(227,3) .OR. sf(230,3) .OR.  &
           sf(463,3)) .AND. l_apply_diag
   sv10 = (sf(210,3) .OR. sf(226,3) .OR. sf(227,3) .OR. sf(230,3) .OR.  &
           sf(463,3)) .AND. l_apply_diag
   slh = sf(234,3) .AND. l_apply_diag
   sq_t1p5 = ( sf(236,3) .OR. sf(237,3) .OR. sf(245,3)                  &
        .OR. sf(247,3) .OR. sf(248,3) .OR. sf(250,3)                    &
        .OR. l_scrn                                                     &
        .OR. sf(341,3) .OR. sf(342,3)                                   &
        .OR. sf(253,3) .OR. sf(328,3) .OR. sf(329,3)                    &
              ) .AND. l_apply_diag
   sq1p5 = sq_t1p5 .AND. l_apply_diag
   st1p5 = sq_t1p5 .AND. l_apply_diag
   sfme  = sf(224,3) .AND. l_apply_diag
   sz0heff = sf(027,3) .AND. l_apply_diag
!       Set switches for BL diagnostic arrays
! counter gradient term for u
  bl_diag%l_rhogamu      = sf(130,3) .AND. l_apply_diag
! counter gradient term for v
  bl_diag%l_rhogamv      = sf(131,3) .AND. l_apply_diag
! counter gradient term for t
  bl_diag%l_rhogamt      = sf(132,3) .AND. l_apply_diag
! counter gradient term for q
  bl_diag%l_rhogamq      = sf(133,3) .AND. l_apply_diag
! mixing length
  bl_diag%l_elm          = sf(134,3) .AND. l_apply_diag
! production rate of TKE by shear
  bl_diag%l_tke_shr_prod = sf(135,3) .AND. l_apply_diag
! production rate of TKE by buoyancy
  bl_diag%l_tke_boy_prod = sf(136,3) .AND. l_apply_diag
! dissipation rate of TKE
  bl_diag%l_tke_dissp    = sf(137,3) .AND. l_apply_diag
! non-dimensional diffusion coef. for u, v
  bl_diag%l_sm           = sf(138,3) .AND. l_apply_diag
! non-dimensional diffusion coef. for t, q
  bl_diag%l_sh           = sf(139,3) .AND. l_apply_diag
! non-gradient buoyancy flux
  bl_diag%l_wb_ng        = sf(140,3) .AND. l_apply_diag
! cloud fraction used in the TKE schemes
  bl_diag%l_cf_trb       = sf(141,3) .AND. l_apply_diag
! condensed water used in the TKE schemes
  bl_diag%l_ql_trb       = sf(142,3) .AND. l_apply_diag
! standard deviation of the distribution function in the TKE schemes
  bl_diag%l_sgm_trb      = sf(143,3) .AND. l_apply_diag
! Heating increment from turbulence dissipation
  bl_diag%l_dtfric     = sf(188,3) .AND. l_apply_diag
! Top of surface mixed layer (Ksurf profile)
  bl_diag%l_smltop     = sf(356,3) .AND. l_apply_diag
! Top of decoupled stratocu layer
  bl_diag%l_dsctop     = sf(357,3) .AND. l_apply_diag
! BL depth diagnosed from Ri>RiCrit
  bl_diag%l_zhlocal    = sf(358,3) .AND. l_apply_diag
! Height of diagnosis parcel top
  bl_diag%l_zhpar      = sf(359,3) .AND. l_apply_diag
! Decoupled stratocu base height, also needed for TKE diagnostic
  bl_diag%l_dscbase    = ( sf(360,3) .OR. sf(473,3) )                   &
                          .AND. l_apply_diag
! BL cloud base height
  bl_diag%l_cldbase    = sf(361,3) .AND. l_apply_diag
! Entrainment rate
  bl_diag%l_weparm     = sf(362,3) .AND. l_apply_diag
! Entrainment rate for decoupled stratocu
  bl_diag%l_weparm_dsc = sf(363,3) .AND. l_apply_diag
! Obukhov length, also required for gustiness diagnostic (463)
  bl_diag%l_oblen      = ( sf(464,3) .OR. sf(463,3) )                   &
                          .AND. l_apply_diag
! Friction velocity, also required for gustiness diagnostic (463)
  bl_diag%l_ustar      = ( sf(465,3) .OR. sf(463,3) )                   &
                          .AND. l_apply_diag
! Surface buoyancy flux
  bl_diag%l_wbsurf     = sf(467,3) .AND. l_apply_diag
! Gradient Richardson number
  bl_diag%l_gradrich   = sf(468,3) .AND. l_apply_diag
! Convective velocity scale
  bl_diag%l_wstar      = sf(466,3) .AND. l_apply_diag
! Stratification
  bl_diag%l_dbdz       = sf(469,3) .AND. l_apply_diag
! Modulus of shear
  bl_diag%l_dvdzm      = sf(470,3) .AND. l_apply_diag
! Momentum diffusivity
  bl_diag%l_rhokm      = sf(471,3) .AND. l_apply_diag
! Thermal diffusivity
  bl_diag%l_rhokh      = sf(472,3) .AND. l_apply_diag
! Turbulent kinetic energy
  bl_diag%l_tke        = sf(473,3) .AND. l_apply_diag
! x component of orographic stress
  bl_diag%l_ostressx   = sf(474,3) .AND. l_apply_diag
! y component of orographic stress
  bl_diag%l_ostressy   = sf(475,3) .AND. l_apply_diag
! local mixing length for momentum
  bl_diag%l_elm3d      = sf(501,3) .AND. l_apply_diag
! local mixing length for scalars
  bl_diag%l_elh3d      = sf(502,3) .AND. l_apply_diag
! local momentum diffusion coefficient
  bl_diag%l_rhokmloc   = sf(503,3) .AND. l_apply_diag
! local scalar diffusion coefficient
  bl_diag%l_rhokhloc   = sf(504,3) .AND. l_apply_diag
! surface driven momentum diffusion coefficient
  bl_diag%l_rhokmsurf  = sf(505,3) .AND. l_apply_diag
! surface driven scalar diffusion coefficient
  bl_diag%l_rhokhsurf  = sf(506,3) .AND. l_apply_diag
! stratocu-top-driven momentum diffusion coefficient
  bl_diag%l_rhokmsc    = sf(507,3) .AND. l_apply_diag
! stratocu-top-driven scalar diffusion coefficient
  bl_diag%l_rhokhsc    = sf(508,3) .AND. l_apply_diag
! weighting applied to 1D BL scheme in Smag blending
  bl_diag%l_weight1d   = sf(509,3) .AND. l_apply_diag

!       Allocate space for those BL diagnostic arrays required and zero
!       the elements explicitly
  IF (bl_diag%l_oblen) THEN
    ALLOCATE(bl_diag%oblen(                                             &
             pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
    bl_diag%oblen(:,:) = 0.0
  END IF
  IF (bl_diag%l_ustar) THEN
    ALLOCATE(bl_diag%ustar(                                             &
             pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
    bl_diag%ustar(:,:) = 0.0
  END IF
  IF (bl_diag%l_wbsurf) THEN
    ALLOCATE(bl_diag%wbsurf(                                            &
             pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
    bl_diag%wbsurf(:,:) = 0.0
  END IF
  IF (bl_diag%l_gradrich) THEN
    ALLOCATE(bl_diag%gradrich(pdims%i_start:pdims%i_end,                &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%gradrich(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_wstar) THEN
    ALLOCATE(bl_diag%wstar(                                             &
             pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
    bl_diag%wstar(:,:) = 0.0
  END IF
  IF (bl_diag%l_dbdz) THEN
    ALLOCATE(bl_diag%dbdz(pdims%i_start:pdims%i_end,                    &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%dbdz(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_dvdzm) THEN
    ALLOCATE(bl_diag%dvdzm(pdims%i_start:pdims%i_end,                   &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%dvdzm(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_rhokm) THEN
    ALLOCATE(bl_diag%rhokm(pdims%i_start:pdims%i_end,                   &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%rhokm(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_rhokh) THEN
    ALLOCATE(bl_diag%rhokh(pdims%i_start:pdims%i_end,                   &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%rhokh(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_tke) THEN
    ALLOCATE(bl_diag%tke(pdims%i_start:pdims%i_end,                     &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%tke(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_ostressx) THEN
    ALLOCATE(bl_diag%ostressx(pdims%i_start:pdims%i_end,                &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%ostressx(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_ostressy) THEN
    ALLOCATE(bl_diag%ostressy(pdims%i_start:pdims%i_end,                &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%ostressy(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_smltop) THEN
    ALLOCATE(bl_diag%smltop(                                            &
             pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
    bl_diag%smltop(:,:) = 0.0
  END IF
  IF (bl_diag%l_dsctop) THEN
    ALLOCATE(bl_diag%dsctop(                                            &
             pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
    bl_diag%dsctop(:,:) = 0.0
  END IF
  IF (bl_diag%l_zhlocal) THEN
    ALLOCATE(bl_diag%zhlocal(                                           &
             pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
    bl_diag%zhlocal(:,:) = 0.0
  END IF
  IF (bl_diag%l_zhpar) THEN
    ALLOCATE(bl_diag%zhpar(                                             &
             pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
    bl_diag%zhpar(:,:) = 0.0
  END IF
  IF (bl_diag%l_dscbase) THEN
    ALLOCATE(bl_diag%dscbase(                                           &
             pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
    bl_diag%dscbase(:,:) = 0.0
  END IF
  IF (bl_diag%l_cldbase) THEN
    ALLOCATE(bl_diag%cldbase(                                           &
             pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
    bl_diag%cldbase(:,:) = 0.0
  END IF
  IF (bl_diag%l_weparm) THEN
    ALLOCATE(bl_diag%weparm(                                            &
             pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
    bl_diag%weparm(:,:) = 0.0
  END IF
  IF (bl_diag%l_weparm_dsc) THEN
    ALLOCATE(bl_diag%weparm_dsc(                                        &
             pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
    bl_diag%weparm_dsc(:,:) = 0.0
  END IF
  IF (bl_diag%l_dtfric) THEN
    ALLOCATE(bl_diag%dtfric(pdims%i_start:pdims%i_end,                 &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%dtfric(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_elm3d) THEN
    ALLOCATE(bl_diag%elm3d(pdims%i_start:pdims%i_end,                 &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%elm3d(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_elh3d) THEN
    ALLOCATE(bl_diag%elh3d(pdims%i_start:pdims%i_end,                 &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%elh3d(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_rhokmloc) THEN
    ALLOCATE(bl_diag%rhokmloc(pdims%i_start:pdims%i_end,              &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%rhokmloc(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_rhokhloc) THEN
    ALLOCATE(bl_diag%rhokhloc(pdims%i_start:pdims%i_end,              &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%rhokhloc(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_rhokmsurf) THEN
    ALLOCATE(bl_diag%rhokmsurf(pdims%i_start:pdims%i_end,             &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%rhokmsurf(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_rhokhsurf) THEN
    ALLOCATE(bl_diag%rhokhsurf(pdims%i_start:pdims%i_end,             &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%rhokhsurf(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_rhokmsc) THEN
    ALLOCATE(bl_diag%rhokmsc(pdims%i_start:pdims%i_end,               &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%rhokmsc(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_rhokhsc) THEN
    ALLOCATE(bl_diag%rhokhsc(pdims%i_start:pdims%i_end,               &
                             pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%rhokhsc(:,:,:) = 0.0
  END IF
  IF (bl_diag%l_weight1d) THEN
    ALLOCATE(bl_diag%weight1d(pdims%i_start:pdims%i_end,              &
                              pdims%j_start:pdims%j_end,bl_levels))
    bl_diag%weight1d(:,:,:) = 0.0
  END IF

! Code to allocate unity arrays when not
! used (for portability)
  IF (.NOT. bl_diag%l_oblen) THEN
    ALLOCATE(bl_diag%oblen(1,1))
    bl_diag%oblen(:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_ustar) THEN
    ALLOCATE(bl_diag%ustar(1,1))
    bl_diag%ustar(:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_wbsurf) THEN
    ALLOCATE(bl_diag%wbsurf(1,1))
    bl_diag%wbsurf(:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_gradrich) THEN
    ALLOCATE(bl_diag%gradrich(1,1,1))
    bl_diag%gradrich(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_wstar) THEN
    ALLOCATE(bl_diag%wstar(1,1))
    bl_diag%wstar(:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_dbdz) THEN
    ALLOCATE(bl_diag%dbdz(1,1,1))
    bl_diag%dbdz(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_dvdzm) THEN
    ALLOCATE(bl_diag%dvdzm(1,1,1))
    bl_diag%dvdzm(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_rhokm) THEN
    ALLOCATE(bl_diag%rhokm(1,1,1))
    bl_diag%rhokm(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_rhokh) THEN
    ALLOCATE(bl_diag%rhokh(1,1,1))
    bl_diag%rhokh(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_tke) THEN
    ALLOCATE(bl_diag%tke(1,1,1))
    bl_diag%tke(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_ostressx) THEN
    ALLOCATE(bl_diag%ostressx(1,1,1))
    bl_diag%ostressx(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_ostressy) THEN
    ALLOCATE(bl_diag%ostressy(1,1,1))
    bl_diag%ostressy(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_smltop) THEN
    ALLOCATE(bl_diag%smltop(1,1))
    bl_diag%smltop(:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_dsctop) THEN
    ALLOCATE(bl_diag%dsctop(1,1))
    bl_diag%dsctop(:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_zhlocal) THEN
    ALLOCATE(bl_diag%zhlocal(1,1))
    bl_diag%zhlocal(:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_zhpar) THEN
    ALLOCATE(bl_diag%zhpar(1,1))
    bl_diag%zhpar(:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_dscbase) THEN
    ALLOCATE(bl_diag%dscbase(1,1))
    bl_diag%dscbase(:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_cldbase) THEN
    ALLOCATE(bl_diag%cldbase(1,1))
    bl_diag%cldbase(:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_weparm) THEN
    ALLOCATE(bl_diag%weparm(1,1))
    bl_diag%weparm(:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_weparm_dsc) THEN
    ALLOCATE(bl_diag%weparm_dsc(1,1))
    bl_diag%weparm_dsc(:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_dtfric) THEN
    ALLOCATE(bl_diag%dtfric(1,1,1))
    bl_diag%dtfric(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_elm3d) THEN
    ALLOCATE(bl_diag%elm3d(1,1,1))
    bl_diag%elm3d(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_elh3d) THEN
    ALLOCATE(bl_diag%elh3d(1,1,1))
    bl_diag%elh3d(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_rhokmloc) THEN
    ALLOCATE(bl_diag%rhokmloc(1,1,1))
    bl_diag%rhokmloc(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_rhokhloc) THEN
    ALLOCATE(bl_diag%rhokhloc(1,1,1))
    bl_diag%rhokhloc(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_rhokmsurf) THEN
    ALLOCATE(bl_diag%rhokmsurf(1,1,1))
    bl_diag%rhokmsurf(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_rhokmsc) THEN
    ALLOCATE(bl_diag%rhokmsc(1,1,1))
    bl_diag%rhokmsc(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_rhokhsurf) THEN
    ALLOCATE(bl_diag%rhokhsurf(1,1,1))
    bl_diag%rhokhsurf(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_rhokhsc) THEN
    ALLOCATE(bl_diag%rhokhsc(1,1,1))
    bl_diag%rhokhsc(:,:,:) = 0.0
  END IF
  IF (.NOT. bl_diag%l_weight1d) THEN
    ALLOCATE(bl_diag%weight1d(1,1,1))
    bl_diag%weight1d(:,:,:) = 0.0
  END IF


    IF (bl_diag%l_rhogamu) THEN
      ALLOCATE(bl_diag%rhogamu(pdims%i_start:pdims%i_end,              &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%rhogamu(1,1,1))
    END IF

    IF (bl_diag%l_rhogamv) THEN
      ALLOCATE(bl_diag%rhogamv(pdims%i_start:pdims%i_end,              &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%rhogamv(1,1,1))
    END IF

    IF (bl_diag%l_rhogamt) THEN
      ALLOCATE(bl_diag%rhogamt(pdims%i_start:pdims%i_end,              &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%rhogamt(1,1,1))
    END IF

    IF (bl_diag%l_rhogamq) THEN
      ALLOCATE(bl_diag%rhogamq(pdims%i_start:pdims%i_end,              &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%rhogamq(1,1,1))
    END IF

    IF (bl_diag%l_elm) THEN
      ALLOCATE(bl_diag%elm(pdims%i_start:pdims%i_end,                  &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%elm(1,1,1))
    END IF

    IF (bl_diag%l_tke_shr_prod) THEN
      ALLOCATE(bl_diag%tke_shr_prod(pdims%i_start:pdims%i_end,         &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%tke_shr_prod(1,1,1))
    END IF

    IF (bl_diag%l_tke_boy_prod) THEN
      ALLOCATE(bl_diag%tke_boy_prod(pdims%i_start:pdims%i_end,         &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%tke_boy_prod(1,1,1))
    END IF

    IF (bl_diag%l_tke_dissp) THEN
      ALLOCATE(bl_diag%tke_dissp(pdims%i_start:pdims%i_end,            &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%tke_dissp(1,1,1))
    END IF

    IF (bl_diag%l_sm) THEN
      ALLOCATE(bl_diag%sm(pdims%i_start:pdims%i_end,                   &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%sm(1,1,1))
    END IF

    IF (bl_diag%l_sh) THEN
      ALLOCATE(bl_diag%sh(pdims%i_start:pdims%i_end,                   &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%sh(1,1,1))
    END IF

    IF (bl_diag%l_wb_ng) THEN
      ALLOCATE(bl_diag%wb_ng(pdims%i_start:pdims%i_end,                &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%wb_ng(1,1,1))
    END IF

    IF (bl_diag%l_cf_trb) THEN
      ALLOCATE(bl_diag%cf_trb(pdims%i_start:pdims%i_end,               &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%cf_trb(1,1,1))
    END IF

    IF (bl_diag%l_ql_trb) THEN
      ALLOCATE(bl_diag%ql_trb(pdims%i_start:pdims%i_end,               &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%ql_trb(1,1,1))
    END IF

    IF (bl_diag%l_sgm_trb) THEN
      ALLOCATE(bl_diag%sgm_trb(pdims%i_start:pdims%i_end,              &
                            pdims%j_start:pdims%j_end,bl_levels))
    ELSE
      ALLOCATE(bl_diag%sgm_trb(1,1,1))
    END IF

  bl_diag%rhogamu(:,:,:) = 0.0
  bl_diag%rhogamv(:,:,:) = 0.0
  bl_diag%rhogamt(:,:,:) = 0.0
  bl_diag%rhogamq(:,:,:) = 0.0
  bl_diag%elm(:,:,:) = 0.0
  bl_diag%tke_shr_prod(:,:,:) = 0.0
  bl_diag%tke_boy_prod(:,:,:) = 0.0
  bl_diag%tke_dissp(:,:,:) = 0.0
  bl_diag%sm(:,:,:) = 0.0
  bl_diag%sh(:,:,:) = 0.0
  bl_diag%wb_ng(:,:,:) = 0.0
  bl_diag%cf_trb(:,:,:) = 0.0
  bl_diag%ql_trb(:,:,:) = 0.0
  bl_diag%sgm_trb(:,:,:) = 0.0
! ----------------------------------------------------------------------
! Section BL.1 Calculate T at old time level.
! Modified to use latest values to avoid time-level inconsistencies
! with cloud data.
! ---------------------------------------------------------------------
   DO k = 1, bl_levels
      DO j = tdims%j_start, tdims%j_end
         DO i = tdims%i_start, tdims%i_end
            t(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
         END DO
      END DO
   END DO
! ----------------------------------------------------------------------
! Calculate the land height over fractional land points:
! ----------------------------------------------------------------------
  l=0
  DO j = pdims%j_start, pdims%j_end
     DO i = pdims%i_start, pdims%i_end
       z_land(i,j) = 0.0

       IF(land_sea_mask(i,j))THEN
         l=l+1

         IF(l_ctile.AND.fland(l) >  0.0.AND.fland(l) <  1.0)THEN
           z_land(i,j) = r_theta_levels(i,j,0) - earth_radius
           IF(z_land(i,j) <  0.0)z_land(i,j) = 0.0
         END IF

       END IF
     END DO
  END DO
! ----------------------------------------------------------------------
! Section BL.2a Call boundary_layer scheme.
! ----------------------------------------------------------------------
! DEPENDS ON: bdy_layr
   CALL bdy_layr(                                                       &
! IN  parameters for iterative SISL scheme
  numcycles, cycleno,                                                   &
! IN values defining field dimensions and subset to be processed :
  ntiles,land_points,nice_use,bl_levels,                                &
! IN values defining vertical grid of model atmosphere :
  p,p_layer_centres(1,1,1),rho_rsq,rho_wet,rho_dry,rho_tq, z_uv, z_tq,  &
! IN U and V momentum fields.
  u_p, v_p, u_1_px, v_1_px, u_0_px, v_0_px,                             &
! IN soil/vegetation/land surface data :
  land_index,dst_levels,dsm_levels,canopy,catch,catch_snow,hcon,smvccl, &
  smvcst,smvcwt,sthf,sthu,sil_orog_land,ho2r2_orog,sd_orog,             &
! IN for dust scheme
  soil_layer_moisture,                                                  &
! IN sea/sea-ice data :
  ice_fract_cat, k_sice,                                                &
! IN cloud data :
  bulk_cloud_fraction,                                                  &
  q(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,1:bl_levels),   &
  qcf(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,1:bl_levels), &
  qcl(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,1:bl_levels), &
  t,                                                                    &
! IN everything not covered so far :
  co2_mmr, photosynth_act_rad, p_star,                                  &
  rad_hr,micro_tends,l_mixing_ratio,zh_prev,ddmfx,nlcl,zhpar,zlcl,      &
! IN Variables for: prescribed surface flux forcing
  flux_e, flux_h, l_flux_bc, l_spec_z0, z0m_scm, z0h_scm,               &
! IN variables required for CLASSIC/mineral dust schemes
  l_aero_classic,l_dust,l_dust_diag,soil_clay,soil_sand,                &
  dust_mrel1,dust_mrel2,dust_mrel3,dust_mrel4,dust_mrel5,dust_mrel6,    &
! IN additional variables for JULES
  snow_tile,z0_tile,z0h_tile_bare,lw_down,                              &
  sw_tile,t_soil,ti,t_surf,tstar_tile,                                  &
  co2_3d,l_co2_interactive,l_phenol,l_triffid,asteps_since_triffid,     &
  cs,frac,canht_ft,lai_ft,fland,flandg,z_land,albsoil,cos_zenith_angle, &
! IN stash flags :-
  sfme,  slh, sq1p5, st1p5, su10, sv10, sz0heff,                        &
! SCM Diagnostics (dummy values in full UM) and STASH
  nscmdpkgs,l_scmdiags, BL_diag,                                        &
! INOUT data
  gs,z0msea,w, etadot,tstar_sea,tstar_sice_cat,                         &
  zh,ntpar,l_shallow,cumulus,error_code,                                &
! INOUT additional variables for JULES
  g_leaf_acc,npp_ft_acc,resp_w_ft_acc,resp_s_acc, l_q10,                &
! INOUT variables for TKE based turbulence schemes
  e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                         &
! OUT variables for message passing
  flandfac, fseafac, rhokm_land, rhokm_ssi, cdr10m, tau_fd_x, tau_fd_y, &
  rhogamu, rhogamv, f_ngstress,                                         &
! OUT  diagnostic not requiring stash flags :
  e_sea, fqt, fqt_tile, epot_tile, ftl,ftl_tile, h_sea, rhokh, rhokm,   &
  rib_gb, vshr, zht, shallowc,cu_over_orog,                             &
  bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6,          &
  bl_type_7, z0m_eff_gb,z0h_eff_gb,                                     &
! OUT diagnostic requiring stash flags :
  fme,                                                                  &
! OUT diagnostics required for soil moisture nudging scheme :
  wt_ext,ra,                                                            &
! OUT data required for tracer mixing :
  rho_aresist,aresist,resist_b,ntml,kent, we_lim, t_frac, zrzi,         &
  kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,                     &
!OUT variables required for mineral dust scheme
  r_b_dust,dust_flux,dust_emiss_frac,                                   &
  u_s_t_tile,u_s_t_dry_tile,u_s_std_tile,                               &
! OUT variables required in IMP_SOLVER
  alpha1,ashtf,bq_gb,bt_gb,dtrdz_charney_grid,rdz_charney_grid,         &
  dtrdz_u,dtrdz_v,rdz_u,rdz_v,                                          &
  fraca,rhokh_tile,smc,chr1p5m,resfs,z0hssi,z0mssi, z1_tq, ustargbm,    &
! OUT additional variables for JULES
 le_tile,radnet_sice,radnet_tile,rib_tile,rho_aresist_tile,aresist_tile,&
  resist_b_tile,alpha1_sice,ashtf_tile,fqt_ice,                         &
  ftl_ice,resft,rhokh_sice,z0h_tile,z0m_gb,z0m_tile,chr1p5m_sice,       &
  g_leaf,gpp_ft,npp_ft, resp_p_ft,resp_s,resp_s_tot,resp_w_ft,          &
  gc,canhc_tile,wt_ext_tile,flake,tile_index,tile_pts,tile_frac,fsmc,   &
  rib_ssi,vshr_land,vshr_ssi,tstar_land,tstar_ssi,                      &
  rhokh_mix,dtstar_tile,dtstar,hcons,emis_tile,emis_soil,               &
! OUT data required elsewhere in um system :
  gpp,npp,resp_p,t1_sd,q1_sd,ntdsc, nbdsc, wstar,wthvs,uw0,vw0          &
  )

END IF                    ! on error code = 0

! end of routine NI_bl_ctl
IF (lhook) CALL dr_hook('NI_BL_CTL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE ni_bl_ctl
