! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: boundary layer

SUBROUTINE diagnostics_bl(                                              &
! IN levels / grids / switches
  boundary_layer_levels, land_points, dsm_levels                        &
, dim_cs1, dim_cs2, l_murk                                              &
, rhc_row_length, rhc_rows                                              &
, sin_theta_longitude, cos_theta_longitude                              &
, land_index, ntiles,npft,nice,nice_use                                 &
, l_dust,l_dust_diag, sq_t1p5                                           &
! IN fields for diagnostics
, aerosol, rhcrit, plsp, cca_2d                                         &
, ls_rain, ls_snow, conv_rain, conv_snow                                &
, p_star, rhcpt, ntml, cumulus, rho1, qcf                               &
, t_incr_diagnostic, q_incr_diagnostic                                  &
, qcl_incr_diagnostic, qcf_incr_diagnostic                              &
, cf_incr_diagnostic                                                    &
, cfl_incr_diagnostic, cff_incr_diagnostic                              &
, u_incr_diagnostic,v_incr_diagnostic                                   &
, t1p5m, ml_depth, u10m, v10m, q1p5m                                    &
, e_sea, h_sea, ei                                                      &
, sea_ice_htf, sice_mlt_htf                                             &
, snomlt_surf_htf, bl_top                                               &
, bl_type_1,bl_type_2,bl_type_3,bl_type_4                               &
, bl_type_5,bl_type_6,bl_type_7                                         &
, fqt, ftl, z0m_gb, z0m_eff, z0h_eff                                    &
, rib, surf_latent_heat_flux, taux, tauy                                &
, wind_mixing_energy, t_soil                                            &
, surf_ht_flux, surf_ht_flux_land,surf_ht_flux_sice                     &
, rib_ssi,ftl_ssi,e_ssi,ei_sice, vshr_land,vshr_ssi                     &
, taux_land,taux_ssi,tauy_land,tauy_ssi                                 &
, radnet_sice,flandg, land_sea_mask                                     &
, sil_orog_land,ho2r2_orog,gs,gpp,npp,resp_p                            &
, ecan_tile,esoil_tile,gpp_ft,ftl_tile                                  &
, npp_ft,resp_p_ft,resp_s,resp_s_tot,cs                                 &
, rib_tile,es,ecan,fsmc,radnet_tile                                     &
, tstar_tile,canopy,catch,z0m_tile,g_leaf                               &
, t1p5m_tile,q1p5m_tile,le_tile,ei_tile,olr                             &
, epot_tile,tile_frac, co2_flux_tot, land_co2, dust_flux                &
, dust_emiss_frac,u_s_t_tile,u_s_t_dry_tile                             &
, u_s_std_tile,drydep2, bl_diag                                         &
! variables required for soil moisture nudging scheme macro
, rhokh,resfs,chr1p5m,alpha1,ra,wt_ext                                  &
, lai_ft,canht_ft,gc                                                    &
! MGS extra bl vars for UKCA
, rhokh_mix, rho_aresist, aresist, resist_b, r_b_dust                   &
, dtrdz_charney_grid, kent, we_lim, t_frac, zrzi, kent_dsc              &
, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,                               &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
! INOUT stash workspace
  stashwork)

! Purpose:
!  Calculates diagnostics generated from boundary layer routines
!  (UM section 3).

! Method:
! Required level lists and logical switches are determined by the
! calling routine from STASH requests and STASHflags.
! Intercepted arrays - calculated unconditionally - and diagnostic
! arrays - dependent on STASHflags - are input from the boundary
! layer routines called previously. Each diagnostic is simply
! copied into the STASHwork array to be passed on to STASH for
! output processing.

!  Diagnostics currently available

! STASH item (all section 3 )
! ----------------------------------------------------------------------
! 181 temperature increment across bl + PC2 routines (model levels)
! 182 humidity    increment across bl + PC2 routines (model levels)
! 183 liq con qCL increment across bl + PC2 routines (model levels)
! 184 ice con qCF increment across bl routines (model levels)
! 185 u wind      increment across bl routines (model levels)
! 186 v wind      increment across bl routines (model levels)
! 192 total  cloud fraction increment across bl routines (model levels)
! 193 liquid cloud fraction increment across bl routines (model levels)
! 194 ice    cloud fraction increment across bl routines (model levels)
! 189 liq temp Tl increment across bl routines (model levels)
! 190 tot hum  qT increment across bl routines (model levels)
! ------------------------------------------------------------
! 209 10m u wind (native grid = 'C')
! 225 10m u wind 'B' grid
!    Note: simple linear horizontal interpolation from 'C' to 'B' grid
! 210 10m v wind (native grid = 'C')
! 226 10m v wind 'B' grid
!    Note: simple linear horizontal interpolation from 'C' to 'B' grid
! 227 10m wind speed 'B' grid
! 254 1.5m liquid water temperature = tl calculated in imp_solver
! 255 1.5m total water = qt  (kg water/ kg air) calculated in imp_solver
! 245 1.5m relative humidity [%]
!    Note: (tl,qt at 1.5m) converted to (t,q at 1.5m) using cloud
!          scheme, assuming level 1 rhcrit. No negative rh allowed.
! 236 1.5m temperature
! 237 1.5m specific humidity
! 248 1.5m fog fraction
! 250 1.5m dewpoint
! 251 Silhouette orographic roughness field (A/S)
! 252 Half of peak to trough height of sub-grid orography
! 253 1.5m mist fraction
! 216 heat fluxes on model (rho_) levels
!    Note: o/p 1->boundary_levels but level 1 is not meaningful - set
!          same values as surface field 217.
! 217 surface heat flux
! 234 surface latent heat flux
! 219 surface and model level u wind stress
!    Note: not enabled at vn5.1 since o/p should be on (theta_) levels
!          0->boundary_levels-1, but would be labelled as
!          1->boundary_levels    under current system.
! 220 surface and model level v wind stress
!    Note: as 219.
! 221 magnitude of wind stress on B-grid
! 238 soil temperature on soil levels
!    Note: mdi=-1.e30 for sea points
! 025 mixed layer depth
! 224 wind mixing energy
! 305-310 boundary layer types
! 222 total moisture flux profile on model levels
! 223 total surface moisture flux
! 241 total surface moisture flux per timestep
! 026 effective roughness length for momentum
! 208 RIB - measure of surface stability
! 304 turbulent cloud depth
! 247 visibility
! 232 surface evaporation weighted by leads
! 228 sensible heat flux over open sea
! 229 soil evaporation
! 231 sublimation of sea-ice (accumulation over timestep)
! 201 melting of bottom of sea-ice (GBM)
! 235 melting of top of sea-ice (GBM)
! 256 heat flux through sea ice on categories
! 257 heat flux melting surface  sea ice on categories
! 258 heat flux due to melting of snow
! 202 heat flux from surface to deep soil level 1
! 296 evap. from soil surface
! 297 evap. from canopy
! 298 sublim. from surface
! 259 canopy conductance
! 261 gross primary productivity
! 262 net primary productivity
! 263 plant respiration
! 293 soil respiration
! 287 canopy evap. on tiles
! 288 transpiration + soil evap. on tiles
! 290 sensible heat flux on tiles
! 294 bulk Richardson number on tiles
! 314 net radiation on tiles
! 316 surface temperature on tiles
! 317 tile fractions
! 318 Leaf area indices on vegetated tiles
! 319 Canopy height on vegetated tiles
! 321 canopy water on tiles
! 322 canopy capacity of tiles
! 324 snow adjusted roughness length of tiles
! 328 1.5m temperature over tiles
! 329 1.5m specific humidity over tiles
! 330 latent heat flux on tiles
! 331 sublimation on tiles
! 341 1.5m land mean temperature over tiles
! 342 1.5m land mean specific humidity over tiles
! 50 exchange coefficient for moisture
! 51 combined soil/stomatol/aerodynamic resistance to evaporation
!    on tiles
! 52 ratio of coefficients required for 1.5m T calculation
!    on tiles
! 53 grad of sat humidity w.r. to temp between surface and level 1
!    on tiles
! 54 aerodynamic resistance (s/m)
! 55 cumulative evaporation from soil due to transpiration
! 56 combined soil/stomatol/aerodynamic resistance to evaporation
!    aggregated over surface tiles 
! 57 ratio of coefficients required for 1.5m T calculation
!    aggregated over surface tiles 
! 58 grad of sat humidity w.r. to temp between surface and level 1
!    aggregated over surface tiles 
! 289 gross primary productivity on PFTs
! 291 net primary productivity on PFTs
! 292 plant respiration on PFTs
! 313 soil moisture availability on PFTs
! 325 leaf turnover rate of PFTs
! 332 TOA outgoing LW radiation 
! 334 land potential evaporation rate
! 335 potential evaporation rates on land tiles
! 356 height of surface mixed layer top
! 357 height of dec stratocumulus layer top
! 358 boundary layer depth diagnosed from critical Ri
! 359 height of diagnostic parcel top
! 462 Stomatal conductance
! 464 Obukhov length
! 465 Friction velocity
! 467 Surface buoyancy flux
! 468 Gradient Richardson number
! 466 Convective velocity scale
! 469 Vertical buoyancy gradient
! 470 Modulus of wind shear
! 471 BL Momentum diffusion
! 472 BL heat diffusion
! 473 Turbulent kinetic energy
! 474 x component of orographic stress
! 475 y component of orographic stress
! 476 Combined boundary layer type diagnostic
! 501 Mixing length for momentum
! 502 Mixing length for heat and moisture
! 503 Km diffusion coeff from local scheme
! 504 Kh diffusion coeff from local scheme
! 505 Km diffusion coeff for surface-driven turb
! 506 Kh diffusion coeff for surface-driven turb
! 507 Km diffusion coeff for cloud-top-driven turb
! 508 Kh diffusion coeff for cloud-top-driven turb
! ------------------------------------------------------------
! MGS Extra BL variables needed in STASH for UKCA.
! N.B. ZH is already passed in and renamed 'ml_depth'
!  60 rhokh_mix
!  61 RHO_ARESIST (RHOSTAR*CD_STD*VSHR)
!  62 ARESIST [ 1/(CD_STD*VSHR) ]
!  63 RESIST_B (1/CH-1/CD_STD)/VSHR
!  64 DTRDZ_CHARNEY_GRID
!  65 GRID-LEVEL OF SML INVERSION (kent)
!  66 Rho * entrainment rate
!  67 Fraction of the timestep
!  68 zrzi
!  69 GRID-LEVEL OF DSC INVERSION (kent_dsc)
!  70 Rho * entrainment rate (dsc)
!  71 Fraction of the timestep
!  72 zrzi
!  73 ZHSC  Top of decoupled layer
!  74 Surface layer resist for dust div1
!  75 Surface layer resist for dust div2
!  76 Surface layer resist for dust div3
!  77 Surface layer resist for dust div4
!  78 Surface layer resist for dust div5
!  79 Surface layer resist for dust div6
! ------------------------------------------------------------
! TKE based turbulent closure model
! 130 counter gradient term of taux
! 131 counter gradient term of tauy
! 132 counter gradient term of ftl
! 133 counter gradient term of fqw
! 134 mixing length
! 135 production rate of TKE by shear
! 136 production rate of TKE by buoyancy
! 137 dissipation rate of TKE
! 138 non-dimensional diffusion coefficient for u and v (SM)
! 139 non-dimensional diffusion coefficient for t and q (SH)
! 140 non-gradient buoyancy flux
! 141 cloud fraction used in the TKE schemes
! 142 condensed water used in the TKE schemes
! 143 standard deviation of the distribution function
!     in the TKE schemes
! ------------------------------------------------------------

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

USE dynamics_grid_mod, ONLY: l_vatpoles

USE atm_fields_bounds_mod, ONLY:                                        &
    udims, udims_s, vdims, vdims_s, tdims, tdims_s, qdims, pdims,pdims_s
USE atmos_constants_mod, ONLY: vkman, cp
USE mask_compression, ONLY: expand_from_mask
USE bl_diags_mod, ONLY : strnewbldiag
USE swapable_field_mod, ONLY : swapable_field_pointer_type
USE visbty_constants_mod, ONLY:                                         &
     n_vis_thresh, vis_thresh, calc_prob_of_vis
USE water_constants_mod, ONLY: lc
USE conversions_mod, ONLY: pi
USE dust_parameters_mod, ONLY: ndiv, ndivh
USE switches, ONLY: l_aggregate
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY : ereport
USE proc_info_mod, ONLY: global_row_length, global_rows, at_extremity,  &
      model_domain, gc_proc_row_group
USE field_types, ONLY: fld_type_u, fld_type_v
USE domain_params, ONLY: mt_global
USE timestep_mod, ONLY: timestep
USE missing_data_mod, ONLY: rmdi
USE beta_precip_mod, ONLY: beta_precip
USE u_to_p_mod, ONLY: u_to_p
USE uc_to_ub_mod, ONLY: uc_to_ub
USE v_to_p_mod, ONLY: v_to_p
USE vc_to_vb_mod, ONLY: vc_to_vb
USE Submodel_Mod
use cable_data_mod!, only :cable

IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

LOGICAL, INTENT(IN) ::                                                  &
  sq_t1p5                                                               &
, l_murk

INTEGER, INTENT(IN) ::                                                  &
  boundary_layer_levels                                                 &
                        ! number of boundary layer levels
, land_points                                                           &
, dsm_levels                                                            &
, ntiles                                                                &
                   ! number of land surface tiles
, npft                                                                  &
                   ! number of plant funcional types
, nice                                                                  & 
                   ! number of seaice categories 
, nice_use         ! number of seaice categories fully used in  
                   !  surface exchange 

LOGICAL, INTENT(IN) ::                                                  &
  land_sea_mask(pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end)! land sea mask

INTEGER, INTENT(IN) ::                                                  &
  ntml (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                        ! Height of diagnosed BL top

LOGICAL, INTENT(IN) ::                                                  &
 cumulus (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                        ! Logical indicator of convection

INTEGER, INTENT(IN) ::                                                  &
  dim_cs1, dim_cs2        ! soil carbon dimensions

! Primary Arrays used in all models
REAL, INTENT(IN) ::                                                     &
  p_star(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

REAL, INTENT(IN) ::                                                     &
  aerosol     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end)                               &
, rhcrit(qdims%k_end)                                                   &
, rho1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)             &
                                 ! Air density at level 1/ kg m-3
, qcf(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,              &
      qdims%k_end)                                                      &
, t_incr_diagnostic(tdims%i_start:tdims%i_end,                          &
                    tdims%j_start:tdims%j_end,tdims%k_end)              &
, q_incr_diagnostic(qdims%i_start:qdims%i_end,                          &
                    qdims%j_start:qdims%j_end,qdims%k_end)              &
, qcl_incr_diagnostic(qdims%i_start:qdims%i_end,                        &
                     qdims%j_start:qdims%j_end,qdims%k_end)             &
, qcf_incr_diagnostic(qdims%i_start:qdims%i_end,                        &
                     qdims%j_start:qdims%j_end,qdims%k_end)             &
, cf_incr_diagnostic(qdims%i_start:qdims%i_end,                         &
                     qdims%j_start:qdims%j_end,qdims%k_end)             &
, cfl_incr_diagnostic(qdims%i_start:qdims%i_end,                        &
                      qdims%j_start:qdims%j_end,qdims%k_end)            &
, cff_incr_diagnostic(qdims%i_start:qdims%i_end,                        &
                      qdims%j_start:qdims%j_end,qdims%k_end)            &
! _diagnostic contain u,v increments before imp_solver
, u_incr_diagnostic(udims%i_start:udims%i_end,                          &
    udims%j_start:udims%j_end,udims%k_start:udims%k_end)                &
, v_incr_diagnostic(vdims%i_start:vdims%i_end,                          &
    vdims%j_start:vdims%j_end,vdims%k_start:vdims%k_end)                &
, t1p5m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)            &
                           ! IN tl at 1.5m but also
                           ! OUT workspace as T at 1.5m
, u10m(udims%i_start:udims%i_end,udims%j_start:udims%j_end)             &
, v10m(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)             &
, q1p5m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)            &
                           ! IN qt at 1.5m (kg water/ kg air)
                           ! OUT workspace as q at 1.5m
!                                ! (kg water/ kg air)
, ftl_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)          &
, surf_latent_heat_flux(pdims%i_start:pdims%i_end,                      &
                        pdims%j_start:pdims%j_end)                      &
, ml_depth(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)         &
, taux(udims%i_start:udims%i_end,udims%j_start:udims%j_end,             &
                       boundary_layer_levels)                           &
, tauy(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,             &
                         boundary_layer_levels)                         &
, taux_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end)        &
, taux_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end)         &
, tauy_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)        &
, tauy_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)         &
, radnet_sice(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,nice_use)                       &
, flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)           &
, wind_mixing_energy(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end)                         &
, t_soil(land_points,dsm_levels)                                        &
, surf_ht_flux(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end)                               &
, surf_ht_flux_land(pdims%i_start:pdims%i_end,                          &
                    pdims%j_start:pdims%j_end)                          &
, surf_ht_flux_sice(pdims%i_start:pdims%i_end,                          &
                    pdims%j_start:pdims%j_end,nice_use)                 &
!*APL*DIAGS
, bl_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)           &
, bl_type_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)        &
                            ! IN Indicator set to 1.0 if stable
!                                 !     b.l. diagnosed, 0.0 otherwise.
, bl_type_2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)        &
                            ! IN Indicator set to 1.0 if Sc over
!                                 !     stable surface layer diagnosed,
!                                 !     0.0 otherwise.
, bl_type_3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)        &
                            ! IN Indicator set to 1.0 if well
!                                 !     mixed b.l. diagnosed,
!                                 !     0.0 otherwise.
, bl_type_4(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)        &
                            ! IN Indicator set to 1.0 if
!                                 !     decoupled Sc layer (not over
!                                 !     cumulus) diagnosed,
!                                 !     0.0 otherwise.
, bl_type_5(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)        &
                            ! IN Indicator set to 1.0 if
!                                 !     decoupled Sc layer over cumulus
!                                 !     diagnosed, 0.0 otherwise.
, bl_type_6(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)        &
                            ! IN Indicator set to 1.0 if a
!                                 !     cumulus capped b.l. diagnosed,
!                                 !     0.0 otherwise.
, bl_type_7(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)        &
                            ! IN Indicator set to 1.0 if a
!                                 !     shear-dominated b.l. diagnosed,
!                                 !     0.0 otherwise.
, ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
                      boundary_layer_levels)                            &
, fqt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
                      boundary_layer_levels)                            &
, z0m_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)           &
, z0m_eff(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)          &
, z0h_eff(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)          &
, rib(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)              &
, rib_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)          &
, e_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)            &
, e_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)            &
, h_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)            &
, ei(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)               &
, ei_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use) &
, sea_ice_htf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              nice)                                                     &
, sice_mlt_htf(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,nice)                          &
, snomlt_surf_htf(pdims%i_start:pdims%i_end,                            &
                  pdims%j_start:pdims%j_end)                            &
, plsp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)             &
, cca_2d(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)           &
, ls_rain(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)          &
, ls_snow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)          &
, conv_rain(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)        &
, conv_snow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)        &
, vshr_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)        &
, vshr_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! MGS new bl vars for use by UKCA, intent(in):
REAL, INTENT(IN) ::                                                     &
  rhokh_mix (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             boundary_layer_levels)                                     &
, rho_aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
, aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)          &
, resist_b(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)         &
, r_b_dust(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           ndiv)                                                        &
,  dtrdz_charney_grid(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end,                        &
                      boundary_layer_levels)                            &
, we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)         &
, t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)         &
, zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)           &
,we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)      &
,t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)      &
, zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)       &
, zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
INTEGER, INTENT(IN) ::                                                  &
  kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
  kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

!     Declaration of new BL diagnostics.
TYPE (strnewbldiag) :: bl_diag

REAL, INTENT(IN) ::                                                     &
  sil_orog_land(land_points)                                            &
                                 ! IN silhouette orog roughness
, ho2r2_orog(land_points)                                               &
                                 ! IN 1/2 pk to trough ht of orog
, gs(land_points)                                                       &
                                 ! IN Canopy conductance
, gpp(land_points)                                                      &
                                 ! IN Gross primary productivity
, npp(land_points)                                                      &
                                 ! IN Net primary productivity
, resp_p(land_points)                                                   &
                                 ! IN Plant respiration
, ecan_tile(land_points,ntiles)                                         &
                                 ! IN Canopy evaporation on tiles
, esoil_tile(land_points,ntiles)                                        &
                                 ! IN Soil evap. on tiles
!CABLE: kdcorbin, 11/10 - changed from npft
, gpp_ft(land_points,ntiles)                                              &
                                 ! IN GPP on PFTs
, ftl_tile(land_points,ntiles)                                          &
                                 ! IN Sensible heat flux on tiles
, npp_ft(land_points,ntiles)                                              &
                                 ! IN NPP on PFTs
, resp_p_ft(land_points,ntiles)                                           &
                                 ! IN RESP_P on PFTs
, resp_s(land_points,dim_cs1)                                           &
                                 ! IN Soil respiration
, resp_s_tot(dim_cs2)                                                   &
                                 ! IN Soil respiration
, cs(land_points,dim_cs1)                                               &
                                 ! IN Soil carbon
, rib_tile(land_points,ntiles)                                          &
                                 ! IN RIB on tiles
, es(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)               &
                                 ! IN Evap from soil surface
, ecan(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)             &
                                 ! IN Evap from canopy
, fsmc(land_points,npft)                                                &
                                 ! IN Soil moisture availability
!                                      !    factor on PFTs
, radnet_tile(land_points,ntiles)                                       &
                                 ! IN Net radiation on tiles
, tstar_tile(land_points,ntiles)                                        &
                                 ! IN Tile surface temperatures
, lai_ft(land_points,npft)                                              &
                                 ! IN LAI of PFTs
, canht_ft(land_points,npft)                                            &
                                 ! IN Canopy height of PFTs
, gc(land_points,ntiles)                                                &
                                 ! IN Stomatal conductance on tile
, canopy(land_points,ntiles)                                            &
                                 ! IN Canopy water on tiles
, catch(land_points,ntiles)                                             &
                                 ! IN Canopy capacity on tiles
, z0m_tile(land_points,ntiles)                                          &
                                 ! IN Roughness length on tiles
, g_leaf(land_points,ntiles)                                              &
                                  ! IN Leaf turnover rate for PFTs
, t1p5m_tile(land_points,ntiles)                                        &
                                 ! IN 1.5m temperature over tiles
, q1p5m_tile(land_points,ntiles)                                        &
                                 ! IN 1.5m specific humidity
!                                      !    over tiles
, le_tile(land_points,ntiles)                                           &
                                 ! IN Surface latent heat flux
!                                      !    on tiles
, ei_tile(land_points,ntiles)                                           &
                                 ! IN Sublimation on tiles
, olr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)              &
                                 ! IN TOA outgoing LW radiation
, epot_tile(land_points,ntiles)                                         &
                                 ! IN Potential evap on tiles
,tile_frac(land_points,ntiles)                                          &
                        ! IN fractional coverage for each
                        !    surface tile
,co2_flux_tot(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
                                 ! IN total CO2 flux
, land_co2(land_points)             ! IN terrestrial CO2 flux

LOGICAL, INTENT(IN) ::                                                  &
  l_dust                                                                &
         !IN switch for prognostic mineral dust
, l_dust_diag
         !IN switch for diagnostic mineral dust lifting

! Variables for STASH macro for soil moisture nudging scheme
REAL, INTENT(IN) ::                                                     &
 rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,             &
       boundary_layer_levels)                                           &
,resfs(land_points, ntiles)                                             &
,chr1p5m(land_points, ntiles)                                           &
,alpha1(land_points, ntiles)                                            &
,ra(land_points)                                                        &
,wt_ext(land_points,dsm_levels)

INTEGER, INTENT(IN) ::                                                  &
 land_index(land_points)

! Variables required for conversion of 1.5m TL and QT
INTEGER, INTENT(IN) ::                                                  &
 rhc_row_length, rhc_rows ! rhcpt dimensions
!  rhcpt is point-varying rhcrit and can be a 1-d parameter or a
!  3-d diagnostic field, but only level 1 is required here
REAL, INTENT(IN) ::                                                     &
 rhcpt(rhc_row_length, rhc_rows, 1)

! mineral dust variables

REAL, INTENT(IN) ::                                                     &
  dust_flux(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            ndiv)                                                       &
                                  !dust emission flux
, dust_emiss_frac(land_points,ntiles)                                   &
                          ! OUT fraction of tile can emit dust
, u_s_t_tile(land_points,ntiles,ndivh)                                  &
                                     !OUT threshold frict. vel
, u_s_t_dry_tile(land_points,ntiles,ndivh)                              &
                                         !OUT dry soil value
, u_s_std_tile(land_points,ntiles)                                      &
                                  !OUT friction velocity
, drydep2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          ndiv)                    !dep from 2nd level

REAL, INTENT(IN) :: cos_theta_longitude(tdims%i_start:tdims%i_end,      &
                                        tdims%j_start:tdims%j_end)
REAL, INTENT(IN) :: sin_theta_longitude(tdims%i_start:tdims%i_end,      &
                                        tdims%j_start:tdims%j_end)

! Diagnostics info
REAL, INTENT(INOUT) ::                                                  &
 stashwork(*)        ! STASH workspace
!
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
!----------------------------------------------------------------------
! Local variables

!  Local parameters and other physical constants------------------------
REAL lcrcp                        ! Derived parameter.
PARAMETER ( lcrcp=lc/cp )         ! Lat ht of condensation/Cp.
LOGICAL, PARAMETER :: l_mixing_ratio=.FALSE. ! Use mixing ratios

! Tempory variable
REAL, ALLOCATABLE :: gc_tmp(:,:)

INTEGER                                                                 &
  i, j, k, l, n                                                         &
                   ! loop indices
,    icode                                                              &
                          ! Return code  =0 Normal exit  >1 Error
 ,item                                                                  &
                  ! STASH item
 ,sect            ! STASH section

INTEGER :: i_field    ! counter for swap_bounds_mv

PARAMETER( sect = 3 ) ! for boundary layer

CHARACTER(LEN=80)  cmessage

CHARACTER(LEN=*) routinename
PARAMETER ( routinename='diagnostics_bl')

INTEGER                                                                 &
  im_index        ! internal model index

     INTEGER :: land_points_dum
!                Dummy number of land points returned from the
!                routine expand_from_mask: the variable land_points
!                is already set on input

INTEGER                                                                 &
  pp_code_bl_types(7)

INTEGER                                                                 &
 pslevel                                                                &
               !  loop counter for pseudolevels
,pslevel_out   !  index for pseudolevels sent to STASH

LOGICAL                                                                 &
 plltile(ntiles)                                                        &
                    ! pseudolevel list for surface types
,pllpft(npft)                                                           &
                    ! pseudolevel list for PFTs
,pllnice(nice)                                                          &
,pllbl(npft)
              ! pseudolevel list for BL diagnostics with 3rd dim=3

REAL                                                                    &
  interp_data(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
, interp_data_3(pdims%i_end*pdims%j_end*pdims%k_end)                    &
                                              ! work array
, interp_data_bl(pdims%i_start:pdims%i_end,                             &
                 pdims%j_start:pdims%j_end,npft)                        &
! work array, for diagnostics with 3rd dim = 3. Arbitrarily use
!  npft, which equals 5.
, qcl1p5m(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)          &
, beta_ls_rain(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
                                 ! Scattering in LS Rain.
, beta_ls_snow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)     &
                                 ! Scattering in LS Snow.
, beta_c_rain(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
                                 ! Scattering in Conv Rain
, beta_c_snow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)      &
                                 ! Scattering in Conv Snow
, vis(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)              &
, vis_no_precip(pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end)                              &
, vis_ls_precip(pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end)                              &
, vis_c_precip(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end)                               &
, vis_pr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         1,n_vis_thresh)                                                &
, vis_threshold(pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,1,n_vis_thresh)               &
                                               !FOG_FR works for n
                                               ! levels, we want 1
, pvis(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,             &
       n_vis_thresh)                                                    &
, u10mb(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end)            &
                             !10m  u-wind     B-grid
, v10mb(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end)            &
                             !10m  v-wind     B-grid
, ws10mb(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end)           &
                             !10m  wind speed B-grid
, ws10m_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)          &
                             !10m  wind speed P-grid
, tauxb(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end,            &
        boundary_layer_levels)                                          &
                                                !x-stress   B-grid
, tauyb(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end,            &
        boundary_layer_levels)                                          &
                                                !y-stress   B-grid
, taub(udims%i_start:udims%i_end,vdims%j_start:vdims%j_end,             &
       boundary_layer_levels)                                           &
                                                !stress mag B-grid
, sea_ice_htf_gb(pdims%i_start:pdims%i_end,                             &
                 pdims%j_start:pdims%j_end)                             &
                                   ! Gridbox mean flux through ice
, sice_mlt_htf_gb(pdims%i_start:pdims%i_end,                            &
                  pdims%j_start:pdims%j_end)                            &
                                   ! Gridbox_mean ice surface melt
, surf_ht_flux_sice_gb(pdims%i_start:pdims%i_end,                       & 
                       pdims%j_start:pdims%j_end)                       & 
                                   ! Gridbox mean sea ice sfc flux 
, ei_sice_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)       & 
                                   ! Gridbox mean sublimation 
, radnet_sice_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)   &
                                   ! Gridbox mean net sea ice rad
, epot_land(land_points)                                                &
                                 ! Land mean Potential evap
, cs_tot(dim_cs2)                                                       &
                                 ! diagnosed total soil carbon
, t1p5m_land(land_points)                                               &
                           ! Land mean tl at 1.5m but also
, q1p5m_land(land_points)                                               &
                           ! Land mean qt at 1.5m
, bl_type_comb(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end)                               & 
                                   ! Combined boundary layer
                                   ! type diagnostic.
, resfs_land(land_points)                                               &
                                   ! Land mean resfs
, chr1p5m_land(land_points)                                             &
                                   ! Land mean chr1p5m
, alpha1_land(land_points)
                                   ! Land mean alpha1

TYPE(swapable_field_pointer_type) :: fields_to_swap(4)   ! for multivar
                                                         ! swap_bounds

! Variables required for calculating the maximum wind gust
! or wind speed diagnostics
REAL, ALLOCATABLE :: gust_wind(:,:)
REAL, ALLOCATABLE :: std_dev(:,:)
REAL, ALLOCATABLE, TARGET :: u10m_halo(:,:)
REAL, ALLOCATABLE, TARGET :: v10m_halo(:,:)
REAL, ALLOCATABLE :: u10m_p(:,:)
REAL, ALLOCATABLE :: v10m_p(:,:)
REAL, ALLOCATABLE, TARGET :: taux_halo(:,:)
REAL, ALLOCATABLE, TARGET :: tauy_halo(:,:)
REAL, ALLOCATABLE :: tmp_p(:,:)
REAL, ALLOCATABLE :: ustar_imp(:,:)

REAL, ALLOCATABLE :: work3(:,:,:)

REAL, PARAMETER :: tol_obukhov = 1.0e-10
!       Tolerance for testing the Obukhov length: the FORTRAN intrinsic
!       cannot be used here because it is already used as a standard
!       variable in the UM

! The following are declared as arrays rather than scalars to maintain
! consistency of argument type, since they are passed to a general
! purpose subroutine in which they are array arguments
REAL :: mag_vector_np(1)
REAL :: mag_vector_sp(1)
REAL :: dir_vector_np(1)
REAL :: dir_vector_sp(1)

REAL :: SCALE

! Tunable parameters used in the calculation of the maximum
! wind gust. See http://www-nwp/~frpz/gust_diag/max_gust.html
REAL, PARAMETER :: c_ugn=4.0
REAL, PARAMETER :: gust_const=2.29

INTEGER :: i_g,j_g
INTEGER :: y,z
INTEGER, PARAMETER :: pnorth=1
INTEGER, PARAMETER :: peast =2
INTEGER, PARAMETER :: psouth=3
INTEGER, PARAMETER :: pwest =4
INTEGER, PARAMETER :: nodomain = -1

INTEGER :: n_rows

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!CABLE:
! tiled snow diagnostics 
real, dimension(land_points, ntiles, dsm_levels) ::                            &
   ftsoil, fsmcl, fsthf
 
! tiled snow diagnostics 
real, dimension(land_points, ntiles, msn_cable) ::                            &
   fsnow_depth3L, fsnow_mass3L,  &
   fsnow_temp3L, fsnow_rho3L

real, dimension(land_points, ntiles) ::                            &
   fsnow_rho1L, fsnow_age, fsnow_flag3L

IF (lhook) CALL dr_hook('DIAGNOSTICS_BL',zhook_in,zhook_handle)

n_rows=vdims%j_end-vdims%j_start+1

icode = 0 ! Initialise error status
im_index = internal_model_index(atmos_im)

! initialize bl_types pp codes :
pp_code_bl_types(1)      = 305
pp_code_bl_types(2)      = 306
pp_code_bl_types(3)      = 307
pp_code_bl_types(4)      = 308
pp_code_bl_types(5)      = 309
pp_code_bl_types(6)      = 310
pp_code_bl_types(7)      = 340

! Allocate 3d work array if calculating +/- increments for cfl,cff,qcf.
IF (sf(172,sect) .OR. sf(173,sect) .OR.                                 &
    sf(176,sect) .OR. sf(177,sect) .OR.                                 &
    sf(178,sect) .OR. sf(179,sect)) THEN
  ALLOCATE ( work3(qdims%i_start:qdims%i_end,                           &
                   qdims%j_start:qdims%j_end,1:qdims%k_end) )
END IF

! ----------------------------------------------------------------------
! Section 1.  Diagnostic Calculation and output :
! ----------------------------------------------------------------------
! Do any communications here, but only if diagnostic is needed

IF (icode <= 0 .AND. ( sf(230,3) .OR. sf(463,3) )) THEN
! 10m winds needed on p grid
  ALLOCATE(u10m_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
  ALLOCATE(v10m_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
  ALLOCATE(u10m_halo(udims_s%i_start:udims_s%i_end,                     &
                     udims_s%j_start:udims_s%j_end))
  ALLOCATE(v10m_halo(vdims_s%i_start:vdims_s%i_end,                     &
                     vdims_s%j_start:vdims_s%j_end))

  u10m_halo(:,:) = 0.0
  v10m_halo(:,:) = 0.0

  DO y=udims%j_start,udims%j_end
    DO z=udims%i_start,udims%i_end
      u10m_halo(z,y)=u10m(z,y)
    END DO
  END DO
  DO y=vdims%j_start,vdims%j_end
    DO z=vdims%i_start,vdims%i_end
      v10m_halo(z,y)=v10m(z,y)
    END DO
  END DO

! Update halos for u10m_halo and v10m_halo
  i_field = 0
  i_field = i_field + 1
  fields_to_swap(i_field) % field_2d    => u10m_halo(:,:)
  fields_to_swap(i_field) % field_type  =  fld_type_u
  fields_to_swap(i_field) % levels      =  1
  fields_to_swap(i_field) % rows        =  udims%j_end
  fields_to_swap(i_field) % vector      =  .TRUE.

  i_field = i_field + 1
  fields_to_swap(i_field) % field_2d    => v10m_halo(:,:)
  fields_to_swap(i_field) % field_type  =  fld_type_v
  fields_to_swap(i_field) % levels      =  1
  fields_to_swap(i_field) % rows        =  n_rows
  fields_to_swap(i_field) % vector      =  .TRUE.

  IF (sf(463,3)) THEN
!       x-stress:
    ALLOCATE(taux_halo(udims_s%i_start:udims_s%i_end,                   &
                       udims_s%j_start:udims_s%j_end))
    taux_halo(:,:)  = 0.0
    taux_halo(udims%i_start:udims%i_end,udims%j_start:udims%j_end)      &
    = taux(udims%i_start:udims%i_end,udims%j_start:udims%j_end,1)
!       Set up for message passing
    i_field = i_field + 1
    fields_to_swap(i_field) % field_2d    => taux_halo(:,:)
    fields_to_swap(i_field) % field_type  =  fld_type_u
    fields_to_swap(i_field) % levels      =  1
    fields_to_swap(i_field) % rows        =  udims%j_end
    fields_to_swap(i_field) % vector      =  .TRUE.

!       y-stress:
    ALLOCATE(tauy_halo(vdims_s%i_start:vdims_s%i_end,                   &
                       vdims_s%j_start:vdims_s%j_end))
    tauy_halo(:,:) = 0.0
    tauy_halo(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)      &
    = tauy(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,1)
    i_field = i_field + 1
    fields_to_swap(i_field) % field_2d    => tauy_halo(:,:)
    fields_to_swap(i_field) % field_type  =  fld_type_v
    fields_to_swap(i_field) % levels      =  1
    fields_to_swap(i_field) % rows        =  n_rows
    fields_to_swap(i_field) % vector      =  .TRUE.

  END IF

! DEPENDS ON: swap_bounds_2d_mv
  CALL swap_bounds_2d_mv( fields_to_swap, i_field,                      &
                        pdims%i_end, pdims_s%halo_i, pdims_s%halo_j)

! interpolate u and v to p grid.
CALL u_to_p(u10m_halo,                                                  &
                  udims_s%i_start,udims_s%i_end,                        &
                  udims_s%j_start,udims_s%j_end,                        &
                  pdims%i_start,pdims%i_end,                            &
                  pdims%j_start,pdims%j_end,                            &
                  1,                                                    &
                  model_domain,at_extremity,u10m_p)


CALL v_to_p(v10m_halo,                                                  &
                  vdims_s%i_start,vdims_s%i_end,                        &
                  vdims_s%j_start,vdims_s%j_end,                        &
                  pdims%i_start,pdims%i_end,                            &
                  pdims%j_start,pdims%j_end,                            &
                  1,                                                    &
                  model_domain,at_extremity,v10m_p)

IF (.NOT. l_vatpoles) THEN
IF (model_domain  ==  mt_global) THEN
! Overwrite values of U_P, V_P at the poles with the magnitude of
! the vector wind.
! DEPENDS ON: polar_vector_wind_n
    CALL polar_vector_wind_n(                                           &
                       v10m_halo,                                       &
                       sin_theta_longitude,                             &
                       cos_theta_longitude,vdims%i_end,                 &
                       n_rows, 1 , mag_vector_np,                       &
                       dir_vector_np, mag_vector_sp,                    &
                       dir_vector_sp,                                   &
                       vdims_s%halo_i,vdims_s%halo_j, global_row_length,&
                       gc_proc_row_group, at_extremity)

    IF (at_extremity(psouth) ) THEN
        DO i=pdims%i_start,pdims%i_end
          v10m_p(i,1) = mag_vector_sp(1)
          u10m_p(i,1) = 0.0
        END DO
    END IF

    IF (at_extremity(pnorth) ) THEN
      DO i=pdims%i_start,pdims%i_end
        v10m_p(i,pdims%j_end) = mag_vector_np(1)
        u10m_p(i,pdims%j_end) = 0.0
      END DO
    END IF

END IF
END IF ! vatpoles

DEALLOCATE(u10m_halo)
DEALLOCATE(v10m_halo)

  IF (sf(463,3)) THEN

    ALLOCATE(tmp_p(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end))
    ALLOCATE(ustar_imp(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end))

CALL u_to_p(taux_halo,                                                  &
                  udims_s%i_start,udims_s%i_end,                        &
                  udims_s%j_start,udims_s%j_end,                        &
                  pdims%i_start,pdims%i_end,                            &
                  pdims%j_start,pdims%j_end,                            &
                  1,                                                    &
                  model_domain,at_extremity,tmp_p)

!       Zero the zonal stress for consistency with the winds
IF (.NOT. l_vatpoles) THEN
  IF (model_domain  ==  mt_global) THEN

    IF (at_extremity(psouth) ) THEN
      tmp_p(pdims%i_start:pdims%i_end,1) = 0.0
    END IF
    IF (at_extremity(pnorth) ) THEN
      tmp_p(pdims%i_start:pdims%i_end,pdims%j_end) = 0.0
    END IF

  END IF
END IF ! vatpoles
!       Begin to form the magnitude of the stress.
    ustar_imp(:,:) = tmp_p(:,:) * tmp_p(:,:)

CALL v_to_p(tauy_halo,                                                  &
                  vdims_s%i_start,vdims_s%i_end,                        &
                  vdims_s%j_start,vdims_s%j_end,                        &
                  pdims%i_start,pdims%i_end,                            &
                  pdims%j_start,pdims%j_end,                            &
                  1,                                                    &
                  model_domain,at_extremity,tmp_p)

IF (.NOT. l_vatpoles) THEN    
  IF (model_domain  ==  mt_global) THEN
!         Set the stress at the poles as for the winds.
! DEPENDS ON: polar_vector_wind_n
      CALL polar_vector_wind_n(                                         &
                         tauy_halo,                                     &
                         sin_theta_longitude,                           &
                         cos_theta_longitude, vdims%i_end,              &
                         n_rows, 1 , mag_vector_np,                     &
                         dir_vector_np, mag_vector_sp,                  &
                         dir_vector_sp,                                 &
                        vdims_s%halo_i,vdims_s%halo_j,global_row_length,&
                         gc_proc_row_group, at_extremity)

      IF (at_extremity(psouth) ) THEN
        tmp_p(pdims%i_start:pdims%i_end,1) = mag_vector_sp(1)
      END IF
      IF (at_extremity(pnorth) ) THEN
        tmp_p(pdims%i_start:pdims%i_end,pdims%j_end) = mag_vector_np(1)
      END IF

  END IF
END IF ! vatpoles
  ustar_imp(:,:) = ustar_imp(:,:) + tmp_p(:,:) * tmp_p(:,:)


  DEALLOCATE(taux_halo)
  DEALLOCATE(tauy_halo)
  DEALLOCATE(tmp_p)

END IF

END IF ! sf(230,3) or sf(463,3)
!
! Interpolation to b-grid also requires swap_bounds
!
IF (icode <= 0 .AND. (sf(225,3) .OR. sf(227,3))) THEN
         ! Horizontal interpolation from 'C' to 'B' grid
         CALL uC_to_uB(u10m,                                            &
              pdims%i_end,udims%j_end,n_rows,1,                         &
              udims_s%halo_i,udims_s%halo_j,                            &
              u10mB)
END IF

      IF (icode <= 0 .and. (sf(226,3) .or. sf(227,3))) THEN
       IF (l_vatpoles) THEN
        IF (model_domain == mt_global .AND. .NOT. sf(225,3) ) THEN
          ! In the vatpoles case, v on B grid 10m is 
          !   dependent on u on B grid 10m
          icode = -1        ! Warning
          Cmessage='v on B grid 10m not possible when '                 &
               //'u on B grid 10m not active. Diagnostic aborted'

   CALL Ereport(Routinename,icode,Cmessage)

        ELSE
          CALL vC_to_vB(v10m, udims%j_end,                               &
               vdims%i_end,n_rows,1,vdims_s%halo_i,vdims_s%halo_j,       &
               global_row_length, v10mB ,u10mB )
        END IF
       ELSE      
          CALL vC_to_vB(v10m, udims%j_end,                               &
               vdims%i_end,n_rows,1,vdims_s%halo_i,vdims_s%halo_j,       &
               global_row_length, v10mB )
       END IF ! vatpoles
      END IF

      IF (icode <= 0 .AND. sf(221,3)) THEN
         ! Horizontal interpolation from 'C' to 'B' grid
         CALL uc_to_ub(taux,                                            &
              pdims%i_end,udims%j_end,n_rows,boundary_layer_levels,     &
              udims_s%halo_i,udims_s%halo_j, tauxb)
         IF (l_vatpoles) THEN
         CALL vC_to_vB(tauy, udims%j_end,                               &
              vdims%i_end,n_rows,boundary_layer_levels,                 &
              udims_s%halo_i,udims_s%halo_j, global_row_length,         &
              tauyB, tauxB )
         ELSE
         CALL vC_to_vB(tauy, udims%j_end,                               &
              vdims%i_end,n_rows,boundary_layer_levels,                 &
              udims_s%halo_i,udims_s%halo_j, global_row_length, tauyB )
         END IF ! vatpoles
      END IF

! End of communication section

! ----------------------------------------------------------------------
!   Copy diagnostic information to STASHwork for STASH processing

! ----------------------------------------------------------------------
! DIAG.03181 Copy T   INC: bdy layer + PC2 condensation to stashwork
! ----------------------------------------------------------------------
item = 181
! Diag03181_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

! Diag03181_do1:
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        l = i+(j-1)*tdims%i_end+(k-1)*tdims%j_end*tdims%i_end
        interp_data_3(l) = t_incr_diagnostic(i, j, k)
      END DO
    END DO
  END DO  ! Diag03181_do1

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
       interp_data_3,                                                   &
       tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="diag_bl  : error in copydiag_3d(T   INC: bdy lay+PC2)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag03181_if1

! ----------------------------------------------------------------------
! DIAG.03182 Copy q   INC: bdy layer + PC2 condensation to stashwork
! ----------------------------------------------------------------------
item = 182
! Diag03182_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

! Diag03182_do1:
  DO k = 1, qdims%k_end
    DO j = qdims%j_start,qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        l = i+(j-1)*qdims%i_end+(k-1)*qdims%j_end*qdims%i_end
        interp_data_3(l)  =   q_incr_diagnostic(i, j, k)
      END DO
    END DO
  END DO  ! Diag03182_do1

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
       interp_data_3,                                                   &
       qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="diag_bl  : error in copydiag_3d(q   INC: bdy lay+PC2)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag03182_if1

! ----------------------------------------------------------------------
! DIAG.03183 Copy qCL INC: bdy layer + PC2 condensation to stashwork
! ----------------------------------------------------------------------
item = 183
! Diag03183_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

! Diag03183_do1:
  DO k = 1, qdims%k_end
    DO j = qdims%j_start,qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        l = i+(j-1)*qdims%i_end+(k-1)*qdims%j_end*qdims%i_end
        interp_data_3(l)  = qcl_incr_diagnostic(i, j, k)
      END DO
    END DO
  END DO  ! Diag03183_do1

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
       interp_data_3,                                                   &
       qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="diag_bl  : error in copydiag_3d(qCL INC: bdy lay+PC2)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag03183_if1

! ----------------------------------------------------------------------
! DIAG.03170 Positive part of 3183
! ----------------------------------------------------------------------
item = 170
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  DO k = 1, qdims%k_end
    DO j = qdims%j_start,qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        l = i+(j-1)*qdims%i_end+(k-1)*qdims%j_end*qdims%i_end
        interp_data_3(l)  = MAX(0.0,qcl_incr_diagnostic(i, j, k))
      END DO
    END DO
  END DO

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
       interp_data_3,                                                   &
       qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="diag_bl  : error in copydiag_3d(qCL INC: bdy lay+PC2)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! DIAG.03171 Negative part of 3183
! ----------------------------------------------------------------------
item = 171
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  DO k = 1, qdims%k_end
    DO j = qdims%j_start,qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        l = i+(j-1)*qdims%i_end+(k-1)*qdims%j_end*qdims%i_end
        interp_data_3(l)  = MIN(0.0,qcl_incr_diagnostic(i, j, k))
      END DO
    END DO
  END DO

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
       interp_data_3,                                                   &
       qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="diag_bl  : error in copydiag_3d(qCL INC: bdy lay+PC2)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! DIAG.03184 Copy qCF INC: bdy layer to stashwork
! ----------------------------------------------------------------------
item = 184
! Diag03184_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN
! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
       qcf_incr_diagnostic,                                             &
       qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="diag_bl  : error in copydiag_3d(Qcf INC: bdy layer)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag03184_if1

! ----------------------------------------------------------------------
! DIAG.03172 Positive part of 3184
! ----------------------------------------------------------------------
item = 172
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  DO k = 1, qdims%k_end
    DO j = qdims%j_start,qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        work3(i,j,k)=MAX(0.0,qcf_incr_diagnostic(i,j,k))
      END DO
    END DO
  END DO

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
       work3,                                                           &
       qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="diag_bl  : error in copydiag_3d(Qcf INC: bdy layer)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! DIAG.03173 Negative part of 3184
! ----------------------------------------------------------------------
item = 173
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

  DO k = 1, qdims%k_end
    DO j = qdims%j_start,qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        work3(i,j,k)=MIN(0.0,qcf_incr_diagnostic(i,j,k))
      END DO
    END DO
  END DO

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
       work3,                                                           &
       qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="diag_bl  : error in copydiag_3d(Qcf INC: bdy layer)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF

! ----------------------------------------------------------------------
! DIAG.03185 Copy U Wind INC: bdy layer to stashwork
! ----------------------------------------------------------------------
item = 185  ! u wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork(si(item,sect,im_index)),                 &
        u_incr_diagnostic,                                              &
        pdims%i_end,udims%j_end,udims%k_end,0,0,0,0, at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
     cmessage=": error in copydiag_3d(item 185)"//cmessage
     CALL ereport(routinename,icode,cmessage)
   END IF

END IF  !  sf(item,sect)

item = 186  ! v wind increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork(si(item,sect,im_index)),                 &
        v_incr_diagnostic,                                              &
        vdims%i_end,n_rows,vdims%k_end,0,0,0,0, at_extremity,           &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
     cmessage=": error in copydiag_3d(item 186)"//cmessage
     CALL ereport(routinename,icode,cmessage)
   END IF

END IF  !  sf(item,sect)

item = 192  ! total cloud fraction increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork(si(item,sect,im_index)),                 &
        cf_incr_diagnostic,                                             &
        qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
     cmessage=": error in copydiag_3d(item 192)"//cmessage
     CALL ereport(routinename,icode,cmessage)
   END IF

END IF  !  sf(item,sect)

item = 193  ! liquid cloud fraction increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork(si(item,sect,im_index)),                 &
        cfl_incr_diagnostic,                                            &
        qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
     cmessage=": error in copydiag_3d(item 193)"//cmessage
     CALL ereport(routinename,icode,cmessage)
   END IF

END IF  !  sf(item,sect)

item = 176  ! liquid cloud fraction increment: positive
IF (icode <= 0 .AND. sf(item,sect)) THEN

  DO k = 1, qdims%k_end
    DO j = qdims%j_start,qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        work3(i,j,k)=MAX(0.0,cfl_incr_diagnostic(i,j,k))
      END DO
    END DO
  END DO

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        work3,                                                          &
        qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 176)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF !  sf(item,sect)

item = 177  ! liquid cloud fraction increment: negative
IF (icode <= 0 .AND. sf(item,sect)) THEN

  DO k = 1, qdims%k_end
    DO j = qdims%j_start,qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        work3(i,j,k)=MIN(0.0,cfl_incr_diagnostic(i,j,k))
      END DO
    END DO
  END DO

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        work3,                                                          &
        qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 177)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF !  sf(item,sect)

item = 194  ! ice cloud fraction increment
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d (stashwork(si(item,sect,im_index)),                 &
        cff_incr_diagnostic,                                            &
        qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode >  0) THEN
     cmessage=": error in copydiag_3d(item 194)"//cmessage
     CALL ereport(routinename,icode,cmessage)
   END IF

END IF  !  sf(item,sect)

item = 178  ! ice cloud fraction increment: positive
IF (icode <= 0 .AND. sf(item,sect)) THEN

  DO k = 1, qdims%k_end
    DO j = qdims%j_start,qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        work3(i,j,k)=MAX(0.0,cff_incr_diagnostic(i,j,k))
      END DO
    END DO
  END DO

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        work3,                                                          &
        qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 178)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF !  sf(item,sect)

item = 179  ! ice cloud fraction increment: negative
IF (icode <= 0 .AND. sf(item,sect)) THEN

  DO k = 1, qdims%k_end
    DO j = qdims%j_start,qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        work3(i,j,k)=MIN(0.0,cff_incr_diagnostic(i,j,k))
      END DO
    END DO
  END DO

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        work3,                                                          &
        qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,      &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
    cmessage=": error in copydiag_3d(item 179)"//cmessage
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF !  sf(item,sect)

! ----------------------------------------------------------------------
! DIAG.03188 TL increment from frictional dissipation
! ----------------------------------------------------------------------
item = 188

IF (icode  <=  0  .AND.  sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
       bl_diag%dtfric,                                                  &
       pdims%i_end,pdims%j_end,                                         &
       boundary_layer_levels,0,0,0,0,at_extremity,                      &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="diag_bl  : error in copydiag_3d(Tl fric INC: bdy lay)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! sf(item,sect)
! ----------------------------------------------------------------------
! DIAG.03189 Copy Tl(Liq) INC: bdy layer to stashwork
! ----------------------------------------------------------------------
item = 189
! Diag03189_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

! Diag03189_do1:
  DO k = 1, qdims%k_end
    DO j = qdims%j_start,qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        l = i+(j-1)*qdims%i_end+(k-1)*qdims%j_end*qdims%i_end
        interp_data_3(l) = t_incr_diagnostic(i, j, k)                   &
                         - lcrcp * qcl_incr_diagnostic(i, j, k)
      END DO
    END DO
  END DO  ! Diag03189_do1

! Diag03189_if2:
  IF (tdims%k_end   >   qdims%k_end)  THEN

! Diag03189_do2:
    DO k = (qdims%k_end + 1), tdims%k_end
      DO j = qdims%j_start,qdims%j_end
        DO i = qdims%i_start, qdims%i_end
          l = i+(j-1)*qdims%i_end+(k-1)*qdims%j_end*qdims%i_end
          interp_data_3(l) = t_incr_diagnostic(i, j, k)
        END DO
      END DO
    END DO  ! Diag03189_do2

  END IF  ! Diag03189_if2

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
       interp_data_3,                                                   &
       tdims%i_end,tdims%j_end,tdims%k_end,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="diag_bl  : error in copydiag_3d(Tl(Liq) INC: bdy lay)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag03189_if1

! ----------------------------------------------------------------------
! DIAG.03190 Copy qT(Liq) INC: bdy layer to stashwork
! ----------------------------------------------------------------------
item = 190
! Diag03190_if1:
IF (icode  <=  0  .AND.  sf(item,sect)) THEN

! Diag03190_do1:
  DO k = 1, qdims%k_end
    DO j = qdims%j_start,qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        l = i+(j-1)*qdims%i_end+(k-1)*qdims%j_end*qdims%i_end
        interp_data_3(l)  =   q_incr_diagnostic(i, j, k)                &
                          +   qcl_incr_diagnostic(i, j, k)
      END DO
    END DO
  END DO  ! Diag03190_do1

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
       interp_data_3,                                                   &
       qdims%i_end,qdims%j_end,qdims%k_end,0,0,0,0, at_extremity,       &
       stlist(1,stindex(1,item,sect,im_index)),len_stlist,              &
       stash_levels,num_stash_levels+1,                                 &
       atmos_im,sect,item,                                              &
       icode,cmessage)

  IF (icode >  0) THEN
    cmessage="diag_bl  : error in copydiag_3d(QT(Liq) INC: bdy lay)"
    CALL ereport(routinename,icode,cmessage)
  END IF
END IF  ! Diag03190_if1

! ----------------------------------------------------------------------
!  10m x Wind
! ----------------------------------------------------------------------
! Item 225 u10m

IF (icode <= 0 .AND. sf(225,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(225,3,im_index)), u10mb,                  &
        pdims%i_end,n_rows,0,0,0,0, at_extremity,                       &
        atmos_im,3,225,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(u10mb) "
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

IF (sf(209,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(209,3,im_index)),u10m,                    &
        pdims%i_end,udims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,209,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(u10m) "
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF ! sf(209,3)

! ----------------------------------------------------------------------
!  10m y Wind
! ----------------------------------------------------------------------
! Item 226 v10m

IF (icode <= 0 .AND. sf(226,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(226,3,im_index)), v10mb,                  &
        pdims%i_end,n_rows,0,0,0,0, at_extremity,                       &
        atmos_im,3,226,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(v10mb) "
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

IF (sf(210,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(210,3,im_index)),v10m,                    &
        vdims%i_end,n_rows,0,0,0,0, at_extremity,                       &
        atmos_im,3,210,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(v10m) "
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF ! sf(210,3)

! ----------------------------------------------------------------------
!  10m Wind Speed on B-grid
! ----------------------------------------------------------------------
! Item 227

IF (icode <= 0 .AND. sf(227,3)) THEN

   ! Calculate wind speed
   DO j = vdims%j_start,vdims%j_end
     DO i = udims%i_start,udims%i_end
       ws10mb(i,j) = SQRT(u10mb(i,j)*u10mb(i,j)                         &
                         +v10mb(i,j)*v10mb(i,j))
     END DO !i
   END DO !j

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(227,3,im_index)), ws10mb,                 &
        pdims%i_end,n_rows,0,0,0,0, at_extremity,                       &
        atmos_im,3,227,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(10m wnd spd B) "
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF ! sf(227,3)

! ----------------------------------------------------------------------
!  10m Wind Speed on C-grid P points for ocean
! ----------------------------------------------------------------------
! Item 230

IF (icode <= 0 .AND. sf(230,3)) THEN

! Calculate 10m wind speed at p C-grid points
  DO j = pdims%j_start,pdims%j_end
    DO i = pdims%i_start,pdims%i_end
      ws10m_p(i,j) = SQRT(u10m_p(i,j)*u10m_p(i,j)                       &
                         +v10m_p(i,j)*v10m_p(i,j))
    END DO !i
  END DO !j
! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(230,3,im_index)),ws10m_p,                 &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,230,                                                 &
        icode,cmessage)
   IF (icode > 0) THEN
      cmessage=": error in copydiag(ws10m_p C-grid  ) "
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF ! sf(230,3)

! ----------------------------------------------------------------------
!   1.5m liquid temperature
! ----------------------------------------------------------------------

! 1.5m liquid temperature
item =   254  ! tl at 1.5m
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,sect,im_index)),                     &
        t1p5m,                                                          &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 254)"//cmessage
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!   1.5m  total water
! ----------------------------------------------------------------------

! 1.5m total water
item =   255  ! qt at 1.5m
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,sect,im_index)),                     &
        q1p5m,                                                          &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 255)"//cmessage
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!  conversion of 1.5m liquid temperature, total water
! ----------------------------------------------------------------------
IF (icode <= 0 .AND. sq_t1p5) THEN

! Use the diagnostic cloud scheme
! Convert 1.5m TL and QT to T and q, assuming p_star and rhcrit at
!  model level 1 approximates to 1.5m.

! DEPENDS ON: ls_cld
  CALL ls_cld(                                                          &
   p_star, rhcpt,  1, boundary_layer_levels                             &
  ,rhc_row_length, rhc_rows, ntml, cumulus                              &
                                           ! in
  ,l_mixing_ratio                                                       &
                                           ! in
  ,t1p5m                                                                &
                                           ! in/out
  ,interp_data_3                                                        &
                                           ! work_cf
  ,q1p5m                                                                &
                                           ! in/out
  ,qcf                                                                  &
                                           ! in
  ,qcl1p5m                                                              &
  ,interp_data_3(2*pdims%i_end*pdims%j_end+1)                           &
                                           ! work_cfl
  ,interp_data_3(3*pdims%i_end*pdims%j_end+1)                           &
                                           ! work_cff
  ,icode )

END IF ! on 1.5m STASHflags

! ----------------------------------------------------------------------
!  1.5m Relative humidity
! ----------------------------------------------------------------------
! Item 245: Calculate relative humidity at 1.5m from q1p5m and t1p5m
!   Note that surface pressure is used to determine saturation
!   humidity (instead of p at 1.5m), but this is an extremely close
!   approximation.

item =   245  ! relative humidity at 1.5m
IF (icode <= 0 .AND. sf(item,sect)) THEN

!  Find humidity saturation at 1.5m, store in interp_data work.
!  Q1.5m is always specific humidity (conversion in IMP_SOLVER)
!  so we want the specific qsat
! DEPENDS ON: qsat_mix
   CALL qsat_mix(interp_data,t1p5m,p_star,pdims%i_end*pdims%j_end       &
                 ,.FALSE.)

   DO j = pdims%j_start,pdims%j_end
   DO i = pdims%i_start,pdims%i_end

!  Supersaturation (>100%) can occur with mixed phase scheme but
!  negative humidity is removed from the diagnostic:
      interp_data(i,j) = MAX( 0. , q1p5m(i,j) ) * 100.                  &
                                              / interp_data(i,j)

   END DO ! i
   END DO ! j

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,sect,im_index)),                     &
        interp_data,                                                    &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 245)"//cmessage
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
!  1.5m Temperature
! ----------------------------------------------------------------------
! Item 236 T1p5m

IF (icode <= 0 .AND. sf(236,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(236,3,im_index)),t1p5m,                   &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,236,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 236)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF
! Top of surface mixed layer (Ksurf profile)
item=356
IF (icode <= 0 .AND. sf(item,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,3,im_index)),bl_diag%smltop,         &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,item,                                                &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 356)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

!  Top of decoupled stratocu layer
item=357
IF (icode <= 0 .AND. sf(item,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,3,im_index)),bl_diag%dsctop,         &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,item,                                                &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 357)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! BL depth diagnosed as Ri>RiCrit
item=358
IF (icode <= 0 .AND. sf(item,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,3,im_index)),bl_diag%zhlocal,        &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,item,                                                &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 358)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Height of diagnosis parcel top
item=359
IF (icode <= 0 .AND. sf(item,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,3,im_index)),bl_diag%zhpar,          &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,item,                                                &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 359)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Height of decoupled layer base (=SML top if doesn't exist)
item=360
IF (icode <= 0 .AND. sf(item,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,3,im_index)),bl_diag%dscbase,        &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,item,                                                &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 360)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Height of stratocumulus cloud base
item=361
IF (icode <= 0 .AND. sf(item,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,3,im_index)),bl_diag%cldbase,        &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,item,                                                &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 361)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Parametrized entrainment rate for surface-based mixed layer
item=362
IF (icode <= 0 .AND. sf(item,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,3,im_index)),bl_diag%weparm,         &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,item,                                                &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 362)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Parametrized entrainment rate for decoupled mixed layer, SML if not
item=363
IF (icode <= 0 .AND. sf(item,3)) THEN

! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,3,im_index)),bl_diag%weparm_dsc,      &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,item,                                                &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 363)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Obukhov length
item=464
IF (icode <= 0 .AND. sf(item,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,3,im_index)),bl_diag%oblen,          &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,item,                                                &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 464)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Friction velocity
item=465
IF (icode <= 0 .AND. sf(item,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,3,im_index)),bl_diag%ustar,          &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,item,                                                &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 465)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Surface buoyancy flux
item=467
IF (icode <= 0 .AND. sf(item,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,3,im_index)),bl_diag%wbsurf,         &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,item,                                                &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 467)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Gradient Richardson number
item =   468
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%gradrich,                                               &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 468)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! Convective velocity scale
item=466
IF (icode <= 0 .AND. sf(item,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,3,im_index)),bl_diag%wstar,          &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,item,                                                &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 466)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Stratification
item =   469
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%dbdz,                                                   &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 469)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! Modulus of shear
item =   470
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%dvdzm,                                                  &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 470)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! Momentum diffusivity
item =   471
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%rhokm,                                                  &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 471)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! Heat diffusivity
item =   472
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%rhokh,                                                  &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 472)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! Turbulent kinetic energy
item =   473
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%tke,                                                    &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 473)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

!     Orographic stress (x component)
item =   474
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%ostressx,                                               &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 474)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

!     Orographic stress (y component)
item =   475
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%ostressy,                                               &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 475)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)
! Mixing length for momentum
item =   501
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%elm3d,                                                  &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 501)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! Mixing length for heat and moisture
item =   502
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%elh3d,                                                  &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 402)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! Km diffusion coeff from local scheme
item =   503
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%rhokmloc,                                               &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 503)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF !  sf(item,sect)

! Kh diffusion coeff from local scheme
item =   504
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%rhokhloc,                                               &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 504)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! Km diffusion coeff for surface-driven turb
item =   505
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%rhokmsurf,                                              &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 505)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! Kh diffusion coeff for surface-driven turb
item =   506
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%rhokhsurf,                                              &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 507)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

 ! Km diffusion coeff for cloud-top-driven turb
item =   507
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%rhokmsc,                                                &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 506)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! Kh diffusion coeff for cloud-top-driven turb
item =   508
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%rhokhsc,                                                &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 508)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! weighting applied to 1D BL scheme in Smag blending
item =   509
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%weight1d,                                               &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 508)"//cmessage
  END IF

END IF  !  sf(item,sect)

! -----------------------------------------------------------------
! Gust Diagnostic Calculation
! -----------------------------------------------------------------
item=463
IF (icode <= 0 .AND. sf(item,3)) THEN

  ALLOCATE(gust_wind(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end))
  ALLOCATE(std_dev(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end))

!       Calculation of the friction velocity; ustar currently holds
!       the magnitude of the stress.
  ustar_imp(:,:) = SQRT( SQRT(ustar_imp(:,:)) / rho1(:,:) )

!       Include the stability dependence. If the Obukhov length were
!       close to 0 the diagnostic would be ill-defined and is left
!       as 0.0.
  std_dev(:,:)     = 0.0
  WHERE( bl_diag%oblen > tol_obukhov )
    std_dev = gust_const * ustar_imp
  endwhere
  WHERE( bl_diag%oblen < -tol_obukhov )
    std_dev = gust_const * ustar_imp *                                  &
      (1.0 - ml_depth / (24.0 * bl_diag%oblen) )**(1.0/3.0)
  endwhere

  gust_wind(:,:) = SQRT( u10m_p(:,:) * u10m_p(:,:) +                    &
                         v10m_p(:,:) * v10m_p(:,:) )
  gust_wind(:,:) = gust_wind(:,:) + std_dev(:,:) * (1.0/vkman) *        &
    LOG( (5.0 * EXP(vkman * c_ugn) + z0m_eff(:,:) ) /                   &
         (5.0 + z0m_eff(:,:) ) )

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(item,3,im_index)),gust_wind,             &
         pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                 &
         atmos_im,3,item,                                               &
         icode,cmessage)
   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 463)"
      CALL ereport(routinename,icode,cmessage)
   END IF

   DEALLOCATE(gust_wind)
   DEALLOCATE(std_dev)
   DEALLOCATE(ustar_imp)

END IF    ! End of gust calculation

IF (icode <= 0 .AND. ( sf(230,3) .OR. sf(463,3) )) THEN
! deallocate these as they are no longer needed
   DEALLOCATE(u10m_p)
   DEALLOCATE(v10m_p)
END IF

! Silhouette Orographic Roughness (A/S)
! ----------------------------------------------------------------------
! Item 251

IF (sf(251,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(251,3,im_index)),                                    &
       sil_orog_land,                                                   &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

! Orographic Roughness Peak to Trough Height
! ----------------------------------------------------------------------
! Item 252

IF (sf(252,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(252,3,im_index)),                                    &
       ho2r2_orog,                                                      &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

! ----------------------------------------------------------------------
!  1.5m Specific Humidity
! ----------------------------------------------------------------------
! Item 237 q1p5m

IF (icode <= 0 .AND. sf(237,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(237,3,im_index)),q1p5m,                   &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,237,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 237)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!  1.5m Fog Fraction and Mist Fraction
! ----------------------------------------------------------------------
!  Items 248 and 253

IF (icode <= 0 .AND. (sf(248,3) .OR. sf(253,3)) ) THEN
  DO i = pdims%i_start,pdims%i_end
    DO j = pdims%j_start,pdims%j_end
      DO k=1,n_vis_thresh
        vis_threshold(i,j,1,k)=vis_thresh(k)
      END DO
    END DO
  END DO
! DEPENDS ON: fog_fr
  CALL fog_fr (p_star, rhcrit,1,                                        &
       pdims%i_end*pdims%j_end,                                         &
       t1p5m, aerosol, l_murk,                                          &
       q1p5m, qcl1p5m, qcf,                                             &
       vis_threshold,pvis,n_vis_thresh)

  IF (icode <= 0 .AND. sf(248,3)) THEN
! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(248,3,im_index)), pvis(1,1,1),           &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,248,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 248)"
      CALL ereport(routinename,icode,cmessage)
   END IF
  END IF     ! sf(248,3)

  IF (icode <= 0 .AND. sf(253,3)) THEN
! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(253,3,im_index)), pvis(1,1,2),           &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,253,                                                 &
        icode,cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 253)"
      CALL ereport(routinename,icode,cmessage)
    END IF
  END IF    ! sf(253,3)

END IF

! ----------------------------------------------------------------------
!  1.5m Dewpoint
! ----------------------------------------------------------------------
! Item 250

IF (icode <= 0 .AND. sf(250,3)) THEN

! DEPENDS ON: dewpnt
   CALL dewpnt (q1p5m, p_star, t1p5m,                                   &
               pdims%i_end*pdims%j_end, interp_data)

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(250,3,im_index)), interp_data,            &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,250,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 250)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!  Visibility in precipitation diagnostics
! ----------------------------------------------------------------------
! Items 281, 282, 283, 284 and 285.

 IF (icode <= 0 .AND. (sf(281,3) .OR. sf(282,3) .OR. sf(283,3)          &
                       .OR. sf(284,3) .OR. sf(285,3)) ) THEN

!     Calculate scattering coefficients due to precipitation

  CALL beta_precip(ls_rain, ls_snow                                     &
                 ,conv_rain, conv_snow, qcf(1,1,1)                      &
                 ,rho1, t1p5m, p_star                                   &
                 ,plsp,cca_2d,.FALSE.,.TRUE.                            &
                 ,pdims%i_end*pdims%j_end,pdims%i_end*pdims%j_end       &
                 ,1                                                     &
                 ,beta_ls_rain, beta_ls_snow                            &
                 ,beta_c_rain, beta_c_snow, icode)

END IF

IF (icode <= 0 .AND. (sf(282,3) .OR. sf(283,3)) ) THEN

!     Calculate screen level probability of vis less than thresholds

! DEPENDS ON: calc_vis_prob
  CALL calc_vis_prob(p_star,                                            &
       rhcrit,1,                                                        &
       pdims%i_end*pdims%j_end,pdims%i_end*pdims%j_end,                 &
       t1p5m,aerosol,l_murk,                                            &
       q1p5m,qcl1p5m,qcf,                                               &
       vis_thresh,n_vis_thresh,                                         &
       plsp,cca_2d,.FALSE.,                                             &
       beta_ls_rain,beta_ls_snow,                                       &
       beta_c_rain,beta_c_snow,                                         &
       vis_pr,                                                          &
       icode)
END IF

! Item 282 Prob vis<1000 m

IF(icode <= 0 .AND. sf(282,3)) THEN
! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(282,3,im_index)), vis_pr(1,1,1,1),        &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,282,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 282)"
      CALL ereport(routinename,icode,cmessage)
   END IF
END IF

! Item 283 Prob vis<5000 m

IF(icode <= 0 .AND. sf(283,3)) THEN
! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(283,3,im_index)), vis_pr(1,1,1,2),        &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,283,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 283)"
      CALL ereport(routinename,icode,cmessage)
   END IF
END IF

IF (icode <= 0 .AND.                                                    &
    (sf(281,3) .OR. sf(284,3) .OR. sf(285,3)) ) THEN

!     Visibility at 1.5 m including precipitation

! DEPENDS ON: visbty
         CALL visbty(                                                   &
           p_star, T1p5m, q1p5m, Qcl1p5m, Qcf(1,1,1)                    &
                                                      !INPUT
           ,Aerosol, calc_prob_of_vis, RHcrit(1), L_murk                &
                                                      !INPUT
           , pdims%i_end*pdims%j_end                                    &
                                                      !INPUT
           ,Vis_no_precip)                            !OUTPUT
! DEPENDS ON: vis_precip
   CALL vis_precip(vis_no_precip                                        &
               ,plsp,cca_2d,.FALSE.                                     &
               ,beta_ls_rain, beta_ls_snow                              &
               ,beta_c_rain, beta_c_snow                                &
               ,pdims%i_end*pdims%j_end,pdims%i_end*pdims%j_end,1       &
               ,vis,vis_ls_precip,vis_c_precip                          &
               ,icode)
END IF


! Item 281 Median Visibility including precip

IF(icode <= 0 .AND. sf(281,3)) THEN
! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(281,3,im_index)), vis,                    &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,281,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 281)"
      CALL ereport(routinename,icode,cmessage)
   END IF
END IF

! Item 284 Visibility in ls precip

IF(icode <= 0 .AND. sf(284,3)) THEN
! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(284,3,im_index)), vis_ls_precip,          &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,284,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 284)"
      CALL ereport(routinename,icode,cmessage)
   END IF
END IF

! Item 285 Visibility in Conv precip

IF(icode <= 0 .AND. sf(285,3)) THEN
! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(285,3,im_index)), vis_c_precip,           &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,285,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 285)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!   Sensible Heat flux: surface and model levels
! ----------------------------------------------------------------------
! Items 216,217 array ftl combines surface and model level values
!         (,,1) is the surface value.
!         Since model level fields are on rho_levels and the surface
!         is a theta_level, these are output as 2 separate diagnostics.

item =   216  ! T_L flux profile on rho_levels
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        ftl,                                                            &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 216)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! surface = ftl(,,1)
item =   217  ! T_L flux at surface
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,sect,im_index)),                     &
        ftl,                                                            &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 217)"//cmessage
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!   Sensible Moisture flux: surface and model levels
! ----------------------------------------------------------------------
! Items 222,223,241 array fqt combines surface and model level values
!         (,,1) is the surface value.
!         Since model level fields are on rho_levels and the surface
!         is a theta_level, these are output as separate diagnostics.

item =   222  ! q_T flux profile on rho_levels
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        fqt,                                                            &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 222)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! surface = fqt(,,1)
item =   223  ! T_L flux at surface
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,sect,im_index)),                     &
        fqt,                                                            &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 223)"//cmessage
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF  !  sf(item,sect)

item =   241  ! total surface moisture flux per timestep
IF (icode <= 0 .AND. sf(item,sect)) THEN

! Convert from rate to timestep accumulation explicitly
   DO j = pdims%j_start,pdims%j_end
     DO i = pdims%i_start,pdims%i_end
       interp_data(i,j) = fqt(i,j,1) * timestep
     END DO
   END DO

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(item,sect,im_index)),                     &
        interp_data,                                                    &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,item,                                             &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 241)"//cmessage
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF  !  sf(item,sect)

! ----------------------------------------------------------------------
!   Latent Heat flux.
! ----------------------------------------------------------------------
! Item 234 surf_latent_heat_flux

IF (icode <= 0 .AND. sf(234,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(234,3,im_index)),                         &
        surf_latent_heat_flux,                                          &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,234,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="error in copydiag(item 234)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!  U component of wind stress.
! ----------------------------------------------------------------------
! Item 219 taux

IF (icode <= 0 .AND. sf(219,3)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d(stashwork(si(219,3,im_index)),taux,                 &
        pdims%i_end,udims%j_end,boundary_layer_levels,0,0,0,0,          &
        at_extremity,                                                   &
        stlist(1,stindex(1,219,3,im_index)),len_stlist,                 &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,3,219,                                                 &
        icode,cmessage)


   IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 219)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!  V component of wind stress.
! ----------------------------------------------------------------------
! Item 220 tauy

IF (icode <= 0 .AND. sf(220,3)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d(stashwork(si(220,3,im_index)),tauy,                 &
        vdims%i_end,n_rows,boundary_layer_levels,0,0,0,0,               &
        at_extremity,                                                   &
        stlist(1,stindex(1,220,3,im_index)),len_stlist,                 &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,3,220,                                                 &
        icode,cmessage)


   IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 220)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!  Magnitude of wind stress on B-grid
! ----------------------------------------------------------------------
! Item 221 tauB

IF (icode <= 0 .AND. sf(221,3)) THEN

   ! Calculate stress magnitude
   DO k = 1,boundary_layer_levels
     DO j = vdims%j_start,vdims%j_end
       DO i = udims%i_start,udims%i_end
         taub(i,j,k) = SQRT(tauxb(i,j,k)*tauxb(i,j,k)                   &
                           +tauyb(i,j,k)*tauyb(i,j,k))
       END DO !i
     END DO !j
   END DO !k

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d(stashwork(si(221,3,im_index)),taub,                 &
        pdims%i_end,n_rows,boundary_layer_levels,0,0,0,0,               &
        at_extremity,                                                   &
        stlist(1,stindex(1,221,3,im_index)),len_stlist,                 &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,3,221,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag_3d(item 221)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
! Section 2.2.9  Soil temperature
! ----------------------------------------------------------------------
! Item 238  t_soil

IF (icode <= 0 .AND. sf(238,3)) THEN

   DO k = 1, dsm_levels
      DO j= pdims%j_start,pdims%j_end
         DO i = pdims%i_start,pdims%i_end
           l = i+(j-1)*pdims%i_end+(k-1)*pdims%j_end*pdims%i_end
           interp_data_3(l) = rmdi
         END DO
      END DO
      DO i = 1, land_points
         l = land_index(i) + (k-1)*pdims%j_end*pdims%i_end
         interp_data_3(l) = t_soil(i,k)
      END DO
   END DO

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d(stashwork(si(238,3,im_index)),interp_data_3,        &
        pdims%j_end,pdims%i_end,dsm_levels,0,0,0,0,                     &
        at_extremity,                                                   &
        stlist(1,stindex(1,238,3,im_index)),len_stlist,                 &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,3,238,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 238)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!  Mixed Layer depth.
! ----------------------------------------------------------------------
! Item 025 ml_depth

IF (icode <= 0 .AND. sf(025,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(025,3,im_index)),                         &
        ml_depth,                                                       &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,025,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 025)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!   wind mixing energy.
! ----------------------------------------------------------------------
! Item 224 wind_mixing_energy

IF (icode <= 0 .AND. sf(224,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(224,3,im_index)),                         &
        wind_mixing_energy,                                             &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,224,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 224)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!  Boundary Layer top.
! ----------------------------------------------------------------------
! Item 304 bl_top

IF (icode <= 0 .AND. sf(304,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(304,3,im_index)),                         &
        bl_top,                                                         &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,304,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 304)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!  Boundary Layer Type (all 6)
! ----------------------------------------------------------------------

! Item 305...310 bl_type  1..6

IF (icode <= 0 .AND. sf(pp_code_bl_types(1),3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(                                                       &
        stashwork(si(pp_code_bl_types(1),3,im_index)),                  &
        bl_type_1,                                                      &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,pp_code_bl_types(1),                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(bl_type_1)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

IF (icode <= 0 .AND. sf(pp_code_bl_types(2),3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(                                                       &
        stashwork(si(pp_code_bl_types(2),3,im_index)),                  &
        bl_type_2,                                                      &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,pp_code_bl_types(2),                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(bl_type_2)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

IF (icode <= 0 .AND. sf(pp_code_bl_types(3),3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(                                                       &
        stashwork(si(pp_code_bl_types(3),3,im_index)),                  &
        bl_type_3,                                                      &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,pp_code_bl_types(3),                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(bl_type_3)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

IF (icode <= 0 .AND. sf(pp_code_bl_types(4),3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(                                                       &
        stashwork(si(pp_code_bl_types(4),3,im_index)),                  &
        bl_type_4,                                                      &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,pp_code_bl_types(4),                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(bl_type_4)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

IF (icode <= 0 .AND. sf(pp_code_bl_types(5),3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(                                                       &
        stashwork(si(pp_code_bl_types(5),3,im_index)),                  &
        bl_type_5,                                                      &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,pp_code_bl_types(5),                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(bl_type_5)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

IF (icode <= 0 .AND. sf(pp_code_bl_types(6),3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(                                                       &
        stashwork(si(pp_code_bl_types(6),3,im_index)),                  &
        bl_type_6,                                                      &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,pp_code_bl_types(6),                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(bl_type_6)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

IF (icode <= 0 .AND. sf(pp_code_bl_types(7),3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(                                                       &
        stashwork(si(pp_code_bl_types(7),3,im_index)),                  &
        bl_type_7,                                                      &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,pp_code_bl_types(7),                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(bl_type_7)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Assign combined boundary layer type diagnostic

item=476
IF (icode <= 0 .AND. sf(item,3)) THEN

  DO j=pdims%j_start,pdims%j_end
    DO i=pdims%i_start,pdims%i_end

    bl_type_comb(i,j)=                                                  &
    1.0*bl_type_1(i,j) + 2.0*bl_type_2(i,j) + 3.0*bl_type_3(i,j)        &
  + 4.0*bl_type_4(i,j) + 5.0*bl_type_5(i,j) + 6.0*bl_type_6(i,j)        &
  + 7.0*bl_type_7(i,j)

    END DO
  END DO

! DEPENDS ON: copydiag
  CALL copydiag(stashwork(si(item,3,im_index)),bl_type_comb,            &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,item,                                                &
        icode,cmessage)

  IF (icode  >   0) THEN
     cmessage=": error in copydiag(item 476)"
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  ! item 476

! ----------------------------------------------------------------------
!   Effective roughness length for momentum
! ----------------------------------------------------------------------
! Item 026 z0m_eff

IF (icode <= 0 .AND. sf(026,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(026,3,im_index)),                         &
        z0m_eff,                                                        &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,026,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 026)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!   Effective roughness length for heat
! ----------------------------------------------------------------------
! Item 027 z0h_eff

IF (icode <= 0 .AND. sf(027,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(027,3,im_index)),                         &
        z0h_eff,                                                        &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,027,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 027)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!   Grid-box mean Vegetative roughness length for momentum
! ----------------------------------------------------------------------
! Item 028 z0m_gb

IF (icode <= 0 .AND. sf(028,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(028,3,im_index)),                         &
        z0m_gb,                                                         &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,028,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 028)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
!   RIB - measure of surface stability
! ----------------------------------------------------------------------
! Item 208 rib

IF (icode <= 0 .AND. sf(208,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(208,3,im_index)),                         &
        rib,                                                            &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,208,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 208)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
! visibility diagnostics
! ----------------------------------------------------------------------
! Item 247 visibility

IF (icode <= 0 .AND. sf(247,3)) THEN

! DEPENDS ON: visbty
         CALL visbty(                                                   &
           p_star, T1p5m, q1p5m, Qcl1p5m, Qcf(1,1,1)                    &
                                                     !INPUT
           ,Aerosol, calc_prob_of_vis, RHcrit(1), L_murk                &
                                                     !INPUT
           ,pdims%i_end*pdims%j_end                                     &
                                                 !INPUT
           ,interp_data)                         !OUTPUT

! DEPENDS ON: copydiag
         CALL copydiag(STASHwork(si(247,3,im_index)),                   &
              interp_data,                                              &
              pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,            &
              atmos_im,3,247,                                           &
              icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 247)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Surface evaporation weighted by leads
! ----------------------------------------------------------------------
! Item 232 Surface evaporation weighted by leads

IF (icode <= 0 .AND. sf(232,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(232,3,im_index)),                         &
        e_sea,                                                          &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,232,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 232)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Sensible heat flux over open sea
! ----------------------------------------------------------------------
! Item 228 Sensible heat flux over open sea

IF (icode <= 0 .AND. sf(228,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(228,3,im_index)),                         &
        h_sea,                                                          &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,228,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 228)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Evaporation amount from soil surface
! ----------------------------------------------------------------------
! Item 229 soil evaporation

IF (icode <= 0 .AND. sf(229,3)) THEN

! Convert from rate to timestep accumulation explicitly

   DO j = pdims%j_start,pdims%j_end
     DO i = pdims%i_start,pdims%i_end
       interp_data(i,j) = es(i,j) * timestep
     END DO
   END DO

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(229,3,im_index)),                         &
        interp_data,                                                    &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,229,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 229)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Sublimation ( if sea-ice included in ocean)
! ----------------------------------------------------------------------
! Item 231 Sublimation

IF (icode <= 0 .AND. sf(231,3)) THEN

! Convert from rate to timestep accumulation explicitly
   DO j = pdims%j_start,pdims%j_end
     DO i = pdims%i_start,pdims%i_end
       interp_data(i,j) = ei(i,j) * timestep
     END DO
   END DO

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(231,3,im_index)),                         &
        interp_data,                                                    &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,231,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 231)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Heat flux through sea-ice ( if sea-ice included in ocean)
! ----------------------------------------------------------------------
! Item 201 Melting

IF (icode <= 0 .AND. sf(201,3)) THEN
   sea_ice_htf_gb(:,:)=0.0
   DO j=pdims%j_start,pdims%j_end
     DO i=pdims%i_start,pdims%i_end
       DO k=1,nice
         sea_ice_htf_gb(i,j)=sea_ice_htf_gb(i,j)+                       &
                                     sea_ice_htf(i,j,k)
       END DO
     END DO
   END DO

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(201,3,im_index)),                         &
        sea_ice_htf_gb,                                                 &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,201,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 201)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
! 3d ice category heat flux through sea ice
! ----------------------------------------------------------------------
! Item 256  sea_ice_htf

IF (icode <= 0 .AND. sf(256,3)) THEN

! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(nice,len_stlist,                                 &
       stlist(1,stindex(1,256,3,im_index)),                             &
       pllnice,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 256 = sea_ice_htf)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,nice
    IF (pllnice(pslevel)) THEN
      pslevel_out=pslevel_out+1
! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(256,3,im_index)+(pslevel_out-1)        &
        *pdims%i_end*pdims%j_end),                                      &
        sea_ice_htf(1,1,pslevel),                                       &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,256,                                                 &
        icode,cmessage)
      IF (icode  >   0) THEN
          cmessage=": error in copydiag(item 256 = sea_ice_htf)"
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF
  END DO
END IF

! Melting of top of sea-ice ( if sea-ice included in ocean)
! ----------------------------------------------------------------------
! Item 235 Melting

IF (icode <= 0 .AND. sf(235,3)) THEN
   sice_mlt_htf_gb(:,:)=0.0
   DO j=pdims%j_start,pdims%j_end
     DO i=pdims%i_start,pdims%i_end
       DO k=1,nice
         sice_mlt_htf_gb(i,j)=sice_mlt_htf_gb(i,j)+                     &
                                  sice_mlt_htf(i,j,k)
       END DO
     END DO
   END DO

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(235,3,im_index)),                         &
        sice_mlt_htf_gb,                                                &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,235,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 235)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
! 3d ice category sea ice surface melt heat flux
! ----------------------------------------------------------------------
! Item 257  sice_mlt_htf

IF (icode <= 0 .AND. sf(257,3)) THEN

! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(nice,len_stlist,                                 &
       stlist(1,stindex(1,257,3,im_index)),                             &
       pllnice,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 257 = sice_mlt_htf)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,nice
    IF (pllnice(pslevel)) THEN
      pslevel_out=pslevel_out+1
! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(257,3,im_index)+(pslevel_out-1)        &
        *pdims%i_end*pdims%j_end),                                      &
        sice_mlt_htf(1,1,pslevel),                                      &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,257,                                                 &
        icode,cmessage)
      IF (icode  >   0) THEN
          cmessage=": error in copydiag(item 257 = sice_mlt_htf)"
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF
  END DO
END IF

! Heat flux due to melting of snow (Watts per sq metre).
! ----------------------------------------------------------------------
! Item 258 Melting

IF (icode <= 0 .AND. sf(258,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(258,3,im_index)),                         &
        snomlt_surf_htf,                                                &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,258,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 258)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Heat flux from surface to deep soil level 1 (Watts per sq metre).
! ----------------------------------------------------------------------
! Item 202 Surf soil flux

IF (sf(202,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(202,3,im_index)),                         &
        surf_ht_flux,                                                   &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,202,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=                                                         &
       "imp_ctl  : error in copydiag(item 202 = surf_ht_flux)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Evaporation from soil surface : rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 296 Soil evap

IF (sf(296,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(296,3,im_index)),                         &
        es,                                                             &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,296,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=                                                         &
       "imp_ctl  : error in copydiag(item 296 = es)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Evaporation from canopy : rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 297 Canopy evap

IF (sf(297,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(297,3,im_index)),                         &
        ecan,                                                           &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,297,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=                                                         &
       "imp_ctl  : error in copydiag(item 297 = ecan)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Surface sublimation : rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 298 Surf sublim

IF (sf(298,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(298,3,im_index)),                         &
        ei,                                                             &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,298,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=                                                         &
       "imp_ctl  : error in copydiag(item 298 = ei)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! TOA outgoing longwave radiation (Watts per sq metre).
! ----------------------------------------------------------------------
! Item 332 TOA outgoing LW

IF (sf(332,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(332,3,im_index)),                         &
        olr,                                                            &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,332,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=                                                         &
       "imp_ctl  : error in copydiag(item 332 = olr)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Land Mean Potential Evaporation: rate (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 334 Land potential evap

IF (sf(334,3)) THEN

   DO l=1,land_points
     epot_land(l)=0.0
     DO n=1,ntiles
       epot_land(l)=epot_land(l)+tile_frac(l,n)*epot_tile(l,n)
     END DO
   END DO

   CALL expand_from_mask (                                              &
      stashwork(si(334,3,im_index)),                                    &
       epot_land,                                                       &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)

END IF

! Potential evaporation on tiles
! ----------------------------------------------------------------------
! Item 335 Tiled potential evap

IF (sf(335,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,335,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 335 = epot_tile)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(335,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),epot_tile(1,pslevel_out),          &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Canopy conductance (m/s)
! ----------------------------------------------------------------------
! Item 259 Canopy conductance

IF (sf(259,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(259,3,im_index)),                                    &
       gs,                                                              &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

! Gross primary productivity (Kg C/m2/s)
! ----------------------------------------------------------------------
! Item 261 GPP

IF (sf(261,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(261,3,im_index)),                                    &
       gpp,                                                             &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

! Net primary productivity (Kg C/m2/s)
! ----------------------------------------------------------------------
! Item 262 NPP

IF (sf(262,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(262,3,im_index)),                                    &
       npp,                                                             &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

! Plant respiration (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 263 Plant respiration

IF (sf(263,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(263,3,im_index)),                                    &
       resp_p,                                                          &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

! Soil respiration (Kg C/m2/s)
! ----------------------------------------------------------------------
! Item 293 Soil respiration

IF (sf(293,3)) THEN
  IF (dim_cs1  ==  4) THEN
    CALL expand_from_mask (                                             &
        stashwork(si(293,3,im_index)),                                  &
         resp_s_tot,                                                    &
         land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
  ELSE
    CALL expand_from_mask (                                             &
        stashwork(si(293,3,im_index)),                                  &
         resp_s,                                                        &
         land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
  END IF
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Soil carbon (Kg C m-2)
!--------------------------------------
! ITEM 320: TOTAL SOIL CARBON CONTENT

IF (sf(320,3)) THEN
  IF (dim_cs1  ==  4) THEN
    DO i=1,land_points
      cs_tot(i) = cs(i,1) + cs(i,2) + cs(i,3) + cs(i,4)
    END DO
    CALL expand_from_mask (                                             &
        stashwork(si(320,3,im_index)),                                  &
         cs_tot,                                                        &
         land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
  ELSE
    CALL expand_from_mask (                                             &
        stashwork(si(320,3,im_index)),                                  &
         cs,                                                            &
         land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
  END IF
END IF

! Soil carbon (Kg C m-2)
!----------------------------------------------------
! ITEMS 477-480: INDIVIDUAL POOL SOIL CARBON CONTENT

! DPM
IF (sf(477,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(477,3,im_index)),                                    &
       cs(:,1),                                                         &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

! RPM
IF (sf(478,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(478,3,im_index)),                                    &
       cs(:,2),                                                         &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

! BIO
IF (sf(479,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(479,3,im_index)),                                    &
       cs(:,3),                                                         &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

!HUM
IF (sf(480,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(480,3,im_index)),                                    &
       cs(:,4),                                                         &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

! Soil respiration (Kg C m-2 s-1)
!-------------------------------------------------------------
! ITEMS 481-484: INDIVIDUAL POOL SOIL RESPIRATION

! DPM
IF (sf(481,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(481,3,im_index)),                                    &
       resp_s(:,1),                                                     &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

! RPM
IF (sf(482,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(482,3,im_index)),                                    &
       resp_s(:,2),                                                     &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

! BIO
IF (sf(483,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(483,3,im_index)),                                    &
       resp_s(:,3),                                                     &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

!HUM
IF (sf(484,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(484,3,im_index)),                                    &
       resp_s(:,4),                                                     &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

! Canopy evaporation on tiles
! ----------------------------------------------------------------------
! Item 287 Tiled canopy evap

IF (sf(287,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,287,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 287 = ecan_tile)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(287,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),ecan_tile(1,pslevel_out),          &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Transpiration + soil evaporation on tiles
! ----------------------------------------------------------------------
! Item 288 Tiled soil evap

IF (sf(288,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,288,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 288 = esoil_tile)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(288,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),esoil_tile(1,pslevel_out),         &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Surface sensible heat flux on tiles (Watts per sq metre).
! ----------------------------------------------------------------------
! Item 290 Tiled surface heat flux

IF (sf(290,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,290,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 290 = ftl_tile)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(290,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),ftl_tile(1,pslevel_out),           &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Bulk Richardson number on tiles
! ----------------------------------------------------------------------
! Item 294 Tiled RIB

IF (sf(294,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,294,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 294 = rib_tile)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(294,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),rib_tile(1,pslevel_out),           &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Surface net radiation on tiles (Watts per sq metre).
! ----------------------------------------------------------------------
! Item 314 Tiled net radiation

IF (sf(314,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,314,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 314 = radnet_tile)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(314,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),radnet_tile(1,pslevel_out),        &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Surface temperature on tiles
! ----------------------------------------------------------------------
! Item 316 Tiled TSTAR

IF (sf(316,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,316,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 316 = tstar_tile)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(316,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),tstar_tile(1,pslevel_out),         &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Surface tile fractions
! ----------------------------------------------------------------------
! Item 317 Surface fractions

IF (sf(317,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
    stlist(1,stindex(1,317,3,im_index)),                                &
    plltile,stash_pseudo_levels,num_stash_pseudo,                       &
    icode,cmessage)
  IF (icode >  0) THEN
    cmessage=                                                           &
    "imp_ctl  : error in set_pseudo_list (item 317 = tile_frac)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out = pslevel_out+1
      CALL expand_from_mask(                                            &
        stashwork(si(317,3,im_index)+(pslevel_out-1)                    &
        *pdims%i_end*pdims%j_end),tile_frac(1,pslevel_out),             &
        land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! LAI on tiles
! ----------------------------------------------------------------------
! Item 318 LAI on tiles

IF (sf(318,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(npft,len_stlist,                                 &
    stlist(1,stindex(1,318,3,im_index)),                                &
    pllpft,stash_pseudo_levels,num_stash_pseudo,                        &
    icode,cmessage)
  IF (icode >  0) THEN
    cmessage=                                                           &
    "imp_ctl  : error in set_pseudo_list (item 318 = lai_ft)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,npft
    IF (pllpft(pslevel)) THEN
      pslevel_out = pslevel_out+1
      CALL expand_from_mask(                                            &
        stashwork(si(318,3,im_index)+(pslevel_out-1)                    &
        *pdims%i_end*pdims%j_end),lai_ft(1,pslevel_out),                &
        land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Canopy height on tiles
! ----------------------------------------------------------------------
! Item 319 Canopy height on tiles

IF (sf(319,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(npft,len_stlist,                                 &
    stlist(1,stindex(1,319,3,im_index)),                                &
    pllpft,stash_pseudo_levels,num_stash_pseudo,                        &
    icode,cmessage)
  IF (icode >  0) THEN
    cmessage=                                                           &
    "imp_ctl  : error in set_pseudo_list (item 319 = canht_ft)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,npft
    IF (pllpft(pslevel)) THEN
      pslevel_out = pslevel_out+1
      CALL expand_from_mask(                                            &
        stashwork(si(319,3,im_index)+(pslevel_out-1)                    &
        *pdims%i_end*pdims%j_end),canht_ft(1,pslevel_out),              &
        land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! "Stomatal" conductance (m/s) on PFTs
! ----------------------------------------------------------------------
! Item 462 Stomatal conductance on PFTs
! Just consider PFTs as only the stomatal conductance over
! vegetated tiles is needed.
IF (sf(462,3)) THEN
  ALLOCATE(gc_tmp(land_points, npft))
  IF (l_aggregate) THEN
    gc_tmp(:,:) = rmdi
    gc_tmp(:,1) = gc(:,1)
  ELSE IF (ntiles > npft) THEN
    gc_tmp(:,:) = gc(:,1:npft)
  ELSE
    gc_tmp(:,:) = rmdi
  END IF
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(npft,len_stlist,                                 &
    stlist(1,stindex(1,462,3,im_index)),                                &
    pllpft,stash_pseudo_levels,num_stash_pseudo,                        &
    icode,cmessage)
  IF (icode >  0) THEN
    cmessage=                                                           &
    "imp_ctl  : error in set_pseudo_list (item 462 = gc)"
      CALL ereport(routinename,icode,cmessage)
  END IF
 pslevel_out=0
  DO pslevel=1,npft
    IF (pllpft(pslevel)) THEN
      pslevel_out = pslevel_out+1
      CALL expand_from_mask(                                            &
        stashwork(si(462,3,im_index)+(pslevel_out-1)                    &
        *pdims%i_end*pdims%j_end),gc_tmp(1,pslevel_out),                &
        land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
  DEALLOCATE(gc_tmp)
END IF

! Canopy water on tiles (Kg/m2).
! ----------------------------------------------------------------------
! Item 321 Tiled canopy water

IF (sf(321,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,321,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 321 = canopy)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(321,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),canopy(1,pslevel_out),             &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Canopy capacity on tiles (Kg/m2).
! ----------------------------------------------------------------------
! Item 322 Tiled canopy capacity

IF (sf(322,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,322,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 322 = catch)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(322,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),catch(1,pslevel_out),              &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Roughness length on tiles (M).
! ----------------------------------------------------------------------
! Item 324 Tiled z0

IF (sf(324,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,324,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 324 = z0m_tile)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(324,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),z0m_tile(1,pslevel_out),           &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! 1.5m temperature over tiles
! ----------------------------------------------------------------------
! Item 328 Tiled T at 1.5m

IF (sf(328,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,328,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 328 = t1p5m_tile)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(328,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),t1p5m_tile(1,pslevel_out),         &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! 1.5m specific humidity over tiles
! ----------------------------------------------------------------------
! Item 329 Tiled Q at 1.5 m

IF (sf(329,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,329,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 329 = q1p5m_tile)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(329,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),q1p5m_tile(1,pslevel_out),         &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Surface latent heat flux on tiles (Watts per sq metre).
! ----------------------------------------------------------------------
! Item 330 Tiled latent heat flux

IF (sf(330,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,330,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 330 = le_tile)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(330,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),le_tile(1,pslevel_out),            &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Sublimation on tiles
! ----------------------------------------------------------------------
! Item 331 Tiled sublimation

IF (sf(331,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,331,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 331 = ei_tile)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(331,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),ei_tile(1,pslevel_out),            &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Gross primary productiviy on plant functional types (Kg C/m2/s)
! ----------------------------------------------------------------------
! Item 289 GPP on tiles (CABLE)

IF (sf(289,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                                 &
       stlist(1,stindex(1,289,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                     &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 289 = gpp_ft)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(289,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),gpp_ft(1,pslevel_out),             &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Net primary productiviy on plant functional types (Kg C/m2/s)
! ----------------------------------------------------------------------
! Item 291 NPP on tiles (CABLE)

IF (sf(291,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                                 &
       stlist(1,stindex(1,291,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                     &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 291 = npp_ft)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(291,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),npp_ft(1,pslevel_out),             &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Plant respiration on plant functional types (Kg/m2/s)
! ----------------------------------------------------------------------
! Item 292 Plant respiration on tiles (CABLE) 

IF (sf(292,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                                 &
       stlist(1,stindex(1,292,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                     &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 292 = resp_p_ft)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(292,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),resp_p_ft(1,pslevel_out),          &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Soil moisture availability factor on plant functional types
! ----------------------------------------------------------------------
! Item 313 FSMC on PFTs

IF (sf(313,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(npft,len_stlist,                                 &
       stlist(1,stindex(1,313,3,im_index)),                             &
       pllpft,stash_pseudo_levels,num_stash_pseudo,                     &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 313 = fsmc)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,npft
    IF (pllpft(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(313,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),fsmc(1,pslevel_out),               &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Leaf turnover rate on plant functional types
! ----------------------------------------------------------------------
! Item 325 Leaf turnover rate on tiles (CABLE)       

IF (sf(325,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                                 &
       stlist(1,stindex(1,325,3,im_index)),                             &
       plltile,stash_pseudo_levels,num_stash_pseudo,                     &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 325 = g_leaf)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(325,3,im_index)+(pslevel_out-1)                  &
           *pdims%i_end*pdims%j_end),g_leaf(1,pslevel_out),             &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! CO2 total flux to atmosphere
! ----------------------------------------------------------------------
! Item 327 CO2 flux

  IF (sf(327,3)) THEN

! DEPENDS ON: copydiag
     CALL copydiag(stashwork(si(327,3,im_index)),                       &
          co2_flux_tot,                                                 &
          pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                &
          atmos_im,3,327,                                               &
          icode,cmessage)

     IF (icode  >   0) THEN
        cmessage=                                                       &
         "imp_ctl  : error in copydiag(item 327 = co2_flux_tot)"
       CALL ereport(routinename,icode,cmessage)
     END IF

  END IF

! CO2 land surface flux
! ----------------------------------------------------------------------
! Item 326 Land CO2

  IF (sf(326,3)) THEN
    CALL expand_from_mask (                                             &
        stashwork(si(326,3,im_index)),                                  &
         land_co2,                                                      &
         land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
  END IF

! Land heat flux from surface to level 1 (land mean) (W/m2)
! ----------------------------------------------------------------------
! Item 337

IF (sf(337,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(337,3,im_index)),                         &
        surf_ht_flux_land,                                              &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,337,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=                                                         &
       "imp_ctl  : error in copydiag(item 337=surf_ht_flux_land)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
! Net surface sea-ice heat flux (sea mean) (W/m2)
! ----------------------------------------------------------------------
! Item 338

IF (sf(338,3)) THEN

  surf_ht_flux_sice_gb(:,:)=0.0
  DO j=pdims%j_start,pdims%j_end
    DO i=pdims%i_start,pdims%i_end
      DO k=1,nice
        surf_ht_flux_sice_gb(i,j)=surf_ht_flux_sice_gb(i,j)+            &
                                  surf_ht_flux_sice(i,j,k)
      END DO
    END DO
  END DO

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(338,3,im_index)),                         &
        surf_ht_flux_sice_gb,                                           &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,338,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=                                                         &
       "imp_ctl  : error in copydiag(item 338=surf_ht_flux_sice_gb)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
! 3d ice category surface heat flux
! ----------------------------------------------------------------------
! Item 510  surf_ht_flux_sice

      If (icode <= 0 .and. sf(510,3)) Then

! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(nice,len_stlist,                                 &
       stlist(1,stindex(1,510,3,im_index)),                             &
       pllnice,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 510 = surf_ht_flux_sice)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,nice
    IF (pllnice(pslevel)) THEN
      pslevel_out=pslevel_out+1
! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(510,3,im_index)+(pslevel_out-1)        &
        *pdims%i_end*pdims%j_end),                                      &
        surf_ht_flux_sice(1,1,pslevel),                                 &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,510,                                                 &
        icode,cmessage)
      IF (icode  >   0) THEN
          cmessage=": error in copydiag(item 510 = surf_ht_flux_sice)"
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF
  END DO
END IF

! ----------------------------------------------------------------------
!   RIB - measure of surface stability
! ----------------------------------------------------------------------
! Item 339 rib over meaned over open sea and sea-ice

IF (sf(339,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(339,3,im_index)),                         &
        rib_ssi,                                                        &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,339,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(rib_ssi)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
! Land Mean 1.5m air temperature
! ----------------------------------------------------------------------
! Item 341 Land mean 1.5m air temperature

IF (sf(341,3)) THEN

   DO l=1,land_points
     t1p5m_land(l)=0.0
     DO n=1,ntiles
       t1p5m_land(l)=t1p5m_land(l)+tile_frac(l,n)*t1p5m_tile(l,n)
     END DO
   END DO

   CALL expand_from_mask (                                              &
      stashwork(si(341,3,im_index)),                                    &
       t1p5m_land,                                                      &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)

END IF

! ----------------------------------------------------------------------
! Land Mean 1.5m Specific Humidity
! ----------------------------------------------------------------------
! Item 342 Land mean 1.5m specific humidity

IF (sf(342,3)) THEN

   DO l=1,land_points
     q1p5m_land(l)=0.0
     DO n=1,ntiles
       q1p5m_land(l)=q1p5m_land(l)+tile_frac(l,n)*q1p5m_tile(l,n)
     END DO
   END DO

   CALL expand_from_mask (                                              &
      stashwork(si(342,3,im_index)),                                    &
       q1p5m_land,                                                      &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)

END IF
! Surface sensible heat flux  meaned over open sea and sea-ice
! ----------------------------------------------------------------------
! Item 343
! surface = ftl(,,1)

IF (sf(343,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(343,3,im_index)),                         &
        ftl_ssi,                                                        &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,343,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(ftl_ssi)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF
! Surface moisture flux  meaned over open sea and sea-ice
! ----------------------------------------------------------------------
! Item 347
! surface = fqw(,,1)
IF (sf(347,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(347,3,im_index)),                         &
        e_ssi,                                                          &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,347,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(e_ssi)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
! Sublimation over sea-ice per category (multipied by category ice fraction)
! ----------------------------------------------------------------------
! Item 509  ei_sice  (only available if nice=nice_use)

IF (icode <= 0 .and. sf(509,3)) THEN

! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(nice,len_stlist,                                 &
       stlist(1,stindex(1,509,3,im_index)),                             &
       pllnice,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 509 = ei_sice)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,nice
    IF (pllnice(pslevel)) THEN
      pslevel_out=pslevel_out+1
! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(509,3,im_index)+(pslevel_out-1)        &
        *pdims%i_end*pdims%j_end),                                      &
        ei_sice(1,1,pslevel),                                           &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,509,                                                 &
        icode,cmessage)
      IF (icode  >   0) THEN
          cmessage=": error in copydiag(item 509 = ei_sice)"
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF
  END DO
END IF

! Sublimation over sea-ice meaned over sea portion of gridbox
! ----------------------------------------------------------------------
! Item 353

IF (sf(353,3)) THEN

  ei_sice_gb(:,:)=0.0
  DO j=pdims%j_start,pdims%j_end
    DO i=pdims%i_start,pdims%i_end
      DO k=1,nice_use
        ei_sice_gb(I,J)=ei_sice_gb(I,J)+ei_sice(I,J,K)
      END DO
    END DO
  END DO
       
! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(353,3,im_index)),                         &
        ei_sice_gb,                                                     &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,353,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(ei_sice_gb)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
! Sea-ice surface net radiation (sea mean) (W/m2)
! ----------------------------------------------------------------------
! Item 381

IF (sf(381,3)) THEN

  radnet_sice_gb(:,:)=0.0
  DO J=pdims%j_start,pdims%j_end
    DO I=pdims%i_start,pdims%i_end
      DO k=1,nice_use
        radnet_sice_gb(I,J)=radnet_sice_gb(I,J)+                        &
                            radnet_sice(I,J,K)
      END DO
    END DO
  END DO


! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(381,3,im_index)),                         &
        radnet_sice_gb,                                                 &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,381,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=                                                         &
       "imp_ctl  : error in copydiag(item 381=radnet_sice_gb)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
! Land fraction
! ----------------------------------------------------------------------
! Item 395

IF (sf(395,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(395,3,im_index)),                         &
        flandg,                                                         &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,395,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=                                                         &
       "imp_ctl  : error in copydiag(item 395=flandg)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! ----------------------------------------------------------------------
! Wind shear over land
! ----------------------------------------------------------------------
! Item 389

IF (sf(389,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(389,3,im_index)),                         &
        vshr_land,                                                      &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,389,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(vshr_land)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Wind shear over sea
! ----------------------------------------------------------------------
! Item 390

IF (sf(390,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(390,3,im_index)),                         &
        vshr_ssi,                                                       &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,390,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(vshr_ssi)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! X-component of surface stress
! ----------------------------------------------------------------------
! Item 460

IF (sf(460,3)) THEN
! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(460,3,im_index)),                         &
      taux(udims%i_start:udims%i_end,udims%j_start:udims%j_end,1),      &
        pdims%i_end,udims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,460,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(taux_surf)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Y-component of surface stress
! ----------------------------------------------------------------------
! Item 461

IF (sf(461,3)) THEN
! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(461,3,im_index)),                         &
      tauy(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,1),      &
        vdims%i_end,n_rows,0,0,0,0, at_extremity,                       &
        atmos_im,3,461,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(tauy_surf)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! X-component of stress over land
! ----------------------------------------------------------------------
! Item 391

IF (sf(391,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(391,3,im_index)),                         &
        taux_land,                                                      &
        pdims%i_end,udims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,391,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(taux_land)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! X-component of stress over sea
! ----------------------------------------------------------------------
! Item 392

IF (sf(392,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(392,3,im_index)),                         &
        taux_ssi,                                                       &
        pdims%i_end,udims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,3,392,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(taux_ssi)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Y-component of stress over land
! ----------------------------------------------------------------------
! Item 393

IF (sf(393,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(393,3,im_index)),                         &
        tauy_land,                                                      &
        vdims%i_end,n_rows,0,0,0,0, at_extremity,                       &
        atmos_im,3,393,                                                 &
        icode,cmessage)


   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(tauy_land)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Y-component of stress over sea
! ----------------------------------------------------------------------
! Item 394

IF (sf(394,3)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(394,3,im_index)),                         &
        tauy_ssi,                                                       &
        vdims%i_end,n_rows,0,0,0,0, at_extremity,                       &
        atmos_im,3,394,                                                 &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(tauy_ssi)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Diagnostics required for soil moisture nudging scheme macro

! Item 50, RHOKH, exchange coefficients for moisture

IF( sf(50,3)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d(stashwork(si(50,3,im_index)),rhokh,                 &
        pdims%i_end,pdims%j_end,boundary_layer_levels,0,0,0,0,          &
        at_extremity,                                                   &
        stlist(1,stindex(1,50,3,im_index)),len_stlist,                  &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,3,50,                                                  &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(RHOKH)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Item 51, RESFS, Combined soil, stomatol and aerodynamic
!                 resistance to evaporation on tiles
IF (sf(51,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,51,3,im_index)),                              &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 51 = resfs)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(51,3,im_index)+(pslevel_out-1)                   &
           *pdims%i_end*pdims%j_end),resfs(1,pslevel_out),              &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Item 52, CHR1P5M, Ratio of coefficients required for calculation
!                 of 1.5m temperature on tiles
IF (sf(52,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,52,3,im_index)),                              &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 52 = chr1p5m)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(52,3,im_index)+(pslevel_out-1)                   &
           *pdims%i_end*pdims%j_end),chr1p5m(1,pslevel_out),            &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Item 53, ALPHA1, Gradient of saturated humidity with respect to
!                  tempearture between surface and model level 1 on
!                  tiles
IF (sf(53,3)) THEN
! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(ntiles,len_stlist,                               &
       stlist(1,stindex(1,53,3,im_index)),                              &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
  IF (icode >  0) THEN
      cmessage=                                                         &
    "imp_ctl  : error in set_pseudo_list(item 53 = alpha1)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out=0
  DO pslevel=1,ntiles
    IF (plltile(pslevel)) THEN
      pslevel_out=pslevel_out+1
      CALL expand_from_mask (                                           &
          stashwork(si(53,3,im_index)+(pslevel_out-1)                   &
           *pdims%i_end*pdims%j_end),alpha1(1,pslevel_out),             &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
    END IF
  END DO
END IF

! Item 54, RA Aerodyanmic resiatance (s/m)
! ----------------------------------------------------------------------

IF (sf(54,3)) THEN
  CALL expand_from_mask (                                               &
      stashwork(si(54,3,im_index)),                                     &
       ra,                                                              &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
END IF

! Item 55, WT_EXT
! Culmulative transpiration from soil (needed for soil moisture nudging)
! ----------------------------------------------------------------------

IF (sf(55,3)) THEN

   DO k = 1, dsm_levels
      DO j= pdims%j_start,pdims%j_end
         DO i = pdims%i_start,pdims%i_end
           l = i+(j-1)*pdims%i_end+ (k-1)*pdims%i_end*pdims%j_end
           interp_data_3(l) = rmdi
         END DO
      END DO
      DO i = 1, land_points
         l = land_index(i) + (k-1)*pdims%i_end*pdims%j_end
         interp_data_3(l) = wt_ext(i,k)
      END DO
   END DO

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d(stashwork(si(55,3,im_index)),interp_data_3,         &
        pdims%i_end,pdims%j_end,dsm_levels,0,0,0,0,                     &
        at_extremity,                                                   &
        stlist(1,stindex(1,55,3,im_index)),len_stlist,                  &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,3,55,                                                  &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage=": error in copydiag(item 55)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Item 56, resfs_land: Land mean combined soil, stomatol and aerodynamic
!                      resistance to evaporation

IF (sf(56,3)) THEN

   DO l=1,land_points
      resfs_land(l)=0.0
      DO n=1,ntiles
         resfs_land(l)=resfs_land(l)+tile_frac(l,n)*resfs(l,n)
      END DO
   END DO

   CALL expand_from_mask (                                              &
      stashwork(si(56,3,im_index)),                                     &
       resfs_land,                                                      &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)

END IF

! Item 57, chr1p5m_land: Land mean ratio of coefficients required for 
!                        calculation of 1.5m temperature

IF (sf(57,3)) THEN

   DO l=1,land_points
      chr1p5m_land(l)=0.0
      DO n=1,ntiles
         chr1p5m_land(l)=chr1p5m_land(l)+tile_frac(l,n)*chr1p5m(l,n)
      END DO
   END DO

   CALL expand_from_mask (                                              &
      stashwork(si(57,3,im_index)),                                     &
       chr1p5m_land,                                                    &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)

END IF

! Item 58, alpha1_land: Land mean Gradient of saturated humidity with 
!                       respect to, tempearture between surface and 
!                       model level 1

IF (sf(58,3)) THEN

   DO l=1,land_points
      alpha1_land(l)=0.0
      DO n=1,ntiles
         alpha1_land(l)=alpha1_land(l)+tile_frac(l,n)*alpha1(l,n)
      END DO
   END DO

   CALL expand_from_mask (                                              &
      stashwork(si(58,3,im_index)),                                     &
       alpha1_land,                                                     &
       land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)

END IF

! MGS. Extra boundary layer diagnostics required for UKCA model
! -------------------------------------------------------------------

! Item 60, RHOKH_MIX

IF( sf(60,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d(stashwork(si(60,sect,im_index)),                    &
        rhokh_mix,                                                      &
        pdims%i_end,pdims%j_end,boundary_layer_levels,0,0,0,0,          &
        at_extremity,                                                   &
        stlist(1,stindex(1,60,sect,im_index)),len_stlist,               &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,60,                                               &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(RHOKH_MIX)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Item 61, RHO_ARESIST

IF (sf(61,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(61,sect,im_index)),                       &
        rho_aresist,                                                    &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,61,                                               &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(rho_aresist)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Item 62, ARESIST [ 1/(CD_STD*VSHR) ]

IF (sf(62,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(62,sect,im_index)),                       &
        aresist,                                                        &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,62,                                               &
        icode,cmessage)

   IF (icode > 0) THEN
      cmessage="imp_ctl  : error in copydiag(aresist)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Item 63, RESIST_B

IF (sf(63,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(63,sect,im_index)),                       &
        resist_b,                                                       &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,63,                                               &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(resist_b)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Item 64, DTRDZ_CHARNEY_GRID

IF( sf(64,sect)) THEN

! DEPENDS ON: copydiag_3d
   CALL copydiag_3d(stashwork(si(64,sect,im_index)),                    &
        dtrdz_charney_grid,                                             &
        pdims%i_end,pdims%j_end,boundary_layer_levels,0,0,0,0,          &
        at_extremity,                                                   &
        stlist(1,stindex(1,64,sect,im_index)),len_stlist,               &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,64,                                               &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(dtrdz_charney_grid)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Item 65, GRID-LEVEL OF SML INVERSION

IF (sf(65,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(65,sect,im_index)),                       &
        kent,                                                           &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,65,                                               &
        icode,cmessage)

   IF (icode  >   0) THEN
      cmessage="imp_ctl  : error in copydiag(tauy_ssi)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! IF any of diagnsotics 66-68, 70-72 are selected, set interp_data_bl
! to zero
IF (sf(66,sect) .OR. sf(67,sect) .OR. sf(68,sect) .OR.                  &
    sf(70,sect) .OR. sf(71,sect) .OR. sf(72,sect)) THEN
  DO k = 1, npft
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        interp_data_bl(i,j,k) = 0.0
      END DO
    END DO
  END DO

END IF

! Item 66, Rho * entrainment rate  (we_lim)

IF (sf(66,sect)) THEN

  DO k = 1, 3
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        interp_data_bl(i,j,k) = we_lim(i,j,k)
      END DO
    END DO
  END DO

! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(npft,len_stlist,                                 &
       stlist(1,stindex(1,66,sect,im_index)),                           &
       pllbl,stash_pseudo_levels,num_stash_pseudo,                      &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage =                                                          &
    "imp_ctl  : error in set_pseudo_list(item 66 = we_lim)"
      CALL ereport(routinename,icode,cmessage)
  END IF
  pslevel_out = 0
  DO pslevel = 1, npft
    IF (pllbl(pslevel)) THEN
      pslevel_out = pslevel_out + 1
! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(66,sect,im_index)+(pslevel_out-1)      &
        *pdims%i_end*pdims%j_end),                                      &
        interp_data_bl(1,1,pslevel),                                    &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,66,                                               &
        icode,cmessage)
      IF (icode  >   0) THEN
        cmessage=": error in copydiag(item 66)"
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF
  END DO
END IF

! Item 67, Fraction of the timestep (t_frac)

IF (sf(67,sect)) THEN

  DO k = 1, 3
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        interp_data_bl(i,j,k) = t_frac(i,j,k)
      END DO
    END DO
  END DO

! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(npft,len_stlist,                                 &
       stlist(1,stindex(1,67,sect,im_index)),                           &
       pllbl,stash_pseudo_levels,num_stash_pseudo,                      &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage =                                                          &
    "imp_ctl  : error in set_pseudo_list(item 67 = t_frac)"
      CALL ereport(routinename,icode,cmessage)
  END IF

  pslevel_out = 0
  DO pslevel = 1, npft
    IF (pllbl(pslevel)) THEN
      pslevel_out = pslevel_out + 1
! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(67,sect,im_index)+(pslevel_out-1)      &
        *pdims%i_end*pdims%j_end),                                      &
        interp_data_bl(1,1,pslevel),                                    &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,67,                                               &
        icode,cmessage)
      IF (icode  >   0) THEN
        cmessage=": error in copydiag(item 67)"
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF
  END DO
END IF

! Item 68, zrzi

IF( sf(68,sect)) THEN

  DO k = 1, 3
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        interp_data_bl(i,j,k) = zrzi(i,j,k)
      END DO
    END DO
  END DO

! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(npft,len_stlist,                                 &
       stlist(1,stindex(1,68,sect,im_index)),                           &
       pllbl,stash_pseudo_levels,num_stash_pseudo,                      &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage =                                                          &
    "imp_ctl  : error in set_pseudo_list(item 68 = zrzi)"
      CALL ereport(routinename,icode,cmessage)
  END IF

  pslevel_out = 0
  DO pslevel = 1, npft
    IF (pllbl(pslevel)) THEN
      pslevel_out = pslevel_out + 1
! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(68,sect,im_index)+(pslevel_out-1)      &
        *pdims%i_end*pdims%j_end),                                      &
        interp_data_bl(1,1,pslevel),                                    &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,68,                                               &
        icode,cmessage)
      IF (icode  >   0) THEN
        cmessage=": error in copydiag(item 68)"
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF
  END DO
END IF

! Item 69, GRID-LEVEL OF DSC INVERSION (kent_dsc)

IF (sf(69,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(69,sect,im_index)),                       &
        kent_dsc,                                                       &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,69,                                               &
        icode,cmessage)

   IF (icode > 0) THEN
      cmessage="imp_ctl  : error in copydiag(kent_dsc)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Item 70, Rho * entrainment rate  (we_lim_dsc)

IF( sf(70,sect)) THEN

  DO k = 1, 3
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        interp_data_bl(i,j,k) = we_lim_dsc(i,j,k)
      END DO
    END DO
  END DO

! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(npft,len_stlist,                                 &
       stlist(1,stindex(1,70,sect,im_index)),                           &
       pllbl,stash_pseudo_levels,num_stash_pseudo,                      &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage =                                                          &
    "imp_ctl  : error in set_pseudo_list(item 70 = we_lim_dsc)"
      CALL ereport(routinename,icode,cmessage)
  END IF

  pslevel_out = 0
  DO pslevel = 1, npft
    IF (pllbl(pslevel)) THEN
      pslevel_out = pslevel_out + 1
! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(70,sect,im_index)+(pslevel_out-1)      &
        *pdims%i_end*pdims%j_end),                                      &
        interp_data_bl(1,1,pslevel),                                    &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,70,                                               &
        icode,cmessage)
      IF (icode  >   0) THEN
        cmessage=": error in copydiag(item 70)"
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF
  END DO
END IF

! Item 71, Fraction of the timestep (t_frac_dsc)

IF (sf(71,sect)) THEN

  DO k = 1, 3
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        interp_data_bl(i,j,k) = t_frac_dsc(i,j,k)
      END DO
    END DO
  END DO

! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(npft,len_stlist,                                 &
       stlist(1,stindex(1,71,sect,im_index)),                           &
       pllbl,stash_pseudo_levels,num_stash_pseudo,                      &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage =                                                          &
    "imp_ctl  : error in set_pseudo_list(item 71 = t_frac_dsc)"
      CALL ereport(routinename,icode,cmessage)
  END IF

  pslevel_out = 0
  DO pslevel = 1, npft
    IF (pllbl(pslevel)) THEN
      pslevel_out = pslevel_out + 1
! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(71,sect,im_index)+(pslevel_out-1)      &
        *pdims%i_end*pdims%j_end),                                      &
        interp_data_bl(1,1,pslevel),                                    &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,71,                                               &
        icode,cmessage)
      IF (icode > 0) THEN
        cmessage=": error in copydiag(item 71)"
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF
  END DO
END IF

! Item 72, zrzi_dsc

IF (sf(72,sect)) THEN

  DO k = 1, 3
    DO j = pdims%j_start,pdims%j_end
      DO i = pdims%i_start,pdims%i_end
        interp_data_bl(i,j,k) = zrzi_dsc(i,j,k)
      END DO
    END DO
  END DO

! DEPENDS ON: set_pseudo_list
  CALL set_pseudo_list(npft,len_stlist,                                 &
       stlist(1,stindex(1,72,sect,im_index)),                           &
       pllbl,stash_pseudo_levels,num_stash_pseudo,                      &
       icode,cmessage)
  IF (icode >  0) THEN
    cmessage =                                                          &
    "imp_ctl  : error in set_pseudo_list(item 72 = zrzi_dsc)"
      CALL ereport(routinename,icode,cmessage)
  END IF

  pslevel_out = 0
  DO pslevel = 1, npft
    IF (pllbl(pslevel)) THEN
      pslevel_out = pslevel_out + 1
! DEPENDS ON: copydiag
      CALL copydiag(stashwork(si(72,sect,im_index)+(pslevel_out-1)      &
        *pdims%i_end*pdims%j_end),                                      &
        interp_data_bl(1,1,pslevel),                                    &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,72,                                               &
        icode,cmessage)
      IF (icode > 0) THEN
        cmessage=": error in copydiag(item 72)"
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF
  END DO
END IF

! Item 73, ZHSC  Height of top of decoupled layer

IF (sf(73,sect)) THEN

! DEPENDS ON: copydiag
   CALL copydiag(stashwork(si(73,sect,im_index)),                       &
        zhsc,                                                           &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,73,                                               &
        icode,cmessage)

   IF (icode > 0) THEN
      cmessage="imp_ctl  : error in copydiag(zhsc)"
      CALL ereport(routinename,icode,cmessage)
   END IF

END IF

! Mineral Dust Diagnostics
! -------------------------------------------------------------------

! Variables for prognostic and diagnostic dust lifting

IF (l_dust .OR. l_dust_diag) THEN

! Items 74-79 Surface layer resist for dust

  DO k = 1, ndiv
    item = 73 + k
    IF (sf(item,sect)) THEN
      DO j = pdims%j_start,pdims%j_end
        DO i = pdims%i_start,pdims%i_end
          interp_data(i,j) = r_b_dust(i,j,k)
        END DO
      END DO

! DEPENDS ON: copydiag
    CALL copydiag(stashwork(si(item,sect,im_index)),                    &
        interp_data,                                                    &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,item,                                             &
        icode,cmessage)

      IF (icode > 0) THEN
        cmessage="imp_ctl  : error in copydiag(r_b_dust div  )"
        CALL ereport(routinename,icode,cmessage)
      END IF

    END IF
  END DO

  item = 400 !Dust emission fraction on tiles
  IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: set_pseudo_list
    CALL set_pseudo_list(ntiles,len_stlist,                           &
     stlist(1,stindex(1,item,sect,im_index)),                         &
     plltile,stash_pseudo_levels,num_stash_pseudo,                    &
     icode,cmessage)
    IF (icode >  0) THEN
      cmessage=                                                       &
      "IMP_CTL  : ERROR IN SET_PSEUDO_LIST(ITEM 400)"
      CALL ereport(routinename,icode,cmessage)
    END IF

    pslevel_out=0
    DO pslevel=1,ntiles
      IF (plltile(pslevel)) THEN
        pslevel_out=pslevel_out+1
        CALL expand_from_mask (                                       &
         stashwork(si(item,sect,im_index)+(pslevel_out-1)             &
         *pdims%i_end*pdims%j_end),dust_emiss_frac(1,pslevel_out),       &
         land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
      END IF
    END DO !NTILES
  END IF

  DO i = 1,ndiv
    item = 400+i  ! dust emission flux
    IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
      CALL copydiag (stashwork(si(item,sect,im_index)),                 &
        dust_flux(pdims%i_start:pdims%i_end,                            &
                  pdims%j_start:pdims%j_end,i),                         &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,item,                                             &
        icode,cmessage)
      IF (icode >  0) THEN
        cmessage=": ERROR IN COPYDIAG_3D(ITEM 401-6)"//cmessage
        CALL ereport(routinename,icode,cmessage)
      END IF

    END IF  !  SF(ITEM,SECT)
  END DO !NDIV

  DO i = 1,ndivh
    item = 410+i !threshold friction velocity on tiles
    IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: set_pseudo_list
      CALL set_pseudo_list(ntiles,len_stlist,                           &
       stlist(1,stindex(1,item,sect,im_index)),                         &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
      IF (icode >  0) THEN
        cmessage=                                                       &
        "IMP_CTL  : ERROR IN SET_PSEUDO_LIST(ITEM 411-419)"
        CALL ereport(routinename,icode,cmessage)
      END IF

      pslevel_out=0
      DO pslevel=1,ntiles
        IF (plltile(pslevel)) THEN
          pslevel_out=pslevel_out+1
          CALL expand_from_mask (                                       &
           stashwork(si(item,sect,im_index)+(pslevel_out-1)             &
           *pdims%i_end*pdims%j_end),u_s_t_tile(1,pslevel_out,i),       &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
        END IF
      END DO !NTILES
    END IF
  END DO !NDIV

  DO i = 1,ndivh
    item = 420+i !dry threshold friction velocity on tiles
    IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: set_pseudo_list
      CALL set_pseudo_list(ntiles,len_stlist,                           &
       stlist(1,stindex(1,item,sect,im_index)),                         &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
      IF (icode >  0) THEN
        cmessage=                                                       &
        "IMP_CTL  : ERROR IN SET_PSEUDO_LIST(ITEM 421-429)"
        CALL ereport(routinename,icode,cmessage)
      END IF

      pslevel_out=0
      DO pslevel=1,ntiles
        IF (plltile(pslevel)) THEN
          pslevel_out=pslevel_out+1
          CALL expand_from_mask (                                       &
           stashwork(si(item,sect,im_index)+(pslevel_out-1)             &
        *pdims%i_end*pdims%j_end),u_s_t_dry_tile(1,pslevel_out,i),      &
           land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
         END IF
       END DO !NTILES
     END IF
   END DO !NDIV

   item = 430 ! friction velocity on tiles
   IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: set_pseudo_list
     CALL set_pseudo_list(ntiles,len_stlist,                            &
       stlist(1,stindex(1,item,sect,im_index)),                         &
       plltile,stash_pseudo_levels,num_stash_pseudo,                    &
       icode,cmessage)
     IF (icode >  0) THEN
       cmessage=                                                        &
       "IMP_CTL  : ERROR IN SET_PSEUDO_LIST(ITEM 431-436)"
       CALL ereport(routinename,icode,cmessage)
     END IF

     pslevel_out=0
     DO pslevel=1,ntiles
       IF (plltile(pslevel)) THEN
         pslevel_out=pslevel_out+1
         CALL expand_from_mask (                                        &
          stashwork(si(item,sect,im_index)+(pslevel_out-1)              &
          *pdims%i_end*pdims%j_end),u_s_std_tile(1,pslevel_out),        &
          land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
       END IF
     END DO
   END IF
END IF !L_DUST

! Variables for prognostic dust only

IF (l_dust) THEN
   DO i = 1,ndiv
    item = 450+i  ! dry deposition from 2nd layer
    IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag
      CALL copydiag (stashwork(si(item,sect,im_index)),                 &
        drydep2(pdims%i_start:pdims%i_end,                              &
                pdims%j_start:pdims%j_end,i),                           &
        pdims%i_end,pdims%j_end,0,0,0,0, at_extremity,                  &
        atmos_im,sect,item,                                             &
        icode,cmessage)
      IF (icode >  0) THEN
        cmessage=": ERROR IN COPYDIAG_3D(ITEM 451-6)"//cmessage
        CALL ereport(routinename,icode,cmessage)
      END IF

    END IF  !  SF(ITEM,SECT)
  END DO !NDIV

END IF !L_DUST

! counter gradient term of taux
item =   130
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%rhogamu,                                                &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 130)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! counter gradient term of tauy
item =   131
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%rhogamv,                                                &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 131)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! counter gradient term of ftl
item =   132
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%rhogamt,                                                &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 132)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! counter gradient term of fqw
item =   133
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%rhogamq,                                                &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 133)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! mixing length
item =   134
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%elm,                                                    &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 134)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! production rate of TKE by shear
item =   135
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%tke_shr_prod,                                           &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 135)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! production rate of TKE by buoyancy
item =   136
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%tke_boy_prod,                                           &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 136)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! dissipation rate of TKE
item =   137
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%tke_dissp,                                              &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 137)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! non-dimensional diffusion coefficient for u and v (SM)
item =   138
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%sm,                                                     &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 138)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! non-dimensional diffusion coefficient for t and q (SH)
item =   139
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%sh,                                                     &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 139)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! non-gradient buoyancy flux
item =   140
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%wb_ng,                                                  &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 140)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! cloud fraction used in the TKE schemes
item =   141
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%cf_trb,                                                 &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 141)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

! condensed water used in the TKE schemes
item =   142
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%ql_trb,                                                 &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 142)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

!standard deviation of the distribution function in the TKE schemes
item =   143
IF (icode <= 0 .AND. sf(item,sect)) THEN

! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork(si(item,sect,im_index)),                  &
        bl_diag%sgm_trb,                                                &
        pdims%i_end,pdims%j_end,boundary_layer_levels,                  &
        0,0,0,0, at_extremity,                                          &
        stlist(1,stindex(1,item,sect,im_index)),len_stlist,             &
        stash_levels,num_stash_levels+1,                                &
        atmos_im,sect,item,                                             &
        icode,cmessage)

  IF (icode >  0) THEN
      cmessage=": error in copydiag_3d(item 143)"//cmessage
      CALL ereport(routinename,icode,cmessage)
  END IF

END IF  !  sf(item,sect)

IF (sf(172,sect) .OR. sf(173,sect) .OR.                                 &
    sf(176,sect) .OR. sf(177,sect) .OR.                                 &
    sf(178,sect) .OR. sf(179,sect)) THEN
  DEALLOCATE ( work3 )

END IF

!if CABLE:
call reshape_cable_diags()

! tiled soil diagnostics 
ftsoil = tsoil_tile_diag
fsmcl = SMCL_TILE_diag    
fsthf = STHF_TILE_diag   

! tiled soil temp, ftsoil (511-516)
DO i = 1,dsm_levels
   item = 510+i 
   IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: set_pseudo_list
      CALL set_pseudo_list(ntiles,len_stlist,                           &
            stlist(1,stindex(1,item,sect,im_index)),                         &
            plltile,stash_pseudo_levels,num_stash_pseudo,                    &
            icode,cmessage)
      IF (icode >  0) THEN
         cmessage=                                                       &
         "CABLE  : ERROR IN SET_PSEUDO_LIST(ITEM 511-516)"
         CALL ereport(routinename,icode,cmessage)
      END IF

      pslevel_out=0
      DO pslevel=1,ntiles
         IF (plltile(pslevel)) THEN
            pslevel_out=pslevel_out+1
            CALL expand_from_mask (                                       &
                  stashwork(si(item,sect,im_index)+(pslevel_out-1)             &
                  *pdims%i_end*pdims%j_end),ftsoil(1,pslevel_out,i),      &
                  land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
         END IF
      END DO !NTILES
   END IF
END DO !dsm_levels

! tiled soil moisture, fsmcl (517-522)
DO i = 1,dsm_levels
   item = 516+i 
   IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: set_pseudo_list
      CALL set_pseudo_list(ntiles,len_stlist,                           &
            stlist(1,stindex(1,item,sect,im_index)),                         &
            plltile,stash_pseudo_levels,num_stash_pseudo,                    &
            icode,cmessage)
      IF (icode >  0) THEN
         cmessage=                                                       &
         "CABLE  : ERROR IN SET_PSEUDO_LIST(ITEM 517-522)"
         CALL ereport(routinename,icode,cmessage)
      END IF

      pslevel_out=0
      DO pslevel=1,ntiles
         IF (plltile(pslevel)) THEN
            pslevel_out=pslevel_out+1
            CALL expand_from_mask (                                       &
                  stashwork(si(item,sect,im_index)+(pslevel_out-1)             &
                  *pdims%i_end*pdims%j_end),fsmcl(1,pslevel_out,i),      &
                  land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
            END IF
         END DO !NTILES
   END IF
END DO !dsm_levels

! tiled frozen soil moisture, fsthf (523-528)
DO i = 1,dsm_levels
   item = 522+i 
   IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: set_pseudo_list
      CALL set_pseudo_list(ntiles,len_stlist,                           &
            stlist(1,stindex(1,item,sect,im_index)),                         &
            plltile,stash_pseudo_levels,num_stash_pseudo,                    &
            icode,cmessage)
      IF (icode >  0) THEN
         cmessage=                                                       &
         "CABLE  : ERROR IN SET_PSEUDO_LIST(ITEM 523-528)"
         CALL ereport(routinename,icode,cmessage)
      END IF

      pslevel_out=0
      DO pslevel=1,ntiles
         IF (plltile(pslevel)) THEN
            pslevel_out=pslevel_out+1
            CALL expand_from_mask (                                       &
                  stashwork(si(item,sect,im_index)+(pslevel_out-1)             &
                  *pdims%i_end*pdims%j_end),fsthf(1,pslevel_out,i),      &
                  land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
         END IF
      END DO !NTILES
   END IF
END DO !dsm_levels

! tiled snow diagnostics 
fsnow_depth3L  = SNOW_DEPTH3L_diag
fsnow_mass3L   = SNOW_MASS3L_diag  
fsnow_temp3L   = SNOW_TMP3L_diag   
fsnow_rho3L    = SNOW_RHO3L_diag   
fsnow_rho1L    = SNOW_RHO1L_diag   
fsnow_age      = SNOW_AGE_diag     
fsnow_flag3L   = SNOW_FLG3L_diag   

! tiled, 3 layer snow depth, fsnow_depth3L (529-531)
DO i = 1, msn_cable 
   item = 528+i 
   IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: set_pseudo_list
      CALL set_pseudo_list(ntiles,len_stlist,                           &
            stlist(1,stindex(1,item,sect,im_index)),                         &
            plltile,stash_pseudo_levels,num_stash_pseudo,                    &
            icode,cmessage)
      IF (icode >  0) THEN
         cmessage=                                                       &
         "CABLE  : ERROR IN SET_PSEUDO_LIST(ITEM 529-531)"
         CALL ereport(routinename,icode,cmessage)
      END IF

      pslevel_out=0
      DO pslevel=1,ntiles
         IF (plltile(pslevel)) THEN
            pslevel_out=pslevel_out+1
            CALL expand_from_mask (                                       &
                  stashwork(si(item,sect,im_index)+(pslevel_out-1)             &
                  *pdims%i_end*pdims%j_end),fsnow_depth3L(1,pslevel_out,i),    &
                  land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
         END IF
      END DO !NTILES
   END IF
END DO !msn_cable

! tiled, 3 layer snow mass, fsnow_mass3L (532-534)
DO i = 1, msn_cable 
   item = 531+i 
   IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: set_pseudo_list
      CALL set_pseudo_list(ntiles,len_stlist,                           &
            stlist(1,stindex(1,item,sect,im_index)),                         &
            plltile,stash_pseudo_levels,num_stash_pseudo,                    &
            icode,cmessage)
      IF (icode >  0) THEN
         cmessage=                                                       &
         "CABLE  : ERROR IN SET_PSEUDO_LIST(ITEM 532-534)"
         CALL ereport(routinename,icode,cmessage)
      END IF

      pslevel_out=0
      DO pslevel=1,ntiles
         IF (plltile(pslevel)) THEN
            pslevel_out=pslevel_out+1
            CALL expand_from_mask (                                       &
                  stashwork(si(item,sect,im_index)+(pslevel_out-1)             &
                  *pdims%i_end*pdims%j_end),fsnow_mass3L(1,pslevel_out,i),    &
                  land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
         END IF
      END DO !NTILES
   END IF
END DO !msn_cable

! tiled, 3 layer snow temp, fsnow_temp3L (535-537)
DO i = 1, msn_cable 
   item = 534+i 
   IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: set_pseudo_list
      CALL set_pseudo_list(ntiles,len_stlist,                           &
            stlist(1,stindex(1,item,sect,im_index)),                         &
            plltile,stash_pseudo_levels,num_stash_pseudo,                    &
            icode,cmessage)
      IF (icode >  0) THEN
         cmessage=                                                       &
         "CABLE  : ERROR IN SET_PSEUDO_LIST(ITEM 535-537)"
         CALL ereport(routinename,icode,cmessage)
      END IF

      pslevel_out=0
      DO pslevel=1,ntiles
         IF (plltile(pslevel)) THEN
            pslevel_out=pslevel_out+1
            CALL expand_from_mask (                                       &
                  stashwork(si(item,sect,im_index)+(pslevel_out-1)             &
                  *pdims%i_end*pdims%j_end),fsnow_temp3L(1,pslevel_out,i),    &
                  land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
         END IF
      END DO !NTILES
   END IF
END DO !msn_cable

! tiled, 3 layer snow rho, fsnow_rho3L (538-540)
DO i = 1, msn_cable 
   item = 537+i 
   IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: set_pseudo_list
      CALL set_pseudo_list(ntiles,len_stlist,                           &
            stlist(1,stindex(1,item,sect,im_index)),                         &
            plltile,stash_pseudo_levels,num_stash_pseudo,                    &
            icode,cmessage)
      IF (icode >  0) THEN
         cmessage=                                                       &
         "CABLE  : ERROR IN SET_PSEUDO_LIST(ITEM 538-540)"
         CALL ereport(routinename,icode,cmessage)
      END IF

      pslevel_out=0
      DO pslevel=1,ntiles
         IF (plltile(pslevel)) THEN
            pslevel_out=pslevel_out+1
            CALL expand_from_mask (                                       &
                  stashwork(si(item,sect,im_index)+(pslevel_out-1)             &
                  *pdims%i_end*pdims%j_end),fsnow_rho3L(1,pslevel_out,i),    &
                  land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
         END IF
      END DO !NTILES
   END IF
END DO !msn_cable

! tiled, snow rho, fsnow_rho1L (541)
item = 541
IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: set_pseudo_list
   CALL set_pseudo_list(ntiles,len_stlist,                           &
         stlist(1,stindex(1,item,sect,im_index)),                         &
         plltile,stash_pseudo_levels,num_stash_pseudo,                    &
         icode,cmessage)
   IF (icode >  0) THEN
      cmessage=                                                       &
      "CABLE  : ERROR IN SET_PSEUDO_LIST(ITEM 541)"
      CALL ereport(routinename,icode,cmessage)
   END IF

   pslevel_out=0
   DO pslevel=1,ntiles
      IF (plltile(pslevel)) THEN
         pslevel_out=pslevel_out+1
         CALL expand_from_mask (                                       &
               stashwork(si(item,sect,im_index)+(pslevel_out-1)             &
               *pdims%i_end*pdims%j_end),fsnow_rho1L(1,pslevel_out),    &
               land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
      END IF
   END DO !NTILES
END IF

! tiled, snow age, fsnow_age (542)
item = 542
IF (icode <= 0 .AND. sf(item,sect)) THEN
! DEPENDS ON: set_pseudo_list
   CALL set_pseudo_list(ntiles,len_stlist,                           &
         stlist(1,stindex(1,item,sect,im_index)),                         &
         plltile,stash_pseudo_levels,num_stash_pseudo,                    &
         icode,cmessage)
   IF (icode >  0) THEN
      cmessage=                                                       &
      "CABLE  : ERROR IN SET_PSEUDO_LIST(ITEM 542)"
      CALL ereport(routinename,icode,cmessage)
   END IF

   pslevel_out=0
   DO pslevel=1,ntiles
      IF (plltile(pslevel)) THEN
         pslevel_out=pslevel_out+1
         CALL expand_from_mask (                                       &
               stashwork(si(item,sect,im_index)+(pslevel_out-1)             &
               *pdims%i_end*pdims%j_end),fsnow_age(1,pslevel_out),    &
               land_sea_mask,pdims%i_end*pdims%j_end,land_points_dum)
      END IF
   END DO !NTILES
END IF

!End if CABLE:


! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('DIAGNOSTICS_BL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE diagnostics_bl
