
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE BDY_LAYR-----------------------------------------------
!
!  Purpose: Calculate turbulent fluxes of heat, moisture and momentum
!           between (a) surface and atmosphere, (b) atmospheric levels
!           within the boundary layer, and/or the effects of these
!           fluxes on the primary model variables.  The flux of heat
!           into and through the soil is also modelled.  Numerous
!           related diagnostics are also calculated.
!
! Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 24.
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE bdy_layr (                                                   &
! IN  parameters for iterative SISL scheme
 numcycles, cycleno,                                                    &
! IN values defining field dimensions and subset to be processed :
 ntiles,land_pts,nice_use,bl_levels,                                    &
! IN values defining vertical grid of model atmosphere :
 p,p_theta_levels, rho_rsq, rho_wet, rho_dry, rho_tq, z_uv, z_tq,       &
! IN U, V and W momentum fields.
 u_p, v_p, u_1_px, v_1_px, u_0_px, v_0_px,                              &
! IN soil/vegetation/land surface data :
 land_index,st_levels,sm_levels,canopy,catch,catch_snow,hcon,smvccl,    &
 smvcst,smvcwt,sthf,sthu,sil_orog_land,ho2r2_orog,sd_orog,              &
! IN for dust scheme
 soil_layer_moisture,                                                   &
! IN sea/sea-ice data :
 ice_fract_cat, k_sice,                                                 &
! IN cloud data :
 cf,q,qcf,qcl,t,                                                        &
! IN everything not covered so far :
 co2_mmr,photosynth_act_rad,pstar,                                      &
 rad_hr,micro_tends,lq_mix_bl,zh_prev,ddmfx,nlcl,zhpar,z_lcl,           &
! IN SCM variables
 flux_e, flux_h, l_flux_bc, l_spec_z0, z0m_scm, z0h_scm,                &
! IN variables required for CLASSIC aerosol/mineral dust schemes
 l_aero_classic,l_dust,l_dust_diag,soil_clay,soil_sand,                 &
 dust_mrel1,dust_mrel2,dust_mrel3,dust_mrel4,dust_mrel5,dust_mrel6,     &
! IN additional variables for JULES
 snow_tile,z0_tile,z0h_tile_bare,lw_down,                               &
 sw_tile,t_soil,ti,tstar,tstar_tile,                                    &
 co2_3d,l_co2_interactive,l_phenol,l_triffid,asteps_since_triffid,      &
 cs,frac,canht_ft,lai_ft,fland,flandg,z_land,albsoil,cos_zenith_angle,  &
! IN STASH flags
 sfme,slh,sq1p5,st1p5,su10,sv10,sz0heff,                                &
! SCM Diagnostics (dummy values in full UM) and stash diags
 nSCMDpkgs,L_SCMDiags, BL_diag,                                         &
! INOUT data :
 Gs,z0msea,w, etadot,tstar_sea,tstar_sice_cat,                          &
 zh, ntpar,l_shallow,cumulus,error,                                     &
! INOUT additional variables for JULES
 g_leaf_acc,npp_ft_acc,resp_w_ft_acc,resp_s_acc,l_q10,                  &
! INOUT variables on TKE based turbulence schemes
 e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                          &
! OUT variables for message passing
 flandfac, fseafac, rhokm_land, rhokm_ssi, cdr10m, tau_fd_x, tau_fd_y,  &
 rhogamu, rhogamv, f_ngstress,                                          &
! OUT Diagnostic not requiring STASH flags :
 e_sea,fqw,fqw_tile,epot_tile, ftl,ftl_tile,h_sea,rhokh,rhokm,          &
 rib_gb,vshr,zht, shallowc,cu_over_orog,                                &
 bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6,           &
 bl_type_7, z0m_eff_gb,z0h_eff_gb,                                      &
! OUT diagnostic requiring STASH flags :
 fme,                                                                   &
! OUT diagnostics required for soil moisture nudging scheme :
 wt_ext,ra,                                                             &
! OUT data required for tracer mixing :
 rho_aresist,aresist,resist_b, ntml, kent, we_lim, t_frac, zrzi,        &
 kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,                      &
!OUT variables required for mineral dust scheme
 r_b_dust,dust_flux,dust_emiss_frac,                                    &
 u_s_t_tile,u_s_t_dry_tile,u_s_std_tile,                                &
! OUT variables required in IMP_SOLVER
 alpha1,ashtf,bq_gb,bt_gb, dtrdz_charney_grid,rdz_charney_grid,         &
 dtrdz_u,dtrdz_v,rdz_u,rdz_v,                                           &
 fraca,rhokh_tile,smc,chr1p5m,resfs,z0hssi,z0mssi, z1_tq,uStarGBM,      &
! OUT additional variables for JULES
 le_tile,radnet_sice,radnet_tile,rib_tile,rho_aresist_tile,aresist_tile,&
 resist_b_tile,alpha1_sice,ashtf_tile,fqw_ice,                          &
 ftl_ice,resft,rhokh_sice,z0h_tile,z0m_gb,z0m_tile,chr1p5m_sice,        &
 g_leaf,gpp_ft,npp_ft, resp_p_ft,resp_s,resp_s_tot,resp_w_ft,           &
 gc,canhc_tile,wt_ext_tile,flake, tile_index,tile_pts,tile_frac,fsmc,   &
 rib_ssi,vshr_land,vshr_ssi,tstar_land,tstar_ssi,                       &
 rhokh_mix,dtstar_tile,dtstar,hcons,emis_tile,emis_soil,                &
! OUT data required elsewhere in UM system :
 gpp,npp,resp_p,t1_sd,q1_sd,ntdsc,nbdsc,wstar,wthvs,uw0,vw0             &
 )

  USE atm_fields_bounds_mod, ONLY:                                      &
   udims, vdims, tdims, qdims, pdims, pdims_s, pdims_l, wdims
  USE atmos_constants_mod, ONLY : cp
  USE water_constants_mod, ONLY: lc
  USE atm_step_local, ONLY: dim_cs1, dim_cs2, land_pts_trif, npft_trif, &
       co2_dim_len,co2_dim_row
  USE bl_diags_mod, ONLY : strnewbldiag
  USE bl_option_mod, ONLY: on
  USE nstypes, ONLY: ntype, npft
  USE switches, ONLY: IScrnTDiag, i_modiscopt
  USE earth_constants_mod, ONLY: g
  USE dust_parameters_mod, ONLY: ndiv, ndivh, ndivl
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  use cable_data_mod, only : cable_control3

  IMPLICIT NONE

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
!---------------------------------------------------------------------
!  Inputs :-
!---------------------------------------------------------------------
! (a) Defining horizontal grid and subset thereof to be processed.
  INTEGER, INTENT(IN) ::                                                &
   numcycles,                                                           &
                  ! Number of cycles (iterations) for iterative SISL.
   cycleno
                  ! Iteration no

! Switch for calculating exchange coeffs from latest values.
  LOGICAL, INTENT(IN) ::                                                &
   lq_mix_bl

  INTEGER, INTENT(IN) ::                                                &
   ntiles,                                                              &
                                   ! IN No. of land-surface tiles
   land_pts,                                                            &
                                   ! IN No.of land points in whole grid.
   nice_use,                                                            &
                                   ! IN No. of sea ice categories used fully
                                   !    in surface calculations
   bl_levels,                                                           &
                                   ! IN Max. no. of "boundary" levels
   nlcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! IN No of levels to LCL

! (b) Defining vertical grid of model atmosphere.
  REAL, INTENT(IN) ::                                                   &
    p(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,      &
      pdims_s%k_start:bl_levels+1),                                     &
    p_theta_levels(tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end, bl_levels+1),             &
    rho_rsq(pdims_s%i_start:pdims_s%i_end,                              &
            pdims_s%j_start:pdims_s%j_end,pdims_s%k_start:bl_levels+1), &
                                         ! IN Density * R**2
    rho_wet(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels+1),                                               &
                        ! wet density on rho levels (kg/m3)
    rho_dry(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
                pdims%k_start:bl_levels+1),                             &
                        ! dry density on rho levels (kg/m3)
  rho_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                ! density on TQ (ie. theta) levels;
                                !    used in RHOKM so wet density
  z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                    ! Z_uv(*,K) is height of half
                                    ! level k-1/2.
  z_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                ! Z_tq(*,K) is height of full level k.
    u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
    v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
  u_1_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),  &
  v_1_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),  &
  u_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),  &
  v_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

! (c) Soil/vegetation/land surface parameters (mostly constant).
  INTEGER, INTENT(IN) ::                                                &
   land_index(land_pts)      ! IN LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.

  INTEGER, INTENT(IN) ::                                                &
   st_levels,                                                           &
                                   ! IN No. of deep soil temp. levels
   sm_levels                   ! IN No. of soil moisture levels

  REAL, INTENT(IN) ::                                                   &
   canopy(land_pts,ntiles),                                             &
                                   ! IN Surface/canopy water for
                                   !    snow-free land tiles (kg/m2)
   catch(land_pts,ntiles),                                              &
                                   ! IN Surface/canopy water capacity
                                   !    of snow-free land tiles (kg/m2).
   catch_snow(land_pts,ntiles),                                         &
                                   ! IN Snow interception capacity of
                                   !    NLT tile (kg/m2).
   hcon(land_pts),                                                      &
                                 ! IN Soil thermal conductivity
!                                     (W/m/K).
   smvccl(land_pts,sm_levels),                                          &
                                 ! IN Critical volumetric SMC (m3/m3
!                                     of soil).
   smvcst(land_pts,sm_levels),                                          &
                                 ! IN Volumetric saturation point
!                                     (m3/m3 of soil).
   smvcwt(land_pts,sm_levels),                                          &
                                 ! IN Volumetric wilting point (m3/m3
!                                     of soil).
   sthf(land_pts,sm_levels),                                            &
                                 ! IN Frozen soil moisture content of
!                                     each layer as a fraction of
!                                     saturation.
   sthu(land_pts,sm_levels),                                            &
                                 ! IN Unfrozen soil moisture content
!                                     of each layer as a fraction of
!                                     saturation.
   sil_orog_land(land_pts),                                             &
                                 ! IN Silhouette area of unresolved
!                                     orography per unit horizontal area
!                                     on land points only.
   ho2r2_orog(land_pts),                                                &
                                 ! IN Standard Deviation of orography.
!                                     equivilent to peak to trough
!                                     height of unresolved orography
!                                     devided by 2SQRT(2) on land
!                                     points only (m)
   sd_orog(land_pts),                                                   &
                                 ! IN Standard Deviation of unresolved 
!                                     orography on land points only (m)
   soil_layer_moisture(land_pts,sm_levels),                             &
                                 !IN soil moisture per layer (kg m-2)
   zh_prev(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                  ! IN boundary layer height from previous timestep
   ddmfx(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                  ! IN Convective downdraught mass-flux at cloud base
    zhpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                 ! IN Height of top of initial
                                 !     parcel ascent
    z_lcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! IN Height of LCL

! (d) Sea/sea-ice data.
  REAL, INTENT(IN) ::                                                   &
   ice_fract_cat(pdims%i_start:pdims%i_end,                             &
                 pdims%j_start:pdims%j_end,nice_use),                   &
                                  ! IN Fraction of gridbox covered by
!                                      category sea-ice (decimal fraction).
                       ! If nice_use=1, this is the sum of the categories
   k_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use)
                                   ! IN sea ice surface layer effective 
!                                  !    conductivity (W/m2/K)
! (e) Cloud data.
  REAL, INTENT(IN) ::                                                   &
   cf(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,bl_levels),   &
                                          ! IN Cloud fraction (decimal).
   qcf(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,bl_levels),  &
                                          ! IN Cloud ice (kg per kg air)
   qcl(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,bl_levels),  &
                                          ! IN Cloud liquid water
   q(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,bl_levels),    &
                                          ! IN specific humidity
   t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                          ! IN temperature

! (f) Atmospheric + any other data not covered so far, incl control.
  REAL, INTENT(IN) ::                                                   &
   co2_mmr,                                                             &
                                   ! IN CO2 Mass Mixing Ratio
   photosynth_act_rad(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end),                       &
                                   ! IN Net downward shortwave radiation
                                   !    in band 1 (w/m2).
   pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! IN Surface pressure (Pascals).
   rad_hr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          2,bl_levels),                                                 &
                                    ! IN (LW,SW) rad heating rate (K/s)
    micro_tends(qdims%i_start:qdims%i_end,qdims%j_start:qdims%j_end,    &
                2, bl_levels),                                          &
                           ! Tendencies from microphys within BL levels
                           ! (TL, K/s; QW, kg/kg/s)
    soil_clay (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
    soil_sand (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
    dust_mrel1 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
    dust_mrel2 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
    dust_mrel3 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
    dust_mrel4 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
    dust_mrel5 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
    dust_mrel6 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
   flux_e(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                   ! IN Surf. lat. heat flux   (W/m^2)
   flux_h(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                   ! IN Surf. sens. heat flux  (W/m^2)
   tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
                                   ! IN Surface temperature (K).
   tstar_tile(land_pts,ntiles),                                         &
                                   ! IN Surface tile temperatures
   z_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                   ! IN    Land height (m).
   z0m_scm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
   z0h_scm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                   ! IN   Fixed sea-surface roughness
                                   !      lengths for momentum and
                                   !      scalars (m, from SCM namelist)
   t_soil(land_pts,sm_levels),                                          &
                                   ! IN Soil temperatures (K).
   ti(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                   ! IN Sea-ice surface layer

  LOGICAL, INTENT(IN) ::                                                &
   l_flux_bc,                                                           &
                                   ! IN Switches for prescribed surface
   l_spec_z0,                                                           &
                                 !    fluxes and roughness lengths
   l_aero_classic,                                                      &
                                 ! IN Switch for CLASSIC aerosol scheme
   l_dust,                                                              &
                                 ! IN Switch for prognostic mineral dust
   l_dust_diag
                                 ! IN Switch for diagnostic mineral dust
                                 !    lifting

! Additional JULES variables
  INTEGER, INTENT(IN) ::                                                &
   asteps_since_triffid
                                   ! IN Number of atmospheric
                                   !    timesteps since last call
                                   !    to TRIFFID.
  REAL, INTENT(IN) ::                                                   &
   snow_tile(land_pts,ntiles),                                          &
                                   ! IN Lying snow on tiles (kg/m2)
   z0_tile(land_pts,ntiles),                                            &
                                   ! IN Tile roughness lengths (m).
   z0h_tile_bare(land_pts,ntiles),                                      &
                                   ! IN Tile thermal roughness 
                                   ! lengths (m) without snow.
   lw_down(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                   ! IN Surface downward LW radiation
                                   !    (W/m2).
   sw_tile(land_pts,ntiles),                                            &
                                   ! IN Surface net SW radiation on
                                   !    land tiles (W/m2).
   co2_3d(co2_dim_len,co2_dim_row),                                     &
                                   ! IN 3D CO2 field if required.
   cs(land_pts,dim_cs1),                                                &
                              ! IN Soil carbon (kg C/m2).
   frac(land_pts,ntype),                                                &
                                   ! IN Fractions of surface types.
   fland(land_pts),                                                     &
                                   ! IN Land fraction on land tiles.
   flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end), &
                                   ! IN Land fraction on all tiles.
   canht_ft(land_pts,npft),                                             &
                                   ! IN Canopy height (m)
   lai_ft(land_pts,npft),                                               &
                                   ! IN Leaf area index
   albsoil(land_pts),                                                   &
                                   ! Soil albedo.
   cos_zenith_angle(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                   ! Cosine of the zenith angle

  LOGICAL, INTENT(IN) ::                                                &
   l_co2_interactive,                                                   &
                                   ! IN Switch for 3D CO2 field
   l_phenol,                                                            &
                                   ! IN Indicates whether phenology
                                   !    in use
   l_triffid,                                                           &
                                   ! IN Indicates whether TRIFFID
                                   !    in use.
   sfme,                                                                &
               ! IN Flag for FME (q.v.).
   sz0heff,                                                             &
               ! IN Flag for Z0H_EFF
   slh,                                                                 &
               ! IN Flag for LATENT_HEAT (q.v.)
   sq1p5,                                                               &
               ! IN Flag for Q1P5M (q.v.)
   st1p5,                                                               &
               ! IN Flag for T1P5M (q.v.)
   su10,                                                                &
               ! IN Flag for U10M (q.v.)
   sv10    ! IN Flag for V10M (q.v.)

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER, INTENT(IN) ::                                                &
   nSCMDpkgs             ! No of SCM diagnostics packages

  LOGICAL, INTENT(IN) ::                                                &
   L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages
!     Declaration of new BL diagnostics.
  TYPE (strnewbldiag) :: BL_diag
!---------------------------------------------------------------------
!  In/outs :-
!---------------------------------------------------------------------
  REAL, INTENT(INOUT) ::                                                &
   Gs(land_pts),                                                        &
                                   ! INOUT "Stomatal" conductance to
                                   !        evaporation (m/s).
   tstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),      &
                                   ! INOUT  Open sea sfc temperature (K).
   tstar_sice_cat(tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,nice_use),                  &
                                   ! INOUT Sea-ice sfc temperature (K).
    z0msea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                   ! INOUT Sea-surface roughness
                                   !       length for momentum (m).
                                   !       NB: same storage is used
                                   !       for Z0V, so the intent is
                                   !       IN for land points.
    w(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,0:bl_levels), &
    etadot(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,         &
           0:bl_levels),                                                &
    zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! INOUT Height above surface of top of
!                                   boundary layer (metres).
  INTEGER, INTENT(INOUT) ::                                             &
   error                    ! OUT 0 - AOK;
                                !     1 to 7  - bad grid definition
                                !     detected;
! INOUT variables on TKE based turbulence schemes
  REAL, INTENT(INOUT) ::                                                &
   e_trb(pdims_l%i_start:pdims_l%i_end,pdims_l%j_start:pdims_l%j_end,   &
           bl_levels),                                                  &
                    ! TKE defined on theta levels K-1
   tsq_trb(pdims_l%i_start:pdims_l%i_end,pdims_l%j_start:pdims_l%j_end, &
           bl_levels),                                                  &
                    ! Self covariance of liquid potential temperature
                    ! (thetal'**2) defined on theta levels K-1
   qsq_trb(pdims_l%i_start:pdims_l%i_end,pdims_l%j_start:pdims_l%j_end, &
           bl_levels),                                                  &
                    ! Self covariance of total water
                    ! (qw'**2) defined on theta levels K-1
   cov_trb(pdims_l%i_start:pdims_l%i_end,pdims_l%j_start:pdims_l%j_end, &
           bl_levels),                                                  &
                    ! Correlation between thetal and qw
                    ! (thetal'qw') defined on theta levels K-1
   zhpar_shcu(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                    ! Height of mixed layer used to evaluate
                    ! the non-gradient buoyancy flux

  LOGICAL, INTENT(INOUT) ::                                             &
    cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                   ! INOUT Logical switch for trade Cu
    l_shallow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                           ! INOUT Flag to indicate shallow convection
   l_q10                       ! INOUT Indicates Q10 for soil resp'n

  INTEGER, INTENT(INOUT) ::                                             &
   ntpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! INOUT Top level of initial parcel
                                 !  ascent. Used in convection scheme.

! Additional JULES variables
  REAL, INTENT(INOUT) ::                                                &
   g_leaf_acc(land_pts,npft),                                           &
                                   ! INOUT Accumulated G_LEAF
   npp_ft_acc(land_pts_trif,npft_trif),                                 &
                                   ! INOUT Accumulated NPP_FT
   resp_w_ft_acc(land_pts_trif,npft_trif),                              &
                                   ! INOUT Accum RESP_W_FT
   resp_s_acc(land_pts_trif,dim_cs1) ! INOUT Accumulated RESP_S
!---------------------------------------------------------------------
!  Outputs :-
!---------------------------------------------------------------------
!  (a) Calculated anyway (use STASH space from higher level) :-
  REAL, INTENT(OUT) ::                                                  &
   e_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                   ! OUT Evaporation from sea times
!                                     leads fraction. Zero over land.
!                                     (kg per square metre per sec).
   fqw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                   ! OUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
   fqw_tile(land_pts,ntiles),                                           &
                                   ! OUT Surface FQW for land tiles
   epot_tile(land_pts,ntiles),                                          &
                                   ! OUT Local EPOT for land tiles.
   ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                   ! OUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
   ftl_tile(land_pts,ntiles),                                           &
                                   ! OUT Surface FTL for land tiles
   h_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                   ! OUT Surface sensible heat flux over
!                                     sea times leads fraction. (W/m2)
   rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                   ! OUT Exchange coeffs for moisture.
   rhokm(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,   &
         bl_levels),                                                    &
                                   ! OUT Exchange coefficients for
                                   !     momentum on P-grid
   rib_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                   ! OUT Mean bulk Richardson number for
!                                     lowest layer.
   rib_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                   ! OUT Sea mean bulk Richardson no.
!                                        for lowest layer.
   vshr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                   ! OUT Magnitude of surface-to-lowest
!                                     atm level wind shear (m per s).
   vshr_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                   ! OUT Magnitude of land sfc-to-lowest
!                                     atm level wind shear (m per s).
   vshr_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                   ! OUT Mag. of mean sea sfc-to-lowest
!                                     atm level wind shear (m per s).
    zht(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                   ! OUT Max height of turb mixing
   bl_type_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if stable
                                   !     b.l. diagnosed, 0.0 otherwise.
   bl_type_2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if Sc over
                                   !     stable surface layer diagnosed,
                                   !     0.0 otherwise.
   bl_type_3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if well
                                   !     mixed b.l. diagnosed,
                                   !     0.0 otherwise.
   bl_type_4(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if
                                   !     decoupled Sc layer (not over
                                   !     cumulus) diagnosed,
                                   !     0.0 otherwise.
   bl_type_5(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if
                                   !     decoupled Sc layer over cumulus
                                   !     diagnosed, 0.0 otherwise.
   bl_type_6(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if a
                                   !     cumulus capped b.l. diagnosed,
                                   !     0.0 otherwise.
   bl_type_7(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if a
                                   !     Shear-dominated unstable b.l.
                                   !     diagnosed, 0.0 otherwise.
   rho_aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                   ! OUT RHOSTAR*CD_STD*VSHR
!                                        for CLASSIC aerosol scheme
   aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                   ! OUT 1/(CD_STD*VSHR)
!                                        for CLASSIC aerosol scheme
   resist_b(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                   ! OUT (1/CH-1/(CD_STD)/VSHR
!                                         for CLASSIC aerosol scheme
    we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                                    ! OUT rho*entrainment rate implied b
                                    !     placing of subsidence
    zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                                    ! OUT (z-z_base)/(z_i-z_base)
    t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                                    ! OUT a fraction of the timestep
    we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),  &
                                    ! OUT rho*entrainment rate implied b
                                    !     placing of subsidence
    zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),    &
                                    ! OUT (z-z_base)/(z_i-z_base)
    t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),  &
                                    ! OUT a fraction of the timestep
    zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! OUT Top of decoupled layer
    r_b_dust(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,ndiv), &
                                      !OUT surface layer resist for dust
    wstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                   ! OUT Convective velocity scale (m/s)
    wthvs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                   ! OUT surface flux of thv (Km/s)
    shallowc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                   ! OUT Shallow Cu diagnostic
                                   !   Indicator set to 1.0 if shallow,
                                   !   0.0 if not shallow or not cumulus
    cu_over_orog(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),  &
                                   ! OUT Indicator for cumulus
                                   !     over steep orography
                                   !   Indicator set to 1.0 if true,
                                   !   0.0 if false. Exclusive.
   wt_ext(land_pts,sm_levels),                                          &
                                  !OUT cumulative fraction of transp'n
   ra(land_pts),                                                        &
                                  !OUT Aerodynamic resistance (s/m)
   tstar_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                   ! OUT   Land mean sfc temperature (K)
   tstar_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                   ! OUT Sea mean sfc temperature (K).
!                                      temperature (K).

REAL, INTENT(OUT) ::                                                    &
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
          pdims_s%j_start:pdims_s%j_end ,bl_levels),                    &
  rhogamv(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end ,bl_levels),                    &
  f_ngstress(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end,2:bl_levels)

  INTEGER, INTENT(OUT) ::                                               &
   ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                   ! OUT Top level for turb mixing in
!                                           any decoupled Sc layer
   nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                   ! OUT Bottom level of any decoupled
                                   !     turbulently-mixed Sc layer.
    kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                   ! OUT grid-level of SML inversion
    kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                   ! OUT grid-level of DSC inversion
   ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! OUT Number of model levels in the
                                 ! surface-based turbulently mixed layer.

!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-
  REAL, INTENT(OUT) ::                                                  &
   fme(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                     ! OUT Wind mixing "power" (W per sq m).

!-2 Genuinely output, needed by other atmospheric routines :-
  REAL, INTENT(OUT) ::                                                  &
   gpp(land_pts),                                                       &
                                 ! OUT Gross primary productivity
                                 !     (kg C/m2/s).
   npp(land_pts),                                                       &
                                 ! OUT Net primary productivity
                                 !     (kg C/m2/s).
   resp_p(land_pts),                                                    &
                                 ! OUT Plant respiration (kg C/m2/s).
   t1_sd(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                 ! OUT Standard deviation of turbulent
!                                   fluctuations of layer 1 temperature;
!                                   for use in initiating convection.
   q1_sd(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                 ! OUT Standard deviation of turbulent
!                                   fluctuations of layer 1 humidity;
!                                   for use in initiating convection.
   z0m_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
   z0h_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! OUT Effective grid-box roughness
!                                   lengths for momentum and for
!                                   heat, moisture
  REAL, INTENT(OUT) ::                                                  &
    uw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                        !OUT U-component of surface wind stress (P-grid)
    vw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                        ! V-component of surface wind stress (P-grid)

! OUT variables for IMP_SOLVER (that used to be local arrays)
  REAL, INTENT(OUT) ::                                                  &
   alpha1(land_pts,ntiles),                                             &
                                  ! OUT Mean gradient of saturated
                                  !     specific humidity with respect
                                  !     to temperature between the
                                  !     bottom model layer and tile
                                  !     surfaces
   ashtf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use), &
                                  ! OUT Coefficient to calculate surface
!                                 heat flux into sea-ice.
   dtrdz_charney_grid(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end,bl_levels),             &
                                  ! OUT dt/(rho*r*r*dz) for scalar
                                  !     flux divergence
   rdz_charney_grid(pdims%i_start:pdims%i_end,                          &
                    pdims%j_start:pdims%j_end,bl_levels),               &
                                ! OUT RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
   dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                           ! OUT dt/(rho*r*r*dz) for
   dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                           ! OUT U,V flux divergence
   bq_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                ! OUT grid-box mean buoyancy parameter
                                !     on p,T,q-levels (full levels).
   bt_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                ! OUT grid-box mean buoyancy parameter
                                !     on p,T,q-levels (full levels).
   rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
          2:bl_levels),                                                 &
                                ! OUT RDZ (K > 1) on UV-grid.
   rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
          2:bl_levels),                                                 &
                                ! OUT RDZ (K > 1) on UV-grid.
   fraca(land_pts,ntiles),                                              &
                                ! OUT Fraction of surface moisture
                                !     flux with only aerodynamic
                                !     resistance for snow-free land
                                !     tiles.
   rhokh_tile(land_pts,ntiles),                                         &
                                ! OUT Surface exchange coefficients
                                !     for land tiles
   smc(land_pts),                                                       &
                                ! OUT Available moisture in the
                                !     soil profile (mm).
   chr1p5m(land_pts,ntiles),                                            &
                                ! OUT Ratio of coefffs for
                                !     calculation of 1.5m temp for
                                !     land tiles.
   resfs(land_pts,ntiles),                                              &
                                ! OUT Combined soil, stomatal
                                !     and aerodynamic resistance
                                !     factor for fraction (1-FRACA)
                                !     of snow-free land tiles.
   z0hssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
   z0mssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                ! OUT Roughness lengths over sea (m).
   z1_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                         ! OUT Height of lowest theta level.

  REAL, INTENT(OUT) :: uStarGBM(pdims%i_start:pdims%i_end,              &
                                pdims%j_start:pdims%j_end)
        ! GBM surface friction velocity for diagnosis of decoupling

! Additional JULES variables
  INTEGER, INTENT(OUT) ::                                               &
   tile_index(land_pts,ntype),                                          &
                                   ! OUT Index of tile points
   tile_pts(ntype)             ! OUT Number of tile points

  REAL, INTENT(OUT) ::                                                  &
   le_tile(land_pts,ntiles),                                            &
                                   ! OUT Surface latent heat flux for
                                   !     land tiles
   radnet_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,     &
               nice_use),                                               &
                                   ! OUT Surface net radiation on
                                   !     sea-ice (W/m2)
   radnet_tile(land_pts,ntiles),                                        &
                                   ! OUT Surface net radiation on
                                   !     land tiles (W/m2)
   rib_tile(land_pts,ntiles),                                           &
                                   ! OUT RIB for land tiles.
   rho_aresist_tile(land_pts,ntiles),                                   &
                                   ! OUT RHOSTAR*CD_STD*VSHR on land
                                   !     tiles for CLASSIC aerosol scheme
   aresist_tile(land_pts,ntiles),                                       &
                                   ! OUT 1/(CD_STD*VSHR) on land tiles
                                   !     for CLASSIC aerosol scheme
   resist_b_tile(land_pts,ntiles),                                      &
                                   ! OUT (1/CH-1/CD_STD)/VSHR on land
                                   !     tiles for CLASSIC aerosol scheme
   alpha1_sice(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,nice_use),                     &
                                   ! OUT ALPHA1 for sea-ice.
   ashtf_tile(land_pts,ntiles),                                         &
                                   !OUT Coefficient to calculate
                                   !     surface heat flux into land
                                   !     tiles.
   fqw_ice(pdims%i_start:pdims%i_end,                                   &
           pdims%j_start:pdims%j_end,nice_use)                 ,        &
                                   ! OUT Surface FQW for sea-ice
   ftl_ice(pdims%i_start:pdims%i_end,                                   &
           pdims%j_start:pdims%j_end,nice_use),                         &
                                   ! OUT Surface FTL for sea-ice
   resft(land_pts,ntiles),                                              &
                                   ! OUT Total resistance factor.
                                   !     FRACA+(1-FRACA)*RESFS for
                                   !     snow-free land, 1 for snow.
   rhokh_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                   ! OUT Surface exchange coefficients
                                   !     for sea and sea-ice
   z0h_tile(land_pts,ntiles),                                           &
                                   ! OUT Tile roughness lengths for heat
                                   !     and moisture (m).
   z0m_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                   ! OUT Gridbox mean Roughness length
                                   !      for momentum (m).
   z0m_tile(land_pts,ntiles),                                           &
                                   ! OUT Tile roughness lengths for
                                   !     momentum.
  chr1p5m_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                                   ! OUT CHR1P5M for sea and sea-ice
                                   !     (leads ignored).
   g_leaf(land_pts,npft),                                               &
                                   ! OUT Leaf turnover rate (/360days).
   gpp_ft(land_pts,npft),                                               &
                                   ! OUT Gross primary productivity
                                   !     on PFTs (kg C/m2/s).
   npp_ft(land_pts,npft),                                               &
                                   ! OUT Net primary productivity
                                   !     (kg C/m2/s).
   resp_p_ft(land_pts,npft),                                            &
                                   ! OUT Plant respiration on PFTs
                                   !     (kg C/m2/s).
   resp_s(land_pts,dim_cs1),                                            &
                               ! OUT Soil respiration (kg C/m2/s).
   resp_s_tot(dim_cs2),                                                 &
                              ! OUT Total soil respiration
                              ! (kg C/m2/s).
   resp_w_ft(land_pts,npft),                                            &
                                   ! OUT Wood maintenance respiration
                                   !     (kg C/m2/s).
   gc(land_pts,ntiles),                                                 &
                                   ! OUT "Stomatal" conductance to
                                   !      evaporation for land tiles
                                   !      (m/s).
   canhc_tile(land_pts,ntiles),                                         &
                                   ! OUT Areal heat capacity of canopy
                                   !    for land tiles (J/K/m2).
   wt_ext_tile(land_pts,sm_levels,ntiles),                              &
                                   ! OUT Fraction of evapotranspiration
                                   !    which is extracted from each
                                   !    soil layer by each tile.
   flake(land_pts,ntiles),                                              &
                                   ! OUT Lake fraction.
   tile_frac(land_pts,ntiles),                                          &
                                   ! OUT Tile fractions including
                                   !     snow cover in the ice tile.
   fsmc(land_pts,npft)
                                   ! OUT Moisture availability factor.

! Additional variables for JULES
  REAL, INTENT(OUT) ::                                                  &
   rhokh_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                   ! Exchange coeffs for moisture.
   dtstar_tile(land_pts,ntiles),                                        &
                                   ! Change in TSTAR over timestep
                                   ! for land tiles
   dtstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use),&
                                   ! Change is TSTAR over timestep
                                   ! for sea-ice
   hcons(land_pts),                                                     &
                                   ! Soil thermal conductivity
                                   ! including water and ice
   emis_tile(land_pts,ntiles),                                          &
                                   ! Emissivity for land tiles
   emis_soil(land_pts)
                                   ! Emissivity of underlying soil
!-----------------------------------------------------------------------
!  Workspace :-
  INTEGER ::                                                            &
  idiv,                                                                 &
            ! loop counter, mineral dust divisions
  m !loop counter

  REAL ::                                                               &
   qw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),   &
   tl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),   &
   a_dqsdt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels),                                                  &
                                ! Saturated lapse rate factor
                                ! on p,T,q-levels (full levels).
   a_qs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels), &
                                ! Saturated lapse rate factor
                                ! on p,T,q-levels (full levels).
   bq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),   &
                                ! A buoyancy parameter for clear air
                                ! on p,T,q-levels (full levels).
  bq_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                ! A buoyancy parameter for cloudy air
                                ! on p,T,q-levels (full levels).
   bt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),   &
                                ! A buoyancy parameter for clear air
                                ! on p,T,q-levels (full levels).
  bt_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                ! A buoyancy parameter for cloudy air
                                ! on p,T,q-levels (full levels).
  deltap(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                ! Difference in pressure between levels
   dqsdt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                ! Derivative of q_SAT w.r.t. T
   dzl_charney(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,     &
               bl_levels),                                              &
                                ! DZL(,K) is depth in m of theta level
!                                 K, i.e. distance from boundary
!                                 K-1/2 to boundary K+1/2.
   fb_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                ! Surface flux buoyancy over density
                                ! (m^2/s^3)
  p_half(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                ! P_HALF(*,K) is pressure at half
                                ! level k-1/2.
   rdz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                ! RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
   rho_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels+1),                                                 &
                                ! OUT density on UV (ie. rho) levels;
                                !    used in RHOKH so dry density if
                                !    Lq_mix_bl is true
   rho_dry_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              bl_levels),                                               &
                                ! OUT density on TQ (ie. theta) levels;
                                !    used in non-turb flux integration
                                !    so dry density if Lq_mix_bl is true
   u_s(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                ! Surface friction velocity (m/s)
   recip_l_mo_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end), &
                                 ! Reciprocal of the surface
                                 ! Obukhov length over the sea
                                 ! (m-1).
  z1_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                ! Height of lowest u,v level.
  h_blend_orog(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                ! Blending height used as part of
!                                 effective roughness scheme
  rhostar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                ! Surface air density
  dust_flux(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,ndiv),  &
                                       !dust production flux(kg m-2 s-1)
  dust_flux_tile(land_pts,ntiles,ndiv),                                 &
                                            !production flux from tiles
  dust_emiss_frac(land_pts,ntiles),                                     &
                                   ! dust emiss frac on each tile
  cd_std_dust(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                   !  Bulk transfer coef. for
                              ! momentum, excluding orographic effects
  u_s_std_tile(land_pts,ntiles),                                        &
                                     ! Surface friction velocity
  clay_land(land_pts),                                                  &
                           ! soil clay fraction on land pts
  sand_land(land_pts),                                                  &
                           ! soil sand fraction on land pts
  pstar_land(land_pts),                                                 &
                            ! surface pressure on land pts
  rhostar_land(land_pts),                                               &
                              ! surface air density on land pts
  mrel_land(land_pts,ndivl),                                            &
                                ! soil size fraction on land pts
  u_s_t_tile(land_pts,ntiles,ndivh),                                    &
                                !threshold friction vel on tiles for dust
  u_s_t_dry_tile(land_pts,ntiles,ndivh),                                &
                                !dry threshold friction velocity on tiles
   rib_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                   !     Land mean bulk Richardson no.
!                                        for lowest layer.

! Dummy field for ozone concentration into JULES 
! until real field can be passed down  
  REAL :: o3_dummy(land_pts) 

  REAL ::                                                               &
  rholem,                                                               &
               ! surface density in LEM
  tv1_sd,                                                               &
               ! virt T standard deviation (approx)
  w_m     ! velocity scale

  REAL, ALLOCATABLE :: z1_uv_top(:,:)
               ! Height of top of lowest uv layer above the surface
  REAL, ALLOCATABLE :: z1_tq_top(:,:)
               ! Height of top of lowest Tq layer above the surface
  INTEGER ::                                                            &
   i,j,k,l,n

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('BDY_LAYR',zhook_in,zhook_handle)
!-----------------------------------------------------------------------
! 0. Verify grid/subset definitions.  Arakawa 'B' grid with P-rows at
!    extremes is assumed.  Extreme-most P-rows are ignored; extreme-
!    most UV-rows are used only for interpolation and are not updated.
!-----------------------------------------------------------------------
  IF ( bl_levels <  1 .OR. st_levels <  1 .OR. sm_levels <  1           &
   .OR. pdims%j_end <  1 ) THEN
    error = 1
    GO TO 9999
  ELSE IF ( land_pts >  pdims%j_end*pdims%i_end ) THEN
    error = 7
    GO TO 9999
  END IF

  IF (l_flux_bc) THEN
        ! For specified surface fluxes impose uniform TSTAR,
        ! calculated in CONV_DIAG to be consistent with fluxes
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        tstar_land(i,j) = tstar(i,j)
        tstar_sea(i,j)  = tstar(i,j)
        tstar_sice_cat(i,j,:) = tstar(i,j)
        tstar_ssi(i,j)  = tstar(i,j)
      END DO
    END DO
  END IF
!-----------------------------------------------------------------------
!     Allocate arrays for conservative discretization of the
!     surface layer.
  IF (i_modiscopt == on) THEN
    ALLOCATE(z1_uv_top(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end))
    ALLOCATE(z1_tq_top(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end))
  ELSE
    ALLOCATE(z1_uv_top(1,1))
    ALLOCATE(z1_tq_top(1,1))
  END IF

  ! these vars declared for the first time here are passed to CABLE as tl_1, qw_1
  call cable_control3( tl, qw )

! DEPENDS ON: bdy_expl1
  CALL bdy_expl1 (                                                      &
! IN values defining vertical grid of model atmosphere :
   bl_levels,                                                           &
   p,p_theta_levels, rho_rsq, rho_wet, rho_dry,                         &
! IN cloud data :
   cf,q,qcf,qcl,t,                                                      &
! IN everything not covered so far :
   pstar,lq_mix_bl,                                                     &
! OUT
   dtrdz_charney_grid,rdz_charney_grid,dtrdz_u,dtrdz_v,                 &
   rdz_u,rdz_v,rho_uv,rho_dry_tq,dzl_charney,rdz,                       &
   z1_tq,z1_tq_top,z1_uv,z1_uv_top,                                     &
   p_half,deltap,qw,tl,                                                 &
   bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt                   &
    )

! Initialise JULES dummy variable for o3 concentration
    o3_dummy(:) = 0.

! DEPENDS ON: sf_expl_l
    CALL sf_expl_l (                                                    &
! IN values defining field dimensions and subset to be processed :
     land_pts, nice_use,                                                &
! IN  parameters for iterative SISL scheme
     numcycles, cycleno,                                                &
! IN parameters required from boundary-layer scheme :
     bq_gb,bt_gb,z1_uv,z1_uv_top,z1_tq,z1_tq_top,qw,tl,                 &
! IN soil/vegetation/land surface data :
     land_index,ntiles,sm_levels,canopy,catch,catch_snow,hcon,          &
     ho2r2_orog,fland,flandg,                                           &
     snow_tile,sil_orog_land,smvccl,smvcst,smvcwt,sthf,sthu,z0_tile,    &
     z0h_tile_bare,                                                     &
! IN sea/sea-ice data :
     ice_fract_cat, k_sice,                                             &
! IN everything not covered so far :
     pstar,lw_down,sw_tile,zh,ddmfx,                                    &
     co2_mmr,co2_3d,l_co2_interactive,l_phenol,l_triffid,               &
     asteps_since_triffid,cs,frac,canht_ft,photosynth_act_rad,lai_ft,   &
     lq_mix_bl,t_soil,ti,tstar,tstar_sea,tstar_sice_cat,                &
     tstar_tile,z_land,albsoil,cos_zenith_angle,                        &
     l_aero_classic,l_dust,l_dust_diag,soil_clay,o3_dummy,              &
! IN idealised and SCM things
     l_spec_z0, z0m_scm, z0h_scm,                                       &
! IN variables for message passing
     u_1_px, v_1_px, u_0_px, v_0_px,                                    &
! IN STASH flags :-
     sfme,sq1p5,st1p5,su10,sv10,sz0heff,                                &
! INOUT data :
     l_q10,z0msea,Gs,g_leaf_acc,npp_ft_acc,resp_w_ft_acc,resp_s_acc,    &
! OUT Diagnostic not requiring STASH flags :
     recip_l_mo_sea,e_sea,fqw,ftl,ftl_tile,le_tile,h_sea,               &
     radnet_sice,radnet_tile,rhokm,rib_gb,rib_tile,                     &
! OUT variables for message passing
     flandfac, fseafac,rhokm_land, rhokm_ssi, cdr10m,                   &
! OUT diagnostic requiring STASH flags :
     fme,                                                               &
! OUT diagnostics required for soil moisture nudging scheme :
     wt_ext,ra,                                                         &
! OUT data required for tracer mixing :
     rho_aresist,aresist,resist_b,                                      &
     rho_aresist_tile,aresist_tile,resist_b_tile,                       &
     r_b_dust,cd_std_dust,u_s_std_tile,                                 &
! OUT data required elsewhere in UM system :
     fb_surf,u_s,t1_sd,q1_sd,                                           &
! OUT data required elsewhere in boundary layer or surface code
     alpha1,alpha1_sice,ashtf,ashtf_tile,fqw_tile,epot_tile,            &
     fqw_ice,ftl_ice,fraca,rhostar,resfs,resft,                         &
     rhokh,rhokh_tile,rhokh_sice,rhokh_mix,dtstar_tile,dtstar,          &
     h_blend_orog,z0hssi,z0h_tile,z0h_eff_gb,z0m_gb,z0mssi,z0m_tile,    &
     z0m_eff_gb,chr1p5m,chr1p5m_sice,smc,hcons,vshr,vshr_land,vshr_ssi, &
     gpp,npp,resp_p,g_leaf,gpp_ft,npp_ft,                               &
     resp_p_ft,resp_s,resp_s_tot,resp_w_ft,                             &
     gc,canhc_tile,wt_ext_tile,flake,                                   &
     tile_index,tile_pts,tile_frac,fsmc,emis_tile,emis_soil             &
     )

!     Store the surface friction velocity if required for the diagnosis
!     of decoupling.
  IF (IScrnTDiag == IP_ScrnDecpl2) uStarGBM = u_s

  IF (l_flux_bc) THEN
        ! Surface fluxes calculated in SFEXPL are substituted with
        ! forcing values.  NOTE: Surface calculation also made
        ! explicit (time weight set to zero).
    DO i = pdims%i_start, pdims%i_end
      DO j = pdims%j_start, pdims%j_end
         !..Converts Fluxes from W/m^2 to rho*K/s
        rholem = rhostar(i,j)

         !..If comparing against LES with rho ne rhostar then match
         ! w'theta' rather than rho*wtheta

         ! RHOLEM = 1.0
        fqw(i,j,1)   = (rhostar(i,j)*flux_e(i,j))/(lc*rholem)
        ftl(i,j,1)   = (rhostar(i,j)*flux_h(i,j))/(cp*rholem)

        fb_surf(i,j) = g * ( bt_gb(i,j,1)*ftl(i,j,1) +                  &
                             bq_gb(i,j,1)*fqw(i,j,1) ) /rhostar(i,j)
        IF ( fb_surf(i,j)  >   0.0) THEN
          w_m        = ( 0.25*zh(i,j)*fb_surf(i,j) +                    &
                         u_s(i,j)*u_s(i,j)*u_s(i,j) ) ** (1.0/3.0)
          t1_sd(i,j) = 1.93 * ftl(i,j,1) / (rhostar(i,j) * w_m)
          q1_sd(i,j) = 1.93 * fqw(i,j,1) / (rhostar(i,j) * w_m)
          tv1_sd     = t(i,j,1) * ( bt_gb(i,j,1)*t1_sd(i,j) +           &
                                    bq_gb(i,j,1)*q1_sd(i,j) )
          t1_sd(i,j) = MAX ( 0.0 , t1_sd(i,j) )
          q1_sd(i,j) = MAX ( 0.0 , q1_sd(i,j) )
          IF (tv1_sd  <=  0.0) THEN
            t1_sd(i,j) = 0.0
            q1_sd(i,j) = 0.0
          END IF
        ELSE
          t1_sd(i,j) = 0.0
          q1_sd(i,j) = 0.0
        END IF
      END DO    ! J
    END DO    ! I
    DO i = 1, land_pts
      DO l = 1, ntiles
        fqw_tile(i,l) = fqw(1,1,1)
        ftl_tile(i,l) = ftl(1,1,1)
      END DO ! L
    END DO ! I
  END IF

! NOTE: If the arguments in the subroutine "bdy_expl2" are going to
!       be changed, the same changes must be made both in
!       bdy_expl2.F90 and bdy_expl2_1a.F90.

! DEPENDS ON: bdy_expl2
  CALL bdy_expl2 (                                                      &
! IN values defining vertical grid of model atmosphere :
   bl_levels,p_theta_levels,land_pts,land_index,                        &
! IN U, V and W momentum fields.
   u_p,v_p,                                                             &
   u_0_px(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
   v_0_px(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
! IN variables for TKE scheme
 pstar,p_half,                                                          &
! IN from other part of explicit boundary layer code
   rho_uv,rho_tq,rho_dry_tq,dzl_charney,rdz,rdz_charney_grid,           &
   z_tq,z_uv,rhostar,bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt,&
   recip_l_mo_sea,                                                      &
   flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
   rib_gb, sil_orog_land,z0m_eff_gb,                                    &
! IN cloud/moisture data :
   cf,q,qcf,qcl,t,qw,tl,                                                &
! IN everything not covered so far :
   rad_hr,micro_tends,fb_surf,u_s,h_blend_orog,                         &
   lq_mix_bl, zh_prev, nlcl,zhpar,z_lcl,ho2r2_orog,sd_orog,             &
! SCM Diagnostics (dummy values in full UM) & stash diag
   nSCMDpkgs,L_SCMDiags,BL_diag,                                        &
! INOUT variables
   zh,ntpar,l_shallow,cumulus,fqw,ftl,rhokh,rhokm,w,etadot,t1_sd,q1_sd, &
! INOUT variables on TKE based turbulence schemes
   e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                        &
! OUT New variables for message passing
   tau_fd_x, tau_fd_y, f_ngstress,rhogamu, rhogamv,                     &
! OUT Diagnostic not requiring STASH flags :
   zht,shallowc,cu_over_orog,                                           &
   bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6,         &
   bl_type_7,                                                           &
! OUT data required for tracer mixing :
   ntml, kent, we_lim, t_frac, zrzi,                                    &
   kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,                          &
! OUT data required elsewhere in UM system :
   zhsc,ntdsc,nbdsc,wstar,wthvs,uw0,vw0                                 &
    )

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      rib_land(i,j)=0.0
      rib_ssi(i,j)=0.0
    END DO
  END DO

  DO n = 1, ntiles
    DO k = 1, tile_pts(n)
      l = tile_index(k,n)
      j=(land_index(l)-1)/pdims%i_end + 1
      i = land_index(l) - (j-1)*pdims%i_end
      rib_land(i,j)=rib_land(i,j) +                                     &
        rib_tile(l,n)*tile_frac(l,n)
    END DO
  END DO

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF(flandg(i,j) <  1.0)                                            &
        rib_ssi(i,j)=(rib_gb(i,j)-rib_land(i,j)*flandg(i,j))            &
          /(1.0-flandg(i,j))
    END DO
  END DO
!-----------------------------------------------------------------------
! Mineral dust production
!-----------------------------------------------------------------------
  IF (l_dust .OR. l_dust_diag) THEN
!initialisation
    dust_flux(:,:,:)=0.0
    dust_flux_tile(:,:,:) = 0.
    u_s_t_tile(:,:,:) = 0.
    u_s_t_dry_tile(:,:,:) = 0.

!put fields into land arrays
    DO l = 1, land_pts
      j = (land_index(l)-1)/pdims%i_end + 1
      i = land_index(l) - (j-1)*pdims%i_end
      pstar_land(l) = pstar(i,j)
      rhostar_land(l) = rhostar(i,j)
      sand_land(l) = soil_sand(i,j)
      clay_land(l) = soil_clay(i,j)
    END DO !LAND_PTS

    DO l = 1, land_pts
      j = (land_index(l)-1)/pdims%i_end + 1
      i = land_index(l) - (j-1)*pdims%i_end
      mrel_land(l,1) = dust_mrel1(i,j)
      mrel_land(l,2) = dust_mrel2(i,j)
      mrel_land(l,3) = dust_mrel3(i,j)
      mrel_land(l,4) = dust_mrel4(i,j)
      mrel_land(l,5) = dust_mrel5(i,j)
      mrel_land(l,6) = dust_mrel6(i,j)
    END DO !LAND_PTS

! DEPENDS ON: dust_srce
    CALL dust_srce(                                                     &
! IN arguments
         land_pts,ntiles,sm_levels,tile_pts,tile_index,fland,           &
         tstar_tile,rhostar_land,soil_layer_moisture,snow_tile,         &
         u_s_std_tile,mrel_land,clay_land,sand_land,ho2r2_orog,         &
! OUT arguments
         dust_flux_tile,u_s_t_tile,u_s_t_dry_tile                       &
         )

! Get the fraction within each tile which is bare soil, for the purpose
! of dust emission:
! DEPENDS ON: dust_calc_emiss_frac
    CALL dust_calc_emiss_frac(                                          &
  land_pts,ntiles,tile_pts,tile_index,frac,lai_ft,dust_emiss_frac       &
    )

! Produce a total dust flux over all tiles, by looping through tiles and 
! multiplying the flux on that tile by the dust_emiss_frac, and summing.
    DO idiv = 1, ndiv
      DO m = 1, ntiles
        DO n = 1, tile_pts(m)
          l = tile_index(n,m)
          j = (land_index(l)-1)/pdims%i_end + 1
          i = land_index(l) - (j-1)*pdims%i_end
          dust_flux(i,j,idiv) = dust_flux(i,j,idiv) +                   &
           dust_flux_tile(l,m,idiv)*dust_emiss_frac(l,m)
        END DO !TILE_PTS
      END DO !NTILES
    END DO !NDIV

  END IF !L_DUST .OR. L_DUST_DIAG

  DEALLOCATE(z1_uv_top)
  DEALLOCATE(z1_tq_top)

9999 CONTINUE  ! Branch for error exit.

  IF (lhook) CALL dr_hook('BDY_LAYR',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE bdy_layr
