! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module containing variables local to atmos_physics2
!
MODULE atmos_physics2_alloc_mod

USE atm_fields_bounds_mod, ONLY: pdims,pdims_s,udims,udims_s,           &
                                 vdims,vdims_s,rkmdims
USE bl_diags_mod, ONLY: strnewbldiag
!
IMPLICIT NONE
SAVE
!
! Description:
!   Contains all variables in atmos_physics2 calculated using
!   start-of-timestep quantities, and which therefore do not need
!   re-calculating at the second ENDGAME semi-lagrangian cycle
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! things from communication section at top
! for surface current interpolation 
REAL, ALLOCATABLE, TARGET :: uhalo(:,:), vhalo(:,:)
! winds on p-grid
REAL, ALLOCATABLE :: u_p(:,:,:),v_p(:,:,:)
! Land fraction on land tiles.
REAL, ALLOCATABLE :: fland(:)
! Land fraction on all points
REAL, ALLOCATABLE, TARGET :: flandg(:,:)
! extended halo winds for coastal tiling
REAL, ALLOCATABLE, TARGET :: u_1_px(:,:), v_1_px(:,:),                  &
                             u_0_px(:,:), v_0_px(:,:)
!
! things from ni_bl_ctl argument list
! out variables required in IMP_SOLVER
REAL, ALLOCATABLE :: alpha1_sice(:,:,:),ashtf(:,:,:),bq_gb(:,:,:),      &
 bt_gb(:,:,:),dtrdz_charney_grid(:,:,:),rdz_charney_grid(:,:,:),        &
 dtrdz_u(:,:,:),dtrdz_v(:,:,:),rdz_u(:,:,:),rdz_v(:,:,:),               &
 z1_tq(:,:),ustargbm(:,:)
! out diagnostics (done after implicit solver)
! e_sea, fqt, ftl, h_sea
REAL, ALLOCATABLE :: rib_gb(:,:),vshr(:,:),zht(:,:),shallowc(:,:),      &
 cu_over_orog(:,:),bl_type_1(:,:),bl_type_2(:,:),bl_type_3(:,:),        &
 bl_type_4(:,:),bl_type_5(:,:),bl_type_6(:,:),bl_type_7(:,:),           &
 z0m_eff_gb(:,:),z0h_eff_gb(:,:),fme(:,:)
! OUT diagnostics required for soil moisture nudging scheme :
REAL, ALLOCATABLE :: wt_ext(:,:),ra(:)
! out data required for tracer mixing :
REAL, ALLOCATABLE :: rho_aresist(:,:),aresist(:,:),resist_b(:,:)
!OUT variables required for mineral dust scheme
REAL, ALLOCATABLE :: r_b_dust(:,:,:),dust_flux(:,:,:),                  &
 dust_emiss_frac(:,:),u_s_t_tile(:,:,:),u_s_t_dry_tile(:,:,:),          &
 u_s_std_tile(:,:),we_lim(:,:,:),t_frac(:,:,:),zrzi(:,:,:),             &
 we_lim_dsc(:,:,:),t_frac_dsc(:,:,:),zrzi_dsc(:,:,:),zhsc(:,:)
INTEGER, ALLOCATABLE :: kent(:,:),kent_dsc(:,:)
! OUT additional variables for MOSES II
! ftl_tile, le_tile, radnet_sice, radnet_tile, fqt_tile, epot_tile
! fqt_ice, ftl_ice
REAL, ALLOCATABLE :: rib_tile(:,:),rho_aresist_tile(:,:),               &
 aresist_tile(:,:),resist_b_tile(:,:),alpha1(:,:),ashtf_tile(:,:),      &
 fraca(:,:),resfs(:,:),resft(:,:),rhokh_tile(:,:),rhokh_sice(:,:),      &
 rhokpm(:,:),rhokpm_pot(:,:),rhokpm_sice(:,:),z0hssi(:,:),z0h_tile(:,:),&
 z0m_gb(:,:),z0mssi(:,:),z0m_tile(:,:),chr1p5m(:,:),CHR1P5M_SICE(:,:),  &
 smc(:),gpp(:),npp(:),resp_p(:),g_leaf(:,:),gpp_ft(:,:),npp_ft(:,:),    &
 resp_p_ft(:,:),resp_s(:,:),resp_s_tot(:), resp_w_ft(:,:),gc(:,:),      &
 canhc_tile(:,:),wt_ext_tile(:,:,:),flake(:,:),tile_frac(:,:),fsmc(:,:),&
 cs_ch4(:),rib_ssi(:,:),vshr_land(:,:),vshr_ssi(:,:)
INTEGER, ALLOCATABLE, save :: tile_index(:,:),tile_pts(:)
! OUT additional variables for JULES
! dtstar_tile
! dtstar
REAL, ALLOCATABLE :: rhokh_mix(:,:),hcons(:),emis_tile(:,:),emis_soil(:)
! bottom block
! t1_sd, q1_sd, ntml, cumulus, nbdsc, ntdsc
! ntpar, nlcl, zhpar, zlcl, l_shallow, wstar, wthvs, delthvu
REAL, ALLOCATABLE :: uw0(:,:),vw0(:,:)
!     Declaration of new BL diagnostics.
TYPE (Strnewbldiag) :: BL_diag
! and from bdy_expl3
REAL, ALLOCATABLE :: rhokm_u(:,:,:),rhokm_v(:,:,:),cdr10m_u(:,:),       &
 cdr10m_v(:,:),flandg_u(:,:),flandg_v(:,:)
! conditional arrays, only needed with endgame saving
! for saving conv_diag
REAL, ALLOCATABLE :: conv_diag_reals(:,:,:)
INTEGER, ALLOCATABLE :: conv_diag_ints(:,:,:)
LOGICAL, ALLOCATABLE :: conv_diag_logs(:,:,:)
! for saving ni_bl_ctl
REAL, ALLOCATABLE :: bl_ctl_2d(:,:,:),bl_ctl_3d(:,:,:,:),               &
 rhokm_save(:,:,:),tile_save(:,:,:),land_save(:,:),sice_save(:,:,:,:)
INTEGER, ALLOCATABLE :: bl_ctl_int2d(:,:,:)
LOGICAL, ALLOCATABLE :: bl_ctl_log2d(:,:,:)
! for saving bdy_expl3
REAL, ALLOCATABLE :: bdy_expl3_u3d(:,:,:),bdy_expl3_v3d(:,:,:),         &
                     bdy_expl3_u2d(:,:,:), bdy_expl3_v2d(:,:,:)
!
CONTAINS
SUBROUTINE atmos_physics2_alloc(land_points,ntiles,ndiv,ndivh,npft,     &
             ntype,dim_cs1,dim_cs2,dsm_levels,bl_levels,nice_use)
!
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE switches, ONLY: IScrnTDiag
USE bl_option_mod, ONLY: l_quick_ap2
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
!
INTEGER :: land_points,bl_levels,ntiles,ndiv,ndivh,npft,ntype,dim_cs1,  &
           dim_cs2,dsm_levels,nice_use
!
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!
IF (lhook) CALL dr_hook('atmos_physics2_alloc',zhook_in,zhook_handle)
!
! things from communication section
ALLOCATE(uhalo(udims_s%i_start:udims_s%i_end,                           &
               udims_s%j_start:udims_s%j_end))
ALLOCATE(vhalo(vdims_s%i_start:vdims_s%i_end,                           &
               vdims_s%j_start:vdims_s%j_end))
ALLOCATE(u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             pdims%k_start:pdims%k_end))
ALLOCATE(v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             pdims%k_start:pdims%k_end))
ALLOCATE(fland(land_points)) ! Land fraction on land tiles.
ALLOCATE(flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end))
ALLOCATE(u_1_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end))
ALLOCATE(v_1_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end))
ALLOCATE(u_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end))
ALLOCATE(v_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end))
!
! things from ni_bl_ctl argument list
! out variables required in IMP_SOLVER
ALLOCATE(alpha1_sice(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end,nice_use))
ALLOCATE(ashtf(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,nice_use))
ALLOCATE(bq_gb(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,bl_levels))
ALLOCATE(bt_gb(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,bl_levels))
ALLOCATE(dtrdz_charney_grid(pdims%i_start:pdims%i_end,                  &
                            pdims%j_start:pdims%j_end,bl_levels))
ALLOCATE(rdz_charney_grid(pdims%i_start:pdims%i_end,                    &
                          pdims%j_start:pdims%j_end,bl_levels))
ALLOCATE(dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,   &
                 bl_levels))
ALLOCATE(dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,   &
                 bl_levels))
ALLOCATE(rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,     &
         2:bl_levels))
ALLOCATE(rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,     &
         2:bl_levels))
ALLOCATE(z1_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
IF (IScrnTDiag == IP_ScrnDecpl2) THEN
   ALLOCATE(uStarGBM(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
   ! Friction velocity for the diagnosis of decoupling
ELSE
   ALLOCATE(uStarGBM(1,1))
END IF
! out diagnostics (done after implicit solver)
ALLOCATE(rib_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
ALLOCATE(vshr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
ALLOCATE(zht(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
ALLOCATE(shallowc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! indicator of shallow cumulus
ALLOCATE(cu_over_orog(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! indicator of cumulus over steep orography
ALLOCATE(bl_type_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! stable bl indicator
ALLOCATE(bl_type_2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! decoupled scu over stable bl indicator
ALLOCATE(bl_type_3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! well mixed bl indicator
ALLOCATE(bl_type_4(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! decoupled scu over well mixed bl indicator
ALLOCATE(bl_type_5(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! decoupled scu over cumulus bl indicator
ALLOCATE(bl_type_6(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! cumulus bl indicator
ALLOCATE(bl_type_7(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! shear driven bl indicator
ALLOCATE(z0m_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
ALLOCATE(z0h_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
ALLOCATE(fme(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
! OUT diagnostics required for soil moisture nudging scheme :
ALLOCATE(wt_ext(land_points,dsm_levels)) ! cumulative fract of transp'n
ALLOCATE(ra(land_points))                ! aerodynamic resiatance (s/m)
! out data required for CLASSIC aerosol mixing and deposition :
ALLOCATE(rho_aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
ALLOCATE(aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
ALLOCATE(resist_b(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
!OUT variables required for mineral dust scheme
ALLOCATE(r_b_dust(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,ndiv))
         ! surface layer resist for dust
ALLOCATE(dust_flux(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,ndiv))
         ! dust emissions (kg m-2 s-1)
ALLOCATE(dust_emiss_frac(land_points,ntiles))
         ! OUT fraction of tile can emit dust
ALLOCATE(u_s_t_tile(land_points,ntiles,ndivh))
         !OUT threshold frict. vel
ALLOCATE(u_s_t_dry_tile(land_points,ntiles,ndivh))
         !OUT dry soil value
ALLOCATE(u_s_std_tile(land_points,ntiles))
         !OUT friction velocity
ALLOCATE(kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         !  grid-level of SML inversion
ALLOCATE(kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         !  grid-level of DSC inversion
ALLOCATE(we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3))
         !  rho*entrainment rate implied by placing of subsidence
ALLOCATE(zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3))
         !  (z-z_base)/(z_i-z_base)
ALLOCATE(t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3))
         !  a fraction of the timestep
ALLOCATE(we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3))
ALLOCATE(zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3))
ALLOCATE(t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3))
ALLOCATE(zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         !  Top of decoupled layer
! OUT additional variables for MOSES II
ALLOCATE(rib_tile(land_points,ntiles))
         ! RIB for land tiles
ALLOCATE(rho_aresist_tile(land_points,ntiles))
         ! RHOSTAR*CD_STD*VSHR on land tiles
ALLOCATE(aresist_tile(land_points,ntiles))
         ! 1/(CD_STD*VSHR) on land tiles for CLASSIC aerosol scheme
ALLOCATE(resist_b_tile(land_points,ntiles))
         ! (1/CH-1/CD_STD)/VSHR on land tiles for CLASSIC aerosol scheme
ALLOCATE(alpha1(land_points,ntiles))
         ! Mean gradient of saturated specific humidity with respect
         ! to temperature between the bottom model layer and tile surfaces
ALLOCATE(ashtf_tile(land_points,ntiles))
         !Coefficient to calculate surface heat flux into land tiles.
ALLOCATE(fraca(land_points,ntiles))
         ! Fraction of surface moisture flux with only aerodynamic
         ! resistance for snow-free land tiles.
ALLOCATE(resfs(land_points,ntiles))
         ! Combined soil, stomatal and aerodynamic resistance
         ! factor for fraction (1-FRACA) of snow-free land tiles.
ALLOCATE(resft(land_points,ntiles))
         ! Total resistance factor FRACA+(1-FRACA)*RESFS for
         !     snow-free land, 1 for snow.
ALLOCATE(rhokh_tile(land_points,ntiles))
         ! Surface exchange coefficients for land tiles
ALLOCATE(rhokh_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! Surface exchange coefficients for sea and sea-ice
ALLOCATE(rhokpm(land_points,ntiles))
         ! Surface exchange coefficient.
ALLOCATE(rhokpm_pot(land_points,ntiles))
         ! Potential evaporation exchange coefficient.
ALLOCATE(rhokpm_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! Sea-ice surface exchange coeff.
ALLOCATE(z0hssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
ALLOCATE(z0h_tile(land_points,ntiles))
          ! Tile roughness lengths for heat and moisture (m).
ALLOCATE(z0m_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! Gridbox mean roughness length for momentum (m)
ALLOCATE(z0mssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! Roughness lengths over sea (m)
ALLOCATE(z0m_tile(land_points,ntiles))
         ! Tile roughness lengths for momentum.
ALLOCATE(chr1p5m(land_points,ntiles))
         ! Ratio of coefffs for calculation of 1.5m temp for land tiles.
ALLOCATE(chr1p5m_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! CHR1P5M for sea and sea-ice (leads ignored).
ALLOCATE(smc(land_points))
ALLOCATE(gpp(land_points))
         ! Gross primary productivity (kg C/m2/s).
ALLOCATE(npp(land_points))
         ! Net primary productivity (kg C/m2/s)
ALLOCATE(resp_p(land_points))
         ! Plant respiration (kg C/m2/s)
!CABLE:kdcorbin, 11/10 - changed from NPFT
ALLOCATE(g_leaf(land_points,ntiles))
         ! Leaf turnover rate (/360days)
ALLOCATE(gpp_ft(land_points,ntiles))
         ! Gross primary productivity on PFTs (kg C/m2/s)
ALLOCATE(npp_ft(land_points,ntiles))
         ! Net primary productivity (kg C/m2/s)
ALLOCATE(resp_p_ft(land_points,ntiles))
         ! Plant respiration on PFTs (kg C/m2/s)
ALLOCATE(resp_s(land_points,dim_cs1))
         ! Soil respiration (kg C/m2/s)
ALLOCATE(resp_s_tot(dim_cs2))
         ! OUT total RESP_S over pools
ALLOCATE(resp_w_ft(land_points,npft))
         ! Wood maintenance respiration (kg C/m2/s)
ALLOCATE(gc(land_points,ntiles))
         ! "Stomatal" conductance to evaporation for land tiles (m/s)
ALLOCATE(canhc_tile(land_points,ntiles))
         ! Areal heat capacity of canopy for land tiles (J/K/m2)
ALLOCATE(wt_ext_tile(land_points,dsm_levels,ntiles))
         ! Fraction of evapotranspiration which is extracted from each
         !    soil layer by each tile.
ALLOCATE(flake(land_points,ntiles))
         ! Lake fraction
ALLOCATE(tile_index(land_points,ntype))
ALLOCATE(tile_pts(ntype))
ALLOCATE(tile_frac(land_points,ntiles))
         ! Tile fractions including snow cover in the ice tile
ALLOCATE(fsmc(land_points,npft))
         ! Soil moisture availability factor
ALLOCATE(cs_ch4(land_points))
         ! Effective soil carbon for CH4 wetlands calc (kg C/m2)
ALLOCATE(rib_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! Rib over sea part of grid box
ALLOCATE(vshr_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! Vshr over land part of grid box
ALLOCATE(vshr_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! Vshr over sea part of grid box
! OUT additional variables for JULES
ALLOCATE(rhokh_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
ALLOCATE(hcons(land_points))
ALLOCATE(emis_tile(land_points,ntiles))
ALLOCATE(emis_soil(land_points))
! bottom block
ALLOCATE(uw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! U-component of surface wind stress (P-grid)
ALLOCATE(vw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
         ! V-component of surface wind stress (P-grid)
! and from bdy_expl3
ALLOCATE(cdr10m_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end))
ALLOCATE(cdr10m_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end))
ALLOCATE(flandg_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end))
         ! Land frac (on U-grid)
ALLOCATE(flandg_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end))
         ! Land frac (on V-grid)
ALLOCATE(rhokm_u(udims%i_start:udims%i_end,                             &
                 udims%j_start:udims%j_end, bl_levels))
ALLOCATE(rhokm_v(vdims%i_start:vdims%i_end,                             &
                 vdims%j_start:vdims%j_end, bl_levels))
!
! if we are optimising the 2nd ENDGAME cycle, then we need to allocate
! these arrays to save output
IF (l_quick_ap2) THEN
! save conv_diag output
   ALLOCATE(conv_diag_reals(pdims%i_start:pdims%i_end,                  &
                            pdims%j_start:pdims%j_end,13)) ! reals
   ALLOCATE(conv_diag_ints(pdims%i_start:pdims%i_end,                   &
                           pdims%j_start:pdims%j_end,4))   ! integers
   ALLOCATE(conv_diag_logs(pdims%i_start:pdims%i_end,                   &
                           pdims%j_start:pdims%j_end,4))   ! logicals
! save bl_ctl output
   ALLOCATE(bl_ctl_2d(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end,6))        ! 2d reals
   ALLOCATE(bl_ctl_int2d(pdims%i_start:pdims%i_end,                     &
                         pdims%j_start:pdims%j_end,1))     ! 2d integers
   ALLOCATE(bl_ctl_log2d(pdims%i_start:pdims%i_end,                     &
                         pdims%j_start:pdims%j_end,2))     ! 2d logicals
   ALLOCATE(bl_ctl_3d(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end,bl_levels,3)) !3d reals
   ALLOCATE(rhokm_save(rkmdims%i_start:rkmdims%i_end,                   &
                       rkmdims%j_start:rkmdims%j_end,                   &
                       rkmdims%k_start:rkmdims%k_end))
   ALLOCATE(tile_save(land_points,ntiles,6))        ! land reals on tiles
   ALLOCATE(land_save(land_points,1))               ! land reals
   ALLOCATE(sice_save(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end, nice_use, 4)) ! sea-ice fields
! save bdy_expl3 output
   ALLOCATE(bdy_expl3_u3d(udims%i_start:udims%i_end,                    &
                          udims%j_start:udims%j_end,bl_levels))
            ! 3d reals on u-points
   ALLOCATE(bdy_expl3_v3d(vdims%i_start:vdims%i_end,                    &
                          vdims%j_start:vdims%j_end,bl_levels))
            ! 3d reals on v-points
   ALLOCATE(bdy_expl3_u2d(udims%i_start:udims%i_end,                    &
                          udims%j_start:udims%j_end,2))
            ! 2d reals on u-points
   ALLOCATE(bdy_expl3_v2d(vdims%i_start:vdims%i_end,                    &
                          vdims%j_start:vdims%j_end,2))
            ! 2d reals on v-points
END IF !l_quick_ap2
!
IF (lhook) CALL dr_hook('atmos_physics2_alloc',zhook_out,zhook_handle)
!
RETURN
END SUBROUTINE atmos_physics2_alloc
END MODULE atmos_physics2_alloc_mod
