! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!!  SUBROUTINE IMP_SOLVER---------------------------------------------
!!!
!!!  Purpose: implicit solver for diffusion equation
!!!           split from bdy_layr routine
!!!
!!!
!!! Programming standard : unified model documentation paper No 3
!!!
!!!  Documentation: UMDP 24.
!!!
!!!---------------------------------------------------------------------

!    Arguments :-
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Boundary Layer
SUBROUTINE imp_solver (                                                 &

! IN MPP variables
 halo_i, halo_j, offx , offy , row_length, rows, n_rows,                &
 global_row_length, proc_row_group, at_extremity,                       &
 n_proc, n_procx, n_procy, neighbour,                                   &

! IN values defining field dimensions and subset to be processed :
 ntiles, land_pts, nice, nice_use,                                      &

! IN values defining vertical grid of model atmosphere :
 model_domain,                                                          &
 bl_levels,                                                             &
 GAMMA,                                                                 &

! IN U and V momentum fields.
 u, v,                                                                  &
! IN Non turbulent increments to momentum
!  (New dynamics only).
 du_nt,dv_nt,                                                           &

! IN soil/vegetation/land surface data :
 land_mask,land_index,                                                  &
 st_levels,sm_levels,tile_frac,canopy,                                  &
 fland,flandg,                                                          &

! IN sea/sea-ice data :
 di, ice_fract, di_ncat, ice_fract_ncat, k_sice, u_0, v_0,              &

! IN cloud data :
 q,qcf,qcl,qcf_latest,qcl_latest, t,                                    &

! IN everything not covered so far :
 rho_wet_theta,pstar,timestep,l_sice_heatflux,                          &
 l_sice_multilayers,                                                    &

! IN variables from BDY_LAYR (that used to be local arrays)
 alpha1,ashtf,bq_gb,bt_gb,dtrdz_charney_grid,rdz_charney_grid,          &
 dtrdz_u,dtrdz_v,rdz_u,rdz_v,                                           &
 fraca,rhokh_tile,smc,                                                  &
 chr1p5m,resfs,z0hssi,z0mssi,cdr10m_u,cdr10m_v,z1_tq,zh,                &
!ajm extra variable added
 rhokm_u,rhokm_v,                                                       &
! needed for new BL numerical solver
 bl_type_1,bl_type_2,bl_type_3,bl_type_4,                               &
 bl_type_5,bl_type_6,bl_type_7,                                         &

! IN additional variables for MOSES II
 tile_pts,tile_index,                                                   &
 canhc_tile,flake,wt_ext_tile,                                          &
 lw_down,sw_tile,alpha1_sice,ashtf_tile,                                &
 fqw_ice,ftl_ice,resft,rhokh_sice,rhokpm,rhokpm_pot,rhokpm_sice,        &
 z0h_tile,z0m_tile,chr1p5m_sice,                                        &
 flandg_u,flandg_v,                                                     &

! IN additional variables used to calculate atmospheric cooling
!    at the screen level
 l_co2_interactive, co2_mmr, co2_3d,                                    &
 rho1, f3_at_p, uStarGBM,                                               &

! IN additional variables for JULES
 rhokh_mix_dummy,dtstar_tile,dtstar,hcons,emis_tile,emis_soil,          &

! INOUT data :
 t_soil,ti,ti_gb,tstar,                                                 &
 tstar_land,tstar_sea,tstar_sice,tstar_sice_cat,tstar_ssi,              &
 tstar_tile, t_latest,q_latest,                                         &

! INOUT Diagnostics started in BDY_LAYR not requiring STASH flags :
 e_sea,fqw,fqw_tile,epot_tile,ftl,ftl_tile,h_sea,rhokh,                 &
 taux,tauy,                                                             &
 taux_land,taux_ssi,tauy_land,tauy_ssi,                                 &
 TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans,                                &

! INOUT additional variables for MOSES II
 snow_tile,le_tile,radnet_sice,radnet_tile,olr,                         &

! OUT Increments to U and V momentum fields and Tl qw
!  (New dynamics only).
 du,dv,                                                                 &

! OUT Diagnostic not requiring STASH flags :
 rhokh_mix,sea_ice_htf,                                                 &

! OUT diagnostics requiring STASH flags :
 sice_mlt_htf,snomlt_surf_htf,latent_heat,                              &
 q1p5m,t1p5m,u10m,v10m,                                                 &

! (IN) STASH flags :-
 simlt,smlt,slh,sq1p5,st1p5,su10,sv10,                                  &
 l_ftl, l_fqw, l_taux, l_tauy,                                          &

! OUT additional variables for MOSES II
 esoil_tile,surf_ht_flux,surf_ht_flux_land,surf_ht_flux_sice,           &
 surf_htf_tile,ei_tile,                                                 &
 q1p5m_tile,t1p5m_tile,ecan_tile,melt_tile,                             &
 e_ssi,ei_sice,ftl_ssi,                                                 &
 error,                                                                 &

! OUT data required elsewhere in UM system :
 ecan,ei,es,ext,snowmelt,                                               &
 BL_diag,lq_mix_bl,                                                     &
! IN SCM LOGICAL
 l_flux_bc,                                                             &
 ltimer                                                                 &
 )

  USE dynamics_grid_mod, ONLY: l_vatpoles

  USE nstypes
  USE atm_fields_bounds_mod, ONLY:                                      &
   udims, vdims, udims_s, vdims_s, udims_l, vdims_l, pdims

  USE atmos_constants_mod, ONLY: cp
  
  USE level_heights_mod, ONLY:                                          &
    r_theta_levels, r_rho_levels, r_at_u, r_at_v

  USE bl_option_mod, ONLY:                                              &
      fric_heating, l_us_blsol, puns,pstb, on

  USE bl_diags_mod, ONLY :                                              &
      strnewbldiag

  USE swapable_field_mod, ONLY :                                        &
      swapable_field_pointer_type

  USE earth_constants_mod, ONLY:g 

  USE water_constants_mod, ONLY: lc, lf
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  USE Field_Types
  USE domain_params
  USE u_to_p_mod, ONLY: u_to_p
  USE v_to_p_mod, ONLY: v_to_p
  IMPLICIT NONE

!  Inputs :-

! (a) Defining horizontal grid and subset thereof to be processed.

  INTEGER ::                                                            &
! cjj additions - MPP variables.
    row_length,                                                         &
                   ! Local number of points on a row
    rows,                                                               &
                   ! Local number of rows in a theta field
    n_rows,                                                             &
                   ! Local number of rows in a v field
    halo_i,                                                             &
                   ! Size of halo in i direction.
    halo_j,                                                             &
                   ! Size of halo in j direction.
    offx ,                                                              &
                   ! Size of small halo in i
    offy ,                                                              &
                   ! Size of small halo in j.
    global_row_length,                                                  &
                           ! number of points on a row
    proc_row_group,                                                     &
                       ! Group id for processors on the same row
    n_proc,                                                             &
                   ! Total number of processors
    n_procx,                                                            &
                   ! Number of processors in longitude
    n_procy,                                                            &
                   ! Number of processors in latitude
    neighbour(4)   ! Array with the Ids of the four neighbours in
                       !   the horizontal plane

  LOGICAL ::                                                            &
    at_extremity(4)  ! Indicates if this processor is at north,
                         !   south, east or west of the processor grid

! Parameters
  INTEGER ::                                                            &
     PNorth,                                                            &
                      ! North processor address in the neighbor array
     PEast,                                                             &
                      ! East processor address in the neighbor array
     PSouth,                                                            &
                      ! South processor address in the neighbor array
     PWest,                                                             &
                      ! West processor address in the neighbor array
     NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
  PARAMETER (                                                           &
     PNorth   = 1,                                                      &
     PEast    = 2,                                                      &
     PSouth   = 3,                                                      &
     PWest    = 4,                                                      &
     NoDomain = -1)

  INTEGER ::                                                            &
   ntiles,                                                              &
                                 ! IN number of land tiles
   land_pts,                                                            &
                                 ! IN No.of land points in whole grid.
   nice,                                                                &
                                 ! IN No.of sea ice categories
   nice_use                      ! IN No.of sea ice categories used fully
                                 !    in surface exchange

! (b) Defining vertical grid of model atmosphere.

  INTEGER ::                                                            &
   bl_levels,                                                           &
                                   ! IN Max. no. of "boundary" levels
!                                     allowed.Assumed <= 30 for dim-
!                                     sioning of GAMMA in common deck
!                                     C_GAMMA used in SF_EXCH and KMKH
   model_domain

  REAL ::                                                               &
    bl_type_1(row_length,rows),                                         &
                                  ! IN Indicator set to 1.0 if stable
                                  !     b.l. diagnosed, 0.0 otherwise.
    bl_type_2(row_length,rows),                                         &
                                  ! IN Indicator set to 1.0 if Sc over
                                  !     stable surface layer diagnosed,
                                  !     0.0 otherwise.
    bl_type_3(row_length,rows),                                         &
                                  ! IN Indicator set to 1.0 if well
                                  !     mixed b.l. diagnosed,
                                  !     0.0 otherwise.
    bl_type_4(row_length,rows),                                         &
                                  ! IN Indicator set to 1.0 if
                                  !     decoupled Sc layer (not over
                                  !     cumulus) diagnosed,
                                  !     0.0 otherwise.
    bl_type_5(row_length,rows),                                         &
                                  ! IN Indicator set to 1.0 if
                                  !     decoupled Sc layer over cumulus
                                  !     diagnosed, 0.0 otherwise.
    bl_type_6(row_length,rows),                                         &
                                  ! IN Indicator set to 1.0 if a
                                  !     cumulus capped b.l. diagnosed,
                                  !     0.0 otherwise.
    bl_type_7(row_length,rows)
                                  ! IN Indicator set to 1.0 if a
                                  !     shear-dominated b.l. diagnosed,
                                  !     0.0 otherwise.

  LOGICAL ::                                                            &
   lq_mix_bl              ! TRUE if mixing ratios used in
                              ! boundary layer code

  REAL ::                                                               &
    GAMMA(bl_levels),                                                   &
    du_nt(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,  &
          bl_levels),                                                   &
                            ! non-turbulent increment to u wind field
    dv_nt(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,  &
          bl_levels),                                                   &
                            ! non-turbulen increment to v wind field
    u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,      &
        bl_levels),                                                     &
    v(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,      &
        bl_levels)


! (c) Soil/vegetation/land surface parameters (mostly constant).

  LOGICAL ::                                                            &
   land_mask(row_length,rows)  ! IN T if land, F elsewhere.



  INTEGER ::                                                            &
   land_index(land_pts)      ! IN LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.

  INTEGER ::                                                            &
   st_levels,                                                           &
                                   ! IN No. of deep soil temp. levels
   sm_levels                   ! IN No. of soil moisture levels

  REAL ::                                                               &
   tile_frac(land_pts,ntiles),                                          &
                                ! IN fractional coverage for each
                                !    surface tile
   canopy(land_pts,ntiles),                                             &
                                ! IN Surface/canopy water (kg/m2)
   fland(land_pts),                                                     &
                                ! IN Land fraction on land tiles.
   flandg(row_length,rows)
                                ! IN Land fraction on all tiles.

! (d) Sea/sea-ice data.

  REAL ::                                                               &
   di(row_length,rows),                                                 &
                                   ! IN "Equivalent thickness" of
                                   !   sea-ice(m).
   ice_fract(row_length,rows),                                          &
                                   ! IN Fraction of gridbox covered by
                                   !  sea-ice (decimal fraction).
   di_ncat(row_length,rows,nice),                                       &
                                   ! IN "Equivalent thickness" of
                                   !   sea-ice on categories(m).
   ice_fract_ncat(row_length,rows,nice),                                &
                                   ! IN Fraction of gridbox
                                   !   covered by sea-ice category.
   k_sice(row_length, rows, nice),                                      &
                                   ! IN sea ice effective conductivity in
                                   !  sfc layer on categories (W/m2/K)
   u_0(udims%i_start:udims%i_end,udims%j_start:udims%j_end),            &
                                   ! IN W'ly component of surface
                                   !  current (m/s).
   v_0(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                   ! IN S'ly component of surface
                                   !  current (m/s).

! (e) Cloud data.

  REAL ::                                                               &
   qcf(row_length,rows,bl_levels),                                      &
                                       ! IN Cloud ice (kg per kg air)
   qcl(row_length,rows,bl_levels),                                      &
                                       ! IN Cloud liquid water
   q(1,row_length,rows,bl_levels),                                      &
                                       ! IN specific humidity
   t(row_length,rows,bl_levels),                                        &
                                       ! IN temperature
                          ! Latest estimates to time level n+1 values
   qcf_latest(row_length,rows,bl_levels),                               &
                                ! IN Cloud ice (kg per kg air)
   qcl_latest(row_length,rows,bl_levels) ! IN Cloud liquid water

! (f) Atmospheric + any other data not covered so far, incl control.

  REAL ::                                                               &
   pstar(row_length,rows),                                              &
                                   ! IN Surface pressure (Pascals).
   rho_wet_theta(row_length,rows,bl_levels),                            &
                                   ! IN wet density on theta levels
   timestep                    ! IN Timestep (seconds).

  LOGICAL ::                                                            &
   l_sice_heatflux,                                                     &
                               ! IN T: semi-implicit sea-ice temp
   l_sice_multilayers 
                               ! IN T: coupled to multilayers sea ice model


! IN variables from BDY_LAYR (that used to be local arrays)
  REAL ::                                                               &
   alpha1(land_pts,ntiles),                                             &
                               ! IN Mean gradient of saturated
!                                 specific humidity with
!                                 respect to temperature between
!                                 the bottom model layer and the
!                                 tile surfaces.
   ashtf(row_length,rows,nice_use),                                     &
                                  ! IN Coefficient to calculate surface
!                                 heat flux into sea-ice.
   dtrdz_charney_grid(row_length,rows,bl_levels),                       &
                                ! IN -g.dt/dp for model layers.
   rdz_charney_grid(row_length,rows,bl_levels),                         &
                                ! IN RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
   dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                          ! IN
   dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                           ! IN
!                                 -g.dt/dp for model wind layers.
   bq_gb(row_length,rows,bl_levels),                                    &
                                ! IN grid-box mean buoyancy parameter
                                !     on p,T,q-levels (full levels).
   bt_gb(row_length,rows,bl_levels),                                    &
                                ! IN grid-box mean buoyancy parameter
                                !     on p,T,q-levels (full levels).
   rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
          2:bl_levels),                                                 &
                                ! IN 1/(Z_U(K)-Z_U(K-1)) m^{-1}
   rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
          2:bl_levels),                                                 &
                                ! IN 1/(Z_V(K)-Z_V(K-1)) m^{-1}
   fraca(land_pts,ntiles),                                              &
                                ! IN Fraction of surface
                                !            moisture flux with only
                                !            aerodynamic resistance.
   rhokh_tile(land_pts,ntiles),                                         &
                                   ! IN
!                                 Tile surface exchange coefficients
!                                 for heat
   smc(land_pts),                                                       &
                               ! IN Soil moisture content in root depth
!                                  (kg/m2).
   chr1p5m(land_pts,ntiles),                                            &
                                ! IN Ratio of coefficients reqd for
!                                 calculation of 1.5 m T.
   resfs(land_pts,ntiles),                                              &
                              ! IN Combined soil, stomatal
!                                 and aerodynamicresistance
!                                 factor = PSIS/(1+RS/RA) for
!                                 fraction (1-FRACA)
   z0hssi(row_length,rows),                                             &
   z0mssi(row_length,rows),                                             &
                                ! IN Roughness lengths over sea
   cdr10m_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
                                ! IN Ratio of CD's reqd for calculation
   cdr10m_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),       &
                                ! IN Ratio of CD's reqd for calculation
   z1_tq(row_length,rows),                                              &
                                ! IN Height of lowest theta level.
   zh(row_length,rows),                                                 &
                                ! IN BL depth (m)
   rhokm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
            bl_levels),                                                 &
                                ! IN Exchange coefficients for u
   rhokm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
            bl_levels)
                                ! IN Exchange coefficients for v

! IN Additional variables for screen-level diagnostics
  LOGICAL, INTENT(IN) :: l_co2_interactive
                                ! Flag for interactive 3-D CO2
  REAL, INTENT(IN)    :: co2_mmr
                                ! Initial or fixed mass mixing ratio
                                ! of CO2
  REAL, INTENT(IN)    :: co2_3d                                         &
                      (1-offx :row_length+offx , 1-offy :rows+offy )
                                ! 3-D field of CO2
  REAL, INTENT(IN)    :: rho1(row_length, rows)
                                ! Density on lowest level
  REAL, INTENT(IN)    :: f3_at_p(row_length, rows)
                                ! Coriolis parameter
  REAL, INTENT(IN)    :: uStarGBM(row_length,rows)
                                ! GBM surface friction velocity



! IN additional variables for MOSES II

  INTEGER ::                                                            &
   tile_pts(ntype),                                                     &
                                 ! IN Number of tile points.
   tile_index(land_pts,ntype)
                                 ! IN Index of tile points.

  REAL ::                                                               &
   canhc_tile(land_pts,ntiles),                                         &
                                ! IN Areal heat capacity of canopy
                                !    for land tiles (J/K/m2).
   flake(land_pts,ntiles),                                              &
                                ! IN Lake fraction.
   wt_ext_tile(land_pts,sm_levels,ntiles),                              &
                                ! IN Fraction of evapotranspiration
                                !    which is extracted from each
                                !    soil layer by each tile.
   lw_down(row_length,rows),                                            &
                                ! IN Surface downward LW radiation
                                !    (W/m2).
   sw_tile(land_pts,ntiles),                                            &
                                ! IN Surface net SW radiation on land
                                !    tiles (W/m2).
   alpha1_sice(row_length,rows,nice_use),                               &
                                ! IN ALPHA1 for sea-ice.
   ashtf_tile(land_pts,ntiles),                                         &
                                ! IN Coefficient to calculate
                                !    surface heat flux into land
                                !    tiles.
   fqw_ice(row_length,rows,nice_use),                                   &
                                ! IN Surface FQW for sea-ice
   ftl_ice(row_length,rows,nice_use),                                   &
                                ! IN Surface FTL for sea-ice
   resft(land_pts,ntiles),                                              &
                                ! IN Total resistance factor.
                                !    FRACA+(1-FRACA)*RESFS for
                                !    snow-free land, 1 for snow.
   rhokh_sice(row_length,rows),                                         &
                                ! IN Surface exchange coefficients
                                !    for sea and sea-ice
   rhokpm(land_pts,ntiles),                                             &
                                ! IN Land surface exchange coeff.
   rhokpm_pot(land_pts,ntiles),                                         &
                                ! IN Land surface exchange coeff.
!                                    for potential evaporation.
   rhokpm_sice(row_length,rows),                                        &
                                ! IN Sea-ice surface exchange coeff.
   z0h_tile(land_pts,ntiles),                                           &
                                ! IN Tile roughness lengths for heat
                                !    and moisture (m).
   z0m_tile(land_pts,ntiles),                                           &
                                ! IN Tile roughness lengths for
                                !    momentum.
   chr1p5m_sice(row_length,rows),                                       &
                                ! IN CHR1P5M for sea and sea-ice
                                !    (leads ignored).
   flandg_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
                                ! IN Land frac (on U-grid, with 1st
                                !    and last rows undefined or, at
                                !    present, set to "missing data")
   flandg_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                ! IN Land frac (on V-grid, with 1st
                                !    and last rows undefined or, at
                                !    present, set to "missing data")


! Additional variables for JULES

  REAL ::                                                               &
   rhokh_mix_dummy(row_length,rows),                                    &
                                   ! IN  Exchange coeffs for moisture.
   dtstar_tile(land_pts,ntiles),                                        &
                                   ! IN  Change in TSTAR over timestep
                                   !     for land tiles
   dtstar(row_length,rows,nice_use),                                    &
                                   ! IN  Change is TSTAR over timestep
                                   !     for sea-ice
   hcons(land_pts),                                                     &
                                   ! IN  Soil thermal conductivity
                                   !     including water and ice
   emis_tile(land_pts,ntiles),                                          &
                                   ! IN  Emissivity for land tiles
   emis_soil(land_pts)
                                   ! IN  Emissivity of underlying soil

!  STASH flags :-

  LOGICAL ::                                                            &
   simlt,                                                               &
               ! IN Flag for SICE_MLT_HTF (q.v.)
   smlt,                                                                &
               ! IN Flag for SNOMLT_SURF_HTF (q.v.)
   slh,                                                                 &
               ! IN Flag for LATENT_HEAT (q.v.)
   sq1p5,                                                               &
               ! IN Flag for Q1P5M (q.v.)
   st1p5,                                                               &
               ! IN Flag for T1P5M (q.v.)
   su10,                                                                &
               ! IN Flag for U10M (q.v.)
   sv10,                                                                &
               ! IN Flag for V10M (q.v.)
   l_ftl,                                                               &
   l_fqw,                                                               &
   l_taux,                                                              &
   l_tauy

!  In/outs :-

  REAL ::                                                               &
   q_latest(row_length,rows,bl_levels),                                 &
                                            ! IN specific humidity
   t_latest(row_length,rows,bl_levels),                                 &
                                            ! IN temperature
   t_soil(land_pts,sm_levels),                                          &
                                  ! INOUT Soil temperatures (K).
   ti(row_length,rows,nice),                                            &
                                  ! INOUT Sea-ice category surface layer
                                  !    temperature (K)
                                  ! (IN only if l_sice_multilayers=T)
   ti_gb(row_length,rows),                                              &
                                  ! INOUT Sea-ice sfc layer temperature 
                                  !   (ice mean) (K)
                                  ! (IN only if l_sice_multilayers=T)
   tstar(row_length,rows),                                              &
                                  ! INOUT Surface temperature (K).
   tstar_land(row_length,rows),                                         &
                                  ! OUT   Land mean sfc temperature (K)
   tstar_sea(row_length,rows),                                          &
                                  ! IN    Open sea sfc temperature (K).
   tstar_sice(row_length,rows),                                         &
                                  ! OUT Sea-ice sfc temperature (K)
                                  ! (ice mean over categories)
   tstar_sice_cat(row_length,rows,nice_use),                            &
                                  ! INOUT Sea-ice sfc temperature (K).
   tstar_ssi(row_length,rows),                                          &
                                  ! INOUT Sea mean sfc temperature (K).
   tstar_tile(land_pts,ntiles)
                                  ! INOUT Surface tile temperature


!  In/Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-

  REAL ::                                                               &
   e_sea(row_length,rows),                                              &
                                ! INOUT Evaporation from sea times
!                                     leads fraction. Zero over land.
!                                     (kg per square metre per sec).
   fqw(row_length,rows,bl_levels),                                      &
                                ! INOUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
   fqw_tile(land_pts,ntiles),                                           &
                                ! INOUT surface tile moisture flux
   epot_tile(land_pts,ntiles),                                          &
                                ! INOUT surface tile potential
                                !       evaporation
   e_ssi(row_length,rows),                                              &
                                ! OUT   Surface FQW for mean sea.
   ei_sice(row_length,rows,nice_use),                                   &
                                ! OUT   Sea-ice sumblimation
                                !       (sea mean).
   ftl(row_length,rows,bl_levels),                                      &
                                ! INOUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
   ftl_tile(land_pts,ntiles),                                           &
                                ! INOUT surface tile heat flux
   ftl_ssi(row_length,rows),                                            &
                                ! OUT sea mean surface heat flux
   h_sea(row_length,rows),                                              &
                                ! INOUT Surface sensible heat flux over
!                                     sea times leads fraction. (W/m2)
   rhokh(row_length,rows,bl_levels),                                    &
                                ! INOUT Exchange coeffs for moisture.
   taux(udims%i_start:udims%i_end,udims%j_start:udims%j_end,            &
         bl_levels),                                                    &
                                ! INOUT W'ly component of surface wind
!                                     stress (N/sq m).(On UV-grid with
!                                     first and last rows undefined or
!                                     at present, set to missing data
   taux_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end),      &
                                 ! INOUT W'ly component of land sfc wind
                                 !     stress (N/sq m). (On U-grid
                                 !     with first and last rows
                                 !     undefined or, at present,
                                 !     set to missing data
   taux_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
                                 ! INOUT W'ly compt of mean sea sfc wind
                                 !     stress (N/sq m). (On U-grid
                                 !     with first and last rows
                                 !     undefined or, at present,
                                 !     set to missing data
   tauy(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,            &
         bl_levels),                                                    &
                                ! INOUT S'ly component of surface wind
!                                     stress (N/sq m).  On UV-grid;
!                                     comments as per TAUX.
   tauy_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),      &
                                 ! INOUT S'ly component of land sfc wind
                                 !     stress (N/sq m).  On V-grid;
                                 !     comments as per TAUX.
   tauy_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                 ! INOUT S'ly compt of mean sea sfc wind
                                 !     stress (N/sq m).  On V-grid;
                                 !     comments as per TAUX.

  REAL, INTENT(INOUT) :: TScrnDcl_SSI(row_length,rows)
                            !    Decoupled screen-level temperature
                            !    over sea or sea-ice
  REAL, INTENT(INOUT) :: TScrnDcl_TILE(land_pts,ntiles)
                            !    Decoupled screen-level temperature
                            !    over land tiles
  REAL, INTENT(INOUT) :: tStbTrans(row_length,rows)
                            !    Time since transition to stable
                            !    conditions


! INOUT additional variables for MOSES II
  REAL ::                                                               &
   snow_tile(land_pts,ntiles),                                          &
                                ! INOUT Snow on tiles (kg/m2).
   le_tile(land_pts,ntiles),                                            &
                                ! INOUT Surface latent heat flux for
                                !       land tiles (W/m2).
   radnet_sice(row_length,rows,nice_use),                               &
                                ! INOUT Sea-ice surface net radiation.
   radnet_tile(land_pts,ntiles),                                        &
                                ! INOUT Tile surface net radiation.
   olr(row_length,rows)     ! IN    TOA - surface upward LW on
                                !       last radiation timestep
                                ! OUT   Corrected TOA outward LW


!  Outputs :-

! (New dynamics only)
  REAL ::                                                               &
   du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,      &
        bl_levels),                                                     &
                                ! OUT BL increment to u wind field
   dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,      &
        bl_levels)          ! OUT BL increment to u wind field


!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-

  REAL ::                                                               &
   rhokh_mix(row_length,rows,bl_levels),                                &
                                ! OUT Exchange coeffs for moisture.
! for use in tracer mixing routines
   sea_ice_htf(row_length,rows,nice)
                                ! OUT Heat flux through sea-ice
!                                     (W/m2, positive downwards).

!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-

  REAL ::                                                               &
   sice_mlt_htf(row_length,rows,nice),                                  &
                                ! OUT Heat flux due to melting of sea-
!                                   ice (Watts per sq metre).
   snomlt_surf_htf(row_length,rows),                                    &
                                ! OUT Heat flux required for surface
!                                   melting of snow (W/m2).
   latent_heat(row_length,rows),                                        &
                                    ! OUT Surface latent heat flux, +ve
!                                   upwards (Watts per sq m).
   q1p5m(row_length,rows),                                              &
                                ! OUT Q at 1.5 m (kg water per kg air).
   t1p5m(row_length,rows),                                              &
                                ! OUT T at 1.5 m (K).
   u10m(udims%i_start:udims%i_end,udims%j_start:udims%j_end),           &
                                ! OUT U at 10 m (m per s).
   v10m(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)  
                                ! OUT V at 10 m (m per s).

!-2 Genuinely output, needed by other atmospheric routines :-

  REAL ::                                                               &
   ei(row_length,rows),                                                 &
                                ! OUT Sublimation from lying snow or
!                                   sea-ice (kg/m2/s).
   ecan(row_length,rows),                                               &
                                ! OUT Gridbox mean evaporation from
!                                   canopy/surface store (kg/m2/s).
!                                   Zero over sea.
   es(row_length,rows),                                                 &
                                ! OUT Surface evapotranspiration from
                                !     soil moisture store (kg/m2/s).
   ext(land_pts,sm_levels),                                             &
                                ! OUT Extraction of water from each
                                !     soil layer (kg/m2/s).
   snowmelt(row_length,rows)! OUT Snowmelt (kg/m2/s).


! OUT additional variables for MOSES II
  REAL ::                                                               &
   esoil_tile(land_pts,ntiles),                                         &
                                ! OUT Evaporation from bare soil (kg/m2)
   surf_ht_flux(row_length,rows),                                       &
                                ! OUT Net downward heat flux at surface
!                                     over land and sea-ice fraction of
!                                     gridbox (W/m2).
   surf_ht_flux_land(row_length,rows),                                  &
                                ! OUT Net downward heat flux at
                                !     surface over land
                                !     fraction of gridbox (W/m2).
   surf_ht_flux_sice(row_length,rows,nice),                             &
                                ! OUT Net category downward heat flux at
                                !     surface over sea-ice
                                !     fraction of gridbox (W/m2).
   surf_htf_tile(land_pts,ntiles),                                      &
                                ! OUT Net downward surface heat flux
                                !     on tiles (W/m2)
   ei_tile(land_pts,ntiles),                                            &
                                ! OUT EI for land tiles
   q1p5m_tile(land_pts,ntiles),                                         &
                                ! OUT Q1P5M over land tiles.
   t1p5m_tile(land_pts,ntiles),                                         &
                                ! OUT T1P5M over land tiles.
   ecan_tile(land_pts,ntiles),                                          &
                                ! OUT ECAN for land tiles
   melt_tile(land_pts,ntiles)
                                ! OUT Snowmelt on tiles (kg/m2/s).

  INTEGER ::                                                            &
   error                    ! OUT 0 - AOK;
                                ! 1 to 7  - bad grid definition detected;


  LOGICAL :: ltimer               ! Logical switch for TIMER diags
  LOGICAL :: l_flux_bc            ! Logical for prescribing surface
                                   ! flux forcing in the SCM

!     Declaration of new BL diagnostics.
  TYPE (strnewbldiag) :: BL_diag

!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-


! Derived local parameters.

  REAL :: lcrcp,ls,lsrcp
  REAL :: pnonl,p1,p2,i1,e1,e2 ! parameters for new BL solver
  REAL :: sqrt2                ! SQRT(2.)

  PARAMETER (                                                           &
   lcrcp=lc/cp,                                                         &
                             ! Evaporation-to-dT conversion factor.
   ls=lf+lc,                                                            &
                             ! Latent heat of sublimation.
   lsrcp=ls/cp                                                          &
                             ! Sublimation-to-dT conversion factor.
    )

  REAL ::                                                               &
    gamma1(row_length,rows),                                            &
    gamma2(row_length,rows)

!-----------------------------------------------------------------------

!  Workspace :-

  INTEGER ::                                                            &
   nblyr(row_length,rows)           ! number of levels in boundary
                                        !  layer, for frictional heating

  REAL ::                                                               &
   qw(row_length, rows, bl_levels),                                     &
                                        ! LOCAL total water
   tl(row_length, rows, bl_levels),                                     &
                                        ! LOCAL liquid water temperature
   dqw(row_length,rows,bl_levels),                                      &
                                        ! LOCAL BL increment to q field
   dtl(row_length,rows,bl_levels),                                      &
                                        ! LOCAL BL increment to T field
! DU_STAR, DV_STAR: 1st stage Temporary BL incr to u, v wind
! components from new (stable) BL solver.
   du_star(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end, &
           bl_levels),                                                  &
   dv_star(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end, &
           bl_levels),                                                  &
   dqw_nt(row_length,rows,bl_levels),                                   &
                                        ! OUT NT incr to qw
   dtl_nt(row_length,rows,bl_levels),                                   &
                                        ! OUT NT incr to TL
   ct_ctq(row_length,rows,bl_levels),                                   &
                                        ! LOCAL Coefficient in T and q
                                        !       tri-diagonal implicit
                                        !       matrix
   cq_cm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                        ! LOCAL Coefficient in U
                                        !       tri-diagonal implicit
                                        !       matrix
   cq_cm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                        ! LOCAL Coefficient in V
                                        !       tri-diagonal implicit
                                        !       matrix

   e_land(row_length,rows),                                             &
                                        ! LOCAL FQW over mean land
   ei_land(row_length,rows),                                            &
                                        ! LOCAL EI over mean land
   ftl_land(row_length,rows),                                           &
                                        ! LOCAL FTL over mean land
   fric_heating_blyr(row_length,rows),                                  &
                                        ! LOCAL Frictional heating rate in
                                        ! surface layer
   weight1, weight2, weight3,                                           &
                                        ! LOCAL Weights for f_buoy interp
   ftl_m, fqw_m, f_buoy_m,                                              &
                                        ! LOCAL Fluxes interpd to theta-levs
   dissip_mol,                                                          &
                                        ! LOCAL Molecular dissipation rate
   dissip_0_int,                                                        &
                                        ! LOCAL Integral to rho-level 1
   z_tq,                                                                &
                                        ! LOCAL Height of theta-levels
   z_blyr,                                                              &
                                        ! LOCAL Height of surface layer
   fric_heating_inc,                                                    &
                                        ! LOCAL heating rate
   dz_k, dz_kp1

! Arrays below are needed for frictional dissipation heating source

  REAL, DIMENSION(:,:,:), ALLOCATABLE, TARGET ::                        &
    dissip_u,dissip_v,dissip_u_p,dissip_v_p

! The three set of arrays below are needed by the uncond stable
! BL numerical solver

  REAL, DIMENSION(:,:,:), ALLOCATABLE ::                                &
    dqw1, dtl1, ctctq1

! Following arrays: stable & non-oscillatory solver fluxes for 1st
!                   stage (predictor)

  REAL, DIMENSION(:,:,:), ALLOCATABLE ::                                &
    fqw_star, ftl_star, taux_star, tauy_star


! Automatic arrays for land and sea surface stress diagnostics (used
! as inputs to the ocean when the coupled model is selected)

  REAL, DIMENSION(udims%i_start:udims%i_end,udims%j_start:udims%j_end)::&
    taux_land_star, taux_ssi_star

  REAL, DIMENSION(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)::&
    tauy_land_star, tauy_ssi_star

  LOGICAL :: l_correct

  INTEGER :: i,j,k,l,n
  INTEGER :: i_field ! for swap_bounds_mv

  TYPE(swapable_field_pointer_type) :: fields_to_swap(2)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('IMP_SOLVER',zhook_in,zhook_handle)

! Use standard BL solver
  IF ( .NOT. l_us_blsol ) THEN

! DEPENDS ON: bdy_impl1
    CALL bdy_impl1 (                                                    &
     halo_i,halo_j,offx ,offy ,row_length,rows,n_rows,bl_levels,        &
     global_row_length,n_proc, n_procy, proc_row_group,at_extremity,    &
     q,qcl,qcf,q_latest,qcl_latest,qcf_latest,                          &
     t,t_latest,                                                        &
     dtrdz_charney_grid,dtrdz_u,dtrdz_v,                                &
     rhokh,rhokm_u,rhokm_v,                                             &
     rdz_charney_grid,rdz_u,rdz_v,GAMMA,                                &
     du_nt,dv_nt,                                                       &
     fqw,ftl,taux,tauy,                                                 &
     qw,tl,                                                             &
     ct_ctq,dqw,dtl,cq_cm_u,cq_cm_v,du,dv,                              &
     ltimer                                                             &
      )

! NOTE ON SEA ICE
! The new schemes added to JULES at UM8.2 have not been copied into
! sf_impl.  So if using sf_impl, it is assumed that nice=nice_use=1

! DEPENDS ON: sf_impl
    CALL sf_impl (                                                      &

! IN values defining field dimensions and subset to be processed :
     offx ,offy ,row_length,rows,n_rows,land_pts,                       &

! IN soil/vegetation/land surface data :
     land_index,land_mask,nice,                                         &
     ntiles,tile_index,tile_pts,sm_levels,                              &
     canhc_tile,canopy,flake,smc,                                       &
     tile_frac,wt_ext_tile,                                             &
     fland,flandg,                                                      &

! IN sea/sea-ice data :
     di,ice_fract,di_ncat,ice_fract_ncat,u_0,v_0,                       &

! IN everything not covered so far :
     pstar,lw_down,sw_tile,timestep,                                    &
     t_soil,qw,tl,u,v,                                                  &
     rhokm_u,rhokm_v,GAMMA(1),                                          &
     alpha1,alpha1_sice,ashtf,ashtf_tile,                               &
     dtrdz_charney_grid,du,dv,                                          &
     fqw_tile,epot_tile,fqw_ice,ftl_ice,                                &
     fraca,resfs,resft,rhokh,rhokh_tile,rhokh_sice,rhokpm,rhokpm_pot,   &
     rhokpm_sice,                                                       &
     dtstar_tile,dtstar,                                                &
     z1_tq,z0hssi,z0mssi,z0h_tile,z0m_tile,                             &
     cdr10m_u,cdr10m_v,chr1p5m,chr1p5m_sice,ct_ctq,dqw,dtl,             &
     cq_cm_u,cq_cm_v,                                                   &
     flandg_u,flandg_v,l_sice_heatflux,                                 &
     emis_tile,emis_soil,                                               &

! IN additional variables used to calculate atmospheric cooling
!    at the screen level
     l_co2_interactive, co2_mmr, co2_3d,                                &
     rho1, f3_at_p, ustargbm,                                           &

! IN STASH flags :-
     simlt,smlt,slh,sq1p5,st1p5,su10,sv10,                              &

! INOUT data :
     ti,ti_gb,tstar,                                                    &
     tstar_land,tstar_sea,tstar_sice_cat,tstar_ssi,                     &
     tstar_tile,snow_tile,                                              &
     le_tile,radnet_sice,radnet_tile,                                   &
     e_sea,fqw,ftl,ftl_tile,h_sea,olr,                                  &
     taux,tauy,                                                         &
     taux_land,taux_ssi,tauy_land,tauy_ssi,                             &
     TScrnDcl_SSI,TScrnDcl_TILE,tStbTrans,                              &

! OUT Diagnostic not requiring STASH flags :
     ecan,ei_tile,esoil_tile,                                           &
     sea_ice_htf,surf_ht_flux,surf_ht_flux_land,surf_ht_flux_sice,      &
     surf_htf_tile,                                                     &

! OUT diagnostic requiring STASH flags :
     sice_mlt_htf,snomlt_surf_htf,latent_heat,                          &
     q1p5m,q1p5m_tile,t1p5m,t1p5m_tile,u10m,v10m,                       &

! OUT data required elsewhere in UM system :
     ecan_tile,ei,es,ext,snowmelt,melt_tile,rhokh_mix,                  &
     error,                                                             &

! LOGICAL LTIMER
     lq_mix_bl,                                                         &
     l_flux_bc,                                                         &
     ltimer                                                             &
     )

    tstar_sice(:,:) = tstar_sice_cat(:,:,1)   ! nice_use=1 in this case

! DEPENDS ON: bdy_impl2
    CALL bdy_impl2 (                                                    &

! IN values defining field dimensions and subset to be processed :
     offx ,offy ,row_length,rows,n_rows,                                &

! IN values defining vertical grid of model atmosphere :
     bl_levels,                                                         &

! IN Diagnostic switches
       l_ftl, l_fqw, l_taux, l_tauy,                                    &
! IN data :
     GAMMA,                                                             &
     rhokh,rhokm_u,rhokm_v,rdz_charney_grid,rdz_u,rdz_v,                &

! INOUT data :
     qw,tl,fqw,ftl,taux,tauy,                                           &
     du,dv,ct_ctq,dqw,dtl,cq_cm_u,cq_cm_v,                              &

! OUT data :
     t_latest,q_latest,rhokh_mix,                                       &

! LOGICAL LTIMER
     ltimer                                                             &
      )

  ELSE

! Use new unconditionally stable and non-oscillatory BL solver

    ALLOCATE (dqw1(row_length,rows,bl_levels))
    ALLOCATE (dtl1(row_length,rows,bl_levels))
    ALLOCATE (ctctq1(row_length,rows,bl_levels))

! Following arrays: stable & non-oscillatory solver fluxes for 1st
!                   stage (predictor)

    ALLOCATE (fqw_star(row_length,rows,bl_levels))
    ALLOCATE (ftl_star(row_length,rows,bl_levels))
    ALLOCATE (taux_star(udims%i_start:udims%i_end,                      &
                         udims%j_start:udims%j_end,bl_levels))
    ALLOCATE (tauy_star(vdims%i_start:vdims%i_end,                      &
                         vdims%j_start:vdims%j_end,bl_levels))

!----------------------------------------------------------------------
! Compute 1st stage solution (predictor).
!----------------------------------------------------------------------

! First compute the scheme coefficients for the 1st stage. Make
! coefficients dependent on the BL type for achieving better balance
! between stability-accuracy: stable BL can be strongly nonlinear and
! stiff and thus numerically unstable, so choose a large P value.
! Unstable BL are weakly nonlinear so the solver should be able to cope
! with small P.

!----------------------------------------------------------------------
    sqrt2 = SQRT(2.)

    DO j = 1, rows
      DO i = 1, row_length
        p1=bl_type_1(i,j)*pstb+(1.-bl_type_1(i,j))*puns
        p2=bl_type_2(i,j)*pstb+(1.-bl_type_2(i,j))*puns
        pnonl=MAX(p1,p2)
        i1 = (1.+1./sqrt2)*(1.+pnonl)
        e1 = (1.+1./sqrt2)*( pnonl + (1./sqrt2) +                       &
                            SQRT(pnonl*(sqrt2-1.)+0.5) )
        gamma1(i,j) = i1
        gamma2(i,j) = i1 - e1
      END DO
    END DO

    l_correct = .FALSE.

! DEPENDS ON: BDY_IMPL3
    CALL bdy_impl3 (                                                    &
! IN levels/switches
     bl_levels, l_correct,                                              &
! IN fields
     q,qcl,qcf,q_latest,qcl_latest,qcf_latest,t,t_latest,               &
     dtrdz_charney_grid,dtrdz_u,dtrdz_v,                                &
     rhokh(1,1,2),rhokm_u(udims%i_start,udims%j_start,2),               &
     rhokm_v(vdims%i_start,vdims%j_start,2),                            &
     rdz_charney_grid,rdz_u,rdz_v,gamma1,gamma2,GAMMA,                  &
     du_nt,dv_nt,                                                       &
! INOUT fields
     fqw,ftl,taux,tauy,du,dv,                                           &
! OUT fields
     dqw_nt,dtl_nt,qw,tl,dqw1,dtl1,ct_ctq,ctctq1,dqw,dtl,cq_cm_u,cq_cm_v&
      )

! DEPENDS ON: SF_IMPL2
    CALL sf_impl2 (                                               &
! IN values defining field dimensions and subset to be processed :
     land_pts,land_index,nice,nice_use,ntiles,tile_index,tile_pts,      &
     sm_levels,canhc_tile,canopy,flake,smc,tile_frac,wt_ext_tile,       &
     fland,flandg,lq_mix_bl, l_flux_bc,                                 &
! IN sea/sea-ice data :
     ice_fract,ice_fract_ncat,k_sice,u_0,v_0,                           &
! IN everything not covered so far :
     pstar,lw_down,sw_tile,                                             &
     t_soil,qw,tl,u,v,rhokm_u,rhokm_v,GAMMA(1),                         &
     gamma1,gamma2,alpha1,alpha1_sice,ashtf,ashtf_tile,                 &
     dtrdz_charney_grid,du,dv,                                          &
     fraca,resfs,resft,rhokh,rhokh_tile,rhokh_sice,z1_tq,               &
     z0hssi,z0mssi,z0h_tile,z0m_tile,cdr10m_u,cdr10m_v,                 &
     chr1p5m,chr1p5m_sice,ctctq1,                                       &
     dqw1,dtl1,du_star,dv_star,cq_cm_u,cq_cm_v,                         &
     l_correct,flandg_u,flandg_v,                           &
     emis_tile,ti,tstar_sea,snow_tile,                                  &
! IN variables used to calculate cooling at the screen level 
     l_co2_interactive, co2_mmr, co2_3d, rho1, f3_at_p, ustargbm,       & 
! IN STASH flags :-
     simlt,smlt,slh,sq1p5,st1p5,su10,sv10,                              &
! INOUT data :
     epot_tile,fqw_ice,ftl_ice,dtstar_tile,dtstar,                      &
     tstar_sice_cat,tstar_ssi,tstar_tile,radnet_sice,fqw_tile,          &
     fqw,ftl,ftl_tile,olr,taux_land,taux_ssi,tauy_land,tauy_ssi,        &
     TScrnDcl_SSI,TScrnDcl_TILE,tStbTrans,                              &
! OUT Diagnostic not requiring STASH flags :
     ecan,ei_tile,esoil_tile,sea_ice_htf,surf_ht_flux,                  &
     surf_ht_flux_land,surf_ht_flux_sice,surf_htf_tile,                 &
! OUT diagnostic requiring STASH flags :
     sice_mlt_htf,snomlt_surf_htf,latent_heat,                          &
     q1p5m,q1p5m_tile,t1p5m,t1p5m_tile,u10m,v10m,                       &
! OUT data required elsewhere in UM system :
     tstar,tstar_land,tstar_sice,le_tile,radnet_tile,e_sea,h_sea,       &
     taux,tauy,taux_land_star,tauy_land_star,taux_ssi_star,             &
     tauy_ssi_star,ecan_tile,ei,ei_sice,es,ext,snowmelt,melt_tile,      &
     rhokh_mix,error                                                    &
     )

! DEPENDS ON: BDY_IMPL4
    CALL bdy_impl4 (                                                    &
! IN levels, switches
     bl_levels,  l_correct, l_ftl, l_fqw, l_taux, l_tauy,               &
! IN data :
     gamma1,gamma2,                                                     &
     rhokm_u(udims%i_start,udims%j_start,2),                            &
     rhokm_v(vdims%i_start,vdims%j_start,2),                            &
     rdz_charney_grid,dtrdz_charney_grid,rdz_u,rdz_v,                   &
     ct_ctq,cq_cm_u,cq_cm_v,dqw_nt,dtl_nt,                              &
! INOUT data :
     qw,tl,fqw,ftl,taux,tauy,fqw_star,ftl_star,taux_star,tauy_star,     &
     du,dv,du_star,dv_star,dqw,dtl,rhokh(1,1,2),                        &
! OUT data
     t_latest,q_latest,rhokh_mix                                        &
      )

!----------------------------------------------------------------------
! Compute 2nd stage (final) solution (corrector).
!----------------------------------------------------------------------

! First compute the scheme coefficients for the 2nd stage. Make
! coefficients dependent on the BL type for achieving better balance
! between stability-accuracy: stable BL can be strongly nonlinear and
! stiff and thus numerically unstable, so choose a large P value.
! Unstable BL are weakly nonlinear so the solver should be able to cope
! with small P.

!----------------------------------------------------------------------
    DO j = 1, rows
      DO i = 1, row_length
        p1=bl_type_1(i,j)*pstb+(1.-bl_type_1(i,j))*puns
        p2=bl_type_2(i,j)*pstb+(1.-bl_type_2(i,j))*puns
        pnonl=MAX(p1,p2)
        i1=(1.+1./sqrt2)*(1.+pnonl)
        e2=(1.+1./sqrt2)*( pnonl+(1./sqrt2) -                           &
                          SQRT(pnonl*(sqrt2-1.)+0.5))
        gamma1(i,j) = i1
        gamma2(i,j) = i1 - e2
      END DO
    END DO

    l_correct = .TRUE.
! DEPENDS ON: BDY_IMPL3
    CALL bdy_impl3 (                                                    &
! IN levels/switches
     bl_levels,l_correct,                                               &
! IN fields
     q,qcl,qcf,q_latest,qcl_latest,qcf_latest,t,t_latest,               &
     dtrdz_charney_grid,dtrdz_u,dtrdz_v,                                &
     rhokh(1,1,2),rhokm_u(udims%i_start,udims%j_start,2),               &
     rhokm_v(vdims%i_start,vdims%j_start,2),                            &
     rdz_charney_grid,rdz_u,rdz_v,gamma1,gamma2,GAMMA,                  &
     du_nt,dv_nt,                                                       &
! INOUT fields
     fqw,ftl,taux,tauy,du,dv,                                           &
     dqw_nt,dtl_nt,qw,tl,dqw1,dtl1,ct_ctq,ctctq1,dqw,dtl,cq_cm_u,cq_cm_v&
      )

! DEPENDS ON: SF_IMPL2
    CALL sf_impl2 (                                                     &
! IN values defining field dimensions and subset to be processed :
     land_pts,land_index,nice,nice_use,ntiles,tile_index,tile_pts,      &
     sm_levels,canhc_tile,canopy,flake,smc,tile_frac,wt_ext_tile,       &
     fland,flandg,lq_mix_bl,l_flux_bc,                                  &
! IN sea/sea-ice data :
     ice_fract,ice_fract_ncat,k_sice,u_0,v_0,                           &
! IN everything not covered so far :
     pstar,lw_down,sw_tile,                                             &
     t_soil,qw,tl,u,v,rhokm_u,rhokm_v,GAMMA(1),                         &
     gamma1,gamma2,alpha1,alpha1_sice,ashtf,ashtf_tile,                 &
     dtrdz_charney_grid,du,dv,                                          &
     fraca,resfs,resft,rhokh,rhokh_tile,rhokh_sice,z1_tq,               &
     z0hssi,z0mssi,z0h_tile,z0m_tile,cdr10m_u,cdr10m_v,                 &
     chr1p5m,chr1p5m_sice,ctctq1,                                       &
     dqw1,dtl1,du_star,dv_star,cq_cm_u,cq_cm_v,                         &
     l_correct,flandg_u,flandg_v,                           &
     emis_tile,ti,tstar_sea,snow_tile,                                  &
! IN variables used to calculate cooling at the screen level
     l_co2_interactive, co2_mmr, co2_3d, rho1, f3_at_p, ustargbm,       &
! IN STASH flags :-
     simlt,smlt,slh,sq1p5,st1p5,su10,sv10,                              &
! INOUT data :
     epot_tile,fqw_ice,ftl_ice,dtstar_tile,dtstar,                      &
     tstar_sice_cat,tstar_ssi,tstar_tile,radnet_sice,fqw_tile,          &
     fqw,ftl,ftl_tile,olr,taux_land,taux_ssi,tauy_land,tauy_ssi,        &
     TScrnDcl_SSI,TScrnDcl_TILE,tStbTrans,                              &
! OUT Diagnostic not requiring STASH flags :
     ecan,ei_tile,esoil_tile,sea_ice_htf,surf_ht_flux,                  &
     surf_ht_flux_land,surf_ht_flux_sice,surf_htf_tile,                 &
! OUT diagnostic requiring STASH flags :
     sice_mlt_htf,snomlt_surf_htf,latent_heat,                          &
     q1p5m,q1p5m_tile,t1p5m,t1p5m_tile,u10m,v10m,                       &
! OUT data required elsewhere in UM system :
     tstar,tstar_land,tstar_sice,le_tile,radnet_tile,e_sea,h_sea,       &
     taux,tauy,taux_land_star,tauy_land_star,taux_ssi_star,             &
     tauy_ssi_star,ecan_tile,ei,ei_sice,es,ext,snowmelt,melt_tile,      &
     rhokh_mix,error                                                    &
     )

! DEPENDS ON: BDY_IMPL4
    CALL bdy_impl4 (                                                    &
! IN levels, switches
     bl_levels, l_correct, l_ftl, l_fqw, l_taux, l_tauy,                &
! IN data :
     gamma1,gamma2,                                                     &
     rhokm_u(udims%i_start,udims%j_start,2),                            &
     rhokm_v(vdims%i_start,vdims%j_start,2),                            &
     rdz_charney_grid,dtrdz_charney_grid,rdz_u,rdz_v,                   &
     ct_ctq,cq_cm_u,cq_cm_v,dqw_nt,dtl_nt,                              &
! INOUT data :
     qw,tl,fqw,ftl,taux,tauy,fqw_star,ftl_star,taux_star,tauy_star,     &
     du,dv,du_star,dv_star,dqw,dtl,rhokh(1,1,2),                        &
! OUT data, NB these are really tl and qt on exit!
     t_latest,q_latest,rhokh_mix                                        &
      )

    DEALLOCATE (dqw1)
    DEALLOCATE (dtl1)
    DEALLOCATE (ctctq1)

! Following arrays: stable & non-oscillatory solver fluxes for 1st
!                   stage (predictor)

    DEALLOCATE (fqw_star)
    DEALLOCATE (ftl_star)
    DEALLOCATE (taux_star)
    DEALLOCATE (tauy_star)

  END IF ! IF L_us_blsol

    DO j = 1, rows
      DO i = 1, row_length
        ftl_land(i,j)=0.0
        ftl_ssi(i,j)=0.0
        e_land(i,j)=0.0
        e_ssi(i,j)=0.0
        ei_land(i,j)=0.0
      END DO
    END DO

    DO n = 1, ntiles
      DO k = 1, tile_pts(n)
        l = tile_index(k,n)
        j=(land_index(l)-1)/row_length + 1
        i = land_index(l) - (j-1)*row_length
        ftl_land(i,j)=ftl_land(i,j) +                                   &
          ftl_tile(l,n)*tile_frac(l,n)
        e_land(i,j)=e_land(i,j) +                                       &
          fqw_tile(l,n)*tile_frac(l,n)
        ei_land(i,j)=ei_land(i,j) +                                     &
          ei_tile(l,n)*tile_frac(l,n)
      END DO
    END DO

    DO j = 1, rows
      DO i = 1, row_length
        IF(flandg(i,j) <  1.0)THEN
          ftl_ssi(i,j)=(ftl(i,j,1)-ftl_land(i,j)*flandg(i,j))           &
            /(1.0-flandg(i,j))
          e_ssi(i,j)=(fqw(i,j,1)-e_land(i,j)*flandg(i,j))               &
            /(1.0-flandg(i,j))
        END IF
      END DO
    END DO

  IF (sq1p5 .AND. lq_mix_bl) THEN
!----------------------------------------------------------------------
! Convert 1.5m mixing ratio to specific humidity
!-----------------------------------------------
! Not immediately clear what to do about QCL and QCF at 1.5m so assume
! the same as level 1.  An alternative would be to assume zero.
!----------------------------------------------------------------------
    DO j = 1, rows
      DO i = 1, row_length
        q1p5m(i,j)=q1p5m(i,j)/(1.0+q1p5m(i,j)+qcl(i,j,1)+qcf(i,j,1))
      END DO
    END DO

    DO n = 1, ntiles
      DO k = 1, tile_pts(n)
        l = tile_index(k,n)
        j=(land_index(l)-1)/row_length + 1
        i = land_index(l) - (j-1)*row_length

        q1p5m_tile(l,n) = q1p5m_tile(l,n)/ (1.0 + q1p5m_tile(l,n)       &
                                             +qcl(i,j,1)+qcf(i,j,1) )
      END DO
    END DO
  END IF

  IF (.NOT. l_sice_multilayers) THEN

!-----------------------------------------------------------------------
! Update sea-ice surface layer temperature, if not coupled to multilayer
! sea ice model.
! This routine should appear further up the code tree, but it is unclear
! at present how this would impact on the diagnostics,
! so put it here for now
!-----------------------------------------------------------------------

! DEPENDS ON: sice_htf
  CALL sice_htf(                                                        &
   row_length,rows,flandg,simlt,nice,nice_use,                          &
   di_ncat,ice_fract,ice_fract_ncat,surf_ht_flux_sice,                  &
   tstar_sice_cat,timestep,                                             &
   ti,ti_gb,sice_mlt_htf,sea_ice_htf,l_sice_heatflux,                   &
   ltimer)
  
  ELSE
    sea_ice_htf(:,:,:) =0.0   ! for safety
  ENDIF

! Convert sea and sea-ice fluxes to be fraction of grid-box
! (as required by sea and sea-ice modellers)

  DO n = 1, nice
    DO j = 1, rows
      DO i = 1, row_length
        surf_ht_flux_sice(i,j,n)=ice_fract_ncat(i,j,n)*                 &
                            surf_ht_flux_sice(i,j,n)
        sice_mlt_htf(i,j,n)=ice_fract_ncat(i,j,n)*sice_mlt_htf(i,j,n)
        sea_ice_htf(i,j,n)=ice_fract_ncat(i,j,n)*sea_ice_htf(i,j,n)
      END DO
    END DO
  END DO

  IF ( fric_heating == on ) THEN

    ALLOCATE (dissip_u(udims_s%i_start:udims_s%i_end,                   &
                       udims_s%j_start:udims_s%j_end,bl_levels))
    ALLOCATE (dissip_v(vdims_s%i_start:vdims_s%i_end,                   &
                       vdims_s%j_start:vdims_s%j_end,bl_levels))
    ALLOCATE (dissip_u_p(row_length,rows,bl_levels))
    ALLOCATE (dissip_v_p(row_length,rows,bl_levels))
!----------------------------------------------------------------------
! Add heating increment from turbulence dissipation
!--------------------------------------------------
! Calculate dissipation rate on theta-levels,
! recalling that TAU(K) is defined on theta_level(K-1)
!----------------------------------------------------------------------

!  Calculate as TAU * DU/DZ

    DO k = 1, bl_levels-1
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          dissip_u(i,j,k) = taux(i,j,k+1) * rdz_u(i,j,k+1) *            &
                      (u(i,j,k+1)+du(i,j,k+1)-u(i,j,k)-du(i,j,k))
        END DO
      END DO
    END DO

    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
          ! Also need to account for dissipation between the surface and
          ! uv-level 1 in the theta-level(K=1) heating term.
          ! Total theta-level 1 dissipation is the cell-average between
          ! the surface and rho-level 2.
          ! First calculate integral to rho-level 1:
        dissip_0_int = taux(i,j,1) * ( u(i,j,1)+du(i,j,1)-u_0(i,j) )
!!           Note integration implies multiplication by zrho1-zth0
!!           hence no division by zrho1-zth0
          ! Then take cell average up to rho-level 2:
        dissip_u(i,j,1) = ( dissip_0_int + dissip_u(i,j,1)*             &

                        (r_at_u(i,j,2)-r_at_u(i,j,1)) )                 &
                   / ( r_at_u(i,j,2)-                                   &
          0.5*(r_theta_levels(i,j,0)+r_theta_levels(i+1,j,0)) )




          ! stress equals zero on top theta-level
        dissip_u(i,j,bl_levels) = 0.0

      END DO
    END DO

    DO k = 1, bl_levels-1
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          dissip_v(i,j,k) = tauy(i,j,k+1) * rdz_v(i,j,k+1) *            &
                      (v(i,j,k+1)+dv(i,j,k+1)-v(i,j,k)-dv(i,j,k))
        END DO
      END DO
    END DO

    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
          ! Also need to account for dissipation between the surface and
          ! uv-level 1 in the theta-level(K=1) heating term
          ! Total theta-level 1 dissipation is the cell-average between
          ! the surface and rho-level 2.
          ! First calculate integral to rho-level 1:
        dissip_0_int = tauy(i,j,1) * ( v(i,j,1)+dv(i,j,1)-v_0(i,j) )
        dissip_v(i,j,1) = ( dissip_0_int + dissip_v(i,j,1)*             &

                        (r_at_v(i,j,2)-r_at_v(i,j,1)) )                 &
                   / ( r_at_v(i,j,2)-                                   &
          0.5*(r_theta_levels(i,j,0)+r_theta_levels(i,j+1,0)) )




          ! stress equals zero on top theta-level
        dissip_v(i,j,bl_levels) = 0.0

      END DO
    END DO



    i_field = 1
    fields_to_swap(i_field) % field       => dissip_u(:,:,:)
    fields_to_swap(i_field) % field_type  =  fld_type_u
    fields_to_swap(i_field) % levels      =  bl_levels
    fields_to_swap(i_field) % rows        =  rows
    fields_to_swap(i_field) % vector      =  .TRUE.

    i_field = i_field + 1
    fields_to_swap(i_field) % field       => dissip_v(:,:,:)
    fields_to_swap(i_field) % field_type  =  fld_type_v
    fields_to_swap(i_field) % levels      =  bl_levels
    fields_to_swap(i_field) % rows        =  n_rows
    fields_to_swap(i_field) % vector      =  .TRUE.

! DEPENDS ON: swap_bounds_mv
    CALL swap_bounds_mv(fields_to_swap, i_field, row_length,            &
                      offx ,          offy )

        ! Interpolate u and v components to p grid.
      CALL u_to_p(dissip_u,                                             &
                        udims_s%i_start,udims_s%i_end,                  & 
                        udims_s%j_start,udims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        bl_levels,                                      &
                        model_domain,at_extremity,dissip_u_p)

      CALL v_to_p(dissip_v,                                             &
                        vdims_s%i_start,vdims_s%i_end,                  & 
                        vdims_s%j_start,vdims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        bl_levels,                                      &
                        model_domain,at_extremity,dissip_v_p)        


    IF(model_domain  ==  mt_global) THEN

IF (.NOT. l_vatpoles) THEN

! Set values at the poles to zero
        IF (at_extremity(PSouth) ) THEN
          DO k = 1, bl_levels
            DO i = 1, row_length
              dissip_v_p(i,1,k) = 0.0
              dissip_u_p(i,1,k) = 0.0
            END DO
          END DO
        END IF
        IF (at_extremity(PNorth) ) THEN
          DO k = 1, bl_levels
            DO i = 1, row_length
              dissip_v_p(i,rows,k) = 0.0
              dissip_u_p(i,rows,k) = 0.0
            END DO
          END DO
        END IF

END IF ! vatpoles

    END IF



!-----------------------------------------------------------------------
! First, estimate molecular dissipation rate by assuming steady state
! subgrid KE, so that     dissip_mol = dissip_rke + f_buoy
! where f_buoy is the buoyancy flux interpolated to theta levels
!-----------------------------------------------------------------------
! Then convert dissipation rate to heating rate,
! noting that dissipation rates are zero on BL_LEVELS,
! and redistribute within the boundary layer to account
! for lack of BL mixing of these increments
!-----------------------------------------------------------------------
    DO j = 1, rows
      DO i = 1, row_length
        nblyr(i,j) = 1
        k = 1

        weight1 = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,0)
        weight2 = r_theta_levels(i,j,k) - r_theta_levels(i,j,0)
        weight3 = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k)
        ftl_m = weight2 * ftl(i,j,k+1) + weight3 * ftl(i,j,k)
        fqw_m = weight2 * fqw(i,j,k+1) + weight3 * fqw(i,j,k)
        f_buoy_m = g*( bt_gb(i,j,k)*(ftl_m/cp) +                        &
                       bq_gb(i,j,k)*fqw_m )/weight1

        dissip_mol = dissip_u_p(i,j,k)+dissip_v_p(i,j,k)                &
                              + f_buoy_m
        fric_heating_inc = MAX (0.0, timestep * dissip_mol              &
                                       / ( cp*rho_wet_theta(i,j,k) ) )

          ! Save level 1 heating increment for redistribution over
          ! boundary layer
        fric_heating_blyr(i,j) = fric_heating_inc *                     &
                         (r_rho_levels(i,j,2)-r_theta_levels(i,j,0))

      END DO
    END DO
    DO k = 2, bl_levels-1
      DO j = 1, rows
        DO i = 1, row_length

          weight1 = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
          weight2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
          weight3 = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k)
          ftl_m = weight2 * ftl(i,j,k+1) + weight3 * ftl(i,j,k)
          fqw_m = weight2 * fqw(i,j,k+1) + weight3 * fqw(i,j,k)

          f_buoy_m = g*( bt_gb(i,j,k)*(ftl_m/cp) +                      &
                         bq_gb(i,j,k)*fqw_m )/weight1

          dissip_mol = dissip_u_p(i,j,k)+dissip_v_p(i,j,k)              &
                              + f_buoy_m
          fric_heating_inc = MAX (0.0, timestep * dissip_mol            &
                                       / ( cp*rho_wet_theta(i,j,k) ) )

          z_tq = r_theta_levels(i,j,k)-r_theta_levels(i,j,0)
          IF ( z_tq <= zh(i,j) ) THEN
              !------------------------------------------------------
              ! Sum increments over boundary layer to avoid
              ! adding large increments in level 1
              !------------------------------------------------------
            nblyr(i,j) = k
            fric_heating_blyr(i,j) = fric_heating_blyr(i,j) +           &
                         fric_heating_inc *                             &
                         (r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))
          ELSE
            t_latest(i,j,k) = t_latest(i,j,k) + fric_heating_inc
            IF (BL_diag%l_dtfric) THEN
              BL_diag%dTfric(i,j,k) = fric_heating_inc
            END IF
          END IF

        END DO
      END DO
    END DO
!-----------------------------------------------------------------------
! Redistribute heating within the boundary layer to account
! for lack of BL mixing of these increments
!-----------------------------------------------------------------------
    DO j = 1, rows
      DO i = 1, row_length

        z_blyr = r_rho_levels(i,j,nblyr(i,j)+1)                         &
                                  - r_theta_levels(i,j,0)
        fric_heating_blyr(i,j) = fric_heating_blyr(i,j) / z_blyr

        DO k = 1, nblyr(i,j)

            ! Linearly decrease heating rate across surface layer
          z_tq = r_theta_levels(i,j,k)-r_theta_levels(i,j,0)
          fric_heating_inc = 2.0 * fric_heating_blyr(i,j) *             &
                                (1.0-z_tq/z_blyr)

          t_latest(i,j,k) = t_latest(i,j,k) + fric_heating_inc

          IF (BL_diag%l_dtfric) THEN
            BL_diag%dTfric(i,j,k) = fric_heating_inc
          END IF

        END DO
      END DO
    END DO

    DEALLOCATE (dissip_v_p)
    DEALLOCATE (dissip_u_p)
    DEALLOCATE (dissip_v)
    DEALLOCATE (dissip_u)

  END IF

  IF (lhook) CALL dr_hook('IMP_SOLVER',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE imp_solver
