! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!!! Subroutine IMPS_INTCT -------------------------------------------
!!!
!!! Purpose : Intermediate control level to call requested version of
!!!           IMP_SOLVER with the appropriate arguments.
!!!
!!! Programming standard : unified model documentation paper No 3
!!!
!!! System components covered : P24
!!!
!!! System task : P0
!!!
!!!END -----------------------------------------------------------------

!    Arguments :-
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Boundary Layer
SUBROUTINE imps_intct (                                                 &

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
 r_u,r_v,                                                               &

! IN soil/vegetation/land surface data :
 land_mask,land_index,                                                  &
 st_levels,sm_levels,                                                   &

! IN sea/sea-ice data :
 di,ice_fract,di_ncat,ice_fract_ncat,k_sice,u_0,v_0,                    &

! IN cloud data :
 q,qcf,qcl,q_conv,qcf_conv,qcl_conv,qcf_latest,qcl_latest,              &
 t,t_conv,                                                              &

! IN everything not covered so far :
 rho_wet_theta,pstar,timestep,l_sice_heatflux,                          &
 l_sice_multilayers,                                                    &

! IN variables from BDY_LAYR (that used to be local arrays)
 alpha1_sice,ashtf,bq_gb,bt_gb,dtrdz_charney_grid,                      &
 rdz_charney_grid,dtrdz_u,dtrdz_v,rdz_u,rdz_v,                          &
 cdr10m_u,cdr10m_v,z1_tq,zh,                                            &
!ajm extra variable added
 rhokm_u,rhokm_v,                                                       &

! IN variables for new BL solver
 bl_type_1,bl_type_2,bl_type_3,bl_type_4,                               &
 bl_type_5,bl_type_6,bl_type_7,                                         &

! IN additional variables for MOSES II
 tile_pts,tile_index,tile_frac,canopy,                                  &
 alpha1,fraca,rhokh_tile,smc,chr1p5m,                                   &
 resfs,z0hssi,z0mssi,canhc_tile,flake,                                  &
 wt_ext_tile,lw_down,sw_tile,ashtf_tile,                                &
 fqw_ice,ftl_ice,resft,rhokh_sice,rhokpm,rhokpm_pot,rhokpm_sice,        &
 z0h_tile,z0m_tile,chr1p5m_sice,                                        &
 fland,flandg,flandg_u,flandg_v,tstar_sea,                              &

! IN additional variables used to calculate atmospheric cooling
!    at the screen level
 l_co2_interactive, co2_mmr, co2_3d,                                    &
 rho1, f3_at_p, uStarGBM,                                               &

! IN additional variables for JULES
 rhokh_mix_dummy,dtstar_tile,dtstar,hcons,emis_tile,emis_soil,          &

! INOUT data :
 t_soil, ti, ti_gb, tstar, t_latest,q_latest,                           &

! INOUT Diagnostics started in BDY_LAYR not requiring STASH flags :
 e_sea,fqw,ftl,h_sea,rhokh,taux,tauy,                                   &

! INOUT additional variables for MOSES II
 tstar_tile,fqw_tile,epot_tile,ftl_tile,                                &
 snow_tile,le_tile,radnet_sice,radnet_tile,olr,                         &
 tstar_sice,tstar_sice_cat,tstar_ssi,                                   &
 taux_land,taux_ssi,tauy_land,tauy_ssi,                                 &
 TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans,                                &

! OUT Increments to U and V momentum fields and Tl qw
!  (New dynamics only).
 du,dv,                                                                 &

! OUT Diagnostic not requiring STASH flags :
 rhokh_mix,sea_ice_htf,surf_ht_flux,                                    &
 surf_ht_flux_land,surf_ht_flux_sice,                                   &
 surf_htf_tile,                                                         &

! OUT diagnostics requiring STASH flags :
 sice_mlt_htf,snomlt_surf_htf,latent_heat,                              &
 q1p5m,t1p5m,u10m,v10m,                                                 &

! (IN) STASH flags :-
 simlt,smlt,slh,sq1p5,st1p5,su10,sv10,l_ftl,l_fqw,l_taux,l_tauy,        &

! OUT additional variables for MOSES II
 esoil_tile,es,ei_tile,                                                 &
 q1p5m_tile,t1p5m_tile,ecan_tile,melt_tile,                             &
 tstar_land,e_ssi,ei_sice,ftl_ssi,                                      &
 error,                                                                 &

! OUT data required elsewhere in UM system :
 ecan,ei,ext,snowmelt,                                                  &
 BL_diag,lq_mix_bl,                                                     &

! SCM namelist logical
 l_flux_bc,                                                             &
 ltimer                                                                 &
 )

  USE nstypes
  USE atm_fields_bounds_mod, ONLY:                                      &
   udims, vdims, udims_s, vdims_s, udims_l, vdims_l
  USE atmos_constants_mod, ONLY: cp


  USE water_constants_mod, ONLY: lc, lf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE bl_diags_mod, ONLY :                                              &
      strnewbldiag

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
                                 ! IN No. of sea ice catagories
   nice_use                      ! IN No. of sea ice categories used
                                 !    fully in surface exchange

! (b) Defining vertical grid of model atmosphere.

  INTEGER ::                                                            &
   bl_levels,                                                           &
                                   ! IN Max. no. of "boundary" levels
!                                     allowed.Assumed <= 30 for dim-
!                                     sioning of GAMMA in common deck
!                                     C_GAMMA used in SF_EXCH and KMKH
   model_domain

! SCM logical for prescribing surface flux forcing
  LOGICAL ::                                                            &
    l_flux_bc

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
    r_u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,    &
          bl_levels),                                                   &
                            ! non-turbulent increment to u wind field
    r_v(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,    &
          bl_levels),                                                   &
                            ! non-turbulen increment to v wind field
    u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,      &
       bl_levels),                                                      &
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


! (d) Sea/sea-ice data.

  REAL ::                                                               &
   di(row_length,rows),                                                 &
                                   ! IN "Equivalent thickness" of
                                   !   sea-ice(m).
   ice_fract(row_length,rows),                                          &
                                   ! IN Fraction of gridbox covered by
                                   !   sea-ice (decimal fraction).
   di_ncat(row_length,rows,nice),                                       &
                                   ! IN "Equivalent thickness"
                                   !   of ice category in grid box (m).
   ice_fract_ncat(row_length,rows,nice),                                &
                                   ! IN Fraction of ice category
                                   !   in gridbox (decimal fraction).
   k_sice(row_length, rows, nice),                                      &
                                   ! IN sea ice effective conductivity in
                                   !  sfc layer on categories (W/m2/K)
   u_0(udims%i_start:udims%i_end,udims%j_start:udims%j_end),            &
                                   ! IN W'ly component of surface
                                   !   current (m/s).
   v_0(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                   ! IN S'ly component of surface
                                   !   current (m/s).

! (e) Cloud data.

  REAL ::                                                               &
   qcf(1-halo_i:row_length+halo_i,                                      &
       1-halo_j:rows+halo_j,bl_levels),                                 &
                                          ! IN Cloud ice (kg per kg air)
   qcl(1-halo_i:row_length+halo_i,                                      &
       1-halo_j:rows+halo_j,bl_levels),                                 &
                                          ! IN Cloud liquid water
   q(1-halo_i:row_length+halo_i,                                        &
     1-halo_j:rows+halo_j,bl_levels),                                   &
                                          ! IN specific humidity
   t(row_length,rows,bl_levels),                                        &
                                          ! IN temperature
                          ! Latest estimates to time level n+1 values
   qcf_latest(row_length,rows,bl_levels),                               &
                                ! IN Cloud ice (kg per kg air)
   qcl_latest(row_length,rows,bl_levels),                               &
                                             ! IN Cloud liquid water
   t_conv(row_length, rows, bl_levels),                                 &
   q_conv(row_length, rows, bl_levels),                                 &
   qcl_conv(row_length, rows, bl_levels),                               &
   qcf_conv(row_length, rows, bl_levels)

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
   alpha1_sice(row_length,rows,nice_use),                               &
                                  ! IN ALPHA1 for sea-ice.
                                  !    gridbox ALPHA1 for MOSES I
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
                                          !  IN RDZ (K > 1) on UV-grid.
   rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
          2:bl_levels),                                                 &
                                !  IN RDZ (K > 1) on UV-grid.
   cdr10m_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
                                ! IN Ratio of CD's reqd for calculation
   cdr10m_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),       &
                                ! IN Ratio of CD's reqd for calculation
   z1_tq(row_length,rows),                                              &
                                ! IN Height of lowest theta level.
   zh(row_length,rows)      ! IN BL depth (m)

  LOGICAL :: ltimer               ! Logical switch for TIMER diags

! IN additional variables for MOSES II

  INTEGER ::                                                            &
   tile_pts(ntype),                                                     &
                                 ! IN Number of tile points.
   tile_index(land_pts,ntype)
                                 ! IN Index of tile points.

  REAL ::                                                               &
   tile_frac(land_pts,ntiles),                                          &
                                ! IN fractional coverage for each
                                !    surface tile
   canopy(land_pts,ntiles),                                             &
                                ! IN Surface/canopy water (kg/m2)
   alpha1(land_pts,ntiles),                                             &
                               ! IN Mean gradient of saturated
!                                 specific humidity with
!                                 respect to temperature between
!                                 the bottom model layer and the
!                                 tile surfaces.
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
   fland(land_pts),                                                     &
                                ! IN Land fraction on land points.
   flandg(row_length,rows),                                             &
                                ! IN Land fraction on all points.
   flandg_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
                                ! IN Land frac (on U-grid, with 1st
                                !    and last rows undefined or, at
                                !    present, set to "missing data")
   flandg_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),       &
                                ! IN Land frac (on V-grid, with 1st
                                !    and last rows undefined or, at
                                !    present, set to "missing data")
   tstar_sea(row_length,rows)   ! IN Open sea sfc temperature (K).

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

!     Declaration of new BL diagnostics.
  TYPE (strnewbldiag) :: BL_diag

!  In/outs :-

  REAL ::                                                               &
   q_latest(row_length,rows,bl_levels),                                 &
                                            ! IN specific humidity
   t_latest(row_length,rows,bl_levels),                                 &
                                            ! IN temperature
   t_soil(land_pts,sm_levels),                                          &
                                ! INOUT Soil temperatures (K).
   ti(row_length,rows,nice),                                            &
                                ! INOUT Sea-ice surface layer
                                !    temperature on categories (K)
                                ! (IN only if l_sice_multilayers=T) 
   ti_gb(row_length,rows),                                              &
                                ! INOUT Sea-ice surface layer temp (ice mean)
                                ! (IN only if l_sice_multilayers=T) 
   tstar(row_length,rows)       ! INOUT Surface temperature (K).


! INOUT additional variables for MOSES II
  REAL ::                                                               &
   tstar_tile(land_pts,ntiles),                                         &
                                ! INOUT Surface tile temperature
   fqw_tile(land_pts,ntiles),                                           &
                                ! INOUT surface tile moisture flux
   epot_tile(land_pts,ntiles),                                          &
                                ! INOUT surface tile potential
                                !       evaporation
   ftl_tile(land_pts,ntiles),                                           &
                                ! INOUT surface tile heat flux
   snow_tile(land_pts,ntiles),                                          &
                                ! INOUT Snow on tiles (kg/m2).
   le_tile(land_pts,ntiles),                                            &
                                ! INOUT Surface latent heat flux for
                                !       land tiles (W/m2).
   radnet_sice(row_length,rows,nice_use),                               &
                                ! INOUT Sea-ice surface net radiation.
   radnet_tile(land_pts,ntiles),                                        &
                                ! INOUT Tile surface net radiation.
   olr(row_length,rows),                                                &
                                ! IN    TOA - surface upward LW on
                                !       last radiation timestep
                                ! OUT   Corrected TOA outward LW
   tstar_sice(row_length,rows),                                         &
                                ! OUT Sea-ice sfc temperature (K)
                                !    (ice mean over categories)
   tstar_sice_cat(row_length,rows,nice_use),                            &
                                ! INOUT Sea-ice sfc temperature (K).
   tstar_ssi(row_length,rows)
                                ! INOUT Sea mean sfc temperature (K).

  REAL, INTENT(INOUT) :: TScrnDcl_SSI(row_length,rows)
                            !    Decoupled screen-level temperature
                            !    over sea or sea-ice
  REAL, INTENT(INOUT) :: TScrnDcl_TILE(land_pts,ntiles)
                            !    Decoupled screen-level temperature
                            !    over land tiles
  REAL, INTENT(INOUT) :: tStbTrans(row_length,rows)
                            !    Time since the transition


!  Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-

!  (a) Calculated anyway (use STASH space from higher level) :-

  REAL ::                                                               &
   e_sea(row_length,rows),                                              &
                                ! OUT Evaporation from sea times
!                                     leads fraction. Zero over land.
!                                     (kg per square metre per sec).
   fqw(row_length,rows,bl_levels),                                      &
                                ! OUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
   ftl(row_length,rows,bl_levels),                                      &
                                ! OUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
   h_sea(row_length,rows),                                              &
                                ! OUT Surface sensible heat flux over
!                                     sea times leads fraction. (W/m2)
   rhokh(row_length,rows,bl_levels),                                    &
                                ! OUT Exchange coeffs for moisture.
   rhokh_mix(row_length,rows,bl_levels),                                &
                                ! OUT Exchange coeffs for moisture.
! for use in tracer mixing routines
   rhokm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
            bl_levels),                                                 &
                                ! OUT Exchange coefficients for u
   rhokm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
            bl_levels),                                                 &
                                ! OUT Exchange coefficients for v
   sea_ice_htf(row_length,rows,nice),                                   &
                                ! OUT Heat flux through sea-ice
!                                     (W/m2, positive downwards).
   surf_ht_flux(row_length,rows),                                       &
                                ! OUT Net downward heat flux at
!                                     surface over land or sea-ice
!                                     fraction of gridbox (W/m2).
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
   taux(udims%i_start:udims%i_end,udims%j_start:udims%j_end,            &
         bl_levels),                                                    &
                                ! OUT W'ly component of surface wind
!                                     stress (N/sq m).(On UV-grid with
!                                     first and last rows undefined or
!                                     at present, set to missing data
   tauy(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,            &
         bl_levels)
                                ! OUT S'ly component of surface wind
!                                     stress (N/sq m).  On UV-grid;
!                                     comments as per TAUX.

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
   ext(land_pts,sm_levels),                                             &
                                ! OUT Extraction of water from each
                                !     soil layer (kg/m2/s).
   snowmelt(row_length,rows),                                           &
                                ! OUT Snowmelt (kg/m2/s).
   du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,      &
        bl_levels),                                                     &
                                ! OUT BL increment to u wind field
   dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,      &
        bl_levels)          ! OUT BL increment to u wind field


! OUT additional variables for MOSES II
  REAL ::                                                               &
   esoil_tile(land_pts,ntiles),                                         &
                                ! OUT Evaporation from bare soil (kg/m2)
   es(row_length,rows),                                                 &
                                ! OUT Surface evapotranspiration from
                                !     soil moisture store (kg/m2/s).
   ei_tile(land_pts,ntiles),                                            &
                                ! OUT EI for land tiles
   q1p5m_tile(land_pts,ntiles),                                         &
                                ! OUT Q1P5M over land tiles.
   t1p5m_tile(land_pts,ntiles),                                         &
                                ! OUT T1P5M over land tiles.
   ecan_tile(land_pts,ntiles),                                          &
                                ! OUT ECAN for land tiles
   melt_tile(land_pts,ntiles),                                          &
                                ! OUT Snowmelt on tiles (kg/m2/s).
   tstar_land(row_length,rows),                                         &
                                ! OUT Land mean sfc temperature (K)
   e_ssi(row_length,rows),                                              &
                                ! OUT Mean sea moisture heat flux.
   ei_sice(row_length,rows,nice_use),                                   &
                                ! OUT Sea-ice sublimation rate
                                !     (sea mean).
   ftl_ssi(row_length,rows),                                            &
                                ! OUT Mean sea surface heat flux.
   taux_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end),      &
                                ! OUT W'ly component of land sfc wind
                                !     stress (N/sq m). (On U-grid
                                !     with first and last rows
                                !     undefined or, at present,
                                !     set to missing data
   taux_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
                                ! OUT W'ly compt of mean sea sfc wind
                                !     stress (N/sq m). (On U-grid
                                !     with first and last rows
                                !     undefined or, at present,
                                !     set to missing data
   tauy_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),      &
                                ! OUT S'ly componentt of land sfc wind
                                !     stress (N/sq m).  On V-grid;
                                !     comments as per TAUX.
   tauy_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                ! OUT S'ly compt of mean sea sfc wind
                                !     stress (N/sq m).  On V-grid;
                                !     comments as per TAUX.

  INTEGER ::                                                            &
   error            ! OUT 0 - AOK;
                        !     1 to 7  - bad grid definition detected;


!-----------------------------------------------------------------------
! Derived local parameters.

  REAL :: lcrcp,ls,lsrcp

  PARAMETER (                                                           &
   lcrcp=lc/cp,                                                         &
                             ! Evaporation-to-dT conversion factor.
   ls=lf+lc,                                                            &
                             ! Latent heat of sublimation.
   lsrcp=ls/cp                                                          &
                             ! Sublimation-to-dT conversion factor.
    )

  INTEGER :: i, j, k

!  Workspace :-

  REAL ::                                                               &
    du_nt(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,  &
          bl_levels),                                                   &
                            ! non-turbulent increment to u wind field
    dv_nt(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,  &
          bl_levels)    ! non-turbulen increment to v wind field

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------


  IF (lhook) CALL dr_hook('IMPS_INTCT',zhook_in,zhook_handle)

! DEPENDS ON: imp_solver
    CALL imp_solver (                                                   &

! IN MPP variables
     halo_i, halo_j, offx , offy , row_length, rows, n_rows,            &
     global_row_length, proc_row_group, at_extremity,                   &
     n_proc, n_procx, n_procy, neighbour,                               &

! IN values defining field dimensions and subset to be processed :
     ntiles,land_pts,nice,nice_use,                                     &

! IN values defining vertical grid of model atmosphere :
     model_domain,                                                      &
     bl_levels,                                                         &
     GAMMA,                                                             &

! IN U and V momentum fields.
     u, v,                                                              &
! IN Non turbulent increments to momentum
!  (New dynamics only).
     r_u, r_v,                                                          &

! IN soil/vegetation/land surface data :
     land_mask,land_index,                                              &
     st_levels,sm_levels,tile_frac,canopy,                              &
     fland,flandg,                                                      &

! IN sea/sea-ice data :
     di,ice_fract,di_ncat,ice_fract_ncat,k_sice,u_0,v_0,                &

! IN cloud data :
     q(1:row_length,1:rows,1:bl_levels),                                &
     qcf(1:row_length,1:rows,1:bl_levels),                              &
     qcl(1:row_length,1:rows,1:bl_levels),                              &
     qcf_latest,qcl_latest, t,                                          &

! IN everything not covered so far :
     rho_wet_theta,pstar,timestep,l_sice_heatflux,                      &
     l_sice_multilayers,                                                &

! IN variables from BDY_LAYR (that used to be local arrays)
     alpha1,ashtf,bq_gb,bt_gb,dtrdz_charney_grid,rdz_charney_grid,      &
     dtrdz_u,dtrdz_v,rdz_u,rdz_v,                                       &
     fraca,rhokh_tile,smc,                                              &
     chr1p5m,resfs,z0hssi,z0mssi,cdr10m_u,cdr10m_v,z1_tq,zh,            &
!ajm extra variable added
     rhokm_u,rhokm_v,                                                   &
! needed for new BL numerical solver
     bl_type_1,bl_type_2,bl_type_3,bl_type_4,                           &
     bl_type_5,bl_type_6,bl_type_7,                                     &

! IN additional variables for MOSES II
     tile_pts,tile_index,                                               &
     canhc_tile,flake,wt_ext_tile,                                      &
     lw_down,sw_tile,alpha1_sice,ashtf_tile,                            &
     fqw_ice,ftl_ice,resft,rhokh_sice,rhokpm,rhokpm_pot,rhokpm_sice,    &
     z0h_tile,z0m_tile,chr1p5m_sice,                                    &
     flandg_u,flandg_v,                                                 &

! IN additional variables used to calculate atmospheric cooling
!    at the screen level
     l_co2_interactive, co2_mmr, co2_3d,                                &
     rho1, f3_at_p, uStarGBM,                                           &

! IN additional variables for JULES
     rhokh_mix_dummy,dtstar_tile,dtstar,hcons,emis_tile,emis_soil,      &

! INOUT data :
!(Note ti and ti_gb are IN only if l_sice_multilayers=T)
     t_soil,ti,ti_gb,tstar,                                             &
     tstar_land,tstar_sea,tstar_sice,tstar_sice_cat,tstar_ssi,          &
     tstar_tile, t_latest,q_latest,                                     &

! INOUT Diagnostics started in BDY_LAYR not requiring STASH flags :
     e_sea,fqw,fqw_tile,epot_tile,ftl,ftl_tile,h_sea,rhokh,             &
     taux,tauy,                                                         &
     taux_land,taux_ssi,tauy_land,tauy_ssi,                             &
     TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans,                            &

! INOUT additional variables for MOSES II
     snow_tile,le_tile,radnet_sice,radnet_tile,olr,                     &

! OUT Increments to U and V momentum fields and Tl qw
!  (New dynamics only).
     du,dv,                                                             &

! OUT Diagnostic not requiring STASH flags :
     rhokh_mix,sea_ice_htf,                                             &

! OUT diagnostics requiring STASH flags :
     sice_mlt_htf,snomlt_surf_htf,latent_heat,                          &
     q1p5m,t1p5m,u10m,v10m,                                             &

! (IN) STASH flags :-
     simlt,smlt,slh,sq1p5,st1p5,su10,sv10,                              &
     l_ftl, l_fqw, l_taux, l_tauy,                                      &

! OUT additional variables for MOSES II
     esoil_tile,surf_ht_flux,surf_ht_flux_land,surf_ht_flux_sice,       &
     surf_htf_tile,ei_tile,                                             &
     q1p5m_tile,t1p5m_tile,ecan_tile,melt_tile,                         &
     e_ssi,ei_sice,ftl_ssi,                                             &
     error,                                                             &

! OUT data required elsewhere in UM system :
     ecan,ei,es,ext,snowmelt,                                           &
     BL_diag,lq_mix_bl,                                                 &
! IN SCM namelist data
     l_flux_bc,                                                         &
     ltimer                                                             &
     )

  IF (lhook) CALL dr_hook('IMPS_INTCT',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE imps_intct
