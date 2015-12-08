! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Longwave interface to the radiance code.
!
! Purpose:
!   This routine prepares the call to the ES radiance
!   scheme in the longwave.
!
! Method:
!   Principally, this routine transfers arrays into the correct formats.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
! Description of code:
!   FORTRAN 90
!
!- ---------------------------------------------------------------------
      subroutine r2_lwrad3z(ierr                                        &
!                       Gaseous mixing ratios
     &  , h2o, co2, o3                                                  &
     &  , co2_dim1, co2_dim2, co2_3d, l_co2_3d                          &
! chemical greenhouse gas fields
     &  , ngrgas, grgas_field                                           &
     &  , n2o_mix_ratio, ch4_mix_ratio                                  &
     &  , cfc11_mix_ratio, cfc12_mix_ratio, cfc113_mix_ratio            &
     &  , cfc114_mix_ratio                                              &
     &  , hcfc22_mix_ratio, hfc125_mix_ratio, hfc134a_mix_ratio         &
!                       Thermodynamic variables
        , tac, t_rad_surf, t_rad_land, t_rad_sice, t_rad_sea            &
        , l_ctile, pstar                                                &
     &  , p_layer_boundaries                                            &
     &  , p_layer_centres                                               &
     &  , height_theta                                                  &
     &  , height_rho                                                    &
!                       Options for COSP
        , L_cosp                                                        &
!                       Options for treating clouds
     &  , global_cloud_top, l_inhom_cloud, inhom_cloud                  &
     &  , dp_corr_strat, dp_corr_conv                                   &
!                       Stratiform cloud fields
        , l_pc2, lca_area, lca_bulk, lccwc1, lccwc2, n_drop_pot         &
!                       Convective cloud fields
     &  , cca, cccwp, ccw, lcbase, ccb, cct                             &
!                       Surface fields
     &  , land, flandg, ice_fraction                                    &
        , lying_snow, emis_land                                         &
!                       Solar fields
     &  , coszin, lit, scs                                              &
!                       Aerosol fields
     &  , l_climat_aerosol, l_clim_aero_hgt, L_HadGEM1_Clim_Aero        &
     &  , bl_depth, n_levels_bl                                         &
     &  , l_dust, l_use_dust, dust_dim1, dust_dim2                      &
     &  , dust_1, dust_2, dust_3, dust_4, dust_5, dust_6                &
     &  , l_use_biogenic, biogenic_dim1, biogenic_dim2, biogenic        &
     &  , l_sulpc_so2, l_use_sulpc_direct, l_use_sulpc_indirect         &
     &  , sulp_dim1,sulp_dim2                                           &
     &  , accum_sulphate, aitken_sulphate, diss_sulphate                &
     &  , sea_salt_film, sea_salt_jet, l_use_seasalt_indirect           &
     &  , l_use_seasalt_direct, salt_dim_a, salt_dim_b                  &
     &  , l_soot, l_use_soot_direct, soot_dim1, soot_dim2               &
     &  , fresh_soot, aged_soot                                         &
     &  , l_biomass, l_use_bmass_direct, bmass_dim1, bmass_dim2         &
     &  , fresh_bmass, aged_bmass, cloud_bmass, l_use_bmass_indirect    &
     &  , l_ocff, l_use_ocff_direct, ocff_dim1, ocff_dim2               &
     &  , fresh_ocff, aged_ocff, cloud_ocff, l_use_ocff_indirect        &
     &  , l_nitrate, l_use_nitrate_direct, nitrate_dim1, nitrate_dim2   &
     &  , accum_nitrate, diss_nitrate, l_use_nitrate_indirect           &
     &  , l_use_arcl, arcl_dim1, arcl_dim2, n_arcl_species              &
     &  , n_arcl_compnts, i_arcl_compnts, arcl                          &
     &  , aero_meso, l_murk_rad, Ntot_land, Ntot_sea                    &
     &  , l_ukca_radaer, ukca_radaer, ukca_dim1, ukca_dim2              &
     &  , ukca_mmr, ukca_cvl, ukca_dry, ukca_wet                        &
     &  , ukca_rho, ukca_vol, ukca_wtv, ukca_nbr                        & 
!                       Time
     &   , PREVIOUS_TIME                                                &
!                       Grid-dependent arrays
     &  , true_latitude                                                 &
!                       Level of tropopause
     &  , trindx                                                        &
!                       Spectrum
     &  , lw_spectrum                                                   &
!                       Algorithmic options
     &  , lw_control                                                    &
        , pts, l_mod_k_flux, l_scale_inc, list                          &
!                       Satellite viewing geometry
     &  , n_viewing_direction                                           &
     &  , viewing_direction1, viewing_direction2                        &
     &  , n_viewing_level, viewing_level                                &
!                       Diagnostics
     &  , n_channel, map_channel                                        &
     &  , LW_diag, row_list, col_list                                   &
!                       Physical dimensions
     &  , n_points, nlevs, n_layer, nclds                               &
     &  , nwet, nozone, row_length, rows, nd_field                      &
     &  , nd_field_flux_diag, nd_field_rad_diag                         &
     &  , nd_profile, nd_layer, nd_column, n_cca_lev, nd_channel        &
     &  , nd_flux_profile, nd_radiance_profile                          &
     &  , nd_viewing_level, nd_direction                                &
     &  , nd_cloud_component, nd_cloud_type                             &
     &  , nd_brdf_basis_fnc, nd_brdf_trunc                              &
     &  , nd_point_tile, nd_tile, id_ct, n_ukca_mode, n_ukca_cpnt       &
!                       Output fields
     &  , olr, lw_down, top_absorption, lwsea, lwout                    &
           ! Variables needed to calculate layer masses
     &   , rho_r2, r_rho_levels, r_theta_levels                         &
     &   , q, qcl, qcf, qcf2, qrain, qgraup                             &
     &   , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio        &
!           COSP input arguments
         , cosp_gbx                                                     &
     &   )
!
!
!
!     Modules included.
      USE rad_pcf
      USE conversions_mod, ONLY: recip_pi_over_180
      use dec_spec
      use control_struc
      use tileid3z
      use lwrdiag_mod, only: strlwdiag
      use mcica_mod
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ukca_radaer_struct_mod
      USE ereport_mod, ONLY : ereport
      USE arcl_mod, ONLY: npd_arcl_species, npd_arcl_compnts
      USE cosp_constants_mod, ONLY: i_cvcliq, i_cvcice, i_lscliq, i_lscice
      USE cosp_types_mod, ONLY: cosp_gridbox
      USE rad_input_mod, ONLY: sc
      IMPLICIT NONE

! Solar constant

!
!     Dummy arguments
!
      integer                                                           &
                !, intent(out)
     &    ierr
!           Error flag
!
!     Dimensions of arrays:
      integer, intent(in) :: row_length
!                              length of rows on each domain
      integer, intent(in) :: rows
!                              number of rows in the domain
      integer                                                           &
                !, intent(in)
     &    nd_field                                                      &
!           Field size in calling program
     &  , nd_field_flux_diag                                            &
!           Field size for flux diagnostics
     &  , nd_field_rad_diag                                             &
!           Field size for radiance diagnostics
     &  , nd_profile                                                    &
!           Size of array of grid-points
     &  , nd_layer                                                      &
!           Array sizes for layers
     &  , nd_column
!           Number of columns per point
      integer                                                           &
                !, intent(in)
     &    nd_channel                                                    &
!           Size allocated for spectral channels
     &  , nd_flux_profile                                               &
!           Size allocated for grid-points where fluxes are
!           diagnosed
     &  , nd_radiance_profile                                           &
!           Number of points where radiances are required
     &  , nd_viewing_level                                              &
!           Number of levels where radiances are required
     &  , nd_direction                                                  &
!           Number of directions in which radiances are calculated
     &  , nd_cloud_component                                            &
!           Number of components permitted in clouds
     &  , nd_cloud_type                                                 &
!           Number of permitted types of cloud
     &  , nd_brdf_basis_fnc                                             &
!           Maximum permitted number of BRDF basis functions
     &  , nd_brdf_trunc                                                 &
!           Maximum permitted order of truncation for BRDFs
     &  , nd_point_tile                                                 &
!           Size allocated for points where the surface is tiled
     &  , nd_tile
!           Size allocated for surface tiles
!
!     Actual sizes used:
      integer                                                           &
                !, intent(in)
     &    n_points                                                      &
!           Number of points
     &  , nwet                                                          &
!           Number of wet levels
     &  , nozone                                                        &
!           Number of levels with ozone
     &  , nlevs                                                         &
!           Number of layers in the main model
     &  , n_layer                                                       &
!           Number of layers seen in the radiation scheme
     &  , nclds                                                         &
!           Number of cloudy levels
     &  , n_levels_bl                                                   &
!           Number of layers occupied by boundary-layer aerosol
!           if l_clim_aero_hgt is false.
     &  , n_cca_lev                                                     &
     &  , n_ukca_mode                                                   &
!           Number of aerosol modes in UKCA_RADAER
     &  , n_ukca_cpnt
!           Number of aerosol components in UKCA_RADAER
!
!     Spectral data:
      type (spectrum) lw_spectrum
!
!     Control data:
      type (control_option) lw_control
!
      logical                                                           &
     &    l_scale_inc
!           Flag for scaling of heating rates to increments
!
!     Gaseous mixing ratios:
      real                                                              &
                !, intent(in)
     &    h2o(nd_field, nwet)                                           &
!           Mass mixing ratio of water
     &  , co2                                                           &
!           Mass mixing ratio of CO2
     &  , o3(nd_field, nozone)                                          &
!           Mass mixing ratios of ozone
     &  , n2o_mix_ratio                                                 &
!           Mass mixing ratio of nitrous oxide
     &  , ch4_mix_ratio                                                 &
!           Mass mixing ratio of methane
     &  , cfc11_mix_ratio                                               &
!           Mass mixing ratio of CFC11
     &  , cfc12_mix_ratio                                               &
!           Mass mixing ratio of CFC12
     &  , cfc113_mix_ratio                                              &
!           Mass mixing ratio of CFC113
     &  , cfc114_mix_ratio                                              &
!           Mass mixing ratio of CFC114
     &  , hcfc22_mix_ratio                                              &
!           Mass mixing ratio of HCFC22
     &  , hfc125_mix_ratio                                              &
!           Mass mixing ratio of HFC125
     &  , hfc134a_mix_ratio
!           Mass mixing ratio of HFC134a
      integer, intent(IN) :: ngrgas
      real, intent(IN) :: grgas_field(nd_field, nlevs, ngrgas)
!
!     General atmospheric properties:
      real                                                              &
                !, intent(in)
     &    p_layer_boundaries(nd_field,0:nlevs)                          &
!             pressure at boundaries of layers
     &  , p_layer_centres(nd_field,0:nlevs)                             &
!             pressure at centres of layers
     &  , height_theta(nd_field, 0:nlevs)                               &
     &  , height_rho(nd_field, nlevs)                                   &
     &  , tac(nd_field, nlevs)
!           Temperatures at centres of layers

!     Incident solar radiation:
      REAL, INTENT(IN) ::                                               &
     &    coszin(nd_field)                                              &
!           Cosines of zenith angle
     &  , scs                                                           &
!           Scaling of solar incident field
     &  , lit(nd_field)
!           Fraction of time point is lit

      Logical, intent(in)::                                             &
     &     l_mcr_qcf2                                                   &
                          ! Use second ice category
     &,    l_mcr_qrain                                                  &
                          ! Use prognostic rain
     &,    l_mcr_qgraup                                                 &
                          ! Use graupel
     &,    l_mixing_ratio ! Use mixing ratios in layer mass calculation
!     Flag for COSP
      LOGICAL,INTENT(IN) :: L_cosp
!
!     Options for treating clouds
      integer                                                           &
                !, intent(in)
     &    global_cloud_top
!           Global topmost cloudy layer
      logical                                                           &
                !, intent(in)
     &    l_inhom_cloud
!           Flag to use scaling factors for inhomogeneous cloud
      real                                                              &
                !, intent(in)
     &    inhom_cloud(nd_cloud_component)                               &
!           Scaling factors for inhomogeneous cloud
     &  , dp_corr_strat                                                 &
!           Decorrelation pressure scale for large scale cloud
     &  , dp_corr_conv
!           Decorrelation pressure scale for convective cloud
!
!     Properties of stratiform clouds:
      logical                                                           &
                !, intent(in)
     &    l_pc2
!           Flag to use PC2 cloud scheme
      real                                                              &
                !, intent(in)
     &    lccwc1(nd_field, nclds+1/(nclds+1))                           &
!           Liquid water contents (these are not used directly in
!           the radiation: the total condensed water content is
!           repartitioned using focwwil).
     &  , lccwc2(nd_field, nclds+1/(nclds+1))                           &
!           Ice water contents (these are not used directly in
!           The radiation: the total condensed water content is
!           Repartitioned using focwwil).
     &  , lca_area(nd_field, nclds+1/(nclds+1))                         &
!           Area fractions of layer clouds outside convective towers
     &  , lca_bulk(nd_field, nclds+1/(nclds+1))
!           Bulk fractions of layer clouds outside convective towers

      REAL, INTENT(IN) :: n_drop_pot(nd_field, nclds)

!     Properties of convective clouds:
      integer                                                           &
                !, intent(in)
     &    ccb(nd_field)                                                 &
!           Base of convective cloud
     &  , lcbase(nd_field)                                              &
!           Base of convective cloud (corrected)
     &  , cct(nd_field)
!           Top of convective cloud
      real                                                              &
                !, intent(in)
     &    cccwp(nd_field)                                               &
!           Water path of convective cloud
     &  , ccw(nd_field, nwet)                                           &
!           Convective cloud water
     &  , cca(nd_field,n_cca_lev)
!           Fraction of grid-box covered by convective cloud

!     Aerosols:
      logical                                                           &
                !, intent(in)
     &    l_climat_aerosol                                              &
!           Flag for climatological aerosol
     &  , l_clim_aero_hgt                                               &
!           Flag to use the depth of the boundary layer to set
!           the climatological aerosol
     &   , L_HadGEM1_Clim_Aero                                          &
!           Flag to use HadGEM1 setting for climatological aerosols
     &  , l_murk_rad
!           Flag for mesoscale model aerosol
      logical                                                           &
                !, intent(in)
     &    l_sulpc_so2                                                   &
!           Sulphur cycle available for effecting fluxes or diagnostics
     &  , l_use_sulpc_direct                                            &
!           Flag to use sulphur cycle for direct effect
     &  , l_use_sulpc_indirect                                          &
!           Flag to use sulphur cycle for indirect effect
     &  , l_dust                                                        &
!           Dust is available for effecting fluxes or diagnostics
     &  , l_use_dust                                                    &
!           Flag to use direct rad effect of mineral dust
     &  , l_use_biogenic                                                &
!           Flag to use biogenic for direct effect
     &  , l_soot                                                        &
!            Soot is available for effecting fluxes or diagnostics
     &  , l_use_soot_direct                                             &
!           Use direct rad. effect of soot aerosol
     &  , l_biomass                                                     &
!            Biomass is available for effecting fluxes or diagnostics
     &  , l_use_bmass_direct                                            &
!           Flag to use direct rad. effect of biomass smoke
     &  , l_use_bmass_indirect                                          &
!           Flag to use indirect effect of biomass smoke
     &  , l_ocff                                                        &
!            OCFF is available for effecting fluxes or diagnostics
     &  , l_use_ocff_direct                                             &
!           Flag to use direct rad. effect of ocff
     &  , l_use_ocff_indirect                                           &
!           Flag to use indirect effect of ocff
     &  , l_use_seasalt_indirect                                        &
!           Flag to use sea-salt for indirect effect
     &  , l_use_seasalt_direct                                          &
!           Flag to use sea-salt for direct effect
     &  , l_nitrate                                                     &
!            Nitrate is available for effecting fluxes or diagnostics
     &  , l_use_nitrate_direct                                          &
!           Flag to use nitrate for direct effect
     &  , l_use_nitrate_indirect
!           Flag to use nitrate for indirect effect

      integer                                                           &
                !,intent (in)
     &    sulp_dim1,sulp_dim2                                           &
!           Dimensions for _sulphate arrays, (P_FIELD,P_LEVELS or 1,1)
     &  , dust_dim1, dust_dim2                                          &
!           Dimensions for mineral dust arrays (p_field,p_levels or 1,1)
     &  , biogenic_dim1, biogenic_dim2                                  &
!           dimensions for biogenic array passed down to
!           r2_set_aerosol_field if direct effect required.
     &  , soot_dim1, soot_dim2                                          &
!           dimensions for soot arrays (P_FIELD,P_LEVELS or 1,1)
     &  , bmass_dim1, bmass_dim2                                        &
!           Dimensions for biomass arrays (P_FIELD,P_LEVELS or 1,1)
     &  , ocff_dim1, ocff_dim2                                          &
!           Dimensions for ocff arrays (P_FIELD,P_LEVELS or 1,1)
     &  , nitrate_dim1, nitrate_dim2                                    &
!           Dimensions for nitrate arrays (P_FIELD,P_LEVELS or 1,1)
     &  , salt_dim_a, salt_dim_b                                        &
!           dimensions for salt arrays on input (salt_dim_a=p_field
!           and salt_dim_b=p_levels, or else 1,1)
     &  , salt_dim_ind_a, salt_dim_ind_b                                &
!           dimensions for sea-salt arrays passed down to
!           r2_set_cloud_field if indirect effect required.
     &  , salt_dim_dir_a, salt_dim_dir_b
!           dimensions for sea-salt arrays passed down to
!           r2_set_aerosol_field if direct effect required.
      real                                                              &
                !, intent(in)
     &    accum_sulphate(sulp_dim1, sulp_dim2)                          &
!           Mass mixing ratio of accumulation mode aerosol
     &  , aitken_sulphate(sulp_dim1, sulp_dim2)                         &
!           Mass mixing ratio of aitken mode aerosol
     &  , diss_sulphate(sulp_dim1, sulp_dim2)                           &
!           Mixing ratio of dissolved sulphate
     &  , dust_1(dust_dim1, dust_dim2)                                  &
!           Mass mixing ratio of div1 dust
     &  , dust_2(dust_dim1, dust_dim2)                                  &
!           Mass mixing ratio of div2 dust
     &  , dust_3(dust_dim1, dust_dim2)                                  &
!           Mass mixing ratio of div3 dust
     &  , dust_4(dust_dim1, dust_dim2)                                  &
!           Mass mixing ratio of div4 dust
     &  , dust_5(dust_dim1, dust_dim2)                                  &
!           Mass mixing ratio of div5 dust
     &  , dust_6(dust_dim1, dust_dim2)                                  &
!           Mass mixing ratio of div6 dust
     &  , biogenic(biogenic_dim1, biogenic_dim2)                        &
!           Mixing ratios of biogenic aerosol
     &  , sea_salt_film(salt_dim_a, salt_dim_b)                         &
!           Number concentration of film-mode sea-salt aerosol
     &  , sea_salt_jet(salt_dim_a, salt_dim_b)                          &
!           Number concentration of jet-mode sea-salt aerosol
     &  , fresh_soot(soot_dim1, soot_dim2)                              &
!           Mixing ratios of fresh soot
     &  , aged_soot(soot_dim1, soot_dim2)                               &
!           Mixing ratios of aged soot
     &  , fresh_bmass(bmass_dim1, bmass_dim2)                           &
!           Mass mixing ratio of fresh biomass smoke
     &  , aged_bmass(bmass_dim1, bmass_dim2)                            &
!           Mass mixing ratio of aged biomass smoke
     &  , cloud_bmass(bmass_dim1, bmass_dim2)                           &
!           Mass mixing ratio of in-cloud biomass smoke
     &  , fresh_ocff(ocff_dim1, ocff_dim2)                              &
!           Mass mixing ratio of fresh fossil-fuel organic carbon aer
     &  , aged_ocff(ocff_dim1, ocff_dim2)                               &
!           Mass mixing ratio of aged fossil-fuel organic carbon aer
     &  , cloud_ocff(ocff_dim1, ocff_dim2)                              &
!           Mass mixing ratio of in-cloud fossil-fuel org carbon aer
     &  , accum_nitrate(nitrate_dim1, nitrate_dim2)                     &
!           Mass mixing ratio of accumulation nitrate aerosol
     &  , diss_nitrate(nitrate_dim1, nitrate_dim2)                      &
!           Mass mixing ratio of dissolved nitrate aerosol
     &  , aero_meso(nd_field, nlevs)
!           Mixing ratio of 'urban' aerosol of mesoscale model
!
!     Aerosol climatology for NWP

      ! Number of requested species within the climatology
      integer n_arcl_species
      
      ! Corresponding number of requested components
      integer n_arcl_compnts
      
      ! Model switch for each species
      logical l_use_arcl(NPD_ARCL_SPECIES)
      
      ! Index of each component
      integer i_arcl_compnts(NPD_ARCL_COMPNTS)
      
      ! Array dimensions
      integer                                                           &
     &        arcl_dim1                                                 &
     &   ,    arcl_dim2
     
      ! Mass-mixing ratios 
      real                                                              &
     &        arcl(arcl_dim1, arcl_dim2, n_arcl_compnts)
!
!     UKCA_RADAER: 
 
      ! Model switch 
      LOGICAL l_ukca_radaer 
       
      ! UKCA_RADAER structure 
      TYPE (ukca_radaer_struct) :: ukca_radaer 
       
      ! Dimensions 
      INTEGER ukca_dim1 
      INTEGER ukca_dim2 
       
      ! Component mass mixing ratios and volumes 
      REAL ukca_mmr(ukca_dim1, ukca_dim2, n_ukca_cpnt) 
      REAL ukca_cvl(ukca_dim1, ukca_dim2, n_ukca_cpnt) 
       
      ! Dry and wet modal diameters 
      REAL ukca_dry(ukca_dim1, ukca_dim2, n_ukca_mode) 
      REAL ukca_wet(ukca_dim1, ukca_dim2, n_ukca_mode) 
       
      ! Modal densities, volumes, volume of water, and number conc 
      REAL ukca_rho(ukca_dim1, ukca_dim2, n_ukca_mode) 
      REAL ukca_vol(ukca_dim1, ukca_dim2, n_ukca_mode) 
      REAL ukca_wtv(ukca_dim1, ukca_dim2, n_ukca_mode) 
      REAL ukca_nbr(ukca_dim1, ukca_dim2, n_ukca_mode) 
!
!     Carbon cycle:
      logical   l_co2_3d    !controls use of 3d co2 field
      integer                                                           &
                !, intent(in)
     &    co2_dim1, co2_dim2
!           Dimensions for CO2 array, (P_FIELD,P_LEVELS or 1,1)
      real                                                              &
                !, intent(in)
     &    co2_3d(co2_dim1, co2_dim2)
!           Mass mixing ratio of carbon dioxide
!     Surface fields:
      logical                                                           &
                !, intent(in)
     &    land(nd_field)                                                &
!           Land mask (true if land fraction >0.5)
     &  , l_ctile
!           coastal tiling switch
      real                                                              &
                !, intent(in)
     &    flandg(nd_field)
!           land fraction in grid box
      real                                                              &
                !, intent(in)
     &    pstar(nd_field)                                               &
!           Surface pressures
        , t_rad_surf(nd_field)                                          &
!           Effective radiative temperature over whole grid-box
        , t_rad_land(nd_field)                                          &
!           Effective radiative temperature over land
        , t_rad_sice(nd_field)                                          &
!           Effective radiative temperature over sea-ice
        , t_rad_sea(nd_field)                                           &
!           Radiative temperature over open sea
     &  , ice_fraction(nd_field)                                        &
!           Sea ice fraction of sea portion of grid box
     &  , lying_snow(nd_field)                                          &
!           Mass loading of lying snow
        , emis_land(nd_field)                                           &    
!           Mean land emissivity in a gridbox if l_um_jules is true 
     &  , bl_depth(nd_field)
!           depth of the boundary layer
!
      real                                                              &
                !, intent(in)
     &     Ntot_land                                                    &
!           number of droplets over land / m-3
     &   , Ntot_sea
!           number of droplets over sea / m-3
!
      REAL, INTENT(IN) :: true_latitude(nd_field)
      INTEGER, INTENT(IN) :: PREVIOUS_TIME(7)

!                       Level of tropopause
      integer, intent(in) ::                                            &
     &    trindx(nd_field)
!           The layer boundary of the tropopause
!
!     increment of time:
      real                                                              &
                !, intent(in)
     &    pts
!           Time increment
!
!     Use modulus of fluxes to remove negative effective extinctions
      logical, intent(in) :: l_mod_k_flux

      integer, intent(in) :: list(nd_field)
!                              list of points where radiation is to be
!                              calculated
!
!     Satellite viewing geometry
      INTEGER, Intent(IN) :: n_channel
!           Number of channels calculated simultaneously in one
!           call to the radiation code
      INTEGER, Intent(IN) :: map_channel(lw_spectrum%npd_band)
!           Mapping of bands in the spectral file to channels in the
!           diagnostic output
      INTEGER, Intent(IN) :: n_viewing_direction
!           Number of viewing directions
      INTEGER, Intent(IN) :: n_viewing_level
!           Number of levels where the radiance is calculated
      REAL, Intent(IN) :: viewing_direction1(nd_field, nd_direction)
!           Satellite viewing direction (Zenith Angle)
      REAL, Intent(IN) :: viewing_direction2(nd_field, nd_direction)
!           Satellite viewing directions (Azimuthal Angle)
      REAL, Intent(IN) :: viewing_level(nd_viewing_level)
!           Levels where radiances are calculated
!
!
!     Information for the calculation of layer masses
      Real, intent(in)::                                                &
     &  rho_r2(nd_field,nlevs)                                          &
                                ! Air density*radius of earth**2 / kg m-1
     &, r_rho_levels(nd_field,nlevs)                                    &
                                      ! Height of rho levels / m
     &, r_theta_levels(nd_field,0:nlevs)                                &
                                           ! Height of theta levels / m
     &, q(nd_field,nwet)                                                &
                                ! Water vapour mixing ratio / kg kg-1
     &, qcl(nd_field,nwet)                                              &
                                ! Liquid water mixing ratio / kg kg-1
     &, qcf(nd_field,nwet)                                              &
                                ! Ice mixing ratio / kg kg-1
     &, qcf2(nd_field,nwet)                                             &
                                ! Second ice category mr / kg kg-1
     &, qrain(nd_field,nwet)                                            &
                                ! Rain mixing ratio / kg kg-1
     &, qgraup(nd_field,nwet)  ! Graupel mixing ratio / kg kg-1

!     Calculated fluxes:
      real                                                              &
                !, intent(out)
     &    olr(nd_field)                                                 &
!           Net outgoing radiation
     &  , lw_down(nd_field)                                             &
!           Downwards surface flux
     &  , top_absorption(nd_field)                                      &
!           Absorption in the extra radiative layer at the top
!           of the model
     &  , lwout(nd_field, nlevs+1)                                      &
!           Net downward fluxes or heating rates
     &  , lwsea(nd_field)
!           Sea-surface components of flux
!
!
!
!     Diagnostics:
      type (strlwdiag) :: lw_diag
!
      integer, intent(in) :: row_list(nd_field)
!                              list of row indices of points
!                              to be treated
      integer, intent(in) :: col_list(nd_field)
!                              list of column indices of points
!                              to be treated
!
!     Structure with COSP arguments
      TYPE(cosp_gridbox),INTENT(OUT) :: cosp_gbx
!
!
!     Local variables.
!
!     Locally allocated dimensions
      integer                                                           &
     &    nd_source_coeff                                               &
!           Size allocated for source coefficients
     &  , nd_layer_clr                                                  &
!           Size allocated for totally clear layers
     &  , nd_2sg_profile                                                &
!           Size allocated for profiles in two-stream or Gaussian
!           runs
     &  , nd_region                                                     &
!           Size allocated for cloudy regions
     &  , nd_overlap_coeff                                              &
!           Size allocated for cloudy overlap coefficients
     &  , nd_max_order                                                  &
!           Size allocated for orders of spherical harmonics
     &  , nd_sph_coeff                                                  &
!           Size allocated for coefficients of spherical harmonics
     &  , id_ct
!           Top level in arrays of cloud properties
!
      integer                                                           &
     &    i                                                             &
!           Loop variable
     &  , j                                                             &
!           Loop variable
     &  , k                                                             &
!           Loop variable
     &  , l                                                             &
!           Loop variable
     &  , ll                                                            &
!           Loop variable
     &  , lll                                                           &
!           Loop variable
     &  , ic
!           Loop variable
      logical                                                           &
     &    l_clear
!           Calculate clear-sky fields
!     Flags for processes actually enabled.
      integer                                                           &
     &    i_solver_clear                                                &
!           Solver for clear-sky fluxes
     &  , i_gas_overlap(lw_spectrum%npd_band)                           &
!           Overlaps in each band
     &  , n_order_forward
!           Order of term used to define the
!           forward scattering fraction
!
!     General atmospheric properties:
      real                                                              &
     &    d_mass(nd_profile, nd_layer)                                  &
!           Mass thicknesses of layers
     &  , p(nd_profile, nd_layer)                                       &
!           Pressure field
     &  , t(nd_profile, nd_layer)                                       &
!           Temperature field
     &  , t_bdy(nd_profile, 0: nd_layer)                                &
!           Temperature field at boundaries
     &  , gas_mix_ratio(nd_profile, nd_layer, lw_spectrum%npd_species)
!           Mass fractions of gases
!
      real :: layer_heat_capacity(nd_profile, nd_layer)
!           Specific heat capacity of layer * d_mass

      real :: alat(nd_profile)
!           Latitude in degrees

!     Surface fields:
      Logical :: Land_g(nd_profile)
!           Gathered land-surface mask
      Real :: ice_fraction_g(nd_profile)
!           Gathered ice fraction
      integer                                                           &
     &    n_brdf_basis_fnc
!           Number of BRDF basis functions
      real                                                              &
     &    rho_alb(nd_profile, nd_brdf_basis_fnc, lw_spectrum%npd_band)  &
!           Weights for BRDF basis functions
     &  , f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                  &
     &      , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)                     &
!           Array of BRDF basis terms
     &  , t_surface(nd_profile)                                         &
!           Gathered temperature of surface
     &  , flandg_g(nd_profile)
!           Gathered land fraction
!
!     Arrays related to tiling of the surface
      logical                                                           &
     &    l_rad_tile
!           Local to allow tiling within the radiation scheme: this is
!           used to handle all surface heterogeneities, whether or not
!           explicit coastal tiling is required (i.e. even when a grid
!           box may be only land or sea, this facility is still used to
!           represent the split between open sea and sea-ice.)
      integer                                                           &
     &    n_point_tile                                                  &
!           Number of points to tile
     &  , n_tile                                                        &
!           Number of tiles used
     &  , list_tile(nd_point_tile)                                      &
!           List of points with surface tiling indexed over the list
!           of points where radiances are calculated
     &  , list_tile_outer(nd_point_tile)                                &
!           List of points with surface tiling indexed over the full
!           list of points where this routine has been called
     &  , index_tile(npd_tile_type)
!           The indexing number of tiles of the given type
      real                                                              &
     &    rho_alb_tile(nd_point_tile, nd_brdf_basis_fnc                 &
     &      , nd_tile, lw_spectrum%npd_band)                            &
!           Weights for the basis functions of the BRDFs
!           at the tiled points
     &  , t_tile(nd_point_tile, nd_tile)                                &
!           Local surface temperatures on individual tiles
     &  , frac_tile(nd_point_tile, nd_tile)
!           Fractions of ech tiled grid-point occupied by tiles of
!           the relevant types
!
!     Cloudy properties:
      integer                                                           &
     &    n_condensed                                                   &
!           Number of condensed phases
     &  , type_condensed(nd_cloud_component)                            &
!           Types of condensed components
     &  , i_condensed_param(nd_cloud_component)                         &
!           Parametrization schemes for components
     &  , condensed_n_phf(nd_cloud_component)                           &
!           Number of terms in the phase functions
     &  , n_cloud_top_global                                            &
!           Inverted global topmost cloudy layer
     &  , n_cloud_type
!           Number of types of clouds
      real                                                              &
     &    condensed_param_list(lw_spectrum%npd_cloud_parameter          &
     &      , nd_cloud_component, lw_spectrum%npd_band)                 &
!           Parameters for condensed phases
     &  , condensed_dim_char(nd_profile, id_ct: nd_layer                &
     &      , nd_cloud_component)                                       &
!           Characteristic dimensions of condensed species
     &  , condensed_mix_ratio(nd_profile, id_ct: nd_layer               &
     &      , nd_cloud_component)                                       &
!           Mass fractions of liquid water
     &  , w_cloud(nd_profile, id_ct: nd_layer)                          &
!           Cloud amounts
     &  , frac_cloud(nd_profile, id_ct: nd_layer, nd_cloud_type)        &
!           Fractions of different types of cloud
     &  , tot_cloud_cover(nd_profile)                                   &
!           Total cloud cover
     &  , condensed_min_dim(nd_cloud_component)                         &
!           Minimum dimensions of condensed components
     &  , condensed_max_dim(nd_cloud_component)
!           Maximum dimensions of condensed components
!
!     Properties of aerosols:
      logical                                                           &
     &    l_aerosol                                                     &
!           Flag to enable direct aerosol effects within
!           the radiation code
     &  , l_aerosol_ccn
!           Flag to define CCN from aerosols
      real                                                              &
     &    aerosol_mix_ratio(nd_profile, nd_layer                        &
     &      , lw_spectrum%npd_aerosol_mixratio)
!           Mixing ratios of aerosols
      integer                                                           &
          aerosol_mr_type_index(lw_spectrum%npd_aerosol_mixratio)       &
!             Index relating aerosol_mix_ratio aerosols to aerosols in
!             the spectral information
        , aerosol_mr_source(lw_spectrum%npd_aerosol_mixratio)
!             Scheme/source of the aerosol data, to determine use in
!             changing radiative fluxes and use in diagnostics
!
!     Properties of ukca aerosols 
      REAL ukca_mix_ratio(nd_profile, nd_layer, n_ukca_cpnt) 
      REAL ukca_comp_vol(nd_profile, nd_layer, n_ukca_cpnt) 
      REAL ukca_dry_diam(nd_profile, nd_layer, n_ukca_mode) 
      REAL ukca_wet_diam(nd_profile, nd_layer, n_ukca_mode) 
      REAL ukca_modal_rho(nd_profile, nd_layer, n_ukca_mode) 
      REAL ukca_modal_vol(nd_profile, nd_layer, n_ukca_mode) 
      REAL ukca_modal_wtv(nd_profile, nd_layer, n_ukca_mode) 
      REAL ukca_modal_nbr(nd_profile, nd_layer, n_ukca_mode) 
!
!     Solar fields:
      REAL ::                                                           &
          zen_0(nd_profile)                                             &
!           Secants of zenith angle (two-stream) or cosines of
!           the zenith angle (spherical harmonics)
        , solar_irrad(nd_profile)                                       &
!           Normally incident solar irradiance
        , solar_tail_flux(nd_profile)
!
!     Local algorithmic variables:
      real                                                              &
     &    euler_factor
!           Weighting applied to the last term of alternating series
      integer                                                           &
     &    i_scatter_method_band(lw_spectrum%npd_band)
!           Method of treating scattering in the calculation
!           of optical propeties in each band
!
!     Gathered viewing directions:
      real                                                              &
     &    viewing_direction_g(nd_radiance_profile, nd_direction, 2)
!           Satellite viewing directions
!
!     Fluxes:
      real                                                              &
     &    flux_direct(nd_flux_profile, 0: nd_layer, nd_channel)         &
!           Direct flux
     &  , flux_direct_clear(nd_flux_profile, 0: nd_layer, nd_channel)   &
!           Clear-sky direct flux
     &  , flux_net(nd_flux_profile, 0: nd_layer, nd_channel)            &
!           Downward/net flux
     &  , flux_net_clear(nd_flux_profile, 0: nd_layer, nd_channel)      &
!           Clear-sky downward/net flux
     &  , flux_up(nd_flux_profile, 0: nd_layer, nd_channel)             &
!           Upward flux
     &  , flux_up_clear(nd_flux_profile, 0: nd_layer, nd_channel)
!           Clear-sky upward flux
!
!     Local fields for cloud diagnostics:
      real                                                              &
     &    cloud_absorptivity_g(nd_flux_profile, nd_layer)               &
!           Mean absorption coefficient in clouds weighted by the
!           cloud amount and the clear-sky flux
     &  , cloud_weight_absorptivity_g(nd_flux_profile, nd_layer)        &
!           Weighting factor for absorption in clouds: the product
!           of the cloud amount and the clear-sky flux
     &  , ls_cloud_absorptivity_g(nd_flux_profile, nd_layer)            &
!           Mean absorption coefficient in layer clouds weighted by the
!           cloud amount and the clear-sky flux
     &  , ls_cloud_weight_absorptivity_g(nd_flux_profile, nd_layer)     &
!           Weighting factor for absorption in layer clouds: the product
!           of the cloud amount and the clear-sky flux
     &  , cnv_cloud_absorptivity_g(nd_flux_profile, nd_layer)           &
!           Mean absorption coefficient in conv.clouds weighted by the
!           cloud amount and the clear-sky flux
     &  , cnv_cloud_weight_absorptivity_g(nd_flux_profile, nd_layer)
!           Weighting factor for absorption in conv.clouds: the product
!           of the cloud amount and the clear-sky flux

!     Switches for cloud diagnostics:
       LOGICAL :: l_ls_cloud_absorptivity, l_cnv_cloud_absorptivity

!
!     local arrays to fill diagnostics:
      real :: total_cloud_cover_g(nd_profile)
!             cloud fraction at gathered points
!     Tiled fluxes:
      real                                                              &
     &    flux_up_tile(nd_point_tile, nd_tile, nd_channel)
!           Upward fluxes at tiled surface points
!
!     Radiances:
      real                                                              &
     &    radiance(nd_radiance_profile, nd_viewing_level                &
     &      , nd_direction, nd_channel)
!           Radiances calculated
!
!     Fields required for call to radiation code but not used
      INTEGER, Parameter :: nd_j_profile=1
      INTEGER :: i_gas
!
!
!     Auxiliary variables:
      real                                                              &
     &    dacon                                                         &
!           Difference in A's
     &  , dbcon                                                         &
!           Difference in B's
     &  , weight_band(lw_spectrum%npd_band)                             &
!           Weighting factors for bands
     &  , nullmmr
!           Null mass mixing ratio
      parameter(nullmmr=0.0)
!
!     Dummy fields for radiation code
      logical                                                           &
     &     l_dummy
      real                                                              &
     &     dummy                                                        &
     &    ,dummy1d(1)                                                   &
     &    ,dummy2d(1,1)                                                 &
     &    ,dummy3d(1,1,1)
!
!     Local arrays for Aerosol Optical Depth diagnostics
      real                                                              &
     &     aod_sulphate(nd_profile, lw_spectrum%npd_aod_wavel),         &
     &     aod_dust(nd_profile, lw_spectrum%npd_aod_wavel),             &
     &     aod_seasalt(nd_profile, lw_spectrum%npd_aod_wavel),          &
     &     aod_soot(nd_profile, lw_spectrum%npd_aod_wavel),             &
     &     aod_biomass(nd_profile, lw_spectrum%npd_aod_wavel),          &
     &     aod_biogenic(nd_profile, lw_spectrum%npd_aod_wavel),         &
     &     aod_ocff(nd_profile, lw_spectrum%npd_aod_wavel),             &
     &     aod_delta(nd_profile, lw_spectrum%npd_aod_wavel),            &
     &     aod_nitrate(nd_profile, lw_spectrum%npd_aod_wavel),          &
     &     aod_total_radn(nd_profile, lw_spectrum%npd_aod_wavel),       &
     &     angst_total_radn(nd_profile, lw_spectrum%npd_aod_wavel),     &
     &     aod_prog_sulphate(nd_profile, lw_spectrum%npd_aod_wavel),    &
     &     aod_prog_dust(nd_profile, lw_spectrum%npd_aod_wavel),        &
     &     aod_prog_seasalt(nd_profile, lw_spectrum%npd_aod_wavel),     &
     &     aod_prog_soot(nd_profile, lw_spectrum%npd_aod_wavel),        &
     &     aod_prog_biomass(nd_profile, lw_spectrum%npd_aod_wavel),     &
     &     aod_prog_ocff(nd_profile, lw_spectrum%npd_aod_wavel),        &
     &     aod_prog_nitrate(nd_profile, lw_spectrum%npd_aod_wavel)
!
!     Local arrays for ukca aerosol optical depth diagnostics 
      REAL                                                              & 
           aod_ukca_ait_sol(nd_profile, lw_spectrum%npd_aod_wavel),     & 
           aod_ukca_acc_sol(nd_profile, lw_spectrum%npd_aod_wavel),     & 
           aod_ukca_cor_sol(nd_profile, lw_spectrum%npd_aod_wavel),     & 
           aod_ukca_ait_ins(nd_profile, lw_spectrum%npd_aod_wavel),     & 
           aod_ukca_acc_ins(nd_profile, lw_spectrum%npd_aod_wavel),     & 
           aod_ukca_cor_ins(nd_profile, lw_spectrum%npd_aod_wavel) 

      CHARACTER (LEN=*), PARAMETER :: RoutineName = 'r2_lwrad3z'
      CHARACTER (LEN=240) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
!
      IF (lhook) CALL dr_hook('R2_LWRAD3Z',zhook_in,zhook_handle)
!
!     Initialize the error flag for the radiation code.
      ierr=i_normal
      condensed_dim_char(:,:,:) = 0.0

!     Initialise logicals which are used later on.
      l_ls_cloud_absorptivity  = .FALSE.
      l_cnv_cloud_absorptivity = .FALSE.

!     Becuase the radiation code may be invoked on subsamples of
!     the same segment with different selections of points there
!     may be ocassions when this routine is called only to fill
!     diagnostics with zeros.
      if (n_points >  0) then
!
!       Set the logical flag for dummy diagnostics not available from
!       the lower code in the long-wave to .FALSE..
        l_dummy=.false.
!
!
!       Set clear-sky calculations.
        l_clear = lw_diag%l_flux_up_clear                               &
         .OR. lw_diag%l_flux_down_clear                                 &
         .OR. lw_diag%l_clear_olr                                       &
         .OR. lw_diag%l_surf_down_clr                                   &
         .OR. lw_diag%l_clear_hr                                        &
         .OR.(lw_diag%l_cloud_absorptivity.AND.                         &
              lw_diag%l_cloud_weight_absorptivity)                      &
         .OR.(lw_diag%l_ls_cloud_absorptivity.AND.                      &
              lw_diag%l_ls_cloud_weight_absorptivity)                   &
         .OR.(lw_diag%l_cnv_cloud_absorptivity.AND.                     &
              lw_diag%l_cnv_cloud_weight_absorptivity).OR.l_cosp

!
!       Set dynamic array sizes and dependent options.
! DEPENDS ON: r2_set_option
        call r2_set_option(ierr                                         &
     &    , n_layer, lw_spectrum                                        &
     &    , lw_control%isolir                                           &
     &    , lw_control%i_gas_overlap, i_gas_overlap                     &
     &    , lw_control%l_aerosol, l_climat_aerosol                      &
     &    , l_use_sulpc_direct, l_use_soot_direct, l_use_biogenic       &
     &    , l_use_dust, l_use_bmass_direct, l_use_ocff_direct           &
     &    , l_use_nitrate_direct                                        &
     &    , l_use_seasalt_direct, l_murk_rad, l_aerosol                 &
     &    , l_use_sulpc_indirect, l_aerosol_ccn, n_arcl_species         &
     &    , lw_control%l_global_cloud_top                               &
     &    , global_cloud_top, n_cloud_top_global                        &
     &    , l_clear                                                     &
     &    , lw_control%i_angular_integration                            &
     &    , lw_control%i_solver, i_solver_clear                         &
     &    , lw_control%l_rescale, n_order_forward                       &
     &    , lw_control%i_truncation, lw_control%ls_global_trunc         &
     &    , lw_control%l_euler_trnf, euler_factor                       &
     &    , lw_control%i_scatter_method, i_scatter_method_band          &
     &    , l_rad_tile                                                  &
     &    , weight_band                                                 &
     &    , nd_overlap_coeff, nd_2sg_profile, nd_layer_clr              &
     &    , nd_source_coeff, nd_max_order, nd_sph_coeff                 &
     &    , nd_region                                                   &
     &    , nd_profile                                                  &
     &    )
!
!       If radiances are to be calculated the viewing directions
!       are to be copied into the gathered array.
        if ( (lw_control%i_angular_integration ==                       &
     &        ip_spherical_harmonic).and.                               &
     &       (lw_control%i_sph_mode == ip_sph_mode_rad) ) then
          do ll=1, n_points
            l=list(ll)
            viewing_direction_g(ll, 1, 1)=viewing_direction1(l, 1)
            viewing_direction_g(ll, 1, 2)=viewing_direction2(l, 1)
          enddo
        endif
!
!
!
!       Set the mixing ratios of gases.
! DEPENDS ON: r2_set_gas_mix_ratio
        call r2_set_gas_mix_ratio(ierr                                  &
     &    , n_points, nlevs, n_layer, nwet, nozone                      &
     &    , list, lw_control%l_extra_top                                &
     &    , lw_spectrum%n_absorb, lw_spectrum%type_absorb               &
     &    , lw_control%l_n2o, lw_control%l_ch4, lw_control%l_cfc11      &
     &    , lw_control%l_cfc12, .false.                                 &
     &    , lw_control%l_cfc113, lw_control%l_cfc114                    &
     &    , lw_control%l_hcfc22                                         &
     &    , lw_control%l_hfc125, lw_control%l_hfc134a                   &
     &    , h2o, co2, o3, n2o_mix_ratio, ch4_mix_ratio                  &
     &    , cfc11_mix_ratio, cfc12_mix_ratio, nullmmr                   &
     &    , cfc113_mix_ratio, cfc114_mix_ratio, hcfc22_mix_ratio        &
     &    , hfc125_mix_ratio, hfc134a_mix_ratio                         &
     &    , gas_mix_ratio                                               &
     &    , co2_dim1, co2_dim2, co2_3d, l_co2_3d                        &
     &    , nd_field, nd_profile, nd_layer, lw_spectrum%npd_species, 1  &
! chemical greenhouse gas fields
     &    , ngrgas, grgas_field                                         &
     &    )
        if (ierr /= i_normal) then 
          cmessage = 'Error following call to r2_set_gas_mix_ratio'
          GO TO 9999
        end if
!
!
!       Calculate pressures and temperatures.
! DEPENDS ON: r2_set_thermodynamic
        call r2_set_thermodynamic(n_points,nlevs,n_layer,nwet,list      &
     &    , lw_control%l_extra_top, .true.                              &
          , pstar                                                       &
     &    , p_layer_boundaries                                          &
     &    , p_layer_centres                                             &
     &    , height_theta                                                &
     &    , height_rho                                                  &
     &    , tac                                                         &
     &    , rho_r2, r_rho_levels, r_theta_levels                        &
     &    , q, qcl, qcf, qcf2, qrain, qgraup                            &
     &    , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio       &
          , p, t, t_bdy, d_mass                                         &
     &    , layer_heat_capacity                                         &
     &    , nd_field, nd_profile, nd_layer                              &
     &    )
!
!
!       Gathering of land fields and surface temperatures
        DO l=1, n_points
          land_g(l)=land(list(l))
          flandg_g(l)=flandg(list(l))
          t_surface(l) = t_rad_surf(list(l))
        END DO
!
! DEPENDS ON: r2_set_surface_field_lw
        CALL r2_set_surface_field_lw(ierr                               &
          , n_points, list                                              &
          , lw_spectrum%n_band, lw_control%ls_brdf_trunc                &
          , l_ctile, flandg, emis_land                                  &
          , n_brdf_basis_fnc, f_brdf, rho_alb                           &
          , land, ice_fraction, t_rad_land, t_rad_sice, t_rad_sea       &
          , l_rad_tile, n_point_tile, n_tile, list_tile, rho_alb_tile   &
          , frac_tile, t_tile, list_tile_outer, index_tile              &
          , nd_field, nd_profile, lw_spectrum%npd_band                  &
          , nd_brdf_basis_fnc                                           &
          , nd_brdf_trunc, nd_point_tile, nd_tile                       &
          )
!
!
!       set sea-salt array dimensions.
        if (l_use_seasalt_direct) then
           salt_dim_dir_a=salt_dim_a
           salt_dim_dir_b=salt_dim_b
        else
           salt_dim_dir_a=1
           salt_dim_dir_b=1
        endif
!
      Do l=1, n_points
        alat(l) = Recip_Pi_Over_180 * true_latitude(list(l))
      Enddo      
        
!
!       Set the mixing ratios of aerosols.
        if (l_aerosol.or.l_aerosol_ccn) then
! DEPENDS ON: r2_set_aerosol_field
          call r2_set_aerosol_field(ierr                                &
     &      , n_points, nlevs, n_layer, lw_spectrum%n_aerosol           &
     &      , lw_spectrum%n_aerosol_mr, lw_spectrum%type_aerosol        &
     &      , list, lw_control%l_extra_top                              &
     &      , l_climat_aerosol, l_clim_aero_hgt, L_HadGEM1_Clim_Aero    &
     &      , bl_depth, t, n_levels_bl, l_murk_rad, aero_meso           &
     &      , l_dust, l_use_dust, dust_dim1, dust_dim2                  &
     &      , dust_1, dust_2, dust_3, dust_4, dust_5, dust_6            &
     &      , l_use_biogenic, biogenic_dim1, biogenic_dim2, biogenic    &
     &      , l_sulpc_so2, l_use_sulpc_direct, sulp_dim1, sulp_dim2     &
     &      , accum_sulphate, aitken_sulphate                           &
     &      , l_use_seasalt_direct, salt_dim_dir_a, salt_dim_dir_b      &
     &      , sea_salt_film, sea_salt_jet, p                            &
     &      , l_soot, l_use_soot_direct, soot_dim1, soot_dim2           &
     &      , fresh_soot, aged_soot                                     &
     &      , l_biomass, l_use_bmass_direct, bmass_dim1, bmass_dim2     &
     &      , fresh_bmass, aged_bmass                                   &
     &      , l_ocff, l_use_ocff_direct, ocff_dim1, ocff_dim2           &
     &      , fresh_ocff, aged_ocff                                     &
     &      , l_nitrate,l_use_nitrate_direct, nitrate_dim1, nitrate_dim2&
     &      , accum_nitrate                                             &
     &      , n_arcl_species, n_arcl_compnts, i_arcl_compnts            &
     &      , l_use_arcl, arcl_dim1, arcl_dim2, arcl                    &
     &      , land, lying_snow, pstar                                   &
     &      , p_layer_boundaries                                        &
     &      , trindx, alat, PREVIOUS_TIME                               &
     &      , aerosol_mix_ratio, aerosol_mr_source,aerosol_mr_type_index&
     &      , nd_field, nd_profile, nd_layer                            &
     &      , lw_spectrum%npd_aerosol_species                           &
     &      , lw_spectrum%npd_aerosol_mixratio, 1                       &
     &      )
           if (ierr /= i_normal) then 
             cmessage = 'Error following call to r2_set_aerosol_field'
             GO TO 9999
           end if
        endif
!
!     Set the mixing ratios of UKCA aerosols 
      IF (l_ukca_radaer) THEN 
! DEPENDS ON: ukca_radaer_set_aerosol_field 
         CALL ukca_radaer_set_aerosol_field(                            & 
              list, lw_control%l_extra_top, n_layer, n_points           & 
            , ukca_dim1, ukca_dim2                                      & 
            , ukca_mmr, ukca_cvl                                        & 
            , ukca_dry, ukca_wet, ukca_rho, ukca_vol, ukca_wtv          & 
            , ukca_nbr                                                  & 
            , ukca_mix_ratio, ukca_comp_vol                             & 
            , ukca_dry_diam, ukca_wet_diam, ukca_modal_rho              & 
            , ukca_modal_vol, ukca_modal_wtv                            & 
            , ukca_modal_nbr                                            & 
            , nd_field, nd_profile, nd_layer                            & 
            , n_ukca_cpnt, n_ukca_mode                                  & 
            ) 
      END IF 
! 
!
!
!
!       Assign the properties of clouds. A dummy array must be passed
!       for the microphysical diagnostics since they are not available
!       through STASH in the long-wave.
!
! DEPENDS ON: r2_set_cloud_parametrization
        call r2_set_cloud_parametrization(ierr, lw_spectrum%n_band      &
     &    , lw_control%i_st_water, lw_control%i_cnv_water               &
     &    , lw_control%i_st_ice, lw_control%i_cnv_ice                   &
     &    , lw_spectrum%l_drop_type                                     &
     &    , lw_spectrum%i_drop_parametrization                          &
     &    , lw_spectrum%n_drop_phf_term, lw_spectrum%drop_parameter_list&
     &    , lw_spectrum%drop_parm_min_dim, lw_spectrum%drop_parm_max_dim&
     &    , lw_spectrum%l_ice_type                                      &
     &    , lw_spectrum%i_ice_parametrization                           &
     &    , lw_spectrum%n_ice_phf_term,lw_spectrum%ice_parameter_list   &
     &    , lw_spectrum%ice_parm_min_dim, lw_spectrum%ice_parm_max_dim  &
     &    , i_condensed_param, condensed_n_phf, condensed_param_list    &
     &    , condensed_min_dim, condensed_max_dim                        &
     &    , lw_spectrum%npd_band, lw_spectrum%npd_drop_type             &
     &    , lw_spectrum%npd_ice_type, lw_spectrum%npd_cloud_parameter   &
     &    , nd_cloud_component                                          &
     &    )
        if (ierr /= i_normal) then 
          cmessage = 'Error following call to r2_set_cloud_parametrization'
          GO TO 9999
        end if
!
!
!       set sea-salt array dimensions.
        if (l_use_seasalt_indirect) then
           salt_dim_ind_a=salt_dim_a
           salt_dim_ind_b=salt_dim_b
        else
           salt_dim_ind_a=1
           salt_dim_ind_b=1
        endif
!
!
! DEPENDS ON: r2_set_cloud_field
        call r2_set_cloud_field(n_points, nlevs, n_layer, nclds         &
     &    , list                                                        &
     &    , p, t, d_mass, alat                                          &
     &    , ccb, cct, cca, cccwp, ccw, lcbase                           &
          , lccwc1, lccwc2, lca_area, lca_bulk, n_drop_pot              &
     &    , l_pc2, lw_control%l_microphysics, l_aerosol_ccn             &
     &    , sea_salt_film, sea_salt_jet                                 &
     &    , l_use_seasalt_indirect, salt_dim_ind_a, salt_dim_ind_b      &
     &    , l_use_biogenic, biogenic, biogenic_dim1, biogenic_dim2      &
     &    , sulp_dim1, sulp_dim2, accum_sulphate, diss_sulphate         &
     &    , aitken_sulphate, l_use_bmass_indirect                       &
     &    , bmass_dim1, bmass_dim2, aged_bmass, cloud_bmass             &
     &    , l_use_ocff_indirect, ocff_dim1, ocff_dim2                   &
     &    , aged_ocff, cloud_ocff                                       &
     &    , l_use_nitrate_indirect, nitrate_dim1, nitrate_dim2          &
     &    , accum_nitrate, diss_nitrate                                 &
     &    , lying_snow                                                  &
     &    , land_g, flandg_g                                            &
     &    , lw_control%i_cloud_representation, i_condensed_param        &
     &    , condensed_min_dim, condensed_max_dim                        &
     &    , n_condensed, type_condensed                                 &
     &    , w_cloud, n_cloud_type, frac_cloud                           &
     &    , lw_control%l_local_cnv_partition                            &
     &    , condensed_mix_ratio, condensed_dim_char                     &
!                       Microphysical diagnostics are not available
!                       in this spectral region.
     &    , dummy2d, .false., dummy2d, .false.                          &
     &    , dummy2d, .false., dummy2d, .false.                          &
     &    , dummy2d, .false.                                            &
     &    , dummy2d, .false., dummy2d, .false.                          &
     &    , dummy2d, .false., dummy2d, .false.                          &
     &    , dummy2d, .false., dummy2d, .false.                          &
     &    , col_list, row_list, row_length, rows                        &
     &    , nd_field, nd_profile, nd_layer                              &
     &    , lw_spectrum%npd_aerosol_species, 1, id_ct, id_ct            &
     &    , n_cca_lev, Ntot_land, Ntot_sea                              &
     &    )
!

          IF (lw_control%i_cloud == ip_cloud_mcica) THEN 
            ALLOCATE(subcol_k(lw_spectrum%npd_band                      & 
                      ,lw_spectrum%npd_esft_term)) 
            ALLOCATE(subcol_reorder(subcol_need)) 
            DO ll=1, lw_spectrum%npd_esft_term 
              DO l=1, lw_spectrum%n_band
                subcol_k(l,ll)=lw_subcol_k(l,ll) 
              ENDDO 
            ENDDO 
            DO l=1,subcol_need 
              subcol_reorder(l)=lw_subcol_reorder(l) 
            ENDDO 
! DEPENDS ON: mcica_order 
            CALL mcica_order(ierr                                       & 
!                       General spectral properties 
             , lw_control%first_band, lw_spectrum%n_band                & 
!                       Gaseous absorption 
             , lw_spectrum%i_band_esft, lw_spectrum%index_absorb        & 
!                       Dimensions of arrays 
             , lw_spectrum%npd_esft_term, lw_spectrum%npd_band          & 
             , lw_spectrum%npd_species                                  & 
         ) 
        
            ALLOCATE(c_sub(nd_profile, id_ct:nd_layer, tot_subcol_gen   &
              , n_cloud_type)) 
            ALLOCATE(frac_cloudy(nd_profile)) 
            DO l=1,n_points 
              i=((row_list(l)-1)*row_length)+col_list(l)
              frac_cloudy(l)=frac_cloudy_full(i)
              tot_cloud_cover(l)=REAL(ncldy(i))/REAL(tot_subcol_gen)
            ENDDO 
            IF (ALLOCATED(cic_sub_full)) THEN
              DO lll=1,subcol_need
                DO ll=id_ct,n_layer
                  DO l=1,n_points
                    i=((row_list(l)-1)*row_length)+col_list(l)
                    c_sub(l,ll,lll,1)=clw_sub_full(i,ll,lll)
                    c_sub(l,ll,lll,2)=cic_sub_full(i,ll,lll)
                  END DO
                END DO
              END DO
            ELSE
              DO lll=1,subcol_need
                DO ll=id_ct,n_layer
                  DO l=1,n_points
                    i=((row_list(l)-1)*row_length)+col_list(l)
                    c_sub(l,ll,lll,1)=clw_sub_full(i,ll,lll)
                    c_sub(l,ll,lll,2)=clw_sub_full(i,ll,lll)
                  END DO
                END DO
              END DO
            END IF
          END IF

          IF (lw_control%i_cloud/=IP_cloud_mcica) THEN
            IF (ALLOCATED(cloud_inhom_param_full)) THEN
              ALLOCATE(cloud_inhom_param(nd_profile, id_ct:nd_layer))
              DO ll=id_ct,n_layer
                DO l=1,n_points
                  i=((row_list(l)-1)*row_length)+col_list(l)
                  cloud_inhom_param(l,ll)=cloud_inhom_param_full(i,ll)
                END DO
              END DO
            END IF
          END IF
!
!
!
!       Check the consistency of cloud diagnostics
!
        if (LW_diag%L_cloud_absorptivity) then
           if (.not.LW_diag%L_cloud_weight_absorptivity) then
              cmessage =                                                &
                '*** Error: The cloud absorptivity ' //                 &
                'may be diagnosed only in conjunction ' //              &
                'with the corresponding weights.'
              ierr=i_err_fatal
              GO TO 9999
           endif
        endif
!
        if (LW_diag%L_ls_cloud_absorptivity) then
           if (.not.LW_diag%L_ls_cloud_weight_absorptivity) then
              cmessage =                                                &
                '*** Error: The layer cloud absorptivity ' //           &
                'may be diagnosed only in conjunction ' //              &
                'with the corresponding weights.'
              ierr=i_err_fatal
              GO TO 9999
           endif
        endif
!
        if (LW_diag%L_cnv_cloud_absorptivity) then
           if (.not.LW_diag%L_cnv_cloud_weight_absorptivity) then
              cmessage =                                                &
                '*** Error: The conv. cloud absorptivity ' //           &
                'may be diagnosed only in conjunction ' //              &
                'with the corresponding weights.'
              ierr=i_err_fatal
              GO TO 9999
           endif
        endif
!
!    Set flags for cloud diagnostics
        IF (LW_diag%L_ls_cloud_absorptivity  .OR. L_cosp)   &
          l_ls_cloud_absorptivity = .TRUE.
        IF (LW_diag%L_cnv_cloud_absorptivity .OR. L_cosp)   &
          l_cnv_cloud_absorptivity = .TRUE.
!
!
        IF (lw_control%l_solar_tail_flux) THEN
!         Set the incident solar flux.
          l=n_points
          WHERE (coszin(list(1:l)) > 0.0)
            solar_irrad(1:l)=scs*sc*lit(list(1:l))
            zen_0(1:l)=1.0/coszin(list(1:l))
          ELSEWHERE
            solar_irrad(1:l)=0.0
            zen_0(1:l)=1.0
          END WHERE
        END IF

!
! DEPENDS ON: radiance_calc
        call radiance_calc(ierr                                         &
!                       Logical flags for processes
     &    , lw_control%l_rayleigh, l_aerosol                            &
     &    , lw_control%l_gas, lw_control%l_continuum                    &
     &    , lw_control%l_cloud, lw_control%l_drop, lw_control%l_ice     &
!                       Angular integration
     &    , lw_control%i_angular_integration, lw_control%l_rescale      &
     &    , n_order_forward                                             &
     &    , lw_control%i_2stream, lw_control%n_order_gauss              &
     &    , lw_control%i_truncation, lw_control%ls_global_trunc         &
     &    , lw_control%ms_min, lw_control%ms_max                        &
     &    , lw_control%accuracy_adaptive, euler_factor                  &
          , lw_control%l_henyey_greenstein_pf, .FALSE.                  &
          , lw_control%ls_brdf_trunc                                    &
     &    , lw_control%i_sph_algorithm, lw_control%n_order_phase_solar  &
     &    , n_viewing_direction, viewing_direction_g                    &
     &    , n_viewing_level, viewing_level                              &
     &    , lw_control%i_sph_mode                                       &
!                       Treatment of scattering
     &    , i_scatter_method_band                                       &
!                       Options for treating clouds
     &    , lw_control%l_global_cloud_top, n_cloud_top_global           &
     &    , l_inhom_cloud, inhom_cloud                                  &
!                       Options for solver
     &    , lw_control%i_solver                                         &
!                       Properties of diagnostics
     &    , map_channel                                                 &
!                       General spectral properties
     &    , lw_spectrum%n_band                                          &
          , lw_control%first_band, lw_spectrum%n_band                   &
     &    , weight_band, lw_spectrum%l_present(14)                      &
     &    , lw_spectrum%n_band_exclude, lw_spectrum%index_exclude       &
!                       General atmospheric properties
     &    , n_points, n_layer                                           &
     &    , p, t, t_surface, t_bdy, d_mass                              &
!                       Spectral region
     &    , lw_control%isolir                                           &
!                       Solar fields
     &    , zen_0, solar_irrad, lw_spectrum%solar_flux_band             &
          , lw_spectrum%solar_flux_band_ses                             &
          , lw_control%l_solar_tail_flux, solar_tail_flux               &
     &    , lw_spectrum%rayleigh_coefficient                            &
!                       Infra-red fields
     &    , lw_spectrum%n_deg_fit, lw_spectrum%thermal_coefficient      &
     &    , lw_spectrum%t_ref_planck, lw_control%l_ir_source_quad       &
!                       Gaseous absorption
     &    , i_gas_overlap, i_gas, gas_mix_ratio                         &
     &    , lw_spectrum%n_band_absorb, lw_spectrum%index_absorb         &
     &    , lw_spectrum%i_band_esft                                     &
     &    , lw_spectrum%w_esft, lw_spectrum%k_esft                      &
     &    , lw_spectrum%i_scale_esft, lw_spectrum%i_scale_fnc           &
     &    , lw_spectrum%scale_vector                                    &
     &    , lw_spectrum%p_reference, lw_spectrum%t_reference            &
     &    , l_mod_k_flux, lw_spectrum%mix_gas_band                      &
          , lw_spectrum%n_mix_gas, lw_spectrum%index_mix_gas            &
          , lw_spectrum%f_mix, lw_spectrum%i_band_esft_ses              &
          , lw_spectrum%k_esft_ses                                      &
          , lw_spectrum%k_mix_gas, lw_spectrum%w_esft_ses               &
!                       Doppler broadening
     &    , lw_spectrum%l_doppler_present                               &
     &    , lw_spectrum%doppler_correction                              &
!                       Surface fields
     &    , n_brdf_basis_fnc, rho_alb, f_brdf                           &
!                       Tiling options for heterogeneous surfaces
     &    , l_rad_tile, n_point_tile, n_tile, list_tile, rho_alb_tile   &
     &    , frac_tile, t_tile                                           &
!                       Continuum absorption
     &    , lw_spectrum%n_band_continuum                                &
     &    , lw_spectrum%index_continuum, lw_spectrum%index_water        &
     &    , lw_spectrum%k_continuum, lw_spectrum%i_scale_fnc_cont       &
     &    , lw_spectrum%scale_continuum                                 &
     &    , lw_spectrum%p_ref_continuum, lw_spectrum%t_ref_continuum    &
          , lw_spectrum%k_continuum_ses, lw_spectrum%k_h2oc             &
!                       Properties of aerosols
     &    , lw_spectrum%n_aerosol, lw_spectrum%n_aerosol_mr             &
     &    , aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index &
     &    , l_use_arcl                                                  &
     &    , lw_spectrum%aerosol_absorption                              &
     &    , lw_spectrum%aerosol_scattering                              &
     &    , lw_spectrum%n_aerosol_phf_term                              &
     &    , lw_spectrum%aerosol_phase_fnc                               &
     &    , lw_spectrum%i_aerosol_parametrization                       &
     &    , lw_spectrum%nhumidity, lw_spectrum%humidities               &
     &    , lw_spectrum%type_aerosol                                    &
!                       Optical depth of aerosols
     &    , lw_spectrum%n_aod_wavel,     lw_spectrum%aod_wavel          &
     &    , LW_diag%L_aod_sulphate,      aod_sulphate                   &
     &    , LW_diag%L_aod_dust,          aod_dust                       &
     &    , LW_diag%L_aod_seasalt,       aod_seasalt                    &
     &    , LW_diag%L_aod_soot,          aod_soot                       &
     &    , LW_diag%L_aod_biomass,       aod_biomass                    &
     &    , LW_diag%L_aod_biogenic,      aod_biogenic                   &
     &    , LW_diag%L_aod_ocff,          aod_ocff                       &
     &    , LW_diag%L_aod_delta,         aod_delta                      &
     &    , LW_diag%L_aod_nitrate,       aod_nitrate                    &
     &    , LW_diag%L_aod_total_radn,    aod_total_radn                 &
     &    , LW_diag%L_angst_total_radn,  angst_total_radn               &
     &    , LW_diag%L_aod_prog_sulphate, aod_prog_sulphate              &
     &    , LW_diag%L_aod_prog_dust,     aod_prog_dust                  &
     &    , LW_diag%L_aod_prog_seasalt,  aod_prog_seasalt               &
     &    , LW_diag%L_aod_prog_soot,     aod_prog_soot                  &
     &    , LW_diag%L_aod_prog_biomass,  aod_prog_biomass               &
     &    , LW_diag%L_aod_prog_ocff,     aod_prog_ocff                  &
     &    , LW_diag%L_aod_prog_nitrate,  aod_prog_nitrate               &
     &    , lw_spectrum%aod_absorption, lw_spectrum%aod_scattering      &
     &    , lw_spectrum%i_aod_type                                      &
!                       Properties of UKCA aerosols 
     &    , l_ukca_radaer, ukca_radaer, ukca_mix_ratio, ukca_comp_vol   & 
     &    , ukca_dry_diam, ukca_wet_diam, ukca_modal_rho, ukca_modal_vol& 
     &    , ukca_modal_wtv, ukca_modal_nbr                              & 
!                       Optical depth of UKCA aerosols 
     &    , LW_diag%L_aod_ukca_ait_sol, aod_ukca_ait_sol                & 
     &    , LW_diag%L_aod_ukca_acc_sol, aod_ukca_acc_sol                & 
     &    , LW_diag%L_aod_ukca_cor_sol, aod_ukca_cor_sol                & 
     &    , LW_diag%L_aod_ukca_ait_ins, aod_ukca_ait_ins                & 
     &    , LW_diag%L_aod_ukca_acc_ins, aod_ukca_acc_ins                &
     &    , LW_diag%L_aod_ukca_cor_ins, aod_ukca_cor_ins                &
!                       Properties of clouds
     &    , n_condensed, type_condensed                                 &
     &    , lw_control%i_cloud, lw_control%i_cloud_representation       &
     &    , w_cloud, n_cloud_type, frac_cloud, tot_cloud_cover          &
     &    , condensed_mix_ratio, condensed_dim_char                     &
     &    , i_condensed_param, condensed_n_phf, condensed_param_list    &
     &    , dp_corr_strat, dp_corr_conv                                 &
!                       Calculated Fluxes
     &    , flux_direct, flux_net, flux_up                              &
     &    , dummy3d, dummy3d, dummy3d                                   &
     &    , l_dummy, l_dummy, l_dummy                                   &
!                       Calculated Radiances
     &    , radiance, dummy2D                                           &
!                       Options for clear-sky fluxes
     &    , l_clear, i_solver_clear                                     &
!                       Clear-sky fluxes calculated
     &    , flux_direct_clear, flux_net_clear, flux_up_clear            &
!                       Special Surface Fluxes
     &    , .false., dummy1d, dummy1d, dummy1d, dummy1d, dummy1d        &
          , .FALSE., dummy1d, .FALSE., dummy1d                          &
!                       Tiled Surface Fluxes
     &    , flux_up_tile, dummy3d                                       &
!                       Arrays for diagnostics specific to the UM
     &    , l_dummy, dummy2d, dummy2d                                   &
     &    , l_dummy, dummy2d, dummy2d                                   &
     &    , l_dummy, dummy2d, dummy2d                                   &
     &    , LW_diag%L_cloud_absorptivity, cloud_absorptivity_g          &
     &    , cloud_weight_absorptivity_g                                 &
          , L_ls_cloud_absorptivity, ls_cloud_absorptivity_g            &
     &    , ls_cloud_weight_absorptivity_g                              &
          , L_cnv_cloud_absorptivity, cnv_cloud_absorptivity_g          &
     &    , cnv_cloud_weight_absorptivity_g                             &
!                       Dimensions of arrays
     &    , nd_profile, nd_layer, nd_column, nd_layer_clr, id_ct        &
     &    , nd_2sg_profile, nd_flux_profile, nd_radiance_profile        &
     &    , nd_j_profile, nd_channel                                    &
     &    , lw_spectrum%npd_band, lw_spectrum%npd_species               &
     &    , lw_spectrum%npd_esft_term, lw_spectrum%npd_scale_variable   &
     &    , lw_spectrum%npd_continuum                                   &
     &    , lw_spectrum%npd_aerosol_species                             & 
     &    , lw_spectrum%npd_aerosol_mixratio                            & 
     &    , lw_spectrum%npd_humidities                                  &
     &    , lw_spectrum%npd_cloud_parameter                             &
     &    , lw_spectrum%npd_thermal_coeff, nd_source_coeff              &
     &    , nd_brdf_basis_fnc, nd_brdf_trunc, lw_spectrum%npd_aod_wavel &
     &    , lw_spectrum%npd_phase_term, nd_max_order, nd_sph_coeff      &
     &    , nd_direction, nd_viewing_level                              &
     &    , nd_region, nd_cloud_type, nd_cloud_component                &
     &    , nd_overlap_coeff, nd_point_tile, nd_tile                    &
          , lw_spectrum%npd_tmp, lw_spectrum%npd_pre                    &
          , lw_spectrum%npd_mix, lw_spectrum%npd_band_mix_gas           &
     &    , n_ukca_mode, n_ukca_cpnt, lw_spectrum%npd_exclude           &
     &    )
        if (ierr /= i_normal) then 
          cmessage = 'Error following call to radiance_calc'
          GO TO 9999
        end if
!
!
      endif
!
!     Assignment of diagnostics:
!
!     Processing depends on whether the code has been invoked to
!     calculate radiances or fluxes.
      if ( (lw_control%i_angular_integration == ip_two_stream).or.      &
     &     ( (lw_control%i_angular_integration ==                       &
     &        ip_spherical_harmonic).and.                               &
     &        (lw_control%i_sph_mode == ip_sph_mode_flux) ) ) then


!       Total upward & downward all-sky & clear-sky flux at all levels.
!       (Note at this point "flux_net" holds the downward flux.)
        IF (LW_diag%l_flux_up) THEN
          DO i=1, nlevs+1
            DO l=1, n_points
              LW_diag%flux_up(col_list(l), row_list(l),i)               &
                = flux_up(l,n_layer+1-i, 1)
            END DO
          END DO
        END IF
        IF (LW_diag%l_flux_down) THEN
          DO i=1, nlevs+1
            DO l=1, n_points
              LW_diag%flux_down(col_list(l), row_list(l),i)             &
                = flux_net(l,n_layer+1-i, 1)
            END DO
          END DO
        END IF
        IF (LW_diag%l_flux_up_clear) THEN
          DO i=1, nlevs+1
            DO l=1, n_points
              LW_diag%flux_up_clear(col_list(l), row_list(l),i)         &
                = flux_up_clear(l,n_layer+1-i, 1)
            END DO
          END DO
        END IF
        IF (LW_diag%l_flux_down_clear) THEN
          DO i=1, nlevs+1
            DO l=1, n_points
              LW_diag%flux_down_clear(col_list(l), row_list(l),i)       &
                = flux_net_clear(l,n_layer+1-i, 1)
            END DO
          END DO
        END IF


!       Convert downward fluxes to net fluxes.
        do i=0, n_layer
          do l=1, n_points
            flux_net(l, i, 1)=flux_net(l, i, 1)-flux_up(l, i, 1)
          enddo
        enddo
        if (l_clear) then
          do i=0, n_layer
            do l=1, n_points
              flux_net_clear(l, i, 1)                                   &
     &          =flux_net_clear(l, i, 1)-flux_up_clear(l, i, 1)
            enddo
          enddo
        endif
!
!
!       OLR:
!
        IF ( lw_control%l_solar_tail_flux ) THEN
!         Solar tail flux needs to be removed for consistency
!         with normal OLR
          do l=1, n_points
            olr(list(l))=-(flux_net(l, 0, 1)-solar_tail_flux(l)/zen_0(l))
          enddo
          if (LW_diag%L_clear_olr) then
            do l=1, n_points
              LW_diag%clear_olr(col_list(l), row_list(l))               &
                =-(flux_net_clear(l, 0, 1)-solar_tail_flux(l)/zen_0(l))
            enddo
          endif
        ELSE
          DO l=1, n_points
            olr(list(l))=-flux_net(l, 0, 1)
          END DO
          IF (LW_diag%L_clear_olr) THEN
            DO l=1, n_points
              LW_diag%clear_olr(col_list(l), row_list(l))               &
                =-flux_net_clear(l, 0, 1)
            END DO
          END IF
        END IF
!
!
!       Total cloud cover:
!
        if (LW_diag%l_total_cloud_cover) then
          if ( (lw_control%i_cloud == ip_cloud_mix_max).or.             &
               (lw_control%i_cloud == ip_cloud_mix_random).or.          &
               (lw_control%i_cloud == ip_cloud_triple).or.              &
               (lw_control%i_cloud == ip_cloud_part_corr).or.           &
               (lw_control%i_cloud == ip_cloud_part_corr_cnv).or.       &
               (lw_control%i_cloud == ip_cloud_mcica) ) then
            do l=1, n_points
              LW_diag%total_cloud_cover(col_list(l), row_list(l))       &
                =tot_cloud_cover(l)
            enddo
          else
! DEPENDS ON: r2_calc_total_cloud_cover
            call r2_calc_total_cloud_cover(n_points, nclds, nclds       &
              , lw_control%i_cloud, w_cloud, total_cloud_cover_g        &
              , nd_profile, nd_layer)
            do l=1, n_points
              LW_diag%total_cloud_cover(col_list(l), row_list(l))       &
                =total_cloud_cover_g(l)
            enddo
          endif
        endif
!
!
!       Net flux at the tropopause:
!
        if (LW_diag%l_net_flux_trop) then
          do l=1, n_points
            LW_diag%net_flux_trop(col_list(l), row_list(l))             &
     &        =flux_net(l, n_layer+1-trindx(list(l)), 1)
          enddo
        endif
!
!
!       Downward flux at the tropopause:
!
        if (LW_diag%l_down_flux_trop) then
          do l=1, n_points
            LW_diag%down_flux_trop(col_list(l), row_list(l))            &
     &        =flux_net(l, n_layer+1-trindx(list(l)), 1)                &
     &        +flux_up(l, n_layer+1-trindx(list(l)), 1)
          enddo
        endif
!
!
!       Downward flux at the surface:
!
        DO l=1, n_points
          lw_down(list(l)) = flux_net(l, n_layer, 1)                    &
                           + flux_up(l, n_layer, 1)
        END DO
!
!
!       Clear-sky downward flux at the surface:
!
        if (LW_diag%l_surf_down_clr) then
          do l=1, n_points
            LW_diag%surf_down_clr(col_list(l), row_list(l))             &
     &        =flux_net_clear(l, n_layer, 1)                            &
     &        +flux_up_clear(l, n_layer, 1)
          enddo
        endif
!
!
!       Cloud absorptivity diagnostics
!
        if (LW_diag%L_cloud_absorptivity) then
          do i=1, nclds
            do l=1, n_points
              LW_diag%cloud_absorptivity(col_list(l),row_list(l),i)     &
                 =cloud_absorptivity_g(l, n_layer+1-i)
               LW_diag%cloud_weight_absorptivity(col_list(l),           &
                 row_list(l), i)                                        &
                 =cloud_weight_absorptivity_g(l, n_layer+1-i)
            enddo
          enddo
        endif
!
        if (LW_diag%L_ls_cloud_absorptivity) then
          do i=1, nclds
            do l=1, n_points
              LW_diag%ls_cloud_absorptivity(col_list(l),row_list(l),i)  &
                 =ls_cloud_absorptivity_g(l, n_layer+1-i)
              LW_diag%ls_cloud_weight_absorptivity(col_list(l),         &
                 row_list(l), i)                                        &
                 =ls_cloud_weight_absorptivity_g(l, n_layer+1-i)
            enddo
          enddo
        endif
!
        if (LW_diag%L_cnv_cloud_absorptivity) then
          do i=1, nclds
            do l=1, n_points
              LW_diag%cnv_cloud_absorptivity(col_list(l),row_list(l),i) &
                 =cnv_cloud_absorptivity_g(l, n_layer+1-i)
              LW_diag%cnv_cloud_weight_absorptivity(col_list(l),        &
                 row_list(l), i)                                        &
                 =cnv_cloud_weight_absorptivity_g(l, n_layer+1-i)
            enddo
          enddo
        endif
!
!       Total cloud fraction on model levels
!
        if (LW_diag%L_total_cloud_on_levels) then
          do i=1, nclds
            do l=1, n_points
              LW_diag%total_cloud_on_levels(col_list(l),row_list(l),i)  &
                 =w_cloud(l, n_layer+1-i)
            enddo
          enddo
        endif        
!
! ######################################################
! CLOUD WATER MIXING RATIOS
! ######################################################
!
! LARGE-SCALE cloud water GRIDBOX MEAN mixing ratio (LIQUID) 
! ============================================================== 
  IF (LW_diag%L_ls_qcl_rad) THEN 
    IF (lw_control%i_cloud == ip_cloud_mcica) THEN 
      DO I=1, nclds 
        DO L=1, n_points 
          lll=MIN( subcol_need,                                         &
                   ncldy((row_list(l)-1)*row_length+col_list(l)) )
          IF (lll > 0) THEN 
            LW_diag%ls_qcl_rad(col_list(L),row_list(L),I) =             & 
              SUM(frac_cloud(l,n_layer+1-i,ip_cloud_type_sw)            & 
              *condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_st_water)     & 
              *c_sub(L,n_layer+1-I,1:lll,ip_cloud_type_sw))             &
              *tot_cloud_cover(l)/REAL(lll)
          END IF 
        END DO 
      END DO 
    ELSE 
      DO I=1, nclds 
        DO L=1, n_points 
          LW_diag%ls_qcl_rad(col_list(L),row_list(L),I)                 & 
             = CONDENSED_MIX_RATIO(L,N_LAYER+1-I,IP_CLCMP_ST_WATER)     & 
             * FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_SW)               & 
             * W_CLOUD(L,N_LAYER+1-I) 
        END DO 
      END DO 
    END IF 
  END IF

! 
! LARGE-SCALE cloud water GRIDBOX MEAN mixing ratio (ICE) 
! ============================================================== 
  IF (LW_diag%L_ls_qcf_rad) THEN 
    IF (lw_control%i_cloud == ip_cloud_mcica) THEN 
      DO I=1, nclds 
        DO L=1, n_points 
          lll=MIN( subcol_need,                                         &
                   ncldy((row_list(l)-1)*row_length+col_list(l)) )
          IF (lll > 0) THEN 
            LW_diag%ls_qcf_rad(col_list(L),row_list(L),I) =             & 
              SUM(frac_cloud(l,n_layer+1-i,ip_cloud_type_si)            & 
              *condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_st_ice)       & 
              *c_sub(L,n_layer+1-I,1:lll,ip_cloud_type_si))             &
              *tot_cloud_cover(l)/REAL(lll)
          END IF 
        END DO 
      END DO 
    ELSE 
      DO I=1, nclds 
        DO L=1, n_points 
          LW_diag%ls_qcf_rad(col_list(L),row_list(L),I)                 & 
             = CONDENSED_MIX_RATIO(L,N_LAYER+1-I,IP_CLCMP_ST_ICE)       & 
             * FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_SI)               & 
             * W_CLOUD(L,N_LAYER+1-I) 
        END DO 
      END DO 
    END IF 
  END IF 
 
! 
! CONVECTIVE cloud water GRIDBOX MEAN mixing ratio (LIQUID) 
! ============================================================== 
  IF (LW_diag%L_cc_qcl_rad) THEN 
    DO I=1, nclds 
      DO L=1, n_points 
        LW_diag%cc_qcl_rad(col_list(L),row_list(L),I)                   & 
           = CONDENSED_MIX_RATIO(L,N_LAYER+1-I,IP_CLCMP_CNV_WATER)      & 
           * FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_CW)                 & 
           * W_CLOUD(L,N_LAYER+1-I) 
      END DO 
    END DO 
  END IF 

! 
! CONVECTIVE cloud water GRIDBOX MEAN mixing ratio (ICE) 
! ============================================================== 
 
  IF (LW_diag%L_cc_qcf_rad) THEN 
    DO I=1, nclds 
      DO L=1, n_points 
        LW_diag%cc_qcf_rad(col_list(L),row_list(L),I)                   & 
           = CONDENSED_MIX_RATIO(L,N_LAYER+1-I,IP_CLCMP_CNV_ICE)        & 
           * FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_CI)                 & 
           * W_CLOUD(L,N_LAYER+1-I) 
      END DO 
    END DO 
  END IF 
!
! ######################################################
! CLOUD FRACTIONS
! ######################################################
!
! LARGE-SCALE cloud GRIDBOX FRACTION seen by radiation. (LIQUID) 
! ============================================================== 
  IF (LW_diag%L_ls_cl_rad) THEN 
    DO I=1, nclds 
      DO L=1, n_points 
        LW_diag%ls_cl_rad(col_list(L),row_list(L),I)                    & 
           = FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_SW)                 & 
           * W_CLOUD(L,N_LAYER+1-I) 
      END DO 
    END DO 
  END IF 
 
! 
! LARGE-SCALE cloud GRIDBOX fraction seen by radiation. (ICE) 
! ============================================================== 
  IF (LW_diag%L_ls_cf_rad) THEN 
    DO I=1, nclds 
      DO L=1, n_points 
        LW_diag%ls_cf_rad(col_list(L),row_list(L),I)                    & 
           = FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_SI)                 & 
           * W_CLOUD(L,N_LAYER+1-I) 
      END DO 
    END DO 
  END IF 
 
! 
! CONVECTIVE cloud GRIDBOX fraction seen by radiation. (LIQUID) 
! ============================================================== 
  IF (LW_diag%L_cc_cl_rad) THEN 
    DO I=1, nclds 
      DO L=1, n_points 
        LW_diag%cc_cl_rad(col_list(L),row_list(L),I)                    & 
           = FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_CW)                 & 
           * W_CLOUD(L,N_LAYER+1-I) 
      END DO 
    END DO 
  END IF 
 
! 
! CONVECTIVE cloud GRIDBOX FRACTION seen by radiation. (ICE) 
! ============================================================== 
  IF (LW_diag%L_cc_cf_rad) THEN 
    DO I=1, nclds 
      DO L=1, n_points 
        LW_diag%cc_cf_rad(col_list(L),row_list(L), I)                   & 
           = FRAC_CLOUD(L,N_LAYER+1-I,IP_CLOUD_TYPE_CI)                 & 
           * W_CLOUD(L,N_LAYER+1-I) 
      END DO 
    END DO 
  END IF  

! ######################################################
! COSP arguments
! ######################################################
!
  IF (l_cosp) THEN
    DO i=1, nclds
      DO l=1, n_points
        j = (row_list(l)-1)*row_length + col_list(l)
!       LARGE-SCALE cloud water GRIDBOX MEAN mixing ratio (LIQUID)
        cosp_gbx%mr_hydro(j,i,i_lscliq) =                                      &
           condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_st_water)                &
           * frac_cloud(l,n_layer+1-i,ip_cloud_type_sw) * w_cloud(l,n_layer+1-i)
!       LARGE-SCALE cloud water GRIDBOX MEAN mixing ratio (ICE)
        cosp_gbx%mr_hydro(j,i,i_lscice) =                                      &
           condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_st_ice)                  &
           * frac_cloud(l,n_layer+1-i,ip_cloud_type_si) * w_cloud(l,n_layer+1-i)
!       CONVECTIVE cloud water GRIDBOX MEAN mixing ratio (LIQUID)
        cosp_gbx%mr_hydro(j,i,i_cvcliq) =                                      &
           condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_cnv_water)               &
           * frac_cloud(l,n_layer+1-i,ip_cloud_type_cw) * w_cloud(l,n_layer+1-i)
!       CONVECTIVE cloud water GRIDBOX MEAN mixing ratio (ICE)
        cosp_gbx%mr_hydro(j,i,i_cvcice) =                                      &
           condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_cnv_ice)                 &
           * frac_cloud(l,n_layer+1-i,ip_cloud_type_ci) * w_cloud(l,n_layer+1-i)
!       TOTAL cloud GRIDBOX fraction seen by radiation
        cosp_gbx%tca(j,i) = w_cloud(l,n_layer+1-i)
!       CONVECTIVE cloud GRIDBOX fraction seen by radiation
        cosp_gbx%cca(j,i) =                                                    &
           (frac_cloud(l,n_layer+1-i,ip_cloud_type_cw) +                       &
           frac_cloud(l,n_layer+1-i,ip_cloud_type_ci)) * w_cloud(l,n_layer+1-i)
!       LARGE-SCALE cloud water effective radius
        cosp_gbx%reff(j,i,i_lscliq) =                                          &
           condensed_dim_char(l,n_layer+1-i,ip_clcmp_st_water)
        cosp_gbx%reff(j,i,i_lscice) =                                          &
           condensed_dim_char(l,n_layer+1-i,ip_clcmp_st_ice)
!       CONVECTIVE cloud water effective radius
        cosp_gbx%reff(j,i,i_cvcliq) =                                          &
           condensed_dim_char(l,n_layer+1-i,ip_clcmp_cnv_water)
        cosp_gbx%reff(j,i,i_cvcice) =                                          &
           condensed_dim_char(l,n_layer+1-i,ip_clcmp_cnv_ice)
!       Large-Scale and convective cloud emissivities
        IF (ls_cloud_weight_absorptivity_g(l, n_layer+1-i) > 1.0E-20) THEN
            cosp_gbx%dem_s(j,i) = 1.0-EXP(-1.666*d_mass(l,n_layer+1-i)         &
                * ls_cloud_absorptivity_g(l, n_layer+1-i)                      &
                / ls_cloud_weight_absorptivity_g(l, n_layer+1-i))
        END IF
        IF (cnv_cloud_weight_absorptivity_g(l, n_layer+1-i) > 1.0E-20) THEN
            cosp_gbx%dem_c(j,i) = 1.0-EXP(-1.666*d_mass(l,n_layer+1-i)         &
                * cnv_cloud_absorptivity_g(l, n_layer+1-i)                     &
                / cnv_cloud_weight_absorptivity_g(l, n_layer+1-i))
        END IF
      END DO
    END DO
  END IF

!
!     Aerosol Optical Depth
!
      if (LW_diag%l_aod_sulphate) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_sulphate(col_list(l), row_list(l), i)           &
             = aod_sulphate(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_dust) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_dust(col_list(l), row_list(l), i)               &
             = aod_dust(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_seasalt) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_seasalt(col_list(l), row_list(l), i)            &
             = aod_seasalt(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_soot) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_soot(col_list(l), row_list(l), i)               &
             = aod_soot(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_biomass) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_biomass(col_list(l), row_list(l), i)            &
             = aod_biomass(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_biogenic) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_biogenic(col_list(l), row_list(l), i)           &
             = aod_biogenic(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_ocff) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_ocff(col_list(l), row_list(l), i)               &
             = aod_ocff(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_delta) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_delta(col_list(l), row_list(l), i)              &
             = aod_delta(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_nitrate) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_nitrate(col_list(l), row_list(l), i)            &
             = aod_nitrate(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_total_radn) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_total_radn(col_list(l), row_list(l), i)         &
             = aod_total_radn(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_angst_total_radn) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%angst_total_radn(col_list(l), row_list(l), i)         &
             = angst_total_radn(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_prog_sulphate) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_prog_sulphate(col_list(l), row_list(l), i)      &
             = aod_prog_sulphate(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_prog_dust) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_prog_dust(col_list(l), row_list(l), i)          &
             = aod_prog_dust(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_prog_seasalt) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_prog_seasalt(col_list(l), row_list(l), i)          &
             = aod_prog_seasalt(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_prog_soot) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_prog_soot(col_list(l), row_list(l), i)          &
             = aod_prog_soot(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_prog_biomass) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_prog_biomass(col_list(l), row_list(l), i)       &
             = aod_prog_biomass(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_prog_ocff) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_prog_ocff(col_list(l), row_list(l), i)          &
             = aod_prog_ocff(l, i)
          enddo
        enddo
      endif
      if (LW_diag%l_aod_prog_nitrate) then
        do i=1, lw_spectrum%n_aod_wavel
          do l=1, n_points
            LW_diag%aod_prog_nitrate(col_list(l), row_list(l), i)       &
             = aod_prog_nitrate(l, i)
          enddo
        enddo
      endif

! 
!    UKCA aerosol optical depth 
! 
      IF (LW_diag%L_aod_ukca_ait_sol) THEN 
        DO i=1, lw_spectrum%n_aod_wavel
          DO l=1, n_points
            LW_diag%aod_ukca_ait_sol(col_list(l), row_list(l), i)       & 
             = aod_ukca_ait_sol(l, i) 
          END DO 
        END DO 
      END IF 
      IF (LW_diag%L_aod_ukca_acc_sol) THEN 
        DO i=1, lw_spectrum%n_aod_wavel
          DO l=1, n_points
            LW_diag%aod_ukca_acc_sol(col_list(l), row_list(l), i)       & 
             = aod_ukca_acc_sol(l, i) 
          END DO 
        END DO 
      END IF 
      IF (LW_diag%L_aod_ukca_cor_sol) THEN 
        DO i=1, lw_spectrum%n_aod_wavel
          DO l=1, n_points
            LW_diag%aod_ukca_cor_sol(col_list(l), row_list(l), i)       & 
             = aod_ukca_cor_sol(l, i) 
          END DO 
        END DO 
      END IF 
      IF (LW_diag%L_aod_ukca_ait_ins) THEN 
        DO i=1, lw_spectrum%n_aod_wavel
          DO l=1, n_points
            LW_diag%aod_ukca_ait_ins(col_list(l), row_list(l), i)       & 
             = aod_ukca_ait_ins(l, i) 
          END DO 
        END DO 
      END IF 
      IF (LW_diag%L_aod_ukca_acc_ins) THEN 
        DO i=1, lw_spectrum%n_aod_wavel
          DO l=1, n_points
            LW_diag%aod_ukca_acc_ins(col_list(l), row_list(l), i)       & 
             = aod_ukca_acc_ins(l, i) 
          END DO 
        END DO 
      END IF 
      IF (LW_diag%L_aod_ukca_cor_ins) THEN 
        DO i=1, lw_spectrum%n_aod_wavel
          DO l=1, n_points
            LW_diag%aod_ukca_cor_ins(col_list(l), row_list(l), i)       & 
             = aod_ukca_cor_ins(l, i) 
         END DO 
        END DO
      END IF
!
!
!
!       Output arrays:
!
!       Convert the fluxes to increments in the heating rate except at
!       the surface: there, the net downward flux is assigned to LWOUT.
        do k=nlevs, 1, -1
          if (l_scale_inc) then
            do l=1, n_points
              lwout(list(l), k+1)=(flux_net(l, n_layer-k, 1)            &
     &          -flux_net(l, n_layer+1-k, 1))                           &
     &          *pts/layer_heat_capacity(l, n_layer+1-k)
            enddo
          else
            do l=1, n_points
              lwout(list(l), k+1)=(flux_net(l, n_layer-k, 1)            &
     &          -flux_net(l, n_layer+1-k, 1))                           &
     &          /layer_heat_capacity(l, n_layer+1-k)
            enddo
          endif
          if (LW_diag%l_clear_hr) then
!           The factor of PTS is included here to yield a rate from an
!           increment.
            do l=1, n_points
              LW_diag%clear_hr(col_list(l), row_list(l), k)             &
     &          =(flux_net_clear(l, n_layer-k, 1)                       &
     &          -flux_net_clear(l, n_layer+1-k, 1))                     &
     &          /layer_heat_capacity(l, n_layer+1-k)
            enddo
          endif
        enddo
!
        if (lw_control%l_extra_top) then
!         calculate the radiation absorbed in the extra layer
!         above the top of the rest of the model.
          do l=1, n_points
            top_absorption(list(l))=flux_net(l, 0, 1)                   &
     &        -flux_net(l, n_layer-nlevs, 1)
          enddo
        endif
!
        do l=1, n_points
          lwout(list(l), 1)=flux_net(l, n_layer, 1)
        enddo
!
!       Separate the contributions over open sea and sea-ice.
!       Fluxes returned from the radiation code itself are not
!       weighted by the fraction of the tile, but here are converted
!       to grid-box mean values. This split is possible only if the
!       ocean surface has been tiled.
!
        if (l_rad_tile) then
!
!         The variable flandg is set even if coastal tiling is not
!         used, so fairly generic code can be written.
!
!         It is simplest to zero LWsea at all points and reset
!         elsewhere.
          lwsea(list(1:n_points)) = 0.0
!
          do ll=1, n_points
            l=list(ll)
            if ( (flandg(l) < TINY(flandg)) .AND.                       &
     &           (ice_fraction(l) < TINY(ice_fraction) )                &
     &         ) then
!             This point is open sea with no sea ice.
              lwsea(l)=lwout(l, 1)
              lwout(l, 1)=0.0
            endif
          enddo
!
!         Tiled points will have both land and sea. Note that the
!         channel index of flux_up_tile is hard-wired to 1 because
!         we don't envisage calling the code in other cases.
          do lll=1, n_point_tile
            ll=list_tile(lll)
            l=list_tile_outer(lll)
            lwsea(l)=(1.0-ice_fraction(l))*(flux_net(ll, n_layer, 1)    &
     &        +flux_up(ll, n_layer, 1)                                  &
     &        -flux_up_tile(lll, index_tile(ip_ocean_tile), 1))
            lwout(l, 1)=lwout(l, 1)-(1.0-flandg(l))*lwsea(l)
          enddo
!
!         The remaining points are entirely land points and lwout
!         need not be altered.
!
        else
!
!         Without radiative tiling we must assume that fluxes are
!         uniform across the grid-box.
          Where (flandg(list(1:n_points)) < 1.0-TINY(flandg))
            lwsea(list(1:n_points))=(1.0-ice_fraction(list(1:n_points)))&
     &        *flux_net(1:n_points, n_layer, 1)
            lwout(list(1:n_points), 1)=lwout(list(1:n_points), 1)       &
     &        -(1.0-flandg(list(1:n_points)))*lwsea(list(1:n_points))
          Endwhere
!
        endif
!
      else
!
!       Radiances are being calculated.
!
        if (lw_diag%l_toa_radiance) then
          do ic=1, n_channel
            do l=1, n_points
              lw_diag%toa_radiance(col_list(l), row_list(l), ic)        &
     &          =radiance(l, 1, 1, ic)
            enddo
          enddo
        endif
!
      endif
!
!
      IF (ALLOCATED(cloud_inhom_param)) DEALLOCATE(cloud_inhom_param)
      IF (lw_control%i_cloud == ip_cloud_mcica) THEN 
        DEALLOCATE(first_subcol_k) 
        DEALLOCATE(first_subcol_band) 
        DEALLOCATE(subcol_band) 
        DEALLOCATE(frac_cloudy) 
        DEALLOCATE(c_sub) 
        DEALLOCATE(subcol_reorder) 
        DEALLOCATE(subcol_k) 
      END IF 

 9999 CONTINUE
! Check error condition
      IF (ierr /= i_normal) THEN

        CALL ereport(RoutineName, ierr, cmessage)
      END IF

      IF (lhook) CALL dr_hook('R2_LWRAD3Z',zhook_out,zhook_handle)
      END SUBROUTINE r2_lwrad3z
