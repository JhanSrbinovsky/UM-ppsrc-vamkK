! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Shortwave interface to the radiance code.
!
! Purpose:
!   This subroutine interface the radiance version of the ES code
!   in the shortwave.
!
! Method:
!   Principally, arrays are transferred to the appropriate formats.
!   separate subroutines are called for each physical process.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------
      SUBROUTINE r2_swrad3z(ierr                                        &
!                       Mixing ratios
        , h2o, co2, o3, o2_mix_ratio                                    &
        , co2_dim1, co2_dim2, co2_3d, l_co2_3d                          &
        , n2o_mix_ratio, ch4_mix_ratio                                  &
!                       Pressure fields
        , pstar                                                         &
        , p_layer_boundaries                                            &
        , p_layer_centres                                               &
!                       Temperatures
        , tac                                                           &
!                       Options for COSP
        , l_cosp                                                        &
!                       Options for treating clouds
        , global_cloud_top                                              &
        , l_inhom_cloud, inhom_cloud, dp_corr_strat, dp_corr_conv       &
!                       Stratiform cloud fields
        , l_pc2, lca_area, lca_bulk, lccwc1, lccwc2, n_drop_pot         &
!                       Convective cloud fields
        , cca, cccwp, ccw, lcbase, ccb, cct                             &
!                       Surface fields
        , land_albedo, l_ctile, l_use_spec_sea, flandg, sea_ice_albedo  &
        , open_sea_albedo, ice_fraction, land, land0p5, lying_snow      &
!                       Solar fields
        , coszin, lit, list, scs                                        &
!                       Aerosol fields
        , l_climat_aerosol, l_clim_aero_hgt, L_HadGEM1_Clim_Aero        &
        , bl_depth, n_levels_bl                                         &
        , l_dust, l_use_dust, dust_dim1, dust_dim2                      &
        , dust_1, dust_2, dust_3, dust_4, dust_5, dust_6                &
        , l_use_biogenic, biogenic_dim1, biogenic_dim2, biogenic        &
        , l_sulpc_so2, l_use_sulpc_direct, l_use_sulpc_indirect         &
        , sulp_dim1, sulp_dim2                                          &
        , accum_sulphate, aitken_sulphate, diss_sulphate                &
        , sea_salt_film, sea_salt_jet, l_use_seasalt_indirect           &
        , l_use_seasalt_direct, salt_dim_a, salt_dim_b                  &
        , l_soot, l_use_soot_direct, soot_dim1, soot_dim2               &
        , fresh_soot, aged_soot                                         &
        , l_biomass, l_use_bmass_direct, bmass_dim1, bmass_dim2         &
        , fresh_bmass, aged_bmass, cloud_bmass, l_use_bmass_indirect    &
        , l_ocff, l_use_ocff_direct, ocff_dim1, ocff_dim2               &
        , fresh_ocff, aged_ocff, cloud_ocff, l_use_ocff_indirect        &
        , l_nitrate, l_use_nitrate_direct, nitrate_dim1, nitrate_dim2   &
        , accum_nitrate, diss_nitrate, l_use_nitrate_indirect           &
        , l_use_arcl, arcl_dim1, arcl_dim2, n_arcl_species              &
        , n_arcl_compnts, i_arcl_compnts, arcl                          &
        , aero_meso, l_murk_rad, Ntot_land, Ntot_sea                    &
        , l_ukca_radaer, ukca_radaer, ukca_dim1, ukca_dim2              &
        , ukca_mmr, ukca_cvl, ukca_dry, ukca_wet                        &
        , ukca_rho, ukca_vol, ukca_wtv, ukca_nbr                        & 
!                       Time
        , PREVIOUS_TIME                                                 &
!                       Grid-dependent arrays
        , true_latitude                                                 &
!                       Level of tropopause
        , trindx                                                        &
!                       Spectrum
        , sw_spectrum                                                   &
!                       Algorithmic options
        , sw_control                                                    &
        , pts, l_mod_k_flux, l_scale_inc                                &
!                       Satellite viewing geometry
        , n_viewing_direction,viewing_direction1,viewing_direction2     &
        , n_viewing_level, viewing_level                                &
!                       diagnostics
        , n_channel, map_channel                                        &
        , sw_diag, row_list, col_list                                   &
!                       Physical dimensions
        , nlit, n_points, nlevs, n_layer, nclds                         &
        , nwet, nozone, row_length, rows                                &
        , nd_field, nd_field_flux_diag, nd_field_rad_diag               &
        , nd_profile, nd_layer, nd_column                               &
        , n_cca_lev, nd_channel, nd_flux_profile                        &
        , nd_radiance_profile, nd_viewing_level, nd_direction           &
        , nd_cloud_component, nd_cloud_type                             &
        , nd_brdf_basis_fnc, nd_brdf_trunc                              &
        , nd_point_tile, nd_tile, id_ct, n_ukca_mode, n_ukca_cpnt       &
!                       Output
        , surf_down_sw                                                  &
        , flux_below_690nm_surf, l_flux_below_690nm_surf                &
        , netsw, top_absorption, swsea, swout                           &
          ! Variables needed to calculate layer masses
        , rho_r2, r_rho_levels, r_theta_levels                          &
        , q, qcl, qcf, qcf2, qrain, qgraup                              &
        , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio         &
! COSP input arguments
        , cosp_gbx                                                      &
           )
!
!
!
!     Modules included.
      USE rad_pcf
      USE conversions_mod, ONLY: recip_pi_over_180
      use dec_spec
      use control_struc
      use tileid3z
      USE solinc_data, ONLY: sol_bearing, lg_orog_corr, orog_corr,      &
     &                       lg_f_orog, f_orog, L_orog
      USE swrdiag_mod, ONLY: strswdiag
      USE mcica_mod
      USE ukca_radaer_struct_mod
      USE yomhook,     ONLY: lhook, dr_hook
      USE parkind1,    ONLY: jprb, jpim
      USE cosp_types_mod, ONLY: cosp_gridbox
      USE ereport_mod, ONLY : ereport
      USE arcl_mod,    ONLY: npd_arcl_species, npd_arcl_compnts
      USE rad_input_mod, ONLY: sc
      USE missing_data_mod, ONLY: rmdi
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
!           Size of array of profiles
     &  , nd_layer                                                      &
!           Array sizes for layers
     &  , nd_column
!           Number of columns per point
      integer                                                           &
                !, intent(in)
     &    nd_channel                                                    &
!           Size allocated for diagnostic spectral bands
     &  , nd_flux_profile                                               &
!           Size allocated for output flux profiles
     &  , nd_radiance_profile                                           &
!           Number of profiles where radiances are required
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
     &  , nd_tile                                                       &
!           Size allocated for surface tiles
     &  , id_ct
!           Top level in arrays of cloud properties
!
!     Actual sizes used:
      integer                                                           &
                !, intent(in)
     &    n_points                                                      &
!           Number of points to be diagnosed including unlit points
     &  , nwet                                                          &
!           Number of wet levels
     &  , nozone                                                        &
!           Number of levels with ozone
     &  , nlevs                                                         &
!           number of layers in the main model
     &  , n_layer                                                       &
!           number of layers seen in the radiation scheme
     &  , nclds                                                         &
!           Number of cloudy levels
     &  , n_levels_bl                                                   &
!           Number of layers occupied by boundary-layer aerosol
!           if l_clim_aero_hgt is false.
     &  , n_cca_lev                                                     &
!           Number of convective cloud levels
     &  , n_ukca_mode                                                   &
!           Number of aerosol modes in UKCA_RADAER
     &  , n_ukca_cpnt
!           Number of aerosol components in UKCA_RADAER
!
!     Spectral data:
      type (spectrum) :: sw_spectrum
!
!     Controlling options:
      type (control_option) :: sw_control
!
      logical                                                           &
     &    l_scale_inc
!           Flag for scaling of heating rates to increments
!
!     Gaseous mixing ratios
      real                                                              &
                !, intent(in)
     &    h2o(nd_field, nwet)                                           &
!           Mass mixing ratio of water
     &  , co2                                                           &
!           Mass mixing ratio of co2
     &  , o3(nd_field, nozone)                                          &
!           Mass mixing ratios of ozone
     &  , o2_mix_ratio                                                  &
!           Mass mixing ratio of oxygen
     &  , n2o_mix_ratio                                                 &
!           Mass mixing ratio of nitrous oxide
     &  , ch4_mix_ratio
!           Mass mixing ratio of methane
!
!     General atmospheric properties:
      real                                                              &
                !, intent(in)
     &    pstar(nd_field)                                               &
!           Surface pressures
     &  , p_layer_boundaries(nd_field,0:nlevs)                          &
!            pressure at boundaries of layers
     &  , p_layer_centres(nd_field,0:nlevs)                             &
!            pressure at centres of layers
     &  , tac(nd_field, nlevs)
!           Temperatures at centres of layers
!
!     Incident solar radiation:
      integer                                                           &
                !, intent(in)
     &    nlit                                                          &
!           Number of lit points
     &  , list(nd_field)
!           List of lit points
      real                                                              &
                !, intent(in)
     &    coszin(nd_field)                                              &
!           Cosines of zenith angle
     &  , scs                                                           &
!           Scaling of solar incident field
     &  , lit(nd_field)
!           Fraction of time point is lit
!
      LOGICAL                                                           &
                !, intent(in)
           l_mcr_qcf2                                                   &
                          ! Use second ice category
     &,    l_mcr_qrain                                                  &
                          ! Use prognostic rain
     &,    l_mcr_qgraup                                                 &
                          ! Use graupel
     &,    l_mixing_ratio ! Use mixing ratios in layer mass calculation

!     Flag for COSP
      LOGICAL,INTENT(IN) :: l_cosp
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
!           Nominal liquid water contents
     &  , lccwc2(nd_field, nclds+1/(nclds+1))                           &
!           Nominal ice water contents
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
!           Fraction of convective cloud

!     Aerosols:
      logical                                                           &
                !, intent(in)
     &    l_climat_aerosol                                              &
!           Flag for climatological aerosol
     &  , l_clim_aero_hgt                                               &
!           flag to use the depth of the boundary layer to set
!           the climatological aerosol
     &  , L_HadGEM1_Clim_Aero                                           &
!           Flag to use HadGEM1 setting for climatological aerosols
     &  , l_murk_rad
!           flag for mesoscale model aerosol
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
                !, intent(in)
     &    sulp_dim1,sulp_dim2                                           &
!           Dimensions for _sulphate arrays, (P_FIELD,P_LEVELS or 1,1)
     &  , dust_dim1, dust_dim2                                          &
!           Dimensions for mineral dust arrays (p_field,p_levels or 1,1)
     &  , biogenic_dim1, biogenic_dim2                                  &
!           dimensions for biogenic array passed down to
!           r2_set_aerosol_field if direct effect required.
     &  , soot_dim1, soot_dim2                                          &
!           Dimensions for soot arrays (P_FIELD,P_LEVELS or 1,1)
     &  , bmass_dim1, bmass_dim2                                        &
!           Dimensions for biomass arrays (P_FIELD,P_LEVELS or 1,1)
     &  , ocff_dim1, ocff_dim2                                          &
!           Dimensions for ocff arrays (P_FIELD,P_LEVELS or 1,1)
     &  , nitrate_dim1, nitrate_dim2                                    &
!           Dimensions for nitrate arrays (P_FIELD,P_LEVELS or 1,1)
     &  , salt_dim_a, salt_dim_b
!           dimensions for salt arrays on input (salt_dim_a=p_field
!           and salt_dim_b=p_levels, or else 1,1)
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
!             number concentration of film-mode sea-salt aerosol
     &  , sea_salt_jet(salt_dim_a, salt_dim_b)                          &
!             number concentration of jet-mode sea-salt aerosol
     &  , fresh_soot(soot_dim1, soot_dim2)                              &
!           Soot mixing ratios
     &  , aged_soot(soot_dim1, soot_dim2)                               &
!           Soot mixing ratios
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
     &  , bl_depth(nd_field)                                            &
!           depth of the boundary layer
     &  , aero_meso(nd_field, nlevs)
!           mixing ratio of 'urban' aerosol of mesoscale model
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
      logical                                                           &
     &    l_co2_3d
!           Controls use of 3D CO2 field
      integer                                                           &
                !, intent(in)
     &    co2_dim1, co2_dim2
!           Dimensions for CO2 array, (P_FIELD,P_LEVELS or 1,1)
      real                                                              &
                !, intent(in)
     &    co2_3d(co2_dim1, co2_dim2)
!           Mass mixing ratio of carbon dioxide

!     Properties of the surface:
      Logical, Intent(IN) :: Land(nd_field)
!                             Land mask
      Logical, Intent(IN) :: Land0p5(nd_field)
!                             Land mask (TRUE if land fraction > 0.5)
      Logical, Intent(IN) :: L_flux_below_690nm_surf
!                             Flag to calculate flux at wavelengths
!                             shorter than 690 nm at the surface:
!                             This may be required as a diagnostic
!                             or for MOSES
      Logical, Intent(IN) :: L_Ctile
!                             Switch for coastal tiling
      Logical, Intent(IN) :: L_use_spec_sea
!                             Switch for spectrally dep. sea albedos
      real                                                              &
                !, intent(in)
     &    ice_fraction(nd_field)                                        &
!             fraction of sea ice in sea portion of grid box
     &  , land_albedo(nd_field,4)                                       &
!             land surface albedo fields
!             (*,1) - direct beam visible
!             (*,2) - diffuse visible
!             (*,3) - direct beam near-ir
!             (*,4) - diffuse near-ir
     &  , flandg(nd_field)                                              &
!             land fraction in grid box
     &  , sea_ice_albedo(nd_field,4)                                    &
!             sea ice albedo fields
!             (*,1) - direct beam visible
!             (*,2) - diffuse visible
!             (*,3) - direct beam near-ir
!             (*,4) - diffuse near-ir
     &  , open_sea_albedo(nd_field, 2)                                  &
!           Surface albedo field of open sea
!             (direct and diffuse components)
     &  , lying_snow(nd_field)
!           Mass loading of lying snow
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
!
!                       level of tropopause
      integer                                                           &
     &    trindx(nd_field)
!           The layer boundary of the tropopause
!
!     Increment of time:
      real                                                              &
                !, intent(in)
     &    pts
!           Time increment
!
!     Use modulus of fluxes to remove negative effective extinctions
      logical, intent(in) :: l_mod_k_flux

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
!
!     Satellite viewing geometry
      INTEGER, Intent(IN) :: n_channel
!           Number of channels calculated simultaneously in one
!           call to the radiation code
      INTEGER, Intent(IN) :: map_channel(sw_spectrum%npd_band)
!           Mapping of bands in the spectral file to channels in the
!           diagnostic output
      INTEGER, Intent(IN) :: n_viewing_direction
!           Number of viewing directions
      INTEGER, Intent(IN) :: n_viewing_level
!           Number of levels where the radiance is calculated
      REAL, Intent(IN) :: viewing_direction1(nd_field, nd_direction)
!           Satellite viewing directions (Zenith Angle)
      REAL, Intent(IN) :: viewing_direction2(nd_field, nd_direction)
!           Satellite viewing directions (Azimuthal Angle)
      REAL, Intent(IN) :: viewing_level(nd_viewing_level)
!           Levels where radiances are calculated
!
!
!     Calculated fluxes:
      real                                                              &
                !, intent(out)
     &    swout(nd_field, nlevs+2)                                      &
!           Net downward fluxes
     &  , swsea(nd_field)                                               &
!           Sea-surface components of flux
!           weighted by (open sea)/(total sea) fraction
     &  , netsw(nd_field)                                               &
!           Net absorbed shortwave radiation
     &  , flux_below_690nm_surf(nd_field_flux_diag)                     &
!           Net surface flux below 690nm (at points where there
!           is sea-ice this is weighted by the fraction of open sea.)
     &  , surf_down_sw(nd_field_flux_diag, 4)                           &
!             surface downward shortwave radiation components
!             (*,1) - direct beam visible
!             (*,2) - diffuse visible
!             (*,3) - direct beam near-ir
!             (*,4) - diffuse near-ir
     &  , top_absorption(nd_field)
!             radiative absorption above the top of the atmosphere
!             as seen in the main model
!
!
!     Diagnostics:
!
!     Definition of the diagnostic structure
      type (strswdiag) :: sw_diag
!
      integer, intent(in) :: row_list(nd_field)
!                              list of row indices of lit points
      integer, intent(in) :: col_list(nd_field)
!                              list of column indices of lit points

! Structure with COSP inputs
      TYPE(cosp_gridbox),INTENT(OUT) :: cosp_gbx
!
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
     &  , nd_sph_coeff
!           Size allocated for coefficients of spherical harmonics
!
      integer                                                           &
     &    salt_dim_ind_a, salt_dim_ind_b                                &
!           dimensions for sea-salt arrays passed down to
!           r2_set_cloud_field if indirect effect required.
     &  , salt_dim_dir_a, salt_dim_dir_b
!           dimensions for sea-salt arrays passed down to
!           r2_set_aerosol_field if direct effect required.
!
      integer                                                           &
     &    i                                                             &
!           Loop variable
     &  , j                                                             &
!           Loop variable
     &  , l                                                             &
!           Loop variable
     &  , ll                                                            &
!           Loop variable
     &  , lll                                                           &
!           Loop variable
        , k                                                             &
!           Loop variable
     &  , ic
!           Loop variable
      logical                                                           &
     &    l_clear
!           Calculate clear-sky fields
!     Flags for processes actually enabled.
      logical                                                           &
     &    l_aerosol                                                     &
!           Flag to enable direct aerosol effects within
!           the radiation code
     &  , l_aerosol_ccn
!           Local flag to use aerosols to determine ccn
      integer                                                           &
     &    i_solver_clear                                                &
!           Solver for clear-sky fluxes
     &  , i_gas_overlap(sw_spectrum%npd_band)                           &
!           Overlaps in each band
     &  , i_scatter_method_band(sw_spectrum%npd_band)                   &
!           treatment of scattering in each band
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
     &  , gas_mix_ratio(nd_profile, nd_layer, sw_spectrum%npd_species)  &
!           Mass fractions of gases
     &  , nullmmr
!           Null mass mixing ratio
      parameter(                                                        &
     &     nullmmr=0.0                                                  &
     &  )
!
      real :: layer_heat_capacity(nd_profile, nd_layer)
!           Specific heat capacity of layer * d_mass

      REAL :: alat(nd_profile)
!           Latitude in degrees

!     Cloudy properties:
      integer                                                           &
     &    n_condensed                                                   &
!           Number of condensed phases
     &  , type_condensed(nd_cloud_component)                            &
!           Types of condensed components
     &  , i_condensed_param(nd_cloud_component)                         &
!           Parametrization schemes for components
     &  , condensed_n_phf(nd_cloud_component)                           &
!           Number of terms in the phase function
     &  , n_cloud_top_global                                            &
!           Inverted global topmost cloudy layer
     &  , n_cloud_type                                                  &
!           Number of types of clouds
     &  , i_cloud_tmp
!           Cloud Overlap used by r2_cloud_level_diag
      real                                                              &
     &    condensed_param_list(sw_spectrum%npd_cloud_parameter          &
     &      , nd_cloud_component, sw_spectrum%npd_band)                 &
!           Parameters for condensed phases
     &  , condensed_dim_char(nd_profile, id_ct: nd_layer                &
     &      , nd_cloud_component)                                       &
!           Characteristic dimensions of condensed species
     &  , condensed_mix_ratio(nd_profile, id_ct: nd_layer               &
     &      , nd_cloud_component)                                       &
!           Mass fractions of condensed species
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
      real                                                              &
     &    aerosol_mix_ratio(nd_profile, nd_layer                        &
     &      , sw_spectrum%npd_aerosol_mixratio)
!           Mixing ratios of aerosols
      integer                                                           &
          aerosol_mr_type_index(sw_spectrum%npd_aerosol_mixratio)       &
!             Index relating aerosol_mix_ratio aerosols to aerosols in
!             the spectral information
        , aerosol_mr_source(sw_spectrum%npd_aerosol_mixratio)
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
!     Surface properties:
      integer                                                           &
     &    n_brdf_basis_fnc
!           Number of BRDF basis functions
      LOGICAL :: land0p5_g(nd_profile)
!           Gathered land mask (true if land fraction >0.5)
      REAL :: flandg_g(nd_profile)
!           Gathered land fraction
      real                                                              &
     &    rho_alb(nd_profile, nd_brdf_basis_fnc, sw_spectrum%npd_band)  &
!           Weights for BRDF basis functions
     &  , f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                  &
     &      , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)
!           Array of BRDF basis terms
!     Variables related to tiling of the surface
      LOGICAL :: l_rad_tile
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
!           List of points with surface tiling
     &  , list_tile_outer(nd_point_tile)                                &
!           List of points with surface tiling indexed over the full
!           list of points where this routine has been called
     &  , index_tile(npd_tile_type)
!           The indexing number of tiles of the given type
      real                                                              &
     &    rho_alb_tile(nd_point_tile, nd_brdf_basis_fnc                 &
     &      , nd_tile, sw_spectrum%npd_band)                            &
!           Weights for the basis functions of the BRDFs
!           at the tiled points
     &  , t_tile(nd_point_tile, nd_tile)
!           Local surface temperatures on individual tiles
      REAL :: frac_solid
!       Total solid fraction in a grid-box
!
!     Local algorithmic variables:
      real                                                              &
     &    euler_factor
!           Weighting applied to the last term of alternating series
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
!           Net/downward flux
     &  , flux_net_clear(nd_flux_profile, 0: nd_layer, nd_channel)      &
!           Clear-sky net/downward total flux
     &  , flux_up(nd_flux_profile, 0: nd_layer, nd_channel)             &
!           Upward flux
     &  , flux_up_clear(nd_flux_profile, 0: nd_layer, nd_channel)
!           Clear-sky upward flux
!
!     UV-Fluxes
!
      real                                                              &
     &    uv_flux_direct(nd_flux_profile, 0: nd_layer, nd_channel)      &
!           direct UV-flux
     &  , uv_flux_up(nd_flux_profile, 0: nd_layer, nd_channel)          &
!           upward UV-flux
     &  , uv_flux_net(nd_flux_profile, 0: nd_layer, nd_channel)
!           net UV-flux
!     Radiances:
      real                                                              &
     &    radiance(nd_radiance_profile, nd_viewing_level                &
     &      , nd_direction, nd_channel)
!           Radiances calculated
!
!     Arrays for use with diagnostics:
      real                                                              &
     &    weight_690nm(sw_spectrum%npd_band)                            &
!           Weights for each band for region below 690 nm
     &  , weight_uv(sw_spectrum%npd_band)
!           Weight for each band for regin in UV
      real                                                              &
     &    flux_direct_blue_surf(nd_flux_profile)                        &
!           Total direct blue flux at the surface
     &  , flux_down_blue_surf(nd_flux_profile)                          &
!           Total downward blue flux at the surface
        , flux_up_blue_surf(nd_flux_profile)                            &
!           Upward blue flux at the surface
        , flux_down_uv_surf(nd_flux_profile)                            &
!           Total downward UV flux at the surface
        , flux_down_clr_uv_surf(nd_flux_profile)
!           Clear-sky downward UV flux at the surface
!
!     Local fields for cloud diagnostics:
      real                                                              &
     &    cloud_extinction_g(nd_flux_profile, nd_layer)                 &
!           Mean extinction coefficient in clouds weighted by the
!           cloud amount and the clear-sky flux (gathered array)
     &  , cloud_weight_extinction_g(nd_flux_profile, nd_layer)          &
!           Weighting factor for extinction in clouds: the product
!           of the cloud amount and the clear-sky direct flux
!           (gathered array)
     &  , ls_cloud_extinction_g(nd_flux_profile, nd_layer)              &
!           Mean extinction coefficient in clouds weighted by the
!           cloud amount and the clear-sky flux (gathered array)
     &  , ls_cloud_weight_extinction_g(nd_flux_profile, nd_layer)       &
!           Weighting factor for extinction in clouds: the product
!           of the cloud amount and the clear-sky direct flux
!           (gathered array)
     &  , cnv_cloud_extinction_g(nd_flux_profile, nd_layer)             &
!           Mean extinction coefficient in clouds weighted by the
!           cloud amount and the clear-sky flux (gathered array)
     &  , cnv_cloud_weight_extinction_g(nd_flux_profile, nd_layer)
!           Weighting factor for extinction in clouds: the product
!           of the cloud amount and the clear-sky direct flux
!           (gathered array)

!     Small real number used in COSP diagnostics calculations
      REAL :: epsreal

!     Tiled surface fluxes
      real                                                              &
     &    flux_up_tile(nd_point_tile, nd_tile, nd_channel)              &
!           Spectrally integrated upward flux on tiles
     &  , flux_up_blue_tile(nd_point_tile, nd_tile, nd_channel)
!           Upward blue flux on the tiles
!
!     Temporary field associated with the orography correction.
      real                                                              &
     &     swout_temp(nd_field)
!
!     Fields required for call to radiation code but not used
      INTEGER, Parameter :: nd_j_profile=1
      INTEGER :: i_gas
!
!     Auxiliary variables:
      real                                                              &
     &    weight_band(sw_spectrum%npd_band)
!           Weighting factors for bands
!
!     Variables required for compatibility with subroutines:
      real                                                              &
     &     dummy                                                        &
     &    ,dummy1d(1)                                                   &
     &    ,dummy2d(1,1)                                                 &
     &    ,dummy3d(1,1,1)

!     Switches for cloud diagnostics:
      LOGICAL :: l_ls_cloud_extinction, l_cnv_cloud_extinction

      CHARACTER (LEN=*), PARAMETER :: RoutineName = 'r2_swrad3z'
      CHARACTER (LEN=240) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
!
      IF (lhook) CALL dr_hook('R2_SWRAD3Z',zhook_in,zhook_handle)

!     Initialize the error flag for the radiation code.
      ierr=i_normal
      condensed_dim_char(:,:,:) = 0.0
      
!     Initialise logicals which are used later on.
      l_ls_cloud_extinction  = .FALSE.
      l_cnv_cloud_extinction = .FALSE.


      if (nlit >  0) then
!
!       Initializations for diagnostics depending on bands
        IF (sw_spectrum%weight_blue(1) == rmdi) THEN
! DEPENDS ON: r2_set_690nm_weight
          call r2_set_690nm_weight(sw_spectrum%n_band                   &
     &      , sw_spectrum%l_present                                     &
     &      , sw_spectrum%n_band_exclude                                &
     &      , sw_spectrum%index_exclude                                 &
     &      , sw_spectrum%wave_length_short                             &
     &      , sw_spectrum%wave_length_long                              &
     &      , weight_690nm                                              &
     &      , sw_spectrum%npd_band, sw_spectrum%npd_exclude             &
     &      , sw_spectrum%npd_type                                      &
     &      )
        ELSE
          weight_690nm(1:sw_spectrum%n_band) =                          &
            sw_spectrum%weight_blue(1:sw_spectrum%n_band)
        END IF

        IF (sw_diag%l_uvflux_direct .OR.                                &
            sw_diag%l_uvflux_up .OR.                                    &
            sw_diag%l_uvflux_net .OR.                                   &
            sw_diag%l_surf_uv .OR.                                      &
            sw_diag%l_surf_uv_clr) THEN

! DEPENDS ON: r2_set_uv_weight
          call r2_set_uv_weight(sw_spectrum%n_band                      &
     &      , sw_spectrum%l_present                                     &
     &      , sw_spectrum%n_band_exclude                                &
     &      , sw_spectrum%index_exclude                                 &
     &      , sw_spectrum%wave_length_short                             &
     &      , sw_spectrum%wave_length_long                              &
     &      , weight_uv                                                 &
     &      , sw_spectrum%npd_band, sw_spectrum%npd_exclude             &
     &      , sw_spectrum%npd_type                                      &
     &      )
        end if
!
!       Set clear-sky calculations.
        l_clear = sw_diag%l_flux_up_clear                               &
         .OR. sw_diag%l_flux_down_clear                                 &
         .OR. sw_diag%l_solar_out_clear                                 &
         .OR. sw_diag%l_surf_down_clr                                   &
         .OR. sw_diag%l_surf_up_clr                                     &
         .OR. sw_diag%l_clear_hr                                        &
         .OR.(sw_diag%l_cloud_extinction.AND.                           &
              sw_diag%l_cloud_weight_extinction)                        &
         .OR.(sw_diag%l_ls_cloud_extinction.AND.                        &
              sw_diag%l_ls_cloud_weight_extinction)                     &
         .OR.(sw_diag%l_cnv_cloud_extinction.AND.                       &
              sw_diag%l_cnv_cloud_weight_extinction).OR.l_cosp
!
!       Set dynamic array sizes and dependent options.
! DEPENDS ON: r2_set_option
        call r2_set_option(ierr                                         &
     &    , n_layer, sw_spectrum                                        &
     &    , sw_control%isolir                                           &
     &    , sw_control%i_gas_overlap, i_gas_overlap                     &
     &    , sw_control%l_aerosol, l_climat_aerosol                      &
     &    , l_use_sulpc_direct, l_use_soot_direct, l_use_biogenic       &
     &    , l_use_dust, l_use_bmass_direct, l_use_ocff_direct           &
     &    , l_use_nitrate_direct                                        &
     &    , l_use_seasalt_direct, l_murk_rad, l_aerosol                 &
     &    , l_use_sulpc_indirect, l_aerosol_ccn, n_arcl_species         &
     &    , sw_control%l_global_cloud_top                               &
     &    , global_cloud_top, n_cloud_top_global                        &
          , l_clear.OR.sw_diag%l_surf_uv_clr                            &
     &    , sw_control%i_angular_integration                            &
     &    , sw_control%i_solver, i_solver_clear                         &
     &    , sw_control%l_rescale, n_order_forward                       &
     &    , sw_control%i_truncation, sw_control%ls_global_trunc         &
     &    , sw_control%l_euler_trnf, euler_factor                       &
     &    , sw_control%i_scatter_method, i_scatter_method_band          &
     &    , l_rad_tile                                                  &
     &    , weight_band                                                 &
     &    , nd_overlap_coeff, nd_2sg_profile, nd_layer_clr              &
     &    , nd_source_coeff, nd_max_order, nd_sph_coeff                 &
     &    , nd_region                                                   &
     &    , nd_profile                                                  &
     &    )
        if (ierr /= i_normal) then 
          cmessage = 'Error following call to r2_set_option_sw'
          GO TO 9999
        end if
!
!       If radiances are to be calculated the viewing directions
!       are to be copied into the gathered array.
        if ( (sw_control%i_angular_integration ==                       &
     &        ip_spherical_harmonic).and.                               &
     &       (sw_control%i_sph_mode == ip_sph_mode_rad) ) then
          do ll=1, nlit
            l=list(ll)
            viewing_direction_g(ll, 1, 1)=viewing_direction1(l, 1)
            viewing_direction_g(ll, 1, 2)=viewing_direction2(l, 1)
          end do
        end if
!
!
!
!       Set the properties of the surface
! DEPENDS ON: r2_set_surface_field_sw
        call r2_set_surface_field_sw(ierr                               &
     &    , sw_spectrum%n_band, sw_control%ls_brdf_trunc                &
     &    , nlit, list                                                  &
     &    , l_ctile, l_use_spec_sea                                     &
     &    , land, land0p5, open_sea_albedo                              &
     &    , sea_ice_albedo                                              &
     &    , flandg, ice_fraction                                        &
     &    , land_albedo, weight_690nm                                   &
     &    , land0p5_g, flandg_g                                         &
     &    , n_brdf_basis_fnc, f_brdf, rho_alb                           &
     &    , l_rad_tile, n_point_tile, n_tile, list_tile, rho_alb_tile   &
     &    , list_tile_outer, index_tile                                 &
     &    , nd_field, nd_profile, sw_spectrum%npd_band                  &
     &    , nd_brdf_basis_fnc, nd_brdf_trunc                            &
     &    , nd_point_tile, nd_tile                                      &
     &    )
! and write out any surface diagnostics 
        IF (sw_diag%l_direct_albedo) THEN 
          DO k=1, sw_spectrum%n_band 
            DO ll=1, nlit 
              !l=list(ll) 
              sw_diag%direct_albedo(col_list(ll),row_list(ll),k)        &
                    = rho_alb(ll,ip_surf_alb_dir,k) 
            END DO 
          END DO 
        END IF 
        IF (sw_diag%l_diffuse_albedo) THEN 
          DO k=1, sw_spectrum%n_band 
            DO ll=1, nlit 
              !l=list(ll) 
              sw_diag%diffuse_albedo(col_list(ll),row_list(ll),k)       &
                    = rho_alb(ll,ip_surf_alb_diff,k) 
            END DO 
          END DO 
        END IF 

!
!       Set the mixing ratios of gases.
! DEPENDS ON: r2_set_gas_mix_ratio
        call r2_set_gas_mix_ratio(ierr                                  &
     &    , nlit, nlevs, n_layer, nwet, nozone                          &
     &    , list, sw_control%l_extra_top                                &
     &    , sw_spectrum%n_absorb, sw_spectrum%type_absorb               &
     &    , sw_control%l_n2o, sw_control%l_ch4                          &
     &    , .false., .false., sw_control%l_o2                           &
     &    , .false., .false., .false., .false., .false.                 &
     &    , h2o, co2, o3, n2o_mix_ratio, ch4_mix_ratio                  &
     &    , nullmmr, nullmmr, o2_mix_ratio                              &
     &    , nullmmr, nullmmr, nullmmr, nullmmr, nullmmr                 &
     &    , gas_mix_ratio                                               &
     &    , co2_dim1, co2_dim2, co2_3d, l_co2_3d                        &
     &    , nd_field, nd_profile, nd_layer, sw_spectrum%npd_species, 1  &
! Note that ozone is passed directly and all the other fields are 
! ignored here, so for the number of greenhouse gases we supply zero.
     &    ,0, dummy3d                                                   &
     &    )
        if (ierr /= i_normal) then 
          cmessage = 'Error following call to r2_set_gas_mix_ratio'
          GO TO 9999
        end if
!
!       Set the thermodynamic properties of the atmosphere.
! DEPENDS ON: r2_set_thermodynamic
        call r2_set_thermodynamic(nlit, nlevs, n_layer, nwet, list      &
     &    , sw_control%l_extra_top, .false.                             &
          , pstar                                                       &
     &    , p_layer_boundaries                                          &
     &    , p_layer_centres                                             &
     &    , dummy2d                                                     &
     &    , dummy2d, tac                                                &
     &    , rho_r2, r_rho_levels, r_theta_levels                        &
     &    , q, qcl, qcf, qcf2, qrain, qgraup                            &
     &    , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup, l_mixing_ratio       &
          , p, t, dummy2d, d_mass                                       &
     &    , layer_heat_capacity                                         &
     &    , nd_field, nd_profile, nd_layer                              &
     &    )
!
!
!       set sea-salt array dimensions.
        if (l_use_seasalt_direct) then
           salt_dim_dir_a=salt_dim_a
           salt_dim_dir_b=salt_dim_b
        else
           salt_dim_dir_a=1
           salt_dim_dir_b=1
        end if
!
        DO l=1, nlit
          alat(l) = Recip_Pi_Over_180 * true_latitude(list(l))
        END DO

!
!       Set the mixing ratios of aerosols.
        if (l_aerosol.or.l_aerosol_ccn) then
! DEPENDS ON: r2_set_aerosol_field
          call r2_set_aerosol_field(ierr                                &
     &      , nlit, nlevs, n_layer, sw_spectrum%n_aerosol               &
     &      , sw_spectrum%n_aerosol_mr, sw_spectrum%type_aerosol        &
     &      , list, sw_control%l_extra_top                              &
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
     &      , land0p5, lying_snow, pstar                                &
     &      , p_layer_boundaries, trindx, alat, PREVIOUS_TIME           &
     &      , aerosol_mix_ratio, aerosol_mr_source,aerosol_mr_type_index&
     &      , nd_field, nd_profile, nd_layer                            &
     &      , sw_spectrum%npd_aerosol_species                           &
     &      , sw_spectrum%npd_aerosol_mixratio, 1                       &
     &      )
           if (ierr /= i_normal) then 
             cmessage = 'Error following call to r2_set_aerosol_field'
             GO TO 9999
           end if
        end if
!
!       Set mixing ratios of UKCA aerosols 
        IF (l_ukca_radaer) THEN 
! DEPENDS ON: ukca_radaer_set_aerosol_field 
         CALL ukca_radaer_set_aerosol_field(                            & 
              list, sw_control%l_extra_top, n_layer, nlit               & 
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
!       Assign the properties of clouds.
!
!       check the consistency of cloud diagnostics.
        if (sw_diag%re_conv_flag) then
          if (.not.sw_diag%wgt_conv_flag) then
            cmessage =                                                  &
              '*** error: microphysical diagnostics for convective ' // &
              'cloud must include the cloud weighting.'
            ierr=i_err_fatal
            GO TO 9999
          end if
        end if
!
        if ( (sw_diag%re_strat_flag).or.(sw_diag%lwp_strat_flag) ) then
          if (.not.sw_diag%wgt_strat_flag) then
            cmessage =                                                  &
              '*** error: microphysical diagnostics for stratiform ' // &
              'cloud must include the cloud weighting.'
            ierr=i_err_fatal
            GO TO 9999
          end if
        end if
!
        if (sw_diag%l_cloud_extinction) then
          if (.not.sw_diag%l_cloud_weight_extinction) then
            cmessage =                                                  &
              '*** error: the cloud extinction ' //                     &
              'may be diagnosed only in conjunction ' //                &
              'with the corresponding weights.'
            ierr=i_err_fatal
            GO TO 9999
          end if
        end if
!
        if (sw_diag%l_ls_cloud_extinction) then
          if (.not.sw_diag%l_ls_cloud_weight_extinction) then
            cmessage =                                                  &
              '*** error: the layer cloud extinction ' //               &
              'may be diagnosed only in conjunction ' //                &
              'with the corresponding weights.'
            ierr=i_err_fatal
            GO TO 9999
          end if
        end if
!
        if (sw_diag%l_cnv_cloud_extinction) then
          if (.not.sw_diag%l_cnv_cloud_weight_extinction) then
            cmessage =                                                  &
              '*** error: the conv. cloud extinction ' //               &
              'may be diagnosed only in conjunction ' //                &
              'with the corresponding weights.'
            ierr=i_err_fatal
            GO TO 9999
          end if
        end if

!
!    Set flags for cloud diagnostics
        IF (SW_diag%L_ls_cloud_extinction  .OR. L_cosp)  &
          l_ls_cloud_extinction = .TRUE.
        IF (SW_diag%L_cnv_cloud_extinction .OR. L_cosp)  &
          l_cnv_cloud_extinction = .TRUE.

!
! DEPENDS ON: r2_set_cloud_parametrization
        call r2_set_cloud_parametrization(ierr, sw_spectrum%n_band      &
     &    , sw_control%i_st_water, sw_control%i_cnv_water               &
     &    , sw_control%i_st_ice, sw_control%i_cnv_ice                   &
     &    , sw_spectrum%l_drop_type                                     &
     &    , sw_spectrum%i_drop_parametrization                          &
     &    , sw_spectrum%n_drop_phf_term, sw_spectrum%drop_parameter_list&
     &    , sw_spectrum%drop_parm_min_dim                               &
     &    , sw_spectrum%drop_parm_max_dim, sw_spectrum%l_ice_type       &
     &    , sw_spectrum%i_ice_parametrization                           &
     &    , sw_spectrum%n_ice_phf_term, sw_spectrum%ice_parameter_list  &
     &    , sw_spectrum%ice_parm_min_dim, sw_spectrum%ice_parm_max_dim  &
     &    , i_condensed_param, condensed_n_phf, condensed_param_list    &
     &    , condensed_min_dim, condensed_max_dim                        &
     &    , sw_spectrum%npd_band                                        &
     &    , sw_spectrum%npd_drop_type, sw_spectrum%npd_ice_type         &
     &    , sw_spectrum%npd_cloud_parameter                             &
     &    , nd_cloud_component                                          &
     &    )
        if (ierr /= i_normal) then 
          cmessage = 'Error following call to r2_set_cloud_parametrization'
          GO TO 9999
        end if
!
!       set sea-salt array dimensions.
        if (l_use_seasalt_indirect) then
          salt_dim_ind_a=salt_dim_a
          salt_dim_ind_b=salt_dim_b
        else
          salt_dim_ind_a=1
          salt_dim_ind_b=1
        end if
!
!
! DEPENDS ON: r2_set_cloud_field
        call r2_set_cloud_field(nlit, nlevs, n_layer, nclds             &
     &    , list                                                        &
     &    , p, t, d_mass, alat                                          &
     &    , ccb, cct, cca, cccwp, ccw, lcbase                           &
          , lccwc1, lccwc2, lca_area, lca_bulk, n_drop_pot              &
          , l_pc2, sw_control%l_microphysics, l_aerosol_ccn             &
     &    , sea_salt_film, sea_salt_jet                                 &
     &    , l_use_seasalt_indirect, salt_dim_ind_a, salt_dim_ind_b      &
     &    , l_use_biogenic, biogenic, biogenic_dim1, biogenic_dim2      &
     &    , sulp_dim1, sulp_dim2, accum_sulphate, diss_sulphate         &
     &    , aitken_sulphate, L_USE_BMASS_INDIRECT                       &
     &    , bmass_dim1, bmass_dim2, aged_bmass, cloud_bmass             &
     &    , l_use_ocff_indirect, ocff_dim1, ocff_dim2                   &
     &    , aged_ocff, cloud_ocff                                       &
     &    , l_use_nitrate_indirect, nitrate_dim1, nitrate_dim2          &
     &    , accum_nitrate, diss_nitrate                                 &
     &    , lying_snow                                                  &
     &    , land0p5_g, flandg_g                                         &
     &    , sw_control%i_cloud_representation, i_condensed_param        &
     &    , condensed_min_dim, condensed_max_dim                        &
     &    , n_condensed, type_condensed                                 &
     &    , w_cloud, n_cloud_type, frac_cloud                           &
     &    , sw_control%l_local_cnv_partition                            &
     &    , condensed_mix_ratio, condensed_dim_char                     &
     &    , sw_diag%re_conv, sw_diag%re_conv_flag                       &
     &    , sw_diag%re_strat, sw_diag%re_strat_flag                     &
     &    , sw_diag%wgt_conv, sw_diag%wgt_conv_flag                     &
     &    , sw_diag%wgt_strat, sw_diag%wgt_strat_flag                   &
     &    , sw_diag%lwp_strat, sw_diag%lwp_strat_flag                   &
     &    , sw_diag%ntot_diag, sw_diag%ntot_diag_flag                   &
     &    , sw_diag%strat_lwc_diag, sw_diag%strat_lwc_diag_flag         &
     &    , sw_diag%so4_ccn_diag, sw_diag%so4_ccn_diag_flag             &
     &    , sw_diag%cond_samp_wgt, sw_diag%cond_samp_wgt_flag           &
     &    , sw_diag%nc_diag, sw_diag%nc_diag_flag                       &
     &    , sw_diag%nc_weight, sw_diag%nc_weight_flag                   &
     &    , col_list, row_list, row_length, rows                        &
     &    , nd_field, nd_profile, nd_layer                              &
     &    , sw_spectrum%npd_aerosol_species, 1, id_ct, id_ct            &
     &    , n_cca_lev, Ntot_land, Ntot_sea                              &
     &    )

        select case (sw_control%i_cloud)
        case (ip_cloud_mcica)
           i_cloud_tmp=ip_cloud_mix_max
        case (ip_cloud_part_corr)
           i_cloud_tmp=ip_cloud_mix_max
        case (ip_cloud_part_corr_cnv)
           i_cloud_tmp=ip_cloud_triple
        case default
           i_cloud_tmp=sw_control%i_cloud
        end select

        if (sw_diag%weighted_re_flag.and.                               &
     &      sw_diag%sum_weight_re_flag) then
! DEPENDS ON: r2_cloud_level_diag
          call r2_cloud_level_diag(ierr, nlit, n_layer, nclds           &
     &      , list                                                      &
     &      , i_cloud_tmp, sw_control%i_cloud_representation            &
     &      , t, w_cloud, frac_cloud, .true.                            &
     &      , condensed_mix_ratio, condensed_dim_char                   &
     &      , sw_diag%weighted_re_flag, sw_diag%weighted_re             &
     &      , sw_diag%sum_weight_re                                     &
     &      , col_list, row_list, row_length, rows                      &
     &      , nd_field, nd_profile, nd_layer, 1, id_ct, id_ct           &
     &      )
          if (ierr /= i_normal) then 
            cmessage = 'Error following call to r2_cloud_level_diag'
            GO TO 9999
          end if
        end if
!
        if (sw_diag%wgtd_warm_re_flag.and.                              &
     &      sw_diag%sum_wgt_warm_re_flag) then
          call r2_cloud_level_diag(ierr, nlit, n_layer, nclds           &
     &      , list                                                      &
     &      , i_cloud_tmp, sw_control%i_cloud_representation            &
     &      , t, w_cloud, frac_cloud, .false.                           &
     &      , condensed_mix_ratio, condensed_dim_char                   &
     &      , sw_diag%wgtd_warm_re_flag, sw_diag%weighted_warm_re       &
     &      , sw_diag%sum_weight_warm_re                                &
     &      , col_list, row_list, row_length, rows                      &
     &      , nd_field, nd_profile, nd_layer, 1, id_ct, id_ct           &
     &      )
          if (ierr /= i_normal) then 
            cmessage = 'Error following call to r2_cloud_level_diag'
            GO TO 9999
          end if
        end if
!


        IF (sw_control%i_cloud==IP_cloud_mcica) THEN
          ALLOCATE(subcol_k(sw_spectrum%npd_band                        &
                    ,sw_spectrum%npd_esft_term))
          ALLOCATE(subcol_reorder(subcol_need))
          DO ll=1, sw_spectrum%npd_esft_term
            DO l=1, sw_spectrum%n_band
              subcol_k(l,ll)=sw_subcol_k(l,ll)
            END DO
          END DO
          DO l=1,subcol_need
            subcol_reorder(l)=sw_subcol_reorder(l)
          END DO
! DEPENDS ON: mcica_order
          CALL mcica_order(ierr                                         &
!                       General spectral properties
             , sw_control%first_band, sw_spectrum%n_band                &
!                       Gaseous absorption
             , sw_spectrum%i_band_esft, sw_spectrum%index_absorb        &
!                       Dimensions of arrays
             , sw_spectrum%npd_esft_term, sw_spectrum%npd_band          &
             , sw_spectrum%npd_species                                  &
         )

          ALLOCATE(c_sub(nd_profile, id_ct:nd_layer, tot_subcol_gen     &
            , n_cloud_type))
          ALLOCATE(frac_cloudy(nd_profile))
          DO l=1,nlit
            i=((row_list(l)-1)*row_length)+col_list(l)
            frac_cloudy(l)=frac_cloudy_full(i)
          END DO
          IF (ALLOCATED(cic_sub_full)) THEN
            DO lll=1,subcol_need
              DO ll=id_ct,n_layer
                DO l=1,nlit
                  i=((row_list(l)-1)*row_length)+col_list(l)
                  c_sub(l,ll,lll,1)=clw_sub_full(i,ll,lll)
                  c_sub(l,ll,lll,2)=cic_sub_full(i,ll,lll)
                END DO
              END DO
            END DO
          ELSE
            DO lll=1,subcol_need
              DO ll=id_ct,n_layer
                DO l=1,nlit
                  i=((row_list(l)-1)*row_length)+col_list(l)
                  c_sub(l,ll,lll,1)=clw_sub_full(i,ll,lll)
                  c_sub(l,ll,lll,2)=clw_sub_full(i,ll,lll)
                END DO
              END DO
            END DO
          END IF
        END IF

        IF (sw_control%i_cloud/=IP_cloud_mcica) THEN
          IF (ALLOCATED(cloud_inhom_param_full)) THEN
            ALLOCATE(cloud_inhom_param(nd_profile, id_ct:nd_layer))
            DO ll=id_ct,n_layer
              DO l=1,nlit
                i=((row_list(l)-1)*row_length)+col_list(l)
                cloud_inhom_param(l,ll)=cloud_inhom_param_full(i,ll)
              END DO
            END DO
          END IF
        END IF
        

!
        if (L_orog) allocate(lg_orog_corr(nlit))
!
!       Set the incident solar flux.
        do l=1, nlit
          solar_irrad(l)=scs*sc*lit(list(l))

!         Gather the orography correction factor
!         into lit points.
          if (L_orog) then
             lg_orog_corr(l)=orog_corr(col_list(l),row_list(l))
          end if

        end do
        if (sw_control%i_angular_integration == ip_two_stream) then
          do l=1, nlit
            zen_0(l)=1.0/coszin(list(l))
          end do
        else if (sw_control%i_angular_integration ==                    &
     &           ip_spherical_harmonic) then
          do l=1, nlit
            zen_0(l)=coszin(list(l))
          end do
        end if
!
!
!
!
! DEPENDS ON: radiance_calc
        call radiance_calc(ierr                                         &
!                       Logical flags for processes
     &    , sw_control%l_rayleigh, l_aerosol                            &
     &    , sw_control%l_gas, sw_control%l_continuum                    &
     &    , sw_control%l_cloud, sw_control%l_drop, sw_control%l_ice     &
!                       Angular integration
     &    , sw_control%i_angular_integration, sw_control%l_rescale      &
     &    , n_order_forward                                             &
     &    , sw_control%i_2stream, sw_control%n_order_gauss              &
     &    , sw_control%i_truncation, sw_control%ls_global_trunc         &
     &    , sw_control%ms_min, sw_control%ms_max                        &
     &    , sw_control%accuracy_adaptive, euler_factor                  &
          , sw_control%l_henyey_greenstein_pf, .FALSE.                  &
          , sw_control%ls_brdf_trunc                                    &
     &    , sw_control%i_sph_algorithm, sw_control%n_order_phase_solar  &
     &    , n_viewing_direction, viewing_direction_g                    &
     &    , n_viewing_level, viewing_level                              &
     &    , sw_control%i_sph_mode                                       &
!                       Treatment of scattering
     &    , i_scatter_method_band                                       &
!                       Options for treating clouds
     &    , sw_control%l_global_cloud_top, n_cloud_top_global           &
     &    , l_inhom_cloud, inhom_cloud                                  &
!                       Options for solver
     &    , sw_control%i_solver                                         &
!                       Properties of diagnostics
     &    , map_channel                                                 &
!                       General spectral properties
     &    , sw_spectrum%n_band                                          &
          , sw_control%first_band, sw_spectrum%n_band                   &
     &    , weight_band, sw_spectrum%l_present(14)                      &
     &    , sw_spectrum%n_band_exclude, sw_spectrum%index_exclude       &
!                       General atmospheric properties
     &    , nlit, n_layer                                               &
     &    , p, t, dummy, dummy, d_mass                                  &
!                       Spectral region
     &    , sw_control%isolir                                           &
!                       Solar fields
     &    , zen_0, solar_irrad, sw_spectrum%solar_flux_band             &
          , sw_spectrum%solar_flux_band_ses, .FALSE., solar_tail_flux   &
     &    , sw_spectrum%rayleigh_coefficient                            &
!                       Infra-red fields
     &    , sw_spectrum%n_deg_fit, sw_spectrum%thermal_coefficient      &
     &    , sw_spectrum%t_ref_planck, .false.                           &
!                       Gaseous absorption
     &    , i_gas_overlap, i_gas, gas_mix_ratio                         &
     &    , sw_spectrum%n_band_absorb, sw_spectrum%index_absorb         &
     &    , sw_spectrum%i_band_esft                                     &
     &    , sw_spectrum%w_esft, sw_spectrum%k_esft                      &
     &    , sw_spectrum%i_scale_esft, sw_spectrum%i_scale_fnc           &
     &    , sw_spectrum%scale_vector                                    &
     &    , sw_spectrum%p_reference, sw_spectrum%t_reference            &
     &    , l_mod_k_flux, sw_spectrum%mix_gas_band                      &
          , sw_spectrum%n_mix_gas, sw_spectrum%index_mix_gas            &
          , sw_spectrum%f_mix, sw_spectrum%i_band_esft_ses              &
          , sw_spectrum%k_esft_ses                                      &
          , sw_spectrum%k_mix_gas, sw_spectrum%w_esft_ses               &
!                       Doppler broadening
     &    , sw_spectrum%l_doppler_present                               &
     &    , sw_spectrum%doppler_correction                              &
!                       Surface fields
     &    , n_brdf_basis_fnc, rho_alb, f_brdf                           &
!                       Tiling options for heterogeneous surfaces
     &    , l_rad_tile, n_point_tile, n_tile, list_tile                 &
     &    , rho_alb_tile, dummy2d, t_tile                               &
!                       Continuum absorption
     &    , sw_spectrum%n_band_continuum                                &
     &    , sw_spectrum%index_continuum, sw_spectrum%index_water        &
     &    , sw_spectrum%k_continuum, sw_spectrum%i_scale_fnc_cont       &
     &    , sw_spectrum%scale_continuum                                 &
     &    , sw_spectrum%p_ref_continuum, sw_spectrum%t_ref_continuum    &
          , sw_spectrum%k_continuum_ses, sw_spectrum%k_h2oc             &
!                       Properties of aerosols
     &    , sw_spectrum%n_aerosol, sw_spectrum%n_aerosol_mr             &
     &    , aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index &
     &    , l_use_arcl                                                  &
     &    , sw_spectrum%aerosol_absorption                              &
     &    , sw_spectrum%aerosol_scattering                              &
     &    , sw_spectrum%n_aerosol_phf_term                              &
     &    , sw_spectrum%aerosol_phase_fnc                               &
     &    , sw_spectrum%i_aerosol_parametrization                       &
     &    , sw_spectrum%nhumidity, sw_spectrum%humidities               &
     &    , sw_spectrum%type_aerosol                                    &
!                       Aerosol optical depth
!                       (computed in the LW: dummy values for SW)
     &    , sw_spectrum%n_aod_wavel, sw_spectrum%aod_wavel              & 
     &    , .false., dummy2d, .false., dummy2d, .false., dummy2d        &
     &    , .false., dummy2d, .false., dummy2d, .false., dummy2d        &
     &    , .false., dummy2d, .false., dummy2d, .false., dummy2d        &
     &    , .false., dummy2d, .false., dummy2d, .false., dummy2d        &
     &    , .false., dummy2d, .false., dummy2d, .false., dummy2d        &
     &    , .false., dummy2d, .false., dummy2d, .false., dummy2d        &    
     &    , sw_spectrum%aod_absorption                                  &
     &    , sw_spectrum%aod_scattering, sw_spectrum%i_aod_type          &
!                       Properties of UKCA aerosols 
     &    , l_ukca_radaer, ukca_radaer, ukca_mix_ratio, ukca_comp_vol   & 
     &    , ukca_dry_diam, ukca_wet_diam, ukca_modal_rho, ukca_modal_vol& 
     &    , ukca_modal_wtv, ukca_modal_nbr                              & 
!                       Optical depth of UKCA aerosols
!                       (computed in the LW: dummy values for SW)
     &    , .false., dummy2d, .false., dummy2d, .false., dummy2d        &
     &    , .false., dummy2d, .false., dummy2d, .false., dummy2d        &
!                       Properties of clouds
     &    , n_condensed, type_condensed                                 &
     &    , sw_control%i_cloud, sw_control%i_cloud_representation       &
     &    , w_cloud, n_cloud_type, frac_cloud, tot_cloud_cover          &
     &    , condensed_mix_ratio, condensed_dim_char                     &
     &    , i_condensed_param, condensed_n_phf, condensed_param_list    &
     &    , dp_corr_strat, dp_corr_conv                                 &
!                       Calculated Fluxes
     &    , flux_direct, flux_net, flux_up                              &
     &    , uv_flux_direct, uv_flux_net, uv_flux_up                     &
     &    , SW_diag%L_uvflux_direct                                     &
     &    , SW_diag%L_uvflux_net,SW_diag%L_uvflux_up                    &
!                       Calculated Radiances
     &    , radiance                                                    &
!                       Photolysis
     &    , dummy2d                                                     &
!                       Options for clear-sky fluxes
     &    , l_clear, i_solver_clear                                     &
!                       Clear-sky fluxes calculated
     &    , flux_direct_clear, flux_net_clear, flux_up_clear            &
!                       Special Surface Fluxes
     &    , l_flux_below_690nm_surf, weight_690nm,weight_uv             &
     &    , flux_direct_blue_surf, flux_down_blue_surf                  &
     &    , flux_up_blue_surf                                           &
          , SW_diag%l_surf_uv, flux_down_uv_surf                        &
          , SW_diag%l_surf_uv_clr, flux_down_clr_uv_surf                &
!                       Tiled Surface Fluxes
     &    , flux_up_tile, flux_up_blue_tile                             &
!                       Arrays for diagnostics specific to the UM
     &    , SW_diag%L_cloud_extinction, cloud_extinction_g              &
     &    , cloud_weight_extinction_g                                   &
          , l_ls_cloud_extinction, ls_cloud_extinction_g                &
     &    , ls_cloud_weight_extinction_g                                &
          , l_cnv_cloud_extinction, cnv_cloud_extinction_g              &
     &    , cnv_cloud_weight_extinction_g                               &
     &    , .false., dummy2d, dummy2d, .false., dummy2d, dummy2d        &
     &    , .false., dummy2d, dummy2d                                   &
!                       Dimensions of arrays
     &    , nd_profile, nd_layer, nd_column, nd_layer_clr, id_ct        &
     &    , nd_2sg_profile, nd_flux_profile, nd_radiance_profile        &
     &    , nd_j_profile, nd_channel                                    &
     &    , sw_spectrum%npd_band, sw_spectrum%npd_species               &
     &    , sw_spectrum%npd_esft_term                                   &
     &    , sw_spectrum%npd_scale_variable                              &
     &    , sw_spectrum%npd_continuum                                   &
     &    , sw_spectrum%npd_aerosol_species                             &
     &    , sw_spectrum%npd_aerosol_mixratio                            &
     &    , sw_spectrum%npd_humidities                                  &
     &    , sw_spectrum%npd_cloud_parameter                             &
     &    , sw_spectrum%npd_thermal_coeff, nd_source_coeff              &
     &    , nd_brdf_basis_fnc, nd_brdf_trunc                            &
     &    , sw_spectrum%npd_aod_wavel                                   &
     &    , sw_spectrum%npd_phase_term                                  &
     &    , nd_max_order, nd_sph_coeff                                  &
     &    , nd_direction, nd_viewing_level                              &
     &    , nd_region, nd_cloud_type, nd_cloud_component                &
     &    , nd_overlap_coeff                                            &
     &    , nd_point_tile, nd_tile                                      &
          , sw_spectrum%npd_tmp, sw_spectrum%npd_pre                    &
          , sw_spectrum%npd_mix, sw_spectrum%npd_band_mix_gas           &
     &    , n_ukca_mode, n_ukca_cpnt, sw_spectrum%npd_exclude           &
     &    )
        if (ierr /= i_normal) then 
          cmessage = 'Error following call to radiance_calc'
          GO TO 9999
        end if
!
!
      end if

!
!     Prepare the output arrays:
!
!     Processing depends on whether the code has been invoked to
!     calculate radiances or fluxes.
      if ( (sw_control%i_angular_integration == ip_two_stream).or.      &
     &     ( (sw_control%i_angular_integration ==                       &
     &        ip_spherical_harmonic).and.                               &
     &        (sw_control%i_sph_mode == ip_sph_mode_flux) ) ) then


!       Total upward & downward all-sky & clear-sky flux at all levels:
!       (Note at this point "flux_net" holds the downward flux.)
        IF (SW_diag%l_flux_up) THEN
          DO i=1, nlevs+1
            DO l=1, nlit
              SW_diag%flux_up(col_list(l), row_list(l),i)               &
                = flux_up(l,n_layer+1-i, 1)
            END DO
          END DO
        END IF
        IF (SW_diag%l_flux_down) THEN
          DO i=1, nlevs+1
            DO l=1, nlit
              SW_diag%flux_down(col_list(l), row_list(l),i)             &
                = flux_net(l,n_layer+1-i, 1)
            END DO
          END DO
        END IF
        IF (SW_diag%l_flux_up_clear) THEN
          DO i=1, nlevs+1
            DO l=1, nlit
              SW_diag%flux_up_clear(col_list(l), row_list(l),i)         &
                = flux_up_clear(l,n_layer+1-i, 1)
            END DO
          END DO
        END IF
        IF (SW_diag%l_flux_down_clear) THEN
          DO i=1, nlevs+1
            DO l=1, nlit
              SW_diag%flux_down_clear(col_list(l), row_list(l),i)       &
                = flux_net_clear(l,n_layer+1-i, 1)
            END DO
          END DO
        END IF


!       Convert downward fluxes to net fluxes.
        do i=0, n_layer
          do l=1, nlit
            flux_net(l, i, 1)=flux_net(l, i, 1)-flux_up(l, i, 1)
          end do
        end do
        if (l_clear) then
          do i=0, n_layer
            do l=1, nlit
              flux_net_clear(l, i, 1)                                   &
     &          =flux_net_clear(l, i, 1)-flux_up_clear(l, i, 1)
            end do
          end do
        end if
!
!
!       Scatter the net downward flux at each level into SWOUT.
        do i=1, nlevs+1
          do l=1, nlit
            swout(list(l), i)=flux_net(l, n_layer+1-i, 1)
          end do
        end do
!
!
!       Net shortwave radiation absorbed by the planet
!       (i. e. earth and atmosphere together):
!
        do l=1, nlit
          netsw(list(l))=swout(list(l), nlevs+1)
        end do
!
        if (sw_control%l_extra_top) then
!         calculate the radiation absorbed in the extra layer
!         above the top of the rest of the model.
          do l=1, nlit
            top_absorption(list(l))=flux_net(l, 0, 1)                   &
     &        -flux_net(l, n_layer-nlevs, 1)
          end do
        end if
!
!
!       Extra direct SW flux reaching the surface due to the
!       orography correction:

        if (L_orog) then
           allocate(lg_f_orog(nlit))
        
           lg_f_orog = flux_direct(1:nlit, n_layer,1) *                 &
     &             (lg_orog_corr - 1.0)/lg_orog_corr
        
           do l=1, nlit
             f_orog(col_list(l),row_list(l)) = lg_f_orog(l)
           end do
        end if
!
!       Assignment of diagnostics
!
!       Note: purely diagnostic quantities allocated dynamically in
!       RAD_CTL2 are zeroed there and need to be filled only at lit
!       points.
!
!
!       Outgoing solar radiation at TOA:
!
        if (SW_diag%l_solar_out_toa) then
          do l=1, nlit
            SW_diag%solar_out_toa(col_list(l), row_list(l))             &
              = flux_up(l, 0, 1)
          end do
        end if
!
!
!       Clear-sky outgoing solar radiation at TOA:
!
        if (SW_diag%l_solar_out_clear) then
          do l=1, nlit
            SW_diag%solar_out_clear(col_list(l), row_list(l))           &
              = flux_up_clear(l, 0, 1)
          end do
        end if
!
!
!       Surface flux below 690nm.
!
        if (l_flux_below_690nm_surf) then
!
          if (l_rad_tile) then
            do ll=1, nlit
              l=list(ll)
              flux_below_690nm_surf(l)                                        &
                = flux_down_blue_surf(ll)-flux_up_blue_surf(ll)
              if ( (flandg(l) < TINY(flandg)) .AND.                           &
                   (ice_fraction(l) < TINY(ice_fraction)) ) then
!               This point is open sea with no sea ice.
                if (SW_diag%L_FlxSeaBelow690nmSurf) then
                  SW_diag%FlxSeaBelow690nmSurf(col_list(ll),row_list(ll))     &
                    = flux_below_690nm_surf(l)
                end if
              else
                if (SW_diag%L_FlxSolBelow690nmSurf) then
                  SW_diag%FlxSolBelow690nmSurf(col_list(ll),row_list(ll))     &
                    = flux_below_690nm_surf(l)
                end if
              end if
            end do

!           Tiled points will have both land and sea. Note that the
!           channel index of flux_up_tile is hard-wired to 1 because
!           we don't envisage calling the code in other cases.
            do lll=1, n_point_tile
              ll=list_tile(lll)
              l=list_tile_outer(lll)
              flux_below_690nm_surf(l)                                        &
                =(1.0-ice_fraction(l))*(flux_down_blue_surf(ll)               &
                -flux_up_blue_tile(lll, index_tile(ip_ocean_tile), 1))
              if (SW_diag%L_FlxSeaBelow690nmSurf) then
                SW_diag%FlxSeaBelow690nmSurf(col_list(ll), row_list(ll))      &
                  =flux_below_690nm_surf(l)
              end if
              if (SW_diag%L_FlxSolBelow690nmSurf) then
                frac_solid=flandg(l)+(1.0-flandg(l))*ice_fraction(l)
                if (frac_solid > 0.0) then
                  SW_diag%FlxSolBelow690nmSurf(col_list(ll),row_list(ll))     &
                    =(SW_diag%FlxSolBelow690nmSurf(col_list(ll),row_list(ll)) &
                    -(1.0-flandg(l))*flux_below_690nm_surf(l))/frac_solid
                end if
              end if
            end do
          else
            do ll=1, nlit
              l=list(ll)
              if (land(l)) then
                flux_below_690nm_surf(l)                                &
     &            =flux_down_blue_surf(ll)-flux_up_blue_surf(ll)
              else
                flux_below_690nm_surf(l)                                &
     &            =(flux_down_blue_surf(ll)-flux_up_blue_surf(ll))      &
     &            *(1.0-ice_fraction(l))
              end if
            end do
          end if
!
        end if
!
!
!       Orography correction to direct SW flux:
!
        if (SW_diag%L_orog_corr .AND. L_orog) then
           do l=1, nlit
              SW_diag%orog_corr(col_list(l), row_list(l))               &
     &          =lg_orog_corr(l)
           end do
        end if
!
!       We need to extend the core of the
!       radiation code to keep direct blue fluxes to get these
!       diagnostics right. For now we pass the total to the
!       diffuse.
          surf_down_sw(list(1:nlit), 1)                                 &
     &      =flux_direct_blue_surf(1:nlit)
          surf_down_sw(list(1:nlit), 2)                                 &
     &      =flux_down_blue_surf(1:nlit)-flux_direct_blue_surf(1:nlit)
          surf_down_sw(list(1:nlit), 3)                                 &
     &      =flux_direct(1:nlit, n_layer, 1)                            &
     &      -flux_direct_blue_surf(1:nlit)
          surf_down_sw(list(1:nlit), 4)                                 &
     &      =flux_net(1:nlit, n_layer, 1)+flux_up(1:nlit, n_layer, 1)   &
     &      -flux_down_blue_surf(1:nlit)                                &
     &      -flux_direct(1:nlit, n_layer, 1)                            &
     &      +flux_direct_blue_surf(1:nlit)
!
!
!       Downward flux at the surface:
!
        if (SW_diag%l_surface_down_flux) then
          do l=1, nlit
            SW_diag%surface_down_flux(col_list(l), row_list(l))         &
     &        =flux_net(l, n_layer, 1)+flux_up(l, n_layer, 1)
          end do
        end if
!
!
!       Clear-sky downward flux at the surface:
!
        if (SW_diag%l_surf_down_clr) then
          do l=1, nlit
            SW_diag%surf_down_clr(col_list(l), row_list(l))             &
     &        =flux_net_clear(l, n_layer, 1)                            &
     &        +flux_up_clear(l, n_layer, 1)
          end do
        end if
!
!
!       Clear-sky upward flux at the surface:
!
        if (SW_diag%l_surf_up_clr) then
          do l=1, nlit
            SW_diag%surf_up_clr(col_list(l), row_list(l))               &
     &        =flux_up_clear(l, n_layer, 1)
          end do
        end if
!
!
!       Net flux at the tropopause:
!
        if (SW_diag%l_net_flux_trop) then
          do l=1, nlit
            SW_diag%net_flux_trop(col_list(l), row_list(l))             &
     &        =flux_net(l, n_layer+1-trindx(list(l)), 1)
          end do
        end if
!
!
!       Upward flux at the tropopause:
!
        if (SW_diag%l_up_flux_trop) then
          do l=1, nlit
            SW_diag%up_flux_trop(col_list(l), row_list(l))              &
     &        =flux_up(l, n_layer+1-trindx(list(l)), 1)
          end do
        end if

!       Direct and diffuse downward flux
        IF (SW_diag%l_flux_direct) THEN
          DO i=1, nlevs+1
            DO l=1, nlit
              SW_diag%flux_direct(col_list(l), row_list(l),i)           &
                = flux_direct(l,n_layer+1-i,1)
            END DO
          END DO
        END IF
        IF (SW_diag%l_flux_diffuse) THEN
          DO i=1, nlevs+1
            DO l=1, nlit
              SW_diag%flux_diffuse(col_list(l), row_list(l),i)          & 
                = flux_net(l,n_layer+1-i,1) + flux_up(l,n_layer+1-i,1)  &
                - flux_direct(l,n_layer+1-i,1)
            END DO
          END DO
        END IF     

!       UV-Fluxes
        IF (SW_diag%l_uvflux_direct) THEN
          DO i=1, nlevs+1
            DO l=1, nlit
              SW_diag%uvflux_direct(col_list(l), row_list(l),i)         &
                = uv_flux_direct(l,n_layer+1-i,1)
            END DO
          END DO
        END IF
        IF (SW_diag%l_uvflux_up) THEN
          DO i=1,nlevs+1
            DO l=1, nlit
              SW_diag%uvflux_up(col_list(l), row_list(l),i)             &
                = uv_flux_up(l,n_layer+1-i,1)
            END DO
          END DO
        END IF
        IF (SW_diag%l_uvflux_net) THEN
          DO i=1,nlevs+1
            DO l=1, nlit
              SW_diag%uvflux_net(col_list(l), row_list(l),i)            &
                = uv_flux_net(l,n_layer+1-i,1)
            END DO
          END DO
        END IF
        IF (SW_diag%l_surf_uv) THEN
          DO l=1, nlit
            SW_diag%surf_uv(col_list(l), row_list(l))                   &
              = flux_down_uv_surf(l)
          END DO
        END IF
        IF (SW_diag%l_surf_uv_clr) THEN
          DO l=1, nlit
            SW_diag%surf_uv_clr(col_list(l), row_list(l))               &
              = flux_down_clr_uv_surf(l)
          END DO
        END IF
!
!    Cloud extinction diagnostics
!
      if (SW_diag%L_cloud_extinction) then
        do i=1, nclds
         do l=1, nlit
            SW_diag%cloud_extinction                                    &
     &         (col_list(l), row_list(l), i)                            &
     &         =cloud_extinction_g(l, n_layer+1-i)
             SW_diag%cloud_weight_extinction                            &
     &         (col_list(l), row_list(l), i)                            &
     &         =cloud_weight_extinction_g(l, n_layer+1-i)
         end do
        end do
      end if
!
      if (SW_diag%L_ls_cloud_extinction) then
        do i=1, nclds
         do l=1, nlit
            SW_diag%ls_cloud_extinction                                 &
     &         (col_list(l), row_list(l), i)                            &
     &         =ls_cloud_extinction_g(l, n_layer+1-i)
            SW_diag%ls_cloud_weight_extinction                          &
     &         (col_list(l), row_list(l), i)                            &
     &         =ls_cloud_weight_extinction_g(l, n_layer+1-i)
         end do
        end do
      end if
!
      if (SW_diag%L_cnv_cloud_extinction) then
        do i=1, nclds
         do l=1, nlit
            SW_diag%cnv_cloud_extinction                                &
     &         (col_list(l), row_list(l), i)                            &
     &         =cnv_cloud_extinction_g(l, n_layer+1-i)
            SW_diag%cnv_cloud_weight_extinction                         &
     &         (col_list(l), row_list(l), i)                            &
     &         =cnv_cloud_weight_extinction_g(l, n_layer+1-i)
         end do
        end do
      end if

! ######################################################
! COSP arguments
! ######################################################
      IF (l_cosp) THEN
        epsreal = EPSILON(1.0)
        DO i=1, nclds
          DO l=1, nlit
            j = (row_list(l)-1)*row_length + col_list(l)
!           Large-scale cloud optical depth
            IF (ls_cloud_weight_extinction_g(l, n_layer+1-i) > epsreal) THEN
              cosp_gbx%dtau_s(j,i) = d_mass(l,n_layer+1-i)                     &
                 * ls_cloud_extinction_g(l, n_layer+1-i)                       &
                 / ls_cloud_weight_extinction_g(l, n_layer+1-i)
            END IF
!           Convective cloud optical depth
            IF (cnv_cloud_weight_extinction_g(l, n_layer+1-i) > epsreal) THEN
              cosp_gbx%dtau_c(j,i) = d_mass(l,n_layer+1-i)                     &
                 * cnv_cloud_extinction_g(l, n_layer+1-i)                      &
                 / cnv_cloud_weight_extinction_g(l, n_layer+1-i)
              END IF
          END DO
        END DO
      END IF

!
!
!
!
!       Final processing of output fields
!
        if (L_orog) then
           do l = 1, nlit
              swout_temp(list(l))=swout(list(l),1)
              swout(list(l),1)=swout(list(l),1) - lg_f_orog(l)
           end do
           deallocate(lg_f_orog)
           deallocate(lg_orog_corr)
        end if

!       Convert the fluxes to increments.
        do i=nlevs, 1, -1
!
          if (l_scale_inc) then
            do l=1, nlit
              swout(list(l),i+1)=(swout(list(l),i+1)-swout(list(l),i))  &
     &          *pts/layer_heat_capacity(l, n_layer+1-i)
            end do
          else
            do l=1, nlit
              swout(list(l),i+1)=(swout(list(l),i+1)-swout(list(l),i))  &
     &          /layer_heat_capacity(l, n_layer+1-i)
            end do
          end if
!
          if (SW_diag%l_clear_hr) then
            do l=1, nlit
              SW_diag%clear_hr(col_list(l), row_list(l), i)             &
     &          =(flux_net_clear(l, n_layer-i, 1)                       &
     &          -flux_net_clear(l, n_layer+1-i, 1))                     &
     &          /layer_heat_capacity(l, n_layer+1-i)
            end do
          end if
!
        end do
!
        if (L_orog) then
           do l = 1, nlit
             swout(list(l),1)=swout_temp(list(l))
           end do
        end if
!
!
!
!       Separate the contributions over open sea and sea ice.
!       Fluxes returned from the radiation code itself are not
!       weighted by the fraction of the tile, but here are converted
!       to grid-box mean values. This can only be done if the surface
!       has been tiled.
!
        if (l_rad_tile) then
!
!         The variable flandg is set even if coastal tiling is not
!         used, so fairly generic code can be written. Note that
!         swsea is initialized in the calling routine, so its
!         setting at land points is implicit.
!
          do ll=1, nlit
            l=list(ll)
            if ( (flandg(l) < TINY(flandg)) .AND.                       &
     &           (ice_fraction(l) < TINY(ice_fraction) )                &
     &         ) then
!             This point is open sea with no sea ice.
              swsea(l)=swout(l, 1)
              swout(l, 1)=0.0
            end if
          end do
!
!         Tiled points will have both land and sea. Note that the
!         channel index of flux_up_tile is hard-wired to 1 because
!         we don't envisage calling the code in other cases.
          do lll=1, n_point_tile
            ll=list_tile(lll)
            l=list_tile_outer(lll)
            swsea(l)=(1.0-ice_fraction(l))*(flux_net(ll, n_layer, 1)    &
     &        +flux_up(ll, n_layer, 1)                                  &
     &        -flux_up_tile(lll, index_tile(ip_ocean_tile), 1))
             swout(l, 1)=swout(l, 1)-(1.0-flandg(l))*swsea(l)
          end do
!
!         The remaining points are entirely land points and swout
!         need not be altered.
!
        else
!
!         Without radiative tiling we must assume that fluxes are
!         uniform across the grid-box.
          Where (flandg(list(1:nlit)) < 1.0-TINY(flandg))
            swsea(list(1:nlit))=(1.0-ice_fraction(list(1:nlit)))        &
     &        *flux_net(1:nlit, n_layer, 1)
            swout(list(1:nlit), 1)=swout(list(1:nlit), 1)               &
     &        -(1.0-flandg(list(1:nlit)))*swsea(list(1:nlit))
          Endwhere
!
        end if
!
!
!       Only output from the main call should be scaled back by the
!       zenith angle.
        if (l_scale_inc) then
!         Divide by cosine of solar zenith angle to provide values for
!         upper routines. This applies only to SWOUT.
          do i=1, nlevs+1
            do l=1, n_points
              swout(l, i)=swout(l, i)/(coszin(l)*lit(l)+tiny(1.0))
            end do
          end do
        end if
!
!       Process surface fluxes
          do i=1, 4
            surf_down_sw(1:n_points, i) = surf_down_sw(1:n_points, i)   &
     &        / (coszin(1:n_points) * lit(1:n_points) + TINY(coszin) )
          end do
!
      else
!
!       Radiances are being calculated. This code is written on the
!       assumption that only one radiance is calculated at a point:
!       this can of course be generalized.
!
        if (SW_diag%l_toa_radiance) then
          do ic=1, n_channel
            do l=1, nlit
              SW_diag%toa_radiance(col_list(l), row_list(l), ic)        &
     &          =radiance(l, 1, 1, ic)
            end do
          end do
        end if
!
      end if
!
!
      IF (ALLOCATED(cloud_inhom_param)) DEALLOCATE(cloud_inhom_param)
      IF (sw_control%i_cloud==IP_cloud_mcica) THEN
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

      IF (lhook) CALL dr_hook('R2_SWRAD3Z',zhook_out,zhook_handle)
      END SUBROUTINE r2_swrad3z
