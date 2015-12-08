! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the radiance field.
!
! Method:
!   Properties independent of the spectral bands are set.
!   a loop over bands is then entered. Grey optical properties
!   are set and an appropriate subroutine is called to treat
!   the gaseous overlaps. The final radiances are assigned.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiance Core
!
!- ---------------------------------------------------------------------
SUBROUTINE radiance_calc(ierr                                           &
!               Logical flags for processes
    , l_rayleigh, l_aerosol, l_gas, l_continuum                         &
    , l_cloud, l_drop, l_ice                                            &
!               Angular integration
    , i_angular_integration, l_rescale, n_order_forward, i_2stream      &
    , n_order_gauss, i_truncation, ls_global_trunc, ms_min, ms_max      &
    , accuracy_adaptive, euler_factor                                   &
    , l_henyey_greenstein_pf, l_lanczos, ls_brdf_trunc                  &
    , i_sph_algorithm, n_order_phase_solar                              &
    , n_direction, direction                                            &
    , n_viewing_level, viewing_level, i_sph_mode                        &
!               Treatment of scattering
    , i_scatter_method                                                  &
!               Options for treating clouds
    , l_global_cloud_top, n_cloud_top_global                            &
    , l_inhom_cloud, inhom_cloud                                        &
!               Options for solver
    , i_solver                                                          &
!               Properties of diagnostics
    , map_channel                                                       &
!               General spectral properties
    , n_band, i_first_band, i_last_band, weight_band                    &
    , l_exclude, n_band_exclude, index_exclude                          &
!               General atmospheric properties
    , n_profile, n_layer                                                &
    , p, t, t_ground, t_level, d_mass                                   &
!               Spectral region
    , isolir                                                            &
!               Solar fields
    , zen_0, solar_irrad, solar_flux_band, solar_flux_band_ses          &
    , l_solar_tail_flux, solar_tail_flux, rayleigh_coefficient          &
!               Infra-red fields
    , n_deg_fit, thermal_coefficient, t_ref_planck                      &
    , l_ir_source_quad                                                  &
!               Gaseous absorption
    , i_gas_overlap, i_gas                                              &
    , gas_mix_ratio, n_band_absorb, index_absorb                        &
    , i_band_esft, w_esft, k_esft, i_scale_esft                         &
    , i_scale_fnc, scale_vector                                         &
    , p_reference, t_reference, l_mod_k_flux                            &
    , mix_gas_band, n_mix_gas, index_mix_gas, f_mix                     &
    , i_band_esft_ses, k_esft_ses, k_mix_gas, w_esft_ses                &
!               Doppler broadening
    , l_doppler, doppler_correction                                     &
!               Surface fields
    , n_brdf_basis_fnc, rho_alb, f_brdf                                 &
!               Tiling options for heterogeneous surfaces
    , l_tile, n_point_tile, n_tile, list_tile, rho_alb_tile             &
    , frac_tile, t_tile                                                 &
!               Continuum absorption
    , n_band_continuum, index_continuum, index_water                    &
    , k_continuum, i_scale_fnc_cont, scale_continuum                    &
    , p_ref_continuum, t_ref_continuum, k_continuum_ses, k_h2oc         &
!               Properties of aerosols
    , n_aerosol, n_aerosol_mr, aerosol_mix_ratio                        &
    , aerosol_mr_source, aerosol_mr_type_index, l_use_arcl              &
    , aerosol_absorption, aerosol_scattering                            &
    , n_aerosol_phf_term, aerosol_phase_fnc                             &
    , i_aerosol_parametrization, nhumidity, humidities                  &
    , type_aerosol                                                      &
!               Aerosol optical depth
    , n_aod_wavel, aod_wavel                                            &
    , l_aod_sulphate, aod_sulphate, l_aod_dust,      aod_dust           &
    , l_aod_seasalt,  aod_seasalt,  l_aod_soot,      aod_soot           &
    , l_aod_biomass,  aod_biomass,  l_aod_biogenic,  aod_biogenic       &
    , l_aod_ocff,     aod_ocff,     l_aod_delta,     aod_delta          &
    , l_aod_nitrate,  aod_nitrate,  l_aod_total_radn,aod_total_radn     &
    , l_angst_total_radn,   angst_total_radn                            &
    , l_aod_prog_sulphate,  aod_prog_sulphate,   l_aod_prog_dust        &
    , aod_prog_dust,        l_aod_prog_seasalt,  aod_prog_seasalt       &
    , l_aod_prog_soot,      aod_prog_soot,       l_aod_prog_biomass     &
    , aod_prog_biomass,     l_aod_prog_ocff,     aod_prog_ocff          &
    , l_aod_prog_nitrate,   aod_prog_nitrate                            &
    , aod_absorption, aod_scattering, i_aod_type                        &
!               Properties of UKCA aerosols
    , l_ukca_radaer, ukca_radaer, ukca_mix_ratio, ukca_comp_vol         &
    , ukca_dry_diam, ukca_wet_diam, ukca_modal_rho, ukca_modal_vol      &
    , ukca_modal_wtv, ukca_modal_nbr                                    &
!               Optical depth of UKCA aerosols
    , l_aod_ukca_ait_sol, aod_ukca_ait_sol                              &
    , l_aod_ukca_acc_sol, aod_ukca_acc_sol                              &
    , l_aod_ukca_cor_sol, aod_ukca_cor_sol                              &
    , l_aod_ukca_ait_ins, aod_ukca_ait_ins                              &
    , l_aod_ukca_acc_ins, aod_ukca_acc_ins                              &
    , l_aod_ukca_cor_ins, aod_ukca_cor_ins                              &
!               Properties of clouds
    , n_condensed, type_condensed                                       &
    , i_cloud, i_cloud_representation, w_cloud                          &
    , n_cloud_type, frac_cloud, tot_cloud_cover                         &
    , condensed_mix_ratio, condensed_dim_char                           &
    , i_condensed_param, condensed_n_phf, condensed_param_list          &
    , dp_corr_strat, dp_corr_conv                                       &
!               Calculated Fluxes or Radiances
    , flux_direct, flux_down, flux_up                                   &
    , uv_flux_direct, uv_flux_down, uv_flux_up                          &
    , l_uvflux_direct, l_uvflux_down, l_uvflux_up                       &
    , radiance, photolysis                                              &
!               Options for clear-sky fluxes
    , l_clear, i_solver_clear                                           &
!               Clear-sky fluxes calculated
    , flux_direct_clear, flux_down_clear, flux_up_clear                 &
!               Special Surface Fluxes
    , l_blue_flux_surf, weight_blue,weight_uv                           &
    , flux_direct_blue_surf, flux_down_blue_surf, flux_up_blue_surf     &
    , l_surf_uv, flux_down_uv_surf                                      &
    , l_surf_uv_clr, flux_down_clr_uv_surf                              &
!               Tiled Surface Fluxes
    , flux_up_tile, flux_up_blue_tile                                   &
!               Arrays for diagnostics specific to the UM
    , l_cloud_extinction, cloud_extinction                              &
    , cloud_weight_extinction                                           &
    , l_ls_cloud_extinction, ls_cloud_extinction                        &
    , ls_cloud_weight_extinction                                        &
    , l_cnv_cloud_extinction, cnv_cloud_extinction                      &
    , cnv_cloud_weight_extinction                                       &
    , l_cloud_absorptivity, cloud_absorptivity                          &
    , cloud_weight_absorptivity                                         &
    , l_ls_cloud_absorptivity, ls_cloud_absorptivity                    &
    , ls_cloud_weight_absorptivity                                      &
    , l_cnv_cloud_absorptivity, cnv_cloud_absorptivity                  &
    , cnv_cloud_weight_absorptivity                                     &
!               Dimensions of arrays
    , nd_profile, nd_layer, nd_column, nd_layer_clr, id_ct              &
    , nd_2sg_profile, nd_flux_profile, nd_radiance_profile              &
    , nd_j_profile, nd_channel, nd_band                                 &
    , nd_species, nd_esft_term, nd_scale_variable                       &
    , nd_continuum, nd_aerosol_species, nd_aerosol_mixratio             &
    , nd_humidities, nd_cloud_parameter, nd_thermal_coeff               &
    , nd_source_coeff, nd_brdf_basis_fnc, nd_brdf_trunc                 &
    , nd_aod_wavel, nd_phase_term, nd_max_order, nd_sph_coeff           &
    , nd_direction, nd_viewing_level                                    &
    , nd_region, nd_cloud_type, nd_cloud_component                      &
    , nd_overlap_coeff, nd_point_tile, nd_tile                          &
    , nd_tmp, nd_pre, nd_mix, nd_band_mix_gas                           &
    , n_ukca_mode, n_ukca_cpnt, nd_exclude                              &
    )
   USE arcl_mod, ONLY: npd_arcl_species


  USE realtype_rd, ONLY: RealK
  USE rad_pcf
  USE def_ss_prop
  USE ukca_radaer_struct_mod
  USE dust_parameters_mod, ONLY: l_twobin_dust
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE mcica_mod, ONLY: cloud_inhom_param
  USE ereport_mod, ONLY: ereport

  USE aodtype_mod, ONLY: ip_type_allaod, ip_type_sulphate, ip_type_dust,&
                         ip_type_seasalt, ip_type_soot, ip_type_biomass,&
                         ip_type_biogenic, ip_type_ocff, ip_type_delta, &
                         ip_type_nitrate, ip_type_twobdust
  IMPLICIT NONE


! Sizes of dummy arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for atmospheric profiles
    , nd_layer                                                          &
!       Size allocated for atmospheric layers
    , nd_layer_clr                                                      &
!       Size allocated for totally clear layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_2sg_profile                                                    &
!       Size allocated for profiles of fluxes
    , nd_flux_profile                                                   &
!       Size allocated for profiles of output fluxes
    , nd_radiance_profile                                               &
!       Size allocated for profiles of radiances
    , nd_j_profile                                                      &
!       Size allocated for profiles of mean radiances
    , nd_channel                                                        &
!       Size allocated for channels of output
    , nd_band                                                           &
!       Size allocated for bands in spectral computation
    , nd_species                                                        &
!       Size allocated for gaseous species
    , nd_continuum                                                      &
!       Size allocated for types of continua
    , nd_aerosol_species                                                &
!       Size allocated for aerosol species in spectral information
    , nd_aerosol_mixratio                                               &
!       Size allocated for aerosol species in mixing array
    , nd_humidities                                                     &
!       Size allocated for humidities
    , nd_esft_term                                                      &
!       Size allocated for ESFT terms
    , nd_scale_variable                                                 &
!       Size allocated for variables in scaling functions
    , nd_cloud_parameter                                                &
!       Size allocated for cloud parameters
    , nd_thermal_coeff                                                  &
!       Size allocated for thermal coefficients
    , nd_source_coeff                                                   &
!       Size allocated for two-stream source coefficients
    , nd_brdf_basis_fnc                                                 &
!       Size allowed for BRDF basis functions
    , nd_brdf_trunc                                                     &
!       Size allowed for truncation of BRDF basis functions
    , nd_column                                                         &
!       Size allocated for columns at each grid-point
    , nd_aod_wavel                                                      &
!       Number of wavelengths for the Aerosol Optical Depth
    , nd_phase_term                                                     &
!       Size allocated for terms in the phase function
!       supplied to the routine
    , nd_max_order                                                      &
!       Size allocated for polar orders
    , nd_direction                                                      &
!       Size allocated for viewing directions at each point
    , nd_viewing_level                                                  &
!       Size allocated for levels where the radiance
!       may be calculated
    , nd_region                                                         &
!       Size allocated for cloudy regions
    , nd_cloud_type                                                     &
!       Size allocated for types of clouds
    , nd_cloud_component                                                &
!       Size allocated for components in clouds
    , nd_overlap_coeff                                                  &
!       Size allocated for overlap coefficients
    , nd_sph_coeff                                                      &
!       Size allocated for arrays of spherical coefficients
!       used in determining radiances
    , nd_point_tile                                                     &
!       Size allocated for points with surface tiling
    , nd_tile                                                           &
!       Size allocated for the number of tiles
    , nd_tmp                                                            &
!       Number of reference temperatures for k-distribution
    , nd_pre                                                            &
!       Number of reference pressures for k-distribution
    , nd_mix                                                            &
!       Number of interpolation parameters for mixture species
    , nd_band_mix_gas                                                   &
!       Number of bands where mixed species exist
    , n_ukca_mode                                                       &
!       Number of aerosol modes in UKCA_RADAER
    , n_ukca_cpnt                                                       &
!       Number of aerosol components in UKCA_RADAER
    , nd_exclude
!       Size allocated for excluded bands

  LOGICAL, INTENT(IN) ::                                                &
      l_exclude
!       Presence of excluded bands
  INTEGER, INTENT(IN) ::                                                &
      n_band_exclude
!       Number of excluded bands
  INTEGER, INTENT(IN) ::                                                &
      index_exclude(nd_exclude, nd_band)
!       Index of excluded bands from each band


! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag

! General logical switches:
  LOGICAL, INTENT(IN) ::                                                &
      l_clear                                                           &
!       Calculate clear-sky fluxes
    , l_ir_source_quad                                                  &
!       Use a quadratic source function
    , l_rescale                                                         &
!       Flag for delta-rescaling
    , l_henyey_greenstein_pf                                            &
!       Use Henyey-Greenstein phase functions
    , l_lanczos
!       Flag to use Lanczos smoothing of solar phf

! Parameters controlling algorithms:
! Representation of clouds:
  INTEGER, INTENT(IN) ::                                                &
      i_cloud
!       Cloud scheme used
! Numerical algorithms:
  INTEGER, INTENT(IN) ::                                                &
      map_channel(nd_band)
!       Mapping of actual bands to the output channels
  INTEGER, INTENT(IN) ::                                                &
      isolir                                                            &
!       Visible or IR
    , i_solver                                                          &
!       Solver used
    , i_solver_clear                                                    &
!       Clear solver used
    , i_2stream                                                         &
!       Two-stream scheme
    , i_angular_integration                                             &
!       Angular integration scheme
    , n_order_gauss                                                     &
!       Order of Gaussian integration
    , i_truncation                                                      &
!       Type of spherical truncation
    , ls_global_trunc                                                   &
!       Truncating order of spherical harmonics
    , ms_min                                                            &
!       Lowest azimuthal order calculated
    , ms_max                                                            &
!       Highest azimuthal order calculated
    , n_order_forward                                                   &
!       Order of the term used to `define' the forward scattering
!       fraction.
    , i_sph_mode                                                        &
!       Mode in which the spherical harmonic solver is being used
    , i_sph_algorithm
!       Algorithm used for spherical harmonic calculation
  REAL (RealK), INTENT(IN) ::                                           &
      accuracy_adaptive                                                 &
!       Accuracy for adaptive truncation
    , euler_factor
!       Factor applied to the last term of an alternating series
  INTEGER, INTENT(INOUT) ::                                             &
      ls_brdf_trunc
!       Order of truncation applied to BRDFs
!       (This will be reset to 0 if a Lambertian surface is assumed)

! Specification of the viewing geometry
  INTEGER, INTENT(IN) ::                                                &
      n_direction                                                       &
!       Number of directions at which to calculate radiances
    , n_viewing_level
!       Number of levels where the radiance is required
  REAL (RealK), INTENT(IN) ::                                           &
      direction(nd_radiance_profile, nd_direction, 2)                   &
!       Directions in which to calculate radiances
    , viewing_level(nd_viewing_level)
!       List of levels where the radiance is required

! Range of spectral bands:
  INTEGER, INTENT(IN) ::                                                &
      i_first_band                                                      &
!       First band
    , i_last_band
!       Last band

! General properties of spectrum and aerosol data:
  INTEGER, INTENT(IN) ::                                                &
      n_band                                                            &
!       Number of spectral bands
    , n_aerosol                                                         &
!       Number of aerosol species in spectral information
    , n_aerosol_mr
!       Number of aerosol species in aerosol_mix_ratio array

! Solar fields:
  REAL (RealK), INTENT(IN) ::                                           &
      solar_irrad(nd_profile)                                           &
!       Incident solar radiation
    , solar_flux_band(nd_band)                                          &
!       Normalized flux in each spectral band
    , solar_flux_band_ses(nd_esft_term, nd_band)                        &
!       Normalised flux for each k-term
    , zen_0(nd_profile)
!       Secant (two-stream) or cosine (spherical harmonics)
!       of solar zenith angle
  LOGICAL, INTENT(IN) ::                                                &
      l_solar_tail_flux
!       Flag for adding solar tail flux to LW region
  REAL (RealK), INTENT(OUT) ::                                          &
      solar_tail_flux(nd_profile)
!       Solar tail flux considered in LW region

! Atmospheric profiles:
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer
!       Number of layers
  REAL (RealK), INTENT(IN) ::                                           &
      p(nd_profile, nd_layer)                                           &
!       Pressure
    , t(nd_profile, nd_layer)                                           &
!       Temperature
    , t_ground(nd_profile)                                              &
!       Temperature of ground
    , t_level(nd_profile, 0: nd_layer)                                  &
!       Temperature on levels
    , d_mass(nd_profile, nd_layer)                                      &
!       Mass thickness of each layer
    , gas_mix_ratio(nd_profile, nd_layer, nd_species)
!       Gaseous mass mixing ratios

! Surface properties:
  INTEGER, INTENT(IN) ::                                                &
      n_brdf_basis_fnc
!       Number of BRDF basis functions
  REAL (RealK), INTENT(INOUT) ::                                        &
      rho_alb(nd_profile, nd_brdf_basis_fnc, nd_band)                   &
!       Weights of the basis functions
    , f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                      &
        , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)
!       Array of BRDF basis terms

! Arrays related to tiling of the surface
  LOGICAL, INTENT(IN) ::                                                &
      l_tile
!       Local to allow tiling options
  INTEGER, INTENT(IN) ::                                                &
      n_point_tile                                                      &
!       Number of points to tile
    , n_tile                                                            &
!       Number of tiles used
    , list_tile(nd_point_tile)
!       List of points with surface tiling
  REAL (RealK), INTENT(INOUT) ::                                        &
      rho_alb_tile(nd_point_tile, nd_brdf_basis_fnc                     &
        , nd_tile, nd_band)                                             &
!       Weights for the basis functions of the BRDFs
!       at the tiled points
    , frac_tile(nd_point_tile, nd_tile)                                 &
!       Fraction of tiled grid-points occupied by each tile
    , t_tile(nd_point_tile, nd_tile)
!       Local surface temperatures on individual tiles

! Rayleigh scattering:
  LOGICAL, INTENT(IN) ::                                                &
      l_rayleigh
!       Include rayleigh scattering in the calculation.
  REAL (RealK), INTENT(IN) ::                                           &
       rayleigh_coefficient(nd_band)
!       Rayleigh coefficients

! fields for gaseous absorption:
  LOGICAL, INTENT(IN) ::                                                &
      l_gas
!       Include gas absorption in the calculation
! gaseous overlaps:
  INTEGER, INTENT(IN) ::                                                &
      i_gas_overlap(nd_band)                                            &
!       Gas overlap assumption
    , i_gas
!       Gas to be considered (one gas only)
! ESFTs:
  INTEGER, INTENT(IN) ::                                                &
      n_band_absorb(nd_band)                                            &
!       Number of absorbers in band
    , index_absorb(nd_species, nd_band)                                 &
!       List of absorbers in bands
    , i_band_esft(nd_band, nd_species)                                  &
!       Number of terms in band
    , i_scale_esft(nd_band, nd_species)                                 &
!       Type of esft scaling
    , i_scale_fnc(nd_band, nd_species)
!       Type of scaling function
  REAL (RealK), INTENT(IN) ::                                           &
      w_esft(nd_esft_term, nd_band, nd_species)                         &
!       Weights for ESFT
    , k_esft(nd_esft_term, nd_band, nd_species)                         &
!       Exponential ESFT terms
    , scale_vector(nd_scale_variable, nd_esft_term, nd_band             &
          , nd_species)                                                 &
!       Absorber scaling parameters
    , p_reference(nd_species, nd_band)                                  &
!       Reference scaling pressure
    , t_reference(nd_species, nd_band)
!       Reference scaling temperature

! Use modulus of fluxes to remove negative effective extinctions
  LOGICAL, INTENT(IN) :: l_mod_k_flux

  INTEGER, INTENT(IN) ::                                                &
      mix_gas_band(nd_band)                                             &
!       Sequence band number (not real band number) of mixed species
    , n_mix_gas(nd_band)                                                &
!       Index of band where mixed species occurs
    , index_mix_gas(2, nd_band_mix_gas)                                 &
!       Index of mixed species
    , i_band_esft_ses(nd_band)
!       Number of terms in band

  REAL (RealK), INTENT(IN) ::                                           &
      k_esft_ses(nd_pre, nd_tmp, nd_species, nd_esft_term, nd_band)     &
!       Absorption coefficients
    , w_esft_ses(nd_esft_term, nd_band)                                 &
!       k-term weights
    , k_mix_gas(nd_pre, nd_tmp, nd_mix, nd_esft_term, nd_band_mix_gas)  &
    , f_mix(nd_band)

! Spectral data for the continuum:
  LOGICAL, INTENT(IN) ::                                                &
      l_continuum
!       Include continuum absorption in the calculation
  INTEGER, INTENT(IN) ::                                                &
      n_band_continuum(nd_band)                                         &
!       Number of continua in bands
    , index_continuum(nd_band, nd_continuum)                            &
!       Indices of continua
    , index_water                                                       &
!       Index of water
    , i_scale_fnc_cont(nd_band, nd_continuum)
!       Type of scaling function for continuum
  REAL (RealK), INTENT(IN) ::                                           &
      k_continuum(nd_band, nd_continuum)                                &
!       Continuum extinction coefficients
    , k_continuum_ses(nd_esft_term, nd_tmp, nd_band, nd_continuum)      &
!       Continuum extinction coefficients (for SES2)
    , k_h2oc(nd_pre, nd_tmp, nd_esft_term, nd_band)                     &
!       H2O continuum Absorption coefficients
    , scale_continuum(nd_scale_variable, nd_band, nd_continuum)         &
!       Continuum scaling parameters
    , p_ref_continuum(nd_continuum, nd_band)                            &
!       Continuum reference pressure
    , t_ref_continuum(nd_continuum, nd_band)
!       Continuum reference temperature


! General cloud fields:
  LOGICAL, INTENT(IN) ::                                                &
      l_cloud
!       Clouds are required in the calculation
  REAL (RealK), INTENT(IN) ::                                           &
      w_cloud(nd_profile, id_ct: nd_layer)                              &
!       Amount of cloud
    , frac_cloud(nd_profile, id_ct: nd_layer, nd_cloud_type)            &
!       Fractions of different types of cloud
    , dp_corr_strat                                                     &
!       Decorrelation pressure scale for large scale cloud
    , dp_corr_conv
!       Decorrelation pressure scale for convective cloud
  REAL (RealK), INTENT(INOUT) ::                                        &
      tot_cloud_cover(nd_profile)
!       Total cloud cover


! Fields for microphysical quantities:
  LOGICAL, INTENT(IN) ::                                                &
      l_drop                                                            &
!       Include droplets in the calculation
    , l_ice                                                             &
!       Include ice in the calculation
    , l_global_cloud_top                                                &
!       Use a global value for the top of clouds
!       (This is for use in a GCM where the code must be
!       bit-reproducible across different configurations of PEs).
    , l_inhom_cloud
!       Use scaling factors for inhomogeneous cloud
  INTEGER, INTENT(IN) ::                                                &
      n_condensed                                                       &
!       Number of condensed components in clouds
    , type_condensed(nd_cloud_component)                                &
!       Types of condensed components
    , i_condensed_param(nd_cloud_component)                             &
!       Parametrization schemes for components
    , condensed_n_phf(nd_cloud_component)                               &
!       Number of terms in the phase function
    , i_cloud_representation                                            &
!       Representation of mixing rule chosen
    , n_cloud_type                                                      &
!       Number of types of cloud
    , n_cloud_top_global
!       Global cloud top

  REAL (RealK), INTENT(INOUT) ::                                        &
      condensed_mix_ratio(nd_profile, id_ct: nd_layer                   &
        , nd_cloud_component)
!       Mixing ratios of condensed components

  REAL (RealK), INTENT(IN) ::                                           &
      condensed_dim_char(nd_profile, id_ct: nd_layer                    &
        , nd_cloud_component)                                           &
!       Characteristic dimensions of condensed components
    , condensed_param_list(nd_cloud_parameter                           &
        , nd_cloud_component, nd_band)                                  &
!       Coefficients in parametrizations of condensed phases
    , inhom_cloud(nd_cloud_component)
!       Scaling factors for inhomogeneous cloud



! Fields for aerosols:
  LOGICAL, INTENT(IN) ::                                                &
      l_aerosol                                                         &
!       Include aerosols in the calculation
    , l_use_arcl(npd_arcl_species)
!       Logical for each aerosol climatology species
  INTEGER, INTENT(IN) ::                                                &
      i_aerosol_parametrization(nd_aerosol_species)                     &
!       Parametrization flags for aerosol
    , n_aerosol_phf_term(nd_aerosol_species)
!       Number of terms in the phase function of aerosols
  INTEGER, INTENT(IN) ::                                                &
      nhumidity(nd_aerosol_species)                                     &
!       Number of humidities
    , aerosol_mr_type_index(nd_aerosol_mixratio)                        &
!       Index relating aerosol_mix_ratio aerosols to aerosols in
!       the spectral information
    , aerosol_mr_source(nd_aerosol_mixratio)
!       Scheme/source of the aerosol data, to determine use in
!       changing radiative fluxes and use in diagnostics

  REAL (RealK), INTENT(IN) ::                                           &
      aerosol_mix_ratio(nd_profile, nd_layer, nd_aerosol_mixratio)
!       Mixing ratio of aerosols
  REAL (RealK), INTENT(IN) ::                                           &
      aerosol_absorption(nd_humidities, nd_aerosol_species, nd_band)    &
!       Absorption by aerosols
    , aerosol_scattering(nd_humidities, nd_aerosol_species, nd_band)    &
!       Scattering by aerosols
    , aerosol_phase_fnc(nd_humidities, nd_phase_term                    &
          , nd_aerosol_species, nd_band)                                &
!       Phase function of aerosols
    , humidities(nd_humidities, nd_aerosol_species)
!       Humidities for species



! Fields for Aerosol Optical Depth

  LOGICAL, INTENT(IN) ::                                                &
      l_aod_sulphate,                                                   &
      l_aod_dust,                                                       &
      l_aod_seasalt,                                                    &
      l_aod_soot,                                                       &
      l_aod_biomass,                                                    &
      l_aod_biogenic,                                                   &
      l_aod_ocff,                                                       &
      l_aod_delta,                                                      &
      l_aod_nitrate,                                                    &
      l_aod_total_radn,                                                 &
      l_angst_total_radn,                                               &
      l_aod_prog_sulphate,                                              &
      l_aod_prog_dust,                                                  &
      l_aod_prog_seasalt,                                               &
      l_aod_prog_soot,                                                  &
      l_aod_prog_biomass,                                               &
      l_aod_prog_ocff,                                                  &
      l_aod_prog_nitrate

  INTEGER, INTENT(IN) ::                                                &
      n_aod_wavel
  REAL (RealK), INTENT(IN) ::                                           &
      aod_wavel(nd_aod_wavel),                                          &
      aod_absorption(nd_humidities, nd_aerosol_species, nd_aod_wavel),  &
      aod_scattering(nd_humidities, nd_aerosol_species, nd_aod_wavel)
  INTEGER, INTENT(IN) ::                                                &
      i_aod_type(nd_aerosol_species)
  INTEGER, INTENT(IN) ::                                                &
      type_aerosol(nd_aerosol_species)
  REAL (RealK), INTENT(OUT) ::                                          &
      aod_sulphate(nd_profile, nd_aod_wavel),                           &
      aod_dust(nd_profile, nd_aod_wavel),                               &
      aod_seasalt(nd_profile, nd_aod_wavel),                            &
      aod_soot(nd_profile, nd_aod_wavel),                               &
      aod_biomass(nd_profile, nd_aod_wavel),                            &
      aod_biogenic(nd_profile, nd_aod_wavel),                           &
      aod_ocff(nd_profile, nd_aod_wavel),                               &
      aod_delta(nd_profile, nd_aod_wavel),                              &
      aod_nitrate(nd_profile, nd_aod_wavel),                            &
      aod_total_radn(nd_profile, nd_aod_wavel),                         &
      angst_total_radn(nd_profile, nd_aod_wavel),                       &
      aod_prog_sulphate(nd_profile, nd_aod_wavel),                      &
      aod_prog_dust(nd_profile, nd_aod_wavel),                          &
      aod_prog_seasalt(nd_profile, nd_aod_wavel),                       &
      aod_prog_soot(nd_profile, nd_aod_wavel),                          &
      aod_prog_biomass(nd_profile, nd_aod_wavel),                       &
      aod_prog_ocff(nd_profile, nd_aod_wavel),                          &
      aod_prog_nitrate(nd_profile, nd_aod_wavel)

! Fields for UKCA aerosols
  LOGICAL, INTENT(IN) ::                                                &
      l_ukca_radaer
!       Model switch
  TYPE (ukca_radaer_struct), INTENT(IN) :: ukca_radaer
!       UKCA_RADAER structure
  REAL (RealK), INTENT(IN) ::                                           &
      ukca_mix_ratio(nd_profile, nd_layer, n_ukca_cpnt),                &
      ukca_comp_vol(nd_profile, nd_layer, n_ukca_cpnt),                 &
      ukca_dry_diam(nd_profile, nd_layer, n_ukca_mode),                 &
      ukca_wet_diam(nd_profile, nd_layer, n_ukca_mode),                 &
      ukca_modal_rho(nd_profile, nd_layer, n_ukca_mode),                &
      ukca_modal_vol(nd_profile, nd_layer, n_ukca_mode),                &
      ukca_modal_wtv(nd_profile, nd_layer, n_ukca_mode),                &
      ukca_modal_nbr(nd_profile, nd_layer, n_ukca_mode)
!       UKCA aerosol component mass-mixing ratios, modal dry and
!       wet diameters, modal densities, volumes, and volumes of
!       water.
  LOGICAL, INTENT(IN) ::                                                &
      l_aod_ukca_ait_sol,                                               &
      l_aod_ukca_acc_sol,                                               &
      l_aod_ukca_cor_sol,                                               &
      l_aod_ukca_ait_ins,                                               &
      l_aod_ukca_acc_ins,                                               &
      l_aod_ukca_cor_ins
!       Switch for UKCA modal optical depth diagnostics.
  REAL (RealK), INTENT(OUT) ::                                          &
      aod_ukca_ait_sol(nd_profile, nd_aod_wavel),                       &
      aod_ukca_acc_sol(nd_profile, nd_aod_wavel),                       &
      aod_ukca_cor_sol(nd_profile, nd_aod_wavel),                       &
      aod_ukca_ait_ins(nd_profile, nd_aod_wavel),                       &
      aod_ukca_acc_ins(nd_profile, nd_aod_wavel),                       &
      aod_ukca_cor_ins(nd_profile, nd_aod_wavel)
!       UKCA modal optical depth diagnostics.

! Fitting of the Planckian function:
  INTEGER, INTENT(IN) ::                                                &
      n_deg_fit
!       Degree of thermal fitting fnc.
  REAL (RealK), INTENT(IN) ::                                           &
      thermal_coefficient(0: nd_thermal_coeff-1, nd_band)               &
!       Coefficients of source terms
    , t_ref_planck
!       Planckian reference temperature

! Doppler broadening
  LOGICAL, INTENT(IN) ::                                                &
      l_doppler(nd_species)
!       Flags to activate doppler corrections
  REAL (RealK), INTENT(IN) ::                                           &
      doppler_correction(nd_species)
!       Doppler broadening term
  REAL (RealK), INTENT(IN) ::                                           &
      weight_band(nd_band)
!       Weighting function for bands

! Control of scattering:
  INTEGER, INTENT(IN) ::                                                &
      i_scatter_method(nd_band)
!       Method of treating scattering in each band

! Fluxes or radiances calculated:
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_flux_profile, 0: nd_layer, nd_channel)             &
!       Direct flux
    , flux_down(nd_flux_profile, 0: nd_layer, nd_channel)               &
!       Downward flux
    , flux_up(nd_flux_profile, 0: nd_layer, nd_channel)                 &
!       Upward flux
    , flux_direct_clear(nd_2sg_profile, 0: nd_layer, nd_channel)        &
!       Clear direct flux
    , flux_down_clear(nd_2sg_profile, 0: nd_layer, nd_channel)          &
!       Clear downward flux
    , flux_up_clear(nd_2sg_profile, 0: nd_layer, nd_channel)            &
!       Clear upward flux
    , uv_flux_direct(nd_flux_profile, 0: nd_layer, nd_channel)          &
!       Direct UV-flux
    , uv_flux_up(nd_flux_profile, 0: nd_layer, nd_channel)              &
!       Upwards UV-flux
    , uv_flux_down(nd_flux_profile, 0: nd_layer, nd_channel)            &
!       Downward UV-Flux
    , radiance(nd_radiance_profile, nd_viewing_level                    &
        , nd_direction, nd_channel)                                     &
!       Calculated radiances
    , photolysis(nd_j_profile, nd_viewing_level, nd_channel)
!       Rate of photolysis
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_up_tile(nd_point_tile, nd_tile, nd_channel)                  &
!       Upward fluxes at tiled surface points
    , flux_up_blue_tile(nd_point_tile, nd_tile, nd_channel)
!       Upward blue fluxes at tiled surface points

! Special Diagnostics:
  LOGICAL, INTENT(IN) ::                                                &
      l_blue_flux_surf                                                  &
!       Flag to calculate the blue surface fluxes
     , l_uvflux_direct                                                  &
!       Flag for direct ultraviolet fluxes
     , l_uvflux_up                                                      &
!       Flag for upward ultraviolet fluxes
     , l_uvflux_down                                                    &
!       Flag for downward ultraviolet fluxes
     , l_surf_uv, l_surf_uv_clr
!       Flags for downward surface UV fluxes
  REAL (RealK), INTENT(IN) ::                                           &
      weight_blue(nd_band)                                              &
!       Weights for each band for blue fluxes
    , weight_uv(nd_band)
!       Weights for each band for the UV-interval
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct_blue_surf(nd_flux_profile)                            &
!       Direct blue flux at the surface
    , flux_down_blue_surf(nd_flux_profile)                              &
!       Total downward blue flux at the surface
    , flux_up_blue_surf(nd_flux_profile)                                &
!       Upward blue flux at the surface
    , flux_down_uv_surf(nd_flux_profile)                                &
!       Total downward UV flux at the surface
    , flux_down_clr_uv_surf(nd_flux_profile)
!       Clear-sky downward UV flux at the surface

! Arrays specific to the Unified Model

! Switches for diagnostics:
  LOGICAL, INTENT(IN) ::                                                &
      l_cloud_absorptivity,                                             &
!       FLAG TO CALCULATE ABSORPTIVITY OF CLOUDS
!       (ONLY INFRA-RED)
      l_cloud_extinction,                                               &
!       FLAG TO CALCULATE EXTINCTION OF CLOUDS
!       (ONLY SOLAR)
      l_ls_cloud_absorptivity,                                          &
!       FLAG TO CALCULATE ABSORPTIVITY OF LAYER CLOUDS
!       (ONLY INFRA-RED)
      l_ls_cloud_extinction,                                            &
!       FLAG TO CALCULATE EXTINCTION OF LAYER CLOUDS
!       (ONLY SOLAR)
      l_cnv_cloud_absorptivity,                                         &
!       FLAG TO CALCULATE ABSORPTIVITY OF CONV.CLOUDS
!       (ONLY INFRA-RED)
      l_cnv_cloud_extinction
!       FLAG TO CALCULATE EXTINCTION OF CONV.CLOUDS
!       (ONLY SOLAR)

! Diagnostics for clouds
  REAL (RealK), INTENT(OUT) ::                                          &
      cloud_absorptivity(nd_profile, nd_layer),                         &
!       Absorptivity of cloud weighted by cloud fraction
!       and upward clear-sky infra-red flux.
      cloud_weight_absorptivity(nd_profile, nd_layer),                  &
!       Weights to be applied to absorptivies.
      ls_cloud_absorptivity(nd_profile, nd_layer),                      &
!       Absorptivity of layer cloud weighted by cloud fraction
!       and upward clear-sky infra-red flux.
      ls_cloud_weight_absorptivity(nd_profile, nd_layer),               &
!       Weights to be applied to layer cld. absorptivies.
      cnv_cloud_absorptivity(nd_profile, nd_layer),                     &
!       Absorptivity of conv.cloud weighted by cloud fraction
!       and upward clear-sky infra-red flux.
      cnv_cloud_weight_absorptivity(nd_profile, nd_layer),              &
!       Weights to be applied to conv.cld absorptivies.
      cloud_extinction(nd_profile, nd_layer),                           &
!       Absorptivity of cloud weighted by cloud fraction
!       and downward clear-sky solar flux.
      cloud_weight_extinction(nd_profile, nd_layer),                    &
!       Weights to be applied to extinctions.
      ls_cloud_extinction(nd_profile, nd_layer),                        &
!       Absorptivity of layer cloud weighted by cloud fraction
!       and downward clear-sky solar flux.
      ls_cloud_weight_extinction(nd_profile, nd_layer),                 &
!       Weights to be applied to layer cld. extinctions.
      cnv_cloud_extinction(nd_profile, nd_layer),                       &
!       Absorptivity of conv.cloud weighted by cloud fraction
!       and downward clear-sky solar flux.
      cnv_cloud_weight_extinction(nd_profile, nd_layer)
!       Weights to be applied to conv. cld. extinctions.



! Local arguments.
! General pointers:
  INTEGER                                                               &
      i_top                                                             &
!       Top level of profiles
    , i_band                                                            &
!       Spectral band
    , n_gas                                                             &
!       Number of active gases
    , i_gas_band                                                        &
!       Single variable for gas in band
    , n_continuum                                                       &
!       Number of continua in band
    , i_continuum                                                       &
!       Continuum number
    , i_continuum_pointer(nd_continuum)                                 &
!       Pointers to continua
    , i_pointer_water
!       Pointer to water vapour

! Additional variables for angular integration:
  LOGICAL                                                               &
      l_solar_phf                                                       &
!       Logical to specify a separate treatment of the singly
!       scattered solar beam
    , l_rescale_solar_phf
!       Logical to apply rescaling to the singly scattered
!       solar phase function
  INTEGER                                                               &
      n_order_phase                                                     &
!       Order of Legendre polynomials of the phase function
    , n_order_phase_solar
!       Order of Legendre polynomials of the phase function
!       retained for the singly scattered solar beam

! Pointers to the contents of layers:
  INTEGER                                                               &
      n_cloud_top                                                       &
!       Topmost cloudy layer
    , n_region                                                          &
!       Number of cloudy regions
    , n_cloud_profile(id_ct: nd_layer)                                  &
!       Number of cloudy profiles
    , i_cloud_profile(nd_profile, id_ct: nd_layer)
!       Profiles containing clouds

! Pointers to types of clouds:
  LOGICAL                                                               &
      l_cloud_cmp(nd_cloud_component)
!       Logical switches to `include' components
  INTEGER                                                               &
      i_phase_cmp(nd_cloud_component)                                   &
!       Phases of components
    , i_cloud_type(nd_cloud_component)                                  &
!       Pypes of cloud to which each component contributes
    , type_region(nd_region)                                            &
!       The types of the regions
    , k_clr                                                             &
!       Index of clear-sky region
    , i_region_cloud(nd_cloud_type)
!       Regions in which particular type of cloud fall

! Fractional coverage of different regions:
  REAL (RealK) ::                                                       &
      frac_region(nd_profile, id_ct: nd_layer, nd_region)
!       Fraction of total cloud occupied by specific regions

! Pointer to table of humidities:
  INTEGER                                                               &
      i_humidity_pointer(nd_profile, nd_layer)
!       Pointer to look-up table for aerosols

! Controlling variables:
  INTEGER                                                               &
      i                                                                 &
!       Loop variable
    , j                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , l, ll
!       Loop variable

! Logical switches:
  LOGICAL                                                               &
      l_gas_band                                                        &
!       Flag to `include' gaseous absorption in a particular band
    , l_moist_aerosol                                                   &
!       Flag for moist aerosol
    , l_aerosol_density                                                 &
!       Flag for calculation of atmospheric density for aerosols
    , l_clear_band
!       Flag to calculate clear-sky fluxes for this band

  REAL (RealK) ::                                                       &
      solar_irrad_band(nd_profile)                                      &
!       Solar irradiance in the band
    , solar_irrad_band_ses(nd_profile, nd_esft_term)
!       Incident solar flux for each k-term
  REAL (RealK) ::                                                       &
      gas_frac_rescaled(nd_profile, nd_layer, nd_species)               &
!       Rescaled gas mixing ratios
    , gas_mix_amt(nd_profile, nd_layer)                                 &
!       Mixed gas mixing ratio
    , amount_continuum(nd_profile, nd_layer, nd_continuum)              &
!       Amounts of continua
    , k_continuum_mono(nd_continuum)
!       Monochromatic continuum components

! Thermal arrays:
  REAL (RealK) ::                                                       &
      planck_flux_band(nd_profile, 0: nd_layer)                         &
!       Planckian flux in band at edges of layers
    , diff_planck_band(nd_profile, nd_layer)                            &
!       Difference in the Planckian flux across layers
    , diff_planck_band_2(nd_profile, nd_layer)                          &
!       Twice the 2nd difference of in the Planckian flux across
!       layers
    , planck_flux_ground(nd_profile)
!       Planckian flux at the surface temperature

! Surface BRDF terms
  LOGICAL                                                               &
      l_diff_alb
!       Flag to calculate diffuse albedos
  REAL (RealK) ::                                                       &
      brdf_sol(nd_profile, nd_brdf_basis_fnc, nd_direction)             &
!       The BRDF evaluated for scattering from the solar
!       beam into the viewing direction
    , brdf_hemi(nd_profile, nd_brdf_basis_fnc, nd_direction)            &
!       The BRDF evaluated for scattering from isotropic
!       radiation into the viewing direction
    , diffuse_alb_basis(nd_brdf_basis_fnc)
!       The diffuse albedo of isotropic radiation for each
!       basis function

! Atmospheric densities:
  REAL (RealK) ::                                                       &
      density(nd_profile, nd_layer)                                     &
!       Overall density
    , molar_density_water(nd_profile, nd_layer)                         &
!       Molar density of water
    , molar_density_frn(nd_profile, nd_layer)
!       Molar density of foreign species

! Fields for moist aerosols:
  REAL (RealK) ::                                                       &
      delta_humidity                                                    &
!       Increment in look-up table for hum.
    , mean_rel_humidity(nd_profile, nd_layer)
!       Mean relative humidity of layers


! Fields for UKCA aerosols
  REAL                                                                  &
      ukca_modal_mixr(nd_profile, nd_layer, n_ukca_mode)
!         Total mass-mixing ratio for each mode
  REAL                                                                  &
      ukca_modal_number(nd_profile, nd_layer, n_ukca_mode)
!         Modal number concentrations
  REAL                                                                  &
      ukca_absorption(nd_profile, nd_layer, n_ukca_mode, nd_band)       &
!         Waveband-averaged UKCA aerosol absorption
    , ukca_scattering(nd_profile, nd_layer, n_ukca_mode, nd_band)       &
!         Waveband-averaged UKCA aerosol scattering
    , ukca_asymmetry(nd_profile, nd_layer,  n_ukca_mode, nd_band)
!         Waveband-averaged UKCA aerosol asymmetry

! Fundamental optical properties of layers:
  TYPE(str_ss_prop) :: ss_prop
!   Single scattering properties of the atmosphere

  REAL (RealK) ::                                                       &
      k_esft_layer(nd_profile, nd_esft_term, nd_layer, nd_species)      &
!       Exponential ESFT terms at actual pressure layer
    , k_mix_gas_layer(nd_profile, nd_esft_term, nd_layer)               &
!       Exponential ESFT terms at actual pressure layer
    , k_contm_layer(nd_profile, nd_esft_term, nd_layer, nd_continuum)
!       Continuum absorption coefficients at layer pressure

! Local variables for spherical harmonic integration
  INTEGER                                                               &
      ls_max_order                                                      &
!       Maximum order of terms required
    , ls_local_trunc(0: nd_max_order)                                   &
!       Actual truncation for each particular value of m
    , ms_trunc(0: nd_max_order)                                         &
!       Maximum azimuthal quantum number for each order
    , ia_sph_mm(0: nd_max_order)
!       Address of spherical coefficient of (m, m) for each m
  REAL (RealK) ::                                                       &
      cg_coeff(nd_sph_coeff)                                            &
!       Clebsch-gordon coefficients
    , uplm_zero(nd_sph_coeff)                                           &
!       Upsilon terms
    , uplm_sol(nd_radiance_profile, nd_sph_coeff)                       &
!       Upsilon terms for solar radiation
    , cos_sol_view(nd_radiance_profile, nd_direction)
!       Cosines of the angles between the solar direction
!       and the viewing direction

! Specification of the grid for radiances:
  INTEGER                                                               &
      i_rad_layer(nd_viewing_level)
!       Layers in which to intercept radiances
  REAL (RealK) ::                                                       &
      frac_rad_layer(nd_viewing_level)
!       Fractions below the tops of the layers

  REAL (RealK) ::                                                       &
      i_direct(nd_radiance_profile, 0: nd_layer)
!       Direct solar irradiance on levels (not split by
!       diagnostic bands or returned, but retained for
!       future use)
  REAL (RealK) ::                                                       &
      planck_radiance_band(nd_radiance_profile, nd_viewing_level)
!       Planckian radiance in the current band
  LOGICAL                                                               &
      l_initial
!       Flag to run the routine incrementing diagnostics in
!       its initializing mode

! Coefficients for the transfer of energy between
! Partially cloudy layers:
  REAL (RealK) ::                                                       &
      cloud_overlap(nd_profile, id_ct-1: nd_layer                       &
        , nd_overlap_coeff)                                             &
!       Coefficients defining overlapping options for clouds:
!       these also depend on the solver selected.
    , w_free(nd_profile, id_ct: nd_layer)
!       Clear-sky fraction

! Cloud geometry
  INTEGER                                                               &
      n_column_cld(nd_profile)                                          &
!       Number of columns in each profile (including those of
!       zero width)
    , n_column_slv(nd_profile)                                          &
!       Number of columns to be solved in each profile
    , list_column_slv(nd_profile, nd_column)                            &
!       List of columns requiring an actual solution
    , i_clm_lyr_chn(nd_profile, nd_column)                              &
!       Layer in the current column to change
    , i_clm_cld_typ(nd_profile, nd_column)
!       Type of cloud to introduce in the changed layer
  REAL (RealK) ::                                                       &
      area_column(nd_profile, nd_column)
!       Areas of columns

! Local variables for tiled fluxes:
  REAL (RealK) ::                                                       &
      planck_flux_tile(nd_point_tile, nd_tile)
!       Local Planckian fluxes on surface tiles

! Secondary arrays for diagnostics
  REAL (RealK) ::                                                       &
      cloud_absorptivity_band(nd_profile, nd_layer)                     &
!       Absorptivity of cloud in a particular band
    , cloud_extinction_band(nd_profile, nd_layer)                       &
!       Absorptivity of cloud in a particular band
    , ls_cloud_absorptivity_band(nd_profile, nd_layer)                  &
!       Absorptivity of cloud in a particular band
    , ls_cloud_extinction_band(nd_profile, nd_layer)                    &
!       Absorptivity of cloud in a particular band
    , cnv_cloud_absorptivity_band(nd_profile, nd_layer)                 &
!       Absorptivity of cloud in a particular band
    , cnv_cloud_extinction_band(nd_profile, nd_layer)                   &
!       Absorptivity of cloud in a particular band
    , flux_direct_clear_band(nd_2sg_profile, 0: nd_layer)               &
!       Clear direct flux in a particular band
    , flux_up_clear_band(nd_2sg_profile, 0: nd_layer)
!       Clear upward flux in a particular band

  INTEGER ::                                                            &
      jp(nd_profile, nd_layer)                                          &
!       Index for pressure interpolation of absorption coefficient
    , jph2oc(nd_profile, nd_layer)                                      &
!       Same as JP but for water vapour pressure
    , jt(nd_profile,nd_layer)                                           &
!       Index for temperature interpolation of absorption coeff
    , jtt(nd_profile, nd_layer)                                         &
!       Index of reference temperature at level i+1
!       such that the actual temperature is between JTT and JTT+1
    , jto2c(nd_profile, nd_layer)                                       &
!       Index of o2 continuum  reference temperature at level I
!       such that the actual temperature is between JTO2C and JTO2C+1
    , jtswo3(nd_profile, nd_layer)
!       Index of sw o3 reference temp

  REAL (RealK) ::                                                       &
      fac00(nd_profile, nd_layer), fac01(nd_profile, nd_layer)          &
    , fac10(nd_profile, nd_layer), fac11(nd_profile, nd_layer)          &
!       Multiplication factors for P & T interpolation
    , fac00c(nd_profile, nd_layer)                                      &
    , fac01c(nd_profile, nd_layer)                                      &
    , fac10c(nd_profile, nd_layer)                                      &
    , fac11c(nd_profile, nd_layer)                                      &
!       Multiplication factors for H2O cont P & T interpolation
    , facc00(nd_profile, nd_layer)                                      &
    , facc01(nd_profile, nd_layer)
!       Multiplication factors for O2 continuum T interpolation

  LOGICAL :: l_grey_cont
!       Flag to add continuum in grey_opt_prop


! Functions called:
  LOGICAL, EXTERNAL :: l_cloud_density
!       Flag for calculation of atmospheric densities for clouds


  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  CHARACTER (LEN=256)           :: cmessage
  CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'radiance_calc'


  IF (lhook) CALL dr_hook('RADIANCE_CALC',zhook_in,zhook_handle)

! Initial determination of flags and switches:

  IF (i_angular_integration == ip_two_stream) THEN

!   Only one term in the phase function is required.
    n_order_phase=1

    l_solar_phf=.FALSE.
    l_rescale_solar_phf=.FALSE.

  ELSE IF (i_angular_integration == ip_spherical_harmonic) THEN

!   Set limits on ranges of harmonics and set pointers to arrays.
! DEPENDS ON: set_truncation
    CALL set_truncation(ierr                                            &
      , i_truncation, ls_global_trunc                                   &
      , ls_max_order, ls_local_trunc                                    &
      , ms_min, ms_max, ms_trunc                                        &
      , ia_sph_mm, n_order_phase                                        &
      , nd_max_order                                                    &
      )

!   Determine whether special treatment of the solar
!   beam is required.
    l_solar_phf=(isolir == ip_solar).AND.                               &
                (i_sph_algorithm == ip_sph_reduced_iter)
    l_rescale_solar_phf=l_rescale.AND.l_solar_phf
!   Calculate the solar scattering angles if treating the
!   solar beam separately.
    IF (l_solar_phf) THEN
! DEPENDS ON: sol_scat_cos
      CALL sol_scat_cos(n_profile, n_direction                          &
        , zen_0, direction, cos_sol_view                                &
        , nd_profile, nd_direction)
    END IF

!   Calculate Clebsch-Gordan coefficients once and for all.
! DEPENDS ON: calc_cg_coeff
    CALL calc_cg_coeff(ls_max_order                                     &
      , ia_sph_mm, ms_min, ms_trunc                                     &
      , cg_coeff                                                        &
      , nd_max_order, nd_sph_coeff)

!   Calculate spherical harmonics at polar angles of pi/2 for
!   use in Marshak's boundary conditions.
! DEPENDS ON: calc_uplm_zero
    CALL calc_uplm_zero(ms_min, ms_max, ia_sph_mm                       &
      , ls_local_trunc, uplm_zero                                       &
      , nd_max_order, nd_sph_coeff)

    IF (isolir == ip_solar) THEN
!     Calculate the spherical harmonics of the solar direction.
! DEPENDS ON: calc_uplm_sol
      CALL calc_uplm_sol(n_profile, ms_min, ms_max, ia_sph_mm           &
        , ls_local_trunc, zen_0, uplm_sol                               &
        , nd_profile, nd_max_order, nd_sph_coeff)
    END IF

    IF (i_sph_algorithm == ip_sph_reduced_iter) THEN
!     Calcuate some arrays of terms for the BRDF.
! DEPENDS ON: calc_brdf
      CALL calc_brdf(isolir, ms_min, ms_max, ia_sph_mm                  &
        , uplm_sol, uplm_zero                                           &
        , n_brdf_basis_fnc, ls_brdf_trunc, f_brdf                       &
        , n_profile, n_direction, direction                             &
        , brdf_sol, brdf_hemi                                           &
        , nd_profile, nd_radiance_profile, nd_direction                 &
        , nd_max_order, nd_sph_coeff                                    &
        , nd_brdf_basis_fnc, nd_brdf_trunc)
    END IF

!   For the calculation of equivalent extinction in the IR
!   we need the diffuse albedo for each basis function.
    l_diff_alb=.FALSE.
    DO i_band=1, n_band
      l_diff_alb=l_diff_alb.OR.                                         &
        (i_gas_overlap(i_band) == ip_overlap_k_eqv)
    END DO
    IF ( (isolir == ip_infra_red).AND.l_diff_alb ) THEN
! DEPENDS ON: diff_albedo_basis
      CALL diff_albedo_basis(n_brdf_basis_fnc                           &
        , ls_brdf_trunc, f_brdf                                         &
        , uplm_zero(ia_sph_mm(0))                                       &
        , diffuse_alb_basis                                             &
        , nd_brdf_basis_fnc, nd_brdf_trunc, nd_sph_coeff                &
        )
    END IF

!   Determine which layers will be required to give radiances.
! DEPENDS ON: set_rad_layer
    CALL set_rad_layer(ierr                                             &
      , n_layer, n_viewing_level, viewing_level                         &
      , i_rad_layer, frac_rad_layer                                     &
      , nd_viewing_level                                                &
      )


  END IF

! Set the top level of the profiles. This is currently reatined
! for historical reasons.
  i_top=1



! Initial calculations for aerosols:

! Set the spectrally independent properties of moist aerosols.
  l_moist_aerosol=.FALSE.
  DO j=1, n_aerosol
    l_moist_aerosol=l_moist_aerosol.OR.                                 &
      (i_aerosol_parametrization(j)                                     &
       == ip_aerosol_param_moist).OR.                                   &
      (i_aerosol_parametrization(j)                                     &
       == ip_aerosol_param_phf_moist)
  END DO

  IF (l_moist_aerosol) THEN
! DEPENDS ON: set_moist_aerosol_properties
    CALL set_moist_aerosol_properties(ierr                              &
      , n_profile, n_layer                                              &
      , n_aerosol, i_aerosol_parametrization, nhumidity                 &
      , gas_mix_ratio(1, 1, index_water), t, p, w_cloud                 &
      , delta_humidity, mean_rel_humidity, i_humidity_pointer           &
      , nd_profile, nd_layer, id_ct, nd_aerosol_species                 &
      )
  END IF


! Check whether the densities will be needed for
! unparametrized aerosols.
  l_aerosol_density=.FALSE.
  IF (l_aerosol) THEN
    DO j=1, n_aerosol
      l_aerosol_density=l_aerosol_density.OR.                           &
        (i_aerosol_parametrization(j) ==                                &
         ip_aerosol_param_moist)                                        &
         .OR.(i_aerosol_parametrization(j) ==                           &
         ip_aerosol_unparametrized)
    END DO
  END IF



! Initial calculations for clouds:

  IF (l_cloud) THEN

!   Set pointers to the types of cloud.
! DEPENDS ON: set_cloud_pointer
    CALL set_cloud_pointer(ierr                                         &
      , n_condensed, type_condensed, i_cloud_representation             &
      , l_drop, l_ice                                                   &
      , i_phase_cmp, i_cloud_type, l_cloud_cmp                          &
      , nd_cloud_component                                              &
      )


!   Set the geometry of the clouds.
! DEPENDS ON: set_cloud_geometry
    CALL set_cloud_geometry(n_profile, n_layer                          &
      , l_global_cloud_top, n_cloud_top_global, w_cloud                 &
      , n_cloud_top, n_cloud_profile, i_cloud_profile                   &
      , nd_profile, nd_layer, id_ct                                     &
      )

!   Scale the condensed water contents to simulate
!   inhomogeneities in the clouds.
    IF ( (l_inhom_cloud) .AND. ( i_cloud /= ip_cloud_mcica ) ) THEN
      IF ( ALLOCATED(cloud_inhom_param) ) THEN
        DO k=1,nd_cloud_component
          DO j=n_cloud_top,n_layer
            DO i=1,n_profile
              condensed_mix_ratio(i,j,k)=cloud_inhom_param(i,j)         &
               *condensed_mix_ratio(i,j,k)
            END DO
          END DO           
        END DO
      ELSE
        DO i=1,nd_cloud_component
          condensed_mix_ratio(1:n_profile,n_cloud_top:n_layer,i)        &
           =inhom_cloud(i) *                                            &
            condensed_mix_ratio(1:n_profile,n_cloud_top:n_layer,i)
        END DO
      END IF
    END IF

    k_clr=1
    IF ( (i_cloud == ip_cloud_triple).OR.                               &
         (i_cloud == ip_cloud_part_corr_cnv) ) THEN
!     Aggregate clouds into regions for solving.
!     Three regions are used with this option. Additionally,
!     flag the clear-sky region.
      n_region=3
      type_region(1)=ip_region_clear
      type_region(2)=ip_region_strat
      type_region(3)=ip_region_conv
! DEPENDS ON: aggregate_cloud
      CALL aggregate_cloud(ierr                                         &
        , n_profile, n_layer, n_cloud_top                               &
        , i_cloud, i_cloud_representation                               &
        , n_cloud_type, frac_cloud                                      &
        , i_region_cloud, frac_region                                   &
        , nd_profile, nd_layer, nd_cloud_type, nd_region                &
        , id_ct                                                         &
        )
    ELSE IF ( (i_cloud == ip_cloud_mix_max).OR.                         &
              (i_cloud == ip_cloud_mix_random).OR.                      &
              (i_cloud == ip_cloud_part_corr) ) THEN
!     There will be only one cloudy region.
      n_region=2
      type_region(1)=ip_region_clear
      type_region(2)=ip_region_strat
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          frac_region(l, i, 2)=1.0e+00_RealK
        END DO
      END DO
    END IF

!   Calculate energy transfer coefficients in a mixed column,
!   or split the atmosphere into columns with a column model:
    IF ( (i_cloud == ip_cloud_mix_max).OR.                              &
         (i_cloud == ip_cloud_mix_random).OR.                           &
         (i_cloud == ip_cloud_triple).OR.                               &
         (i_cloud == ip_cloud_part_corr).OR.                            &
         (i_cloud == ip_cloud_part_corr_cnv) ) THEN

! DEPENDS ON: overlap_coupled
      CALL overlap_coupled(n_profile, n_layer, n_cloud_top              &
        , w_cloud, w_free, n_region, type_region, frac_region, p        &
        , i_cloud                                                       &
        , cloud_overlap                                                 &
        , nd_profile, nd_layer, nd_overlap_coeff, nd_region             &
        , id_ct, dp_corr_strat, dp_corr_conv, tot_cloud_cover           &
        )

    ELSE IF (i_cloud == ip_cloud_column_max) THEN

! DEPENDS ON: cloud_maxcs_split
        CALL cloud_maxcs_split(ierr, n_profile, n_layer, n_cloud_top    &
          , w_cloud, frac_cloud                                         &
          , n_cloud_type                                                &
          , n_column_cld, n_column_slv, list_column_slv                 &
          , i_clm_lyr_chn, i_clm_cld_typ, area_column                   &
          , nd_profile, nd_layer, id_ct, nd_column, nd_cloud_type       &
          )

    END IF

  ELSE

    n_cloud_top=n_layer+1

  END IF


! Calculate the atmospheric densities:
  IF ( l_continuum                                                      &
        .OR.l_aerosol_density                                           &
        .OR.(l_cloud                                                    &
! DEPENDS ON: l_cloud_density
        .AND.l_cloud_density(n_condensed, i_phase_cmp, l_cloud_cmp      &
                          , i_condensed_param, nd_cloud_component       &
                          ) ) ) THEN

!   Set the pointer for water vapour to a legal value: this must
!   be done for cases where water vapour is not included in the
!   spectral file, but densities are needed for aerosols.
    i_pointer_water=MAX(index_water, 1)

! DEPENDS ON: calculate_density
    CALL calculate_density(n_profile, n_layer, l_continuum              &
      , gas_mix_ratio(1, 1, i_pointer_water)                            &
      , p, t, i_top                                                     &
      , density, molar_density_water, molar_density_frn                 &
      , nd_profile, nd_layer                                            &
      )
  END IF

! Calculate temperature and pressure interpolation factor for ESFT
  IF (ANY(i_scale_fnc(i_first_band:i_last_band,1)==ip_scale_ses2)) THEN
! DEPENDS ON: inter_pt
    CALL inter_pt(nd_profile, nd_layer                                  &
      , n_profile, n_layer, gas_mix_ratio(1,1,index_water)              &
      , p, t, fac00, fac01, fac10, fac11                                &
      , fac00c, fac01c, fac10c, fac11c                                  &
      , facc00, facc01, jp, jph2oc, jt, jtt, jto2c, jtswo3)
  END IF

! Initialise solar_tail_flux if necessary
  IF (l_solar_tail_flux ) solar_tail_flux=0.0

! Check that there is enough information in the case of spherical
! harmonics. This check is rather late in the logical order of
! things, but we had to wait for certain other calculations to be
! made.
  IF (i_angular_integration == ip_spherical_harmonic) THEN
! DEPENDS ON: check_phf_term
    CALL check_phf_term(ierr                                            &
      , l_aerosol, n_aerosol, i_aerosol_parametrization                 &
      , n_aerosol_phf_term                                              &
      , l_cloud, n_condensed, i_condensed_param, i_phase_cmp            &
      , condensed_n_phf                                                 &
      , n_order_phase, l_henyey_greenstein_pf                           &
      , l_rescale, n_order_forward                                      &
      , l_solar_phf, n_order_phase_solar                                &
      , nd_aerosol_species, nd_cloud_component                          &
      )
  END IF


! Initialization of cloud diagnostics for the unified model.
  IF (l_cloud_extinction) THEN
     cloud_extinction        = 0.0
     cloud_weight_extinction = 0.0
  END IF

  IF (l_ls_cloud_extinction) THEN
     ls_cloud_extinction        = 0.0
     ls_cloud_weight_extinction = 0.0
  END IF

  IF (l_cnv_cloud_extinction) THEN
     cnv_cloud_extinction        = 0.0
     cnv_cloud_weight_extinction = 0.0
  END IF

  IF (l_cloud_absorptivity) THEN
     cloud_absorptivity        = 0.0
     cloud_weight_absorptivity = 0.0
  END IF

  IF (l_ls_cloud_absorptivity) THEN
     ls_cloud_absorptivity        = 0.0
     ls_cloud_weight_absorptivity = 0.0
  END IF

  IF (l_cnv_cloud_absorptivity) THEN
     cnv_cloud_absorptivity        = 0.0
     cnv_cloud_weight_absorptivity = 0.0
  END IF


! For the CLASSIC aerosols, all AOD calculations are dealt with via a
! wrapper subroutine
  IF ( l_aod_sulphate .OR. l_aod_dust .OR. l_aod_seasalt .OR.           &
       l_aod_soot .OR. l_aod_biomass .OR. l_aod_biogenic .OR.           &
       l_aod_ocff .OR. l_aod_delta .OR. l_aod_nitrate .OR.              &
       l_aod_total_radn .OR. l_angst_total_radn .OR.                    &
       l_aod_prog_sulphate .OR. l_aod_prog_dust .OR.                    &
       l_aod_prog_seasalt .OR. l_aod_prog_soot .OR.                     &
       l_aod_prog_biomass .OR. l_aod_prog_ocff .OR.                     &
       l_aod_prog_nitrate ) THEN
! DEPENDS ON: compute_all_aod
    CALL compute_all_aod(ierr,                                          &
      n_aerosol,      n_aerosol_mr,                                     &
      nd_aerosol_species, nd_aerosol_mixratio,                          &
      n_aod_wavel,    nd_aod_wavel,  aod_wavel,                         &
      n_profile,      nd_profile,                                       &
      n_layer, 1,     nd_layer,                                         &
      nd_humidities,  nd_humidities,                                    &
      type_aerosol,   i_aod_type,                                       &
      aerosol_mix_ratio, d_mass,                                        &
      aerosol_mr_source, aerosol_mr_type_index,                         &
      i_aerosol_parametrization, aod_absorption, aod_scattering,        &
      i_humidity_pointer, mean_rel_humidity,                            &
      humidities,     delta_humidity,    l_use_arcl,                    &
      l_aod_sulphate, aod_sulphate, l_aod_dust,      aod_dust,          &
      l_aod_seasalt,  aod_seasalt,  l_aod_soot,      aod_soot,          &
      l_aod_biomass,  aod_biomass,  l_aod_biogenic,  aod_biogenic,      &
      l_aod_ocff,     aod_ocff,     l_aod_delta,     aod_delta,         &
      l_aod_nitrate,  aod_nitrate,  l_aod_total_radn,aod_total_radn,    &
      l_angst_total_radn,   angst_total_radn,                           &
      l_aod_prog_sulphate,  aod_prog_sulphate,   l_aod_prog_dust,       &
      aod_prog_dust,        l_aod_prog_seasalt,  aod_prog_seasalt,      &
      l_aod_prog_soot,      aod_prog_soot,       l_aod_prog_biomass,    &
      aod_prog_biomass,     l_aod_prog_ocff,     aod_prog_ocff,         &
      l_aod_prog_nitrate,   aod_prog_nitrate)
  END IF


! For UKCA aerosols, perform the waveband averaging of optical
! properties.
! Note: this is done outside the loop over bands in order to
!       support excluded bands.

  IF (l_ukca_radaer) THEN


       ! Calculations shared by ukca_radaer_band_average() and
       ! ukca_radaer_compute_aod().

! DEPENDS ON: ukca_radaer_prepare
       CALL ukca_radaer_prepare(                                        &
         ! Actual array dimensions
         n_profile, n_layer, n_ukca_mode, n_ukca_cpnt,                  &
         ! UKCA_RADAER structure
         ukca_radaer,                                                   &
         ! Component mass-mixing ratios
         ukca_mix_ratio,                                                &
         ! Modal mass-mixing ratios
         ukca_modal_mixr,                                               &
         ! Input modal number concentrations
         ukca_modal_nbr,                                                &
         ! Output modal number concentrations
         ukca_modal_number,                                             &
         ! Pressure and temperature
         p, t,                                                          &
         ! Fixed array dimensions
         nd_profile, nd_layer)


       ! Compute the band-averaged optical properties for UKCA aerosols

! DEPENDS ON: ukca_radaer_band_average
       CALL ukca_radaer_band_average(                                   &
         ! Spectral information
         n_band, isolir, l_exclude, n_band_exclude, index_exclude,      &
         ! Actual array dimensions
         n_profile, n_layer, n_ukca_mode, n_ukca_cpnt,                  &
         ! UKCA_RADAER structure
         ukca_radaer,                                                   &
         ! Modal mass-mixing ratios
         ukca_modal_mixr,                                               &
         ! Modal number concentrations
         ukca_modal_number,                                             &
         ! Modal diameters from UKCA module
         ukca_dry_diam, ukca_wet_diam,                                  &
         ! Other inputs from UKCA module
         ukca_comp_vol, ukca_modal_vol, ukca_modal_rho,                 &
         ukca_modal_wtv,                                                &
         ! Band-averaged optical properties (outputs)
         ukca_absorption, ukca_scattering, ukca_asymmetry,              &
         ! Fixed array dimensions
         nd_profile, nd_layer, nd_band, nd_exclude)

  END IF ! l_ukca_radaer

  IF (l_aod_ukca_ait_sol) THEN

! DEPENDS ON: ukca_radaer_compute_aod
       CALL ukca_radaer_compute_aod(                                    &
              ! UKCA model switch
              l_ukca_radaer,                                            &
              ! Actual array dimension
              n_profile, n_layer, n_ukca_mode,                          &
              n_ukca_cpnt, n_aod_wavel,                                 &
              ! UKCA_RADAER structure
              ukca_radaer,                                              &
              ! Modal diameters from UKCA module
              ukca_dry_diam, ukca_wet_diam,                             &
              ! Mass thickness of layers
              d_mass,                                                   &
              ! Component volumes
              ukca_comp_vol,                                            &
              ! Modal volumes, densities, and water content
              ukca_modal_vol, ukca_modal_rho, ukca_modal_wtv,           &
              ! Modal mass-mixing ratios
              ukca_modal_mixr,                                          &
              ! Modal number concentrations
              ukca_modal_number,                                        &
              ! Type selection
              ip_ukca_mode_aitken, .TRUE.,                              &
              ! Modal extinction aerosol opt depth (output)
              aod_ukca_ait_sol,                                         &
              ! Fixed array dimensions
              nd_profile, nd_layer, nd_aod_wavel)

  END IF ! l_aod_ukca_ait_sol

  IF (l_aod_ukca_acc_sol) THEN

! DEPENDS ON: ukca_radaer_compute_aod
       CALL ukca_radaer_compute_aod(                                    &
              ! UKCA model switch
              l_ukca_radaer,                                            &
              ! Actual array dimension
              n_profile, n_layer, n_ukca_mode,                          &
              n_ukca_cpnt, n_aod_wavel,                                 &
              ! UKCA_RADAER structure
              ukca_radaer,                                              &
              ! Modal diameters from UKCA module
              ukca_dry_diam, ukca_wet_diam,                             &
              ! Mass thickness of layers
              d_mass,                                                   &
              ! Component volumes
              ukca_comp_vol,                                            &
              ! Modal volumes, densities, and water content
              ukca_modal_vol, ukca_modal_rho, ukca_modal_wtv,           &
              ! Modal mass-mixing ratios
              ukca_modal_mixr,                                          &
              ! Modal number concentrations
              ukca_modal_number,                                        &
              ! Type selection
              ip_ukca_mode_accum, .TRUE.,                               &
              ! Modal extinction aerosol opt depth (output)
              aod_ukca_acc_sol,                                         &
              ! Fixed array dimensions
              nd_profile, nd_layer, nd_aod_wavel)

  END IF ! l_aod_ukca_acc_sol

  IF (l_aod_ukca_cor_sol) THEN

! DEPENDS ON: ukca_radaer_compute_aod
       CALL ukca_radaer_compute_aod(                                    &
              ! UKCA model switch
              l_ukca_radaer,                                            &
              ! Actual array dimension
              n_profile, n_layer, n_ukca_mode,                          &
              n_ukca_cpnt, n_aod_wavel,                                 &
              ! UKCA_RADAER structure
              ukca_radaer,                                              &
              ! Modal diameters from UKCA module
              ukca_dry_diam, ukca_wet_diam,                             &
              ! Mass thickness of layers
              d_mass,                                                   &
              ! Component volumes
              ukca_comp_vol,                                            &
              ! Modal volumes, densities, and water content
              ukca_modal_vol, ukca_modal_rho, ukca_modal_wtv,           &
              ! Modal mass-mixing ratios
              ukca_modal_mixr,                                          &
              ! Modal number concentrations
              ukca_modal_number,                                        &
              ! Type selection
              ip_ukca_mode_coarse, .TRUE.,                              &
              ! Modal extinction aerosol opt depth (output)
              aod_ukca_cor_sol,                                         &
              ! Fixed array dimensions
              nd_profile, nd_layer, nd_aod_wavel)

  END IF ! l_aod_ukca_cor_sol

  IF (l_aod_ukca_ait_ins) THEN

! DEPENDS ON: ukca_radaer_compute_aod
       CALL ukca_radaer_compute_aod(                                    &
              ! UKCA model switch
              l_ukca_radaer,                                            &
              ! Actual array dimension
              n_profile, n_layer, n_ukca_mode,                          &
              n_ukca_cpnt, n_aod_wavel,                                 &
              ! UKCA_RADAER structure
              ukca_radaer,                                              &
              ! Modal diameters from UKCA module
              ukca_dry_diam, ukca_wet_diam,                             &
              ! Mass thickness of layers
              d_mass,                                                   &
              ! Component volumes
              ukca_comp_vol,                                            &
              ! Modal volumes, densities, and water content
              ukca_modal_vol, ukca_modal_rho, ukca_modal_wtv,           &
              ! Modal mass-mixing ratios
              ukca_modal_mixr,                                          &
              ! Modal number concentrations
              ukca_modal_number,                                        &
              ! Type selection
              ip_ukca_mode_aitken, .FALSE.,                             &
              ! Modal extinction aerosol opt depth (output)
              aod_ukca_ait_ins,                                         &
              ! Fixed array dimensions
              nd_profile, nd_layer, nd_aod_wavel)

  END IF ! l_aod_ukca_ait_ins

  IF (l_aod_ukca_acc_ins) THEN

! DEPENDS ON: ukca_radaer_compute_aod
       CALL ukca_radaer_compute_aod(                                    &
              ! UKCA model switch
              l_ukca_radaer,                                            &
              ! Actual array dimension
              n_profile, n_layer, n_ukca_mode,                          &
              n_ukca_cpnt, n_aod_wavel,                                 &
              ! UKCA_RADAER structure
              ukca_radaer,                                              &
              ! Modal diameters from UKCA module
              ukca_dry_diam, ukca_wet_diam,                             &
              ! Mass thickness of layers
              d_mass,                                                   &
              ! Component volumes
              ukca_comp_vol,                                            &
              ! Modal volumes, densities, and water content
              ukca_modal_vol, ukca_modal_rho, ukca_modal_wtv,           &
              ! Type selection
              ip_ukca_mode_accum, .FALSE.,                              &
              ! Modal mass-mixing ratios
              ukca_modal_mixr,                                          &
              ! Modal number concentrations
              ukca_modal_number,                                        &
              ! Modal extinction aerosol opt depth (output)
              aod_ukca_acc_ins,                                         &
              ! Fixed array dimensions
              nd_profile, nd_layer, nd_aod_wavel)

  END IF ! l_aod_ukca_acc_ins

  IF (l_aod_ukca_cor_ins) THEN

! DEPENDS ON: ukca_radaer_compute_aod
       CALL ukca_radaer_compute_aod(                                    &
              ! UKCA model switch
              l_ukca_radaer,                                            &
              ! Actual array dimension
              n_profile, n_layer, n_ukca_mode,                          &
              n_ukca_cpnt, n_aod_wavel,                                 &
              ! UKCA_RADAER structure
              ukca_radaer,                                              &
              ! Modal diameters from UKCA module
              ukca_dry_diam, ukca_wet_diam,                             &
              ! Mass thickness of layers
              d_mass,                                                   &
              ! Component volumes
              ukca_comp_vol,                                            &
              ! Modal volumes, densities, and water content
              ukca_modal_vol, ukca_modal_rho, ukca_modal_wtv,           &
              ! Modal mass-mixing ratios
              ukca_modal_mixr,                                          &
              ! Modal number concentrations
              ukca_modal_number,                                        &
              ! Type selection
              ip_ukca_mode_coarse, .FALSE.,                             &
              ! Modal extinction aerosol opt depth (output)
              aod_ukca_cor_ins,                                         &
              ! Fixed array dimensions
              nd_profile, nd_layer, nd_aod_wavel)

  END IF ! l_aod_ukca_cor_ins




! Solve the equation of transfer in each band and
! increment the fluxes.

  DO i_band=i_first_band, i_last_band

!   Set the flag to initialize the diagnostic arrays.
    IF (i_band == i_first_band) THEN
      l_initial=.TRUE.
    ELSE
      l_initial=(map_channel(i_band) > map_channel(i_band-1))
    END IF

!   Initialise fluxes to calculate band increments for
!   cloud diagnostics
    IF (l_cloud_extinction.OR.l_ls_cloud_extinction                     &
                          .OR.l_cnv_cloud_extinction) THEN
      IF (l_initial) THEN
        flux_direct_clear_band = 0.0
      ELSE
        flux_direct_clear_band                                          &
          = flux_direct_clear(:,:,map_channel(i_band))
      END IF
    END IF
    IF (l_cloud_absorptivity.OR.l_ls_cloud_absorptivity                 &
                            .OR.l_cnv_cloud_absorptivity) THEN
      IF (l_initial) THEN
        flux_up_clear_band = 0.0
      ELSE
        flux_up_clear_band                                              &
          = flux_up_clear(:,:,map_channel(i_band))
      END IF
    END IF

!   Determine whether clear-sky fluxes are required for this band
    l_clear_band = l_clear
    IF ((isolir == ip_solar).AND.(i_band == i_first_band)) THEN
!     As currently implemented UV flux must be in band 1
      l_clear_band = l_clear_band .OR. l_surf_uv_clr
    END IF

!   Determine whether gaseous absorption is included in this band.
    IF ( (l_gas).AND.(n_band_absorb(i_band) > 0) ) THEN

!     Note: I_GAS_BAND is used extensively below since nested
!     array elements in a subroutine call (see later) can
!     confuse some compilers.

!     Normally the number of gases in the calculation will be
!     as in the spectral file, but particular options may result
!     in the omission of some gases.

      n_gas=n_band_absorb(i_band)

      IF (i_gas_overlap(i_band) == ip_overlap_single) THEN

!       There will be no gaseous absorption in this band
!       unless the selected gas appears.
        n_gas=0

        DO i=1, n_band_absorb(i_band)
          IF (index_absorb(i, i_band) == i_gas) n_gas=1
        END DO

      END IF


      IF (n_gas > 0) THEN

!       Set the flag for gaseous absorption in the band.
        l_gas_band=.TRUE.

        DO j=1, n_gas

          i_gas_band=index_absorb(j, i_band)

!         Reset the pointer if there is just one gas.

          IF (i_gas_overlap(i_band) == ip_overlap_single)               &
            THEN
!           Only the selected gas is active in the band.
            i_gas_band=i_gas

          END IF

          IF (i_scale_esft(i_band, i_gas_band)                          &
              == ip_scale_band) THEN
!           Rescale the amount of gas for this band now.
! DEPENDS ON: scale_absorb
            CALL scale_absorb(ierr, n_profile, n_layer                  &
              , gas_mix_ratio(1, 1, i_gas_band), p, t                   &
              , i_top                                                   &
              , gas_frac_rescaled(1, 1, i_gas_band)                     &
              , i_scale_fnc(i_band, i_gas_band)                         &
              , p_reference(i_gas_band, i_band)                         &
              , t_reference(i_gas_band, i_band)                         &
              , scale_vector(1, 1, i_band, i_gas_band)                  &
              , 1, i_band                                               &
              , l_doppler(i_gas_band)                                   &
              , doppler_correction(i_gas_band)                          &
              , nd_profile, nd_layer                                    &
              , nd_scale_variable                                       &
              )

          ELSE IF (i_scale_esft(i_band, i_gas_band)                     &
              == ip_scale_null) THEN
!           Copy across the unscaled array.
            DO i=i_top, n_layer
              DO l=1, n_profile
                gas_frac_rescaled(l, i, i_gas_band)                     &
                  =gas_mix_ratio(l, i, i_gas_band)
              END DO
            END DO
          END IF
        END DO
      ELSE
        l_gas_band=.FALSE.
      END IF

    ELSE
      l_gas_band=.FALSE.
    END IF



!   Rescale amounts of continua.
    l_grey_cont = .FALSE.
    IF (l_continuum) THEN
      n_continuum=n_band_continuum(i_band)
      DO i=1, n_continuum
        i_continuum_pointer(i)=index_continuum(i_band, i)
        i_continuum=i_continuum_pointer(i)
        IF (i_scale_fnc_cont(i_band,i_continuum) == ip_scale_ses2) THEN
          k_continuum_mono(i_continuum)=0.0
! DEPENDS ON: ses_rescale_contm
          CALL ses_rescale_contm(nd_profile, nd_layer                   &
            , i_continuum, n_profile, n_layer                           &
            , p, t, gas_mix_ratio(1,1,index_water)                      &
            , amount_continuum(1, 1, i)                                 &
            )
        ELSE
          l_grey_cont = .TRUE.
          k_continuum_mono(i_continuum)=k_continuum(i_band, i_continuum)
! DEPENDS ON: rescale_continuum
          CALL rescale_continuum(n_profile, n_layer, i_continuum        &
            , p, t, i_top                                               &
            , density, molar_density_water, molar_density_frn           &
            , gas_mix_ratio(1, 1, index_water)                          &
            , amount_continuum(1, 1, i_continuum)                       &
            , i_scale_fnc_cont(i_band, i_continuum)                     &
            , p_ref_continuum(i_continuum, i_band)                      &
            , t_ref_continuum(i_continuum, i_band)                      &
            , scale_continuum(1, i_band, i_continuum)                   &
            , nd_profile, nd_layer                                      &
            , nd_scale_variable                                         &
            )
        END IF
      END DO
    END IF

!   Allocate the single scattering propeties.
    ALLOCATE(ss_prop%k_grey_tot_clr                                     &
      (nd_profile, 1:nd_layer_clr))
    ALLOCATE(ss_prop%k_ext_scat_clr                                     &
      (nd_profile, 1:nd_layer_clr))
    ALLOCATE(ss_prop%phase_fnc_clr                                      &
      (nd_profile, 1:nd_layer_clr, nd_max_order))
    ALLOCATE(ss_prop%forward_scatter_clr                                &
      (nd_profile, 1:nd_layer_clr))
    ALLOCATE(ss_prop%phase_fnc_solar_clr                                &
      (nd_radiance_profile, 1:nd_layer_clr, nd_direction))
    ALLOCATE(ss_prop%forward_solar_clr                                  &
      (nd_radiance_profile, 1:nd_layer_clr))

    ALLOCATE(ss_prop%tau_clr                                            &
      (nd_profile, 1: nd_layer_clr))
    ALLOCATE(ss_prop%omega_clr                                          &
      (nd_profile, 1: nd_layer_clr))

    ALLOCATE(ss_prop%k_grey_tot                                         &
      (nd_profile, id_ct: nd_layer,                                     &
       0: nd_cloud_type))
    ALLOCATE(ss_prop%k_ext_scat                                         &
      (nd_profile, id_ct: nd_layer,                                     &
       0: nd_cloud_type))
    ALLOCATE(ss_prop%phase_fnc                                          &
      (nd_profile, id_ct: nd_layer, nd_max_order,                       &
       0: nd_cloud_type))
    ALLOCATE(ss_prop%forward_scatter                                    &
      (nd_profile, id_ct: nd_layer,                                     &
       0: nd_cloud_type))
    ALLOCATE(ss_prop%phase_fnc_solar                                    &
      (nd_radiance_profile, id_ct: nd_layer, nd_direction,              &
       0: nd_cloud_type))
    ALLOCATE(ss_prop%forward_solar                                      &
      (nd_radiance_profile, id_ct: nd_layer,                            &
       0: nd_cloud_type))

    ALLOCATE(ss_prop%tau                                                &
      (nd_profile, id_ct: nd_layer, 0: nd_cloud_type))
    ALLOCATE(ss_prop%omega                                              &
      (nd_profile, id_ct: nd_layer, 0: nd_cloud_type))

    ALLOCATE(ss_prop%k_ext_tot_cloud_comp                               &
      (nd_profile, id_ct: nd_layer,                                     &
       0: nd_cloud_component))
    ALLOCATE(ss_prop%k_ext_scat_cloud_comp                              &
      (nd_profile, id_ct: nd_layer,                                     &
       0: nd_cloud_component))
    ALLOCATE(ss_prop%phase_fnc_cloud_comp                               &
      (nd_profile, id_ct: nd_layer, nd_max_order,                       &
       0: nd_cloud_component))
    ALLOCATE(ss_prop%forward_scatter_cloud_comp                         &
      (nd_profile, id_ct: nd_layer,                                     &
       0: nd_cloud_component))
    ALLOCATE(ss_prop%phase_fnc_solar_cloud_comp                         &
      (nd_radiance_profile, id_ct: nd_layer, nd_direction,              &
       0: nd_cloud_component))
    ALLOCATE(ss_prop%forward_solar_cloud_comp                           &
      (nd_radiance_profile, id_ct: nd_layer,                            &
       0: nd_cloud_component))

    ALLOCATE(ss_prop%phase_fnc_no_cloud                                 &
      (nd_profile, id_ct: nd_layer, nd_max_order))
    ALLOCATE(ss_prop%forward_scatter_no_cloud                           &
      (nd_profile, id_ct: nd_layer))

!   Calculate the grey extinction within the band.

! DEPENDS ON: grey_opt_prop
    CALL grey_opt_prop(ierr                                             &
      , n_profile, n_layer, p, t, density                               &
      , n_order_phase, l_rescale, n_order_forward                       &
      , l_henyey_greenstein_pf, l_solar_phf, l_lanczos                  &
      , n_order_phase_solar, n_direction, cos_sol_view                  &
      , l_rayleigh, rayleigh_coefficient(i_band)                        &
      , l_grey_cont, n_continuum, i_continuum_pointer                   &
      , k_continuum_mono, amount_continuum                              &
      , l_aerosol, n_aerosol, n_aerosol_mr, aerosol_mix_ratio           &
      , aerosol_mr_source, aerosol_mr_type_index                        &
      , i_aerosol_parametrization                                       &
      , i_humidity_pointer, humidities, delta_humidity                  &
      , mean_rel_humidity                                               &
      , aerosol_absorption(1, 1, i_band)                                &
      , aerosol_scattering(1, 1, i_band)                                &
      , aerosol_phase_fnc(1, 1, 1, i_band)                              &
      , l_ukca_radaer, n_ukca_mode, ukca_modal_mixr                     &
      , ukca_absorption(1, 1, 1, i_band)                                &
      , ukca_scattering(1, 1, 1, i_band)                                &
      , ukca_asymmetry(1, 1, 1, i_band)                                 &
      , l_cloud, i_cloud, n_cloud_profile, i_cloud_profile              &
      , n_cloud_top, n_condensed, l_cloud_cmp, i_phase_cmp              &
      , i_condensed_param, condensed_n_phf                              &
      , condensed_param_list(1, 1, i_band)                              &
      , condensed_mix_ratio, condensed_dim_char                         &
      , n_cloud_type, i_cloud_type                                      &
      , ss_prop                                                         &
      , frac_cloud                                                      &
      , l_cloud_extinction, cloud_extinction_band                       &
      , l_cloud_absorptivity, cloud_absorptivity_band                   &
      , l_ls_cloud_extinction, ls_cloud_extinction_band                 &
      , l_ls_cloud_absorptivity, ls_cloud_absorptivity_band             &
      , l_cnv_cloud_extinction, cnv_cloud_extinction_band               &
      , l_cnv_cloud_absorptivity, cnv_cloud_absorptivity_band           &
      , nd_profile, nd_radiance_profile, nd_layer                       &
      , nd_layer_clr, id_ct                                             &
      , nd_continuum, nd_aerosol_species, nd_aerosol_mixratio           &
      , nd_humidities, nd_cloud_parameter, nd_cloud_component           &
      , nd_cloud_type, nd_phase_term, nd_max_order, nd_direction        &
      , n_ukca_mode                                                     &
      )


    IF ( (i_angular_integration == ip_two_stream).OR.                   &
         (i_angular_integration == ip_spherical_harmonic) ) THEN

!     Rescale the phase function and calculate the scattering
!     fractions. (These are grey and may be calculated outside
!     a loop over gases).

      IF (l_rescale) THEN

!       Rescale clear-sky phase function:

!       The section above clouds.
! DEPENDS ON: rescale_phase_fnc
        CALL rescale_phase_fnc(n_profile, 1, n_cloud_top-1              &
          , n_direction, cos_sol_view                                   &
          , n_order_phase                                               &
          , ss_prop%phase_fnc_clr, ss_prop%forward_scatter_clr          &
          , ss_prop%forward_solar_clr                                   &
          , l_rescale_solar_phf, n_order_phase_solar                    &
          , ss_prop%phase_fnc_solar_clr                                 &
          , nd_profile, nd_radiance_profile, nd_layer_clr, 1            &
          , nd_direction, nd_max_order                                  &
          )
!       The section including clouds.
        CALL rescale_phase_fnc(n_profile, n_cloud_top                   &
          , n_layer, n_direction, cos_sol_view                          &
          , n_order_phase                                               &
          , ss_prop%phase_fnc(1, id_ct, 1, 0)                           &
          , ss_prop%forward_scatter(1, id_ct, 0)                        &
          , ss_prop%forward_solar(1, id_ct, 0)                          &
          , l_rescale_solar_phf, n_order_phase_solar                    &
          , ss_prop%phase_fnc_solar(1, id_ct, 1, 0)                     &
          , nd_profile, nd_radiance_profile, nd_layer, id_ct            &
          , nd_direction, nd_max_order                                  &
          )


        IF (l_cloud .AND. i_cloud /= ip_cloud_mcica) THEN

!         Rescale cloudy phase functions:
          DO k=1, n_cloud_type
            CALL rescale_phase_fnc(n_profile, n_cloud_top               &
              , n_layer, n_direction, cos_sol_view                      &
              , n_order_phase                                           &
              , ss_prop%phase_fnc(1, id_ct, 1, k)                       &
              , ss_prop%forward_scatter(1, id_ct, k)                    &
              , ss_prop%forward_solar(1, id_ct, k)                      &
              , l_rescale_solar_phf, n_order_phase_solar                &
              , ss_prop%phase_fnc_solar(1, id_ct, 1, k)                 &
              , nd_profile, nd_radiance_profile, nd_layer, id_ct        &
              , nd_direction, nd_max_order                              &
              )
          END DO

        END IF

      END IF

    END IF




!   Preliminary calculations for source terms:

    IF (i_gas_overlap(i_band) == ip_overlap_mix_ses2) THEN

!     Interpolate absorption coefficients onto model grid
! DEPENDS ON: inter_k
      CALL inter_k(n_profile, n_layer, n_band_absorb                    &
        , mix_gas_band, n_mix_gas(i_band), index_mix_gas                &
        , i_band_esft_ses(i_band), i_band                               &
        , n_continuum, k_continuum_ses, l_continuum, index_continuum    &
        , k_h2oc(1,1,1,i_band), fac00, fac01, fac10, fac11              &
        , fac00c, fac01c, fac10c, fac11c, facc00, facc01                &
        , jp, jph2oc, jt, jtt, jto2c, jtswo3                            &
        , k_esft_ses(1,1,1,1,i_band), k_mix_gas                         &
        , f_mix(i_band), gas_mix_ratio, gas_mix_amt                     &
        , k_esft_layer, k_mix_gas_layer, k_contm_layer                  &
!       Dimensions
        , nd_profile, nd_layer                                          &
        , nd_band, nd_species, nd_continuum                             &
        , nd_esft_term, nd_mix, nd_tmp, nd_pre                          &
        , nd_band_mix_gas                                               &
        )

      IF ((isolir == ip_solar) .OR. l_solar_tail_flux) THEN
!       Convert normalized band fluxes to actual energy fluxes.
        DO k=1,i_band_esft_ses(i_band)
          DO l=1, n_profile
            solar_irrad_band_ses(l,k)=solar_irrad(l)                    &
              *solar_flux_band_ses(k, i_band)
          END DO
        END DO
      END IF

      IF (l_solar_tail_flux .AND. isolir == ip_infra_red) THEN
!       Calculate solar tail flux and add it to the solar region
!       for diagnostic output
        DO k=1,i_band_esft_ses(i_band)
          DO l=1, n_profile
            solar_tail_flux(l)=solar_tail_flux(l)                       &
              +w_esft_ses(k,i_band)*solar_irrad_band_ses(l,k)
          END DO
        END DO
      END IF

    ELSE

      IF (isolir == ip_solar) THEN
!       Convert normalized band fluxes to actual energy fluxes.
        DO l=1, n_profile
          solar_irrad_band(l)=solar_irrad(l)                            &
            *solar_flux_band(i_band)
        END DO
      END IF

    END IF


    IF (isolir == ip_infra_red) THEN

!     Calculate the change in the thermal source function
!     across each layer for the infra-red part of the spectrum.

! DEPENDS ON: diff_planck_source
      CALL diff_planck_source(n_profile, n_layer                        &
        , n_deg_fit, thermal_coefficient(0, i_band)                     &
        , t_ref_planck, t_level, t_ground                               &
        , planck_flux_band, diff_planck_band                            &
        , planck_flux_ground                                            &
        , l_ir_source_quad, t, diff_planck_band_2                       &
        , i_angular_integration                                         &
        , n_viewing_level, i_rad_layer, frac_rad_layer                  &
        , planck_radiance_band                                          &
        , l_tile, n_point_tile, n_tile, list_tile                       &
        , frac_tile, t_tile, planck_flux_tile                           &
        , nd_profile, nd_layer, nd_thermal_coeff                        &
        , nd_radiance_profile, nd_viewing_level                         &
        , nd_point_tile, nd_tile                                        &
        )

    END IF




!   Call a solver appropriate to the presence of gases and
!   the overlap assumed:

    IF (.NOT.l_gas_band) THEN

!     There is no gaseous absorption. Solve for the
!     radiances directly.

! DEPENDS ON: solve_band_without_gas
      CALL solve_band_without_gas(ierr                                  &
!                 Atmospheric properties
        , n_profile, n_layer, d_mass                                    &
!                 Angular integration
        , i_angular_integration, i_2stream                              &
        , n_order_phase, l_rescale, n_order_gauss                       &
        , ms_min, ms_max, i_truncation, ls_local_trunc                  &
        , accuracy_adaptive, euler_factor, i_sph_algorithm              &
        , i_sph_mode                                                    &
!                 Precalculated angular arrays
        , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                      &
!                 Treatment of scattering
        , i_scatter_method(i_band)                                      &
!                 Options for solver
        , i_solver                                                      &
!                 Spectral region
        , isolir                                                        &
!                 Solar properties
        , zen_0, solar_irrad_band                                       &
!                 Infra-red properties
        , planck_flux_band(1, 0), planck_flux_band(1, n_layer)          &
        , diff_planck_band, l_ir_source_quad, diff_planck_band_2        &
!                 Surface properties
        , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb(1, 1, i_band)        &
        , f_brdf, brdf_sol, brdf_hemi                                   &
        , planck_flux_ground                                            &
!                 Tiling of the surface
        , l_tile, n_point_tile, n_tile, list_tile                       &
        , rho_alb_tile(1, 1, 1, i_band), planck_flux_tile               &
!                 Optical Properties
        , ss_prop                                                       &
!                 Cloudy properties
        , l_cloud, i_cloud                                              &
!                 Cloudy geometry
        , n_cloud_top                                                   &
        , n_cloud_type, frac_cloud                                      &
        , n_region, k_clr, i_region_cloud, frac_region                  &
        , w_free, w_cloud, cloud_overlap                                &
        , n_column_slv, list_column_slv                                 &
        , i_clm_lyr_chn, i_clm_cld_typ, area_column                     &
!                 Levels for calculating radiances
        , n_viewing_level, i_rad_layer, frac_rad_layer                  &
!                 Viewing Geometry
        , n_direction, direction                                        &
!                 Weighting factor for the band
        , weight_band(i_band), l_initial                                &
!                 Calculated fluxes
        , flux_direct(1, 0, map_channel(i_band))                        &
        , flux_down(1, 0, map_channel(i_band))                          &
        , flux_up(1, 0, map_channel(i_band))                            &
!                 Radiances
        , i_direct, radiance(1, 1, 1, map_channel(i_band))              &
!                 Rate of photolysis
        , photolysis(1, 1, map_channel(i_band))                         &
!                 Flags for clear-sky fluxes
        , l_clear_band, i_solver_clear                                  &
!                 Calculated clear-sky fluxes
        , flux_direct_clear(1, 0, map_channel(i_band))                  &
        , flux_down_clear(1, 0, map_channel(i_band))                    &
        , flux_up_clear(1, 0, map_channel(i_band))                      &
!                 Tiled Surface Fluxes
        , flux_up_tile(1, 1, map_channel(i_band))                       &
        , flux_up_blue_tile(1, 1, map_channel(i_band))                  &
!                 Special Surface Fluxes
        , l_blue_flux_surf, weight_blue(i_band)                         &
        , flux_direct_blue_surf                                         &
        , flux_down_blue_surf, flux_up_blue_surf                        &
!                 Dimensions of arrays
        , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column          &
        , nd_flux_profile, nd_radiance_profile, nd_j_profile            &
        , nd_cloud_type, nd_region, nd_overlap_coeff                    &
        , nd_max_order, nd_sph_coeff                                    &
        , nd_brdf_basis_fnc, nd_brdf_trunc                              &
        , nd_viewing_level, nd_direction                                &
        , nd_source_coeff, nd_point_tile, nd_tile                       &
        )


    ELSE

!     Gases are included.

!     Treat the gaseous overlaps as directed by
!     the overlap switch.

      IF (i_gas_overlap(i_band) == ip_overlap_single) THEN

! DEPENDS ON: solve_band_one_gas
        CALL solve_band_one_gas(ierr                                    &
!                 Atmospheric properties
          , n_profile, n_layer, i_top, p, t, d_mass                     &
!                 Angular integration
          , i_angular_integration, i_2stream                            &
          , n_order_phase, l_rescale, n_order_gauss                     &
          , ms_min, ms_max, i_truncation, ls_local_trunc                &
          , accuracy_adaptive, euler_factor                             &
          , i_sph_algorithm, i_sph_mode                                 &
!                 Precalculated angular arrays
          , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                    &
!                 Treatment of scattering
          , i_scatter_method(i_band)                                    &
!                 Options for solver
          , i_solver                                                    &
!                 Gaseous properties
          , i_band, i_gas                                               &
          , i_band_esft, i_scale_esft, i_scale_fnc                      &
          , k_esft, w_esft, scale_vector                                &
          , p_reference, t_reference                                    &
          , gas_mix_ratio, gas_frac_rescaled                            &
          , l_doppler, doppler_correction                               &
!                 Spectral region
          , isolir                                                      &
!                 Solar properties
          , zen_0, solar_irrad_band                                     &
!                 Infra-red properties
          , planck_flux_band(1, 0)                                      &
          , planck_flux_band(1, n_layer)                                &
          , diff_planck_band                                            &
          , l_ir_source_quad, diff_planck_band_2                        &
!                 Surface properties
          , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb(1, 1, i_band)      &
          , f_brdf, brdf_sol, brdf_hemi                                 &
          , planck_flux_ground                                          &
!                 Tiling of the surface
          , l_tile, n_point_tile, n_tile, list_tile                     &
          , rho_alb_tile(1, 1, 1, i_band)                               &
          , planck_flux_tile                                            &
!                 Optical Properties
          , ss_prop                                                     &
!                 Cloudy properties
          , l_cloud, i_cloud                                            &
!                 Cloud geometry
          , n_cloud_top                                                 &
          , n_cloud_type, frac_cloud                                    &
          , n_region, k_clr, i_region_cloud, frac_region                &
          , w_free, w_cloud, cloud_overlap                              &
          , n_column_slv, list_column_slv                               &
          , i_clm_lyr_chn, i_clm_cld_typ, area_column                   &
!                 Levels for calculating radiances
          , n_viewing_level, i_rad_layer, frac_rad_layer                &
!                 Viewing Geometry
          , n_direction, direction                                      &
!                 Weighting factor for the band
          , weight_band(i_band), l_initial                              &
!                 Fluxes calculated
          , flux_direct(1, 0, map_channel(i_band))                      &
          , flux_down(1, 0, map_channel(i_band))                        &
          , flux_up(1, 0, map_channel(i_band))                          &
!                 Radiances
          , i_direct, radiance(1, 1, 1, map_channel(i_band))            &
!                 Rate of photolysis
          , photolysis(1, 1, map_channel(i_band))                       &
!                 Flags for clear-sky calculations
          , l_clear_band, i_solver_clear                                &
!                 Clear-sky fluxes
          , flux_direct_clear(1, 0, map_channel(i_band))                &
          , flux_down_clear(1, 0, map_channel(i_band))                  &
          , flux_up_clear(1, 0, map_channel(i_band))                    &
!                 Tiled Surface Fluxes
          , flux_up_tile(1, 1, map_channel(i_band))                     &
          , flux_up_blue_tile(1, 1, map_channel(i_band))                &
!                 Special Surface Fluxes
          , l_blue_flux_surf, weight_blue(i_band)                       &
          , flux_direct_blue_surf                                       &
          , flux_down_blue_surf, flux_up_blue_surf                      &
!                 Dimensions of arrays
          , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column        &
          , nd_flux_profile, nd_radiance_profile, nd_j_profile          &
          , nd_band, nd_species                                         &
          , nd_esft_term, nd_scale_variable                             &
          , nd_cloud_type, nd_region, nd_overlap_coeff                  &
          , nd_max_order, nd_sph_coeff                                  &
          , nd_brdf_basis_fnc, nd_brdf_trunc                            &
          , nd_viewing_level, nd_direction                              &
          , nd_source_coeff, nd_point_tile, nd_tile                     &
          )

      ELSE IF (i_gas_overlap(i_band) == ip_overlap_random) THEN

! DEPENDS ON: solve_band_random_overlap
        CALL solve_band_random_overlap(ierr                             &
!                 Atmospheric properties
          , n_profile, n_layer, i_top, p, t, d_mass                     &
!                 Angular integration
          , i_angular_integration, i_2stream                            &
          , n_order_phase, l_rescale, n_order_gauss                     &
          , ms_min, ms_max, i_truncation, ls_local_trunc                &
          , accuracy_adaptive, euler_factor                             &
          , i_sph_algorithm, i_sph_mode                                 &
!                 Precalculated angular arrays
          , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                    &
!                 Treatment of scattering
          , i_scatter_method(i_band)                                    &
!                 Options for solver
          , i_solver                                                    &
!                 Gaseous properties
          , i_band, n_gas                                               &
          , index_absorb, i_band_esft, i_scale_esft, i_scale_fnc        &
          , k_esft, w_esft, scale_vector                                &
          , p_reference, t_reference                                    &
          , gas_mix_ratio, gas_frac_rescaled                            &
          , l_doppler, doppler_correction                               &
!                 Spectral region
          , isolir                                                      &
!                 Solar properties
          , zen_0, solar_irrad_band                                     &
!                 Infra-red properties
          , planck_flux_band(1, 0)                                      &
          , planck_flux_band(1, n_layer)                                &
          , diff_planck_band                                            &
          , l_ir_source_quad, diff_planck_band_2                        &
!                 Surface properties
          , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb(1, 1, i_band)      &
          , f_brdf, brdf_sol, brdf_hemi                                 &
          , planck_flux_ground                                          &
!                 Tiling of the surface
          , l_tile, n_point_tile, n_tile, list_tile                     &
          , rho_alb_tile(1, 1, 1, i_band)                               &
          , planck_flux_tile                                            &
!                 Optical Properties
          , ss_prop                                                     &
!                 Cloudy properties
          , l_cloud, i_cloud                                            &
!                 Cloud geometry
          , n_cloud_top                                                 &
          , n_cloud_type, frac_cloud                                    &
          , n_region, k_clr, i_region_cloud, frac_region                &
          , w_free, w_cloud, cloud_overlap                              &
          , n_column_slv, list_column_slv                               &
          , i_clm_lyr_chn, i_clm_cld_typ, area_column                   &
!                 Levels for calculating radiances
          , n_viewing_level, i_rad_layer, frac_rad_layer                &
!                 Viewing Geometry
          , n_direction, direction                                      &
!                 Weighting factor for the band
          , weight_band(i_band), l_initial                              &
!                 Fluxes calculated
          , flux_direct(1, 0, map_channel(i_band))                      &
          , flux_down(1, 0, map_channel(i_band))                        &
          , flux_up(1, 0, map_channel(i_band))                          &
!                 Radiances
          , i_direct, radiance(1, 1, 1, map_channel(i_band))            &
!                 Rate of photolysis
          , photolysis(1, 1, map_channel(i_band))                       &
!                 Flags for clear-sky calculations
          , l_clear_band, i_solver_clear                                &
!                 Clear-sky fluxes
          , flux_direct_clear(1, 0, map_channel(i_band))                &
          , flux_down_clear(1, 0, map_channel(i_band))                  &
          , flux_up_clear(1, 0, map_channel(i_band))                    &
!                 Tiled Surface Fluxes
          , flux_up_tile(1, 1, map_channel(i_band))                     &
          , flux_up_blue_tile(1, 1, map_channel(i_band))                &
!                 Special Surface Fluxes
          , l_blue_flux_surf, weight_blue(i_band)                       &
          , flux_direct_blue_surf                                       &
          , flux_down_blue_surf, flux_up_blue_surf                      &
!                 Dimensions of arrays
          , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column        &
          , nd_flux_profile, nd_radiance_profile, nd_j_profile          &
          , nd_band, nd_species                                         &
          , nd_esft_term, nd_scale_variable                             &
          , nd_cloud_type, nd_region, nd_overlap_coeff                  &
          , nd_max_order, nd_sph_coeff                                  &
          , nd_brdf_basis_fnc, nd_brdf_trunc, nd_viewing_level          &
          , nd_direction, nd_source_coeff, nd_point_tile, nd_tile       &
          )

      ELSE IF (i_gas_overlap(i_band) == ip_overlap_k_eqv) THEN

! DEPENDS ON: solve_band_k_eqv
        CALL solve_band_k_eqv(ierr                                      &
!                 Atmospheric properties
          , n_profile, n_layer, i_top, p, t, d_mass                     &
!                 Angular integration
          , i_angular_integration, i_2stream                            &
          , n_order_phase, l_rescale, n_order_gauss                     &
          , ms_min, ms_max, i_truncation, ls_local_trunc                &
          , accuracy_adaptive, euler_factor                             &
          , i_sph_algorithm, i_sph_mode                                 &
!                 Precalculated angular arrays
          , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                    &
!                 Treatment of scattering
          , i_scatter_method(i_band)                                    &
!                 Options for solver
          , i_solver                                                    &
!                 Gaseous properties
          , i_band, n_gas                                               &
          , index_absorb, i_band_esft, i_scale_esft, i_scale_fnc        &
          , k_esft, w_esft, scale_vector                                &
          , p_reference, t_reference, l_mod_k_flux                      &
          , gas_mix_ratio, gas_frac_rescaled                            &
          , l_doppler, doppler_correction                               &
!                 Spectral region
          , isolir                                                      &
!                 Solar properties
          , zen_0, solar_irrad_band                                     &
!                 Infra-red properties
          , planck_flux_band                                            &
          , diff_planck_band                                            &
          , l_ir_source_quad, diff_planck_band_2                        &
!                 Surface properties
          , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb(1, 1, i_band)      &
          , f_brdf, brdf_sol, brdf_hemi                                 &
          , diffuse_alb_basis                                           &
          , planck_flux_ground                                          &
!                 Tiling of the surface
          , l_tile, n_point_tile, n_tile, list_tile                     &
          , rho_alb_tile(1, 1, 1, i_band)                               &
          , planck_flux_tile                                            &
!                 Optical Properties
          , ss_prop                                                     &
!                 Cloudy properties
          , l_cloud, i_cloud                                            &
!                 Cloud geometry
          , n_cloud_top                                                 &
          , n_cloud_type, frac_cloud                                    &
          , n_region, k_clr, i_region_cloud, frac_region                &
          , w_free, w_cloud, cloud_overlap                              &
          , n_column_slv, list_column_slv                               &
          , i_clm_lyr_chn, i_clm_cld_typ, area_column                   &
!                 Additional variables required for mcica
          , l_cloud_cmp, n_condensed, n_cloud_profile, i_cloud_profile  &
          , i_cloud_type, nd_cloud_component                            &
!                 Levels for calculating radiances
          , n_viewing_level, i_rad_layer, frac_rad_layer                &
!                 Viewing Geometry
          , n_direction, direction                                      &
!                 Weighting factor for the band
          , weight_band(i_band), l_initial                              &
!                 Fluxes calculated
          , flux_direct(1, 0, map_channel(i_band))                      &
          , flux_down(1, 0, map_channel(i_band))                        &
          , flux_up(1, 0, map_channel(i_band))                          &
!                 Radiances
          , i_direct, radiance(1, 1, 1, map_channel(i_band))            &
!                 Rate of photolysis
          , photolysis(1, 1, map_channel(i_band))                       &
!                 Flags for clear-sky calculations
          , l_clear_band, i_solver_clear                                &
!                 Clear-sky fluxes calculated
          , flux_direct_clear(1, 0, map_channel(i_band))                &
          , flux_down_clear(1, 0, map_channel(i_band))                  &
          , flux_up_clear(1, 0, map_channel(i_band))                    &
!                 Tiled Surface Fluxes
          , flux_up_tile(1, 1, map_channel(i_band))                     &
          , flux_up_blue_tile(1, 1, map_channel(i_band))                &
!                 Special Surface Fluxes
          , l_blue_flux_surf, weight_blue(i_band)                       &
          , flux_direct_blue_surf                                       &
          , flux_down_blue_surf, flux_up_blue_surf                      &
!                 Dimensions of arrays
          , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column        &
          , nd_flux_profile, nd_radiance_profile, nd_j_profile          &
          , nd_band, nd_species                                         &
          , nd_esft_term, nd_scale_variable                             &
          , nd_cloud_type, nd_region, nd_overlap_coeff                  &
          , nd_max_order, nd_sph_coeff                                  &
          , nd_brdf_basis_fnc, nd_brdf_trunc                            &
          , nd_viewing_level, nd_direction                              &
          , nd_source_coeff, nd_point_tile, nd_tile                     &
          )

      ELSE IF (i_gas_overlap(i_band) == ip_overlap_mix_ses2) THEN

! DEPENDS ON: solve_band_ses
        CALL solve_band_ses(ierr                                        &
!                 Atmospheric properties
          , n_profile, n_layer, d_mass                                  &
!                 Angular integration
          , i_angular_integration, i_2stream                            &
          , n_order_phase, l_rescale, n_order_gauss                     &
          , ms_min, ms_max, i_truncation, ls_local_trunc                &
          , accuracy_adaptive, euler_factor                             &
          , i_sph_algorithm, i_sph_mode                                 &
!                 Precalculated angular arrays
          , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                    &
!                 Treatment of scattering
          , i_scatter_method(i_band)                                    &
!                 Options for solver
          , i_solver                                                    &
!                 Gaseous properties
          , i_band, n_gas, index_absorb                                 &
          , i_band_esft_ses, k_esft_layer, w_esft_ses                   &
          , k_mix_gas_layer, n_mix_gas(i_band)                          &
          , gas_mix_ratio, gas_mix_amt                                  &
!                 Continuum absorption
          , k_contm_layer, l_continuum                                  &
          , n_continuum, amount_continuum                               &
!                 Spectral region
          , isolir                                                      &
!                 Solar properties
          , zen_0, solar_irrad_band_ses, l_solar_tail_flux              &
!                 Infra-red properties
          , planck_flux_band(1, 0)                                      &
          , planck_flux_band(1, n_layer)                                &
          , diff_planck_band                                            &
          , l_ir_source_quad, diff_planck_band_2                        &
!                 Surface properties
          , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb(1, 1, i_band)      &
          , f_brdf, brdf_sol, brdf_hemi                                 &
          , planck_flux_ground                                          &
!                 Tiling of the surface
          , l_tile, n_point_tile, n_tile, list_tile                     &
          , rho_alb_tile(1, 1, 1, i_band)                               &
          , planck_flux_tile                                            &
!                 Optical Properties
          , ss_prop                                                     &
!                 Cloudy properties
          , l_cloud, i_cloud                                            &
!                 Cloud geometry
          , n_cloud_top                                                 &
          , n_cloud_type, frac_cloud                                    &
          , n_region, k_clr, i_region_cloud, frac_region                &
          , w_free, w_cloud, cloud_overlap                              &
          , n_column_slv, list_column_slv                               &
          , i_clm_lyr_chn, i_clm_cld_typ, area_column                   &
!                 Additional variables required for mcica
          , l_cloud_cmp, n_condensed, n_cloud_profile, i_cloud_profile  &
          , i_cloud_type, nd_cloud_component                            &
!                 Levels for calculating radiances
          , n_viewing_level, i_rad_layer, frac_rad_layer                &
!                 Viewing Geometry
          , n_direction, direction                                      &
!                 Weighting factor for the band
          , weight_band(i_band), l_initial                              &
!                 Fluxes calculated
          , flux_direct(1, 0, map_channel(i_band))                      &
          , flux_down(1, 0, map_channel(i_band))                        &
          , flux_up(1, 0, map_channel(i_band))                          &
!                 Radiances
          , i_direct, radiance(1, 1, 1, map_channel(i_band))            &
!                 Rate of photolysis
          , photolysis(1, 1, map_channel(i_band))                       &
!                 Flags for clear-sky calculations
          , l_clear_band, i_solver_clear                                &
!                 Clear-sky fluxes
          , flux_direct_clear(1, 0, map_channel(i_band))                &
          , flux_down_clear(1, 0, map_channel(i_band))                  &
          , flux_up_clear(1, 0, map_channel(i_band))                    &
!                 Tiled Surface Fluxes
          , flux_up_tile(1, 1, map_channel(i_band))                     &
          , flux_up_blue_tile(1, 1, map_channel(i_band))                &
!                 Special Surface Fluxes
          , l_blue_flux_surf, weight_blue(i_band)                       &
          , flux_direct_blue_surf                                       &
          , flux_down_blue_surf, flux_up_blue_surf                      &
!                 Dimensions of arrays
          , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column        &
          , nd_flux_profile, nd_radiance_profile, nd_j_profile          &
          , nd_band, nd_species                                         &
          , nd_esft_term, nd_continuum                                  &
          , nd_cloud_type, nd_region, nd_overlap_coeff                  &
          , nd_max_order, nd_sph_coeff                                  &
          , nd_brdf_basis_fnc, nd_brdf_trunc                            &
          , nd_viewing_level, nd_direction                              &
          , nd_source_coeff, nd_point_tile, nd_tile                     &
          )

      ELSE
        cmessage = '*** Error: An appropriate gaseous overlap ' //      &
          'has not been specified, even though gaseous ' //             &
          'absorption is to be included.'
        ierr=i_err_fatal
        CALL ereport(RoutineName, ierr, cmessage)
      END IF
    END IF

    IF ((isolir == ip_solar).AND.(i_band == i_first_band)) THEN

!     Calculate UV-Fluxes
!     (Note: as presently coded these diagnostics will only work if
!      the weight is non-zero for band 1 only. A more complete
!      treatment would require passing WEIGHT_UV through to the
!      routine AUGMENT_RADIANCE as done for WEIGHT_BLUE.)

      IF (l_uvflux_direct) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            uv_flux_direct(l,i,1)=                                      &
                    weight_uv(i_band)*flux_direct(l,i,1)
          END DO
        END DO
      END IF
      IF (l_uvflux_up) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            uv_flux_up(l,i,1)=                                          &
                    weight_uv(i_band)*flux_up(l,i,1)
          END DO
        END DO
      END IF
      IF (l_uvflux_down) THEN
        DO i=0, n_layer
          DO l=1, n_profile
            uv_flux_down(l,i,1)=                                        &
                    weight_uv(i_band)*flux_down(l,i,1)
          END DO
        END DO
      END IF
      IF (l_surf_uv) THEN
        DO l=1, n_profile
          flux_down_uv_surf(l)=                                         &
            weight_uv(i_band)*flux_down(l,n_layer,1)
        END DO
      END IF
      IF (l_surf_uv_clr) THEN
        DO l=1, n_profile
          flux_down_clr_uv_surf(l)=                                     &
            weight_uv(i_band)*flux_down_clear(l,n_layer,1)
        END DO
      END IF

    END IF


!   Spectral diagnostics for clouds:

    IF (l_cloud_extinction) THEN

!     Increment the arrays of diagnostics. The extinction
!     calculated in this band is a mean value weighted with the
!     fractions of individual types of cloud (which sum to 1).
!     Here it is weighted with the clear-sky direct solar
!     flux in the band at the top of the current
!     layer and the total amount of cloud in the grid-box.
!     This definition has the advantage of convenience, but there
!     appears to be no optimal definition of an average extinction.

      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO ll=1, n_cloud_profile(i)
           l=i_cloud_profile(ll, i)
           cloud_weight_extinction(l, i)                                &
              =cloud_weight_extinction(l, i) + w_cloud(l, i)*           &
              (flux_direct_clear(l, i-1,map_channel(i_band))            &
              -flux_direct_clear_band(l, i-1))
           cloud_extinction(l, i)                                       &
              =cloud_extinction(l, i) + w_cloud(l, i)*                  &
              (flux_direct_clear(l, i-1,map_channel(i_band))            &
              -flux_direct_clear_band(l, i-1))                          &
              *cloud_extinction_band(l, i)
        END DO
      END DO

    END IF

    IF (l_cloud_absorptivity) THEN

!     Increment the arrays of diagnostics. The absorptivity
!     calculated in this band is a mean value weighted with the
!     fractions of individual types of cloud (which sum to 1).
!     Here it is weighted with the modulus of the clear_sky
!     differential flux in the band at the top of the current
!     layer and the total amount of cloud in the grid-box, as the
!     diagnostic is a measure of the effect of introducing an
!     infinitesimal layer of layer at the top of the current
!     layer into a clear atmosphere on the upward flux at the
!     top of the cloud.

      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO ll=1, n_cloud_profile(i)
           l=i_cloud_profile(ll, i)
           cloud_weight_absorptivity(l, i)                              &
              =cloud_weight_absorptivity(l, i) + w_cloud(l, i)*         &
              ABS(flux_up_clear(l, i-1,map_channel(i_band))             &
              -flux_up_clear_band(l, i-1))
           cloud_absorptivity(l, i)                                     &
              =cloud_absorptivity(l, i) + w_cloud(l, i)*                &
              ABS(flux_up_clear(l, i-1,map_channel(i_band))             &
              -flux_up_clear_band(l, i-1))                              &
              *cloud_absorptivity_band(l, i)
        END DO
      END DO

    END IF
    IF (l_ls_cloud_extinction) THEN

!     Increment the arrays of diagnostics. The extinction
!     calculated in this band is a mean value weighted with the
!     fractions of individual types of cloud (which sum to 1).
!     Here it is weighted with the clear-sky direct solar
!     flux in the band at the top of the current
!     layer and the total amount of cloud in the grid-box.
!     This definition has the advantage of convenience, but there
!     appears to be no optimal definition of an average extinction.

      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO ll=1, n_cloud_profile(i)
           l=i_cloud_profile(ll, i)
           ls_cloud_weight_extinction(l, i)                             &
              =ls_cloud_weight_extinction(l, i)                         &
              +w_cloud(l, i)*(frac_cloud(l,i,1)+frac_cloud(l,i,2))      &
              *(flux_direct_clear(l, i-1,map_channel(i_band))           &
               -flux_direct_clear_band(l, i-1))
           ls_cloud_extinction(l, i)                                    &
              =ls_cloud_extinction(l, i)                                &
              +w_cloud(l, i)*(frac_cloud(l,i,1)+frac_cloud(l,i,2))      &
              *(flux_direct_clear(l, i-1,map_channel(i_band))           &
               -flux_direct_clear_band(l, i-1))                         &
              *ls_cloud_extinction_band(l, i)
        END DO
      END DO

    END IF

    IF (l_ls_cloud_absorptivity) THEN

!     Increment the arrays of diagnostics. The absorptivity
!     calculated in this band is a mean value weighted with the
!     fractions of individual types of cloud (which sum to 1).
!     Here it is weighted with the modulus of the clear_sky
!     differential flux in the band at the top of the current
!     layer and the total amount of cloud in the grid-box as the
!     diagnostic is a measure of the effect of introducing an
!     infinitesimal layer of layer at the top of the current
!     layer into a clear atmosphere on the upward flux at the
!     top of the cloud.

      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO ll=1, n_cloud_profile(i)
           l=i_cloud_profile(ll, i)
           ls_cloud_weight_absorptivity(l, i)                           &
              =ls_cloud_weight_absorptivity(l, i)                       &
              +w_cloud(l, i)*(frac_cloud(l,i,1)+frac_cloud(l,i,2))      &
              *ABS(flux_up_clear(l, i-1,map_channel(i_band))            &
              -flux_up_clear_band(l, i-1))
           ls_cloud_absorptivity(l, i)                                  &
              =ls_cloud_absorptivity(l, i)                              &
              +w_cloud(l, i)*(frac_cloud(l,i,1)+frac_cloud(l,i,2))      &
              *ABS(flux_up_clear(l, i-1,map_channel(i_band))            &
              -flux_up_clear_band(l, i-1))                              &
              *ls_cloud_absorptivity_band(l, i)
        END DO
      END DO

    END IF
    IF (l_cnv_cloud_extinction) THEN

!     Increment the arrays of diagnostics. The extinction
!     calculated in this band is a mean value weighted with the
!     fractions of individual types of cloud (which sum to 1).
!     Here it is weighted with the clear-sky direct solar
!     flux in the band at the top of the current
!     layer and the total amount of cloud in the grid-box.
!     This definition has the advantage of convenience, but there
!     appears to be no optimal definition of an average extinction.

      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO ll=1, n_cloud_profile(i)
           l=i_cloud_profile(ll, i)
           cnv_cloud_weight_extinction(l, i)                            &
              =cnv_cloud_weight_extinction(l, i)                        &
              +w_cloud(l, i)*(frac_cloud(l,i,3)+frac_cloud(l,i,4))      &
              *(flux_direct_clear(l, i-1,map_channel(i_band))           &
               -flux_direct_clear_band(l, i-1))
           cnv_cloud_extinction(l, i)                                   &
              =cnv_cloud_extinction(l, i)                               &
              +w_cloud(l, i)*(frac_cloud(l,i,3)+frac_cloud(l,i,4))      &
              *(flux_direct_clear(l, i-1,map_channel(i_band))           &
               -flux_direct_clear_band(l, i-1))                         &
              *cnv_cloud_extinction_band(l, i)
        END DO
      END DO

    END IF

    IF (l_cnv_cloud_absorptivity) THEN

!     Increment the arrays of diagnostics. The absorptivity
!     calculated in this band is a mean value weighted with the
!     fractions of individual types of cloud (which sum to 1).
!     Here it is weighted with the modulus of the clear_sky
!     differential flux in the band at the top of the current
!     layer and the total amount of cloud in the grid-box as the
!     diagnostic is a measure of the effect of introducing an
!     infinitesimal layer of layer at the top of the current
!     layer into a clear atmosphere on the upward flux at the
!     top of the cloud.

      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO ll=1, n_cloud_profile(i)
           l=i_cloud_profile(ll, i)
           cnv_cloud_weight_absorptivity(l, i)                          &
              =cnv_cloud_weight_absorptivity(l, i)                      &
              +w_cloud(l, i)*(frac_cloud(l,i,3)+frac_cloud(l,i,4))      &
              *ABS(flux_up_clear(l, i-1,map_channel(i_band))            &
              -flux_up_clear_band(l, i-1))
           cnv_cloud_absorptivity(l, i)                                 &
              =cnv_cloud_absorptivity(l, i)                             &
              +w_cloud(l, i)*(frac_cloud(l,i,3)+frac_cloud(l,i,4))      &
              *ABS(flux_up_clear(l, i-1,map_channel(i_band))            &
              -flux_up_clear_band(l, i-1))                              &
              *cnv_cloud_absorptivity_band(l, i)
        END DO
      END DO

    END IF

!   Deallocate the single scattering propeties (in reverse order).
    DEALLOCATE(ss_prop%forward_scatter_no_cloud)
    DEALLOCATE(ss_prop%phase_fnc_no_cloud)

    DEALLOCATE(ss_prop%forward_solar_cloud_comp)
    DEALLOCATE(ss_prop%phase_fnc_solar_cloud_comp)
    DEALLOCATE(ss_prop%forward_scatter_cloud_comp)
    DEALLOCATE(ss_prop%phase_fnc_cloud_comp)
    DEALLOCATE(ss_prop%k_ext_scat_cloud_comp)
    DEALLOCATE(ss_prop%k_ext_tot_cloud_comp)

    DEALLOCATE(ss_prop%omega)
    DEALLOCATE(ss_prop%tau)

    DEALLOCATE(ss_prop%forward_solar)
    DEALLOCATE(ss_prop%phase_fnc_solar)
    DEALLOCATE(ss_prop%forward_scatter)
    DEALLOCATE(ss_prop%phase_fnc)
    DEALLOCATE(ss_prop%k_ext_scat)
    DEALLOCATE(ss_prop%k_grey_tot)

    DEALLOCATE(ss_prop%omega_clr)
    DEALLOCATE(ss_prop%tau_clr)

    DEALLOCATE(ss_prop%forward_solar_clr)
    DEALLOCATE(ss_prop%phase_fnc_solar_clr)
    DEALLOCATE(ss_prop%forward_scatter_clr)
    DEALLOCATE(ss_prop%phase_fnc_clr)
    DEALLOCATE(ss_prop%k_ext_scat_clr)
    DEALLOCATE(ss_prop%k_grey_tot_clr)


!   Make any adjustments to fluxes and radiances to convert
!   to actual values. This is done inside the loop over bands
!   to allow for division of the output fluxes between
!   separate diagnostic bands.
    IF (isolir == ip_infra_red) THEN
! DEPENDS ON: adjust_ir_radiance
      CALL adjust_ir_radiance(n_profile, n_layer, n_viewing_level       &
        , n_direction, i_angular_integration, i_sph_mode                &
        , planck_flux_band, planck_radiance_band                        &
        , flux_down(1, 0, map_channel(i_band))                          &
        , flux_up(1, 0, map_channel(i_band))                            &
        , radiance(1, 1, 1, map_channel(i_band))                        &
        , l_clear_band                                                  &
        , flux_down_clear(1, 0, map_channel(i_band))                    &
        , flux_up_clear(1, 0, map_channel(i_band))                      &
        , nd_2sg_profile, nd_flux_profile, nd_radiance_profile          &
        , nd_layer, nd_direction, nd_viewing_level                      &
        )
    END IF


  END DO ! i_band


  IF (lhook) CALL dr_hook('RADIANCE_CALC',zhook_out,zhook_handle)

END SUBROUTINE radiance_calc
