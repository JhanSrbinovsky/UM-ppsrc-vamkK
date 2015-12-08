! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate the fluxes within the band with one gas.
!
! Method:
!   Monochromatic calculations are performed for each ESFT term
!   and the results are summed.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiance Core
!
!- ---------------------------------------------------------------------
SUBROUTINE solve_band_one_gas(ierr                                      &
!                 Atmospheric Column
    , n_profile, n_layer, i_top, p, t, d_mass                           &
!                 Angular Integration
    , i_angular_integration, i_2stream                                  &
    , n_order_phase, l_rescale, n_order_gauss                           &
    , ms_min, ms_max, i_truncation, ls_local_trunc                      &
    , accuracy_adaptive, euler_factor                                   &
    , i_sph_algorithm, i_sph_mode                                       &
!                 Precalculated angular arrays
    , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                          &
!                 Treatment of Scattering
    , i_scatter_method                                                  &
!                 Options for Solver
    , i_solver                                                          &
!                 Gaseous Properties
    , i_band, i_gas                                                     &
    , i_band_esft, i_scale_esft, i_scale_fnc                            &
    , k_esft, w_esft, scale_vector                                      &
    , p_reference, t_reference                                          &
    , gas_mix_ratio, gas_frac_rescaled                                  &
    , l_doppler, doppler_correction                                     &
!                 Spectral Region
    , isolir                                                            &
!                 Solar Properties
    , zen_0, solar_irrad                                                &
!                 Infra-red Properties
    , planck_flux_top, planck_flux_bottom                               &
    , diff_planck_band                                                  &
    , l_ir_source_quad, diff_planck_band_2                              &
!                 Surface Properties
    , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                          &
    , f_brdf, brdf_sol, brdf_hemi                                       &
    , planck_flux_ground                                                &
!                   Tiling of the surface
    , l_tile, n_point_tile, n_tile, list_tile, rho_alb_tile             &
    , planck_flux_tile                                                  &
!                 Optical Properties
    , ss_prop                                                           &
!                 Cloudy Properties
    , l_cloud, i_cloud                                                  &
!                 Cloud Geometry
    , n_cloud_top                                                       &
    , n_cloud_type, frac_cloud                                          &
    , n_region, k_clr, i_region_cloud, frac_region                      &
    , w_free, w_cloud, cloud_overlap                                    &
    , n_column_slv, list_column_slv                                     &
    , i_clm_lyr_chn, i_clm_cld_typ, area_column                         &
!                   Levels for calculating radiances
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
!                   Viewing Geometry
    , n_direction, direction                                            &
!                   Weighting factor for the band
    , weight_band, l_initial                                            &
!                 Calculated Fluxes
    , flux_direct, flux_down, flux_up                                   &
!                   Calculated radiances
    , i_direct, radiance                                                &
!                   Calculated rate of photolysis
    , photolysis                                                        &
!                 Flags for Clear-sky Fluxes
    , l_clear, i_solver_clear                                           &
!                 Clear-sky Fluxes
    , flux_direct_clear, flux_down_clear, flux_up_clear                 &
!                   Tiled Surface Fluxes
    , flux_up_tile, flux_up_blue_tile                                   &
!                   Special Surface Fluxes
    , l_blue_flux_surf, weight_blue                                     &
    , flux_direct_blue_surf                                             &
    , flux_down_blue_surf, flux_up_blue_surf                            &
!                 Dimensions of Arrays
    , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column              &
    , nd_flux_profile, nd_radiance_profile, nd_j_profile                &
    , nd_band, nd_species                                               &
    , nd_esft_term, nd_scale_variable                                   &
    , nd_cloud_type, nd_region, nd_overlap_coeff                        &
    , nd_max_order, nd_sph_coeff                                        &
    , nd_brdf_basis_fnc, nd_brdf_trunc, nd_viewing_level                &
    , nd_direction, nd_source_coeff                                     &
    , nd_point_tile, nd_tile                                            &
    )


  USE realtype_rd, ONLY: RealK
  USE def_ss_prop
  USE rad_pcf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_layer                                                          &
!       Size allocated for layers
    , nd_layer_clr                                                      &
!       Size allocated for totally clear layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_flux_profile                                                   &
!       Size allocated for profiles in arrays of fluxes
    , nd_radiance_profile                                               &
!       Size allocated for profiles in arrays of radiances
    , nd_j_profile                                                      &
!       Size allocated for profiles in arrays of mean radiances
    , nd_column                                                         &
!       Size allocated for sub-columns per point
    , nd_band                                                           &
!       Size allocated for bands
    , nd_species                                                        &
!       Size allocated for species
    , nd_esft_term                                                      &
!       Size allocated for ESFT variables
    , nd_scale_variable                                                 &
!       Size allocated for scaling variables
    , nd_cloud_type                                                     &
!       Size allocated for cloud types
    , nd_region                                                         &
!       Size allocated for cloudy regions
    , nd_overlap_coeff                                                  &
!       Size allocated for cloudy overlap coefficients
    , nd_max_order                                                      &
!       Size allocated for orders of spherical harmonics
    , nd_sph_coeff                                                      &
!       Size allocated for coefficients of spherical harmonics
    , nd_brdf_basis_fnc                                                 &
!       Size allowed for BRDF basis functions
    , nd_brdf_trunc                                                     &
!       Size allowed for orders of BRDFs
    , nd_viewing_level                                                  &
!       Size allocated for levels where radiances are calculated
    , nd_direction                                                      &
!       Size allocated for viewing directions
    , nd_source_coeff                                                   &
!       Size allocated for source coefficients
    , nd_point_tile                                                     &
!       Size allocated for points where the surface is tiled
    , nd_tile
!       Size allocated for surface tiles


! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag

!                 Atmospheric column
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer                                                           &
!       Number of layers
    , i_top
!       Top of vertical grid
  REAL (RealK), INTENT(IN) ::                                           &
      p(nd_profile, nd_layer)                                           &
!       Pressure
    , t(nd_profile, nd_layer)                                           &
!       Temperature
    , d_mass(nd_profile, nd_layer)
!       Mass thickness of each layer

!                 Angular integration
  INTEGER, INTENT(IN) ::                                                &
      i_angular_integration                                             &
!       Angular integration scheme
    , i_2stream                                                         &
!       Two-stream scheme
    , n_order_phase                                                     &
!       Maximum order of terms in the phase function used in
!       the direct calculation of spherical harmonics
    , n_order_gauss                                                     &
!       Order of gaussian integration
    , ms_min                                                            &
!       Lowest azimuthal order used
    , ms_max                                                            &
!       Highest azimuthal order used
    , i_truncation                                                      &
!       Type of truncation used
    , ia_sph_mm(0: nd_max_order)                                        &
!       Address of spherical coefficient of (m, m) for each m
    , ls_local_trunc(0: nd_max_order)                                   &
!       Orders of truncation at each azimuthal order
    , i_sph_mode                                                        &
!       Mode in which the spherical harmonic code is used
    , i_sph_algorithm
!       Algorithm used for spherical harmonic calculation
  LOGICAL, INTENT(IN) ::                                                &
      l_rescale
!       Rescale optical properties
  REAL (RealK) ::                                                       &
      cg_coeff(nd_sph_coeff)                                            &
!       Clebsch-Gordan coefficients
    , uplm_zero(nd_sph_coeff)                                           &
!       Values of spherical harmonics at polar angles pi/2
    , uplm_sol(nd_radiance_profile, nd_sph_coeff)                       &
!       Values of spherical harmonics in the solar direction
    , accuracy_adaptive                                                 &
!       Accuracy for adaptive truncation
    , euler_factor
!       Factor applied to the last term of an alternating series
  REAL (RealK), INTENT(IN) ::                                           &
      weight_band
!       Weighting factor for the current band
  LOGICAL, INTENT(INOUT) ::                                             &
      l_initial
!       Flag to initialize diagnostics

!                 Treatment of scattering
  INTEGER, INTENT(IN) ::                                                &
      i_scatter_method
!       Method of treating scattering

!                 Options for solver
  INTEGER, INTENT(IN) ::                                                &
      i_solver
!       Solver used

!                 Gaseous properties
  INTEGER, INTENT(IN) ::                                                &
      i_band                                                            &
!       Band being considered
    , i_gas                                                             &
!       Gas being considered
    , i_band_esft(nd_band, nd_species)                                  &
!       Number of terms in band
    , i_scale_esft(nd_band, nd_species)                                 &
!       Type of ESFT scaling
    , i_scale_fnc(nd_band, nd_species)
!       Type of scaling function
  LOGICAL, INTENT(IN) ::                                                &
      l_doppler(nd_species)
!       Doppler broadening included
  REAL (RealK), INTENT(IN) ::                                           &
      k_esft(nd_esft_term, nd_band, nd_species)                         &
!       Exponential ESFT terms
    , w_esft(nd_esft_term, nd_band, nd_species)                         &
!       Weights for ESFT
    , scale_vector(nd_scale_variable, nd_esft_term, nd_band             &
        , nd_species)                                                   &
!       Absorber scaling parameters
    , p_reference(nd_species, nd_band)                                  &
!       Reference scaling pressure
    , t_reference(nd_species, nd_band)                                  &
!       Reference scaling temperature
    , gas_mix_ratio(nd_profile, nd_layer, nd_species)                   &
!       Gas mass mixing ratios
    , doppler_correction(nd_species)
!       Doppler broadening terms
  REAL (RealK), INTENT(INOUT) ::                                        &
      gas_frac_rescaled(nd_profile, nd_layer, nd_species)
!       Rescaled gas mass fractions

!                 Spectral region
  INTEGER, INTENT(IN) ::                                                &
      isolir
!       Spectral region

!                 Solar properties
  REAL (RealK), INTENT(IN) ::                                           &
       zen_0(nd_profile)                                                &
!       Secant (two-stream) or cosine (spherical harmonics)
!       of the solar zenith angle
    , solar_irrad(nd_profile)
!       Incident solar irradiance in band

!                 Infra-red properties
  LOGICAL, INTENT(IN) ::                                                &
      l_ir_source_quad
!       Use a quadratic source function
  REAL (RealK), INTENT(IN) ::                                           &
      planck_flux_top(nd_profile)                                       &
!       Planckian flux at the top of the layer
    , planck_flux_bottom(nd_profile)                                    &
!       Planckian source at the bottom of the layer
    , diff_planck_band(nd_profile, nd_layer)                            &
!       Thermal source function
    , diff_planck_band_2(nd_profile, nd_layer)
!       Twice second difference of Planckian in band

!                 Surface properties
  REAL (RealK), INTENT(IN) ::                                           &
      planck_flux_ground(nd_profile)
!       Thermal source function at ground
  INTEGER, INTENT(IN) ::                                                &
      ls_brdf_trunc                                                     &
!       Order of truncation of BRDFs
    , n_brdf_basis_fnc
!       Number of BRDF basis functions
  REAL (RealK), INTENT(IN) ::                                           &
      rho_alb(nd_profile, nd_brdf_basis_fnc)                            &
!       Weights of the basis functions
    , f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                      &
        , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)                         &
!       Array of BRDF basis terms
    , brdf_sol(nd_profile, nd_brdf_basis_fnc, nd_direction)             &
!       The BRDF evaluated for scattering from the solar
!       beam into the viewing direction
    , brdf_hemi(nd_profile, nd_brdf_basis_fnc, nd_direction)
!       The BRDF evaluated for scattering from isotropic
!       radiation into the viewing direction

! Variables related to tiling of the surface
  LOGICAL, INTENT(IN) ::                                                &
      l_tile
!       Logical to allow invoke options
  INTEGER, INTENT(IN) ::                                                &
      n_point_tile                                                      &
!       Number of points to tile
    , n_tile                                                            &
!       Number of tiles used
    , list_tile(nd_point_tile)
!       List of points with surface tiling
  REAL (RealK), INTENT(IN) ::                                           &
      rho_alb_tile(nd_point_tile, nd_brdf_basis_fnc, nd_tile)           &
!       Weights for the basis functions of the BRDFs
!       at the tiled points
    , planck_flux_tile(nd_point_tile, nd_tile)
!       Local Planckian fluxes on surface tiles

!                 Optical properties
  TYPE(STR_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere

!                 Cloudy properties
  LOGICAL, INTENT(IN) ::                                                &
      l_cloud
!       Clouds required
  INTEGER, INTENT(IN) ::                                                &
      i_cloud
!       Cloud scheme used

!                 Cloud geometry
  INTEGER, INTENT(IN) ::                                                &
      n_cloud_top                                                       &
!       Top cloudy layer
    , n_cloud_type                                                      &
!       Number of types of clouds
    , n_region                                                          &
!       Number of cloudy regions
    , k_clr                                                             &
!       Index of clear-sky region
    , i_region_cloud(nd_cloud_type)
!       Regions in which types of clouds fall

! Cloud geometry
  INTEGER, INTENT(IN) ::                                                &
      n_column_slv(nd_profile)                                          &
!       Number of columns to be solved in each profile
    , list_column_slv(nd_profile, nd_column)                            &
!       List of columns requiring an actual solution
    , i_clm_lyr_chn(nd_profile, nd_column)                              &
!       Layer in the current column to change
    , i_clm_cld_typ(nd_profile, nd_column)
!       Type of cloud to introduce in the changed layer
  REAL (RealK), INTENT(IN) ::                                           &
      w_free(nd_profile, id_ct: nd_layer)                               &
!       Clear-sky fraction
    , w_cloud(nd_profile, id_ct: nd_layer)                              &
!       Cloudy fraction
    , frac_cloud(nd_profile, id_ct: nd_layer, nd_cloud_type)            &
!       Fractions of types of clouds
    , cloud_overlap(nd_profile, id_ct-1: nd_layer, nd_overlap_coeff)    &
!       Coefficients for transfer for energy at interfaces
    , area_column(nd_profile, nd_column)                                &
!       Areas of columns
    , frac_region(nd_profile, id_ct: nd_layer, nd_region)
!       Fractions of total cloud occupied by each region


  INTEGER, INTENT(IN) ::                                                &
      n_viewing_level                                                   &
!       Number of levels where radiances are calculated
    , i_rad_layer(nd_viewing_level)
!       Layers in which radiances are calculated
  REAL (RealK), INTENT(IN) ::                                           &
      frac_rad_layer(nd_viewing_level)
!       Fractions below the tops of the layers

!                   Viewing Geometry
  INTEGER, INTENT(IN) ::                                                &
      n_direction
!       Number of viewing directions
  REAL (RealK), INTENT(IN) ::                                           &
      direction(nd_radiance_profile, nd_direction, 2)
!       Viewing directions

!                 Calculated fluxes
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct(nd_flux_profile, 0: nd_layer)                         &
!       Direct flux
    , flux_down(nd_flux_profile, 0: nd_layer)                           &
!       Total downward flux
    , flux_up(nd_flux_profile, 0: nd_layer)
!       Upward flux

!                   Calculated radiances
  REAL (RealK), INTENT(INOUT) ::                                        &
      i_direct(nd_radiance_profile, 0: nd_layer)                        &
!       Direct solar irradiance on levels
    , radiance(nd_radiance_profile,  nd_viewing_level                   &
        , nd_direction)
!       Radiances

!                   Calculated mean radiances
  REAL (RealK), INTENT(INOUT) ::                                        &
      photolysis(nd_j_profile,  nd_viewing_level)
!       Mean rate of photolysis

!                 Flags for clear-sky calculations
  LOGICAL, INTENT(IN) ::                                                &
      l_clear
!       Calculate net clear-sky properties
  INTEGER, INTENT(IN) ::                                                &
      i_solver_clear
!       Clear solver used

!                 Clear-sky fluxes calculated
  REAL (RealK), INTENT(OUT) ::                                          &
      flux_direct_clear(nd_profile, 0: nd_layer)                        &
!       Clear-sky direct flux
    , flux_down_clear(nd_profile, 0: nd_layer)                          &
!       Clear-sky total downward flux
    , flux_up_clear(nd_profile, 0: nd_layer)                            &
!       Clear-sky upward flux
    , flux_up_tile(nd_point_tile, nd_tile)                              &
!       Upward fluxes at tiled surface points
    , flux_up_blue_tile(nd_point_tile, nd_tile)
!       Upward blue fluxes at tiled surface points

!                 Special Diagnostics:
  LOGICAL, INTENT(IN) ::                                                &
      l_blue_flux_surf
!       Flag to calculate the blue flux at the surface
  REAL (RealK), INTENT(IN) ::                                           &
      weight_blue
!       Weights for blue fluxes in this band
  REAL (RealK), INTENT(INOUT) ::                                        &
      flux_direct_blue_surf(nd_flux_profile)                            &
!       Direct blue flux at the surface
    , flux_down_blue_surf(nd_flux_profile)                              &
!       Total downward blue flux at the surface
    , flux_up_blue_surf(nd_flux_profile)
!       Upward blue flux at the surface



! Local variables.
  INTEGER                                                               &
      l
!       Loop variable
  INTEGER                                                               &
      i_gas_pointer(nd_species)                                         &
!       Pointer array for monochromatic ESFTs
    , iex
!       Index of ESFT term
  REAL (RealK) ::                                                       &
      k_esft_mono(nd_species)                                           &
!       ESFT monochromatic exponents
    , k_gas_abs(nd_profile, nd_layer)                                   &
!       Gaseous absorptive extinction
    , d_planck_flux_surface(nd_profile)                                 &
!       Ground source function
    , flux_inc_direct(nd_profile)                                       &
!       Incident direct flux
    , flux_inc_down(nd_profile)                                         &
!       Incident downward flux
    , dummy_ke(nd_profile, nd_layer)
!       Dummy array (not used)

! Monochromatic incrementing radiances:
  REAL (RealK) ::                                                       &
      flux_direct_part(nd_profile, 0: nd_layer)                         &
!       Partial direct flux
    , flux_total_part(nd_profile, 2*nd_layer+2)                         &
!       Partial total flux
    , flux_direct_clear_part(nd_profile, 0: nd_layer)                   &
!       Partial clear-sky direct flux
    , flux_total_clear_part(nd_profile, 2*nd_layer+2)
!       Partial clear-sky total flux
  REAL (RealK) ::                                                       &
      i_direct_part(nd_radiance_profile, 0: nd_layer)                   &
!       Partial solar irradiances
    , radiance_part(nd_radiance_profile, nd_viewing_level               &
        , nd_direction)
!       Partial radiances
  REAL (RealK) ::                                                       &
      photolysis_part(nd_j_profile, nd_viewing_level)
!       Partial rate of photolysis
  REAL (RealK) ::                                                       &
      weight_incr                                                       &
!       Weight applied to increments
    , weight_blue_incr
!       Weight applied to blue increments

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('SOLVE_BAND_ONE_GAS',zhook_in,zhook_handle)

! The ESFT terms for the first gas in the band alone are used.
  i_gas_pointer(1)=i_gas
  DO iex=1, i_band_esft(i_band, i_gas)

!   Rescale for each ESFT term if that is required.
    IF (i_scale_esft(i_band, i_gas) == ip_scale_term) THEN
! DEPENDS ON: scale_absorb
      CALL scale_absorb(ierr, n_profile, n_layer                        &
        , gas_mix_ratio(1, 1, i_gas), p, t                              &
        , i_top                                                         &
        , gas_frac_rescaled(1, 1, i_gas)                                &
        , i_scale_fnc(i_band, i_gas)                                    &
        , p_reference(i_gas, i_band)                                    &
        , t_reference(i_gas, i_band)                                    &
        , scale_vector(1, iex, i_band, i_gas)                           &
        , iex, i_band                                                   &
        , l_doppler(i_gas), doppler_correction(i_gas)                   &
        , nd_profile, nd_layer                                          &
        , nd_scale_variable                                             &
        )
    END IF

!   Set the appropriate boundary terms for the total
!   upward and downward fluxes.

    IF ( (i_angular_integration == ip_two_stream).OR.                   &
         (i_angular_integration == ip_ir_gauss) ) THEN

      IF (isolir == ip_solar) THEN
!       Solar region.
        DO l=1, n_profile
          d_planck_flux_surface(l)=0.0e+00_RealK
          flux_inc_down(l)=solar_irrad(l)/zen_0(l)
          flux_inc_direct(l)=solar_irrad(l)/zen_0(l)
        END DO
      ELSE IF (isolir == ip_infra_red) THEN
!       Infra-red region.
        DO l=1, n_profile
          flux_inc_direct(l)=0.0e+00_RealK
          flux_inc_down(l)=-planck_flux_top(l)
          d_planck_flux_surface(l)                                      &
            =(1.0e+00_RealK-rho_alb(l, ip_surf_alb_diff))               &
            *(planck_flux_ground(l)-planck_flux_bottom(l))
        END DO
      END IF

    ELSE IF (i_angular_integration == ip_spherical_harmonic) THEN

      IF (isolir == ip_solar) THEN
        DO l=1, n_profile
          i_direct_part(l, 0)=solar_irrad(l)
          flux_inc_down(l)=0.0e+00_RealK
        END DO
      ELSE
        DO l=1, n_profile
          flux_inc_down(l)=-planck_flux_top(l)
          d_planck_flux_surface(l)                                      &
            =planck_flux_ground(l)-planck_flux_bottom(l)
        END DO
      END IF

    END IF

!   Assign the monochromatic absorption coefficient.
    k_esft_mono(i_gas)=k_esft(iex, i_band, i_gas)

! DEPENDS ON: gas_optical_properties
    CALL gas_optical_properties(n_profile, n_layer                      &
      , 1, i_gas_pointer, k_esft_mono                                   &
      , gas_frac_rescaled                                               &
      , k_gas_abs                                                       &
      , nd_profile, nd_layer, nd_species                                &
      )


! DEPENDS ON: monochromatic_radiance
    CALL monochromatic_radiance(ierr                                    &
!                 Atmospheric properties
      , n_profile, n_layer, d_mass                                      &
!                 Angular integration
      , i_angular_integration, i_2stream                                &
      , l_rescale, n_order_gauss                                        &
      , n_order_phase, ms_min, ms_max, i_truncation, ls_local_trunc     &
      , accuracy_adaptive, euler_factor                                 &
      , i_sph_algorithm, i_sph_mode                                     &
!                   Precalculated angular arrays
      , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                        &
!                 Treatment of scattering
      , i_scatter_method                                                &
!                 Options for solver
      , i_solver                                                        &
!                 Gaseous propreties
      , k_gas_abs                                                       &
!                 Options for equivalent extinction
      , .FALSE., dummy_ke                                               &
!                 Spectral region
      , isolir                                                          &
!                 Infra-red properties
      , diff_planck_band                                                &
      , l_ir_source_quad, diff_planck_band_2                            &
!                 Conditions at TOA
      , zen_0, flux_inc_direct, flux_inc_down                           &
      , i_direct_part                                                   &
!                 Surface properties
      , d_planck_flux_surface                                           &
      , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                        &
      , f_brdf, brdf_sol, brdf_hemi                                     &
!                 Optical properties
      , ss_prop                                                         &
!                 Cloudy properties
      , l_cloud, i_cloud                                                &
!                 Cloud geometry
      , n_cloud_top                                                     &
      , n_cloud_type, frac_cloud                                        &
      , n_region, k_clr, i_region_cloud, frac_region                    &
      , w_free, w_cloud, cloud_overlap                                  &
      , n_column_slv, list_column_slv                                   &
      , i_clm_lyr_chn, i_clm_cld_typ, area_column                       &
!                   Levels for calculating radiances
      , n_viewing_level, i_rad_layer, frac_rad_layer                    &
!                   Viewing Geometry
      , n_direction, direction                                          &
!                 Calculated flxues
      , flux_direct_part, flux_total_part                               &
!                   Calculated Radiances
      , radiance_part                                                   &
!                   Calculated rates of photolysis
      , photolysis_part                                                 &
!                 Flags for clear-sky calculations
      , l_clear, i_solver_clear                                         &
!                 Clear-sky fluxes calculated
      , flux_direct_clear_part, flux_total_clear_part                   &
!                 Dimensions of arrays
      , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column            &
      , nd_flux_profile, nd_radiance_profile, nd_j_profile              &
      , nd_cloud_type, nd_region, nd_overlap_coeff                      &
      , nd_max_order, nd_sph_coeff                                      &
      , nd_brdf_basis_fnc, nd_brdf_trunc, nd_viewing_level              &
      , nd_direction, nd_source_coeff                                   &
      )


!   Increment the radiances within the band. Each increment
!   represents a single k-term within a band weighted with
!   its own weighting factor, hence for each increment the
!   weighting is the product of these two factors: similarly
!   for the blue flux.
    weight_incr=weight_band*w_esft(iex, i_band,  i_gas)
    IF (l_blue_flux_surf)                                               &
      weight_blue_incr=weight_blue*w_esft(iex, i_band,  i_gas)
! DEPENDS ON: augment_radiance
    CALL augment_radiance(n_profile, n_layer                            &
      , i_angular_integration, i_sph_mode                               &
      , n_viewing_level, n_direction                                    &
      , isolir, l_clear, l_initial, weight_incr                         &
      , l_blue_flux_surf, weight_blue_incr                              &
!                   Actual radiances
      , flux_direct, flux_down, flux_up                                 &
      , flux_direct_blue_surf                                           &
      , flux_down_blue_surf, flux_up_blue_surf                          &
      , i_direct, radiance, photolysis                                  &
      , flux_direct_clear, flux_down_clear, flux_up_clear               &
!                   Increments to radiances
      , flux_direct_part, flux_total_part                               &
      , i_direct_part, radiance_part, photolysis_part                   &
      , flux_direct_clear_part, flux_total_clear_part                   &
!                   Dimensions
      , nd_flux_profile, nd_radiance_profile, nd_j_profile              &
      , nd_layer, nd_viewing_level, nd_direction                        &
      )

!   Add in the increments from surface tiles
    IF (l_tile) THEN
! DEPENDS ON: augment_tiled_radiance
      CALL augment_tiled_radiance(ierr                                  &
        , n_point_tile, n_tile, list_tile                               &
        , i_angular_integration, isolir, l_initial                      &
        , weight_incr, l_blue_flux_surf, weight_blue_incr               &
!                   Surface characteristics
        , rho_alb_tile                                                  &
!                   Actual radiances
        , flux_up_tile, flux_up_blue_tile                               &
!                   Increments to radiances
        , flux_direct_part(1, n_layer)                                  &
        , flux_total_part(1, 2*n_layer+2)                               &
        , planck_flux_tile, planck_flux_bottom                          &
!                   Dimensions
        , nd_flux_profile, nd_point_tile, nd_tile                       &
        , nd_brdf_basis_fnc                                             &
        )
    END IF

!   After the first call to these routines quantities should be
!   incremented rather than initialized, until the flag is reset.
    l_initial=.FALSE.

  END DO


  IF (lhook) CALL dr_hook('SOLVE_BAND_ONE_GAS',zhook_out,zhook_handle)

END SUBROUTINE solve_band_one_gas
