! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to sample the cloud field for a single k-term
!
! Method:
!   Returns the fluxes for a single k-term. A number of cloudy
!   sub-columns are sampled (dependent on the configuration of McICA)
!   and a single clear-sky calculation is performed. The results are
!   meaned making use of the fraction of cloudy sub-columns.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiance Core
!
!- ---------------------------------------------------------------------
SUBROUTINE mcica_sample(ierr                                            &
!                 Atmospheric Propertries
    , n_profile, n_layer, d_mass                                        &
!                 Angular Integration
    , i_angular_integration, i_2stream                                  &
    , l_rescale, n_order_gauss                                          &
    , n_order_phase, ms_min, ms_max, i_truncation, ls_local_trunc       &
    , accuracy_adaptive, euler_factor, i_sph_algorithm                  &
    , i_sph_mode                                                        &
!                 Precalculated angular arrays
    , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                          &
!                 Treatment of Scattering
    , i_scatter_method                                                  &
!                 Options for Solver
    , i_solver                                                          &
!                 Gaseous Properties
    , k_gas_abs                                                         &
!                 Options for Equivalent Extinction
    , l_scale_solar, adjust_solar_ke                                    &
!                 Spectral Region
    , isolir                                                            &
!                 Infra-red Properties
    , diff_planck, l_ir_source_quad, diff_planck_2                      &
!                 Conditions at TOA
    , zen_0, flux_inc_direct, flux_inc_down                             &
    , i_direct                                                          &
!                 Surface Properties
    , d_planck_flux_surface                                             &
    , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                          &
    , f_brdf, brdf_sol, brdf_hemi                                       &
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
!                 Additional variables required for mcica
    , l_cloud_cmp, n_condensed, n_cloud_profile, i_cloud_profile        &
    , i_cloud_type, nd_cloud_component, iex, i_band                     &
!                 Levels for calculating radiances
    , n_viewing_level, i_rad_layer, frac_rad_layer                      &
!                 Viewing geometry
    , n_direction, direction                                            &
!                 Calculated Fluxes
    , flux_direct, flux_total                                           &
!                 Calculated Radiances
    , radiance                                                          &
!                 Calculated mean radiances
    , j_radiance                                                        &
!                 Flags for Clear-sky Calculation
    , l_clear, i_solver_clear                                           &
!                 Clear-sky Fluxes Calculated
    , flux_direct_clear, flux_total_clear                               &
!                 Dimensions of Arrays
    , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column              &
    , nd_flux_profile, nd_radiance_profile, nd_j_profile                &
    , nd_cloud_type, nd_region, nd_overlap_coeff                        &
    , nd_max_order, nd_sph_coeff                                        &
    , nd_brdf_basis_fnc, nd_brdf_trunc, nd_viewing_level                &
    , nd_direction, nd_source_coeff                                     &
    )


  USE realtype_rd, ONLY: RealK
  USE def_ss_prop
  USE rad_pcf, ONLY: ip_solar
  USE mcica_mod, ONLY: first_subcol_k, index_subcol, subcol_need,       &
    subcol_reorder, subcol_k, l_avg_phase_fnc, c_sub, frac_cloudy
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Maximum number of profiles
    , nd_layer                                                          &
!       Maximum number of layers
    , nd_layer_clr                                                      &
!       Maximum number of completely clear layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_flux_profile                                                   &
!       Maximum number of profiles in arrays of fluxes
    , nd_radiance_profile                                               &
!       Maximum number of profiles in arrays of radiances
    , nd_j_profile                                                      &
!       Maximum number of profiles in arrays of mean radiances
    , nd_column                                                         &
!       Number of columns per point
    , nd_cloud_type                                                     &
!       Maximum number of types of cloud
    , nd_region                                                         &
!       Maximum number of cloudy regions
    , nd_overlap_coeff                                                  &
!       Maximum number of overlap coefficients
    , nd_max_order                                                      &
!       Maximum order of spherical harmonics used
    , nd_sph_coeff                                                      &
!       Allocated size for spherical coefficients
    , nd_brdf_basis_fnc                                                 &
!       Size allowed for BRDF basis functions
    , nd_brdf_trunc                                                     &
!       Size allowed for orders of BRDFs
    , nd_viewing_level                                                  &
!       Allocated size for levels where radiances are calculated
    , nd_direction                                                      &
!       Allocated size for viewing directions
    , nd_source_coeff
!       Size allocated for source coefficients


! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag

!                 Atmospheric properties
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer
!       Number of layers
  REAL (RealK), INTENT(IN) ::                                           &
      d_mass(nd_profile, nd_layer)
!       Mass thickness of each layer

!                 Angular integration
  INTEGER, INTENT(IN) ::                                                &
      i_angular_integration                                             &
!       Angular integration scheme
    , i_2stream                                                         &
!       Two-stream scheme
    , n_order_gauss                                                     &
!       Order of Gaussian integration
    , n_order_phase                                                     &
!       Highest order retained in the phase function
    , i_truncation                                                      &
!       Type of spherical truncation adopted
    , ms_min                                                            &
!       Lowest azimuthal order calculated
    , ms_max                                                            &
!       Highest azimuthal order calculated
    , ia_sph_mm(0: nd_max_order)                                        &
!       Address of spherical coefficient for (m, m) for each m
    , ls_local_trunc(0: nd_max_order)                                   &
!       Orders of truncation at each azimuthal order
    , i_sph_mode                                                        &
!       Mode in which teh spherical harmonic solver is being used
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
    , uplm_sol(nd_radiance_profile, nd_sph_coeff)
!       Values of spherical harmonics in the solar direction
  REAL (RealK), INTENT(IN) ::                                           &
      accuracy_adaptive                                                 &
!       Accuracy for adaptive truncation
    , euler_factor
!       Factor applied to the last term of an alternating series

!                 Treatment of scattering
  INTEGER, INTENT(IN) ::                                                &
      i_scatter_method

!                 Options for solver
  INTEGER, INTENT(IN) ::                                                &
      i_solver
!       Solver used

!                 Gaseous properties
  REAL (RealK), INTENT(IN) ::                                           &
      k_gas_abs(nd_profile, nd_layer)
!       Gaseous absorptive extinctions

!                 Variables for equivalent extinction
  LOGICAL, INTENT(IN) ::                                                &
      l_scale_solar
!       Apply scaling to solar flux
  REAL (RealK), INTENT(IN) ::                                           &
      adjust_solar_ke(nd_profile, nd_layer)
!       Adjustment of solar beam with equivalent extinction

!                 Spectral region
  INTEGER, INTENT(IN) ::                                                &
      isolir
!       Visible or IR

!                 Infra-red properties
  LOGICAL, INTENT(IN) ::                                                &
      l_ir_source_quad
!       Flag for quadratic IR-source
  REAL (RealK), INTENT(IN) ::                                           &
      diff_planck(nd_profile, nd_layer)                                 &
!       DIfferences in the Planckian function across layers
    , diff_planck_2(nd_profile, nd_layer)
!       Twice the second differences of Planckian source function

!                 Conditions at TOA
  REAL (RealK), INTENT(IN) ::                                           &
      zen_0(nd_profile)                                                 &
!       Secants (two-stream) or cosines (spherical harmonics)
!       of the solar zenith angles
    , flux_inc_direct(nd_profile)                                       &
!       Incident direct flux
    , flux_inc_down(nd_profile)
!       Incident downward flux
  REAL (RealK), INTENT(INOUT) ::                                        &
      i_direct(nd_radiance_profile, 0: nd_layer)
!       Direct radiance (the first row contains the incident
!       solar radiance: the other rows are calculated)

!                 Surface properties
  REAL (RealK), INTENT(IN) ::                                           &
      d_planck_flux_surface(nd_profile)
!       Differential Planckian flux from the surface
  INTEGER, INTENT(IN) ::                                                &
      ls_brdf_trunc                                                     &
!       Order of trunation of BRDFs
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
!       Topmost cloudy layer
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
      w_cloud(nd_profile, id_ct: nd_layer)                              &
!       Cloudy fraction
    , frac_cloud(nd_profile, id_ct: nd_layer, nd_cloud_type)            &
!       Fractions of different types of cloud
    , w_free(nd_profile, id_ct: nd_layer)                               &
!       Clear-sky fraction
    , cloud_overlap(nd_profile, id_ct-1: nd_layer                       &
        , nd_overlap_coeff)                                             &
!       Coefficients for energy transfer at interfaces
    , area_column(nd_profile, nd_column)                                &
!       Areas of columns
    , frac_region(nd_profile, id_ct: nd_layer, nd_region)
!       Fractions of total cloud occupied by each region


!                 Levels where radiance are calculated
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

!                 Calculated Fluxes
  REAL (RealK), INTENT(INOUT) ::                                        &
      flux_direct(nd_flux_profile, 0: nd_layer)                         &
!       Direct flux
    , flux_total(nd_flux_profile, 2*nd_layer+2)
!       Total flux

!                 Calculated radiances
  REAL (RealK), INTENT(INOUT) ::                                        &
      radiance(nd_radiance_profile, nd_viewing_level, nd_direction)
!       Radiances
!                 Calculated mean radiances
  REAL (RealK), INTENT(INOUT) ::                                        &
      j_radiance(nd_j_profile, nd_viewing_level)
!       Mean radiances

!                 Flags for clear-sky calculations
  LOGICAL, INTENT(IN) ::                                                &
      l_clear
!       Calculate clear-sky properties
  INTEGER, INTENT(IN) ::                                                &
      i_solver_clear
!       Clear solver used

!                 Clear-sky fluxes calculated
  REAL (RealK), INTENT(INOUT) ::                                        &
      flux_direct_clear(nd_flux_profile, 0: nd_layer)                   &
!       Clear-sky direct flux
    , flux_total_clear(nd_flux_profile, 2*nd_layer+2)
!       Clear-sky total flux

! Variables required for McICA
  INTEGER, INTENT(IN) ::                                                &
      n_condensed                                                       &
!       number of condensed components
    , n_cloud_profile(id_ct: nd_layer)                                  &
!       number of cloudy profiles in each layer
    , i_cloud_profile(nd_profile, id_ct: nd_layer)                      &
!       profiles containing clouds
    , nd_cloud_component                                                &
!       size allocated for components of clouds
    , i_cloud_type(nd_cloud_component)                                  &
!       types of cloud to which each component contributes
    , iex                                                               &
!       Index of ESFT term
    , i_band
!       Band being considered

  LOGICAL, INTENT(IN) ::                                                &
      l_cloud_cmp(nd_cloud_component)
!       flags to activate cloudy components



! Local variables.
  INTEGER :: i, j, l, ll, ls, k, m
!       Loop variables

  REAL (RealK) ::                                                       &
      i_direct_subcol(nd_radiance_profile, 0: nd_layer)                 &
!       Partial solar irradiances
    , flux_direct_subcol(nd_flux_profile, 0: nd_layer)                  &
!       Partial direct flux
    , flux_total_subcol(nd_flux_profile, 2*nd_layer+2)                  &
!       Partial total flux
    , radiance_subcol(nd_radiance_profile, nd_viewing_level             &
        , nd_direction)                                                 &
!       Partial radiances
    , photolysis_subcol(nd_j_profile, nd_viewing_level)
!       Partial rate of photolysis

  REAL (RealK) ::                                                       &
      subcol_k_inv
!       1.0/subcol_k

  LOGICAL ::                                                            &
      l_clear_calc
!       flag for calculating clear sky conrtribution to McICA

  REAL (RealK), PARAMETER :: eps=EPSILON(ss_prop%k_ext_scat)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('MCICA_SAMPLE',zhook_in,zhook_handle)

  DO m=first_subcol_k(i_band,iex), first_subcol_k(i_band,iex+1)-1

    index_subcol=MOD(m,subcol_need)
    IF (index_subcol == 0)  index_subcol=subcol_need
    index_subcol=subcol_reorder(index_subcol)
    subcol_k_inv=1.0e+00_RealK/subcol_k(i_band,iex)


    IF (m==first_subcol_k(i_band,iex))THEN
      l_clear_calc=.TRUE.
    ELSE
      l_clear_calc=.FALSE.
    END IF


!   Loop over the condensed components, calculating their optical
!   properties by simply scaling the values previously calculated
!   for the average cloud condensate amounts.
    IF (l_avg_phase_fnc) THEN
      DO k=1, n_condensed
        IF (l_cloud_cmp(k)) THEN
          DO i=n_cloud_top, n_layer
!CDIR NODEP
            DO ll=1, n_cloud_profile(i)
              l=i_cloud_profile(ll, i)
              ss_prop%k_grey_tot(l, i, i_cloud_type(k))                 &
                =ss_prop%k_grey_tot(l, i, 0)                            &
                +ss_prop%k_ext_tot_cloud_comp(l, i, k)                  &
                *c_sub(l, i, index_subcol, i_cloud_type(k))
              ss_prop%k_ext_scat(l, i, i_cloud_type(k))                 &
                =ss_prop%k_ext_scat(l, i, 0)                            &
                +ss_prop%k_ext_scat_cloud_comp(l, i, k)                 &
                *c_sub(l, i, index_subcol, i_cloud_type(k))
            END DO
          END DO
        END IF
      END DO
    ELSE
      IF (l_rescale) THEN
        DO k=1, n_condensed
          IF (l_cloud_cmp(k)) THEN
            DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                l=i_cloud_profile(ll, i)
                ss_prop%k_grey_tot(l, i, i_cloud_type(k))               &
                  =ss_prop%k_grey_tot(l, i, 0)                          &
                  +ss_prop%k_ext_tot_cloud_comp(l, i, k)                &
                  *c_sub(l, i, index_subcol, i_cloud_type(k))
                ss_prop%k_ext_scat(l, i, i_cloud_type(k))               &
                  =ss_prop%k_ext_scat(l, i, 0)                          &
                  +ss_prop%k_ext_scat_cloud_comp(l, i, k)               &
                  *c_sub(l, i, index_subcol, i_cloud_type(k))
                ss_prop%forward_scatter(l, i, i_cloud_type(k))          &
                  =ss_prop%forward_scatter_no_cloud(l, i)               &
                  +ss_prop%forward_scatter_cloud_comp(l, i, k)          &
                  *c_sub(l, i, index_subcol, i_cloud_type(k))
                DO ls=1, n_order_phase
                  ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))          &
                    =(ss_prop%phase_fnc_no_cloud(l, i, ls)              &
                    +ss_prop%phase_fnc_cloud_comp(l, i, ls, k)          &
                    *c_sub(l, i, index_subcol, i_cloud_type(k))         &
                    -ss_prop%forward_scatter(l, i, i_cloud_type(k)))    &
                    /MAX(ss_prop%k_ext_scat(l, i, i_cloud_type(k))      &
                    -ss_prop%forward_scatter(l, i, i_cloud_type(k))     &
                    ,eps)
                END DO
                ss_prop%forward_scatter(l, i, i_cloud_type(k))          &
                  =ss_prop%forward_scatter(l, i, i_cloud_type(k))       &
                  /MAX(ss_prop%k_ext_scat(l, i, i_cloud_type(k))        &
                  ,eps)
              END DO
            END DO
          END IF
        END DO
      ELSE
        DO k=1, n_condensed
          IF (l_cloud_cmp(k)) THEN
            DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                l=i_cloud_profile(ll, i)
                ss_prop%k_grey_tot(l, i, i_cloud_type(k))               &
                  =ss_prop%k_grey_tot(l, i, 0)                          &
                  +ss_prop%k_ext_tot_cloud_comp(l, i, k)                &
                  *c_sub(l, i, index_subcol, i_cloud_type(k))
                ss_prop%k_ext_scat(l, i, i_cloud_type(k))               &
                  =ss_prop%k_ext_scat(l, i, 0)                          &
                  +ss_prop%k_ext_scat_cloud_comp(l, i, k)               &
                  *c_sub(l, i, index_subcol, i_cloud_type(k))
                DO ls=1, n_order_phase
                  ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))          &
                    =(ss_prop%phase_fnc_no_cloud(l, i, ls)              &
                    +ss_prop%phase_fnc_cloud_comp(l, i, ls, k)          &
                    *c_sub(l, i, index_subcol, i_cloud_type(k)))        &
                    /MAX(ss_prop%k_ext_scat(l, i, i_cloud_type(k))      &
                    ,eps)
                END DO
              END DO
            END DO
          END IF
        END DO
      END IF
    END IF


! DEPENDS ON: monochromatic_radiance
    CALL monochromatic_radiance(ierr                                    &
!             atmospheric properties
      , n_profile, n_layer, d_mass                                      &
!             angular integration
      , i_angular_integration, i_2stream                                &
      , l_rescale, n_order_gauss                                        &
      , n_order_phase, ms_min, ms_max, i_truncation                     &
      , ls_local_trunc                                                  &
      , accuracy_adaptive, euler_factor                                 &
      , i_sph_algorithm, i_sph_mode                                     &
!               precalculated angular arrays
      , ia_sph_mm, cg_coeff, uplm_zero, uplm_sol                        &
!             treatment of scattering
      , i_scatter_method                                                &
!             options for solver
      , i_solver                                                        &
!             gaseous propreties
      , k_gas_abs                                                       &
!             options for equivalent extinction
      , l_scale_solar, adjust_solar_ke                                  &
!             spectral region
      , isolir                                                          &
!             infra-red properties
      , diff_planck, l_ir_source_quad, diff_planck_2                    &
!             conditions at toa
      , zen_0, flux_inc_direct, flux_inc_down                           &
      , i_direct_subcol                                                 &
!             surface properties
      , d_planck_flux_surface                                           &
      , ls_brdf_trunc, n_brdf_basis_fnc, rho_alb                        &
      , f_brdf, brdf_sol, brdf_hemi                                     &
!             optical properties
      , ss_prop                                                         &
!             cloudy properties
      , l_cloud, i_cloud                                                &
!             cloud geometry
      , n_cloud_top                                                     &
      , n_cloud_type, frac_cloud                                        &
      , n_region, k_clr, i_region_cloud, frac_region                    &
      , w_free, w_cloud, cloud_overlap                                  &
      , n_column_slv, list_column_slv                                   &
      , i_clm_lyr_chn, i_clm_cld_typ, area_column                       &
!               levels for the calculation of radiances
      , n_viewing_level, i_rad_layer, frac_rad_layer                    &
!               viewing geometry
      , n_direction, direction                                          &
!               calculated fluxes
      , flux_direct_subcol, flux_total_subcol                           &
!               calculated radiances
      , radiance_subcol                                                 &
!               calculated rate of photolysis
      , photolysis_subcol                                               &
!             flags for clear-sky calculations
      , l_clear_calc, i_solver_clear                                    &
!             clear-sky fluxes calculated
      , flux_direct_clear, flux_total_clear                             &
!             planckian function
      , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_column            &
      , nd_flux_profile, nd_radiance_profile, nd_j_profile              &
      , nd_cloud_type, nd_region, nd_overlap_coeff                      &
      , nd_max_order, nd_sph_coeff                                      &
      , nd_brdf_basis_fnc, nd_brdf_trunc, nd_viewing_level              &
      , nd_direction, nd_source_coeff                                   &
      )

    IF (m == first_subcol_k(i_band,iex)) THEN
      IF (isolir == ip_solar) THEN
        DO k=0,n_layer
          DO j=1,n_profile
            flux_direct(j,k)=flux_direct_subcol(j,k)
          END DO
        END DO
      END IF

      DO k=1,2*n_layer+2
        DO j=1,n_profile
          flux_total(j,k)=flux_total_subcol(j,k)
        END DO
      END DO
    ELSE
      IF (isolir == ip_solar) THEN
        DO k=0,n_layer
          DO j=1,n_profile
            flux_direct(j,k)=flux_direct(j,k)+flux_direct_subcol(j,k)
          END DO
        END DO
      END IF

      DO k=1,2*n_layer+2
        DO j=1,n_profile
          flux_total(j,k)=flux_total(j,k)+flux_total_subcol(j,k)
        END DO
      END DO
    END IF

  END DO


  IF (isolir == ip_solar) THEN
    DO k=0,n_layer
      DO j=1,n_profile
        flux_direct(j,k)=flux_direct(j,k)*subcol_k_inv
        flux_direct(j,k)=frac_cloudy(j)*flux_direct(j,k)                &
          + ((1.0_RealK-frac_cloudy(j))*flux_direct_clear(j,k))
      END DO
    END DO
  END IF

  DO k=1,2*n_layer+2
    DO j=1,n_profile
      flux_total(j,k)=flux_total(j,k)*subcol_k_inv
      flux_total(j,k)=(frac_cloudy(j)*flux_total(j,k))                  &
        + ((1.0_RealK-frac_cloudy(j))*flux_total_clear(j,k))
    END DO
  END DO

  IF (lhook) CALL dr_hook('MCICA_SAMPLE',zhook_out,zhook_handle)

END SUBROUTINE mcica_sample
