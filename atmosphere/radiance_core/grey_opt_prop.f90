! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to calculate grey optical properties.
!
! Method:
!   For each activated optical process, excluding gaseous
!   absorption, increments are calculated for the total and
!   scattering extinctions, and the products of the asymmetry
!   factor and the forward scattering factor in clear and
!   cloudy regions. These increments are summed, and the grey
!   total and scattering extinctions and the asymmetry and forward
!   scattering factors are thus calculated.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiance Core
!
!- ---------------------------------------------------------------------
SUBROUTINE grey_opt_prop(ierr                                           &
    , n_profile, n_layer, p, t, density                                 &
    , n_order_phase, l_rescale, n_order_forward                         &
    , l_henyey_greenstein_pf, l_solar_phf, l_lanczos                    &
    , n_order_phase_solar, n_direction, cos_sol_view                    &
    , l_rayleigh, rayleigh_coeff                                        &
    , l_continuum, n_continuum, i_continuum_pointer, k_continuum        &
    , amount_continuum                                                  &
    , l_aerosol, n_aerosol, n_aerosol_mr, aerosol_mix_ratio             &
    , aerosol_mr_source, aerosol_mr_type_index                          &
    , i_aerosol_parametrization                                         &
    , i_humidity_pointer, humidities, delta_humidity                    &
    , mean_rel_humidity                                                 &
    , aerosol_absorption, aerosol_scattering, aerosol_phase_fnc         &
    , l_ukca_radaer, n_ukca_mode, ukca_modal_mixr                       &
    , ukca_absorption, ukca_scattering, ukca_asymmetry                  &
    , l_cloud, i_cloud, n_cloud_profile, i_cloud_profile                &
    , n_cloud_top, n_condensed, l_cloud_cmp, i_phase_cmp                &
    , i_condensed_param, condensed_n_phf, condensed_param_list          &
    , condensed_mix_ratio, condensed_dim_char                           &
    , n_cloud_type, i_cloud_type                                        &
    , ss_prop                                                           &
    , frac_cloud                                                        &
    , l_cloud_extinction, cloud_extinction                              &
    , l_cloud_absorptivity, cloud_absorptivity                          &
    , l_ls_cloud_extinction, ls_cloud_extinction                        &
    , l_ls_cloud_absorptivity, ls_cloud_absorptivity                    &
    , l_cnv_cloud_extinction, cnv_cloud_extinction                      &
    , l_cnv_cloud_absorptivity, cnv_cloud_absorptivity                  &
    , nd_profile, nd_radiance_profile, nd_layer                         &
    , nd_layer_clr, id_ct                                               &
    , nd_continuum, nd_aerosol_species, nd_aerosol_mixratio             &
    , nd_humidities, nd_cloud_parameter, nd_cloud_component             &
    , nd_cloud_type, nd_phase_term, nd_max_order, nd_direction          &
    , nd_ukca_mode                                                      &
    )


  USE realtype_rd, ONLY: RealK
  USE def_ss_prop
  USE rad_pcf
  USE mcica_mod, ONLY: l_avg_phase_fnc
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of dummy arrays
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_radiance_profile                                               &
!       Size allocated for profiles of quantities specifically
!       used in calulating radiances
    , nd_layer                                                          &
!       Size allocated for layers
    , nd_layer_clr                                                      &
!       Size allocated for completely clear layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_direction                                                      &
!       Size allocated for viewing directions
    , nd_aerosol_species                                                &
!       Size allocated for aerosols in spectral information
    , nd_aerosol_mixratio                                               &
!       Size allocated for aerosols in aerosol_mix_ratio array
    , nd_humidities                                                     &
!       Size allocated for humidities
    , nd_continuum                                                      &
!       Size allocated for continua
    , nd_phase_term                                                     &
!       Size allocated for terms in the phase function
    , nd_max_order                                                      &
!       Size allocated for the order of the calculation
    , nd_cloud_parameter                                                &
!       Size allocated for cloud parameters
    , nd_cloud_component                                                &
!       Size allocated for components of clouds
    , nd_cloud_type                                                     &
!       Size allocated for types of clouds
    , nd_ukca_mode
!       Size allocated for UKCA aerosols


! Dummy arguments.
  INTEGER, INTENT(INOUT) ::                                             &
      ierr
!       Error flag


! Basic atmospheric properties:
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
    , density(nd_profile, nd_layer)
!       Density at levels


! Optical switches:
  LOGICAL, INTENT(IN) ::                                                &
      l_rescale                                                         &
!       Delta-rescaling required
    , l_henyey_greenstein_pf                                            &
!       Flag to use a Henyey-Greenstein phase function
    , l_solar_phf                                                       &
!       Flag to use an extended phase function for solar
!       radiation
    , l_lanczos
!       Flag to use Lanczos smoothing of solar phf
  INTEGER, INTENT(IN) ::                                                &
      n_order_phase                                                     &
!       Order of terms in the phase function
    , n_order_phase_solar                                               &
!       Order of truncation of the solar beam
    , n_order_forward
!       Order used in forming the forward scattering parameter

! Directional information
  INTEGER, INTENT(IN) ::                                                &
      n_direction
!       Number of viewing directions
  REAL (RealK), INTENT(IN) ::                                           &
      cos_sol_view(nd_radiance_profile, nd_direction)
!       Cosines of the angles between the solar direction
!       and the viewing direction


! Rayleigh scattering:
  LOGICAL, INTENT(IN) ::                                                &
      l_rayleigh
!       Rayleigh scattering activated

  REAL (RealK), INTENT(IN) ::                                           &
      rayleigh_coeff
!       Rayleigh coefficient


! Continuum processes:
  LOGICAL, INTENT(IN) ::                                                &
      l_continuum
!       Continuum absorption activated

  INTEGER, INTENT(IN) ::                                                &
      n_continuum                                                       &
!       Number of continua
    , i_continuum_pointer(nd_continuum)
!       Pointers to active continua

  REAL (RealK), INTENT(IN) ::                                           &
      k_continuum(nd_continuum)                                         &
!       Continuum extinction
    , amount_continuum(nd_profile, nd_layer, nd_continuum)
!       Amounts for continua


! Properties of aerosols:
  LOGICAL, INTENT(IN) ::                                                &
      l_aerosol
!       Aerosols activated

  INTEGER, INTENT(IN) ::                                                &
      n_aerosol                                                         &
!       Number of aerosol species in spectral information
    , n_aerosol_mr                                                      &
!       Number of aerosol species in aerosol_mix_ratio array
    , aerosol_mr_type_index(nd_aerosol_mixratio)                        &
!       Index relating aerosol_mix_ratio aerosols to aerosols in
!       the spectral information
    , aerosol_mr_source(nd_aerosol_mixratio)                            &
!       Scheme/source of the aerosol data, to determine use in
!       changing radiative fluxes and use in diagnostics
    , i_aerosol_parametrization(nd_aerosol_species)                     &
!       Parametrizations of aerosols
    , i_humidity_pointer(nd_profile,  nd_layer)
!       Pointer to aerosol look-up table

  REAL (RealK), INTENT(IN) ::                                           &
      aerosol_mix_ratio(nd_profile, nd_layer, nd_aerosol_mixratio)      &
!       Number densty of aerosols
    , aerosol_absorption(nd_humidities, nd_aerosol_species)             &
!       Aerosol absorption in band for a mixing ratio of unity
    , aerosol_scattering(nd_humidities, nd_aerosol_species)             &
!       Aerosol scattering in band for a mixing ratio of unity
    , aerosol_phase_fnc(nd_humidities                                   &
        , nd_phase_term, nd_aerosol_species)                            &
!       Aerosol phase function in band
    , humidities(nd_humidities, nd_aerosol_species)                     &
!       Array of humidities
    , delta_humidity                                                    &
!       Increment in humidity
    , mean_rel_humidity(nd_profile, nd_layer)
!       Mixing ratio of water vapour


! Properties of UKCA aerosols:
  LOGICAL, INTENT(IN) ::                                                &
      l_ukca_radaer
!       UKCA aerosols activated

  INTEGER, INTENT(IN) ::                                                &
      n_ukca_mode
!       Actual number of UKCA aerosol modes

  REAL (RealK), INTENT(IN) ::                                           &
      ukca_modal_mixr(nd_profile, nd_layer, nd_ukca_mode)               &
!       Modal mass mixing ratio
    , ukca_absorption(nd_profile, nd_layer, nd_ukca_mode)               &
!       Waveband-averaged specific coefficient for absorption
    , ukca_scattering(nd_profile, nd_layer, nd_ukca_mode)               &
!       Waveband-averaged specific coefficient for scattering
    , ukca_asymmetry(nd_profile, nd_layer, nd_ukca_mode)
!       Waveband-averaged asymmetry factor


! Properties of clouds:
  LOGICAL, INTENT(IN) ::                                                &
      l_cloud
!       Clouds activated

  INTEGER, INTENT(IN) ::                                                &
      i_cloud
!       Cloud scheme used

! Geometry of clouds:
  INTEGER, INTENT(IN) ::                                                &
      n_cloud_top                                                       &
!       Topmost cloudy layer
    , n_cloud_type                                                      &
!       Number of types of clouds
    , n_cloud_profile(id_ct: nd_layer)                                  &
!       Number of cloudy profiles in each layer
    , i_cloud_profile(nd_profile, id_ct: nd_layer)                      &
!       Profiles containing clouds
    , i_cloud_type(nd_cloud_component)
!       Types of cloud to which each component contributes

! Microphysical quantities:
  INTEGER, INTENT(IN) ::                                                &
      n_condensed                                                       &
!       Number of condensed components
    , i_phase_cmp(nd_cloud_component)                                   &
!       Phases of cloudy components
    , i_condensed_param(nd_cloud_component)                             &
!       Parametrization schemes for cloudy components
    , condensed_n_phf(nd_cloud_component)
!       Number of terms in the phase function for each
!       cloudy component

  LOGICAL, INTENT(IN) ::                                                &
      l_cloud_cmp(nd_cloud_component)
!       Flags to activate cloudy components

  REAL (RealK), INTENT(IN) ::                                           &
      condensed_param_list(nd_cloud_parameter                           &
        , nd_cloud_component)                                           &
!       Coefficients in parametrization schemes
    , condensed_mix_ratio(nd_profile, id_ct: nd_layer                   &
        , nd_cloud_component)                                           &
!       Mixing ratios of cloudy components
    , condensed_dim_char(nd_profile, id_ct: nd_layer                    &
        , nd_cloud_component)
!       Characteristic dimensions of cloudy components


! Variables required for extra diagnostic calculations.
  LOGICAL, INTENT(IN) ::                                                &
      l_cloud_extinction                                                &
!       Flag for explicit calculation of extinction
    , l_cloud_absorptivity                                              &
!       Flag for explicit calculation of absorptivity
    , l_ls_cloud_extinction                                             &
!       Flag for explicit calculation of extinction
    , l_ls_cloud_absorptivity                                           &
!       Flag for explicit calculation of absorptivity
    , l_cnv_cloud_extinction                                            &
!       Flag for explicit calculation of extinction
    , l_cnv_cloud_absorptivity
!       Flag for explicit calculation of absorptivity

  REAL (RealK), INTENT(IN) ::                                           &
      frac_cloud(nd_profile, id_ct:nd_layer, nd_cloud_type)
!       Fractions of each type of cloud


! Optical properties
  TYPE(str_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere

  REAL (RealK), INTENT(OUT) ::                                          &
      cloud_extinction(nd_profile, nd_layer)                            &
!       Mean cloud extinction
    , cloud_absorptivity(nd_profile, nd_layer)                          &
!       Mean cloud absorptivity
    , ls_cloud_extinction(nd_profile, nd_layer)                         &
!       Mean cloud extinction
    , ls_cloud_absorptivity(nd_profile, nd_layer)                       &
!       Mean cloud absorptivity
    , cnv_cloud_extinction(nd_profile, nd_layer)                        &
!       Mean cloud extinction
    , cnv_cloud_absorptivity(nd_profile, nd_layer)
!       Mean cloud absorptivity



! Local variables.
  INTEGER ::                                                            &
      i_continuum                                                       &
!       Temporary continuum `index'
    , l                                                                 &
!       Loop variable
    , ll                                                                &
!       Loop variable
    , i                                                                 &
!       Loop variable
    , id                                                                &
!       Loop variable
    , j                                                                 &
!       Loop variable
    , k                                                                 &
!       Loop variable
    , ls                                                                &
!       Loop variable
    , n_index                                                           &
!       Number of indices satisfying the test
    , indx(nd_profile)
!       Indices satifying the test

! Temporary variable for the divisions
  REAL (RealK) :: tmp_inv(nd_profile)

  REAL (RealK), PARAMETER :: tiny_k=TINY(ss_prop%k_ext_scat)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('GREY_OPT_PROP',zhook_in,zhook_handle)

! If using a separate solar phase function that must be initialized.
  IF (l_solar_phf) THEN
    DO id=1, n_direction
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          ss_prop%phase_fnc_solar_clr(l, i, id)=0.0_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          ss_prop%phase_fnc_solar(l, i, id, 0)=0.0_RealK
        END DO
      END DO
    END DO
    DO i=1, n_cloud_top-1
      DO l=1, n_profile
        ss_prop%forward_solar_clr(l, i)=0.0_RealK
      END DO
    END DO
    DO i=n_cloud_top, n_layer
      DO l=1, n_profile
        ss_prop%forward_solar(l, i, 0)=0.0_RealK
      END DO
    END DO
  END IF



! Consider each optical process in turn.

! Rayleigh scattering:

  IF (l_rayleigh) THEN
!   Forward scattering is required only when delta-rescaling
!   is performed.
    IF (l_rescale) THEN
      DO i=1, n_cloud_top-1
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot_clr(l, i)=0.0_RealK
          ss_prop%k_ext_scat_clr(l, i)=rayleigh_coeff
          ss_prop%phase_fnc_clr(l, i, 1)=0.0_RealK
          ss_prop%forward_scatter_clr(l, i)=0.0_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot(l, i, 0)=0.0_RealK
          ss_prop%k_ext_scat(l, i, 0)=rayleigh_coeff
          ss_prop%phase_fnc(l, i, 1, 0)=0.0_RealK
          ss_prop%forward_scatter(l, i, 0)=0.0_RealK
        END DO
      END DO
    ELSE
      DO i=1, n_cloud_top-1
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot_clr(l, i)=0.0_RealK
          ss_prop%k_ext_scat_clr(l, i)=rayleigh_coeff
          ss_prop%phase_fnc_clr(l, i, 1)=0.0_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot(l, i, 0)=0.0_RealK
          ss_prop%k_ext_scat(l, i, 0)=rayleigh_coeff
          ss_prop%phase_fnc(l, i, 1, 0)=0.0_RealK
        END DO
      END DO
    END IF

!   Only the second Lengendre polynomial contributes.
    IF (n_order_phase >= 2) THEN
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          ss_prop%phase_fnc_clr(l, i, 2)=rayleigh_coeff*1.0e-01_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          ss_prop%phase_fnc(l, i, 2, 0)=rayleigh_coeff*1.0e-01_RealK
        END DO
      END DO
      DO ls=3, n_order_phase
        DO i=1, n_cloud_top-1
          DO l=1, n_profile
            ss_prop%phase_fnc_clr(l, i, ls)=0.0_RealK
          END DO
        END DO
        DO i=n_cloud_top, n_layer
          DO l=1, n_profile
            ss_prop%phase_fnc(l, i, ls, 0)=0.0_RealK
          END DO
        END DO
      END DO
    END IF

!   No formal rescaling is applied to the phase function for
!   Rayleigh scattering, as only g_2 is non-zero.

    IF (l_solar_phf) THEN

      DO id=1, n_direction
        DO i=1, n_cloud_top-1
          DO l=1, n_profile
            ss_prop%phase_fnc_solar_clr(l, i, id)                       &
              =ss_prop%phase_fnc_solar_clr(l, i, id)                    &
              +rayleigh_coeff                                           &
              *0.75_RealK*(1.0_RealK+cos_sol_view(l, id)**2)
          END DO
        END DO
        DO i=n_cloud_top, n_layer
          DO l=1, n_profile
            ss_prop%phase_fnc_solar(l, i, id, 0)                        &
              =ss_prop%phase_fnc_solar(l, i, id, 0)                     &
              +rayleigh_coeff                                           &
              *0.75_RealK*(1.0_RealK+cos_sol_view(l, id)**2)
          END DO
        END DO
      END DO

    END IF

  ELSE

    IF (l_rescale) THEN
      DO i=1, n_cloud_top-1
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot_clr(l, i)=0.0_RealK
          ss_prop%k_ext_scat_clr(l, i)=0.0_RealK
          ss_prop%phase_fnc_clr(l, i, 1)=0.0_RealK
          ss_prop%forward_scatter_clr(l, i)=0.0_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot(l, i, 0)=0.0_RealK
          ss_prop%k_ext_scat(l, i, 0)=0.0_RealK
          ss_prop%phase_fnc(l, i, 1, 0)=0.0_RealK
          ss_prop%forward_scatter(l, i, 0)=0.0_RealK
        END DO
      END DO
    ELSE
      DO i=1, n_cloud_top-1
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot_clr(l, i)=0.0_RealK
          ss_prop%k_ext_scat_clr(l, i)=0.0_RealK
          ss_prop%phase_fnc_clr(l, i, 1)=0.0_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
!CDIR NODEP
        DO l=1, n_profile
          ss_prop%k_grey_tot(l, i, 0)=0.0_RealK
          ss_prop%k_ext_scat(l, i, 0)=0.0_RealK
          ss_prop%phase_fnc(l, i, 1, 0)=0.0_RealK
        END DO
      END DO
    END IF
    DO ls=2, n_order_phase
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          ss_prop%phase_fnc_clr(l, i, ls)=0.0_RealK
        END DO
      END DO
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          ss_prop%phase_fnc(l, i, ls, 0)=0.0_RealK
        END DO
      END DO
    END DO

  END IF

  IF (l_aerosol) THEN
!   Include the effects of aerosol.
!   Above clouds.
! DEPENDS ON: opt_prop_aerosol
    CALL opt_prop_aerosol(ierr                                          &
      , n_profile, 1, n_cloud_top-1                                     &
      , n_order_phase, l_rescale, n_order_forward                       &
      , l_henyey_greenstein_pf                                          &
      , n_aerosol, n_aerosol_mr, aerosol_mix_ratio                      &
      , aerosol_mr_source, aerosol_mr_type_index                        &
      , i_aerosol_parametrization                                       &
      , i_humidity_pointer, humidities, delta_humidity                  &
      , mean_rel_humidity                                               &
      , aerosol_absorption, aerosol_scattering, aerosol_phase_fnc       &
      , l_solar_phf, n_order_phase_solar, n_direction, cos_sol_view     &
      , ss_prop%k_grey_tot_clr, ss_prop%k_ext_scat_clr                  &
      , ss_prop%phase_fnc_clr, ss_prop%forward_scatter_clr              &
      , ss_prop%forward_solar_clr, ss_prop%phase_fnc_solar_clr          &
      , nd_profile, nd_radiance_profile, nd_layer                       &
      , 1, nd_layer_clr                                                 &
      , nd_aerosol_species, nd_aerosol_mixratio, nd_humidities          &
      , nd_phase_term, nd_max_order, nd_direction                       &
      )
!   Within clouds:
    CALL opt_prop_aerosol(ierr                                          &
      , n_profile, n_cloud_top, n_layer                                 &
      , n_order_phase, l_rescale, n_order_forward                       &
      , l_henyey_greenstein_pf                                          &
      , n_aerosol, n_aerosol_mr, aerosol_mix_ratio                      &
      , aerosol_mr_source, aerosol_mr_type_index                        &
      , i_aerosol_parametrization                                       &
      , i_humidity_pointer, humidities, delta_humidity                  &
      , mean_rel_humidity                                               &
      , aerosol_absorption, aerosol_scattering, aerosol_phase_fnc       &
      , l_solar_phf, n_order_phase_solar, n_direction, cos_sol_view     &
      , ss_prop%k_grey_tot(1, id_ct, 0)                                 &
      , ss_prop%k_ext_scat(1, id_ct, 0)                                 &
      , ss_prop%phase_fnc(1, id_ct, 1, 0)                               &
      , ss_prop%forward_scatter(1, id_ct, 0)                            &
      , ss_prop%forward_solar(1, id_ct, 0)                              &
      , ss_prop%phase_fnc_solar(1, id_ct, 1, 0)                         &
      , nd_profile, nd_radiance_profile, nd_layer                       &
      , id_ct, nd_layer                                                 &
      , nd_aerosol_species, nd_aerosol_mixratio, nd_humidities          &
      , nd_phase_term, nd_max_order, nd_direction                       &
      )
  END IF

  IF (l_ukca_radaer) THEN
!   Include the effects of UKCA aerosols.
!   Above clouds.
! DEPENDS ON: opt_prop_ukca_aerosol
    CALL opt_prop_ukca_aerosol(ierr                                     &
      , n_profile, 1, n_cloud_top-1                                     &
      , n_order_phase, l_rescale, n_order_forward                       &
      , n_ukca_mode, ukca_modal_mixr                                    &
      , ukca_absorption, ukca_scattering, ukca_asymmetry                &
      , ss_prop%k_grey_tot_clr, ss_prop%k_ext_scat_clr                  &
      , ss_prop%phase_fnc_clr, ss_prop%forward_scatter_clr              &
      , nd_profile, nd_layer, 1, nd_layer_clr, nd_ukca_mode             &
      , nd_max_order                                                    &
      )
!   Within clouds:
    CALL opt_prop_ukca_aerosol(ierr                                     &
      , n_profile, n_cloud_top, n_layer                                 &
      , n_order_phase, l_rescale, n_order_forward                       &
      , n_ukca_mode, ukca_modal_mixr                                    &
      , ukca_absorption, ukca_scattering, ukca_asymmetry                &
      , ss_prop%k_grey_tot(1, id_ct, 0)                                 &
      , ss_prop%k_ext_scat(1, id_ct, 0)                                 &
      , ss_prop%phase_fnc(1, id_ct, 1, 0)                               &
      , ss_prop%forward_scatter(1, id_ct, 0)                            &
      , nd_profile, nd_layer, id_ct, nd_layer, nd_ukca_mode             &
      , nd_max_order                                                    &
      )
  END IF

  IF (l_continuum) THEN
!   Include continuum absorption.
    DO j=1, n_continuum
      i_continuum=i_continuum_pointer(j)
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          ss_prop%k_grey_tot_clr(l, i)=ss_prop%k_grey_tot_clr(l, i)     &
            +k_continuum(i_continuum)                                   &
            *amount_continuum(l, i, i_continuum)
        END DO
      END DO
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          ss_prop%k_grey_tot(l, i, 0)=ss_prop%k_grey_tot(l, i, 0)       &
            +k_continuum(i_continuum)                                   &
            *amount_continuum(l, i, i_continuum)
        END DO
      END DO
    END DO
  END IF


! Add the scattering on to the total extinction. The final clear-sky
! phase function not calculated here since the product of the phase
! function and scattering is also needed to calculate the cloudy
! phase function.
  DO i=1, n_cloud_top-1
    DO l=1, n_profile
      ss_prop%k_grey_tot_clr(l, i)=ss_prop%k_grey_tot_clr(l, i)         &
        +ss_prop%k_ext_scat_clr(l, i)
    END DO
  END DO
  DO i=n_cloud_top, n_layer
    DO l=1, n_profile
      ss_prop%k_grey_tot(l, i, 0)=ss_prop%k_grey_tot(l, i, 0)           &
        +ss_prop%k_ext_scat(l, i, 0)
    END DO
  END DO


! If there are no clouds calculate the final optical properties
! and return to the calling routine.
  IF (.NOT.l_cloud) THEN

    IF (l_rescale) THEN
      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          IF (ss_prop%k_ext_scat_clr(l, i) > tiny_k) THEN
            tmp_inv(l)=1.0_RealK/ss_prop%k_ext_scat_clr(l, i)
            ss_prop%forward_scatter_clr(l, i)                           &
              =ss_prop%forward_scatter_clr(l, i)*tmp_inv(l)
            DO ls=1, n_order_phase
              ss_prop%phase_fnc_clr(l, i, ls)                           &
                =ss_prop%phase_fnc_clr(l, i, ls)*tmp_inv(l)
            END DO
          END IF
        END DO
      END DO
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          IF (ss_prop%k_ext_scat(l, i, 0) > tiny_k) THEN
            tmp_inv(l)=1.0_RealK/ss_prop%k_ext_scat(l, i, 0)
            ss_prop%forward_scatter(l, i, 0)                            &
              =ss_prop%forward_scatter(l, i, 0)*tmp_inv(l)
            DO ls=1, n_order_phase
              ss_prop%phase_fnc(l, i, ls, 0)                            &
                =ss_prop%phase_fnc(l, i, ls, 0)*tmp_inv(l)
            END DO
          END IF
        END DO
      END DO
    ELSE
      DO ls=1, n_order_phase
        DO i=1, n_cloud_top-1
          DO l=1, n_profile
            IF (ss_prop%k_ext_scat_clr(l, i) > tiny_k) THEN
              ss_prop%phase_fnc_clr(l, i, ls)                           &
                =ss_prop%phase_fnc_clr(l, i, ls)                        &
                /ss_prop%k_ext_scat_clr(l, i)
            END IF
          END DO
        END DO
        DO i=n_cloud_top, n_layer
          DO l=1, n_profile
            IF (ss_prop%k_ext_scat(l, i, 0) > tiny_k) THEN
              ss_prop%phase_fnc(l, i, ls, 0)                            &
                =ss_prop%phase_fnc(l, i, ls, 0)                         &
                /ss_prop%k_ext_scat(l, i, 0)
            END IF
          END DO
        END DO
      END DO
    END IF

    IF (l_solar_phf) THEN

      DO i=1, n_cloud_top-1
        DO l=1, n_profile
          IF (ss_prop%k_ext_scat_clr(l, i) > tiny_k) THEN
            ss_prop%forward_solar_clr(l, i)                             &
              =ss_prop%forward_solar_clr(l, i)                          &
              /ss_prop%k_ext_scat_clr(l, i)
          END IF
        END DO
      END DO
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          IF (ss_prop%k_ext_scat(l, i, 0) > tiny_k) THEN
            ss_prop%forward_solar(l, i, 0)                              &
              =ss_prop%forward_solar(l, i, 0)                           &
              /ss_prop%k_ext_scat(l, i, 0)
          END IF
        END DO
      END DO

      DO id=1, n_direction
        DO i=1, n_cloud_top-1
          DO l=1, n_profile
            IF (ss_prop%k_ext_scat_clr(l, i) > tiny_k) THEN
              ss_prop%phase_fnc_solar_clr(l, i, id)                     &
                =ss_prop%phase_fnc_solar_clr(l, i, id)                  &
                /ss_prop%k_ext_scat_clr(l, i)
            END IF
          END DO
        END DO
        DO i=n_cloud_top, n_layer
          DO l=1, n_profile
            IF (ss_prop%k_ext_scat(l, i, 0) > tiny_k) THEN
              ss_prop%phase_fnc_solar(l, i, id, 0)                      &
                =ss_prop%phase_fnc_solar(l, i, id, 0)                   &
                /ss_prop%k_ext_scat(l, i, 0)
            END IF
          END DO
        END DO
      END DO

    END IF

    IF (lhook) CALL dr_hook('GREY_OPT_PROP',zhook_out,zhook_handle)
    RETURN

  END IF



! Addition of cloudy properties:


! Add in background contibutions:


! All the processes occurring outside clouds also occur within them.
  DO k=1, n_cloud_type
    DO i=n_cloud_top, n_layer
      DO l=1, n_profile
        ss_prop%k_grey_tot(l, i, k)=ss_prop%k_grey_tot(l, i, 0)
        ss_prop%k_ext_scat(l, i, k)=ss_prop%k_ext_scat(l, i, 0)
        ss_prop%forward_scatter(l, i, k)                                &
          =ss_prop%forward_scatter(l, i, 0)
      END DO
      DO ls=1, n_order_phase
        DO l=1, n_profile
          ss_prop%phase_fnc(l, i, ls, k)                                &
            =ss_prop%phase_fnc(l, i, ls, 0)
        END DO
      END DO
    END DO
!   If using a separate solar phase function that must be initialized.
    IF (l_solar_phf) THEN
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          ss_prop%forward_solar(l, i, k)                                &
            =ss_prop%forward_solar(l, i, 0)
        END DO
      END DO
      DO id=1, n_direction
        DO i=n_cloud_top, n_layer
          DO l=1, n_profile
            ss_prop%phase_fnc_solar(l, i, id, k)                        &
              =ss_prop%phase_fnc_solar(l, i, id, 0)
          END DO
        END DO
      END DO
    END IF
  END DO

! For use with McICA, save the clear-sky phase function (actually the
! product of the phase function and the scattering) before it is
! divided by the mean scatering.
  IF (i_cloud == ip_cloud_mcica .AND. .NOT.l_avg_phase_fnc) THEN
    DO i=n_cloud_top, n_layer
      DO ls=1, n_order_phase
        DO l=1, n_profile
          ss_prop%phase_fnc_no_cloud(l, i, ls)                          &
            =ss_prop%phase_fnc(l, i, ls, 0)
        END DO
      END DO
    END DO

    IF (l_rescale) THEN
      DO i=n_cloud_top, n_layer
        DO l=1, n_profile
          ss_prop%forward_scatter_no_cloud(l, i)                        &
            =ss_prop%forward_scatter(l, i, 0)
        END DO
      END DO
    END IF
  END IF


! Initialize arrays for diagnostic use.
  IF (l_cloud_extinction) THEN
     DO i=1, n_layer
        DO l=1, n_profile
           cloud_extinction(l, i)=0.0_RealK
        END DO
     END DO
  END IF

  IF (l_cloud_absorptivity) THEN
     DO i=1, n_layer
        DO l=1, n_profile
           cloud_absorptivity(l, i)=0.0_RealK
        END DO
     END DO
  END IF

  IF (l_ls_cloud_extinction) THEN
     DO i=1, n_layer
        DO l=1, n_profile
           ls_cloud_extinction(l, i)=0.0_RealK
        END DO
     END DO
  END IF

  IF (l_ls_cloud_absorptivity) THEN
     DO i=1, n_layer
        DO l=1, n_profile
           ls_cloud_absorptivity(l, i)=0.0_RealK
        END DO
     END DO
  END IF

  IF (l_cnv_cloud_extinction) THEN
     DO i=1, n_layer
        DO l=1, n_profile
           cnv_cloud_extinction(l, i)=0.0_RealK
        END DO
     END DO
  END IF

  IF (l_cnv_cloud_absorptivity) THEN
     DO i=1, n_layer
        DO l=1, n_profile
           cnv_cloud_absorptivity(l, i)=0.0_RealK
        END DO
     END DO
  END IF



! Add on the terms representing processes within clouds.

! Loop over the condensed components, calculating their optical
! properties and then assign them to the arrays for the types of
! cloud.

  DO k=1, n_condensed

!   Flags for dealing with components were set in the subroutine
!   set_cloud_pointer. we now determine whether the component is
!   to be included and calculate its optical properties according
!   to the phase of the component. these contributions are added
!   to the arrays for the selected type of cloud.

    IF (l_cloud_cmp(k)) THEN

      IF (i_phase_cmp(k) == ip_phase_water) THEN

!       Include scattering by water droplets.

! DEPENDS ON: opt_prop_water_cloud
        CALL opt_prop_water_cloud(ierr                                  &
          , n_profile, n_layer, n_cloud_top                             &
          , n_cloud_profile, i_cloud_profile                            &
          , n_order_phase, l_rescale, n_order_forward                   &
          , l_henyey_greenstein_pf, l_solar_phf, l_lanczos              &
          , n_order_phase_solar, n_direction, cos_sol_view              &
          , i_condensed_param(k)                                        &
          , condensed_param_list(1, k)                                  &
          , condensed_mix_ratio(1, id_ct, k)                            &
          , condensed_dim_char(1, id_ct, k)                             &
          , ss_prop%k_ext_tot_cloud_comp(1, id_ct, k)                   &
          , ss_prop%k_ext_scat_cloud_comp(1, id_ct, k)                  &
          , ss_prop%phase_fnc_cloud_comp(1, id_ct, 1, k)                &
          , ss_prop%forward_scatter_cloud_comp(1, id_ct, k)             &
          , ss_prop%forward_solar_cloud_comp(1, id_ct, k)               &
          , ss_prop%phase_fnc_solar_cloud_comp(1, id_ct, 1, k)          &
          , nd_profile, nd_radiance_profile, nd_layer, id_ct            &
          , nd_direction, nd_phase_term, nd_max_order                   &
          , nd_cloud_parameter                                          &
          )

      ELSE IF (i_phase_cmp(k) == ip_phase_ice) THEN

!       Include scattering by ice crystals.

! DEPENDS ON: opt_prop_ice_cloud
        CALL opt_prop_ice_cloud(ierr                                    &
          , n_profile, n_layer, n_cloud_top                             &
          , n_cloud_profile, i_cloud_profile                            &
          , n_order_phase, l_rescale, n_order_forward                   &
          , l_henyey_greenstein_pf, l_solar_phf, l_lanczos              &
          , n_order_phase_solar, n_direction, cos_sol_view              &
          , i_condensed_param(k)                                        &
          , condensed_param_list(1, k)                                  &
          , condensed_mix_ratio(1, id_ct, k)                            &
          , condensed_dim_char(1, id_ct, k)                             &
          , t, density                                                  &
          , ss_prop%k_ext_tot_cloud_comp(1, id_ct, k)                   &
          , ss_prop%k_ext_scat_cloud_comp(1, id_ct, k)                  &
          , ss_prop%phase_fnc_cloud_comp(1, id_ct, 1, k)                &
          , ss_prop%forward_scatter_cloud_comp(1, id_ct, k)             &
          , ss_prop%forward_solar_cloud_comp(1, id_ct, k)               &
          , ss_prop%phase_fnc_solar_cloud_comp(1, id_ct, 1, k)          &
          , nd_profile, nd_radiance_profile, nd_layer, id_ct            &
          , nd_direction                                                &
          , nd_phase_term, nd_max_order, nd_cloud_parameter             &
          )

      END IF


!     Increment the arrays of optical properties.

      IF (i_cloud /= ip_cloud_mcica .OR. l_avg_phase_fnc) THEN
      IF (l_rescale) THEN
        DO i=n_cloud_top, n_layer
!CDIR NODEP
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            ss_prop%k_grey_tot(l, i, i_cloud_type(k))                   &
              =ss_prop%k_grey_tot(l, i, i_cloud_type(k))                &
              +ss_prop%k_ext_tot_cloud_comp(l, i, k)
            ss_prop%k_ext_scat(l, i, i_cloud_type(k))                   &
              =ss_prop%k_ext_scat(l, i, i_cloud_type(k))                &
              +ss_prop%k_ext_scat_cloud_comp(l, i, k)
            DO ls=1, n_order_phase
              ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))              &
                =ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))           &
                +ss_prop%phase_fnc_cloud_comp(l, i, ls, k)
            END DO
            ss_prop%forward_scatter(l, i, i_cloud_type(k))              &
              =ss_prop%forward_scatter(l, i, i_cloud_type(k))           &
              +ss_prop%forward_scatter_cloud_comp(l, i, k)
          END DO
        END DO
      ELSE
        DO i=n_cloud_top, n_layer
!CDIR NODEP
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            ss_prop%k_grey_tot(l, i, i_cloud_type(k))                   &
              =ss_prop%k_grey_tot(l, i, i_cloud_type(k))                &
              +ss_prop%k_ext_tot_cloud_comp(l, i, k)
            ss_prop%k_ext_scat(l, i, i_cloud_type(k))                   &
              =ss_prop%k_ext_scat(l, i, i_cloud_type(k))                &
              +ss_prop%k_ext_scat_cloud_comp(l, i, k)
            DO ls=1, n_order_phase
              ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))              &
                =ss_prop%phase_fnc(l, i, ls, i_cloud_type(k))           &
                +ss_prop%phase_fnc_cloud_comp(l, i, ls, k)
            END DO
          END DO
        END DO
      END IF
      IF (l_solar_phf) THEN
        DO i=n_cloud_top, n_layer
!CDIR NODEP
          DO ll=1, n_cloud_profile(i)
            l=i_cloud_profile(ll, i)
            ss_prop%forward_solar(l, i, i_cloud_type(k))                &
              =ss_prop%forward_solar(l, i, i_cloud_type(k))             &
              +ss_prop%forward_solar_cloud_comp(l, i, k)
          END DO
        END DO
        DO i=n_cloud_top, n_layer
!CDIR NODEP
          DO id=1, n_direction
            DO ll=1, n_cloud_profile(i)
              l=i_cloud_profile(ll, i)
              ss_prop%phase_fnc_solar(l, i, id, i_cloud_type(k))        &
                =ss_prop%phase_fnc_solar(l, i, id, i_cloud_type(k))     &
                +ss_prop%phase_fnc_solar_cloud_comp(l, i, id, k)
            END DO
          END DO
        END DO
      END IF
      END IF


!     Extra calculations for diagnostics.

      IF (l_cloud_extinction) THEN
         DO i=n_cloud_top, n_layer
!CDIR NODEP
            DO ll=1, n_cloud_profile(i)
               l=i_cloud_profile(ll, i)
               cloud_extinction(l, i)                                   &
                  =cloud_extinction(l, i)                               &
                  +ss_prop%k_ext_tot_cloud_comp(l, i, k)                &
                  *frac_cloud(l, i, i_cloud_type(k))
            END DO
         END DO
      END IF


      IF (l_cloud_absorptivity) THEN
         DO i=n_cloud_top, n_layer
!CDIR NODEP
            DO ll=1, n_cloud_profile(i)
               l=i_cloud_profile(ll, i)
               cloud_absorptivity(l, i)                                 &
                  =cloud_absorptivity(l, i)                             &
                  +(ss_prop%k_ext_tot_cloud_comp(l, i, k)               &
                  -ss_prop%k_ext_scat_cloud_comp(l, i, k))              &
                  *frac_cloud(l, i, i_cloud_type(k))
            END DO
         END DO
      END IF

      IF ((i_cloud_type(k) == ip_cloud_type_sw).OR.                     &
          (i_cloud_type(k) == ip_cloud_type_si)) THEN

        IF (l_ls_cloud_extinction) THEN
           DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                 l=i_cloud_profile(ll, i)
                 ls_cloud_extinction(l, i)                              &
                    =ls_cloud_extinction(l, i)                          &
                    +ss_prop%k_ext_tot_cloud_comp(l, i, k)              &
                    *frac_cloud(l, i, i_cloud_type(k))
              END DO
           END DO
        END IF


        IF (l_ls_cloud_absorptivity) THEN
           DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                 l=i_cloud_profile(ll, i)
                 ls_cloud_absorptivity(l, i)                            &
                    =ls_cloud_absorptivity(l, i)                        &
                    +(ss_prop%k_ext_tot_cloud_comp(l, i, k)             &
                    -ss_prop%k_ext_scat_cloud_comp(l, i, k))            &
                    *frac_cloud(l, i, i_cloud_type(k))
              END DO
           END DO
        END IF

      ELSE  !  Cloud is of convective type

        IF (l_cnv_cloud_extinction) THEN
           DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                 l=i_cloud_profile(ll, i)
                 cnv_cloud_extinction(l, i)                             &
                    =cnv_cloud_extinction(l, i)                         &
                    +ss_prop%k_ext_tot_cloud_comp(l, i, k)              &
                    *frac_cloud(l, i, i_cloud_type(k))
              END DO
           END DO
        END IF


        IF (l_cnv_cloud_absorptivity) THEN
           DO i=n_cloud_top, n_layer
!CDIR NODEP
              DO ll=1, n_cloud_profile(i)
                 l=i_cloud_profile(ll, i)
                 cnv_cloud_absorptivity(l, i)                           &
                    =cnv_cloud_absorptivity(l, i)                       &
                    +(ss_prop%k_ext_tot_cloud_comp(l, i, k)             &
                    -ss_prop%k_ext_scat_cloud_comp(l, i, k))            &
                    *frac_cloud(l, i, i_cloud_type(k))
              END DO
           END DO
        END IF
      END IF

    END IF

  END DO



! Calculate the final optical properties.
! The scattering was included in the free total extinction earlier,
! but we have yet to divide the product of the phase function and
! the scattering by the mean scattering.

  DO i=1, n_cloud_top-1

    n_index=0
    DO l=1, n_profile
      IF (ss_prop%k_ext_scat_clr(l, i) > tiny_k) THEN
        n_index=n_index+1
        indx(n_index)=l
      END IF
    END DO

    IF (l_rescale) THEN
!CDIR NODEP
      DO k=1, n_index
        tmp_inv(k)=1.0_RealK/ss_prop%k_ext_scat_clr(indx(k), i)
        ss_prop%forward_scatter_clr(indx(k), i)                         &
          =ss_prop%forward_scatter_clr(indx(k), i)*tmp_inv(k)
        DO ls=1, n_order_phase
          ss_prop%phase_fnc_clr(indx(k), i, ls)                         &
            =ss_prop%phase_fnc_clr(indx(k), i, ls)*tmp_inv(k)
        END DO
      END DO
    ELSE
      DO ls=1, n_order_phase
!CDIR NODEP
        DO k=1, n_index
          ss_prop%phase_fnc_clr(indx(k), i, ls)                         &
            =ss_prop%phase_fnc_clr(indx(k), i, ls)                      &
            /ss_prop%k_ext_scat_clr(indx(k), i)
        END DO
      END DO
    END IF

    IF (l_solar_phf) THEN
      DO k=1, n_index
        ss_prop%forward_solar_clr(indx(k), i)                           &
          =ss_prop%forward_solar_clr(indx(k), i)                        &
          /ss_prop%k_ext_scat_clr(indx(k), i)
      END DO
      DO id=1, n_direction
        DO k=1, n_index
          ss_prop%phase_fnc_solar_clr(indx(k), i, id)                   &
            =ss_prop%phase_fnc_solar_clr(indx(k), i, id)                &
            /ss_prop%k_ext_scat_clr(indx(k), i)
        END DO
      END DO
    END IF

  END DO

  DO i=n_cloud_top, n_layer

    n_index=0
    DO l=1, n_profile
      IF (ss_prop%k_ext_scat(l, i, 0) > tiny_k) THEN
        n_index=n_index+1
        indx(n_index)=l
      END IF
    END DO

    IF (l_rescale) THEN
!CDIR NODEP
      DO k=1, n_index
        tmp_inv(k)=1.0_RealK/ss_prop%k_ext_scat(indx(k), i, 0)
        ss_prop%forward_scatter(indx(k), i, 0)                          &
          =ss_prop%forward_scatter(indx(k), i, 0)*tmp_inv(k)
        DO ls=1, n_order_phase
          ss_prop%phase_fnc(indx(k), i, ls, 0)                          &
            =ss_prop%phase_fnc(indx(k), i, ls, 0)*tmp_inv(k)
        END DO
      END DO
    ELSE
      DO ls=1, n_order_phase
!CDIR NODEP
        DO k=1, n_index
          ss_prop%phase_fnc(indx(k), i, ls, 0)                          &
            =ss_prop%phase_fnc(indx(k), i, ls, 0)                       &
            /ss_prop%k_ext_scat(indx(k), i, 0)
        END DO
      END DO
    END IF

    IF (l_solar_phf) THEN
      DO k=1, n_index
        ss_prop%forward_solar(indx(k), i, 0)                            &
          =ss_prop%forward_solar(indx(k), i, 0)                         &
          /ss_prop%k_ext_scat(indx(k), i, 0)
      END DO
      DO id=1, n_direction
        DO k=1, n_index
          ss_prop%phase_fnc_solar(indx(k), i, id, 0)                    &
            =ss_prop%phase_fnc_solar(indx(k), i, id, 0)                 &
            /ss_prop%k_ext_scat(indx(k), i, 0)
        END DO
      END DO
    END IF

  END DO

  IF (i_cloud /= ip_cloud_mcica .OR. l_avg_phase_fnc) THEN

! Repeat for clouds.
  DO k=1, n_cloud_type
    DO i=n_cloud_top, n_layer
      n_index=0
      DO l=1, n_profile
        IF (ss_prop%k_ext_scat(l, i, k) > tiny_k) THEN
          n_index=n_index+1
          indx(n_index)=l
        END IF
      END DO

      IF (l_rescale) THEN
!CDIR NODEP
        DO j=1, n_index
          tmp_inv(j)=1.0_RealK/ss_prop%k_ext_scat(indx(j), i, k)
          ss_prop%forward_scatter(indx(j), i, k)                        &
            =ss_prop%forward_scatter(indx(j), i, k)*tmp_inv(j)
          DO ls=1, n_order_phase
            ss_prop%phase_fnc(indx(j), i, ls, k)                        &
              =ss_prop%phase_fnc(indx(j), i, ls, k)*tmp_inv(j)
          END DO
        END DO
      ELSE
        DO ls=1, n_order_phase
!CDIR NODEP
          DO j=1, n_index
            ss_prop%phase_fnc(indx(j), i, ls, k)                        &
              =ss_prop%phase_fnc(indx(j), i, ls, k)                     &
              /ss_prop%k_ext_scat(indx(j), i, k)
          END DO
        END DO
      END IF

      IF (l_solar_phf) THEN
        DO j=1, n_index
          ss_prop%forward_solar(indx(j), i, k)                          &
            =ss_prop%forward_solar(indx(j), i, k)                       &
            /ss_prop%k_ext_scat(indx(j), i, k)
        END DO
        DO id=1, n_direction
          DO j=1, n_index
            ss_prop%phase_fnc_solar(indx(j), i, id, k)                  &
              =ss_prop%phase_fnc_solar(indx(j), i, id, k)               &
              /ss_prop%k_ext_scat(indx(j), i, k)
          END DO
        END DO
      END IF

    END DO
  END DO

  END IF


  IF (lhook) CALL dr_hook('GREY_OPT_PROP',zhook_out,zhook_handle)

END SUBROUTINE grey_opt_prop
