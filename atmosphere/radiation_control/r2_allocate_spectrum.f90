! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to allocate a reduced spectrum.

! Purpose:
!   To allocate space for the arrays of the reduced spectrum.

! Method:
!   The sizes of the arrays have already been selected. Here we
!   set aside the necessary space.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of code:
!   FORTRAN 90

!- ---------------------------------------------------------------------
      SUBROUTINE r2_allocate_spectrum(spectrum_a, l_ses_spec)




!     Modules used:
      USE dec_spec
      USE missing_data_mod, ONLY: rmdi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE


!     Define the spectrum.
      TYPE (spectrum) spectrum_a

      LOGICAL l_ses_spec
!          Control ses2 radiation

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     The dimensions for the arrays have been set in the calling
!     routine. Space for each array is now allocated and initialised
!     to zero.

      IF (lhook) CALL dr_hook('R2_ALLOCATE_SPECTRUM',zhook_in,zhook_handle)

      ALLOCATE(spectrum_a%l_present(                                    &
          0: spectrum_a%npd_type                                        &
        ))
      spectrum_a%l_present(:)=.FALSE.

      ALLOCATE(spectrum_a%wave_length_short(                            &
          spectrum_a%npd_band                                           &
        ))
      spectrum_a%wave_length_short(:)=0.0

      ALLOCATE(spectrum_a%wave_length_long(                             &
          spectrum_a%npd_band                                           &
        ))
      spectrum_a%wave_length_long(:)=0.0

      ALLOCATE(spectrum_a%n_band_exclude(                               &
          spectrum_a%npd_band                                           &
        ))
      spectrum_a%n_band_exclude(:)=0

      ALLOCATE(spectrum_a%index_exclude(                                &
          spectrum_a%npd_exclude                                        &
        , spectrum_a%npd_band                                           &
        ))
      spectrum_a%index_exclude(:,:)=0

      ALLOCATE(spectrum_a%solar_flux_band(                              &
          spectrum_a%npd_band                                           &
        ))
      spectrum_a%solar_flux_band(:)=0.0

      ALLOCATE(spectrum_a%weight_blue(                                  &
          spectrum_a%npd_band                                           &
        ))
      spectrum_a%weight_blue(:)=rmdi

      ALLOCATE(spectrum_a%rayleigh_coefficient(                         &
          spectrum_a%npd_band                                           &
        ))
      spectrum_a%rayleigh_coefficient(:)=0.0

      ALLOCATE(spectrum_a%n_band_absorb(                                &
          spectrum_a%npd_band                                           &
        ))
      spectrum_a%n_band_absorb(:)=0

      ALLOCATE(spectrum_a%index_absorb(                                 &
          spectrum_a%npd_species                                        &
        , spectrum_a%npd_band                                           &
        ))
      spectrum_a%index_absorb(:,:)=0

      ALLOCATE(spectrum_a%type_absorb(                                  &
          spectrum_a%npd_species                                        &
        ))
      spectrum_a%type_absorb(:)=0

      ALLOCATE(spectrum_a%i_band_esft(                                  &
          spectrum_a%npd_band                                           &
        , spectrum_a%npd_species                                        &
        ))
      spectrum_a%i_band_esft(:,:)=0

      ALLOCATE(spectrum_a%i_scale_esft(                                 &
          spectrum_a%npd_band                                           &
        , spectrum_a%npd_species                                        &
        ))
      spectrum_a%i_scale_esft(:,:)=0

      ALLOCATE(spectrum_a%i_scale_fnc(                                  &
          spectrum_a%npd_band                                           &
        , spectrum_a%npd_species                                        &
        ))
      spectrum_a%i_scale_fnc(:,:)=0

      ALLOCATE(spectrum_a%k_esft(                                       &
          spectrum_a%npd_esft_term                                      &
        , spectrum_a%npd_band                                           &
        , spectrum_a%npd_species                                        &
        ))
      spectrum_a%k_esft(:,:,:)=0.0

      ALLOCATE(spectrum_a%w_esft(                                       &
          spectrum_a%npd_esft_term                                      &
        , spectrum_a%npd_band                                           &
        , spectrum_a%npd_species                                        &
        ))
      spectrum_a%w_esft(:,:,:)=0.0

      ALLOCATE(spectrum_a%scale_vector(                                 &
          spectrum_a%npd_scale_variable                                 &
        , spectrum_a%npd_esft_term                                      &
        , spectrum_a%npd_band                                           &
        , spectrum_a%npd_species                                        &
        ))
      spectrum_a%scale_vector(:,:,:,:)=0.0

      ALLOCATE(spectrum_a%p_reference(                                  &
          spectrum_a%npd_species                                        &
        , spectrum_a%npd_band                                           &
        ))
      spectrum_a%p_reference(:,:)=0.0

      ALLOCATE(spectrum_a%t_reference(                                  &
          spectrum_a%npd_species                                        &
        , spectrum_a%npd_band                                           &
        ))
      spectrum_a%t_reference(:,:)=0.0

      ALLOCATE(spectrum_a%thermal_coefficient(                          &
          0: spectrum_a%npd_thermal_coeff-1                             &
        , spectrum_a%npd_band                                           &
        ))
      spectrum_a%thermal_coefficient(:,:)=0.0

      ALLOCATE( spectrum_a%i_spec_surface(                              &
          spectrum_a%npd_surface                                        &
        ))
      spectrum_a%i_spec_surface(:)=0

      ALLOCATE( spectrum_a%n_dir_albedo_fit(                            &
          spectrum_a%npd_surface                                        &
        ))
      spectrum_a%n_dir_albedo_fit(:)=0

      ALLOCATE( spectrum_a%l_surface(                                   &
          spectrum_a%npd_surface                                        &
        ))
      spectrum_a%l_surface(:)=.FALSE.

      ALLOCATE( spectrum_a%surface_albedo(                              &
          spectrum_a%npd_band                                           &
        , spectrum_a%npd_surface                                        &
        ))
      spectrum_a%surface_albedo(:,:)=0.0

      ALLOCATE( spectrum_a%direct_albedo_parm(                          &
          0:spectrum_a%npd_albedo_parm                                  &
        , spectrum_a%npd_band                                           &
        , spectrum_a%npd_surface                                        &
       ))
      spectrum_a%direct_albedo_parm(:,:,:)=0.0

      ALLOCATE( spectrum_a%emissivity_ground(                           &
          spectrum_a%npd_band                                           &
        , spectrum_a%npd_surface                                        &
       ))
      spectrum_a%emissivity_ground(:,:)=0.0

      ALLOCATE(spectrum_a%n_band_continuum(                             &
          spectrum_a%npd_band                                           &
        ))
      spectrum_a%n_band_continuum(:)=0

      ALLOCATE(spectrum_a%index_continuum(                              &
          spectrum_a%npd_band                                           &
        , spectrum_a%npd_continuum                                      &
        ))
      spectrum_a%index_continuum(:,:)=0

      ALLOCATE(spectrum_a%i_scale_fnc_cont(                             &
          spectrum_a%npd_band                                           &
        , spectrum_a%npd_continuum                                      &
        ))
      spectrum_a%i_scale_fnc_cont(:,:)=0

      ALLOCATE(spectrum_a%k_continuum(                                  &
          spectrum_a%npd_band                                           &
        , spectrum_a%npd_continuum                                      &
        ))
      spectrum_a%k_continuum(:,:)=0.0

      ALLOCATE(spectrum_a%scale_continuum(                              &
          spectrum_a%npd_scale_variable                                 &
        , spectrum_a%npd_band                                           &
        , spectrum_a%npd_continuum                                      &
        ))
      spectrum_a%scale_continuum(:,:,:)=0.0

      ALLOCATE(spectrum_a%p_ref_continuum(                              &
          spectrum_a%npd_continuum                                      &
        , spectrum_a%npd_band                                           &
        ))
      spectrum_a%p_ref_continuum(:,:)=0.0

      ALLOCATE(spectrum_a%t_ref_continuum(                              &
          spectrum_a%npd_continuum                                      &
        , spectrum_a%npd_band                                           &
        ))
      spectrum_a%t_ref_continuum(:,:)=0.0

      ALLOCATE(spectrum_a%i_drop_parametrization(                       &
          spectrum_a%npd_drop_type                                      &
        ))
      spectrum_a%i_drop_parametrization(:)=0

      ALLOCATE(spectrum_a%l_drop_type(                                  &
          spectrum_a%npd_drop_type                                      &
        ))
      spectrum_a%l_drop_type(:)=.FALSE.

      ALLOCATE(spectrum_a%n_drop_phf_term(                              &
          spectrum_a%npd_drop_type                                      &
        ))
      spectrum_a%n_drop_phf_term(:)=0

      ALLOCATE(spectrum_a%drop_parameter_list(                          &
          spectrum_a%npd_cloud_parameter                                &
        , spectrum_a%npd_band                                           &
        , spectrum_a%npd_drop_type                                      &
        ))
      spectrum_a%drop_parameter_list(:,:,:)=0.0

      ALLOCATE(spectrum_a%drop_parm_min_dim(                            &
          spectrum_a%npd_drop_type                                      &
        ))
      spectrum_a%drop_parm_min_dim(:)=0.0

      ALLOCATE(spectrum_a%drop_parm_max_dim(                            &
          spectrum_a%npd_drop_type                                      &
        ))
      spectrum_a%drop_parm_max_dim(:)=0.0

      ALLOCATE(spectrum_a%type_aerosol(                                 &
          spectrum_a%npd_aerosol_species                                &
        ))
      spectrum_a%type_aerosol(:)=0

      ALLOCATE(spectrum_a%i_aerosol_parametrization(                    &
          spectrum_a%npd_aerosol_species                                &
        ))
      spectrum_a%i_aerosol_parametrization(:)=0

      ALLOCATE(spectrum_a%n_aerosol_phf_term(                           &
          spectrum_a%npd_aerosol_species                                &
        ))
      spectrum_a%n_aerosol_phf_term(:)=0

      ALLOCATE(spectrum_a%nhumidity(                                    &
          spectrum_a%npd_aerosol_species                                &
        ))
      spectrum_a%nhumidity(:)=0

      ALLOCATE(spectrum_a%l_aerosol_species(                            &
          spectrum_a%npd_aerosol_species                                &
        ))
      spectrum_a%l_aerosol_species(:)=.FALSE.

      ALLOCATE(spectrum_a%aerosol_absorption(                           &
          spectrum_a%npd_humidities                                     &
        , spectrum_a%npd_aerosol_species                                &
        , spectrum_a%npd_band                                           &
        ))
      spectrum_a%aerosol_absorption(:,:,:)=0.0

      ALLOCATE(spectrum_a%aerosol_scattering(                           &
          spectrum_a%npd_humidities                                     &
        , spectrum_a%npd_aerosol_species                                &
        , spectrum_a%npd_band                                           &
        ))
      spectrum_a%aerosol_scattering(:,:,:)=0.0

      ALLOCATE(spectrum_a%aerosol_phase_fnc(                            &
          spectrum_a%npd_humidities                                     &
        , spectrum_a%npd_phase_term                                     &
        , spectrum_a%npd_aerosol_species                                &
        , spectrum_a%npd_band                                           &
        ))
      spectrum_a%aerosol_phase_fnc(:,:,:,:)=0.0

      ALLOCATE(spectrum_a%humidities(                                   &
          spectrum_a%npd_humidities                                     &
        , spectrum_a%npd_aerosol_species                                &
        ))
      spectrum_a%humidities(:,:)=0.0

      ALLOCATE(spectrum_a%i_ice_parametrization(                        &
          spectrum_a%npd_ice_type                                       &
        ))
      spectrum_a%i_ice_parametrization(:)=0

      ALLOCATE(spectrum_a%l_ice_type(                                   &
          spectrum_a%npd_ice_type                                       &
        ))
      spectrum_a%l_ice_type(:)=.FALSE.

      ALLOCATE(spectrum_a%n_ice_phf_term(                               &
          spectrum_a%npd_ice_type                                       &
        ))
      spectrum_a%n_ice_phf_term(:)=0

      ALLOCATE(spectrum_a%ice_parameter_list(                           &
          spectrum_a%npd_cloud_parameter                                &
        , spectrum_a%npd_band                                           &
        , spectrum_a%npd_ice_type                                       &
        ))
      spectrum_a%ice_parameter_list(:,:,:)=0.0

      ALLOCATE(spectrum_a%ice_parm_min_dim(                             &
          spectrum_a%npd_ice_type                                       &
        ))
      spectrum_a%ice_parm_min_dim(:)=0.0

      ALLOCATE(spectrum_a%ice_parm_max_dim(                             &
          spectrum_a%npd_ice_type                                       &
        ))
      spectrum_a%ice_parm_max_dim(:)=0.0

      ALLOCATE(spectrum_a%l_doppler_present(                            &
          spectrum_a%npd_species                                        &
        ))
      spectrum_a%l_doppler_present(:)=.FALSE.

      ALLOCATE(spectrum_a%doppler_correction(                           &
          spectrum_a%npd_species                                        &
        ))
      spectrum_a%doppler_correction(:)=0.0

      ALLOCATE(spectrum_a%aod_wavel(                                    &
          spectrum_a%npd_aod_wavel                                      &
        ))
      spectrum_a%aod_wavel(:)=0.0

      ALLOCATE(spectrum_a%i_aod_type(                                   &
          spectrum_a%npd_aerosol_species                                &
        ))
      spectrum_a%i_aod_type(:)=0

      ALLOCATE(spectrum_a%aod_absorption(                               &
          spectrum_a%npd_humidities                                     &
        , spectrum_a%npd_aerosol_species                                &
        , spectrum_a%npd_aod_wavel                                      &
        ))
      spectrum_a%aod_absorption(:,:,:)=0.0

      ALLOCATE(spectrum_a%aod_scattering(                               &
          spectrum_a%npd_humidities                                     &
        , spectrum_a%npd_aerosol_species                                &
        , spectrum_a%npd_aod_wavel                                      &
        ))
      spectrum_a%aod_scattering(:,:,:)=0.0


      IF ( l_ses_spec ) THEN

!     Some spectral variables in ses2 have different demonsion from ES and
!     these are renamed to avoid conflict.
!     SES2 also needs several new spectral variables as allocated below

         ALLOCATE(spectrum_a%solar_flux_band_ses(                       &
            spectrum_a%npd_esft_term,                                   &
            spectrum_a%npd_band                                         &
        ))
         spectrum_a%solar_flux_band_ses(:, :)=0.0

         ALLOCATE(spectrum_a%i_band_esft_ses(                           &
             spectrum_a%npd_band                                        &
          ))
         spectrum_a%i_band_esft_ses(:)=0

         ALLOCATE(spectrum_a%k_esft_ses(                                &
             spectrum_a%npd_pre                                         &
           , spectrum_a%npd_tmp                                         &
           , spectrum_a%npd_species                                     &
           , spectrum_a%npd_esft_term                                   &
           , spectrum_a%npd_band                                        &
           ))
!        spectrum_a%k_esft_ses(:,:,:,:,:)=0.0

         ALLOCATE(spectrum_a%w_esft_ses(                                &
             spectrum_a%npd_esft_term                                   &
           , spectrum_a%npd_band                                        &
           ))
!        spectrum_a%w_esft_ses(:,:)=0.0

         ALLOCATE(spectrum_a%k_continuum_ses(                           &
             spectrum_a%npd_esft_term                                   &
           , spectrum_a%npd_tmp                                         &
           , spectrum_a%npd_band                                        &
           , spectrum_a%npd_continuum                                   &
           ))
!        spectrum_a%k_continuum_ses(:,:,:,:)=0.0

         ALLOCATE(spectrum_a%k_h2oc(                                    &
             spectrum_a%npd_pre                                         &
           , spectrum_a%npd_tmp                                         &
           , spectrum_a%npd_esft_term                                   &
           , spectrum_a%npd_band                                        &
           ))
!        spectrum_a%k_h2oc(:,:,:,:)=0.0

         ALLOCATE(spectrum_a%n_mix_gas(                                 &
             spectrum_a%npd_band                                        &
           ))
         spectrum_a%n_mix_gas(:)=0

         ALLOCATE(spectrum_a%index_mix_gas(                             &
             2, spectrum_a%npd_band_mix_gas                             &
           ))
         spectrum_a%index_mix_gas(:,:)=0

         ALLOCATE(spectrum_a%mix_gas_band(                              &
             spectrum_a%npd_band                                        &
           ))
         spectrum_a%mix_gas_band(:)=0

         ALLOCATE(spectrum_a%num_mix(                                   &
             spectrum_a%npd_band                                        &
           ))
         spectrum_a%num_mix(:)=0

         ALLOCATE(spectrum_a%num_ref_p(                                 &
             spectrum_a%npd_band                                        &
           ))
         spectrum_a%num_ref_p(:)=0

         ALLOCATE(spectrum_a%k_mix_gas(                                 &
             spectrum_a%npd_pre                                         &
           , spectrum_a%npd_tmp                                         &
           , spectrum_a%npd_mix                                         &
           , spectrum_a%npd_esft_term                                   &
           , spectrum_a%npd_band_mix_gas                                &
           ))
!        spectrum_a%k_mix_gas(:,:,:,:,:)=0.0

         ALLOCATE(spectrum_a%f_mix(spectrum_a%npd_band))
         spectrum_a%f_mix(:)=0.0

         ALLOCATE( spectrum_a%plk_frac(                                 &
             spectrum_a%npd_esft_term                                   &
           , spectrum_a%npd_band                                        &
           ))
!        spectrum_a%plk_frac(:,:)=0.0

         ALLOCATE( spectrum_a%planckb_ref(                              &
             161, spectrum_a%npd_band                                   &
           ))
!        spectrum_a%planckb_ref(:,:)=0.0

      ELSE

         ALLOCATE(spectrum_a%solar_flux_band_ses(1,1))
         ALLOCATE(spectrum_a%i_band_esft_ses(1)      )
         ALLOCATE(spectrum_a%k_esft_ses(1,1,1,1,1)   )
         ALLOCATE(spectrum_a%w_esft_ses(1,1)         )
         ALLOCATE(spectrum_a%k_continuum_ses(1,1,1,1))
         ALLOCATE(spectrum_a%k_h2oc(1,1,1,1)         )
         ALLOCATE(spectrum_a%n_mix_gas(1)            )
         ALLOCATE(spectrum_a%index_mix_gas(1,1)      )
         ALLOCATE(spectrum_a%mix_gas_band(1)         )
         ALLOCATE(spectrum_a%num_mix(1)              )
         ALLOCATE(spectrum_a%num_ref_p(1)            )
         ALLOCATE(spectrum_a%k_mix_gas(1,1,1,1,1)    )
         ALLOCATE(spectrum_a%f_mix(1)                )
         ALLOCATE(spectrum_a%plk_frac(1,1)           )
         ALLOCATE(spectrum_a%planckb_ref(1,1)        )

      END IF

      IF (lhook) CALL dr_hook('R2_ALLOCATE_SPECTRUM',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE r2_allocate_spectrum
