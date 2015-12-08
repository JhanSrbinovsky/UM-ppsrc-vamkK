! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to transfer spectrum to reduced array.

! Purpose:
!       Spectral data from the large dynamically allocated array
!       are transferred to the reduced array.

! Method:
!       Elements are copied across.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.

!- ---------------------------------------------------------------------
      SUBROUTINE r2_compress_spectrum(                                  &
!                       Original spectrum
          l_present                                                     &
        , n_band, wave_length_short , wave_length_long                  &
        , n_band_exclude, index_exclude                                 &
        , solar_flux_band, weight_blue, rayleigh_coefficient            &
        , n_absorb, n_band_absorb, index_absorb, type_absorb            &
        , l_retain_absorb, n_absorb_retain, index_absorb_retain         &
        , compressed_index, i_band_esft, k_esft, w_esft, i_scale_esft   &
        , i_scale_fnc, scale_vector, p_reference, t_reference           &
        , n_deg_fit, thermal_coefficient, t_ref_planck                  &
        , i_spec_surface, l_surface, surface_albedo                     &
        , n_dir_albedo_fit, direct_albedo_parm, emissivity_ground       &
        , n_band_continuum, index_continuum, index_water                &
        , k_continuum, i_scale_fnc_cont, scale_continuum                &
        , p_ref_continuum, t_ref_continuum                              &
        , l_drop_type, i_drop_parametrization                           &
        , n_drop_phf_term, drop_parameter_list                          &
        , drop_parm_min_dim, drop_parm_max_dim                          &
        , l_ice_type, i_ice_parametrization                             &
        , n_ice_phf_term, ice_parameter_list                            &
        , ice_parm_min_dim, ice_parm_max_dim                            &
        , n_aerosol, type_aerosol                                       &
        , n_aerosol_retain, index_aerosol_retain                        &
        , l_aerosol_species, aerosol_absorption                         &
        , aerosol_scattering, n_aerosol_phf_term, aerosol_phase_fnc     &
        , nhumidity, humidities, i_aerosol_parametrization              &
        , l_doppler_present, doppler_correction, l_use_aod              &
        , n_aod_wavel, aod_wavel                                        &
        , aod_absorption, aod_scattering, i_aod_type                    &
!                       Reduced spectral array
        , spectrum_rd                                                   &
        )



!     Modules used:
      USE rad_pcf
      USE dec_spec
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE


!     ------------------------------------------------------------------
!     Declaration of initial spectrum.
!     ------------------------------------------------------------------
!     ------------------------------------------------------------------
!     MODULE SETTING MAXIMUM DIMENSIONS OF ARRAYS IN THE RADIATION CODE.

      INTEGER,PARAMETER :: NPD_TYPE      = 16 ! number of types of data
      INTEGER,PARAMETER :: NPD_BAND      = 20 ! number of spectral bands
      INTEGER,PARAMETER :: NPD_EXCLUDE   = 2  ! nubmer of excluded bands
      INTEGER,PARAMETER :: NPD_SPECIES   = 11 ! number of gaseous species
      INTEGER,PARAMETER :: NPD_ESFT_TERM = 16 ! number of esft terms
      INTEGER,PARAMETER :: NPD_SCALE_FNC = 5  ! number of scaling funcs
      INTEGER,PARAMETER :: NPD_SOLAR_SRC = 6  ! number of solar source
                                              ! functions
      INTEGER,PARAMETER :: NPD_MIX       = 9  ! number of eta for mixture
                                              ! absorbing species
      INTEGER,PARAMETER :: NPD_BAND_MIX_GAS=3 ! number of bands where
                                              ! mixed species exist
      INTEGER,PARAMETER :: NPD_PRE       = 59 ! number of reference
                                              ! pressure for esft k
      INTEGER,PARAMETER :: NPD_TMP       = 5  ! number of reference
                                              ! temperature for esft k

      ! number of scaling variables
      INTEGER,PARAMETER :: NPD_SCALE_VARIABLE=4

      INTEGER,PARAMETER :: NPD_SURFACE     = 2 ! no of surface types
      INTEGER,PARAMETER :: NPD_ALBEDO_PARM = 4 ! no of albedo parameters
      INTEGER,PARAMETER :: NPD_CONTINUUM   = 4 ! no of continua
      INTEGER,PARAMETER :: NPD_DROP_TYPE   = 6 ! no of drop types
      INTEGER,PARAMETER :: NPD_ICE_TYPE    =16 ! no of ice crystal types


      ! max no of cloud parameters
      INTEGER,PARAMETER :: NPD_CLOUD_PARAMETER=504

      INTEGER, PARAMETER :: NPD_AEROSOL_SPECIES=50

      ! maximum no of humidities
      INTEGER,PARAMETER :: NPD_HUMIDITIES=21

      ! no of thermal coefficients
      INTEGER,PARAMETER :: NPD_THERMAL_COEFF=9

      ! number of wavelengths for aerosol optical depth
      INTEGER,PARAMETER :: NPD_AOD_WAVEL=6

      ! Size allocated for terms in phase functions
      INTEGER,PARAMETER :: NPD_PHASE_TERM=101
! MXSIZE3A end
! SPDEC3A contains declarations for spectral file in two-stream
! radiation code.

! general fields:

      LOGICAL:: L_PRESENT(0:NPD_TYPE) ! flag for types of data present

! properties of the spectral bands:

      INTEGER :: N_BAND ! number of spectral bands

      REAL ::  WAVE_LENGTH_SHORT(NPD_BAND) ! shorter wavelength limits
      REAL ::  WAVE_LENGTH_LONG(NPD_BAND)  ! longer wavelength limits

      ! exclusion of specific bands from parts of the spectrum:

      ! number of excluded bands within each spectral band
      INTEGER :: N_BAND_EXCLUDE(NPD_BAND)

      ! indices of excluded bands
      INTEGER :: INDEX_EXCLUDE(NPD_EXCLUDE,NPD_BAND)

! fields for the solar flux:

      ! fraction of the incident solar flux in each band
      REAL :: SOLAR_FLUX_BAND(NPD_BAND)
      REAL :: weight_blue(npd_band)

      ! fields for Rayleigh scattering:

      REAL ::  RAYLEIGH_COEFFICIENT(NPD_BAND) ! Rayleigh coefficients

      ! fields for gaseous absorption:

      INTEGER :: N_ABSORB ! number of absorbers

      ! number of absorbers in each band
      INTEGER :: N_BAND_ABSORB(NPD_BAND)

      ! list of absorbers in each band
      INTEGER ::  INDEX_ABSORB(NPD_SPECIES, NPD_BAND)

      ! types of each gas in the spectral file
      INTEGER ::  TYPE_ABSORB(NPD_SPECIES)

      ! number of esft terms in band for each gas
      INTEGER ::  I_BAND_ESFT(NPD_BAND, NPD_SPECIES)

      ! type of esft scaling
      INTEGER ::  I_SCALE_ESFT(NPD_BAND, NPD_SPECIES)

      ! type of scaling function
      INTEGER :: I_SCALE_FNC(NPD_BAND, NPD_SPECIES)

      ! esft exponents
      REAL :: K_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)

      ! esft weights
      REAL :: W_ESFT(NPD_ESFT_TERM, NPD_BAND, NPD_SPECIES)

      ! scaling parameters for each absorber and term
      REAL :: SCALE_VECTOR(NPD_SCALE_VARIABLE, NPD_ESFT_TERM, NPD_BAND, &
     &  NPD_SPECIES)

      ! reference pressure for scaling function
      REAL :: P_REFERENCE(NPD_SPECIES, NPD_BAND)

      ! reference temperature for scaling function
      REAL :: T_REFERENCE(NPD_SPECIES, NPD_BAND)

      ! representation of the Planckian:

      INTEGER :: N_DEG_FIT ! degree of thermal polynomial

      ! coefficients in polynomial fit to source function

      REAL :: THERMAL_COEFFICIENT(0: NPD_THERMAL_COEFF-1, NPD_BAND)

      REAL :: T_REF_PLANCK ! planckian reference temperature

! surface properties:

      ! method of specifying properties of surface

      INTEGER :: I_SPEC_SURFACE(NPD_SURFACE)

      ! number of parameters fitting the direct albedo
      INTEGER :: N_DIR_ALBEDO_FIT(NPD_SURFACE)

      LOGICAL :: L_SURFACE(NPD_SURFACE) ! surface types included

      REAL :: SURFACE_ALBEDO(NPD_BAND, NPD_SURFACE) ! surface albedos

! coefficients for fitting direct albedo

      REAL :: DIRECT_ALBEDO_PARM(0:NPD_ALBEDO_PARM,NPD_BAND,NPD_SURFACE)

      ! surface emissivities
      REAL ::   EMISSIVITY_GROUND(NPD_BAND, NPD_SURFACE)


! fields for continua:

      ! number of continua in each band
      INTEGER :: N_BAND_CONTINUUM(NPD_BAND)

      ! list of continua continuua in each band
      INTEGER :: INDEX_CONTINUUM(NPD_BAND, NPD_CONTINUUM)

      INTEGER :: INDEX_WATER ! index of water vapour

      ! type of scaling function for continuum
      INTEGER :: I_SCALE_FNC_CONT(NPD_BAND, NPD_CONTINUUM)

      ! grey extinction coefficients for continuum
      REAL :: K_CONTINUUM(NPD_BAND, NPD_CONTINUUM)

      ! scaling parameters for continuum
      REAL :: SCALE_CONTINUUM(NPD_SCALE_VARIABLE,NPD_BAND,NPD_CONTINUUM)

      ! reference pressure for scaling of continuum
      REAL :: P_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)

      ! reference temperature for scaling of continuum
      REAL :: T_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)

! fields for water droplets:

      ! parametrization type of droplets

      INTEGER :: I_DROP_PARAMETRIZATION(NPD_DROP_TYPE)
      INTEGER :: N_DROP_PHF_TERM(NPD_DROP_TYPE)
!       Number of terms in the phase function for droplets

      ! types of droplet present

      LOGICAL :: L_DROP_TYPE(NPD_DROP_TYPE)

      ! parameters used to fit optical properties of clouds
      REAL :: DROP_PARAMETER_LIST(NPD_CLOUD_PARAMETER,NPD_BAND,         &
     &  NPD_DROP_TYPE)

      ! minimum dimension permissible in the parametrization
      REAL :: DROP_PARM_MIN_DIM(NPD_DROP_TYPE)

      ! maximum dimension permissible in the parametrization
      REAL :: DROP_PARM_MAX_DIM(NPD_DROP_TYPE)

! fields for aerosols:

      INTEGER :: N_AEROSOL ! number of species of aerosol in spectrum

      INTEGER :: N_AEROSOL_MR ! number of aerosols in mixing ratio array

      ! types of aerosols
      INTEGER :: TYPE_AEROSOL(NPD_AEROSOL_SPECIES)

      ! parametrization of aerosols
      INTEGER :: I_AEROSOL_PARAMETRIZATION(NPD_AEROSOL_SPECIES)
      INTEGER :: N_AEROSOL_PHF_TERM(NPD_AEROSOL_SPECIES)
!       Number of terms in the phase function

      ! numbers of humidities
      INTEGER :: NHUMIDITY(NPD_AEROSOL_SPECIES)

      ! aerosol species included

      LOGICAL :: L_AEROSOL_SPECIES(NPD_AEROSOL_SPECIES)

      ! absorption by aerosols
      REAL :: AEROSOL_ABSORPTION(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,   &
     &  NPD_BAND)

      ! scattering by aerosols
      REAL :: AEROSOL_SCATTERING(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,   &
     &  NPD_BAND)

      REAL :: AEROSOL_PHASE_FNC(NPD_HUMIDITIES, NPD_PHASE_TERM          &
     &     ,  NPD_AEROSOL_SPECIES, NPD_BAND)
      REAL :: AEROSOL_ASYMMETRY(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES     &
     &     ,  NPD_BAND)
!       Asymmetry of aerosols

      ! humidities for components
      REAL :: HUMIDITIES(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES)

! fields for ice crystals:

      ! types of parametrization of ice crystals
      INTEGER :: I_ICE_PARAMETRIZATION(NPD_ICE_TYPE)
      INTEGER :: N_ICE_PHF_TERM(NPD_ICE_TYPE)
!       number of terms in the phase function for ice crystals

      ! types of ice crystal present
      LOGICAL :: L_ICE_TYPE(NPD_ICE_TYPE)

      ! parameters used to fit single scattering of ice crystals
      REAL :: ICE_PARAMETER_LIST(NPD_CLOUD_PARAMETER,NPD_BAND,          &
     &   NPD_ICE_TYPE)

      ! minimum dimension permissible in the parametrization
      REAL :: ICE_PARM_MIN_DIM(NPD_ICE_TYPE)

      ! maximum dimension permissible in the parametrization
      REAL :: ICE_PARM_MAX_DIM(NPD_ICE_TYPE)

! fields for doppler broadening:

      ! flag for doppler broadening for each species
      LOGICAL :: L_DOPPLER_PRESENT(NPD_SPECIES)

      ! doppler correction terms
      REAL :: DOPPLER_CORRECTION(NPD_SPECIES)

! fields for aerosol optical depth

      ! number of wavelengths at which the AOD is computed
      INTEGER :: N_AOD_WAVEL

      ! wavelengths at which the AOD is computed
      REAL :: AOD_WAVEL(NPD_AOD_WAVEL)

      ! Monochromatic specific absorption and scattering
      ! coefficients for aerosols
      REAL :: AOD_ABSORPTION(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,       &
     &  NPD_AOD_WAVEL)
      REAL :: AOD_SCATTERING(NPD_HUMIDITIES, NPD_AEROSOL_SPECIES,       &
     &  NPD_AOD_WAVEL)

      ! Aerosol type for each aerosol component
      INTEGER :: I_AOD_TYPE(NPD_AEROSOL_SPECIES)

! SPDEC3A end

!     ------------------------------------------------------------------
!     Auxiliary variables used to select parts of the initial spectrum
!     ------------------------------------------------------------------

      LOGICAL, INTENT(IN) :: l_retain_absorb(npd_species)
!           Flags for the retention of gases in the spectral file
      LOGICAL, INTENT(IN) :: l_use_aod
!           Aerosol optical depth diagnostics are requested.
      INTEGER, INTENT(IN) :: n_absorb_retain
!           Number of absorbers to be retained
      INTEGER, INTENT(IN) :: index_absorb_retain(npd_species)
!           Indices of absorbers to be retained
      INTEGER, INTENT(IN) :: compressed_index(npd_species)
!           Mapping from old to new indices of absorbers
      INTEGER, INTENT(IN) :: n_aerosol_retain
!           Number of aerosols in the initial spectrum to be used
!           In the calculation
      INTEGER, INTENT(IN) :: index_aerosol_retain(npd_aerosol_species)
!           Indices of the retained aerosols


!     ------------------------------------------------------------------
!     Declaration of reduced spectrum.
!     ------------------------------------------------------------------

      TYPE (spectrum) :: spectrum_rd



!     Local variables
      INTEGER                                                           &
          i                                                             &
!           Loop variable
        , j                                                             &
!           Loop variable
        , k                                                             &
!           Loop variable
        , l                                                             &
!           Loop variable
        , ls                                                            &
!           Loop variable
        , n_parameter                                                   &
!           Number of parameters in scheme
        , i_initial                                                     &
!           Indexing number in initial spectrum
        , i_species                                                     &
!           Species of gas
        , i_continuum
!           Type of continuum

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      INTEGER  n_scale_variable(0: npd_scale_fnc)
!           Number of scaling variables

      DATA n_scale_variable(ip_scale_power_law)/2/
      DATA n_scale_variable(ip_scale_fnc_null)/0/
      DATA n_scale_variable(ip_scale_power_quad)/3/
      DATA n_scale_variable(ip_scale_doppler_quad)/4/
      DATA n_scale_variable(ip_scale_wenyi)/0/


      IF (lhook) CALL dr_hook('R2_COMPRESS_SPECTRUM',zhook_in,zhook_handle)


!     Initailize all blocks of the compressed spectrum to .FALSE.
      DO i=0, spectrum_rd%npd_type
        spectrum_rd%l_present(i)=.FALSE.
      END DO


!     Proceed through each block of the spectral file transferring
!     the data from the input array to the reduced array.


!     Block 0:

      IF (l_present(0)) THEN
        spectrum_rd%l_present(0)=.TRUE.
        spectrum_rd%n_band=n_band
        spectrum_rd%n_absorb=n_absorb_retain
        spectrum_rd%n_aerosol=n_aerosol_retain
        DO i=1, n_absorb_retain
          spectrum_rd%type_absorb(i)                                    &
            =type_absorb(index_absorb_retain(i))
        END DO
        DO i=1, n_aerosol_retain
          spectrum_rd%type_aerosol(i)                                   &
            =type_aerosol(index_aerosol_retain(i))
        END DO
      END IF

!     Block 1:
      IF (l_present(1)) THEN
        spectrum_rd%l_present(1)=.TRUE.
        DO i=1, n_band
          spectrum_rd%wave_length_short(i)=wave_length_short(i)
          spectrum_rd%wave_length_long(i)=wave_length_long(i)
        END DO
      END IF

!     Block 2:
      IF (l_present(2)) THEN
        spectrum_rd%l_present(2)=.TRUE.
        DO i=1, n_band
          spectrum_rd%solar_flux_band(i)=solar_flux_band(i)
          spectrum_rd%weight_blue(i)=weight_blue(i)
        END DO
      END IF

!     Block 3:
      IF (l_present(3)) THEN
        spectrum_rd%l_present(3)=.TRUE.
        DO i=1, n_band
          spectrum_rd%rayleigh_coefficient(i)=rayleigh_coefficient(i)
        END DO
      END IF

!     Block 4:
      IF (l_present(4)) THEN
        spectrum_rd%l_present(4)=.TRUE.
        DO i=1, n_band
          spectrum_rd%n_band_absorb(i)=0
          DO j=1, n_band_absorb(i)
            IF (l_retain_absorb(index_absorb(j, i))) THEN
              spectrum_rd%n_band_absorb(i)                              &
                =spectrum_rd%n_band_absorb(i)+1
              spectrum_rd%index_absorb(spectrum_rd%n_band_absorb(i), i) &
                =compressed_index(index_absorb(j, i))
            END IF
          END DO
        END DO
      END IF

!     Block 5:
      IF (l_present(5)) THEN
        spectrum_rd%l_present(5)=.TRUE.

        DO i=1, n_band
          DO j=1, spectrum_rd%n_band_absorb(i)

            i_species=spectrum_rd%index_absorb(j, i)
            i_initial=index_absorb_retain(i_species)

            spectrum_rd%i_band_esft(i, i_species)                       &
              =i_band_esft(i, i_initial)
            spectrum_rd%i_scale_esft(i, i_species)                      &
              =i_scale_esft(i, i_initial)
            spectrum_rd%i_scale_fnc(i, i_species)                       &
              =i_scale_fnc(i, i_initial)
            spectrum_rd%p_reference(i_species, i)                       &
              =p_reference(i_initial, i)
            spectrum_rd%t_reference(i_species, i)                       &
              =t_reference(i_initial, i)

            DO k=1, i_band_esft(i, i_initial)
              spectrum_rd%k_esft(k, i, i_species)                       &
                =k_esft(k, i, i_initial)
              spectrum_rd%w_esft(k, i, i_species)                       &
                =w_esft(k, i, i_initial)
              DO l=1, n_scale_variable(i_scale_fnc(i, i_initial))
                spectrum_rd%scale_vector(l, k, i, i_species)            &
                  =scale_vector(l, k, i, i_initial)
              END DO
            END DO

          END DO
        END DO
      END IF

!     Block 6:
      IF (l_present(6)) THEN
        spectrum_rd%l_present(6)=.TRUE.
        spectrum_rd%n_deg_fit=n_deg_fit
        spectrum_rd%t_ref_planck=t_ref_planck
        DO i=1, n_band
          DO j=0, n_deg_fit
            spectrum_rd%thermal_coefficient(j, i)                       &
              =thermal_coefficient(j, i)
          END DO
        END DO
      END IF

!     Block 8:
      IF (l_present(8)) THEN
        spectrum_rd%l_present(8)=.TRUE.
        DO i=1, n_band
          spectrum_rd%n_band_continuum(i)=n_band_continuum(i)
          DO j=1, n_band_continuum(i)
            spectrum_rd%index_continuum(i, j)=index_continuum(i, j)
          END DO
        END DO

        spectrum_rd%index_water=0
        DO i=1, n_absorb_retain
          IF (index_absorb_retain(i) == index_water) THEN
            spectrum_rd%index_water=i
          END IF
        END DO

      END IF

!     Block 9:
      IF (l_present(9)) THEN
        spectrum_rd%l_present(9)=.TRUE.
        DO i=1, n_band
          DO j=1, n_band_continuum(i)
            i_continuum=index_continuum(i, j)
            spectrum_rd%i_scale_fnc_cont(i, i_continuum)                &
              =i_scale_fnc_cont(i, i_continuum)
            spectrum_rd%p_ref_continuum(i_continuum, i)                 &
              =p_ref_continuum(i_continuum, i)
            spectrum_rd%t_ref_continuum(i_continuum, i)                 &
              =t_ref_continuum(i_continuum, i)

            spectrum_rd%k_continuum(i, i_continuum)                     &
              =k_continuum(i, i_continuum)
            DO l=1, n_scale_variable(i_scale_fnc_cont                   &
                    (i, i_continuum))
              spectrum_rd%scale_continuum(l, i, i_continuum)            &
                =scale_continuum(l, i, i_continuum)
            END DO
          END DO
        END DO
      END IF

!     Block 10:
      IF (l_present(10)) THEN
        spectrum_rd%l_present(10)=.TRUE.
        DO i=1, spectrum_rd%npd_drop_type
          IF (l_drop_type(i)) THEN
            spectrum_rd%l_drop_type(i)=.TRUE.
            spectrum_rd%i_drop_parametrization(i)                       &
              =i_drop_parametrization(i)
            spectrum_rd%n_drop_phf_term(i)                              &
              =n_drop_phf_term(i)
            spectrum_rd%drop_parm_min_dim(i)=drop_parm_min_dim(i)
            spectrum_rd%drop_parm_max_dim(i)=drop_parm_max_dim(i)
            IF (i_drop_parametrization(i) == ip_slingo_schrecker) THEN
              n_parameter=6
            ELSE IF (i_drop_parametrization(i) ==                       &
                      ip_slingo_schr_phf) THEN
              n_parameter=2*n_drop_phf_term(i)+4
            ELSE IF (i_drop_parametrization(i) ==                       &
                     ip_ackerman_stephens) THEN
              n_parameter=9
            ELSE IF (i_drop_parametrization(i) == ip_drop_pade_2) THEN
              n_parameter=16
            END IF

            DO j=1, n_parameter
              DO k=1, n_band
                spectrum_rd%drop_parameter_list(j, k, i)                &
                  =drop_parameter_list(j, k, i)
              END DO
            END DO
          ELSE
            spectrum_rd%l_drop_type(i)=.FALSE.
          END IF
        END DO
      END IF

!     Block 11:
      IF (l_present(11)) THEN
        spectrum_rd%l_present(11)=.TRUE.
        DO i=1, n_aerosol_retain
          i_initial=index_aerosol_retain(i)
          IF (l_aerosol_species(i_initial)) THEN
            spectrum_rd%l_aerosol_species(i)=.TRUE.
            spectrum_rd%i_aerosol_parametrization(i)                    &
              =i_aerosol_parametrization(i_initial)

! *** Changed by JCT 26/05/04 ****

! Consider Full Phase Function:

            IF ( (i_aerosol_parametrization(i_initial) ==               &
                  ip_aerosol_param_phf_dry).OR.                         &
                 (i_aerosol_parametrization(i_initial) ==               &
                ip_aerosol_param_phf_moist) ) THEN
              spectrum_rd%n_aerosol_phf_term(i)                         &
                =n_aerosol_phf_term(i_initial)

! Consider only the asymmetry

            ELSE IF ( (i_aerosol_parametrization(i_initial) ==          &
                  ip_aerosol_param_dry).OR.                             &
                 (i_aerosol_parametrization(i_initial) ==               &
                ip_aerosol_param_moist) ) THEN
              n_aerosol_phf_term(i_initial)=1
              spectrum_rd%n_aerosol_phf_term(i)=1
            END IF

! Read in Asymmetry or phase function for dry aerosols

            IF ( (i_aerosol_parametrization(i_initial) ==               &
                  ip_aerosol_param_dry).OR.                             &
                 (i_aerosol_parametrization(i_initial) ==               &
                ip_aerosol_param_phf_dry) ) THEN
              spectrum_rd%nhumidity(i)=0
              DO k=1, n_band
                spectrum_rd%aerosol_absorption(1, i, k)                 &
                  =aerosol_absorption(1, i_initial, k)
                spectrum_rd%aerosol_scattering(1, i, k)                 &
                  =aerosol_scattering(1, i_initial, k)
                DO ls=1, n_aerosol_phf_term(i_initial)
                  spectrum_rd%aerosol_phase_fnc(1, ls, i, k)            &
                    =aerosol_phase_fnc(1, ls, i_initial, k)
                END DO
              END DO

! Read in Asymmetry or phase function for wet aerosols

            ELSE IF ( (i_aerosol_parametrization(i_initial) ==          &
                     ip_aerosol_param_moist).OR.                        &
                      (i_aerosol_parametrization(i_initial) ==          &
                     ip_aerosol_param_phf_moist) ) THEN
              spectrum_rd%index_water=index_water
              spectrum_rd%nhumidity(i)=nhumidity(i_initial)
              DO j=1, nhumidity(i_initial)
                spectrum_rd%humidities(j, i)=humidities(j, i_initial)
                DO k=1, n_band
                  spectrum_rd%aerosol_absorption(j, i, k)               &
                    =aerosol_absorption(j, i_initial, k)
                  spectrum_rd%aerosol_scattering(j, i, k)               &
                    =aerosol_scattering(j, i_initial, k)
                  DO ls=1, n_aerosol_phf_term(i_initial)
                    spectrum_rd%aerosol_phase_fnc(j, ls, i, k)          &
                      =aerosol_phase_fnc(j, ls, i_initial, k)
                  END DO
                END DO
              END DO
            END IF

          ELSE
            spectrum_rd%l_aerosol_species(i)=.FALSE.
          END IF
        END DO
      END IF

! ********************************

!     Block 12:
      IF (l_present(12)) THEN
        spectrum_rd%l_present(12)=.TRUE.
        DO i=1, spectrum_rd%npd_ice_type
          IF (l_ice_type(i)) THEN
            spectrum_rd%l_ice_type(i)=.TRUE.
            spectrum_rd%ice_parm_min_dim(i)=ice_parm_min_dim(i)
            spectrum_rd%ice_parm_max_dim(i)=ice_parm_max_dim(i)

            spectrum_rd%i_ice_parametrization(i)                        &
              =i_ice_parametrization(i)
            spectrum_rd%n_ice_phf_term(i)                               &
              =n_ice_phf_term(i)
            IF (i_ice_parametrization(i) ==                             &
                ip_slingo_schrecker_ice) THEN
              n_parameter=6
            ELSE IF (i_ice_parametrization(i) ==                        &
                     ip_sun_shine_vn2_vis) THEN
              n_parameter=6
            ELSE IF (i_ice_parametrization(i) ==                        &
                     ip_sun_shine_vn2_ir) THEN
              n_parameter=0
            ELSE IF (i_ice_parametrization(i) == ip_ice_adt) THEN
              n_parameter=30
            ELSE IF (i_ice_parametrization(i) == ip_ice_adt_10) THEN
              n_parameter=36
            ELSE IF (i_ice_parametrization(i) == ip_ice_fu_phf) THEN
              n_parameter=5*n_ice_phf_term(i)+9
            ELSE IF (i_ice_parametrization(i) == ip_ice_agg_de) THEN
              n_parameter=14
            ELSE IF (i_ice_parametrization(i) == ip_ice_t_iwc) THEN
              n_parameter=9
            END IF

            DO j=1, n_parameter
              DO k=1, n_band
                spectrum_rd%ice_parameter_list(j, k, i)                 &
                  =ice_parameter_list(j, k, i)
              END DO
            END DO

          ELSE
            spectrum_rd%l_ice_type(i)=.FALSE.
          END IF
        END DO
      END IF

!     Block 13:
      IF (l_present(13)) THEN
        spectrum_rd%l_present(13)=.TRUE.
        DO i=1, n_absorb
          IF (l_retain_absorb(i)) THEN
            spectrum_rd%l_doppler_present(compressed_index(i))          &
              =l_doppler_present(i)
            IF (l_doppler_present(i))                                   &
              spectrum_rd%doppler_correction(compressed_index(i))       &
                =doppler_correction(i)
          END IF
        END DO
      ELSE
        DO i=1, n_absorb_retain
          spectrum_rd%l_doppler_present(i)=.FALSE.
        END DO
      END IF


!     Block 14:
      IF (l_present(14)) THEN
        spectrum_rd%l_present(14)=.TRUE.
        DO i=1, n_band
          spectrum_rd%n_band_exclude(i)=n_band_exclude(i)
          DO j=1, n_band_exclude(i)
            spectrum_rd%index_exclude(j, i)=index_exclude(j, i)
          END DO
        END DO
      END IF


!     BLOCK 15 (we rely on the work done for block 11)
      IF (l_present(15).AND.l_use_aod) THEN
         spectrum_rd%l_present(15)=.TRUE.
         spectrum_rd%n_aod_wavel = n_aod_wavel
         DO k=1, n_aod_wavel
           spectrum_rd%aod_wavel(k) = aod_wavel(k)
         END DO

         DO i=1, n_aerosol_retain
            i_initial=index_aerosol_retain(i)
            IF (l_aerosol_species(i_initial)) THEN
               spectrum_rd%i_aod_type(i) = i_aod_type(i_initial)
               IF (i_aerosol_parametrization(i_initial)                 &
                    ==  ip_aerosol_param_dry) THEN
                  DO k=1, n_aod_wavel
                    spectrum_rd%aod_absorption(1, i, k) =               &
                     aod_absorption(1, i_initial, k)
                    spectrum_rd%aod_scattering(1, i, k) =               &
                     aod_scattering(1, i_initial, k)
                  END DO
               ELSE IF (i_aerosol_parametrization(i_initial)            &
                         ==  ip_aerosol_param_moist) THEN
                  DO j=1, nhumidity(i_initial)
                    DO k=1, n_aod_wavel
                      spectrum_rd%aod_absorption(j, i, k) =             &
                       aod_absorption(j, i_initial, k)
                      spectrum_rd%aod_scattering(j, i, k) =             &
                       aod_scattering(j, i_initial, k)
                    END DO
                  END DO
               END IF
            END IF
         END DO
      END IF


      IF (lhook) CALL dr_hook('R2_COMPRESS_SPECTRUM',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE r2_compress_spectrum
