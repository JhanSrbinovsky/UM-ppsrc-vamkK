! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!***********************************************************************
!* INITIALISATION OF SPECTRAL PARAMETERS REQUIRED IN EDWARDS SLINGO    *
!*     LW RADIATION ROUTINE                                            *
!*                                                                     *
!*                     SZA       JAN. 1997                             *
!***********************************************************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
      SUBROUTINE ses_inilw(ierr, cmessage                               &
        , spectral_file                                                 &
        , l_ch4, l_n2o, l_cfc11, l_cfc12                                &
        , l_cfc113, l_cfc114, l_hcfc22, l_hfc125, l_hfc134a             &
        , l_climat_aerosol                                              &
        , l_dust, l_use_arcldust                                        &
        , l_sulpc_so2, l_use_arclsulp                                   &
        , l_soot, l_use_arclblck                                        &
        , l_biomass, l_use_arclbiom                                     &
        , l_use_seasalt_direct, l_use_arclsslt                          &
        , l_ocff, l_use_arclocff                                        &
        , l_use_biogenic, l_use_arcldlta                                &
        , l_nitrate                                                     &
        , l_murk_rad, l_use_aod                                         &
        , l_gas, l_continuum, l_drop, l_aerosol, l_ice                  &
        , i_solar_src, lw_spectrum                                      &
        )

!     Modules included.
      USE dec_spec
      USE rad_pcf
      USE filenamelength_mod, ONLY:                                     &
          filenamelength
      USE dust_parameters_mod, ONLY: l_twobin_dust

!
! Common deks required by Edwards & Slingo radiation code

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE ereport_mod, ONLY : ereport
      USE gasid3a_mod, ONLY: npd_gases, ip_h2o, ip_co2, ip_o3, ip_n2o,  &
                             ip_co, ip_ch4, ip_o2, ip_no, ip_so2, ip_no2,&
                             ip_nh3, ip_hno3, ip_n2, ip_cfc11, ip_cfc12,&
                             ip_cfc113, ip_hcfc22, ip_hfc125, ip_hfc134a,&
                             ip_cfc114
      USE stdio3a_mod, ONLY: iu_stdin, iu_stdout, iu_err
      IMPLICIT NONE
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

!     Dummy arguments
      CHARACTER  (LEN=filenamelength), INTENT(IN) :: spectral_file
!           Name of file containing the spectral data
      LOGICAL, INTENT(IN) :: l_ch4
!           Absorption by methane is included
      LOGICAL, INTENT(IN) :: l_n2o
!           Absorption by nitrous oxide is included
      LOGICAL, INTENT(IN) :: l_cfc11
!           Absorption by cfc11 is included
      LOGICAL, INTENT(IN) :: l_cfc12
!           Absorption by cfc12 is included
      LOGICAL, INTENT(IN) :: l_cfc113
!           Absorption by cfc113 is included
      LOGICAL, INTENT(IN) :: l_cfc114
!           Absorption by cfc114 is included
      LOGICAL, INTENT(IN) :: l_hcfc22
!           Absorption by hcfc22 is included
      LOGICAL, INTENT(IN) :: l_hfc125
!           Absorption by hfc125 is included
      LOGICAL, INTENT(IN) :: l_hfc134a
!           Absorption by hfc134a is included
      LOGICAL, INTENT(IN) :: l_climat_aerosol
!           Climatological aerosols are to be included
      LOGICAL, INTENT(IN) :: l_sulpc_so2
!           The direct effects of sulphate aerosols are
!           to be included
      LOGICAL, INTENT(IN) :: l_use_arclsulp
!           The direct effects of sulphate aerosols from NWP climatology
!           are to be included
      LOGICAL, INTENT(IN) :: l_dust
!           Dust aerosols are to be included, either for their
!           direct effect or on radiation or for radiative diagnostics
      LOGICAL, INTENT(IN) :: l_use_arcldust
!           Use the direct radiative effects of mineral dust from
!           NWP climatology in the SW
      LOGICAL, INTENT(IN) :: l_biomass
!           Biomass burning aerosols are to be included, either for their
!           direct effect or on radiation or for radiative diagnostics
      LOGICAL, INTENT(IN) :: l_use_arclbiom
!           Use the direct radiative effects of biomass from
!           NWP climatology in the SW
      LOGICAL, INTENT(IN) :: l_soot
!           Soot aerosols are to be included, either for their
!           direct effect or on radiation or for radiative diagnostics
      LOGICAL, INTENT(IN) :: l_use_arclblck
!           Use the direct radiative effects of black-carbon from NWP
!           climatology in the SW
      LOGICAL, INTENT(IN) :: l_use_seasalt_direct
!           Direct effect of Sea-salt aerosols are to be included
      LOGICAL, INTENT(IN) :: l_use_arclsslt
!           Include the direct radiative effect of sea-salt from NWP
!           climatology in the SW
      LOGICAL, INTENT(IN) :: l_use_biogenic
!           Use the biogenic aerosol direct effect in the SW
      LOGICAL, INTENT(IN) :: l_ocff
!           Fossil-fuel OC aerosols are to be included, either for their
!           direct effect or on radiation or for radiative diagnostics
      LOGICAL, INTENT(IN) :: l_use_arclocff
!           Include the direct radiative effect of fossil-fuel organic
!           carbon aerosol from NWP climatology in the SW
      LOGICAL, INTENT(IN) :: l_use_arcldlta
!           Include the direct radiative effect of delta aerosol
!           from NWP climatology in the SW
      LOGICAL, INTENT(IN) :: l_nitrate
!           Nitrate aerosols are to be included, either for their
!           direct effect or on radiation or for radiative diagnostics
      LOGICAL, INTENT(IN) :: l_murk_rad
!           Mesoscale aerosols are to be included
      LOGICAL, INTENT(IN) :: l_use_aod
!           At least one of the aerosol optical depth diagnostics
!           is requested

!     Flags for physical processes
      LOGICAL, INTENT(INOUT) :: l_gas
!           Flag to include gaseous absorption
      LOGICAL, INTENT(INOUT) :: l_continuum
!           Flag to include continuum absorption
      LOGICAL, INTENT(INOUT) :: l_drop
!           Flag to include scattering by droplets
      LOGICAL, INTENT(INOUT) :: l_aerosol
!           Flag to include scattering by aerosols
      LOGICAL, INTENT(INOUT) :: l_ice
!           Flag to include scattering by ice crystals

      INTEGER, INTENT(IN) :: i_solar_src
!           Index of solar source function used
      INTEGER, INTENT(OUT) :: ierr
!           Error flag
      CHARACTER  (LEN=80), INTENT(OUT) :: cmessage

!     Define the reduced LW spectrum.
      TYPE (spectrum) :: lw_spectrum

!     Declare the initial spectrum for dynamic and define
!     an appropriate namelist.

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

      LOGICAL:: l_present(0:npd_type) ! flag for types of data present

! properties of the spectral bands:

      INTEGER :: n_band ! number of spectral bands

      REAL ::  wave_length_short(npd_band) ! shorter wavelength limits
      REAL ::  wave_length_long(npd_band)  ! longer wavelength limits

      ! exclusion of specific bands from parts of the spectrum:

      ! number of excluded bands within each spectral band
      INTEGER :: n_band_exclude(npd_band)

      ! indices of excluded bands
      INTEGER :: index_exclude(npd_exclude,npd_band)

! fields for the solar flux:

      ! fraction of the incident solar flux in each band
      REAL :: solar_flux_band(npd_esft_term, npd_band, npd_solar_src)

      ! fields for Rayleigh scattering:

      REAL ::  rayleigh_coefficient(npd_band) ! Rayleigh coefficients

      ! fields for gaseous absorption:

      INTEGER :: n_absorb ! number of absorbers

      ! number of absorbers in each band
      INTEGER :: n_band_absorb(npd_band)

      ! list of absorbers in each band
      INTEGER ::  index_absorb(npd_species, npd_band)

      ! types of each gas in the spectral file
      INTEGER ::  type_absorb(npd_species)

      ! Number of mixed gases in a band
      INTEGER ::  n_mix_gas(npd_band)

      ! Index  OF Mixed GASOUS ABSORBERS IN EACH BAND
      INTEGER ::  index_mix_gas(2, npd_band_mix_gas)

      ! sequence number (band) for mixed species
      INTEGER ::  mix_gas_band(npd_band)
      ! number of esft terms in band for each gas
      INTEGER ::  i_band_esft(npd_band)

      ! NUMBER OF REFERENCE PRESSURES
      INTEGER ::  num_ref_p(npd_band)


      ! NUMBER OF BINARY PARAMETER
      INTEGER ::  num_mix(npd_band)

      ! type of esft scaling
!     INTEGER ::  I_SCALE_ESFT(NPD_BAND, NPD_SPECIES)

      ! type of scaling function
!     INTEGER :: I_SCALE_FNC(NPD_BAND, NPD_SPECIES)

      ! esft exponents
      REAL :: k_esft(npd_pre, npd_tmp, npd_mix, npd_esft_term, npd_band)

      ! esft weights
      REAL :: w_esft(npd_esft_term, npd_band)

      ! Absorption coefficient for mixture of binary species
      REAL :: k_mix_gas(npd_pre,npd_tmp,npd_mix, npd_esft_term          &
              , npd_band_mix_gas)

      ! Mixing ratio of mixed absorber amount
      REAL :: f_mix(npd_band)
      ! scaling parameters for each absorber and term
!     REAL :: SCALE_VECTOR(NPD_SCALE_VARIABLE, NPD_ESFT_TERM, NPD_BAND, &
!    &  NPD_SPECIES)

      ! reference pressure for scaling function
!     REAL :: P_REFERENCE(NPD_SPECIES, NPD_BAND)

      ! reference temperature for scaling function
!     REAL :: T_REFERENCE(NPD_SPECIES, NPD_BAND)

      ! representation of the Planckian:

      ! PLANCKIAN FRACTION FUNCTION
      REAL :: plk_frac(npd_esft_term, npd_band)

      ! Planckian function at 161 reference temperatures
      REAL :: planckb_ref(161, npd_band)

      INTEGER :: N_DEG_FIT ! degree of thermal polynomial

      ! coefficients in polynomial fit to source function

      REAL :: THERMAL_COEFFICIENT(0: NPD_THERMAL_COEFF-1, NPD_BAND)

      REAL :: T_REF_PLANCK ! planckian reference temperature

! surface properties:

      ! method of specifying properties of surface

      INTEGER :: i_spec_surface(npd_surface)

      ! number of parameters fitting the direct albedo
      INTEGER :: n_dir_albedo_fit(npd_surface)

      LOGICAL :: l_surface(npd_surface) ! surface types included

      REAL :: surface_albedo(npd_band, npd_surface) ! surface albedos

! coefficients for fitting direct albedo

      REAL :: direct_albedo_parm(0:npd_albedo_parm,npd_band,npd_surface)

      ! surface emissivities
      REAL ::   emissivity_ground(npd_band, npd_surface)


! fields for continua:

      ! number of continua in each band
      INTEGER :: n_band_continuum(npd_band)

      ! list of continua continuua in each band
      INTEGER :: index_continuum(npd_band, npd_continuum)

      INTEGER :: index_water ! index of water vapour

      ! type of scaling function for continuum
!     INTEGER :: I_SCALE_FNC_CONT(NPD_BAND, NPD_CONTINUUM)

      ! grey extinction coefficients for continuum
      REAL :: k_continuum(npd_esft_term, 5, npd_band, npd_continuum)

      ! Extinction for water vapour
      REAL :: k_h2oc(npd_pre, npd_tmp, npd_esft_term,npd_band)

      ! scaling parameters for continuum
!     REAL :: SCALE_CONTINUUM(NPD_SCALE_VARIABLE,NPD_BAND,NPD_CONTINUUM)

      ! reference pressure for scaling of continuum
!     REAL :: P_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)

      ! reference temperature for scaling of continuum
!     REAL :: T_REF_CONTINUUM(NPD_CONTINUUM, NPD_BAND)

! fields for water droplets:

      ! parametrization type of droplets

      INTEGER :: i_drop_parametrization(npd_drop_type)
      INTEGER :: n_drop_phf_term(npd_drop_type)
!       Number of terms in the phase function for droplets

      ! types of droplet present

      LOGICAL :: l_drop_type(npd_drop_type)

      ! parameters used to fit optical properties of clouds
      REAL :: drop_parameter_list(npd_cloud_parameter,npd_band,         &
        npd_drop_type)

      ! minimum dimension permissible in the parametrization
      REAL :: drop_parm_min_dim(npd_drop_type)

      ! maximum dimension permissible in the parametrization
      REAL :: drop_parm_max_dim(npd_drop_type)

! fields for aerosols:

      INTEGER :: n_aerosol ! number of species of aerosol in spectral info

      INTEGER :: n_aerosol_mr ! number of aerosols in mixing ratio array

      ! types of aerosols
      INTEGER :: type_aerosol(npd_aerosol_species)

      ! parametrization of aerosols
      INTEGER :: i_aerosol_parametrization(npd_aerosol_species)
      INTEGER :: n_aerosol_phf_term(npd_aerosol_species)
!       Number of terms in the phase function

      ! numbers of humidities
      INTEGER :: nhumidity(npd_aerosol_species)

      ! aerosol species included

      LOGICAL :: l_aerosol_species(npd_aerosol_species)

      ! absorption by aerosols
      REAL :: aerosol_absorption(0:npd_humidities, npd_aerosol_species, &
        npd_band)

      ! scattering by aerosols
      REAL :: aerosol_scattering(0:npd_humidities, npd_aerosol_species, &
        npd_band)
      ! aerosol optical depth ratio
      REAL :: aerosol_asymmetry(0:npd_humidities, npd_aerosol_species   &
           ,  npd_band)
!       Asymmetry of aerosols
      REAL :: aerosol_phase_fnc(0: npd_humidities, npd_phase_term       &
           ,  npd_aerosol_species, npd_band)

      ! humidities for components
      REAL :: humidities(npd_humidities, npd_aerosol_species)

! fields for ice crystals:

      ! types of parametrization of ice crystals
      INTEGER :: i_ice_parametrization(npd_ice_type)
      INTEGER :: n_ice_phf_term(npd_ice_type)
!       number of terms in the phase function for ice crystals

      ! types of ice crystal present
      LOGICAL :: l_ice_type(npd_ice_type)

      ! parameters used to fit single scattering of ice crystals
      REAL :: ice_parameter_list(npd_cloud_parameter,npd_band,          &
         npd_ice_type)

      ! minimum dimension permissible in the parametrization
      REAL :: ice_parm_min_dim(npd_ice_type)

      ! maximum dimension permissible in the parametrization
      REAL :: ice_parm_max_dim(npd_ice_type)

! fields for doppler broadening:

      ! flag for doppler broadening for each species
      LOGICAL :: l_doppler_present(npd_species)

      ! doppler correction terms
      REAL :: doppler_correction(npd_species)

! fields for aerosol optical depth

      ! number of wavelengths at which the AOD is computed
      INTEGER :: n_aod_wavel

      ! wavelengths at which the AOD is computed
      REAL :: aod_wavel(npd_aod_wavel)

      ! Monochromatic specific absorption and scattering
      ! coefficients for aerosols
      REAL :: aod_absorption(npd_humidities, npd_aerosol_species,       &
        npd_aod_wavel)
      REAL :: aod_scattering(npd_humidities, npd_aerosol_species,       &
        npd_aod_wavel)

      ! Aerosol type for each aerosol component
      INTEGER :: i_aod_type(npd_aerosol_species)

! SPDEC3A end


! Local variables

      LOGICAL                                                           &
          l_retain_absorb(npd_species)                                  &
!           Flag set to .TRUE. if the absorber is to be retained
        , l_gas_included(npd_gases)
!           Logical to test for actual gases included
      INTEGER                                                           &
          n_absorb_retain                                               &
!           Number of absorbers to retain
        , index_absorb_retain(npd_species)                              &
!           Indices of absorbers to be retained
        , n_aerosol_retain                                              &
!           Number of aerosols in the spectral file to be retained
!           For the radiative calculation
        , index_aerosol_retain(npd_aerosol_species)                     &
!           Indexing numbers of the retained aerosols
        , compressed_index(npd_species)                                 &
!           Mapping from old to new indices of absorbers
        , n_aerosol_found                                               &
!           Number of aerosols for the current group of processes
!           Found in the spectral file
        , n_aerosol_mixratio                                            &
!           Number of aerosols to be allocated in the aerosol mixing
!           ratio arrays
        , n_parameter_water                                             &
!           Number of parameters for water cloud optical scheme
        , n_parameter_ice
!           Number of parameters for ice cloud optical scheme

!     Other local variables.

      INTEGER                                                           &
          i, j, k, l, m, n, ii                                          &
        , ls, i_initial
!           Loop variable

      LOGICAL                                                           &
          l_asymmetry
!           Flag to check whether aerosol_asymmetry or
!           aerosol_phase_fnc is read in from spectral file

      CHARACTER (LEN=*), PARAMETER ::  RoutineName = 'ses_inilw'

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
 
      DATA drop_parm_min_dim/npd_drop_type*3.5e-07/
      DATA drop_parm_max_dim/npd_drop_type*3.7e-05/
      DATA ice_parm_min_dim/npd_ice_type*10.0e-06/
      DATA ice_parm_max_dim/npd_ice_type*14.5e-05/


      IF (lhook) CALL dr_hook('SES_INILW',zhook_in,zhook_handle)

!     Each block is initialized as missing:
      l_present=.FALSE.

      aerosol_asymmetry = rmdi
      l_asymmetry = .FALSE.

      nhumidity=0

! Initialise the number of terms in the phase function for the case
! where these are not present in the spectral file
      n_ice_phf_term=1
      n_drop_phf_term=1

      ierr = i_normal


! call in_spectrum here for read in spectral data

      PRINT*,'SES file: ',spectral_file
! DEPENDS ON: SES_SPECTRUM
      CALL ses_spectrum(ierr                                            &
         , spectral_file                                                &
         , l_present                                                    &
         , n_band                                                       &
         , n_absorb, n_band_absorb, index_absorb, type_absorb           &
         , n_mix_gas, index_mix_gas, mix_gas_band                       &
         , n_aerosol, type_aerosol                                      &
         , n_band_continuum, index_continuum, index_water               &
         , i_band_esft, num_ref_p, num_mix                              &
         , f_mix, k_esft, k_mix_gas, w_esft                             &
         , wave_length_short , wave_length_long                         &
         , n_band_exclude, index_exclude                                &
         , solar_flux_band, planckb_ref, plk_frac                       &
         , n_deg_fit, thermal_coefficient, t_ref_planck                 &
         , rayleigh_coefficient                                         &
         , i_spec_surface, surface_albedo, n_dir_albedo_fit             &
         , direct_albedo_parm, emissivity_ground, l_surface             &
         , k_h2oc, k_continuum                                          &
         , n_parameter_water, n_parameter_ice                           &
         , l_drop_type, i_drop_parametrization                          &
         , n_drop_phf_term, drop_parameter_list                         &
         , drop_parm_min_dim, drop_parm_max_dim                         &
         , l_ice_type, i_ice_parametrization, ice_parameter_list        &
         , l_aerosol_species, aerosol_absorption                        &
         , aerosol_scattering, aerosol_asymmetry                        &
         , nhumidity, humidities, i_aerosol_parametrization             &
         , n_aod_wavel, i_aod_type, aod_absorption, aod_scattering      &
         , l_doppler_present, doppler_correction                        &
         )

      IF (ierr /= i_normal) THEN
        cmessage = 'Error on exit from ses_spectrum'
        GO TO 9999
      END IF


!     Set reduced dimensions, either from the sizes of the fixed arrays
!     or from the arrays read in.

      lw_spectrum%npd_type=npd_type
      lw_spectrum%npd_pre=npd_pre
      lw_spectrum%npd_tmp =npd_tmp
      lw_spectrum%npd_mix=npd_mix
      lw_spectrum%npd_band_mix_gas=npd_band_mix_gas
      lw_spectrum%npd_band=MAX(n_band, 1)
      lw_spectrum%npd_species=MAX(n_absorb, 1)
      lw_spectrum%npd_surface=npd_surface
      lw_spectrum%npd_continuum=npd_continuum
      lw_spectrum%npd_thermal_coeff=n_deg_fit+1
      lw_spectrum%npd_cloud_parameter=npd_cloud_parameter
      lw_spectrum%npd_aod_wavel=1
      lw_spectrum%npd_phase_term=1


!     Search the spectrum to find maximum dimensions.

      lw_spectrum%npd_exclude=1
      IF (l_present(14)) THEN
        DO i=1, n_band
          lw_spectrum%npd_exclude=MAX(lw_spectrum%npd_exclude           &
            , n_band_exclude(i))
        END DO
      END IF

!     Search the spectrum to find those gases to be retained.
!     Water vapour, carbon dioxide and ozone are included
!     if present, but a warning is printed if they are
!     not included.
      DO i=1, npd_gases
        l_gas_included(i)=.FALSE.
      END DO
      n_absorb_retain=0

      DO i=1, n_absorb

        l_retain_absorb(i)=.FALSE.
        compressed_index(i)=0

        IF ( (type_absorb(i) == ip_h2o).OR.                             &
             (type_absorb(i) == ip_co2).OR.                             &
             (type_absorb(i) == ip_o3).OR.                              &
             ( (type_absorb(i) == ip_ch4).AND.l_ch4 ).OR.               &
             ( (type_absorb(i) == ip_n2o).AND.l_n2o ).OR.               &
             ( (type_absorb(i) == ip_cfc11).AND.l_cfc11 ).OR.           &
             ( (type_absorb(i) == ip_cfc12).AND.l_cfc12 ).OR.           &
             ( (type_absorb(i) == ip_cfc113).AND.l_cfc113 ).OR.         &
             ( (type_absorb(i) == ip_cfc114).AND.l_cfc114 ).OR.         &
             ( (type_absorb(i) == ip_hcfc22).AND.l_hcfc22 ).OR.         &
             ( (type_absorb(i) == ip_hfc125).AND.l_hfc125 ).OR.         &
             ( (type_absorb(i) == ip_hfc134a).AND.l_hfc134a ) ) THEN
          n_absorb_retain=n_absorb_retain+1
          index_absorb_retain(n_absorb_retain)=i
          compressed_index(i)=n_absorb_retain
          l_retain_absorb(i)=.TRUE.
          l_gas_included(type_absorb(i))=.TRUE.
        END IF

      END DO


!     Print warning messages if those gases normally expected
!     are not present. This is done only for the call advancing
!     the integration. Subsequent calls are diagnostic and may
!     be made only for a limited spectral range when it is
!     not appropriate to include all gases.

        IF (.NOT.l_gas_included(ip_h2o)) THEN
          WRITE(iu_err, '(/a, /a)')                                     &
            '*** warning: water vapour is not included in the '         &
            , 'longwave spectral file.'
        END IF

        IF (.NOT.l_gas_included(ip_co2)) THEN
          WRITE(iu_err, '(/a, /a)')                                     &
            '*** warning: carbon dioxide is not included in the '       &
            , 'longwave spectral file.'
        END IF

        IF (.NOT.l_gas_included(ip_o3)) THEN
          WRITE(iu_err, '(/a, /a)')                                     &
            '*** warning: ozone is not included in the '                &
            , 'longwave spectral file.'
        END IF

        IF ((.NOT.l_gas_included(ip_ch4)).AND.l_ch4) THEN
          cmessage = 'Methane is not included in the LW '               &
            //'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          GO TO 9999
        END IF

        IF ((.NOT.l_gas_included(ip_n2o)).AND.l_n2o) THEN
          cmessage = 'N2O is not included in the LW '                   &
            //'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          GO TO 9999
        END IF

        IF ((.NOT.l_gas_included(ip_cfc11)).AND.l_cfc11) THEN
          cmessage = 'CFC11 is not included in the LW '                 &
            //'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          GO TO 9999
        END IF

        IF ((.NOT.l_gas_included(ip_cfc12)).AND.l_cfc12) THEN
          cmessage = 'CFC12 is not included in the LW '                 &
            //'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          GO TO 9999
        END IF

        IF ((.NOT.l_gas_included(ip_cfc113)).AND.l_cfc113) THEN
          cmessage = 'CFC113 is not included in the LW '                &
            //'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          GO TO 9999
        END IF

        IF ((.NOT.l_gas_included(ip_cfc114)).AND.l_cfc114) THEN
          cmessage = 'CFC114 is not included in the LW '                &
            //'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          GO TO 9999
        END IF

        IF ((.NOT.l_gas_included(ip_hcfc22)).AND.l_hcfc22) THEN
          cmessage = 'HCFC22 is not included in the LW '                &
            //'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          GO TO 9999
        END IF

        IF ((.NOT.l_gas_included(ip_hfc125)).AND.l_hfc125) THEN
          cmessage = 'HFC125 is not included in the LW '                &
            //'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          GO TO 9999
        END IF

        IF ((.NOT.l_gas_included(ip_hfc134a)).AND.l_hfc134a) THEN
          cmessage = 'HFC134A is not included in the LW '               &
            //'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          GO TO 9999
        END IF


!     Set an appropriate reduced dimension.
      lw_spectrum%npd_species=MAX(n_absorb_retain, 1)

      lw_spectrum%npd_esft_term=1
      IF (l_present(5)) THEN
        DO i=1, n_band
          DO j=1, n_band_absorb(i)
            IF (l_retain_absorb(index_absorb(j, i)))                    &
              lw_spectrum%npd_esft_term                                 &
                =MAX(lw_spectrum%npd_esft_term, i_band_esft(i))
          END DO
        END DO
      END IF

      lw_spectrum%npd_drop_type=1
      IF (l_present(10)) THEN
        DO i=1, npd_drop_type
          IF (l_drop_type(i)) THEN
            lw_spectrum%npd_drop_type                                   &
              =MAX(lw_spectrum%npd_drop_type, i)
          END IF
        END DO
      END IF

      lw_spectrum%npd_ice_type=1
      IF (l_present(12)) THEN
        DO i=1, npd_ice_type
          IF (l_ice_type(i)) THEN
            lw_spectrum%npd_ice_type                                    &
              =MAX(lw_spectrum%npd_ice_type, i)
          END IF
        END DO
      END IF



!     Aerosols must be treated carefully to allow for various
!     different combinations without requiring the spectral file
!     to be too constrained. Only those required will be retained.

!     Basic initialization to safe values.
      lw_spectrum%npd_humidities=1
      n_aerosol_retain=0
      n_aerosol_mixratio=0

!     Only check for presence of aerosols if they are turned on.
      IF (l_aerosol) THEN

!     Check the spectral file for climatological aerosols
      IF (l_climat_aerosol .OR. l_murk_rad) THEN

        IF (l_present(11)) THEN

!         Search for the aerosols required for this scheme.
          n_aerosol_found=0
          DO i=1, n_aerosol

            IF ( (type_aerosol(i) == ip_water_soluble).OR.              &
                 (type_aerosol(i) == ip_dust_like).OR.                  &
                 (type_aerosol(i) == ip_oceanic).OR.                    &
                 (type_aerosol(i) == ip_soot).OR.                       &
                 (type_aerosol(i) == ip_sulphuric) ) THEN
              n_aerosol_retain=n_aerosol_retain+1
              index_aerosol_retain(n_aerosol_retain)=i
              n_aerosol_found=n_aerosol_found+1
              n_aerosol_mixratio=n_aerosol_mixratio+1
            END IF
          END DO

          IF (n_aerosol_found /= 5) THEN

            ierr=i_err_fatal
            cmessage='The LW spectral file lacks some '                 &
              //'climatological aerosols.'
            GO TO 9999

          END IF
        ELSE

          ierr=i_err_fatal
          cmessage='LW spectral file contains no aerosol data.'
          GO TO 9999

        END IF

      END IF


      IF (l_dust .OR. l_use_arcldust) THEN

         IF (l_present(11)) THEN ! aerosol block present in spec file

! If running with 2 bin dust, retain the spectrum for it, even if the climatology is used
! for radiation, so the spectrum is there for AODs
           IF (l_dust .AND. l_twobin_dust) THEN

!           Search for the aerosols required for this scheme.
            n_aerosol_found=0
            DO i=1, n_aerosol

               IF ( (type_aerosol(i) == ip_twobindust_1).OR.            &
                    (type_aerosol(i) == ip_twobindust_2) ) THEN
                  n_aerosol_retain=n_aerosol_retain+1
                  index_aerosol_retain(n_aerosol_retain)=i
                  n_aerosol_found=n_aerosol_found+1
                  n_aerosol_mixratio=n_aerosol_mixratio+1
               END IF

            END DO

            IF (n_aerosol_found /= 2) THEN

               ierr=i_err_fatal
               cmessage='The LW Spectral file lacks some '              &
                  //'mineral dust aerosol for the two-bin scheme.'
               GO TO 9999

            END IF
           
           END IF ! ends tests on two-bin scheme


           IF ( l_use_arcldust .OR.                                     &
               ( l_dust .AND. .NOT.l_twobin_dust ) ) THEN

!           Search for the aerosols required for this scheme.
            n_aerosol_found=0
            DO i=1, n_aerosol

               IF ( (type_aerosol(i) == ip_dust_1).OR.                  &
                    (type_aerosol(i) == ip_dust_2).OR.                  &
                    (type_aerosol(i) == ip_dust_3).OR.                  &
                    (type_aerosol(i) == ip_dust_4).OR.                  &
                    (type_aerosol(i) == ip_dust_5).OR.                  &
                    (type_aerosol(i) == ip_dust_6) ) THEN
                  n_aerosol_retain=n_aerosol_retain+1
                  index_aerosol_retain(n_aerosol_retain)=i
                  n_aerosol_found=n_aerosol_found+1
                  IF ( l_use_arcldust .AND.                                     &
                     l_dust .AND. .NOT.l_twobin_dust ) THEN
                    ! if both 6 bin dust and aersol clim dust is on
                    ! then both need to be stored in aerosol mix ratio
                    n_aerosol_mixratio=n_aerosol_mixratio+2
                  ELSE 
                    n_aerosol_mixratio=n_aerosol_mixratio+1
                  END IF
               END IF

            END DO

            IF (n_aerosol_found /= 6) THEN

               ierr=i_err_fatal
               cmessage='The LW Spectral file lacks some '              &
                  //'mineral dust aerosol for the six-bin scheme.'
               GO TO 9999

            END IF

           END IF ! ends test on 6 bin dust being required
           
         ELSE

            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            GO TO 9999

         END IF

      END IF

!     Check the spectral file for soot aerosols.
      IF (l_soot .OR. l_use_arclblck) THEN
        IF (l_present(11)) THEN
!         Search for the aerosols required for this scheme.
          n_aerosol_found=0
          DO i=1, n_aerosol
            IF ( (type_aerosol(i) == ip_fresh_soot).OR.                 &
                 (type_aerosol(i) == ip_aged_soot)) THEN
              n_aerosol_retain=n_aerosol_retain+1
              index_aerosol_retain(n_aerosol_retain)=i
              n_aerosol_found=n_aerosol_found+1
              IF ( l_soot .AND. l_use_arclblck ) THEN
                ! if both prog soot and aersol clim soot is on
                ! then both need to be stored in aerosol mix ratio
                n_aerosol_mixratio=n_aerosol_mixratio+2
              ELSE 
                n_aerosol_mixratio=n_aerosol_mixratio+1
              END IF
            END IF
          END DO

          IF (n_aerosol_found /= 2) THEN

            ierr=i_err_fatal
            cmessage='The LW spectral file lacks some '                 &
              //'soot aerosol data.'
            GO TO 9999

          END IF

        ELSE

          ierr=i_err_fatal
          cmessage='LW spectral file contains no soot data.'
          GO TO 9999

        END IF

      END IF

!     Check the spectral file for biomass aerosol modes.
!     (Only required for the direct effect).

      IF (l_biomass .OR. l_use_arclbiom) THEN

         IF (l_present(11)) THEN ! aerosol block present in spec file

!           Search for the aerosols required for this scheme.
            n_aerosol_found=0
            DO i=1, n_aerosol

               IF ( (type_aerosol(i) == ip_biomass_1).OR.               &
                    (type_aerosol(i) == ip_biomass_2) ) THEN
                  n_aerosol_retain=n_aerosol_retain+1
                  index_aerosol_retain(n_aerosol_retain)=i
                  n_aerosol_found=n_aerosol_found+1
                  IF ( l_biomass .AND. l_use_arclbiom ) THEN
                    ! if both prog biom and aersol clim biom is on
                    ! then both need to be stored in aerosol mix ratio
                    n_aerosol_mixratio=n_aerosol_mixratio+2
                  ELSE 
                    n_aerosol_mixratio=n_aerosol_mixratio+1
                  END IF
               END IF

            END DO

            IF (n_aerosol_found /= 2) THEN

               ierr=i_err_fatal
               cmessage='The LW Spectral file lacks some '              &
                  //'biomass aerosol.'
               GO TO 9999

            END IF

         ELSE

            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            GO TO 9999

         END IF

      END IF

!     Check the spectral file for sulphate aerosols. (These are
!     required only for the direct effect).

      IF (l_sulpc_so2 .OR. l_use_arclsulp) THEN

        IF (l_present(11)) THEN

!         Search for the aerosols required for this scheme.
          n_aerosol_found=0
          DO i=1, n_aerosol

            IF ( (type_aerosol(i) == ip_accum_sulphate).OR.             &
                 (type_aerosol(i) == ip_aitken_sulphate) ) THEN
              n_aerosol_retain=n_aerosol_retain+1
              index_aerosol_retain(n_aerosol_retain)=i
              n_aerosol_found=n_aerosol_found+1
              IF ( l_sulpc_so2 .AND. l_use_arclsulp ) THEN
                ! if both prog sulphate and aersol clim sulphate is on
                ! then both need to be stored in aerosol mix ratio
                n_aerosol_mixratio=n_aerosol_mixratio+2
              ELSE 
                n_aerosol_mixratio=n_aerosol_mixratio+1
              END IF
            END IF
          END DO

          IF (n_aerosol_found /= 2) THEN

            ierr=i_err_fatal
            cmessage='The LW spectral file lacks some '                 &
              //'sulphate aerosols.'
            GO TO 9999

          END IF
        ELSE

          ierr=i_err_fatal
          cmessage='LW spectral file contains no aerosol data.'
          GO TO 9999

        END IF

      END IF


!     Check the spectral file for sea-salt aerosols. (These are
!     required only for the direct effect).

      IF (l_use_seasalt_direct .OR. l_use_arclsslt) THEN

        IF (l_present(11)) THEN

!         Search for the aerosols required for this scheme.
          n_aerosol_found=0
          DO i=1, n_aerosol

            IF ( (type_aerosol(i) == ip_seasalt_film).OR.               &
                 (type_aerosol(i) == ip_seasalt_jet) ) THEN
              n_aerosol_retain=n_aerosol_retain+1
              index_aerosol_retain(n_aerosol_retain)=i
              n_aerosol_found=n_aerosol_found+1
              IF ( l_use_seasalt_direct .AND. l_use_arclsslt ) THEN
                ! if both prog seasalt and aersol clim seasalt is on
                ! then both need to be stored in aerosol mix ratio
                n_aerosol_mixratio=n_aerosol_mixratio+2
              ELSE 
                n_aerosol_mixratio=n_aerosol_mixratio+1
              END IF
            END IF
          END DO

          IF (n_aerosol_found /= 2) THEN

            ierr=i_err_fatal
            cmessage='The LW spectral file lacks some '                 &
              //'seasalt aerosols.'
            GO TO 9999

          END IF
        ELSE

          ierr=i_err_fatal
          cmessage='LW spectral file contains no aerosol data.'
          GO TO 9999

        END IF

      END IF

!     Check the spectral file for biogenic aerosol. (direct effect
!     only).

      IF (l_use_biogenic) THEN

         IF (l_present(11)) THEN ! aerosol block present in spec file

!           Search for the aerosol required for this scheme
            n_aerosol_found=0
            DO i=1, n_aerosol

               IF ( type_aerosol(i) == ip_biogenic ) THEN
                 n_aerosol_retain=n_aerosol_retain+1
                 index_aerosol_retain(n_aerosol_retain)=i
                 n_aerosol_found=n_aerosol_found+1
                 n_aerosol_mixratio=n_aerosol_mixratio+1
               END IF

            END DO

            IF (n_aerosol_found /= 1) THEN

               ierr=i_err_fatal
               cmessage='The LW Spectral file lacks some '              &
                  //'biogenic aerosol.'
               GO TO 9999

            END IF

          ELSE

            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            GO TO 9999

         END IF

      END IF

!     Check the spectral file for fossil-fuel organic carbon aerosol.
!     (Only required for the direct effect).

      IF (l_ocff .OR. l_use_arclocff) THEN

         IF (l_present(11)) THEN ! aerosol block present in spec file

!           Search for the aerosols required for this scheme.
            n_aerosol_found=0
            DO i=1, n_aerosol

               IF ( (type_aerosol(i) == ip_ocff_fresh).OR.              &
                    (type_aerosol(i) == ip_ocff_aged) ) THEN
                  n_aerosol_retain=n_aerosol_retain+1
                  index_aerosol_retain(n_aerosol_retain)=i
                  n_aerosol_found=n_aerosol_found+1
                  IF ( l_ocff .AND. l_use_arclocff ) THEN
                    ! if both prog seasalt and aersol clim seasalt is on
                    ! then both need to be stored in aerosol mix ratio
                    n_aerosol_mixratio=n_aerosol_mixratio+2
                  ELSE 
                    n_aerosol_mixratio=n_aerosol_mixratio+1
                  END IF
               END IF

            END DO

            IF (n_aerosol_found /= 2) THEN

               ierr=i_err_fatal
               cmessage='The LW Spectral file lacks some '              &
                  //'fossil fuel org.carb. aerosol.'
               GO TO 9999

            END IF

         ELSE

            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            GO TO 9999

         END IF

      END IF

!     Check the spectral file for delta aerosol. (direct effect
!     only).

      IF (l_use_arcldlta) THEN

         IF (l_present(11)) THEN ! aerosol block present in spec file

!           Search for the aerosol required for this scheme
            n_aerosol_found=0
            DO i=1, n_aerosol

               IF ( type_aerosol(i) == ip_delta ) THEN
                 n_aerosol_retain=n_aerosol_retain+1
                 index_aerosol_retain(n_aerosol_retain)=i
                 n_aerosol_found=n_aerosol_found+1
                 n_aerosol_mixratio=n_aerosol_mixratio+1
               END IF

            END DO

            IF (n_aerosol_found /= 1) THEN

               ierr=i_err_fatal
               cmessage='The LW Spectral file lacks some '              &
                  //'delta aerosol.'
               GO TO 9999

            END IF

          ELSE

            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            GO TO 9999

         END IF

      END IF

!     Check the spectral file for nitrate aerosol. (direct effect
!     only).

      IF (l_nitrate) THEN

         IF (l_present(11)) THEN ! aerosol block present in spec file

!           Search for the aerosol required for this scheme
            n_aerosol_found=0
            DO i=1, n_aerosol

               IF ( type_aerosol(i) == ip_nitrate ) THEN
                 n_aerosol_retain=n_aerosol_retain+1
                 index_aerosol_retain(n_aerosol_retain)=i
                 n_aerosol_found=n_aerosol_found+1
                 n_aerosol_mixratio=n_aerosol_mixratio+1
               END IF

            END DO

            IF (n_aerosol_found /= 1) THEN

               ierr=i_err_fatal
               cmessage='The LW Spectral file lacks some '              &
                  //'nitrate aerosol.'
               GO TO 9999

            END IF

          ELSE

            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            GO TO 9999

         END IF

      END IF



      END IF ! l_aerosol

!     Set an appropriate reduced dimension.
      lw_spectrum%npd_aerosol_species=MAX(n_aerosol_retain, 1)
      lw_spectrum%n_aerosol_mr=n_aerosol_mixratio
      lw_spectrum%npd_aerosol_mixratio=MAX(n_aerosol_mixratio, 1)

!     Set the allowed number of humidities from the number of
!     retained aerosols.

      IF (l_present(11)) THEN
        DO i=1, n_aerosol_retain
          IF ( (i_aerosol_parametrization(index_aerosol_retain(i)) ==   &
                ip_aerosol_param_phf_dry).OR.                           &
               (i_aerosol_parametrization(index_aerosol_retain(i)) ==   &
                ip_aerosol_param_phf_moist) ) THEN
            lw_spectrum%npd_phase_term=MAX(lw_spectrum%npd_phase_term   &
              , n_aerosol_phf_term(index_aerosol_retain(i)))
          END IF
          IF ( (i_aerosol_parametrization(index_aerosol_retain(i)) ==   &
                ip_aerosol_param_moist).OR.                             &
               (i_aerosol_parametrization(index_aerosol_retain(i)) ==   &
                ip_aerosol_param_phf_moist) ) THEN
            lw_spectrum%npd_humidities=MAX(lw_spectrum%npd_humidities   &
              , nhumidity(index_aerosol_retain(i)))
          END IF
        END DO

! Check whether aerosol_phase_fnc or aerosol_asymmetry is specified.

        IF ( ANY(aerosol_asymmetry /= rmdi) ) THEN
           aerosol_phase_fnc(:,1,:,:)=aerosol_asymmetry
        END IF


      END IF


!     Check that block 15 is present if the aerosol optical depth
!     was requested.
      IF (l_use_aod.AND.l_aerosol) THEN
        IF(.NOT. l_present(15)) THEN
          ierr = i_err_fatal
          cmessage='Block 15 needed in the LW spectral file.'
          GO TO 9999
        END IF
       ! Check that the number of wavelengths in the spectral
       ! file is not larger than the number set in the include
       ! file MXSIZE3A.
        IF(n_aod_wavel  >   npd_aod_wavel) THEN
          ierr = i_err_fatal
          cmessage='Increase NPD_AOD_WAVEL in MXSIZE3A.'
          GO TO 9999
        END IF
        lw_spectrum%npd_aod_wavel = npd_aod_wavel
      END IF


!     Allocate space for the LW spectrum.
! DEPENDS ON: r2_allocate_spectrum
      CALL r2_allocate_spectrum(lw_spectrum, .TRUE.)

!     Transfer the large namelist to the reduced spectrum.



!     Initailize all blocks of the compressed spectrum to .FALSE.
      DO i=0, lw_spectrum%npd_type
        lw_spectrum%l_present(i)=.FALSE.
      END DO


!     Proceed through each block of the spectral file transferring
!     the data from the input array to the reduced array.


!     Block 0:

      IF (l_present(0)) THEN
        lw_spectrum%l_present(0)=.TRUE.
        lw_spectrum%n_band=n_band
        lw_spectrum%n_absorb=n_absorb_retain
        lw_spectrum%n_aerosol=n_aerosol_retain
        DO i=1, n_absorb_retain
          lw_spectrum%type_absorb(i)                                    &
            =type_absorb(index_absorb_retain(i))
        END DO
        DO i=1, n_aerosol_retain
          lw_spectrum%type_aerosol(i)                                   &
            =type_aerosol(index_aerosol_retain(i))
        END DO
       END IF

!     Block 1:
      IF (l_present(1)) THEN
        lw_spectrum%l_present(1)=.TRUE.
        DO i=1, n_band
          lw_spectrum%wave_length_short(i)=wave_length_short(i)
          lw_spectrum%wave_length_long(i)=wave_length_long(i)
        END DO
      END IF

!     Block 2:
      IF (l_present(2)) THEN
        lw_spectrum%l_present(2)=.TRUE.
        DO i=1, n_band
           DO k=1,i_band_esft(i)
              lw_spectrum%solar_flux_band_ses(k,i)=                     &
                solar_flux_band(k,i,i_solar_src)
           END DO
        END DO
      END IF

!     Block 3:

!     Block 4:
      IF (l_present(4)) THEN
        lw_spectrum%l_present(4)=.TRUE.
        DO i=1, n_band
          lw_spectrum%n_band_absorb(i)=0
          DO j=1, n_band_absorb(i)
            IF (l_retain_absorb(index_absorb(j, i))) THEN
              lw_spectrum%n_band_absorb(i)                              &
                =lw_spectrum%n_band_absorb(i)+1
              lw_spectrum%index_absorb(lw_spectrum%n_band_absorb(i), i) &
                =compressed_index(index_absorb(j, i))
            END IF
          END DO
        END DO
      END IF

!     Block 5:
      IF (l_present(5)) THEN
        lw_spectrum%l_present(5)=.TRUE.
        lw_spectrum%i_scale_fnc=ip_scale_ses2
        lw_spectrum%i_scale_esft=ip_scale_null

        DO i=1, n_band

          lw_spectrum%num_ref_p(i)=num_ref_p(i)
          lw_spectrum%i_band_esft_ses(i)=i_band_esft(i)
          lw_spectrum%num_mix(i)=num_mix(i)
          lw_spectrum%f_mix(i)=f_mix(i)
          DO j=1, i_band_esft(i)
             lw_spectrum%w_esft_ses(j,i)=w_esft(j,i)
          END DO

          DO j=1, lw_spectrum%n_band_absorb(i)
             DO k=1, i_band_esft(i)
                DO l=1, num_ref_p(i)
                   DO n=1, npd_tmp
                      lw_spectrum%k_esft_ses(l, n, j, k, i)                 &
                           =k_esft(l, n, j, k, i)
                   END DO
                END DO
             END DO
           END DO

        END DO
      END IF



!     Block 6:
      IF (l_present(6)) THEN
        lw_spectrum%l_present(6)=.TRUE.
        lw_spectrum%n_deg_fit=n_deg_fit
        lw_spectrum%t_ref_planck=t_ref_planck
        DO i=1, n_band
          DO j=0, n_deg_fit
            lw_spectrum%thermal_coefficient(j, i)                       &
              =thermal_coefficient(j, i)
          END DO
          DO j=1, i_band_esft(i)
             lw_spectrum%plk_frac(j, i)                                 &
               =plk_frac(j, i)
          END DO

          DO j=1,161
             lw_spectrum%planckb_ref(j, i)=planckb_ref(j, i)
          END DO
        END DO
      END IF

!     Block 8:
      IF (l_present(8)) THEN
        lw_spectrum%l_present(8)=.TRUE.
        DO i=1, n_band
          lw_spectrum%n_band_continuum(i)=n_band_continuum(i)
          DO j=1, n_band_continuum(i)
            lw_spectrum%index_continuum(i, j)=index_continuum(i, j)
          END DO
        END DO

        lw_spectrum%index_water=0
        DO i=1, n_absorb_retain
          IF (index_absorb_retain(i) == index_water) THEN
            lw_spectrum%index_water=i
          END IF
        END DO

      END IF

!     Block 9:
      IF (l_present(9)) THEN
        lw_spectrum%l_present(9)=.TRUE.
        lw_spectrum%i_scale_fnc_cont=ip_scale_ses2
        DO i=1, n_band
           DO j=1, n_band_continuum(i)
              IF( index_continuum(i, j) == 1) THEN

! Cotinuum for water vapour

                 DO k=1, i_band_esft(i)
                    DO n=1, 21
                       DO l=1, npd_tmp
                          lw_spectrum%k_h2oc(n, l, k, i)                &
                             =k_h2oc(n, l, k, i)
                       END DO
                    END DO
                 END DO
              ELSE

! Continuum for oxygen

                 IF ( i  <=  6 ) THEN
                    DO k=1, i_band_esft(i)
                       lw_spectrum%k_continuum_ses(k, 1, i, j)          &
                             =k_continuum(k, 1, i, j)
                    END DO
                 ELSE
                    DO k=1, i_band_esft(i)
                       DO l=1, 3
                          lw_spectrum%k_continuum_ses(k, l, i, j)       &
                             =k_continuum(k, l, i, j)
                       END DO
                    END DO
                 END IF

              END IF
           END DO
        END DO
      END IF

!     Block 10:
      IF (l_present(10)) THEN
        lw_spectrum%l_present(10)=.TRUE.
        DO i=1, lw_spectrum%npd_drop_type
          IF (l_drop_type(i)) THEN
            lw_spectrum%l_drop_type(i)=.TRUE.
            lw_spectrum%i_drop_parametrization(i)                       &
              =i_drop_parametrization(i)
            lw_spectrum%n_drop_phf_term(i)                              &
              =n_drop_phf_term(i)
            lw_spectrum%drop_parm_min_dim(i)=drop_parm_min_dim(i)
            lw_spectrum%drop_parm_max_dim(i)=drop_parm_max_dim(i)
            IF (i_drop_parametrization(i) == ip_slingo_schrecker) THEN
              n_parameter_water=6
            ELSE IF (i_drop_parametrization(i) ==                       &
                      ip_slingo_schr_phf) THEN
              n_parameter_water=2*n_drop_phf_term(i)+4
            ELSE IF (i_drop_parametrization(i) ==                       &
                     ip_ackerman_stephens) THEN
              n_parameter_water=9
            ELSE IF (i_drop_parametrization(i) == ip_drop_pade_2) THEN
              n_parameter_water=16
            ELSE IF (i_drop_parametrization(i) == ip_drop_stamnes) THEN
              n_parameter_water=15
            END IF

            DO j=1, n_parameter_water
              DO k=1, n_band
                lw_spectrum%drop_parameter_list(j, k, i)                &
                  =drop_parameter_list(j, k, i)
              END DO
            END DO
          ELSE
            lw_spectrum%l_drop_type(i)=.FALSE.
          END IF
        END DO
      END IF

!     Block 11:
      IF (l_present(11)) THEN
        lw_spectrum%l_present(11)=.TRUE.
        DO i=1, n_aerosol_retain
          i_initial=index_aerosol_retain(i)
          IF (l_aerosol_species(i_initial)) THEN
            lw_spectrum%l_aerosol_species(i)=.TRUE.
            lw_spectrum%i_aerosol_parametrization(i)                    &
              =i_aerosol_parametrization(i_initial)

!**** Changed by JCT 26/05/04 ****

! Consider Full Phase Function:

            IF ( (i_aerosol_parametrization(i_initial) ==               &
                  ip_aerosol_param_phf_dry).OR.                         &
                 (i_aerosol_parametrization(i_initial) ==               &
                ip_aerosol_param_phf_moist) ) THEN
              lw_spectrum%n_aerosol_phf_term(i)                         &
                =n_aerosol_phf_term(i_initial)

! Consider only the asymmetry

            ELSE IF ( (i_aerosol_parametrization(i_initial) ==          &
                  ip_aerosol_param_dry).OR.                             &
                 (i_aerosol_parametrization(i_initial) ==               &
                ip_aerosol_param_moist) ) THEN
              n_aerosol_phf_term(i_initial)=1
              lw_spectrum%n_aerosol_phf_term(i)=1
            END IF

! Read in Asymmetry or phase function for dry aerosols

            IF ( (i_aerosol_parametrization(i_initial) ==               &
                  ip_aerosol_param_dry).OR.                             &
                 (i_aerosol_parametrization(i_initial) ==               &
                ip_aerosol_param_phf_dry) ) THEN
              lw_spectrum%nhumidity(i)=1

! Initialize humidities here for vectorization purpose
! in grey_extinction

              DO k=1, lw_spectrum%npd_humidities
                 lw_spectrum%humidities(k, i)=humidities(k, i)
              END DO


              DO k=1, n_band
                lw_spectrum%aerosol_absorption(1, i, k)                 &
                  =aerosol_absorption(1, i_initial, k)
                lw_spectrum%aerosol_scattering(1, i, k)                 &
                  =aerosol_scattering(1, i_initial, k)
                DO ls=1, n_aerosol_phf_term(i_initial)
                  lw_spectrum%aerosol_phase_fnc(1, ls, i, k)            &
                    =aerosol_phase_fnc(1, ls, i_initial, k)
                END DO
              END DO

! Read in Asymmetry or phase function for wet aerosols

            ELSE IF ( (i_aerosol_parametrization(i_initial) ==          &
                     ip_aerosol_param_moist).OR.                        &
                      (i_aerosol_parametrization(i_initial) ==          &
                     ip_aerosol_param_phf_moist) ) THEN
              lw_spectrum%index_water=index_water
              lw_spectrum%nhumidity(i)=nhumidity(i_initial)
              DO j=1, nhumidity(i_initial)
                lw_spectrum%humidities(j, i)=humidities(j, i_initial)
                DO k=1, n_band
                  lw_spectrum%aerosol_absorption(j, i, k)               &
                    =aerosol_absorption(j, i_initial, k)
                  lw_spectrum%aerosol_scattering(j, i, k)               &
                    =aerosol_scattering(j, i_initial, k)
                  DO ls=1, n_aerosol_phf_term(i_initial)
                    lw_spectrum%aerosol_phase_fnc(j, ls, i, k)          &
                      =aerosol_phase_fnc(j, ls, i_initial, k)
                  END DO
                END DO
              END DO
            END IF

          ELSE
            lw_spectrum%l_aerosol_species(i)=.FALSE.
          END IF
        END DO
      END IF

!*********************************

!     Block 12:
      IF (l_present(12)) THEN
        lw_spectrum%l_present(12)=.TRUE.
        DO i=1, lw_spectrum%npd_ice_type
          IF (l_ice_type(i)) THEN
            lw_spectrum%l_ice_type(i)=.TRUE.
            lw_spectrum%ice_parm_min_dim(i)=ice_parm_min_dim(i)
            lw_spectrum%ice_parm_max_dim(i)=ice_parm_max_dim(i)

            lw_spectrum%i_ice_parametrization(i)                        &
              =i_ice_parametrization(i)
            lw_spectrum%n_ice_phf_term(i)                               &
              =n_ice_phf_term(i)
            IF (i_ice_parametrization(i) ==                             &
                ip_slingo_schrecker_ice) THEN
              n_parameter_ice=6
            ELSE IF (i_ice_parametrization(i) ==                        &
                     ip_sun_shine_vn2_vis) THEN
              n_parameter_ice=6
            ELSE IF (i_ice_parametrization(i) ==                        &
                     ip_sun_shine_vn2_ir) THEN
              n_parameter_ice=0
            ELSE IF (i_ice_parametrization(i) == ip_ice_adt) THEN
              n_parameter_ice=30
            ELSE IF (i_ice_parametrization(i) == ip_ice_adt_10) THEN
              n_parameter_ice=36
            ELSE IF (i_ice_parametrization(i) == ip_ice_fu_phf) THEN
              n_parameter_ice=5*n_ice_phf_term(i)+9
            ELSE IF (i_ice_parametrization(i) == ip_ice_sun_fu) THEN
              n_parameter_ice=11
            ELSE IF (i_ice_parametrization(i) == ip_ice_chou_vis) THEN
              n_parameter_ice=6
            ELSE IF (i_ice_parametrization(i) == ip_ice_agg_de_sun) THEN
              n_parameter_ice=11
            ELSE IF (i_ice_parametrization(i) == ip_ice_t_iwc) THEN
              n_parameter_ice=9
            END IF

            DO j=1, n_parameter_ice
              DO k=1, n_band
                lw_spectrum%ice_parameter_list(j, k, i)                 &
                  =ice_parameter_list(j, k, i)
              END DO
            END DO

          ELSE
            lw_spectrum%l_ice_type(i)=.FALSE.
          END IF
        END DO
      END IF

!     Block 13:
      IF (l_present(13)) THEN
        lw_spectrum%l_present(13)=.TRUE.
        DO i=1, n_absorb
          IF (l_retain_absorb(i)) THEN
            lw_spectrum%l_doppler_present(compressed_index(i))          &
              =l_doppler_present(i)
            IF (l_doppler_present(i))                                   &
              lw_spectrum%doppler_correction(compressed_index(i))       &
                =doppler_correction(i)
          END IF
        END DO
      ELSE
        DO i=1, n_absorb_retain
          lw_spectrum%l_doppler_present(i)=.FALSE.
        END DO
      END IF


!     Block 14:
      IF (l_present(14)) THEN
        lw_spectrum%l_present(14)=.TRUE.
        DO i=1, n_band
          lw_spectrum%n_band_exclude(i)=n_band_exclude(i)
          DO j=1, n_band_exclude(i)
            lw_spectrum%index_exclude(j, i)=index_exclude(j, i)
          END DO
        END DO
      END IF


!     BLOCK 15 (we rely on the work done for block 11)
      IF (l_present(15).AND.l_use_aod) THEN
         lw_spectrum%l_present(15)=.TRUE.
         lw_spectrum%n_aod_wavel = n_aod_wavel
         lw_spectrum%aod_wavel = (/ 380., 440., 550., 670., 865., 1020. /) 
! wl above are in nanometers, convert to SI:
         lw_spectrum%aod_wavel = lw_spectrum%aod_wavel * 1.0e-9

         DO i=1, n_aerosol_retain
            i_initial=index_aerosol_retain(i)
            IF (l_aerosol_species(i_initial)) THEN
               lw_spectrum%i_aod_type(i) = i_aod_type(i_initial)
               IF (i_aerosol_parametrization(i_initial)                 &
                    ==  ip_aerosol_param_dry) THEN
                  DO k=1, n_aod_wavel
                    lw_spectrum%aod_absorption(1, i, k) =               &
                     aod_absorption(1, i_initial, k)
                    lw_spectrum%aod_scattering(1, i, k) =               &
                     aod_scattering(1, i_initial, k)
                  END DO
               ELSE IF (i_aerosol_parametrization(i_initial)            &
                         ==  ip_aerosol_param_moist) THEN
                  DO j=1, nhumidity(i_initial)
                    DO k=1, n_aod_wavel
                      lw_spectrum%aod_absorption(j, i, k) =             &
                       aod_absorption(j, i_initial, k)
                      lw_spectrum%aod_scattering(j, i, k) =             &
                       aod_scattering(j, i_initial, k)
                    END DO
                  END DO
               END IF
            END IF
         END DO
      END IF


!     Block 16 in ses2 for mixture species in a band
      IF (l_present(16)) THEN
         lw_spectrum%l_present(16)=.TRUE.
         ii=0
         DO i=1, n_band

            lw_spectrum%n_mix_gas(i)=n_mix_gas(i)
            lw_spectrum%mix_gas_band(i)=mix_gas_band(i)

            IF ( lw_spectrum%n_mix_gas(i)  /=  0 ) THEN
               ii=ii+1
               DO j=1, n_mix_gas(i)
                  lw_spectrum%index_mix_gas(j,ii)=index_mix_gas(j,ii)
               END DO

               DO k=1, i_band_esft(i)
                  DO l=1, npd_mix
                     DO n=1, num_ref_p(i)
                        DO m=1, npd_tmp
                           lw_spectrum%k_mix_gas(n,m,l,k,ii)            &
                              =k_mix_gas(n,m,l,k,ii)
                        END DO
                     END DO
                  END DO
               END DO
            END IF
         END DO

      END IF

! Check error condition
 9999 IF (ierr /= 0) THEN

        CALL ereport(RoutineName, ierr, cmessage)
      END IF

      IF (lhook) CALL dr_hook('SES_INILW',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ses_inilw
