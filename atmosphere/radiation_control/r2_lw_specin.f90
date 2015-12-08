! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to read a longwave spectral namelist.

! Purpose:
!   To read a longwave namelist into a spectral array.

! Method:
!   The spectrum is read into the dynamically allocated array
!   and then reduced to a more manageable size.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!- ---------------------------------------------------------------------
      SUBROUTINE r2_lw_specin(ierr, cmessage                            &
        , lw_spectral_file                                              &
        , l_ch4, l_n2o, l_cfc11, l_cfc12                                &
        , l_cfc113, l_hcfc22, l_hfc125, l_hfc134a                       &
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
        , lw_spectrum                                                   &
        )


!     Modules used:
      USE rad_pcf
      USE dec_spec
      USE filenamelength_mod, ONLY: filenamelength
      USE dust_parameters_mod, ONLY: l_twobin_dust
      USE ereport_mod, ONLY: ereport
      
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE gasid3a_mod, ONLY: npd_gases, ip_h2o, ip_co2, ip_o3, ip_n2o,   &
                             ip_co, ip_ch4, ip_o2, ip_no, ip_so2, ip_no2,&
                             ip_nh3, ip_hno3, ip_n2, ip_cfc11, ip_cfc12, &
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
      CHARACTER  (LEN=filenamelength), INTENT(IN) :: lw_spectral_file
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
      LOGICAL, INTENT(IN) :: l_hcfc22
!           Absorption by hcfc22 is included
      LOGICAL, INTENT(IN) :: l_hfc125
!           Absorption by hfc125 is included
      LOGICAL, INTENT(IN) :: l_hfc134a
!           Absorption by hfc134a is included
      LOGICAL, INTENT(IN) :: l_climat_aerosol
!           Climatological aerosols are to be included
      LOGICAL, INTENT(IN) :: l_sulpc_so2
!           Sulphur cycle is on, so sulphate aerosols can be included, either
!           for their direct effect or on radiation or for radiative diagnostics
      LOGICAL, INTENT(IN) :: l_use_arclsulp
!           The direct effects of sulphate aerosols from NWP climatology
!           are to be included
      LOGICAL, INTENT(IN) :: l_dust
!           Dust aerosols are to be included, either for their
!           direct effect or on radiation or for radiative diagnostics
      LOGICAL, INTENT(IN) :: l_use_arcldust
!           Use the direct radiative effects of mineral dust from
!           NWP climatology in the LW
      LOGICAL, INTENT(IN) :: l_biomass
!           Biomass burning aerosols are to be included, either for their
!           direct effect or on radiation or for radiative diagnostics
      LOGICAL, INTENT(IN) :: l_use_arclbiom
!           Use the direct radiative effects of biomass from
!           NWP climatology in the LW
      LOGICAL, INTENT(IN) :: l_soot
!           Soot aerosols are to be included, either for their
!           direct effect or on radiation or for radiative diagnostics
      LOGICAL, INTENT(IN) :: l_use_arclblck
!           Use the direct radiative effects of black-carbon from NWP
!           climatology in the LW
      LOGICAL, INTENT(IN) :: l_use_seasalt_direct
!           Direct effect of Sea-salt aerosols are to be included
      LOGICAL, INTENT(IN) :: l_use_arclsslt
!           Include the direct radiative effect of sea-salt from NWP
!           climatology in the LW
      LOGICAL, INTENT(IN) :: l_use_biogenic
!           Use the biogenic aerosol direct effect in the LW
      LOGICAL, INTENT(IN) :: l_ocff
!           Fossil-fuel OC aerosols are to be included, either for their
!           direct effect or on radiation or for radiative diagnostics
      LOGICAL, INTENT(IN) :: l_use_arclocff
!           Include the direct radiative effect of fossil-fuel organic
!           carbon aerosol from NWP climatology in the LW
      LOGICAL, INTENT(IN) :: l_use_arcldlta
!           Include the direct radiative effect of delta aerosol
!           from NWP climatology in the LW
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

      INTEGER, INTENT(OUT) :: ierr
!           Error flag
      CHARACTER  (LEN=80), INTENT(OUT) :: cmessage

!     Define the reduced LW spectrum.
      TYPE (spectrum) :: lw_spectrum


!     Local variables.

!     Declare the initial spectrum for the namelist.

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
! LWSP3A declares elements of a spectral file as a namelist.

      NAMELIST/R2LWSP/                                                  &
        ! Blocks Present
     &  L_PRESENT,                                                      &
        ! Block 0
     &  N_BAND, N_ABSORB, N_AEROSOL, TYPE_ABSORB, TYPE_AEROSOL,         &
        ! Block 1
     &  WAVE_LENGTH_SHORT, WAVE_LENGTH_LONG,                            &
        ! Block 2
        solar_flux_band, weight_blue,                                   &
        ! Block 3
     &  RAYLEIGH_COEFFICIENT,                                           &
        ! Block 4
     &  N_BAND_ABSORB, INDEX_ABSORB,                                    &
        ! Block 5
     &  I_BAND_ESFT, I_SCALE_ESFT, I_SCALE_FNC,                         &
     &  P_REFERENCE, T_REFERENCE, K_ESFT, W_ESFT, SCALE_VECTOR,         &
        ! Block 6
     &  N_DEG_FIT, T_REF_PLANCK, THERMAL_COEFFICIENT,                   &
        ! Block 7
     &  I_SPEC_SURFACE, N_DIR_ALBEDO_FIT, L_SURFACE,                    &
     &  SURFACE_ALBEDO, DIRECT_ALBEDO_PARM, EMISSIVITY_GROUND,          &
        ! Block 8
     &  N_BAND_CONTINUUM, INDEX_CONTINUUM, INDEX_WATER,                 &
        ! Block 9
     &  I_SCALE_FNC_CONT, P_REF_CONTINUUM, T_REF_CONTINUUM,             &
     &  K_CONTINUUM, SCALE_CONTINUUM,                                   &
        ! Block 10
     &  I_DROP_PARAMETRIZATION, L_DROP_TYPE, DROP_PARAMETER_LIST,       &
     &  N_DROP_PHF_TERM, DROP_PARM_MIN_DIM, DROP_PARM_MAX_DIM,          &
        ! Block 11
     &  I_AEROSOL_PARAMETRIZATION, NHUMIDITY, L_AEROSOL_SPECIES,        &
     &  AEROSOL_ABSORPTION, AEROSOL_SCATTERING, AEROSOL_ASYMMETRY,      &
     &  N_AEROSOL_PHF_TERM, AEROSOL_PHASE_FNC,                          &
     &  HUMIDITIES,                                                     &
        ! Block 12
     &  I_ICE_PARAMETRIZATION, L_ICE_TYPE, ICE_PARAMETER_LIST,          &
     &  N_ICE_PHF_TERM,ICE_PARM_MIN_DIM, ICE_PARM_MAX_DIM,              &
        ! Block 13
     &  L_DOPPLER_PRESENT, DOPPLER_CORRECTION,                          &
        ! Block 14
     &  N_BAND_EXCLUDE, INDEX_EXCLUDE,                                  &
        ! Block 15
     &  N_AOD_WAVEL, AOD_WAVEL,                                         &
     &  AOD_ABSORPTION, AOD_SCATTERING, I_AOD_TYPE

! LWSP3A end


!     Radiative variables for reducing the spectrum


      INTEGER                                                           &
          ios
!           Status of I/O

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
        , n_aerosol_mixratio
!           Number of aerosols to be allocated in the aerosol mixing
!           ratio arrays

!     Other local variables.

      INTEGER                                                           &
          i, j, k
!           Loop variable

      CHARACTER                                                         &
           ch_ios*5
!           Character string for iostat error

      CHARACTER (LEN=12), PARAMETER :: RoutineName='r2_lw_specin'
      
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_in,zhook_handle)

!     Each block is initialized as missing:
      l_present(:)   = .FALSE.
      l_surface(:)   = .FALSE.
      l_drop_type(:) = .FALSE.
      l_aerosol_species(:) = .FALSE.
      l_ice_type(:)        = .FALSE.
      l_doppler_present(:) = .FALSE.


! Some spectral files specify the aerosol_asymmetry and others specify
! the aerosol phase-function. We want to read in both types of spectral
! file and hence we need to find out whether the namelists contains
! one or the other. In order to do this we initialise aerosol_asymmetry
! and aerosol_phase_fnc to the missing data flag and check which one
! has been allocate proper values. If aerosol_phase_fnc is specified
! nothing else needs to be done, otherwise we copy aerosol_asymmetry
! into aerosol_phase_fnc.

      aerosol_asymmetry = rmdi
      nhumidity=1

! Initialise the number of terms in the phase function for the case
! where these are not present in the spectral file
      n_ice_phf_term=1
      n_drop_phf_term=1

! Initialise the AOD wavelengths in case these are not present in the
! spectral file
      aod_wavel(1:6) = (/ 380., 440., 550., 670., 865., 1020. /) 
! wl above are in nanometers, convert to SI:
      aod_wavel(1:6) = aod_wavel(1:6) * 1.0e-9

! Initialise weight_blue in case it is not present
      weight_blue = rmdi

!     Read the longwave spectrum as a namelist.
      WRITE(6,'(a, a)') 'Namelist file: ',TRIM(ADJUSTL(lw_spectral_file))
      OPEN(UNIT=80, FILE=lw_spectral_file, IOSTAT=ios)
      IF (ios /= 0) THEN
        ierr=i_err_io
        WRITE(ch_ios, '(i5)') ios
        cmessage='error opening longwave spectral file.'                &
          //' iostat='//ch_ios
        IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
        RETURN
      END IF
      READ(80, r2lwsp)
      CLOSE(80)


!     Test for minimal requisite information.
      IF ( .NOT.(l_present(0).AND.                                      &
                 l_present(6) ) ) THEN
        cmessage='longwave spectrum is deficient.'
        ierr=i_err_fatal
        IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
        RETURN
      END IF

      IF (l_gas) THEN
        IF (.NOT.l_present(5)) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** warning: the lw spectrum contains no '                 &
            //'gaseous absorption data.'
          ierr = i_warning
          CALL ereport(RoutineName, ierr, cmessage)  
          l_gas=.FALSE.
        END IF
      END IF

      IF (l_continuum) THEN
        IF (.NOT.l_present(9)) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** warning: the lw spectrum contains no '                 &
            //'continuum absorption data.'
          ierr = i_warning
          CALL ereport(RoutineName, ierr, cmessage)  
          l_continuum=.FALSE.
        END IF
      END IF

      IF (l_drop) THEN
        IF (.NOT.l_present(10)) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** warning: the lw spectrum contains no '                 &
            //'data for water droplets.'
          ierr = i_warning
          CALL ereport(RoutineName, ierr, cmessage)  
          l_drop=.FALSE.
        END IF
      END IF

      IF (l_aerosol) THEN
        IF (.NOT.l_present(11)) THEN
!         No warning given as this is expected for the incremental
!         radiative timestepping scheme
          l_aerosol=.FALSE.
        END IF
      END IF

      IF (l_ice) THEN
        IF (.NOT.l_present(12)) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** warning: the lw spectrum contains no '                 &
            //'data for ice crystals.'
          ierr = i_warning
          CALL ereport(RoutineName, ierr, cmessage)  
          l_ice=.FALSE.
        END IF
      END IF

!     Set reduced dimensions, either from the sizes of the fixed arrays
!     or from the arrays read in.

      lw_spectrum%npd_type=npd_type
      lw_spectrum%npd_band=MAX(n_band, 1)
      lw_spectrum%npd_species=MAX(n_absorb, 1)
      lw_spectrum%npd_scale_variable=npd_scale_variable
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
          lw_spectrum%npd_exclude                                       &
            =MAX(lw_spectrum%npd_exclude, n_band_exclude(i))
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
          WRITE(cmessage, '(a)')                                        &
            '*** warning: water vapour is not included in the '         &
            // 'longwave spectral file.'
          ierr = i_warning
          CALL ereport(RoutineName, ierr, cmessage)  
        END IF

        IF (.NOT.l_gas_included(ip_co2)) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** warning: carbon dioxide is not included in the '       &
            // 'longwave spectral file.'
          ierr = i_warning
          CALL ereport(RoutineName, ierr, cmessage)  
        END IF

        IF (.NOT.l_gas_included(ip_o3)) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** warning: ozone is not included in the '                &
            // 'longwave spectral file.'
          ierr = i_warning
          CALL ereport(RoutineName, ierr, cmessage)  
        END IF

        IF ((.NOT.l_gas_included(ip_ch4)).AND.l_ch4) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** error: methane is not included in the longwave '       &
            // 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN
        END IF

        IF ((.NOT.l_gas_included(ip_n2o)).AND.l_n2o) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** error: nitrous oxide is not included in the longwave ' &
            // 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN
        END IF

        IF ((.NOT.l_gas_included(ip_cfc11)).AND.l_cfc11) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** error: cfc11 is not included in the longwave '         &
            // 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN
        END IF

        IF ((.NOT.l_gas_included(ip_cfc12)).AND.l_cfc12) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** error: cfc12 is not included in the longwave '         &
            // 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN
        END IF

        IF ((.NOT.l_gas_included(ip_cfc113)).AND.l_cfc113) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** error: cfc113 is not included in the longwave '        &
            // 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN
        END IF

        IF ((.NOT.l_gas_included(ip_hcfc22)).AND.l_hcfc22) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** error: hcfc22 is not included in the longwave '        &
            // 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN
        END IF

        IF ((.NOT.l_gas_included(ip_hfc125)).AND.l_hfc125) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** error: hfc125 is not included in the longwave '        &
            // 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN
        END IF

        IF ((.NOT.l_gas_included(ip_hfc134a)).AND.l_hfc134a) THEN
          WRITE(cmessage, '(a)')                                        &
            '*** error: hfc134a is not included in the longwave '       &
            // 'spectral file, but was requested in the run.'
          ierr=i_err_fatal
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN
        END IF


!     Set an appropriate reduced dimension.
      lw_spectrum%npd_species=MAX(n_absorb_retain, 1)

      lw_spectrum%npd_esft_term=1
      IF (l_present(5)) THEN
        DO i=1, n_band
          DO j=1, n_band_absorb(i)
            IF (l_retain_absorb(index_absorb(j, i)))                    &
              lw_spectrum%npd_esft_term                                 &
                =MAX(lw_spectrum%npd_esft_term                          &
                , i_band_esft(i, index_absorb(j, i)))
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
            cmessage='the lw spectral file lacks some '                 &
              //'climatological aerosols.'
            IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
            RETURN

          END IF
        ELSE

          ierr=i_err_fatal
          cmessage='lw spectral file contains no aerosol data.'
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN

        END IF

      END IF ! ends land/sea climatological aerosol, or murk in radiation

!     Check the spectral file for mineral dust classes.
!     (Only required for the direct effect).

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
                  //'mineral dust aerosol for two-bin scheme.'
               IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
               RETURN

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
                  //'mineral dust aerosol for six-bin scheme.'
               IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
               RETURN

            END IF

           END IF ! ends test on 6 bin dust being required

         ELSE

            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
            RETURN

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
            cmessage='the lw spectral file lacks some '                 &
              //'soot aerosol data.'
            IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
            RETURN

          END IF

        ELSE


          ierr=i_err_fatal
          cmessage='lw spectral file contains no soot data.'
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN

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
               IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
               RETURN

            END IF

         ELSE

            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
            RETURN

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
            cmessage='the lw spectral file lacks some '                 &
              //'sulphate aerosols.'
              IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
              RETURN

          END IF
        ELSE

          ierr=i_err_fatal
          cmessage='lw spectral file contains no aerosol data.'
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN

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
            cmessage='the lw spectral file lacks some '                 &
              //'seasalt aerosols.'
              IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
              RETURN

          END IF
        ELSE

          ierr=i_err_fatal
          cmessage='lw spectral file contains no aerosol data.'
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN

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
                  IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
                  RETURN

            END IF

          ELSE

            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
            RETURN

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
                  IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
                  RETURN

            END IF

         ELSE

            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
            RETURN

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
                  IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
                  RETURN

            END IF

          ELSE

            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
            RETURN

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
                  IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
                  RETURN

            END IF

          ELSE

            ierr=i_err_fatal
            cmessage='LW Spectral file contains no aerosol data.'
            IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
            RETURN

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

        IF (ANY(aerosol_asymmetry /= rmdi)) THEN
           aerosol_phase_fnc(:,1,:,:)=aerosol_asymmetry
        END IF

      END IF


!     Check that block 15 is present if the aerosol optical depth
!     was requested.
      IF (l_use_aod.AND.l_aerosol) THEN
        IF(.NOT. l_present(15)) THEN
          ierr = i_err_fatal
          cmessage='Block 15 needed in the LW spectral file.'
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN
        END IF
       ! Check that the number of wavelengths in the spectral
       ! file is not larger than the number set in the include
       ! file MXSIZE3A.
        IF(n_aod_wavel  >   npd_aod_wavel) THEN
          ierr = i_err_fatal
          cmessage='Increase NPD_AOD_WAVEL in MXSIZE3A.'
          IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
          RETURN
        END IF
        lw_spectrum%npd_aod_wavel = npd_aod_wavel
      END IF


!     Allocate space for the LW spectrum.
! DEPENDS ON: r2_allocate_spectrum
      CALL r2_allocate_spectrum(lw_spectrum, .FALSE.)

!     Transfer the large namelist to the reduced spectrum.


! DEPENDS ON: r2_compress_spectrum
      CALL r2_compress_spectrum(                                        &
!                       Spectral array in namelist
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
        , lw_spectrum                                                   &
        )

        IF (lhook) CALL dr_hook('R2_LW_SPECIN',zhook_out,zhook_handle)
        RETURN
      END SUBROUTINE r2_lw_specin
