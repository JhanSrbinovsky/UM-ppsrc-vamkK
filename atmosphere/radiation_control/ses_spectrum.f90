! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!     ******************************************************************
!     *                                                                *
!     *            MODULE TO READ IN THE SPECTRAL DATA FILE.           *
!     *            TO RADIATION CODE (INTERPOLATION VERSION)           *
!     *                                                                *
!     ******************************************************************
!     ******************************************************************
!     *                                                                *
!     *     TOP LEVEL SUBROUTINE TO CONTROL READING OF THE SPECTRAL    *
!     *                        DATA FILE.                              *
!     *                                                                *
!     ******************************************************************
!
!     Code Owner: See Unified Model Code Owners HTML page
!     This file belongs in section: Radiation Control

SUBROUTINE ses_spectrum(ierr                                            &
         , file_spectral                                                &
         , l_present                                                    &
         , n_band                                                       &
         , n_absorb, n_band_absorb, index_absorb, type_absorb           &
         , n_mix_gas, index_mix_gas, mix_gas_band                       &
         , n_aerosol, type_aerosol                                      &
         , n_band_continuum, index_continuum, index_water               &
         , i_band_esft, num_ref_p, num_mix                              &
         , f_mix, k_esft, k_mix_gas,  w_esft                            &
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


      USE rad_pcf, ONLY: i_normal, i_err_fatal, i_err_exist
      USE filenamelength_mod, ONLY:                                     &
          filenamelength
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport

      USE gasid3a_mod, ONLY: npd_gases, ip_h2o, ip_co2, ip_o3, ip_n2o,   &
                             ip_co, ip_ch4, ip_o2, ip_no, ip_so2, ip_no2,&
                             ip_nh3, ip_hno3, ip_n2, ip_cfc11, ip_cfc12, &
                             ip_cfc113, ip_hcfc22, ip_hfc125, ip_hfc134a,&
                             ip_cfc114
      IMPLICIT NONE
          
!     INCLUDE HEADER FILES.
!     ------------------------------------------------------------------
!     MODULE SETTING UNIT NUMBERS FOR READING SPECTRAL FILE.

      INTEGER                                                           &
           iu_spc                                                       &
!             UNIT NUMBER TO READ SPECTRAL FILE
         , iu_spc1                                                      &
!             Unit number to read gas absorption coefficients
         , iu_spc2                                                      &
!             Unit number to read minor gas absorption coefficients
         , iu_spc3

      PARAMETER(                                                        &
           iu_spc=91, iu_spc1=iu_spc+1, iu_spc2=iu_spc+2                &
         , iu_spc3=iu_spc2+1                                            &
         )

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

!     DUMMY VARIABLES.
      CHARACTER                                                         &
                !, INTENT(IN)
           file_spectral*(*)
!             SPECTRAL FILE

      INTEGER, INTENT(INOUT) :: ierr  ! Error flag

      INTEGER                                                           &
              !, INTENT(OUT)
           n_parameter_water                                            &
!             Number of parameters for water cloud optical property scheme
         , n_parameter_ice
!             Number of parameters for ice cloud optical property scheme

!     LOCAL VARIABLES.

      CHARACTER (LEN=*), PARAMETER ::  RoutineName = 'ses_spectrum'
      CHARACTER (LEN=256) :: cmessage

      CHARACTER                                                         &
           line*80                                                      &
!             LINE READ FROM FILE
         , char_dum*80                                                  &
!             DUMMY CHARCATER VARIABLE
         , char_in*1
!             input character

      CHARACTER(LEN=filenamelength) :: spectral_k   
!             spectral absorption data


      INTEGER                                                           &
           ios                                                          &
!             IO VARIABLE
         , i_type                                                       &
!             TYPE OF BLOCK READ IN
         , i_subtype                                                    &
!             SUBTYPE OF BLOCK
         , i_version                                                    &
!             VERSION FOR TYPE AND SUBTYPE
         , i
!             LOOP VARIABLE
      LOGICAL                                                           &
           l_exist                                                      &
!             EXISTENCE FLAG FOR FILE
         , l_open
!             OPEN FLAG FOR FILE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('SES_SPECTRUM',zhook_in,zhook_handle)

! Initialize array

      DO i=1, npd_drop_type
         i_drop_parametrization(i) = 0
      END DO
      n_parameter_water=0
      n_parameter_ice=0

!     It is important to know which gas is water vapour for some
!     applications: here we initialize index_water to 0 to
!     produce an error in the code if it is not reset to a legal value.
!     This should guard against omitting water vapour when it is needed.
      index_water=0
      i=INDEX(file_spectral, ' ') - 1
      spectral_k = file_spectral(1:i) // '_k'

!     CHECK THAT THE SPECTRAL FILE EXISTS.
      INQUIRE(FILE=file_spectral, EXIST=l_exist)
      IF (.NOT.l_exist) THEN
         cmessage = file_spectral(1:i) // ' FILE DOES NOT EXIST.'
         ierr=i_err_exist
         GO TO 9999
      END IF
      INQUIRE(FILE=spectral_k, EXIST=l_exist)
      IF (.NOT.l_exist) THEN
         cmessage = spectral_k // ' FILE DOES NOT EXIST.'
         ierr=i_err_exist
         GO TO 9999
      END IF
      INQUIRE(iu_spc, OPENED=l_open)
      IF (l_open) THEN
         cmessage = 'UNIT IS ALREADY CONNECTED.'
         ierr=i_err_fatal
         GO TO 9999
      END IF

!     OPEN THE FILE FOR READING
      OPEN(UNIT=iu_spc, FILE=file_spectral, IOSTAT=ios                  &
         , STATUS='OLD')
      IF (ios /= 0) THEN
         cmessage = file_spectral(1:i) // ' FILE COULD NOT BE OPENED.'
         ierr=i_err_fatal
         GO TO 9999
      END IF
      OPEN(UNIT=iu_spc1, FILE=spectral_k, IOSTAT=ios                    &
         , STATUS='OLD')
      IF (ios /= 0) THEN
         cmessage = spectral_k // ' FILE COULD NOT BE OPENED.'
         ierr=i_err_fatal
         GO TO 9999
      END IF

! Skip over head of the file containing absorption coefficients

110   READ(iu_spc1, '(a1)') char_in
      IF ( char_in  /=  '*' ) GO TO 110

!     INITIALIZE THE ELEMENTS OF THE ARRAYS L_PRESENT
!     L_SURFACE AND L_AERSOL_SPECIES TO FALSE.
      DO i=0, npd_type
         l_present(i)=.FALSE.
      END DO
      DO i=1, npd_surface
         l_surface(i)=.FALSE.
      END DO
      DO i=1, npd_aerosol_species
         l_aerosol_species(i)=.FALSE.
      END DO
      DO i=1, npd_drop_type
         l_drop_type(i)=.FALSE.
      END DO


!     READ THROUGH THE FILE PROCESSING THE BLOCKS OF DATA AS THEY
!     ARE ENCOUNTERED.
!     EACH LINE IS READ INTO AN INTERNAL FILE AND THEN PROCESSED.
100   READ(iu_spc, '(A80)', END=101) line
!     LOCATE A BLOCK HEADER.
      IF (line(1:6) == '*BLOCK') THEN
         READ(line, FMT='(A15, I4, 12X, I4, 12X, I4)', IOSTAT=ios)      &
            char_dum, i_type, i_subtype, i_version
         IF (ios /= 0) THEN
            cmessage = 'BLOCK HEADER IS INCORRECT.'
            ierr=i_err_fatal
            GO TO 9999
         END IF

!        READ IN THE REST OF THE BLOCK.
         CALL ses_block
         IF (ierr /= i_normal) THEN
            GO TO 9999
         END IF

!        READ IN THE TERMINATION STATEMENT
         READ(iu_spc, '(A4)') char_dum
         IF (char_dum(1:4) /= '*END') THEN
            WRITE(cmessage, '(2X, A37, I5, A13, I5, A12, I5)')          &
               'BLOCK INCORRECTLY TERMINATED. TYPE = ', i_type          &
               , ': SUB-TYPE = ', i_subtype                             &
               , ': VERSION = ', i_version
            ierr=i_err_fatal
            GO TO 9999
         END IF

!        THE BLOCK HAS BEEN PROPERLY READ: RECORD THE DATA AS PRESENT.
         l_present(i_type)=.TRUE.

      END IF

!     READ NEXT LINE.
      GO TO 100


!     CLOSE THE FILE.
101   CLOSE(iu_spc)
      CLOSE(iu_spc1)

! Check error condition
 9999 IF (ierr /= i_normal) THEN
        CALL ereport(RoutineName, ierr, cmessage)
      END IF

      IF (lhook) CALL dr_hook('SES_SPECTRUM',zhook_out,zhook_handle)
      RETURN


CONTAINS

!     *****************************************************************
!     *                                                               *
!     *      SUBROUTINE TO CONTROL THE READING IN OF BLOCKS OF DATA.  *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block

!     Depending on the value of i_type, the appropriate subroutine
!     is called.
      IF (i_type == 0) THEN
         IF (i_subtype == 0) THEN
            IF (i_version == 1) THEN
               CALL ses_block_0_0_1
            END IF
         END IF
      ELSE IF (i_type == 1) THEN
         IF (i_subtype == 0) THEN
            IF (i_version == 0) THEN
               CALL ses_block_1_0_0
            END IF
         END IF
      ELSE IF (i_type == 2) THEN
         IF (i_subtype == 1) THEN
            IF (i_version == 0) THEN
               CALL ses_block_2_1_0
            END IF
         END IF
      ELSE IF (i_type == 3) THEN
         IF (i_subtype == 0) THEN
            IF (i_version == 0) THEN
               CALL ses_block_3_0_0
            END IF
         END IF
      ELSE IF (i_type == 4) THEN
         IF (i_subtype == 0) THEN
            IF (i_version == 1) THEN
               CALL ses_block_4_0_1
            END IF
         END IF
      ELSE IF (i_type == 5) THEN
         IF (i_subtype == 1) THEN
            IF (i_version == 0) THEN
               CALL ses_block_5_1_0
            END IF
         END IF
      ELSE IF (i_type == 6) THEN
         IF (i_subtype == 0) THEN
           CALL read_block_6_0_0
         ELSE IF (i_subtype == 1) THEN
           CALL ses_block_6_1_0
         END IF
      ELSE IF (i_type == 8) THEN
         IF (i_subtype == 0) THEN
            IF (i_version == 0) THEN
               CALL ses_block_8_0_0
            END IF
         END IF
      ELSE IF (i_type == 9) THEN
         IF (i_subtype == 1) THEN
            IF (i_version == 0) THEN
               CALL ses_block_9_1_0
            END IF
         END IF
      ELSE IF (i_type == 10) THEN
         IF (i_subtype == 0) THEN
            IF (i_version == 2) THEN
               CALL ses_block_10_0_2
            END IF
         END IF
      ELSE IF (i_type == 11) THEN
         IF (i_subtype == 0) THEN

!           This is for water insolable aerosols
            IF (i_version == 1) THEN
               CALL ses_block_11_0_1
            ELSE IF (i_version == 2) THEN
               CALL ses_block_11_0_2
            END IF
         ELSE IF (i_subtype == 1) THEN

!           This is for moisture aerosols
            IF (i_version == 0) THEN
               CALL ses_block_11_1_0
            ELSE IF (i_version == 2) THEN
               CALL ses_block_11_1_2
            END IF
         END IF
      ELSE IF (i_type == 12) THEN
         IF (i_subtype == 0) THEN
            IF (i_version == 0) THEN
               CALL ses_block_12_0_0
            ELSE IF (i_version == 2) THEN
               CALL read_block_12_0_2
            END IF
         END IF
      ELSE IF (i_type == 14) THEN
         IF (i_subtype == 0) THEN
            IF (i_version == 0) THEN
               CALL ses_block_14_0_0
            END IF
         END IF
      ELSE IF (i_type == 15) THEN
         IF (i_subtype == 0) THEN
            CALL ses_block_15_0_0
         ELSE IF(i_subtype == 1 ) THEN
            CALL ses_block_15_1_0
         END IF
      ELSE IF (i_type == 16) THEN
         IF (i_subtype == 0) THEN
            IF (i_version == 0) THEN
               CALL ses_block_16_0_0
            END IF
         END IF
      ELSE
!        The value of i_type does not correspond to a supported type
!        of block.
         WRITE(cmessage, '(A18, I4, A32)') 'VALUE OF I_TYPE = ',i_type  &
            , ': THIS VALUE IS NOT IMPLEMENTED.'
         ierr=i_err_fatal
         RETURN
      END IF

      END SUBROUTINE ses_block

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 0.       *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_0_0_1

!     Local variables.
      CHARACTER :: chdum
!             Dummy character
      INTEGER :: idum
!             Dummy integer

!     Skip over the header.
      READ(iu_spc, *)

!     Read in the number of spectral intervals, the number of
!     gaseous absorbers and the number of aerosols.
      READ(iu_spc, '(27X, I5)', IOSTAT=ios) n_band
      IF (ios /= 0) THEN
         cmessage =  '*** Error in subroutine ses_block_0_0_1. ' //     &
           'Number of bands could not be read.'
         ierr=i_err_fatal
         RETURN
      END IF
      IF (n_band >  npd_band) THEN
         cmessage =  '*** Error in subroutine ses_block_0_0_1. ' //     &
           'Number of bands exceeds maximum permitted number. ' //      &
           'Increase npd_band and recompile.'
         ierr=i_err_fatal
         RETURN
      END IF
      READ(iu_spc, '(36X, I5)', IOSTAT=ios) n_absorb
      IF (ios /= 0) THEN
         cmessage =  '*** Error in subroutine ses_block_0_0_1. ' //     &
           'Number of absorbers could not be read.'
         ierr=i_err_fatal
         RETURN
      END IF
      IF (n_absorb >  npd_species) THEN
         cmessage =  '*** Error in subroutine ses_block_0_0_1. ' //     &
           'Number of gaseous absorbers exceeds maximum permitted ' //  &
           'number. Increase npd_species and recompile.'
         ierr=i_err_fatal
         RETURN
      END IF
      READ(iu_spc, '(27X, I5)', IOSTAT=ios) n_aerosol
      IF (ios /= 0) THEN
         cmessage =  '*** Error in subroutine ses_block_0_0_1. ' //     &
           'Number of aerosols could not be read.'
         ierr=i_err_fatal
         RETURN
      END IF
      IF (n_aerosol >  npd_aerosol_species) THEN
         cmessage =  '*** Error in subroutine ses_block_0_0_1. ' //     &
           'Number of aerosols exceeds maximum permitted number. ' //   &
           'Increase npd_aerosol_species and recompile.'
         ierr=i_err_fatal
         RETURN
      END IF

!     Read over the headers and the list of absorbers.
      READ(iu_spc, '(/)')
      DO i=1, n_absorb
         READ(iu_spc, '(I5, 7X, I5, 7X, A)') idum, type_absorb(i), chdum
      END DO

!     Read over the headers and the list of aerosols.
      READ(iu_spc, '(/)')
      DO i=1, n_aerosol
         READ(iu_spc, '(I5, 7X, I5, 7X, A)')                            &
            idum, type_aerosol(i), chdum
      END DO

      END SUBROUTINE ses_block_0_0_1

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 1.       *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_1_0_0

!     Local variables.
      INTEGER :: idum
!             Dummy integer

!     Skip over the headers.
      READ(iu_spc, '(//)')

!     Read in the limits on the intervals in the spectrum
      DO i=1, n_band
         READ(iu_spc, FMT='(I5, 7X, 1PE16.9, 4X, 1PE16.9)'              &
            , IOSTAT=ios)                                               &
            idum, wave_length_short(i), wave_length_long(i)
         IF (ios /= 0) THEN
            cmessage =  '*** Error in subroutine ses_block_1_0_0. ' //  &
              'Limits on wavelength could not be read.'
            ierr=i_err_fatal
            RETURN
         END IF
      END DO

      END SUBROUTINE ses_block_1_0_0

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 2.       *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_2_1_0

!     Local variables.
      INTEGER :: is, k
!       Loop variables
      INTEGER :: i_band
!       Dummy variable
      INTEGER :: n_term
!       Number of terms


!     Read in the limits on the intervals in the spectrum
      DO is=1, npd_solar_src
!       Skip over the headers.
        READ(iu_spc, '(/)')
        DO i=1, n_band
          READ(iu_spc, FMT='(I5, 7X, I5)', IOSTAT=ios) i_band, n_term

          READ(iu_spc, FMT='(4(1X,1PE16.9))',IOSTAT=ios)                &
            (solar_flux_band(k, i, is),k=1, n_term)
          IF (ios /= 0) THEN
            cmessage = '*** Error in subroutine ses_block_2_1_0. ' //   &
              'Solar spectral data is not correct.'
            ierr=i_err_fatal
            RETURN
          END IF
        END DO
      END DO

      END SUBROUTINE ses_block_2_1_0

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 3.       *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_3_0_0

!     Local variables.
      INTEGER :: idum
!             Dummy integer

!     Skip over the headers.
      READ(iu_spc, '(//)')

!     Read in the limits on the intervals in the spectrum
      DO i=1, n_band
         READ(iu_spc, FMT='(I5, 7X, 1PE16.9)', IOSTAT=ios)              &
            idum, rayleigh_coefficient(i)
         IF (ios /= 0) THEN
            cmessage =  '*** Error in subroutine ses_block_3_0_0. ' //  &
              'Rayleigh scattering data is not correct.'
            ierr=i_err_fatal
            RETURN
         END IF
      END DO

      END SUBROUTINE ses_block_3_0_0

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 4.       *
!     *          INDEX NUMBERS OF ABSORBERS IN EACH BAND.             *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_4_0_1

!     Local variables.
      INTEGER                                                           &
           idum                                                         &
!             Dummy integer
         , j
!             Loop variable


!     Skip over the headers.
      READ(iu_spc, '(////)')

!     Read in the list of absorbers in each band.
      DO i=1, n_band
         READ(iu_spc, FMT='(I5, 7X, I5)', IOSTAT=ios)                   &
            idum, n_band_absorb(i)
         IF (ios /= 0) THEN
            cmessage =                                                  &
               '*** Error in subroutine ses_block_4_0_1. ' //           &
               'The list of absorbers is not correct.'
            ierr=i_err_fatal
            RETURN
         END IF
         IF (n_band_absorb(i) >  0) THEN
            READ(iu_spc, FMT='(5X, 8(2X, I3))', IOSTAT=ios)             &
               ( index_absorb(j, i), j=1, n_band_absorb(i) )
         END IF
         IF (ios /= 0) THEN
            cmessage =                                                  &
               '*** Error in subroutine ses_block_4_0_1. ' //           &
               'The index of absorbers is not correct.'
            ierr=i_err_fatal
            RETURN
         END IF
      END DO

      END SUBROUTINE ses_block_4_0_1

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 5.       *
!     *            EXPONENTIAL SUM FITTING TERMS.                     *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_5_1_0

!     Local variables.
      INTEGER                                                           &
           idum_band                                                    &
!             Dummy integer
         , number_term                                                  &
!             Number of terms
         , j, k                                                         &
!             Loop variables
         , it, ip, i_gas, n_gas, i_term                                 &
         , index_band(npd_band)

      READ(iu_spc, '(/)')

!     Read in the number of ESFT terms in each band.
      DO i = 1 , n_band
         READ(iu_spc, FMT='(4I5,   1x,1PE16.9 )', IOSTAT=ios)           &
            index_band(i), i_band_esft(i), num_ref_p(i), num_mix(i)     &
          , f_mix(i)
         IF (ios /= 0) THEN
            cmessage = '*** 1st error in subroutine ses_block_5_1_0.'// &
              ' ESFT data is not consistent with summary.'
            ierr=i_err_fatal
            RETURN
         END IF

         IF(i_band_esft(i) >  npd_esft_term) THEN
            cmessage = '*** Error in subroutine ses_block_5_1_0. ' //   &
              'Too many esft terms have been given. ' //                &
              'Increase npd_esft_term and recompile.'
            ierr=i_err_fatal
            RETURN
         END IF
!        For each band read the values of the coefficients.
         READ(iu_spc, FMT='(4(1X,1PE16.9))', IOSTAT=ios)                &
             (w_esft(k, i),k=1, i_band_esft(i))

         IF (ios /= 0) THEN
            cmessage = '*** 2nd error in subroutine ses_block_5_1_0.'// &
              ' ESFT data is not consistent with summary.'
            ierr=i_err_fatal
            RETURN
         END IF
      END DO

      DO i=1, n_band
!        Skip over the headers.
         READ(iu_spc1, '(/)')
         DO i_gas=1, n_band_absorb(i)
            DO i_term=1, i_band_esft(index_band(i))
               DO ip=1, num_ref_p(i)
                  READ(iu_spc1, FMT='(5(1PE12.6,1x))', IOSTAT=ios)      &
                  (k_esft(ip,it,i_gas,i_term,i),it=1, npd_tmp)
                  IF (ios /= 0) THEN
                    WRITE(cmessage, '(A45, 4i4)')                       &
                      '*** 4th error in subroutine ses_block_5_1_0. ',  &
                      i,i_term,i_gas,ip
                    ierr=i_err_fatal
                    RETURN
                  END IF
               END DO
            END DO
         END DO
      END DO

      END SUBROUTINE ses_block_5_1_0

!     ******************************************************************
!     *                                                                *
!     *      Subroutine to read the contents of block 6.               *
!     *      Data for the thermal source function in each band.        *
!     *                                                                *
!     ******************************************************************
      SUBROUTINE read_block_6_0_0

!     Local variables.
      INTEGER ::                                                        &
          k                                                             &
!           Loop variable
        , i_band
!           Number of band

      READ(iu_spc, '(/, 23x, i5, 26x, 1pe16.9)') n_deg_fit, t_ref_planck
      IF (n_deg_fit > npd_thermal_coeff-1) THEN
        cmessage = '*** Error in subroutine read_block_6_0_0. ' //      &
          'The degree of the polynomial fit is too high. ' //           &
          'Increase npd_thermal_coeff and recompile.'
        ierr=i_err_fatal
        RETURN
      ENDIF

      READ(iu_spc, '(/)')
      DO i=1, n_band
        READ(iu_spc, '(i5, 7x, (t13, 3(1pe16.9, 4x)))', iostat=ios)     &
          i_band, (thermal_coefficient(k, i), k=0, n_deg_fit)
        IF (ios /= 0) THEN
          cmessage = '*** Error in subroutine read_block_6_0_0. ' //    &
           'The data for the thermal source function could not be read.'
          ierr=i_err_fatal
          RETURN
        ENDIF
      ENDDO

      END SUBROUTINE read_block_6_0_0

!     ******************************************************************
!     *                                                                *
!     *      SUBROUTINE TO READ THE CONTENTS OF BLOCK 6.               *
!     *         COEFFICIENTS FOR SOURCE FUNCTION IN EACH BAND.         *
!     *                                                                *
!     ******************************************************************
      SUBROUTINE ses_block_6_1_0

!     LOCAL VARIABLES.
      INTEGER                                                           &
           k                                                            &
!             LOOP VARIABLE
         , i_band                                                       &
!             NUMBER OF BAND
         , n_term
!             NUMBER OF FRACTION TERM

      READ(iu_spc, '(/)')
      DO i=1, n_band
         READ(iu_spc, '(I5, 4(7X, I5))', IOSTAT=ios)                    &
            i_band, n_term
         READ(iu_spc, '(4(1X,1PE16.9))', IOSTAT=ios)                    &
            (plk_frac(k, i),k=1, n_term)
         IF (ios /= 0) THEN
            cmessage = '*** Error in subroutine ses_block_6_1_0. ' //   &
              'The Planck fraction data could not be read.'
            ierr=i_err_fatal
            RETURN
         END IF
         READ(iu_spc1, *)
         READ(iu_spc1, '(5(1PE12.6,1x))', IOSTAT=ios)                   &
             (planckb_ref(k, i), k=1,161)
         IF(ios /= 0) THEN
            cmessage = '*** Error in subroutine ses_block_6_1_0. ' //   &
              'The reference Planck data error'
            ierr=i_err_fatal
            RETURN
         END IF
      END DO

      END SUBROUTINE ses_block_6_1_0

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 8.       *
!     *           INDEX NUMBERS OF CONTINUA IN EACH BAND.             *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_8_0_0

!     LOCAL VARIABLES.
      INTEGER                                                           &
           idum                                                         &
!             DUMMY INTEGER
         , j
!             LOOP VARIABLE

!     Skip over the headers.
      READ(iu_spc, '(////)')

!     Read in the limits on the intervals in the spectrum
      DO i=1, n_band
         READ(iu_spc, FMT='(I5, 7X, I5)', IOSTAT=ios)                   &
            idum, n_band_continuum(i)
         IF (ios /= 0) THEN
            cmessage = '*** Error in subroutine ses_block_8_0_0. ' //   &
              'The list of continua is not correct.'
            ierr=i_err_fatal
            RETURN
         END IF
         IF (n_band_continuum(i) >  npd_continuum) THEN
            cmessage = '*** Error in subroutine ses_block_8_0_0. ' //   &
              'There are too many continua: ' //                        &
              'increase npd_continuum and recompile.'
            ierr=i_err_fatal
            RETURN
         END IF
         IF (n_band_continuum(i) >  0) THEN
            READ(iu_spc, FMT='(5X, 4(2X, I3))', IOSTAT=ios)             &
               ( index_continuum(i, j), j=1, n_band_continuum(i) )
         END IF
         IF (ios /= 0) THEN
            cmessage = '*** Error in subroutine ses_block_8_0_0. ' //   &
              'The list of continua is not correct.'
            ierr=i_err_fatal
            RETURN
         END IF
      END DO

!     Read in the indices of gases forming the continuum species.
      READ(iu_spc, '(/, 22X, I5)') index_water

      END SUBROUTINE ses_block_8_0_0

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 9.       *
!     *            EXPONENTIAL SUM FITTING TERMS.                     *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_9_1_0

!     LOCAL VARIABLES.
      INTEGER                                                           &
           idum_band                                                    &
!             DUMMY INTEGER
         , idum_continuum                                               &
!             DUMMY INTEGER
         , n_term                                                       &
!             DUMMY INTEGER
         , j, ip, it                                                    &
!             LOOP VARIABLE
         , i_tmp, num_t                                                 &
         , i_term

!     Skip over the headers.
      READ(iu_spc, '(/)')

      DO i=1, n_band
         DO j=1, n_band_continuum(i)
            READ(iu_spc, *, IOSTAT=ios)                                 &
               idum_band, idum_continuum, n_term, num_t
            IF (ios /= 0) THEN
               cmessage = '*** 1st error in subroutine ses_block_9_1_0.'&
                // ' Continuum type or band could not be read.'
               ierr=i_err_fatal
               RETURN
            END IF

            IF (idum_continuum == 1) THEN
               READ(iu_spc1, '(/)')
               DO i_term=1, n_term
                  DO ip=1, 21
                     READ(iu_spc1, '(5(1PE12.6,1x))', IOSTAT=ios)       &
                     (k_h2oc(ip,it,i_term,i),it=1, npd_tmp)
                     IF (ios /= 0) THEN
                       WRITE(cmessage, '(A45, 3i4)')                    &
                         '*** 2nd error in subroutine ses_block_9_1_0. '&
                         ,i,i_term,ip
                       ierr=i_err_fatal
                       RETURN
                     END IF
                  END DO
               END DO
            ELSE
               DO i_term=1, n_term
                  READ(iu_spc, '(6X, 5(1x, 1PE16.9))', IOSTAT=ios)      &
                    (k_continuum(i_term, i_tmp, idum_band, j)           &
                     , i_tmp=1, num_t)
               END DO
            END IF
         END DO
      END DO

      END SUBROUTINE ses_block_9_1_0

!     *****************************************************************
!     *                                                               *
!     *     Subroutine to read in the contents of block type 10.      *
!     *           Parametrized droplet scattering data.               *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_10_0_2

      USE rad_pcf, ONLY:                                                &
        ip_slingo_schrecker, ip_ackerman_stephens, ip_drop_pade_2

!     Local variables.
      INTEGER ::                                                        &
          i_drop                                                        &
!           Type of droplet
        , i_parametrization_drop                                        &
!           Dummy index of parameter scheme
        , k                                                             &
!           Loop variable
        , i_dummy
!           Dummy reading variable

!     Read the headers.
      READ(iu_spc, '(/, 27x, i5, /, 34x, i5, 27x, i5)')                 &
        i_drop, i_parametrization_drop, i_dummy
      IF (i_drop > npd_drop_type) THEN
        cmessage = '*** Error in subroutine ses_block_10_0_2. ' //      &
          'The indexing number of a droplet exceeds the maximum ' //    &
          'permitted value: increase npd_drop_type and recompile.'
        ierr=i_err_fatal
        RETURN
      END IF

      IF ( (i_parametrization_drop == ip_slingo_schrecker).OR.          &
           (i_parametrization_drop == ip_ackerman_stephens).OR.         &
           (i_parametrization_drop == ip_drop_pade_2) ) THEN
!       Data are parametrized.
        READ(iu_spc, '(42x, i5)') n_drop_phf_term(i_drop)
        READ(iu_spc, '(39x, 1pe12.5, 4x, 1pe12.5)')                     &
          drop_parm_min_dim(i_drop), drop_parm_max_dim(i_drop)

        n_parameter_water=i_dummy

        DO i=1, n_band
          READ(iu_spc, FMT='(/, (4(4x, 1pe12.5)))', IOSTAT=ios)         &
            (drop_parameter_list(k, i, i_drop), k=1, n_parameter_water)
!         For each band read the values of the parameters.
          IF (ios /= 0) THEN
            cmessage = '*** Error in subroutine ses_block_10_0_2. ' //  &
              'Data for droplets are not in the correct format.'
            ierr=i_err_fatal
            RETURN
          END IF
        END DO
      ELSE
!       Illegal parametrization scheme encountered.
        cmessage = '*** Error in subroutine ses_block_10_0_2. ' //      &
          'An unknown parametrization scheme has been specified.'
        ierr=i_err_fatal
        RETURN
      END IF

!     Record the presence of the drop type and the index
!     of the parametrization
      l_drop_type(i_drop)=.TRUE.
      i_drop_parametrization(i_drop)=i_parametrization_drop

      END SUBROUTINE ses_block_10_0_2

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 11.      *
!     *                AEROSOL SCATTERING DATA.                       *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_11_0_1

      USE rad_pcf, ONLY: ip_aerosol_param_dry

!     Local variables.
      INTEGER                                                           &
           i_species                                                    &
!             Index of species
         , idum
!             Dummy reading variable

      READ(iu_spc, '(/, 19X, I5, //)') i_species

!     Initialize humidities here for vectorization purpose
!     in grey_extinction
      DO i=1, npd_humidities
         humidities(i, i_species) = REAL(i-1)
      END DO

!     Read in the scattering parameters for each band
      DO i=1, n_band
         READ(iu_spc, FMT='(I5, 3(4X, 1PE16.9))', IOSTAT=ios) idum      &
            , aerosol_absorption(1, i_species, i)                       &
            , aerosol_scattering(1, i_species, i)                       &
            , aerosol_asymmetry(1, i_species, i)

!        Initialize zero component for vectorization
         aerosol_absorption(0,i_species, i)=0.0
         aerosol_scattering(0,i_species, i)=0.0
         aerosol_asymmetry(0,i_species, i)=0.0

         IF (ios /= 0) THEN
            cmessage = '*** Error in subroutine ses_block_11_0_1. ' //  &
             'Dry aerosol scattering data is not in the correct format.'
            ierr=i_err_fatal
            RETURN
         END IF
      END DO
      nhumidity(i_species)=1
      i_aerosol_parametrization(i_species)=ip_aerosol_param_dry

!     After successful reading the presence of this species is specified.
      l_aerosol_species(i_species)=.TRUE.

      END SUBROUTINE ses_block_11_0_1

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 11.      *
!     *                AEROSOL SCATTERING DATA.                       *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_11_0_2

      USE rad_pcf, ONLY: ip_aerosol_param_dry

!     Local variables.
      INTEGER                                                           &
           i_species                                                    &
!             Index of species
         , idum
!             Dummy reading variable

      READ(iu_spc, '(/, 19X, I5, ///)') i_species

!     Read in the scattering parameters for each band
      DO i=1, n_band
         READ(iu_spc, FMT='(I5, 3(4X, 1PE16.9))', IOSTAT=ios) idum      &
            , aerosol_absorption(1, i_species, i)                       &
            , aerosol_scattering(1, i_species, i)                       &
            , aerosol_asymmetry(1, i_species, i)

         IF (ios /= 0) THEN
            cmessage = '*** Error in subroutine ses_block_11_0_2. ' //  &
             'Dry aerosol scattering data is not in the correct format.'
            ierr=i_err_fatal
            RETURN
         END IF
      END DO
      nhumidity(i_species)=1
      i_aerosol_parametrization(i_species)=ip_aerosol_param_dry

!     After successful reading the presence of this species is specified.
      l_aerosol_species(i_species)=.TRUE.

      END SUBROUTINE ses_block_11_0_2

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 11.      *
!     *              MOIST AEROSOL SCATTERING DATA.                   *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_11_1_0

      USE rad_pcf, ONLY: ip_aerosol_param_moist

!     Local variables.
      INTEGER                                                           &
           k                                                            &
!             Loop variable
         , i_species
!             Index of component

      READ(iu_spc, '(/, 19X, I5 )') i_species
      READ(iu_spc, '(28X, I3)') nhumidity(i_species)

!     Read in the scattering parameters for each band
      DO i=1, n_band
         READ(iu_spc, '(//)')
         aerosol_absorption(0,i_species, i)=0.0
         aerosol_scattering(0,i_species, i)=0.0
         aerosol_asymmetry(0,i_species, i)=0.0

         DO k=1, nhumidity(i_species)
            READ(iu_spc, '(4(4X, 1PE12.5))', IOSTAT=ios)                &
                 humidities(k, i_species)                               &
               , aerosol_absorption(k, i_species, i)                    &
               , aerosol_scattering(k, i_species, i)                    &
               , aerosol_asymmetry(k, i_species, i)
            IF (ios /= 0) THEN
              cmessage = '*** Error in subroutine ses_block_11_1_0. '// &
                'Moist aerosol scattering data is ' //                  &
                'not in the correct format.'
              ierr=i_err_fatal
              RETURN
            END IF
         END DO
      END DO

      i_aerosol_parametrization(i_species)=ip_aerosol_param_moist

!     After successful reading the presence of this species is specified
      l_aerosol_species(i_species)=.TRUE.

      END SUBROUTINE ses_block_11_1_0

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 11.      *
!     *              MOIST AEROSOL SCATTERING DATA.                   *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_11_1_2

      USE rad_pcf, ONLY: ip_aerosol_param_moist

!     Local variables.
      INTEGER                                                           &
           k                                                            &
!             Loop variable
         , i_species
!             Index of component

      READ(iu_spc, '(/, 19X, I5 )') i_species
      READ(iu_spc, '(28X, I3,/)') nhumidity(i_species)

!     Read in the scattering parameters for each band
      DO i=1, n_band
         READ(iu_spc, '(//)')
         aerosol_absorption(0,i_species, i)=0.0
         aerosol_scattering(0,i_species, i)=0.0
         aerosol_asymmetry(0,i_species, i)=0.0
         DO k=1, nhumidity(i_species)
            READ(iu_spc, '(3x, 1pe16.9, 3(4x, 1pe16.9))', IOSTAT=ios)   &
                 humidities(k, i_species)                               &
               , aerosol_absorption(k, i_species, i)                    &
               , aerosol_scattering(k, i_species, i)                    &
               , aerosol_asymmetry(k, i_species, i)
            IF (ios /= 0) THEN
              cmessage = '*** Error in subroutine ses_block_11_1_1. '// &
                'Moist aerosol scattering data is ' //                  &
                'not in the correct format.'
               ierr=i_err_fatal
               RETURN
            END IF
         END DO
      END DO

      i_aerosol_parametrization(i_species)=ip_aerosol_param_moist

!     After successful reading the presence of this species is specified.
      l_aerosol_species(i_species)=.TRUE.

      END SUBROUTINE ses_block_11_1_2

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 12.      *
!     *            PARAMETRIZED ICE CRYSTAL SCATTERING DATA.          *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_12_0_0

      USE rad_pcf, ONLY:                                                &
        ip_ice_sun_fu, ip_ice_chou_vis, ip_ice_agg_de_sun

!     Local variables.
      INTEGER                                                           &
           i_ice                                                        &
!             Type of ice crystal
         , i_parametrization_ice                                        &
!             Dummy index of parameter scheme
         , k
!             Loop variable

!     Read the headers.
      READ(iu_spc, '(/, 31X, I5, /, 34X, I5, 27X, I5)')                 &
         i_ice, i_parametrization_ice, n_parameter_ice

      IF ( (i_parametrization_ice == ip_ice_sun_fu)                     &
         .OR.(i_parametrization_ice == ip_ice_chou_vis)                 &
         .OR.(i_parametrization_ice == ip_ice_agg_de_sun) ) THEN
!        Data is parametrized.
         DO i=1, n_band
            READ(iu_spc, *)
            READ(iu_spc,  *, IOSTAT=ios)                                &
               (ice_parameter_list(k, i, i_ice), k=1, n_parameter_ice)
!           For each band read the values of the parameters.
            IF (ios /= 0) THEN
               cmessage = '*** Error in subroutine ses_block_12_0_0. '//&
                 'Data for ice crystals are not in the correct format.'
               ierr=i_err_fatal
               RETURN
            END IF
         END DO
      ELSE
!        Illegal parametrization scheme encountered.
         cmessage = '*** Error in subroutine ses_block_12_0_0. ' //     &
           'An unknown parametrization scheme has been specified.'
         ierr=i_err_fatal
         RETURN
      END IF

!     Record the presence of the ice crystal type and the index
!     of the parametrization
      l_ice_type(i_ice)=.TRUE.
      i_ice_parametrization(i_ice)=i_parametrization_ice

      END SUBROUTINE ses_block_12_0_0

!     *****************************************************************
!     *                                                               *
!     *     Subroutine to read in the contents of block type 12.      *
!     *          Parametrized ice crystal scattering data.            *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE read_block_12_0_2

      USE rad_pcf, ONLY:                                                &
        ip_slingo_schrecker_ice, ip_sun_shine_vn2_vis,                  &
        ip_sun_shine_vn2_ir, ip_ice_adt, ip_ice_adt_10,                 &
        ip_ice_fu_solar, ip_ice_fu_ir,  ip_slingo_schr_ice_phf,         &
        ip_ice_fu_phf

!     Local variables.
      INTEGER                                                           &
          i_ice                                                         &
!           Type of ice crystal
        , i_parametrization_ice                                         &
!           Dummy index of parameter scheme
        , n_parameter_ice                                               &
!           Number of parameters
        , k                                                             &
!           Loop variable
        , i_dummy
!           Dummy reading variable


!     Read the headers.
      READ(iu_spc, '(/, 31x, i5, /, 34x, i5, 27x, i5)')                 &
        i_ice, i_parametrization_ice, i_dummy
      IF (i_ice > npd_ice_type) THEN
        cmessage = '*** Error in subroutine read_block_12_0_2. ' //     &
          'The indexing number of an ice crystal exceeds the maximum'// &
          ' permitted value: increase npd_ice_type and recompile.'
        ierr=i_err_fatal
        RETURN
      ENDIF

      IF ( (i_parametrization_ice == IP_slingo_schrecker_ice).OR.       &
           (i_parametrization_ice == IP_ice_adt).OR.                    &
           (i_parametrization_ice == IP_ice_adt_10).OR.                 &
           (i_parametrization_ice == IP_ice_fu_solar).OR.               &
           (i_parametrization_ice == IP_ice_fu_ir).OR.                  &
           (i_parametrization_ice == IP_slingo_schr_ice_phf).OR.        &
           (i_parametrization_ice == IP_ice_fu_phf).OR.                 &
           (i_parametrization_ice == IP_sun_shine_vn2_vis).OR.          &
           (i_parametrization_ice == IP_sun_shine_vn2_ir) ) THEN
!       Data are parametrized.
        n_parameter_ice=i_dummy

        READ(iu_spc, '(42x, i5)') n_ice_phf_term(i_ice)
        READ(iu_spc, '(39x, 1pe12.5, 4x, 1pe12.5)')                     &
          ice_parm_min_dim(i_ice), ice_parm_max_dim(i_ice)

        DO i=1, n_band
          READ(iu_spc, '()')
          READ(iu_spc, *, iostat=ios)                                   &
            (ice_parameter_list(k, i, i_ice), k=1, n_parameter_ice)
!         For each band read the values of the parameters.
          IF (ios /= 0) THEN
            cmessage = '*** Error in subroutine read_block_12_0_2. ' // &
              'Data for ice crystals are not in the correct format.'
            ierr=i_err_fatal
            RETURN
          ENDIF
        ENDDO
      ELSE
!       Illegal parametrization scheme encountered.
        cmessage = '*** Error in subroutine read_block_12_0_2. ' //     &
          'An unknown parametrization scheme has been specified.'
        ierr=i_err_fatal
        RETURN
      ENDIF

!     Record the presence of the ice crystal type and the index
!     of the parametrization
      l_ice_type(i_ice)=.true.
      i_ice_parametrization(i_ice)=i_parametrization_ice

      END SUBROUTINE read_block_12_0_2

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 14.      *
!     *        BANDS EXCLUDED FROM REGIONS OF THE SPECTRUM.           *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_14_0_0

!     Local variables.
      INTEGER                                                           &
           idum                                                         &
!             Dummy integer
         , j
!             Loop variable

!     Skip over the headers.
      READ(iu_spc, '(//)')

!     Read in the list of excluded bands for each band in turn.
      DO i=1, n_band
         READ(iu_spc, FMT='(I5, 7X, I5)', IOSTAT=ios)                   &
            idum, n_band_exclude(i)
         IF (ios /= 0) THEN
            cmessage = '*** Error in subroutine ses_block_14_0_0. ' //  &
              'The list of excluded bands is not correct.'
            ierr=i_err_fatal
            RETURN
         END IF
         IF (n_band_exclude(i) >  0) THEN
            READ(iu_spc, '(14X, 8(3X, I5))', IOSTAT=ios)                &
               (index_exclude(j, i), j=1, n_band_exclude(i) )
         END IF
         IF (ios /= 0) THEN
            cmessage = '*** Error in subroutine ses_block_14_0_0. ' //  &
              'The index of excluded bands is not correct.'
            ierr=i_err_fatal
            RETURN
         END IF
      END DO

      END SUBROUTINE ses_block_14_0_0

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 15.      *
!     *            mixed absorbing gas coefficients                   *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_15_0_0

!     Local variables.
      INTEGER                                                           &
           idum_band                                                    &
!             Dummy integer
         , j                                                            &
!             AOD type
         , k
!             Loop variable

      READ(iu_spc, '(/14X, I5)') n_aod_wavel
      READ(iu_spc, '(19x,I5, 15x, I5)') i, j
      i_aod_type(i) = j

      READ(iu_spc, *)
      DO k = 1 , n_aod_wavel
         READ(iu_spc, FMT='(I5,2(4X, 1PE12.5))', IOSTAT=ios)            &
            idum_band, aod_absorption(1, i,k), aod_scattering(1, i, k)
         IF (ios /= 0) THEN
            cmessage = '*** Error in subroutine ses_block_15_0_0. ' //  &
              'AOD data is not consistent with summary.'
            ierr=i_err_fatal
            RETURN
         END IF
      END DO

      END SUBROUTINE ses_block_15_0_0

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 15.      *
!     *     Moisture aerosols                                         *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_15_1_0

!     Local variables.
      INTEGER                                                           &
           idum_band                                                    &
!             Dummy integer
         , j                                                            &
!             AOD type
         , k
!             Loop variable

      READ(iu_spc, '(/14X, I5)') n_aod_wavel
      READ(iu_spc, '(19x,I5, 15x, I5)') i, j
      i_aod_type(i) = j

      DO k = 1 , n_aod_wavel

        READ(iu_spc, '(/)')
        DO j=1, 21
          READ(iu_spc, '(2(4X, 1PE12.5))', IOSTAT=ios)                  &
            aod_absorption(j, i,k), aod_scattering(j, i, k)
          IF (ios /= 0) THEN
            cmessage = '*** Error in subroutine ses_block_15_1_0. ' //  &
              'AOD data is not consistent with summary.'
            ierr=i_err_fatal
            RETURN
          END IF
        END DO

      END DO

      END SUBROUTINE ses_block_15_1_0

!     *****************************************************************
!     *                                                               *
!     *     SUBROUTINE TO READ IN THE CONTENTS OF BLOCK TYPE 15.      *
!     *            mixed absorbing gas coefficients                   *
!     *                                                               *
!     *****************************************************************
      SUBROUTINE ses_block_16_0_0

!     Local variables.
      INTEGER                                                           &
           idum_band                                                    &
!             Dummy integer
         , j, k, i_band, i_mix, it, ip, i_term
!             Loop variables

      CHARACTER char_dum*80

      READ(iu_spc, '(/)')

!     Read in the number of ESFT terms in each band.
      i = 0
      mix_gas_band(1:npd_band)=0
      DO i_band = 1 , n_band
         READ(iu_spc, FMT='(I5, 7X, I5)', IOSTAT=ios)                   &
            idum_band, n_mix_gas(i_band)
         IF (ios /= 0) THEN
            cmessage = '*** Error in subroutine ses_block_16_0_0. ' //  &
              'ESFT data is not consistent with summary.'
            ierr=i_err_fatal
            RETURN
         END IF

         IF ( n_mix_gas(i_band)  /=  0 ) THEN
            i=i+1

            mix_gas_band(i_band) = i

            READ(iu_spc,*)(index_mix_gas(j,i),j=1,n_mix_gas(i_band))
            READ(iu_spc1, '(/)')

            DO i_term = 1, i_band_esft(i_band)
               DO i_mix=1, 9
                  DO ip = 1, num_ref_p(i_band)
                     READ(iu_spc1, *, err=1000, IOSTAT=ios)             &
                     (k_mix_gas(ip,it,i_mix,i_term,i),it=1, npd_tmp)
                  END DO
               END DO
            END DO

1000        IF (ios /= 0) THEN
              cmessage = '*** Error in subroutine ses_block_16_0_0. '// &
                'The error is in block data file'
               PRINT*,'ip= ',ip,' i_mix= ',i_mix,' i_term= ',i_term,    &
                    ' i_band= ',i_band
               ierr=i_err_fatal
               RETURN
            END IF
         END IF
      END DO

      END SUBROUTINE ses_block_16_0_0

END SUBROUTINE ses_spectrum
