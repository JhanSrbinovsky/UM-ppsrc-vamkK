! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!   This module defines the elements of the structure defining
!   algorithmic control of the radiation code.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90

MODULE dec_spec

IMPLICIT NONE

  TYPE spectrum


!   Sizes of arrays:

    INTEGER :: npd_band
!     Number of spectral bands
    INTEGER :: npd_exclude
!     Numer of excluded bands
    INTEGER :: npd_esft_term
!     Number of esft terms
    INTEGER :: npd_pre
!     Number of reference pressure for esft
    INTEGER :: npd_tmp
!     Number of reference temperature for esft
    INTEGER :: npd_mix
!     Number of interpolation parameters for mixture species
    INTEGER :: npd_band_mix_gas
!     Number of bands where mixed species exist
    INTEGER :: npd_type
!     Number of data types
    INTEGER :: npd_species
!     Number of gaseous species
    INTEGER :: npd_scale_fnc
!     Number of scaling functions
    INTEGER :: npd_scale_variable
!     Number of scaling variables
    INTEGER :: npd_surface
!     Number of surface types in LW spectrum
    INTEGER :: npd_albedo_parm
!     Number of albedo parameters in LW spectrum
    INTEGER :: npd_continuum
!     Number of continua
    INTEGER :: npd_drop_type
!     Number of drop types
    INTEGER :: npd_ice_type
!     Number of ice crystal types
    INTEGER :: npd_aerosol_species
!     Number of aerosol species in spectral information
    INTEGER :: npd_aerosol_mixratio
!     Number of aerosol species in mixing ratio information
    INTEGER :: npd_thermal_coeff
!     Number of thermal coefficients
    INTEGER :: npd_cloud_parameter
!     Number of cloud parameters
    INTEGER :: npd_humidities
!     Number of humidities
    INTEGER :: npd_aod_wavel
!     Number of wavelengths for aod in LW spectrum
    INTEGER :: npd_phase_term
!     Number of terms in the phase function

! General fields:

    LOGICAL, POINTER ::    l_present(:)
!     Flag for types of data present

! Properties of the spectral bands:

    INTEGER :: n_band
!     Number of spectral bands

    REAL, POINTER :: wave_length_short(:)
!     Shorter wavelength limits
    REAL, POINTER :: wave_length_long(:)
!     Longer wavelength limits



! Exclusion of specific bands from parts of the spectrum:

    INTEGER, POINTER :: n_band_exclude(:)
!     Number of excluded bands within each spectral band
    INTEGER, POINTER :: index_exclude(:, :)
!     Indices of excluded bands



! Fields for the solar flux:

    REAL, POINTER :: solar_flux_band(:)
!     Fraction of the incident solar flux in each band
    REAL, POINTER :: solar_flux_band_ses(:, :)
!     Fraction of the incident solar flux in each g interval in each band
    REAL, POINTER :: weight_blue(:)
!     Fraction of the surface flux designated as "blue" in each band


! Fields for rayleigh scattering:

    REAL, POINTER :: rayleigh_coefficient(:)
!     Rayleigh coefficients



! Fields for gaseous absorption:

    INTEGER :: n_absorb
!     Number of absorbers
    INTEGER, POINTER :: n_band_absorb(:)
!     Number of absorbers in each band
    INTEGER, POINTER :: index_absorb(:, :)
!     List of absorbers in each band
    INTEGER, POINTER :: type_absorb(:)
!     Types of each gas in the spectral file
    INTEGER, POINTER :: n_mix_gas(:)
!     Number of mixed gases in a band
    INTEGER, POINTER :: index_mix_gas(:, :)
!     Index of mixed absorbers in each band
    INTEGER, POINTER :: num_mix(:)
!     Number of binary parameter for interpolation of absorption
!     coefficient for mixture of two species
    INTEGER, POINTER :: mix_gas_band(:)
!     Sequence band number (not real band number) of mixed species
    INTEGER, POINTER :: num_ref_p(:)
!     Number of refrence pressures
    INTEGER, POINTER :: i_band_esft(:, :)
    INTEGER, POINTER :: i_band_esft_ses(:)
!     Number of esft terms in band for each gas
    INTEGER, POINTER :: i_scale_esft(:, :)
!     Type of esft scaling
    INTEGER, POINTER :: i_scale_fnc(:, :)
!     Type of scaling function

    REAL, POINTER :: k_esft(:, :, :)
!     ESFT exponents
    REAL, POINTER :: w_esft(:, :, :)
!     ESFT weights
    REAL, POINTER :: scale_vector(:, :, :, :)
!     Scaling parameters for each absorber and term
    REAL, POINTER :: p_reference(:, :)
!     Reference pressure for scaling function
    REAL, POINTER :: t_reference(:, :)
!     Reference temperature for scaling function



! Representation of the Planckian:

    INTEGER :: n_deg_fit
!     Degree of thermal polynomial

    REAL, POINTER :: thermal_coefficient(:, :)
!     Coefficients in polynomial fit to source function
    REAL :: t_ref_planck
!     Planckian reference temperature


! Surface properties

    INTEGER, POINTER :: i_spec_surface(:)
!     Method of specifying properties of surface
    INTEGER, POINTER :: n_dir_albedo_fit(:)
!     Number of parameters fitting the direct albedo

    LOGICAL, POINTER :: l_surface(:)
!     Surface types included

    REAL, POINTER :: surface_albedo(:,:)
!     Surface albedos
    REAL, POINTER :: direct_albedo_parm(:,:,:)
!     Coefficients for fitting the direct albedo
    REAL, POINTER :: emissivity_ground(:,:)
!     Surface emissivities


! Fields for continua:

    INTEGER, POINTER :: n_band_continuum(:)
!     Number of continua in each band
    INTEGER, POINTER :: index_continuum(:, :)
!     List of continua continuua in each band
    INTEGER :: index_water
!     Index of water vapour
    INTEGER, POINTER :: i_scale_fnc_cont(:, :)
!     Type of scaling function for continuum

    REAL, POINTER :: k_continuum(:, :)
!     Grey extinction coefficients for continuum
    REAL, POINTER :: scale_continuum(:, :, :)
!     Scaling parameters for continuum
    REAL, POINTER :: p_ref_continuum(:, :)
!     Reference pressure for scaling of continuum
    REAL, POINTER :: t_ref_continuum(:, :)
!     Reference temperature for scaling of continuum



! Fields for water droplets:

    INTEGER, POINTER :: i_drop_parametrization(:)
!     Parametrization type of droplets

    LOGICAL, POINTER :: l_drop_type(:)
!     Types of droplet present

    INTEGER, POINTER :: n_drop_phf_term(:)
!     Number of moments in the phase function for droplets
    REAL, POINTER :: drop_parameter_list(:, :, :)
!     Parameters used to fit optical properties of clouds
    REAL, POINTER :: drop_parm_min_dim(:)
!     Minimum dimension permissible in the parametrization
    REAL, POINTER :: drop_parm_max_dim(:)
!     Maximum dimension permissible in the parametrization



! Fields for aerosols:

    INTEGER :: n_aerosol
!     Number of species of aerosol in spectral information
    INTEGER, POINTER :: type_aerosol(:)
!     Types of aerosols
    INTEGER, POINTER :: i_aerosol_parametrization(:)
!     Parametrization of aerosols
    INTEGER, POINTER :: n_aerosol_phf_term(:)
!     Number of terms in the phase function
    INTEGER, POINTER :: nhumidity(:)
!     Numbers of humidities
    INTEGER :: n_aerosol_mr
!     Number of species of aerosol mixing ratios

    LOGICAL, POINTER :: l_aerosol_species(:)
!     Aerosol species included

    REAL, POINTER :: aerosol_absorption(:, :, :)
!     Absorption by aerosols
    REAL, POINTER :: aerosol_scattering(:, :, :)
!     Scattering by aerosols
    REAL, POINTER :: aerosol_phase_fnc(:, :, :, :)
!     Phase function of aerosols
    REAL, POINTER :: humidities(:, :)
!     Humidities for components



! Fields for ice crystals:

    INTEGER, POINTER :: i_ice_parametrization(:)
!     Types of parametrization of ice crystals

    LOGICAL, POINTER :: l_ice_type(:)
!     Types of ice crystal present

    INTEGER, POINTER :: n_ice_phf_term(:)
!     Number of moments in the phase function for ice crystals
    REAL, POINTER :: ice_parameter_list(:, :, :)
!     Parameters used to fit single scattering of ice crystals
    REAL, POINTER :: ice_parm_min_dim(:)
!     Minimum dimension permissible in the parametrization
    REAL, POINTER :: ice_parm_max_dim(:)
!     Maximum dimension permissible in the parametrization



! Fields for Doppler broadening:

    LOGICAL, POINTER :: l_doppler_present(:)
!     Flag for Doppler broadening for each species

    REAL, POINTER :: doppler_correction(:)
!     Doppler correction terms



! Fields for aerosol optical depth:

    INTEGER :: n_aod_wavel
!     number of wavelengths
    INTEGER, POINTER :: i_aod_type(:)
!     relationship between aerosol component and type
    REAL, POINTER :: aod_wavel(:)
!     wavelengths for the aod
    REAL, POINTER :: aod_absorption(:, :, :)
!     monochromatic specific absorption coefficient
    REAL, POINTER :: aod_scattering(:, :, :)
!     monochromatic specific scattering coefficient

! Large arrays from ses2 defined here for efficient debugging

    REAL, POINTER :: k_esft_ses(:, :, :, :, :)
    REAL, POINTER :: w_esft_ses(:, :)
    REAL, POINTER :: k_mix_gas(:, :, :, :, :)
!     Absorption coefficients for mixture species
    REAL, POINTER :: f_mix(:)
!     Mixing ratio of mixed absorber amount
    REAL, POINTER :: k_continuum_ses(:, :, :, :)
    REAL, POINTER :: k_h2oc(:, :, :, :)
!     Absorption coefficient for water vapour continuum
    REAL, POINTER :: plk_frac(:, :)
!     Planckian fraction function
    REAL, POINTER :: planckb_ref(:, :)
!     Planckian function at 161 reference temperatures
  END TYPE spectrum



END MODULE dec_spec
