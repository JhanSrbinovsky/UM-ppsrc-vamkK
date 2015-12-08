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

MODULE control_struc

  USE filenamelength_mod, ONLY:                                       &
      filenamelength

  IMPLICIT NONE

  TYPE control_option

!   Name of spectral file
    CHARACTER  (LEN=filenamelength) :: spectral_file

!   Range of spectral bands
    INTEGER :: first_band
!     First band to use in the calculation
    INTEGER :: last_band
!     Last band to use in the calculation

!   Miscallaneous options
    INTEGER :: isolir
!     Spectral region
    INTEGER :: i_solar_src
!     Index of solar source function
    INTEGER :: i_scatter_method
!     Method of treating scattering
    LOGICAL :: l_ir_source_quad
!     Flag to use a quadratic source function in the IR
    LOGICAL :: l_extra_top
!     Flag to insert an extra layer into radiation above the
!     top of the model (this is sometimes desirable to ensure
!     proper radiative heating in the top resolved layer).
    LOGICAL :: l_rad_deg
!     Flag to apply spatial degradation in radiation: in this case
!     radiation quantities are interpolated to all grid-points,
!     whereas subsampling refers to selecting a portion of the area
!     on the PE and returning radiative quentatities only at those
!     points
    LOGICAL :: l_subsample
!     Flag to apply spatial subsampling (for satellite footprints)
!     in radiation.

!   Properties of clouds:
    INTEGER :: i_cloud
!     Cloud scheme
    INTEGER :: i_cloud_representation
!     Representation of clouds
    INTEGER :: i_st_water
!     Type of water droplet in stratiform clouds
    INTEGER :: i_cnv_water
!     Type of water droplet in convective clouds
    INTEGER :: i_st_ice
!     Type of ice crystal in stratiform clouds
    INTEGER :: i_cnv_ice
!     Type of ice crystal in convective clouds
    LOGICAL :: l_local_cnv_partition
!     Flag to partition convective clouds between water and
!     ice using the local temperature
    LOGICAL :: l_global_cloud_top
!     Flag to use a global value for the topmost cloudy layer
!     (This is used to obtained bit-reproducible results
!     across different configurations of PEs on MPP systems)
    INTEGER :: i_fsd
!     method for obtaining values of cloud water content variability
    INTEGER :: i_inhom
!     method for treating of cloud water content variability
    INTEGER :: i_overlap
!     method for treating of cloud vertcial overlap

!   Physical processes
    LOGICAL :: l_microphysics
!     Flag for microphysics
    LOGICAL :: l_gas
!     Flag for gaseous absorption
    LOGICAL :: l_rayleigh
!     Flag for Rayleigh scattering
    LOGICAL :: l_continuum
!     Flag for the continuum
    LOGICAL :: l_cloud
!     Flag for clouds
    LOGICAL :: l_drop
!     Flag for droplets
    LOGICAL :: l_ice
!     Flag for ice crystals
    LOGICAL :: l_aerosol
!     Flag for aerosols
    LOGICAL :: l_aerosol_ccn
!     Flag for aerosols as CCN
    LOGICAL :: l_solar_tail_flux
!     Flag for adding solar tail flux to LW ragion

!   Gaseous absorption:
    INTEGER :: i_gas_overlap
!     Treatment of gaseous overlaps
    LOGICAL :: l_o2
!     Flag for absorption by oxygen
    LOGICAL :: l_n2o
!     Flag for absorption by nitrous oxide
    LOGICAL :: l_ch4
!     Flag for absorption by methane
    LOGICAL :: l_cfc11
!     Flag for absorption by CFC11
    LOGICAL :: l_cfc12
!     Flag for absorption by CFC12
    LOGICAL :: l_cfc113
!     Flag for absorption by CFC113
    LOGICAL :: l_cfc114
!     Flag for absorption by CFC114
    LOGICAL :: l_hcfc22
!     Flag for absorption by HCFC22
    LOGICAL :: l_hfc125
!     Flag for absorption by HFC125
    LOGICAL :: l_hfc134a
!     Flag for absorption by HFC134A

!   Angular integration:
    INTEGER :: i_angular_integration
!     Method of angular integration
    INTEGER :: i_2stream
!     Two-stream scheme
    INTEGER :: i_solver
!     Two-stream solver
    INTEGER :: n_order_gauss
!     Order of Gaussian quadrature
    INTEGER :: i_truncation
!     Type of truncation for spherical harmonics
    INTEGER :: i_sph_algorithm
!     Algorithm used for spherical harmonic calculations
    INTEGER :: n_order_phase_solar
!     Order of truncation of the solar phase function
    INTEGER :: ls_global_trunc
!     Global order of truncation
    INTEGER :: ms_min
!     Minimum azimuthal order
    INTEGER :: ms_max
!     Maximum azimuthal order
    INTEGER :: ls_brdf_trunc
!     Order of truncation of BRDFs
    REAL :: accuracy_adaptive
!     Accuracy for adaptive truncation
    LOGICAL :: l_rescale
!     Flag for rescaling

    LOGICAL :: l_henyey_greenstein_pf
!     Flag to use Henyey-Greenstein phase functions
    INTEGER :: i_sph_mode
!     Mode of operation of spherical harmonic code
    LOGICAL :: l_euler_trnf
!     Flag to apply Euler's transformation to alternating series



!   Satellite Data:

!   Current provisions are largely intended for geostationary
!   satellites. Note that in accordance with the UM's conventions
!   these must be given in SI units; i.e. angels will be in radians
!   etc.

    LOGICAL :: l_geostationary
!     Flag to signal that a geostationary satellite is assumed.
    CHARACTER  (LEN=80)  :: sat_desc
!     String for description of satellite
    REAL :: sat_hgt
!     Height of the orbit above the Earth's surface
    REAL :: sat_lon
!     Longitude of the (geostationary) satellite
    REAL :: sat_lat
!     Latitude of the (geostationary) satellite (in practice, for
!     a geostationary satellite this must be 0.0)

!   Viewing domain:
    REAL :: max_view_lon
!     Maximum longitude of viewing domain
    REAL :: min_view_lon
!     Minimum longitude of viewing domain
    REAL :: max_view_lat
!     Maximum latitude of viewing domain
    REAL :: min_view_lat
!     Minimum latitude of viewing domain

    LOGICAL :: l_layer
!     Flag for properties in layers
    LOGICAL :: l_cloud_layer
!     Flag for cloudy properties in layers
    LOGICAL :: l_2_stream_correct
!     Flag for corrections to 2-stream scheme

  END TYPE control_option

END MODULE control_struc
