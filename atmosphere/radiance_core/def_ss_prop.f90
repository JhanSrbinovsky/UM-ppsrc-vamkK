! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to define the elements of single scattering properties
!
! Description:
!   This module defines the components of the structure of single
!   scattering propeties,as used throughout the radiation code: 
!   forward scattering properties and raw extinctions are 
!   included for convenience. 
!
!   The model for storage involves a split between the clear-sky
!   properties, dimensioned over (profiles, cloud-free layers)
!   and those in potentially cloudy regions, dimensioned over
!   (profiles, top-most cloudy layer: bottom layer, 0: cloud types)
!   A zero in the third entry refers to clear-sky properties.
!
!   A reordering of this structure to store the clear-sky properties
!   over (profiles, layers) and the cloudy properties over
!   (profiles, top-most cloudy layer: bottom layer, 1:cloud types)
!   would be possible. Which is more convenient depends on the
!   algorithms for cloud overlap.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiance Core
!
!- ---------------------------------------------------------------------
MODULE def_ss_prop

  USE realtype_rd

  IMPLICIT NONE


  TYPE STR_ss_prop

    REAL (RealK), POINTER :: k_grey_tot_clr(:, :)              => NULL()
!     Grey extinction in clear-sky region above clouds
    REAL (RealK), POINTER :: k_ext_scat_clr(:, :)              => NULL()
!     Scattering in clear-sky region above clouds
    REAL (RealK), POINTER :: phase_fnc_clr(:, :, :)            => NULL()
!     Phase function in clear-sky region above clouds
!     (Held as moments)
    REAL (RealK), POINTER :: forward_scatter_clr(:, :)         => NULL()
!     Forward scattering in clear-sky region above clouds
    REAL (RealK), POINTER :: phase_fnc_solar_clr(:, :, :)      => NULL()
!     Solar phase function in clear-sky region above clouds
!     (Held as the actual phase function in the viewing direction)
    REAL (RealK), POINTER :: forward_solar_clr(:, :)           => NULL()
!     Solar forward scattering in clear-sky region above clouds

    REAL (RealK), POINTER :: tau_clr(:, :)                     => NULL()
!     Optical depth in clear-sky region above clouds
    REAL (RealK), POINTER :: omega_clr(:, :)                   => NULL()
!     Albedo of single scattering in clear-sky region above clouds

    REAL (RealK), POINTER :: k_grey_tot(:, :, :)               => NULL()
!     Grey extinction in potentially cloudy regions
    REAL (RealK), POINTER :: k_ext_scat(:, :, :)               => NULL()
!     Scattering in potentially cloudy regions
    REAL (RealK), POINTER :: phase_fnc(:, :, :, :)             => NULL()
!     Phase function in potentially cloudy regions
!     (Held as moments)
    REAL (RealK), POINTER :: forward_scatter(:, :, :)          => NULL()
!     Forward scattering in potentially cloudy regions
    REAL (RealK), POINTER :: phase_fnc_solar(:, :, :, :)       => NULL()
!     Solar phase function in potentially cloudy regions
!     (Held as the actual phase function in the viewing direction)
    REAL (RealK), POINTER :: forward_solar(:, :, :)            => NULL()
!     Solar forward scattering in potentially cloudy regions

    REAL (RealK), POINTER :: tau(:, :, :)                      => NULL()
!     Optical depth in potentially cloudy regions
    REAL (RealK), POINTER :: omega(:, :, :)                    => NULL()
!     Albedo of single scattering in potentially cloudy regions

    REAL (RealK), POINTER :: k_ext_tot_cloud_comp(:, :, :)     => NULL()
!       Total extinction of cloudy component
    REAL (RealK), POINTER :: k_ext_scat_cloud_comp(:, :, :)    => NULL()
!       Scattering extinction of cloudy component
    REAL (RealK), POINTER :: phase_fnc_cloud_comp(:, :, :, :)  => NULL()
!       Phase function of cloudy components
    REAL (RealK), POINTER :: forward_scatter_cloud_comp(:, :, :) => NULL()
!       Forward scattering of cloudy component
    REAL (RealK), POINTER :: phase_fnc_solar_cloud_comp(:, :, :, :) => NULL()
!       Solar phase function of cloudy components
    REAL (RealK), POINTER :: forward_solar_cloud_comp(:, :, :) => NULL()
!       Solar forward scattering in the cloudy component

!   For McICA:
    REAL (RealK), POINTER :: phase_fnc_no_cloud(:, :, :)       => NULL()
!       Cloudy phase_fnc before cloud contribution is added
    REAL (RealK), POINTER :: forward_scatter_no_cloud(:, :)    => NULL()
!       Cloudy forward_scatter before cloud contribution is added

  END TYPE STR_ss_prop


END MODULE def_ss_prop
