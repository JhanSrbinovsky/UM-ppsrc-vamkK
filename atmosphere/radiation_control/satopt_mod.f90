! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     Module defining satellite options for the radiation code.
!     The variables here parallel those in the controlling structure
!     because elements of a structure cannot be set in a namelist.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE sat_opt_mod

IMPLICIT NONE

      LOGICAL :: l_subsample
!       Flag to signal that subsmapling for the footprint is required
      LOGICAL :: l_geostationary
!       Flag to signal that a geostationary satellite is assumed.
      CHARACTER  (LEN=80)  :: sat_desc
!       String for description of satellite
      REAL :: sat_hgt
!       Height of the orbit above the Earth's surface
      REAL :: sat_lon
!       Longitude of the (geostationary) satellite
      REAL :: sat_lat
!       Latitude of the (geostationary) satellite (in practice, for
!       a geostationary staellite this must be 0.0)
!
!     Viewing domain:
      REAL :: max_view_lon
!       Maximum longitude of viewing domain
      REAL :: min_view_lon
!       Minimum longitude of viewing domain
      REAL :: max_view_lat
!       Maximum latitude of viewing domain
      REAL :: min_view_lat
!       Minimum latitude of viewing domain
!
!     ------------------------------------------------------------------
END MODULE sat_opt_mod
