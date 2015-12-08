! *****************************copyright*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************copyright*******************************

!+ Subroutine to set dimensions for the radiation code.

! Purpose:
!   To set the dimensions of arrays for use in the radiation code
!   depending on the options chosen.

! Method:
!   This routine covers a range of disparate requirements, working by
!   IF-tests.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of code:
!   Fortran 90

!- ---------------------------------------------------------------------
SUBROUTINE r2_set_rad_dim &
(row_length, rows, lat, lon, l_solar, sindec, seconds_since_midnight, &
 l_geostationary, sat_hgt, sat_lat, sat_lon, &
 i_cloud, i_angular_integration, i_sph_mode, &
 nd_field, n_points, model_levels, cloud_levels, &
 nd_cloud_type, nd_cloud_component, &
 l_extra_top, id_ct, n_rad_layers, nd_column, &
 nd_band, n_band, map_channel, n_channel, nd_channel, &
 nd_viewing_level, n_viewing_level, viewing_level, &
 nd_direction, n_view_direction, view_direction, &
 nd_brdf_basis_fnc, nd_brdf_trunc, &
 nd_profile, nd_flux_profile, nd_radiance_profile, &
 nd_field_flux_diag, nd_field_rad_diag, &
 l_ctile, nd_tile, nd_point_tile &
)



  USE rad_pcf
  USE earth_constants_mod, ONLY: earth_radius
  USE conversions_mod, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY : ereport
  IMPLICIT NONE


! Dummy arguments

  INTEGER, INTENT(IN) :: row_length
!   Length of EW rows
  INTEGER, INTENT(IN) :: rows
!   Number of NS rows

  REAL, INTENT(IN) :: lat(row_length, rows)
!   Latitudes of grid-points
  REAL, INTENT(IN) :: lon(row_length, rows)
!   Longitudes of grid-points
  LOGICAL, INTENT(IN) :: l_solar
!   Flag for solar region of the spectrum
  REAL, INTENT(IN) :: sindec
!   Sine of solar declination
  REAL, INTENT(IN) :: seconds_since_midnight
!   Seconds since midnight


  INTEGER, INTENT(IN) :: i_cloud
!   Option for cloud geomtery
  INTEGER, INTENT(IN) :: i_angular_integration
!   Type of angular integration
  INTEGER, INTENT(IN) :: i_sph_mode
!   Mode of operation of the spherical harmonic code
  LOGICAL, INTENT(IN) :: l_extra_top
!   Flag to include an extra layer in radiation to model the part
!   of the atmosphere above the top of the model

  LOGICAL, INTENT(IN) :: l_ctile
!   Flag for coastal tiling

  INTEGER, INTENT(IN) :: nd_field
!   Horizontal size of input fields
  INTEGER, INTENT(IN) :: nd_band
!   Size allocated for spectral bands
  INTEGER, INTENT(IN) :: n_band
!   Number of spectral bands used
  INTEGER, INTENT(IN) :: n_points
!   Number of grid points to operate on
  INTEGER, INTENT(IN) :: model_levels
!   Number of theta levels in the model
  INTEGER, INTENT(IN) :: cloud_levels
!   Number of potentially cloudy layers in the model

  LOGICAL, INTENT(IN) :: l_geostationary
!   Flag for a geostationary orbit
  REAL, INTENT(IN) :: sat_hgt
!   Height of the satellite's orbit
  REAL, INTENT(IN) :: sat_lon
!   Longitude of the (geostationary) satellite
  REAL, INTENT(IN) :: sat_lat
!   Laitude of the (geostationary) satellite

  INTEGER, INTENT(OUT) :: n_channel
!   Number of channels in output
  INTEGER, INTENT(OUT) :: map_channel(nd_band)
!   Mapping of spectral bands to output channels

  INTEGER, INTENT(OUT) :: n_rad_layers
!   Number of layers used in radiation
  INTEGER, INTENT(OUT) :: nd_cloud_type
!   Size to be allocated for types of cloud
  INTEGER, INTENT(OUT) :: nd_cloud_component
!   Size to be allocated for condensed components
  INTEGER, INTENT(OUT) :: nd_column
!   Size to be allocated for subcolumns in a grid-box
  INTEGER, INTENT(OUT) :: id_ct
!   Topmost layer allowed to contain cloud
  INTEGER, INTENT(OUT) :: nd_channel
!   Size to be allocated for satellite channels available from
!   one call to the radiation scheme
  INTEGER, INTENT(OUT) :: nd_viewing_level
!   Size to be allocated for viewing levels, where fluxes or radiances
!   will be calculated
  INTEGER, INTENT(OUT) :: nd_direction
!   Size to be allocated for viewing directions
  INTEGER, INTENT(OUT) :: nd_profile
!   Size to be allocated for points where fields are required
!   for both flux and radiance calculations
  INTEGER, INTENT(OUT) :: nd_radiance_profile
!   Size to be allocated for points, where radiances
!   will be calculated
  INTEGER, INTENT(OUT) :: nd_flux_profile
!   Size to be allocated for points, where fluxes
!   will be calculated
  INTEGER, INTENT(OUT) :: nd_field_flux_diag
!   Size to be allocated for flux diagnostics
  INTEGER, INTENT(OUT) :: nd_field_rad_diag
!   Size to be allocated for radiance diagnostics

  INTEGER, INTENT(OUT) :: nd_brdf_basis_fnc
!   Size to be allocated for basis functions for BRDFs
  INTEGER, INTENT(OUT) :: nd_brdf_trunc
!   Size to be allocated for order of truncation fo be applied to
!   BRDFs

  INTEGER, INTENT(OUT) :: nd_tile
!   Size to be allocated for number of surface tiles within the
!   radiation scheme
  INTEGER, INTENT(OUT) :: nd_point_tile
!   Size to be allocated for number of points where surface tiling will
!   be applied

  INTEGER, INTENT(OUT) :: n_viewing_level
!   Actual number of viewing levels
  REAL, INTENT(OUT) :: viewing_level(model_levels+1)
!   Viewing levels (this hard-wired dimension is a temporary measure)

  INTEGER, INTENT(OUT) :: n_view_direction
!   Actual number of viewing directions
  REAL, INTENT(OUT) :: view_direction(nd_field, 1, 2)
!   Viewing directions (this hard-wired dimension is a temporary measure)
!   The second dimension allows just one viewing direction



! Local variables
  INTEGER :: i
!   Loop variable
  INTEGER :: j
!   Loop variable
  INTEGER :: l
!   Loop variable

  REAL :: r_sat(3)
!   Cartesian coordinates of the satellite
  REAL :: r_obs(3)
!   Cartesian coordinates of observed point
  REAL :: r_diff(3)
!   Cartesian coordinates of separation of the points
  REAL :: r_sun(3)
!    Cartesian coordinates of the sun
  REAL :: r_hor(3)
!    Projection of r_diff on local horizontal plane
  REAL :: r_hor_s(3)
!    Projection of r_sun on local horizontal plane
  REAL :: mag_rdiff
!   Magnitude of the separation
  REAL :: mag_robs
!   Magnitude of r_obs
  REAL :: mag_rhor
!    Magnitude or r_hor
  REAL :: mag_rhor_s
!    Magnitude of r_hor_s
  REAL :: solar_zen
!    Solar zenith angle
  REAL :: h_angle
!    Hour Angle




  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('R2_SET_RAD_DIM',zhook_in,zhook_handle)



! Set the actual size of arrays in the radiation code:
! for some architectures (such as that of Cray vector
! machines) on odd size is preferred to avoid memory
! bank conflicts.
  nd_profile=2*(n_points/2)+1

! Set the number of layers seen in the radiation code.
! This may optionally be 1 greater than the number used
! in the rest of the model to avoid spurious effects
! resulting from the upper boundary (principally in
! stratospheric configurations).
  IF (l_extra_top) THEN
    n_rad_layers=model_levels+1
  ELSE
    n_rad_layers=model_levels
  END IF

! Standard sizes for cloud arrays.
  nd_cloud_type      = 4
  nd_cloud_component = 4
  SELECT CASE(i_cloud)
    CASE(ip_cloud_column_max)
      nd_column = 3 * cloud_levels + 2
    CASE DEFAULT
      nd_column = 1
  END SELECT


  IF (   (i_angular_integration == ip_two_stream) .OR. &
       ( (i_angular_integration == ip_spherical_harmonic) .AND. &
         (i_sph_mode == ip_sph_mode_flux) ) ) THEN

    SELECT CASE(i_angular_integration)
      CASE(ip_two_stream)
        nd_viewing_level     = 1
        n_viewing_level      = 0
        nd_radiance_profile  = 1
      CASE(ip_spherical_harmonic)
        nd_viewing_level     = model_levels+1
        n_viewing_level      = model_levels+1
        nd_radiance_profile  = nd_profile
    END SELECT
    nd_direction       = 1
    nd_brdf_basis_fnc  = 2
    nd_brdf_trunc      = 1
    nd_field_flux_diag = nd_field
    nd_flux_profile    = nd_profile
    nd_field_rad_diag  = 1

    DO i=1, model_levels + 1
      viewing_level(i) = REAL(i-1)
    END DO

!   Dummy initailization.
    n_view_direction = 1

!   Set the mapping of bands in the spectral file onto channels for
!   the diagnostics. When calculating fluxes all bands are normally
!   mapped into a single channel.
    nd_channel            = 1
    n_channel             = 1
    map_channel(1:n_band) = 1

  ELSE

!   Set up space for only one viewing direction initially.

    nd_direction     = 1
    n_view_direction = 1

!   Set the mapping of bands in the spectral file onto channels for
!   the diagnostics. We assume for now that each band in the spectral
!   file corresponds to a particular channel.
    nd_channel = n_band
    n_channel  = n_band
    DO i=1, n_band
      map_channel(i)=i
    END DO

    SELECT CASE(l_geostationary)

      CASE(.TRUE.)

! Define the Cartesian position of the satellite.

        r_sat(1) = ( earth_radius + sat_hgt ) * &
                     COS(sat_lat) * COS(sat_lon)
        r_sat(2) = ( earth_radius + sat_hgt ) * &
                     COS(sat_lat) * SIN(sat_lon)
        r_sat(3) = ( earth_radius + sat_hgt ) * &
                     SIN(sat_lat)


        DO j=1, rows
          DO i=1, row_length

            l = i+(j-1)*row_length


! Calculate zenith angle


! Define point in satelite footprint the satelite is lookin at:

            r_obs(1) = earth_radius * &
                         COS(lat(i, j)) * COS(lon(i, j))
            r_obs(2) = earth_radius * &
                         COS(lat(i, j)) * SIN(lon(i, j))
            r_obs(3) = earth_radius * &
                         SIN(lat(i, j))

            r_diff(1) = r_sat(1) - r_obs(1)
            r_diff(2) = r_sat(2) - r_obs(2)
            r_diff(3) = r_sat(3) - r_obs(3)

            mag_rdiff = SQRT( r_diff(1) * r_diff(1) + &
                              r_diff(2) * r_diff(2) + &
                              r_diff(3) * r_diff(3) )

            mag_robs  = SQRT( r_obs(1) * r_obs(1) + &
                              r_obs(2) * r_obs(2) + &
                              r_obs(3) * r_obs(3) )


            view_direction(l, 1, 1) = ( r_obs(1) * r_diff(1) + &
                                        r_obs(2) * r_diff(2) + &
                                        r_obs(3) * r_diff(3) ) / &
                                      ( earth_radius * mag_rdiff )


! Calculate Azimuthal angle

            IF (l_solar) THEN

! Define the normalised vector pointing to the sun

              h_angle = 2.0*pi*(0.5-seconds_since_midnight/86400.0)

              r_sun(1)=COS(sindec)*COS(h_angle)
              r_sun(2)=COS(sindec)*SIN(h_angle)
              r_sun(3)=SIN(sindec)

! Projection of r_diff on local horizontal plane

              r_hor(1)=r_diff(1)-view_direction(l,1,1) &
                 *r_obs(1)*(mag_rdiff/mag_robs)
              r_hor(2)=r_diff(2)-view_direction(l,1,1) &
                 *r_obs(2)*(mag_rdiff/mag_robs)
              r_hor(3)=r_diff(3)-view_direction(l,1,1) &
                 *r_obs(3)*(mag_rdiff/mag_robs)


! Projection of unity vector r_sun on horizontal plane


              solar_zen=(r_obs(1)*r_sun(1)+r_obs(2)*r_sun(2)  &
                       +r_obs(3)*r_sun(3))/mag_robs

! 1.0 refers to magnitude of unit vector pointing to the sun.
              r_hor_s(1)=r_sun(1)-solar_zen*r_obs(1)*(1.0/mag_robs)
              r_hor_s(2)=r_sun(2)-solar_zen*r_obs(2)*(1.0/mag_robs)
              r_hor_s(3)=r_sun(3)-solar_zen*r_obs(3)*(1.0/mag_robs)


! Calculate the azimuth angle and store in view_direction(l,1,2)


              mag_rhor  = SQRT( r_hor(1) * r_hor(1) + &
                                r_hor(2) * r_hor(2) + &
                                r_hor(3) * r_hor(3) )

              mag_rhor_s  = SQRT( r_hor_s(1) * r_hor_s(1) + &
                                  r_hor_s(2) * r_hor_s(2) + &
                                  r_hor_s(3) * r_hor_s(3) )

              view_direction(l,1,2)= - (r_hor(1)*r_hor_s(1) &
                                      + r_hor(2)*r_hor_s(2) &
                                      + r_hor(3)*r_hor_s(3)) &
                                      /mag_rhor/mag_rhor_s

              view_direction(l,1,2)=ACOS(view_direction(l,1,2))

            ELSE

!             For infra-red simulations we set the azimuth to 0.0.
              view_direction(l, 1, 2) = 0.0

            END IF

          END DO
        END DO

      CASE(.FALSE.)


        CALL ereport("rddim", 21, &
          "Only geostationary satellites are available so far")

    END SELECT

!   There should be only one viewing level: that at the top of the
!   atmosphere.
    nd_viewing_level  = 1
    n_viewing_level  = 1
    viewing_level(1) = 0.0

    nd_field_flux_diag  = 1
    nd_flux_profile     = 1
    nd_radiance_profile = nd_profile
    nd_brdf_basis_fnc   = 2
    nd_brdf_trunc       = 0
    nd_field_flux_diag  = 1
    nd_field_rad_diag   = nd_field

  END IF

! Remaining common sizes.
  nd_point_tile = n_points
  IF (l_ctile) THEN
    nd_tile     = 3
  ELSE
    nd_tile     = 2
  END IF
  id_ct         = n_rad_layers+1-cloud_levels

  IF (lhook) CALL dr_hook('R2_SET_RAD_DIM',zhook_out,zhook_handle)

END SUBROUTINE r2_set_rad_dim
