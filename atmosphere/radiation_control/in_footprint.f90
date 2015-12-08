! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Routine to determine whether a grid-point is within the
!          footprint seen by a satellite.

! Method:  This routine is subject to considerable alteration and
!          will need rewriting for each application. A satellite is
!          identified only by the number of the call.

! Code description:
!   Language: Fortran 90

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

FUNCTION in_footprint(i_call, l_sw, &
  l_geostationary, min_view_lon, max_view_lon, &
  min_view_lat, max_view_lat, &
  pt_lat, pt_lon, time_0) &
RESULT (l_in)

! Declarations

  USE conversions_mod, ONLY: pi
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE ereport_mod, ONLY : ereport
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: i_call
!   Number of call to radiation code
  LOGICAL, INTENT(IN) :: l_sw
!   Flag for Shortwave region
  REAL, INTENT(IN) :: pt_lon
!   Longitude of grid-point
  REAL, INTENT(IN) :: pt_lat
!   Latitude of grid-point
  REAL, INTENT(IN) :: time_0
!   Time at beginning of current timestep

  LOGICAL, INTENT(IN) :: l_geostationary
!   Flag for geostationary satellite
  REAL, INTENT(IN) :: min_view_lon
!   Minimum viewing longitude
  REAL, INTENT(IN) :: max_view_lon
!   Maximum viewing longitude
  REAL, INTENT(IN) :: min_view_lat
!   Minimum viewing latitude
  REAL, INTENT(IN) :: max_view_lat
!   Maximum viewing latitude

  LOGICAL :: l_in
!   Returned value of the function: true if the point is
!   in the footprint



! Local variables
  INTEGER :: error_code
!   Error code passed to termintaing routine in case of fatal errors
  REAL :: sat_lon
!   Longitude of viewing satellite
  REAL :: sat_lat
!   Latitude of viewing satellite

  REAL, PARAMETER :: period =  11000.0

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
!   Time at beginning of current timestep



! SW calls first:
  IF (lhook) CALL dr_hook('IN_FOOTPRINT',zhook_in,zhook_handle)
  SELECT CASE(l_sw)


    CASE(.TRUE.)

      SELECT CASE(i_call)

        CASE(1)
!         This should not be called by the main call.
          error_code = 1

          CALL ereport('in_footprint', error_code, &
            'Foot-printing in the main call is not permitted.')
        CASE(2)

          IF (l_geostationary) THEN
!           Simply test the preset window.
!           Latitudes first:
            l_in = (pt_lat >= min_view_lat) .AND. &
                   (pt_lat <= max_view_lat)
!           Longitudes, allowing for different conventions
            l_in = l_in .AND. &
                 ( ( (pt_lon >= min_view_lon) .AND. &
                     (pt_lon <= max_view_lon) ) .OR. &
                   ( (pt_lon - 2.0 * pi >= min_view_lon) .AND. &
                     (pt_lon - 2.0 * pi <= max_view_lon) ) )
          ELSE
!           This code is purely here for illustration.
            sat_lon = 2.0 * pi * REAL( time_0 / period - &
                        REAL(INT(time_0 / period) ) )
            sat_lat = pi * SIN (sat_lon)
            IF ( (ABS(pt_lon - sat_lon) < 0.2) .AND. &
                 (ABS(pt_lat - sat_lat) < 0.2) ) THEN
              l_in = .TRUE.
            ELSE
              l_in = .FALSE.
            END IF
          END IF

        CASE DEFAULT
!         Code is not available for this call.
          error_code=10

          CALL ereport('in_footprint', error_code, &
            'No footprinting is available for this call.')

      END SELECT


    CASE(.FALSE.)

      SELECT CASE(i_call)

        CASE(1)
!         This should not be called by the main call.
          error_code = 1

          CALL ereport('in_footprint', error_code, &
            'Foot-printing in the main call is not permitted.')
        CASE(2)

          IF (l_geostationary) THEN
!           Simply test the preset window.
!           Latitudes first:
            l_in = (pt_lat >= min_view_lat) .AND. &
                   (pt_lat <= max_view_lat)
!           Longitudes, allowing for different conventions
            l_in = l_in .AND. &
                 ( ( (pt_lon >= min_view_lon) .AND. &
                     (pt_lon <= max_view_lon) ) .OR. &
                   ( (pt_lon - 2.0 * pi >= min_view_lon) .AND. &
                     (pt_lon - 2.0 * pi <= max_view_lon) ) )
          ELSE
!           This code is purely here for illustration.
            sat_lon = 2.0 * pi * REAL( time_0 / period - &
                        REAL(INT(time_0 / period) ) )
            sat_lat = pi * SIN (sat_lon)
            IF ( (ABS(pt_lon - sat_lon) < 0.2) .AND. &
                 (ABS(pt_lat - sat_lat) < 0.2) ) THEN
              l_in = .TRUE.
            ELSE
              l_in = .FALSE.
            END IF
          END IF

        CASE DEFAULT
!         Code is not available for this call.
          error_code=10

          CALL ereport('in_footprint', error_code, &
            'No footprinting is available for this call.')

      END SELECT


    END SELECT
    IF (lhook) CALL dr_hook('IN_FOOTPRINT',zhook_out,zhook_handle)
    RETURN



END FUNCTION in_footprint
