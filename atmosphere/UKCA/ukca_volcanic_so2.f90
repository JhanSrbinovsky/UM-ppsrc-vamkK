! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject 
! to the terms and conditions set out therein.
! [Met Office Ref SC138] 
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!   To perform emissions of SO2 from explosive volcanic eruptions
!   into the stratosphere. See www.start.or.th for the source of the data.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      SUBROUTINE ukca_volcanic_so2                                      &
                (so2_mmr, mass, latitude, longitude, dellat, dellon,    &
                 row_length, rows, model_levels,                        &
                 year, month, day, timestep, tropopause_height,         &
                 r_theta_levels)

      USE earth_constants_mod, ONLY: earth_radius

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE parcons_mod, ONLY: g, pi, circ, zpi, rad, deg, r
      IMPLICIT NONE


! input parameters:

      INTEGER, INTENT(IN) :: row_length, rows, model_levels

      INTEGER, INTENT(IN) :: year, month, day
      REAL, INTENT(IN) :: timestep, dellat, dellon

! input fields
      REAL, INTENT(IN) :: longitude(row_length, rows)
      REAL, INTENT(IN) :: latitude(row_length, rows)
      REAL, INTENT(IN) :: tropopause_height(row_length, rows)
      REAL, INTENT(IN) :: mass(row_length, rows, model_levels)
      REAL, INTENT(IN) :: r_theta_levels(row_length,rows,model_levels)

! IN/OUT field
! Mass mixing ratio of SO2 in kg(S)/kg(air)
      REAL, INTENT(INOUT) :: so2_mmr(row_length, rows, model_levels)

! local variables
      REAL :: time

! The following is from www.start.or.th/aerosol/emiss2.html
! References: Evans and Kerr, 1983, Krueger, 1990, Bluth et al, 1992,
! Doiron et al, 1991

      INTEGER, PARAMETER :: nerupt = 6 ! no. of volcanic eruptions considered

      REAL :: lat_volc_rad(nerupt)
      REAL :: lon_volc_rad(nerupt)

! Latitude of eruptive volcanoes:
! Mt St Helens, El Chichon, Nevado dei Ruiz, Pinatubo, Cerro Hudson, Fantasy
! 1980/5/18   ,1982/3/28-4/4, 1985/11/13,   1991/6/15, 1991/8/12-8/15 1979/9/1
! latitude / degrees of volcanoes
      REAL, PARAMETER :: lat_volc(nerupt) =                             &
                      (/  46.,   17.3,  4.89,  15.1,  -46.,  0./)
! longitude / degrees of volcanoes
      REAL, PARAMETER :: lon_volc(nerupt) =                             &
                      (/ 238.,  266.8,284.63, 120.4,  287.,  0./)
! year of eruption
      INTEGER, PARAMETER :: year_volc(nerupt) =                         &
                      (/1980,  1982, 1985, 1991,1991, 1979 /)
! day of year when eruption started
      INTEGER, PARAMETER :: start_volc(nerupt) =                        &
                      (/139 ,    87,  255,  166, 223, 241 /)
! day of year when eruption ended
      INTEGER, PARAMETER :: end_volc(nerupt)   =                        &
                      (/139 ,    94,  255,  166, 227, 244 /)
! magnitude in kg(SO2) injected into the stratosphere
      REAL, PARAMETER :: magnitude(nerupt)  =                           &
                      (/0.6e9, 13.4e9, .66e9, 20.e9, 1.5e9, 10.e9/)
! top of eruptive plume, in m. The data is a guestimate, based on the
! SPARC aerosol climatology
      REAL, PARAMETER :: height_top(nerupt) =                           &
                      (/22000.,32000.,26000.,29000.,22000.,30000./)
! Pinatubo height changed from 28 to 29 km as per Oman et al (2006).

      INTEGER, SAVE :: lat_index(nerupt) = (/-1,-1,-1,-1,-1,-1/)
      INTEGER, SAVE :: lon_index(nerupt) = (/-1,-1,-1,-1,-1,-1/)

      REAL, SAVE :: emission(nerupt) ! emission in kg(S)/timestep

      LOGICAL, SAVE :: first = .TRUE.

      REAL, PARAMETER :: day_in_secs = 86400.

! local variables
      REAL :: day_of_year
      INTEGER :: i,j,k,l1,l2

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! determine locations on grid where volcanoes are located
      IF (lhook) CALL dr_hook('UKCA_VOLCANIC_SO2',zhook_in,zhook_handle)
      IF (first) THEN

! revert from degrees to radians for latitude and longitude of volcanoes

        lat_volc_rad = lat_volc * rad
        lon_volc_rad = lon_volc * rad

        DO k=1,nerupt
          DO i=1,row_length
            DO j=1,rows
              IF ((ABS(latitude(i,j)-lat_volc_rad(k))<0.5*dellat) .AND.  &
                 (ABS(longitude(i,j)-lon_volc_rad(k))<0.5*dellon)) THEN
                lat_index(k) = j
                lon_index(k) = i
              END IF
            END DO
          END DO
        END DO

        emission = magnitude / (end_volc - start_volc + 1.)             &
                     * timestep / day_in_secs
! above is as applied to UKCA SO2 which has units of kgSO2/kgair

        first = .FALSE.
      END IF !first

      DO k=1,nerupt
! check whether in processor domain of eruptive volcano
        IF ((lat_index(k) > 0) .AND. (year == year_volc(k))) THEN
          day_of_year = (month-1)*30 + day
! check whether in eruptive phase
          IF ((day_of_year >= start_volc(k)) .AND.                      &
              (day_of_year <= end_volc(k))) THEN
            i = lon_index(k)
            j = lat_index(k)
            l1 = 1
            DO WHILE (r_theta_levels(i,j,l1) - earth_radius <= 19000.)
! changed from tropopause to 19km to match Oman et al 2006.
              l1 = l1 + 1
            END DO
            l2 = model_levels
            DO WHILE (r_theta_levels(i,j,l2) - earth_radius >           &
                height_top(k))
              l2 = l2 - 1
            END DO
            IF (l2 < l1) l1 = l2
! express emission as increase in average increase in mass mixing ratio
            so2_mmr(i,j,l1:l2) = so2_mmr(i,j,l1:l2) +                   &
                emission(k)/SUM(mass(i,j,l1:l2))
          END IF
        END IF
      END DO

      IF (lhook) CALL dr_hook('UKCA_VOLCANIC_SO2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ukca_volcanic_so2
