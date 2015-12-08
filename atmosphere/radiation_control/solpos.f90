! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!   Subroutine SOLPOS   ----------------------------------------------
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Radiation Control

!   Purpose :
!    Calculations of the earth's orbit described in the first page of
!    the "Calculation of incoming insolation" section of UMDP 23, i.e.
!    from the day of the year (and, in forecast mode, whether it is a
!    leap year) and the orbital "constants" (which vary over
!    "Milankovitch" timescales) it calculates the sin of the solar
!    declination and the inverse-square scaling factor for the solar
!    "constant".  It is thus intrinsically scalar.  The FORTRAN code
!    present depends on whether *DEF CAL360 is set during UPDATE: this
!    replaces the Julian calendar with the climate-mode 360-day calendar

!   Programming standard:
!      Written in FORTRAN90 to comply with Version 7.2 of the UMDP3
!       dated the 5/2/98.
!   Logical components covered : P233

!   Project task :

!   External documentation: P23

!   ------------------------------------------------------------------

      SUBROUTINE solpos (day, year, lcal360, l_sec_var, l_eqt           &
                       , eqt, sindec, scs )


      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE yearlen_mod, ONLY: tropyearlength
      IMPLICIT NONE

      LOGICAL, INTENT(IN) :: lcal360    !True if 360 day calendar in use
      LOGICAL, INTENT(IN) :: l_eqt      !True, include equation of time
      LOGICAL, INTENT(IN) :: l_sec_var  !True, include secular variation
!                                       !      of the Earth's orbit
      INTEGER, INTENT(IN) :: day        ! Day-number in the year
      INTEGER, INTENT(IN) :: year       ! Calendar year

      REAL, INTENT(OUT)   :: sindec     ! Sin(solar declination)
      REAL, INTENT(OUT)   :: scs        ! Solar constant scaling factor
      REAL, INTENT(OUT)   :: eqt        ! The equation of time,
!                                       !  specified as an hour angle
!                                       !  in radians.

!  This routine has no dynamically allocated work areas and no
!   significant structure.  It calls the intrinsic functions FLOAT, SIN
!   & COS, and the user subroutine ORBPRM
!  The user deck <astron/astron.h> is included.


!     Mathematical constants:
      REAL, PARAMETER :: twopi = 2. * pi

!     Parameters of the Earth's orbit:
      REAL :: e                        ! Eccentricity of the orbit
      REAL :: gamma                    ! Supplement of the longitude
!                                      !  of the perihelion
      REAL :: tau0                     ! Time of the perihelion
!                                      !  passage in days
      REAL :: oblq                     ! Obliquity of the orbit
      REAL :: diny                     ! Length of the calendar year
!                                      !  (in whole days)

!     Derived orbital constants:
      REAL :: e1
      REAL :: e2
      REAL :: e3
      REAL :: e4
      REAL :: y

      REAL :: m       ! Mean anomaly: positional angle of a "mean" Earth
!                     !  rotating around the sun with a constant angular
!                     !  speed equal to 2pi/T and counted counterclock-
!                     !  wise from the perihelion
      REAL :: v       ! True anomaly: positional angle of Earth in its
!                     !  orbit, counted countercloskwise from perihelion

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('SOLPOS',zhook_in,zhook_handle)

!     Determine the orbital parameters for the calendar selected.
! DEPENDS ON: orbprm
      CALL orbprm(l_sec_var, year, lcal360                              &
                , e, gamma, oblq, tau0, diny )

!     Calculate the mean anomaly at 12Z on the current day.
!     The 0.5 accounts for the time in days since mid-night.
!     The references are to Smart 1944 (and UMDP23)
!     Eq 67 p. 113 and n=2pi/orbital period     (Eq 3.1.1)

      IF (lcal360) THEN
        m = (twopi / diny)           * (FLOAT(day) - tau0 - .5)
      ELSE
        m = (twopi / tropyearlength) * (FLOAT(day) - tau0 - .5)
      END IF

!       Calculate the coefficients in the equation of the centre and
!        thence derive the true anomaly.
      e1 = e * ( 2. - .25 * e*e )
      e2 = 1.25 * e*e
      e3 = e*e*e * 13. / 12.

!       True anomaly, equation 87 in Smart on p. 120 (UMDP23 Eq 3.1.2)
      v  = m + e1*SIN(m) + e2*SIN(2.*m) + e3*SIN(3.*m)

!       Solar constant scaling factor (UMDP23 Eq 3.1.4)
      e4  = ( (1. + e*e*.5) / (1. - e*e) )**2
      scs = e4 * ( 1. + e * COS(v) ) **2

!       sin(solar declination) (UMDP23 Eq 3.1.5)
!       The solar declination is related to
!        the true longitude of the earth (lambda) by:
!        sindec = sin(obliquity) * sin(lambda)
!       Lambda is counted counterclockwise from the vernal equinox
!        and is related to v (the true anomaly) through
!        lambda = v + (longitude of perihelion)

      sindec = SIN(oblq) * SIN (v - gamma)


!     Calculate the equation of time as given by equation (29)
!     on page 149 of Smart (1944). (Recall the factor of 5/4 in the
!     definition of E2).
      IF (l_eqt) THEN
        y   = ( TAN ( 0.5*oblq ) )**2
        eqt = y * SIN(2.0 * ( m - gamma ))                              &
              - 2.0*e * SIN(m)                                          &
              + 4.0*e*y * SIN(m) * COS(2.0*( m - gamma ))               &
              - 0.5*y*y * SIN(4.0*( m - gamma ))                        &
              - e2 * SIN(2.0*m)
      ELSE
        eqt=0.0e+00
      END IF

      IF (lhook) CALL dr_hook('SOLPOS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE solpos
