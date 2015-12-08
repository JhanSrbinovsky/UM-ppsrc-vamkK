! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to calculate the parameters of the Earth's orbit.

!  Purpose:
!  This routine returns the parameters of the Earth's orbit
!   (the eccentricity, obliquity, supplement of the longitude of
!    perihelion) the time of the parihelion passage in days, and
!   Length of the calendar year (in whole days).

!  Method:
!  For long runs there may be an interest in running with secular
!   variations in the astronomy. The orbital constants have
!   been derived from A. L. Berger 1978, J. Atm. Sci, Volume 35
!   2362-2367. A copy of which can be found in the Met Office library.
!  For short current runs, or long control runs it is preferrable
!   not to allow the astronomy to vary, so fixed values are used.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of Code:
!   FORTRAN90 complying with UMDP3 Version 7.2 from 5/2/98

!- ---------------------------------------------------------------------
      SUBROUTINE orbprm(l_sec_var, year, lcal360                        &
                      , e, gamma, oblq, tau0, diny )

      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE yearlen_mod, ONLY: tropyearlength
      IMPLICIT NONE

      INTEGER, INTENT (IN) :: year       ! Calendar year
      LOGICAL, INTENT (IN) :: l_sec_var  ! Include secular variations
!                                        !  of the orbit
      LOGICAL, INTENT (IN) :: lcal360    ! Use a calendar of 360 days

!     Parameters of the Earth's orbit:

      REAL, INTENT(OUT) :: e             ! Eccentricity of the orbit
      REAL, INTENT(OUT) :: gamma         ! Supplement of the longitude
!                                        !  of the perihelion
      REAL, INTENT(OUT) :: oblq          ! Obliquity of the orbit
      REAL, INTENT(OUT) :: tau0          ! Time of the perihelion
!                                        !  passage in days
      REAL, INTENT(OUT) :: diny          ! Length of the calendar year
!                                        !  (in whole days)

!     Local Variables for use within ORBPRM

      REAL :: year_offset                ! Offset of the year from the
!                                        !  reference year when default
!                                        !  values apply
      REAL :: ecn_sn                     ! Eccentricity multiplied by
!                                        !  the sine of the longitude
!                                        !  of the perihelion
      REAL :: ecn_cn                     ! Eccentricity multiplied by
!                                        !  the cosine of the longitude
!                                        !  of the perihelion
      REAL :: lph_fixed_ve               ! Longitude of the perihelion
!                                        !  relative to a fixed vernal
!                                        !  equinox
      REAL :: gn_prcs                    ! General precession
      REAL :: date_ve                    ! Date of the vernal equinox
!                                        !  in days into the year
      REAL :: no_leap_days               ! The number of leap days,
!                                        !  used to calculate DATE_VE
      REAL :: mean_anom_ve               ! Mean anomaly at the vernal
!                                        !  equinox

!     Synthetic constants

      REAL :: beta
      REAL :: ee1
      REAL :: ee2
      REAL :: ee3

      INTEGER :: i                       ! Loop variable

!     Mathematical constants:
      REAL, PARAMETER :: twopi = 2. * pi

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!     Astronomical Parameters:
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! ----------------------- Header file ASTRON  -------------------------
! Description: Parameters of the Earth's orbit.
!----------------------------------------------------------------------
!
!     Default values of the orbital elements
!      (currently those for the epoch J2000 which is 1.5d Jan. 2000):
!     The Eccentricity and Longitue of perhelion are recommended by NASA
!      see (http://ssd.jpl.nasa.gov/elem_planets.html)
!     The Obliquity value comes from the Astronomical Almanac for 1984
!      page S26 and is used on several webpages e.g.
!      nedwww.ipac.caltech.edu/help/calc_doc.txt
!      www.stargazing.net/kepler/astrovba2.html
!      http://edhs1.gsfc.nasa.gov/waisdata/docsw/txt/tp4450505.txt
!
!     The data in the series expansions are adapted from Berger 1978.
!      Andre' L. Berger Journ. of Atm. Sci. Volume 35 p. 2362-2367,
!      and is available from the Met Office library.
!
!     ! Eccentricity of the orbit
      Real, Parameter    :: E_DFLT         = 1.6710222E-02
!
!     ! Longitude of the perihelion in radians
      Real, Parameter    :: LPH_DFLT       = 102.94719*PI/180.0
!
!     ! Obliquity of the orbit - corresponds to 23.43929111 degrees
      Real, Parameter    :: OBLQ_DFLT      = 0.409092804
!
!     ! Reference year for setting the date of the vernal equinox
      Integer, Parameter :: YEAR_REF_VE    = 2000
!
!     ! Date of the vernal equinox in days after the start of the year
!     !  This date is for the year 2000.
      Real, Parameter    :: DATE_VE_DFLT   = 79.3159
!
!     The final parameter required is the time of the perihelion
!     passage, TAU0. For a pure Keplerian orbit, with a specified
!     eccentricity and longitude of the perihelion, this can be
!     deduced from the date of the vernal equinox (as is specified
!     in AMIP-2, for example). In practice it is somewhat more
!     complicated to calculate the time of the perihelion.
!     For simplicity, a mean value for the years 1995-2005 is used
!     here: note that the range of TAU0 in this period is from
!     1.0 to 3.75 and that there is no simple relationship with leap
!     years.
!
!     ! Time of the perihelion passage in days
      Real, Parameter    :: TAU0_DFLT      = 2.667
!
!     ------------------------------------------------------------------
!
!     The parameters used to calculate secular variations of the orbital
!     elements are taken from A. L. Berger (1978), J. Atm. Sci., vol 35,
!     p. 2362.These have been converted so that:
!     amplitudes are in radians,
!     angular frequencies in radians per year and
!     phases in radians
!     Phases have also been converted to be taken relative to
!     J2000 for consistency with the time for the default values of the
!     orbital parameters.
!     Berger's numbers (with time correction) differ slightly from the
!     default values above.
!
!     The obliquity and longitude of the perihelion have been adjusted
!     to agree with the NASA values for J2000, but adjustment of the
!     eccentricity to NASA values is not so easy and has not been done.
!
!
!     ! Reference year     YEAR_REF
      Integer, Parameter :: YEAR_REF       = 2000
!
!   -----------------------------------------------------------------
!     Obliquity (Table 1) from the first 24 terms:
!     (enough for deviations of less than 0.002 degrees)
!
!     ! Constant term in the obliquity: from the Astr. Almanac for 1984
!     !  The following value corresponds to 23.320870 degrees at J2000
      Real, Parameter    :: OBLQ_CNST      = 0.40702597
!
!     ! Number of terms retained in the series for the obliquity
      Integer, Parameter :: N_TERM_OBQ     = 24
!
!     ! Amplitude
      Real :: A(N_TERM_OBQ)
!     ! Angular frequency
      Real :: F(N_TERM_OBQ)
!     ! Phase in the series
      Real :: D(N_TERM_OBQ)
!
!   -----------------------------------------------------------------
!     Eccentricity and longitude of the fixed perihelion (Table 4):
!
!     ! Number of terms retained in the series for the
!     !  eccentricty and longitude of the perihelion
      Integer, Parameter :: N_TERM_ECN_LPH = 19
!
!     ! Amplitude
      Real :: M(N_TERM_ECN_LPH)
!     ! Angular frequency
      Real :: G(N_TERM_ECN_LPH)
!     ! Phase in the series
      Real :: B(N_TERM_ECN_LPH)
!
!   ------------------------------------------------------------------72
!     General Precession (Table 5):
!
!     ! Linear rate of precession!
!     ! The value corresponds to 50.439273 seconds per year -Berger 1979
      Real, Parameter :: LIN_RATE_GN_PRCS  = 2.44536496E-04
!
!     ! Constant offset to general precession (in seconds pre year),
!     ! corrected for 50 years difference in reference time.
      Real, Parameter :: GN_PRCS_CNST      = 7.14372244E-02
!
!     ! Number of terms kept in the series for the general precession
      Integer, Parameter :: N_TERM_GN_PRCS = 10
!
!     ! Amplitude
      Real :: C(N_TERM_GN_PRCS)
!     ! Angular frequency
      Real :: H(N_TERM_GN_PRCS)
!     ! Phase in the series
      Real :: R(N_TERM_GN_PRCS)
!
!   -----------------------------------------------------------------
!    Table 1
!
      DATA A/                                                           &
     &    -1.19372E-02, -4.15640E-03, -3.05103E-03, -2.00849E-03        &
     &  , -1.51146E-03,  1.49778E-03, -7.88065E-04, -5.62917E-04        &
     &  ,  4.90244E-04, -3.28170E-04,  1.20767E-04,  1.09471E-04        &
     &  , -1.02587E-04, -7.58733E-05,  7.46128E-05,  7.11222E-05        &
     &  , -5.68686E-05,  4.97904E-05,  3.14644E-05,  2.83616E-05        &
     &  , -2.66163E-05, -2.63254E-05,  2.50164E-05,  2.46285E-05/
      DATA F/                                                           &
     &     1.5324946E-04,  1.5814864E-04,  1.1719011E-04                &
     &  ,  1.5506174E-04,  2.1733392E-04,  1.5016256E-04                &
     &  ,  2.1170962E-04,  1.5633636E-04,  1.4835028E-04                &
     &  ,  2.0692488E-04,  2.1252514E-04,  2.2999289E-04                &
     &  ,  3.0649899E-04,  3.1139817E-04,  4.8991877E-06                &
     &  ,  3.6059331E-05,  2.7043965E-04,  1.8122966E-06                &
     &  ,  6.4084427E-05,  3.0341210E-04,  3.0831127E-04                &
     &  ,  3.7058338E-04,  2.2211866E-04,  4.0958519E-05/
      DATA D/                                                           &
     &     4.4041E+00,  4.9093E+00,  2.2451E+00,  5.1167E+00            &
     &  ,  2.7912E-01,  4.6115E+00,  5.3935E+00,  4.1966E+00            &
     &  ,  3.8990E+00,  4.7014E+00,  5.5397E+00,  5.5896E+00            &
     &  ,  2.5251E+00,  3.0303E+00,  5.0517E-01,  2.1589E+00            &
     &  ,  3.6608E-01,  7.1253E-01,  2.1582E+00,  2.7325E+00            &
     &  ,  3.2376E+00,  4.6833E+00,  9.7121E-01,  2.6640E+00/
!
!   -----------------------------------------------------------------
!    Table 4
!
      DATA M/                                                           &
     &     1.8607980E-02,  1.6275220E-02, -1.3006600E-02                &
     &  ,  9.8882900E-03, -3.3670000E-03,  3.3307700E-03                &
     &  , -2.3540000E-03,  1.4001500E-03,  1.0070000E-03                &
     &  ,  8.5700000E-04,  6.4990000E-04,  5.9900000E-04                &
     &  ,  3.7800000E-04, -3.3700000E-04,  2.7600000E-04                &
     &  ,  1.8200000E-04, -1.7400000E-04, -1.2400000E-04                &
     &  ,  1.2500000E-05/
      DATA G/                                                           &
     &     2.0397105E-05,  3.5614854E-05,  8.6574454E-05                &
     &  ,  8.3487563E-05,  8.1675266E-05,  2.5205846E-05                &
     &  ,  8.8386751E-05,  1.2710243E-04,  3.0830121E-05                &
     &  ,  7.8588375E-05,  1.4860417E-05,  8.0400672E-05                &
     &  ,  8.9661345E-05,  3.0014587E-05,  9.1473642E-05                &
     &  ,  8.4481533E-05,  2.9990579E-05,  8.9290274E-05                &
     &  ,  3.2378912E-06/
      DATA B/                                                           &
     &     5.0053E-01,  3.3839E+00,  5.3852E+00,  5.5925E+00            &
     &  ,  4.8800E+00,  1.5230E+00,  6.0977E+00,  2.2481E+00            &
     &  ,  2.6918E+00,  5.0874E+00,  2.0054E+00,  5.8001E+00            &
     &  ,  5.1778E+00,  2.5455E+00,  5.8903E+00,  2.6587E+00            &
     &  ,  2.2151E+00,  3.6812E+00,  1.2585E+00/
!
!   -----------------------------------------------------------------
!    Table 5
!
      DATA C/                                                           &
     &     3.58327E-02,  1.23877E-02,  9.80662E-03, -9.56853E-03        &
     &  ,  6.01280E-03,  4.62449E-03, -4.51725E-03,  4.22942E-03        &
     &  ,  2.93967E-03, -2.40482E-03/
      DATA H/                                                           &
     &     1.5324946E-04,  1.5814864E-04,  1.1719011E-04,  3.0868911E-06&
     &  ,  1.5506174E-04,  1.5217749E-05,  1.5016256E-04,  2.1733392E-04&
     &  ,  4.8087409E-06,  1.8122966E-06/
      DATA R/                                                           &
     &     4.4041E+00,  4.9093E+00,  2.2451E+00,  6.0756E+00            &
     &  ,  5.1167E+00,  2.8833E+00,  4.6115E+00,  2.7912E-01            &
     &  ,  1.0225E+00,  7.1253E-01/
!

!     The length of the calendar year may be set for a 360-day calendar
!      (as is often used in climate runs),
!      or for a real Gregorian calendar which has 365 days in
!      non-leap years and 366 in leap years.

      IF (lhook) CALL dr_hook('ORBPRM',zhook_in,zhook_handle)
      IF (lcal360) THEN

        diny=360.0

      ELSE
!      Is this a leap year?
        IF (MOD(year,4)    ==  0 .AND.                                  &
           (MOD(year,400)  ==  0 .OR. MOD(year,100)  /=  0)) THEN

          diny = 366.0

!      Is this a normal year?
        ELSE

          diny = 365.0

        END IF
      END IF

!     The orbital elements are normally set to default values, but
!     secular variations may be required in some longer climate runs.

      IF (l_sec_var) THEN

        year_offset = REAL( year - year_ref )

!       Obliquity: (Equation 1 from Berger 1978)

        oblq = oblq_cnst
        DO i=1, n_term_obq
          oblq = oblq+a(i)*COS(f(i)*year_offset+d(i))
        END DO

!       Eccentricity: this is better computed from its components
!       than directly from the series.(Equation (4) of Berger 1978).

        ecn_sn = m(1) * SIN (g(1) * year_offset + b(1))
        ecn_cn = m(1) * COS (g(1) * year_offset + b(1))

        DO i=2, n_term_ecn_lph
          ecn_sn = ecn_sn + m(i) * SIN (g(i) * year_offset + b(i))
          ecn_cn = ecn_cn + m(i) * COS (g(i) * year_offset + b(i))
        END DO
        e = SQRT(ecn_sn*ecn_sn+ecn_cn*ecn_cn)

!       We now obtain the longitude of the perihelion relative to the
!       fixed equinox.

        lph_fixed_ve = ATAN2 (ecn_sn,ecn_cn)

!       The longitude of perihelion and
!        the supplement of the longitude of the perihelion relative to
!        the actual vernal equinox requires the general precession.

!      General Precession.
        gn_prcs = lin_rate_gn_prcs * year_offset + gn_prcs_cnst
        DO i=1, n_term_gn_prcs
          gn_prcs = gn_prcs + c(i) * SIN (h(i) * year_offset + r(i))
        END DO

!      Supplement of the longitude of the perihelion
        gamma = pi - lph_fixed_ve - gn_prcs

!       Time of perihelion: The time at which an object is at perihelion
!        (its closest distance to the sun).
!       The time of perihelion is inferred from the date of
!        the vernal equinox using the Gregorian calendar.

!      Calculate the date of the vernal equinox.
!       First we need to:
!        Calculate the no of leap days between year & year_ref_ve.
!        This needs to be corrected when using the Gregorian calendar.
!         by adding (DINY-366.0) when the year_ref_ve is a leap year or
!         by adding (DINY-365.0) when the year_ref_ve is a normal year.
!        This correction is done when the DATE_VE is calculated below!

!        In the calculation of NO_LEAP_DAYS below, the divisions of type
!         'YEAR'/x (where x is 4, 100 or 400) are integer computations.
!         These integers are then subtracted and the resulting integer
!         is then converted to a real.

        no_leap_days = ( tropyearlength - 365.0)                        &
          * REAL( year     - year_ref_ve     )                          &
          - REAL( year/4   - year_ref_ve/4   )                          &
          + REAL( year/100 - year_ref_ve/100 )                          &
          - REAL( year/400 - year_ref_ve/400 )

!      Now we can calculate the date of the vernal equinox!
!      Because the date of the vernal equinox is varying with the year,
!      we have to keep track of its position in the sky.
!      In order to accomodate a time varying vernal equinox when using
!      a 360-day year, we still have to calculate the difference in
!      the vernal equinox depending on leap years, normal years and the
!      difference between the length of the tropical year and the
!      "normal" year and then we adjust this by multiplying the
!      DATE_VE by 360/(length of tropical year).

!      Is a 360 day calendar being used?

        IF (lcal360) THEN
          date_ve = date_ve_dflt + no_leap_days
          date_ve = date_ve * diny / tropyearlength


!      Is a 365 day calendar being used?
        ELSE

!        Is the epoch reference year a leap year?

          IF (MOD(year_ref_ve,4)    ==  0 .AND.                         &
             (MOD(year_ref_ve,400)  ==  0 .OR.                          &
              MOD(year_ref_ve,100)  /=  0)) THEN

            date_ve = date_ve_dflt + (no_leap_days + (diny - 366.0))

!        Is the epoch reference year a normal year?

          ELSE

            date_ve = date_ve_dflt + (no_leap_days + (diny - 365.0))

          END IF
        END IF

        beta = SQRT(1.0e+00-e*e)
        ee1  = (0.5*e + 0.125*e*e*e)*(1.0 + beta)
        ee2  = -0.25*e*e* (0.5 + beta)
        ee3  = 0.125*e*e*e*((1.0/3.0) + beta)
        mean_anom_ve = gamma - 2.0e+00 * (                              &
            ee1 * SIN (gamma)                                           &
          + ee2 * SIN (2.0 * gamma)                                     &
          + ee3 * SIN (3.0 * gamma)                                     &
          )

        tau0 = date_ve - mean_anom_ve * tropyearlength/(twopi)

      ELSE

        e     = e_dflt
        oblq  = oblq_dflt
        gamma = pi - lph_dflt
        tau0  = tau0_dflt

      END IF

!     If using a 360-day calendar the time of the perihelion is
!     adjusted.
      IF (lcal360) THEN
        tau0 = tau0*(360.0/tropyearlength)+0.71
      END IF



      IF (lhook) CALL dr_hook('ORBPRM',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE orbprm
