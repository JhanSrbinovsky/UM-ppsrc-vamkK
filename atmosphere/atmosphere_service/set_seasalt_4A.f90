! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Subroutine to calculate sea-salt aerosol number concentration.
!
! Purpose:
!   Sea-salt aerosol number concentration in 2 modes is calculated
!   from windspeed, height above surface and land fraction.
!
! Method:
!   10m windspeed is used to
!   calculate sea-salt aerosol number concentration in 2 modes
!   following O'Dowd et al. (1999). An exponential decay with
!   height of these concentrations is then assumed, and they are
!   scaled depending on the fraction of open water in the gridbox.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Atmosphere Service
!
! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.
!
!- ---------------------------------------------------------------------

      MODULE set_seasalt_mod

      IMPLICIT NONE
      CONTAINS

      SUBROUTINE set_seasalt_4A(                                        &
           windspeed_10m,                                               &
           height,                                                      &
           land_fract,                                                  &
           sea_ice,                                                     &
           row_length,                                                  &
           num_rows,                                                    &
           num_levels,                                                  &
           salt_dim1, salt_dim2, salt_dim3,                             &
           n_levels_bl,                                                 &
           sea_salt_film,                                               &
           sea_salt_jet                                                 &
         )

      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

      INTEGER                                                           &
           row_length,                                                  &
!             number of points in a row
           num_rows,                                                    &
!             number of rows
           num_levels,                                                  &
!             number of levels
           n_levels_bl,                                                 &
!             number of boundary-layer levels
           salt_dim1, salt_dim2, salt_dim3
!             Dimensions of sea-salt arrays
!
      REAL :: windspeed_10m(row_length, num_rows)
!             10m windspeed (ms-1)
      REAL                                                              &
           land_fract(row_length, num_rows),                            &
!             land fraction
           sea_ice(row_length, num_rows),                               &
!             sea-ice fraction
           height(row_length, num_rows, num_levels)
!             height of theta level centres above sea surface (m)
!
      REAL                                                              &
           sea_salt_film(salt_dim1, salt_dim2, salt_dim3),              &
!             number concentration of sea-salt aerosol derived
!             from film-mode droplets (m-3)
           sea_salt_jet(salt_dim1, salt_dim2, salt_dim3) 
!             number concentration of sea-salt aerosol derived
!             from jet-mode droplets (m-3)
!
!     Local variables:
!
      INTEGER :: i, j, k
!             Loop counters

      REAL, PARAMETER :: film_scale_height = 900.0
!             Scale height for film-mode aerosol (m)
      REAL, PARAMETER :: jet_scale_height = 900.0
!             Scale height for jet-mode aerosol (m)

      REAL :: p2
      REAL :: p3

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('SET_SEASALT_4A',zhook_in,zhook_handle)

      p2 = EXP(10.0/film_scale_height)
      p3 = EXP(10.0/jet_scale_height)

! Parameters used: P2, P3
!$OMP  PARALLEL DO DEFAULT(NONE) SHARED(num_levels, num_rows,           &
!$OMP& row_length, land_fract, height, sea_salt_film,                   &
!$OMP& sea_salt_jet, n_levels_bl, sea_ice, windspeed_10m,P2,P3)         &
!$OMP& PRIVATE(i, j, k)                                                 &
!$OMP& SCHEDULE(STATIC)
      DO k = 1, num_levels
        DO j = 1, num_rows
          DO i = 1, row_length

            IF (land_fract(i, j) <  1.0 .AND. sea_ice(i, J) <  1.0      &
                                        .AND. k <= n_levels_bl) THEN
               IF (windspeed_10m(i, j)  <   2.0) THEN
                 sea_salt_film(i, j, k)                                 &
                             =3.856E06*(1.0-EXP(-0.736*windspeed_10m(i, j)))
                 sea_salt_jet(i, j, k)                                  &
                             =0.671E06*(1.0-EXP(-1.351*windspeed_10m(i, j)))
               ELSE IF (windspeed_10m(i, j)  >   17.5) THEN
                 sea_salt_film(i, j, k)=150.0e06*(1.0-(97.874           &
                                          *EXP(-0.313*windspeed_10m(i, j))))
                 sea_salt_jet(i, j, k)=3.6e06*(1.0-(103.926             &
                                          *EXP(-0.353*windspeed_10m(i, j))))
               ELSE
                 sea_salt_film(i, j, k)                                 &
                                =10.0**((0.095*windspeed_10m(i, j))+6.283)
                 sea_salt_jet(i, j, k)                                  &
                                =10.0**((0.0422*windspeed_10m(i, j))+5.7122)
               ENDIF

               sea_salt_film(i, j, k)=sea_salt_film(i, j, k)            &
                        *p2*EXP(-height(i, j, k)/film_scale_height)     &
                        *(1.0-land_fract(i, j))*(1.0-sea_ice(i, j))

               sea_salt_jet(i, j, k)=sea_salt_jet(i, j, k)              &
                        *p3*EXP(-height(i, j, k)/jet_scale_height)      &
                        *(1.0-land_fract(i, j))*(1.0-sea_ice(i, j))

            ELSE

               sea_salt_film(i, j, k)=0.0
               sea_salt_jet (i, j, k)=0.0

            ENDIF

          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

      IF (lhook) CALL dr_hook('SET_SEASALT_4A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SET_SEASALT_4A

      END MODULE
