! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to calculate weights for ultraviolet fluxes.

! Purpose:
!   Weights to calculate the flux in a wavelength interval [a,b].
!   This is a straightforward extension of the subroutine
!   'R2_SET_690NM_WEIGHT'.

! Method:
!   Straightforward. The flux is assumed to be linearly distributed
!   across bands.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.

!----------------------------------------------------------------------
      SUBROUTINE r2_set_uv_weight(n_band                                &
         , l_present                                                    &
         , n_band_exclude, index_exclude                                &
         , wave_length_short, wave_length_long                          &
         , weight_uv                                                    &
         , npd_band_sw, npd_exclude_sw, npd_type_sw                     &
         )



      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!
!    Define the upper and lower boundaries of the interval in which the
!    fluxes should be calculated.
!
! Method:
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:
!
! Global variables (#include statements etc):
!

!
!    Define the upper and lower boundaries of the interval in which the
!    fluxes should be calculated.
!
        Real, Parameter :: UV_INTERVAL_SHORT = 2.0E-07
        Real, Parameter :: UV_INTERVAL_LONG  = 3.2E-07

!     DUMMY VARIABLES:

!     DIMENSIONS OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
           npd_band_sw                                                  &
!             MAXIMUM NUMBER OF SPECTRAL BANDS
         , npd_exclude_sw                                               &
!             MAXIMUM NUMBER OF EXCLUDED REGIONS
         , npd_type_sw
!             MAXIMUM NUMBER OF TYPES OF SPECTRAL DATA

!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
           n_band                                                       &
!             NUMBER OF SPECTRAL BANDS
         , n_band_exclude(npd_band_sw)                                  &
!             NUMBER OF EXCLUDED REGIONS IN BANDS
         , index_exclude(npd_exclude_sw, npd_band_sw)
!             INDICES OF EXCLUDED REGIONS IN BANDS

      LOGICAL                                                           &
                !, INTENT(IN)
           l_present(0: npd_type_sw)
!             FLAG FOR TYPES OF SPECTRAL DATA PRESENT
      REAL                                                              &
                !, INTENT(IN)
           wave_length_short(npd_band_sw)                               &
!             SHORT WAVELENGTH LIMITS OF BANDS
         , wave_length_long(npd_band_sw)
!             LONG WAVELENGTH LIMITS OF BANDS


!     WEIGHTS SET.
      REAL                                                              &
                !, INTENT(OUT)
           weight_uv(npd_band_sw)
!             WEIGHTS APPLYING TO EACH BAND

!     LOCAL VARIABLES.
      INTEGER                                                           &
           i                                                            &
!             LOOP VARIABLE
         , j
!             LOOP VARIABLE

      REAL                                                              &
           total_energy_range                                           &
!             TOTAL RANGE OF ENERGIES COVERED BY BAND
         , energy_range_uv
!             RANGE OF ENERGIES IN BAND IN THE WAVELENGHT INTERVAL

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



      IF (lhook) CALL dr_hook('R2_SET_UV_WEIGHT',zhook_in,zhook_handle)
      weight_uv=0.0e+00

      DO i=1, n_band

         IF (wave_length_long(i) <  uv_interval_short)                  &
                                         weight_uv(i)=0.0e+00

         IF (wave_length_short(i) >  uv_interval_long)                  &
                                         weight_uv(i)=0.0e+00

         IF (wave_length_short(i) <  uv_interval_short) THEN

            total_energy_range=1.0e+00/wave_length_short(i)             &
                              -1.0e+00/wave_length_long(i)

            IF (wave_length_long(i) >  uv_interval_long) THEN
               energy_range_uv=1.0/uv_interval_short                    &
                              -1.0/uv_interval_long
               IF (l_present(14)) THEN
!                REMOVE CONTRIBUTIONS FROM EXCLUDED BANDS.
                 DO j=1, n_band_exclude(i)
                   IF (wave_length_short(index_exclude(j, i)) <         &
                      uv_interval_short) THEN
                      energy_range_uv=energy_range_uv                   &
                       -1.0e+00/uv_interval_short                       &
                       +1.0e+00/wave_length_long(index_exclude(j, i))
                   ELSE IF (wave_length_long(index_exclude(j, i)) >     &
                      uv_interval_long) THEN
                      energy_range_uv=energy_range_uv                   &
                       -1.0e+00/wave_length_short(index_exclude(j, i))  &
                       +1.0e+00/uv_interval_long
                   ELSE
                      energy_range_uv=energy_range_uv                   &
                       -1.0e+00/wave_length_short(index_exclude(j, i))  &
                       +1.0e+00/wave_length_long(index_exclude(j, i))
                   END IF
                   total_energy_range=total_energy_range                &
                      -1.0e+00/wave_length_short(index_exclude(j, i))   &
                      +1.0e+00/wave_length_long(index_exclude(j, i))
                 END DO
               END IF
              weight_uv(i)=energy_range_uv/total_energy_range
            END IF

            IF ((wave_length_long(i) <  uv_interval_long).AND.          &
                (wave_length_long(i) >  uv_interval_short)) THEN
               energy_range_uv=1.0e+00/uv_interval_short                &
                              -1.0e+00/wave_length_long(i)
               IF (l_present(14)) THEN
!                REMOVE CONTRIBUTIONS FROM EXCLUDED BANDS.
                 DO j=1, n_band_exclude(i)
                   IF (wave_length_short(index_exclude(j, i)) >         &
                      uv_interval_short) THEN
                      energy_range_uv=energy_range_uv                   &
                       -1.0e+00/wave_length_short(index_exclude(j, i))  &
                       +1.0e+00/wave_length_long(index_exclude(j, i))
                   ELSE IF (wave_length_short(index_exclude(j, i)) <    &
                      uv_interval_short) THEN
                      energy_range_uv=energy_range_uv                   &
                       -1.0e+00/uv_interval_short                       &
                       +1.0e+00/wave_length_long(index_exclude(j, i))
                   END IF
                   total_energy_range=total_energy_range                &
                      -1.0e+00/wave_length_short(index_exclude(j, i))   &
                      +1.0e+00/wave_length_long(index_exclude(j, i))
                 END DO
               END IF
               weight_uv(i)=energy_range_uv/total_energy_range
            END IF

         ELSE

            total_energy_range=1.0e+00/wave_length_short(i)             &
                              -1.0e+00/wave_length_long(i)

            IF (wave_length_long(i) <= uv_interval_long)                &
                                           weight_uv(i)=1.0

            IF ((wave_length_short(i) <  uv_interval_long).AND.         &
                (wave_length_long(i) >  uv_interval_long)) THEN
               energy_range_uv=1.0e+00/wave_length_short(i)             &
                              -1.0e+00/uv_interval_long

               IF (l_present(14)) THEN
!                REMOVE CONTRIBUTIONS FROM EXCLUDED BANDS.
                 DO j=1, n_band_exclude(i)
                   IF (wave_length_long(index_exclude(j, i)) <          &
                      uv_interval_long) THEN
                      energy_range_uv=energy_range_uv                   &
                       -1.0e+00/wave_length_short(index_exclude(j, i))  &
                       +1.0e+00/wave_length_long(index_exclude(j, i))
                   ELSE IF (wave_length_short(index_exclude(j, i)) <    &
                       uv_interval_long) THEN
                       energy_range_uv=energy_range_uv                  &
                        -1.0e+00/wave_length_short(index_exclude(j, i)) &
                        +1.0e+00/uv_interval_long
                   END IF
                   total_energy_range=total_energy_range                &
                      -1.0e+00/wave_length_short(index_exclude(j, i))   &
                      +1.0e+00/wave_length_long(index_exclude(j, i))
                 END DO
               END IF
               weight_uv(i)=energy_range_uv/total_energy_range
            END IF

         END IF

      END DO

      IF (lhook) CALL dr_hook('R2_SET_UV_WEIGHT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE r2_set_uv_weight
