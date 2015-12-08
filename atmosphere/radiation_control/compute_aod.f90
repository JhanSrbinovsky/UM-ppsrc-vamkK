! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to compute the aerosol optical depth (AOD).

! Method:
!       The optical depth is the vertical integral of the aerosol
!       specific extinction weighted by the aerosol mass.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of Code:
!   FORTRAN 77

!-----------------------------------------------------------------------
      SUBROUTINE compute_aod (                                          &
! actual and fixed array dimensions
         n_aerosol, n_aerosol_mr, npd_aerosol, npd_aerosol_mr,          &
         n_aod_wavel, npd_aod_wavel,                                    &
         n_profile, npd_profile,                                        &
         n_layer, first_layer, npd_layer,                               &
         n_humidities, npd_humidities,                                  &
! variables with intent in
         type_aerosol,                                                  &
         i_aod_type, i_type_wanted, l_radn_on, l_prognostic,            &
         aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,   &
         d_mass, i_aerosol_parametrization,                             &
         aod_absorption, aod_scattering,                                &
         i_humidity_pointer, mean_rh,                                   &
         humidities, delta_humidity,                                    &
! variable with intent out
         aod                                                            &
       )

      USE rad_pcf
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE aodtype_mod, ONLY: ip_type_allaod, ip_type_sulphate,          &
                             ip_type_dust, ip_type_seasalt,             &
                             ip_type_soot, ip_type_biomass,             &
                             ip_type_biogenic, ip_type_ocff,            &
                             ip_type_delta, ip_type_nitrate,            &
                             ip_type_twobdust
      IMPLICIT NONE

! declaration of arguments

! array dimensions (actual and fixed)

      INTEGER n_aerosol
!  number of aerosol components in spectral information
      INTEGER n_aerosol_mr
!  number of aerosol components in mixing ration information
      INTEGER npd_aerosol 
!  number of aerosol components in the spectral information
      INTEGER npd_aerosol_mr 
!  number of aerosol components in aerosol_mix_ratio

!  number of wavelengths at which the AOD is to be computed
      INTEGER n_aod_wavel
      INTEGER npd_aod_wavel
!  number of grid-boxes
      INTEGER n_profile
      INTEGER npd_profile
!  number of vertical layers
      INTEGER n_layer
      INTEGER first_layer
      INTEGER npd_layer
!  number of humidities used in moist aerosol parameterisation
      INTEGER n_humidities
      INTEGER npd_humidities

! arguments with intent in

!  array giving the type of each aerosol component (see AERCMP3A)
      INTEGER type_aerosol(npd_aerosol)
!  array connecting an aerosol component to an aerosol type
      INTEGER i_aod_type(npd_aerosol)
!  aerosol type to be considered in this call
      INTEGER i_type_wanted
!  aerosol component mass mixing ratio
      REAL aerosol_mix_ratio(npd_profile, first_layer:npd_layer,        &
                             npd_aerosol_mr)
! Index relating aerosol_mix_ratio aerosols to aerosols in
! the spectral information
      INTEGER aerosol_mr_type_index(npd_aerosol_mr)
! Scheme/source of the aerosol data, to determine use in
! changing radiative fluxes and use in diagnostics
      INTEGER aerosol_mr_source(npd_aerosol_mr)
!  mass thickness of vertical layers
      REAL d_mass(npd_profile, npd_layer)
!  aerosol parameterisation (dry or moist)
      INTEGER i_aerosol_parametrization(npd_aerosol)
!  aerosol specific coefficients for absorption and scattering
!  (monochromatic)
      REAL aod_absorption(npd_humidities, npd_aerosol, npd_aod_wavel)
      REAL aod_scattering(npd_humidities, npd_aerosol, npd_aod_wavel)
!  mean relative humidities, and indices to the look-up tables
!  it may be the grid-box mean or clear-sky mean relative humidity,
!  depending on calculations made in FLUX_CALC
      INTEGER i_humidity_pointer(npd_profile, npd_layer)
      REAL mean_rh(npd_profile, npd_layer)
      REAL humidities(npd_humidities, npd_aerosol)
      REAL delta_humidity
!  logical determining if the aod to be calculated is only from 
!  radiatively active aerosol mixing ratio data
      LOGICAL l_radn_on
!  logical determining if the aod to be calculated is only from 
!  prognostic aerosol mixing ratio data, as opposed to climatologies
      LOGICAL l_prognostic

! arguments with intent out

!  computed aerosol optical depth
      REAL aod(npd_profile, npd_aod_wavel)


! local variables

!  loop indices
      INTEGER i, j, jj, j_mr, jj_mr, k, l
!  number of aerosol components included in the current aerosol type
      INTEGER n_aer_in_type
!  indices of those components (at most, all components in the same
!                               type)
      INTEGER n_aer_in_type_index(npd_aerosol_mr)
!  aerosol extinction (absorption + scattering)
      REAL extinction
!  variables needed for the interpolation on humidities
      INTEGER i_pointer
      REAL weight_upper
      REAL weight_lower

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('COMPUTE_AOD',zhook_in,zhook_handle)
      DO k = 1, n_aod_wavel
        DO l = 1, n_profile
          aod(l, k) = 0.0
        END DO
      END DO

      ! FIND WHICH AEROSOL COMPONENT BELONGS TO THE CURRENT TYPE
      ! DO NOT CONSIDER CLIMATOLOGICAL AEROSOLS (TYPE_AEROSOL IS
      ! SMALLER THAN 10) UNLESS i_type_wanted is ALLAOD
      n_aer_in_type = 0
      DO j_mr = 1, n_aerosol_mr
        ! convert the j_mr into a j for spectral information:
        j=aerosol_mr_type_index(j_mr)
        ! test to see if the aerosol data is required
        IF( i_type_wanted == IP_TYPE_ALLAOD .OR.                        &
           (type_aerosol(j) >= 10 .AND.                                 &
            i_aod_type(j) == i_type_wanted) ) THEN
          IF (l_radn_on .AND. .NOT.l_prognostic) THEN
          ! Only want aerosols that are radiatively active
          ! and do not care if they are prognostic or not
            IF (aerosol_mr_source(j_mr) == ip_aersrc_cusack_ron .or.    &
             aerosol_mr_source(j_mr) == ip_aersrc_classic_ron  .or.     &
             aerosol_mr_source(j_mr) == ip_aersrc_arcl_ron) THEN
               n_aer_in_type = n_aer_in_type + 1
               n_aer_in_type_index(n_aer_in_type) = j_mr
            END IF
          ELSE IF (.NOT.l_radn_on .AND. l_prognostic) THEN
          ! only want prognostic aerosol mixing ratio data
            IF (aerosol_mr_source(j_mr) == ip_aersrc_classic_ron .or.   &
             aerosol_mr_source(j_mr) == ip_aersrc_classic_roff) THEN
               n_aer_in_type = n_aer_in_type + 1
               n_aer_in_type_index(n_aer_in_type) = j_mr
              END IF
          ELSE IF (l_radn_on .AND. l_prognostic) THEN
            IF (aerosol_mr_source(j_mr) == ip_aersrc_classic_ron) THEN
              n_aer_in_type = n_aer_in_type + 1
              n_aer_in_type_index(n_aer_in_type) = j_mr
            END IF
          ELSE 
            n_aer_in_type = n_aer_in_type + 1
            n_aer_in_type_index(n_aer_in_type) = j_mr
          END IF
        END IF
      END DO

      ! COMPUTE THE AOD IF AT LEAST ONE COMPONENT MATCHES
      ! (IF NO MATCH, THEN THE AOD WILL REMAIN AT ZERO)
      IF (n_aer_in_type  >   0) THEN

        DO jj = 1, n_aer_in_type
          j_mr = n_aer_in_type_index(jj)
          j=aerosol_mr_type_index(j_mr)

          IF(i_aerosol_parametrization(j)  ==                           &
             ip_aerosol_param_dry) THEN

            ! Non-hygroscopic aerosol

            DO k = 1, n_aod_wavel
              extinction = aod_absorption(1, j, k) +                    &
                           aod_scattering(1, j, k)
              DO i = 1, n_layer
                DO l = 1, n_profile
                  aod(l, k) = aod(l, k) +                               &
                    aerosol_mix_ratio(l, i, j_mr) *                     &
                    d_mass(l, i) * extinction
                END DO ! L
              END DO ! I
            END DO ! K

          ELSE IF(i_aerosol_parametrization(j)  ==                      &
                  ip_aerosol_param_moist) THEN

            ! Hygroscopic aerosol
            ! interpolation on the mean relative humidity

            DO k = 1, n_aod_wavel
              DO i = 1, n_layer
                DO l = 1, n_profile
                  i_pointer = i_humidity_pointer(l, i)
                  weight_upper = ( mean_rh(l, i)                        &
                         - humidities(i_pointer, j))                    &
                         / delta_humidity
                  weight_lower = 1.00e+00 - weight_upper

                  extinction =                                          &
                    (aod_absorption(i_pointer, j, k) +                  &
                     aod_scattering(i_pointer, j, k))                   &
                    * weight_lower + weight_upper *                     &
                    (aod_absorption(i_pointer+1,j,k) +                  &
                     aod_scattering(i_pointer+1,j,k))

                  aod(l, k) = aod(l, k) +                               &
                    aerosol_mix_ratio(l, i, j_mr) *                     &
                    d_mass(l, i) * extinction
                END DO ! L
              END DO ! I
            END DO ! K
          END IF

        END DO ! JJ

      END IF ! N_AER_IN_TYPE > 0

      IF (lhook) CALL dr_hook('COMPUTE_AOD',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE compute_aod
