! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! A wrapper subroutine to calculate the aerosol optical depth diagnostic
! a number of different times for different aerosols

SUBROUTINE compute_all_aod(ierr,                                        &
       n_aerosol,      n_aerosol_mr,                                    &
       nd_aerosol_species, nd_aerosol_mixratio,                         &
       n_aod_wavel,    nd_aod_wavel,    aod_wavel,                      &
       n_profile,      nd_profile,                                      &
       n_layer,        first_layer,     nd_layer,                       &
       n_humidities,   nd_humidities,                                   &
       type_aerosol,   i_aod_type,                                      &
       aerosol_mix_ratio, d_mass,                                       &
       aerosol_mr_source, aerosol_mr_type_index,                        &
       i_aerosol_parametrization, aod_absorption, aod_scattering,       &
       i_humidity_pointer, mean_rel_humidity,                           &
       humidities,     delta_humidity,   l_use_arcl,                    &
       l_aod_sulphate, aod_sulphate, l_aod_dust,      aod_dust,         &
       l_aod_seasalt,  aod_seasalt,  l_aod_soot,      aod_soot,         &
       l_aod_biomass,  aod_biomass,  l_aod_biogenic,  aod_biogenic,     &
       l_aod_ocff,     aod_ocff,     l_aod_delta,     aod_delta,        &
       l_aod_nitrate,  aod_nitrate,  l_aod_total_radn,aod_total_radn,   &
       l_angst_total_radn,   angst_total_radn,                          &
       l_aod_prog_sulphate,  aod_prog_sulphate,   l_aod_prog_dust,      &
       aod_prog_dust,        l_aod_prog_seasalt,  aod_prog_seasalt,     &
       l_aod_prog_soot,      aod_prog_soot,       l_aod_prog_biomass,   &
       aod_prog_biomass,     l_aod_prog_ocff,     aod_prog_ocff,        &
       l_aod_prog_nitrate,   aod_prog_nitrate)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE dust_parameters_mod, ONLY: l_twobin_dust
USE arcl_mod, ONLY: npd_arcl_species, ip_arcl_dust

USE aodtype_mod, ONLY: ip_type_allaod, ip_type_sulphate, ip_type_dust,   &
                       ip_type_seasalt, ip_type_soot, ip_type_biomass,   &
                       ip_type_biogenic, ip_type_ocff, ip_type_delta,    &
                       ip_type_nitrate, ip_type_twobdust
IMPLICIT NONE
!
! Description:
!   Acts as a wrapper subroutine so that compute_aod can be called multiple
!   times for different aerosol types.
!
! Method:
!   Calls the compute_aod subroutine multiple times, for different aerosol
!   types.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: 2: Long wave radiation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!
! Declarations:
!   These are of the form:-
!     INTEGER, INTENT(IN) :: ExampleVariable      !Description of variable
!
! Global variables (#include statements etc):

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!     Dummy arguments.
 INTEGER, INTENT(INOUT) :: IERR      ! Error flag

!     General properties of input data
 INTEGER, INTENT(IN) :: n_aod_wavel  ! Number of wavelengths for aod calcualtions
 INTEGER, INTENT(IN) :: n_aerosol    ! Number of aerosols in spectral information
 INTEGER, INTENT(IN) :: n_aerosol_mr ! Number of aerosols in aerosol_mix_ratio array
 INTEGER, INTENT(IN) :: n_profile    ! Number of Atmospheric profiles
 INTEGER, INTENT(IN) :: n_layer      ! Number of Atmospheric layers
 INTEGER, INTENT(IN) :: first_layer  ! First atmosphere layer
 INTEGER, INTENT(IN) :: n_humidities ! Number of humidities used in computing aod
!
!     Size allocated for...
 INTEGER, INTENT(IN) :: nd_profile   ! ... atmospheric profiles
 INTEGER, INTENT(IN) :: nd_layer     ! ... atmospheric layers
 INTEGER, INTENT(IN) :: nd_aerosol_species  ! ... aerosol species in spectral info
 INTEGER, INTENT(IN) :: nd_aerosol_mixratio ! ... aerosol species in mixing array
 INTEGER, INTENT(IN) :: nd_aod_wavel        ! ... the AOD wavelengths
 INTEGER, INTENT(IN) :: nd_humidities       ! ... humidities for computing AOD

!     Input data about aerosol properties
 INTEGER, INTENT(IN) :: type_aerosol(nd_aerosol_species)
                                         ! Aerosol type in spectral information
 INTEGER, INTENT(IN) :: i_aod_type(nd_aerosol_species)   
                                         ! Map between aerosol type and aod type
 INTEGER, INTENT(IN) :: i_aerosol_parametrization(nd_aerosol_species)
                                             ! Parametrization flags for aerosol
 INTEGER, INTENT(IN) :: aerosol_mr_type_index(nd_aerosol_mixratio)
                         ! Index relating aerosol_mix_ratio aerosols to aerosols in
                         ! the spectral information
 INTEGER, INTENT(IN) :: aerosol_mr_source(nd_aerosol_mixratio)
                         ! Scheme/source of the aerosol data, to determine use in
                         ! changing radiative fluxes and use in diagnostics


 REAL, INTENT(IN) :: aod_absorption(nd_humidities, nd_aerosol_species,  &
                        nd_aod_wavel)    ! Aerosol absorption array for AOD
 REAL, INTENT(IN) :: aod_scattering(nd_humidities, nd_aerosol_species,  &
                        nd_aod_wavel)    ! Aerosol scattering array for AOD
 REAL, INTENT(IN) :: humidities(nd_humidities, nd_aerosol_species)
                                         ! Humidities for aerpsol species
 LOGICAL, INTENT(IN) :: l_use_arcl(npd_arcl_species)
                                         ! Logical for aerosol clim. species

! logical switches for the calculation of differnt AODs:
 LOGICAL, INTENT(IN) :: l_aod_sulphate 
 LOGICAL, INTENT(IN) :: l_aod_dust
 LOGICAL, INTENT(IN) :: l_aod_seasalt
 LOGICAL, INTENT(IN) :: l_aod_soot
 LOGICAL, INTENT(IN) :: l_aod_biomass
 LOGICAL, INTENT(IN) :: l_aod_biogenic
 LOGICAL, INTENT(IN) :: l_aod_ocff
 LOGICAL, INTENT(IN) :: l_aod_delta
 LOGICAL, INTENT(IN) :: l_aod_nitrate
 LOGICAL, INTENT(IN) :: l_aod_total_radn
 LOGICAL, INTENT(IN) :: l_angst_total_radn
 LOGICAL, INTENT(IN) :: l_aod_prog_sulphate
 LOGICAL, INTENT(IN) :: l_aod_prog_dust
 LOGICAL, INTENT(IN) :: l_aod_prog_seasalt
 LOGICAL, INTENT(IN) :: l_aod_prog_soot
 LOGICAL, INTENT(IN) :: l_aod_prog_biomass
 LOGICAL, INTENT(IN) :: l_aod_prog_ocff
 LOGICAL, INTENT(IN) :: l_aod_prog_nitrate

!     Input data itself:
 INTEGER, INTENT(IN):: i_humidity_pointer(nd_profile, nd_layer)
                                         ! Pointer to look-up table for aerosols
 REAL, INTENT(IN) ::  aerosol_mix_ratio(nd_profile, nd_layer,         &
                        nd_aerosol_mixratio) ! Mixing ratios of aerosols
 REAL, INTENT(IN) ::  d_mass(nd_profile, nd_layer) ! Mass thickness of each layer
 REAL, INTENT(IN) ::  mean_rel_humidity(nd_profile, nd_layer)  
                                             ! Mean relative humidity of layers
 REAL, INTENT(IN) ::  delta_humidity         ! Increment in look-up table for hum.

 REAL, INTENT(IN) ::  aod_wavel(nd_aod_wavel) ! Wavelengths for AOD calculations

! AOD data to be calculated:
 REAL, INTENT(INOUT) ::aod_sulphate(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_dust(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_seasalt(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_soot(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_biomass(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_biogenic(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_ocff(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_delta(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_nitrate(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_total_radn(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::angst_total_radn(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_prog_sulphate(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_prog_dust(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_prog_seasalt(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_prog_soot(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_prog_biomass(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_prog_ocff(nd_profile, nd_aod_wavel)
 REAL, INTENT(INOUT) ::aod_prog_nitrate(nd_profile, nd_aod_wavel)

!  loop indices
 INTEGER k, l

!<Declaration order for all the above; INTEGER, REAL, LOGICAL, CHARACTER>

! Function definitions

! End of header
IF (lhook) CALL dr_hook('compute_all_aod',zhook_in,zhook_handle)


!     For each aerosol type, compute the optical depth if
!     it was requested.
!
IF (l_aod_sulphate) THEN
! DEPENDS ON: compute_aod
  CALL COMPUTE_AOD(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_sulphate, .true., .false.,        &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_sulphate)
END IF
IF (l_aod_dust) THEN
  IF (l_twobin_dust .and. .not.l_use_arcl(ip_arcl_dust) ) THEN
     ! two-bin prognostic dust in radiation
! DEPENDS ON: compute_aod
    CALL compute_aod(                                                   &
      n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio, &
      n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                 &
      n_layer, first_layer, nd_layer,                                   &
      n_humidities, nd_humidities,                                      &
      type_aerosol, i_aod_type, ip_type_twobdust, .true., .false.,      &
      aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,      &
      d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,&
      i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,&
      aod_dust)
  ELSE  ! not two-bin dust
! DEPENDS ON: compute_aod
    CALL compute_aod(                                                   &
      n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio, &
      n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                 &
      n_layer, first_layer, nd_layer,                                   &
      n_humidities, nd_humidities,                                      &
      type_aerosol, i_aod_type, ip_type_dust, .true., .false.,          &
      aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,      &
      d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,&
      i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,&
      aod_dust)
  END IF   ! ends not two-bin dust
END IF
IF (l_aod_seasalt) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_seasalt, .true., .false.,         &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_seasalt)
END IF
IF (l_aod_soot) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_soot, .true., .false.,            &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_soot)
END IF
IF (l_aod_biomass) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_biomass, .true., .false.,         &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_biomass)
END IF
IF (l_aod_biogenic) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_biogenic, .true., .false.,        &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_biogenic)
END IF
IF (l_aod_ocff) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_ocff, .true., .false.,            &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_ocff)
END IF
IF (l_aod_delta) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_delta, .true., .false.,           &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_delta)
END IF
IF (l_aod_nitrate) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_nitrate, .true., .false.,         &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_nitrate)
END IF
IF (l_aod_total_radn .or. l_angst_total_radn) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_allaod, .true., .false.,          &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_total_radn)
  IF (l_angst_total_radn) THEN
    ! calculate the angstrom coefficient for different aod wavelengths
    ! angst=alog(aot_w1/aot_w2)/alog(w2/w1)

    IF (n_aod_wavel > 1 ) THEN
      ! for now, set the aod wavelengths hardcoded here:

      ! first wavelength, non-centred gradient:
      k = 1
      DO l = 1, n_profile
        angst_total_radn(l, k) = log(aod_total_radn(l,k)/aod_total_radn(l,k+1))/&
                                 log(aod_wavel(k+1)/aod_wavel(k))
      END DO
      ! centre wavelengths, centred gradients:
      DO k=2, n_aod_wavel-1
        DO l = 1, n_profile
         angst_total_radn(l,k)=log(aod_total_radn(l,k-1)/aod_total_radn(l,k+1))/&
                               log(aod_wavel(k+1)/aod_wavel(k-1))
        END DO
      END DO
      ! final wavelength, non-centred gradient:
      k = n_aod_wavel
      DO l = 1, n_profile
        angst_total_radn(l, k) = log(aod_total_radn(l,k-1)/aod_total_radn(l,k))/&
                                 log(aod_wavel(k)/aod_wavel(k-1))
      END DO

    ELSE
      DO k = 1, n_aod_wavel
        DO l = 1, n_profile
          angst_total_radn(l, k) = 0.0
        END DO
      END DO
    END IF

  END IF
END IF
IF (l_aod_prog_sulphate) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_sulphate, .false., .true.,        &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_prog_sulphate)
END IF
IF (l_aod_prog_dust) THEN
  IF (l_twobin_dust) THEN
     ! two-bin prognostic dust
! DEPENDS ON: compute_aod
    CALL compute_aod(                                                   &
      n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio, &
      n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                 &
      n_layer, first_layer, nd_layer,                                   &
      n_humidities, nd_humidities,                                      &
      type_aerosol, i_aod_type, ip_type_twobdust, .false., .true.,      &
      aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,      &
      d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,&
      i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,&
      aod_prog_dust)
  ELSE  ! not two-bin dust
! DEPENDS ON: compute_aod
    CALL compute_aod(                                                   &
      n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio, &
      n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                 &
      n_layer, first_layer, nd_layer,                                   &
      n_humidities, nd_humidities,                                      &
      type_aerosol, i_aod_type, ip_type_dust, .false., .true.,          &
      aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,      &
      d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,&
      i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,&
      aod_prog_dust)
  END IF   ! ends not two-bin dust
END IF
IF (l_aod_prog_seasalt) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_seasalt, .false., .true.,         &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_prog_seasalt)
END IF
IF (l_aod_prog_soot) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_soot, .false., .true.,            &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_prog_soot)
END IF
IF (l_aod_prog_biomass) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_biomass, .false., .true.,         &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_prog_biomass)
END IF
IF (l_aod_prog_ocff) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_ocff, .false., .true.,            &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_prog_ocff)
END IF
IF (l_aod_prog_nitrate) THEN
! DEPENDS ON: compute_aod
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_nitrate, .false., .true.,         &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_prog_nitrate)
END IF

IF (lhook) CALL dr_hook('compute_all_aod',zhook_out,zhook_handle)

RETURN
END SUBROUTINE compute_all_aod
