! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : LOGIC

MODULE s_logic

  USE s_maxdim,  ONLY:                                                        &
    mx_rw_lng, mx_rw

  USE scm_cntl_mod,  ONLY: scm_nml

  USE scm_utils, ONLY:                                                        &
    imdi, rw_lng, rw, zhook_in, zhook_out, jprb, lhook, dr_hook

  USE s_main_force, ONLY:                                                     &
    dev_test, test, altdat, geoforce, geoinit, noforce, obs, obs_surf, stats  &
  , local_time, ancyc, l_flux_bc, l_spec_z0, l_geo_centred, l_qpos_for        &
  , grafdump_day, grafdump_days, grafdump_step, prindump_day                  &
  , prindump_days, prindump_obs, prindump_step, prinstat                      &
  , radcloud_fixed, tapein, tapeout                                           &
  , scm_land_sea_mask => land_sea_mask &
  , scm_land_ice_mask => land_ice_mask &
  , scm_soil_mask     => soil_mask     &
  , scm_cumulus       => cumulus        

   IMPLICIT NONE

!=============================================================================
!
! Description:
!   Allows SCM to read in LOGIC namelist from forcing file, <scm_nml>
!
! Method:
!   Namelist LOGIC is defined in this module and read in by contained
!   subroutine read_logic.  Scalar variables are read directly to
!   s_main_force.  Array variables are declared as private single rank arrays 
!   of fixed length (arrays cannot be varying length in namelists) and
!   reshaped before being transferred to arrays of the correct size/shape in
!   s_main_force.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!=============================================================================
!
! Declarations:

  LOGICAL, PRIVATE ::                      &
    land_sea_mask (mx_rw_lng*mx_rw) = .false. &! land point flag (any type)
  , land_ice_mask (mx_rw_lng*mx_rw) = .false. &! land point flag (ice)
  , soil_mask     (mx_rw_lng*mx_rw) = .false. &! land point flag (soil)
  , cumulus       (mx_rw_lng*mx_rw) = .false.  ! bl convection flag


  NAMELIST/logic/                                                             &
    dev_test, test, altdat, geoforce, geoinit, noforce, obs, obs_surf, stats  &
  , local_time, ancyc, l_flux_bc, l_spec_z0, l_geo_centred, l_qpos_for        &
  , grafdump_day, grafdump_days, grafdump_step, prindump_obs, prinstat        &
  , prindump_day, prindump_days, prindump_step, radcloud_fixed, tapein        &
  , tapeout, land_sea_mask, land_ice_mask, soil_mask, cumulus

  PRIVATE :: logic

!=============================================================================
CONTAINS

  SUBROUTINE read_logic

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    INTEGER :: istatus
    INTEGER :: icode

    CHARACTER(LEN=11), PARAMETER :: routinename='read_logic'

    !=========================================================================

    IF (lhook) CALL dr_hook('READ_LOGIC',zhook_in,zhook_handle)

    OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus)

    IF (istatus /= 0) THEN
      icode = 500
      WRITE(6,*) ' Error opening file on unit 10 from '//routinename
      WRITE(6,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
      WRITE(6,*) ' IOstat =', istatus
    END IF

    READ  (10, LOGIC)
    CLOSE (10)


    scm_land_sea_mask = RESHAPE( land_sea_mask, (/rw_lng,rw/) )

    scm_land_ice_mask = RESHAPE( land_ice_mask, (/rw_lng,rw/) )

    scm_soil_mask     = RESHAPE( soil_mask,     (/rw_lng,rw/) )

    scm_cumulus       = RESHAPE( cumulus,       (/rw_lng,rw/) )

    IF (lhook) CALL dr_hook('READ_LOGIC',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE read_logic

!=============================================================================
END MODULE s_logic
