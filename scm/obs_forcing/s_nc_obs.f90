! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : NC_OBS

MODULE s_nc_obs

  USE scm_cntl_mod, ONLY: scm_nml

  USE scm_utils, ONLY:                                                        &
      zhook_in, zhook_out, jprb, lhook, dr_hook

  USE s_main_force, ONLY: netcdf_file, obs_t0_prf, source

  IMPLICIT NONE

!=============================================================================
!
! Description:
!   Allows SCM to read in NC_OBS namelist from forcing file, scm_nml.
!   NC_OBS contains information relating to the use of NetCDF driver files to
!   force the SCM
!
! Method:
!   Namelist NC_OBS is defined in this module and read in by contained
!   subroutine read_nc_obs.  Scalar variables are read directly to
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

  !---------------------------------------------------------------------------
  ! Define namelist
  !---------------------------------------------------------------------------
  NAMELIST/NC_OBS/                                                            &
   netcdf_file, obs_t0_prf, source

  PRIVATE :: nc_obs

!=============================================================================
CONTAINS

   SUBROUTINE read_nc_obs

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    INTEGER :: istatus
    INTEGER :: icode

    CHARACTER(LEN=11), PARAMETER :: routinename='read_nc_obs'

    IF (lhook) CALL dr_hook('READ_NC_OBS',zhook_in,zhook_handle)

    ! Set defaults
    netcdf_file = ''
    obs_t0_prf  = 1

    OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus)

    IF (istatus /= 0) THEN
      icode = 500
      WRITE(6,*) ' Error opening file on unit 10 from '//routinename
      WRITE(6,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
      WRITE(6,*) ' IOstat =', istatus
    END IF

    READ  (10, nc_obs)
    CLOSE (10)

    IF (lhook) CALL dr_hook('READ_NC_OBS',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE read_nc_obs

!=============================================================================
END MODULE s_nc_obs
