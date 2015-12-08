! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : INGWD

MODULE s_ingwd

  USE scm_cntl_mod, ONLY: scm_nml

  USE scm_utils, ONLY:                                  &
    rmdi, zhook_in, zhook_out                           &
  , jprb, lhook, dr_hook

  USE s_main_force, ONLY:                               &
    sd_orog_land                                        &
  , orog_grad_xx_land                                   &
  , orog_grad_xy_land                                   &
  , orog_grad_yy_land

  IMPLICIT NONE

!=============================================================================
!
! Description:
!   Allows SCM to read in INGWD namelist from forcing file, <scm_nml>
!
! Method:
!   Namelist INGWD is defined in this module and read in by contained
!   subroutine read_ingwd.  Scalar variables are read directly to
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

  NAMELIST/ingwd/                                                             &
    sd_orog_land, orog_grad_xx_land,  orog_grad_xy_land,  orog_grad_yy_land

  PRIVATE :: ingwd
!=============================================================================
CONTAINS

  SUBROUTINE read_ingwd

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    INTEGER :: istatus
    INTEGER :: icode

    CHARACTER(LEN=11), PARAMETER :: routinename='read_ingwd'

    !=========================================================================

    IF (lhook) CALL dr_hook('READ_INGWD',zhook_in,zhook_handle)

    OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus)

    IF (istatus /= 0) THEN
      icode = 500
      WRITE(6,*) ' Error opening file on unit 10 from '//routinename
      WRITE(6,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
      WRITE(6,*) ' IOstat =', istatus
    END IF

    READ  (10, INGWD)
    CLOSE (10)

    IF (lhook) CALL dr_hook('READ_INGWD',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE read_ingwd

!=============================================================================
END MODULE s_ingwd
