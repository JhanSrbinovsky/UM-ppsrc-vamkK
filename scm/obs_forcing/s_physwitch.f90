! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : PHYSWITCH

MODULE s_physwitch

  USE scm_cntl_mod, ONLY: scm_nml

  USE scm_utils,  ONLY:                                                       &
      zhook_in, zhook_out, jprb, lhook, dr_hook

  USE s_main_force, ONLY: conv_mode

  IMPLICIT NONE

!=============================================================================
!
! Description:
!   Allows SCM to read in PHYSWITCH namelist from forcing file, scm_nml.
!   PHYSWITCH contains information on the different physics modes.
!
! Method:
!   Namelist PHYSWITCH is defined in this module and read in by contained
!   subroutine read_physwitch.  Scalar variables are read directly to
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

  NAMELIST/physwitch/ conv_mode

  PRIVATE :: physwitch

!=============================================================================
CONTAINS

  SUBROUTINE read_physwitch

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    INTEGER :: istatus
    INTEGER :: icode

    CHARACTER(LEN=11), PARAMETER :: routinename='read_physwitch'

    IF (lhook) CALL dr_hook('READ_PHYSWITCH',zhook_in,zhook_handle)
 
    OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus)

    IF (istatus /= 0) THEN
      icode = 500
      WRITE(6,*) ' Error opening file on unit 10 from '//routinename
      WRITE(6,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
      WRITE(6,*) ' IOstat =', istatus
    END IF

    READ  (10, physwitch)
    CLOSE (10)

    IF (lhook) CALL dr_hook('READ_PHYSWITCH',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE read_physwitch

!=============================================================================
END MODULE s_physwitch
