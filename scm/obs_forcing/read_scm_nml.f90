! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads in forcing file 'namelist.scm' and passes data back to scm_main

SUBROUTINE read_scm_nml

  USE scm_utils, ONLY:                                                        &
    nlnd, zhook_in, zhook_out, jprb, lhook, dr_hook

  ! Module holding forcing file data for run
  USE s_main_force,  ONLY:  geoforce, radcloud_fixed

  ! Modules to read in SCM sub-namelists from forcing file
  USE s_logic,     ONLY: read_logic
  USE s_inprof,    ONLY: read_inprof
  USE s_indata,    ONLY: read_indata
  USE s_rundata,   ONLY: read_rundata
  USE s_inobsfor,  ONLY: read_inobsfor
  USE s_injules,   ONLY: read_injules
  USE s_ingwd,     ONLY: read_ingwd
  USE s_ingeofor,  ONLY: read_ingeofor
  USE s_radcloud,  ONLY: read_radcloud
  USE s_physwitch, ONLY: read_physwitch

  IMPLICIT NONE

!-----------------------------------------------------------------------------
!
! Description:
!   Calls rountines to reads in SCM sub-namelists from scm forcing file
!   (default:namelist.scm)
!
! Method:
!   Fortran namelists cannot be declared with non-constant dimensions. So
!   namelist variables are read using their respective subroutines. The
!   reading subroutines read in the namelists, resize the arrays and
!   interpolate vertical profiles as required if the forcing namelist is
!   a different resolution to that of the SCM run.
!
!   After calling this routine, all namelist variables can be accessed via
!   s_main_force module. 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!
!-----------------------------------------------------------------------------
! Declarations:

  ! Dr Hook
  !==============================
  REAL(KIND=jprb) :: zhook_handle


  IF (lhook) CALL dr_hook('READ_SCM_NML',zhook_in,zhook_handle)

  CALL read_logic

  CALL read_inprof

  CALL read_indata

  CALL read_rundata

  CALL read_inobsfor

  CALL read_physwitch

  IF (nlnd > 0) THEN
    CALL read_injules
    CALL read_ingwd
  END IF

  IF (geoforce) THEN 
    CALL read_ingeofor
  END IF

  IF (radcloud_fixed) THEN 
    CALL read_radcloud
  END IF

  IF (lhook) CALL dr_hook('READ_SCM_NML',zhook_out,zhook_handle)


  RETURN
END SUBROUTINE read_scm_nml
