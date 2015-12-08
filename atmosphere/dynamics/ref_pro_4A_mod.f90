! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


MODULE ref_pro_mod
! Description: Contains the horizontal ENDGame grid.
!  
!
! Method:
!
! Documentation: ENDGame formulation version 1.01
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.



! Horizontal coordinates

IMPLICIT NONE

REAL, SAVE  , ALLOCATABLE ::                                          &
        exner_ref_pro(:,:,:)                                          & 
      , thetav_ref_pro(:,:,:)                                         &
      , rho_ref_pro(:,:,:)                                            &
      , thetav_ref_eta(:,:,:)                                         &
      , rho_ref_eta(:,:,:)



LOGICAL, SAVE ::  reference_profile_changed      = .FALSE.



CONTAINS

SUBROUTINE init_ref_pro()

USE atm_fields_bounds_mod, ONLY: pdims_s,tdims_s
USE ereport_mod, ONLY: ereport
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

INTEGER ierr
CHARACTER(LEN=*) RoutineName
PARAMETER (   RoutineName='ref_pro_4A')
CHARACTER(LEN=256)       :: message
INTEGER                  :: errorstatus
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('ref_pro_4A',zhook_in,zhook_handle)

ALLOCATE                                                              &
       (exner_ref_pro(   pdims_s%i_start:pdims_s%i_end,               &
                         pdims_s%j_start:pdims_s%j_end,               &
                         pdims_s%k_start-1:pdims_s%k_end+1)           &
      , thetav_ref_pro(  tdims_s%i_start:tdims_s%i_end,               &
                         tdims_s%j_start:tdims_s%j_end,               &
                         tdims_s%k_start:tdims_s%k_end)               &
      , rho_ref_pro(     pdims_s%i_start:pdims_s%i_end,               &
                         pdims_s%j_start:pdims_s%j_end,               &
                         pdims_s%k_start:pdims_s%k_end)               &
      , thetav_ref_eta(  pdims_s%i_start:pdims_s%i_end,               &
                         pdims_s%j_start:pdims_s%j_end,               &
                         pdims_s%k_start:pdims_s%k_end)               &
      , rho_ref_eta(     tdims_s%i_start:tdims_s%i_end,               &
                         tdims_s%j_start:tdims_s%j_end,               &
                         tdims_s%k_start:tdims_s%k_end),STAT=ierr)

IF(ierr.ne.0) THEN
  errorstatus = abs(ierr)
  message = "Allocation failed"
  CALL ereport(RoutineName, errorstatus, message)
ENDIF

IF (lhook) CALL dr_hook('ref_pro_4A',zhook_out,zhook_handle)

END SUBROUTINE


END MODULE
