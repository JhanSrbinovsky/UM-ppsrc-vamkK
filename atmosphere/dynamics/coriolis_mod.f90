MODULE coriolis_mod
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Description: Contains the ENDGame coriolis coefficients
!  
!
! Method:
!
! Documentation: ENDGame formulation version 1.01
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Atmosphere Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

IMPLICIT NONE

! Horizontal coordinates


REAL,SAVE, ALLOCATABLE ::                                             &
! Variable Coriolis terms
        f1_star (:,:,:),                                              &
        f2_star (:,:,:),                                              &
        f3_star (:,:,:),                                              &
! Fixed Coriolis terms
        f1_comp (:,:),                                                &
        f2_comp (:,:),                                                &
        f3_comp (:,:)

CONTAINS

SUBROUTINE init_coriolis()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod
USE atm_fields_bounds_mod

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER ierr

IF (lhook) CALL dr_hook('INIT_CORIOLIS',zhook_in,zhook_handle)

ALLOCATE (    f1_star (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  & 
                       pdims_s%k_start:pdims_s%k_end),                 &
              f2_star (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  & 
                       pdims_s%k_start:pdims_s%k_end),                 &
              f3_star (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  & 
                       pdims_s%k_start:pdims_s%k_end),                 &
              f1_comp (pdims%i_start:pdims%i_end,                      &
                       pdims%j_start:pdims%j_end),                     &
              f2_comp (pdims%i_start:pdims%i_end,                      &
                       pdims%j_start:pdims%j_end),                     &
              f3_comp (pdims%i_start:pdims%i_end,                      &
                       pdims%j_start:pdims%j_end),                     &
              STAT=ierr)

IF (ierr/=0) CALL Ereport("init_coriolis",ierr, "Unable to allocate.")

! initialisation needed here because we do not swapbound this field
! in the global model
f1_star = 0.0

IF (lhook) CALL dr_hook('INIT_CORIOLIS',zhook_out,zhook_handle)

END SUBROUTINE init_coriolis
END MODULE coriolis_mod
