MODULE helmholtz_const_matrix_mod
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Description: Contains the ENDGame departure points
!  
!
! Method:
!
! Documentation: ENDGame formulation version 1.01
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Atmosphere Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

USE um_types, ONLY  : real64,real32
IMPLICIT NONE

! Horizontal coordinates



INTEGER, PARAMETER :: real_eg_hlm_kind=real32


REAL (KIND=real_eg_hlm_kind) , SAVE, ALLOCATABLE ::                    &
         Hlm_Lp(:,:,:),                                                &
         Hlm_Ln(:,:,:),                                                &
         Hlm_Ls(:,:,:),                                                &
         Hlm_Le(:,:,:),                                                &
         Hlm_Lw(:,:,:),                                                &
         Hlm_Lu(:,:,:),                                                &
         Hlm_Ld(:,:,:),                                                &
         Hlm_Ek(:,:,:),                                                &
         Hlm_Fk(:,:,:),                                                &
         Hlm_Ck(:,:,:),                                                &
         Hu_k(:,:,:), Hd_k(:,:,:)

CONTAINS

SUBROUTINE init_helmholtz_const_matrix()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE ereport_mod
USE atm_fields_bounds_mod

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER ierr

IF (lhook) CALL dr_hook('INIT_HELMHOLTZ_CONST_MATRIX',zhook_in,zhook_handle)

ALLOCATE (   Hlm_Lp   (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  & 
                       pdims_s%k_start:pdims_s%k_end),                 &
             Hlm_Ln   (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  & 
                       pdims_s%k_start:pdims_s%k_end),                 &
             Hlm_Ls   (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  & 
                       pdims_s%k_start:pdims_s%k_end),                 &
             Hlm_Le   (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  & 
                       pdims_s%k_start:pdims_s%k_end),                 &
             Hlm_Lw   (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  & 
                       pdims_s%k_start:pdims_s%k_end),                 &
             Hlm_Lu   (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  & 
                       pdims_s%k_start:pdims_s%k_end),                 &
             Hlm_Ld   (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  & 
                       pdims_s%k_start:pdims_s%k_end),                 &
             Hlm_Ek   (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  & 
                       pdims_s%k_start:pdims_s%k_end),                 &
             Hlm_Fk   (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  & 
                       pdims_s%k_start:pdims_s%k_end),                 &
             Hlm_Ck   (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  & 
                       0:pdims_s%k_end),                               &
             Hu_k     (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  &
                       pdims_s%k_start:pdims_s%k_end),                 &
             Hd_k     (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  &
                       pdims_s%k_start:pdims_s%k_end),                 &
                STAT=ierr)

IF (ierr/=0) CALL Ereport("init_helmholtz_const_matrix",ierr,        &
                            "Unable to allocate.")

IF (lhook) CALL dr_hook('INIT_HELMHOLTZ_CONST_MATRIX',zhook_out,zhook_handle)

END SUBROUTINE init_helmholtz_const_matrix
END MODULE helmholtz_const_matrix_mod
