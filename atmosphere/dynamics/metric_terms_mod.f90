MODULE metric_terms_mod
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Description: Contains the ENDGame metric terms
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
                         h1_p    (:,:,:),                             &
                         h1_xi1_u(:,:,:),                             &
                         h1_xi2_v(:,:,:),                             &
                         h1_p_eta(:,:,:),                             &
                         h2_p    (:,:,:),                             &
                         h2_xi1_u(:,:,:),                             &
                         h2_xi2_v(:,:,:),                             &
                         h2_p_eta(:,:,:),                             &
                         h3_p    (:,:,:),                             &
                         h3_xi1_u(:,:,:),                             &
                         h3_xi2_v(:,:,:),                             &
                         h3_p_eta(:,:,:),                             &
                         deta_xi3(:,:,:),                             &
                         deta_xi3_theta(:,:,:),                       &
                         deta_xi3_u(:,:,:),                           &
                         deta_xi3_v(:,:,:),                           &
                         dxi1_xi3(:,:,:),                             &
                         dxi2_xi3(:,:,:)

CONTAINS

SUBROUTINE init_metric_terms()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE ereport_mod
USE atm_fields_bounds_mod

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER ierr

IF (lhook) CALL dr_hook('INIT_METRIC_TERMS',zhook_in,zhook_handle)

ALLOCATE (  h1_p  (pdims_s%i_start:pdims_s%i_end,                     &
                   pdims_s%j_start:pdims_s%j_end,                     &
                   pdims_s%k_start:pdims_s%k_end),                    &
            h2_p  (pdims_s%i_start:pdims_s%i_end,                     &
                   pdims_s%j_start:pdims_s%j_end,                     &
                   pdims_s%k_start:pdims_s%k_end),                    &
            h3_p  (pdims_s%i_start:pdims_s%i_end,                     &
                   pdims_s%j_start:pdims_s%j_end,                     &
                   pdims_s%k_start:pdims_s%k_end),                    &
          h1_p_eta(tdims_s%i_start:tdims_s%i_end,                     &
                   tdims_s%j_start:tdims_s%j_end,                     &
                   tdims_s%k_start:tdims_s%k_end),                    &
          h2_p_eta(tdims_s%i_start:tdims_s%i_end,                     &
                   tdims_s%j_start:tdims_s%j_end,                     &
                   tdims_s%k_start:tdims_s%k_end),                    &
          h3_p_eta(tdims_s%i_start:tdims_s%i_end,                     &
                   tdims_s%j_start:tdims_s%j_end,                     &
                   tdims_s%k_start:tdims_s%k_end),                    &
          h1_xi1_u(udims_s%i_start:udims_s%i_end,                     &
                   udims_s%j_start:udims_s%j_end,                     &
                   udims_s%k_start:udims_s%k_end),                    &
          h2_xi1_u(udims_s%i_start:udims_s%i_end,                     &
                   udims_s%j_start:udims_s%j_end,                     &
                   udims_s%k_start:udims_s%k_end),                    &
          h3_xi1_u(udims_s%i_start:udims_s%i_end,                     &
                   udims_s%j_start:udims_s%j_end,                     &
                   udims_s%k_start:udims_s%k_end),                    &
          h1_xi2_v(vdims_s%i_start:vdims_s%i_end,                     &
                   vdims_s%j_start:vdims_s%j_end,                     &
                   vdims_s%k_start:vdims_s%k_end),                    &
          h2_xi2_v(vdims_s%i_start:vdims_s%i_end,                     &
                   vdims_s%j_start:vdims_s%j_end,                     &
                   vdims_s%k_start:vdims_s%k_end),                    &
          h3_xi2_v(vdims_s%i_start:vdims_s%i_end,                     &
                   vdims_s%j_start:vdims_s%j_end,                     &
                   vdims_s%k_start:vdims_s%k_end),                    &
          deta_xi3(pdims_s%i_start:pdims_s%i_end,                     &
                   pdims_s%j_start:pdims_s%j_end,                     &
                   pdims_s%k_start:pdims_s%k_end),                    &
    deta_xi3_theta(tdims_s%i_start:tdims_s%i_end,                     &
                   tdims_s%j_start:tdims_s%j_end,                     &
                   tdims_s%k_start:tdims_s%k_end),                    &
          dxi1_xi3(tdims_s%i_start:tdims_s%i_end,                     &
                   tdims_s%j_start:tdims_s%j_end,                     &
                   tdims_s%k_start:tdims_s%k_end),                    &
          dxi2_xi3(tdims_s%i_start:tdims_s%i_end,                     &
                   tdims_s%j_start:tdims_s%j_end,                     &
                   tdims_s%k_start:tdims_s%k_end),                    &
        deta_xi3_u(udims_s%i_start:udims_s%i_end,                     &
                   udims_s%j_start:udims_s%j_end,                     &
                   udims_s%k_start:udims_s%k_end),                    &
        deta_xi3_v(vdims_s%i_start:vdims_s%i_end,                     &
                   vdims_s%j_start:vdims_s%j_end,                     &
                   vdims_s%k_start:vdims_s%k_end),                    &
                STAT=ierr)

IF (ierr/=0) CALL Ereport("init_metric_terms",ierr, "Unable to allocate.")

IF (lhook) CALL dr_hook('INIT_METRIC_TERMS',zhook_out,zhook_handle)

END SUBROUTINE init_metric_terms
END MODULE metric_terms_mod
