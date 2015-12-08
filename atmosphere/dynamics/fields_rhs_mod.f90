MODULE fields_rhs_mod
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

IMPLICIT NONE

! Horizontal coordinates


REAL,SAVE, ALLOCATABLE ::                                             &
   r_u      (:,:,:),&
   r_v      (:,:,:),&
   r_w      (:,:,:),&
   r_u_p2   (:,:,:),&
   r_v_p2   (:,:,:),&
   r_w_p2   (:,:,:),&
   r_u_p2_n (:,:,:),&
   r_v_p2_n (:,:,:),&
   r_w_p2_n (:,:,:),&
   r_theta  (:,:,:),&
   r_rho    (:,:,:),&
   r_m_v    (:,:,:),&
   r_m_cl   (:,:,:),&
   r_m_cf   (:,:,:),&
   r_m_r    (:,:,:),&
   r_m_gr   (:,:,:),&
   r_m_cf2  (:,:,:),&
   r_u_d    (:,:,:),&
   r_v_d    (:,:,:),&
   r_w_d    (:,:,:),&
   r_theta_d(:,:,:),&
   r_rho_d  (:,:,:),&
   r_m_v_d  (:,:,:),&
   r_m_cl_d (:,:,:),&
   r_m_cf_d (:,:,:),&
   r_m_r_d  (:,:,:),&
   r_m_gr_d (:,:,:),&
   r_m_cf2_d(:,:,:),&
   r_p_p    (:,:,:),&
   r_p_p_d  (:,:,:),&
   s_u      (:,:,:),&
   s_v      (:,:,:),&
   s_w      (:,:,:),&
   s_thetav (:,:,:),&
   s_m_v    (:,:,:),&
   s_m_cl   (:,:,:),&
   s_m_cf   (:,:,:),&
   s_m_r    (:,:,:),&
   s_m_gr   (:,:,:),&
   s_m_cf2  (:,:,:),&
   r_u_skeb (:,:,:),&
   r_v_skeb (:,:,:)

CONTAINS

SUBROUTINE init_r_skeb()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod
USE atm_fields_bounds_mod

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER ierr

IF (lhook) CALL dr_hook('FIELDS_RHS_MOD:INIT_R_SKEB',zhook_in,zhook_handle)

ALLOCATE (    r_u_skeb(udims_s%i_start:udims_s%i_end,                     &
                       udims_s%j_start:udims_s%j_end,                     &
                       udims_s%k_start:udims_s%k_end),STAT=ierr)

IF (ierr/=0) CALL Ereport("init_r_skeb",ierr, "Unable to allocate1.")

ALLOCATE (    r_v_skeb(vdims_s%i_start:vdims_s%i_end,                     &
                       vdims_s%j_start:vdims_s%j_end,                     &
                       vdims_s%k_start:vdims_s%k_end),STAT=ierr)

IF (ierr/=0) CALL Ereport("init_r_skeb",ierr, "Unable to allocate2.")

r_u_skeb = 0.0
r_v_skeb = 0.0

IF (lhook) CALL dr_hook('FIELDS_RHS_MOD:INIT_R_SKEB',zhook_out,zhook_handle)

END SUBROUTINE init_r_skeb

SUBROUTINE destroy_r_skeb()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE ereport_mod
USE atm_fields_bounds_mod

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('FIELDS_RHS_MOD:DESTROY_R_SKEB',zhook_in,zhook_handle)

DEALLOCATE ( r_v_skeb )
DEALLOCATE ( r_u_skeb )

IF (lhook) CALL dr_hook('FIELDS_RHS_MOD:DESTROY_R_SKEB',zhook_out,zhook_handle)

END SUBROUTINE destroy_r_skeb


SUBROUTINE init_fields_rhs(l_skeb2)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE ereport_mod
USE atm_fields_bounds_mod

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

LOGICAL, INTENT (IN) :: l_skeb2

INTEGER ierr

IF (lhook) CALL dr_hook('FIELDS_RHS_MOD:INIT_FIELDS_RHS',zhook_in,zhook_handle)

ALLOCATE (         r_u(udims_s%i_start:udims_s%i_end,                     &
                       udims_s%j_start:udims_s%j_end,                     &
                       udims_s%k_start:udims_s%k_end),                    &
                 r_u_d(udims_s%i_start:udims_s%i_end,                     &
                       udims_s%j_start:udims_s%j_end,                     &
                       udims_s%k_start:udims_s%k_end),                    &
                r_u_p2(udims_s%i_start:udims_s%i_end,                     &
                       udims_s%j_start:udims_s%j_end,                     &
                       udims_s%k_start:udims_s%k_end),                    &
              r_u_p2_n(udims_s%i_start:udims_s%i_end,                     &
                       udims_s%j_start:udims_s%j_end,                     &
                       udims_s%k_start:udims_s%k_end),                    &
                   s_u(udims_s%i_start:udims_s%i_end,                     &
                       udims_s%j_start:udims_s%j_end,                     &
                       udims_s%k_start:udims_s%k_end),                    &
                   r_v(vdims_s%i_start:vdims_s%i_end,                     &
                       vdims_s%j_start:vdims_s%j_end,                     &
                       vdims_s%k_start:vdims_s%k_end),                    &
                 r_v_d(vdims_s%i_start:vdims_s%i_end,                     &
                       vdims_s%j_start:vdims_s%j_end,                     &
                       vdims_s%k_start:vdims_s%k_end),                    &
                r_v_p2(vdims_s%i_start:vdims_s%i_end,                     &
                       vdims_s%j_start:vdims_s%j_end,                     &
                       vdims_s%k_start:vdims_s%k_end),                    &
              r_v_p2_n(vdims_s%i_start:vdims_s%i_end,                     &
                       vdims_s%j_start:vdims_s%j_end,                     &
                       vdims_s%k_start:vdims_s%k_end),                    &
                   s_v(vdims_s%i_start:vdims_s%i_end,                     &
                       vdims_s%j_start:vdims_s%j_end,                     &
                       vdims_s%k_start:vdims_s%k_end),                    &
                   r_w(wdims%i_start:wdims%i_end,                         &
                       wdims%j_start:wdims%j_end,                         &
                       wdims%k_start:wdims%k_end),                        &
                 r_w_d(wdims%i_start:wdims%i_end,                         &
                       wdims%j_start:wdims%j_end,                         &
                       wdims%k_start:wdims%k_end),                        &
                r_w_p2(wdims%i_start:wdims%i_end,                         &
                       wdims%j_start:wdims%j_end,                         &
                       wdims%k_start:wdims%k_end),                        &
              r_w_p2_n(wdims%i_start:wdims%i_end,                         &
                       wdims%j_start:wdims%j_end,                         &
                       wdims%k_start:wdims%k_end),                        &
                   s_w(wdims%i_start:wdims%i_end,                         &
                       wdims%j_start:wdims%j_end,                         &
                       wdims%k_start:wdims%k_end),                        &
               r_theta(tdims_s%i_start:tdims_s%i_end,                     &
                       tdims_s%j_start:tdims_s%j_end,                     &
                       tdims_s%k_start:tdims_s%k_end),                    &
             r_theta_d(tdims_s%i_start:tdims_s%i_end,                     &
                       tdims_s%j_start:tdims_s%j_end,                     &
                       tdims_s%k_start:tdims_s%k_end),                    &
              s_thetav(tdims_s%i_start:tdims_s%i_end,                     &
                       tdims_s%j_start:tdims_s%j_end,                     &
                       tdims_s%k_start:tdims_s%k_end),                    &
                 r_rho(pdims_s%i_start:pdims_s%i_end,                     &
                       pdims_s%j_start:pdims_s%j_end,                     &
                       pdims_s%k_start:pdims_s%k_end),                    &
               r_rho_d(pdims_s%i_start:pdims_s%i_end,                     &
                       pdims_s%j_start:pdims_s%j_end,                     &
                       pdims_s%k_start:pdims_s%k_end),                    &
                 r_p_p(pdims_s%i_start:pdims_s%i_end,                     &
                       pdims_s%j_start:pdims_s%j_end,                     &
                       pdims_s%k_start:pdims_s%k_end),                    &
               r_p_p_d(pdims_s%i_start:pdims_s%i_end,                     &
                       pdims_s%j_start:pdims_s%j_end,                     &
                       pdims_s%k_start:pdims_s%k_end),                    &
                 r_m_v(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
                r_m_cl(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
                r_m_cf(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
               r_m_cf2(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
                r_m_gr(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
                 r_m_r(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
               r_m_v_d(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
              r_m_cl_d(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
              r_m_cf_d(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
             r_m_cf2_d(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
              r_m_gr_d(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
               r_m_r_d(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
                 s_m_v(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
                s_m_cl(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
                s_m_cf(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
               s_m_cf2(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
                s_m_gr(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),                    &
                 s_m_r(qdims_s%i_start:qdims_s%i_end,                     &
                       qdims_s%j_start:qdims_s%j_end,                     &
                       qdims_s%k_start:qdims_s%k_end),STAT=ierr)

IF (ierr/=0) CALL Ereport("init_fields_rhs",ierr, "Unable to allocate.")

IF (l_skeb2) CALL init_r_skeb()

! These fields may not be initialised with certain physics switches 
! but are then used
! in eg_helm_rhs_star (for example) without the flags being tested! 
! Therefore we do initialise them here to zero:

s_m_r  (:,:,:) = 0.0
s_m_gr (:,:,:) = 0.0
s_m_cf2(:,:,:) = 0.0
s_m_cf (:,:,:) = 0.0
s_m_cl (:,:,:) = 0.0

IF (lhook) CALL dr_hook('FIELDS_RHS_MOD:INIT_FIELDS_RHS',zhook_out,zhook_handle)

END SUBROUTINE init_fields_rhs
END MODULE fields_rhs_mod
