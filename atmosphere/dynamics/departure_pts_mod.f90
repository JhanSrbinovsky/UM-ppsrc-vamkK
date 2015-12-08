MODULE departure_pts_mod
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


REAL,SAVE, ALLOCATABLE ::                                             &
  depart_xi1_u  (:,:,:), depart_xi2_u  (:,:,:), depart_xi3_u  (:,:,:),&
  depart_xi1_v  (:,:,:), depart_xi2_v  (:,:,:), depart_xi3_v  (:,:,:),&
  depart_xi1_w  (:,:,:), depart_xi2_w  (:,:,:), depart_xi3_w  (:,:,:),&
  depart_xi1_rho(:,:,:), depart_xi2_rho(:,:,:), depart_xi3_rho(:,:,:)

REAL,SAVE, ALLOCATABLE ::                                             &
  depart_xi1_w_disp(:,:,:),  depart_xi2_w_disp(:,:,:)


REAL, ALLOCATABLE ::                                                  &
     depart_xi1_rho_halo(:,:,:),                                      &
     depart_xi2_rho_halo(:,:,:),                                      &
     depart_xi3_rho_halo(:,:,:)


CONTAINS

SUBROUTINE init_depart_pts()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod
USE atm_fields_bounds_mod

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER ierr

IF (lhook) CALL dr_hook('INIT_DEPART_PTS',zhook_in,zhook_handle)

ALLOCATE (depart_xi1_u(udims%i_start:udims%i_end,                     &
                       udims%j_start:udims%j_end,                     &
                       udims%k_start:udims%k_end),                    &
          depart_xi2_u(udims%i_start:udims%i_end,                     &
                       udims%j_start:udims%j_end,                     &
                       udims%k_start:udims%k_end),                    &
          depart_xi3_u(udims%i_start:udims%i_end,                     &
                       udims%j_start:udims%j_end,                     &
                       udims%k_start:udims%k_end),                    &
          depart_xi1_v(vdims%i_start:vdims%i_end,                     &
                       vdims%j_start:vdims%j_end,                     &
                       vdims%k_start:vdims%k_end),                    &
          depart_xi2_v(vdims%i_start:vdims%i_end,                     &
                       vdims%j_start:vdims%j_end,                     &
                       vdims%k_start:vdims%k_end),                    &
          depart_xi3_v(vdims%i_start:vdims%i_end,                     &
                       vdims%j_start:vdims%j_end,                     &
                       vdims%k_start:vdims%k_end),                    &
          depart_xi1_w(wdims%i_start:wdims%i_end,                     &
                       wdims%j_start:wdims%j_end,                     &
                       wdims%k_start:wdims%k_end),                    &
          depart_xi2_w(wdims%i_start:wdims%i_end,                     &
                       wdims%j_start:wdims%j_end,                     &
                       wdims%k_start:wdims%k_end),                    &
          depart_xi3_w(wdims%i_start:wdims%i_end,                     &
                       wdims%j_start:wdims%j_end,                     &
                       wdims%k_start:wdims%k_end),                    &
        depart_xi1_rho(pdims%i_start:pdims%i_end,                     &
                       pdims%j_start:pdims%j_end,                     &
                       pdims%k_start:pdims%k_end),                    &
        depart_xi2_rho(pdims%i_start:pdims%i_end,                     &
                       pdims%j_start:pdims%j_end,                     &
                       pdims%k_start:pdims%k_end),                    &
        depart_xi3_rho(pdims%i_start:pdims%i_end,                     &
                       pdims%j_start:pdims%j_end,                     &
                       pdims%k_start:pdims%k_end),                    &
        STAT=ierr)

IF (ierr/=0) CALL Ereport("init_depart_pts",ierr, "Unable to allocate.")

ALLOCATE (depart_xi1_w_disp(wdims  %i_start:wdims  %i_end,             &
                            wdims  %j_start:wdims  %j_end,             &
                            wdims  %k_start:wdims  %k_end),            &
          depart_xi2_w_disp(wdims  %i_start:wdims  %i_end,             &
                            wdims  %j_start:wdims  %j_end,             &
                            wdims  %k_start:wdims  %k_end),STAT=ierr)

IF (ierr/=0) CALL Ereport("init_depart_pts",ierr, "Unable to allocate.")

ALLOCATE(                                                             &
         depart_xi1_rho_halo(pdims_s%i_start:pdims_s%i_end,           &
                             pdims_s%j_start:pdims_s%j_end,           &
                             pdims  %k_start:pdims  %k_end))

ALLOCATE(                                                             &
         depart_xi2_rho_halo(pdims_s%i_start:pdims_s%i_end,           &
                             pdims_s%j_start:pdims_s%j_end,           &
                             pdims  %k_start:pdims  %k_end))

ALLOCATE(                                                             &
         depart_xi3_rho_halo(pdims_s%i_start:pdims_s%i_end,           &
                             pdims_s%j_start:pdims_s%j_end,           &
                             pdims  %k_start:pdims  %k_end))

IF (lhook) CALL dr_hook('INIT_DEPART_PTS',zhook_out,zhook_handle)

END SUBROUTINE init_depart_pts

SUBROUTINE reset_dpt_pts()

USE horiz_grid_mod
USE atm_fields_bounds_mod, ONLY : udims,vdims,pdims
USE level_heights_mod,     ONLY : eta_rho_levels,eta_theta_levels
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER i,j,k
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('RESET_DPT_PTS',zhook_in,zhook_handle)

  DO k = 1, pdims%k_end
     DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
           depart_xi1_u(i,j,k) = xi1_u(i)
           depart_xi2_u(i,j,k) = xi2_p(j)
           depart_xi3_u(i,j,k) = eta_rho_levels(k)
        END DO
     END DO
  END DO
  DO k = 1, pdims%k_end
     DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
             depart_xi1_v(i,j,k) = xi1_p(i)
             depart_xi2_v(i,j,k) = xi2_v(j)
             depart_xi3_v(i,j,k) = eta_rho_levels(k)
        END DO
     END DO
  END DO
  DO k = 0, pdims%k_end
     DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
           depart_xi1_w(i,j,k) = xi1_p(i)
           depart_xi2_w(i,j,k) = xi2_p(j)
           depart_xi3_w(i,j,k) = eta_theta_levels(k)
        END DO
     END DO
  END DO
  DO k = 1, pdims%k_end
     DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
           depart_xi1_rho(i,j,k) = xi1_p(i)
           depart_xi2_rho(i,j,k) = xi2_p(j)
           depart_xi3_rho(i,j,k) = eta_rho_levels(k)
        END DO
     END DO
  END DO

  IF (lhook) CALL dr_hook('RESET_DPT_PTS',zhook_out,zhook_handle)

END SUBROUTINE


END MODULE departure_pts_mod
