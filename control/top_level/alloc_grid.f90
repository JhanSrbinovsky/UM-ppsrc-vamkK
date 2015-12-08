! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE alloc_grid_mod
IMPLICIT NONE
CONTAINS
  SUBROUTINE alloc_grid()

  USE parkind1, ONLY: jpim, jprb       !DrHook
  USE yomhook,  ONLY: lhook, dr_hook   !DrHook
  USE horiz_grid_mod
  USE atm_fields_bounds_mod

  USE um_parvars,    ONLY : halo_i, halo_j
  USE proc_info_mod, ONLY: datastart=>l_datastart,model_domain,         &
                           global_row_length,global_rows

  IMPLICIT NONE
!
! Description: Allocate arrays common to the stretched and uniform
!              grids.
!  
! Method: 
!
! Documentation: ENDGame formulation version 1.01
!  
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

  INTEGER                       :: ErrorStatus, alloc_stat

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('ALLOC_GRID',zhook_in,zhook_handle)

  ALLOCATE( xi1_p(pdims_l%i_start:pdims_l%i_end),                       &
            xi1_u(udims_l%i_start:udims_l%i_end),                       &
            xi2_p(pdims_l%j_start:pdims_l%j_end),                       &
            xi2_v(vdims_l%j_start:vdims_l%j_end))

  ALLOCATE( glob_xi1_p(1-halo_i:global_row_length+halo_i),              &
            glob_xi1_u(-halo_i:global_row_length+halo_i),               &
            glob_xi2_p(1-halo_j:global_rows+halo_j),                    &
            glob_xi2_v(-halo_j:global_rows+halo_j) )

  ALLOCATE( glob_dxi1_p(1-halo_i:global_row_length+halo_i),             &
            glob_dxi1_u(-halo_i:global_row_length+halo_i),              &
            glob_dxi2_p(1-halo_j:global_rows+halo_j),                   &
            glob_dxi2_v(-halo_j:global_rows+halo_j) )

  ALLOCATE( glob_rdxi1_p(1-halo_i:global_row_length+halo_i),            &
            glob_rdxi1_u(-halo_i:global_row_length+halo_i),             &
            glob_rdxi2_p(1-halo_j:global_rows+halo_j),                  &
            glob_rdxi2_v(-halo_j:global_rows+halo_j) )

  ALLOCATE( csxi1_p(pdims_l%i_start:pdims_l%i_end),                     &
            csxi1_u(udims_l%i_start:udims_l%i_end),                     &
            csxi2_p(pdims_l%j_start:pdims_l%j_end),                     &
            csxi2_v(vdims_l%j_start:vdims_l%j_end),                     &
            snxi1_p(pdims_l%i_start:pdims_l%i_end),                     &
            snxi1_u(udims_l%i_start:udims_l%i_end),                     &
            snxi2_p(pdims_l%j_start:pdims_l%j_end),                     &
            snxi2_v(vdims_l%j_start:vdims_l%j_end),                     &
            phi_at_p(pdims_s%i_start:pdims_s%i_end,                     &
                     pdims_s%j_start:pdims_s%j_end,                     &
                     pdims_s%k_start:pdims_s%k_end),                    &
            phi_at_u(udims_s%i_start:udims_s%i_end,                     &
                     udims_s%j_start:udims_s%j_end,                     &
                     udims_s%k_start:udims_s%k_end),                    &
            phi_at_v(vdims_s%i_start:vdims_s%i_end,                     &
                     vdims_s%j_start:vdims_s%j_end,                     &
                     vdims_s%k_start:vdims_s%k_end),                    &
            phi_at_eta(tdims_s%i_start:tdims_s%i_end,                   &
                       tdims_s%j_start:tdims_s%j_end,                   &
                       tdims_s%k_start:tdims_s%k_end) )

  ALLOCATE( intw_u2p(pdims%i_start:pdims%i_end,2),                      &
            intw_v2p(pdims%j_start:pdims%j_end,2),                      &
            intw_p2u(udims%i_start:pdims%i_end,2),                      &
            intw_p2v(vdims%j_start:vdims%j_end+1,2),                    & 
            intw_rho2w(wdims%k_end,2),                                  &
            intw_w2rho(wdims%k_end,2) )

  IF (lhook) CALL dr_hook('ALLOC_GRID',zhook_out,zhook_handle)

  END SUBROUTINE alloc_grid
END MODULE alloc_grid_mod
