MODULE horiz_grid_mod
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
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

IMPLICIT NONE

! Horizontal coordinates


REAL,SAVE, ALLOCATABLE ::                                             &
  xi1_p(:),                                                           &
  xi1_u(:),                                                           &
  xi2_p(:),                                                           &
  xi2_v(:)

REAL,SAVE, ALLOCATABLE ::                                             &
  glob_xi1_p(:),                                                      &
  glob_xi1_u(:),                                                      &
  glob_xi2_p(:),                                                      &
  glob_xi2_v(:)

REAL,SAVE, ALLOCATABLE ::                                             &
  glob_dxi1_p(:),                                                     &
  glob_dxi1_u(:),                                                     &
  glob_dxi2_p(:),                                                     &
  glob_dxi2_v(:)

REAL,SAVE, ALLOCATABLE ::                                             &
  glob_rdxi1_p(:),                                                    &
  glob_rdxi1_u(:),                                                    &
  glob_rdxi2_p(:),                                                    &
  glob_rdxi2_v(:)


! Horizontal & Vertical coordinates

REAL,SAVE, ALLOCATABLE ::                                             &
  csxi1_p(:),                                                         &
  csxi1_u(:),                                                         &
  csxi2_p(:),                                                         &
  csxi2_v(:)

REAL,SAVE, ALLOCATABLE ::                                             &
  snxi1_p(:),                                                         &
  snxi1_u(:),                                                         &
  snxi2_p(:),                                                         &
  snxi2_v(:)

REAL, SAVE, ALLOCATABLE ::                                            &
  phi_at_p(:,:,:),                                                    &
  phi_at_u(:,:,:),                                                    &
  phi_at_v(:,:,:),                                                    &
  phi_at_eta(:,:,:)

!
! Linear Interpolation weights
!

REAL, SAVE, ALLOCATABLE ::                                            &
  intw_u2p(:,:),                                                      &
  intw_v2p(:,:),                                                      & 
  intw_p2u(:,:),                                                      &
  intw_p2v(:,:),                                                      &
  intw_rho2w(:,:),                                                    &
  intw_w2rho(:,:)

REAL  base_xi1,base_xi2,delta_xi1, delta_xi2

! Variable resolution parameters
INTEGER :: Nxi1L, Nxi1V, Nxi2L, Nxi2V
REAL    :: delta_xi1_H, delta_xi1_L, delta_xi2_H, delta_xi2_L

END MODULE
