! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description: Contains the horizontal ENDGame grid.
!  
! Method:
!
! Documentation: ENDGame formulation version 1.01
!
! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

MODULE interp_grid_const_mod
IMPLICIT NONE

! Constants derived from the horizontal grid for use in interpolation

! sig and tau are the positions of the points which would be at
! -1 & 2 on a uniform grid.

REAL,SAVE, TARGET, ALLOCATABLE ::                                     &
  sig_xi1_p(:), sig_xi1_u(:),                                         &
  sig_xi2_p(:), sig_xi2_v(:)

REAL,SAVE, TARGET, ALLOCATABLE ::                                     &
  tau_xi1_p(:), tau_xi1_u(:),                                         &
  tau_xi2_p(:), tau_xi2_v(:)

! Divisors for the cubic lagrange interpolants q1,  q2 & q3.
! q4 is for the quintic in the vertical

REAL,SAVE, TARGET, ALLOCATABLE ::                                     &
  q1_p(:,:), q1_u(:,:), q2_p(:,:), q2_v(:,:),                        &
  q3_p(:,:), q3_w(:,:), q5_p(:,:), q5_w(:,:)

END MODULE
