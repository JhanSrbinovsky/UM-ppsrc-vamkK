! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Horizontal grid for variable resolution

MODULE vrhoriz_grid_mod

USE Atmos_Max_Sizes
USE UM_ParParams
IMPLICIT NONE
!
! Description:
!   Contains the horizontal grid namelist for variable resolution, shared
!   between MakeBC and the reconfiguration
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
REAL :: lambda_input_p(row_length_max)      ! longitude points on p-grid
REAL :: lambda_input_u(row_length_max)      ! longitude points on u-grid
REAL :: phi_input_p(rows_max)               ! latitude points on p-grid
REAL :: phi_input_v(rows_max)               ! latitude points on v-grid

NAMELIST /HORIZGRID/                                                       &
  lambda_input_p, lambda_input_u,                                          &
  phi_input_p, phi_input_v

END MODULE vrhoriz_grid_mod
