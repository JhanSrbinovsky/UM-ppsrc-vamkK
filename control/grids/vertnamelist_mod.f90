! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Vertical grid information

MODULE vertnamelist_mod

USE Atmos_Max_Sizes, ONLY : model_levels_max

USE Control_Max_Sizes
IMPLICIT NONE
!
! Description:
!   Contains the vertical grid namelist for variable resolution, shared
!   between LBC generation, MakeBC, the SCM and the reconfiguration
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
! This includes model_levels_max - note that cmaxsize provides model_levels_max
! to the reconfiguration, and a large number of contents of cmaxsize to the 
! SCM via this module

INTEGER :: first_constant_r_rho_level     
 ! Lowest level where rho level has constant radius (i.e. not terrain-following)

REAL    :: z_top_of_model                 ! Top of model (metres)
REAL    :: eta_theta (model_levels_max+1) ! Theta levels (as fraction of whole)
REAL    :: eta_rho (model_levels_max)     ! Rho levels (as fraction of whole)

NAMELIST /VERTLEVS/                                                       &
  first_constant_r_rho_level,                                             &
  z_top_of_model, eta_theta, eta_rho

END MODULE vertnamelist_mod
