! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  reg_grid_mod

      MODULE reg_grid_mod
      IMPLICIT NONE

! Description:  regular grid constants

! Method:   These are calculated in set_var_grid via setcona. 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

      REAL  :: base_lambda
      REAL  :: base_phi
      REAL  :: delta_lambda
      REAL  :: delta_phi
      REAL  :: recip_delta_lambda
      REAL  :: recip_delta_phi
      REAL  :: f_plane_rad
      REAL  :: ff_plane_rad
      REAL  :: lat_rot_NP_deg
      REAL  :: long_rot_NP_deg
      REAL  :: lat_rot_NP
      REAL  :: long_rot_NP
 
      END MODULE reg_grid_mod
