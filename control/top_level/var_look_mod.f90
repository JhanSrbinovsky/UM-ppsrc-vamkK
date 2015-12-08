! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  var_look_mod

      MODULE var_look_mod

      IMPLICIT NONE

! Description:  End values for variable resolution grid
!
! Method:   These are calculated in set_var_grid via setcona. 
!           
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

! Set of look-up tables for searching on variable grids
      INTEGER,  ALLOCATABLE :: look_lam(:)     ! lambda p search
      INTEGER,  ALLOCATABLE :: look_phi(:)     ! phi p search

      INTEGER :: halo_lam              ! halo for lamp look-up table
      INTEGER :: halo_phi              ! halo for phip look-up table
      INTEGER :: max_look              ! max size of look-up table

      REAL ::  recip_dlam                 ! smallest delta_lambda p
      REAL ::  recip_dphi                 ! smallest delta_phi p

      END MODULE var_look_mod
