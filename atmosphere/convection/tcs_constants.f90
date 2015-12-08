! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! 
! Module holding collection of physical constants and derived 
! combinations thereof.
!
MODULE tcs_constants

USE earth_constants_mod, ONLY: g, earth_radius

USE atmos_constants_mod, ONLY:             &
    r, cp, repsilon, c_virtual, kappa, pref, rv,                       &
    recip_epsilon, recip_kappa

USE water_constants_mod, ONLY: lc, lf

USE cv_diag_param_mod, ONLY:                                           &
   a_bolton, b_bolton, c_bolton, d_bolton

  IMPLICIT NONE
  !
  ! Description:
  !   This module holds a collection of physical constants and derived 
  ! combinations thereof.
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.1 programming standards.
  !

  ! Model constants:

  !--------------------------------------
  ! Combined Physical constants
  !--------------------------------------
  REAL, PARAMETER ::                                          &
     lc_o_cp     = lc / cp                                    &
     , gamma_dry = g/cp                                       &
     , ra2       = 1. / (earth_radius*earth_radius)           &
     , cv        = cp - R

END MODULE tcs_constants
