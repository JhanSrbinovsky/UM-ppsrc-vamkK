! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! 
! Module holding commonly derived constants used in convection
!
MODULE cv_derived_constants_mod

USE earth_constants_mod, ONLY: g, earth_radius

USE atmos_constants_mod, ONLY:             &
    r, cp

USE water_constants_mod, ONLY: lc, lf


IMPLICIT NONE
SAVE

!-------------------------------------------------------------------
! Description:
! This module calculate commonly used constants derived from 
! UM constants
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!-------------------------------------------------------------------

REAL, PARAMETER ::           &
  ls = lc+lf                 & ! Latent heat of sublimation
 ,lsrcp = ls/cp              &
 ,lcrcp = lc/cp              &
 ,lfrcp = lf/cp              &
 ,gamma_dry = g/cp           & ! dry adiabatic lapse rate
 ,cv        = cp - R           ! specific heat of dry air at constant volume

REAL, PARAMETER ::                                  &
  ra2 = 1.0 / (earth_radius*earth_radius) 


END MODULE cv_derived_constants_mod
