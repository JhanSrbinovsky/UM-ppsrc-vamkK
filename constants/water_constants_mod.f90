! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Water related physical constants

MODULE water_constants_mod

! Description:
!        Water related physical constants

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.

IMPLICIT NONE

! tfs, temperature at which sea water freezes
REAL, PARAMETER :: tfs                 = 271.35
! tm, temperature at which fresh water freezes and ice melts
REAL, PARAMETER :: tm                  = 273.15

! density of pure water (kg/m3)
REAL, PARAMETER :: rho_water           = 1000.0
! density of sea water (kg/m3)
REAL, PARAMETER :: rhosea              = 1026.0

! latent heat of condensation of water at 0degc
REAL, PARAMETER :: lc                  = 2.501E6
! latent heat of fusion of water at 0degc
REAL, PARAMETER :: lf                  = 0.334E6

END MODULE water_constants_mod
