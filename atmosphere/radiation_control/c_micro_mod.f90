! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
MODULE c_micro_mod

IMPLICIT NONE
!
! Description:
!
!  Contains various cloud droplet parameters, defined for
!  land and sea areas.
!
!  NTOT_* is the total number concentration (m-3) of cloud droplets;
!  KPARAM_* is the ratio of the cubes of the volume-mean radius
!                                           and the effective radius;
!  DCONRE_* is the effective radius (m) for deep convective clouds;
!  DEEP_CONVECTION_LIMIT_* is the threshold depth (m) bewteen shallow
!                                          and deep convective cloud.
!
! Current Code Owner: Andy Jones
!
!
!     REAL,PARAMETER:: NTOT_LAND is set in UMUI
!     REAL,PARAMETER:: NTOT_SEA is set in UMUI
      REAL,PARAMETER:: KPARAM_LAND = 0.67
      REAL,PARAMETER:: KPARAM_SEA = 0.80
      REAL,PARAMETER:: DCONRE_LAND = 9.5E-06
      REAL,PARAMETER:: DCONRE_SEA = 16.0E-06
      REAL,PARAMETER:: DEEP_CONVECTION_LIMIT_LAND = 500.0
      REAL,PARAMETER:: DEEP_CONVECTION_LIMIT_SEA = 1500.0
!
! Maximum Temp for homogenous nucleation (deg C)
! Thomo has moved back into mphys_ice_mod
! To call this in your code use the following syntax:
! USE mphys_ice_mod, ONLY: thomo

END MODULE c_micro_mod
