! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE grdtypes_mod

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

IMPLICIT NONE
! PARAMETERs defining the GRID_TYPE characteristics

! General
      INTEGER, PARAMETER :: gt_unset = -1
      ! Any value which is unset

! MODEL_TYPE
      INTEGER, PARAMETER :: gt_atmos = 1
      ! Atmosphere field

      INTEGER, PARAMETER :: gt_ocean = 2
      ! Ocean field

      INTEGER, PARAMETER :: gt_wave = 3
      ! Wave field

! CONTENT
      INTEGER, PARAMETER :: gt_thetamass = 1
      ! Contains theta or mass points

      INTEGER, PARAMETER :: gt_velocity = 2
      ! Contains velocity (B grid U,V) points

      INTEGER, PARAMETER :: gt_U_C = 3
      ! Contains U points on C grid

      INTEGER, PARAMETER :: gt_V_C = 4
      ! Contains V points on C grid

      INTEGER, PARAMETER :: gt_hybrid = 5
      ! Points on none of the above
      INTEGER, PARAMETER :: gt_river = 6
      ! River routing grid

! COVERAGE
      INTEGER, PARAMETER :: gt_allpts = 1
      ! All points

      INTEGER, PARAMETER :: gt_land = 2
      ! Land points

      INTEGER, PARAMETER :: gt_sea = 3
      ! Sea points

! DOMAIN
      INTEGER, PARAMETER :: gt_full = 1
      ! Full field

      INTEGER, PARAMETER :: gt_zonal = 2
      ! Zonal field

      INTEGER, PARAMETER :: gt_meridional = 3
      ! Meridional field

      INTEGER, PARAMETER :: gt_ozone = 4
      ! Ozone field

      INTEGER, PARAMETER :: gt_scalar = 5
      ! Single point

      INTEGER, PARAMETER :: gt_compressed = 6
      ! Compressed points

      INTEGER, PARAMETER :: gt_LBC = 7
      ! Lateral Boundary Condition Field

! CYCLIC

      INTEGER, PARAMETER :: gt_nocyclic = 1
      ! No cyclic columns

      INTEGER, PARAMETER :: gt_optcyclic = 2
      ! Optional cyclic columns

      INTEGER, PARAMETER :: gt_cyclic = 3
      ! Includes cyclic columns

END MODULE grdtypes_mod
