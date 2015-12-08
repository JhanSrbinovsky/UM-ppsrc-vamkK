! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Module to define identifiers for surface tiles.

MODULE tileid3z

IMPLICIT NONE

! Description:
!   This module defines identifiers for different surface types
!   as used in the radiation scheme.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90

!- End of header



  INTEGER, PARAMETER :: npd_tile_type  = 3
!   Number of types of tile for which space is allocated
  INTEGER, PARAMETER :: ip_ocean_tile  = 1
!   Identifier for open sea
  INTEGER, PARAMETER :: ip_seaice_tile = 2
!   Idenitifer for ice
  INTEGER, PARAMETER :: ip_land_tile   = 3
!   Identifer for land



END MODULE tileid3z
