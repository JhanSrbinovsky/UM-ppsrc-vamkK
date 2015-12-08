! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!+ data module for switches/options concerned with the land surface.
  ! Description:
  !   Module containing runtime options/data used by the land surface.
  !
  ! Method:
  !   Switches and associated data values used by the land-surface scheme
  !   are defined here and assigned default values. These may be overridden
  !   by namelist input.
  !
  !   Any routine wishing to use these options may do so with the 'USE'
  !   statement.
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Boundary Layer
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !

MODULE land_surf_mod

  ! Declarations:

  IMPLICIT NONE

!=======================================================================
! The module contains switches that exist in the JULES modules switches
! in the section #if !defined(UM_RUN) so are not available to the UM
!=======================================================================

!=======================================================================
! RUN_BLVEG namelists items
!=======================================================================

INTEGER :: ilayers=10    ! For radiation model

NAMELIST/RUN_BLVEG/ilayers

END MODULE land_surf_mod
