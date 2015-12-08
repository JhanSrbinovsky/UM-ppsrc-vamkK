! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module: RIV_CONCERNS -------------------------------------------------
!
!  Purpose: Location of persistent data structures needed by 
!  the river routing routines to coordinate parallel regridding 
!  between atmosphere and river grids.
!  As mentioned data structures below are persistent until the end 
!  of a UM run (i.e. allocated but not deallocated, any future
!  modifications should take this into consideration)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: River Routing


MODULE RIV_CONCERNS

  USE regrid_types

  IMPLICIT NONE 
  
  TYPE (CONCERN), POINTER :: runoff_send_concern(:) => NULL()
  ! information for other process's need to regrid 
  ! runoff to trips grid

  TYPE (CONCERN), POINTER :: runoff_recv_concern(:) => NULL()
  ! information this process needs to regrid 
  ! runoff to trips grid
  
  TYPE (CONCERN_MAX), POINTER :: riverout_send_concern(:) => NULL()
  ! information for other process's need to regrid 
  ! river outflow to atmos (pressure) grid

  TYPE (CONCERN_MAX), POINTER :: riverout_recv_concern(:) => NULL()
  ! information this process needs to regrid 
  ! river outflow to atmos (pressure grid)

  TYPE (CONTRIBUTION_INFO), POINTER :: riverout_contribution(:) => NULL()
  ! convenience data structure for storing regridding information for 
  ! a single target grid point for river outflow regridding


  TYPE (CONCERN_MAX), POINTER :: outflow_send_concern(:) => NULL()
  ! information for other process's need to regrid 
  ! outflow to atmos (pressure) grid 
  
  TYPE (CONCERN_MAX), POINTER :: outflow_recv_concern(:) => NULL()
  ! information this process needs to regrid 
  ! outflow to atmos (pressure) grid 

  TYPE (CONTRIBUTION_INFO), POINTER :: outflow_contribution(:) => NULL()
  ! convenience data structure for storing regridding information for 
  ! a single target grid point for river outflow regridding

END MODULE RIV_CONCERNS
