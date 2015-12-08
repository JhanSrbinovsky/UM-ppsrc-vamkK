! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
! Description:
!   Module containing runtime options/data used by the river routing scheme
!
! Method:
!   Switches and associated data values used by the river routing scheme
!   are defined here and assigned default values. These may be overridden
!   by namelist input.
!
!   Any routine wishing to use these options may do so with the 'USE'
!   statement.
!
! Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: River Routing

MODULE river_inputs_mod

USE missing_data_mod, ONLY: rmdi,imdi

IMPLICIT NONE

!===========================================================================
! REAL values set from RUN_RIVERS namelist
!===========================================================================

REAL :: river_vel   = rmdi      ! River velocity

REAL :: river_mcoef = rmdi      ! Meandering coefficient

REAL :: river_step  = rmdi      ! River routing timestep

!===========================================================================
! LOGICAL options set from RUN_RIVERS namelist
!===========================================================================

LOGICAL :: l_rivers = .FALSE.   ! Use global river routing scheme 
  
LOGICAL :: l_inland = .FALSE.   ! Control rerouting of inland basin water


!===========================================================================
! INTEGER options set from RUN_RIVERS namelist
!===========================================================================
INTEGER :: i_river_vn = imdi

!----------------------------------------------------------------------

! Define the run_rivers namelist

NAMELIST /RUN_RIVERS/ river_vel, river_mcoef, river_step, l_rivers, &
                      l_inland, i_river_vn

END MODULE river_inputs_mod
