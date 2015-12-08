! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
! Description:
!   Module containing runtime options/data used by the murk scheme
!
! Method:
!   Switches and associated data values used by the murk scheme
!   are defined here and assigned default values. These may be overridden
!   by namelist input.
!
!   Any routine wishing to use these options may do so with the 'USE'
!   statement.
!
! Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Atmosphere Service

MODULE murk_inputs_mod

IMPLICIT NONE

!===========================================================================
! LOGICAL options set from run_murk namelist
!===========================================================================

LOGICAL :: l_murk        = .FALSE.  ! Use total aerosol field

LOGICAL :: l_murk_advect = .FALSE.  ! Aerosol advection

LOGICAL :: l_murk_source = .FALSE.  ! Aerosol source & sink terms

LOGICAL :: l_murk_lbc    = .FALSE.  ! Murk aerosol lbcs active

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------

! Define the run_murk namelist

NAMELIST /run_murk/ l_murk, l_murk_advect, l_murk_source, l_murk_lbc

!===========================================================================
! LOGICAL options not set in namelist
!===========================================================================

LOGICAL :: l_murk_bdry   = .FALSE.  ! UK Mes boundary model

LOGICAL :: l_murk_rad    = .FALSE.  ! Include radiative effects of aerosol

END MODULE murk_inputs_mod
