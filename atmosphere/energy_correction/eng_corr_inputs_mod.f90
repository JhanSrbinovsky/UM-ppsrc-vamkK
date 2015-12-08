! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
! Description:
!   Module containing runtime options/data used by the energy correction scheme
!
! Method:
!   Switches and associated data values used by the energy correction scheme
!   are defined here and assigned default values. These may be overridden
!   by namelist input.
!
!   Any routine wishing to use these options may do so with the 'USE'
!   statement.
!
! Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Energy Correction

MODULE eng_corr_inputs_mod

USE missing_data_mod, ONLY: imdi

IMPLICIT NONE

!===========================================================================
! LOGICAL options set from run_eng_corr namelist
!===========================================================================

LOGICAL :: l_emcorr    = .FALSE.  ! Turns on energy correction code

LOGICAL :: lmass_corr  = .FALSE.  ! Apply mass correction

LOGICAL :: lemq_print  = .FALSE.  ! Print additional info from em code

!===========================================================================
! INTEGER values set from run_eng_corr namelist
!===========================================================================

INTEGER :: a_energysteps = imdi   ! Number of model timesteps per energy 
                                  ! correction period.

!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------

! Define the run_eng_corr namelist

NAMELIST /run_eng_corr/ l_emcorr, lmass_corr, lemq_print, a_energysteps

!===========================================================================
! LOGICAL options not set in namelist
!===========================================================================

LOGICAL :: lqt_corr    = .FALSE.  ! Apply total moisture correction

LOGICAL :: lenergy     = .FALSE.  ! Timestep to cal energy correction

LOGICAL :: lflux_reset = .FALSE.  ! Timestep to reset flux array in D1

END MODULE eng_corr_inputs_mod
