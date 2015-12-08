! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Global bypass module for switches concerned with microphysics

MODULE mphys_bypass_mod

  USE domain_params
  IMPLICIT NONE
  SAVE

!  Description:
!   Module containing logical switches used by the microphysics code
!   that are not part of the run_precip or run_cloud namelists

!  Method:
!   Switches are initialised to false and then set in atm_step to
!   the corresponding switch used in cntlatm.h and read in from the
!   UMUI. The module may then be used directly where the switches
!   are needed within the microphysics code.

!   This module is only for switches not in the run_cloud or
!   run_precip namelists. The switches below can usually be found in
!   the nlstcatm namelist, which is currently not in a module, but 
!   remains as an include file. New microphysics or cloud scheme 
!   logicals should be put in the run_precip or run_cloud namelist 
!   and not here. 

!   Once someone places the cntlatm namelsist into a module, or it is
!   made redundant and the switches placed in physics-specific modules
!   then this module can be removed and the code altered to point to the
!   new variable locations.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large-scale preciptation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

!-----------------------------------------------------
! 1. Microphysics Switches not part of run_precip
!-----------------------------------------------------

LOGICAL :: l_crystals     = .FALSE.

!-----------------------------------------------------
! 2. Dimensions and top of model; 
! required in microphysics
!-----------------------------------------------------

REAL :: mphys_mod_top
REAL :: mp_dell
REAL :: mp_delp

END MODULE mphys_bypass_mod
