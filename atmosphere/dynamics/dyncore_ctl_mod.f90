! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for dynamics.

! Description:
!   Module containing input ctl switches/settings as used by the
!   dynamical core
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: top_level

! Method:
!   Switches are initialised to false and read in from the
!   UMUI. The module may then be used directly where the switches
!   are needed within the dynamics code.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE dyncore_ctl_mod

IMPLICIT NONE

! Suarez-Held variables
REAL, ALLOCATABLE :: SuHe_level_weight(:)
REAL, ALLOCATABLE :: friction_level(:)
REAL, ALLOCATABLE :: frictional_timescale(:)

! ----------------------

END MODULE  dyncore_ctl_mod
