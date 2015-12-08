! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
MODULE tprofile_mod

IMPLICIT NONE
! Description: COMDECK containing temperature profile types
!  for use in idealised problems
!
      INTEGER, PARAMETER :: tp_dthetadz=1
      INTEGER, PARAMETER :: tp_isothermal=2
      INTEGER, PARAMETER :: tp_BruntV=3
      INTEGER, PARAMETER :: tp_BV_isoth=4
      INTEGER, PARAMETER :: tp_dyn_core=6
      INTEGER, PARAMETER :: tp_dyn_core_lam=7
      INTEGER, PARAMETER :: tp_namelist=9
      INTEGER, PARAMETER :: tp_dump=10

END MODULE tprofile_mod
