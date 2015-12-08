! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
MODULE c_idl_force_options_mod

IMPLICIT NONE
! Start C_IDL_FORCE_OPTIONS

! Description: Include file containing idealised forcing options
!
      INTEGER, PARAMETER :: no_forcing      = 0
      INTEGER, PARAMETER :: force_increment = 1
      INTEGER, PARAMETER :: force_relax     = 2
      INTEGER, PARAMETER :: force_reset     = 3

! End C_IDL_FORCE_OPTIONS

END MODULE c_idl_force_options_mod
