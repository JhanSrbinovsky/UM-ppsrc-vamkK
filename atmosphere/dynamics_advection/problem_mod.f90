! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
MODULE problem_mod

IMPLICIT NONE
! Description: COMDECK containing problem_number
!  for use in setting problem types
!
      INTEGER, PARAMETER :: standard=0
      INTEGER, PARAMETER :: monsoon=1
      INTEGER, PARAMETER :: dynamical_core=2
      INTEGER, PARAMETER :: idealised_problem=3
      INTEGER, PARAMETER :: standard_namelist=4

END MODULE problem_mod
