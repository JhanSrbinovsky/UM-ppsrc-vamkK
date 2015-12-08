! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
MODULE gflat_mod

IMPLICIT NONE
! Description: COMDECK containing vertical coordinate flattening types
!  for use in idealised problems

      INTEGER, PARAMETER :: gflat_old=1
      INTEGER, PARAMETER :: gflat_linear=0
      INTEGER, PARAMETER :: gflat_linear1=2
      INTEGER, PARAMETER :: gflat_quadratic=3

END MODULE gflat_mod
