! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP.
!
! parameters enumerating different model types (essentially captring
! their boundary conditions in a convenient form factor)

MODULE domain_params
  IMPLICIT NONE

  ! Small execs don't use this, so we have a placeholder value
  INTEGER,PARAMETER :: mt_smexe         = 0

  ! Global models (i.e. cyclic EW, and NS overpole - a sphere)
  INTEGER,PARAMETER :: mt_global        = 1

  ! No boundary conditions (i.e. the world has an edge)
  INTEGER,PARAMETER :: mt_lam           = 2

  ! models that wrap EW
  INTEGER,PARAMETER :: mt_cyclic_lam    = 3

  ! models that wrap EW and NS (i.e. torus like)
  INTEGER,PARAMETER :: mt_bi_cyclic_lam = 4

  ! Single column (i.e. one horizontal point)
  INTEGER,PARAMETER :: mt_single_column = 5

END MODULE domain_params
