! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Climate Diagnostics
MODULE interpor_mod

IMPLICIT NONE
! ----------------------- include file: INTERPOR -----------------------
! Description: Hold magic numbers for order of vertical interpolation,
!              where prime use is for atmosphere physics/dynamics
!              diagnostics.
!              interp_order should be set to a _default value where
!              used generally in the code, but can be set explicitly
!              locally where required.
!
      INTEGER,PARAMETER:: interp_order_linear       = 1
      INTEGER,PARAMETER:: interp_order_linear_noex  = 2
      INTEGER,PARAMETER:: interp_order_cubic        = 3
      INTEGER,PARAMETER:: interp_order_quintic      = 5

! INTERPOR end

END MODULE interpor_mod
