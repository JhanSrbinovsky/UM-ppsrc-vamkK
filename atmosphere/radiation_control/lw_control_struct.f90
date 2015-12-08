! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!   This module defines the controlling structure for LW calculations.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90

MODULE lw_control_struct

  USE max_calls
  USE control_struc
  USE missing_data_mod  ! Note, rmdi needed for glue_rad

  IMPLICIT NONE

  TYPE (control_option), SAVE :: lw_control(npd_lwcall)

  INTEGER, SAVE :: n_lwcall = imdi
!   The number of LW calls to the radiation

END MODULE lw_control_struct
