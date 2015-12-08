! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!   This module declares the controlling structure for SW radiation.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90

MODULE sw_control_struct

  USE max_calls
  USE control_struc
  USE missing_data_mod  ! Note, rmdi needed for glue_rad

  IMPLICIT NONE

  TYPE (control_option), SAVE :: sw_control(npd_swcall)

  INTEGER, SAVE :: n_swcall = imdi
!   The number of SW calls to the radiation

END MODULE sw_control_struct
