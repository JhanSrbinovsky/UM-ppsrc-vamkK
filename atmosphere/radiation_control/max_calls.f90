! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module to set the maximum number of calls to the radiation.
!
! Description:
!   This module declares the maximum number of calls to the
!   radiation code permitted on a single timestep.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
!- End of header

MODULE max_calls

IMPLICIT NONE

!   Size allocated for arrays concerned with the SW call
  INTEGER, PARAMETER :: npd_swcall=5

!   Size allocated for arrays concerned with the LW call
  INTEGER, PARAMETER :: npd_lwcall=5

END MODULE max_calls
