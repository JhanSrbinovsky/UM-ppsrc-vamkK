! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!=====================================================================
! FUNCTION DayNew
! PURPOSE:-           To calculate SIN of argument (in eqn. 12
!                     in SCM doc.) required in calculation of
!                     mean or SD of variable at day relative to winter
!                     solstice
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
!=====================================================================

FUNCTION daynew (at, bt, itd)

!---------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER :: &
    itd       ! In Dayno. relative to winter solstice

  REAL ::    &
    at       &! In Constants for calculating annual cycle
  , bt

!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------
  REAL ::    &
    arg      &! Argument
  , daynew    ! SIN of argument

  arg = at * float(itd) + bt
  daynew = SIN(arg)

  RETURN

END FUNCTION daynew

!=====================================================================

