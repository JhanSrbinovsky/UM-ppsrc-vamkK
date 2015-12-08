! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Module for allocation of Filename length.

MODULE filenamelength_mod

IMPLICIT NONE

! Description:
!   The module sets the max length of a filename in the UM.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

INTEGER, PARAMETER      :: filenamelength = 256
INTEGER, PARAMETER      :: datawnamelength = 256
INTEGER, PARAMETER      :: runidnamelength = 20

!-------------------------------------------------------------------

END MODULE filenamelength_mod
