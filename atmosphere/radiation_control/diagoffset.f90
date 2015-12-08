! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:
!   Set the offset off the item numbers for the forcing (and radiance)
!   diagnostics in STASH.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.

MODULE diagoffset

IMPLICIT NONE

! Offsets for diagnostics from successive calls to radiation code
  INTEGER, PARAMETER :: diagnostic_offset = 200

END MODULE diagoffset
