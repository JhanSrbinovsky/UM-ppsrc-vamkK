! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!   Module containing JRA pressure levels

!  Part of the Nudged model (see nudging_main.F90)

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
MODULE nudging_jra_plevel_def

IMPLICIT NONE

INTEGER,PARAMETER :: jra_nplevs = 23

! JRA fixed pressure levels
! Note last level is 0: set to 0.0000001 to stop code crashing
REAL,PARAMETER :: jra_pressure(1:jra_nplevs) = (/                      &
 100000, 92500, 85000, 70000, 60000,                                   &
  50000, 40000, 30000, 25000, 20000,                                   &
  15000, 10000,  7000,  5000,  3000,                                   &
   2000,  1000,   700,   500,   300,                                   &
    200,   100,     1                                                  &
 /)

END MODULE nudging_jra_plevel_def

