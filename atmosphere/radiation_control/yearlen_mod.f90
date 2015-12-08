! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
MODULE yearlen_mod

IMPLICIT NONE
!
!
! ----------------------- Header file YEARLEN  -------------------------
! Description: Parameter of the length of the tropical year
!----------------------------------------------------------------------
!     ! Number of days in the tropical year
!     ! The tropical year is defined as the mean interval between two
!     ! successive passages of the sun through the vernal equinox.
!     ! This value is 365.2424 and will be used here.
!     ! The value 365.2422 most often used is the mean length of
!     ! the year starting at different points on the ellipse.
      Real, Parameter    :: TropYearLength = 365.2424

END MODULE yearlen_mod
