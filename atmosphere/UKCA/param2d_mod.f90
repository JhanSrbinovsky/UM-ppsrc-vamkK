! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
MODULE param2d_mod

IMPLICIT NONE
!
!     parameter param2d
!     Parameters for 2-D model/photolysis
!
!     Parameter    Description
!     ---------    -----------
!     nolat        Number of 2-D model latitudes
!     nolev        Number of 2-D model levels
!     nlphot       Number of 2-D model sub-levels used in the calculation
!                  of photolysis rates (approx. 1km spacing).
!     ntphot       Number of times per day that data from the 2-D model
!                  photolysis scheme is stored (time 3 is noon).
!     jpphio       Fortran i/o unit to write out and read in photolysis
!                  rates
!
!     ----------------------------------------------------------------
!
      INTEGER, PARAMETER :: nolat  = 19
      INTEGER, PARAMETER :: nolev  = 17
      INTEGER, PARAMETER :: nlphot = 51      
      INTEGER, PARAMETER :: ntphot = 3
      INTEGER, PARAMETER :: jpphio = 84
      INTEGER, PARAMETER :: jpphin = 58      

END MODULE param2d_mod
