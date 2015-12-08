! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************
! Description: Surface types settings

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer

! Code Description:
!   Language: FORTRAN 90

MODULE ukca_nstypes

IMPLICIT NONE

! **********************************************************************
! This is a copy of nstypes_MOSES and should not be proliferated. Use
! the module nstypes instead, which in JULES has these defined as
! variables
! **********************************************************************

! Number of non-vegetation surface types
INTEGER, PARAMETER :: nnvg  = 4

! Number of plant functional types.
INTEGER, PARAMETER :: npft  = 5

! Number of surface types.
INTEGER, PARAMETER :: ntype = 9

! Index of the surface type 'Soil'
INTEGER, PARAMETER :: soil  = 8

! Index of the surface type 'Lake'
INTEGER, PARAMETER :: lake  = 7

! Land surface types :
!     1 - Broadleaf Tree
!     2 - Needleleaf Tree
!     3 - C3 Grass
!     4 - C4 Grass
!     5 - Shrub
!     6 - Urban
!     7 - Water
!     8 - Soil
!     9 - Ice

END MODULE ukca_nstypes
