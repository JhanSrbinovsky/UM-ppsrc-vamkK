! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ List of UKMO grib related magic numbers

MODULE rcf_ukmocodes_mod

! Description:
!    Magic numbers for grib codes, as defined by
!    the UKMO versions of 'Table 2' 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: Fortran 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

! UKMO codes used in the reconfiguration

! GEMS related fields without codes in table 210
! Now stored in UKMO table no 141
! GEMS NO 
INTEGER, PARAMETER      :: ukmocode_no                  = 1

! GEMS HNO3 
INTEGER, PARAMETER      :: ukmocode_hno3                = 2

! GEMS PAN 
INTEGER, PARAMETER      :: ukmocode_pan                 = 3

! GEMS C2H6 
INTEGER, PARAMETER      :: ukmocode_c2h6                = 4

! GEMS C3H8 
INTEGER, PARAMETER      :: ukmocode_c3h8                = 5

! Dust in UM bins converted from GEMS 3 bin scheme
INTEGER, PARAMETER      :: ukmocode_umdust1             = 6
INTEGER, PARAMETER      :: ukmocode_umdust2             = 7
INTEGER, PARAMETER      :: ukmocode_umdust3             = 8
INTEGER, PARAMETER      :: ukmocode_umdust4             = 9
INTEGER, PARAMETER      :: ukmocode_umdust5             = 10
INTEGER, PARAMETER      :: ukmocode_umdust6             = 11

! Sulfate - UM uses 2 modes (plus dissolved). converted from single mode in GEMS
INTEGER, PARAMETER      :: ukmocode_umso4aitk           = 12
INTEGER, PARAMETER      :: ukmocode_umso4accu           = 13

! SO2 - note that in mmr(S) not mmr
INTEGER, PARAMETER      :: ukmocode_so2                 = 14


END MODULE rcf_ukmocodes_mod
