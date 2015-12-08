! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Module containing constants relating to ECMWF fixed & hybrid
!  pressure levels for old (60 level) version.
!  Based on the routine Rcf_ecmwf_60level_def

!  Method:
!   The levels used in the ECMWF model are defined by two numbers for
!   each level k, a(k) and b(k). The height of the level above the
!   surface at a grid point j in terms of pressure is given by

!         P(j,k) = a(k) + b(k) * Pstar(j)

!   where all pressure measurements are in Pa (not hPa). The following
!   parameters define the half levels used in the ECMWF model. For
!   for the full levels mean average the two surrounding half levels.

!  Part of the Nudged model (see nudging_main.F90)

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
MODULE nudging_ecmwf_60level_def

IMPLICIT NONE

!Number of model levels in ecmwf model
INTEGER,PARAMETER :: ecmwf_nlevs  = 60
INTEGER,PARAMETER :: ecmwf_nplevs = 23

! ECMWF fixed pressure levels
REAL,PARAMETER :: data_pressure(1:ecmwf_nplevs) = (/                   &
 100000, 92500, 85000, 77500, 70000,                                   &
  60000, 50000, 40000, 30000, 25000,                                   &
  20000, 15000, 10000,  7000,  5000,                                   &
   3000,  2000,  1000,   700,   500,                                   &
    300,   200,  100                                                   &
 /)

!A(k) parameters for ECMWF half levels
REAL,PARAMETER :: ak(0:ecmwf_nlevs) =      (/    0.000000,             &
     20.000000,    38.425343,    63.647804,     95.636963,             &
    134.483307,   180.584351,   234.779053,    298.495789,             &
    373.971924,   464.618134,   575.651001,    713.218079,             &
    883.660522,  1094.834717,  1356.474609,   1680.640259,             &
   2082.273926,  2579.888672,  3196.421631,   3960.291504,             &
   4906.708496,  6018.019531,  7306.631348,   8765.053711,             &
  10376.126953, 12077.446289, 13775.325195,  15379.805664,             &
  16819.474609, 18045.183594, 19027.695312,  19755.109375,             &
  20222.205078, 20429.863281, 20384.480469,  20097.402344,             &
  19584.330078, 18864.750000, 17961.357422,  16899.468750,             &
  15706.447266, 14411.124023, 13043.218750,  11632.758789,             &
  10209.500977,  8802.356445,  7438.803223,   6144.314941,             &
   4941.778320,  3850.913330,  2887.696533,   2063.779785,             &
   1385.912598,   855.361755,   467.333588,    210.393890,             &
     65.889244,     7.367743,     0.000000,      0.000000              &
   /)

!B(k) parameters for ECMWF half levels
REAL,PARAMETER :: bk(0:ecmwf_nlevs) =(/ 0.00000000,                    &
    0.00000000, 0.00000000, 0.00000000, 0.00000000,                    &
    0.00000000, 0.00000000, 0.00000000, 0.00000000,                    &
    0.00000000, 0.00000000, 0.00000000, 0.00000000,                    &
    0.00000000, 0.00000000, 0.00000000, 0.00000000,                    &
    0.00000000, 0.00000000, 0.00000000, 0.00000000,                    &
    0.00000000, 0.00000000, 0.00000000, 0.00007582,                    &
    0.00046139, 0.00181516, 0.00508112, 0.01114291,                    &
    0.02067788, 0.03412116, 0.05169041, 0.07353383,                    &
    0.09967469, 0.13002251, 0.16438432, 0.20247594,                    &
    0.24393314, 0.28832296, 0.33515489, 0.38389215,                    &
    0.43396294, 0.48477158, 0.53570992, 0.58616841,                    &
    0.63554746, 0.68326861, 0.72878581, 0.77159661,                    &
    0.81125343, 0.84737492, 0.87965691, 0.90788388,                    &
    0.93194032, 0.95182151, 0.96764523, 0.97966272,                    &
    0.98827010, 0.99401945, 0.99763012, 1.00000000 /)

END MODULE nudging_ecmwf_60level_def
