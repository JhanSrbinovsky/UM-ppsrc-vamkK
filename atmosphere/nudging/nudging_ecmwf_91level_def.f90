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
!  pressure levels for (new) 91 level version.
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
MODULE nudging_ecmwf_91level_def

IMPLICIT NONE

!Number of model levels in ecmwf model
INTEGER,PARAMETER :: ecmwf_nlevs  = 91
INTEGER,PARAMETER :: ecmwf_nplevs = 23

! ECMWF fixed pressure levels
REAL,PARAMETER :: data_pressure(1:ecmwf_nplevs) =(/                    &
 100000, 92500, 85000, 77500, 70000,                                   &
  60000, 50000, 40000, 30000, 25000,                                   &
  20000, 15000, 10000,  7000,  5000,                                   &
   3000,  2000,  1000,   700,   500,                                   &
    300,   200,  100                                                   &
 /)

!A(k) parameters for ECMWF half levels
REAL,PARAMETER :: ak(0:ecmwf_nlevs) =(/                                &
     0.000000,     2.000040,     3.980832,     7.387186,               &
    12.908319,    21.413612,    33.952858,    51.746601,               &
    76.167656,   108.715561,   150.986023,   204.637451,               &
   271.356506,   352.824493,   450.685791,   566.519226,               &
   701.813354,   857.945801,  1036.166504,  1237.585449,               &
  1463.163940,  1713.709595,  1989.874390,  2292.155518,               &
  2620.898438,  2976.302246,  3358.425781,  3767.196045,               &
  4202.416504,  4663.776367,  5150.859863,  5663.156250,               &
  6199.839355,  6759.727051,  7341.469727,  7942.926270,               &
  8564.624023,  9208.305664,  9873.560547, 10558.881836,               &
 11262.484375, 11982.662109, 12713.897461, 13453.225586,               &
 14192.009766, 14922.685547, 15638.053711, 16329.560547,               &
 16990.623047, 17613.281250, 18191.029297, 18716.968750,               &
 19184.544922, 19587.513672, 19919.796875, 20175.394531,               &
 20348.916016, 20434.158203, 20426.218750, 20319.011719,               &
 20107.031250, 19785.357422, 19348.775391, 18798.822266,               &
 18141.296875, 17385.595703, 16544.585938, 15633.566406,               &
 14665.645508, 13653.219727, 12608.383789, 11543.166992,               &
 10471.310547,  9405.222656,  8356.252930,  7335.164551,               &
  6353.920898,  5422.802734,  4550.215820,  3743.464355,               &
  3010.146973,  2356.202637,  1784.854614,  1297.656128,               &
   895.193542,   576.314148,   336.772369,   162.043427,               &
    54.208336,     6.575628,     0.003160,     0.000000                &
  /)

!      !B(k) parameters for ECMWF half levels
REAL,PARAMETER :: bk(0:ecmwf_nlevs) =(/0.000000,                       &
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,                  &
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,                  &
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,                  &
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,                  &
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,                  &
    0.000000, 0.000000, 0.000000, 0.000000, 0.000000,                  &
    0.000000, 0.000000, 0.000000, 0.000000, 0.000014,                  &
    0.000055, 0.000131, 0.000279, 0.000548, 0.001000,                  &
    0.001701, 0.002765, 0.004267, 0.006322, 0.009035,                  &
    0.012508, 0.016860, 0.022189, 0.028610, 0.036227,                  &
    0.045146, 0.055474, 0.067316, 0.080777, 0.095964,                  &
    0.112979, 0.131935, 0.152934, 0.176091, 0.201520,                  &
    0.229315, 0.259554, 0.291993, 0.326329, 0.362203,                  &
    0.399205, 0.436906, 0.475016, 0.513280, 0.551458,                  &
    0.589317, 0.626559, 0.662934, 0.698224, 0.732224,                  &
    0.764679, 0.795385, 0.824185, 0.850950, 0.875518,                  &
    0.897767, 0.917651, 0.935157, 0.950274, 0.963007,                  &
    0.973466, 0.982238, 0.989153, 0.994204, 0.997630,                  &
    1.000000                                                           &
 /)

END MODULE nudging_ecmwf_91level_def
