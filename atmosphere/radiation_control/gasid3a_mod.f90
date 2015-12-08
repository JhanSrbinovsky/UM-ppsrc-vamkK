! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
MODULE gasid3a_mod

IMPLICIT NONE
! GASID3A defines indexing numbers of gaseous absorbing species for
! two-stream radiation code.
! the numbering 1-12 corresponds to lowtran 7.

      INTEGER,PARAMETER:: NPD_GASES=20 ! Number of indexed gases

      INTEGER,PARAMETER:: IP_H2O=1
      INTEGER,PARAMETER:: IP_CO2=2
      INTEGER,PARAMETER:: IP_O3=3
      INTEGER,PARAMETER:: IP_N2O=4
      INTEGER,PARAMETER:: IP_CO=5
      INTEGER,PARAMETER:: IP_CH4=6
      INTEGER,PARAMETER:: IP_O2=7
      INTEGER,PARAMETER:: IP_NO=8
      INTEGER,PARAMETER:: IP_SO2=9
      INTEGER,PARAMETER:: IP_NO2=10
      INTEGER,PARAMETER:: IP_NH3=11
      INTEGER,PARAMETER:: IP_HNO3=12
      INTEGER,PARAMETER:: IP_N2=13
      INTEGER,PARAMETER:: IP_CFC11=14
      INTEGER,PARAMETER:: IP_CFC12=15
      INTEGER,PARAMETER:: IP_CFC113=16
      INTEGER,PARAMETER:: IP_HCFC22=17
      INTEGER,PARAMETER:: IP_HFC125=18
      INTEGER,PARAMETER:: IP_HFC134A=19
      INTEGER,PARAMETER:: IP_CFC114=20

! GASID3A end

END MODULE gasid3a_mod
