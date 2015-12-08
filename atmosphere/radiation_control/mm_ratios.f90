! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!  Define the reference mixing ratios for the radiatively active
!  gases for the diagnostic call in order to calculate the
!  radiative forcings.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.

MODULE mm_ratios

IMPLICIT NONE

!     Values taken from adutv by WJI 31/3/04.
      REAL, PARAMETER :: co2_mmr_d = 4.34800e-04
      REAL, PARAMETER :: n2o_mix_ratio_d = 4.205e-07
      REAL, PARAMETER :: ch4_mix_ratio_d = 4.461e-07
      REAL, PARAMETER :: o2_mmr_d = 0.2314
      REAL, PARAMETER :: cfc11_mix_ratio_d = 0.0
      REAL, PARAMETER :: cfc12_mix_ratio_d = 0.0
      REAL, PARAMETER :: c113mmr_d    = 0.0
      REAL, PARAMETER :: hcfc22mmr_d  = 0.0
      REAL, PARAMETER :: hfc125mmr_d  = 0.0
      REAL, PARAMETER :: hfc134ammr_d = 0.0
END MODULE mm_ratios
