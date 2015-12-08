! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! UKCA_RADAER look-up tables.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
! There are four look-up tables, for accumulation and coarse 
! mode aerosols (first dimension) and for SW and LW spectra 
! (second dimension).
!
MODULE ukca_radaer_lut
  
  USE ukca_radaer_tlut_mod
  
  IMPLICIT NONE
  SAVE
  
  INTEGER, PARAMETER :: npd_ukca_lut_mode  = 2
  INTEGER, PARAMETER :: IP_UKCA_LUT_ACCUM  = 1
  INTEGER, PARAMETER :: IP_UKCA_LUT_COARSE = 2
  
  INTEGER, PARAMETER :: npd_ukca_lut_spectrum = 2
!                 those two values must be consistent with
!                 the parameters in SPCRG3A
  INTEGER, PARAMETER :: IP_UKCA_LUT_SW        = 1
  INTEGER, PARAMETER :: IP_UKCA_LUT_LW        = 2
  
  TYPE (ukca_radaer_tlut), &
    DIMENSION(npd_ukca_lut_mode, npd_ukca_lut_spectrum) :: ukca_lut
  
END MODULE ukca_radaer_lut
