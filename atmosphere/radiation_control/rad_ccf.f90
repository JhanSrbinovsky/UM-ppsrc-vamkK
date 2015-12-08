! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Radiation Constant Configuration File

MODULE rad_ccf

! Description:
!   Module containing settings of standard constants used within the
!   radiation scheme. This replaces a number of separate constant
!   configuration files. The original file names are indicated at the
!   head of each section.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: Fortran 95
!   This code is written to UMDP3 v8.3 programming standards.

  USE conversions_mod, ONLY: pi
  USE atmos_constants_mod, ONLY: r, repsilon, c_virtual

  IMPLICIT NONE

! diff_elsasser_ccf, elsass3a
! ------------------------------------------------------------------
! Module to set diffusivity for elsasser's scheme.
  REAL, PARAMETER :: elsasser_factor = 1.66e+00
!   Diffusivity factor for elsasser's scheme

! ------------------------------------------------------------------
! diff_keqv_ucf, diffke3a
! ------------------------------------------------------------------
! Module to set the diffusivity factor for use with equivalent
! extinction
  REAL, PARAMETER :: diffusivity_factor_minor = 1.66e+00
!   Minor diffusivity factor

! ------------------------------------------------------------------
! physical_constants_0_ccf, phycn03a
! ------------------------------------------------------------------
! Module setting physical constants.
  REAL, PARAMETER :: mol_weight_air = 28.966e-03
!   Molar weight of dry air
  REAL, PARAMETER :: n2_mass_frac   = 0.781e+00
!   Mass fraction of nitrogen

! ------------------------------------------------------------------
! physical_constants_pp_ccf
! ------------------------------------------------------------------
! Module setting physical constants.
  REAL, PARAMETER :: n_avogadro = 6.022045e+23
!   Avogadro's number
  REAL, PARAMETER :: k_boltzmann = 1.380662e-23
!   Boltzmann's constant

END MODULE rad_ccf
