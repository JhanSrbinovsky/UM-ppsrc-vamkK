! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!   Contains:
!   Parameters used to increase the true surface albedo linearly from
!   a minimum value, alphaM at melting point, to a "cold-ice" value,
!   alphaC, over a temperature range, dTice
!   Mass Mixing Ratios of minor Gases N2O,CH4,CFC11,CFC12,O2
!   CFC113,HCFC22,HFC125, and HFC134A

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP

MODULE rad_com_mod

USE missing_data_mod, ONLY: imdi, rmdi

IMPLICIT NONE

! Number of columns used for internal sampling by the ISCCP simulator
INTEGER :: is_ncol  = imdi



! Albedo of sea-ice at melting point (TM) if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at melting point (TM) if l_ssice_albedo
REAL :: alpham = rmdi  ! "M" for "melting"

! Namelist input for alpham if l_ssice_albedo
! - assigned to alpham in readlsta.F90
REAL :: ssalpham = rmdi  ! "M" for "melting"

! Albedo of sea-ice at and below TM-DTICE if .not.l_ssice_albedo, or
! Albedo of snow on sea-ice at and below TM-DTICE if l_ssice_albedo
REAL :: alphac = rmdi   ! "C" for "cold"

! Namelist input for alphac if l_ssice_albedo
! - assigned to alphac in readlsta.F90
REAL :: ssalphac = rmdi  ! "C" for "cold"

! Albedo of snow-free sea-ice if l_ssice_albedo
REAL :: alphab = rmdi  ! "B" for "bare"

! Temperature range in which albedo of sea-ice, if .not.l_ssice_albedo,
! or of snow on sea-ice, if l_ssice_albedo, varies between its limits
REAL :: dtice= rmdi

! Namelist input for dtice if l_ssice_albedo
! - assignd to dtice in readlsta.F90
REAL :: ssdtice= rmdi

! Temperature range below TM over which meltponds form if
! l_sice_meltponds and l_ssice_albedo
REAL :: dt_bare= rmdi

! Increment to albedo for each degree temperature rises above
! TM-DT_BARE. Only used if l_sice_meltponds and l_ssice_albedo
REAL :: dalb_bare_wet= rmdi

! Fraction of SW radiation that penetrates seaice and scatters
! back causing an addition to the albedo. Only active if l_ssice_albedo
! and l_sice_scattering
REAL :: pen_rad_frac= rmdi

! attenutation parameter for SW in seaice which controls the
! additional albedo due to internal scattering
REAL :: sw_beta= rmdi
REAL :: n2ommr = rmdi ! N2O mmr (Mass Mixing Ratio)
REAL :: ch4mmr = rmdi ! CH4 mmr
REAL :: c11mmr = rmdi ! CFC11 mmr
REAL :: c12mmr = rmdi ! CFC12 mmr
REAL :: o2mmr  = rmdi ! O2 mmr
REAL :: c113mmr    = rmdi ! CFC113 mmr
REAL :: c114mmr    = rmdi ! CFC114 mmr
REAL :: hcfc22mmr  = rmdi ! HCFC22 mmr
REAL :: hfc125mmr  = rmdi ! HFC125 mmr
REAL :: hfc134ammr = rmdi ! HFC134A mmr


END MODULE rad_com_mod
