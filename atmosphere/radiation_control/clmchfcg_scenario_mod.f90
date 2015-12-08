! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   A module containing tunable parameters used for histories/scenarios 
!   of climate change forcings (clm ch fcg).
!
MODULE clmchfcg_scenario_mod

!
! Description:
!   This module contains declarations for tunable parameters
!   used for histories/scenarios of climate change forcings.
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Radiation Control
!
! Code description:
!   Language: Fortran 90
!   This code is written to UM programming standards UMDP 3 vn8.2
!
USE missing_data_mod, ONLY: imdi, rmdi

IMPLICIT NONE

  ! Number of well-mixed greenhouse gases
  INTEGER, PARAMETER :: Nwmghg   = 10

  ! Number of sulphate loading patterns
  INTEGER, PARAMETER :: Nsulpat  = 2

  ! Maximum length of scenarios
  INTEGER, PARAMETER :: LenScen  = 50

  ! Number of such scenarios, made up of:
  INTEGER, PARAMETER :: Nclmfcgs = Nwmghg + Nsulpat

  !  Indices indicating which scenario corresponds to which forcing:
  INTEGER, PARAMETER :: S_CO2     = 1
  INTEGER, PARAMETER :: S_CH4     = 2
  INTEGER, PARAMETER :: S_N2O     = 3
  INTEGER, PARAMETER :: S_CFC11   = 4
  INTEGER, PARAMETER :: S_CFC12   = 5
  INTEGER, PARAMETER :: S_SO4     = 6
  INTEGER, PARAMETER :: S_CFC113  = 8
  INTEGER, PARAMETER :: S_HCFC22  = 9
  INTEGER, PARAMETER :: S_HFC125  = 10
  INTEGER, PARAMETER :: S_HFC134A = 11
  INTEGER, PARAMETER :: S_CFC114  = 12

  !  Carbon dioxide (CO2), methane (CH4), nitrous oxide (N2O),
  !  trichlorofluoromethane (CCl3F, "CFC-11"),
  !  dichlorodifluoromethane (CCl2F2, "CFC-12"), and then the first
  !  HadCM2-style anthropogenic sulphate loading pattern - these
  !  come at the end as their number in principle may vary.

  ! Switch for time varying GHG
  LOGICAL     :: L_ClmChFcg = .FALSE. 

  ! Years at which a rate or level is specified
  INTEGER     :: Clim_Fcg_Years(LenScen,Nclmfcgs) = imdi

  ! Number of such years, for each forcing
  INTEGER     :: Clim_Fcg_NYears(Nclmfcgs) = imdi

  ! Values, or rates of increase, for the designated years.
  !  See GAS_CALC (in Section 70) or the umui panels for details.
  REAL        :: Clim_Fcg_Levls(LenScen,Nclmfcgs) = rmdi
  REAL        :: Clim_Fcg_Rates(LenScen,Nclmfcgs) = rmdi

  NAMELIST / clmchfcg / Clim_Fcg_NYears, Clim_Fcg_Years,            &
                        Clim_Fcg_Levls, Clim_Fcg_Rates, L_ClmChFcg

END MODULE clmchfcg_scenario_mod
