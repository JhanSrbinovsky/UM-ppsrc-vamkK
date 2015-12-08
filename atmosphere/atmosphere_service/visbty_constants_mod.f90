! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   A module containing constants used in visibility diagnosis
!
MODULE visbty_constants_mod

USE conversions_mod, ONLY: pi
USE  missing_data_mod, ONLY: rmdi
!
! Description:
!   This module contains declarations for constants used to diagnose
!   horizontal visibility.
!  
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Atmosphere Service
!
! Code description:
!   Language: Fortran 90
!   This code is written to UM programming standards version 8.2.
!

IMPLICIT NONE

!
! General constants (these should be promoted to a more general module)
!
REAL, PARAMETER :: OneThird   = 1.0/3.0
REAL, PARAMETER :: FourThirds = 4*OneThird

!
! Information about the diagnostics themselves
!      
INTEGER, PARAMETER :: n_vis_thresh = 2
REAL, PARAMETER    :: visfog = 1000.0      ! Visibility defining fog
REAL, PARAMETER    :: vismist = 5000.0     ! Visibility defining mist
REAL, PARAMETER    :: vis_thresh(n_vis_thresh)=(/ visfog, vismist /)
REAL               :: calc_prob_of_vis = rmdi
                    ! tunable parameter: the cumulative prob value at 
                    ! which vis is estimated, set in RUN_BL namelist

!
! Parameters used in visibility calculations
!
REAL, PARAMETER :: LiminalContrast = 0.05
REAL, PARAMETER :: LnLiminalContrast = -2.99573
                                        ! Natural log of Liminal contrast
REAL, PARAMETER :: RecipVisAir = 1.0E-5 ! Reciprocal of the clean air 
                                        ! visibility (100km)

!
! Microphysical parameters used in visibility calculations
!
REAL, PARAMETER :: N0 = 200.0E7         ! Standard number density of murk (/m3)
REAL, PARAMETER :: B0= 0.14             ! Activation parameter
REAL, PARAMETER :: radius0 = 0.11E-6    ! Radius of standard murk particle (m)
REAL, PARAMETER :: rho_aerosol = 1700.0 ! Density of murk (kg/m3)
REAL, PARAMETER :: rho_air = 1.0        ! Density of air (kg/m3)

REAL, PARAMETER :: beta0  = 1.5 * Pi    ! Scattering coefficient normalisation

REAL, PARAMETER :: a0 = 1.2E-9          ! Constant involving surface energy of water
REAL, PARAMETER :: m0 = FourThirds * Pi * radius0 * radius0 * radius0           &
                        * (rho_aerosol/rho_air) * N0
                                        ! Standard aerosol mass mixing ratio (kg/kg)

REAL, PARAMETER :: power = 1.0/6.0      ! Murk particle radius/mass loading power

REAL, PARAMETER :: VisFactor = -LnLiminalContrast / Beta0
                                        ! transformation to visibility 
                                        ! ( = ln(liminal contrast) / Beta0 )

REAL, PARAMETER :: aero0 = 0.1          ! Minimum allowed murk aerosol
REAL, PARAMETER :: aeromax = 200.0      ! Maximum allowed murk aerosol

END MODULE visbty_constants_mod
