! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Gravity Wave Drag
MODULE c_gwave_mod

IMPLICIT NONE
!
!  Description: This comdeck defines the constants for the 4A and 5A
!               versions of the orographic gravity wave drag scheme. These are
!               tuneable parameters but are unlikely to be changed.

! Orographic drag scheme parameters:
! All schemes
      ! Number of standard deviations above the mean orography of top
      ! of sub-grid mountains
      REAL,PARAMETER     ::  nsigma       = 2.5

! 4A scheme only
      ! Switch to determine which wave saturation test is used
      INTEGER,PARAMETER :: Amplitude_saturation = 1
      INTEGER,PARAMETER :: Stress_saturation    = 0

! 5A scheme only
      ! Fixed value for group velocity angle (used when l_nonhydro true 
      ! and l_dynbeta false)
      REAL,PARAMETER     ::  beta_fix     = 100.   
      ! Fraction of vert wavelength to deposit acceleration over 
      ! (when l_smooth true)  
      REAL,PARAMETER     ::  frac_wl      = 1.0      
      ! Minimum value for local vert wavelength (U/N) (m)                  
      REAL,PARAMETER     ::  lambdaz_min  = 250. 
      ! Maximum value for local vert wavelength (U/N) (m)     
      REAL,PARAMETER     ::  lambdaz_max  = 3000.  
      ! Threshold for brunt-vaisala frequency squared (ensq) 
      ! to define a neutral layer (used in gw_block and gw_wave)
      REAL, PARAMETER    ::  nsq_neutral  = 1.E-6  
      ! Threshold for % agreement for Zav convergence in gw_block
      REAL, PARAMETER    ::  zav_converge = 0.05
      ! Num of iterations to allow Zav to converge in gw_block
      INTEGER, PARAMETER ::  zav_iterate  = 5      
!

END MODULE c_gwave_mod
