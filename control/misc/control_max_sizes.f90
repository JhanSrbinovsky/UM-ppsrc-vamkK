! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc.

! CMAXSIZE maximum sizes for dimensioning arrays
! of model constants whose sizes are configuration dependent. This
! allows constants to be read in from a NAMELIST file and maintain
! the flexibility of dynamic allocation for primary variables. The
! maximum sizes should agree with the maximum sizes implicit in the
! front-end User Interface.

MODULE Control_Max_Sizes

  USE Atmos_Max_Sizes, ONLY : model_levels_max
  IMPLICIT NONE

  ! Max levels in boundary layer
  INTEGER, PARAMETER :: max_bl_levels = model_levels_max
  
  ! Max size of alpha_Cd
  INTEGER, PARAMETER :: max_number_alpha_cds = max_bl_levels
  
  ! Max no. of levels for pvort output
  INTEGER, PARAMETER :: MAX_REQ_THPV_LEVS = model_levels_max

  ! Max no. 1-2-1 rows in polar filter
  INTEGER, PARAMETER ::  max_121_rows =  32
  ! 0 is used for horizontal diffusion pointer

  ! Max no. of levels (from top) to apply upper level diffusion
  INTEGER, PARAMETER ::  max_updiff_levels = 10

  ! Max size of any sponge zones
  INTEGER, PARAMETER ::  max_sponge_width = 10

  ! Max size of look-up tables for searches
  INTEGER, PARAMETER ::  max_look = 2048

  ! Max no. of atmos interface areas
  INTEGER, PARAMETER :: MAX_N_INTF_A =  8

  ! Max no. of points in LBC      
  INTEGER, PARAMETER :: max_intf_lbcrow_length = 5000
  INTEGER, PARAMETER :: max_intf_lbcrows = 5000

  ! Max no. of atmos interface levels
  INTEGER, PARAMETER :: MAX_INTF_LEVELS = model_levels_max

  ! Maximum number of physics segments
  INTEGER, PARAMETER :: MAX_NO_OF_SEGS = 200

END MODULE Control_Max_Sizes
