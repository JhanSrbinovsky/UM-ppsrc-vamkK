! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!+ Global data module for switches/options concerned with convection.
 
MODULE cv_param_mod

  ! Description:
  !   Module containing parameters used by the convection code.
  !
  ! Method:
  !   Default values have been declared where appropriate.
  !   
  !   Any routine wishing to use these options may do so with the 'Use' 
  !   statement.
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 v7.4 programming standards.
  !
  ! Declarations:

  IMPLICIT NONE
  SAVE

  !===========================================================================
  ! Integer parameters to remove use of magic numbers
  !===========================================================================

  ! Basis for 2d convective cloud amount
  !   Total Condensed Water(TCW) original 4A
  !   LES mb/wsc scalings (Grant and Lock, 2004) (shallow)
  !   Surface rain from respective convective cloud type (Dp and Md-level)
  Integer, Parameter :: cca2d_total_condensed_water = 0
  Integer, Parameter :: cca2d_grant_lock            = 1 ! Shallow cnv (CCRad)
  Integer, Parameter :: cca2d_srf_precip_rate       = 1 ! mid/deep cnv (CCRad)

  Integer, Parameter :: total_condensed_water = 0
  Integer, Parameter :: grant_lock            = 1 ! Shallow cnv (CCRad)
  Integer, Parameter :: srf_precip            = 1 ! mid/deep cnv (CCRad)

  ! Convective cloud decay
  Integer, Parameter :: rad_decay_off           = 0
  Integer, Parameter :: rad_decay_full_timestep = 1
  Integer, Parameter :: rad_decay_conv_substep  = 2 ! CCRad only


  ! Convective cloud decay timescale
  Integer, Parameter :: cld_life_constant = 0
  Integer, Parameter :: cld_life_func_hgt = 1 ! CCRad only

  ! For all the options below except 3), the anvil base is at the freezing
  ! Level
  Integer, Parameter :: anv_pressure      = 0
  Integer, Parameter :: anv_height        = 1
  Integer, Parameter :: anv_model_levels  = 2
  Integer, Parameter :: anv_limited_pressure_depth = 3


  
  ! Parameters used to relate cca_2d of cld to its precipitation rate
  Real, Parameter    :: a_land = 0.3   
  Real, Parameter    :: a_sea  = 0.3   
  Real, Parameter    :: b_land = 0.025 
  Real, Parameter    :: b_sea  = 0.025 

  ! Application of convective cloud anvils
  Real, parameter :: deep_dp = 50000.0  ! Min. depth for anvil criteria(Pa)

  ! Critical depth of cloud for the formation of
  ! convective precipitation over sea (m)
  REAL, PARAMETER :: critdsea = 1.5E3

  ! critical depth of cloud for the formation of convective
  ! precipitation over land (m)
  REAL, PARAMETER :: critdlnd = 4.0E3

  ! critical depth of a glaciated cloud for the formation of
  ! convective precipitation (m)
  REAL, PARAMETER :: critdice = 1.0E3

  ! Parcel ascent in diagnosis, cloud water condensate
  REAL, PARAMETER :: qlcrit = 1.0e-3   ! critical cloud water
  
  ! Timestep frequency for calling convection - hardwired to call every timestep
  INTEGER,PARAMETER :: a_conv_step = 1 

!============================================================================
! Parcel ascent
!============================================================================

  ! coefficients used in calculation of entrainment rate
  REAL, PARAMETER :: ae1     = 1.0        ! Not used 
  REAL, PARAMETER :: ae2     = 1.5

  ! minimum excess buoyancy to continue parcel ascent (K)

  REAL,PARAMETER:: xsbmin = 0.2       ! Used in mid scheme

  ! initial excess potential temperature (k) and mixing ratio
  ! (kg/kg) for deep convection
  REAL, PARAMETER :: thpixs_deep= 0.2
  REAL, PARAMETER :: qpixs_deep =0.0

  ! initial excess potential temperature (k) and mixing ratio
  ! (kg/kg) for shallow convection
  REAL, PARAMETER :: thpixs_shallow = 0.2
  REAL, PARAMETER :: qpixs_shallow  = 0.0

  ! initial excess potential temperature (k) and mixing ratio
  ! (kg/kg) for mid-level convection
  REAL, PARAMETER :: thpixs_mid= 0.2
  REAL, PARAMETER :: qpixs_mid =0.0

  ! Minimum parcel buoyancy/layer thickness (K/Pa)
  REAL, PARAMETER :: mparb = 1.0            ! Used in mid scheme

  ! Constants to determine initial convective mass flux from parcel buoyancy
  ! Deep convection
  REAL, PARAMETER :: c_deep = 5.17E-4       ! No longer used
  REAL, PARAMETER :: d_deep = 0.0           ! No longer used

  ! Shallow convection
  REAL, PARAMETER :: c_shallow = 5.17E-4    ! No longer used
  REAL, PARAMETER :: d_shallow = 0.0        ! No longer used

  ! Mid convection
  REAL, PARAMETER :: c_mid = 5.17E-4
  REAL, PARAMETER :: d_mid = 0.0

  ! limits on the initial convective parcel perturbations
  REAL, PARAMETER :: max_diag_thpert  =  2.0
  REAL, PARAMETER :: max_dp_thpert    =  2.0
  REAL, PARAMETER :: min_dp_thpert    = -2.0
  REAL, PARAMETER :: max_sh_thpert    =  2.0
  REAL, PARAMETER :: min_sh_thpert    = -2.0
  REAL, PARAMETER :: max_dp_qpert_fac =  0.2
  REAL, PARAMETER :: max_sh_qpert_fac =  0.2

  ! mparfl = 1E-3 * minimum parcel buoyancy * mass flux parameter c
  REAL, PARAMETER :: mparfl = 1.0E-3 * 1.0 * 3.33E-4

  ! Difference in potential temperature between levels above which the 
  ! atmosphere is assumed to be too stable to convect (K)
  REAL, PARAMETER :: delthst = 0.5

!============================================================================
! Convective closure
!============================================================================

  ! Coefficient relating sub-cloud convective velocity scale to cumulus
  ! mass flux for shallow convection
  REAL, PARAMETER :: c_mass=0.03

  ! Tuneable factor used in denominator of W_CAPE timescale calculation 
  REAL, PARAMETER :: wcape_fac = 3.0

  ! Parameter governing the speed of transition from parametrized to 
  ! resolved convection, for use with sh_grey option
  REAL :: beta_cu = 0.15 
!============================================================================
! Downdraught and evaporation below cloud base calculations
!============================================================================

  ! Coefficients used in calculation of downdraught entrainment
  ! rates
  REAL, PARAMETER :: ddcoef1 = 1.8E6
  REAL, PARAMETER :: ddcoef2 = 3.0

  ! Thickness level used in calculation of mixing detrainment for
  ! downdraught  (pa)
  REAL, PARAMETER :: det_lyr = 10000.0

  ! exponents used in calculation of evaporation of liquid
  REAL, PARAMETER :: p_lq1 = 0.52
  REAL, PARAMETER :: p_lq2 = 0.67

  ! exponents used in calculation of evaporation of ice
  REAL, PARAMETER :: p_ice1 = 0.55
  REAL, PARAMETER :: p_ice2 = 0.76

  ! exponents and constants associated with density term in
  ! evaporation of liquid
  REAL, PARAMETER :: rho_lqp1 = 0.26
  REAL, PARAMETER :: rho_lqp2 = 0.59
  REAL, PARAMETER :: rho_lqa  = 108.80
  REAL, PARAMETER :: rho_lqb  = 830.73

  ! exponents and constants associated with density term in
  ! evaporation of ice
  REAL, PARAMETER :: rho_icp1 = 0.28
  REAL, PARAMETER :: rho_icp2 = 0.63
  REAL, PARAMETER :: rho_icea = 1569.52
  REAL, PARAMETER :: rho_iceb = 32069.02

  ! constants used in quadratic formula for evaporation of liquid
  REAL, PARAMETER :: lq_a = 2.008E-9
  REAL, PARAMETER :: lq_b = -1.385E-6
  REAL, PARAMETER :: lq_c = 2.424E-4

  ! constants used in quadratic formula for evaporation of ice
  REAL, PARAMETER :: ice_a = -5.2E-9
  REAL, PARAMETER :: ice_b = 2.5332E-6
  REAL, PARAMETER :: ice_c = -2.911E-4

  ! Fractional cloud area when no dd, used in evaporation calc
  REAL, PARAMETER :: cldarea = 1.0

  ! fractional cloud area of dd
  REAL, PARAMETER :: ddcldfra = 0.5 

  ! downdraught precipitation transfer efficiency factor
  REAL, PARAMETER :: ddptef = 2.0 

!============================================================================
! Convective momentum transport (CMT) calculations
!============================================================================

  ! Coefficient for "pressure term" relative to "shear term" as found in 
  ! Kershaw & Gregory 1997

  REAL, PARAMETER :: cpress_term = 0.7   


END MODULE cv_param_mod
