! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! 
! Module to hold the parameters which are used in the tcs warm rain
! routines
!
MODULE tcs_parameters_warm

  IMPLICIT NONE
  !
  ! Description:
  !   Module to hold the parameters which are used in the tcs warm rain
  !   routines
  !
  ! Method:
  !   Most of this is a straight forward declaration of parameters, 
  !   except for the derived type "similarity_coefficient" which 
  !   defines a complete set of parameters to be used by different
  !   types of convection (e.g. shallow, congestus). 
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !   Language: Fortran 90.
  !   This code is written to UMDP3 version 8.1 programming standards.
  !

  INTEGER :: icong_options=1 !options to select methods in the code

! Used in:  calc_scales

  REAL ::                   &
    mb_a1  = 0.03           &  ! Cloud base mass flux parameter
  , mb_a2  = 0.2            &  ! Cloud base mass flux parameter
  , wup_a1 = 1.1               ! Cloud base updraught velocity parameter
  

! Used in:  cgconv_precip

  REAL :: ql_t=0.001 
                               ! threshold for precipitation
  REAL :: epsilon_rain=0.5 
                               ! coefficient for surface
                               ! precipitation calculation.
  REAL :: beta_rain=0.5    
                               ! coefficient for precip
                               ! production.
  REAL :: gamma_rain=0.5   
                               ! coefficient for rainfall
                               ! rate at base of inversion.
  REAL :: kappa_rain=1.0   
                               ! coefficient for rainfall
                               ! rate at base of inversion.

! Used in:  cgconv_cloudbase

  REAL :: wql_cb_a1=0.4
                               ! coefficient for cloudbase
                               ! liquid water flux
  REAL :: jump_cb_a1=0.7
                               ! cloud base jump coefficient 
                               ! wthv_- = jump_cb_a1*mb*dthv
  REAL :: wthv_cb_a1=0.16
                               ! coefficient for cloudbase 
                               ! thetav flux
                               ! wthv_- = wthv_cb_a1*wthv_surf
  REAL :: wh_cb_a1=1.0
                               ! coefficient for cloudbase 
                               ! moist static energy flux
                               ! wh_- = wthv_cb + wh_cb_a1
                               !        *(l/cp-.61*TH)*wqv_surf 
  REAL :: sat_a1=0.4
                               ! Saturation condition coefficient 

! Used in:  cgconv_inversion

  REAL :: beta_cld = 1.1
  REAL :: jump_inv_a1=0.05 !0.19
                               ! inversion jump coefficient 
                               ! wthv_- = jump_inv_a1*mb*dthv_inv
  ! Working note: The value for jump_inv_a2 (moist static energy flux at
  ! inversion) needs to be checked out further from the LEM simulations.
  REAL :: jump_inv_a2=1.0
                               ! inversion jump coefficient 
                               ! wh_- = jump_inv_a2*mb*dh_inv
  ! Working note: The value for jump_inv_a3 (wthetav at
  ! inversion) needs to be checked out further from the LEM simulations.
  REAL :: jump_inv_a3=.385
                               ! inversion jump coefficient 
                               ! wthetav_- = jump_inv_a3*mb*dthetav_inv
  REAL :: mb_inv_a1=0.35
                               ! ratio of inversion mass
                               ! flux to cloudbase mass flux
                               ! m_inv=mb_inv_a1*mb
  REAL :: sat_inv_a1=0.85
                               ! Inversion Saturation
                               ! condition coefficient
  REAL :: q_min=1.e-8
                               ! Minimum allowed value of water vapour
                               ! mixing ratio at inversion base

! Used in:  cong_similarity

  TYPE :: similarity_coefficient
     
     INTEGER :: conv_type
     REAL :: K(2)
     REAL :: Fng(3)
     REAL :: B(3)
     REAL :: g(5)
     REAL :: fql(3)
     REAL :: gql(2)
     REAL :: f0(2)
     REAL :: f1(2)
     REAL :: fth(3)
     REAL :: fw(6)
     REAL :: pc2_d(3)

  END TYPE similarity_coefficient
  
  ! Coefficients to be used for shallow non-precipitating convection
  ! , i.e. conv_type=1
  TYPE (similarity_coefficient), TARGET ::              &
     sim_coeff_nonp_sh = similarity_coefficient(                   &
     1,                                                            & !type
     (/0.115,    10.                                           /), & !K
     (/0.3,       1.75,     10.                                /), & !Fng
     (/1.4,       0.6,      10.                                /), & !B
     (/1.45,      0.665,     4.67/36.,  1.33/6.,   6.0         /), & !g
     (/0.6,       0.6,       5.0                               /), & !fql
     (/0.6,       4.5                                          /), & !gql
     (/0.473,     0.64                                         /), & !f0
     (/9.0,       1.5                                          /), & !f1
     (/2.25,      1.75,      4.0                               /), & !fth
     (/-0.4273/3., 0.9752, 0.06329,1e-8, 13.819/3., -5.5197/3. /), & !fw
     (/0.8       , 0.0005    ,0.0005                           /)  & !pc2_detr
     )

  ! Coefficients to be used for shallow precipitating convection
  ! , i.e. conv_type=2
  ! Currently this is set to be the same as non-precipitating
  TYPE (similarity_coefficient), TARGET ::              &
     sim_coeff_prec_sh = similarity_coefficient(                   &
     2,                                                            & !type
     (/0.115,    10.                                           /), & !K
     (/0.3,       1.75,     10.                                /), & !Fng
     (/1.4,       0.6,      10.                                /), & !B
     (/1.45,      0.665,     4.67/36.,  1.33/6.,   6.0         /), & !g
     (/0.6,       0.6,       5.0                               /), & !fql
     (/0.6,       4.5                                          /), & !gql
     (/0.473,     0.64                                         /), & !f0
     (/9.0,       1.5                                          /), & !f1
     (/2.25,      1.75,      4.0                               /), & !fth
     (/-0.4273/3., 0.9752, 0.06329,1e-8, 13.819/3., -5.5197/3. /), & !fw
     (/0.8,       0.0005,    0.0005                            /)  & !pc2_detr
     )

  
  ! Coefficients to be used for warm congestus convection
  ! , i.e. conv_type=3
  TYPE (similarity_coefficient), TARGET ::              &
     sim_coeff_warm_cg = similarity_coefficient(                   &
     3,                                                            & !type
     (/0.08,     20.                                           /), & !K
     (/0.28,      1.75,     20.                                /), & !Fng
     (/3.,        1.333,    10.                                /), & !B
     (/1.45,      0.665,     4.67/36.,  1.33/6.,   6.0         /), & !g
     (/1.0,       0.6,       5.0                               /), & !fql
     (/0.6,       4.5                                          /), & !gql
     (/0.473,     0.64                                         /), & !f0
     (/9.0,       1.5                                          /), & !f1
     (/2.25,      1.75,      4.0                               /), & !fth
     (/-0.4273,   0.9752,    0.06329,1e-8, 13.819, -5.5197     /), & !fw
     (/0.8,       0.0005,    0.0005                            /)  & !pc2_detr
     )

  ! Coefficients to be used for mixed-phase congestus convection
  ! , i.e. conv_type=4
  ! Currently this is set to be the same as warm  
  TYPE (similarity_coefficient), TARGET ::              &
     sim_coeff_cold_cg = similarity_coefficient(                   &
     4,                                                            & !type
     (/0.08,     20.                                           /), & !K
     (/0.28,      1.75,     20.                                /), & !Fng
     (/3.,        1.333,    10.                                /), & !B
     (/1.45,      0.665,     4.67/36.,  1.33/6.,   6.0         /), & !g
     (/1.0,       0.6,       5.0                               /), & !fql
     (/0.6,       4.5                                          /), & !gql
     (/0.473,     0.64                                         /), & !f0
     (/9.0,       1.5                                          /), & !f1
     (/2.25,      1.75,      4.0                               /), & !fth
     (/-0.4273,   0.9752,    0.06329,1e-8, 13.819, -5.5197     /), & !fw
     (/0.8,       0.0005,    0.0005                            /)  & !pc2_detr
     )


! Used in: cgtconv_base_stress
  REAL, PARAMETER ::            &
       beta_cmt  = 0.04         &  ! Coefficient for CMT calculation
       ,delta_cmt = 2.3         &  ! Coefficient for CMT calculation
       ,gamma_cmt = 1.63           ! Coefficient for CMT calculation

! Used in:  cgconv_turb_fluxes

  REAL, PARAMETER :: wthvl_factor = 1.

! Used in: calc_pc2
  
  REAL, PARAMETER :: Aepsilon = .03 ! fraction of TKE production available for entrainment

! Used in: tcs_warm_mod

  REAL :: asmooth_inv = 2.0  ! Smoothing parameter for inversion 
  REAL :: bsmooth_inv = 1.0  ! Smoothing parameter for inversion 

END MODULE tcs_parameters_warm
