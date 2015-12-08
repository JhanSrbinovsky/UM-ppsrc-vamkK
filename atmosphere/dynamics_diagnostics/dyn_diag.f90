! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate various diagnostics related to dynamics variables.
!
! Subroutine Interface:
      SUBROUTINE Dyn_diag(                                              &
! Primary data: in
     & exner_rho_levels,rho,u,v,w                                       &
     &,exner_theta_levels                                               &
     &,theta                                                            &
! Grid sizes and definition: in
     &,rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &,global_rows,global_row_length                                    &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,eta_theta_levels,eta_rho_levels                                  &
     &,Model_domain                                                     &
! Grid coordinates: in
     &,delta_lambda,delta_phi                                           &
     &,ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole      &
! Pre-calculated grid associated arrays: in
     &,r_at_u,r_at_v                                                    &
     &, r_theta_levels, r_rho_levels, sec_v_latitude                    &
     &, tan_v_latitude, sec_theta_latitude, f3_at_v                     &
     &,rot_coeff1,rot_coeff2                                            &
! Time information: in
     &,forecast_hrs                                                     &
! Theta levels for output arrays
     &, desired_theta                                                   &
! Pressure levels for output arrays: in
     &,ucomB_press,vcomB_press                                          &
     &,ucomp_press,vcomp_press,wcomp_press                              &
     &,testd_press                                                      &
     &,p_height,theta_height,rho_height                                 &
     &,w_height,u_height,v_height                                       &
     &,pv_press                                                         &
! Model levels for output arrays: in
     &,ucomB_model,vcomB_model                                          &
     &,testd_model                                                      &
     &,Htheta_model,Hrho_model                                          &
! Flags to request each diagnostic output field: in
! wind related diagnostics
     &,qucomB_m,qvcomB_m                                                &
     &,qucomB_p,qvcomB_p                                                &
     &,qucomp_p,qvcomp_p,qwcomp_p                                       &
     &,qu50mB_h,qv50mB_h                                                &
     &,qu50m_h,qv50m_h                                                  &
! PV related diagnostics
     &,qpotn_vort_theta,qtheta_potn_vort,qtheta_pv_points               &
     &,qpv_mod_levs,qpv_theta_mlev,qpotn_vort_press                     &
! test fields
     &,qdia1,qdia2,qdia3,qdia4                                          &
! flux diagnostics
     &,qrhow,qrhouw,qrhovw,qrhow_up                                     &
     &,qrhow_down,qrhowc_up,qrhowc_down                                 &
! height and height level diagnostics
     &,qHtheta_ml,qHrho_ml                                              &
     &,qpress_h,qtheta_h,qrho_h                                         &
     &,qu_h,qv_h,qw_h                                                   &
! other diagnostics
     &,spec_w,qtrue_density                                             &
! Flags for wind rotation (lam grid): in
     &,rot_uvcomB_p                                                     &
! Diagnostics lengths: in
     &,ucomB_m_levs,vcomB_m_levs                                        &
     &,ucomB_p_levs,vcomB_p_levs                                        &
     &,ucomp_p_levs,vcomp_p_levs,wcomp_p_levs                           &
     &,pv_theta_levs,pv_press_levs                                      &
     &,testd_p_levs,testd_m_levs                                        &
     &,Htheta_m_levs,Hrho_m_levs                                        &
     &,p_h_levs,theta_h_levs,rho_h_levs,w_h_levs,u_h_levs,v_h_levs      &
! Diagnostic arrays: out
! wind related diagnostics
     &,ucomB_m,vcomB_m                                                  &
     &,ucomB_p,vcomB_p                                                  &
     &,ucomp_p,vcomp_p,wcomp_p                                          &
     &,u50mB_h,v50mB_h                                                  &
     &,u50m_h,v50m_h                                                    &
! PV related diagnostics
     &,potn_vort_theta,theta_potn_vort,theta_pv_points                  &
     &,pv_mod_levs,pv_theta_mlev,potn_vort_press                        &
! test fields
     &,testdiag1,testdiag2,testdiag3,testdiag4                          &
! flux diagnostics
     &,rhow,rhouw,rhovw,rhow_up                                         &
     &,rhow_down,rhow_convup,rhow_convdown                              &
! height and height level diagnostics
     &,height_theta_ml,height_rho_ml                                    &
     &,press_h,theta_h,rho_h                                            &
     &,ucomp_h,vcomp_h,wcomp_h                                          &
! other diagnostics
     &,spec_3D,true_density)

      USE dynamics_grid_mod, ONLY: l_vatpoles

      USE atm_fields_bounds_mod, ONLY:                                  &
          udims, vdims, wdims, udims_s, vdims_s, wdims_s,               &
          tdims_s, pdims_s, udims_l, vdims_l

      USE atmos_constants_mod, ONLY: kappa, p_zero, recip_kappa

      USE swapable_field_mod, ONLY:                                     &
          swapable_field_pointer_type

      USE earth_constants_mod, ONLY: g, earth_radius

      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE interpor_mod, ONLY: interp_order_linear, interp_order_cubic
      USE p_to_u_mod, ONLY: p_to_u
      USE p_to_v_mod, ONLY: p_to_v
      USE uc_to_ub_mod, ONLY: uc_to_ub
      USE vert_interp2_mod, ONLY: vert_interp2
      USE vert_interp_mdi_mod, ONLY: vert_interp_mdi
      USE vert_interp_mdi2_mod, ONLY: vert_interp_mdi2
      USE w_eqtoll_mod, ONLY: w_eqtoll
      USE vc_to_vb_mod, ONLY: vc_to_vb
      IMPLICIT NONE
!
! Description:
!   Calculate dynamics-related diagnostics - held in STASH section 15 -
!   which may include interpolation onto pressure surfaces. Diagnostics
!   currently supported:
!   [All diagnostics on native 'C' grid unless specified.]
!   STASH item
!     2 u component of wind on model levels      on 'B' grid
!     3 v component of wind on model levels      on 'B' grid
!   201 u component of wind on pressure surfaces on 'B' grid
!   202 v component of wind on pressure surfaces on 'B' grid
!   243 u component of wind on pressure surfaces
!   244 v component of wind on pressure surfaces
!   242 w component of wind on pressure surfaces
!   245 u component of wind at 50m height
!   246 v component of wind at 50m height
!   212 u component of wind at 50m height on 'B' grid
!   213 v component of wind at 50m height on 'B' grid
!   214 potential vorticity on theta levels
!   229 potential vorticity on pressure levels
!   215 theta on potential vorticity = +/-2 surface
!   216 theta at potential vorticity points
!   217 potential vorticity on model levels
!   218 potential vorticity on model theta grid and theta levels
!   231 test analytic field on v grid - single level
!   232 test analytic field on p grid - single level
!   233 test analytic field on p grid - pressure levels
!   234 test analytic field on p grid - model levels
!   260 mass flux (rhow) on model levels
!   261 momentum flux (rhouw) on model levels
!   262 momentum flux (rhovw) on model levels
!   263 upward mass flux (rhow, w> 0m/s) on model levels
!   264 downward mass flux (rhow, w<0m/s) on model levels
!   265 upward convective mass flux (rhow, w>1m/s) on model levels
!   266 downward convective mass flux (rhow, w<-1m/s) on model levels
!   101 Height from sea level on theta model levels
!   102 Height from sea level on rho model levels
!   108 Pressure on height (from sea) levels
!   119 Potential temperature on height (from sea) levels
!   127 Rho (density) on height (from sea) levels
!   142 W component of wind on height (from sea) levels
!   143 U component of wind on height (from sea) levels
!   144 V component of wind on height (from sea) levels
!   270 real part of spectra (model levels)
!   271 true unscaled density on model levels
!
! Method:
!   Required level lists and logical switches are determined by the
!   calling routine from STASH requests and STASHflags.
!   Primary model data is input, and each diagnostic is calculated - if
!   its flag is set - in simple sequential order. Where the extraction
!   of the diagnostic quantity requires further calculation, a lower
!   level diagnostic-specific routine is called.
!
! STASH items 002,003     : u,v wind components on model levs 'B' grid:
!    Simple horizontal interpolation from u,v to uv 'B' position.
!
! STASH items 201,202     : u,v wind components on p surfaces 'B' grid:
! 1. Simple horizontal interpolation of p onto u,v staggering.
! 2. Perform vertical interpolation of u,v onto each o/p p surface.
! 3. Simple horizontal interpolation from u,v to uv 'B' position.
! 4. (Lam only) Rotate winds from native eq lat-long to standard grid
!
! STASH item 214: pv on theta levels:
! 1. Interpolate theta onto PV points.
! 2. Calculate PV on model levels.
! 3. Perform vertical interpolation of PV onto theta levels.

! STASH items 243,244,242 : u,v,w wind components on pressure surfaces:
! 1. Simple horizontal interpolation of p onto u,v staggering.
! 2. Perform vertical interpolation of u,v onto each o/p p surface.
!
! STASH items 212,213 : u,v wind components at 50m height 'B' grid:
! STASH items 244,245 : u,v wind components at 50m height:
! 1. Re-use p_at_u/v arrays to get height surface 50m above orography
! 2. Perform vertical interpolation of u,v onto height surface
! 3. Simple horizontal interpolation from u,v to uv 'B' position.
!
! STASH items 231,232,233,234: Test diagnostics 1-4:
! Call TestDiag routine. See UM DOc Paper D7.
!
! STASH item 215: theta on pv=+/-2 surface:
! 1. Interpolate theta onto PV points.
! 2. Calculate PV on model levels.
! 3. Set mod_PV = |PV|
! 4. Perform vertical interpolation of theta onto PV=2 surface
!
! STASH item 216: theta at pv points:
! 1. Interpolate theta onto PV points.
!
! STASH item 217: pv on model levels:
! 1. Calculate PV on model levels.
!
! STASH item 218: pv on model theta levels:
! 1. Calculate PV on model levels, using 'calc_pv_at_theta'.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Diagnostics
!
! Code Description:
!   Language: FORTRAN 90
!
! Declarations:
! Global variables (*CALLed COMDECKs etc...):
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
! Model_domain meaningful names

! Subroutine arguments
!   Scalar arguments with intent(in):
! Grid sizes:
      INTEGER                                                           &
     & rows,n_rows,row_length,model_levels,wet_levels,bl_levels         &
     &, global_rows,global_row_length                                   &
     &,theta_field_size,u_field_size,v_field_size                       &
     &,Model_domain
! Grid coordinates: in
      REAL                                                              &
     & delta_lambda,delta_phi                                           &
     &,ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole
! Time information: in
      REAL                                                              &
     & forecast_hrs       ! T+forecast_hrs: hours after analysis time
                          ! UM6.5 - MODEL_ANALYSIS_HRS changed to REAL -
                          ! requires FORECAST_HRS changed to REAL also  
! Diagnostics lengths: in
      INTEGER                                                           &
     & ucomB_m_levs                                                     & 
                          ! No of levels for output u on m 'B' grid
     &,vcomB_m_levs                                                     & 
                          ! No of levels for output v on m 'B' grid
     &,ucomB_p_levs                                                     & 
                          ! No of levels for output u on p 'B' grid
     &,vcomB_p_levs                                                     & 
                          ! No of levels for output v on p 'B' grid
     &,ucomp_p_levs                                                     &
                          ! No of levels for output of u on p
     &,vcomp_p_levs                                                     &
                          ! No of levels for output of v on p
     &,wcomp_p_levs                                                     &
                          ! No of levels for output of w on p
     &, pv_theta_levs                                                   &
                         !No of levels for output of pv on theta
     &, pv_press_levs                                                   &
                            !No of levels for output of pv on pressure
     &,testd_p_levs                                                     &
                          ! No of levels for output of testdiag3
     &,testd_m_levs                                                     &
                          ! No of levels for output of testdiag4
     &,Htheta_m_levs                                                    &
                              ! Num of levs for output H on theta lev
     &,Hrho_m_levs                                                      &
                              ! Num of levs for output H on rho lev
     &,p_h_levs                                                         &
                              ! Num of levs for output p on H levs
     &,theta_h_levs                                                     &
                              ! Num of levs for output theta on H levs
     &,rho_h_levs                                                       &
                              ! Num of levs for output rho on H levs
     &,w_h_levs                                                         &
                              ! Num of levs for output w on H levs
     &,u_h_levs                                                         &
                              ! Num of levs for output u on H levs
     &,v_h_levs               ! Num of levs for output v on H levs
! Flags to request each diagnostic output field: IN
      LOGICAL                                                           &
     & qucomB_m                                                         & 
                     ! Flag for U wind on model    levels  'B' grid
     &,qvcomB_m                                                         & 
                     ! Flag for V wind on model    levels  'B' grid
     &,qucomB_p                                                         & 
                     ! Flag for U wind on pressure levels  'B' grid
     &,qvcomB_p                                                         & 
                     ! Flag for V wind on pressure levels  'B' grid
     &,qucomp_p                                                         &
                     ! Flag for U wind component on pressure levels
     &,qvcomp_p                                                         &
                     ! Flag for V wind component on pressure levels
     &,qwcomp_p                                                         &
                     ! Flag for W wind component on pressure levels
     &,qu50mB_h                                                         & 
                     ! Flag for U wind comp at 50m height  'B' grid
     &,qv50mB_h                                                         & 
                     ! Flag for V wind comp at 50m height  'B' grid
     &,qu50m_h                                                          &
                     ! Flag for U wind component at 50m height
     &,qv50m_h                                                          &
                     ! Flag for V wind component at 50m height
     &, qpotn_vort_theta                                                &
                            !Flag for pv on theta levels
     &, qpotn_vort_press                                                &
                            !Flag for pv on pressure levels
     &,qtheta_potn_vort                                                 &
                           !Flag for theta on pv=+/-2 surface
     &,qtheta_pv_points                                                 &
                           !Flag for theta at pv points
     &,qpv_mod_levs                                                     &
                           !Flag for pv on model levels
     &,qpv_theta_mlev                                                   &
                           !Flag for pv on model theta points and levels
     &,qdia1                                                            &
                     ! Flag for test diagnostic 1
     &,qdia2                                                            &
                     ! Flag for test diagnostic 2
     &,qdia3                                                            &
                     ! Flag for test diagnostic 3
     &,qdia4                                                            &
                     ! Flag for test diagnostic 4
     &,qrhow                                                            &
                     ! Flag for rhow
     &,qrhouw                                                           &
                     ! Flag for rhouw
     &,qrhovw                                                           &
                     ! Flag for rhovw
     &,qrhow_up                                                         &
                     ! Flag for rhow_up
     &,qrhow_down                                                       &
                     ! Flag for rhow_down
     &,qrhowc_up                                                        &
                     ! Flag for rhow_convup
     &,qrhowc_down                                                      &
                     ! Flag for rhow_convdown
     &,qHtheta_ml                                                       &
                                 ! Flag for height on theta levs
     &,qHrho_ml                                                         &
                                 ! Flag for height on rho levs
     &,qpress_H                                                         &
                               ! Flag for press on height levs
     &,qtheta_H                                                         &
                               ! Flag for theta on height levs
     &,qrho_H                                                           &
                               ! Flag for rho on height levs
     &,qw_H                                                             &
                               ! Flag for w on height levs
     &,qu_H                                                             &
                               ! Flag for u on height levs
     &,qv_H                                                             &
                               ! Flag for v on height levs
     &,spec_w                                                           &
                               ! Flag for real part of spectra          
     &,qtrue_density           ! flag for true unscaled density
! Flags for wind rotation (lam grid)  (rotate if .T. ) : in
      LOGICAL, INTENT(IN) ::                                            &
     & rot_uvcomB_p          ! u,v (B grid) on p levels
!   Array  arguments with intent(in):
! Primary data: IN
      REAL, TARGET ::                                                   &
         u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end, &
           udims_s%k_start:udims_s%k_end)                               &
        ,v(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end, &
           vdims_s%k_start:vdims_s%k_end)                               &
        ,rho(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end,                             &
             pdims_s%k_start:pdims_s%k_end)                             &    
        ,exner_rho_levels(pdims_s%i_start:pdims_s%i_end,                &
                          pdims_s%j_start:pdims_s%j_end,                &
                          pdims_s%k_start:pdims_s%k_end)                &
        ,theta(tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end)                           
     
      REAL                                                              &
         w(wdims_s%i_start:wdims_s%i_end,wdims_s%j_start:wdims_s%j_end, &
           wdims_s%k_start:wdims_s%k_end)                               &      
        ,exner_theta_levels(tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end)              &
! Vertical grid definition: IN
      ,eta_theta_levels(0:model_levels)                                 &
                                        ! vertical grid for theta vars
      ,eta_rho_levels    (model_levels)                                 &
                                        ! vertical grid for rho   vars
! Pre-calculated grid associated arrays: IN
      ,r_at_u (udims_l%i_start:udims_l%i_end,udims_l%j_start:udims_l%j_end, &
           udims_l%k_start:udims_l%k_end)                               &
      ,r_at_v (vdims_l%i_start:vdims_l%i_end,vdims_l%j_start:vdims_l%j_end, &
           vdims_l%k_start:vdims_l%k_end)                               &
      , r_theta_levels(1-halo_i:row_length+halo_i,                      &
                         1-halo_j:rows+halo_j,0:model_levels)           &
      , r_rho_levels(1-halo_i:row_length+halo_i,                        &
                       1-halo_j:rows+halo_j, model_levels)              &
      , sec_v_latitude (vdims_s%i_start:vdims_s%i_end,                  &
                        vdims_s%j_start:vdims_s%j_end)                  &
      , tan_v_latitude(vdims%i_start:vdims%i_end,                       &
                       vdims%j_start:vdims%j_end)                       &
      , sec_theta_latitude (1-offx:row_length+offx,1-offy:rows+offy)    &
      , f3_at_v (vdims_s%i_start:vdims_s%i_end,                         &
                 vdims_s%j_start:vdims_s%j_end)                         &
      ,rot_coeff1(udims%i_start:udims%i_end,                            &
                  vdims%j_start:vdims%j_end)                            &
                                      ! for lam wind rotations
      ,rot_coeff2(udims%i_start:udims%i_end,                            &
                  vdims%j_start:vdims%j_end)                            &
                                      !     (on B grid)
! Pressure levels (units mb) for output arrays: IN
     &,ucomB_press(ucomB_p_levs)                                        & 
                                        ! for u wind  'B' grid  press
     &,vcomB_press(vcomB_p_levs)                                        & 
                                        ! for v wind  'B' grid  press
     &,ucomp_press(ucomp_p_levs)                                        &
                                        ! for u wind
     &,vcomp_press(vcomp_p_levs)                                        &
                                        ! for v wind
     &,wcomp_press(wcomp_p_levs)                                        &
                                        ! for w wind
     &,testd_press(testd_p_levs)                                        &
                                        ! for test diagnostic3
     &,p_height(p_h_levs)                                               &
                                        ! for diagnostic 108 (press)
     &,theta_height(theta_h_levs)                                       &
                                        ! for diagnostic 119 (theta)
     &,rho_height(rho_h_levs)                                           &
                                        ! for diagnostic 127 (rho)
     &,w_height(w_h_levs)                                               &
                                        ! for diagnostic 142 (w)
     &,u_height(u_h_levs)                                               &
                                        ! for diagnostic 143 (u)
     &,v_height(v_h_levs)                                               &
                                        ! for diagnostic 144 (v)
     &, desired_theta(pv_theta_levs)                                    &
                                      ! for potential vorticity
                                      !  on theta levels
     &, pv_press(pv_press_levs)                                         &
                                      ! for potential vorticity
                                      !  on pressure levels
! Model    levels for output arrays: IN
     &,testd_model(testd_m_levs)        ! for test diagnostic4

      INTEGER                                                           &
     & ucomB_model(ucomB_m_levs)                                        & 
                                        ! for u wind  'B' grid  model
     &,vcomB_model(vcomB_m_levs)                                        & 
                                        ! for v wind  'B' grid  model
     &,Htheta_model(Htheta_m_levs)                                      &
                                        ! for diagnostic 101
     &,Hrho_model(Hrho_m_levs)          ! for diagnostic 102
!   Scalar arguments with intent(InOut):

!   Array  arguments with intent(InOut):

!   Scalar arguments with intent(out):

!   Array  arguments with intent(out):
! Diagnostic arrays: OUT
      REAL                                                                &
       ucomB_m(udims%i_start:udims%i_end,                                 &
               vdims%j_start:vdims%j_end,ucomB_m_levs)                    & 
                                               ! u 'B' grid at mod levs
      ,vcomB_m(udims%i_start:udims%i_end,                                 &
               vdims%j_start:vdims%j_end,vcomB_m_levs)                    & 
                                               ! v 'B' grid at mod levs
      ,ucomB_p(udims%i_start:udims%i_end,                                 &
               vdims%j_start:vdims%j_end,ucomB_p_levs)                    & 
                                               ! u 'B' grid at pressures
      ,vcomB_p(udims%i_start:udims%i_end,                                 &
               vdims%j_start:vdims%j_end,vcomB_p_levs)                    & 
                                               ! v 'B' grid at pressures
      ,ucomp_p(udims%i_start:udims%i_end,                                 &
               udims%j_start:udims%j_end,ucomp_p_levs)                    &
                                               ! u at selected pressures
      ,vcomp_p(vdims%i_start:vdims%i_end,                                 &
               vdims%j_start:vdims%j_end,ucomp_p_levs)                    &
                                               ! v at selected pressures
      ,wcomp_p(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,       &
                                         wcomp_p_levs)                    &
                                               ! w at selected pressures
      ,u50mB_h(udims%i_start:udims%i_end,                                 &
               vdims%j_start:vdims%j_end)                                 & 
                                               ! u at 50m ht 'B' grid
      ,v50mB_h(udims%i_start:udims%i_end,                                 &
               vdims%j_start:vdims%j_end)                                 & 
                                               ! v at 50m ht 'B' grid
      ,u50m_h(udims%i_start:udims%i_end,                                  &
              udims%j_start:udims%j_end)                                  &
                                               ! u at 50m height
      ,v50m_h(vdims%i_start:vdims%i_end,                                  &
              vdims%j_start:vdims%j_end)                                  &
                                               ! v at 50m height
      , potn_vort_theta(udims%i_start:udims%i_end,                        &
                        vdims%j_start:vdims%j_end, pv_theta_levs)         &
!                                  pv on theta levels
      , potn_vort_press(udims%i_start:udims%i_end,                        &
                        vdims%j_start:vdims%j_end, pv_press_levs)         &
                                               ! pv on pressure levels
      ,theta_potn_vort(udims%i_start:udims%i_end,                         &
                       vdims%j_start:vdims%j_end)                         &
                                               ! theta on pv+/-2
      ,theta_pv_points(udims%i_start:udims%i_end,                         &
                       vdims%j_start:vdims%j_end,model_levels)            &
                                                       ! theta at pv poi
      ,pv_mod_levs(udims%i_start:udims%i_end,                             &
                   vdims%j_start:vdims%j_end,model_levels)                &
                                                   ! pv on model levels
      ,pv_theta_mlev(row_length,rows,model_levels)                        &
                                              ! pv on model theta levels
     &,testdiag1(row_length,n_rows)                                     &
                                               ! testdiag1 -single level
     &,testdiag2(row_length,  rows)                                     &
                                               ! testdiag2 -single level
     &,testdiag3(row_length,rows,testd_p_levs)                          &
                                               ! testdiag3 -press levels
     &,testdiag4(row_length,rows,testd_m_levs)                          &
                                               ! testdiag4 -model levels
     &,rhow(row_length,rows,model_levels)                               &
                                   ! mass flux on model levels
     &,rhow_up(row_length,rows,model_levels)                            &
                                   ! up mass flux on model levels
     &,rhow_down(row_length,rows,model_levels)                          &
                                   ! down mass flux on model levels
     &,rhow_convup(row_length,rows,model_levels)                        &
                                   ! up conv mass flux on model levels
     &,rhow_convdown(row_length,rows,model_levels)                      &
                                   ! down conv mass flux on model levels
     &,rhouw(row_length,rows,model_levels)                              &
                                             ! momentum flux
     &,rhovw(row_length,rows,model_levels)                              &
                                             ! momentum flux
! stash work arrays for height on theta and rho model levels
      ,height_theta_ml(row_length,rows,Htheta_m_levs)                   &
      ,height_rho_ml(row_length,rows,Hrho_m_levs)                       &
! stash work arrays for height level diagnostics
      ,press_h(row_length,rows,p_h_levs)                                &
      ,theta_h(row_length,rows,theta_h_levs)                            &
      ,rho_h(row_length,rows,rho_h_levs)                                &
      ,wcomp_h(row_length,rows,w_h_levs)                                &
      ,ucomp_h(udims%i_start:udims%i_end,                               &
               udims%j_start:udims%j_end,u_h_levs)                      &
      ,vcomp_h(vdims%i_start:vdims%i_end,                               &
               vdims%j_start:vdims%j_end,v_h_levs)                      &
      ,true_density(row_length,rows,model_levels) ! unscaled density


! Local parameters:
      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='Dyn_diag')
      REAL                                                              &
     & z_50m                    ! height=50m above earth's surface
      PARAMETER (                                                       &
     & z_50m = 50.                                                      &
     &)

! Local scalars:
!   ErrorStatus
      INTEGER      ErrorStatus                                          &
                                        ! Error flag (0 = OK)
     &,i,j,k                                                            &
                                        ! loop counters
     &,kk                                                               &
                                        ! k level index
     &,interp_order    !  order of vertical interpolation

      INTEGER :: i_field  ! counter for swap_bounds_mv

      CHARACTER(LEN=256)                                                     &
     & CMessage         ! Error message if return code >0

      LOGICAL :: LAM              ! T: limited area model
      LOGICAL :: GLOBAL           ! T: global model

      REAL                                                              &
     & dummy                                                            &
                                ! dummy argument - not referenced
     &,desired_potn_vort                                                &
                                ! value of pv surface for theta
     &,pressure_pa                                                      &
                                ! pressure in pascals
     &,pressure_ex                                                      &
                                ! exner pressure
     &,weight                                                           &
                           ! weight for calculating w_on_rho
     &,desired_r              ! height to interpolate

! Local dynamic arrays:
      REAL                                                              &
       exner_at_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,  &
                    udims%k_start:udims%k_end),                         &
                                                    ! exner at u points
       exner_at_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,  &
                    vdims%k_start:vdims%k_end)                          &
                                                    ! exner at v points      
      ,exner_at_w(wdims%i_start:wdims%i_end,                      &
                  wdims%j_start:wdims%j_end,1:wdims%k_end)              &
                                                  !pressure at w points
      ,exner_at_pv(udims%i_start:udims%i_end,                           &
                   vdims%j_start:vdims%j_end, model_levels)             &
                                               !pressure at pv points
![Note: exner_at_u,exner_at_v are optionally re-used as height points.]
      ,theta_at_PV(udims%i_start:udims%i_end,                           &
                   vdims%j_start:vdims%j_end, model_levels)             &
      ,T (1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
      ,PV(udims%i_start:udims%i_end,                                    &
          vdims%j_start:vdims%j_end, model_levels)                      &
      ,mod_PV(udims%i_start:udims%i_end,                                &
              vdims%j_start:vdims%j_end, model_levels)                  &
                                                 ! for pv=-2 diag
      ,w_on_rho(row_length,rows,model_levels)                           &
                                          ! vertical velocity on rho pts
      ,u_on_rho(row_length,rows,model_levels)                           &
                                              ! u velocity on rho pts
      ,v_on_rho(row_length,rows,model_levels)                           &
                                              ! v velocity on rho pts
      ,true_rho(row_length,rows,model_levels) ! rho/(r_rho_levels^2)

     REAL :: work_u(udims%i_start:udims%i_end,                          &
                       udims%j_start:udims%j_end)  
                            ! workspace, used for u & v
             
     REAL :: work_v(vdims%i_start:vdims%i_end,                          &
                       vdims%j_start:vdims%j_end)  
                            ! workspace, used for u & v


     REAL :: work_1(udims%i_start:udims%i_end,                          &
              vdims%j_start:vdims%j_end)  ! workspace, used for u & v
     REAL :: work_2(udims%i_start:udims%i_end,                          &
              vdims%j_start:vdims%j_end)  ! workspace, used for u & v


! Spectral variables

      Real                                                              &
     & spec_3D(row_length,rows,model_levels)                            &
     &,spec_2D(1-offx:row_length+offx,1-offy:rows+offy)

      Integer                                                           &
     & klev                                                             &
     &,gath_proc

      Type(swapable_field_pointer_type) :: fields_to_swap(5) ! mv swapbounds

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!- End of header


! ----------------------------------------------------------------------
! Section 0.  Initialisation.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('DYN_DIAG',zhook_in,zhook_handle)

! Set Error code to zero
      ErrorStatus = 0

! Set order of vertical interpolation
      interp_order = interp_order_linear

! Determine whether limited area model
      IF (model_domain  ==  mt_lam .OR.                                 &
     &    model_domain  ==  mt_cyclic_lam .or.                          &
     &    model_domain  ==  mt_bi_cyclic_lam) THEN
        lam = .true.
      ELSE
        lam = .false.
      END IF

      IF (model_domain  ==  mt_global) THEN
        global = .TRUE.
      END IF

!-----------------------------------------------
!    do all swapbounds calls needed
!-----------------------------------------------
      i_field = 0
      IF(qucomp_p .OR. qucomB_p .OR. qpotn_vort_press) THEN
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => exner_rho_levels(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.
      END IF

      IF(qpotn_vort_theta .OR. qtheta_potn_vort .OR.                    &
     &   qtheta_pv_points .OR. qpotn_vort_press .OR. qpv_mod_levs .OR.  &
     &   qpv_theta_mlev) THEN  
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => theta(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  tdims_s%k_len
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.
      END IF

      IF(qucomB_m .OR. qpotn_vort_press .OR. qtheta_potn_vort .OR.      &
     &   qpv_mod_levs .OR. qrhouw .OR. qpotn_vort_theta .OR.            &
     &   qpv_theta_mlev) THEN  
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => u(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_u
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .TRUE.
      END IF

      IF(qvcomB_m .OR. qpotn_vort_press .OR. qtheta_potn_vort .OR.      &
     &   qpv_mod_levs .OR. qrhovw .OR.qpotn_vort_theta .OR.             &
     &   qpv_theta_mlev) THEN  
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => v(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_v
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  n_rows
        fields_to_swap(i_field) % vector      =  .TRUE.
      END IF

      IF (qpotn_vort_press .OR. qtheta_potn_vort .OR.                   &
     &    qpv_mod_levs .OR. qpotn_vort_theta .OR. qpv_theta_mlev) THEN
        i_field = i_field + 1
        fields_to_swap(i_field) % field       => rho(:,:,:)
        fields_to_swap(i_field) % field_type  =  fld_type_p
        fields_to_swap(i_field) % levels      =  model_levels
        fields_to_swap(i_field) % rows        =  rows
        fields_to_swap(i_field) % vector      =  .FALSE.
      END IF

      IF (i_field > 0) THEN
! DEPENDS ON: swap_bounds_mv
        CALL swap_bounds_mv(fields_to_swap, i_field, row_length,        &
                            offx,           offy)
      END IF

      IF(qucomp_p.OR.qucomB_p) THEN
!  Calculate exner at u points. Store in exner_at_u

      CALL p_to_u(exner_rho_levels,                               &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   udims%k_start, udims%k_end,                    &
                   exner_at_u)

      END IF ! on STASHflag

      IF(qvcomp_p.OR.qvcomB_p) THEN
!  Calculate exner at v points. Store in exner_at_v

      CALL p_to_v(exner_rho_levels,                               &   
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   vdims%i_start,vdims%i_end,                     &
                   vdims%j_start,vdims%j_end,                     &
                   vdims%k_start, vdims%k_end,                    &
                   exner_at_v)  

      END IF ! on STASHflag

      IF(qwcomp_p) THEN
!  Calculate exner at w points. Store in exner_at_w

      DO k = 1, wdims%k_end
        DO j = wdims%j_start, wdims%j_end
          DO i = wdims%i_start, wdims%i_end
             exner_at_w(i,j,k) = exner_theta_levels(i,j,k)
           END DO ! i
        END DO ! j
      END DO ! k

      END IF ! on STASHflag

      IF(qpotn_vort_theta.OR.qtheta_potn_vort.OR.                       &
     &   qtheta_pv_points) THEN

! Calculate theta at PV points. Store in theta_at_pv
! first interpolate theta to rho levels. Use linear interpolation.
! Store in T as work space.

        DO j = 1-offy, rows+offy
          DO i = 1-offx, row_length+offx
             T(i,j,1) = theta(i,j,1)
          END DO
        END DO

        DO k = 2, model_levels
          DO j = 1-offy, rows+offy
            DO i = 1-offx, row_length+offx
                    T(i,j,k) = (theta(i,j,k) *                          &
     &                        (r_rho_levels(i,j,k) -                    &
     &                         r_theta_levels(i,j,k-1) ) +              &
     &                        theta(i,j,k-1) *                          &
     &                        (r_theta_levels(i,j,k) -                  &
     &                         r_rho_levels(i,j,k) ) ) /                &
     &                        (r_theta_levels(i,j,k) -                  &
     &                         r_theta_levels(i,j,k-1) )
            END DO
          END DO
        END DO

        DO k = 1, model_levels
          DO j = vdims%j_start,vdims%j_end
            DO i = udims%i_start,udims%i_end
               theta_at_pv(i,j,k) = .25*(T(i+1,j,k) + T(i,j,k) +        &
     &                                   T(i+1,j+1,k) + T(i,j+1,k) )
            END DO
          END DO
        END DO
      END IF ! on STASHflag

      IF(qpotn_vort_press) THEN

! Interpolate exner onto PV points.

        DO k = 1, model_levels
          DO j = vdims%j_start,vdims%j_end
            DO i = udims%i_start,udims%i_end
               exner_at_pv(i,j,k) = .25*(exner_rho_levels(i+1,j,k) +    &
     &                                   exner_rho_levels(i,j,k) +      &
     &                                   exner_rho_levels(i+1,j+1,k) +  &
     &                                   exner_rho_levels(i,j+1,k) )
            END DO
          END DO
        END DO
      END IF ! on STASHflag

! STASH items 002,003     : u,v wind components on model levs 'B' grid
! ----------------------------------------------------------------------

      IF(qucomB_m) THEN

        DO  k=1,ucomB_m_levs
! Perform simple horizontal interpolation from 'C' to 'B' grid
! Halos already populated, so interpolate directly:

          kk = ucomB_model(k) ! selected model level

           DO j=udims%j_start, udims%j_end
             DO i=udims%i_start, udims%i_end
               work_u(i,j) = u(i,j,kk)
             END DO ! i
           END DO ! j


          CALL  uC_to_uB(work_u,                                        &
            row_length,rows,n_rows,1,offx,offy,                         &
            ucomB_m(:,:,k))

        END DO  ! k model levels loop

      END IF ! on STASHflag

      IF(qvcomB_m) THEN

        IF (l_vatpoles) THEN

        IF (global .AND. .NOT. qucomB_m) THEN
          ! In the vatpoles case, v on B grid model levels is 
          !   dependent on u on B grid model levels
          ErrorStatus = -1        ! Warning
          Cmessage='v on B grid model levels not possible when '   &
               //'u on B grid model levels not active. Diagnostic aborted'

          CALL Ereport(Routinename,ErrorStatus,Cmessage)

        ELSE

        DO  k=1,vcomB_m_levs

! Perform simple horizontal interpolation from 'C' to 'B' grid

! Perform simple horizontal interpolation from 'C' to 'B' grid
! Halos already populated, so interpolate directly:

          kk = vcomB_model(k) ! selected model level

           DO j=vdims%j_start, vdims%j_end
             DO i=vdims%i_start, vdims%i_end
               work_v(i,j) = v(i,j,kk)
             END DO ! i
           END DO ! j



            ! In the vatpoles case, u on the B grid is required by v 
            ! on the B grid to compute polar values, so the model 
            ! levels must match.
            IF(ucomB_model(k) == vcomB_model(k)) THEN

          CALL  vC_to_vB(work_v,                                            &
            rows,row_length,n_rows,1,offx,offy,                             &
            global_row_length,                                              &
            vcomB_m(udims%i_start,vdims%j_start,k)                          &
            ,ucomB_m(udims%i_start,vdims%j_start,k)                         &
           )    

            ELSE
              ErrorStatus = -1        ! Warning
              Cmessage='u and v on B grid mode;s levels - '   &
               //'requested pressure levels do not match, diagnostic aborted'

              CALL Ereport(Routinename,ErrorStatus,Cmessage)
              EXIT                    ! jump out of levels loop            
            END IF


        END DO  ! k model levels loop

        END IF  ! global .AND. .NOT. qucomB_m

        ELSE  ! vatpoles

        DO  k=1,vcomB_m_levs

! Perform simple horizontal interpolation from 'C' to 'B' grid
! Halos already populated, so interpolate directly:

          kk = vcomB_model(k) ! selected model level

           DO j=vdims%j_start, vdims%j_end
             DO i=vdims%i_start, vdims%i_end
               work_v(i,j) = v(i,j,kk)
             END DO ! i
           END DO ! j


          CALL  vC_to_vB(work_v,                                            &
            rows,row_length,n_rows,1,offx,offy,                             &
            global_row_length,                                              &
            vcomB_m(udims%i_start,vdims%j_start,k)                          &
           )


        END DO  ! k model levels loop

        END IF  ! vatpoles

      END IF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 201,202     : u,v wind components on p surfaces 'B' grid
! ----------------------------------------------------------------------

      IF(qucomB_p) THEN

        DO  k=1,ucomB_p_levs

          pressure_pa = ucomB_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
          CALL vert_interp2 (u, row_length, rows, model_levels          &
                                ,pressure_ex                            &
                                ,offx, offy, 0, 0                       &
                                ,exner_at_u, interp_order               &
                                ,work_u )

! Perform simple horizontal interpolation from 'C' to 'B' grid

          CALL  uC_to_uB(work_u,                                        &
            row_length,rows,n_rows,1,offx,offy,                         &
            ucomB_p(:,:,k))

        END DO  ! k pressure levels loop

      END IF ! on STASHflag

      IF(qvcomB_p) THEN

        IF (l_vatpoles) THEN

        IF (global .AND. .NOT. qucomB_p) THEN
          ! In the vatpoles case, v on B grid pressure levels is 
          !   dependent on u on B grid pressure levels
          ErrorStatus = -1        ! Warning
          Cmessage='v on B grid pressure levels not possible when '   &
               //'u on B grid pressure levels not active. Diagnostic aborted'

          CALL Ereport(Routinename,ErrorStatus,Cmessage)

        ELSE
          DO  k=1,vcomB_p_levs

            pressure_pa = vcomB_press(k)*100.0   ! convert to Pascals
            pressure_ex = ( pressure_pa /p_zero )**kappa
            CALL vert_interp2 (v, row_length, n_rows, model_levels      &
                                ,pressure_ex                            &
                                ,offx, offy, 0, 0                       &
                                ,exner_at_v, interp_order               &
                                ,work_v )


! Perform simple horizontal interpolation from 'C' to 'B' grid

            ! In the vatpoles case, u on the B grid is required by v 
            ! on the B grid to compute polar values, so the pressure 
            ! levels must match.
            IF(ucomB_press(k) == vcomB_press(k)) THEN

          CALL  vC_to_vB(work_v,                                        &
            rows,row_length,n_rows,1,offx,offy,                         &
            global_row_length,                                          &
            vcomB_p(udims%i_start,vdims%j_start,k)                          &
            ,ucomB_p(udims%i_start,vdims%j_start,k)                         &
           )    


            ELSE
              ErrorStatus = -1        ! Warning
              Cmessage='u and v on B grid pressure levels - '   &
               //'requested pressure levels do not match, diagnostic aborted'

              CALL Ereport(Routinename,ErrorStatus,Cmessage)
              EXIT                    ! jump out of levels loop            
            END IF

! Rotate winds from model to standard lat-long grid
            IF (lam .AND. rot_uvcomB_p) THEN
! First check valid requests: implicit assumption that u components
! and v components are requested for the same pressure levels
              IF(qucomB_p .AND. ucomB_press(k) == vcomB_press(k)) THEN

                DO j=vdims%j_start,vdims%j_end
                  DO i=udims%i_start,udims%i_end
                    work_1(i,j) = ucomB_p(i,j,k)
                    work_2(i,j) = vcomB_p(i,j,k)
                  END DO
                END DO

! Rotation calculation on B grid
                CALL W_EqtoLL(rot_coeff1,rot_coeff2,work_1,work_2,        &
                 ucomB_p(udims%i_start,vdims%j_start,k),                  &
                 vcomB_p(udims%i_start,vdims%j_start,k),v_field_size,.TRUE.)

              ELSE

                ErrorStatus = -1        ! Warning
                Cmessage='wind diagnostics cannot be rotated: u and v '   &
                   //'requested components on pressure levels must match'

                CALL Ereport(Routinename,ErrorStatus,Cmessage)
                EXIT                    ! jump out of levels loop

              END IF       ! Check valid request

            END IF  ! (lam wind rotation)

          END DO  ! k pressure levels loop

        END IF  ! global .AND. .NOT. qucomB_p

        ELSE  ! vatpoles

          DO  k=1,vcomB_p_levs

            pressure_pa = vcomB_press(k)*100.0   ! convert to Pascals
            pressure_ex = ( pressure_pa /p_zero )**kappa
            CALL vert_interp2 (v, row_length, n_rows, model_levels      &
                                ,pressure_ex                            &
                                ,offx, offy, 0, 0                       &
                                ,exner_at_v, interp_order               &
                                ,work_v )


! Perform simple horizontal interpolation from 'C' to 'B' grid

          CALL  vC_to_vB(work_v,                                        &
            rows,row_length,n_rows,1,offx,offy,                         &
            global_row_length,                                          &
            vcomB_p(udims%i_start,vdims%j_start,k)                          &
           )    


! Rotate winds from model to standard lat-long grid
            IF (lam .AND. rot_uvcomB_p) THEN
! First check valid requests: implicit assumption that u components
! and v components are requested for the same pressure levels
              IF(qucomB_p .AND. ucomB_press(k) == vcomB_press(k)) THEN

                DO j=vdims%j_start,vdims%j_end
                  DO i=udims%i_start,udims%i_end
                    work_1(i,j) = ucomB_p(i,j,k)
                    work_2(i,j) = vcomB_p(i,j,k)
                  END DO
                END DO

! Rotation calculation on B grid
                CALL W_EqtoLL(rot_coeff1,rot_coeff2,work_1,work_2,        &
                 ucomB_p(udims%i_start,vdims%j_start,k),                  &
                 vcomB_p(udims%i_start,vdims%j_start,k),v_field_size,.TRUE.)

              ELSE

                ErrorStatus = -1        ! Warning
                Cmessage='wind diagnostics cannot be rotated: u and v '   &
                   //'requested components on pressure levels must match'

                CALL Ereport(Routinename,ErrorStatus,Cmessage)
                EXIT                    ! jump out of levels loop

              END IF       ! Check valid request

            END IF  ! (lam wind rotation)

          END DO  ! k pressure levels loop

        END IF  ! vatpoles

      END IF ! on STASHflag

! ----------------------------------------------------------------------
! Calculate PV for use with STASH items 229,214,215,217
! ----------------------------------------------------------------------

      IF (qpotn_vort_press .OR. qtheta_potn_vort .OR.                   &
     &    qpv_mod_levs .or. qpotn_vort_theta) THEN

! Calculate PV on model levels.
 
! DEPENDS ON: calc_pv
        Call Calc_PV                                                    &
     &            (u, v, theta, rho,                                    &
     &             r_theta_levels, r_rho_levels,                        &
     &             r_at_u, r_at_v,                                      &
     &             sec_v_latitude, tan_v_latitude,                      &
     &             sec_theta_latitude, f3_at_v,                         &
     &             delta_lambda, delta_phi,                             &
     &             row_length, rows, n_rows, model_levels,              &
     &             offx, offy, halo_i, halo_j,                          &
     &             at_extremity,                                        &
     &             PV)

      END IF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 229 : potential vorticity on pressure levels
! ----------------------------------------------------------------------

      IF (qpotn_vort_press) THEN

        DO k = 1, pv_press_levs

          pressure_pa = pv_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
          Call vert_interp2 (pv, row_length, n_rows,                    &
                             model_levels,                              &
                             pressure_ex,                               &
                             0, 0, 0, 0,                                &
                             exner_at_pv, interp_order_cubic,           &
                             potn_vort_press(:,:,k) )


        END DO  ! k pressure levels loop

      END IF ! on STASHflag

! ----------------------------------------------------------------------
! STASH item 214 : potential vorticity on theta levels
! ----------------------------------------------------------------------

      IF (qpotn_vort_theta) THEN

        DO k = 1, pv_theta_levs

          Call vert_interp_mdi2 (PV, row_length, n_rows,                &
                                 model_levels,                          &
                                 desired_theta(k),                      &
                                 0, 0, 0, 0,                            &
                                 theta_at_pv, interp_order_cubic,       &
                                 rmdi,                                  &
                                 potn_vort_theta(:,:,k) )

        END DO  ! k theta levels loop

      END IF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 243,244,242 : u,v,w wind components on pressure surfaces
! ----------------------------------------------------------------------

      IF(qucomp_p) THEN

        DO  k=1,ucomp_p_levs

          pressure_pa = ucomp_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
          CALL vert_interp2 (u, row_length, rows, model_levels          &
                                ,pressure_ex                            &
                                ,offx, offy, 0, 0                       &
                                ,exner_at_u, interp_order               &
                                ,ucomp_p(:,:,k) )

        END DO  ! k pressure levels loop

      END IF ! on STASHflag

      IF(qvcomp_p) THEN
        DO  k=1,vcomp_p_levs

          pressure_pa = vcomp_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
          CALL vert_interp2 (v, row_length, n_rows, model_levels        &
                                ,pressure_ex                            &
                                ,offx, offy, 0, 0                       &
                                ,exner_at_v, interp_order               &
                                ,vcomp_p(:,:,k) )

        END DO  ! k pressure levels loop

      END IF ! on STASHflag

      IF(qwcomp_p) THEN
        DO  k=1,wcomp_p_levs

          pressure_pa = wcomp_press(k)*100.0   ! convert to Pascals
          pressure_ex = ( pressure_pa /p_zero )**kappa
          CALL vert_interp2 (w(:,:,1:), row_length, rows, model_levels  &
                                ,pressure_ex                            &
                                ,offx, offy, 0, 0                       &
                                ,exner_at_w, interp_order               &
                                ,wcomp_p(:,:,k) )


        END DO  ! k pressure levels loop

      END IF ! on STASHflag

! ----------------------------------------------------------------------
! STASH items 212,213 : u,v wind components at 50m height 'B' grid
! STASH items 245,246 : u,v wind components at 50m height
! ----------------------------------------------------------------------

! Restrict interpolation to boundary levels, which should always
! greatly exceed levels in the vicinity of 50m.


      IF(qu50mB_h.OR.qu50m_h) THEN
! Generate height field above orography for lower (ie boundary) levels
! at u pts: (re-use array exner_at_u for workspace)
        DO k=1,bl_levels
          DO j=udims%j_start, udims%j_end
            DO i=udims%i_start, udims%i_end
              exner_at_u(i,j,k)= r_at_u(i,j,k) - r_at_u(i,j,1)
            END DO ! i
          END DO ! j
        END DO ! k

         CALL vert_interp_mdi (u, row_length, rows,                     &
                               bl_levels, z_50m,                        &
                               offx, offy,                              &
                               0, 0,                                    &
                               exner_at_u, interp_order,                &
                               rmdi, work_u )
         IF(qu50m_h) THEN
           DO j=udims%j_start, udims%j_end
             DO i=udims%i_start, udims%i_end
               u50m_h(i,j) = work_u(i,j)
             END DO ! i
           END DO ! j

         END IF ! on STASHflag qu50m_h

       IF(qu50mB_h) THEN
! Perform simple horizontal interpolation from 'C' to 'B' grid

          CALL  uC_to_uB(work_u,                                        &
            row_length,rows,n_rows,1,offx,offy,                         &
            u50mB_h(udims%i_start,vdims%j_start))

       END IF ! on STASHflag qu50mB_h

      END IF ! on STASHflags  qu50mB_h.OR.qu50m_h

      IF(qv50mB_h.OR.qv50m_h) THEN

! Generate height field above orography for lower (ie boundary) levels
! at v pts: (re-use array exner_at_v for workspace)
        DO k=1,bl_levels
          DO j=vdims%j_start, vdims%j_end
            DO i=vdims%i_start, vdims%i_end
              exner_at_v(i,j,k)= r_at_v(i,j,k) - r_at_v(i,j,1)
            END DO ! i
          END DO ! j
        END DO ! k

         CALL vert_interp_mdi (v, row_length, n_rows,                   &
                               bl_levels, z_50m,                        &
                               offx, offy,                              &
                               0, 0,                                    &
                               exner_at_v, interp_order,                &
                               rmdi, work_v)
         IF(qv50m_h) THEN
           DO j=vdims%j_start, vdims%j_end
             DO i=vdims%i_start, vdims%i_end
               v50m_h(i,j) = work_v(i,j)
             END DO ! i
           END DO ! j

         END IF ! on STASHflag qv50m_h

        IF(qv50mB_h) THEN
! Perform simple horizontal interpolation from 'C' to 'B' grid

          IF (l_vatpoles) THEN

          IF (global .AND. .NOT. qu50mB_h) THEN
            ! In the vatpoles case, v on B grid 50m is 
            !   dependent on u on B grid 50m
            ErrorStatus = -1        ! Warning
            Cmessage='v on B grid 50m not possible when '   &
               //'u on B grid 50m not active. Diagnostic aborted'

            CALL Ereport(Routinename,ErrorStatus,Cmessage)

          ELSE

            CALL vC_to_vB(work_v, rows,                                &
                          row_length,n_rows,1,offx,offy,               &
                          global_row_length,                           &
                          v50mB_h(udims%i_start,vdims%j_start)         &
                        , u50mB_h (udims%i_start,vdims%j_start)      &
                          )
          END IF
          ELSE  ! vatpoles

            CALL vC_to_vB(work_v, rows,                                &
                          row_length,n_rows,1,offx,offy,               &
                          global_row_length,                           &
                          v50mB_h(udims%i_start,vdims%j_start)         &
                          )
          END IF  ! vatpoles

         END IF ! on STASHflag qv50mB_h

      END IF ! on STASHflags qv50mB_h.OR.qv50m_h


! ----------------------------------------------------------------------
! STASH items 231,232,233,234: Test diagnostics 1-4
! ----------------------------------------------------------------------

      IF (qdia1.OR.qdia2.OR.qdia3.OR.qdia4) THEN

! DEPENDS ON: testdiag
        CALL TestDiag(                                                  &
     &  theta_field_size,v_field_size,rows,n_rows,row_length            &
     & ,ew_space,ns_space,first_lat,first_long,phi_pole,lambda_pole     &
     & ,lam                                                             &
     & ,testd_press,testd_p_levs                                        &
     & ,testd_model,testd_m_levs,forecast_hrs                           &
     & ,testdiag1,testdiag2,testdiag3,testdiag4                         &
     & ,qdia1,qdia2,qdia3,qdia4)

      END IF ! on STASHflag

! ---------------------------------------------------------------------
! STASH item 215 : theta on potential vorticity = +/-2 surface
! ---------------------------------------------------------------------

      IF (qtheta_potn_vort) THEN
        desired_potn_vort = 0.0000020    ! set pv surface to pv=2x10^-6

! Take the absolute value of pv so that interpolation is to pv = +/- 2

        DO k = 1, model_levels
          DO j =  vdims%j_start, vdims%j_end
            DO i = udims%i_start, udims%i_end
              mod_PV(i,j,k) = ABS(PV(i,j,k))
            END DO
          END DO
        END DO

! Interpolate theta onto pv=2 surface

        Call vert_interp_mdi2 (theta_at_pv, row_length, n_rows,         &
                               model_levels,                            &
                               desired_potn_vort,                       &
                               0, 0, 0, 0,                              &
                               mod_PV, interp_order_linear,             &
                               rmdi,                                    &
              theta_potn_vort(udims%i_start,vdims%j_start) )

      END IF ! on STASHflag

! ---------------------------------------------------------------------
! STASH item 216 : theta at potential vorticity points
! ---------------------------------------------------------------------

      IF (qtheta_pv_points) THEN

         DO k = 1, model_levels
            DO j = vdims%j_start, vdims%j_end
               DO i = udims%i_start, udims%i_end
                  theta_pv_points(i,j,k) = theta_at_pv(i,j,k)
               END DO
            END DO
         END DO

      END IF ! on STASHflag

! ---------------------------------------------------------------------
! STASH item 217 : potential vorticity on model levels
! ---------------------------------------------------------------------

      IF (qpv_mod_levs) THEN

         DO k = 1, model_levels
            DO j = vdims%j_start, vdims%j_end
               DO i = udims%i_start, udims%i_end
                  pv_mod_levs(i,j,k) = PV(i,j,k)
               END DO
            END DO
         END DO

      END IF ! on STASHflag

! ---------------------------------------------------------------------
! STASH item 218 : potential vorticity on model theta points and levels
! ---------------------------------------------------------------------

      IF (qpv_theta_mlev) THEN

! DEPENDS ON: calc_pv_at_theta
        Call Calc_PV_at_theta(u, v, theta, rho,                         &
     &                        r_theta_levels, r_rho_levels,             &
     &                        r_at_u, r_at_v,                           &
     &                        sec_v_latitude, tan_v_latitude,           &
     &                        sec_theta_latitude, f3_at_v,              &
     &                        delta_lambda, delta_phi,                  &
     &                        model_domain,                             &
     &                        pv_theta_mlev)

      END IF ! on STASHflag
 
      ! True rho is required for many of the rho STASH options.
      IF (qrhow .OR. qrhouw .OR. qrhovw .OR. qrhow_up .OR.              &
     &     qrhow_down .OR. qrhowc_up .OR. qrhowc_down .OR.              &
           qrho_h .OR. qtrue_density) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              true_rho(i,j,k)=rho(i,j,k)/                               &
     &                        (r_rho_levels(i,j,k)*r_rho_levels(i,j,k))
            END DO
          END DO
        END DO
      END IF

! ----------------------------------------------------------------------
! STASH items 260-266 : mass and momentum fluxes
! 260 = mass flux = rhow
! 261 = rhouw
! 262 = rhovw
! 263 = upward mass flux = rhow with w >0m/s
! 264 = downward mass flux = rhow with w <0m/s
! 265 = upward convective mass flux = rhow with w> 1m/s
! 266 = upward convective mass flux = rhow with w < -1m/s
! ----------------------------------------------------------------------
      IF (qrhow .OR. qrhouw .OR. qrhovw .OR. qrhow_up .OR.              &
     &     qrhow_down .OR. qrhowc_up .OR. qrhowc_down) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              weight = (r_theta_levels(i,j,k) - r_rho_levels(i,j,k))/   &
     &                 (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
              w_on_rho(i,j,k)= weight        * w(i,j,k-1) +             &
     &                        (1.0 - weight) * w(i,j,k)
              rhow(i,j,k) = true_rho(i,j,k)*w_on_rho(i,j,k)
              IF (qrhouw) THEN
                u_on_rho(i,j,k)=0.5*(u(i-1,j,k)+u(i,j,k))
                rhouw(i,j,k)=rhow(i,j,k)*u_on_rho(i,j,k)
              END IF
              IF (qrhovw) THEN
                v_on_rho(i,j,k)=0.5*(v(i,j-1,k)+v(i,j,k))
                rhovw(i,j,k)=rhow(i,j,k)*v_on_rho(i,j,k)
              END IF
              IF (qrhow_up) THEN
                IF (w_on_rho(i,j,k)  >   0.0) THEN
                  rhow_up(i,j,k)=rhow(i,j,k)
                ELSE
                  rhow_up(i,j,k)=0.0
                END IF
              END IF
              IF (qrhow_down) THEN
                IF (w_on_rho(i,j,k)  <   0.0) THEN
                  rhow_down(i,j,k) = rhow(i,j,k)
                ELSE
                  rhow_down(i,j,k)=0.0
                END IF
              END IF
              IF (qrhowc_up) THEN
                IF (w_on_rho(i,j,k)  >   1.0) THEN
                  rhow_convup(i,j,k)=rhow(i,j,k)
                ELSE
                  rhow_convup(i,j,k)=0.0
                END IF
              END IF
              IF (qrhowc_down) THEN
                IF (w_on_rho(i,j,k)  <   -1.0) THEN
                  rhow_convdown(i,j,k) = rhow(i,j,k)
                ELSE
                  rhow_convdown(i,j,k)=0.0
                END IF
              END IF
            END DO
          END DO
        END DO

      END IF ! on STASHflag
!-----------------------------------------------------------------------
! STASH items 101, 102, 108, 119, 127, 142, 143, 144
!-----------------------------------------------------------------------
! Height on theta-levels diagnostic

      IF ( qHtheta_ml ) THEN
         DO k=1, Htheta_m_levs
            DO j=1, rows
               DO i=1, row_length
                  kk = Htheta_model(k)
                  height_theta_ml(i,j,k)                                &
     &                        = r_theta_levels(i,j,kk) - Earth_Radius
               END DO
            END DO
         END DO
      END IF

! Height on rho-levels diagnostic

      IF ( qHrho_ml ) THEN
         DO k=1, Hrho_m_levs
            DO j=1, rows
               DO i=1, row_length
                  kk = Hrho_model(k)
                  height_rho_ml(i,j,k)                                  &
     &                 = r_rho_levels(i,j,kk) - Earth_Radius
               END DO
            END DO
         END DO
      END IF

! pressure on height-levels diagnostic

      IF ( qpress_H ) THEN

         DO k=1, p_h_levs

            desired_r = p_height(k) + Earth_Radius

            Call vert_interp_mdi (exner_rho_levels,row_length,rows,     &
     &           model_levels,desired_r,offx,offy,halo_i,halo_j,        &
     &           r_rho_levels,interp_order,rmdi,press_h(:,:,k) )

! Convert to standard pressure

            DO j=1, rows
               DO i=1, row_length
                  IF ( press_h(i,j,k)  /=  rmdi ) THEN
                     press_h(i,j,k) = p_zero *                          &
                          press_h(i,j,k)**recip_kappa
                  END IF
               END DO
            END DO

         END DO  ! end k-loop

      END IF

! potential temperature on height-levels diagnostic

      IF ( qtheta_H ) THEN
         DO k=1,theta_h_levs
            desired_r = theta_height(k) + Earth_Radius
            Call vert_interp_mdi ( theta,row_length,rows,                     &
              tdims_s%k_len,desired_r,offx,offy,halo_i,halo_j,                &
              r_theta_levels(1-halo_i,1-halo_j,tdims_s%k_start),interp_order, &
              rmdi,theta_h(:,:,k) )
         END DO
      END IF

! density (rho) on height-levels diagnostic

      IF ( qrho_H ) THEN

         DO k=1,rho_h_levs

            desired_r = rho_height(k) + Earth_Radius

            Call vert_interp_mdi (true_rho,row_length,rows,model_levels,     &
     &           desired_r,0,0,halo_i,halo_j,r_rho_levels,        &
     &           interp_order,rmdi,rho_h(:,:,k))

         END DO

      END IF

! w on height-levels diagnostic

      IF ( qw_H ) THEN
         DO k=1, w_h_levs
            desired_r = w_height(k) + Earth_Radius
            Call vert_interp_mdi (w,row_length,rows,wdims_s%k_len,      &
     &           desired_r,offx,offy,halo_i,halo_j,r_theta_levels,      &
     &           interp_order,rmdi,wcomp_h(:,:,k) )
         END DO
      END IF

! u on height-levels diagnostic

      IF ( qu_H ) THEN
         DO k=1, u_h_levs
            desired_r = u_height(k) + Earth_Radius
            Call vert_interp_mdi (u,row_length,rows,model_levels,       &
     &           desired_r,offx,offy,halo_i,halo_j,                     &
     &           r_at_u,interp_order,rmdi,ucomp_h(:,:,k))
         END DO
      END IF

! v on height-levels diagnostic

      IF ( qv_H ) THEN
         DO k=1, v_h_levs
            desired_r = v_height(k) + Earth_Radius
            Call vert_interp_mdi (v,row_length,n_rows,model_levels,     &
     &           desired_r,offx,offy,halo_i,halo_j,r_at_v,              &
     &           interp_order,rmdi,vcomp_h(:,:,k))
         END DO
      END IF

      IF (spec_w) THEN
        IF (model_domain  ==  4) THEN
          gath_proc=0
          DO klev=1,model_levels
! DEPENDS ON: calc_spectra
            CALL CALC_SPECTRA(w(wdims_s%i_start,wdims_s%j_start,klev),  &
                                                    spec_2D,            &
                   row_length+2*offx,rows+2*offy,                       &
                   global_row_length,global_rows,                       &
                   fld_type_p,halo_type_single,                         &
                   gath_proc)
            DO j = 1,rows
              DO i = 1,row_length
                spec_3D(i,j,klev)=spec_2D(i,j)
              END DO
            END DO
          END DO
        ELSE
          write(6,*)'The spectra is not set up for this domain'
        END IF
      END IF

! -----------------------------------------------------
! stash item 271  : true density on model_levels
! -----------------------------------------------------
      IF (qtrue_density) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              true_density(i,j,k)=true_rho(i,j,k)
            END DO
          END DO
        END DO
      END IF

! ---------------------------------------------------------------------
! Check error condition
      IF(ErrorStatus >  0) THEN

         CALL Ereport(RoutineName,ErrorStatus,Cmessage)
      END IF

      IF (lhook) CALL dr_hook('DYN_DIAG',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Dyn_diag
