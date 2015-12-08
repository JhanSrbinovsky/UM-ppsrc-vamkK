MODULE g_wave_5a_mod
USE dynamics_grid_mod, ONLY: l_vatpoles
USE Field_Types
USE domain_params
IMPLICIT NONE      
CONTAINS
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calls components of version 5A of gravity wave drag scheme.
!
      SUBROUTINE g_wave_5a(                                     &
        theta, u, v, row_length, rows, nrows,u_rows,v_rows,     & 
        off_x, off_y,                                           &
        global_row_length,n_proc, n_procy, proc_row_group,      &
        at_extremity,  model_domain,                            &
        levels, rho,                                            &
        delta_lambda, delta_phi, true_latitude,                 &
        sd_orog_land, orog_grad_xx_land, orog_grad_xy_land,     &
        orog_grad_yy_land,                                      &
        land_index, land_points, timestep,                      &
        kay, frc, r_u, r_v, T_inc, l_taus_scale, l_fix_gwsatn,  &
        l_gwd_40km, sat_scheme, fsat, Gsharp, fbcd,             &
        l_smooth, l_nonhydro, l_dynbeta, l_gw_heating,          &
! diagnostics
        stress_ud,     stress_ud_on,     stress_ud_p_on,        &
                                         points_stress_ud,      &
        stress_vd,     stress_vd_on,     points_stress_vd,      &
        stress_ud_satn,stress_ud_satn_on,points_stress_ud_satn, &
        stress_vd_satn,stress_vd_satn_on,points_stress_vd_satn, &
        stress_ud_wake,stress_ud_wake_on,points_stress_ud_wake, & 
        stress_vd_wake,stress_vd_wake_on,points_stress_vd_wake, &
        du_dt_satn,    du_dt_satn_on,    du_dt_satn_p_on,       & 
                                         points_du_dt_satn,     &
        dv_dt_satn,    dv_dt_satn_on,    points_dv_dt_satn,     &
        du_dt_wake,    du_dt_wake_on,    points_du_dt_wake,     &
        dv_dt_wake,    dv_dt_wake_on,    points_dv_dt_wake,     &
        u_s_d,         u_s_d_on,         points_u_s_d,          &
        v_s_d,         v_s_d_on,         points_v_s_d,          &
        nsq_s_d,       nsq_s_d_on,       points_nsq_s_d,        &
        fr_d,          fr_d_on,          points_fr_d,           &
        bld_d,         bld_d_on,         points_bld_d,          &
        bldt_d,        bldt_d_on,        points_bldt_d,         &
        num_lim_d,     num_lim_d_on,     points_num_lim_d,      &
        num_fac_d,     num_fac_d_on,     points_num_fac_d,      &
        tausx_d,       tausx_d_on,       points_tausx_d,        &
        tausy_d,       tausy_d_on,       points_tausy_d,        &
        taus_scale_d,  taus_scale_d_on,  points_taus_scale_d,   &
        orog_slope_d,  orog_slope_d_on,  points_orog_slope_d,   &
        orog_anis_d,   orog_anis_d_on,   points_orog_anis_d,    &
        orog_dir_d,    orog_dir_d_on,    points_orog_dir_d,     &
        iret)

  USE atm_fields_bounds_mod, ONLY:               &
    udims, vdims, wdims, tdims, pdims,           &
    udims_s, vdims_s, wdims_s, tdims_s, pdims_s, array_dims   

! Model level heights from centre of Earth
  USE level_heights_mod, ONLY: &
    r_theta_levels,            &  ! Radii on theta levels (m)
    r_rho_levels                  ! Radii on rho levels (m)
  USE earth_constants_mod, ONLY: earth_radius
  
  USE eg_v_at_poles_mod, ONLY: eg_v_at_poles  
  
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE Field_Types
  USE p_to_u_mod, ONLY: p_to_u
  USE p_to_v_mod, ONLY: p_to_v
  USE polar_row_mean_mod, ONLY: polar_row_mean
  USE u_to_p_mod, ONLY: u_to_p
  USE v_to_p_mod, ONLY: v_to_p
  IMPLICIT NONE      
!
! Description:
! 1) Interpolate winds to theta points
!    gather data for land points only
! 2) Call routine to set-up variables for the gravity wave and 
!    flow blocking routines
! 3) Calculate stress profiles due to different components of the
!    scheme. Calculate associated wind increments.
! 4) Interpolate acceleration to wind points and update winds
! 5) Gather diagnostics from land points and interpolate from p to
!    u,v staggering on 'c' grid
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Gravity Wave Drag
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.


!
! SUBROUTINE ARGUMENTS
!
!----------------------
! Intent(in) variables 
!----------------------   

  INTEGER, INTENT(IN)  ::        &!block for intent(in)         
    row_length,                  &! number of points per row
    rows,                        &! number of rows on theta grid
    nrows,                       &! number of rows on v grid
    u_rows,                      &! rows on u grid
    v_rows,                      &! rows on v grid
    model_domain,                &! MPP switch - global model or not.
    off_x,                       &! small x-halo
    off_y,                       &! small y-halo
    global_row_length,           &! number of points on a row
    proc_row_group,              &! Group id for processors on the same row
    n_proc,                      &! Total number of processors
    n_procy,                     &! Number of processors in latitude 
    levels,                      &! number of model levels
    land_points,                 &! number of land points
    land_index(rows*row_length)   ! index for land points

  INTEGER, INTENT(IN) ::         &!intent(in) 4a scheme only
    sat_scheme                    !Switch to determine whether to use
                                  !amplitude or stress based saturation test

  INTEGER, INTENT(IN) ::         &! intent(in) for diags
    points_stress_ud,            &
    points_stress_vd,            &
    points_stress_ud_satn,       &
    points_stress_vd_satn,       &
    points_stress_ud_wake,       &
    points_stress_vd_wake,       &
    points_du_dt_satn,           &
    points_dv_dt_satn,           &
    points_du_dt_wake,           &
    points_dv_dt_wake,           &
    points_u_s_d,                &
    points_v_s_d,                &
    points_nsq_s_d,              &
    points_fr_d,                 &
    points_bld_d,                &
    points_bldt_d,               &
    points_num_lim_d,            &
    points_num_fac_d,            &
    points_tausx_d,              &
    points_tausy_d,              &
    points_taus_scale_d,         &
    points_orog_slope_d,         &
    points_orog_anis_d,          &
    points_orog_dir_d

  REAL, INTENT(IN)  ::                            &!block for intent(in)       
    theta(row_length,rows,levels),                &!primary theta field (K)
    rho(pdims_s%i_start:pdims_s%i_end,            &!density *r*r (kg/m)
        pdims_s%j_start:pdims_s%j_end,            &!
        pdims_s%k_start:pdims_s%k_end),           &!
    delta_lambda,                                 &!spacing between pts - i dir
    delta_phi,                                    &!spacing between pts - j dir
    true_latitude(row_length, rows),              &!latitude   
    timestep,                                     &!timestep (s)
    kay,                                          &!sfc stress const [4A only]
    frc,                                          &!critical Froude number
    fsat,                                         &!crit Fr num - wave breaking
    sd_orog_land(land_points),                    &!standard deviation orog(m)
    orog_grad_xx_land(land_points),               &!(dh/dx)**2 gradient orog
    orog_grad_xy_land(land_points),               &!(dh/dx)(dh/dy) grad orog
    orog_grad_yy_land(land_points)                 !(dh/dy)**2 gradient orog


  REAL, INTENT(IN) ::    &!block for intent(in) 5a scheme only
    Gsharp,              &!function of mtn sharpness
    fbcd                  !flow blocking drag coefficient

  LOGICAL, INTENT(IN) :: &!intent(in):
    at_extremity(4)       ! indicates if this processor is at north,
                          ! south, east or west of the processor grid

  LOGICAL, INTENT(IN) :: &!intent(in) 4a scheme only
    l_taus_scale,        &!if true - sfc stress depends on the low level Fr num
    l_fix_gwsatn,        &!if true - invoke minor bug fixes in gwsatn    
    l_gwd_40km            !if true - don't apply GWD above 40km  
      
  LOGICAL, INTENT(IN) :: &!intent(in) 5a scheme only
    l_dynbeta,           &!if true:dynamically adjusting beta (group vel angle)
    l_nonhydro,          &!if true:use nonhydro scheme
    l_smooth,            &!if true:lambda_z smoothing of acc
    l_gw_heating(3)       !Calculate heating tendency
                          !l_gw_heating(1)=.true.: due to flow blocking drag
                          !l_gw_heating(2)=.true.: due to gravity wave drag

! Below are the stash flags for calculating diagnostics:
  LOGICAL, INTENT(IN) :: &!intent(in) for diags  
    stress_ud_on,        &!u stress
    stress_ud_p_on,      &!u stress p points
    stress_vd_on,        &!v stress
    stress_ud_satn_on,   &!u satn stress
    stress_vd_satn_on,   &!v satn stress
    stress_ud_wake_on,   &!u wake stress
    stress_vd_wake_on,   &!v wake stress
    du_dt_satn_on,       &!u accel (saturation)
    du_dt_satn_p_on,     &!u accel (saturation) on p points
    dv_dt_satn_on,       &!v accel (saturation)
    du_dt_wake_on,       &!u accel blocked flow
    dv_dt_wake_on,       &!v accel blocked flow
    u_s_d_on,            &!u_s_d diag switch
    v_s_d_on,            &!v_s_d diag switch
    nsq_s_d_on,          &!nsq_s_d diag switch
    fr_d_on,             &!fr_d switch
    bld_d_on,            &!bld_d switch
    bldt_d_on,           &!bldt_d switch
    num_lim_d_on,        &!num_lim_d switch (4a scheme only)
    num_fac_d_on,        &!num_fac_d switch (4a scheme only)
    tausx_d_on,          &!tausx_d switch
    tausy_d_on,          &!tausy_d switch
    taus_scale_d_on,     &!taus_scale_d switch (4a scheme only)
    orog_slope_d_on,     &!orog_slope_d switch (5a scheme only)
    orog_anis_d_on,      &!orog_anis_d switch (5a scheme only)
    orog_dir_d_on         !orog_dir_d switch (5a scheme only)

!----------------------
! Intent(inout) variables 
!----------------------
  REAL, INTENT(INOUT) ::                &!block for intent(inout)   
    u(udims_s%i_start:udims_s%i_end,    & ! primary u field (ms**-1)
      udims_s%j_start:udims_s%j_end,    &
      udims_s%k_start:udims_s%k_end),   &
    v(vdims_s%i_start:vdims_s%i_end,    & ! primary v field (ms**-1)
      vdims_s%j_start:vdims_s%j_end,    &
      vdims_s%k_start:vdims_s%k_end),   &
    T_inc(tdims%i_start:tdims%i_end,    & !Temperature increment 
          tdims%j_start:tdims%j_end,    &
                      1:tdims%k_end)   

  REAL, INTENT(INOUT) ::                &!intent(inout):
    r_u(udims_s%i_start:udims_s%i_end,  & !u wind increment diagnostic
        udims_s%j_start:udims_s%j_end,  &
        udims_s%k_start:udims_s%k_end), &
    r_v(vdims_s%i_start:vdims_s%i_end,  & !v wind increment diagnostic
        vdims_s%j_start:vdims_s%j_end,  &
        vdims_s%k_start:vdims_s%k_end)    

!----------------------
! Intent(out) variables 
!----------------------
  INTEGER, INTENT(OUT)  ::          &!block for intent(in)      
    iret                             ! return code : iret=0 normal exit

  REAL, INTENT(OUT) ::                         &!block for intent(out)
    stress_ud(row_length,u_rows,0:levels),     &!u total stress
    stress_vd(row_length,v_rows,0:levels),     &!v total stress
    stress_ud_satn(row_length,u_rows,0:levels),&!u satn stress
    stress_vd_satn(row_length,v_rows,0:levels),&!v satn stress
    stress_ud_wake(row_length,u_rows,0:levels),&!u wake stress
    stress_vd_wake(row_length,v_rows,0:levels),&!v wake stress
    du_dt_satn(row_length,u_rows,levels),      &!u accel (saturation)
    dv_dt_satn(row_length,v_rows,levels),      &!v accel (saturation)
    du_dt_wake(row_length,u_rows,levels),      &!u accel (blocked flow)
    dv_dt_wake(row_length,v_rows,levels)        !v accel (blocked flow)

  REAL              ::            &!block for 4a only (do not use INOUT)
    num_lim_d (row_length,  rows),&!% of time numerical limiter used (4a only) 
    num_fac_d (row_length,  rows),&!% redn. of flow-blocking stress (4a only) 
    taus_scale_d(row_length,rows)  !Factor surface stress scaled by (4a only) 
 
  REAL, INTENT(OUT) ::            &!block for intent(out)
    u_s_d  (row_length,  rows),   &!u_s diag at theta pts
    v_s_d  (row_length,  rows),   &!v_s diag at theta pts
    nsq_s_d(row_length,  rows),   &!n_s(not nsq) diag at theta pts
    fr_d   (row_length,  rows),   &!Fr diag at theta pts
    bld_d  (row_length,  rows),   &!blocked layer depth at theta pts
    bldt_d (row_length,  rows),   &!% of time blocked layer diagnosed
                                  !after numerical limiter invoked
    tausx_d(row_length,u_rows),   &!x-component of surface stress
    tausy_d(row_length,v_rows),   &!y-component of surface stress
    orog_slope_d(row_length,rows),&!Slope of sub-grid orography (5a only) 
    orog_anis_d(row_length,rows), &!Anisotropy of sub-grid orography (5a only) 
    orog_dir_d(row_length,rows)    !Orientation of sub-grid orography (5a only)


!--------------------------------------------------------------------
! Local variables
!--------------------------------------------------------------------

  INTEGER  ::           &! block for local variables
    k_top(land_points), &! first model level above mountain tops
    i,j,k,l              ! loop counters in routine

  INTEGER, PARAMETER :: &! parameters for MPP
    pnorth   = 1,       &! north processor address in the neighbor array
    peast    = 2,       &! east processor address in the neighbor array
    psouth   = 3,       &! south processor address in the neighbor array
    pwest    = 4,       &! west processor address in the neighbor array
    nodomain = -1        ! value in neighbor array if the domain has
                         !  no neighbor in this direction. otherwise
                         !  the value will be the tid of the neighbor

  REAL      ::                                  &!block for work arrays
    work_u(row_length,rows,levels),             &       
    work_v(row_length,rows,levels),             &      
    work_halo(wdims_s%i_start:wdims_s%i_end,    & 
              wdims_s%j_start:wdims_s%j_end,    &
              wdims_s%k_start:wdims_s%k_end),   &
    work_on_u_grid(udims%i_start:udims%i_end,   & 
                   udims%j_start:udims%j_end,   &
                   udims%k_start:udims%k_end),  &
    work_on_v_grid(vdims%i_start:vdims%i_end,   & 
                   vdims%j_start:vdims%j_end,   &
                   vdims%k_start:vdims%k_end)

  REAL ::               &!currently local variables
    slope(land_points), &!sso params - slope
    anis(land_points),  &!sso params - anisotropy
    mtdir(land_points), &!sso params - angle of major axis 
    mt_high(land_points) !sso params - height

  REAL    ::                        &!block for local variables 
    up_land(land_points,levels),    &!interpolated u on theta grid
    vp_land(land_points,levels),    &!interpolated V on theta grid
    theta_land(land_points,levels), &!land theta field on theta levels
    r_rho_levels_land               &!land field of heights of
    (land_points,levels),           &!rho levels above z=0
    r_theta_levels_land             &!land field of heights of
    (land_points,levels),           &!theta levels above z=0
    rho_land(land_points,levels),   &!density at land points
    du_dt(land_points,levels),      &!total mountain drag du/dt on land/theta
    dv_dt(land_points,levels),      &!total mountain drag dv/dt on land/theta
    dt_dt(land_points,levels),      &!temperature tendency on land/theta
    nsq(land_points,levels),        &!buoyancy frequency
    latitude_land(land_points),     &!true latitude on land points     
    banis(land_points),             &!const (function of anisotropy)
    canis(land_points),             &!const (function of anisotropy)
    ulow(land_points),              &!u averaged from 0.5h to h  
    vlow(land_points),              &!v averaged from 0.5h to h
    modu(land_points),              &!modulus of horizontal wind(sqrt(u^2+v^2))
    nlow(land_points),              &!N bulk averaged from 0.5h to h
    psilow(land_points),            &!psi averaged from 0.5h to h
    rholow(land_points),            &!rho averaged from 0.5h to h
    psi1(land_points),              &!atan(vlow/ulow)
    zb(land_points)                  !depth of flow blocking layer

  REAL, PARAMETER ::                        &
    recip_a2=1./(earth_radius*earth_radius)  !1/(radius of earth)^2

  LOGICAL   ::          &
     l_drag(land_points) !whether point has a non-zero stress or not

!--------------------------------------------------------------------
! Local diagnostic variables
!--------------------------------------------------------------------     
!
! Land points arrays below are for the total GWD stress and
! its 4 individual components. (x and y components for each)
!
  REAL  ::                                                   &
    stress_ud_land      ( points_stress_ud      , 0:levels ),&
    stress_vd_land      ( points_stress_vd      , 0:levels ),&
    stress_ud_satn_land ( points_stress_ud_satn , 0:levels ),&
    stress_vd_satn_land ( points_stress_vd_satn , 0:levels ),&
    stress_ud_wake_land ( points_stress_ud_wake , 0:levels ),&
    stress_vd_wake_land ( points_stress_vd_wake , 0:levels )
!
!  Land point arrays below are for the 4 individual components of
!  the GWD wind increment.  (x and y components for each)
!
  REAL  ::                                         &                           
    du_dt_satn_land ( points_du_dt_satn , levels ),&               
    dv_dt_satn_land ( points_dv_dt_satn , levels ),&                
    du_dt_wake_land ( points_du_dt_wake , levels ),&                
    dv_dt_wake_land ( points_dv_dt_wake , levels )

!
!  Land point arrays below are for the 9 GWD 'surface' diagnostics.
!
  REAL  ::                                 &                                  
    u_s_d_land(points_u_s_d),              &
    v_s_d_land(points_v_s_d),              &
    nsq_s_d_land(points_nsq_s_d),          &            
    fr_d_land(points_fr_d),                &                    
    bld_d_land(points_bld_d),              &                     
    bldt_d_land(points_bldt_d),            &                                  
    tausx_d_land(points_tausx_d),          &                         
    tausy_d_land(points_tausy_d)

  TYPE (array_dims)  :: gwd_udims          &  ! gwd 0 to top udims for diags
                       ,gwd_vdims          &  ! gwd 0 to top vdims for diags
                       ,gwd_udims1         &  ! gwd single level udims
                       ,gwd_vdims1            ! gwd single level vdims


  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


!------------------------------------------------------------------
! set up local dims for gwd diags as required by eg_v_at_poles
! -----------------------------------------------------------------

      gwd_udims                   = udims
      gwd_udims%k_start           = 0
      gwd_vdims                   = vdims
      gwd_vdims%k_start           = 0
      
      gwd_udims1                  = udims
      gwd_udims1%k_start          = 1
      gwd_udims1%k_end            = 1
      gwd_vdims1                  = vdims
      gwd_vdims1%k_start          = 1
      gwd_vdims1%k_end            = 1

!------------------------------------------------------------------
!l    1.1 interpolate winds to p/theta-grid
!------------------------------------------------------------------

    IF (lhook) CALL dr_hook('G_WAVE_5A',zhook_in,zhook_handle)
      CALL u_to_p(u,                                                    &
                        udims_s%i_start,udims_s%i_end,                  & 
                        udims_s%j_start,udims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        levels,                                         &
                        model_domain,at_extremity,work_u)
!
      CALL v_to_p(v,                                                    &
                        vdims_s%i_start,vdims_s%i_end,                  & 
                        vdims_s%j_start,vdims_s%j_end,                  &
                        pdims%i_start,pdims%i_end,                      & 
                        pdims%j_start,pdims%j_end,                      &
                        levels,                                         &
                        model_domain,at_extremity,work_v)


IF (.NOT. l_vatpoles) THEN
! set polar values of u and v to zero.
      IF (model_domain  ==  mt_global) THEN
        IF (at_extremity(psouth)) THEN
          DO k = 1, levels
            DO i = 1, row_length
               work_v(i,pdims%j_start,k) = 0.0
               work_u(i,pdims%j_start,k) = 0.0
            END DO
          END DO
        END IF
        IF (at_extremity(pnorth)) THEN
          DO k = 1, levels
            DO i = 1, row_length
               work_v(i,pdims%j_end,k) = 0.0
               work_u(i,pdims%j_end,k) = 0.0
            END DO
          END DO
        END IF
      END IF  ! model_domain  ==  mt_global
ELSE
   ! No special conditions for vatpoles grid as p points not on the pole
END IF ! vatpoles

!------------------------------------------------------------------
!l    1.2  gather winds at land points
!------------------------------------------------------------------

    DO k=1,levels
      DO l=1,land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          up_land(l,k) =work_u(i,j,k)
          vp_land(l,k) =work_v(i,j,k)
      END DO
    END DO

!------------------------------------------------------------------
!l    1.3  gather theta, rho and heights at land points
!------------------------------------------------------------------

    IF (land_points  >   0) THEN
      DO k=1,levels
        DO l=1,land_points
          j = (land_index(l)-1)/row_length + 1
          i = land_index(l) - (j-1)*row_length
          r_rho_levels_land(l,k)  = r_rho_levels(i,j,k) -               &
                                    r_theta_levels(i,j,0) 
          r_theta_levels_land(l,k)= r_theta_levels(i,j,k) -             &
                                    r_theta_levels(i,j,0)
          rho_land(l,k)           = rho(i,j,k)*recip_a2
          theta_land(l,k)         = theta(i,j,k)
          latitude_land(l)       = true_latitude(i,j)  
        END DO
      END DO

!DEPENDS ON: gw_setup
      CALL GW_SETUP(levels,land_points,up_land,vp_land,rho_land,theta_land, &
                    nsq,ulow,vlow,rholow,nlow,psilow,psi1,modu,             &
                    r_rho_levels_land,r_theta_levels_land,                  &
                    sd_orog_land,orog_grad_xx_land,orog_grad_xy_land,       &
                    orog_grad_yy_land,mt_high,slope,anis,banis,canis,mtdir, &
                    k_top,l_drag,                                           &
!diagnostics
                    u_s_d_land,u_s_d_on,points_u_s_d,                       &
                    v_s_d_land,v_s_d_on,points_v_s_d,                       &
                    nsq_s_d_land, nsq_s_d_on, points_nsq_s_d)

!DEPENDS ON: gw_block
      CALL GW_BLOCK(levels,land_points,timestep,up_land,vp_land,rho_land,   &
                   nsq,ulow,vlow,rholow,psilow,modu,r_rho_levels_land,      &
                   r_theta_levels_land,mt_high,sd_orog_land,                &
                   slope,anis,mtdir,zb,banis,                               &
                   canis,du_dt,dv_dt,dt_dt,fbcd,frc,l_drag,l_gw_heating(1), &
!diagnostics
                   du_dt_wake_land,points_du_dt_wake,du_dt_wake_on,         &
                   dv_dt_wake_land,points_dv_dt_wake,dv_dt_wake_on,         &
                   stress_ud_land,stress_ud_on,points_stress_ud,            &
                   stress_ud_p_on,                                          &
                   stress_vd_land,stress_vd_on,points_stress_vd,            &
                   stress_ud_wake_land,points_stress_ud_wake,               &
                   stress_ud_wake_on,                                       &
                   stress_vd_wake_land,points_stress_vd_wake,               &
                   stress_vd_wake_on,                                       &
                   fr_d_land,fr_d_on,points_fr_d,                           &
                   bld_d_land,bld_d_on,points_bld_d,                        & 
                   bldt_d_land,bldt_d_on, points_bldt_d,                    &
                   tausx_d_land,tausx_d_on,points_tausx_d,                  &
                   tausy_d_land,tausy_d_on,points_tausy_d)

!------------------------------------------------------------------
!l    3. calculate stress profile and accelerations,
!l       CALL gw_vert
!------------------------------------------------------------------
!DEPENDS ON: gw_wave
      CALL GW_WAVE(levels,land_points,up_land,vp_land,rho_land,nsq,ulow,vlow, &
                   rholow,psi1,psilow,nlow,modu,k_top,r_rho_levels_land,      &
                   r_theta_levels_land,delta_lambda,delta_phi,latitude_land,  &
                   mt_high,sd_orog_land,slope,zb,banis,canis,du_dt,dv_dt,     &
                   dt_dt,timestep,l_dynbeta,l_nonhydro,l_smooth,fsat,Gsharp,  &
                   l_drag,l_gw_heating(2),                                    &
!diagnostics
                   du_dt_satn_land,points_du_dt_satn,du_dt_satn_on,           &
                   du_dt_satn_p_on,                                           &
                   dv_dt_satn_land,points_dv_dt_satn,dv_dt_satn_on,           &
                   stress_ud_land ,points_stress_ud ,stress_ud_on,            &
                   stress_ud_p_on,                                            &
                   stress_vd_land ,points_stress_vd ,stress_vd_on,            &
                   stress_ud_satn_land,points_stress_ud_satn,                 &
                   stress_ud_satn_on,                                         &
                   stress_vd_satn_land,points_stress_vd_satn,                 &
                   stress_vd_satn_on,                                         &
                   tausx_d_land   , tausx_d_on   , points_tausx_d  ,          &
                   tausy_d_land   , tausy_d_on   , points_tausy_d)
    END IF ! on land_points > 0

!------------------------------------------------------------------
!l    4. scatter accelerations to full area, interpolate to uv-grid
!l       and update winds
!------------------------------------------------------------------

! initialise work array: required to ensure zero values generated by
!   non-land points. [note that it is not necessary to re-initialise
!   the work array for each subsequent diagnostic since only land mask
!   points will differ, and these will be overwritten.]

 
    DO k = 0, levels
      DO j = wdims_s%j_start,wdims_s%j_end
        DO i = wdims_s%i_start,wdims_s%i_end
          work_halo(i,j,k)=0.
        END DO !i
      END DO  !j
    END DO   !k


! expand from land points and interpolate to 'c' u,v grid
    DO k=1,levels
      DO l=1,land_points
         j=(land_index(l)-1)/row_length + 1
         i=land_index(l) - (j-1)*row_length
         work_halo(i,j,k) = du_dt(l,k)
      END DO ! l
    END DO ! k

!DEPENDS ON: swap_bounds
     CALL swap_bounds(                                                 &
          work_halo(pdims_s%i_start,pdims_s%j_start,1),                &
          row_length, rows,                                            &
          levels, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_u(work_halo(pdims_s%i_start,pdims_s%j_start,1),   &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   1,levels,work_on_u_grid)

    DO k=1,levels
      DO l=1,land_points
         j=(land_index(l)-1)/row_length + 1
         i=land_index(l) - (j-1)*row_length
         work_halo(i,j,k) = dv_dt(l,k)
      END DO ! l
    END DO ! k

!DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                &
          work_halo(pdims_s%i_start,pdims_s%j_start,1),                &
          row_length, rows,                                            &
          levels, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_v(work_halo(pdims_s%i_start,pdims_s%j_start,1),   &   
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   vdims%i_start,vdims%i_end,                     &
                   vdims%j_start,vdims%j_end,                     &
                   1,levels,work_on_v_grid)

IF (l_vatpoles) THEN
! correct dv/dt at poles
      IF( model_domain == mt_global ) THEN
        IF( at_extremity(psouth) ) THEN


          CALL eg_v_at_poles(work_on_u_grid, work_on_v_grid, 1.0, &
                             udims%j_start, vdims%j_start,        &
                             udims,vdims)
        
        END IF
        
        IF( at_extremity(pnorth) ) THEN


          CALL eg_v_at_poles(work_on_u_grid, work_on_v_grid, -1.0,&
                             udims%j_end, vdims%j_end,            &
                             udims,vdims)
        END IF
      END IF                       
END IF ! vatpoles            
    DO k=1,levels
      DO j=udims%j_start,udims%j_end 
        DO i=udims%i_start,udims%i_end 
          r_u(i,j,k)=r_u(i,j,k)+ timestep*work_on_u_grid(i,j,k)
        END DO
      END DO
    END DO

    DO k=1,levels
      DO j=vdims%j_start,vdims%j_end 
        DO i=vdims%i_start,vdims%i_end 
          r_v(i,j,k)=r_v(i,j,k)+ timestep*work_on_v_grid(i,j,k)
        END DO
      END DO
    END DO

!------------------------------------------------------------------
!l    5. scatter temperature increments to full area,
!l       and update 
!------------------------------------------------------------------

    IF ( l_gw_heating(1) .OR. l_gw_heating(2)) THEN !Apply heating tendency

     DO k=1,levels
       DO l=1,land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          T_inc(i,j,k)=T_inc(i,j,k)+ timestep*dt_dt(l,k)
       END DO ! l
     END DO ! k

    END IF ! l_gw_heating

!------------------------------------------------------------------
!l    6. gather diagnostics from land points and interpolate from p to
!l       u,v staggering on 'c' grid.
!l       note that stress_ud_on,stress_vd_on are on theta levels
!l       whereas remaining diagnostics are on rho levels.
!------------------------------------------------------------------

    IF (stress_ud_on .OR. stress_ud_p_on) THEN  ! stashflag for this diagnostic
! expand from land points and interpolate to 'c' u,v grid
      DO k=0,levels
        DO l=1,land_points
           j=(land_index(l)-1)/row_length + 1
           i=land_index(l) - (j-1)*row_length
           work_halo(i,j,k) = stress_ud_land(l,k)
        END DO ! l
      END DO ! k

!DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                &
          work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
          row_length, rows,                                            &
          levels+1, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_u(work_halo,                                      &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   0,levels,stress_ud)


    END IF ! (stress_ud_on .OR. stress_ud_p_on) on stashflag

    IF (stress_vd_on) THEN  ! stashflag for this diagnostic
! expand from land points and interpolate to 'c' u,v grid
      DO k=0,levels
        DO l=1,land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          work_halo(i,j,k) = stress_vd_land(l,k)
        END DO ! l
      END DO ! k

!DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                &
          work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
          row_length, rows,                                            &
          levels+1, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_v(work_halo,                                      &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   vdims%i_start,vdims%i_end,                     &
                   vdims%j_start,vdims%j_end,                     &
                   0,levels,stress_vd)

IF (l_vatpoles) THEN
! correct stress at poles
      IF( model_domain == mt_global ) THEN
        IF( at_extremity(psouth) ) THEN


          CALL eg_v_at_poles(stress_ud, stress_vd, 1.0,           &
                             udims%j_start, vdims%j_start,        &
                             gwd_udims,gwd_vdims)
        
        END IF
        
        IF( at_extremity(pnorth) ) THEN


          CALL eg_v_at_poles(stress_ud, stress_vd, -1.0,          &
                             udims%j_end, vdims%j_end,            &
                             gwd_udims,gwd_vdims)
        
        END IF           
      END IF
END IF ! vatpoles

    END IF ! (stress_vd_on) on stashflag

    IF (stress_ud_satn_on) THEN  ! stashflag for this diagnostic
! expand from land points and interpolate to 'c' u,v grid
      DO k=0,levels
        DO l=1,land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          work_halo(i,j,k) = stress_ud_satn_land(l,k)
        END DO ! l
      END DO ! k

!DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                &
          work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
          row_length, rows,                                            &
          levels+1, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_u(work_halo,                                      &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   0,levels,stress_ud_satn)


    END IF ! (stress_ud_satn_on) on stashflag

    IF (stress_vd_satn_on) THEN  ! stashflag for this diagnostic
! expand from land points and interpolate to 'c' u,v grid
      DO k=0,levels
        DO l=1,land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          work_halo(i,j,k) = stress_vd_satn_land(l,k)
        END DO ! l
      END DO ! k

!DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                &
          work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
          row_length, rows,                                            &
          levels+1, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_v(work_halo,                                      &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   vdims%i_start,vdims%i_end,                     &
                   vdims%j_start,vdims%j_end,                     &
                   0,levels,stress_vd_satn)

IF (l_vatpoles) THEN
! correct stress at poles
      IF( model_domain == mt_global ) THEN
        IF( at_extremity(psouth) ) THEN


          CALL eg_v_at_poles(stress_ud_satn, stress_vd_satn, 1.0, &
                             udims%j_start, vdims%j_start,        &
                             gwd_udims,gwd_vdims)
        
        END IF
        
        IF( at_extremity(pnorth) ) THEN


          CALL eg_v_at_poles(stress_ud_satn, stress_vd_satn, -1.0,&
                             udims%j_end, vdims%j_end,            &
                             gwd_udims,gwd_vdims)
        
        END IF           
      END IF
END IF ! vatpoles



   END IF ! (stress_vd_satn_on) on stashflag

   IF (stress_ud_wake_on) THEN  ! stashflag for this diagnostic
! expand from land points and interpolate to 'c' u,v grid
     DO k=0,levels
       DO l=1,land_points
         j=(land_index(l)-1)/row_length + 1
         i=land_index(l) - (j-1)*row_length
         work_halo(i,j,k) = stress_ud_wake_land(l,k)
       END DO ! l
     END DO ! k

!DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                &
          work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
          row_length, rows,                                            &
          levels+1, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_u(work_halo,                                      &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   0,levels,stress_ud_wake)

   END IF ! (stress_ud_wake_on) on stashflag

   IF (stress_vd_wake_on) THEN  ! stashflag for this diagnostic
! expand from land points and interpolate to 'c' u,v grid
     DO k=0,levels
       DO l=1,land_points
         j=(land_index(l)-1)/row_length + 1
         i=land_index(l) - (j-1)*row_length
         work_halo(i,j,k) = stress_vd_wake_land(l,k)
       END DO ! l
     END DO ! k

!DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                &
          work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
          row_length, rows,                                            &
          levels+1, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_v(work_halo,                                      &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   vdims%i_start,vdims%i_end,                     &
                   vdims%j_start,vdims%j_end,                     &
                   0,levels,stress_vd_wake)

IF (l_vatpoles) THEN
! correct stress at poles
      IF( model_domain == mt_global ) THEN
        IF( at_extremity(psouth) ) THEN


          CALL eg_v_at_poles(stress_ud_wake, stress_vd_wake, 1.0, &
                             udims%j_start, vdims%j_start,        &
                             gwd_udims,gwd_vdims)
        
        END IF
        
        IF( at_extremity(pnorth) ) THEN


          CALL eg_v_at_poles(stress_ud_wake, stress_vd_wake, -1.0,&
                             udims%j_end, vdims%j_end,            &
                             gwd_udims,gwd_vdims)
        
        END IF           
      END IF
END IF ! vatpoles

   END IF ! (stress_vd_wake_on) on stashflag

   IF (du_dt_satn_on .OR. du_dt_satn_p_on) THEN ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
      DO k=1,levels
        DO l=1,land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          work_halo(i,j,k) = du_dt_satn_land(l,k)
        END DO ! l
      END DO ! k

!DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                &
          work_halo(pdims_s%i_start,pdims_s%j_start,1),                &
          row_length, rows,                                            &
          levels, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_u(work_halo(pdims_s%i_start,pdims_s%j_start,1),   & 
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   1,levels,du_dt_satn)

    END IF ! (du_dt_satn_on .OR. du_dt_satn_p_on) on stashflag

    IF (dv_dt_satn_on) THEN ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
     DO k=1,levels
       DO l=1,land_points
         j=(land_index(l)-1)/row_length + 1
         i=land_index(l) - (j-1)*row_length
         work_halo(i,j,k) = dv_dt_satn_land(l,k)
       END DO ! l
     END DO ! k

!DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                &
          work_halo(pdims_s%i_start,pdims_s%j_start,1),                &
          row_length, rows,                                            &
          levels, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_v(work_halo(pdims_s%i_start,pdims_s%j_start,1),   &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   vdims%i_start,vdims%i_end,                     &
                   vdims%j_start,vdims%j_end,                     &
                   1,levels,dv_dt_satn)
IF (l_vatpoles) THEN
! correct stress at poles
      IF( model_domain == mt_global ) THEN
        IF( at_extremity(psouth) ) THEN


          CALL eg_v_at_poles(du_dt_satn, dv_dt_satn, 1.0,         &
                             udims%j_start, vdims%j_start,        &
                             udims,vdims)
        
        END IF
        
        IF( at_extremity(pnorth) ) THEN


          CALL eg_v_at_poles(du_dt_satn, dv_dt_satn, -1.0,        &
                             udims%j_end, vdims%j_end,            &
                             udims,vdims)
        
        END IF           
      END IF
END IF ! vatpoles        


   END IF ! (dv_dt_satn_on) on stashflag

   IF (du_dt_wake_on) THEN  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
     DO k=1,levels
       DO l=1,land_points
         j=(land_index(l)-1)/row_length + 1
         i=land_index(l) - (j-1)*row_length
         work_halo(i,j,k) = du_dt_wake_land(l,k)
       END DO ! l
     END DO ! k

!DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                &
          work_halo(pdims_s%i_start,pdims_s%j_start,1),                &
          row_length, rows,                                            &
          levels, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_u(work_halo(pdims_s%i_start,pdims_s%j_start,1),   &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   1,levels,du_dt_wake)

   END IF ! (du_dt_wake_on) on stashflag

   IF (dv_dt_wake_on) THEN  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
     DO k=1,levels
       DO l=1,land_points
         j=(land_index(l)-1)/row_length + 1
         i=land_index(l) - (j-1)*row_length
         work_halo(i,j,k) = dv_dt_wake_land(l,k)
       END DO ! l
     END DO ! k

!DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                &
          work_halo(pdims_s%i_start,pdims_s%j_start,1),                &
          row_length, rows,                                            &
          levels, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_v(work_halo(pdims_s%i_start,pdims_s%j_start,1),   &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   vdims%i_start,vdims%i_end,                     &
                   vdims%j_start,vdims%j_end,                     &
                   1,levels,dv_dt_wake)
IF (l_vatpoles) THEN
! correct stress at poles
      IF( model_domain == mt_global ) THEN
        IF( at_extremity(psouth) ) THEN


          CALL eg_v_at_poles(du_dt_wake, dv_dt_wake, 1.0,         &
                             udims%j_start, vdims%j_start,        &
                             udims,vdims)
        
        END IF
        
        IF( at_extremity(pnorth) ) THEN


          CALL eg_v_at_poles(du_dt_wake, dv_dt_wake, -1.0,        &
                             udims%j_end, vdims%j_end,            &
                             udims,vdims)
        
        END IF           
      END IF
END IF ! vatpoles
   END IF ! (dv_dt_wake_on) on stashflag

   IF (u_s_d_on) THEN  ! stashflag for this diagnostic
! expand from land points
     DO j=1,rows
       DO i=1,row_length
         u_s_d(i,j) = 0.0
       END DO
     END DO
     DO l=1,land_points
       j=(land_index(l)-1)/row_length + 1
       i=land_index(l) - (j-1)*row_length
       u_s_d(i,j) = u_s_d_land(l)
     END DO ! l
   END IF ! (u_s_d_on) on stashflag

   IF (v_s_d_on) THEN  ! stashflag for this diagnostic
! expand from land points
     DO j=1,rows
       DO i=1,row_length
         v_s_d(i,j) = 0.0
       END DO
     END DO
     DO l=1,land_points
       j=(land_index(l)-1)/row_length + 1
       i=land_index(l) - (j-1)*row_length
       v_s_d(i,j) = v_s_d_land(l)
     END DO ! l
   END IF ! (v_s_d_on) on stashflag

   IF (nsq_s_d_on) THEN  ! stashflag for this diagnostic
! expand from land points
     DO j=1,rows
       DO i=1,row_length
         nsq_s_d(i,j) = 0.0
       END DO
     END DO
     DO l=1,land_points
       j=(land_index(l)-1)/row_length + 1
       i=land_index(l) - (j-1)*row_length
       nsq_s_d(i,j) = nsq_s_d_land(l)
     END DO ! l
   END IF ! (nsq_s_d_on) on stashflag

   IF (fr_d_on) THEN  ! stashflag for this diagnostic
! expand from land points
     DO j=1,rows
       DO i=1,row_length
         fr_d(i,j) = 0.0
       END DO
     END DO
     DO l=1,land_points
       j=(land_index(l)-1)/row_length + 1
       i=land_index(l) - (j-1)*row_length
       fr_d(i,j) = fr_d_land(l)
     END DO ! l
   END IF ! (fr_d_on) on stashflag

   IF (bld_d_on) THEN  ! stashflag for this diagnostic
! expand from land points
     DO j=1,rows
       DO i=1,row_length
         bld_d(i,j) = 0.0
       END DO
     END DO
     DO l=1,land_points
       j=(land_index(l)-1)/row_length + 1
       i=land_index(l) - (j-1)*row_length
       bld_d(i,j) = bld_d_land(l)
     END DO ! l
   END IF ! (bld_d_on) on stashflag

   IF (bldt_d_on) THEN  ! stashflag for this diagnostic
! expand from land points
     DO j=1,rows
       DO i=1,row_length
         bldt_d(i,j) = 0.0
       END DO
     END DO
     DO l=1,land_points
       j=(land_index(l)-1)/row_length + 1
       i=land_index(l) - (j-1)*row_length
       bldt_d(i,j) = bldt_d_land(l)
     END DO ! l
   END IF ! (bldt_d_on) on stashflag
 
   IF (tausx_d_on) THEN  ! stashflag for this diagnostic
! expand from land points and interpolate to 'c' u,v grid
     DO l=1,land_points
       j=(land_index(l)-1)/row_length + 1
       i=land_index(l) - (j-1)*row_length
       work_halo(i,j,0) = tausx_d_land(l)
     END DO ! l
     
!DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                &
          work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
          row_length, rows,                                            &
          1, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_u(work_halo(pdims_s%i_start,pdims_s%j_start,0),   &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   1,1,tausx_d(udims%i_start,udims%j_start))

   END IF ! (tausx_d_on) on stashflag

   IF (tausy_d_on) THEN  ! stashflag for this diagnostic
! expand from land points and interpolate to 'c' u,v grid
     DO l=1,land_points
       j=(land_index(l)-1)/row_length + 1
       i=land_index(l) - (j-1)*row_length
       work_halo(i,j,0) = tausy_d_land(l)
     END DO ! l

!DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                &
          work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
          row_length, rows,                                            &
          1, off_x, off_y, fld_type_p, .false.)
          
      CALL p_to_v(work_halo(pdims_s%i_start,pdims_s%j_start,0),   &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   vdims%i_start,vdims%i_end,                     &
                   vdims%j_start,vdims%j_end,                     &
                   1,1,tausy_d(vdims%i_start,vdims%j_start))
IF (l_vatpoles) THEN
! correct stress at poles
      IF( model_domain == mt_global ) THEN
        IF( at_extremity(psouth) ) THEN


          CALL eg_v_at_poles(tausx_d,tausy_d, 1.0,                &
                             udims%j_start, vdims%j_start,        &
                             gwd_udims1,gwd_vdims1)
        
        END IF
        
        IF( at_extremity(pnorth) ) THEN


          CALL eg_v_at_poles(tausx_d,tausy_d, -1.0,               &
                             udims%j_end, vdims%j_end,            &
                             gwd_udims1,gwd_vdims1)
        
        END IF           
      END IF
END IF ! vatpoles

   END IF ! (tausy_d_on) on stashflag

   IF (orog_slope_d_on) THEN  ! stashflag for this diagnostic
! expand from land points
     DO j=1,rows
       DO i=1,row_length
         orog_slope_d(i,j) = 0.0
       END DO
     END DO
     DO l=1,land_points
       j=(land_index(l)-1)/row_length + 1
       i=land_index(l) - (j-1)*row_length
       orog_slope_d(i,j) = slope(l)
     END DO ! l
   END IF ! (orog_slope_d_on) on stashflag

   IF (orog_anis_d_on) THEN  ! stashflag for this diagnostic
! expand from land points
     DO j=1,rows
       DO i=1,row_length
         orog_anis_d(i,j) = 0.0
       END DO
     END DO
     DO l=1,land_points
       j=(land_index(l)-1)/row_length + 1
       i=land_index(l) - (j-1)*row_length
       orog_anis_d(i,j) = anis(l)
     END DO ! l
   END IF ! (orog_anis_d_on) on stashflag

   IF (orog_dir_d_on) THEN  ! stashflag for this diagnostic
! expand from land points
     DO j=1,rows
       DO i=1,row_length
         orog_dir_d(i,j) = 0.0
       END DO
     END DO
     DO l=1,land_points
       j=(land_index(l)-1)/row_length + 1
       i=land_index(l) - (j-1)*row_length
       orog_dir_d(i,j) = mtdir(l)
     END DO ! l
   END IF ! (orog_dir_d_on) on stashflag

   iret=0
   IF (lhook) CALL dr_hook('G_WAVE_5A',zhook_out,zhook_handle)
   RETURN

END SUBROUTINE g_wave_5a
END MODULE
