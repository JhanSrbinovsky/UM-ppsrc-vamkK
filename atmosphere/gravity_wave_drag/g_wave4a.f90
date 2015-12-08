MODULE g_wave_4a_mod
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
! Calls components of version 4A of gravity wave drag scheme.
!
      SUBROUTINE g_wave_4a(                                             &
     &  theta,u,v,row_length,rows,nrows,u_rows,v_rows,off_x,off_y,      &
     &  global_row_length,n_proc, n_procy, proc_row_group,at_extremity, &
     &  model_domain,                                                   &
     &  levels,rho,                                                     &
     &  delta_lambda,delta_phi,true_latitude,                           &
     &  sd_orog_land,                                                   &
     &  orog_grad_xx_land,orog_grad_xy_land,orog_grad_yy_land,          &
     &  land_index,land_points,timestep,kay,frc,r_u,r_v,T_inc,          &
     &  l_taus_scale,l_fix_gwsatn,l_gwd_40km,                           &
     &  sat_scheme,fsat,                                                &
     &  Gsharp,fbcd,l_smooth,l_nonhydro,l_dynbeta,l_gw_heating,         &
! diagnostics
     &  stress_ud     ,stress_ud_on     , stress_ud_p_on,               &
     &                                    points_stress_ud,             &
     &  stress_vd     ,stress_vd_on     , points_stress_vd    ,         &
     &  stress_ud_satn,stress_ud_satn_on,points_stress_ud_satn,         &
     &  stress_vd_satn,stress_vd_satn_on,points_stress_vd_satn,         &
     &  stress_ud_wake,stress_ud_wake_on,points_stress_ud_wake,         &
     &  stress_vd_wake,stress_vd_wake_on,points_stress_vd_wake,         &
     &  du_dt_satn    ,du_dt_satn_on    ,du_dt_satn_p_on,               &
     &                                   points_du_dt_satn,             &
     &  dv_dt_satn    ,dv_dt_satn_on    ,points_dv_dt_satn    ,         &
     &  du_dt_wake    ,du_dt_wake_on    ,points_du_dt_wake    ,         &
     &  dv_dt_wake    ,dv_dt_wake_on    ,points_dv_dt_wake    ,         &
     &  u_s_d         ,u_s_d_on         ,points_u_s_d         ,         &
     &  v_s_d         ,v_s_d_on         ,points_v_s_d         ,         &
     &  nsq_s_d       ,nsq_s_d_on       ,points_nsq_s_d       ,         &
     &  fr_d          ,fr_d_on          ,points_fr_d          ,         &
     &  bld_d         ,bld_d_on         ,points_bld_d         ,         &
     &  bldt_d        ,bldt_d_on        ,points_bldt_d        ,         &
     &  num_lim_d     ,num_lim_d_on     ,points_num_lim_d     ,         &
     &  num_fac_d     ,num_fac_d_on     ,points_num_fac_d     ,         &
     &  tausx_d       ,tausx_d_on       ,points_tausx_d       ,         &
     &  tausy_d       ,tausy_d_on       ,points_tausy_d       ,         &
     &  taus_scale_d  ,taus_scale_d_on  ,points_taus_scale_d  ,         &
     &  orog_slope_d  ,orog_slope_d_on  ,points_orog_slope_d  ,         &
     &  orog_anis_d   ,orog_anis_d_on   ,points_orog_anis_d   ,         &
     &  orog_dir_d    ,orog_dir_d_on    ,points_orog_dir_d    ,         &
     &  iret)
      USE atm_fields_bounds_mod, ONLY:              &
        udims, vdims, wdims, tdims, pdims,          &
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
! 2) Call surface stress routine
! 3) Calculate stress profiles due to different components of the
!    scheme. Calculate associated wind increments.
! 4) Interpolate acceleration to wind points and update winds
! 5) Gather diagnostics from land points and interpolate from p to
!    u,v staggering on 'c' grid
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Gravity Wave Drag
!
!   language: fortran 77 + common extensions.
!   this code is written to umdp3 v6 programming standards.
!
! suitable for single column use, with calls to: uv_to_p removed
!                                                p_to_uv removed
! suitable for rotated grids
!
!
! SUBROUTINE ARGUMENTS
!
      INTEGER                                                           &
                           !,intent(in):
     & row_length          ! number of points per row

      INTEGER                                                           &
                           !,intent(in):
     & rows                ! number of rows on theta grid

      INTEGER                                                           &
                           !,intent(in):
     & nrows               ! number of rows on v grid

      INTEGER                                                           &
                           !,intent(in):
     & u_rows                                                           &
                           ! rows on u grid
     &,v_rows              ! rows on v grid

      INTEGER                                                           &
                           !,intent(in):
     & model_domain        ! MPP switch - global model or not.


      INTEGER             &!,intent(in):
     & off_x              &! small x-halo
     &,off_y              &! small y-halo
     &,global_row_length  &! number of points on a row
     &,proc_row_group     &! Group id for processors on the same row
     &,n_proc             &! Total number of processors
     &,n_procy             ! Number of processors in latitude 

      INTEGER                                                           &
                           !,intent(in):
     & levels              ! number of model levels

      INTEGER                                                           &
                           !,intent(in):
     & land_points         ! number of land points

      INTEGER                                                           &
                           !,intent(in):
     & land_index(rows*row_length)
!                          ! index for land points

      INTEGER                                                           &
                           !,intent(in):
     &  iret               ! return code : iret=0 normal exit

!
! The integers below are set to size land_points if the corresponding
! diagnostic is called or to 1 if it is not. These are set in GWD_CTL2
!
      INTEGER                                                           &
                           !,intent(in):
     &  points_stress_ud                                                &
     &, points_stress_vd                                                &
     &, points_stress_ud_satn                                           &
     &, points_stress_vd_satn                                           &
     &, points_stress_ud_wake                                           &
     &, points_stress_vd_wake                                           &
     &, points_du_dt_satn                                               &
     &, points_dv_dt_satn                                               &
     &, points_du_dt_wake                                               &
     &, points_dv_dt_wake                                               &
     &, points_u_s_d                                                    &
     &, points_v_s_d                                                    &
     &, points_nsq_s_d                                                  &
     &, points_fr_d                                                     &
     &, points_bld_d                                                    &
     &, points_bldt_d                                                   &
     &, points_num_lim_d                                                &
     &, points_num_fac_d                                                &
     &, points_tausx_d                                                  &
     &, points_tausy_d                                                  &
     &, points_taus_scale_d                                             &
     &, points_orog_slope_d                                             &
     &, points_orog_anis_d                                              &
     &, points_orog_dir_d

      REAL                                                              &
                           !,intent(in):
     & theta(row_length,rows,levels)
                           ! primary theta field (K)

      REAL                                                              &   
                           !,intent(in):
        rho(pdims_s%i_start:pdims_s%i_end,  & !  density *r*r (kg/m)
            pdims_s%j_start:pdims_s%j_end,  & !
            pdims_s%k_start:pdims_s%k_end)    

      REAL                                                              &
                           !block for intent(inout)   
        u(udims_s%i_start:udims_s%i_end,    & ! primary u field (ms**-1)
     &    udims_s%j_start:udims_s%j_end,    &
     &    udims_s%k_start:udims_s%k_end)    &
     &, v(vdims_s%i_start:vdims_s%i_end,    & ! primary v field (ms**-1)
     &    vdims_s%j_start:vdims_s%j_end,    &
     &    vdims_s%k_start:vdims_s%k_end)    

      REAL  ::                         &!block for intent(in) 
       delta_lambda,                   &!spacing between pts - i dir
       delta_phi,                      &!spacing between pts - j dir
       true_latitude(row_length, rows) !latitude   
      REAL                                                              &
                           !,intent(in):
     & sd_orog_land(land_points)
!                          ! standard deviation of orography (m)

      REAL                                                              &
                           !,intent(in):
     & orog_grad_xx_land(land_points)
!                          ! dh/dx squared gradient orography

      REAL                                                              &
                           !,intent(in):
     & orog_grad_xy_land(land_points)
!                          ! (dh/dx)(dh/dy) gradient orography

      REAL                                                              &
                           !,intent(in):
     & orog_grad_yy_land(land_points)
                           ! dh/dy squared gradient orography

      REAL                                                              &
                           !,intent(in):
     & timestep            ! timestep (s)

      REAL                                                              &
                           !,intent(in):
     & kay                 ! surface stress constant ( m**-1)

      REAL                                                              &
                           !,intent(in):
     & frc                 ! critical Froude number

      REAL                                                              &
                            !,intent(in):
     & fsat                 ! Froude number used to scale critical
!                           ! wave amplitude for breaking

!
!  Start of full field diagnostic arrays. This space is allocated
! in GWD_CTL2 only if a diagnostic is called.
!
      REAL                                    &!intent(inout):
     &  r_u(udims_s%i_start:udims_s%i_end,    & !u wind increment diagnostic
     &      udims_s%j_start:udims_s%j_end,    &
     &      udims_s%k_start:udims_s%k_end)    &
     &, r_v(vdims_s%i_start:vdims_s%i_end,    & !v wind increment diagnostic
     &      vdims_s%j_start:vdims_s%j_end,    &
     &      vdims_s%k_start:vdims_s%k_end)    &
     &, T_inc(tdims%i_start:tdims%i_end,      & !Temperature increment 
              tdims%j_start:tdims%j_end,      &
                            1:tdims%k_end)


      REAL                                                              &
                           !,intent(out):
     & stress_ud(row_length, u_rows,0:levels)                            &
                                                  !u   total stress
     &,stress_vd(row_length, v_rows,0:levels)      !v   total stress

      REAL                                                              &
                           !,intent(out):
     & stress_ud_satn(row_length, u_rows,0:levels)                       &
                                                  !u   satn  stress
     &,stress_vd_satn(row_length, v_rows,0:levels) !v   satn  stress

      REAL                                                              &
                           !,intent(out):
     & stress_ud_wake(row_length, u_rows,0:levels)                       &
                                                  !u   wake  stress
     &,stress_vd_wake(row_length, v_rows,0:levels) !v   wake  stress

      REAL                                                              &
                           !,intent(out):
     & du_dt_satn(row_length, u_rows,levels)                             &
                                            !u acceln (saturation)
     &,dv_dt_satn(row_length, v_rows,levels) !v acceln (saturation)

      REAL                                                              &
                           !,intent(out):
     & du_dt_wake(row_length, u_rows,levels)                             &
                                            !u acceln (blocked flow)
     &,dv_dt_wake(row_length, v_rows,levels) !v acceln (blocked flow)

      REAL                                                              &
                           !,intent(out):
     & u_s_d  (row_length,  rows)                                       &
                                     ! u_s  diag at theta pts
     &,v_s_d  (row_length,  rows)    ! v_s  diag at theta pts

      REAL                                                              &
                           !,intent(out):
     & nsq_s_d(row_length,  rows)    ! nsq_s diag at theta pts

      REAL                                                              &
                           !,intent(out):
     & fr_d   (row_length,  rows)    ! Fr diag at theta pts

      REAL                                                              &
                           !,intent(out):
     & bld_d  (row_length,  rows)    ! blocked layer depth at theta pts

      REAL                                                              &
                           !,intent(out):
     & bldt_d (row_length,  rows)    ! % of time blocked layer diagnosed

      REAL                                                              &
                           !,intent(out):
     & num_lim_d (row_length,  rows) ! % of time numerical
                                     ! limiter invoked

      REAL                                                              &
                           !,intent(out):
     & num_fac_d (row_length,  rows) ! % redn. of flow-blocking stress
                                     ! after numerical limiter invoked

      REAL                                                              &
                           !,intent(out):
     & tausx_d(row_length, u_rows)                                       &
                                   !x-component of surface stress
     &,tausy_d(row_length, v_rows)  !y-component of surface stress

      REAL                                                              &
                           !,intent(out):
     & taus_scale_d(row_length,rows) ! Factor surface stress scaled by
!                                    ! if Froude no. dependence is on.

      REAL                                                              &
                           !,intent(out):
     & orog_slope_d(row_length,rows)                                    & 
                           ! orog slope (unavailable: 5a only)
     &,orog_anis_d(row_length,rows)                                     & 
                           ! orog anisotropy (unavailable: 5a only)
     &,orog_dir_d(row_length,rows)                                   
                           ! orog orientation (unavailable:5a only)


      LOGICAL                                                           &
                           !,intent(in):
     &  at_extremity(4)    ! indicates if this processor is at north,
                           ! south, east or west of the processor grid

      LOGICAL                                                           &
                           !,intent(in):
     &  l_taus_scale                                                    &
                           ! if true then surface stress is made to
!                          ! depend on the low level Froude number
     &, l_fix_gwsatn                                                    &
!                          ! if true then invoke minor bug fixes in     
!                          ! gwsatn
     &, l_gwd_40km         ! if true then don't apply GWD above 40km 

      INTEGER                                                           &
                            !,intent(in):
     & sat_scheme           ! Switch to determine whether to use
!                           ! amplitude or stress based saturation test


!
! END of full field diagnostic arrays
! Below are the stash flags for calculating diagnostics:
!
      LOGICAL                                                           &
                           !,intent(in):
     & stress_ud_on                                                     &
                           !u      stress
     &,stress_ud_p_on                                                   &
                           !u on press      stress
     &,stress_vd_on                                                     &
                           !v      stress
     &,stress_ud_satn_on                                                &
                           !u satn stress
     &,stress_vd_satn_on                                                &
                           !v satn stress
     &,stress_ud_wake_on                                                &
                           !u wake stress
     &,stress_vd_wake_on                                                &
                           !v wake stress
     &,du_dt_satn_on                                                    &
                           !u accel (saturation)
     &,du_dt_satn_p_on                                                  &
                           !u accel (saturation)
     &,dv_dt_satn_on                                                    &
                           !v accel (saturation)
     &,du_dt_wake_on                                                    &
                           !u accel blocked flow
     &,dv_dt_wake_on                                                    &
                           !v accel blocked flow
     &,u_s_d_on                                                         &
                           !u_s_d   diag switch
     &,v_s_d_on                                                         &
                           !v_s_d   diag switch
     &,nsq_s_d_on                                                       &
                           !nsq_s_d diag switch
     &,fr_d_on                                                          &
                           !fr_d    switch
     &,bld_d_on                                                         &
                           !bld_d   switch
     &,bldt_d_on                                                        &
                           !bldt_d   switch
     &,num_lim_d_on                                                     &
                           !num_lim_d switch
     &,num_fac_d_on                                                     &
                           !num_fac_d switch
     &,tausx_d_on                                                       &
                           !tausx_d switch
     &,tausy_d_on                                                       &
                           !tausy_d switch
     &,taus_scale_d_on                                                  &
                           !taus_scale_d switch
     &,orog_slope_d_on                                                  &
                           !orog_slope_d switch (5a scheme only)
     &,orog_anis_d_on                                                   &
                           !orog_anis_d switch (5a scheme only)
     &,orog_dir_d_on       !orog_dir_d switch (5a scheme only)

!-------------------------------------------------------
! New variables for the 5A scheme(not used in this scheme) 
!-------------------------------------------------------                     
  LOGICAL   ::          &!intent(in):
    l_dynbeta,          &!if true -dynamically adjusting beta (group vel angle)
    l_nonhydro,         &!if true - use nonhydro scheme
    l_smooth,           &!if true - lambda_z smoothing of acc
    l_gw_heating(3)      !if true - calculate heating tendency       


  REAL ::       &!block for intent(in) 
    Gsharp,     &!function of mtn sharpness
    fbcd         !flow blocking drag coefficient

!--------------------------------------------------------------------
! LOCAL DYNAMIC ARRAYS:
!--------------------------------------------------------------------

      INTEGER                                                           &
     & k_top(land_points)                                               &
                                     ! model level at mountain tops -
!                                    ! exact definition given in gwsurf
     &,k_top_max                     ! max(k_top)

! parameters for MPP
      INTEGER                                                           &
     &   pnorth,                                                        &
                      ! north processor address in the neighbor array
     &   peast,                                                         &
                      ! east processor address in the neighbor array
     &   psouth,                                                        &
                      ! south processor address in the neighbor array
     &   pwest,                                                         &
                      ! west processor address in the neighbor array
     &   nodomain     ! value in neighbor array if the domain has
                      !  no neighbor in this direction. otherwise
                      !  the value will be the tid of the neighbor
      PARAMETER (                                                       &
     &   pnorth   = 1,                                                  &
     &   peast    = 2,                                                  &
     &   psouth   = 3,                                                  &
     &   pwest    = 4,                                                  &
     &   nodomain = -1)

      INTEGER                                                           &
     & i,j,k,l                 ! loop counters in routine

! Work arrays
      REAL      ::                                &!block for work arrays
     & work_u(row_length,rows,levels)             &       
     &,work_v(row_length,rows,levels)             &      
     &,work_halo(wdims_s%i_start:wdims_s%i_end,   & 
     &         wdims_s%j_start:wdims_s%j_end,     &
     &         wdims_s%k_start:wdims_s%k_end),    &
     & work_on_u_grid(udims%i_start:udims%i_end,  & 
     &                udims%j_start:udims%j_end,  &
     &                udims%k_start:udims%k_end), &
     & work_on_v_grid(vdims%i_start:vdims%i_end,  & 
     &                vdims%j_start:vdims%j_end,  &
     &                vdims%k_start:vdims%k_end) 

      REAL                                                              &
     & up_land(land_points,levels)                                      &
                                     ! interpolated u on theta grid
     &,vp_land(land_points,levels)   ! interpolated V on theta grid

      REAL                                                              &
     & theta_land(land_points,levels)! land theta field on theta levels

      REAL                                                              &
     & r_rho_levels_land(land_points,levels)
!                                    !  land field of heights of
!                                    !  rho levels above z=0

      REAL                                                              &
     & r_theta_levels_land(land_points,levels)!)
!                                    !  land field of heights of
!                                    !  theta levels above z=0

      REAL                                                              &
     & rho_land(land_points,levels)  ! density at land points

      REAL                                                              &
     & s_x_lin_stress(land_points)                                      & 
                                     ! 'surface'  x_lin_stress land pnts
     &,s_y_lin_stress(land_points)   ! 'surface'  y_lin_stress land pnts

      REAL                                                              &
     & s_x_wake_stress(land_points)                                     & 
                                     ! 'surface' x_wake_stress land pts
     &,s_y_wake_stress(land_points)  ! 'surface' y_wake_stress land pts

      REAL                                                              &
     & s_x_orog(land_points)                                            & 
                                     ! 'surface' x_orog on land points
     &,s_y_orog(land_points)         ! 'surface' y_orog on land points

      REAL                                                              &
     & du_dt(land_points,levels)                                        &
                                     ! total GWD du/dt on land/theta
     &,dv_dt(land_points,levels)     ! total GWD dv/dt on land/theta

      REAL                                                              &
     & lift(land_points)             ! depth of blocked layer

      REAL                                                              &
     & fr(land_points)               ! low level froude number

      REAL                                                              &
     & rho_s(land_points)            ! low level density

!
! Land points arrays below are for the total GWD stress and
! its 4 individual components. (x and y components for each)
!
      REAL                                                              &
     & stress_ud_land      ( points_stress_ud      , 0:levels )         &
     &,stress_vd_land      ( points_stress_vd      , 0:levels )         &
     &,stress_ud_satn_land ( points_stress_ud_satn , 0:levels )         &
     &,stress_vd_satn_land ( points_stress_vd_satn , 0:levels )         &
     &,stress_ud_wake_land ( points_stress_ud_wake , 0:levels )         &
     &,stress_vd_wake_land ( points_stress_vd_wake , 0:levels )
!
!  Land point arrays below are for the 4 individual components of
!  the GWD wind increment.  (x and y components for each)
!
      REAL                                                              &
     & du_dt_satn_land ( points_du_dt_satn , levels )                   &
     &,dv_dt_satn_land ( points_dv_dt_satn , levels )                   &
     &,du_dt_wake_land ( points_du_dt_wake , levels )                   &
     &,dv_dt_wake_land ( points_dv_dt_wake , levels )

!
!  Land point arrays below are for the 9 GWD 'surface' diagnostics.
!
      REAL                                                              &
     & u_s_d_land     ( points_u_s_d     )                              &
     &,v_s_d_land     ( points_v_s_d     )                              &
     &,nsq_s_d_land   ( points_nsq_s_d   )                              &
     &,fr_d_land      ( points_fr_d      )                              &
     &,bld_d_land     ( points_bld_d     )                              &
     &,bldt_d_land    ( points_bldt_d    )                              &
     &,num_lim_d_land ( points_num_lim_d )                              &
     &,num_fac_d_land ( points_num_fac_d )                              &
     &,tausx_d_land   ( points_tausx_d   )                              &
     &,tausy_d_land   ( points_tausy_d   )                              &
     &,taus_scale_d_land ( points_taus_scale_d )

      REAL                                                              &
     & recip_a2       ! 1/(radius of earth)^2

      PARAMETER                                                         &
     & (recip_a2=1./(earth_radius*earth_radius) )

      LOGICAL                                                           &
     & l_drag(land_points)           ! whether point has a non-zero
!                                    ! stress or not

      TYPE (array_dims)  :: gwd_udims  &  ! gwd 0 to top udims for diags
                           ,gwd_vdims  &  ! gwd 0 to top vdims for diags
                           ,gwd_udims1 &  ! gwd single level udims
                           ,gwd_vdims1    ! gwd single level vdims


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
!
!------------------------------------------------------------------
!l    1.1 interpolate winds to p/theta-grid
!------------------------------------------------------------------

      IF (lhook) CALL dr_hook('G_WAVE_4A',zhook_in,zhook_handle)
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

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(up_land,       &
!$OMP& vp_land, work_u, work_v, land_index, land_points, levels,        &
!$OMP& row_length) PRIVATE(k,l,i,j)
      Do k=1,levels
        Do l=1,land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          up_land(l,k) =work_u(i,j,k)
          vp_land(l,k) =work_v(i,j,k)
        End do
      End do
!$OMP END PARALLEL DO


!------------------------------------------------------------------
!l    1.3  gather theta, rho and heights at land points
!------------------------------------------------------------------

      If (land_points  >   0) Then

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(rho_land,      &
!$OMP& r_rho_levels_land, r_theta_levels_land, theta_land, theta, rho,  &
!$OMP& r_theta_levels, row_length, land_index, land_points, levels,     &
!$OMP& r_rho_levels) PRIVATE(i,j,k,l)
       Do k=1,levels
        Do l=1,land_points
          j=(land_index(l)-1)/row_length + 1
          i=land_index(l) - (j-1)*row_length
          r_rho_levels_land(l,k)  = r_rho_levels(i,j,k) -               &
     &                              r_theta_levels(i,j,0)
          r_theta_levels_land(l,k)= r_theta_levels(i,j,k) -             &
     &                              r_theta_levels(i,j,0)
          rho_land(l,k)           = rho(i,j,k)*recip_a2
          theta_land(l,k)         = theta(i,j,k)
        End do
      End do
!$OMP END PARALLEL DO

!------------------------------------------------------------------
!l    2. calculate anisotropic 'surface' stress,CALL gw_surf
!------------------------------------------------------------------

! DEPENDS ON: gw_surf
      CALL gw_surf(                                                     &
     &            r_theta_levels_land,                                  &
     &            rho_land,theta_land,up_land,vp_land,timestep,         &
     &            sd_orog_land,orog_grad_xx_land,orog_grad_xy_land,     &
     &            orog_grad_yy_land,s_x_lin_stress,s_y_lin_stress,      &
     &            s_x_wake_stress,s_y_wake_stress,                      &
     &            s_x_orog,s_y_orog,levels,land_points,kay,rho_s,       &
     &            l_taus_scale, k_top,k_top_max,lift,l_drag,fr,frc,     &
     &            u_s_d_land     , u_s_d_on     , points_u_s_d    ,     &
     &            v_s_d_land     , v_s_d_on     , points_v_s_d    ,     &
     &            nsq_s_d_land   , nsq_s_d_on   , points_nsq_s_d  ,     &
     &            fr_d_land      , fr_d_on      , points_fr_d     ,     &
     &            bld_d_land     , bld_d_on     , points_bld_d    ,     &
     &            bldt_d_land    , bldt_d_on    , points_bldt_d   ,     &
     &            num_lim_d_land , num_lim_d_on , points_num_lim_d,     &
     &            num_fac_d_land , num_fac_d_on , points_num_fac_d,     &
     &            tausx_d_land   , tausx_d_on   , points_tausx_d  ,     &
     &            tausy_d_land   , tausy_d_on   , points_tausy_d  ,     &
     &            taus_scale_d_land    , taus_scale_d_on          ,     &
     &                                    points_taus_scale_d      )


!------------------------------------------------------------------
!l    3. calculate stress profile and accelerations,
!l       CALL gw_vert
!------------------------------------------------------------------

! DEPENDS ON: gw_vert
      CALL gw_vert(                                                     &
     &   rho_land,r_rho_levels_land,r_theta_levels_land,                &
     &   theta_land,up_land,vp_land,levels,land_points,                 &
     &   kay,sd_orog_land,s_x_lin_stress,s_y_lin_stress,                &
     &   s_x_wake_stress,s_y_wake_stress,                               &
     &   s_x_orog,s_y_orog,du_dt,dv_dt,                                 &
     &   k_top,k_top_max,lift,l_drag,fr,rho_s,l_fix_gwsatn,l_gwd_40km,  &
     &   sat_scheme,fsat,                                               &
     &   stress_ud_land ,points_stress_ud ,stress_ud_on,stress_ud_p_on, &
     &   stress_vd_land ,points_stress_vd ,stress_vd_on,                &
     &   stress_ud_satn_land,points_stress_ud_satn,stress_ud_satn_on,   &
     &   stress_vd_satn_land,points_stress_vd_satn,stress_vd_satn_on,   &
     &   stress_ud_wake_land,points_stress_ud_wake,stress_ud_wake_on,   &
     &   stress_vd_wake_land,points_stress_vd_wake,stress_vd_wake_on,   &
     &  du_dt_satn_land,points_du_dt_satn,du_dt_satn_on,du_dt_satn_p_on,&
     &   dv_dt_satn_land,points_dv_dt_satn,dv_dt_satn_on,               &
     &   du_dt_wake_land,points_du_dt_wake,du_dt_wake_on,               &
     &   dv_dt_wake_land,points_dv_dt_wake,dv_dt_wake_on )

      End If ! on land_points > 0

!------------------------------------------------------------------
!l    4. scatter accelerations to full area, interpolate to uv-grid
!l       and update winds
!------------------------------------------------------------------

! initialise work array: required to ensure zero values generated by
!   non-land points. [note that it is not necessary to re-initialise
!   the work array for each subsequent diagnostic since only land mask
!   points will differ, and these will be overwritten.]

!$OMP  PARALLEL DEFAULT(NONE) SHARED(work_halo, du_dt, row_length,      &
!$OMP& levels, wdims_s, land_index, land_points) PRIVATE(i,j,k,l)

!$OMP DO SCHEDULE(STATIC)
      Do k=0,levels
       Do j=wdims_s%j_start,wdims_s%j_end
        Do i=wdims_s%i_start,wdims_s%i_end
          work_halo(i,j,k)=0.
        End do !i
       End do  !j
      End do   !k
!$OMP END DO

! expand from land points and interpolate to 'c' u,v grid

!$OMP DO SCHEDULE(STATIC)
      Do k=1,levels

        Do l=1,land_points
           j=(land_index(l)-1)/row_length + 1
           i=land_index(l) - (j-1)*row_length

           work_halo(i,j,k) = du_dt(l,k)

        End do ! l

      End do ! k
!$OMP END DO

!$OMP END PARALLEL

! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,1),                &
     &     row_length, rows,                                            &
     &     levels, off_x, off_y, fld_type_p, .false.)


      CALL p_to_u(work_halo(pdims_s%i_start,pdims_s%j_start,1),   &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   1,levels,work_on_u_grid)


!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(land_points,   &
!$OMP& land_index,work_halo,row_length,dv_dt,levels) PRIVATE(i,j,k,l)
      Do k=1,levels

        Do l=1,land_points
           j=(land_index(l)-1)/row_length + 1
           i=land_index(l) - (j-1)*row_length

           work_halo(i,j,k) = dv_dt(l,k)

        End do ! l

      End do ! k
!$OMP END PARALLEL DO 


! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,1),                &
     &     row_length, rows,                                            &
     &     levels, off_x, off_y, fld_type_p, .false.)


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

!$OMP  PARALLEL DEFAULT(NONE) SHARED(r_u,r_v,timestep,work_on_u_grid,   &
!$OMP& levels, vdims, udims, work_on_v_grid) PRIVATE(k,j,i)

!$OMP DO SCHEDULE(STATIC)
      Do k=1,levels
        Do j=udims%j_start,udims%j_end
          Do i=udims%i_start,udims%i_end
            r_u(i,j,k)=r_u(i,j,k)+ timestep*work_on_u_grid(i,j,k)
          End do
        End do
      End do
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
      Do k=1,levels
        Do j=vdims%j_start,vdims%j_end
          Do i=vdims%i_start,vdims%i_end
            r_v(i,j,k)=r_v(i,j,k)+ timestep*work_on_v_grid(i,j,k)
          End do
        End do
      End do
!$OMP END DO

!$OMP END PARALLEL 

!------------------------------------------------------------------
!l    5. gather diagnostics from land points and interpolate from p to
!l       u,v staggering on 'c' grid.
!l       note that stress_ud_on,stress_vd_on are on theta levels
!l       whereas remaining diagnostics are on rho levels.
!------------------------------------------------------------------


      If (stress_ud_on .or. stress_ud_p_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=0,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = stress_ud_land(l,k)

          End do ! l

        End do ! k


! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
     &     row_length, rows,                                            &
     &     levels+1, off_x, off_y, fld_type_p, .false.)


      CALL p_to_u(work_halo,                                      &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   0,levels,stress_ud)

      EndIf ! on stashflag

      If (stress_vd_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=0,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = stress_vd_land(l,k)

          End do ! l

        End do ! k


! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
     &     row_length, rows,                                            &
     &     levels+1, off_x, off_y, fld_type_p, .false.)


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

      EndIf ! on stashflag


      If (stress_ud_satn_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid

        Do k=0,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = stress_ud_satn_land(l,k)

          End do ! l

        End do ! k


! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
     &     row_length, rows,                                            &
     &     levels+1, off_x, off_y, fld_type_p, .false.)


      CALL p_to_u(work_halo,                                      &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   0,levels,stress_ud_satn)


      EndIf ! on stashflag

      If (stress_vd_satn_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=0,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = stress_vd_satn_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
     &     row_length, rows,                                            &
     &     levels+1, off_x, off_y, fld_type_p, .false.)


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

      EndIf ! on stashflag


      If (stress_ud_wake_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=0,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = stress_ud_wake_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
     &     row_length, rows,                                            &
     &     levels+1, off_x, off_y, fld_type_p, .false.)


      CALL p_to_u(work_halo,                                      &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   0,levels,stress_ud_wake)

      EndIf ! on stashflag

      If (stress_vd_wake_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=0,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = stress_vd_wake_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
     &     row_length, rows,                                            &
     &     levels+1, off_x, off_y, fld_type_p, .false.)


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

      EndIf ! on stashflag


      If (du_dt_satn_on .or. du_dt_satn_p_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=1,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = du_dt_satn_land(l,k)

          End do ! l

        End do ! k


! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,1),                &
     &     row_length, rows,                                            &
     &     levels, off_x, off_y, fld_type_p, .false.)


      CALL p_to_u(work_halo(pdims_s%i_start,pdims_s%j_start,1),   & 
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   1,levels,du_dt_satn)


      EndIf ! on stashflag

      If (dv_dt_satn_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=1,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = dv_dt_satn_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,1),                &
     &     row_length, rows,                                            &
     &     levels, off_x, off_y, fld_type_p, .false.)


      CALL p_to_v(                                                &
              work_halo(pdims_s%i_start,pdims_s%j_start,1),       &
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

      EndIf ! on stashflag


      If (du_dt_wake_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=1,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = du_dt_wake_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,1),                &
     &     row_length, rows,                                            &
     &     levels, off_x, off_y, fld_type_p, .false.)


      CALL p_to_u(work_halo(pdims_s%i_start,pdims_s%j_start,1),   &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   1,levels,du_dt_wake)

      EndIf ! on stashflag


      If (dv_dt_wake_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid
        Do k=1,levels

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,k) = dv_dt_wake_land(l,k)

          End do ! l

        End do ! k

! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,1),                &
     &     row_length, rows,                                            &
     &     levels, off_x, off_y, fld_type_p, .false.)


      CALL p_to_v(                                                &
                  work_halo(pdims_s%i_start,pdims_s%j_start,1),   &
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
      EndIf ! on stashflag


      If (u_s_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             u_s_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            u_s_d(i,j) = u_s_d_land(l)
         End Do ! l

      EndIf ! on stashflag

      If (v_s_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             v_s_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            v_s_d(i,j) = v_s_d_land(l)
         End Do ! l

      EndIf ! on stashflag

      If (nsq_s_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             nsq_s_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            nsq_s_d(i,j) = nsq_s_d_land(l)
         End Do ! l

      EndIf ! on stashflag


      If (fr_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             fr_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            fr_d(i,j) = fr_d_land(l)
         End Do ! l

      EndIf ! on stashflag

      If (bld_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             bld_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            bld_d(i,j) = bld_d_land(l)
         End Do ! l

      EndIf ! on stashflag


      If (bldt_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             bldt_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            bldt_d(i,j) = bldt_d_land(l)
         End Do ! l

      EndIf ! on stashflag


      If (num_lim_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             num_lim_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            num_lim_d(i,j) = num_lim_d_land(l)
         End Do ! l

      EndIf ! on stashflag


      If (num_fac_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             num_fac_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            num_fac_d(i,j) = num_fac_d_land(l)
         End Do ! l

      EndIf ! on stashflag

      If (tausx_d_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid

        Do l=1,land_points
           j=(land_index(l)-1)/row_length + 1
           i=land_index(l) - (j-1)*row_length

           work_halo(i,j,0) = tausx_d_land(l)

        End do ! l

! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
     &     row_length, rows,                                            &
     &     1, off_x, off_y, fld_type_p, .false.)


      CALL p_to_u(work_halo(pdims_s%i_start,pdims_s%j_start,0),   &
                   pdims_s%i_start,pdims_s%i_end,                 &
                   pdims_s%j_start,pdims_s%j_end,                 &
                   udims%i_start,udims%i_end,                     &
                   udims%j_start,udims%j_end,                     &
                   1,1,tausx_d(udims%i_start,udims%j_start))

      EndIf ! on stashflag

      If (tausy_d_on) Then  ! stashflag for this diagnostic

! expand from land points and interpolate to 'c' u,v grid

          Do l=1,land_points
             j=(land_index(l)-1)/row_length + 1
             i=land_index(l) - (j-1)*row_length

             work_halo(i,j,0) = tausy_d_land(l)

          End do ! l

! DEPENDS ON: swap_bounds
      CALL swap_bounds(                                                 &
     &     work_halo(pdims_s%i_start,pdims_s%j_start,0),                &
     &     row_length, rows,                                            &
     &     1, off_x, off_y, fld_type_p, .false.)


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

      EndIf ! on stashflag

      If (taus_scale_d_on) Then  ! stashflag for this diagnostic

! expand from land points

         Do j=1,rows
           Do i=1,row_length
             taus_scale_d(i,j) = 0.0
           End do
         End do

         Do l=1,land_points
            j=(land_index(l)-1)/row_length + 1
            i=land_index(l) - (j-1)*row_length
            taus_scale_d(i,j) = taus_scale_d_land(l)
         End Do ! l

      EndIf ! on stashflag


      iret=0

      IF (lhook) CALL dr_hook('G_WAVE_4A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE g_wave_4a
END MODULE
