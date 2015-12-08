! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Declares variables for turbulent diffusion input 
!
!  Code Owner: See Unified Model Code Owner's HTML page
!  This file belongs in section: Top Level

MODULE turb_diff_mod

! Subgrid turbulence scheme input.

USE atmos_max_sizes, ONLY : model_levels_max
USE Control_Max_Sizes, ONLY : max_121_rows, max_sponge_width,           &
    max_updiff_levels

IMPLICIT NONE

LOGICAL :: L_diff_active
LOGICAL :: L_subfilter_horiz
LOGICAL :: L_subfilter_vert
LOGICAL :: L_subfilter_blend  ! blend diffusion coeffs

INTEGER :: turb_startlev_horiz   ! 1st lev for horiz subgrid turb
INTEGER :: turb_endlev_horiz     ! last lev for horiz subgrid turb
INTEGER :: turb_startlev_vert    ! 1st lev for vert subgrid turb
INTEGER :: turb_endlev_vert      ! last lev for vert subgrid turb

REAL :: diff_factor
REAL :: mix_factor

      LOGICAL :: L_print_w
      LOGICAL :: L_print_div
      LOGICAL :: L_diag_print_ops     ! diagnostic prints for ops
      LOGICAL :: L_print_pe     ! print diagnostics on all pe's if true
      LOGICAL :: L_print_shear        ! wind diagnostic prints
      LOGICAL :: L_print_max_wind     ! wind diagnostic prints
      LOGICAL :: L_diag_L2norms  ! l2norm diagnostic prints
      LOGICAL :: L_diag_L2helm   ! l2norm diagnostic prints from solver
      LOGICAL :: L_diag_wind     ! wind diagnostic prints
      LOGICAL :: L_diag_noise    ! w diagnostic prints
      LOGICAL :: L_flush6        ! if T then flush buffers on failure
      LOGICAL :: L_diag_print  ! Print diagnostics
      LOGICAL :: L_print_lapse ! Print lapse_rate diagnostics
      LOGICAL :: L_print_wmax ! Print max w diagnostic
      LOGICAL :: L_print_theta1 ! Print level 1 theta diagnostic

      LOGICAL :: L_diffusion
      LOGICAL :: L_upper_ramp   ! ramp upper-level diffusion
      LOGICAL :: L_adjust_theta ! activate convective adjustment
      LOGICAL :: L_vdiff_uv    ! activate targeted diffusion uv
      LOGICAL :: L_filter    ! activate polar filter or diffusion
      LOGICAL :: L_filter_incs    ! activate polar filter of incs
      LOGICAL :: L_pftheta   ! activate polar filter of theta
      LOGICAL :: L_pfuv      ! activate polar filter of u,v
      LOGICAL :: L_pfw       ! activate polar filter of w
      LOGICAL :: L_pofil_new  ! activate new polar filter
      LOGICAL :: L_pfexner  ! activate polar filter of Exner pressure
      LOGICAL :: L_pofil_hadgem2 ! run with HadGEM2 polar filter setting   
      LOGICAL :: L_diff_exner  ! activate diffusion of Exner pressure
      LOGICAL :: L_pfcomb   ! combined polar filter/diffusion active
      LOGICAL :: L_sponge   ! activate lateral boundariessponge zones
      LOGICAL :: L_pfincs    ! activate polar filter of incs
      LOGICAL :: L_diff_thermo ! horiz. diffusion of theta
      LOGICAL :: L_diff_wind   ! horiz. diffusion of u, v
      LOGICAL :: L_diff_w      ! horiz. diffusion of w
      LOGICAL :: L_diff_incs   ! horiz. diffusion of increments
      LOGICAL :: L_diff_auto   ! UM calculates diffusion parameters
      LOGICAL :: L_tardiff_q     ! activate targeted diffusion q
      LOGICAL :: L_diff_ctl      ! general diffusion control
! L_diff_ctl == L_diffusion .or. L_cdiffusion .or. L_vertical_diffusion
!               .or. L_divdamp .or.  L_tardiff_q .or. L_diag_print
      LOGICAL :: L_cdiffusion

      INTEGER :: print_step    ! To control diagnostic printing interval
      INTEGER :: diag_interval ! diagnostic printing sampling frequency
      INTEGER :: norm_lev_start ! start level for norm diagnostics
      INTEGER :: norm_lev_end   ! end level for norm diagnostics
      INTEGER :: first_norm_print ! first timestep for norm printing
      INTEGER :: dom_w_in    ! define start for block printing
      INTEGER :: dom_e_in    ! define end for block printing
      INTEGER :: dom_s_in    ! define start for block printing
      INTEGER :: dom_n_in    ! define end for block printing
      INTEGER :: blockx_in    ! define size for block printing
      INTEGER :: blocky_in    ! define size for  block printing

      INTEGER  :: diffusion_order_thermo(model_levels_max)
      INTEGER  :: diffusion_order_wind(model_levels_max)
      INTEGER  :: diffusion_order_w(model_levels_max)
      INTEGER  :: diffusion_order_q(model_levels_max)
      INTEGER :: u_begin(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: u_end(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: v_begin(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: v_end(0:max_121_rows) ! Sweep control on 121 filter
      INTEGER :: u_sweeps(max_121_rows) ! Sweep control on 121 filter
      INTEGER :: v_sweeps(max_121_rows) ! Sweep control on 121 filter
      INTEGER :: max_sweeps ! Max sweeps wanted for 121 filter
      INTEGER :: global_u_filter ! Sweep control on 121 filter
      INTEGER :: global_v_filter ! Sweep control on 121 filter
      INTEGER :: diff_order_thermo   ! diffusion order for theta
      INTEGER :: diff_order_wind ! diffusion order for winds
      INTEGER :: diff_timescale_thermo  ! diffusion timescale for theta
      INTEGER :: diff_timescale_wind    ! diffusion timescale for wind
      INTEGER :: vdiffuv_timescale ! diffusion e-folding timesteps
      INTEGER :: vdiffuv_start ! start level  targeted diffusion
      INTEGER :: vdiffuv_end ! end level targeted diffusion
      INTEGER :: adjust_theta_start ! start level convective adjustment
      INTEGER :: adjust_theta_end   ! end levelconvective adjustment
      INTEGER :: top_filt_start ! start level upper-level diffusion
      INTEGER :: top_filt_end   ! end level upper-level diffusion
      INTEGER :: sponge_ew   ! left/right boundaries sponge zone width
      INTEGER :: sponge_ns   ! north/south boundaries sponge zone width
      INTEGER :: sponge_power ! sponge zone weighting order
      INTEGER :: tardiffq_test ! test level test w targetted diffusion
      INTEGER :: tardiffq_start ! start level test w targetted diffusion
      INTEGER :: tardiffq_end ! end level test w targetted diffusion

      !  level - assume surfaces are horizontal
      INTEGER  :: horizontal_level
      INTEGER  :: tar_horizontal  ! steep slope test targeted diffusion

      REAL :: diffusion_coefficient_thermo(model_levels_max)
      REAL :: diffusion_coefficient_wind(model_levels_max)
      REAL :: diffusion_coefficient_w(model_levels_max)
      REAL :: diffusion_coefficient_q(model_levels_max)
      REAL :: vdiffuv_test ! test to activate shear diffusion of u,v
      REAL :: vdiffuv_factor ! vertical diffusion coeff for shear diff
      REAL :: scale_ratio ! Pass control on 121 filter
      REAL :: ref_lat_deg ! Reference lat for auto diffusion
      REAL :: top_diff    !  upper-level diffusion coefficient
      REAL :: up_diff_scale    !  upper-level diffusion ramping factor
      REAL :: adjust_lapse_min !  min dtheta/dz in vertical adjustment

      REAL :: up_diff(max_updiff_levels) ! upper-level diffusion coeff
      REAL :: sponge_wts_ew(max_sponge_width) ! sponge weights
      REAL :: sponge_wts_ns(max_sponge_width) ! sponge weights

      REAL :: diff_coeff_ref ! EW diffusion coefficient at polar cap
      REAL :: diff_coeff_thermo  ! NS theta diffusion coeff
      REAL :: diff_coeff_wind    ! NS u,v diffusion coeff
      REAL :: diff_coeff_phi    ! North-South diffusion coefficient
!   reference latitudes for filtering and diffusion
      REAL :: polar_cap  ! Apply 1-2-1 filter polewards
      REAL :: tardiffq_factor ! targeted diffusion coefficient
      REAL :: w_print_limit ! w Threshold for diagnostic printing
      REAL :: w_conv_limit  ! w Threshold for limiting convection

      ! Divergence damping control options:
      LOGICAL :: L_divdamp

      REAL :: div_damp_coefficient(model_levels_max)

      ! Polar filter control options:
      LOGICAL :: L_polar_filter
      ! T: use polar filter to filter increment
      LOGICAL :: L_polar_filter_incs

      REAL :: polar_filter_north_lat_limit
      REAL :: polar_filter_south_lat_limit
      REAL :: polar_filter_coefficient

      ! amount in radians to increment start latitude by per sweep
      REAL :: polar_filter_step_per_sweep

      ! max latitude at which filter can star
      REAL :: polar_filter_lat_limit

      ! number of sweeps of filter to do
      INTEGER :: polar_filter_n_sweeps

      ! Vertical Diffusion control options:
      LOGICAL :: L_vertical_diffusion
      LOGICAL :: L_ramp

      INTEGER :: level_start_wind
      INTEGER :: level_stop_wind
      INTEGER :: level_start_theta
      INTEGER :: level_stop_theta
      INTEGER :: level_start_q
      INTEGER :: level_stop_q

      REAL :: vert_diffusion_coeff_wind
      REAL :: vert_diffusion_coeff_theta
      REAL :: vert_diffusion_coeff_q
      REAL :: ramp_lat_radians

      ! Moisture resetting control options:
      LOGICAL :: l_qpos              ! logical to run qpos code
      LOGICAL :: l_qpos_diag_pr      ! Diagnostic print switch

      INTEGER :: q_pos_method        ! Algorithm choice
      INTEGER :: q_pos_tracer_method ! Algorithm choice for tracers

      REAL :: qpos_diag_limit        ! Limit for diagnostic prints
      REAL :: qlimit                 ! lowest allowed value of q

NAMELIST/RUN_Diffusion/L_diffusion,diffusion_order_thermo,              &
        diffusion_order_wind,diffusion_order_w,diffusion_order_q,       &
        diffusion_coefficient_thermo,diffusion_coefficient_wind,        &
        diffusion_coefficient_w,diffusion_coefficient_q,                &
        diff_order_thermo, diff_timescale_thermo,                       &
        diff_order_wind, diff_timescale_wind,                           &
        horizontal_level, tar_horizontal,                               &
        L_cdiffusion,L_ramp, ramp_lat_radians,                          &
        L_divdamp,div_damp_coefficient,                                 &
        L_polar_filter,L_polar_filter_incs,                             &
        polar_filter_north_lat_limit,polar_filter_south_lat_limit,      &
        polar_filter_coefficient,polar_filter_step_per_sweep,           &
        polar_filter_lat_limit,polar_filter_n_sweeps,                   &
        l_qpos,  q_pos_method, q_pos_tracer_method, qlimit,             &
        l_qpos_diag_pr, qpos_diag_limit,                                &
        L_diff_ctl, L_tardiff_q, w_conv_limit,                          &
        tardiffq_factor, tardiffq_test, tardiffq_start, tardiffq_end,   &
        L_subfilter_horiz, L_subfilter_vert, L_subfilter_blend,         &
        diff_factor, mix_factor, turb_startlev_horiz, turb_endlev_horiz,&
        turb_startlev_vert, turb_endlev_vert,                           &
        L_vertical_diffusion,                                           &
        L_pftheta, L_pfuv, L_pfw, L_pfincs, L_pfexner, L_pofil_hadgem2, &
        L_diff_incs, L_diff_thermo, L_diff_wind, L_diff_w,              &
        level_start_wind, level_stop_wind,                              &
        level_start_theta, level_stop_theta,                            &
        level_start_q, level_stop_q,                                    &
        L_pofil_new, L_diff_auto,                                       &
        diff_coeff_ref, polar_cap, scale_ratio, ref_lat_deg,            &
        max_sweeps, L_upper_ramp, up_diff_scale, top_diff,              &
        top_filt_start, top_filt_end, L_vdiff_uv, vdiffuv_timescale,    &
        vdiffuv_test, vdiffuv_factor, vdiffuv_start, vdiffuv_end,       &
        L_adjust_theta, adjust_theta_start, adjust_theta_end,           &
        L_sponge, sponge_power, sponge_ew, sponge_ns,                   &
        vert_diffusion_coeff_wind, vert_diffusion_coeff_theta,          &
        vert_diffusion_coeff_q,                                         &
        L_diag_print, L_diag_print_ops, L_print_pe,                     &
        L_print_w, L_print_wmax, L_print_lapse, L_print_theta1,         &
        L_print_div, L_diag_wind, L_print_shear, L_print_max_wind,      &
        L_diag_noise, L_diag_L2norms, L_diag_L2helm,                    &
        norm_lev_start, norm_lev_end, first_norm_print,                 &
        dom_w_in, dom_e_in, dom_s_in, dom_n_in, blockx_in, blocky_in,   &
        print_step, diag_interval, w_print_limit, L_Flush6

END MODULE turb_diff_mod
