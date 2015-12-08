! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine NI_DIFF_CTL
!
      Subroutine NI_DIFF_CTL(                                           &
     &                     L_diffusion, L_cdiffusion,                   &
     &                     L_vertical_diffusion, L_divdamp,             &
     &                     L_ramp, ramp_lat_radians,                    &
     &                     L_Backwards, Ltimer,                         &
     &                     timestep, pos_timestep, neg_timestep,        &
     &                     theta, w, q, qcl,qcf,                        &
     &                     qcf2, qrain, qgraup,                         &
     &                     mix_v, mix_cl, mix_cf,                       &
     &                     mix_cf2, mix_rain, mix_graup,                &
     &                     L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,       &
     &                     u, v, rho, Exner, w_adv,                     &
     &                     r_theta_levels, r_rho_levels, r_at_u, r_at_v,&
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     sec_theta_latitude,                          &
     &                     cos_theta_latitude, sec_v_latitude,          &
     &                     FV_sec_theta_latitude, cos_v_latitude,       &
     &                     sin_theta_latitude, sin_v_latitude,          &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     at_extremity, gc_proc_row_group,             &
     &                     mype, nproc, nproc_x, nproc_y, neighbour,    &
     &                     delta_lambda, delta_phi, L_regular,          &
     &                     lambda_p, phi_p, lambda_u, phi_v,            &
     &                     recip_dlamp, recip_dphip,                    &
     &                     recip_dlamu, recip_dphiv,                    &
     &                     rows, row_length, n_rows,                    &
     &                     model_levels, bl_levels, model_domain,      &
     &                     global_row_length, global_rows,              &
     &                     diffusion_coefficient_thermo,                &
     &                     diffusion_coefficient_w,                     &
     &                     diffusion_coefficient_q,                     &
     &                     diffusion_coefficient_wind,                  &
     &                     diffusion_order_thermo,                      &
     &                     diffusion_order_w,                           &
     &                     diffusion_order_q,                           &
     &                     diffusion_order_wind,                        &
     &                     horizontal_level, tar_horizontal,            &
     &                     level_start_wind, level_stop_wind,           &
     &                     level_start_q, level_stop_q,                 &
     &                     level_start_theta, level_stop_theta,         &
     &                     L_tardiff_q, w_conv_limit, tardiffq_factor,  &
     &                     tardiffq_test, tardiffq_start, tardiffq_end, &
     &                     L_diag_w, w_local_mask,                      &
     &                     L_adjust_theta,                              &
     &                     adjust_theta_start, adjust_theta_end,        &
     &                     L_vdiff_uv, vdiffuv_test,                    &
     &                     vdiffuv_factor, vdiffuv_start, vdiffuv_end,  &
     &                     vert_diffusion_coeff_wind,                   &
     &                     vert_diffusion_coeff_q,                      &
     &                     vert_diffusion_coeff_theta,                  &
     &                     div_damp_coefficient,                        &
     &                     theta_star, R_w, q_star,                     &
     &                     qcl_star, qcf_star,                          &
     &                     qcf2_star, qrain_star, qgraup_star,          &
     &                     R_u, R_v,                                    &
     &                     L_mix_ratio)

! Purpose:
!          Subroutine to interface to diffusion code
!
! Method:
!          Is described in ;
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
      IMPLICIT NONE

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, bl_levels                                                      &
                         ! number of bl levels.
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, offx                                                            &
                     ! Size of small halo in i
     &, offy                                                            &
                     ! Size of small halo in j.
     &, mype                                                            &
                     ! processor Id
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, gc_proc_row_group                                               &
     &, nproc                                                           &
                  ! Total number of processors
     &, nproc_x                                                         &
                   ! Number of processors in longitude
     &, nproc_y                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)                                                    &
                             ! Array with the Ids of the four neighbours
                             ! in the horizontal plane
     &, i,j,k

      REAL                                                              &
     & timestep                                                         &
                            ! atmosphere model timestep
     &,pos_timestep                                                     &
                    ! = +timestep.
     &,neg_timestep                                                     &
                    ! = -timestep.
     &,  w_conv_limit                                                   &
     &, vdiffuv_test 

      Real, Intent (InOut) ::                                           &
                              ! primary model variables
     &  u(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
     &, v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)     &
     &, w(1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels)     &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &          0:model_levels)                                         &
     &, rho(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)     &
     &, Exner(1-offx:row_length+offx, 1-offy:rows+offy,                 &
     &        model_levels + 1)                                         &
     &, q (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &     model_levels)                                                  &
     &, qcl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &     model_levels)                                                  &
     &, qcf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &     model_levels)                                                  &
     &, qcf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &     model_levels)                                                  &
     &, qrain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &     model_levels)                                                  &
     &, qgraup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &     model_levels)                                                  &
     &, theta (1-offx:row_length+offx, 1-offy:rows+offy,                &
     &         model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_at_u (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          model_levels)                                           &
     &, r_at_v (1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,     &
     &          model_levels)                                           &
     &, eta_theta_levels (0:model_levels)                               &
     &, eta_rho_levels (model_levels)                                   &
     &, cos_theta_latitude (1-offx:row_length+offx,                     &
     &                    1-offy:rows+offy)                             &
     &, sec_theta_latitude (1-offx:row_length+offx,                     &
     &                    1-offy:rows+offy)                             &
     &, sec_v_latitude (1-offx:row_length+offx, 1-offy:n_rows+offy)     &
     &, cos_v_latitude (1-offx:row_length+offx, 1-offy:n_rows+offy)     &
     &, FV_sec_theta_latitude (1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy)                        &
     &, sin_theta_latitude(row_length, rows)                            &
     &, sin_v_latitude(row_length, n_rows)

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, diffusion_order(model_levels)

      Real                                                              &
     &  diffusion_coefficient(model_levels)

      Logical                                                           &
     &  L_diffusion                                                     &
     &, L_cdiffusion                                                    &
     &, L_regular                                                       &
                        !  Variable resolution switch
     &, L_vertical_diffusion                                            &
     &, L_ramp                                                          &
     &, L_divdamp                                                       &
     &, L_Backwards                                                     &
     &, Ltimer                                                          &
     &, L_mix_ratio                                                     &
     &, L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup                           &
     &, L_vdiff_uv                                                      &
     &, L_adjust_theta                                                  &
     &, L_tardiff_q                                                     &
     &, L_diag_w

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
           !VarRes horizontal co-ordinate information
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, lambda_u(1-halo_i:row_length+halo_i)                            &
     &, phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)       &
     &, recip_dlamp(1-halo_i:row_length+halo_i)                         &
     &, recip_dphip(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)   &
     &, recip_dlamu(1-halo_i:row_length+halo_i)                         &
     &, recip_dphiv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

      INTEGER                                                           &
     & diffusion_order_thermo(model_levels)                             &
                                              ! value * del^2: order=2
     &,diffusion_order_wind(model_levels)                               &
                                              ! gives del^4 diffusion
     &,diffusion_order_w(model_levels-1)                                &
                                              !
     &,diffusion_order_q(model_levels)          !

      REAL                                                              &
     & diffusion_coefficient_thermo(model_levels)                       &
     &,diffusion_coefficient_wind(model_levels)                         &
     &,diffusion_coefficient_w(model_levels-1)                          &
     &,diffusion_coefficient_q(model_levels)                              &
     &, tardiffq_factor       ! targeted diffusion coefficient

      Integer                                                           &
     & horizontal_level                                                 &
     &, tar_horizontal    !  steep slope test targeted diffusion

      Integer                                                           &
     & level_start_wind, level_stop_wind                                &
     &,level_start_theta, level_stop_theta                              &
     &,level_start_q, level_stop_q                                      &
     &, adjust_theta_start                                              &
     &, adjust_theta_end                                                &
     &, vdiffuv_start                                                   &
     &, vdiffuv_end                                                     &
     &, tardiffq_test                                                   &
     &, tardiffq_start                                                  &
     &, tardiffq_end

      Real                                                              &
     & vert_diffusion_coeff_wind                                        &
     &,vert_diffusion_coeff_theta                                       &
     &,vert_diffusion_coeff_q                                           &
     &, ramp_lat_radians                                                &
     &, vdiffuv_factor


      Real                                                              &
     & div_damp_coefficient(model_levels)

!   variables with intent INOUT

      Real, Intent (InOut) ::                                           &
     &  R_u(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &        model_levels)                                             &
     &, R_v(1-offx:row_length+offx, 1-offy:n_rows+offy,                 &
     &        model_levels)                                             &
     &, R_w(row_length, rows, model_levels)                             &
     &, q_star(1-offx:row_length+offx,                                  &
     &           1-offy:rows+offy, model_levels)                          &
     &, qcl_star(1-offx:row_length+offx,                                &
     &           1-offy:rows+offy, model_levels)                          &
     &, qcf_star(1-offx:row_length+offx,                                &
     &           1-offy:rows+offy, model_levels)                          &
     &, qcf2_star(1-offx:row_length+offx,                               &
     &           1-offy:rows+offy, model_levels)                          &
     &, qrain_star(1-offx:row_length+offx,                              &
     &           1-offy:rows+offy, model_levels)                          &
     &, qgraup_star(1-offx:row_length+offx,                             &
     &           1-offy:rows+offy, model_levels)                          &
     &, theta_star(1-offx:row_length+offx,                              &
     &               1-offy:rows+offy, model_levels)                    &
     &, w_local_mask(row_length,rows)                                   &
     &, mix_v (1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, model_levels)                        &
     &, mix_cl(1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, model_levels)                        &
     &, mix_cf(1-halo_i:row_length+halo_i,                              &
     &         1-halo_j:rows+halo_j, model_levels)                        &
     &, mix_cf2(1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:rows+halo_j, model_levels)                       &
     &, mix_rain(1-halo_i:row_length+halo_i,                            &
     &           1-halo_j:rows+halo_j, model_levels)                      &
     &, mix_graup(1-halo_i:row_length+halo_i,                           &
     &            1-halo_j:rows+halo_j, model_levels)

!    Local arrays for when using mixing ratios
      REAL, DIMENSION (:,:,:), ALLOCATABLE ::                           &
     &  mix_v_star, mix_cl_star, mix_cf_star                            &
     &, mix_cf2_star, mix_rain_star, mix_graup_star

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



      IF (lhook) CALL dr_hook('NI_DIFF_CTL',zhook_in,zhook_handle)
      IF(L_mix_ratio)then
        allocate ( mix_v_star (1-offx:row_length+offx,                  &
     &                         1-offy:rows+offy,model_levels) )
        allocate ( mix_cl_star (1-offx:row_length+offx,                 &
     &                         1-offy:rows+offy,model_levels) )
        allocate ( mix_cf_star (1-offx:row_length+offx,                 &
     &                         1-offy:rows+offy,model_levels) )
        IF(L_mcr_qcf2)then
          allocate ( mix_cf2_star (1-offx:row_length+offx,              &
     &                             1-offy:rows+offy,model_levels) )
        else
          allocate ( mix_cf2_star (1,1,1) )
        endif
        IF(L_mcr_qrain)then
          allocate ( mix_rain_star(1-offx:row_length+offx,              &
     &                             1-offy:rows+offy,model_levels) )
        else
          allocate ( mix_rain_star (1,1,1) )
        endif
        IF(L_mcr_qgraup)then
          allocate ( mix_graup_star(1-offx:row_length+offx,             &
     &                             1-offy:rows+offy,model_levels) )
        else
          allocate ( mix_graup_star (1,1,1) )
        endif


!  convert q, qcl,qcf _star to mix_v, mix_cl,mix_cf _star
! q_star actually holds q_dash
! DEPENDS ON: q_to_mix
        call q_to_mix (row_length, rows, model_levels,                    &
     &                 offx,offy,                                       &
     &                 q_star, qcl_star, qcf_star,                      &
     &                 qcf2_star,qrain_star, qgraup_star,               &
     &                 L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
     &                 mix_v_star, mix_cl_star, mix_cf_star,            &
     &                 mix_cf2_star, mix_rain_star, mix_graup_star      &
     &                 )
! DEPENDS ON: diff_divdamp_ctl
          Call Diff_Divdamp_Ctl(                                        &
     &                     L_diffusion, L_cdiffusion,                   &
     &                     L_vertical_diffusion, L_divdamp,             &
     &                     L_ramp, ramp_lat_radians,                    &
     &                     L_Backwards,                                 &
     &                     timestep, pos_timestep, neg_timestep,        &
     &                     theta, w, mix_v,                             &
     &                     u, v, rho, Exner, w_adv,                     &
     &                     r_theta_levels, r_rho_levels, r_at_u, r_at_v,&
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     sec_theta_latitude,                          &
     &                     cos_theta_latitude, sec_v_latitude,          &
     &                     FV_sec_theta_latitude, cos_v_latitude,       &
     &                     sin_theta_latitude, sin_v_latitude,          &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     at_extremity, gc_proc_row_group,             &
     &                     mype, nproc, nproc_x, nproc_y, neighbour,    &
     &                     delta_lambda, delta_phi, L_regular,          &
     &                     lambda_p, phi_p, lambda_u, phi_v,            &
     &                     recip_dlamp, recip_dphip,                    &
     &                     recip_dlamu, recip_dphiv,                    &
     &                     rows, row_length, n_rows,                    &
     &                     model_levels, bl_levels, model_domain,      &
     &                     global_row_length,global_rows,               &
     &                     diffusion_coefficient_thermo,                &
     &                     diffusion_coefficient_w,                     &
     &                     diffusion_coefficient_q,                     &
     &                     diffusion_coefficient_wind,                  &
     &                     diffusion_order_thermo,                      &
     &                     diffusion_order_w,                           &
     &                     diffusion_order_q,                           &
     &                     diffusion_order_wind,                        &
     &                     horizontal_level, tar_horizontal,            &
     &                     level_start_wind, level_stop_wind,           &
     &                     level_start_q, level_stop_q,                 &
     &                     level_start_theta, level_stop_theta,         &
     &                     L_tardiff_q, w_conv_limit, tardiffq_factor,  &
     &                     tardiffq_test, tardiffq_start, tardiffq_end, &
     &                     L_diag_w, w_local_mask,                      &
     &                     L_adjust_theta,                              &
     &                     adjust_theta_start, adjust_theta_end,        &
     &                     L_vdiff_uv, vdiffuv_test,                    &
     &                     vdiffuv_factor, vdiffuv_start, vdiffuv_end,  &
     &                     vert_diffusion_coeff_wind,                   &
     &                     vert_diffusion_coeff_q,                      &
     &                     vert_diffusion_coeff_theta,                  &
     &                     div_damp_coefficient,                        &
     &                     theta_star, R_w, mix_v_star, R_u, R_v        &
     &                     )
!  convert mix_v, mix_cl,mix_cf to q, qcl,qcf

! DEPENDS ON: mix_to_q
        call mix_to_q                                                   &
     &                  (row_length, rows, model_levels,                  &
     &                   offx, offy,                                    &
     &                   mix_v_star, mix_cl_star, mix_cf_star,          &
     &                   mix_cf2_star, mix_rain_star, mix_graup_star,   &
     &                   L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,         &
     &                   q_star, qcl_star, qcf_star,                    &
     &                   qcf2_star,qrain_star, qgraup_star              &
     &                   )
        deallocate( mix_v_star)
        deallocate( mix_cl_star)
        deallocate( mix_cf_star)
        deallocate( mix_cf2_star)
        deallocate( mix_rain_star)
        deallocate( mix_graup_star)

      else  ! L_mix_ratio=.false.
! DEPENDS ON: diff_divdamp_ctl
          Call Diff_Divdamp_Ctl(                                        &
     &                     L_diffusion, L_cdiffusion,                   &
     &                     L_vertical_diffusion, L_divdamp,             &
     &                     L_ramp, ramp_lat_radians,                    &
     &                     L_Backwards,                                 &
     &                     timestep, pos_timestep, neg_timestep,        &
     &                     theta, w, q,                                 &
     &                     u, v, rho, Exner, w_adv,                     &
     &                     r_theta_levels, r_rho_levels, r_at_u, r_at_v,&
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     sec_theta_latitude,                          &
     &                     cos_theta_latitude, sec_v_latitude,          &
     &                     FV_sec_theta_latitude, cos_v_latitude,       &
     &                     sin_theta_latitude, sin_v_latitude,          &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     at_extremity, gc_proc_row_group,             &
     &                     mype, nproc, nproc_x, nproc_y, neighbour,    &
     &                     delta_lambda, delta_phi, L_regular,          &
     &                     lambda_p, phi_p, lambda_u, phi_v,            &
     &                     recip_dlamp, recip_dphip,                    &
     &                     recip_dlamu, recip_dphiv,                    &
     &                     rows, row_length, n_rows,                    &
     &                     model_levels, bl_levels, model_domain,      &
     &                     global_row_length, global_rows,              &
     &                     diffusion_coefficient_thermo,                &
     &                     diffusion_coefficient_w,                     &
     &                     diffusion_coefficient_q,                     &
     &                     diffusion_coefficient_wind,                  &
     &                     diffusion_order_thermo,                      &
     &                     diffusion_order_w,                           &
     &                     diffusion_order_q,                           &
     &                     diffusion_order_wind,                        &
     &                     horizontal_level, tar_horizontal,            &
     &                     level_start_wind, level_stop_wind,           &
     &                     level_start_q, level_stop_q,                 &
     &                     level_start_theta, level_stop_theta,         &
     &                     L_tardiff_q, w_conv_limit, tardiffq_factor,  &
     &                     tardiffq_test, tardiffq_start, tardiffq_end, &
     &                     L_diag_w, w_local_mask,                      &
     &                     L_adjust_theta,                              &
     &                     adjust_theta_start, adjust_theta_end,        &
     &                     L_vdiff_uv, vdiffuv_test,                    &
     &                     vdiffuv_factor, vdiffuv_start, vdiffuv_end,  &
     &                     vert_diffusion_coeff_wind,                   &
     &                     vert_diffusion_coeff_q,                      &
     &                     vert_diffusion_coeff_theta,                  &
     &                     div_damp_coefficient,                        &
     &                     theta_star, R_w, q_star, R_u, R_v            &
     &                     )
      endif  !L_mix_ratio

      IF (lhook) CALL dr_hook('NI_DIFF_CTL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE NI_DIFF_CTL
