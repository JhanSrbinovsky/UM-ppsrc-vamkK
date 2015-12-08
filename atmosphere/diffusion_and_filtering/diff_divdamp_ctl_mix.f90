! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Diff_Divdamp_Ctl_mix
          Subroutine Diff_Divdamp_Ctl_mix(                              &
     &                     L_diffusion, L_cdiffusion,L_subfilter_horiz, &
     &                     L_vertical_diffusion, L_divdamp,             &
     &                     L_ramp, ramp_lat_radians,                    &
     &                     L_Backwards, Ltimer,                         &
     &                     timestep, pos_timestep, neg_timestep,        &
     &                     theta, w, mix_v,                             &
     &                     u, v, rho, w_adv,                            &
     &                     r_theta_levels, r_rho_levels, r_at_u, r_at_v,&
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     sec_theta_latitude,                          &
     &                     cos_theta_latitude, sec_v_latitude,          &
     &                     FV_sec_theta_latitude, cos_v_latitude,       &
     &                     sin_theta_latitude, sin_v_latitude, pi,      &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     at_extremity, gc_proc_row_group,             &
     &                     mype, nproc, nproc_x, nproc_y, neighbour,    &
     &                     delta_lambda, delta_phi, L_regular,          &
     &                     lambda_p, phi_p, lambda_u, phi_v,            &
     &                     recip_dlamp, recip_dphi_p,                   &                     &
     &                     rows, row_length, n_rows,                    &
     &                     model_levels, wet_levels, model_domain,      &
     &                     global_row_length, global_rows,              &
     &                     diffusion_coefficient_thermo,                &
     &                     diffusion_coefficient_w,                     &
     &                     diffusion_coefficient_q,                     &
     &                     diffusion_coefficient_wind,                  &
     &                     diffusion_order_thermo,                      &
     &                     diffusion_order_w,                           &
     &                     diffusion_order_q,                           &
     &                     diffusion_order_wind,                        &
     &                     visc_m, visc_h,                              &     
     &                     horizontal_level, tar_horizontal,            &
     &                     level_start_wind, level_stop_wind,           &
     &                     level_start_q, level_stop_q,                 &
     &                     level_start_theta, level_stop_theta,         &
     &                     vert_diffusion_coeff_wind,                   &
     &                     vert_diffusion_coeff_q,                      &
     &                     vert_diffusion_coeff_theta,                  &
     &                     div_damp_coefficient,                        &
     &                     L_tardiff_q, w_conv_limit, tardiffq_factor,  &
     &                     tardiffq_test, tardiffq_start, tardiffq_end, &
     &                     L_diag_w, w_local_mask,                      &
     &                     theta_star, R_w, mix_v_star, R_u, R_v        &
     &                     )

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
      Implicit None

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, wet_levels                                                      &
                         ! number of wet levels.
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
     &,  global_rows                                                    &
     &, gc_proc_row_group                                               &
     &, nproc                                                           &
                  ! Total number of processors
     &, nproc_x                                                         &
                   ! Number of processors in longitude
     &, nproc_y                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)         ! Array with the Ids of the four neighbours
                             ! in the horizontal plane

      REAL                                                              &
     & timestep                                                         &
                            ! atmosphere model timestep
     &,pos_timestep                                                     &
                    ! = +timestep.
     &,neg_timestep                                                     &
                    ! = -timestep.
     &,  Pi                                                             &
     &,  w_conv_limit

      Real                                                              &
           ! primary model variables
     &  u(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
     &, v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)     &
     &, w(1-offx:row_length+offx, 1-offy:rows+offy, 0:model_levels)     &
     &, rho(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)     &
     &, mix_v (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &     wet_levels)                                                  &
     &, theta (1-offx:row_length+offx, 1-offy:rows+offy,                &
     &         model_levels)                                            &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &          0:model_levels)

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
     &  diffusion_coefficient(model_levels)                             &
     &, tardiffq_factor       ! targeted diffusion coefficient

      Logical                                                           &
     &  L_diffusion                                                     &
     &, L_cdiffusion                                                    &
     &, L_tardiff_q                                                     &
     &, L_ramp                                                          &
     &, L_vertical_diffusion                                            &
     &, L_divdamp                                                       &
     &, L_Backwards                                                     &
     &, Ltimer                                                          &
     &, L_diag_w                                                        &
     &, L_subfilter_horiz    ! subgrid turbulence scheme

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
            ! VarRes horizontal co-ordinate information
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
     &,diffusion_order_q(wet_levels)          !

      REAL                                                              &
     & diffusion_coefficient_thermo(model_levels)                       &
     &,diffusion_coefficient_wind(model_levels)                         &
     &,diffusion_coefficient_w(model_levels-1)                          &
     &,diffusion_coefficient_q(wet_levels)
     &,visc_m(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j          &
     &                                             , model_levels)      &
!            ! diffusion coefficient for momentum from subgrid
!            ! turbulence scheme
     &,visc_h(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j          &
     &                                             , model_levels)      &
!            ! diffusion coefficient for heat and moisture from
!            ! subgrid turbulence scheme
! Interpolation weights
     &,w1(1-offx:row_length+offx)                                       &
     &,w2(1-offx:row_length+offx)                                       &
     &,w3(1-offx:row_length,1-offy:rows+offy)                           &
     &,w4(1-offx:row_length,1-offy:rows+offy)                           &
     &,weight1, weight2, weight3 ! vertical weights

      REAL, ALLOCATABLE :: diff_coeff_u     (:,:,:)
      REAL, ALLOCATABLE :: diff_coeff_v     (:,:,:)
      REAL, ALLOCATABLE :: diff_coeff_th    (:,:,:)
      REAL, ALLOCATABLE :: diff_coeff_centre(:,:,:)
                             ! visc_m interpolated onto appropriate
                             ! points for use in the turb_diff_*
                             ! diffusion routines.
      REAL, ALLOCATABLE :: w_diff_coeff_u   (:,:,:)
      REAL, ALLOCATABLE :: w_diff_coeff_v   (:,:,:)
                             ! horizontal diffusion coefficients for
                             ! diffusion of w

! local
       Integer                                                          &
     & i,j,k

      Integer                                                           &
     & horizontal_level                                                 &
     &, tar_horizontal

      Integer                                                           &
     & level_start_wind, level_stop_wind                                &
     &,level_start_theta, level_stop_theta                              &
     &,level_start_q, level_stop_q                                      &
     &, tardiffq_test                                                   &
     &, tardiffq_start                                                  &
     &, tardiffq_end

      Real                                                              &
     & vert_diffusion_coeff_wind                                        &
     &,vert_diffusion_coeff_theta                                       &
     &,vert_diffusion_coeff_q                                           &
     &, ramp_lat_radians

      Real                                                              &
     & div_damp_coefficient(model_levels)

!   variables with intent INOUT

      Real                                                              &
     &  R_u(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &        model_levels)                                             &
     &, R_v(1-offx:row_length+offx, 1-offy:n_rows+offy,                 &
     &        model_levels)                                             &
     &, R_w(row_length, rows, model_levels)                             &
     &, mix_v_star(1-offx:row_length+offx,                              &
     &           1-offy:rows+offy, wet_levels)                          &
     &, theta_star(1-offx:row_length+offx,                              &
     &               1-offy:rows+offy, model_levels)                    &
     &, w_local_mask(row_length,rows)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Function & Subroutine calls:
      External                                                          &
     &  h_diff_theta                                                    &
     &, h_diff_u                                                        &
     &, h_diff_v                                                        &
     &, h_diff_q                                                        &
     &, vert_diff_u                                                     &
     &, vert_diff_v                                                     &
     &, vert_diff_theta                                                 &
     &, vert_diff_q                                                     &
     &, div_damp


! ----------------------------------------------------------------------
! Section 1.0  Horizontal Diffusion
! ----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('DIFF_DIVDAMP_CTL_MIX',zhook_in,zhook_handle)
      If ( L_subfilter_horiz) Then
        ALLOCATE (diff_coeff_u(1-offx:row_length, rows                  &
     &                                        , model_levels))
        ALLOCATE (diff_coeff_v(row_length, 1-offy:n_rows+offy           &
     &                                        , model_levels))
        ALLOCATE (diff_coeff_th(1-offx:row_length+offx, 1-offy:rows+offy&
     &                                        , model_levels))
        ALLOCATE (diff_coeff_centre(1-offx:row_length, 1-offy:rows,     &
     &                                          model_levels))
        ALLOCATE (w_diff_coeff_u(1-offx:row_length, rows                &
     &                                        , model_levels))
        ALLOCATE (w_diff_coeff_v(row_length, 1-offy:n_rows+offy         &
     &                                        , model_levels))
      Else
        ALLOCATE (diff_coeff_u(1,1,1))
        ALLOCATE (diff_coeff_v(1,1,1))
        ALLOCATE (diff_coeff_th(1,1,1))
        ALLOCATE (diff_coeff_centre(1,1,1))
        ALLOCATE (w_diff_coeff_u(1,1,1))
        ALLOCATE (w_diff_coeff_v(1,1,1))
      End If

      If ( L_diffusion ) Then

!       Run with a positive timestep if integrating backwards.
        IF (L_Backwards) timestep = pos_timestep

! DEPENDS ON: h_diff_theta
        Call h_diff_theta(                                              &
     &                     theta,                                       &
     &                     r_theta_levels,                              &
     &                     offx, offy, halo_i, halo_j, offx, offy,      &
     &                     at_extremity, gc_proc_row_group,             &
     &                     nproc, nproc_x, nproc_y, neighbour,          &
     &                     delta_lambda, delta_phi,                     &
     &                     timestep, rows, row_length,                  &
     &                     model_levels, model_domain,                  &
     &                     global_row_length,                           &
     &                     diffusion_coefficient_thermo,                &
     &                     diffusion_order_thermo, theta_star,          &
     &                     horizontal_level)

! diffusion of w not called, code left in case needs to be reinstated
!        Call h_diff_theta(
!     &                     w(1-offx, 1-offy, 1),
!     &                     r_theta_levels,
!     &                     offx, offy, halo_i, halo_j, 0, 0,
!     &                     at_extremity, gc_proc_row_group,
!     &                     nproc, nproc_x, nproc_y, neighbour,
!     &                     delta_lambda, delta_phi,
!     &                     timestep, rows, row_length,
!     &                     model_levels-1, model_domain,
!     &                     global_row_length,
!     &                     diffusion_coefficient_w,
!     &                     diffusion_order_w, R_w,
!     &                     horizontal_level)

! DEPENDS ON: h_diff_q
        Call h_diff_q(                                                  &
     &                     mix_v,                                       &
     &                     r_theta_levels,                              &
     &                     offx, offy, halo_i, halo_j, offx, offy,      &
     &                     at_extremity, gc_proc_row_group,             &
     &                     nproc, nproc_x, nproc_y, neighbour,          &
     &                     delta_lambda, delta_phi,                     &
     &                     timestep, rows, row_length,                  &
     &                     model_levels,                                &
     &                     wet_levels, model_domain,                    &
     &                     global_row_length,                           &
     &                     diffusion_coefficient_q,                     &
     &                     diffusion_order_q, mix_v_star,               &
     &                     horizontal_level)

! DEPENDS ON: h_diff_u
        Call h_diff_u(                                                  &
     &                 u, r_at_u, r_theta_levels,                       &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 at_extremity, gc_proc_row_group,                 &
     &                 nproc, nproc_x, nproc_y, neighbour,              &
     &                 delta_lambda, delta_phi,                         &
     &                 timestep, rows, row_length,                      &
     &                 model_levels, model_domain,                      &
     &                 diffusion_coefficient_wind,                      &
     &                 diffusion_order_wind, R_u,                       &
     &                 horizontal_level)

! DEPENDS ON: h_diff_v
        Call h_diff_v(                                                  &
     &                 v, r_at_v, r_theta_levels, rows,                 &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 at_extremity, gc_proc_row_group,                 &
     &                 nproc, nproc_x, nproc_y, neighbour,              &
     &                 delta_lambda, delta_phi,                         &
     &                 timestep, n_rows, row_length,                    &
     &                 model_levels, model_domain,                      &
     &                 diffusion_coefficient_wind,                      &
     &                 diffusion_order_wind, R_v,                       &
     &                 horizontal_level)

!       Go back to negative timestep if integrating backwards.
        IF (L_Backwards) timestep = neg_timestep

      End If        !  L_diffusion
! If L_cdiffusion .true. then call conservative diffusion
      If ( L_cdiffusion ) Then

!       Run with a positive timestep if integrating backwards.
        IF (L_Backwards) timestep = pos_timestep

! DEPENDS ON: h_cdiff_theta
        Call h_cdiff_theta(                                             &
     &                     theta,                                       &
     &                     r_theta_levels, r_rho_levels,                &
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     sec_theta_latitude, cos_v_latitude,          &
     &                     offx, offy, halo_i, halo_j, offx, offy,      &
     &                     at_extremity, gc_proc_row_group,             &
     &                     nproc, nproc_x, nproc_y, neighbour,          &
     &                     delta_lambda, delta_phi,                     &
     &                     timestep, rows, n_rows, row_length,          &
     &                     model_levels, model_domain,                  &
     &                     global_row_length,                           &
     &                     diffusion_coefficient_thermo,                &
     &                     diffusion_order_thermo, theta_star,          &
     &                     horizontal_level)

! DEPENDS ON: h_cdiff_q
        Call h_cdiff_q(                                                 &
     &                     mix_v,                                       &
     &                     r_theta_levels, r_rho_levels,                &
     &                     eta_theta_levels, eta_rho_levels,            &
     &                     sec_theta_latitude, cos_v_latitude,          &
     &                     offx, offy, halo_i, halo_j, offx, offy,      &
     &                     at_extremity, gc_proc_row_group,             &
     &                     nproc, nproc_x, nproc_y, neighbour,          &
     &                     delta_lambda, delta_phi,                     &
     &                     timestep, rows, n_rows, row_length,          &
     &                     model_levels,                                &
     &                     wet_levels, model_domain,                    &
     &                     global_row_length,                           &
     &                     diffusion_coefficient_q,                     &
     &                     diffusion_order_q, mix_v_star,               &
     &                     horizontal_level)

! DEPENDS ON: h_cdiff_u
        Call h_cdiff_u(                                                 &
     &                 u, r_at_u, r_theta_levels,                       &
     &                 eta_theta_levels, eta_rho_levels,                &
     &                 sec_theta_latitude, cos_v_latitude,              &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 at_extremity, gc_proc_row_group,                 &
     &                 nproc, nproc_x, nproc_y, neighbour,              &
     &                 delta_lambda, delta_phi,                         &
     &                 timestep, rows, n_rows, row_length,              &
     &                 model_levels, model_domain,                      &
     &                 diffusion_coefficient_wind,                      &
     &                 diffusion_order_wind, R_u,                       &
     &                 horizontal_level)

! DEPENDS ON: h_cdiff_v
        Call h_cdiff_v(                                                 &
     &                 v, r_at_v, r_theta_levels, rows,                 &
     &                 eta_theta_levels, eta_rho_levels,                &
     &                 cos_theta_latitude, sec_v_latitude,              &
     &                 offx, offy, halo_i, halo_j,                      &
     &                 at_extremity, gc_proc_row_group,                 &
     &                 nproc, nproc_x, nproc_y, neighbour,              &
     &                 delta_lambda, delta_phi,                         &
     &                 timestep, n_rows, row_length,                    &
     &                 model_levels, model_domain,                      &
     &                 diffusion_coefficient_wind,                      &
     &                 diffusion_order_wind, R_v,                       &
     &                 horizontal_level)

!       Go back to negative timestep if integrating backwards.
        IF (L_Backwards) timestep = neg_timestep

      End If        !  L_cdiffusion

! Test to apply targetted diffusion on q if w exceeds prescribed limit
        IF( L_tardiff_q) then
          if( tar_horizontal > 0 ) then
! DEPENDS ON: tardiff_q_wss
            Call tardiff_q_wss(                                         &
     &                         mix_v, w_adv,                            &
     &                         r_theta_levels, r_rho_levels,            &
     &                         sec_theta_latitude,                      &
     &                         cos_theta_latitude, cos_v_latitude,      &
     &                         offx, offy, halo_i, halo_j, offx, offy,  &
     &                         at_extremity, gc_proc_row_group,         &
     &                         delta_lambda, delta_phi,                 &
     &                         timestep, rows, n_rows, row_length,      &
     &                         model_levels, wet_levels,                &
     &                         model_domain, global_row_length,         &
     &                         mix_v_star, w_conv_limit, tar_horizontal,&
     &                         tardiffq_factor, tardiffq_test,          &
     &                         tardiffq_start, tardiffq_end,            &
     &                         L_diag_w, w_local_mask )
          else  !  tar_horizontal = 0
! DEPENDS ON: tardiff_q_w
          Call tardiff_q_w(                                             &
     &                    mix_v, w_adv,                                 &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    sec_theta_latitude,                           &
     &                    cos_theta_latitude, cos_v_latitude,           &
     &                    offx, offy, halo_i, halo_j, offx, offy,       &
     &                    at_extremity, gc_proc_row_group,              &
     &                    delta_lambda, delta_phi,                      &
     &                    timestep, rows, n_rows, row_length,           &
     &                    model_levels, wet_levels,                     &
     &                    model_domain, global_row_length,              &
     &                    mix_v_star, w_conv_limit,                     &
     &                    tardiffq_factor, tardiffq_test,               &
     &                    tardiffq_start, tardiffq_end,                 &
     &                    L_diag_w, w_local_mask)
          endif !  tar_horizontal > 0
        ENDIF !   L_tardiff_q

       If ( L_subfilter_horiz) Then
!        Run with a positive timestep if integrating backwards.
         If (L_Backwards) timestep = pos_timestep
!
!Interpolate visc_m onto appropriate points for the diffusion routines.
! At this point the diffusion coefficients are on w (theta) points.
!
! Horizontal weights
!
      If (L_regular) then  ! regular grid

         Do k = 2, model_levels

           Do j = 1, rows
             Do i = 1-offx,row_length
               diff_coeff_u(i,j,k) = 0.5*(                              &
     &                   visc_h(i+1,j,k) + visc_h(i,j,k) )
               w_diff_coeff_u(i,j,k) = 0.5*(                            &
     &                   visc_m(i+1,j,k) + visc_m(i,j,k) )
             End Do
           End Do

           Do j = 1-offy, n_rows+offy
             Do i = 1, row_length
               diff_coeff_v(i,j,k) = 0.5*(                              &
     &                   visc_h(i,j+1,k) + visc_h(i,j,k) )
               w_diff_coeff_v(i,j,k) = 0.5*(                            &
     &                   visc_m(i,j+1,k) + visc_m(i,j,k) )
             End Do
           End Do

           Do j = 1-offy, rows + offy
             Do i = 1-offx, row_length+offx
               weight1 =                                                &
     &           r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
               weight2 =                                                &
     &           r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
               weight3 =                                                &
     &           r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
               diff_coeff_th(i,j,k) = ( weight3*visc_m(i,j,k) +         &
     &               weight2*visc_m(i,j,k-1) )/weight1
             End Do
           End Do

           Do j = 1-offy, rows
             Do i = 1-offx, row_length
               diff_coeff_centre(i,j,k) = 0.25*(                         &
     &             diff_coeff_th(i,j,k) + diff_coeff_th(i+1,j,k) +       &
     &             diff_coeff_th(i,j+1,k) + diff_coeff_th(i+1,j+1,k)  )
             End Do
           End Do

         End Do

         Do j = 1, rows
           Do i = 1-offx, row_length
             diff_coeff_u(i,j,1) = 0.5*(                                &
     &                 visc_h(i+1,j,1) + visc_h(i,j,1) )
             w_diff_coeff_u(i,j,1) = 0.5*(                              &
     &                 visc_m(i+1,j,1) + visc_m(i,j,1) )
           End Do
         End Do

         Do j = 1-offy, n_rows+offy
           Do i = 1, row_length
             diff_coeff_v(i,j,1) = 0.5*(                                &
     &                   visc_h(i,j+1,1) + visc_h(i,j,1) )
             w_diff_coeff_v(i,j,1) = 0.5*(                              &
     &                   visc_m(i,j+1,1) + visc_m(i,j,1) )
           End Do
         End Do

         Do j = 1-offy, rows + offy
           Do i = 1-offx, row_length+offx
             diff_coeff_th(i,j,1) = visc_m(i,j,1)
           End Do
         End Do

         Do j = 1-offy, rows
           Do i = 1-offx, row_length
             diff_coeff_centre(i,j,1) = 0.25* (                         &
     &         visc_m(i,j,1) + visc_m(i,j+1,1) +                        &
     &         visc_m(i+1,j,1) + visc_m(i+1,j+1,1)  )
           End Do
         End Do

      Else ! variable resolution      

         Do i = 1-offx, row_length+offx
           w1(i) = (lambda_p(i+1) - lambda_u(i))*recip_dlamp(i)
           w2(i) = (lambda_u(i) - lambda_p(i))*recip_dlamp(i)
         End Do

         Do j = 1-offy, n_rows+offy
           Do i = 1-offx, row_length
             w3(i,j) = (phi_v(i,j) - phi_p(i,j))*recip_dphip(i,j)
             w4(i,j) = (phi_p(i,j+1)- phi_v(i,j))*recip_dphip(i,j)           
           End Do
         End Do

         Do i = 1, row_length
           w3(i,rows + offy) = w3(i,n_rows + offy)
           w4(i,rows + offy) = w4(i,n_rows + offy)
         End Do

         Do k = 2, model_levels

           Do j = 1, rows
             
             Do i = 1-offx,row_length
               diff_coeff_u(i,j,k) = w1(i)*visc_h(i,j,k)                &
     &                             + w2(i)*visc_h(i+1,j,k)
               w_diff_coeff_u(i,j,k) = w1(i)*visc_m(i,j,k)              &
     &                               + w2(i)*visc_m(i+1,j,k)     
             End Do
           End Do

           Do j = 1-offy, n_rows+offy
             Do i = 1, row_length
               diff_coeff_v(i,j,k) = w3(i,j)*visc_h(i,j,k)              &
     &                             + w4(i,j)*visc_h(i,j+1,k)
               w_diff_coeff_v(i,j,k) = w3(i,j)*visc_m(i,j,k)            &
     &                               + w4(i,j)*visc_m(i,j+1,k)  
             End Do
           End Do

           Do j = 1-offy, rows + offy
             Do i = 1-offx, row_length+offx
               weight1 =                                                &
     &           r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
               weight2 =                                                &
     &           r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
               weight3 =                                                &
     &           r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1)
               diff_coeff_th(i,j,k) = ( weight3*visc_m(i,j,k) +         &
     &               weight2*visc_m(i,j,k-1) )/weight1
             End Do
           End Do

           Do j = 1-offy, rows
             Do i = 1-offx, row_length
               diff_coeff_centre(i,j,k) =                               &
     &            w4(i,j)*( w1(i)*diff_coeff_th(i,j,k)                  &
     &                  + w2(i)*diff_coeff_th(i+1,j,k)  )               &
     &          + w3(i,j)*( w1(i)*diff_coeff_th(i,j+1,k)                &
                        + w2(i)*diff_coeff_th(i+1,j+1,k)  )
             End Do
           End Do

         End Do

         Do j = 1, rows
           Do i = 1-offx, row_length
             diff_coeff_u(i,j,1) = w1(i)*visc_h(i,j,1)                  &
     &                           + w2(i)*visc_h(i+1,j,1)
             w_diff_coeff_u(i,j,1) = w1(i)*visc_m(i,j,1)                &
     &                             + w2(i)*visc_m(i+1,j,1)
           End Do
         End Do

         Do j = 1-offy, n_rows+offy
           Do i = 1, row_length
             diff_coeff_v(i,j,1) = w3(i,j)*visc_h(i,j,1)                &
     &                           + w4(i,j)*visc_h(i,j+1,1)
             w_diff_coeff_v(i,j,1) = w3(i,j)*visc_m(i,j,1)              &
     &                             + w4(i,j)*visc_m(i,j+1,1) 
           End Do
         End Do

         Do j = 1-offy, rows + offy
           Do i = 1-offx, row_length+offx
             diff_coeff_th(i,j,1) = visc_m(i,j,1)
           End Do
         End Do

         Do j = 1-offy, rows
           Do i = 1-offx, row_length
             diff_coeff_centre(i,j,1) =                                 &
     &         w4(i,j)*( w1(i)*visc_m(i,j,1) + w2(i)*visc_m(i+1,j,1) )  &
     &        +w3(i,j)*( w1(i)*visc_m(i,j+1,1) +                        & 
     &                   w2(i)* visc_m(i+1,j+1,1))
           Enddo
         Enddo

      End If 
!
!  Call diffusion routines with (3D) coefficients from the subgrid
!  turbulence scheme.
!
         If (L_regular) then  ! regular grid
! DEPENDS ON: turb_diff_th
           Call turb_diff_th(                                           &
     &                    theta,                                        &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    sec_theta_latitude,                           &
     &                    cos_v_latitude, sec_v_latitude,               &
     &                    offx, offy, halo_i, halo_j,                   &
     &                    delta_lambda, delta_phi, timestep,            &
     &                    rows, n_rows, row_length, model_levels,       &
     &                    model_levels, diff_coeff_u, diff_coeff_v,     &
     &                    theta_star)
         Else ! variable resolution
! DEPENDS ON: turb_diff_th_VarRes
           Call turb_diff_th_VarRes(                                    &
     &                        theta,                                    &
     &                        r_theta_levels, r_rho_levels,             &
     &                        sec_theta_latitude,                       &
     &                        cos_v_latitude, sec_v_latitude,           &
     &                        offx, offy, halo_i, halo_j,               &
     &                        recip_dlamp, recip_dphip,                 &
     &                        recip_dlamu, recip_dphiv,                 &
     &                        timestep,                                 &
     &                        rows, n_rows, row_length, model_levels,   &
     &                        model_levels, diff_coeff_u, diff_coeff_v, &
     &                        theta_star)
         End If

         If (L_regular) then  ! regular grid
! DEPENDS ON: turb_diff_q
           Call turb_diff_q(                                            &
     &                    moist,                                        &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    sec_theta_latitude,                           &
     &                    cos_v_latitude, sec_v_latitude,               &
     &                    offx, offy, halo_i, halo_j,                   &
     &                    delta_lambda, delta_phi, timestep,            &
     &                    rows, n_rows, row_length,                     &
     &                    model_levels, wet_levels, wet_levels,         &
     &                    diff_coeff_u, diff_coeff_v,                   &
     &                    moist_star)
         Else ! variable resolution
! DEPENDS ON: turb_diff_q_VarRes
           Call turb_diff_q_VarRes(                                     &
     &                       moist,                                     &
     &                       r_theta_levels, r_rho_levels,              &
     &                       sec_theta_latitude,                        &
     &                       cos_v_latitude, sec_v_latitude,            &
     &                       offx, offy, halo_i, halo_j,                &
     &                       recip_dlamp, recip_dphip,                  &
     &                       recip_dlamu, recip_dphiv,                  &
     &                       timestep,                                  &
     &                       rows, n_rows, row_length,                  &
     &                       model_levels, wet_levels, wet_levels,      &
     &                       diff_coeff_u, diff_coeff_v,                &
     &                       moist_star)
         End If

         If (L_regular) then  ! regular grid
! DEPENDS ON: turb_diff_u
           Call turb_diff_u(                                            &
     &                    u, r_at_u,                                    &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    sec_theta_latitude,                           &
     &                    cos_v_latitude, sec_v_latitude,               &
     &                    offx, offy, halo_i, halo_j,                   &
     &                    delta_lambda, delta_phi, timestep,            &
     &                    rows, n_rows, row_length,                     &
     &                    model_levels, model_levels,                   &
     &                    diff_coeff_th, diff_coeff_centre,             &
     &                    R_u )
         Else ! variable resolution
! DEPENDS ON: turb_diff_u_VarRes
           Call turb_diff_u_VarRes(                                     &
     &                    u, r_at_u,                                    &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    sec_theta_latitude,                           &
     &                    cos_v_latitude, sec_v_latitude,               &
     &                    offx, offy, halo_i, halo_j,                   &
     &                    recip_dlamp, recip_dphip,                     &
     &                    recip_dlamu, recip_dphiv,                     &
     &                    timestep,                                     &
     &                    rows, n_rows, row_length,                     &
     &                    model_levels, model_levels,                   &
     &                    diff_coeff_th, diff_coeff_centre, R_u )
         End If

         If (L_regular) then  ! regular grid
! DEPENDS ON: turb_diff_v
           Call turb_diff_v(                                            &
     &                    v, r_at_v,                                    &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    sec_v_latitude,                               &
     &                    cos_theta_latitude, sec_theta_latitude,       &
     &                    offx, offy, halo_i, halo_j,                   &
     &                    delta_lambda, delta_phi, timestep,            &
     &                    rows, n_rows, row_length,                     &
     &                    model_levels, model_levels,                   &
     &                    diff_coeff_th, diff_coeff_centre,             &
     &                    R_v )
         Else ! variable resolution
! DEPENDS ON: turb_diff_v_VarRes
           Call turb_diff_v_VarRes(                                     &
     &                       v, r_at_v,                                 &
     &                       r_theta_levels, r_rho_levels,              &
     &                       sec_v_latitude,                            &
     &                       cos_theta_latitude, sec_theta_latitude,    &
     &                       offx, offy, halo_i, halo_j,                &
     &                       recip_dlamp, recip_dphip,                  &
     &                       recip_dlamu, recip_dphiv,                  &
     &                       timestep,                                  &
     &                       rows, n_rows, row_length,                  &
     &                       model_levels, model_levels,                &
     &                       diff_coeff_th, diff_coeff_centre, R_v )
         End If

         If (L_regular) then  ! regular grid
! DEPENDS ON: turb_diff_w
           Call turb_diff_w(                                            &
     &                   w,                                             &
     &                   r_theta_levels, r_rho_levels,                  &
     &                   sec_theta_latitude, cos_v_latitude,            &
     &                   sec_v_latitude,                                &
     &                   offx, offy, halo_i, halo_j,                    &
     &                   delta_lambda, delta_phi,                       &
     &                   timestep, rows, n_rows, row_length,            &
     &                   model_levels, model_levels - 1,                &
     &                   w_diff_coeff_u, w_diff_coeff_v,                &
     &                   R_w )
         Else ! variable resolution
! DEPENDS ON: turb_diff_w_VarRes
           Call turb_diff_w_VarRes(                                     &
     &                       w,                                         &
     &                       r_theta_levels, r_rho_levels,              &
     &                       sec_theta_latitude,                        &
     &                       cos_v_latitude, sec_v_latitude,            &
     &                       offx, offy, halo_i, halo_j,                &
     &                       recip_dlamp, recip_dphip,                  &
     &                       recip_dlamu, recip_dphiv,                  &
     &                       timestep, rows, n_rows, row_length,        &
     &                       model_levels, model_levels-1,              &
     &                       w_diff_coeff_u, w_diff_coeff_v,            & 
     &                       R_w)
         End If

!        Go back to negative timestep if integrating backwards.
         IF (L_Backwards) timestep = neg_timestep
       End If   !L_subfilter_horiz
+

! ----------------------------------------------------------------------
! Section 2.0  Vertical Diffusion
! ----------------------------------------------------------------------

      IF (L_vertical_diffusion) then

!       Run with a positive timestep if integrating backwards.
        IF (L_Backwards) timestep = pos_timestep

        IF (vert_diffusion_coeff_wind  >   tiny(1.0)) then
! DEPENDS ON: vert_diff_u
          Call vert_diff_u(                                             &
     &                     u, r_theta_levels, r_at_u,                   &
     &                     sin_theta_latitude, Pi, L_ramp,              &
     &                     ramp_lat_radians,                            &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     at_extremity,                                &
     &                     timestep, rows, row_length,                  &
     &                     model_levels, model_domain,                  &
     &                     level_start_wind, level_stop_wind,           &
     &                     vert_diffusion_coeff_wind,                   &
     &                     R_u )

! DEPENDS ON: vert_diff_v
          Call vert_diff_v(                                             &
     &                     v, r_theta_levels, r_at_v,                   &
     &                     sin_theta_latitude, Pi, L_ramp,              &
     &                     ramp_lat_radians,                            &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     at_extremity,                                &
     &                     timestep, rows, n_rows, row_length,          &
     &                     model_levels, model_domain,                  &
     &                     level_start_wind, level_stop_wind,           &
     &                     vert_diffusion_coeff_wind,                   &
     &                     R_v )
        Endif  ! vert_diffusion_coeff_wind

        IF (vert_diffusion_coeff_q  >   tiny(1.0)) then
! DEPENDS ON: vert_diff_q
          Call vert_diff_q(                                             &
     &                     mix_v, r_theta_levels, r_rho_levels,         &
     &                     sin_theta_latitude, Pi, L_ramp,              &
     &                     ramp_lat_radians,                            &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     offx, offy,                                  &
     &                     at_extremity,                                &
     &                     timestep, rows, row_length,                  &
     &                     model_levels, wet_levels, model_domain,      &
     &                     level_start_q, level_stop_q,                 &
     &                     vert_diffusion_coeff_q,                      &
     &                     mix_v_star )
        Endif  ! vert_diffusion_coeff_q

        IF (vert_diffusion_coeff_theta  >   tiny(1.0)) then
! DEPENDS ON: vert_diff_theta
          Call vert_diff_theta(                                         &
     &                     theta, r_theta_levels, r_rho_levels,         &
     &                     sin_theta_latitude, Pi, L_ramp,              &
     &                     ramp_lat_radians,                            &
     &                     offx, offy, halo_i, halo_j,                  &
     &                     at_extremity,                                &
     &                     timestep, rows, row_length,                  &
     &                     model_levels, model_domain,                  &
     &                     level_start_theta, level_stop_theta,         &
     &                     vert_diffusion_coeff_theta,                  &
     &                     theta_star)
        Endif  ! vert_diffusion_coeff_q

!       Go back to negative timestep if integrating backwards.
        IF (L_Backwards) timestep = neg_timestep

      END IF     !  (L_vertical_diffusion)

! ----------------------------------------------------------------------
! Section 3.0  Divergence damping
! ----------------------------------------------------------------------

      If (L_divdamp) Then

!       Run with a positive timestep if integrating backwards.
        IF (L_Backwards) timestep = pos_timestep

! Call divergence damping code
! DEPENDS ON: div_damp
        Call div_damp(                                                  &
     &                u, v, rho,                                        &
     &                r_at_u, r_at_v, r_rho_levels,                     &
     &                FV_sec_theta_latitude, cos_v_latitude,            &
     &                offx, offy, halo_i, halo_j,                       &
     &                mype, at_extremity, gc_proc_row_group,            &
     &                nproc, nproc_x, nproc_y, neighbour,               &
     &                delta_lambda, delta_phi,                          &
     &                timestep, rows, n_rows, row_length,               &
     &                model_levels, model_domain,                       &
     &                global_row_length, div_damp_coefficient,          &
     &                R_u, R_v)

!       Go back to negative timestep if integrating backwards.
        IF (L_Backwards) timestep = neg_timestep

      End If         !  L_divdamp

       If (L_subfilter_horiz) then
         DEALLOCATE (diff_coeff_u)
         DEALLOCATE (diff_coeff_v)
         DEALLOCATE (diff_coeff_th)
         DEALLOCATE (diff_coeff_centre)
         DEALLOCATE (w_diff_coeff_u)
         DEALLOCATE (w_diff_coeff_v)
       End If

!   call to w_div_diag removed from diffctl (now called from atmstep)

      IF (lhook) CALL dr_hook('DIFF_DIVDAMP_CTL_MIX',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Diff_Divdamp_Ctl_mix
