! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  subroutine update_rho
      subroutine update_rho(                                            &
     &                      timestep, NumCycles, CycleNo, rows, n_rows, &
     &                      row_length, model_levels, wet_model_levels, &
     &                      model_domain, first_constant_r_rho_level,   &
     &                      alpha_1, alpha_2, rims_to_do,               &
     &                      nproc, gc_proc_row_group,                   &
     &                      L_regular, at_extremity, global_row_length, &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      cos_v_latitude, delta_lambda, delta_phi,    &
     &                      dlambda_p, dphi_p,                          &
     &                      recip_dlambda_u, recip_dphi_v,              &
     &                      wt_lambda_p, wt_lambda_u,                   &
     &                      wt_phi_p, wt_phi_v,                         &
     &                      R_u, R_v, R_w,                              &
     &                      u, v, w, rho, rho_np1,                      &
     &                      vap, cl, cf, cf2, rain, graup,              &
     &                      vap_star, cl_star, cf_star,                 &
     &                      cf2_star, rain_star, graup_star,            &
     &                      vap_np1, cl_np1, cf_np1,                    &
     &                      cf2_np1, rain_np1, graup_np1,               &
     &                      L_mcr_cf2, L_mcr_rain, L_mcr_graup,         &
     &                      r_theta_levels, r_rho_levels,               &
     &                      eta_theta_levels, eta_rho_levels,           &
     &                      FV_sec_theta_latitude,                      &
     &                      wet_to_dry_n, wet_to_dry_np1,               &
     &                      rho_n, inc_rho, L_do_increment,             &
     &                      L_new_tdisc, L_mix_ratio, L_dry )

! Purpose:
!         subroutine to update rho, to tidy ATM_STEP
!
! Method:
!          Is described in ;
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Use swapable_field_mod, Only: &
          swapable_field_pointer_type

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
              ! model dimensions
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held.
     &, off_x                                                           &
     &, off_y                                                           &
     &, halo_i                                                          &
                     ! Size of halo in i.
     &, halo_j                                                          &
                     ! Size of halo in j.
     &, nproc                                                           &
                    ! Total number of processors
     &, gc_proc_row_group                                               &
                          ! Group id for processors on the same row
     &, global_row_length                                               &
     &, rims_to_do                                                      &
                        ! rim size of lbc weights = 1
     &, NumCycles                                                       &
     &, CycleNo

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_regular                                                       &
                    ! false if variable resolution
     &, L_mcr_cf2, L_mcr_rain, L_mcr_graup                              &
     &, L_mix_ratio                                                     &
     &, L_dry                                                           &
     &, L_new_tdisc

      Integer                                                           &
     &  first_constant_r_rho_level ! first rho level on which r
                                   ! is constant.

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
     &  timestep                                                        &
     &, alpha_1                                                         &
     &, alpha_2

      Real                                                              &
     &  vap (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &       wet_model_levels)                                          &
     &, cl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &       wet_model_levels)                                          &
     &, cf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &       wet_model_levels)                                          &
     &, cl_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &          wet_model_levels)                                       &
     &, cf_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &          wet_model_levels)                                       &
     &, vap_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &          wet_model_levels)                                       &
     &, cf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &       wet_model_levels)                                          &
     &, rain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &       wet_model_levels)                                          &
     &, graup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &       wet_model_levels)                                          &
     &, cf2_star(1-off_x:row_length+off_x,                              &
     &            1-off_y:rows+off_y, wet_model_levels)                 &
     &, rain_star(1-off_x:row_length+off_x,                             &
     &             1-off_y:rows+off_y, wet_model_levels)                &
     &, graup_star(1-off_x:row_length+off_x,                            &
     &             1-off_y:rows+off_y, wet_model_levels)                &
     &, vap_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &          wet_model_levels)                                       &
     &, cl_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          wet_model_levels)                                       &
     &, cf_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &          wet_model_levels)                                       &
     &, cf2_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &          wet_model_levels)                                       &
     &, rain_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &          wet_model_levels)                                       &
     &, graup_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &          wet_model_levels)                                       &
     &, rho_n (1-off_x:row_length+off_x,                                &
     &         1-off_y:rows+off_y, model_levels)

      Real, Intent(Out) ::                                              &
     &  inc_rho(1-off_x:row_length+off_x,                               &
     &             1-off_y:rows+off_y, model_levels)                    &
     &, wet_to_dry_n (1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
     &          model_levels)                                           &
     &, wet_to_dry_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,   &
     &          model_levels)


      Real                                                              &
     &  u(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)   &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y, model_levels) &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y, 0:model_levels) &
     &, R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels) &
     &, R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,             &
     &        model_levels)                                             &
     &, R_w(row_length, rows, model_levels)

      Real                                                              &
           ! trigonometric arrays.
     &  FV_sec_theta_latitude (1-off_x:row_length+off_x,                &
     &                         1-off_y:rows+off_y)                      &
     &, cos_v_latitude (1-off_x:row_length+off_x,                       &
     &                  1-off_y:n_rows+off_y)

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  eta_theta_levels (0:model_levels)                               &
     &, eta_rho_levels (model_levels)                                   &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
           !VarRes horizontal co-ordinate information
     &  dlambda_p(1-halo_i:row_length+halo_i)                           &
     &, dphi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)        &
     &, recip_dlambda_u(1-halo_i : row_length + halo_i)                 &
     &, recip_dphi_v( 1-halo_i : row_length + halo_i                    &
     &,               1-halo_j : n_rows+halo_j )                        &
     &, wt_lambda_p(1-halo_i:row_length+halo_i)                         &
     &, wt_lambda_u(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)      &
     &, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

      Logical :: L_do_increment

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  rho (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &       model_levels)                                              &
     &, rho_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &       model_levels)

! Local Variables.

      Logical                                                           &
     &  L_poles          ! true for etadot at poles

      Integer                                                           &
     &  i, j, k                                                         &
                      ! Loop indices
     &, i_start, i_stop                                                 &
     &, j_start, j_stop                                                 &
     &, j_begin, j_end                                                  &
     &, dim_u, dim_v
      Integer :: i_field

! Local arrays

      Real,Target ::                                                    &
     &  u1(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)  &
     &, v1(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y, model_levels)&
     &, w1(row_length, rows, 0:model_levels)                            &
     &, moist (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &          wet_model_levels)                                       &
     &, moist_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &              wet_model_levels)                                   &
     &, dry_to_wet (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &              model_levels)

      Real, Dimension(:,:,:), Allocatable::                             &
     &  moist_np1

      Type(swapable_field_pointer_type) :: fields_to_swap(3)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! ----------------------------------------------------------------------
! Section 0. set up moisture loading terms
!      Fields can be mixing ratios or specific quantities
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('UPDATE_RHO',zhook_in,zhook_handle)
      moist= vap + cl + cf
      moist_star= vap_star + cl_star + cf_star

      If(L_mcr_cf2)then
        moist      = moist + cf2
        moist_star = moist_star + cf2_star
      endif
      If(L_mcr_rain)then
        moist      = moist + rain
        moist_star = moist_star + rain_star
      endif
      If(L_mcr_graup)then
        moist      = moist + graup
        moist_star = moist_star + graup_star
      endif

      If ( CycleNo > 1 .and. L_new_tdisc ) Then
        Allocate( moist_np1 (1-off_x:row_length+off_x,                  &
     &            1-off_y:rows+off_y, wet_model_levels) )
        moist_np1 = vap_np1 + cl_np1 + cf_np1
        If ( L_mcr_cf2 ) moist_np1 = moist_np1 + cf2_np1
        If ( L_mcr_rain ) moist_np1 = moist_np1 + rain_np1
        If ( L_mcr_graup ) moist_np1 = moist_np1 + graup_np1
      Else
        Allocate( moist_np1 (1,1,1) )
      End If


        Do k = 1, model_levels
          Do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
              u1(i,j,k) = u(i,j,k) + alpha_1 * R_u(i,j,k)
            End Do
          End Do
        End Do
        Do k = 1, model_levels
          Do j = 1-off_y, n_rows+off_y
            Do i = 1-off_x, row_length+off_x
              v1(i,j,k) = v(i,j,k) + alpha_1 * R_v(i,j,k)
            End Do
          End Do
        End Do

        Do k = 1, model_levels - 1
          Do j = 1, rows
            Do i = 1, row_length
              w1(i,j,k) = w(i,j,k) + alpha_2 * R_w(i,j,k)
            End Do
          End Do
        End Do  !  k = 1, model_levels - 1
        Do k = 0, model_levels, model_levels
          Do j = 1, rows
            Do i = 1, row_length
              w1(i,j,k) = 0.
            End Do
          End Do
        End Do  !  k = 0, model_levels, model_levels

! ----------------------------------------------------------------------
! Section 1. set up pointers for both Flux_Rho and Flux_Rho_part
! ----------------------------------------------------------------------

      j_begin = 1
      j_end = rows
      j_start = 1
      j_stop = rows
      i_start = 1
      i_stop = row_length
      L_poles = .false.
 
      If (model_domain  ==  mt_Global) then
        L_poles = .true.
        If (at_extremity(PSouth)) j_begin = 2
        If (at_extremity(PNorth)) j_end = rows-1
      Else If ( model_domain  ==  mt_LAM ) Then 
        If (at_extremity(PSouth)) Then 
          j_start = 2
          j_begin = 2
        End If !  at_extremity(PSouth)
        If (at_extremity(PNorth)) Then
          j_stop = rows - 1
          j_end = rows - 1
        End If  ! at_extremity(PNorth)
        If (at_extremity(PWest)) i_start = 2
        If (at_extremity(PEast)) i_stop = row_length - 1
      Else If ( model_domain  ==  mt_cyclic_LAM ) Then
        If (at_extremity(PSouth)) j_begin = 2
        If (at_extremity(PNorth)) j_end = rows-1        
      End If  !  model_domain  ==  mt_Global

! ----------------------------------------------------------------------
! Section 2. set up terms needed in Flux_Rho row which differ for
!            mixing ratios or specific quantities
! ----------------------------------------------------------------------

! DEPENDS ON: flux_rho_part
        Call Flux_Rho_part(                                             &
     &                     moist, moist_star, moist_np1,                &
     &                     r_theta_levels, r_rho_levels,                &
     &                     rows, row_length,                            &
     &                     model_levels, wet_model_levels,              &
     &                     halo_i, halo_j, off_x, off_y,                &
     &                     alpha_1, j_begin, j_end,                     &
     &                     wet_to_dry_n, wet_to_dry_np1, dry_to_wet,    &
     &                     rho, rho_np1, CycleNo, L_new_tdisc,          &
     &                     L_mix_ratio, L_dry )

! ----------------------------------------------------------------------
! Section 3. Calculate update to rho
! ----------------------------------------------------------------------

! DEPENDS ON: flux_rho
        Call Flux_Rho(                                                  &
     &                 u1, v1, w1,                                      &
     &                 r_theta_levels, r_rho_levels,                    &
     &                 eta_theta_levels, eta_rho_levels,                &
     &                 FV_sec_theta_latitude,                           &
     &                 cos_v_latitude, delta_lambda, delta_phi,         &
     &                 dlambda_p, dphi_p,                               &
     &                 recip_dlambda_u, recip_dphi_v,                   &
     &                 wt_lambda_p, wt_lambda_u,                        &
     &                 wt_phi_p, wt_phi_v,                              &
     &                 timestep, NumCycles, CycleNo,                    &
     &                 rows, n_rows, row_length,                        &
     &                 model_levels, wet_model_levels, model_domain,    &
     &                 first_constant_r_rho_level, rims_to_do,          &
     &                 halo_i, halo_j, nproc, gc_proc_row_group,        &
     &                 L_regular, at_extremity, global_row_length,      &
     &                 off_x, off_y, i_start, i_stop,                   &
     &                 j_start, j_stop, j_begin, j_end,                 &
     &                 wet_to_dry_n, dry_to_wet,                        &
     &                 rho, rho_np1, L_poles, L_new_tdisc )

      Deallocate ( moist_np1 )

      If (L_do_increment) Then    ! calculate increment over timestep

        Do k = 1, model_levels
          Do j = 1, rows
            Do i = 1, row_length
            inc_rho(i,j,k) = rho(i,j,k) - rho_n(i,j,k)
            End Do
          End Do
        End Do  !  k = 1, model_levels

      End If  !  L_do_increment

      IF (lhook) CALL dr_hook('UPDATE_RHO',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE update_rho

