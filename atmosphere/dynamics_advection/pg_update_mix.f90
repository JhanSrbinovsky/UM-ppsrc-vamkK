! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Pg_update_mix
!
      Subroutine Pg_update_mix(                                         &
                               mix_v, mix_cl, mix_cf,                   &
                     mix_cf2, mix_rain, mix_graup,                      &
                               mix_v_star, mix_cl_star, mix_cf_star,    &
                       mix_cf2_star, mix_rain_star, mix_graup_star,     &
                       L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup,           &
                               theta_star, theta,                       &
                               exner,                                   &
                               exner_theta_levels,                      &
                               r_theta_levels, r_rho_levels,            &
                               eta_theta_levels, eta_rho_levels,        &
                               sec_theta_latitude,                      &
                               delta_lambda, delta_phi,                 &
                               lambda_p, phi_p, lambda_u, phi_v,        &
                               timestep, alpha_3, alpha_4,              &
                               Cp,                                      &
                               row_length, rows, n_rows,                &
                               model_levels, wet_model_levels,          &
                               model_domain,                            &
                               first_constant_r_rho_level,              &
                               halo_i, halo_j, off_x, off_y,            &
                               L_regular, at_extremity,                 &
                               CycleNo, L_new_tdisc,                    &
                               R_u, R_v, R_w)

! Purpose:
!          Performs update to pressure gradient term at arrival point
!          to include estimate of theta change between time-levels
!          n and n+1.
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
!          and
!
!          A semi-Implicit scheme for the Unified Model.
!          F.R. Division working paper No 154.
!          M. J. P. Cullen and T. Davies.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE atmos_constants_mod, ONLY: recip_epsilon
      USE domain_params
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, wet_model_levels                                                &
                         ! number of model levels where moisture is held
     &, first_constant_r_rho_level                                      &
                                   ! first rho level on which r
                                   ! is constant.
     &, halo_i                                                          &
                     ! Size of halo in i.
     &, halo_j                                                          &
                     ! Size of halo in j.
     &, off_x                                                           &
                     ! Size of small halo in i
     &, off_y                                                           &
                     ! Size of small halo in j.
     &, CycleNo

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_regular                                                       &
                     ! Variable resolution if false
     &, L_mcr_qcf2, L_mcr_qrain, L_mcr_qgraup                           &
     &, L_new_tdisc

      Real                                                              &
           ! physical constants
        Cp

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
        delta_lambda                                                    &
                         ! grid-length in lambda direction
      , delta_phi        ! grid-length in phi direction

      Real                                                              &
           !VarRes horizontal co-ordinate information
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, lambda_u(1-halo_i:row_length+halo_i)                            &
     &, phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

      Real                                                              &
           ! time-weighting parameters, see WP 154.
     &  alpha_3                                                         &
     &, alpha_4                                                         &
     &, timestep

      Real                                                              &
           ! trigonometric functions
     &  sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)

      Real                                                              &
           ! primary model variables
     &  mix_v  (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_model_levels)                                       &
     &, mix_cl (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_model_levels)                                       &
     &, mix_cf (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,       &
     &          wet_model_levels)                                       &
     &, mix_cf2 (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,      &
     &       wet_model_levels)                                          &
     &, mix_rain (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,     &
     &       wet_model_levels)                                          &
     &, mix_graup (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,    &
     &       wet_model_levels)                                          &
     &, mix_v_star  (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &               wet_model_levels)                                  &
     &, mix_cl_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &               wet_model_levels)                                  &
     &, mix_cf_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,      &
     &               wet_model_levels)                                  &
     &, mix_cf2_star (1-off_x:row_length+off_x,                         &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_rain_star (1-off_x:row_length+off_x,                        &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, mix_graup_star (1-off_x:row_length+off_x,                       &
     &               1-off_y:rows+off_y,wet_model_levels)               &
     &, theta (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)                                            &
     &, theta_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &              model_levels)                                       &
     &, exner (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)                                            &
     &, exner_theta_levels(1-off_x:row_length+off_x, 1-off_y:rows+off_y,&
     &                     model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, eta_theta_levels (0:model_levels)                               &
     &, eta_rho_levels (model_levels)

! Arguments with Intent INOUT. ie: Input variables changed on output.

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
           ! See WP154 and WP162 for definitions.
     &  R_u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &        model_levels)                                             &
     &, R_v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,            &
     &       model_levels)                                              &
     &, R_w (row_length, rows, model_levels)

! Parameters
      Integer                                                           &
     &   PNorth,                                                        &
                      ! North processor address in the neighbor array
     &   PEast,                                                         &
                      ! East processor address in the neighbor array
     &   PSouth,                                                        &
                      ! South processor address in the neighbor array
     &   PWest,                                                         &
                      ! West processor address in the neighbor array
     &   NoDomain     ! Value in neighbor array if the domain has
                      !  no neighbor in this direction. Otherwise
                      !  the value will be the tid of the neighbor
      Parameter (                                                       &
     &   PNorth   = 1,                                                  &
     &   PEast    = 2,                                                  &
     &   PSouth   = 3,                                                  &
     &   PWest    = 4,                                                  &
     &   NoDomain = -1)

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                  ! Loop indices
     &, j0, j1

      Real                                                              &
        recip_delta_lambda                                              &
      , recip_delta_phi

      Real                                                              &
     &  Yw (row_length+1, rows+1, 0:model_levels)                       &
     &, thetav_inc (row_length+1, rows+1, model_levels)                 &
     &, thetav_inc_rho_levels (row_length+1, rows+1, model_levels)      &
     &, ave_vpg(row_length+1, rows+1, first_constant_r_rho_level)       &
     &, recip_dlambda_p(2-halo_i:row_length+halo_i)                     &
     &, recip_dphi_p(2-halo_i:row_length+halo_i, 2-halo_j:rows+halo_j)  &
     &, wt_lambda_p(2-halo_i:row_length+halo_i)                         &
     &, wt_lambda_u(2-halo_i:row_length+halo_i)                         &
     &, wt_phi_p(2-halo_i:row_length+halo_i, 2-halo_j:rows+halo_j)      &
     &, wt_phi_v(2-halo_i:row_length+halo_i, 2-halo_j:n_rows+halo_j)    &
     &, moist (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
     &       wet_model_levels)                                          &
     &, moist_star (1-off_x:row_length+off_x,                           &
     &               1-off_y:rows+off_y,wet_model_levels)

      Real                                                              &
     &  interp2 (row_length, rows)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! No External Routines:

! ----------------------------------------------------------------------
! Section 0.  Calculate thetav_inc and thetav_inc on rho levels.
! ----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('PG_UPDATE_MIX',zhook_in,zhook_handle)
      moist =1 + mix_v + mix_cl + mix_cf
      moist_star = 1 + mix_v_star + mix_cl_star + mix_cf_star
      If(L_mcr_qcf2)then
        moist      = moist + mix_cf2
        moist_star = moist_star + mix_cf2_star
      endif
      If(L_mcr_qrain)then
        moist      = moist + mix_rain
        moist_star = moist_star + mix_rain_star
      endif
      If(L_mcr_qgraup)then
        moist      = moist + mix_graup
        moist_star = moist_star + mix_graup_star
      endif

      j0 = 1
      j1 = rows
      If (model_domain  /=  mt_bi_cyclic_lam) then
      If (at_extremity(PSouth)) j0 = 2
      If (at_extremity(PNorth)) j1 = rows-1
      Endif
      if (L_regular ) then
        recip_delta_lambda = 1./delta_lambda
        recip_delta_phi = 1./delta_phi
      else  ! variable resolution
        do i = 2-halo_i, row_length+halo_i
          recip_dlambda_p(i) =  1.0 /  (lambda_p(i) - lambda_p(i-1))
          wt_lambda_p(i) = lambda_p(i) - lambda_u(i-1)
          wt_lambda_u(i) = lambda_u(i) - lambda_p(i)
        enddo
        do j = 2-halo_j, rows+halo_j
          do i = 1-halo_i, row_length+halo_i
            recip_dphi_p(i,j) = 1.0 / ( phi_p(i,j) - phi_p(i,j-1))
            wt_phi_p(i,j) = phi_p(i,j) - phi_v(i,j-1)
          enddo
        enddo
        do j = 2-halo_j, n_rows+halo_j
          do i = 1-halo_i, row_length+halo_i
            wt_phi_v(i,j) = phi_v(i,j) - phi_p(i,j)
          enddo
        enddo
      endif ! L_regular


      Do k = 1, wet_model_levels
        Do j = 1, j1+1
          Do i = 1, row_length+1
              thetav_inc(i,j,k) =  ( theta_star(i,j,k) *                &
     &        ( 1. +   recip_epsilon * mix_v_star(i,j,k)  ) /           &
     &        moist_star(i,j,k) )                                       &
     &                          -  ( theta(i,j,k) *                     &
     &        ( 1. +   recip_epsilon * mix_v(i,j,k)  ) /                &
     &        moist(i,j,k) )
          End Do
        End Do
      End Do
      Do k = 1 + wet_model_levels, model_levels
        Do j = 1, j1+1
          Do i = 1, row_length+1
            thetav_inc(i,j,k) = theta_star(i,j,k) - theta(i,j,k)
          End Do
        End Do
      End Do

! calculate thetav_inc * vertical derivative at all points.
! store in Yw to save re-calculation.
      Do k = 1, model_levels - 1
        Do j = 1, j1+1
          Do i = 1, row_length+1
            Yw(i,j,k) = Cp * thetav_inc(i,j,k) *                        &
     &                 (exner(i,j,k+1) - exner(i,j,k)) /                &
     &                 (r_rho_levels(i,j,k+1) -                         &
     &                  r_rho_levels(i,j,k))
          End Do
        End Do
      End Do
      k = 0
        Do j = 1, j1+1
          Do i = 1, row_length+1
            Yw(i,j,k) = 0
          End Do
        End Do

! If dynamics are iterated no vpg term update is necessary for u, v
! as this has taken place in u, v SL advection routine slfullw2a.dk
      If ( CycleNo == 1 .OR. (.NOT. L_new_tdisc) ) Then

        Do k = 1, first_constant_r_rho_level - 1
! ave_vpg term at rho levels
          Do j = 1, j1+1
            Do i = 1, row_length+1
              ave_vpg(i,j,k) = ((r_rho_levels(i,j,k)                    &
     &                           - r_theta_levels(i,j,k-1))             &
     &                           * Yw(i,j,k) +                          &
     &                          (r_theta_levels(i,j,k)                  &
     &                           - r_rho_levels(i,j,k))                 &
     &                           * Yw(i,j,k-1) )                        &
     &                        / ( r_theta_levels(i,j,k) -               &
     &                            r_theta_levels(i,j,k-1) )
            End Do
          End Do
        End Do

        k = 1
        Do j = 1, j1+1
          Do i = 1, row_length+1
            thetav_inc_rho_levels(i,j,k) = thetav_inc(i,j,k)
          End Do
        End Do
        Do k = 2, model_levels
        Do j = 1, j1+1
          Do i = 1, row_length+1
! thetav_inc at rho levels
            thetav_inc_rho_levels(i,j,k) = ((r_rho_levels(i,j,k)        &
     &                                   - r_theta_levels(i,j,k-1))     &
     &                                  * thetav_inc(i,j,k) +           &
     &                                  (r_theta_levels(i,j,k)          &
     &                                   - r_rho_levels(i,j,k))         &
     &                                  * thetav_inc(i,j,k-1) )         &
     &                                 / ( r_theta_levels(i,j,k) -      &
     &                                     r_theta_levels(i,j,k-1) )
          End Do
        End Do
        End Do

! ----------------------------------------------------------------------
! Section 1.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 1.1  Calculate new values of R_u.
! ----------------------------------------------------------------------

        Do k = 1, model_levels

! Calculate pressure gradient term and store in interp2.

        If ( k >= first_constant_r_rho_level ) Then
        if (L_regular ) then
          Do j = j0, j1
            Do i = 1, row_length
              interp2(i,j) =  Cp * (thetav_inc_rho_levels(i,j,k) +      &
     &                              thetav_inc_rho_levels(i+1,j,k))     &
     &                     * (exner(i+1,j,k) - exner(i,j,k)) *          &
     &                          recip_delta_lambda                      &
     &                     * sec_theta_latitude(i,j)/                   &
     &                       (r_rho_levels(i,j,k) +                     &
     &                        r_rho_levels(i+1,j,k))
            End Do
          End Do

        else  ! variable resolution

          Do j = j0, j1
            Do i = 1, row_length
              interp2(i,j) =  Cp * (thetav_inc_rho_levels(i,j,k) *      &
     &                                            wt_lambda_p(i+1) +    &
     &                              thetav_inc_rho_levels(i+1,j,k) *    &
     &                                            wt_lambda_u(i) ) *    &
     &                             (exner(i+1,j,k) - exner(i,j,k)) *    &
     &                                        recip_dlambda_p(i+1) *    &
     &                                     sec_theta_latitude(i,j) /    &
     &                    (r_rho_levels(i,j,k) *  wt_lambda_p(i+1) +    &
     &                     r_rho_levels(i+1,j,k) * wt_lambda_u(i) )
            End Do
          End Do

        endif ! L_regular

        Else ! k < first_constant_r_rho_level

        if (L_regular ) then

          Do j = j0, j1
            Do i = 1, row_length
              interp2(i,j) = ( Cp * (thetav_inc_rho_levels(i,j,k) +     &
     &                               thetav_inc_rho_levels(i+1,j,k))    &
     &                         * (exner(i+1,j,k) - exner(i,j,k))        &
     &                         - (r_rho_levels(i+1,j,k) -               &
     &                            r_rho_levels(i,j,k) ) *               &
     &                         (ave_vpg(i+1,j,k) + ave_vpg(i,j,k)) )    &
     &                     * recip_delta_lambda                         &
     &                     * sec_theta_latitude(i,j)/                   &
     &                       (r_rho_levels(i,j,k) +                     &
     &                        r_rho_levels(i+1,j,k))
            End Do
          End Do

        else  ! variable resolution

          Do j = j0, j1
            Do i = 1, row_length
              interp2(i,j) = ( Cp * (thetav_inc_rho_levels(i,j,k) *     &
     &                                            wt_lambda_p(i+1) +    &
     &                               thetav_inc_rho_levels(i+1,j,k) *   &
     &                                             wt_lambda_u(i) ) *   &
     &                              (exner(i+1,j,k) - exner(i,j,k)) -   &
     &               (r_rho_levels(i+1,j,k) - r_rho_levels(i,j,k) ) *   &
     &                         (ave_vpg(i,j,k) * wt_lambda_p(i+1) +     &
     &                          ave_vpg(i+1,j,k) * wt_lambda_u(i) ) ) * &
     &                                           recip_dlambda_p(i+1) * &
     &                                        sec_theta_latitude(i,j) / &
     &                        (r_rho_levels(i,j,k) * wt_lambda_p(i+1) + &
     &                         r_rho_levels(i+1,j,k) * wt_lambda_u(i) )
            End Do
          End Do

        endif ! L_regular

        End If !  k >= first_constant_r_rho_level


! Modify R_u from the calculated terms.
        Do j = j0, j1
          Do i = 1, row_length
            R_u(i,j,k) = R_u(i,j,k)- interp2(i,j) * alpha_3 * timestep
          End Do
        End Do

! End loop over levels
        End Do

! ----------------------------------------------------------------------
! Section 1.2  Calculate new values of R_v.
! ----------------------------------------------------------------------

        Do k = 1, model_levels

        If ( k < first_constant_r_rho_level ) Then
        if (L_regular ) then
          Do j = 1, n_rows
            Do i = 1, row_length
              R_v(i,j,k) = R_v(i,j,k) - alpha_3 * timestep *            &
     &                    ( Cp * (thetav_inc_rho_levels(i,j+1,k) +      &
     &                            thetav_inc_rho_levels(i,j,k))         &
     &                         * (exner(i,j+1,k) - exner(i,j,k))        &
     &                         - (r_rho_levels(i,j+1,k) -               &
     &                            r_rho_levels(i,j,k) ) *               &
     &                         (ave_vpg(i,j+1,k) + ave_vpg(i,j,k)) )    &
     &                     * recip_delta_phi /                          &
     &                       (r_rho_levels(i,j,k) +                     &
     &                        r_rho_levels(i,j+1,k))
            End Do
          End Do

        else  ! variable resolution

          Do j = 1, n_rows
            Do i = 1, row_length
              R_v(i,j,k) = R_v(i,j,k) - alpha_3 * timestep *            &
     &                          ( Cp * (thetav_inc_rho_levels(i,j,k) *  &
     &                                             wt_phi_p(i,j+1)   +  &
     &                                thetav_inc_rho_levels(i,j+1,k) *  &
     &                                               wt_phi_v(i,j) ) *  &
     &                               (exner(i,j+1,k) - exner(i,j,k)) -  &
     &                (r_rho_levels(i,j+1,k) -  r_rho_levels(i,j,k)) *  &
     &                             (ave_vpg(i,j,k) * wt_phi_p(i,j+1) +  &
     &                          ave_vpg(i,j+1,k) * wt_phi_v(i,j) ) ) *  &
     &                                           recip_dphi_p(i,j+1) /  &
     &                        (r_rho_levels(i,j,k) * wt_phi_p(i,j+1) +  &
     &                         r_rho_levels(i,j+1,k) * wt_phi_v(i,j) )
            End Do
          End Do

        endif ! L_regular

        Else !  k >= first_constant_r_rho_level

        if (L_regular ) then

          Do j = 1, n_rows
            Do i = 1, row_length
              R_v(i,j,k) = R_v(i,j,k) - alpha_3 * timestep *            &
     &                     Cp * (thetav_inc_rho_levels(i,j,k) +         &
     &                           thetav_inc_rho_levels(i,j+1,k))        &
     &                     * (exner(i,j+1,k) - exner(i,j,k)) *          &
     &                          recip_delta_phi /                       &
     &                       (r_rho_levels(i,j+1,k) +                   &
     &                        r_rho_levels(i,j,k))
            End Do
          End Do
        else  ! variable resolution

          Do j = 1, n_rows
            Do i = 1, row_length
              R_v(i,j,k) = R_v(i,j,k) - alpha_3 * timestep *            &
     &                             Cp * (thetav_inc_rho_levels(i,j,k) * &
     &                                                wt_phi_p(i,j+1) + &
     &                                 thetav_inc_rho_levels(i,j+1,k) * &
     &                                                wt_phi_v(i,j) ) * &
     &                                (exner(i,j+1,k) - exner(i,j,k)) * &
     &                                            recip_dphi_p(i,j+1) / &
     &                        ( r_rho_levels(i,j,k) * wt_phi_p(i,j+1) + &
     &                          r_rho_levels(i,j+1,k) * wt_phi_v(i,j) )
            End Do
          End Do

        endif ! L_regular

        End If !   k < first_constant_r_rho_level


! End loop over levels
        End Do

      End If ! CycleNo == 1 .OR. .NOT. L_new_tdisc

! ----------------------------------------------------------------------
! Section 1.3  Calculate new values of R_w.
! ----------------------------------------------------------------------

      Do k = 1, model_levels-1
        Do j = 1, rows
          Do i = 1, row_length
            R_w(i,j,k) = R_w(i,j,k)                                     &
     &                   - alpha_4 * timestep * Yw(i,j,k)
          End Do
        End Do
      End Do

! End of routine.
      IF (lhook) CALL dr_hook('PG_UPDATE_MIX',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Pg_update_mix

