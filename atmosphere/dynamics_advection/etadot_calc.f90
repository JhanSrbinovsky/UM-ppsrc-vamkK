! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Etadot_Calc.

      Subroutine Etadot_Calc(                                           &
     &                       r_theta_levels, r_rho_levels,              &
     &                       eta_theta_levels, eta_rho_levels,          &
     &                       u, v, w,                                   &
     &                       sec_theta_latitude,                        &
     &                       row_length, rows, n_rows, model_levels,    &
     &                       delta_lambda, delta_phi,                   &
     &                       dlambda_p, dphi_p,                         &
     &                       wt_lambda_p, wt_lambda_u,                  &
     &                       wt_phi_p, wt_phi_v,                        &
     &                       model_domain, first_constant_r_rho_level,  &
     &                       proc_row_group, at_extremity, g_row_length,&
     &                       off_x, off_y, halo_i, halo_j,              &
     &                       halo_i_wind, halo_j_wind,                  &
     &                       halo_i_w, halo_j_w, halo_i_eta, halo_j_eta,&
     &                       i_start, i_stop, j_start, j_stop,          &
     &                       j_begin, j_end, L_regular, L_poles,        &
     &                       l_s_poles, l_n_poles, etadot)

! Purpose:
!          Calculates etadot from given u, v and w fields.
!
! Method:
!          Is described in ;
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!


      USE global_2d_sums_mod, ONLY : global_2d_sums
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Integer                                                           &
              ! model dimensions
     &  row_length                                                      &
                         ! number of points on a row
     &, g_row_length                                                    &
                         ! global number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels
     &, off_x                                                           &
     &, off_y                                                           &
     &, halo_i                                                          &
     &, halo_j                                                          &
     &, halo_i_wind                                                     &
                     ! Size of halo in i for u, v
     &, halo_j_wind                                                     &
                     ! Size of halo in j for u, v
     &, halo_i_w                                                        &
                     ! Size of halo in i for w
     &, halo_j_w                                                        &
                     ! Size of halo in j for w
     &, halo_i_eta                                                      &
                     ! Size of halo in i for w
     &, halo_j_eta                                                      &
                     ! Size of halo in j for w
! Pointers for domain calculation start
     &, i_start                                                         &
     &, i_stop                                                          &
     &, j_start                                                         &
     &, j_stop                                                          &
     &, j_begin                                                         &
     &, j_end                                                           &
! Pointers for domain calculation end
     &, proc_row_group   ! Group id for processors on the same row

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_poles                                                         &
                         ! true for etadot at poles
     &, L_regular        ! true for regular resolution


      Integer                                                           &
     &  first_constant_r_rho_level ! first rho level on which r
                                   ! is constant.

      Real                                                              &
           ! horizontal co-ordinate information
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
            ! VarRes horizontal co-ordinate information
     &  dlambda_p(1-halo_i:row_length+halo_i)                           &
     &, dphi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)        &
     &, wt_lambda_p(1-halo_i:row_length+halo_i)                         &
     &, wt_lambda_u(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)      &
     &, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

      Real                                                              &
           ! trigonometric functions
     &  sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)

      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

! Primary Arrays

      Real                                                              &
     &  u(1-halo_i_wind:row_length+halo_i_wind,                         &
     &    1-halo_j_wind:rows+halo_j_wind, model_levels)                 &
     &, v(1-halo_i_wind:row_length+halo_i_wind,                         &
     &    1-halo_j_wind:n_rows+halo_j_wind, model_levels)               &
     &, w(1-halo_i_w:row_length+halo_i_w, 1-halo_j_w:rows+halo_j_w,     &
     &     0:model_levels)

! Arguments with Intent OUT. ie: Output variables.

      Real                                                              &
     &  etadot(1-halo_i_eta:row_length+halo_i_eta,                      &
               1-halo_j_eta:rows+halo_j_eta, 0:model_levels)            &
! Optional output variables which can be used in flux_rho
     &, l_s_poles(row_length, first_constant_r_rho_level-1)             &
     &, l_n_poles(row_length, first_constant_r_rho_level-1)

! Local Variables.

! scalars

      Integer                                                           &
     &  i, j, k   ! Loop indices

      Real                                                              &
     &  Recip_delta_lambda                                              &
     &, Recip_delta_phi                                                 &
     &, weight1                                                         &
     &, weight2

      Integer info

! arrays

      Real                                                              &
     &  u_on_theta(0:row_length, rows)                                  &
     &, v_on_theta(row_length, 0:rows)                                  &
     &, term(0:row_length, 0:rows)                                      &
     &, sum_s(model_levels)                                             &
     &, sum_n(model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! No External Routines:

! Functions: None

! ----------------------------------------------------------------------
! Section 1.   Interpolate u and v to w points.
!              Not done at top level.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('ETADOT_CALC',zhook_in,zhook_handle)
      If ( L_regular ) Then
        recip_delta_lambda = 1./ delta_lambda
        recip_delta_phi = 1./ delta_phi
      end If ! L_regular
      
!  Initialise to zero everywhere
      etadot = 0.0

      Do k = 1, first_constant_r_rho_level-1

        weight2 = (eta_rho_levels(k+1) - eta_theta_levels(k) )/         &
     &            (eta_rho_levels(k+1) - eta_rho_levels(k) )
        weight1 = (eta_theta_levels(k) - eta_rho_levels(k))/            &
     &            (eta_rho_levels(k+1) - eta_rho_levels(k) )

!  calculate u on theta levels.

        Do j = j_begin, j_end
          Do i = i_start - 1, i_stop
            u_on_theta (i,j) = weight2 * u(i,j,k) +                     &
     &                         weight1 * u(i,j,k+1)
          End Do
        End Do

!  calculate v on theta levels.
        Do j = j_begin - 1, j_end
          Do i = i_start, i_stop
            v_on_theta (i,j) = weight2 * v(i,j,k) +                     &
     &                         weight1 * v(i,j,k+1)
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2.   Calculate u.grad(r)
! ----------------------------------------------------------------------

! Calculate u dr/d lambda.
        If ( L_regular ) Then

        Do j = j_begin, j_end
          Do i = i_start - 1, i_stop
            term(i,j) = u_on_theta(i,j) *                               &
     &                 (r_theta_levels(i+1,j,k) -                       &
     &                  r_theta_levels(i,j,k) ) * recip_delta_lambda    &
     &                  * 2.0 / (r_theta_levels(i+1,j,k) +              &
     &                           r_theta_levels(i,j,k) )
          End Do
        End Do

        else  ! variable resolution

        Do j = j_begin, j_end
          Do i = i_start - 1, i_stop
            term(i,j) = u_on_theta(i,j) *                               &
     &           ( r_theta_levels(i+1,j,k) -  r_theta_levels(i,j,k) ) / &
     &                   ( ( wt_lambda_p(i+1) * r_theta_levels(i,j,k) + &
     &            (1.0 - wt_lambda_p(i+1)) * r_theta_levels(i+1,j,k) ) *&
     &                                                  dlambda_p(i) )
          End Do
        End Do

        endIf ! L_regular

! average quantity and store in etadot
        If ( L_regular ) Then
        Do j = j_begin, j_end
          Do i = i_start, i_stop
            etadot(i,j,k) = (term(i,j) + term(i-1,j) ) * 0.5 *          &
     &                        sec_theta_latitude(i,j)
          End Do
        End Do
        else  ! variable resolution
        Do j = j_begin, j_end
          Do i = i_start, i_stop
             etadot(i,j,k) = (  wt_lambda_u(i) * term(i-1,j) +          &
    &                        (1.0 - wt_lambda_u(i)) * term(i,j) ) *     &
    &                         sec_theta_latitude(i,j)
          End Do
        End Do
        endIf ! L_regular

! Calculate v dr/d phi.

        If ( L_regular ) Then
        Do j = j_begin-1, j_end 
          Do i = i_start, i_stop
            term(i,j) = v_on_theta(i,j) *                               &
     &                 (r_theta_levels(i,j+1,k) -                       &
     &                  r_theta_levels(i,j,k) ) * recip_delta_phi       &
     &                    * 2.0 / (r_theta_levels(i,j+1,k) +            &
     &                             r_theta_levels(i,j,k) )
          End Do
        End Do
        else  ! variable resolution
        Do j = j_begin-1, j_end
          Do i = i_start, i_stop
            term(i,j) = v_on_theta(i,j) *                               &
     &             ( r_theta_levels(i,j+1,k) - r_theta_levels(i,j,k) ) /&
     &                     ( ( wt_phi_p(i,j+1) * r_theta_levels(i,j,k) +&
     &             (1.0 - wt_phi_p(i,j+1)) * r_theta_levels(i,j+1,k) ) *&
     &                         dphi_p(i,j+1))
          End Do
        End Do
        endIf ! L_regular

! average quantity and add on to u term
        If ( L_regular ) Then
        Do j = j_begin, j_end
          Do i = i_start, i_stop
            etadot(i,j,k) = (term(i,j) + term(i,j-1)) * 0.5 +           &
     &                        etadot(i,j,k)
          End Do
        End Do
        else  ! variable resolution
        Do j = j_begin, j_end
          Do i = i_start, i_stop
             etadot(i,j,k) =  wt_phi_v(i,j) * term(i,j-1) +             &
     &                      (1.0 - wt_phi_v(i,j)) * term(i,j) +         &
     &                                        etadot(i,j,k)
          End Do
        End Do
        endIf ! L_regular

! northern and southern boundary updates.
! save values for summing.

        If (L_poles) Then

          If(at_extremity(PSouth))then
            Do i = 1, row_length
              l_s_poles(i,k) = term(i,1)
            End Do
          End If

          If(at_extremity(PNorth))then
            Do i =1, row_length
              l_n_poles(i,k) = term(i,n_rows)
            End Do
          End If

        End If  !  L_poles

      End Do   !  k = 1, first_constant_r_rho_level - 1

! if some levels need polar meaning
      If (first_constant_r_rho_level > 1 .and. L_poles ) Then

! calculate average value at poles
        If (at_extremity(PSouth)) Then

          CALL global_2d_sums(l_s_poles, row_length, 1, 0, 0,           &
                              first_constant_r_rho_level-1,             &
                              sum_s, proc_row_group)

          Do k = 1, first_constant_r_rho_level-1
            sum_s(k) = sum_s(k) / g_row_length
            Do i = 1, row_length
              etadot(i,1,k) = sum_s(k) * 2.0
            End Do
          End Do

        End If  !  at_extremity(PSouth)

        If (at_extremity(PNorth)) Then

          CALL global_2d_sums(l_n_poles, row_length, 1, 0, 0,           &
                              first_constant_r_rho_level-1,             &
                              sum_n, proc_row_group)

          Do k = 1, first_constant_r_rho_level-1
            sum_n(k) = sum_n(k) / g_row_length
            Do i = 1, row_length
              etadot(i,rows,k) = sum_n(k) * 2.0
            End Do
          End Do

        End If  !  at_extremity(PNorth)

      End If   !  first_constant_r_rho_level > 1 .and. L_poles

! calculate final etadot value

      Do k = 1, first_constant_r_rho_level-1
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            etadot(i,j,k) = (w(i,j,k) - etadot(i,j,k))                  &
     &                       * (eta_rho_levels(k+1) -                   &
     &                          eta_rho_levels(k))                      &
     &                       / (r_rho_levels(i,j,k+1) -                 &
     &                          r_rho_levels(i,j,k))
          End Do
        End Do
      End Do  ! k = 1, first_constant_r_rho_level-1

      Do k = first_constant_r_rho_level, model_levels - 1
! grad r is zero so just perform vertical re-scaling.
        Do j = j_start, j_stop
          Do i = i_start, i_stop
            etadot(i,j,k) = w(i,j,k)                                    &
     &                       * (eta_rho_levels(k+1) -                   &
     &                          eta_rho_levels(k))                      &
     &                       / (r_rho_levels(i,j,k+1) -                 &
     &                          r_rho_levels(i,j,k))
          End Do
        End Do
      End Do  !  k = first_constant_r_rho_level, model_levels - 1

! etadot has already been set to zero at top and bottom

      IF (lhook) CALL dr_hook('ETADOT_CALC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Etadot_Calc
