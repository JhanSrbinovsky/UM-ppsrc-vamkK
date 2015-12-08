! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Flux_rho_mix
!

      Subroutine Flux_rho_mix(                                          &
     &                    u1, v1, w1, moist, moist_star, moist_np1,     &
     &                    r_theta_levels, r_rho_levels,                 &
     &                    eta_theta_levels, eta_rho_levels,             &
     &                    FV_sec_theta_latitude,                        &
     &                    cos_v_latitude, delta_lambda, delta_phi,      &
     &                    lambda_p, phi_p, lambda_u, phi_v,             &
     &                    wt_lambda_u, wt_phi_v,                        &
     &                    recip_dlambda_u, recip_dphi_v,                &
     &                    timestep, NumCycles, CycleNo,                 &
     &                    rows, n_rows, row_length,                     &
     &                    model_levels, wet_model_levels, model_domain, &
     &                    first_constant_r_rho_level, rims_to_do,       &
     &                    halo_i, halo_j, n_proc, proc_row_group,       &
     &                    L_regular, at_extremity, g_row_length,        &
     &                    off_x, off_y,halo_i_wind, halo_j_wind,        &
     &                    wet_to_dry_n, wet_to_dry_np1, alpha_1,        &
     &                    rho, rho_np1, L_new_tdisc )
! Purpose:
!          Calculates density rho at new time level using the flux form
!          of the continuity equation.
!
! Method:
!          Is described in ;
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


      USE global_2d_sums_mod, ONLY : global_2d_sums
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

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
     &, wet_model_levels                                                &
                         ! number of model levels where moisture
                         ! variables are held.
     &, halo_i                                                          &
                     ! Size of halo in i.
     &, halo_j                                                          &
                     ! Size of halo in j.
     &, off_x                                                           &
     &, off_y                                                           &
     &, halo_i_wind                                                     &
                     ! Size of halo in i for wind fields.
     &, halo_j_wind                                                     &
                     ! Size of halo in j for wind fields.
     &, n_proc                                                          &
                     ! Total number of processors
     &, proc_row_group                                                  &
                       ! Group id for processors on the same row
     &, rims_to_do                                                      &
                        !  lbc weights = 1 zone
     &, NumCycles                                                       &
     &, CycleNo

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_regular                                                       &
                    ! false if variable resolution
     &, L_new_tdisc

      Integer                                                           &
     &  first_constant_r_rho_level ! first rho level on which r
                                   ! is constant.

      Integer                                                           &
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
     &  timestep                                                        &
     &, alpha_1

      Real                                                              &
     &  u1(1-halo_i_wind:row_length+halo_i_wind,                        &
     &     1-halo_j_wind:rows+halo_j_wind, model_levels)                &
     &, v1(1-halo_i_wind:row_length+halo_i_wind,                        &
     &     1-halo_j_wind:n_rows+halo_j_wind, model_levels)              &
     &, w1(1-halo_i_wind:row_length+halo_i_wind,                        &
     &     1-halo_j_wind:rows+halo_j_wind, 0:model_levels)              &
     &, moist (1-halo_i_wind:row_length+halo_i_wind,                    &
     &     1-halo_j_wind:rows+halo_j_wind, wet_model_levels)            &
     &, moist_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &          wet_model_levels)                                       &
     &, moist_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &          wet_model_levels)                                       &
     &, wet_to_dry_n (1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
     &          model_levels)                                           &
     &, wet_to_dry_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,   &
     &          model_levels)

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
           ! VarRes horizontal co-ordinate spacing.
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, phi_p(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, lambda_u(1-halo_i:row_length+halo_i)                            &
     &, phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)       &
     &, wt_lambda_u(1-halo_i:row_length+halo_i)                         &
     &, wt_phi_v(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)    &
     &, recip_dlambda_u(1-halo_i:row_length+halo_i)                     &
     &, recip_dphi_v(1-halo_i:row_length+halo_i,1-halo_j:n_rows+halo_j)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  rho (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &       model_levels)                                              &
     &, rho_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &       model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                     ! Loop indices
     &, j0, j1                                                          &
     &, i_start, i_end                                                  &
     &, j_start, j_end                                                  &
     &, levels

      Integer info

      Real                                                              &
     &  weight1                                                         &
     &, weight2                                                         &
     &, weight3                                                         &
     &, temp

      Real                                                              &
     &  recip_delta_lambda                                              &
     &, recip_delta_phi

! Local arrays
      Real                                                              &
     &  dry_to_wet(0:row_length+1, 0:rows+1, wet_model_levels+1)        &
     &, increment(row_length, rows, model_levels)                       &
     &, Dr_Deta_rho_levels(0:row_length+1, 0:rows+1, model_levels)      &
     &, etadot1 (1-halo_i_wind:row_length+halo_i_wind,                  &
     &           1-halo_j_wind:rows+halo_j_wind, 0:model_levels)        &
     &, l_n_poles(row_length,model_levels)                              &
     &, l_s_poles(row_length,model_levels)                              &
     &, sum_n(model_levels)                                             &
     &, sum_s(model_levels)

      Real                                                              &
     &  interp1(0:row_length+1, 0:rows+1)                               &
     &, interp2(row_length, rows)                                       &
     &, lambda_pmiu(1-halo_i:row_length+halo_i)                         &
     &, lambda_umip(1-halo_i:row_length+halo_i)                         &
     &, phi_pmiv(1-halo_j:rows+halo_j)                                  &
     &, phi_vmip(1-halo_j:n_rows+halo_j)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External Routines:
      External Etadot_Calc

! ----------------------------------------------------------------------
! Section 1.   Convert density to dry density.
!              Calculate etadot1.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('FLUX_RHO_MIX',zhook_in,zhook_handle)
      j0 = 1
      j1 = rows
      If (model_domain  /=  mt_bi_cyclic_lam) then
      If (at_extremity(PSouth)) j0 = 2
      If (at_extremity(PNorth)) j1 = rows-1
      Endif
      j_start = 1
      j_end = rows
      i_start = 1
      i_end = row_length
      If (model_domain  ==  mt_LAM) Then
        If (at_extremity(PSouth)) j_start = 2
        If (at_extremity(PNorth)) j_end = rows - 1
        If (at_extremity(PWest)) i_start = 2
        If (at_extremity(PEast)) i_end = row_length - 1
      End If

      If ( L_regular ) Then
        recip_delta_lambda = 1./ delta_lambda
        recip_delta_phi = 1./ delta_phi
      else  ! variable resolution
! polar value for recip_delta_phi
      If (at_extremity(PSouth))                                         &
     &   recip_delta_phi = 1.0 / (phi_v(1,2) - phi_v(1,1))
      If (at_extremity(PNorth))                                         &
     &   recip_delta_phi = 1.0 / (phi_v(1,n_rows) - phi_v(1,n_rows-1))
        Do j = 1-off_y, rows+off_y
          phi_pmiv(j) = phi_p(1,j) - phi_v(1,j-1)
        enddo
        Do j = 1-off_y, n_rows+off_y
          phi_pmiv(j) = phi_p(1,j) - phi_v(1,j-1)
          phi_vmip(j) = phi_v(1,j) - phi_p(1,j)
        enddo
        Do i = 1-off_x, row_length+off_x
          lambda_pmiu(i) = lambda_p(i) - lambda_u(i-1)
          lambda_umip(i) = lambda_u(i) - lambda_p(i)
        enddo
      endIf ! L_regular


      If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then

      k = 1
        Do j = j0-1, j1+1
          Do i = 0, row_length+1
            dry_to_wet(i,j,k) = 1. + moist_star(i,j,k)
           wet_to_dry_np1(i,j,k) = 1. /dry_to_wet(i,j,k)
! rho holds rho_moist, convert to rho_dry
           wet_to_dry_n(i,j,k) = 1. /(1. + moist(i,j,k))
           rho(i,j,k) = rho(i,j,k)*wet_to_dry_n(i,j,k)
          End Do
        End Do

      Do k = 2, wet_model_levels
        Do j = j0-1, j1+1
          Do i = 0, row_length+1
            weight2 = r_rho_levels(i,j,k)                               &
     &                - r_theta_levels(i,j,k-1)
            weight1 = r_theta_levels(i,j,k)                             &
     &                - r_rho_levels(i,j,k)
            weight3 = r_theta_levels(i,j,k)                             &
     &                - r_theta_levels(i,j,k-1)
            dry_to_wet(i,j,k) = ( weight2 * (1. + moist_star(i,j,k)) +  &
     &                            weight1 * (1. + moist_star(i,j,k-1)) )&
     &                          / weight3
           wet_to_dry_np1(i,j,k) = 1. /dry_to_wet(i,j,k)
            temp = ( weight2 * (1. + moist(i,j,k)) +                    &
     &               weight1 * (1. + moist(i,j,k-1)) )                  &
     &             / weight3
! rho holds rho_moist, convert to rho_dry
           wet_to_dry_n(i,j,k) = 1. /temp
           rho(i,j,k) = rho(i,j,k) * wet_to_dry_n(i,j,k)

          End Do
        End Do
      End Do

      k = wet_model_levels + 1
      if ( k <= model_levels ) then
        Do j = j0-1, j1+1
          Do i = 0, row_length+1
            weight2 = r_rho_levels(i,j,k)                               &
     &                - r_theta_levels(i,j,k-1)
            weight1 = r_theta_levels(i,j,k)                             &
     &                - r_rho_levels(i,j,k)
            weight3 = r_theta_levels(i,j,k)                             &
     &                - r_theta_levels(i,j,k-1)

            dry_to_wet(i,j,k) = ( weight2  +                            &
     &                            weight1 * (1. + moist_star(i,j,k-1)) )&
     &                          / weight3
           wet_to_dry_np1(i,j,k) = 1. /dry_to_wet(i,j,k)
            temp = ( weight2  +                                         &
     &               weight1 * (1. + moist(i,j,k-1)) )                  &
     &             / weight3
! rho holds rho_moist, convert to rho_dry
           wet_to_dry_n(i,j,k) = 1. /temp
           rho(i,j,k) = rho(i,j,k) * wet_to_dry_n(i,j,k)
          End Do
        End Do

      endif !  k <= model_levels

      Else
        k=1
        Do j = j0-1, j1+1
          Do i = 0, row_length+1
            dry_to_wet(i,j,k)    =  1.+ moist_star(i,j,k)
            wet_to_dry_np1(i,j,k)=  1./dry_to_wet(i,j,k)
! rho holds rho_moist, convert to rho_dry
            wet_to_dry_n(i,j,k)  =  1./(1. + moist(i,j,k))
            rho(i,j,k)= rho(i,j,k)*wet_to_dry_n(i,j,k)
! compute t-averaged dry rho using last iteration's
! rho^{n+1} estimate
            rho_np1(i,j,k)= (1.-alpha_1)*rho(i,j,k)                     &
     &                    + alpha_1*rho_np1(i,j,k)/(1.+moist_np1(i,j,k))
          End Do
        End Do

        Do k = 2, wet_model_levels
          Do j = j0-1, j1+1
            Do i = 0, row_length+1
              weight2 = r_rho_levels(i,j,k)                             &
     &                - r_theta_levels(i,j,k-1)
              weight1 = r_theta_levels(i,j,k)                           &
     &                - r_rho_levels(i,j,k)
              weight3 = r_theta_levels(i,j,k)                           &
     &                - r_theta_levels(i,j,k-1)
              dry_to_wet(i,j,k) = ( weight2*(1.+moist_star(i,j,k)) +    &
     &                            weight1*(1.+moist_star(i,j,k-1)) )    &
     &                          / weight3
              wet_to_dry_np1(i,j,k) = 1./dry_to_wet(i,j,k)
              temp = ( weight2 * (1. + moist(i,j,k)) +                  &
     &                 weight1 * (1. + moist(i,j,k-1)) )                &
     &               / weight3
! rho holds rho_moist, convert to rho_dry
              wet_to_dry_n(i,j,k) = 1. /temp
              rho(i,j,k) = rho(i,j,k) * wet_to_dry_n(i,j,k)
! compute t-averaged dry rho using last iteration's
! rho^{n+1} estimate
              temp = ( weight2 * (1. + moist_np1(i,j,k)) +              &
     &                 weight1 * (1. + moist_np1(i,j,k-1)) )            &
     &               / weight3
              rho_np1(i,j,k) = (1.-alpha_1)*rho(i,j,k)                  &
     &                       + alpha_1*rho_np1(i,j,k)/temp
            End Do
          End Do
        End Do

        k = wet_model_levels+1
        If ( k <= model_levels ) Then
          Do j = j0-1, j1+1
            Do i = 0, row_length+1

              weight2 = r_rho_levels(i,j,k)                             &
     &                - r_theta_levels(i,j,k-1)
              weight1 = r_theta_levels(i,j,k)                           &
     &                - r_rho_levels(i,j,k)
              weight3 = r_theta_levels(i,j,k)                           &
     &                - r_theta_levels(i,j,k-1)
              dry_to_wet(i,j,k) = ( weight2  +                          &
     &                              weight1*(1.+moist_star(i,j,k-1)) )  &
     &                            / weight3
              wet_to_dry_np1(i,j,k) = 1./dry_to_wet(i,j,k)
              temp = ( weight2  +                                       &
     &                 weight1 * (1. + moist(i,j,k-1)) )                &
     &               / weight3
! rho holds rho_moist, convert to rho_dry
              wet_to_dry_n(i,j,k) = 1. /temp
              rho(i,j,k) = rho(i,j,k) * wet_to_dry_n(i,j,k)
! compute t-averaged dry rho using last iteration's
! rho^{n+1} estimate
              temp = ( weight2 + weight1*(1.+moist_np1(i,j,k-1)) )      &
     &               / weight3
              rho_np1(i,j,k) = (1.-alpha_1)*rho(i,j,k)                  &
     &                       + alpha_1*rho_np1(i,j,k)/temp
            End Do
          End Do
        End If

        Do k=wet_model_levels+1, model_levels
          Do j = j0-1, j1+1
            Do i = 0, row_length+1
              rho_np1(i,j,k) = (1.-alpha_1)*rho(i,j,k)                  &
     &                       + alpha_1*rho_np1(i,j,k)
            End Do
          End Do
        End Do

      End If ! .NOT. L_new_tdisc


! DEPENDS ON: etadot_calc
      Call Etadot_Calc (r_theta_levels, r_rho_levels,                   &
     &                  eta_theta_levels, eta_rho_levels,               &
     &                  u1, v1, w1,                                     &
     &                  FV_sec_theta_latitude,                          &
     &                  row_length, rows, n_rows, model_levels,         &
     &                  delta_lambda, delta_phi,                        &
     &                  lambda_p, phi_p, lambda_u, phi_v,               &
     &                  wt_lambda_u, wt_phi_v,                          &
     &                  model_domain, first_constant_r_rho_level,       &
     &                  halo_i, halo_j, n_proc, proc_row_group,         &
     &                  at_extremity, g_row_length,                     &
     &                  off_x,off_y,halo_i_wind, halo_j_wind,           &
     &                  L_regular, rims_to_do, etadot1)

! ----------------------------------------------------------------------
! Section 2.   Calculate horizontal part of increment.
! ----------------------------------------------------------------------

      Do k = 1, model_levels

! calculate D(r)/D(eta) about rho levels

        Do j = j0-1, j1+1
          Do i = 0, row_length+1
            Dr_Deta_rho_levels(i,j,k) = (r_theta_levels(i,j,k) -        &
     &                                   r_theta_levels(i,j,k-1))/      &
     &                                  (eta_theta_levels(k) -          &
     &                                   eta_theta_levels(k-1))
          End Do
        End Do

! ----------------------------------------------------------------------
! Section 2.1  Calculate lambda direction term.
! ----------------------------------------------------------------------

! Calculate averaged term at u points
        If ( L_regular ) Then
        If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then
        Do j = j0, j1
          Do i = 0, row_length
            interp1(i,j) = ( rho(i+1,j,k) * Dr_Deta_rho_levels(i+1,j,k) &
     &                       + rho(i,j,k) * Dr_Deta_rho_levels(i,j,k) ) &
     &                       / ( r_rho_levels(i+1,j,k) +                &
     &                           r_rho_levels(i,j,k) )
          End Do
        End Do
        Else
          Do j = j0, j1
            Do i = 0, row_length
              interp1(i,j)=(rho_np1(i+1,j,k)*Dr_Deta_rho_levels(i+1,j,k)&
     &                    + rho_np1(i,j,k)* Dr_Deta_rho_levels(i,j,k) ) &
     &                    / ( r_rho_levels(i+1,j,k) +                   &
     &                        r_rho_levels(i,j,k) )
            End Do
          End Do
        End If

        else  ! variable resolution
        If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then
        Do j = j0, j1
          Do i = 0, row_length
            interp1(i,j) = ( lambda_pmiu(i+1) * rho(i,j,k) *            &
     &                           Dr_Deta_rho_levels(i,j,k) +            &
     &                       lambda_umip(i) * rho(i+1,j,k) *            &
     &                         Dr_Deta_rho_levels(i+1,j,k) ) /          &
     &              ( lambda_pmiu(i+1) * r_rho_levels(i,j,k) +          &
     &                lambda_umip(i) * r_rho_levels(i+1,j,k) )
          End Do
        End Do
        Else
        Do j = j0, j1
          Do i = 0, row_length
            interp1(i,j) = ( lambda_pmiu(i+1) * rho_np1(i,j,k) *        &
     &                           Dr_Deta_rho_levels(i,j,k) +            &
     &                       lambda_umip(i) * rho_np1(i+1,j,k) *        &
     &                         Dr_Deta_rho_levels(i+1,j,k) ) /          &
     &              ( lambda_pmiu(i+1) * r_rho_levels(i,j,k) +          &
     &                lambda_umip(i) * r_rho_levels(i+1,j,k) )
          End Do
        End Do
        End If

        endIf ! L_regular
! Calculate derivative at rho points

        If ( L_regular ) Then
        Do j = j0, j1
          Do i = 1, row_length
            increment(i,j,k) = ( u1(i,j,k) * interp1(i,j) -             &
     &                           u1(i-1,j,k) * interp1(i-1,j) ) *       &
     &                          recip_delta_lambda *                    &
     &                          FV_sec_theta_latitude(i,j)
          End Do
        End Do
        else  ! variable resolution
        Do j = j0, j1
          Do i = 1, row_length
            increment(i,j,k) = ( u1(i,j,k) * interp1(i,j) -             &
     &                           u1(i-1,j,k) * interp1(i-1,j) ) *       &
     &                                     recip_dlambda_u(i-1) *       &
     &                                 FV_sec_theta_latitude(i,j)
          End Do
        End Do
        endIf ! L_regular

! set values at northern and southern boundaries to zero.

      If (model_domain  /=  mt_bi_cyclic_lam) then
        If (at_extremity(PSouth)) Then
          Do i = 1, row_length
            increment(i,1,k) = 0.
          End Do
        End If
        If (at_extremity(PNorth)) Then
          Do i = 1, row_length
            increment(i,rows,k) = 0.
          End Do
        End If
      Endif

! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! ----------------------------------------------------------------------

! Calculate averaged term at v points

        If ( L_regular ) Then
      If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then
        Do j = j0-1, j1
          Do i = 1, row_length
            interp1(i,j) = ( rho(i,j,k) * Dr_Deta_rho_levels(i,j,k) +   &
     &                      rho(i,j+1,k) * Dr_Deta_rho_levels(i,j+1,k)) &
     &                       / ( r_rho_levels(i,j,k) +                  &
     &                           r_rho_levels(i,j+1,k) )
          End Do
        End Do
      Else
        Do j = j0-1, j1
          Do i = 1, row_length
            interp1(i,j)=(rho_np1(i,j,k)*Dr_Deta_rho_levels(i,j,k) +    &
     &                    rho_np1(i,j+1,k)*Dr_Deta_rho_levels(i,j+1,k)) &
     &                    / ( r_rho_levels(i,j,k) +                     &
     &                        r_rho_levels(i,j+1,k) )
          End Do
        End Do
      End If ! CycleNo == 1

        else  ! variable resolution
      If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then
        Do j = j0-1, j1
          Do i = 1, row_length
            interp1(i,j) =  ( phi_pmiv(j+1) * rho(i,j,k) *              &
     &                         Dr_Deta_rho_levels(i,j,k) +              &
     &                        phi_vmip(j) * rho(i,j+1,k) *              &
     &                         Dr_Deta_rho_levels(i,j+1,k) ) /          &
     &                ( phi_pmiv(j+1) * r_rho_levels(i,j,k) +           &
     &                  phi_vmip(j) * r_rho_levels(i,j+1,k) )
          End Do
        End Do
      Else
        Do j = j0-1, j1
          Do i = 1, row_length
            interp1(i,j) =  ( phi_pmiv(j+1) * rho_np1(i,j,k) *          &
     &                         Dr_Deta_rho_levels(i,j,k) +              &
     &                        phi_vmip(j) * rho_np1(i,j+1,k) *          &
     &                         Dr_Deta_rho_levels(i,j+1,k) ) /          &
     &                ( phi_pmiv(j+1) * r_rho_levels(i,j,k) +           &
     &                  phi_vmip(j) * r_rho_levels(i,j+1,k) )
          End Do
        End Do
      End If ! CycleNo == 1
        endIf ! L_regular
! Calculate derivative at rho points

        If ( L_regular ) Then
        Do j = j0, j1
          Do i = 1, row_length
            increment(i,j,k) = increment(i,j,k) +                       &
     &                         ( v1(i,j,k) * interp1(i,j) *             &
     &                           cos_v_latitude(i,j) -                  &
     &                           v1(i,j-1,k) * interp1(i,j-1) *         &
     &                           cos_v_latitude(i,j-1) ) *              &
     &                          recip_delta_phi *                       &
     &                          FV_sec_theta_latitude(i,j)
          End Do
        End Do
        else  ! variable resolution
        Do j = j0, j1
          Do i = 1, row_length
            increment(i,j,k) = increment(i,j,k) +                       &
     &                         ( v1(i,j,k) * interp1(i,j) *             &
     &                                cos_v_latitude(i,j) -             &
     &                           v1(i,j-1,k) * interp1(i,j-1) *         &
     &                                cos_v_latitude(i,j-1) ) *         &
     &                                    recip_dphi_v(i,j-1) *         &
     &                           FV_sec_theta_latitude(i,j)
          End Do
        End Do
        endIf ! L_regular

! calculate values at poles if a global model.

        If (model_domain  ==  mt_Global) Then

! save values at poles for polar calculation
          If(at_extremity(PSouth))then
            Do i=1,row_length
              l_s_poles(i,k) = v1(i,1,k) * interp1(i,1)
            End Do
           End If
          If(at_extremity(PNorth))then
            Do i=1,row_length
              l_n_poles(i,k) = v1(i,n_rows,k) * interp1(i,n_rows)
            End Do
          End If
        End If

      End Do  !  k = 1, model_levels

! calculate values at poles if a global model.

      If (model_domain == mt_Global) Then


        If (at_extremity(PSouth)) Then

          CALL global_2d_sums(l_s_poles, row_length, 1, 0, 0,           &
                              model_levels, sum_s, proc_row_group)

          Do k = 1, model_levels
            sum_s(k) = sum_s(k) * recip_delta_phi                       &
     &               * cos_v_latitude(1,1)                              &
     &               * FV_sec_theta_latitude(1,1) / g_row_length
            Do i = 1, row_length
              increment(i,1,k) = sum_s(k)
            End Do
          End Do
        End If

        If (at_extremity(PNorth)) Then
          CALL global_2d_sums(l_n_poles, row_length, 1, 0, 0,           &
                              model_levels, sum_n, proc_row_group)

          Do k = 1, model_levels
            sum_n(k) = - sum_n(k) * recip_delta_phi                     &
     &               * cos_v_latitude(1,n_rows)                         &
     &               * FV_sec_theta_latitude(1,rows) / g_row_length
            Do i = 1, row_length
              increment(i,rows,k) = sum_n(k)
            End Do
          End Do
        End If

      End If

! ----------------------------------------------------------------------
! Section 3.   Calculate vertical part of increment and add on to
!              horizontal part.
! ----------------------------------------------------------------------

      If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then

        k = 1
! calculate term inside vertical derivative
! interp1 holds value at lower level, interp2 value at upper level
      Do j = 1, rows
        Do i = 1, row_length
          interp1(i,j) = 0.
          interp2(i,j) =  ( rho(i,j,k) *                                &
     &                     (r_rho_levels(i,j,k+1) -                     &
     &                      r_theta_levels(i,j,k) ) +                   &
     &                      rho(i,j,k+1) *                              &
     &                     (r_theta_levels(i,j,k) -                     &
     &                      r_rho_levels(i,j,k) ) ) /                   &
     &                     (r_rho_levels(i,j,k+1) -                     &
     &                      r_rho_levels(i,j,k) ) * etadot1(i,j,k) *    &
     &                   (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)) /&
     &                   (eta_rho_levels(k+1) - eta_rho_levels(k))
          increment(i,j,k) = increment(i,j,k) + (interp2(i,j) -         &
     &                                           interp1(i,j)) /        &
     &                       (eta_theta_levels(k) -                     &
     &                        eta_theta_levels(k-1))

        End Do
      End Do

! calculate term inside vertical derivative
! interp1 holds value at lower level, interp2 value at upper level
      Do k = 2, model_levels - 1
        Do j = 1, rows
          Do i = 1, row_length
            interp1(i,j) = interp2(i,j)
            interp2(i,j) = ( rho(i,j,k) *                               &
     &                     (r_rho_levels(i,j,k+1) -                     &
     &                      r_theta_levels(i,j,k) ) +                   &
     &                      rho(i,j,k+1) *                              &
     &                     (r_theta_levels(i,j,k) -                     &
     &                      r_rho_levels(i,j,k) ) ) /                   &
     &                     (r_rho_levels(i,j,k+1) -                     &
     &                      r_rho_levels(i,j,k) ) * etadot1(i,j,k) *    &
     &                   (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))/ &
     &                   (eta_rho_levels(k+1) - eta_rho_levels(k))
            increment(i,j,k) = increment(i,j,k) + (interp2(i,j) -       &
     &                                             interp1(i,j)) /      &
     &                         (eta_theta_levels(k) -                   &
     &                          eta_theta_levels(k-1))

          End Do
        End Do
      End Do

      Else

      k = 1
! calculate term inside vertical derivative
! interp1 holds value at lower level, interp2 value at upper level
      Do j = 1, rows
        Do i = 1, row_length
          interp1(i,j) = 0.
          interp2(i,j) =  ( rho_np1(i,j,k) *                            &
     &                     (r_rho_levels(i,j,k+1) -                     &
     &                      r_theta_levels(i,j,k) ) +                   &
     &                      rho_np1(i,j,k+1) *                          &
     &                     (r_theta_levels(i,j,k) -                     &
     &                      r_rho_levels(i,j,k) ) ) /                   &
     &                     (r_rho_levels(i,j,k+1) -                     &
     &                      r_rho_levels(i,j,k) ) * etadot1(i,j,k) *    &
     &                   (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)) /&
     &                   (eta_rho_levels(k+1) - eta_rho_levels(k))
          increment(i,j,k) = increment(i,j,k) + (interp2(i,j) -         &
     &                                           interp1(i,j)) /        &
     &                       (eta_theta_levels(k) -                     &
     &                        eta_theta_levels(k-1))
        End Do
      End Do

! calculate term inside vertical derivative
! interp1 holds value at lower level, interp2 value at upper level
      Do k = 2, model_levels - 1
        Do j = 1, rows
          Do i = 1, row_length
            interp1(i,j) = interp2(i,j)
            interp2(i,j) =( rho_np1(i,j,k) *                            &
     &                     (r_rho_levels(i,j,k+1) -                     &
     &                      r_theta_levels(i,j,k) ) +                   &
     &                      rho_np1(i,j,k+1) *                          &
     &                     (r_theta_levels(i,j,k) -                     &
     &                      r_rho_levels(i,j,k) ) ) /                   &
     &                     (r_rho_levels(i,j,k+1) -                     &
     &                      r_rho_levels(i,j,k) ) * etadot1(i,j,k) *    &
     &                   (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)) /&
     &                   (eta_rho_levels(k+1) - eta_rho_levels(k))
            increment(i,j,k) = increment(i,j,k) + (interp2(i,j) -       &
     &                                             interp1(i,j)) /      &
     &                         (eta_theta_levels(k) -                   &
     &                          eta_theta_levels(k-1))

          End Do
        End Do
      End Do

      End If ! CycleNo == 1 .OR. .NOT. L_new_tdisc

      k = model_levels
! calculate term inside vertical derivative
! interp1 holds value at lower level, interp2 value at upper level
      Do j = 1, rows
        Do i = 1, row_length
          interp1(i,j) = interp2(i,j)
          interp2(i,j) = 0.
          increment(i,j,k) = increment(i,j,k) + (interp2(i,j) -         &
     &                                           interp1(i,j)) /        &
     &                       (eta_theta_levels(k) -                     &
     &                        eta_theta_levels(k-1))

        End Do
      End Do

! ----------------------------------------------------------------------
! Section 4.   Update dry density.
! ----------------------------------------------------------------------

! not updated on LAM boundaries
      If ( CycleNo == NumCycles ) Then
        Do k = 1, model_levels
          Do j = j_start, j_end
            Do i = i_start, i_end
              rho(i,j,k) = rho(i,j,k) - timestep * increment(i,j,k) /   &
     &                                  Dr_Deta_rho_levels(i,j,k)
            End Do
          End Do
        End Do
      Else
        Do k = 1, model_levels
          Do j = j_start, j_end
            Do i = i_start, i_end
               rho_np1(i,j,k) = rho(i,j,k) - timestep*increment(i,j,k)/ &
     &                                       Dr_Deta_rho_levels(i,j,k)
            End Do
          End Do
        End Do
      End If

! ----------------------------------------------------------------------
! Section 5.   Convert dry density to wet density.
! ----------------------------------------------------------------------
      If ( CycleNo == NumCycles ) Then
        levels = min( model_levels, wet_model_levels + 1 )
        Do k = 1, levels
          Do j = 1, rows
            Do i = 1, row_length
              rho(i,j,k) = rho(i,j,k) * dry_to_wet(i,j,k)
            End Do
          End Do
        End Do !  k = 1, levels
      Else
        levels = min( model_levels, wet_model_levels + 1 )
        Do k = 1, levels
          Do j = 1, rows
            Do i = 1, row_length
              rho_np1(i,j,k) = rho_np1(i,j,k) * dry_to_wet(i,j,k)
            End Do
          End Do
          Do j = j0-1, j1+1
            Do i = 0, row_length+1
              rho(i,j,k) = rho(i,j,k)/wet_to_dry_n(i,j,k)
            End Do
          End Do
        End Do
      End If ! If ( CycleNo == NumCycles )

! End of routine
      IF (lhook) CALL dr_hook('FLUX_RHO_MIX',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Flux_rho_mix

