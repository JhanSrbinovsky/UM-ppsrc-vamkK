! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Flux_rho_part
!

      Subroutine Flux_rho_part(                                         &
     &                         vap, vap_star, vap_np1,                  &
     &                         r_theta_levels, r_rho_levels,            &
     &                         rows, row_length,                        &
     &                         model_levels, wet_model_levels,          &
     &                         halo_i, halo_j, off_x, off_y,            &
     &                         alpha_1, j_begin, j_end,                 &
     &                         wet_to_dry_n, wet_to_dry_np1,            &
     &                         dry_to_wet,                              &
     &                         rho, rho_np1, CycleNo, L_new_tdisc,      &
     &                         L_mix_ratio, L_dry )

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
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
              ! model dimensions
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
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
     &, j_begin, j_end                                                  &
     &, CycleNo

      Logical                                                           &
     &  L_mix_ratio                                                     &
     &, L_dry                                                           &
     &, L_new_tdisc

      Real                                                              &
     &  alpha_1

      Real                                                              &
     &  vap (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
     &       wet_model_levels)                                          &
     &, wet_to_dry_n (1-off_x:row_length+off_x, 1-off_y:rows+off_y,     &
     &          model_levels)                                           &
     &, wet_to_dry_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,   &
     &          model_levels)                                           &
     &, vap_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &          wet_model_levels)                                       &
     &, vap_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &          wet_model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  rho (1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &       model_levels)                                              &
     &, rho_np1 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,          &
     &       model_levels)

      Real, Intent(Out) ::                                              &
     &  dry_to_wet (1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
     &              model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k      ! Loop indices

      Real                                                              &
     &  weight1                                                         &
     &, weight2                                                         &
     &, weight3                                                         &
     &, temp

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 1.   Convert density to dry density.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('FLUX_RHO_PART',zhook_in,zhook_handle)
      If (L_mix_ratio) then      

      If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then

        k = 1
        Do j = j_begin-1, j_end+1
          Do i = 0, row_length+1
            dry_to_wet(i,j,k) = 1. + vap_star(i,j,k)
            wet_to_dry_np1(i,j,k) = 1. /dry_to_wet(i,j,k)
! rho holds rho_moist, convert to rho_dry
            wet_to_dry_n(i,j,k) = 1. /(1. + vap(i,j,k))
            rho(i,j,k) = rho(i,j,k) * wet_to_dry_n(i,j,k)
          End Do
        End Do

      Do k = 2, wet_model_levels
        Do j = j_begin-1, j_end+1
          Do i = 0, row_length+1
            weight2 = r_rho_levels(i,j,k)                               &
     &                - r_theta_levels(i,j,k-1)
            weight1 = r_theta_levels(i,j,k)                             &
     &                - r_rho_levels(i,j,k)
            weight3 = r_theta_levels(i,j,k)                             &
     &                - r_theta_levels(i,j,k-1)
            dry_to_wet(i,j,k) = ( weight2 * (1. + vap_star(i,j,k)) +    &
     &                            weight1 * (1. + vap_star(i,j,k-1)) )  &
     &                          / weight3
            wet_to_dry_np1(i,j,k) = 1. /dry_to_wet(i,j,k)
            temp = ( weight2 * (1. + vap(i,j,k)) +                      &
     &               weight1 * (1. + vap(i,j,k-1)) ) / weight3
! rho holds rho_moist, convert to rho_dry
           wet_to_dry_n(i,j,k) = 1. / temp
           rho(i,j,k) = rho(i,j,k) * wet_to_dry_n(i,j,k)
          End Do
        End Do
      End Do !  k = 2, wet_model_levels

      k = wet_model_levels + 1
      if ( k <= model_levels ) then
        Do j = j_begin-1, j_end+1
          Do i = 0, row_length+1
            weight2 = r_rho_levels(i,j,k)                               &
     &                - r_theta_levels(i,j,k-1)
            weight1 = r_theta_levels(i,j,k)                             &
     &                - r_rho_levels(i,j,k)
            weight3 = r_theta_levels(i,j,k)                             &
     &                - r_theta_levels(i,j,k-1)
            dry_to_wet(i,j,k) = ( weight2  +                            &
     &                            weight1 * (1. + vap_star(i,j,k-1)) )  &
     &                          / weight3
            wet_to_dry_np1(i,j,k) = 1. /dry_to_wet(i,j,k)
            temp = ( weight2  +                                         &
     &               weight1 * (1. + vap(i,j,k-1)) )                    &
     &             / weight3
! rho holds rho_moist, convert to rho_dry
            wet_to_dry_n(i,j,k) = 1. / temp
            rho(i,j,k) = rho(i,j,k) * wet_to_dry_n(i,j,k)
          End Do
        End Do

      endif !  k <= model_levels

      Else  !  CycleNo > 1  .AND. L_new_tdisc

        k = 1
        Do j = j_begin-1, j_end+1
          Do i = 0, row_length+1
            dry_to_wet(i,j,k)    =  1.+ vap_star(i,j,k)
            wet_to_dry_np1(i,j,k)=  1./ dry_to_wet(i,j,k)
! rho holds rho_moist, convert to rho_dry
            wet_to_dry_n(i,j,k)  =  1./(1. + vap(i,j,k))
            rho(i,j,k) = rho(i,j,k) * wet_to_dry_n(i,j,k)
! compute t-averaged dry rho using last iteration's
! rho^{n+1} estimate
            rho_np1(i,j,k) = (1.- alpha_1) * rho(i,j,k) +               &
     &                     alpha_1*rho_np1(i,j,k)/(1. + vap_np1(i,j,k))
          End Do
        End Do

        Do k = 2, wet_model_levels
          Do j = j_begin-1, j_end+1
            Do i = 0, row_length+1
              weight2 = r_rho_levels(i,j,k)                             &
     &                - r_theta_levels(i,j,k-1)
              weight1 = r_theta_levels(i,j,k)                           &
     &                - r_rho_levels(i,j,k)
              weight3 = r_theta_levels(i,j,k)                           &
     &                - r_theta_levels(i,j,k-1)
              dry_to_wet(i,j,k) = ( weight2*(1.+ vap_star(i,j,k)) +     &
     &                              weight1*(1.+ vap_star(i,j,k-1)) )   &
     &                            / weight3
              wet_to_dry_np1(i,j,k) = 1. / dry_to_wet(i,j,k)
              temp = ( weight2 * (1. + vap(i,j,k)) +                    &
     &                 weight1 * (1. + vap(i,j,k-1)) ) / weight3
! rho holds rho_moist, convert to rho_dry
              wet_to_dry_n(i,j,k) = 1. / temp
              rho(i,j,k) = rho(i,j,k) * wet_to_dry_n(i,j,k)
! compute t-averaged dry rho using last iteration's
! rho^{n+1} estimate
              temp = ( weight2 * (1. + vap_np1(i,j,k)) +                &
     &                 weight1 * (1. + vap_np1(i,j,k-1)) )              &
     &               / weight3
              rho_np1(i,j,k) = (1.-alpha_1) * rho(i,j,k) +              &
     &                         alpha_1 * rho_np1(i,j,k) / temp
            End Do
          End Do
        End Do !  k = 2, wet_model_levels

        k = wet_model_levels+1
        If ( k <= model_levels ) Then
          Do j = j_begin-1, j_end+1
            Do i = 0, row_length+1

              weight2 = r_rho_levels(i,j,k)                             &
     &                - r_theta_levels(i,j,k-1)
              weight1 = r_theta_levels(i,j,k)                           &
     &                - r_rho_levels(i,j,k)
              weight3 = r_theta_levels(i,j,k)                           &
     &                - r_theta_levels(i,j,k-1)
              dry_to_wet(i,j,k) = ( weight2  +                          &
     &                              weight1*(1. + vap_star(i,j,k-1)) )  &
     &                            / weight3
              wet_to_dry_np1(i,j,k) = 1. / dry_to_wet(i,j,k)
              temp = ( weight2  +                                       &
     &                 weight1 * (1. + vap(i,j,k-1)) )                  &
     &               / weight3
! rho holds rho_moist, convert to rho_dry
              wet_to_dry_n(i,j,k) = 1. / temp
              rho(i,j,k) = rho(i,j,k) * wet_to_dry_n(i,j,k)
! compute t-averaged dry rho using last iteration's
! rho^{n+1} estimate
              temp = ( weight2 + weight1 * (1. + vap_np1(i,j,k-1)) )    &
     &               / weight3
              rho_np1(i,j,k) = (1.-alpha_1) * rho(i,j,k) +              &
     &                         alpha_1 * rho_np1(i,j,k) / temp
            End Do
          End Do
        End If ! k <= model_levels

        Do k = wet_model_levels+1, model_levels
          Do j = j_begin-1, j_end+1
            Do i = 0, row_length+1
              rho_np1(i,j,k) = (1.-alpha_1) * rho(i,j,k) +              &
     &                              alpha_1 * rho_np1(i,j,k)
            End Do
          End Do
        End Do

      End If ! .NOT. L_new_tdisc
        
      Else If (L_dry) then  !  dry run      

        Do k = 1, model_levels
          Do j = j_begin-1, j_end+1
            Do i = 0, row_length+1
              dry_to_wet(i,j,k) =  1.0
              wet_to_dry_np1(i,j,k) =  dry_to_wet(i,j,k)
              wet_to_dry_n(i,j,k)=  1.0
            End Do
          End Do
        End Do  !  k = 1, model_levels

        If ( CycleNo > 1 ) then

! compute t-averaged dry rho using last iteration's
! rho^{n+1} estimate
          Do k = 1, model_levels
            Do j = j_begin-1, j_end+1
              Do i = 0, row_length+1
                rho_np1(i,j,k) = (1.0 - alpha_1) * rho(i,j,k) +         &
     &                               alpha_1 * rho_np1(i,j,k)
              End Do
            End Do
          End Do  !  k = 1, model_levels

        End If ! CycleNo > 1

      Else ! L_mix_ratio is false - specific quantities      

      If ( CycleNo == 1 .OR. .NOT. L_new_tdisc ) Then

        k = 1
        Do j = j_begin-1, j_end+1
          Do i = 0, row_length+1
            dry_to_wet(i,j,k)    =  1. /(1. - vap_star(i,j,k))
            wet_to_dry_np1(i,j,k)=  1. - vap_star(i,j,k)
            wet_to_dry_n(i,j,k)  =  1. - vap(i,j,k)
            rho(i,j,k) = rho(i,j,k) * wet_to_dry_n(i,j,k)
          End Do
        End Do

      Do k = 2, wet_model_levels
        Do j = j_begin-1, j_end+1
          Do i = 0, row_length+1
            weight2 = r_rho_levels(i,j,k)                               &
     &                - r_theta_levels(i,j,k-1)
            weight1 = r_theta_levels(i,j,k)                             &
     &                - r_rho_levels(i,j,k)
            weight3 = r_theta_levels(i,j,k)                             &
     &                - r_theta_levels(i,j,k-1)
            dry_to_wet(i,j,k) =                                         &
     &                ( weight2 * (1. - vap_star(i,j,k)  ) +            &
     &                  weight1 * (1. - vap_star(i,j,k-1)) )            &
     &                          / weight3
            wet_to_dry_np1(i,j,k)=  dry_to_wet(i,j,k)
            dry_to_wet(i,j,k) = 1./ dry_to_wet(i,j,k)
            temp = ( weight2 * (1. - vap(i,j,k)  ) +                    &
     &               weight1 * (1. - vap(i,j,k-1)) )                    &
     &             / weight3
            rho(i,j,k) = rho(i,j,k) * temp
            wet_to_dry_n(i,j,k)=  temp
          End Do
        End Do
      End Do

      k = wet_model_levels + 1
      if ( k <= model_levels ) then

        Do j = j_begin-1, j_end+1
          Do i = 0, row_length+1
            weight2 = r_rho_levels(i,j,k)                               &
     &                - r_theta_levels(i,j,k-1)
            weight1 = r_theta_levels(i,j,k)                             &
     &                - r_rho_levels(i,j,k)
            weight3 = r_theta_levels(i,j,k)                             &
     &                - r_theta_levels(i,j,k-1)

            dry_to_wet(i,j,k) = ( weight2  +                            &
     &                            weight1 * (1. - vap_star(i,j,k-1)) )  &
     &                          / weight3
            wet_to_dry_np1(i,j,k)=  dry_to_wet(i,j,k)
            dry_to_wet(i,j,k) = 1./ dry_to_wet(i,j,k)
            temp = ( weight2  +                                         &
     &               weight1 * (1. - vap(i,j,k-1)) )                    &
     &             / weight3
            rho(i,j,k) = rho(i,j,k) * temp
            wet_to_dry_n(i,j,k)=  temp
          End Do
        End Do

      endif !  k <= model_levels

      Else  !  CycleNo > 1  .AND. L_new_tdisc

        k = 1
        Do j = j_begin-1, j_end+1
          Do i = 0, row_length+1
            dry_to_wet(i,j,k)    =  1. /(1. - vap_star(i,j,k))
            wet_to_dry_np1(i,j,k)=  1. - vap_star(i,j,k)
            wet_to_dry_n(i,j,k)  =  1. - vap(i,j,k)
! compute dry density at timelevel n
            rho(i,j,k)= rho(i,j,k)*wet_to_dry_n(i,j,k)
! compute t-averaged dry rho using last iteration's
! rho^{n+1} estimate
            rho_np1(i,j,k)= (1.-alpha_1)*rho(i,j,k)                     &
     &                    + alpha_1*rho_np1(i,j,k)*(1.-vap_np1(i,j,k))
          End Do
        End Do

        Do k = 2, wet_model_levels
          Do j = j_begin-1, j_end+1
            Do i = 0, row_length+1
              weight2 = r_rho_levels(i,j,k)                             &
     &                - r_theta_levels(i,j,k-1)
              weight1 = r_theta_levels(i,j,k)                           &
     &                - r_rho_levels(i,j,k)
              weight3 = r_theta_levels(i,j,k)                           &
     &                - r_theta_levels(i,j,k-1)
              dry_to_wet(i,j,k) =                                       &
     &                ( weight2 * (1. - vap_star(i,j,k)  ) +            &
     &                  weight1 * (1. - vap_star(i,j,k-1)) )            &
     &                          / weight3
              wet_to_dry_np1(i,j,k)=  dry_to_wet(i,j,k)
              dry_to_wet(i,j,k) = 1./ dry_to_wet(i,j,k)
              temp = ( weight2 * (1. - vap(i,j,k)  ) +                  &
     &                 weight1 * (1. - vap(i,j,k-1)) )                  &
     &               / weight3
! compute dry density at timelevel n
              rho(i,j,k) = rho(i,j,k)*temp
! compute t-averaged dry rho using last iteration's
! rho^{n+1} estimate
              rho_np1(i,j,k) = (1.-alpha_1)*rho(i,j,k)                  &
     &               + alpha_1*rho_np1(i,j,k)                           &
     &               * ( weight2 * (1. - vap_np1(i,j,k)  ) +            &
     &                  weight1 * (1. - vap_np1(i,j,k-1)) )/weight3
              wet_to_dry_n(i,j,k)  =  temp
            End Do
          End Do
        End Do  !  k = 2, wet_model_levels

        k = wet_model_levels+1
        If ( k <= model_levels ) Then
          Do j = j_begin-1, j_end+1
            Do i = 0, row_length+1

              weight2 = r_rho_levels(i,j,k)                             &
     &                - r_theta_levels(i,j,k-1)
              weight1 = r_theta_levels(i,j,k)                           &
     &                - r_rho_levels(i,j,k)
              weight3 = r_theta_levels(i,j,k)                           &
     &                - r_theta_levels(i,j,k-1)

              dry_to_wet(i,j,k) = ( weight2  +                          &
     &                            weight1 * (1. - vap_star(i,j,k-1)) )  &
     &                          / weight3
              wet_to_dry_np1(i,j,k)=  dry_to_wet(i,j,k)
              dry_to_wet(i,j,k) = 1./ dry_to_wet(i,j,k)
              temp = ( weight2  +                                       &
     &               weight1 * (1. - vap(i,j,k-1)) )                    &
     &             / weight3
! compute dry density at timelevel n
              rho(i,j,k) = rho(i,j,k)*temp
! compute t-averaged dry rho using last iteration's
! rho^{n+1} estimate
              rho_np1(i,j,k) = (1.-alpha_1)*rho(i,j,k)                  &
     &               + alpha_1*rho_np1(i,j,k)                           &
     &               * (weight2 +weight1*(1.-vap_np1(i,j,k-1)))/weight3
              wet_to_dry_n(i,j,k)  =  temp
            End Do
          End Do
        End If !  k <= model_levels

        Do k = wet_model_levels+1, model_levels
          Do j = j_begin-1, j_end+1
            Do i = 0, row_length+1
              rho_np1(i,j,k) = (1.-alpha_1) * rho(i,j,k) +              &
     &                              alpha_1 * rho_np1(i,j,k)
            End Do
          End Do
        End Do

      End If ! .NOT. L_new_tdisc

      EndIf ! L_mix_ratio     

      IF (lhook) CALL dr_hook('FLUX_RHO_PART',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Flux_rho_part

