! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_calc_exner
!
      Subroutine IDL_calc_exner(                                        &
     &                      row_length, rows, model_levels              &
     &,                     wet_levels, me, halo_i, halo_j              &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     theta, q, exner_rho_levels                  &
     &,                     exner_theta_levels                          &
     &,                     height_domain                               &
! Profile settings
     &,                     theta_ref, q_ref, exner_ref                 &
     &,                     p_surface, theta_surface                    &
!  Options
     &,                     L_code_test)

! Purpose: Calculates hydrostatic exner pressure from
!          virtual potential temperature.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
! Code Description:
!   Language: FORTRAN 77 + extensions
!   This code is written to UMDP3 programming standards.
!
      USE atmos_constants_mod, ONLY:                                    &
          r, cp, kappa, repsilon, p_zero 

      USE earth_constants_mod, ONLY: g, earth_radius

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

      Integer, Intent(In) ::                                            &
     &  row_length                                                      &
                           ! Number of points on a row
     &, rows                                                            &
                           ! Number of rows in a theta field
     &, model_levels                                                    &
                           ! Number of model levels
     &, wet_levels                                                      &
                           ! Number of model wet levels
     &, halo_i                                                          &
                           ! Size of halo in i direction.
     &, halo_j             ! Size of halo in j direction.

      Real, Intent(In) ::                                               &
     &  theta_surface                                                   &
                           ! Surface temperature
     &, p_surface                                                       &
                           ! Surface pressure
     &, height_domain                                                   &
                           ! Height of top of domain
     &, theta_ref(model_levels)                                         &
                                    ! Reference theta profile
     &, q_ref(model_levels)         ! Reference q profile

      Logical, Intent(In) ::                                            &
     &  L_code_test        ! User switch

      Integer, Intent(In) ::                                            &
     &  me                 ! My processor number


      ! Vertical co-ordinate information
      Real, Intent(InOut) ::                                            &
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, model_levels)                &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

      ! Output Arrays from this routine
      Real, Intent(InOut) ::                                            &
     &  theta(1-halo_i:row_length+halo_i                                &
     &,       1-halo_j:rows+halo_j, model_levels)                       &
     &, q(1-halo_i:row_length+halo_i                                    &
     &,       1-halo_j:rows+halo_j, wet_levels)                         &
     &, exner_rho_levels(1-halo_i:row_length+halo_i                     &
     &,                  1-halo_j:rows+halo_j, model_levels+1)          &
     &, exner_theta_levels(1-halo_i:row_length+halo_i                   &
     &,                    1-halo_j:rows+halo_j, model_levels)          &
     &, exner_ref(model_levels + 1) ! Reference exner profile

! Local variables

      Integer                                                           &
     &  i, j, k, k2          ! Loop indices

      Real                                                              &
     &  weight                                                          &
     &, re_epsilon                                                      &
     &, exner_surface                                                   &
     &, theta_q_surface                                                 &
     &, theta_ref_at_rho1                                               &
     &, z_at_theta, z_at_theta_p1                                       &
     &, z_at_rho,   z_at_rho_p1                                         &
     &, z_at_orog

      ! Work arrays
      Real                                                              &
     &  theta_temp(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)    &
     &, theta_orog(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)    &
     &, theta_rholev1(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j) &
     &, exner_orog(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)    &
     &, work1(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, thetav_ref(model_levels)

      ! Error reporting
      Character (Len=*),  Parameter :: RoutineName='idl_calc_exner'
      Character (Len=256)           :: Cmessage
      Integer                       :: ErrorStatus

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------
!
! Section 1  set profile from namelist data
!
!-------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_CALC_EXNER',zhook_in,zhook_handle)

      ! Calculate exner pressure at surface
      re_epsilon =  1.0 / repsilon  - 1.0
      exner_surface = (p_surface/p_zero) ** kappa


      ! Set up reference thetav profile (i.e. no orography)
      Do k = 1, model_levels
        thetav_ref(k) = theta_ref(k) * (1.0 + re_epsilon * q_ref(k))
      End Do

!-----------------------------------------------------------------
! Generate pressure. (From hydrostatic relation for
!                     Exner pressure)
!-------------------------------------------------------------------
!  First set up reference exner profile (i.e. no orography) on rho levs

      ! Calculate thetav_ref on first rho level
      theta_q_surface = theta_surface *                                 &
     &                    (1.0 + re_epsilon * q_ref(1))
      weight = (eta_rho_levels(1)*height_domain)                        &
     &                 /(eta_theta_levels(1)*height_domain)
      theta_ref_at_rho1 = theta_q_surface + weight*                     &
     &                    (thetav_ref(1) - theta_q_surface)

      ! Calculate exner_ref on first rho levels using average theta
      exner_ref(1) = exner_surface - g *                                &
     &               eta_rho_levels(1) * height_domain /                &
     &               (Cp * 0.5*(theta_ref_at_rho1+theta_q_surface) )

      Do k = 2,model_levels
        exner_ref(k)= exner_ref(k-1) - g * height_domain *              &
     &                 (eta_rho_levels(k) - eta_rho_levels(k-1))        &
     &                                     / ( Cp * thetav_ref(k-1))
      End Do

      k = model_levels + 1
      exner_ref(k) = exner_ref(k-1) - 2.0 * g * height_domain *         &
     &              (eta_theta_levels(k-1) - eta_rho_levels(k-1))       &
     &                                   / ( Cp * thetav_ref(k-1) )

! Find exner pressure at r_theta_levels(i,j,0), i.e. on the orography
! and put in exner_orog(i,j)

      Do k = 1, model_levels
        z_at_rho    = eta_rho_levels(k)*height_domain
        IF (k == model_levels) THEN
          ! Next level above is the same height difference from level below.
          z_at_rho_p1 = (2*eta_rho_levels(k)-eta_rho_levels(k-1))*height_domain
        ELSE
          z_at_rho_p1 = eta_rho_levels(k+1)*height_domain
        END IF
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            z_at_orog  = r_theta_levels(i,j,0) - Earth_radius
            If (k == 1 .and. z_at_orog < z_at_rho) Then
              ! Set exner_orog if below rho level 1
              weight = z_at_orog/z_at_rho
              exner_orog(i,j) = exner_surface + weight*                 &
     &                       (exner_ref(1) - exner_surface)
            Else If (z_at_orog >= z_at_rho .and.                        &
     &               z_at_orog < z_at_rho_p1) Then
              ! Set exner_orog if above rho level 1
              weight = (z_at_orog - z_at_rho)                           &
     &                 /(z_at_rho_p1 - z_at_rho)
              exner_orog(i,j) = exner_ref(k) + weight*                  &
     &                        (exner_ref(k+1) - exner_ref(k))
            End If
          End Do
        End Do
      End Do

! Find theta at r_theta_levels(i,j,0), i.e. on the orography
! and put in theta_orog(i,j).

      Do k = 1, model_levels - 1
        z_at_theta  = eta_theta_levels(k)*height_domain
        z_at_theta_p1 = eta_theta_levels(k+1)*height_domain
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            z_at_orog  = r_theta_levels(i,j,0) - Earth_radius
            If (k == 1 .and. z_at_orog < z_at_theta) Then
              ! Set theta_orog if below theta level 1
              weight = z_at_orog/z_at_theta
              theta_orog(i,j) = theta_q_surface + weight*               &
     &                       (thetav_ref(1) - theta_q_surface)
            Else If (z_at_orog >= z_at_theta .and.                      &
     &               z_at_orog < z_at_theta_p1) Then
              weight = (z_at_orog - z_at_theta)                         &
     &                 /(z_at_theta_p1 - z_at_theta)
              theta_orog(i,j) = thetav_ref(k) + weight*                 &
     &                        (thetav_ref(k+1) - thetav_ref(k))
            End If
          End Do
        End Do
      End Do

! Find theta at r_rho_levels(i,j,1), i.e. first rho levels
! and put in theta_rholev1(i,j)

      Do k = 1, model_levels - 1
        z_at_theta  = eta_theta_levels(k)*height_domain
        z_at_theta_p1 = eta_theta_levels(k+1)*height_domain
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            z_at_rho = r_rho_levels(i,j,1) - Earth_radius
            If (k == 1 .and. z_at_rho < z_at_theta) Then
              ! Set theta_rholev1 if below rho level 1
              weight = z_at_rho/z_at_theta
              theta_rholev1(i,j) = theta_q_surface + weight*            &
     &                       (thetav_ref(1) - theta_q_surface)
            Else If (z_at_rho >= z_at_theta .and.                       &
     &               z_at_rho < z_at_theta_p1) Then
              weight = (z_at_rho - z_at_theta)                          &
     &                 /(z_at_theta_p1 - z_at_theta)
              theta_rholev1(i,j) = thetav_ref(k) + weight*              &
     &                        (thetav_ref(k+1) - thetav_ref(k))
            End If
          End Do
        End Do
      End Do

      ! Put theta midway between orography and first rho level
      ! into theta_temp(i,j)
      ! theta_orog holds theta on orography, r_theta_levels(i,j,0)
      ! theta_rholev1 holds theta on 1st rho level, r_rho_levels(i,j,1)

      Do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          theta_temp(i,j) = 0.5 * (theta_orog(i,j)+theta_rholev1(i,j))
        End Do
      End Do


      Do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          exner_rho_levels(i,j,1) = exner_orog(i,j) - g *               &
     &                     (r_rho_levels(i,j,1) - r_theta_levels(i,j,0))&
     &                               / ( Cp * theta_temp(i,j) )
        end do
      end do

      do k = 2,model_levels
        do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            exner_rho_levels(i,j,k)= exner_rho_levels(i,j,k-1)-         &
     &             g * (r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1))    &
     &     / ( Cp * theta(i,j,k-1) *(1.0 + re_epsilon * q(i,j,k-1)) )
          end do
        end do
      end do

      k = model_levels + 1
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          exner_rho_levels(i,j,k) = exner_rho_levels(i,j,k-1) -         &
     &     2.0 * g * (r_theta_levels(i,j,k-1) - r_rho_levels(i,j,k-1))  &
     &     / ( Cp * theta(i,j,k-1) * (1.0 + re_epsilon * q(i,j,k-1)) )
        end do
      end do

      Do k = 1, model_levels - 1
        do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            weight = (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))/   &
     &               (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
            exner_theta_levels(i,j,k) =                                 &
     &                         weight * exner_rho_levels(i,j,k) +       &
     &                 (1.0 - weight) * exner_rho_levels(i,j,k+1)
          End Do
        End Do
      End Do

      k =  model_levels
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          exner_theta_levels(i,j,k) =  0.5 * exner_rho_levels(i,j,k) +  &
     &                                 0.5 * exner_rho_levels(i,j,k+1)
        End Do
      End Do

      IF (lhook) CALL dr_hook('IDL_CALC_EXNER',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_calc_exner

      !  End subroutine IDL_calc_exner

