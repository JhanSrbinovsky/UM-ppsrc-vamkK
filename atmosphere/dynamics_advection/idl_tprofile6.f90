! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine IDL_tprofile6

      Subroutine IDL_tprofile6(                                         &
     &                      row_length, rows, model_levels              &
     &,                     me, off_x, off_y, halo_i, halo_j            &
     &,                     cos_theta_latitude                          &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     theta, exner_rho_levels, exner_theta_levels &
!  Profiles for fixed lbcs and sponge zones
     &,                     theta_ref, exner_ref                        &
!  Grid information
     &,                     height_domain, big_layers                   &
! Profile settings
     &,                     p_surface, theta_surface                    &
! Dynamical core settings
     &,                     SuHe_pole_equ_deltaT, SuHe_static_stab      &
     &,                     L_SH_Williamson                             &
!  Options
     &,                     L_code_test)

! Purpose:
!          Sets up initial data for idealised problems.
!          Zonally symmetric initial temperature profile
!          for global dynamical core
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
!
      USE atmos_constants_mod, ONLY:                                    &
          r, cp, kappa, p_zero, recip_kappa 

      USE earth_constants_mod, ONLY: g, earth_radius
       
      USE conversions_mod, ONLY: pi

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, model_levels                                                    &
                         ! number of model levels
     &, big_layers                                                      &
                          ! number of isothermal layers
     &, halo_i                                                          &
                             ! Size of halo in i direction.
     &, halo_j                                                          &
                             ! Size of halo in j direction.
     &, off_x                                                           &
                             ! Size of small halo in i
     &, off_y                ! Size of small halo in j.

      Real                                                              &
     &  theta_surface                                                   &
     &, p_surface                                                       &
     &, height_domain

      Logical                                                           &
     & L_code_test    ! user switch

      Integer                                                           &
     &  me         ! My processor number


      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i                       &
     &,                1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i                         &
     &,              1-halo_j:rows+halo_j, model_levels)                &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)                                    &
     &, cos_theta_latitude(1-off_x:row_length+off_x                     &
     &,                     1-off_y:rows+off_y)

!Output Arrays from this routine
      Real                                                              &
     &  theta(1-halo_i:row_length+halo_i                                &
     &,       1-halo_j:rows+halo_j, model_levels)                       &
     &, exner_rho_levels(1-halo_i:row_length+halo_i                     &
     &,                  1-halo_j:rows+halo_j, model_levels+1)          &
     &, exner_theta_levels(1-halo_i:row_length+halo_i                   &
     &,                    1-halo_j:rows+halo_j, model_levels)          &
     &, theta_ref(model_levels)                                         &
                                 !theta profile for use in sponge & lbcs
     &, exner_ref(model_levels + 1)  ! Exner profile for use in lbcs

! Suarez Held variables
      Real                                                              &
     &  SuHe_pole_equ_deltaT                                            &
     &, SuHe_static_stab

      Logical                                                           &
     &  L_SH_williamson

! local variables
      Integer                                                           &
     &  i, j, k

      Real                                                              &
        exner_surface                                                   &
      , temp1                                                           &
      , temp2                                                           &
      , weight                                                          &
      , sigma_ref(model_levels)                                         &
      , sigma_to_kappa(model_levels)

! Williamson variables1. 1. parameters.

      Real                                                              &
     &  p_d                                                             &
     &, p_pl                                                            &
     &, lapse_rate_d                                                    &
     &, lapse_rate_i                                                    &
     &, delta_phi_0                                                     &
     &, A_will                                                          &
     &, phi_0                                                           &
     &, T_0

! 2. derived variables
      Real                                                              &
     &  p_i                                                             &
     &, p_eq                                                            &
     &, power_d                                                         &
     &, power_i                                                         &
     &, latitude                                                        &
     &, p_lim                                                           &
     &, minusgoverRTref

      Real                                                              &
     &  p_theta_levels(1-off_x:row_length+off_x, 1-off_y:rows+off_y     &
     &,    model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------
! Section 1  Zonally symmetric initial data based on dynamical core
!              set up (Held-Suarez)
!-------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_TPROFILE6',zhook_in,zhook_handle)
      exner_surface = (p_surface/p_zero) ** kappa
      minusgoverRTref = - g / ( R * 250.0)

      if(me  ==  0)then
        print*,'** Zonally symmetric profile for dynamical core ***'
      Endif   !(me  ==  0)

      do k= 1, model_levels
        sigma_ref(k) = exp( minusgoverRTref *                           &
     &                       (r_theta_levels(1,1,k) -  Earth_radius ) )
        sigma_to_kappa(k) =  sigma_ref(k) ** kappa
      end do

      Do k = 1, model_levels
        temp1 = (315. -                                                 &
     &                      SuHe_static_stab * log(sigma_ref(k)))       &
     &                    * sigma_to_kappa(k)
        do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
            temp2 = (315. - SuHe_pole_equ_deltaT *                      &
     &                   (1.0 - cos_theta_latitude(i,j)                 &
     &                        * cos_theta_latitude(i,j))                &
     &                    - SuHe_static_stab * log(sigma_ref(k))        &
     &                    * cos_theta_latitude(i,j)                     &
     &                    * cos_theta_latitude(i,j))                    &
     &                    * sigma_to_kappa(k)

            theta(i,j,k) = max(200.0, temp2)
! NB  theta currently holds temperature
          End Do
        End Do
        theta_ref(k) = max(200.0, temp1)
      End Do
! ----------------------------------------------------------------------
! Section 1.6.1  Calculate pressure profile (overwriting dump)
! NB  theta currently holds temperature
! ----------------------------------------------------------------------
      do j = 1-off_y, rows+off_y
        Do i = 1-off_x, row_length+off_x
          exner_rho_levels(i,j,1) = exp( g *                            &
     &                ( r_theta_levels(i,j,0) - r_rho_levels(i,j,1) )   &
     &                               / ( Cp * theta(i,j,1) ) )
        end do
      end do

      Do k = 2, model_levels
        do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
            exner_rho_levels(i,j,k) = exner_rho_levels(i,j,k-1) *       &
     &          exp ( g * (r_rho_levels(i,j,k-1) - r_rho_levels(i,j,k)) &
     &                                     / ( Cp * theta(i,j,k-1)))
          End Do
        End Do
      End Do

      k = model_levels
      do j = 1-off_y, rows+off_y
        Do i = 1-off_x, row_length+off_x
          exner_rho_levels(i,j,k+1) = exner_rho_levels(i,j,k) *         &
     &      exp (2.0 * g * (r_rho_levels(i,j,k) - r_theta_levels(i,j,k))&
     &                                     / ( Cp * theta(i,j,k)))
        End Do
      End Do

      Do k = 1, model_levels - 1
        do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
            weight = (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))/   &
     &                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
            exner_theta_levels(i,j,k) =   weight *                      &
     &                                   exner_rho_levels(i,j,k) +      &
     &                  (1.0 - weight) * exner_rho_levels(i,j,k+1)
          End Do
        End Do
      End Do
      k = model_levels
      do j = 1-off_y, rows+off_y
        Do i = 1-off_x, row_length+off_x
          exner_theta_levels(i,j,k) =   0.5 *                           &
     &            (exner_rho_levels(i,j,k) + exner_rho_levels(i,j,k+1))
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 1.6.2  More realistical stratosphere
! ----------------------------------------------------------------------
      If (L_SH_Williamson .and. height_domain  <   45000.0)then
! Williamson profile for more realistic stratosphere
!  only realistic for altitudes up to 2hPa (about 45km)
! calculate p at theta levels
        do k = 1, model_levels
          do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
              p_theta_levels(i,j,k) = p_zero *                          &
     &                        exner_theta_levels(i,j,k) ** recip_Kappa
            end do
          end do
        end do

        p_d = 10000.0                 ! 100 hpa
        p_pl = 200.0                  ! limited to 2hpa
        delta_phi_0 = pi / 12.0       ! 15 degrees converted to radians
        A_will = 2.65 / delta_phi_0
        phi_0 = 60.0 * pi / 180.      ! 60 degrees converted to radians
        T_0 = 200.0                   ! temperature at 100 hpa in
                                      !Held-Suarez
        lapse_rate_d = 0.002          ! 2 K /km
        lapse_rate_i = -0.003345      ! 3.345 K /km

        p_eq = p_d
        power_d = R * lapse_rate_d / g
        power_i = R * lapse_rate_i / g

        Do k = 1, model_levels
          do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
!  p_theta_levels  constant along longitude
              p_lim = p_theta_levels(i,j,k)
! restrict Williamson temperature to not change with height above 2hpa
! as values get very unphysical.
              If (p_lim  <   200.0) p_lim =200.0

              If (p_lim  <=  p_d) Then
                temp1 = T_0 * (p_lim / p_d) ** power_d

                latitude = acos(cos_theta_latitude(i,j))

                p_i = p_eq - (p_eq - p_pl) * 0.5 *                      &
     &                       (1+tanh(A_will*(abs(latitude)-phi_0)))

                If (p_lim  <=  p_i) Then
                  temp2 = T_0 *                                         &
     &                    ((p_lim / p_i) ** power_i - 1.0)
                  temp1 = temp1 + temp2
                End If
!  temp1 holds temperature. Converted to potential temperature later
                theta(i,j,k) =  temp1
              EndIf !p_lim  <=  p_d
            End Do  ! i = 1, row_length
          End Do  !j = 1, rows
        End Do  ! k = 1, model_levels

! ----------------------------------------------------------------------
! Section 1.6.3  Re Calculate pressure profile (overwriting dump)
! NB  theta currently holds temperature
! ----------------------------------------------------------------------
        do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
            exner_rho_levels(i,j,1) = exp( g *                          &
     &                ( r_theta_levels(i,j,0) - r_rho_levels(i,j,1) )   &
     &                               / ( Cp * theta(i,j,1) ) )
          end do
        end do

        Do k = 2, model_levels
          do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
              exner_rho_levels(i,j,k) = exner_rho_levels(i,j,k-1) *     &
     &          exp ( g * (r_rho_levels(i,j,k-1) - r_rho_levels(i,j,k)) &
     &                                     / ( Cp * theta(i,j,k-1)))
            End Do
          End Do
        End Do

        k = model_levels
        do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
            exner_rho_levels(i,j,k+1) = exner_rho_levels(i,j,k) *       &
     &      exp (2.0 * g * (r_rho_levels(i,j,k) - r_theta_levels(i,j,k))&
     &                                     / ( Cp * theta(i,j,k)))
          End Do
        End Do

        Do k = 1, model_levels - 1
          do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
              weight = (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))/ &
     &                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
              exner_theta_levels(i,j,k) =   weight *                    &
     &                                   exner_rho_levels(i,j,k) +      &
     &                  (1.0 - weight) * exner_rho_levels(i,j,k+1)
            End Do
          End Do
        End Do
        k = model_levels
        do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
            exner_theta_levels(i,j,k) =   0.5 *                         &
     &            (exner_rho_levels(i,j,k) + exner_rho_levels(i,j,k+1))
          End Do
        End Do

      EndIf !L_SH_Williamson .and. height_domain  <   45000.0

! ----------------------------------------------------------------------
! Section 1.6.4  Convert temperature to potential temperature
! ----------------------------------------------------------------------
      Do k = 1, model_levels
        do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
            theta(i,j,k) =  theta(i,j,k) / exner_theta_levels(i,j,k)
          End Do
        End Do
      End Do

!  exner_ref required to set lbcs - theta_ref contains temperature
      exner_ref(1) =  exp(- g * height_domain *                         &
     &                       eta_rho_levels(1) / ( Cp * theta_ref(1)))

      do k = 2,model_levels
        exner_ref(k)= exner_ref(k-1) *                                  &
     &                                 exp(g * height_domain *          &
     &                 (eta_rho_levels(k-1) - eta_rho_levels(k))        &
     &                                     / ( Cp * theta_ref(k-1)))
      end do

      k = model_levels + 1
      exner_ref(k) = exner_ref(k-1) * exp( g * height_domain *          &
     &           2.0 * (eta_rho_levels(k-1) - eta_theta_levels(k-1))    &
     &                                   / ( Cp * theta_ref(k-1)))

!  theta_ref contains temperature - convert to potential temperature
      do k = 1, model_levels
        theta_ref(k)= theta_ref(k) / exner_ref(k)
      end do

      IF (lhook) CALL dr_hook('IDL_TPROFILE6',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_tprofile6

!  End subroutine IDL_tprofile6

