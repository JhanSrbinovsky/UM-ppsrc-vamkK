! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************


! Subroutine IDL_tprofile4

      Subroutine IDL_tprofile4(                                         &
     &                      row_length, rows, model_levels              &
     &,                     me, halo_i, halo_j                          &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     theta, exner_rho_levels, exner_theta_levels &
!  Profiles for fixed lbcs and sponge zones
     &,                     theta_ref, exner_ref                        &
!  Grid information
     &,                     height_domain, big_layers                   &
! Profile settings
     &,                     p_surface, theta_surface, Brunt_Vaisala     &
!  Options
     &,                     L_code_test)

! Purpose:
!          Sets up initial data for idealised problems.
!          Initial temperature profile at all points are set
!          to the same everywhere
!          Constant Brunt-Vaisala frequency
!              PLUS isothermal big_layers
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
          r, cp, kappa, p_zero 

      USE earth_constants_mod, ONLY: g, earth_radius
      
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
     &, halo_j               ! Size of halo in j direction.

      Real                                                              &
     &  theta_surface                                                   &
     &, p_surface                                                       &
     &, Brunt_Vaisala                                                   &
     &, height_domain

      Logical                                                           &
     &  L_code_test    ! user switch

      Integer                                                           &
     &  me         ! My processor number


      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, model_levels)                &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

! Primary Arrays used in all models
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

! local variables
      Integer                                                           &
     &  i, j, k

      Real                                                              &
     &  weight                                                          &
     &, exner_surface                                                   &
     &, delta_z                                                         &
     &, temp                                                            &
     &, BV_squared_over_g

      Real                                                              &
     &  exner_theta_ref(model_levels)                                   &
                                       ! Exner ref profile
     &, theta_temp(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)    &
     &, work1(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, work2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------
! Section 1  constant Brunt-Vaisala frequency
!        ** PLUS isothermal big_layers ***
!-------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_TPROFILE4',zhook_in,zhook_handle)
      BV_squared_over_g =  Brunt_Vaisala * Brunt_Vaisala / g

      exner_surface = (p_surface/p_zero) ** kappa

      if(me  ==  0)then
        print*,'** Constant Brunt-Vaisala frequency  profile ***'
        print*,'** PLUS isothermal big_layers ***'
      Endif   !(me  ==  0)

! Set theta field up to model_levels - big_layers.

      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          theta(i,j,1) = theta_surface *                                &
     &                     exp ( BV_squared_over_g *                    &
     &              (r_theta_levels(i,j,1) - Earth_radius ))
        end do
      end do

      do k = 2, model_levels - big_layers
        do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            theta(i,j,k) = theta(i,j,k-1) *                             &
     &                     exp ( BV_squared_over_g *                    &
     &             (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)))
          end do
        end do
      end do
!-----------------------------------------------------------------
! Generate pressure. (From hydrostatic relation for Exner pressure
!                     work1 contains Exner_surface
!-------------------------------------------------------------------
!  Need theta_temp, theta midway between orography and z = 0
!   Hence the factor 0.5 multiplying deltaz in exp function
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          theta_temp(i,j) = theta_surface *                             &
     &                     exp ( BV_squared_over_g * 0.5 *              &
     &              (r_theta_levels(i,j,0) - Earth_radius ))
        end do
      end do

      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          work1(i,j) = exner_surface - g *                              &
     &                         (r_theta_levels(i,j,0) - Earth_radius )  &
     &                               / ( Cp * theta_temp(i,j) )
        end do
      end do

!  Need theta_temp, theta midway between rho level and orography
!   Hence the factor 0.5 multiplying deltaz in exp function
!   and the plus sign since it is relative to z=0.
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          theta_temp(i,j) = theta_surface *                             &
     &                     exp ( BV_squared_over_g * ( 0.5 *            &
     &              (r_rho_levels(i,j,1) + r_theta_levels(i,j,0))       &
     &                                   -  Earth_radius ) )
        end do
      end do

      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          exner_rho_levels(i,j,1) = work1(i,j) - g *                    &
     &                    (r_rho_levels(i,j,1) - r_theta_levels(i,j,0)) &
     &                               / ( Cp * theta_temp(i,j) )
        end do
      end do

      do k = 2, model_levels - big_layers + 1
        do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            exner_rho_levels(i,j,k)= exner_rho_levels(i,j,k-1)-         &
     &            g * (r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1))     &
     &                                     / ( Cp * theta(i,j,k-1))
          end do
        end do
      end do

      ! Calculate exner_theta_levels up to model_levels-big_layers-1
      Do k = 1, model_levels - big_layers - 1
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            weight = (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))/   &
     &                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
            exner_theta_levels(i,j,k) =   weight *                      &
     &                                   exner_rho_levels(i,j,k) +      &
     &                  (1.0 - weight) * exner_rho_levels(i,j,k+1)
          End Do
        End Do
      End Do

!     need to calculate Exner pressure on model_levels - big_layers
!     to get temperature  at start of big_layers
      k = model_levels - big_layers
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          delta_z = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
          weight = (r_rho_levels(i,j,k+1 ) - r_theta_levels(i,j,k))     &
     &                     /  delta_z

          exner_theta_levels(i,j,k) =  exner_rho_levels(i,j,k+1)*       &
     &                                weight +                          &
     &                                exner_rho_levels(i,j,k) *         &
     &                                (1.0 - weight)
        End Do
      End Do

!     now calculate  temperature at  model_levels - big_layers
!       store in work2 work array

      k= model_levels - big_layers
      do j = 1-halo_j, rows+halo_j
        do i = 1-halo_i, row_length+halo_i
          work2(i,j) = exner_theta_levels(i,j,k) * theta(i,j,k)
        end do
      end do

!     and calculate exner pressure on  rest of big_layers
!     now temperature (work2) is known. Only needed if big_layers > 1
      if(big_layers  >   1)then

        do k = model_levels - big_layers + 1, model_levels
          do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              delta_z = r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1)
              weight = (r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1))  &
     &                     /  delta_z

              exner_rho_levels(i,j,k) =  exner_rho_levels(i,j,k-1) *    &
     &            ( Cp * work2(i,j) -                                   &
     &               (1.0 - weight) * g * delta_z )/                    &
     &            ( Cp * work2(i,j) + weight * g * delta_z )
            end do
          end do
        end do

      endif   !      (big_layers  >   1)then

      k = model_levels + 1
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          delta_z = r_theta_levels(i,j,k-1) - r_rho_levels(i,j,k-1)

          exner_rho_levels(i,j,k) =  exner_rho_levels(i,j,k-1) *        &
     &    ( Cp * work2(i,j) - g * delta_z )/                            &
     &    ( Cp * work2(i,j) + g * delta_z )
        end do
      end do

!***************************************************
! Now calculate exner_theta_levels for big-layers
      Do k = model_levels - big_layers + 1, model_levels
        do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            if(k  ==  model_levels)then
             weight = 0.5
            else
              delta_z = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
              weight = (r_rho_levels(i,j,k+1 ) - r_theta_levels(i,j,k)) &
     &                     /  delta_z
            endif    !(k  ==  model_levels)
            exner_theta_levels(i,j,k) =  exner_rho_levels(i,j,k+1) *    &
     &                                weight +                          &
     &                                exner_rho_levels(i,j,k) *         &
     &                                (1.0 - weight)
          End Do
        End Do
      End Do

!*****************************************************
!     now calculate  thetas on  big_layers

      Do k = model_levels - big_layers + 1, model_levels
        do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            theta(i,j,k) = work2(i,j) /                                 &
     &                      exner_theta_levels(i,j,k)
          end do
        end do
      End do

! Calculate theta_ref profile to be bit-reproducible for all processors
!  This is needed for sponge.  Assumes undisturbed (no-orography) flow

      theta_ref(1) =  theta_surface *                                   &
     &                    exp ( BV_squared_over_g *                     &
     &                       eta_theta_levels(1) * height_domain)
      do k = 2, model_levels - big_layers
        theta_ref(k) = theta_ref(k-1) *                                 &
     &                  exp ( BV_squared_over_g * height_domain *       &
     &                (eta_theta_levels(k) - eta_theta_levels(k-1)))
      end do

!  Need theta_temp, theta midway between rho level and orography
!   Hence the factor 0.5 multiplying deltaz in exp function
      temp = theta_surface *                                            &
     &                     exp ( BV_squared_over_g * 0.5 *              &
     &                        eta_rho_levels(1) * height_domain )

!   Use profile array to store reference exner-pressure
      exner_ref(1) = exner_surface - g *                                &
     &                    eta_rho_levels(1) * height_domain             &
     &                               / ( Cp * temp )

      do k = 2, model_levels - big_layers + 1
        exner_ref(k)= exner_ref(k-1) - g * height_domain *              &
     &                 (eta_rho_levels(k) - eta_rho_levels(k-1))        &
     &                                     / ( Cp * theta_ref(k-1))
      end do

!     need to calculate Exner pressure on model_levels - big_layers
!     to get temperature  at start of big_layers
      k = model_levels - big_layers
      delta_z = eta_rho_levels(k+1) - eta_rho_levels(k)
      weight = (eta_rho_levels(k+1) - eta_theta_levels(k))              &
     &                     /  delta_z

      exner_theta_ref(k) =  exner_ref(k+1) * weight +                   &
     &                         exner_ref(k) * (1.0 - weight)

!     now calculate  temperature at  model_levels - big_layers

      k= model_levels - big_layers
      temp = exner_theta_ref(k) * theta_ref(k)

!     and calculate exner pressure on  rest of big_layers
!     now temperature is known. Only needed if big_layers > 1

      if(big_layers  >   1)then

        do k = model_levels - big_layers + 1, model_levels
          delta_z = eta_rho_levels(k) - eta_rho_levels(k-1)
          weight = (eta_rho_levels(k) - eta_theta_levels(k-1))          &
     &                     /  delta_z
          delta_z = delta_z * height_domain

          exner_ref(k) =  exner_ref(k-1) *                              &
     &    ( Cp * temp - (1.0 - weight) * g * delta_z ) /                &
     &    ( Cp * temp +        weight  * g * delta_z )
        end do

      endif   !      (big_layers  >   1)then

      k = model_levels + 1
      delta_z = (eta_theta_levels(k-1) - eta_rho_levels(k-1)) *         &
     &                                  height_domain

      exner_ref(k) =  exner_ref(k-1) *                                  &
     &     ( Cp * temp - g * delta_z ) / ( Cp * temp + g * delta_z )

!***************************************************
! Now calculate exner_theta_levels for big-layers
      Do k = model_levels - big_layers + 1, model_levels
        if(k  ==  model_levels)then
          weight = 0.5
        else
          delta_z = eta_rho_levels(k+1) - eta_rho_levels(k)
          weight = (eta_rho_levels(k+1) - eta_theta_levels(k))          &
     &                     /  delta_z
        endif    !(k  ==  model_levels)
          exner_theta_ref(k) =  exner_ref(k+1) * weight +               &
     &                           exner_ref(k) * (1.0 - weight)
      End Do

!*****************************************************
!     now calculate  thetas on  big_layers

      Do k = model_levels - big_layers + 1, model_levels
        theta_ref(k) = temp / exner_theta_ref(k)
      End do

      IF (lhook) CALL dr_hook('IDL_TPROFILE4',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_tprofile4
!  End subroutine IDL_tprofile4

