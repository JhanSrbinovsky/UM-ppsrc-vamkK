! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine IDL_tprofile2

      Subroutine IDL_tprofile2(                                         &
     &                      row_length, rows, model_levels              &
     &,                     me, halo_i, halo_j                          &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     theta, exner_rho_levels, exner_theta_levels &
     &,                     theta_ref, exner_ref                        &
     &,                     height_domain                               &
     &,                     p_surface, theta_surface                    &
     &,                     L_code_test)

! Purpose:
!          Sets up initial data for idealised problems.
!          Initial isothermal temperature profile at all points
!          is set to the same everywhere
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE earth_constants_mod, ONLY: g, earth_radius
      USE atmos_constants_mod, ONLY:                                    &
          r, cp, kappa, p_zero 

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
     &, halo_i                                                          &
                             ! Size of halo in i direction.
     &, halo_j               ! Size of halo in j direction.

      Real                                                              &
     &  theta_surface                                                   &
     &, p_surface                                                       &
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

! Output Arrays from this routine
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
     &, temp                                                            &
     &, goverCpT                                                        &
     &, delta_z                                                         &
     &, exner_surface

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------
! Section 1  ISOTHERMAL   temperature constant
!        (constant) temperature =  theta_surface
!-------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_TPROFILE2',zhook_in,zhook_handle)
      if(me  ==  0)then
        print*,'***** ISOTHERMAL  profile *****'
      Endif   !(me  ==  0)

      goverCpT = g/(Cp * theta_surface)
      exner_surface = (p_surface/p_zero) ** kappa

      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          delta_z = 0.5 * (r_rho_levels(i,j,1) - Earth_radius)
          exner_rho_levels(i,j,1) =  exner_surface *                    &
     &                              (1.0 - goverCpT * delta_z ) /       &
     &                               (1.0 + goverCpT * delta_z )
        end do
      end do

      do k = 2, model_levels
        do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            weight = (r_theta_levels(i,j,k-1) - r_rho_levels(i,j,k-1))
            temp = (r_rho_levels(i,j,k) - r_theta_levels(i,j,k-1))
             exner_rho_levels(i,j,k) =  exner_rho_levels(i,j,k-1) *     &
     &                              ( 1.0 -  temp * goverCpT  ) /       &
     &                              ( 1.0 +  weight * goverCpT  )
          end do
        end do
      end do

      k = model_levels + 1
! weight = 0.5 cancels factor of 2 in delta_z
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          delta_z =  r_theta_levels(i,j,k-1) - r_rho_levels(i,j,k-1)
          exner_rho_levels(i,j,k) =  exner_rho_levels(i,j,k-1) *        &
     &                            ( 1.0 -  goverCpT * delta_z ) /       &
     &                            ( 1.0 +  goverCpT * delta_z )
        end do
      end do

      Do k = 1, model_levels - 1
        do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            weight = (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))/   &
     &                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
            exner_theta_levels(i,j,k) =   weight *                      &
     &                                   exner_rho_levels(i,j,k) +      &
     &                  (1.0 - weight) * exner_rho_levels(i,j,k+1)
            theta(i,j,k) =   theta_surface /                            &
     &                                   exner_theta_levels(i,j,k)
          End Do
        End Do
      End Do

      k = model_levels
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          exner_theta_levels(i,j,k) =  0.5 * exner_rho_levels(i,j,k) +  &
     &                                  0.5 * exner_rho_levels(i,j,k+1)
          theta(i,j,k) =   theta_surface /                              &
     &                                   exner_theta_levels(i,j,k)
        End Do
      End Do

!  exner_ref required to set lbcs
      delta_z = 0.5 * height_domain * eta_rho_levels(1)
      exner_ref(1) =  (1.0 - goverCpT * delta_z ) /                     &
     &                               (1.0 + goverCpT * delta_z )

      do k = 2, model_levels
        weight = (eta_theta_levels(k-1) - eta_rho_levels(k-1)) *        &
     &                                              height_domain
        temp = (eta_rho_levels(k) - eta_theta_levels(k-1)) *            &
     &                                              height_domain
        exner_ref(k) =  exner_ref(k-1) *                                &
     &                              ( 1.0 -  temp * goverCpT  ) /       &
     &                              ( 1.0 +  weight * goverCpT  )
      end do

      k = model_levels + 1
! weight = 0.5 cancels factor of 2 in delta_z
      delta_z =  (eta_theta_levels(k-1) - eta_rho_levels(k-1)) *        &
     &                                              height_domain
      exner_ref(k) =  exner_ref(k-1) *                                  &
     &                            ( 1.0 -  goverCpT * delta_z ) /       &
     &                            ( 1.0 +  goverCpT * delta_z )

!  theta_ref required to set lbcs and for sponge region
      Do k = 1, model_levels - 1
        weight = (eta_rho_levels(k+1) - eta_theta_levels(k)) /          &
     &             (eta_rho_levels(k+1) - eta_rho_levels(k) )
        temp =              weight * exner_ref(k) +                     &
     &                   (1.0 - weight) * exner_ref(k+1)
        theta_ref(k) =   theta_surface /temp
      End Do

      k = model_levels
      temp =  0.5 * (exner_ref(k) + exner_ref(k+1) )
      theta_ref(k) =   theta_surface / temp

      IF (lhook) CALL dr_hook('IDL_TPROFILE2',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_tprofile2
!  End subroutine IDL_tprofile2

