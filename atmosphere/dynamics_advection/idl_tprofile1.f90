! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine IDL_tprofile1

      Subroutine IDL_tprofile1(                                         &
     &                      row_length, rows, model_levels              &
     &,                     me, halo_i, halo_j                          &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     theta, exner_rho_levels, exner_theta_levels &
!  Profiles for fixed lbcs and sponge zones
     &,                     theta_ref, exner_ref                        &
!  Grid information
     &,                     height_domain, p_surface                    &
! Profile settings
     &,                     theta_surface, dtheta_dz1, height_dz1       &
!  Options
     &,                     L_constant_dz                               &
     &,                     L_code_test)

! Purpose:
!          Sets up initial data for idealised problems.
!          Initial temperature profile at all points is set
!          to the same everywhere
!          constant dtheta/dz.
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
     &, height_domain                                                   &
     &, dtheta_dz1(3)                                                   &
                        !Allows different values of dtheta_dz to be set
     &, height_dz1(2)   ! at different heights specified by the
                        ! height_dz variable.
      Logical                                                           &
     &  L_constant_dz                                                   &
                       ! Sets constant dtheta_dz from dtheta_dz1(1)
     &, L_code_test    ! user switch

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

! local variables
      Integer                                                           &
     &  i, j, k

      Real                                                              &
     & weight                                                           &
     &, exner_surface                                                   &
     &, r_theta_height_k, r_theta_height_km1

      Real                                                              &
     &  theta_temp(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)    &
     &, work1(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, dtheta_dz(model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------
! Section 1  constant static stability  dtheta_dz
!-------------------------------------------------------------------
      IF (lhook) CALL dr_hook('IDL_TPROFILE1',zhook_in,zhook_handle)
      exner_surface = (p_surface/p_zero) ** kappa

      If (me  ==  0) Then
        If (L_constant_dz) Then
          If (dtheta_dz1(1) == 0.0) Then
            Write (Unit=6,Fmt='(A21,A14,F6.2,A2)')                      &
     &                  '   Isentropic profile',                        &
     &                  '   with theta = ',theta_surface,' K'
          Else
            Write (6,*) '   Constant static stability profile',         &
     &                  '   with surface theta = ',theta_surface,' K',  &
     &                  '   and dtheta/dz = ',dtheta_dz1(1),' K/m'
          End If
        Else
          Print*,'***** constant static stability layers *****'
        End If    !(L_constant_dz)
      Endif   !(me  ==  0)

      Do k = 1, model_levels
        r_theta_height_km1 = r_theta_levels(1,1,k-1)-Earth_radius
        r_theta_height_k   = r_theta_levels(1,1,k)-Earth_radius
        If (L_constant_dz) Then
          dtheta_dz(k) = dtheta_dz1(1)
        Else
! set up a profile if required
          If ((r_theta_levels(1,1,k) - Earth_radius)                    &
     &           <= height_dz1(1)) Then

            dtheta_dz(k) = dtheta_dz1(1)

          Else If (r_theta_height_k   > height_dz1(1)                   &
     &       .and. r_theta_height_km1 < height_dz1(1))                  &
     &        Then  ! height_dz1 is inbetween levels

            dtheta_dz(k) =                                              &
     &        (dtheta_dz1(1)*(height_dz1(1)-r_theta_height_km1)         &
     &       + dtheta_dz1(2)*(r_theta_height_k-height_dz1(1)))          &
     &        /(r_theta_height_k - r_theta_height_km1)

          Else If ((r_theta_levels(1,1,k) - Earth_radius)               &
     &               <= height_dz1(2)) Then

            dtheta_dz(k) = dtheta_dz1(2)

          Else If (r_theta_height_k   > height_dz1(2)                   &
     &       .and. r_theta_height_km1 < height_dz1(2))                  &
     &        Then  ! height_dz2 is inbetween levels

            dtheta_dz(k) =                                              &
     &        (dtheta_dz1(2)*(height_dz1(2)-r_theta_height_km1)         &
     &       + dtheta_dz1(3)*(r_theta_height_k-height_dz1(2)))          &
     &        /(r_theta_height_k - r_theta_height_km1)

          Else

            dtheta_dz(k) = dtheta_dz1(3)

          End If

        End If    !(L_constant_dz)
      End Do

! Set theta field.
!  First calculate theta_temp, theta on orography
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          theta_temp(i,j) = theta_surface + dtheta_dz(1) *              &
     &                      (r_theta_levels(i,j,0) - Earth_radius)
        end do
      end do

      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          theta(i,j,1) = theta_temp(i,j) + dtheta_dz(1) *               &
     &               (r_theta_levels(i,j,1) - r_theta_levels(i,j,0))
        end do
      end do

      do k = 2, model_levels
        do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            theta(i,j,k) = theta(i,j,k-1) + dtheta_dz(k) *              &
     &             (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1))
          end do
        end do
      end do

! Calculate theta_ref profile to be bit-reproducible for all processors
!  This is needed for sponge.  Assumes undisturbed (no-orography) flow

      theta_ref(1) = theta_surface + dtheta_dz(1) * height_domain *     &
     &                       eta_theta_levels(1)
      do k = 2, model_levels
        theta_ref(k) = theta_ref(k-1) + dtheta_dz(k) * height_domain *  &
     &                (eta_theta_levels(k) - eta_theta_levels(k-1))
      end do

!-----------------------------------------------------------------
! Generate pressure. (From hydrostatic relation for Exner pressure)
!                     work1 is exner_surface
!-------------------------------------------------------------------
!  Need theta_temp, theta midway between orography and z = 0
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          theta_temp(i,j) = theta_surface + dtheta_dz(1) * 0.5 *        &
     &                    (r_theta_levels(i,j,0) - Earth_radius)
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
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          theta_temp(i,j) = theta_surface + dtheta_dz(1) * 0.5 *        &
     &                    (r_rho_levels(i,j,1) - r_theta_levels(i,j,0))
        end do
      end do

      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          exner_rho_levels(i,j,1) = work1(i,j) - g *                    &
     &                    (r_rho_levels(i,j,1) - r_theta_levels(i,j,0)) &
     &                               / ( Cp * theta_temp(i,j) )
        end do
      end do

      do k = 2,model_levels
        do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            exner_rho_levels(i,j,k)= exner_rho_levels(i,j,k-1)-         &
     &            g * (r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1))     &
     &                                     / ( Cp * theta(i,j,k-1))
          end do
        end do
      end do

      k = model_levels + 1
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          exner_rho_levels(i,j,k) = exner_rho_levels(i,j,k-1) -         &
     &    2.0 * g * (r_theta_levels(i,j,k-1) - r_rho_levels(i,j,k-1))   &
     &                                   / ( Cp * theta(i,j,k-1) )
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
          End Do
        End Do
      End Do

      k =  model_levels
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          exner_theta_levels(i,j,k) =  0.5 * exner_rho_levels(i,j,k) +  &
     &                                  0.5 * exner_rho_levels(i,j,k+1)
        End Do
      End Do

!-----------------------------------------------------------------
! Generate exner_ref required to set lbcs
!-------------------------------------------------------------------

      exner_ref(1) = exner_surface - g *                                &
     &                             height_domain * eta_rho_levels(1)    &
     &                               / ( Cp * theta_ref(1) )

      do k = 2,model_levels
        exner_ref(k)= exner_ref(k-1)- g * height_domain *               &
     &                 (eta_rho_levels(k) - eta_rho_levels(k-1))        &
     &                                     / ( Cp * theta_ref(k-1))
      end do

      k = model_levels + 1
      exner_ref(k) = exner_ref(k-1) - 2.0 * g * height_domain *         &
     &              (eta_theta_levels(k-1) - eta_rho_levels(k-1))       &
     &                                   / ( Cp * theta_ref(k-1) )

      IF (lhook) CALL dr_hook('IDL_TPROFILE1',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_tprofile1

!  End subroutine IDL_tprofile1

