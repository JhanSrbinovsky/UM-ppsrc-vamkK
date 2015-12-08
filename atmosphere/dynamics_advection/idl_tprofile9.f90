! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine IDL_tprofile9

      Subroutine IDL_tprofile9(                                         &
     &                      row_length, rows, model_levels              &
     &,                     me, halo_i, halo_j                          &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     theta, exner_rho_levels, exner_theta_levels &
!  Profiles for fixed lbcs and sponge zones
     &,                     theta_ref, exner_ref                        &
!  Grid information
     &,                     height_domain, big_layers                   &
     &,                     surface_type                                &
! Profile settings
     &,                     zprofile_data, tprofile_data                &
     &,                     p_surface, theta_surface, num_profile_data  &
     &,                     zprofile_orog, idl_interp_option, hf        &
!  Options
     &,                     L_code_test)

! Purpose:
!          Sets up initial data for idealised problems.
!          Input from NAMELIST input
!
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
      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE

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
     &, height_domain                                                   &
     &, zprofile_orog                                                   &
     &, hf                                                              &
     &, zprofile_data(100)                                              &
                           ! Heights for theta profiles
     &, tprofile_data(100) ! Data values for theta profile

      Integer                                                           &
     &  surface_type                                                    &
                           ! idealised orography type
     &, num_profile_data                                                &
                           ! number of values in z,q,tprofile_data
     &, idl_interp_option  ! Profile interpolation option

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
     &  i, j, k, k2

      Real                                                              &
     &  weight                                                          &
     &, exner_surface                                                   &
     &, theta_ref_at_rho1                                               &
     &, z_at_theta, z_at_theta_p1                                       &
     &, z_at_rho,   z_at_rho_p1                                         &
     &, z_at_orog                                                       &
     &, eta_model, out1, out2, hs                                       &
     &, eta_profile(num_profile_data)


      Real                                                              &
     &  theta_temp(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)    &
     &, theta_orog(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)    &
     &, theta_rholev1(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j) &
     &, exner_orog(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)    &
     &, work1(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)

      ! Error reporting
      Character (Len=*),  Parameter :: RoutineName='idl_tprofile9'
      Character (Len=256)           :: Cmessage
      Integer                       :: ErrorStatus

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------
! Section 1  set profile from namelist data
!
!-------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_TPROFILE9',zhook_in,zhook_handle)
      exner_surface = (p_surface/p_zero) ** kappa

      ! Check to make sure the namelist profile data extends
      ! to the top of the model.
      Do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          If (r_theta_levels(i,j,model_levels) - Earth_radius           &
     &         >   zprofile_data(num_profile_data)) Then
            Write(Cmessage,*)                                           &
     &        'Idealised namelist vertical profile data'                &
     &        //'does not extend to the top of the model.'              &
     &        //'Please modify the namelist data.'
            ErrorStatus = 1

            Call Ereport( RoutineName, ErrorStatus, Cmessage )
          End If
        End Do
      End Do

      ! Interpolate theta from namelist profile to model levels
      Do k = 1, model_levels
        Do k2 = 1, num_profile_data-1
          Do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i

              z_at_theta = r_theta_levels(i,j,k) - Earth_radius

              If (z_at_theta  >   zprofile_data(k2) .and.               &
     &            z_at_theta  <=  zprofile_data(k2+1)) Then

                weight = (z_at_theta - zprofile_data(k2))               &
     &                 /(zprofile_data(k2+1) - zprofile_data(k2))

                theta(i,j,k) = tprofile_data(k2) + weight*              &
     &                   (tprofile_data(k2+1) - tprofile_data(k2))
              End If

            End Do
          End Do
        End Do
      End Do


      !-----------------------------------------------------------------
      ! Alternative interpolation of initial profile over orography
      !-----------------------------------------------------------------
      !
      !  idl_interp_option = 1: constant on height levels
      !                         (default above, no need to modify)
      !  idl_interp_option = 2: hybrid height everywhere up to a
      !                         specified height "hf" (the same as model
      !                         levels are defined). hs = height_domain
      !  idl_interp_option = 3: as option 2 but only when orography is
      !                         less than input profile orography.
      !                         hs = zprofile_orog
      !
      ! Sets up an eta coord for each level of the initial profile data
      ! and an eta coordinate for each model column, and interpolates
      ! in eta space if the model level height is less than "hf" and
      ! the model orography height is less than "hs"
      !-----------------------------------------------------------------

      If (idl_interp_option == 2 .or. idl_interp_option == 3) Then

        If (idl_interp_option == 2) hs = height_domain
        If (idl_interp_option == 3) hs = zprofile_orog

        ! Set up an eta coord for each level of the initial profile data
        eta_profile(1) = 0.0
        Do k2 = 2, num_profile_data
          eta_profile(k2) =  (zprofile_data(k2) - zprofile_orog)        &
     &                      /(hf - zprofile_orog)
        End Do

        ! Interpolate in eta space and overwrite theta where appropriate
        Do k = 1, model_levels
          Do k2 = 1, num_profile_data - 1
            Do j = 1-halo_j, rows+halo_j
              Do i = 1-halo_i, row_length+halo_i
                z_at_theta = r_theta_levels(i,j,k) - Earth_radius
                z_at_orog  = r_theta_levels(i,j,0) - Earth_radius
                eta_model = (z_at_theta - z_at_orog)/(hf - z_at_orog)

                If ( (z_at_orog <= hs ) .and. (z_at_theta < hf) .and.   &
     &               (eta_model > eta_profile(k2))  .and.               &
     &               (eta_model <= eta_profile(k2+1)) ) Then

                  weight = (eta_model - eta_profile(k2))                &
     &                     /(eta_profile(k2+1) - eta_profile(k2))
                  theta(i,j,k) = tprofile_data(k2) + weight*            &
     &                      (tprofile_data(k2+1) - tprofile_data(k2))
                End If
              End Do
            End Do
          End Do
        End Do

      End If ! on idl_interp_option = 2 or 3

      ! Set up reference theta profile (i.e. no orography)
      Do k = 1, model_levels
        z_at_theta = eta_theta_levels(k)*height_domain
        Do k2 = 1, num_profile_data-1
          If (z_at_theta  >   zprofile_data(k2) .and.                   &
     &        z_at_theta  <=  zprofile_data(k2+1)) Then

            weight = (z_at_theta - zprofile_data(k2))                   &
     &                 /(zprofile_data(k2+1) - zprofile_data(k2))

            theta_ref(k) = tprofile_data(k2) + weight*                  &
     &                   (tprofile_data(k2+1) - tprofile_data(k2))
          End If
        End Do
      End Do


!-----------------------------------------------------------------
! Generate pressure. (From hydrostatic relation for
!                     Exner pressure)
!-------------------------------------------------------------------
!  First set up reference exner profile (i.e. no orography) on rho levs

      ! Calculate theta_ref on first rho level
      weight = (eta_rho_levels(1)*height_domain)                        &
     &                 /(eta_theta_levels(1)*height_domain)
      theta_ref_at_rho1 = theta_surface + weight*                       &
     &                    (theta_ref(1) - theta_surface)

      ! Calculate exner_ref on first rho levels using average theta
      exner_ref(1) = exner_surface - g *                                &
     &               eta_rho_levels(1) * height_domain /                &
     &               (Cp * 0.5*(theta_ref_at_rho1+theta_surface) )

      Do k = 2,model_levels
        exner_ref(k)= exner_ref(k-1) - g * height_domain *              &
     &                 (eta_rho_levels(k) - eta_rho_levels(k-1))        &
     &                                     / ( Cp * theta_ref(k-1))
      End Do

      k = model_levels + 1
      exner_ref(k) = exner_ref(k-1) - 2.0 * g * height_domain *         &
     &              (eta_theta_levels(k-1) - eta_rho_levels(k-1))       &
     &                                   / ( Cp * theta_ref(k-1) )

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
            If (k  ==  1 .and. z_at_orog  <   z_at_rho) Then
              ! Set exner_orog if below rho level 1
              weight = z_at_orog/z_at_rho
              exner_orog(i,j) = exner_surface + weight*                 &
     &                       (exner_ref(1) - exner_surface)
            Else If (z_at_orog  >=  z_at_rho .and.                      &
     &               z_at_orog  <   z_at_rho_p1) Then
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
! Use theta_ref profile rather than tprofile_data.

      Do k = 1, model_levels - 1
        z_at_theta  = eta_theta_levels(k)*height_domain
        z_at_theta_p1 = eta_theta_levels(k+1)*height_domain
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            z_at_orog  = r_theta_levels(i,j,0) - Earth_radius
            If (k  ==  1 .and. z_at_orog  <   z_at_theta) Then
              ! Set theta_orog if below theta level 1
              weight = z_at_orog/z_at_theta
              theta_orog(i,j) = theta_surface + weight*                 &
     &                       (theta_ref(1) - theta_surface)
            Else If (z_at_orog  >=  z_at_theta .and.                    &
     &               z_at_orog  <   z_at_theta_p1) Then
              weight = (z_at_orog - z_at_theta)                         &
     &                 /(z_at_theta_p1 - z_at_theta)
              theta_orog(i,j) = theta_ref(k) + weight*                  &
     &                        (theta_ref(k+1) - theta_ref(k))
            End If
          End Do
        End Do
      End Do

! Find theta at r_rho_levels(i,j,1), i.e. first rho levels
! and put in theta_rholev1(i,j)
! Use theta_ref profile rather than tprofile_data.

      Do k = 1, model_levels - 1
        z_at_theta  = eta_theta_levels(k)*height_domain
        z_at_theta_p1 = eta_theta_levels(k+1)*height_domain
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            z_at_rho  = r_rho_levels(i,j,1) - Earth_radius
            If (k  ==  1 .and. z_at_rho  <   z_at_theta) Then
              ! Set theta_rholev1 if below rho level 1
              weight = z_at_rho/z_at_theta
              theta_rholev1(i,j) = theta_surface + weight*              &
     &                       (theta_ref(1) - theta_surface)
            Else If (z_at_rho  >=  z_at_theta .and.                     &
     &               z_at_rho  <   z_at_theta_p1) Then
              weight = (z_at_rho - z_at_theta)                          &
     &                 /(z_at_theta_p1 - z_at_theta)
              theta_rholev1(i,j) = theta_ref(k) + weight*               &
     &                        (theta_ref(k+1) - theta_ref(k))
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


      do j = 1-halo_j, rows+halo_j
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
     &                                     / ( Cp * theta(i,j,k-1))
          end do
        end do
      end do

      k = model_levels + 1
      do j = 1-halo_j, rows+halo_j
        Do i = 1-halo_i, row_length+halo_i
          exner_rho_levels(i,j,k) = exner_rho_levels(i,j,k-1) -         &
     &     2.0 * g * (r_theta_levels(i,j,k-1) - r_rho_levels(i,j,k-1))  &
     &                                   / ( Cp * theta(i,j,k-1) )
        end do
      end do

      Do k = 1, model_levels - 1
        do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            weight = (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))/   &
     &               (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
            exner_theta_levels(i,j,k) =   weight *                      &
     &                                  exner_rho_levels(i,j,k) +       &
     &                 (1.0 - weight) * exner_rho_levels(i,j,k+1)
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


      IF (lhook) CALL dr_hook('IDL_TPROFILE9',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_tprofile9

      !  End subroutine IDL_tprofile9

