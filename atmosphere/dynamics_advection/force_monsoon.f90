! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
! Subroutine Force_Monsoon

      Subroutine Force_Monsoon                                          &
     &                        (row_length                               &
     &,                        rows, model_levels, timestep             &
     &,                        delta_lambda, delta_phi, model_domain    &
     &,                        lambda_half_width, phi_half_width        &
     &,                        p_max, p_top, p_bottom, p_arbitrary      &
     &,                        lambda_heat_centre, phi_heat_centre      &
     &,                        max_heat_per_day, newtonian_timescale    &
     &,                        off_x, off_y, l_datastart                &
     &,                        at_extremity                             &
     &,                        exner_theta_levels, p, theta )

! Purpose:
!          Provides temperature forcing term for Idealised Monsoon
!          experiments.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.


      USE earth_constants_mod, ONLY: g
      USE atmos_constants_mod, ONLY: r, kappa, p_zero, recip_kappa

      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None



      Integer                                                           &
     &  model_domain

      Integer                                                           &
     &  row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, model_levels                                                    &
                         ! number of model levels
     &, off_x                                                           &
                     ! Size of small halo in i
     &, off_y                                                           &
                     ! Size of small halo in j.
     &, l_datastart(3)

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

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

      Real                                                              &
     &  timestep

      Real                                                              &
           ! horizontal co-ordinate information
     &  delta_lambda                                                    &
     &, delta_phi

      Real                                                              &
     &  p(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)   &
     &, theta (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)                                            &
     &, exner_theta_levels(1-off_x:row_length+off_x, 1-off_y:rows+off_y,&
     &         model_levels)

! Monsoon variables

      Real                                                              &
     &  lambda_half_width                                               &
     &, phi_half_width                                                  &
     &, p_max                                                           &
     &, p_top                                                           &
     &, p_bottom                                                        &
     &, p_arbitrary                                                     &
     &, lambda_heat_centre                                              &
     &, phi_heat_centre                                                 &
     &, max_heat_per_day                                                &
     &, newtonian_timescale

! local variables

      Integer                                                           &
     &  i, j, k, gi, gj, j0, j1

      Real                                                              &
     &  Tau_over_kappa                                                  &
     &, p_at_theta                                                      &
     &, lambda                                                          &
     &, lambda_term                                                     &
     &, phi                                                             &
     &, phi_term                                                        &
     &, theta_ref                                                       &
     &, recip_lambda_half_width                                         &
     &, recip_phi_half_width                                            &
     &, b_term                                                          &
     &, heating_function                                                &
     &, newtonian_cooling                                               &
     &, pressure_term                                                   &
     &, Base_phi                                                        &
     &, max_heat_per_second

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! No External routines

! ----------------------------------------------------------------------
! Section 1.  Calculate heat source and Newtonian cooling and add on.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('FORCE_MONSOON',zhook_in,zhook_handle)
      max_heat_per_second = max_heat_per_day / (24. * 3600.)

      If (model_domain  ==  1 ) Then
        Base_phi = - pi / 2
      Else
        Base_phi = 0.
      End If

      pressure_term = 1. / ( (p_max - p_top) * (p_max - p_bottom) *     &
     &                       (p_max - p_arbitrary) )
      recip_lambda_half_width = 1. / lambda_half_width
      recip_phi_half_width = 1. / phi_half_width
      Tau_over_kappa = 620. * g / (R * kappa)

      j0 = 1
      j1 = rows
      If (at_extremity(PSouth)) j0 = 2
      If (at_extremity(PNorth)) j1 = rows-1

      Do k = 1, model_levels
        Do j = j0, j1
          Do i = 1, row_length
            gi = l_datastart(1)+i-1
            gj = l_datastart(2)+j-1

            p_at_theta = exner_theta_levels(i,j,k) ** recip_kappa       &
     &                   * p_zero
            theta_ref = (p_zero/p_at_theta)**kappa * Tau_over_kappa     &
     &                      + 300. - Tau_over_kappa

            lambda = (gi-1) * delta_lambda
            lambda_term = ( (lambda - lambda_heat_centre) *             &
     &                       recip_lambda_half_width ) **2
            phi = Base_phi + (gj-1) * delta_phi
            phi_term = ( (phi - phi_heat_centre) *                      &
     &                       recip_phi_half_width ) **2

            b_term = max_heat_per_second *                              &
     &               exp( - lambda_term) * exp( - phi_term)

            If (p_at_theta  <   p_bottom .and.                          &
     &          p_at_theta  >   p_top ) Then
              heating_function = b_term * (p_at_theta - p_top) *        &
     &                                    (p_at_theta - p_bottom) *     &
     &                                    (p_at_theta - p_arbitrary) *  &
     &                                    pressure_term
            Else
              heating_function = 0.
            End If

            newtonian_cooling = newtonian_timescale *                   &
     &                          (theta(i,j,k) - theta_ref)

            theta(i,j,k) = theta(i,j,k) + timestep *                    &
     &                     (heating_function - newtonian_cooling)

          End Do
        End Do

      End Do

      If (at_extremity(PSouth)) Then
        gj = l_datastart(2)
        j = 1
        Do k = 1, model_levels
          Do i = 1, row_length

            p_at_theta = exner_theta_levels(i,j,k) ** recip_kappa       &
     &                   * p_zero
            theta_ref = (p_zero/p_at_theta)**kappa * Tau_over_kappa     &
     &                      + 300. - Tau_over_kappa

            lambda = 0.
            lambda_term = ( (lambda - lambda_heat_centre) *             &
     &                       recip_lambda_half_width ) **2
            phi = Base_phi + (gj-1) * delta_phi
            phi_term = ( (phi - phi_heat_centre) *                      &
     &                       recip_phi_half_width ) **2

            b_term = max_heat_per_second *                              &
     &               exp( - lambda_term) * exp( - phi_term)

            If (p_at_theta  <   p_bottom .and.                          &
     &          p_at_theta  >   p_top ) Then
              heating_function = b_term * (p_at_theta - p_top) *        &
     &                                    (p_at_theta - p_bottom) *     &
     &                                    (p_at_theta - p_arbitrary) *  &
     &                                    pressure_term
            Else
              heating_function = 0.
            End If

            newtonian_cooling = newtonian_timescale *                   &
     &                          (theta(i,j,k) - theta_ref)

            theta(i,j,k) = theta(i,j,k) + timestep *                    &
     &                     (heating_function - newtonian_cooling)

          End Do
        End Do

      End If

      If (at_extremity(PNorth)) Then
        j = rows
        gj = l_datastart(2) + j - 1
        Do k = 1, model_levels
          Do i = 1, row_length

            p_at_theta = exner_theta_levels(i,j,k) ** recip_kappa       &
     &                   * p_zero
            theta_ref = (p_zero/p_at_theta)**kappa * Tau_over_kappa     &
     &                      + 300. - Tau_over_kappa

            lambda = 0.
            lambda_term = ( (lambda - lambda_heat_centre) *             &
     &                       recip_lambda_half_width ) **2
            phi = Base_phi + (gj-1) * delta_phi
            phi_term = ( (phi - phi_heat_centre) *                      &
     &                       recip_phi_half_width ) **2

            b_term = max_heat_per_second *                              &
     &               exp( - lambda_term) * exp( - phi_term)

            If (p_at_theta  <   p_bottom .and.                          &
     &          p_at_theta  >   p_top ) Then
              heating_function = b_term * (p_at_theta - p_top) *        &
     &                                    (p_at_theta - p_bottom) *     &
     &                                    (p_at_theta - p_arbitrary) *  &
     &                                    pressure_term
            Else
              heating_function = 0.
            End If

            newtonian_cooling = newtonian_timescale *                   &
     &                          (theta(i,j,k) - theta_ref)

            theta(i,j,k) = theta(i,j,k) + timestep *                    &
     &                     (heating_function - newtonian_cooling)

          End Do
        End Do

      End If

      IF (lhook) CALL dr_hook('FORCE_MONSOON',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Force_Monsoon

