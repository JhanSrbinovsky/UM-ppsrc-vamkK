! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
      Subroutine eta_vert_weights_level (                               &
     &                        eta_in,                                   &
     &                        r_in,                                     &
     &                        check_bottom_levels,                      &
     &                        interp_vertical_search_tol,               &
     &                        first_flat_level_in,                      &
     &                        dim_i_in, dim_j_in, dim_k_in,             &
     &                        dim_i_out, dim_j_out, dim_k_out,          &
     &                        high_order_scheme, monotone_scheme,       &
     &                        model_domain, L_high, L_mono,             &
     &                        L_conserv,                                &
     &                        halo_i, halo_j,                           &
     &                        r_out,                                    &
     &                        coeff_z, coeff_z_lin,                     &
     &                        k_out, i_out, j_out,                      &
     &                        weight_lambda, weight_phi, k)

! Purpose: Calculates vertical interpolation weights and level which
!          is below the desired point.
! This version of vert_weights calculates vertical interpolation
! weights appropriate to doing interpolation in eta coordinate,
! rather than r.
! Eta coordinates for departure points are computed by assuming that
! r_out(i,j,k) is the trilinear interpolant of r in (lambda,phi,eta)
! coordinates. This is easily inverted to obtain an interpolation
! weight for the eta coordinate. This weight is then used to linearly
! interpolate the eta coordinate to the departure point level.
! Ideally the departure point routine should return eta coordinates
! of the departure points, rather than r coordinates. This version of
! vert_weights is simply for testing the use of the eta coordinate in
! the interpolation subroutines.
!
! NB. This approach doesn't give the full accuracy of working
!     wholly in the eta coord system, as the location of departure
!     points in eta coords has to be obtained by interpolation.
!     One result of this is that trilinear interpolation using this
!     routine is just the same as with the standard r-coordinate
!     routine.
!
!
! Method:
!          Is described in ;
!          The proposed semi-Lagrangian advection scheme for the
!          semi-Implicit Unified Model integration scheme.
!          F.R. Division working paper No 162.
!          Mark H. Mawson
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      USE um_types,   ONLY: integer32
      USE yomhook,    ONLY: lhook, dr_hook
      USE parkind1,   ONLY: jprb, jpim
      USE highos_mod, ONLY: cubicLagrange, quinticLagrange,             &
                            ECMWF_quasiCubic, ECMWF_mono_quasiCubic,    &
                            hCubic_vLin, hQuasiCubic_vQuintic,          &
                            hCubic_vQuintic
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  dim_i_in                                                        &
                    ! Dimension of Data_in in i direction.
     &, dim_j_in                                                        &
                    ! Dimension of Data_in in j direction.
     &, dim_k_in                                                        &
                    ! Dimension of Data_in in k direction.
     &, dim_i_out                                                       &
                    ! Dimension of Data_out in i direction.
     &, dim_j_out                                                       &
                    ! Dimension of Data_out in j direction.
     &, dim_k_out                                                       &
                    ! Dimension of Data_out in k direction.
     &, halo_i                                                          &
                    ! Size of halo in i direction.
     &, halo_j      ! Size of halo in j direction.

      Integer                                                           &
     &  high_order_scheme                                               &
                           ! a code saying which high order scheme to
                           ! use.
     &, monotone_scheme                                                 &
                        ! a code saying which monotone scheme to use.
     &, interp_vertical_search_tol                                      &
                                   !number of levels either side of
                                   ! default level to search.
     &, check_bottom_levels ! used in interpolation code, and is
                            ! the number of levels to check to see
                            ! if the departure point lies inside the
                            ! orography.

      Logical                                                           &
     &  L_high                                                          &
                       ! True, if high order interpolation required.
     &, L_mono                                                          &
                       ! True, if interpolation required to be monotone.
     &, L_conserv      ! True, if interpolation to be monotone and
                       !       conservative.

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, first_flat_level_in

! start of SWW mods
      Integer M,MM




      Real  r_in (1-halo_i:dim_i_in+halo_i,                             &
     &            1-halo_j:dim_j_in+halo_j,                             &
     &            dim_k_in) 

      real r_here_tmp(dim_i_out), r_here_plus_tmp(dim_i_out)
! end of mods
                                                      ! Vertical
                                                      ! co-ordinate
                                                      ! of input data.
      Real  eta_in(dim_k_in)                          ! eta levels of
                                                      ! input data.

      Real                                                              &
     &  r_out (dim_i_out, dim_j_out, dim_k_out)                         &
                                                      ! Vertical
                                                      ! co-ordinate
                                                      ! of output data.
     &, array_r_here(dim_i_out, dim_j_out)

! Arguments with Intent IN/OUT.
      Real                                                              &
     &  coeff_z(dim_i_out, dim_j_out, dim_k_out, -2:3)                  &
     &, coeff_z_lin(dim_i_out, dim_j_out, dim_k_out, 0:1)               &
     &, weight_lambda (dim_i_out, dim_j_out, dim_k_out)                 &
     &, weight_phi (dim_i_out, dim_j_out, dim_k_out)

      INTEGER (KIND=integer32) ::                                       &
     &  k_out (dim_i_out, dim_j_out, dim_k_out)

      INTEGER (KIND=integer32) ::                                       &
     &  i_out (dim_i_out, dim_j_out, dim_k_out)                         &
     &, j_out (dim_i_out, dim_j_out, dim_k_out)
!
      Logical                                                           &
     &  L_continue(dim_i_out, dim_j_out)                                &
     &, L_continue_up(dim_i_out, dim_j_out)                             &
     &, L_continue_down(dim_i_out, dim_j_out)                           &
     &, L_at_lower_lim                                                  &
     &, L_at_upper_lim

      Integer                                                           &
     &  i, j, k, index                                                  &
                       ! Loop indices
     &, count, count_down, count_up                                     &
     &, lower_limit                                                     &
     &, upper_limit

      Real                                                              &
     &  r_here_minus2                                                   &
     &, r_here_minus                                                    &
     &, r_here                                                          &
     &, r_here_plus                                                     &
     &, r_here_plus2                                                    &
     &, r_here_plus3                                                    &
     &, numer_minus                                                     &
     &, numer                                                           &
     &, numer_plus                                                      &
     &, numer_plus2                                                     &
     &, weight_eta                                                      &
                      ! Interpolation weight for eta coordinate.
     &, eta_out                                                         &
                      ! eta coordinate of a departure point.
     &, numer_minus2                                                    &
                      ! Needed to make code for quintic look similar
     &, numer_plus3   ! to code for other orders of interpolation.

! arrays

      Real                                                              &
     &  coeff_a                                                         &
     &, coeff_b                                                         &
     &, coeff_c                                                         &
     &, coeff_d

! Description: COMDECK containing the allowed
!              monotone scheme options
!
      INTEGER                                                           &
     &     triLinear                                                    &
     &,    mono_quasiCubic

      PARAMETER(                                                        &
     &     triLinear       = 1                                          &
     &,    mono_quasiCubic = 2 )

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! No External Routines:

! Functions: None

! ----------------------------------------------------------------------
!  Section 1.   Find model levels just below departure point.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('ETA_VERT_WEIGHTS_LEVEL',zhook_in,zhook_handle)

! calculate horizontal interpolation weights

! start of SWW mods




       M=1
       MM=1


! Find k point.
        If (interp_vertical_search_tol  >   dim_k_out/2 ) Then
! search all levels as time saving is pretty minimal

! Set minimum value to level one.
          Do j = 1, dim_j_out
            Do i = 1, dim_i_out
              k_out(i,j,k) = 1
            End Do
          End Do

! Find level which is just below r_out value
          Do index = 2, dim_k_in-1

            Do j = 1, dim_j_out
              Do i = 1, dim_i_out




            coeff_d = weight_lambda(i,j,k)*                             &
     &                       weight_phi(i,j,k)
            coeff_c = weight_phi(i,j,k) - coeff_d
            coeff_b = weight_lambda(i,j,k) - coeff_d
            coeff_a = 1. - weight_lambda(i,j,k) - coeff_c
                r_here = coeff_a *                                      &
     &                    r_in (i_out(i,j,k),j_out(i,j,k),index)        &
     &                     + coeff_b *                                  &
     &                    r_in (i_out(i,j,k)+1,j_out(i,j,k),index)      &
     &                     + coeff_c *                                  &
     &                    r_in (i_out(i,j,k),j_out(i,j,k)+M,index)      &
     &                     + coeff_d *                                  &
     &                    r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,index)
                If (r_out(i,j,k)  >   r_here ) k_out(i,j,k) = index

              End Do
            End Do



          End Do

! set flag to say that there exists at least one point on each level
! for which k_out is 1 or dim_k-1
          L_at_lower_lim = .true.
          L_at_upper_lim = .true.

        Else

! use search over restricted levels.

! Find level which is just below r_out value, min possible is
! max(1,k- interp_vertical_search_tol), max possible is
! min(dim_k_in, k+interp_vertical_search_tol)-1
          If( k  <=   check_bottom_levels ) Then
            lower_limit = max(2,k - check_bottom_levels)
            upper_limit = min(dim_k_in-1,                               &
     &                         k + check_bottom_levels)
          Else
            lower_limit = max(2, k - interp_vertical_search_tol)
            upper_limit = min(dim_k_in-1,                               &
     &                          k + interp_vertical_search_tol)
          End If

! set flag to say that there exists at least one point on each level
! for which k_out is 1 or dim_k-1 to .true.
          If (lower_limit  >   2) Then
            L_at_lower_lim = .false.
          Else
            L_at_lower_lim = .true.
          End If
          If (upper_limit  <   dim_k_in-1) Then
            L_at_upper_lim = .false.
          Else
            L_at_upper_lim = .true.
          End If

! level 1 Only performs upward search
          If (k  ==  1 ) Then
            count = 1
            index = k + 1
            Do j = 1, dim_j_out
              Do i = 1, dim_i_out
                k_out(i,j,k) = k
                L_continue(i,j) = .true.
              End Do
            End Do

            Do index = 2, upper_limit
              If (count  >   0) Then
                count = 0
                If (index  <   first_flat_level_in ) Then

                  Do j = 1, dim_j_out
                  Do i = 1, dim_i_out





                    If (L_continue(i,j)) Then
! upward search
            coeff_d = weight_lambda(i,j,k)*                             &
     &                       weight_phi(i,j,k)
            coeff_c = weight_phi(i,j,k) - coeff_d
            coeff_b = weight_lambda(i,j,k)  - coeff_d
            coeff_a = 1. - weight_lambda(i,j,k) - coeff_c

                      r_here = coeff_a *                                &
     &               r_in (i_out(i,j,k),j_out(i,j,k),index)             &
     &               + coeff_b *                                        &
     &               r_in (i_out(i,j,k)+1,j_out(i,j,k),index)           &
     &               + coeff_c *                                        &
     &               r_in (i_out(i,j,k),j_out(i,j,k)+M,index)           &
     &               + coeff_d *                                        &
     &               r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,index)

                     r_here_tmp(i) = r_here
                    end if
                  End Do
                  Do i = 1, dim_i_out
                    If (L_continue(i,j)) Then

                      IF &

                          (r_out(i,j,k)  >   r_here_tmp(i) ) &



                          Then
                        k_out(i,j,k) = index
                        count = count + 1
!                     Else
!                       L_continue(i,j) = .false.
                      End If
                    End If

                  End Do
                  End Do



                Else
! r_here is now constant for all points so
                  r_here = r_in (1,MM,index)

                  Do j = 1, dim_j_out
                  Do i = 1, dim_i_out





                    If (L_continue(i,j)) Then
! upward search
                      If (r_out(i,j,k)  >   r_here ) Then
                        k_out(i,j,k) = index
                        count = count + 1
                      Else
                        L_continue(i,j) = .false.
                      End If
                    End If

                  End Do
                End Do



                End If

                If (index  ==  2 .and.                                  &
     &              count  ==  dim_j_out*dim_i_out) then
! all points at level 1 have index 2 or more
                  L_at_lower_lim = .false.
                End If

              End If
            End Do

! AT upper limit Only performs downward search
          Else If (k  >=  dim_k_in ) Then

            Do j = 1, dim_j_out
              Do i = 1, dim_i_out
                k_out(i,j,k) = dim_k_in - 1
                L_continue(i,j) = .true.
              End Do
            End Do

            count = 1

            Do index = dim_k_in-1, lower_limit, -1
              If (count  >   0) Then
                count = 0
                If (index  <   first_flat_level_in ) Then

                Do j = 1, dim_j_out
                Do i = 1, dim_i_out





                  If (L_continue(i,j)) Then
! downward search
            coeff_d = weight_lambda(i,j,k)*                             &
     &                       weight_phi(i,j,k)
            coeff_c = weight_phi(i,j,k) - coeff_d
            coeff_b = weight_lambda(i,j,k) - coeff_d
            coeff_a = 1. - weight_lambda(i,j,k) - coeff_c

                    r_here = coeff_a *                                  &
     &               r_in (i_out(i,j,k),j_out(i,j,k),index)             &
     &               + coeff_b *                                        &
     &               r_in (i_out(i,j,k)+1,j_out(i,j,k),index)           &
     &               + coeff_c *                                        &
     &               r_in (i_out(i,j,k),j_out(i,j,k)+M,index)           &
     &               + coeff_d *                                        &
     &               r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,index)

                     r_here_tmp(i) = r_here
                  End If
                End Do
                Do i = 1, dim_i_out
                  If (L_continue(i,j)) Then

                    IF &

                        (r_out(i,j,k)  <   r_here_tmp(i) ) &



                        Then
                      k_out(i,j,k) = index - 1
                      count = count + 1
                    Else
                      L_continue(i,j) = .false.
                    End If
                  End If

                End Do
              End Do



                Else
! r_here is now constant for all points so
                r_here = r_in (1,MM,index)
                Do j = 1, dim_j_out
                Do i = 1, dim_i_out
                  If (L_continue(i,j)) Then
! downward search
                    If (r_out(i,j,k)  <   r_here ) Then
                      k_out(i,j,k) = index - 1
                      count = count + 1
                    Else
                      L_continue(i,j) = .false.
                    End If
                  End If
                End Do
                End Do

                End If

                If (index  ==  dim_k_in-1                               &
     &            .and. count  ==  dim_j_out*dim_i_out) then
! all points at level dim_k_in-1 have index dim_k_in-2 or less
                  L_at_upper_lim = .false.
                End If

              End If

            End Do


          Else
! Calculate level k value

            count_up = 0
            count_down = 0
            index = k
            If (index  <   first_flat_level_in ) Then


            Do j = 1, dim_j_out
            Do i = 1, dim_i_out





            coeff_d = weight_lambda(i,j,k) * weight_phi(i,j,k)
            coeff_c = weight_phi(i,j,k) - coeff_d
            coeff_b = weight_lambda(i,j,k) - coeff_d
            coeff_a = 1. - weight_lambda(i,j,k) - coeff_c

!J      array_r_here(i,j) = coeff_a *                                   &
            r_here = coeff_a *                                   &
     &               r_in (i_out(i,j,k),j_out(i,j,k),index)             &
     &               + coeff_b *                                        &
     &               r_in (i_out(i,j,k)+1,j_out(i,j,k),index)           &
     &               + coeff_c *                                        &
     &               r_in (i_out(i,j,k),j_out(i,j,k)+M,index)           &
     &               + coeff_d *                                        &
     &               r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,index)
!J          End Do

!J          End Do
!J          Do j = 1, dim_j_out
!J          Do i = 1, dim_i_out
              L_continue_up(i,j) = .false.
              L_continue_down(i,j) = .false.
!J            r_here=array_r_here(i,j)
              If (r_out(i,j,k)  >   r_here ) Then
                k_out(i,j,k) = k
                count_up = count_up + 1
                L_continue_up(i,j) = .true.
              Else
                k_out(i,j,k) = k-1
                count_down = count_down + 1
                L_continue_down(i,j) = .true.
              End If
            End Do
            End Do
            Else
! r_here is now constant for all points so
            r_here = r_in (1,MM,index)
            Do j = 1, dim_j_out
            Do i = 1, dim_i_out
              L_continue_up(i,j) = .false.
              L_continue_down(i,j) = .false.
              If (r_out(i,j,k)  >   r_here ) Then
                k_out(i,j,k) = k
                count_up = count_up + 1
                L_continue_up(i,j) = .true.
              Else
                k_out(i,j,k) = k-1
                count_down = count_down + 1
                L_continue_down(i,j) = .true.
              End If
            End Do
            End Do
            End If

            Do index = k+1, upper_limit
            If (count_up  >   0) Then
              count_up = 0
              If (index  <   first_flat_level_in ) Then
                Do j = 1, dim_j_out
                Do i = 1, dim_i_out
                  If (L_continue_up(i,j)) Then
! upward search
            coeff_d = weight_lambda(i,j,k) * weight_phi(i,j,k)
            coeff_c = weight_phi(i,j,k) - coeff_d
            coeff_b = weight_lambda(i,j,k) - coeff_d
            coeff_a = 1. - weight_lambda(i,j,k) - coeff_c

                    r_here = coeff_a *                                  &
     &               r_in (i_out(i,j,k),j_out(i,j,k),index)             &
     &               + coeff_b *                                        &
     &               r_in (i_out(i,j,k)+1,j_out(i,j,k),index)           &
     &               + coeff_c *                                        &
     &               r_in (i_out(i,j,k),j_out(i,j,k)+M,index)           &
     &               + coeff_d *                                        &
     &               r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,index)
                    If (r_out(i,j,k)  >   r_here ) Then
                      k_out(i,j,k) = index
                      count_up = count_up + 1
                    Else
                      L_continue_up(i,j) = .false.
                    End If
                  End If
                End Do
                End Do
              Else
! r_here is now constant for all points so
                r_here = r_in (1,MM,index)
                Do j = 1, dim_j_out
                Do i = 1, dim_i_out
                  If (L_continue_up(i,j)) Then
! upward search
                    If (r_out(i,j,k)  >   r_here ) Then
                      k_out(i,j,k) = index
                      count_up = count_up + 1
                    Else
                      L_continue_up(i,j) = .false.
                    End If
                  End If
                End Do
              End Do
              End If

              If (index  ==  dim_k_in-1 .and. count_up   /=  0) Then
                L_at_upper_lim = .true.
              End If
            End If

            End Do

            Do index = k-1, lower_limit, -1
            If (count_down  >   0) Then
              count_down = 0
              If (index  <   first_flat_level_in ) Then
                Do j = 1, dim_j_out
                Do i = 1, dim_i_out
                  If (L_continue_down(i,j)) Then
! downward search
            coeff_d = weight_lambda(i,j,k) * weight_phi(i,j,k)
            coeff_c = weight_phi(i,j,k) - coeff_d
            coeff_b = weight_lambda(i,j,k) - coeff_d
            coeff_a = 1. - weight_lambda(i,j,k) - coeff_c


                    r_here = coeff_a *                                  &
     &               r_in (i_out(i,j,k),j_out(i,j,k),index)             &
     &               + coeff_b *                                        &
     &               r_in (i_out(i,j,k)+1,j_out(i,j,k),index)           &
     &               + coeff_c *                                        &
     &               r_in (i_out(i,j,k),j_out(i,j,k)+M,index)           &
     &               + coeff_d *                                        &
     &               r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,index)
                    If (r_out(i,j,k)  <   r_here ) Then
                      k_out(i,j,k) = index - 1
                      count_down = count_down + 1
                    Else
                      L_continue_down(i,j) = .false.
                    End If
                  End If
                End Do
              End Do
              Else
! r_here is now constant for all points so
                r_here = r_in (1,MM,index)
                Do j = 1, dim_j_out
                Do i = 1, dim_i_out
                  If (L_continue_down(i,j)) Then
! downward search
                    If (r_out(i,j,k)  <   r_here ) Then
                      k_out(i,j,k) = index - 1
                      count_down = count_down + 1
                    Else
                      L_continue_down(i,j) = .false.
                    End If
                  End If
                End Do
              End Do
              End If

              If (index  ==  2 .and. count_down   /=  0) Then
                L_at_lower_lim = .true.
              End If

            End If

            End Do

          End If

        End If

! ----------------------------------------------------------------------
! Section 4.   Calculate vertical interpolation weights
! ----------------------------------------------------------------------

        If ( (high_order_scheme  ==  quinticLagrange .or.               &
     &        high_order_scheme  ==  hQuasiCubic_vQuintic .or.          &
     &        high_order_scheme  ==  hCubic_vQuintic) .and. L_High )    &
     &      Then
! quintic interpolation

! Calculate interpolation coefficients

          Do j = 1, dim_j_out
            Do i = 1, dim_i_out

            coeff_d = weight_lambda(i,j,k) * weight_phi(i,j,k)
            coeff_c = weight_phi(i,j,k) - coeff_d
            coeff_b = weight_lambda(i,j,k) - coeff_d
            coeff_a = 1. - weight_lambda(i,j,k) - coeff_c

              If (k_out(i,j,k)  ==  1 .or. k_out(i,j,k)  ==             &
     &              dim_k_in-1) Then
! use linear to interpolate between
! bottom most levels or top most levels.

                r_here = coeff_a                                        &
     &         * r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k))          &
     &                  + coeff_b                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k))        &
     &                  + coeff_c                                       &
     &         * r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k))        &
     &                  + coeff_d                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k))

                r_here_plus = coeff_a                                   &
     &        *  r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+1)        &
     &                  + coeff_b                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+1)      &
     &                  + coeff_c                                       &
     &        *  r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k)+1)      &
     &                  + coeff_d                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k)+1)

! Obtain an eta coordinate for the departure point.
                weight_eta = ( r_out(i,j,k) - r_here ) /                &
     &                       ( r_here_plus - r_here )

                eta_out = eta_in(k_out(i,j,k)) + weight_eta *           &
     &                    ( eta_in(k_out(i,j,k)+1) -                    &
     &                      eta_in(k_out(i,j,k)) )


                coeff_z (i,j,k,-2) = 0.
                coeff_z (i,j,k,-1) = 0.

! cjs Compute coeffs using eta coord system.
! cjs Linear interpolation weight already calculated.
                coeff_z(i,j,k,0) = 1 - weight_eta
                coeff_z(i,j,k,1) = weight_eta

                coeff_z (i,j,k,2) = 0.
                coeff_z (i,j,k,3) = 0.

              Else If (k_out(i,j,k)  ==  2 .or. k_out(i,j,k)  ==        &
     &                 dim_k_in-2) Then

! use cubic interpolation.

! cjs 110200 It is necessary to retain the calculation of r_here and
!     r_here_plus, in order to calculate an eta coordinate for the
!     departure point.
!     This will not be the case if departure point scheme returns
!     an eta coordinate rather than an r coodinate.

                r_here = coeff_a                                        &
     &         * r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k))          &
     &                  + coeff_b                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k))        &
     &                  + coeff_c                                       &
     &         * r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k))        &
     &                  + coeff_d                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k))

                r_here_plus = coeff_a *                                 &
     &           r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+1)        &
     &                  + coeff_b                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+1)      &
     &                  + coeff_c                                       &
     &        *  r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k)+1)      &
     &                  + coeff_d                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k)+1)


!  Obtain an eta coordinate for the departure point.
                weight_eta = ( r_out(i,j,k) - r_here ) /                &
     &                       ( r_here_plus - r_here )

                eta_out = eta_in(k_out(i,j,k)) + weight_eta *           &
     &                    ( eta_in(k_out(i,j,k)+1) -                    &
     &                      eta_in(k_out(i,j,k)) )

! cjs Re-use the r coordinate variables to hold the corresponding
!     eta coordinate values. Interpolation in eta simplifies the
!     algortihm : the knots for vertical interpolation are now the
!     tabulated eta levels, rather than the linearly interpolated r
!     values.

                r_here_minus = eta_in(k_out(i,j,k) - 1)
                r_here       = eta_in(k_out(i,j,k))
                r_here_plus  = eta_in(k_out(i,j,k) + 1)
                r_here_plus2 = eta_in(k_out(i,j,k) + 2)


! Compute coordinate differences in eta.

                numer_minus = eta_out - r_here_minus
                numer       = eta_out - r_here
                numer_plus  = eta_out - r_here_plus
                numer_plus2 = eta_out - r_here_plus2

                coeff_z(i,j,k,-2)= 0.

                coeff_z(i,j,k,-1)= (numer*numer_plus * numer_plus2 ) /  &
     &                            ( (r_here_minus - r_here) *           &
     &                              (r_here_minus - r_here_plus )*      &
     &                              (r_here_minus - r_here_plus2 ) )

                coeff_z(i,j,k,0) = (numer_minus * numer_plus *          &
     &                             numer_plus2 ) /                      &
     &                            ( (r_here - r_here_minus) *           &
     &                              (r_here - r_here_plus )*            &
     &                              (r_here - r_here_plus2 ) )

                coeff_z(i,j,k,1) = (numer_minus * numer *               &
     &                             numer_plus2 ) /                      &
     &                            ( (r_here_plus - r_here_minus) *      &
     &                              (r_here_plus - r_here )*            &
     &                              (r_here_plus - r_here_plus2 ) )

                coeff_z(i,j,k,2) = (numer_minus * numer *               &
     &                             numer_plus ) /                       &
     &                            ( (r_here_plus2 - r_here_minus) *     &
     &                              (r_here_plus2 - r_here )*           &
     &                              (r_here_plus2 - r_here_plus ) )

                coeff_z(i,j,k,3)= 0.

              Else
! quintic interpolation


                r_here = coeff_a                                        &
     &         * r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k))          &
     &                  + coeff_b                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k))        &
     &                  + coeff_c                                       &
     &         * r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k))        &
     &                  + coeff_d                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k))

                r_here_plus = coeff_a *                                 &
     &           r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+1)        &
     &                  + coeff_b                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+1)      &
     &                  + coeff_c                                       &
     &        *  r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k)+1)      &
     &                  + coeff_d                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k)+1)


!  Obtain an eta coordinate for the departure point.
                weight_eta = ( r_out(i,j,k) - r_here ) /                &
     &                       ( r_here_plus - r_here )

                eta_out = eta_in(k_out(i,j,k)) + weight_eta *           &
     &                    ( eta_in(k_out(i,j,k)+1) -                    &
     &                      eta_in(k_out(i,j,k)) )

! Re-use the r coordinate variables to hold the corresponding
!     eta coordinate values.

                r_here_minus2 = eta_in(k_out(i,j,k) - 2)
                r_here_minus  = eta_in(k_out(i,j,k) - 1)
                r_here        = eta_in(k_out(i,j,k))
                r_here_plus   = eta_in(k_out(i,j,k) + 1)
                r_here_plus2  = eta_in(k_out(i,j,k) + 2)
                r_here_plus3  = eta_in(k_out(i,j,k) + 3)

!  Compute coordinate differences in eta.

                numer_minus2 = eta_out - r_here_minus2
                numer_minus  = eta_out - r_here_minus
                numer        = eta_out - r_here
                numer_plus   = eta_out - r_here_plus
                numer_plus2  = eta_out - r_here_plus2
                numer_plus3  = eta_out - r_here_plus3

                coeff_z(i,j,k,-2)= ( numer_minus * numer * numer_plus * &
     &                               numer_plus2 * numer_plus3 ) /      &
     &                             ( (r_here_minus2 - r_here_minus) *   &
     &                               (r_here_minus2 - r_here ) *        &
     &                               (r_here_minus2 - r_here_plus ) *   &
     &                               (r_here_minus2 - r_here_plus2 ) *  &
     &                               (r_here_minus2 - r_here_plus3 ) )

                coeff_z(i,j,k,-1)= ( numer_minus2 * numer * numer_plus *&
     &                               numer_plus2 * numer_plus3 ) /      &
     &                             ( (r_here_minus - r_here_minus2) *   &
     &                               (r_here_minus - r_here )*          &
     &                               (r_here_minus - r_here_plus )*     &
     &                               (r_here_minus - r_here_plus2 )*    &
     &                               (r_here_minus - r_here_plus3 ) )

                coeff_z(i,j,k,0) = ( numer_minus2 * numer_minus *       &
     &                               numer_plus * numer_plus2 *         &
     &                               numer_plus3 )/                     &
     &                             ( (r_here - r_here_minus2) *         &
     &                               (r_here - r_here_minus )*          &
     &                               (r_here - r_here_plus )*           &
     &                               (r_here - r_here_plus2 )*          &
     &                               (r_here - r_here_plus3 ) )

                coeff_z(i,j,k,1) = ( numer_minus2 * numer_minus *       &
     &                               numer * numer_plus2 *              &
     &                               numer_plus3 ) /                    &
     &                             ( (r_here_plus - r_here_minus2) *    &
     &                               (r_here_plus - r_here_minus )*     &
     &                               (r_here_plus - r_here )*           &
     &                               (r_here_plus - r_here_plus2 )*     &
     &                               (r_here_plus - r_here_plus3 ) )

                coeff_z(i,j,k,2) = ( numer_minus2 * numer_minus *       &
     &                               numer * numer_plus *               &
     &                               numer_plus3 ) /                    &
     &                             ( (r_here_plus2 - r_here_minus2) *   &
     &                               (r_here_plus2 - r_here_minus )*    &
     &                               (r_here_plus2 - r_here )*          &
     &                               (r_here_plus2 - r_here_plus )*     &
     &                               (r_here_plus2 - r_here_plus3 ) )

                coeff_z(i,j,k,3)= ( numer_minus2 * numer_minus *        &
     &                              numer * numer_plus *                &
     &                              numer_plus2 ) /                     &
     &                             ( (r_here_plus3 - r_here_minus2) *   &
     &                               (r_here_plus3 - r_here_minus )*    &
     &                               (r_here_plus3 - r_here )*          &
     &                               (r_here_plus3 - r_here_plus )*     &
     &                               (r_here_plus3 - r_here_plus2 ) )
              End If

            End Do
          End Do

        Else If ( ( ( high_order_scheme  ==  cubicLagrange              &
     &          .or.  high_order_scheme  ==  ECMWF_quasiCubic           &
     &          .or.  high_order_scheme  ==  ECMWF_mono_quasiCubic )    &
     &          .and. L_High )                                          &
     &          .or. ( monotone_scheme  ==  mono_quasiCubic             &
     &         .and. L_Mono ) ) Then

! cubic interpolation

! Calculate interpolation coefficients

          If (L_at_upper_lim .or. L_at_lower_lim ) Then
            Do j = 1, dim_j_out
            Do i = 1, dim_i_out

            coeff_d = weight_lambda(i,j,k) * weight_phi(i,j,k)
            coeff_c = weight_phi(i,j,k) - coeff_d
            coeff_b = weight_lambda(i,j,k) - coeff_d
            coeff_a = 1. - weight_lambda(i,j,k) - coeff_c

              If (k_out(i,j,k)  ==  1 .or. k_out(i,j,k)  ==             &
     &              dim_k_in-1) Then
! use linear to interpolate between
! bottom most levels or top most levels.

                r_here = coeff_a                                        &
     &         * r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k))          &
     &                  + coeff_b                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k))        &
     &                  + coeff_c                                       &
     &         * r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k))        &
     &                  + coeff_d                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k))

                r_here_plus = coeff_a                                   &
     &        *  r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+1)        &
     &                  + coeff_b                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+1)      &
     &                  + coeff_c                                       &
     &        *  r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k)+1)      &
     &                  + coeff_d                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k)+1)

!  Obtain an eta coordinate for the departure point.
                weight_eta = ( r_out(i,j,k) - r_here ) /                &
     &                       ( r_here_plus - r_here )

                eta_out = eta_in(k_out(i,j,k)) + weight_eta *           &
     &                    ( eta_in(k_out(i,j,k)+1) -                    &
     &                      eta_in(k_out(i,j,k)) )

                coeff_z (i,j,k,-1) = 0.

!  Compute coeffs using eta coord system.

!  Linear interpolation weight already calculated.
                coeff_z(i,j,k,0) = 1 - weight_eta
                coeff_z(i,j,k,1) = weight_eta

                coeff_z (i,j,k,2) = 0.

              Else ! If not near boundary.
! use cubic interpolation.

                r_here = coeff_a                                        &
     &         * r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k))          &
     &                  + coeff_b                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k))        &
     &                  + coeff_c                                       &
     &         * r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k))        &
     &                  + coeff_d                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k))

                r_here_plus = coeff_a *                                 &
     &           r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+1)        &
     &                  + coeff_b                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+1)      &
     &                  + coeff_c                                       &
     &        *  r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k)+1)      &
     &                  + coeff_d                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k)+1)

!  Obtain an eta coordinate for the departure point.
                weight_eta = ( r_out(i,j,k) - r_here ) /                &
     &                       ( r_here_plus - r_here )

                eta_out = eta_in(k_out(i,j,k)) + weight_eta *           &
     &                    ( eta_in(k_out(i,j,k)+1) -                    &
     &                      eta_in(k_out(i,j,k)) )

! cjs Re-use the r coordinate variables to hold the corresponding
!     eta coordinate values. Interpolation in eta simplifies the
!     algortihm : the knots for vertical interpolation are now the
!     tabulated eta levels, rather than the linearly interpolated r
!     values.

                r_here_minus = eta_in(k_out(i,j,k) - 1)
                r_here       = eta_in(k_out(i,j,k))
                r_here_plus  = eta_in(k_out(i,j,k) + 1)
                r_here_plus2 = eta_in(k_out(i,j,k) + 2)


! cjs Compute coordinate differences in eta.

                numer_minus = eta_out - r_here_minus
                numer       = eta_out - r_here
                numer_plus  = eta_out - r_here_plus
                numer_plus2 = eta_out - r_here_plus2


                coeff_z(i,j,k,-1)= (numer * numer_plus *numer_plus2) /  &
     &                            ( (r_here_minus - r_here) *           &
     &                              (r_here_minus - r_here_plus )*      &
     &                              (r_here_minus - r_here_plus2 ) )

                coeff_z(i,j,k,0) = (numer_minus * numer_plus *          &
     &                             numer_plus2 ) /                      &
     &                            ( (r_here - r_here_minus) *           &
     &                              (r_here - r_here_plus )*            &
     &                              (r_here - r_here_plus2 ) )

                coeff_z(i,j,k,1) = (numer_minus * numer *               &
     &                             numer_plus2 ) /                      &
     &                            ( (r_here_plus - r_here_minus) *      &
     &                              (r_here_plus - r_here )*            &
     &                              (r_here_plus - r_here_plus2 ) )

                coeff_z(i,j,k,2) = (numer_minus * numer *               &
     &                             numer_plus ) /                       &
     &                            ( (r_here_plus2 - r_here_minus) *     &
     &                              (r_here_plus2 - r_here )*           &
     &                              (r_here_plus2 - r_here_plus ) )

              End If ! Check on being close to boundary.

            End Do
          End Do

          Else  ! If not at upper/lower limit of search.

            Do j = 1, dim_j_out
              Do i = 1, dim_i_out
! use cubic interpolation.
            coeff_d = weight_lambda(i,j,k) * weight_phi(i,j,k)
            coeff_c = weight_phi(i,j,k) - coeff_d
            coeff_b = weight_lambda(i,j,k) - coeff_d
            coeff_a = 1. - weight_lambda(i,j,k) - coeff_c

                r_here = coeff_a                                        &
     &         * r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k))          &
     &                  + coeff_b                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k))        &
     &                  + coeff_c                                       &
     &         * r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k))        &
     &                  + coeff_d                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k))

                r_here_plus = coeff_a *                                 &
     &           r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+1)        &
     &                  + coeff_b                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+1)      &
     &                  + coeff_c                                       &
     &        *  r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k)+1)      &
     &                  + coeff_d                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k)+1)

!  Obtain an eta coordinate for the departure point.
                weight_eta = ( r_out(i,j,k) - r_here ) /                &
     &                       ( r_here_plus - r_here )

                eta_out = eta_in(k_out(i,j,k)) + weight_eta *           &
     &                    ( eta_in(k_out(i,j,k)+1) -                    &
     &                      eta_in(k_out(i,j,k)) )

! cjs Re-use the r coordinate variables to hold the corresponding
!     eta coordinate values. Interpolation in eta simplifies the
!     algortihm : the knots for vertical interpolation are now the
!     tabulated eta levels, rather than the linearly interpolated r
!     values.

                r_here_minus = eta_in(k_out(i,j,k) - 1)
                r_here       = eta_in(k_out(i,j,k))
                r_here_plus  = eta_in(k_out(i,j,k) + 1)
                r_here_plus2 = eta_in(k_out(i,j,k) + 2)


! cjs Compute coordinate differences in eta.

                numer_minus = eta_out - r_here_minus
                numer       = eta_out - r_here
                numer_plus  = eta_out - r_here_plus
                numer_plus2 = eta_out - r_here_plus2


                coeff_z(i,j,k,-1)= (numer * numer_plus *numer_plus2) /  &
     &                            ( (r_here_minus - r_here) *           &
     &                              (r_here_minus - r_here_plus )*      &
     &                              (r_here_minus - r_here_plus2 ) )

                coeff_z(i,j,k,0) = (numer_minus * numer_plus *          &
     &                             numer_plus2 ) /                      &
     &                            ( (r_here - r_here_minus) *           &
     &                              (r_here - r_here_plus )*            &
     &                              (r_here - r_here_plus2 ) )

                coeff_z(i,j,k,1) = (numer_minus * numer *               &
     &                             numer_plus2 ) /                      &
     &                            ( (r_here_plus - r_here_minus) *      &
     &                              (r_here_plus - r_here )*            &
     &                              (r_here_plus - r_here_plus2 ) )

                coeff_z(i,j,k,2) = (numer_minus * numer *               &
     &                             numer_plus ) /                       &
     &                            ( (r_here_plus2 - r_here_minus) *     &
     &                              (r_here_plus2 - r_here )*           &
     &                              (r_here_plus2 - r_here_plus ) )
              End Do
            End Do

          End If

        End If

      If ( ( high_order_scheme  ==  hCubic_vLin .and. L_High )          &
     &        .or. ( monotone_scheme  ==  triLinear .and. L_Mono )      &
     &        .or. L_conserv ) Then

! linear interpolation
          Do j = 1, dim_j_out
            Do i = 1, dim_i_out

! use linear to interpolate between
! bottom most levels or top most levels.

            coeff_d = weight_lambda(i,j,k) * weight_phi(i,j,k)
            coeff_c = weight_phi(i,j,k) - coeff_d
            coeff_b = weight_lambda(i,j,k) - coeff_d
            coeff_a = 1. - weight_lambda(i,j,k) - coeff_c

              r_here = coeff_a                                          &
     &         * r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k))          &
     &                  + coeff_b                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k))        &
     &                  + coeff_c                                       &
     &         * r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k))        &
     &                  + coeff_d                                       &
     &         * r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k))

              r_here_plus = coeff_a                                     &
     &        *  r_in (i_out(i,j,k),j_out(i,j,k),k_out(i,j,k)+1)        &
     &                  + coeff_b                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k),k_out(i,j,k)+1)      &
     &                  + coeff_c                                       &
     &        *  r_in (i_out(i,j,k),j_out(i,j,k)+M,k_out(i,j,k)+1)      &
     &                  + coeff_d                                       &
     &        *  r_in (i_out(i,j,k)+1,j_out(i,j,k)+M,k_out(i,j,k)+1)

                r_here_tmp(i) = r_here
                r_here_plus_tmp(i) = r_here_plus
              end do 
            Do i = 1, dim_i_out
                weight_eta = ( r_out(i,j,k) - r_here_tmp(i) ) /                &
     &                       ( r_here_plus_tmp(i) - r_here_tmp(i) )

              If (weight_eta < 0.0) Then
                coeff_z_lin (i,j,k,0) = 1.0 
                coeff_z_lin (i,j,k,1) = 0.0
              Else If(weight_eta > 1.0) Then
                coeff_z_lin (i,j,k,0) = 0.0 
                coeff_z_lin (i,j,k,1) = 1.0
              Else
              coeff_z_lin (i,j,k,0) = 1 - weight_eta
              coeff_z_lin (i,j,k,1) = weight_eta
              end if

            End Do
          End Do

        End If

! End of routine.
        IF (lhook) CALL dr_hook('ETA_VERT_WEIGHTS_LEVEL',zhook_out,zhook_handle)
        return
      END SUBROUTINE eta_vert_weights_level

