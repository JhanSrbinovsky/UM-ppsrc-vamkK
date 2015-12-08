! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine eta_vert_weights_e
      Subroutine eta_vert_weights_e (                                   &
     &                        eta_in,                                   &
     &                        r_in,                                     &
     &                        check_bottom_levels,                      &
     &                        interp_vertical_search_tol,               &
     &                        first_flat_level_in,                      &
     &                        dim_i_in, dim_j_in, dim_k_in,             &
     &                        dim_k_out, dim_e_out,                     &
     &                        high_order_scheme, monotone_scheme,       &
     &                        model_domain, L_high, L_mono,             &
     &                        L_conserv,                                &
     &                        halo_i, halo_j,                           &
     &                        r_out, k_level,                           &
     &                        coeff_z, coeff_z_lin,                     &
     &                        k_out, i_out, j_out,                      &
     &                        weight_lambda, weight_phi )

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
! Purpose: Calculates vertical interpolation weights and level which
!          is below the desired point.
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
     &, dim_k_out                                                       &
                    ! Dimension of departure point Data in k direction.
     &, dim_e_out                                                       &
                    ! Number of points to do.
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
     &, first_flat_level_in                                             &
                            ! level at which model levels go flat
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
     &  model_domain     ! holds integer code for model domain

      Real                                                              &
     &  r_in (1-halo_i:dim_i_in+halo_i,                                 &
     &        1-halo_j:dim_j_in+halo_j, dim_k_in)                       &
                                                      ! Vertical
                                                      ! co-ordinate
                                                      ! of input data.
     &, eta_in(dim_k_in)                              ! eta levels of
                                                      ! input data.

      Real                                                              &
     &  r_out (dim_e_out)       ! Vertical co-ordinate
                                ! of output data.

      Integer                                                           &
     &  k_level (dim_e_out)     ! vertical level of arrival point

! Arguments with Intent IN/OUT.
      Real                                                              &
     &  coeff_z(dim_e_out, -2:3)                                        &
     &, coeff_z_lin(dim_e_out, 0:1)                                     &
     &, weight_lambda (dim_e_out)                                       &
     &, weight_phi (dim_e_out)

      INTEGER (KIND=integer32) ::                                       &
     &  k_out (dim_e_out)

      INTEGER (KIND=integer32) ::                                       &
     &  i_out (dim_e_out)                                               &
     &, j_out (dim_e_out)

! Local Variables.

! scalars

      Logical                                                           &
     &  L_continue

      Integer                                                           &
     &  i, index                                                        &
                 ! Loop indices
     &, lower_limit, upper_limit

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
     &  coeff_a(dim_e_out)                                              &
     &, coeff_b(dim_e_out)                                              &
     &, coeff_c(dim_e_out)                                              &
     &, coeff_d(dim_e_out)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Description: COMDECK containing the allowed
!              monotone scheme options
!
      INTEGER                                                           &
     &     triLinear                                                    &
     &,    mono_quasiCubic

      PARAMETER(                                                        &
     &     triLinear       = 1                                          &
     &,    mono_quasiCubic = 2 )

! No External Routines:

! Functions: None

! ----------------------------------------------------------------------
!  Section 1.   Find model levels just below departure point.
! ----------------------------------------------------------------------

! calculate horizontal interpolation weights
      IF (lhook) CALL dr_hook('ETA_VERT_WEIGHTS_E',zhook_in,zhook_handle)
      
!$OMP  PARALLEL DO SCHEDULE(STATIC) SHARED(coeff_d, coeff_c, coeff_b,   &
!$OMP& coeff_a, weight_lambda, weight_phi) PRIVATE(i)                   
      Do i = 1, dim_e_out
        coeff_d(i) = weight_lambda(i) * weight_phi(i)
        coeff_c(i) = weight_phi(i) - coeff_d(i)
        coeff_b(i) = weight_lambda(i) - coeff_d(i)
        coeff_a(i) = 1. - weight_lambda(i) - coeff_c(i)
      End Do
!$OMP END PARALLEL DO 

! Find k point.
      If (interp_vertical_search_tol  >   dim_k_out/2 ) Then
! search all levels as time saving is pretty minimal


! Set minimum value to level one.

        Do i = 1, dim_e_out
          k_out(i) = 1
        End Do
        
! Find level which is just below r_out value

        Do index = 2, dim_k_in-1
          Do i = 1, dim_e_out
            r_here = coeff_a(i) *                                       &
     &                   r_in (i_out(i),j_out(i),index)                 &
     &                 + coeff_b(i) *                                   &
     &                   r_in (i_out(i)+1,j_out(i),index)               &
     &                 + coeff_c(i) *                                   &
     &                   r_in (i_out(i),j_out(i)+1,index)               &
     &                 + coeff_d(i) *                                   &
     &                   r_in (i_out(i)+1,j_out(i)+1,index)
            If (r_out(i)  >   r_here ) k_out(i) = index
          End Do
        End Do

      Else

! use search over restricted levels.
! Find level which is just below r_out value, min possible is
! max(1,k- interp_vertical_search_tol), max possible is
! min(dim_k_in, k+interp_vertical_search_tol)-1

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(lower_limit,&
!$OMP& upper_limit, index, l_continue, r_here, i)
        Do i = 1, dim_e_out

          If (k_level(i)  <=  check_bottom_levels) Then
            lower_limit = max(1,k_level(i) - check_bottom_levels)
            upper_limit = min(dim_k_in,                                 &
     &                            k_level(i) + check_bottom_levels)
          Else
            lower_limit = max(1,k_level(i) -                            &
     &                            interp_vertical_search_tol)
            upper_limit = min(dim_k_in,                                 &
     &                          k_level(i) + interp_vertical_search_tol)
          End If

! level 1 Only performs upward search
          If (k_level(i)  ==  1) Then
            k_out(i) = k_level(i)
            index = k_level(i)
            L_continue = .true.
! upward search
            Do while (index  <   upper_limit-1 .and.                    &
     &                L_continue)
              index = index + 1
              If (index >= first_flat_level_in) Then
                r_here=r_in (i_out(i),j_out(i),index)
              Else
                r_here = coeff_a(i) *                                   &
     &                   r_in (i_out(i),j_out(i),index)                 &
     &                 + coeff_b(i) *                                   &
     &                   r_in (i_out(i)+1,j_out(i),index)               &
     &                 + coeff_c(i) *                                   &
     &                   r_in (i_out(i),j_out(i)+1,index)               &
     &                 + coeff_d(i) *                                   &
     &                   r_in (i_out(i)+1,j_out(i)+1,index)
              End If  
              If (r_out(i)  >   r_here ) Then
                k_out(i) = index
              Else
                L_continue = .false.
              End If
            End Do
          Else If (k_level(i)  >=  dim_k_in   ) Then
! AT upper limit Only performs downward search
            k_out(i) = upper_limit - 1
            index = upper_limit - 1
            L_continue = .true.
! downward search
            Do while (index  >   lower_limit+1 .and.                    &
     &                L_continue)
              index = index - 1
              If (index >= first_flat_level_in) Then
                r_here=r_in (i_out(i),j_out(i),index)
              Else
                r_here = coeff_a(i) *                                   &
     &                   r_in (i_out(i),j_out(i),index)                 &
     &                 + coeff_b(i) *                                   &
     &                   r_in (i_out(i)+1,j_out(i),index)               &
     &                 + coeff_c(i) *                                   &
     &                   r_in (i_out(i),j_out(i)+1,index)               &
     &                 + coeff_d(i) *                                   &
     &                   r_in (i_out(i)+1,j_out(i)+1,index)
              End If
              If (r_out(i)  <   r_here ) Then
                k_out(i) = index - 1
              Else
                L_continue = .false.
              End If
            End Do

          Else
! Calculate level k value
            If (k_level(i) >= first_flat_level_in) Then
              r_here=r_in (i_out(i),j_out(i),k_level(i))
            Else
              r_here = coeff_a(i) *                                     &
     &                 r_in (i_out(i),j_out(i),k_level(i))              &
     &               + coeff_b(i) *                                     &
     &                 r_in (i_out(i)+1,j_out(i),k_level(i))            &
     &               + coeff_c(i) *                                     &
     &                 r_in (i_out(i),j_out(i)+1,k_level(i))            &
     &               + coeff_d(i) *                                     &
     &                 r_in (i_out(i)+1,j_out(i)+1,k_level(i))
            End If
            L_continue = .true.
            index = k_level(i)
            If (r_out(i)  >   r_here ) Then
              k_out(i) = k_level(i)
! upward search
              Do while (index  <   upper_limit-1 .and.                  &
     &                  L_continue)
                index = index + 1
                If (index >= first_flat_level_in) Then
                  r_here=r_in (i_out(i),j_out(i),index)
                Else
                  r_here = coeff_a(i) *                                 &
     &                   r_in (i_out(i),j_out(i),index)                 &
     &                   + coeff_b(i) *                                 &
     &                   r_in (i_out(i)+1,j_out(i),index)               &
     &                   + coeff_c(i) *                                 &
     &                   r_in (i_out(i),j_out(i)+1,index)               &
     &                   + coeff_d(i) *                                 &
     &                   r_in (i_out(i)+1,j_out(i)+1,index)
                End If
                If (r_out(i)  >   r_here ) Then
                  k_out(i) = index
                Else
                  L_continue = .false.
                End If
              End Do

            Else
              k_out(i) = k_level(i)-1
! downward search
              Do while (index  >   lower_limit+1 .and.                  &
     &                  L_continue)
                index = index - 1
                If (index >= first_flat_level_in) Then
                  r_here=r_in (i_out(i),j_out(i),index)
                Else
                  r_here = coeff_a(i) *                                 &
     &                   r_in (i_out(i),j_out(i),index)                 &
     &                   + coeff_b(i) *                                 &
     &                   r_in (i_out(i)+1,j_out(i),index)               &
     &                   + coeff_c(i) *                                 &
     &                   r_in (i_out(i),j_out(i)+1,index)               &
     &                   + coeff_d(i) *                                 &
     &                   r_in (i_out(i)+1,j_out(i)+1,index)
                End If
                If (r_out(i)  <   r_here ) Then
                  k_out(i) = index - 1
                Else
                  L_continue = .false.
                End If
              End Do

            End If

          End If

        End Do       !i
!$OMP END PARALLEL DO 

      End If     ! restricted

! ----------------------------------------------------------------------
! Section 4.   Calculate vertical interpolation weights
! ----------------------------------------------------------------------

        If ( ( high_order_scheme  ==  quinticLagrange .or.              &
     &         high_order_scheme  ==  hQuasiCubic_vQuintic .or.         &
     &         high_order_scheme  ==  hCubic_vQuintic ) .and. L_High )  &
     &      Then
! quintic interpolation

! Calculate interpolation coefficients

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i, r_here,  &
!$OMP& r_here_plus, weight_eta, eta_out, numer_minus, numer, numer_plus,&
!$OMP& numer_plus2, r_here_plus2, r_here_minus, numer_plus3,            &
!$OMP& r_here_plus3, r_here_minus2, numer_minus2)
        Do i = 1, dim_e_out

          If (k_out(i)  ==  1 .or. k_out(i)  ==  dim_k_in-1) Then
! use linear to interpolate between
! bottom most levels or top most levels.

            r_here = coeff_a(i)                                         &
     &               * r_in (i_out(i),j_out(i),k_out(i))                &
     &             + coeff_b(i)                                         &
     &               * r_in (i_out(i)+1,j_out(i),k_out(i))              &
     &             + coeff_c(i)                                         &
     &               * r_in (i_out(i),j_out(i)+1,k_out(i))              &
     &             + coeff_d(i)                                         &
     &               * r_in (i_out(i)+1,j_out(i)+1,k_out(i))

            r_here_plus = coeff_a(i)                                    &
     &                    *  r_in (i_out(i),j_out(i),k_out(i)+1)        &
     &                  + coeff_b(i)                                    &
     &                    *  r_in (i_out(i)+1,j_out(i),k_out(i)+1)      &
     &                  + coeff_c(i)                                    &
     &                    *  r_in (i_out(i),j_out(i)+1,k_out(i)+1)      &
     &                  + coeff_d(i)                                    &
     &                    *  r_in (i_out(i)+1,j_out(i)+1,k_out(i)+1)

! Obtain an eta coordinate for the departure point.
            weight_eta = ( r_out(i) - r_here ) /                        &
     &                       ( r_here_plus - r_here )

            coeff_z (i,-2) = 0.
            coeff_z (i,-1) = 0.

            coeff_z (i,0) = 1 - weight_eta
            coeff_z (i,1) = weight_eta

            coeff_z (i,2) = 0.
            coeff_z (i,3) = 0.

          Else If (k_out(i)  ==  2 .or. k_out(i)  ==                    &
     &              dim_k_in-2) Then

! use cubic interpolation.

! cjs 110200 It is necessary to retain the calculation of r_here and
!     r_here_plus, in order to calculate an eta coordinate for the
!     departure point.
!     This will not be the case if departure point scheme returns
!     an eta coordinate rather than an r coodinate.

            r_here = coeff_a(i)                                         &
     &         * r_in (i_out(i),j_out(i),k_out(i))                      &
     &                  + coeff_b(i)                                    &
     &         * r_in (i_out(i)+1,j_out(i),k_out(i))                    &
     &                  + coeff_c(i)                                    &
     &         * r_in (i_out(i),j_out(i)+1,k_out(i))                    &
     &                  + coeff_d(i)                                    &
     &         * r_in (i_out(i)+1,j_out(i)+1,k_out(i))

            r_here_plus = coeff_a(i) *                                  &
     &           r_in (i_out(i),j_out(i),k_out(i)+1)                    &
     &                  + coeff_b(i)                                    &
     &        *  r_in (i_out(i)+1,j_out(i),k_out(i)+1)                  &
     &                  + coeff_c(i)                                    &
     &        *  r_in (i_out(i),j_out(i)+1,k_out(i)+1)                  &
     &                  + coeff_d(i)                                    &
     &        *  r_in (i_out(i)+1,j_out(i)+1,k_out(i)+1)


!  Obtain an eta coordinate for the departure point.
            weight_eta = ( r_out(i) - r_here ) /                        &
     &                   ( r_here_plus - r_here )

! cjs Eta_out calculated by assuming eta varies linearly between
!     neighbouring grid levels.
            eta_out = eta_in(k_out(i)) + weight_eta *                   &
     &                ( eta_in(k_out(i)+1) - eta_in(k_out(i)) )

! cjs Re-use the r coordinate variables to hold the corresponding
!     eta coordinate values. Interpolation in eta simplifies the
!     algortihm : the knots for vertical interpolation are now the
!     tabulated eta levels, rather than the linearly interpolated r
!     values.

            r_here_minus = eta_in(k_out(i) - 1)
            r_here       = eta_in(k_out(i))
            r_here_plus  = eta_in(k_out(i) + 1)
            r_here_plus2 = eta_in(k_out(i) + 2)

            numer_minus = eta_out - r_here_minus
            numer       = eta_out - r_here
            numer_plus  = eta_out - r_here_plus
            numer_plus2 = eta_out - r_here_plus2

            coeff_z(i,-2)= 0.
            coeff_z(i,-1)= (numer * numer_plus * numer_plus2 ) /        &
     &                            ( (r_here_minus - r_here) *           &
     &                              (r_here_minus - r_here_plus )*      &
     &                              (r_here_minus - r_here_plus2 ) )
            coeff_z(i,0) = (numer_minus * numer_plus *                  &
     &                             numer_plus2 ) /                      &
     &                            ( (r_here - r_here_minus) *           &
     &                              (r_here - r_here_plus )*            &
     &                              (r_here - r_here_plus2 ) )
            coeff_z(i,1) = (numer_minus * numer *                       &
     &                             numer_plus2 ) /                      &
     &                            ( (r_here_plus - r_here_minus) *      &
     &                              (r_here_plus - r_here )*            &
     &                              (r_here_plus - r_here_plus2 ) )
            coeff_z(i,2) = (numer_minus * numer *                       &
     &                             numer_plus ) /                       &
     &                            ( (r_here_plus2 - r_here_minus) *     &
     &                              (r_here_plus2 - r_here )*           &
     &                              (r_here_plus2 - r_here_plus ) )
            coeff_z(i,3)= 0.

          Else
! quintic interpolation

            r_here = coeff_a(i)                                         &
     &         * r_in (i_out(i),j_out(i),k_out(i))                      &
     &                  + coeff_b(i)                                    &
     &         * r_in (i_out(i)+1,j_out(i),k_out(i))                    &
     &                  + coeff_c(i)                                    &
     &         * r_in (i_out(i),j_out(i)+1,k_out(i))                    &
     &                  + coeff_d(i)                                    &
     &         * r_in (i_out(i)+1,j_out(i)+1,k_out(i))

            r_here_plus = coeff_a(i) *                                  &
     &           r_in (i_out(i),j_out(i),k_out(i)+1)                    &
     &                  + coeff_b(i)                                    &
     &        *  r_in (i_out(i)+1,j_out(i),k_out(i)+1)                  &
     &                  + coeff_c(i)                                    &
     &        *  r_in (i_out(i),j_out(i)+1,k_out(i)+1)                  &
     &                  + coeff_d(i)                                    &
     &        *  r_in (i_out(i)+1,j_out(i)+1,k_out(i)+1)

!  Obtain an eta coordinate for the departure point.
            weight_eta = ( r_out(i) - r_here ) /                        &
     &                   ( r_here_plus - r_here )

            eta_out = eta_in(k_out(i)) + weight_eta *                   &
     &                ( eta_in(k_out(i)+1) - eta_in(k_out(i)) )

! cjs Re-use the r coordinate variables to hold the corresponding
!     eta coordinate values.

            r_here_minus2 = eta_in(k_out(i) - 2)
            r_here_minus  = eta_in(k_out(i) - 1)
            r_here        = eta_in(k_out(i))
            r_here_plus   = eta_in(k_out(i) + 1)
            r_here_plus2  = eta_in(k_out(i) + 2)
            r_here_plus3  = eta_in(k_out(i) + 3)

! cjs Compute coordinate differences in eta.

            numer_minus2 = eta_out - r_here_minus2
            numer_minus  = eta_out - r_here_minus
            numer        = eta_out - r_here
            numer_plus   = eta_out - r_here_plus
            numer_plus2  = eta_out - r_here_plus2
            numer_plus3  = eta_out - r_here_plus3

            coeff_z(i,-2)= ( numer_minus * numer * numer_plus *         &
     &                               numer_plus2 * numer_plus3 ) /      &
     &                             ( (r_here_minus2 - r_here_minus) *   &
     &                               (r_here_minus2 - r_here ) *        &
     &                               (r_here_minus2 - r_here_plus ) *   &
     &                               (r_here_minus2 - r_here_plus2 ) *  &
     &                               (r_here_minus2 - r_here_plus3 ) )

            coeff_z(i,-1)= ( numer_minus2 * numer * numer_plus *        &
     &                               numer_plus2 * numer_plus3 ) /      &
     &                             ( (r_here_minus - r_here_minus2) *   &
     &                               (r_here_minus - r_here )*          &
     &                               (r_here_minus - r_here_plus )*     &
     &                               (r_here_minus - r_here_plus2 )*    &
     &                               (r_here_minus - r_here_plus3 ) )

            coeff_z(i,0) = ( numer_minus2 * numer_minus *               &
     &                               numer_plus * numer_plus2 *         &
     &                               numer_plus3 )/                     &
     &                             ( (r_here - r_here_minus2) *         &
     &                               (r_here - r_here_minus )*          &
     &                               (r_here - r_here_plus )*           &
     &                               (r_here - r_here_plus2 )*          &
     &                               (r_here - r_here_plus3 ) )

            coeff_z(i,1) = ( numer_minus2 * numer_minus *               &
     &                               numer * numer_plus2 *              &
     &                               numer_plus3 ) /                    &
     &                             ( (r_here_plus - r_here_minus2) *    &
     &                               (r_here_plus - r_here_minus )*     &
     &                               (r_here_plus - r_here )*           &
     &                               (r_here_plus - r_here_plus2 )*     &
     &                               (r_here_plus - r_here_plus3 ) )

            coeff_z(i,2) = ( numer_minus2 * numer_minus *               &
     &                               numer * numer_plus *               &
     &                               numer_plus3 ) /                    &
     &                             ( (r_here_plus2 - r_here_minus2) *   &
     &                               (r_here_plus2 - r_here_minus )*    &
     &                               (r_here_plus2 - r_here )*          &
     &                               (r_here_plus2 - r_here_plus )*     &
     &                               (r_here_plus2 - r_here_plus3 ) )

            coeff_z(i,3)= ( numer_minus2 * numer_minus *                &
     &                              numer * numer_plus *                &
     &                              numer_plus2 ) /                     &
     &                             ( (r_here_plus3 - r_here_minus2) *   &
     &                               (r_here_plus3 - r_here_minus )*    &
     &                               (r_here_plus3 - r_here )*          &
     &                               (r_here_plus3 - r_here_plus )*     &
     &                               (r_here_plus3 - r_here_plus2 ) )
          End If

        End Do
!$OMP END PARALLEL DO 

        Else If ( ( ( high_order_scheme  ==  cubicLagrange              &
     &          .or.  high_order_scheme  ==  ECMWF_quasiCubic           &
     &          .or.  high_order_scheme  ==  ECMWF_mono_quasiCubic )    &
     &          .and. L_High )                                          &
     &          .or. ( monotone_scheme  ==  mono_quasiCubic             &
     &         .and. L_Mono ) ) Then

! cubic interpolation

! Calculate interpolation coefficients

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i, r_here,  &
!$OMP& r_here_plus2, r_here_minus, numer_minus, numer, numer_plus,      &
!$OMP& numer_plus2, eta_out, r_here_plus, weight_eta)
        Do i = 1, dim_e_out

          If (k_out(i)  ==  1 .or. k_out(i)  ==  dim_k_in-1) Then
! use linear to interpolate between
! bottom most levels or top most levels.

            r_here = coeff_a(i)                                         &
     &         * r_in (i_out(i),j_out(i),k_out(i))                      &
     &                  + coeff_b(i)                                    &
     &         * r_in (i_out(i)+1,j_out(i),k_out(i))                    &
     &                  + coeff_c(i)                                    &
     &         * r_in (i_out(i),j_out(i)+1,k_out(i))                    &
     &                  + coeff_d(i)                                    &
     &         * r_in (i_out(i)+1,j_out(i)+1,k_out(i))

            r_here_plus = coeff_a(i)                                    &
     &        *  r_in (i_out(i),j_out(i),k_out(i)+1)                    &
     &                  + coeff_b(i)                                    &
     &        *  r_in (i_out(i)+1,j_out(i),k_out(i)+1)                  &
     &                  + coeff_c(i)                                    &
     &        *  r_in (i_out(i),j_out(i)+1,k_out(i)+1)                  &
     &                  + coeff_d(i)                                    &
     &        *  r_in (i_out(i)+1,j_out(i)+1,k_out(i)+1)

!  Obtain an eta coordinate for the departure point.
            weight_eta = ( r_out(i) - r_here ) /                        &
     &                   ( r_here_plus - r_here )

            coeff_z (i,-1) = 0.

            coeff_z (i,0) = 1 - weight_eta

            coeff_z (i,1) = weight_eta

            coeff_z (i,2) = 0.

          Else
! use cubic interpolation.

            r_here = coeff_a(i)                                         &
     &         * r_in (i_out(i),j_out(i),k_out(i))                      &
     &                  + coeff_b(i)                                    &
     &         * r_in (i_out(i)+1,j_out(i),k_out(i))                    &
     &                  + coeff_c(i)                                    &
     &         * r_in (i_out(i),j_out(i)+1,k_out(i))                    &
     &                  + coeff_d(i)                                    &
     &         * r_in (i_out(i)+1,j_out(i)+1,k_out(i))

            r_here_plus = coeff_a(i) *                                  &
     &           r_in (i_out(i),j_out(i),k_out(i)+1)                    &
     &                  + coeff_b(i)                                    &
     &        *  r_in (i_out(i)+1,j_out(i),k_out(i)+1)                  &
     &                  + coeff_c(i)                                    &
     &        *  r_in (i_out(i),j_out(i)+1,k_out(i)+1)                  &
     &                  + coeff_d(i)                                    &
     &        *  r_in (i_out(i)+1,j_out(i)+1,k_out(i)+1)

!  Obtain an eta coordinate for the departure point.
            weight_eta = ( r_out(i) - r_here ) /                        &
     &                   ( r_here_plus - r_here )

            eta_out = eta_in(k_out(i)) + weight_eta *                   &
     &                ( eta_in(k_out(i)+1) - eta_in(k_out(i)) )

! cjs Re-use the r coordinate variables to hold the corresponding
!     eta coordinate values. Interpolation in eta simplifies the
!     algortihm : the knots for vertical interpolation are now the
!     tabulated eta levels, rather than the linearly interpolated r
!     values.

            r_here_minus = eta_in(k_out(i) - 1)
            r_here       = eta_in(k_out(i))
            r_here_plus  = eta_in(k_out(i) + 1)
            r_here_plus2 = eta_in(k_out(i) + 2)


! cjs Compute coordinate differences in eta.

            numer_minus = eta_out - r_here_minus
            numer       = eta_out - r_here
            numer_plus  = eta_out - r_here_plus
            numer_plus2 = eta_out - r_here_plus2


            coeff_z(i,-1)= (numer * numer_plus * numer_plus2 ) /        &
     &                            ( (r_here_minus - r_here) *           &
     &                              (r_here_minus - r_here_plus )*      &
     &                              (r_here_minus - r_here_plus2 ) )
            coeff_z(i,0) = (numer_minus * numer_plus *                  &
     &                             numer_plus2 ) /                      &
     &                            ( (r_here - r_here_minus) *           &
     &                              (r_here - r_here_plus )*            &
     &                              (r_here - r_here_plus2 ) )
            coeff_z(i,1) = (numer_minus * numer *                       &
     &                             numer_plus2 ) /                      &
     &                            ( (r_here_plus - r_here_minus) *      &
     &                              (r_here_plus - r_here )*            &
     &                              (r_here_plus - r_here_plus2 ) )
            coeff_z(i,2) = (numer_minus * numer *                       &
     &                             numer_plus ) /                       &
     &                            ( (r_here_plus2 - r_here_minus) *     &
     &                              (r_here_plus2 - r_here )*           &
     &                              (r_here_plus2 - r_here_plus ) )

          End If

        End Do
!$OMP END PARALLEL DO 

      End If

      If ( ( high_order_scheme  ==  hCubic_vLin .and. L_High )          &
     &      .or. (monotone_scheme  ==  triLinear .and. L_Mono )         &
     &      .or. L_conserv ) Then

! linear interpolation
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(r_here,    &
!$OMP& r_here_plus, weight_eta, i)
        Do i = 1, dim_e_out

! use linear to interpolate between
! bottom most levels or top most levels.

          r_here = coeff_a(i) * r_in (i_out(i),j_out(i),k_out(i))       &
     &           + coeff_b(i) * r_in (i_out(i)+1,j_out(i),k_out(i))     &
     &           + coeff_c(i) * r_in (i_out(i),j_out(i)+1,k_out(i))     &
     &           + coeff_d(i) * r_in (i_out(i)+1,j_out(i)+1,k_out(i))

          r_here_plus = coeff_a(i)                                      &
     &                *  r_in (i_out(i),j_out(i),k_out(i)+1)            &
     &                + coeff_b(i)                                      &
     &                *  r_in (i_out(i)+1,j_out(i),k_out(i)+1)          &
     &                + coeff_c(i)                                      &
     &                *  r_in (i_out(i),j_out(i)+1,k_out(i)+1)          &
     &                + coeff_d(i)                                      &
     &                *  r_in (i_out(i)+1,j_out(i)+1,k_out(i)+1)

!  Obtain an eta coordinate for the departure point.
          weight_eta = ( r_out(i) - r_here ) /                          &
     &                 ( r_here_plus - r_here )

          If (weight_eta<0.0) Then
            coeff_z_lin (i,0) = 1.0
            coeff_z_lin (i,1) = 0.0
          Else If (weight_eta>1.0) Then 
            coeff_z_lin (i,0) = 0.0
            coeff_z_lin (i,1) = 1.0
          Else
            coeff_z_lin (i,0) = 1 - weight_eta
            coeff_z_lin (i,1) = weight_eta
          End If

        End Do
!$OMP   END PARALLEL DO 
      End If

! End of routine.
      IF (lhook) CALL dr_hook('ETA_VERT_WEIGHTS_E',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE eta_vert_weights_e

