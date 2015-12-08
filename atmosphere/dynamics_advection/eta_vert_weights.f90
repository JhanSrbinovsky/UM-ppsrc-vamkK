! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine eta_vert_weights
      Subroutine eta_vert_weights (                                     &
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
     &                        weight_lambda, weight_phi )

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
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE um_types, ONLY: integer32

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

      Real                                                              &
     &  r_in (1-halo_i:dim_i_in+halo_i,                                 &
     &        1-halo_j:dim_j_in+halo_j, dim_k_in)                       &
                                                      ! Vertical
                                                      ! co-ordinate
                                                      ! of input data.
     &, eta_in(dim_k_in)                              ! eta levels of
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

! Local Variables.

! scalars

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

      INTEGER (KIND=integer32) :: j_out_copy(dim_i_out,dim_j_out,dim_k_out)

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
! Start of SWW mods

      IF (lhook) CALL dr_hook('ETA_VERT_WEIGHTS',zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(dim_k_out,        &
!$OMP&     dim_j_out,dim_i_out,j_out_copy,j_out,halo_j,halo_i)          &
!$OMP&     SCHEDULE(STATIC)                                     
      do k=1, dim_k_out
      Do j = 1, dim_j_out
        Do i = 1, dim_i_out


             j_out_copy (i,j,k)=j_out(i,j,k)

        End Do
      End Do
      End Do
!$OMP END PARALLEL DO

! end of mods

! ----------------------------------------------------------------------
!  Section 1.   Find model levels just below departure point.
! ----------------------------------------------------------------------

!$OMP PARALLEL DO SCHEDULE(DYNAMIC,2) DEFAULT(NONE) PRIVATE(k)          &
!$OMP& SHARED(eta_in,r_in, check_bottom_levels,L_conserv,               &
!$OMP&        interp_vertical_search_tol,first_flat_level_in,r_out,     &
!$OMP&        dim_i_in, dim_j_in, dim_k_in,dim_i_out, dim_j_out,        &
!$OMP&        dim_k_out,high_order_scheme, monotone_scheme,             &
!$OMP&        model_domain, L_high, L_mono,halo_i, halo_j,coeff_z,      &
!$OMP&        coeff_z_lin,k_out, i_out, j_out_copy,weight_lambda,       &
!$OMP&        weight_phi)
      do k=1,dim_k_out   
! DEPENDS ON: eta_vert_weights_level
       call eta_vert_weights_level (                                    &
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
     &                        k_out, i_out, j_out_copy,                 &
     &                        weight_lambda, weight_phi,k )
      End Do ! end loop over k
!$OMP END PARALLEL DO

! End of routine.
      IF (lhook) CALL dr_hook('ETA_VERT_WEIGHTS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE eta_vert_weights
      
