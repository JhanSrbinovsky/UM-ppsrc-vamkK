! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine  IDL_Generate_grid

      Subroutine IDL_Generate_grid(                                     &
     &                      row_length, rows, n_rows, model_levels      &
     &,                     boundary_layer_levels                       &
     &,                     first_constant_rho_level                    &
     &,                     halo_i, halo_j, me                          &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     r_at_u, r_at_v                              &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     a_ixsts,len_a_ixsts, a_spsts,len_a_spsts    &
!  Grid information
     &,                     grid_number, grid_flat, height_domain       &
     &,                     h_o, z_orog_print                           &
     &,                     L_code_test)

! Purpose:
!        Generate r_theta_levels, r_rho_levels given
!        eta_theta_levels, eta_rho_levels and height_domain
!        r_theta_levels(i,j,0) contains orography on input
!        h_o is used to generate reference grid for printing
!        Grid type ( grid_number) and flattening (grid_flat)
!        are required (user options)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE dynamics_grid_mod, ONLY: l_vatpoles

      USE earth_constants_mod, ONLY: earth_radius
      
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
      USE gflat_mod, ONLY: gflat_old, gflat_linear, gflat_linear1,      &
                           gflat_quadratic
      IMPLICIT NONE

      Real                                                              &
        height_domain                                                   &
      , h_o

      Integer                                                           &
        grid_number                                                     &
      , grid_flat

      Integer                                                           &
        me         ! My processor number
      Logical                                                           &
        L_code_test ! User switch


      Integer                                                           &
        row_length                                                      &
                         ! number of points on a row
      , rows                                                            &
                         ! number of rows in a theta field
      , n_rows                                                          &
                         ! number of rows in a v field
      , model_levels                                                    &
                         ! number of model levels
      , boundary_layer_levels                                           &
                                ! number of  boundary_layer_levels
      , first_constant_rho_level                                        &
      , halo_i                                                          &
                             ! Size of halo in i direction.
      , halo_j               ! Size of halo in j direction.

      Real                                                              &
           ! vertical co-ordinate information
        r_theta_levels(1-halo_i:row_length+halo_i,                      &
                       1-halo_j:rows+halo_j,0:model_levels)             &
      , r_rho_levels(1-halo_i:row_length+halo_i,                        &
                     1-halo_j:rows+halo_j, model_levels)                &
      , r_at_u(1-halo_i:row_length+halo_i,                              &
                     1-halo_j:rows+halo_j, model_levels)                &
      , r_at_v(1-halo_i:row_length+halo_i,                              &
                     1-halo_j:n_rows+halo_j, model_levels)              &
      , eta_theta_levels(0:model_levels)                                &
      , eta_rho_levels(model_levels)                                    &
      , z_orog_print(0:model_levels)

      ! Stash indexing arrays
      Integer len_a_ixsts
      Integer len_a_spsts
      Integer a_ixsts(len_a_ixsts)     ! stash index array
      Real    a_spsts(len_a_spsts)     ! atmos stash array

! local variables

      Integer                                                           &
        i,                                                              &
        j,                                                              &
        k,                                                              &
        fcrl

      Real                                                              &
        height_domain_rho                                               &
      , z_ref_theta(model_levels)                                       &
      , z_ref_rho(model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! Description: COMDECK containing vertical grid types
!  for use in idealised problems
!
      INTEGER, PARAMETER :: vert_regular=1
      INTEGER, PARAMETER :: vert_quadratic_theta=21
      INTEGER, PARAMETER :: vert_bi_quadratic=22
      INTEGER, PARAMETER :: vert_quadratic_uvtheta=23
      INTEGER, PARAMETER :: vert_schar=3
      INTEGER, PARAMETER :: vert_dwd=4
      INTEGER, PARAMETER :: vert_stretch_plus_regular=5
      INTEGER, PARAMETER :: vert_quad_stretch_thin=6
      INTEGER, PARAMETER :: vert_regular_thin=7
      INTEGER, PARAMETER :: vert_geometric_theta=8
      INTEGER, PARAMETER :: vert_dump=10
      INTEGER, PARAMETER :: vert_idl_um_grid=11

! No External routines

! ----------------------------------------------------------------------
! Section 0.  Initialise Data fields.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_GENERATE_GRID',zhook_in,zhook_handle)

!   fcrl shorthand for first_constant_rho_level
      fcrl = first_constant_rho_level

! ----------------------------------------------------------------------
! Section 1. Set up vertical co-ordinate arrays relative to surface
!            On entry r_theta_levels(level0) contains orographic
!            height relative to Earth's surface
!            Add Earth-radius to  r_theta_levels and r_rho_levels
!            after section 3
! ----------------------------------------------------------------------
! height_domain is either set for problem or changed in vertical grid

!      If (me == 0) Then
!        Write (6,*) ' '
!        Write (6,*) ' VERTICAL GRID'
!      End If
      z_orog_print(0) =  h_o
      Do k = 1, model_levels
        z_ref_theta(k) = eta_theta_levels(k) * height_domain
        z_ref_rho(k) = eta_rho_levels(k) * height_domain
        z_orog_print(k) = z_ref_theta(k)
      End Do

      if(grid_number  ==  vert_dump) then
        if(me  ==  0)then
          print*,'grid_number = ',grid_number
          print*,'Vertical grid unchanged - as in input dump '
          print*,' BUT    WARNING '
          print*,' grid_flat must be set in your NAMELIST input'
          print*,' For quadratic flattening set grid_flat=3'
        Endif   !(me  ==  0)
      endif  !   grid_number  ==  vert_dump

! Set stash level arrays to the new level heights (zero orography)
! This sets blev and bhlev correctly in the PP header.
! At present the C(k) values (brlev and bhrlev) are set to zero
! as there are various different options in this routine.
! To plot on true height above msl, output height field on model levs
! z(i,j,k) = zsea(k) + C(k)*zorog(i,j)

      If (grid_number  /=  vert_dump) then
        ! Set theta height at surface to 0.0
        a_spsts(a_ixsts(3)) = 0.0
        a_spsts(a_ixsts(4)) = 0.0
        Do k = 1, model_levels
          ! Set rho levels
          a_spsts(a_ixsts(1) + k-1) = z_ref_rho(k)
          ! Set C(k) to zero at present
          a_spsts(a_ixsts(2) + k-1) = 0.0

          ! Set theta levels
          a_spsts(a_ixsts(3) + k) = z_ref_theta(k)
          ! Set C(k) to zero at present
          a_spsts(a_ixsts(4) + k) = 0.0
        End Do
      End If

      if(grid_flat  ==  gflat_linear) then

        if(me  ==  0)then
          print*,'** Linear flattening of grid between surface '
          print*,' and first_constant_rho_level grid_flat = ',grid_flat
        Endif   !(me  ==  0)

        if(grid_number  <=  vert_quadratic_theta)then

! For constant levels set r to be a constant on the level
        Do k = first_constant_rho_level, model_levels
          z_orog_print(k) =  z_ref_theta(k)
          do j = 1-halo_j, rows+halo_j
            do i = 1-halo_i, row_length+halo_i
              r_theta_levels(i,j,k) = z_ref_theta(k)
              r_rho_levels(i,j,k) = z_ref_rho(k)
            End Do
          End Do
        End Do

! From surface to first_constant_rho_level  use linear relaxation
! for flattening
        Do k = 1, first_constant_rho_level - 1
          z_orog_print(k) =  z_ref_theta(k) +                           &
                    h_o * (1.0 - eta_theta_levels(k) /                  &
                    eta_rho_levels(first_constant_rho_level))
          do j = 1-halo_j, rows+halo_j
            do i = 1-halo_i, row_length+halo_i
              r_rho_levels(i,j,k) = z_ref_rho(k) +                      &
                    r_theta_levels(i,j,0) * (1.0 - eta_rho_levels(k) /  &
                    eta_rho_levels(first_constant_rho_level))
              r_theta_levels(i,j,k) = z_ref_theta(k)  +                 &
                  r_theta_levels(i,j,0) * (1.0 - eta_theta_levels(k) /  &
                    eta_rho_levels(first_constant_rho_level))
            End Do
          End Do
        End Do

        elseif(grid_number  ==  vert_bi_quadratic) then
          height_domain_rho = height_domain * 0.5 *                     &
               (eta_theta_levels(model_levels) +                        &
                eta_theta_levels(model_levels - 1))

          do k=1,model_levels
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                if(k <  first_constant_rho_level)then

                  r_theta_levels(i,j,k) = r_theta_levels(i,j,0) +       &
                                      eta_theta_levels(k) *             &
                            (height_domain - r_theta_levels(i,j,0) /    &
                                              eta_theta_levels(fcrl))

                  r_rho_levels(i,j,k)= r_theta_levels(i,j,0) +          &
                                       eta_theta_levels(k) *            &
                         (height_domain_rho - r_theta_levels(i,j,0) /   &
                          eta_theta_levels(fcrl) )
                else

                  r_theta_levels(i,j,k) = eta_theta_levels(k) *         &
                                     height_domain

                  r_rho_levels(i,j,k) = eta_theta_levels(k) *           &
                                    height_domain_rho

                end if
              end do
            end do
          end do

        else     !     grid_number  ==  vert_quadratic_uvtheta

          do k=1,model_levels
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                if(k <  first_constant_rho_level)then

                  r_theta_levels(i,j,k) = r_theta_levels(i,j,0) +       &
                                         eta_theta_levels(k) *          &
                              (height_domain - r_theta_levels(i,j,0) /  &
                                               eta_theta_levels(fcrl))

                  r_rho_levels(i,j,k)= r_theta_levels(i,j,0) +          &
                                         eta_rho_levels(k) *            &
                              (height_domain - r_theta_levels(i,j,0) /  &
                                                eta_rho_levels(fcrl) )
                else

                  r_theta_levels(i,j,k) = eta_theta_levels(k) *         &
                                     height_domain

                  r_rho_levels(i,j,k) = eta_rho_levels(k) *             &
                                    height_domain

                end if
              end do
            end do
          end do

        endif         !    grid_number  <=  vert_quadratic_theta

      elseif (grid_flat  ==  gflat_old) then

        if(me  ==  0)then
          print*,'** Linear flattening of grid between boundary layer '
          print*,' and first_constant_rho_level grid_flat = ',grid_flat
        Endif   !(me  ==  0)

        if(grid_number  <=  vert_quadratic_theta)then

! For boundary layer levels set depth to be constant.
          Do k = 1, boundary_layer_levels
            z_orog_print(k) =  h_o  + z_ref_theta(k)
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                r_theta_levels(i,j,k) = r_theta_levels(i,j,0)           &
                                       + z_ref_theta(k)
                r_rho_levels(i,j,k) = r_theta_levels(i,j,0)             &
                                     + z_ref_rho(k)
              End Do
            End Do
          End Do
! For constant levels set r to be a constant on the level
          Do k = first_constant_rho_level, model_levels
            z_orog_print(k) =  z_ref_theta(k)
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                r_theta_levels(i,j,k) = z_ref_theta(k)
                r_rho_levels(i,j,k) = z_ref_rho(k)
              End Do
            End Do
          End Do
! For intermediate levels use linear relaxation to constant value.
          Do k = boundary_layer_levels+1, first_constant_rho_level-1
            z_orog_print(k) =                                           &
                ( r_rho_levels(1,1,first_constant_rho_level) -          &
                          z_orog_print(boundary_layer_levels) ) *       &
                ( eta_theta_levels(k) -                                 &
                  eta_theta_levels(boundary_layer_levels) ) /           &
                ( eta_rho_levels(first_constant_rho_level) -            &
                  eta_theta_levels(boundary_layer_levels) )             &
                 +  z_orog_print(boundary_layer_levels)
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                r_rho_levels(i,j,k) =                                   &
                  ( r_rho_levels(i,j,first_constant_rho_level) -        &
                    r_theta_levels(i,j,boundary_layer_levels) ) *       &
                  ( eta_rho_levels(k) -                                 &
                    eta_theta_levels(boundary_layer_levels) ) /         &
                  ( eta_rho_levels(first_constant_rho_level) -          &
                    eta_theta_levels(boundary_layer_levels) )           &
                   +  r_theta_levels(i,j,boundary_layer_levels)
                r_theta_levels(i,j,k) =                                 &
                  ( r_rho_levels(i,j,first_constant_rho_level) -        &
                    r_theta_levels(i,j,boundary_layer_levels) ) *       &
                  ( eta_theta_levels(k) -                               &
                    eta_theta_levels(boundary_layer_levels) ) /         &
                  ( eta_rho_levels(first_constant_rho_level) -          &
                    eta_theta_levels(boundary_layer_levels) )           &
                   +  r_theta_levels(i,j,boundary_layer_levels)
              End Do
            End Do
          End Do

        elseif(grid_number  ==  vert_bi_quadratic) then

          height_domain_rho = height_domain * 0.5 *                     &
               (eta_theta_levels(model_levels) +                        &
                eta_theta_levels(model_levels - 1))

          do k = 1, first_constant_rho_level - 1
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                r_theta_levels(i,j,k) = r_theta_levels(i,j,0) +         &
                                        eta_theta_levels(k) *           &
                          (height_domain - r_theta_levels(i,j,0) /      &
                                           eta_theta_levels(fcrl))

                r_rho_levels(i,j,k) = r_theta_levels(i,j,0) +           &
                                           eta_theta_levels(k) *        &
                         (height_domain_rho - r_theta_levels(i,j,0) /   &
                                             eta_theta_levels(fcrl) )
              end do
            end do
          end do

          do k = first_constant_rho_level, model_levels
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                r_theta_levels(i,j,k) = eta_theta_levels(k) *           &
                                     height_domain

                r_rho_levels(i,j,k) = eta_theta_levels(k) *             &
                                    height_domain_rho
              end do
            end do
          end do

        else     !      grid_number  ==  vert_quadratic_uvtheta

          do k = 1, first_constant_rho_level -1
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i

                r_theta_levels(i,j,k) = r_theta_levels(i,j,0) +         &
                                        eta_theta_levels(k) *           &
                          (height_domain - r_theta_levels(i,j,0) /      &
                                           eta_theta_levels(fcrl))

                r_rho_levels(i,j,k)= r_theta_levels(i,j,0) +            &
                                       eta_rho_levels(k) *              &
                          (height_domain - r_theta_levels(i,j,0) /      &
                                            eta_rho_levels(fcrl) )
              end do
            end do
          end do

          do k = first_constant_rho_level, model_levels
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i

                r_theta_levels(i,j,k) = eta_theta_levels(k) *           &
                                    height_domain
                r_rho_levels(i,j,k) = eta_rho_levels(k) *               &
                                    height_domain
              end do
            end do
          end do

        endif         !     grid_number  <=  vert_quadratic_theta

      elseif (grid_flat  ==  gflat_linear1) then
        if(me  ==  0)then
          print*,'** Linear flattening of grid over orography '
          print*,'** starting from level 1.  grid_flat = ',grid_flat
        Endif   !(me  ==  0)

        if(grid_number  <=  vert_quadratic_theta)then

! Only first level is fully terrain-following
          z_orog_print(1) =  h_o  + z_ref_theta(1)

          do j = 1-halo_j, rows+halo_j
            do i = 1-halo_i, row_length+halo_i
              r_theta_levels(i,j,1) = r_theta_levels(i,j,0)             &
                                       + z_ref_theta(1)
              r_rho_levels(i,j,1) = r_theta_levels(i,j,0)               &
                                     + z_ref_rho(1)
            End Do
          End Do

! For constant levels set r to be a constant on the level
          Do k = first_constant_rho_level, model_levels
            z_orog_print(k) =  z_ref_theta(k)
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                r_theta_levels(i,j,k) = z_ref_theta(k)
                r_rho_levels(i,j,k) = z_ref_rho(k)
              End Do
            End Do
          End Do
! From surface to first_constant_rho_level  use linear relaxation
! for flattening
          Do k = 2, first_constant_rho_level-1
            z_orog_print(k) =                                           &
                ( r_rho_levels(1,1,first_constant_rho_level) -          &
                                           z_orog_print(1) ) *          &
                ( eta_theta_levels(k) - eta_theta_levels(1) ) /         &
                ( eta_rho_levels(first_constant_rho_level) -            &
                  eta_theta_levels(1) ) +  z_orog_print(1)
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                r_rho_levels(i,j,k) =                                   &
                  ( r_rho_levels(i,j,first_constant_rho_level) -        &
                    r_theta_levels(i,j,1) ) *                           &
                  ( eta_rho_levels(k) - eta_theta_levels(1) ) /         &
                  ( eta_rho_levels(first_constant_rho_level) -          &
                    eta_theta_levels(1) ) +  r_theta_levels(i,j,1)
                r_theta_levels(i,j,k) =                                 &
                  ( r_rho_levels(i,j,first_constant_rho_level) -        &
                    r_theta_levels(i,j,1) ) *                           &
                  ( eta_theta_levels(k) - eta_theta_levels(1) ) /       &
                  ( eta_rho_levels(first_constant_rho_level) -          &
                    eta_theta_levels(1) ) +  r_theta_levels(i,j,1)
              End Do
            End Do
          End Do

        elseif(grid_number  ==  vert_bi_quadratic) then

          height_domain_rho = height_domain * 0.5 *                     &
               (eta_theta_levels(model_levels) +                        &
                eta_theta_levels(model_levels - 1))

          do k=1,model_levels
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                if(k <  first_constant_rho_level)then

                  r_theta_levels(i,j,k) = r_theta_levels(i,j,0) +       &
                                         eta_theta_levels(k) *          &
                            (height_domain - r_theta_levels(i,j,0) /    &
                                             eta_theta_levels(fcrl))

                  r_rho_levels(i,j,k) = r_theta_levels(i,j,0) +         &
                                         eta_theta_levels(k) *          &
                         (height_domain_rho- r_theta_levels(i,j,0) /    &
                                            eta_theta_levels(fcrl) )
                else

                  r_theta_levels(i,j,k) = eta_theta_levels(k) *         &
                                     height_domain

                  r_rho_levels(i,j,k) = eta_theta_levels(k) *           &
                                    height_domain_rho

                end if
              end do
            end do
          end do

        else     !     grid_number  ==  vert_quadratic_uvtheta

          do k=1,model_levels
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                if(k <  first_constant_rho_level)then

                  r_theta_levels(i,j,k) = r_theta_levels(i,j,0) +       &
                                          eta_theta_levels(k) *         &
                         (height_domain - r_theta_levels(i,j,0) /       &
                                          eta_theta_levels(fcrl))

                  r_rho_levels(i,j,k) = r_theta_levels(i,j,0) +         &
                                          eta_rho_levels(k) *           &
                           (height_domain - r_theta_levels(i,j,0) /     &
                                             eta_rho_levels(fcrl) )
                else

                  r_theta_levels(i,j,k) = eta_theta_levels(k) *         &
                                     height_domain
                  r_rho_levels(i,j,k) = eta_rho_levels(k) *             &
                                    height_domain
                end if
              end do
            end do
          end do

        endif    !    grid_number  <=  vert_quadratic_theta

      elseif (grid_flat  ==  gflat_quadratic) then

        If (me == 0) Then
          Write (Unit=6,Fmt='(A47,A36,I2,A1)')                          &
                '   Quadratic flattening of grid over orography ',      &
                '   starting from surface. (grid_flat = ',grid_flat,')'
        End If

        if(grid_number  <=  vert_quadratic_theta)then

! From surface to first_constant_rho_level use quadratic relaxation
! for flattening

          Do k = 1, first_constant_rho_level - 1
            z_orog_print(k) =  z_ref_theta(k) +                         &
                    h_o * (1.0 - eta_theta_levels(k) /                  &
                    eta_rho_levels(first_constant_rho_level))**2
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                r_rho_levels(i,j,k) = z_ref_rho(k)  +                   &
                     r_theta_levels(i,j,0) * (1.0 - eta_rho_levels(k) / &
                    eta_rho_levels(first_constant_rho_level))**2
                r_theta_levels(i,j,k) = z_ref_theta(k)  +               &
                   r_theta_levels(i,j,0) * (1.0 - eta_theta_levels(k) / &
                    eta_rho_levels(first_constant_rho_level))**2
              End Do
            End Do
          End Do

! For constant levels set r to be a constant on the level
          Do k = first_constant_rho_level, model_levels
            z_orog_print(k) =   z_ref_theta(k)
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                r_theta_levels(i,j,k) = z_ref_theta(k)
                r_rho_levels(i,j,k) = z_ref_rho(k)
              End Do
            End Do
          End Do

        elseif(grid_number  ==  vert_bi_quadratic) then
          height_domain_rho = height_domain * 0.5 *                     &
               (eta_theta_levels(model_levels) +                        &
                eta_theta_levels(model_levels - 1))

          do k =1,model_levels
            do j = 1-halo_j, rows+halo_j
              do i = 1-halo_i, row_length+halo_i
                if(k <  first_constant_rho_level)then

                   r_theta_levels(i,j,k) = r_theta_levels(i,j,0) +      &
                                           eta_theta_levels(k) *        &
                            (height_domain - r_theta_levels(i,j,0) /    &
                                             eta_theta_levels(fcrl))

                   r_rho_levels(i,j,k)= r_theta_levels(i,j,0) +         &
                                        eta_theta_levels(k) *           &
                         (height_domain_rho - r_theta_levels(i,j,0) /   &
                                             eta_theta_levels(fcrl) )
                 else

                   r_theta_levels(i,j,k) = eta_theta_levels(k) *        &
                                     height_domain
                   r_rho_levels(i,j,k) = eta_theta_levels(k) *          &
                                    height_domain_rho
                 end if
               end do
             end do
           end do

         else     !      grid_number  ==  vert_quadratic_uvtheta

           do k=1,model_levels
             do j = 1-halo_j, rows+halo_j
               do i = 1-halo_i, row_length+halo_i
                 if(k <  first_constant_rho_level)then

                   r_theta_levels(i,j,k) = r_theta_levels(i,j,0) +      &
                                           eta_theta_levels(k) *        &
                            (height_domain - r_theta_levels(i,j,0) /    &
                                             eta_theta_levels(fcrl))

                   r_rho_levels(i,j,k) = r_theta_levels(i,j,0) +        &
                                         eta_rho_levels(k) *            &
                             (height_domain - r_theta_levels(i,j,0) /   &
                                               eta_rho_levels(fcrl) )
                 else

                   r_theta_levels(i,j,k) = eta_theta_levels(k) *        &
                                     height_domain
                   r_rho_levels(i,j,k) = eta_rho_levels(k) *            &
                                    height_domain
                 end if
               end do
             end do
           end do

        endif         !    grid_number  <   vert_quadratic_theta

      else

        print*,' **  grid flattening option ',grid_flat,' NOT SUPPORTED'

      endif       !  on grid_flat

! ----------------------------------------------------------------------
! Section 2. Add Earth-radius to  r_theta_levels and r_rho_levels
!            calculate r_at_u points and r_at_v points on rho levels
! ----------------------------------------------------------------------
      do j = 1-halo_j, rows+halo_j
        do i = 1-halo_i, row_length+halo_i
          r_theta_levels(i,j,0) = r_theta_levels(i,j,0) +               &
                                  Earth_radius
        End Do
      End Do

      Do k = 1, model_levels
        do j = 1-halo_j, rows+halo_j
          do i = 1-halo_i, row_length+halo_i
            r_theta_levels(i,j,k) = r_theta_levels(i,j,k) +             &
                                    Earth_radius
            r_rho_levels(i,j,k)   = r_rho_levels(i,j,k) +               &
                                    Earth_radius
          End Do
        End Do
      End Do
! call swap_bounds to set haloes
! DEPENDS ON: swap_bounds
      call Swap_Bounds                                                  &
                        (r_theta_levels,                                &
                         row_length, rows, model_levels+1,              &
                         halo_i, halo_j, fld_type_p, .false.)
! DEPENDS ON: swap_bounds
      call Swap_Bounds                                                  &
                        (r_rho_levels,                                  &
                         row_length, rows, model_levels,                &
                         halo_i, halo_j, fld_type_p, .false.)

      IF (l_vatpoles) THEN
        DO k = 1, model_levels
          DO j = 1-halo_j, rows+halo_j
            DO i = 1-halo_i, row_length+halo_i-1
              r_at_u(i,j,k) = .5 * (r_rho_levels(i,j,k) +                 &
                                         r_rho_levels(i-1,j,k) )
            END DO
          END DO
        END DO
      ELSE
        DO k = 1, model_levels
          DO j = 1-halo_j, rows+halo_j
            DO i = 1-halo_i, row_length+halo_i-1
              r_at_u(i,j,k) = .5 * (r_rho_levels(i,j,k) +                 &
                                         r_rho_levels(i+1,j,k) )
            END DO
          END DO
        END DO
      END IF ! vatpoles

      IF (l_vatpoles) THEN
        DO k = 1, model_levels
          DO j = 1-halo_j, n_rows+halo_j-1
            DO i = 1-halo_i, row_length+halo_i
              r_at_v(i,j,k) = .5 * (r_rho_levels(i,j,k) +                 &
                                       r_rho_levels(i,j-1,k) )
            END DO
          END DO
        END DO
      ELSE
        DO k = 1, model_levels
          DO j = 1-halo_j, n_rows+halo_j-1
            DO i = 1-halo_i, row_length+halo_i
              r_at_v(i,j,k) = .5 * (r_rho_levels(i,j,k) +                 &
                                       r_rho_levels(i,j+1,k) )
            END DO
          END DO
        END DO
      END IF ! vatpoles

! call swap_bounds to set extra points
! DEPENDS ON: swap_bounds
      call Swap_Bounds                                                  &
                        (r_at_u,                                        &
                         row_length, rows, model_levels,                &
                         halo_i, halo_j, fld_type_u, .false.)
! DEPENDS ON: swap_bounds
      call Swap_Bounds                                                  &
                        (r_at_v,                                        &
                         row_length, n_rows, model_levels,              &
                         halo_i, halo_j, fld_type_v, .false.)

      IF (lhook) CALL dr_hook('IDL_GENERATE_GRID',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_Generate_grid

!  End Subroutine IDL_Generate_grid
