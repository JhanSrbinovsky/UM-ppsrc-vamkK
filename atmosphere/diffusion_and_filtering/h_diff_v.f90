! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine H_diff_v

      Subroutine H_diff_v(                                              &
     &                      v, r_at_v, ground,rows,                     &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      at_extremity, proc_row_group,               &
     &                      n_proc, n_procx, n_procy, neighbour,        &
     &                      delta_lambda, delta_phi,                    &
     &                      timestep, n_rows, row_length,               &
     &                      model_levels, model_domain,                 &
     &                      diffusion_coefficient,                      &
     &                      diffusion_order,                            &
     &                      R_v,                                        &
     &                      horizontal_level)

! Purpose:
!          Calculates horizontal diffusion increment to v.
!
! Method:
!          Is described in ;
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, n_rows                                                          &
                         ! number of rows.
     &, rows                                                            &
     &, model_levels                                                    &
                         ! number of model levels.
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, proc_row_group                                                  &
     &, n_proc                                                          &
                   ! Total number of processors
     &, n_procx                                                         &
                   ! Number of processors in longitude
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, neighbour(4)         ! Array with the Ids of the four neighbours
                             ! in the horizontal plane

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, diffusion_order(model_levels)

      Real                                                              &
     &  timestep

      Real                                                              &
     &  diffusion_coefficient(model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_at_v (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:n_rows+halo_j, model_levels)                   &
     &, ground (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:rows+halo_j)

      Real                                                              &
     &  v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,              &
     &     model_levels)

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  R_v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,            &
     &       model_levels)

      Integer                                                           &
     &  horizontal_level  ! level at which steep slope test no
!                         ! longer operates

! Local Variables.

      Integer                                                           &
     &  i, j, k, ka, order      ! Loop indices

      Integer :: kk                     ! level reference

      Real                                                              &
     &  scalar1                                                         &
     &, scalar2                                                         &
     &, sign                                                            &
     &, lower_surface

! Local arrays

      Real                                                              &
     &  lambda_term(row_length, n_rows, model_levels)                   &
     &, phi_term(row_length, n_rows, model_levels)                      &
     &, timestep_over_r_squared(row_length, n_rows, model_levels)       &
     &, field(1-off_x:row_length+off_x,                                 &
     &        1-off_y:n_rows+off_y, model_levels)

      Integer                                                           &
     &  active_diffusion_order(model_levels)                            &
     &, n_active_levels                                                 &
     &, max_diffusion_order

      Real                                                              &
     &  active_diff_coeff(model_levels)
      Integer                                                           &
     &  mask_i(1-off_x:row_length, n_rows, model_levels)                &
     &, mask_j(row_length, 1-off_y:n_rows, model_levels)                &
     &, level_base                                                      &
     &, level_flat

      Integer :: indx(model_levels)    ! index of levels
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle




! No External Routines:

      IF (lhook) CALL dr_hook('H_DIFF_V',zhook_in,zhook_handle)
      scalar1 = 1. / (delta_lambda * delta_lambda)
      scalar2 = 1. / (delta_phi * delta_phi)

! ----------------------------------------------------------------------
! Section 1.   Copy v into field for those model levels with a
!              positive diffusion coefficient.
!              Copy coeffs and diffusion order for all active levels.
! ----------------------------------------------------------------------

      n_active_levels = 0
      max_diffusion_order = 0
      Do k = 1, model_levels
        If (diffusion_coefficient(k)  >   0.0 .and.                     &
     &      diffusion_order(k)  >   0 )Then
          n_active_levels = n_active_levels + 1
          indx(n_active_levels) = k
          active_diffusion_order(n_active_levels) = diffusion_order(k)
          active_diff_coeff(n_active_levels) = diffusion_coefficient(k)
          If (diffusion_order(k)  >   max_diffusion_order )             &
     &        max_diffusion_order = diffusion_order(k)
        End If
      End Do

      level_base = model_levels - n_active_levels
!  Last model-level upon which diffusion is inactive
! ----------------------------------------------------------------------
! Section 2.   Loop over the order to produce order * del squared
! ----------------------------------------------------------------------
!  Horizontal_level should be the first flat level and as such
!  does not need testing for a sloping surface.
!   i.e if level_flat < 2   then
!       diffusion behaves as original (along eta surfaces)
! To switch off slope test completetely then
! either hard-wire level_flat = 1 or ensure that the argument
!    horizontal level <= wet_model_levels - n_active_levels + 1
!      level_flat = 1
      level_flat = n_active_levels -                                   &
     &              model_levels + horizontal_level
      If (level_flat  <=  0) level_flat = 1

!$OMP  PARALLEL DO DEFAULT(NONE)  SCHEDULE(STATIC)                      &
!$OMP& PRIVATE(k,kk,j,i)                                                &
!$OMP& SHARED(off_x, off_y, n_rows, row_length, field, v, r_at_v,       &
!$OMP&        timestep, timestep_over_r_squared, mask_i, mask_j,        &
!$OMP&        n_active_levels, indx)
      Do k = 1, n_active_levels

        kk = indx(k)

        Do j = 1-off_y, n_rows+off_y
          Do i = 1-off_x, row_length+off_x
            field(i,j,k) = v(i,j,kk)
          End Do
        End Do

        Do j = 1, n_rows
          Do i = 1, row_length
            timestep_over_r_squared(i,j,k) =                            &
                            timestep /                                  &
                         ( r_at_v(i,j,kk) * r_at_v(i,j,kk) )
          End Do
        End Do

        Do j = 1, n_rows
          Do i = 0, row_length
           mask_i(i,j,k)=1
          End Do
        End Do

        Do j = 0, n_rows
          Do i = 1, row_length
           mask_j(i,j,k)=1
          End Do
        End Do
      End Do
!$OMP END PARALLEL DO

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& PRIVATE(ka, k, j, i, lower_surface)                              &
!$OMP& SHARED(level_flat, level_base, n_rows, row_length, ground,       &
!$OMP&        r_at_v, mask_i, mask_j)
      Do ka = 1, level_flat - 1
         k = ka + level_base
        Do j = 1, n_rows
          Do i = 0, row_length
            If (k == 1) Then
              lower_surface=( ground(i,j)+ground(i,j+1) )/2.
            Else
              lower_surface=r_at_v(i,j,k-1)
            End If
            If (r_at_v(i+1,j,k)  <   r_at_v(i,j,k+1) .AND.              &
                r_at_v(i+1,j,k)  >   lower_surface) Then
              mask_i(i,j,ka)=1
            Else
              mask_i(i,j,ka)=0
            End If
          End Do
        End Do

        Do j = 0, n_rows
          Do i = 1, row_length
            If (k == 1) Then
              lower_surface=( ground(i,j)+ground(i,j+1) )/2.
            Else
              lower_surface=r_at_v(i,j,k-1)
            End If
            If (r_at_v(i,j+1,k)  <   r_at_v(i,j,k+1) .AND.              &
                r_at_v(i,j+1,k)  >   lower_surface) Then
              mask_j(i,j,ka)=1
            Else
              mask_j(i,j,ka)=0
            End If
          End Do
        End Do
      End Do
!$OMP END PARALLEL DO


      Do order = 1, max_diffusion_order

! ----------------------------------------------------------------------
! Section 2.1  Calculate lambda direction term.
! ----------------------------------------------------------------------

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& PRIVATE(k,j,i)                                                   &
!$OMP& SHARED(order, active_diffusion_order, n_rows, row_length,        &
!$OMP&        lambda_term, phi_term, field, active_diff_coeff, mask_i,  &
!$OMP&        mask_j, scalar1, scalar2, timestep_over_r_squared,        &
!$OMP&        at_extremity, rows, model_domain, n_active_levels)
        Do k = 1, n_active_levels
! Only do calculations for levels whose diffusion order is less than
! the current value
          If (order  <=  active_diffusion_order(k)) Then
            Do j = 1, n_rows
              Do i = 1, row_length
                lambda_term(i,j,k) = (field(i+1,j,k) - field(i,j,k)) *  &
     &                            mask_i(i,j,k)*active_diff_coeff(k)    &
     &                         - (field(i,j,k) - field(i-1,j,k)) *      &
     &                            mask_i(i-1,j,k)*active_diff_coeff(k)

! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! ----------------------------------------------------------------------
                phi_term(i,j,k) = (field(i,j+1,k) - field(i,j,k)) *     &
     &                         mask_j(i,j,k)*active_diff_coeff(k)       &
     &                      - (field(i,j,k) - field(i,j-1,k)) *         &
     &                         mask_j(i,j-1,k)*active_diff_coeff(k)
              End Do
            End Do

! ----------------------------------------------------------------------
! Section 2.3   Calculate new variable.
! ----------------------------------------------------------------------
            Do j = 1, n_rows
              Do i = 1, row_length
                field(i,j,k) = timestep_over_r_squared(i,j,k) *         &
     &                     (lambda_term(i,j,k) * scalar1 +              &
     &                      phi_term(i,j,k) * scalar2 )
              End Do
            End Do
          End If
        End Do
!$OMP END PARALLEL DO

        If (order  /=  max_diffusion_order ) then
! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                   field, row_length, n_rows, n_active_levels,    &
     &                   off_x, off_y, fld_type_v, .TRUE.)
        End If

! End loop over order
      End Do

! ----------------------------------------------------------------------
! Section 3.   Diffusion increment is field * appropriate sign.
! ----------------------------------------------------------------------

! use previously calculated index to calculate which level of field
! maps to which level of input data.
      Do k = 1, n_active_levels
        kk = indx(k)
        sign = (-1) ** (diffusion_order(kk)-1)
        Do j = 1, n_rows
          Do i = 1, row_length
            R_v(i,j,kk) = R_v(i,j,kk) + field(i,j,k) * sign
          End Do
        End Do
      End Do

! End of routine
      IF (lhook) CALL dr_hook('H_DIFF_V',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE H_diff_v

