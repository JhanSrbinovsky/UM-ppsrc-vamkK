! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine h_diff_q
      subroutine h_diff_q                                               &
     &                     (q, r_theta_levels,                          &
     &                      off_x, off_y, halo_i, halo_j,               &
     &                      halo_i_star, halo_j_star,                   &
     &                      at_extremity, proc_row_group,               &
     &                      n_proc, n_procx, n_procy, neighbour,        &
     &                      delta_lambda, delta_phi,                    &
     &                      timestep, rows, row_length,                 &
     &                      model_levels, wet_model_levels,             &
     &                      model_domain, global_row_length,            &
     &                      diffusion_coefficient,                      &
     &                      diffusion_order,                            &
     &                      q_star,                                     &
     &                      horizontal_level)

! Purpose:
!          Calculates horizontal diffusion increment to q.
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


      USE global_2d_sums_mod, ONLY: global_2d_sums
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
     &, rows                                                            &
                         ! number of rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, wet_model_levels                                                &
                             ! number of model levels.
     &, halo_i_star                                                     &
                      ! Size of halo in i direction for star field.
     &, halo_j_star                                                     &
                      ! Size of halo in j direction for star field.
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, global_row_length                                               &
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
     &, diffusion_order(wet_model_levels)

      Real                                                              &
     &  timestep

      Real                                                              &
     &  diffusion_coefficient(wet_model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)

      Real                                                              &
     &  q (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
     &         wet_model_levels)

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  q_star (1-halo_i_star:row_length+halo_i_star,                   &
     &          1-halo_j_star:rows+halo_j_star, wet_model_levels)

      Integer                                                           &
     &  horizontal_level   ! level at which steep slope test no
!                          ! longer operates

! Local Variables.

      Integer                                                           &
     &  i, j, k, ka, order, j0, j1                                      &
                                        ! Loop indices
     &, info

      Integer :: kk                     ! level reference

      Real                                                              &
     &  scalar1                                                         &
     &, scalar2                                                         &
     &, sign

! Local arrays

      Real                                                              &
     &  lambda_term(row_length, rows, model_levels)                     &
     &, phi_term(row_length, rows, model_levels)                        &
     &, timestep_over_r_squared(row_length, rows, model_levels)         &
     &, field(1-off_x:row_length+off_x,                                 &
     &        1-off_y:rows+off_y, model_levels)                         &
     &, l_s_poles(row_length, model_levels)                             &
     &, l_n_poles(row_length, model_levels)                             &
     &, sum_s(model_levels)                                             &
     &, sum_n(model_levels)

      Integer                                                           &

     &  mask_i(1-off_x:row_length, rows, model_levels)                  &
     &, mask_j(row_length, 1-off_y:rows, model_levels)

      Integer :: indx(model_levels)       ! index of levels

      Integer                                                           &
     &  active_diffusion_order(model_levels)                            &
     &, n_active_levels                                                 &
     &, max_diffusion_order                                             &
     &, level_base                                                      &
     &, level_flat

      Real                                                              &
     &  active_diff_coeff(model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! No External Routines:

      IF (lhook) CALL dr_hook('H_DIFF_Q',zhook_in,zhook_handle)
      scalar1 = 1. / (delta_lambda * delta_lambda)
      scalar2 = 1. / (delta_phi * delta_phi)
      j0 = 1
      j1 = rows
      If (model_domain  ==  1) Then
        If (at_extremity(PSouth)) j0 = 2
        If (at_extremity(PNorth)) j1 = rows-1
      End If

! ----------------------------------------------------------------------
! Section 1.   Copy q into field for those model levels with a
!              positive diffusion coefficient.
!              Copy coeffs and diffusion order for all active levels.
! ----------------------------------------------------------------------

      n_active_levels = 0
      max_diffusion_order = 0
      Do k = 1, wet_model_levels
        If (diffusion_coefficient(k)  >   0.0 .AND.                     &
            diffusion_order(k)  >   0 ) Then
          n_active_levels = n_active_levels + 1
          indx(n_active_levels) = k
          active_diffusion_order(n_active_levels) = diffusion_order(k)
          active_diff_coeff(n_active_levels) = diffusion_coefficient(k)
          If (diffusion_order(k)  >   max_diffusion_order )             &
              max_diffusion_order = diffusion_order(k)
        End If
      End Do

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& PRIVATE(k,kk,j,i)                                                &
!$OMP& SHARED(off_x, off_y, rows, row_length, field, q, r_theta_levels, &
!$OMP&        timestep, timestep_over_r_squared, mask_i, mask_j,        &
!$OMP&        n_active_levels, indx)
      Do k = 1, n_active_levels

        kk = indx(k)
        Do j = 1-off_y, rows+off_y
          Do i = 1-off_x, row_length+off_x
            field(i,j,k) = q(i,j,kk)
          End Do
        End Do

        Do j = 1, rows
          Do i = 1, row_length
            timestep_over_r_squared(i,j,k) =  timestep /                &
                                            ( r_theta_levels(i,j,kk) *  &
                                              r_theta_levels(i,j,kk) )
           End Do
        End Do

        Do j = 1, rows 
          Do i = 0, row_length 
            mask_i(i,j,k)=1 
          End Do 
        End Do 

        Do j = 0, rows 
          Do i = 1, row_length 
           mask_j(i,j,k)=1 
          End Do 
        End Do 
      End Do 
!$OMP END PARALLEL DO 

      level_base = wet_model_levels - n_active_levels
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
     &              wet_model_levels + horizontal_level
       If (level_flat  <=  0) level_flat = 1

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& PRIVATE(k,ka,j,i)                                                &
!$OMP& SHARED(level_flat, level_base, rows, row_length, r_theta_levels, &
!$OMP&        mask_i, mask_j)
      Do ka = 1, level_flat - 1
        k = ka + level_base
        Do j = 1, rows
          Do i = 0, row_length
           If (r_theta_levels(i+1,j,k) <  r_theta_levels(i,j,k+1) .AND. &
               r_theta_levels(i+1,j,k) >  r_theta_levels(i,j,k-1)) Then
             mask_i(i,j,ka)=1
           Else
             mask_i(i,j,ka)=0
           Endif
          End Do
        End Do

        Do j = 0, rows
          Do i = 1, row_length
           If (r_theta_levels(i,j+1,k) <  r_theta_levels(i,j,k+1) .AND. &
               r_theta_levels(i,j+1,k) >  r_theta_levels(i,j,k-1)) THEN
             mask_j(i,j,ka)=1
           Else
             mask_j(i,j,ka)=0
           Endif
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
!$OMP& SHARED(n_active_levels, active_diffusion_order, lambda_term,     &
!$OMP&        field, mask_i, mask_j, active_diff_coeff, j0, j1, rows,   &
!$OMP&        row_length, model_domain, at_extremity, phi_term,         &
!$OMP&        l_s_poles, l_n_poles, order)
        Do k = 1, n_active_levels
! Only do calculations for levels whose diffusion order is less than
! the current value
          If (order  <=  active_diffusion_order(k)) Then
            Do j = 1, rows
              Do i = 1, row_length
                lambda_term(i,j,k) = (field(i+1,j,k) - field(i,j,k)) *  &
     &                            mask_i(i,j,k)*active_diff_coeff(k)    &
     &                         - (field(i,j,k) - field(i-1,j,k)) *      &
     &                            mask_i(i-1,j,k)*active_diff_coeff(k)
              End Do
            End Do

! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! ----------------------------------------------------------------------

            Do j = j0, j1
              Do i = 1, row_length
                phi_term(i,j,k) = (field(i,j+1,k) - field(i,j,k)) *     &
     &                         mask_j(i,j,k)*active_diff_coeff(k)       &
     &                      - (field(i,j,k) - field(i,j-1,k)) *         &
     &                         mask_j(i,j-1,k)*active_diff_coeff(k)
              End Do
            End Do

            If (model_domain  ==  1) Then
              If(at_extremity(PSouth))then
                Do i = 1, row_length
                  l_s_poles(i,k) = (field(i,2,k) - field(i,1,k)) *      &
     &                             active_diff_coeff(k)
                End Do
              End If
              If(at_extremity(PNorth))then
                Do i = 1, row_length
                  l_n_poles(i,k)=(field(i,rows-1,k) - field(i,rows,k))* &
     &                           active_diff_coeff(k)
                End Do
              End If
            End If

          End If ! On if (order  <=  active_diffusion_order(k))
        End Do ! end loop over model levels
!$OMP END PARALLEL DO

        If (model_domain  ==  1) Then
          If (at_extremity(PSouth)) Then
            CALL global_2d_sums(l_s_poles, row_length, 1, 0, 0,         &
                                n_active_levels, sum_s,                 &
                                proc_row_group)

            Do k = 1, n_active_levels
              sum_s(k) = sum_s(k) * 2. / global_row_length
              Do i = 1, row_length
                phi_term(i,1,k) = sum_s(k)
              End Do
            End Do
          End If

          If (at_extremity(PNorth)) Then
            CALL global_2d_sums(l_n_poles, row_length, 1, 0, 0,         &
                                n_active_levels, sum_n,                 &
                                proc_row_group)

            Do k = 1, n_active_levels
              sum_n(k) = sum_n(k) * 2. / global_row_length
              Do i = 1, row_length
                phi_term(i,rows,k) = sum_n(k)
              End Do
            End Do
          End If
        End If

! ----------------------------------------------------------------------
! Section 2.3   Calculate new variable.
! ----------------------------------------------------------------------

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& PRIVATE(k,j,i)                                                   &
!$OMP& SHARED(n_active_levels, order, active_diffusion_order, field,    &
!$OMP&        rows, row_length, timestep_over_r_squared, lambda_term,   &
!$OMP&        scalar1, scalar2, phi_term)
        Do k = 1, n_active_levels
          If (order  <=  active_diffusion_order(k)) Then
            Do j = 1, rows
              Do i = 1, row_length
                field(i,j,k) = timestep_over_r_squared(i,j,k) *         &
     &                       (lambda_term(i,j,k) * scalar1 +            &
     &                        phi_term(i,j,k) * scalar2 )
              End Do
            End Do
          End If
        End Do
!$OMP END PARALLEL DO

        If (order  /=  max_diffusion_order ) Then
! DEPENDS ON: swap_bounds
          call Swap_Bounds(                                             &
     &                   field, row_length, rows, n_active_levels,      &
     &                   off_x, off_y, fld_type_p, .false.)
        End If

! End loop over order
      End Do

! ----------------------------------------------------------------------
! Section 3.   Diffusion increment is field * appropriate sign.
! ----------------------------------------------------------------------

! use previously calculated index to assertain which level of field
! maps to which level of input data.
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& PRIVATE(k,j,i,kk,sign)                                           &
!$OMP& SHARED(n_active_levels, indx, diffusion_order, rows, row_length, &
!$OMP&        q_star, field)
      Do k = 1, n_active_levels
        kk = indx(k)
        sign = (-1) ** (diffusion_order(kk)-1)
        Do j = 1, rows
          Do i = 1, row_length
            q_star(i,j,kk) = q_star(i,j,kk) +                           &
                                field(i,j,k) * sign
          End Do
        End Do
      End Do
!$OMP END PARALLEL DO

! End of routine
      IF (lhook) CALL dr_hook('H_DIFF_Q',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE h_diff_q

