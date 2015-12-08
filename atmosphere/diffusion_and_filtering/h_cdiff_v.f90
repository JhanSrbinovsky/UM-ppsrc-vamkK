! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine H_cdiff_v

      Subroutine H_cdiff_v(                                             &
     &                      v, r_at_v, r_theta_levels, rows,            &
     &                      eta_theta_levels, eta_rho_levels,           &
     &                      cos_theta_latitude, sec_v_latitude,         &
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
!          Calculates horizontal diffusion increment to v
!            in conservative form
! Method:
!          Is described in ;
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
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
     &, diffusion_order(model_levels)                                   &
     &, horizontal_level  ! level at which steep slope test no
!                         ! longer operates

      Real                                                              &
     &  timestep

      Real                                                              &
     &  diffusion_coefficient(model_levels)

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_at_v (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:n_rows+halo_j, model_levels)                   &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, eta_theta_levels (0:model_levels)                               &
     &, eta_rho_levels (model_levels)

      Real                                                              &
           !  trigonometric functions
     &  cos_theta_latitude (1-off_x:row_length+off_x,                   &
     &                    1-off_y:rows+off_y)                           &
     &, sec_v_latitude (1-off_x:row_length+off_x,                       &
     &                    1-off_y:n_rows+off_y)

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

! Local Variables.

      Integer                                                           &
     &  i, j, k, ka, order                                              &
                                ! Loop indices
     &, n_active_levels                                                 &
     &, max_diffusion_order                                             &
     &, level_base                                                      &
     &, level_flat

      Real                                                              &
     &  scalar1                                                         &
     &, scalar2                                                         &
     &, sign                                                            &
     &, recip_delta_phi                                                 &
     &, lower_surface

! Local arrays

      Real                                                              &
     &  lambda_term(row_length, n_rows, model_levels)                   &
     &, phi_term(row_length, n_rows, model_levels)                      &
     &, timestep_over_r_squared(row_length, n_rows, model_levels)       &
     &, field(1-off_x:row_length+off_x,                                 &
     &        1-off_y:n_rows+off_y, model_levels)                       &
     &, r_theta_at_v(1-off_x:row_length+off_x,                          &
     &        1-off_y:n_rows+off_y, 0:model_levels)                     &
     &, drdeta(1-off_x:row_length+off_x,                                &
     &        1-off_y:rows+off_y, model_levels)                         &
     &, temp(1-off_x:row_length+off_x, 1-off_y:rows+off_y)              &
     &, active_diff_coeff(model_levels)

      Integer                                                           &
     &  active_diffusion_order(model_levels)                            &
     &, mask_i(1-off_x:row_length, n_rows, model_levels)                &
     &, mask_j(row_length, 1-off_y:n_rows, model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! No External Routines:

      IF (lhook) CALL dr_hook('H_CDIFF_V',zhook_in,zhook_handle)
      scalar1 = 1. / (delta_lambda * delta_lambda)
      recip_delta_phi = 1.0/ delta_phi
      scalar2 = recip_delta_phi * recip_delta_phi
        Do k = 0, model_levels
         Do j = 1-off_y, n_rows+off_y
           Do i = 1-off_x, row_length+off_x
              r_theta_at_v (i,j,k) = .5 * (r_theta_levels(i,j,k) +      &
     &                               r_theta_levels(i,j+1,k) )
            End Do
          End Do
        End Do

! call swap_bounds to set halo points

! DEPENDS ON: swap_bounds
        call Swap_Bounds                                                &
     &                  (r_theta_at_v,                                  &
     &                   row_length, n_rows, model_levels + 1,          &
     &                   off_x, off_y, fld_type_v, .false.)

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
          active_diffusion_order(n_active_levels) = diffusion_order(k)
          active_diff_coeff(n_active_levels) = diffusion_coefficient(k)
          Do j = 1-off_y, n_rows+off_y
            Do i = 1-off_x, row_length+off_x
              field(i,j,n_active_levels) = v(i,j,k)
            End Do
          End Do
! calculate D(r)/D(eta) about v levels
          if(k  <   model_levels)then
          Do j = 1-off_y, n_rows+off_y
            Do i = 1-off_x, row_length+off_x
            drdeta(i,j,n_active_levels) = (r_theta_at_v(i,j,k) -        &
     &                                     r_theta_at_v(i,j,k-1) )  /   &
     &                     (eta_theta_levels(k) - eta_theta_levels(k-1))
           End Do
         End Do
          else   ! (k = model_levels)
! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at top level and then no need to take special action in later loops
          Do j = 1-off_y, n_rows+off_y
            Do i = 1-off_x, row_length+off_x
            drdeta(i,j,n_active_levels) = 1.0
           End Do
         End Do
         endif          !(k  <   model_levels)
          Do j = 1, n_rows
            Do i = 1, row_length
              timestep_over_r_squared(i,j,n_active_levels) =            &
     &                        timestep /                                &
     &                     ( r_at_v(i,j,k) * r_at_v(i,j,k) *            &
     &                                drdeta(i,j,n_active_levels) )
            End Do
          End Do
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
       if (level_flat  <=  0) level_flat = 1

       Do k = 1, n_active_levels
        Do j = 1, n_rows
          Do i = 0, row_length
           mask_i(i,j,k)=1
          Enddo
        Enddo
        Do j = 0, n_rows
          Do i = 1, row_length
           mask_j(i,j,k)=1
          Enddo
        Enddo
      Enddo

      Do ka = 1, level_flat - 1
         k = ka + level_base
        Do j = 1, n_rows
          Do i = 0, row_length
            if(k == 1)then
              lower_surface=r_theta_at_v(i,j,0)
            else
              lower_surface=r_at_v(i,j,k-1)
            endif
            if(r_at_v(i+1,j,k)  <   r_at_v(i,j,k+1)                     &
     &        .and.                                                     &
     &        r_at_v(i+1,j,k)  >   lower_surface)then
              mask_i(i,j,ka)=1
            else
              mask_i(i,j,ka)=0
            endif
          Enddo
        Enddo
        Do j = 0, n_rows
          Do i = 1, row_length
            if(k == 1)then
              lower_surface=r_theta_at_v(i,j,0)
            else
              lower_surface=r_at_v(i,j,k-1)
            endif
            if(r_at_v(i,j+1,k)  <   r_at_v(i,j,k+1)                     &
     &        .and.                                                     &
     &        r_at_v(i,j+1,k)  >   lower_surface)then
              mask_j(i,j,ka)=1
            else
              mask_j(i,j,ka)=0
            endif
          Enddo
        Enddo
      Enddo


      Do order = 1, max_diffusion_order

! ----------------------------------------------------------------------
! Section 2.1  Calculate lambda direction term.
! ----------------------------------------------------------------------

        Do k = 1, n_active_levels
! Only do calculations for levels whose diffusion order is less than
! the current value
          If (order  <=  active_diffusion_order(k)) Then
            Do j = 1, n_rows
              Do i = 0, row_length
                temp(i,j) = (field(i+1,j,k) * drdeta(i+1,j,k) -         &
     &                       field(i  ,j,k) * drdeta(i  ,j,k) ) *       &
     &                          mask_i(i,j,k) * active_diff_coeff(k)
              End Do
              Do i = 1, row_length
                lambda_term(i,j,k) = temp(i,j) - temp(i-1,j)
              End Do
            End Do

! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! ----------------------------------------------------------------------

            Do j = 0, n_rows
              Do i = 1, row_length
                temp(i,j) = (field(i,j+1,k) * drdeta(i,j+1,k) -         &
     &                       field(i,j,  k) * drdeta(i,j,  k) ) *       &
     &                         mask_j(i,j,k) * active_diff_coeff(k) *   &
     &                         cos_theta_latitude(i,j+1)
              End Do
            End Do
            Do j = 1, n_rows
             Do i = 1, row_length
                phi_term(i,j,k) = (temp(i,j) - temp(i,j-1)) *           &
     &                             sec_v_latitude(i,j)
              End Do
            End Do


          End If ! On if (order  <=  active_diffusion_order(k))
        End Do ! end loop over model levels

! ----------------------------------------------------------------------
! Section 2.3   Calculate new variable.
! ----------------------------------------------------------------------

        Do k = 1, n_active_levels
          If (order  <=  active_diffusion_order(k)) Then
            Do j = 1, n_rows
              Do i = 1, row_length
                field(i,j,k) = timestep_over_r_squared(i,j,k) *         &
     &                     (lambda_term(i,j,k) * scalar1 +              &
     &                      phi_term(i,j,k) * scalar2 )
              End Do
            End Do
          End If
        End Do

        If (order  /=  max_diffusion_order ) then
! DEPENDS ON: swap_bounds
          call Swap_Bounds(                                             &
     &                   field, row_length, n_rows, n_active_levels,    &
     &                   off_x, off_y, fld_type_v, .true.)
        End If

! End loop over order
      End Do

! ----------------------------------------------------------------------
! Section 3.   Diffusion increment is field * appropriate sign.
! ----------------------------------------------------------------------

! use ka to calculate which level of field maps to which level of
! input data.
      ka = 0
      Do k = 1, model_levels
        If (diffusion_coefficient(k)  >   0.0 .and.                     &
     &      diffusion_order(k)  >   0 )Then
          ka = ka + 1
          sign = (-1) ** (diffusion_order(k)-1)
          Do j = 1, n_rows
            Do i = 1, row_length
              R_v(i,j,k) = R_v(i,j,k) + field(i,j,ka) * sign
            End Do
          End Do
        End If
      End Do

! End of routine H_cdiff_v
      IF (lhook) CALL dr_hook('H_CDIFF_V',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE H_cdiff_v

