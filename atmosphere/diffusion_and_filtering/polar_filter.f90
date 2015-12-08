! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine Polar_filter

      Subroutine Polar_filter(                                          &
     &                      theta, u, v, w,                             &
     &                      row_length, rows, n_rows, model_levels,     &
     &                      global_row_length, delta_lambda,            &
     &                      cos_theta_longitude, sin_theta_longitude,   &
     &                      sin_theta_latitude, sin_v_latitude,         &
     &                      north_lat_limit, south_lat_limit,           &
     &                      filt_coeff,                                 &

     &                      n_sweeps, step_per_sweep,                   &
     &                      polar_start_lat_limit,                      &

     &                      off_x, off_y, halo_i, halo_j,               &
     &                      n_proc, n_procx, n_procy, l_datastart,      &
     &                      neighbour, at_extremity, proc_row_group)

! Purpose:
!          Filter fields on rows adjacent to pole and
!          set polar v components to wave 1.
!
! Method:
!          Using 1-2-1 Shapiro filter to zap 2-grid waves
!          T. Davies.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      Use swapable_field_mod, Only :                                    &
          swapable_field_pointer_type


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Field_Types
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels.
     &, global_row_length                                               &
     &, off_x                                                           &
                   ! Size of small halo in i
     &, off_y                                                           &
                   ! Size of small halo in j.
     &, halo_i                                                          &
                   ! Size of halo in i direction.
     &, halo_j                                                          &
                   ! Size of halo in j direction.
     &, l_datastart(3)       ! First gridpoints held by this processor

      Real                                                              &
     &  delta_lambda

      Real                                                              &
     &  cos_theta_longitude(row_length,rows)                            &
                                             ! cos of longitude
     &, sin_theta_longitude(row_length,rows)                            &
                                             ! sin of longitude
     &, sin_theta_latitude (row_length, rows)                           &
     &, sin_v_latitude (row_length, n_rows)

      Real                                                              &
           ! limits for filtering (in radians)
     &  north_lat_limit                                                 &
                         ! latitude above which 1-2-1 filter is applied
     &, south_lat_limit                                                 &
                         ! latitude below which 1-2-1 filter is applied
     &, filt_coeff                                                      &
                         !   filter coefficient
     &, step_per_sweep                                                  &
                         ! amount in radians to increment start latitude
                         ! by per sweep
     &, polar_start_lat_limit ! max latitude at which filter can start

      Integer                                                           &
     &  n_sweeps        ! number of sweeps of filter to do

      Integer                                                           &
     &  n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, proc_row_group                                                  &
     &, neighbour(4)

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

! Arguments with Intent IN/OUT.

      Real, Target ::                                                   &
     &  u (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &     model_levels)                                                &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)                                               &
     &, w(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      0:model_levels)                                             &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                     ! Loop indices
     &, j_start, j_end, i_sweep, gi

      Integer :: i_field     ! counter for swapbounds

      Real                                                              &
     &  u_filt(row_length,rows)                                         &
     &, v_filt(row_length,n_rows)                                       &
     &, w_filt(row_length,rows)                                         &
     &, theta_filt(row_length,rows)                                     &
     &, latitude_theta(rows)                                            &
     &, latitude_v(n_rows)

      Real                                                              &
     &  scalar                                                          &
                         !   filter coefficient smoother
     &, local_north_lat_limit                                           &
     &, local_south_lat_limit

      Real                                                              &
     &  mag_vector_np (model_levels)                                    &
     &, dir_vector_np (model_levels)                                    &
     &, mag_vector_sp (model_levels)                                    &
     &, dir_vector_sp (model_levels)

      Type(swapable_field_pointer_type) :: fields_to_swap(4)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! ----------------------------------------------------------------------
! Section 0. Calulate latitude for each row.
! ----------------------------------------------------------------------
! calculate latitude of each row

      IF (lhook) CALL dr_hook('POLAR_FILTER',zhook_in,zhook_handle)
      Do j = 1, rows
        latitude_theta(j) = asin(sin_theta_latitude(1,j))
      End Do
      Do j = 1, n_rows
        latitude_v(j) = asin(sin_v_latitude(1,j))
      End Do


      Do i_sweep = 1, n_sweeps

        local_north_lat_limit = north_lat_limit + (i_sweep - 1) *       &
     &                          step_per_sweep
        If (local_north_lat_limit  >   polar_start_lat_limit ) Then
          local_north_lat_limit = polar_start_lat_limit
        End If
        local_south_lat_limit = south_lat_limit - (i_sweep - 1) *       &
     &                          step_per_sweep
        If (local_south_lat_limit  <   - polar_start_lat_limit ) Then
          local_south_lat_limit = - polar_start_lat_limit
        End If

! ----------------------------------------------------------------------
! Section 1. 1-2-1 filter desired area near north pole which lies on
!            this processor
! ----------------------------------------------------------------------

! calculate which rows to filter on this processor, default is none
      j_start = 0
      j_end = -1
      Do j = 1, rows
        If (latitude_theta(j)  >   local_north_lat_limit                &
     &       .and. j_start  ==  0) Then
          j_start = j
          j_end = rows
        End If
      End Do

! on north polar processor do not filter polar row
      If (at_extremity(PNorth) .and. j_end  ==  rows) Then
        j_end = rows - 1
      End If

      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = 1, row_length
            u_filt(i,j) = filt_coeff * (u(i-1,j,k) - u(i,j,k) * 2.      &
     &                                    + u(i+1,j,k))
            theta_filt(i,j) = filt_coeff*(theta(i-1,j,k)                &
     &                                    - theta(i,j,k) * 2.           &
     &                                    + theta(i+1,j,k))
          End Do
        End Do
        Do j = j_start, j_end
!         scalar=1.0
          Do i = 1, row_length
            theta(i,j,k) = theta(i,j,k) + theta_filt(i,j)
            u(i,j,k) = u(i,j,k) + u_filt(i,j)
          End Do
!            scalar=scalar+0.1
        End Do
      End Do
      Do k = 1, model_levels-1
        Do j = j_start, j_end
          Do i = 1, row_length
            w_filt(i,j) = filt_coeff*(w(i-1,j,k) - w(i,j,k) * 2.        &
     &                                + w(i+1,j,k))
          End Do
        End Do
        Do j = j_start, j_end
!         scalar=1.0
          Do i = 1, row_length
            w(i,j,k) = w(i,j,k) + w_filt(i,j)
          End Do
!            scalar=scalar+0.1
        End Do
      End Do

! calculate which v rows to filter on this processor, default is none
      j_start = 0
      j_end = -1
      Do j = 1, n_rows
        If (latitude_v(j)  >   local_north_lat_limit                    &
     &      .and. j_start  ==  0) Then
          j_start = j
          j_end = n_rows
        End If
      End Do
      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = 1, row_length
            v_filt(i,j) = filt_coeff * (v(i-1,j,k) - v(i,j,k) * 2.      &
     &                                  + v(i+1,j,k))
          End Do
        End Do
        Do j = j_start, j_end
!            scalar=1.0
          Do i = 1, row_length
             v(i,j,k) = v(i,j,k) + v_filt(i,j)
          End Do
!            scalar=scalar+0.1
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 2. 1-2-1 filter desired area near south pole which lies on
!            this processor
! ----------------------------------------------------------------------

! calculate which rows to filter on this processor
      If (latitude_theta(1)  <   local_south_lat_limit) Then
! at least one row to filter, default to all
        j_start = 1
        j_end = rows
        Do j = 1, rows
          If (latitude_theta(j)  >=  local_south_lat_limit              &
     &        .and. j_end  ==  rows) Then
! limit filtering to previous row
            j_end = j-1
          End If
        End Do
      Else
! no rows to filter
        j_start = 1
        j_end = 0
      End If

! on north polar processor do not filter polar row
      If (at_extremity(PSouth) .and. j_start  ==  1) Then
        j_start = 2
      End If

      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = 1, row_length
            u_filt(i,j) = filt_coeff * (u(i-1,j,k) - u(i,j,k) * 2.      &
     &                                    + u(i+1,j,k))
            theta_filt(i,j) = filt_coeff*(theta(i-1,j,k)                &
     &                                    - theta(i,j,k) * 2.           &
     &                                    + theta(i+1,j,k))
          End Do
        End Do
        Do j = j_start, j_end
!         scalar=1.0
          Do i = 1, row_length
            theta(i,j,k) = theta(i,j,k) + theta_filt(i,j)
            u(i,j,k) = u(i,j,k) + u_filt(i,j)
          End Do
!            scalar=scalar+0.1
        End Do
      End Do

      Do k = 1, model_levels-1
        Do j = j_start, j_end
          Do i = 1, row_length
            w_filt(i,j) = filt_coeff*(w(i-1,j,k) - w(i,j,k) * 2.        &
     &                                + w(i+1,j,k))
          End Do
        End Do
        Do j = j_start, j_end
!         scalar=1.0
          Do i = 1, row_length
            w(i,j,k) = w(i,j,k) + w_filt(i,j)
          End Do
!            scalar=scalar+0.1
        End Do
      End Do

! calculate which rows to filter on this processor
      If (latitude_v(1)  <   local_south_lat_limit) Then
! at least one row to filter, default to all
        j_start = 1
        j_end = n_rows
        Do j = 1, n_rows
          If (latitude_v(j)  >=  local_south_lat_limit                  &
     &        .and. j_end  ==  n_rows) Then
! limit filtering to previous row
            j_end = j-1
          End If
        End Do
      Else
! no rows to filter
        j_start = 1
        j_end = 0
      End If

      Do k = 1, model_levels
        Do j = j_start, j_end
          Do i = 1, row_length
            v_filt(i,j) = filt_coeff * (v(i-1,j,k) - v(i,j,k) * 2.      &
     &                                  + v(i+1,j,k))
          End Do
        End Do
        Do j = j_start, j_end
!            scalar=1.0
          Do i = 1, row_length
             v(i,j,k) = v(i,j,k) + v_filt(i,j)
          End Do
!            scalar=scalar+0.1
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 3. swop haloes for those fields that have been altered
! ----------------------------------------------------------------------

        If (i_sweep  ==  n_sweeps) Then
! DEPENDS ON: polar_vector_wind_n
          Call Polar_vector_wind_n(                                     &
     &                       v,                                         &
     &                       sin_theta_longitude,                       &
     &                       cos_theta_longitude, row_length,           &
     &                       n_rows, model_levels, mag_vector_np,       &
     &                       dir_vector_np, mag_vector_sp,              &
     &                       dir_vector_sp,                             &
     &                       off_x, off_y, global_row_length,           &
     &                       proc_row_group, at_extremity)

          If (at_extremity(PSouth) ) Then
            Do k = 1,model_levels
              Do i = 1, row_length
                gi = l_datastart(1) + i - 1
                u(i,1,k) = - mag_vector_sp(k) * sin ( (gi-.5)*          &
     &                       delta_lambda - dir_vector_sp(k))
              End Do
            End Do
          End If
          If (at_extremity(PNorth) ) Then
            Do k = 1,model_levels
              Do i = 1, row_length
                gi = l_datastart(1) + i - 1
                u(i,rows,k) = mag_vector_np(k) *                        &
     &                            sin ( (gi-.5)*delta_lambda -          &
     &                                dir_vector_np(k))
              End Do
            End Do
          End If

        End If

      i_field = 0
      i_field = i_field + 1
      fields_to_swap(i_field) % field      => u(:,:,:)
      fields_to_swap(i_field) % field_type = fld_type_u
      fields_to_swap(i_field) % levels     = model_levels
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .TRUE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field      => v(:,:,:)
      fields_to_swap(i_field) % field_type = fld_type_v
      fields_to_swap(i_field) % levels     = model_levels
      fields_to_swap(i_field) % rows       = n_rows
      fields_to_swap(i_field) % vector     = .TRUE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field      => theta(:,:,:)
      fields_to_swap(i_field) % field_type = fld_type_p
      fields_to_swap(i_field) % levels     = model_levels
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.

      i_field = i_field + 1
      fields_to_swap(i_field) % field      => w(:,:,1:)
      fields_to_swap(i_field) % field_type = fld_type_p
      fields_to_swap(i_field) % levels     = model_levels - 1
      fields_to_swap(i_field) % rows       = rows
      fields_to_swap(i_field) % vector     = .FALSE.

! DEPENDS ON: swap_bounds_mv
      Call swap_bounds_mv(fields_to_swap, i_field, row_length,         &
                          off_x, off_y)

      End Do ! end loop over number of sweeps

! End of routine
      IF (lhook) CALL dr_hook('POLAR_FILTER',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Polar_filter

