! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine pofil_v

      Subroutine pofil_v(                                               &
     &                    v, r_at_v, r_theta_levels, r_rho_levels,      &
     &                    off_x, off_y, halo_i, halo_j,                 &
     &                    sec_v_latitude, cos_theta_latitude,           &
     &                    l_datastart, global_row_length,               &
     &                    delta_lambda, delta_phi,                      &
     &                    me, at_extremity, proc_row_group, n_procy,    &
     &                    rows, n_rows, row_length, model_levels,       &
     &                    max_filter_rows, v_begin, v_end,              &
     &                    v_sweeps, global_v_filter,                    &
     &                    horizontal_level,                             &
     &                    diff_order_wind, diff_coeff_wind,             &
     &                    diff_coeff_v, L_diff_wind )

! Purpose:
!          Filter/diffusion of a v-type
!          Based on conservative diffusion operator.
!          Applied as multiple EW 1-2-1 filter near poles.
!          E-W filter of a v-type field based on E-W first order
!          conservative diffusion operator.
!          Applied as multiple 1-2-1 filter near poles.
!
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
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_diff_wind

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                        ! number of v-rows
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
     &, n_procy                                                         &
                   ! Number of processors in latitude
     &, l_datastart(3)                                                  &
                             ! First gridpoints held by this processor
     &, global_row_length                                               &
     &, me        ! processor id

      Integer                                                           &
     &  max_filter_rows                                                 &
                             !  array dimension for v_begin etc.
     &, v_begin(0:max_filter_rows)                                      &
     &, v_end(0:max_filter_rows)                                        &
     &, v_sweeps(max_filter_rows)                                       &
     &, global_v_filter                                                 &
     &, horizontal_level                                                &
                           ! level at which steep slope test no
!                               ! longer operates
     &, diff_order_wind   ! diffusion order

      Real                                                              &
     &  delta_lambda                                                    &
     &, delta_phi                                                       &
     &, diff_coeff_wind   ! NS diffusion coefficient

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_at_v (1-halo_i:row_length+halo_i,                             &
     &          1-halo_j:n_rows+halo_j, model_levels)                   &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)               &
     &, r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)

      Real                                                              &
           !  trigonometric functions
     &  sec_v_latitude (1-off_x:row_length+off_x,                       &
     &                  1-off_y:n_rows+off_y)                           &
     &, cos_theta_latitude (1-off_x:row_length+off_x,                   &
     &                      1-off_y:rows+off_y)                         &
!   EW diffusion coefficient on v-rows
     &, diff_coeff_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  v (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,              &
     &     model_levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                          ! Loop indices
     &, order                                                           &
     &, j_begin, j_end                                                  &
                          ! Loop bounds
     &, j_store                                                         &
     &, gi                                                              &
     &, level_flat                                                      &
     &, i_filter                                                        &
     &, i_sweep

      Logical                                                           &
     &  L_cycle

      Real                                                              &
     &  sign

! Local arrays

      Real                                                              &
     &  field (1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,          &
     &     model_levels)                                                &
     &, r_theta_at_v(1-off_x:row_length+off_x,                          &
     &               1-off_y:n_rows+off_y, 0:model_levels)              &
     &, delta_z(1-off_x:row_length+off_x,                               &
     &          1-off_y:n_rows+off_y, model_levels )                    &
     &, recip_r_squared_delz(1-off_x:row_length+off_x,                  &
     &                       1-off_y:n_rows+off_y, model_levels )       &
     &, temp(1-off_x:row_length, 1-off_y:n_rows)                        &
     &, diff_lambda(1-off_x:row_length, n_rows)                         &
     &, lambda_term(row_length, n_rows)                                 &
     &, phi_term(row_length, n_rows)                                    &
     &, mask_i(1-off_x:row_length, n_rows, model_levels)                &
     &, mask_j(row_length, 1-off_y:n_rows, model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Set up loop bounds and pointers
!            v_begin and v_end for each potential sweep are set SETCON
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('POFIL_V',zhook_in,zhook_handle)

! To switch off slope test completetely then
! either hard-wire level_flat = 1 or ensure that the argument
!    horizontal level <= wet_model_levels - - model_levels + 1
!      level_flat = 1
      level_flat =  horizontal_level
      if (level_flat  <=  0) level_flat = 1

      if( L_diff_wind ) then

        sign =  -1.0 ** ( diff_order_wind - 1 )

        Do k = 1, model_levels
          Do j = 1-off_y, n_rows+off_y
            Do i = 1-off_x, row_length+off_x
              field(i,j,k) = v(i,j,k)
            Enddo
          EndDo
        EndDo ! k = 1, model_levels

      endif !  L_diff_wind

! ----------------------------------------------------------------------
! Section 2.   Set r field  and calculate delta_z
! ----------------------------------------------------------------------

      Do k = 0, model_levels
        Do j = 1-off_y, n_rows+off_y
          Do i = 1-off_x, row_length+off_x
            r_theta_at_v(i,j,k) = .5 * (r_theta_levels(i+1,j,k) +       &
     &                                  r_theta_levels(i  ,j,k) )
          End Do
        End Do
      End Do ! k = 0, model_levels

! call swap_bounds to set halo points

! DEPENDS ON: swap_bounds
      call Swap_Bounds                                                  &
     &                (r_theta_at_v,                                    &
     &                 row_length, n_rows, model_levels + 1,            &
     &                 off_x, off_y, fld_type_v, .false.)

! calculate D(r)/D(eta) about theta levels
      Do k = 1, model_levels - 1
        Do j = 1-off_y, n_rows+off_y
          Do i = 1-off_x, row_length+off_x
            delta_z(i,j,k) = r_theta_at_v(i,j,k) -                      &
     &                         r_theta_at_v(i,j,k-1)
            recip_r_squared_delz(i,j,k) = 1.0 /                         &
     &                                 ( r_at_v(i,j,k) * r_at_v(i,j,k) *&
     &                                                  delta_z(i,j,k) )
            End Do
          End Do
      End Do   !  k = 1, model_levels - 1
! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at top level and then no need to take special action in later loops
      k = model_levels
      Do j = 1-off_y, n_rows+off_y
        Do i = 1-off_x, row_length+off_x
          delta_z(i,j,k) = 1.0
          recip_r_squared_delz(i,j,k) = 1.0 /                           &
     &                                 ( r_at_v(i,j,k) * r_at_v(i,j,k) )
        End Do
      End Do

! ----------------------------------------------------------------------
! Section 3.   Switch off theta diffusion at steep slopes
! ----------------------------------------------------------------------
!  Horizontal_level should be the first flat level and as such
!  does not need testing for a sloping surface.
!   i.e if level_flat < 2   then
!       diffusion behaves as original (along eta surfaces)

      Do k = 1, model_levels
        Do j = 1, n_rows
          Do i = 1-off_x, row_length
            mask_i(i,j,k) = 1.0
          Enddo
        Enddo
        Do j = 0, n_rows
          Do i = 1, row_length
            mask_j(i,j,k) = 1.0
          Enddo
        Enddo
      Enddo  ! 1, model_levels

      Do j = 1, n_rows
        Do i = 1-off_x, row_length
          if( r_at_v(i,j,1) < r_theta_at_v(i+1,j,0) .or.                &
     &         r_at_v(i+1,j,1) < r_theta_at_v(i,j,0) )                  &
     &        mask_i(i,j,1) = 0.0
        Enddo
      Enddo
      Do j = 0, n_rows
        Do i = 1, row_length
          if( r_at_v(i,j,1) < r_theta_at_v(i,j+1,0) .or.                &
     &         r_at_v(i,j+1,1) < r_theta_at_v(i,j,0) )                  &
     &        mask_j(i,j,1) = 0.0
        Enddo
      Enddo
      Do k = 2, level_flat - 1
        Do j = 1, n_rows
          Do i = 1-off_x, row_length
            if( r_at_v(i,j,k) < r_at_v(i+1,j,k-1) .or.                  &
     &           r_at_v(i+1,j,k) < r_at_v(i,j,k-1) )                    &
     &          mask_i(i,j,k) = 0.0
          Enddo
        Enddo
        Do j = 0, n_rows
          Do i = 1, row_length
            if( r_at_v(i,j,k) < r_at_v(i,j+1,k-1) .or.                  &
     &           r_at_v(i,j+1,k) < r_at_v(i,j,k-1) )                    &
     &         mask_j(i,j,k) = 0.0
          Enddo
        Enddo
      Enddo  ! k = 2, level_flat - 1

! ----------------------------------------------------------------------
! Section 4.0 Diffusion
! ----------------------------------------------------------------------

      if( L_diff_wind ) then

        j_begin = v_begin(0)
        j_end = v_end(0)

        Do order = 1, diff_order_wind

          Do k = 1, model_levels

! ----------------------------------------------------------------------
! Section 4.1  Calculate lambda direction term.
! ----------------------------------------------------------------------

            Do j = j_begin, j_end
              Do i = 1-off_x, row_length
                temp(i,j) = (field(i+1,j,k) * delta_z(i+1,j,k) -        &
     &                        field(i  ,j,k) * delta_z(i  ,j,k) ) *     &
     &                          mask_i(i,j,k) * diff_coeff_v(i,j) *     &
     &                 0.25 * ( r_at_v(i,j,k) + r_at_v(i+1,j,k) ) *     &
     &                        ( r_at_v(i,j,k) + r_at_v(i+1,j,k) )
              End Do
              Do i = 1, row_length
                lambda_term(i,j) = temp(i,j) - temp(i-1,j)
              End Do
            End Do

! ----------------------------------------------------------------------
! Section 4.2  Calculate phi direction term.
! ----------------------------------------------------------------------

            Do j = j_begin-off_y, j_end
              Do i = 1, row_length
                temp(i,j) = (field(i,j+1,k) * delta_z(i,j+1,k) -        &
     &                        field(i,j,  k) * delta_z(i,j, k) ) *      &
     &                              mask_j(i,j,k) * diff_coeff_wind *   &
     &                r_rho_levels(i,j+1,k) * r_rho_levels(i,j+1,k) *   &
     &                                  cos_theta_latitude(i,j+1)
              End Do
            End Do
            Do j = j_begin, j_end
              Do i = 1, row_length
                phi_term(i,j) = ( temp(i,j) - temp(i,j-1) ) *           &
     &                               sec_v_latitude(i,j)
              End Do
            End Do

! ----------------------------------------------------------------------
! Section 4.3   Calculate new variable.
! ----------------------------------------------------------------------

            Do j = j_begin, j_end
              Do i = 1, row_length
                field(i,j,k) = recip_r_squared_delz(i,j,k) *            &
     &                         ( lambda_term(i,j) + phi_term(i,j) )
              End Do
            End Do  ! j = j_begin, j_end

          End Do  ! k = 1, model_levels

! DEPENDS ON: swap_bounds
          call Swap_Bounds(                                             &
     &                     field, row_length, n_rows, model_levels,     &
     &                     off_x, off_y, fld_type_v, .true.)

        EndDo !  order = 1, diff_order_wind

        Do k = 1, model_levels
          Do j = j_begin, j_end
            Do i = 1, row_length
              v(i,j,k) = v(i,j,k) + field(i,j,k) * sign
            End Do
          End Do  ! j = j_begin, j_end
        End Do ! k = 1, model_levels

!   to get correct results there should be a swapbounds call for v
!   here. This is not corrected so HadGEM1 gets the same results.

      endif !  L_diff_wind

! ----------------------------------------------------------------------
! Section 5.0  Polar Filtering
! ----------------------------------------------------------------------

      DO i_filter = 1, global_v_filter

        j_begin = v_begin(i_filter)
        j_end = v_end(i_filter)
        L_cycle = .false.

!  Need to filter hemispheres separately
        if( n_procy == 1 ) L_cycle = .true.

        i_sweep = 1
        Do  ! Sweeping loop  i_sweep = 1, v_sweeps(i_filter)

! ----------------------------------------------------------------------
! Section 5.1   Calculate new variable.
! ----------------------------------------------------------------------
          Do k = 1, model_levels

! There are 2 factors of 1/1 multiplied to give 1/16
            Do j = j_begin, j_end
              Do i = 1-off_x, row_length
                temp(i,j) = (v(i+1,j,k) * delta_z(i+1,j,k) -            &
     &                         v(i  ,j,k) * delta_z(i  ,j,k) ) *        &
     &                                  mask_i(i,j,k) * 0.0625 *        &
     &                     ( r_at_v(i,j,k) + r_at_v(i+1,j,k) ) *        &
     &                     ( r_at_v(i,j,k) + r_at_v(i+1,j,k) )
              End Do

              Do i = 1, row_length
                v(i,j,k) = v(i,j,k) + ( temp(i,j) - temp(i-1,j) ) *     &
     &                                 recip_r_squared_delz(i,j,k)
              End Do
           End Do  ! j = j_begin, j_end

          End Do  ! k = 1, model_levels

          if(n_procy == 1 ) then
! Northern hemisphere needs to be done if n_procy=1
! loop bounds are equvalent of those used so far for S. Hem filter
            If ( L_cycle ) then
              j_store = j_begin
              j_begin = n_rows - j_end + 1
              j_end = n_rows - j_store + 1
              L_cycle = .false.
        CYCLE
            else  ! L_cycle = .false.
! Reset Southern hemisphere pointers for next sweep
              j_begin = v_begin(i_filter)
              j_end = v_end(i_filter)
              L_cycle = .true.
            endif ! L_cycle)
          endif  ! n_procy == 1

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                     v, row_length, n_rows, model_levels,         &
     &                     off_x, off_y, fld_type_v, .true.)

          i_sweep = i_sweep + 1
          if ( i_sweep > v_sweeps(i_filter) ) EXIT
        CYCLE
          EndDo ! sweeping  loop

      EndDO   !  i_filter = 1, global_v_filter


      IF (lhook) CALL dr_hook('POFIL_V',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE pofil_v

