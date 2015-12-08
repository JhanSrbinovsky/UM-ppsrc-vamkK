! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine pofil_new
      Subroutine pofil_new(                                             &
     &                     field, fld_type, off_u, off_v,               &
     &                     base_f, base_h, base_i, base_j,              &
     &                     levels, model_levels, active_levels,         &
     &                     metric_layer,                                &
     &                     in_rows, cos_rows, rows, row_length,         &
     &                     r_levels, r_half_levels, r_midi, r_midj,     &
     &                     off_x, off_y, halo_x, halo_y,                &
     &                     halo_ri, halo_rj, halo_hi, halo_hj,          &
     &                     halo_ii, halo_ij, halo_ji, halo_jj,          &
     &                     sec_latitude, cos_latitude, grad_theta1,     &
     &                     model_domain, at_extremity, n_procy,         &
     &                     max_filter_rows, global_filter,              &
     &                     sweeps, begin, end, horizontal_level,        &
     &                     diff_coeff_phi, diff_coeff, L_diff, L_vector,&
     &                     L_pofil_hadgem2 )

! Purpose:
!          Filter/diffusion based on first order
!          conservative diffusion operator.
!          Applied as multiple EW 1-2-1 filter near poles.
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

      Logical, Intent(In) ::                                            &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_diff                                                          &
                         ! true if diffusion active
     &, L_vector         ! true if vector field
      Logical, Intent(In) :: L_pofil_hadgem2       ! setting for hadgem2

      Integer, Intent(In) ::                                            &
     &  max_filter_rows                                                 &
                         ! max dimension for sweeping arrays
     &, global_filter                                                   &
                        ! max number of filter sweeps 0=diffusion only
     &, row_length                                                      &
                          ! number of point on a row.
     &, rows                                                            &
                          ! number of rows.
     &, in_rows                                                         &
                          ! number of rows in field
     &, cos_rows                                                        &
                          ! number of rows for cos lat
     &, levels                                                          &
                          ! number of levels in field
     &, active_levels                                                   &
                          ! number of levels to be filtered
     &, model_levels                                                    &
                          ! number of model levels.
     &, metric_layer                                                    &
                        ! uppermost non-constant layer for metric terms
     &, base_f                                                          &
                      ! full levels base dim; 0 for theta, 1 for rho
     &, base_h                                                          &
                      ! half levels base dim; 0 for theta, 1 for rho
     &, base_i                                                          &
                      ! start level for midi levels
     &, base_j                                                          &
                      ! start level for midj levels
     &, off_u                                                           &
                      ! 1 for u points, 0 for v,p points
     &, off_v                                                           &
                      ! 1 for v points, 0 for u,p points
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, halo_x                                                          &
                      ! Field halo in i direction.
     &, halo_y                                                          &
                      ! Field halo in j direction.
     &, halo_hi                                                         &
                      ! Size of halo in i for r on half-levels
     &, halo_hj                                                         &
                      ! Size of halo in j for r on half-levels
     &, halo_ri                                                         &
                      ! Size of  halo_i for field r levels
     &, halo_rj                                                         &
                      ! Size of  halo_j for field r levels
     &, halo_ii                                                         &
                      ! Size of  halo_i for r at mid i points
     &, halo_ij                                                         &
                      ! Size of  halo_j for r at mid i points
     &, halo_ji                                                         &
                      ! Size of  halo_i for r at mid j points
     &, halo_jj                                                         &
                      ! Size of  halo_j for r at mid j points
     &, model_domain                                                    &
     &, fld_type                                                        &
                      ! field type (p=1, u=2 or v=3)
     &, n_procy    ! Number of processors in latitude

      Integer, Intent(In) ::                                            &
     &  sweeps(max_filter_rows)                                         &
     &, begin(0:max_filter_rows)                                        &
     &, end(0:max_filter_rows)                                          &
     &, horizontal_level   ! level at which steep slope test no
!                               ! longer operates

      Real, Intent(In) ::                                               &
     &  diff_coeff_phi                                                  &
                         ! NS diffusion coefficient
     &, grad_theta1      ! level 1 theta gradient test value

      Real, Intent(In) ::                                               &
     &  r_levels (1-halo_ri:row_length+halo_ri,                         &
     &            1-halo_rj:in_rows+halo_rj, base_f:model_levels)       &
     &, r_half_levels (1-halo_hi:row_length+halo_hi,                    &
     &                 1-halo_hj:in_rows+halo_hj, base_h:model_levels)  &
     &, r_midi (1-halo_ii:row_length+halo_ii,                           &
     &          1-halo_ij:in_rows+halo_ij, base_i:model_levels)         &
     &, r_midj (1-halo_ji:row_length+halo_ji,                           &
     &          1-halo_jj:cos_rows+halo_jj, base_j:model_levels)

      Real, Intent(In) ::                                               &
     &  sec_latitude (1-off_x:row_length+off_x,                         &
     &                  1-off_y:in_rows+off_y)                          &
     &, cos_latitude (1-off_x:row_length+off_x,                         &
     &                      1-off_y:cos_rows+off_y)                     &
!   EW diffusion coefficient
     &, diff_coeff(1-off_x:row_length+off_x, 1-off_y:in_rows+off_y)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real, Intent(InOut) ::                                            &
     &  field (1-halo_x:row_length+halo_x, 1-halo_y:in_rows+halo_y,     &
     &         levels)

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                          ! Loop indices
     &, k_start                                                         &
                          ! Loop bound
     &, j_begin, j_end                                                  &
                          ! Loop bounds
     &, j_store                                                         &
     &, level_flat                                                      &
     &, i_filter                                                        &
     &, i_sweep                                                         &
     &, num_pass                                                        &
     &, count                                                           &
     &, i_cold(100)                                                     &
     &, j_cold(100)

      Logical                                                           &
     &  L_cycle                                                         &
     &, L_combine

! Local arrays

      Real                                                              &
     &  delta_z(1-off_x:row_length+off_x, 1-off_y:in_rows+off_y,        &
     &                                               metric_layer )     &
     &, recip_r_squared_delz( 1-off_x:row_length+off_x,                 &
     &                        1-off_y:in_rows+off_y, metric_layer )     &
     &, temp( 1-off_x:row_length, 1-off_y:in_rows)                      &
     &, mask_i( 1-off_x:row_length, in_rows, metric_layer )             &
     &, mask_j( row_length, 1-off_y:in_rows, metric_layer )

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Set metric terms and calculate D(r)/D(eta)
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('POFIL_NEW',zhook_in,zhook_handle)

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) SHARED(off_x, delta_z, &
!$OMP& recip_r_squared_delz, r_levels, in_rows, off_y, row_length,       &
!$OMP& metric_layer, r_half_levels, base_h) PRIVATE(i, j, k)
      Do k = 1, metric_layer
!   No metric terms needed above first constant rho-level
!   since they cancel on constant r-surfaces

          Do j = 1-off_y, in_rows+off_y
            Do i = 1-off_x, row_length+off_x
              delta_z(i,j,k) = r_half_levels(i,j,k+base_h) -            &
     &                         r_half_levels(i,j,k+base_h-1)
              recip_r_squared_delz(i,j,k) = 1.0 / ( r_levels(i,j,k) *   &
     &                                              r_levels(i,j,k) *   &
     &                                              delta_z(i,j,k) )
            End Do
          End Do

      End Do   !  k = 1, metric_layer
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Section 3.   Switch off diffusion at steep slopes
! ----------------------------------------------------------------------
!  level_flat is the uppermost non-flat level and as such
!  does not need testing for a sloping surface.

      Do k = 1, metric_layer
        Do j = 1, in_rows
          Do i = 1-off_x, row_length
            mask_i(i,j,k) = 1.0
          Enddo
        Enddo
        Do j = 0, in_rows
          Do i = 1, row_length
            mask_j(i,j,k) = 1.0
          Enddo
        Enddo
      Enddo  ! 1, metric_layer

      if (.NOT. L_pofil_hadgem2) then
!   turn off the check below for Hadgem2 runs

      if( horizontal_level > 0 ) then

        level_flat = min(horizontal_level, metric_layer)

        k_start = 1
        if ( base_f == 1 ) then
! for rho levels, bottom level can only be tested against surface
          Do j = 1, in_rows
            Do i = 1-off_x, row_length
              if( r_levels(i,j,1) < r_half_levels(i+1,j,0) .or.         &
     &            r_levels(i+1,j,1) < r_half_levels(i,j,0) )            &
     &           mask_i(i,j,1) = 0.0
            Enddo
          Enddo

          Do j = 0, in_rows
            Do i = 1, row_length
              if( r_levels(i,j,1) < r_half_levels(i,j+1,0) .or.         &
     &            r_levels(i,j+1,1) < r_half_levels(i,j,0) )            &
     &            mask_j(i,j,1) = 0.0
            Enddo
          Enddo

          k_start = 2
        endif !  base_f == 1

        Do k = k_start, level_flat
          Do j = 1, in_rows
            Do i = 1-off_x, row_length
              if( r_levels(i,j,k) < r_levels(i+1,j,k-1) .or.            &
     &             r_levels(i+1,j,k) < r_levels(i,j,k-1) )              &
     &            mask_i(i,j,k) = 0.0
            Enddo
          Enddo
         Do j = 0, in_rows
            Do i = 1, row_length
              if( r_levels(i,j,k) < r_levels(i,j+1,k-1) .or.            &
     &             r_levels(i,j+1,k) < r_levels(i,j,k-1) )              &
     &           mask_j(i,j,k) = 0.0
            Enddo
          Enddo
        Enddo  ! k = k_start, level_flat

      endif ! horizontal_level > 0

      endif ! L_pofil_hadgem2

      if ( grad_theta1 > 1.0 ) then
        count = 0

        do k = 1, 100
          i_cold(k) = 0
          j_cold(k) = 0
        enddo ! k = 1, 100

! test for too cold level 1 theta' look for large horizontal gradient
! re-set mask to 1 to activate diffusion for steep slopes

        Do j = 1, in_rows
          Do i = 1-off_x, row_length
            if( abs(field(i+1,j,1) - field(i,j,1)) > grad_theta1 ) then
              mask_i(i,j,1) = 1.0
              count = count + 1
              if( count < 101 ) then
                i_cold(count) = i
                j_cold(count) = j
              endif ! count < 101
            endif ! abs(field(i+1,j,1) - field(i,j,1)) > grad_theta1
          Enddo
        Enddo
        Do j = 0, in_rows
          Do i = 1, row_length
            if( abs(field(i,j+1,1) - field(i,j,1)) > grad_theta1 ) then
              mask_j(i,j,1) = 1.0
              count = count + 1
              if( count < 101 ) then
                i_cold(count) = i
                j_cold(count) = j
              endif ! count < 101
            endif ! abs(field(i,j+1,1) - field(i,j,1)) > grad_theta1
          Enddo
        Enddo

        write(6,*) ' '
        write(6,*)'Test of bottom level theta gradient'
        write(6,*) count,' points have temperature difference > '       &
     &         , grad_theta1
        do k = 1, count
         write(6,*) k ,' at i,j = ',i_cold(k),j_cold(k)
        enddo ! k = 1, count

      endif ! grad_theta1 > 1.0

! ----------------------------------------------------------------------
! Section 4.0  Filtering and diffusion
! ----------------------------------------------------------------------

      L_combine = .false.
      num_pass = global_filter
      j_begin = begin(1)
      j_end = end(1)
! Combine EW diffusion with first sweep of polar filter
      if ( L_diff ) then
        L_combine = .true.
        if( global_filter < 1 ) num_pass = 1
        if ( j_begin < 0 ) then
          j_begin = begin(0)
          j_end = end(0)
        elseif ( j_end < end(0) ) then
          j_end = end(0)
        elseif ( j_begin > begin(0) ) then
          j_begin = begin(0)
        endif !  j_begin < 0
      endif ! L_diff

      DO i_filter = 1, num_pass

        if( i_filter > 1 ) then
          L_combine = .false.
          j_begin = begin(i_filter)
          j_end = end(i_filter)
        endif ! i_filter > 1

       L_cycle = .false.

!  If n_procy = 1, need to filter hemispheres separately
        if( n_procy == 1 ) L_cycle = .true.

        i_sweep = 1

        Do  ! Sweeping loop  i_sweep = 1, sweeps(i_filter)

! ----------------------------------------------------------------------
! Section 4.1   EW filtering and diffusion
! ----------------------------------------------------------------------

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(active_levels,  &
!$OMP& j_begin, j_end, row_length, metric_layer, field, delta_z, r_midi, &
!$OMP& diff_coeff, mask_i, recip_r_squared_delz, off_u) PRIVATE(temp, i, &
!$OMP& j, k)
        Do k = 1, active_levels

          if ( k <= metric_layer ) then
!   No mask or metric terms since levels are horizontal

            Do j = j_begin, j_end
              Do i = 0, row_length
                temp(i,j) = ( field(i+1,j,k) * delta_z(i+1,j,k) -       &
     &                        field(i  ,j,k) * delta_z(i  ,j,k) ) *     &
     &                         mask_i(i,j,k) * diff_coeff(i,j) *        &
     &                   r_midi(i+off_u,j,k) * r_midi(i+off_u,j,k)
              End Do

              Do i = 1, row_length
                field(i,j,k) = field(i,j,k) + (temp(i,j) - temp(i-1,j)) &
     &                                   * recip_r_squared_delz(i,j,k)
              End Do
           End Do  ! j = j_begin, j_end

          else  ! k > metric_layer

            Do j = j_begin, j_end

              Do i = 0, row_length
                temp(i,j) = ( field(i+1,j,k) - field(i,j,k) ) *         &
     &                                      diff_coeff(i,j)
      
              End Do

              Do i = 1, row_length
                field(i,j,k) = field(i,j,k) + temp(i,j) - temp(i-1,j)
              End Do

            End Do  ! j = j_begin, j_end

          endif ! k <= metric_layer

        End Do  ! k = 1, active_levels
!$OMP END PARALLEL DO

! If n_procy=1, Northern hemisphere needs to be done after S. Hem
          If ( L_cycle ) then
            j_store = j_begin
            if ( i_sweep == 1 .and. L_combine ) then
              j_begin = j_end + 1
            else
              j_begin = in_rows - j_end + 1
            endif ! i_sweep == 1 .and. L_combine
            j_end = in_rows - j_store + 1
            L_cycle = .false.
        CYCLE
          endif ! L_cycle)

! Reset pointers for next sweep
! either because N Hem has just been done or
! 1st sweep was combined filter and diffusion
          j_begin = begin(i_filter)
          j_end = end(i_filter)
          if( n_procy == 1 ) L_cycle = .true.

! DEPENDS ON: swap_bounds
          Call Swap_Bounds(                                             &
     &                     field, row_length, in_rows, levels,          &
     &                     halo_x, halo_y, fld_type, L_vector)

          i_sweep = i_sweep + 1
          if ( i_sweep > sweeps(i_filter) ) EXIT
        CYCLE
          EndDo ! sweeping  loop

      EndDO   !  i_filter = 1, num_pass

! ----------------------------------------------------------------------
! Section 4.2   NS diffusion
! ----------------------------------------------------------------------

      if( L_diff ) then
        j_begin = begin(0)
        j_end = end(0)
        if( n_procy == 1 ) j_end = in_rows - j_begin + 1

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) SHARED(active_levels,  &
!$OMP& j_begin, j_end, row_length, metric_layer, field, delta_z, r_midj, &
!$OMP& diff_coeff_phi, mask_j, recip_r_squared_delz, cos_latitude,       &
!$OMP& sec_latitude, fld_type,  model_domain, off_v, at_extremity)       &
!$OMP& PRIVATE(temp, i, j, k)
      Do k = 1, active_levels

        if ( k <= metric_layer ) then
!   No mask or metric terms since levels are horizontal

          Do j = j_begin-1, j_end
            Do i = 1, row_length
              temp(i,j) = (field(i,j+1,k) * delta_z(i,j+1,k) -          &
     &                     field(i,j,  k) * delta_z(i,j, k) ) *         &
     &                      mask_j(i,j,k) * diff_coeff_phi *            &
     &                r_midj(i,j+off_v,k) * r_midj(i,j+off_v,k) *       &
     &                                  cos_latitude(i,j+off_v)
            End Do
          End Do

          Do j = j_begin, j_end
            Do i = 1, row_length
              field(i,j,k) = field(i,j,k) + (temp(i,j) - temp(i,j-1)) * &
     &                                              sec_latitude(i,j) * &
     &                                      recip_r_squared_delz(i,j,k)
            End Do
          End Do

          If ( model_domain == mt_global .and.                          &
     &             fld_type == fld_type_v ) Then
            If (at_extremity(PSouth) ) Then
              j = j_begin-1
              Do i = 1, row_length
                field(i,j,k) = field(i,j,k) + temp(i,j) *               &
     &                                sec_latitude(i,j) *               &
     &                                recip_r_squared_delz(i,j,k)
              End Do
            endIf !at_extremity(PSouth)
            If (at_extremity(PNorth) ) Then
              j = j_end + 1
              Do i = 1, row_length
                field(i,j,k) = field(i,j,k) - temp(i,j-1) *             &
     &                                  sec_latitude(i,j) *             &
     &                                  recip_r_squared_delz(i,j,k)
              End Do
            endIf !at_extremity(PNorth)
          endIf ! model_domain = mt_global and fld_type = fld_type_v

        else  ! k > metric_layer
!   No mask or metric terms since levels are horizontal

          Do j = j_begin-1, j_end
            Do i = 1, row_length
              temp(i,j) = (field(i,j+1,k) - field(i,j,k) ) *            &
     &                     diff_coeff_phi * cos_latitude(i,j+off_v)
            End Do
          End Do

          Do j = j_begin, j_end
            Do i = 1, row_length
              field(i,j,k) = field(i,j,k) + (temp(i,j) - temp(i,j-1)) * &
     &                                              sec_latitude(i,j)
            End Do
          End Do

          If ( model_domain == mt_global .and.                          &
     &             fld_type == fld_type_v ) Then
            If (at_extremity(PSouth) ) Then
              j = j_begin-1
              Do i = 1, row_length
                field(i,j,k) = field(i,j,k) + temp(i,j) *               &
     &                                           sec_latitude(i,j)
              End Do
            endIf !at_extremity(PSouth)
            If (at_extremity(PNorth) ) Then
              j = j_end + 1
              Do i = 1, row_length
                field(i,j,k) = field(i,j,k) - temp(i,j-1) *             &
     &                                        sec_latitude(i,j)
              End Do
            endIf !at_extremity(PNorth)
          endIf ! model_domain = mt_global and fld_type = fld_type_v

        endif ! k <= metric_layer

      End Do  ! k = 1, active_levels
!$OMP END PARALLEL DO

! DEPENDS ON: swap_bounds
        call Swap_Bounds(                                               &
     &                   field, row_length, in_rows, levels,            &
     &                   halo_x, halo_y, fld_type, L_vector)

      endif ! L_diff

      IF (lhook) CALL dr_hook('POFIL_NEW',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE pofil_new

