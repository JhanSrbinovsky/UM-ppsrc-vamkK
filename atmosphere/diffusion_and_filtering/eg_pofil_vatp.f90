! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine eg_pofil_vatp
!
! Purpose:
!          Filter/diffusion based on first order
!          conservative diffusion operator.
!          Applied as multiple EW 1-2-1 filter near poles.
!
! Method:
!          For any field. Pointers used to input grid
!          This version tested for v-at-the-poles/ENDGAME
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: DIFfusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!
      MODULE eg_pofil_vatp_mod
      CONTAINS
      SUBROUTINE eg_pofil_vatp(                                         &
                            field, fld_type, off_u, off_v,              &
                            base_f, base_h, base_i, base_j,             &
                            levels, model_levels, active_levels,        &
                            metric_layer,                               &
                            in_rows, cos_rows, rows, row_length,        &
                            r_levels, r_half_levels, r_midi, r_midj,    &
                            off_x, off_y, halo_x, halo_y,               &
                            halo_ri, halo_rj, halo_hi, halo_hj,         &
                            halo_ii, halo_ij, halo_ji, halo_jj,         &
                            sec_latitude, cos_latitude,                 &
                            n_procy, max_filter_rows, global_filter,    &
                            sweeps, begin, end, horizontal_level,       &
                            diff_coeff_phi, diff_coeff, L_diff,         &
                            L_vector, L_pofil_hadgem2,                  &
                            csxi2_on, csxi2_off, v_shift,L_u_ns_diff)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      LOGICAL, Intent(In) ::                                            &
        L_diff,                                                         &
                         ! true IF diffusion active
        L_vector         ! true IF vector field

      LOGICAL, Intent(In) :: L_pofil_hadgem2       ! setting for hadgem2

      LOGICAL, Intent(In) :: L_u_ns_diff ! true IF diffusing u

      INTEGER, Intent(In) ::                                            &
        max_filter_rows,                                                &
                         ! max dimension for sweeping arrays
        global_filter,                                                  &
                        ! max number of filter sweeps 0=diffusion only
        row_length,                                                     &
                          ! number of point on a row.
        rows,                                                           &
                          ! number of rows.
        in_rows,                                                        &
                          ! number of rows in field
        cos_rows,                                                       &
                          ! number of rows for cos lat
        levels,                                                         &
                          ! number of levels in field
        active_levels,                                                  &
                          ! number of levels to be filtered
        model_levels,                                                   &
                          ! number of model levels.
        metric_layer,                                                   &
                        ! uppermost non-constant layer for metric terms
        base_f,                                                         &
                      ! full levels base dim; 0 for theta, 1 for rho
        base_h,                                                         &
                      ! half levels base dim; 0 for theta, 1 for rho
        base_i,                                                         &
                      ! start level for midi levels
        base_j,                                                         &
                      ! start level for midj levels
        off_u,                                                          &
                      ! 1 for u points, 0 for v,p points
        off_v,                                                          &
                      ! 1 for v points, 0 for u,p points
        off_x,                                                          &
                      ! Size of small halo in i
        off_y,                                                          &
                      ! Size of small halo in j.
        halo_x,                                                         &
                      ! Field halo in i direction.
        halo_y,                                                         &
                      ! Field halo in j direction.
        halo_hi,                                                        &
                      ! Size of halo in i for r on half-levels
        halo_hj,                                                        &
                      ! Size of halo in j for r on half-levels
        halo_ri,                                                        &
                      ! Size of  halo_i for field r levels
        halo_rj,                                                        &
                      ! Size of  halo_j for field r levels
        halo_ii,                                                        &
                      ! Size of  halo_i for r at mid i points
        halo_ij,                                                        &
                      ! Size of  halo_j for r at mid i points
        halo_ji,                                                        &
                      ! Size of  halo_i for r at mid j points
        halo_jj,                                                        &
                      ! Size of  halo_j for r at mid j points
        fld_type,                                                       &
                      ! field type (p=1, u=2 or v=3)
        n_procy    ! Number of processors in latitude

      INTEGER, Intent(In) ::                                            &
        sweeps(max_filter_rows),                                        &
        begin(0:max_filter_rows),                                       &
        end(0:max_filter_rows),                                         &
        horizontal_level   ! level at which steep slope test no
!                               ! longer operates

!       REAL, Intent(In) ::  diff_coeff_phi   ! NS diffusion coefficient

      REAL, Intent(In) ::                                               &
        r_levels (1-halo_ri:row_length+halo_ri,                         &
                  1-halo_rj:in_rows+halo_rj, base_f:model_levels),      &
        r_half_levels (1-halo_hi:row_length+halo_hi,                    &
                       1-halo_hj:in_rows+halo_hj, base_h:model_levels), &
        r_midi (1-halo_ii:row_length+halo_ii,                           &
                1-halo_ij:in_rows+halo_ij, base_i:model_levels),        &
        r_midj (1-halo_ji:row_length+halo_ji,                           &
                1-halo_jj:cos_rows+halo_jj, base_j:model_levels)

      REAL, Intent(In) ::                                               &
        sec_latitude (1-off_x:row_length+off_x,                         &
                      1-off_y:in_rows+off_y),                           &
        cos_latitude (1-off_x:row_length+off_x,                         &
                      1-off_y:cos_rows+off_y),                          &
!   EW diffusion coefficient
        diff_coeff    (1-off_x:row_length+off_x,                        &
                       1-off_y:in_rows+off_y),                          &
! NS diffusion coefficient
        diff_coeff_phi(1-off_x:row_length+off_x,                        &
                       1-off_y:in_rows+off_y)


      INTEGER, Intent(In) :: v_shift
      REAL, Intent(In)    :: csxi2_on(v_shift-halo_rj:                  &
                                      in_rows-1+v_shift+halo_rj),       &
                             csxi2_off(1-v_shift-halo_rj:               &
                                       cos_rows-v_shift+halo_rj)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      REAL, Intent(InOut) ::                                            &
        field (1-halo_x:row_length+halo_x, 1-halo_y:in_rows+halo_y,     &
               levels)

! Local Variables.

      INTEGER ::                                                        &
        i, j, k,                                                        &
                          ! Loop indices
        k_start,                                                        &
                          ! Loop bound
        j_begin, j_end,                                                 &
                          ! Loop bounds
        j_store,                                                        &
        level_flat,                                                     &
        i_filter,                                                       &
        i_sweep,                                                        &
        num_pass

      LOGICAL :: L_cycle, L_combine

! Local arrays
      REAL  :: csphi_on(v_shift-halo_rj:in_rows-1+v_shift+halo_rj),     &
               csphi_off(1-v_shift-halo_rj:cos_rows-v_shift+halo_rj),   &
               u_fac(v_shift-halo_rj:in_rows-1+v_shift+halo_rj)

      REAL                                                              &
        delta_z(1-off_x:row_length+off_x, 1-off_y:in_rows+off_y,        &
                                                     metric_layer ),    &
        recip_r_squared_delz( 1-off_x:row_length+off_x,                 &
                              1-off_y:in_rows+off_y, metric_layer ),    &
!      &, temp( 1-off_x:row_length, 1-off_y:in_rows)                    &
! ENDGame fix 
        temp( 1-off_x:row_length, 0-off_y:in_rows),                     &
        mask_i( 1-off_x:row_length, in_rows, metric_layer ),            &
        mask_j( row_length, 1-off_y:in_rows, metric_layer )

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Set metric terms and calculate D(r)/D(eta)
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('EG_POFIL_VATP',zhook_in,zhook_handle)

!   No metric terms needed above first constant rho-level
!   since they cancel on constant r-surfaces
      DO k = 1, metric_layer
        DO j = 1-off_y, in_rows+off_y
          DO i = 1-off_x, row_length+off_x
            delta_z(i,j,k) = r_half_levels(i,j,k+base_h) -            &
                             r_half_levels(i,j,k+base_h-1)
            recip_r_squared_delz(i,j,k) = 1.0 / ( r_levels(i,j,k) *   &
                                                  r_levels(i,j,k) *   &
                                                   delta_z(i,j,k) )
          END DO
        END DO
      END DO   !  k = 1, metric_layer

! ----------------------------------------------------------------------
! Section 3.   Switch off diffusion at steep slopes
! ----------------------------------------------------------------------
!  level_flat is the uppermost non-flat level and as such
!  DOes not need testing for a sloping surface.

      DO k = 1, metric_layer
        DO j = 1, in_rows
          DO i = 1-off_x, row_length
            mask_i(i,j,k) = 1.0
          END DO
        END DO
        DO j = 0, in_rows
          DO i = 1, row_length
            mask_j(i,j,k) = 1.0
          END DO
        END DO
      END DO  ! 1, metric_layer

      IF (.NOT. L_pofil_hadgem2) THEN
!   turn off the check below for Hadgem2 runs

        IF( horizontal_level > 0 ) THEN

          level_flat = min(horizontal_level, metric_layer)

          k_start = 1
          IF ( base_f == 1 ) THEN
! for rho levels, bottom level can only be tested against surface
            DO j = 1, in_rows
              DO i = 1-off_x, row_length
                IF( r_levels(i,j,1) < r_half_levels(i+1,j,0) .OR.       &
                    r_levels(i+1,j,1) < r_half_levels(i,j,0) )          &
                mask_i(i,j,1) = 0.0
              END DO
            END DO
            DO j = 0, in_rows
              DO i = 1, row_length
                IF( r_levels(i,j,1) < r_half_levels(i,j+1,0) .OR.       &
                    r_levels(i,j+1,1) < r_half_levels(i,j,0) )          &
                mask_j(i,j,1) = 0.0
              END DO
            END DO
            k_start = 2
          END IF !  base_f == 1

          DO k = k_start, level_flat
            DO j = 1, in_rows
              DO i = 1-off_x, row_length
                IF( r_levels(i,j,k) < r_levels(i+1,j,k-1) .OR.          &
                    r_levels(i+1,j,k) < r_levels(i,j,k-1) )             &
                mask_i(i,j,k) = 0.0
              END DO
            END DO
            DO j = 0, in_rows
              DO i = 1, row_length
                IF( r_levels(i,j,k) < r_levels(i,j+1,k-1) .OR.          &
                    r_levels(i,j+1,k) < r_levels(i,j,k-1) )             &
                mask_j(i,j,k) = 0.0
              END DO
            END DO
          END DO  ! k = k_start, level_flat

        END IF ! horizontal_level > 0

      END IF ! L_pofil_hadgem2

! ----------------------------------------------------------------------
! Section 4.0  Filtering and diffusion
! ----------------------------------------------------------------------

      L_combine = .false.
      num_pass = global_filter
      j_begin = begin(1)
      j_end = end(1)
! Combine EW diffusion with first sweep of polar filter
      IF ( L_diff ) THEN
        L_combine = .true.
        IF( global_filter < 1 ) num_pass = 1
        IF ( j_begin < 0 ) THEN
          j_begin = begin(0)
          j_end = end(0)
        ELSEIF ( j_end < end(0) ) THEN
          j_end = end(0)
        ELSEIF ( j_begin > begin(0) ) THEN
          j_begin = begin(0)
        END IF !  j_begin < 0
      END IF ! L_diff

      DO i_filter = 1, num_pass

        IF( i_filter > 1 ) THEN
          L_combine = .false.
          j_begin = begin(i_filter)
          j_end = end(i_filter)
        END IF ! i_filter > 1

        L_cycle = .false.
!  IF global and n_procy = 1, need to filter hemispheres separately
        IF ( n_procy == 1) L_cycle = .TRUE.

        i_sweep = 1
 
! ----------------------------------------------------------------------
! Section 4.1   EW filtering and diffusion
! ----------------------------------------------------------------------
        DO  ! Sweeping loop  i_sweep = 1, sweeps(i_filter)

          DO k = 1, active_levels
            IF ( k <= metric_layer ) THEN
!   No mask or metric terms since levels are horizontal
              DO j = j_begin, j_end
                DO i = 0, row_length
!                   temp(i,j) = ( field(i+1,j,k) * delta_z(i+1,j,k) -     &
!                                 field(i  ,j,k) * delta_z(i  ,j,k) )     &
                  temp(i,j) = ( field  (i+1,j,k) -  field(i ,j,k) )     &
                             *( delta_z(i+1,j,k) + delta_z(i,j,k) )*0.5 &
                                 *mask_i(i,j,k) * diff_coeff(i,j) *     &
                            r_midi(i+off_u,j,k) * r_midi(i+off_u,j,k)
                END DO

                DO i = 1, row_length
                  field(i,j,k) = field(i,j,k)                           &
                               + (temp(i,j) - temp(i-1,j))              &
                                * recip_r_squared_delz(i,j,k) 
                END DO
              END DO  ! j = j_begin, j_end
 
            ELSE  ! k > metric_layer
 
              DO j = j_begin, j_end
                DO i = 0, row_length
                  temp(i,j) = ( field(i+1,j,k) - field(i,j,k) ) *       &
                                diff_coeff(i,j)                    
                END DO
                DO i = 1, row_length
                  field(i,j,k) = field(i,j,k) + (temp(i,j) - temp(i-1,j)) 
                END DO
              END DO  ! j = j_begin, j_end
            END IF ! k <= metric_layer
          END DO  ! k = 1, active_levels

! IF n_procy=1, Northern hemisphere needs to be DOne after S. Hem
          IF ( L_cycle ) THEN
            j_store = j_begin
            IF ( i_sweep == 1 .AND. L_combine ) THEN
              j_begin = j_end + 1
            ELSE
              j_begin = in_rows - j_end + 1
            END IF ! i_sweep == 1 .AND. L_combine
            j_end = in_rows - j_store + 1
            L_cycle = .false.
        CYCLE
          END IF ! L_cycle)

! Reset pointers for next sweep
! either because N Hem has just been DOne or
! 1st sweep was combined filter and diffusion
          j_begin = begin(i_filter)
          j_end = end(i_filter)
          IF( n_procy == 1 ) L_cycle = .true.

! DEPENDS ON: swap_bounds
          CALL Swap_Bounds(                                             &
                           field, row_length, in_rows, levels,          &
                           halo_x, halo_y, fld_type, L_vector)

          i_sweep = i_sweep + 1
          IF ( i_sweep > sweeps(i_filter) ) EXIT
        CYCLE
        END DO ! sweeping  loop

      END DO   !  i_filter = 1, num_pass

! ----------------------------------------------------------------------
! Section 4.2   NS diffusion
! ----------------------------------------------------------------------

      IF( L_diff ) THEN
        j_begin = begin(0)
        j_end = end(0)
        IF( n_procy == 1 ) j_end = in_rows - j_begin + 1

        ! settings for diffusion of u
        IF ( L_u_ns_diff ) Then
          DO k = 1, active_levels
            DO j = j_begin-1, j_end+1
              DO i = 1, row_length
                field(i,j,k) = field(i,j,k)/csxi2_on(j+v_shift-1)
              END DO
            END DO
          END DO

          DO j = v_shift-halo_rj,in_rows-1+v_shift+halo_rj
            csphi_on(j) = csxi2_on(j)*csxi2_on(j)
            u_fac(j) = csxi2_on(j)
          ENDDO
          DO j =1-v_shift-halo_rj,cos_rows-v_shift+halo_rj
            csphi_off(j) = csxi2_off(j)*csxi2_off(j)*csxi2_off(j)
          ENDDO 
        ELSE
          DO j = v_shift-halo_rj,in_rows-1+v_shift+halo_rj
            csphi_on(j) = csxi2_on(j)
            u_fac(j) = 1.0
          ENDDO
          DO j =1-v_shift-halo_rj,cos_rows-v_shift+halo_rj
            csphi_off(j) = csxi2_off(j)
          ENDDO        
        END IF ! L_u_ns_diff

        DO k = 1, active_levels

          IF ( k <= metric_layer ) THEN
!   No mask or metric terms since levels are horizontal

            DO j = j_begin-1, j_end
              DO i = 1, row_length
                temp(i,j) = ( field  (i,j+1,k) - field  (i,j,k) )       &
                           *( delta_z(i,j+1,k) + delta_z(i,j,k) )*0.5   &
                             *mask_j(i,j,k) * diff_coeff_phi(i,j)       &
                             *r_midj(i,j+off_v,k) * r_midj(i,j+off_v,k) &
                             *csphi_off(j)
              END DO
            END DO

            DO j = j_begin, j_end            
              DO i = 1, row_length              
                field(i,j,k) = field(i,j,k)*u_fac(j+v_shift-1)          &
                             + (temp(i,j) - temp(i,j-1))                &
                             *1.0/csphi_on(j+v_shift-1)                 &
                             *recip_r_squared_delz(i,j,k)
              END DO
            END DO

          ELSE  ! k > metric_layer
!   No mask or metric terms since levels are horizontal

            DO j = j_begin-1, j_end
              DO i = 1, row_length
                temp(i,j) = (field(i,j+1,k) - field(i,j,k) ) *          &
                             diff_coeff_phi(i,j) * csphi_off(j)
              END DO
            END DO

            DO j = j_begin, j_end
              DO i = 1, row_length
                field(i,j,k) = field(i,j,k)*u_fac(j+v_shift-1)          & 
                             + (temp(i,j) - temp(i,j-1))                &
                               *1.0/csphi_on(j+v_shift-1)
              END DO
            END DO

          END IF ! k <= metric_layer

        END DO  ! k = 1, active_levels

! DEPENDS ON: swap_bounds
        CALL Swap_Bounds(                                               &
                         field, row_length, in_rows, levels,            &
                         halo_x, halo_y, fld_type, L_vector)

      END IF ! L_diff

      IF (lhook) CALL dr_hook('EG_POFIL_VATP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE eg_pofil_vatp
      END MODULE eg_pofil_vatp_mod
