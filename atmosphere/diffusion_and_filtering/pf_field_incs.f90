! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering

MODULE pf_field_incs
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE halo_exchange, ONLY : &
      swap_bounds,          &
      swap_bounds_NS,       &
      swap_bounds_EW
  USE UM_ParParams
  IMPLICIT NONE

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

  CONTAINS

  SUBROUTINE pfwthinc(                              &
      field, r_theta_levels, r_rho_levels,          &
      off_x, off_y, halo_i, halo_j,                 &
      sec_theta_latitude, cos_v_latitude,           &
      me, n_procy,                                  &
      rows, n_rows, row_length, model_levels,       &
      max_filter_rows, u_begin, u_end,              &
      u_sweeps, global_u_filter,                    &
      levels, horizontal_level,                     &
      diff_coeff_phi, diff_coeff_u, L_diff )
    
! Purpose:
!          Filter/diffusion of a theta-type field ( theta or w )
!          Based on conservative diffusion operator.
!          Applied as multiple EW 1-2-1 filter near poles.
    

    IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

    LOGICAL, INTENT(IN) :: L_diff   ! Horizontal diffusion active

    INTEGER, INTENT(IN) ::  &
        row_length,   &             ! number of point on a row.
        rows,         &             ! number of rows.
        n_rows,       &             ! number of v-rows.
        model_levels, &             ! number of model levels.
        halo_i,       &             ! Size of halo in i direction.
        halo_j,       &             ! Size of halo in j direction.
        off_x,        &             ! Size of small halo in i
        off_y,        &             ! Size of small halo in j.
        n_procy,      &             ! Number of processors in latitude
        me        ! processor id

    INTEGER, INTENT(IN) ::          &
        max_filter_rows,            &!  array dimension for u_begin etc.
        u_begin(0:max_filter_rows), &
        u_end(0:max_filter_rows),   &
        u_sweeps(max_filter_rows),  &
        global_u_filter,            &
        levels,                     &! levels active
        horizontal_level             ! level at which steep slope test no
!                                ! longer operates

    REAL, INTENT(IN) ::       &! vertical co-ordinate arrays.
        r_theta_levels (1-halo_i:row_length+halo_i,                     &
        1-halo_j:rows+halo_j, 0:model_levels),                          &
        r_rho_levels (1-halo_i:row_length+halo_i,                       &
        1-halo_j:rows+halo_j, model_levels)

    REAL, INTENT(IN) ::       &!  trigonometric functions
        sec_theta_latitude(1-off_x:row_length+off_x,1-off_y:rows+off_y),&
        cos_v_latitude(1-off_x:row_length+off_x,1-off_y:n_rows+off_y),  &
!   EW  diffusion coefficient
        diff_coeff_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

    REAL, INTENT(INOUT) ::                                              &
        field (1-off_x:row_length+off_x, 1-off_y:rows+off_y,levels)

    REAL, INTENT(IN) :: diff_coeff_phi   ! NS diffusion coefficient

! Local Variables.

    INTEGER ::                    &
        i, j, k,                  &! Loop indices
        j_begin, j_end,           &! Loop bounds
        j_store,                  &
        level_flat,               &
        i_filter,                 &
        i_sweep,                  &
        times,                    &
        num_pass

    LOGICAL :: L_cycle, L_combine

    LOGICAL :: Active ! true if we are actually participating in this operation
    !                 ! and have polar rows to work on

! Local arrays

    REAL ::                       &
        delta_z(1-off_x:row_length+off_x,1-off_y:rows+off_y, levels ), &
        recip_r_squared_delz                                           &
        (1-off_x:row_length+off_x,1-off_y:rows+off_y, levels ),        &
        temp(1-off_x:row_length, 1-off_y:rows),                        &
        mask_i(1-off_x:row_length, rows, levels),                      &
        mask_j(row_length, 1-off_y:rows, levels)

    REAL(KIND=jprb)               :: zhook_handle


! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Set up loop bounds and pointers
!             u_begin and u_end for each potential sweep are set up
!             in SETCON (as u_begin and u_end since on same latitude)
! ----------------------------------------------------------------------
! To switch off slope test completetely then
! either hard-wire level_flat = 1 or ensure that the argument
!    horizontal level <= wet_model_levels - model_levels + 1
!      level_flat = 1
    IF (lhook) CALL dr_hook('PF_FIELD_INCS:PFWTHINC', &
        zhook_in,zhook_handle)
    level_flat =  horizontal_level
    IF (level_flat  <=  0) level_flat = 1

! ----------------------------------------------------------------------
! Section 2.   Set field to be filtered and calculate delta_z
!              Diffusion coefficient = 0.25 set in mask array
! ----------------------------------------------------------------------

! calculate D(r)/D(eta) about theta levels
    DO k = 1, levels
      IF ( k < model_levels ) THEN
        DO j = 1-off_y, rows+off_y
          DO i = 1-off_x, row_length+off_x
            delta_z(i,j,k) = r_rho_levels(i,j,k+1) -                 &
                r_rho_levels(i,j,k )
            recip_r_squared_delz(i,j,k) = 1.0 /                      &
                ( r_theta_levels(i,j,k) * r_theta_levels(i,j,k) *    &
                delta_z(i,j,k) )
          END DO
        END DO
      ELSE ! k = model_levels
! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at top level and then no need to take special action in later loops
        DO j = 1-off_y, rows+off_y
          DO i = 1-off_x, row_length+off_x
            delta_z(i,j,k) = 1.0
            recip_r_squared_delz(i,j,k) = 1.0 /                      &
                ( r_theta_levels(i,j,k) * r_theta_levels(i,j,k) )
          END DO
        END DO
      END IF !  k < model_levels
    END DO   !  k = 1, levels

! ----------------------------------------------------------------------
! Section 3.   Switch off diffusion at steep slopes
! ----------------------------------------------------------------------
!  Horizontal_level should be the first flat level and as such
!  does not need testing for a sloping surface.
!   i.e if level_flat < 2   then
!       diffusion behaves as original (along eta surfaces)

    DO k = 1, levels
      DO j = 1, rows
        DO i = 0, row_length
          mask_i(i,j,k) = 1.0
        END DO
      END DO
      DO j = 0, rows
        DO i = 1, row_length
          mask_j(i,j,k) = 1.0
        END DO
      END DO
    END DO ! k = 1, levels

    DO k = 1, level_flat - 1
      DO j = 1, rows
        DO i = 0, row_length
          IF ( r_theta_levels(i,j,k) < r_theta_levels(i+1,j,k-1) .OR.  &
              r_theta_levels(i+1,j,k) < r_theta_levels(i,j,k-1) )      &
              mask_i(i,j,k) = 0.0
        END DO
      END DO
      DO j = 0, rows
        DO i = 1, row_length
          IF ( r_theta_levels(i,j,k) < r_theta_levels(i,j+1,k-1) .OR.  &
              r_theta_levels(i,j+1,k) < r_theta_levels(i,j,k-1) )      &
              mask_j(i,j,k) = 0.0
        END DO
      END DO
    END DO    ! k = 1, level_flat - 1

! ----------------------------------------------------------------------
! Section 4.0  Polar Filtering and diffusion
!              j loops will no-op for filtering
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
!        write(6,*)'pofil_wth_incs global_u_filter = ', global_u_filter

    L_combine = .FALSE.
    num_pass = global_u_filter
    j_begin = u_begin(1)
    j_end = u_end(1)

! Combine EW diffusion with first sweep of polar filter
    IF ( L_diff ) THEN
      L_combine = .TRUE.
      IF ( global_u_filter < 1 ) num_pass = 1
      IF ( j_begin < 0 ) THEN
        j_begin = u_begin(0)
        j_end   = u_end(0)
      ELSE IF ( j_end < u_end(0) ) THEN
        j_end   = u_end(0)
      ELSE IF ( j_begin > u_begin(0) ) THEN
        j_begin = u_begin(0)
      END IF !  j_begin < 0
    END IF ! L_diff

    DO i_filter = 1, num_pass

      IF ( i_filter > 1 ) THEN
        L_combine = .FALSE.
        j_begin = u_begin(i_filter)
        j_end = u_end(i_filter)
      END IF !i_filter > 1


      L_cycle = .FALSE.
!  if n_procy = 1, need to filter hemispheres separately
      IF ( n_procy == 1 ) L_cycle = .TRUE.

      i_sweep = 1
      DO 

        ! We will only halo exchange if we participate.
        active=.FALSE.
        IF (j_end-j_begin >= 0) active=.TRUE.

! ----------------------------------------------------------------------
! Section 4.1   EW filtering and diffusion
! ----------------------------------------------------------------------
        DO k = 1, levels
! There are 2 factors of 1/4 multiplied to give 1/16
          DO j = j_begin, j_end
            DO i = 1-off_x, row_length
              temp(i,j) = (field(i+1,j,k) * delta_z(i+1,j,k) -              &
                  field(i,j,k) * delta_z(i,j,k) ) *                         &
                  mask_i(i,j,k) * diff_coeff_u(i,j) *                       &
                  0.25 * (r_theta_levels(i,j,k) + r_theta_levels(i+1,j,k)) *&
                  (r_theta_levels(i,j,k) + r_theta_levels(i+1,j,k))
            END DO

            DO i = 1, row_length
              field(i,j,k) = field(i,j,k) +                                 &
                  recip_r_squared_delz(i,j,k) * ( temp(i,j) - temp(i-1,j) )
            END DO

          END DO  ! j = j_begin, j_end
        END DO  ! k = 1, levels

! If n_procy=1, Northern hemisphere needs to be done after S. Hem
        IF ( L_cycle ) THEN
          j_store = j_begin
          IF ( i_sweep == 1 .AND. L_combine ) THEN
            j_begin = j_end + 1
          ELSE
            j_begin = rows - j_end + 1
          END IF ! i_sweep == 1 .AND. L_combine
          j_end = rows - j_store + 1
          L_cycle = .FALSE.
          CYCLE
        END IF ! L_cycle)

! Reset pointers for next sweep
! either because N Hem has just been done or
! 1st sweep was combined filter and diffusion
        j_begin = u_begin(i_filter)
        j_end = u_end(i_filter)
        IF ( n_procy == 1 ) L_cycle = .TRUE.

        IF (active)                                  & 
            CALL Swap_Bounds_EW(                     &
            field, row_length, rows, levels,         &
            off_x, off_y)

        i_sweep = i_sweep + 1
        IF ( i_sweep > u_sweeps(i_filter) ) EXIT
        CYCLE
      END DO ! sweeping  loop  i_sweep = 1, u_sweeps(i_filter)        
    END DO   !  i_filter = 1, num_pass

! ----------------------------------------------------------------------
! Section 4.2  NS diffusion
! ----------------------------------------------------------------------



    IF ( L_diff ) THEN
! 'Field' isn't completely boundary swapped after the above
!  we need j+1 (North) elements next
      CALL Swap_Bounds(                          &
          field, row_length, rows, levels,          &
          off_x, off_y, fld_type_p, .FALSE.)

      j_begin = u_begin(0)
      j_end = u_end(0)
      IF ( n_procy == 1 ) j_end = rows - j_begin + 1
!          write(6,*)' NS diffusion  diff_coeff_phi = ', diff_coeff_phi
!               ,' j_begin = ', j_begin,' j_end = ', j_end

      DO k = 1, levels

        DO j = j_begin - 1, j_end
          DO i = 1, row_length
            temp(i,j) = (field(i,j+1,k) * delta_z(i,j+1,k) -               &
                field(i,j,  k) * delta_z(i,j,  k) ) *                      &
                mask_j(i,j,k) * diff_coeff_phi *                           &
                0.25 * (r_theta_levels(i,j,k) + r_theta_levels(i,j+1,k)) * &
                (r_theta_levels(i,j,k) + r_theta_levels(i,j+1,k)) *        &
                cos_v_latitude(i,j)
          END DO
        END DO

        DO j = j_begin, j_end
          DO i = 1, row_length
            field(i,j,k) = field(i,j,k) + sec_theta_latitude(i,j) *   &
                recip_r_squared_delz(i,j,k) *                         &
                (temp(i,j) - temp(i,j-1))
          END DO
        END DO  ! j = j_begin, j_end

      END DO  ! k = 1, levels
    END IF ! L_diff
! Do a full halo exchange to the field is halo swaped at exit. 
    CALL Swap_Bounds(                                  &
        field, row_length, rows, levels,               &
        off_x, off_y, fld_type_p, .FALSE.)



    IF (lhook) CALL dr_hook('PF_FIELD_INCS:PFWTHINC', &
        zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE pfwthinc


  SUBROUTINE pfuinc(                                  &
      R_u, r_at_u,                                    &
      r_theta_levels, r_rho_levels,                   &
      off_x, off_y, halo_i, halo_j,                   &
      sec_theta_latitude, cos_v_latitude,             &
      me, n_procy,                                    &
      rows, n_rows, row_length, model_levels,         &
      max_filter_rows, u_begin, u_end,                &
      u_sweeps, global_u_filter,                      &
      horizontal_level,                               &
      diff_coeff_phi, diff_coeff_u, L_diff )

! Purpose:
!          Filter/diffusion of a u-type field
!          Based on conservative diffusion operator.
!          Applied as multiple EW 1-2-1 filter near poles.
!


    IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

    LOGICAL                                                           &
        L_diff   ! Horizontal diffusion of increments
    
    INTEGER            &
        row_length,    &! number of point on a row.
        rows,          &! number of rows.
        n_rows,        &! number of v-rows.
        model_levels,  &! number of model levels.
        halo_i,        &! Size of halo in i direction.
        halo_j,        &! Size of halo in j direction.
        off_x,         &! Size of small halo in i
        off_y,         &! Size of small halo in j.
        n_procy,       &! Number of processors in latitude
        me              ! processor id

    INTEGER                        &
        max_filter_rows,           &!  array dimension for u_begin etc.
        u_begin(0:max_filter_rows),&
        u_end(0:max_filter_rows),  &
        u_sweeps(max_filter_rows), &
        global_u_filter,           &
        horizontal_level            ! level at which steep slope test no
    !                               ! longer operates

    REAL                           &
        diff_coeff_phi              ! NS diffusion coefficient

    REAL                       &! vertical co-ordinate arrays.
        r_at_u (1-halo_i:row_length+halo_i,          &
        1-halo_j:rows+halo_j, model_levels),         &
        r_theta_levels (1-halo_i:row_length+halo_i,  &
        1-halo_j:rows+halo_j, 0:model_levels),       &
        r_rho_levels (1-halo_i:row_length+halo_i,    &
        1-halo_j:rows+halo_j, model_levels)

    REAL                       &!  trigonometric functions
        sec_theta_latitude(1-off_x:row_length+off_x,  &
        1-off_y:rows+off_y),                          &
        cos_v_latitude(1-off_x:row_length+off_x,      &
        1-off_y:n_rows+off_y),                        &
    !                           ! EW diffusion coefficient on u rows
        diff_coeff_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.
    
    REAL                                                              &
        R_u (1-off_x:row_length+off_x,                                &
        1-off_y:rows+off_y, model_levels)

! Local Variables.

    INTEGER                                                           &
        i, j, k,       &! Loop indices
        j_begin, j_end,&! Loop bounds
        j_store,       &
        level_flat,    &
        i_filter,      &
        i_sweep,       &
        num_pass

    LOGICAL                                           &
        L_cycle,                                      &
        L_combine

! Local arrays

    REAL                                              &
        r_theta_at_u(1-off_x:row_length+off_x,        &
        1-off_y:rows+off_y, 0:model_levels),          &
        delta_z(1-off_x:row_length+off_x,             &
        1-off_y:rows+off_y, model_levels ),           &
        recip_r_squared_delz(1-off_x:row_length+off_x,&
        1-off_y:rows+off_y, model_levels ),           &
        temp(row_length+off_x, 1-off_y:rows),         &
        mask_i(row_length+off_x, rows, model_levels), &
        mask_j(row_length, 1-off_y:rows, model_levels)

    LOGICAL :: Active ! true if we are actually participating in this operation
    !                 ! and have polar rows to work on
    
    REAL(KIND=jprb)               :: zhook_handle
    
    

! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Set up loop bounds and pointers
!             u_begin and u_end for each potential sweep are set up
! ----------------------------------------------------------------------

    IF (lhook) CALL dr_hook('PF_FIELD_INCS:PFUINC', &
        zhook_in,zhook_handle)
! To switch off slope test completetely then
! either hard-wire level_flat = 1 or ensure that the argument
!    horizontal level <= wet_model_levels - - model_levels + 1
!      level_flat = 1
    level_flat =  horizontal_level
    IF (level_flat <= 0) level_flat = 1

! ----------------------------------------------------------------------
! Section 2.   Set r values and calculate delta_z
! ----------------------------------------------------------------------

    DO k = 0, model_levels
      DO j = 1-off_y, rows+off_y
        DO i = 1-off_x, row_length+off_x
          r_theta_at_u(i,j,k) = .5 * (r_theta_levels(i+1,j,k) +       &
              r_theta_levels(i  ,j,k) )
        END DO
      END DO
    END DO ! k = 0, model_levels

! calculate D(r)/D(eta) about theta levels
    DO k = 1, model_levels - 1
      DO j = 1-off_y, rows+off_y
        DO i = 1-off_x, row_length+off_x
          delta_z(i,j,k) = r_theta_at_u(i,j,k) - r_theta_at_u(i,j,k-1)
          recip_r_squared_delz(i,j,k) = 1.0 /          &
              ( r_at_u(i,j,k) * r_at_u(i,j,k) *        &
              delta_z(i,j,k) )
        END DO
      END DO
    END DO   !  k = 1, model_levels - 1
! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at top level and then no need to take special action in later loops
    k = model_levels
    DO j = 1-off_y, rows+off_y
      DO i = 1-off_x, row_length+off_x
        delta_z(i,j,k) = 1.0
        recip_r_squared_delz(i,j,k) = 1.0 /            &
            ( r_at_u(i,j,k) * r_at_u(i,j,k) )
      END DO
    END DO

! ----------------------------------------------------------------------
! Section 3.   Switch off theta diffusion at steep slopes
! ----------------------------------------------------------------------
!  Horizontal_level should be the first flat level and as such
!  does not need testing for a sloping surface.
!   i.e if level_flat < 2   then
!       diffusion behaves as original (along eta surfaces)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length+1
          mask_i(i,j,k) = 1.0
        END DO
      END DO
      DO j = 0, rows
        DO i = 1, row_length
          mask_j(i,j,k) = 1.0
        END DO
      END DO
    END DO    ! k = 1, model_levels

    DO j = 1, rows
      DO i = 1, row_length+1
        if( r_at_u(i-1,j,1)   <  r_theta_at_u(i,j,0) .OR.           &
            r_at_u(i,j,1) <  r_theta_at_u(i-1,j,0) )                &
            mask_i(i,j,1) = 0.0
      END DO
    END DO
    DO j = 0, rows
      DO i = 1, row_length
        if( r_at_u(i,j,1) < r_theta_at_u(i,j+1,0) .OR.              &
            r_at_u(i,j+1,1) < r_theta_at_u(i,j,0) )                 &
            mask_j(i,j,1) = 0.0
      END DO
    END DO
    DO k = 2, level_flat - 1
      DO j = 1, rows
        DO i = 1, row_length + 1
          if( r_at_u(i-1,j,k) < r_at_u(i,j,k-1) .OR.                  &
              r_at_u(i,j,k) < r_at_u(i-1,j,k-1) )                     &
              mask_i(i,j,k) = 0.0
        END DO
      END DO
      DO j = 0, rows
        DO i = 1, row_length
          if( r_at_u(i,j,k) < r_at_u(i,j+1,k-1) .OR.                  &
              r_at_u(i,j+1,k) < r_at_u(i,j,k-1) )                     &
              mask_j(i,j,k) = 0.0
        END DO
      END DO
    END DO    ! k = 2, level_flat - 1

! ----------------------------------------------------------------------
! Section 4.0  Polar filtering and diffusion
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
!        write(6,*)'pofil_u_incs global_u_filter = ', global_u_filter

    L_combine = .FALSE.
    num_pass = global_u_filter
    j_begin = u_begin(1)
    j_end = u_end(1)
! Combine EW diffusion with first sweep of polar filter
    IF ( L_diff ) THEN
      L_combine = .TRUE.
      IF ( global_u_filter < 1 ) num_pass = 1
      IF ( j_begin < 0 ) THEN
        j_begin = u_begin(0)
        j_end = u_end(0)
      ELSE IF ( j_end < u_end(0) ) THEN
        j_end = u_end(0)
      ELSE IF ( j_begin > u_begin(0) ) THEN
        j_begin = u_begin(0)
      END IF !  j_begin < 0
    END IF ! L_diff

    DO i_filter = 1, num_pass
      
      IF ( i_filter > 1 ) THEN
        L_combine = .FALSE.
        j_begin = u_begin(i_filter)
        j_end = u_end(i_filter)
      END IF !i_filter > 1



      L_cycle = .FALSE.
!  if n_procy = 1, need to filter hemispheres separately
      IF ( n_procy == 1 ) L_cycle = .TRUE.
      
      i_sweep = 1
      DO  

        ! We only need to halo exchange if our j range is finite
        active=.FALSE.
        IF (j_end-j_begin >= 0) active=.TRUE.

! ----------------------------------------------------------------------
! Section 4.1   EW filtering and diffusion
! ----------------------------------------------------------------------
        DO k = 1, model_levels

          DO j = j_begin, j_end

            DO i = 1, row_length+1
              temp(i,j) = (R_u(i  ,j,k) * delta_z(i  ,j,k) -  &
                  R_u(i-1,j,k) * delta_z(i-1,j,k) ) *         &
                  mask_i(i,j,k) * diff_coeff_u(i,j) *         &
                  r_rho_levels(i,j,k) * r_rho_levels(i,j,k)
            END DO

            DO i = 1, row_length
              R_u(i,j,k) = R_u(i,j,k) + ( temp(i+1,j) - temp(i,j) ) * &
                  recip_r_squared_delz(i,j,k)
            END DO

          END DO  ! j = j_begin, j_end

        END DO  ! k = 1, model_levels

! If n_procy=1, Northern hemisphere needs to be done after S. Hem
        IF ( L_cycle ) THEN
          j_store = j_begin
          IF ( i_sweep == 1 .AND. L_combine ) THEN
            j_begin = j_end + 1
          ELSE
            j_begin = rows - j_end + 1
          END IF ! i_sweep == 1 .AND. L_combine
          j_end = rows - j_store + 1
          L_cycle = .FALSE.
          CYCLE
        END IF ! L_cycle)

! Reset pointers for next sweep
! either because N Hem has just been done or
! 1st sweep was combined filter and diffusion
        j_begin = u_begin(i_filter)
        j_end = u_end(i_filter)
        IF ( n_procy == 1 ) L_cycle = .TRUE.

        IF (Active)                                     &
            CALL Swap_Bounds_EW(                         &
            R_u, row_length, rows, model_levels,        &
            off_x, off_y)

        i_sweep = i_sweep + 1
        IF ( i_sweep > u_sweeps(i_filter) ) EXIT
        CYCLE
      END DO ! sweeping  loop  i_sweep = 1, u_sweeps(i_filter)

    END DO   !  i_filter = 1, num_pass


! ----------------------------------------------------------------------
! Section 4.2  NS diffusion
! ----------------------------------------------------------------------

    IF ( L_diff ) THEN
      CALL Swap_Bounds  (                            &
          R_u, row_length, rows, model_levels,          &
          off_x, off_y, fld_type_u, .TRUE.)

      j_begin = u_begin(0)
      j_end = u_end(0)
      IF ( n_procy == 1 ) j_end = rows - j_begin + 1
!          write(6,*)' NS diffusion  diff_coeff_phi = ', diff_coeff_phi
!     &          ,' j_begin = ', j_begin,' j_end = ', j_end

      DO k = 1, model_levels

        DO j = j_begin - 1, j_end
          DO i = 1, row_length
            temp(i,j) = (R_u(i,j+1,k) * delta_z(i,j+1,k) -    &
                R_u(i,j,  k) * delta_z(i,j, k) ) *            &
                mask_j(i,j,k) * diff_coeff_phi *              &
                0.25 * ( r_at_u(i,j,k) + r_at_u(i,j+1,k) ) *  &
                ( r_at_u(i,j,k) + r_at_u(i,j+1,k) ) *         &
                cos_v_latitude(i,j)
          END DO
        END DO
        DO j = j_begin, j_end
          DO i = 1, row_length
            R_u(i,j,k) = R_u(i,j,k) + recip_r_squared_delz(i,j,k) * &
                ( temp(i,j) - temp(i,j-1) ) *                       &
                sec_theta_latitude(i,j)
          END DO
        END DO

      END DO  ! k = 1, model_levels
    END IF ! L_diff

! Ensure fully swapped at exit
    CALL Swap_Bounds(                           &
        R_u, row_length, rows, model_levels,    &
        off_x, off_y, fld_type_u, .TRUE.)
    
    IF (lhook) CALL dr_hook('PF_FIELD_INCS:PFUINC', &
        zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE pfuinc


  SUBROUTINE pfvinc(                             &               
      R_v, r_at_v, r_theta_levels, r_rho_levels, &
      off_x, off_y, halo_i, halo_j,              &
      sec_v_latitude, cos_theta_latitude,        &
      me, n_procy, at_extremity,                 &
      rows, n_rows, row_length, model_levels,    &
      max_filter_rows, v_begin, v_end,           &
      v_sweeps, global_v_filter,                 &
      horizontal_level, model_domain,            &
      diff_coeff_phi, diff_coeff_v, L_diff )

! Purpose:
!          Filter/diffusion of a v-type field
!          Based on conservative diffusion operator.
!          Applied as multiple EW 1-2-1 filter near poles.


    IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

    LOGICAL, INTENT(IN) :: &
        at_extremity(4),   & ! Indicates if this processor is at north,
!                    ! south, east or west of the processor grid
        L_diff

    INTEGER, INTENT(IN) ::                                             &
        row_length,   &             ! number of point on a row.
        rows,         &             ! number of rows.
        n_rows,       &             ! number of v-rows.
        model_levels, &             ! number of model levels.
        halo_i,       &             ! Size of halo in i direction.
        halo_j,       &             ! Size of halo in j direction.
        off_x,        &             ! Size of small halo in i
        off_y,        &             ! Size of small halo in j.
        model_domain, &
        n_procy,      &             ! Number of processors in latitude
        me        ! processor id   

    INTEGER, INTENT(IN) ::          &
        max_filter_rows,            &!  array dimension for u_begin etc.
        v_begin(0:max_filter_rows), &
        v_end(0:max_filter_rows),   &
        v_sweeps(max_filter_rows),  &
        global_v_filter,            &
        horizontal_level             ! level at which steep slope test no
!                                ! longer operates

    REAL, INTENT(IN) :: diff_coeff_phi   ! NS diffusion coefficient

    REAL, INTENT(IN) ::       &! vertical co-ordinate arrays.
        r_at_v (1-halo_i:row_length+halo_i,                             &
        1-halo_j:n_rows+halo_j, model_levels),                          &
        r_theta_levels (1-halo_i:row_length+halo_i,                     &
        1-halo_j:rows+halo_j, 0:model_levels),                          &
        r_rho_levels (1-halo_i:row_length+halo_i,                       &
        1-halo_j:rows+halo_j, model_levels)

    REAL, INTENT(IN) ::       &!  trigonometric functions
        sec_v_latitude(1-off_x:row_length+off_x,1-off_y:n_rows+off_y),&
        cos_theta_latitude(1-off_x:row_length+off_x,1-off_y:rows+off_y),    &
!   EW  diffusion coefficient
        diff_coeff_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y)

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

    REAL, INTENT(INOUT) :: &
        R_v(1-off_x:row_length+off_x,1-off_y:n_rows+off_y, model_levels)

! Local Variables.

    INTEGER ::                    &
        i, j, k,                  &! Loop indices
        j_begin, j_end,           &! Loop bounds
        j_store,                  &
        level_flat,               &
        i_filter,                 &
        i_sweep,                  &
        num_pass

    LOGICAL :: L_cycle, L_combine

    LOGICAL :: Active ! true if we are actually participating in this operation
!                 ! and have polar rows to work on

! Local arrays

    REAL ::                                               &
        r_theta_at_v(1-off_x:row_length+off_x,            &
        1-off_y:n_rows+off_y, 0:model_levels),            &
        delta_z(1-off_x:row_length+off_x,                 &
        1-off_y:n_rows+off_y, model_levels ),             &
        recip_r_squared_delz(1-off_x:row_length+off_x,    &
        1-off_y:n_rows+off_y, model_levels ),             &
        temp(1-off_x:row_length, 1-off_y:n_rows),         &
        mask_i(1-off_x:row_length, n_rows, model_levels), &
        mask_j(row_length, 1-off_y:n_rows, model_levels)

    REAL(KIND=jprb)               :: zhook_handle


! No External Routines:

! ----------------------------------------------------------------------
! Section 1.   Set up loop bounds and pointers
!            v_begin and v_end for each potential sweep are set SETCON
! ----------------------------------------------------------------------
! To switch off slope test completetely THEN
! either hard-wire level_flat = 1 or ensure that the argument
!    horizontal level <= wet_model_levels - - model_levels + 1
!      level_flat = 1
    IF (lhook) CALL dr_hook('PF_FIELD_INCS:PFVINC', &
        zhook_in,zhook_handle)
    level_flat =  horizontal_level
    IF (level_flat  <=  0) level_flat = 1

! ----------------------------------------------------------------------
! Section 2.   Set r field  and calculate delta_z
! ----------------------------------------------------------------------

    DO k = 0, model_levels
      DO j = 1-off_y, n_rows+off_y
        DO i = 1-off_x, row_length+off_x
          r_theta_at_v(i,j,k) = .5 * (r_theta_levels(i,j+1,k) +       &
              r_theta_levels(i,j  ,k) )
        END DO
      END DO
    END DO ! k = 0, model_levels

! calculate D(r)/D(eta) about theta levels
    DO k = 1, model_levels - 1
      DO j = 1-off_y, n_rows+off_y
        DO i = 1-off_x, row_length+off_x
          delta_z(i,j,k) = r_theta_at_v(i,j,k) -     &
              r_theta_at_v(i,j,k-1)
          recip_r_squared_delz(i,j,k) = 1.0 /        &
              ( r_at_v(i,j,k) * r_at_v(i,j,k) *      &
              delta_z(i,j,k) )
        END DO
      END DO
    END DO   !  k = 1, model_levels - 1
! By setting dr/deta=1 we effectively cancel out constant dr/deta
! at top level and THEN no need to take special action in later loops
    k = model_levels
    DO j = 1-off_y, n_rows+off_y
      DO i = 1-off_x, row_length+off_x
        delta_z(i,j,k) = 1.0
        recip_r_squared_delz(i,j,k) = 1.0 /          &
            ( r_at_v(i,j,k) * r_at_v(i,j,k) )
      END DO
    END DO

! ----------------------------------------------------------------------
! Section 3.   Switch off theta diffusion at steep slopes
! ----------------------------------------------------------------------
!  Horizontal_level should be the first flat level and as such
!  does not need testing for a sloping surface.
!   i.e if level_flat < 2   THEN
!       diffusion behaves as original (along eta surfaces)

    DO k = 1, model_levels
      DO j = 1, n_rows
        DO i = 1-off_x, row_length
          mask_i(i,j,k) = 1.0
        END DO
      END DO
      DO j = 0, n_rows
        DO i = 1, row_length
          mask_j(i,j,k) = 1.0
        END DO
      END DO
    END DO  ! 1, model_levels

    DO j = 1, n_rows
      DO i = 1-off_x, row_length
        IF ( r_at_v(i,j,1) < r_theta_at_v(i+1,j,0) .OR. &
            r_at_v(i+1,j,1) < r_theta_at_v(i,j,0) )     &
            mask_i(i,j,1) = 0.0
      END DO
    END DO
    DO j = 0, n_rows
      DO i = 1, row_length
        IF ( r_at_v(i,j,1) < r_theta_at_v(i,j+1,0) .OR. &
            r_at_v(i,j+1,1) < r_theta_at_v(i,j,0) )     &
            mask_j(i,j,1) = 0.0
      END DO
    END DO
    DO k = 2, level_flat - 1
      DO j = 1, n_rows
        DO i = 1-off_x, row_length
          IF ( r_at_v(i,j,k) < r_at_v(i+1,j,k-1) .OR.   &
              r_at_v(i+1,j,k) < r_at_v(i,j,k-1) )       &
              mask_i(i,j,k) = 0.0
        END DO
      END DO
      DO j = 0, n_rows
        DO i = 1, row_length
          IF ( r_at_v(i,j,k) < r_at_v(i,j+1,k-1) .OR.   &
              r_at_v(i,j+1,k) < r_at_v(i,j,k-1) )       &
              mask_j(i,j,k) = 0.0
        END DO
      END DO
    END DO  ! k = 2, level_flat - 1

! ----------------------------------------------------------------------
! Section 4.0  Polar Filtering and diffusion
! ----------------------------------------------------------------------

!        write(6,*)'pofil_v_incs global_v_filter = ', global_v_filter

    L_combine = .FALSE.
    num_pass = global_v_filter
    j_begin = v_begin(1)
    j_end = v_end(1)
! Combine EW diffusion with first sweep of polar filter
    IF ( L_diff ) THEN
      L_combine = .TRUE.
      IF ( global_v_filter < 1 ) num_pass = 1
      IF ( j_begin < 0 ) THEN
        j_begin = v_begin(0)
        j_end = v_end(0)
      ELSE IF ( j_end < v_end(0) ) THEN
        j_end = v_end(0)
      ELSE IF ( j_begin > v_begin(0) ) THEN
        j_begin = v_begin(0)
      END IF !  j_begin < 0
    END IF ! L_diff

    DO i_filter = 1, num_pass

      IF ( i_filter > 1 ) THEN
        L_combine = .FALSE.
        j_begin = v_begin(i_filter)
        j_end = v_end(i_filter)
      END IF !i_filter > 1


      L_cycle = .FALSE.
!  If n_procy = 1, need to filter hemispheres separately
      IF ( n_procy == 1 ) L_cycle = .TRUE.

      i_sweep = 1

      DO  
        ! We only need to halo exchange if our j range is finite
        active=.FALSE.
        IF (j_end-j_begin >= 0) active=.TRUE.

! ----------------------------------------------------------------------
! Section 4.1   EW Filtering and diffusion
! ----------------------------------------------------------------------
        DO k = 1, model_levels

          DO j = j_begin, j_end
            DO i = 1-off_x, row_length
              temp(i,j) = (R_v(i+1,j,k) * delta_z(i+1,j,k) -   &
                  R_v(i  ,j,k) * delta_z(i  ,j,k) ) *          &
                  mask_i(i,j,k) * diff_coeff_v(i,j) *          &
                  0.25 * ( r_at_v(i,j,k) + r_at_v(i+1,j,k) ) * &
                  ( r_at_v(i,j,k) + r_at_v(i+1,j,k) )
            END DO

            DO i = 1, row_length
              R_v(i,j,k) = R_v(i,j,k) + ( temp(i,j) - temp(i-1,j) ) * &
                  recip_r_squared_delz(i,j,k)
            END DO
          END DO  ! j = j_begin, j_end

        END DO  ! k = 1, model_levels

! If n_procy=1, Northern hemisphere needs to be done after S. Hem
        IF ( L_cycle ) THEN
          j_store = j_begin
          IF ( i_sweep == 1 .AND. L_combine ) THEN
            j_begin = j_end + 1
          ELSE
            j_begin = n_rows - j_end + 1
          END IF ! i_sweep == 1 .AND. L_combine
          j_end = n_rows - j_store + 1
          L_cycle = .FALSE.
          CYCLE
        END IF ! L_cycle)

! Reset pointers for next sweep
! either because N Hem has just been done or
! 1st sweep was combined filter and diffusion
        j_begin = v_begin(i_filter)
        j_end = v_end(i_filter)
        IF ( n_procy == 1 ) L_cycle = .TRUE.

        IF (active)                                   &
            CALL Swap_Bounds_EW(                      &
            R_v, row_length, n_rows, model_levels,    &
            off_x, off_y)

        i_sweep = i_sweep + 1
        IF ( i_sweep > v_sweeps(i_filter) ) EXIT
        CYCLE
      END DO ! sweeping  loop

    END DO   !  i_filter = 1, num_pass

! ----------------------------------------------------------------------
! Section 4.2   NS diffusion
! ----------------------------------------------------------------------

    IF ( L_diff ) THEN
! Because the swaps above only updated the east halo and the code 
! below needs the N halo....
      CALL Swap_Bounds(                            &
          R_v, row_length, n_rows, model_levels,      &
          off_x, off_y, fld_type_v, .TRUE.)

      j_begin = v_begin(0)
      j_end = v_end(0)
      IF ( n_procy == 1 ) j_end = n_rows - j_begin + 1
!          write(6,*)' NS diffusion  diff_coeff_phi = ', diff_coeff_phi
!          ,' j_begin = ', j_begin,' j_end = ', j_end

      DO k = 1, model_levels

        DO j = j_begin-off_y, j_end
          DO i = 1, row_length
            temp(i,j) = (R_v(i,j+1,k) * delta_z(i,j+1,k) -     &
                R_v(i,j,  k) * delta_z(i,j, k) ) *             &
                mask_j(i,j,k) * diff_coeff_phi *               &
                r_rho_levels(i,j+1,k) * r_rho_levels(i,j+1,k) *&
                cos_theta_latitude(i,j+1)
          END DO
        END DO
        DO j = j_begin, j_end
          DO i = 1, row_length
            R_v(i,j,k) = R_v(i,j,k) + sec_v_latitude(i,j) *    &
                ( temp(i,j) - temp(i,j-1) ) *          &
                recip_r_squared_delz(i,j,k)
          END DO
        END DO

        IF ( model_domain == mt_global ) THEN
          IF (at_extremity(PSouth) ) THEN
            j = j_begin-1
            DO i = 1, row_length
              R_v(i,j,k) = R_v(i,j,k) + recip_r_squared_delz(i,j,k) * &
                  temp(i,j) * sec_v_latitude(i,j)
            END DO
          END IF !at_extremity(PSouth)
          IF (at_extremity(PNorth) ) THEN
            j = j_end + 1
            DO i = 1, row_length
              R_v(i,j,k) = R_v(i,j,k) - recip_r_squared_delz(i,j,k) * &
                  temp(i,j-1) * sec_v_latitude(i,j)
            END DO
          END IF !at_extremity(PNorth)
        END IF ! model_domain == mt_global )

      END DO  ! k = 1, model_levels

    END IF ! L_diff

! Ensure the field is fully swapped on exit
    CALL Swap_Bounds(                                  &
        R_v, row_length, n_rows, model_levels,         &
        off_x, off_y, fld_type_v, .TRUE.)

    IF (lhook) CALL dr_hook('PF_FIELD_INCS:PFVINC', &
        zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE pfvinc

END MODULE pf_field_incs
