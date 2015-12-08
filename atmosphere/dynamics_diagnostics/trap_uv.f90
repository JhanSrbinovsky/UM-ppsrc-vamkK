! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine trap_uv

      SUBROUTINE trap_uv(                                               &
                         u, v, u_adv, v_adv, max_wind,                  &
                         rows, n_rows, row_length, model_levels,        &
                         off_x, off_y, halo_i, halo_j,                  &
                         trap_option )
      USE proc_info_mod

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

! Purpose:
!          Diagnostic routine for trapping when wind components
!             exceed maximum

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Diagnostics

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      INTEGER , INTENT(IN) ::                                           &
        row_length,                                                     &
                         ! number of point on a row.
        rows,                                                           &
                         ! number of rows.
        n_rows,                                                         &
                         ! number of rows.
        model_levels,                                                   &
                         ! number of model levels.
        off_x,                                                          &
                         ! Size of small halo in i
        off_y,                                                          &
                         ! Size of small halo in j.
        halo_i,                                                         &
                         ! Size of halo for wind components.
        halo_j,                                                         &
                         ! Size of halo for wind component.
        trap_option      ! 0 = reset, no prints
                         ! 1 = reset + print message
                         ! 2 = NO reset but print details max/min

      REAL , INTENT(INOUT) ::                                           &
        u (1-off_x:row_length+off_x,                                    &
           1-off_y:rows+off_y, model_levels),                           &
        v (1-off_x:row_length+off_x,                                    &
           1-off_y:n_rows+off_y, model_levels),                         &
        u_adv (1-halo_i:row_length+halo_i,                              &
               1-halo_j:rows+halo_j, model_levels),                     &
        v_adv (1-halo_i:row_length+halo_i,                              &
               1-halo_j:n_rows+halo_j, model_levels)

      REAL , INTENT(IN) :: max_wind

! Local Variables.

      INTEGER                                                           &
        i, j, k,                                                        &
                               ! Loop counters
        gi, gj,                                                         &
                               ! global pointers
        ki, ku, kv,                                                     &
                               ! pointers for arrays summed over pe's
        ic,                                                             &
        info,                                                           &
        itests,                                                         &
        i_stop_u,                                                       &
        j_start_u, j_stop_u,                                            &
                               ! Loop indices
        level_p,                                                        &
                               ! print level of max/min
        lambda_p,                                                       &
                               ! print % domain longitude of max/min
        phi_p                  ! print % domain latitude of  max/min

! Local arrays
      INTEGER                                                           &
        itest(4),                                                       &
        sumi( 8 * model_levels )

      REAL                                                              &
        max_u_pe(model_levels),                                         &
        min_u_pe(model_levels),                                         &
        max_v_pe(model_levels),                                         &
        min_v_pe(model_levels),                                         &
        max_real(2*model_levels),                                       &
        min_real(2*model_levels)

      INTEGER                                                           &
        i_max_u(model_levels),                                          &
        j_max_u(model_levels),                                          &
        i_min_u(model_levels),                                          &
        j_min_u(model_levels),                                          &
        i_max_v(model_levels),                                          &
        j_max_v(model_levels),                                          &
        i_min_v(model_levels),                                          &
        j_min_v(model_levels)

      INTEGER  :: min_indices(2), max_indices(2)

      REAL                                                              &
        recip_row_length,                                               &
        recip_rows,                                                     &
        w_p                     ! print max/min

      CHARACTER(LEN=8) North, East
      PARAMETER ( East=' % East ')
      PARAMETER ( North=' % North')
      CHARACTER(LEN=*) mperspr
      PARAMETER ( mperspr=' m/s at ')
      CHARACTER(LEN=*) prulim
      PARAMETER ( prulim=' winds u and u_adv limited to ')
      CHARACTER(LEN=*) prvlim
      PARAMETER ( prvlim=' winds v and v_adv limited to ')

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 0. Initialise values
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('TRAP_UV',zhook_in,zhook_handle)

      j_start_u = 1
      j_stop_u = rows
      i_stop_u = row_length
      IF ( model_domain  ==  mt_global ) THEN
        IF (at_extremity(PSouth)) j_start_u = 2
        IF (at_extremity(PNorth)) j_stop_u = rows - 1
      END IF ! model_domain  ==  mt_global
      IF ( model_domain  ==  mt_LAM ) THEN
        IF (at_extremity(PEast)) i_stop_u = row_length - 1
      END IF ! model_domain  ==  mt_LAM

      itest = 0   ! Integer switch array > 0 if u,v > u/max_wind anywhere
                  !                       or if u,v < -u/max_wind anywhere

!-------------------------------------------------------------------
! 1.1  Test horizontal wind components
!-------------------------------------------------------------------

      ki = 1

      DO  k = 1, model_levels

        max_indices = MAXLOC(u_adv(1:i_stop_u,j_start_u:j_stop_u,k))
        i_max_u(k) = max_indices(1)
        j_max_u(k) = max_indices(2)
        max_u_pe(k) = u_adv(max_indices(1),max_indices(2),k)
        min_indices = MINLOC(u_adv(1:i_stop_u,j_start_u:j_stop_u,k))
        i_min_u(k) = min_indices(1)
        j_min_u(k) = min_indices(2)
        min_u_pe(k) = u_adv(min_indices(1),min_indices(2),k)

        IF ( trap_option < 2 ) THEN

          IF ( max_u_pe(k) > max_wind ) THEN
            ic = 0
            DO j = j_start_u, j_stop_u
              DO i = 1, i_stop_u
                IF ( u_adv(i,j,k) > max_wind ) THEN
                  ic = ic + 1
                  u_adv(i,j,k) = max_wind
                END IF ! u_adv(i,j,k) > max_wind
                IF ( u(i,j,k) > max_wind ) THEN
                  ic = ic + 1
                  u(i,j,k) = max_wind
                END IF ! u(i,j,k) > max_wind
              END DO
            END DO
            itest(ki) = itest(ki) + ic
          END IF ! max_u_pe(k) > max_wind

          IF ( min_u_pe(k) < -max_wind ) THEN
            ic = 0
            DO j = j_start_u, j_stop_u
              DO i = 1, i_stop_u
                IF ( u_adv(i,j,k) < -max_wind ) THEN
                  ic = ic + 1
                  u_adv(i,j,k) = -max_wind
                END IF ! u_adv(i,j,k) < -max_wind
                IF ( u(i,j,k) < -max_wind ) THEN
                  ic = ic + 1
                  u(i,j,k) = -max_wind
                END IF ! u(i,j,k) < -max_wind
              END DO
            END DO
            itest(ki+1) = itest(ki+1) + ic
          END IF ! min_u_pe(k) < -max_wind

        END IF ! trap_option < 2

      END DO  !  k = 1, model_levels

      ki = 3

      DO k = 1, model_levels

        max_indices = MAXLOC(v_adv(1:row_length,1:n_rows,k))
        i_max_v(k) = max_indices(1)
        j_max_v(k) = max_indices(2)
        max_v_pe(k) = v_adv(max_indices(1),max_indices(2),k)
        min_indices = MINLOC(v_adv(1:row_length,1:n_rows,k))
        i_min_v(k) = min_indices(1)
        j_min_v(k) = min_indices(2)
        min_v_pe(k) = v_adv(min_indices(1),min_indices(2),k)

        IF ( trap_option < 2 ) THEN

          IF ( max_v_pe(k) > max_wind ) THEN
            ic = 0
            DO j = 1, n_rows
              DO i = 1, row_length
                IF ( v_adv(i,j,k) > max_wind ) THEN
                  ic = ic + 1
                  v_adv(i,j,k) = max_wind
                END IF !  v_adv(i,j,k) > max_wind
                IF ( v(i,j,k) > max_wind ) THEN
                  ic = ic + 1
                  v(i,j,k) = max_wind
                END IF !  v(i,j,k) > max_wind
              END DO
            END DO
            itest(ki) = itest(ki) + ic
          END IF ! max_v_pe(k) > max_wind
          IF ( min_v_pe(k) < -max_wind ) THEN
            ic = 0
            DO j = 1, n_rows
              DO i = 1, row_length
                IF ( v_adv(i,j,k) < -max_wind ) THEN
                  ic = ic + 1
                  v_adv(i,j,k) = -max_wind
                END IF ! v_adv(i,j,k) < -max_wind
                IF ( v(i,j,k) < -max_wind ) THEN
                  ic = ic + 1
                  v(i,j,k) = -max_wind
                END IF ! v(i,j,k) < -max_wind
              END DO
            END DO
            itest(ki+1) = itest(ki+1) + ic
          END IF ! min_v_pe(k) < -max_wind

        END IF ! trap_option < 2

      END DO  !  k = 1, model_levels

! ----------------------------------------------------------------------
! Section 2. sum itest over all processors to see if any
!             u,v exceed thresholds
! ----------------------------------------------------------------------

      CALL gc_isum(4, n_proc, info, itest)

! ----------------------------------------------------------------------
! Section 3. For trap_option = 2 find max/mins over all processors
!            All required fields done at same time
! ----------------------------------------------------------------------

      IF ( trap_option > 1 ) THEN

! Copy local max/mins  into ?_max/min which will hold
!            global ?_max/min after CALL to gc_rmax

        ku = 0
        DO k = 1, model_levels

          ku = ku + 1
          kv = ku + model_levels
          max_real(ku) = max_u_pe(k)
          min_real(ku) = min_u_pe(k)
          max_real(kv) = max_v_pe(k)
          min_real(kv) = min_v_pe(k)

        END DO ! k = 1, model_levels

        CALL gc_rmax(kv, n_proc, info, max_real)
        CALL gc_rmin(kv, n_proc, info, min_real)

! ----------------------------------------------------------------------
! Section 4. Now locations of max/mins
! ----------------------------------------------------------------------

        recip_row_length = 1.0 / real(global_row_length)
        recip_rows = 1.0 / real(global_rows)

! ----------------------------------------------------------------------
! Section 4.1  Obtain max and mins from max_real, min_real
!              and fill isum, sumr arrays for summing over pe's
! ----------------------------------------------------------------------

        sumi = 0
!  Re-Initialise pointers in arrays for summing over pe's
        ki = -7
        ku = 0
        gi = l_datastart(1) - 1
        gj = l_datastart(2) - 1
        DO   k = 1, model_levels
          ku = ku + 1
          kv = ku + model_levels
          ki = ki + 8
          IF ( max_u_pe(k) >= max_real(ku)) THEN
            ! max is on this processor for this level
            i = gi + i_max_u(k)
            j = gj + j_max_u(k)
            sumi(ki) = nint( (i-1) * recip_row_length * 100.0 )
            sumi(ki+1) = nint( (j-1) * recip_rows * 100.0 )
          END IF  ! max_u_pe(k) >= max_real(ku)
          IF ( min_u_pe(k) <= min_real(ku)) THEN
            ! min is on this processor for this level
            i = gi + i_min_u(k)
            j = gj + j_min_u(k)
            sumi(ki+2) = nint( (i-1) * recip_row_length * 100.0 )
            sumi(ki+3) = nint( (j-1) * recip_rows * 100.0 )
          END IF  ! min_u_pe(k) <= min_real(ku)

          IF ( max_v_pe(k) >= max_real(kv)) THEN
            ! max is on this processor for this level
            i = gi + i_max_v(k)
            j = gj + j_max_v(k)
            sumi(ki+4) = nint( (i-1) * recip_row_length * 100.0 )
            sumi(ki+5) = nint( (j-1) * recip_rows * 100.0 )
          END IF  ! max_v_pe(k) >= max_real(kv)
          IF ( min_v_pe(k) <= min_real(kv)) THEN
            ! min is on this processor for this level
            i = gi + i_min_v(k)
            j = gj + j_min_v(k)
            sumi(ki+6) = nint( (i-1) * recip_row_length * 100.0 )
            sumi(ki+7) = nint( (j-1) * recip_rows * 100.0 )
          END IF  ! min_v_pe(k) <= min_real(kv)
        END DO  ! k = 1, model_levels

        ki = ki + 7

        CALL gc_isum(ki, n_proc, info, sumi)

      END IF ! trap_option > 1

! ----------------------------------------------------------------------
! Section 5  Summing over pe's to obtain sums and location of max/mins
! ----------------------------------------------------------------------

      itests = itest(1) + itest(2)

      IF ( itests > 0 ) THEN
! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
                         u, row_length, rows, model_levels,             &
                         off_x, off_y, fld_type_u, .TRUE.)
! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
                         u_adv, row_length, rows, model_levels,         &
                         halo_i, halo_j, fld_type_u, .TRUE.)
        IF ( trap_option == 1 ) THEN
          IF ( itest(1) > 0 )                                           &
            WRITE(6,'("Westerly", A, F7.1, A ,I6," points")')           &
                                  prulim, max_wind, mperspr, itest(1)
          IF ( itest(2) > 0 )                                           &
            WRITE(6,'("Easterly", A, F7.1, A ,I6," points")')           &
                                  prulim, max_wind, mperspr, itest(2)
        END IF ! trap_option == 1

      END IF ! itests > 0

      itests = itest(3) + itest(4)

      IF ( itests > 0 ) THEN
! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
                         v, row_length, n_rows, model_levels,           &
                         off_x, off_y, fld_type_v, .TRUE.)
! DEPENDS ON: swap_bounds
        CALL swap_bounds(                                               &
                         v_adv, row_length, n_rows, model_levels,       &
                         halo_i, halo_j, fld_type_v, .TRUE.)
        IF ( trap_option == 1 ) THEN
          IF ( itest(3) > 0 )                                           &
            WRITE(6,'("Southerly", A, F7.1, A ,I6," points")')          &
                                  prvlim, max_wind, mperspr, itest(3)
          IF ( itest(4) > 0 )                                           &
            WRITE(6,'("Northerly", A, F7.1, A ,I6," points")')          &
                                  prvlim, max_wind, mperspr, itest(4)
        END IF ! trap_option == 1

      END IF ! itests > 0

      IF ( trap_option == 2 ) THEN

        w_p = 0.0
        ku = 0
        ki = -7
        DO   k = 1, model_levels
          ku = ku + 1
          ki = ki + 8
          IF ( w_p < max_real(ku) ) THEN
            w_p = max_real(ku)
            level_p = k
            lambda_p = sumi(ki)
            phi_p = sumi(ki+1)
          END IF ! w_p < max_real(ku)
        END DO   !  k = 1, model_levels
        WRITE(6,'("Maximum westerly wind = ", F7.1, A, "level ", I4)')  &
                                                w_p, mperspr, level_p
        WRITE(6,'("Location in domain is ", I4, A8, I4, A8 )')          &
                                         lambda_p, East, phi_p, North

        w_p = 1000.0
        ki = -5
        ku = 0
        DO k = 1, model_levels
          ku = ku + 1
          ki = ki + 8
          IF ( w_p > min_real(ku) ) THEN
            w_p = min_real(ku)
            level_p = k
            lambda_p = sumi(ki)
            phi_p = sumi(ki+1)
          END IF ! w_p > min_real(ku)
        END DO   !  k = 1, model_levels
        w_p = -w_p
        WRITE(6,'("Maximum  easterly wind= ", F7.1, A, "level ", I4)')  &
                                            w_p, mperspr, level_p
        WRITE(6,'("Location in domain is ", I4, A8, I4, A8 )')          &
                                     lambda_p, East, phi_p, North

        w_p = 0.0
        kv = model_levels
        ki = -3
        DO k = 1, model_levels
          kv = kv + 1
          ki = ki + 8
          IF ( w_p < max_real(kv) ) THEN
            w_p = max_real(kv)
            level_p = k
            lambda_p = sumi(ki)
            phi_p = sumi(ki+1)
          END IF ! w_p < max_real(kv)
        END DO   !  k = 1, model_levels
        WRITE(6,'("Maximum southerly wind = ", F7.1, A, "level ", I4)') &
                                              w_p, mperspr, level_p
        WRITE(6,'("Location in domain is ", I4, A8, I4, A8 )')          &
                                       lambda_p, East, phi_p, North

        w_p = 1000.0
        kv = model_levels
        ki = - 1
        DO k = 1, model_levels
          kv = kv + 1
          ki = ki + 8
          IF ( w_p > min_real(kv) ) THEN
            w_p = min_real(kv)
            level_p = k
            lambda_p = sumi(ki)
            phi_p = sumi(ki+1)
          END IF ! w_p > min_real(ku)
        END DO   !  k = 1, model_levels
        w_p = -w_p
        WRITE(6,'("Maximum  Northerly wind= ", F7.1, A, "level ", I4)') &
                                              w_p, mperspr, level_p
        WRITE(6,'("Location in domain is ", I4, A8, I4, A8 )')          &
                                       lambda_p, East, phi_p, North

      END IF ! trap_option == 2

      IF (lhook) CALL dr_hook('TRAP_UV',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE trap_uv
