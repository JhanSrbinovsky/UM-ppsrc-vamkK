! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine trap_theta

      SUBROUTINE trap_theta(                                            &
                            theta, theta_star, max_inc,                 &
                            rows, row_length, model_levels,             &
                            off_x, off_y, halo_i, halo_j,               &
                            trap_option )
      USE proc_info_mod

! Purpose:
!          Diagnostic routine for trapping when diabatic
!          heating/cooling is too strong

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Diagnostics

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE UM_ParParams
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      INTEGER , INTENT(IN) ::                                           &
        row_length,                                                     &
                         ! number of point on a row.
        rows,                                                           &
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
        theta (1-off_x:row_length+off_x,                                &
               1-off_y:rows+off_y, model_levels),                       &
        theta_star (1-off_x:row_length+off_x,                           &
                    1-off_y:rows+off_y, model_levels)

      REAL , INTENT(IN) :: max_inc

! Local Variables.

      INTEGER                                                           &
        i, j, k,                                                        &
                            ! Loop counters
        gi, gj,                                                         &
                            ! global pointers
        ki, kr,                                                         &
                            ! pointers for arrays summed over pe's
        ic,                                                             &
        info,                                                           &
        itests,                                                         &
        j_start, j_stop,                                                &
                            ! Loop indices
        level_p,                                                        &
                            ! print height of max/min
        lambda_p,                                                       &
                            ! print % domain longitude of max/min
        phi_p               ! print % domain latitude of  max/min

! Local arrays

      REAL  inc(row_length, rows)

      REAL                                                              &
        max_inc_pe(model_levels),                                       &
        min_inc_pe(model_levels),                                       &
        max_real(model_levels),                                         &
        min_real(model_levels)

      INTEGER                                                           &
        i_max(model_levels),                                            &
        j_max(model_levels),                                            &
        i_min(model_levels),                                            &
        j_min(model_levels),                                            &
        sumi( 4 * model_levels )

      INTEGER :: min_indices(2), max_indices(2), itest(2)

      REAL                                                              &
        recip_row_length,                                               &
        recip_rows,                                                     &
        w_p                     ! print max/min

      CHARACTER(LEN=8) North, East
      PARAMETER ( North=' % North')
      PARAMETER ( East=' % East ')
      CHARACTER(LEN=*) maxprint
      PARAMETER ( maxprint='Maximum diabatic ')
      CHARACTER(LEN=*) resetpr
      PARAMETER ( resetpr=' points by reseting theta_star ')


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Section 0. Initialise values
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('TRAP_THETA',zhook_in,zhook_handle)

      j_start = 1
      j_stop = rows
      IF ( model_domain  ==  mt_global ) THEN
        IF (at_extremity(PSouth)) j_start = 2
        IF (at_extremity(PNorth)) j_stop = rows - 1
      END IF ! model_domain  ==  mt_global

      itest = 0   ! Integer switch array > 0 if inc > max_inc anywhere
                  !                       or if inc < -max_inc anywhere
      inc = 0.0

!-------------------------------------------------------------------
! 1.1  Test increments
!-------------------------------------------------------------------

      ki = 1
      DO  k = 1, model_levels

        DO j = j_start, j_stop
          DO i = 1, row_length
            inc(i,j) = theta_star(i,j,k) - theta(i,j,k)
          END DO
        END DO

        max_indices = MAXLOC(inc(1:row_length,1:rows))
        i_max(k) = max_indices(1)
        j_max(k) = max_indices(2)
        max_inc_pe(k) = inc(max_indices(1),max_indices(2))
        min_indices = MINLOC(inc(1:row_length,1:rows))
        i_min(k) = min_indices(1)
        j_min(k) = min_indices(2)
        min_inc_pe(k) = inc(min_indices(1),min_indices(2))

        IF ( trap_option < 2 ) THEN

          IF ( max_inc_pe(k) > max_inc ) THEN
            ic = 0
            DO j = j_start, j_stop
              DO i = 1, row_length
                IF ( inc(i,j) > max_inc ) THEN
                  ic = ic + 1
                  theta_star(i,j,k) = theta(i,j,k) + max_inc
                END IF ! inc(i,j) > max_inc
              END DO
            END DO
            itest(ki) = itest(ki) + ic
          END IF ! max_inc_pe(k) > max_inc
          IF ( min_inc_pe(k) < -max_inc ) THEN
            ic = 0
            DO j = j_start, j_stop
              DO i = 1, row_length
                IF ( inc(i,j) < -max_inc ) THEN
                  ic = ic + 1
                  theta_star(i,j,k) = theta(i,j,k) - max_inc
                END IF ! inc(i,j) < -max_inc
              END DO
            END DO
            itest(ki+1) = itest(ki+1) + ic
          END IF ! min_inc_pe(k) > -max_inc

        END IF ! trap_option < 2

      END DO  !  k = 1, model_levels

! ----------------------------------------------------------------------
! Section 2. sum itest over all processors to see if
!            increments exceed thresholds
! ----------------------------------------------------------------------

      CALL GC_ISUM(2, n_proc, info, itest)

! ----------------------------------------------------------------------
! Section 3. For trap_option = 2 find max/mins over all processors
!            All required fields done at same time
! ----------------------------------------------------------------------

      IF ( trap_option > 1 ) THEN

! Copy local max/mins  into ?_max/min which will hold
!            global ?_max/min after CALL to gc_rmax

        kr = 0

        DO k = 1, model_levels
          kr = kr + 1
          max_real(kr) = max_inc_pe(k)
          min_real(kr) = min_inc_pe(k)
        END DO ! k = 1, model_levels

        CALL gc_rmax(kr, n_proc, info, max_real)
        CALL gc_rmin(kr, n_proc, info, min_real)

! ----------------------------------------------------------------------
! Section 4. For trap_option = 2 find locations of max/mins
! ----------------------------------------------------------------------

          sumi = 0
          recip_row_length = 1.0 / real(global_row_length)
          recip_rows = 1.0 / real(global_rows)

! ----------------------------------------------------------------------
! Section 4.1  Obtain max and mins from max_real, min_real
!              and fill sumi arrays for summing over pe's
! ----------------------------------------------------------------------

          kr = 0
          ki = -3
          DO   k = 1, model_levels
            kr = kr + 1
            ki = ki + 4
            IF ( max_inc_pe(k) >= max_real(kr)) THEN
              ! max is on this processor for this level
              gi = l_datastart(1) + i_max(k) - 1
              gj = l_datastart(2) + j_max(k) - 1
              sumi(ki) = nint( (gi-1) * recip_row_length * 100.0 )
              sumi(ki+1) = nint( (gj-1) * recip_rows * 100.0 )
            END IF  ! max_inc_pe(k) >= max_real(kr)
            IF ( min_inc_pe(k) <= min_real(kr)) THEN
              ! min is on this processor for this level
              gi = l_datastart(1) + i_min(k) - 1
              gj = l_datastart(2) + j_min(k) - 1
              sumi(ki+2) = nint( (gi-1) * recip_row_length * 100.0 )
              sumi(ki+3) = nint( (gj-1) * recip_rows * 100.0 )
            END IF  ! min_inc_pe(k) <= min_real(kr)
          END DO  ! k = 1, model_levels

          CALL gc_isum(ki, n_proc, info, sumi)

        END IF ! trap_option > 1

! ----------------------------------------------------------------------
! Section 5  Summing over pe's to obtain sums and location of max/mins
! ----------------------------------------------------------------------

      itests = itest(1) + itest(2)

      IF ( itests > 0 ) THEN
! DEPENDS ON: Swap_bounds
        CALL Swap_bounds(                                               &
                         theta_star, row_length, rows, model_levels,    &
                         off_x, off_y, fld_type_p, .FALSE.)

        IF ( trap_option == 1 ) THEN
          IF ( itest(1) > 0 )                                           &
            WRITE(6,'(A, "heating limited at", I6, A )')                &
                                         maxprint, itest(1), resetpr
          IF ( itest(2) > 0 )                                           &
            WRITE(6,'(A, "cooling limited at", I6, A )')                &
                                         maxprint, itest(2), resetpr
        END IF ! trap_option == 1

      END IF ! itests > 0

      IF ( trap_option == 2 ) THEN

        w_p = 0.0
        kr = 0
        ki = -3
        Do   k = 1, model_levels
          kr = kr + 1
          ki = ki + 4
          IF ( w_p < max_real(kr) ) THEN
            w_p = max_real(kr)
            level_p = k
            lambda_p = sumi(ki)
            phi_p = sumi(ki+1)
          END IF ! w_p < max_real(kr)
        END DO   !  k = 1, model_levels
        WRITE(6,'(A, "heating = ", F7.1, " at level ", I4)')            &
                                            maxprint, w_p, level_p
        WRITE(6,'("Location in domain is ", I4, A8, I4, A8 )')          &
                                      lambda_p, East, phi_p, North

        w_p = 1000.0
        kr = 0
        ki = -1
        Do   k = 1, model_levels
          kr = kr + 1
          ki = ki + 4
          IF ( w_p > min_real(kr) ) THEN
            w_p = min_real(kr)
            level_p = k
            lambda_p = sumi(ki)
            phi_p = sumi(ki+1)
          END IF ! w_p > min_real(kr)
        END DO   !  k = 1, model_levels
        WRITE(6,'(A, "cooling = ", F7.1, " at level ", I4)')            &
                                               maxprint, w_p, level_p
        WRITE(6,'("Location in domain is ", I4, A8, I4, A8 )')          &
                                         lambda_p, East, phi_p, North

      END IF ! trap_option == 2

      IF (lhook) CALL dr_hook('TRAP_THETA',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Trap_theta
