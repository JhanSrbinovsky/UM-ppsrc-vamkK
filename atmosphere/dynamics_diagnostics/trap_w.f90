! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine Trap_w

      SUBROUTINE trap_w(                                                &
                        w, w_adv, Cw_max, Cw_test_lev,                  &
                        rows, row_length, model_levels,                 &
                        r_theta_levels, r_rho_levels, timestep,         &
                        off_x, off_y, halo_i, halo_j,                   &
                        trap_option )
      USE proc_info_mod

! Purpose:
!          Diagnostic routine for trapping when w component
!             exceeds maximum Courant number.

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
                         ! Size of halo for wind component.
        halo_j,                                                         &
                         ! Size of halo for wind component.
        Cw_test_lev,                                                    &
                         ! lowest level to test Cw_max (default=0)
        trap_option      ! 0 = reset, no prints
                         ! 1 = reset + print message
                         ! 2 = NO reset but print details max/min

      REAL , INTENT(IN) ::                                              &
                           ! vertical co-ordinate arrays.
        r_theta_levels (1-halo_i:row_length+halo_i,                     &
                        1-halo_j:rows+halo_j, 0:model_levels),          &
        r_rho_levels (1-halo_i:row_length+halo_i,                       &
                      1-halo_j:rows+halo_j, model_levels)

      REAL , INTENT(INOUT) ::                                           &
       w (1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
          0:model_levels),                                              &
       w_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
              0:model_levels)

      REAL , INTENT(IN) ::  Cw_max
      REAL , INTENT(IN) ::  timestep

! Local Variables.

      INTEGER                                                           &
        i, j, k,                                                        &
                          ! Loop counters
        gi, gj,                                                         &
                          ! global pointers
        ki, kr, kv,                                                     &
                          ! pointers for arrays summed over pe's
        ic,                                                             &
        info,                                                           &
        itests,                                                         &
        j_start, j_stop,                                                &
                          ! Loop indices
        level_p,                                                        &
                          ! print level of max/min
        lambda_p,                                                       &
                          ! print % domain longitude of max/min
        phi_p             ! print % domain latitude of  max/min

! Local arrays

      REAL                                                              &
        w_test(row_length, rows),                                       &
        Cw(row_length, rows)

      REAL                                                              &
        max_Cw_pe(model_levels-1),                                      &
        min_Cw_pe(model_levels-1),                                      &
        max_real(model_levels-1),                                       &
        min_real(model_levels-1)

      INTEGER                                                           &
        i_max_w(model_levels-1),                                        &
        j_max_w(model_levels-1),                                        &
        i_min_w(model_levels-1),                                        &
        j_min_w(model_levels-1),                                        &
        sumi( 4 * model_levels - 4 ),                                   &
        ic_xy(row_length, rows)

      INTEGER  :: min_indices(2), max_indices(2)
      INTEGER  :: itest(3)

      REAL                                                              &
        delta_z,                                                        &
        recip_row_length,                                               &
        recip_rows,                                                     &
        w_p                  ! print max/min

      CHARACTER(LEN=8) North, East
      PARAMETER ( East=' % East ')
      PARAMETER ( North=' % North')
      CHARACTER(LEN=*) maxprint
      PARAMETER ( maxprint=' Courant number = ')
      CHARACTER(LEN=*) resetpr
      PARAMETER ( resetpr=' points by reseting w and w_adv')
      CHARACTER(LEN=*) resetcol
      PARAMETER ( resetcol=' columns by reseting w and w_adv')


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Section 0. Initialise values
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('TRAP_W',zhook_in,zhook_handle)

      j_start = 1
      j_stop = rows
      IF ( model_domain  ==  mt_global ) THEN
        IF (at_extremity(PSouth)) j_start = 2
        IF (at_extremity(PNorth)) j_stop = rows - 1
      END IF ! model_domain  ==  mt_global

      kr = 0

      itest = 0   ! Integer switch array > 0 if Cw > Cw_max anywhere
                  !                       or if Cw < -Cw_max anywhere

      Cw = 0.0
      ic_xy = 0

!-------------------------------------------------------------------
! 1.1  Test vertical wind components
!-------------------------------------------------------------------

      ki = 1
      DO  k = model_levels - 1, 1, -1

        DO j = j_start, j_stop
          DO i = 1, row_length
            delta_z = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
            Cw(i,j) = w_adv(i,j,k) * timestep / delta_z
            w_test(i,j) = Cw_max * delta_z / timestep
          END DO
        END DO

        max_indices = MAXLOC(Cw(1:row_length,1:rows))
        i_max_w(k) = max_indices(1)
        j_max_w(k) = max_indices(2)
        max_Cw_pe(k) = Cw(max_indices(1),max_indices(2))
        min_indices = MINLOC(Cw(1:row_length,1:rows))
        i_min_w(k) = min_indices(1)
        j_min_w(k) = min_indices(2)
        min_Cw_pe(k) = Cw(min_indices(1),min_indices(2))

        IF ( trap_option < 2 ) THEN

          IF ( max_Cw_pe(k) > Cw_max ) THEN
            ic = 0
            DO j = j_start, j_stop
              DO i = 1, row_length
                IF ( w_adv(i,j,k) > w_test(i,j) ) THEN
                  IF ((k>Cw_test_lev).OR.(ic_xy(i,j)>0)) THEN
                    ic_xy(i,j) = 1
                    ic = ic + 1
                    w_adv(i,j,k) = w_test(i,j)
                  END IF
                END IF ! w_adv(i,j,k) > w_test(i,j)
                IF ( w(i,j,k) > w_test(i,j) ) THEN
                  IF ((k>Cw_test_lev).OR.(ic_xy(i,j)>0)) THEN
                    ic_xy(i,j) = 1
                    ic = ic + 1
                    w(i,j,k) = w_test(i,j)
                  END IF
                END IF ! w(i,j,k) > w_test(i,j)
              END DO
            END DO
            itest(ki) = itest(ki) + ic
          END IF ! max_Cw_pe(k) > Cw_max
          IF ( min_Cw_pe(k) < -Cw_max ) THEN
            ic = 0
            DO j = j_start, j_stop
              DO i = 1, row_length
                w_test(i,j) = -1.0 * w_test(i,j)
                IF ( w_adv(i,j,k) < w_test(i,j) ) THEN
                  IF ((k>Cw_test_lev).OR.(ic_xy(i,j)>0)) THEN
                    ic_xy(i,j) = 1
                    ic = ic + 1
                    w_adv(i,j,k) = w_test(i,j)
                  END IF
                END IF ! w_adv(i,j,k) < w_test(i,j)
                IF ( w(i,j,k) < w_test(i,j) ) THEN
                  IF ((k>Cw_test_lev).OR.(ic_xy(i,j)>0)) THEN
                    ic_xy(i,j) = 1
                    ic = ic + 1
                    w(i,j,k) = w_test(i,j)
                  END IF
                END IF ! w(i,j,k) < w_test(i,j)
              END DO
            END DO
            itest(ki+1) = itest(ki+1) + ic
          END IF ! min_Cw_pe(k) > -Cw_max

        END IF ! trap_option < 2

      END DO  !  k = model_levels - 1, 1, -1

! Sum up number of columns acted upon and store in itest(3)
      DO j = j_start, j_stop
        DO i = 1, row_length
          itest(ki+2) = itest(ki+2) + ic_xy(i,j)
        END DO
      END DO
! ----------------------------------------------------------------------
! Section 2. sum itest over all processors to see if any
!            of u,v,w exceed thresholds
! ----------------------------------------------------------------------

      CALL gc_isum(3, n_proc, info, itest)

! ----------------------------------------------------------------------
! Section 3. For trap_option = 2 find max/mins over all processors
!            All required fields done at same time
! ----------------------------------------------------------------------

      IF ( trap_option > 1 ) THEN

! Copy local max/mins  into ?_max/min which will hold
!            global ?_max/min after CALL to gc_rmax

        kr = 0

        DO k = 1, model_levels - 1
          kr = kr + 1
          max_real(kr) = max_Cw_pe(k)
          min_real(kr) = min_Cw_pe(k)
        END DO ! k = 1, model_levels - 1

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
          gi = l_datastart(1) - 1
          gj = l_datastart(2) - 1
          DO   k = 1, model_levels - 1
            kr = kr + 1
            ki = ki + 4
            IF ( max_Cw_pe(k) >= max_real(kr)) THEN
              ! max is on this processor for this level
              i = gi + i_max_w(k)
              j = gj + j_max_w(k)
              sumi(ki) = nint( (i-1) * recip_row_length * 100.0 )
              sumi(ki+1) = nint( (j-1) * recip_rows * 100.0 )
            END IF  ! max_Cw_pe(k) >= max_real(kr)
            IF ( min_Cw_pe(k) <= min_real(kr)) THEN
              ! min is on this processor for this level
              i = gi + i_min_w(k)
              j = gj + j_min_w(k)
              sumi(ki+2) = nint( (i-1) * recip_row_length * 100.0 )
              sumi(ki+3) = nint( (j-1) * recip_rows * 100.0 )
            END IF  ! min_Cw_pe(k) <= min_real(kr)
          END DO  ! k = 1, model_levels - 1

          ki = ki + 3

          CALL gc_isum(ki, n_proc, info, sumi)

        END IF ! trap_option > 1

! ----------------------------------------------------------------------
! Section 5  Summing over pe's to obtain sums and location of max/mins
! ----------------------------------------------------------------------

      itests = itest(1) + itest(2)

      IF ( itests > 0 ) THEN
! DEPENDS ON: Swap_bounds
        CALL Swap_bounds(                                               &
                         w(1-off_x, 1-off_y, 1),                        &
                         row_length, rows, model_levels - 1,            &
                         off_x, off_y, fld_type_p, .FALSE.)
! DEPENDS ON: Swap_bounds
        CALL Swap_bounds(                                               &
                         w_adv(1-halo_i, 1-halo_j, 1),                  &
                         row_length, rows, model_levels - 1,            &
                         halo_i, halo_j, fld_type_p, .FALSE.)

        IF ( trap_option == 1 ) THEN
          IF ( itest(1) > 0 )                                           &
            WRITE(6,'("Maximum upward", A,"limited at", I6, A)')        &
                                         maxprint, itest(1), resetpr
          IF ( itest(2) > 0 )                                           &
            WRITE(6,'("Maximum downward", A," limited at", I6, A)')     &
                                         maxprint, itest(2), resetpr
          IF ( itest(3) > 0 )                                           &
            WRITE(6,'("Vertical", A," limited in", I6, A)')             &
                                maxprint, itest(3), resetcol
        END IF ! trap_option == 1

      END IF ! itests > 0

      IF ( trap_option == 2 ) THEN

        w_p = 0.0
        kr = 0
        ki = -3
        DO   k = 1, model_levels - 1
          kr = kr + 1
          ki = ki + 4
          IF ( w_p < max_real(kr) ) THEN
            w_p = max_real(kr)
            level_p = k
            lambda_p = sumi(ki)
            phi_p = sumi(ki+1)
          END IF ! w_p < max_real(kr)
        END DO   !  k = 1, model_levels - 1
        WRITE(6,'("Maximum upward", A, F7.1," at level ", I4)')         &
                                             maxprint, w_p, level_p
        WRITE(6,'("Location in domain is ", I4, A8, I4, A8 )')          &
                   lambda_p, East, phi_p, North

        w_p = 1000.0
        kr = 0
        ki = -1
        DO   k = 1, model_levels - 1
          kr = kr + 1
          ki = ki + 4
          IF ( w_p > min_real(kr) ) THEN
            w_p = min_real(kr)
            level_p = k
            lambda_p = sumi(ki)
            phi_p = sumi(ki+1)
          END IF ! w_p > min_real(kr)
        END DO   !  k = 1, model_levels - 1
        WRITE(6,'("Maximum downward", A, F7.1," at level ", I4)')       &
                                             maxprint, w_p, level_p
        WRITE(6,'("Location in domain is ", I4, A8, I4, A8 )')          &
                   lambda_p, East, phi_p, North

      END IF ! trap_option == 2

      IF (lhook) CALL dr_hook('TRAP_W',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE trap_w
