! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine trap_array

      SUBROUTINE trap_array(                                            &
                            array, max_val,                             &
                            rows, row_length, levels,                   &
                            halo_x, halo_y, halo_i, halo_j,             &
                            fld_type, trap_option )
      USE proc_info_mod

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

! Purpose:
!          Diagnostic routine for trapping array values

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Diagnostics

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      USE UM_ParParams
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      INTEGER , INTENT(IN) ::                                           &
       row_length,                                                      &
                       ! number of point on a row.
       rows,                                                            &
                       ! number of rows.
       levels,                                                          &
                       ! number of model levels.
       halo_x,                                                          &
                       ! Size of small halo in i
       halo_y,                                                          &
                       ! Size of small halo in j.
       halo_i,                                                          &
                       ! Size of halo for wind components.
       halo_j,                                                          &
                       ! Size of halo for wind component.
       fld_type,                                                        &
                       ! field type for swap_bounds
       trap_option     ! 0 = reset, no prints
                       ! 1 = reset + print message
                       ! 2 = NO reset but print details max/min

      REAL , INTENT(INOUT) ::  array(1-halo_x:row_length+halo_x,        &
                                     1-halo_y:rows+halo_y, levels)

      REAL , INTENT(IN) :: max_val

! Local Variables.

      INTEGER                                                           &
       i, j, k,                                                         &
                       ! Loop counters
       gi, gj,                                                          &
                       ! global pointers
       ki, kr,                                                          &
                       ! pointers for arrays summed over pe's
       ic,                                                              &
       info,                                                            &
       itests,                                                          &
       j_start, j_stop,                                                 &
                       ! Loop indices
       level_p,                                                         &
                       ! print height of max/min
       lambda_p,                                                        &
                       ! print % domain longitude of max/min
       phi_p           ! print % domain latitude of  max/min

! Local arrays

      REAL                                                              &
        max_val_pe(levels),                                             &
        min_val_pe(levels),                                             &
        max_real(levels),                                               &
        min_real(levels)

      INTEGER                                                           &
        i_max(levels),                                                  &
        j_max(levels),                                                  &
        i_min(levels),                                                  &
        j_min(levels),                                                  &
        sumi( 4 * levels )

      INTEGER :: min_indices(2), max_indices(2), itest(2)

      REAL                                                              &
        recip_row_length,                                               &
        recip_rows,                                                     &
        w_p                      ! print max/min

      CHARACTER(LEN=8) North
      PARAMETER ( North=' % North')
      CHARACTER(LEN=8) East
      PARAMETER ( East=' % East ')
      CHARACTER(LEN=*) maxprint
      PARAMETER ( maxprint='Maximum array value')
      CHARACTER(LEN=*) resetpr
      PARAMETER ( resetpr=' points by reseting array ')

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 0. Initialise values
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('TRAP_ARRAY',zhook_in,zhook_handle)

      j_start = 1
      j_stop = rows
      IF ( model_domain  ==  mt_global ) THEN
        IF (at_extremity(PSouth)) j_start = 2
        IF (at_extremity(PNorth)) j_stop = rows - 1
      END IF ! model_domain  ==  mt_global

      itest = 0       ! Switch  > 0 if array > max_val anywhere
                      !          or if array < -max_val anywhere

!-------------------------------------------------------------------
! 1.1  Test valrements
!-------------------------------------------------------------------

      ki = 1
      DO  k = 1, levels

        max_indices = MAXLOC(array(1:row_length,1:rows,k))
        i_max(k) = max_indices(1)
        j_max(k) = max_indices(2)
        max_val_pe(k) = array(max_indices(1),max_indices(2),k)
        min_indices = MINLOC(array(1:row_length,1:rows,k))
        i_min(k) = min_indices(1)
        j_min(k) = min_indices(2)
        min_val_pe(k) = array(min_indices(1),min_indices(2),k)

        IF ( trap_option < 2 ) THEN

          IF ( max_val_pe(k) > max_val ) THEN
            ic = 0
            DO j = j_start, j_stop
              DO i = 1, row_length
                IF ( array(i,j,k) > max_val ) THEN
                  ic = ic + 1
                  array(i,j,k) = max_val
                END IF ! array(i,j,k) > max_val
              END DO
            END DO
            itest(ki) = itest(ki) + ic
          END IF ! max_val_pe(k) > max_val
          IF ( min_val_pe(k) < -max_val ) THEN
            ic = 0
            DO j = j_start, j_stop
              DO i = 1, row_length
                IF ( array(i,j,k) < -max_val ) THEN
                  ic = ic + 1
                  array(i,j,k) =  - max_val
                END IF ! array(i,j,k) < -max_val
              END DO
            END DO
            itest(ki+1) = itest(ki+1) + ic
          END IF ! min_val_pe(k) > -max_val

        END IF ! trap_option < 2

      END DO  !  k = 1, levels

! ----------------------------------------------------------------------
! Section 2. sum itest over all processors to see if
!            valrements exceed thresholds
! ----------------------------------------------------------------------

      CALL gc_isum(2, n_proc, info, itest)

! ----------------------------------------------------------------------
! Section 3. For trap_option = 2 find max/mins over all processors
!            All required fields done at same time
! ----------------------------------------------------------------------

      IF ( trap_option > 1 ) THEN

! Copy local max/mins  into ?_max/min which will hold
!            global ?_max/min after CALL to gc_rmax

        kr = 0

        DO k = 1, levels
          kr = kr + 1
          max_real(kr) = max_val_pe(k)
          min_real(kr) = min_val_pe(k)
        END DO ! k = 1, levels

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
          DO   k = 1, levels
            kr = kr + 1
            ki = ki + 4
            IF ( max_val_pe(k) >= max_real(kr)) THEN
              ! max is on this processor for this level
              gi = l_datastart(1) + i_max(k) - 1
              gj = l_datastart(2) + j_max(k) - 1
              sumi(ki) = nint( (gi-1) * recip_row_length * 100.0 )
              sumi(ki+1) = nint( (gj-1) * recip_rows * 100.0 )
            END IF  ! max_val_pe(k) >= max_real(kr)
            IF ( min_val_pe(k) <= min_real(kr)) THEN
              ! min is on this processor for this level
              gi = l_datastart(1) + i_min(k) - 1
              gj = l_datastart(2) + j_min(k) - 1
              sumi(ki+2) = nint( (gi-1) * recip_row_length * 100.0 )
              sumi(ki+3) = nint( (gj-1) * recip_rows * 100.0 )
            END IF  ! min_val_pe(k) <= min_real(kr)
          END DO  ! k = 1, levels

          CALL gc_isum(ki, n_proc, info, sumi)

        END IF ! trap_option > 1

! ----------------------------------------------------------------------
! Section 5  Summing over pe's to obtain sums and location of max/mins
! ----------------------------------------------------------------------

      itests = itest(1) + itest(2)

      IF ( itests > 0 ) THEN
! DEPENDS ON: Swap_bounds
        CALL Swap_bounds(                                               &
                         array, row_length, rows, levels,               &
                         halo_x, halo_y, fld_type, .FALSE.)

        IF ( trap_option == 1 ) THEN
          IF ( itest(1) > 0 )                                           &
            WRITE(6,'(A, " limited at", I6, A)')                        &
                                        maxprint, itest(1), resetpr
          IF ( itest(2) > 0 )                                           &
            WRITE(6,'(A, " limited at", I6, A)')                        &
                                        maxprint, itest(2), resetpr
        END IF ! trap_option == 1

      END IF ! itests > 0

      IF ( trap_option == 2 ) THEN

        w_p = 0.0
        kr = 0
        ki = -3
        DO   k = 1, levels
          kr = kr + 1
          ki = ki + 4
          IF ( w_p < max_real(kr) ) THEN
            w_p = max_real(kr)
            level_p = k
            lambda_p = sumi(ki)
            phi_p = sumi(ki+1)
          END IF ! w_p < max_real(kr)
        END DO   !  k = 1, levels
        WRITE(6,'(A," = ", F7.1," at level ", I4)')                     &
                     maxprint, w_p, level_p
        WRITE(6,'("Location in domain is ", I4, A8, I4, A8 )')          &
                     lambda_p, East, phi_p, North

        w_p = 1000.0
        kr = 0
        ki = -1
        DO   k = 1, levels
          kr = kr + 1
          ki = ki + 4
          IF ( w_p > min_real(kr) ) THEN
            w_p = min_real(kr)
            level_p = k
            lambda_p = sumi(ki)
            phi_p = sumi(ki+1)
          END IF ! w_p > min_real(kr)
        END DO   !  k = 1, levels
        WRITE(6,'(A," = ", F7.1," at level ", I4)')                     &
                     maxprint, w_p, level_p
        WRITE(6,'("Location in domain is ", I4, A8, I4, A8 )')          &
                     lambda_p, East, phi_p, North

      END IF ! trap_option == 2

      IF (lhook) CALL dr_hook('TRAP_ARRAY',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Trap_array
