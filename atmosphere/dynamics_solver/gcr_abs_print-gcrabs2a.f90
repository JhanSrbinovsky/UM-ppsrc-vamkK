! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine GCR_abs_print

      subroutine GCR_abs_print(                                         &
     &                           Error, off_x, off_y,                   &
     &                            gc_proc_row_group,                    &
     &                            global_row_length, global_rows,       &
     &                            g_rows, g_row_length,                 &
     &                            n_proc, n_procx, n_procy, me,         &
     &                            row_length, rows,                     &
     &                            model_levels, model_domain,           &
     &                            at_extremity)

! Purpose:
!          Diagnostic print of absolute error norms
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver

      USE global_2d_sums_mod, ONLY : global_2d_sums
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      USE UM_ParVars, ONLY : g_datastart
      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.
      INTEGER                                                           &
     &  row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, gc_proc_row_group                                               &
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, me                                                              &
     &, off_x                                                           &
     &, off_y                                                           &
     &, model_domain

      LOGICAL                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


      INTEGER                                                           &
     &  g_rows(0:n_proc-1)                                              &
     &, g_row_length(0:n_proc-1)                                        

      REAL                                                              &
     &  Error (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)

! Local variables
      REAL                                                              &
     &  error_sum(row_length, rows)                                     &
     &, mean(global_rows)                                               &
     &, rbuf(global_rows)                                               &
     &, max_val(model_levels)

      INTEGER d, i, j, k, istat, info, gj                               &
     &,       i_start, i_end, j_start, j_end

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!    External Routines:

! ----------------------------------------------------------------------
! Section 1.   Calculate error norm over rows and levels.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GCR_ABS_PRINT',zhook_in,zhook_handle)

      i_start = 1
      i_end = row_length
      j_start = 1
      j_end = rows
      IF (model_domain  ==  mt_lam) THEN
        IF(at_extremity(PSouth)) j_start = 2
        IF(at_extremity(PNorth)) j_end = rows-1
        IF(at_extremity(PEast)) i_end = row_length-1
        IF(at_extremity(PWest)) i_start = 2
      END IF
      IF (model_domain  ==  mt_cyclic_LAM) THEN
        IF(at_extremity(PSouth)) j_start = 2
        IF(at_extremity(PNorth)) j_end = rows-1
      END IF

      DO k = 1, model_levels
        max_val(k) = 0.0
      END DO
      k = 1
      DO j = 1, rows
        DO i = 1, row_length
          error_sum(i,j) = 0.0
        END DO
      END DO
      DO k = 1, model_levels
        DO j = j_start, j_end
          DO i = i_start, i_end
            error_sum(i,j) = error_sum(i,j) + error(i,j,k)
            IF (error(i,j,k)  >   max_val(k) )                          &
     &        max_val(k) = error(i,j,k)
          END DO
        END DO
      END DO

      CALL global_2d_sums(error_sum, row_length, 1, 0, 0, rows,         &
                          mean, gc_proc_row_group)

      IF (model_domain  /=  mt_lam) THEN
        DO j = 1, rows
          mean(j) = mean(j) / (model_levels*global_row_length)
        END DO
      ELSE
        DO j = 1, rows
          mean(j) = mean(j) / (model_levels*(global_row_length-2))
        END DO
      END IF

! Gather all data to processor 0 from its column and form
! global rows fields

      CALL gc_imax(model_levels, n_proc, info, max_val)

! step 1 send data to processor zero
      IF (me  /=  0) THEN
        DO j = 1, rows
          rbuf(j) = mean(j)
        END DO

        CALL gc_rsend(100*me, g_rows(me),                               &
     &                  0, info, rbuf, rbuf)
      END IF

! step 2 processor zero receives data and puts in correct location

      IF (me  ==  0) THEN
        DO d = 1, n_procx*n_procy-1
          CALL gc_rrecv(100*d, g_rows(d),                               &
     &                  d, info, rbuf, rbuf)
          DO j = 1, g_rows(d)
            gj = g_datastart(2,d) + j - 1
            mean(gj) = rbuf(j)
          END DO
        END DO

        WRITE(6,'(A)') ' mean approximate absolute density change left '
        DO j = 1, global_rows
          WRITE(*,'(1x,I3,1x,E10.2)') j, mean(j)
        END DO
        WRITE(6,'(A)') ' max value left per level '
        DO k = 1, model_levels
            WRITE(*,'(1x,I3,1x,E10.2)') k, max_val(k)
        END DO

      END IF

      IF (lhook) CALL dr_hook('GCR_ABS_PRINT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_abs_print
