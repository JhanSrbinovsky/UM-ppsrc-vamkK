! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine GCR_error_print

      subroutine GCR_error_print(                                       &
     &                            Error, gc_proc_row_group,             &
     &                            global_row_length, global_rows,       &
     &                            g_rows, g_row_length,                 &
     &                            n_proc, n_procx, n_procy, me,         &
     &                            row_length, rows, model_levels,       &
     &                            model_domain, at_extremity,           &
     &                            init_error_mean,switch)

! Purpose:
!          Diagnostic print of residual error norms
!
! Method:
!          Is described in ;
!
!          Documentation yet to be written
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver



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
     &, switch                                                          &
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, gc_proc_row_group                                               &
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, me                                                              &
     &, model_domain

      INTEGER                                                           &
     &  g_rows(0:n_proc-1)                                              &
     &, g_row_length(0:n_proc-1)                                        

      REAL                                                              &
     &  Error (row_length, rows, model_levels)

      LOGICAL                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid


! Local variables
      REAL                                                              &
     &  init_Error_mean (global_rows)                                   &
     &, error_sq(row_length, rows)                                      &
     &, mean(global_rows)                                               &
     &, rbuf(global_rows)

      REAL                                                              &
     &  non_zero

      INTEGER d, i, j, k, istat, info, gj                               &
     &,  i_start, i_end                                                 &
     &,  j_start, j_end

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!    External Routines:
      External                                                          &
     &  gcg_rvecsumr, gc_rsend, gc_rrecv

! ----------------------------------------------------------------------
! Section 1.   Calculate error norms
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GCR_ERROR_PRINT',zhook_in,zhook_handle)
      non_zero = 1.e-10
      i_start = 1
      i_end = row_length
      j_start = 1
      j_end = rows
      IF (model_domain  ==  mt_lam) THEN
        IF(at_extremity(PSouth)) j_start = 2
        IF(at_extremity(PNorth)) j_end = rows-1
        IF(at_extremity(PEast)) i_end = row_length-1
        IF(at_extremity(PWest)) i_start = 2
      ELSE IF (model_domain  ==  mt_cyclic_lam) THEN
        IF(at_extremity(PSouth)) j_start = 2
        IF(at_extremity(PNorth)) j_end = rows-1
      END IF

! If switch = 0 then set up initial error mean
! else calculate row error means.

      DO j = 1, rows
        DO i = 1, row_length
          error_sq(i,j) = 0.0
        END DO
      END DO
      DO k = 1, model_levels
        DO j = j_start, j_end
          DO i = i_start, i_end
            error_sq(i,j) = error_sq(i,j) +                             &
     &                      error(i,j,k) * error(i,j,k)
          END DO
        END DO
      END DO

      CALL gcg_rvecsumr(row_length,row_length,1,                        &
     &   rows,error_sq,gc_proc_row_group,istat,mean)

      IF (model_domain  /=  mt_lam) THEN
        DO j = 1, rows
          mean(j) = SQRT (mean(j) ) / (model_levels*global_row_length)
        END DO
      ELSE
        DO j = 1, rows
          mean(j) = SQRT (mean(j) ) /                                   &
     &              (model_levels*(global_row_length-2))
        END DO
      END IF

! Gather all data to processor 0 from its column and form
! global rows fields

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

        IF (switch  ==  0) THEN
          DO j = 1, global_rows
            init_Error_mean(j) = max(mean(j),non_zero)
          END DO
        ELSE
          WRITE(6,'(A)') ' End of convergence tolerances achieved '
          DO j = 1, global_rows
            mean(j) = mean(j) / init_Error_mean(j)
            WRITE(*,'(1x,I3,1x,E10.2)') j, mean(j)
          END DO
        END IF

      END IF

      IF (lhook) CALL dr_hook('GCR_ERROR_PRINT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GCR_error_print
