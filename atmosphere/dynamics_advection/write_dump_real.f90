! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine write_dump_real

      Subroutine write_dump_real                                        &
     &                  (row_length, rows, n_levels,                    &
     &                   global_row_length, global_rows,                &
     &                   l_datastart, off_x, off_y,                     &
     &                   me, n_proc, n_procx, n_procy,                  &
     &                   g_rows, g_row_length, dun,                     &
     &                   field)

! Purpose:
!          Gathers a global real field and writes it out to a file.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars, ONLY : g_datastart
      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Integer                                                           &
     &  row_length                                                      &
                        ! number of points on a row
     &, rows                                                            &
                        ! number of rows of data
     &, n_levels                                                        &
     &, l_datastart(3)                                                  &
     &, off_x                                                           &
     &, off_y                                                           &
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, me                                                              &
     &, n_proc                                                          &
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, g_rows (0:n_proc-1)                                             &
     &, g_row_length (0:n_proc-1)                                       &
     &, dun

! Arguments with Intent OUT. ie: Output

      Real                                                              &
     &  field(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &        n_levels)

! local variables
      Real                                                              &
     &  global_real (global_row_length, global_rows)                    &
     &, rbuf(global_row_length, global_rows)

      Integer                                                           &
     & i, j, k, gi, gj, d, info

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 1. Collect all data onto processor zero then write out level
!            by level.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('WRITE_DUMP_REAL',zhook_in,zhook_handle)
      Do k = 1, n_levels

! step 1 send data to processor zero
        If (me  /=  0) then
          Do j = 1, rows
            Do i = 1, row_length
              rbuf(i,j) = field(i,j,k)
            End Do
          End Do

          call gc_rsend(100*me, global_rows*global_row_length,          &
     &                  0, info, rbuf, rbuf)
        Else
          Do j = 1, rows
            Do i = 1, row_length
              global_real(i,j) = field(i,j,k)
            End Do
          End Do
        End If

! step 2 processor zero receives data and puts in correct location

        If (me  ==  0) then
          Do d = 1, n_procx*n_procy-1
            call gc_rrecv(100*d, global_rows*global_row_length,         &
     &                    d, info, rbuf, rbuf)
            Do j = 1, g_rows(d)
              gj = g_datastart(2,d) + j - 1
              Do i = 1, g_row_length(d)
                gi = g_datastart(1,d) + i - 1
                global_real(gi,gj) = rbuf(i,j)
              End Do
            End Do
          End Do
          Write(dun) global_real
        End If

      End Do

! end of routine

      IF (lhook) CALL dr_hook('WRITE_DUMP_REAL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE write_dump_real

