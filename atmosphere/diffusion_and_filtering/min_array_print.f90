! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! Subroutine Min_array_print

      Subroutine Min_array_print(                                       &
     &                           field,                                 &
     &                           row_length, rows, levels,              &
     &                           halo_i, halo_j,                        &
     &                           start_level, end_level,                &
     &                           global_row_length, global_rows,        &
     &                           l_datastart, me, n_proc )

! Purpose:
!          Print value and location of minimum value of an array
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!

      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook

      Implicit None

      Integer , Intent(IN) ::                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, levels                                                          &
                         ! number of levels in field
     &, start_level                                                     &
                         ! start level for finding min
     &, end_level                                                       &
                         ! end level for finding min
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, l_datastart(2)                                                  &
                             ! First gridpoints held by this processor
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, n_proc                                                          &
                         ! Total number of processors
     &, me               ! This processor number

      Real , Intent(IN) ::                                              &
     &       field (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,   &
     &              levels)

! Local Variables.

      Integer                                                           &
     &  k, ki, kr                                                       &
                    ! Loop counters
     &, num_levels                                                      &
     &, gi, gj                                                          &
                    ! global pointers
     &, pe                                                              &
                    ! processor identifier
     &, info

      Real                                                              &
     &  sumr(2,levels)                                                  &
                              ! array  for summing integers over pe's
     &, min_pe(levels)                                                  &
                              ! min wind each level on this pe
     &, min_all(levels)       ! array  for finding min over pe's
 
      Integer                                                           &
     &  sumi(3,levels)                                                  &
                            ! array  for summing integers over pe's
     &, i_min(levels)                                                   &
                            ! i pointers for prints
     &, j_min(levels)       ! j pointers for prints
 
      Real                                                              &
     &  recip_row_length                                                &
     &, recip_rows                                                      &
     &, lambda                                                          &
     &, phi

      CHARACTER(LEN=8)                                                       &
     &  l_string                                                        &
     &, p_string

      Integer, dimension(2) :: min_indices

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 0. Initialise values
! ----------------------------------------------------------------------

        IF (lhook) CALL dr_hook('MIN_ARRAY_PRINT',zhook_in,zhook_handle)

        recip_row_length = 1.0 / real(global_row_length)
        recip_rows = 1.0 / real(global_rows)
        p_string ='% North'
        l_string ='% East'
        num_levels = end_level - start_level + 1

        Do k = 1, num_levels
          min_indices = minloc(field(1:row_length,1:rows,k))
          i_min(k) = min_indices(1)
          j_min(k) = min_indices(2)
          min_pe(k) = field(min_indices(1),min_indices(2),k)       
          min_all(k) = min_pe(k)       
        End Do !  k = 1, num_levels
 
        call gc_rmin(num_levels, n_proc, info, min_all)
        
        sumi = 0
        sumr = 0.0
        ki = 0
        kr = 0
        Do k = 1, num_levels
          ki = ki + 3
          kr = kr + 2
          if( min_pe(k) <= min_all(k) )then
            sumi(1,k) = me ! min is on this processor for this level
            gi = l_datastart(1) + i_min(k) - 1
            gj = l_datastart(2) + j_min(k) - 1
            sumi(2,k) = gi
            sumi(3,k) = gj
            sumr(1,k) = (gi-1) * recip_row_length * 100.0
            sumr(2,k) = (gj-1) * recip_rows * 100.0
          endif ! min_pe(k) <= min_all(k)
        End Do !  k = 1, num_levels

        call gc_isum(ki, n_proc, info, sumi)
        call gc_rsumr(kr, n_proc, info, sumr)
        
        write(6,*) '                         global position'
        write(6,*) 'level   min    proc   i   ,  j '
        Do   k = 1, num_levels
          pe = sumi(1,k)
          gi = sumi(2,k)
          gj = sumi(3,k)
          lambda = sumr(1,k)
          phi    = sumr(2,k)
          write(6,'(I4, E12.4, 3I5, F6.1, A7, F6.1, A7)')               &
     &                    k, min_all(k), pe, gi, gj,                    &
     &                    lambda, l_string, phi, p_string
          end Do   !  k = 1, num_levels

        IF (lhook) CALL dr_hook('MIN_ARRAY_PRINT',zhook_out,zhook_handle)   
        RETURN
      END SUBROUTINE MIN_ARRAY_PRINT

