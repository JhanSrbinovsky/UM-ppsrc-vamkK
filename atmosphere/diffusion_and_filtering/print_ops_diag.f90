! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT*****************************
!
! Subroutine Print_ops_diag

      Subroutine Print_ops_diag(                                        &
     &                     w, theta, row_length, rows,                  &
     &                     model_levels, model_domain,                  &
     &                     global_row_length, global_rows,              &
     &                     off_x, off_y, halo_i, halo_j,                &
     &                     me, n_proc, proc_row_group,                  &
     &                     timestep_number, print_step, diag_interval,  &
     &                     rpemax, rpemin, ipesum,                      &
     &                     L_print_pe, L_print_wmax, L_print_theta1,    &
     &                     max_w_run, min_theta1_run,                   &
     &                     time_w_max, time_theta1_min )

! Purpose:
!          Diagnostic routine for w and level 1 theta
!
! Method:
!          Is described in ;
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! History:
! Version   Date      Comment
! ----     -------    -------
!  6.2   28/02/05     Original code   Terry Davies
!   Cut down version of print_diag for operational use
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.

      Integer , Intent(IN) ::                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, proc_row_group                                                  &
     &, n_proc                                                          &
                   ! Total number of processors
     &, me                                                              &
     &, model_domain                                                    &
                         ! holds integer code for model domain
     &, timestep_number                                                 &
     &, print_step                                                      &
     &, diag_interval                                                   &
! the following 3 sizes are set in SETCONA
     &, rpemax                                                          & 
                  ! total size of array needed to find max over pe's
     &, rpemin                                                          & 
                  ! total size of array needed to find min over pe's
     &, ipesum    ! total size of integer array needed to sum over pe's

      Integer , Intent(INOUT) ::                                        &
     &  time_theta1_min                                                 &
     &, time_w_max(model_levels)

      Real , Intent(INOUT) ::                                           &
     &  min_theta1_run                                                  &
     &, max_w_run(model_levels)

      Real , Intent(IN) ::                                              &
     &  w(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &          0:model_levels)                                         &
     &, theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &          model_levels)

      Logical , Intent(IN) ::                                           &
     &  L_print_wmax                                                    &
                         ! Print control
     &, L_print_theta1                                                  &
                         ! Print control
     &, L_print_pe       ! print diagnostics on all pe's if true

! Local Variables.

      Integer                                                           &
     &  i, j, k                                                         &
                     ! Loop counters
     &, ki                                                              & 
                ! pointers for arrays summed over pe's
     &, kminr                                                           & 
                   ! pointers for arrays summed over pe's
     &, kmaxr                                                           & 
                   ! pointers for arrays summed over pe's
     &, info                                                            &
     &, pe                                                              &
     &, pe_max                                                          &
     &, level_max                                                       &
     &, level_max_run                                                   &
     &, time_max_run

      Real                                                              &
     &  min_theta_pe                                                    &
                     ! min theta level 1 on this pe
     &, max_w                                                           &
     &, max_run

! Local arrays

      Integer                                                           &
     &  sumi(ipesum) ! array  for summing integers over pe's

      Real                                                              &
     &  max_w_pe(model_levels - 1)                                      &
                                   ! max w each level on this pe
     &, max_real(rpemax)                                                & 
                         ! array  for finding max over pe's
     &, min_real(rpemin) ! array  for finding min over pe's

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! No External Routines:

      IF (lhook) CALL dr_hook('PRINT_OPS_DIAG',zhook_in,zhook_handle)
      If( mod( timestep_number , diag_interval )  ==  0) then

! ----------------------------------------------------------------------
! Section 1  !  Initialise counters
! ----------------------------------------------------------------------
      kminr = 0
      kmaxr = 0

! ----------------------------------------------------------------------
! Section 2  Find theta min on level 1
! ----------------------------------------------------------------------
      if ( L_print_theta1 ) then
! copy into min_real(1) which will hold global min later
        min_theta_pe = minval(theta(1:row_length,1:rows,1))
        kminr = kminr + 1
        min_real(kminr) = min_theta_pe
      endif ! L_print_theta1

! ----------------------------------------------------------------------
! Section 3   find max w for each level
!               w(model_levels) = 0.0 everywhere so omit
! ----------------------------------------------------------------------
      if ( L_print_wmax ) then

        Do   k = 1, model_levels - 1
! Copy local max w into w_max which will hold global w_max after gc_rmax
          max_w_pe(k) = maxval(w(1:row_length,1:rows,k))
          kmaxr = kmaxr + 1
          max_real(kmaxr) = max_w_pe(k)
        End Do  ! k = 1, model_levels - 1

      endif ! L_print_wmax

! ----------------------------------------------------------------------
! Section 4. Now find max/mins over all processors
!            All required fields done at same time
! ----------------------------------------------------------------------

      if ( kmaxr > 0 ) then
        call gc_rmax(kmaxr, n_proc, info, max_real)
      endif !  kmaxr > 0
      if ( kminr > 0 ) then
        call gc_rmin(kminr, n_proc, info, min_real)
      endif !  kminr > 0

! ----------------------------------------------------------------------
! Section 5  Obtain max and mins from max_real, min_real
!             and fill sumi arrays for summing over pe's
! ----------------------------------------------------------------------

!  Initialise pointers in arrays for summing over pe's
      ki = 0
!  Re-initialise pointers in max/min arrays
      kminr = 0
      kmaxr = 0

      if ( L_print_theta1 ) then
        kminr = kminr + 1
        if( min_real(kminr) < min_theta1_run )then
            min_theta1_run = min_real(kminr)
            time_theta1_min = timestep_number
        endif ! min_real(kminr) < min_theta1_run
        ki = ki + 1
        sumi(ki) = 0
        if( min_theta_pe <= min_real(kminr) )then
          sumi(ki) = me ! min is on this processor
        endif ! min_theta_pe <= min_real(kminr)

      endif ! L_print_theta1

      if ( L_print_wmax ) then
        Do   k = 1, model_levels - 1
          kmaxr = kmaxr + 1
          if( max_real(kmaxr) > max_w_run(k) )then
             max_w_run(k) =  max_real(kmaxr)
             time_w_max(k) = timestep_number
          endif ! max_real(kmaxr) > max_w_run(k)
          ki = ki + 1
          sumi(ki) = 0
          if ( max_w_pe(k) >= max_real(kmaxr)) then
            sumi(ki) = me ! max is on this processor for this level
          endif  ! max_w_pe(k) <= max_real(kmaxr)
        End do  ! k = 1, model_levels - 1
      endif ! L_print_wmax

! ----------------------------------------------------------------------
!  Printing will be done in section 7
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Section 6  Summing over pe's to obtain sums and location of max/mins
! ----------------------------------------------------------------------

      if( ki > 0 ) call gc_isum(ki, n_proc, info, sumi)

      EndIf ! mod( timestep_number , diag_interval )  ==  0

! ----------------------------------------------------------------------
! Section 7. Print diagnostic information
! ----------------------------------------------------------------------

      If( mod( timestep_number , print_step ) == 0                      &
     &   .and.   (L_print_pe .or. me == 0)       ) then

! Re-initialise pointers
        kminr = 0
        kmaxr = 0
        ki = 0

        if ( L_print_theta1 ) then
          kminr = kminr + 1
          ki = ki + 1
          pe = sumi(ki)
          write(6,*) '  '
         write(6,*)'Minimum theta level 1 for timestep ',timestep_number
          write(6,*) '           This timestep',                        &
     &             '           This run'
          write(6,*) '  Min theta1     proc   ',                        &
     &            '   Min theta1 timestep'
          write(6,990) min_real(kminr), pe,                             &
     &               min_theta1_run, time_theta1_min

        endif ! L_print_theta1

        if ( L_print_wmax ) then
          max_w = -100.0
          max_run = -100.0
          Do   k = 1, model_levels - 1
            kmaxr = kmaxr + 1
            ki = ki + 1
            pe = sumi(ki)
            if ( max_w < max_real(kmaxr) ) then
              max_w = max_real(kmaxr)
              level_max = k
              pe_max = pe
            endif ! max_w < max_real(kmaxr)
            if ( max_run < max_w_run(k) ) then
              max_run = max_w_run(k)
              level_max_run = k
              time_max_run = time_w_max(k)
            endif ! max_run < max_w_run(k)
          endDo   !  k = 1, model_levels - 1
          write(6,*) '  '
          write(6,*) ' Maximum vertical velocity at timestep ',         &
     &            timestep_number,'      Max w this run '
          write(6,*) '   w_max   level  proc  ',                        &
     &           '     run w_max level timestep'
          write(6,991) max_w, level_max, pe_max,                        &
     &               max_run, level_max_run, time_max_run
          endif ! L_print_wmax

      endif !  mod(timestep_number,print_step) == 0

! ----------------------------------------------------------------------
! Section 8. Print formats
! ----------------------------------------------------------------------

 990   FORMAT(1X, F11.3, I8, F11.3, I6)
 991   FORMAT(1X, E11.3, I4, I7, E11.3, I5, I6)

      IF (lhook) CALL dr_hook('PRINT_OPS_DIAG',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Print_ops_diag

