! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine tardiff_q_wss
      subroutine tardiff_q_wss(                                         &
     &                         q, w_adv,                                &
     &                         r_theta_levels, r_rho_levels,            &
     &                         sec_theta_latitude,                      &
     &                         cos_theta_latitude, cos_v_latitude,      &
     &                         off_x, off_y, halo_i, halo_j,            &
     &                         halo_i_star, halo_j_star,                &
     &                         at_extremity, proc_row_group,            &
     &                         delta_lambda, delta_phi,                 &
     &                         timestep, rows, n_rows, row_length,      &
     &                         model_levels, wet_model_levels,          &
     &                         model_domain, global_row_length,         &
     &                         q_star, w_limit, horizontal_level,       &
     &                         factor, test_level,                      &
     &                         start_level, end_level,                  &
     &                         L_diag_w, w_local_mask)

! Purpose:
!          Calculates conservative horizontal diffusion increment to q
!          subject to w > w_limit, and has a steep slope test
!
! Method:
!          Is described in ;
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!


      USE global_2d_sums_mod, ONLY : global_2d_sums
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
! Modules for z-level diffusion
      USE atm_fields_bounds_mod
      USE level_heights_mod,    ONLY : eta_theta_levels, eta_rho_levels 
      USE pofil_zlevel_mod
      USE conversions_mod,      ONLY: pi_over_180
      USE proc_info_mod,        ONLY: me, n_procy, l_datastart, n_proc
      USE PrintStatus_mod,      ONLY: printstatus, prstatus_diag

      Implicit None

! Arguments with Intent IN. ie: Input variables.

      Logical                                                           &
     &  at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
     &, L_diag_w         ! diagnostic control

      Integer                                                           &
     &  row_length                                                      &
                         ! number of point on a row.
     &, rows                                                            &
                         ! number of rows.
     &, n_rows                                                          &
                         ! number of v rows.
     &, model_levels                                                    &
                         ! number of model levels.
     &, wet_model_levels                                                &
                             ! number of model levels.
     &, halo_i_star                                                     &
                      ! Size of halo in i direction for star field.
     &, halo_j_star                                                     &
                      ! Size of halo in j direction for star field.
     &, halo_i                                                          &
                      ! Size of halo in i direction.
     &, halo_j                                                          &
                      ! Size of halo in j direction.
     &, off_x                                                           &
                      ! Size of small halo in i
     &, off_y                                                           &
                      ! Size of small halo in j.
     &, global_row_length                                               &
     &, proc_row_group

      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, test_level                                                      &
     &, start_level                                                     &
     &, end_level                                                       &
     &, horizontal_level   !  upper test level for steep slopes

      Real                                                              &
     &  timestep                                                        &
     &, w_limit                                                         &
                 ! Vertical velocity test value
     &, factor   ! effective diffusion coefficient

      Real                                                              &
           ! vertical co-ordinate arrays.
     &  r_theta_levels (1-halo_i:row_length+halo_i,                     &
     &                  1-halo_j:rows+halo_j, 0:model_levels)           &
     &, r_rho_levels (1-halo_i:row_length+halo_i,                       &
     &                1-halo_j:rows+halo_j, model_levels)

      Real                                                              &
           !  trigonometric functions
     &  sec_theta_latitude (1-off_x:row_length+off_x,                   &
     &                    1-off_y:rows+off_y)                           &
     &, cos_theta_latitude (1-off_x:row_length+off_x,                   &
     &                    1-off_y:rows+off_y)                           &
     &, cos_v_latitude (1-off_x:row_length+off_x,                       &
     &                    1-off_y:n_rows+off_y)

      Real                                                              &
     &  q (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j             &
     &,         wet_model_levels)                                       &
     &, w_adv (1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j         &
     &,        0:model_levels)

      Real                                                              &
           ! horizontal co-ordinate spacing.
     &  delta_lambda                                                    &
     &, delta_phi

! Arguments with Intent IN/OUT. ie: Input variables changed on Output.

      Real                                                              &
     &  q_star (1-halo_i_star:row_length+halo_i_star,                   &
     &          1-halo_j_star:rows+halo_j_star, wet_model_levels)       &
     &, w_local_mask(row_length,rows)

! Local Variables.

      Integer                                                           &
     &  i, j, k, ka                                                     &
                                  ! Loop indices
     &, j_start, j_stop                                                 &
                              ! Loop indices
     &, info                                                            &
     &, active_levels                                                   &
     &, level_flat

      Real                                                              &
     &  recip_delta_phi

! Local arrays

      Real                                                              &
     &  lambda_term(row_length, rows, model_levels )                    &
     &, phi_term(row_length, rows, model_levels )                       &
     &, delta_z(1-off_x:row_length+off_x,                               &
     &          1-off_y:rows+off_y, model_levels )                      &
     &, temp(1-off_x:row_length+off_x, 1-off_y:rows+off_y)              &
     &, l_s_poles(row_length, model_levels )                            &
     &, l_n_poles(row_length, model_levels )                            &
     &, sum_s(model_levels )                                            &
     &, sum_n(model_levels )                                            &
     &, diffc_i(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, diffc_j(1-off_x:row_length+off_x, 1-off_y:rows+off_y)           &
     &, mask_i(1-off_x:row_length, rows, model_levels)                  &
     &, mask_j(row_length, 1-off_y:rows, model_levels)                  &
     &, w_max(1-off_x:row_length+off_x, 1-off_y:rows+off_y)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      INTEGER ::                                                        &
                           ! array size for u_begin etc
        j_q_sweeps_begin(0:2),                                          &
                           ! row pointers for 1-2-1 filter
        j_q_sweeps_end(0:2),                                            &
        q_sweeps(2)

      REAL :: xi1_p(1-halo_i:row_length+halo_i),                        &
              xi1_u(1-halo_i:row_length+halo_i),                        &
              xi2_p(1-halo_j:rows+halo_j),                              &
              xi2_v(1-halo_j:n_rows+1+halo_j)

! z level test arrays and diagnostic
      REAL ::                                                           &
        diff_temp_q(1-off_x:row_length+off_x, 1-off_y:rows+off_y,       &
                    wet_model_levels) 

      REAL :: flux

      REAL :: csphi_theta(1-halo_j:rows+halo_j),                        &
              csphi_v(1-halo_j:n_rows+1+halo_j)

      REAL,DIMENSION(:,:,:),ALLOCATABLE ::                              &
     &  r_theta_at_u                                                    &
     &, r_theta_at_v


      Integer :: istat
      Real    :: k_sum(row_length,rows),  k_sum_v(row_length,n_rows) 
      Real    :: row_sum(rows) 

      INTEGER :: gj, gi, count
      REAL    :: base_lambda, base_phi

! No External Routines:

      IF (lhook) CALL dr_hook('TARDIFF_Q_WSS',zhook_in,zhook_handle)
      active_levels = end_level - start_level + 1
      recip_delta_phi = 1.0/ delta_phi

      j_start = 1
      j_stop = rows
      If (model_domain  ==  1) Then
        If (at_extremity(PSouth)) j_start = 2
        If (at_extremity(PNorth)) j_stop = rows - 1
      End If

! initialise w_local_mask to zero

      If (L_diag_w) Then
        Do j = 1,rows
          Do i = 1, row_length
            w_local_mask(i,j) = 0.0
          End Do
        End Do
      EndIf !  L_diag_w

! ----------------------------------------------------------------------
! Section 1.   Calculate delta_z and diffusion coefficients
! ----------------------------------------------------------------------
! Originally calculate D(r)/D(eta) about q levels BUT
! instead calculate delta_z about q levels
! We can omit the d(eta) part since this is constant and can be
! brought out through the derivatives and cancelled


      Do k = 1, active_levels
        ka = k + start_level - 1
        If(k < model_levels) then
          Do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
               delta_z(i,j,k) = (r_rho_levels(i,j,ka+1) -               &
     &                           r_rho_levels(i,j,ka) )
            End Do
          End Do
        Else ! k = model_levels
! can use any constant value for delta_z since it will cancel
          Do j = 1-off_y, rows+off_y
            Do i = 1-off_x, row_length+off_x
              delta_z(i,j,k) = 1.0
            End Do
          End Do
        EndIf ! k < model_levels
      End Do   !  k = 1, active_levels


  count = 0

!  Initialise diffusion coefficients to zero
      Do j = 1-off_y, rows+off_y
        Do i = 1-off_x, row_length+off_x
          diffc_i(i,j) = 0.0
          diffc_j(i,j) =   0.0
          w_max(i,j) = 0.0
        Enddo  ! i = 1-off_x, row_length+off_x
      Enddo    ! j = 1-off_y, rows+off_y

!  Check to see if w > w_print_limit at each point
       Do k =  test_level, end_level
         Do j = 1-off_y, rows+off_y
           Do i = 1-off_x, row_length+off_x
!  Find if vertical velocity above threshold at this point
!  start at level 5 since delta_z small near surface
              if( w_max(i,j) < w_adv(i,j,k) ) then
                  w_max(i,j) = w_adv(i,j,k)
              endif ! w_max(i,j) < w_adv(i,j,k
            Enddo  ! i = 1-off_x, row_length+off_x
          Enddo    ! j = 1-off_y, rows+off_y
        EndDo  !   k =  test_level, end_level

         Do j = 1-off_y, rows+off_y
           Do i = 1-off_x, row_length+off_x
             if( w_max(i,j) > w_limit)then
               diffc_i(i,j) = factor
               diffc_j(i,j) = factor
             endif ! w_max(i,j) > w_limit
           Enddo  ! i = 1-off_x, row_length+off_x
         Enddo    ! j = 1-off_y, rows+off_y
         
       IF (printstatus >= prstatus_diag) THEN 
       ! Compute diagnostic and output
         Do j = 1, rows
           Do i = 1, row_length
             if( w_max(i,j) > w_limit)then
               count = count + 1
             endif ! w_max(i,j) > w_limit
           Enddo  ! i = 1-off_x, row_length+off_x
         Enddo    ! j = 1-off_y, rows+off_y
        CALL gc_isum(1,n_proc,istat,count)        
        IF ( me == 0) WRITE(6,*) 'targeted diffusion in  ',count,' columns'
       END IF


        Do j = 1, rows+off_y
          Do i = 1, row_length+off_x
           if( w_max(i,j) > w_max(i-1,j))then
             diffc_i(i-1,j) = diffc_i(i,j)
           endif ! w_max(i,j) > w_max(i-1,j)
           if( w_max(i,j) > w_max(i,j-1))then
             diffc_j(i,j-1) = diffc_j(i,j)
           endif ! w_max(i,j) > w_max(i,j-1)
          Enddo  !i = 1, row_length+off_x
        Enddo    !j = 1, rows+off_y

      If (L_diag_w) Then
        Do j = 1,rows
          Do i = 1, row_length
            if( w_max(i,j) > w_limit)then
                w_local_mask(i,j) = 1.0
            endif ! w_max(i,j) > w_limit
          End Do
        End Do
      EndIf !  L_diag_w


! ----------------------------------------------------------------------
! Section 2.   Switch off theta diffusion at steep slopes
! ----------------------------------------------------------------------
!  Horizontal_level should be the first flat level and as such
!  does not need testing for a sloping surface.
!   i.e if level_flat < 2   then
!       diffusion behaves as original (along eta surfaces)
! To switch off slope test completetely then
! either hard-wire level_flat = 1 or ensure that the argument
!    horizontal level <= wet_model_levels - - model_levels + 1
!      level_flat = 1
       level_flat =  horizontal_level
       if (level_flat  <=  0) level_flat = 1

       Do k = 1, active_levels
        Do j = 1, rows
          Do i = 0, row_length
           mask_i(i,j,k) = diffc_i(i,j)
          Enddo
        Enddo
        Do j = 0, rows
          Do i = 1, row_length
           mask_j(i,j,k) = diffc_j(i,j)
          Enddo
        Enddo
      Enddo   !  k = 1, active_levels

      Do k = start_level, level_flat - 1
         ka = k - start_level + 1
        Do j = 1, rows
          Do i = 0, row_length
            if( r_theta_levels(i,j,k) < r_theta_levels(i+1,j,k-1) .or.  &
     &          r_theta_levels(i+1,j,k) < r_theta_levels(i,j,k-1)) then
              mask_i(i,j,ka) = 0.0
           endif
          Enddo
        Enddo
        Do j = 0, rows
          Do i = 1, row_length
            if( r_theta_levels(i,j,k) < r_theta_levels(i,j+1,k-1) .or.  &
     &          r_theta_levels(i,j+1,k) < r_theta_levels(i,j,k-1)) then
              mask_j(i,j,ka) = 0.0
           endif
          Enddo
        Enddo
      Enddo    ! ka = start_level, level_flat - 1

! ----------------------------------------------------------------------
! Section 2.1  Calculate lambda direction term.
! ----------------------------------------------------------------------

        Do k = 1, active_levels
           ka = k + start_level - 1
            Do j = 1, rows
              Do i = 1-off_x, row_length
                temp(i,j) = (q(i+1,j,ka) * delta_z(i+1,j,k) -           &
     &                        q(i  ,j,ka) * delta_z(i ,j,k) ) *         &
     &                                          mask_i(i,j,k) *         &
     &             r_theta_levels(i,j,ka) * r_theta_levels(i,j,ka)
              End Do
              Do i = 1, row_length
                lambda_term(i,j,k) = temp(i,j) - temp(i-1,j)
              End Do
            End Do

! ----------------------------------------------------------------------
! Section 2.2  Calculate phi direction term.
! ----------------------------------------------------------------------

             Do j = j_start-1, j_stop
              Do i = 1, row_length
                temp(i,j) = (q(i,j+1,ka) * delta_z(i,j+1,k) -           &
     &                        q(i,j,  ka) * delta_z(i,j,  k) ) *        &
     &                     mask_j(i,j,k) * cos_v_latitude(i,j) *        &
     &               r_theta_levels(i,j,ka) * r_theta_levels(i,j,ka)
              End Do
             End Do
             Do j = j_start, j_stop
              Do i = 1, row_length
                phi_term(i,j,k) = (temp(i,j) - temp(i,j-1)) *           &
     &                              sec_theta_latitude(i,j)
              End Do
             End Do

        If ( model_domain  ==  mt_Global ) Then

              If(at_extremity(PSouth))then
                Do i = 1, row_length
             l_s_poles(i,k) = (q(i,2,ka) * delta_z(i,2,k) -             &
     &                          q(i,1,ka) * delta_z(i,1,k)) *           &
     &                      mask_j(i,1,k) * cos_v_latitude(i,1) *       &
     &             r_theta_levels(i,1,ka) * r_theta_levels(i,1,ka)
                End Do
              End If

              If(at_extremity(PNorth))then
                Do i = 1, row_length
            l_n_poles(i,k) = (q(i,rows-1,ka) * delta_z(i,rows-1,k) -    &
     &                         q(i,rows,ka) * delta_z(i,rows,k)) *      &
     &                   mask_j(i,rows-1,k) * cos_v_latitude(i,rows-1) *&
     &         r_theta_levels(i,rows-1,ka) * r_theta_levels(i,rows-1,ka)
                End Do
              End If

            End If !model_domain  ==  mt_Global

        EndDo ! k = 1, active_levels

        If ( model_domain  ==  mt_Global ) Then

          If (at_extremity(PSouth)) Then
            CALL global_2d_sums(l_s_poles, row_length, 1, 0, 0,         &
                                active_levels, sum_s,                   &
                                proc_row_group)
            Do k = 1, active_levels
              sum_s(k) = sum_s(k) * 8. * recip_delta_phi /              &
     &                                   global_row_length
              Do i = 1, row_length
                phi_term(i,1,k) = sum_s(k)
              End Do
            End Do   ! k = 1, active_levels
          End If  ! at_extremity(PSouth)

          If (at_extremity(PNorth)) Then
            CALL global_2d_sums(l_n_poles, row_length, 1, 0, 0,         &
                                active_levels, sum_n,                   &
                                proc_row_group)
            Do k = 1, active_levels
              sum_n(k) = sum_n(k) * 8. * recip_delta_phi /              &
     &                                   global_row_length
              Do i = 1, row_length
                phi_term(i,rows,k) = sum_n(k)
              End Do
            End Do  ! k = 1, active_levels
          End If ! at_extremity(PNorth)

        End If !model_domain  ==  mt_Global

! ----------------------------------------------------------------------
! Section 3.   Diffusion for q
! ----------------------------------------------------------------------

      Do k = 1, active_levels
         ka = k + start_level - 1
          Do j = 1, rows
            Do i = 1, row_length
              q_star(i,j,ka) = q_star(i,j,ka) +                         &
     &                        ( lambda_term(i,j,k) + phi_term(i,j,k) ) /&
     &             (r_theta_levels(i,j,ka) * r_theta_levels(i,j,ka) *   &
     &                                              delta_z(i,j,k) )
            End Do
          End Do
      End Do   ! k = 1, active_levels


      IF (lhook) CALL dr_hook('TARDIFF_Q_WSS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE tardiff_q_wss

