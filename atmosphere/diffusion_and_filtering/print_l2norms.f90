! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! subroutine print_l2norms
      subroutine print_l2norms(                                         &
     &                         exner_rho_levels, rho,                   &
     &                         u, v, w,                                 &
     &                         u_adv, v_adv, w_adv,                     &
     &                         theta, q, qcl, qcf,                      &
     &                         R_u, R_v, R_w, theta_star,               &
     &                         q_star, qcl_star, qcf_star,              &
     &                         row_length, rows, n_rows,                &
     &                         model_levels, wet_levels,                &
     &                         start_level, end_level,                  &
     &                         offx, offy, halo_i, halo_j,              &
     &                         me, nprocx, nprocy,                      &
     &                         proc_row_group, proc_col_group,          &
     &                         global_row_length, global_rows,          &
     &                         at_extremity, datastart, model_domain,   &
     &                         L_full, L_do_halos, L_print_pe )

! Purpose:
!          To calculate and print l2norms of tendencies or star fields
!
! Method:
!          Is described in ;
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Diffusion and Filtering
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.
!


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

! Arguments with Intent IN. ie: Input variables.
      INTEGER, Intent(In) ::                                            &
     &  row_length                                                      &
                          ! in: no of points per local row
     &, rows                                                            &
                          ! in: no of local (theta) rows
     &, n_rows                                                          &
                          ! in: no of local (v) rows
     &, model_levels                                                    &
                          ! in: no of model levels
     &, wet_levels                                                      &
                          ! in: no of moist-levels
     &, start_level                                                     &
                          ! start level for norm calculation
     &, end_level                                                       &
                          ! end level for norm calculation
     &, offx                                                            &
                          ! standard halo size in east-west
     &, offy                                                            &
                          ! standard halo size in North-South
     &, halo_i                                                          &
                          ! extended halo size in East-West
     &, halo_j                                                          &
                          ! extended halo size in North-South
     &, model_domain                                                    &
     &, me                                                              &
                          ! pe id
     &, nprocx                                                          &
     &, nprocy                                                          &
     &, proc_row_group                                                  &
     &, proc_col_group                                                  &
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, datastart(3)

      Real, intent(in) ::                                               &
     &  u(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)       &
     &, v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels)     &
     &, w(1-offx:row_length+offx, 1-offy:rows+offy,                     &
     &      0:model_levels)                                             &
     &, u_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &          model_levels)                                           &
     &, v_adv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,       &
     &          model_levels)                                           &
     &, w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &          0:model_levels)                                         &
     &, theta(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)   &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &      wet_levels)                                                 &
     &, qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,           &
     &        wet_levels)                                               &
     &, rho(1-offx:row_length+offx, 1-offy:rows+offy, model_levels)     &
     &, exner_rho_levels(1-offx:row_length+offx,                        &
     &                   1-offy:rows+offy, model_levels+1)

      Real, intent(in) ::                                               &
     &  R_u(1-offx:row_length+offx, 1-offy:rows+offy,                   &
     &        model_levels)                                             &
     &, R_v(1-offx:row_length+offx, 1-offy:n_rows+offy,                 &
     &        model_levels)                                             &
     &, R_w(row_length, rows, model_levels - 1)                         &
     &, theta_star(1-offx:row_length+offx,                              &
     &               1-offy:rows+offy, model_levels)                    &
     &, q_star(1-offx:row_length+offx,                                  &
     &           1-offy:rows+offy, wet_levels)                          &
     &, qcl_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)                        &
     &, qcf_star(1-offx:row_length+offx,                                &
     &             1-offy:rows+offy, wet_levels)

      Logical                                                           &
     &  L_do_halos                                                      &
     &, L_full                                                          &
                  ! T if norms of prognostics, F if norms of increments
     &, L_print_pe      ! true if  printing on all pe's

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

!    local variables

      Integer :: first_wet
      Integer :: last_wet
      Integer :: p_end_level
      Integer :: w_end_level

      Real ::  Two_norm
      Logical  :: L_do_w

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('PRINT_L2NORMS',zhook_in,zhook_handle)

      If ( start_level > wet_levels ) Then
        first_wet = wet_levels
        last_wet = wet_levels
      Else If ( end_level > wet_levels ) Then
        first_wet = start_level
        last_wet = wet_levels
      Else
        first_wet = start_level
        last_wet = end_level
      End If ! start_level > wet_levels
      
      p_end_level = end_level
      if ( end_level == model_levels ) p_end_level = model_levels + 1

      L_do_w = .true.
      w_end_level = end_level
      If ( end_level == model_levels ) Then
        w_end_level = model_levels - 1
        If ( start_level == model_levels ) L_do_w = .false.
      End If ! end_level == model_levels
      
      If (L_full) Then   !  Norms of prognostic variables

! DEPENDS ON: two_norm_levels
        Call Two_Norm_levels(                                           &
     &                     exner_rho_levels, row_length, rows,          &
     &                     model_levels + 1, start_level, p_end_level,  &
     &                     model_domain, offx, offy, 0, 0,              &
     &                     at_extremity, nprocx, nprocy,                &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     datastart, .false., .true., 1, Two_Norm )
        if ( L_print_pe .or. me == 0 ) then
          write(6,*)' L2 Norm of prognostic variables '
          write(6,*)'Levels',start_level,' to', p_end_level ,           &
     &           ' Exner pressure Two_Norm = ' , Two_Norm
        endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
        Call Two_Norm_levels(                                           &
     &                       rho, row_length, rows, model_levels,       &
     &                       start_level, end_level,                    &
     &                       model_domain, offx, offy, 0, 0,            &
     &                       at_extremity, nprocx, nprocy,              &
     &                       global_row_length, global_rows,            &
     &                       proc_col_group, proc_row_group,            &
     &                       datastart, .false., .true., 1, Two_Norm )
        if ( L_print_pe .or. me == 0 ) then
          write(6,*)'Levels',start_level,' to',end_level,               &
     &         ' rho Two_Norm = ' , Two_Norm
        endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
        Call Two_Norm_levels(                                           &
     &                     theta, row_length, rows, model_levels,       &
     &                     start_level, end_level,                      &
     &                     model_domain, offx, offy, 0, 0,              &
     &                     at_extremity, nprocx, nprocy,                &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     datastart, .false., .true., 1, Two_Norm )
        if ( L_print_pe .or. me == 0 ) then
          write(6,*)'Levels',start_level,' to',end_level,               &
     &         ' theta Two_Norm = ' , Two_Norm
        endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
        Call Two_Norm_levels(                                           &
     &                     q, row_length, rows, wet_levels,             &
     &                     first_wet, last_wet,                         &
     &                     model_domain, halo_i, halo_j, 0, 0,          &
     &                     at_extremity, nprocx, nprocy,                &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     datastart, .false., .true., 1, Two_Norm )
        if ( L_print_pe .or. me == 0 ) then
          write(6,*)'Levels',first_wet,' to',last_wet,                  &
     &         ' q Two_Norm = ' , Two_Norm
        endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
        Call Two_Norm_levels(                                           &
     &                     u, row_length, rows, model_levels,           &
     &                     start_level, end_level,                      &
     &                     model_domain, offx, offy, 1, 0,              &
     &                     at_extremity, nprocx, nprocy,                &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     datastart, .false., .true., 1, Two_Norm )
        if ( L_print_pe .or. me == 0 ) then
          write(6,*)'Levels',start_level,' to',end_level,               &
     &           ' u Two_Norm = ' , Two_Norm
        endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
        Call Two_Norm_levels(                                           &
     &                     v, row_length, n_rows, model_levels,         &
     &                     start_level, end_level,                      &
     &                     model_domain, offx, offy, 0, 1,              &
     &                     at_extremity, nprocx, nprocy,                &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     datastart, .false., .true., 1, Two_Norm )
        if ( L_print_pe .or. me == 0 ) then
          write(6,*)'Levels',start_level,' to',end_level,               &
     &           ' v Two_Norm = ' , Two_Norm
        endif ! L_print_pe .or. me == 0

        If (L_do_w) Then 
! DEPENDS ON: two_norm_levels
        Call Two_Norm_levels(                                           &
     &                     w, row_length, rows, model_levels + 1,       &
     &                     start_level, p_end_level,                    &
     &                     model_domain, offx, offy, 0, 0,              &
     &                     at_extremity, nprocx, nprocy,                &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     datastart, .false., .true., 1, Two_Norm )
        if ( L_print_pe .or. me == 0 ) then
          write(6,*)'Levels',start_level,' to', p_end_level,            &
     &           ' w Two_Norm = ' , Two_Norm
        endif ! L_print_pe .or. me == 0
        End If ! L_do_w

      Else   !  Norms of increments or star fields

! DEPENDS ON: two_norm_levels
        Call Two_Norm_levels(                                           &
     &                     theta_star, row_length, rows, model_levels,  &
     &                     start_level, end_level,                      &
     &                     model_domain, offx, offy, 0, 0,              &
     &                     at_extremity, nprocx, nprocy,                &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     datastart, .false., .true., 1, Two_Norm )
        if ( L_print_pe .or. me == 0 ) then
          write(6,*)'Levels',start_level,' to',end_level,               &
     &           ' theta_star Two_Norm = ' , Two_Norm
        endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
        Call Two_Norm_levels(                                           &
     &                     q_star, row_length, rows, wet_levels,        &
     &                     first_wet, last_wet,                         &
     &                     model_domain, offx, offy, 0, 0,              &
     &                     at_extremity, nprocx, nprocy,                &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     datastart, .false., .true., 1, Two_Norm )
        if ( L_print_pe .or. me == 0 ) then
          write(6,*)'Levels',first_wet,' to',last_wet,                  &
     &           ' q_star Two_Norm = ' , Two_Norm
        endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
        Call Two_Norm_levels(                                           &
     &                     R_u, row_length, rows, model_levels,         &
     &                     start_level, end_level,                      &
     &                     model_domain, offx, offy, 1, 0,              &
     &                     at_extremity, nprocx, nprocy,                &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     datastart, .false., .true., 1, Two_Norm )
        if ( L_print_pe .or. me == 0 ) then
          write(6,*)'Levels',start_level,' to',end_level,               &
     &           ' R_u Two_Norm = ' , Two_Norm
        endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
        Call Two_Norm_levels(                                           &
     &                     R_v, row_length, n_rows, model_levels,       &
     &                     start_level, end_level,                      &
     &                     model_domain, offx, offy, 0, 1,              &
     &                     at_extremity, nprocx, nprocy,                &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     datastart, .false., .true., 1, Two_Norm )
        if ( L_print_pe .or. me == 0 ) then
          write(6,*)'Levels',start_level,' to',end_level,               &
     &           ' R_v Two_Norm = ' , Two_Norm
        endif ! L_print_pe .or. me == 0

        If (L_do_w) Then 
! DEPENDS ON: two_norm_levels
        Call Two_Norm_levels(                                           &
     &                     R_w, row_length, rows, model_levels,         &
     &                     start_level, w_end_level,                    &
     &                     model_domain, 0, 0, 0, 0,                    &
     &                     at_extremity, nprocx, nprocy,                &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     datastart, .false., .true., 1, Two_Norm )
        if ( L_print_pe .or. me == 0 ) then
          write(6,*)'Levels',start_level,' to', w_end_level,            &
     &           ' R_w Two_Norm = ' , Two_Norm
        endif ! L_print_pe .or. me == 0
        End If ! L_do_w

      End If !  L_full

      IF (lhook) CALL dr_hook('PRINT_L2NORMS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE print_l2norms

