! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! subroutine print_L2nhelm
      subroutine print_L2nhelm(                                         &
     &                         HM_Cxx1, HM_Cxx2, HM_Cyy1, HM_Cyy2,      &
     &                         HM_Czz, HM_Cz, HM_Cxz, HM_Cyz,           &
     &                         HM_Cxp, HM_Cyp, HM_C2n,                  &
     &                         HM_Cxy1, HM_Cxy2, HM_Cyx1, HM_Cyx2,      &
     &                         HM_C3, HM_C4, HM_C5, HM_RHS,             &
     &                         row_length, rows, model_levels,          &
     &                         first_constant_r_rho_level_m1,           &
     &                         start_level, end_level,                  &
     &                         off_x, off_y, me, n_procx, n_procy,      &
     &                         proc_row_group, proc_col_group,          &
     &                         global_row_length, global_rows,          &
     &                         at_extremity, l_datastart, model_domain, &
     &                         L_do_halos, L_do_rims, rims_to_do,       &
     &                         L_print_pe )

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
     &, model_levels                                                    &
                          ! in: no of model levels
     &, off_x                                                           &
                           ! standard halo size in east-west
     &, off_y                                                           &
                           ! standard halo size in North-South
     &, rims_to_do                                                      &
     &, model_domain                                                    &
     &, me                                                              &
                          ! processor id
     &, n_procx                                                         &
     &, n_procy                                                         &
     &, proc_row_group                                                  &
     &, proc_col_group                                                  &
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, first_constant_r_rho_level_m1                                   &
                                      ! value used to dimension
                                      ! arrays, max of (1 and
                                      ! first_constant_r_rho_level)
     &, start_level                                                     &
                          ! start level for norm calculation
     &, end_level                                                       &
                          ! end level for norm calculation
     &, l_datastart(3)

      Real, Intent(In) ::                                               &
     & HM_Cxx1(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cxx2(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cxy1(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cxy2(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cyy1(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cyy2(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cyx1(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cyx2(1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Czz (1-off_x:row_length+off_x,1-off_y:rows+off_y,model_levels)&
     &,HM_Cz (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)&
     &,HM_C2n (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         first_constant_r_rho_level_m1)                           &
     &,HM_C3 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)&
     &,HM_C4 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,model_levels)&
     &,HM_C5 (1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cxz (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cyz (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cxp (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_Cyp (1-off_x:row_length+off_x, 1-off_y:rows+off_y,           &
     &         first_constant_r_rho_level_m1)                           &
     &, HM_RHS(row_length, rows, model_levels)

      Logical                                                           &
     &  L_do_halos                                                      &
     &, L_do_rims                                                       &
     &, L_print_pe      ! true if  printing on all pe's

      Logical                                                           &
     &  at_extremity(4)  ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid

!    local variables
      Integer :: first_level
      Integer :: last_level

      Real ::  Two_norm
      Logical ::  L_do

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('PRINT_L2NHELM',zhook_in,zhook_handle)

!    set pointers for arrays with levels < first_constant_r_rho_level
      L_do = .true.
      first_level = start_level
      last_level = end_level
      If ( start_level > first_constant_r_rho_level_m1) Then
        L_do = .false.
      Else If ( end_level > first_constant_r_rho_level_m1) Then
        first_level = start_level
        last_level = first_constant_r_rho_level_m1
      End If ! start_level > first_constant_r_rho_level_m1

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Cxx1, row_length, rows, model_levels,     &
     &                     start_level, end_level,                      &
     &                     model_domain, off_x, off_y, 1, 0,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'   ***   L2_Norms of solver coefficients  ****  '
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_Cxx1 Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Cxx2, row_length, rows, model_levels,     &
     &                     start_level, end_level,                      &
     &                     model_domain, off_x, off_y, 1, 0,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_Cxx2 Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Cyy1, row_length, rows, model_levels,     &
     &                     start_level, end_level,                      &
     &                     model_domain, off_x, off_y, 0, 1,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_Cyy1 Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Cyy2, row_length, rows, model_levels,     &
     &                     start_level, end_level,                      &
     &                     model_domain, off_x, off_y, 0, 1,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_Cyy2 Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Czz, row_length, rows, model_levels,      &
     &                     start_level, end_level,                      &
     &                     model_domain, off_x, off_y, 0, 0,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_Czz Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Cz, row_length, rows, model_levels,       &
     &                     start_level, end_level,                      &
     &                     model_domain, off_x, off_y, 0, 0,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_Cz Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

      If (L_do) Then
! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Cxz, row_length, rows,                    &
     &                     first_constant_r_rho_level_m1,               &
     &                     first_level, last_level,                     &
     &                     model_domain, off_x, off_y, 1, 0,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',first_level,' to',last_level,                &
     &         ' HM_Cxz Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Cyz, row_length, rows,                    &
     &                     first_constant_r_rho_level_m1,               &
     &                     first_level, last_level,                     &
     &                     model_domain, off_x, off_y, 0, 1,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',first_level,' to',last_level,                &
     &         ' HM_Cyz Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Cxp, row_length, rows,                    &
     &                     first_constant_r_rho_level_m1,               &
     &                     first_level, last_level,                     &
     &                     model_domain, off_x, off_y, 1, 0,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',first_level,' to',last_level,                &
     &         ' HM_Cxp Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Cyp, row_length, rows,                    &
     &                     first_constant_r_rho_level_m1,               &
     &                     first_level, last_level,                     &
     &                     model_domain, off_x, off_y, 0, 1,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',first_level,' to',last_level,                &
     &         ' HM_Cyp Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0
      End If !  L_do

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Cxy1, row_length, rows, model_levels,     &
     &                     start_level, end_level,                      &
     &                     model_domain, off_x, off_y, 1, 0,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_Cxy1 Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Cxy2, row_length, rows, model_levels,     &
     &                     start_level, end_level,                      &
     &                     model_domain, off_x, off_y, 0, 1,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_Cxy2 Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Cyx1, row_length, rows, model_levels,     &
     &                     start_level, end_level,                      &
     &                     model_domain, off_x, off_y, 0, 1,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_Cyx1 Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_Cyx2, row_length, rows, model_levels,     &
     &                     start_level, end_level,                      &
     &                     model_domain, off_x, off_y, 1, 0,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_Cyx2 Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

      If (L_do) Then
! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_C2n, row_length, rows,                    &
     &                     first_constant_r_rho_level_m1,               &
     &                     first_level, last_level,                     &
     &                     model_domain, off_x, off_y, 0, 0,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',first_level,' to',last_level,                &
     &         ' HM_C2n Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0
      End If !  L_do

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_C3, row_length, rows, model_levels,       &
     &                     start_level, end_level,                      &
     &                     model_domain, off_x, off_y, 0, 0,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_C3 Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_C4, row_length, rows, model_levels,       &
     &                     start_level, end_level,                      &
     &                     model_domain, off_x, off_y, 0, 0,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_C4 Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

      If (L_do) Then
! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_C5, row_length, rows,                     &
     &                     first_constant_r_rho_level_m1,               &
     &                     first_level, last_level,                     &
     &                     model_domain, off_x, off_y, 0, 0,            &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_C5 Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0
      End If !  L_do

! DEPENDS ON: two_norm_levels
      Call Two_Norm_levels(                                             &
     &                     HM_RHS, row_length, rows, model_levels,      &
     &                     start_level, end_level,                      &
     &                     model_domain, 0, 0, 0, 0,                    &
     &                     at_extremity, n_procx, n_procy,              &
     &                     global_row_length, global_rows,              &
     &                     proc_col_group, proc_row_group,              &
     &                     l_datastart, L_do_halos,                     &
     &                     L_do_rims, rims_to_do, Two_Norm )
      if ( L_print_pe .or. me == 0 ) then
        write(6,*)'Levels',start_level,' to',end_level,                 &
     &         ' HM_RHS Two_Norm = ' , Two_Norm
      endif ! L_print_pe .or. me == 0

      IF (lhook) CALL dr_hook('PRINT_L2NHELM',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE print_L2nhelm

