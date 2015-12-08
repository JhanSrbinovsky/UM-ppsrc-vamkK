! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ******************************COPYRIGHT******************************
!
! subroutine array_l2_print
      SUBROUTINE array_l2_print(                                        &
                                exner_rho_levels, rho,                  &
                                u, v, w,                                &
                                u_adv, v_adv, w_adv,                    &
                                theta, q, qcl, qcf,                     &
                                R_u, R_v, R_w, theta_star,              &
                                q_star, qcl_star, qcf_star,             &
                                row_length, rows, n_rows, rims,         &
                                model_levels, wet_levels,               &
                                start_level, end_level,                 &
                                offx, offy, halo_i, halo_j, me,         &
                                L_do_halos, L_do_rims,                  &
                                L_full, L_print_pe )

! Purpose:
!          To calculate and print l2norms of tendencies or star fields

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Diagnostics

! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

! Arguments with Intent IN. ie: Input variables.
      INTEGER, INTENT(IN) ::                                            &
       row_length,                                                      &
                         ! in: no of points per local row
       rows,                                                            &
                        ! in: no of local (theta) rows
       n_rows,                                                          &
                        ! in: no of local (v) rows
       model_levels,                                                    &
                        ! in: no of model levels
       wet_levels,                                                      &
                        ! in: no of moist-levels
       start_level,                                                     &
                        ! start level for norm calculation
       end_level,                                                       &
                        ! end level for norm calculation
       offx,                                                            &
                        ! standard halo size in east-west
       offy,                                                            &
                        ! standard halo size in North-South
       halo_i,                                                          &
                        ! extended halo size in East-West
       halo_j,                                                          &
                        ! extended halo size in North-South
       me,                                                              &
                        ! this processor
       rims             ! rim size for LAMs

      REAL, INTENT(IN) ::                                               &
       u(1-offx:row_length+offx, 1-offy:rows+offy, model_levels),       &
       v(1-offx:row_length+offx, 1-offy:n_rows+offy, model_levels),     &
       w(1-offx:row_length+offx, 1-offy:rows+offy,                      &
          0:model_levels),                                              &
       u_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
              model_levels),                                            &
       v_adv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,        &
              model_levels),                                            &
       w_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,          &
              0:model_levels),                                          &
       theta(1-offx:row_length+offx, 1-offy:rows+offy, model_levels),   &
       q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,              &
          wet_levels),                                                  &
       qcl(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
            wet_levels),                                                &
       qcf(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,            &
            wet_levels),                                                &
       rho(1-offx:row_length+offx, 1-offy:rows+offy, model_levels),     &
       exner_rho_levels(1-offx:row_length+offx,                         &
                        1-offy:rows+offy, model_levels+1)

      REAL, INTENT(IN) ::                                               &
        R_u(1-offx:row_length+offx, 1-offy:rows+offy,                   &
            model_levels),                                              &
        R_v(1-offx:row_length+offx, 1-offy:n_rows+offy,                 &
            model_levels),                                              &
        R_w(row_length, rows, model_levels - 1),                        &
        theta_star(1-offx:row_length+offx,                              &
                   1-offy:rows+offy, model_levels),                     &
        q_star(1-offx:row_length+offx,                                  &
               1-offy:rows+offy, wet_levels),                           &
        qcl_star(1-offx:row_length+offx,                                &
                 1-offy:rows+offy, wet_levels),                         &
        qcf_star(1-offx:row_length+offx,                                &
                 1-offy:rows+offy, wet_levels)

      LOGICAL                                                           &
       L_do_halos,                                                      &
       L_do_rims,                                                       &
       L_full,                                                          &
                ! T if norms of prognostics, F if norms of increments
       L_print_pe      ! true if  printing on all pe's

!    local variables

      INTEGER :: first_wet
      INTEGER :: last_wet
      INTEGER :: p_end_level
      INTEGER :: w_end_level

      REAL ::  Two_norm

      LOGICAL  :: L_do_w

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('ARRAY_L2_PRINT',zhook_in,zhook_handle)

      IF ( start_level > wet_levels ) THEN
        first_wet = wet_levels
        last_wet = wet_levels
      ELSE IF ( end_level > wet_levels ) THEN
        first_wet = start_level
        last_wet = wet_levels
      ELSE
        first_wet = start_level
        last_wet = end_level
      END IF ! start_level > wet_levels

      p_end_level = end_level
      IF ( end_level == model_levels ) p_end_level = model_levels + 1

      L_do_w = .true.
      w_end_level = end_level
      IF ( end_level == model_levels ) THEN
        w_end_level = model_levels - 1
        IF ( start_level == model_levels ) L_do_w = .false.
      END IF ! end_level == model_levels

      IF (L_full) THEN   !  Norms of prognostic variables

! DEPENDS ON: array_l2norm
        CALL array_l2norm(                                              &
                           exner_rho_levels, row_length, rows,          &
                           model_levels + 1, start_level, p_end_level,  &
                           offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                           rims, rims, Two_Norm )
        IF ( L_print_pe .OR. me == 0 ) THEN
          WRITE(6,*)' L2 Norm of prognostic variables '
          WRITE(6,*)'Levels',start_level,' to', p_end_level ,           &
                    ' Exner pressure Two_Norm = ' , Two_Norm
        END IF ! L_print_pe .OR. me == 0

! DEPENDS ON: array_l2norm
        CALL array_l2norm(                                              &
                           rho, row_length, rows, model_levels,         &
                           start_level, end_level,                      &
                           offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                           rims, rims, Two_Norm )
        IF ( L_print_pe .OR. me == 0 ) THEN
          WRITE(6,*)'Levels',start_level,' to',end_level,               &
                    ' rho Two_Norm = ' , Two_Norm
        END IF ! L_print_pe .OR. me == 0

! DEPENDS ON: array_l2norm
        CALL array_l2norm(                                              &
                           theta, row_length, rows, model_levels,       &
                           start_level, end_level,                      &
                           offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                           rims, rims, Two_Norm )
        IF ( L_print_pe .OR. me == 0 ) THEN
          WRITE(6,*)'Levels',start_level,' to',end_level,               &
                    ' theta Two_Norm = ' , Two_Norm
        END IF ! L_print_pe .OR. me == 0

! DEPENDS ON: array_l2norm
        CALL array_l2norm(                                              &
                           q, row_length, rows, wet_levels,             &
                           first_wet, last_wet,                         &
                           halo_i, halo_j, 0, 0, L_do_halos, L_do_rims, &
                           rims, rims, Two_Norm )
        IF ( L_print_pe .OR. me == 0 ) THEN
          WRITE(6,*)'Levels',first_wet,' to',last_wet,                  &
                    ' q Two_Norm = ' , Two_Norm
        END IF ! L_print_pe .OR. me == 0

! DEPENDS ON: array_l2norm
        CALL array_l2norm(                                              &
                           u, row_length, rows, model_levels,           &
                           start_level, end_level,                      &
                           offx, offy, 1, 0, L_do_halos, L_do_rims,     &
                           rims, rims, Two_Norm )
        IF ( L_print_pe .OR. me == 0 ) THEN
          WRITE(6,*)'Levels',start_level,' to',end_level,               &
                    ' u Two_Norm = ' , Two_Norm
        END IF ! L_print_pe .OR. me == 0

! DEPENDS ON: array_l2norm
        CALL array_l2norm(                                              &
                           v, row_length, n_rows, model_levels,         &
                           start_level, end_level,                      &
                           offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                           rims, rims, Two_Norm )
        IF ( L_print_pe .OR. me == 0 ) THEN
          WRITE(6,*)'Levels',start_level,' to',end_level,               &
                    ' v Two_Norm = ' , Two_Norm
        END IF ! L_print_pe .OR. me == 0

        IF (L_do_w) THEN
! DEPENDS ON: array_l2norm
          CALL array_l2norm(                                            &
                            w, row_length, rows, model_levels + 1,      &
                            start_level, p_end_level,                   &
                            offx, offy, 0, 0, L_do_halos, L_do_rims,    &
                            rims, rims, Two_Norm )
          IF ( L_print_pe .OR. me == 0 ) THEN
            WRITE(6,*)'Levels',start_level,' to', p_end_level,          &
                      ' w Two_Norm = ' , Two_Norm
          END IF ! L_print_pe .OR. me == 0
        END IF ! L_do_w

      ELSE   !  Norms of increments or star fields

! DEPENDS ON: array_l2norm
        CALL array_l2norm(                                              &
                          theta_star, row_length, rows, model_levels,   &
                          start_level, end_level,                       &
                          offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                          rims, rims, Two_Norm )
        IF ( L_print_pe .OR. me == 0 ) THEN
          WRITE(6,*)'Levels',start_level,' to',end_level,               &
                    ' theta_star Two_Norm = ' , Two_Norm
        END IF ! L_print_pe .OR. me == 0

! DEPENDS ON: array_l2norm
        CALL array_l2norm(                                              &
                          q_star, row_length, rows, wet_levels,         &
                          first_wet, last_wet,                          &
                          offx, offy, 0, 0, L_do_halos, L_do_rims,      &
                          rims, rims, Two_Norm )
        IF ( L_print_pe .OR. me == 0 ) THEN
          WRITE(6,*)'Levels',first_wet,' to',last_wet,                  &
                    ' q_star Two_Norm = ' , Two_Norm
        END IF ! L_print_pe .OR. me == 0

! DEPENDS ON: array_l2norm
        CALL array_l2norm(                                              &
                           R_u, row_length, rows, model_levels,         &
                           start_level, end_level,                      &
                           offx, offy, 1, 0, L_do_halos, L_do_rims,     &
                           rims, rims, Two_Norm )
        IF ( L_print_pe .OR. me == 0 ) THEN
          WRITE(6,*)'Levels',start_level,' to',end_level,               &
                    ' R_u Two_Norm = ' , Two_Norm
        END IF ! L_print_pe .OR. me == 0

! DEPENDS ON: array_l2norm
        CALL array_l2norm(                                              &
                           R_v, row_length, n_rows, model_levels,       &
                           start_level, end_level,                      &
                           offx, offy, 0, 0, L_do_halos, L_do_rims,     &
                           rims, rims, Two_Norm )
        IF ( L_print_pe .OR. me == 0 ) THEN
          WRITE(6,*)'Levels',start_level,' to',end_level,               &
                    ' R_v Two_Norm = ' , Two_Norm
        END IF ! L_print_pe .OR. me == 0

        IF (L_do_w) THEN
! DEPENDS ON: array_l2norm
          CALL array_l2norm(                                            &
                            R_w, row_length, rows, model_levels,        &
                            start_level, w_end_level,                   &
                            0, 0, 0, 0, L_do_halos, L_do_rims,          &
                            rims, rims, Two_Norm )
          IF ( L_print_pe .OR. me == 0 ) THEN
            WRITE(6,*)'Levels',start_level,' to', w_end_level,          &
                      ' R_w Two_Norm = ' , Two_Norm
          END IF ! L_print_pe .OR. me == 0
        END IF ! L_do_w

      END IF !  L_full

      IF (lhook) CALL dr_hook('ARRAY_L2_PRINT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE array_l2_print
