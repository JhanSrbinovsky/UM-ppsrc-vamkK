! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to fill in missing data

! Method:

!   Radiative fluxes may not have been calculated at all
!   points: we now fill in as required. This part of the
!   code was originally located in RAD_CTL2 (v6.1 and below)
!   but has been move into a subroutine in order to make
!   RAD_CTL2 more readable.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!-----------------------------------------------------------------------

      SUBROUTINE fill_missing_data_isccp(                               &
        off_x, off_y, row_length, rows,                                 &
        first_row,last_row,                                             &
        first_data_interp, j_lw, es_space_interp,                       &
        l_complete_north, l_complete_south, l_complete_deg)

      USE lw_diag_mod, ONLY: lw_diag
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE


! VARIABLES WITH INTENT IN

      INTEGER                                                           &
        off_x                                                           &
      , off_y                                                           &
      , row_length                                                      &
      , rows                                                            &
      , model_levels

      INTEGER                                                           &
        first_row                                                       &
      , last_row                                                        &
      , first_data_interp                                               &
      , j_lw

      REAL                                                              &
         es_space_interp(4, row_length, rows)

      LOGICAL                                                           &
        l_complete_north                                                &
      , l_complete_south                                                &
      , l_complete_deg                                                  &
      , l_extra_top

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!   Isccp diagnostics

      IF (lhook) CALL dr_hook('FILL_MISSING_DATA_ISCCP',zhook_in,zhook_handle)
        IF ( lw_diag(j_lw)%l_isccp_weights ) THEN
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(                                               &
            l_complete_north, l_complete_south, l_complete_deg,         &
            row_length, rows, off_x, off_y, first_row, last_row,        &
            first_data_interp, es_space_interp, 1,                      &
            lw_diag(j_lw)%isccp_weights                                 &
            )
        END IF

        IF ( lw_diag(j_lw)%l_isccp_cf ) THEN
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(                                               &
            l_complete_north, l_complete_south, l_complete_deg,         &
            row_length, rows, off_x, off_y, first_row, last_row,        &
            first_data_interp, es_space_interp, 7,                      &
            lw_diag(j_lw)%isccp_cf                                      &
            )
        END IF

        IF ( lw_diag(j_lw)%l_isccp_cf_tau_0_to_p3 ) THEN
! DEPENDS ON: rad3d_inp
           CALL rad3d_inp(                                              &
             l_complete_north, l_complete_south, l_complete_deg,        &
             row_length, rows, off_x, off_y, first_row, last_row,       &
             first_data_interp, es_space_interp, 7,                     &
             lw_diag(j_lw)%isccp_cf_tau_0_to_p3                         &
             )
        END IF

        IF ( lw_diag(j_lw)%l_isccp_cf_tau_p3_to_1p3 ) THEN
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(                                               &
             l_complete_north, l_complete_south, l_complete_deg,        &
             row_length, rows, off_x, off_y, first_row, last_row,       &
             first_data_interp, es_space_interp, 7,                     &
             lw_diag(j_lw)%isccp_cf_tau_p3_to_1p3                       &
             )
        END IF

        IF ( lw_diag(j_lw)%l_isccp_cf_tau_1p3_to_3p6 ) THEN
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(                                               &
             l_complete_north, l_complete_south, l_complete_deg,        &
             row_length, rows, off_x, off_y, first_row, last_row,       &
             first_data_interp, es_space_interp, 7,                     &
             lw_diag(j_lw)%isccp_cf_tau_1p3_to_3p6                      &
             )
        END IF

        IF ( lw_diag(j_lw)%l_isccp_cf_tau_3p6_to_9p4 ) THEN
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(                                               &
             l_complete_north, l_complete_south, l_complete_deg,        &
             row_length, rows, off_x, off_y, first_row, last_row,       &
             first_data_interp, es_space_interp, 7,                     &
             lw_diag(j_lw)%isccp_cf_tau_3p6_to_9p4                      &
             )
        END IF

        IF ( lw_diag(j_lw)%l_isccp_cf_tau_9p4_to_23 ) THEN
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(                                               &
            l_complete_north, l_complete_south, l_complete_deg,         &
            row_length, rows, off_x, off_y, first_row, last_row,        &
            first_data_interp, es_space_interp, 7,                      &
            lw_diag(j_lw)%isccp_cf_tau_9p4_to_23                        &
            )
        END IF

        IF ( lw_diag(j_lw)%l_isccp_cf_tau_23_to_60 ) THEN
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(                                               &
             l_complete_north, l_complete_south, l_complete_deg,        &
             row_length, rows, off_x, off_y, first_row, last_row,       &
             first_data_interp, es_space_interp, 7,                     &
             lw_diag(j_lw)%isccp_cf_tau_23_to_60                        &
             )
        END IF

        IF ( lw_diag(j_lw)%l_isccp_cf_tau_ge_60 ) THEN
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(                                               &
            l_complete_north, l_complete_south, l_complete_deg,         &
            row_length, rows, off_x, off_y, first_row, last_row,        &
            first_data_interp, es_space_interp, 7,                      &
            lw_diag(j_lw)%isccp_cf_tau_ge_60                            &
            )
        END IF

        IF ( lw_diag(j_lw)%l_meanalbedocld ) THEN
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(                                               &
            l_complete_north, l_complete_south, l_complete_deg,         &
            row_length, rows, off_x, off_y, first_row, last_row,        &
            first_data_interp, es_space_interp, 1,                      &
            lw_diag(j_lw)%meanalbedocld                                 &
            )
        END IF

        IF ( lw_diag(j_lw)%l_meantaucld ) THEN
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(                                               &
            l_complete_north, l_complete_south, l_complete_deg,         &
            row_length, rows, off_x, off_y, first_row, last_row,        &
            first_data_interp, es_space_interp, 1,                      &
            lw_diag(j_lw)%meantaucld                                    &
            )
        END IF

        IF ( lw_diag(j_lw)%l_meanptop ) THEN
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(                                               &
            l_complete_north, l_complete_south, l_complete_deg,         &
            row_length, rows, off_x, off_y, first_row, last_row,        &
            first_data_interp, es_space_interp, 1,                      &
            lw_diag(j_lw)%meanptop                                      &
            )
        END IF

        IF ( lw_diag(j_lw)%l_totalcldarea ) THEN
! DEPENDS ON: rad3d_inp
          CALL rad3d_inp(                                               &
            l_complete_north, l_complete_south, l_complete_deg,         &
            row_length, rows, off_x, off_y, first_row, last_row,        &
            first_data_interp, es_space_interp, 1,                      &
            lw_diag(j_lw)%totalcldarea                                  &
            )
        END IF

      IF (lhook) CALL dr_hook('FILL_MISSING_DATA_ISCCP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE fill_missing_data_isccp
