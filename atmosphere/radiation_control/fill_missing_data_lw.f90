! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to fill in missing data
!
! Method:
!
!   Radiative fluxes may not have been calculated at all
!   points: we now fill in as required. This part of the
!   code was originally located in RAD_CTL2 (v6.1 and below)
!   but has been move into a subroutine in order to make
!   RAD_CTL2 more readable.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
! Description of Code:
!   FORTRAN 77/90  with extensions listed in documentation.
!
!-----------------------------------------------------------------------
!
      Subroutine fill_missing_data_lw(                                  &
     &  off_x, off_y, row_length, rows, model_levels,                   &
     &  cloud_levels, n_aod_wavel,                                      &
     &  first_row,last_row,                                             &
     &  first_data_interp, ES_space_interp,                             &
     &  L_complete_North, L_complete_South, L_complete_deg,             &
     &  n_channel, j_lw, l_extra_top,                                   &
     &  LW_incs, OLR, lw_down, LWsea, top_absorption )

      Use lw_diag_mod, Only: LW_diag
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

!
! VARIABLES WITH INTENT IN
!
      Integer                                                           &
     &  off_x                                                           &
     &, off_y                                                           &
     &, row_length                                                      &
     &, rows                                                            &
     &, model_levels                                                    &
     &, cloud_levels                                                    &
     &, n_aod_wavel

      Integer                                                           &
     &  first_row                                                       &
     &, last_row                                                        &
     &, first_data_interp                                               &
     &, n_channel                                                       &
     &, j_lw

      Real                                                              &
     &   ES_SPACE_INTERP(4, row_length, rows)

      Logical                                                           &
     &  L_complete_North                                                &
     &, L_complete_South                                                &
     &, L_complete_deg                                                  &
     &, l_extra_top
!
! VARIABLES WITH INTENT IN/OUT
!
      Real                                                              &
     &  LW_incs(row_length, rows, 0:model_levels)                       &
     &, LWsea(row_length, rows)                                         &
     &, OLR(row_length, rows)                                           &
     &, lw_down(row_length, rows)                                       &
     &, top_absorption(row_length, rows)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('FILL_MISSING_DATA_LW',zhook_in,zhook_handle)

!
! Primary Fields:
!
! DEPENDS ON: rad3d_inp
        Call rad3d_inp(                                                 &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp,                         &
     &      model_levels+1,                                             &
     &      LW_incs                                                     &
     &      )
! DEPENDS ON: rad3d_inp
        Call rad3d_inp(                                                 &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 1,                      &
     &      OLR                                                         &
     &      )
! DEPENDS ON: rad3d_inp
        Call rad3d_inp(                                                 &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 1,                      &
     &      lw_down                                                     &
     &      )
! DEPENDS ON: rad3d_inp
        Call rad3d_inp(                                                 &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 1,                      &
     &      LWsea                                                       &
     &      )
        If (l_extra_top) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 1,                      &
     &      top_absorption                                              &
     &      )
        Endif
!
! LW Diagnostics:
!
        IF ( lw_diag(j_lw)%l_flux_up ) THEN
           CALL rad3d_inp(                                              &
              l_complete_north, l_complete_south, l_complete_deg,       &
              row_length, rows, off_x, off_y, first_row, last_row,      &
              first_data_interp, es_space_interp,                       &
              model_levels+1, lw_diag(j_lw)%flux_up )
        END IF
        IF ( lw_diag(j_lw)%l_flux_down ) THEN
           CALL rad3d_inp(                                              &
              l_complete_north, l_complete_south, l_complete_deg,       &
              row_length, rows, off_x, off_y, first_row, last_row,      &
              first_data_interp, es_space_interp,                       &
              model_levels+1, lw_diag(j_lw)%flux_down )
        END IF
        IF ( lw_diag(j_lw)%l_flux_up_clear ) THEN
           CALL rad3d_inp(                                              &
              l_complete_north, l_complete_south, l_complete_deg,       &
              row_length, rows, off_x, off_y, first_row, last_row,      &
              first_data_interp, es_space_interp,                       &
              model_levels+1, lw_diag(j_lw)%flux_up_clear )
        END IF
        IF ( lw_diag(j_lw)%l_flux_down_clear ) THEN
           CALL rad3d_inp(                                              &
              l_complete_north, l_complete_south, l_complete_deg,       &
              row_length, rows, off_x, off_y, first_row, last_row,      &
              first_data_interp, es_space_interp,                       &
              model_levels+1, lw_diag(j_lw)%flux_down_clear )
        END IF
        If ( LW_diag(j_lw)%L_total_cloud_cover ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 1,                      &
     &      LW_diag(j_lw)%total_cloud_cover                             &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_clear_olr ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 1,                      &
     &      LW_diag(j_lw)%clear_olr                                     &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_surf_down_clr ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 1,                      &
     &      LW_diag(j_lw)%surf_down_clr                                 &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_clear_hr ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp,                         &
     &      model_levels,                                               &
     &      LW_diag(j_lw)%clear_hr                                      &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_net_flux_trop ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 1,                      &
     &      LW_diag(j_lw)%net_flux_trop                                 &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_down_flux_trop ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, 1,                      &
     &      LW_diag(j_lw)%down_flux_trop                                &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_total_cloud_on_levels ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, cloud_levels,           &
     &      LW_diag(j_lw)%total_cloud_on_levels                         &
     &      )
        Endif
!       
! Grid-box mean cloud diagnostics as seen by radiation:
!
          If ( LW_diag(j_lw)%L_ls_qcl_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag(j_lw)%ls_qcl_rad                                  &
     &        )
          Endif
          If ( LW_diag(j_lw)%L_ls_qcf_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag(j_lw)%ls_qcf_rad                                  &
     &        )
          Endif
          If ( LW_diag(j_lw)%L_cc_qcl_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag(j_lw)%cc_qcl_rad                                  &
     &        )
          Endif
          If ( LW_diag(j_lw)%L_cc_qcf_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag(j_lw)%cc_qcf_rad                                  &
     &        )
          Endif
          If ( LW_diag(j_lw)%L_ls_cl_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag(j_lw)%ls_cl_rad                                   &
     &        )
          Endif
          If ( LW_diag(j_lw)%L_ls_cf_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag(j_lw)%ls_cf_rad                                   &
     &        )
          Endif
          If ( LW_diag(j_lw)%L_cc_cl_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag(j_lw)%cc_cl_rad                                   &
     &        )
          Endif
          If ( LW_diag(j_lw)%L_cc_cf_rad ) Then
! DEPENDS ON: rad3d_inp
            Call rad3d_inp(                                             &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, model_levels,         &
     &        LW_diag(j_lw)%cc_cf_rad                                   &
     &        )
          Endif 
!
! Radiances
!
        If ( LW_diag(j_lw)%L_toa_radiance ) THEN
! DEPENDS ON: rad3d_inp
           Call rad3d_inp(                                              &
     &        L_complete_North, L_complete_South, L_complete_deg,       &
     &        row_length, rows, off_x, off_y, first_row, last_row,      &
     &        first_data_interp, ES_space_interp, n_channel,            &
     &        LW_diag(j_lw)%toa_radiance                                &
     &       )
        Endif
!
!   Absorptivity diagnostics:
!
        If ( LW_diag(j_lw)%L_cloud_absorptivity ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, cloud_levels,           &
     &      LW_diag(j_lw)%cloud_absorptivity                            &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_cloud_weight_absorptivity ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, cloud_levels,           &
     &      LW_diag(j_lw)%cloud_weight_absorptivity                     &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_ls_cloud_absorptivity ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, cloud_levels,           &
     &      LW_diag(j_lw)%ls_cloud_absorptivity                         &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_ls_cloud_weight_absorptivity ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, cloud_levels,           &
     &      LW_diag(j_lw)%ls_cloud_weight_absorptivity                  &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_cnv_cloud_absorptivity ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, cloud_levels,           &
     &      LW_diag(j_lw)%cnv_cloud_absorptivity                        &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_cnv_cloud_weight_absorptivity ) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, cloud_levels,           &
     &      LW_diag(j_lw)%cnv_cloud_weight_absorptivity                 &
     &      )
        Endif

! Aerosol optical depth diagnostics

        If ( LW_diag(j_lw)%L_aod_sulphate) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_sulphate                                  &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_dust) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_dust                                      &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_seasalt) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_seasalt                                   &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_soot) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_soot                                      &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_biomass) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_biomass                                   &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_biogenic) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_biogenic                                  &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_ocff) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_ocff                                      &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_delta) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_delta                                     &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_nitrate) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_nitrate                                   &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_prog_sulphate) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_prog_sulphate                             &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_prog_dust) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_prog_dust                                 &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_prog_seasalt) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_prog_seasalt                              &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_prog_soot) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_prog_soot                                 &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_prog_biomass) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_prog_biomass                              &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_prog_ocff) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_prog_ocff                                 &
     &      )
        Endif
        If ( LW_diag(j_lw)%L_aod_prog_nitrate) THEN
! DEPENDS ON: rad3d_inp
          Call rad3d_inp(                                               &
     &      L_complete_North, L_complete_South, L_complete_deg,         &
     &      row_length, rows, off_x, off_y, first_row, last_row,        &
     &      first_data_interp, ES_space_interp, n_aod_wavel,            &
     &      LW_diag(j_lw)%aod_prog_nitrate                              &
     &      )
        Endif

! UKCA aerosol optical depth diagnostics 
 
        IF (LW_diag(j_lw)%L_aod_ukca_ait_sol) THEN 
! DEPENDS ON: rad3d_inp 
          CALL rad3d_inp(                                               & 
            L_complete_North, L_complete_South, L_complete_deg,         & 
            row_length, rows, off_x, off_y, first_row, last_row,        & 
            first_data_interp, ES_space_interp, n_aod_wavel,            & 
            LW_diag(j_lw)%aod_ukca_ait_sol                              & 
            ) 
        END IF 
        IF (LW_diag(j_lw)%L_aod_ukca_acc_sol) THEN 
! DEPENDS ON: rad3d_inp 
          CALL rad3d_inp(                                               & 
            L_complete_North, L_complete_South, L_complete_deg,         & 
            row_length, rows, off_x, off_y, first_row, last_row,        & 
            first_data_interp, ES_space_interp, n_aod_wavel,            & 
            LW_diag(j_lw)%aod_ukca_acc_sol                              & 
            ) 
        END IF 
        IF (LW_diag(j_lw)%L_aod_ukca_cor_sol) THEN 
! DEPENDS ON: rad3d_inp 
          CALL rad3d_inp(                                               & 
            L_complete_North, L_complete_South, L_complete_deg,         & 
            row_length, rows, off_x, off_y, first_row, last_row,        & 
            first_data_interp, ES_space_interp, n_aod_wavel,            & 
            LW_diag(j_lw)%aod_ukca_cor_sol                              & 
            ) 
        END IF 
        IF (LW_diag(j_lw)%L_aod_ukca_ait_ins) THEN 
! DEPENDS ON: rad3d_inp 
          CALL rad3d_inp(                                               & 
            L_complete_North, L_complete_South, L_complete_deg,         & 
            row_length, rows, off_x, off_y, first_row, last_row,        & 
            first_data_interp, ES_space_interp, n_aod_wavel,            & 
            LW_diag(j_lw)%aod_ukca_ait_ins                              & 
            ) 
        END IF 
        IF (LW_diag(j_lw)%L_aod_ukca_acc_ins) THEN 
! DEPENDS ON: rad3d_inp 
          CALL rad3d_inp(                                               & 
            L_complete_North, L_complete_South, L_complete_deg,         & 
            row_length, rows, off_x, off_y, first_row, last_row,        & 
            first_data_interp, ES_space_interp, n_aod_wavel,            & 
            LW_diag(j_lw)%aod_ukca_acc_ins                              & 
            ) 
        END IF 
        IF (LW_diag(j_lw)%L_aod_ukca_cor_ins) THEN 
! DEPENDS ON: rad3d_inp 
          CALL rad3d_inp(                                               & 
            L_complete_North, L_complete_South, L_complete_deg,         & 
            row_length, rows, off_x, off_y, first_row, last_row,        & 
            first_data_interp, ES_space_interp, n_aod_wavel,            & 
            LW_diag(j_lw)%aod_ukca_cor_ins                              & 
            ) 
        END IF 

      IF (lhook) CALL dr_hook('FILL_MISSING_DATA_LW',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE fill_missing_data_lw
