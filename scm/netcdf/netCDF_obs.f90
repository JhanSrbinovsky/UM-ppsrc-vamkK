! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM to read in netcdf data for forcing

MODULE netcdf_obs

  USE netcdf_arr

  IMPLICIT NONE

!-----------------------------------------------------------------------------
!
! Description:
!   Controls extraction of forcing information from netCDF driver files and
!   passes them back to scm_main
!
! Method:
!   General procedure is to:
!   1) Allocate arrays to hold data read directly from the source netcdf file
!   2) Read in data from netCDF files into arrays created in 1)
!   3) Allocate arrays dimensioned to current SCM resolutions
!   4) Convert/calculate (if required) data in netCDF file to units/fields
!      required by SCM
!   5) Interpolate netCDF data onto SCM grid resolution
!   6) Deallocate dynamic arrays before returning to scm_main
!
!   Each netcdf source will require a module with routines specific
!   to itself.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!-----------------------------------------------------------------------------
! Declarations

  INTEGER, PARAMETER :: &
    TWPICE = 1          &
  , RACMO  = 2

CONTAINS
!==============================================================================

  SUBROUTINE get_netcdf_obs

    USE scm_utils, ONLY:                                                      &
        zhook_in, zhook_out, jprb, lhook, dr_hook

    USE s_main_force, ONLY:                                                   &
        year_init, month_init, day_init, hour_init, min_init, sec_init        &
      , obs_pd, obs_top, obs_bot,  smi_opt, l_vertadv, albsoil, smci          &
      , canopy_gbi, lat, long, gridbox_area, orog, smcli, t_deep_soili        &
      , snodepi, tstari, z0mseai, z0m_scm, z0h_scm, frac_typ, canopy          &
      , z0_tile, tstar_tile, ozone, theta, qi, p_in, ui, vi, wi               &
      , tstar_forcing, flux_h, flux_e, t_inc, q_star, u_inc, v_inc, w_inc     &
      , t_bg, q_bg, u_bg, v_bg, w_bg, land_sea_mask, soil_type, source        &
      , netcdf_file

    USE twpice_netcdf
    USE racmo_netcdf


    IMPLICIT NONE


    !----------------------------------------------
    ! Local variables
    !----------------------------------------------
    INTEGER        :: k, cnt ! Counters
    CHARACTER(LEN=20)  :: nc_list(100) ! Register which fields in netcdf
                                   ! data have been used.

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('GET_NETCDF_OBS',zhook_in,zhook_handle)

    SELECT CASE (source)

    CASE (TWPICE)
      CALL alloc_TWPICE_nc
      CALL read_TWPICE_nc
      CALL alloc_nc_arr
      CALL TWPICE_nc_var_to_scm
      CALL TWPICE_nc_grid_to_scm
      CALL dealloc_TWPICE_nc

    CASE (RACMO)
      CALL alloc_RACMO_nc
      CALL alloc_nc_arr
      CALL read_RACMO_nc
      CALL RACMO_nc_var_to_scm
      CALL RACMO_nc_grid_to_scm
      CALL dealloc_RACMO_nc

      PRINT*, &
      ' Interactive vertical advection enabled by RACMO'
      l_vertadv = .TRUE.
      ! Currently no vertical adv T/q incs in driver files

    CASE default
      PRINT*, 'Unknown NetCDF source'
    END SELECT

    PRINT *,                                            &
         ' Data sourced from NetCDF File:'              &
         , ' '//TRIM(netcdf_file)                       &
         , '========================================'// &
         '======================================='

    cnt = 1

    IF (yeari_nc /= IMDI) THEN
      year_init = yeari_nc
      nc_list(cnt) = 'year_init'
      cnt = cnt + 1
    END IF

    IF (monthi_nc /= IMDI) THEN
      month_init = monthi_nc
      nc_list(cnt) = 'month_init'
      cnt = cnt + 1
    END IF

    IF (dayi_nc /= IMDI) THEN
      day_init = dayi_nc
      nc_list(cnt) = 'day_init'
      cnt = cnt + 1
    END IF

    IF (houri_nc /= IMDI) THEN
      hour_init  = houri_nc
      nc_list(cnt) = 'hour_init'
      cnt = cnt + 1
    END IF

    IF (mini_nc /= IMDI)  THEN
      min_init   = mini_nc
      nc_list(cnt) = 'min_init'
      cnt = cnt + 1
    END IF

    IF (seci_nc /= IMDI)  THEN
      sec_init = seci_nc
      nc_list(cnt) = 'sec_init'
      cnt = cnt + 1
    END IF

    IF (obs_pd_nc /= RMDI)  THEN
      obs_pd = obs_pd_nc
      nc_list(cnt) = 'obs_pd'
      cnt = cnt + 1
!      nc_list(cnt) = 'ichgf'
!      cnt = cnt + 1
    END IF

    IF (obs_top_nc /= RMDI)  THEN
      obs_top = obs_top_nc
      nc_list(cnt) = 'obs_top'
      cnt = cnt + 1
    END IF

    IF (obs_bot_nc /= RMDI)  THEN
      obs_bot = obs_bot_nc
      nc_list(cnt) = 'obs_bot'
      cnt = cnt + 1
    END IF

    IF (lat_nc /= RMDI) THEN
      lat = lat_nc
      nc_list(cnt) = 'lat'
      cnt = cnt + 1
    END IF

    IF (long_nc /= RMDI) THEN
      long = long_nc
      nc_list(cnt) = 'long'
      cnt = cnt + 1
    END IF

    IF (tstari_nc /= RMDI) THEN
      tstari = tstari_nc
      nc_list(cnt) = 'tstari'
      cnt = cnt + 1
    END IF


    IF (lsm_nc /= RMDI) THEN
      IF (lsm_nc > 0.5) THEN
        land_sea_mask = .TRUE.

        nc_list(cnt) = 'land_sea_mask'
        cnt = cnt + 1

        IF (l_smci) THEN
          WHERE (smci_nc /= RMDI) smci = smci_nc
          nc_list(cnt) = 'smci'
          cnt = cnt + 1
        END IF


        IF (canopy_gbi_nc(1) /= RMDI) canopy_gbi = canopy_gbi_nc
        nc_list(cnt) = 'canopy_gbi'
        cnt = cnt + 1


        IF (l_smcli) THEN
          WHERE (smcli_nc /= RMDI) smcli = smcli_nc
          smi_opt = 0
          nc_list(cnt) = 'smi_opt'
          cnt = cnt + 1
          nc_list(cnt) = 'smcli'
          cnt = cnt + 1
        END IF

        IF (l_t_deep_soili) THEN
          WHERE (t_deep_soili_nc /= RMDI) t_deep_soili = t_deep_soili_nc
          nc_list(cnt) = 't_deep_soili'
          cnt = cnt + 1
        END IF

        IF (l_frac_typ) THEN
          WHERE (frac_typ_nc /= RMDI) frac_typ = frac_typ_nc
          nc_list(cnt) = 'frac_typ'
          cnt = cnt + 1
        END IF

        IF (l_canopy) THEN
          WHERE (canopy_nc /= RMDI) canopy = canopy_nc
          nc_list(cnt) = 'canopy'
          cnt = cnt + 1
        END IF

        IF (l_z0_tile) THEN
          WHERE (z0_tile_nc /= RMDI) z0_tile = z0_tile_nc
          nc_list(cnt) = 'z0_tile'
          cnt = cnt + 1
        END IF

        IF (l_tstar_tile) THEN
          WHERE (tstar_tile_nc /= RMDI) tstar_tile = tstar_tile_nc
          nc_list(cnt) = 'tstar_tile'
          cnt = cnt + 1
        END IF

        IF (sndep_nc /= RMDI) THEN
          snodepi = sndep_nc
          nc_list(cnt) = 'snodepi'
          cnt = cnt + 1
        END IF

        IF (orog_nc /= RMDI) THEN
          orog = orog_nc
          nc_list(cnt) = 'orog'
          cnt = cnt + 1
        END IF

        IF (albsoil_nc /= RMDI) THEN
          albsoil = albsoil_nc
          nc_list(cnt) = 'albsoil'
          cnt = cnt + 1
        END IF

        IF (soil_type_nc /= IMDI) THEN
          soil_type = soil_type_nc
          nc_list(cnt) = 'soil_type'
          cnt = cnt + 1
        END IF

      ELSE
        land_sea_mask = .FALSE.
        nc_list(cnt) = 'land_sea_mask'
        cnt = cnt + 1

        IF (z0m_nc /= RMDI) THEN
          z0m_scm = z0m_nc
          z0mseai = z0m_nc
          nc_list(cnt) = 'z0m_scm'
          cnt = cnt + 1
          nc_list(cnt) = 'z0mseai'
          cnt = cnt + 1
        END IF

        IF (z0h_nc /= RMDI)  THEN
          z0h_scm = z0h_nc
          nc_list(cnt) = 'z0h_scm'
          cnt = cnt + 1
        END IF

      END IF
    END IF


    IF (grdbx_area_nc /= RMDI) THEN
      gridbox_area = grdbx_area_nc
      nc_list(cnt) = 'gridbox_area'
      cnt = cnt + 1
    END IF

    IF (l_ozone) THEN
      WHERE (o3_nc /= RMDI) ozone = o3_nc
      nc_list(cnt) = 'ozone'
      cnt = cnt + 1
    END IF

    IF (l_theta) THEN
      WHERE (thi_nc /= RMDI) theta = thi_nc
      nc_list(cnt) = 'theta'
      cnt = cnt + 1
    END IF

    IF (l_qi) THEN
      WHERE (qi_nc /= RMDI) qi = qi_nc
      nc_list(cnt) = 'qi'
      cnt = cnt + 1
    END IF

    IF (l_pi) THEN
      WHERE (pi_rh_nc /= RMDI) p_in = pi_rh_nc
      nc_list(cnt) = 'p_in'
      cnt = cnt + 1
    END IF

    IF (l_ui) THEN
      WHERE (ui_nc /= RMDI) ui = ui_nc
      DO k=2, nmod_lv
        WHERE ((ui(:,:,k) == rmdi).AND.(ui(:,:,k-1) /= rmdi))                &
          ui(:,:,k) = ui(:,:,k-1)
      END DO
      nc_list(cnt) = 'ui'
      cnt = cnt + 1
    END IF

    IF (l_vi) THEN
      WHERE (vi_nc /= RMDI) vi = vi_nc
      DO k=2, nmod_lv
        WHERE ((vi(:,:,k) == rmdi).AND.(vi(:,:,k-1) /= rmdi))                &
          vi(:,:,k) = vi(:,:,k-1)
      END DO
      nc_list(cnt) = 'vi'
      cnt = cnt + 1
    END IF

    IF (l_wi) THEN
      WHERE (wi_nc /= RMDI) wi = wi_nc
      nc_list(cnt) = 'wi'
      cnt = cnt + 1
    END IF

    IF (nc_obsfor) THEN

      IF (l_lhflx) THEN
        WHERE (lhflx_nc /= RMDI) flux_h = lhflx_nc
        nc_list(cnt) = 'flux_h'
        cnt = cnt + 1
      END IF

      IF (l_shflx) THEN
        WHERE (shflx_nc /= RMDI) flux_e = shflx_nc
        nc_list(cnt) = 'flux_e'
        cnt = cnt + 1
      END IF

      IF (l_tstar_for) THEN
        WHERE (tstar_nc /= RMDI) tstar_forcing = tstar_nc
        nc_list(cnt) = 'tstar_forcing'
        cnt = cnt + 1
      END IF

      IF (l_u_inc) THEN
        WHERE (u_h_inc_nc /= RMDI) u_inc = u_h_inc_nc
        IF (.NOT. l_vertadv) THEN
          WHERE (u_v_inc_nc /= RMDI) u_inc = u_inc  + u_v_inc_nc
        END IF
        DO k=2, nmod_lv
          WHERE ((u_inc(:,:,:,k) == rmdi).AND.(u_inc(:,:,:,k-1) /= rmdi))    &
            u_inc(:,:,:,k) = u_inc(:,:,:,k-1)
        END DO
        nc_list(cnt) = 'u_inc'
        cnt = cnt + 1
      END IF

      IF (l_v_inc) THEN
        WHERE (v_h_inc_nc /= RMDI) v_inc = v_h_inc_nc
        IF (.NOT. l_vertadv) THEN
          WHERE (v_v_inc_nc /= RMDI) v_inc = v_inc  + v_v_inc_nc
        END IF
        DO k=2, nmod_lv
          WHERE ((v_inc(:,:,:,k) == rmdi).AND.(v_inc(:,:,:,k-1) /= rmdi))    &
            v_inc(:,:,:,k) = v_inc(:,:,:,k-1)
        END DO
        nc_list(cnt) = 'v_inc'
        cnt = cnt + 1
      END IF

      IF (l_w_inc) THEN
        WHERE (w_inc_nc /= RMDI) w_inc = w_inc_nc
        nc_list(cnt) = 'w_inc'
        cnt = cnt + 1
      END IF

      IF (l_t_inc) THEN
        WHERE (t_h_inc_nc /= RMDI) t_inc = t_h_inc_nc
        IF (.NOT. l_vertadv) THEN
          WHERE (t_v_inc_nc /= RMDI) t_inc = t_inc  + t_v_inc_nc
        END IF
        nc_list(cnt) = 't_inc'
        cnt = cnt + 1
      END IF

      IF (l_q_inc) THEN
        WHERE (q_h_inc_nc /= RMDI) q_star = q_h_inc_nc
        IF (.NOT. l_vertadv) THEN
          WHERE (q_v_inc_nc /= RMDI) q_star = q_star + q_v_inc_nc
        END IF
        nc_list(cnt) = 'q_star'
        cnt = cnt + 1
      END IF

      IF (l_t_bg) THEN
        WHERE (t_bg_nc /= RMDI) t_bg = t_bg_nc
        nc_list(cnt) = 't_bg'
        cnt = cnt + 1
      END IF

      IF (l_q_bg) THEN
        WHERE (q_bg_nc /= RMDI) q_bg = q_bg_nc
        nc_list(cnt) = 'q_bg'
        cnt = cnt + 1
      END IF

      IF (l_u_bg) THEN
        WHERE (u_bg_nc /= RMDI) u_bg = u_bg_nc
        DO k=2, nmod_lv
          WHERE ((u_bg(:,:,:,k) == rmdi).AND.(u_bg(:,:,:,k-1) /= rmdi))    &
            u_bg(:,:,:,k) = u_bg(:,:,:,k-1)
        END DO
        nc_list(cnt) = 'u_bg'
        cnt = cnt + 1
      END IF

      IF (l_v_bg) THEN
        WHERE (v_bg_nc /= RMDI) v_bg = v_bg_nc
        DO k=2, nmod_lv
          WHERE ((v_bg(:,:,:,k) == rmdi).AND.(v_bg(:,:,:,k-1) /= rmdi))    &
            v_bg(:,:,:,k) = v_bg(:,:,:,k-1)
        END DO
        nc_list(cnt) = 'v_bg'
        cnt = cnt + 1
      END IF

      IF (l_w_bg) THEN
        WHERE (w_bg_nc /= RMDI) w_bg = w_bg_nc
        nc_list(cnt) = 'w_bg'
        cnt = cnt + 1
      END IF

    END IF

    IF (cnt > 1) THEN
      WRITE(6,'(("  ",A15,3("  ",A15)))') nc_list(1:cnt-1)
      PRINT *, ' '
      PRINT *, 'Note: The above variables will be overwritten'
      PRINT *, '      if specified in namelist.scm'
    ELSE
      WRITE (6,'(A50)')  'No data used'
    END IF

    CALL dealloc_nc_arr

    IF (lhook) CALL dr_hook('GET_NETCDF_OBS',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE get_netcdf_obs

!==============================================================================
END MODULE netcdf_obs
