! *********************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *********************************COPYRIGHT*********************************
!
! Contains routines to read in netcdf driver format created in RACMO format
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: SCM

MODULE racmo_netcdf

  USE netcdf
  USE s_main_force,  ONLY: netcdf_file
  USE scm_utils
  USE soil_param,    ONLY: dzsoil_jules => dzsoil
  USE s_interp_mod,  ONLY: interp1d
  USE id_arr_netcdf, ONLY:                                                    &
      status,     ncid,       dimid_time, dimid_nlev, dimid_nlevp1            &
    , dimid_slev, nt_nc,      nk_nc,      nk_half_nc, ns_nc                   &
    , xtype,      coor_lena,  coor_lenb,  attnum                              &

    ! Variable ids
    !
    , id_date,    id_time,    id_sec,     id_lat,     id_lon                  &
    , id_ps,      id_u,       id_v,       id_t,       id_q                    &
    , id_tadv,    id_qadv,    id_omega,   id_z0m,     id_z0h                  &
    , id_tsskin,  id_qskin,   id_flx_H,   id_flx_E,   id_orog                 &
    , id_alb,     id_lsm,     id_sndep,   id_s_thk,   id_sm                   &
    , id_st,      id_hvc,     id_lvc,     id_ug,      id_vg                   &
    , id_uadv,    id_vadv                                                     &

    ! Variables not currently used, but may be utilised in future
    !============================================================
    ! id_OSST,    id_sice_t,  id_sice_cf, id_rho_snw, id_alb_snw, id_tsnow,
    ! id_uadv,    id_vadv,    id_cf,      id_ql,      id_qi

    ! Receiving arrays
    , file_coor_par_a, file_coor_par_b, file_z_half_prf, file_z_full_prf      &
    , file_p_half,     file_p_full,     file_flx_h,      file_flx_e           &
    , file_t_skin,     file_q_skin,     file_p_srf,      file_q               &
    , file_qadv,       file_t,          file_tadv,       file_u               &
    , file_v,          file_w,          file_omega,      file_exner_full      &
    , file_dzsoil,     file_st,         file_sm,         file_uadv            &
    , file_vadv,       file_date,       file_secs,       file_time            &
    , file_lat,        file_long,       file_lsm,        file_orog            &
    , file_z0h,        file_alb,        file_z_snow,     file_hvc             &
    , file_lvc,        file_z0m

  USE netcdf_arr
  USE global_scmop, ONLY: incdf

  USE atmos_constants_mod,  ONLY: r, kappa, p_zero, recip_kappa
  USE earth_constants_mod,  ONLY: g

  IMPLICIT NONE


  REAL, PRIVATE, ALLOCATABLE :: rdum2d(:,:)

CONTAINS
!==============================================================================

  SUBROUTINE alloc_racmo_nc

    IMPLICIT NONE

    ! Dr Hook
    !=============================================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('ALLOC_RACMO_NC',zhook_in,zhook_handle)

    status = NF90_OPEN(TRIM(ADJUSTL(netcdf_file)), nf90_nowrite, ncid)

    IF (status /= nf90_noerr) THEN
      PRINT*,'Error opening netcdf file: '//TRIM(ADJUSTL(netcdf_file))
    END IF

    ! Get NetCDF file dimension ID numbers
    status = NF90_INQ_DIMID (ncid, 'time',   dimid_time)
    status = NF90_INQ_DIMID (ncid, 'nlev',   dimid_nlev)
    status = NF90_INQ_DIMID (ncid, 'nlevp1', dimid_nlevp1)
    status = NF90_INQ_DIMID (ncid, 'nlevs',  dimid_slev)


    ! Get dimension sizes
    status = NF90_INQUIRE_DIMENSION (ncid, dimid_time,   len=nt_nc)
    status = NF90_INQUIRE_DIMENSION (ncid, dimid_nlev,   len=nk_nc)
    status = NF90_INQUIRE_DIMENSION (ncid, dimid_nlevp1, len=nk_half_nc)
    status = NF90_INQUIRE_DIMENSION (ncid, dimid_slev,   len=ns_nc)

    status = NF90_INQUIRE_ATTRIBUTE (ncid, nf90_global, 'coor_par_a',&
                                     xtype, coor_lena, attnum)
    status = NF90_INQUIRE_ATTRIBUTE (ncid, nf90_global, 'coor_par_a',&
                                     xtype, coor_lenb, attnum)


    ! Allocate dynamic arrays
    !----------------------------------------------
    ! Rank 1 arrays
    ALLOCATE                         &
      ( file_coor_par_a (coor_lena)  &
      , file_coor_par_b (coor_lenb)  &
      , file_z_half_prf (nk_half_nc) &
      , file_z_full_prf (nk_nc)      &
      , file_dzsoil     (ns_nc)      &! Soil layer thickness (m)
      , file_date       (nt_nc)      &
      , file_secs       (nt_nc)      &
      , file_time       (nt_nc)      &
      , file_lat        (nt_nc)      &
      , file_long       (nt_nc)      &
      , file_lsm        (nt_nc)      &
      , file_orog       (nt_nc)      &
      , file_p_srf      (nt_nc)      &! Surface pressure (Pa)
      , file_t_skin     (nt_nc)      &! Skin temperature (K)
      , file_q_skin     (nt_nc)      &! Skin reservoir content (m)
      , file_flx_h      (nt_nc)      &! Srf. sensible heat flux (W/m2)
      , file_flx_e      (nt_nc)      &! Srf. latent heat flux (W/m2)
      , file_z0h        (nt_nc)      &
      , file_z0m        (nt_nc)      &
      , file_alb        (nt_nc)      &
      , file_z_snow     (nt_nc)      &
      , file_hvc        (nt_nc)      &
      , file_lvc        (nt_nc) )


    ! Rank 2 arrays
    ALLOCATE                               &
      ( file_p_half     (nt_nc,nk_half_nc) &!
      , file_p_full     (nt_nc,nk_nc)      &!
      , file_st         (nt_nc,ns_nc)      &! Soil temperture profile (K)
      , file_sm         (nt_nc,ns_nc)      &! Soil moisture profile
      , file_u          (nt_nc,nk_nc)      &! Zonal wind (m/s)
      , file_v          (nt_nc,nk_nc)      &! Meridional wind (m/s)
      , file_w          (nt_nc,nk_nc)      &! vertical wind (m/s)
      , file_t          (nt_nc,nk_nc)      &! Temperture (K)
      , file_q          (nt_nc,nk_nc)      &! H2O mixing ratio (vapour)(kg/kg)
      , file_tadv       (nt_nc,nk_nc)      &! Adv. T tendency (K/s)
      , file_qadv       (nt_nc,nk_nc)      &! Adv. q tendency (kg/kg/s)
      , file_uadv       (nt_nc,nk_nc)      &! Adv. u tendency (m/s)/s
      , file_vadv       (nt_nc,nk_nc)      &! Adv. u tendency (m/s)/s
      , file_omega      (nt_nc,nk_nc)      &! Vert. pressure velocity (Pa/s)
      , file_exner_full (nt_nc,nk_nc) )

    status = NF90_CLOSE(ncid)

    ! Initialise Arrays
    file_p_full     = rmdi
    file_p_half     = rmdi
    file_z_full_prf = rmdi
    file_z_half_prf = rmdi
    file_exner_full = rmdi

    file_date   = imdi
    file_secs   = imdi
    file_time   = imdi

    file_dzsoil = rmdi
    file_lat    = rmdi
    file_long   = rmdi
    file_lsm    = rmdi
    file_orog   = rmdi
    file_p_srf  = rmdi
    file_t_skin = rmdi
    file_q_skin = rmdi
    file_flx_h  = rmdi
    file_flx_e  = rmdi
    file_z0h    = rmdi
    file_z0m    = rmdi
    file_alb    = rmdi
    file_z_snow = rmdi
    file_hvc    = rmdi
    file_lvc    = rmdi
    file_u      = rmdi
    file_v      = rmdi
    file_w      = rmdi
    file_t      = rmdi
    file_q      = rmdi
    file_tadv   = rmdi
    file_qadv   = rmdi
    file_omega  = rmdi
    file_sm     = rmdi
    file_st     = rmdi
    file_uadv   = rmdi
    file_vadv   = rmdi

    IF (lhook) CALL dr_hook('ALLOC_RACMO_NC',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE alloc_racmo_nc

!==============================================================================

  SUBROUTINE dealloc_racmo_nc

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('DEALLOC_RACMO_NC',zhook_in,zhook_handle)

    DEALLOCATE                                                               &
      ( file_coor_par_a, file_coor_par_b, file_z_half_prf, file_z_full_prf   &
      , file_date,       file_secs,       file_time,       file_lat          &
      , file_long,       file_lsm,        file_dzsoil                        &
      , file_orog,       file_p_srf,      file_t_skin,     file_q_skin       &
      , file_flx_h ,     file_flx_e,      file_z0h,        file_z0m          &
      , file_alb,        file_z_snow,     file_hvc,        file_lvc          &
      , file_p_half,     file_p_full,     file_st,         file_sm           &
      , file_u,          file_v,          file_w,          file_t            &
      , file_q,          file_tadv,       file_qadv,       file_uadv         &
      , file_vadv,       file_omega,      file_exner_full )

    IF (lhook) CALL dr_hook('DEALLOC_RACMO_NC',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE dealloc_racmo_nc

!=============================================================================

   SUBROUTINE read_racmo_nc

    IMPLICIT NONE

    REAL :: zsoil      (sm_lv+1)
    REAL :: file_zsoil (sm_lv+1)
    REAL :: wat2dep    (sm_lv+1)

    INTEGER :: i,k,kk

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('READ_RACMO_NC',zhook_in,zhook_handle)

    ! Open netcdf file
    status = NF90_OPEN (TRIM(ADJUSTL(netcdf_file)), nf90_noWrite, ncid)

    IF (status /= nf90_noerr) THEN
      PRINT*,'Error opening netcdf file: '//TRIM(ADJUSTL(netcdf_file))
    END IF

    ! Get id of data to extract
    ! NOTE: These names are specific to RACMO netcdf file
    status = NF90_INQ_VARID (ncid, 'date',           id_date)
    status = NF90_INQ_VARID (ncid, 'time',           id_time)
    status = NF90_INQ_VARID (ncid, 'second',         id_sec)
    status = NF90_INQ_VARID (ncid, 'lat',            id_lat)
    status = NF90_INQ_VARID (ncid, 'lon',            id_lon)
    status = NF90_INQ_VARID (ncid, 'ps',             id_ps)
    status = NF90_INQ_VARID (ncid, 'u',              id_u)
    status = NF90_INQ_VARID (ncid, 'v',              id_v)
    status = NF90_INQ_VARID (ncid, 't',              id_t)
    status = NF90_INQ_VARID (ncid, 'q',              id_q)
    status = NF90_INQ_VARID (ncid, 'tadv',           id_tadv)
    status = NF90_INQ_VARID (ncid, 'qadv',           id_qadv)
    status = NF90_INQ_VARID (ncid, 'uadv',           id_uadv)
    status = NF90_INQ_VARID (ncid, 'vadv',           id_vadv)
    status = NF90_INQ_VARID (ncid, 'omega',          id_omega)
    status = NF90_INQ_VARID (ncid, 'mom_rough',      id_z0m)
    status = NF90_INQ_VARID (ncid, 'heat_rough',     id_z0h)
    status = NF90_INQ_VARID (ncid, 't_skin',         id_tsskin)
    status = NF90_INQ_VARID (ncid, 'q_skin',         id_qskin)
    status = NF90_INQ_VARID (ncid, 'sfc_sens_flx',   id_flx_H)
    status = NF90_INQ_VARID (ncid, 'sfc_lat_flx',    id_flx_E)
    status = NF90_INQ_VARID (ncid, 'orog',           id_orog)
    status = NF90_INQ_VARID (ncid, 'albedo',         id_alb)
    status = NF90_INQ_VARID (ncid, 'lsm',            id_lsm)
    status = NF90_INQ_VARID (ncid, 'snow',           id_sndep)
    status = NF90_INQ_VARID (ncid, 'q_soil',         id_sm)
    status = NF90_INQ_VARID (ncid, 't_soil',         id_st)
    status = NF90_INQ_VARID (ncid, 'ug',             id_ug)
    status = NF90_INQ_VARID (ncid, 'vg',             id_vg)
    status = NF90_INQ_VARID (ncid, 'h_soil',         id_s_thk)
    status = NF90_INQ_VARID (ncid, 'high_veg_cover', id_hvc)
    status = NF90_INQ_VARID (ncid, 'low_veg_cover',  id_lvc)

    ! Variables not currently used, but may be utilised in future
    !============================================================
    ! status = NF90_INQ_VARID (ncid, 'open_sst',       id_OSST)
    ! status = NF90_INQ_VARID (ncid, 't_sea_ice',      id_sice_t)
    ! status = NF90_INQ_VARID (ncid, 'sea_ice_frct',   id_sice_cf)
    ! status = NF90_INQ_VARID (ncid, 'density_snow',   id_rho_snw)
    ! status = NF90_INQ_VARID (ncid, 'albedo_snow',    id_alb_snw)
    ! status = NF90_INQ_VARID (ncid, 't_snow',         id_tsnow)
    ! status = NF90_INQ_VARID (ncid, 'cloud_fraction', id_cf)
    ! status = NF90_INQ_VARID (ncid, 'ql',             id_ql)
    ! status = NF90_INQ_VARID (ncid, 'qi',             id_qi)


    !-------------------------------------------------------------------------
    ! Get obs data from netcdf file
    !-------------------------------------------------------------------------

    ! Extract start date/time at time obs_t0
    !----------------------------------------------
    ! RACMO stores date as a integer in form yyyymmdd
    status    = NF90_GET_VAR (ncid, id_date, file_date)
    yeari_nc  = INT(file_date(obs_t0)*0.0001)
    monthi_nc = MOD(INT(file_date(obs_t0)*0.01),100)
    dayi_nc   = MOD(file_date(obs_t0), 100)


    ! Extract start time
    ! RACMO stores time as seconds into day
    status    = NF90_GET_VAR (ncid, id_sec, file_secs)
    houri_nc  = INT(file_secs(obs_t0)/3600)
    mini_nc   = INT((file_secs(obs_t0) - houri_nc*3600)/60)
    seci_nc   = file_secs(obs_t0) - houri_nc*3600 - mini_nc*60


    ! Calculate observation period
    status    = NF90_GET_VAR (ncid, id_time, file_time)
    obs_pd_nc = file_time(obs_t0+1) - file_time(obs_t0)


    ! Get surface timeseries of:
    !----------------------------------------------
    ! Sensible heat flux (W/m2)
    status = NF90_GET_VAR (ncid, id_flx_h, file_flx_h)

    ! Latent heat flux (W/m2)
    status = NF90_GET_VAR (ncid, id_flx_e, file_flx_e)

    ! Skin temperature (K)
    status = NF90_GET_VAR (ncid, id_tsskin, file_t_skin)

    ! Mean Pressure (mb)
    status = NF90_GET_VAR (ncid, id_ps, file_p_srf)


    ! Get set-up variables at time obs_t0
    !----------------------------------------------

    ! Latitude
    status   = NF90_GET_VAR (ncid, id_lat, file_lat)
    lat_nc   = file_lat(obs_t0)

    ! Longitude
    status   = NF90_GET_VAR (ncid, id_lon, file_long)
    long_nc  = file_long(obs_t0)

    ! Surface temperature
    tstari_nc = file_t_skin(obs_t0)

    ! Roughness length (heat)
    status   = NF90_GET_VAR (ncid, id_z0h, file_z0h)
    z0h_nc   = file_z0h(obs_t0)

    ! Roughness length (momentum)
    status   = NF90_GET_VAR (ncid, id_z0h, file_z0m)
    z0m_nc   = file_z0m(obs_t0)

    ! Land mask/fraction
    status   = NF90_GET_VAR (ncid, id_lsm, file_lsm)
    lsm_nc   = file_lsm(obs_t0)


    ! Class land fractions > 0.5 as land-points.
    !----------------------------------------------
    ! In future, may use this fraction to be able to run coastal points
    IF (lsm_nc > 0.5) THEN

      ! Soil albedo
      status     = NF90_GET_VAR (ncid, id_alb, file_alb)
      albsoil_nc = file_alb(obs_t0)

      ! Canopy water content
      ! Convert skin reservoir content from metres to mm of water
      status           = NF90_GET_VAR (ncid, id_qskin, file_q_skin)
      canopy_nc(:,:)   = file_q_skin(obs_t0)*1000.0
      canopy_gbi_nc(:) = file_q_skin(obs_t0)*1000.0


      ! Orographic height
      status   = NF90_GET_VAR (ncid, id_qskin, file_orog)
      orog_nc  = file_orog(obs_t0)

      ! Snow depth
      status   = NF90_GET_VAR (ncid, id_sndep, file_z_snow)
      sndep_nc = file_z_snow(obs_t0)
      sndep_nc = MAX(sndep_nc,0.0)

      ! Tile roughness lengths
      z0_tile_nc(:,:)    = z0m_nc

      ! Tile thermal roughness lengths
      z0h_tile_nc(:,:)    = z0h_nc

      ! Tile surface temperature
      tstar_tile_nc(:,:) = tstari_nc


      ! Fractional type of land tiles
      !--------------------------------------------
      ! Spread high vegetation cover across tiles 1:2 which are trees
      status             = NF90_GET_VAR (ncid, id_hvc, file_hvc)
      frac_typ_nc(:,1:2) = file_hvc(obs_t0)/2.0

      ! Spread low vegetation cover across tiles 3:5 which are grasses/shrubs
      status  = nf90_get_var (ncid, id_lvc, file_lvc)
      frac_typ_nc(:,3:5) = file_lvc(obs_t0)/3.0

      ! Set ice tile fraction to zero
      frac_typ_nc(:,9)   = 0.0

      ! Spread remaining fraction over other tiles
      frac_typ_nc(:,6:8) = (1.0 - SUM(frac_typ_nc(:,1:5)))/3.0

      ! Top up trees with any discrepency so sum of fractional tiles equal 1
      frac_typ_nc(:,1) = 1.0 - SUM(frac_typ_nc(:,2:8))

      ! Set logicals to register variables obtained from netcdf file
      l_frac_typ   = .TRUE.
      l_canopy     = .TRUE.
      l_z0_tile    = .TRUE.
      l_tstar_tile = .TRUE.

    END IF


    ! Soil moisture/temperature
    IF (lsm_nc > 0.5) THEN

      ! Get soil profile
      !--------------------------------------------
      ALLOCATE( rdum2d(ns_nc, nt_nc) )

      ! Soil layer thickness (m)
      status = NF90_GET_VAR (ncid, id_s_thk, file_dzsoil)

      ! Soil temperature
      status = nf90_get_var (ncid, id_st, rdum2d)
      file_st(:,1:ns_nc) = TRANSPOSE(rdum2d)

      ! Soil moisture
      status = nf90_get_var (ncid, id_sm, rdum2d)
      file_sm(:,1:ns_nc) = TRANSPOSE(rdum2d)

      DEALLOCATE( rdum2d )


      ! Calculate soil depths
      !--------------------------------------------
      zsoil(1)      = 0.0  ! UM soil depths
      file_zsoil(1) = 0.0  ! RACMO soil depths

        DO k=2, sm_lv+1
          zsoil(k) = zsoil(k-1) + dzsoil_jules(k-1)
        END DO
 
        DO k=2, ns_nc+1
          file_zsoil(k) = file_zsoil(k-1) + file_dzsoil(k-1)
        END DO

      ! Soil moisture profiles
      !--------------------------------------------
      ! Calculated UM total water(/m^2) in soil to depth of level
      ! RACMO provides water content in layer /m2,
      wat2dep  = 0.0
      DO k=2, sm_lv+1
        IF (zsoil(k) > file_zsoil(ns_nc+1)) THEN
          ! UM soil level is deeper than deepest RACMO soil level
          ! Sum water content in RACMO, assume water content /m2 is
          ! is constant after deepest RACMO soil level
          wat2dep(k) = SUM(file_dzsoil(:)*file_sm(obs_t0,:))               &
                     + (zsoil(k)-file_zsoil(ns_nc+1))*file_sm(obs_t0,ns_nc)

        ELSE IF (zsoil(k) < file_zsoil(2)) THEN
          ! UM level is shallower than shallowest RACMO level
          wat2dep(k) = zsoil(k)*file_sm(obs_t0,1)
        ELSE

          DO kk=2, ns_nc+1
            IF ( (zsoil(k) >   file_zsoil(kk-1)) .AND.                     &
                 (zsoil(k) <=  file_zsoil(kk))) THEN

              wat2dep(k) = SUM(file_dzsoil(1:kk-2)*file_sm(obs_t0,1:kk-2)) &
                         + (zsoil(k) - file_zsoil(kk-1))                   &
                         * file_sm(obs_t0,kk)
            END IF
          END DO

        END IF
      END DO

        DO k=1, sm_lv
          DO i=1, nlnd
            smcli_nc(i,k) = (wat2dep(k+1) - wat2dep(k)) / dzsoil_jules(k)
          END DO
        END DO

      DO i=1, nlnd
        smci_nc(i) = wat2dep(sm_lv+1)
      END DO



      ! Soil temperature profiles
      !--------------------------------------------
      ! Assuming RACMO soil temperatures are soil temperature at a soil
      ! depth rather than temperature of soil layer.

      DO k=2, st_lv+1
        DO i=1, nlnd
          IF (zsoil(k) > file_zsoil(ns_nc+1)) then
            ! UM soil level is deeper than deepest RACMO soil level
            ! Set Temperature to be same as lowest obs
            t_deep_soili_nc(i,k-1) = file_st(obs_t0,ns_nc)

          ELSE IF (zsoil(k) < file_zsoil(2)) then
            ! UM level is shallower than shallowest RACMO level
            t_deep_soili_nc(i,k-1) = tstari_nc                                &
                                   + (zsoil(k)/file_dzsoil(1))                &
                                   * (file_st(obs_t0,1) - tstari_nc)
          ELSE

            DO kk=2, ns_nc+1
              IF ( (zsoil(k) >   file_zsoil(kk-1)) .AND.                      &
                   (zsoil(k) <=  file_zsoil(kk))) THEN


                t_deep_soili_nc(i,k-1) = file_st(obs_t0,kk-1)                 &
                                + (zsoil(k) - file_zsoil(kk-1))               &
                                * (file_st(obs_t0,kk) - file_st(obs_t0,kk-1)) &
                                / (file_zsoil(kk) - file_zsoil(kk-1))

              END IF

            END DO
          END IF
        END DO
      END DO

      l_smcli        = .TRUE.
      l_smci         = .TRUE.
      l_t_deep_soili = .TRUE.

    END IF

    ! Get pressure co-ord levels
    !===========================
    status = NF90_GET_ATT (ncid, nf90_Global, 'coor_par_a', file_coor_par_a)
    status = NF90_GET_ATT (ncid, nf90_Global, 'coor_par_b', file_coor_par_b)

    ALLOCATE( rdum2d(nk_nc, nt_nc) )

    ! Get atmospheric profiles
    !=========================
    ! Temperature (K)
    status = NF90_GET_VAR (ncid, id_t, rdum2d)
    file_t(:,1:nk_nc) = TRANSPOSE(rdum2d)

    ! q mixing ratio (kg/kg)
    status = NF90_GET_VAR (ncid, id_q, rdum2d)
    file_q(:,1:nk_nc) = TRANSPOSE(rdum2d)

    ! U-wind (m/s)
    status = NF90_GET_VAR (ncid, id_u, rdum2d)
    file_u(:,1:nk_nc) = TRANSPOSE(rdum2d)

    ! V-wind (m/s)
    status = NF90_GET_VAR (ncid, id_v, rdum2d)
    file_v(:,1:nk_nc) = TRANSPOSE(rdum2d)

    ! omega (Pa/s)
    status = NF90_GET_VAR (ncid, id_omega, rdum2d)
    file_omega(:,1:nk_nc) = TRANSPOSE(rdum2d)


    ! Get forcing tendencies
    !=======================
    ! Horizontal advection + diffusion of temperature
    status = NF90_GET_VAR (ncid, id_tadv, rdum2d)
    file_tadv(:,1:nk_nc) = TRANSPOSE(rdum2d)

    ! Horizontal advection + diffusion of mixing ratio (q)
    status = NF90_GET_VAR (ncid, id_qadv, rdum2d)
    file_qadv(:,1:nk_nc) = TRANSPOSE(rdum2d)

    ! Horizontal advection + diffusion of U
    status = NF90_GET_VAR (ncid, id_uadv, rdum2d)
    file_uadv(:,1:nk_nc) = TRANSPOSE(rdum2d)

    ! Horizontal advection + diffusion of V
    status = NF90_GET_VAR (ncid, id_vadv, rdum2d)
    file_vadv(:,1:nk_nc) = TRANSPOSE(rdum2d)

    status = nf90_close(ncid)

    DEALLOCATE(rdum2d)


    ! Invert Vertical axis, since these were read in order of
    ! increasing pressure rather than in increasing height.
    file_coor_par_a(:) = file_coor_par_a(coor_lena:1:-1)
    file_coor_par_b(:) = file_coor_par_b(coor_lenb:1:-1)

    file_t (:,1:nk_nc) = file_t(:,nk_nc:1:-1)
    file_q (:,1:nk_nc) = file_q(:,nk_nc:1:-1)
    file_u (:,1:nk_nc) = file_u(:,nk_nc:1:-1)
    file_v (:,1:nk_nc) = file_v(:,nk_nc:1:-1)

    file_omega (:,1:nk_nc) = file_omega (:,nk_nc:1:-1)
    file_tadv  (:,1:nk_nc) = file_tadv  (:,nk_nc:1:-1)
    file_qadv  (:,1:nk_nc) = file_qadv  (:,nk_nc:1:-1)
    file_uadv  (:,1:nk_nc) = file_uadv  (:,nk_nc:1:-1)
    file_vadv  (:,1:nk_nc) = file_vadv  (:,nk_nc:1:-1)

    IF (lhook) CALL dr_hook('READ_RACMO_NC',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE read_racmo_nc

!==============================================================================

  SUBROUTINE racmo_nc_var_to_scm

    IMPLICIT NONE

    INTEGER :: i,j,k

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('RACMO_NC_VAR_TO_SCM',zhook_in,zhook_handle)

    DO i=1, nt_nc
      file_p_half(i,:) = file_coor_par_a(:)                        &
                       + file_coor_par_b(:) * file_p_srf(i)
    END DO

    DO k=1, nk_nc
      file_p_full(:,k) = 0.5*(file_p_half(:,k)                     &
                       +      file_p_half(:,k+1))
    END DO

    DO k=1, ns_nc
    END DO



    !=================================
    ! Add Surface values to the arrays
    !=================================
    !===========================================
    ! Scale/Convert to variables required by SCM
    !===========================================
    !======================
    ! Surface Timeseries
    !======================
    !======================
    ! Atmospheric Profiles
    !======================
    WHERE (file_q    == -0.0) file_q    = 0.0 ! Remove -0.0 values
    WHERE (file_qadv == -0.0) file_qadv = 0.0 ! Remove -0.0 values


    ! Mixing ratio ---> Specific humidity
    file_q(:,:)    =  file_q(:,:)    / (file_q(:,:)    + 1)
    file_qadv(:,:) =  file_qadv(:,:) / (file_qadv(:,:) + 1)

    ! Omega(Pa/s) --> w-wind(m/s) using -(Omega*R*T)/(P*g)
    ! Assuming surface pressure does not change

    DO j=1, nt_nc
      file_w(j,:) = - R * file_t(j,:) * file_omega(j,:)              &
                  / (file_p_full(j,:) * g)
    END DO

    WHERE (file_w == -0.0) file_w = 0.0 ! Stop creation of -0.0 values

    ! Forcing rates convert from 1/s to 1/day
    file_qadv(:,:) =  file_qadv(:,:)*86400.0
    file_tadv(:,:) =  file_tadv(:,:)*86400.0
    file_uadv(:,:) =  file_uadv(:,:)*86400.0
    file_vadv(:,:) =  file_vadv(:,:)*86400.0

    ! Surface fluxes in RACMO netcdf file are defined as a downward flux
    ! Reverse signs for SCM
    file_flx_h(:) =  -1.0 * file_flx_h(:)
    file_flx_e(:) =  -1.0 * file_flx_e(:)

    ! Calculate z-profile
    !====================

    !=================================================================!
    ! ASSUMING z-heights of forcing data does not vary significantly, !
    !          Use initial profile for interpolation on UM grid.      !
    !=================================================================!

    ! Calculate height on half levels
    ! Assume element 1 is surface since it is same pressure as p_srf

    ! Half Levels
    file_z_half_prf(1)   = 0.0

    DO k=2, nk_half_nc
      file_z_half_prf(k) = -(R/g) * file_t(obs_t0,k-1)        &
                         * (LOG(  file_p_half(obs_t0,k)       &
                                / file_p_half(obs_t0,k-1)))   &
                         + file_z_half_prf(k-1)
    END DO

    ! Full Levels

    file_z_full_prf(1)  = - (R/g) * file_t(obs_t0,1)          &
                        * (LOG(  file_p_full(obs_t0,1)        &
                               / file_p_srf(obs_t0)))

    DO k=2, nk_nc
      file_z_full_prf(k) = -(R/g)                             &
                         * ((  file_t(obs_t0,k)               &
                             + file_t(obs_t0,k-1))/2.0)       &
                         * (LOG(  file_p_full(obs_t0,k)       &
                                / file_p_full(obs_t0,k-1)))   &
                         + file_z_full_prf(k-1)
    END DO

    obs_top_nc = 32000.0
    obs_bot_nc = file_z_full_prf(1)

    IF (lhook) CALL dr_hook('RACMO_NC_VAR_TO_SCM',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE RACMO_nc_var_to_scm

!==============================================================================

  SUBROUTINE RACMO_nc_grid_to_scm

    IMPLICIT NONE

    INTEGER :: i,j,k,kk
    REAL, ALLOCATABLE :: exner_th(:,:,:), exner_rh(:,:,:)

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('RACMO_NC_GRID_TO_SCM',zhook_in,zhook_handle)

    ALLOCATE( exner_th(rw_lng,rw,nmod_lv)   &
            , exner_rh(rw_lng,rw,nmod_lv+1))

    DO j=1, rw
      DO i=1, rw_lng
        shflx_nc(i,j,:) = file_flx_H  (obs_t0:nt_nc) ! Sen. heat flux (W/m2)
        lhflx_nc(i,j,:) = file_flx_E  (obs_t0:nt_nc) ! Lat. heat flux (W/m2)
        tstar_nc(i,j,:) = file_t_skin (obs_t0:nt_nc) ! (K)
      END DO
    END DO

    l_shflx     = .TRUE.
    l_lhflx     = .TRUE.
    l_tstar_for = .TRUE.

    ! No ozone provided
    l_ozone     = .FALSE.

    file_exner_full(:,:) = (file_p_full(:,:)/p_zero)**(kappa)



    DO j=1, rw
      DO i=1, rw_lng
        CALL interp1d(file_z_full_prf, file_exner_full(obs_t0,:),  &
                      z_th, exner_th(i,j,:))
        CALL interp1d(file_z_full_prf, file_exner_full(obs_t0,:),  &
                      z_rh, exner_rh(i,j,:))

        CALL interp1d(file_z_full_prf, file_t(obs_t0,:),      &
                      z_rh,ti_rh_nc(i,j,:))

        CALL interp1d(file_z_full_prf, file_t(obs_t0,:),      &
                      z_th,ti_th_nc(i,j,:))

        CALL interp1d(file_z_full_prf, file_q(obs_t0,:),      &
                      z_th,qi_nc(i,j,:)        )

        CALL interp1d(file_z_full_prf, file_u(obs_t0,:),      &
                      z_rh(1:nmod_lv), ui_nc(i,j,:))

        CALL interp1d(file_z_full_prf, file_v(obs_t0,:),      &
                      z_rh(1:nmod_lv), vi_nc(i,j,:))

        CALL interp1d(file_z_full_prf, file_w(obs_t0,:),      &
                      z_th,wi_nc(i,j,1:nmod_lv))

      END DO
    END DO

   ! Calculate File  theta and p
   WHERE (exner_rh /= rmdi) pi_rh_nc = p_zero * exner_rh ** recip_kappa
   WHERE ((exner_th /= rmdi) .AND. &
          (ti_th_nc /= rmdi)) thi_nc = ti_th_nc / exner_th

    l_qi = .TRUE.
    l_pi = .TRUE.
    l_ui = .TRUE.
    l_vi = .TRUE.
    l_wi = .TRUE.
    l_theta = .TRUE.

    DO kk=obs_t0, nt_nc
      DO j=1, rw
        DO i=1, rw_lng
          ! obs forcing
          CALL interp1d(file_z_full_prf, file_t(kk,:),                 &
                        z_th, t_bg_nc(i,j,kk-obs_t0+1,:))
          CALL interp1d(file_z_full_prf, file_q(kk,:),                 &
                        z_th, q_bg_nc(i,j,kk-obs_t0+1,:))
          CALL interp1d(file_z_full_prf, file_u(kk,:),                 &
                        z_rh(1:nmod_lv),                               &
                        u_bg_nc(i,j,kk-obs_t0+1,:))
          CALL interp1d(file_z_full_prf, file_v(kk,:),                 &
                        z_rh(1:nmod_lv),                               &
                        v_bg_nc(i,j,kk-obs_t0+1,:))




          CALL interp1d(file_z_full_prf, file_tadv(kk,:),              &
                        z_th, t_h_inc_nc(i,j,kk-obs_t0+1,:))

          CALL interp1d(file_z_full_prf, file_qadv(kk,:),              &
                        z_th, q_h_inc_nc(i,j,kk-obs_t0+1,:))

          CALL interp1d(file_z_full_prf, file_uadv(kk,:),              &
                        z_rh(1:nmod_lv),                               &
                        u_h_inc_nc(i,j,kk-obs_t0+1,:))

          CALL interp1d(file_z_full_prf, file_vadv(kk,:),              &
                        z_rh(1:nmod_lv),                               &
                        v_h_inc_nc(i,j,kk-obs_t0+1,:))

          CALL interp1d(file_z_full_prf, file_w(kk,:),                 &
                        z_th,                                          &
                        w_bg_nc(i,j,kk-obs_t0+1,1:nmod_lv))

        END DO
      END DO
    END DO

    l_t_bg   = .TRUE.
    l_q_bg   = .TRUE.
    l_u_bg   = .TRUE.
    l_v_bg   = .TRUE.
    l_w_bg   = .TRUE.
    l_t_inc  = .TRUE.
    l_q_inc  = .TRUE.
    l_u_inc  = .TRUE.
    l_v_inc  = .TRUE.


    ! Set t_v_inc_nc = 0 as file_tadv is horizontal
    ! advection/diffusion

    t_v_inc_nc(:,:,:,:) = 0.0
    q_v_inc_nc(:,:,:,:) = 0.0
    u_v_inc_nc(:,:,:,:) = 0.0
    v_v_inc_nc(:,:,:,:) = 0.0

    ! No LS tendency + diffusion of W in driver file, so calculate
    ! w_inc from background w field
    DO kk=obs_t0, nt_nc-1
      w_inc_nc(:,:,kk-obs_t0+1,:) = (86400.0/obs_pd_nc)          &
                                  * ( w_bg_nc(:,:,kk-obs_t0+2,:) &
                                    - w_bg_nc(:,:,kk-obs_t0+1,:))
    END DO

    ! Because there is no w-field information after end of driver
    ! file, dw/dt cannot be calculated from the field on the
    ! last timestep, assume that it does not change, so:
    w_inc_nc(:,:,nt_nc-obs_t0+1,:) = w_inc_nc(:,:,nt_nc-obs_t0,:)

    ! Set surface level of w fields to 0
    wi_nc(:,:,0)      = 0.0
    w_inc_nc(:,:,:,0) = 0.0
    w_bg_nc(:,:,:,0)  = 0.0


    l_w_inc  = .TRUE.
    print*,' No LS W-tendency in driver file, Forcing calculated from w-field.'

    DEALLOCATE(exner_th, exner_rh)


    IF (lhook) CALL dr_hook('RACMO_NC_GRID_TO_SCM',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE RACMO_nc_grid_to_scm

!==============================================================================
END MODULE RACMO_netCDF
