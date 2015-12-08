
MODULE TWPICE_netcdf


! Description:
!
!
!
!
!
! *********************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *********************************COPYRIGHT*********************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model


  ! Modules
  USE netcdf
  USE scm_utils
  USE s_interp_mod, ONLY: interp1d
  USE s_main_force, ONLY: netcdf_file
  USE netCDF_arr
  USE global_SCMop, ONLY: incdf
  USE MCC_data
  USE TWPICE_data

  USE earth_constants_mod,  ONLY: g

  IMPLICIT NONE

  REAL ::           &
    file_lat        &
  , file_long       &
  , file_albsoil

  REAL, ALLOCATABLE, DIMENSION(:,:) :: &
    file_t          &! Temperature
  , file_q          &! H2o mixing ratio
  , file_u          &! Zonal wind
  , file_v          &! Meridional wind
  , file_omega      &! Rate of change of pressure (dp/dt)
  , file_divt       &! Hor.  adv. of s (dry static temperature)
  , file_divq       &! Hor.  adv. of q
  , file_vertdivt   &! Vert. adv. of s (dry static temperature)
  , file_vertdivq   &! Vert. adv. of q
  , file_w          &! Vertical wind         (Calculated)
  , file_th          ! Potential temperature (Calculated)

  REAL, ALLOCATABLE, DIMENSION(:)   :: &
    file_time       &! Time levels
  , file_year       &
  , file_month      &
  , file_day        &
  , file_hour       &
  , file_min        &
  , file_z_prf      &! z of obs levels (assumed fixed at initial run time)
  , file_layer_thk  &! Layer thickness of obs level (Using mean T of layer)
  , file_flx_H      &! Srf. sensible heat flux
  , file_flx_E      &! Srf. latent   heat flux
  , file_tsskin     &! Srf. Skin temperature
  , file_ts_air     &! Srf. Air temperature
  , file_p_srf      &! Srf. Mean Pressure
  , file_omega_srf  &! Srf. Rate of change of pressure (dp/dt)
  , file_qs_srf     &! Srf. Saturation mixing ratio
  , file_u_srf      &! Srf. Zonal wind
  , file_v_srf      &! Srf. Meridional wind
  , file_rhs_air    &! Srf. Air Relative humidity
  , file_w_srf      &! Srf. Vertical velocity (Calculated)
  , file_th_srf     &! Srf. Potential temperature (Calculated)
  , file_q_srf      &! Srf. Mixing Ratio (Calculated)
  , file_p           ! Pressure levels


  INTEGER(incdf) :: &
    status          &! NetCDF Error status
  , ncid            &! NetCDF File id
  , time_dimid      &! Time  dimension id
  , nlev_dimid      &! Nlevs dimension id
  , slev_dimid      &! n soil levels dimension id
  , xtype           &
  , nt_nc           &! Number of time levels in file
  , nk_in           &! Number of levels in file
  , attnum          &!
  , coor_lena       &!
  , coor_lenb       &!
  , nk_half_in

  ! Variable Ids for:
  INTEGER(incdf) :: &
    time_id         &! Time (days from 2005-12-31, 00:00:00)
  , date_id         &!
  , year_id         &!
  , mon_id          &!
  , day_id          &!
  , hour_id         &!
  , min_id          &!
  , sec_id          &!
  , lat_id          &! Latitude
  , lon_id          &! Longitude
  , p_in_id         &! Pressure
  , t_id            &! Temperature
  , q_id            &! H2o vapour Mixing ratio
  , ql_id           &! H2o Liquid Mixing Ratio (kg/kg)
  , qi_id           &! H2o Ice    Mixing Ratio (kg/kg)
  , u_id            &! Zonal wind
  , v_id            &! Meridional wind
  , omega_id        &! Rate of change of pressure (dp/dt)
  , hdivt_id        &! Hor.  adv. of s (dry static temperature)
  , hdivq_id        &! Hor.  adv. of q
  , vdivt_id        &! Vert. adv. of s (dry static temperature)
  , vdivq_id        &! Vert. adv. of q
  , alb_id          &! Albedo
  , alb_snow_id     &! Snow Albedo
  , flx_H_id        &! Srf. sensible heat flux
  , flx_E_id        &! Srf. latent   heat flux
  , open_sst_id     &!
  , orog_id         &!
  , t_soil_id       &!
  , q_soil_id       &!
  , z0m_id          &!
  , z0h_id          &!
  , tsskin_id       &! Srf. Skin temperature
  , ts_air_id       &! Srf. Air temperature
  , ps_id           &! Srf. Mean Pressure
  , omegas_id       &! Srf. Rate of change of pressure (dp/dt)
  , qs_srf_id       &! Srf. Saturation mixing ratio
  , u_srf_id        &! Srf. Zonal wind
  , v_srf_id        &! Srf. Meridional wind
  , land_sea_id     &!
  , rhs_air_id       ! Srf. Air Relative humidity

  REAL, PRIVATE, ALLOCATABLE :: dummy(:,:)

CONTAINS
!============================================================================

  SUBROUTINE alloc_twpice_nc

    IMPLICIT NONE

    status = NF90_OPEN(TRIM(ADJUSTL(netcdf_file)), nf90_noWrite, ncid)

    IF (status /= nf90_noerr) THEN
      PRINT*,'Error opening netcdf file: '//TRIM(ADJUSTL(netcdf_file))
    END IF

    ! Get NetCDF file Dimension ID numbers
    status = NF90_INQ_DIMID (ncid, 'time', time_dimid)
    status = NF90_INQ_DIMID (ncid, 'lev',  nlev_dimid)


    ! Get Dimension sizes
    status = NF90_INQUIRE_DIMENSION (ncid, time_dimid, len=nt_nc)
    status = NF90_INQUIRE_DIMENSION (ncid, nlev_dimid, len=nk_in)


    ! Allocate arrays to receive obs data 0 element to hold srf face values
    ! 1D arrays
    ALLOCATE                   &
      ( file_flx_h     (nt_nc) &
      , file_flx_e     (nt_nc) &
      , file_tsskin    (nt_nc) &
      , file_time      (nt_nc) &
      , file_year      (nt_nc) &
      , file_month     (nt_nc) &
      , file_day       (nt_nc) &
      , file_hour      (nt_nc) &
      , file_min       (nt_nc) &
      , file_ts_air    (nt_nc) &
      , file_th_srf    (nt_nc) &
      , file_p_srf     (nt_nc) &
      , file_omega_srf (nt_nc) &
      , file_qs_srf    (nt_nc) &
      , file_q_srf     (nt_nc) &
      , file_rhs_air   (nt_nc) &
      , file_u_srf     (nt_nc) &
      , file_v_srf     (nt_nc) &
      , file_w_srf     (nt_nc) &
      , file_p       (0:nk_in) & ! srf calc'd depends on starting point
      , file_z_prf   (0:nk_in) & ! add srf, i.e. 0m
      , file_layer_thk (nk_in) )

    ! 2D Arrays
    ! (extra k level so to add srf value to array for interpolation)
    ALLOCATE                           &
      ( file_divT     (nt_nc, 0:nk_in) &! Assume srf value = lowest obs
      , file_vertdivT (nt_nc, 0:nk_in) &! Assume srf value = lowest obs
      , file_divq     (nt_nc, 0:nk_in) &! Assume srf value = lowest obs
      , file_vertdivq (nt_nc, 0:nk_in) &! Assume srf value = lowest obs
      , file_u        (nt_nc, 0:nk_in) &
      , file_v        (nt_nc, 0:nk_in) &
      , file_w        (nt_nc, 0:nk_in) &
      , file_t        (nt_nc, 0:nk_in) &
      , file_th       (nt_nc, 0:nk_in) &
      , file_q        (nt_nc, 0:nk_in) &
      , file_omega    (nt_nc, 0:nk_in) )

    status = NF90_CLOSE(ncid)

    ! Initialise Arrays
    file_flx_h     = RMDI  ;  file_flx_e  = RMDI  ;  file_tsskin = RMDI
    file_time      = RMDI  ;  file_year   = RMDI  ;  file_month  = RMDI
    file_day       = RMDI  ;  file_hour   = RMDI  ;  file_min    = RMDI
    file_ts_air    = RMDI  ;  file_th_srf = RMDI  ;  file_p_srf  = RMDI
    file_omega_srf = RMDI  ;  file_qs_srf = RMDI  ;  file_q_srf  = RMDI
    file_rhs_air   = RMDI  ;  file_u_srf  = RMDI  ;  file_v_srf  = RMDI
    file_w_srf     = RMDI  ;  file_p      = RMDI  ;  file_z_prf  = RMDI
    file_layer_thk = RMDI


    file_divt      = RMDI  ;  file_vertdivt = RMDI ; file_divq   = RMDI
    file_vertdivq  = RMDI  ;  file_u        = RMDI ; file_v      = RMDI
    file_w         = RMDI  ;  file_t        = RMDI ; file_th     = RMDI
    file_q         = RMDI  ;  file_omega    = RMDI


  END SUBROUTINE alloc_twpice_nc

!============================================================================

  SUBROUTINE dealloc_twpice_nc

    IMPLICIT NONE

    DEALLOCATE(                                                      &
      ! 1D Arrays
      file_flx_H,     file_flx_E,     file_tsskin,    file_time,     &
      file_year,      file_month,     file_day,       file_hour,     &
      file_min,       file_ts_air,    file_th_srf,    file_p_srf,    &
      file_omega_srf, file_qs_srf,    file_q_srf,     file_rhs_air,  &
      file_u_srf,     file_v_srf,     file_w_srf,     file_p,        &
      file_z_prf,     file_layer_thk,                                &

      ! 2D Arrays
      file_divT,      file_vertdivT,  file_divq,      file_vertdivq, &
      file_u,         file_v,         file_w,         file_t,        &
      file_th,        file_q,         file_omega)

  END SUBROUTINE dealloc_twpice_nc

!============================================================================

 SUBROUTINE read_TWPICE_nc

      IMPLICIT NONE

      status = nf90_open(TRIM(ADJUSTL(netcdf_file)), nf90_noWrite, ncid)

      IF (status /= nf90_noerr) THEN
        PRINT*,'Error opening netcdf file: '//TRIM(ADJUSTL(netcdf_file))
      END IF

      ! Setup variables
      status = NF90_INQ_VARID (ncid, 'time',      time_id)
      status = NF90_INQ_VARID (ncid, 'year',      year_id)
      status = NF90_INQ_VARID (ncid, 'month',      mon_id)
      status = NF90_INQ_VARID (ncid, 'day',        day_id)
      status = NF90_INQ_VARID (ncid, 'hour',      hour_id)
      status = NF90_INQ_VARID (ncid, 'minute',     min_id)
      status = NF90_INQ_VARID (ncid, 'lat',        lat_id)
      status = NF90_INQ_VARID (ncid, 'lon',        lon_id)

      ! Forcing variables
      status = NF90_INQ_VARID (ncid, 'lev',       p_in_id)
      status = NF90_INQ_VARID (ncid, 'T',            t_id)
      status = NF90_INQ_VARID (ncid, 'q',            q_id)
      status = NF90_INQ_VARID (ncid, 'u',            u_id)
      status = NF90_INQ_VARID (ncid, 'v',            v_id)
      status = NF90_INQ_VARID (ncid, 'omega',    omega_id)
      status = NF90_INQ_VARID (ncid, 'divT',     hdivt_id)
      status = NF90_INQ_VARID (ncid, 'divq',     hdivq_id)
      status = NF90_INQ_VARID (ncid, 'vertdivT', vdivt_id)
      status = NF90_INQ_VARID (ncid, 'vertdivq', vdivq_id)
      status = NF90_INQ_VARID (ncid, 'alb',        alb_id)
      status = NF90_INQ_VARID (ncid, 'shflx',    flx_H_id)
      status = NF90_INQ_VARID (ncid, 'lhflx',    flx_E_id)
      status = NF90_INQ_VARID (ncid, 'Tskin',   tsskin_id)
      status = NF90_INQ_VARID (ncid, 'Tsair',   ts_air_id)
      status = NF90_INQ_VARID (ncid, 'Ps',          ps_id)
      status = NF90_INQ_VARID (ncid, 'omegas',  omegas_id)
      status = NF90_INQ_VARID (ncid, 'qs',      qs_srf_id)
      status = NF90_INQ_VARID (ncid, 'RHsair', rhs_air_id)
      status = NF90_INQ_VARID (ncid, 'usrf',     u_srf_id)
      status = NF90_INQ_VARID (ncid, 'vsrf',     v_srf_id)


      !=======================
      ! Get obs data from file
      !=======================

      ! Get set-up variables
      !=======================
      status = NF90_GET_VAR (ncid, lat_id,  file_lat)
      status = NF90_GET_VAR (ncid, lon_id,  file_long)
      status = NF90_GET_VAR (ncid, alb_id,  file_albsoil)
      status = NF90_GET_VAR (ncid, year_id, file_year)
      status = NF90_GET_VAR (ncid, mon_id,  file_month)
      status = NF90_GET_VAR (ncid, day_id,  file_day)
      status = NF90_GET_VAR (ncid, hour_id, file_hour)
      status = NF90_GET_VAR (ncid, min_id,  file_min)



      ! Get surface timeseries of:
      !===========================
      status = NF90_GET_VAR (ncid, flx_H_id,   file_flx_H)     ! (W/m2)
      status = NF90_GET_VAR (ncid, flx_E_id,   file_flx_E)     ! (W/m2)
      status = NF90_GET_VAR (ncid, tsskin_id,  file_tsskin)    ! (K)
      status = NF90_GET_VAR (ncid, qs_srf_id,  file_qs_srf)    ! (kg/kg)
      status = NF90_GET_VAR (ncid, rhs_air_id, file_rhs_air)   ! (%)
      status = NF90_GET_VAR (ncid, ts_air_id,  file_ts_air)    ! (K)
      status = NF90_GET_VAR (ncid, ps_id,      file_p_srf)     ! (mb)
      status = NF90_GET_VAR (ncid, u_srf_id,   file_u_srf)     ! (m/s)
      status = NF90_GET_VAR (ncid, v_srf_id,   file_v_srf)     ! (m/s)
      status = NF90_GET_VAR (ncid, omegas_id,  file_omega_srf) ! (mb/hr)



      ! Get time levels
      !=========================

      status = NF90_GET_VAR (ncid, time_id, file_time) ! Time (days since
                                                       ! 2005-12-31, 00:00:00)


      obs_pd_nc = (file_time(2) - file_time(1))*86400.0
      yeari_nc  = file_year(obs_t0)
      monthi_nc = file_month(obs_t0)
      dayi_nc   = file_day(obs_t0)
      houri_nc  = file_hour(obs_t0)
      mini_nc   = file_min(obs_t0)

      albsoil_nc = file_albsoil
      lat_nc     = file_lat
      long_nc    = file_long
      tstari_nc  = file_tsskin(obs_t0)

!      ndayin = INT((  time_in(nt_nc)                            &
!                    - time_in(obs_t0)))

!      nminin = INT((  time_in(nt_nc)                            &
!                    - time_in(obs_t0)                           &
!                    - ndayin)*24.0*60.0)

!      nsecin = INT((((  time_in(nt_nc)                          &
!                      - time_in(obs_t0)                         &
!                      - ndayin)*24.0*60.0)                      &
!                      - nminin)*60.0)


      ! Get pressure co-ord levels
      !===========================
      status = NF90_GET_VAR (ncid, p_in_id, file_p(1:nk_in)) ! Pressure(mb)

      ALLOCATE( dummy(nk_in, nt_nc) )

      ! Get atmospheric profiles
      !=========================

      ! Temperature (K)
      status = NF90_GET_VAR (ncid, T_id, dummy)
      file_t(:,1:nk_in)  = TRANSPOSE(dummy)

      ! q mixing ratio (g/kg)
      status = NF90_GET_VAR (ncid, q_id, dummy)
      file_q(:,1:nk_in) = TRANSPOSE(dummy)

      ! U-wind (m/s)
      status = NF90_GET_VAR (ncid, u_id, dummy)
      file_u(:,1:nk_in) = TRANSPOSE(dummy)

      ! V-wind (m/s)
      status = NF90_GET_VAR (ncid, v_id, dummy)
      file_v(:,1:nk_in) = TRANSPOSE(dummy)

      ! omega (mb/hr)
      status = NF90_GET_VAR (ncid, omega_id, dummy)
      file_omega(:,1:nk_in) = TRANSPOSE(dummy)


      ! Get Forcing tendencies
      !=======================

      ! Horizontal temperature advection
      status = NF90_GET_VAR (ncid, hdivt_id, dummy)
      file_divT(:,1:nk_in) = TRANSPOSE(dummy)

      ! Vertical dry static energy(s) advection
      status = NF90_GET_VAR (ncid, vdivt_id, dummy)
      file_vertdivT(:,1:nk_in) = TRANSPOSE(dummy)

      ! Horizontal mixing ratio (q) advection
      status = NF90_GET_VAR (ncid, hdivq_id, dummy)
      file_divq(:,1:nk_in) = TRANSPOSE(dummy)

      ! Vertical mixing ratio (q) advection
      status = NF90_GET_VAR (ncid, vdivq_id, dummy)
      file_vertdivq(:,1:nk_in) = TRANSPOSE(dummy)

      DEALLOCATE(dummy)

      status = nf90_close(ncid)

      ! Invert Vertical axis, since these were read in order of
      ! increasing pressure rather than in increasing height.

      file_p (1:nk_in)   = file_p(  nk_in:1:-1)
      file_t (:,1:nk_in) = file_t(:,nk_in:1:-1)
      file_q (:,1:nk_in) = file_q(:,nk_in:1:-1)
      file_u (:,1:nk_in) = file_u(:,nk_in:1:-1)
      file_v (:,1:nk_in) = file_v(:,nk_in:1:-1)

      file_omega    (:,1:nk_in) = file_omega    (:,nk_in:1:-1)
      file_divT     (:,1:nk_in) = file_divT     (:,nk_in:1:-1)
      file_vertdivT (:,1:nk_in) = file_vertdivT (:,nk_in:1:-1)
      file_divq     (:,1:nk_in) = file_divq     (:,nk_in:1:-1)
      file_vertdivq (:,1:nk_in) = file_vertdivq (:,nk_in:1:-1)

    END SUBROUTINE read_TWPICE_nc

!==============================================================================

!==============================================================================

    SUBROUTINE TWPICE_nc_var_to_scm

      IMPLICIT NONE

      INTEGER :: i,j,k

      !=================================
      ! Add Surface values to the arrays
      !=================================
      file_p(0)          = file_p_srf(obs_t0)         ! Pressure (Pa)
      file_t(:,0)        = file_ts_air                ! Temperature (K)
      file_q(:,0)        = (file_rhs_air(:)/100.0) &  ! q mix rat. (g/kg)
                         * file_qs_srf(:)*1000.0
      file_u(:,0)        = file_u_srf(:)              ! U-wind (m/s)
      file_v(:,0)        = file_v_srf(:)              ! V-wind (m/s)

      file_omega(:,0)    = file_omega_srf(:)          ! omega (mb/hr)
      file_divT(:,0)     = file_divT(:,1)     ! Hor. temperature adv.
      file_vertdivT(:,0) = file_vertdivT(:,1) ! Ver. dry static energy(s) adv.
      file_divq(:,0)     = file_divq(:,1)     ! Hor. mixing ratio (q) adv.
      file_vertdivq(:,0) = file_vertdivq(:,1) ! Ver. mixing ratio (q) adv.


      !===========================================
      ! Scale/Convert to variables required by SCM
      !===========================================

      ! Pressure (Pa)
      file_p(:) = file_p(:) * 100.0

      ! Mixing ratio: g/kg --> kg/kg
      file_q(:,:)        = 1.0e-3 * file_q(:,:)         ! Note: File comments
      file_divq(:,:)     = 1.0e-3 * file_divq(:,:)      ! may be wrong, file
      file_vertdivq(:,:) = 1.0e-3 * file_vertdivq(:,:)  ! specifies kg/kg but
                                                        ! appear to be in g/kg.
      WHERE (file_q == -0.0)        file_q = 0.0        ! Remove -0.0 values
      WHERE (file_divq == -0.0)     file_divq = 0.0     ! Remove -0.0 values
      WHERE (file_vertdivq == -0.0) file_vertdivq = 0.0 ! Remove -0.0 values

      ! Mixing ratio ---> Specific humidity
      file_q(:,:)        =  file_q(:,:)        / (file_q(:,:)        + 1)
      file_divq(:,:)     =  file_divq(:,:)     / (file_divq(:,:)     + 1)
      file_vertdivq(:,:) =  file_vertdivq(:,:) / (file_vertdivq(:,:) + 1)

      ! Omega(mb/hr) --> w-wind(m/s) using -Omega*R*T/P
      ! Assuming surface pressure does not change
      DO j=1, nt_nc
        file_w(j,:) = - R * file_t(j,:)                               &
                    * file_omega(j,:)*(1.0/36.0)                      &
                    / (file_p(:) * g )
      END DO

      WHERE (file_w == -0.0) file_w = 0.0 ! Stop creation of -0.0 values

      ! Temperature  --> Theta
      ! Assumeing surface pressure does not change
      DO j=1, nt_nc
        file_th(j,:) = file_t(j,:)                                    &
                      * ((1.0e+5/file_p(:))**(R/cp))
      END DO


      ! Calculate z-profile
      !====================

      !=================================================================!
      ! ASSUMING z-heights of forcing data does not vary significantly, !
      !          Use initial profile for interpolation on UM grid.      !
      !=================================================================!

      ! Set Surface height to 0m
      file_z_prf(0)    = 0.0

      DO k=1, nk_in
        file_layer_thk(k) = -(R/g)                                       &
                          * ((  file_t(obs_t0,k)                         &
                            + file_t(obs_t0,k-1))/2.0)                   &
                            * (LOG(  file_p(k)                           &
                                   / file_p(k-1)))

        file_z_prf(k)     = file_z_prf(k-1) + file_layer_thk(k)
      END DO

    END SUBROUTINE TWPICE_nc_var_to_scm

!==============================================================================

    SUBROUTINE TWPICE_nc_grid_to_scm

      IMPLICIT NONE

      INTEGER :: i,j,k,kk

      !============================================================
      ! Ozone
      ! Ozone is dealt with separately as it is not containing in the
      ! TWPice netcdf files

!!$      ! Forcing data on UM levels
!!$      ! Initial profiles
!!$      real, intent(inout) :: init_th(:)
!!$      real, intent(inout) :: init_p(:)
!!$      real, intent(inout) :: init_q(:)
!!$      real, intent(inout) :: init_u(:)
!!$      real, intent(inout) :: init_v(:)
!!$      real, intent(inout) :: init_w(:)
!!$      real, intent(inout) :: ozone(:)
!!$
!!$      ! obs forcing
!!$      real, intent(inout) :: t_h_inc(:,:)
!!$      real, intent(inout) :: t_v_inc(:,:)
!!$      real, intent(inout) :: q_h_inc(:,:)
!!$      real, intent(inout) :: q_v_inc(:,:)
!!$      real, intent(inout) :: u_inc(:,:)
!!$      real, intent(inout) :: v_inc(:,:)
!!$      real, intent(inout) :: w_inc(:,:)

       DO j=1, rw
         DO i=1, rw_lng
          shflx_nc(i,j,:) = file_flx_H  (obs_t0:nt_nc) ! Sen. heat flux (W/m2)
          lhflx_nc(i,j,:) = file_flx_E  (obs_t0:nt_nc) ! Lat. heat flux (W/m2)
          tstar_nc(i,j,:) = file_tsskin (obs_t0:nt_nc) ! (K)
          CALL interp1d(twp_data_o3_z, twp_data_o3, z_th, o3_nc(i,j,:))
        END DO
      END DO

      IF (MAXVAL(z_th) > MAXVAL(twp_data_o3_z)) THEN
        CALL alloc_mcc(nmod_lv)
        CALL interp1d(mcc_trp_z, mcc_trp_o3, z_th, mcc_o3)
        DO k=1, nmod_lv
          DO j=1, rw
            DO i=1, rw_lng
              IF (o3_nc(i,j,k) == RMDI) o3_nc(i,j,k) = mcc_o3(k)
            END DO
          END DO
        END DO
        CALL dealloc_mcc
      END IF

      DO j=1, rw
        DO i=1, rw_lng

          CALL interp1d(file_z_prf, file_p(:),             &
                        z_rh,       pi_rh_nc(i,j,:))
          CALL interp1d(file_z_prf, file_q(obs_t0,:),      &
                        z_th,       qi_nc(i,j,:))
          CALL interp1d(file_z_prf, file_th(obs_t0,:),     &
                        z_th,       thi_nc(i,j,:))
          CALL interp1d(file_z_prf, file_u(obs_t0,:),      &
                        z_rh(1:nmod_lv), ui_nc(i,j,:))
          CALL interp1d(file_z_prf, file_v(obs_t0,:),      &
                        z_rh(1:nmod_lv), vi_nc(i,j,:))
          CALL interp1d(file_z_prf, file_w(obs_t0,:),      &
                        z_th,       wi_nc(i,j,1:nmod_lv))
        END DO
      END DO

      DO kk=obs_t0, nt_nc
        DO j=1, rw
          DO i=1, rw_lng
            ! obs forcing
            CALL interp1d(file_z_prf, file_divt(kk,:),              &
                          z_th, t_h_inc_nc(i,j,kk-obs_t0+1,:))

            CALL interp1d(file_z_prf, file_divq(kk,:),              &
                          z_th, q_h_inc_nc(i,j,kk-obs_t0+1,:))

            CALL interp1d(file_z_prf, file_vertdivt(kk,:),          &
                          z_th, t_v_inc_nc(i,j,kk-obs_t0+1,:))

            CALL interp1d(file_z_prf, file_vertdivq(kk,:),          &
                          z_th, q_v_inc_nc(i,j,kk-obs_t0+1,:))

            CALL interp1d(file_z_prf, file_u(kk,:),                 &
                          z_rh(1:nmod_lv),                          &
                          u_h_inc_nc(i,j,kk-obs_t0+1,:))

            CALL interp1d(file_z_prf, file_v(kk,:),                 &
                          z_rh(1:nmod_lv),                          &
                          v_h_inc_nc(i,j,kk-obs_t0+1,:))

            CALL interp1d(file_z_prf, file_w(kk,:),                 &
                          z_th,                                     &
                          w_inc_nc(i,j,kk-obs_t0+1,:))
          END DO
        END DO
      END DO


      ! UM resolution is outside current obs/forcing data
      ! Use McClatchly routines to fill in missing domain

      IF (MAXVAL(z_rh) > MAXVAL(file_z_prf)) THEN

        CALL alloc_mcc(nmod_lv)

        CALL interp1d(mcc_trp_z, mcc_trp_p,  z_rh, mcc_rh_p)
        CALL interp1d(mcc_trp_z, mcc_trp_th, z_th, mcc_th)
        CALL interp1d(mcc_trp_z, mcc_trp_q,  z_th, mcc_q)

        DO k=1, nmod_lv
          DO j=1, rw
            DO i=1, rw_lng
              IF (thi_nc (i,j,k) == RMDI) thi_nc(i,j,k) = mcc_th(k)
            END DO
          END DO
        END DO

        DO k=1, nmod_lv+1
          DO j=1, rw
            DO i=1, rw_lng
              IF (pi_rh_nc(i,j,k) == RMDI) pi_rh_nc(i,j,k) = mcc_rh_p(k)
            END DO
          END DO
        END DO

        DO k=1, nwet_lv
          DO j=1, rw
            DO i=1, rw_lng
              IF (qi_nc(i,j,k) == RMDI) qi_nc(i,j,k) = mcc_q(k)
            END DO
          END DO
        END DO

        CALL dealloc_mcc


        DO k=1, nmod_lv
          DO j=1, rw
            DO i=1, rw_lng
              IF (ui_nc(i,j,k) == RMDI) ui_nc(i,j,k) = ui_nc(i,j,k-1)
              IF (vi_nc(i,j,k) == RMDI) vi_nc(i,j,k) = vi_nc(i,j,k-1)
              IF (wi_nc(i,j,k) == RMDI) wi_nc(i,j,k) = wi_nc(i,j,k-1)
            END DO
          END DO
        END DO

        DO k=1, nmod_lv
          DO kk=1, nobs
            DO j=1, rw
              DO i=1, rw_lng
                IF (t_h_inc_nc(i,j,kk,k) == RMDI)                             &
                    t_h_inc_nc(i,j,kk,k) =  t_h_inc_nc(i,j,kk,k-1)

                IF (q_h_inc_nc(i,j,kk,k) == RMDI)                             &
                    q_h_inc_nc(i,j,kk,k) =  q_h_inc_nc(i,j,kk,k-1)

                IF (t_v_inc_nc(i,j,kk,k) == RMDI)                             &
                    t_v_inc_nc(i,j,kk,k) =  t_v_inc_nc(i,j,kk,k-1)

                IF (q_v_inc_nc(i,j,kk,k) == RMDI)                             &
                    q_v_inc_nc(i,j,kk,k) =  q_v_inc_nc(i,j,kk,k-1)

                IF (u_h_inc_nc(i,j,kk,k) == RMDI)                             &
                    u_h_inc_nc(i,j,kk,k) =  u_h_inc_nc(i,j,kk,k-1)

                IF (v_h_inc_nc(i,j,kk,k) == RMDI)                             &
                    v_h_inc_nc(i,j,kk,k) =  v_h_inc_nc(i,j,kk,k-1)

                IF (w_inc_nc(i,j,kk,k) == RMDI)                               &
                    w_inc_nc(i,j,kk,k) =  w_inc_nc(i,j,kk,k-1)
              END DO
            END DO
          END DO
        END DO

      END IF

    END SUBROUTINE TWPICE_nc_grid_to_scm                                      !

!==============================================================================

END MODULE TWPICE_netCDF

