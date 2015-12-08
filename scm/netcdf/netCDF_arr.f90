! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE netCDF_arr

  USE scm_utils, ONLY:                                                  &
      nobs, nmod_lv, nwet_lv, rw, rw_lng, nlnd, ntile, ntype            &
    , sm_lv, st_lv, rmdi, imdi, zhook_in, zhook_out, jprb, lhook        &
    , dr_hook


  IMPLICIT NONE

!-----------------------------------------------------------------------------
!
! Description:
!   Declares registers and arrays to hold forcing data from netcdf file
!
! Method:
!   Registers indicate what data has been extracted from a specified
!   netCDF file.  The forcing data/profiles which are held in this
!   module (_nc) are allocated/deallocated by subroutines in this module.
!  
!   NetCDF data will have been converted into units suitable for the SCM
!   and interpolated on levels specified by the current SCM run.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!-----------------------------------------------------------------------------
! Declarations:

  INTEGER ::             &
    yeari_nc      = IMDI &! Initial year
  , monthi_nc     = IMDI &! Initial month
  , dayi_nc       = IMDI &! Initial day
  , houri_nc      = IMDI &! Initial hour
  , mini_nc       = IMDI &! Initial minute
  , seci_nc       = IMDI &! Initial second
  , soil_type_nc  = IMDI  ! Soil type

  REAL ::                &
    grdbx_area_nc = RMDI &! Gridbox area
  , lat_nc        = RMDI &! Latitude
  , long_nc       = RMDI &! Longitude
  , albsoil_nc    = RMDI &! Soil albedo
  , orog_nc       = RMDI &! Orographic height (m)
  , obs_pd_nc     = RMDI &! Observation period (s)
  , obs_top_nc    = RMDI &! Highest obs height (m)
  , obs_bot_nc    = RMDI &! Lowest  obs height (m)
  , z0h_nc        = RMDI &! Surface roughness length (heat)
  , z0m_nc        = RMDI &! Surface roughness length (momentum)
  , sndep_nc      = RMDI &! Initial snow depth (m)
  , tstari_nc     = RMDI &! Initial surface temperature (K)
  , lsm_nc        = RMDI  ! Land-sea mask (fraction)


  ! Logicals to register what obs data has been
  ! extracted from netCDF file. Logicals are for ...
  LOGICAL ::                 &
    l_qi           = .FALSE. &!... init spec. humid profile
  , l_pi           = .FALSE. &!... init pressure profile
  , l_ui           = .FALSE. &!... init zonal wind
  , l_vi           = .FALSE. &!... init merid wind
  , l_wi           = .FALSE. &!... init vert wind
  , l_t_bg         = .FALSE. &!... background temperature
  , l_q_bg         = .FALSE. &!... background spec. humid
  , l_u_bg         = .FALSE. &!... background zonal wind
  , l_v_bg         = .FALSE. &!... background merid wind
  , l_w_bg         = .FALSE. &!... background vert wind
  , l_u_inc        = .FALSE. &!... zonal wind forcing incs
  , l_v_inc        = .FALSE. &!... merid wind forcing incs
  , l_w_inc        = .FALSE. &!... vert  wind forcing incs
  , l_t_inc        = .FALSE. &!... temperature forcing incs
  , l_q_inc        = .FALSE. &!... spec humid forcing incs
  , l_ozone        = .FALSE. &!... ozone profile
  , l_theta        = .FALSE. &!... init potential temperature
  , l_shflx        = .FALSE. &!... surface sensible heat flux
  , l_lhflx        = .FALSE. &!... surface latent heat flux
  , l_tstar_for    = .FALSE. &!... surface temperture
  , l_smci         = .FALSE. &!... init soil moisture content
  , l_smcli        = .FALSE. &!... init soil moisture content in layers
  , l_frac_typ     = .FALSE. &!... fraction of each tile type
  , l_canopy       = .FALSE. &!... surface canopy water content
  , l_z0_tile      = .FALSE. &!... surface tile roughness lengths
  , l_z0h_tile     = .FALSE. &!... surface tile thermal roughness lengths
  , l_tstar_tile   = .FALSE. &!... surface tile temperatures
  , l_t_deep_soili = .FALSE.  !... soil temperature profile


  ! Rank 1 allocatable arrays
  REAL, ALLOCATABLE ::    &
    canopy_gbi_nc   (:)   &! Init gridbox mean canopy water content (kg/m2)
  , smci_nc         (:)    ! Init soil moisture content (kg/m2)

  ! Rank 2 allocatable arrays
  REAL, ALLOCATABLE ::    &
    smcli_nc        (:,:) &! Init soil moisture content in layers
  , t_deep_soili_nc (:,:) &! Deep soil temperature profile
  , canopy_nc       (:,:) &! Canopy water content on tiles
  , z0_tile_nc      (:,:) &! Tile surface roughness lengths
  , z0h_tile_nc     (:,:) &! Tile surface roughness lengths
  , frac_typ_nc     (:,:) &! Fraction of each tile type
  , tstar_tile_nc   (:,:)  ! Surface tile temperature

  ! Rank 3 allocatable arrays
  REAL, ALLOCATABLE :: &
    tstar_nc (:,:,:)   &! Surface temperature timeseries
  , lhflx_nc (:,:,:)   &! Latent heat flux timeseries
  , shflx_nc (:,:,:)   &! Sensible heat flux timeseries
  , o3_nc    (:,:,:)   &! Ozone profile  
  , thi_nc   (:,:,:)   &! Potential temperature profile
  , ti_rh_nc (:,:,:)   &! Init. temperature on rho-grid
  , ti_th_nc (:,:,:)   &! Init. temperature on theta-grid
  , pi_rh_nc (:,:,:)   &! Init. pressure on rho-grid
  , pi_th_nc (:,:,:)   &! Init. pressure on theta-grid
  , qi_nc    (:,:,:)   &! Init. specific humidity profile
  , ui_nc    (:,:,:)   &! Init. zonal wind profile
  , vi_nc    (:,:,:)   &! Init. merid wind profile
  , wi_nc    (:,:,:)    ! Init. vert. wind profile

  ! Rank 4 allocatable arrays
  ! Obs. profiles of ...
  REAL, ALLOCATABLE ::   &
    t_bg_nc    (:,:,:,:) &!... background temperature
  , q_bg_nc    (:,:,:,:) &!... background specific humidity
  , u_bg_nc    (:,:,:,:) &!... background zonal wind 
  , v_bg_nc    (:,:,:,:) &!... background merid wind 
  , w_bg_nc    (:,:,:,:) &!... background vert. wind 
  , t_h_inc_nc (:,:,:,:) &!... temperature inc due to horiz adv
  , t_v_inc_nc (:,:,:,:) &!... temperature inc due to vert  adv
  , q_h_inc_nc (:,:,:,:) &!... spec. humid inc due to horiz adv
  , q_v_inc_nc (:,:,:,:) &!... spec. humid inc due to vert  adv
  , u_h_inc_nc (:,:,:,:) &!... zonal wind  inc due to horiz adv
  , u_v_inc_nc (:,:,:,:) &!... zonal wind  inc due to vert  adv
  , v_h_inc_nc (:,:,:,:) &!... merid wind  inc due to horiz adv
  , v_v_inc_nc (:,:,:,:) &!... merid wind  inc due to vert  adv
  , w_inc_nc   (:,:,:,:)  !... vert  wind  inc

CONTAINS
!-----------------------------------------------------------------------

  SUBROUTINE alloc_nc_arr

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('ALLOC_NC_ARR',zhook_in,zhook_handle)

    ! Land Points
    IF (nlnd > 0) THEN
      ALLOCATE                              &
        ( smci_nc         (nlnd)            &
        , canopy_gbi_nc   (nlnd)            &
        , smcli_nc        (nlnd, sm_lv)     &
        , t_deep_soili_nc (nlnd, st_lv)     &
        , canopy_nc       (nlnd, ntile)     &
        , z0_tile_nc      (nlnd, ntile)     &
        , z0h_tile_nc     (nlnd, ntile)     &
        , frac_typ_nc     (nlnd, ntype)     &
        , tstar_tile_nc   (nlnd, ntile) )
    END IF

    ! Initial Profiles
    ALLOCATE                                &
      ( pi_rh_nc (rw_lng, rw,   nmod_lv+1)  &
      , pi_th_nc (rw_lng, rw,   nmod_lv  )  &
      , thi_nc   (rw_lng, rw,   nmod_lv  )  &  
      , ti_rh_nc (rw_lng, rw,   nmod_lv+1)  &  
      , ti_th_nc (rw_lng, rw,   nmod_lv  )  &   
      , qi_nc    (rw_lng, rw,   nwet_lv  )  &
      , o3_nc    (rw_lng, rw,   nmod_lv  )  &
      , ui_nc    (rw_lng, rw,   nmod_lv  )  &
      , vi_nc    (rw_lng, rw,   nmod_lv  )  &
      , wi_nc    (rw_lng, rw, 0:nmod_lv  ) )

    ! Obs forcing
    ALLOCATE                                      &
      ( tstar_nc   (rw_lng, rw, nobs)             &
      , lhflx_nc   (rw_lng, rw, nobs)             &
      , shflx_nc   (rw_lng, rw, nobs) )

    ALLOCATE                                      &
      ( t_h_inc_nc (rw_lng, rw, nobs,   nmod_lv)  &
      , t_v_inc_nc (rw_lng, rw, nobs,   nmod_lv)  &
      , u_h_inc_nc (rw_lng, rw, nobs,   nmod_lv)  &
      , u_v_inc_nc (rw_lng, rw, nobs,   nmod_lv)  &
      , v_h_inc_nc (rw_lng, rw, nobs,   nmod_lv)  &
      , v_v_inc_nc (rw_lng, rw, nobs,   nmod_lv)  &
      , q_h_inc_nc (rw_lng, rw, nobs,   nwet_lv)  &
      , q_v_inc_nc (rw_lng, rw, nobs,   nwet_lv)  &
      , w_inc_nc   (rw_lng, rw, nobs, 0:nmod_lv) )

    ! Background obs. fields
    ALLOCATE                                      &
      ( t_bg_nc    (rw_lng, rw, nobs,   nmod_lv)  &
      , u_bg_nc    (rw_lng, rw, nobs,   nmod_lv)  &
      , v_bg_nc    (rw_lng, rw, nobs,   nmod_lv)  &
      , q_bg_nc    (rw_lng, rw, nobs,   nwet_lv)  &
      , w_bg_nc    (rw_lng, rw, nobs, 0:nmod_lv) )


    ! Initialise with Missing data indicators
    pi_rh_nc   = RMDI
    pi_th_nc   = RMDI
    ti_rh_nc   = RMDI
    ti_th_nc   = RMDI
    thi_nc     = RMDI
    qi_nc      = RMDI
    o3_nc      = RMDI
    ui_nc      = RMDI
    vi_nc      = RMDI
    wi_nc      = RMDI
    tstar_nc   = RMDI
    lhflx_nc   = RMDI
    shflx_nc   = RMDI
    t_h_inc_nc = RMDI
    t_v_inc_nc = RMDI
    q_h_inc_nc = RMDI
    q_v_inc_nc = RMDI
    u_h_inc_nc = RMDI 
    u_v_inc_nc = RMDI
    v_h_inc_nc = RMDI
    v_v_inc_nc = RMDI  
    w_inc_nc   = RMDI
    t_bg_nc    = RMDI
    q_bg_nc    = RMDI
    u_bg_nc    = RMDI
    v_bg_nc    = RMDI
    w_bg_nc    = RMDI

    IF (nlnd > 0) THEN
      smci_nc         = RMDI  
      smcli_nc        = RMDI
      canopy_gbi_nc   = RMDI
      t_deep_soili_nc = RMDI
      canopy_nc       = RMDI
      z0_tile_nc      = RMDI
      z0h_tile_nc     = RMDI
      tstar_tile_nc   = RMDI
      frac_typ_nc     = RMDI
    END IF

    IF (lhook) CALL dr_hook('ALLOC_NC_ARR',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE alloc_nc_arr

!-----------------------------------------------------------------------

  SUBROUTINE dealloc_nc_arr

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('DEALLOC_NC_ARR',zhook_in,zhook_handle)

    ! Land points
    IF (nlnd > 0) THEN
      DEALLOCATE                                                       &
        ( smci_nc,    canopy_gbi_nc,  smcli_nc,     t_deep_soili_nc    &
        , canopy_nc,  z0_tile_nc,     z0h_tile_nc,  frac_typ_nc        &
        ,  tstar_tile_nc )
    END IF

    ! Initial Profiles
    DEALLOCATE                                                         &
      ( pi_rh_nc,  pi_th_nc,  thi_nc,  ti_rh_nc,  ti_th_nc             &
      , qi_nc,     o3_nc,     ui_nc,   vi_nc,     wi_nc )

    ! Obs forcing
    DEALLOCATE                                                         &
      ( tstar_nc,    lhflx_nc,    shflx_nc                             &
      , t_h_inc_nc,  t_v_inc_nc,  u_h_inc_nc,  u_v_inc_nc              &
      , v_h_inc_nc,  v_v_inc_nc,  q_h_inc_nc,  q_v_inc_nc,  w_inc_nc )

    ! Background Obs fields
    DEALLOCATE                                                         &
      ( t_bg_nc,  u_bg_nc,  v_bg_nc,  q_bg_nc,  w_bg_nc )


    IF (lhook) CALL dr_hook('DEALLOC_NC_ARR',zhook_out,zhook_handle)

    RETURN  
  END SUBROUTINE dealloc_nc_arr

!=============================================================================
END MODULE netCDF_arr
