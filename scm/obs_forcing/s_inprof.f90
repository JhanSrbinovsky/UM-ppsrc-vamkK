! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : INPROF

MODULE s_inprof

  USE scm_cntl_mod, ONLY: scm_nml

  USE ancil_info,  ONLY: nsmax

  USE s_maxdim,  ONLY:                                      &
    mx_rw_lng, mx_rw, mx_nlnd, mx_ntile, mx_nsmax           &
  , mx_mod_lv, mx_wet_lv, mx_st_lv, mx_tr_lv, mx_tr_vars

  USE scm_utils, ONLY:                                      &
    rmdi, imdi, rw_lng, rw, nlnd, ntile, nml_nmod_lv        &
  , nmod_lv, nwet_lv, st_lv, ntr_lv, ntr_var                &
  , interpolate, z_th, z_rh, nml_z_th, nml_z_rh             &
  , zhook_in, zhook_out, jprb, lhook, dr_hook

  USE s_interp_mod, ONLY: interp1d

  USE s_main_force, ONLY:                   &
    kill_interp, obs_bot, obs_top           &
  , nml_inprof_thetal                       &
  , scm_iccti           => iccti            &
  , scm_iccbi           => iccbi            &
  , scm_di              => di               &
  , scm_u_0             => u_0              &
  , scm_v_0             => v_0              &
  , scm_ccai            => ccai             &
  , scm_snodepi         => snodepi          &
  , scm_tstari          => tstari           &
  , scm_z0mseai         => z0mseai          &
  , scm_z0m_scm         => z0m_scm          &
  , scm_z0h_scm         => z0h_scm          &
  , scm_smci            => smci             &
  , scm_canopy_gbi      => canopy_gbi       &
  , scm_sil_orog_land   => sil_orog_land    &
  , scm_ho2r2_orog      => ho2r2_orog       &
  , scm_ice_fract       => ice_fract        &
  , scm_t_deep_soili    => t_deep_soili     &
  , scm_theta           => theta            &
  , scm_ui              => ui               &
  , scm_vi              => vi               &
  , scm_qi              => qi               &
  , scm_wi              => wi               &
  , scm_w_advi          => w_advi           &
  , scm_p_in            => p_in             &
  , scm_i_snowdepth     => i_snowdepth      &
  , scm_i_snow_grnd     => i_snow_grnd      &
  , scm_i_rho_snow_grnd => i_rho_snow_grnd  &
  , scm_nsnow           => nsnow            &
  , scm_i_tsnowlayer    => i_tsnowlayer     &
  , scm_i_sice          => i_sice           &
  , scm_i_sliq          => i_sliq           &
  , scm_i_ds            => i_ds             &
  , scm_i_rgrainl       => i_rgrainl

  USE atmos_constants_mod,  ONLY: kappa, p_zero, recip_kappa

  IMPLICIT NONE

!=============================================================================
!
! Description:
!   Allows SCM to read in INPROF namelist from forcing file, scm_nml.
!   INPROF contains information relating to initial model primary variable
!   profiles (UMDP 1) plus T, but QCL and QCF are initialised by subroutine
!   INITQLCF conditions.
!
! Method:
!   Namelist INPROF is defined in this module and read in by contained
!   subroutine read_inprof.  Scalar variables are read directly to
!   s_main_force.  Array variables are declared as private single rank arrays
!   of fixed length (arrays cannot be varying length in namelists) and
!   reshaped/interpolated before being transferred to arrays of the
!   correct size/shape in s_main_force.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.
!
!=============================================================================
!
! Declarations:

  REAL :: z_tom_nml = rmdi

  !--------------------------------------------------
  ! Fixed arrays to read in namelist before reshaping
  !--------------------------------------------------

  INTEGER, PRIVATE ::                &
    iccti (mx_rw_lng*mx_rw)   = imdi &! Convective cloud top  (theta levs)
  , iccbi (mx_rw_lng*mx_rw)   = imdi  ! Convective cloud base (theta levs)

  REAL, PRIVATE ::                   &
    eta_th_nml (mx_mod_lv)    = rmdi &
  , eta_rh_nml (mx_mod_lv+1)  = rmdi

  REAL, PRIVATE ::                   &
    canopy_gbi (mx_nlnd)      = rmdi &
  , smci       (mx_nlnd)      = rmdi

  REAL, PRIVATE ::                   &
    di      (mx_rw_lng*mx_rw) = rmdi &! Equivalent thickness of sea-ice (m)
  , u_0     (mx_rw_lng*mx_rw) = rmdi &! Comp. of srf current, westerly (m/s)
  , v_0     (mx_rw_lng*mx_rw) = rmdi &! Comp. of srf current, southernly (m/s)
  , snodepi (mx_rw_lng*mx_rw) = rmdi &! Init. snow depth (kg/m2)
  , tstari  (mx_rw_lng*mx_rw) = rmdi  ! Init. surface temperature (K)

  REAL, PRIVATE ::                   &
    z0mseai (mx_rw_lng*mx_rw) = rmdi &! ... Momentum (sea-point) (m)
  , z0m_scm (mx_rw_lng*mx_rw) = rmdi &! ... Momentum (fixed) (m)
  , z0h_scm (mx_rw_lng*mx_rw) = rmdi  ! ... Heat (fixed) (m)

  REAL, PRIVATE ::                                            &
    sil_orog_land   (mx_rw_lng*mx_rw)  = rmdi                 &
                    ! Silhouette area of unresolved orography
                    ! per unit horizontal area (land points)

  , ho2r2_orog      (mx_rw_lng*mx_rw)  = rmdi                 &
                    ! SD of orography equivalent to peak to
                    ! trough height of unresolved orography
                    ! divided by 2SQRT(2) (land points) (m)

  , ice_fract       (mx_rw_lng*mx_rw)  = rmdi                 &
                    ! Gridbox fraction covered by sea ice

  , t_deep_soili    (mx_nlnd*mx_st_lv) = rmdi                 &
                    ! Init soil temp. profile (K)

  , i_snowdepth     (mx_nlnd*mx_ntile) = rmdi                 &
                    ! Init. snow depth (m)

  , i_snow_grnd     (mx_nlnd*mx_ntile) = rmdi                 &
                    ! Under canopy snow store (kg/m2)

  , i_rho_snow_grnd (mx_nlnd*mx_ntile) = rmdi                 &
                    ! Init. snow density ()

  , nsnow           (mx_nlnd*mx_ntile) = rmdi                 &
                    ! Number of snow levels (Real)

  , i_tsnowlayer    (mx_nlnd*mx_ntile*mx_nsmax) = rmdi        &
                    ! Init. snow layer temperature (K)

  , i_sice          (mx_nlnd*mx_ntile*mx_nsmax) = rmdi        &
                    ! Init. ice in snow pack (kg/m2)

  , i_sliq          (mx_nlnd*mx_ntile*mx_nsmax) = rmdi        &
                    ! Init. liquid in snow pack (kg/m2)

  , i_ds            (mx_nlnd*mx_ntile*mx_nsmax) = rmdi        &
                    ! Init. snow layer thickness (m)

  , i_rgrainl       (mx_nlnd*mx_ntile*mx_nsmax) = rmdi
                    ! Snow grain size in each layer (microns)


  REAL, PRIVATE ::                                 &
    theta  (mx_rw_lng*mx_rw*mx_mod_lv)      = rmdi &! Pot. temperature (K)
  , ui     (mx_rw_lng*mx_rw*mx_mod_lv)      = rmdi &! Zonal  wind (m/s)
  , vi     (mx_rw_lng*mx_rw*mx_mod_lv)      = rmdi &! Merid. wind (m/s)
  , qi     (mx_rw_lng*mx_rw*mx_wet_lv)      = rmdi &! Spec. humid (kg/kg)
  , wi     (mx_rw_lng*mx_rw*(mx_mod_lv+1))  = rmdi &! Vert.  wind (m/s)
  , w_advi (mx_rw_lng*mx_rw*(mx_mod_lv+1))  = rmdi &
  , p_in   (mx_rw_lng*mx_rw*(mx_mod_lv+1))  = rmdi &! Pressure (Pa)
  , ccai   (mx_rw_lng*mx_rw)                = rmdi


  !---------------------------------------------------------------------------
  ! Allocatable arrays:
  ! Used when intepolating namelist forcing profiles from namelist resolution
  ! to SCM resolution
  !---------------------------------------------------------------------------
  REAL,ALLOCATABLE ::     &
    nml_theta    (:,:,:)  &! Pot. temperature (K)
  , nml_ui       (:,:,:)  &! Zonal  wind (m/s)
  , nml_vi       (:,:,:)  &! Merid. wind (m/s)
  , nml_qi       (:,:,:)  &! Spec. humid (kg/kg)
  , nml_wi       (:,:,:)  &! Vert.  wind (m/s)
  , nml_w_advi   (:,:,:)  &
  , nml_p_in     (:,:,:)   ! Pressure (Pa

  REAL, ALLOCATABLE ::    &
    nml_ti_th    (:,:,:)  &
  , nml_exner_rh (:,:,:)  &
  , nml_exner_th (:,:,:)  &
  , exner_rh     (:,:,:)  &
  , exner_th     (:,:,:)  &
  , ti_rh        (:,:,:)  &
  , ti_th        (:,:,:)


  !---------------------------------------------------------------------------
  ! Define namelist
  !---------------------------------------------------------------------------
  NAMELIST/INPROF/                                                            &
    iccbi, iccti, di, u_0, v_0, ccai, snodepi, tstari, z0mseai, z0m_scm       &
  , z0h_scm, sil_orog_land, ho2r2_orog, ice_fract, t_deep_soili, canopy_gbi   &
  , qi, smci, theta, ui, vi, wi, p_in, w_advi, eta_th_nml                     &
  , eta_rh_nml, z_tom_nml, kill_interp, nml_inprof_thetal, i_snowdepth        &
  , i_snow_grnd, i_rho_snow_grnd, nsnow, i_tsnowlayer, i_sice, i_sliq, i_ds   &
  , i_rgrainl


  PRIVATE :: inprof

!=============================================================================
CONTAINS

  SUBROUTINE read_inprof

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    INTEGER :: i,j,k
    INTEGER :: istatus
    INTEGER :: icode
    CHARACTER(LEN=11), PARAMETER :: routinename='read_inprof'

    IF (lhook) CALL dr_hook('READ_INPROF',zhook_in,zhook_handle)

    OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus)

    IF (istatus /= 0) THEN
      icode = 500
      WRITE(6,*) ' Error opening file on unit 10 from '//routinename
      WRITE(6,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
      WRITE(6,*) ' IOstat =', istatus
    END IF

    READ  (10, inprof)
    CLOSE (10)

    IF (nml_nmod_lv /= IMDI) THEN

      nml_z_th = z_tom_nml*RESHAPE(eta_th_nml, (/nml_nmod_lv/))
      nml_z_rh = z_tom_nml*RESHAPE(eta_rh_nml, (/nml_nmod_lv+1/))

      ! Some profile data has been specified in namelist.scm
      ! Does The namelist need interpolating to UM levels ?
      !=====================================================
      IF (nml_nmod_lv /= nmod_lv) THEN
        interpolate = .TRUE.
      ELSE
        ! Allow user to override interpolation in case levels
        ! are flagged as different because of bit level differences
        ! Suppose it should really user specified., default .true.?
        IF (kill_interp) THEN
          interpolate = .FALSE.
        ELSE
          DO k=1, nmod_lv
            IF (z_th(k) /= nml_z_th(k)) interpolate = .TRUE.
            IF (z_rh(k) /= nml_z_th(k)) interpolate = .TRUE.
          END DO

          IF (z_rh(nmod_lv+1) /= nml_z_th(nmod_lv+1)) THEN
            interpolate = .TRUE.
          END IF
        END IF ! kill_interp
      END IF
    END IF

    ! By default, use namelist obs range for obs_bot/obs_top.
    ! This maybe overwritten later when reading in inobsfor
    ! obs_bot = nml_z_rh(1)
    ! obs_top = nml_z_rh(nml_nmod_lv+1)

    ! Single-level/surface fields
    IF (iccti(1) /= imdi) scm_iccti = RESHAPE( iccti, (/rw_lng,rw/) )
    IF (iccbi(1) /= imdi) scm_iccbi = RESHAPE( iccbi, (/rw_lng,rw/) )
    IF (di(1)    /= rmdi) scm_di    = RESHAPE( di,    (/rw_lng,rw/) )
    IF (u_0(1)   /= rmdi) scm_u_0   = RESHAPE( u_0,   (/rw_lng,rw/) )
    IF (v_0(1)   /= rmdi) scm_v_0   = RESHAPE( v_0,   (/rw_lng,rw/) )
    IF (ccai(1)  /= rmdi) scm_ccai  = RESHAPE( ccai,  (/rw_lng,rw/) )

    IF (snodepi(1) /= rmdi) scm_snodepi = RESHAPE( snodepi, (/rw_lng,rw/) )
    IF (tstari(1)  /= rmdi) scm_tstari  = RESHAPE( tstari,  (/rw_lng,rw/) )
    IF (z0mseai(1) /= rmdi) scm_z0mseai = RESHAPE( z0mseai, (/rw_lng,rw/) )
    IF (z0m_scm(1) /= rmdi) scm_z0m_scm = RESHAPE( z0m_scm, (/rw_lng,rw/) )
    IF (z0h_scm(1) /= rmdi) scm_z0h_scm = RESHAPE( z0h_scm, (/rw_lng,rw/) )

    IF (sil_orog_land(1) /= rmdi)                                             &
      scm_sil_orog_land = RESHAPE( sil_orog_land, (/rw_lng,rw/) )

    IF (ho2r2_orog(1) /= rmdi)                                                &
      scm_ho2r2_orog = RESHAPE( ho2r2_orog, (/rw_lng,rw/) )

    IF (ice_fract(1) /= rmdi)                                                 &
      scm_ice_fract = RESHAPE( ice_fract, (/rw_lng,rw/) )


    IF (nlnd > 0) THEN

      IF (canopy_gbi(1) /= rmdi)                                              &
        scm_canopy_gbi = RESHAPE( canopy_gbi, (/nlnd/) )

      IF (smci(1) /= rmdi)                                                    &
        scm_smci = RESHAPE( smci, (/nlnd/) )

      IF (t_deep_soili(1) /= rmdi)                                            &
        scm_t_deep_soili = RESHAPE( t_deep_soili, (/nlnd,st_lv/) )

      IF (i_snowdepth(1) /= rmdi)                                             &
        scm_i_snowdepth = RESHAPE( i_snowdepth, (/nlnd,ntile/) )

      IF (i_snow_grnd(1) /= rmdi)                                             &
        scm_i_snow_grnd = RESHAPE( i_snow_grnd, (/nlnd,ntile/) )

      IF (i_rho_snow_grnd(1) /= rmdi)                                         &
        scm_i_rho_snow_grnd = RESHAPE( i_rho_snow_grnd, (/nlnd,ntile/) )

      IF (nsnow(1) /= rmdi)                                                   &
        scm_nsnow = RESHAPE( nsnow, (/nlnd,ntile/) )

      IF (i_tsnowlayer(1) /= rmdi)                                            &
        scm_i_tsnowlayer = RESHAPE( i_tsnowlayer, (/nlnd,ntile,nsmax/) )

      IF (i_sice(1) /= rmdi)                                                  &
        scm_i_sice = RESHAPE( i_sice, (/nlnd,ntile,nsmax/) )

      IF (i_sliq(1) /= rmdi)                                                  &
        scm_i_sliq = RESHAPE( i_sliq, (/nlnd,ntile,nsmax/) )

      IF (i_ds(1) /= rmdi)                                                    &
        scm_i_ds = RESHAPE( i_ds, (/nlnd,ntile,nsmax/) )

      IF (i_rgrainl(1) /= rmdi)                                               &
        scm_i_rgrainl = RESHAPE( i_rgrainl, (/nlnd,ntile,nsmax/) )

    END IF


    ! Multi-level fields which may require interpolation
    IF (interpolate) THEN

      CALL alloc_inprof
      CALL interp_inprof
      CALL dealloc_inprof

    ELSE

      IF (theta(1)  /= rmdi)                                                  &
        scm_theta  = RESHAPE( theta,  (/rw_lng,rw,nmod_lv/) )

      IF (ui(1)     /= rmdi)                                                  &
        scm_ui     = RESHAPE( ui,     (/rw_lng,rw,nmod_lv/) )

      IF (vi(1)     /= rmdi)                                                  &
        scm_vi     = RESHAPE( vi,     (/rw_lng,rw,nmod_lv/) )

      IF (qi(1)     /= rmdi)                                                  &
        scm_qi     = RESHAPE( qi,     (/rw_lng,rw,nwet_lv/) )

      IF (wi(1)     /= rmdi)                                                  &
        scm_wi     = RESHAPE( wi,     (/rw_lng,rw,nmod_lv+1/) )

      IF (w_advi(1) /= rmdi)                                                  &
        scm_w_advi = RESHAPE( w_advi, (/rw_lng,rw,nmod_lv+1/) )

      IF (p_in(1)   /= rmdi)                                                  &
        scm_p_in   = RESHAPE( p_in,   (/rw_lng,rw,nmod_lv+1/) )

    END IF

    IF (lhook) CALL dr_hook('READ_INPROF',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE read_inprof

!=============================================================================

  SUBROUTINE interp_inprof

    IMPLICIT NONE

    ! Local variables
    INTEGER :: i,j,k

    ! Dummy arrays
    REAL, ALLOCATABLE :: &
      rdum1d  (:)        &
    , rdum1db (:)

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('INTERP_INPROF',zhook_in,zhook_handle)

    !-------------------------
    ! Initial pressure profile
    !-------------------------
    IF (p_in(1) /= rmdi) THEN
      nml_p_in = RESHAPE( p_in, (/rw_lng,rw,nml_nmod_lv+1/) )

      ! Calculate exner for namelist pressure levels

      DO k=1, nml_nmod_lv+1
        nml_exner_rh(:,:,k) = (nml_p_in(:,:,k)/p_zero)**(kappa)
      END DO

      ! Linearly interpolate for exner on UM theta levels
      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_rh(:), nml_exner_rh(i,j,:)                                &
            , nml_z_th(:), nml_exner_th(i,j,:) )

          CALL interp1d                                                       &
            ( nml_z_rh(:), nml_exner_rh(i,j,:)                                &
            , z_rh(:),     exner_rh(i,j,:) )

          CALL interp1d                                                       &
            ( nml_z_rh(:), nml_exner_rh(i,j,:)                                &
            , z_th(:),     exner_th(i,j,:) )
        END DO
      END DO

      ! Calc p_in_out/Theta_out from exner_rh/exner_th
      WHERE (exner_rh /= rmdi) scm_p_in = p_zero*exner_rh**recip_kappa

    END IF

    !-------------------------------------
    ! Initial theta profile
    !-------------------------------------
    IF ((theta(1) /= rmdi) .AND.                                              &
        (p_in(1)  /= rmdi)) THEN

      ! Calculate initial T profile
      nml_theta = RESHAPE(theta, (/rw_lng,rw,nml_nmod_lv/))

      ! Convert nml_theta to nml_ti
      nml_ti_th(:,:,:) = nml_theta(:,:,:) * nml_exner_th(:,:,:)

      ! Get T on rh/th levels
      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th(:), nml_ti_th(i,j,:)                                   &
            , z_th(:),     ti_th(i,j,:) )

          CALL interp1d                                                       &
            ( nml_z_th(:), nml_ti_th(i,j,:)                                   &
            , z_rh(:),     ti_rh(i,j,:) )
        END DO
      END DO

      WHERE ((exner_th /= rmdi) .AND. (ti_th /= rmdi))                        &
        scm_theta = ti_th / exner_th

    END IF


    !----------------------------------
    ! Initial specific humidity profile
    !----------------------------------
    IF (qi(1) /= rmdi) THEN
      nml_qi = RESHAPE( qi, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( nml_z_th(:), nml_qi(i,j,1:nml_nmod_lv)                          &
            , z_th,        scm_qi(i,j,1:nwet_lv) )
        END DO
      END DO

    END IF

    !--------------------------------
    ! Initial wind profiles, u, v, w
    !--------------------------------
    ! For obs, force wind to be zero at z=0
    ALLOCATE( rdum1d (0:nml_nmod_lv)                                          &
            , rdum1db(0:nml_nmod_lv) )

    rdum1d(0)  = 0.0
    rdum1db(0) = 0.0

    ! Initial u-wind
    IF (ui(1) /= rmdi) THEN
      nml_ui = RESHAPE( ui, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          rdum1d(1:nml_nmod_lv)  = nml_z_rh(1:nml_nmod_lv)
          rdum1db(1:nml_nmod_lv) = nml_ui(i,j,1:nml_nmod_lv)

          CALL interp1d                                                       &
            ( rdum1d,          rdum1db                                        &
            , z_rh(1:nmod_lv), scm_ui(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Set missing data at top of profile to be equal to highest level
      ! which does contain data.
      DO k=1, nmod_lv
        WHERE ((scm_ui(:,:,k)   == rmdi) .AND.                                &
               (scm_ui(:,:,k-1) /= rmdi))                                     &
          scm_ui(:,:,k) = scm_ui(:,:,k-1)
      END DO

    END IF


    ! Initial v-wind
    IF (vi(1) /= rmdi) THEN
      nml_vi = RESHAPE( vi, (/rw_lng,rw,nml_nmod_lv/) )

      DO j=1, rw
        DO i=1, rw_lng
          rdum1d(1:nml_nmod_lv)  = nml_z_rh(1:nml_nmod_lv)
          rdum1db(1:nml_nmod_lv) = nml_vi(i,j,1:nml_nmod_lv)

          CALL interp1d                                                       &
            ( rdum1d,          rdum1db                                        &
            , z_rh(1:nmod_lv), scm_vi(i,j,1:nmod_lv) )
        END DO
      END DO

      ! Set missing data at top of profile to be equal to highest level
      ! which does contain data.
      DO k=2, nmod_lv
        WHERE ((scm_vi(:,:,k)   == rmdi) .AND.                                &
               (scm_vi(:,:,k-1) /= rmdi) )                                    &
          scm_vi(:,:,k) = scm_vi(:,:,k-1)
      END DO

    END IF


    ! Initial w-wind
    IF ( wi(1) /= rmdi) THEN

      nml_wi = RESHAPE( wi, (/rw_lng,rw,nml_nmod_lv+1/) )

      rdum1d(1:nml_nmod_lv) = nml_z_th(1:nml_nmod_lv)
      rdum1d(0) = 0.0

      DO j=1, rw
        DO i=1, rw_lng
          CALL interp1d                                                       &
            ( rdum1d, nml_wi(i,j,:)                                           &
            , z_th,   scm_wi(i,j,1:nmod_lv))
        END DO
      END DO

      ! Set missing data at top of profile to be equal to 0.0
      DO k=2, nmod_lv
        WHERE ((scm_wi(:,:,k)   == rmdi) .AND.                                &
               (scm_wi(:,:,k-1) /= rmdi))                                     &
          scm_wi(:,:,k) = 0.0
      END DO

    END IF

    DEALLOCATE( rdum1d, rdum1db )

    IF (lhook) CALL dr_hook('INTERP_INPROF',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE interp_inprof

!=============================================================================

  SUBROUTINE alloc_inprof

    IMPLICIT NONE

    ! Dr Hook
    !=============================================================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('ALLOC_INPROF',zhook_in,zhook_handle)

    ALLOCATE                                    &
      ( nml_theta  (rw_lng,rw,  nml_nmod_lv)    &
      , nml_ui     (rw_lng,rw,  nml_nmod_lv)    &
      , nml_vi     (rw_lng,rw,  nml_nmod_lv)    &
      , nml_qi     (rw_lng,rw,  nml_nmod_lv)    &
      , nml_wi     (rw_lng,rw,0:nml_nmod_lv)    &
      , nml_w_advi (rw_lng,rw,0:nml_nmod_lv)    &
      , nml_p_in   (rw_lng,rw,  nml_nmod_lv+1) )

    ALLOCATE                                    &
      ( nml_ti_th    (rw_lng,rw,nml_nmod_lv)    &
      , nml_exner_rh (rw_lng,rw,nml_nmod_lv+1)  &
      , nml_exner_th (rw_lng,rw,nml_nmod_lv)    &
      , exner_rh     (rw_lng,rw,    nmod_lv+1)  &
      , exner_th     (rw_lng,rw,    nmod_lv)    &
      , ti_rh        (rw_lng,rw,    nmod_lv+1)  &
      , ti_th        (rw_lng,rw,    nmod_lv) )

    IF (lhook) CALL dr_hook('ALLOC_INPROF',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE alloc_inprof

!=============================================================================

  SUBROUTINE dealloc_inprof

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('DEALLOC_INPROF',zhook_in,zhook_handle)

    DEALLOCATE         &
      ( ti_rh          &
      , ti_th          &
      , exner_th       &
      , exner_rh       &
      , nml_exner_th   &
      , nml_exner_rh   &
      , nml_ti_th      )

    DEALLOCATE         &
      ( nml_p_in       &
      , nml_w_advi     &
      , nml_wi         &
      , nml_qi         &
      , nml_vi         &
      , nml_ui         &
      , nml_theta )

    IF (lhook) CALL dr_hook('DEALLOC_INPROF',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE dealloc_inprof

!=============================================================================
END MODULE s_inprof
