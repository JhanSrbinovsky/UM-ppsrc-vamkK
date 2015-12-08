! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module for SCM namelist : INOBSFOR

MODULE s_inobsfor

  USE scm_cntl_mod, ONLY: scm_nml

  USE s_maxdim,   ONLY:                                               &
    mx_rw_lng, mx_rw, mx_nlnd                                         &
  , mx_mod_lv, mx_wet_lv, mx_st_lv, mx_tr_lv, mx_tr_vars, mx_nobs

  USE scm_utils,  ONLY:                                               &
    rmdi, imdi, rw_lng, rw, nml_nmod_lv, nmod_lv, nwet_lv, nobs       &
  , interpolate, z_th, z_rh, nml_z_th, nml_z_rh                       &
  , zhook_in, zhook_out, jprb, lhook, dr_hook                         &
  , old_nml, old_rlx, old_vertadv

  USE s_interp_mod, ONLY: interp1d

  USE s_main_force, ONLY:                         &
    l_vertadv                                     &
  , obs_pd, obs_top, obs_bot                      &
  , rlx_t,  rlx_q,  rlx_u,  rlx_v,  rlx_w         &
  , plev_t, plev_q, plev_u, plev_v, plev_w        &
  , tau_t,  tau_q,  tau_u,  tau_v,  tau_w         &
  , scm_flux_h        => flux_h                   &
  , scm_flux_e        => flux_e                   &
  , scm_tstar_forcing => tstar_forcing            &
  , scm_q_star        => q_star                   &
  , scm_t_inc         => t_inc                    &
  , scm_u_inc         => u_inc                    &
  , scm_v_inc         => v_inc                    &
  , scm_w_inc         => w_inc                    &
  , scm_q_bg          => q_bg                     &
  , scm_t_bg          => t_bg                     &
  , scm_u_bg          => u_bg                     &
  , scm_v_bg          => v_bg                     &
  , scm_w_bg          => w_bg


  IMPLICIT NONE

!=============================================================================
!
! Description:
!   Allows SCM to read in INOBSFOR namelist from forcing file, scm_nml.
!   INOBSFOR contains data required if observational forcing is chosen
!   (obs logical)
!
! Method:
!   Namelist INOBSFOR is defined in this module and read in by contained
!   subroutine read_inobsfor.  Scalar variables are read directly to
!   s_main_force.  Array variables are declared as private single rank arrays 
!   of fixed length (arrays cannot be varying length in namelists) and
!   reshaped/interpolated before being transferred to arrays of the correct
!   size/shape in s_main_force.
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

  !--------------------------------------------------
  ! Fixed arrays to read in namelist before reshaping
  !--------------------------------------------------
  REAL, PRIVATE ::                                               &
    flux_h        (mx_rw_lng*mx_rw*mx_nobs)               = rmdi &
  , flux_e        (mx_rw_lng*mx_rw*mx_nobs)               = rmdi &
  , tstar_forcing (mx_rw_lng*mx_rw*mx_nobs)               = rmdi &
  , t_bg          (mx_rw_lng*mx_rw*mx_nobs*mx_mod_lv)     = rmdi &
  , u_bg          (mx_rw_lng*mx_rw*mx_nobs*mx_mod_lv)     = rmdi &
  , v_bg          (mx_rw_lng*mx_rw*mx_nobs*mx_mod_lv)     = rmdi &
  , q_bg          (mx_rw_lng*mx_rw*mx_nobs*mx_wet_lv)     = rmdi &
  , w_bg          (mx_rw_lng*mx_rw*mx_nobs*(mx_mod_lv+1)) = rmdi &
  , t_inc         (mx_rw_lng*mx_rw*mx_nobs*mx_mod_lv)     = rmdi &
  , u_inc         (mx_rw_lng*mx_rw*mx_nobs*mx_mod_lv)     = rmdi &
  , v_inc         (mx_rw_lng*mx_rw*mx_nobs*mx_mod_lv)     = rmdi &
  , q_star        (mx_rw_lng*mx_rw*mx_nobs*mx_wet_lv)     = rmdi &
  , w_inc         (mx_rw_lng*mx_rw*mx_nobs*(mx_mod_lv+1)) = rmdi


  !---------------------------------------------------------------------------
  ! Allocatable arrays:
  ! Used when intepolating namelist forcing profiles from namelist resolution
  ! to SCM resolution
  !---------------------------------------------------------------------------
  REAL, ALLOCATABLE ::     &
    nml_q_star  (:,:,:,:)  &! ... Spec. Humid  (kg/kg)/day
  , nml_t_inc   (:,:,:,:)  &! ... Temperature  K/day
  , nml_u_inc   (:,:,:,:)  &! ... Zonal wind   (m/s)/day
  , nml_v_inc   (:,:,:,:)  &! ... Merid wind   (m/s)/day
  , nml_w_inc   (:,:,:,:)  &! ... Vert. wind   (m/s)/day
  , nml_q_bg    (:,:,:,:)  &! ... Temperature  (K)
  , nml_t_bg    (:,:,:,:)  &! ... Spec. humid. (kg/kg)
  , nml_u_bg    (:,:,:,:)  &! ... Zonal wind   (m/s)
  , nml_v_bg    (:,:,:,:)  &! ... Merid wind   (m/s)
  , nml_w_bg    (:,:,:,:)   ! ... Vert. wind   (m/s)

  !---------------------------------------------------------------------------
  ! Define namelist
  !---------------------------------------------------------------------------
  NAMELIST/INOBSFOR/                                                          &
    l_vertadv, old_nml, old_rlx, old_vertadv                                  &
  , obs_pd, obs_top, obs_bot, rlx_t, rlx_q, rlx_u, rlx_v, rlx_w               &
  , plev_t, plev_q, plev_u, plev_v, plev_w, tau_t, tau_q, tau_u, tau_v, tau_w &
  , q_star, t_inc, u_inc, v_inc, w_inc, q_bg, t_bg, u_bg, v_bg, w_bg          &
  , flux_h, flux_e, tstar_forcing

  PRIVATE :: inobsfor

!=============================================================================
CONTAINS

   SUBROUTINE read_inobsfor

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    ! Local variables
    INTEGER :: i,j,k,kk
    INTEGER :: istatus
    INTEGER :: icode
    CHARACTER(LEN=13), PARAMETER :: routinename='read_inobsfor'

    IF (lhook) CALL dr_hook('READ_INOBSFOR',zhook_in,zhook_handle)

    OPEN(10, FILE=TRIM(ADJUSTL(scm_nml)), IOSTAT=istatus)

    IF (istatus /= 0) THEN
      icode = 500
      WRITE(6,*) ' Error opening file on unit 10 from '//routinename
      WRITE(6,*) ' Filename = '//TRIM(ADJUSTL(scm_nml))
      WRITE(6,*) ' IOstat =', istatus
    END IF

    READ  (10, INOBSFOR)
    CLOSE (10)


    IF (flux_e(1) /= rmdi)                                                    &
      scm_flux_e = RESHAPE( flux_e ,(/rw_lng,rw,nobs/))

    IF (flux_h(1) /= rmdi)                                                    &
      scm_flux_h = RESHAPE( flux_h ,(/rw_lng,rw,nobs/))

    IF (tstar_forcing(1) /= rmdi)                                             &
      scm_tstar_forcing = RESHAPE( tstar_forcing ,(/rw_lng,rw,nobs/))



    IF (interpolate) THEN

      CALL alloc_inobsfor
      CALL interp_inobsfor
      CALL dealloc_inobsfor

    ELSE
      ! Reshape INOBSFOR variables
      !===========================
      IF (nobs /= IMDI) THEN

        IF (t_inc(1)  /= rmdi)                                                &
          scm_t_inc  = RESHAPE( t_inc  ,(/rw_lng,rw,nobs,nmod_lv/) )

        IF (q_star(1) /= rmdi)                                                &
          scm_q_star = RESHAPE( q_star ,(/rw_lng,rw,nobs,nwet_lv/) )

        IF (u_inc(1)  /= rmdi)                                                &
          scm_u_inc  = RESHAPE( u_inc  ,(/rw_lng,rw,nobs,nmod_lv/) )

        IF (v_inc(1)  /= rmdi)                                                &
          scm_v_inc  = RESHAPE( v_inc  ,(/rw_lng,rw,nobs,nmod_lv/) )

        IF (w_inc(1)  /= rmdi)                                                &
          scm_w_inc  = RESHAPE( w_inc  ,(/rw_lng,rw,nobs,nmod_lv+1/) )




        IF (t_bg(1)   /= rmdi)                                                &
          scm_t_bg   = RESHAPE( t_bg   ,(/rw_lng,rw,nobs,nmod_lv/) )

        IF (q_bg(1)   /= rmdi)                                                &
          scm_q_bg   = RESHAPE( q_bg   ,(/rw_lng,rw,nobs,nwet_lv/) )

        IF (u_bg(1)   /= rmdi)                                                &
          scm_u_bg   = RESHAPE( u_bg   ,(/rw_lng,rw,nobs,nmod_lv/) )

        IF (v_bg(1)   /= rmdi)                                                &
          scm_v_bg   = RESHAPE( v_bg   ,(/rw_lng,rw,nobs,nmod_lv/) )

        IF (w_bg(1)   /= rmdi)                                                &
          scm_w_bg   = RESHAPE( w_bg   ,(/rw_lng,rw,nobs,nmod_lv+1/) )

      END IF        ! nobs
    END IF        ! interpolate

    IF (lhook) CALL dr_hook('READ_INOBSFOR',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE read_inobsfor

!-----------------------------------------------------------------------------

  SUBROUTINE interp_inobsfor

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    ! Local variables
    INTEGER :: i,j,k,kk

    ! Dummy arrays
    REAL, ALLOCATABLE :: &
      rdum1d  (:)        &
    , rdum1db (:)

    IF (lhook) CALL dr_hook('INTERP_INOBSFOR',zhook_in,zhook_handle)

    ! For obs, force wind to be zero at z=0
    ALLOCATE( rdum1d (0:nml_nmod_lv) &
            , rdum1db(0:nml_nmod_lv) )

    rdum1d(0)  = 0.0
    rdum1db(0) = 0.0


    !=========================================================================
    ! Forcing fields
    !=========================================================================

    !-------------------------------------------------------------------------
    ! q_star - Specific humidity forcing (kg/kg)/day
    !-------------------------------------------------------------------------
    IF (q_star(1) /= rmdi) THEN

      nml_q_star = RESHAPE( q_star, (/rw_lng,rw,nobs,nml_nmod_lv/) )

      DO kk=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            CALL interp1d                                                     &
              ( nml_z_th(:), nml_q_star(i,j,kk,:)                             &
              , z_th,        scm_q_star(i,j,kk,:) )
          END DO
        END DO
      END DO

      ! Assume q is constant from lowest obs level to level 1
      DO k=1, nmod_lv
        IF (z_th(k) <= nml_z_th(1))                                           &
          scm_q_star(:,:,:,k) = nml_q_star(:,:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! t_inc - Temperature forcing (K/day)
    !-------------------------------------------------------------------------
    IF (t_inc(1) /= rmdi) THEN

      nml_t_inc = RESHAPE( t_inc, (/rw_lng,rw,nobs,nml_nmod_lv/) )

      DO kk=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            CALL interp1d                                                     &
              ( nml_z_th(:), nml_t_inc(i,j,kk,:)                              &
              , z_th,        scm_t_inc(i,j,kk,:) )
          END DO
        END DO
      END DO

      ! Assume T is constant from lowest obs level to level 1
      DO k=1,nmod_lv
        IF (z_th(k) <= nml_z_th(1))                                           &
          scm_t_inc(:,:,:,k) = nml_t_inc(:,:,:,1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! u_inc - Zonal wind forcing (m/s)/day
    !-------------------------------------------------------------------------
    IF (u_inc(1) /= rmdi) THEN
      nml_u_inc = RESHAPE( u_inc, (/rw_lng,rw,nobs,nml_nmod_lv/) )

      DO kk=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            rdum1d(1:nml_nmod_lv)  = nml_z_rh(1:nml_nmod_lv)
            rdum1db(1:nml_nmod_lv) = nml_u_inc(i,j,kk,1:nml_nmod_lv)

            CALL interp1d                                                     &
              ( rdum1d,          rdum1db                                      &
              , z_rh(1:nmod_lv), scm_u_inc(i,j,kk,1:nmod_lv) )
          END DO
        END DO
      END DO

      ! Set missing data at top of profile to be equal to highest level
      ! which does contain data.
      DO k=2, nmod_lv
        WHERE ((scm_u_inc(:,:,:,k)   == rmdi) .AND.                           &
               (scm_u_inc(:,:,:,k-1) /= rmdi))                                &
          scm_u_inc(:,:,:,k) = scm_u_inc(:,:,:,k-1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! v_inc - Meridional wind forcing (m/s)/day
    !-------------------------------------------------------------------------
    IF (v_inc(1) /= rmdi) THEN
      nml_v_inc = RESHAPE( v_inc, (/rw_lng,rw,nobs,nml_nmod_lv/) )

      DO kk=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            rdum1d(1:nml_nmod_lv)  = nml_z_rh(1:nml_nmod_lv)
            rdum1db(1:nml_nmod_lv) = nml_v_inc(i,j,kk,1:nml_nmod_lv)

            CALL interp1d                                                     &
              ( rdum1d,          rdum1db                                      &
              , z_rh(1:nmod_lv), scm_v_inc(i,j,kk,1:nmod_lv) )
          END DO
        END DO
      END DO

      ! Set missing data at top of profile to be equal to highest level
      ! which does contain data.
      DO k=2, nmod_lv
        WHERE ((scm_v_inc(:,:,:,k)   == rmdi) .AND.                           &
               (scm_v_inc(:,:,:,k-1) /= rmdi))                                &
          scm_v_inc(:,:,:,k) = scm_v_inc(:,:,:,k-1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! w_inc - Vertical velocity forcing (m/s)/day
    !-------------------------------------------------------------------------
    IF (w_inc(1) /= rmdi) THEN
      nml_w_inc = RESHAPE( w_inc, (/rw_lng,rw,nobs,nml_nmod_lv+1/) )

      rdum1d(1:nml_nmod_lv) = nml_z_th(1:nml_nmod_lv)

      DO kk=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            CALL interp1d                                                     &
              ( rdum1d, nml_w_inc(i,j,kk,:)                                   &
              , z_th,   scm_w_inc(i,j,kk,1:nmod_lv) )
          END DO
        END DO
      END DO

      ! Set missing data at top of profile to be equal 0.0
      DO k=1, nmod_lv
        WHERE ((scm_w_inc(:,:,:,k)   == rmdi) .AND.                           &
               (scm_w_inc(:,:,:,k-1) /= rmdi))                                &
          scm_w_inc(:,:,:,k) = 0.0
      END DO
    END IF


    !=========================================================================
    ! Background fields
    !=========================================================================

    !-------------------------------------------------------------------------
    ! q_bg - Specific humidity (kg/kg)
    !-------------------------------------------------------------------------
    IF (q_bg(1) /= rmdi) THEN
      nml_q_bg = RESHAPE( q_bg, (/rw_lng,rw,nobs,nml_nmod_lv/) )

      DO kk=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            CALL interp1d                                                     &
              ( nml_z_th(:), nml_q_bg(i,j,kk,:)                               &
              , z_th,        scm_q_bg(i,j,kk,:) )
          END DO
        END DO
      END DO

    END IF


    !-------------------------------------------------------------------------
    ! t_bg - Temperature (K)
    !-------------------------------------------------------------------------
    IF (t_bg(1) /= rmdi) THEN
      nml_t_bg = RESHAPE( t_bg, (/rw_lng,rw,nobs,nml_nmod_lv/) )

      DO kk=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            CALL interp1d                                                     &
              ( nml_z_th(:), nml_t_bg(i,j,kk,:)                               &
              , z_th,        scm_t_bg(i,j,kk,:) )
          END DO
        END DO
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! u_bg - Zonal wind (m/s)
    !-------------------------------------------------------------------------
    IF (u_bg(1) /= rmdi) THEN
      nml_u_bg = RESHAPE( u_bg, (/rw_lng,rw,nobs,nml_nmod_lv/) )

      rdum1d(1:nml_nmod_lv)  = nml_z_rh(1:nml_nmod_lv)

      DO kk=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            rdum1db(1:nml_nmod_lv) = nml_u_bg(i,j,kk,1:nml_nmod_lv)

            CALL interp1d                                                     &
              ( rdum1d,          rdum1db                                      &
              , z_rh(1:nmod_lv), scm_u_bg(i,j,kk,1:nmod_lv) )
          END DO
        END DO
      END DO

      ! Set missing data at top of profile to be equal to highest level
      ! which does contain data.
      DO k=2, nmod_lv
        WHERE ((scm_u_bg(:,:,:,k)   == rmdi) .AND.                            &
               (scm_u_bg(:,:,:,k-1) /= rmdi))                                 &
          scm_u_bg(:,:,:,k) = scm_u_bg(:,:,:,k-1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! v_bg - Meridonal wind (m/s)
    !-------------------------------------------------------------------------
    IF (v_bg(1) /= rmdi) THEN
      nml_v_bg = RESHAPE( v_bg, (/rw_lng,rw,nobs,nml_nmod_lv/) )

      rdum1d(1:nml_nmod_lv)  = nml_z_rh(1:nml_nmod_lv)

      DO kk=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            rdum1db(1:nml_nmod_lv) = nml_v_bg(i,j,kk,1:nml_nmod_lv)

            CALL interp1d                                                     &
              ( rdum1d,          rdum1db                                      &
              , z_rh(1:nmod_lv), scm_v_bg(i,j,kk,1:nmod_lv) )
          END DO
        END DO
      END DO

      ! Set missing data at top of profile to be equal to highest level
      ! which does contain data.
      DO k=2, nmod_lv
        WHERE ((scm_v_bg(:,:,:,k)   == rmdi) .AND.                            &
               (scm_v_bg(:,:,:,k-1) /= rmdi))                                 &
          scm_v_bg(:,:,:,k) = scm_v_bg(:,:,:,k-1)
      END DO
    END IF


    !-------------------------------------------------------------------------
    ! w_bg - Vertical wind (m/s)
    !-------------------------------------------------------------------------
    IF (w_bg(1) /= rmdi) THEN
      nml_w_bg = RESHAPE( w_bg, (/rw_lng,rw,nobs,nml_nmod_lv+1/) )

      rdum1d(1:nml_nmod_lv) = nml_z_th(1:nml_nmod_lv)

      DO kk=1, nobs
        DO j=1, rw
          DO i=1, rw_lng
            CALL interp1d                                                     &
              ( rdum1d, nml_w_bg(i,j,kk,:)                                    &
              , z_th,   scm_w_bg(i,j,kk,1:nmod_lv) )
          END DO
        END DO
      END DO

      ! Set missing data at top of profile to be equal to highest level
      ! which does contain data.
      DO k=1, nmod_lv
        WHERE ((scm_w_bg(:,:,:,k)   == rmdi) .AND.                            &
               (scm_w_bg(:,:,:,k-1) /= rmdi))                                 &
          scm_w_bg(:,:,:,k) = scm_w_bg(:,:,:,k-1)
      END DO
    END IF

    DEALLOCATE( rdum1d, rdum1db )

    IF (lhook) CALL dr_hook('INTERP_INOBSFOR',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE interp_inobsfor

!-----------------------------------------------------------------------------

  SUBROUTINE alloc_inobsfor

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('ALLOC_INOBSFOR',zhook_in,zhook_handle)

    ALLOCATE                                        &
      ( nml_q_star (rw_lng,rw,nobs,  nml_nmod_lv)   &
      , nml_t_inc  (rw_lng,rw,nobs,  nml_nmod_lv)   &
      , nml_u_inc  (rw_lng,rw,nobs,  nml_nmod_lv)   &
      , nml_v_inc  (rw_lng,rw,nobs,  nml_nmod_lv)   &
      , nml_w_inc  (rw_lng,rw,nobs,0:nml_nmod_lv)   &

      , nml_q_bg   (rw_lng,rw,nobs,  nml_nmod_lv)   &
      , nml_t_bg   (rw_lng,rw,nobs,  nml_nmod_lv)   &
      , nml_u_bg   (rw_lng,rw,nobs,  nml_nmod_lv)   &
      , nml_v_bg   (rw_lng,rw,nobs,  nml_nmod_lv)   &
      , nml_w_bg   (rw_lng,rw,nobs,0:nml_nmod_lv) )

    IF (lhook) CALL dr_hook('ALLOC_INOBSFOR',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE alloc_inobsfor

!-----------------------------------------------------------------------------

  SUBROUTINE dealloc_inobsfor

    IMPLICIT NONE

    ! Dr Hook
    !==============================
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('DEALLOC_INOBSFOR',zhook_in,zhook_handle)

    DEALLOCATE        &
      (  nml_w_bg     &
      ,  nml_v_bg     &
      ,  nml_u_bg     &
      ,  nml_t_bg     &
      ,  nml_q_bg     &

      ,  nml_w_inc    &
      ,  nml_v_inc    &
      ,  nml_u_inc    &
      ,  nml_t_inc    &
      ,  nml_q_star )

    IF (lhook) CALL dr_hook('DEALLOC_INOBSFOR',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE dealloc_inobsfor

!=============================================================================
END MODULE s_inobsfor
