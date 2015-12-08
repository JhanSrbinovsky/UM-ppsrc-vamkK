! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Calculate and output some SCM diagnostics

SUBROUTINE dgnstcs_scm_main                                                   &
  ( row_length, rows, land_points, model_levels, wet_levels, tr_levels        &
  , tr_vars, sm_levels                                                        &
  , st_levels, ntype, rho, timestep, u, v, T                                  &
  , theta, q, qcl, qcf, layer_cloud, cca, ccw, t_deep_soil, p_star, tstar     &
  , smc, canopy_gb, snodep, zh, z0msea, smcl, sthu, sthf, gs, lw_incs         &
  , photosynth_act_rad, tstar_tile, aerosol, free_tracers                     &
  , p_theta_levels, p                                                         &
  , iccb, icct, w, w_adv, area_cloud_fraction, bulk_cloud_fraction            &
  , cloud_fraction_liquid, cloud_fraction_frozen, cclwp, nSCMDpkgs            &
  , L_SCMDiags )

  USE scm_utils, ONLY:                                                        &
      zhook_in, zhook_out, jprb, lhook, dr_hook

  USE cv_run_mod, ONLY: l_ccrad, l_3d_cca

  USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels
  
  IMPLICIT NONE

! Description:
!   Essentially just a lot of calls to SCMoutput

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
!   Language: Fortran 90

! All parameters are Intent In

! Variables passed in through the parameter list...
  INTEGER ::      &
    row_length    &
  , rows          &
  , land_points   &
  , model_levels  &
  , wet_levels    &
  , tr_levels     &
  , tr_vars       &
  , sm_levels     &
  , st_levels     &
  , ntype


  REAL ::                                   &
    rho(row_length,rows,model_levels)       &! Density*r^2 (kg/m)
  , u(row_length,rows,model_levels)         &! Zonal wind (m/s)
  , v(row_length,rows,model_levels)         &! Meridional wind (m/s)
  , w(row_length,rows,0:model_levels)       &! Vertical velocity (m/s)
  , w_adv(row_length,rows,0:model_levels)   &! Advective vertical velocity
                                             ! (m/s)
  , t(row_length,rows,model_levels)         &! Temperature(K)
  , theta(row_length,rows,model_levels)     &! Potential temperature (K)
  , q(row_length,rows,wet_levels)           &! Specific humidity (kg/kg)
  , qcf(row_length,rows,wet_levels)         &! Cloud ice content (kg/kg)
  , qcl(row_length,rows,wet_levels)         &! Cloud water content(kg/kg)
  , layer_cloud(row_length,rows,wet_levels) &! Layer cloud amount
                                             ! (decimal fraction)
  , cca(row_length,rows,model_levels)       &! Convective cloud amount
  , ccw(row_length,rows,wet_levels)          ! Convective Cloud Water passed
                                             ! to radiation scheme (kg/kg)

  REAL ::                                             &
    area_cloud_fraction(row_length,rows,wet_levels)   &! Area cloud amount
                                                       ! (decimal fraction)
  , bulk_cloud_fraction(row_length,rows,wet_levels)   &! Cloud amount
                                                       ! (decimal fraction)
  , cloud_fraction_liquid(row_length,rows,wet_levels) &! Liquid cloud amount
                                                       ! (decimal fraction)
  , cloud_fraction_frozen(row_length,rows,wet_levels) &! Frozen cloud amount
                                                       ! (decimal fraction)

  , cclwp(row_length,rows)             &! condensed water path (kg/m^2)
  , t_deep_soil(land_points,st_levels) &! Deep soil temperatures (K)
  , p_star(row_length,rows)            &! Pressure at earth's surface
  , tstar(row_length,rows)             &! Surface temperature (K)
  , smc(land_points)                   &! Soil moisture content(kg/m^2)
  , smcl(land_points,sm_levels)        &! Soil moisture content in layers
                                        ! (kg/m^2)
  , sthf(land_points,sm_levels)        &! Frozen soil moisture content of
                                        ! each layer as a fraction of
                                        ! saturation.
  , canopy_gb(land_points)             &! Canopy water content (kg/m^2)
  , snodep(row_length,rows)            &! Snow depth (kg/m^2)
  , p(row_length,rows,model_levels+1)   ! Pressure on rho levels

  REAL ::                                         &
    p_theta_levels(row_length,rows,model_levels)  &! Pressure on theta levels
  , zh(row_length,rows)                           &! Height above surface of top
                                                   ! of boundary layer (m)
  , z0msea(row_length,rows)            &! Sea surface roughness length
  , sthu(land_points,sm_levels)        &! Unfrozen soil moisture content
                                        ! of each layer as a fraction of
                                        ! saturation. (kg/m^2)
  , gs(row_length*rows)                &! Stomatal conductance
  , LW_incs(row_length,rows,0:model_levels) &
  , photosynth_act_rad(row_length,rows)     &! Net downward shortwave
                                             ! radiation in band 1 (w/m^2).
  , tstar_tile(row_length*rows,ntype)       &! Surface tile temperature
  , aerosol(row_length,rows,wet_levels)     &! Aerosol values ; only used if
                                             ! l_murk=.true. ; default .false.
  , free_tracers(row_length,rows,tr_levels,tr_vars)    &
  , timestep

  INTEGER ::                  &
    iccb(row_length, rows)    &! Convective cloud base and top
  , icct(row_length, rows)     ! at levels 1 to model_levels

  INTEGER :: nSCMDpkgs      ! No of SCM diagnostics packages
  LOGICAL :: L_SCMDiags(nSCMDpkgs)  ! Diagnostics packages logicals

! Include parameters necessary for calls to SCMoutput...
! Start of include file: s_scmop.h
! Description:
!  Declares and defines some parameters necessary for calling SCMoutput
!
!

! Integers to represent the different time profiles. All must
! be non-negative and less than "only_radsteps".

  INTEGER, PARAMETER :: &
    t_inst        = 1   &! Give the instantaneous value
  , t_avg         = 2   &! Construct the average value
  , t_max         = 3   &! " maximum value
  , t_min         = 4   &! " minimum value
  , t_acc         = 5   &! " accumulated value
  , t_div         = 7   &! " average value divided
                         !   by another diagnostic
  , t_mult        = 8   &! " average value multiplied
                         !   by another diagnostic
  , t_acc_div     = 9   &! " accumulated value divided
                         !   by another diagnostic
  , t_acc_mult    = 10  &! " accumulated value multiplied
                         !   by another diagnostic
  , t_const       = 11  &! The value is constant.
  , only_radsteps = 100  ! When added to one of the above parameters,
                         ! flags that the diagnostic is only available
                         ! on radiation timesteps

! Integers to represent the different domain profiles
  INTEGER, PARAMETER :: &
    d_sl      = 1       &
  , d_soilt   = 2       &
  , d_bl      = 3       &
  , d_wet     = 4       &
  , d_all     = 5       &
  , d_soilm   = 6       &
  , d_tile    = 7       &
  , d_vis     = 9       &
  , d_point   = 13      &
  , d_allxtra = 14      &
  , d_land    = 15      &
  , d_cloud   = 16

! Statement function to encode a stream number into an integer
  INTEGER :: &
    Stream   &
  , strm

  Stream(strm) = 2**(strm-1)

! The default streams for diagnostics to go to will be 1,2,3,4,5 and 6.
! The following should thus be equal to:
!
! Stream(1) [2^0=1] + Stream(2) [2^1=2]  + Stream(3) [2^2=4]
! Stream(4) [2^3=8] + Stream(5) [2^4=16] + Stream(6) [2^5=32]
! Total = 63
!
! where Stream() is the statement function defined above.
! default is 63 (all)

  INTEGER, PARAMETER :: &
    default_streams = 63

! Integers to represent the different diagnostics packages
  INTEGER, PARAMETER :: &
    SCMDiag_gen   = 1   & ! General diagnostics
  , SCMDiag_rad   = 2   & ! Radiation
  , SCMDiag_bl    = 3   & ! Boundary layer
  , SCMDiag_surf  = 4   & ! Surface
  , SCMDiag_land  = 5   & ! Land points only
  , SCMDiag_sea   = 6   & ! Sea points only
  , SCMDiag_lsp   = 7   & ! Large scale precip
  , SCMDiag_conv  = 8   & ! Convection
  , SCMDiag_lscld = 9   & ! Large scale cloud
  , SCMDiag_pc2   = 10  & ! PC2
  , SCMDiag_forc  = 11  & ! Forcing
  , SCMDiag_incs  = 12  & ! Increments
  , SCMDiag_gwd   = 13    ! Gravity Wave Drag

! End of include file: s_scmop.h

! Local variables...

  CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'Dgnstcs_Scm_Main'
  CHARACTER (LEN=2)        :: cdiv    ! Character string for tracer counter

  REAL ::                            &
    sum_p_col(row_length,rows)       &! Sums of pressures on theta levels
  , coltemp(row_length,rows)         &! Average pressure-weighted
                                      ! column temperature
  , colq(row_length,rows)            &! Average pressure-weighted
                                      ! column humidity
  , accum_pptn(row_length,rows)      &! Accumulated precipitation
  , pptn_rate(row_length,rows)       &! Precipitation rate
  , ccb_m(row_length,rows)           &! Convective cloud base (m)
  , ccb_Pa(row_length,rows)          &! Convective cloud base (Pa)
  , cct_m(row_length,rows)           &! Convective cloud top  (m)
  , cct_Pa(row_length,rows)          &! Convective cloud top  (Pa)
  , rh(row_length,rows,wet_levels,2) &! Relative humidity (over land and water)
  , qst                               ! A saturation mixing ratio


  REAL ::                                  &
    z_theta(row_length,rows,model_levels)  &! height above surface (m) (Theta)
  , z_rho(row_length,rows,model_levels)    &! height above surface (m) (Rho)
  , rho_only(row_length,rows,model_levels)  ! Actual density (kg/m3)

  REAL :: a2out(row_length,rows,model_levels)
  INTEGER :: i,j,k,n

  ! Dr Hook
  !==============================
  REAL(KIND=jprb) :: zhook_handle

  IF (lhook) CALL dr_hook('DGNSTCS_SCM_MAIN',zhook_in,zhook_handle)

! Store diagnostics...

!-----------------------------------------------------------------------
!     SCM General Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_gen)) THEN

! Stash 0,002
! DEPENDS ON: scmoutput
    CALL scmoutput(u, 'u'                                                     &
      , 'Zonal wind', 'm/s'                                                   &
      , t_avg, d_all, default_streams, '', RoutineName)

! Stash 0,003
! DEPENDS ON: scmoutput
    CALL scmoutput(v, 'v'                                                     &
      , 'Meridional wind', 'm/s'                                              &
      , t_avg, d_all, default_streams, '', RoutineName)

! Stash 30,111
! DEPENDS ON: scmoutput
    CALL scmoutput(T, 'T'                                                     &
      , 'Temperature', 'K'                                                    &
      , t_avg, d_all, default_streams, '', RoutineName)

! Stash 0,004
! DEPENDS ON: scmoutput
    CALL scmoutput(theta, 'theta'                                             &
      , 'Potential temperature', 'K'                                          &
      , t_avg, d_all, default_streams, '', RoutineName)

! Stash 0,010
! DEPENDS ON: scmoutput
    CALL scmoutput(q, 'q'                                                     &
      , 'Specific humidity', 'kg/kg'                                          &
      , t_avg, d_wet, default_streams, '', RoutineName)

    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          a2out(i,j,k) = w(i,j,k)
        END DO
      END DO
    END DO

! Stash 0,150
! DEPENDS ON: scmoutput
    CALL scmoutput(a2out, 'w'                                                 &
      , 'Vertical velocity', 'm/s'                                            &
      , t_avg, d_all, default_streams, '', RoutineName)

! Stash 0,407
! DEPENDS ON: scmoutput
    CALL scmoutput(p, 'p_rho'                                                 &
      , 'Pressure on rho levels', 'Pa', t_avg, d_all, default_streams, ''     &
      , RoutineName)

! Stash 0,408
! DEPENDS ON: scmoutput
    CALL scmoutput(p_theta_levels, 'p_theta'                                  &
      , 'Pressure on theta levels', 'Pa'                                      &
      , t_avg, d_all, default_streams, '', RoutineName)

    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
! Height above model surface as calculation in ni_conv_ctl
          z_theta(i,j,k)  = r_theta_levels(i,j,k) - r_theta_levels(i,j,0)
          z_rho(i,j,k)    = r_rho_levels(i,j,k)   - r_theta_levels(i,j,0)
          rho_only(i,j,k) = rho(i,j,k)                                        &
                          / (r_rho_levels(i,j,k) * r_rho_levels(i,j,k))
        END DO ! i
      END DO ! j
    END DO ! k

! Stash 15,101
! DEPENDS ON: scmoutput
    CALL scmoutput(z_theta, 'h_theta'                                         &
      , 'Height of model theta levels', 'm'                                   &
      , t_avg, d_all, default_streams, '', RoutineName)

! Stash 15,102
! DEPENDS ON: scmoutput
    CALL scmoutput(z_rho, 'h_rho'                                             &
      , 'Height of model rho levels', 'm'                                     &
      , t_avg, d_all, default_streams, '', RoutineName)

! Stash 0,253
! DEPENDS ON: scmoutput
    CALL scmoutput(rho, 'rho_r2'                                              &
      , 'Density *r*r after timestep', 'kg/m'                                 &
      , t_avg, d_all, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(rho_only, 'rho_only'                                       &
      , 'Density after timestep', 'kg/m3'                                     &
      , t_avg, d_all, default_streams, '', RoutineName)

! Calculate relative humidities on theta levels
    DO j=1, rows
      DO i=1, row_length
        DO k=1, wet_levels

! DEPENDS ON: qsat
          CALL qsat (qst, T(i,j,k), p_theta_levels(i,j,k), 1)
          rh(i,j,k,1) = q(i,j,k) / qst

! DEPENDS ON: qsat_wat
          CALL qsat_wat (qst, T(i,j,k), p_theta_levels(i,j,k), 1)
          rh(i,j,k,2) = q(i,j,k) / qst

        END DO ! k
      END DO ! i
    END DO ! j

! Stash 30,113
! DEPENDS ON: scmoutput
    CALL scmoutput(rh(1,1,1,1), 'rh'                                          &
      , 'Relative humidity', '%'                                              &
      , t_avg, d_wet, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(rh(1,1,1,2), 'rh2'                                         &
      , 'Relative humidity over liquid water', '%'                            &
      , t_avg, d_wet, default_streams, '', RoutineName)

! Calculate pressure-weighted average temperature and humidity
    i = 1
    j = 1
    colq(i,j) = 0.0
    coltemp(i,j) = 0.0
    sum_p_col(i,j) = 0.0
    DO k=1, model_levels
      sum_p_col(i,j) = sum_p_col(i,j) + p_theta_levels(i,j,k)
      coltemp(i,j)   = coltemp(i,j)   + T(i,j,k)*p_theta_levels(i,j,k)
      colq(i,j)      = colq(i,j)      + Q(i,j,k)*p_theta_levels(i,j,k)
    END DO ! k

! DEPENDS ON: scmoutput
    CALL scmoutput(sum_p_col, 'sum_p_col'                                     &
      , 'Sum of theta-level pressures', 'Pa'                                  &
      , t_acc, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(coltemp, 'tatmos'                                          &
      , 'Mean column temperature', 'K'                                        &
      , t_acc_div, d_sl, default_streams, 'sum_p_col', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(colq, 'qatmos'                                             &
      , 'Mean column specific humidity', 'kg/kg'                              &
      , t_acc_div, d_sl, default_streams, 'sum_p_col', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(aerosol, 'aerosol'                                         &
      , 'Murk aerosol concentration', 'micro-g/kg'                            &
      , t_avg, d_wet, default_streams, '', RoutineName)

    IF (tr_vars > 0) THEN

      DO n=1,tr_vars 

        a2out(:,:,:) = 0.0
        DO k=1, tr_levels
          a2out(:,:,k) = free_tracers(:,:,k,n)
        END DO
        WRITE(cdiv,'(I2.2)') n

        ! DEPENDS ON: scmoutput
        CALL SCMoutput(a2out,                                                 &
           'tracer'//cdiv,'Concentration of tracer '//cdiv,'kg/kg',           &
           t_avg,d_all,default_streams,'',RoutineName)

      END DO

    END IF  ! test on tr_vars

  END IF ! L_SCMDiags(SCMDiag_gen)

!-----------------------------------------------------------------------
!     SCM General OR Large Scale Cloud Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_lscld)) THEN

! Stash 0,254
! DEPENDS ON: scmoutput
    CALL scmoutput(qcl, 'qcl'                                                 &
      , 'QCL cloud water content', 'kg/kg'                                    &
      , t_avg, d_wet, default_streams, '', RoutineName)

! Stash 0,012
! DEPENDS ON: scmoutput
    CALL scmoutput(qcf, 'qcf'                                                 &
      , 'QCF cloud ice content', 'kg/kg'                                      &
      , t_avg, d_wet, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(layer_cloud, 'layer_cloud'                                 &
      , 'Layer cloud amount', 'Fraction'                                      &
      , t_avg, d_wet, default_streams, '', RoutineName)

! Stash 0, 265
! DEPENDS ON: scmoutput
    CALL scmoutput(area_cloud_fraction, 'acf'                                 &
      , 'Area cloud fraction', 'Fraction'                                     &
      , t_avg, d_wet, default_streams, '', RoutineName)

! Stash 0, 266
! DEPENDS ON: scmoutput
    CALL scmoutput(bulk_cloud_fraction, 'bcf'                                 &
      , 'Bulk cloud fraction', 'Fraction'                                     &
      , t_avg, d_wet, default_streams, '', RoutineName)

! Stash 0, 267
! DEPENDS ON: scmoutput
    CALL scmoutput(cloud_fraction_liquid, 'cfl'                               &
      , 'Liquid cloud fraction', 'Fraction'                                   &
      , t_avg, d_wet, default_streams, '', RoutineName)

! Stash 0, 268
! DEPENDS ON: scmoutput
    CALL scmoutput(cloud_fraction_frozen, 'cff'                               &
      , 'Frozen cloud fraction', 'Fraction'                                   &
      , t_avg, d_wet, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_lscld)

!-----------------------------------------------------------------------
!     SCM General OR Convection Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_conv)) THEN

! Stash 0,016
! Condensed water path convective cloud passed to radiation
! DEPENDS ON: scmoutput
    CALL scmoutput(cclwp, 'cclwp2rad'                                         &
      , 'cclwp passed to Radiation', 'kg/m^2'                                 &
      , t_avg, d_sl, default_streams, '', RoutineName)


    IF (l_3d_cca) THEN
! Stash 0,211
! Convective Cloud Amount passed to Radiation
! DEPENDS ON: scmoutput
      CALL SCMoutput(cca                                                      &
        , 'cca2rad','CCA passed to radiation', 'Fraction'                     &
        , t_avg, d_all, default_streams, '', RoutineName)
    ELSE
! Stash 0,211
! Convective Cloud Amount passed to Radiation
! DEPENDS ON: scmoutput
      CALL SCMoutput(cca(1,1,1)                                               &
        , 'cca2rad','CCA passed to radiation', 'Fraction'                     &
        , t_avg, d_sl, default_streams, '', RoutineName)
    END IF

! Stash 0,212
! Convective Cloud Water passed to Radiation
    IF (l_ccrad) THEN

!DEPENDS ON: scmoutput
      CALL scmoutput(ccw, 'ccw2rad'                                           &
        , 'CCW passed to radiation', 'kg/kg'                                  &
        , t_avg, d_wet, default_streams, '', RoutineName)
    END IF

    DO k=1,model_levels
      DO j=1,rows
        DO i=1,row_length
          a2out(i,j,k) = w_adv(i,j,k)
        END DO
      END DO
    END DO

! Stash 0,258
! DEPENDS ON: scmoutput
    CALL scmoutput(a2out, 'w_adv'                                             &
      , 'Advective vertical velocity', 'm/s'                                  &
      , t_avg, d_all, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_conv)

!-----------------------------------------------------------------------
!     SCM General OR Surface Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_surf)) THEN

! Stash 0,409
! DEPENDS ON: scmoutput
    CALL scmoutput(p_star, 'pstar'                                            &
      , 'Surface pressure', 'Pa'                                              &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 0,024
! DEPENDS ON: scmoutput
    CALL scmoutput(tstar, 'tstar'                                             &
      , 'Surface temperature', 'K'                                            &
      , t_avg, d_sl, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_surf)

!-----------------------------------------------------------------------
!     SCM General OR Radiation Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_rad)) THEN

! Stash 2,232
! DEPENDS ON: scmoutput
    CALL scmoutput(lw_incs, 'lwrate_day'                                      &
      , 'LW heating rate', 'K/day'                                            &
      , t_mult, d_allxtra, default_streams, 'ntspday', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(lw_incs, 'lwrate_step'                                     &
      , 'LW heating rate', 'K/timestep'                                       &
      , t_avg, d_allxtra, default_streams, 'ntspday', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(lw_incs, 'lwrate_acc'                                      &
      , 'Accumulated LW heating rate over dumping period', 'K/ts*dumps'       &
      , t_acc, d_allxtra, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_gen) .OR. L_SCMDiags(SCMDiag_rad)

!-----------------------------------------------------------------------
!     SCM Surface OR Radiation Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_surf) .OR. L_SCMDiags(SCMDiag_rad)) THEN

! surf_radflux has been removed since it was calculated only for
! sea ice and did not contain the full flux.

! DEPENDS ON: scmoutput
    CALL scmoutput(photosynth_act_rad, 'down_surf_sw_b1'                      &
      , 'Downward SW radn in band 1', 'W/m^2'                                 &
      , t_avg, d_sl, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_surf) .OR. L_SCMDiags(SCMDiag_rad)

!-----------------------------------------------------------------------
!     SCM Convection Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_conv)) THEN

! Convert ccb as a model level to ccb in meters and Pascals
    DO j=1, rows
      DO i=1, row_length
        IF (iccb(i,j) /= 0) THEN
          ccb_m(i,j)  = r_theta_levels(i,j,iccb(i,j)) - r_theta_levels(i,j,0)
          cct_m(i,j)  = r_theta_levels(i,j,icct(i,j)) - r_theta_levels(i,j,0)
          ccb_pa(i,j) = p_theta_levels(i,j,iccb(i,j))
          cct_pa(i,j) = p_theta_levels(i,j,icct(i,j))
        ELSE
          ccb_m(i,j)  = 0.0
          cct_m(i,j)  = 0.0
          ccb_pa(i,j) = 0.0
          cct_pa(i,j) = 0.0
        END IF
      END DO ! i
    END DO ! j

! DEPENDS ON: scmoutput
    CALL scmoutput                                                            &
       ( ccb_m*cca(1:row_length,1:rows, 1), 'ccb_z'                           &
       , 'Convective cloud base height (weighted average)', 'm'               &
       , t_div, d_sl, default_streams, 'cca2rad', RoutineName )

! DEPENDS ON: scmoutput
    CALL scmoutput                                                            &
       ( ccb_pa*cca(1:row_length,1:rows, 1), 'ccb_pa'                         &
       , 'Convective cloud base pressure (weighted average)', 'Pa'            &
       , t_div, d_sl, default_streams, 'cca2rad', RoutineName )

! DEPENDS ON: scmoutput
    CALL scmoutput                                                            &
       ( cct_m*cca(1:row_length,1:rows,1), 'cct_z'                            &
       , 'Convective cloud top height (weighted average)', 'm'                &
       , t_div, d_sl, default_streams, 'cca2rad', RoutineName )

! DEPENDS ON: scmoutput
    CALL scmoutput                                                            &
       ( cct_pa*cca(1:row_length,1:rows,1), 'cct_pa'                          &
       , 'Convective cloud top pressure (weighted average)', 'Pa'             &
       , t_div, d_sl, default_streams, 'cca2rad', RoutineName )

  END IF ! L_SCMDiags(SCMDiag_conv)

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_bl)) THEN

! DEPENDS ON: scmoutput
    CALL scmoutput(zh, 'bl_depth'                                             &
      , 'Boundary layer depth (zh)', 'm'                                      &
      , t_avg, d_sl, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_bl)

!-----------------------------------------------------------------------
!     SCM Sea Points Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_sea)) THEN

! DEPENDS ON: scmoutput
    CALL scmoutput(z0msea, 'sea_roughness'                                    &
      , 'Sea surface roughness length', 'm'                                   &
      , t_avg, d_sl, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_sea)

!-----------------------------------------------------------------------
!     SCM Land Points Diagnostics Package
!-----------------------------------------------------------------------

  IF (L_SCMDiags(SCMDiag_land) .AND. land_points > 0) THEN

! Stash 0,020
! DEPENDS ON: scmoutput
    CALL scmoutput(t_deep_soil, 'tsoildeep'                                   &
      , 'Deep soil temperatures', 'K'                                         &
      , t_avg, d_soilt, default_streams, '', RoutineName)

! Stash 8,208
! DEPENDS ON: scmoutput
    CALL scmoutput(smc, 'soilmoist'                                           &
      , 'Soil moisture content', 'kg/m^2'                                     &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 8,223
! DEPENDS ON: scmoutput
    CALL scmoutput(smcl, 'soilmoistlay'                                       &
      , 'Layer soil moisture content', 'kg/m^2'                               &
      , t_avg, d_soilm, default_streams, '', RoutineName)

! Stash 0,214
! DEPENDS ON: scmoutput
    CALL scmoutput(sthu, 'soilmoistunfroz'                                    &
      , 'Unfrozen soil moisture content', 'kg/m^2'                            &
      , t_avg, d_soilm, default_streams, '', RoutineName)

! Stash 0,215
! DEPENDS ON: scmoutput
    CALL scmoutput(sthf, 'soilmoistfroz'                                      &
      , 'Frozen soil moisture content', 'kg/m^2'                              &
      , t_avg, d_soilm, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(gs, 'stoma_cond'                                           &
      , 'Stomatal conductance', 'm/s'                                         &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 0,022
! DEPENDS ON: scmoutput
    CALL scmoutput(canopy_gb, 'canopy_gb'                                     &
      , 'Canopy water content', 'kg/m^2'                                      &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 0,023
! DEPENDS ON: scmoutput
    CALL scmoutput(snodep, 'snowdepth'                                        &
      , 'Snow depth', 'kg/m^2'                                                &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 3,316
! DEPENDS ON: scmoutput
    CALL scmoutput(tstar_tile, 'tstartl'                                      &
      , 'Tile surface temp', 'K'                                              &
      , t_avg, d_tile, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_land)/land_points

  IF (lhook) CALL dr_hook('DGNSTCS_SCM_MAIN',zhook_out,zhook_handle)

  RETURN
END SUBROUTINE dgnstcs_scm_main

