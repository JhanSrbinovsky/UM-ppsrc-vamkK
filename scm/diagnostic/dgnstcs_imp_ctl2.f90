! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Calculate and output some SCM diagnostics

SUBROUTINE dgnstcs_imp_ctl2                                                   &
  ! In
  ( row_length, rows, rhc_row_length, rhc_rows                                &
  , model_levels, wet_levels, bl_levels, cloud_levels, off_x, off_y           &
  , rhcrit, L_murk, combined_cloud, nclds, cumulus, ntml                      &
  , plsp, conv_rain, conv_snow, ls_rain, ls_snow, n_rows, R_u, R_v            &
  , u_incr_diagnostic, v_incr_diagnostic, zh, zht, bl_type_1, bl_type_2       &
  , bl_type_3, bl_type_4, bl_type_5, bl_type_6, bl_type_7, bl_alltypes        &
  , fqt, ftl, T, q, cf, cfl, cff, T_earliest, q_earliest, qcl_earliest        &
  , qcf_earliest, cf_earliest, cfl_earliest, cff_earliest, p_theta_levels     &
  , LWP, IWP, z0m, z0h_eff_gb, z0m_eff_gb, sea_ice_htf, sice_mlt_htf          &
  , taux, tauy, area_cloud_fraction, p_star, u10m, v10m                       &
  , latent_heat, cca_2d, rho1, qcl, qcf, surf_ht_flux_sice, rhcpt             &
  , surf_ht_flux_gb, aerosol, nSCMDpkgs, L_SCMDiags, BL_diag                  &

  ! InOut
  , q1p5m, t1p5m )

  USE atmos_constants_mod, ONLY: cp

  USE conversions_mod, ONLY: recip_pi_over_180, pi

  USE bl_diags_mod, ONLY:                                                     &
      strnewbldiag

  USE bl_option_mod, ONLY:                                                    &
      Fric_heating
  USE water_constants_mod, ONLY: lc
  
  USE vertnamelist_mod, ONLY:  first_constant_r_rho_level,                    &
      z_top_of_model, eta_theta, eta_rho, vertlevs

  USE level_heights_mod, ONLY: eta_theta_levels

  USE visbty_constants_mod, ONLY: n_vis_thresh, vis_thresh, calc_prob_of_vis 

  USE rad_pcf, ONLY: ip_cloud_mix_max, ip_cloud_mix_random

  USE beta_precip_mod, ONLY: beta_precip
  IMPLICIT NONE

! Description:
! Much of what is in here is (regrettably) duplication of code
! in diagnostics_lscld and diagnostics_bl. These routines can't
! be called in the SCM because of explicit STASH calls, but we
! still want the diagnostics.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

! Code Description:
! Language: Fortran77 with some bits of Fortran90

! ALL PARAMETERS ARE INTENT In APART FROM T1P5M AND Q1P5M
! WHICH ARE InOut.

  INTEGER ::          &
    row_length        &
  , rows              &
  , rhc_row_length    &
  , rhc_rows          &
  , model_levels      &
  , wet_levels        &
  , bl_levels         &
  , cloud_levels      &
  , off_x             &! MPP associated variables, expect both to be
  , off_y             &! zero since this is the SCM
  , nclds             &
  , n_rows

  LOGICAL ::                     &
    L_murk                       &! Switch for (visibility) aerosol
  , cumulus(row_length,rows)     &
  , ntml(row_length,rows)         ! Height of diagnosed BL top

  REAL ::                            &
    rhcrit(wet_levels)               &
  , plsp(row_length,rows)            &
  , conv_rain(row_length,rows)       &
  , conv_snow(row_length,rows)       &
  , ls_rain(row_length,rows)         &
  , ls_snow(row_length,rows)         &
  , R_u( 1-off_x:row_length+off_x    &
       , 1-off_y:rows+off_y          &
       , model_levels)               &
  , R_v( 1-off_x:row_length+off_x    &
       , 1-off_y:n_rows+off_y        &
       , model_levels)

  REAL ::                                               &
    u_incr_diagnostic(row_length, rows, model_levels)   &
  , v_incr_diagnostic(row_length, n_rows, model_levels) &
  , zh(row_length,rows)                                 &
  , zht(row_length,rows)

  REAL ::                        &
    bl_type_1(row_length,rows)   &! In Indicator set to 1.0 if stable b.l.
                                  !    diagnosed, 0.0 otherwise.
  , bl_type_2(row_length,rows)   &! In Indicator set to 1.0 if Sc over
                                  !    stable surface layer diagnosed,
                                  !    0.0 otherwise.
  , bl_type_3(row_length,rows)   &! In Indicator set to 1.0 if well
                                  !    mixed b.l. diagnosed, 0.0 otherwise.
  , bl_type_4(row_length,rows)   &! In Indicator set to 1.0 if decoupled
                                  !    Sc layer (not over cumulus) diagnosed,
                                  !    0.0 otherwise.
  , bl_type_5(row_length,rows)   &! In Indicator set to 1.0 if decoupled
                                  !    Sc layer over cumulus diagnosed,
                                  !    0.0 otherwise.
  , bl_type_6(row_length,rows)   &! In Indicator set to 1.0 if a cumulus
                                  !    capped b.l. diagnosed, 0.0 otherwise.
  , bl_type_7(row_length,rows)   &! In Indicator set to 1.0 if a
                                  !    shear-dominated b.l. diagnosed, 0.0
                                  !    otherwise.
  , bl_alltypes(row_length,rows)

  REAL ::                          &
    fqt(row_length,rows,bl_levels) &
  , ftl(row_length,rows,bl_levels) &
  , z0m(row_length,rows)           &! In Roughness length for mom
  , z0m_eff_gb(row_length,rows)    &! In Orographic roughness length
  , z0h_eff_gb(row_length,rows)    &! In Roughness length for heat

! variables needed for PC2 diagnostics
  , T(row_length, rows, model_levels)                           &
  , q(row_length, rows, wet_levels)                             &
  , cf(row_length, rows, wet_levels)                            &
  , cfl(row_length, rows, wet_levels)                           &
  , cff(row_length, rows, wet_levels)                           &

! _earliest arrays contain fields at start of imp_ctl
  , T_earliest(row_length, rows, model_levels)                  &
  , q_earliest(row_length, rows, wet_levels)                    &
  , qcl_earliest(row_length, rows, wet_levels)                  &
  , qcf_earliest(row_length, rows, wet_levels)                  &

! _earliest values contain the values of temperature, water contents
! and cloud fractions before the boundary layer call.
  , cf_earliest(row_length, rows, wet_levels)                   &
  , cfl_earliest(row_length, rows, wet_levels)                  &
  , cff_earliest(row_length, rows, wet_levels)                  &
  , p_theta_levels(row_length, rows, model_levels)              &

  , LWP(row_length, rows)              &! Liquid water path  (kg/m2)
  , IWP(row_length, rows)              &! Ice water path     (kg/m2)
  , sea_ice_htf(row_length,rows)       &
  , sice_mlt_htf(row_length,rows)      &
  , taux(row_length,rows,bl_levels)    &
  , tauy(row_length,rows,bl_levels)    &

  , area_cloud_fraction(row_length,rows,wet_levels) &

  , p_star(row_length,rows)            &
  , u10m(row_length,rows)              &
  , v10m(row_length,rows)              &
  , latent_heat(row_length,rows)       &
  , cca_2d(row_length,rows)            &
  , rho1(row_length,rows)              &! Air density at level 1/kg/m^3
  , qcl(row_length,rows,wet_levels)    &
  , qcf(row_length,rows,wet_levels)    &
  , q1p5m(row_length,rows)             &! InOut
  , t1p5m(row_length,rows)             &! InOut
  , surf_ht_flux_sice(row_length,rows) &
  , surf_ht_flux_gb(row_length,rows)   &

  , rhcpt(rhc_row_length,rhc_rows,wet_levels)       &
  , aerosol( 1-off_x:row_length+off_x               &
           , 1-off_y:rows+off_y                     &
           , model_levels)                          &
  , combined_cloud(row_length,rows,wet_levels)

! Logicals for SCM Diagnostics packages
  INTEGER ::               &
    nSCMDpkgs               ! No of SCM diagnostics packages

  LOGICAL ::               &
    L_SCMDiags(nSCMDpkgs)   ! Logicals for SCM diagnostics packages

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

! Declaration of new BL diagnostics.
  TYPE (Strnewbldiag) :: BL_diag

  ! Local variables...
  CHARACTER(LEN=*), PARAMETER :: RoutineName = 'Dgnstcs_Imp_Ctl2'

  INTEGER ::             &
    i                    &
  , j                    &
  , k                    &
  , error_code_ignored   &
  , kinvert               ! Vertical index for inverted arrays

  ! Work array
  REAL ::                                       &
    interp_data(row_length*rows*model_levels)   &
  , work_2d(row_length,rows)                    &
  , work_3d_w(row_length,rows,wet_levels)       &
  , work_3d(row_length,rows,model_levels)

  ! Global variables

  ! Local parameters and other physical constants
  REAL, PARAMETER :: LCRCP = LC/CP    ! Derived parameter.
                                      ! Lat ht of condensation/Cp.

  ! Variables used to calculate visibility stuff
  REAL ::                                          &
    Beta_LS_Rain(row_length,rows)                  &! Scattering in LS Rain.
  , Beta_LS_Snow(row_length,rows)                  &! Scattering in LS Snow.
  , Beta_C_Rain(row_length,rows)                   &! Scattering in Conv Rain
  , Beta_C_Snow(row_length,rows)                   &! Scattering in Conv Snow
  , Vis_Threshold(row_length,rows,1,n_vis_thresh)  &! FOG_FR works for n levels,
                                                    ! we want 1
  , PVis(row_length,rows,n_vis_thresh)

  REAL ::                                   &
    v_1p5m(row_length,rows)                 &! Visibility overall at 1.5m
  , v_no_precip(row_length,rows)            &! Visibility without precip.
  , v_ls_precip(row_length,rows)            &! Visibility in LS Precip.
  , v_c_precip(row_length,rows)             &! Visibility in Conv Precip.
  , v_probs(row_length,rows,1,n_vis_thresh)  ! Vis probs at first level
                                             ! (decimal fraction)

  ! Total cloud amounts.
  ! (X,X,1) - max. overlap. (X,X,2) - random overlap.
  REAL :: tot_cloud(row_length,rows,2)

  ! More output diagnostics
  REAL ::                             &
    layer_cloud1p5m(row_length,rows)  &! Layer cloud at 1.5m
                                       ! (decimal fraction)
  , qcl1p5m(row_length,rows)           ! Cloud liquid water content at
                                       ! 1.5m (kg/kg of air)

  ! Variables used in the calculation of the low_cloud, med_cloud
  ! and high_cloud diagnostics.
  INTEGER, PARAMETER :: num_cloud_types = 3

  REAL ::                           &
    h_split(num_cloud_types+1)      &
  , cloud_bound(num_cloud_types+1)

  INTEGER ::       &
    kk             &
  , level          &
  , low_bot_level  &
  , low_top_level  &
  , med_bot_level  &
  , med_top_level  &
  , high_bot_level &
  , high_top_level


  REAL ::                         &
    low_cloud(row_length,rows)    &
  , med_cloud(row_length,rows)    &
  , high_cloud(row_length,rows)   &
  , td1p5m(row_length,rows)       &! 1.5m dewpoint temperature
  , rh1p5m(row_length,rows)       &! 1.5m relative humidity
  , rhw1p5m(row_length,rows)      &! 1.5m relative humidity over liquid water
  , wspd10m(row_length,rows)      &! 10m wind speed
  , wdrn10m(row_length,rows)      &! 10m wind direction
  , bl_gust(row_length,rows)      &! Shear wind gust
  , alpha                         &! 'Angle' of wind
  , qst(1,1)                       ! A saturation mixing ratio

  ! Parameters
  LOGICAL, PARAMETER ::           &
    l_mixing_ratio = .FALSE.       ! Use mixing rations


!-----------------------------------------------------------------------
!     SCM Large Scale Cloud Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_lscld)) THEN

! To calculate low, medium and high cloud amounts, first need to set up
! variables low_bot_level, low_top_level, med_bot_level, etc. In global
! and limited-area models this is done during initialisation.

! Default settings for h_split as set in routine READLSTA
! (which is not called in the SCM)
!   low   : middle :  high   cloud model levels
! (1)->(2):(2)->(3):(3)->(4)

    h_split(1) =   111.0    ! ICAO 1000mb height (m)
    h_split(2) =  1949.0    ! ICAO  800mb height (m)
    h_split(3) =  5574.0    ! ICAO  500mb height (m)
    h_split(4) = 13608.0    ! ICAO  150mb height (m)

    ! Set up cloud type boundaries for low/medium/high cloud as in
    ! routine SETCONA (which is not called in the SCM). See this
    ! routine for more comments on what's going on here.
    DO kk=1, num_cloud_types+1

      level = 1
      DO WHILE ((eta_theta_levels(level)*z_top_of_model <= h_split(kk)) &
                 .AND. (level <= model_levels))
        level=level+1
      END DO

      IF (level >  model_levels) THEN
        PRINT*,'ERROR in ni_imp_ctl: level>model_levels',kk,level,model_levels
        level = model_levels
      END IF

      cloud_bound(kk) = level

    END DO

    low_bot_level  = cloud_bound(1)
    low_top_level  = cloud_bound(2) - 1
    med_bot_level  = cloud_bound(2)
    med_top_level  = cloud_bound(3) - 1
    high_bot_level = cloud_bound(3)
    high_top_level = cloud_bound(4) - 1

    IF (low_top_level >  cloud_levels) THEN
      PRINT*,'NI_Imp_Ctl ERROR: no of cloud levels less than '//  &
             'top of low', low_top_level, cloud_levels
    END IF

    IF (med_top_level >  cloud_levels) THEN
      PRINT*,'NI_Imp_Ctl ERROR: no of cloud levels less than '//  &
             'top of med', med_top_level, cloud_levels
    END IF

    IF (high_top_level >  cloud_levels) THEN
      PRINT*,'NI_Imp_Ctl ERROR: no of cloud levels less than '//  &
             'top of high', high_top_level, cloud_levels
    END IF

!---------------------------------------------------------------------
! In non-SCM runs the following stuff would be done in diagnostics_lscld

    ! Diagnostic 203: Low Cloud Amount
    ! Initialize cloud amount to lowest level value.
    DO j=1, rows
      DO i=1, row_length
        low_cloud(i,j)=area_cloud_fraction(i,j,LOW_BOT_LEVEL)
      END DO
    END DO

    ! Cloud amount is calculated under maximum overlap assumption.
    ! Limit cloud amount to maximum of 1.
    DO k=low_bot_level+1, low_top_level
      DO j=1, rows
        DO i=1, row_length
          low_cloud(i,j) = MAX(low_cloud(i,j),area_cloud_fraction(i,j,k))
          low_cloud(i,j) = MIN(low_cloud(i,j),1.0)
        END DO
      END DO
    END DO

! DEPENDS ON: scmoutput
    CALL scmoutput(low_cloud, 'lowcld'                                        &
      , 'Low cloud fraction', 'Fraction'                                      &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Diagnostic 204: Medium Cloud Amount
! Initialize cloud amount to lowest level value.
    DO j=1, rows
      DO i=1, row_length
        med_cloud(i,j)=area_cloud_fraction(i,j,med_bot_level)
      END DO
    END DO

! Cloud amount is calculated under maximum overlap assumption.
! Limit cloud amount to maximum of 1.
    DO k=med_bot_level+1, med_top_level
      DO j=1, rows
        DO i=1, row_length
          med_cloud(i,j) = MAX(med_cloud(i,j),area_cloud_fraction(i,j,k))
          med_cloud(i,j) = MIN(med_cloud(i,j),1.0)
        END DO
      END DO
    END DO

! DEPENDS ON: scmoutput
    CALL scmoutput(med_cloud, 'medcld'                                        &
      , 'Medium cloud fraction', 'Fraction'                                   &
      , t_avg, d_sl, default_streams, '', RoutineName)

    ! Diagnostic 205: High Cloud Amount
    ! Initialize cloud amount to lowest level value.
    DO j=1, rows
      DO i=1, row_length
        high_cloud(i,j)=area_cloud_fraction(i,j,high_bot_level)
      END DO
    END DO

    ! Cloud amount is calculated under maximum overlap assumption.
    ! Limit cloud amount to maximum of 1.
    DO k=high_bot_level+1, high_top_level
      DO j=1, rows
        DO i=1, row_length
          high_cloud(i,j) = MAX(high_cloud(i,j),area_cloud_fraction(i,j,k))
          high_cloud(i,j) = MIN(high_cloud(i,j),1.0)
        END DO
      END DO
    END DO

! DEPENDS ON: scmoutput
    CALL scmoutput(high_cloud, 'highcld'                                      &
      , 'High cloud fraction', 'Fraction'                                     &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Total Cloud Amount RANDOM Overlap
! DEPENDS ON: r2_calc_total_cloud_cover
    CALL r2_calc_total_cloud_cover                                            &
      ( row_length*rows, wet_levels, nclds, IP_CLOUD_MIX_RANDOM               &
      , combined_cloud(1,1,1), tot_cloud(1,1,2), row_length*rows, wet_levels )

! Total Cloud Amount MAX/RANDOM Overlap
! DEPENDS ON: r2_calc_total_cloud_cover
    CALL r2_calc_total_cloud_cover                                            &
      ( row_length*rows, wet_levels, nclds, IP_CLOUD_MIX_MAX                  &
      , combined_cloud(1,1,1), tot_cloud(1,1,1), row_length*rows              &
      , wet_levels )

! DEPENDS ON: scmoutput
    CALL scmoutput(tot_cloud(1,1,2), 'tcarndm'                                &
      , 'Total cloud amount (random overlap)', 'Fraction'                     &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(tot_cloud(1,1,1), 'tcamxrn'                                &
      , 'Total cloud amount (max. random)', 'Fraction'                        &
      , t_avg, d_sl, default_streams, '', RoutineName)

    ! Note that stash equivalents 9,18[123] are identical to 3,18[123]
    ! and so done below stash 9,226
    DO k=1, wet_levels
      DO j=1, rows
        DO i=1, row_length
          IF (area_cloud_fraction(i, j, k) <= 0.0)  THEN
            work_3d_w(i, j, k) = 0.0
          ELSE
            work_3d_w(i, j, k) = 1.0
          END IF
        END DO ! i
      END DO ! j
    END DO ! k

!
! DEPENDS ON: scmoutput
    CALL scmoutput( work_3d_w, 'lyrcldfreq'                                   &
      , 'Layer cloud frequency indicator', 'Indicator'                        &
      , t_avg, d_wet, default_streams, '', RoutineName)

! Stash 9,228
! DEPENDS ON: scmoutput
    CALL scmoutput(rhcpt, 'rhcpt'                                             &
      , 'Critical relative humidity', '%'                                     &
      , t_avg, d_wet, default_streams, '', RoutineName)

! Stash 9,229
    DO k=1, wet_levels
! DEPENDS ON: qsat
      CALL qsat (work_2d, T(1,1,K), p_theta_levels(1,1,k), row_length*rows)

      DO j=1, rows
        DO i=1, row_length

          ! relative humidity in per cent
          work_3d_w(i, j, k) = q(i, j, k) / work_2d(i, j) *100.0

          ! Supersaturation (>100%) can occur with mixed phase
          ! scheme but negative humidity is removed from the
          ! diagnostic:
          IF ( work_3d_w(i, j, k) < 0.0) THEN
            work_3d_w(i, j, k) = 0.0
          END IF

        END DO ! i
      END DO ! j
    END DO ! k

!
! DEPENDS ON: scmoutput
    CALL scmoutput(work_3d_w, 'rhaftcld'                                      &
      , 'Relative humidity after main cloud', '%'                             &
      , t_avg, d_wet, default_streams, '', RoutineName)

! Stash 9,231
    DO k=1, wet_levels

! NB: Convention in Sect 70 (Radiation) is to invert levels,
!     1 at top. Combined_cloud is calculated this way but
!     re-inverted for STASH.

      kinvert = wet_levels+1-k

      DO j=1, rows
        DO i=1, row_length
          work_3d_w(i, j, k) = combined_cloud(i,j,kinvert)
        END DO ! i
      END DO ! j
    END DO ! k

! DEPENDS ON: scmoutput
    CALL scmoutput(work_3d_w, 'combca_ls'                                     &
      , 'Combined cloud amount in each layer', 'Fraction'                     &
      , t_avg, d_wet, default_streams, '', RoutineName)


! End of stuff that would be done in diagnostics_lscld               !
!---------------------------------------------------------------------


  END IF ! L_SCMDiags(SCMDiag_lscld)


!-----------------------------------------------------------------------
!     SCM Surface Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_surf)) THEN

! Stash 3,255
! DEPENDS ON: scmoutput
    CALL SCMoutput(q1p5m, 'qt1p5m'                                            &
      , '1.5m total water kg water/kg air', 'kgH20/kgAIR'                     &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 3,254
! DEPENDS ON: scmoutput
    CALL SCMoutput(t1p5m, 'tl1p5m'                                            &
      , '1.5m liquid temperature', 'K'                                        &
      , t_avg, d_sl, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_surf)

!-----------------------------------------------------------------------
!     SCM Boundary Layer OR Large Scale Cloud Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_bl) .OR. L_SCMDiags(SCMDiag_lscld)) THEN

! DEPENDS ON: scmoutput
    CALL SCMoutput(lwp, 'LWP'                                                 &
      , 'Liquid water path', 'kg/m2'                                          &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL SCMoutput(iwp, 'IWP'                                                 &
      , 'Ice water path', 'kg/m2'                                             &
      , t_avg, d_sl, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_bl) .OR. L_SCMDiags(SCMDiag_lscld)


!---------------------------------------------------------------------
! In non-SCM runs the following stuff would be done in diagnostics_bl

  ! This call to ls_cld converts 1.5m TL and QT to T and q,
  ! assuming p_star and rhcrit at model level 1 approximates to
  ! 1.5m. This should be done even if ssfm=.false. It also
  ! calculates layer_cloud1p5m and qcl1p5m.
! DEPENDS ON: ls_cld
  CALL ls_cld                                                                 &
    ( p_star, rhcpt, 1, bl_levels, rhc_row_length, rhc_rows                   &
    , ntml, cumulus                                                           &
    , l_mixing_ratio, t1p5m, layer_cloud1p5m, q1p5m, qcf, qcl1p5m             &
    , interp_data(2*row_length*rows+1), interp_data(3*row_length*rows+1)      &
    , error_code_ignored)

!-----------------------------------------------------------------------
!     SCM Large Scale Cloud OR Surface Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_lscld) .OR. L_SCMDiags(SCMDiag_surf)) THEN

! DEPENDS ON: scmoutput
    CALL scmoutput(layer_cloud1p5m, 'lca1p5m'                                 &
      , '1.5m layer cloud amount', 'Fraction'                                 &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(qcl1p5m, 'qcl1p5m'                                         &
       , '1.5m cloud water', 'kg/kg'                                          &
       , t_avg, d_sl, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_lscld) .OR. L_SCMDiags(SCMDiag_surf)


!-----------------------------------------------------------------------
!     SCM Surface Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_surf)) THEN

    ! 1.5m Fog Fraction and Mist Fraction
    DO i=1, row_length
      DO j=1, rows
        DO k=1, n_vis_thresh
          vis_threshold(i,j,1,k) = vis_thresh(k)
        END DO
      END DO
    END DO

! DEPENDS ON: fog_fr
    CALL fog_fr                                                               &
      ( p_star, rhcrit, 1, row_length*rows, t1p5m, aerosol, L_murk, q1p5m     &
      , qcl, qcf, vis_threshold, PVis, n_vis_thresh )

     ! Calculate scattering coefficients due to precipitation
    CALL beta_precip                                                          &
      ( ls_rain, ls_snow, conv_rain,conv_snow,qcf(1,1,1), rho1, T1p5m, p_star &
      , plsp, cca_2d, .FALSE., .TRUE., row_length*rows, row_length*rows, 1    &
      , beta_ls_rain, beta_ls_snow, beta_c_rain, beta_c_snow                  &
      , error_code_ignored )

     ! Calculate screen level probability of visibility
     ! less than thresholds
! DEPENDS ON: calc_vis_prob
    CALL calc_vis_prob                                                        &
      ( p_star, rhcrit, 1, row_length*rows, row_length*rows, t1p5m, aerosol   &
      , l_murk, q1p5m, qcl, qcf, vis_thresh, n_vis_thresh, plsp, cca_2d       &
      , .FALSE., beta_ls_rain, beta_ls_snow, beta_c_rain, beta_c_snow         &
      , v_probs, error_code_ignored )

     ! Visibility at 1.5 m including precipitation
! DEPENDS ON: visbty
    CALL visbty                                                               &
      ( p_star, T1p5m, q1p5m, Qcl, Qcf, Aerosol, calc_prob_of_vis, RHcrit     &
      , L_murk, row_length*rows, v_no_precip )

! DEPENDS ON: vis_precip
    CALL vis_precip                                                           &
      ( v_no_precip, plsp, cca_2d, .FALSE., Beta_ls_Rain, Beta_ls_Snow        &
      , Beta_C_Rain, Beta_C_Snow, row_length*rows, row_length*rows, 1         &
      , v_1p5m, v_ls_precip, v_C_Precip, error_code_ignored)

! DEPENDS ON: scmoutput
    CALL scmoutput(v_probs(1,1,1,1), 'pfog1p5m'                               &
      , 'Probability of fog at 1.5m', '-'                                     &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(v_probs(1,1,1,2), 'pmist1p5m'                              &
      , 'Probability of mist at 1.5m', '-'                                    &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(v_probs, 'pvisthresh'                                      &
      , 'Probability of vis<=threshold', '-'                                  &
      , t_avg, d_vis, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(v_no_precip, 'visnop1p5m'                                  &
      , '1.5m visibility outside precip', 'm'                                 &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(v_ls_precip, 'vislsp1p5m'                                  &
      , '1.5m visibility in LS precip', 'm'                                   &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(v_c_precip, 'viscp1pm5'                                    &
      , '1.5m visibility in conv precip', 'm'                                 &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(v_1p5m, 'vis1p5m'                                          &
      , '1.5m visibility', 'm'                                                &
      , t_avg, d_sl, default_streams, '', RoutineName)

    ! Calculate relative humidities at screen level
    DO j=1, rows
      DO i=1, row_length

! DEPENDS ON: qsat
        CALL qsat (qst(1,1), t1p5m(i,j), p_star(i,j), 1)
        rh1p5m(i,j) = q1p5m(i,j)/qst(1,1)*100.0
! DEPENDS ON: qsat_wat
        CALL qsat_wat (qst(1,1), t1p5m(i,j), p_star(i,j), 1)
        rhw1p5m(i,j) = q1p5m(i,j)/qst(1,1)*100.0

      END DO
    END DO

! DEPENDS ON: scmoutput
   CALL scmoutput(rh1p5m, 'rh1p5m'                                            &
     , 'Relative humidity at 1.5m', '%'                                       &
     , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
   CALL scmoutput(rhw1p5m, 'rhw1p5m'                                          &
     , 'Relative humidity wrt H2O at 1.5m', 'W/m2'                            &
     , t_avg, d_sl, default_streams, '', RoutineName)

! Calculate dewpoint

! DEPENDS ON: dewpnt
   CALL dewpnt (q1p5m, p_star, t1p5m, row_length*rows, td1p5m)

! DEPENDS ON: scmoutput
   CALL scmoutput(td1p5m, 'td1p5m'                                            &
     , '1.5m dewpoint temperature', 'W/m2'                                    &
     , t_avg, d_sl, default_streams, '', RoutineName)

     ! Calculate the magnitude and direction (by meteorological
     ! convention) of the wind at 10m from the 10m u and v
     ! components
     DO j=1, rows
        DO i=1, row_length
           wspd10m(i,j) = SQRT(u10m(i,j)**2+v10m(i,j)**2)
           IF (wspd10m(i,j) == 0.0) THEN
              wdrn10m(i,j) = 0.0
           ELSE
              alpha = ASIN(v10m(i,j)/wspd10m(i,j))*                           &
                   recip_pi_over_180
              IF (u10m(i,j) >  0.0) THEN
                 wdrn10m(i,j) = 270 - alpha
              ELSE
                 wdrn10m(i,j) = 90 + alpha
              END IF
           END IF
        END DO
     END DO

! DEPENDS ON: scmoutput
     CALL scmoutput(wspd10m, 'wspd10m'                                        &
       , '10m wind speed', 'm/s'                                              &
       , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
     CALL scmoutput(wdrn10m, 'wdrn10m'                                        &
       ,'10m wind direction', 'degs'                                          &
       , t_avg, d_sl, default_streams, '', RoutineName)

     ! Boundary layer gust diagnostic
     DO j=1, rows
        DO i=1, row_length
           bl_gust(i,j) = wspd10m(i,j)+2.0*2.5                                &
                *SQRT(SQRT(taux(i,j,1)**2+tauy(i,j,1)**2))
        END DO
     END DO

! DEPENDS ON: scmoutput
     CALL scmoutput(bl_gust, 'gust10m'                                        &
       , '10m Gust', 'm/s'                                                    &
       , t_avg, d_sl, default_streams, '', RoutineName)


  END IF ! L_SCMDiags(SCMDiag_surf)


!-----------------------------------------------------------------------
!     SCM Sea Points Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_sea)) THEN

! DEPENDS ON: scmoutput
    CALL scmoutput(sice_mlt_htf, 'sice_mlt_htf'                               &
      , 'Heat flux due to melting sea ice', 'W/m2'                            &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(sea_ice_htf, 'sea_ice_htf'                                 &
      , 'Heat flux through sea ice', 'W/m2'                                   &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(surf_ht_flux_sice, 'surf_ht_flux_si'                       &
      , 'Net downward heat flux into sea frctn', 'W/m2'                       &
      , t_avg, d_sl, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_sea)

!-----------------------------------------------------------------------
!     SCM Surface Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_surf)) THEN

! DEPENDS ON: scmoutput
    CALL scmoutput(latent_heat, 'lat_ht'                                      &
      , 'Surface latent heat flux', 'W/m2'                                    &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(surf_ht_flux_gb, 'surf_ht_flux'                            &
      , 'Net downward heat flux at surface', 'W/m2'                           &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(z0m, 'z0m'                                                 &
      , 'Roughness length for momentum', 'm'                                  &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(z0h_eff_gb, 'z0h'                                          &
      , 'Effective roughness length for heat', 'm'                            &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(z0m_eff_gb, 'z0m_eff'                                      &
      , 'Effective roughness length for momentum', 'm'                        &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(q1p5m, 'q1p5m'                                             &
      , '1.5m specific humidity', 'kg/kg'                                     &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(t1p5m, 't1p5m'                                             &
      , '1.5m temperature', 'K'                                               &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(t1p5m, 't1p5m_max'                                         &
      , 'Max 1.5m temperature', 'K'                                           &
      , t_max, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(t1p5m, 't1p5m_min'                                         &
      , 'Min 1.5m temperature', 'K'                                           &
      , t_min, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(u10m, 'u10m'                                               &
      , 'Zonal 10m wind', 'm/s'                                               &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(v10m, 'v10m'                                               &
      , 'Meridional 10m wind', 'm/s'                                          &
      , t_avg, d_sl, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_surf)
!
!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_bl)) THEN

! Stash 3, 025
! DEPENDS ON: scmoutput
    CALL scmoutput(zh, 'zh'                                                   &
      , 'Boundary layer depth after B.layer', 'm'                             &
      , t_avg, d_sl, default_streams, '', RoutineName)

!       Since model level fields are on rho_levels and the surface
!       is a theta_level, these are output as separate diagnostics

! Stash 3, 217
    DO j=1, rows
      DO i=1, row_length
        work_2d(i, j) = ftl(i, j, 1)
      END DO ! i
    END DO ! j

! DEPENDS ON: scmoutput
    CALL scmoutput(work_2d, 'ftl_surf'                                        &
      , 'Surface sensible heat flux from B.layer', 'W/m2'                     &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 3, 223
    DO j=1, rows
      DO i=1, row_length
        work_2d(i,j) = fqt(i,j,1)
      END DO ! i
    END DO ! j

! DEPENDS ON: scmoutput
    CALL scmoutput(work_2d, 'fqt_surf'                                        &
      , 'Surface sensible moisture flux from B.layer', 'kg/m2/s'              &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 3,222
! DEPENDS ON: scmoutput
    CALL scmoutput(fqt, 'fqt_bl'                                              &
      , 'Sensible moisture flux', 'kg/m2/s'                                   &
      , t_avg, d_bl, default_streams, '', RoutineName)

! Stash 3,473
! DEPENDS ON: scmoutput
    CALL scmoutput(BL_diag%tke, 'TKE'                                         &
      , 'BL diagnostic of Turbulent Kinetic Energy', 'm2/s'                   &
      , t_avg, d_bl, default_streams, '', RoutineName)

  IF (Fric_heating /= 0) THEN
! Stash 3,188
! DEPENDS ON: scmoutput
    CALL scmoutput(BL_diag%dTfric, 'dt_fric'                                  &
      , 'T increment from turbulence dissipation', 'K'                        &
      , t_avg, d_bl, default_streams, '', RoutineName)
  END IF

! Stash 3,304
! DEPENDS ON: scmoutput
    CALL scmoutput(zht, 'zht'                                                 &
      , 'Turbulent mixing height after B.layer', 'm'                          &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 3,305
! DEPENDS ON: scmoutput
    CALL scmoutput(bl_type_1, 'bl_type_1'                                     &
      , 'Boundary layer type: stable', 'Indicator'                            &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 3,306
! DEPENDS ON: scmoutput
    CALL scmoutput(bl_type_2, 'bl_type_2'                                     &
      , 'Boundary layer type: Sc over stable', 'Indicator'                    &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 3,307
! DEPENDS ON: scmoutput
    CALL scmoutput(bl_type_3, 'bl_type_3'                                     &
      , 'Boundary layer type: well mixed', 'Indicator'                        &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 3,308
! DEPENDS ON: scmoutput
    CALL scmoutput(bl_type_4, 'bl_type_4'                                     &
      , 'Boundary layer type: decoup Sc not over Cu', 'Indicator'             &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 3,309
! DEPENDS ON: scmoutput
    CALL scmoutput(bl_type_5, 'bl_type_5'                                     &
      , 'Boundary layer type: decoup Sc over Cu', 'Indicator'                 &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 3,310
! DEPENDS ON: scmoutput
    CALL scmoutput(bl_type_6, 'bl_type_6'                                     &
      , 'Boundary layer type: cumulus capped', 'Indicator'                    &
      , t_avg, d_sl, default_streams, '', RoutineName)

! Stash 3,340
! DEPENDS ON: scmoutput
    CALL scmoutput(bl_type_7, 'bl_type_7'                                     &
      , 'Boundary layer type: shear driven', 'Indicator'                      &
      , t_avg, d_sl, default_streams, '', RoutineName)

! DEPENDS ON: scmoutput
    CALL scmoutput(bl_alltypes, 'bl_alltypes'                                 &
      , 'Boundary layer types', ''                                            &
      , t_inst, d_sl, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_bl)
!
!-----------------------------------------------------------------------
!     SCM Boundary Layer OR Increments Diagnostics Package
!-----------------------------------------------------------------------
  IF (L_SCMDiags(SCMDiag_bl) .OR. L_SCMDiags(SCMDiag_incs)) THEN

! Stash 3,181 and 9,181
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,k) = T(i,j,k) - T_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

! DEPENDS ON: scmoutput
    CALL scmoutput(work_3d, 'dt_pc2blls'                                      &
      , 'Temperature inc PC2+bdy layer+ls cld', 'K'                           &
      ,  t_avg, d_all, default_streams, '', RoutineName)

! Stash 3,182 and 9,182
    DO k=1, wet_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d_w(i,j,k) = q(i,j,k) - q_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

! DEPENDS ON: scmoutput
    CALL scmoutput(work_3d_w, 'dq_pc2blls'                                    &
      , 'Specific humidity inc PC2+bdy layer+ls cld', 'kg/kg'                 &
      , t_avg, d_wet, default_streams, '', RoutineName)

! Stash 3,183 and 9,183
    DO k=1, wet_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d_w(i,j,k) = qcl(i,j,k) - qcl_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

! DEPENDS ON: scmoutput
    CALL scmoutput(work_3d_w, 'dqcl_pc2blls'                                  &
      , 'QCL increment PC2+bdy layer+ls cld', 'kg/kg'                         &
      , t_avg, d_wet, default_streams,  '',  RoutineName)

! Stash 3,184
    DO k=1, wet_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d_w(i,j,k) = qcf(i,j,k) - qcf_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

! DEPENDS ON: scmoutput
    CALL scmoutput(work_3d_w, 'dqcf_pc2blls'                                  &
      , 'QCF increment PC2+bdy layer+ls cld', 'kg/kg'                         &
      , t_avg, d_wet, default_streams, '', RoutineName)

! Stash 3,185
! DEPENDS ON: scmoutput
    CALL scmoutput(u_incr_diagnostic, 'du_bl'                                 &
      , 'U wind increment bdy layer', 'm/s'                                   &
      , t_avg, d_all, default_streams, '', RoutineName)

! Stash 3,186
! DEPENDS ON: scmoutput
    CALL scmoutput(v_incr_diagnostic, 'dv_bl'                                 &
      , 'V wind increment bdy layer', 'm/s'                                   &
      , t_avg, d_all, default_streams, '', RoutineName)

! Stash 3,192
    DO k=1, wet_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d_w(i,j,k) = cf(i,j,k) - cf_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

! DEPENDS ON: scmoutput
    CALL scmoutput(work_3d_w, 'dbcf_bl'                                       &
      , 'Bulk cloud fraction increment bdy layer', 'Fraction'                 &
      ,  t_avg, d_wet, default_streams, '', RoutineName)

! Stash 3,193
    DO k=1, wet_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d_w(i,j,k) = cfl(i,j,k) - cfl_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

! DEPENDS ON: scmoutput
    CALL scmoutput(work_3d_w, 'dcfl_bl'                                       &
      , 'Liquid cloud fraction increment bdy layer', 'Fraction'               &
      , t_avg, d_wet, default_streams, '', RoutineName)

! Stash 3,194
    DO k=1, wet_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d_w(i,j,k) = cff(i,j,k) - cff_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

! DEPENDS ON: scmoutput
    CALL scmoutput(work_3d_w, 'dcff_bl'                                       &
      , 'Frozen cloud fraction increment bdy layer', 'Fraction'               &
      , t_avg, d_wet, default_streams, '', RoutineName)

! Stash 3,189
    DO k=1, wet_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d(i,j,k) = T(i,j,k) -   T_earliest(i,j,k)                     &
               - LCRCP * ( qcl(i,j,k) - qcl_earliest(i,j,k) )
        END DO ! i
      END DO ! j
    END DO ! k

    IF (model_levels > wet_levels)  THEN
      DO k=(wet_levels + 1), model_levels
        DO j=1, rows
          DO i=1, row_length
            work_3d(i,j,k) = T(i,j,k) - T_earliest(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
    END IF ! model_levels > wet_levels

! DEPENDS ON: scmoutput
    CALL scmoutput(work_3d, 'dlwt_bl'                                         &
      , 'Liquid water temp increment bdy layer', 'K'                          &
      , t_avg, d_all, default_streams, '', RoutineName)

! Stash 3,190
    DO k=1, wet_levels
      DO j=1, rows
        DO i=1, row_length
          work_3d_w(i,j,k) =   q(i,j,k) -   q_earliest(i,j,k)                 &
                           + qcl(i,j,k) - qcl_earliest(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

! DEPENDS ON: scmoutput
    CALL scmoutput(work_3d_w, 'tqinc_bl'                                      &
      , 'Total (liquid) water increment bdy layer', 'kg/kg'                   &
      , t_avg, d_all, default_streams, '', RoutineName)

  END IF ! L_SCMDiags(SCMDiag_bl) .OR. L_SCMDiags(SCMDiag_incs)


! End of stuff that would be done in diagnostics_bl
!---------------------------------------------------------------------

  RETURN

END SUBROUTINE dgnstcs_imp_ctl2

