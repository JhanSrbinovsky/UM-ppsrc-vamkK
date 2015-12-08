! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Advection through falling of ice
! and rain
! Subroutine Interface:
MODULE lsp_fall_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_fall(                                                    &
  points, timestep,                                                     &
                                          ! Number of points and tstep
  qcf_cry_0, qcf_agg_0, frac_agg,                                       &
                                          ! Ice contents
  qrain_0,qgraup_0, t,                                                  &
                                          ! Rain and graupel contents
  snow_agg,snow_cry,rainrate,grauprate,                                 &
                                             ! Sedimentation into layer
  snowt_agg,snowt_cry,rainratet,graupratet,                             &
                                                 ! Sedim. out of layer
  vf_agg, vf_cry, vf_rain, vf_graup,                                    &
                                          ! Fall speeds of hydrometeors
  area_clear,area_ice,cfice,cficei,                                     &
                                          ! Cloud fraction information
  frac_ice_above,                                                       &
                                            ! at start of microphy. ts
  cf, cfl, cff,                                                         &
                                          ! Current cloud fractions for
                                          ! updating
  rho, rhor, tcgi, tcgci,                                               &
                                          ! Parametrization information
  corr, dhi, dhir, rainfrac,                                            &
                                          ! Parametrization information
  d_qcf_cry_dt,d_qcf_agg_dt,                                            &
                                          ! Mass transfer diagnostics
  d_qrain_dt,d_qgraup_dt,                                               &
                                          ! Mass transfer diagnostics
  iterations, one_over_tsi,                                             &
                                          ! 1/(timestep*iterations)
  cftransfer,cfftransfer,                                               &
                                          ! Cloud transfer diagnostics
  uk, vk, ukp1, vkp1,                                                   &
                                          ! Winds for calc of wind-shear
  r_theta_levels_c, fv_cos_theta_latitude_c                             &
                                          ! Grid info
  )

  ! Microphysics and Clouds modules
  USE mphys_ice_mod,        ONLY: m0
  USE mphys_constants_mod,  ONLY: cx, constp
  USE mphys_inputs_mod,     ONLY: l_rainfall_as, l_psd, l_warm_new,     &
                                  l_mcr_qgraup, l_mcr_qcf2,             &
                                  l_mcr_qrain 
  USE cloud_inputs_mod,     ONLY: pc2_falliceshear_method,              &
                                  cff_spread_rate
  USE mphys_bypass_mod,     ONLY: l_crystals, mp_dell, mp_delp
  USE pc2_constants_mod,    ONLY: wind_shear_factor,                    &
                                  original_but_wrong, ignore_shear,     &
                                  real_shear

  ! Dr Hook Modules
  USE yomhook,              ONLY: lhook, dr_hook
  USE parkind1,             ONLY: jprb, jpim

  ! Large scale precip modules
  USE lsp_moments_mod,      ONLY: lsp_moments
  USE lsp_sedim_eulexp_mod, ONLY: lsp_sedim_eulexp

  IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of ice particle and
!   raindrop fall

! Method:
!   Calculate particle fall speeds and solve the advection equation
!   for mixing ratios following Rotstayns method.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! Both small and large ice, and raindrops, will fall, hence require
! advection downwards. We include both code for the two-ice prognostics
! and the single ice prognostic that is diagnostically split. Although
! the advection methods are very similar, they are different enough
! (in their calculation of fall-speed of the ice from above) to
! use different branches of code.

! Subroutine Arguments

  INTEGER, INTENT(IN) ::                                                &
    points,                                                             &
                          ! Number of points to calculate
    iterations
                          ! Number of microphysics iterations
                          ! (for rate diagnostic)

  REAL, INTENT(IN) ::                                                   &
    timestep,                                                           &
                          ! Timestep / s
    frac_agg(points),                                                   &
                          ! Fraction of ice mass that is aggregates
    snow_cry(points),                                                   &
                          ! Crystal flux into layer / kg m-2 s-1
    snow_agg(points),                                                   &
                          ! Aggregate flux into layer / kg m-2 s-1
    rainrate(points),                                                   &
                          ! Rain flux into layer / kg m-2 s-1
    grauprate(points),                                                  &
                          ! Graupel flux into layer / kg m-2 s-1
    area_clear(points),                                                 &
                          ! Fraction of gridbox with clear sky
    area_ice(points),                                                   &
                          ! Frac of gridbox with ice but not liquid
    cfice(points),                                                      &
                          ! Fraction of gridbox with ice cloud
    cficei(points),                                                     &
                          ! 1/Fraction of gridbox with ice cloud
    rainfrac(points),                                                   &
                          ! Rain fraction in grid-box
    frac_ice_above(points),                                             &
                               ! Ice cloud fraction in layer above
    rho(points),                                                        &
                          ! Air density / kg m-3
    rhor(points),                                                       &
                          ! 1 / Air density / m3 kg-1
    dhi(points),                                                        &
                          ! Timestep / thickness of model layer / s m-1
    dhir(points),                                                       &
                          ! 1/dhi / m s-1
    tcgi(points),                                                       &
                          ! 1/tcg (no units)
    tcgci(points),                                                      &
                          ! 1/tcgc (no units)
    corr(points),                                                       &
                          ! Air density fall speed correction (no units)
    uk(points),                                                         &
                          ! U wind on level k
    vk(points),                                                         &
                          ! V wind on level k
    ukp1(points),                                                       &
                          ! U wind on level k+1
    vkp1(points),                                                       &
                          ! V wind on level k+1
    r_theta_levels_c(points),                                           &
                          ! Distance from centre of Earth and ... 
    fv_cos_theta_latitude_c(points),                                    &
                          ! ... grid info for working out gridbox size.
    one_over_tsi          
                          ! 1/timestep(*iterations)

  REAL, INTENT(INOUT) ::                                                &
    qcf_cry_0(points),                                                  &
                          ! Ice crystal mixing ratio / kg kg-1
    qcf_agg_0(points),                                                  &
                          ! Ice aggregate mixing ratio / kg kg-1
    qrain_0(points),                                                    &
                          ! Rain mixing ratio / kg kg-1
    qgraup_0(points),                                                   &
                          ! Graupel mixing ratio / kg kg-1
    t(points),                                                          &
                          ! Temperature / K
    snowt_cry(points),                                                  &
                          ! Snowfall rate out of this layer / kg m-2 s-1
    snowt_agg(points),                                                  &
                            ! for crystals and aggregates
    rainratet(points),                                                  &
                          ! Rain rate out of this layer / kg m-2 s-1
    graupratet(points),                                                 &
                          ! Graupel rate out of this layer / kg m-2 s-1
    vf_cry(points),                                                     &
                          ! On input: Fall speed of hydrometeors
    vf_agg(points),                                                     &
                                    ! entering the current layer / m s-1
    vf_rain(points),                                                    &
                          ! On output: Fall speed of hydrometeors
    vf_graup(points),                                                   &
                                    ! leaving the current layer / m s-1
    cf(points),                                                         &
                          ! Current cloud fraction
    cfl(points),                                                        &
                          ! Current liquid cloud fraction
    cff(points)       ! Current ice cloud fraction

  REAL, INTENT(INOUT) ::                                                &
    d_qcf_cry_dt(points),                                               &
                             ! Rate of change of crystal, aggregate,
    d_qcf_agg_dt(points),                                               &
                               ! rain and graupel mixing ratios due
    d_qrain_dt(points),                                                 &
                               ! to sedimentation / kg kg-1 s-1
    d_qgraup_dt(points),                                                &
    cftransfer(points),                                                 &
                           ! Cloud fraction increment this tstep
    cfftransfer(points)! Ice cloud fraction inc this tstep

! Local Variables

  INTEGER ::                                                            &
    i                 ! Loop counter for points

  REAL ::                                                               &
    qcf_cry(points),                                                    &
                          ! Working copy of qcf_cry_0 / kg kg-1
    qcf_agg(points),                                                    &
                          ! Working copy of qcf_agg_0 / kg kg-1
    qrain(points),                                                      &
                          ! Working copy of q_rain / kg kg-1
    qgraup(points),                                                     &
                          ! Working copy of q_graup / kg kg-1
    fqi_agg(points),                                                    &
                          ! Fall speed of aggregates out of the
                          ! current layer / m s-1
    fqi_cry(points),                                                    &
                          ! Fall speed of crystals out of the
                          ! current layer / m s-1
    fqirqi2_agg(points),                                                &
                            ! Fraction of aggregate and crystal mass
    fqirqi2_cry(points),                                                &
                              ! that remains in the current layer
                              ! after sedimentation
    fqirqi_agg,                                                         &
                          ! Flux of aggregates out of layer / m s-1
    fqirqi_cry,                                                         &
                          ! Flux of crystals out of layer / m s-1
    fqi_rain(points),                                                   &
                          ! Bulk fall speed of rain and graupel
    fqi_graup(points),                                                  &
                            ! out of the current layer / m s-1
    mixratio_fromabove(points),                                         &
                                   ! Estimate of the mixing ratio
                          ! of ice in the layer above / kg kg-1
    temp3(points),                                                      &
                          ! Fraction of layer ice falls through
    overhang,                                                           &
                          ! Ice cloud fraction that overhangs
                          ! the current layer from above
    deltacf(points),                                                    &
                          ! Change in cloud fraction across timestep
    deltacff(points),                                                   &
                          ! Change in ice cloud fraction across tstep
    m_bi_di(points),                                                    &
                          ! bi+di'th moment of the ice particle dist.
                          ! Change in ice cloud fraction across tstep
    lamr(points),                                                       &
                          ! Lambda for rain.
    lamrh1(points),                                                     &
                          ! Lambda for rain + h1r
    lamrh2(points),                                                     &
                          ! Lambda for rain + h2r
    save_vf_agg(points),                                                &
    save_vf_cry(points),                                                &

    ! For calculating lateral displacement of falling ice due to shear.
    !------------------------------------------------------------------ 
    mwfv, dudz, dvdz, shear, horiz_scale, lateral_disp, cff_perimeter

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('LSP_FALL',zhook_in,zhook_handle)

      !-----------------------------------------------
      ! Calculate moment of size distribution if appropriate
      !-----------------------------------------------

  IF (l_psd) THEN
        ! Calculate the bi+di (cx(82)) moment of the
        ! ice particle size distribution
    CALL lsp_moments(points,rho,t,qcf_agg_0,cficei,cx(82),m_bi_di)
  END IF

  DO i = 1, points
        !-----------------------------------------------
        ! Copy input water contents (allows sequential updates)
        !-----------------------------------------------
    qcf_cry(i) = qcf_cry_0(i)
    qcf_agg(i) = qcf_agg_0(i)
    qrain(i)   = qrain_0(i)
    qgraup(i)  = qgraup_0(i)

    ! Save the fall speeds coming in from layer above before 
    ! they get overwritten with fall speed passed to layer below.
    ! Will be used in updating falling ice cloud fraction
    save_vf_agg(i)=vf_agg(i)
    save_vf_cry(i)=vf_cry(i)

    IF (qcf_agg(i) >  m0) THEN
          !-----------------------------------------------
          ! Estimate fall speed out of this layer
          !-----------------------------------------------
      IF (l_psd) THEN
            ! Use the generic PSD
            ! constp(82) = ci*ai
        fqi_agg(i) = constp(82) * corr(i) * m_bi_di(i)                  &
                   / (rho(i) * qcf_agg(i) * cficei(i))
      ELSE
        fqi_agg(i) = constp(24) * corr(i) *                             &
        (rho(i)*qcf_agg(i)*constp(25)*tcgi(i)*cficei(i))**cx(23)
      END IF  ! l_psd

    ELSE
          !-----------------------------------------------
          ! Ice content is small so it is numerically best to
          ! assume there is no ice here, so the fall speed is zero.
          !-----------------------------------------------
      fqi_agg(i)=0.0
    END IF  ! qcf_agg gt 0

    IF (qcf_cry(i) >  m0 .AND. l_crystals) THEN
         !-----------------------------------------------
         ! Estimate fall speed out of this layer
         !-----------------------------------------------
      fqi_cry(i) = constp(4) * corr(i)*                                 &
      (rho(i)*qcf_cry(i)*constp(5)*tcgci(i)*cficei(i))**cx(3)
    ELSE
          !-----------------------------------------------
          ! Ice content is small so it is numerically best to
          ! assume there is no ice here, so the fall speed is zero.
          !-----------------------------------------------
      fqi_cry(i)=0.0
    END IF  ! qcf_cry gt 0

  END DO  ! Points

  IF (.NOT. l_mcr_qcf2) THEN

    DO i = 1, points
          !-----------------------------------------------
          ! Calculate fall speed from above as an average of the
          ! crystal and aggregate fall speeds.
          !-----------------------------------------------
      mixratio_fromabove(i) = (snow_agg(i) + snow_cry(i))               &
                               * dhi(i)*rhor(i)
      IF ((qcf_cry(i)+qcf_agg(i)+mixratio_fromabove(i)) >  m0) THEN

            !-----------------------------------------------
            ! Make a linear combination of fall speeds between this
            ! layer and the layer above to aid numerical solution
            !-----------------------------------------------
        fqi_cry(i)= ( fqi_cry(i)*(qcf_cry(i)+qcf_agg(i))                &
                      + vf_agg(i)*mixratio_fromabove(i) )               &
                    / (qcf_cry(i)+qcf_agg(i)+mixratio_fromabove(i))
        fqi_agg(i)= ( fqi_agg(i)*(qcf_cry(i)+qcf_agg(i))                &
                      + vf_agg(i)*mixratio_fromabove(i) )               &
                    / (qcf_cry(i)+qcf_agg(i)+mixratio_fromabove(i))

            !-----------------------------------------------
            ! Fall speed of ice to pass to layer below
            !-----------------------------------------------
        vf_agg(i) = fqi_cry(i) * (1.0-frac_agg(i))                      &
                  + fqi_agg(i) *      frac_agg(i)

      ELSE  ! qcf gt m0
            ! Little ice so set fall speed to zero
        vf_agg(i)=0.0

      END IF  ! qcf gt m0

          !-----------------------------------------------
          ! Solve for fraction of ice that remains in the same layer
          !-----------------------------------------------
      fqirqi2_agg(i) = EXP(-fqi_agg(i)*dhi(i))
      fqirqi2_cry(i) = EXP(-fqi_cry(i)*dhi(i))

      IF (fqi_agg(i)  >   0.0) THEN
            !-----------------------------------------------
            ! Advect aggregates
            !-----------------------------------------------
        fqirqi_agg = snow_agg(i) + dhir(i) *                            &
                    (rho(i)*qcf_agg(i)-snow_agg(i)/fqi_agg(i))          &
                    * (1.0-fqirqi2_agg(i))
        qcf_agg(i) = snow_agg(i) * rhor(i) / fqi_agg(i)                 &
                    * (1.0-fqirqi2_agg(i))                              &
                    + qcf_agg(i)*fqirqi2_agg(i)
      ELSE
            !-----------------------------------------------
            ! No fall of ice out of the layer
            !-----------------------------------------------
        fqirqi_agg = 0.0
        qcf_agg(i) = snow_agg(i)*rhor(i)*dhi(i)
      END IF  ! fqi_agg gt 0

      IF (fqi_cry(i)  >   0.0 .AND. l_crystals) THEN
            !-----------------------------------------------
            ! Advect crystals
            !-----------------------------------------------
        fqirqi_cry = snow_cry(i) + dhir(i) *                            &
                    (rho(i)*qcf_cry(i)-snow_cry(i)/fqi_cry(i))          &
                    * (1.0-fqirqi2_cry(i))
        qcf_cry(i) = snow_cry(i) * rhor(i) / fqi_cry(i)                 &
                    * (1.0-fqirqi2_cry(i))                              &
                    + qcf_cry(i)*fqirqi2_cry(i)
      ELSE
            !-----------------------------------------------
            ! No fall of ice out of the layer
            !-----------------------------------------------
        fqirqi_cry = 0.0
        qcf_cry(i) = snow_cry(i)*rhor(i)*dhi(i)
      END IF  ! fqi_cry gt 0

          ! --------------------------------------------------
          ! Snow is used to save flux out of layer
          ! --------------------------------------------------
      snowt_cry(i) = snowt_cry(i) + fqirqi_cry/iterations
      snowt_agg(i) = snowt_agg(i) + fqirqi_agg/iterations

    END DO  ! Points

  ELSE  ! l_mcr_qcf2
        !--------------------------------------------------
        ! Prognostic Ice Crystal and Snow Aggregate Sedimentation
        !--------------------------------------------------
        ! Ice Crystals
    IF (l_crystals) THEN
      CALL lsp_sedim_eulexp(                                            &
        iterations, points, m0, dhi, dhir, rho, rhor,                   &
        snow_cry, fqi_cry, qcf_cry, vf_cry,                             &
        snowt_cry)
    END IF  ! l_crystals

        ! Aggregates
    CALL lsp_sedim_eulexp(                                              &
      iterations, points, m0, dhi, dhir, rho, rhor,                     &
      snow_agg, fqi_agg, qcf_agg, vf_agg,                               &
      snowt_agg)

  END IF ! l_mcr_qcf2

      !--------------------------------------------------
      ! Rain Sedimentation
      !--------------------------------------------------
  IF (l_mcr_qrain) THEN

    DO i = 1, points

      IF (qrain(i)  >   m0) THEN
            !--------------------------------------------------
            ! Estimate fall speed out of this layer (FQI_RAIN)
            !--------------------------------------------------

            !Calculate lambda for rain:
            ! assuming rain is a mixing ratio (kg kg-1)

        IF (l_rainfall_as) THEN

          IF (l_warm_new) THEN
            ! flux should be scaled by rain fraction
            lamr(i) = (1.0/(qrain(i)*rho(i)*constp(50)/rainfrac(i)))    &
                                                             **cx(52)
          ELSE
            ! original (wrong) method which ignored rain fraction
            lamr(i) = (1.0/(qrain(i)*rho(i)*constp(50)))**cx(52)
          END IF

              !Calculate additional properties

              ! lamda +h1r
          lamrh1(i) = lamr(i)+cx(56)
              ! lamda +h2r
          lamrh2(i) = lamr(i)+cx(57)

          fqi_rain(i) = ((lamr(i)**(cx(45)))/constp(53))*corr(i)*       &
             ( (constp(54) / lamrh1(i)**cx(59) ) +                      &
               (constp(55) / lamrh2(i)**cx(60) ) )

          IF (l_warm_new .AND. fqi_rain(i) < 0.001) THEN
            ! minimum fall speed allowed for rain
            fqi_rain(i) = 0.001
          END IF

        ELSE

          IF (l_warm_new) THEN
            ! flux should be scaled by rain fraction
            fqi_rain(i) = constp(41) * corr(i) / 6.0 *                  &
            (rho(i) * qrain(i) * constp(50)/rainfrac(i)) ** cx(51)
          ELSE
            ! original (wrong) method which ignored rain fraction
            fqi_rain(i) = constp(41) * corr(i) / 6.0 *                  &
            (rho(i) * qrain(i) * constp(50)) ** cx(51)
          END IF

        END IF !L_rainfall_as

      ELSE
            ! Fall speed is set to zero
        fqi_rain(i) = 0.0
      END IF  ! qrain gt m0

    END DO  ! Points

    CALL lsp_sedim_eulexp(                                              &
      iterations, points, m0, dhi, dhir, rho, rhor,                     &
      rainrate, fqi_rain, qrain, vf_rain,                               &
      rainratet)

  END IF  ! l_mcr_qrain

      !--------------------------------------------------
      ! Graupel Sedimentation
      !--------------------------------------------------
  IF (l_mcr_qgraup) THEN

    DO i = 1, points

      IF (qgraup(i)  >   m0) THEN
            !--------------------------------------------------
            ! Estimate fall speed out of this layer (FQI_GRAUP)
            !--------------------------------------------------
        fqi_graup(i) = constp(64) * corr(i) *                           &
          (rho(i) * qgraup(i) * constp(65)) ** cx(63)
      ELSE
            ! Fall speed is set to zero
        fqi_graup(i) = 0.0
      END IF  ! qgraup gt m0

    END DO  ! Points

    CALL lsp_sedim_eulexp(                                              &
       iterations, points, m0, dhi, dhir, rho, rhor,                    &
       grauprate, fqi_graup, qgraup, vf_graup,                          &
       graupratet)

  END IF  !  l_mcr_qgraup

  DO i = 1, points

      !-----------------------------------------------
      ! Update cloud fractions
      !-----------------------------------------------

      IF ((qcf_cry(i)+qcf_agg(i))  >   0.0) THEN

        IF (pc2_falliceshear_method == original_but_wrong) THEN

          !----------------------------------------------------------
          ! Calculate fraction of a layer the ice has fallen
          !----------------------------------------------------------
          ! This incorrectly uses the fall velocities in THIS layer
          ! rather than those from the layer ABOVE.
          !----------------------------------------------------------
          IF (l_mcr_qcf2) THEN
            temp3(i) = dhi(i) * (vf_agg(i)*qcf_agg(i)                   &
                       + vf_cry(i)*qcf_cry(i))/(qcf_cry(i)+qcf_agg(i))
          ELSE
            temp3(i) = dhi(i) * vf_agg(i)
          END IF
          !----------------------------------------------------------
          ! Calculate the amount of cloud overhang between levels
          !----------------------------------------------------------
          overhang = MAX(frac_ice_above(i)-cff(i),0.0)
          IF (temp3(i)  >   0.0) THEN
            overhang = overhang + wind_shear_factor / temp3(i)          &
                       * timestep
          END IF

          !----------------------------------------------------------
          ! Physially limit the amount of overhang. Note there is
          ! no limit on the fraction of ice in the layer above
          !----------------------------------------------------------
          overhang = MIN(overhang,1.0-cff(i))
          ! Now limit the fall out quantity
          temp3(i) = MIN(MAX(temp3(i),0.0),1.0)

          !----------------------------------------------------------
          ! Calculate change in ice cloud fraction
          !----------------------------------------------------------
            
          temp3(i) = temp3(i) * overhang
          deltacff(i) = temp3(i)

          IF (temp3(i)  <=  0.0) THEN
            !--------------------------------------------------------
            ! Total cloud fraction will be reduced
            !--------------------------------------------------------
            deltacf(i) = temp3(i)*area_ice(i)*cficei(i) !Random o'lap
          ELSE IF (cfice(i)  <   1.0) THEN
            !--------------------------------------------------------
            ! Total cloud fraction will be increased
            !--------------------------------------------------------
!           deltacf(i) = temp3(i)*area_clear(i) /(1.0-cfice(i))
!                                                      !Random o'lap
            deltacf(i) = (MIN(temp3(i),area_clear(i))) !Minimum o'lap
          END IF

        ELSE ! (pc2_falliceshear_method /= original_but_wrong)

          !--------------------------------------------------------
          ! Calculate temp3 or "fraction ice fallen":  
          ! fraction of depth of this layer the ice from the layer  
          ! above has fallen. Using fall speed INTO layer. 
          !--------------------------------------------------------

          IF (l_mcr_qcf2) THEN
            ! Find mass-weighted fall velocity (mwfv) for combined 
            ! ice categories.
            mwfv = (save_vf_agg(i)*qcf_agg(i)+save_vf_cry(i)*qcf_cry(i))&
                 / (qcf_cry(i)+qcf_agg(i))
          ELSE
            ! Using single ice category
            mwfv = save_vf_agg(i)
          END IF

          temp3(i) = mwfv * dhi(i) 

          ! Ensure temp3 is positive
          ! but allow "fraction fallen" to be > 1. 
          temp3(i)            = MAX(temp3(i),0.0)

          !------------------------------------------------------
          ! Calculate the amount of cloud overhang between levels
          !------------------------------------------------------
          overhang = MAX(frac_ice_above(i)-cff(i),0.0)

          IF (pc2_falliceshear_method == real_shear) THEN
            ! Increase the overhang depending on the vertical
            ! shear of the model wind.

            ! Magnitude of vertical shear of the horizontal wind.
            ! |dU/dz| = SQRT( dudz^2 + dvdz^2 )
            dudz = ( ukp1(i) - uk(i) )
            dvdz = ( vkp1(i) - vk(i) )
            shear = SQRT( (dudz*dudz) + (dvdz*dvdz) )

            ! The horizontal scale is taken as the square root 
            ! of the area of the grid box.

            horiz_scale = SQRT (   r_theta_levels_c(i) * mp_dell         &
                                 * r_theta_levels_c(i) * mp_delp         &
                                 * fv_cos_theta_latitude_c(i)     )

            ! Calculate the horizontal distance (in metres) the ice
            ! has moved across
            lateral_disp = shear * timestep

            ! Convert the lateral displacement of the falling ice
            ! cloud fraction to an increase in ice cloud fraction
            ! overhang by considering the size of the grid-box.
            overhang = overhang + ( lateral_disp / horiz_scale )

          END IF ! pc2_falliceshear_method == real shear

          !-----------------------------------------------
          ! Calculate change in ice cloud fraction
          !-----------------------------------------------
          ! The overhanging cloud gets advected down a 
          ! certain fraction of the depth of the layer. Now assume the 
          ! cloud fills the whole depth of the layer and
          ! reduce the lateral extent while conserving cloud volume.

          deltacff(i)=min(temp3(i) * overhang, 1.0-cff(i))

          ! Augment the change in ice cloud fraction to account
          ! for the lateral spreading out of ice cloud (e.g. cirrus).
          ! This will increase CFF while keeping IWC the same.
          !
          ! Cloud can only spread out from its edges, so work out the 
          ! perimeter of the cloud edge as a function of cloud fraction.
          cff_perimeter=-(2.0*cff(i)*cff(i))+(2.0*cff(i))
          deltacff(i)=deltacff(i)+(cff_spread_rate*cff_perimeter*timestep)
          deltacff(i)=min(deltacff(i), 1.0-cff(i)) 

          IF (cfice(i) < 1.0) THEN
            !-----------------------------------------------
            ! Total cloud fraction will be increased
            !-----------------------------------------------
            deltacf(i) = (MIN(deltacff(i),area_clear(i))) !Minimum o'lap
          ELSE IF (cfice(i) == 1.0) THEN 
            deltacf(i) = 0.0
          END IF

        END IF ! pc2_falliceshear_method

      ELSE  ! qcf gt 0
        ! Set ice cloud fraction to zero and total cloud
        ! fraction to the liquid cloud fraction
        deltacff(i) = -cff(i)
        deltacf(i)  = (cfl(i) - cf(i))
      END IF  ! qcf gt 0


        !-----------------------------------------------
        ! Form transfer diagnostics
        !-----------------------------------------------

    d_qcf_cry_dt(i) = d_qcf_cry_dt(i) +                                 &
      (qcf_cry(i) - qcf_cry_0(i)) * one_over_tsi
    d_qcf_agg_dt(i) = d_qcf_agg_dt(i) +                                 &
      (qcf_agg(i) - qcf_agg_0(i)) * one_over_tsi
    d_qrain_dt(i) = d_qrain_dt(i) +                                     &
      (qrain(i) - qrain_0(i)) * one_over_tsi
    d_qgraup_dt(i) = d_qgraup_dt(i) +                                   &
      (qgraup(i) - qgraup_0(i)) * one_over_tsi

        !-----------------------------------------------
        ! Update water contents
        !-----------------------------------------------

      qcf_cry_0(i) = qcf_cry(i)
      qcf_agg_0(i) = qcf_agg(i)
      qrain_0(i)   = qrain(i)
      qgraup_0(i)  = qgraup(i)

        !-----------------------------------------------
        ! Update cloud fractions
        !-----------------------------------------------

      cf(i)  = cf(i)  + deltacf(i)
      cff(i) = cff(i) + deltacff(i)

      cftransfer(i)  = cftransfer(i)  + deltacf(i)  * one_over_tsi
      cfftransfer(i) = cfftransfer(i) + deltacff(i) * one_over_tsi

  END DO  ! Points

  IF (lhook) CALL dr_hook('LSP_FALL',zhook_out,zhook_handle)

  RETURN
END SUBROUTINE lsp_fall
END MODULE lsp_fall_mod
