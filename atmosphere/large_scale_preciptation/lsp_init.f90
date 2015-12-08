! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Initialisation of variables
! Subroutine Interface:
MODULE lsp_init_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_init(                                                    &
  points,                                                               &
                                     ! Number of points
  timestepfixed, timestep,                                              &
                                     ! Timesteps
  t, p, cttemp,                                                         &
                                     ! Temperature and pressure
  deltaz,rhodz,rhodz_dry,rhodz_moist,                                   &
                                           ! Air density information
  rho, rhor,                                                            &
  dhi, dhir, dhilsiterr,                                                &
                                     ! Tstep and layer thick. ratios
  q, qcl, qcf, qcf2, qrain, qgraup, rainrate,                           &
                                                   ! Water contents
  qcf_agg, qcf_cry, qcf_tot, frac_agg, frac_cry_dep, frac_agg_dep,      &
  qs, qsl, esi, esw,                                                    &
                                     ! Saturation water quantities
  psdep, psaut, psacw, psacr,                                           &
                                     ! Aggregate transfer rates
  psaci, psmlt, psmltevp,psfall,                                        &
  pifrw, pifrr, piprm, pidep, piacw,                                    &
                                     ! Ice transfer rates
  piacr, pimlt, pimltevp,pifall,                                        &
  praut, pracw, prevp, prfall,                                          &
                                     ! Rain transfer rates
  plset, plevpset,                                                      &
                                     ! Droplet settling transfers
  pgaut, pgacw, pgacs, pgmlt,                                           &
                                     ! Graupel transfer rates
  pgfall,                                                               &
  cf_transfer_diag, cfl_transfer_diag,                                  &
  cff_transfer_diag, rf_transfer_diag,                                  &
  snow_cry, snow_agg, snowt_cry,                                        &
                                     ! Precipitation rates
  snowt_agg, rainratet, graupratet,                                     &
  lheat_correc_liq, lheat_correc_ice,                                   &
                                           ! Latent heat corrections
  corr, corr2, rocor,                                                   &
                                     ! Fall speed and diffusivity
                                     ! corrections
  tcg,tcgi,tcgc,tcgci,tcgg,tcggi,                                       &
                                           ! Ice PSD intercepts
  lsiter, niter_bs )


  ! Microphysics modules
  USE mphys_ice_mod,       ONLY: t_scaling, m0, qcf0, t_agg_min
  USE mphys_constants_mod, ONLY: cx
  USE mphys_inputs_mod,    ONLY: l_cry_agg_dep, l_psd, l_mcr_qrain,     &
                                 l_mcr_qcf2, l_mcr_qgraup
  USE mphys_bypass_mod,    ONLY: l_crystals

  ! Mathematical modules
  USE vectlib_mod,   ONLY: powr_v, exp_v

  ! General and Atmospheric Modules
  USE conversions_mod,      ONLY: zerodegc
  USE water_constants_mod,  ONLY: lc, lf
  USE atmos_constants_mod,  ONLY:                                       &
                            cp, r, repsilon, c_virtual, recip_epsilon
  USE um_input_control_mod, ONLY: l_mr_physics1

  ! Dr Hook Modules
  USE yomhook,             ONLY: lhook, dr_hook
  USE parkind1,            ONLY: jprb, jpim

  IMPLICIT NONE

! Purpose:
!   Perform initialization of variables required for the
!   microphysics scheme.

! Method:
!   Set variables to their initial values.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation


! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.



! Subroutine Arguments


  INTEGER, INTENT(IN) ::                                                &
    points,                                                             &
                       ! Number of points
    lsiter,                                                             &
                       ! Iterations of microphysics (inside column)
    niter_bs
                       ! Iterations of microphysics (outside column)

  REAL, INTENT(IN) ::                                                   &
    timestepfixed,                                                      &
                       ! Physics timestep / s
    t(points),                                                          &
                       ! Temperature / K
    p(points),                                                          &
                       ! Pressure / N m-2
    q(points),                                                          &
                       ! Vapour content / kg kg-1
    qcl(points),                                                        &
                       ! Liquid water content / kg kg-1
    qcf(points),                                                        &
                       ! Ice aggregate content / kg kg-1
    qcf2(points),                                                       &
                       ! Ice crystal content / kg kg-1
    cttemp(points),                                                     &
                       ! Cloud top temperature / K
    rainrate(points),                                                   &
                         ! Rainrate into layer / kg m-2 s-1
    deltaz(points),                                                     &
                       ! Layer thickness / m
    rhodz(points),                                                      &
                       ! Air density * deltaz based on hydrostatic
                       ! assumption / kg m-2
    rhodz_dry(points),                                                  &
                           ! Air density*deltaz for dry air based on
                           ! non-hydrostatic assumption / kg m-2
    rhodz_moist(points)
                           ! Air density*deltaz for moist air based on


  REAL, INTENT(INOUT) ::                                                &
    qgraup(points),                                                     &
                       ! Graupel content / kg kg-1
    snow_cry(points),                                                   &
                          ! Ice crystal precip into layer / kg m-2 s-1
    snow_agg(points)  ! Ice agg. precip into layer / kg m-2 s-1

  REAL, INTENT(INOUT) ::                                                &
        ! Microphysical process rate diagnostics / kg kg-1 s-1
    psdep(points),                                                      &
                       ! Deposition of vapour to snow aggregates
    psaut(points),                                                      &
                       ! Autoconversion of aggregates from crystals
    psacw(points),                                                      &
                       ! Accretion of liq. water by snow aggregates
    psacr(points),                                                      &
                       ! Collection of rain by snow aggregates
    psaci(points),                                                      &
                       ! Collection of ice crystals by aggregates
    psmlt(points),                                                      &
                       ! Melting of snow aggregates
    psmltevp(points),                                                   &
                          ! Evaporation of melting aggregates
    psfall(points),                                                     &
                       ! Fall-out of snow aggregates

    praut(points),                                                      &
                       ! Autoconversion of cloud drops to rain
    pracw(points),                                                      &
                       ! Accretion of liq. water by rain
    prevp(points),                                                      &
                       ! Evaporation of rain
    prfall(points),                                                     &
                       ! Fall-out of rain

    plset(points),                                                      &
                       ! Droplet settling of liquid water
    plevpset(points),                                                   &
                          ! Evaporated settled droplets

    pgaut(points),                                                      &
                       ! Autoconversion of graupel from aggregates
    pgacw(points),                                                      &
                       ! Accretion of liq. water by graupel
    pgacs(points),                                                      &
                       ! Collection of snow aggregates by graupel
    pgmlt(points),                                                      &
                       ! Melting of graupel
    pgfall(points),                                                     &
                       ! Fall-out of graupel

    pifrw(points),                                                      &
                       ! Homogeneous freezing nucleation
    pifrr(points),                                                      &
                       ! Homogeneous freezing nucleation of rain
    piprm(points),                                                      &
                       ! Heterogeneous (primary) nucleation
    pidep(points),                                                      &
                       ! Deposition of vapour to ice crystals
    piacw(points),                                                      &
                       ! Accretion of liq. water by ice crystals
    piacr(points),                                                      &
                       ! Collection of rain by ice crystals
    pimlt(points),                                                      &
                       ! Melting of ice crystals
    pimltevp(points),                                                   &
                          ! Evaporation of melting ice crystals
    pifall(points),                                                     &
                       ! Fall-out of ice crystals
    frac_agg(points)
                          ! Fraction of ice that is aggregates

  REAL, INTENT(OUT) ::                                                  &
    timestep,                                                           &
                          ! Timestep of each iteration / s
    snowt_cry(points),                                                  &
                          ! Ice crystal precip out of layer / kg m-2 s-1
    snowt_agg(points),                                                  &
                          ! Ice agg. precip out of layer / kg m-2 s-1
    rainratet(points),                                                  &
                          ! Rain precip rate out of layer / kg m-2 s-1
    graupratet(points),                                                 &
                          ! Graupel precip rate out of layer/ kg m-2 s-1
    qcf_agg(points),                                                    &
                          ! Ice aggregate mixing ratio / kg kg-1
    qcf_cry(points),                                                    &
                          ! Ice crystal mixing ratio / kg kg-1
    qcf_tot(points),                                                    &
                          ! Total ice content / kg kg-1
    frac_cry_dep(points),                                               &
                          ! Fraction of supersaturation that can be
                          ! removed by crystals
    frac_agg_dep(points),                                               &
                          ! Fraction of supersaturation that can be
                          !removed by aggregates
    qrain(points),                                                      &
                          ! Rain water content / kg kg-1
    qs(points),                                                         &
                          ! Saturated humidity wrt ice / kg kg-1
    qsl(points),                                                        &
                          ! Saturated humidity wrt liquid / kg kg-1
    rho(points),                                                        &
                       ! Air density / kg m-3
    rhor(points),                                                       &
                       ! 1 / air density / m kg-1
    esi(points),                                                        &
                       ! Vapour pressure wrt ice / N m-2
    esw(points),                                                        &
                       ! Vapour pressure wrt liquid / N m-2
    lheat_correc_liq(points),                                           &
                                 ! Latent heat correction for
                                 ! liquid (no units)
    lheat_correc_ice(points),                                           &
                                 ! Latent heat correction for
                                 ! ice (no units)
    dhi(points),                                                        &
                          ! Timestep / deltaz / s m-1
    dhir(points),                                                       &
                          ! deltaz / timestep / m s-1
    dhilsiterr(points),                                                 &
                          ! deltaz / (timestep*lsiter) / m s-1
    corr(points),                                                       &
                          ! Fall speed correction factor (no units)
    corr2(points),                                                      &
                          ! Diffusivity correction factor (no units)
    rocor(points),                                                      &
                          ! sqrt(rho*corr*corr2) (no units)
    tcg(points),                                                        &
                          ! Temperature dependent aggregate PSD
                          ! intercept factor (no units)
    tcgi(points),                                                       &
                          ! 1 / tcg (no units)
    tcgc(points),                                                       &
                          ! Temperature dependent crystal PSD
                          ! intercept factor (no units)
    tcgci(points),                                                      &
                          ! 1 / tcgc (no units)
    tcgg(points),                                                       &
                          ! Temperature dependent graupel PSD
                          ! intercept factor (no units)
    tcggi(points),                                                      &
                          ! 1 / tcgg (no units)
    cf_transfer_diag(points),                                           &
                          ! Diagnostic for change in cloud frac
    cfl_transfer_diag(points),                                          &
                          ! Diagnostic for change in liq cloud frac
    cff_transfer_diag(points),                                          &
                          ! Diagnostic for change in ice cloud frac
    rf_transfer_diag(points)
                          ! Diagnostic for change in rain fraction

! Local Variables

  INTEGER ::                                                            &
    i              ! Loop counter

  REAL ::                                                               &
    qc_tot(points) ! Total condensate / kg kg-1

  REAL, PARAMETER ::                                                    &
    one_over_zerodegc = 1.0/zerodegc,                                   &
    rho1 = 1.0
                                     ! Ref. air density / kg m-3

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


      !-----------------------------------------------
      ! Set up short timestep
      !-----------------------------------------------

  IF (lhook) CALL dr_hook('LSP_INIT',zhook_in,zhook_handle)

  timestep = timestepfixed/(lsiter*niter_bs)

  DO i = 1, points

      !-----------------------------------------------
      ! Set fluxes to zero
      !-----------------------------------------------
    snowt_cry(i)  = 0.0
    snowt_agg(i)  = 0.0
    rainratet(i)  = 0.0
    graupratet(i) = 0.0

      !-----------------------------------------------
      ! Initialize transfer diagnostics to zero
      !-----------------------------------------------
    pifrw(i) = 0.0
    pifrr(i) = 0.0
    piprm(i) = 0.0
    pidep(i) = 0.0
    piacw(i) = 0.0
    piacr(i) = 0.0
    pimlt(i) = 0.0
    pimltevp(i) = 0.0
    pifall(i)= 0.0

    psdep(i) = 0.0
    psaut(i) = 0.0
    psacw(i) = 0.0
    psacr(i) = 0.0
    psaci(i) = 0.0
    psmlt(i) = 0.0
    psmltevp(i) = 0.0
    psfall(i)= 0.0

    praut(i) = 0.0
    pracw(i) = 0.0
    prevp(i) = 0.0
    prfall(i)= 0.0

    plset(i) = 0.0
    plevpset(i) = 0.0

    pgaut(i) = 0.0
    pgacw(i) = 0.0
    pgacs(i) = 0.0
    pgmlt(i) = 0.0
    pgfall(i)= 0.0

    cf_transfer_diag(i) = 0.0
    cfl_transfer_diag(i)= 0.0
    cff_transfer_diag(i)= 0.0
    rf_transfer_diag(i) = 0.0

  END DO

      !-----------------------------------------------
      ! Set mixing ratios of graupel to zero if not used.
      ! If not done then (with high optimisation) this can be any
      ! old rubbish from memory and you get d_qgraup_dt being
      ! calculated as nan-nan in later routines, which causes a crash
      !-----------------------------------------------
  IF (.NOT. l_mcr_qgraup) THEN

    DO i = 1, points

      qgraup(i)   = 0.0

    END DO

  END IF  ! .not. l_mcr_qgraup

  IF (l_mcr_qcf2) THEN
        ! If l_mcr_qcf2 is true then there are two ice/snow
        ! prognostics (qcf and qcf2), so copy them to
        ! qcf_cry and qcf_agg

    DO i = 1, points

      qcf_cry(i) = qcf2(i)
      qcf_agg(i) = qcf(i)

    END DO

  ELSE IF (.NOT. l_crystals) THEN
        ! All ice is placed in a single ice category (aggregates).

    DO i = 1, points

      frac_agg(i) = 1.0
      qcf_cry(i)  = 0.0
      qcf_agg(i)  = qcf(i)
      snow_cry(i) = 0.0

    END DO

  ELSE
        ! Split the one ice/snow
        ! prognostic (qcf) diagnostically into ice crystals
        ! (qcf_cry) and snow aggregates (qcf_agg)

    DO i = 1, points

      frac_agg(i) = -t_scaling*MAX((t(i)-cttemp(i)),0.0)                &
                    *MAX(qcf(i)*qcf0,0.0)
    END DO
    CALL exp_v(points,frac_agg,frac_agg)
    DO i = 1, points
      frac_agg(i) = MAX(1.0-frac_agg(i) , 0.0)
          ! Allocate ice content to crystals and aggregates
      qcf_cry(i) = qcf(i) * (1.0-frac_agg(i))
      qcf_agg(i) = qcf(i) * frac_agg(i)

          ! Assume falling snow is partitioned into crystals and
          ! aggregates. Snow_agg contains total snow on input
      snow_cry(i) = snow_agg(i) * (1.0-frac_agg(i))
      snow_agg(i) = snow_agg(i) * frac_agg(i)

    END DO

  END IF ! l_mcr_qcf2

      !-----------------------------------------------
      ! Calculate total ice content
      !-----------------------------------------------
  DO i = 1, points

    qcf_tot(i) = qcf_cry(i) + qcf_agg(i)

  END DO

      !-----------------------------------------------
      ! Calculate the fraction of ice in the crystals and
      ! aggregates category that is allowed to remove
      ! supersaturation
      !-----------------------------------------------
  IF (l_cry_agg_dep) THEN

    IF (l_mcr_qcf2) THEN

          ! Use the partition given by the prognostic ice categories.
          ! Here we add a small amount of ice to the crystal category
          ! to ensure no problems when the ice contents are zero.

      DO i = 1, points

        frac_cry_dep(i) = (qcf_cry(i)+m0)/(MAX(qcf_tot(i),0.0)+m0)
        frac_agg_dep(i) = qcf_agg(i)/(MAX(qcf_tot(i),0.0)+m0)

      END DO

    ELSE  ! l_mcr_qcf2
          ! Use the diagnostic partition function

      DO i = 1, points

        frac_cry_dep(i) = 1.0 - frac_agg(i)
        frac_agg_dep(i) = frac_agg(i)

      END DO

    END IF

  ELSE
        ! Set the fractions to 1 to maintain bit reproducibility

    DO i = 1, points

      frac_cry_dep(i)=1.0
      frac_agg_dep(i)=1.0

    END DO

  END IF

      !-----------------------------------------------
      ! If rain is a diagnostic, convert flux (kg m-2 s-1)
      ! to mass (kg kg-1)
      !-----------------------------------------------
  IF (.NOT. l_mcr_qrain) THEN

  ! Rain is a diagnostic quantity

    DO i = 1, points

      IF (rainrate(i)  >   0.0) THEN

        IF (l_mr_physics1) THEN

          ! Mixing ratio formulation
          qrain(i) = rainrate(i) * timestepfixed / rhodz_dry(i)

        ELSE ! l_mr_physics1

          ! Specific humidity formulation
          qrain(i) = rainrate(i) * timestepfixed / rhodz_moist(i)

        END IF  ! l_mr_physics1

      ELSE    ! rainrate > 0

        qrain(i) = 0.0

      END IF  ! rainrate > 0

    END DO ! points

  END IF  ! .not. l_mcr_qrain

      !-----------------------------------------------
      ! Calculate saturation specific humidities
      !-----------------------------------------------
      ! Qsat with respect to ice
! DEPENDS ON: qsat_mix
  CALL qsat_mix(qs,t,p,points,l_mr_physics1)

      ! Qsat with respect to liquid water
! DEPENDS ON: qsat_wat_mix
  CALL qsat_wat_mix(qsl,t,p,points,l_mr_physics1)

  DO i = 1, points

      !-----------------------------------------------
      ! Calculate saturation vapour pressures
      !-----------------------------------------------
    esi(i) = qs (i) * p(i) * recip_epsilon
    esw(i) = qsl(i) * p(i) * recip_epsilon

      !-----------------------------------------------
      ! Calculate density of air
      !-----------------------------------------------
    IF (l_mr_physics1) THEN

      ! rho is the dry density
        rho(i) = rhodz_dry(i) / deltaz(i)

    ELSE

      ! rho is the moist density
        rho(i) = rhodz_moist(i) / deltaz(i)

    END IF  ! l_mr_physics1

        ! Calculate the inverse of the air density
    rhor(i) = 1.0 / rho(i)

        !-----------------------------------------------
        ! Estimate latent heat correction to rate of evaporation
        !-----------------------------------------------
    lheat_correc_liq(i) = 1.0/(1.0+repsilon*lc**2*qsl(i)                &
                             /(cp*r*t(i)**2))
    lheat_correc_ice(i) = 1.0/(1.0+repsilon*(lc+lf)**2*qs(i)            &
                             /(cp*r*t(i)**2))

        !-----------------------------------------------
        ! Calculate CFL timestep divided by level separation
        !-----------------------------------------------

        ! Use the formulation based on the heights of the levels
    dhi(i) = timestep/deltaz(i)

        ! Calculate the inverse
    dhir(i)       = 1.0/dhi(i)
        ! Calculate the inverse for the long timestep
    dhilsiterr(i) = 1.0/(dhi(i)*lsiter)

  END DO ! points

      !-----------------------------------------------
      ! Correction factors due to air density and temperature
      !-----------------------------------------------

  DO i = 1, points

        ! Correction of fall speeds
    IF (l_mr_physics1) THEN

      corr(i) = rho1*deltaz(i) / rhodz_moist(i)

    ELSE

      corr(i) = rho1*rhor(i)

    END IF

        ! Correction factor in viscosity etc. due to temperature

    corr2(i) = (t(i) * one_over_zerodegc)

  END DO

  CALL powr_v( points, corr,  0.4, corr  )
  CALL powr_v( points, corr2, 1.5, corr2 )

  DO i = 1, points

    corr2(i) = corr2(i)*(393.0/(t(i)+120.0))

    ! Combined correction factor

    IF (l_mr_physics1) THEN

      rocor(i) = rhodz_moist(i) / deltaz(i) * corr(i) * corr2(i)

    ELSE

      rocor(i) = rho(i)*corr(i)*corr2(i)

    END IF

  END DO

  CALL powr_v( points, rocor, 0.5, rocor )

      !-----------------------------------------------
      ! Calculate ice particle size distributions
      !-----------------------------------------------

  DO i = 1, points

        ! Calculate a temperature factor for N0crystals
    tcgc(i) = -cx(12)*MAX(t(i)-zerodegc,t_agg_min)
        ! Calculate a temperature factor for N0aggregates
    tcg(i)  = -cx(32)*MAX(t(i)-zerodegc,t_agg_min)

        ! Define inverse of TCG values
    tcgci(i) = -tcgc(i)
    tcgi(i)  = -tcg(i)

  END DO

  CALL exp_v( points, tcg,   tcg   )
  CALL exp_v( points, tcgc,  tcgc  )
  CALL exp_v( points, tcgi,  tcgi  )
  CALL exp_v( points, tcgci, tcgci )

      !-----------------------------------------------
      ! Calculate graupel size distributions
      !-----------------------------------------------
  DO i = 1, points

    tcgg(i)  = 1.0
    tcggi(i) = 1.0/tcgg(i)

  END DO

  IF (lhook) CALL dr_hook('LSP_INIT',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_init
END MODULE lsp_init_mod
