! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Iterative melting interface
MODULE lsp_it_fall_melt_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_it_fall_melt(                                            &
  points, timestep,                                                     &
                                          ! Number of points and tstep
  t, p, q, q_ice, qsl,                                                  &
                                          ! Temperature and vapour
  qcf_cry, qcf_agg, frac_agg,                                           &
                                          ! Ice contents
  qrain, qgraup,                                                        &
                                          ! Rain and graupel contents
  snow_agg,snow_cry,rainrate,grauprate,                                 &
                                             ! Sedimentation into layer
  snowt_agg,snowt_cry,rainratet,graupratet,                             &
                                                 ! Sedim. out of layer
  vf_agg, vf_cry, vf_rain, vf_graup,                                    &
                                          ! Fall speeds of hydrometeors
  area_liq, area_mix, area_clear,                                       &
                                          ! Cloud fraction information
  area_ice,cfice,cficei,                                                &
                                          ! at start of microphy. ts
  frac_ice_above,                                                       &
  cf, cfl, cff,                                                         &
                                          ! Current cloud fractions for
                                          ! updating
  rainfrac, rain_liq, rain_mix,                                         &
                                          ! Rain fraction information
  rain_ice, rain_clear,                                                 &
  rho, rhor,                                                            &
  tcg,tcgi,tcgc,tcgci,tcgg,tcggi,                                       &
                                          ! Parametrization information
  corr, corr2, rocor, dhi, dhir,                                        &
                                          ! Parametrization information
  lfrcp,                                                                &
                                          ! Microphysical information
  l_psd,                                                                &

  pifall, psfall, prfall, pgfall,                                       &
                                          ! Mass transfer diagnostics
  pimlt, psmlt, pgmlt,                                                  &
                                          ! Mass transfer diagnostics
  iterations, n_melt_iterations,                                        &
                                          ! Iterations of microphysics
  cftransfer,cfftransfer,rftransfer,                                    &
                                          ! Cloud transfer diagnostics
  uk, vk, ukp1, vkp1,                                                   &
                                          ! Winds for calc of wind-shear
  r_theta_levels_c, fv_cos_theta_latitude_c                             &
  )

  ! Microphysics Modules
  USE mphys_ice_mod,        ONLY: m0
  USE mphys_bypass_mod,     ONLY: l_crystals
  USE mphys_inputs_mod,     ONLY: l_mcr_qgraup
  
  ! Dr Hook Modules
  USE yomhook,              ONLY: lhook, dr_hook
  USE parkind1,             ONLY: jprb, jpim

  ! Large scale precipitation modules
  USE lsp_fall_mod,    ONLY: lsp_fall
  USE lsp_melting_mod, ONLY: lsp_melting
  IMPLICIT NONE

! Purpose:
!   Call fall out and melting terms in an iterative way

! Method:
!   Iterate over the fall-out and melting terms for a prescribed
!   number of iterations.
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
                          ! Number of points to calculate
    iterations,                                                         &
                          ! Number of microphysics iterations
    n_melt_iterations
                          ! Number of iterative melting iterations

  REAL, INTENT(IN) ::                                                   &
    timestep,                                                           &
                          ! Timestep / s
    q(points),                                                          &
                          ! Vapour content / kg kg-1
    q_ice(points),                                                      &
                          ! Vapour content in ice only region / kg kg-1
    qsl(points),                                                        &
                          ! Sat. vapour content wrt liquid / kg kg-1
    p(points),                                                          &
                          ! Pressure / N m-2
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
    area_liq(points),                                                   &
                          ! Frac of gridbox with liq cloud but not ice
    area_mix(points),                                                   &
                          ! Fraction of gridbox with ice and liq cloud
    area_clear(points),                                                 &
                          ! Fraction of gridbox with clear sky
    area_ice(points),                                                   &
                          ! Frac of gridbox with ice cloud but not liq
    cfice(points),                                                      &
                          ! Fraction of gridbox with ice cloud
    cficei(points),                                                     &
                          ! 1/Fraction of gridbox with ice cloud
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
!    &, m0                ! Seed ice water / kg kg-1 (in c_lspmic)
!    &, wind_shear_factor ! Degree of overhang of ice fraction between
                          ! two model layers (in c_lspmic)
    tcg(points),                                                        &
                          ! T dependent func in ice agg. size dist'n
    tcgi(points),                                                       &
                          ! 1/tcg (no units)
    tcgc(points),                                                       &
                          ! T dependent func in ice crystal size dist'n
    tcgci(points),                                                      &
                          ! 1/tcgc (no units)
    tcgg(points),                                                       &
                          ! T dependent func in graupel size dist'n
    tcggi(points),                                                      &
                          ! 1/tcgg (no units)
    corr(points),                                                       &
                          ! Air density fall speed correction (no units)
    corr2(points),                                                      &
                          ! Temperature correcn to diffusivity (no units
    rocor(points),                                                      &
                          ! sqrt(rho*corr*corr2)
    uk(points),                                                         &
                          ! U wind on level k
    vk(points),                                                         &
                          ! V wind on level k
    ukp1(points),                                                       &
                          ! U wind on level k+1
    vkp1(points),                                                       &
                          ! V wind on level k+1
    r_theta_levels_c(points),                                           &
                          ! Distance from centre of Earth
    fv_cos_theta_latitude_c(points),                                    &
                          ! Finite vol cos of latitude
    lfrcp
                          ! Latent heat of fusion
                          !   / heat capacity of air / K

  REAL, INTENT(INOUT) ::                                                &
    t(points),                                                          &
                          ! Temperature / K
    qcf_cry(points),                                                    &
                          ! Ice crystal mixing ratio / kg kg-1
    qcf_agg(points),                                                    &
                          ! Ice aggregate mixing ratio / kg kg-1
    qrain(points),                                                      &
                          ! Rain mixing ratio / kg kg-1
    qgraup(points),                                                     &
                          ! Graupel mixing ratio / kg kg-1
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
    cff(points)           ! Current ice cloud fraction

  REAL, INTENT(INOUT) ::                                                &
    rainfrac(points),                                                   &
                          ! Fraction of gridbox with rain
    rain_liq(points),                                                   &
                          ! Fraction of gridbox with rain and liq cloud
    rain_mix(points),                                                   &
                          ! Frac of gbox with rain and mixed phase cloud
    rain_ice(points),                                                   &
                          ! Fraction of gridbox with rain and ice cloud
    rain_clear(points),                                                 &
                          ! Fraction of gridbox with rain but no cloud
    pifall(points),                                                     &
                          ! Rate of change of crystal, aggregate,
    psfall(points),                                                     &
                          ! rain and graupel mixing ratios due
    prfall(points),                                                     &
                          ! to sedimentation / kg kg-1 s-1
    pgfall(points),                                                     &
    pimlt(points),                                                      &
                          ! Melting rate of crystals / kg kg-1 s-1
    psmlt(points),                                                      &
                          ! Melting rate of snow / kg kg-1 s-1
    pgmlt(points),                                                      &
                          ! Melting rate of graupel / kg kg-1 s-1
    cftransfer(points),                                                 &
                          ! Cloud fraction increment this tstep
    cfftransfer(points),                                                &
                           ! Ice cloud fraction inc this tstep
    rftransfer(points)! Rain fraction increment this tstep

  LOGICAL, INTENT(IN) ::                                                &
    l_psd
                          ! Use generic ice particle size distribution

! Local Variables

  INTEGER ::                                                            &
    i,                                                                  &
                          ! Loop counter for points
    k,                                                                  &
                          ! Loop counter for melting iterations
    total_iterations  ! Net number of iterations

  REAL ::                                                               &
    qcft(points),                                                       &
                          ! Total ice content / kg kg-1
    dhi_it(points),                                                     &
                          ! Iterated tstep / layer thick. / s m-1
    dhir_it(points),                                                    &
                          ! Layer thick. / iterated tstep / m s-1
    timestep_it,                                                        &
                          ! Iterated timestep / s
    one_over_tsi          ! 1/(timestep*iterations)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


      !-----------------------------------------------
      ! Set up short timestep
      !-----------------------------------------------

  IF (lhook) CALL dr_hook('LSP_IT_FALL_MELT',zhook_in,zhook_handle)

  timestep_it = timestep / n_melt_iterations
  total_iterations = iterations * n_melt_iterations
  one_over_tsi = 1.0/(timestep_it*total_iterations)

  DO i = 1, points
        ! Update dhi values for short timestep
    dhi_it(i) = dhi(i) / n_melt_iterations
    dhir_it(i) = dhir(i) * n_melt_iterations
  END DO

      !-----------------------------------------------
      ! Call fall and melting for the iterated points
      !-----------------------------------------------
  DO k = 1, n_melt_iterations

    CALL lsp_fall(points, timestep_it,                                  &
                     qcf_cry, qcf_agg, frac_agg, qrain, qgraup, t,      &
                     snow_agg, snow_cry, rainrate, grauprate,           &
                     snowt_agg, snowt_cry, rainratet, graupratet,       &
                     vf_agg, vf_cry, vf_rain, vf_graup,                 &
                     area_clear, area_ice, cfice, cficei,               &
                     frac_ice_above, cf, cfl, cff,                      &
                     rho, rhor, tcgi, tcgci,                            &
                     corr, dhi_it, dhir_it, rainfrac,                   &
                     pifall, psfall, prfall, pgfall,                    &
                     total_iterations, one_over_tsi,                    &
                     cftransfer, cfftransfer,                           &
                     uk, vk, ukp1, vkp1,                                &
                     r_theta_levels_c, fv_cos_theta_latitude_c          &
                    )

    DO i = 1, points
          ! Calculate total ice content
      qcft(i) = qcf_cry(i) + qcf_agg(i)
    END DO

    IF (l_crystals) THEN
      CALL lsp_melting(points, timestep_it,                             &
                      q, q_ice, qcf_cry, qcft, qrain, qsl, t, p,        &
                      area_liq, area_mix, area_ice, area_clear,         &
                      cfice, cficei, frac_ice_above,                    &
                      rainfrac, rain_liq, rain_mix,                     &
                      rain_ice, rain_clear, cf, cff,                    &
                      rho, rhor, m0, tcg, tcgi, corr2, rocor,           &
                      lfrcp, 0,                                         &
                      .FALSE.,                                          &
                      pimlt, one_over_tsi,                              &
                      cftransfer, cfftransfer, rftransfer               &
                     )
    END IF  ! l_crystals

    DO i = 1, points
          ! Calculate total ice content
      qcft(i) = qcf_cry(i) + qcf_agg(i)
    END DO

    CALL lsp_melting(points, timestep_it,                               &
                      q, q_ice, qcf_agg, qcft, qrain, qsl, t, p,        &
                      area_liq, area_mix, area_ice, area_clear,         &
                      cfice, cficei, frac_ice_above,                    &
                      rainfrac, rain_liq, rain_mix,                     &
                      rain_ice, rain_clear, cf, cff,                    &
                      rho, rhor, m0, tcgc, tcgci, corr2, rocor,         &
                      lfrcp, 1,                                         &
                      l_psd,                                            &
                      psmlt, one_over_tsi,                              &
                      cftransfer, cfftransfer, rftransfer               &
                     )

    IF (l_mcr_qgraup) THEN
          ! Graupel does not update cloud fractions so there is no need
          ! to update qcft (it is not used)
      CALL lsp_melting(points, timestep_it,                             &
                        q, q_ice, qgraup, qcft, qrain, qsl, t, p,       &
                        area_liq, area_mix, area_ice, area_clear,       &
                        cfice, cficei, frac_ice_above,                  &
                        rainfrac, rain_liq, rain_mix,                   &
                        rain_ice, rain_clear, cf, cff,                  &
                        rho, rhor, m0, tcgg, tcggi, corr2, rocor,       &
                        lfrcp, 3,                                       &
                        .FALSE.,                                        &
                        pgmlt, one_over_tsi,                            &
                        cftransfer, cfftransfer, rftransfer             &
                       )
    END IF  ! l_mcr_qgraup

  END DO  ! n_melt_iterations

  IF (lhook) CALL dr_hook('LSP_IT_FALL_MELT',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_it_fall_melt
END MODULE lsp_it_fall_melt_mod
