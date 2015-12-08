! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Melting of ice particles
! Subroutine Interface:
MODULE lsp_melting_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_melting(                                                 &
  points, timestep,                                                     &
                                          ! Number of points and tstep
  q, q_ice, qcf, qcft, qrain, qsl,                                      &
                                          ! Water contents
  t, p,                                                                 &
                                          ! Temperature and pressure
  area_liq, area_mix, area_ice, area_clear,                             &
                                          ! Partition information
  cfice, cficei, frac_ice_above,                                        &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  rainfrac, rain_liq, rain_mix,                                         &
                                          ! Rain fraction information
  rain_ice, rain_clear,                                                 &
  cf, cff,                                                              &
                                          ! Current cloud fractions for
!                                         ! updating
  rho, rhor, m0, tcg, tcgi,                                             &
                                          ! Parametrization information
  corr2, rocor,                                                         &
  lfrcp, ice_type,                                                      &
                                          ! Microphysical information
  l_psd,                                                                &
                                          ! Code options
  ptransfer,                                                            &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                         &
                                          ! 1/(timestep*iterations)
  cftransfer,cfftransfer,                                               &
                                          ! Cloud transfer diagnostics
  rftransfer                                                            &
                                          ! Rain fraction transfer
  )

  ! Microphysics modules
  USE mphys_constants_mod, ONLY: cx, constp, ice_type_offset
  USE lsp_dif_mod,         ONLY: tw1, tw2, tw3, tw4, tw5

  ! General atmosphere modules
  USE conversions_mod,     ONLY: zerodegc
  
  ! Dr Hook Modules
  USE yomhook,             ONLY: lhook, dr_hook
  USE parkind1,            ONLY: jprb, jpim

  ! Large scale precipitation modules
  USE lsp_moments_mod, ONLY: lsp_moments
  
  IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of melting of ice particles

! Method:
!   Solve implicitly the microphysical transfer equation for a
!   specified distribution of ice particles melting at an
!   approximated wet-bulb temperature
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
    ice_type          ! Type of ice (0 - crystals, 1 - aggregates
                          !              3 - graupel)

  REAL, INTENT(IN) ::                                                   &
    timestep,                                                           &
                          ! Timestep / s
    q(points),                                                          &
                          ! Gridbox mean vapour content / kg kg-1
    q_ice(points),                                                      &
                          ! Vapour content in ice partition / kg kg-1
    qcft(points),                                                       &
                          ! Total ice water content / kg kg-1
    qsl(points),                                                        &
                          ! Saturated mixing ratio wrt liquid / kg kg-1
    p(points),                                                          &
                          ! Pressure / N m-2
    area_liq(points),                                                   &
                          ! Fraction of gridbox with liquid-only cloud
    area_mix(points),                                                   &
                          ! Fraction of gridbox with mixed phase cloud
    area_ice(points),                                                   &
                          ! Fraction of gridbox with ice-only cloud
    area_clear(points),                                                 &
                          ! Frac of gridbox with no cloud
    cfice(points),                                                      &
                          ! Fraction of gridbox with ice cloud
    cficei(points),                                                     &
                          ! 1/Fraction of gridbox with ice cloud
    frac_ice_above(points),                                             &
                               ! Ice fraction in layer above this layer
    rho(points),                                                        &
                          ! Air density / kg m-3
    rhor(points),                                                       &
                          ! 1/Air density / m3 kg-1
    corr2(points),                                                      &
                          ! Temperature dependency of diffusion params
    rocor(points),                                                      &
                          ! Combination of corr and corr2
    m0,                                                                 &
                          ! Seed ice water content / kg kg-1
    tcg(points),                                                        &
                          ! T dependent function in ice size dist'n
    tcgi(points),                                                       &
                          ! 1/tcg (no units)
    lfrcp,                                                              &
                          ! Latent heat of fusion
                          ! / heat capacity of air / K
    one_over_tsi          ! 1/(timestep*iterations)

  REAL, INTENT(INOUT) ::                                                &
    qcf(points),                                                        &
                          ! Ice water content of this category / kg kg-1
    qrain(points),                                                      &
                          ! Rain mixing ratio / kg kg-1
    t(points),                                                          &
                          ! Temperature / K
    cf(points),                                                         &
                          ! Current cloud fraction
    cff(points),                                                        &
                          ! Current ice cloud fraction
    rainfrac(points),                                                   &
                          ! Rain fraction
    rain_liq(points),                                                   &
                          ! Overlap fraction of rain and liquid
    rain_mix(points),                                                   &
                          ! Overlap fraction of rain and mixed phase
    rain_ice(points),                                                   &
                          ! Overlap fraction of rain and ice
    rain_clear(points),                                                 &
                          ! Overlap fraction of rain and clear sky
    ptransfer(points),                                                  &
                          ! Mass rimed in this timestep / kg kg-1
    rftransfer(points),                                                 &
                          ! Transfer of rain fraction this timestep
    cftransfer(points),                                                 &
                           ! Cumulative cloud fraction increment
    cfftransfer(points)! Cumulative liquid cloud fraction increment

  LOGICAL, INTENT(IN) ::                                                &
    l_psd
                          ! Use generic ice particle size distribution

! Local Variables

  INTEGER ::                                                            &
    i,                                                                  &
                          ! Loop counter for points
    cry_offset        ! Index offset for ice crystals

  REAL ::                                                               &
    tempw(points),                                                      &
                          ! Average supersaturation in gridbox
    temp7(points),                                                      &
                          ! Wet bulb temperature minus 0 C / deg C
    pr02(points),                                                       &
                          ! Temporary in melting mass calculation
    dpr(points),                                                        &
                          ! Mass melted this timestep / kg kg-1
    delta_cf(points),                                                   &
                          ! Cloud fraction transfered (no units)
    delta_cff(points),                                                  &
                          ! Ice cloud fraction transfered (no units)
    rf_final(points),                                                   &
                          ! Final rain fraction (no units)
    m_1(points),                                                        &
                          ! 1st moment of generic particle size dist.
    m_0p5_dp3(points)
                          ! 1+(di+1)/2 moment of generic PSD

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


      !-----------------------------------------------
      ! Select appropriate ice parametrization (see c_lspsiz)
      !-----------------------------------------------

  IF (lhook) CALL dr_hook('LSP_MELTING',zhook_in,zhook_handle)

  cry_offset=ice_type*ice_type_offset

      ! Use the generic ice particle size distribution
      ! Calculate the 1st (cx(84)) moment and the 1+0.5(di+1) (cx(85))
      ! moment of the ice particle size distribution
  IF (l_psd) THEN

    CALL lsp_moments(points,rho,t,qcf,cficei,cx(84),m_1)
    CALL lsp_moments(points,rho,t,qcf,cficei,cx(85),m_0p5_dp3)

  END IF

  DO i = 1, points

    IF (qcf(i) >  m0 .AND. t(i) >  zerodegc) THEN

          !-----------------------------------------------
          ! Calculate the average wet-bulb temperature (tempw)
          !-----------------------------------------------
          ! First calculate the average supersaturation
      IF (ice_type  ==  3) THEN
            ! Graupel does not use ice cloud fraction
        tempw(i) = MAX(qsl(i)-q(i),0.0)
      ELSE
            ! Use the full ice cloud fraction information
        tempw(i)=area_ice(i)*MAX(qsl(i)-q_ice(i),0.0)*cficei(i)
      END IF

          ! Now calculate the average amount (temp7) by which the
          ! wet-bulb temperature exceeds 0 degrees C.
      temp7(i) = t(i) - zerodegc - tempw(i)                             &
               * (tw1 + tw2 * (p(i)-tw3) - tw4*(t(i)-tw5) )
      temp7(i) = MAX(temp7(i),0.0)

          !-----------------------------------------------
          ! Calculate melting rate
          !-----------------------------------------------
      IF (l_psd) THEN
            ! Use generic particle size distribution
            ! constp(84) = ventilation coefficient 1
            ! constp(85) = ventilation coefficient 2
            ! constp(90) = 2 pi axial_ratio_correction
            !              * air_conductivity / Lf
        dpr(i) = rhor(i) * constp(90) * timestep *                      &
                 ( constp(84)*m_1(i)*corr2(i) +                         &
                   constp(85)*rocor(i)*m_0p5_dp3(i) )
      ELSE
            ! Use particle size distribution based on intercepts
        IF (ice_type  ==  3) THEN
              ! Graupel does not use ice cloud fraction
          pr02(i) = rho(i)*qcf(i)*constp(5+cry_offset)*tcgi(i)
        ELSE
              ! Use the full ice cloud information
          pr02(i) = rho(i)*qcf(i)*cficei(i) *                           &
                    constp(5+cry_offset)*tcgi(i)
        END IF  ! ice_type==3

        dpr(i)  = tcg(i) * constp(14+cry_offset) * timestep *           &
                  (constp(7+cry_offset) * corr2(i) *                    &
                  pr02(i)**cx(4+cry_offset)                             &
                + constp(8+cry_offset) * rocor(i) *                     &
                  pr02(i)**cx(5+cry_offset)) * rhor(i)

      END IF  ! l_psd

          !-----------------------------------------------
          ! Calculate transfer
          !-----------------------------------------------
          ! Solve implicitly in terms of temperature
      dpr(i) = temp7(i) * (1.0-1.0/(1.0+dpr(i)*lfrcp))/lfrcp

          ! Limit mass transfer to the mass available
      dpr(i) = MIN(dpr(i),qcf(i))

      ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

          !------------------------------------------------
          ! Adjust ice and rain contents
          !------------------------------------------------
      qcf(i)   = qcf(i)   - dpr(i)
      qrain(i) = qrain(i) + dpr(i)
      t(i)     = t(i)     - dpr(i) * lfrcp

          !------------------------------------------------
          ! Update rain fractions
          !------------------------------------------------
      IF (dpr(i)  >   0.0) THEN
        rf_final(i) = MAX(MAX(rainfrac(i),cff(i)),frac_ice_above(i))

        rftransfer(i) = rftransfer(i) + (rf_final(i) - rainfrac(i))

        rainfrac(i)= rf_final(i)
        rain_liq(i)= MIN(area_liq(i),rainfrac(i))
        rain_mix(i)= MIN(area_mix(i),rainfrac(i)-rain_liq(i))
        rain_ice(i)=                                                    &
             MIN(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
        rain_clear(i)=                                                  &
               rainfrac(i) - rain_liq(i)-rain_mix(i)-rain_ice(i)
        rain_clear(i) = MIN(area_clear(i),rain_clear(i))

      END IF  ! dpr  >   0

          !------------------------------------------------
          ! Update cloud fractions
          !------------------------------------------------
      IF (ice_type /= 3) THEN
            ! ( ice_type /=3 as graupel does not change
            !   cloud fractions )


            ! Assume that ice cloud fraction is reduced in
            ! proportion to the mass removed
        delta_cff(i) = -cff(i) * dpr(i) / qcft(i)

            ! Assume a random overlap between the removed ice
            ! cloud and the liquid cloud
        delta_cf(i)  = delta_cff(i) * area_ice(i) * cficei(i)

            ! Calculate transfer rate diagnostics
        cftransfer(i)  = cftransfer(i)                                  &
                         + delta_cf(i)  * one_over_tsi
        cfftransfer(i) = cfftransfer(i)                                 &
                         + delta_cff(i) * one_over_tsi

        cf(i)  = cf (i) + delta_cf (i)
        cff(i) = cff(i) + delta_cff(i)

      END IF  ! ice_type /= 3

    END IF  ! qcf  >   m0 etc.

  END DO  ! Points

  IF (lhook) CALL dr_hook('LSP_MELTING',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_melting
END MODULE lsp_melting_mod
