! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Riming of ice particles
! Subroutine Interface:
MODULE lsp_tidy_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_tidy(                                                    &
  points, iterations,                                                   &
                                          ! Number of points and tstep
  one_over_tsi,                                                         &
                                          ! 1/(timestep*iterations)
  q, qcl, qcf, qcf2, qrain, t,                                          &
                                          ! Water contents and temp
  area_liq, area_mix, area_ice,                                         &
                                          ! Cloud fraction information
  cfice, cficei,                                                        &
                                          ! at start of microphysics ts
  cf, cfl, cff,                                                         &
                                          ! Current cloud fractions for
                                          ! updating
  rainfrac, rain_liq, rain_mix,                                         &
                                          ! Rain fractions
  rain_ice, rain_clear,                                                 &
  q_ice, qs, qsl, snow_agg, snow_cry,                                   &
                                          ! Other water contents
  rho, rhor, p,                                                         &
                                          ! Other model prognostics
  cttemp,dhi,dhilsiterr,frac_ice_above,                                 &
                                             ! Other information
  lcrcp, lfrcp, lsrcp,                                                  &
                                          ! Microphysical information
  psdep, pidep, psmlt, pimlt, prevp,                                    &
                                          ! Mass transfer diagnostics
  cftransfer,cfltransfer,cfftransfer,                                   &
                                          ! Cloud transfer diagnostics
  rftransfer                                                            &
                                          ! Rain transfer diagnostics
  )

  ! Microphysics modules
  USE mphys_ice_mod,    ONLY: qcfmin
  USE lsp_dif_mod,      ONLY: tw1, tw2, tw3, tw4, tw5
  USE mphys_inputs_mod, ONLY: l_it_melting, l_mcr_qcf2

  ! General atmosphere modules
  USE conversions_mod, ONLY: zerodegc

  ! Dr Hook Modules
  USE yomhook,         ONLY: lhook, dr_hook
  USE parkind1,        ONLY: jprb, jpim

  IMPLICIT NONE

! Purpose:
!   Tidy up small numerical values after large-scale precipitation.
!   Ideally, this code would be unnecessary, but in reality there
!   are always going to be small values left over numerically that
!   should be reset.

! Method:
!   1. Evaporate rain amounts
!   2. Evaporate small ice amounts
!   3. Reset cloud fractions
!   4. Melt small snow amounts
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
    iterations        ! Number of iterations

  REAL, INTENT(IN) ::                                                   &
    one_over_tsi,                                                       &
                          ! 1/(timestep*iterations)
    qcl(points),                                                        &
                          ! Liquid water content / kg kg-1
    area_liq(points),                                                   &
                          ! Fraction of gridbox with liquid-only cloud
    area_mix(points),                                                   &
                          ! Fraction of gridbox with mixed phase cloud
    area_ice(points),                                                   &
                          ! Fraction of gridbox with ice-only cloud
    cfice(points),                                                      &
                          ! Fraction of gridbox with ice cloud
    cficei(points),                                                     &
                          ! 1/Fraction of gridbox with ice cloud
    q_ice(points),                                                      &
                          ! Vapour content in ice cloud / kg kg-1
    qs(points),                                                         &
                          ! Saturated water content wrt ice / kg kg-1
    qsl(points),                                                        &
                          ! Saturated water content wrt liquid / kg kg-1
    rho(points),                                                        &
                          ! Air density / kg m-3
    rhor(points),                                                       &
                          ! 1 / air density / m3 kg-1
    p(points),                                                          &
                          ! Air pressure / N m-2
    dhi(points),                                                        &
                          ! Timestep/layer thickness / s m-1
    dhilsiterr(points),                                                 &
                          ! 1/(dhi*iterations) / m s-1
    frac_ice_above(points),                                             &
                               ! Ice cloud fraction in layer above
    lcrcp,                                                              &
                          ! Latent heat of condensation/cP / K
    lfrcp,                                                              &
                          ! Latent heat of fusion/cP / K
    lsrcp             ! Latent heat of sublimation/cP / K

  REAL, INTENT(INOUT) ::                                                &
    q(points),                                                          &
                          ! Vapour mixing ratio / kg kg-1
    qcf(points),                                                        &
                          ! Aggregates mixing ratio / kg kg-1
    qcf2(points),                                                       &
                          ! Crystals mixing ratio / kg kg-1
    qrain(points),                                                      &
                          ! Rain mixing ratio / kg kg-1
    t(points),                                                          &
                          ! Temperature / K
    cttemp(points),                                                     &
                          ! Ice-cloud top temperature / K
    cf(points),                                                         &
                          ! Current cloud fraction
    cfl(points),                                                        &
                          ! Current liquid cloud fraction
    cff(points),                                                        &
                          ! Current ice cloud fraction
    rainfrac(points),                                                   &
                          ! Rain fraction (no units)
    rain_liq(points),                                                   &
                          ! Overlap of rain with liquid (no units)
    rain_mix(points),                                                   &
                          ! Overlap of rain with mixed phase region
    rain_ice(points),                                                   &
                          ! Overlap of rain with ice
    rain_clear(points),                                                 &
                          ! Overlap of rain with clear sky
    snow_agg(points),                                                   &
                          ! Aggregate snowfall rate / kg m-2 s-1
    snow_cry(points),                                                   &
                          ! Crystal snowfall rate / kg m-2 s-1
    psdep(points),                                                      &
                          ! Deposition of aggregates diag. / kg kg s-1
    pidep(points),                                                      &
                          ! Deposition of crystals diag. / kg kg s-1
    psmlt(points),                                                      &
                          ! Melting of aggregates diagnostic / kg kg s-1
    pimlt(points),                                                      &
                          ! Melting of crystals diagnostic / kg kg s-1
    prevp(points),                                                      &
                          ! Evaporation of rain diagnostic / kg kg s-1
    cftransfer(points),                                                 &
                           ! Rate of change of bulk cloud frac / s-1
    cfltransfer(points),                                                &
                           ! Rate of change of liquid cloud frac / s-1
    cfftransfer(points),                                                &
                           ! Rate of change of ice cloud frac / s-1
    rftransfer(points)! Rate of change of rain fraction / s-1


! Local Variables


  REAL, PARAMETER :: smallnum = 2.2e-14
                          ! Small value used in if tests

  INTEGER ::                                                            &
    i                 ! Loop counter for points

  REAL ::                                                               &
    dpr(points),                                                        &
                          ! Mass transfer / kg kg-1
    dpr2(points),                                                       &
                          ! Equivalent to dpr, but with units kg m-2 s-1
    temp7(points),                                                      &
                          ! Wet bulb temperature / deg C
    tempw(points),                                                      &
                          ! Saturation excess / kg kg-1
    cfnew(points),                                                      &
                          ! New cloud fraction
    cflnew(points),                                                     &
                          ! New liquid cloud fraction
    cffnew(points),                                                     &
                          ! New ice cloud fraction
    t_melt            ! Emergency melting temperature

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('LSP_TIDY',zhook_in,zhook_handle)
  DO i = 1, points

!    IF ((qrain(i)  <=  qcfmin .AND. qcl(i) <= 0.0)                      &
    IF ((qrain(i)  <=  qcfmin .AND. qcl(i) < 1.0E-16)                   &
         .OR. qrain(i)  <   0.0) THEN
          !-----------------------------------------------
          ! 1. If there is a tiny rain amount and no liquid then
          ! evaporate the rain
          !-----------------------------------------------
      dpr(i) = qrain(i)

          ! Update prognostics
      q(i)   = q(i) + dpr(i)
      t(i)   = t(i) - dpr(i) * lcrcp
      qrain(i) = 0.0

          ! Update rain fractions
      rainfrac(i) = 0.0
      rain_liq(i) = 0.0
      rain_mix(i) = 0.0
      rain_ice(i) = 0.0
      rain_clear(i) = 0.0

          ! Update evaporation rate
      rftransfer(i) = rftransfer(i)                                     &
                    - rainfrac(i) * one_over_tsi
      prevp(i) = prevp(i) + dpr(i) * one_over_tsi

    END IF  ! rain lt qcfmin etc.

        !-----------------------------------------------
        ! 2. Evaporate small ice amounts
        !-----------------------------------------------
    IF (l_mcr_qcf2) THEN
          ! Evaporate small aggregate and crystal amounts separately

          !-----------------------------------------------
          ! 2a. Aggregates
          !-----------------------------------------------
      IF ((qcf(i)  <   qcfmin) .AND.                                    &
        (t(i)  >   zerodegc .OR.                                        &
        (q_ice(i)  <=  qs(i) .AND. area_mix(i)  <  smallnum)            &
        .OR. qcf(i) <  0.0)) THEN
            ! Ice is very small and T>0 and ice is not growing by
            ! deposition, so evaporate it.
        dpr(i) = qcf(i)

            ! Update prognostics

        q(i)   = q(i) + dpr(i)
        t(i)   = t(i) - lsrcp * dpr(i)
        qcf(i) = 0.0


            ! Update deposition rate
        psdep(i) = psdep(i) - dpr(i) * one_over_tsi

      END IF  ! qcf < qcfmin etc

          !-----------------------------------------------
          ! 2b. Crystals
          !-----------------------------------------------
      IF ((qcf2(i)  <   qcfmin) .AND.                                   &
        (t(i)  >   zerodegc .OR.                                        &
        (q_ice(i)  <=  qs(i) .AND. area_mix(i)  <  smallnum)            &
        .OR. qcf2(i)  <   0.0)) THEN
            ! Ice is very small and T>0 and ice is not growing by
            ! deposition, so evaporate it.
        dpr(i) = qcf2(i)

            ! Update prognostics
        q(i)   = q(i) + dpr(i)
        t(i)   = t(i) - lsrcp * dpr(i)
        qcf2(i)= 0.0

            ! Update deposition rate
        pidep(i) = pidep(i) - dpr(i) * one_over_tsi

      END IF  ! qcf < 0 etc.

      IF (qcf(i)  <=  0.0 .AND. qcf2(i)  <=  0.0) THEN

            !----------------------------------------------
            ! 2c. Update ice cloud amounts
            !----------------------------------------------

        cfftransfer(i) = -cff(i) * one_over_tsi
        cftransfer(i)  = ( cfl(i) - cff(i) ) * one_over_tsi
        cff(i)         = 0.0
        cf(i)          = cfl(i)


            !----------------------------------------------
            ! 2d. Update cloud top temperature
            !----------------------------------------------
        cttemp(i) = t(i)
      END IF  ! qcf eq 0 etc

    ELSE  ! l_mcr_qcf2

          !----------------------------------------------
          ! 2e. One prognostic ice content is active
          !----------------------------------------------
      IF ((qcf(i)  <   qcfmin) .AND.                                    &
        (t(i)  >   zerodegc .OR.                                        &
        (q_ice(i)  <=  qs(i) .AND. area_mix(i)  <  smallnum)            &
        .OR. qcf(i) <= 0.0) )  THEN
            ! Ice is very small and T>0 and ice is not growing by
            ! deposition, so evaporate it.
        dpr(i) = qcf(i)

            ! Update prognostics

        q(i)   = q(i) + dpr(i)
        t(i)   = t(i) - lsrcp * dpr(i)
        qcf(i) = 0.0

            ! Update deposition rate
        psdep(i) = psdep(i) - dpr(i) * one_over_tsi

            ! Update ice cloud amounts
        cfftransfer(i) = -cff(i) * one_over_tsi
        cftransfer(i)  = (cfl(i)-cff(i)) * one_over_tsi

        cff(i) = 0.0
        cf(i)  = cfl(i)

            ! Update cloud top temperature
        cttemp(i) = t(i)

      END IF  ! qcf lt qcfmin etc

    END IF  ! l_mcr_qcf2

        !------------------------------------------------
        ! 3. Limit cloud fractions to physically reasonable values
        !------------------------------------------------
        ! It isn't clear that this ought to be a parallel calculation
        ! but it is allowed to be for the present time.

          ! Calculate new cloud fraction values
      cffnew(i) = MAX(MIN(cff(i),1.0),0.0)
      cflnew(i) = MAX(MIN(cfl(i),1.0),0.0)
      cfnew (i) = MAX( MAX(cflnew(i),cffnew(i)) , cf(i) )
      cfnew (i) = MIN( MIN(cflnew(i)+cffnew(i),1.0), cfnew(i) )

          ! Calculate transfer rates
      cfftransfer(i) = cfftransfer(i)                                   &
                     + (cffnew(i) - cff(i)) * one_over_tsi
      cfltransfer(i) = cfltransfer(i)                                   &
                     + (cflnew(i) - cfl(i)) * one_over_tsi
      cftransfer(i)  = cftransfer(i)                                    &
                     + (cfnew(i)  - cf (i)) * one_over_tsi

          ! Update cloud fractions

      cff(i) = cffnew(i)
      cfl(i) = cflnew(i)
      cf(i)  = cfnew(i)


        !------------------------------------------------
        ! 4. Emergency melting of snow to avoid excess snow at
        !    warm temperatures
        !------------------------------------------------
    IF (l_it_melting) THEN
          ! Use a warmer temperature if iterative melting is active
      t_melt = zerodegc + 2.0
    ELSE
      t_melt = zerodegc
    END IF  ! l_it_melting

        !------------------------------------------------
        ! Melt snow_agg first
        !------------------------------------------------
    IF (snow_agg(i)  >   0.0 .AND. t(i)  >   t_melt) THEN

          ! Numerical approximation of wet bulb temperature excess
      tempw(i) = area_ice(i) * MAX(qsl(i)-q_ice(i),0.0)*cficei(i)
      temp7(i) = t(i) - zerodegc                                        &
                 -tempw(i)*(tw1+tw2*(p(i)-tw3)-tw4*(t(i)-tw5))
      temp7(i) = MAX(temp7(i),0.0)

          ! Calculate transfer rate
      dpr(i)  = temp7(i) / lfrcp ! Rate based on Tw excess
      dpr2(i) = dpr(i)*rho(i)*dhilsiterr(i)

          ! Limit to the amount of snow available
      dpr(i)  = MIN(dpr(i) , snow_agg(i)                                &
                           * dhi(i)*iterations*rhor(i) )
      dpr2(i) = MIN(dpr2(i), snow_agg(i))

          ! Add to melting rate
      psmlt(i) = psmlt(i) + dpr(i) * one_over_tsi

          ! Update values of snow and rain

      snow_agg(i) = snow_agg(i) - dpr2(i)
      qrain(i)    = qrain(i)    + dpr(i)
      t(i)        = t(i)        - dpr(i) * lfrcp

    END IF  !  snow_agg gt 0 etc

        !------------------------------------------------
        ! Melt snow_cry next
        !------------------------------------------------
    IF (l_mcr_qcf2 .AND.                                                &
        snow_cry(i)  >   0.0 .AND. t(i)  >   t_melt) THEN

          ! Numerical approximation of wet bulb temperature excess
      tempw(i) = area_ice(i) * MAX(qsl(i)-q_ice(i),0.0)*cficei(i)
      temp7(i) = t(i) - zerodegc                                        &
                 - tempw(i)*(tw1+tw2*(p(i)-tw3)-tw4*(t(i)-tw5))
      temp7(i) = MAX(temp7(i),0.0)

          ! Calculate transfer rate
      dpr(i)  = temp7(i) / lfrcp ! Rate based on Tw excess
      dpr2(i) = dpr(i)*rho(i)*dhilsiterr(i)

          ! Limit to the amount of snow available
      dpr(i)  = MIN(dpr(i) , snow_cry(i)                                &
                           * dhi(i)*iterations*rhor(i) )
      dpr2(i) = MIN(dpr2(i), snow_cry(i))

          ! Add to melting rate
      pimlt(i) = pimlt(i) + dpr(i) * one_over_tsi

          ! Rain fraction will take on the value of the ice fraction
      rftransfer(i) = rftransfer(i) +                                   &
        MAX( MAX(frac_ice_above(i),cff(i))-rainfrac(i) , 0.0)           &
        * one_over_tsi

          ! Update values of snow and rain

      snow_cry(i) = snow_cry(i) - dpr2(i)
      qrain(i)    = qrain(i)    + dpr(i)
      t(i)        = t(i)        - dpr(i) * lfrcp

            ! Update rain fractions
      rainfrac(i) = MAX(rainfrac(i),cfice(i))
      rain_liq(i) = MIN(area_liq(i),rainfrac(i))
      rain_mix(i) = MIN(area_mix(i),rainfrac(i)-rain_liq(i))
      rain_ice(i) =                                                     &
             MIN(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
      rain_clear(i) =                                                   &
             rainfrac(i)-rain_liq(i)-rain_mix(i)-rain_ice(i)

    END IF  ! snow_cry lt 0 etc

  END DO  ! Points

  IF (lhook) CALL dr_hook('LSP_TIDY',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_tidy
END MODULE lsp_tidy_mod
