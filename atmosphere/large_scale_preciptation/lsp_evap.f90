! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Evaporation of rain
! Subroutine Interface:
MODULE lsp_evap_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_evap(                                                    &
  points, timestep,                                                     &
                                          ! Number of points and tstep
  p, q, qrain, t, q_ice, q_clear,                                       &
                                          ! Water contents and temp
  area_liq, area_mix,                                                   &
                                          ! Cloud fraction partitions
  area_ice, area_clear,                                                 &
  rainfrac, rain_liq, rain_mix,                                         &
                                          ! Rain fractions for updating
  rain_ice, rain_clear,                                                 &
  rho, corr, corr2,rocor,dhilsiterr,                                    &
                                          ! Parametrization information
  lcrcp, lheat_correc_liq, qsl, esw,                                    &
  ptransfer, rftransfer,                                                &
                                          ! Mass transfer diagnostic
  one_over_tsi                                                          &
                                          ! 1/(timestep*iterations)
  )

  ! Microphysics modules

  USE lsp_dif_mod,         ONLY: apb4, apb5, apb6
  USE mphys_ice_mod,       ONLY: qcfmin, m0
  USE mphys_constants_mod, ONLY: cx, constp, rho_q_veloc, lam_evap_enh, &
                                 max_as_enh

  USE mphys_inputs_mod,    ONLY: l_rainfall_as, l_warm_new, l_mcr_qrain


  ! Dr Hook Modules
  USE yomhook,             ONLY: lhook, dr_hook
  USE parkind1,            ONLY: jprb, jpim
 
  IMPLICIT NONE

! Purpose:
!   Update rain and water vapour due to evaporation

! Method:
!   Integrate evaporation rates of a single raindrop over the
!   raindrop size distribution
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

  INTEGER, INTENT(IN) ::                                                &
    points
                          ! Number of points to calculate

  REAL, INTENT(IN) ::                                                   &
    timestep,                                                           &
                          ! Timestep / s
    p(points),                                                          &
                          ! Air pressure / N m-2
    q_ice(points),                                                      &
                          ! Local vapour in ice partition / kg kg-1
    q_clear(points),                                                    &
                          ! Local vapour in clear partition / kg kg-1
    area_liq(points),                                                   &
                          ! Fraction of gridbox with liquid but no ice
    area_mix(points),                                                   &
                          ! Fraction of gridbox with liquid and ice
    area_ice(points),                                                   &
                          ! Fraction of gridbox with liquid and ice
    area_clear(points),                                                 &
                          ! Fraction of gridbox with no condensate
    lcrcp,                                                              &
                          ! Latent heat of condensation/cP / K
    dhilsiterr(points),                                                 &
                          ! layer thickness/timestep / m s-1
    corr(points),                                                       &
                          ! Fall speed correction factor due to
                          !   air density changes (no units)
    corr2(points),                                                      &
                          ! Air diffusivity correction (no units)
    rho(points),                                                        &
                          ! Air density / kg m-3
    lheat_correc_liq(points),                                           &
                                  ! Correction of subsaturation due
                                  ! to latent heat
    qsl(points),                                                        &
                          ! Saturated spec humidity wrt liquid / kg kg-1
    esw(points),                                                        &
                          ! Saturated vapour pressure wrt liquid / N m-2
    rocor(points),                                                      &
                          ! Air density and viscosity correction factor
                          !   (no units)
    one_over_tsi
                          ! 1/(timestep*iterations)

  REAL, INTENT(INOUT) ::                                                &
    qrain(points),                                                      &
                          ! Rain water content / kg kg-1
    q(points),                                                          &
                          ! Vapour content / kg kg-1
    t(points),                                                          &
                          ! Temperature / K
    ptransfer(points),                                                  &
                          ! Evaporation rate / kg kg-1 s-1
    rftransfer(points),                                                 &
                          ! Rate of change of rain fraction / s-1
    rainfrac(points),                                                   &
                          ! Rain fraction (no units)
    rain_liq(points),                                                   &
                          ! Overlap of rain with liquid (no units)
    rain_mix(points),                                                   &
                          ! Overlap of rain with mixed phase region
    rain_ice(points),                                                   &
                          ! Overlap of rain with ice
    rain_clear(points)! Overlap of rain with clear sky

!     Real, Intent(Out) ::

! Local Variables

  INTEGER ::                                                            &
    i                 ! Loop counter

  REAL ::                                                               &
    dpr(points),                                                        &
                          ! Amount of mass evaporated / kg kg-1
    pr04(points),                                                       &
                          ! Temporary in evaporation rate calc.
    lamr1(points),                                                      &
                          ! Powers of drop size distribution slope
    lamr2(points),                                                      &
                          !   lambda.
    temp7(points),                                                      &
                          ! Subsaturation in gridbox / kg kg-1
    rainfracnew(points),                                                &
                          ! Updated rain fraction

! Local variables for Abel & Shipway-style diagnostic rain evaporation
! Enhancement

    prelamnew(points),                                                  &
                          ! Precursor to new lambda calculated for A/S
    lamnew(points),                                                     &
                          ! Lambda for Abel and Shipway calculations
    standard_as_veloc(points),                                          &
                          ! Standard velocity for A/S
                          ! determined from rho_q_veloc input
    true_as_veloc(points),                                              &
                          ! Actual Abel and Shipway velocity
                          ! determined from lamnew
    evap_enh_factor(points)
                          ! Factor with which to increase (enhance)
                          ! the evaporation rate

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!End of header and all declarations

  IF (lhook) CALL dr_hook('LSP_EVAP',zhook_in,zhook_handle)
  DO i = 1, points
        !-----------------------------------------------
        ! If there is only a small amount of rain present, evaporate
        ! it completely to avoid numerical problems.
        !-----------------------------------------------
    IF (qrain(i)  <   qcfmin .AND. .NOT. l_warm_new .OR.                &
        qrain(i) < m0 .AND. l_warm_new) THEN
      ! original code allowed significant evaporation of rain in-cloud,
      ! which is unphysical. l_warm_new prevents this by lowering
      ! the threshold to define a "small amount"

          ! Evaporate all this rain
      dpr(i) = qrain(i)

        t(i)   = t(i) - lcrcp * dpr(i)
        q(i)   = q(i) + dpr(i)
        qrain(i) = 0.0

          ! Store evaporation rate
      ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

          ! Update rain fractions
      rftransfer(i)  = rftransfer(i)                                    &
                     -rainfrac(i) * one_over_tsi

        rainfrac(i)  = 0.0
        rain_liq(i)  = 0.0
        rain_mix(i)  = 0.0
        rain_ice(i)  = 0.0
        rain_clear(i)= 0.0

    END IF  ! qrain lt qcfmin

    IF (qrain(i)  >   0.0 .AND. rainfrac(i)  >   0.0) THEN
          !-----------------------------------------------
          ! Calculate evaporation parameters
          !-----------------------------------------------
      pr04(i) = ((apb4-apb5*t(i))*esw(i)+apb6*p(i)*t(i)**3)

          !-----------------------------------------------
          ! Calculate powers of drop size distribution slopes
          !-----------------------------------------------
      IF (l_mcr_qrain) THEN
            ! Rain is a mixing ratio - use fall speeds as parametrized
        lamr1(i) = qrain(i) * constp(50) * rho(i) / (rainfrac(i))
        lamr2(i) = lamr1(i) ** (cx(47)*cx(52))
        lamr1(i) = lamr1(i) ** (cx(49)*cx(52))
      ELSE
            ! Rain is a diagnostic quantity
            ! - use fall speeds to allow sedimentation into next layer
        lamr1(i) = qrain(i) * rho(i) * dhilsiterr(i) /                  &
                   (constp(42) * corr(i) * rainfrac(i))
        lamr2(i) = lamr1(i) ** (cx(47)*cx(48))
        lamr1(i) = lamr1(i) ** (cx(49)*cx(48))
      END IF  ! l_mcr_qrain

          !-----------------------------------------------
          ! Calculate evaporation rate
          !-----------------------------------------------
      dpr(i) = constp(46) * t(i)**2 * esw(i) * timestep
      dpr(i) = dpr(i) * ( (constp(47) * corr2(i) * lamr2(i))            &
                 + (constp(48) * rocor(i) * lamr1(i)) )


          !-----------------------------------------------
          ! Simulate  Abel & Shipway Section for diagnostic
          ! rain only
          !-----------------------------------------------

          ! Idea is that droplets are all small for big
          ! lambda and so for
          ! diagnostic rain they won't fall as far before
          ! evaporating.

      IF (.NOT. l_mcr_qrain .AND. l_rainfall_as) THEN

            ! Calculate new lambda for rain based on A/S (07)

        prelamnew(i) = qrain(i)*rho(i)*dhilsiterr(i)                    &
                   /(rainfrac(i)*constp(42)*corr(i))

        lamnew(i) =  1.0 / (prelamnew(i)**cx(48))

        IF (lamnew(i)  >  lam_evap_enh) THEN

              !Drops are smaller than mean of those at the
              !point at where the Abel and Shipway curve
              !diverges from the standard UM curve.

              !-----------------------------------------------------
              ! Calculate standard velocity
              !-----------------------------------------------------
          standard_as_veloc(i) = rho_q_veloc / (rho(i)*qrain(i))

              !-----------------------------------------------------
              ! Calculate actual Abel and Shipway Bulk Velocity
              !-----------------------------------------------------
          true_as_veloc(i) = (constp(57)/(rho(i)*qrain(i)))*            &
                             (lamnew(i)**cx(46)) *(                     &
         ( constp(54) / ((lamnew(i)+cx(56))**cx(59) ) ) +               &
         ( constp(55) / ((lamnew(i)+cx(57))**cx(60) ) )   )

              !-----------------------------------------------------
              ! Calculate evaporation enhancement factor
              ! and
              ! Limit enhancement factor to maximum value
              ! (Set in module)
              !-----------------------------------------------------
          evap_enh_factor(i) =                                          &
            MIN((standard_as_veloc(i)/true_as_veloc(i)), max_as_enh)

              !-----------------------------------------------------
              ! Finally, enhance evaporation rate
              !-----------------------------------------------------
          dpr(i) = dpr(i)*evap_enh_factor(i)

        END IF !lamnew(i) .GT. Lam_evap_enh

      END IF ! .not. l_mcr_qrain and l_rainfall_as

          !-----------------------------------------------
          ! Calculate transfers
          !-----------------------------------------------
          ! Calculate gridbox mean supersaturation
      temp7(i) = (q_ice(i) - qsl(i)) * rain_ice(i)                      &
             +(q_clear(i) - qsl(i)) * rain_clear(i)

          ! Limit on the gridbox mean supersaturation
      dpr(i) = dpr(i) * MAX(-temp7(i)*lheat_correc_liq(i),0.0)          &
             /(qsl(i) * rho(i) * pr04(i) + dpr(i))

          ! Limit on the amount of rain available
      IF (rain_mix(i) == 0.0 .AND. rain_liq(i) == 0.0) THEN
            ! No overlap between rain and liquid cloud, so
            ! (rain_ice + rain_clear) / rainfrac must be 1.
        dpr(i) = MIN( dpr(i), qrain(i))
      ELSE
        dpr(i) = MIN( dpr(i), qrain(i) *                                &
                (rain_ice(i) + rain_clear(i)) / rainfrac(i))
      END IF

          !-----------------------------------------------
          ! Store process rate (kg/kg/s)
          !-----------------------------------------------
      ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

          !-----------------------------------------------
          ! Update values of rain, vapour and temperature
          !-----------------------------------------------
        qrain(i) = qrain(i) - dpr(i)
        q(i)     = q(i)     + dpr(i)
        t(i)     = t(i)     - dpr(i) * lcrcp

          !-----------------------------------------------
          ! Calculate new rain fraction
          !-----------------------------------------------
!         These are commented out to ensure that rain fraction converges
!         to a 0 or 1 behaviour when RHcrit tends to 1.
!          rainfracnew(i) = rainfrac(i)*qrain(i)/(qrain(i)+dpr(i))
!          rftransfer(i)  = rftransfer(i)
!     &       + (rainfracnew(i) - rainfrac(i)) / (timestep*iterations)

          !-----------------------------------------------
          ! Update rain fractions
          !-----------------------------------------------
!         These are commented out to ensure that rain fraction converges
!         to a 0 or 1 behaviour when RHcrit tends to 1.
!       
!            rainfrac(i)  = rainfracnew(i)
!            rain_liq(i)  = min(area_liq(i) , rainfrac(i))
!            rain_mix(i)  = min(area_mix(i) , rainfrac(i)-rain_liq(i))
!            rain_ice(i)  =
!     &          min(area_ice(i) , rainfrac(i)-rain_liq(i)-rain_mix(i))
!            rain_clear(i)= rainfrac(i)-rain_liq(i)
!     &                    -rain_mix(i)-rain_ice(i)
!         

    END IF  !  qrain gt 0 etc.

  END DO  ! Points

  IF (lhook) CALL dr_hook('LSP_EVAP',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_evap
END MODULE lsp_evap_mod
