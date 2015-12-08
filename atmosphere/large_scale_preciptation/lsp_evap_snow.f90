! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Evaporation of melting snow
! Subroutine Interface:
MODULE lsp_evap_snow_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_evap_snow(                                               &
  points, timestep,                                                     &
                                          ! Number of points and tstep
  q, q_ice, qcf, qcft, t, p,                                            &
                                          ! Water contents, temp, pres
  esw, qsl,                                                             &
                                          ! Saturated quantities
  area_ice, cficei,                                                     &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  cf, cff,                                                              &
                                          ! Current cloud fractions for
                                          ! updating
  rho, tcg, tcgi,                                                       &
                                          ! Parametrization information
  corr2, rocor, lheat_correc_liq,                                       &
  lsrcp, ice_type,                                                      &
                                          ! Microphysical information
  l_psd,                                                                &
                                          ! Code options
  ptransfer,                                                            &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                         &
                                          ! 1/(timestep*iterations)
  cftransfer,cfftransfer                                                &
                                          ! Cloud transfer diagnostics
  )

  ! Microphysics Modules
  USE lsp_dif_mod,         ONLY: apb4, apb5, apb6
  USE mphys_ice_mod,       ONLY: m0
  USE mphys_constants_mod, ONLY: cx, constp, ice_type_offset

  ! General Atmosphere Modules
  USE conversions_mod,     ONLY: zerodegc

  ! Dr Hook Modules
  USE yomhook,             ONLY: lhook, dr_hook
  USE parkind1,            ONLY: jprb, jpim

  ! Large scale precipitation modules
  USE lsp_moments_mod,     ONLY: lsp_moments
  
  IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of sublimation of melting
!   snow

! Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles in a distribution of vapour
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! Source for vapour. Sink for ice.

! Subroutine Arguments

  INTEGER, INTENT(IN) ::                                                &
    points,                                                             &
                          ! Number of points to calculate
    ice_type          ! Type of ice (0 - crystals, 1 - aggregates)

  REAL, INTENT(IN) ::                                                   &
    timestep,                                                           &
                          ! Timestep / s
    p(points),                                                          &
                          ! Air pressure / N m-2
    q_ice(points),                                                      &
                          ! Vapour content in ice partition / kg kg-1
    esw(points),                                                        &
                          ! Saturated vapour pres. over liquid / N m-2
    qsl(points),                                                        &
                          ! Saturated humidity wrt liquid / kg kg-1
    area_ice(points),                                                   &
                          ! Fraction of gridbox with ice but no liquid
    cficei(points),                                                     &
                          ! 1/Fraction of gridbox with ice cloud
    rho(points),                                                        &
                          ! Air density / kg m-3
    tcg(points),                                                        &
                          ! T dependent function in ice size dist'n
    tcgi(points),                                                       &
                          ! 1/tcg (no units)
    corr2(points),                                                      &
                          ! Temperature correction factor (no units)
    rocor(points),                                                      &
                          ! Combined fall and corr2 correction factor
    lheat_correc_liq(points),                                           &
                                 ! Liquid latent heat correction factor
    lsrcp,                                                              &
                          ! Latent heat of sublimation
                          ! /heat capacity of air (cP) / K
    one_over_tsi          ! 1/(timestep*iterations)

  REAL, INTENT(INOUT) ::                                                &
    q(points),                                                          &
                          ! Vapour content / kg kg-1
    qcf(points),                                                        &
                          ! Ice water content in ice category to be
!                           updated    / kg kg-1
    qcft(points),                                                       &
                          ! Ice water in all ice categories
!                           (for cloud fraction calculations)
    t(points),                                                          &
                          ! Temperature / K
    cf(points),                                                         &
                          ! Current cloud fraction
    cff(points),                                                        &
                          ! Current ice cloud fraction
    ptransfer(points)  ! Mass deposited in this timestep / kg kg-1

  REAL, INTENT(INOUT) ::                                                &
    cftransfer(points),                                                 &
                           ! Cloud fraction increment this tstep
    cfftransfer(points)! Ice cloud fraction inc this tstep

  LOGICAL, INTENT(IN) ::                                                &
    l_psd
                          ! Use generic ice particle size distribution

! Local Variables

  INTEGER ::                                                            &
    i,                                                                  &
                          ! Loop counter for points
    cry_offset        ! Index offset for ice crystals

  REAL ::                                                               &
    pr02(points),                                                       &
                          ! Temporary in calculation of PSD slopes
    pr04(points),                                                       &
                          ! Temporary in calculation of PSD slopes
    dpr(points),                                                        &
                          ! Temporary in calculating ice transfers
    deltacf(points),                                                    &
                          ! Change in cf across timestep
    tempw(points),                                                      &
                          ! Available subsaturation for evaporation
    m_1(points),                                                        &
                          ! 1st moment of the generic ice size distribtn
    m_0p5_dp3(points)
                          ! 1+(di+1)/2 moment of the generic ice PSD

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


      !-----------------------------------------------
      ! Select appropriate ice parametrization (see c_lspsiz)
      !-----------------------------------------------

  IF (lhook) CALL dr_hook('LSP_EVAP_SNOW',zhook_in,zhook_handle)

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
          ! Diffusional growth parameters
          !-----------------------------------------------
      pr04(i) = ((apb4-apb5*t(i))*esw(i)+apb6*p(i)*t(i)**3)

          !-----------------------------------------------
          ! Calculate transfer rates
          !-----------------------------------------------
      IF (l_psd) THEN
            ! Use generic particle size distribution
            ! constp(83) = 2 pi axial_ratio_correction
            ! constp(84) = ventilation coefficient 1
            ! constp(85) = ventilation coefficient 2
            !        * Sc^(1/3)*ci^0.5/viscosity0^0.5
        dpr(i) = constp(83) * timestep * t(i)**2 * esw(i)               &
                            * (constp(84)*m_1(i)*corr2(i)               &
                            +  constp(85)*rocor(i)*m_0p5_dp3(i))        &
                            / (qsl(i) * rho(i) * pr04(i))

      ELSE
            ! Use particle size distribution based on intercepts
        pr02(i)=rho(i)*qcf(i)*cficei(i)*constp(5+cry_offset)*tcgi(i)
        dpr(i)=tcg(i)*constp(6+cry_offset)*t(i)**2*esw(i)*timestep*     &
        (constp(7+cry_offset)*corr2(i)*pr02(i)**cx(4+cry_offset)        &
         +constp(8+cry_offset)*rocor(i)*pr02(i)**cx(5+cry_offset))      &
         /(qsl(i)*rho(i)*pr04(i))

      END IF  ! l_psd

          !-----------------------------------------------
          ! Limit transfers
          !-----------------------------------------------
          ! tempw is the subsaturation that is available
      tempw(i) = area_ice(i) * (qsl(i) - q_ice(i))
      dpr(i) = dpr(i) * tempw(i)
      dpr(i) = MAX(MIN(dpr(i),tempw(i)*lheat_correc_liq(i)),0.0)

          ! Limit on the amount of ice available
      dpr(i) = MIN(dpr(i),qcf(i))

          !-----------------------------------------------
          ! Store process rate
          !-----------------------------------------------
      ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

          !-----------------------------------------------
          ! Store and update ice cloud fraction and total fractions
          !-----------------------------------------------

      deltacf(i)     = ( -cff(i) * dpr(i)/qcft(i) )
      cftransfer(i)  = cftransfer(i) + deltacf(i) * one_over_tsi
      cfftransfer(i) = cfftransfer(i)+ deltacf(i) * one_over_tsi

      cf(i)  = cf(i)  + deltacf(i)
      cff(i) = cff(i) + deltacf(i)


          !-----------------------------------------------
          ! Update values of ice and vapour
          !-----------------------------------------------

      qcf(i) = qcf(i) - dpr(i)
      q(i)   = q(i)   + dpr(i)
      t(i)   = t(i)   - dpr(i)*lsrcp

    END IF  ! qcf gt m0 etc.

  END DO  ! Points

  IF (lhook) CALL dr_hook('LSP_EVAP_SNOW',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_evap_snow
END MODULE lsp_evap_snow_mod
