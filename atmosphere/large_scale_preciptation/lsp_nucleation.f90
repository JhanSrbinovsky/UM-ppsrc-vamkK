! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Nucleation of ice particles
MODULE lsp_nucleation_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_nucleation(                                              &
  points,                                                               &
                                          ! Number of points and tstep
  q, qcl, qrain, qcf, t,                                                &
                                          ! Water contents, temperature
  qs, qsl,                                                              &
                                          ! Saturated quantities
  cfliq, cfice,                                                         &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  area_liq, area_mix, cf, cfl, cff, rainfrac,                           &
                                          ! Current cloud and rain 
                                          ! fractions for updating
  rain_liq, rain_mix,                                                   &
                                          ! Overlaps of cloud with rain
  rhor, lheat_correc_ice,                                               &
                                          ! Parametrization information
  lfrcp, lsrcp,                                                         &
                                          ! Microphysical information
  hettransfer, homtransfer, homtransfer2,                               &
                                          ! Mass transfer diagnostics
  one_over_tsi,                                                         &
                                          ! 1/(timestep*iterations)
  cftransfer, cfltransfer, cfftransfer, rf_transfer_diag                &
                                          ! Cloud transfer diagnostics
  )

  ! Microphysics Modules
  USE mphys_ice_mod,    ONLY: thomo,  m0
  USE mphys_inputs_mod, ONLY: tnuc

  ! General atmosphere modules
  USE conversions_mod,  ONLY: zerodegc

  ! Dr Hook modules
  USE yomhook,         ONLY: lhook, dr_hook
  USE parkind1,        ONLY: jprb, jpim

  IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of homogeneous and
!   heterogeneous ice nucleation

! Method:
!   Homogeneous nucleation converts all liquid to ice with a temperature
!   threshold. Heterogeneous nucleation converts a seed amount of liquid
!   or vapour to ice.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! Subroutine Arguments

  INTEGER, INTENT(IN) ::                                                &
    points
                          ! Number of points to calculate
  REAL, INTENT(IN) ::                                                   &
    qs(points),                                                         &
                          ! Saturated humidity wrt ice / kg kg-1
    qsl(points),                                                        &
                          ! Saturated humidity wrt liquid / kg kg-1
    cfliq(points),                                                      &
                          ! Fraction of gridbox with liquid cloud
    cfice(points),                                                      &
                          ! Fraction of gridbox with ice cloud
    area_liq(points),                                                   &
                          ! Fraction of gridbox with only liquid cloud
    area_mix(points),                                                   &
                          ! Fraction of gridbox with mixed phase cloud
    rain_liq(points),                                                   &
                          ! Frac. of gridbox with rain and liquid cld.
    rain_mix(points),                                                   &
                          ! Frac. of grbx. with rain and mixed ph. cld
    rhor(points),                                                       &
                          ! 1/Air density / kg m-3
    lheat_correc_ice(points),                                           &
                                 ! Ice latent heat correction factor
!    &, m0                ! Seed ice water content / kg kg-1 c_lspmic
    lfrcp,                                                              &
                          ! Latent heat of fusion
                          ! /heat capacity of air (cP) / K
    lsrcp,                                                              &
                          ! Latent heat of sublimation/cP / K
    one_over_tsi          ! 1/(timestep*iterations)

  REAL, INTENT(INOUT) ::                                                &
    q(points),                                                          &
                          ! Vapour content / kg kg-1
    qcl(points),                                                        &
                          ! Liquid water content / kg kg-1
    qrain(points),                                                      &
                          ! Rain mixing ratio / kg kg-1
    qcf(points),                                                        &
                          ! Ice water content in ice category to be
!                           updated    / kg kg-1
    t(points),                                                          &
                          ! Temperature / K
    cf(points),                                                         &
                          ! Current cloud fraction
    cfl(points),                                                        &
                          ! Current liquid cloud fraction
    cff(points),                                                        &
                          ! Current ice cloud fraction
    rainfrac(points),                                                   &
                          ! Current rain fraction
    homtransfer(points),                                                &
                          ! Mass homog. nucleated this ts / kg kg-1
    homtransfer2(points),                                               &
                          ! Mass homog. (rain)
                          !  nucleated this ts / kg kg-1
    hettransfer(points)   
                          ! Mass het. nucleated this ts / kg kg-1

  REAL, INTENT(INOUT) ::                                                &
    cftransfer(points),                                                 &
                             ! Cloud fraction increment this tstep
    cfltransfer(points),                                                &
                             ! Liquid cloud fraction inc this tstep
    cfftransfer(points),                                                &
                             ! Ice cloud fraction inc this tstep
    rf_transfer_diag(points)
                             ! Rain fraction inc this tstep

! Local Variables

  INTEGER ::                                                            &
    i                 ! Loop counter for points

  REAL ::                                                               &
    dqi(points),                                                        &
                          ! Temporary in calculating ice transfers
    dqil(points),                                                       &
                          ! Temporary in calculating liquid transfers
    rhnuc             ! Nucleation relative humidity threshold

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('LSP_NUCLEATION',zhook_in,zhook_handle)
  DO i = 1, points

        !-----------------------------------------------
        ! Homogeneous nucleation - freeze all liquid if T < threshold
        !-----------------------------------------------
    IF (t(i) <  (zerodegc+thomo) .AND. qcl(i) >  0.0) THEN

          !-----------------------------------------------
          ! Update cloud fractions
          !-----------------------------------------------
      cfftransfer(i) = cfftransfer(i) + (cf(i) - cff(i) )* one_over_tsi
      cfltransfer(i) = cfltransfer(i) -  cfl(i)* one_over_tsi

      cff(i) = cf(i)
      cfl(i) = 0.0

          !-----------------------------------------------
          ! Update water contents
          !-----------------------------------------------
      homtransfer(i) = homtransfer(i) + qcl(i)*one_over_tsi

      qcf(i) = qcf(i) + qcl(i)
      t(i)   = t(i)   + lfrcp * qcl(i)
      qcl(i) = 0.0

    END IF ! T lt 0+thomo etc.

        !-----------------------------------------------
        ! Homogeneous nucleation - freeze all rain if T < threshold
        !-----------------------------------------------
    IF (t(i) <  (zerodegc+thomo) .AND. qrain(i) >  0.0) THEN

          !-----------------------------------------------
          ! Update rain fractions
          !-----------------------------------------------
      cfftransfer(i)      = cfftransfer(i) +                            &
                          ( cf(i) - cff(i) ) * one_over_tsi

      rf_transfer_diag(i) = rf_transfer_diag(i) -                       &
                            rainfrac(i) * one_over_tsi

      cff(i)      = cf(i)
      rainfrac(i) = 0.0

          !-----------------------------------------------
          ! Update water contents
          !-----------------------------------------------
      homtransfer2(i) = homtransfer2(i) + qrain(i)*one_over_tsi

      qcf(i)   = qcf(i) + qrain(i)
      t(i)     = t(i)   + lfrcp * qrain(i)
      qrain(i) = 0.0

    END IF ! T lt 0+thomo etc.

        !-----------------------------------------------
        ! Heterogeneous nucleation for T < tnuc
        !-----------------------------------------------
    IF (t(i) <  (zerodegc+tnuc) .AND. area_liq(i) >  0.0                &
        .AND. cfliq(i)  >   0.0 ) THEN

          !-----------------------------------------------
          ! Calculate number and mixing ratio of nucleated crystals
          !-----------------------------------------------
      dqi(i)=MIN(0.01*EXP(-0.6*(t(i)-zerodegc)),1.0e5)
      dqi(i)=m0 * dqi(i) * rhor(i)

          !-----------------------------------------------
          ! How much moisture is available for ice formation
          !-----------------------------------------------
          ! Firstly, calculate the threshold relative humidity
      rhnuc=(188.92+2.81*(t(i)-zerodegc)                                &
                  +0.013336*(t(i)-zerodegc)**2)*0.01
      rhnuc=MIN(rhnuc,1.0)-0.1
      rhnuc=MAX(qsl(i)*rhnuc,qs(i))

          ! Next calculate the available water
      dqil(i)=(qcl(i)/cfliq(i)+qsl(i)-rhnuc)
      dqi(i) = MAX(area_liq(i)*                                         &
                   MIN(dqi(i),dqil(i)*lheat_correc_ice(i)),0.0)
      qcf(i) = qcf(i)+dqi(i)

          !-----------------------------------------------
          ! Store nucleation rate
          !-----------------------------------------------
      hettransfer(i) = hettransfer(i)+dqi(i)*one_over_tsi

          !-----------------------------------------------
          ! Calculate mass transfers
          !-----------------------------------------------

            ! Firstly take mass from liquid water

        IF (area_mix(i) == 0.0) THEN ! area_liq/cfliq must be 1

          dqil(i) = MIN(dqi(i),qcl(i))

        ELSE

          dqil(i) = MIN(dqi(i),qcl(i)*area_liq(i)/cfliq(i))

        END IF

      qcl(i)  = qcl(i)-dqil(i)
      t(i)    = t(i)+lfrcp*dqil(i)
            ! If more ice is needed then the mass comes from vapour
      dqi(i)  = dqi(i)-dqil(i)
      t(i)    = t(i)+lsrcp*dqi(i)
      q(i)    = q(i)-dqi(i)

          !-----------------------------------------------
          ! Udate cloud fractions.
          !-----------------------------------------------

      cfftransfer(i) = cfftransfer(i) + (cf(i) - cff(i))*one_over_tsi
      cff(i)      = cf(i)


    END IF  ! On temperature threshold

  END DO  ! Points

  IF (lhook) CALL dr_hook('LSP_NUCLEATION',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_nucleation
END MODULE lsp_nucleation_mod
