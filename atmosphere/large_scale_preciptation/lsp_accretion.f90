! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Accretion of droplets by raindrops
! Subroutine Interface:
MODULE lsp_accretion_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_accretion(                                               &
  points, timestep,                                                     &
                                          ! Number of points and tstep
  qcl, qrain,                                                           &
                                          ! Water contents
!    &, area_liq, area_mix, area_ice      ! Partition information
  cfliq,                                                                &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  rainfrac, rain_liq, rain_mix,                                         &
                                          ! Rain fraction information
!    &, rain_ice, rain_clear
!    &, cf, cfl                           ! Current cloud fractions for
!                                         ! updating
  rho, corr, dhilsiterr,                                                &
                                          ! Parametrization information
  ptransfer,                                                            &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                         &
                                          ! 1/(timestep*iterations)
!    &, cftransfer,cfltransfer            ! Cloud transfer diagnostics
!    &, rftransfer                        ! Rain fraction transfer
  r_theta_levels_c, fv_cos_theta_latitude_c                             &
  )


  ! Microphysics Modules
  USE mphys_constants_mod, ONLY: cx, constp, acc_pref, acc_qc,          &
                                 acc_qr, l_inhomog
  USE mphys_inputs_mod,    ONLY: l_rainfall_as, l_warm_new, c_r_correl, &
                                 l_mcr_qrain
  USE mphys_ice_mod, ONLY: qcfmin
  USE mphys_bypass_mod, ONLY: mp_dell, mp_delp

  ! Dr Hook Modules
  USE yomhook,             ONLY: lhook, dr_hook
  USE parkind1,            ONLY: jprb, jpim

  IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of accretion of cloud
!   droplets by rain drops

! Method:
!   Solve implicitly the microphysical transfer equation for a
!   specified distribution of raindrops sweeping out a volume
!   of air uniformally inhabited by cloud water droplets.
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
    timestep,                                                           &
                          ! Timestep / s
    one_over_tsi,                                                       &
                          ! 1/(timestep*iterations)
    rain_liq(points),                                                   &
                          ! Overlap fraction of rain and liquid
    rain_mix(points),                                                   &
                          ! Overlap fraction of rain and mixed phase
    rainfrac(points),                                                   &
                          ! Rain fraction
    cfliq(points),                                                      &
                          ! Liquid cloud fraction at start of timestep
!    &, cf(points)        ! Current cloud fraction
!    &, cfl(points)       ! Current liquid cloud fraction
    rho(points),                                                        &
                          ! Air density / kg m-3
    corr(points),                                                       &
                          ! Fall velocity correction factor (no units)
    r_theta_levels_c(points),                                           &
                          ! Distance from centre of Earth and ... 
    fv_cos_theta_latitude_c(points),                                    &
                          ! ... grid info for working out gridbox size.
    dhilsiterr(points)! Depth of layer / timestep  / m s-1

  REAL, INTENT(INOUT) ::                                                &
    qcl(points),                                                        &
                          ! Liquid water content / kg kg-1
    qrain(points),                                                      &
                          ! Rain mixing ratio / kg kg-1
!    &, rain_ice(points)  ! Overlap fraction of rain and ice
!    &, rain_clear(points)! Overlap fraction of rain and clear sky
    ptransfer(points) ! Mass rimed in this timestep / kg kg-1
!    &, rftransfer(points)! Transfer of rain fraction this timestep
!    &, cftransfer(points) ! Cumulative cloud fraction increment
!    &, cfltransfer(points)! Cumulative liquid cloud fraction increment

! Local Variables

  INTEGER ::                                                            &
    i                 ! Loop counter for points

  REAL ::                                                               &
    dpr(points),                                                        &
                          ! Transfer of mixing ratio  / kg kg-1
    temp1(points),                                                      &
                          ! Temporary in rate calculations
    lambda(points),                                                     &
                          ! Temporary lambda for Abel/Shipway
    lambda_h1(points),                                                  &
                          ! Temporary lambda + h1r for Abel/Shipway
    lambda_h2(points),                                                  &
                          ! Temporary lambda + h2r for Abel/Shipway
    temp7(points),                                                      &
                          ! Rain free liquid cloud fraction (no units)
    qclnew(points)    ! Updated value of liquid water / kg kg-1

  REAL ::                                                               &
    fsd_qc,                                                             &
            ! fractional standard dev of cloud water content
    fsd_qr,                                                             &
            ! fractional standard dev of rain water content
    bias,                                                               &
            ! accretion rate bias
    x_in_km ! grid-box size in km

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('LSP_ACCRETION',zhook_in,zhook_handle)

  IF (l_warm_new) THEN

  !---------------------------------------------
  ! Use new accretion scheme
  !---------------------------------------------

    DO i = 1, points

      IF (rainfrac(i) >  0.0 .AND. qrain(i) >  0.0                      &
          .AND. cfliq(i) >  0.0) THEN

        IF (l_inhomog) THEN

          ! calculate bias E factor based on Boutle et al 2012 QJ
          x_in_km = 0.001*SQRT (   r_theta_levels_c(i) * mp_dell        &
                                 * r_theta_levels_c(i) * mp_delp        &
                                 * fv_cos_theta_latitude_c(i)     )

          IF ( cfliq(i) < 1.0 ) THEN
            fsd_qc=(0.45-0.25*cfliq(i))*(((x_in_km*cfliq(i))**0.333)    &
                   *((0.06*x_in_km*cfliq(i))**1.5+1.0)**(-0.17))
          ELSE
            fsd_qc = 0.11*(((x_in_km*cfliq(i))**0.333)                  &
                     *((0.06*x_in_km*cfliq(i))**1.5+1.0)**(-0.17))
          END IF

          fsd_qr = (1.1-0.8*rainfrac(i))                                &
                   *(((x_in_km*rainfrac(i))**0.333)                     &
                   *((0.11*x_in_km*rainfrac(i))**1.14+1.0)**(-0.22))

          fsd_qc = fsd_qc*1.414
          fsd_qr = fsd_qr*1.414

          bias = ((1+fsd_qc**2)**(-0.5*acc_qc))*                        &
                 ((1+fsd_qc**2)**(0.5*acc_qc**2))*                      &
                 ((1+fsd_qr**2)**(-0.5*acc_qr))*                        &
                 ((1+fsd_qr**2)**(0.5*acc_qr**2))*                      &
                 EXP(c_r_correl*acc_qc*acc_qr*                          &
                 SQRT(LOG(1+fsd_qc**2)*LOG(1+fsd_qr**2)))

        ELSE ! no inhomog param

           bias = 1.0

        END IF

        dpr(i) = acc_pref * bias *                                      &
                 ( (qcl(i)/cfliq(i))**acc_qc) *                         &
                 ( (qrain(i)/rainfrac(i))**acc_qr )                     &
                 * timestep * (rain_liq(i)+rain_mix(i))
        dpr(i) = MAX(MIN(dpr(i),                                        &
                 qcl(i)*(rain_liq(i)+rain_mix(i))/cfliq(i)-qcfmin),0.0)

        ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

        !------------------------------------------------
        ! Adjust liquid and rain contents
        !------------------------------------------------
        qrain(i) = qrain(i) + dpr(i)
        qcl(i)   = qcl(i)   - dpr(i)

      END IF  ! rainfrac > 0.0 etc.

    END DO  ! Points

  ELSE ! original accretion

    DO i = 1, points

      IF (rainfrac(i) >  0.0 .AND. qrain(i) >  0.0                      &
          .AND. cfliq(i) >  0.0) THEN

        IF (l_mcr_qrain) THEN
                ! Use the prognostic rain formulation
          temp1(i) = qrain(i)*rho(i)*constp(50)/(rainfrac(i))
        ELSE
                ! Use the diagnostic rain formulation
          temp1(i) = qrain(i)*rho(i)*dhilsiterr(i)                      &
                   /(rainfrac(i)*constp(42)*corr(i))
  
        END IF  ! l_mcr_qrain
  
            !-----------------------------------------------
            ! Calculate the new local value of qcl
            !-----------------------------------------------
        IF (l_rainfall_as) THEN
  
          IF (l_mcr_qrain) THEN
  
               !Prognostic Abel and Shipway version
  
               !First need to calculate lambda:
               ! lambda = (1 / (rain fraction temp1)) to power of cx(52)
  
            lambda(i) = ( 1.0 / (rainfrac(i)*temp1(i)) )**cx(52)
  
            lambda_h1(i) = lambda(i)+cx(56)
            lambda_h2(i) = lambda(i)+cx(57)
  
            qclnew(i)= qcl(i) / (                                       &
                       cfliq(i)+(cfliq(i)*timestep*corr(i)*  (          &
                       (constp(51) * (lambda(i)**cx(46))) /             &
                        (lambda_h1(i)**cx(61))     +                    &
                       (constp(52) * (lambda(i)**cx(46))) /             &
                        (lambda_h2(i)**cx(62))     )  )                 &
                               )
  
  
          ELSE
  
               !Diagnostic Abel and Shipway version
               !First need lambda, which for this case is
               !lambda = (1 /(rain fraction temp1)) to power of cx(48)
  
            lambda(i) = ( 1.0 / (rainfrac(i)*temp1(i)) )**cx(48)
  
            lambda_h1(i) = lambda(i)+cx(56)
            lambda_h2(i) = lambda(i)+cx(57)
  
               !Now we have lambda, should be the same as above.
               !Keeping prognostic and diagnostic separate for now so
               !I can change one without affecting the other.
  
            qclnew(i)= qcl(i) / (                                       &
                      cfliq(i)+(cfliq(i)*timestep*corr(i)*  (           &
                      (constp(51) * (lambda(i)**cx(46))) /              &
                       (lambda_h1(i)**cx(61))     +                     &
                      (constp(52) * (lambda(i)**cx(46))) /              &
                       (lambda_h2(i)**cx(62))     )  )                  &
                              )
  
  
          END IF ! l_mcr_qrain
  
        ELSE
  
               ! Use Sachinananda and Zrnic (1986) rain velocity formula
  
          IF (l_mcr_qrain) THEN
            qclnew(i)=qcl(i)/((cfliq(i)+cfliq(i)*constp(49)*corr(i)     &
                                     *timestep*temp1(i)**cx(53)))
          ELSE
            qclnew(i)=qcl(i)/((cfliq(i)+cfliq(i)*constp(49)*corr(i)     &
                                     *timestep*temp1(i)**cx(50)))
          END IF ! l_mcr_qrain
  
        END IF ! L_rainfall_as
  
            !-----------------------------------------------
            ! Convert qclnew to a gridbox mean
            !-----------------------------------------------
            ! temp7 is the rain free region of liquid cloud
        temp7(i) = MAX(cfliq(i)-rain_liq(i)-rain_mix(i),0.0)
        qclnew(i) = qclnew(i)*(rain_liq(i)+rain_mix(i))                 &
                   +qcl(i)/cfliq(i)*temp7(i)
  
            !-----------------------------------------------
            ! Calculate transfer
            !-----------------------------------------------
        dpr(i) = qcl(i) - qclnew(i)

        ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

          !------------------------------------------------
          ! Adjust liquid and rain contents
          !------------------------------------------------
        qrain(i) = qrain(i) + dpr(i)
        qcl(i)   = qcl(i)   - dpr(i)

          !------------------------------------------------
          ! Update cloud fractions
          !------------------------------------------------
!         These are commented out since there is currently no
!         cloud fraction update associated with the accretion term.
!          cf_transfer_rate(i)  = 0.0 / (timestep*iterations)
!          cfl_transfer_rate(i) = 0.0 / (timestep*iterations)


!            cftransfer(i)  = cftransfer(i)  + cf_transfer_rate(i)
!            cfltransfer(i) = cfltransfer(i) + cfl_transfer_rate(i)

!          
!            cf(i)  = cf(i)  +cf_transfer_rate(i) *timestep*iterations
!            cfl(i) = cfl(i) +cfl_transfer_rate(i)*timestep*iterations
!            

          !-----------------------------------------------
          ! Update rain fractions
          !-----------------------------------------------
!         These are commented out since there is currently no
!         rain fraction update associated with the accretion term.
!          rf_final(i) = rainfrac(i)
!          rftransfer(i) = rftransfer(i) + (rf_final(i) - rainfrac(i))

!         
!          rainfrac(i)= rf_final(i)
!          rain_liq(i)= min(area_liq(i),rainfrac(i))
!          rain_mix(i)= min(area_mix(i),rainfrac(i)-rain_liq(i))
!          rain_ice(i)=
!     &         min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
!          rain_clear(i)=
!     &         rainfrac(i) - rain_liq(i)-rain_mix(i)-rain_ice(i)
!          

      END IF  ! rainfrac gt 0 etc.

    END DO  ! Points

  END IF ! l_warm_new

  IF (lhook) CALL dr_hook('LSP_ACCRETION',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_accretion
END MODULE lsp_accretion_mod
