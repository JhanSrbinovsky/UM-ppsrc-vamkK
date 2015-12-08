! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Capture of one ice category by
!  another.
MODULE lsp_collection_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_collection(                                              &
  points, timestep,                                                     &
                                          ! Number of points and tstep
  qcf1, qcf2, t,                                                        &
                                          ! Water contents and temp
  area_mix, area_ice,                                                   &
  cficei,                                                               &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
!    &, cf, cff                           ! Current cloud fractions for
!                                         ! updating
  rho, rhor, m0, tcg1, tcg1i, tcg2, tcg2i,                              &
                                          ! Parametrization information
  corr, ice_type1, ice_type2,                                           &
                                          ! Microphysical information
  l_psd_2 ,  l_use_area,                                                &
  l_no_t_check,                                                         &
                                          ! Code options
  ptransfer,                                                            &
                                          ! Mass transfer diagnostic
  one_over_tsi                                                          &
                                          ! 1/(timestep*iterations)
!    &, cftransfer, cfftransfer           ! Cloud transfer diagnostics
  )

  ! Microphysics modules
  USE mphys_constants_mod, ONLY: cx, constp, ice_type_offset

  ! General atmosphere modules
  USE conversions_mod,     ONLY: zerodegc

  ! Dr Hook Modules
  USE yomhook,             ONLY: lhook, dr_hook
  USE parkind1,            ONLY: jprb, jpim

  ! Large scale precipitation modules
  USE lsp_moments_mod,     ONLY: lsp_moments

  IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of collisions between different
!   categories of ice particles

!  Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles (species 1) sweeping out a different
!   specified distribution of ice particles (species 2) and adding
!   combining the mass into species 1.
!   Note that the current formulation is limited to the collection of
!   snow by graupel or the collection of ice by snow. If we want to do
!   graupel to ice then we need to look again at the prefactors
!   constp(16 to 19 +ice_offset1) because they depend on both the
!   species, not just species 1.
!   The commented out code referring to cloud fraction updates is to
!   highlight where in the subroutine changes would need to be made
!   if it was later decided to update cloud fractions in some way.
!   We note that we are only allowing use of the generic ice particle
!   size distribution for second species, intending its use for the
!   collection of aggregates (generic psd) by graupel and *not* the
!   collection of crystals by aggregates (this would go against the
!   idea that the generic psd represents both the crystals and the
!   aggregates). The logic in the call to the routine should prevent
!   any inconsistencies since l_mcr_qcf2 and l_psd should not both be
!   true.
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
    ice_type1,                                                          &
    ice_type2
                          ! Type of ice (0 - crystals, 1 - aggregates
                          !              3 - graupel)

  REAL, INTENT(IN) ::                                                   &
    timestep,                                                           &
                          ! Timestep / s
    area_mix(points),                                                   &
                          ! Fraction of gridbox with mixed phase cloud
    area_ice(points),                                                   &
                          ! Fraction of gridbox with ice-only cloud
    cficei(points),                                                     &
                          ! 1/Fraction of gridbox with ice cloud
!    &, cf(points)        ! Current cloud fraction
!    &, cff(points)       ! Current ice cloud fraction
    rho(points),                                                        &
                          ! Air density / kg m-3
    rhor(points),                                                       &
                          ! 1/Air density / m3 kg-1
    m0,                                                                 &
                          ! Seed ice water content / kg kg-1
    tcg1(points),                                                       &
                          ! T dependent function in ice size dist'n
    tcg1i(points),                                                      &
                          ! 1/tcg (no units)
    tcg2(points),                                                       &
                          ! T dependent function in ice size dist'n
    tcg2i(points),                                                      &
                          ! 1/tcg (no units)
    corr(points),                                                       &
                          ! Fall velocity correction factor (no units)
    one_over_tsi,                                                       &
                          ! 1/(timestep*iterations)
    t(points)
                          ! Temperature (k)

  REAL, INTENT(INOUT) ::                                                &
    qcf1(points),                                                       &
                          ! Ice water content of 1st species / kg kg-1
    qcf2(points),                                                       &
                          ! Ice water content of 2nd species / kg kg-1
    ptransfer(points)
                          ! Mass rimed in this timestep / kg kg-1
!    &, cftransfer(points) ! Cumulative cloud fraction increment
!    &, cfftransfer(points)! Cumulative ice cloud fraction increment

  LOGICAL, INTENT(IN) ::                                                &
    l_psd_2,                                                            &
                          ! Use generic ice particle size distribution
    l_use_area,                                                         &
                          ! Use ice area amount in calculating
                          ! gridbox mean transfer rate
    l_no_t_check
                          ! Do not check that temperature is less
                          ! than 0 deg C before proceeding

! Local Variables

  INTEGER ::                                                            &
    i,                                                                  &
                          ! Loop counter for points
    ice_offset1,                                                        &
                          ! Index offset for first ice species
    ice_offset2
                          ! Index offset for second ice species

  REAL ::                                                               &
    dqi(points),                                                        &
                          ! Transfer of mixing ratio  / kg kg-1
    vi1(points),                                                        &
                          ! Average fall speed of first ice
                          ! particle species  / m s-1
    vi2(points),                                                        &
                          ! Average fall speed of second ice
                          ! particle species  / m s-1
    fv1(points),                                                        &
                          ! Average fall speed difference between
                          ! ice species  / m s-1
    lami1(points),                                                      &
                          ! Reciprocal of slope parameter in
                          ! first ice particle size distribution  / m
    lami2(points),                                                      &
                          ! Reciprocal of slope parameter in
                          ! second ice particle size distribution  / m
    lamfac1(points),                                                    &
                          ! Combination of lamr1 and lamr2
    collision_eff(points),                                              &
                          ! Collision efficiency / no units
    m_0(points), m_1(points), m_2(points),                              &
                          ! zero, 1st and 2nd moments of the ice PSD
    m_bi_di(points)
                          ! bi+di moment of the generic ice size distn

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


      !-----------------------------------------------
      ! Select appropriate ice parametrization (see c_lspsiz)
      !-----------------------------------------------

  IF (lhook) CALL dr_hook('LSP_COLLECTION',zhook_in,zhook_handle)

  ice_offset1=ice_type1*ice_type_offset
  ice_offset2=ice_type2*ice_type_offset

      !-----------------------------------------------
      ! Calculate moments of size distribution if appropriate
      !-----------------------------------------------
  IF (l_psd_2) THEN
        ! Calculate the 0th, 1st, 2nd and bi+di (cx(82)) moments of the
        ! ice particle size distribution
    CALL lsp_moments(points,rho,t,qcf2,cficei,0.0,m_0)
    CALL lsp_moments(points,rho,t,qcf2,cficei,1.0,m_1)
    CALL lsp_moments(points,rho,t,qcf2,cficei,2.0,m_2)
    CALL lsp_moments(points,rho,t,qcf2,cficei,cx(82),m_bi_di)
  END IF  ! l_psd_2

  DO i = 1, points

    IF (qcf1(i)  >   m0 .AND. qcf2(i)  >   m0                           &
        .AND. ( t(i) < zerodegc .OR. l_no_t_check)                      &
        .AND. ( (area_ice(i)+area_mix(i))  >   0.0                      &
               .OR. (.NOT. l_use_area) ) ) THEN

          !-----------------------------------------------
          ! Calculate first ice mass-weighted fallspeed
          !-----------------------------------------------
          ! Use size distribution based on intercepts
      vi1(i) = constp(4+ice_offset1) * corr(i) *                        &
             ( rho(i) * qcf1(i) * cficei(i)                             &
           * constp(5+ice_offset1) * tcg1i(i))**cx(3+ice_offset1)

          !-----------------------------------------------
          ! Calculate second ice mass-weighted fallspeed
          !-----------------------------------------------
      IF (l_psd_2) THEN
            ! Use generic ice size distribution
        vi2(i) = constp(82) * corr(i) * m_bi_di(i)                      &
                   / (rho(i) * qcf1(i) * cficei(i))
      ELSE
            ! Use size distribution based on intercepts
        vi2(i) = constp(4+ice_offset2) * corr(i) *                      &
               ( rho(i) * qcf2(i) * cficei(i)                           &
             * constp(5+ice_offset2) * tcg2i(i))**cx(3+ice_offset2)
      END IF  ! l_psd_2

          !-----------------------------------------------
          ! Estimate the mean absolute differences in velocities
          !-----------------------------------------------
      fv1(i) = MAX(ABS(vi1(i)-vi2(i)),(vi1(i)+vi2(i))/8.0)

          !-----------------------------------------------
          ! Calculate reciprocal of lambda for first ice distribution
          !-----------------------------------------------
      lami1(i) = (rho(i) * qcf1(i) * cficei(i)                          &
        *constp(5+ice_offset1)*tcg1i(i))**(-cx(7+ice_offset1))

          !-----------------------------------------------
          ! Calculate reciprocal of lambda for second ice distribution
          !-----------------------------------------------
      IF (.NOT. l_psd_2) THEN
        lami2(i) = (rho(i) * qcf2(i) * cficei(i)                        &
          *constp(5+ice_offset2)*tcg2i(i))**(-cx(7+ice_offset2))
      END IF  ! l_psd_2

          !------------------------------------------------
          ! Calculate transfer
          !------------------------------------------------
      collision_eff(i) = 0.02*EXP(0.08*(t(i)-zerodegc))

      IF (l_psd_2) THEN
            ! Use the generic ice particle size distribution
            ! constp(80)=pi**2/24 x1g ag (80 = 20+ice_offset for graup)
            ! constp(91)=gamma(bi+3+x4i)
            ! constp(92)=2 gamma(bi+2+x4i)
            ! constp(93)=gamma(bi+1+x4i)
        dqi(i) = collision_eff(i) * fv1(i) *  timestep * rhor(i) *      &
                 tcg1(i) * constp(20+ice_offset1) *                     &
                 lami1(i)**(-cx(11+ice_offset1))  *                     &
              (constp(91) * m_0(i) * lami1(i)**cx(13+ice_offset1) +     &
               constp(92) * m_1(i) * lami1(i)**cx(14+ice_offset1) +     &
               constp(93) * m_2(i) * lami1(i)**cx(15+ice_offset1))
      ELSE
            ! Use the distribution defined by intercepts
        lamfac1(i) =                                                    &
           constp(16+ice_offset1)*(lami2(i)**cx(8+ice_offset2)          &
                                 * lami1(i)**cx(13+ice_offset1)) +      &
           constp(17+ice_offset1)*(lami2(i)**cx(9+ice_offset2)          &
                                 * lami1(i)**cx(14+ice_offset1)) +      &
           constp(18+ice_offset1)*(lami2(i)**cx(10+ice_offset2)         &
                                 * lami1(i)**cx(15+ice_offset1))

        dqi(i) = collision_eff(i) * tcg1(i) * tcg2(i)                   &
                 * constp(19+ice_offset1)                               &
                 * lami1(i)**(-cx(11+ice_offset1))                      &
                 * lami2(i)**(-cx(11+ice_offset2))                      &
                 * fv1(i) * lamfac1(i) * timestep * rhor(i)

      END IF  ! l_psd_2

      IF (l_use_area) THEN
            ! The calculations above have used in-cloud quantities to
            ! obtain the transfer rates. Multiply the transfer by the
            ! ice cloud amount to get a gridbox mean.
        dqi(i) = dqi(i) * (area_mix(i) + area_ice(i))
      END IF

          !------------------------------------------------
          ! Limit transfer to the mass of species 2 that is available
          !------------------------------------------------
      dqi(i) = MIN(dqi(i),qcf2(i))

      ptransfer(i) = ptransfer(i) + dqi(i) * one_over_tsi

          !------------------------------------------------
          ! Adjust ice species contents
          !------------------------------------------------
      qcf1(i)  = qcf1(i)   + dqi(i)
      qcf2(i)  = qcf2(i)   - dqi(i)

          !------------------------------------------------
          ! Update cloud fractions
          !------------------------------------------------
!         These are commented out since there is currently no
!         cloud fraction update associated with the collision terms.
!          cf_transfer_rate(i)  = 0.0 / (timestep*iterations)
!          cff_transfer_rate(i) = 0.0 / (timestep*iterations)

!            cftransfer(i)  = cftransfer(i)  + cf_transfer_rate(i)
!            cfftransfer(i) = cfftransfer(i) + cff_transfer_rate(i)

!            cf(i)  = cf(i)  +cf_transfer_rate(i) *timestep*iterations
!            cff(i) = cff(i) +cff_transfer_rate(i)*timestep*iterations


    END IF ! qcf1(i) >  m0 etc.

  END DO  ! Points

  IF (lhook) CALL dr_hook('LSP_COLLECTION',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_collection
END MODULE lsp_collection_mod
