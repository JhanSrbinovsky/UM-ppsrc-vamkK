! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Riming of ice particles
! Subroutine Interface:
MODULE lsp_riming_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_riming(                                                  &
  points, timestep,                                                     &
                                          ! Number of points and tstep
  qcl, qcf, t,                                                          &
                                          ! Water contents and temp
  area_liq, area_mix, cfliq, cficei,                                    &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
!    &, cf, cfl, cff                      ! Current cloud fractions for
!                                         ! updating
  rho, m0, tcg, tcgi, corr,                                             &
                                          ! Parametrization information
  lfrcp , ice_type,                                                     &
                                          ! Microphysical information
  l_psd,                                                                &
                                          ! Code options
  ptransfer,                                                            &
                                          ! Mass transfer diagnostic
  one_over_tsi                                                          &
                                          ! 1/(timestep*iterations)
!    &, cftransfer,cfltransfer,cfftransfer! Cloud transfer diagnostics
  )

  ! Microphysics modules
  USE mphys_constants_mod, ONLY: cx, constp, ice_type_offset

  ! General atmosphere modules
  USE conversions_mod,     ONLY: zerodegc

  ! Dr Hook modules
  USE yomhook,             ONLY: lhook, dr_hook
  USE parkind1,            ONLY: jprb, jpim

  ! Large scale precipitation modules
  USE lsp_moments_mod, ONLY: lsp_moments
  IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of ice particle riming

! Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles sweeping out liquid water.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! Riming is a source term for ice crystals
! and a sink term for supercooled water.
! Occurs if there are ice crystals and liquid water drops present
! in a mixed phase overlap region. The Numerical solution is implicit
! in liquid water.
! Riming does not, in this formulation, update cloud fractions
! which is why these variables have been commented out.

! Subroutine Arguments

  INTEGER, INTENT(IN) ::                                                &
    points,                                                             &
                          ! Number of points to calculate
    ice_type          ! Type of ice (0 - crystals, 1 - aggregates
                          !              3 - graupel)

  REAL, INTENT(IN) ::                                                   &
    timestep,                                                           &
                          ! Timestep / s
    area_liq(points),                                                   &
                          ! Fraction of gridbox with liquid-only cloud
    area_mix(points),                                                   &
                          ! Fraction of gridbox with mixed phase cloud
    cfliq(points),                                                      &
                          ! Fraction of gridbox with liquid cloud
    cficei(points),                                                     &
                          ! 1/Fraction of gridbox with ice cloud
!    &, cf(points)        ! Current cloud fraction
!    &, cfl(points)       ! Current liquid cloud fraction
!    &, cff(points)       ! Current ice cloud fraction
    rho(points),                                                        &
                          ! Air density / kg m-3
    m0,                                                                 &
                          ! Seed ice water content / kg kg-1
    tcg(points),                                                        &
                          ! T dependent function in ice size dist'n
    tcgi(points),                                                       &
                          ! 1/tcg (no units)
    corr(points),                                                       &
                          ! Fall velocity correction factor (no units)
    lfrcp,                                                              &
                          ! Latent heat of fusion
                          ! / heat capacity of air / K
    one_over_tsi
                          ! 1/(timestep*iterations)

  REAL, INTENT(INOUT) ::                                                &
    qcl(points),                                                        &
                          ! Liquid water content / kg kg-1
    qcf(points),                                                        &
                          ! Ice water content    / kg kg-1
    t(points),                                                          &
                          ! Temperature / K
    ptransfer(points) ! Mass rimed in this timestep / kg kg-1

!      Real, Intent(InOut) ::
!    &, cf_transfer_rate(points) ! Cloud fraction increment this tstep
!    &, cfl_transfer_rate(points)! Liquid cloud fraction inc this tstep
!    &, cff_transfer_rate(points)! Ice cloud fraction inc this tstep

  LOGICAL, INTENT(IN) ::                                                &
    l_psd
                          ! Use generic ice particle size distribution

! Local Variables

  INTEGER ::                                                            &
    i,                                                                  &
                          ! Loop counter for points
    cry_offset        ! Index offset for ice crystals

  REAL ::                                                               &
    qclnew(points),                                                     &
                          ! For value of qcl after riming  / kg kg-1
    dqi(points),                                                        &
                          ! Amount of ice rimed  / kg kg-1
    m_2_di(points)
                          ! 2+DI moment of particle size distribution

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


      !-----------------------------------------------
      ! Select appropriate ice parametrization (see c_lspsiz)
      !-----------------------------------------------
  IF (lhook) CALL dr_hook('LSP_RIMING',zhook_in,zhook_handle)

  cry_offset=ice_type*ice_type_offset

  IF (l_psd) THEN
        ! Use the generic ice particle size distribution
        ! Calculate the 2+di (cx(81)) moment of the
        ! ice particle size distribution.

    CALL lsp_moments(points,rho,t,qcf,cficei,cx(81),m_2_di)

  END IF

  DO i = 1, points

    IF (qcf(i) >  m0 .AND. qcl(i) >  0.0 .AND. t(i)  <   zerodegc       &
        .AND. area_mix(i) >  0.0 .AND. cfliq(i) >  0.0) THEN

          !-----------------------------------------------
          ! Calculate water content of mixed phase region
          !-----------------------------------------------
      IF (l_psd) THEN

            ! Calculate the riming rate using the generic PSD
            ! constp(81) = (pi/4) ci
        qclnew(i) = qcl(i) /                                            &
                  (cfliq(i)+cfliq(i)*constp(81)                         &
                  *corr(i)*timestep*m_2_di(i))

      ELSE
            ! Use the defined gamma distribution

        qclnew(i) = qcl(i) /                                            &
                  (cfliq(i)+cfliq(i)*constp(9+cry_offset)*tcg(i)        &
                  *corr(i)*timestep*(rho(i)*qcf(i)*cficei(i)            &
                  *constp(5+cry_offset)*tcgi(i))**cx(6+cry_offset))
      END IF

          !-----------------------------------------------
          ! Convert to new grid box total water content
          !-----------------------------------------------
      qclnew(i)=qcl(i)*area_liq(i)/cfliq(i)+qclnew(i)*area_mix(i)
      dqi(i)=(qcl(i)-qclnew(i))


          !-----------------------------------------------
          ! Store process rate / kg kg-1 s-1 and cloud fraction changes
          !-----------------------------------------------
      ptransfer(i) = ptransfer(i) + dqi(i) * one_over_tsi

!         There are no cloud fraction updates associated
!         with the riming term. They will have already been set to
!         zero on input to this subroutine, which is why these
!         are commented out.
!          cf_transfer_rate(i)  = 0.0 / (timestep*iterations)
!          cfl_transfer_rate(i) = 0.0 / (timestep*iterations)
!          cff_transfer_rate(i) = 0.0 / (timestep*iterations)

          !-----------------------------------------------
          ! Update water contents
          !-----------------------------------------------
      qcf(i) = qcf(i) + dqi(i)
      t(i)   = t(i) + lfrcp * dqi(i)
      qcl(i) = qclnew(i)

          !-----------------------------------------------
          ! Update cloud fractions
          !-----------------------------------------------
!          IF (ice_type /= 3) THEN
!            These are commented out since there is currently no
!            cloud fraction update associated with the riming term.

!            cf(i)  = cf(i) + cf_transfer_rate(i) *timestep*iterations
!            cfl(i) = cfl(i)+ cfl_transfer_rate(i)*timestep*iterations
!            cff(i) = cff(i)+ cff_transfer_rate(i)*timestep*iterations

!            cftransfer(i)  = cftransfer(i)  + cf_transfer_rate(i)
!            cfltransfer(i) = cfltransfer(i) + cfl_transfer_rate(i)
!            cfftransfer(i) = cfftransfer(i) + cff_transfer_rate(i)

!          END IF  ! ice_type /= 3

    END IF ! qcf(i) >  m0 etc.

  END DO  ! Points

  IF (lhook) CALL dr_hook('LSP_RIMING',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_riming
END MODULE lsp_riming_mod
