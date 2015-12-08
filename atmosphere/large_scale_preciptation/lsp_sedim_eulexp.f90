! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Microphysics hydrometeor Eulerian sedimentation scheme
MODULE lsp_sedim_eulexp_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_sedim_eulexp(                                            &
  lsiter,points,m0,dhi,dhir,rho,rhor,                                   &
  flux_fromabove, fallspeed_thislayer,                                  &
  mixratio_thislayer, fallspeed_fromabove,                              &
  total_flux_out)


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

! Description:
!   First order Eulerian advection scheme for hydrometeor sedimentation
!   with an exponential-based limiter to ensure stability for CFL>1
!   used by routine LSP_ICE for ice crystals, snow, rain and graupel.

! Method:
!   Based on method described in Rotstayn (1997)(QJRMS, 123, 1227-1282)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.


! Subroutine arguments

      ! Intent (In)
  INTEGER ::                                                            &
    lsiter,                                                             &
                     ! number of iterations of microphysics param
    points       ! number of points to process

  REAL ::                                                               &
    m0,                                                                 &
                           ! Small mass (kg/kg) defined in c_lspmic
    dhi(points),                                                        &
                           ! CFL limit (s m-1)
    dhir(points),                                                       &
                           ! 1.0/DHI (m s-1)
    rho(points),                                                        &
                           ! Air density (kg m-3)
    rhor(points),                                                       &
                           ! 1.0/Rho
    flux_fromabove(points),                                             &
    fallspeed_thislayer(points)

      ! Intent (InOut)
  REAL ::                                                               &
    mixratio_thislayer(points),                                         &
    fallspeed_fromabove(points)

      ! Intent (Out)
  REAL ::                                                               &
    total_flux_out(points)

! Local variables

  REAL ::                                                               &
    mixratio_fromabove(points),                                         &
                                   ! Mixing Ratio from above
    flux_out(points),                                                   &
                                   ! Temporary flux out of layer
    expfactor(points)          ! Exponential Factor

  INTEGER :: i                    ! Loop counter

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------

  IF (lhook) CALL dr_hook('LSP_SEDIM_EULEXP',zhook_in,zhook_handle)
  DO i = 1, points

        ! -------------------------------------------------------------
        ! Adjust fall speeds to make a linear combination of the fall
        ! speed from the layer above. This ensures that the fall speed
        ! in this layer will not be calculated as zero even though
        ! there is mass falling into it from above.
        ! -------------------------------------------------------------

    mixratio_fromabove(i) = flux_fromabove(i)*dhi(i)*rhor(i)

    IF (mixratio_thislayer(i) + mixratio_fromabove(i)  >   m0) THEN

      fallspeed_thislayer(i) =                                          &
            (fallspeed_thislayer(i)*mixratio_thislayer(i)               &
            +fallspeed_fromabove(i)*mixratio_fromabove(i))              &
            /(mixratio_thislayer(i)+mixratio_fromabove(i))

    ELSE

      fallspeed_thislayer(i) = 0.0

    END IF

        ! -------------------------------------------------------------
        ! Eulerian solution with exponential limiter
        ! -------------------------------------------------------------

    IF (fallspeed_thislayer(i)  >   0.0) THEN

      expfactor(i) = EXP(-1.0*fallspeed_thislayer(i)*dhi(i))

          ! Calculate flux out of this layer

      flux_out(i) = flux_fromabove(i)+dhir(i)                           &
                   *(rho(i)*mixratio_thislayer(i)-flux_fromabove(i)     &
                   /fallspeed_thislayer(i))*(1.0-expfactor(i))

          ! Calculate mass (kg/kg) that remains in this layer

      mixratio_thislayer(i) = flux_fromabove(i)*rhor(i)                 &
                       / fallspeed_thislayer(i)                         &
                       * (1.0-expfactor(i)) + mixratio_thislayer(i)     &
                       * expfactor(i)

    ELSE

          ! No fall out of the layer.
          ! Set MixingRatio to be the amount of mass falling in.
          ! FallSpeed can only be zero if MixingRatio_ThisLayer LE M0
          ! and MixingRatio_FromAbove LE M0
          ! so this is slightly inconsistent.
      flux_out(i)       = 0.0
      mixratio_thislayer(i) = flux_fromabove(i)*rhor(i)*dhi(i)

    END IF

        ! No need to compute fall speed out of the layer in this method
        ! -------------------------------------------------------------
        !  Total_Flux_Out is the flux out of this layer that falls
        !  into the next layer down
        ! -------------------------------------------------------------

    total_flux_out(i) = total_flux_out(i) + flux_out(i)/lsiter

        ! -------------------------------------------------------------
        ! Store fall speed in this layer to be the fallspeed from above
        ! for the next layer down
        ! -------------------------------------------------------------

    fallspeed_fromabove(i) = fallspeed_thislayer(i)

  END DO ! on loop over points

  IF (lhook) CALL dr_hook('LSP_SEDIM_EULEXP',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_sedim_eulexp
END MODULE lsp_sedim_eulexp_mod
