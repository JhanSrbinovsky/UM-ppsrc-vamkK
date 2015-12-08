! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Autoconversion of snow.
MODULE lsp_snow_autoc_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_snow_autoc(                                              &
  points, timestep,                                                     &
                                          ! Number of points and tstep
  qcf_cry, qcf_agg, t, cttemp,                                          &
                                          ! Water contents and temp
  m0, t_scaling, qcf0,                                                  &
                                          ! Parametrization information
  ptransfer,                                                            &
                                          ! Mass transfer diagnostic
  one_over_tsi                                                          &
                                          ! 1/(timestep*iterations)
  )


  ! Dr Hook Modules
  USE yomhook,          ONLY: lhook, dr_hook
  USE parkind1,         ONLY: jprb, jpim

  IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of the autoconversion of
!   cloud ice to snow aggregates

!  Method:
!   Transfer some mass from ice crystals to snow aggregates depending
!   on the temperature and cloud top temperature.
!   Simple explicit Kessler type param. of autoconversion
!   with the autoconversion limit and rate set to emulate the
!   split-ice scheme.
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
                          ! Number of points to process

  REAL, INTENT(IN) ::                                                   &
    timestep,                                                           &
                          ! Timestep / s
    m0,                                                                 &
                          ! Seed ice mass / kg kg-1
    t_scaling,                                                          &
                          ! Scaling temperature / K
    qcf0,                                                               &
                          ! Prescribed ice content / kg kg-1
    t(points),                                                          &
                          ! Temperature / K
    cttemp(points),                                                     &
                          ! Cloud top temperature / K
    one_over_tsi
                          ! 1/(timestep*iterations)

  REAL, INTENT(INOUT) ::                                                &
    qcf_cry(points),                                                    &
                          ! Ice water content of ice crystals / kg kg-1
    qcf_agg(points),                                                    &
                          ! Ice water content of snow aggs. / kg kg-1
    ptransfer(points)
                          ! Autoconversion rate / kg kg-1 s-1

! Local Variables

  INTEGER ::                                                            &
    i


  REAL ::                                                               &
    dpr(points),                                                        &
                          ! Transfer amount from ice to snow / kg kg-1
    qcfautolim(points),                                                 &
                          ! Autoconversion limit / kg kg-1
    qcfautorate(points),                                                &
                          ! Rate of transfer / s-1
    qc(points)
                          ! Ice remaining after autoconversion
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('LSP_SNOW_AUTOC',zhook_in,zhook_handle)
  DO i = 1, points

    IF (qcf_cry(i) > m0) THEN

          !-----------------------------------------------
          ! Set autoconversion limit to emulate split-ice scheme
          !-----------------------------------------------
      qcfautolim(i) = (qcf_agg(i)+qcf_cry(i))                           &
                 *MAX(EXP(-t_scaling * MAX((t(i)-cttemp(i)),0.0)        &
                 *MAX(qcf_agg(i)+qcf_cry(i),0.0)*qcf0) , 0.0)

          !-----------------------------------------------
          ! Set rate to emulate spilt-ice scheme, i.e. infinite
          !-----------------------------------------------
      qcfautorate(i) = 1.0/timestep

      qc(i)  = MIN(qcfautolim(i) , qcf_cry(i))
      dpr(i) = MIN(qcfautorate(i) * timestep * (qcf_cry(i)-qc(i)),      &
                    qcf_cry(i) - qc(i))

          !-----------------------------------------------
          ! Store process rate (kg kg-1 s-1)
          !-----------------------------------------------
      ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

          !-----------------------------------------------
          ! Update ice/snow variables
          !-----------------------------------------------
      qcf_cry(i) = qcf_cry(i) - dpr(i)
      qcf_agg(i) = qcf_agg(i) + dpr(i)

          ! No cloud fraction updating is needed

    END IF  ! qcf_cry > 0

  END DO  ! Points

  IF (lhook) CALL dr_hook('LSP_SNOW_AUTOC',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_snow_autoc
END MODULE lsp_snow_autoc_mod
