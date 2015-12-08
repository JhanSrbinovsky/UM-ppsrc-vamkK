! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Large-scale precipitation scheme. Autoconversion of graupel.
! Subroutine Interface:
MODULE lsp_graup_autoc_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE lsp_graup_autoc(                                             &
  points, timestep,                                                     &
                                          ! Number of points and tstep
  qcf_agg, qgraup, t, rho,                                              &
                                          ! Water contents and temp
  psacw, psdep, ptransfer,                                              &
                                          ! Mass transfer diagnostic
  one_over_tsi                                                          &
                                          ! 1/(timestep*iterations)
  )

! Microphysics modules
  USE lsp_autoc_consts_mod,ONLY: auto_graup_qcf_thresh,                 &
                                 auto_graup_t_thresh, auto_graup_coeff 

! General modules
  USE conversions_mod,     ONLY: zerodegc
  
! Dr Hook Modules
  USE yomhook,             ONLY: lhook, dr_hook
  USE parkind1,            ONLY: jprb, jpim
 
  IMPLICIT NONE

! Purpose:
!   Update cloud prognostics as a result of the autoconversion of
!   snow aggregates to graupel.

!  Method:
!   This is a source term for graupel when snow growth is dominated
!   by riming liquid water cloud and is the only term that can create
!   graupel where there was no graupel before.
!   The conversion rate is proportional to how much the riming rate
!   of aggregates (PSACW) exceeds the rate of growth due to vapour
!   deposition (PSDEP) and collection/accretion of ice crystals (PSACI).
!   The coefficient reduces the conversion rate as the riming
!   aggregates will not immediately increase their density to that of
!   graupel.
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
    one_over_tsi,                                                       &
                          ! 1/(timestep*iterations)
    t(points),                                                          &
                          ! Temperature / K
    rho(points),                                                        &
                          ! Air density / kg kg s-1
    psacw(points),                                                      &
                          ! Riming rate of snow aggs. / kg kg-1 s-1
    psdep(points)
                          ! Deposition rate of snow aggs. / kg kg-1 s-1

  REAL, INTENT(INOUT) ::                                                &
    qgraup(points),                                                     &
                          ! Graupel content / kg kg-1
    qcf_agg(points),                                                    &
                          ! Ice water content of snow aggs. / kg kg-1
    ptransfer(points)
                          ! Autoconversion rate / kg kg-1 s-1

! Local Variables

  INTEGER ::                                                            &
    i


  REAL ::                                                               &
    dqi(points)
                          ! Transfer amount from ice to snow / kg kg-1

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle



  IF (lhook) CALL dr_hook('LSP_GRAUP_AUTOC',zhook_in,zhook_handle)
  DO i = 1, points

    IF (rho(i)*qcf_agg(i)  >   auto_graup_qcf_thresh                    &
        .AND. t(i)  <  (zerodegc+auto_graup_t_thresh)) THEN

          !-----------------------------------------------
          ! Calculate conversion rate / kg kg-1 s-1
          !-----------------------------------------------
      dqi(i) = auto_graup_coeff *                                       &
               MAX(0.0, psacw(i)-psdep(i))

      dqi(i) = MIN(dqi(i)*timestep,qcf_agg(i))

          !-----------------------------------------------
          ! Store process rate / kg kg-1 s-1
          !-----------------------------------------------
      ptransfer(i) = ptransfer(i) + dqi(i) * one_over_tsi

          !-----------------------------------------------
          ! Recalculate water contents
          !-----------------------------------------------
      qgraup(i)  = qgraup(i)  + dqi(i)
      qcf_agg(i) = qcf_agg(i) - dqi(i)

          ! No cloud fraction updates

    END IF !  qcf_agg > threshold etc.

  END DO  ! Points

  IF (lhook) CALL dr_hook('LSP_GRAUP_AUTOC',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lsp_graup_autoc
END MODULE lsp_graup_autoc_mod
