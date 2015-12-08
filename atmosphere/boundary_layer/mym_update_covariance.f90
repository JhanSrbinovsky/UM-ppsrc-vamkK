! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_update_covariance-----------------------------------
!
!  Purpose: To integrate the covariances(tsq, qsq, cov) appeared
!           in the MY model.
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE mym_update_covariance(                                       &
! IN levels
      bl_levels,                                                        &
! IN fields
      qkw, el, dfm, pdt_tsq, pdt_cov, pdt_res,                          &
      pdq_qsq, pdq_cov, pdq_res, pdc_cov, pdc_tsq, pdc_qsq, pdc_res,    &
! INOUT fields
      tsq, qsq, cov)

  USE atm_fields_bounds_mod, ONLY: tdims, tdims_s, tdims_l
  USE timestep_mod, ONLY: timestep
  USE mym_const_mod, ONLY: b2, coef_trbvar_diff
  USE mym_option_mod, ONLY: l_my_extra_level, tke_levels
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

! Intent IN Variables
  INTEGER, INTENT(IN) ::                                                &
     bl_levels
                   ! Max. no. of "boundary" level

  REAL, INTENT(IN) ::                                                   &
     qkw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                   ! sqrt(qke) = sqrt(2TKE)
     el(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! mixing length
     dfm(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,   &
                                                    bl_levels),         &
                   ! diffusion coefficients fot momentum
     pdt_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to tsq in the production term of tsq
     pdt_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to cov in the production term of tsq
     pdt_res(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a residual part in the production term of tsq
     pdq_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to qsq in the production term of qsq
     pdq_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to cov in the production term of qsq
     pdq_res(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a residual part in the production term of qsq
     pdc_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to cov in the production term of cov
     pdc_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to tsq in the production term of cov
     pdc_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to qsq in the production term of cov
     pdc_res(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels)
                   ! a residual part in the production term of cov

  REAL, INTENT(INOUT) ::                                                &
     tsq(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
         bl_levels),                                                    &
                   ! Self covariance of liquid potential temperature
                   ! (thetal'**2) defined on theta levels K-1
     qsq(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
         bl_levels),                                                    &
                   ! Self covariance of total water
                   ! (qw'**2) defined on theta levels K-1
     cov(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
         bl_levels)
                   ! Correlation between thetal and qw
                   ! (thetal'qw') defined on theta levels K-1

! Local Variables
  INTEGER ::                                                            &
     i, j, k, k_start
                  ! loop indexes, etc.
  REAL ::                                                               &
     elem
                  ! work variables

  REAL ::                                                               &
     disp_coef
                  ! coefficients of the prognostic variables in
                  ! dissipation terms

  REAL ::                                                               &
     aa(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
     bb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
     cc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                  ! tri-diagonal matrix elements due to diffusion
     qq_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
     qq_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
     qq_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
     aa_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
     bb_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
     cc_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
     pp_tc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           tke_levels),                                                 &
     aa_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
     bb_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
     cc_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
     pp_qc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           tke_levels),                                                 &
     aa_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
     bb_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
     cc_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
     pp_ct(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           tke_levels),                                                 &
     pp_cq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           tke_levels)
                  ! matrix elements (see the documents for details)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_UPDATE_COVARIANCE',                      &
                                      zhook_in,zhook_handle)

! DEPENDS ON: mym_diff_matcoef
  CALL mym_diff_matcoef(                                                &
        bl_levels,coef_trbvar_diff, dfm, aa, bb, cc)

  IF (l_my_extra_level) THEN
    k_start = 1
  ELSE
    k_start = 2
  END IF

  ! set maxtrix elements
  DO k = k_start, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        disp_coef = 2.0 * qkw(i, j, k) / (b2 * el(i, j, k))
        elem = 1.0 - bb(i, j, k) * timestep                             &
                               + timestep * disp_coef
        bb_tsq(i, j, k) = elem                                          &
                            - 2.0 * pdt_tsq(i, j, k) * timestep
        bb_qsq(i, j, k) = elem                                          &
                            - 2.0 * pdq_qsq(i, j, k) * timestep
        bb_cov(i, j, k) = elem                                          &
                            - 2.0 * pdc_cov(i, j, k) * timestep
        qq_tsq(i, j, k) = tsq(i, j, k)                                  &
                          + 2.0 * pdt_res(i, j, k) * timestep
        qq_qsq(i, j, k) = qsq(i, j, k)                                  &
                          + 2.0 * pdq_res(i, j, k) * timestep
        qq_cov(i, j, k) = cov(i, j, k)                                  &
                          + 2.0 * pdc_res(i, j, k) * timestep

        elem = -aa(i, j, k) * timestep
        aa_tsq(i, j, k) = elem
        aa_qsq(i, j, k) = elem
        aa_cov(i, j, k) = elem

        elem = -cc(i, j, k) * timestep
        cc_tsq(i, j, k) = elem
        cc_qsq(i, j, k) = elem
        cc_cov(i, j, k) = elem

        pp_tc(i, j, k) = - 2.0 * pdt_cov(i, j, k) * timestep
        pp_qc(i, j, k) = - 2.0 * pdq_cov(i, j, k) * timestep
        pp_ct(i, j, k) = - 2.0 * pdc_tsq(i, j, k) * timestep
        pp_cq(i, j, k) = - 2.0 * pdc_qsq(i, j, k) * timestep
      END DO
    END DO
  END DO

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      aa_tsq(i, j, k_start) = 0.0
      aa_qsq(i, j, k_start) = 0.0
      aa_cov(i, j, k_start) = 0.0

      cc_tsq(i, j, tke_levels) = 0.0
      cc_qsq(i, j, tke_levels) = 0.0
      cc_cov(i, j, tke_levels) = 0.0
    END DO
  END DO

  IF (.NOT. l_my_extra_level) THEN
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        bb_tsq(i, j, 1) = 1.0
        bb_qsq(i, j, 1) = 1.0
        bb_cov(i, j, 1) = 1.0
        qq_tsq(i, j, 1) = 0.0
        qq_qsq(i, j, 1) = 0.0
        qq_cov(i, j, 1) = 0.0

        aa_tsq(i, j, 1) = 0.0
        aa_qsq(i, j, 1) = 0.0
        aa_cov(i, j, 1) = 0.0

        cc_tsq(i, j, 1) = 0.0
        cc_qsq(i, j, 1) = 0.0
        cc_cov(i, j, 1) = 0.0

        pp_tc(i, j, 1) = 0.0
        pp_qc(i, j, 1) = 0.0
        pp_ct(i, j, 1) = 0.0
        pp_cq(i, j, 1) = 0.0
      END DO
    END DO
  END IF

  ! Solve the simultaneous equations for tsq, qsq and cov
! DEPENDS ON: mym_solve_simeq
  CALL mym_solve_simeq(                                                 &
! IN levels
        bl_levels,                                                      &
! IN fields
        qq_tsq, qq_qsq, qq_cov, aa_tsq, bb_tsq, cc_tsq, pp_tc,          &
        aa_qsq, bb_qsq, cc_qsq, pp_qc,aa_cov,bb_cov,cc_cov,pp_ct, pp_cq,&
! OUT fields
        tsq, qsq, cov)

  IF (lhook) CALL dr_hook('MYM_UPDATE_COVARIANCE',                      &
                                      zhook_out,zhook_handle)
  RETURN

END SUBROUTINE mym_update_covariance
