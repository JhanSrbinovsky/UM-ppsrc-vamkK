! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_solve_simeq-----------------------------------------
!
!  Purpose: To solve simultaneous equations for tsq, qsq and cov
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE mym_solve_simeq(                                             &
! IN levels
      bl_levels,                                                        &
! IN fields
      qq_tsq, qq_qsq, qq_cov, aa_tsq, bb_tsq, cc_tsq, pp_tc,            &
      aa_qsq, bb_qsq, cc_qsq, pp_qc,aa_cov,bb_cov, cc_cov, pp_ct, pp_cq,&
! OUT fields
      tsq, qsq, cov)

  USE atm_fields_bounds_mod, ONLY: tdims, tdims_l
  USE mym_option_mod, ONLY: tke_levels
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

! Intent IN Variables
  INTEGER, INTENT(IN) ::                                                &
     bl_levels
                   ! Max. no. of "boundary" level

  REAL, INTENT(IN) ::                                                   &
     ! matrix elements (for meanings of each, see the document)
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
            tke_levels),                                                &
     pp_cq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           tke_levels)

  REAL, INTENT(OUT) ::                                                  &
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

! Local variables
  INTEGER ::                                                            &
     i, j, k,                                                           &
     endflag

  REAL ::                                                               &
     ! one-dimensional variables to secure continuous memory accesses
     qq_tsq_k(tke_levels),                                              &
     qq_qsq_k(tke_levels),                                              &
     qq_cov_k(tke_levels),                                              &
     aa_tsq_k(tke_levels),                                              &
     bb_tsq_k(tke_levels),                                              &
     cc_tsq_k(tke_levels),                                              &
     pp_tc_k(tke_levels),                                               &
     aa_qsq_k(tke_levels),                                              &
     bb_qsq_k(tke_levels),                                              &
     cc_qsq_k(tke_levels),                                              &
     pp_qc_k(tke_levels),                                               &
     aa_cov_k(tke_levels),                                              &
     bb_cov_k(tke_levels),                                              &
     cc_cov_k(tke_levels),                                              &
     pp_ct_k(tke_levels),                                               &
     pp_cq_k(tke_levels),                                               &
     tsq_k(tke_levels),                                                 &
     qsq_k(tke_levels),                                                 &
     cov_k(tke_levels)

! Parameters
  INTEGER, PARAMETER ::                                                 &
     max_itr    = 500
               ! the maximum iteration number

  REAL, PARAMETER ::                                                    &
     eps       = 1.0e-15
               ! convergence creteria

  REAL, PARAMETER ::                                                    &
     tsq_scale = 1.0e0,                                                 &
     qsq_scale = 1.0e6,                                                 &
     cov_scale = 1.0e3,                                                 &
     r_tsq_scale = 1.0 / tsq_scale,                                     &
     r_qsq_scale = 1.0 / qsq_scale,                                     &
     r_cov_scale = 1.0 / cov_scale,                                     &
     tc_scale = tsq_scale * r_cov_scale,                                &
     qc_scale = qsq_scale * r_cov_scale,                                &
     ct_scale = cov_scale * r_tsq_scale,                                &
     cq_scale = cov_scale * r_qsq_scale
               ! scaling factors for the matrix elements

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

      ! Copy to 1dim variables to secure continuous memory accesses
      DO k = 1, tke_levels
        qq_tsq_k(k) = qq_tsq(i, j, k) * tsq_scale
        qq_qsq_k(k) = qq_qsq(i, j, k) * qsq_scale
        qq_cov_k(k) = qq_cov(i, j, k) * cov_scale
        aa_tsq_k(k) = aa_tsq(i, j, k)
        bb_tsq_k(k) = bb_tsq(i, j, k)
        cc_tsq_k(k) = cc_tsq(i, j, k)
        pp_tc_k(k)  = pp_tc(i, j, k)  * tc_scale
        aa_qsq_k(k) = aa_qsq(i, j, k)
        bb_qsq_k(k) = bb_qsq(i, j, k)
        cc_qsq_k(k) = cc_qsq(i, j, k)
        pp_qc_k(k)  = pp_qc(i, j, k)  * qc_scale
        aa_cov_k(k) = aa_cov(i, j, k)
        bb_cov_k(k) = bb_cov(i, j, k)
        cc_cov_k(k) = cc_cov(i, j, k)
        pp_ct_k(k)  = pp_ct(i, j, k)  * ct_scale
        pp_cq_k(k)  = pp_cq(i, j, k)  * cq_scale
      END DO

! DEPENDS ON: mym_solve_simeq_bcgstab
      CALL  mym_solve_simeq_bcgstab(                                    &
              max_itr, eps,                                             &
              qq_tsq_k, qq_qsq_k, qq_cov_k,                             &
              aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                    &
              aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                    &
              aa_cov_k, bb_cov_k, cc_cov_k,                             &
              pp_ct_k, pp_cq_k,                                         &
              tsq_k, qsq_k, cov_k, endflag)

      IF (endflag < 0) THEN
      ! if failed to converge, solve eqs. by LU decomposition
! DEPENDS ON: mym_solve_simeq_lud
        CALL mym_solve_simeq_lud(                                       &
              qq_tsq_k, qq_qsq_k, qq_cov_k,                             &
              aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                    &
              aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                    &
              aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,           &
              tsq_k, qsq_k, cov_k)
      END IF

      ! set the values into the original arrays.
      DO k = 1, tke_levels
        tsq(i, j, k) = tsq_k(k) * r_tsq_scale
        qsq(i, j, k) = qsq_k(k) * r_qsq_scale
        cov(i, j, k) = cov_k(k) * r_cov_scale
      END DO

    END DO !loop i = tdims%i_start, tdims%i_end
  END DO   !loop j = tdims%j_start, tdims%j_end
  IF (lhook) CALL dr_hook('MYM_SOLVE_SIMEQ',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE mym_solve_simeq
