! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_simeq_matrix_prod-----------------------------------
!
!  Purpose: To calculate products of a matrix and a vector
!           in solving the simultaneous equations.
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE mym_simeq_matrix_prod(                                       &
      aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                            &
      aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                            &
      aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,                   &
      x_tsq_k, x_qsq_k, x_cov_k,                                        &
      y_tsq_k, y_qsq_k, y_cov_k)
  USE mym_option_mod, ONLY: tke_levels
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

  REAL, INTENT(IN) ::                                                   &
     ! matrix elements (for meanings of each, see the document)
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
     ! vector elements
     x_tsq_k(tke_levels),                                               &
     x_qsq_k(tke_levels),                                               &
     x_cov_k(tke_levels)

  REAL, INTENT(OUT) ::                                                  &
     ! vector elements of products (answers)
     y_tsq_k(tke_levels),                                               &
     y_qsq_k(tke_levels),                                               &
     y_cov_k(tke_levels)

  INTEGER :: k

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_SIMEQ_MATRIX_PROD',                      &
                                     zhook_in,zhook_handle)

  ! y = A * x
  k = 1
  y_tsq_k(k) =      bb_tsq_k(k) * x_tsq_k(k)                            &
                  + cc_tsq_k(k) * x_tsq_k(k + 1)                        &
                  + pp_tc_k(k) * x_cov_k(k)

  y_qsq_k(k) =      bb_qsq_k(k) * x_qsq_k(k)                            &
                  + cc_qsq_k(k) * x_qsq_k(k + 1)                        &
                  + pp_qc_k(k) * x_cov_k(k)

  y_cov_k(k) =     bb_cov_k(k) * x_cov_k(k)                             &
                  + cc_cov_k(k) * x_cov_k(k + 1)                        &
                  + pp_ct_k(k) * x_tsq_k(k)                             &
                  + pp_cq_k(k) * x_qsq_k(k)

  DO k = 2, tke_levels - 1
    y_tsq_k(k) = aa_tsq_k(k) * x_tsq_k(k - 1)                           &
                  + bb_tsq_k(k) * x_tsq_k(k)                            &
                  + cc_tsq_k(k) * x_tsq_k(k + 1)                        &
                  + pp_tc_k(k) * x_cov_k(k)

    y_qsq_k(k) = aa_qsq_k(k) * x_qsq_k(k - 1)                           &
                  + bb_qsq_k(k) * x_qsq_k(k)                            &
                  + cc_qsq_k(k) * x_qsq_k(k + 1)                        &
                  + pp_qc_k(k) * x_cov_k(k)

    y_cov_k(k) = aa_cov_k(k) * x_cov_k(k - 1)                           &
                  + bb_cov_k(k) * x_cov_k(k)                            &
                  + cc_cov_k(k) * x_cov_k(k + 1)                        &
                  + pp_ct_k(k) * x_tsq_k(k)                             &
                  + pp_cq_k(k) * x_qsq_k(k)

  END DO

  k = tke_levels
  y_tsq_k(k) = aa_tsq_k(k) * x_tsq_k(k - 1)                             &
                  + bb_tsq_k(k) * x_tsq_k(k)                            &
                  + pp_tc_k(k) * x_cov_k(k)

  y_qsq_k(k) = aa_qsq_k(k) * x_qsq_k(k - 1)                             &
                  + bb_qsq_k(k) * x_qsq_k(k)                            &
                  + pp_qc_k(k) * x_cov_k(k)

  y_cov_k(k) = aa_cov_k(k) * x_cov_k(k - 1)                             &
                  + bb_cov_k(k) * x_cov_k(k)                            &
                  + pp_ct_k(k) * x_tsq_k(k)                             &
                  + pp_cq_k(k) * x_qsq_k(k)


  IF (lhook) CALL dr_hook('MYM_SIMEQ_MATRIX_PROD',                      &
                                     zhook_out,zhook_handle)
  RETURN

END SUBROUTINE mym_simeq_matrix_prod
