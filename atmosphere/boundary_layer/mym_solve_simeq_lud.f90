! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_solve_simeq_lud-------------------------------------
!
!  Purpose: To solve simultaneous equations by LU decomposition
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE mym_solve_simeq_lud(                                         &
      qq_tsq_k, qq_qsq_k, qq_cov_k,                                     &
      aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                            &
      aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                            &
      aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,                   &
      tsq_k, qsq_k, cov_k)

  USE mym_option_mod, ONLY: tke_levels
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

! intent in variables
  REAL, INTENT(IN) ::                                                   &
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
     pp_cq_k(tke_levels)
             ! matrix elements

  REAL, INTENT(OUT) ::                                                  &
     tsq_k(tke_levels),                                                 &
     qsq_k(tke_levels),                                                 &
     cov_k(tke_levels)
             ! solved tsq, qsq and cov

  INTEGER ::                                                            &
     k, l, m, n,                                                        &
             ! loop indexes
     kpiv
             ! index of a pivot

  REAL ::                                                               &
     wk
            ! work variables

  REAL ::                                                               &
     amat(3 * tke_levels, 3 * tke_levels),                              &
             ! coefficient matrix
     bvec(3 * tke_levels)
             ! vector in the right hand side

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_SOLVE_SIMEQ_LUD',zhook_in,zhook_handle)

  amat(:, :) = 0.0

  DO k = 1, tke_levels
    amat(k, k) = bb_tsq_k(k)
    amat(tke_levels + k, tke_levels + k) = bb_qsq_k(k)
    amat(2 * tke_levels + k, 2 * tke_levels + k)                        &
                                          = bb_cov_k(k)
    bvec(k)    = qq_tsq_k(k)
    bvec(tke_levels + k) = qq_qsq_k(k)
    bvec(2 * tke_levels + k) = qq_cov_k(k)
  END DO

  DO k = 2, tke_levels
    amat(k, k-1) = aa_tsq_k(k)
    amat(tke_levels + k, tke_levels + k - 1) = aa_qsq_k(k)
    amat(2 * tke_levels + k, 2 * tke_levels + k - 1)                    &
                                          = aa_cov_k(k)
  END DO

  DO k = 1, tke_levels - 1
    amat(k, k+1) = cc_tsq_k(k)
    amat(tke_levels + k, tke_levels + k + 1) = cc_qsq_k(k)
    amat(2 * tke_levels + k, 2 * tke_levels + k + 1)                    &
                                          = cc_cov_k(k)
  END DO

  DO k = 1, tke_levels
    amat(k, 2 * tke_levels + k) = pp_tc_k(k)
    amat(tke_levels + k, 2 * tke_levels + k) = pp_qc_k(k)
    amat(2 * tke_levels + k, k) = pp_ct_k(k)
    amat(2 * tke_levels + k, tke_levels + k) = pp_cq_k(k)
  END DO

  n = 3 * tke_levels
  ! main part
  DO k = 1, n
    kpiv = k
    wk  = ABS(amat(k, k))
    DO l = k + 1, n
      IF(ABS(amat(l, k)) > wk) THEN
        kpiv = l
        wk  = ABS(amat(l, k))
      END IF
    END DO

    IF(kpiv /= k) THEN
      DO m = 1, n
        wk       = amat(k, m)
        amat(k, m)    = amat(kpiv, m)
        amat(kpiv, m) = wk
      END DO
      wk   = bvec(k)
      bvec(k) = bvec(kpiv)
      bvec(kpiv) = wk
    END IF

    amat(k, k) = 1.0 / amat(k, k)

    DO l = k + 1, n
      amat(l, k) = amat(l, k) * amat(k, k)
    END DO

    DO m = k + 1, n
      DO l = k+1, n
        amat(l, m) = amat(l, m) - amat(k, m) * amat(l, k)
      END DO
    END DO
  END DO  ! loop k = 1, n

  DO m = 1, n - 1
    DO l = m + 1, n
      bvec(l) = bvec(l) - bvec(m) * amat(l, m)
    END DO
  END DO

  DO m = n, 1, -1
    bvec(m) = bvec(m) * amat(m, m)
    DO l = 1, m - 1
      bvec(l) = bvec(l) - amat(l, m) * bvec(m)
    END DO
  END DO

  DO k = 1, tke_levels
    tsq_k(k) = bvec(k)
    qsq_k(k) = bvec(tke_levels + k)
    cov_k(k) = bvec(2 * tke_levels + k)
  END DO
  IF (lhook) CALL dr_hook('MYM_SOLVE_SIMEQ_LUD',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE mym_solve_simeq_lud
