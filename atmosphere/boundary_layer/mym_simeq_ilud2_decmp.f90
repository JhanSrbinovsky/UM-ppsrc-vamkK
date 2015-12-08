! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_simeq_ilud2_decmp-----------------------------------
!
!  Purpose: To perform the incomplete LU decomposition with fill-in
!           level 2
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE mym_simeq_ilud2_decmp(                                       &
      aa_tsq_k, bb_tsq_k, cc_tsq_k, pp_tc_k,                            &
      aa_qsq_k, bb_qsq_k, cc_qsq_k, pp_qc_k,                            &
      aa_cov_k, bb_cov_k, cc_cov_k, pp_ct_k, pp_cq_k,                   &
      aap_tsq_k, r_bbp_tsq_k, ccp_tsq_k,                                &
      ppp_tc_k, pp1_tc_k,  pp2_tc_k,                                    &
      aap_qsq_k, r_bbp_qsq_k, ccp_qsq_k,                                &
      ppp_qc_k, pp1_qc_k,  pp2_qc_k,                                    &
      aap_cov_k, r_bbp_cov_k, ccp_cov_k,                                &
      ppp_ct_k, ppp_cq_k, pp1_ct_k, pp1_cq_k, pp2_ct_k, pp2_cq_k)

  USE mym_option_mod, ONLY: tke_levels
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

! intent in variables
  REAL, INTENT(IN) ::                                                   &
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
     aap_tsq_k(tke_levels),                                             &
     r_bbp_tsq_k(tke_levels),                                           &
     ccp_tsq_k(tke_levels),                                             &
     ppp_tc_k(tke_levels),                                              &
     pp1_tc_k(tke_levels),                                              &
     pp2_tc_k(tke_levels),                                              &
     aap_qsq_k(tke_levels),                                             &
     r_bbp_qsq_k(tke_levels),                                           &
     ccp_qsq_k(tke_levels),                                             &
     ppp_qc_k(tke_levels),                                              &
     pp1_qc_k(tke_levels),                                              &
     pp2_qc_k(tke_levels),                                              &
     aap_cov_k(tke_levels),                                             &
     r_bbp_cov_k(tke_levels),                                           &
     ccp_cov_k(tke_levels),                                             &
     ppp_ct_k(tke_levels),                                              &
     ppp_cq_k(tke_levels),                                              &
     pp1_ct_k(tke_levels),                                              &
     pp1_cq_k(tke_levels),                                              &
     pp2_ct_k(tke_levels),                                              &
     pp2_cq_k(tke_levels)
             ! matrix elements of the ILU decomposed matrix

  INTEGER :: k
             ! loop indexes

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_SIMEQ_ILUD2_DECMP',zhook_in,zhook_handle)

  aap_tsq_k(1) = aa_tsq_k(1)
  r_bbp_tsq_k(1) = 1.0 / bb_tsq_k(1)
  ccp_tsq_k(1) = cc_tsq_k(1)

  aap_qsq_k(1) = aa_qsq_k(1)
  r_bbp_qsq_k(1) = 1.0 / bb_qsq_k(1)
  ccp_qsq_k(1) = cc_qsq_k(1)

  ppp_tc_k(1) = pp_tc_k(1)
  ppp_qc_k(1) = pp_qc_k(1)

  pp1_tc_k(1) = 0.0
  pp1_qc_k(1) = 0.0
  pp2_tc_k(1) = 0.0
  pp2_qc_k(1) = 0.0

  DO k = 2, tke_levels
    aap_tsq_k(k) = aa_tsq_k(k)
    r_bbp_tsq_k(k) = 1.0 / (bb_tsq_k(k)                                 &
              - aap_tsq_k(k) * ccp_tsq_k(k - 1) * r_bbp_tsq_k(k - 1))
    ccp_tsq_k(k) = cc_tsq_k(k)

    aap_qsq_k(k) = aa_qsq_k(k)
    r_bbp_qsq_k(k) = 1.0 / (bb_qsq_k(k)                                 &
              - aap_qsq_k(k) * ccp_qsq_k(k - 1) * r_bbp_qsq_k(k - 1))
    ccp_qsq_k(k) = cc_qsq_k(k)

    ppp_tc_k(k) = pp_tc_k(k)
    ppp_qc_k(k) = pp_qc_k(k)

    pp1_tc_k(k) = - aap_tsq_k(k) * r_bbp_tsq_k(k - 1) * ppp_tc_k(k - 1)
    pp1_qc_k(k) = - aap_qsq_k(k) * r_bbp_qsq_k(k - 1) * ppp_qc_k(k - 1)

    pp2_tc_k(k) = - aap_tsq_k(k) * r_bbp_tsq_k(k - 1) * pp1_tc_k(k - 1)
    pp2_qc_k(k) = - aap_qsq_k(k) * r_bbp_qsq_k(k - 1) * pp1_qc_k(k - 1)

  END DO

  k = 1

  ppp_ct_k(k) = pp_ct_k(k)
  ppp_cq_k(k) = pp_cq_k(k)

  pp1_ct_k(k) = -ppp_ct_k(k) * r_bbp_tsq_k(k) * ccp_tsq_k(k)
  pp1_cq_k(k) = -ppp_cq_k(k) * r_bbp_qsq_k(k) * ccp_qsq_k(k)

  pp2_ct_k(k) = - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * ccp_tsq_k(k + 1)
  pp2_cq_k(k) = - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * ccp_qsq_k(k + 1)

  aap_cov_k(k) = 0.0

  r_bbp_cov_k(k) = 1.0 / (                                              &
                  bb_cov_k(k)                                           &
                  - ppp_ct_k(k) * r_bbp_tsq_k(k) * ppp_tc_k(k)          &
                  - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * pp1_tc_k(k + 1)  &
                  - pp2_ct_k(k) * r_bbp_tsq_k(k + 2) * pp2_tc_k(k + 2)  &
                  - ppp_cq_k(k) * r_bbp_qsq_k(k) * ppp_qc_k(k)          &
                  - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * pp1_qc_k(k + 1)  &
                  - pp2_cq_k(k) * r_bbp_qsq_k(k + 2) * pp2_qc_k(k + 2))

  ccp_cov_k(k) = cc_cov_k(k)                                            &
                 - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * ppp_tc_k(k + 1)   &
                 - pp2_ct_k(k) * r_bbp_tsq_k(k + 2) * pp1_tc_k(k + 2)   &
                 - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * ppp_qc_k(k + 1)   &
                 - pp2_cq_k(k) * r_bbp_qsq_k(k + 2) * pp1_qc_k(k + 2)

  DO k = 2, tke_levels - 2
    ppp_ct_k(k) = pp_ct_k(k)
    ppp_cq_k(k) = pp_cq_k(k)

    pp1_ct_k(k) = -ppp_ct_k(k) * r_bbp_tsq_k(k) * ccp_tsq_k(k)
    pp1_cq_k(k) = -ppp_cq_k(k) * r_bbp_qsq_k(k) * ccp_qsq_k(k)

    pp2_ct_k(k) = - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * ccp_tsq_k(k + 1)
    pp2_cq_k(k) = - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * ccp_qsq_k(k + 1)

    aap_cov_k(k) = aa_cov_k(k)                                          &
                   - ppp_ct_k(k) * r_bbp_tsq_k(k) * pp1_tc_k(k)         &
                   - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * pp2_tc_k(k + 1) &
                   - ppp_cq_k(k) * r_bbp_qsq_k(k) * pp1_qc_k(k)         &
                   - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * pp2_qc_k(k + 1)

    r_bbp_cov_k(k) = 1.0 / (                                            &
                 bb_cov_k(k)                                            &
                 - ppp_ct_k(k) * r_bbp_tsq_k(k) * ppp_tc_k(k)           &
                 - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * pp1_tc_k(k + 1)   &
                 - pp2_ct_k(k) * r_bbp_tsq_k(k + 2) * pp2_tc_k(k + 2)   &
                 - ppp_cq_k(k) * r_bbp_qsq_k(k) * ppp_qc_k(k)           &
                 - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * pp1_qc_k(k + 1)   &
                 - pp2_cq_k(k) * r_bbp_qsq_k(k + 2) * pp2_qc_k(k + 2)   &
                 - aap_cov_k(k) * r_bbp_cov_k(k - 1) * ccp_cov_k(k - 1))

    ccp_cov_k(k) = cc_cov_k(k)                                          &
                   - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * ppp_tc_k(k + 1) &
                   - pp2_ct_k(k) * r_bbp_tsq_k(k + 2) * pp1_tc_k(k + 2) &
                   - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * ppp_qc_k(k + 1) &
                   - pp2_cq_k(k) * r_bbp_qsq_k(k + 2) * pp1_qc_k(k + 2)
  END DO

  k = tke_levels - 1

  ppp_ct_k(k) = pp_ct_k(k)
  ppp_cq_k(k) = pp_cq_k(k)

  pp1_ct_k(k) = -ppp_ct_k(k) * r_bbp_tsq_k(k) * ccp_tsq_k(k)
  pp1_cq_k(k) = -ppp_cq_k(k) * r_bbp_qsq_k(k) * ccp_qsq_k(k)
  pp2_ct_k(k) = 0.0
  pp2_cq_k(k) = 0.0

  aap_cov_k(k) = aa_cov_k(k)                                            &
                 - ppp_ct_k(k) * r_bbp_tsq_k(k) * pp1_tc_k(k)           &
                 - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * pp2_tc_k(k + 1)   &
                 - ppp_cq_k(k) * r_bbp_qsq_k(k) * pp1_qc_k(k)           &
                 - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * pp2_qc_k(k + 1)

  r_bbp_cov_k(k) = 1.0 / (                                              &
                 bb_cov_k(k)                                            &
                 - ppp_ct_k(k) * r_bbp_tsq_k(k) * ppp_tc_k(k)           &
                 - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * pp1_tc_k(k + 1)   &
                 - ppp_cq_k(k) * r_bbp_qsq_k(k) * ppp_qc_k(k)           &
                 - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * pp1_qc_k(k + 1)   &
                 - aap_cov_k(k) * r_bbp_cov_k(k - 1) * ccp_cov_k(k - 1))

  ccp_cov_k(k) = cc_cov_k(k)                                            &
                 - pp1_ct_k(k) * r_bbp_tsq_k(k + 1) * ppp_tc_k(k + 1)   &
                 - pp1_cq_k(k) * r_bbp_qsq_k(k + 1) * ppp_qc_k(k + 1)

  k = tke_levels

  ppp_ct_k(k) = pp_ct_k(k)
  ppp_cq_k(k) = pp_cq_k(k)

  pp1_ct_k(k) = 0.0
  pp1_cq_k(k) = 0.0
  pp2_ct_k(k) = 0.0
  pp2_cq_k(k) = 0.0

  aap_cov_k(k) = aa_cov_k(k)                                            &
                 - ppp_ct_k(k) * r_bbp_tsq_k(k) * pp1_tc_k(k)           &
                 - ppp_cq_k(k) * r_bbp_qsq_k(k) * pp1_qc_k(k)
  r_bbp_cov_k(k) = 1.0 / (                                              &
                 bb_cov_k(k)                                            &
                 - ppp_ct_k(k) * r_bbp_tsq_k(k) * ppp_tc_k(k)           &
                 - ppp_cq_k(k) * r_bbp_qsq_k(k) * ppp_qc_k(k)           &
                 - aap_cov_k(k) * r_bbp_cov_k(k - 1) * ccp_cov_k(k - 1))
  ccp_cov_k(k) = 0.0

  IF (lhook) CALL dr_hook('MYM_SIMEQ_ILUD2_DECMP',                      &
                                                zhook_out,zhook_handle)
  RETURN

END SUBROUTINE mym_simeq_ilud2_decmp
