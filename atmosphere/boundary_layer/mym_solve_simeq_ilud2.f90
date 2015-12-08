! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_solve_simeq_ilud2-----------------------------------
!
!  Purpose: To solve simultaneous equations of which the coefficient
!           matrix is obtained by imcompelete LU decomposition
!           with fill-in level 2 (ILU(2)) for the original coefficient
!           matrix.
!           ILU(2) decomposition is assumed to have been already done
!           in mym_simeq_ilud2_dcmp.
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE mym_solve_simeq_ilud2(                                       &
      imode,                                                            &
      qq_tsq_k, qq_qsq_k, qq_cov_k,                                     &
      aap_tsq_k, r_bbp_tsq_k, ccp_tsq_k,                                &
      ppp_tc_k, pp1_tc_k,  pp2_tc_k,                                    &
      aap_qsq_k, r_bbp_qsq_k, ccp_qsq_k,                                &
      ppp_qc_k, pp1_qc_k,  pp2_qc_k,                                    &
      aap_cov_k, r_bbp_cov_k, ccp_cov_k,                                &
      ppp_ct_k, ppp_cq_k, pp1_ct_k, pp1_cq_k, pp2_ct_k, pp2_cq_k,       &
      tsq_k, qsq_k, cov_k)

  USE mym_option_mod, ONLY: tke_levels
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

! intent in variables
  INTEGER, INTENT(IN) :: imode
                         ! mode switch for the Matrix
                         ! 0: normal, 1: transposed

  REAL, INTENT(IN) ::                                                   &
     qq_tsq_k(tke_levels),                                              &
     qq_qsq_k(tke_levels),                                              &
     qq_cov_k(tke_levels),                                              &
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
             ! matrix elements of ILU decomposed matrix
             ! See the document for details

  REAL, INTENT(OUT) ::                                                  &
     tsq_k(tke_levels),                                                 &
     qsq_k(tke_levels),                                                 &
     cov_k(tke_levels)
             ! solution vectors

  INTEGER :: k
             ! loop indexes

  INTEGER, PARAMETER ::                                                 &
     normal = 0,                                                        &
     transposed = 1
             ! symbols for the mode

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_SOLVE_SIMEQ_ILUD2',zhook_in,zhook_handle)

  tsq_k(1) = qq_tsq_k(1) * r_bbp_tsq_k(1)
  qsq_k(1) = qq_qsq_k(1) * r_bbp_qsq_k(1)

  IF (imode == normal) THEN
    DO k = 2, tke_levels
      tsq_k(k) = (qq_tsq_k(k) - aap_tsq_k(k) * tsq_k(k - 1))            &
                                                 * r_bbp_tsq_k(k)
      qsq_k(k) = (qq_qsq_k(k) - aap_qsq_k(k) * qsq_k(k - 1))            &
                                                 * r_bbp_qsq_k(k)
    END DO

    k = 1
    cov_k(k) = (qq_cov_k(k) - ppp_ct_k(k) * tsq_k(k)                    &
                            - pp1_ct_k(k) * tsq_k(k + 1)                &
                            - pp2_ct_k(k) * tsq_k(k + 2)                &
                            - ppp_cq_k(k) * qsq_k(k)                    &
                            - pp1_cq_k(k) * qsq_k(k + 1)                &
                            - pp2_cq_k(k) * qsq_k(k + 2))               &
               * r_bbp_cov_k(k)

    DO k = 2, tke_levels - 2
      cov_k(k) = (qq_cov_k(k) - ppp_ct_k(k) * tsq_k(k)                  &
                              - pp1_ct_k(k) * tsq_k(k + 1)              &
                              - pp2_ct_k(k) * tsq_k(k + 2)              &
                              - ppp_cq_k(k) * qsq_k(k)                  &
                              - pp1_cq_k(k) * qsq_k(k + 1)              &
                              - pp2_cq_k(k) * qsq_k(k + 2)              &
                              - aap_cov_k(k) * cov_k(k - 1))            &
               * r_bbp_cov_k(k)
    END DO

    k = tke_levels - 1
    cov_k(k) = (qq_cov_k(k) - ppp_ct_k(k) * tsq_k(k)                    &
                            - pp1_ct_k(k) * tsq_k(k + 1)                &
                            - ppp_cq_k(k) * qsq_k(k)                    &
                            - pp1_cq_k(k) * qsq_k(k + 1)                &
                            - aap_cov_k(k) * cov_k(k - 1))              &
               * r_bbp_cov_k(k)


    k = tke_levels
    cov_k(k) = (qq_cov_k(k) - ppp_ct_k(k) * tsq_k(k)                    &
                            - ppp_cq_k(k) * qsq_k(k)                    &
                            - aap_cov_k(k) * cov_k(k - 1))              &
               * r_bbp_cov_k(k)


    DO k = tke_levels - 1, 1, -1
      cov_k(k) = cov_k(k)                                               &
                        - ccp_cov_k(k) * cov_k(k + 1) * r_bbp_cov_k(k)
    END DO

    k = tke_levels
    qsq_k(k) = qsq_k(k) - (ppp_qc_k(k) * cov_k(k)                       &
                          + pp1_qc_k(k) * cov_k(k - 1)                  &
                          + pp2_qc_k(k) * cov_k(k - 2))                 &
                         * r_bbp_qsq_k(k)
    tsq_k(k) = tsq_k(k) - (ppp_tc_k(k) * cov_k(k)                       &
                          + pp1_tc_k(k) * cov_k(k - 1)                  &
                          + pp2_tc_k(k) * cov_k(k - 2))                 &
                         * r_bbp_tsq_k(k)

    DO k = tke_levels - 1, 3, -1
      qsq_k(k) = qsq_k(k)                                               &
            - (ppp_qc_k(k) * cov_k(k) + ccp_qsq_k(k) * qsq_k(k + 1)     &
               + pp1_qc_k(k) * cov_k(k - 1)                             &
               + pp2_qc_k(k) * cov_k(k - 2))                            &
            * r_bbp_qsq_k(k)
      tsq_k(k) = tsq_k(k)                                               &
            - (ppp_tc_k(k) * cov_k(k) + ccp_tsq_k(k) * tsq_k(k + 1)     &
               + pp1_tc_k(k) * cov_k(k - 1)                             &
               + pp2_tc_k(k) * cov_k(k - 2))                            &
            * r_bbp_tsq_k(k)
    END DO

    k = 2
    qsq_k(k) = qsq_k(k)                                                 &
          - (ppp_qc_k(k) * cov_k(k) + ccp_qsq_k(k) * qsq_k(k + 1)       &
             + pp1_qc_k(k) * cov_k(k - 1))                              &
          * r_bbp_qsq_k(k)
    tsq_k(k) = tsq_k(k)                                                 &
          - (ppp_tc_k(k) * cov_k(k) + ccp_tsq_k(k) * tsq_k(k + 1)       &
             + pp1_tc_k(k) * cov_k(k - 1))                              &
          * r_bbp_tsq_k(k)


    k = 1
    qsq_k(k) = qsq_k(k)                                                 &
          - (ppp_qc_k(k) * cov_k(k) + ccp_qsq_k(k) * qsq_k(k + 1))      &
          * r_bbp_qsq_k(k)
    tsq_k(k) = tsq_k(k)                                                 &
          - (ppp_tc_k(k) * cov_k(k) + ccp_tsq_k(k) * tsq_k(k + 1))      &
          * r_bbp_tsq_k(k)

  ELSE IF (imode == transposed) THEN
    DO k = 2, tke_levels
      tsq_k(k) = (qq_tsq_k(k)                                           &
                   - ccp_tsq_k(k - 1) * tsq_k(k - 1)) * r_bbp_tsq_k(k)
      qsq_k(k) = (qq_qsq_k(k)                                           &
                   - ccp_qsq_k(k - 1) * qsq_k(k - 1)) * r_bbp_qsq_k(k)
    END DO

    k = 1
    cov_k(k) = (qq_cov_k(k) - ppp_tc_k(k) * tsq_k(k)                    &
                            - pp1_tc_k(k + 1) * tsq_k(k + 1)            &
                            - pp2_tc_k(k + 2) * tsq_k(k + 2)            &
                            - ppp_qc_k(k) * qsq_k(k)                    &
                            - pp1_qc_k(k + 1) * qsq_k(k + 1)            &
                            - pp2_qc_k(k + 2) * qsq_k(k + 2))           &
              * r_bbp_cov_k(k)

    DO k = 2, tke_levels - 2
      cov_k(k) = (qq_cov_k(k) - ppp_tc_k(k) * tsq_k(k)                  &
                              - pp1_tc_k(k + 1) * tsq_k(k + 1)          &
                              - pp2_tc_k(k + 2) * tsq_k(k + 2)          &
                              - ppp_qc_k(k) * qsq_k(k)                  &
                              - pp1_qc_k(k + 1) * qsq_k(k + 1)          &
                              - pp2_qc_k(k + 2) * qsq_k(k + 2)          &
                              - ccp_cov_k(k - 1) * cov_k(k - 1))        &
                * r_bbp_cov_k(k)
    END DO

    k = tke_levels - 1
    cov_k(k) = (qq_cov_k(k) - ppp_tc_k(k) * tsq_k(k)                    &
                            - pp1_tc_k(k + 1) * tsq_k(k + 1)            &
                            - ppp_qc_k(k) * qsq_k(k)                    &
                            - pp1_qc_k(k + 1) * qsq_k(k + 1)            &
                            - ccp_cov_k(k - 1) * cov_k(k - 1))          &
              * r_bbp_cov_k(k)


    k = tke_levels
    cov_k(k) = (qq_cov_k(k) - ppp_tc_k(k) * tsq_k(k)                    &
                          - ppp_qc_k(k) * qsq_k(k)                      &
                          - ccp_cov_k(k - 1) * cov_k(k - 1))            &
              * r_bbp_cov_k(k)

    DO k = tke_levels - 1, 1, -1
      cov_k(k) = cov_k(k)                                               &
                   - aap_cov_k(k + 1) * cov_k(k + 1) * r_bbp_cov_k(k)
    END DO

    k = tke_levels
    qsq_k(k) = qsq_k(k) - (ppp_cq_k(k) * cov_k(k)                       &
                          + pp1_cq_k(k - 1) * cov_k(k - 1)              &
                          + pp2_cq_k(k - 2) * cov_k(k - 2))             &
                         * r_bbp_qsq_k(k)
    tsq_k(k) = tsq_k(k) - (ppp_ct_k(k) * cov_k(k)                       &
                          + pp1_ct_k(k - 1) * cov_k(k - 1)              &
                          + pp2_ct_k(k - 2) * cov_k(k - 2))             &
                         * r_bbp_tsq_k(k)

    DO k = tke_levels - 1, 3, -1
      qsq_k(k) = qsq_k(k)                                               &
            - (ppp_cq_k(k) * cov_k(k) + aap_qsq_k(k + 1) * qsq_k(k + 1) &
               + pp1_cq_k(k - 1) * cov_k(k - 1)                         &
               + pp2_cq_k(k - 2) * cov_k(k - 2))                        &
            * r_bbp_qsq_k(k)
      tsq_k(k) = tsq_k(k)                                               &
            - (ppp_ct_k(k) * cov_k(k) + aap_tsq_k(k + 1) * tsq_k(k + 1) &
               + pp1_ct_k(k - 1) * cov_k(k - 1)                         &
               + pp2_ct_k(k - 2) * cov_k(k - 2))                        &
            * r_bbp_tsq_k(k)
    END DO

    k = 2
    qsq_k(k) = qsq_k(k)                                                 &
          - (ppp_cq_k(k) * cov_k(k) + aap_qsq_k(k + 1) * qsq_k(k + 1)   &
             + pp1_cq_k(k - 1) * cov_k(k - 1))                          &
          * r_bbp_qsq_k(k)
    tsq_k(k) = tsq_k(k)                                                 &
          - (ppp_ct_k(k) * cov_k(k) + aap_tsq_k(k + 1) * tsq_k(k + 1)   &
             + pp1_ct_k(k - 1) * cov_k(k - 1))                          &
          * r_bbp_tsq_k(k)


    k = 1
    qsq_k(k) = qsq_k(k)                                                 &
          - (ppp_cq_k(k) * cov_k(k) + aap_qsq_k(k + 1) * qsq_k(k + 1))  &
          * r_bbp_qsq_k(k)
    tsq_k(k) = tsq_k(k)                                                 &
          - (ppp_ct_k(k) * cov_k(k) + aap_tsq_k(k + 1) * tsq_k(k + 1))  &
          * r_bbp_tsq_k(k)

  END IF

  IF (lhook) CALL dr_hook('MYM_SOLVE_SIMEQ_ILUD2',                      &
                                                zhook_out,zhook_handle)
  RETURN

END SUBROUTINE mym_solve_simeq_ilud2
