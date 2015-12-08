! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_turbulence-----------------------------------------
!
!  Purpose: To calculate the diffusion coefficients and counter
!           gradient term for momentum, heat and moisture, and
!           integrate the prognostic variables appearing
!           in the MY model.
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
!
!  This code is based on the code provided by the authors who wrote
!  the following papers.
!   * Nakanishi, M. and H. Niino, 2009: Development of an improved
!      turbulence closure model for the atmospheric boundary layer.
!      J. Meteor. Soc. Japan, 87, 895-912.
!   * Nakanishi, M. and H. Niino, 2006: An improved Mellor-Yamada
!      Level-3 model: Its numerical stability and application to
!      a regional prediction of advection fog.
!      Boundary-Layer Meteor., 119, 397-407.
!   * Nakanishi, M. and H. Niino, 2004: An improved Mellor-Yamada
!      Level-3 model with condensation physics: Its design and
!      verification.
!      Boundary-Layer Meteor., 112, 1-31.
!   * Nakanishi, M., 2001: Improvement of the Mellor-Yamada
!      turbulence closure model based on large-eddy simulation data.
!      Boundary-Layer Meteor., 99, 349-378.
!   The web site publicising their code:
!    http://www.nda.ac.jp/~naka/MYNN/index.html
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE mym_turbulence(                                              &
! IN levels/switches
      bl_levels, levflag, nSCMDpkgs,L_SCMDiags, BL_diag,                &
! IN fields
      t, q, qcl, qcf, z_uv, z_tq,                                       &
      vq, vt, gtr, fqw, ftl, wb_ng,                                     &
      dbdz, dtldz, dqwdz, dvdzm, dudz, dvdz,                            &
      r_mosurf, u_s, fb_surf, pmz, phh,                                 &
! INOUT fields
      qke, tsq, qsq, cov, dfm, dfh,                                     &
! OUT fields
      dfu_cg, dfv_cg, dft_cg, dfq_cg)

  USE atm_fields_bounds_mod, ONLY: tdims, pdims, tdims_l, tdims_s
  USE atmos_constants_mod, ONLY: vkman
  USE conversions_mod, ONLY: pi
  USE mym_const_mod, ONLY: e1c,e2c,e3c,e4c,e5c,a1,a2,c1,b2,qke_max,     &
        coef_trbvar_diff,coef_trbvar_diff_tke,two_thirds,a1_2,          &
        b1,one_third,cc3
  USE mym_option_mod, ONLY:                                             &
        my_lowest_pd_surf, l_my_extra_level, my_z_extra_fact,           &
        l_my_prod_adj, my_prod_adj_fact, tke_levels,                    &
        l_my_lowest_pd_surf_tqc
  USE bl_diags_mod, ONLY: strnewbldiag
  USE earth_constants_mod, ONLY: g
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

! Intent IN Variables
  INTEGER, INTENT(IN) ::                                                &
     bl_levels,                                                         &
                   ! Max. no. of "boundary" levels
     levflag
                   ! to indicate the level of the MY model
                   ! 2: level 2.5
                   ! 3: level 3

  REAL, INTENT(IN) ::                                                   &
     t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                   ! temperature
     q(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                   ! specific humidity
     qcl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                   ! Cloud liquid water
     qcf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                   ! Cloud ice (kg per kg air)
     z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels),                                                   &
                   ! Z_UV(*,K) is height of u level k
     z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                   ! Z_TQ(*,K) is height of theta level k.
     vq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! A buoyancy param on theta level k-1
     vt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! A buoyancy param on theta level k-1
     gtr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                   ! g/thetav on theta level k-1
     fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                   ! Moisture flux between layers
                   !   (kg per square metre per sec).
                   ! FQW(,1) is total water flux
                   ! from surface, 'E'.
                   ! Defined on rho levels.
     ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                   ! FTL(,K) contains net turbulent
                   ! sensible heat flux into layer K
                   ! from below; so FTL(,1) is the
                   ! surface sensible heat, H. (W/m2)
                   ! Defined on rho levels.
     wb_ng(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           tke_levels),                                                 &
                   ! buoyancy flux related to the skewness
                   ! on theta K-1 levels
     dbdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          2:tke_levels),                                                &
                   ! Buoyancy gradient across layer
                   ! interface interpolated to theta levels.
                   ! (:,:,K) repserents the value on theta level K-1
     dtldz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           2:bl_levels),                                                &
                   ! gradient of TL across layer
                   ! interface interpolated to theta levels.
                   ! (:,:,K) repserents the value on theta level K-1
     dqwdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           2:bl_levels),                                                &
                   ! gradient of QW across layer
                   ! interface interpolated to theta levels.
                   ! (:,:,K) repserents the value on theta level K-1
     dvdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           2:bl_levels),                                                &
                   ! Modulus of wind shear at theta levels.
                   ! (:,:,K) repserents the value on theta level K-1
     dudz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          2:bl_levels),                                                 &
                   ! Gradient of u at theta levels.
                   !(:,:,K) repserents the value on theta level K-1
     dvdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          2:bl_levels),                                                 &
                   ! Gradient of v at theta levels.
                   !(:,:,K) repserents the value on theta level K-1
     r_mosurf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
                   ! reciprocal of Monin-Obukhov length
     u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
                   ! Surface friction velocity
     fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),      &
                   ! Surface buoyancy flux
     pmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
                   ! gradient function for momentum at surface
     phh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                   ! gradient function for scalars at surface

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER, INTENT(IN) ::                                                &
     nSCMDpkgs             ! No of SCM diagnostics packages

  LOGICAL, INTENT(IN) ::                                                &
     L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

! Intent INOUT Variables
  REAL, INTENT(INOUT) ::                                                &
     qke(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
         bl_levels),                                                    &
                   ! twice of TKE (denoted to q**2) on theta level K-1
     tsq(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
         bl_levels),                                                    &
                   ! Self covariance of liquid potential temperature
                   ! (thetal'**2) defined on theta levels K-1
     qsq(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
         bl_levels),                                                    &
                   ! Self covariance of total water
                   ! (qw'**2) defined on theta levels K-1
     cov(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
         bl_levels),                                                    &
                   ! Correlation between thetal and qw
                   ! (thetal'qw') defined on theta levels K-1
     dfm(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,   &
         bl_levels),                                                    &
                   ! diffusion coefficient for momentum
                   ! on theta level K-1
     dfh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels)
                   ! diffusion coefficient for scalars
                   ! on theta level K-1

!  Declaration of BL diagnostics.
  TYPE (strnewbldiag), INTENT(INOUT) :: BL_diag

! Intent OUT Variables
  REAL, INTENT(OUT) ::                                                  &
     dfu_cg(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,&
            bl_levels),                                                 &
                   ! counter gradient term for u
                   ! on theta level K-1
     dfv_cg(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,&
            bl_levels),                                                 &
                   ! counter gradient term for v
                   ! on theta level K-1
     dft_cg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
                   ! counter gradient term for TL
                   ! on theta level K-1
     dfq_cg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels)
                   ! counter gradient term for QW
                   ! on theta level K-1

! Local variables
! Scalar
  INTEGER ::                                                            &
     i, j, k, k_start, k_start_cor
                   ! Loop indexes






  REAL ::                                                               &
     e1,                                                                &
                   ! a variable denoted to E1 in the papers
     e3,                                                                &
                   ! a variable denoted to E3 in the papers
     e4,                                                                &
                   ! a variable denoted to E4 in the paper
     q2sq,                                                              &
                   ! qke derived by level 2 scheme
     t2sq,                                                              &
                   ! tsq derived by level 2 scheme
     r2sq,                                                              &
                   ! qsq derived by level 2 scheme
     t3sq,                                                              &
                   ! tsq derived by level 2.5 or 3 scheme
     r3sq,                                                              &
                   ! qsq derived by level 2.5 or 3 scheme
     c3sq,                                                              &
                   ! cov deribed by level 2.5 or 3 scheme
     wden,                                                              &
                   ! work variable
     eden,                                                              &
                   ! work variable
     reden,                                                             &
                   ! reciprocal of eden
     cw,                                                                &
                   ! <w'**2>/qke = 1-(<u'**2> + <v'**2>) / qke
     e6c,                                                               &
                   ! work variable
     coef,                                                              &
                   ! work variable
     elq,                                                               &
                   ! mixing length  times qkw appeared
                   ! in the production term of qke
     elh,                                                               &
                   ! mixing length times qkw appeared
                   ! in the production terms of tsq, qsq and cov.
     phm,                                                               &
                   ! work variable
     enum,                                                              &
                   ! work variable
     disp_coef,                                                         &
     b1l,                                                               &
                   ! B1 (closure constant) times mixing length
     b2l,                                                               &
                   ! work variable
     clow,                                                              &
                   ! lower limit for difference between cov in level 3
                   ! and level 2
     cupp
                   ! upper limit for difference between cov in level 3
                   ! and level 2


  REAL ::                                                               &
     gm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! square of wind shear on theta level K-1
                   ! (a denominator of gradient Richardson number)
     gh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! - buoyancy gradient on theta level K-1
                   ! (a numerator of gradient Richardson number)
     sm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! Non-dimensional diffusion coefficients for
                   ! momentum derived by level 2 scheme
                   ! defined on theta level K-1
     sh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! Non-dimensional diffusion coefficients for
                   ! scalars derived by level 2 scheme
                   ! define on theta level K-1
     qkw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
          tke_levels),                                                  &
                   ! q=sqrt(qke) on theta level K-1
     elsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                   ! square of mixing length
     gmel(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                   ! GM times the mixing length
     ghel(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                   ! GH times the mixing length
     qdiv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                   ! factor for flux correction: sqrt(q3sq/q2sq)
     gamv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                   ! counter gradient correction for production of qke
     gamv_coef(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               tke_levels),                                             &
                   ! coefficient of (c3sq-c2sq) in gamv
     gamt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                   ! counter gradient term for flux of TL
                   ! gamt = gamt_tsq * tsq + gamt_cov * cov + gamt_res
     gamt_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              tke_levels),                                              &
                   ! a linear part to tsq in gamt
     gamt_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              tke_levels),                                              &
                   ! a linear part to cov in gamt
     gamt_res(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              tke_levels),                                              &
                   ! a residual part in gamt
     gamt_factor(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,   &
                 tke_levels),                                           &
                   ! stability factor for gamt
     gamq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                   ! counter gradient term for flux of QW
                   ! gamq = gamq_qsq * qsq + gamq_cov * cov + gamq_res
     gamq_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              tke_levels),                                              &
                   ! a linear part to qsq in gamq
     gamq_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              tke_levels),                                              &
                   ! a linear part to cov in gamq
     gamq_res(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              tke_levels),                                              &
                   ! a residual part in gamq
     gamq_factor(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,   &
                 tke_levels),                                           &
                   ! stability factor for gamq
     pdc_factor(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                tke_levels),                                            &
                   ! stability factor for pdc
     smd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                   ! counter gradient correction for SM
     smd_coef(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              tke_levels),                                              &
                   ! coefficient of (c3sq-c2sq) in smd
     e2(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        tke_levels),                                                    &
                   ! a variable denoted to E2 in the papers
     cu(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! <u'**2> / qke in level 3
     cv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! <v'**2> / qke  in level 3
     cu25(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                   ! <u'**2> / qke  in level 2.5
     cv25(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                   ! <v'**2> / qke in level 2.5
     cw25(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                   ! <w'**2> / qke = 1 - <u'**2> - <v'**2> in level 2.5
     pdk(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
          tke_levels),                                                  &
                   ! production term of qke
     pdt_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to tsq in a production term of tsq
     pdt_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to cov in a production term of tsq
     pdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                   ! production term of tsq excluding a linear part
                   ! of tsq and cov
                   ! (Note that it is the production term itself
                   ! in level 2.5)
     pdq_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to qsq in a production term of qsq
     pdq_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to cov in a production term of qsq
     pdq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                   ! production term of qsq excluding a linear part
                   ! of qsq and cov
                   ! (Note that it is the production term itself
                   ! in level 2.5)
     pdc_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to tsq in a production term of cov
     pdc_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to qsq in a production term of cov
     pdc_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                   ! a linear part to cov in a production term of cov
     pdc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                   ! production term of cov excluding a linear part
                   ! of cov, tsq and qsq
                   ! (Note that it is the production term itself
                   ! in level 2.5)
     bp(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! coefficients of qke in the dissipation term
     rp(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! production term of qke
     el(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! mixing length
     q3sq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                   ! qke derived by level 2.5 or level 3
     c2sq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
           tke_levels)
                   ! cov derived by level 2

  REAL, ALLOCATABLE ::                                                  &
     ! These variables are required only when imp_mode /= FULL_IMPL
     ! So usually they are not used.
     ! (That is why they have an "allocatable" attribute.)
     bp_tsq(:, :, :),                                                   &
                   ! coefficients of tsq in the dissipation term
     rp_tsq(:, :, :),                                                   &
                   ! production term of tsq
     bp_qsq(:, :, :),                                                   &
                   ! coefficients of qsq in the dissipation term
     rp_qsq(:, :, :),                                                   &
                   ! production term of qsq
     bp_cov(:, :, :),                                                   &
                   ! coefficients of cov in the dissipation term
     rp_cov(:, :, :)
                   ! production term of cov

  INTEGER, PARAMETER ::                                                 &
     ! Symbols for a switch
     full_impl = 0,                                                     &
     half_impl = 1,                                                     &
     expl      = 2

  INTEGER, PARAMETER ::                                                 &
     imp_mode = full_impl
        ! mode to integrate covariances


  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_TURBULENCE',zhook_in,zhook_handle)
  IF (l_my_extra_level) THEN
    k_start = 1
  ELSE
    k_start = 2
  END IF

! DEPENDS ON: mym_level2
  CALL mym_level2(                                                      &
        bl_levels,dbdz, dvdzm,gm, gh, sm, sh)

! DEPENDS ON: mym_length
  CALL mym_length(                                                      &
        tdims%i_end,tdims%j_end,tdims_l%halo_i,tdims_l%halo_j,bl_levels,&
        qke, z_uv, z_tq, dbdz, r_mosurf, fb_surf, qkw, el)

  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        elsq(i, j, k) = el(i, j, k) ** 2
        q2sq = b1 * elsq(i, j, k)                                       &
             * (sm(i, j, k) * gm(i, j, k) + sh(i, j, k) * gh(i, j, k))
        q3sq(i, j, k) = qkw(i, j, k) ** 2
        gmel(i, j, k) = gm(i, j, k) * elsq(i, j, k)
        ghel(i, j, k) = gh(i, j, k) * elsq(i, j, k)

        ! adjust SM and SH by SQRT(q3sq / q2sq)
        IF ( q3sq(i, j, k) < q2sq ) THEN
          qdiv(i, j, k) = SQRT(q3sq(i, j, k) / q2sq)
          sm(i, j, k) = sm(i, j, k) * qdiv(i, j, k)
          sh(i, j, k) = sh(i, j, k) * qdiv(i, j, k)

          e1   = q3sq(i, j, k)                                          &
               - e1c * ghel(i, j, k) * qdiv(i, j, k) ** 2
          e2(i, j, k)   = q3sq(i, j, k)                                 &
               - e2c * ghel(i, j, k) * qdiv(i, j, k) ** 2
          e3 = e1 + e3c * ghel(i, j, k) * qdiv(i, j, k) ** 2
          e4 = e1 - e4c * ghel(i, j, k) * qdiv(i, j, k) ** 2
          eden = e2(i, j, k) * e4                                       &
               + e3 * e5c * gmel(i, j, k) * qdiv(i, j, k) ** 2
          eden = MAX(eden, 1.0e-20)
          reden = 1.0 / eden
        ELSE
          e1 = q3sq(i, j, k) - e1c * ghel(i, j, k)
          e2(i, j, k) = q3sq(i, j, k) - e2c * ghel(i, j, k)
          e3 = e1 + e3c * ghel(i, j, k)
          e4 = e1 - e4c * ghel(i, j, k)
          eden = e2(i, j, k) * e4 + e3 * e5c * gmel(i, j, k)
          eden = MAX(eden, 1.0e-20)
          reden = 1.0 / eden

          qdiv(i, j, k) = 1.0
          sm(i, j, k) = q3sq(i, j, k) * a1 * (e3 - 3.0 * c1 *e4)        &
                             * reden
          sh(i, j, k) = q3sq(i, j, k)                                   &
               * a2 * (e2(i, j, k) + 3.0 * c1 * e5c * gmel(i, j, k))    &
               * reden
        END IF ! test if q3sq < q2sq
        cu25(i, j, k) =(e2(i, j, k)                                     &
             + 3.0 * c1 * e5c * gmel(i, j, k)                           &
             * qdiv(i, j, k) ** 2) * one_third * reden
        cv25(i, j, k) = cu25(i, j, k)                                   &
             * (e4 - 0.5 * e4c * ghel(i, j, k) * qdiv(i, j, k) ** 2)
        cw25(i, j, k) = cu25(i, j, k) * e1
        cu25(i, j, k) = 1.0 - cv25(i, j, k) - cw25(i, j, k)
      END DO
    END DO
  END DO

  IF ( levflag == 3 ) THEN
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          t2sq = qdiv(i, j, k) * b2 * elsq(i, j, k)                     &
                       * sh(i, j, k) * dtldz(i, j, k) ** 2
          r2sq = qdiv(i, j, k) * b2 * elsq(i, j, k)                     &
                       * sh(i, j, k) * dqwdz(i, j, k) ** 2
          c2sq(i, j, k) = qdiv(i, j, k) * b2 * elsq(i, j, k)            &
                * sh(i, j, k) * dtldz(i, j, k) * dqwdz(i, j, k)
          t3sq = MAX(tsq(i, j, k), 0.0)
          r3sq = MAX(qsq(i, j, k), 0.0)
          c3sq = cov(i, j, k)

          c3sq = SIGN( MIN( ABS(c3sq), SQRT(t3sq*r3sq) ), c3sq )

          t2sq = vt(i, j, k) * t2sq + vq(i, j, k) * c2sq(i, j, k)
          r2sq = vt(i, j, k) * c2sq(i, j, k) + vq(i, j, k) * r2sq
          c2sq(i, j, k) = MAX(vt(i, j, k) * t2sq + vq(i, j, k) * r2sq,  &
                              0.0)
          t3sq = vt(i, j, k) * t3sq + vq(i, j, k) * c3sq
          r3sq = vt(i, j, k) * c3sq + vq(i, j, k) * r3sq
          c3sq = MAX(vt(i, j, k) * t3sq + vq(i, j, k) * r3sq, 0.0)

          !  Limitation on q, instead of L/q
          IF ( q3sq(i, j, k) < -gh(i, j, k) * elsq(i, j, k)) THEN
            q3sq(i, j, k) = -elsq(i, j, k) * gh(i, j, k)
          END IF

          ! Limitation on c3sq (0.12 =< cw =< 0.76)
          ! e2 = q^2 * phi2'
          e2(i, j, k)   = q3sq(i, j, k)                                 &
                              - e2c*ghel(i, j, k) * qdiv(i, j, k)**2
          ! e3 = q^2 * phi3'
          e3   = q3sq(i, j, k) + e3c*ghel(i, j, k) * qdiv(i, j, k)**2
          ! e4 = q^2 * phi4'
          e4   = q3sq(i, j, k) - e4c*ghel(i, j, k) * qdiv(i, j, k)**2
          ! eden = q^4 D'
          eden = e2(i, j, k) * e4                                       &
                          + e3 *e5c*gmel(i, j, k) * qdiv(i, j, k)**2

          ! wden = numerator in the square braket in (10a) in NN2006
          !        times (1-c3) * (g/thetav)**2 * GH
          wden = cc3*gtr(i, j, k) **2                                   &
                   * elsq(i, j, k)**2 / elsq(i, j, k)                   &
                   * qdiv(i, j, k)**2                                   &
                   *( e2(i, j, k)*e4c                                   &
                         - e3c*e5c*gmel(i, j, k) * qdiv(i, j, k)**2 )

          IF ( wden /= 0.0 ) THEN
            clow = q3sq(i, j, k) * ( 0.12-cw25(i, j, k) )*eden/wden
            cupp = q3sq(i, j, k) *( 0.76-cw25(i, j, k) )*eden/wden

            IF ( wden > 0.0 ) THEN
              c3sq  = MIN( MAX( c3sq, c2sq(i, j, k) + clow),            &
                                       c2sq(i, j, k) + cupp)
            ELSE
              c3sq  = MAX( MIN( c3sq, c2sq(i, j, k) + clow),            &
                                       c2sq(i, j, k) + cupp)
            END IF
          END IF

          e1   = e2(i, j, k) + e5c*gmel(i, j, k) * qdiv(i, j, k) ** 2
          eden = MAX( eden, 1.0e-20 )
          reden = 1.0 / eden

          e6c  = 3.0 * a2 *cc3 * gtr(i, j, k)                           &
                                     * elsq(i, j, k) / elsq(i, j, k)

          ! Calculate each term in  Gamma_theta
          coef = - e1 * qdiv(i, j, k) * e6c * reden
          gamt_tsq(i, j, k) = coef * vt(i, j, k)
          gamt_cov(i, j, k) = coef * vq(i, j, k)
          gamt_res(i, j, k) = - coef * t2sq

          ! Calculate each term in  Gamma_q
          gamq_qsq(i, j, k) = coef * vq(i, j, k)
          gamq_cov(i, j, k) = coef * vt(i, j, k)
          gamq_res(i, j, k) = - coef * r2sq

          ! for Sm' and Sh'd(Theta_V)/dz
          smd_coef(i, j, k)  = elsq(i, j, k) * qdiv(i, j, k) * e6c      &
               * gtr(i, j, k) * reden * qdiv(i, j, k) ** 2              &
               * (e3c + e4c) * a1_2
          gamv_coef(i, j, k) = e1 * qdiv(i, j, k) * e6c * gtr(i, j, k)  &
                                * reden
          smd(i, j, k) = smd_coef(i, j, k) * (c3sq - c2sq(i, j, k))
          gamv(i, j, k) = gamv_coef(i, j, k) * (c3sq - c2sq(i, j,k))

          ! For elh (see below), qdiv in Level 3 is reset to 1.0.
          qdiv(i, j, k) = 1.0

          ! Calculate diffusion coefficients
          elq = el(i, j, k) * qkw(i, j, k)
          dfm(i, j, k) = elq * sm(i, j, k)
          dfh(i, j, k) = elq * sh(i, j, k)

        END DO
      END DO
    END DO

    ! Adjustment for Gamma_theta and Gamma_q
    ! After the adjustment, Gamma_theta and Gamma_q are calculated
    IF (l_my_prod_adj .AND.                                             &
          (imp_mode == half_impl .OR. imp_mode == full_impl)) THEN
      DO k = 2, tke_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

            elq = el(i, j, k) * qkw(i, j, k)
            elh = elq * qdiv(i, j, k)
            disp_coef = qkw(i, j, k) / (b2 * el(i, j, k))               &
                     + 0.5 * coef_trbvar_diff * dfm(i, j, k)            &
                        * (2.0 * pi * my_prod_adj_fact(k)               &
                          / (z_uv(i, j, k) - z_uv(i, j, k - 1))) ** 2

            pdt_tsq(i, j, k) = elh * gamt_tsq(i, j, k) * dtldz(i, j, k)
            IF (disp_coef < pdt_tsq(i, j, k)) THEN
              gamt_factor(i, j, k) = disp_coef / pdt_tsq(i, j, k)
            ELSE
              gamt_factor(i, j, k) = 1.0
            END IF

            pdq_qsq(i, j, k) = elh * gamq_qsq(i, j, k) * dqwdz(i, j, k)
            IF (disp_coef < pdq_qsq(i, j, k)) THEN
              gamq_factor(i, j, k) = disp_coef / pdq_qsq(i, j, k)
            ELSE
              gamq_factor(i, j, k) = 1.0
            END IF

            gamt_tsq(i, j, k) = gamt_factor(i, j, k) * gamt_tsq(i, j, k)
            gamt_cov(i, j, k) = gamt_factor(i, j, k) * gamt_cov(i, j, k)
            gamt_res(i, j, k) = gamt_factor(i, j, k) * gamt_res(i, j, k)

            gamq_qsq(i, j, k) = gamq_factor(i, j, k) * gamq_qsq(i, j, k)
            gamq_cov(i, j, k) = gamq_factor(i, j, k) * gamq_cov(i, j, k)
            gamq_res(i, j, k) = gamq_factor(i, j, k) * gamq_res(i, j, k)

            pdc_cov(i, j, k) = elh                                      &
                          * (gamt_cov(i, j, k) * dqwdz(i, j, k)         &
                           + gamq_cov(i, j, k) * dtldz(i, j, k)) * 0.5
            IF (disp_coef < pdc_cov(i, j, k)) THEN
              pdc_factor(i, j, k) = disp_coef / pdc_cov(i, j, k)
            ELSE
              pdc_factor(i, j, k) = 1.0
            END IF
            gamt_tsq(i, j, k) = pdc_factor(i, j, k) * gamt_tsq(i, j, k)
            gamt_cov(i, j, k) = pdc_factor(i, j, k) * gamt_cov(i, j, k)
            gamt_res(i, j, k) = pdc_factor(i, j, k) * gamt_res(i, j, k)

            gamq_qsq(i, j, k) = pdc_factor(i, j, k) * gamq_qsq(i, j, k)
            gamq_cov(i, j, k) = pdc_factor(i, j, k) * gamq_cov(i, j, k)
            gamq_res(i, j, k) = pdc_factor(i, j, k) * gamq_res(i, j, k)

            gamt(i, j, k) = gamt_tsq(i, j, k) * tsq(i, j, k)            &
                          + gamt_cov(i, j, k) * cov(i, j, k)            &
                          + gamt_res(i, j, k)

            gamq(i, j, k) = gamq_qsq(i, j, k) * qsq(i, j, k)            &
                          + gamq_cov(i, j, k) * cov(i, j, k)            &
                          + gamq_res(i, j, k)

          END DO
        END DO
      END DO
    ELSE
      DO k = 2, tke_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            gamt(i, j, k) = gamt_tsq(i, j, k) * tsq(i, j, k)            &
                          + gamt_cov(i, j, k) * cov(i, j, k)            &
                          + gamt_res(i, j, k)

            gamq(i, j, k) = gamq_qsq(i, j, k) * qsq(i, j, k)            &
                          + gamq_cov(i, j, k) * cov(i, j, k)            &
                          + gamq_res(i, j, k)

            gamt_factor(i, j, k) = 1.0
            gamq_factor(i, j, k) = 1.0
            pdc_factor(i, j, k) = 1.0
          END DO
        END DO
      END DO
    END IF ! IF L_MY_PROD_ADJ

    ! Calculate production terms
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end

          elq = el(i, j, k) * qkw(i, j, k)
          elh = elq * qdiv(i, j, k)

          pdk(i, j, k) = elq * (sm(i, j, k) * gm(i, j, k)               &
                                  + sh(i, j, k) * gh(i, j, k))          &
                                  + wb_ng(i,j,k)

          pdt(i, j, k) = elh                                            &
                  * (sh(i, j, k) * dtldz(i, j, k) + gamt_res(i, j, k))  &
                  * dtldz(i, j, k)
          pdt_tsq(i, j, k) = elh * gamt_tsq(i, j, k) * dtldz(i, j, k)
          pdt_cov(i, j, k) = elh * gamt_cov(i, j, k) * dtldz(i, j, k)

          pdq(i, j, k) = elh                                            &
               * (sh(i, j, k) * dqwdz(i, j, k) + gamq_res(i, j, k))     &
               * dqwdz(i, j, k)
          pdq_qsq(i, j, k) = elh * gamq_qsq(i, j, k) * dqwdz(i, j, k)
          pdq_cov(i, j, k) = elh * gamq_cov(i, j, k) * dqwdz(i, j, k)

          pdc(i, j, k) = 0.5 * elh                                      &
               * ((sh(i, j, k) * dtldz(i, j, k)                         &
                       + gamt_res(i, j, k)) * dqwdz(i, j, k)            &
                + (sh(i, j, k) * dqwdz(i, j, k)                         &
                       + gamq_res(i, j, k)) * dtldz(i, j, k))

          pdc_tsq(i, j, k) = elh                                        &
                           * gamt_tsq(i, j, k) * dqwdz(i, j, k) * 0.5
          pdc_qsq(i, j, k) = elh                                        &
                           * gamq_qsq(i, j, k) * dtldz(i, j, k) * 0.5
          pdc_cov(i, j, k) = 0.5 * elh                                  &
                           * (gamt_cov(i, j, k) * dqwdz(i, j, k)        &
                            + gamq_cov(i, j, k) * dtldz(i, j, k))

          dfu_cg(i, j, k) = elq * smd(i, j, k) * dudz(i, j, k)
          dfv_cg(i, j, k) = elq * smd(i, j, k) * dvdz(i, j, k)
          dft_cg(i, j, k) = elq * gamt(i, j, k)
          dfq_cg(i, j, k) = elq * gamq(i, j, k)
        END DO
      END DO
    END DO
  ELSE  ! level 2.5
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          !     In Level 2.5, qdiv is not reset.
          gamt(i, j, k) = 0.0
          gamq(i, j, k) = 0.0
          gamv(i, j, k) = 0.0
          smd(i, j, k) = 0.0
          cu(i, j, k) = cu25(i, j, k)
          cv(i, j, k) = cv25(i, j, k)

          elq = el(i, j, k) * qkw(i, j, k)
          elh = elq * qdiv(i, j, k)

          pdk(i, j, k) = elq                                            &
                      * (sm(i, j, k) * gm(i, j, k)                      &
                       + sh(i, j, k) * gh(i, j, k)) + wb_ng(i,j,k)
          pdt(i, j, k) = elh                                            &
                      * (sh(i, j, k) * dtldz(i, j, k)) * dtldz(i, j, k)
          pdq(i, j, k) = elh                                            &
                      * (sh(i, j, k) * dqwdz(i, j, k)) * dqwdz(i, j, k)
          pdc(i, j, k) = elh                                            &
               * (sh(i, j, k) * dtldz(i, j, k)) * dqwdz(i, j, k) * 0.5  &
               + elh                                                    &
               * (sh(i, j, k) * dqwdz(i, j, k)) * dtldz(i, j, k) * 0.5

          dfm(i, j, k) = elq * sm(i, j, k)
          dfh(i, j, k) = elq * sh(i, j, k)
          dfu_cg(i, j, k) = 0.0
          dfv_cg(i, j, k) = 0.0
          dft_cg(i, j, k) = 0.0
          dfq_cg(i, j, k) = 0.0
        END DO
      END DO
    END DO
  END IF  ! test if levflag == 3

  ! Overwrite production terms by ones calculated with surface fluxes
  IF (my_lowest_pd_surf > 0) THEN
    IF (l_my_extra_level) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          pdk(i, j, 1) = 1.0 * u_s(i, j) ** 3 * pmz(i, j)               &
                      / (vkman * z_tq(i, j, 1) * my_z_extra_fact)
        END DO
      END DO
      IF (l_my_lowest_pd_surf_tqc) THEN
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            phm  = 1.0 / u_s(i, j) * phh(i, j)                          &
               / (vkman * z_tq(i, j, 1) * my_z_extra_fact)
            pdt(i, j, 1) = phm * ftl(i, j, 1) ** 2
            pdq(i, j, 1) = phm * fqw(i, j, 1) ** 2
            pdc(i, j, 1) = phm * ftl(i, j, 1) * fqw(i, j, 1)
          END DO
        END DO
      END IF
    ELSE    ! NOT L_MY_Extra_level
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          pdk(i, j, 2) = 1.0 * u_s(i, j) ** 3 * pmz(i, j)               &
                               / (vkman * z_tq(i, j, 1))
          pdk(i, j, 1) = 0.0
          pdt(i, j, 1) = 0.0
          pdq(i, j, 1) = 0.0
          pdc(i, j, 1) = 0.0
        END DO
      END DO
      IF (l_my_lowest_pd_surf_tqc) THEN
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            phm  = 1.0 / u_s(i, j)* phh(i, j)                           &
                              / (vkman * z_tq(i, j, 1))
            pdt(i, j, 2) = phm * ftl(i, j, 1) ** 2
            pdq(i, j, 2) = phm * fqw(i, j, 1) ** 2
            pdc(i, j, 2) = phm * ftl(i, j, 1) * fqw(i, j, 1)
            pdt_tsq(i, j, 2) = 0.0
            pdt_cov(i, j, 2) = 0.0
            pdq_qsq(i, j, 2) = 0.0
            pdq_cov(i, j, 2) = 0.0
            pdc_tsq(i, j, 2) = 0.0
            pdc_qsq(i, j, 2) = 0.0
            pdc_cov(i, j, 2) = 0.0
          END DO
        END DO
      END IF  ! IF L_MY_lowest_pd_surf_tqc
    END IF ! IF L_MY_EXTRA_LEVEL
  ELSE  ! MY_lowest_pd_surf = 0
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        pdk(i, j, 1) = 0.0
        pdt(i, j, 1) = 0.0
        pdq(i, j, 1) = 0.0
        pdc(i, j, 1) = 0.0
      END DO
    END DO
  END IF  ! IF MY_lowest_pd_surf

  ! for diagnostics
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      gamt(i, j, 1) = 0.0
      gamq(i, j, 1) = 0.0
      gamv(i, j, 1) = 0.0
      smd(i, j, 1) = 0.0
      dfh(i, j, 1) = 0.0
      dfu_cg(i, j, 1) = 0.0
      dfv_cg(i, j, 1) = 0.0
      dft_cg(i, j, 1) = 0.0
      dfq_cg(i, j, 1) = 0.0
      pdt_tsq(i, j, 1) = 0.0
      pdt_cov(i, j, 1) = 0.0
      pdq_qsq(i, j, 1) = 0.0
      pdq_cov(i, j, 1) = 0.0
      pdc_tsq(i, j, 1) = 0.0
      pdc_qsq(i, j, 1) = 0.0
      pdc_cov(i, j, 1) = 0.0
    END DO
  END DO

  IF (BL_diag%l_tke_shr_prod) THEN
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          elq = el(i, j, k) * qkw(i, j, k)
          BL_diag%tke_shr_prod(i, j, k) = el(i, j, k) * qkw(i, j, k)    &
                         * (sm(i, j, k) + smd(i, j, k)) * gm(i, j, k)
        END DO
      END DO
    END DO
  END IF

  IF (BL_diag%l_tke_boy_prod) THEN
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%tke_boy_prod(i, j, k) = el(i, j, k) * qkw(i, j, k)    &
                         * (sh(i, j, k) * gh(i, j, k)                   &
                                      + gamv(i, j, k)) + wb_ng(i,j,k)
        END DO
      END DO
    END DO
  END IF

  IF (BL_diag%l_tke_boy_prod) THEN
    DO k = k_start, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%tke_dissp(i, j, k) = qkw(i, j, k) ** 3                &
                                                / (b1 * el(i, j, k))
        END DO
      END DO
    END DO
  END IF

  IF (levflag == 3) THEN
    ! Integrate the covariances

    IF (imp_mode == full_impl) THEN
! DEPENDS ON: mym_update_covariance
      CALL mym_update_covariance(                                       &
! IN levels
            bl_levels,                                                  &
! IN fields
            qkw, el, dfm, pdt_tsq, pdt_cov, pdt,                        &
            pdq_qsq, pdq_cov, pdq, pdc_cov, pdc_tsq, pdc_qsq, pdc,      &
! INOUT fields
            tsq, qsq, cov)
    ELSE   ! half implict or explicit
      ALLOCATE(bp_tsq(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end, tke_levels))
      ALLOCATE(rp_tsq(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end, tke_levels))
      ALLOCATE(bp_qsq(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end, tke_levels))
      ALLOCATE(rp_qsq(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end, tke_levels))
      ALLOCATE(bp_cov(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end, tke_levels))
      ALLOCATE(rp_cov(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end, tke_levels))

      IF (imp_mode == half_impl) THEN
        DO k = k_start, tke_levels
          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              pdt(i, j, k) = pdt(i, j, k)                               &
                                  + pdt_cov(i, j, k) * cov(i, j, k)

              pdq(i, j, k) = pdq(i, j, k)                               &
                                  + pdq_cov(i, j, k) * cov(i, j, k)

              pdc(i, j, k) = pdc(i, j, k)                               &
                                  + pdc_tsq(i, j, k) * tsq(i, j, k)     &
                                  + pdc_qsq(i, j, k) * qsq(i, j, k)

              b2l = 2.0 * qkw(i, j, k) / (b2 * el(i, j, k))

              bp_tsq(i, j, k) = b2l - 2.0 * pdt_tsq(i, j, k)
              rp_tsq(i, j, k) = 2.0 * pdt(i, j, k)

              bp_qsq(i, j, k) = b2l - 2.0 * pdq_qsq(i, j, k)
              rp_qsq(i, j, k) = 2.0 * pdq(i, j, k)

              bp_cov(i, j, k) = b2l - 2.0 * pdc_cov(i, j, k)
              rp_cov(i, j, k) = 2.0 * pdc(i, j, k)
            END DO
          END DO
        END DO
      ELSE IF (imp_mode == expl) THEN
        DO k = k_start, tke_levels
          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              pdt(i, j, k) = pdt(i, j, k)                               &
                                  + pdt_tsq(i, j, k) * tsq(i, j, k)     &
                                  + pdt_cov(i, j, k) * cov(i, j, k)
              pdq(i, j, k) = pdq(i, j, k)                               &
                                  + pdq_qsq(i, j, k) * qsq(i, j, k)     &
                                  + pdq_cov(i, j, k) * cov(i, j, k)
              pdc(i, j, k) = pdc(i, j, k)                               &
                                  + pdc_cov(i, j, k) * cov(i, j, k)     &
                                  + pdc_tsq(i, j, k) * tsq(i, j, k)     &
                                  + pdc_qsq(i, j, k) * qsq(i, j, k)

              b2l = 2.0 * qkw(i, j, k) / (b2 * el(i, j, k))

              bp_tsq(i, j, k) = b2l
              rp_tsq(i, j, k) = 2.0 * pdt(i, j, k)

              bp_qsq(i, j, k) = b2l
              rp_qsq(i, j, k) = 2.0 * pdq(i, j, k)

              bp_cov(i, j, k) = b2l
              rp_cov(i, j, k) = 2.0 * pdc(i, j, k)
            END DO
          END DO
        END DO
      END IF
! DEPENDS ON: mym_update_fields
      CALL mym_update_fields(                                           &
            bl_levels, coef_trbvar_diff,dfm, rp_tsq, bp_tsq,tsq)

! DEPENDS ON: mym_update_fields
      CALL mym_update_fields(                                           &
            bl_levels, coef_trbvar_diff,dfm, rp_qsq, bp_qsq,qsq)

! DEPENDS ON: mym_update_fields
      CALL mym_update_fields(                                           &
            bl_levels, coef_trbvar_diff,dfm, rp_cov, bp_cov,cov)

      DEALLOCATE(rp_cov)
      DEALLOCATE(bp_cov)
      DEALLOCATE(rp_qsq)
      DEALLOCATE(bp_qsq)
      DEALLOCATE(rp_tsq)
      DEALLOCATE(bp_tsq)

    END IF  ! if imp_mode == FULL_IMPL
  ELSE  ! level 2.5
    ! In level 2.5, tsq, qsq, cov are diagnosed assuming balance between
    ! prodcution and dissipation.
    DO k = k_start, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          IF (qkw(i, j, k) <= 1.0e-4)  THEN
            b2l = 0.0
          ELSE
            b2l = b2 * el(i, j, k) / qkw(i, j, k)
          END IF
          tsq(i, j, k) = b2l * 2.0 * pdt(i, j, k)
          qsq(i, j, k) = b2l * 2.0 * pdq(i, j, k)
          cov(i, j, k) = b2l * 2.0 * pdc(i, j, k)
        END DO
      END DO
    END DO
  END IF

  IF (levflag >= 2) THEN
    ! predict qke
    IF (my_lowest_pd_surf > 0) THEN
      k_start_cor = k_start + 1
    ELSE
      k_start_cor = k_start
    END IF

    IF (levflag == 3 .AND.                                              &
           (imp_mode == half_impl .OR. imp_mode == full_impl)) THEN
      ! add correction terms evaluated with integrated tsq, qsq and cov
      DO k = k_start_cor, tke_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end

            t3sq = MAX(tsq(i, j, k), 0.0)
            r3sq = MAX(qsq(i, j, k), 0.0)
            c3sq = cov(i, j, k)

            c3sq = SIGN( MIN( ABS(c3sq), SQRT(t3sq*r3sq) ), c3sq )

            t3sq = vt(i, j, k) * t3sq + vq(i, j, k) * c3sq
            r3sq = vt(i, j, k) * c3sq + vq(i, j, k) * r3sq
            c3sq = MAX(vt(i, j, k) * t3sq + vq(i, j, k) * r3sq, 0.0)

            elq = el(i, j, k) * qkw(i, j, k)
            smd(i, j, k) = smd_coef(i, j, k) * (c3sq - c2sq(i, j, k))

            pdk(i, j, k) = pdk(i, j, k) + elq                           &
                 * (smd(i, j, k) * gm(i, j, k)                          &
                    + gamv_coef(i, j, k) * (c3sq- c2sq(i, j, k)))
          END DO
        END DO
      END DO
    ELSE
      DO k = k_start_cor, tke_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            pdk(i, j, k) = pdk(i, j, k)                                 &
                  + el(i, j, k) * qkw(i, j, k)                          &
                        * (smd(i, j, k) * gm(i, j, k) + gamv(i, j, k))
          END DO
        END DO
      END DO
    END IF ! if test levflag == 3

    DO k = k_start, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          b1l = b1 * el(i, j, k)
          bp(i, j, k) = 2.0 * qkw(i, j, k) / b1l
          rp(i, j, k) = 2.0 * pdk(i, j, k)
        END DO
      END DO
    END DO

! DEPENDS ON: mym_update_fields
    CALL mym_update_fields(                                             &
          bl_levels, coef_trbvar_diff_tke,dfm, rp, bp, qke)
  ELSE
     ! level 2
     ! diagnose qke
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          b2l = b2 * el(i, j, k)
          qke(i, j, k) = (MAX(b2l * 2.0 * pdk(i, j, k), 0.0))           &
                                                       ** two_thirds
        END DO
      END DO
    END DO
  END IF  ! test if levflag >= 2

  DO k = 1, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qke(i, j, k) = MIN(MAX(qke(i, j, k), 1.e-20), qke_max)
        tsq(i, j, k) = MAX(tsq(i, j, k), 0.0)
        qsq(i, j, k) = MAX(qsq(i, j, k), 0.0)
      END DO
    END DO
  END DO

  DO k = tke_levels + 1, bl_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        qke(i, j, k) = 0.0
        tsq(i, j, k) = 0.0
        qsq(i, j, k) = 0.0
        cov(i, j, k) = 0.0
        dfm(i, j, k) = 0.0
        dfh(i, j, k) = 0.0
        dfu_cg(i, j, k) = 0.0
        dfv_cg(i, j, k) = 0.0
        dft_cg(i, j, k) = 0.0
        dfq_cg(i, j, k) = 0.0
      END DO
    END DO
  END DO


  IF (BL_diag%l_elm) THEN
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%elm(i, j, k) = el(i, j, k)
        END DO
      END DO
    END DO
  END IF

  IF (BL_diag%l_sm) THEN
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%sm(i, j, k) = sm(i, j, k)
        END DO
      END DO
    END DO
  END IF

  IF (BL_diag%l_sh) THEN
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%sh(i, j, k) = sh(i, j, k)
        END DO
      END DO
    END DO
  END IF

  IF (lhook) CALL dr_hook('MYM_TURBULENCE',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE mym_turbulence

