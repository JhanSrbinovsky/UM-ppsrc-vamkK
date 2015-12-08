! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_initialize------------------------------------------
!
!  Purpose: To set initial values of the prognostic variables
!           appeared in MYmodel (E_TRB, TSQ, QSQ, COV) with MY level
!           2 model iteration, that is, assuming balance between
!           production and dissipation.
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
SUBROUTINE mym_initialize(                                              &
! IN levels
      bl_levels,                                                        &
! IN fields
      z_uv, z_tq, dbdz, dvdzm, dtldzm, dqwdzm,                          &
      fqw, ftl, u_s, r_mosurf, fb_surf,                                 &
! INOUT fields
      e_trb, tsq, qsq, cov)

  USE atm_fields_bounds_mod, ONLY: tdims, tdims_l, pdims, tdims_s
  USE atmos_constants_mod, ONLY:                                        &
      vkman, cp
  USE mym_const_mod, ONLY: b1, b2, qke_max, coef_trbvar_diff,           &
        coef_trbvar_diff_tke
  USE mym_option_mod, ONLY:                                             &
        my_lowest_pd_surf, l_my_extra_level, my_z_extra_fact,           &
        tke_levels, l_my_lowest_pd_surf_tqc
  USE earth_constants_mod, ONLY: g
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

! Intent IN Variables
  INTEGER, INTENT(IN) ::                                                &
     bl_levels
                    ! Max. no. of "boundary" levels

  REAL, INTENT(IN) ::                                                   &
     z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels),                                                   &
                    ! Z_UV(*,K) is height of u level k
     z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                    ! Z_TQ(*,K) is height of theta
                    ! level k.
     dbdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          2:tke_levels),                                                &
                    ! Buoyancy gradient across layer
                    ! interface interpolated to theta levels.
                    ! (:,:,K) repserents the value on theta level K-1
     dvdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           2:bl_levels),                                                &
                    ! Modulus of wind shear at theta levels.
                    ! (:,:,K) repserents the value on theta level K-1
     dtldzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            2:bl_levels),                                               &
                    ! gradient of TL across layer
                    ! interface interpolated to theta levels.
                    ! (:,:,K) repserents the value on theta level K-1
     dqwdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            2:bl_levels),                                               &
                    ! gradient of QW across layer
                    ! interface interpolated to theta levels.
                    ! (:,:,K) repserents the value on theta level K-1
     fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                    ! Moisture flux between layers
                    ! (kg per square metre per sec).
                    ! FQW(,1) is total water flux
                    ! from surface, 'E'.
     ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                    ! FTL(,K) contains net turbulent
                    ! sensible heat flux into layer K
                    ! from below; so FTL(,1) is the
                    ! surface sensible heat, H. (W/m2)
     u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
                    ! Surface friction velocity
     r_mosurf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
     fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                    ! Surface flux buoyancy over
                    ! density (m^2/s^3)

! Intent INOUT Variables
  REAL, INTENT(INOUT) ::                                                &
     e_trb(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end, &
            bl_levels),                                                 &
                    ! TKE defined on theta levels K-1
     tsq(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
            bl_levels),                                                 &
                    ! Self covariance of liquid potential temperature
                    ! (thetal'**2) defined on theta levels K-1
     qsq(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
            bl_levels),                                                 &
                    ! Self covariance of total water
                    ! (qw'**2) defined on theta levels K-1
     cov(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
            bl_levels)
                    ! Correlation between thetal and qw
                    ! (thetal'qw') defined on theta levels K-1

! Local variables
  INTEGER ::                                                            &
     i, j, k, ll
                   ! Loop indexes

  REAL ::                                                               &
     phm,                                                               &
                   ! gradient function at the surface
     elq
                   ! mixing length * qkw

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
     el(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                   ! mixing length on theta level K-1
     qkw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                   ! q=sqrt(qke) on theta level K-1
     pmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
                   ! gradient function for momentum minus zeta
                   ! (zeta: height over Monin-Obkhov length)
     phh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                   ! gradient function for scalars


  REAL ::                                                               &
     pdk(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                   ! production terms of qke divided by elq
     pdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                   ! production terms of tsq divided by elq
     pdq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                   ! production terms of qsq divided by elq
     pdc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                   ! production terms of cov divided by elq
     pdk0(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                   ! production terms of qke at the lowest level
     pdt0(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                   ! production terms of tsq at the lowest level
     pdq0(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                   ! production terms of qsq at the lowest level
     pdc0(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                   ! production terms of cov at the lowest level
     aa_qke(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! coefficients of the tri-diagonal eqs. for qke
     bb_qke(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! coefficients of the tri-diagonal eqs. for qke
     cc_qke(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! coefficients of the tri-diagonal eqs. for qke
     aa_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! coefficients of the tri-diagonal eqs. for tsq
     bb_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! coefficients of the tri-diagonal eqs. for tsq
     cc_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! coefficients of the tri-diagonal eqs. for tsq
     aa_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! coefficients of the tri-diagonal eqs. for qsq
     bb_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! coefficients of the tri-diagonal eqs. for qsq
     cc_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! coefficients of the tri-diagonal eqs. for qsq
     aa_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! coefficients of the tri-diagonal eqs. for cov
     bb_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! coefficients of the tri-diagonal eqs. for cov
     cc_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! coefficients of the tri-diagonal eqs. for cov
     aa_oth(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! work variable for aa_tsq, aa_qsq and aa_cov
     bb_oth(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! work variable for bb_tsq, bb_qsq and bb_cov
     cc_oth(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            tke_levels),                                                &
                   ! work variable for cc_tsq, cc_qsq and cc_cov
     qke_nohalo(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                tke_levels),                                            &
                   ! qke without halos
     tsq_nohalo(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                tke_levels),                                            &
                   ! tsq without halos
     qsq_nohalo(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                tke_levels),                                            &
                   ! qsq without halos
     cov_nohalo(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                tke_levels),                                            &
                   ! cov without halos
     dfm(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,   &
          bl_levels)
                   ! diffusion coefficient for momentum
                   ! on theta level K-1

  INTEGER ::                                                            &
     my3_itr_ini
                   ! number of iteration
  INTEGER :: k_start

! Combined constants
  REAL, PARAMETER ::                                                    &
     grcp = g/cp

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_INITIALIZE',zhook_in,zhook_handle)

  my3_itr_ini = tke_levels + 1

  IF (my_lowest_pd_surf == 0) THEN
    l_my_extra_level = .FALSE.
    my_z_extra_fact = 1.0
  END IF

  IF (l_my_extra_level) THEN
    k_start = 1
  ELSE
    k_start = 2
  END IF

  DO k = 1, bl_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        e_trb(i, j, k) = 0.0
        tsq(i, j, k) = 0.0
        qsq(i, j, k) = 0.0
        cov(i, j, k) = 0.0
      END DO
    END DO
  END DO

! DEPENDS ON: mym_level2
  CALL mym_level2(                                                      &
        bl_levels, dbdz, dvdzm, gm, gh, sm, sh)

  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        pdk(i, j, k) = sm(i, j, k) * gm(i, j, k)                        &
                                         + sh(i, j, k) * gh(i, j, k)
        IF (pdk(i, j, k) <= 0.0) THEN
          qke_nohalo(i, j, k) = 0.0
          pdk(i, j, k) = 0.0
          pdt(i, j, k) = 0.0
          pdq(i, j, k) = 0.0
          pdc(i, j, k) = 0.0
        ELSE
          qke_nohalo(i, j, k) = 1.0e-5
          pdt(i, j, k) = sh(i, j, k) * dtldzm(i, j, k) ** 2
          pdq(i, j, k) = sh(i, j, k) * dqwdzm(i, j, k) ** 2
          pdc(i, j, k) = sh(i, j, k) * dtldzm(i, j, k) * dqwdzm(i, j, k)
        END IF
      END DO
    END DO
  END DO

  IF (my_lowest_pd_surf > 0) THEN
! DEPENDS ON: mym_calcphi
    CALL mym_calcphi(                                                   &
          bl_levels, z_tq, r_mosurf, pmz, phh)
    IF (l_my_extra_level) THEN
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          pdk0(i, j) = 1.0 * u_s(i, j) ** 3 * pmz(i, j)                 &
                      / (vkman * z_tq(i, j, 1) * my_z_extra_fact)
        END DO
      END DO

      IF (l_my_lowest_pd_surf_tqc) THEN
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            phm  = 1.0 / u_s(i, j) * phh(i, j)                          &
                         / (vkman * z_tq(i, j, 1) * my_z_extra_fact)
            pdt0(i, j) = phm * ftl(i, j, 1) ** 2
            pdq0(i, j) = phm * fqw(i, j, 1) ** 2
            pdc0(i, j) = phm * ftl(i, j, 1) * fqw(i, j, 1)
          END DO
        END DO
      END IF
    ELSE
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          pdk0(i, j) = 1.0 * u_s(i, j) ** 3 * pmz(i, j)                 &
                      / (vkman * z_tq(i, j, 1))
        END DO
      END DO

      IF (l_my_lowest_pd_surf_tqc) THEN
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            phm  = 1.0 / u_s(i, j)* phh(i, j)                           &
                             / (vkman * z_tq(i, j, 1))
            pdt0(i, j) = phm * ftl(i, j, 1) ** 2
            pdq0(i, j) = phm * fqw(i, j, 1) ** 2
            pdc0(i, j) = phm * ftl(i, j, 1) * fqw(i, j, 1)
          END DO
        END DO
      END IF
    END IF ! IF L_MY_EXTRA_LEVEL
  END IF  ! IF MY_lowest_pd_surf

  DO ll = 1, my3_itr_ini
! DEPENDS ON: mym_length
    CALL mym_length(                                                    &
          tdims%i_end, tdims%j_end, 0, 0, bl_levels,                    &
          qke_nohalo, z_uv, z_tq, dbdz, r_mosurf, fb_surf,              &
          qkw, el)

    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          IF (qke_nohalo(i, j, k) <= 0.0) THEN
            qkw(i, j, k) = 0.0
          END IF
          dfm(i, j, k) = sm(i, j, k) * qkw(i, j, k) * el(i, j, k)
        END DO
      END DO
    END DO

! DEPENDS ON: mym_diff_matcoef
    CALL mym_diff_matcoef(                                              &
          bl_levels, coef_trbvar_diff_tke, dfm,                         &
          aa_qke, bb_qke, cc_qke)

! DEPENDS ON: mym_diff_matcoef
    CALL mym_diff_matcoef(                                              &
          bl_levels, coef_trbvar_diff, dfm,                             &
          aa_oth, bb_oth, cc_oth)

    DO k = k_start, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          IF (bb_qke(i, j, k) == 0.0) THEN
            aa_qke(i, j, k) = 0.0
            bb_qke(i, j, k) = 1.0
            cc_qke(i, j, k) = 0.0
            qke_nohalo(i, j, k) = 0.0

            aa_tsq(i, j, k) = 0.0
            bb_tsq(i, j, k) = 1.0
            cc_tsq(i, j, k) = 0.0
            tsq_nohalo(i, j, k) = 0.0

            aa_qsq(i, j, k) = 0.0
            bb_qsq(i, j, k) = 1.0
            cc_qsq(i, j, k) = 0.0
            qsq_nohalo(i, j, k) = 0.0

            aa_cov(i, j, k) = 0.0
            bb_cov(i, j, k) = 1.0
            cc_cov(i, j, k) = 0.0
            cov_nohalo(i, j, k) = 0.0
          ELSE
            elq = qkw(i, j, k) * el(i, j, k)
            aa_qke(i, j, k) = - aa_qke(i, j, k)
            bb_qke(i, j, k) = - bb_qke(i, j, k)                         &
                              + 2.0 * qkw(i, j, k) / (b1 * el(i, j, k))
            bb_qke(i, j, k) = SIGN(MAX(ABS(bb_qke(i, j, k)), 1.0e-20),  &
                                      bb_qke(i, j, k))

            cc_qke(i, j, k) = - cc_qke(i, j, k)
            qke_nohalo(i, j, k) = 2.0 * elq * pdk(i, j, k)

            aa_oth(i, j, k) = - aa_oth(i, j, k)
            bb_oth(i, j, k) = - bb_oth(i, j, k)                         &
                              + 2.0 * qkw(i, j, k) / (b2 * el(i, j, k))
            bb_oth(i, j, k) = SIGN(MAX(ABS(bb_oth(i, j, k)), 1.0e-20),  &
                                      bb_oth(i, j, k))
            cc_oth(i, j, k) = - cc_oth(i, j, k)

            aa_tsq(i, j, k) = aa_oth(i, j, k)
            bb_tsq(i, j, k) = bb_oth(i, j, k)
            cc_tsq(i, j, k) = cc_oth(i, j, k)
            tsq_nohalo(i, j, k) = 2.0 * elq * pdt(i, j, k)

            aa_qsq(i, j, k) = aa_oth(i, j, k)
            bb_qsq(i, j, k) = bb_oth(i, j, k)
            cc_qsq(i, j, k) = cc_oth(i, j, k)
            qsq_nohalo(i, j, k) = 2.0 * elq * pdq(i, j, k)

            aa_cov(i, j, k) = aa_oth(i, j, k)
            bb_cov(i, j, k) = bb_oth(i, j, k)
            cc_cov(i, j, k) = cc_oth(i, j, k)
            cov_nohalo(i, j, k) = 2.0 * elq * pdc(i, j, k)
          END IF
        END DO
      END DO
    END DO

    IF (my_lowest_pd_surf > 0) THEN
      k = k_start
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          IF (bb_qke(i, j, k) /= 0.0) THEN
            qke_nohalo(i, j, k) = 2.0 * pdk0(i, j)
          END IF
        END DO
      END DO

      IF (l_my_lowest_pd_surf_tqc) THEN
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            IF (bb_qke(i, j, k) /= 0.0) THEN
              tsq_nohalo(i, j, k) = 2.0 * pdt0(i, j)
              qsq_nohalo(i, j, k) = 2.0 * pdq0(i, j)
              cov_nohalo(i, j, k) = 2.0 * pdc0(i, j)
            END IF
          END DO
        END DO
      END IF ! IF L_MY_lowest_pd_surf_tqc
    END IF   ! IF MY_lowest_pd_surf > 0

! DEPENDS ON: mym_implic
    CALL mym_implic(                                                    &
                    tke_levels, k_start, tke_levels,                    &
                    aa_qke, bb_qke, cc_qke, qke_nohalo)

! DEPENDS ON: mym_implic
    CALL mym_implic(                                                    &
                    tke_levels, k_start, tke_levels,                    &
                    aa_tsq, bb_tsq, cc_tsq, tsq_nohalo)

! DEPENDS ON: mym_implic
    CALL mym_implic(                                                    &
                    tke_levels, k_start, tke_levels,                    &
                    aa_qsq, bb_qsq, cc_qsq, qsq_nohalo)

! DEPENDS ON: mym_implic
    CALL mym_implic(                                                    &
                    tke_levels, k_start, tke_levels,                    &
                    aa_cov, bb_cov, cc_cov, cov_nohalo)

  END DO  ! iteration ll = 1, my3_itr_ini

  DO k = k_start, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        e_trb(i, j, k) = 0.5 * MIN(                                     &
                                    MAX(qke_nohalo(i, j, k), 1.0e-20),  &
                                    qke_max)
        tsq(i, j, k) = MAX(tsq_nohalo(i, j, k), 0.0)
        qsq(i, j, k) = MAX(qsq_nohalo(i, j, k), 0.0)
        cov(i, j, k) = cov_nohalo(i, j, k)
      END DO
    END DO
  END DO

  IF (lhook) CALL dr_hook('MYM_INITIALIZE',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE mym_initialize

