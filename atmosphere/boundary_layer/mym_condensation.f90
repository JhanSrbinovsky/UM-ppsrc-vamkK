! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_condensation----------------------------------------
!
!  Purpose: To evaluate buoyancy parameters, cloud fraction, and
!           liquid water content in the MY model.
!           Fluctuation of heat and moisture is assumed to obey the
!           bi-normal distribution function.
!           Note that the cloud fraction and liquid water content
!           diagnosed here is not used in the other processes.
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
SUBROUTINE mym_condensation(                                            &
! IN levels/switches
      bl_levels, levflag, nSCMDpkgs,L_SCMDiags, l_mixing_ratio, BL_diag,&
! IN fields
      qw, tl, t, p_theta_levels, tsq, qsq, cov,                         &
! OUT fields
      vt, vq, q1, cld, ql)
      
  USE atm_fields_bounds_mod, ONLY: tdims, tdims_l
  USE conversions_mod, ONLY: pi
  USE atmos_constants_mod, ONLY:                                        &
      cp, r, repsilon, pref, kappa, c_virtual
  USE water_constants_mod, ONLY: lc, lf
  USE mym_option_mod, ONLY: tke_levels
  USE bl_diags_mod, ONLY: strnewbldiag
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

! Intent In Variables
  INTEGER, INTENT(IN) ::                                                &
     bl_levels,                                                         &
                    ! Max. no. of "boundary" levels
     levflag
                    ! flag to indicate the level of MY
                    ! 2: MY2.5, 3:MY3

  REAL, INTENT(IN) ::                                                   &
     qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
                    ! Total water content
     tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels), &
                    ! Ice/liquid water temperature
     t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                    ! temperature
     p_theta_levels(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
                    bl_levels+1),                                       &
                    ! pressure on theta levels (Pa)
     tsq(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
            bl_levels),                                                 &
                    ! Self covariant of liquid potential temperature
                    ! (thetal'**2) defined on theta levels K-1
     qsq(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
            bl_levels),                                                 &
                    ! Self covariant of total water
                    ! (qw'**2) defined on theta levels K-1
     cov(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,   &
            bl_levels)
                    ! Correlation between thetal and qw
                    ! (thetal'qw') defined on theta levels K-1

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER, INTENT(IN) ::  nSCMDpkgs
                    ! No of SCM diagnostics packages

  LOGICAL, INTENT(IN) ::  L_SCMDiags(nSCMDpkgs)
                    ! Logicals for SCM diagnostics packages

  LOGICAL, INTENT(IN) ::                                                &
     l_mixing_ratio
                    !  flag for qsat_mix
                    !  .true. return qsat as a mixing ratio
                    !  .false. return qsat as a specific humidity

!  Declaration of BL diagnostics.
  TYPE (strnewbldiag), INTENT(INOUT) :: BL_diag

! Intent OUT Variables
  REAL, INTENT(OUT) ::                                                  &
     vt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                    ! Buoyancy parameter (coefficients of <w'thetal'>)
                    ! on theta K-1
                    ! Note that g/thetav is not included.
     vq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                    ! Buoyancy parameter (coefficients of <w'qw'>)
                    ! on theta K-1
                    ! Note that g/thetav is not included.
     q1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                    ! normalized excess water content on theta K-1
     cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                    ! cloud fraction on theta K-1
     ql(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, tke_levels)
                    ! condensed liquid water content

! Local Variables
  INTEGER ::                                                            &
     i, j, k
                    ! loop indexes
  REAL ::                                                               &
     rr2,                                                               &
                    ! 1 / sqrt(2)
     rrp,                                                               &
                    ! 1 / sqrt(2 pi)
     hl,                                                                &
                    ! latent heat
     qsl,                                                               &
                    ! saturated specific humidity
     dqsl,                                                              &
                    ! derivative of saturated specific humidity by
                    ! temperature
     t3sq,                                                              &
                    ! work variable for <thetal'**2>
     r3sq,                                                              &
                    ! work variable for <qw'**2>
     c3sq,                                                              &
                    ! work variable for <thetal'qw'>
     alp_qsl,                                                           &
                    ! alpha * qsl
     eq1,                                                               &
                    ! work variable
     qll,                                                               &
                    ! work variable
     r_exner,                                                           &
                    ! reciprocal of Exner function
     q2p,                                                               &
                    ! work variable
     pt_tmp,                                                            &
                    ! work variable
     qt,                                                                &
                    ! work variable
     rac
                    ! work variable







  REAL ::                                                               &
     rice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                    ! ratio of ice.
                    ! temperature > 0C   : rice =0
                    ! temperature < -36C : rice =1
                    ! Between 0C and 36C, it is linearly interpolated
                    ! with temperature
     hl_ovr_cp(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               tke_levels),                                             &
                    ! latent heat over heat capacity
     exner(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           tke_levels),                                                 &
                    ! Exner function
     qmq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                    ! Excess of total water from saturated one
     alp(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                    ! coefficient related to saturation
                    ! (alpha in the paper)
     bet(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                    ! coefficient related to saturation
                    ! (beta in the paper)
     sgm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                    ! standard deviation of the bi-normal distribution
                    ! function
     erf_arg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                    ! array to store an argument for error function
     erf_val(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             tke_levels),                                               &
                    ! array to store a result of error function
     qsw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         tke_levels),                                                   &
                    ! saturated specific ratio for water
     qsi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels)
                    ! saturated specific ratio for ice

  REAL, PARAMETER ::                                                    &
     e0cw = 6.11e2,                                                     &
     tetn1w = 17.27,                                                    &
     tetn2w = 273.15,                                                   &
     tetn3w = 35.85,                                                    &
     e0ci = 6.11e2,                                                     &
     tetn1i = 21.875,                                                   &
     tetn2i = 273.15,                                                   &
     tetn3i = 7.65,                                                     &
                    ! coefficients in the Tetens' Formula
     ttriple = 273.16,                                                  &
                    ! a triple point of water
     temp_zero = 273.15,                                                &
                    ! Kelvin for O degree Celsius
     temp_ice = 237.15
                    ! Below this temperature, all of condensed water
                    ! should be ice. -36C
  REAL, PARAMETER ::                                                    &
     my_sgm_min_fct = 0.0,                                              &
                    ! factor to set the lower limit for sgm
     my_sgm_max_fct = 1.0
                    ! factor to set the upper limit for sgm

! Derived local parameters.

  REAL, PARAMETER ::                                                    &
     ls = lf+lc,                                                        &
          ! Latent heat of sublimation.
     one_minus_epsilon = 1.0 - repsilon
          ! 1 - repsilon = 0.378






  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_CONDENSATION',zhook_in,zhook_handle)
  rr2 = 1.0 / SQRT(2.0)
  rrp = 1.0 / SQRT(2.0 * pi)

  ! Here, qsw and qsi are saturated vapor pressure.
  ! Using the Teten's formula instead of the subroutine "qmix"
  ! because the saturated vapor pressure on liquid water is necessary
  ! even in sub-zero temperature.
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        qsw(i, j, k) = e0cw * EXP(tetn1w *                              &
              (tl(i, j, k - 1)  - tetn2w)                               &
                    / (tl(i, j, k - 1) - tetn3w) )
        qsi(i, j, k) = e0ci * EXP(tetn1i *                              &
              (tl(i, j, k - 1) - tetn2i)                                &
             / (tl(i, j, k - 1) - tetn3i) )
      END DO
    END DO
  END DO

  ! convert to mixing ratio or specific humidity
  IF (l_mixing_ratio) THEN
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qsw(i, j, k) = repsilon * qsw(i, j, k)                        &
                                    / p_theta_levels(i, j, k - 1)
          qsi(i, j, k) = repsilon * qsi(i, j, k)                        &
                                    / p_theta_levels(i, j, k - 1)
        END DO
      END DO
    END DO
  ELSE
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qsw(i, j, k) = repsilon * qsw(i, j, k)                        &
                          / (p_theta_levels(i, j, k - 1)                &
                               - one_minus_epsilon * qsw(i, j, k))
          qsi(i, j, k) = repsilon * qsi(i, j, k)                        &
                          / (p_theta_levels(i, j, k - 1)                &
                               - one_minus_epsilon * qsi(i, j, k))
        END DO
      END DO
    END DO
  END IF

  ! Calculate sgm
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (tl(i, j, k -1) >= ttriple) THEN
          rice(i, j, k) = 0.0
        ELSE IF (tl(i, j, k - 1) < temp_ice) THEN
          rice(i, j, k) = 1.0
        ELSE
          rice(i, j, k) = (ttriple - tl(i, j, k - 1))                   &
                                         / (ttriple - temp_ice)
        END IF

        hl = (1.0 - rice(i, j, k)) * lc + rice(i, j, k) * ls
        qsl = (1.0 - rice(i, j, k)) * qsw(i, j, k)                      &
                + rice(i, j, k) * qsi(i, j, k)

        hl_ovr_cp(i, j, k) = hl / cp

        dqsl = qsl * repsilon * hl / (r * tl(i, j, k - 1) **2)

        exner(i, j, k) = (p_theta_levels(i, j, k - 1) / pref) ** kappa
        qmq(i, j, k) = qw(i, j, k - 1) - qsl
        alp(i, j, k) = 1.0 /(1.0 + dqsl * hl_ovr_cp(i, j, k))
        bet(i, j, k) = dqsl * exner(i, j, k)

        t3sq = MAX(tsq(i, j, k), 0.0)
        r3sq = MAX(qsq(i, j, k), 0.0)
        c3sq = cov(i, j, k)
        c3sq = SIGN(MIN(ABS(c3sq), SQRT(t3sq * r3sq)), c3sq)

        r3sq = r3sq + bet(i, j, k) ** 2 * t3sq                          &
             -2.0 * bet(i, j, k) * c3sq
        alp_qsl = MIN(alp(i, j, k) * qsl, qw(i, j, k - 1))
        sgm(i, j, k) = MAX(                                             &
                        MIN(0.5 * alp(i, j, k) * SQRT(MAX(r3sq, 0.0)),  &
                        my_sgm_max_fct * alp_qsl),                      &
                        my_sgm_min_fct * alp_qsl, 1.e-10)
      END DO
    END DO
  END DO

  IF (levflag /= 3) THEN
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        sgm(i, j, 2) = sgm(i, j, 3)
      END DO
    END DO
  END IF

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      erf_arg(i, j, 1) = 0.0
      cld(i, j, 1) = 0.0
      ql(i, j, 1) = 0.0
      sgm(i, j, 1) = 0.0
      q1(i, j, 1) = 0.0
    END DO
  END DO
  !
  ! Preparation to calculate values of the err function
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        q1(i, j, k) = 0.5 * alp(i, j, k)                                &
                               * qmq(i, j, k) / sgm(i, j, k)
        erf_arg(i, j, k) = q1(i, j, k) * rr2
      END DO
    END DO
  END DO

! DEPENDS ON: mym_errfunc
  CALL mym_errfunc(tdims%i_end*tdims%j_end*tke_levels, erf_arg, erf_val)

  ! Calculate the buoyancy parameters vt and vq
  DO k = 2, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        cld(i, j, k) = 0.5 * (1.0 + erf_val(i, j, k))
        IF(ABS(q1(i, j, k)) > 10.0 ) THEN
          eq1 = 0.0
        ELSE
          eq1  = rrp * EXP(- 0.5 * q1(i, j, k) ** 2)
        END IF
        ! qll = ql / (2 * sgm)
        qll  = MAX(cld(i, j, k) * q1(i, j, k) + eq1, 0.0)

        IF(qw(i, j, k) < 1.e-10) THEN
          ql(i, j, k) = 0.0
        ELSE
          ql(i, j, k) = MAX(                                            &
               2. * sgm(i, j, k) * qll, 0.0)
        END IF
        ! To avoid negative QV (for safety)
        ql(i, j, k) = MIN(ql(i, j, k), qw(i, j, k - 1) * 0.5)

        r_exner = 1.0 / exner(i, j, k)
        q2p  = hl_ovr_cp(i, j, k) * r_exner
        pt_tmp = t(i, j, k - 1) * r_exner
        qt = 1.0 + c_virtual * qw(i, j, k - 1)                          &
                 - (1.0 + c_virtual) * ql(i, j, k)
        rac = alp(i, j, k) * (cld(i, j, k)                              &
             - qll * eq1) * (q2p * qt - (1.0 + c_virtual) * pt_tmp)

        vt (i, j, k) = qt - rac * bet(i, j, k)
        vq (i, j, k) = c_virtual * pt_tmp + rac
      END DO
    END DO
  END DO

  IF (BL_diag%l_cf_trb) THEN
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%cf_trb(i, j, k) = cld(i, j, k)
        END DO
      END DO
    END DO
  END IF

  IF (BL_diag%l_ql_trb) THEN
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%ql_trb(i, j, k) = ql(i, j, k)
        END DO
      END DO
    END DO
  END IF

  IF (BL_diag%l_sgm_trb) THEN
    DO k = 2, tke_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%sgm_trb(i, j, k) = sgm(i, j, k)
        END DO
      END DO
    END DO
  END IF


  IF (lhook) CALL dr_hook('MYM_CONDENSATION',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE mym_condensation

