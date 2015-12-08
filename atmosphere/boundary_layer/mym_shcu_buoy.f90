! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE mym_shcu_buoy-------------------------------------------
!
!  Purpose: To evaluate buoyancy flux in shallow convection
!           corresponding to the skewness of the distribution function
!           in the MY model.
!
!  Programming standard : unified model documentation paper No 3
!
!  Documentation: UMDP 024A
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE mym_shcu_buoy(                                               &
! IN levels/switches
                     bl_levels,nSCMDpkgs, L_SCMDiags, l_mixing_ratio,   &
                     BL_diag,                                           &
! IN fields
                     fb_surf, ustar, pstar,                             &
                     z_tq, z_uv, p_theta_levels, p_rho_levels,          &
                     u_p, v_p, t, q, qcl, qcf, q1, frac_gauss,          &
! INOUT / OUT fields
                     zhpar,frac, wb_ng)
                     
  USE atm_fields_bounds_mod, ONLY: tdims, pdims
  USE conversions_mod, ONLY: pi
  USE timestep_mod, ONLY: timestep
  USE mym_option_mod, ONLY: tke_levels, wb_ng_max, shcu_levels
  USE atmos_constants_mod, ONLY:                                        &
      cp, r, repsilon, pref, kappa, c_virtual, recip_kappa
  USE water_constants_mod, ONLY: lc, lf, tm    
  USE bl_diags_mod, ONLY: strnewbldiag
  USE earth_constants_mod, ONLY: g
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  IMPLICIT NONE

  INTEGER, INTENT(IN) ::                                                &
     bl_levels
                 ! Max. no. of "boundary" levels

  REAL, INTENT(IN) ::                                                   &
     fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),      &
                ! buoyancy flux at the surface
     ustar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                ! surface friction velocity
     pstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                ! surface pressure (Pa)

  REAL, INTENT(IN) ::                                                   &
     z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                ! height of theta levels
     z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          bl_levels),                                                   &
                ! height of p levels
     p_theta_levels(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
                    bl_levels + 1),                                     &
                ! pressure at theta levels (Pa)
     p_rho_levels(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  bl_levels)
                ! pressure at rho levels (Pa)

  REAL, INTENT(IN) ::                                                   &
     u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                ! U at pressure points
     v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                ! V at pressure points
     t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels), &
                ! temperature at theta levels
     q(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels), &
                ! specific humidity at theta levels
     qcl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                ! liquid water content at theta levels
     qcf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                ! frozen water content at theta levels
     q1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),&
                ! normalized excess of water
                ! (:,:,K) is located at theta level K-1
     frac_gauss(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                tke_levels)
                ! cloud fraction derived with Gaussian distribution
                ! function
                ! (:,:,K) is located at theta level K-1

! Additional variables for SCM diagnostics which are dummy in full UM
  INTEGER, INTENT(IN) ::                                                &
     nSCMDpkgs             ! No of SCM diagnostics packages
  LOGICAL, INTENT(IN) ::                                                &
     L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

  LOGICAL, INTENT(IN) ::                                                &
     l_mixing_ratio
                 !  flag for qsat_mix
                 !  .true. return qsat as a mixing ratio
                 !  .false. return qsat as a specific humidity

  REAL, INTENT(INOUT) ::                                                &
     zhpar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
               ! boundary layer height evaluated with Richardson Number

!  Declaration of BL diagnostics.
  TYPE (strnewbldiag), INTENT(INOUT) :: BL_diag

  REAL, INTENT(OUT) ::                                                  &
     frac(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          tke_levels),                                                  &
                ! cloud fraction including that by convection
                ! (:,:,K) is located at theta level K-1
     wb_ng(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           tke_levels)
                ! Non-gradint buoyancy flux due to the skewness
                ! (:,:,K) is located at theta level K-1

! local variables









  INTEGER :: i, j, k,                                                   &
     k_par(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                   ! level for start of parcel ascent
     ktpar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                   ! highest theta level below inversion (at ZHPAR)
     k_neut(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                   ! level of neutral parcel buoyancy
     ktinv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                   ! top level of inversion
     k_lcl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                   ! level of lifting condensation
     interp_inv,                                                        &
                   ! flag to interpolated inversion heights (1=yes)
     topbl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                   ! 1 => top of bl reached
                   ! 2 => max allowable height reached
  REAL ::                                                               &
     exner(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           shcu_levels),                                                &
                   ! sigma^cappa
     th(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        shcu_levels),                                                   &
                   ! potential temperature
     thl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         shcu_levels),                                                  &
                   ! liquid water potential tempeature
     tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        shcu_levels),                                                   &
                   ! liquid water tempeature
     qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        shcu_levels),                                                   &
                   ! total water spec humidity
     thvl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          shcu_levels),                                                 &
                   ! virtual thl
     THv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         shcu_levels),                                                  &
                   ! virtual th
     thv_par(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             shcu_levels),                                              &
                   ! virtual th for parcel
     qc_par(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            shcu_levels),                                               &
                   ! parcel liquid water
     dthvdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            shcu_levels),                                               &
                   ! gradient of THV at rho levels
     dthvdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             shcu_levels)
                   ! gradient of THV at theta levels

  REAL ::                                                               &
     thl_par(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),      &
                   ! parcel thl
     qw_par(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                   ! parcel qw
     sl_par(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                   ! parcel static energy
     th_ref(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                   ! reference theta for parcel ascent
     th_par_kp1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),   &
                   ! parcel theta at level below
     p_lcl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                   ! pressure of LCL
     thv_pert(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
                   ! threshold for parcel thv
     z_lcl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                   ! height of LCL
     zh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
                   ! boundary layer depth (from RI)
     zhpar_old(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),    &
                   ! height of cloud-top on previous timestep
     zhpar_max(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),    &
                   ! Maximum allowed height for cloud-top
                   ! (to limit growth rate of boundary layer)
     w_star(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                   ! sub-cloud layer velocity scale (m/s)
     cape(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                   ! convective available potential energy (m2/s2)
     dbdz_inv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
                   ! buoyancy gradient across inversion
     dz_inv_cu(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),    &
                   ! inversion thickness above Cu
     qsat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                   ! saturated water mixing ratio
     t_ref(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                   ! reference temperature

  REAL ::                                                               &
     virt_factor,                                                       &
                   ! Vfac = 1+0.61qv - qcl - qcf
     z_surf,                                                            &
                   ! height of surface layer
     w_s,                                                               &
                   ! velocity scale
     thv_sd,                                                            &
                   ! standard deviation of thv in surface layer
     dqsatdt,qsatfac,                                                   &
                   ! saturation coefficients in buoyancy parameters
     qc_env,                                                            &
                   ! environment liquid water
     vap_press,                                                         &
                   ! Vapour pressure.
     t_lcl,                                                             &
                   ! temperature of LCL
     th_par,                                                            &
                   ! theta of parcel
     t_par,                                                             &
                   ! temperature of parcel
     dpar_bydz,                                                         &
                   ! parcel thv gradient
     denv_bydz,                                                         &
                   ! environement thv gradient
     gamma_fa,                                                          &
                   ! free atmospheric lapse rate
     gamma_cld,                                                         &
                   ! cloud layer lapse rate
     wb_scale,                                                          &
                   ! buoyancy flux scaling (m2/s3)
     w_cld,                                                             &
                   ! cloud layer velocity scale (m/s)
     z_cld,                                                             &
                   ! cloud layer depth (m)
     dz_inv_cu_rec,                                                     &
                   ! reconstructed inversion thickness above Cu
     vscalsq_incld,                                                     &
                   ! incloud squared velocity scale
     m_base,                                                            &
                   ! cloud base mass flux (m/s)
     z_pr,ze_pr,                                                        &
                   ! scaled height
     zpr_top,                                                           &
                   ! inversion top in scaled coordinate
     f_ng,                                                              &
                   ! non-gradient shape functions
     fnn,                                                               &
                   ! entrainment factor gN
     z0,z1,z2,z3,                                                       &
                   ! heights for polynomial interpolation
     d0,d1,d2,d3,                                                       &
                   ! values for polynomial interpolation
     a2,a3,xi,                                                          &
                   ! work variables for polynomial interpolation
     a_poly,b_poly,c_poly,                                              &
                   ! coefficients in polynomial interpolation
     ri,                                                                &
                   ! Richardson number
     grid_int,                                                          &
                   ! THV integral over inversion
     zhdisc,                                                            &
                   ! height of subgrid interpolated inversion
     weight1, weight2, weight3,                                         &
                   ! weight factors in interpolating
     lrcp_c,                                                            &
                   ! Latent heat over heat capacity
     l_heat,                                                            &
                   ! Latent heat
     frcu
                   ! cloud fraction due to convection

  LOGICAL ::                                                            &
     topinv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                   ! indicates top of inversion being reached
     topprof(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),      &
                   ! indicates top of ascent being reached
     above_lcl
                   ! indicates being above the LCL

  REAL, PARAMETER ::                                                    &
     a_parcel=0.2,                                                      &
     b_parcel=3.26,                                                     &
     max_t_grad=1.0e-3,                                                 &
     ric=0.25

! Derived local parameters.

  REAL, PARAMETER ::                                                    &
     lcrcp=lc/cp,                                                       &
                              ! Evaporation-to-dT conversion factor.
     ls=lf+lc,                                                          &
                              ! Latent heat of sublimation.
     lsrcp=ls/cp,                                                       &
                              ! Sublimation-to-dT conversion factor.
     grcp=g/cp

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('MYM_SHCU_BUOY',zhook_in,zhook_handle)
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      zhpar_old(i, j) = zhpar(i, j)
      ! Limit boundary layer growth rate to 0.14 m/s
      ! (approx 500m/hour)
      zhpar_max(i,j) = MIN( z_tq(i, j, shcu_levels-1),                  &
                            zhpar_old(i, j)+timestep*0.14 )
      zh(i, j) = 0.0
    END DO
  END DO
  DO k = 1, shcu_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ! initialise cumulus cloud fraction to zero
        exner(i, j, k) = (p_theta_levels(i, j, k) / pref) ** kappa
        th(i, j, k) = t(i, j, k) / exner(i, j, k)
        thl(i, j, k) = th(i, j, k)                                      &
             - (lcrcp*qcl(i,j,k) + lsrcp*qcf(i,j,k)) / exner(i, j, k)
        qw(i,j,k) = q(i,j,k) + qcl(i,j,k) + qcf(i,j,k)
        thvl(i,j,k)= thl(i,j,k) * ( 1. + c_virtual*qw(i,j,k) )
        virt_factor = 1. + c_virtual*q(i,j,k)                           &
                                - qcl(i,j,k) - qcf(i,j,k)
        THv(i,j,k) = th(i,j,k) * virt_factor
        thv_par(i,j,k) = THv(i,j,k)   ! default for stable bls
        wb_ng(i,j,k) = 0.0
        frac(i,j,k) = frac_gauss(i,j,k)
        tl(i,j,k) = t(i,j,k) - lcrcp*qcl(i,j,k) - lsrcp*qcf(i,j,k)
      END DO
    END DO
  END DO

  DO k = 2, shcu_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dthvdz(i,j,k) = THv(i,j,k) - THv(i,j,k-1)
      END DO
    END DO
  END DO

  DO k = 3, shcu_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        weight1 = z_uv(i,j,k) - z_uv(i,j,k-1)
        weight2 = z_tq(i,j,k-1)- z_uv(i,j,k-1)
        weight3 = z_uv(i,j,k) - z_tq(i,j,k-1)
        dthvdzm(i,j,k) = (weight2 * dthvdz(i,j,k)                       &
                         + weight3 * dthvdz(i,j,k-1)) / weight1
      END DO
    END DO
  END DO

  k = 2
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dthvdzm(i,j,k) = dthvdz(i,j,k)
    END DO
  END DO

  DO k = 2, shcu_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        ri = (u_p(i,j,k)-u_p(i,j,k-1))**2                               &
                         +(v_p(i,j,k)-v_p(i,j,k-1))**2
        ri = (g*(z_uv(i,j,k)-z_uv(i,j,k-1))                             &
                 *dthvdzm(i,j,k)/THv(i,j,k)) / MAX( 1.e-14, ri )
        IF ( ri > ric .AND. zh(i,j) == 0.0 ) THEN
          zh(i,j)=z_uv(i,j,k)
        END IF
        qc_par(i,j,k) = 0.0
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
! 1. Set up parcel
!-----------------------------------------------------------------------
! Start parcel ascent from grid-level above top of surface layer, taken
! to be at a height, z_surf, given by 0.1*ZH
!-----------------------------------------------------------------------
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      k_par(i,j) = 1
      zhpar(i,j) = zh(i,j)  ! initialise to bl depth (from RI)
      k_lcl(i,j) = 1
      IF (fb_surf(i,j) >= 0.0) THEN
        z_surf = 0.1 * zh(i,j)
        DO WHILE ( z_uv(i,j,k_par(i,j)) < z_surf .AND.                  &
                    ! not reached Z_SURF
                    thvl(i,j,k_par(i,j)+1) <= thvl(i,j,k_par(i,j)) )
                    ! not reached inversion
          k_par(i,j) = k_par(i,j) + 1
        END DO
        w_s = ( fb_surf(i,j)*zh(i,j) + ustar(i,j)**3 )**(1.0/3.0)
        thv_sd = 1.93 * fb_surf(i,j) * THv(i,j,k_par(i,j))              &
                      / ( g * w_s )
        thl_par(i,j) = thl(i,j,k_par(i,j))
        qw_par(i,j)  = qw(i,j,k_par(i,j))
        sl_par(i,j) = tl(i,j,k_par(i,j))                                &
                            + grcp * z_tq(i,j,k_par(i,j))
!-----------------------------------------------------------------------
! Calculate temperature and pressure of lifting condensation level
! using approximations from Bolton (1980)
!-----------------------------------------------------------------------
        vap_press = q(i,j,k_par(i,j)) *                                 &
              p_theta_levels(i,j,k_par(i,j)) / ( 100.0*repsilon )
        IF (vap_press >= 0.0) THEN
          t_lcl = 55.0 + 2840.0 / ( 3.5*ALOG(t(i,j,k_par(i,j)))         &
                     - ALOG(vap_press) - 4.805 )
          p_lcl(i,j) =  p_theta_levels(i,j,k_par(i,j)) *                &
               ( t_lcl / t(i,j,k_par(i,j)) )**(recip_kappa)
        ELSE
          p_lcl(i,j) = pstar(i,j)
        END IF
         ! K_LCL is model level BELOW the lifting condensation level
        k_lcl(i,j) = 1
        DO k = 2, shcu_levels
          IF (p_rho_levels(i,j,k) > p_lcl(i,j)) THEN
            k_lcl(i,j) = k - 1
          END IF
        END DO
        z_lcl(i,j) = z_uv(i,j,k_lcl(i,j)+1)                             &
              + ( z_uv(i,j,k_lcl(i,j))-z_uv(i,j,k_lcl(i,j)+1) )         &
              * ( p_rho_levels(i,j,k_lcl(i,j)+1) - p_lcl(i,j))          &
              / ( p_rho_levels(i,j,k_lcl(i,j)+1)                        &
                               - p_rho_levels(i,j,k_lcl(i,j)) )
        z_lcl(i,j) = MAX( z_uv(i,j,1), z_lcl(i,j) )
!-----------------------------------------------------------------------
! Threshold on parcel buoyancy for ascent, THV_PERT, is related to
! standard deviation of thv in surface layer
!-----------------------------------------------------------------------
        thv_pert(i,j)= MAX( a_parcel,                                   &
              MIN( max_t_grad*zh(i,j), b_parcel*thv_sd ) )

        th_ref(i,j) = thl_par(i,j)
        th_par_kp1(i,j) = thl_par(i,j)
      ELSE
        ! dummy
        th_ref(i,j) = thl(i,j,1)
        z_lcl(i,j) = z_uv(i, j, 1)
      END IF   ! test on unstable
    END DO
  END DO
!-----------------------------------------------------------------------
! 2  Parcel ascent:
!-----------------------------------------------------------------------
! Lift parcel conserving its THL and QW.
! Calculate parcel QC by linearising q_sat about the parcel's
! temperature extrapolated up to the next grid-level

  DO k = 1, shcu_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        t_ref(i,j) = th_ref(i,j)*exner(i,j,k)
      END DO
    END DO

    !DEPENDS ON: qsat_mix
    CALL qsat_mix(qsat,t_ref,p_theta_levels(1,1,k),                     &
                  tdims%i_end*tdims%j_end,l_mixing_ratio)

    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (fb_surf(i,j) > 0.0) THEN
          IF (t_ref(i,j) > tm) THEN
            lrcp_c = lcrcp
            l_heat = lc
          ELSE
            lrcp_c = lsrcp
            l_heat = ls
          END IF

          dqsatdt = repsilon * l_heat * qsat(i,j)/(r*t_ref(i,j)**2)
          qsatfac = 1./(1.+(lrcp_c)*dqsatdt)
          qc_par(i,j,k)  = MAX( 0.0,                                    &
                 qsatfac*( qw_par(i,j) - qsat(i,j)                      &
                - (thl_par(i,j)-th_ref(i,j))                            &
                     *exner(i,j,k)*dqsatdt ) )
          qc_env  = MAX( 0.0, qsatfac*( qw(i,j,k) - qsat(i,j)           &
                - (tl(i,j,k)-t_ref(i,j)) *dqsatdt ) )
          qc_par(i,j,k)  = qc_par(i,j,k) + qcl(i,j,k) + qcf(i,j,k)      &
                               - qc_env
          t_par = sl_par(i,j) - grcp * z_tq(i,j,k)                      &
                    + lrcp_c * qc_par(i,j,k)
          ! recalculate if signs of T_REF and T_PAR are different
          IF (t_ref(i,j) <= tm .AND. t_par > tm) THEN
            lrcp_c = lcrcp
            qsatfac = 1./(1.+(lrcp_c)*dqsatdt)
            qc_par(i,j,k)  = MAX( 0.0,                                  &
                   qsatfac*( qw_par(i,j) - qsat(i,j)                    &
                  - (sl_par(i,j)-grcp*z_tq(i,j,k)-t_ref(i,j))           &
                       *dqsatdt ) )
            qc_par(i,j,k)  = qc_par(i,j,k) + qcl(i,j,k) + qcf(i,j,k)    &
                                 - qc_env
            t_par = sl_par(i,j) - grcp * z_tq(i,j,k)                    &
                      + lrcp_c * qc_par(i,j,k)
          END IF
          th_par = t_par / exner(i,j,k)
          thv_par(i,j,k) = th_par *                                     &
                           (1.+c_virtual*qw_par(i,j)                    &
                             -(1.+c_virtual)*qc_par(i,j,k))
          IF (k > 1 .AND. k < shcu_levels - 1) THEN
          ! extrapolate reference TH gradient up to next grid-level
            z_pr      = (z_tq(i,j,k+1)-z_tq(i,j,k))                     &
                                 /(z_tq(i,j,k)-z_tq(i,j,k-1))
            th_ref(i,j) = th_par*(1.0+z_pr)                             &
                                 - th_par_kp1(i,j)*z_pr
            th_par_kp1(i,j) = th_par
          END IF
        END IF   ! test on unstable
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
! 3 Identify layer boundaries
!-----------------------------------------------------------------------
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      topbl(i,j) = 0
      topprof(i,j) = .FALSE.
      topinv(i,j)= .FALSE.
      ktpar(i,j) = 1
      k_neut(i,j) = 1
      ktinv(i,j) = 1
      dbdz_inv(i,j) = 0.003
                    ! start with a weak minimum inversion lapse rate
                    ! (~1.e-4 s^-2, converted from K/m to s^-2 later)
    END DO
  END DO

  DO k = 2, shcu_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        IF (fb_surf(i,j) > 0.0) THEN
          !------------------------------------------------------------
          ! Set flag to true when level BELOW is above the lcl
          ! and above LCL transition zone
          !------------------------------------------------------------
          above_lcl = k-1 > k_lcl(i,j) + 1                              &
                           .AND. z_tq(i,j,k-1) > 1.1*z_lcl(i,j)
          !-------------------------------------------------------------
          ! Calculate vertical gradients in parcel and environment THV
          !-------------------------------------------------------------
          dpar_bydz = (thv_par(i,j,k) - thv_par(i,j,k-1)) /             &
                (z_tq(i,j,k) - z_tq(i,j,k-1))
          denv_bydz = (THv(i,j,k) - THv(i,j,k-1)) /                     &
                (z_tq(i,j,k) - z_tq(i,j,k-1))
          !-------------------------------------------------------------
          ! Find top of inversion - where parcel has minimum buoyancy
          !-------------------------------------------------------------
          IF ( topbl(i,j) > 0 .AND. .NOT. topinv(i,j) ) THEN
            dbdz_inv(i,j) = MAX( dbdz_inv(i,j), denv_bydz )
            IF ( k-1 > ktpar(i,j)+2 .AND. (                             &
                 ! Inversion at least two grid-levels thick
                           denv_bydz <= dpar_bydz .OR.                  &
                 ! => at a parcel buoyancy minimum
                  z_uv(i,j,k) > zhpar(i,j)+MIN(1000., 0.5*zhpar(i,j))   &
                 )) THEN
                 ! restrict inversion thickness < 1/2 bl depth and 1km
              topinv(i,j) = .TRUE.
              ktinv(i,j) = k-1
            END IF
          END IF
          !-------------------------------------------------------------
          ! Find base of inversion - where parcel has maximum buoyancy
          !                          or is negatively buoyant
          !-------------------------------------------------------------
          IF ( .NOT. topprof(i,j) .AND. k > k_par(i,j) .AND.            &
               ((thv_par(i,j,k)-THv(i,j,k)                              &
                          <= - thv_pert(i,j)) .OR.                      &
                 k > shcu_levels - 1 )) THEN
            topprof(i,j) = .TRUE.
            k_neut(i,j) = k-1
          END IF

          IF ( topbl(i,j) == 0 .AND. k > k_par(i,j) .AND.               &
               (  ( thv_par(i,j,k)-THv(i,j,k)                           &
                            <= - thv_pert(i,j)) .OR.                    &
!                      plume non buoyant

                ( above_lcl .AND. (denv_bydz > 1.25*dpar_bydz) )        &

!                      or environmental virtual temperature gradient
!                      significantly larger than parcel gradient
!                      above lifting condensation level

                       )) THEN

            topbl(i,j) = 1
            ktpar(i,j) = k-1   ! marks most buoyant theta-level
                             ! (just below inversion)
            zhpar(i,j)    = z_uv(i,j,k)
            dbdz_inv(i,j) = MAX( dbdz_inv(i,j), denv_bydz )
          END IF

          IF ( topbl(i,j) == 0 .AND.                                    &
               (z_tq(i,j,k-1) >= zhpar_max(i,j)                         &
                                            .OR. k == shcu_levels)) THEN
!                      gone above maximum allowed height
            topbl(i,j) = 2
            ktpar(i,j) = k-2
            dbdz_inv(i,j) = MAX( dbdz_inv(i,j), denv_bydz )
          END IF
        END IF   ! test on unstable
      END DO
    END DO
  END DO

!-----------------------------------------------------------------------
! 3.1 Interpolate inversion base and top between grid-levels
!-----------------------------------------------------------------------
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      IF ( ktpar(i,j) > 1 ) THEN
          !-----------------------------------------------------
          ! parcel rose successfully
          !-----------------------------------------------------
        zhpar(i,j)    = z_uv(i,j,ktpar(i,j)+1)

        ! to determine if interpolation of the inversion is performed
        IF (topbl(i,j) == 2) THEN
            ! Stopped at max allowable height
          interp_inv= 0
          zhpar(i,j)  = zhpar_max(i,j)
          k = ktpar(i,j)
        ELSE
          interp_inv=1
          !-------------------------------------------------------
          ! First interpolate inversion base (max buoyancy excess)
          !-------------------------------------------------------
          !-----------------------------------------------------------
          ! interpolate height by fitting a cubic to 3 parcel excesses,
          ! at Z1 (top of cloud layer) and the two grid-levels above,
          ! and matching cloud layer gradient (D1-D0) at Z1.
          !-----------------------------------------------------------
          k = ktpar(i,j)+2
          z3=z_tq(i,j,k)  -z_tq(i,j,k-2)
          d3=thv_par(i,j,k)-THv(i,j,k)
          z2=z_tq(i,j,k-1)-z_tq(i,j,k-2)
          d2=thv_par(i,j,k-1)-THv(i,j,k-1)
          z1=z_tq(i,j,k-2)-z_tq(i,j,k-2)
          d1=thv_par(i,j,k-2)-THv(i,j,k-2)
          z0=z_tq(i,j,k-3)-z_tq(i,j,k-2)
          d0=thv_par(i,j,k-3)-THv(i,j,k-3)
          c_poly = (d1-d0)/(z1-z0)
          a2= d2 - d1 - c_poly*z2
          a3= d3 - d1 - c_poly*z3
          b_poly = (a3-a2*z3**3/z2**3)/(z3*z3*(1.0-z3/z2))
          a_poly = (a2-b_poly*z2*z2)/z2**3

          xi=b_poly*b_poly-3.*a_poly*c_poly
          IF (a_poly /= 0.0 .AND. xi > 0.0) THEN
            ! ZHPAR is then the height where the above
            ! polynomial has zero gradient
            zhpar(i,j) = z_tq(i,j,k-2)-(b_poly+SQRT(xi))                &
                                           /(3.*a_poly)
            zhpar(i,j) = MAX( MIN( zhpar(i,j), z_tq(i,j,k) ),           &
                                          z_tq(i,j,k-2) )
            IF ( zhpar(i,j) > z_tq(i,j,ktpar(i,j)+1) ) THEN
              ktpar(i,j)=ktpar(i,j)+1
            END IF
          END IF
          k = ktpar(i,j)
          denv_bydz = (THv(i,j,k+1) - THv(i,j,k)) /                     &
                          (z_tq(i,j,k+1) - z_tq(i,j,k))
        END IF
        IF ( interp_inv == 1 ) THEN
          !-----------------------------------------------------
          ! Now interpolate inversion top
          !-----------------------------------------------------
          IF ( ktinv(i,j) > ktpar(i,j)+1 ) THEN
            k = ktinv(i,j)+1
            dpar_bydz = (thv_par(i,j,k) - thv_par(i,j,k-1)) /           &
                  (z_tq(i,j,k) - z_tq(i,j,k-1))
            denv_bydz = (THv(i,j,k) - THv(i,j,k-1)) /                   &
                  (z_tq(i,j,k) - z_tq(i,j,k-1))
            IF (denv_bydz < dpar_bydz) THEN
            !-----------------------------------------------------------
            ! interpolate height by fitting a parabola to parcel
            ! excesses and finding the height of its minimum
            !-----------------------------------------------------------
              z1=z_tq(i,j,k)
              d1=thv_par(i,j,k)-THv(i,j,k)
              z2=z_tq(i,j,k-1)
              d2=thv_par(i,j,k-1)-THv(i,j,k-1)
              z3=z_tq(i,j,k-2)
              d3=thv_par(i,j,k-2)-THv(i,j,k-2)
              xi=z2**2-z3**2
              b_poly=( d1-d3 - (d2-d3)*(z1**2-z3**2)/xi ) /             &
                    ( z1-z3 - (z2-z3)*(z1**2-z3**2)/xi )
              a_poly=(d2 - d3 - b_poly*(z2-z3) )/xi
            END IF
          END IF   ! inversion top grid-level 2 levels above parcel top
        END IF   ! interp_inv flag
      END IF   ! parcel rose
    END DO
  END DO
!-----------------------------------------------------------------------
! 4. Integrate parcel excess buoyancy
!-----------------------------------------------------------------------
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      cape(i,j) = 0.0
    END DO
  END DO
  DO k = 2, shcu_levels - 1
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        IF (k > k_lcl(i,j) .AND. k <=  k_neut(i,j)-1) THEN
          cape(i,j) = cape(i,j) + (thv_par(i,j,k) - THv(i,j,k))         &
                      * (z_uv(i,j,k+1)-z_uv(i,j,k)) / THv(i,j,k)
        END IF
      END DO
    END DO
  END DO
!-----------------------------------------------------------------------
! 6. Calculate non-gradient fluxes and velocity scales
!-----------------------------------------------------------------------
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dz_inv_cu(i,j) = 0.0
      frcu = 0.0
      IF (fb_surf(i,j) > 0.0) THEN
        w_star(i,j) = ( fb_surf(i,j)*zhpar(i,j) )**(1.0/3.0)
                                     ! dry bl scale
        dbdz_inv(i,j) = g*dbdz_inv(i,j)/THv(i,j,ktpar(i,j))
                                     ! convert to buoyancy units
        dz_inv_cu(i,j) = 0.2*zhpar(i,j)
                                     ! default for no CAPE
      END IF

      IF (cape(i,j) > 0.0 .AND. zhpar(i,j) - z_lcl(i,j) > 0.0) THEN
        k = k_lcl(i,j)
        ! calculate velocity scales
        w_star(i,j) = ( fb_surf(i,j)*z_lcl(i,j) )**(1.0/3.0)
        m_base    = 0.04*w_star(i,j)
        cape(i,j) = g * cape(i,j)
        w_cld     = ( m_base * cape(i,j) )**(1./3.)
        z_cld     = zhpar(i,j) - z_lcl(i,j)
        ! calculate fluxes at LCL
        wb_scale  = ( w_cld**3/z_cld ) * SQRT( m_base/w_cld )

        !----------------------------------------------------------
        ! Estimate inversion thickness.
        !----------------------------------------------------------
        vscalsq_incld = 2.0*cape(i,j)
        dz_inv_cu(i,j)  = SQRT( vscalsq_incld/dbdz_inv(i,j) )
        !
        ! If inversion is unresolved (less than 3 grid-levels thick)
        ! then use profile reconstruction
        !
        IF ( ktpar(i,j) <= shcu_levels - 4 ) THEN
          IF ( dz_inv_cu(i,j)                                           &
               < z_tq(i,j,ktpar(i,j)+3) - z_tq(i,j,ktpar(i,j)) ) THEN
          !
          ! First interpolate to find height of discontinuous inversion
          !
            k = ktpar(i,j)
            gamma_cld = (THv(i,j,k)-THv(i,j,k-1))                       &
                                     /(z_tq(i,j,k)-z_tq(i,j,k-1))
            IF (k-2 > k_lcl(i,j)) THEN
              gamma_cld = MIN( gamma_cld,                               &
                           ( THv(i,j,k-1)-THv(i,j,k-2) )                &
                              /(   z_tq(i,j,k-1)-  z_tq(i,j,k-2) ) )
            END IF
            gamma_cld = MAX(0.0, gamma_cld)
            gamma_fa = (THv(i,j,k+4)-THv(i,j,k+3))                      &
                                  /(z_tq(i,j,k+4)-z_tq(i,j,k+3))
            gamma_fa = MAX(0.0, gamma_fa)
            ! Integrate thv over the inversion grid-levels
            grid_int =  (THv(i,j,k+1)-THv(i,j,k))                       &
                                 *(z_uv(i,j,k+2)-z_uv(i,j,k+1))         &
                   + (THv(i,j,k+2)-THv(i,j,k))                          &
                                 *(z_uv(i,j,k+3)-z_uv(i,j,k+2))         &
                   + (THv(i,j,k+3)-THv(i,j,k))                          &
                                 *( z_tq(i,j,k+3)-z_uv(i,j,k+3))

            c_poly = (THv(i,j,k+3)-THv(i,j,k))                          &
                          *(z_tq(i,j,k+3)-z_tq(i,j,k))                  &
                       - 0.5*gamma_fa*(z_tq(i,j,k+3)-z_tq(i,j,k))**2    &
                       - grid_int
            b_poly = -(THv(i,j,k+3)-THv(i,j,k)                          &
                             -gamma_fa*(z_tq(i,j,k+3)-z_tq(i,j,k)))
            a_poly = 0.5*(gamma_cld-gamma_fa)
            xi     = b_poly*b_poly-4.*a_poly*c_poly

            IF (xi >= 0.0.AND.                                          &
                        ( a_poly  /= 0.0 .OR. b_poly  /= 0.0 )) THEN
              IF (a_poly == 0.0) THEN
                dz_inv_cu_rec = -c_poly/b_poly
              ELSE
                dz_inv_cu_rec = (-b_poly-SQRT(xi))/(2.*a_poly)
              END IF
              zhdisc  = z_tq(i,j,k)+dz_inv_cu_rec
              !
              ! Now calculate inversion stability given Dz=V^2/DB
              !
              c_poly = -vscalsq_incld*THv(i,j,k+1)/g
              b_poly = THv(i,j,k+3)-gamma_fa *(z_tq(i,j,k+3)-zhdisc)    &
                    -THv(i,j,k)  -gamma_cld*(zhdisc  -z_tq(i,j,k))
              a_poly = 0.5*(gamma_cld+gamma_fa)
              xi=b_poly*b_poly-4.*a_poly*c_poly

              IF (xi >= 0.0 .AND.                                       &
                       ( a_poly /= 0.0 .OR. b_poly /= 0.0 )) THEN
                IF (a_poly == 0.0) THEN
                  dz_inv_cu_rec = -c_poly/b_poly
                ELSE
                  dz_inv_cu_rec = (-b_poly+SQRT(xi))/(2.*a_poly)
                END IF
                dz_inv_cu_rec = MIN( dz_inv_cu_rec,                     &
                                 2.0*(zhdisc-z_tq(i,j,ktpar(i,j))) )
                IF (dz_inv_cu_rec <= dz_inv_cu(i,j)) THEN
                  dz_inv_cu(i,j) = dz_inv_cu_rec
                END IF
              END IF  ! interpolation for DZ_INV_CU successful
            END IF    ! interpolation for ZHDISC successful

          END IF  ! inversion not resolved
        END IF  ! if ktpar(i,j) <= shcu_levels - 4

        zpr_top     = 1. + MIN(1., dz_inv_cu(i,j)/z_cld )
        DO k = 1, shcu_levels-1
            ! Z_PR=0 at cloud-base, 1 at cloud-top
          z_pr = ( z_uv(i,j,k+1) - z_lcl(i,j) )/ z_cld
          IF (z_pr > 0.0) THEN

!   Non-gradient function for WB

            f_ng = 0.0
            IF ( z_pr <= 0.9 ) THEN
                  ! function with gradient=0 at z=0.9
                  !                    f=0,1 at z=0,0.9
              ze_pr = z_pr/0.9
              f_ng  = 0.5 * SQRT(ze_pr) * (3.0-ze_pr)
            ELSE IF (z_pr <= zpr_top) THEN
              ze_pr = (z_pr-0.9)/(zpr_top-0.9)  ! from 0 to 1
              f_ng  = 0.5 * (1.+COS(pi*ze_pr))
            END IF
            fnn = 0.5 * (1.0 + TANH(0.8 * (q1(i,j,k+1) + 0.5)))
            wb_ng(i,j,k+1) = MIN((1.-fnn)*3.7*f_ng*wb_scale, wb_ng_max)
          END IF   ! if Z_PR > 0

!   Cloud fraction enhancement and sigma_s calculation (for ql)
!   (on Z rather than ZE levels)

          z_pr = ( z_tq(i,j,k) - z_lcl(i,j) )/ z_cld
            ! Z_PR=0 at cloud-base, 1 at cloud-top

          IF (z_pr > 0.0) THEN
            f_ng = 0.0
            IF ( z_pr <= 0.9 ) THEN
              f_ng = 1.0+3.0*EXP(-5.0*z_pr)     ! =4 at cloud-base
            ELSE IF ( z_pr < zpr_top ) THEN
              ze_pr = (z_pr-0.9)/(zpr_top-0.9)  ! from 0 to 1
              f_ng = 0.5*(1.+COS(pi*ze_pr))
            END IF
            frcu = 0.5*f_ng*MIN(0.5,m_base/w_cld)
          END IF   ! Z_PR > 0
          frac(i,j,k+1)    = MAX( frac_gauss(i,j,k+1), frcu)
        END DO   ! loop over K
      END IF ! Test on CAPE
    END DO
  END DO

  DO k = shcu_levels + 1, tke_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        frac(i, j, k) = 0.0
        wb_ng(i, j, k) = 0.0
      END DO
    END DO
  END DO

  IF (BL_diag%l_wb_ng) THEN
    DO k = 2, shcu_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%wb_ng(i, j, k) = wb_ng(i, j, k)
        END DO
      END DO
    END DO
  END IF

  IF (lhook) CALL dr_hook('MYM_SHCU_BUOY',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE mym_shcu_buoy
