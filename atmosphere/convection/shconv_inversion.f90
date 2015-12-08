! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Values of fields at the inversion level
!

      SUBROUTINE shconv_inversion( n_sh, ntra, imethod_precip, l_tracer &
     &,                     mb,wstar_up,wstar_dn,zcld,zlcl              &
     &,                     wthetav_surf, dthetav_cb, dthetav_cld       &
     &,                     theta_below, theta_mid, theta_above         &
     &,                     r_below,r_mid, r_above                      &
     &,                     p_below, p_mid, p_above                     &
     &,                     p_inv_base, p_inv_top, dp_cld               &
     &,                     exner_below, exner_mid,exner_above          &
     &,                     exner_inv_base                              &
     &,                   tracer_below,tracer_mid                       &
     &,                   wqr_inv, precip_product_inv_base              &
     &,                   wtheta_inv_base,wq_inv,wql_inv,wthetav_inv    &
     &,                   wthetal_inv,wqt_inv_base,wtracer_inv)

!
! Purpose:
!   Calculate values of theta, q etc at the inversion level. These
!   values are then used to calculate the fluxes at the inversion.
!   The inversion base level is defined as a uv level in the scheme.
!
!   Called by Shallow_conv (5A version)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!

USE cv_diag_param_mod, ONLY:                                          &
   a_bolton, b_bolton, c_bolton, d_bolton

USE earth_constants_mod, ONLY: g
USE atmos_constants_mod, ONLY:                                        &
          cp, r, repsilon, pref, kappa, c_virtual, rv, recip_kappa

USE water_constants_mod, ONLY: lc

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!

      integer, intent(in) ::                                            &
     &   n_sh                                                           &
                         ! No. of shallow convection points
     &,  ntra                                                           &
                         ! Number of tracer fields
     &, imethod_precip   ! switch for precip calculation

      logical, intent(in) ::                                            &
     &  L_tracer         ! true - tracers present.

      real, intent(in) ::                                               &
     &  mb(n_sh)                                                        &
                               ! cloud base mass flux
     &, wstar_up(n_sh)                                                  &
                               ! wstar_up  (m/s)
     &, wstar_dn(n_sh)                                                  &
                               ! wstar_dn  (m/s)
     &, zcld(n_sh)                                                      &
                               ! convective cloud depth (m)
     &, zlcl(n_sh)                                                      &
                               ! lifting condensation level (m)
     &, wthetav_surf(n_sh)                                              &
                               ! wthetav flux at the surface
     &, dthetav_cb(n_sh)                                                &
                               ! dthetav across cloud base (k)
     &, dthetav_cld(n_sh)                                               &
                               !  dthetav across cloud
     &, theta_below(n_sh)                                               &
                               ! theta on level below inversion base (K)
     &, theta_mid(n_sh)                                                 &
                               ! theta on level above inversion base (K)
     &, theta_above(n_sh)                                               &
                               ! theta on level above inversion top (K)
     &, r_below(n_sh)                                                   &
                               ! r on level below inversion base (kg/kg)
     &, r_mid(n_sh)                                                     &
                               ! r on level above inversion base (kg/kg)
     &, r_above(n_sh)                                                   &
                               ! r on level above inversion top (kg/kg)
     &, p_below(n_sh)                                                   &
                         ! pressure on level below inversion base (N/m2)
     &, p_mid(n_sh)                                                     &
                         ! pressure on level above inversion base (N/m2)
     &, p_above(n_sh)                                                   &
                         ! pressure on level above inversion top (N/m2)
     &, p_inv_base(n_sh)                                                &
                               ! pressure at inversion base(N/m2)
     &, p_inv_top(n_sh)                                                 &
                               ! pressure at inversion top (N/m2)
     &, dp_cld(n_sh)                                                    &
                               !  dp across cloud
     &, exner_below(n_sh)                                               &
                               ! exner on level below inversion base
     &, exner_mid(n_sh)                                                 &
                               ! exner on level above inversion base
     &, exner_above(n_sh)                                               &
                               ! exner on level above inversion top
     &, exner_inv_base(n_sh)                                            &
                               ! exner at inversion base
     &, wqr_inv(n_sh)                                                   &
                               ! rain water flux (kg/kg m/s)
     &, precip_product_inv_base(n_sh)                                   &
                                       ! precip production from inversion
                                  ! (kg/m2/s) / density ?
     &, tracer_below(n_sh,ntra)                                         &
                                 ! tracer level below inversion (kg/kg)
     &, tracer_mid(n_sh,ntra)  ! tracer level above inversion (kg/kg)
!
! Arguments with intent INOUT:
!
!          NONE

!
! Arguments with intent OUT:
!
      real, intent(out) ::                                              &
     &  wthetal_inv(n_sh)                                               & 
                                  !  w'thetal' at inversion (Km/s)
     &, wtheta_inv_base(n_sh)                                           & 
                                       !  w'theta' at inversion (Km/s)
     &, wq_inv(n_sh)                                                    & 
                                  !  w'q'  at inversion (kg/kg m/s)
     &, wqt_inv_base(n_sh)                                              & 
                                       !  w'qt'  at inversion (kg/kg m/s)
     &, wql_inv(n_sh)                                                   & 
                                  !  w'ql at inversion (kg/kg m/s)
     &, wthetav_inv(n_sh)                                               & 
                                  !  w'thetav' at inversion (Km/s)
     &, wtracer_inv(n_sh,ntra)    !  w'tracer' at inversion (kgm/kg/s)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
      real                                                              &
     & lc_o_cp                 ! lc/cp

      real                                                              &
     &  dthetav_inv(n_sh)                                               &
                              ! thetav across inversion (k)
     &, thetav_mid(n_sh)                                                &
     &, thetav_below(n_sh)                                              &
     &, thetav_above(n_sh)                                              &
     &, t_mid(n_sh)                                                     &
     &, t_above(n_sh)                                                   &
     &, t_below(n_sh)                                                   &
     &, r_inv_base(n_sh)                                                &
                              ! water vapour at inversion (mixing
                              !  ratio, kg/kg)
     &, t_inv_base(n_sh)                                                &
                                   ! T at inversion (K)
     &, rsat_inv_base(n_sh)                                             &
                                   ! rsat at inversion (kg/kg)
     &, drsatdt_inv_base(n_sh)                                          &
                                   ! drsat/dT at inversion (kg/kg/K)
     &, t_lcl_below(n_sh)                                               &
                              ! T at LCL for level below inversion base
     &, t_lcl_mid(n_sh)                                                 &
                              ! T at LCL for level above inversion base
     &, t_lcl_above(n_sh)                                               &
                              ! T at LCL for level above inversion top
     &, t_lcl_inv(n_sh)                                                 &
                              ! T at LCL for  inversion base
     &, pstar_below(n_sh)                                               &
                              ! p at LCL for level below inversion base
     &, pstar_mid(n_sh)                                                 &
                              ! p at LCL for level above inversion base
     &, pstar_above(n_sh)                                               &
                              ! p at LCL for level above inversion top
     &, mb_new1(n_sh)                                                   &
                              ! possible modified mb
     &, root_new(n_sh)                                                  &
                              ! possible modified root(mb/wsc)
     &, dthetav_inv2(n_sh)                                              &
                              !  dthetav across inversion
     &, dthetav_invnew(n_sh)                                            &
                              !  dthetav across inversion
     &, r_rsat_inv_base(n_sh)                                           &
                              ! r-rsat at inversion
     &, dpstar_inv(n_sh)                                                &
                              ! dpstar across inversion
     &, dp_inv(n_sh)                                                    &
                              ! dp across inversion
     &, rsat_lcl(n_sh)                                                  &
                              ! rsat at lcl
     &, dthetav_dpstar_top(n_sh)                                        &
                                  ! dthetav/dp* at top of inversion
     &, dthetav_dpstar_base(n_sh)                                       &
                                  ! dthetav/dp* at base of inversion
     &, dthetav_dpstar_layer(n_sh)                                      &
                                   ! dthetav/dp* across inversion
     &, beta_inv(n_sh)                                                  &
                              ! beta for the inversion layer
     &, dp_lower(n_sh)                                                  &
                              ! dp across lower part of inversion
     &, dp_upper(n_sh)                                                  &
                              ! dp across upper part of inversion
     &, thetav_inv_base(n_sh)                                           &
                              ! thetav and base of inversion (k)
     &, thetav_inv_top(n_sh)                                            &
                              ! thetav and top of inversion (k)
     &, theta_inv_base(n_sh)                                            &
                              ! thetav and base of inversion (k)
     &, pstar_inv_base(n_sh)  ! pstar at inversion base

! temporary variables

      real                                                              &
     &  term_a                                                          &
     &, term_b                                                          &
     &, term_c                                                          &
     &, vap_press                                                       &
     &, div_term2, divisor_term                                         &
     &,r_min

      real                                                              &
     &  beta_cld

       parameter (beta_cld = 1.1)       ! could be 1.2

! Loop counters
!

      integer :: i, ktra
      
      Logical :: l_3lev_inv     ! true for 3 level inversion cal
                                ! current default is false.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
       
!-----------------------------------------------------------------------
! 1.0 Initialise arrays
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SHCONV_INVERSION',zhook_in,zhook_handle)

      l_3lev_inv = .false.

      r_min = 1.0e-8        ! should use user value to be consistent
                            ! Quick fix.
      Do i = 1,n_sh

        wtheta_inv_base(i) = 0.0
        wthetav_inv(i) = 0.0
        wthetal_inv(i) = 0.0
        wq_inv(i)     = 0.0
        wqt_inv_base(i)     = 0.0
        wql_inv(i)    = 0.0

      End do

      If (l_tracer) then
        Do ktra=1,ntra
          Do i = 1,n_sh
            wtracer_inv(i,ktra) = 0.0
          End Do
        End Do
      End If

! set up constants to be used in routine

      lc_o_cp = lc/cp

!-----------------------------------------------------------------------
! values at levels above and below inversion base and top
!-----------------------------------------------------------------------

      Do i = 1,n_sh
! For mixing ratio
        thetav_below(i) = theta_below(i)*                               &
     &                        (1.+r_below(i)/repsilon)/(1.+r_below(i))
        thetav_mid(i)   = theta_mid(i)*                                 &
     &                        (1.+r_mid(i)/repsilon)/(1.+r_mid(i))
        thetav_above(i) = theta_above(i)*                               &
     &                        (1.+r_above(i)/repsilon)/(1.+r_above(i))

        t_below(i) = theta_below(i)*exner_below(i)
        t_mid(i)   = theta_mid(i)  *exner_mid(i)
        t_above(i) = theta_above(i)*exner_above(i)

      End Do


! use modified mb for all inversion calculations

      Do i = 1,n_sh
         ! Problem with wthetav_surf = 0.0
         If (wthetav_surf(i) == 0) then
             mb_new1(i) = 0.2*wstar_dn(i)
         Else
             mb_new1(i) = 0.2*wthetav_surf(i)/                          &
     &           (max(dthetav_cb(i),wthetav_surf(i)/wstar_dn(i)))
         End If
         root_new(i) = sqrt(mb_new1(i)/wstar_up(i))
      End Do

!=======================================================================
! calculation of q and th at inversion
!=======================================================================

!-----------------------------------------------------------------------
! Calculate of r and theta at inversion
! Require dp*/dp to do this  dtheta/dp* and dr/dp* are ~ constant
! through inversion
!
! Need p* above and below inversion base and top
!
! Calculate temperature and pressure of lifting condensation level
!  using approximations from Bolton (1980)
!-----------------------------------------------------------------------

      If (l_3lev_inv) then
       Do i = 1,n_sh
        If (r_mid(i) >  0.0) then
          vap_press = 0.01*r_mid(i)*p_mid(i)/(repsilon+r_mid(i))
        Else
          vap_press = 0.01*r_min*p_mid(i)/(repsilon+r_min)
        End If
        t_lcl_mid(i) =  a_bolton + b_bolton / ( c_bolton*log(t_mid(i))  &
     &                                    - log(vap_press) - d_bolton )
        pstar_mid(i) = p_mid(i) *                                       &
                           ( t_lcl_mid(i) / t_mid(i) )**recip_kappa

        If (r_below(i) >  0.0) then
          vap_press = 0.01*r_below(i)*p_below(i)/(repsilon+r_below(i))
        Else
          vap_press = 0.01*r_min*p_below(i)/(repsilon+r_min)
        End If
        t_lcl_below(i) = a_bolton + b_bolton/( c_bolton*log(t_below(i)) &
     &                                    - log(vap_press) - d_bolton )
        pstar_below(i) = p_below(i) *                                   &
                         ( t_lcl_below(i) / t_below(i) )**recip_kappa

        If (r_above(i) >  0.0) then
          vap_press = 0.01*r_above(i)*p_above(i)/(repsilon+r_above(i))
        Else
          vap_press = 0.01*r_min*p_above(i)/(repsilon+r_min)
        End If
        t_lcl_above(i) = a_bolton + b_bolton/( c_bolton*log(t_above(i)) &
     &                                    - log(vap_press) - d_bolton )
        pstar_above(i) = p_above(i) *                                   &
                         ( t_lcl_above(i) / t_above(i) )**recip_kappa
       End do

! gradients of dthetav/dp* across inversion base and inversion top
!  and across the whole layer

       Do i = 1,n_sh
        dthetav_dpstar_top(i) = (thetav_above(i)-thetav_mid(i))/        &
     &                       (pstar_above(i)-pstar_mid(i))

        dthetav_dpstar_base(i) = (thetav_mid(i)-thetav_below(i))/       &
     &                       (pstar_mid(i)-pstar_below(i))

        dthetav_dpstar_layer(i) = (thetav_above(i)-thetav_below(i))/    &
     &                       (pstar_above(i)-pstar_below(i))

! pressure across inversion layer
        dp_inv(i) = p_inv_top(i)-p_inv_base(i)


! dthetav across the inversion layer

        dthetav_inv(i) = (0.16+0.55*zcld(i)/zlcl(i))                    &
     &               *wthetav_surf(i)/(0.2*mb(i))


! dp* across inversion layer

        dpstar_inv(i) = dthetav_inv(i)/dthetav_dpstar_layer(i)

        beta_inv(i)  = dpstar_inv(i)/dp_inv(i)

        dp_lower(i) = p_mid(i) - p_inv_base(i)
        dp_upper(i) = p_inv_top(i) - p_mid(i)
       End do

! calculate thetav at inversion top and inversion base using interpolation
! from value at inversion middle.

       Do i = 1,n_sh

        thetav_inv_base(i) = thetav_mid(i)                              &
     &          - beta_inv(i)*dthetav_dpstar_base(i)*dp_lower(i)

        thetav_inv_top(i) = thetav_mid(i)                               &
     &          + beta_inv(i)*dthetav_dpstar_top(i)*dp_upper(i)

        pstar_inv_base(i) = pstar_mid(i) - beta_inv(i)*dp_lower(i)

        dthetav_inv2(i) = thetav_inv_top(i) - thetav_inv_base(i)
        If (pstar_inv_base(i) < 0.0) then
           write(6,*) ' Problem pstar_inv_base ',i,pstar_inv_base(i),&
   &thetav_inv_base(i),thetav_inv_top(i),beta_inv(i),dp_lower(i),dp_upper(i),&
   &dthetav_dpstar_base(i),dthetav_dpstar_top(i)
     write(6,*) ' above ',t_lcl_above(i),thetav_above(i),pstar_above(i),&
    &t_above(i),p_above(i)

     write(6,*) ' below ',t_lcl_below(i),thetav_below(i),pstar_below(i),&
    &t_below(i),p_below(i)
     write(6,*) ' mid ',t_lcl_mid(i),thetav_mid(i),pstar_mid(i),&
    &t_mid(i),p_mid(i)
        End If

! There are points with problems 
!  dthetav_inv2   -ve
!  dthetav_inv2 >> dthetav_inv
! For these points need to use another method ?


       End do

      Else   ! 2 level inversion

        Do i = 1,n_sh


          If (r_mid(i) >  0.0) then
            vap_press = 0.01*r_mid(i)*p_mid(i)/(repsilon+r_mid(i))
          Else
            vap_press = 0.01*r_min*p_mid(i)/(repsilon+r_min)
          End If
          t_lcl_mid(i) = a_bolton + b_bolton /( c_bolton*log(t_mid(i)) &
     &                                    - log(vap_press) - d_bolton )
          pstar_mid(i) = p_mid(i) *                                    &
                           ( t_lcl_mid(i) / t_mid(i) )**recip_kappa

        End do

        Do i = 1,n_sh
          dp_inv(i) = p_mid(i)-p_below(i)
          dthetav_invnew(i) = (0.16+0.55*zcld(i)/zlcl(i))               &
     &               *wthetav_surf(i)/(0.2*mb(i))

! restriction on dthetav_cld to avoid problems with strange profiles

          If (dthetav_cld(i).gt.0.01) then
            dpstar_inv(i)     = beta_cld*dthetav_invnew(i)*dp_cld(i)    &
     &                              /dthetav_cld(i)

          Else
            dpstar_inv(i)     = beta_cld*dthetav_invnew(i)*dp_cld(i)    &
     &                              /0.01
          End If

          pstar_inv_base(i)  = pstar_mid(i) - dpstar_inv(i)
          thetav_inv_base(i) = thetav_mid(i) - dthetav_invnew(i)
!
! value of thetav across inversion
!
          dthetav_inv2(i) = thetav_mid(i) - thetav_below(i)

        End Do
      End If  ! inversion cal type


!---------------------------------------------------------------------
! Calculate of r-rsat at inversion
!---------------------------------------------------------------------

      Do i = 1,n_sh
          t_lcl_inv(i) = thetav_inv_base(i)*                            &
     &                           ((pstar_inv_base(i)/pref)**kappa)
      End Do

! DEPENDS ON: qsat_mix
      Call qsat_mix(rsat_lcl,t_lcl_inv,pstar_inv_base,n_sh,.true.)


      Do i = 1,n_sh
        r_inv_base(i) = rsat_lcl(i)/                                    &
     &             (1.+c_virtual*(lc/(rv*t_lcl_inv(i)))*rsat_lcl(i))

        theta_inv_base(i) = thetav_inv_base(i)                          &
     &                                /(1.+c_virtual*r_inv_base(i))
        t_inv_base(i) = theta_inv_base(i)*exner_inv_base(i)
      End Do


! DEPENDS ON: qsat_mix
      call qsat_mix(rsat_inv_base,t_inv_base,p_inv_base,n_sh,.true.)

      Do i = 1,n_sh

        drsatdt_inv_base(i) = Lc *rsat_inv_base(i)/                     &
     &                           (Rv * t_inv_base(i)*t_inv_base(i))
        r_rsat_inv_base(i) = r_inv_base(i) - rsat_inv_base(i)

      End do


!-----------------------------------------------------------------------
!  Calculate fluxes at inversion uses values calculated above
!-----------------------------------------------------------------------

!----------------------------------------------------------------------
! Solution for fluxes at inversion
!----------------------------------------------------------------------


      Do i = 1,n_sh

! 0.19 from LES fits could be regarded as something which could be
! varied (data left by Alan imply values 0.14-0.19)
! precipitation adds extra terms to A and B

       If (imethod_precip == 1 ) then
        term_a = -0.19*mb_new1(i)*dthetav_inv2(i) +                     &
     &   (lc_o_cp/exner_inv_base(i)                                     &
     &            -0.61*theta_inv_base(i))*precip_product_inv_base(i)

        term_b = thetav_inv_base(i)*root_new(i)*(wstar_up(i)**3)/       &
     &                       (g*zcld(i))                                &
     &                                 +theta_inv_base(i)*wqr_inv(i)
       Else
        term_a = -0.19*mb_new1(i)*dthetav_inv2(i)

        term_b = thetav_inv_base(i)*root_new(i)*(wstar_up(i)**3)/       &
     &                       (g*zcld(i))
       End If

        term_c = 0.35*mb_new1(i)*r_rsat_inv_base(i)

! Corrected version
        div_term2    = lc_o_cp/exner_inv_base(i)                        &
     &                               -(1.+c_virtual)*theta_inv_base(i)
        divisor_term = (1.+c_virtual*t_inv_base(i)*drsatdt_inv_base(i))

!
! viscosity based value
!
        wthetav_inv(i) = term_a
!
! Liquid water flux correct for mixing ratio
!
        wql_inv(i) =(term_b-term_a)/div_term2

        wqt_inv_base(i) = ((1.+lc_o_cp*drsatdt_inv_base(i))*            &
     &    wql_inv(i)/0.85-term_c+                                       &
     &    term_a*exner_inv_base(i)*drsatdt_inv_base(i))/divisor_term

        wthetal_inv(i) = term_a                                         &
     &                     -c_virtual*theta_inv_base(i)*wqt_inv_base(i)

!
! above implies
!
        wq_inv(i)=wqt_inv_base(i) - wql_inv(i)

        wtheta_inv_base(i) = wthetal_inv(i)                             &
     &                        +lc_o_cp*wql_inv(i)/exner_inv_base(i)

      End do

! Note at base of inversion expect wtheta & wthetal to be negative and
! wq, wqt to be positive

!-----------------------------------------------------------------------
! Tracer flux at inversion - is this accurate (Do something slightly
! different for dtheta and dq across inversion ?)
!
      If (l_tracer) then
       Do ktra = 1,ntra
         Do i = 1,n_sh
           wtracer_inv(i,ktra) = -0.19*mb_new1(i)*                      &
     &                     (tracer_mid(i,ktra) - tracer_below(i,ktra))
         End Do
       End Do
      End If

!-----------------------------------------------------------------------
!    End Subroutine
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SHCONV_INVERSION',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE shconv_inversion
