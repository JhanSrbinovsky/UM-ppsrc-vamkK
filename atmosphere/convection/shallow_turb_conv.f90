! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Shallow convection scheme - turbulence based version
!

      SUBROUTINE SHALLOW_TURB_CONV(nlev,n_sh,ntra,trlev                 &
     &,                      ntml,ntpar,ishall_precip                   &
     &,                      L_calc_dxek,L_q_interact,L_tracer          &
     &,                      bland                                      &
     &,                      timestep                                   &
     &,                      delthvu,ql_ad,uw0,vw0                      &
     &,                      wstar_dn,wthvs,zlcl,zlcl_uv,ztop_uv        &
     &,                      pstar,p_layer_centres,p_layer_boundaries   &
     &,                      exner_layer_centres                        &
     &,                      exner_layer_boundaries                     &
     &,                      z_theta,z_rho, r_rho, r_theta              &
     &,                      rho, rho_theta, r2rho, r2rho_th            &
     &,                      dr_across_rh,dr_across_th                  &
     &,                      q_mix, theta, u,v                          &
     &,                      qse                                        &
             ! Inout fields (PC2 ones not used)
     &,                      bulk_cf,cf_frozen,cf_liquid                &
     &,                      qcf,qcl, tracer                            &
             ! Output fields
     &,                      cape_out,cclwp,ccw,cca                     &
     &,                      dbcfbydt,dcffbydt,dcflbydt                 &
     &,                      dqbydt,dqcfbydt,dqclbydt                   &
     &,                      dthbydt,dubydt,dvbydt,dtrabydt             &
     &,                      iccb,icct,lcca,lcbase,lctop                &
     &,                      rain,snow, up_flux                         &
     &,                      uw_shall,vw_shall                          &
     &,                      wqt,wthetal,wthetav,wql                    &
     &,                      wstar_up, mb, mb_new                       &
     &,                      tcw, cca_2d)


!
! Purpose:
!   Shallow convection scheme - works on points diagnosed as shallow in
!   subroutine CONV_DIAG.
! 
!   Interaction with PC2.
!   ---------------------
!   This version of the shallow scheme returns a non-zero CCA
!   array for shallow points and a ccw array. All PC2 increments are
!   zero from this scheme. The CCA and CCW are used later in IMP_CTL2
!   to overwrite the large scale cloud liquid fraction with cca and 
!   qcl with cca*ccw for where shallow clouds are diagnosed (ie cca non
!   zero).
! 
!
!   Called by GLUE_CONV4A.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN77 + common extensions
!  This code is written to UMDP3 v6 programming standards
!
USE water_constants_mod, ONLY: lc
USE conversions_mod, ONLY: zerodegc

USE atmos_constants_mod, ONLY:                                          &
   cp, r, repsilon, c_virtual, rv, recip_epsilon

USE cv_derived_constants_mod, ONLY:                                     &
   lcrcp, gamma_dry, ra2, cv  

USE cv_run_mod, ONLY:                                                   &
          i_convection_vn,i_convection_vn_6a, l_mom

USE cv_stash_flg_mod, ONLY:                                             &
   flg_uw_shall, flg_vw_shall              

USE earth_constants_mod, ONLY: g, earth_radius
      
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------
!
! Arguments with intent IN:
!

      Integer, intent(in) ::                                            &
     &  nlev                                                            &
                   ! No. of model layers

     &, n_sh                                                            &
                   ! No. of shallow convection points

     &, ntra                                                            &
                   ! No. of tracer fields

     &, trlev      ! No. of model levels on which tracers are included

      Integer, intent(in) ::                                            &
     &  ntml(n_sh)                                                      &
                      ! Top level of surface mixed layer defined
                      ! relative to theta,q grid

     &, ntpar(n_sh)                                                     &   
                      ! Top level of initial parcel ascent in BL
                      ! scheme defined relative to theta,q grid

     &, ishall_precip
                      ! 0 no shallow precip, 1 scheme with shall precip
     
      Logical, intent(in) ::                                            &
     &  L_calc_dxek                                                     &
                      ! Switch for calculation of condensate increment
                      ! (PC2)
     &, L_q_interact                                                    &
                      ! Switch allows overwriting parcel variables when
                      ! calculating condensate incr. (PC2)

     &, L_tracer      ! Switch for inclusion of tracers


! Not used at present but left in case required in future
      Logical, intent(in) :: bland(n_sh) ! Land/sea mask

      Real, intent(in) ::                                               &
     &  timestep           ! Model timestep (s)

      Real, intent(in) ::                                               &
     &  delthvu(n_sh)                                                   &
                           ! Integral of undilute parcel buoyancy over
                           ! convective cloud layer (Km)
     &, ql_ad(n_sh)                                                     &
                           ! adiabatic liquid water content at inversion
                           ! (kg/kg)
     &, uw0(n_sh)                                                       &
                           ! U-comp of surface stress (N/m2)

     &, vw0(n_sh)                                                       &
                           ! V-comp of surface stress (N/m2)

     &, wstar_dn(n_sh)                                                  &
                           ! Mixed layer convective velocity scale

     &, wthvs(n_sh)                                                     &
                           ! Surface flux of THV (Km/s)


     &, zlcl(n_sh)                                                      &
                           ! Lifting condensation level - accurate
                           ! height (m), not a model level height.

     &, zlcl_uv(n_sh)                                                   &
                           ! Lifting condensation level defined for
                           ! the uv grid (m)

     &, ztop_uv(n_sh)      ! Top of cloud layer defined for the uv
                           ! grid (m)

      Real, intent(in)    ::                                            &
     &  pstar(n_sh)                                                     &
                      ! Surface pressure (Pa)

     &, p_layer_centres(n_sh,0:nlev)                                    &
                                            ! Pressure (Pa)

     &, p_layer_boundaries(n_sh,0:nlev)                                 &
                                          ! Pressure at half level above
                                          ! p_layer_centres (Pa)
     &, exner_layer_centres(n_sh,0:nlev)                                &
                                            !Exner

     &, exner_layer_boundaries(n_sh,0:nlev) ! Exner at half level above
                                            ! exner_layer_centres
      Real, intent(in)    ::                                            &
     &  z_theta(n_sh,nlev)                                              &
                                     ! height of theta levels(m)
     &, z_rho(n_sh,nlev)                                                &
                                     ! height of rho levels (m)
     &, r_rho(n_sh,nlev)                                                &
                                     ! radius of rho levels (m)
     &, r_theta(n_sh,0:nlev)                                            &
                             ! radius of theta levels (+ surface) (m)
     &, rho(n_sh,nlev)                                                  &
                             ! density on rho levels kg/m3
     &, rho_theta(n_sh,nlev)                                            &
                             ! density on theta levels kg/m3
     &, r2rho(n_sh,nlev)                                                &
                                     ! r2*rho rho levels (kg/m)
     &, r2rho_th(n_sh,nlev)                                             &
                                     ! r2*rho theta levels (kg/m)
     &, dr_across_rh(n_sh,nlev)                                         &
                                     ! dr across rho levels  (m)
     &, dr_across_th(n_sh,nlev)      ! dr across theta levels (m)

      Real, intent(in)    ::                                            &
     &  q_mix(n_sh,nlev)                                                &
                           ! q as mixing ratio(kg/kg)


     &, theta(n_sh,nlev)                                                &
                           ! Model potential temperature (K)

     &, u(n_sh,nlev)                                                    &
                           ! Model U field (m/s)

     &, v(n_sh,nlev)                                                    &
                           ! Model V field (m/s)

     &, qse(n_sh,nlev)     ! Saturation mixing ratio of cloud
                           ! environment (kg/kg)

!
! Arguments with intent INOUT:
!
! Variables used only by PC2

      Real, intent(inout) ::                                            &
     &  bulk_cf(n_sh,nlev)                                              &
                             ! Bulk total cloud volume

     &, cf_frozen(n_sh,nlev)                                            &
                             ! Frozen water cloud volume

     &, cf_liquid(n_sh,nlev)                                            &
                             ! Liq water cloud volume

     &, qcf(n_sh,nlev)                                                  &
                             ! Ice condensate mix ratio (kg/kg)

     &, qcl(n_sh,nlev)       ! Liq condensate mix ratio (kg/kg)

! Variable only used if LTRACER

      Real, intent(inout)    :: tracer(n_sh,trlev,ntra)
                                   !Model tracer fields (kg/kg)

!
! Arguments with intent OUT:
!

      Real, intent(out) :: &
        cape_out(n_sh)     & ! Saved convective available potential energy
                             ! for diagnostic output (J/kg) NOT used in 5A
      , cclwp(n_sh)        & ! Condensed water path (kg/m2)
      , ccw(n_sh,nlev)     & ! Convective cloud liquid water on model
                             ! levels (g/kg)
      , cca(n_sh,nlev)       ! Convective cloud amount on model levels 
                             ! (fraction)

      Real, intent(out) ::                                              &
     &  dbcfbydt(n_sh,nlev)                                             &
                            ! Increments to total cld volume due to
                            ! convection(/s)

     &, dcffbydt(n_sh,nlev)                                             &
                            ! Increments to ice cloud volume due to
                            !  convection(/s)

     &, dcflbydt(n_sh,nlev)                                             &
                            ! Increments to liq cloud volume due to
                            ! convection(/s)

     &, dqbydt(n_sh,nlev)                                               &
                            ! Increments to q (water vap mixing ratio)
                            ! due to convection (kg/kg/s)

     &, dqcfbydt(n_sh,nlev)                                             &
                            ! Increments to ice condensate due to
                            ! convection (kg/kg/s)

     &, dqclbydt(n_sh,nlev)                                             &
                            ! Increments to liq condensate due to
                            ! convection (kg/kg/s)

     &, dthbydt(n_sh,nlev)                                              &
                            ! Increments to potential temp. due to
                            ! convection (K/s)

     &, dubydt(n_sh,nlev+1)                                             &
                            ! Increments to U due to CMT (m/s2)

     &, dvbydt(n_sh,nlev+1)                                             &
                            ! Increments to V due to CMT (m/s2)

     &, dtrabydt(n_sh,nlev,ntra) !Increment to tracer due to
                                 ! convection (kg/kg/s)

      Integer, intent(out) ::                                           &
     &  iccb(n_sh)                                                      &
                     ! Convective cloud base level (m)
     &, icct(n_sh)   ! Convective cloud top level (m)

      Real, intent(out) ::                                              &
     &  lcca(n_sh)   ! Lowest conv. cloud amt. (%)

      Integer, intent(out) ::                                           &
     &  lcbase(n_sh)                                                    &
                     ! Lowest conv. cloud base level (m)
     &, lctop(n_sh)  ! Lowest conv. cloud top level (m)

      Real, intent(out) ::                                              &
     &  rain(n_sh)                                                      &
                     ! Surface convective rainfall (kg/m2/s)
     &, snow(n_sh)   ! Surface convective snowfall (kg/m2/s)

! 5A scheme does not have explicit up and down draughts

      Real, intent(out) ::                                              &
     &  up_flux(n_sh,nlev)  ! mass flux (Pa/s)

      Real, intent(out) ::                                              &
     &  uw_shall(n_sh,nlev)                                             &
                            ! X-comp. of stress from shallow convection
                            !(kg/m/s2)
     &, vw_shall(n_sh,nlev)                                             &
                            ! Y-comp. of stress from shallow convection
                            !(kg/m/s2)
     &, wqt(n_sh,nlev)                                                  & 
                            ! w'qt' flux  (m/s kg/kg)
     &, wthetal(n_sh,nlev)                                              & 
                            ! w'thetal' flux  (Km/s)
     &, wthetav(n_sh,nlev)                                              & 
                            ! w'thetav' flux  (Km/s)
     &, wql(n_sh,nlev)      ! w'ql'   (kg m/s/kg)

      Real, intent(out) ::                                              &
     &  wstar_up(n_sh)                                                  &
                        ! cumulus layer convective velocity scale(m/s)
     &, mb(n_sh)                                                        &
                        ! Cloud base mass flux (m/s)
     &, mb_new(n_sh)    ! revised mb for incloud calculations (m/s)

      Real, intent(out) ::                                              &
     &  tcw(n_sh)                                                       &
                     ! Total condensed water(kg/m2/s)

     &, cca_2d(n_sh) ! 2D convective cloud amount (%)

!
!-----------------------------------------------------------------------
! Variables defined loCally
!-----------------------------------------------------------------------

      CHARACTER(LEN=*), Parameter ::  RoutineName = 'shallow_turb_conv'

      Real ::                                                           &
     &   wup(n_sh,nlev)                                                 &
                             ! ensemble vertical velocity (m/s)
     &,  wtheta(n_sh,nlev)                                              &
                             ! wtheta flux
     &,  wq(n_sh,nlev)                                                  &
                             ! wq flux
     &,  qlup(n_sh,nlev)     ! ensemble qlup   (kg/kg)

      Real ::                                                           &
     &   wtracer(n_sh,trlev,ntra) ! wtracer flux on uv levels (kgm/kg/s)
!
! Limit nlev loop to those levels actually required using ntpar
! diagnosed in conv_diag
!
      Integer :: ntpar_max         ! max ntpar value
      Integer :: ntml_max          ! max ntml value

!
! Local compressed arrays  - on cloud levels
!
      Real ::                                                           &
     &  theta_cld(n_sh,nlev)                                            &
                               ! theta
     &, q_cld(n_sh,nlev)                                                &
                               ! q - mixing ratio
     &, exner_cld(n_sh,nlev)                                            &
                               ! exner on theta levels in cloud
     &, exner_uvcld(n_sh,nlev)                                          &
                               ! exner pressure on uv levels in cloud
     &, p_cld(n_sh,nlev)                                                &
                               ! pressure on theta levels in cloud
     &, z_rho_cld(n_sh,nlev)                                            &
                               ! height of rho levels
     &, rho_cld(n_sh,nlev)                                              &
                               ! density on  rho levels
     &, tracer_cld(n_sh,trlev,ntra)  ! tracer in cloud
!
! values evaluated for cloud levels
!
      Real ::                                                           &
     &  eta(n_sh,nlev)                                                  &
                            ! non-dimension height for cloud levels
     &, eta_half(n_sh,nlev)                                             &
                            ! non-dimension height for cloud rho levels
     &, wup_cld(n_sh,nlev)                                              &
                               ! cumulus ensemble vertical velocity
     &, wup_h_cld(n_sh,nlev)                                            &
                               ! cumulus ensemble vertical velocity
     &, ql_up_cld(n_sh,nlev)                                            &
                               ! cumulus ensemble ql
     &, dwql_dz(n_sh,nlev)                                              &
                               ! dwql/dz all level
     &, qsat_cld(n_sh,nlev)                                             &
                               ! qsat in cloud
     &, tc_cld(n_sh,nlev)                                               &
                               ! temperature in degrees Celcius in cloud
     &, t_cld(n_sh,nlev)                                                &
                               ! temperature in degrees Celcius in cloud
     &, dqsatdt_cld(n_sh,nlev)                                          &
                               ! dqsat/dt
     &, mf_cld(n_sh,nlev)                                               &
                               ! mass flux in cld (theta levels)
     &, mf_h_cld(n_sh,nlev)                                             &
                               ! mass flux in cld (rho levels)
     &, dthetavdz(n_sh,nlev)                                            &
                               ! dthetav/dz
     &, kdtracerdz(n_sh,nlev,ntra)                                      &
                                    ! K*dtracer/dz in cloud
     &, wq_cld(n_sh,nlev)                                               &
                               ! wq in cloud
     &, wqt_cld(n_sh,nlev)                                              &
                               ! wqt in cloud
     &, wql_cld(n_sh,nlev)                                              &
                               ! wql in cloud with cb & inv values
     &, wtheta_cld(n_sh,nlev)                                           &
                               ! wtheta in cloud
     &, wthetal_cld(n_sh,nlev)                                          &
                               ! wthetal in cloud
     &, wthetav_cld(n_sh,nlev)                                          &
                               ! wthetav in cloud
     &, ztheta_cld(n_sh,nlev)                                           &
                               ! height of model theta levels in cloud
     &, thetav_cld(n_sh,nlev)                                           &
                                    ! thetav on cloud levels
     &, r2rho_cld(n_sh,nlev)                                            &
                                    ! rho on uv cloud levels
     &, r2rho_theta_cld(n_sh,nlev)                                      &
                                    ! rho on theta cloud levels
     &, p_uvcld(n_sh,nlev)                                              &
                                    ! p on uv cloud levels
     &, u_cld(n_sh,nlev)                                                &
                                    ! u on cloud levels
     &, v_cld(n_sh,nlev)                                                &
                                    ! v on cloud levels
     &, uw_cld(n_sh,nlev)                                               &
                                    ! uw on cloud levels
     &, vw_cld(n_sh,nlev)                                               &
                                    ! vw on cloud levels
     &, dr_across_rh_cld(n_sh,nlev)                                     &
                                     ! thickness of rho cloud levels
     &, dr_across_th_cld(n_sh,nlev)                                     &
                                     ! thickness of theta cloud levels
     &, wqr_cld(n_sh,nlev)                                              &
                                     ! flux of rain water (kg m /kg/s)
     &, precip_product_hcld(n_sh,nlev)                                  &
                                       !precipitation production
                                       ! (kg/m3/s) uv levels
     &, precip_product_cld(n_sh,nlev)                                   &
                                      ! precipitation production
                                      ! (kg/m3/s) theta levels
     &, wtracer_cld(n_sh,trlev,ntra) ! wtracer in cloud values(kgm/kg/s)

!
! Flux gradients on all model levels
!
      Real ::                                                           &
     &  dwthetal_dz(n_sh,nlev)                                          &
                                ! dwthetal/dz
     &, dwqt_dz(n_sh,nlev)                                              &
                                ! dwqt/dz
     &, dwthetav_dz(n_sh,nlev)                                          &
                                ! dwthetav/dz
     &, dwtracer_dz(n_sh,trlev,ntra)  ! dwtracer/dz

!
! Arrays for storing values for all points
!
      Real ::                                                           &
     &  zcld(n_sh)                                                      &
                                ! Depth of cloud layer (m)
     &, zcld_uv(n_sh)                                                   &
                                ! Depth of cloud layer (m) CMT cal
     &, mb_o_wsc(n_sh)                                                  &
                                ! mb/wstar_up
     &, root_mb_o_wsc(n_sh)                                             &
                                ! sqrt of above
     &, wstar_up3(n_sh)                                                 &
                                ! wstar_up**3 * root_mb_o_wsc
     &, wsc_o_mb(n_sh)          ! Convective velocity scale /mb
!
! Cloud base fluxes
!
      Real ::                                                           &
     &  wthetav_cb(n_sh)                                                &
                                ! wthetav at cloud base
     &, wtheta_plus_cb(n_sh)                                            &
                                ! wtheta at cloud base (in cloud value)
     &, wthetal_cb(n_sh)                                                &
                                ! wthetal at cloud base
     &, wq_plus_cb(n_sh)                                                &
                                ! wq at cloud base (in cloud value)
     &, wql_cb(n_sh)                                                    &
                                ! wql at cloud base
     &, wqt_cb(n_sh)                                                    &
                                ! wqt at cloud base
     &, wtracer_cb(n_sh,ntra)   ! wtracer at cloud base

!
! Arrays for cloud base calculations
!
      Real ::                                                           &
     &  theta_plus(n_sh)                                                &
                                ! theta at plus side of cloud base
     &, theta_minus(n_sh)                                               &
                                ! theta at negative side of cloud base
     &, exner_plus(n_sh)                                                &
                                ! exner at plus side of cloud base
     &, exner_minus(n_sh)                                               &
                                ! exner at negative side of cloud base
     &, exner_cb(n_sh)                                                  &
                                ! exner at cloud base
     &, q_plus(n_sh)                                                    &
                                ! q at plus side of cloud base
     &, q_minus(n_sh)                                                   &
                                ! q at negative side of cloud base
     &, qsat_plus(n_sh)                                                 &
                                ! qsat at plus side of cloud base
     &, qsat_minus(n_sh)                                                &
                                ! qsat at negative side of cloud base
     &, wup_cb2(n_sh)                                                   &
                                ! wup at cloud base **2
     &, dthetal_cb(n_sh)                                                &
                                ! dthetal across cloud base
     &, dthetav_cb(n_sh)                                                &
                                ! dthetav across cloud base
     &, dqt_cb(n_sh)                                                    &
                                ! dqt across cloud base
     &, dz_cb(n_sh)                                                     &
                                ! dz across cloud base
     &, dtracer_cb(n_sh,ntra)   ! dtracer across cloud base

      Real ::                                                           &
     &  du_cb(n_sh)                                                     &
                           ! du across cloud base
     &, dv_cb(n_sh)        ! dv across cloud base
!
! Further arrays for 3 level cloud base transition region
!
!   Extra level eg theta_plus1
!
      Real ::                                                           &
     &  theta_plus1(n_sh)                                               &
                             ! theta at level above plus side
     &, exner_plus1(n_sh)                                               &
                             ! exner at level above plus side
     &, exner_uvplus(n_sh)                                              &
                             ! exner at uv level above cb
     &, q_plus1(n_sh)                                                   &
                             ! q at level above plus side cloud base
     &, qsat_plus1(n_sh)                                                &
                             ! qsat at level above plus side cloud base
     &, z_plus1(n_sh)                                                   &
                             ! height of level above plus side
     &, z_plus (n_sh)                                                   &
                             ! height of level plus side
     &, p_plus1(n_sh)                                                   &
                             ! pressure of level above plus side
     &, p_plus (n_sh)                                                   &
                             ! pressure of level plus side
     &, p_uvplus(n_sh)       ! pressure of uv level above cloud base

      Real ::                                                           &
     & zlcl_ntml(n_sh)       ! Height of nearest model uv level to the
                             ! lifting condensation level.
!
! Arrays for across cloud
!
      Integer ::                                                        &
     &  ncld_thlev(n_sh)       ! number of theta levels in cloud
      Real ::                                                           &
     &  dqt_cld(n_sh)                                                   &
                               ! qt across cloud
     &, dthetal_cld(n_sh)                                               &
                               ! thetal across cloud
     &, dp_cld(n_sh)                                                    &
                               ! p across cloud
     &, dthetavdz_ad_lcl(n_sh)                                          &
                               ! dthetav/dz for a moist adiabat at lcl
     &, dthetav_cld(n_sh)      ! needed for one option of inv

!
! Base of inversion fluxes
!
      Real ::                                                           &
     & wtheta_inv(n_sh)                                                 &
                              ! wtheta  across inversion
     &,wthetal_inv(n_sh)                                                &
                              ! wthetal across inversion
     &,wq_inv(n_sh)                                                     &
                              ! flux of wq across inversion
     &,wqt_inv(n_sh)                                                    &
                              ! flux of wqt across inversion
     &,wql_inv(n_sh)                                                    &
                              ! flux of wql across inversion
     &,wthetav_inv(n_sh)                                                &
                              ! wthetav across inversion
     &,wqr_inv(n_sh)                                                    &
                              ! wqr flux at the base of inversion
     &,precip_product_inv(n_sh)                                         &
                                ! integral of precip production for
                                ! inversion
     &,wtracer_inv(n_sh,ntra)   ! wtracer flux across inversion
!
! Arrays for across inversion calculation
!
      Real ::                                                           &
     &  theta_below(n_sh)                                               &
                               ! theta on level below inversion base
     &, theta_above(n_sh)                                               &
                               ! theta on level above inversion base
     &, theta_above2(n_sh)                                              &
                               ! theta on level above inversion top
     &, q_below(n_sh)                                                   &
                               ! q on level below inversion base
     &, q_above(n_sh)                                                   &
                               ! q on level above inversion base
     &, q_above2(n_sh)                                                  &
                               ! q on level above inversion top
     &, exner_below(n_sh)                                               &
                               ! exner on level below inversion base
     &, exner_above(n_sh)                                               &
                               ! exner on level above inversion base
     &, exner_above2(n_sh)                                              &
                               ! exner on level above inversion top
     &, exner_inv(n_sh)                                                 &
                               ! exner at inversion
     &, p_below(n_sh)                                                   &
                           ! pressure on level below inversion base
     &, p_above(n_sh)                                                   &
                           ! pressure on level above inversion base
     &, p_above2(n_sh)                                                  &
                           ! pressure on level above inversion top
     &, p_inv(n_sh)                                                     &
                           ! pressure at inversion base
     &, p_inv_top(n_sh)                                                 &
                           ! pressure at inversion top
     &, dz_inv(n_sh)                                                    &
                           ! depth of inversion (m)
     &, rho_inv(n_sh)                                                   &
                           ! density at inversion (kg/m3)
     &, tracer_below(n_sh,ntra)                                         &
                                 ! tracers on theta level below inv.
     &, tracer_above(n_sh,ntra)  ! tracers on theta level above inv.
!
! Arrays for functions of height - theta levels
!

      Real ::                                                           &
     & g_func(n_sh,nlev)                                                &
                                ! g function
     &,k_func(n_sh,nlev)                                                &
                                ! k function
     &,f0_func(n_sh,nlev)                                               &
                                ! f0 function
     &,f1_func(n_sh,nlev)                                               &
                                ! f1 function
     &,fw_func(n_sh,nlev)                                               &
                                ! fw function
     &,ftheta_func(n_sh,nlev)   ! ftheta function
!
! Arrays for functions of height - uv levels
!
      Real ::                                                           &
     & fql_func(n_sh,nlev)                                              &
                                ! fql function
     &,g_hfunc(n_sh,nlev)                                               &
                                ! g function
     &,fw_hfunc(n_sh,nlev)                                              &
                                ! fw function
     &,fng_hfunc(n_sh,nlev)                                             &
                                ! fng function
     &,k_hfunc(n_sh,nlev)                                               &
                                ! k function
     &,b_hfunc(n_sh,nlev)                                               &
                                ! b function
     &,gql_func(n_sh,nlev)      ! gql function
!
! Arrays for CMT calculation
!
      Integer ::                                                        &
     &  nlcl_uv(n_sh)                                                   &
                           ! Level index for LCL
     &, ntop_uv(n_sh)      ! Level index for top of layer

      Real ::                                                           &
     &  uw(n_sh,nlev)                                                   &
                              ! U- comp stress profile (m2s-2)
     &, vw(n_sh,nlev)                                                   &
                              ! V-comp stress profile (m2s-2)
     &, p_lcl(n_sh)           ! p at LCL

!
! Arrays for CMT calculation - possible other calculations
!
      Real ::                                                           &
     & dr_across_th1(n_sh,nlev)       ! dr across theta levels

!
! Arrays for precipitation calculations
!
      Real ::                                                           &
     &  wqr(n_sh,nlev)                                                  &
                                  ! flux of rain water (kg m /kg/s)
     &, precip_product(n_sh,nlev) ! precipitation production/density
                                  ! (/s)
!
! Arrays for tracer checks
!
      Real ::                                                           &
     &  mtracer(n_sh,ntra)                                              &
                               ! column integral of tracers inc
     &, ttracer(n_sh,ntra)                                              &
                               ! column integral of tracers
     &, masscol(n_sh)     ! column integral of mass
!
!   required by check on -ve q
!
      Real :: qMinInColumn(n_sh)     ! Minimum value for q in column
                                      ! (kg/kg)
      Real :: tempx(n_sh)            ! work array

      Real, parameter :: QMIN = 1.0E-8 ! Global minimum allowed Q
!
! temporary storeage arrays
!
      Real ::                                                           &
     &  wup2_h, wup2, dmass                                             &
     &, factor                                                          &
     &, t_lcl_temp, temp1, temp2
!
! Loop counters
!

      Integer :: i,k,ktra,klev,ktop  ,j
      Integer :: max_cldlev, ilev, ibase, inv, max_cldtrlev

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
!
!  Shallow scheme assumes levels ntml+1 to ntpar are in cloud
!  The parametrisation scheme assumes the atmosphere is hydrostatic and
!  Boussinesq. (This is contrary to the assumptions made in the model
!  dynamics.)
!
!  Diagram of shallow cloud levels
!
!
!       ----------------------------  ntpar + 1
!
!       ++++++++++++++++++++++++++++   Inversion on uv level (ntpar+1)
!
!       ----------------------------  ntpar
!
!       + + + + + + + + + + + + + +
!
!       ----------------------------
!          "          "
!          "          "               Centre of cloud
!          "          "
!       + + + + + + + + + + + + + +
!
!       ----------------------------  ntml + 1   (plus level)
!
!       ++++++++++++++++++++++++++++  cloud base on uv level (ntml+1)
!
!       ----------------------------  ntml       (minus level)
!
!
!  Key
!  ----
!
!  -------   theta levels
!
!  + + + +   uv levels
!
!
!-----------------------------------------------------------------------
! 1.0 Initialise arrays and variables
!-----------------------------------------------------------------------
!
! Initialise some PC2 variables for safety.
!
      IF (lhook) CALL dr_hook('SHALLOW_TURB_CONV',zhook_in,zhook_handle)
      Do k = 1,nlev
        Do i = 1,n_sh
          dqclbydt(i,k) = 0.0
          dqcfbydt(i,k) = 0.0
          dbcfbydt(i,k) = 0.0
          dcflbydt(i,k) = 0.0
          dcffbydt(i,k) = 0.0
        End Do
      End Do

!-----------------------------------------------------------------------
! 1.1  Initialise precipitation, dth/dt, dq/dt, du/dt, dv/dt, tracer
!      increment arrays and cca
!-----------------------------------------------------------------------

      Do k = 1,nlev
        Do i = 1,n_sh
           dthbydt(i,k) = 0.0
           dqbydt(i,k)  = 0.0
        End Do
      End Do

      If (L_mom) then
        Do k = 1,nlev+1
          Do i = 1,n_sh
            dubydt(i,k) = 0.0
            dvbydt(i,k) = 0.0
          End Do
        End Do
      End If  ! L_mom

! Initialise tracer increments if tracers present.

      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 1,trlev
            Do i = 1,n_sh
              dtrabydt(i,k,ktra) = 0.0
            End Do
          End Do
        End Do
      End If  ! L_tracer

!-----------------------------------------------------------------------
! 1.2  Initialise diagnostic arrays selected by STASH flags
!-----------------------------------------------------------------------
      Do k = 1,nlev
        Do i = 1,n_sh
          ccw(i,k)    = 0.0     !
        End Do
      End Do

      Do i = 1,n_sh
         rain(i)    = 0.0     ! by definition of shallow scheme
         snow(i)    = 0.0     ! by definition of shallow scheme
         cclwp(i) =  0.0
         lcca(i) = 0.0
         cape_out(i) = 0.0
         tcw(i) =0.0
       End Do


      If (L_mom) then
        If (flg_uw_shall) then
          Do k = 1,nlev
            Do i = 1,n_sh
              uw_shall(i,k) = 0.0
            End Do
          End Do
        End If
        If (flg_vw_shall) then
          Do k = 1,nlev
            Do i = 1,n_sh
              vw_shall(i,k) = 0.0
            End Do
          End Do
        End If
      End If  ! L_mom


! diagnostics on full levels (i.e. all model levels)
! calculations on just cloud levels  - need 2 arrays for each ?

      Do k = 1,nlev
        Do i = 1,n_sh
           uw(i,k)   = 0.0
           vw(i,k)   = 0.0
           wtheta(i,k) = 0.0
           wq(i,k)     = 0.0
           wqt(i,k)     = 0.0
           wqr(i,k)     = 0.0
           wql(i,k)     = 0.0
           wthetav(i,k) = 0.0
           wthetal(i,k) = 0.0
           dwql_dz(i,k) = 0.0
           up_flux(i,k) = 0.0
           wup(i,k)     = 0.0
           qlup(i,k)    = 0.0
           precip_product(i,k) = 0.0
           eta(i,k)      = 0.0
           eta_half(i,k) = 0.0
        End Do
      End Do

      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 1,trlev
            Do i = 1,n_sh
              wtracer(i,k,ktra) = 0.0
            End Do
          End Do
        End Do
      End If  ! L_tracer
!
! work out maximum level values - used later to reduce loops
!
      ntpar_max = 0      !
      ntml_max  = 0
      max_cldlev = 0      ! maximum number of cloud levels

      Do i = 1,n_sh

        If (ntpar(i) >  ntpar_max) then
          ntpar_max=ntpar(i)
        End If
        If (ntml(i) >  ntml_max) then
          ntml_max=ntml(i)
        End If
!
! number of cloud levels
!
        ncld_thlev(i) = ntpar(i) - ntml(i)+1

        If (ncld_thlev(i) >  max_cldlev) then
          max_cldlev = ncld_thlev(i)
        End If

      End do

! tracers may be on less levels

      If (l_tracer) then
        max_cldtrlev = max_cldlev
        If (trlev <  max_cldtrlev) max_cldtrlev=trlev-1
      End If

      ntpar_max = ntpar_max+1    ! 1 more than value

!-----------------------------------------------------------------------
! Distance between rho levels (not thickness of levels in the case
! of k=1)

      Do k=1,nlev-1
        Do i=1,n_sh
          dr_across_th1(i,k) = r_rho(i,k+1) - r_rho(i,k)
        End do
      End do

!-----------------------------------------------------------------------
! 2.0 Calculate quantities required at all gridpoints
!-----------------------------------------------------------------------
! wstar_dn (from BL) sub cloud layer velocity scale
! wstar_up  cumulus layer convective velocity scale
! wstar_up = (mb CAPE)**1/3
! mb = 0.04*wstar_dn - cloud base mass flux
! zlcl_uv  - convection cloud base  from conv_diag
! ztop_uv  - base of inversion  from conv_diag
! zcld - depth of shallow cloud
! delthuv  - from conv_diag some form of dz
!-----------------------------------------------------------------------
      Do i = 1,n_sh

        mb(i)  = 0.04 * wstar_dn(i)
! Check this may wish to use rho and dz rather than delthuv
        wstar_up(i) = (delthvu(i) * mb(i) * G / (theta(i,ntml(i))       &
     &               * (1.0 + C_VIRTUAL * q_mix(i,ntml(i)))))**0.3333

! do the 2 methods give the same answer ? No
        zcld_uv(i) = ztop_uv(i) - zlcl_uv(i)
        zlcl_ntml(i) = z_rho(i,ntml(i)+1)


! correct version for theta levels
        zcld(i) = z_rho(i,ntpar(i)+1) - z_rho(i,ntml(i)+1)

!
! pre-calculate frequently used expressions
!
        wsc_o_mb(i) = wstar_up(i)/mb(i)

        mb_o_wsc(i) = mb(i)/wstar_up(i)

        root_mb_o_wsc(i) = sqrt(mb_o_wsc(i))

!
! Shallow convective cloud amounts etc
! Taken from mods for HadGEM1
!
        iccb(i) = ntml(i)+1        !
        lcbase(i) = iccb(i)
        icct(i) = ntpar(i)+1       !
        lctop(i) = icct(i)
        cca_2d(i) = 2.0*mb_o_wsc(i)
        lcca(i) = cca_2d(i)
!
! Ensemble vertical velocity at cloud base
!
        wup_cb2(i) = 1.2 * wstar_dn(i) * wstar_dn(i)

!
! The convective parametrisation scheme assumes that the enviromental
! air has qcl=0 and qcf=0. i.e thetal=theta and qt=q (this may not be
! the case in the model).


!
! values across cloud
!
        dthetal_cld(i) = theta(i,ntpar(i)+1)-theta(i,ntml(i))

        dqt_cld(i)= q_mix(i,ntpar(i)+1) - q_mix(i,ntml(i))


        dp_cld(i) = p_layer_boundaries(i,ntpar(i))                      &
     &                       - p_layer_boundaries(i,ntml(i))

        dthetav_cld(i) =    theta(i,ntpar(i)+1)*                        &
     &      (1.+q_mix(i,ntpar(i)+1)/repsilon)/(1.+q_mix(i,ntpar(i)+1))  &
     &                  -   theta(i,ntml(i))*                           &
     &      (1.+q_mix(i,ntml(i))/repsilon)/(1.+q_mix(i,ntml(i)))

      End Do

!-----------------------------------------------------------------------
! 3.0  Cloud base calculations
!-----------------------------------------------------------------------
!
!  Cloud base arrays
!
      Do i= 1, n_sh
!
! Values on upper side of cloud base
!
        exner_plus(i) = exner_layer_centres(i,ntml(i)+1)
        p_plus(i) = p_layer_centres(i,ntml(i)+1)
        theta_plus(i) = theta(i,ntml(i)+1)
        q_plus(i)     = q_mix(i,ntml(i)+1)
        z_plus(i)     = z_theta(i,ntml(i)+1)
        qsat_plus(i) = qse(i,ntml(i)+1)
!
! Values on theta level above + side of cloud base
!
        exner_plus1(i) = exner_layer_centres(i,ntml(i)+2)
        p_plus1(i) = p_layer_centres(i,ntml(i)+2)
        theta_plus1(i) = theta(i,ntml(i)+2)
        q_plus1(i)     = q_mix(i,ntml(i)+2)
        z_plus1(i)     = z_theta(i,ntml(i)+2)

!
! Saturation mixing ratio at cloud base  +
!
        qsat_plus1(i) = qse(i,ntml(i)+2)
!
! Values on lower side of cloud base
!
        exner_minus(i) = exner_layer_centres(i,ntml(i))
        theta_minus(i) = theta(i,ntml(i))
        q_minus(i)     = q_mix(i,ntml(i))
        qsat_minus(i)     = qse(i,ntml(i))
!
! Values at cloud base
!
        exner_cb(i) = exner_layer_boundaries(i,ntml(i))
!
! Values at uv level above cloud base
!
        exner_uvplus(i) = exner_layer_boundaries(i,ntml(i)+1)
        p_uvplus(i)     = p_layer_boundaries(i,ntml(i)+1)
!
! values across cloud base
!


        dz_cb(i) = z_rho(i,ntml(i)+1) - z_rho(i,ntml(i))
!
! dthetav/dz along a moist adiabat at LCL ntml level below cloud base ?
!
        t_lcl_temp = theta(i,ntml(i))*exner_layer_centres(i,ntml(i))
        temp1 = qse(i,ntml(i))*lc/(rv*t_lcl_temp*t_lcl_temp)

        temp2=-gamma_dry*(1.0+t_lcl_temp*temp1*recip_epsilon)            &
                                     /(1.+lcrcp*temp1)

        dthetavdz_ad_lcl(i) = temp2*(1.+0.61*theta(i,ntml(i))*temp1)+    &
             gamma_dry +0.61*g*qse(i,ntml(i))*theta(i,ntml(i))/(R*t_lcl_temp)

       End Do

       If (L_tracer) then
         Do ktra = 1,ntra
           Do i = 1,n_sh
             dtracer_cb(i,ktra) = tracer(i,ntml(i)+1,ktra) -            &
     &                                          tracer(i,ntml(i),ktra)
           End Do
         End Do
       End If  ! L_tracer
!
! Cloud base fluxes
!

! DEPENDS ON: shconv_cloudbase
        Call shconv_cloudbase( n_sh, ntra, l_tracer, mb                 &
     &,                     theta_plus,q_plus,qsat_plus,exner_plus      &
     &,                     theta_minus,q_minus,exner_minus             &
     &,                     exner_cb, wthvs, dtracer_cb                 &
     &,                     dthetal_cb, dqt_cb,dthetav_cb               &
     &,                     wtheta_plus_cb, wthetal_cb,wq_plus_cb       &
     &,                     wqt_cb,wql_cb                               &
     &,                     wthetav_cb,wtracer_cb)

!-----------------------------------------------------------------------
! 5.0 Calculate eta - non-dimensional height
!-----------------------------------------------------------------------
!  eta = (z -zcbase)/(zcld)
!-----------------------------------------------------------------------
! Compression to convective cloud levels - Loop order may be wrong?
!           Yet to decide most efficient way to do this.
!-----------------------------------------------------------------------
! field on theta levels in cloud
! Note fill arrays even if above cloud top


      Do k=1,max_cldlev+1
        Do i = 1,n_sh
          klev = ntml(i)+k
          ktop = ntpar(i)+1
          If (k <= ktop) then
            eta(i,k)   = (z_theta(i,klev) - z_rho(i,ntml(i)+1))/zcld(i)
          End If
          q_cld(i,k)      = q_mix(i,klev)
          qsat_cld(i,k)   = qse(i,klev)
          theta_cld(i,k)  = theta(i,klev)
          thetav_cld(i,k) = theta(i,klev)*(1+q_mix(i,klev)/repsilon)    &
     &                             /(1.+q_mix(i,klev))
          exner_cld(i,k)  = exner_layer_centres(i,klev)
          p_cld(i,k)  = p_layer_centres(i,klev)
          tc_cld(i,k) = theta(i,klev)*exner_cld(i,k)-ZeroDegC
          t_cld(i,k) = theta(i,klev)*exner_cld(i,k)
          ztheta_cld(i,k)   = z_theta(i,klev)
          r2rho_theta_cld(i,k)   = r2rho_th(i,klev)

        End do
      End do

      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 1,max_cldtrlev+1
            Do i = 1,n_sh
              klev = ntml(i)+k
              tracer_cld(i,k,ktra) = tracer(i,klev,ktra)
            End Do
          End Do
        End Do
      End If  ! L_tracer

! Fields on uv levels in cloud including cloud base and top


      Do k=1,max_cldlev+1
        Do i = 1,n_sh
           klev = ntml(i)+k
           ktop = ntpar(i)+1
           If (k <  ktop) then
            eta_half(i,k)= (z_rho(i,klev) - z_rho(i,ntml(i)+1))/zcld(i)
           End If
           exner_uvcld(i,k)  = exner_layer_boundaries(i,klev-1)
           p_uvcld(i,k)      = p_layer_boundaries(i,klev-1)
           z_rho_cld(i,k)   = z_rho(i,klev)
           r2rho_cld(i,k)   = r2rho(i,klev)
           rho_cld(i,k)     = rho(i,klev)
           If (k == ktop) then
             eta_half(i,k) = 1.0
           End If
         End do
       End do
!-----------------------------------------------------------------------
!    Calculation of precipitation
!-----------------------------------------------------------------------
      If (ishall_precip == 1) then
       Do i = 1,n_sh
         dz_inv(i) = z_theta(i,ntpar(i)+1) -z_theta(i,ntpar(i))
         rho_inv(i) = rho(i,ntpar(i)+1)
       End do

! DEPENDS ON: shconv_precip
       Call shconv_precip( n_sh, nlev, max_cldlev                       &
     &,               zcld, mb, wstar_up, ql_ad, dz_inv, rho_inv        &
     &,               eta_half,eta                                      &
     &,               rain, wqr_inv, precip_product_inv                 &
     &,               wqr_cld, precip_product_cld,precip_product_hcld)


      End If

!-----------------------------------------------------------------------
! 4.0  Work out fluxes at inversion
!-----------------------------------------------------------------------
! Calculations use theta and q at either side of the inversion.
! The inverson is defined as the ntpar+1 uv level.
!-----------------------------------------------------------------------
! extract values either side of inversion required for inversion
! calculations.

       Do i=1,n_sh
         q_below(i)  = q_mix(i,ntpar(i))
         q_above(i)  = q_mix(i,ntpar(i)+1)
         q_above2(i) = q_mix(i,ntpar(i)+2)
         theta_below(i)  = theta(i,ntpar(i))
         theta_above(i)  = theta(i,ntpar(i)+1)
         theta_above2(i) = theta(i,ntpar(i)+2)
         exner_below(i)  = exner_layer_centres(i,ntpar(i))
         exner_above(i)  = exner_layer_centres(i,ntpar(i)+1)
         exner_above2(i) = exner_layer_centres(i,ntpar(i)+2)
         exner_inv(i)  = exner_layer_boundaries(i,ntpar(i))
         p_below(i)  = p_layer_centres(i,ntpar(i))
         p_above(i)  = p_layer_centres(i,ntpar(i)+1)
         p_above2(i) = p_layer_centres(i,ntpar(i)+2)
         p_inv(i)     = p_layer_boundaries(i,ntpar(i))
         p_inv_top(i) = p_layer_boundaries(i,ntpar(i)+1)
         p_lcl(i)     = p_layer_boundaries(i,ntml(i))
       End do

! Note currently no check that tracer levels are > ntpar +1

       If (L_tracer) then
         Do ktra = 1,ntra
           Do i = 1,n_sh
             tracer_below(i,ktra) = tracer(i,ntpar(i),ktra)
             tracer_above(i,ktra) = tracer(i,ntpar(i)+1,ktra)
           End Do
         End Do
       End If  ! L_tracer
!
! Calculate fluxes at inversion
!

! DEPENDS ON: shconv_inversion
       Call shconv_inversion( n_sh, ntra, ishall_precip, l_tracer      &
     &,                     mb, wstar_up,wstar_dn,zcld,zlcl_ntml        &
     &,                     wthvs, dthetav_cb,dthetav_cld               &
     &,                     theta_below, theta_above, theta_above2      &
     &,                     q_below,q_above, q_above2                   &
     &,                     p_below, p_above, p_above2                  &
     &,                     p_inv,p_inv_top,dp_cld                      &
     &,                     exner_below, exner_above,exner_above2       &
     &,                     exner_inv                                   &
     &,                     tracer_below,tracer_above                   &
     &,                  wqr_inv, precip_product_inv                    &
     &,                  wtheta_inv,wq_inv,wql_inv,wthetav_inv          &
     &,                  wthetal_inv,wqt_inv,wtracer_inv)


!-----------------------------------------------------------------------
! 5.1 Calculate LES functions of ETA - non-dimensional height
!-----------------------------------------------------------------------
!
!  functions required at theta levels
!
! fw, G for Mass flux and vertical velocity
! ftheta, f0 and f1 for qlup

! DEPENDS ON: les_shall_func_thlev
      Call les_shall_func_thlev(n_sh,max_cldlev,nlev,eta,wsc_o_mb,1     &
     &,     fw_func, g_func, k_func, f0_func, f1_func                   &
     &,     ftheta_func)

!  functions required at uv levels
!
!   K, fng, B for turbulent transports.
! fw, G for Mass flux and vertical velocity

! DEPENDS ON: les_shall_func_rhlev
      Call les_shall_func_rhlev(n_sh,max_cldlev,nlev,eta_half,1         &
     &,     k_hfunc, fng_hfunc, b_hfunc, fw_hfunc, g_hfunc, fql_func    &
     &,     gql_func )


!-----------------------------------------------------------------------
! 6.0 Calculate ensemble vertical velocity and mass flux
!     On eta_half levels or full levels ?
!     At present doing both.
!-----------------------------------------------------------------------

      Do i = 1,n_sh
        wup_cb2(i) = 1.2*wstar_dn(i)*wstar_dn(i)

! Also reset mb and root_mb_o_wsc for use in all incloud cal etc
! Problem wthvs can be zero from BL in full model
        If (wthvs(i) == 0.0) THen
          mb_new(i) = 0.2*wstar_dn(i)
        Else
! added condition on dthetav_cb to avoid problems
          mb_new(i) = 0.2*wthvs(i)/                                     &
     &                    (max(dthetav_cb(i),wthvs(i)/wstar_dn(i)))
        End If

        mb_o_wsc(i) = mb_new(i)/wstar_up(i)

        root_mb_o_wsc(i) = sqrt(mb_o_wsc(i))

        wstar_up3(i)=root_mb_o_wsc(i) *(wstar_up(i)**3)

      End do


      Do k=1,max_cldlev
        do i = 1,n_sh
          wup2 = wup_cb2(i)+wstar_up(i)*wstar_up(i)*fw_func(i,k)
          wup2_h = wup_cb2(i)+wstar_up(i)*wstar_up(i)*fw_hfunc(i,k)
          mf_cld(i,k) = (mb_new(i)*wup_cb2(i)                           &
     &                          + wstar_up3(i)*g_func(i,k)) /wup2
          mf_h_cld(i,k)= (mb_new(i)*wup_cb2(i)                          &
     &                           + wstar_up3(i)*g_hfunc(i,k))/wup2_h
          wup_cld(i,k) = sqrt(wup2)
          wup_h_cld(i,k) = sqrt(wup2_h)
        End do
      End do


!-----------------------------------------------------------------------
! 7.0 Calculate qlup (ie CCW) (all convective cloud levels)
!-----------------------------------------------------------------------
! This routine is purely diagnostic, calculating qlup (CCW)

! DEPENDS ON: shconv_qlup
      Call shconv_qlup( n_sh, max_cldlev                                &
     &,                     mb_new,  mb_o_wsc, wstar_up, zcld           &
     &,                     wthetal_cb,wqt_cb                           &
     &,                     dthetal_cb,  dqt_cb                         &
     &,                     dthetal_cld, dqt_cld                        &
     &,                     theta_cld, q_cld, exner_cld                 &
     &,                     ftheta_func,f0_func,f1_func                 &
     &,                     ql_up_cld)

!-----------------------------------------------------------------------
! 8.0 Calculate liquid water flux and dqsat/dT in cloud
!-----------------------------------------------------------------------

! DEPENDS ON: shallow_wql
      Call shallow_wql(n_sh, nlev, max_cldlev                           &
     &,                ncld_thlev                                       &
     &,                zcld, wstar_up,root_mb_o_wsc                     &
     &,                wql_cb, wql_inv                                  &
     &,                qsat_cld, q_cld, theta_cld                       &
     &,                t_cld,    fql_func                               &
     &,                wql_cld                                          &
     &,                dqsatdt_cld)


      If (ishall_precip == 1) then
! add on wqr to wql calculated so far to take account of precip

        Do k=1,max_cldlev
          do i = 1,n_sh
            wql_cld(i,k) = wql_cld(i,k) + wqr_cld(i,k)
          End do
        End do

      End If

!-----------------------------------------------------------------------
! 9.0 turbulent transports  - implicit Fluxes on half levels
!                             therefore need cloud base and cloud top
!                             values.
!
!  Requires mass_flux, w_up, k_func, fng_func, b_func
!           wthetal_cb, wqt_cb
!
! Note
!  w'h' = cp* w'theta' +(L/exner)*w'q'
!       = cp*w'thetal' + (L/exner)*w'qt'
!-----------------------------------------------------------------------
! Problems here because of levels required for this calculation



! 9.1 Full K and B functions on eta_half

! Directive added to stop unrolling of outer k loop making code optimisation
! dependent on max_cldlev which varies with PE.
!CDIR NOUNROLL
      Do k =1,max_cldlev
        Do i=1,n_sh
          k_hfunc(i,k) = mf_h_cld(i,k)*(wup_h_cld(i,k)/wstar_up(i))     &
     &                        *zcld(i)*k_hfunc(i,k)
        End do
      End do


! 9.2 Evaluation of thetav at t+1


! DEPENDS ON: shallow_grad_h
       Call shallow_grad_h(n_sh, max_cldlev, nlev                       &
     &,                      timestep, ncld_thlev                       &
     &,                      z_rho_cld, ztheta_cld, r2rho_cld           &
     &,                      r2rho_theta_cld                            &
     &,                      thetav_cld, k_hfunc                        &
     &,                      wthetav_inv                                &
     &,                      dthetavdz )


! Tracer at t+1
! loop over tracers calculating dtracer/dz

      If (L_tracer) then

        Do ktra = 1,ntra

! Note kdtracerdz dimensioned on nlev so subroutine works
! DEPENDS ON: shallow_grad_h
          Call shallow_grad_h(n_sh, max_cldtrlev, nlev                  &
     &,                      timestep, ncld_thlev                       &
     &,                      z_rho_cld, ztheta_cld, r2rho_cld           &
     &,                      r2rho_theta_cld                            &
     &,                      tracer_cld(1,1,ktra), k_hfunc              &
     &,                      wtracer_inv(1,ktra)                        &
     &,                      kdtracerdz(1,1,ktra) )
        End Do

      End If  ! L_tracer
!--------------------------------------------------------------------
! 9.3 Evaluation wtheta and wq, wthetav, wqt, wthetal fluxes in cloud


! DEPENDS ON: shconv_turb_fluxes
       Call shconv_turb_fluxes(n_sh, ntra, max_cldlev, nlev             &
     &,                     max_cldtrlev,trlev                          &
     &,                     ishall_precip, l_tracer                     &
     &,                      dthetavdz, kdtracerdz                      &
     &,                      fng_hfunc, gql_func                        &
     &,                      wql_cld, mf_cld, q_cld,qsat_cld            &
     &,                      dqsatdt_cld, theta_cld, exner_cld          &
     &,                      rho_cld,precip_product_hcld,wup_h_cld      &
     &,                      wthvs, zcld, wstar_up, wtracer_cb          &
     &,                      wtheta_cld, wq_cld , wthetav_cld           &
     &,                      wthetal_cld,wqt_cld, wtracer_cld)



! 9.4 store various fluxes back on full model levels for diagnostic
!     output and for calculating gradients

! Cloud base and inversion fluxes

       Do i = 1,n_sh
         ibase = ntml(i)+1
         inv   = ntpar(i) + 1
         wq(i,ibase) = wq_plus_cb(i)
         wq(i,inv)   = wq_inv(i)
         wqt(i,ibase) = wqt_cb(i)
         wqt(i,inv)   = wqt_inv(i)
         wtheta(i,ibase) = wtheta_plus_cb(i)
         wtheta(i,inv)   = wtheta_inv(i)
         wthetav(i,ibase) = wthetav_cb(i)
         wthetav(i,inv)   = wthetav_inv(i)
         wthetal(i,ibase) = wthetal_cb(i)
         wthetal(i,inv)   = wthetal_inv(i)
         wql(i,inv)  = wql_inv(i)
         wql(i,ibase)= wql_cb(i)
         wqr(i,inv)  = wqr_inv(i)
      End Do

      If (L_tracer) then
        Do ktra = 1,ntra
          Do i = 1,n_sh
            ibase = ntml(i)+1
            inv   = ntpar(i) + 1
            wtracer(i,ibase,ktra) = wtracer_cb(i,ktra)
            wtracer(i,inv,ktra)   = wtracer_inv(i,ktra)
          End Do
        End Do
      End If

!
! Expand in cloud values on uv levels
!
      Do k = 2,max_cldlev
        Do i = 1,n_sh
          inv = ntpar(i)+1
          ilev = ntml(i)+k
          If (ilev <  inv) then
             wq(i,ilev)      = wq_cld(i,k)
             wqt(i,ilev)     = wqt_cld(i,k)
             wtheta(i,ilev)  = wtheta_cld(i,k)
             wthetal(i,ilev) = wthetal_cld(i,k)
             wthetav(i,ilev) = wthetav_cld(i,k)
             wql(i,ilev)     = wql_cld(i,k)
             wqr(i,ilev)     = wqr_cld(i,k)
          End If
        End Do
      End Do

      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 2,max_cldtrlev
            Do i = 1,n_sh
              inv = ntpar(i)+1
              ilev = ntml(i)+k
              If (ilev <  inv) then
                wtracer(i,ilev,ktra) = wtracer_cld(i,k,ktra)
              End If
            End Do
          End Do
        End Do
      End If
!
! theta levels
!
      Do k = 1,max_cldlev
        Do i = 1,n_sh
           inv = ntpar(i)+1
           ilev = ntml(i)+k
           If (ilev <= inv) then
             wup(i,ilev) = wup_cld(i,k)
             ccw(i,ilev) = ql_up_cld(i,k)
             up_flux(i,ilev) = mf_cld(i,k)
             precip_product(i,ilev) = precip_product_cld(i,k)/          &
     &                                     rho_theta(i,ilev)
           End If
         End Do

      End Do

! Fluxes below cloud base
! Know fluxes at cloud base for w'q' and w'theta' assume they
! go linearly to zero at surface below cloud base
! At bottom of cloud base transition region
!  w'q' = w'qt'_cb    & w'theta' = w'thetal'_cb

       Do k=1,ntml_max
         Do i=1,n_sh
           if (k <= ntml(i)) then

             factor  = z_rho(i,k)/z_rho(i,ntml(i)+1)
             wqt(i,k)     = factor*wqt_cb(i)
             wthetal(i,k) = factor*wthetal_cb(i)

           End if
         End do
       End do

      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 1,ntml_max
            Do i = 1,n_sh
              If (k <= ntml(i)) then
                factor  = z_rho(i,k)/z_rho(i,ntml(i)+1)
                wtracer(i,k,ktra) = factor*wtracer_cb(i,ktra)
              End If
            End Do
          End Do
        End Do

      End If

! 9.5 calculate gradient of fluxes in cloud
!  Only require wqt and wthetal


! DEPENDS ON: sh_grad_flux
       Call sh_grad_flux (n_sh,ntra, nlev,trlev,ntpar_max,l_tracer      &
     &,                  r2rho,r2rho_th,dr_across_th1                   &
     &,                  wthetav, wthetal, wqt, wtracer                 &
     &,                  dwthetav_dz, dwthetal_dz, dwqt_dz, dwtracer_dz)

!
! Increments to model
! method for calculating increments to theta and rv
!
      If (ishall_precip == 1) then
       Do k = 1,nlev
         Do i = 1,n_sh

           dthbydt(i,k) = -dwthetal_dz(i,k)                             &
                + lcrcp*precip_product(i,k)/exner_layer_centres(i,k)

           dqbydt(i,k)  = -dwqt_dz(i,k) - precip_product(i,k)
         End Do
       End Do
      Else
       Do k = 1,nlev
         Do i = 1,n_sh
           dthbydt(i,k) = -dwthetal_dz(i,k)
           dqbydt(i,k)  = -dwqt_dz(i,k)

         End Do
       End Do
      End If

! tracer increments


      If (L_tracer) then
        Do ktra = 1,ntra
          Do k = 1,trlev
            Do i = 1,n_sh
              dtrabydt(i,k,ktra) = -dwtracer_dz(i,k,ktra)
            End Do
          End Do
        End Do


      End If

!
! Convective cloud liquid water path
!
       Do i = 1,n_sh
         cclwp(i) = 0.0
       End Do

       Do k=1,nlev-1
         Do i=1,n_sh

           dmass = r2rho_th(i,k)*dr_across_th(i,k)*ra2
           cclwp(i) = cclwp(i) + ccw(i,k)*dmass

! Is this still correct given up_flux nolonger upward ?

           tcw(i) = tcw(i) - ccw(i,k)*up_flux(i,k)*dmass

         End do
       End do
!-----------------------------------------------------------------------
! 10.0 Convective Momentum Transport as existing scheme (if L_mom = .T.)
!-----------------------------------------------------------------------
! But recoded to use wup and mass_flux from above plus height directly
!-----------------------------------------------------------------------
!       Because of the definition of nlcl, the pressure of the top of
!       the mixed layer is phalf_uv(nlcl,*)
!
! eta levels for this calculation different from above as cloud
!  base and top taken to be theta levels for this calculation
!-----------------------------------------------------------------------
! Initialize arrays required for Convective Momentum Transport(CMT)
!-----------------------------------------------------------------------

      If (L_mom) then

        Do i = 1,n_sh
          nlcl_uv(i)    = ntml(i) + 1
          ntop_uv(i)    = ntpar(i) + 1
        End Do

!-----------------------------------------------------------------------
! 10.1 compression to convective cloud levels - U & V
!-----------------------------------------------------------------------
! Note fill arrays even if above cloud top

        Do k=1,max_cldlev+1
          Do i = 1,n_sh
            klev = ntml(i)+k
            u_cld(i,k)   = u(i,klev)
            v_cld(i,k)   = v(i,klev)
            dr_across_rh_cld(i,k) = dr_across_rh(i,klev)
            dr_across_th_cld(i,k) = dr_across_th(i,klev)
            r2rho_cld(i,k)       = r2rho(i,klev)
            r2rho_theta_cld(i,k) = r2rho_th(i,klev)
          End do
        End do

!-----------------------------------------------------------------------
! 10.2  Convective Momentum Transport (if L_mom = .T.)
!-----------------------------------------------------------------------

! Altered value of wsc_o_mb ? here as well need original ?
! what values of mf_cld etc do I pass in ? those
! on theta, uv levels relative to zcld_th or zcld_cmt?
! need to try different ones to see which work best.

! DEPENDS ON: shtconv_grad_stress
          Call shtconv_grad_stress(n_sh, nlev,max_cldlev                &
     &,                        ncld_thlev                               &
     &,                        timestep,mb_new,wstar_up                 &
     &,                        zcld_uv                                  &
     &,                        mf_cld,wup_cld,k_func                    &
     &,                        u_cld,v_cld                              &
     &,    r2rho_cld,r2rho_theta_cld,dr_across_rh_cld,dr_across_th_cld  &
                                   ! OUT
     &,                        uw_cld,vw_cld)


!
! Expand in cloud values on uv levels
!
         Do k = 1,max_cldlev
           Do i = 1,n_sh
              inv = ntpar(i)+1
              ilev = ntml(i)+k
              If (ilev <  inv) then
                uw(i,ilev)      = uw_cld(i,k)
                vw(i,ilev)      = vw_cld(i,k)
              End If
            End Do
          End Do

!
!  change in winds across cloud base
!
          Do i = 1,n_sh
            du_cb(i) = u(i,ntml(i)+1) - u(i,ntml(i))
            dv_cb(i) = v(i,ntml(i)+1) - v(i,ntml(i))
          End Do

! DEPENDS ON: shtconv_base_stress
          Call shtconv_base_stress(n_sh,nlev,ntml,ntpar,ntpar_max       &
     &,                            timestep, mb,wsc_o_mb                &
     &,                            uw0,vw0, du_cb,dv_cb                 &
     &,                            rho_theta,z_rho,z_theta              &
     &,                            flg_uw_shall,flg_vw_shall            &
                                   ! INOUT
     &,                            uw,vw                                &
                                   ! OUT
     &,                            uw_shall,vw_shall)


! DEPENDS ON: shtconv_cmt_incr
          Call shtconv_cmt_incr(n_sh,nlev, ntpar_max                    &
     &,                          r2rho, r2rho_th                        &
     &,                          dr_across_rh                           &
     &,                          uw ,vw                                 &
                                !OUT
     &,                          dubydt,dvbydt)


      End If ! L_mom

!-----------------------------------------------------------------------
! 11.0  Check for negative or very small water vapour after convection
!-----------------------------------------------------------------------
! Trys to prevent negative q by adjusting the vertical integral of dq/dt
! in a conservative way.
!-----------------------------------------------------------------------
! Find minimum value of q in column

        Do i = 1,n_sh
          qMinInColumn(i) = q_mix(i,nlev)
        End Do
        Do k = 1,nlev-1
          Do i = 1,n_sh
            If (q_mix(i,k)  <   qMinInColumn(i)) then
              qMinInColumn(i) = q_mix(i,k)
            End If
          End Do
        End Do

!
! Ensure Q does not go below global allowed minimum (QMIN)
!
        Do i = 1,n_sh
          qMinInColumn(i)=MAX(QMIN,qMinInColumn(i))
        End Do

!
! Apply an artificial upwards flux from k-1 level to ensure Q
! remains above minimum value in the column.
! Only safe if increment below is positive (definitely a problem if do
! below ntml where all increments are negative).
!

        Do k = nlev,2,-1
          Do i = 1,n_sh

            tempx(i)=q_mix(i,k) + dqbydt(i,k) * timestep

            If (tempx(i)  <   qMinInColumn(i)) then

! safe to correct level below has a positive increment > negative
              if (dqbydt(i,k-1) >= 0.0                                  &
     &                .and.dqbydt(i,k-1) >  (-1.*dqbydt(i,k))) then
              dqbydt(i,k-1) = dqbydt(i,k-1) -                           &
     &       ((qMinInColumn(i) - q_mix(i,k)) / timestep-dqbydt(i,k))    &
     &               * (r2rho_th(i,k)*dr_across_th(i,k))                &
     &               / (r2rho_th(i,k-1)*dr_across_th(i,k-1))

              dqbydt(i,k) = (qMinInColumn(i) - q_mix(i,k)) / timestep

              Else
! unsafe to correct ?
! try checking on level 2 below
                if (dqbydt(i,k-2) >= 0.0) then

                   dqbydt(i,k-1) = dqbydt(i,k-1) -                      &
     &       ((qMinInColumn(i) - q_mix(i,k)) / timestep-dqbydt(i,k))    &
     &               * (r2rho_th(i,k)*dr_across_th(i,k))                &
     &               / (r2rho_th(i,k-1)*dr_across_th(i,k-1))

              dqbydt(i,k) = (qMinInColumn(i) - q_mix(i,k)) / timestep

                else
                   write(6,*) 'Problem negative q ',i,k,                &
     &                        tempx(i),dqbydt(i,k), q_mix(i,k),         &
     &                        (dqbydt(i,j),j=1,ntpar(i)+1)

                End If
              End If

            End If   ! resultant less than qmin
          End Do     ! n_sh loop
        End Do       ! nlev loop

! check q on bottom level

        k=1
        Do i = 1,n_sh

            tempx(i)=q_mix(i,k) + dqbydt(i,k) * timestep

            If (tempx(i)  <   qMinInColumn(i)) then
               write(6,*) 'Problem negative q k=1 ',i,tempx(i),         &
     &             dqbydt(i,k),q_mix(i,k),(dqbydt(i,j),j=1,ntpar(i)+1)

            End If
        End Do ! n_sh loop

!-----------------------------------------------------------------------
! SCM diagnostics - Note problems if Call convection more then once per
!                 timestep
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! 12.0 Convective cloud amount 3d - no anvil
!-----------------------------------------------------------------------
! Initialise output array

      Do k = 1,nlev
        Do i = 1,n_sh
          cca(i,k) = 0.0
        End Do
      End Do


! Whether PC2 or not set 3D cca values for convective cloud
! The new shallow turbulence scheme provides diagnostic CCA and CCW 
! to PC2

      Do k = 1,ntpar_max
        Do i = 1,n_sh
          If (k >= iccb(i) .and. k <= icct(i)) then  ! correct know bug 
            cca(i,k) = cca_2d(i)
          End If
        End Do    
      End Do  
      

!-----------------------------------------------------------------------
! 13.0  End Subroutine
!-----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('SHALLOW_TURB_CONV',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SHALLOW_TURB_CONV
