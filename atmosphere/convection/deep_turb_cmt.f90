! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
SUBROUTINE deep_turb_cmt (n_dp, nterm, nlev, deep_cmt_opt,               &
                          ntml, kterm, cu_term, freeze_lev,              &
                          timestep,                                      &
                          uw0, vw0, mb, wcld, wstar, zlcl,               &
                          mass_flux,                                     &
                          r_rho, r_theta, z_rho, z_theta,                &
                          rho, rho_theta,                                &
                          r2rho, r2rho_th, dr_across_th, dr_across_rh,   &
                          u, v,                                          &
                          dubydt, dvbydt, uw_diag, vw_diag)

! Modules

  USE conv_type_defs      ! numbers indicating type of convection
  USE tcs_cmt_params_dp   ! parameters from CRM fits for deep CMT

  USE earth_constants_mod, ONLY: g
  
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   This routine calculates convective momentum transport for deep convection
!   using turbulence ideas. This version is designed for use with the mass flux
!   convection scheme.
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

   INTEGER, INTENT(IN) :: &
     n_dp                 & ! number of deep columns
    ,nterm                & ! Number of deep columns which actually convected
    ,nlev                 & ! Number of model levels for calculations
    ,deep_cmt_opt           ! 3 or 4 to indicate treatment of non-local part

   INTEGER, INTENT(IN) :: &
     ntml(n_dp)           & ! Cloud base level information
    ,kterm(n_dp)          & ! cloud top level information
    ,cu_term(n_dp)        & ! Location of deep points which actually convected
    ,freeze_lev(n_dp)       ! First theta level with T < 273.15K

   REAL, INTENT(IN) ::    &
     timestep               ! Convection timestep (s)

   REAL, INTENT(IN) ::    &
     uw0(n_dp)            & ! surface uw stress (N/m2)
    ,vw0(n_dp)            & ! surface vw stress (N/m2)
    ,mb(n_dp)             & ! Cloud base mass flux (Pa/s)
    ,wcld(n_dp)           & ! convective velocity scale (m/s)
    ,wstar(n_dp)          & ! convective sub-cloud velocity scale (m/s)
    ,zlcl(n_dp)           & ! Exact height of LCL (m) ? (or model level height)
    ,mass_flux(n_dp,nlev) & ! Mass flux profile (on theta levels)
    ,r_rho(n_dp,0:nlev)   & ! Radius of rho levels (m)
    ,r_theta(n_dp,0:nlev) & ! Radius of theta levels (m)
    ,z_rho(n_dp,nlev)     & ! Height above surface of rho levels (m)
    ,z_theta(n_dp,nlev)   & ! Height above surface of theta levels (m)
    ,rho(n_dp,nlev)       & ! rho levels (kg/m3)
    ,rho_theta(n_dp,nlev) & ! rho  theta levels (kg/m3)
    ,r2rho(n_dp,nlev)     & ! r*r*rho  rho levels (kg/m)
    ,r2rho_th(n_dp,nlev)  & ! r*r*rho  theta levels (kg/m)
    ,dr_across_th(n_dp,nlev) & ! thickness of theta levels (m)
    ,dr_across_rh(n_dp,nlev) & ! thickness of rho levels (m)
    ,u(n_dp,nlev)            & ! U component of wind (m/s)
    ,v(n_dp,nlev)              ! V component of wind (m/s)

   REAL, INTENT(OUT) ::   &
     dubydt(n_dp,nlev+1)  & ! dU/dt due to deep CMT (m/s/s)
    ,dvbydt(n_dp,nlev+1)  & ! dV/dt due to deep CMT (m/s/s)
    ,uw_diag(n_dp,nlev)   & ! uw stress profile on theta levels (N/m2)
    ,vw_diag(n_dp,nlev)     ! vw stress profile on theta levels (N/m2)


! Local declarations:
   INTEGER  ::           &
      i,k,ii             &  ! loop counter
     ,klev               &  ! level in cloud
     ,ilev               &  ! level in cloud
     ,itop               &  ! top level
     ,ntop_max           &  ! maximum cloud top level
     ,max_cldlev         &  ! maximum number of cloud levels
     ,max_ntml              ! maximum number of below cloud levels

   INTEGER  ::           &
     ntml_uv(n_dp)         ! lcl for UV calculations


   REAL ::               &
     zp                  & !
    ,eta_val             & !
    , factor

   REAL ::               &
     uw(n_dp,nlev)       & ! uw stress profile on theta levels (N/m2)
    ,vw(n_dp,nlev)         ! vw stress profile on theta levels (N/m2)

! arrays compressed to deep points which actually convected

   INTEGER ::           &
     ncld_thlev(nterm)  & ! number of cloud levels
    ,nstart(nterm)        ! level for ust

   REAL ::              &
     uw_cb(nterm)       & ! uw at cloud base (N/m2)
    ,vw_cb(nterm)       & ! vw at cloud base (N/m2)
    ,uw0_term(nterm)    & ! uw surface stress compressed to nterm (N/m2)
    ,vw0_term(nterm)    & ! vw surface stress compressed to nterm (N/m2)
    ,mb_term(nterm)     & ! cloud base mass flux in (m/s)
    ,zcld(nterm)        & ! cloud depth (m)
    ,zlcl_term(nterm)   & ! lifting condensation level (m)
    ,dz_cb(nterm)       & ! depth of cloud base layer (m)
    ,wup_cb(nterm)      & ! wup/wcld at cloud base
    ,fng_nlclp1(nterm)  & ! value of non-gradient function at lcl+1
    ,du_start(nterm)    & ! du across cloud base
    ,dv_start(nterm)    & ! dv across cloud base
    ,wup_peak(nterm)    & ! peak value of wup/wcld
    ,zcld_freeze(nterm) & ! depth of cloud from frezing level to top
    ,zfreeze(nterm)     & ! height of freezing level
    ,zsurf(nterm)         ! depth of surface layer

   REAL ::                        &
     uw_cld(nterm,nlev)           & ! uw on cloud levels (N/m2)
    ,vw_cld(nterm,nlev)           & ! vw on cloud levels (N/m2)
    ,u_cld(nterm,nlev)            & ! u on cloud levels
    ,v_cld(nterm,nlev)            & ! v on cloud levels
    ,uth_cld(nterm,nlev)          & ! u on theta cloud levels
    ,vth_cld(nterm,nlev)          & ! v on theta cloud levels
    ,dr_across_rh_cld(nterm,nlev) & ! thickness of rho levels
    ,dr_across_th_cld(nterm,nlev) & ! thickness of theta levels
    ,r2rho_cld(nterm,nlev)        & ! r2*rho on rho levels
    ,r2rho_theta_cld(nterm,nlev)  & ! r2*rho on theta levels
    ,eta_rho(nterm,nlev)          & ! eta (non-dimensional in cloud depth)
                                    ! on rho levels
    ,eta_theta(nterm,nlev)        & ! eta (non-dimensional in cloud depth)
                                    ! on theta levels
    ,mf_cld(nterm,nlev)           & ! mass flux on cloud levels
    ,wup_cld(nterm,nlev)          & ! wup on cloud levels
    ,k_func(nterm,nlev)           & ! k function
    ,fng_func(nterm,nlev)           ! non-gradient function

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle


!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('DEEP_TURB_CMT',zhook_in,zhook_handle)

!-------------------------------------------------------------------------------
! Initialise output arrays
!-------------------------------------------------------------------------------

DO k=1,nlev
  DO i=1,n_dp
   dubydt(i,k) = 0.0
   dvbydt(i,k) = 0.0
   uw(i,k)     = 0.0
   vw(i,k)     = 0.0
   uw_diag(i,k)  = 0.0
   vw_diag(i,k)  = 0.0
  END DO
END DO

! set top level output to zero                   
k=nlev+1                 
  Do i=1,n_dp            
    dubydt(i,k) = 0.0            
    dvbydt(i,k) = 0.0            
  End Do  

!-------------------------------------------------------------------------------
! Picture of levels used for calculation
!-------------------------------------------------------------------------------

!   ---------------------  T,q

!   + + + + + + + + + + +  u,v   Cloud top  (kterm+1)          ^
!                                                              |
!   ---------------------  T,q   uw,vw     kterm               |
!                                                              |
!   - - - - - - - - - - -  u,v                                 |
!                                                              |
!   ---------------------  T,q  uw,vw                          |
!                                                              |
!   - - - - - - - - - - -  u,v                                zcld
!                                                              |
!   ---------------------  T,q   uw,vw                         |
!                                                              |
!   - - - - - - - - - - -  u,v                                 |
!                                                              |
!   ---------------------  T,q  uw,vw                          |
!                                                              |
!   - - - - - - - - - - -  u,v                                 |
!                                                              |
!   ---------------------  T,q  uw,vw                          |
!                                                              |
!   + + + + + + + + + + +  u,v   Cloud base  (ntml+1)  zlcl    v uw_cb
!
!   ---------------------  T,q   ntml
!
!   - - - - - - - - - - -  u,v
!
!   ---------------------  T,q level 1

!   - - - - - - - - - - -  u,v level 1

!   ---------------------  Surface
!   /////////////////////

! There are Kterm -ntml + 1 wind levels used for in cloud calculations
! (includes cloud base and cloud top winds) to produce uw, vw fluxes
! on  kterm - ntml in cloud levels.
! Below cloud base the flux is assumed to reduce linearly with height from
! the cloud base flux value.

!-------------------------------------------------------------------------------
! From this point the code operates on those points which actually convected
! which for the mass flux scheme is not equal to those first diagnosed as deep
! on this time step.
!-------------------------------------------------------------------------------

! work out maximum level values - used later to reduce loops

ntop_max = 0        ! maximum cloud top level
max_cldlev = 0      ! maximum number of cloud levels
max_ntml = 0        ! maximum below cloud

DO i = 1,nterm
  ii = cu_term(i)
  IF (kterm(ii) >  ntop_max) THEN
    ntop_max=kterm(ii)+1
  END IF
  IF (ntml(ii) >  max_ntml) THEN
    max_ntml = ntml(ii)
  END IF

  ! Number of cloud levels - i.e. number of u,v levels required for in cloud
  ! calculation.

  ncld_thlev(i) = kterm(ii) - ntml(ii)+1

  IF (ncld_thlev(i)>  max_cldlev) THEN
    max_cldlev = ncld_thlev(i)
  END IF

END DO    ! nterm

! Note this should never happen unless problem of convection going to the
! model top.
IF (ntop_max > nlev-1) THEN
  ntop_max = nlev-1
END IF

! Convert cloud base mass flux from Pa/s to m/s (divide by rho*g)

DO i=1,nterm

  ii = cu_term(i)        ! location in full input arrays

  ! Convert cloud base mass flux from Pa/s to m/s (divide by rho*g)
  mb_term(i) = mb(ii)/(g*rho(ii,ntml(ii)+1))

  ntml_uv(i) = ntml(ii) + 1

  ! cloud depth
  zcld(i)  = z_rho(ii,kterm(ii)+1) - zlcl(ii)

  ! compress uw0 and vw0 to just those point which did deep convection
  uw0_term(i) = uw0(ii)
  vw0_term(i) = vw0(ii)
  zlcl_term(i) = zlcl(ii)

END DO

!-------------------------------------------------------------------------------
! Compress to cloud levels and points which actually did deep convection
!-------------------------------------------------------------------------------
!CDIR NOUNROLL
DO k=1,max_cldlev+1
  DO i = 1,nterm
    ii = cu_term(i)
    klev = ntml(ii)+k

    IF (klev >= nlev-1) THEN   ! problem will go outside array
      klev = nlev-1
    END IF

    u_cld(i,k)   = u(ii,klev)
    v_cld(i,k)   = v(ii,klev)
    dr_across_rh_cld(i,k) = dr_across_rh(ii,klev)
    dr_across_th_cld(i,k) = dr_across_th(ii,klev)
    r2rho_cld(i,k)        = r2rho(ii,klev)
    r2rho_theta_cld(i,k)  = r2rho_th(ii,klev)

    ! mass flux converted from Pa/s to m/s

    mf_cld(i,k)  = mass_flux(ii,klev)/(g*rho_theta(ii,klev))

    IF (k <= kterm(ii)+1) THEN

      eta_rho(i,k)   = (z_rho(ii,klev)  - zlcl_term(i))/zcld(i)
      eta_theta(i,k) = (z_theta(ii,klev)- zlcl_term(i))/zcld(i)

      ! visocity function  on theta levels (i.e. uw stress levels)
      !           k(eta) = a*eta

      k_func(i,k) = a_kcmt_deep*eta_theta(i,k)

      ! non-gradient function on theta levels

      fng_func(i,k) = EXP(-eta_theta(i,k))

    ELSE       ! above cloud top

      fng_func(i,k) = 0.0
      k_func(i,k) = 0.0
      eta_rho(i,k)   = 0.0
      eta_theta(i,k) = 0.0

    END IF

  END DO
END DO


! Calculate u & v on theta cloud levels using linear interpolation

IF (deep_cmt_opt == 4) THEN
  DO k=1,max_cldlev+1
    DO i = 1,nterm
      ii = cu_term(i)
      klev = ntml(ii)+k

      IF (klev >= nlev-1) THEN
        klev = nlev-1
      END IF
      factor = (z_theta(ii,klev) - z_rho(ii,klev))/dr_across_th(ii,klev)
      uth_cld(i,k)   = u(ii,klev+1)*factor+(1.0-factor)*u(ii,klev)
      vth_cld(i,k)   = v(ii,klev+1)*factor+(1.0-factor)*v(ii,klev)
    END DO
  END DO
END IF


! wup = wup/wcld
! -------------------------------------------------------------------------
! Scheme assumes a shape for wup
! Peak value at freezing level wup = wcld  from CRM fit (provided freezing
! level is above cloud base).
! Value at cloud base ~ a+b* wstar  (sub-cloud velocity scale)
! Note encountered problems with fit. At many points in a full atmospheric run
! wcld is small or wstar bigger so that wup_cb > wcld and the shape wrong.
! Decided to enforce a condition on wup_cb/wcld <= 0.8.

DO i=1,nterm
  ii = cu_term(i)
  wup_cb(i)   = (a_wup_deep+b_wup_deep*wstar(ii))/wcld(ii)

  IF (wup_cb(i) > 0.8) THEN
  ! Don't use wstar and wcld information just assume a shape.
    wup_cb(i) = 0.8
  END IF

  IF (freeze_lev(ii) > ntml(ii)+1) THEN  ! freezing level above cloud base

  ! Value at freezing level wup_peak = wcld/wcld
    wup_peak(i)    = 1.0
    zcld_freeze(i) = z_rho(ii,kterm(ii)+1) - z_theta(ii,freeze_lev(ii))
    zfreeze(i)     = z_theta(ii,freeze_lev(ii))

  ELSE                                   ! set to cloud base value

  ! peak value assume to be cloud value of wup
    wup_peak(i)    = wup_cb(i)
    zcld_freeze(i) = zcld(i)
    zfreeze(i)     = zlcl_term(i)

  END IF

END DO

!CDIR NOUNROLL
DO k=1,max_cldlev+1
  DO i = 1,nterm
    ii = cu_term(i)     ! location in full array
    klev = ntml(ii)+k

    IF (klev < kterm(ii)+1) THEN    ! work out in cloud wup

      IF (klev < freeze_lev(ii)) THEN  ! below freezing level

        eta_val = (z_theta(ii,klev) - zlcl_term(i))                     &
                                      /(zfreeze(i)-zlcl_term(i))
        factor = wup_cb(i)*wup_cb(i)
        wup_cld(i,k) = (factor + (1.0-factor)*eta_val)**0.5

      ELSE        ! above freezing level

      ! straight line from peak to zero at top

        zp = z_theta(ii,klev) - zfreeze(i)
        wup_cld(i,k) = wup_peak(i)* (1.0-zp/zcld_freeze(i))

      END IF       ! test on freezing level

    ELSE                            ! above cloud top

      wup_cld(i,k) = 0.0

    END IF                          ! test on klev in cloud

  END DO
END DO

!-------------------------------------------------------------------------------
! Gradient part of stress
!-------------------------------------------------------------------------------

! DEPENDS ON: deep_turb_grad_stress
CALL deep_turb_grad_stress(nterm, nlev, max_cldlev, ncld_thlev     &
                      , timestep                                   &
                      , zcld                                       &
                      , mf_cld, wup_cld, k_func                    &
                      , u_cld, v_cld                               &
                      , r2rho_cld, r2rho_theta_cld                 &
                      , dr_across_rh_cld, dr_across_th_cld         &
                      ,uw_cld,vw_cld)


!-------------------------------------------------------------------------------
! Non-gradient part of stress
! Note this code does not calculate the non-gradient stress in exactly the
! same way as the old turbluence based code documented by A Grant.
!-------------------------------------------------------------------------------
IF (deep_cmt_opt == 4) THEN      ! version using m(ust-u)*(1-eta)
                                 ! This version is still experimental.

  ! need level near surface  = 0.1*zlcl
  DO i=1,nterm
    zsurf(i) = 0.1*zlcl_term(i)
  END DO
  k=1
    DO i=1,nterm
      ii = cu_term(i)
      IF (zsurf(i) <= z_theta(ii,k)) THEN
        nstart(i) = k
      END IF
    END DO
!CDIR NOUNROLL
  DO k=2,max_ntml
    DO i=1,nterm
      ii = cu_term(i)
      IF (zsurf(i) <= z_theta(ii,k) .AND. zsurf(i) > z_theta(ii,k-1)) THEN
        nstart(i) = k
      END IF
    END DO
  END DO

!CDIR NOUNROLL
  DO k=1,max_cldlev+1
    DO i = 1,nterm
      ii = cu_term(i)
      klev = ntml(ii)+k
      IF (klev < kterm(ii)+1) THEN    ! only in cloud
        uw_cld(i,k) = uw_cld(i,k)                                       &
                       + mf_cld(i,k)*(u(ii,nstart(i))-uth_cld(i,k))*    &
                          (1.0-eta_theta(i,k))
        vw_cld(i,k) = vw_cld(i,k)                                       &
                       + mf_cld(i,k)*(v(ii,nstart(i))-vth_cld(i,k))*    &
                          (1.0-eta_theta(i,k))
      END IF
    END DO
  END DO
  DO i=1,nterm
    ii = cu_term(i)
    uw_cb(i) = mb_term(i)*(u(ii,nstart(i))-u(ii,ntml(ii)+1))
    vw_cb(i) = mb_term(i)*(v(ii,nstart(i))-v(ii,ntml(ii)+1))
  END DO


ELSE              ! option 3 non-local term uses cloud base stress

! Values across cloud base - at present uses adjacent model levels
! but could be altered.

  DO i=1,nterm
    ii = cu_term(i)
    dz_cb(i) = z_rho(ii,ntml_uv(i)) - z_rho(ii,ntml_uv(i)-1)
    du_start(i) = u(ii,ntml_uv(i)) - u(ii,ntml_uv(i)-1)
    dv_start(i) = v(ii,ntml_uv(i)) - v(ii,ntml_uv(i)-1)

    !  Non gradient function at first in cloud level

    fng_nlclp1(i) = EXP(-eta_rho(i,1))  ! Required at a rho level

  END DO


! Calculate cloud base stresses

! DEPENDS ON: tcs_cb_stress
  CALL tcs_cb_stress (deep_conv, nterm                              &
                     ,timestep                                      &
                     ,uw0_term, vw0_term, mb_term, zlcl_term        &
                     ,du_start, dv_start, dz_cb, fng_nlclp1         &
                     ,uw_cb, vw_cb )

! Add non-gradient part of stress to gradient stress

!CDIR NOUNROLL
  DO k=1,max_cldlev+1
    DO i = 1,nterm

      uw_cld(i,k) = uw_cld(i,k) + uw_cb(i)*fng_func(i,k)
      vw_cld(i,k) = vw_cld(i,k) + vw_cb(i)*fng_func(i,k)

    END DO
  END DO

END IF        ! test on non-local option

!-------------------------------------------------------------------------------
! Expand in cloud values on uv levels

!CDIR NOUNROLL
DO k = 1,max_cldlev+1
  DO i = 1,nterm
    ii = cu_term(i)
    itop = kterm(ii)+1
    ilev = ntml(ii)+k
    IF (ilev < itop) THEN
      uw(ii,ilev)  = uw_cld(i,k)
      vw(ii,ilev)  = vw_cld(i,k)
    END IF
  END DO
END DO

! Below cloud base assume uwcb linearly to zero

!CDIR NOUNROLL
DO k = 1,max_ntml
  DO i = 1,nterm
    ii = cu_term(i)
    ilev = ntml(ii)
    IF (k <= ilev)  THEN
      uw(ii,k)  = uw_cb(i)*z_theta(ii,k)/zlcl_term(i)
      vw(ii,k)  = vw_cb(i)*z_theta(ii,k)/zlcl_term(i)
    END IF
  END DO
END DO

!-------------------------------------------------------------------------------
! Calculate increments to U and V - full arrays


! DEPENDS ON: tcs_cmt_incr
CALL tcs_cmt_incr(n_dp,nlev, ntop_max, kterm             &
                 ,r2rho, r2rho_th                        &
                 ,dr_across_rh                           &
                 ,uw ,vw                                 &
                 !OUT
                 ,dubydt,dvbydt)

!-------------------------------------------------------------------------------
! Copy stress to output arrays multiplying by density on theta levels
! so that diagnostic output has the correct units for stress rather than
! m2/s2.

DO k=1,ntop_max+1
  DO i=1,n_dp
    uw_diag(i,k)=uw(i,k)*rho_theta(i,k)
    vw_diag(i,k)=vw(i,k)*rho_theta(i,k)
  END DO
END DO
!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('DEEP_TURB_CMT',zhook_out,zhook_handle)

!-------------------------------------------------------------------------------
END SUBROUTINE deep_turb_cmt
