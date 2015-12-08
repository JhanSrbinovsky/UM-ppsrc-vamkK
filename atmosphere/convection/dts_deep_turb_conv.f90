! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+  Deep turbulent convection scheme
!

SUBROUTINE deep_turb_conv(                                                    &
         ! Intent IN
           call_number, nbl, nlev, n_cca_lev, ntra, n_dp, trlev               &
         , ntml, ntpar, freeze_lev                                            &
         , l_calc_dxek,l_q_interact,l_tracer                                  &
         , land_mask                                                          &
         , timestep, pstar, p_layer_centres, p_layer_boundaries               &
         , exner_layer_centres, exner_layer_boundaries, z_theta               &
         , z_rho, dr_across_th, dr_across_rh, rho, rho_theta                  &
         , r_rho, r_theta, r2rho, r2rho_th                                    &
         , qse, theta, q, qcl, qcf, bulk_cf, cf_frozen, cf_liquid             &
         , u, v, w                                                            &
         , uw0, vw0, rhowqt, sensible_heat                                    &
         , zlcl, zlcl_uv, wth0, wq0, wstar, wthvs, q1_sd, t1_sd               &
         , w_max                                                              &

         ! Intent INOUT
         , tracer                                                             &

         ! Intent OUT
         , kterm, iccb, icct                                                  &
         , dthbydt, dqbydt, dqclbydt, dqcfbydt, dbcfbydt                      &
         , dcffbydt, dcflbydt, dubydt, dvbydt, dtrabydt                       &
         , rain, snow, rain_3d, snow_3d, up_flux, cca_2d, cca                 &
         , ccw, cclwp                                                         &
         , uw_deep, vw_deep, mb_deep )

USE water_constants_mod, ONLY: lc, lf, tm

USE atmos_constants_mod, ONLY: cp, c_virtual, r

USE cv_run_mod,  ONLY:                                                        &
    l_mom, deep_cmt_opt, cca2d_dp_opt, cca_dp_knob, ccw_dp_knob, l_ccrad,     &
    l_3d_cca

USE dts_fitpars_mod, ONLY:                                                    &
   tcritplume

USE cv_param_mod, ONLY:                                                       &
   a_land, a_sea, b_land, b_sea

USE earth_constants_mod, ONLY: g, earth_radius

! subroutines used

USE dts_tracerflux_mod, ONLY:                                                 &
   dts_tracerflux

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!-----------------------------------------------------------------------
! Description:
!   Deep turbulent convection scheme 
!   
!
!   Called by GLUE_CONV.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection

! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
!-----------------------------------------------------------------------

! Subroutine Arguments

INTEGER, INTENT(IN) :: &
  call_number          & ! Sweep number for convection
 ,nbl                  & ! No. of boundary layer levels
 ,nlev                 & ! No. of model layers
 ,n_cca_lev            & ! No. of convective cloud levels
 ,ntra                 & ! No. of tracer fields
 ,n_dp                 & ! No. of deep convection points
 ,trlev                  ! No. of model levels on which tracers are included
                         ! (may be able to assume this is nlev?)

INTEGER, INTENT(IN) :: &
  ntml(n_dp)           & ! Top level of surface mixed layer defined relative
                         ! to theta,q grid (conv_diag)
 ,ntpar(n_dp)          & ! Top level of initial parcel ascent calculated by
                         !  conv_diag
 ,freeze_lev(n_dp)       ! Level index for freezing level

LOGICAL, INTENT(IN) :: & 
  l_calc_dxek          & !
 ,l_q_interact         & !
 ,l_tracer               ! .true. if tracers present

LOGICAL, INTENT(IN) :: & 
  land_mask(n_dp)        ! Land/sea mask .true. if land

REAL, INTENT(IN)    ::  &
  timestep              & ! Model timestep for convection (s)
 ,pstar(n_dp)             ! Surface pressure (Pa)


! UM prognostics grid info (you may not use all fields)

REAL, INTENT(IN)    ::                &
  exner_layer_centres(n_dp,0:nlev)    & ! Exner at theta levels
 ,exner_layer_boundaries(n_dp,0:nlev) & ! Exner at u,v levels
 ,p_layer_centres(n_dp,0:nlev)        & ! Pressure(Pa) at theta levels
 ,p_layer_boundaries(n_dp,0:nlev)     & ! Pressure  at uv(Pa)
 ,z_theta(n_dp,nlev)                  & ! height of theta levels(m)
 ,z_rho(n_dp,nlev)                    & ! height of rho levels (m)
 ,dr_across_th(n_dp,nlev)             & ! thickness of theta levels (m)
 ,dr_across_rh(n_dp,nlev)             & ! thickness of rho levels (m)
 ,rho(n_dp,nlev)                      & ! density on rho levels (kg/m3)
 ,rho_theta(n_dp,nlev)                & ! density on theta levels (kg/m3)
 ,r_rho(n_dp,nlev)                    & ! Radius rho levels (m)
 ,r_theta(n_dp,0:nlev)                & ! Radius theta levels (m)
 ,r2rho(n_dp,nlev)                    & ! r**2*density rho levels  (kg/m)
 ,r2rho_th(n_dp,nlev)                   ! r**2*density theta levels (kg/m)
      
REAL, INTENT(IN) ::     &
  qse(n_dp,nlev)          ! Saturation mixing ratio of environment (kg/kg)

! prognostics

REAL, INTENT(IN) ::  &
  theta(n_dp,nlev)   &   ! Model potential temperature (K)
 ,q(n_dp,nlev)       &   ! Model mixing ratio  (kg/kg)
 ,qcl(n_dp,nlev)     &   ! Liq condensate mix ratio (kg/kg)
 ,qcf(n_dp,nlev)         ! Ice condensate mix ratio (kg/kg)

! Other PC2 prognostics  (qcl and qcf need to be updated if PC2)

REAL, INTENT(IN) ::    &
  bulk_cf(n_dp,nlev)   &  ! Bulk total cloud volume ( )
 ,cf_frozen(n_dp,nlev) &  ! Frozen water cloud volume
 ,cf_liquid(n_dp,nlev)    ! Liq water cloud volume
     
! winds 
REAL, INTENT(IN)    :: &
  u(n_dp,nlev)         &  !Model U field (m/s) on pgrid horizontally
 ,v(n_dp,nlev)         &  !Model V field (m/s) on pgrid horizontally
 ,w(n_dp,nlev)            !Model W field (m/s)


! From BL scheme caled before convection

REAL, INTENT(IN)    :: &
  uw0(n_dp)            &  ! U-comp of surface stress (N/m2)
 ,vw0(n_dp)            &  ! V-comp of surface stress (N/m2)
 ,rhowqt(n_dp)         &  ! surface moisture flux rho*wqt (kg/m2/s)
 ,sensible_heat(n_dp)  &  ! surface sensible heat flux cp*rho*w'tl'(W/m2)
 ,zlcl(n_dp)           &  ! Lifting condensation level from conv_diag (m)
 ,zlcl_uv(n_dp)        &  ! Lifting condensation level for uv grid (m)
 ,wth0(n_dp)           &  ! ftl(dpi(j))/rho(dpi(j),1)/cp
 ,wq0(n_dp)            &  ! fqt(dpi(j))/rho(dpi(j),1)
 ,wstar(n_dp)          &  ! Convective velocity scale
                          ! (m/s) (sub cloud) COMES from BL
 ,wthvs(n_dp)          &  ! Surface flux of wthv (Km/s)
 ,q1_sd(n_dp)          &  ! Standard deviation of turbulent flucts.
                          ! of layer 1 q (kg/kg)
 ,t1_sd(n_dp)          &  ! Standard deviation of turbulent flucts.
                          ! of layer 1temp. (K)
 ,w_max(n_dp)             ! maximun w in column (m/s)

! Tracers  - no information on type just number ntra 
!            all processed in the same way - i.e. transported

REAL, INTENT(INOUT) ::    &
  tracer(n_dp,trlev,ntra)   !Model tracer fields (kg/kg)


INTEGER, INTENT(OUT) :: &
  kterm(n_dp)           & ! Actual termination level of deep convection  
 ,iccb(n_dp)            & ! Convective cloud base level number 
 ,icct(n_dp)              ! Convective cloud top level number 


REAL, INTENT(OUT) ::    &
  dthbydt(n_dp,nlev)    & ! Increments to potential temp. due to convection(K/s)

 ,dqbydt(n_dp,nlev)     & ! Increments to q due to convection (kg/kg/s)

! If PC2 then also increments to cloud prognostics

 ,dqclbydt(n_dp,nlev)   & ! Increments to liq condensate due to 
                          ! convection (kg/kg/s)
 ,dqcfbydt(n_dp,nlev)   & ! Increments to ice condensate due to
                          ! convection (kg/kg/s)
 ,dbcfbydt(n_dp,nlev)   & ! Increments to total cld volume due to convection(/s)

 ,dcffbydt(n_dp,nlev)   & ! Increments to ice cloud volume due to convection(/s)

 ,dcflbydt(n_dp,nlev)   & ! Increments to liq cloud volume due to convection(/s)

! If CMT
 ,dubydt(n_dp,nlev+1)   & ! Increments to U due to CMT (m/s2)

 ,dvbydt(n_dp,nlev+1)   & ! Increments to V due to CMT (m/s2)

! If tracers 
 ,dtrabydt(n_dp,nlev,ntra) !Increment to tracer due to convection (kg/kg/s)


! Diagnostics
REAL, INTENT(OUT) ::    &
  rain(n_dp)            & ! Surface convective rainfall (kg/m2/s)
 ,snow(n_dp)            & ! Surface convective snowfall (kg/m2/s)
 ,cca_2d(n_dp)            ! 2d convective cloud amount



! Desirable / essential for aviation - convection cloud information
! Currently being set if not PC2

REAL, INTENT(OUT) ::    &
  cca(n_dp,n_cca_lev)   & ! Convective cloud amount on model levels 
 ,ccw(n_dp,nlev)        & ! Convective cloud condensate on model levels (g/kg)
 ,cclwp(n_dp)             ! Condensed water path (kg/m2) Integrated ccw field.


! Desirable for UKCA - some form of rain /snow precipitation profile
!                      CURRENTLY set to zero
REAL, INTENT(OUT) ::    &
  rain_3d(n_dp,nlev)    & ! Convective rainfall flux (kg/m2/s)
 ,snow_3d(n_dp,nlev)      ! Convective snowfall flux (kg/m2/s)


! What else does the scheme produce that we would like as diagnostics?

REAL, INTENT(OUT) ::    &
  up_flux(n_dp,nlev)      ! Updraught mass flux (Pa/s for diag output)
                          ! units here kg/m2/s (*g to get diag output). 

! From CMT scheme

REAL, INTENT(OUT) ::    &
  uw_deep(n_dp,nlev)    & ! X-comp. of stress from deep convection(kg/m/s2)
 ,vw_deep(n_dp,nlev)      ! Y-comp. of stress from deep convection(kg/m/s2)
      
REAL, INTENT(OUT) ::    &
  mb_deep(n_dp)           ! cloud base mass flux in Pa/s 

       

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
INTEGER ::         &
  i,k,j,i_dp        ! loop counters
INTEGER ::         &
  one              & ! integer = 1
 ,icallmicro       & ! ==1 calls micro
 ,iusegradq        & !   
 ,ntparmax         & ! Max ntpar across all conv points
 ,kdo              & ! level below 
 ,kup                ! level above

INTEGER ::           &
  km40(n_dp)         & ! level identifiers
 ,km40m(n_dp)        &
 ,km40p(n_dp)        &
 ,kfrm(n_dp)         &
 ,kfrp(n_dp)         &   
 ,klcl(n_dp)         & ! Theta Level corresponding to zlcl
 ,klclm(n_dp)        & ! Theta Level corresponding to zlcl -1
 ,klclp(n_dp)        & ! Theta Level corresponding to zlcl +1
 ,ktop(n_dp)         & ! Theta Level corresponding to ztop
 ,ktopm(n_dp)        & ! Theta Level corresponding to ztop -1
 ,ktopp(n_dp)        & ! Theta Level corresponding to ztop +1
 ,k_ad(n_dp)         &
 ,k_adm(n_dp)        &
 ,k_adp(n_dp)        &
 ,iconvclass(n_dp)   & ! what type of convection is it?
 ,dts_ntpar(n_dp)    & ! ntpar calc'd from dts_cape program
 ,dts_ntpardum(n_dp)   !dummy variable

REAL ::                   &
  cape_below_fr(n_dp)     & ! cape up to freezing level / J/kg
 ,cape_below_frdum(n_dp)  & ! cape up to freezing level / J/kg
 ,cape_whole_layer(n_dp)  & ! cape for whole layer / J/kg
 ,cape_above_fr(n_dp)     & ! cape above freezing level / J/kg
 ,cin(n_dp)               & ! convective inhibition / J/kg
 ,cinabove(n_dp)          & ! convective inhibition / J/kg
 ,ql_ad(n_dp)             & ! adiabatic liquid water content at
                            ! height of max buoy excess / kg/kg
 ,ql_addum(n_dp)          & ! adiabatic liquid water content dummy
 ,h_ad(n_dp)              & ! height of max buoy excess / m 
 ,h_addum(n_dp)           & ! dummy height of max buoy excess / m 
 ,diffmax(n_dp)           &
 ,diffmaxdum(n_dp)        &
 ,mb(n_dp)                & ! mass flux at cloud base
 ,mfr(n_dp)               & ! mass flux at freezing level
 ,w2lcl(n_dp)             & ! velocity variance at the lcl
 ,wcld(n_dp)              & ! vel scale for conv layer
 ,wfr(n_dp)               & ! vel scale above freez level
 ,wall(n_dp)              & ! vel scale for whole layer
 ,cstar(n_dp)             & ! scale for c-e = mb*ql_ad/(h_ad-zlcl)
 ,sigma(n_dp)             & ! cloud base fractional area
 ,qrainmax(n_dp)          & ! scaling param for max rain amount
 ,qsatsurf(n_dp)          & ! qsat at the surface
 ,qvatfr(n_dp)            & ! q at freezing level
 ,zfr(n_dp)               & ! height of fr lev (intepolated)
 ,zm40(n_dp)              & ! height of -40degC (intepolated)
 ,ztop(n_dp)              & ! z_theta(ntpar)
 ,zlcl_cape(n_dp)         & ! lcl from cape routine -- not as good
 ,zlcl_dum(n_dp)          & ! dummy lcl variable
 ,t40(1)                    ! value -40degC
REAL ::                   &
  temperature(n_dp,nlev)  & ! temperature in K
 ,qsat_moist_ad(n_dp,nlev)& ! qsat along parcel ascent
 ,qsat_moist_addum(n_dp,nlev)& ! dummy qsat along parcel ascent
 ,buoyfl(n_dp,nlev)       & ! this is g/theta_v <w'theta_v'>
 ,wthv(n_dp,nlev)         & ! theta_v flux
 ,wmse(n_dp,nlev)         & ! mse flux 
 ,wqv(n_dp,nlev)            ! qv flux
REAL  ::                        &
  ww(n_dp,nlev)                 & ! velocity variance
 ,wwrho(n_dp,nlev)              & ! velocity variance on rho levs
 ,massfl(n_dp,nlev)             & ! mass flux from vel variance
 ,massfl_rho(n_dp,nlev)         & ! mass flux from vel variance
 ,dthvdt_flux(n_dp,nlev)        & ! theta_v increment due to flux
 ,dqhdt_flux(n_dp,nlev)         & ! qh flux increment due to flux
 ,mse_withflux(n_dp,nlev)       & ! new mse value once flux incs applied
 ,thetav_withflux(n_dp,nlev)    & !new theta_v once flux incs applied
 ,thetav_fluxsmth(n_dp,nlev)    & !new theta_v once flux incs
 ,dthvdt_smth(n_dp,nlev)        & !smoothing increment to theta_v 
                                  ! applied and smoothing applied
 ,wthv_smth(n_dp,nlev)          & !flux implied by smoothing
 ,qhyd(n_dp,nlev)               & ! hydrometeor content
 ,theta_withflux(n_dp,nlev)     & ! new theta once flux incs applied 
 ,q_withflux(n_dp,nlev)         & ! new q once flux incs applied
 ,qsewat(n_dp,nlev)             & ! qsat w r t water
 ,condensation(n_dp,nlev)       & ! condensation rate kg/kg/s
 ,deposition(n_dp,nlev)         & ! deposition rate kg/kg/s
 ,sublimation(n_dp,nlev)        & ! sublimation rate kg/kg/s
 ,evaporation(n_dp,nlev)          ! evaporation rate kg/kg/s

REAL   ::            &
  depint(n_dp)       & ! Integral of deposition
 ,freezint(n_dp)     & ! conversion terms
 ,condint(n_dp)      & ! Integral of condensation
 ,qcfint(n_dp)       & ! Integral of qcf
 ,qclint(n_dp)       & ! Integral of qcl
 ,subint(n_dp)       & ! integral of sublimation
 ,meltint(n_dp)      & ! integral of melting
 ,evapint(n_dp)      &
 ,revpint(n_dp)                 
REAL   ::                     &
  scalefac1(n_dp)             &
 ,scalefac2(n_dp)             &
 ,scalefac3(n_dp)             &               
 ,rainrate(n_dp)              &
 ,snowrate(n_dp)              &
 ,starting_heights(n_dp)      &
 ,pnbabove(n_dp)              & ! level of neutral buoyancy from ascent
                                ! starting at the freezing level
 ,pnb(n_dp)                     ! level of neutral buoyancy from ascent
                                ! starting in the boundary layer
REAL   ::                     &
  melt(n_dp,nlev)             &
 ,freeze(n_dp,nlev)           &
 ,rainprod(n_dp,nlev)         &
 ,snowprod(n_dp,nlev)         &
 ,revp(n_dp,nlev)             &
 ,storethvp(n_dp,nlev)        & ! parcel virtual potential temperature
 ,storethvp_upper(n_dp,nlev)  & ! parcel virtual pot temp
 ,dthvdz_m(n_dp,nlev)         &
 ,dthvdz(n_dp,nlev)           &
 ,dqsedz(n_dp,nlev)           &
 ,qcl_plume(n_dp,nlev)        & ! plume qcl 
 ,qcf_plume(n_dp,nlev)        & ! plume qcf
 ,h1all(n_dp,nlev)

REAL   ::            &
  rtimestep          &  ! 1/timestep         (/s)
 ,delp               &  ! 1/layer depth in pressure (/Pa)
 ,th_excess          &  ! theta excess for parcel
 ,frac               &  !
 ,qlcrit             &  ! ql critical value for setting up plume water
 ,qlmin                 ! minimum water for plume

REAL   ::               &
  thetav(n_dp,nlev)     & ! environment virtual pot temp
 ,mse(n_dp,nlev)          ! moist static energy 

LOGICAL ::      &
  l_increase      !   


! Local variables required for CMT
 INTEGER ::        &
   cu_term(n_dp)      ! Index of columns actually DOing deep convection
                      ! all true for this scheme.

 INTEGER ::        &
   ncmt            &  ! number of deep columns with non-zero wcld value
  ,icall              ! type of dts_capecall

! required by water conservation check
 REAL ::               &
   qMinInColumn(n_dp)  &  ! Minimum value for q in column(kg/kg)
  ,temp1(n_dp)            ! work array
 REAL ::               &
   r_a2                &  ! radius of earth^2
  ,dqt                 &  ! change in total water
  ,dzdt                   ! dz/dt for checking CFL

REAL, PARAMETER :: qmin = 1.0E-8 ! Global minimum allowed Q

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
!  The model grid - as far as convection is concerned 
!-----------------------------------------------------------------------
!
!
!       _________________________
!   
!       - - - - - - - - - - - - -  < Cloud top (ntpar+1)
! ntpar _________________________
! /|!  |    - - - - - - - - - - - - - 

!  |    _________________________
!  |
!  |    - - - - - - - - - - - - - 
!  |    _________________________   theta(n), q(n), 
!  |
!  |    - - - - - - - - - - - - -   rho(n), u(n), v(n) 
!  |    _________________________
! \|/ 
!       - - - - - - - - - - - - -  < Cloud base (ntml+1) (LCL)
! ntml  _________________________   ntml
!
!       - - - - - - - - - - - - -   ntml
!       _________________________  theta(2), q(2), qcl(2),qcf(2), w(2) p_th(2)
!
!       - - - - - - - - - - - - -  rho(2), u(2), v(2), p_rho(1)
!       _________________________  theta(1), q(1), qcl(1),qcf(1), w(1) p_th(1)
!
!       - - - - - - - - - - - - -  rho(1), u(1), v(1)
!       _________________________
!        / / / / / / / / / / / /       Surface    p_th(0), p_rho(0)       
!
!
!
! Note p_th  - short for p_layer_centres
!      p_rho - short for p layer boundaries
!     p_rho numbering is slightly odd but this is how it comes into 
!     convection. The original p field on rho levels has a value on rho 
!     level 1 but is not passed to convection.
!
!  Cloud base defined as the lifting condensation level.
!
!   ntml and ntpar are defined in conv_diag from a parcel ascent.
!-----------------------------------------------------------------------
! parameters are a set of parameters for functions that define various things...
! include parameters.h 

!=============================================================================
!Switches

 IF (lhook) CALL dr_hook('DEEP_TURB_CONV',zhook_in,zhook_handle)
iusegradq = 1  ! 1 for simple q flux, (0  for wthv and wmse method)

icallmicro = 1 ! this should be 1 for normal usage, 0 for debugging

!Constants
one = 1

!================================================================
! INITIALISING FIELDS SECTION

! Ensure theta and q increments are zero at the beginning of this routine

dthbydt(:,:) = 0.0
dqbydt(:,:)  = 0.0


! 1/timestep
rtimestep = 1.0/timestep

!---------------------------------------------------------
! Calculate the temperature from the potential temperature 
! nb isn't this done earlier?
DO k=1,nlev
  DO i_dp=1,n_dp
    temperature(i_dp,k) = theta(i_dp,k)*exner_layer_centres(i_dp,k)
  END DO
END DO

!---------------------------------------------------------------------
! Interpolate to get height of freezing level
! nb may need to think about whether want the level height or the
! actual height -- and which to use in calc_buoy_flux, calc_w_variance

DO i_dp=1,n_dp
  IF(freeze_lev(i_dp) < nlev) THEN 
    IF( (temperature(i_dp,freeze_lev(i_dp))  >= tm .and.   &
         temperature(i_dp,freeze_lev(i_dp)+1) < tm) .or.  &
      ! allowing for some inversion behaviour
        (temperature(i_dp,freeze_lev(i_dp))  <= tm .and.  &
         temperature(i_dp,freeze_lev(i_dp)+1) > tm) ) THEN 

      kdo = freeze_lev(i_dp)
      kup = freeze_lev(i_dp)+1
    ELSE
      kdo = freeze_lev(i_dp)-1
      kup = freeze_lev(i_dp)
    END IF
    IF(kdo <= 0) THEN 
      kdo = 1
      kup = 2
    END IF
    zfr(i_dp) = z_theta(i_dp,kdo) +                       &
                   (z_theta(i_dp,kup)-z_theta(i_dp,kdo))*           & 
                   (temperature(i_dp,kdo)-tm)/                      &
                   (temperature(i_dp,kdo)-temperature(i_dp,kup))
  ELSE

    zfr(i_dp) = z_theta(i_dp,nlev) ! unlikely but protects code

  END IF      ! test on freezing level

! check to see whether the surface is a lower temperature than 0C
! if so set zfr to zero
  IF(temperature(i_dp,1) < tm) THEN
    zfr(i_dp) = 0.0
  END IF
  kfrm(i_dp) = kdo
  kfrp(i_dp) = kup
END DO
        
!-----------------------------------------------------------------------------
! Now do the same for the -40degC level: 
! Temperature trying to locate height of:
t40(1) = tm-38. !nb or should it actually be -38C?

! Locate closest level
l_increase = .false. ! temperature decreases with height

!DEPENDS ON: dts_locate_closest_levels
CALL dts_locate_closest_levels(n_dp,nlev,one,n_dp,l_increase,temperature,t40&
                               ,km40,km40m,km40p)

DO i_dp=1,n_dp
  IF(temperature(i_dp,km40(i_dp)) >= t40(one)) THEN 
    kdo = km40(i_dp)
    kup = km40(i_dp)+1
  ELSE
    kdo = km40(i_dp)-1
    kup = km40(i_dp)
  END IF

  IF(kdo >= 1 .and. kup <= nlev) THEN 
    zm40(i_dp) = z_theta(i_dp,kdo) +                         &
                (z_theta(i_dp,kup)-z_theta(i_dp,kdo))*              &
                (temperature(i_dp,kdo)-t40(one))/                   &
                (temperature(i_dp,kdo)-temperature(i_dp,kup))
  END IF
END DO
!------------------------------------------
! Find level corresponding to zlcl
l_increase = .true. ! z_theta increases with height

!DEPENDS ON: dts_locate_closest_levels
CALL dts_locate_closest_levels(n_dp,nlev,n_dp,n_dp,l_increase,z_theta,zlcl&
                              ,klcl,klclm,klclp)

! calculate qhyd -- take from large scale
qhyd(:,:) = qcf(:,:)+qcl(:,:)

!-------------------------------------------------------------------
! Calculate virtual theta profile -- this now doesn't include hydrometeor 
! loading
DO k=1,nlev
  DO i_dp = 1,n_dp
    thetav(i_dp,k) = theta(i_dp,k)*(1.0+c_virtual*q(i_dp,k)) !nb-qhyd(i_dp,k))
  END DO
END DO
    
!-----------------------------------
! Calculate moist static energy here
DO k=1,nlev
  DO i_dp=1,n_dp
    mse(i_dp,k) = theta(i_dp,k)+q(i_dp,k)*lc/cp
  END DO
END DO

! Now set up a liquid water content for plumes
!nb also calculated in dts_cape, could be made more efficient!

qcl_plume(:,:) = 0.0

! variables for setting up the plume water content
qlcrit = 1.e-3 ! kg/kg
qlmin = 2.e-4  ! kg/kg

DO k=1,nlev
  DO i_dp=1,n_dp
    IF(z_theta(i_dp,k) > zlcl(i_dp)) THEN  
      qcl_plume(i_dp,k) = 0.5*qse(i_dp,k)
      IF(qcl_plume(i_dp,k) > qlcrit) THEN 
        qcl_plume(i_dp,k) = qlcrit
      END IF
      IF(qcl_plume(i_dp,k) < qlmin) THEN 
        qcl_plume(i_dp,k) = qlmin
      END IF
    END IF ! z_theta(i_dp,k) > zlcl(i_dp)
  END DO
END DO

! now determine in a crude way what fraction of the above is in fact ice and
! what fraction liquid: if the environmental temperature is less
! than -10C then turn it all into ice. Could improve on this!

DO k=1,nlev
  DO i_dp=1,n_dp
    IF(temperature(i_dp,k) < tm) THEN 
      frac = (temperature(i_dp,k) - tm)/tcritplume
      IF(frac > 1.0) frac = 1.0
                 
      qcf_plume(i_dp,k) = qcl_plume(i_dp,k)*(frac) 
      qcl_plume(i_dp,k) = qcl_plume(i_dp,k)*(1.0-frac)
    ELSE
      qcf_plume(i_dp,k) = 0.
    END IF
  END DO
END DO

! find the saturation value of this initial parcel

qsewat(:,:) = qse(:,:)
       
!================================================================
! SCALING PARAMETER SECTION
! Calculation of scaling parameters
! In this section we obtain the following parameters:
! cape_lower = Undilute cape up to freezing level
! qsat along a moist adiabat (a by-product of the cape calculation)
! mb = Mass flux at cloud base
! zm40 = height of -40deg
! h_ad = height of maximum buoyancy excess for undilute parcel
! cstar = mb*ql_ad/(h_ad-zlcl)

starting_heights(:) = zlcl(:)/2 ! was 200 m
th_excess = 0.0 ! K
icall=1

! DEPENDS ON: dts_cape
CALL dts_cape(n_dp,nlev,icall,                                             &
              q,theta,thetav,qcl,qcf,qse,                                  &
              p_layer_centres,                                             &
              p_layer_boundaries,exner_layer_centres,z_theta,              &
              rho_theta,zlcl_cape,zlcl,klcl,freeze_lev,starting_heights,   &
              th_excess,ntpar,dts_ntpar,                                   &
              cape_below_fr,cape_whole_layer,cin,                          &
              qsat_moist_ad,ql_ad,h_ad,diffmax,pnb,storethvp)
       
!   write(6,*) ' AFTER dts_cape or from conv_diag ',n_dp
!   write(6,*) ' ntml ',(ntml(i),i=1,n_dp)
!   write(6,*) ' ntpar ',(ntpar(i),i=1,n_dp)
!   write(6,*) ' dts_ntpar ',(dts_ntpar(i),i=1,n_dp)
!   write(6,*) ' cape_below_fr ',(cape_below_fr(i),i=1,n_dp)
!   write(6,*) ' cape_whole_layer ',(cape_whole_layer(i),i=1,n_dp)
!   write(6,*) ' cin ',(cin(i),i=1,n_dp)
!   write(6,*) ' ql_ad ',(ql_ad(i),i=1,n_dp)
!   write(6,*) ' h_ad ',(h_ad(i),i=1,n_dp)
!   write(6,*) ' diffmax ',(diffmax(i),i=1,n_dp)
!   write(6,*) ' pnb ',(pnb(i),i=1,n_dp)
!   DO k=1,20
!     write(6,*) ' qsat_moist ',k,(qsat_moist_ad(i,k),i=1,n_dp)
!     write(6,*) ' storethvp ',k,(storethvp(i,k),i=1,n_dp)     
!   END DO

!---------------------------------------------------
! Work out heights of convection
! Also find the maximum ntpar value so that loops over height can be shortened
ntparmax = 0
DO i_dp = 1,n_dp
  ztop(i_dp) = z_theta(i_dp,dts_ntpar(i_dp))

! Need cloud base and top level numbers output incase cca3d called in mid-level
! scheme 
  icct(i_dp) = dts_ntpar(i_dp)
  iccb(i_dp) = ntml(i_dp)

  IF(dts_ntpar(i_dp) > ntparmax) ntparmax = dts_ntpar(i_dp)
          
END DO
        
ktop(:)  = dts_ntpar(:) ! likely not to be necessary
ktopm(:) = dts_ntpar(:)
ktopp(:) = dts_ntpar(:)+1


l_increase = .true. ! z_theta increases with height
         
!DEPENDS ON: dts_locate_closest_levels
CALL dts_locate_closest_levels(n_dp,nlev,n_dp,n_dp,l_increase, & 
                               z_theta,h_ad,k_ad,k_adm,k_adp)

! Now send up a parcel from the freezing level

starting_heights(:) = zfr(:) ! m
th_excess = 0.5 ! K this is based on LES simulations theta
! excess at freezing level (see frsx/Wave/thvpert.pro)
icall=2

! DEPENDS ON: dts_cape
CALL dts_cape(n_dp,nlev, icall,                                            &
                 q,theta,thetav,qcl,qcf,qse,p_layer_centres,               &
                p_layer_boundaries,exner_layer_centres,z_theta,            &
                rho_theta,zlcl_dum,zlcl,klcl,freeze_lev,starting_heights,  &
                th_excess,ntpar,dts_ntpardum,                              &
                cape_below_frdum,cape_above_fr,cinabove,                   &
                qsat_moist_addum,ql_addum,h_addum,diffmaxdum,              &
                pnbabove,storethvp_upper)

!
! final code
!Calculate a variety of scaling parameters:
! mb, mfr, wcld, wfr, cstar, qrainmax,qsatsurf,qvatfr,sigma

!DEPENDS ON: dts_flux_par
CALL dts_flux_par(n_dp,nlev,klcl,freeze_lev,km40,land_mask,                &
                 q,qse,temperature,dr_across_th,                           &
                 zfr,w_max,zlcl,ztop,timestep,wstar,cape_below_fr,         &
                 cape_above_fr,cape_whole_layer,ql_ad,h_ad,                &
                 mb,mfr,w2lcl,wcld,wfr,wall,cstar,                         &
                 qrainmax,qsatsurf,qvatfr,sigma,scalefac1,scalefac2,scalefac3)

DO i=1,n_dp
  mb_deep(i) = mb(i) *g *rho(i,ntml(i))
END DO  

!=============================================================

! Work out what kind of convection it is -- ie where is the freezing
! level relative to h_ad, zlcl
! Not currently used, but may well need 
!DEPENDS ON: dts_conv_classify
CALL dts_conv_classify(n_dp,klclm,klclp,kfrm,kfrp,k_adm,k_adp              &
                      ,ktopm,ktopp,zlcl,zfr,h_ad,ztop                      &
                      ,iconvclass)


!================================================================
! FLUX CALCULATION SECTION
!================================================================

        
! Calculation of the velocity variance and mass flux         
! There are two methods in this routine -- under a switch: iusetrans
!DEPENDS ON: dts_w_variance
CALL dts_w_variance(iconvclass,n_dp,nlev,dts_ntpar,ntparmax,               &
                    z_theta,z_rho,rho,rho_theta,zlcl,zfr,ztop,             &
                    q,qse, diffmax,mb,mfr,w2lcl,wcld,wfr,wall,             &
                    scalefac1,scalefac2,scalefac3,thetav,storethvp,        &
                    ww,wwrho,massfl,massfl_rho)

up_flux(:,:) = rho_theta(:,:)*massfl(:,:)


! First work out the moist adiabatic thetav gradient along with the
! actual thetav gradient
!DEPENDS ON: dts_dthvdz
CALL dts_dthvdz(n_dp,nlev,dts_ntpar,ntparmax                               &
                ,qse,exner_layer_centres,z_theta,dr_across_rh,thetav       &
                ,storethvp,dthvdz,dthvdz_m,dqsedz)

!DEPENDS ON: dts_wthv
CALL dts_wthv(n_dp,nlev,ntparmax,dts_ntpar,                                &
              z_theta,z_rho,dr_across_rh,dr_across_th,rho,rho_theta,       &
              ztop,zlcl,thetav,dthvdz,dthvdz_m,wwrho,timestep,             &
              wthv,dthvdt_flux)


         
! Add the theta_v increment to theta_v
!DEPENDS ON: dts_update
CALL dts_update(n_dp,nlev,dts_ntpar,ntparmax,                              &
                timestep,thetav,dthvdt_flux,thetav_withflux)



!----------------------------------------------
! Added 8/5/09
! This is a very simple approach in which the water vapour flux is
! calculated explicitly rather than derived from the wthv and
! wmse fluxes. This is an attempt to get around small errors
! in the wthv and wmse fluxes leading to large, unphysical
! increments in the water vapour field
! iusegradq is recommended for the time being

IF(iusegradq == 1) THEN 
           
! Uses a very simple gradient formula for the water vapour flux, and
! derive theta from difference between mse and q flux
  q_withflux(:,:) = q(:,:)


  !DEPENDS ON: dts_qflux
  CALL dts_qflux(n_dp,nlev,dts_ntpar,ntparmax,z_theta,z_rho,rho            &
                  ,rho_theta,dr_across_rh,dr_across_th                     &
                  ,wwrho,massfl,massfl_rho,wstar,theta                     &
                  ,q,qsewat,zlcl,ztop,wq0,w2lcl,                           &
                 timestep,wqv,q_withflux)

  theta_withflux(:,:) = thetav_withflux(:,:)/(1.0+c_virtual*q_withflux(:,:))

ELSE
  !DEPENDS ON: dts_mseflux
  CALL dts_mseflux(n_dp,nlev,dts_ntpar,ntparmax,                           &
                 z_theta,z_rho,rho,rho_theta,dr_across_rh,dr_across_th,    &
                 wwrho,                                                    &
                 wstar,theta,q,mse,zlcl,klcl,ztop,wq0,wth0,timestep,       &
                 wmse,mse_withflux,h1all)
        

  !DEPENDS ON: dts_deduce_thetaandqv
  CALL dts_deduce_thetaandqv(n_dp,nlev,k_ad,q,mse,thetav_withflux,         & 
                 qhyd,mse_withflux,theta,                                  &
                 wmse,wthv,wthv_smth,                                      &
                 q_withflux,theta_withflux)
         

END IF ! iusegradq = 1 or not 1

! Move tracers around in a similar way to mse
! NB this is a guess about how they would behave, and could be improved
! on

IF(l_tracer) THEN

  CALL dts_tracerflux(n_dp,nlev,trlev,ntra,dts_ntpar,ntparmax,                &
                      z_theta,z_rho,dr_across_rh,dr_across_th,rho,            &
                      rho_theta,wwrho,wstar,tracer,zlcl,ztop,timestep,        &
                      dtrabydt)
END IF


! ------------------------------------------------------------------------
! Now work backwards to deduce what the increment due to the fluxes
! should be (this will then have the microphysics terms added to them)
! Includes boundary layer updates

DO k=1,ntparmax
  DO i_dp=1,n_dp
    IF(k <= dts_ntpar(i_dp)) THEN
      dthbydt(i_dp,k)=(theta_withflux(i_dp,k)-theta(i_dp,k))*rtimestep

      ! only apply subsidence term above the lcl
      dqbydt(i_dp,k)=(q_withflux(i_dp,k)-q(i_dp,k))*rtimestep
    END IF
  END DO
END DO

       
!================================================================
! MICROPHYSICS SECTION
!================================================================
! Increments due to phase changes
         
IF(icallmicro == 1) THEN 
           
  !DEPENDS ON: dts_cond_and_dep
  CALL dts_cond_and_dep(n_dp,nlev,ntparmax,dts_ntpar                         &
                        ,zlcl,zfr,rho,rho_theta,dr_across_rh,dr_across_th    &
                        ,qsat_moist_ad,qse,dqsedz,q                          &
                        ,massfl_rho,temperature,z_theta,z_rho                &
                        ,condensation,deposition,condint,depint)
         
! increment theta and q fields for condensation and deposition
  DO i_dp=1,n_dp
    DO k=1,ntparmax
      IF(z_theta(i_dp,k) >= zlcl(i_dp) .and. k <= dts_ntpar(i_dp)) THEN
        dthbydt(i_dp,k)= dthbydt(i_dp,k) +                                   &
               (condensation(i_dp,k)*lc/cp + deposition(i_dp,k)*(lc+lf)/cp)  &
                                 / exner_layer_centres(i_dp,k)
        dqbydt(i_dp,k) = dqbydt(i_dp,k) -                                    &
                            (condensation(i_dp,k)+deposition(i_dp,k))
      END IF
    END DO
  END DO

! -----------------------------------------------------------------
!nb need to think carefully whether the freezing term is required
! given the deposition term above), currently switched off
! inside routine  i.e. freeze=0.0, freezint=0.0 in routine

  freeze(:,:) = 0.0 
  freezint(:) = 0.0

! Don't call routine as currently does nothing will save CPU
!  !DEPENDS ON: dts_freeze
!  CALL dts_freeze(nlev,n_dp,dts_ntpar,ntparmax,z_theta,zlcl,klcl,temperature &
!                ,qrainmax,rho_theta,z_rho,freeze,freezint)
! update for freezing ! COMMENTED out as freeze is currently zero
!  DO k=1,ntparmax
!    DO i_dp=1,n_dp
!      IF(k <= dts_ntpar(i_dp)) THEN
!        dthbydt(i_dp,k)=dthbydt(i_dp,k) +                                  &
!                       (freeze(i_dp,k))*lf/cp/exner_layer_centres(i_dp,k)
!      END IF
!    END DO
!  END DO

  ! Now call the routine that evaluates the PC2 terms
  dbcfbydt(:,:) = 0.0
  ccw(:,:) = 0.0
  cclwp(:) = 0.0
                  
  !DEPENDS ON: dts_pc2
  CALL dts_pc2(n_dp,nlev,dts_ntpar,                                           &
               l_calc_dxek,                                                   &
               timestep,z_theta,zlcl,h_ad,ztop,rho,rho_theta,dr_across_rh,    &
               dr_across_th,temperature,qse,up_flux,q,                        &
               qcl,qcf,cf_frozen,cf_liquid,qcl_plume,qcf_plume,               &
               dcffbydt,dcflbydt,dqcfbydt,dqclbydt,                           &
               depint,freezint,condint,qclint,qcfint)


!  IF (.NOT. l_calc_dxek) THEN   ! not PC2
  ! set CCW to plume values ?

    DO k=1,nlev
      DO i_dp=1,n_dp
        IF (k <= dts_ntpar(i_dp)) THEN
          ccw(i_dp,k) = qcl_plume(i_dp,k)+qcf_plume(i_dp,k)
          delp = p_layer_boundaries(i_dp,k-1) - p_layer_boundaries(i_dp,k)
          cclwp(i_dp) = cclwp(i_dp) + ccw(i_dp,k)*delp/g
        END IF
      END DO 
    END DO 
!  END IF   
  ! ------------------------------------------------------------------------
  !DEPENDS ON: dts_sublimation
  CALL dts_sublimation(n_dp,nlev,dts_ntpar,ntparmax,                        &
                       dr_across_th,rho_theta,q,temperature,qse,            &
                       deposition,depint,freezint,qcfint,                   & 
                       sublimation,subint)

  ! update for sublimation
  DO k=1,ntparmax
    DO i_dp=1,n_dp
      dthbydt(i_dp,k)=dthbydt(i_dp,k) -                                     &
                          sublimation(i_dp,k)*(lc+lf)/cp/                   &
                                               exner_layer_centres(i_dp,k)
      dqbydt(i_dp,k) = dqbydt(i_dp,k) + sublimation(i_dp,k)
    END DO
  END DO

  ! Melt the falling component according to a prescribed profile as a fn
  ! of temperature
  !DEPENDS ON: dts_melt
  CALL dts_melt(nlev,n_dp,dts_ntpar,ntparmax,                               &
                rho_theta,dr_across_th,temperature,                         &
                depint,subint,freezint,qcfint,                              &
                melt,meltint)

  ! update for melting
  DO k=1,ntparmax
    DO i_dp=1,n_dp
      IF(k <= dts_ntpar(i_dp)) THEN
        dthbydt(i_dp,k)=dthbydt(i_dp,k) +                                   &
                       (-melt(i_dp,k))*lf/cp/exner_layer_centres(i_dp,k)
      END IF
    END DO
  END DO

  !DEPENDS ON: dts_evaporation
  CALL dts_evaporation(nlev,n_dp,dts_ntpar,ntparmax                         &
                       ,temperature,rho_theta,z_theta,dr_across_rh          &
                       ,dr_across_th,qsat_moist_ad,q,mb,wstar,wcld          &
                       ,ztop,zfr,zlcl,melt,freeze,condensation              &
                       ,qclint,condint,freezint,meltint                     &
                       ,evaporation,evapint)

  ! update evaporation
  DO k=1,ntparmax
    DO i_dp=1,n_dp
      IF(z_theta(i_dp,k) >= zlcl(i_dp) .and. k <= dts_ntpar(i_dp)) THEN
        dthbydt(i_dp,k)=dthbydt(i_dp,k) -                                   &
                       evaporation(i_dp,k)*(lc/cp)/exner_layer_centres(i_dp,k)
        dqbydt(i_dp,k) = dqbydt(i_dp,k) + evaporation(i_dp,k)
      END IF
    END DO
  END DO
      
! ---------------------------------------------------------------------
! 12. 2/10/07 calculate a rain production term
  rain_3d(:,:) = 0.0
  snow_3d(:,:) = 0.0
  !DEPENDS ON: dts_rainprod
  CALL dts_rainprod(n_dp,nlev,z_rho,z_theta,dr_across_rh,dr_across_th      &
                   ,rho,rho_theta,temperature,q,qse,zlcl,zfr               &
                   ,condensation,evaporation,melt                          &
                   ,freeze,dqclbydt,dqcfbydt,deposition,sublimation        &
                   ,rainprod,snowprod,revp,revpint,rainrate,snowrate)

         

! update for rain evaporation
  dqbydt(:,:) = dqbydt(:,:) + revp(:,:) 
  dthbydt(:,:)= dthbydt(:,:)- revp(:,:)*(lc/cp)             
  rain(:) = rainrate(:)
  snow(:) = snowrate(:)

END IF ! icallmicro==1

!--------------------------------------------------------------------------
! CMT section
!--------------------------------------------------------------------------
! Scale mass flux so units correct for output & for CMT scheme i.e. Pa/s

  up_flux(:,:) = up_flux(:,:)*g


  kterm(:) = dts_ntpar(:)    ! Top of deep convection

  IF (l_mom) THEN

! As input need kterm - level deep convection terminates

    ncmt = 0          ! number of deep columns to deep CMT on  
    DO i=1,n_dp

! check deep convection over at least 3 levels ?
      If (wall(i)  > 0.0 .and. (kterm(i) - ntml(i)) > 1) THEN
        ncmt = ncmt+1  
       cu_term(ncmt) = i         ! index in case some columns diagnosed
                                 ! as deep convective DOn't convect.
                                 ! Assuming all convect for this scheme  
      END IF
    END DO  

    CALL deep_turb_cmt( n_dp, ncmt, nlev, deep_cmt_opt,                 &
                        ntml, kterm, cu_term, freeze_lev,               &
                        timestep,                                       &
                        uw0, vw0, mb, wall, wstar, zlcl_uv,             &
                        up_flux,                                        &
                        r_rho, r_theta, z_rho, z_theta, rho, rho_theta, &
                        r2rho, r2rho_th, dr_across_th, dr_across_rh,    &
                        u, v,                                           &
                        dubydt, dvbydt, uw_deep, vw_deep)

  END IF        ! test on l_mom

!--------------------------------------------------------------------------
! Check for negative q after increments applied i.e. removing all moisture
! Same approach as deep mass flux scheme. Note without this q was going
! negative at higher model levels.
!--------------------------------------------------------------------------
! check on shallow deep
  DO i = 1,n_dp
    IF (dts_ntpar(i) < 5) THEN
      write(6,*) ' deep turb alison ',CALL_number,i,' ntml ntpar ', &
                     klcl(i),dts_ntpar(i),ntml(i),ntpar(i)

    END IF
  END DO

  DO i = 1,n_dp
    qMinInColumn(i) = q(i,nlev)
  END DO
  DO k = 1,nlev-1
    DO i = 1,n_dp
      IF (q(i,k)  <   qMinInColumn(i)) THEN
        qMinInColumn(i) = q(i,k)
      END IF
    END DO
  END DO

! Ensure Q DOes not go below global allowed minimum (QMIN)

  DO i = 1,n_dp
    qMinInColumn(i)=MAX(qmin,qMinInColumn(i))
  END DO

! Apply an artificial upwards flux from k-1 level to ensure Q
! remains above minimum value in the column.

  DO k = nlev,2,-1
    DO i = 1,n_dp
      IF (dqbydt(i,k) /= 0.0) THEN
        temp1(i)=q(i,k) + dqbydt(i,k) * timestep
        IF (temp1(i)  <   qMinInColumn(i)) THEN

          dqbydt(i,k-1) = dqbydt(i,k-1) -                              &
                   ((qMinInColumn(i) - q(i,k)) / timestep-dqbydt(i,k)) &
                    * (r2rho_th(i,k)*dr_across_th(i,k))                &
                    / (r2rho_th(i,k-1)*dr_across_th(i,k-1))

          dqbydt(i,k) = (qMinInColumn(i) - q(i,k)) / timestep
! Only warning print statement not DOing any thing
!          write(6,*) ' negative q deep',i,k,temp1(i),dqbydt(i,k)

        END IF
      END IF
    END DO ! n_dp loop
  END DO  ! nlev

! Check negative q bottom level.

  k=1
    DO i = 1,n_dp
      temp1(i)=q(i,k) + dqbydt(i,k) * timestep
      IF (temp1(i)  <   qMinInColumn(i)) THEN
! This is a warning message - no action possible as no layer below this

        write(6,*) ' negative q deep',i,k,kterm(i),temp1(i),dqbydt(i,k),&
                   w_max(i),mb(i),rain(i),snow(i)
        write(6,*) ' dq ',(dqbydt(i,j),j=1,nlev)
        write(6,*) ' dqcl ', (dqclbydt(i,j),j=1,nlev)
        write(6,*) ' dqcf ', (dqcfbydt(i,j),j=1,nlev)

      END IF
    END DO ! n_dp loop
!--------------------------------------------------------------------------
! Check conservation of moisture by scheme
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! Convective cloud amount section 
! At present the scheme does not return a convective cloud amount.
! The easiest approach to take is to estimate the convective cloud 
! in the same way as the old mass flux scheme.
!--------------------------------------------------------------------------
! cca already initialised to zero

! Set cca_2d based on convective precipitation

  cca_2d(:) = 0.0
  DO i=1, n_dp

    IF ((rain(i) + snow(i)) > 0.0) THEN
      IF(land_mask(i)) THEN
        cca_2d(i)  = a_land + b_land * ALOG(86400.0 * (rain(i)+snow(i)))
      ELSE
        cca_2d(i)  = a_sea  + b_sea  * ALOG(86400.0 * (rain(i)+snow(i)))
      END IF 
    END IF
    IF (icct(i) > iccb(i)) THEN 

      cca_2d(i) = MAX(2.0E-5,cca_2d(i))   ! sets a minimum cloud in case no 
                                          ! precip
    ELSE     ! case of failed deep convection

      iccb(i) = 0
      icct(i) = 0
      cca_2d(i) = 0.0

    END IF  

  END DO 


!NB Probably should not allow a choice here i.e. only have l_ccrad

  cca(:,:) = 0.0         ! initialise to zero

  IF (l_ccrad) THEN


    !---------------------------------------------------------------------
    ! Apply CCA_2D to 3d cloud profile
    !---------------------------------------------------------------------
    ! Apply anvil scheme to deep cloud
    ! DEPENDS ON: CALC_3D_CCA 
      CALL calc_3d_cca(n_dp, n_dp, nlev, n_cca_lev, nbl, iccb, icct     & 
             , p_layer_boundaries, freeze_lev, cca_2d, cca              &
             , z_theta, z_rho)

      ! NOTE: iccb, icct are layer centres (theta levels) at this
      !        point.



  ELSE        ! Not CCRAD  and not PC2 require a cloud

    IF (.NOT.l_calc_dxek) THEN

      ! Assume 3D anvil cloud required. 

      ! DEPENDS ON: CALC_3D_CCA 
      CALL calc_3d_cca(n_dp,n_dp,nlev,n_cca_lev,nbl,iccb,icct,          &
                       p_layer_boundaries,freeze_lev,                   &
                       cca_2d,cca,z_theta,z_rho)

    END IF      

  END IF      ! l_ccrad

!--------------------------------------------------------------------------

  IF (lhook) CALL dr_hook('DEEP_TURB_CONV',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE deep_turb_conv
