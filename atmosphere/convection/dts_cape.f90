! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!
!+  Measure cape and related quantities for deep turbulence scheme
!
SUBROUTINE dts_cape(n_dp,nlev,icall,                                         &
                 q,theta,thetav,qcl,qcf,qse,                                 &
           p_layer_centres,p_layer_boundaries,exner_layer_centres,           &
           z_theta,rho_theta,zlcl_cape,zlcl,klcl,freeze_lev,starting_heights,&
           th_excess,ntpar,dts_ntpar,                                        &
           cape_below_fr,cape_whole_layer,                                   &
           cin,qsat_moist_ad,ql_ad,h_ad,diffmax,pnb,storethvp)

USE atmos_constants_mod, ONLY: r, cp, kappa, c_virtual, rv

USE water_constants_mod, ONLY: lc, lf, tm
USE dts_cntl_mod, ONLY:                                                        &
    dts_qfac

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!-----------------------------------------------------------------------------
!
! Purpose: calculates cape and related quantities for deep turb scheme
!
! Called by: deep_turb_conv
!
! What this subroutine does:
! --------------------------
! -calculates cape and cin via an undilute parcel ascent
! -this is divided into: cape up to the freezing level, and for the
! whole layer
! -retains the option to include or lose water loading

! Inputs: 
! -------
! theta,q,qcl,qcf,p_layer_centres
! p_layer_boundaries,exner_layer_centres,z_theta,rho_theta,freeze_lev

! Outputs: 
! --------
! cape_below_fr
! cape_whole_layer
! cin

! Potential Weaknesses:
! ---------------------
! Variants that could be tried: 
! 1) with/without water loading! 
! 2) allowing different amounts of hydrometeor to precipitate out! 
! 3) different starting level heights

! Other information: 
! ------------------
! - this sub. is called from the SCALING PARAMETER section of the
! deep_turb_conv routine
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP 3 programming standards version vn8.2.
!------------------------------------------------------------------------------

INTEGER, INTENT(in) ::  &
  nlev                  & ! No. of model layers
 ,n_dp                  & ! No. convecting points
 ,freeze_lev(n_dp)      & ! Index Level of freezing 
 ,klcl(n_dp)            & ! Index Level lcl
 ,ntpar(n_dp)

INTEGER, INTENT(in) ::  &
  icall                   ! Call number temporary for debugging    

REAL, INTENT(IN) ::      &
  q(n_dp,nlev)           &  ! Model mixing ratio  (kg/kg)
, theta(n_dp,nlev)       &  ! Model potential temperature (K)
, qcl(n_dp,nlev)         &  ! L-S Liq condensate mix ratio (kg/kg)
, qcf(n_dp,nlev)         &  ! L-S Ice condensate mix ratio (kg/kg)
, starting_heights(n_dp) &  ! starting height of ascent in m
, th_excess              &  ! initial th perturbation 
, qse(n_dp,nlev)            ! Saturation mixing ratio of environment (kg/kg)

REAL, INTENT(IN) ::                &   
  p_layer_centres(n_dp,0:nlev)     & ! Pressure(Pa) at theta levels
, p_layer_boundaries(n_dp,0:nlev)  & ! Pressure(Pa) on rho levels
, exner_layer_centres(n_dp,0:nlev) & ! Exner pressure theta levels
, z_theta(n_dp,nlev)               & ! height of theta levels(m)
, rho_theta(n_dp,nlev)             & ! density on theta levels (kg/m3)
, zlcl(n_dp)                       & ! lifting condensation level (m)
, thetav(n_dp,nlev)                  ! env virtual pot temp   (K)   

INTEGER, INTENT(OUT) ::   &
  dts_ntpar(n_dp)           ! ntpar based on this parcel ascent

REAL, INTENT(OUT) ::           &
  cape_below_fr(n_dp)          & ! cape up to freezing level  (J/kg) 
 ,cape_whole_layer(n_dp)       & ! cape up to level of neutral buoy (J/kg)
 ,cin(n_dp)                    & ! convective inhibition (J/kg)
 ,ql_ad(n_dp)                  & ! adiabatic liquid water content at
                                 ! level of max buoy excess  (kg/kg)
 ,h_ad(n_dp)                   & ! height of max buoy excess  (m) 
 ,qsat_moist_ad(n_dp,nlev)     & ! qsat along parcel ascent (theta levs)
 ,storethvp(n_dp,nlev)         & ! parcel virtual potential temperature
 ,diffmax(n_dp)                & ! maximum buoyancy excess (thv')
 ,pnb(n_dp)                    & ! pressure at level of neutral buoyancy (Pa) 
 ,zlcl_cape(n_dp)                ! give lcl from cape routine
      
!-----------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------------
INTEGER ::    &
  i_dp,k,     & ! loop counters
  n_itime       ! number of points with itime=1

INTEGER ::        &
  klev(n_dp)      & ! levels to start parcel ascent
 ,klevm(n_dp)     &
 ,klevp(n_dp)     &
 ,klev_min        & ! minimum klev across all points
 ,ktop_max        & ! maximum ktop across all points
 ,kpos(n_dp)      & ! level at which buoyancy excess is maximum
 ,itime(n_dp)     & ! counter to determine how 
 ,npnts           & ! used in calls to qsat_mix (set to 1 here)
 ,one

LOGICAL ::    &
  lq_mix      & ! used in qsat_mix: uses specific hum if true
 ,l_itime     & ! no itime values still zero
 ,l_increase    !

REAL ::               &
  delp(n_dp)          & ! pressure diff. between one level and next (Pa)
 ,delp_cen(n_dp)      &
 ,tparcel(n_dp)       & ! temperature of parcel  (K)
 ,thparcel(n_dp)      & ! potential temperature of parcel  (K)
 ,p(n_dp)             & ! pressure of parcel  (Pa)
 ,qvbot(n_dp)         & ! parcel water vapour content at start  (kg/kg)
 ,qlbot(n_dp)         & ! parcel liquid water content at start  (kg/kg)          
 ,qlval(n_dp)         & !liquid water content inside parcel
 ,cape1(n_dp)         & !cape with water loading / J/kg
 ,cape2(n_dp)         & !cape without water loading / J/kg
 ,cin1(n_dp)          & !cin with water loading / J/kg
 ,cin2(n_dp)          & !cin without water loading / J/kg
 ,capelower1(n_dp)    & !cape up to freezing level with water loading
 ,capelower2(n_dp)    & !cape up to freezing level without water loading
 ,smallcape           & ! minimum value of cape before
                        ! negative buoyancy counts as level of neutral buoyancy
                        ! rather than cin
 ,qvparcel(n_dp)        ! water vapour mixing ratio in parcel

REAL ::                 &
  storeqlval(n_dp,nlev) & ! liquid water content along parcel ascent         
 ,diffth(n_dp,nlev)     & ! profile of buoyancy excess (thv')
 ,qcl_plume(n_dp,nlev)    ! Plume liquid water content 

REAL ::          &
  thvparcel1     & ! virtual potential temp of parcel w water loading
 ,thvparcel2     & ! virtual potential temp of parcel without water loading
 ,capecom1       & ! contribution to cape at a given level with water loading 
 ,capecom2       & ! contribution to cape at a given
                   ! level without water loading
 ,vlt            & ! latent heat of deposition or condensation
 ,rs1(n_dp)      & ! sat value on existing level
 ,rs2(n_dp)      & ! sat value on previous levels 
 ,t0,t1,t2       & ! local temp variables
 ,q0(n_dp)       & ! local qsat variables
 ,q1(n_dp)       & ! local qsat variables
 ,q2(n_dp)       & ! local qsat variables
 ,qnext(n_dp)    & ! local qsat variables
 ,alpha          & !
 ,beta           & !
 ,bb1            & !
 ,bb2            & ! local variables
 ,drs0           & !
 ,dt             & !
 ,dq             & !
 ,deltat         & ! local variables
 ,qlcrit         & !
 ,qlmin          & ! maximum liquid water content of parcel -- beyond
 ,tmpthvp(n_dp)  &
 ,smallpressure    ! minimum value of pressure -- set to 1 Pa
                   ! calc'd in dts_pc2 -- not v efficient!)
                   ! that rains out

! extra for SCM output to look at parcel ascent
REAL :: &
  qv_parc(n_dp,nlev)     &
 ,ql_parc(n_dp,nlev)     &
 ,t_parc(n_dp,nlev)


REAL ::     &
  rbb2      &  ! 1/bb2
, rdeltat   &  ! 1/deltat
, dp_o_p    &  ! dp /p
, gdz_o_thv    ! dp/(thv rho)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!======================================================================
 
! Initialise some variables
 IF (lhook) CALL dr_hook('DTS_CAPE',zhook_in,zhook_handle)
  one = 1
  n_itime = 0  
  l_itime = .true.   
  qlcrit = 1.e-3 ! kg/kg
  qlmin  = 2.e-4 ! kg/kg
  deltat   = 0.2 ! K iteration for temperature interval
  rdeltat  = 1.0/0.2 ! K iteration for temperature interval
    
  smallcape= 1.0 ! J/kg a small value for cape (somewhat arbitrary)
  smallpressure = 1.0 ! Pa
  npnts    = 1
  lq_mix   = .true. ! use specific humidity
  dts_ntpar(:) = 1
  itime(:) = 0 ! initialise counter
  cape1(:) = 0.0
  cape2(:) = 0.0
  cin1(:)  = 0.0
  cin2(:)  = 0.0
  qvparcel(:) = 0.0
  pnb(:)   = 0.0 !initialise level of neutral buoyancy
  qlval(:) = 0.0
  delp(:) = 0.0
  delp_cen(:) = 0.0
  qsat_moist_ad(:,:) = 0.0
  cape_below_fr(:) = 0.0
  capelower1(:) = 0.0
  capelower2(:) = 0.0
  cape_whole_layer(:) = 0.0
  storethvp(:,:) = 0.0
  tmpthvp(:) = 0.0
  zlcl_cape(:) = 0.0
  cin(:) = 0.0
  ql_ad(:) = 0.0
  h_ad(:) = 0.0
  kpos(:) = 0.5*(freeze_lev(:)+ntpar(:)) 
  ! something not too unreasonable, in
  ! case it fails to find a level of maximum buoyancy...
  qv_parc(:,:)=0.0
  ql_parc(:,:)=0.0
  t_parc(:,:) =0.0

! Now set up a liquid water content for plumes
DO k=1,nlev
  DO i_dp=1,n_dp
    qcl_plume(i_dp,k) = 0.5*qse(i_dp,k)
    IF(qcl_plume(i_dp,k) > qlcrit) THEN 
      qcl_plume(i_dp,k) = qlcrit
    ELSE IF(qcl_plume(i_dp,k) < qlmin) THEN 
      qcl_plume(i_dp,k) = qlmin
    END IF
  END DO
END DO

! Find the level corresponding to the starting height

 l_increase = .true. ! z_theta increases with height

!DEPENDS ON: dts_locate_closest_levels
call dts_locate_closest_levels(n_dp,nlev,n_dp,n_dp,l_increase,z_theta     &
                               ,starting_heights,klev,klevm,klevp)
      
! now find lowest value of klev = klev_min out of all of the n_dp
! convecting points
 
klev_min = nlev ! start off very high
DO i_dp = 1,n_dp
  IF(klev(i_dp) < klev_min) THEN
    klev_min = klev(i_dp)
  END IF 
END DO
      
ktop_max = nlev-1

!==========================================================================

! Assign initial thermodynamic values
      
DO i_dp = 1,n_dp
  thparcel(i_dp) = theta(i_dp,klev(i_dp))+th_excess
  tparcel(i_dp)  = thparcel(i_dp)*exner_layer_centres(i_dp,klev(i_dp))
  ! use large-scale value 
  qlbot(i_dp)    = qcl(i_dp,klev(i_dp))+qcf(i_dp,klev(i_dp))
  p(i_dp)        =  p_layer_centres(i_dp,klev(i_dp))
END DO ! i_dp
      
! find the saturation value of this initial parcel
! DEPENDS ON: qsat_mix
CALL qsat_mix(q0,tparcel,p,n_dp,lq_mix)
      
DO i_dp = 1,n_dp
  IF(z_theta(i_dp,klev(i_dp)) <= zlcl(i_dp)) THEN 
! if the parcel starts beneath the lcl, then give it the water vapour
            ! value for that level
! nb multiplying by an extra factor of 1.2
!            qvbot(i_dp) = q(i_dp,klev(i_dp))*1.2
! Sensitivity test
    qvbot(i_dp) = q(i_dp,klev(i_dp))*dts_qfac
  ELSE
! if the parcel starts above the lcl, THEN say that it's already saturated
! and give it the plume liquid value
    qvbot(i_dp) = q0(i_dp)
    qlbot(i_dp) = qcl_plume(i_dp,klev(i_dp)) ! i.e. 1 g/kg
  END IF
END DO ! i_dp


!==================================================================
! Start parcel ascent -- increasing height by vertical levels
!==================================================================

DO k = klev_min,ktop_max
       

! create a variable p to shorten lines of code
  p(:) =  p_layer_centres(:,k)
  delp(:) = p_layer_boundaries(:,k)-     &
                             p_layer_boundaries(:,k+1) ! positive quantity

  delp_cen(:) = p_layer_centres(:,k)-p_layer_centres(:,k+1)
! calculate how much water vapour would be needed at this T and p in
! order to be saturated
! for later moist adiabatic ascent, will need gradient as a function of T
  ! DEPENDS ON: qsat_mix
  CALL qsat_mix(q0,tparcel,p,n_dp,lq_mix)
 
  qvparcel(:) = q0(:) ! first suggest saturated

  DO i_dp = 1,n_dp 
          
    IF(k >= klev(i_dp) .and. pnb(i_dp) < smallpressure) THEN

! liquid water in parcel is total water minus saturation value
      qlval(i_dp) = qlbot(i_dp)+qvbot(i_dp)-qvparcel(i_dp) 

      IF (qlval(i_dp) < 0.0) THEN
        qlval(i_dp) = 0.0 ! ensure >= 0
      ELSE IF (qlval(i_dp) > qcl_plume(i_dp,k)) THEN
        qlval(i_dp) = qcl_plume(i_dp,k)
      END IF
! If the parcel was unsaturated, then its value stays that of the initial parcel
      IF(qlval(i_dp) < 1.0e-10) THEN
        qvparcel(i_dp)  = qvbot(i_dp)
        zlcl_cape(i_dp) = z_theta(i_dp,k)
      END IF

      ! evaluation of parcel virt.pot.temp keeping all water loading
      thvparcel1 = thparcel(i_dp)*(1.0+c_virtual*qvparcel(i_dp)-qlval(i_dp))

      ! evaluation of parcel virt.pot.temp losing all liquid water
      thvparcel2 = thparcel(i_dp)*(1.0+c_virtual*qvparcel(i_dp))

      ! contribution to cape from this level:

       gdz_o_thv =  delp(i_dp)/(thetav(i_dp,k)*rho_theta(i_dp,k))
       capecom1 = gdz_o_thv*(thvparcel1-thetav(i_dp,k))

       capecom2 = gdz_o_thv*(thvparcel2-thetav(i_dp,k))

       qv_parc(i_dp,k) = qvparcel(i_dp) 
               
       IF(capecom1 > 0.0) THEN
         cape1(i_dp)=cape1(i_dp)+capecom1 ! w water loading
       END IF
       IF(capecom2 > 0.0) THEN
         cape2(i_dp)=cape2(i_dp)+capecom2 ! no water loading
       END IF

       ! if negative and cape is still small, then count as cin 
       IF(capecom1 < 0.0 .and. cape1(i_dp) <= smallcape) THEN
         cin1(i_dp)=cin1(i_dp)+capecom1
       END IF    
       IF(capecom2 < 0.0 .and. cape2(i_dp) <= smallcape) THEN
         cin2(i_dp)=cin2(i_dp)+capecom2
       END IF    

! define level of neutral buoyancy as the first level at which
! negatively buoyant once cape is greater than a certain value

       IF(capecom2 < 0.0 .and. cape2(i_dp) > smallcape) THEN
         pnb(i_dp)=p(i_dp)
         dts_ntpar(i_dp) = k
       END IF
! store cape up to freezing level separately
       IF(k < freeze_lev(i_dp)) THEN 
         capelower1(i_dp) = cape1(i_dp)
         capelower2(i_dp) = cape2(i_dp)
       END IF
! store various profiles for determination of ql_ad,h_ad further down
       storethvp(i_dp,k) = thvparcel1 ! thparcel(i_dp) 
       qsat_moist_ad(i_dp,k) = q0(i_dp)
       storeqlval(i_dp,k)    = qlval(i_dp)

       ql_parc(i_dp,k) = qlval(i_dp)


     END IF !k >= klev(i_dp) .and. pnb(i_dp) == 0.0
   END DO ! i_dp
             
! determine qsat for the level below
! it's only called once, so having the if(itime==0) test saves time too...
!                 if(itime(i_dp) == 0) THEN 
!                    npnts = 1
!                   ! DEPENDS ON: qsat_mix
!                   CALL qsat_mix(rs1(i_dp),tparcel(i_dp)              &
!                          ,p_layer_centres(i_dp,k-1),npnts,lq_mix)
!                 ELSE
!                    rs1(i_dp) = rs2(i_dp)
!                 END IF
!              END IF !k >= klev(i_dp) .and. pnb(i_dp) == 0.0
!           END DO ! i_dp


! Not Sure this is right ? But avoids qsat in level loop (VERY EXPENSIVE)
! Still points where itime is zero 
   IF (l_itime) THEN    
     ! DEPENDS ON: qsat_mix
     CALL qsat_mix(rs1,tparcel,p_layer_centres(:,k-1),n_dp,lq_mix)
   END IF

   DO i_dp = 1,n_dp
     IF(itime(i_dp) /= 0) THEN 
       IF(k >= klev(i_dp) .and. pnb(i_dp) < smallpressure) THEN
         rs1(i_dp) = rs2(i_dp)

       END IF !k >= klev(i_dp) .and. pnb(i_dp) == 0.0
     END IF
   END DO ! i_dp

   ! DEPENDS ON: qsat_mix
   CALL qsat_mix(qnext,tparcel,p_layer_centres(:,k+1),n_dp,lq_mix)
                 
   ! DEPENDS ON: qsat_mix
   CALL qsat_mix(q1,tparcel-deltat,p_layer_centres(:,k+1),n_dp,lq_mix)

! break out of points loop to DO qsat calculations
   DO i_dp = 1,n_dp
     IF(k >= klev(i_dp) .and. pnb(i_dp) == 0.0) THEN

! determine a local gradient for the qsaturation
! (assume a local linear relationship between T and q for a given p)
 
      IF(qlval(i_dp) > 0.0) THEN ! if it's saturated then do moist ascent
         t0 = tparcel(i_dp)
!        t1 = t0 - deltat   ! not required as use (t0-t1)=deltat
                        
! choose between L_vap and L_fus tm = 273.15 K
         IF(tparcel(i_dp) >= tm) THEN 
            vlt = lc
         ELSE 
            vlt = lc+lf
         END IF
         ! determine gradient of q vs T line  
         alpha = (qnext(i_dp)-q1(i_dp))*rdeltat
         beta = q0(i_dp)-alpha*t0

         dp_o_p = delp_cen(i_dp)/p(i_dp)
                    
         ! Calculation of bb1 including 2nd order term
         bb1 =  -(r*tparcel(i_dp)*dp_o_p)*(1.0-(kappa-1.0)*dp_o_p*0.5)   

!        bb2 = 1.0 + q0(i_dp)
         rbb2 = 1.0/(1.0 + q0(i_dp))
         drs0 = qnext(i_dp)-q0(i_dp)   !   q0(i_dp) - rs1(i_dp)
  
!        dt = (bb1-vlt*drs0/bb2)/(cp+alpha*vlt/bb2)
         dt = (bb1-vlt*drs0*rbb2)/(cp+alpha*vlt*rbb2)
                   
         dq = drs0 + alpha*dt
     
         tparcel(i_dp) = tparcel(i_dp) + dt
         itime(i_dp) = itime(i_dp) + 1
         n_itime = n_itime + 1

! this is the potential temperature at the next height up
         thparcel(i_dp) = tparcel(i_dp)/exner_layer_centres(i_dp,k+1)

       ELSE

!deduce T given new pressure 
         tparcel(i_dp) = thparcel(i_dp)*exner_layer_centres(i_dp,k+1)
       END IF

      t_parc(i_dp,k) = tparcel(i_dp)

    END IF ! (k >= klev(i_dp) .and. pnb(i_dp) == 0.0) 
  END DO ! i_dp

! Reevaluate qsat, this time with new tparcel, but old pressure
  !DEPENDS ON: qsat_mix
  CALL qsat_mix(rs2,tparcel,p,n_dp,lq_mix)


  IF (n_itime == n_dp) THEN
    l_itime = .false.          
  END IF  

END DO ! k
!==============================================================
! End of parcel ascent
!==============================================================
      
! Find the level at which the parcel buoyancy is greatest
! and corresponding level, kpos
  
diffmax(:) = 0.0
DO k=1,nlev
  DO i_dp = 1,n_dp

    ! Below or equal to starting level set to parcel thv (klev+1)
    IF(k <= klev(i_dp)) THEN 
      storethvp(i_dp,k) = storethvp(i_dp,klev(i_dp)+1)

    ! within parcel ascent and below top         
    ELSE if(k >= klev(i_dp) .AND. p_layer_centres(i_dp,k) >= pnb(i_dp)) THEN
               
      tmpthvp(i_dp) = storethvp(i_dp,k) ! just to store value at top
    END IF

    ! Above top or crazy thv value reset to value at top of ascent
    IF(p_layer_centres(i_dp,k) < pnb(i_dp) .OR. &
                                   storethvp(i_dp,k) < 10.0) THEN 
      storethvp(i_dp,k) = tmpthvp(i_dp) ! ensure a value is given
                                        ! right the way up
    END IF

          
    diffth(i_dp,k) = storethvp(i_dp,k)-thetav(i_dp,k)

    IF(diffth(i_dp,k) > diffmax(i_dp)) THEN 
      diffmax(i_dp) = diffth(i_dp,k)
      kpos(i_dp) = k
    END IF

  END DO
END DO ! k
      
      
! now store adiabatic water content, the height of the maximum buoyancy
! excess, and various cape, cin parameters
DO i_dp = 1,n_dp
  ql_ad(i_dp) = storeqlval(i_dp,kpos(i_dp))
  h_ad(i_dp) = z_theta(i_dp,kpos(i_dp))

!nb may change mind about whether or not to include water loading

  cin(i_dp) = cin2(i_dp)                ! without water loading
  cape_whole_layer(i_dp) = cape1(i_dp)  ! nb to decide about water loading
  cape_below_fr(i_dp)= capelower1(i_dp) ! nb to decide about water loading

END DO
IF (lhook) CALL dr_hook('DTS_CAPE',zhook_out,zhook_handle)
RETURN

      
END subroutine dts_cape
