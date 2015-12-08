! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Completes lifting of the parcel from layer k to k+1
!
SUBROUTINE convec2_4a5a(npnts,np_full,nlev,ntra,k,ad_on,sdet_on,new_termc, &
                        start_lev,                                         &
                        l_tracer,l_mom_gk,l_q_interact,l_calc_dxek,        &
                        l_shallow,l_mid,cumulus,                           &
                        bwk,bwkp1,bgmkp1,bland,blowst,                     &
                        timestep,                                          &
                        thek,thekp1,qek,qekp1,qclek,qclekp1,qcfek,qcfekp1, &
                        cflek,cflekp1,cffek,cffekp1,bcfek,bcfekp1,         &
                        uek,uekp1,vek,vekp1,                               &
                        traek,traekp1,trapkp1,                             &
                        qsekp1,dqskp1,pstar,thpkp1,qpkp1,upkp1,vpkp1,      &
                        xsqkp1,rbuoy,qsek,dqsk,thpi,qpi,expi,              &
                        ekp14,ekp34,amdetk, pk,pkp1,exk,exkp1,delexkp1,    &
                        delpk,delpkp1,delp_uv_k, delp_uv_kp1,              &
                        t1_sd,q1_sd,thpixs_v,qpixs_v,xsbmin_v,rbuoy_p_old, &
                        ! In/out
                        binit,bgmk,                                        &
                        thpk,qpk,qclpk,qcfpk,qclpkp1,qcfpkp1,              &
                        upk,vpk,trapk,xpk,flxk,                            &
                        dthek,dqek,dqclek,dqcfek,dcflek,dcffek,dbcfek,     &
                        duek,dvek,dtraek,                                  &
                        tcw,depth,cclwp,cape,dcpbydt,eflux_u_ud,eflux_v_ud,&
                        ! Out
                        iccb,icct,lcbase,lctop,                            &
                        bterm,                                             &
                        prekp1,dthekp1,dqekp1,dqclekp1,dqcfekp1,           &
                        dcflekp1,dcffekp1,dbcfekp1,duekp1,dvekp1,dtraekp1, &
                        cca,ccw,lcca,deltak,flxkp12,max_cfl_c,relh,        &
                        dptot,rbuoy_p_here,the_here,thp_here,qe_here,qp_here)


USE cv_run_mod, ONLY:                                                 &
    cape_opt, cnv_wat_load_opt

USE atmos_constants_mod, ONLY:                                        &
          c_virtual, r

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Completes lifting of the parcel from layer k to k+1.
!
!   Calls subroutines parcel and environ
!  
!   Subroutine parcel calculates an initial mass flux, carries out the
!   detrainment calculation, tests to see if convection is termintating 
!   and calculates the precipitation rate from layer k+1.
!
!   Subroutine environ calculates the effect of convection upon the 
!   large-scale atmosphere.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  npnts                & ! Number of points
 ,np_full              & ! Full vector length
 ,nlev                 & ! Number of model levels for calculations
 ,ntra                 & ! Number of tracer variables
 ,k                    & ! present model layer
 ,ad_on                & ! Flag for adaptive detrainment 
 ,sdet_on              & ! Flag smoothed forced detrainment
 ,new_termc              ! Flag for simplified termination of convection

INTEGER, INTENT(IN) :: &
  start_lev(npnts)       ! Level at which convection is initiated

LOGICAL, INTENT(IN) :: &
  l_tracer             & ! Switch for tracers
 ,l_mom_gk             & ! Gregory-Kershaw convective momentum transport switch 
 ,l_calc_dxek          & ! Switch for calculating condensate inc.
 ,l_q_interact           ! Switch allows overwriting of parcel
                         ! variables when calculating condensate
                         ! increments (will alter results).
LOGICAL, INTENT(IN) ::  &
  l_shallow(npnts)      & ! true in shallow convection
 ,l_mid(npnts)          & ! true if mid convection
 ,cumulus(npnts)          ! true if convection from the boundary layer

LOGICAL, INTENT(IN) :: &
  bwkp1(npnts)         & ! mask for whether condensate is liquid in layer k+1
 ,bwk(npnts)           & ! mask for whether condensate is liquid in layer k
 ,bland(npnts)         & ! Land/sea mask
 ,blowst(npnts)          ! mask for those points at which stability is low
                         ! enough for convection to occur
REAL, INTENT(IN) ::  &
  timestep             ! convection timestep (s)
                       !   (= model timestep if conv. called once)

REAL,INTENT(IN) :: &
  thek(npnts)      & ! potential temperature of cloud environment in layer k (K)
 ,thekp1(npnts)    & ! potential temperature of cloud environment in layer k+1
 ,qek(npnts)       & ! mixing ratio of cloud environment in layer k (kg/kg)
 ,qekp1(npnts)     & ! mixing ratio of cloud environment in layer k+1 (kg/kg)
 ,qclek(npnts)     & ! liquid condensate mixing ratio of cloud environment 
                     ! layer k (kg/kg)
 ,qclekp1(npnts)   & ! liquid condensate mixing ratio of cloud environment 
                     ! layer k+1 (kg/kg)
 ,qcfek(npnts)     & ! frozen condensate mixing ratio of cloud environment 
                     ! layer k (kg/kg)
 ,qcfekp1(npnts)   & ! frozen condensate mixing ratio of cloud environment 
                     ! layer k+1 (kg/kg)
 ,cflek(npnts)     & ! Liquid cloud volume fraction of cloud environment 
                     ! layer k 
 ,cflekp1(npnts)   & ! Liquid cloud volume fraction of cloud environment 
                     ! layer k+1
 ,cffek(npnts)     & ! Frozen cloud volume fraction of cloud environment 
                     ! layer k 
 ,cffekp1(npnts)   & ! Frozen cloud volume fraction of cloud environment 
                     ! layer k+1
 ,bcfek(npnts)     & ! Total cloud volume fraction of cloud environment 
                     ! layer k 
 ,bcfekp1(npnts)   & ! Total cloud volume fraction of cloud environment 
                     ! layer k+1
 ,uek(npnts)       & ! U in environment in layer k (m/s)
 ,uekp1(npnts)     & ! U in environment in layer k+1 (m/s)
 ,vek(npnts)       & ! V in environment in layer k (m/s)
 ,vekp1(npnts)       ! V in environment in layer k+1 (m/s)


REAL,INTENT(IN) ::      &
  traek(np_full,ntra)   & ! tracer content of cloud environment 
                          ! layer k (kg/kg)
 ,traekp1(np_full,ntra) & ! tracer content of cloud environment 
                          ! layer k+1 (kg/kg)
 ,trapkp1(np_full,ntra)   ! parcel tracer content in layer k+1 (kg/kg)

REAL,INTENT(IN) ::      &
  qsekp1(npnts)         & ! saturation mixing ratio of cloud environment
                          ! in layer k+1 (kg/kg)
 ,dqskp1(npnts)         & ! gradient of saturation mixing ratio with 
                          ! potential temperature for the cloud environment
                          ! in layer k+1 (kg/kg)
 ,pstar(npnts)          & ! Surface pressure (Pa)

 ,upkp1(npnts)          & ! parcel U in layer k+1 (m/s)
 ,vpkp1(npnts)          & ! parcel V in layer k+1 (m/s)
 ,rbuoy(npnts)          & ! parcel buoyancy in layer k+1 (K)
 ,qsek(npnts)           & ! saturation mixing ratio of cloud environment 
                          ! in layer k (kg/kg)
 ,dqsk(npnts)           & ! gradient of saturation mixing ratio with
                          ! potential temperature for the cloud environment
                          ! of layer k (kg/kg/K)
 ,thpi(npnts)           & ! Initial parcel potential temperature (K)
 ,qpi(npnts)            & ! Initial parcel mixing ratio (kg/kg)
 ,expi(npnts)             ! Initial parcel Exner pressure


REAL, INTENT(IN) ::  &
  ekp14(npnts)       & ! Entrainment coefficient at level k+1/4 multiplied by 
                       ! appropriate layer thickness
 ,ekp34(npnts)       & ! Entrainment coefficient at level k+3/4 multiplied by 
                       ! appropriate layer thickness
 ,amdetk(npnts)        ! Mixing detrainment coefficient at level k multiplied 
                       ! by appropriate layer thickness

REAL, INTENT(IN) ::  &
  pk(npnts)          & ! pressure at mid-point of layer k (Pa)
 ,pkp1(npnts)        & ! pressure at mid-point of layer k+1   (Pa)
 ,exk(npnts)         & ! Exner ratio at mid-point of layer k
 ,exkp1(npnts)       & ! Exner ratio at mid-point of layer k+1
 ,delexkp1(npnts)    & ! Difference in Exner ratio between mid-points of
                       !  layers k and k+1
 ,delpk(npnts)       & ! Difference in pressure across layer k (Pa)
 ,delpkp1(npnts)     & ! Difference in pressure across layer k+1 (Pa)
 ,delp_uv_k(npnts)   & ! pressure difference across UV layer k (Pa)
 ,delp_uv_kp1(npnts)   ! pressure difference across UV layer K+1 (Pa)

REAL, INTENT(IN) ::   &
  t1_sd(npnts)        & ! Standard deviation of turbulent fluctuations of
                        !  layer 1 temperature (K).
 ,q1_sd(npnts)        & ! Standard deviation of turbulent fluctuations of 
                        !  layer 1  humidity (kg/kg).
 ,thpixs_v(npnts)     & !  parcel excesses of theta
 ,qpixs_v(npnts)      & !  parcel excesses of q
 ,xsbmin_v(npnts)     & !  Threshold buoyancy for forced detrainment
 ,rbuoy_p_old(npnts)    ! buoyancy from previous level

!----------------------------------------------------------------------
! Variables which are input but which are also updated in this routine
!----------------------------------------------------------------------

LOGICAL,INTENT(INOUT) :: & 
  binit(npnts)           &!IN mask for those points at which convection is 
                          !   occuring
                          !OUT mask reset to false if convection terminates in
                          !    layer.
 ,bgmk(npnts)            &!IN Mask for parcels which are saturated in layer k
                          !OUT MASK for parcels which are saturated in layer k+1
 ,bgmkp1(npnts)           !IN mask for parcels which are saturated in layer k+1
                          !OUT - value not wanted, but modified by parcel.

REAL, INTENT(INOUT) :: &  
  thpk(npnts)          & ! IN parcel potential temperature in layer k (K)
                         ! OUT parcel potential temperature in layer k+1 (K)

 ,thpkp1(npnts)        & ! IN parcel potential temperature in layer k+1 (K)
                         ! OUT - value not wanted, but modifed by parcel

 ,qpk(npnts)           & ! IN  parcel mixing ratio in layer k (kg/kg)
                         ! OUT parcel mixing ratio in layer k+1 (kg/kg)

 ,qpkp1(npnts)         & ! IN parcel mixing ratio in layer k+1 (kg/kg)
                         ! OUT - value not wanted, but modified by parcel.

 ,qclpk(npnts)         & ! IN  parcel liquid condensate mixing ratio in 
                         !     layer k (kg/kg)
                         ! OUT parcel liquid condensate mixing ratio in 
                         !     layer k+1 (kg/kg)
 ,qcfpk(npnts)         & ! IN  parcel frozen condensate mixing ratio in 
                         !     layer k (kg/kg)
                         ! OUT parcel frozen condensate mixing ratio in 
                         !     layer k+1 (kg/kg)
 ,qclpkp1(npnts)       & ! IN  parcel liquid condensate mixing ratio in layer
                         !     k+1 after dry ascent only (kg/kg)
                         ! OUT parcel liquid condensate in layer k+1 (kg/kg)

 ,qcfpkp1(npnts)       & ! IN  parcel frozen condensate mixing ratio in layer
                         !     k+1 after dry ascent only (kg/kg)
                         ! OUT parcel frozen condensate in layer k+1 (kg/kg)

 ,xsqkp1(npnts)        & ! IN Excess water in parcel after lifting
                         ! layer k to K+1 (kg/kg)
                         ! OUT - value not wanted, but modified by parcel.
 ,upk(npnts)           & ! IN  parcel U in layer k (m/s)
                         ! OUT parcel U in layer k+1 (m/s)

 ,vpk(npnts)           & ! IN  parcel V in layer k (m/s)
                         ! OUT parcel V in layer k+1 (m/s)

 ,trapk(np_full,ntra)  & ! IN  parcel tracer content in layer k (kg/kg)
                         ! OUT parcel tracer content in layer k+1  (kg/kg)

 ,xpk(npnts)             ! IN parcel cloud water in layer k (kg/kg)
                         ! OUT parcel cloud water in layer k+1 (kg/kg)

! NOTE: XPK is identical with condensate mixing ratio of cloud parcel
!       XPK will be redundant when L_calc_dxek == .True.  because
!       it can be calculated directly from (Qclpk + Qcfpk).

REAL, INTENT(INOUT) :: &  
  flxk(npnts)          & ! IN parcel massflux in layer k (Pa/s) 
                         ! OUT parcel massflux in layer k+1 (Pa/s)
 ,dthek(npnts)         & ! IN Increment to model potential temperature in layer
                         !    k due to convection (may be non-zero due to a 
                         !   previous split final detrainment calculation) (K/s)
                         ! OUT updated Increment to model potential temperature 
                         !     in layer k due to convection (K/s)

 ,dqek(npnts)          & ! IN Increment to model mixing ratio in layer k due 
                         !    to convection (may be non-zero due to a previous
                         !    split final detrainment calculation) (kg/kg/s)
                         ! OUT updated Increment to model mixing ratio
                         !     in layer k due to convection (kg/kg/s)

 ,dqclek(npnts)        & ! IN Increment to model liquid condensate mixing ratio
                         !    in layer k due to convection (may be non-zero due
                         !    to a previous split final detrainment 
                         !     calculation) (kg/kg/s)
                         ! OUT updated Increment to model liquid condensate
                         !     in layer k due to convection (kg/kg/s)

 ,dqcfek(npnts)        & ! IN Increment to model frozen condensate mixing ratio
                         !    in layer k due to convection (may be non-zero due
                         !    to a previous split final detrainment 
                         !     calculation) (kg/kg/s)
                         ! OUT updated Increment to model frozen condensate
                         !     in layer k due to convection (kg/kg/s)

 ,dbcfek(npnts)        & ! IN Increment to model total cloud volume fraction in
                         !    layer k due to convection (may be non-zero due
                         !    to a previous split final detrainment 
                         !     calculation) (/s)
                         ! OUT updated Increment to model total cloud volume 
                         !     in layer k due to convection (/s)

 ,dcflek(npnts)        & ! IN Increment to model liquid cloud volume fraction in
                         !    layer k due to convection (may be non-zero due
                         !    to a previous split final detrainment 
                         !     calculation) (/s)
                         ! OUT updated Increment to model liquid cloud volume 
                         !     in layer k due to convection (/s)

 ,dcffek(npnts)        & ! IN Increment to model frozen cloud volume fraction in
                         !    layer k due to convection (may be non-zero due
                         !    to a previous split final detrainment 
                         !     calculation) (/s)
                         ! OUT updated Increment to model frozen cloud volume 
                         !     in layer k due to convection (/s)

 ,duek(npnts)          & ! IN Increment to model U in layer k due to convection 
                         !  (m/s**2)
                         ! OUT updated Increment to model U in layer k
                         !     due to convection (m/s**2)

 ,dvek(npnts)            ! IN Increment to model V in layer k due to convection 
                         !  (m/s**2)
                         ! OUT updated Increment to model V in layer k
                         !     due to convection (m/s**2)

REAL, INTENT(INOUT) :: &  
  dtraek(np_full,ntra)   ! IN Increment to model tracer in layer k due to
                         !    convection(may be non-zero due to a previous
                         !    final detrainment  calculation) (kg/kg/s) 
                         !     calculation) (/s)
                         ! OUT updated Increment to model tracer in layer k
                         !     due to convection (kg/kg/s)

REAL, INTENT(INOUT) :: &  
  tcw(npnts)           & ! IN total condensed water summed to 
                         !       layer k(kg/m**2/s)
                         ! OUT updated total condensed water summed to 
                         !       layer k+1 (kg/m**2/s)
 ,depth(npnts)         & ! IN depth of convective cloud to layer k (m)
                         ! OUT updated depth of convective cloud to layer k+1(m)
 ,cclwp(npnts)         & ! IN Condensed water path summed to layer k (kg/m**2)
                         ! OUT updated condensed water path summed to 
                         !     layer k+1 (kg/m**2)
 ,cape(npnts)          & ! IN Convective available potential energy up to the
                         !    current convecting layer   
                         ! OUT Convective available potential energy including
                         !     addition due to the CAPE within the current layer
 ,dcpbydt(npnts)       & ! IN  Rate of change of CAPE
                         ! OUT Rate of change of CAPE including contribution
                         !     from current layer
 ,eflux_u_ud(npnts)    & ! IN Eddy flux of momentum to UD at bottom of layer
                         ! OUT Eddy flux of momentum to UD at top of layer
 ,eflux_v_ud(npnts)      !  IN Eddy flux of momentum to UD at bottom of layer
                         ! OUT Eddy flux of momentum to UD at top of layer 



!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------
INTEGER, INTENT(OUT) :: & 
  iccb(npnts)           & ! convective cloud base_level
 ,icct(npnts)           & ! convective cloud top level
 ,lcbase(npnts)         & ! Lowest conv. cloud base level
 ,lctop(npnts)            ! Lowest conv. cloud top level


LOGICAL, INTENT(OUT) :: &  
  bterm(npnts)            ! Mask for parcels which terminate in layer k+1


REAL, INTENT(OUT) :: &  
  prekp1(npnts)      & ! precipitation from parcel as it rises from layer 
                       !  k to k+1 (kg/m**2/s)
 ,dthekp1(npnts)     & ! Increment to model potential temperature in layer 
                       ! k+1 due to convection (K/s)
 ,dqekp1(npnts)      & ! Increment to model mixing ratio in layer k+1 
                       ! due to convection (kg/kg/s)
 ,dqclekp1(npnts)    & ! Increment to model liquid condensate mixing ratio 
                       ! in layer k+1 due to convection (kg/kg/s)
 ,dqcfekp1(npnts)    & ! Increment to model frozen condensate mixing ratio 
                       ! in layer k+1 due to convection (kg/kg/s)
 ,dbcfekp1(npnts)    & ! Increment to model total cloud volume fraction
                       ! in layer k+1 due to convection (/s)
 ,dcflekp1(npnts)    & ! Increment to model liquid cloud volume fraction
                       ! in layer k+1 due to convection (/s)
 ,dcffekp1(npnts)    & ! Increment to model frozen cloud volume fraction
                       ! in layer k+1 due to convection (/s)
 ,duekp1(npnts)      & ! Increment to model U in layer k+1 due to 
                       ! convection (m/s**2)
 ,dvekp1(npnts)        ! Increment to model V in layer k+1 due to 
                       ! convection (m/s**2)


REAL, INTENT(OUT) ::      &  
  dtraekp1(np_full,ntra)    ! Increment to model tracer in layer k+1 due to 
                            ! convection (kg/kg/s)

REAL, INTENT(OUT) ::   &  
  cca(npnts)           & ! convective cloud amount (%)
 ,ccw(npnts)           & ! convective cloud water(g/kg) on model levels
 ,lcca(npnts)          & ! Lowest conv. cloud amount (%)
 ,deltak(npnts)        & ! parcel forced detrainment rate layer k multiplied
                         ! by appropriate layer thickness
 ,flxkp12(npnts)       & ! half level mass flux
 ,max_cfl_c(npnts)     & ! CFL ratio
 ,relh(npnts)          & ! Relative humidity integral (average when convection
                         ! terminates)
 ,dptot(npnts)           ! Delta P integral


REAL, INTENT(OUT) ::  &
  rbuoy_p_here(npnts) & ! buoyancy for SCM diags
 ,the_here(npnts)     & ! TH environ for SCM diags
 ,thp_here(npnts)     & ! TH parcel for SCM diags
 ,qe_here(npnts)      & ! Q environ for SCM diags
 ,qp_here(npnts)        ! Q parcel for SCM diags

!-------------------------------------------------------------------------------
! Local variables

INTEGER ::        & 
  i,ktra           ! loop counters

REAL ::            & 
  thrk(npnts)      & ! parcel detrainment potential temperature in layer k (K)
 ,qrk(npnts)       & ! parcel detrainment mxing ratio in layer k (kg/kg)
 ,xpkp1(npnts)     & ! parcel cloud water in layer k+1 (kg/kg)
 ,flxkp1(npnts)      ! parcel massflux in layer k+1 (Pa/s)

REAL ::       &
 thvp         & ! Virtual tempature of parcel
,thve         & ! Virtual tempature of environment 
,rho          & ! Density required in CAPE calculations
,tmp_dcpbydt    ! Temporary dcpbydt

REAL ::                & 
  dqek_nonpc2(npnts)   & ! increment to qek if PC2 was not in place
 ,dthek_nonpc2(npnts)  & ! increment to thek if PC2 was not in place
 ,dqekp1_nonpc2(npnts) & ! increment to thek if PC2 was not in place
 ,dthekp1_nonpc2(npnts)  ! increment to thekp1 if PC2 was not in place

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('CONVEC2_4A5A',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
!  Complete lifting parcel to layer k+1.
!
!  SUBROUTINE parcel
!
!  UM Documentation paper 27
!  Sections (5),(6),(7),(8),(9)
! ----------------------------------------------------------------------
!
! DEPENDS ON: parcel_4a5a
 CALL parcel_4a5a (k,npnts,nlev,ad_on, new_termc, start_lev,             &
                   l_q_interact, l_shallow,                              &
                   bgmk,bwkp1,blowst, bland,                             &
                   pstar,thekp1,thek,qekp1,qek,qsek,qsekp1,dqsk,dqskp1,  &
                   thpk,qpk,xpk,qclpk,qcfpk,                             &
                   thpi,qpi,expi,rbuoy,rbuoy_p_old,xsbmin_v,             &
                   ekp14,ekp34,amdetk,                                   &
                   delpkp1,pk,pkp1,exk,exkp1,delexkp1,                   &
                   ! In/out     
                   bgmkp1,                                               &
                   thpkp1,qpkp1,qclpkp1,qcfpkp1,flxk,xsqkp1,             &
                   tcw,depth,cclwp,                                      &
                   ! Out            
                   iccb,icct,lcbase,lctop,                               &
                   bterm,                                                &
                   prekp1,thrk,qrk,xpkp1,flxkp1,                         & 
                   cca,ccw,lcca,deltak, flxkp12,                         &
                   rbuoy_p_here,the_here,thp_here,qe_here,qp_here)
!
! ----------------------------------------------------------------------
!  Calculate the effect on the environment (except for the 
!  evaporation of precipitation and changeof phase)
!
!  SUBROUTINE environ
!
!  UM Documenataion paper 27
!  Section (10)
! ----------------------------------------------------------------------
!
! DEPENDS ON: environ
CALL environ (k,npnts,np_full,ntra,sdet_on,                         &
              l_tracer,l_mom_gk,l_calc_dxek,l_q_interact,           &
              bwk, bwkp1, bterm, blowst,                            &
              l_shallow, cumulus, l_mid,                            &
              timestep,                                             &
              thek, qek, qclek, qcfek, bcfek, cflek, cffek,         &
              thekp1,qekp1,qclekp1,qcfekp1,bcfekp1,cflekp1,cffekp1, &
              uek, vek, uekp1, vekp1, traek, traekp1,               &
              thpk, qpk, qclpk, qcfpk, flxk, thrk, qrk,             &
              thpkp1, qpkp1, qclpkp1, qcfpkp1, flxkp1,              &
              upk, vpk, upkp1, vpkp1, trapk, trapkp1,               &
              deltak, ekp14, exk, exkp1, delpk, delpkp1,            &
              delp_uv_k, delp_uv_kp1,                               &
              thpixs_v, qpixs_v, amdetk, t1_sd, q1_sd,              &
              dthek, dqek, dqclek, dqcfek, dbcfek, dcflek, dcffek,  &
              xpk, xpkp1,                                           & 
              duek, dvek, eflux_u_ud, eflux_v_ud, dtraek,           &
              dthekp1, dqekp1, dqclekp1, dqcfekp1,                  &
              dbcfekp1, dcflekp1, dcffekp1,                         &
              dqek_nonpc2, dthek_nonpc2,                            &
              dqekp1_nonpc2, dthekp1_nonpc2,                        &
              duekp1, dvekp1, dtraekp1,                             &
              max_cfl_c )

      
!-----------------------------------------------------------------------
! Reset BINIT where convection has terminated.
!-----------------------------------------------------------------------
DO i=1,npnts
  binit(i) = .NOT.bterm(i)
END DO

! ---------------------------------------------------------------------
!  Calculate contribution to CAPE and rate of change of CAPE  due to
!  the updraught
! ---------------------------------------------------------------------


DO i=1,npnts
  IF (cnv_wat_load_opt == 1) THEN
    thvp=thpk(i)*(1.0+c_virtual*qpk(i)-qclpk(i)-qcfpk(i))
    thve=thek(i)*(1.0+c_virtual*qek(i)-qclek(i)-qcfek(i))
    rho=pk(i)/(r*thve*exk(i))
  ELSE
    thvp=thpk(i)*(1.0+c_virtual*qpk(i))
    thve=thek(i)*(1.0+c_virtual*qek(i))
    rho=pk(i)/(r*thek(i)*exk(i))
  END IF

  cape(i)=cape(i)+(thvp-thve)*delpk(i)/(rho*thve)
  relh(i)=relh(i)+(qek(i)/qsek(i))*delpk(i)
  dptot(i)=dptot(i)+delpk(i)

  IF (cnv_wat_load_opt == 1) THEN
!   Include the effects of water loading
    tmp_dcpbydt=(dthek(i)*(1.0+c_virtual*qek(i)-                &
             qclek(i)-qcfek(i))+thek(i)*                        &
             (c_virtual*dqek(i) - dqclek(i) - dqcfek(i)))*      &
             (delpk(i)/(rho*thve))
  ELSE
!   Water loading is not included and therefore use the original
!   calculation
    tmp_dcpbydt=(dthek_nonpc2(i)*(1.0+c_virtual*qek(i))+        &
             c_virtual*thek(i)*dqek_nonpc2(i))*                 &
             (delpk(i)/(rho*thve))
  END IF

  IF (tmp_dcpbydt  >   0.0) THEN
    dcpbydt(i) = dcpbydt(i) + tmp_dcpbydt
  END IF

  IF (bterm(i)) THEN
    IF (cnv_wat_load_opt == 1) THEN
      thvp=thpkp1(i)*(1.0+c_virtual*qpkp1(i)-qclpkp1(i)-qcfpkp1(i))
      thve=thekp1(i)*(1.0+c_virtual*qekp1(i)-qclekp1(i)-qcfekp1(i))
      rho =pkp1(i)/(r*thve*exkp1(i))
    ELSE
      thvp=thpkp1(i)*(1.0+c_virtual*qpkp1(i))
      thve=thekp1(i)*(1.0+c_virtual*qekp1(i))
      rho =pkp1(i)/(r*thekp1(i)*exkp1(i))
    END IF

    cape(i)=cape(i)+(thvp-thve)*delpkp1(i)/(rho*thve)
    relh(i)=relh(i)+(qekp1(i)/qsekp1(i))*delpkp1(i)
    dptot(i)=dptot(i)+delpkp1(i)

    IF (cnv_wat_load_opt == 1) THEN
!     Include the effects of water loading
      tmp_dcpbydt=(dthekp1(i)*(1.0+c_virtual*qekp1(i)-            &
               qclekp1(i)-qcfekp1(i))+thekp1(i)*                  &
               (c_virtual*dqekp1(i) - dqclekp1(i) - dqcfekp1(i)))*&
               (delpkp1(i)/(rho*thve))
    ELSE
!     Water loading is not included and therefore use the original
!     calculation
      tmp_dcpbydt=(dthekp1_nonpc2(i)*(1.0+c_virtual*qekp1(i))+    &
               c_virtual*thekp1(i)*dqekp1_nonpc2(i))*             &
               (delpkp1(i)/(rho*thve))
    END IF

    IF (tmp_dcpbydt  >   0.0) THEN
      dcpbydt(i) = dcpbydt(i) + tmp_dcpbydt
    END IF
  END IF

END DO


! ---------------------------------------------------------------------
!  Swap parcel values ready for the next part of ascent
!  from layer k+1 to k+2 (looping over number of tracers).
! ---------------------------------------------------------------------

DO i=1,npnts
  thpk(i) = thpkp1(i)
  qpk(i) = qpkp1(i)
  xpk(i) = xpkp1(i)
  qclpk(i) = qclpkp1(i)
  qcfpk(i) = qcfpkp1(i)
  flxk(i) = flxkp1(i)
  bgmk(i) = bgmkp1(i)
END DO

IF(l_mom_gk)THEN
  DO i=1,npnts
    upk(i) = upkp1(i)
    vpk(i) = vpkp1(i)
  END DO
END IF

IF(l_tracer)THEN

  DO ktra = 1,ntra
    DO i=1,npnts
      IF(binit(i))THEN
        trapk(i,ktra) = trapkp1(i,ktra)
      END IF
    END DO
  END DO

END IF

IF (lhook) CALL dr_hook('CONVEC2_4A5A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE convec2_4a5a
