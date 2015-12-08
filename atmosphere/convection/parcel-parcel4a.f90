! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Completes lifting of the parcel from layer k to k+1
!
SUBROUTINE parcel_4a5a (k,npnts,nlev,ad_on, new_termc, start_lev,             &
                        l_q_interact, l_shallow,                              &
                        bgmk,bwkp1,blowst, bland,                             &
                        pstar,thekp1,thek,qekp1,qek,qsek,qsekp1,dqsk,dqskp1,  &
                        thpk,qpk,xpk,qclpk,qcfpk,                             &
                        thpi,qpi,expi,rbuoy,rbuoy_p_old,xsbmin,               &
                        ekp14,ekp34,amdetk,                                   &
                        delpkp1,pk,pkp1,exk,exkp1,delexkp1,                   &
                        ! In/out     
                        bgmkp1,                                               &
                        thpkp1,qpkp1,qclpkp1,qcfpkp1,flxk, xsqkp1,            &
                        tcw,depth,cclwp,                                      &
                        ! Out            
                        iccb,icct,lcbase,lctop,                               &
                        bterm,                                                &
                        prekp1,thrk,qrk,xpkp1,flxkp1,                         & 
                        cca,ccw,lcca, deltak, flxkp12,                        &
                        rbuoy_p_here,the_here,thp_here,qe_here,qp_here)


USE cv_run_mod, ONLY:                                                    &
    r_det

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
! Completes lifting of the parcel from layer k to k+1
! Calls detrain, term_con and cloud_w
! An intial mass flux is Calculated.
! Detrain  - carries out the forced detrainment.
! term_con - tests for any convection which is terminating in layer k+1
! cloud_w  - carries out the cloud microphysics calculation.
!
!  See UM Documentation paper No 27
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
  k                    & ! present model layer
 ,npnts                & ! Number of points
 ,nlev                 & ! Number of model levels for calculations
 ,ad_on                & ! Flag for adaptive detrainment 
 ,new_termc              ! Flag for simplified termination of convection

INTEGER, INTENT(IN) :: &
  start_lev(npnts)       ! Level at which convection is initiated

LOGICAL, INTENT(IN) :: &
  l_q_interact           ! Switch allows overwriting of parcel
                         ! variables when calculating condensate
                         ! increments (will alter results).

LOGICAL, INTENT(IN) ::  &
  l_shallow(npnts)      & ! true in shallow convection
 ,bgmk(npnts)           & ! mask for parcels whcih are saturated in layer k
 ,bwkp1(npnts)          & ! mask for whether condensate is liquid in layer k+1
 ,blowst(npnts)         & ! mask for those points at which stability is low
                          ! enough for convection to occur
 ,bland(npnts)            ! Land/sea mask

REAL,INTENT(IN) :: &
  pstar(npnts)     & ! Surface pressure (Pa)
 ,thekp1(npnts)    & ! potential temperature of cloud environment in layer k+1
 ,thek(npnts)      & ! potential temperature of cloud environment in layer k (K)
 ,qekp1(npnts)     & ! mixing ratio of cloud environment in layer k+1 (kg/kg)
 ,qek(npnts)       & ! mixing ratio of cloud environment in layer k (kg/kg)
 ,qsek(npnts)      & ! saturation mixing ratio of cloud environment 
                     ! in layer k (kg/kg)
 ,qsekp1(npnts)    & ! saturation mixing ratio of cloud environment
                     ! in layer k+1 (kg/kg)
 ,dqsk(npnts)      & ! gradient of saturation mixing ratio with potential 
                     ! temperature for the cloud environment of layer k
                     ! (kg/kg/K)
 ,dqskp1(npnts)      ! gradient of saturation mixing ratio with 
                     ! potential temperature for the cloud environment
                     ! in layer k+1 (kg/kg/K)

REAL, INTENT(IN) :: &  
  thpk(npnts)       & ! parcel potential temperature in layer k (K)
 ,qpk(npnts)        & ! parcel mixing ratio in layer k (kg/kg)
 ,qclpk(npnts)      & ! parcel liquid condensate mixing ratio in 
                      !     layer k (kg/kg)
 ,qcfpk(npnts)      & ! parcel frozen condensate mixing ratio in 
                      !     layer k (kg/kg)
 ,thpi(npnts)       & ! Initial parcel potential temperature (K)
 ,qpi(npnts)        & ! Initial parcel mixing ratio (kg/kg)
 ,expi(npnts)       & ! Initial parcel Exner pressure
 ,rbuoy(npnts)      & ! parcel buoyancy in layer k+1 (K)
 ,rbuoy_p_old(npnts)& ! buoyancy from previous level
 ,xsbmin(npnts)       !  Threshold buoyancy for forced detrainment
                      !    Function of delta P

REAL, INTENT(IN) ::  &
  ekp14(npnts)       & ! Entrainment coefficient at level k+1/4 multiplied by 
                       ! appropriate layer thickness
 ,ekp34(npnts)       & ! Entrainment coefficient at level k+3/4 multiplied by 
                       ! appropriate layer thickness
 ,amdetk(npnts)        ! Mixing detrainment coefficient at level k multiplied 
                       ! by appropriate layer thickness

REAL, INTENT(IN) ::  &
  delpkp1(npnts)     & ! Difference in pressure across layer k+1 (Pa)
 ,pk(npnts)          & ! pressure at mid-point of layer k (Pa)
 ,pkp1(npnts)        & ! pressure at mid-point of layer k+1   (Pa)
 ,exk(npnts)         & ! Exner ratio at mid-point of layer k
 ,exkp1(npnts)       & ! Exner ratio at mid-point of layer k+1
 ,delexkp1(npnts)      ! Difference in Exner ratio between mid-points of
                       !  layers k and k+1

!----------------------------------------------------------------------
! Variables which are input but which are also updated in this routine
!----------------------------------------------------------------------

LOGICAL,INTENT(INOUT) :: & 
  bgmkp1(npnts)           !IN Mask for parcels which are saturated in layer k+1
                          !OUT MASK for parcels which are saturated in layer k+1

REAL, INTENT(INOUT) :: &  
  thpkp1(npnts)        & ! IN parcel potential temperature in layer k+1  after
                         ! entrainment and latent heating (K)
                         ! OUT parcel potential temperature in layer k+1 (K)
                         ! after forced detrainment
 ,qpkp1(npnts)         & ! IN  parcel mixing ratio in layer k+1  after
                         ! entrainment and latent heating (kg/kg)
                         ! OUT parcel mixing ratio in layer k+1 (kg/kg)
                         ! after forced detrainment
 ,qclpkp1(npnts)       & ! IN  parcel liquid condensate mixing ratio in layer
                         !     k+1 after dry ascent only (kg/kg)
                         ! OUT parcel liquid condensate in layer k+1 (kg/kg)
 ,qcfpkp1(npnts)       & ! IN  parcel frozen condensate mixing ratio in layer
                         !     k+1 after dry ascent only (kg/kg)
                         ! OUT parcel frozen condensate in layer k+1 (kg/kg)
 ,xpk(npnts)           & ! IN parcel cloud water in layer k (kg/kg)
                         ! OUT overwritten with qcl+qcf in cloud_w for layer 
                         !     k if PC2
 ,flxk(npnts)          & ! IN parcel massflux in layer k (Pa/s) 
                         ! OUT parcel massflux in layer k (set if convection
                         !     initiated from layer k) (Pa/s)
 ,xsqkp1(npnts)          ! Excess parcel water after lifting from layer k to k+1
                         !   (kg/kg)

REAL, INTENT(INOUT) :: &  
  tcw(npnts)           & ! IN total condensed water summed to 
                         !       layer k(kg/m**2/s)
                         ! OUT updated total condensed water summed to 
                         !       layer k+1 (kg/m**2/s)
 ,depth(npnts)         & ! IN depth of convective cloud to layer k (m)
                         ! OUT updated depth of convective cloud to layer k+1(m)
 ,cclwp(npnts)           ! IN Condensed water path summed to layer k (kg/m**2)
                         ! OUT updated condensed water path summed to 
                         !     layer k+1 (kg/m**2)


!---------------------------------------------------------------------
! Variables which are output
!---------------------------------------------------------------------
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
 ,thrk(npnts)        & ! parcel detrainment potential temperature in layer k (K)
 ,qrk(npnts)         & ! parcel detrainment mxing ratio in layer k (kg/kg)
 ,xpkp1(npnts)       & ! parcel cloud water in layer k+1 (kg/kg)
 ,flxkp1(npnts)        ! parcel massflux in layer k+1 (Pa/s)

REAL, INTENT(OUT) ::   &  
  cca(npnts)           & ! convective cloud amount (%)
 ,ccw(npnts)           & ! convective cloud water(g/kg) on model levels
 ,lcca(npnts)          & ! Lowest conv. cloud amount (%)
 ,deltak(npnts)        & ! parcel forced detrainment rate layer k multiplied
                         ! by appropriate layer thickness
 ,flxkp12(npnts)         ! half level mass flux

REAL, INTENT(OUT) ::  &
  rbuoy_p_here(npnts) & ! buoyancy for SCM diags
 ,the_here(npnts)     & ! TH environ for SCM diags
 ,thp_here(npnts)     & ! TH parcel for SCM diags
 ,qe_here(npnts)      & ! Q environ for SCM diags
 ,qp_here(npnts)        ! Q parcel for SCM diags


!-------------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

INTEGER ::        & 
  i               & ! loop counter 
 ,ndet              ! Compress vector length for the detrainment calculation

INTEGER ::        &
  index1(npnts)     ! Index for compress and expand

LOGICAL ::        &
  bwkp1_c(npnts)  & ! Compressed mask for whether condensate is liquid 
                    ! in layer k+1
 ,bgmk_c(npnts)   & ! Compressed mask for parcels which are saturated
                    ! in layer k
 ,bgmkp1_c(npnts) & ! Compressed mask for parcels which are saturated
                    ! in layer k+1
 ,bdetk(npnts)      ! Mask for points under going forced detrainment


REAL ::            &
  xsbmin_ad(npnts) & ! xsbmin adpative (NOTE, will be different at 
                     ! different points)
 ,thek_c(npnts)    & ! Compressed potential temperature of
                     ! cloud environment in layer k (K)
 ,thekp1_c(npnts)  & ! Compressed potential temperature of
                     ! cloud environment in layer k+1 (K)
 ,qek_c(npnts)     & ! Compressed mixing ratio of cloud
                     ! environment in layer k (kg/kg)
 ,qekp1_c(npnts)   & ! Compressed mixing ratio of cloud
                     ! environment in layer k+1 (kg/kg)
 ,qsek_c(npnts)    & ! Compressed saturation mixing ratio of
                     ! cloud environment in layer k (kg/kg)
 ,dqsk_c(npnts)    & ! Compressed gradient of saturation mixing ratio  with 
                     ! potential temperature for the cloud environment
                     ! of layer k (kg/kg/K)
 ,qsekp1_c(npnts)  & ! Compressed saturation mixing ratio of
                     ! cloud environment in layer k+1 (kg/kg)
 ,dqskp1_c(npnts)    ! Compressed gradient of saturation mixing ratio  with 
                     ! potential temperature for the cloud environment 
                     ! of layer k+1 (kg/kg/K)

REAL ::            &
  thpk_c(npnts)    & ! Compressed parcel potential temperature in layer k (K)
 ,qpk_c(npnts)     & ! Compressed parcel mixing ratio  in layer k  (kg/kg)
 ,thpkp1_c(npnts)  & ! Compressed parcel potential temperature in layer k+1 (K)
 ,qpkp1_c(npnts)   & ! Compressed parcel mixing ratio  in layer k  (kg/kg)
 ,xsqkp1_c(npnts)    ! Excess parcel water after lifting from layer k to k+1
                     ! (kg/kg)
REAL ::            &
  thrk_c(npnts)    & ! Compressed parcel detrainment potential temperature
                     !  in layer k (K)
 ,qrk_c(npnts)     & ! Compressed parcel detrainment mixing ratio
                     ! in layer k (kg/kg)
 ,deltak_c(npnts)  & ! Compressed parcel forced detrainment coefficient in
                     ! layer k multiplied by appropriate layer thickness
 ,ekp14_c(npnts)   & ! Compressed in entrainment coefficient at level k+1/4
                     ! multiplied by appropriate layer thickness
 ,ekp34_c(npnts)     ! Compressed in entrainment coefficient at level k+3/4
                     ! multiplied by appropriate layer thickness

REAL ::            &
  pk_c(npnts)      & ! Compressed pressure at level k (PA)
 ,pkp1_c(npnts)    & ! Compressed pressure at level k+1 (PA)
 ,xsbmin_c(npnts)  & ! Compressed threshold for forced detrainment
                     !  Function of delta P
 ,exk_c(npnts)     & ! Compressed exner function at level k
 ,exkp1_c(npnts)     ! Compressed exner function at level k+1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('PARCEL_4A5A',zhook_in,zhook_handle)
!---------------------------------------------------------------------

DO i=1,npnts

! ---------------------------------------------------------------------
!  Calculate mask for those points under going forced detrainment
!
!  UM Documentation paper 27
!  Section (6), equation (23)
! ---------------------------------------------------------------------

  bdetk(i) = rbuoy(i)  <   xsbmin(i)

  !  Calculate XSBMIN_AD
  !  Use adaptive detrainment if ADAPTIVE ON flag is set to 1
  !  making sure that adaptive min. buoyancy doesn't fall below 0.0
  IF (ad_on  /=  1) THEN
    xsbmin_ad(i) = xsbmin(i)
  ELSE
    xsbmin_ad(i) = r_det * rbuoy_p_old(i) + (1-r_det)*rbuoy(i)
    xsbmin_ad(i) = MAX(xsbmin_ad(i),0.0)
  END IF

  !  Set mask of points for forced detrainment
  bdetk(i) = rbuoy(i)  <   xsbmin_ad(i)

  rbuoy_p_here(i)=rbuoy(i)

  the_here(i)=thek(i)
  qe_here(i)=qek(i)
  thp_here(i)=thpk(i)
  qp_here(i)=qpk(i)

END DO

! ----------------------------------------------------------------------
!  Compress all input arrays for the forced detrainment calculations
! ----------------------------------------------------------------------

ndet = 0
DO i=1,npnts
  IF (bdetk(i)) THEN
    ndet = ndet + 1
    index1(ndet) = i
  END IF
END DO

IF (ndet  /=  0) THEN
  DO i=1,ndet
    thek_c(i)  = thek(index1(i))
    qek_c(i)   = qek(index1(i))
    thpk_c(i)  = thpk(index1(i))
    qpk_c(i)   = qpk(index1(i))
    qsek_c(i)  = qsek(index1(i))
    dqsk_c(i)  = dqsk(index1(i))
    thekp1_c(i)= thekp1(index1(i))
    qekp1_c(i) = qekp1(index1(i))
    thpkp1_c(i)= thpkp1(index1(i))
    qpkp1_c(i) = qpkp1(index1(i))
    qsekp1_c(i)= qsekp1(index1(i))
    dqskp1_c(i)= dqskp1(index1(i))
    xsqkp1_c(i)= xsqkp1(index1(i))
    ekp14_c(i) = ekp14(index1(i))
    ekp34_c(i) = ekp34(index1(i))
    pk_c(i)    = pk(index1(i))
    pkp1_c(i)  = pkp1(index1(i))
    exk_c(i)   = exk(index1(i))
    ! Change xsbmin_c to adaptive version
    xsbmin_c(i)= xsbmin_ad(index1(i))
    exkp1_c(i) = exkp1(index1(i))

    bgmk_c(i)  = bgmk(index1(i))
    bgmkp1_c(i)= bgmkp1(index1(i))
    bwkp1_c(i) = bwkp1(index1(i))
  END DO

! -------------------------------------------------------------------
!  detrainment calculation
!
!  SUBROUTINE detrain
!
!  UM Documentation paper 27
!  Section (6)
! -------------------------------------------------------------------

! DEPENDS ON: detrain_4a5a
   CALL detrain_4a5a (ndet,  bwkp1_c,bgmk_c,                              &
                      thek_c,qek_c,thpk_c,qpk_c,qsek_c,dqsk_c,            &
                      thekp1_c,qekp1_c,qsekp1_c,dqskp1_c,                 &
                      ekp14_c,ekp34_c,pk_c,pkp1_c,exk_c,exkp1_c,xsbmin_c, &
                      bgmkp1_c, thpkp1_c,qpkp1_c,xsqkp1_c,                &
                      deltak_c,thrk_c,qrk_c)

  !-----------------------------------------------------------------------
  ! Decompress/expand output arrays from the detrainment calculations
  !-----------------------------------------------------------------------

  DO i=1,ndet
    thpkp1(index1(i)) = thpkp1_c(i)
    qpkp1(index1(i))  = qpkp1_c(i)
  END DO
  DO i=1,ndet
    xsqkp1(index1(i)) = xsqkp1_c(i)

    bgmkp1(index1(i)) = bgmkp1_c(i)
  END DO

END IF     !  ndet =/ 0

DO i=1,npnts
  deltak(i) = 0.0
  thrk(i) = 0.0
  qrk(i) = 0.0
END DO

DO i=1,ndet
  deltak(index1(i)) = deltak_c(i)
  thrk(index1(i))   = thrk_c(i)
  qrk(index1(i))    = qrk_c(i)
END DO

! ----------------------------------------------------------------------
!   Calculate mass flux at level K+1.
!
!   UM Documentation paper 27
!   Section (2B), equation (10A)
! ----------------------------------------------------------------------

DO i=1,npnts
  flxkp1(i) = flxk(i)*(1.0+ekp14(i))*(1.0+ekp34(i))*(1.0-deltak(i))*      &
                                                      (1.0-amdetk(i))
  flxkp12(i)= flxk(i)*(1.0+ekp14(i))*(1.0-deltak(i))* (1.0-amdetk(i))  

END DO

! ---------------------------------------------------------------------
!  Test for points at which convection terminates in layer k+1
!
!  SUBROUTINE TERM_CON
!
!  UM Documentation paper 27
!  Section (7)
! ---------------------------------------------------------------------

! DEPENDS ON: term_con_4a5a
CALL term_con_4a5a (npnts,nlev,k,new_termc,                                 &
                    bwkp1,                                                  &
                    flxkp1,thekp1,qekp1,thpi,qpi,qsekp1,deltak,             &
                    expi,ekp14,ekp34,pstar,pk,pkp1,xsbmin_ad,               &
                    bterm )

! ----------------------------------------------------------------------
!  Cloud microphysics calculation
!
!  SUBROUTINE cloud_w
!
!  UM Documentation paper 27
!  Section (8), (9)
! ----------------------------------------------------------------------

! DEPENDS ON: cloud_w_4a5a
CALL cloud_w_4a5a (k,npnts,xpkp1,qclpkp1,qcfpkp1,prekp1,xsqkp1,blowst,      &
                   flxkp1,xpk,qclpk,qcfpk,thekp1,qekp1,bwkp1,bland,         &
                   qsekp1,bgmkp1,bterm,cca,iccb,icct,tcw,depth,ekp14,ekp34, &
                   delexkp1,cclwp,delpkp1,ccw,lcca,lcbase,lctop,            &
                   l_shallow,l_q_interact,start_lev)

IF (lhook) CALL dr_hook('PARCEL_4A5A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE parcel_4a5a
