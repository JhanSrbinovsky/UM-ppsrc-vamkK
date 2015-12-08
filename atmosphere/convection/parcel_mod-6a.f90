! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Completes lifting of the parcel from layer k to k+1

MODULE parcel_6a_mod

IMPLICIT NONE

!
! Description:
!   Completes lifting of the parcel from layer k to k+1
!   Calls detrain, term_con and cloud_w
!   Detrain  - carries out the forced detrainment.
!   term_con - tests for any convection which is terminating in layer k+1
!   cloud_w  - carries out the cloud microphysics calculation.
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


CONTAINS

SUBROUTINE parcel_6a (k, npnts, np_full, nlev, ntra, ad_on, new_termc,       &
                      start_lev,                                             &
                      pstar, pk, pkp1, delpkp1, exk, exkp1,                  &
                      thek, thekp1, qek, qekp1,                              &
                      qclek, qclekp1, qcfek, qcfekp1,                        &
                      qsek, qsekp1,                                          &
                      thpk, qpk, qclpk, qcfpk,                               &
                      thpi, qpi, expi,                                       &
                      rbuoyk, rbuoykp1, xsbmin,                              &
                      watldek, watldekp1, watldpk, watldpkp1,                &
                      Qlkp1, Qfkp1, Frezkp1,                                 &
                      ekp14, ekp34, amdetk, flxk,                            &
                      uek, uekp1, vek, vekp1,                                &
                      upk, vpk,                                              &
                      traek, traekp1, trapk,                                 &
                      l_q_interact, l_mom_gk, l_tracer,                      &
                      bgmk, bgmkp1, bwk, bwkp1, blowst, bland,               &
                      ! In/out
                      lcbase, lctop,                                         &
                      thpkp1, qpkp1, qclpkp1, qcfpkp1,                       &
                      tcw, depth, cclwp, lcca,                               &
                      bterm, blatent,                                        &
                      ! Out            
                      iccb, icct,                                            &
                      prekp1, thrk, qrk, deltak, flxkp1, flxkp12,            &
                      cca, ccwkp1,                                           &
                      upkp1, vpkp1,                                          &
                      trapkp1)


USE cv_run_mod, ONLY: r_det
USE water_constants_mod, ONLY: lc, lf
USE cv_derived_constants_mod, ONLY: ls
USE cv_param_mod, ONLY: cpress_term
USE atmos_constants_mod, ONLY: cp, repsilon, kappa, rv
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE detrain_6a_mod
USE cloud_w_6a_mod
!USE term_con_6a_mod   Not currently used, rather using 4a5a versions.

IMPLICIT NONE

!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------

INTEGER,INTENT(IN) :: k             ! present model layer
INTEGER,INTENT(IN) :: npnts         ! Number of points
INTEGER,INTENT(IN) :: np_full       ! Full vector length
INTEGER,INTENT(IN) :: nlev          ! Number of model levels for calculations
INTEGER,INTENT(IN) :: ntra          ! Number of tracer variables
INTEGER,INTENT(IN) :: ad_on         ! Flag for adaptive detrainment 
INTEGER,INTENT(IN) :: new_termc     ! Flag for simplified termination of 
                                    ! convection

INTEGER,INTENT(IN) :: start_lev(npnts) ! Level at which convection is initiated

REAL,INTENT(IN) :: pstar(npnts)     ! Surface pressure (Pa)
REAL,INTENT(IN) :: pk(npnts)        ! pressure at mid-point of layer k (Pa)
REAL,INTENT(IN) :: pkp1(npnts)      ! pressure at mid-point of layer k+1 (Pa)
REAL,INTENT(IN) :: delpkp1(npnts)   ! pressure difference across layer k+1 (Pa)
REAL,INTENT(IN) :: exk(npnts)       ! Exner ratio at mid-point of layer k
REAL,INTENT(IN) :: exkp1(npnts)     ! Exner ratio at mid-point of layer k+1
REAL,INTENT(IN) :: thek(npnts)      ! Env. pot. temperature in layer k (K)
REAL,INTENT(IN) :: thekp1(npnts)    ! Env. pot. temperature in layer k+1 (K)
REAL,INTENT(IN) :: qek(npnts)       ! Env. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: qekp1(npnts)     ! Env. spec. humidity in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qclek(npnts)     ! Env. qcl in layer k (kg/kg)
REAL,INTENT(IN) :: qclekp1(npnts)   ! Env. qcl in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qcfek(npnts)     ! Env. qcf in layer k (kg/kg)
REAL,INTENT(IN) :: qcfekp1(npnts)   ! Env. qcf in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qsek(npnts)      ! Env. saturated specific humidity in 
                                    ! in layer k (kg/kg)
REAL,INTENT(IN) :: qsekp1(npnts)    ! Env. saturated specific humidity in 
                                    ! in layer k+1 (kg/kg)
REAL,INTENT(IN) :: thpk(npnts)      ! Par. pot. temperature in layer k (K)
REAL,INTENT(IN) :: qpk(npnts)       ! Par. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: qclpk(npnts)     ! Par. qcl in layer k (kg/kg)
REAL,INTENT(IN) :: qcfpk(npnts)     ! Par. qcf in layer k (kg/kg)
REAL,INTENT(IN) :: thpi(npnts)      ! Initial parcel potential temperature (K)
REAL,INTENT(IN) :: qpi(npnts)       ! Initial parcel specific humidity (kg/kg)
REAL,INTENT(IN) :: expi(npnts)      ! Initial parcel Exner pressure
REAL,INTENT(IN) :: rbuoyk(npnts)    ! Par. buoyancy in layer k (K)
REAL,INTENT(IN) :: rbuoykp1(npnts)  ! Par. buoyancy in layer k+1 (K)
REAL,INTENT(IN) :: xsbmin(npnts)    ! Threshold buoyancy for forced 
                                    ! detrainment (K)
REAL,INTENT(IN) :: watldek(npnts)   ! Env. water loading in layer k (kg/kg)
REAL,INTENT(IN) :: watldekp1(npnts) ! Env. water loading in layer k+1 (kg/kg)
REAL,INTENT(IN) :: watldpk(npnts)   ! Par. water loading in layer k (kg/kg)
REAL,INTENT(IN) :: watldpkp1(npnts) ! Par. water loading in layer k+1 (kg/kg)
REAL,INTENT(IN) :: Qlkp1(npnts)     ! Amount of condensation to liquid water 
                                    ! in the parcel (kg/kg)
REAL,INTENT(IN) :: Qfkp1(npnts)     ! Amount of deposition to ice water
                                    ! in the parcel (kg/kg)
REAL,INTENT(IN) :: Frezkp1(npnts)   ! Amount of freezing from liquid 
                                    ! to ice water in the parcel (kg/kg)
REAL,INTENT(IN) :: ekp14(npnts)     ! Entrainment coefficient at level k+1/4
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: ekp34(npnts)     ! Entrainment coefficient at level k+3/4 
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: amdetk(npnts)    ! Mixing detrainment coefficient at level k 
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: flxk(npnts)      ! Parcel massflux in layer k (Pa/s) 
REAL,INTENT(IN) :: uek(npnts)       ! Env. U in layer k (m/s)
REAL,INTENT(IN) :: uekp1(npnts)     ! Env. U in layer k+1 (m/s)
REAL,INTENT(IN) :: vek(npnts)       ! Env. V in layer k (m/s)
REAL,INTENT(IN) :: vekp1(npnts)     ! Env. V in layer k+1 (m/s)
REAL,INTENT(IN) :: upk(npnts)       ! Par. U in layer k (m/s)
REAL,INTENT(IN) :: vpk(npnts)       ! Par. V in layer k (m/s)
REAL,INTENT(IN) :: traek(np_full,ntra)    ! Env. tracer content 
                                          ! in layer k (kg/kg)
REAL,INTENT(IN) :: traekp1(np_full,ntra)  ! Env. tracer content 
                                          ! in layer k+1 (kg/kg)
REAL,INTENT(IN) :: trapk(np_full,ntra)    ! Par. tracer content 
                                          ! in layer k (kg/kg)

LOGICAL,INTENT(IN) :: l_q_interact  ! True if PC2 is switched on
LOGICAL,INTENT(IN) :: l_mom_gk      ! Switch for inclusion of Gregory-Kershaw
                                    ! CMT
LOGICAL,INTENT(IN) :: l_tracer      ! Switch for tracers

LOGICAL,INTENT(IN) :: bgmk(npnts)   ! mask for parcels which are saturated 
                                    ! in layer k
LOGICAL,INTENT(IN) :: bgmkp1(npnts) ! Mask for parcels which are saturated 
                                    ! in layer k+1
LOGICAL,INTENT(IN) :: bwk(npnts)    ! mask for parcels which have liquid 
                                    ! condensate in layer k
LOGICAL,INTENT(IN) :: bwkp1(npnts)  ! mask for parcels which have liquid 
                                    ! condensate in layer k+1
LOGICAL,INTENT(IN) :: blowst(npnts) ! mask for those points at which stability
                                    ! is low enough for convection to occur
LOGICAL,INTENT(IN) :: bland(npnts)  ! Land/sea mask

!----------------------------------------------------------------------
! Variables which are input and output
!----------------------------------------------------------------------
INTEGER,INTENT(INOUT) :: lcbase(npnts)! Lowest conv. cloud base level
INTEGER,INTENT(INOUT) :: lctop(npnts) ! Lowest conv. cloud top level

REAL,INTENT(INOUT) :: thpkp1(npnts) ! Par. pot. temperature in layer k+1 (K)
                                    ! IN after entrainment and latent heating
                                    ! OUT after forced detrainment
REAL,INTENT(INOUT) :: qpkp1(npnts)  ! Par. spec. humidity in layer k+1 (kg/kg)
                                    ! IN after entrainment and latent heating
                                    ! OUT after forced detrainment
REAL,INTENT(INOUT) :: qclpkp1(npnts)! Par. qcl in layer k+1 (kg/kg)
                                    ! IN after entrainment and latent heating
                                    ! OUT after forced detrainment
REAL,INTENT(INOUT) :: qcfpkp1(npnts)! Par. qcf in layer k+1 (kg/kg)
                                    ! IN after entrainment and latent heating
                                    ! OUT after forced detrainment
REAL,INTENT(INOUT) :: tcw(npnts)    ! Total condensed water (kg/m**2/s)
                                    ! IN summed to layer k
                                    ! OUT summed to layer k+1
REAL,INTENT(INOUT) :: depth(npnts)  ! Depth of convective cloud (m)
                                    ! IN summed to layer k 
                                    ! OUT summed to layer k+1
REAL,INTENT(INOUT) :: cclwp(npnts)  ! Condensed water path (kg/m**2)
                                    ! IN summed to layer k
                                    ! OUT summed to layer k+1
REAL,INTENT(INOUT) :: lcca(npnts)   ! Lowest conv. cloud amount (%)

LOGICAL,INTENT(INOUT) :: bterm(npnts)   ! Mask for parcels which terminate 
                                        ! in layer k+1
LOGICAL,INTENT(INOUT) :: blatent(npnts) ! Mask for points where latent heat has
                                        ! been released

!---------------------------------------------------------------------
! Variables which are output
!---------------------------------------------------------------------
INTEGER,INTENT(OUT) :: iccb(npnts)  ! convective cloud base_level
INTEGER,INTENT(OUT) :: icct(npnts)  ! convective cloud top level

REAL,INTENT(OUT) :: prekp1(npnts)   ! precipitation from parcel as it rises 
                                    ! from layer k to k+1 (kg/m**2/s)
REAL,INTENT(OUT) :: thrk(npnts)     ! pot. temperature of forced detrained
                                    ! parcel in layer k (K)
REAL,INTENT(OUT) :: qrk(npnts)      ! Specific humidity of forced detrained
                                    ! parcel in layer k (kg/kg)
REAL,INTENT(OUT) :: deltak(npnts)   ! Parcel forced detrainment rate in 
                                    ! layer k multiplied by layer thickness
REAL,INTENT(OUT) :: flxkp1(npnts)   ! parcel massflux in layer k+1 (Pa/s)
REAL,INTENT(OUT) :: flxkp12(npnts)  ! parcel massflux in layer k+1/2 (Pa/s)
REAL,INTENT(OUT) :: cca(npnts)      ! convective cloud amount (%)
REAL,INTENT(OUT) :: ccwkp1(npnts)   ! Total condensate in level k+1 (kg/kg)
REAL,INTENT(OUT) :: upkp1(npnts)    ! Par. U in layer k+1 (m/s)
REAL,INTENT(OUT) :: vpkp1(npnts)    ! Par. V in layer k+1 (m/s)
REAL,INTENT(OUT) :: trapkp1(np_full,ntra) ! parcel tracer content 
                                          ! in layer k+1 (kg/kg)

!-------------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

INTEGER :: i               ! loop counter 
INTEGER :: ndet            ! Compress vector length for the detrainment 
                           ! calculation
INTEGER :: ktra            ! Loop counter for tracers

INTEGER :: index1(npnts)   ! Index for compress and expand

REAL :: pk_c(npnts)        ! Compressed pressure at mid-point of layer k (Pa)
REAL :: pkp1_c(npnts)      ! Compressed pressure at mid-point of layer k+1 (Pa)
REAL :: deltak_c(npnts)    ! Compressed Parcel forced detrainment rate in 
                           ! layer k multiplied by layer thickness
REAL :: exk_c(npnts)       ! Compressed Exner ratio at mid-point of layer k
REAL :: exkp1_c(npnts)     ! Compressed Exner ratio at mid-point of layer k+1
REAL :: thek_c(npnts)      ! Compressed Env. pot. temperature in layer k (K)
REAL :: thekp1_c(npnts)    ! Compressed Env. pot. temperature in layer k+1 (K)
REAL :: qek_c(npnts)       ! Compressed Env. spec. humidity in layer k (kg/kg)
REAL :: qekp1_c(npnts)     ! Compressed Env. spec. humidity in layer k+1 (kg/kg)
REAL :: thpk_c(npnts)      ! Compressed Par. pot. temperature in layer k (K)
REAL :: qpk_c(npnts)       ! Compressed Par. spec. humidity in layer k (kg/kg)
REAL :: thpkp1_c(npnts)    ! Compressed Par. pot. temperature in layer k+1 (K)
REAL :: qpkp1_c(npnts)     ! Compressed Par. spec. humidity in layer k+1 (kg/kg)
REAL :: thrk_c(npnts)      ! Compressed pot. temperature of forced detrained
                           ! parcel in layer k (K)
REAL :: qrk_c(npnts)       ! Compressed specific humidity of forced detrained
                           ! parcel in layer k (kg/kg)
REAL :: watldek_c(npnts)   ! compressed Env. water loading in layer k (kg/kg)
REAL :: watldpk_c(npnts)   ! compressed Par. water loading in layer k (kg/kg)
REAL :: watldekp1_c(npnts) ! compressed Env. water loading in layer k+1 (kg/kg)
REAL :: watldpkp1_c(npnts) ! compressed Par. water loading in layer k+1 (kg/kg)
REAL :: Qlkp1_c(npnts)     ! compressed amount of condensation to liquid water 
                           ! in the parcel (kg/kg)
REAL :: Qfkp1_c(npnts)     ! compressed amount of deposition to ice water
                           ! in the parcel (kg/kg)
REAL :: Frezkp1_c(npnts)   ! compressed amount of freezing from liquid 
                           ! to ice water in the parcel (kg/kg)
REAL :: ekp14_c(npnts)     ! Compressed entrainment coefficient at level k+1/4
                           ! multiplied by appropriate layer thickness
REAL :: ekp34_c(npnts)     ! Compressed entrainment coefficient at level k+3/4 
                           ! multiplied by appropriate layer thickness
REAL :: Factor1            ! Factor used in update calculation
REAL :: Factor2            ! Factor used in update calculation
REAL :: xsbmin_ad(npnts)   ! xsbmin adaptive (NOTE, will be different at 
                           ! different points)
REAL :: xsbmin_c(npnts)    ! Compressed threshold buoyancy for forced 
                           ! detrainment (K)

LOGICAL :: bwk_c(npnts)    ! Compressed mask for parcels which have liquid 
                           ! condensate in layer k
LOGICAL :: bwkp1_c(npnts)  ! Compressed mask for parcels which have liquid 
                           ! condensate in layer k+1
LOGICAL :: bgmk_c(npnts)   ! Compressed mask for parcels which are saturated
                           ! in layer k
LOGICAL :: bgmkp1_c(npnts) ! Compressed mask for parcels which are saturated
                           ! in layer k+1
LOGICAL :: bdetk(npnts)    ! Mask for points under going forced detrainment

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('PARCEL_6A',zhook_in,zhook_handle)
!---------------------------------------------------------------------

DO i=1,npnts

! ---------------------------------------------------------------------
!  Calculate mask for those points under going forced detrainment
!
!  UM Documentation paper 27
!  Section (6), equation (23)
! ---------------------------------------------------------------------

  bdetk(i) = rbuoykp1(i)  <   xsbmin(i)

  !  Calculate XSBMIN_AD
  !  Use adaptive detrainment if ADAPTIVE ON flag is set to 1
  !  making sure that adaptive min. buoyancy doesn't fall below 0.0
  IF (ad_on  ==  1) THEN
    xsbmin_ad(i) = r_det * rbuoyk(i) + (1-r_det)*rbuoykp1(i)
    xsbmin_ad(i) = MAX(xsbmin_ad(i),0.0)
  ELSE
    xsbmin_ad(i) = xsbmin(i)
  END IF

  !  Set mask of points for forced detrainment
  bdetk(i) = rbuoykp1(i)  <   xsbmin_ad(i)

END DO


!Initialise properties of the forced detrained parcel
DO i=1,npnts
  deltak(i) = 0.0
  thrk(i)   = 0.0
  qrk(i)    = 0.0
END DO

!----------------------------------------------------------------------
! Compress all input arrays for the forced detrainment calculations
! NB the compression is in the same order as the argument list
!----------------------------------------------------------------------

ndet = 0
DO i=1,npnts
  IF (bdetk(i)) THEN
    ndet = ndet + 1
    index1(ndet) = i
  END IF
END DO

IF (ndet > 0) THEN
  DO i=1,ndet
    !Compression for INTENT(IN)  
    pk_c(i)        = pk(index1(i))
    pkp1_c(i)      = pkp1(index1(i))
    exk_c(i)       = exk(index1(i))
    exkp1_c(i)     = exkp1(index1(i))
    thek_c(i)      = thek(index1(i))
    thekp1_c(i)    = thekp1(index1(i))
    qek_c(i)       = qek(index1(i))
    qekp1_c(i)     = qekp1(index1(i))
    thpk_c(i)      = thpk(index1(i))
    qpk_c(i)       = qpk(index1(i))
    watldek_c(i)   = watldekp1(index1(i))
    watldpk_c(i)   = watldpkp1(index1(i))
    watldekp1_c(i) = watldekp1(index1(i))
    watldpkp1_c(i) = watldpkp1(index1(i))
    Qlkp1_c(i)     = Qlkp1(index1(i))
    Qfkp1_c(i)     = Qfkp1(index1(i))
    Frezkp1_c(i)   = Frezkp1(index1(i))      
    ekp14_c(i)     = ekp14(index1(i))
    ekp34_c(i)     = ekp34(index1(i))
    xsbmin_c(i)    = xsbmin_ad(index1(i))
    bwk_c(i)       = bwk(index1(i))
    bwkp1_c(i)     = bwkp1(index1(i))
    bgmk_c(i)      = bgmk(index1(i))
    bgmkp1_c(i)    = bgmkp1(index1(i))
    !Compression for INTENT(INOUT)  
    thpkp1_c(i)    = thpkp1(index1(i))
    qpkp1_c(i)     = qpkp1(index1(i))
  END DO

! -------------------------------------------------------------------
!  detrainment calculation
!
!  SUBROUTINE detrain
!
!  UM Documentation paper 27
!  Section (6)
! -------------------------------------------------------------------


   CALL detrain_6a (ndet, pk_c, pkp1_c, exk_c, exkp1_c,                 &
                       thek_c, thekp1_c, qek_c, qekp1_c,                &
                       thpk_c, qpk_c,                                   &
                       watldek_c, watldpk_c, watldekp1_c, watldpkp1_c,  & 
                       Qlkp1_c, Qfkp1_c, Frezkp1_c,                     &
                       ekp14_c, ekp34_c, xsbmin_c,                      &
                       bwk_c, bwkp1_c, bgmk_c, bgmkp1_c,                &
                       thpkp1_c, qpkp1_c,                               &
                       deltak_c, thrk_c, qrk_c)

  !-----------------------------------------------------------------------
  ! Decompress/expand output arrays from the detrainment calculations
  !-----------------------------------------------------------------------

  DO i=1,ndet
    !Decompression for INTENT(INOUT)
    thpkp1(index1(i)) = thpkp1_c(i)
    qpkp1(index1(i))  = qpkp1_c(i)
    !Decompression for INTENT(OUT)
    deltak(index1(i)) = deltak_c(i)
    thrk(index1(i))   = thrk_c(i)
    qrk(index1(i))    = qrk_c(i)
  END DO
  
END IF     !  ndet > 0


! ----------------------------------------------------------------------
! Update parcel properties at level k+1. NB precipitation is calculated
! later.
! ----------------------------------------------------------------------

DO i=1,npnts
  flxkp1(i)  = (1.0-amdetk(i))*(1.0-deltak(i))                        &
               *(1.0+ekp14(i))*(1.0+ekp34(i))*flxk(i) 
END DO

! ---------------------------------------------------------------------
!  Test for points at which convection terminates in layer k+1 due to 
!  mass flux becoming very small or convection reaching the top of
!  the model.
!  Currently includes the old ropey undilute parcel test because it is 
!  still used for shallow convection.
! ---------------------------------------------------------------------

! DEPENDS ON: term_con_4a5a
CALL term_con_4a5a (npnts,nlev,k,new_termc,                           &
               bwkp1,                                                 &
               flxkp1,thekp1,qekp1,thpi,qpi,qsekp1,deltak,            &
               expi,ekp14,ekp34,pstar,pk,pkp1,xsbmin_ad,              &
               bterm )
               
!CALL term_con_6a(npnts,nlev,k,flxkp1,ekp14,ekp34,pstar,bterm)

DO i=1,npnts
  IF (deltak(i) >= 0.95 .OR. bterm(i)) THEN
! ---------------------------------------------------------------------
!   The parcel is terminating so calculate the 
!   the properties at level k+1 for a terminating parcel
! ---------------------------------------------------------------------
    bterm(i)   = .TRUE.
    flxkp1(i)  = 0.0
    flxkp12(i) = 0.0
    thpkp1(i)  = 0.0
    qpkp1(i)   = 0.0
    qclpkp1(i) = 0.0
    qcfpkp1(i) = 0.0
    thrk(i)    = thpk(i)
    qrk(i)     = qpk(i)
    deltak(i)  = 1.0
    IF (l_tracer) THEN
      DO ktra = 1,ntra
         trapkp1(i,ktra) = 0.0
      END DO
    END IF
    IF (l_mom_gk) THEN                         
      upkp1(i)   = 0.0
      vpkp1(i)   = 0.0
    END IF       ! l_mom_gk test


  ELSE
! ---------------------------------------------------------------------
!   The parcel does not terminate if the forced detrainment rate is 
!   not too large and the parcel mass flux is not to small. The parcel
!   properties at level k+1 are therefore calculated.
! ---------------------------------------------------------------------
!    bterm(i)   = .FALSE.
    Factor1     = 1.0/((1.0+ekp14(i))*(1.0+ekp34(i)))
    Factor2     = Factor1/(1.0-deltak(i))
!   flxkp1 is calculated earlier
    flxkp12(i)  = (1.0-amdetk(i))*(1.0-deltak(i))*(1.0+ekp14(i))*flxk(i) 
    thpkp1(i)   = ( thpk(i) - deltak(i)*thrk(i)                           &
                + (1.0-deltak(i))*ekp14(i)*thek(i)                        &
                + (1.0-deltak(i))*(1.0+ekp14(i))*ekp34(i)*thekp1(i) )     &
                * Factor2                                                 &
                + (lc*Qlkp1(i) + ls*Qfkp1(i) + lf*Frezkp1(i))/(cp*exkp1(i))

    qpkp1(i)    = ( qpk(i) - deltak(i)*qrk(i)                             &
                + (1.0-deltak(i))*ekp14(i)*qek(i)                         &
                + (1.0-deltak(i))*(1.0+ekp14(i))*ekp34(i)*qekp1(i) )      &
                * Factor2                                                 &
                - Qlkp1(i) - Qfkp1(i)

    IF (l_q_interact) THEN
    !PC2 so entrainment from the environment
      qclpkp1(i) = ( qclpk(i)                                             &
                 + ekp14(i)*qclek(i)                                      &
                 + (1.0+ekp14(i))*ekp34(i)*qclekp1(i) )                   &
                 * Factor1                                                &
                 + Qlkp1(i) - Frezkp1(i)
      qcfpkp1(i) = ( qcfpk(i)                                             &
                 + ekp14(i)*qcfek(i)                                      &
                 + (1.0+ekp14(i))*ekp34(i)*qcfekp1(i) )                   &
                 * Factor1                                                &
                 + Qfkp1(i) + Frezkp1(i)
    ELSE
    !Not PC2 so no entrainment from the environment
      qclpkp1(i) = qclpk(i)                                               &
                 * Factor1                                                &
                 + Qlkp1(i) - Frezkp1(i)
      qcfpkp1(i) = qcfpk(i)                                               &
                 * Factor2                                                &
                 + Qfkp1(i) + Frezkp1(i)  
    END IF ! l_q_interact

    !Due to floating point arithmetic, qclpkp1 or qcfpkp1 could 
    !be tiny and/or negative when they should be 0.0. This sometimes 
    !happens at phase changes. Therefore reset very small values to zero.
    IF (ABS(qclpkp1(i)) .LT. 1.0E-18) qclpkp1(i) = 0.0
    IF (ABS(qcfpkp1(i)) .LT. 1.0E-18) qcfpkp1(i) = 0.0
    
    blatent(i) = (blatent(i) .OR. Qlkp1(i) > 0.0 .OR. Qfkp1(i) > 0.0)

    IF (l_mom_gk) THEN   ! l_mom_gk set according to type of CMT required
                         ! This code does Gregory-Kershaw CMT
      upkp1(i)   = ( upk(i)                                               &
                 + ekp14(i)*uek(i)                                        &
                 + (1.0+ekp14(i))*ekp34(i)*uekp1(i) )                     &
                 * Factor1                                                &
                 - cpress_term*(uek(i)-uekp1(i))/(1.0+ekp34(i))
                 
      vpkp1(i)   = ( vpk(i)                                               &
                 + ekp14(i)*vek(i)                                        &
                 + (1.0+ekp14(i))*ekp34(i)*vekp1(i) )                     &
                 * Factor1                                                &
                 - cpress_term*(vek(i)-vekp1(i))/(1.0+ekp34(i))

    END IF       ! l_mom_gk test


    IF (l_tracer) THEN
      DO ktra = 1,ntra
         trapkp1(i,ktra) = ( trapk(i,ktra)                                &
                         + ekp14(i)*traek(i,ktra)                         &
                         + (1.0+ekp14(i))*ekp34(i)*traekp1(i,ktra) )      &
                         * Factor1

      END DO
    END IF
  END IF
END DO


! ----------------------------------------------------------------------
!  Cloud microphysics calculation
!
!  SUBROUTINE cloud_w
!
!  UM Documentation paper 27
!  Section (8), (9)
! ----------------------------------------------------------------------

CALL cloud_w_6a (k, npnts, start_lev,                                &
                 flxkp1, qclpk, qcfpk,                               &
                 thekp1, qekp1, qsekp1,                              &
                 ekp14, ekp34, delpkp1,                              &
                 blowst, bwkp1, bland, bterm, l_q_interact,          &
                 lcbase, lctop,                                      &
                 qclpkp1, qcfpkp1,                                   &
                 tcw, depth, cclwp, lcca,                            &
                 iccb, icct, prekp1, cca, ccwkp1)

IF (lhook) CALL dr_hook('PARCEL_6A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE parcel_6a

END MODULE parcel_6a_mod
