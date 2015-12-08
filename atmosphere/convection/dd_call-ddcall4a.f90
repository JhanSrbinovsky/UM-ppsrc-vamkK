! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate initial Downdraught massflux.
!
! Subroutine Interface:
SUBROUTINE dd_call_4a5a(npnts, nterm, kct, nlev, trlev, ntra              &
                  ,iccb, icct, index1                                     &
                  ,l_tracer                                               &
                  ,bwater                                                 &
                  ,exner_layer_centres, exner_layer_boundaries            &
                  ,p_layer_centres, p_layer_boundaries, pstar             &
                  ,recip_pstar, timestep , cca                            &
                  ,thp, qp, the, qe, qse, trap, trae, flx                 &
                  ,precip, dthbydt, dqbydt, dtrabydt                      &
                  ,rain, snow ,rain_3d, snow_3d                           &
                  ,dd_flux, entrain_dwn, detrain_dwn)

USE cv_stash_flg_mod, ONLY:                                               &
  flg_dwn_flx, flg_entr_dwn, flg_detr_dwn  

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! 
! Description: Calculate initial Downdraught massflux.
!            Reset en/detrainment rates for Downdraught
!            Compress/expand variables
!            Initialise downdrought
!            Call downdraught routine
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

INTEGER, INTENT(IN) :: &
  npnts                &  ! Full vector length
 ,nterm                &  ! Compressed vector length for downdraught cal.
 ,kct                  &  ! Convective cloud top layer number 
                          ! set required here equal to kmax_term+1).
 ,nlev                 &  ! Number of model levels
 ,trlev                &  ! Number of tracer levels
 ,ntra                    ! Number of tracers

INTEGER, INTENT(IN) :: &
  iccb(npnts)          &  ! Convective cloud base level (m)
 ,icct(npnts)          &  ! Convective cloud top level (m)
 ,index1(npnts)           ! index of points where DOwndraught possible

LOGICAL, INTENT(IN) :: &
  l_tracer                ! Switch for tracers

LOGICAL, INTENT(IN) ::  &
  bwater(npnts,2:kct+1)   ! Mask for points at which condensate is liquid


REAL, INTENT(IN) ::                      &
  exner_layer_centres(npnts,0:kct+1)     & ! exner pressure
 ,exner_layer_boundaries(npnts,0:kct+1)  & ! exner at half level above
                                           ! exner_layer_centres
 ,p_layer_centres(npnts,0:kct+2)         & ! Pressure (Pa)
 ,p_layer_boundaries(npnts,0:kct+1)      & ! Pressure at half level above
                                           ! p_layer_centres (Pa)
 ,pstar(npnts)                           & ! Surface pressure (Pa)
 ,recip_pstar(npnts)                       ! 1/pstar (Pa)

REAL, INTENT(IN) ::         &
  timestep                    ! timestep

REAL, INTENT(IN) ::      &
  cca(npnts)             & ! 2d convective cloud amount
 ,thp(npnts,kct+1)       & ! Parcel potential temperature (K)
 ,qp(npnts,kct+1)        & ! Parcel mixing ratio (kg/kg)
 ,the(npnts,kct+1)       & ! Model enviromental potential temperature (K)
 ,qe(npnts,kct+1)        & ! Model enviromental mixing ratio (kg/kg)
 ,qse(npnts,kct+1)       & ! Model enviromental qsat (kg/kg)
 ,trap(npnts,nlev,ntra)  & ! Parcel tracer (kg/kg)
 ,trae(npnts,trlev,ntra) & ! Environment tracer (kg/kg)
 ,flx(npnts,kct+1)         ! updraught mass flux (Pa/s)

REAL, INTENT(INOUT) ::      &
  precip(npnts,kct+1)       & ! precipitation added when descending
                              ! from layer k to k-1 (kg/m**2/s)
 ,dthbydt(npnts,kct+1)      & ! increment to model potential temperature (K/s)
 ,dqbydt(npnts,kct+1)       & ! increment to model mixing ratio (kg/kg/s)
 ,dtrabydt(npnts,nlev,ntra) & ! increment to model tracers (kg/kg/s)
 ,rain(npnts)               & ! rainfall at surface (kg/m**2/s)
 ,snow(npnts)               & ! snowfall at surface (kg/m**2/s)
 ,rain_3d(npnts,kct+1)      & ! rainfall flux  (kg/m**2/s)
 ,snow_3d(npnts,kct+1)        ! snowfall flux  (kg/m**2/s)

REAL, INTENT(OUT) ::        &
  entrain_dwn(npnts,kct+1)  & ! fractional entrainment rate for downdraught
                              ! mass flux
 ,detrain_dwn(npnts,kct+1)  & ! fractional detrainment rate for downdraught
                              ! mass flux
 ,dd_flux(npnts,kct+1)        ! Downdraught mass flux (Pa/s)

!-----------------------------------------------------------------------
! Local variables 
!-----------------------------------------------------------------------

INTEGER ::       &
  i,k ,i2,ktra   & ! Loop counters
 ,ndd            & ! Compressed vector length for downdraught calculation
 ,nddon_tmp        ! Number of points with active downdraught


LOGICAL ::            &
  bddi(npnts)         & ! Mask for points where downdraught might occur
 ,bdd_start(npnts)    & ! mask for those points where downdraught is able
                        ! to start from level k
 ,bddwt_k(npnts)      & ! mask for points in downdraught where ppt in 
                        ! layer k is liquid
 ,bddwt_km1(npnts)    & ! mask for points in downdraught where ppt in 
                        ! layer k-1 is liquid
 ,bdd_on(npnts)         ! mask for those points where DD continues from
                        ! layer k+1

REAL ::               &
  flx_strt(npnts)       ! Mass flux at level where downdraught starts (Pa/s)

!-----------------------------------------------------------------------
! Compressed arrays
!-----------------------------------------------------------------------

INTEGER ::        &
  iccb_c(nterm) & ! Compressed cloud base level
 ,kmin(nterm)     ! Freezing level where entrainment rates are increased


LOGICAL ::             &
  bwater_k_c(nterm)    & ! Compressed mask for those points at which condensate
                         ! is water in layer k
 ,bddi_c(nterm)        & ! Compressed mask for points where downdraught may
                         ! initiate
 ,bdd_start_c(nterm)   & ! Compressed mask for those points where downdraught
                         ! is able to start from level k
 ,bddwt_k_c(nterm)     & ! Compressed mask for points in DD where ppt in
                         ! layer k is liquid
 ,bddwt_km1_c(nterm)   & ! Compressed mask for points in DD where ppt in
                         ! layer k-1 is liquid
 ,bdd_on_c(nterm)        ! Compressed mask for points where DD continues from
                         ! layer K+1

REAL ::                &
  exner_km12_c(nterm)  & ! Compressed exner function at layer k
 ,exner_kp12_c(nterm)  & ! Compressed exner function at layer k+1
 ,exner_KM32_c(nterm)  & ! Compressed exner function at layer k-1

 ,pk(nterm)            & ! Pressure of layer k (Pa)
 ,p_km1(nterm)         & ! Pressure of layer k-1 (Pa)
 ,exk(nterm)           & ! exner ratio for layer K
 ,exkm1(nterm)         & ! exner ratio for layer K-1

 ,delpk(nterm)         & ! Pressure difference across layer K  (Pa)
 ,delpkm1(nterm)       & ! Pressure difference across layer K-1 (Pa)

 ,amdetk(nterm)        & ! Mixing detrainment at level k multiplied by 
                         ! appropriate layer thickness
 ,ekm14(nterm)         & ! exner ratio at layer k-1/4
 ,ekm34(nterm)         & ! exner ratio at layer k-3/4

 ,precip_k_c(nterm)    & ! Compressed precipitation added when descending 
                         ! from layer K to K-1 (kg/m**2/s)
 ,q_k_c(nterm)         & ! Compressed parcel mixing ratio of layer K (kg/kg)

 ,th_k_c(nterm)        & ! Compressed parcel potential temperature of
                         ! layer k (K)
 ,tra_k_c(nterm,ntra)  & ! Compressed parcel tracer in layer K (kg/kg)

 ,pstar_c(nterm)       & ! Compressed surface pressure (Pa)

 ,recip_pstar_c(nterm)   ! Reciprocal of comp. pstar array

REAL ::                                 &
  P_layer_centres_c(nterm,0:kct+2)      & ! Pressure (Pa)
 ,P_layer_boundaries_c(nterm,0:kct+1)   & ! Pressure (Pa)
 ,exner_layer_centres_c(nterm,0:kct+1)    ! exner

REAL ::                 &
  dthbydt_K_c(nterm)    & ! Compressed increment to model potential 
                          ! temperature of layer k (K/s)
 ,dthbydt_km1_c(nterm)  & ! Compressed increment to model potential
                          ! temperature of layer k-1 (K/s)
 ,dqbydt_k_c(nterm)     & ! Compressed increment to model mixing ratio
                          ! of layer k (kg/kg/s)
 ,dqbydt_km1_c(nterm)   & ! Compressed increment to model mixing ratio 
                          ! of layer k-1(kg/kg/s)
 ,dtra_k_c(nterm,ntra)  & ! Compressed increment to model tracer of
                          ! layer k (kg/kg/s)
 ,dtra_km1_c(nterm,ntra)  ! Compressed increment to model tracer of 
                          ! layer k-1 (kg/kg/s)

REAL ::                 &
  deltd(nterm)          & ! Cooling necessary to achieve saturation (K)

 ,delqd(nterm)          & ! moistening necessary to achieve saturation (kg/kg)

 ,deltrad(nterm,ntra)   & ! Depletion of environment tracer due to 
                          ! downdraught formation (kg/kg)
 ,qdd_k(nterm)          & ! mixing ratio of downdraught in layer K (kg/kg)

 ,thdd_k(nterm)         & ! Model potential temperature of downdraught 
                          ! in layer k (K)
 ,tradd_k(nterm,ntra)     ! Model tracer of downdraught in layer K (kg/kg)

REAL ::                  &
  flx_dd_k(npnts)        & ! Downdraught initial mass flux (Pa/s)

 ,flx_dd_k_c(nterm)      & ! Compressed downdraught initial mass flux (Pa/s)

 ,qe_k_c(nterm)          & ! Compressed environment mixing ratio of 
                           ! layer k (kg/kg)
 ,qe_km1_c(nterm)        & ! Compressed environment mixing ratio of 
                           ! layer k-1 (kg/kg)
 ,qse_k_c(nterm)         & ! Compressed environment qsat mixing ratio of 
                           ! layer k (kg/kg)
 ,qse_km1_c(nterm)       & ! Compressed environment qsat mixing ratio of 
                           ! layer k-1 (kg/kg)
 ,the_k_c(nterm)         & ! Compressed potential temperature
                           ! of environment in layer k (K)
 ,the_km1_c(nterm)       & ! Compressed potential temperature
                           ! of environment in layer k-1 (K)
 ,trae_k_c(nterm,ntra)   & ! Compressed tracer of environment in layer k
                           ! (kg/kg)
 ,trae_km1_c(nterm,ntra) & ! Compressed tracer of environment in layer k-1
                           ! (kg/kg)
 ,rain_c(nterm)          & ! Compressed surface rainfall (kg/m**2/s)

 ,snow_c(nterm)          & ! Compressed surface snowfall (kg/m**2/s)

 ,flx_ud_k_c(nterm)      & ! updraught mass flux at layer K

 ,rain_env(nterm)        & ! Amount of rainfall passing through environment 
                           ! (kg/m**2/s)
 ,snow_env(nterm)        & ! Amount of snowfall passing through environment 
                           ! (kg/m**2/s)
 ,rain_dd(nterm)         & ! Amount of rainfall passing through
                           ! downdraught (KG/M**2/S)
 ,snow_dd(nterm)         & ! Amount of snowfall passing through
                           ! downdraught (KG/M**2/S)
 ,flx_strt_c(nterm)      & ! Compressed value of flx_strt (Pa/s)

 ,cca_c(nterm)           & ! Compressed convective cloud amount

 ,lr_ud_ref(nterm)         ! precipitation mixing ratio at lowest 
                           ! precipitationg level of UD

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! Compression to Downdraught points  (all levels even above term level)
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('DD_CALL_4A5A',zhook_in,zhook_handle)

ndd=nterm

DO k=0,kct+1
  DO i=1,ndd
    p_layer_centres_c(i,k)     = p_layer_centres(index1(i),k)
    exner_layer_centres_c(i,k) = exner_layer_centres(index1(i),k)
    p_layer_boundaries_c(i,k)  = p_layer_boundaries(index1(i),k)
  END DO
END DO

!----------------------------------------------------------------------
! Iniialise logical arrays as .FALSE.
!-----------------------------------------------------------------------

DO i=1,npnts
  bddi(i)      = .FALSE.
  bdd_start(i) = .FALSE.
  bddwt_k(i)   = .FALSE.
  bddwt_km1(i) = .FALSE.
  bdd_on(i)    = .FALSE.
END DO

DO i2=1,ndd
  bddi(index1(i2)) = .TRUE.    ! all points with downdraught
END DO

!----------------------------------------------------------------------
! Calculate initial downdraught mass flux
!-----------------------------------------------------------------------

! DEPENDS ON: flx_init_4a5a
CALL flx_init_4a5a(npnts, kct, iccb, icct, flx, flx_dd_k, bddi, flx_strt)

!-----------------------------------------------------------------------
! Compress all input arrays for the downdraught calculation down to
! those points where a downdraught is possible when the convection
! terminates in the column.
!-----------------------------------------------------------------------

DO k = kct+1,2,-1

  DO i=1,ndd
    th_k_c(i) = thp(index1(i),k)
    q_k_c(i)  = qp(index1(i),k)
    the_k_c(i)   = the(index1(i),k)
    the_km1_c(i) = the(index1(i),k-1)
    qe_k_c(i)   = qe(index1(i),k)
    qe_km1_c(i) = qe(index1(i),k-1)
    qse_k_c(i)   = qse(index1(i),k)
    qse_km1_c(i) = qse(index1(i),k-1)
    dthbydt_k_c(i)   = dthbydt(index1(i),k)
    dthbydt_km1_c(i) = dthbydt(index1(i),k-1)
    dqbydt_k_c(i)   = dqbydt(index1(i),k)
    dqbydt_km1_c(i) = dqbydt(index1(i),k-1)
    exner_km12_c(i) = exner_layer_boundaries(index1(i),k-1)
    exner_kp12_c(i) = exner_layer_boundaries(index1(i),k)
    exner_km32_c(i) = exner_layer_boundaries(index1(i),k-2)
    precip_k_c(i) = precip(index1(i),k)
    flx_ud_k_c(i) = flx(index1(i),k)
    bwater_k_c(i) = bwater(index1(i),k)
  END DO

  IF (l_tracer) THEN   ! If run has tracers

    DO ktra=1,ntra
      DO i=1,ndd
        tra_k_c(i,ktra)    = trap(index1(i),k,ktra)
        trae_k_c(i,ktra)   = trae(index1(i),k,ktra)
        trae_km1_c(i,ktra) = trae(index1(i),k-1,ktra)
        dtra_k_c(i,ktra)   = dtrabydt(index1(i),k,ktra)
        dtra_km1_c(i,ktra) = dtrabydt(index1(i),k-1,ktra)
      END DO
    END DO

  END IF

  IF (k == kct+1) THEN

    DO i=1,ndd
      flx_dd_k_c(i) = flx_dd_k(index1(i))
      flx_strt_c(i) = flx_strt(index1(i))
      pstar_c(i)      = pstar(index1(i))
      recip_pstar_c(i)= recip_pstar(index1(i))
      iccb_c(i) = iccb(index1(i))
      bddi_c(i) = bddi(index1(i))
      bdd_start_c(i) = bdd_start(index1(i))
      rain_c(i) = rain(index1(i))
      snow_c(i) = snow(index1(i))
      bddwt_k_c(i)   = bddwt_k(index1(i))
      bddwt_km1_c(i) = bddwt_km1(index1(i))
      bdd_on_c(i) = bdd_on(index1(i))
      cca_c(i) = cca(index1(i))
      lr_ud_ref(i) = 0.0
    END DO
  END IF

!----------------------------------------------------------------------
! If below convective cloud base then downdraught not allowed to form
!----------------------------------------------------------------------

  DO i=1,ndd
    IF (k <  iccb_c(i)) bddi_c(i)=.FALSE.
  END DO

!-----------------------------------------------------------------------
! Reset en/detrainment rates for downdraught
!-----------------------------------------------------------------------

! DEPENDS ON: layer_dd_4a5a
  CALL layer_dd_4a5a (ndd,k,kct,the_k_c,the_km1_c,flx_strt_c,           &
                     p_layer_centres_c,p_layer_boundaries_c,            &
                     exner_layer_centres_c,                             &
                     exner_km12_c,exner_kp12_c,                         &
                     exner_km32_c,pstar_c,pk,p_km1,delpk,delpkm1,exk,   &
                     exkm1,amdetk,ekm14,ekm34,kmin,bddi_c,              &
                     recip_pstar_c)

!----------------------------------------------------------------------
! If level k within 150mb of surface then downdraught not allowed to
! form
!----------------------------------------------------------------------

  DO i=1,ndd
    IF (pk(i) >  (pstar_c(i)-15000.0)) bddi_c(i)=.FALSE.
  END DO

!-----------------------------------------------------------------------
! Initialise downdraught
! downdraught not allowed to form from cloud top layer (kct+1)
! or from below cloud base
!-----------------------------------------------------------------------

! DEPENDS ON: dd_init
  CALL dd_init(ndd,nterm,th_k_c,q_k_c,the_k_c,qe_k_c,qse_k_c,pk,exk,  &
                   thdd_k,qdd_k,deltd,delqd,bdd_start_c,k,bddi_c,     &  
                   bdd_on_c,l_tracer,ntra,tra_k_c,                    &  
                   trae_k_c,tradd_k,deltrad)

!-----------------------------------------------------------------------
! Update mask for where downdraught occurs
!-----------------------------------------------------------------------

  DO i=1,ndd
    IF (bdd_start_c(i).or.bdd_on_c(i)) bdd_on_c(i)=.TRUE.
  END DO

!
! If downdraught initiated set diagnostic array
!
!

  If(flg_dwn_flx) THEN
    DO i=1,ndd
      IF(bdd_start_c(i)) dd_flux(index1(i),k)=flx_dd_k(index1(i))
    END DO
  END IF

  nddon_tmp = 0
  DO i=1,ndd
    IF (bdd_on_c(i)) THEN
      nddon_tmp = nddon_tmp+1
    END IF
  END DO

!-----------------------------------------------------------------------
! Call downdraught routine
!-----------------------------------------------------------------------
        
! DEPENDS ON: downd_4a5a
  CALL downd_4a5a(ndd,nterm,k,kct,ntra,nddon_tmp, iccb_c,               &
             l_tracer,bwater_k_c,                                       &
             timestep,the_k_c,the_km1_c,                                &
             qe_k_c,qe_km1_c,qse_km1_c,p_km1,delpk,delpkm1,exk,         &
             exkm1,deltd,delqd,amdetk,ekm14,ekm34,flx_ud_k_c,cca_c,     &
             trae_k_c,trae_km1_c,deltrad,                               &
             bdd_start_c,bddwt_k_c,bddwt_km1_c,bdd_on_c,                &
             thdd_k,qdd_k,dthbydt_k_c,dthbydt_km1_c,dqbydt_k_c,         &
             dqbydt_km1_c,rain_c,snow_c,precip_k_c,rain_env,snow_env,   &
             rain_dd,snow_dd,flx_dd_k_c,lr_ud_ref,                      &
             tradd_k,dtra_k_c,dtra_km1_c)

!-----------------------------------------------------------------------
! Decompress/expand those variables which are to be output
!-----------------------------------------------------------------------

  DO i=1,ndd
    dthbydt(index1(i),k)   = dthbydt_k_c(i)
    dthbydt(index1(i),k-1) = dthbydt_km1_c(i)
    dqbydt(index1(i),k)    = dqbydt_k_c(i)
    dqbydt(index1(i),k-1)  = dqbydt_km1_c(i)
  END DO
!
! Need to check that point would be selected in S.R DOWND or else
! not sensible to set entrainment and detrainment rates  in diagnostics
!
  IF(flg_dwn_flx) THEN
    DO i=1,ndd
      IF(bdd_on_c(i)) THEN
         dd_flux(index1(i),k-1) = flx_dd_k_c(i)
      END IF
    END DO
  END IF

  IF(flg_entr_dwn) THEN
    DO i=1,ndd
      IF(bdd_on_c(i)) THEN
        entrain_dwn(index1(i),k)=(1.0-amdetk(i))*                          &
                                 (ekm14(i)+ekm34(i)*(1.0+ekm14(i)))*       &
                                   dd_flux(index1(i),k)
      END IF
    END DO
  END IF

  IF(flg_detr_dwn) THEN
    DO i=1,ndd
      IF(bdd_on_c(i)) THEN
        detrain_dwn(index1(i),k)=-amdetk(i)*dd_flux(index1(i),k)
      END IF
    END DO
  END IF

  IF (k == 2) THEN
    DO i=1,ndd
      rain(index1(i)) = rain_c(i)
      snow(index1(i)) = snow_c(i)
    END DO
  END IF


  DO i=1, ndd
    precip(index1(i),k)    = precip_k_c(i)
    rain_3d(index1(i),k-1) = rain_3d(index1(i),k-1) + rain_dd(i) + rain_env(i) 
    snow_3d(index1(i),k-1) = snow_3d(index1(i),k-1) + snow_dd(i) + snow_env(i)
  END DO

  IF (l_tracer) THEN       ! Runs with tracers

    DO ktra=1,ntra
      DO i=1,ndd
        dtrabydt(index1(i),k,ktra)   = dtra_k_c(i,ktra)
        dtrabydt(index1(i),k-1,ktra) = dtra_km1_c(i,ktra)
      END DO
    END DO

  END IF

!----------------------------------------------------------------------
!   End of main K loop
!----------------------------------------------------------------------

END DO    ! End of main level loop

IF (lhook) CALL dr_hook('DD_CALL_4A5A',zhook_out,zhook_handle)
RETURN

END SUBROUTINE dd_call_4a5a
