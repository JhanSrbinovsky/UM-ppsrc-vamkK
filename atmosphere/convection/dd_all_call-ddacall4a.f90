! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate initial Downdraught massflux.
!
! Subroutine Interface:
SUBROUTINE dd_all_call_4a5a (npnts, npossdd, klev, nlev, trlev, ntra    &
                       ,kterm, iccb, icct, index1                       &
                       ,l_tracer                                        &
                       ,bwater                                          &
                       ,exner_layer_centres, exner_layer_boundaries     &
                       ,p_layer_centres, p_layer_boundaries, pstar      &
                       ,recip_pstar, timestep , cca                     &
                       ,thp, qp, the, qe, qse, trap, trae, flx          &
                       ,precip, dthbydt, dqbydt, dtrabydt               &
                       ,rain, snow ,rain_3d, snow_3d                    &
                       ,dd_flux, entrain_dwn, detrain_dwn)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE cv_stash_flg_mod, ONLY:                                             &
  flg_dwn_flx, flg_entr_dwn, flg_detr_dwn

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
  npnts                &  ! Vector length of full arrays
 ,npossdd              &  ! Length of array holding index of points where a
                          ! downdraught is possible.
 ,klev                 &  ! Number of levels (may be model levels or a reduced 
                          ! set required here equal to kmax_term+1).
 ,nlev                 &  ! Number of model levels
 ,trlev                &  ! Number of tracer levels
 ,ntra                    ! Number of tracers


INTEGER, INTENT(IN) :: &
  kterm(npnts)         &  ! Convective cloud top layer
 ,iccb(npnts)          &  ! Convective cloud base level (m)
 ,icct(npnts)          &  ! Convective cloud top level (m)
 ,index1(npnts)           ! index of points where downdraught possible

LOGICAL, INTENT(IN) :: &
  l_tracer                ! Switch for tracers

LOGICAL, INTENT(IN) :: &
  bwater(npnts,2:klev+1)  ! Mask for points at which condensate is liquid

REAL, INTENT(IN) ::                      &
  exner_layer_centres(npnts,0:klev+1)    & ! exner pressure
 ,exner_layer_boundaries(npnts,0:klev+1) & ! exner at half level above
                                           ! exner_layer_centres
 ,p_layer_centres(npnts,0:klev+2)        & ! Pressure (Pa)
 ,p_layer_boundaries(npnts,0:klev+1)     & ! Pressure at half level above
                                           ! p_layer_centres (Pa)
 ,pstar(npnts)                           & ! Surface pressure (Pa)
 ,recip_pstar(npnts)                       ! 1/pstar (Pa)

REAL, INTENT(IN) ::    &
  timestep                ! timestep

REAL, INTENT(IN) ::      &
  cca(npnts)             & ! 2d convective cloud amount
 ,thp(npnts,klev+1)      & ! Parcel potential temperature (K)
 ,qp(npnts,klev+1)       & ! Parcel mixing ratio (kg/kg)
 ,the(npnts,klev+1)      & ! Model enviromental potential temperature (K)
 ,qe(npnts,klev+1)       & ! Model enviromental mixing ratio (kg/kg)
 ,qse(npnts,klev+1)      & ! Model enviromental qsat mixing ratio (kg/kg)
 ,trap(npnts,nlev,ntra)  & ! Parcel tracer (kg/kg)
 ,trae(npnts,trlev,ntra) & ! Environment tracer (kg/kg)
 ,flx(npnts,klev+1)        ! updraught mass flux (Pa/s)

REAL, INTENT(INOUT) ::      &
  precip(npnts,klev+1)      & ! precipitation added when descending
                              ! from layer k to k-1 (kg/m**2/s)
 ,dthbydt(npnts,klev+1)     & ! increment to model potential temperature (K/s)
 ,dqbydt(npnts,klev+1)      & ! increment to model mixing ratio (kg/kg/s)
 ,dtrabydt(npnts,nlev,ntra) & ! increment to model tracers (kg/kg/s)
 ,rain(npnts)               & ! rainfall at surface (kg/m**2/s)
 ,snow(npnts)               & ! snowfall at surface (kg/m**2/s)
 ,rain_3d(npnts,klev+1)     & ! rainfall flux  (kg/m**2/s)
 ,snow_3d(npnts,klev+1)       ! snowfall flux  (kg/m**2/s)

REAL, INTENT(OUT) ::        &
  entrain_dwn(npnts,klev+1) & ! fractional entrainment rate for downdraught
                              ! mass flux
 ,detrain_dwn(npnts,klev+1) & ! fractional detrainment rate for downdraught
                              ! mass flux
 ,dd_flux(npnts,klev+1)       ! Downdraught mass flux (Pa/s)

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
  iccb_c(npossdd) & ! Compressed cloud base level
 ,kmin(npossdd)     ! Freezing level where entrainment rates are increased

LOGICAL ::             &
  bwater_k_c(npossdd)  & ! Compressed mask for those points at which condensate
                         ! is water in layer k
 ,bddi_c(npossdd)      & ! Compressed mask for points where downdraught may
                         ! initiate
 ,bdd_start_c(npossdd) & ! Compressed mask for those points where downdraught
                         ! is able to start from level k
 ,bddwt_k_c(npossdd)   & ! Compressed mask for points in DD where ppt in
                         ! layer k is liquid
 ,bddwt_km1_c(npossdd) & ! Compressed mask for points in DD where ppt in
                         ! layer k-1 is liquid
 ,bdd_on_c(npossdd)      ! Compressed mask for points where DD continues from
                         ! layer K+1

REAL ::                 &
  exner_km12_c(npossdd) & ! Compressed exner function at layer k
 ,exner_kp12_c(npossdd) & ! Compressed exner function at layer k+1
 ,exner_KM32_c(npossdd) & ! Compressed exner function at layer k-1

 ,pk(npossdd)           & ! Pressure of layer k (Pa)
 ,p_km1(npossdd)        & ! Pressure of layer k-1 (Pa)
 ,exk(npossdd)          & ! exner ratio for layer K
 ,exkm1(npossdd)        & ! exner ratio for layer K-1

 ,delpk(npossdd)        & ! Pressure difference across layer K  (Pa)
 ,delpkm1(npossdd)      & ! Pressure difference across layer K-1 (Pa)

 ,amdetk(npossdd)       & ! Mixing detrainment at level k multiplied by 
                          ! appropriate layer thickness
 ,ekm14(npossdd)        & ! exner ratio at layer k-1/4
 ,ekm34(npossdd)        & ! exner ratio at layer k-3/4

 ,precip_k_c(npossdd)   & ! Compressed precipitation added when descending 
                          ! from layer K to K-1 (kg/m**2/s)
 ,q_k_c(npossdd)        & ! Compressed parcel mixing ratio of layer K (kg/kg)

 ,th_k_c(npossdd)       & ! Compressed parcel potential temperature of
                          ! layer k (K)
 ,tra_k_c(npossdd,ntra) & ! Compressed parcel tracer in layer K (kg/kg)

 ,pstar_c(npossdd)      & ! Compressed surface pressure (Pa)

 ,recip_pstar_c(npossdd)  ! Reciprocal of comp. pstar array

REAL ::                                    &
  P_layer_centres_c(npossdd,0:klev+2)      & ! Pressure (Pa)
 ,P_layer_boundaries_c(npossdd,0:klev+1)   & ! Pressure (Pa)
 ,exner_layer_centres_c(npossdd,0:klev+1)    ! exner

REAL ::                   &
  dthbydt_K_c(npossdd)    & ! Compressed increment to model potential 
                            ! temperature of layer k (K/s)
 ,dthbydt_km1_c(npossdd)  & ! Compressed increment to model potential
                            ! temperature of layer k-1 (K/s)
 ,dqbydt_k_c(npossdd)     & ! Compressed increment to model mixing ratio
                            ! of layer k (kg/kg/s)
 ,dqbydt_km1_c(npossdd)   & ! Compressed increment to model mixing ratio 
                            ! of layer k-1(kg/kg/s)
 ,dtra_k_c(npossdd,ntra)  & ! Compressed increment to model tracer of
                            ! layer k (kg/kg/s)
 ,dtra_km1_c(npossdd,ntra)  ! Compressed increment to model tracer of 
                            ! layer k-1 (kg/kg/s)

REAL ::                   &
  deltd(npossdd)          & ! Cooling necessary to achieve saturation (K)

 ,delqd(npossdd)          & ! moistening necessary to achieve saturation (K)

 ,deltrad(npossdd,ntra)   & ! Depletion of environment tracer due to 
                            ! downdraught formation (kg/kg)
 ,qdd_k(npossdd)          & ! mixing ratio of downdraught in layer K (kg/kg)

 ,thdd_k(npossdd)         & ! Model potential temperature of downdraught 
                            ! in layer k (K)
 ,tradd_k(npossdd,ntra)     ! Model tracer of downdraught in layer K (kg/kg)

REAL ::                    &
  flx_dd_k(npnts)          & ! Downdraught initial mass flux (Pa/s)
 ,flx_dd_k_c(npossdd)      & ! Compressed downdraught initial mass flux (Pa/s)
 ,qe_k_c(npossdd)          & ! Compressed environment mixing ratio of 
                             ! layer k  (kg/kg)
 ,qe_km1_c(npossdd)        & ! Compressed environment mixing
                             ! ratio of layer k-1 (kg/kg)
 ,qse_k_c(npossdd)         & ! Compressed environment qsat mixing ratio of 
                             ! layer k  (kg/kg)
 ,qse_km1_c(npossdd)       & ! Compressed environment qsat mixing
                             ! ratio of layer k-1 (kg/kg)
 ,the_k_c(npossdd)         & ! Compressed potential temperature
                             ! of environment in layer k (K)
 ,the_km1_c(npossdd)       & ! Compressed potential temperature
                             ! of environment in layer k-1 (K)
 ,trae_k_c(npossdd,ntra)   & ! Compressed tracer of environment in layer k
                             ! (kg/kg)
 ,trae_km1_c(npossdd,ntra) & ! Compressed tracer of environment in layer k-1 
                             ! (kg/kg)
 ,rain_c(npossdd)          & ! Compressed surface rainfall (kg/m**2/s)

 ,snow_c(npossdd)          & ! Compressed surface snowfall (kg/m**2/s)

 ,flx_ud_k_c(npossdd)      & ! updraught mass flux at layer K

 ,rain_env(npossdd)        & ! Amount of rainfall passing through environment 
                             ! (kg/m**2/s)
 ,snow_env(npossdd)        & ! Amount of snowfall passing through environment 
                             ! (kg/m**2/s)
 ,rain_dd(npossdd)         & ! Amount of rainfall passing through
                             ! downdraught (kg/m**2/s)
 ,snow_dd(npossdd)         & ! Amount of snowfall passing through
                             ! downdraught (kg/m**2/s)
 ,flx_strt_c(npossdd)      & ! Compressed value of flx_strt

 ,cca_c(npossdd)           & ! Compressed convective cloud amount

 ,lr_ud_ref(npossdd)         ! precipitation mixing ratio at lowest
                             ! precipitationg level of UD

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0      
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
! Compression to Downdraught points  (all levels even above term level)
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('DD_ALL_CALL_4A5A',zhook_in,zhook_handle)

ndd=npossdd

DO k=0,klev+1
  DO i=1,ndd
    p_layer_centres_c(i,k)     = p_layer_centres(index1(i),k)
    exner_layer_centres_c(i,k) = exner_layer_centres(index1(i),k)
    p_layer_boundaries_c(i,k)  = p_layer_boundaries(index1(i),k)
  END DO
END DO

!----------------------------------------------------------------------
! Initailise logical arrays as .FALSE.
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
! Done on all npnts but uses test on bddi (klev not used).
!-----------------------------------------------------------------------

! DEPENDS ON: flx_init_4a5a
CALL flx_init_4a5a(npnts, klev, iccb, icct, flx, flx_dd_k, bddi, flx_strt)

!-----------------------------------------------------------------------
! Compress all input arrays for the downdraught calculation down to
! those points where a downdraught is possible when the convection
! terminates in the column.
!-----------------------------------------------------------------------
! Initialise various arrays before level loop

DO i=1,ndd
  flx_dd_k_c(i) = flx_dd_k(index1(i))
  flx_strt_c(i) = flx_strt(index1(i))
  pstar_c(i)    = pstar(index1(i))
  recip_pstar_c(i)=recip_pstar(index1(i))
  iccb_c(i) = iccb(index1(i))
  bdd_start_c(i) = bdd_start(index1(i))

  rain_c(i) = rain(index1(i))
  snow_c(i) = snow(index1(i))

  bddwt_k_c(i)   = bddwt_k(index1(i))
  bddwt_km1_c(i) = bddwt_km1(index1(i))
  bdd_on_c(i)    = bdd_on(index1(i))

  cca_c(i) = cca(index1(i))
  lr_ud_ref(i) = 0.0

! Initialise compressed downdraught indicator array to false

  bddi_c(i) = .FALSE.
END DO


!-----------------------------------------------------------------------
! Main level loop working from top Downwards
!-----------------------------------------------------------------------

DO k = klev+1,2,-1

  ! Need at this stage to reset bddi_c to true if reached level where
  ! convection terminated

  DO I=1,ndd
    IF (kterm(index1(i))+1 == k) THEN
      bddi_c(i) = .TRUE.
    END IF

  ! Compress arrays to those points with downdraughts possible in the column
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
        tra_k_c(I,ktra)    = trap(index1(i),k,ktra)
        trae_k_c(I,ktra)   = trae(index1(i),k,ktra)
        trae_km1_c(I,ktra) = trae(index1(i),k-1,ktra)
        dtra_k_c(I,ktra)   = dtrabydt(index1(i),k,ktra)
        dtra_km1_c(I,ktra) = dtrabydt(index1(i),k-1,ktra)
      END DO
    END DO

  END IF


!----------------------------------------------------------------------
! If below convective cloud base downdraught note allowed to form
!----------------------------------------------------------------------

  DO i=1,ndd
    IF (k <  iccb_c(i)) bddi_c(i)=.FALSE.
  END DO

!-----------------------------------------------------------------------
! Reset en/detrainment rates for downdraught
!-----------------------------------------------------------------------
! Possible problem with calculation of kmin in layer_dd when looping
! over levels starting well above the termination level.
! kmin only calculated in layer_dd if k=klev+1

  IF ( k /= klev+1) THEN
    DO i=1,ndd

      IF (kterm(index1(i))+1 == k) THEN

        kmin(i)=klev+2      ! required to ensure amdetk calculated
                            ! correctly in layer_dd
      END IF

    END DO
  END IF

! DEPENDS ON: layer_dd_4a5a
  CALL layer_dd_4a5a (ndd,k,klev,the_k_c,the_km1_c,flx_strt_c,           &
                      P_layer_centres_c,P_layer_boundaries_c,            &
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
! downdraught not allowed to form from cloud top layer (Klev+1)
! or from below cloud base
!-----------------------------------------------------------------------

! DEPENDS ON: dd_init
  CALL dd_init(ndd,npossdd,th_k_c,q_k_c,the_k_c,qe_k_c,qse_k_c,pk,exk,  &
                   thdd_k,qdd_k,deltd,delqd,bdd_start_c,k,bddi_c,       &
                   bdd_on_c,l_tracer,ntra,tra_k_c,                      &
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
  CALL downd_4a5a(ndd,npossdd,k,klev,ntra,nddon_tmp, iccb_c,                 &
                  l_tracer,bwater_k_c,                                       &
                  timestep,the_k_c,the_km1_c,                                &
                  qe_k_c,qe_km1_c,qse_km1_c,p_km1,delpk,delpkm1,exk,         &
                  exkm1,deltd,delqd,amdetk,ekm14,ekm34,flx_ud_k_c,cca_c,     &
                  trae_k_c,trae_km1_c,deltrad,                               &
                  bdd_start_c,bddwt_k_c,bddwt_km1_c,bdd_on_c,                &
                  thdd_k,qdd_k,dthbydt_k_c,dthbydt_km1_c,dqbydt_k_c,         &
                  dqbydt_km1_c,rain_c,snow_c,precip_k_c, rain_env,snow_env,  &
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
        entrain_dwn(index1(i),k)=(1.0-amdetk(i))*                           &
                                 (ekm14(i)+ekm34(i)*(1.0+ekm14(i)))*        &
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

IF (lhook) CALL dr_hook('DD_ALL_CALL_4A5A',zhook_out,zhook_handle)
RETURN

END SUBROUTINE dd_all_call_4a5a
