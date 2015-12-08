! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Downdraught routine
!
! Subroutine Interface: (Argument list order does not yet obey coding standards)
!
SUBROUTINE ddraught_6a (npnts, np_full, k, kct, ntra,                        &
                     l_tracer,                                               &
                     thdd_k,qdd_k,the_k,                                     &
                     the_km1,qe_k,qe_km1,qse_km1,dthbydt_k,dthbydt_km1,      &
                     dqbydt_k,dqbydt_km1,flx_dd_k,p_km1,delpk,               &
                     delpkm1,exk,exkm1,deltd,delqd,amdetk,ekm14,             &
                     ekm34,rain,snow,                                        &
                     bdd_start,bddwt_k,bddwt_km1,                            &
                     bdd_on,b_dd_end,                                        &
                     deltrad, cca, ppn_mix_dd,                               &
                     tradd_k,trae_k,trae_km1,dtrabydt_k,                     &
                     dtrabydt_km1)

USE water_constants_mod, ONLY: tm
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! 
! Description: Downdraught routine
!              Convective downdraught based on parcel theory.
!              Carry out dry descent.
!              Calculate subsaturation.
!              Calculate effect on the environment.
!
! Method: UM documentataion paper 27
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
  npnts                &  ! Vector length
 ,np_full              &  ! Full vector length
 ,k                    &  ! Present model layer
 ,kct                  &  ! Convective cloud top
 ,ntra                    ! Number of tracers

LOGICAL, INTENT(IN) :: &
  l_tracer                ! Switch for tracers


REAL, INTENT(IN) ::     &
  the_km1(npnts)        & ! Potential temperature of enviroment in layer k-1(K)
 ,qe_km1(npnts)         & ! Mixing ratio of enviroment in layer k-1 (kg/kg)
 ,qse_km1(npnts)        & ! qsat Mixing ratio of enviroment in layer k-1 (kg/kg)
 ,p_km1(npnts)          & ! Pressure in layer k-1  (Pa)
 ,delpk(npnts)          & ! Change in pressure across layer k (Pa)
 ,delpkm1(npnts)        & ! Change in pressure across layer k-1  (Pa)
 ,exk(npnts)            & ! Exner ratio in layer k
 ,exkm1(npnts)          & ! Exner ratio in layer k-1
 ,amdetk(npnts)         & ! Mixing detrainment rate
 ,ekm14(npnts)          & ! Exner ratio at layer k-1/4
 ,ekm34(npnts)          & ! Exner ratio at layer k-3/4
 ,deltd(npnts)          & ! Cooling necessary to achieve saturation (K)
 ,delqd(npnts)          & ! moistening necessary to achieve saturation (kg/kg)
 ,cca(npnts)            & ! Convective cloud amount (fraction)  
 ,ppn_mix_dd(npnts)       ! precipitation mixing ratio (kg/kg)      

REAL, INTENT(IN) ::       &
  trae_km1(np_full,ntra)  & ! Tracer content of enviroment in layer k-1 (kg/kg) 
 ,deltrad(npnts,ntra)       ! Depletion of environment tracer due to 
                            ! downdraught formation (kg/kg)

LOGICAL, INTENT(INOUT) ::  &
  bdd_on(npnts)            & ! In  Mask for those points where DD has continued
                             !     from layer k+1 
                             ! Out Mask for those points where DD continues
                             !     to layer k-1 
 ,bddwt_k(npnts)           & ! In Mask for those points in downdraught where 
                             ! precipitation is liquid in layer k 
 ,bddwt_km1(npnts)           ! Mask for those points in downdraught where 
                             ! precipitation is liquid in layer k-1

REAL, INTENT(INOUT) :: &
  thdd_k(npnts)        & ! In  Potential temperature of downdraught in layer k 
                         ! Out Potential temperature reset for next layer (K)  
 ,qdd_k(npnts)         & ! In  Mixing ratio of downdraught in layer k   
                         ! Out Mixing ratio reset for next layer (kg/kg)
 ,tradd_k(np_full,ntra)& ! In  Downdraught tracer content of layer k
                         ! Out tracer content reset for next layer (kg/kg)
 ,the_k(npnts)         & ! In  Potential temperature of environment in layer k
                         ! Out Potential temperature reset for next layer (K)  
 ,qe_k(npnts)          & ! In  Mixing ratio of environment in layer k   
                         ! Out Mixing ratio reset for next layer (kg/kg)
 ,trae_k(np_full,ntra) & ! In  environment tracer content of layer k
                         ! Out tracer content reset for next layer (kg/kg)
 ,flx_dd_k(npnts)      & ! In  Downdraught mass flux of layer k
                         ! Out Downdraught mass flux reset for next layer (Pa/s)
 ,rain(npnts)          & ! In  Amount of rain
                         ! Out Updated amount of rainfall (kg/m**2/s)
 ,snow(npnts)          & ! In  Amount of snow
                         ! Out Updated amount of snowfall (kg/m**2/s)
 ,dthbydt_k(npnts)     & ! In  Increment to potential temperature of layer k
                         ! Out Updated increment potential temperature layer k
                         !           (K/s)
 ,dthbydt_km1(npnts)   & ! In  Increment to potential temperature of layer k-1
                         ! Out Updated increment potential temperature layer k-1
                         !           (K/s)
 ,dqbydt_k(npnts)      & ! In  Increment to mixing ratio of layer k
                         ! Out Updated increment mixing ratio layer k (kg/kg)
 ,dqbydt_km1(npnts)      ! In  Increment to mixing ratio  of layer k-1
                         ! Out Updated increment mixing ratio layer k-1 (kg/kg)

REAL, INTENT(INOUT) ::       &
  dtrabydt_k(np_full,ntra)   & ! In  Increment to tracer content of layer k
                               ! Out Updated increment tracer content layer k 
                               ! (kg/kg)
 ,dtrabydt_km1(np_full,ntra)   ! In  Increment to tracer content of layer k-1
                               ! Out Updated increment tracer content layer k-1
                               ! (kg/kg)

LOGICAL, INTENT(OUT) ::  &
  bdd_start(npnts)       & ! Mask for those points where DD may start 
                           ! in layer k-1
 ,b_dd_end(npnts)          ! Mask for those points where DD is ending in
                           ! layer k-1


!-----------------------------------------------------------------------
! Local variables 
!-----------------------------------------------------------------------

INTEGER ::       &
  i, ktra          ! Loop counters

REAL ::               &
  thdd_km1(npnts)     & ! Potential temperature of downdraught in layer k-1 (K)

 ,qdd_km1(npnts)      & ! Downdraught mixing ratio of layer k-1 (kg/kg)

 ,qsatdd(npnts)       & ! Saturated downdraught mixing ratio of layer k-1
                        ! (kg/kg)
 ,tdd_km1(npnts)      & ! Temperature of downdraught in layer k-1 (K)

 ,thdds(npnts)        & ! Potential temperature of saturated downdraught
                        ! in layer k-1 (K)
 ,qdds(npnts)         & ! Saturated downdraught mixing ratio of layer k-1
                        ! (kg/kg)
 ,flx_dd_km1(npnts)   & ! Downdraught mass flux in layer K-1 (Pa/s)

 ,rain_tmp(npnts)     & ! Liquid precipitation store

 ,snow_tmp(npnts)       ! Snow Store

REAL ::                 &
  tradd_km1(npnts,ntra)   ! Tracer content of downdraught in layer k-1 (kg/kg)

REAL ::                 &
  ekm14_plus1           & ! 1+ekm14
 ,ekm34_plus1           & ! 1+ekm34
 ,dnom                  & ! (1+ekm14)*(1+ekm34)
 ,rdnom                   ! 1/dnom

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Calculate mask for those points in downdraught where precipitation
! is liquid
!
! Store precipitation in layer k in temporary variables.
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('DDRAUGHT_6A',zhook_in,zhook_handle)

DO i=1,npnts
  IF (k  ==  kct .OR. bdd_start(i)) THEN
    bddwt_k(i) = thdd_k(i)*exk(i)  >   TM
  ELSE
    bddwt_k(i) = bddwt_km1(i)
  END IF

  rain_tmp(i) = rain(i)
  snow_tmp(i) = snow(i)

!-----------------------------------------------------------------------
! Dry descent from layer k to k-1
!
! Entrainment calculation
!-----------------------------------------------------------------------

  ekm14_plus1 = 1.0+ekm14(i)
  ekm34_plus1 = 1.0+ekm34(i)
  dnom        = ekm14_plus1 *ekm34_plus1
  rdnom        = 1.0/dnom

  thdd_km1(i) = (thdd_k(i)+(ekm14(i)*the_k(i)) +                           &
                                   ekm14_plus1*ekm34(i)*the_km1(i))*rdnom

  qdd_km1(i) = (qdd_k(i)+(ekm14(i)*qe_k(i)) +                              &
                                   ekm14_plus1*ekm34(i)*qe_km1(i))*rdnom

  !-----------------------------------------------------------------------
  ! Update mass flux and calculate temperature of layer k-1
  !-----------------------------------------------------------------------

  flx_dd_km1(i) = flx_dd_k(i)*dnom*(1.0-amdetk(i))

  tdd_km1(i)    = thdd_km1(i)*exkm1(i)

END DO

!----------------------------------------------------------------------
! Dry descent for tracers
!----------------------------------------------------------------------

IF (l_tracer) THEN

  DO ktra=1,ntra
    DO i=1,npnts
      ekm14_plus1 = 1.0+ekm14(i)
      dnom        = ekm14_plus1 * (1.0+ekm34(i))

      tradd_km1(i,ktra)=(tradd_k(i,ktra)+(ekm14(i)*trae_k(i,ktra))         &
                             + ekm14_plus1*ekm34(i)*trae_km1(i,ktra))/dnom


    END DO
  END DO

END IF

!-----------------------------------------------------------------------
! Calculate subsaturation
! Calculate temperature if brought to saturation
!-----------------------------------------------------------------------

! DEPENDS ON: satcal
CALL satcal(npnts,k,tdd_km1,thdd_km1,p_km1,exkm1,qdd_km1,the_km1,qse_km1,   &
            qdds,thdds) 

DO i=1,npnts
  bddwt_km1(i) = thdds(i)*exkm1(i)  >   TM
END DO

!-----------------------------------------------------------------------
! Calculate change of phase due to downdraught saturation temperature
!-----------------------------------------------------------------------
!
! DEPENDS ON: crs_frzl
CALL crs_frzl(npnts, bddwt_km1, exkm1, flx_dd_km1, thdd_km1, rain, snow)

DO i=1,npnts
  tdd_km1(i) = thdd_km1(i)*exkm1(i)
END DO

!-----------------------------------------------------------------------
! Recalculate subsaturation temperature
!-----------------------------------------------------------------------

! DEPENDS ON: satcal
CALL satcal(npnts,k,tdd_km1,thdd_km1,p_km1,exkm1,qdd_km1,the_km1,qse_km1, &
                 qdds, thdds)

!-----------------------------------------------------------------------
! Calculate moisture subsaturation
!-----------------------------------------------------------------------

! DEPENDS ON: qsat
CALL qsat(qsatdd,tdd_km1,p_km1,npnts)

!-----------------------------------------------------------------------
! Evaporation calculation and adujstment of downdraught temperature
! and moisture
!-----------------------------------------------------------------------

! DEPENDS ON: devap
CALL devap(npnts, bddwt_km1                                            &
           , thdd_k, thdds, qdds, flx_dd_km1, exk, exkm1               &
           , qsatdd, delpkm1, cca, p_km1                               &
           , thdd_km1, qdd_km1, rain, snow )

!-----------------------------------------------------------------------
! Check if parcel still negatively buoyant such that downdraught can
! continue to k-1
!-----------------------------------------------------------------------

! DEPENDS ON: termdd
CALL termdd(npnts, k, bdd_start                                        &
           , thdd_km1, qdd_km1, the_km1, qe_km1, ppn_mix_dd            &
           , b_dd_end, bdd_on)

!-----------------------------------------------------------------------
! Calculate the effect on the environment in layer k
!-----------------------------------------------------------------------

! DEPENDS ON: dd_env
CALL dd_env(npnts, np_full, ntra                                           &
                  ,l_tracer, b_dd_end, bdd_start, bdd_on                   &
                  ,thdd_k, thdd_km1, qdd_k, qdd_km1, the_k, the_km1        &
                  ,qe_k, qe_km1, flx_dd_k, flx_dd_km1, delpk, delpkm1      &
                  ,deltd, delqd, amdetk, ekm14                             &
                  ,tradd_k, tradd_km1, trae_k, trae_km1, deltrad           &
                  ,dthbydt_k, dthbydt_km1, dqbydt_k, dqbydt_km1            &
                  ,dtrabydt_k, dtrabydt_km1 )

!-----------------------------------------------------------------------
! Reset downdraught bit vectors
!-----------------------------------------------------------------------

DO i=1,npnts
  bdd_start(i) = .FALSE.
  IF (.NOT. bdd_on(i)) THEN
    rain(i) = rain_tmp(i)
    snow(i) = snow_tmp(i)
  END IF
  IF (b_dd_end(i)) bdd_on(i) = .FALSE.
END DO

!-----------------------------------------------------------------------
! Switch potential temperature, mixing ratio,  mass flux and 
! tracer ready for calculation at next model layer
!-----------------------------------------------------------------------

IF (k >  2) THEN
  DO i=1,npnts
    IF (bdd_on(i)) THEN
      thdd_k(i)   = thdd_km1(i)
      qdd_k(i)    = qdd_km1(i)
      flx_dd_k(i) = flx_dd_km1(i)
    END IF
  END DO

  IF(l_tracer)THEN

    DO ktra=1,ntra
      DO i=1,npnts
        IF(bdd_on(i))THEN
          tradd_k(i,ktra) = tradd_km1(i,ktra)
        END IF
      END DO
    END DO

  END IF

END IF  ! K > 2

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('DDRAUGHT_6A',zhook_out,zhook_handle)
RETURN

END SUBROUTINE ddraught_6a
