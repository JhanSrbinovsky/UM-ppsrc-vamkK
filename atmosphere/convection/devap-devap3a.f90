! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Evaporation routine
!
! Subroutine Interface:
!
SUBROUTINE devap(npnts, bddwt_km1                                      &
                 ,thdd_k, thdds, qdds, flx_dd_km1, exk, exkm1          &
                 ,qsatdd, delpkm1, cca, pkm1                           &
                 ,thdd_km1, qdd_km1, rain, snow)

USE earth_constants_mod, ONLY: g

USE atmos_constants_mod, ONLY: cp, r

USE cv_param_mod, ONLY: ddcldfra

USE water_constants_mod, ONLY: lc, lf
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! 
! Description: Evaporation routine
!              Carries out evaporation and updates precipitation
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
  npnts                  ! Vector length

LOGICAL, INTENT(IN) :: &
  bddwt_km1(npnts)       ! Mask for those points in downdraught where 
                         ! precipitation is liquid 

REAL, INTENT(IN) ::    &
  thdd_k(npnts)        & ! Potential temperature of downdraught in layer k (K)

 ,thdds(npnts)         & ! Potential temperature of saturated downdraught
                         ! in layer (K)
 ,qdds(npnts)          & ! Saturated downdraught mixing ratio of layer 
                         ! (kg/kg)
 ,flx_dd_km1(npnts)    & ! Downdraught mass flux in layer k-1 (Pa/s)

 ,exk(npnts)           & ! Exner ratio in layer k
 ,exkm1(npnts)         & ! Exner ratio in layer k-1
 ,qsatdd(npnts)        & ! Saturated downdraught mixing ratio (kg/kg)
 ,delpkm1(npnts)       & ! Change in pressure across layer k-1  (Pa)
 ,cca(npnts)           & ! Convective cloud amount (fraction)  
 ,pkm1(npnts)            ! Pressure in layer k-1  (Pa)

REAL, INTENT(INOUT) :: &
  thdd_km1(npnts)      & ! In  Potential temperature of DD in layer k-1(kg/kg)
                         ! Out Potential temperature of DD in layer k-1 after
                         ! evaporation of saturation (K)  
 ,qdd_km1(npnts)       & ! In  Mixing ratio of downdraught in layer k-1 (kg/kg) 
                         ! Out Mixing ratio of DD in layer k-1 after 
                         ! evaporation of saturation (kg/kg)
 ,rain(npnts)          & ! In  Amount of rain
                         ! Out Updated amount of rainfall (kg/m**2/s)
 ,snow(npnts)            ! In  Amount of snow
                         ! Out Updated amount of snowfall (kg/m**2/s)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! Local variables 
!-----------------------------------------------------------------------

INTEGER :: i               ! Loop counter

LOGICAL ::              &
  bevap(npnts)          & ! Mask for those points at which evaporation
                          ! calculation is to be carried out
 ,bsat(npnts)           & ! Mask for those points which are subsaturated
 ,l_rates_adjusted      & ! True if evaporation / sublimation rates have
                          ! been adjusted to give saturation
 ,full_evap_rain(npnts) & ! True if full rain evap in evp
 ,full_evap_snow(npnts)   ! True if full snow evap in evp

REAL ::            &
  tevp(npnts)      & ! Temperature used in evaporation calculation (K)

 ,evap_rain(npnts) & ! Amount of evaporation of rain

 ,sub_snow(npnts)  & ! Amount of snow sublimation

 ,delq(npnts)      & ! Difference in mixing ratio (kg/kg)

 ,delth(npnts)     & ! Increment to downdraught potential temperature in layer
                     ! k-1  due to evaporation (K)
 ,delqe(npnts)     & ! Increment to downdraught mixing ratio in layer k-1
                     !  due to evaporation (kg/kg)
 ,delths(npnts)    & ! Saturated potential temperature minus potential
                     ! temperature of downdraught
 ,factor(npnts)    & ! delths/delth

 ,pincr(npnts)     & ! Increase in precipitation if parcel supersaturates

 ,rho(npnts)         ! Density of air in parcel

!-----------------------------------------------------------------------
! Check if evaporation possible
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('DEVAP',zhook_in,zhook_handle)

DO i=1,npnts
  delq(i) = qsatdd(i)-qdd_km1(i)

  bevap(i) =((rain(i) >  0.0) .OR. (snow(i) >  0.0)) .AND. (delq(i) >  0.0)
  bsat(i) = delq(i)  <   0.0

  !-----------------------------------------------------------------------
  ! Calculate temperature used in calculation of evaporation constants
  ! based on temeprature of parcel after unsaturated descent
  !-----------------------------------------------------------------------

  IF (bevap(i)) THEN
    tevp(i) = ((thdd_k(i)*exk(i))+(thdd_km1(i)*exkm1(i)))*0.5
    rho(i) = pkm1(i) / (r*tevp(i))
  END IF
END DO

!-----------------------------------------------------------------------
! Evaporation calculation - calculate rates for rain and snow
!-----------------------------------------------------------------------

! DEPENDS ON: evp
CALL evp(npnts,1,bevap,ddcldfra,rain,tevp,cca,rho,delq,delpkm1,         &
         pkm1,evap_rain,full_evap_rain)

! DEPENDS ON: evp
CALL evp(npnts,2,bevap,ddcldfra,snow,tevp,cca,rho,delq,delpkm1,         &
         pkm1,sub_snow,full_evap_snow)

DO i=1,npnts
  IF (bevap(i)) THEN

  !-----------------------------------------------------------------------
  ! Adjust evaporation and sublimation rates back to grid box means
  !-----------------------------------------------------------------------

    evap_rain(i) = evap_rain(i) * cca(i) * ddcldfra
    sub_snow(i) = sub_snow(i) * cca(i) * ddcldfra

  !-----------------------------------------------------------------------
  ! Check if parcel supersaturated
  !-----------------------------------------------------------------------

    delth(i) = -((lc*evap_rain(i))+((lc+lf)*sub_snow(i)))*g/        &
                                        (cp*exkm1(i)*flx_dd_km1(i))
    delqe(i) = (evap_rain(i)+sub_snow(i))*g/flx_dd_km1(i)

    delths(i) = thdds(i)-thdd_km1(i)

    l_rates_adjusted = .FALSE.

    IF (delth(i) <  delths(i)) THEN

    !-----------------------------------------------------------------------
    ! Adjust evaporation and sublimation rates to give saturation
    !-----------------------------------------------------------------------

      l_rates_adjusted = .TRUE.
      factor(i) = delths(i)/delth(i)
      delth(i)  = delths(i)
      delqe(i)  = delqe(i)*factor(i)
      evap_rain(i) = evap_rain(i)*factor(i)
      sub_snow(i)  = sub_snow(i)*factor(i)
    END IF

    !-----------------------------------------------------------------------
    ! Update T, q and precipitation
    !-----------------------------------------------------------------------

    rain(i) = rain(i)-evap_rain(i)
    IF (rain(i) <  0.0) rain(i)=0.0
    snow(i) = snow(i)-sub_snow(i)
    IF (snow(i) <  0.0) snow(i)=0.0
    thdd_km1(i) = thdd_km1(i)+ delth(i)
    qdd_km1(i)  = qdd_km1(i) + delqe(i)

! If the call to evp has produced an evaporation rate that evaporates
! all the precipitation, and there was no adjustment of the rate to
! account for saturation, logic dictates that the precipitation rate
! should now be exactly zero. However, because of rounding effects it is
! possible to end up with slightly non-zero rates. To get around this,
! we directly reset to zero:
    IF (.NOT.l_rates_adjusted) THEN
      IF (full_evap_rain(i)) rain(i) = 0.0
      IF (full_evap_snow(i)) snow(i) = 0.0
    END IF

    !-----------------------------------------------------------------------
    ! Parcel is supersaturated before evaporation occurs
    ! Bring parcel to saturation and precipitate water
    !-----------------------------------------------------------------------

  ELSE IF (bsat(i)) THEN
    pincr(i)    = (qdd_km1(i)-qdds(i))*flx_dd_km1(i)/g
    qdd_km1(i)  = qdds(i)
    thdd_km1(i) = thdds(i)
    IF (bddwt_km1(i)) THEN
      rain(i) = rain(i)+pincr(i)
    ELSE
      snow(i) = snow(i)+pincr(i)
    END IF
  END IF
END DO

IF (lhook) CALL dr_hook('DEVAP',zhook_out,zhook_handle)
RETURN
END SUBROUTINE devap

