! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Evaporate rain below cloud base if no downdraught
!
! Subroutine Interface:  
!
SUBROUTINE pevp_bcb (npnts,k,iccb,th,pk,q,qse,delp,rain,snow,         &
                     dthbydt,dqbydt,exk,timestep,cca)


USE cv_param_mod, ONLY: cldarea

USE earth_constants_mod, ONLY: g

USE atmos_constants_mod, ONLY: r, cp

USE water_constants_mod, ONLY: lc, lf

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!-----------------------------------------------------------------------
! Description:  Evaporate rain below cloud base if no downdraught.
!
! Method: UM documentation paper 27
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.3.
!-----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------

INTEGER, INTENT(IN) :: &
  npnts                & ! Vector length
 ,k                      ! present model layer 

INTEGER, INTENT(IN) :: &
  iccb(npnts)            ! convective cloud base layer    

REAL, INTENT(IN) ::    &
  th(npnts)            & ! Potential temperature (K)
 ,pk(npnts)            & ! Pressure (Pa)
 ,q(npnts)             & ! Mixing ratio (kg/kg)
 ,qse(npnts)           & ! Environmental qsat Mixing ratio (kg/kg)
 ,delp(npnts)          & ! Change in pressure across layer k-1 (Pa)
 ,exk(npnts)           & ! Exner ratio of layer k
 ,timestep             & ! Convection timestep (s)
 ,cca(npnts)             ! Convective cloud amount

REAL, INTENT(INOUT) ::   &
  rain(npnts)            & ! IN  Amount of falling rain (kg/m**2/s)
                           ! OUT Updated amount of falling rain
 ,snow(npnts)            & ! IN  Amount of falling snow (kg/m**2/s)
                           ! OUT Updated amount of falling snow
 ,dthbydt(npnts)         & ! IN  Increment to model potential temperature
                           ! OUT Updated Increment to model potential
                           !     temperature (K/s)
 ,dqbydt(npnts)            ! IN  Increment to model mixing ratio (kg/kg/s)
                           ! OUT Updated Increment to model mixing ratio
                           !     after evaporation below cloud base.

!-----------------------------------------------------------------------
! Local variables 
!-----------------------------------------------------------------------

INTEGER ::        &
  i                 ! loop counter

LOGICAL ::        &
  bevap(npnts)          & ! Mask for those points where evaporation occurs
 ,l_rates_adjusted      & ! True if evaporation / sublimation rates have
                          ! been adjusted to give saturation 
 ,full_evap_rain(npnts) & ! True if full rain evap in evp 
 ,full_evap_snow(npnts)   ! True if full snow evap in evp 

REAL ::               &
  t(npnts)            & ! Model temperature (K)
 ,evap_rain(npnts)    & ! Amount of evaporation of rain
 ,sub_snow(npnts)     & ! Amount of snow sublimation
 ,qsate(npnts)        & ! Saturation mixing ratio in environment (kg/kg)
 ,delq(npnts)         & ! Change in mixing ratio across layer k  (kg/kg)
 ,ths(npnts)          & ! Saturated parcel potential temperature (K)
 ,qs(npnts)           & ! Saturated parcel parcel mixing ratio (kg/kg)
 ,dthbydt_evp(npnts)  & ! Increment to potential temperature due to 
                        ! evaporation (K)
 ,dqbydt_evp(npnts)   & ! Increment to mixing ratio due to 
                        ! evaporation (kg/kg)
 ,dthbydt_sat(npnts)  & ! Increment to potential temperature due to 
                        ! saturation (K)
 ,factor(npnts)       & ! dthbydt_sat/dthbydt_evp
 ,rho(npnts)            ! Density of air in parcel

REAL ::               &
  cloud_fraction      & ! cca*cldarea
 ,rtimestep             ! 1/timestep 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('PEVP_BCB',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Evaporate rain in layer k if layer k is below cloud base.
! Calculate moisture sub-saturation
!-----------------------------------------------------------------------

DO i=1,npnts
  t(i) = th(i)*exk(i)
  bevap(i) = .FALSE.
END DO


DO i=1,npnts
  IF (k  <   iccb(i)) THEN
    delq(i) = qse(i)-q(i)

    !-----------------------------------------------------------------------
    ! Check if evaporation possible
    !-----------------------------------------------------------------------

    IF ((rain(i) >  0.0 .OR. snow(i) >  0.0) .AND.                 &
        delq(i)  >   0.0) THEN

      bevap(i) = .TRUE.
      rho(i) = pk(i) / (r*t(i))
    END IF
  END IF
END DO

!-----------------------------------------------------------------------
! Calculate evaporation
!-----------------------------------------------------------------------

! DEPENDS ON: evp
CALL evp (npnts,1,bevap,cldarea,rain,t,cca,rho,delq,delp,      &
            pk,evap_rain,full_evap_rain)

! DEPENDS ON: evp
CALL evp (npnts,2,bevap,cldarea,snow,t,cca,rho,delq,delp,      &
            pk,sub_snow,full_evap_snow)

!-----------------------------------------------------------------------
! Calculate temperature and mixing ratio if layer brought to saturation
! by evaporation and sublimation
!-----------------------------------------------------------------------

! DEPENDS ON: satcal
CALL satcal(npnts,k,t,th,pk,exk,q,th,qse,qs,ths)

rtimestep = 1.0/timestep

DO i=1,npnts
  IF (bevap(i)) THEN
    dthbydt_evp(i) = -((lc*evap_rain(i))+((lc+lf)*sub_snow(i)))*g/       &
                       (cp*exk(i)*delp(i))
    dqbydt_evp(i)  = (evap_rain(i)+sub_snow(i))*g/delp(i)

    dthbydt_sat(i) = (ths(i)-th(i))*rtimestep

    l_rates_adjusted = .FALSE. 

    IF (dthbydt_evp(i) <  dthbydt_sat(i)) THEN

    !---------------------------------------------------------------------
    ! Adjust evaporation and sublimation rates to give saturation
    !---------------------------------------------------------------------

      l_rates_adjusted = .TRUE.
      factor(i) = dthbydt_sat(i)/dthbydt_evp(i)
      dthbydt_evp(i) = dthbydt_sat(i)
      dqbydt_evp(i) = dqbydt_evp(i)*factor(i)
      evap_rain(i) = evap_rain(i)*factor(i)
      sub_snow(i) = sub_snow(i)*factor(i)
    END IF

    !---------------------------------------------------------------------
    ! Update increments and rainfall and adjust back to gridbox means
    !---------------------------------------------------------------------

    cloud_fraction = cca(i)*cldarea
    dthbydt(i) = dthbydt(i)+dthbydt_evp(i)*cloud_fraction
    dqbydt(i)  = dqbydt(i) +dqbydt_evp(i)*cloud_fraction
    rain(i) = rain(i)-evap_rain(i)*cloud_fraction
    snow(i) = snow(i)-sub_snow(i)*cloud_fraction

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

  END IF
END DO

IF (lhook) CALL dr_hook('PEVP_BCB',zhook_out,zhook_handle)
RETURN
END SUBROUTINE pevp_bcb
