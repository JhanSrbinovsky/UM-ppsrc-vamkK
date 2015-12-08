! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ calculates the evaporation rate
!

SUBROUTINE dts_evaporation(nlev,n_dp,dts_ntpar,ntparmax                    &
                          ,temperature,rho_theta,z_theta,dr_across_rh      &
                          ,dr_across_th,qsat_moist_ad,q,mb,wstar,wcld      &
                          ,ztop,zfr,zlcl,melt,freeze,condensation          &
                          ,qclint,condint,freezint,meltint                 &
                          ,evaporation,evapint)

USE atmos_constants_mod, ONLY: cp, r

USE dts_cntl_mod, Only:                                                      &
       dts_gamma_evp, idts_gamma_evp

USE water_constants_mod, ONLY: lc, tm

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!-----------------------------------------------------------------------
! Description:
! calculates the evaporation rate
! calculates gamma (wcld/zcld)*(qsat(env)-qenv)/(1+lv/cp dqsat/dT)
! then renormalises so that it cancels the cloud production term
! ie int(evap) dz = int(condensation-freeze-racw)dz
! note there is no melting term because this is all assumed to go into
! rainfall, which has its own evaporation term 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
!-----------------------------------------------------------------------

! Subroutine Arguments

INTEGER, INTENT(IN) ::  &
  nlev                  & ! No. of model layers
 ,n_dp                  & ! No. convecting points
 ,dts_ntpar(n_dp)       & ! Top level of initial parcel ascent
 ,ntparmax                ! Max ntpar across all conv points

REAL, INTENT(IN)    ::          &
  temperature(n_dp,nlev)        & ! temperature in K
 ,rho_theta(n_dp,nlev)          & ! density on theta levels (kg/m3)
 ,z_theta(n_dp,nlev)            & ! height of theta levels/m
 ,dr_across_rh(n_dp,nlev)       & ! thickness of rho layers (m)
 ,dr_across_th(n_dp,nlev)       & ! thickness of theta layers (m)
 ,qsat_moist_ad(n_dp,nlev)      & ! environmental qsat 
 ,q(n_dp,nlev)                  & ! water vapour  (kg/kg)
 ,mb(n_dp)                      & ! cloud base mass flux
 ,wstar(n_dp)                   & ! Sub cloud velocity scale (m/s)
 ,wcld(n_dp)                    & ! vel scale for conv layer (m/s)
 ,ztop(n_dp)                    & ! z_theta(dts_ntpar)  (m)
 ,zfr(n_dp)                     & ! height of freezing level (m)
 ,zlcl(n_dp)                    & ! Lifting condensation level (m)
 ,melt(n_dp,nlev)               & ! Melting rate (kg/kg/s)
 ,freeze(n_dp,nlev)             & ! Freezing rate (kg/kg/s)
 ,condensation(n_dp,nlev)       & ! condensation rate (kg/kg/s)
 ,meltint(n_dp)                 & ! Column integral of melting
 ,condint(n_dp)                 & ! Column integral of condensation
 ,freezint(n_dp)                & ! Column integral of freezing
 ,qclint(n_dp)                    ! Column integral of qcl

REAL, INTENT(OUT) ::       &
  evaporation(n_dp,nlev)   &      ! evaporation rate (kg/kg/s)
 ,evapint(n_dp)                   ! Column integral of evaporation

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
INTEGER ::         &
  i_dp,k             ! Loop counters

real ::                   &
  dqsedt(n_dp,nlev)       & ! dqsatenv/dT
 ,gamma_evp               & ! 
 ,moistint(n_dp)          &
 ,frac                    &
 ,zsc

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

  IF (lhook) CALL dr_hook('DTS_EVAPORATION',zhook_in,zhook_handle)
  dqsedt(:,:) = 0.0
  evaporation(:,:) = 0.0
        
  evapint(:) = 0.0
  moistint(:) = 0.0

! Or using the Clausius clapeyron relation
   DO k=1,ntparmax
     dqsedt(:,k) = qsat_moist_ad(:,k)*lc/r/temperature(:,k)**2
   END DO ! k

   DO k=1,nlev-1
     DO i_dp=1,n_dp
! Require temperature to be above -40 C
!Finding a tendency to dry out the upper trop, so for the moment
! cheating by allowing evaporation to occur througout the domain

       IF (idts_gamma_evp == 0) THEN ! original code 

         gamma_evp = mb(i_dp)/wstar(i_dp)

       ELSE        ! try fixed values

         gamma_evp = dts_gamma_evp

       END IF
       IF(z_theta(i_dp,k) < ztop(i_dp) .AND.                             &
                   z_theta(i_dp,k) >= zlcl(i_dp)) THEN 

! For the level just above zlcl, only use the fraction of the layer
! above it. Apply a smooth reduction to evap between tm and -40C
! to reflect the expected decrease in liquid water available to
! evaporate between these points (expect it to have frozen)
! give a decaying profile to the evaporation term

         zsc = 1.0 !1000.0/z_theta(i_dp,k) 

         ! decrease linearly to zero up to -40C
         IF(temperature(i_dp,k) < tm .AND. temperature(i_dp,k) > tm-45.0) THEN      
! this is dangerous -- could interact with itself...             
           zsc = zsc*exp(-0.5*(3.0*(tm-temperature(i_dp,k))/40.)**2)
                   
         END IF
         IF(k > 1) THEN 
           IF(z_theta(i_dp,k-1) < zlcl(i_dp)) THEN
             frac = (z_theta(i_dp,k)-zlcl(i_dp))/dr_across_rh(i_dp,k)
           ELSE
             frac = 1.0
           END IF
         ELSE ! k > 1
           frac = 0.0
         END IF

! The qcl_plume term may prove problematic, and may not be correct 30/3/09

         IF(ztop(i_dp) > zlcl(i_dp)) THEN 
           evaporation(i_dp,k) = zsc*frac*gamma_evp*(wcld(i_dp)/ztop(i_dp))  &
                                 *((qsat_moist_ad(i_dp,k)-q(i_dp,k)) /       &
                                   (1.0+(lc/cp)*dqsedt(i_dp,k))  )         

         ELSE
           evaporation(i_dp,k) = zsc*frac*gamma_evp*(wcld(i_dp)/(10000.))    &
                                 *(qsat_moist_ad(i_dp,k)-q(i_dp,k)) /        &
                                   (1.0+(lc/cp)*dqsedt(i_dp,k))
         END IF
                 
         IF(evaporation(i_dp,k) < 0.0) THEN
           evaporation(i_dp,k) = 0.0
         END IF
   ! no evap above -40C level (putting 45 to allow for a much warmer plume)
           IF(temperature(i_dp,k) < tm-45.0) THEN
             evaporation(i_dp,k) = 0.0
           END IF
             evapint(i_dp) = evapint(i_dp) +                                 &
                                evaporation(i_dp,k)*rho_theta(i_dp,k)*       &
                                                    dr_across_th(i_dp,k)
                
         END IF ! T > -40, z >= lcl

       END DO ! i_dp
     END DO ! k
     moistint(:) = condint(:)+meltint(:)-freezint(:)-qclint(:)
        
     DO i_dp=1,n_dp
       IF(evapint(i_dp) > moistint(i_dp)) THEN 
         IF(moistint(i_dp) > 0.0) THEN 
! PROBLEM IN LOOP ORDER
           DO k=1,nlev
             evaporation(i_dp,k) = evaporation(i_dp,k)                       &
                                    *moistint(i_dp)/evapint(i_dp)
             evapint(i_dp) = moistint(i_dp)
           END DO
         ELSE
           evaporation(i_dp,:) = 0.0
           evapint(i_dp) = 0.0
         END IF
       END IF ! evapint > moistint
     END DO
     IF (lhook) CALL dr_hook('DTS_EVAPORATION',zhook_out,zhook_handle)
     RETURN


END SUBROUTINE dts_evaporation
