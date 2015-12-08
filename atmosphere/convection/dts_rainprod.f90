! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! + Calculates the rain production term

SUBROUTINE dts_rainprod(n_dp,nlev,                                          &
                        z_rho,z_theta,dr_across_rh,dr_across_th,            &
                        rho,rho_theta,temperature,q,qse,zlcl,zfr,           &
                        condensation,evaporation,melt,freeze,               &
                        dqclbydt,dqcfbydt,deposition,sublimation,           &
                        rainprod,snowprod,revp,revpint,rainrate,snowrate)

USE dts_fitpars_mod, ONLY:                                                  &
  alfac

USE conversions_mod, ONLY: zerodegc

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

!-----------------------------------------------------------------------
! Description:
!
! Calculates the rain production term
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
!-----------------------------------------------------------------------

! Subroutine Arguments

INTEGER, INTENT(IN) :: &
  n_dp                 & ! number of deep points
 ,nlev                   ! Number of levels   


REAL, INTENT(IN)    ::    &
  temperature(n_dp,nlev)  & ! temperature (K)
 ,q(n_dp,nlev)            & ! water vapour (kg/kg)
 ,qse(n_dp,nlev)          & ! saturation value (kg/kg)
 ,dqclbydt(n_dp,nlev)     & ! rate of change of qcl (kg/kg/s)
 ,dqcfbydt(n_dp,nlev)     & ! rate of change of qcf (kg/kg/s)
 ,condensation(n_dp,nlev) & ! condensation (kg/kg/s)
 ,evaporation(n_dp,nlev)  & ! evaporation  (kg/kg/s)
 ,deposition(n_dp,nlev)   & ! deposition (kg/kg/s)
 ,sublimation(n_dp,nlev)  & ! sublimation (kg/kg/s)
 ,melt(n_dp,nlev)         & ! melting  (kg/kg/s)
 ,freeze(n_dp,nlev)       & ! freezing (kg/kg/s)
 ,z_rho(n_dp,nlev)        & ! height on rho levels
 ,z_theta(n_dp,nlev)      & ! height on theta levels
 ,dr_across_rh(n_dp,nlev) & ! thickness of rho layers (m)
 ,dr_across_th(n_dp,nlev) & ! thickness of theta layers (m)
 ,rho(n_dp,nlev)          & ! density on rho levels (kg/m3)
 ,rho_theta(n_dp,nlev)    & ! density on theta levels (kg/m3)
 ,zfr(n_dp)               & ! height of freezing level (m)
 ,zlcl(n_dp)                ! height of lifting condensation level (m) 

REAL, INTENT(OUT) ::    &
  rainprod(n_dp,nlev)   & ! rain production (kg/kg/s)
 ,snowprod(n_dp,nlev)   & ! snow production (kg/kg/s)
 ,revp(n_dp,nlev)       & ! reevaporation rate (kg/kg/s)
 ,revpint(n_dp)         & ! column integral of reevaporation rate
 ,rainrate(n_dp)        & ! surface rainfall rate 
 ,snowrate(n_dp)          ! surface snowfall rate

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
INTEGER ::         &
  i_dp,k              ! Loop counters

REAL ::              &
  rtot(n_dp)         & ! total prec production
 ,fraintot(n_dp)     &
 ,alpha              &
 ,rh                 & ! relative humidity
 ,frain 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

     
!-----------------------------------------------------------------------

  IF (lhook) CALL dr_hook('DTS_RAINPROD',zhook_in,zhook_handle)

  revpint(:) = 0.0
  rtot(:) = 0.0
  fraintot(:) = 0.0
  revp(:,:) = 0.0

! assume all ice becomes falling hydrometeor (which for the moment is
! presumed to melt...)

! assume falling ice is just what is left over from deposition,
! sublimation etc:
! assume rain production is equal to liquid water production:

  snowprod(:,:) = deposition(:,:)+freeze(:,:)-melt(:,:)-sublimation(:,:) &
                                   -dqcfbydt(:,:)

  rainprod(:,:) = condensation(:,:)-evaporation(:,:)-freeze(:,:)         &
                        +melt(:,:)-dqclbydt(:,:)

    

! ------------------------------------------------------
! Evaporation of rain 
! -------------------
! This is modelled to be proportional to the rain profile multiplied by
! (1-q/qsat). The rain profile is given a fixed functional form
! and is normalised by the rain rate 

! All other levels        
    DO K=1,nlev-1
      DO i_dp=1,n_dp
        IF(z_theta(i_dp,k) <= zfr(i_dp)) THEN
! use generic rainfall shape, which then gets scaled to work out
! evaporation term
          frain = 4.0*(zfr(i_dp)-z_theta(i_dp,k))*                       &
                       (z_theta(i_dp,k)+0.3*zfr(i_dp))/                  &
                       (1.3*zfr(i_dp))**2
              
          fraintot(i_dp) = fraintot(i_dp) + frain                        &
                            *rho_theta(i_dp,k)*dr_across_th(i_dp,k)
                
          alpha = alfac*(1.0-z_theta(i_dp,k)/zfr(i_dp))

! ensure that rain evaporation responds to relative humidity
          rh = q(i_dp,k)/qse(i_dp,k)
          revp(i_dp,k) = alpha*frain*(1.0-rh) 
          IF(revp(i_dp,k) < 0.0) THEN
            revp(i_dp,k) = 0.0
          END IF
         
        END IF ! z <= zfr

        rtot(i_dp) = rtot(i_dp) + rho_theta(i_dp,k)*                        & 
                     (rainprod(i_dp,k)+snowprod(i_dp,k))*dr_across_th(i_dp,k)

      END DO ! i_dp
    END DO ! k
        
        
    DO K=1,nlev-1
! now renormalise by the actual rain rate (without evap)
      DO i_dp=1,n_dp
        IF(rtot(i_dp) > 0.0 .AND. fraintot(i_dp) > 0.0) THEN
          revp(i_dp,k) = revp(i_dp,k)*rtot(i_dp)/fraintot(i_dp)
          rainprod(i_dp,k) = rainprod(i_dp,k)-revp(i_dp,k)                 

! calculate revp integral
          revpint(i_dp) = revpint(i_dp) +  rho_theta(i_dp,k)*  & 
                                    (revp(i_dp,k))*dr_across_th(i_dp,k)
        ELSE
          revp(i_dp,k) = 0.0
          revpint(i_dp) = 0.0
          rainprod(i_dp,k) = 0.0
        END IF
      END DO
    END DO ! k


        
! Calculate the rain rate by integrating the production terms as a
! function of height
    rainrate(:) = 0.0
    snowrate(:) = 0.0
    DO K=1,nlev-1
      DO i_dp=1,n_dp
! nb which levels should it be on?
        snowrate(i_dp) = snowrate(i_dp) + rho_theta(i_dp,k)              &
                                   *snowprod(i_dp,k)*dr_across_th(i_dp,k)
        rainrate(i_dp) = rainrate(i_dp) + rho_theta(i_dp,k)              &
                                   *rainprod(i_dp,k)*dr_across_th(i_dp,k)

      END DO
    END DO

! Deal with negative precipitation case -- evaporated more than was there:
    DO i_dp=1,n_dp
      IF(rainrate(i_dp)+snowrate(i_dp) <= 0.0) THEN

        rainrate(i_dp) = 0.0
        snowrate(i_dp) = 0.0

        ! renormalise evap so that it matches whatever was produced

        IF(revpint(i_dp) > 0.0) THEN
          DO K=1,nlev
            revp(i_dp,k) = revp(i_dp,k)*rtot(i_dp)/revpint(i_dp)
          END DO
        END IF

      END IF
    END DO

!--------------------------------------------------------------
! Simplistic way of converting all to snow if temperature on first     
! model level is less than 0C
! nb there are almost certainly better algorithms for this! 
  DO i_dp=1,n_dp
    IF(temperature(i_dp,1) < ZeroDegC) THEN
      snowrate(i_dp) = snowrate(i_dp) + rainrate(i_dp)
      rainrate(i_dp) = 0.0
    ELSE
      rainrate(i_dp) = rainrate(i_dp) + snowrate(i_dp)
      snowrate(i_dp) = 0.0
    END IF
  END DO
  IF (lhook) CALL dr_hook('DTS_RAINPROD',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE dts_rainprod
