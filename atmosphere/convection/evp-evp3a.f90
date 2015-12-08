! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates the evaporation of precipitation
!
! Subroutine Interface:  
!
SUBROUTINE evp(npnts, iphase, bevap, area_fac                              &
               , precip, tevp, cca, rho, delq, delpkm1, pkm1               &
               , evap, full_evap)

USE earth_constants_mod, ONLY: g

USE cv_param_mod, ONLY:                                                    &
   p_lq1, p_lq2, p_ice1, p_ice2, rho_lqp1, rho_lqp2, rho_lqa, rho_lqb,     &
   rho_icp1, rho_icp2, rho_icea, rho_iceb,                                 &
   lq_a, lq_b, lq_c, ice_a, ice_b, ice_c


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
! 
! Description: Calculates the evaporation of precipitation
!
! Method: UM documentation paper 27
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
INTEGER, INTENT(IN) ::  &
  npnts                 & ! Vector length
 ,iphase                  ! Indication for rain (1), or snow (2)

LOGICAL, INTENT(IN) ::  &
  bevap(npnts)            ! Mask for points where evaporation takes place   

REAL, INTENT(IN) ::     &
  area_fac                ! Fraction of convective cloud amount to give 
                          ! local cloud area
REAL, INTENT(IN) ::     &
  precip(npnts)         & ! Amount of precipitation (kg/m**2/s)

 ,tevp(npnts)           & ! Temperature of layer k (K)

 ,cca(npnts)            & ! Convective cloud amount (fraction)  

 ,rho(npnts)            & ! Density of air

 ,delq(npnts)           & ! Change in humidity mixing ratio across layer k 
                          ! (kg/kg)
 ,delpkm1(npnts)        & ! Change in pressure across layer k-1 (Pa)

 ,pkm1(npnts)             ! Pressure at level k-1 (Pa)

REAL, INTENT(OUT) ::    &
  evap(npnts)             ! Evaporation

LOGICAL, INTENT(OUT) :: &
  full_evap(npnts)        ! True if all precip evaporated

!-----------------------------------------------------------------------
! Local variables 
!-----------------------------------------------------------------------

INTEGER ::       &
  i               ! Loop counter


REAL ::        & 
  econ         & ! Quadratic term
 ,c1           & ! Constant
 ,c2           & ! Constant
 ,lrate        & ! Local rate of recipitation
 ,ca           & ! Lcal cloud area
 , tl1,ti1

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
! Start of routine
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('EVP',zhook_in,zhook_handle)

full_evap(:) = .FALSE.

tl1=0.5*p_lq1
ti1=0.5*p_ice1

IF (iphase == 1) THEN        ! Rain

  DO i=1,npnts
    IF (bevap(i)) THEN
      IF (precip(i)  >   0.0) THEN
        econ = ((lq_a*tevp(i)+lq_b)*tevp(i)+lq_c)* (100000.0/pkm1(i))   
        ca = area_fac*cca(i)
        lrate = precip(i)/ca
        c1 = rho_lqa*ca*(lrate*lrate*rho(i))**tl1
        c2 = rho_lqb*ca*lrate**p_lq2*rho(i)**rho_lqp2
        evap(i) = MIN(econ*(c1+c2)*delq(i)*delpkm1(i)/g,lrate)
        IF (evap(i) == lrate) full_evap(i) = .TRUE.
      ELSE
        evap(i) = 0.0
      END IF
    END IF
  END DO

ELSE IF (iphase == 2) THEN        ! Snow

  DO i=1,npnts
    IF (bevap(i)) THEN
      IF (precip(i)  >   0.0) THEN
        IF(tevp(i) <= 243.58) THEN
          econ = 1.7405e-5*(100000.0/pkm1(i))
        ELSE
          econ = ((ice_a*tevp(i)+ice_b)*tevp(i)+ice_c)*(100000.0/pkm1(i))
        END IF
        ca = area_fac*cca(i)
        lrate = precip(i)/ca
        c1 = rho_icea*ca*(lrate*lrate*rho(i))**ti1
        c2 = rho_iceb*ca*lrate**p_ice2*rho(i)**rho_icp2
        evap(i)=MAX(0.,MIN(econ*(c1+c2)*delq(i)*delpkm1(i)/g,lrate))
        IF (evap(i) == lrate) full_evap(i) = .TRUE.
      ELSE
        evap(i) = 0.0
      END IF
    END IF
  END DO

END IF     ! test on iphase

IF (lhook) CALL dr_hook('EVP',zhook_out,zhook_handle)
RETURN
END SUBROUTINE evp

