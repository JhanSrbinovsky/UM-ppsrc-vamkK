! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Change of phase routine where precip crosses a melting or freezing level
!
! Subroutine Interface:  
!
SUBROUTINE crs_frzl (npnts, bddwt_km1, exkm1, flx_dd_km1                   &
                     , thdd_km1, rain, snow)
USE water_constants_mod, ONLY: lf, tm

USE atmos_constants_mod, ONLY: cp
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE earth_constants_mod, ONLY: g

IMPLICIT NONE

! 
! Description: Change of phase routine where precipitation crosses a melting
!              or freezing level
!
! Method: UM documentation paper 27
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards 8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------
INTEGER, INTENT(IN) :: &
  npnts                  ! Vector length

LOGICAL, INTENT(IN) :: &
  bddwt_km1(npnts)       ! Mask for those points in downdraught where 
                         ! precipitation is liquid in layer k-1 

REAL, INTENT(IN) ::    &
  exkm1(npnts)         & ! Exner ratio in layer k-1
 ,flx_dd_km1(npnts)      ! Downdraught mass flux in layer k-1 (Pa/s)

REAL, INTENT(INOUT) :: &
  thdd_km1(npnts)      & ! In  Potential temperature of DD in layer k-1(K)
                         ! Out Potential temperature of DD in layer k-1 
                         ! updated due to change of phase
 ,rain(npnts)          & ! In  Amount of rain descending from k-1 to k-2
                         ! Out Updated amount of rainfall (kg/m**2/s)
 ,snow(npnts)            ! In  Amount of snow descending from k-1 to k-2
                         ! Out Updated amount of snowfall (kg/m**2/s)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
! Model constants
!----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Local variables 
!-----------------------------------------------------------------------

INTEGER :: i         ! Loop counter

REAL ::            &
  factor           & ! Used in the calculation of the change of phase of 
                     ! falling precipitation
 ,precip_fre       & ! Freezing precipitation 

 ,precip_melt        ! Melting precipitation 

!-----------------------------------------------------------------------
! Add latent heating where precip crosses a melting or freezing level
!
!   UM Documentation paper 27
!   Ssection (11), equation (42)
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('CRS_FRZL',zhook_in,zhook_handle)

DO i=1,npnts

  IF (.NOT.bddwt_km1(i).AND.rain(i) >  0.0.AND.thdd_km1(i)        &
       *exkm1(i) <  tm) THEN
    ! Freeze
    factor = (exkm1(i)*cp*flx_dd_km1(i))/(lf*g)
    precip_fre = (tm/exkm1(i)-thdd_km1(i))* factor
    precip_fre = MIN(rain(i),precip_fre)
    thdd_km1(i) = thdd_km1(i)+precip_fre/factor
    rain(i) = rain(i)-precip_fre
    snow(i) = snow(i)+precip_fre

  ELSE IF (bddwt_km1(i).AND.snow(i) >  0.0) THEN
    ! Melt
    factor = (exkm1(i)*cp*flx_dd_km1(i))/(lf*g)
    precip_melt = (thdd_km1(i)-tm/exkm1(i))*factor
    precip_melt = MIN(snow(i),precip_melt)
    thdd_km1(i) = thdd_km1(i)-precip_melt/factor
    rain(i) = rain(i)+precip_melt
    snow(i) = snow(i)-precip_melt

  END IF

END DO

IF (lhook) CALL dr_hook('CRS_FRZL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE crs_frzl

