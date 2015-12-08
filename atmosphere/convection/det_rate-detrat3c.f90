! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates forced detrainment rate in layer K
!
SUBROUTINE det_rate_4a5a (npnts,bwkp1,bcalc,                        &
                          thrk,xsqr,thpk,thek,thekp1,xsqkp1,thpkp1, &
                          ekp14,ekp34,exk,exkp1,                    &
                          deltak)
USE water_constants_mod, ONLY: lc, lf
USE atmos_constants_mod, ONLY: cp
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description : 
!   Calculates forced detrainment rate in layer K   
!
!   Method    : See Unified Model documentation paper 27.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.3.
!-----------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  npnts                  ! Vector length

LOGICAL, INTENT(IN) :: &
  bwkp1(npnts)         & ! Mask for whether condensate is liquid in layer k+1
 ,bcalc(npnts)           ! Mask for points at which calculations of this
                         ! routine are needed

REAL,INTENT(IN) :: &
  thrk(npnts)      & ! parcel detrainment potential temperature in layer k (K)
 ,xsqr(npnts)      & ! Excess water vapour of the detraining air in layer k
                     ! (kg/kg)
 ,thpk(npnts)      & ! parcel potential temperature in layer k (K)
 ,thek(npnts)      & ! potential temperature of cloud environment in layer k (K)
 ,thekp1(npnts)    & ! potential temperature of cloud environment in layer k+1
 ,xsqkp1(npnts)    & ! Excess water vaopur of the parcel in layer k+1 (kg/kg)
 ,thpkp1(npnts)    & ! Parcel poential temperature in layer k+1 (K)
 ,ekp14(npnts)     & ! Entrainment coefficient at level k+1/4 multiplied by 
                     ! appropriate layer thickness
 ,ekp34(npnts)     & ! Entrainment coefficient at level k+3/4 multiplied by 
                     ! appropriate layer thickness
 ,exk(npnts)       & ! Exner ratio at mid-point of layer k
 ,exkp1(npnts)       ! Exner ratio at mid-point of layer k+1

REAL, INTENT(OUT) :: &
  deltak(npnts)       ! parcel forced detrainment rate layer k multiplied
                      ! by appropriate layer thickness

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER ::        & 
  i                 ! loop counter

REAL ::           &
  el              & ! Latent heat of condensation or (condensation + fusion)
                    !  (J/kg)
 ,epss            & ! (1+EKP14)*(1+EKP34)
 ,denom             ! The denominator part of the deltak eqn

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!---------------------------------------------------------------------

IF (lhook) CALL dr_hook('DET_RATE_4A5A',zhook_in,zhook_handle)

DO i=1,npnts
  IF (bcalc(i)) THEN

  !-----------------------------------------------------------------------
  !   Create a vector of latent heats
  !-----------------------------------------------------------------------

  IF (bwkp1(i)) THEN
    el = lc
  ELSE
    el = lc + lf
  END IF

  !-----------------------------------------------------------------------
  !   Calculate detrainment rates
  !-----------------------------------------------------------------------

  epss = (1.0 + ekp14(i)) * (1.0 + ekp34(i))
  deltak(i) = ekp14(i)*thek(i) + ekp34(i)*(1.0+ekp14(i))*thekp1(i)         &
                 - epss*(thpkp1(i) - el/(exkp1(i)*cp) * xsqkp1(i))

  !      Only perform DELTAK equation if the denominator is not
  !      equal to zero.


  denom = (deltak(i) + thrk(i))*(exk(i)*cp) - el*xsqr(i)
  IF (denom /= 0.0) THEN
    deltak(i) = (deltak(i) + thpk(i)) * (exk(i)*cp) / denom
  ELSE
    deltak(i) = 1.0   ! Detrain all if denominator zero
  END IF


  !----------------------------------------------------------------------
  !  From a theoretical view point DELTAK cannot = 1 . However because of
  !  approximations used in the calculation numerically it may be possible.
  !  Hence if  DELTAK = 1 set it to slightly smaller than 1
  !----------------------------------------------------------------------

    deltak(i) = MIN(0.99999,deltak(i))

  END IF
END DO

IF (lhook) CALL dr_hook('DET_RATE_4A5A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE det_rate_4a5a
