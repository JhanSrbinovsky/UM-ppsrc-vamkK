! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates theta and q of detraining air in layer k
!

SUBROUTINE thetar_4a5a (npnts,bwkp1,thek,qek,qpk,qsek,dqsk,exk,pk,        &
                        bgmk,thrk,qrk,xsqr)

USE atmos_constants_mod, ONLY: c_virtual
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!
! Calculates the potential temperature of the detraining air in layer k
! and also the difference in the water vapour content of the detraining air
! from that of the mean parcel in layer k.
!  See UM Documentation paper No. 27
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  npnts                  ! Number of points

LOGICAL, INTENT(IN) :: &
  bwkp1(npnts)           ! Mask for whether condensate is liquid in layer k+1

REAL, INTENT(IN) :: &
  thek(npnts)       & ! potential temperature of cloud environment in layer k(K)
 ,qek(npnts)        & ! mixing ratio of cloud environment in layer k (kg/kg)
 ,qpk(npnts)        & ! parcel mixing ratio in layer k (kg/kg)
 ,qsek(npnts)       & ! saturation mixing ratio of cloud environment 
                      ! in layer k (kg/kg)
 ,dqsk(npnts)       & ! gradient of saturation mixing ratio with potential 
                      ! temperature for the cloud environment of layer k
                      ! (kg/kg/K)
 ,exk(npnts)        & ! Exner ratio at mid-point of layer k
 ,pk(npnts)           ! pressure at mid-point of layer k (Pa)


!----------------------------------------------------------------------
! Variables which are input/output 
!----------------------------------------------------------------------

LOGICAL, INTENT(INOUT) :: &
  bgmk(npnts)               ! Mask for parcels which are saturated in layer k

!----------------------------------------------------------------------
! Variables which are output 
!----------------------------------------------------------------------

REAL, INTENT(OUT) :: &
  thrk(npnts)        & ! parcel detrainment potential temperature in layer k (K)
 ,qrk(npnts)         & ! parcel detrainment mixing ratio in layer k (kg/kg)
 ,xsqr(npnts)          ! Excess water vapour of detraining air (kg/kg)

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER ::  & 
  i           ! loop counter

REAL ::     &
  tt(npnts)   ! temporary store for temperature for the calculation of
              ! saturated mixing ratio (K)   

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle



!----------------------------------------------------------------------
IF (lhook) CALL dr_hook('THETAR_4A5A',zhook_in,zhook_handle)

DO i=1,npnts

! ----------------------------------------------------------------------
!   Calculate the potential temperature of detraining air
!
!   UM Documentation paper P27
!   section (6), equation(26)
! ----------------------------------------------------------------------

  IF (.NOT.bgmk(i)) THEN
    thrk(i) = thek(i)*(1.0 + c_virtual*qek(i)) /(1.0 + c_virtual*qpk(i))
  ELSE
    thrk(i) = thek(i)*(1.0 + c_virtual*(qek(i)-qsek(i))/          &
                      (1.0 + c_virtual*thek(i)*dqsk(i)))
  END IF               

! ----------------------------------------------------------------------
!   Calculate the mixing ratio of the detraining air. The
!   difference between this and the mixing ratio of the mean
!   parcel in layer k
!
!   The moisture difference is used calculate the
!   COND_DET_K term of equation (30), section (6),
!   UM Documentation paper P27
! ----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Convert potential temperature to temperature and calculate
! pressure of layer k for calculation of saturated mixing ratio.
!-----------------------------------------------------------------------

  tt(i) = thrk(i)*exk(i)
END DO     !i

! DEPENDS ON: qsat
CALL qsat (xsqr,tt,pk,npnts)

DO i=1,npnts
! ----------------------------------------------------------------------
!   Small numerical approximations in the above calculations can mean
!   that the detraining parcel is no longer saturated at THRK. Add a 
!   check to see if the parcel is still saturated, and reset BGMK to
!   FALSE if it is not.
! ---------------------------------------------------------------------
  IF (xsqr(i) >  qpk(i)) THEN
    bgmk(i)=.FALSE.
  END IF

  IF (bgmk(i)) THEN
    qrk(i)  = xsqr(i)
    xsqr(i) = qpk(i) - xsqr(i)
  ELSE
    qrk(i)  = qpk(i)
    xsqr(i) = 0.0
  END IF
END DO

IF (lhook) CALL dr_hook('THETAR_4A5A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE thetar_4a5a
