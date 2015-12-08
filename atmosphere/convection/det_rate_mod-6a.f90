! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Calculates forced detrainment rate in layer K

MODULE det_rate_6a_mod

IMPLICIT NONE

!
! Description:
!   Calculates forced detrainment rate in layer K
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


CONTAINS

SUBROUTINE det_rate_6a (npnts, qek, qekp1, qpk, qpkp1, qrk, Qlkp1, Qfkp1, &
                        ekp14, ekp34, deltak)
USE water_constants_mod, ONLY: lc, lf
USE atmos_constants_mod, ONLY: cp
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Subroutine arguments

!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------

INTEGER, INTENT(IN) :: npnts     ! Number of points

REAL,INTENT(IN) :: qek(npnts)    ! Env. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: qekp1(npnts)  ! Env. spec. humidity in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qpk(npnts)    ! Par. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: qpkp1(npnts)  ! Par. specific humidity in layer k+1 (kg/kg)
REAL,INTENT(IN) :: qrk(npnts)    ! Specific humidity of forced detrained
                                 ! parcel in layer k (kg/kg)
REAL,INTENT(IN) :: Qlkp1(npnts)  ! Amount of condensation to liquid water 
                                 ! in the parcel (kg/kg)
REAL,INTENT(IN) :: Qfkp1(npnts)  ! Amount of deposition to ice water
                                 ! in the parcel (kg/kg)
REAL,INTENT(IN) :: ekp14(npnts)  ! Entrainment coefficient at level k+1/4 
                                 ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: ekp34(npnts)  ! Entrainment coefficient at level k+3/4 
                                 ! multiplied by appropriate layer thickness
                                 
!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------

REAL,INTENT(OUT) :: deltak(npnts)   ! Parcel forced detrainment rate in 
                                    ! layer k multiplied by layer thickness

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER :: i        ! loop counter

REAL :: Factor      ! factor used to calculate the detrainment rate.
REAL :: Denom       ! denominator used to calculate the detrainment rate.
REAL, PARAMETER :: TinyDenom=TINY(Denom) !Smallest allowable denominator
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!---------------------------------------------------------------------

IF (lhook) CALL dr_hook('DET_RATE_6A',zhook_in,zhook_handle)

DO i=1,npnts
  Factor = ekp14(i)*qek(i) + (1.0+ekp14(i))*ekp34(i)*qekp1(i)               &
           - (1.0+ekp14(i))*(1.0+ekp34(i))*(qpkp1(i) + Qlkp1(i) + Qfkp1(i))   

  Denom  = qrk(i) + Factor

  IF (ABS(Denom) <= TinyDenom) THEN
    !If the denominator is zero then set the detrainment rate to 
    !one.
    deltak(i) = 1.0
  ELSE
    deltak(i) = (qpk(i) + Factor)/Denom
  END IF

  !Constrain the detrainment rate to be between zero and almost one.
  deltak(i) = MAX(MIN(1.0,deltak(i)),0.0)

END DO

IF (lhook) CALL dr_hook('DET_RATE_6A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE det_rate_6a
END MODULE det_rate_6a_mod
