! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates theta and q of detraining air in layer k
!
MODULE thetar_6a_mod

IMPLICIT NONE

!
! Description:
!   Calculates the potential temperature and the humidity of the parcel 
!   undergoing forced detrainment in layer k
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

SUBROUTINE thetar_6a(npnts, exk, pk, thek, qek, qpk, watldek, watldpk, &
                     bwk, bgmk, thrk, qrk)

USE water_constants_mod, ONLY: lc, lf
USE cv_derived_constants_mod, ONLY: ls
USE atmos_constants_mod, ONLY: c_virtual, rv

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine arguments

!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------

INTEGER,INTENT(IN) :: npnts        ! Number of points

REAL,INTENT(IN) :: exk(npnts)      ! Exner ratio at mid-point of layer k
REAL,INTENT(IN) :: pk(npnts)       ! pressure at mid-point of layer k (Pa)
REAL,INTENT(IN) :: thek(npnts)     ! Env. p. temperature in layer k (K)
REAL,INTENT(IN) :: qek(npnts)      ! Env. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: qpk(npnts)      ! Par. specific humidity in layer k (kg/kg)
REAL,INTENT(IN) :: watldek(npnts)  ! Env. water loading in layer k (kg/kg)
REAL,INTENT(IN) :: watldpk(npnts)  ! Par. water loading in layer k (kg/kg)

LOGICAL,INTENT(IN) :: bwk(npnts)   ! Mask for whether condensate is      
                                   ! liquid in layer k+1
LOGICAL,INTENT(IN) :: bgmk(npnts)  ! Mask for parcels which are
                                   ! saturated in layer k

!----------------------------------------------------------------------
! Variables which are output 
!----------------------------------------------------------------------

REAL,INTENT(OUT) :: thrk(npnts)    ! Pot. temperature of forced detrained
                                   ! parcel in layer k (K)
REAL,INTENT(OUT) :: qrk(npnts)     ! Specific humidity of forced detrained
                                   ! parcel in layer k (kg/kg)

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER :: i, j ! loop counters

REAL :: trk(npnts) ! estimate of the detrained parcel's temperature (K)
REAL :: qsrk(npnts)! qsat at trk
REAL :: dqsdth  ! Rate of change of qsat with potential temperature
REAL :: el      ! Latent heat of gas-to-whatever-condenses PC2 defn (J/kg)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!----------------------------------------------------------------------
IF (lhook) CALL dr_hook('THETAR_6A',zhook_in,zhook_handle)

DO i=1,npnts
  IF ( bgmk(i) ) THEN
    ! If saturated then an iterative calculation is required
    ! Set the first estimate of the detrained p.temp to that of the 
    ! environment
    thrk(i) = thek(i)
  ELSE
    ! If not saturated the detrained p.temp and humidity can be directly
    ! calculated
    thrk(i) = thek(i)*(1.0 + c_virtual*qek(i) - watldek(i))               &
              /(1.0 + c_virtual*qpk(i) - watldpk(i))
    qrk(i)  = qpk(i)
  END IF
  !update the parcel's temperature
  trk(i)  = thrk(i) * exk(i)
END DO

! DEPENDS ON: qsat
CALL qsat (qsrk, trk, pk, npnts)

DO j=1,3  !Three iterations should be sufficient
  DO i=1,npnts
    IF ( bgmk(i) ) THEN !if saturated
      IF ( bwk(i) ) THEN
        el=lc
      ELSE
        el=ls
      END IF
      
      !Calculate the gradient of qsat. nb qsrk is saturated humidity at trk
      dqsdth = el * qsrk(i) / ( rv * exk(i) * thrk(i) * thrk(i) )

      !Calculate the next estimate of the parcel's p.temp
      thrk(i) = ( thek(i)*(1.0 + c_virtual*qek(i) - watldek(i))           &
                - thrk(i)*c_virtual*(qsrk(i) - thrk(i)*dqsdth) )          &
                / ( 1.0 + c_virtual*thrk(i)*dqsdth - watldpk(i) )
        
      !update the parcel's temperature
      trk(i)  = thrk(i) * exk(i)
    END IF    !if saturated
  END DO      !i

! Assume that the detrained parcel is saturated and update its humidity.
! DEPENDS ON: qsat
  CALL qsat (qsrk, trk, pk, npnts)

END DO  !j

DO i=1,npnts
  IF ( bgmk(i) ) THEN
    ! If saturated update qrk to saturated humidity at trk
    qrk(i) = qsrk(i)
  END IF
END DO

IF (lhook) CALL dr_hook('THETAR_6A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE thetar_6a
END MODULE thetar_6a_mod
