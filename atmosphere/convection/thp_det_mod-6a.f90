! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Calculates theta of the parcel in layer k+1 after forced detrainment
!
MODULE thp_det_6a_mod

IMPLICIT NONE

!
! Description:
!   Calculates potential temperature of the parcel in layer k+1 
!   after forced detrainment in layer k.
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

SUBROUTINE thp_det_6a(npnts, exkp1, pkp1, thekp1, qekp1, watldekp1, watldpkp1, &
                      xsbmin, bwkp1, bgmkp1, qpkp1, thpkp1)

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

INTEGER,INTENT(IN) :: npnts         ! Number of points

REAL,INTENT(IN) :: exkp1(npnts)     ! Exner ratio at mid-point of layer k+1
REAL,INTENT(IN) :: pkp1(npnts)      ! pressure at mid-point of layer k+1 (Pa)
REAL,INTENT(IN) :: thekp1(npnts)    ! Env. p. temperature in layer k+1 (K)
REAL,INTENT(IN) :: qekp1(npnts)     ! Env. spec. humidity in layer k+1 (kg/kg)
REAL,INTENT(IN) :: watldekp1(npnts) ! Env. water loading in layer k+1 (kg/kg)
REAL,INTENT(IN) :: watldpkp1(npnts) ! Par. water loading in layer k+1 (kg/kg)
REAL,INTENT(IN) :: xsbmin(npnts)    ! Threshold buoyancy for forced 
                                    ! detrainment (K)

LOGICAL,INTENT(IN) :: bwkp1(npnts)  ! Mask for whether condensate is 
                                    ! liquid in layer k+1
LOGICAL,INTENT(IN) :: bgmkp1(npnts) ! Mask for parcels which are 
                                    ! saturated in layer k+1

!----------------------------------------------------------------------
! Variables which are input and output
!----------------------------------------------------------------------

REAL,INTENT(INOUT) :: qpkp1(npnts)  ! Par. spec. humidity in layer k+1 (kg/kg)
                                    ! after forced detrainment

!----------------------------------------------------------------------
! Variables which are output 
!----------------------------------------------------------------------

REAL,INTENT(OUT) :: thpkp1(npnts)   ! Par. p. temperature in layer k+1 (K)
                                    ! after forced detrainment

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER :: i, j ! loop counters

REAL :: tpkp1(npnts) ! estimate of the detrained parcel's temperature (K)
REAL :: qspkp1(npnts)! qsat at tpkp1
REAL :: dqsdth  ! Rate of change of qsat with potential temperature
REAL :: el      ! Latent heat of gas-to-whatever-condenses PC2 defn (J/kg)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('THP_DET_6A',zhook_in,zhook_handle)

DO i=1,npnts
  IF ( bgmkp1(i) ) THEN
    ! If saturated then an iterative calculation is required
    ! Set the first estimate of the parcel p.temp to that of the 
    ! environment
    thpkp1(i) = thekp1(i)
    tpkp1(i)  = thpkp1(i) * exkp1(i)
  ELSE
    ! If not saturated the detrained p.temp and humidity can be directly
    ! calculated
    thpkp1(i) = ( thekp1(i)*(1.0 + c_virtual*qekp1(i) - watldekp1(i))     &
                + xsbmin(i) ) / (1.0 + c_virtual*qpkp1(i) - watldpkp1(i))
    !tpkp1 is set to prevent qsat operating on uninitialised data.
    tpkp1(i)  = thpkp1(i) * exkp1(i)
    ! The humidity is unchanged.
  END IF
END DO

! DEPENDS ON: qsat
CALL qsat (qspkp1, tpkp1, pkp1, npnts)

DO j=1,3  !Three iterations should be sufficient
  DO i=1,npnts
    IF ( bgmkp1(i) ) THEN !if saturated
      IF ( bwkp1(i) ) THEN
        el=lc
      ELSE
        el=ls
      END IF
      
      !Calculate the gradient of qsat. nb qspkp1 is saturated humidity at tpkp1
      dqsdth = el * qspkp1(i) / ( rv * exkp1(i) * thpkp1(i) * thpkp1(i) )

      !Calculate the next estimate of the parcel's p.temp
      thpkp1(i) = ( thekp1(i)*(1.0 + c_virtual*qekp1(i) - watldekp1(i))        &
                - thpkp1(i)*c_virtual*(qspkp1(i) - thpkp1(i)*dqsdth) &
                + xsbmin(i) )                                                  &
                / ( 1.0 + c_virtual*thpkp1(i)*dqsdth - watldpkp1(i) )
        
      !update the parcel's temperature
      tpkp1(i)  = thpkp1(i) * exkp1(i)
    END IF    !if saturated
  END DO      !i

! Assume that the detrained parcel is saturated and update its humidity.
! DEPENDS ON: qsat
CALL qsat (qspkp1, tpkp1, pkp1, npnts)

END DO  !j

DO i=1,npnts
  IF ( bgmkp1(i) ) THEN
    ! If saturated update qpkp1 to saturated humidity at tpkp1
    qpkp1(i) = qspkp1(i)
  END IF
END DO

IF (lhook) CALL dr_hook('THP_DET_6A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE thp_det_6a
END MODULE thp_det_6a_mod
