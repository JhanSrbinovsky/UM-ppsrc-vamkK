! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Calculates the amount of water loading for the parcel and the environment

MODULE water_loading_mod

IMPLICIT NONE

!
! Description:
!   Calculates the amount of water loading for the parcel and the environment
!
! Method:
!   Calculates the amount of water loading for the parcel and the environment
!   depending upon the water loading option selected.
!   1) Ignore water loading
!   2) Include all the condensate in the water loading
!   3) Not yet coded: Include only the convective condensate estimated to 
!      remain after precipitation
!
!   The arguments are named in terms of layer k but the routine can also 
!   be applied to layer k+1.
!
! Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


CONTAINS

SUBROUTINE water_loading (npnts, pk, exk, thek, thpk,       &
                          qclek, qcfek, qclpk, qcfpk,       &
                          watldek,watldpk)

USE cv_run_mod, ONLY: cnv_wat_load_opt
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine arguments:
!-----------------------------------------------------------------------

! Arguments with intent IN:

INTEGER, INTENT(IN) :: npnts ! Vector length

REAL,INTENT(IN) :: pk(npnts)    ! Pressure in layer k (Pa)
REAL,INTENT(IN) :: exk(npnts)   ! Exner ratio at mid-point of layer k
REAL,INTENT(IN) :: thek(npnts)  ! environment p. temperature in layer k (K)
REAL,INTENT(IN) :: thpk(npnts)  ! parcel p. temperature in layer k (K)
REAL,INTENT(IN) :: qclek(npnts) ! env. liquid condensate in layer k (kg/kg)
REAL,INTENT(IN) :: qcfek(npnts) ! env. frozen condensate in layer k (kg/kg)
REAL,INTENT(IN) :: qclpk(npnts) ! parcel liquid condensate in layer k (kg/kg)
REAL,INTENT(IN) :: qcfpk(npnts) ! parcel frozen condensate in layer k (kg/kg)

! Arguments with intent OUT:

REAL, INTENT(OUT) :: watldek(npnts) ! Environment water loading in layer k+1
REAL, INTENT(OUT) :: watldpk(npnts) ! Parcel water loading in layer k+1

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER :: i                 ! loop counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!---------------------------------------------------------------------

IF (lhook) CALL dr_hook('DET_RATE',zhook_in,zhook_handle)

SELECT CASE (cnv_wat_load_opt)
CASE (0)
  DO i=1,npnts
    watldek(i) = 0.0
    watldpk(i) = 0.0
  END DO
CASE (1)
  DO i=1,npnts
    watldek(i) = qclek(i) + qcfek(i)
    watldpk(i) = qclpk(i) + qcfpk(i)
  END DO
END SELECT  


IF (lhook) CALL dr_hook('WATER_LOADING',zhook_out,zhook_handle)
RETURN
END SUBROUTINE water_loading
END MODULE water_loading_mod
