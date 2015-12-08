! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates theta of the parcel in layer k+1 after forced detrainment
!
SUBROUTINE thp_det_4a5a (npnts,bgmkp1,bcalc,                            &
                         thekp1,qpkp1,qekp1,qsekp1,dqskp1,xsbmin,       &
                         thpkp1)

USE atmos_constants_mod, ONLY: c_virtual
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
! Calculates potential temperature of the parcel in layer k+1 
! after forced detrainment in layer k.
!  See UM Documentation paper No. 27 Section (6) equation (28)
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
  bgmkp1(npnts)        & ! Mask for parcels which are saturated in layer k+1
 ,bcalc(npnts)           ! Mask for parcels at which calculations of this 
                         ! subroutine are to be carried out

REAL,INTENT(IN) :: &
  thekp1(npnts)    & ! potential temperature of cloud environment in layer k+1
 ,qpkp1(npnts)     & ! parcel mixing ratio in layer k+1 (kg/kg)
 ,qekp1(npnts)     & ! mixing ratio of cloud environment in layer k+1 (kg/kg)
 ,qsekp1(npnts)    & ! saturation mixing ratio of cloud environment
                     ! in layer k+1 (kg/kg)
 ,dqskp1(npnts)    & ! gradient of saturation mixing ratio with potential  
                     ! temperature for the cloud environment in layer k+1
                     ! (kg/kg/K)
 ,xsbmin(npnts)      ! Threshold buoyancy for forced detrainment
                     ! Function of delta P

REAL,INTENT(INOUT) ::  &
  thpkp1(npnts)        ! parcel potential temperature in layer k+1 (K)
                       ! after forced detrainment
!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

INTEGER ::        & 
  i                 ! loop counter

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('THP_DET_4A5A',zhook_in,zhook_handle)
! ---------------------------------------------------------------------
!   No significant structure
! ---------------------------------------------------------------------
!

DO i=1,npnts
  IF (bcalc(i))THEN
    IF (bgmkp1(i)) THEN
      thpkp1(i) = thekp1(i) +                                             &
                (c_virtual*thekp1(i)*(qekp1(i)-qsekp1(i)) + xsbmin(i))    &
                 /( 1.0 + c_virtual*thekp1(i)*dqskp1(i) )

    ELSE
      thpkp1(i) = (thekp1(i)*(1.0 + c_virtual*qekp1(i))+ xsbmin(i))       &
                                   /(1.0 + c_virtual*qpkp1(i))
    END IF
  END IF
END DO

IF (lhook) CALL dr_hook('THP_DET_4A5A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE thp_det_4a5a
