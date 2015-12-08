! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Returns a mask for points at which convection is terminating
!

MODULE term_con_6a_mod

IMPLICIT NONE

!
! Description:
!   Returns a mask for points at which convection is terminating!
! THIS ROUTINE IS CURRENTLY NOT USED BY CONVECTION DESPITE ITS NAME.
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

SUBROUTINE term_con_6a(npnts,nlev,k,flxkp1,ekp14,ekp34,pstar,bterm)

USE cv_param_mod, ONLY: mparfl
USE water_constants_mod, ONLY: lc, lf
USE atmos_constants_mod, ONLY: cp, c_virtual
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER,INTENT(IN) :: npnts         ! Number of points
INTEGER,INTENT(IN) :: nlev          ! Number of model levels for calculations
INTEGER,INTENT(IN) :: k             ! present model layer

REAL,INTENT(IN) :: flxkp1(npnts)    ! parcel massflux in layer k+1 (Pa/s)
REAL,INTENT(IN) :: ekp14(npnts)     ! Entrainment coefficient at level k+1/4
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: ekp34(npnts)     ! Entrainment coefficient at level k+3/4 
                                    ! multiplied by appropriate layer thickness
REAL,INTENT(IN) :: pstar(npnts)     ! Surface pressure (Pa)

!---------------------------------------------------------------------
! Variables which are input and output
!---------------------------------------------------------------------
LOGICAL,INTENT(INOUT) :: bterm(npnts)   ! Mask for parcels which terminate
                                        ! in layer k+1

!-------------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

INTEGER ::   i                 ! loop counter 
REAL :: flxmin

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('TERM_CON_6A',zhook_in,zhook_handle)
!----------------------------------------------------------------------
!  Calculate minimum mass flux below which convection is terminated
!----------------------------------------------------------------------

DO i=1,npnts
  flxmin = mparfl*(1.0+ekp14(i))*(1.0+ekp34(i))*pstar(i)
  IF (.NOT. bterm(i)) THEN
    bterm(i) = (flxkp1(i)  <   flxmin) .OR. ((k+1)  ==  nlev)
  END IF
END DO  ! i

IF (lhook) CALL dr_hook('TERM_CON_6A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE term_con_6a

END MODULE term_con_6a_mod
