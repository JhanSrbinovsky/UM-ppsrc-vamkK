! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Returns a mask for points at which convection is terminating
!
SUBROUTINE term_con_4a5a(npnts,nlev,k,new_termc,                     &
                         bwkp1,                                         &
                         flxkp1,thekp1,qekp1,thpi,qpi,qsekp1,deltak,    &
                         expi,ekp14,ekp34,pstar,pk,pkp1,xsbmin,         &
                         bterm  )

USE cv_run_mod, ONLY:                                              &
    qstice
USE cv_param_mod, ONLY:                                            &
    mparfl
USE water_constants_mod, ONLY: lc, lf
USE atmos_constants_mod, ONLY: cp, c_virtual
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
! Returns a mask for points at which convection is terminating
!  See UM Documentation paper No 27
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in Section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  npnts                & ! Number of points
 ,nlev                 & ! Number of model levels for calculations
 ,k                    & ! present model layer
 ,new_termc              ! Flag for simplified termination of convection

LOGICAL, INTENT(IN) :: &  
  bwkp1(npnts)           ! Mask for whether condensate is liquid in layer k+1

REAL,INTENT(IN) :: &
  flxkp1(npnts)    & ! parcel massflux in layer k+1 (Pa/s)
 ,thekp1(npnts)    & ! potential temperature of cloud environment in layer k+1
 ,qekp1(npnts)     & ! mixing ratio of cloud environment in layer k+1 (kg/kg)
 ,thpi(npnts)      & ! Initial parcel potential temperature (K)
 ,qpi(npnts)       & ! Initial parcel mixing ratio (kg/kg)
 ,qsekp1(npnts)    & ! saturation mixing ratio of cloud environment
                     ! in layer k+1 (kg/kg)
 ,deltak(npnts)    & ! parcel forced detrainment rate layer k multiplied
                     ! by appropriate layer thickness
 ,expi(npnts)      & ! Initial parcel Exner pressure
 ,ekp14(npnts)     & ! Entrainment coefficient at level k+1/4 multiplied by 
                     ! appropriate layer thickness
 ,ekp34(npnts)     & ! Entrainment coefficient at level k+3/4 multiplied by 
                     ! appropriate layer thickness
 ,pstar(npnts)     & ! Surface pressure (Pa)
 ,pk(npnts)        & ! pressure at mid-point of layer k (Pa)
 ,pkp1(npnts)      & ! pressure at mid-point of layer k+1   (Pa)
 ,xsbmin(npnts)      ! Threshold buoyancy for forced detrainment
                     ! Function of delta P

!---------------------------------------------------------------------
! Variables which are input and output
!---------------------------------------------------------------------
LOGICAL, INTENT(INOUT) :: &  
  bterm(npnts)            ! Mask for parcels which terminate in layer k+1

!-------------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

INTEGER ::        & 
  i                 ! loop counter 

REAL ::        &
  el           & ! Latent heat of condensation or (condensation + fusion) (J/kg)
 ,flxmin       & ! Minimum convecive mass flux below which termination
                 ! detrainment occurs (Pa/s)
 ,thvundi      & ! Potential temperature of an undilute parcel in layer k+1
                 ! from the starting layer of convection (K)
 ,thvekp1        ! Virtual potential temperature of environment in layer k+1 (K)


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('TERM_CON_4A5A',zhook_in,zhook_handle)
!----------------------------------------------------------------------
!  Calculate minimum mass flux below which convection is terminated
!----------------------------------------------------------------------

DO i=1,npnts
  flxmin = mparfl*(1.0+ekp14(i))*(1.0+ekp34(i))*pstar(i)

  !-----------------------------------------------------------------------
  !   Create a vector of latent heats
  !-----------------------------------------------------------------------

  IF (bwkp1(i)) THEN
    el = lc
  ELSE
    el = lc + lf
  END IF

  ! ----------------------------------------------------------------------
  !   Parcels are only considered for termination if they are detraining
  !   except at the top model layer, where all convection terminates.
  !
  !   If the parcel has a potential temperature greater than the
  !   potential temperature of an undilute parcel from the starting
  !   layer of convection in layer k+1 then convection is terminated.
  !
  !   UM Documentation paper P27
  !   Section (7), equation (32)
  !
  !   Convection is also terminated if mass flux in layer k+1 is less than 
  !   a minimum value.
  !
  !   UM Documentation paper P27
  !   Section (7), equation (33)
  ! ----------------------------------------------------------------------

  thvundi=( thpi(i) + (el/(expi(i)*cp)) *(qpi(i) - qsekp1(i))       &
           +((lc-el)/(expi(i)*cp))*MAX(0.0,(qpi(i)-qstice))         &
           )*(1.0+c_virtual*qsekp1(i))



  thvekp1 = (thekp1(i)*(1.0+c_virtual*qekp1(i)) + xsbmin(i))


  IF (.NOT. bterm(i)) THEN
    ! Depending on whether option has been chosen in UMUI, either use
    ! original 4a termination condition or Martin Willett's simplified
    ! termination condition (new_termc=1)
    IF (new_termc  ==  1) THEN
      bterm(i) = (flxkp1(i)  <   flxmin) .OR. ((k+1)  ==  nlev)
    ELSE
      bterm(i) = (((flxkp1(i)  <   flxmin) .OR. (thvundi  <   thvekp1))   &
                 .AND. (deltak(i) >  0.0)) .OR. (k+1)  ==  nlev
    END IF
  END IF

END DO  ! i

IF (lhook) CALL dr_hook('TERM_CON_4A5A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE term_con_4a5a
