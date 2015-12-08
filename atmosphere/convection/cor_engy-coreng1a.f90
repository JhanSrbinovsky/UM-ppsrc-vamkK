! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Ensure conservation of moist static energy
!
SUBROUTINE cor_engy_4a(np_field,npnts,ncore,nlev,dthbydt,dqbydt,snow,         &
                       exner_layer_centres, p_layer_boundaries, index4)

USE earth_constants_mod, ONLY: g
USE water_constants_mod, ONLY: lc, lf
USE atmos_constants_mod, ONLY: cp
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
! Adjust the potential temperature increments to ensure the conservation
! of moist static energy
!
!  See UM Documentation paper No 27
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) :: &
  np_field             & ! Length of data 
 ,npnts                & ! Full number of points
 ,ncore                & ! Vector lengths
 ,nlev                   ! Number of model layers


INTEGER, INTENT(IN) :: &
  index4(npnts)          ! index of points

REAL,INTENT(IN) ::      &
  dqbydt(np_field,nlev)   ! Increment to model mixing ratio dur to convection
                          !    (kg/kg/s)

REAL,INTENT(IN) ::      &
  snow(np_field)         ! Snow at surface (kg/m**2/s)


REAL,INTENT(IN) ::                      &
  exner_layer_centres(np_field,0:nlev)  & ! Exner ratio
 ,p_layer_boundaries(np_field,0:nlev)     ! layer boundary (Pa)

!----------------------------------------------------------------------
! Variables which are input but which are also updated in this routine
!----------------------------------------------------------------------

REAL, INTENT(INOUT) ::   &  
  dthbydt(np_field,nlev)   ! IN increments to model potential temperature
                           !    dur to convection (K/s)
                           ! OUT corrected increments to model potential 
                           !    temperature due to convection (K/s) 

!-------------------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

INTEGER ::        & 
  i,k               ! loop counters 

LOGICAL ::        &
  bposer(ncore)   & ! Mask for points in layer k at which increments to model
                    ! potential temperature due to convection are positive
 ,bcorr(ncore)      ! Mask for points at which enthalpy correction in
                    ! necessary

REAL ::         &
  qsum(ncore)   & ! Summation of increments to model mixing ratio due to
                  ! convection in the vertical, weighted according to the
                  ! mass of the layer (kg/m**2/s)
 ,tspos(ncore)  & ! Summation of positive increments to model potential
                  ! temperature due to convection with height,weighted 
                  ! according to the mass of the layer (kg/m**2/s)
 ,tsneg(ncore)  & ! Summation of negative increments to model potential
                  ! temperature due to convection with height,weighted 
                  ! according to the mass of the layer (kg/m**2/s)
 ,terr(ncore)     ! Summation of all increments to model potential
                  ! temperature due to convection with height,weighted 
                  ! according to the mass of the layer (kg/m**2/s)

REAL ::     &
  delpk     & ! Difference in pressure across a layer (Pa)
 ,extempk     ! Exner ratio at the mid-point of layer k

REAL ::  &
  pu,pl

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------
! Model constants
!----------------------------------------------------------------------


IF (lhook) CALL dr_hook('COR_ENGY_4A',zhook_in,zhook_handle)
! ----------------------------------------------------------------------
! Sum up mixing ratio and  +ve and -ve temperature increments
! ----------------------------------------------------------------------

DO i=1,ncore
  qsum (i) = 0.0
  tspos(i) = 0.0
  tsneg(i) = 0.0
END DO

DO k=1,nlev
  DO i=1,ncore

    delpk = - (p_layer_boundaries(index4(i),k)  -                          &
                                        p_layer_boundaries(index4(i),k-1))

    extempk  = exner_layer_centres(index4(i),k)

    qsum(i) = qsum(i) + dqbydt(index4(i),k)*delpk

    IF (dthbydt(index4(i),k)  >   0.0) THEN
      tspos(i) = tspos(i) + dthbydt(index4(i),k)*(cp*delpk*extempk)
    ELSE
      tsneg(i) = tsneg(i) + dthbydt(index4(i),k)*(cp*delpk*extempk)
    END IF

  END DO
END DO

! ----------------------------------------------------------------------
!   Calculate the error and apply the necessary correction
!
!   UM DOCUMENTATION PAPER 27
!   Section (12), Equation (48), (49)
! ----------------------------------------------------------------------

DO i=1,ncore

  terr(i) = lc*qsum(i) - lf*g*snow(index4(i)) + tspos(i) + tsneg(i)

  bposer(i) = terr(i)  >   0.0

  IF (bposer(i) .AND. (tspos(i)  ==  0.0)) THEN
    bposer(i) = .FALSE.
  ELSE IF (.NOT.bposer(i) .AND. (tsneg(i)  ==  0.0)) THEN
    bposer(i) = .TRUE.
  END IF

  bcorr(i) = (tspos(i)  /=  0.0) .OR. (tsneg(i)  /=  0.0)

  IF (bposer(i) .AND. bcorr(i)) THEN
    terr(i) = 1.0 - terr(i)/tspos(i)
  ELSE IF (.NOT.bposer(i) .AND. bcorr(i)) THEN
    terr(i) = 1.0 - terr(i)/tsneg(i)
  END IF

END DO

DO k=1,nlev
!DIR$ IVDEP
  DO i=1,ncore
    IF (bcorr(i) .AND. (( bposer(i) .AND. (dthbydt(index4(i),k)  >   0.0))  &
        .OR. ( .NOT.bposer(i) .AND. (dthbydt(index4(i),k)  <   0.0))))  THEN

       dthbydt(index4(i),k) = dthbydt(index4(i),k)*terr(i)

    END IF
  END DO  ! ncore
END DO ! nlev

IF (lhook) CALL dr_hook('COR_ENGY_4A',zhook_out,zhook_handle)

RETURN
END SUBROUTINE cor_engy_4a
