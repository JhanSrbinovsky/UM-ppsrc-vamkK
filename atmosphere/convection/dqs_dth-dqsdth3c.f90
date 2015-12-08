! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates gradient of saturation mixing ratio
!
! Subroutine Interface: (argument list order does not obey standards
!                        but called in a lot of places so not altered.)
!
SUBROUTINE dqs_dth_4a (dqs,k,thek,qsek,exk,npnts)

USE water_constants_mod, ONLY: lc, lf
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE rv_mod, ONLY: rv
IMPLICIT NONE

! 
! Description: Calculates gradient of saturation mixing ratio with 
!              potential temperature from the Clausius-Clapeyron equation
!
! Method: UM documentation paper 27
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards vn8.2.
!----------------------------------------------------------------------
! Subroutine arguments
!----------------------------------------------------------------------


INTEGER, INTENT(IN) :: &
  npnts                &  ! Vector length
 ,k                       ! Present model layer

REAL, INTENT(IN) ::  &
  thek(npnts)        & ! Environmental potential temperature in layer k

 ,qsek(npnts)        & ! Saturation Mixing ratio of layer k (kg/kg)

 ,exk(npnts)           ! Exner ratio of layer k

REAL, INTENT(OUT) :: &
  dqs(npnts)           ! Gradient of saturation mixing ratio with potential 
                       ! temperature for layer k (kg/kg/K)

!----------------------------------------------------------------------
! Model constants
!----------------------------------------------------------------------


!----------------------------------------------------------------------
! Local variables 
!----------------------------------------------------------------------
INTEGER ::       &
  i                ! Loop counter


REAL ::          &
  el               ! Latent heating of condensation or condensation plus
                   ! heating (J/kg)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!----------------------------------------------------------------------

IF (lhook) CALL dr_hook('DQS_DTH_4A',zhook_in,zhook_handle)

DO i=1,npnts

  !-----------------------------------------------------------------------
  ! Create a vector of latent heats according to whether qsat is with
  ! respect to ice or water  (NOTE equation below should really have 
  ! TM=273.15 not 273.0? Current code is inconstistent.)
  !-----------------------------------------------------------------------

  IF (thek(i)*exk(i)  >   273.) THEN
    el = lc
  ELSE
    el = lc + lf
  END IF

  !-----------------------------------------------------------------------
  ! Calculate d(qsat)/d(theta)
  !-----------------------------------------------------------------------

  dqs(i) = el*qsek(i)/(exk(i)*rv*thek(i)*thek(i))

END DO

IF (lhook) CALL dr_hook('DQS_DTH_4A',zhook_out,zhook_handle)
RETURN
END SUBROUTINE dqs_dth_4a
