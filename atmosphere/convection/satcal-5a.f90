! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate saturated temperature
!
! Subroutine Interface:  
!
SUBROUTINE satcal (npnts, k, t, th, pk, exk, q_k, the_k, qse_k, qs, thdds)

USE atmos_constants_mod, ONLY: cp, rv

USE water_constants_mod, ONLY: lc,  tm

USE cv_derived_constants_mod, ONLY: ls

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! 
! Description: Calculate saturated temperature
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
  t(npnts)           & ! Temperature (K)

 ,th(npnts)          & ! Potential temperature (K)

 ,pk(npnts)          & ! Pressure of layer k (Pa)

 ,exk(npnts)         & ! Exner ratio of layer k

 ,q_k(npnts)         & ! Mixing ratio of layer k (kg/kg)

 ,the_k(npnts)       & ! Environmental potential temperature in layer k
 ,qse_k(npnts)         ! qsat for environment of layer k (kg/kg)

REAL, INTENT(OUT) :: &
  qs(npnts)          & ! Saturated specific  humidity (kg/kg)

 ,thdds(npnts)         ! Saturated environmental potential temperature (K)

!-----------------------------------------------------------------------
! Model constants
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
! Local variables 
!-----------------------------------------------------------------------

INTEGER ::       &
  i, ic            ! Loop counters

REAL ::          &
  l                ! Latent heat

REAL ::          &
  ts(npnts)      & ! Saturated temperature (K)

 ,t_fg(npnts)    & ! Temperature first guess (K)

 ,th_fg(npnts)   & ! Potential temperature first guess (K)

 ,dqbydt(npnts)    ! First guess at mixing ratio increment (kg/kg/s)

REAL ::          &
  cpexk(npnts)     ! cp * exner 

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------
! Set initial first guess temperature and theta - based upon
! environmental temperature in layer k
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook('SATCAL',zhook_in,zhook_handle)

DO i=1,npnts
  th_fg(i) = the_k(i)
  t_fg(i)  = th_fg(i)*exk(i)
  qs(i)    = qse_k(i)        ! first guess qsat from environment value
  cpexk(i) = cp*exk(i)       ! save CPU by calculating before interation 
END DO

!----------------------------------------------------------------------
! Do two iterations to find saturation point due to evaporation
!----------------------------------------------------------------------

DO ic=1,2

  !----------------------------------------------------------------------
  ! Calculate dqsat/dT for first guess temperature
  !----------------------------------------------------------------------

! DEPENDS ON: dqs_dth
!  CALL dqs_dth(dqbydt,k,th_fg,qs,exk,npnts)
! Cheaper to remove call and do inline making use of l

  !----------------------------------------------------------------------
  ! Calculate updated temperature at saturation
  !----------------------------------------------------------------------

  DO i=1,npnts

    IF (t_fg(i) >  tm) THEN
      l=lc
    ELSE
      l=ls     ! lc+lf
    END IF
    !----------------------------------------------------------------------
    ! Calculate dqsat/dT for first guess temperature
    !----------------------------------------------------------------------
    dqbydt(i) = l*qs(i)/(exk(i)*rv*th_fg(i)*th_fg(i))

    thdds(i) = (th(i)*cpexk(i) - l*(qs(i)-q_k(i)-th_fg(i)*dqbydt(i))) /   &
                  (cpexk(i)+l*dqbydt(i))


    !----------------------------------------------------------------------
    ! Calculate temperature at saturation and update first guess
    !----------------------------------------------------------------------

    th_fg(i) = thdds(i)
    t_fg(i) = th_fg(i)*exk(i)

  END DO

  !----------------------------------------------------------------------
  ! Calculate revised saturation mixing ratio at saturation
  !---------------------------------------------------------------------

! DEPENDS ON: qsat
  CALL qsat(qs,t_fg,pk,npnts)

END DO

IF (lhook) CALL dr_hook('SATCAL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE satcal
