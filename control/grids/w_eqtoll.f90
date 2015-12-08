! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Subroutine W_EQTOLL----------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!
!    This file belongs in section: Grids
!    Purpose:  Calculates u and v components of wind on standard
!              latitude-longitude grid by rotating wind
!              components on equatorial latitude-longitude (eq)
!              grid.
!
!    Programming standard : UMDP3 v8.2
!
!   Logical components covered : S133
!
!    Documentation: The transformation formulae are described in
!                   unified model on-line documentation paper S1.
MODULE w_eqtoll_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE w_eqtoll(coeff1,coeff2,u_eq,v_eq,u,v,points,l_mdi)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN) :: points      ! Number of points to be processed

REAL, INTENT(IN)  :: coeff1(points) ! Coefficient of rotation no 1
REAL, INTENT(IN)  :: coeff2(points) ! Coefficient of rotation no 2
REAL, INTENT(IN)  :: u_eq(points)   ! u component of wind on equatorial grid
REAL, INTENT(IN)  :: v_eq(points)   ! v component of wind on equatorial grid
REAL, INTENT(OUT) :: u(points)      ! u component of wind on lat-lon grid
REAL, INTENT(OUT) :: v(points)      ! v component of wind on lat-lon grid

LOGICAL, INTENT(IN) :: l_mdi  ! Shall we check for MDIs?

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! Define local varables
INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!  1. Transform wind components

! Formulae used are from eq (4.14)

IF (lhook) CALL dr_hook('W_EQTOLL',zhook_in,zhook_handle)
IF (l_mdi) THEN
  DO i = 1,points
    IF ( u_eq(i)  ==  rmdi .OR. v_eq(i)  ==  rmdi ) THEN
      u(i) = rmdi
      v(i) = rmdi
    ELSE
      u(i)=coeff1(i)*u_eq(i)+coeff2(i)*v_eq(i)
      v(i)=coeff1(i)*v_eq(i)-coeff2(i)*u_eq(i)
    END IF
  END DO

ELSE

  DO i = 1,points
    u(i)=coeff1(i)*u_eq(i)+coeff2(i)*v_eq(i)
    v(i)=coeff1(i)*v_eq(i)-coeff2(i)*u_eq(i)
  END DO

END IF

IF (lhook) CALL dr_hook('W_EQTOLL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE w_eqtoll
END MODULE w_eqtoll_mod
