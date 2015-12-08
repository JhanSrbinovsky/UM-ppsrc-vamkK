! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Subroutine W_LLTOEQ----------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Grids
!
!    Purpose:  Calculates u and v components of wind on equatorial
!              (eq) latitude longitude grid by rotating wind
!              components on standard latitude-longitude (eq)
!              grid.
!
!    Programming standard : UMDP3 v8.2
!
!    Documentation: The transformation formulae are described in
!                   unified model on-line documentation paper S1.
!      -----------------------------------------------------------------
MODULE w_lltoeq_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE w_lltoeq(coeff1,coeff2,u,v,u_eq,v_eq,points,l_mdi)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, INTENT(IN) :: points !  Number of points to be processed

REAL, INTENT(IN)  :: coeff1(points) ! Coefficient of rotation no 1
REAL, INTENT(IN)  :: coeff2(points) ! Coefficient of rotation no 2
REAL, INTENT(OUT) :: u_eq(points)   ! u component of wind on equatorial grid
REAL, INTENT(OUT) :: v_eq(points)   ! v component of wind on equatorial grid
REAL, INTENT(IN)  :: u(points)      ! u component of wind on lat-lon grid
REAL, INTENT(IN)  :: v(points)      ! v component of wind on lat-lon grid

LOGICAL, INTENT(IN) :: l_mdi     ! Shall we check for MDIs?
! Workspace usage:-----------------------------------------------------
! None
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

! Define local varables:-----------------------------------------------
INTEGER :: i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------

!  1. Transform wind components

! Formulae used are from eq (4.13)

IF (lhook) CALL dr_hook('W_LLTOEQ',zhook_in,zhook_handle)
IF (l_mdi) THEN
  DO i = 1,points
    IF ( u(i)  ==  rmdi .OR. v(i)  ==  rmdi ) THEN
      u_eq(i) = rmdi
      v_eq(i) = rmdi
    ELSE
      u_eq(i)=coeff1(i)*u(i)-coeff2(i)*v(i)
      v_eq(i)=coeff1(i)*v(i)+coeff2(i)*u(i)
    END IF
  END DO

ELSE

  DO i = 1,points
    u_eq(i)=coeff1(i)*u(i)-coeff2(i)*v(i)
    v_eq(i)=coeff1(i)*v(i)+coeff2(i)*u(i)
  END DO

END IF

IF (lhook) CALL dr_hook('W_LLTOEQ',zhook_out,zhook_handle)
RETURN
END SUBROUTINE w_lltoeq
END MODULE w_lltoeq_mod
