! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculates constants used in large-scale precipitation scheme.
!  SUBROUTINE GAMMAF--------------------------------------------------
!   PURPOSE: CALCULATES COMPLETE GAMMAF FUNCTION BY
!   A POLYNOMIAL APPROXIMATION
!
!   Notes: The original routine was only designed for arguments >1.
!          From VN7.6 this has modified this to work with arguments <1
!          following guidance from Ben Shipway.
!
!          NB For negative integers the solution is undefined and
!          a divide by zero will occur, but UMUI settings should prevent
!          this happening.
! --------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation
MODULE gammaf_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE gammaf(y,gam)

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE
  REAL ::                                                               &
                              !, INTENT(IN)
    y
  REAL ::                                                               &
                              !, INTENT(OUT)
    gam
! Gamma function of Y

! LOCAL VARIABLE
  INTEGER :: i,m
  REAL :: gg,g,pare,x

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

! --------------------------------------------------------------------
  IF (lhook) CALL dr_hook('GAMMAF',zhook_in,zhook_handle)
  gg=1.
  m=FLOOR(y)
  x=y-m
  IF (m > 1) THEN
    DO i = 1, m-1
      g=y-i
      gg=gg*g
    END DO
  ELSE IF (m < 1)THEN
    DO i = m, 0
      g=y-i
      gg=gg/g
    END DO
  END IF
  pare=-0.5748646*x+0.9512363*x*x-0.6998588*x*x*x                       &
  +0.4245549*x*x*x*x-0.1010678*x*x*x*x*x+1.
  gam=pare*gg
  IF (lhook) CALL dr_hook('GAMMAF',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE gammaf

! --------------------------------------------------------------------
! --------------------------------------------------------------------

SUBROUTINE gammaf_lower_incomp(s_in,x_in,gamma_out,nn)

  USE yomhook,    ONLY: lhook, dr_hook
  USE parkind1,   ONLY: jprb, jpim

  IMPLICIT NONE

! Input parameters
  REAL, INTENT(IN)    :: s_in,x_in
! Output value of lower incomplete gamma function
  REAL, INTENT(OUT)   :: gamma_out
! Number of elements in summation.
  INTEGER, INTENT(IN) :: nn

! Local variables
  INTEGER :: k
  REAL :: factor, BigGamma, tmpsum

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

! --------------------------------------------------------------------
  IF (lhook) CALL dr_hook('GAMMAF_LOWER_INCOMP',zhook_in,zhook_handle)

  CALL gammaf(s_in, BigGamma)

  factor = (x_in**s_in) * BigGamma * EXP(-x_in)

  tmpsum = 0.0

  DO k = 0, nn
    CALL gammaf(s_in+REAL(k)+1.0, BigGamma)

    tmpsum = tmpsum + ( (x_in**k) / BigGamma )

  END DO

  gamma_out = factor * tmpsum

  IF (lhook) CALL dr_hook('GAMMAF_LOWER_INCOMP',zhook_out,zhook_handle)

  RETURN

END SUBROUTINE gammaf_lower_incomp

END MODULE gammaf_mod
