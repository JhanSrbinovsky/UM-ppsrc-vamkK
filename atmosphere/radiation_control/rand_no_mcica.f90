! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Generate random numbers that are consistent across all machine
!
! Method:
!   Use linear congruential random number generator
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------

MODULE rand_no_mcica


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER, PARAMETER :: m = 86436
INTEGER, PARAMETER :: a = 1093
INTEGER, PARAMETER :: c = 18257

  CONTAINS

  SUBROUTINE mcica_rand_no(IntRand, Out,gi,gk)

! Subroutine arguments
  Integer, INTENT(IN) :: gi,gk
  Integer, INTENT(INOUT) :: IntRand(gi)

  REAL, INTENT(OUT) :: Out(gi,gk)

! Local variables
  INTEGER :: i
  INTEGER :: k
  REAL    :: rm
  REAL    :: temp

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('RAND_NO_MCICA',zhook_in,zhook_handle)

  rm=1.0/REAL(m)

  DO k=1,gk
    DO i=1,gi
      temp=IntRand(i)*a+c
      IntRand(i)=temp-INT(temp*rm)*m
      Out(i,k)=REAL(IntRand(i))*rm
    END DO
  END DO

  IF (lhook) CALL dr_hook('RAND_NO_MCICA',zhook_out,zhook_handle)
  RETURN
  END SUBROUTINE mcica_rand_no
END MODULE rand_no_mcica
