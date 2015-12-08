! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    Purpose of this deck : provide a number of routines equivalent
!    to the random number generation routines of the NAG library.
!    Where possible, the interface of the routines/function
!    matches its equivalent in the NAG library.
!

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model

!
!    Programming standards :
!      Fortran 90, Written to UM coding standards
!      as specified in UMDP 3, vn8.2
!
!    Logical components :
!
!    System task :
!
!----------------------------------------------------------------
!
!     Initialize the random number generator RAN1_JC (and therefore
!     also the functions like g05xxx which use RAN_JC) to give
!     repeatable sequence
!
!     Argument : iseed input, not destroyed after the call
!     Side effect : initialize the /ran_jc/ common block
!
!     This routine *needs* to be called first when using the
!     random routines G05xxE

SUBROUTINE g05cbe(iseed)

  IMPLICIT NONE

! Argument
  INTEGER :: iseed

! Local variables
  INTEGER, PARAMETER :: ntab = 32

  INTEGER :: idum, iv(ntab), iy ! This common block constitutes the
  COMMON /ran_jc/ idum,iv,iy    ! 'memory' of the random timeserie.

  INTEGER :: j                  ! index


! ----------------------
! Zero the arrays of the ran1_jc function
  DO j=1, ntab
    iv(j) = 0
  END DO
  iy = 0

! initialize the ran1_jc function with the user provided value:
  idum = - INT(iseed)       ! idum must be negative when
                            ! initializing RAN1_JC.

  RETURN

END SUBROUTINE g05cbe

