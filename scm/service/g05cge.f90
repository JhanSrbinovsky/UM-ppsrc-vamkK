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
!     Restore the state of the random number generator
!
!     Arguments : they are all input values.
!     They are the values which should be extracted
!     from the memory of the random routines by G05CFE
!     The arguments are not destroyed after the call.
!     Side effect : the /ran_jc/ common block is
!     updated with the values provided as argument.

SUBROUTINE g05cge(idum_in,iv_in,iy_in)
  IMPLICIT NONE

  INTEGER, PARAMETER :: ntab = 32

! Arguments to put into the memory of the random routines
  INTEGER :: idum_in,iv_in(ntab),iy_in

! Local variables
  INTEGER :: idum,iv(ntab),iy

  COMMON /ran_jc/ idum,iv,iy ! This common block constitutes the
                             ! 'memory' of the timeserie.
  INTEGER :: j               ! index

! Copy the values straight into the common
  DO j=1, ntab
    iv(j) = iv_in(j)
  END DO

  idum = idum_in
  iy = iy_in

  RETURN

END SUBROUTINE g05cge

