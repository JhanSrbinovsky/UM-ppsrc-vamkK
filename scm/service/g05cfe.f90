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
!     Save the state of the random number generator
!
!     Arguments : they are all output value. It is the
!     task of the user to store those values somewhere
!     until s/he wants to reinitialise to an identical
!     serie of random numbers with G05CGE, called
!     with the same arguments.
!     Side effect : nil

SUBROUTINE g05cfe(idum_out, iv_out, iy_out)
  IMPLICIT NONE

  INTEGER, PARAMETER :: ntab = 32

! Arguments to extract from the memory of the random routines
  INTEGER :: idum_out,iv_out(ntab),iy_out

! Local variables
  INTEGER :: idum,iv(ntab),iy
  COMMON /ran_jc/ idum,iv,iy ! This common block constitutes the
                             ! 'memory' of the timeserie.
  INTEGER :: j


! Copy the values straight from the common into the _out variables.
  DO j=1, ntab
    iv_out(j) = iv(j)
  END DO

  idum_out = idum
  iy_out = iy

  RETURN

END SUBROUTINE g05cfe

