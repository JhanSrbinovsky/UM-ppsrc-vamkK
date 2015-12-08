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
!     Pseudo-random real numbers, Gaussian distribution of
!     mean m, standard deviation s.
!
!     Also inspired from the 'Numerical Recipies in fortran 77' book.
!
!     Arguments : Input (m,s) ; are not destroyed after the call.
!     calls RAN_JC() in order to access a uniform random number
!     generator over [0,1].
!     Side effect : changes the values of the /ran_jc/ common
!     block
FUNCTION g05dde(m,s)

  IMPLICIT NONE

! Function type
  REAL :: g05dde

! Arguments
  REAL :: m,s         ! (mean, standard deviation) of the distribution.
  REAL :: ran1_jc

!     Returns a normally distributed deviate with m mean and s
!     variance. Uusing ran1_jc() as the source of uniform deviates.

! Local variables
  INTEGER :: iset
  REAL :: fac,gset,rsq,v1,v2
  SAVE iset, gset
  DATA iset/0/

  IF (iset == 0) THEN                    ! We don't have an extra deviate
! DEPENDS ON: ran1_jc
 1  v1 = 2.0*ran1_jc()-1.0                 ! handy, so pick two uniform
! DEPENDS ON: ran1_jc                      ! numbers
    v2 = 2.0*ran1_jc()-1.0                 ! in the square extending from
    rsq = v1**2+v2**2                      ! -1 to +1 in each direction see
    IF (rsq >= 1.0 .OR. rsq == 0.0) GOTO 1 ! if they are in the unit circle.
                                           ! and if they are not, try again.

    fac = SQRT(-2.0*LOG(rsq)/rsq)  ! Now make the Box-Muller
                                   ! transformation to get two normal
                                   ! deviates. Return one and save the
                                   ! other for next time.
    gset = v1*fac
    g05dde = m+s*(v2*fac)
    iset = 1                ! Set flag.
  ELSE                      ! We have an extra deviate handy.
    g05dde = m+s*gset       ! so return it.
    iset = 0                ! and unset the flag.
  ENDIF

  RETURN

END FUNCTION g05dde

