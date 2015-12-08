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
!------------------------------------------------------------------
!     Function inspired from the 'Numerical Recipies in fortran 77'
!     well known book. Generates a random number uniformly distributed
!     over [0,1].
!     It should not be called direct externally from outside this
!     deck. Instead it is accessed through the G05xx
!     routines/functions when they need a random number.
!
!     Arguments : nil
!     Returned value : a random number uniformly distributed on [0,1].
!     Side effects : the status of the /ran_jc/ common block
!                    is updated at each call.

FUNCTION ran1_jc()

  IMPLICIT NONE

  INTEGER :: idum

  INTEGER, PARAMETER :: &
    ia   = 16807        &
  , im   = 2147483647   &
  , iq   = 127773       &
  , ir   = 2836         &
  , ntab = 32           &
  , ndiv = 1 + (im-1)/ntab

  REAL, PARAMETER :: &
    am   = 1.0/im    &
  , eps  = 1.2e-7    &
  , rnmx = 1.0-eps

  REAL :: ran1_jc

! *Minimal* random number generator of Park and Miller with
! Bays-Durham shuffle and added safeguards. Returns a uniform
! random deviate between 0.0  and 1.0 (exclusive of the endpoint values).
!
! Call with idum a negative integer to initialize; thereafter, do
! not alter idum between successive deviates in a sequence.
! RNMX should approximate the largest Floatingvalue that is less than 1.

  INTEGER :: j, k, iv(ntab), iy

  COMMON /ran_jc/ idum,iv,iy ! This common block constitutes the
                             ! 'memory' of the timeserie ;
                             !  it is manipulated by other routines
                             !  for intializing, saving, restoring
                             !  the timeseries.

  IF (idum <= 0.OR.iy == 0) THEN ! Initialize.
    idum = MAX(-idum,1)      ! Be sure to prevent idum = 0

    DO j=ntab+8, 1, -1       ! Load the shuffle table (after 8 warm-ups).
      k = idum/iq
      idum = ia*(idum-k*iq)-ir*k

      IF (idum <  0) THEN
        idum = idum + im
      END IF

      IF (j <= ntab)  THEN
        iv(j) = idum
      END IF

    END DO
     iy = iv(1)
  END IF

  k = idum/IQ                 ! Start here when not initializing.
  idum = ia*(idum-k*iq)-ir*k  ! Compute idum=mod(IA*idum,IM) without
                              ! overflows by Schrage's method.
  IF (idum <  0) THEN
    idum = idum+im
  END IF

  j  = 1 + iy/ndiv            ! Will be in the range 1:NTAB.
  iy = iv(j)                  ! Output previously stored value and
                              ! refill the shuffle table
  iv(j) = idum
  ran1_jc = MIN(am*iy,rnmx)   ! Because users don't expect endpoint
                              ! values.

  RETURN

END FUNCTION ran1_jc

!----------------------------------------------------------------

