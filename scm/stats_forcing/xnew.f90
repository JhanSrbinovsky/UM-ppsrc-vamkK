! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!=====================================================================
! Subroutine Xnew
! Purpose:-           To calculate mean or SD of random variable
!                     at daynumber relative to winter solstice
!                     (eqn. 12 in SCM doc.)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
!=====================================================================

SUBROUTINE xnew(x, xa, xb, row_length, rows, nlevs, xt)

!---------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER ::   &
    row_length &! In dimension of arrays
  , rows       &
  , nlevs       ! model levels

  REAL ::                     &
    x(row_length,rows,nlevs)  &! Out Mean or SD of forcing variable
                               !     at day relative to winter
                               !     solstice
  , xa(row_length,rows,nlevs) &! In  Amplitude of seasonal variation
                               !     of forcing variable.
  , xb(row_length,rows,nlevs) &! In  Mean of seasonal variation
                               !     of forcing variable.
  , xt                         ! In  Sin of argument

!---------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------
  INTEGER ::   &
    i,j,k       ! Loop counters

  DO k=1, nlevs
    DO j=1, rows
      DO i=1, row_length
        x(i,j,k) = xa(i,j,k) * xt + xb(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE xnew

!=====================================================================

