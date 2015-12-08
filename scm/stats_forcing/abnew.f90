! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!=====================================================================
! SUBROUTINE ABnew
! Purpose:-           To calculate amplitude and mean of sinusoidal
!                     distribution for stats. Eqns. 10 and 11
!                     in SCM documentation.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
!=====================================================================

SUBROUTINE abnew(x1, x2, xa, xb, row_length, rows, n)

  IMPLICIT NONE

!---------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------

  INTEGER ::               &
    i, j, k                &! Loop counters
  , row_length, rows, n

  REAL ::                  &
    x1(row_length,rows,n)  &! In SD or mean of forcing variable
                            !    for max. of annual cycle (July)
  , x2(row_length,rows,n)  &! In SD or mean of forcing variable
                            !    for min. of annual cycle (Jan)
  , xa(row_length,rows,n)  &! Out Amplitude of seasonal variation
                            !     of forcing variable
  , xb(row_length,rows,n)   ! Out Mean of seasonal variation
                            !    of forcing variable

  DO i=1, row_length
    DO j=1, rows
      DO k=1, n
        xa(i,j,k) = (x1(i,j,k)-x2(i,j,k))/2.0
        xb(i,j,k) = (x1(i,j,k)+x2(i,j,k))/2.0
      END DO
    END DO
  END DO

  RETURN

END SUBROUTINE abnew

!=====================================================================

