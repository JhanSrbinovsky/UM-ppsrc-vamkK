! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!=====================================================================
! SUBROUTINE ACInit
! PURPOSE:-           To calculate mean and SD of a random variable
!                     eqns 6 and 7 in SCM doc.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Single Column Model
!
!=====================================================================

SUBROUTINE acinit (xbar, xsd, a, cbar, csd, cor, n, row_length, rows)

!---------------------------------------------------------------------
! Arguments
!---------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER ::    &
    n           &! In no of model_levels, rows and columns
  , row_length  &
  , rows

  REAL ::                      &
    a(row_length,rows,n-1)     &! Out term a of eqn. 2.22
  , cbar(row_length,rows,n-1)  &! Out Mean of random variable C
  , cor(row_length,rows)       &! In Vertical correlation coefficient
  , csd(row_length,rows,n-1)   &! Out SD of random variable C
  , xbar(row_length,rows,n)    &! In Mean of forcing variable
  , xsd(row_length,rows,n)      ! In SD of forcing variable
!---------------------------------------------------------------------
!     Local variables
!---------------------------------------------------------------------

  INTEGER ::    &
    i, j, k      ! Loop counters

  DO  k=1 ,n-1
    DO j=1, rows
      DO i=1, row_length
        a(i,j,k) = cor(i,j) * xsd(i,j,k+1) / xsd(i,j,k)
        cbar(i,j,k) = xbar(i,j,k+1) - a(i,j,k) * xbar(i,j,k)
        csd(i,j,k) = SQRT(1.0-cor(i,j)*cor(i,j)) * xsd(i,j,k+1)
      END DO
    END DO
  END DO

  RETURN

END SUBROUTINE acinit

