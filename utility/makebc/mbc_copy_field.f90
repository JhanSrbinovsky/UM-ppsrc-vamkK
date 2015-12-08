! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

SUBROUTINE mbc_copy_field(field_in, levels, rows, row_length, &
                          field_out)
IMPLICIT NONE

! Description : Copy data from a field to another using full sizes
!               Used for optimisation purposes

! Method : A straightforward loop/copy

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MakeBC

! Arguments

INTEGER, INTENT(IN) :: levels            ! Number of levels
INTEGER, INTENT(IN) :: rows              ! Number of rows
INTEGER, INTENT(IN) :: row_length        ! Length of row

REAL, INTENT(IN)  :: field_in(row_length, rows, levels)    ! Input data
REAL, INTENT(OUT) :: field_out(row_length, rows, levels)   ! Output data

! Local variables
INTEGER :: i,j,k     ! Loopers
!-----------------------------------------------------------------------

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                        &
!$OMP&            SHARED(levels, rows, row_length, field_in, field_out) &
!$OMP&            PRIVATE(i,j,k)
DO k = 1, levels
  DO j = 1, rows
    DO i = 1, row_length
      field_out(i,j,k) = field_in(i,j,k)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

RETURN
END SUBROUTINE mbc_copy_field
