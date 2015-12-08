! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ performs polar row averaging for theta fields

Module Rcf_average_polar_mod

!  Subroutine Rcf_average_polar_real    - for real data
!  Subroutine Rcf_average_polar_int     - for integer data
!  Subroutine Rcf_average_polar_logical - for logical data
!
! Description:
!   Polar row averaging perfomed on serial (non distributed) fields.
!   Horizontal interpolation done on a level-by-level basis.
!
! Method:
!   A generic interface is used to determine which averaging method
!   is utilised
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Implicit None

!-------------------------------------------------------------------
! Use a generic (overloaded) interface for different data-types
!-------------------------------------------------------------------
Interface Rcf_average_polar
  Module Procedure Rcf_average_polar_real, Rcf_average_polar_int, &
                   Rcf_average_polar_logical
End Interface

Contains

!------------------------------------------------------------------
! For real data
!------------------------------------------------------------------
Subroutine Rcf_average_polar_real(field, rows, row_length, GLOBAL, &
                                  average_required)

Implicit None

! Arguments
Integer, Intent(In)   :: rows
Integer, Intent(In)   :: row_length
Real, Intent(InOut)   :: field(row_length * rows)
Logical, Intent(In)   :: GLOBAL
Logical, Intent(Out)  :: average_required

! Local variables
Integer           :: i
Integer           :: last_row
Real              :: sum_start
Real              :: sum_end
Real              :: avg_start
Real              :: avg_end

last_row = (rows - 1) * row_length
sum_start = 0.0
sum_end   = 0.0
average_required = .FALSE.

If (GLOBAL) Then
  Do i = 1, row_length

    If ( field( i ) /= field( 1 ) .OR.                              &
         field( last_row + i ) /= field( last_row + 1 ) ) Then
      average_required = .TRUE.
    End If

    sum_start = sum_start + field( i )
    sum_end = sum_end + field( last_row + i )
  End Do

  If ( average_required ) Then
    avg_start = sum_start/Real(row_length)
    avg_end   = sum_end/Real(row_length)

    Do i = 1, row_length
      field(i) = avg_start
      field( last_row + i ) = avg_end
    End Do
  End If
End If

Return
End Subroutine Rcf_Average_polar_real

!------------------------------------------------------------------
! For integer data
!------------------------------------------------------------------
Subroutine Rcf_average_polar_int(field, rows, row_length, GLOBAL, &
                                 average_required)

Implicit None

! Arguments
Integer, Intent(In)        :: rows
Integer, Intent(In)        :: row_length
Integer, Intent(InOut)     :: field(row_length * rows)
Logical, Intent(In)        :: GLOBAL
Logical, Intent(Out)       :: average_required

! Local variables
Integer           :: i
Integer           :: last_row
Integer           :: sum_start
Integer           :: sum_end
Integer           :: avg_start
Integer           :: avg_end

last_row  = (rows - 1) * row_length
sum_start = 0
sum_end   = 0
average_required = .FALSE.

If (GLOBAL) Then
  Do i = 1, row_length

    If ( field( i ) /= field( 1 ) .OR.                              &
         field( last_row + i ) /= field( last_row + 1 ) ) Then
      average_required = .TRUE.
    End If

    sum_start = sum_start + field( i )
    sum_end = sum_end + field( last_row + i )
  End Do

  If ( average_required ) Then
    avg_start = Nint( Real(sum_start)/Real(row_length) )
    avg_end   = Nint( Real(sum_end)  /Real(row_length) )

    Do i = 1, row_length
      field(i) = avg_start
      field( last_row + i ) = avg_end
    End Do
  End If
End If

Return
End Subroutine Rcf_Average_polar_int

!------------------------------------------------------------------
! For logical data
!------------------------------------------------------------------
Subroutine Rcf_average_polar_logical(field, rows, row_length, GLOBAL, &
                                     average_required)

Implicit None

! Arguments
Integer, Intent(In)     :: rows
Integer, Intent(In)     :: row_length
Logical, Intent(InOut)  :: field(row_length * rows)
Logical, Intent(In)     :: GLOBAL
Logical, Intent(Out)    :: average_required

! Local variables
Integer           :: i
Integer           :: last_row
Integer           :: sum_start
Integer           :: sum_end
Logical           :: avg_start
Logical           :: avg_end

last_row = (rows - 1) * row_length
sum_start = 0
sum_end   = 0
average_required = .FALSE.

If (GLOBAL) Then
  Do i = 1, row_length

    If ( field( i ) .neqv. field( 1 ) .OR.                         &
         field( last_row + i ) .neqv. field( last_row + 1 ) ) Then
      average_required = .TRUE.
    End If

    If ( field(i) ) Then
      sum_start = sum_start + 1
    End If

    If ( field( last_row + i) ) Then
      sum_end = sum_end + 1
    End If
  End Do

  If (average_required) Then
    If ( sum_start >= row_length/2 ) Then
      avg_start = .TRUE.
    Else
      avg_start = .FALSE.
    End If

    If ( sum_end >= row_length/2 ) Then
      avg_end = .TRUE.
    Else
      avg_end = .FALSE.
    End If

    Do i = 1, row_length
      field(i) = avg_start
      field( last_row + i ) = avg_end
    End Do

  End If
End If

Return
End Subroutine Rcf_Average_polar_logical

End Module Rcf_Average_Polar_Mod
