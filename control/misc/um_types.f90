! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C70

! Parameters for 32 and 64 bit kinds

MODULE um_types
  IMPLICIT NONE
  ! Precision and range for 64 bit real
  INTEGER, PARAMETER :: prec64  = 15
  INTEGER, PARAMETER :: range64 = 307

  ! Precision and range for 32 bit real
  INTEGER, PARAMETER :: prec32  = 6
  INTEGER, PARAMETER :: range32 = 37

  ! Range for integers
  INTEGER, PARAMETER :: irange64=15
  INTEGER, PARAMETER :: irange32=9

  ! Kind for 64 bit real
  INTEGER, PARAMETER :: real64  = selected_real_kind(prec64,range64)
  ! Kind for 32 bit real
  INTEGER, PARAMETER :: real32  = selected_real_kind(prec32,range32)
  ! Kind for 64 bit integer
  INTEGER, PARAMETER :: integer64 = selected_int_kind(irange64)
  ! Kind for 32 bit integer
  INTEGER, PARAMETER :: integer32 = selected_int_kind(irange32)

  ! Kinds for 64 and 32 bit logicals. Note that there is no
  ! "selected_logical_kind", but using the equivalent integer kind is a
  ! workaround that works on every platform we have tested.
  INTEGER, PARAMETER :: logical64 = integer64
  INTEGER, PARAMETER :: logical32 = integer32

  CONTAINS
    
! Discover the size of a fortran default real
    FUNCTION umFortranRealSize() RESULT (r)
      IMPLICIT NONE
      INTEGER :: r
      REAL    :: a(2)
      CALL um_addr_diff(a(1),a(2),r)
    END FUNCTION umFortranRealSize

! Discover the size of a fortran default integer
    FUNCTION umFortranIntegerSize() RESULT (r)
      IMPLICIT NONE
      INTEGER :: r
      INTEGER :: a(2)
      CALL um_addr_diff(a(1),a(2),r)
    END FUNCTION umFortranIntegerSize

! Discover the size of a pointer (we assume c and fortran are the same)
    FUNCTION umFortranPointerSize() RESULT (r)
      IMPLICIT NONE
      INTEGER, EXTERNAL :: um_addr_size
      INTEGER :: r
      r=um_addr_size()
    END FUNCTION umFortranPointerSize

END MODULE um_types

