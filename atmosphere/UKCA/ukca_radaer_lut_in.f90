! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!+ Reads in a UKCA_RADAER look-up table (namelist format) and copy 
!+ its contents in the structure given in argument.
!
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
SUBROUTINE ukca_radaer_lut_in(icode, cmessage, funit, lut)

USE parkind1, ONLY: jpim, jprb 
USE yomhook,  ONLY: lhook, dr_hook

USE ukca_radaer_tlut_mod

IMPLICIT NONE

!
! Arguments
!

! error indicator (0 is OK, <O is KO)
INTEGER :: icode

! error message
CHARACTER (LEN=80) :: cmessage

! look-up table logical unit number
INTEGER :: funit

! Structure representing a look-up table
TYPE (ukca_radaer_tlut) :: lut

!
! Local variables
!
CHARACTER (len=80) :: filename

!
! Local copy of the contents of a UKCA look-up table.
! Here, fixed array sizes are needed. 
! Note that array indexing starts at 0 in the namelist.
INTEGER, PARAMETER :: NPD_X  = 50, &
                      NPD_NR = 50, &
                      NPD_NI = 50

REAL :: stdev
INTEGER :: n_x
INTEGER :: n_nr
INTEGER :: n_ni
REAL :: x_min
REAL :: x_max
REAL :: nr_min
REAL :: nr_max
REAL :: ni_min
REAL :: ni_max
REAL, DIMENSION(0:NPD_X,0:NPD_NI,0:NPD_NR) :: ukca_absorption
REAL, DIMENSION(0:NPD_X,0:NPD_NI,0:NPD_NR) :: ukca_scattering
REAL, DIMENSION(0:NPD_X,0:NPD_NI,0:NPD_NR) :: ukca_asymmetry
REAL, DIMENSION(0:NPD_X) :: volume_fraction

NAMELIST /UKCANML/                                                      &
          stdev,                                                        &
          n_x, n_nr, n_ni,                                              &
          x_min, x_max,                                                 &
          nr_min, nr_max,                                               &
          ni_min, ni_max,                                               &
          ukca_absorption, ukca_scattering, ukca_asymmetry,             &
          volume_fraction

INTEGER :: ios
INTEGER :: i
INTEGER :: j
INTEGER :: k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('UKCA_RADAER_LUT_IN', zhook_in, zhook_handle)

!
! Get the file name, open the file and perform a namelist read.
!
CALL get_file(funit, filename, 80, ios)
IF (ios /= 0) THEN
  icode = -1
  cmessage = 'Error reading name of UKCA look-up table file.'
  IF (lhook) CALL dr_hook('UKCA_RADAER_LUT_IN', zhook_out, zhook_handle)
  RETURN
END IF

OPEN(unit=funit, file=filename, iostat=ios)
IF (ios /= 0) THEN
  icode = -1
  cmessage='Error opening ' // filename
  IF (lhook) CALL dr_hook('UKCA_RADAER_LUT_IN', zhook_out, zhook_handle)
  RETURN
END IF
READ(funit, UKCANML)
CLOSE(funit)

!
! Check that actual array dimensions do not exceed
! those that were fixed at compile time.
!
IF (n_x > NPD_X) THEN
  icode = -1
  cmessage='Look-up table dimension exceeds built-in limit (X).'
  IF (lhook) CALL dr_hook('UKCA_RADAER_LUT_IN', zhook_out, zhook_handle)
  RETURN
END IF

IF (n_nr > NPD_NR) THEN
  icode = -1
  cmessage='Look-up table dimension exceeds built-in limit (NR).'
  IF (lhook) CALL dr_hook('UKCA_RADAER_LUT_IN', zhook_out, zhook_handle)
  RETURN
END IF

IF (n_ni > NPD_NI) THEN
  icode = -1
  cmessage='Look-up table dimension exceeds built-in limit (NI).'
  IF (lhook) CALL dr_hook('UKCA_RADAER_LUT_IN', zhook_out, zhook_handle)
  RETURN
END IF

! Copy the look-up table into the structure passed as an argument
lut%stdev = stdev

!
! Namelist arrays indices start at 0. For the sake of ease of use, 
! online arrays start at 1 -> simply shift the content by 1 and add 
! 1 to each dimension.
! All calculations within UKCA_RADAER expect array indexing starting 
! at 1.
!
n_x = n_x + 1
n_nr = n_nr + 1
n_ni = n_ni + 1

lut%n_x  = n_x
lut%n_nr = n_nr
lut%n_ni = n_ni

lut%x_min  = x_min
lut%x_max  = x_max
lut%nr_min = nr_min
lut%nr_max = nr_max
lut%incr_nr = (nr_max - nr_min) / float(n_nr-1)
lut%ni_min = ni_min
lut%ni_max = ni_max
lut%incr_ni = (ni_max - ni_min) / float(n_ni-1)

! Allocate the dynamic arrays
ALLOCATE(lut%ukca_absorption(n_x, n_ni, n_nr))
ALLOCATE(lut%ukca_scattering(n_x, n_ni, n_nr))
ALLOCATE(lut%ukca_asymmetry(n_x, n_ni, n_nr))
ALLOCATE(lut%volume_fraction(n_x))

DO k = 1, n_nr
  DO j = 1, n_ni
    DO i = 1, n_x
      lut%ukca_absorption(i, j, k) = ukca_absorption(i-1, j-1, k-1)
      lut%ukca_scattering(i, j, k) = ukca_scattering(i-1, j-1, k-1)
      lut%ukca_asymmetry(i, j, k)  = ukca_asymmetry(i-1, j-1, k-1)
    END DO
  END DO
END DO

DO i = 1, n_x
  lut%volume_fraction(i) = volume_fraction(i-1)
END DO

IF (lhook) CALL dr_hook('UKCA_RADAER_LUT_IN', zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_lut_in
