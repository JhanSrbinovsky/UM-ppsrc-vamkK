! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Given surface pressure calculate ECMWF hybrid p-levels
!  (old 60 level version)

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_VARLOADER

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_calc_ecmwf_60levs(     &
 global_row_length,                       &  ! Length of model row
 global_rows,                             &  ! Number of model rows
 file_levels,                             &  ! Number of file levels
 data_pressure_surf_level,                &  ! Surface pressure
 data_pressure_ecmwf_levels,              &  ! Surface pressure
 debug)                                      ! Debug levels

USE nudging_control                       ! Standard nudging switches
USE nudging_ecmwf_60level_def! ECMWF pressure levels def module

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER, INTENT(IN) :: global_row_length      ! global row length
INTEGER, INTENT(IN) :: global_rows            ! global rows
INTEGER, INTENT(IN) :: file_levels            ! no. ECMWF levels

! Surface pressure from netcdf file
REAL, INTENT(IN):: data_pressure_surf_level                 &
(global_row_length, global_rows)

! Pressure on ECMWF levels
REAL, INTENT(OUT) :: data_pressure_ecmwf_levels               &
(global_row_length, global_rows, file_levels)

INTEGER, INTENT(IN) :: debug                  ! Debug flag

REAL                :: ak_int(1:ecmwf_nlevs)  ! A parameter
REAL                :: bk_int(1:ecmwf_nlevs)  ! B parameter
INTEGER             :: i, j, k                ! Loop counters

INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_CALC_ECMWF_60LEVS',zhook_in,zhook_handle)

!*******************************************************
! End of Header

! Standard Subroutine Entry Comment
IF(debug > 0) THEN
  WRITE(OUT, *)                                                       &
 ' NUDGING_CALC_ECMWF_LEVS: Entering Routine'
END IF

IF(file_levels /= ecmwf_nlevs) THEN
  WRITE(OUT, *)                                                       &
   '  NUDGING_CALC_ECMWF_LEVS:',                                      &
   ' Warning: Number of levels in NetCDF',                            &
  'File does not equal ECMWF Levels'
   RETURN
END IF

DO i=1, ecmwf_nlevs

! Convert the levels to the same order as UM
! (i.e  start at the bottom of the atmosphere)
  j=(ecmwf_nlevs+1) - i

! Average the pressure levels
  ak_int(i) = 0.5*(ak(j-1) +ak(j))
  bk_int(i) = 0.5*(bk(j-1) +bk(j))

END DO

! Determine the pressure levels for ECMWF
DO k=1, file_levels
  data_pressure_ecmwf_levels(:,:,k) =                                 &
  ak_int(k) + bk_int(k)* data_pressure_surf_level(:,:)
END DO

! Standard Subroutine Exit Comment
IF(debug > 0) THEN
  WRITE(OUT, *)                                                       &
 ' NUDGING_CALC_ECMWF_LEVS: Leaving Routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_CALC_ECMWF_60LEVS',zhook_out,zhook_handle)

RETURN

END SUBROUTINE nudging_calc_ecmwf_60levs

