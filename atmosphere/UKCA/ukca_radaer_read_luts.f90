! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!+ Coordinates the input of UKCA_RADAER look-up tables for each
!+ spectra and aerosol modes.
!
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
SUBROUTINE ukca_radaer_read_luts(icode, cmessage)

USE parkind1, ONLY: jpim, jprb 
USE yomhook,  ONLY: lhook, dr_hook
USE ukca_radaer_lut

IMPLICIT NONE

!
! Arguments
!

! error indicator (0 is OK, <0 is KO)
INTEGER :: icode

! error message
CHARACTER (len=80) :: cmessage

! Local variables

! Standard logical unit numbers for the 4 look-up tables
INTEGER, PARAMETER :: unit_sw_ac = 166, &
                      unit_lw_ac = 167, &
                      unit_sw_cr = 168, &
                      unit_lw_cr = 169

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('UKCA_RADAER_READ_LUTS', zhook_in, zhook_handle)

icode = 0

! DEPENDS ON: ukca_radaer_lut_in
CALL ukca_radaer_lut_in(icode, cmessage, unit_sw_ac, &
  ukca_lut(ip_ukca_lut_accum, ip_ukca_lut_sw))
IF (icode /= 0) THEN
  IF (lhook) CALL dr_hook('UKCA_RADAER_READ_LUTS', zhook_out, zhook_handle)
  RETURN
END IF

! DEPENDS ON: ukca_radaer_lut_in
call ukca_radaer_lut_in(icode, cmessage, unit_sw_cr, &
  ukca_lut(ip_ukca_lut_coarse, ip_ukca_lut_sw))
IF (icode /= 0) THEN
  IF (lhook) CALL dr_hook('UKCA_RADAER_READ_LUTS', zhook_out, zhook_handle)
  RETURN
END IF


! DEPENDS ON: ukca_radaer_lut_in
CALL ukca_radaer_lut_in(icode, cmessage, unit_lw_ac, &
  ukca_lut(ip_ukca_lut_accum, ip_ukca_lut_lw))
IF (icode /= 0) THEN
  IF (lhook) CALL dr_hook('UKCA_RADAER_READ_LUTS', zhook_out, zhook_handle)
  RETURN
END IF


! DEPENDS ON: ukca_radaer_lut_in
CALL ukca_radaer_lut_in(icode, cmessage, unit_lw_cr, &
  ukca_lut(ip_ukca_lut_coarse, ip_ukca_lut_lw))
IF (icode /= 0) THEN
  IF (lhook) CALL dr_hook('UKCA_RADAER_READ_LUTS', zhook_out, zhook_handle)
  RETURN
END IF

IF (lhook) CALL dr_hook('UKCA_RADAER_READ_LUTS', zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_read_luts
