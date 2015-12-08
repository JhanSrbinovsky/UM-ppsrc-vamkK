! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!   This module defines the elements required for radiative forcing

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90

!- End of header

MODULE coradoca

USE max_calls
USE rad_input_mod, ONLY: l_forcing

IMPLICIT NONE
INTEGER, PARAMETER :: c2c_size=npd_swcall-1
!     Whether to change each quantity in diagnostic calls to radiation.
!     The arrays will normally be of size 1, but dimensioning them as
!     arrays will greatly simplify any changes to have more than 1 call.

! Defaults for things to change in diagnostic 2nd call to radiation.
LOGICAL ::  c2c_all(c2c_size) = .FALSE.
LOGICAL ::  c2c_wmg(c2c_size) = .FALSE.
LOGICAL ::  c2c_o2(c2c_size) = .FALSE.
LOGICAL ::  c2c_o3(c2c_size) = .FALSE.
LOGICAL ::  c2c_co2(c2c_size) = .FALSE.
LOGICAL ::  c2c_n2o(c2c_size) = .FALSE.
LOGICAL ::  c2c_ch4(c2c_size) = .FALSE.
LOGICAL ::  c2c_cfc11(c2c_size) = .FALSE.
LOGICAL ::  c2c_cfc12(c2c_size) = .FALSE.
LOGICAL ::  c2c_c113(c2c_size) = .FALSE.
LOGICAL ::  c2c_hcfc22(c2c_size) = .FALSE.
LOGICAL ::  c2c_hfc125(c2c_size) = .FALSE.
LOGICAL ::  c2c_hfc134(c2c_size) = .FALSE.
LOGICAL ::  c2c_aerosol(c2c_size) = .FALSE.
LOGICAL ::  c2c_sulpc_d(c2c_size) = .FALSE.
LOGICAL ::  c2c_seas_d(c2c_size) = .FALSE.
LOGICAL ::  c2c_soot_d(c2c_size) = .FALSE.
LOGICAL ::  c2c_bmb_d(c2c_size) = .FALSE.
LOGICAL ::  c2c_ocff_d(c2c_size) = .FALSE.
LOGICAL ::  c2c_nitr_d(c2c_size) = .FALSE.
LOGICAL ::  c2c_dust_d(c2c_size) = .FALSE.
LOGICAL ::  c2c_biog_d(c2c_size) = .FALSE.
LOGICAL ::  c2c_ukca_d(c2c_size) = .FALSE.
LOGICAL ::  c2c_land_s(c2c_size) = .FALSE.


NAMELIST / radfcdia / c2c_o2, c2c_o3, c2c_co2, c2c_n2o, c2c_ch4,  &
  c2c_cfc11, c2c_cfc12, c2c_c113, c2c_hcfc22, c2c_hfc125,         &
  c2c_hfc134, c2c_aerosol, c2c_sulpc_d, c2c_seas_d, c2c_soot_d,   &
  c2c_bmb_d, c2c_ocff_d, c2c_land_s, c2c_all,                     &
  c2c_wmg, c2c_nitr_d, c2c_dust_d, c2c_biog_d, c2c_ukca_d


! ----------------------------------------------------------------
CONTAINS

! Subroutine to set the input values of the control structure.

SUBROUTINE coradoca_defaults

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER :: j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('coradoca_defaults',zhook_in,zhook_handle)

DO j=1, c2c_size
  IF ( .NOT. ( c2c_o2(j) .OR. c2c_o3(j) .OR. c2c_co2(j) .OR.      &
     c2c_n2o(j) .OR. c2c_ch4(j) .OR. c2c_cfc11(j) .OR.            &
     c2c_cfc12(j) .OR. c2c_c113(j) .OR. c2c_hcfc22(j) .OR.        &
     c2c_hfc125(j) .OR. c2c_hfc134(j) .OR. c2c_aerosol(j) .OR.    &
     c2c_sulpc_d(j) .OR. c2c_seas_d(j) .OR. c2c_soot_d(j) .OR.    &
     c2c_bmb_d(j) .OR. c2c_ocff_d(j) .OR.                         &
     c2c_land_s(j) .OR. c2c_all(j) .OR.                           &
     c2c_wmg(j) .OR. c2c_nitr_d(j) .OR. c2c_dust_d(j) .OR.        &
     c2c_biog_d(j) .OR. c2c_ukca_d(j)) .AND. l_forcing ) THEN
!          ErrorStatus = 614
    WRITE(6,'(a,a)') 'warning: forcing switched on but nothing,', &
                     'set to change IN diagnostic CALL!'

  END IF
  IF ( c2c_all(j) ) THEN
! C_O3(J) = .TRUE.  ! Commented out till we have the ancillary
! Best made conditional on the alternate ozone file being specified
! if we don't hardwire the latter's name (in which case it'd always
! be safe to read from that file, & add nothing much to the cost of
! a 2nd call which is going to be made anyway).
    c2c_wmg(j) = .TRUE.
    c2c_aerosol(j) = .TRUE.
    c2c_land_s(j) = .TRUE.
! This is used for surface albedo forcing
  END IF
  IF ( c2c_wmg(j) ) THEN
    c2c_o2(j) = .TRUE.
    c2c_co2(j) = .TRUE.
    c2c_n2o(j) = .TRUE.
    c2c_ch4(j) = .TRUE.
    c2c_cfc11(j) = .TRUE.
    c2c_cfc12(j) = .TRUE.
    c2c_c113(j) = .TRUE.
    c2c_hcfc22(j) = .TRUE.
    c2c_hfc125(j) = .TRUE.
    c2c_hfc134(j) = .TRUE.
  END IF
  IF ( c2c_aerosol(j) ) THEN
    c2c_sulpc_d(j) = .TRUE.
    c2c_seas_d(j) = .TRUE.
    c2c_soot_d(j) = .TRUE.
    c2c_bmb_d(j) = .TRUE.
    c2c_ocff_d(j) = .TRUE.
    c2c_nitr_d(j) = .TRUE.
    c2c_dust_d(j) = .TRUE.
    c2c_biog_d(j) = .TRUE.
    c2c_ukca_d(j) = .TRUE.
  END IF
END DO

IF (lhook) CALL dr_hook('coradoca_defaults',zhook_out,zhook_handle)

END SUBROUTINE coradoca_defaults



END MODULE coradoca
