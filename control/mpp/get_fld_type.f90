! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Function Interface
INTEGER FUNCTION get_fld_type (grid_type_code)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParParams
USE cppxref_mod   ! ppx_XXX items - most of module used
IMPLICIT NONE


! Description:
! Takes a STASH grid type code, and returns which type of
! grid this is - mass or wind grid.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

INTEGER, INTENT(IN)  ::  grid_type_code     ! IN : STASH grid type code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('GET_FLD_TYPE',zhook_in,zhook_handle)
IF ( (grid_type_code  ==  ppx_atm_tall) .OR.                      &
     (grid_type_code  ==  ppx_atm_tland) .OR.                     &
     (grid_type_code  ==  ppx_atm_tsea) .OR.                      &
     (grid_type_code  ==  ppx_atm_tzonal) .OR.                    &
     (grid_type_code  ==  ppx_atm_tmerid) .OR.                    &
     (grid_type_code  ==  ppx_atm_compressed) .OR.                &
     (grid_type_code  ==  ppx_atm_ozone) .OR.                     &
     (grid_type_code  ==  ppx_ocn_tall)    .OR.                   &
     (grid_type_code  ==  ppx_ocn_tfield)  .OR.                   &
     (grid_type_code  ==  ppx_ocn_tzonal)  .OR.                   &
     (grid_type_code  ==  ppx_atm_lbc_theta)  .OR.                &
     (grid_type_code  ==  ppx_atm_lbc_orog)  .OR.                 &
     (grid_type_code  ==  ppx_ocn_lbc_theta)  .OR.                &
     (grid_type_code  ==  ppx_wam_all) .OR.                       &
     (grid_type_code  ==  ppx_ocn_tmerid) ) THEN
  get_fld_type=fld_type_p
ELSE IF                                                           &
   ( (grid_type_code  ==  ppx_atm_cuall) .OR.                     &
     (grid_type_code  ==  ppx_ocn_uall) .OR.                      &
     (grid_type_code  ==  ppx_ocn_cuall) .OR.                     &
     (grid_type_code  ==  ppx_ocn_ufield) .OR.                    &
     (grid_type_code  ==  ppx_ocn_uzonal) .OR.                    &
     (grid_type_code  ==  ppx_atm_lbc_u)  .OR.                    &
     (grid_type_code  ==  ppx_ocn_lbc_u)  .OR.                    &
     (grid_type_code  ==  ppx_ocn_umerid) ) THEN
  get_fld_type=fld_type_u
ELSE IF                                                           &
   ( (grid_type_code  ==  ppx_atm_cvall) .OR.                     &
     (grid_type_code  ==  ppx_atm_lbc_v)  .OR.                    &
     (grid_type_code  ==  ppx_atm_uall ) .OR.                     &
     (grid_type_code  ==  ppx_ocn_cvall) ) THEN
  get_fld_type=fld_type_v
ELSE IF                                                           &
   ( (grid_type_code  ==  ppx_atm_uall) .OR.                      &
     (grid_type_code  ==  ppx_atm_uland) .OR.                     &
     (grid_type_code  ==  ppx_atm_usea) .OR.                      &
     (grid_type_code  ==  ppx_atm_uzonal) .OR.                    &
     (grid_type_code  ==  ppx_atm_umerid) .OR.                    &
     (grid_type_code  ==  ppx_ocn_uall) .OR.                      &
     (grid_type_code  ==  ppx_ocn_ufield) .OR.                    &
     (grid_type_code  ==  ppx_ocn_uzonal) .OR.                    &
     (grid_type_code  ==  ppx_ocn_umerid) ) THEN
  get_fld_type=fld_type_v  ! This is actually the B U/V (velocity)
                           ! grid, but it has the same sizes as
                           ! the C V grid.
ELSE IF                                                           &
     (grid_type_code  ==  ppx_wam_sea) THEN
  get_fld_type=fld_type_comp_wave
ELSE IF                                                           &
     (grid_type_code  ==  ppx_wam_rim) THEN
  get_fld_type=fld_type_rim_wave
ELSE IF                                                           &
   ( (grid_type_code  ==  ppx_atm_river) ) THEN
  get_fld_type=fld_type_r  ! River routing grid

ELSE
  get_fld_type=fld_type_unknown
END IF

IF (lhook) CALL dr_hook('GET_FLD_TYPE',zhook_out,zhook_handle)
RETURN

END FUNCTION get_fld_type

