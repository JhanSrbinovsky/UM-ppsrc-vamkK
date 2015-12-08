! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set variables in module cv_dependent_switch_mod
!
! Subroutine Interface:
SUBROUTINE cv_set_dependent_switches

! Modules used

USE cv_dependent_switch_mod, ONLY:                                         &
  l_var_entrain, l_new_det, l_const_ent, l_rh_dep,                         &
  sh_on, dp_on, cg_on, md_on,                                              &
  mdet_sh_on, mdet_dp_on, mdet_cg_on, mdet_md_on,                          &
  sh_ent_on,  dp_ent_on,  cg_ent_on,  md_ent_on,                           &
  sh_sdet_on, dp_sdet_on, cg_sdet_on, md_sdet_on,                          &
  sh_new_termc, dp_new_termc, cg_new_termc, md_new_termc,                  &
  cape_ts_w, sh_grey

USE cv_run_mod, ONLY:                                                      &
  icvdiag, adapt, termconv, cape_opt, cape_timescale,                      &
  i_convection_vn,i_convection_vn_4a, iconv_shallow, iconv_mid

USE turb_diff_mod, ONLY: l_subfilter_blend

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   This routine overrides default values held in cv_dependent_switch_mod
!   with values depending on convection namelist input. Called from readlsta
!   or scm_shell after convection namelist read in.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------
! Subroutine arguments - NONE at present

! Local variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CV_SET_DEPENDENT_SWITCHES',zhook_in,zhook_handle)

!-------------------------------------------------------------------------------
! 1.0 - Variables used in layer_cn  - 5A scheme only, no impact on 4A code
!-------------------------------------------------------------------------------

IF (icvdiag == 4 .OR. icvdiag == 5) THEN        ! diurnal cycle diagnosis
                           
  l_var_entrain = .true.      ! Use variable entrainment rate 
  l_new_det     = .true.      ! Use new detrainment relationship

END IF 

!-------------------------------------------------------------------------------
! 2.0 - Variables used in glue_conv and below controlling adaptive options
!-------------------------------------------------------------------------------
! Set flags for adaptive scheme depending on values passed in from umui
! flags all initialized to zero, so no changes needed if adapt  ==  0

IF (adapt  ==  1) THEN      !HadGEM1a

  md_on = 1
  dp_on = 1
  mdet_dp_on = 1
  mdet_md_on = 1

ELSE IF(adapt  ==  2) THEN  !convection test code

  md_on = 1
  dp_on = 1
  mdet_dp_on = 1
  mdet_md_on = 1
  md_ent_on=1
  dp_ent_on=1

ELSE IF(adapt  ==  3) THEN  !operational

  mdet_dp_on = 1
  dp_on = 1

ELSE IF(adapt  ==  4) THEN  ! Possible HadGEM3 as option 1 plus shallow

  md_on = 1
  dp_on = 1
  mdet_dp_on = 1
  mdet_md_on = 1
  sh_on = 1
! mdet_sh_on = 1   Mixing detrainment not set on as has no impact
  cg_on = 1         ! Only used in 5A code

ELSE IF(adapt  ==  5) THEN !adapt det + smoothed forced det for deep + mid

  md_on = 1
  dp_on = 1
  mdet_dp_on = 1
  mdet_md_on = 1
  md_sdet_on=1
  dp_sdet_on=1

ELSE IF(adapt  ==  6) THEN ! as 5 + shallow

  sh_on = 1
  md_on = 1
  dp_on = 1
! mdet_sh_on = 1  Mixing detrainment not set on as has no impact 
  mdet_md_on = 1
  mdet_dp_on = 1
  sh_sdet_on=1
  md_sdet_on=1
  dp_sdet_on=1

  cg_on = 1          ! Only used in 5A code
  cg_sdet_on=1       ! Only used in 5A code

ELSE IF(adapt  ==  7) THEN ! adapt det + smoothed forced det of theta, q, qcl and qcf for mid and deep

  md_on = 1
  dp_on = 1
  mdet_md_on = 1
  mdet_dp_on = 1
  md_sdet_on=2
  dp_sdet_on=2

ELSE IF(adapt  ==  8) THEN ! as 7 + shallow + congestus

  sh_on = 1
  md_on = 1
  dp_on = 1
! mdet_sh_on = 1  Mixing detrainment not set on as has no impact 
  mdet_md_on = 1
  mdet_dp_on = 1
  sh_sdet_on=2
  md_sdet_on=2
  dp_sdet_on=2

  cg_on = 1          ! Only used in 5A code
  cg_sdet_on=2       ! Only used in 5A code

END IF

! convection termination conditions

IF(termconv  ==  0) THEN
  
  md_new_termc=0
  dp_new_termc=0

ELSE IF(termconv  ==  1) THEN

  md_new_termc=1
  dp_new_termc=1
  IF (adapt == 4 .or. adapt == 6 .or. adapt == 8) THEN
    sh_new_termc=1          ! use new termination condition
    cg_new_termc=1          ! use new termination condition for congestus
  END IF

END IF 

IF (l_subfilter_blend) THEN
  sh_grey = 1
  iconv_shallow=1
  iconv_mid=0
END IF

IF (i_convection_vn == i_convection_vn_4a) THEN
! Set cape_ts_w at beginning of run - still required because stochastic
! physics uses cape_opt=3 and this uses cape_ts_w which remains fixed while
! cape_timescale alters. cape_ts_w is set to value of cape_timescale set in
! UMUI.

  IF (cape_opt == 3) THEN
    cape_ts_w = cape_timescale
  END IF

END IF

!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('CV_SET_DEPENDENT_SWITCHES',zhook_out,zhook_handle)
RETURN
END SUBROUTINE cv_set_dependent_switches
