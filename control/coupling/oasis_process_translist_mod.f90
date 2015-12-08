MODULE OASIS_process_translist_mod
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************

IMPLICIT NONE

CONTAINS

SUBROUTINE OASIS_process_translist()

!
! Description:
! Process transient list to see what special conditions need to be catered for. 
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Coupling
! 
!==================================================================

USE OASIS_atm_data_mod, ONLY: transient_in, max_transients
USE Field_Types, ONLY : fld_type_u, fld_type_v
USE UM_ParVars, ONLY : glsize, at_extremity, PNorth
USE domain_params, ONLY : mt_global
USE um_input_control_mod, ONLY: model_domain
USE um_types
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER :: i_number

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('OASIS_process_translist',zhook_in,zhook_handle)

DO i_number = 1, max_transients 

  ! Polar meaning only applies to ND, and then only if it's a global 
  ! model, and only for fields on T points and even then we only 
  ! need to set the flag on northern-most PEs.
  IF ( (transient_in(i_number)%grid=="T").AND.             &
       (glsize(2,fld_type_u) > glsize(2,fld_type_v)).AND.  &
       (model_domain == mt_global).AND.                    &
       (at_extremity(PNorth)) ) THEN
       transient_in(:)%polar_mean = .TRUE.
  END IF

END DO

IF (lhook) CALL dr_hook('OASIS_process_translist',zhook_out,zhook_handle)

END SUBROUTINE OASIS_process_translist

END MODULE OASIS_process_translist_mod
