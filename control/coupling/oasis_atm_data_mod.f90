! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************

MODULE oasis_atm_data_mod
  ! Description:
  ! Useful data and arrays for use with OASIS coupling.
  ! This is expected to be applicable to all OASIS3
  ! versions including OASIS3-MCT.
  !
  ! Code Owner: See Unified Model Code Owners HTML page
  ! This file belongs in section: Coupling
  !
  !=====================================================================

  USE um_types
  USE field_types, ONLY : fld_type_u, fld_type_v, fld_type_p

  IMPLICIT NONE

  ! Things needed for general compilation of all models
  INTEGER (KIND=integer64) :: oasis_couple_ts ! Coupling timestep length


END MODULE oasis_atm_data_mod

