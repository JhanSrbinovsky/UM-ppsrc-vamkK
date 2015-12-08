! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads the NLSIZES namelist
!
MODULE rcf_readnl_nlsizes_mod

IMPLICIT NONE

! Description:
!  rcf_readnl_nlsizes reads in the NLSIZES namelist, containing information
!  about the various grid levels/sizes to be used.
!
! Method:
!  Data read in via unit number passed as argument.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 v8.5 programming standards.

CONTAINS

SUBROUTINE rcf_readnl_nlsizes (nft)

! NLSIZES namelist definition
USE nlsizes_namelist_mod   ! All

USE rimtypes, ONLY:                                                      &
  rima_type_norm

USE rcf_lsm_mod, ONLY:                                                   &
  glob_land_out

USE ereport_mod, ONLY:                                                   &
  ereport

USE PrintStatus_mod, ONLY:                                               &
  PrintStatus, PrStatus_Oper

USE UM_parvars, ONLY:                                                    &
  mype

USE rcf_grid_type_mod, ONLY:                                             &
  output_grid

USE missing_data_mod, ONLY:                                              &
  IMDI

USE calc_ntiles_mod, ONLY :                                              &
  calc_ntiles

USE switches, ONLY:                                                      &
  l_aggregate

USE nstypes, ONLY:                                                       &
  npft, nnvg

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)  :: nft  ! unit number

! Local variables
INTEGER              :: ErrorStatus

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'rcf_readnl_nlsizes'
CHARACTER (LEN=80)           :: cmessage


! Set defaults
global_row_length = 0
global_rows       = 0
model_levels      = 0
wet_levels        = 0
st_levels         = 0
bl_levels         = 0
ozone_levels      = 0
tr_ukca           = 0
river_rows        = 180   ! 1 degree default
river_row_length  = 360   ! 1 degree default
rimwidtha(rima_type_norm) = 0
pp_len_inthd      = 46
pp_len_realhd     = 38

! Sea ice improvements still in development, so set to old scheme here.
nice_use          = 1

! Previously set to IMDI in the reconfiguration
a_len2_rowdepc    = IMDI
a_len2_coldepc    = IMDI
a_len_inthd       = IMDI
a_len_realhd      = IMDI
a_len2_levdepc    = IMDI
a_len2_flddepc    = IMDI
a_len_extcnst     = IMDI

! No prior default in the reconfiguration
land_field        = IMDI
nice              = IMDI

! Namelist variables not used by the reconfiguration
cloud_levels      = IMDI
tpps_ozone_levels = IMDI
tr_lbc_vars       = IMDI
tr_lbc_ukca       = IMDI
a_len_cfi1        = IMDI
a_len_cfi2        = IMDI
a_len_cfi3        = IMDI
nancil_lookupsa   = IMDI
n_intf_a          = IMDI
nrim_timesa       = IMDI

! Read the namelists
READ(UNIT = nft, NML = NLSIZES, IOSTAT = ErrorStatus)
IF (ErrorStatus > 0) THEN
  cmessage = 'Fatal error reading NLSIZES namelist'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  WRITE(6, NML = NLSIZES)
END IF

READ(UNIT = nft, NML = ATM_SIZES, IOSTAT = ErrorStatus)
IF (ErrorStatus > 0) THEN
  cmessage = 'Fatal error reading ATM_SIZES namelist'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  WRITE(6, NML = ATM_SIZES)
END IF

ntiles=9
CALL calc_ntiles(l_aggregate,npft,nnvg,ntiles)

! Set output variables as required
output_grid % glob_p_rows       = global_rows
output_grid % glob_p_row_length = global_row_length
output_grid % glob_r_rows       = river_rows
output_grid % glob_r_row_length = river_row_length
output_grid % model_levels      = model_levels
output_grid % wet_levels        = wet_levels
output_grid % ozone_levels      = ozone_levels
output_grid % st_levels         = st_levels
output_grid % bl_levels         = bl_levels

! Default value - not required in the dump
output_grid % cloud_levels      = IMDI

! Set some land/sea maks values
glob_land_out                   = land_field
output_grid % glob_land_field   = land_field

! Set values for tropopause-based ozone.
! In future developments, obtain directly from namelist.
tpps_ozone_levels               = ozone_levels

RETURN
END SUBROUTINE rcf_readnl_nlsizes

END MODULE rcf_readnl_nlsizes_mod
