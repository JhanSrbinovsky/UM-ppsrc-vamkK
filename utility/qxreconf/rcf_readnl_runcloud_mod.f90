! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads the RUN_Cloud namelist
!
MODULE rcf_readnl_runcloud_mod

IMPLICIT NONE

! Description:
!  rcf_readnl_runcloud reads in the RUN_Cloud namelist.
!
! Method:
!  Data read in via unit number passed as argument.
!  Only one element of this namelist is required by the reconfiguration: rhcrit.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 v8.5 programming standards.

CONTAINS

SUBROUTINE rcf_readnl_runcloud (nft)

USE cloud_inputs_mod, ONLY:                                              &
  RUN_Cloud,                                                             &
  rhcrit, l_pc2, l_cld_area

USE rcf_grid_type_mod, ONLY:                                             &
  output_grid

USE atmos_max_sizes, ONLY:                                               &
  model_levels_max

USE ereport_mod, ONLY:                                                   &
  ereport

USE PrintStatus_mod, ONLY:                                               &
  PrintStatus, PrStatus_Oper

USE UM_parvars, ONLY:                                                    &
  mype

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)  :: nft  ! unit number

! Local variables
INTEGER              :: ErrorStatus

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'rcf_readnl_runcloud'
CHARACTER (LEN=80)           :: cmessage


! Defaults are set in the namelist module
! Read the namelist
READ(UNIT = nft, NML = RUN_Cloud, IOSTAT = ErrorStatus)
IF (ErrorStatus > 0) THEN
  cmessage = 'Fatal error reading RUN_Cloud namelist'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  WRITE(6, NML = RUN_Cloud)
END IF

! Allocate space and fill up relevant part of output_grid
ALLOCATE ( output_grid % rhcrit( output_grid % model_levels ) )

output_grid % rhcrit( 1 : output_grid % model_levels ) =                 &
              rhcrit( 1 : output_grid % model_levels )

RETURN
END SUBROUTINE rcf_readnl_runcloud

END MODULE rcf_readnl_runcloud_mod
