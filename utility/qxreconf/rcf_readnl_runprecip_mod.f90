! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Reads the RUN_Precip namelist
!
MODULE rcf_readnl_runprecip_mod

IMPLICIT NONE

! Description:
!  rcf_readnl_runprecip reads in the RUN_Precip namelist.
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

SUBROUTINE rcf_readnl_runprecip (nft)

USE mphys_inputs_mod, ONLY:                                              &
  RUN_Precip, l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain        

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

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'rcf_readnl_runprecip'
CHARACTER (LEN=80)           :: cmessage

! Defaults are set in the namelist module
! Read the namelist
READ(UNIT = nft, NML = RUN_Precip, IOSTAT = ErrorStatus)
IF (ErrorStatus > 0) THEN
  cmessage = 'Fatal error reading RUN_Precip namelist'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

IF ( PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  WRITE(6, NML = RUN_Precip)
END IF

RETURN
END SUBROUTINE rcf_readnl_runprecip

END MODULE rcf_readnl_runprecip_mod
