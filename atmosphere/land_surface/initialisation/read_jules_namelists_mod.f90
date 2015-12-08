! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Wrapper module containing subroutines for reading JULES namelists
!
MODULE read_jules_namelists_mod

! Description:
!  Contains read_jules_<namelist> and read_<urban_namelist> subroutines
!  for reading namelists into JULES during a UM-JULES job.
!
! Method:
!  The unit number holding the namelist is passed as the sole argument 
!  to each file.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Land
!
! Code Description:
!  Language: FORTRAN 95.
!  This code is written to UMDP3 v8.5 programming standards.

  USE check_iostat_mod, ONLY:                                          &
    check_iostat

  USE PrintStatus_mod,  ONLY:                                          &
    PrintStatus, PrStatus_Oper

  USE UM_parvars,       ONLY:                                          &
    mype

  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook

  IMPLICIT NONE

! Local variables common to each subroutine
  INTEGER, PRIVATE             :: errorstatus
  REAL(KIND=jprb), PRIVATE     :: zhook_handle

  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CONTAINS

SUBROUTINE read_jules_csmin (unitnumber)

! Description:
!  Read the JULES_CSMIN namelist

  USE csmin,            ONLY:                                          &
    jules_csmin

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_JULES_CSMIN',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=jules_csmin, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist jules_csmin")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = jules_csmin)
  END IF

  IF (lhook) CALL dr_hook('READ_JULES_CSMIN',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_jules_csmin


SUBROUTINE read_jules_elevate (unitnumber)

! Description:
!  Read the JULES_ELEVATE namelist

  USE c_elevate,        ONLY:                                          &
    jules_elevate

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_JULES_ELEVATE',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=jules_elevate, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist jules_elevate")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = jules_elevate)
  END IF

  IF (lhook) CALL dr_hook('READ_JULES_ELEVATE',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_jules_elevate


SUBROUTINE read_jules_nstypes (unitnumber)

! Description:
!  Read the JULES_NSTYPES namelist

  USE nstypes,          ONLY:                                          &
    jules_nstypes,                                                     &
! namelist variables:
    ntype,            npft,             nnvg

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_JULES_NSTYPES',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=jules_nstypes, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist jules_nstypes")

! Initialise number of surface types from namelist values
  ntype = npft + nnvg

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = jules_nstypes)
  END IF

  IF (lhook) CALL dr_hook('READ_JULES_NSTYPES',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_jules_nstypes


SUBROUTINE read_jules_nvegparm (unitnumber)

! Description:
!  read the JULES_NVEGPARM namelist

  USE nvegparm_io,      ONLY:                                          &
    jules_nvegparm

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_JULES_NVEGPARM',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=jules_nvegparm, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist jules_nvegparm")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = jules_nvegparm)
  END IF

  IF (lhook) CALL dr_hook('READ_JULES_NVEGPARM',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_jules_nvegparm


SUBROUTINE read_jules_pftparm (unitnumber)

! Description:
!  Read the JULES_PFTPARM namelist

  USE pftparm_io,       ONLY:                                          &
    jules_pftparm

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_JULES_PFTPARM',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=jules_pftparm, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist jules_pftparm")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = jules_pftparm)
  END IF

  IF (lhook) CALL dr_hook('READ_JULES_PFTPARM',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_jules_pftparm


SUBROUTINE read_jules_rad_param (unitnumber)

! Description:
!  read the JULES_RAD_PARAM namelist

  USE rad_param,        ONLY:                                          &
    jules_rad_param

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_JULES_RAD_PARAM',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=jules_rad_param, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist jules_rad_param")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = jules_rad_param)
  END IF

  IF (lhook) CALL dr_hook('READ_JULES_RAD_PARAM',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_jules_rad_param


SUBROUTINE read_jules_seed (unitnumber)

! Description:
!  Read the JULES_SEED namelist

  USE seed,             ONLY:                                          &
    jules_seed

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_JULES_SEED',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=jules_seed, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist jules_seed")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = jules_seed)
  END IF

  IF (lhook) CALL dr_hook('READ_JULES_SEED',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_jules_seed


SUBROUTINE read_jules_sigm (unitnumber)

! Description:
!  Read the JULES_SIGM namelist

  USE sigm,             ONLY:                                          &
    jules_sigm

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_JULES_SIGM',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=jules_sigm, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist jules_sigm")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = jules_sigm)
  END IF

  IF (lhook) CALL dr_hook('READ_JULES_SIGM',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_jules_sigm


SUBROUTINE read_jules_snow_param (unitnumber)

! Description:
!  Read the JULES_SNOW_PARAM namelist

  USE snow_param,       ONLY:                                          &
    jules_snow_param

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_JULES_SNOW_PARAM',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=jules_snow_param, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist jules_snow_param")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = jules_snow_param)
  END IF

  IF (lhook) CALL dr_hook('READ_JULES_SNOW_PARAM',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_jules_snow_param


SUBROUTINE read_jules_soil_param (unitnumber)

! Description:
!  Read the JULES_SOIL_PARAM namelist

  USE soil_param,       ONLY:                                          &
    jules_soil_param

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_JULES_SOIL_PARAM',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=jules_soil_param, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist jules_soil_param")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = jules_soil_param)
  END IF

  IF (lhook) CALL dr_hook('READ_JULES_SOIL_PARAM',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_jules_soil_param


SUBROUTINE read_jules_surf_param (unitnumber)

! Description:
!  Read the JULES_SURF_PARAM namelist

  USE surf_param,       ONLY:                                          &
    jules_surf_param

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_JULES_SURF_PARAM',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=jules_surf_param, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist jules_surf_param")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = jules_surf_param)
  END IF

  IF (lhook) CALL dr_hook('READ_JULES_SURF_PARAM',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_jules_surf_param


SUBROUTINE read_jules_switches (unitnumber)

! Description:
!  Read the JULES_SWITCHES namelist

  USE switches,         ONLY:                                          &
    jules_switches

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_JULES_SWITCHES',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=jules_switches, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist jules_switches")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = jules_switches)
  END IF

  IF (lhook) CALL dr_hook('READ_JULES_SWITCHES',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_jules_switches


SUBROUTINE read_jules_triffid (unitnumber)

! Description:
!  Read the JULES_TRIFFID namelist

  USE trif_io,          ONLY:                                          &
    jules_triffid

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_JULES_TRIFFID',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=jules_triffid, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist jules_triffid")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = jules_triffid)
  END IF

  IF (lhook) CALL dr_hook('READ_JULES_TRIFFID',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_jules_triffid


SUBROUTINE read_urban2t_param (unitnumber)

! Description:
!  Read the URBAN2T_PARAM namelist

  USE urban_param,      ONLY:                                          &
    urban2t_param

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_URBAN2T_PARAM',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=urban2t_param, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist urban2t_param")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = urban2t_param)
  END IF

  IF (lhook) CALL dr_hook('READ_URBAN2T_PARAM',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_urban2t_param


SUBROUTINE read_urban_switches (unitnumber)

! Description:
!  Read the URBAN_SWITCHES namelist

  USE switches_urban,   ONLY:                                          &
    urban_switches

  IMPLICIT NONE

! Subroutine arguments
  INTEGER, INTENT(IN) :: unitnumber

  IF (lhook) CALL dr_hook('READ_URBAN_SWITCHES',zhook_in,zhook_handle)

  READ (UNIT=unitnumber, NML=urban_switches, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist urban_switches")

  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE (6, NML = urban_switches)
  END IF

  IF (lhook) CALL dr_hook('READ_URBAN_SWITCHES',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE read_urban_switches

END MODULE read_jules_namelists_mod
