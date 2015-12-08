! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Top level for reading in reconfiguration namelists.

MODULE Rcf_read_namelists_Mod

IMPLICIT NONE

!  Subroutine Rcf_Read_Namelists - read the rcf namelists.
!
! Description:
!   Read the namelists, assigning and freeing Fortran units as
!   required.
!
! Method:
!   The namelists are provided in one file - RECONA and RECONO
!   for Atmos and Ocean respectively.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

SUBROUTINE rcf_read_namelists ( )

USE Submodel_Mod, ONLY : &
    Submodel_Ident,      &
    Atmos_IM

USE Rcf_readnl_nsubmodl_Mod,  ONLY :  Rcf_readnl_nsubmodl
USE Rcf_readnl_nlsizes_Mod,   ONLY :  Rcf_readnl_nlsizes
USE Rcf_readnl_ustsnum_Mod,   ONLY :  Rcf_readnl_ustsnum
USE Rcf_readnl_uancnum_Mod,   ONLY :  Rcf_readnl_uancnum
USE Rcf_readnl_stshcomp_Mod,  ONLY :  Rcf_readnl_stshcomp
USE Rcf_readnl_nlstcatm_Mod,  ONLY :  Rcf_readnl_nlstcatm
USE Rcf_readnl_runukca_Mod,   ONLY :  Rcf_readnl_runukca
USE Rcf_readnl_rundust_Mod,   ONLY :  Rcf_readnl_rundust
USE Rcf_readnl_rungwd_Mod,    ONLY :  Rcf_readnl_rungwd
USE Rcf_readnl_runmurk_Mod,   ONLY :  Rcf_readnl_runmurk
USE Rcf_readnl_runbl_Mod,     ONLY :  Rcf_readnl_runbl
USE Rcf_readnl_runrivers_Mod, ONLY :  Rcf_readnl_runrivers
USE Rcf_readnl_runprecip_Mod, ONLY :  Rcf_readnl_runprecip
USE Rcf_readnl_runcloud_Mod,  ONLY :  Rcf_readnl_runcloud
USE Rcf_readnl_ancilcta_Mod,  ONLY :  Rcf_readnl_ancilcta
USE Rcf_readnl_Recon_Mod,     ONLY :  Rcf_readnl_recon
USE Rcf_readnl_vertical_Mod,  ONLY :  Rcf_readnl_vertical
USE Rcf_readnl_horizont_Mod,  ONLY :  Rcf_readnl_horizont
USE Rcf_readnl_headers_Mod,   ONLY :  Rcf_readnl_headers
USE Rcf_readnl_items_Mod,     ONLY :  Rcf_readnl_items
USE Rcf_readnl_trans_Mod,     ONLY :  Rcf_readnl_trans
USE read_jules_namelists_mod,    ONLY :  read_jules_nstypes,   &
                                         read_jules_switches,  &
                                         read_jules_snow_param,&
                                         read_urban_switches
USE Rcf_readnl_runradiation_Mod,  ONLY :  Rcf_readnl_runradiation
USE Rcf_readnl_runconvection_Mod, ONLY :  Rcf_readnl_runconvection

USE science_fixes_mod
USE check_iostat_mod

USE Ereport_Mod, ONLY :&
    Ereport

USE PrintStatus_mod, ONLY : &
    PrintStatus,            &
    PrStatus_Oper

USE Rcf_FortranIO_Mod, ONLY : &
    Rcf_Get_Unit,             &
   Rcf_Free_Unit,             &
   Max_Filename_Len

USE Rcf_Recon_Mod, ONLY : &
    TRANS

IMPLICIT NONE
! Arguments

! Local Variables/Paramters
CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'Rcf_Read_Namelists'
CHARACTER (LEN=80)                 :: Cmessage
CHARACTER (LEN=Max_Filename_Len)   :: FileName
INTEGER                            :: ErrorStatus
INTEGER                            :: nft
INTEGER                            :: nft_sizes
INTEGER                            :: nft_shared
INTEGER                            :: nft_vertlevs
INTEGER                            :: nft_horizgrid
INTEGER                            :: status
LOGICAL                            :: l_exist

! Get a Fortran Unit for the namelists
CALL Rcf_Get_Unit( nft )
CALL Rcf_Get_Unit( nft_sizes )
CALL Rcf_Get_Unit( nft_shared )
CALL Rcf_Get_Unit( nft_vertlevs )
CALL Rcf_Get_Unit( nft_horizgrid )

! ---------------------------------------------------
! Read namelists common to the UM and reconfiguration
! Part 1 - miscellaneous namelists in file SHARED
! ---------------------------------------------------

! First, find the namelist filename from Env Vars
CALL Fort_Get_Env( 'SHARED_NLIST', 12, FileName,                   &
                    Max_Filename_Len, ErrorStatus )

IF ( ErrorStatus /= 0 ) THEN
  ErrorStatus = 30
  Cmessage = 'Unable to Obtain Shared Namelists Filename from Environment'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

FileName = TRIM( FileName )

! Check file exists - Do all this by hand (not using file_open)
! as wish to do this on *all* PEs together
INQUIRE( file=FileName, exist=l_exist )

IF ( .NOT. l_exist ) THEN
  ErrorStatus = 40
  Cmessage = 'Shared namelists file does not exist!'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! Open the file containing namelists from SHARED file
OPEN( UNIT=nft_shared, FILE=FileName, IOSTAT=status )

IF ( PrintStatus >= PrStatus_Oper ) THEN
  WRITE (6,'(''Shared namelists file: '',a)') FileName
END IF

! ---------------------------------------------------
! read in any temporary UM fixes and warn if required
! ---------------------------------------------------

READ (UNIT=nft_shared, NML=temp_fixes, IOSTAT=errorstatus)
CALL check_iostat(errorstatus, "namelist temp_fixes")
CALL warn_temp_fixes()

! ---------------------------------------------------


CALL read_jules_switches (nft_shared)
CALL read_urban_switches (nft_shared)
CALL rcf_readnl_ustsnum  (nft_shared)

! ---------------------------------------------------
! Read namelists common to the UM and reconfiguration
! Part 2 - model SIZES namelists
! ---------------------------------------------------

! First, find the namelist filename from Env Vars
CALL Fort_Get_Env( 'SIZES_NLIST', 11, FileName,                    &
                    Max_Filename_Len, ErrorStatus )

IF ( ErrorStatus /= 0 ) THEN
  ErrorStatus = 10
  Cmessage = 'Unable to Obtain Sizes Namelists Filename from Environment'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

FileName = TRIM( FileName )

! Check file exists - Do all this by hand (not using file_open)
! as wish to do this on *all* PEs together
INQUIRE( file=FileName, exist=l_exist )

IF ( .NOT. l_exist ) THEN
  ErrorStatus = 20
  Cmessage = 'Sizes namelists file does not exist!'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! Open the file containing namelists from SIZES file
OPEN( UNIT=nft_sizes, FILE=FileName, IOSTAT=status )

IF ( PrintStatus >= PrStatus_Oper ) THEN
  WRITE (6,'(''Sizes namelists file: '',a)') FileName
END IF

CALL rcf_readnl_nsubmodl (nft_sizes)
CALL rcf_readnl_nlsizes  (nft_sizes)

! ------------------------------------------------
! Begin reading reconfiguration-specific namelists
! ------------------------------------------------

! First, find the namelist filename from Env Vars
CALL Fort_Get_Env( 'RCF_NAMELIST', 12, FileName,                   &
                    Max_Filename_Len, ErrorStatus )

IF ( ErrorStatus /= 0 ) THEN
  ErrorStatus = 50
  Cmessage = 'Unable to Obtain Namelist Filename from Environment'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

FileName = TRIM( FileName )

! Check file exists - Do all this by hand (not using file_open)
! as wish to do this on *all* PEs together
INQUIRE( file=FileName, exist=l_exist )

IF ( .NOT. l_exist ) THEN
  ErrorStatus = 60
  Cmessage = 'Namelist file does not exist!'
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! Open the file containing namelists from RECONA/RECONO file
OPEN( UNIT=nft, FILE=FileName, IOSTAT=status )

IF ( PrintStatus >= PrStatus_Oper ) THEN
  WRITE (6,'(''Namelist file: '',a)') FileName
END IF

! RECON must be read after NLSIZES due to using model_levels
CALL rcf_readnl_recon    (nft)
CALL rcf_readnl_uancnum  (nft)

! -----------------------------------------
! Resume reading SHARED and SIZES namelists
! -----------------------------------------

! Submodel dependent namelists
! Submodel_ident is read in from RECON
IF (Submodel_Ident == Atmos_IM) THEN
  CALL rcf_readnl_nlstcatm (nft_shared)
  CALL read_jules_nstypes (nft_sizes)
END IF
CALL rcf_readnl_stshcomp (nft_sizes)
CALL rcf_readnl_runukca (nft_shared)
CALL rcf_readnl_rungwd (nft_shared)
CALL rcf_readnl_runmurk (nft_shared)
CALL rcf_readnl_runconvection (nft_shared)
CALL rcf_readnl_runbl (nft_shared)
CALL rcf_readnl_runrivers (nft_shared)
CALL rcf_readnl_runprecip (nft_shared)
CALL rcf_readnl_runradiation (nft_shared)
CALL rcf_readnl_rundust  (nft_shared)
CALL rcf_readnl_runcloud (nft_shared)
CALL rcf_readnl_ancilcta (nft_shared)
CALL read_jules_snow_param (nft_shared)

! ------------------------
! Vertical Levels Namelist
! ------------------------

IF (Submodel_Ident == Atmos_IM) THEN

! Find the filename containing vertical levels from Env Vars
  CALL Fort_Get_Env('VERT_LEV', 8, FileName,                       &
                     Max_Filename_Len, ErrorStatus)

  IF ( ErrorStatus /= 0 ) THEN
    ErrorStatus = 70
    Cmessage =  &
    'Unable to Obtain Vertical Levels Filename from Environment'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  FileName = TRIM( FileName )

! Check file exists - Do all this by hand (not using file_open)
! as wish to do this on *all* PEs together
  INQUIRE( file=FileName, exist=l_exist )

  IF ( .NOT. l_exist ) THEN
    WRITE (6,'(''Vertical Levels file: '',a)') FileName
    ErrorStatus = 80
    Cmessage = 'Vertical Levels Namelist file does not exist!'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

! Open the file containing vertical levels
  OPEN( UNIT=nft_vertlevs, FILE=FileName, IOSTAT=status )

  IF ( PrintStatus >= PrStatus_Oper ) THEN
    WRITE (6,'(''Vertical Levels file: '',a)') FileName
  END IF

END IF

! For vertical, three namelists - VERTICAL, 
!  VERTLEVS and JULES_SOIL_PARAM - are read in.
CALL rcf_readnl_vertical (nft, nft_vertlevs, nft_shared)
CALL rcf_readnl_horizont (nft, nft_horizgrid)
CALL rcf_readnl_headers  (nft)
CALL rcf_readnl_items    (nft)

IF (TRANS) THEN
  CALL rcf_readnl_trans    (nft)
END IF

CLOSE( Unit=nft )
CLOSE( Unit=nft_sizes )
CLOSE( Unit=nft_shared )
CLOSE( Unit=nft_vertlevs )

CALL Rcf_Free_Unit( nft )
CALL Rcf_Free_Unit( nft_sizes )
CALL Rcf_Free_Unit( nft_shared )
CALL Rcf_Free_Unit( nft_vertlevs )
CALL Rcf_Free_Unit( nft_horizgrid )

RETURN

END SUBROUTINE  Rcf_read_namelists

END MODULE Rcf_read_namelists_Mod
