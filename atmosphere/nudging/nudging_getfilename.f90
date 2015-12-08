! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------
! Description:
!  Converts the model time and date into the approppriate filenames
!  for the files containing the meteorological analyses.
!  Note it is hard wired for the format of ECMWF files obtained from
!  BADC via BDAN. Any changes to the file naming require modification.

!  Part of the Nudged model (see nudging_main.F90)

!  Called from NUDGING_MAIN.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Nudging

! Code description:
!   Language: FORTRAN 90

!------------------------------------------------------------------
SUBROUTINE nudging_getfilename(        &
 i_year,                               &   ! Model year
 i_month,                              &   ! Model month
 i_day,                                &   ! Model day
 i_hour,                               &   ! Model hour
 dataname1,                            &   ! Return first data date string
 dataname2,                            &   ! Return second data date string
 timestep1,                            &   ! Timestep of first data point
 timestep2,                            &   ! Timestep of second data point
 filetype,                             &   ! Naming convention
 debug)                                    ! Debug flag

USE nudging_control

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

USE ereport_mod, ONLY : ereport
IMPLICIT NONE

INTEGER, INTENT(IN) :: i_year     ! Model year
INTEGER, INTENT(IN) :: i_month    ! Model month
INTEGER, INTENT(IN) :: i_day      ! Model day
INTEGER, INTENT(IN) :: i_hour     ! Model hour

CHARACTER(LEN=256), INTENT(OUT) :: dataname1  ! first full filename
CHARACTER(LEN=256), INTENT(OUT) :: dataname2  ! second full filename

INTEGER, INTENT(OUT) :: timestep1  ! first timestep
INTEGER, INTENT(OUT) :: timestep2  ! second timestep
INTEGER, INTENT(IN) :: filetype    ! Naming convention
INTEGER, INTENT(IN) :: debug       ! Debug flag

!**********************************************
! Variables used in writing filename
CHARACTER(LEN=10)  :: date1        ! first date
CHARACTER(LEN=10)  :: date2        ! second date
CHARACTER(LEN=200) :: datafile1    ! first filename
CHARACTER(LEN=200) :: datafile2    ! second filename

INTEGER                   :: errcode              ! error status
INTEGER(KIND=jpim), PARAMETER :: zhook_in = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NUDGING_GETFILENAME',zhook_in,zhook_handle)

!***********************************************************************
! End of Header

! Standard Entry Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                        &
 ' NUDGING_GETFILENAME: Entering Routine'
END IF

! Set timestep that we read from the file
! As only have one timestep per file then this is always trivially 1
! Really kept this here in case we want to have one file per day.
timestep1 = 1
timestep2 = 1

! Calculate the data string
! DEPENDS ON: nudging_getdate
CALL nudging_getdate(          &
  i_year,                      &  ! Model year
  i_month,                     &  ! Model month
  i_day,                       &  ! Model day
  i_hour,                      &  ! Model hour
  date1,                       &  ! Return first data date string
  date2,                       &  ! Return second data date string
  debug)                          ! Miscellanea

! Combine the date with the standard file beginnings and
! endings. Vary these according to the sources of analyses
SELECT CASE(ndg_analysis_source)
  CASE(0)
    datafile1 = file_model_basis // date1 // file_end
    datafile2 = file_model_basis // date2 // file_end
  CASE(1)
    datafile1 = file_pres_basis  // date1 // file_end
    datafile2 = file_pres_basis  // date2 // file_end
  CASE(2)
    datafile1 = file_um_basis  // date1 // file_um_end
    datafile2 = file_um_basis  // date2 // file_um_end
  CASE(3)
    datafile1 = file_jra_basis  // date1 // file_jra_end
    datafile2 = file_jra_basis  // date2 // file_jra_end
  CASE DEFAULT
    nmessage = 'Unknown Analysis data source.'
    errcode = 999

    CALL ereport('NUDGING_GETFILENAME',errcode,nmessage)
END SELECT

! Combine the filenames and filepaths to produce full string
! routines that read the netcdf files
WRITE(dataname1,'(a,a1,a)')TRIM(ndg_datapath),'/',datafile1
WRITE(dataname2,'(a,a1,a)')TRIM(ndg_datapath),'/',datafile2

! Standard Exit Comment
IF(debug > 10) THEN
  WRITE(OUT,*) 'NUDGING_GETFILENAME: datafiles: ',TRIM(dataname1),     &
  ' and ',TRIM(dataname2)
END IF

!********************************************************************

! Standard Exit Comment
IF(debug > 10) THEN
  WRITE(OUT,*)                                                         &
 ' NUDGING_GETFILENAME: Leaving Routine'
END IF

IF (lhook) CALL dr_hook('NUDGING_GETFILENAME',zhook_out,zhook_handle)

RETURN
END SUBROUTINE nudging_getfilename

