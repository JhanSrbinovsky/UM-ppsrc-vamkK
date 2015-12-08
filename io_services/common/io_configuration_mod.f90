! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: C96

! Module defining buffering constants and namelists

MODULE io_configuration_mod

IMPLICIT NONE

INTEGER, PARAMETER :: io_namelist_unit = 127
INTEGER, PARAMETER :: io_timing_on     = 1
INTEGER, PARAMETER :: io_timing_off    = 0

INTEGER :: io_wbuffer_size    = 524288    ! In words. Size of buffer when 
                                         ! buffered I/O is active.
                                         ! Default is 4MB.

INTEGER :: io_rbuffer_size    = 524288    ! In words. Size of buffer when 
                                          ! buffered I/O is active.
                                          ! Default is 4MB.

INTEGER :: io_rbuffer_count    = 4        ! Number of read buffers/stream

INTEGER :: io_rbuffer_update   = 1        ! 
INTEGER :: io_rbuffer_prefetch = 1        ! 

INTEGER :: io_data_alignment = 524288    ! In words. Alignment of data block in
                                         ! dumps/fieldsfiles. Defaults is 4MB.

INTEGER :: io_field_padding  = 512       ! Previously um_sector_size. The 
                                         ! padding between  fields. Now mostly
                                         ! defunct due to buffered I/O.

INTEGER :: io_timing         = io_timing_off

LOGICAL :: io_external_control = .FALSE.

! A namelist will allow us to control these via the UMUI.
NAMELIST /io_control/    & 
    io_wbuffer_size,     &
    io_rbuffer_size,     &
    io_rbuffer_count,    &
    io_rbuffer_update,   &
    io_rbuffer_prefetch, &
    io_data_alignment,   &
    io_field_padding,    &
    io_timing,           &
    io_external_control

CONTAINS
  SUBROUTINE ioLoadConfig()
    USE FilenameLength_mod, ONLY : &
        filenamelength
    USE ereport_mod
    USE check_iostat_mod
    USE application_description, ONLY : isSmallExec
    IMPLICIT NONE
    INTEGER :: error
    INTEGER :: errCode
    CHARACTER (LEN=FileNameLength) :: NamelistFile
    CHARACTER (LEN=132)            :: message

    
!*L --------------------- Comdeck: CENVIR   ----------------------------
!
!   Purpose: COMDECK defining Character enviroment variables used
!            by portable IO to open and close files
!
!
! Type declarations
!
      CHARACTER(LEN=8) FT_ENVIRON(199)  ! Array holding enviroment variables
!                                  for filenames
      INTEGER     LEN_FT_ENVIR(199) ! character length of each variable
!
!
!Common Blocks for character and integer arrays
!
      COMMON/CENVIR/FT_ENVIRON
      COMMON/CLENVIR/LEN_FT_ENVIR
!
    NamelistFile = "dummy file"
    errCode=-10

    CALL FORT_GET_ENV(ft_environ(IO_Namelist_Unit), &
        len_ft_envir(IO_Namelist_Unit),             &
        NameListFile,FileNameLength,error)
    IF (error==0) THEN
      OPEN(IO_Namelist_Unit,FILE=NameListFile,FORM='FORMATTED',&
          STATUS='OLD',IOSTAT=error)
      IF ( error == 0 ) THEN
        ! Read the I/O buffering control namelist too
        READ(IO_Namelist_Unit, io_control, IOSTAT=error)
        CALL check_iostat(error,'NAMELIST - io_control')
        CLOSE(IO_Namelist_Unit)
      ELSE
        ! Stifle warnings from small execs which caused user concern.
        ! Small execs don't need to find the file as default values work
        ! perfectly well.
        IF (.NOT.isSmallExec()) THEN
          WRITE(MESSAGE,'(A,A,A)')'Failed to open IO control file (', &
              TRIM(NameListFile),')'
          CALL ereport('ioLoadConfig',errCode,message)
        END IF
      END IF
    ELSE
      WRITE(MESSAGE,'(A,A)')'Failed to get filename for IO control file ', &
          ' from environment'
      CALL ereport('ioLoadConfig',errCode,message)
    END IF
    
  END SUBROUTINE ioLoadConfig
  
END MODULE io_configuration_mod
