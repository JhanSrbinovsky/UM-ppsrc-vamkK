! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

MODULE io
! A module to contain Fortran interfaces to file operations

! Description:
!  This module provides BUFFOUT,BUFFIN,OPEN,CLOSE,SETPOS routines
!  for the Parallel Unified Model, Reconfiguration, and Small Executables.
!
! Complete list of interfaces
!
!             buffin,buffin64,buffin32
!             buffout,buffout64,buffo32
!
!  Note that routines are selected by array type, and hence the generic
!  (non-lengthed) interface is preferred. The array type supplied MUST
!  reflect its content in order to facilitate correct byteswapping.
!
! Complete list of routines
!
! ioInit()
! ioShutdown()
! checkUnit(unit)
! file_open
!  (unit,env,env_len,read_write,name_in_environ,error,ioLocality)
! file_close(unit,env,env_len,name_in_environ,delete,err)
! buffout64_r(unit,array,numWords,wordsWritten,errorCode)
! buffout64_i(unit,array,numWords,wordsWritten,errorCode)
! buffout64_l(unit,array,numWords,wordsWritten,errorCode)
! buffout64_r2D(unit,array,numWords,wordsWritten,errorCode)
! buffout64_i2D(unit,array,numWords,wordsWritten,errorCode)
! buffout64_l2D(unit,array,numWords,wordsWritten,errorCode)
! buffin64_r(unit,array,numWords,wordsRead,errorCode)
! buffin64_i(unit,array,numWords,wordsRead,errorCode)
! buffin64_l(unit,array,numWords,wordsRead,errorCode)
! buffin64_i2D(unit,array,numWords,wordsRead,errorCode)
! buffin64_l2D(unit,array,numWords,wordsRead,errorCode)
! buffin64_r2D(unit,array,numWords,wordsRead,errorCode)
! buffin32_r(unit,array,numWords,wordsRead,errorCode)
! buffin32_i(unit,array,numWords,wordsRead,errorCode)
! buffin32_i2D(unit,array,numWords,wordsRead,errorCode)
! buffin32_l(unit,array,numWords,wordsRead,errorCode)
! buffin_character(unit,array,numWords,csize,wordsRead,errorCode)
! setpos(unit,ipos,icode)
! setpos8(unit,ipos,icode)
! setpos32(unit,ipos,icode)
! ioFileState(unit,state)
! readbleWords(unit,wordLength)
! set_unit_bcast_flag(unit)
! clear_unit_bcast_flag(unit)
! logical broadcast_read(unit)
!
! PRIVATE:
!       setAttribute(var,mask)
!       clearAttribute(var,mask)
!       LOGICAL testAttribute(var,mask)
!       performOperation(unit)
!       io_ereport
!
!  Note that setpos expresses units of the default integer/real type
!  and there is currently no explicit setpos64, and that 'getpos'
!  functionality is not exposed at the fortran layer - use ioFileState
!  instead to obtain the current file position in bytes.
!
! Method:
!  The C interfaces *_single are called only by the IO module
!  appropriately on tasks designated to perform IO. Calls to C
!  versions should not be made elsewhere. The designation of who
!  performs IO is specified at file open time via the ioLocality
!  optional argument and is persistant until file close (with the
!  exception of buffin, for which broadcasting can be switched on and off)
!
!  Failing to provide the ioLocality argument results in a default
!  behaviour in which IO Servers perform IO if available for writeable
!  files, otherwise rank zero. Rank zero will perform IO for non-writeable
!  files (unless overriden via the IO Server configuration). By default
!  buffin operations are collective and broadcast read data to all tasks.
!
!  Specifying ioLocality=ioAllLocal results in each task performing IO
!  completely independently (take care to write to different files).
!
!  When using any option other than ioAllLocal, calls should be regarded
!  as collective and the call made on all tasks, with correct arguments
!  (e.g. unit numbers, modes etc) even if you think the task
!  will not actually participate. This allows for future optimisations
!  to perform in parallel and/or distribute IO load.
!
!  Some exceptions are made to the above. They are there to support the
!  existing UM code and such use is depricated. Modification to existing UM IO
!  code should attempt to make it collective.

  USE um_types
  USE filenamelength_mod, ONLY :                                               &
      maxFileNameLength => fileNameLength
  USE ios
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE PrintStatus_mod
  USE UM_ParVars, ONLY : mype, nproc
  USE IO_dependencies ! Lower layer dependencies

  IMPLICIT NONE
  PRIVATE


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters indicating the nature of the file open operation
  INTEGER, PUBLIC, PARAMETER :: ioOpenReadOnly        = 0
  INTEGER, PUBLIC, PARAMETER :: ioOpenReadWrite       = 1
  INTEGER, PUBLIC, PARAMETER :: ioOpenWriteOnly       = -99999

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters indicating known file types, for provision via the fileType
! optional argument to file_open (and internally)
  INTEGER, PUBLIC, PARAMETER :: ioFileTypeUnknown = 0
  INTEGER, PUBLIC, PARAMETER :: ioFileTypeUM      = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters indicating how the passed string to file_open/close encodes
! the filename.
  INTEGER, PUBLIC, PARAMETER :: ioNameInEnv       = 0
  INTEGER, PUBLIC, PARAMETER :: ioNameProvided    = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters indicating whether file should be deleted on close.
  INTEGER, PUBLIC, PARAMETER :: ioNoDelete        = 0
  INTEGER, PUBLIC, PARAMETER :: ioDelete          = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters for controlling where IO happens !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER, PUBLIC, PARAMETER :: ioNotKnown=-1
  INTEGER, PUBLIC, PARAMETER :: ioNotSet=0
  INTEGER, PUBLIC, PARAMETER :: ioLocal=1
  INTEGER, PUBLIC, PARAMETER :: ioAllLocal=2
  INTEGER, PUBLIC, PARAMETER :: ioRemote=4
  INTEGER, PUBLIC, PARAMETER :: ioReadReplicate=8
  INTEGER, PUBLIC, PARAMETER :: ioReadOnly=16
  INTEGER, PUBLIC, PARAMETER :: ioOpen=32
  INTEGER, PUBLIC, PARAMETER :: ioWriteOnly=64
!Some shortcuts
   ! Good choice for wr files;
  INTEGER, PUBLIC, PARAMETER :: ioRemoteReplicate=IOR(ioRemote,ioReadReplicate)
   ! Good choice for ro files;
  INTEGER, PUBLIC, PARAMETER :: ioLocalReplicate=IOR (ioLocal, ioReadReplicate)
! The default
  INTEGER, PUBLIC, PARAMETER :: ioDefault=ioRemoteReplicate
!Notes:
! The parameter is supplied at file open and persists until the
! file is closed
!
! ioAllLocal : all tasks will open/write/close the file
!  - It is the programmers responsibility to avoid filename collisions
!  - It is the programmers responsibility to make calls from appropriate PEs
! ioLocal : Rank 0 will open/write/close the file
!  - Other tasks will return with a success flag but not do anything
! ioRemote : Rank 0 will open/write/close the file using a remote IO service
!  - All tasks will return with a success flag
!  - No action taken on ranks that are not 0
! ioReadReplicate : reads will broadcast to all tasks
!  - cannot be or'd with ioAllLocal
!
! ****** A unit number may not be open remotely and locally ****
! ****** Remote IO from all tasks is not supported (yet)    ****

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters for manipulating IO
  INTEGER, PRIVATE, PARAMETER :: ioNoLocation=-1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Parameters for controlling how 'faulty' IO operations are dealt with !
  INTEGER, PUBLIC, PARAMETER :: ioPolicyStrict=27
  INTEGER, PUBLIC, PARAMETER :: ioPolicyConsistent=39
  INTEGER, PUBLIC, PARAMETER :: ioPolicySloppy=4123
  INTEGER, PUBLIC, PARAMETER :: ioPolicyLax=976
  INTEGER, PUBLIC, PARAMETER :: ioPolicyDefault=ioPolicyStrict

!Notes:
! The parameter is interpreted by file_open as
!
! ioPolicyStrict     : It is required that the unit is closed
! ioPolicyConsistent : If the file is already open on this unit,
!                      with the right permissions, position will
!                      be reset to zero
! ioPolicySloppy     : If the file is already open on this unit,
!                      the file position will be reset to 0
! ioPolicyLax      : The file will be opened as requested and
!                      any existing file open on that unit will
!                      be closed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Private storage for the above attributes.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! General parameters
  INTEGER, PARAMETER            :: maxUnit=300
                                    ! should match c_portio.h
  INTEGER, PARAMETER            :: minUnit=10
                                    !to exclude special purpose units

  CHARACTER (LEN=*),                                                           &
      PARAMETER       :: ioNoFileName ='* No Filename available *'
  CHARACTER (LEN=*),                                                           &
      PARAMETER       :: ioNoEnvName  ='* No Environment available *'

! Type for containing file/unit state
  TYPE unitAttributes
     INTEGER :: state
     INTEGER :: policy
     INTEGER :: nextLocation
     INTEGER :: fileType
     LOGICAL :: from_environment
     CHARACTER (LEN=maxFileNameLength) :: fileName
     CHARACTER (LEN=maxFileNameLength) :: env
  END TYPE unitAttributes

! Type for reporting the state of a file, e.g. for enquiry functions
  TYPE, PUBLIC :: fileState
     SEQUENCE
     INTEGER :: unit
     INTEGER :: fileExtent
     INTEGER :: filePosition
  END TYPE fileState

! Length of the above type for transmitting it via e.g. MPI
  INTEGER, PARAMETER :: fileState_len=3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Intefaces and access
!
! Nomenclature for implementations:
!
! buff[in|out]<wordsize>_<array_type><rank>

!read numWords words (of some type) from disk
  INTERFACE buffin
    MODULE PROCEDURE                                                           &
        buffin64_r,                                                            &
        buffin64_i,                                                            &
        buffin64_l,                                                            &
        buffin64_r2D,                                                          &
        buffin64_i2D,                                                          &
        buffin64_l2D,                                                          &
        buffin_character,                                                      &
        buffin32_r,                                                            &
        buffin32_i,                                                            &
        buffin32_i2D,                                                          &
        buffin32_l
  END INTERFACE

!read numWords 64 bit words from disk
  INTERFACE buffin64
    MODULE PROCEDURE                                                           &
        buffin64_r,                                                            &
        buffin64_i,                                                            &
        buffin64_l,                                                            &
        buffin64_r2D,                                                          &
        buffin64_i2D,                                                          &
        buffin64_l2D,                                                          &
        buffin_character
  END INTERFACE

!read numWords 32 bit words from disk
  INTERFACE buffin32
    MODULE PROCEDURE                                                           &
        buffin32_r,                                                            &
        buffin32_i,                                                            &
        buffin32_i2D,                                                          &
        buffin32_l
  END INTERFACE

!write numWords words (of some type) to disk
  INTERFACE buffout
    MODULE PROCEDURE                                                           &
        buffout64_r,                                                           &
        buffout64_i,                                                           &
        buffout64_l,                                                           &
        buffout64_r2D,                                                         &
        buffout64_i2D,                                                         &
        buffout64_l2D,                                                         &
        buffout32_i,                                                           &
        buffout32_r,                                                           &
        buffout32_i2D
  END INTERFACE

!write numWords 64 bit words to disk
  INTERFACE buffout64
    MODULE PROCEDURE                                                           &
        buffout64_r,                                                           &
        buffout64_i,                                                           &
        buffout64_l,                                                           &
        buffout64_r2D,                                                         &
        buffout64_i2D,                                                         &
        buffout64_l2D
  END INTERFACE

!write numWords 32 bit words to disk
  INTERFACE buffo32
    MODULE PROCEDURE                                                           &
        buffout32_i,                                                           &
        buffout32_r,                                                           &
        buffout32_i2D
  END INTERFACE

! Access control
  PRIVATE performOperation
  PRIVATE setAttribute
  PRIVATE clearAttribute
  PRIVATE testAttribute
  PUBLIC  ioInit
  PUBLIC  ioShutdown
  PUBLIC  io_ereport
  PUBLIC  ioDiskSynchronise
  PUBLIC  ioFence
  PUBLIC  setpos
  PUBLIC  setpos8
  PUBLIC  setpos32
  PUBLIC  set_unit_bcast_flag
  PUBLIC  clear_unit_bcast_flag
  PUBLIC  buffin
  PUBLIC  buffout
  PUBLIC  buffo32
  PUBLIC  buffin32
  PUBLIC  buffout64
  PUBLIC  buffin64
  PUBLIC  is_unit_open
  PUBLIC  isRemote
  PUBLIC  file_open
  PUBLIC  file_close
  PUBLIC  io_timestep
  PUBLIC  ioFileState
  PUBLIC  readableWords
  PUBLIC  ReturnFileName
  PUBLIC  isReadOnly

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Global variables (access is not thread safe)

! File States
  TYPE(unitAttributes)               :: fileAttributes(minUnit:maxUnit)

! General purpose message
  CHARACTER (LEN=2*maxFileNameLength):: io_message

! params/vars  for dr_hook
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CONTAINS


!+ Top level initialisation,
!+    load config
!+    initialise C
!+    initialise this module
!+
  SUBROUTINE ioInit()
    USE io_configuration_mod

    IMPLICIT NONE
    LOGICAL, SAVE  :: initted=.FALSE.
    INTEGER        :: unit
    REAL(KIND=jprb):: zhook_handle

    IF (lhook) CALL dr_hook('IO:IOINIT',zhook_in,zhook_handle)

    IF (.NOT.initted) THEN
      CALL ioLoadConfig()
      CALL portioInit(         & 
          io_wbuffer_size,     & 
          io_rbuffer_size,     & 
          io_rbuffer_count,    & 
          io_rbuffer_update,   & 
          io_rbuffer_prefetch, & 
          io_timing)


! Initialise IO Timer
      IF (io_timing==io_timing_on) THEN
        CALL io_timing_init()
      END IF


      DO unit=minUnit,maxUnit
        CALL resetFileAttributes(unit)
      END DO
      IF (printstatus>=prstatus_normal)                                        &
          WRITE(6,'(A)')'IO: Initialised IO'
      initted=.TRUE.

    END IF

    IF (lhook) CALL dr_hook('IO:IOINIT',zhook_out,zhook_handle)

  END SUBROUTINE ioInit

!+ Close down io activity, tidy up and maybe print timing data etc.
!+
  SUBROUTINE ioShutdown()
    USE io_configuration_mod, ONLY :                                           &
        io_timing,                                                             &
        io_timing_on
    IMPLICIT NONE
    REAL(KIND=jprb):: zhook_handle

    IF (lhook) CALL dr_hook('IO:IOSHUTDOWN',zhook_in,zhook_handle)


    ! Write out PE 0 IO timings if required, for 2B this is handled
    ! in the portioshutdown call sequence.
    IF (( mype == 0)  .OR. L_IO_Server ) THEN
      IF (io_timing==io_timing_on) THEN
        WRITE(6,'(A,I5)') 'IO timings for PE ', mype
        CALL io_total_timings()
      END IF
    END IF


    CALL portioShutdown()
    IF (printstatus>=prstatus_diag)                                            &
        WRITE(6,'(A)')'IO: Shutdown IO'

    IF (lhook) CALL dr_hook('IO:IOSHUTDOWN',zhook_out,zhook_handle)

  END SUBROUTINE ioShutdown

!+ Internal subroutine to reset the unit state description to defaults
  SUBROUTINE resetFileAttributes(unit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit ! Fortran unit
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:RESETFILEATTS',zhook_in,zhook_handle)
    fileAttributes(unit)%fileName(:)     =' '
    fileAttributes(unit)%env(:)          =' '
    fileAttributes(unit)%fileName        =ioNoFileName
    fileAttributes(unit)%env             =ioNoEnvName
    fileAttributes(unit)%from_environment=.FALSE.
    fileAttributes(unit)%policy          =ioNotSet
    fileAttributes(unit)%state           =ioNotSet
    fileAttributes(unit)%nextLocation    =ioNoLocation
    fileAttributes(unit)%fileType        =ioFileTypeUnknown
!$OMP FLUSH(fileAttributes)
    IF (lhook) CALL dr_hook('IO:RESETFILEATTS',                                &
        zhook_out,zhook_handle)
  END SUBROUTINE resetFileAttributes


!+ Internal subroutine to wrap ereport and provide more io specific
!  system state
  SUBROUTINE io_ereport(r,code,m,unit)
    USE ereport_mod, ONLY : ereport
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: r
    CHARACTER(LEN=*),INTENT(IN) :: m
    INTEGER, INTENT(IN)         :: code ! code to pass through to ereport
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL                :: unit ! unit responsible for the error
    INTEGER                     :: lcode
    INTEGER                     :: i
    REAL(KIND=jprb)             :: zhook_handle

    IF (lhook) CALL dr_hook('IO:EREPORT',zhook_in,zhook_handle)
    lcode=code ! fix ereport_mod/ereport API mismatch
    WRITE(6,'(A,A)')                                                           &
        '****************** IO Error Report ',                                 &
        '***********************************'
    IF (PRESENT(unit))THEN
      WRITE(6,'(A,I5)')'Unit Generating error=',unit
    END IF
    IF (printstatus>=prstatus_diag.OR.code>0) THEN
      WRITE(6,'(A)')'---File States --------------------------'
      DO i=minUnit,maxUnit
        CALL ioStatusPrint(i)
      END DO
      WRITE(6,'(A)')'---End File States ----------------------'
    ELSE
      WRITE(6,'(A,A)')                                                         &
          '*** File states can be reported by ',                               &
          'setting diagnostic output levels **'
    END IF
    CALL ereport(r,lcode,m)

    IF (lhook) CALL dr_hook('IO:EREPORT',zhook_out,zhook_handle)
  END SUBROUTINE io_ereport


!+ Internal subroutine to report the state of a unit to stdout
  SUBROUTINE ioStatusPrint(unit)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit ! unit to report about
    INTEGER             :: state
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:STATUSPRINT',zhook_in,zhook_handle)

    state=fileAttributes(unit)%state
    IF (testAttribute(state,ioOpen)) THEN
      WRITE(6,'(A,i3,A,A)')'Unit ',unit,' open on filename ',                  &
          TRIM(fileAttributes(unit)%filename)

      IF (fileAttributes(unit)%from_environment) THEN
        WRITE(6,'(A,A)')'  --> Opened from environment variable:',             &
            TRIM(fileAttributes(unit)%env)
      END IF

      WRITE(6,'(A,I3,A,L1,A,L1)') &
          '  --> File Type: ',fileAttributes(unit)%fileType,      &
          ' , Read Only: ',testAttribute(state,ioReadOnly),                    &
          ' , Write Only: ',testAttribute(state,ioWriteOnly)
      WRITE(6,'(A,L1,A,L1,A,L1,A,L1)')                                         &
          '  --> Local: ',testAttribute(state,ioLocal),                        &
          ' AllLocal: ',testAttribute(state,ioAllLocal),                       &
          ' Remote: ',testAttribute(state,ioRemote),                           &
          ' Broadcast: ',testAttribute(state,ioReadReplicate)
    END IF

    IF (lhook) CALL dr_hook('IO:STATUSPRINT',zhook_out,zhook_handle)

  END SUBROUTINE ioStatusPrint


!+ Internal Subroutine to check if the calling process has a valid
!+ unit number.
  SUBROUTINE checkUnit(unit)

    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: unit
    INTEGER               :: abortCode=10

    IF ( unit > maxUnit .OR. unit < minUnit ) THEN
      WRITE(6,'(A,I10,A,i3,A,i3)')                                             &
          'Fortran unit out of allowed range, unit = ',unit,                   &
          ' RANGE = ',minUnit,' - ',maxUnit
      CALL io_ereport ('io:checkUnit',abortCode,                               &
          'Fortran Unit provided is out of range',unit=unit)
    END IF

  END SUBROUTINE checkUnit


!+ Internal function to check if the calling process
!+ will participate in the operation. This enables us
!+ to be able to check parameters sanely
  LOGICAL FUNCTION performOperation(unit,reading,writing)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit    !Fortran unit
    LOGICAL, INTENT(IN),                                                       &
        OPTIONAL        :: reading !Are we reading from the file
    LOGICAL, INTENT(IN),                                                       &
        OPTIONAL        :: writing !Are we writing to the file
    LOGICAL             :: Lreading
    LOGICAL             :: Lwriting
    INTEGER             :: abortCode=11

    Lreading=.FALSE.
    Lwriting=.FALSE.
    IF (PRESENT(reading))Lreading=reading
    IF (PRESENT(writing))Lwriting=writing

    CALL CheckUnit(unit)

    performOperation=.FALSE.
    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      performOperation=.TRUE.
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)                 &
        .AND. mype==0) THEN
      performOperation=.TRUE.
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)                &
        .AND. mype==0) THEN
      performOperation=.TRUE.
    END IF

  ! If we are reading we may participate in a broadcast op.
    IF (Lreading.AND.broadcast_read(unit)) THEN
      performOperation=.TRUE.
    END IF

    !If we are going to be involved in an IO in any way,
    !do some saneness checks
    IF (performOperation) THEN
      !the unit should be recognisably open
      IF (.NOT.testAttribute(fileAttributes(unit)%state,ioOpen))               &
          CALL io_ereport('io:performOperation',                               &
          abortCode,'Operation attempted to a closed unit',unit=unit)
      !We shouldn't  write to a 'read only' defined unit
      IF (Lwriting.AND.                                                        &
          testAttribute(fileAttributes(unit)%state,ioReadOnly))                &
          CALL io_ereport('io:performOperation',                               &
          abortCode,'Operation to write to a ReadOnly unit',unit=unit)
      !We cant read from a 'write only' defined unit
      IF (Lreading.AND.                                                        &
          testAttribute(fileAttributes(unit)%state,ioWriteOnly))               &
          CALL io_ereport('io:performOperation',                               &
          abortCode,'Operation to read from a WriteOnly unit',unit=unit)
    END IF
  END FUNCTION performOperation

!+ Report is a unit is remote or not
  LOGICAL FUNCTION isRemote(unit)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit          !File Unit number
    INTEGER             :: status
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:ISREMOTE',zhook_in,zhook_handle)
    CALL checkUnit(unit)

    isRemote=testAttribute(fileAttributes(unit)%state,ioRemote)

    IF (lhook) CALL dr_hook('IO:ISREMOTE',zhook_out,zhook_handle)

    RETURN
  END FUNCTION isRemote

!+ Report is a unit is open or not
  LOGICAL FUNCTION is_unit_open(unit)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit          !File Unit number
    INTEGER             :: status
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:IS_UNIT_OPEN',zhook_in,zhook_handle)
    CALL ioInit()
    CALL checkUnit(unit)

    is_unit_open=.FALSE.
    IF     (testAttribute(fileAttributes(unit)%state,ioOpen)) THEN
      is_unit_open=.TRUE.
    END IF
    IF (lhook) CALL dr_hook('IO:IS_UNIT_OPEN',zhook_out,zhook_handle)

    RETURN
  END FUNCTION is_unit_open


!+ Open a file on the given unit
  SUBROUTINE file_open (                                                       &
      unit,env,env_len,read_write,name_in_environ,error,                       &
      ioLocality,ioPolicy,allowRemap,fileType)
    USE IOS_Common, ONLY : IOS_local_ro_files
    IMPLICIT NONE
    INTEGER, INTENT(IN)               :: unit
                                       !Unit number for I/O
    CHARACTER(LEN=*),INTENT(IN)       :: env
                                       !Environment name or
                                       !explicit file name
    INTEGER, INTENT(IN), OPTIONAL     :: env_len
                                       !Length of env
    INTEGER, INTENT(IN), OPTIONAL     :: read_write
                                       ! == 0 read only access,
                                       ! /= 0 read+write
    INTEGER, INTENT(IN), OPTIONAL     :: name_in_environ
                                       ! == 0 file name stored in
                                       !environment variable, otherwise
                                       !name specified explicitly
    INTEGER, INTENT(IN), OPTIONAL     :: ioPolicy
                                       !Degree of inconsistency allowed
                                       !between file ops before
                                       !generating errors
    INTEGER, INTENT(IN), OPTIONAL     :: ioLocality
                                       !Where IO on this unit should
                                       !take place (e.g. all processors,
                                       ! remote services)
    LOGICAL, INTENT(IN), OPTIONAL     :: allowRemap
                                       !Specifies that this file has a new
                                       !filename and hence can be remapped
                                       !to a different IO server.
    INTEGER, INTENT(OUT), OPTIONAL    :: error
                                       ! == 0 file OPENED
                                       ! /= 0 file NOT OPENED due to error
    INTEGER, INTENT(IN), OPTIONAL     :: fileType

    ! Local copies of optional args
    INTEGER                           :: Local_fileType
    INTEGER                           :: Local_read_write
    INTEGER                           :: Local_nameinenv
    INTEGER                           :: Local_envlen
    INTEGER                           :: Local_Error

    INTEGER                           :: info
    INTEGER                           :: abortCode=12
    INTEGER                           :: warnCode=-12
    INTEGER                           :: fileNameLength
    INTEGER                           :: newPolicy
    INTEGER                           :: newState
    INTEGER                           :: mymin
    INTEGER                           :: mymax
    INTEGER                           :: mymodel_com
    INTEGER                           :: server
    CHARACTER (LEN=maxFileNameLength) :: fileName
    CHARACTER (LEN=30)                :: serveString
    REAL(KIND=jprb)                   :: zhook_handle

    IF (lhook) CALL dr_hook('IO:FILE_OPEN',zhook_in,zhook_handle)
    CALL ioInit()
    CALL checkUnit(unit)

    Local_fileType=ioFileTypeUnknown
    IF (PRESENT(fileType)) THEN
      Local_fileType=fileType
    END IF


    Local_nameinenv=ioNameProvided
    IF (PRESENT(name_in_environ)) THEN
      IF (name_in_environ==ioNameinEnv)                                        &
          Local_nameinenv=ioNameInEnv
    END IF

    Local_read_write = ioOpenReadWrite
    IF (PRESENT(read_write)) THEN
      IF (read_write/= ioOpenReadOnly .AND. read_write /= ioOpenWriteOnly) THEN
        Local_read_write = ioOpenReadWrite
      ELSE
        Local_read_write = read_write
      END IF
    END IF

    IF (PRESENT(env_len)) THEN
      local_envlen=env_len
    ELSE
      local_envlen=LEN_TRIM(env)
    END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, we perform sanity checks
! - Set params for optional arguments
! - adapt arguments to the IO Server configuration
! - trap inconsistent and illegal requests
! - enforce collective semantics
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Check arguments and set default policy
    IF(.NOT.PRESENT(ioPolicy)) THEN
      newPolicy=ioPolicyDefault
    ELSE
      IF(ioPolicy==ioPolicyStrict.OR.                                          &
          ioPolicy==ioPolicySloppy.OR.                                         &
          ioPolicy==ioPolicyLax) THEN
        newPolicy=ioPolicy
      ELSE
        WRITE(6,'(A,i5,A)')'Policy value of ',ioPolicy,' is invalid'
        CALL io_ereport('io:file_open',abortCode,                              &
            'Invalid Optional Policy Argument',unit=unit)
      END IF
    END IF

    ! Check arguments and set default Locality
    IF(.NOT.PRESENT(ioLocality)) THEN
      newState=ioDefault
    ELSE
      newState=ioLocality
    END IF

    ! If there is no IO server and remote is set,
    ! then set the local flag instead
    IF( testAttribute(newState,ioRemote) .AND. .NOT.l_ios_active())THEN
      CALL setAttribute(newState,ioLocal)
      CALL clearAttribute(newState,ioRemote)
      IF (printstatus>=prstatus_diag)                                          &
          WRITE(6,'(A,A)')'IO: Switching file mode to local because ',         &
          'there is no IO server'
    END IF

    ! Check that user doesn't ask for read broadcast and per task IO
    IF (testAttribute(newState,ioAllLocal) .AND.                               &
        testAttribute(newState,ioReadReplicate)) THEN
      CALL io_ereport('io:file_open',abortCode,                                &
          'Requested per task IO and broadcast :-o',unit=unit)
    END IF

    !Route read only files locally if requested.
    IF ( IOS_local_ro_files ) THEN
      IF (local_read_write==0) THEN
        IF( testAttribute(newState,ioRemote)) THEN
          CALL setAttribute(newState,ioLocal)
          CALL clearAttribute(newState,ioRemote)
          IF (printstatus>=prstatus_diag)                                      &
              WRITE(6,'(A)')'IO: Routing read only file to local'
        END IF
      END IF
    END IF

    ! Clear the broadcast bit if there is only one processor, so that we
    ! avoid calling gcom in serial applications
    IF(testAttribute(newState,ioReadReplicate).AND.                            &
        model_procs==1) THEN
      CALL clearAttribute(newState,ioReadReplicate)
      IF (printstatus>=prstatus_diag)                                          &
          WRITE(6,'(A,A)')'IO: Clearing the broadcast read flag ',             &
          'because there is only one CPU'
    END IF

    IF (printstatus>=prstatus_diag) THEN
      IF (testAttribute(newState,ioReadReplicate)) THEN
        WRITE(6,'(A,I3,A)')'IO: Opening unit ',unit,                           &
            ' with collective(broadcast) semantics'
      ELSE
        WRITE(6,'(A,I3,A)')'IO: Opening unit ',unit,                           &
            ' without collective(broadcast) semantics'
      END IF
      IF (local_read_write==0)                                                 &
          WRITE(6,'(A)')'IO: Read Only mode'
    END IF

    Local_Error=0

    ! Check if the unit is already attached
    IF (is_unit_open(unit)) THEN
      IF (newPolicy==ioPolicyStrict .OR.                                       &
          fileAttributes(unit)%policy==ioPolicyStrict) THEN
        CALL io_ereport('io:file_open',abortCode,                              &
            'Unit already open',unit)
      ELSE
        IF (printstatus>=prstatus_diag) THEN
          CALL io_ereport('io:file_open',warnCode,                             &
              'Unit already open - attempting close',unit)
        END IF
        CALL file_close(unit,env,env_len,name_in_environ,ioNoDelete,Local_Error)
      END IF
    END IF


    ! For small jobs (where its cheap) check for consistency
    ! We also avoid this for small execs as currently the IBM
    ! 32 bit builds link with a gcom 64 bit library, so we will skip it
    ! for the nproc=1 case too.

    IF (.NOT.testAttribute(newstate,ioAllLocal).AND.                           &
        nproc<=256.AND.nproc>1) THEN
      IF (printstatus>=prstatus_diag) THEN
        WRITE(6,'(A)')'IO: Checking consistency of unit open request...'
!DEPENDS ON : um_fort_flush
        CALL um_fort_flush(6,Local_Error)
      END IF
      CALL gc_get_communicator(mymodel_com,Local_Error)

      mymin=unit
      mymax=unit
      CALL GC_imin(1,nproc,Local_Error,mymin)
      CALL GC_imax(1,nproc,Local_Error,mymax)
      IF (mymin/=mymax) CALL io_ereport('file_open',abortCode,                 &
          'unit no. mismatch',unit)

      mymin=local_read_write
      mymax=local_read_write
      CALL GC_imin(1,nproc,Local_Error,mymin)
      CALL GC_imax(1,nproc,Local_Error,mymax)
      IF (mymin/=mymax) CALL io_ereport('file_open',abortCode,                 &
          'read_write flag mismatch',unit)

      mymin=Local_nameinenv
      mymax=Local_nameinenv
      CALL GC_imin(1,nproc,Local_Error,mymin)
      CALL GC_imax(1,nproc,Local_Error,mymax)
      IF (mymin/=mymax) CALL io_ereport('file_open',abortCode,                 &
          'name_in_environ mismatch',unit)

      mymin=newstate
      mymax=newstate
      CALL GC_imin(1,nproc,Local_Error,mymin)
      CALL GC_imax(1,nproc,Local_Error,mymax)
      IF (mymin/=mymax) CALL io_ereport('file_open',abortCode,                 &
          'locale mismatch',unit)

      IF (printstatus>=prstatus_diag) THEN
        WRITE(6,'(A)')'IO: Valid request'
!DEPENDS ON : um_fort_flush
        CALL um_fort_flush(6,Local_Error)
      END IF

    END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Second, we recover the filename from the environment if needed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF (Local_nameinenv == ioNameInEnv) THEN
      CALL fort_get_env(env,local_envlen,filename,maxFileNameLength,info)
      fileNameLength=LEN_TRIM(filename)
      IF (info == -1 ) THEN
        WRITE(io_message,'(A,A,A)')'Environment variable ',TRIM(env),          &
            ' not set.'
        CALL io_ereport('io:file_open',abortCode,io_message,                   &
            unit=unit)
      END IF
    ELSE
      IF (local_envlen>maxFileNameLength) THEN
        WRITE(io_message,'(A,I10,A,I3)')'File name too long (',                &
            local_envlen,') must be less than ',maxFileNameLength
        CALL io_ereport('io:file_open',abortCode,                              &
            io_message,unit=unit)
      END IF
      fileNameLength=local_envlen
      filename  =env(1:fileNameLength)
    END IF
    !Blank the end of the filename string so we can use trim()
    ! with impunity
    filename(fileNameLength+1:)=' '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Third, we actually open the file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Local_Error=0
    serveString=""
    IF      (testAttribute(newState,ioAllLocal)) THEN
      CALL open_single(unit,                                                   &
          TRIM(filename),LEN_TRIM(filename),local_read_write,1,Local_Error)
    ELSE IF (testAttribute(newState,ioLocal) .AND. mype==0) THEN
      CALL open_single(unit,                                                   &
          TRIM(filename),LEN_TRIM(filename),local_read_write,1,Local_Error)
    ELSE IF (testAttribute(newState,ioRemote)) THEN
      server=IOS_Open(unit,                                                    &
          TRIM(filename),local_read_write,Local_filetype,allowRemap)
      IF (server /= -1) THEN
        WRITE(serveString,'(A,I4,A)')' (Server = ',server,')'
      ELSE
        WRITE(serveString,'(A,I4,A)')' (ServerUnknown)'
      END IF
    ELSE IF(mype==0) THEN
      WRITE(6,'(A)')'IO: FILE_OPEN: Bad Request.'
      Local_Error=-999
    ELSE
      Local_Error=0
    END IF

    IF (Local_Error/=0)THEN
      IF ( fileAttributes(unit)%policy==ioPolicyStrict) THEN
        CALL io_ereport('io:file_open',abortCode,                              &
            'An error occured opening a file',unit=unit)
      ELSE
        CALL io_ereport('io:file_open',warnCode,                               &
            'An error occured opening a file',unit=unit)
        Local_Error=0
      END IF
    END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fourth, we set the state of the file according to what
! actually happened, or die if we failed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF (Local_Error==0) THEN
      ! File was opened so set any other attributes
      ! even if its not open local to this rank
      IF (Local_nameinenv == ioNameInEnv) THEN
        fileAttributes(unit)%env=env
        fileAttributes(unit)%env(local_envlen+1:)=' '
        fileAttributes(unit)%filename=filename
        fileAttributes(unit)%filename(fileNameLength+1:)=' '
        fileAttributes(unit)%from_environment=.TRUE.
      ELSE
        fileAttributes(unit)%fileName=env
        fileAttributes(unit)%fileName(local_envlen+1:)=' '
        fileAttributes(unit)%env=ioNoEnvName
        fileAttributes(unit)%from_environment=.FALSE.
      END IF
      fileAttributes(unit)%state=newState
      fileAttributes(unit)%policy=newPolicy
      fileAttributes(unit)%fileType=Local_fileType
      CALL setAttribute(fileAttributes(unit)%state,ioOpen)
      IF (local_read_write==ioReadOnly)                                        &
          CALL setAttribute(fileAttributes(unit)%state,ioReadOnly)
      IF (local_read_write==ioWriteOnly)                                       &
          CALL setAttribute(fileAttributes(unit)%state,ioWriteOnly)


      WRITE(6,'(A,A,A,I3,A)')                                                  &
          'IO: Open: ',TRIM(fileAttributes(unit)%fileName),                    &
          ' on unit ',unit,serveString
      IF (Local_nameinenv == ioNameInEnv) THEN
        WRITE(6,'(A,A)')'IO: from environment variable ',                      &
            TRIM(fileAttributes(unit)%env)
      END IF
!$OMP FLUSH(fileAttributes)
    ELSE
      ! We failed.
      IF (Local_nameinenv == ioNameInEnv) THEN
        WRITE(6,'(A,A,A)')                                                     &
            'IO: Failure to open file specified by environment ',              &
            'variable: ',                                                      &
            TRIM(env)
      END IF
      WRITE(6,'(A,A)') 'IO: Failure to open file: ',TRIM(filename)
      ! Reset the file attributes
      CALL io_ereport('io:file_open',error,'Failed to open file',              &
          unit=unit)
      CALL resetFileAttributes(unit)
    END IF

    IF ( PRESENT(error) ) error = Local_Error

    IF (lhook) CALL dr_hook('IO:FILE_OPEN',zhook_out,zhook_handle)

  END SUBROUTINE file_open


!+ Close the file on the given unit
  SUBROUTINE file_close(unit,env,env_len,name_in_environ,delete,error)

    IMPLICIT NONE
    INTEGER, INTENT(IN):: unit           !Unit number for I/O
    CHARACTER(LEN=*), INTENT(IN)                                               &
                  :: env                 !Environment name or filename
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL  :: env_len             !Length of env
    INTEGER, INTENT(IN), &               ! == 0, do not delete file
        OPTIONAL       :: delete         ! /= 0, delete file
    INTEGER, INTENT(IN), &               ! == 0 file name stored in env var
        OPTIONAL  :: name_in_environ     ! /= 0 file name given
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL  :: error               !return code (0=no error)

! Local copies of optional args
    INTEGER            :: Local_Error
    INTEGER            :: Local_nameinenv
    INTEGER            :: Local_delete
    INTEGER            :: Local_envlen

    INTEGER            :: abortCode=33   !for ereport
    INTEGER            :: warnCode=-34   !for ereport
    REAL(KIND=jprb)    :: zhook_handle
    INTEGER                           :: fileNameLength
    CHARACTER (LEN=maxFileNameLength) :: fileName
    INTEGER                           :: info

    IF (lhook) CALL dr_hook('IO:FILE_CLOSE',zhook_in,zhook_handle)

    Local_nameinenv=ioNameProvided
    IF (PRESENT(name_in_environ)) THEN
      IF (name_in_environ==ioNameInEnv)                                        &
          Local_nameinenv=ioNameInEnv
    END IF

    Local_delete=ioNoDelete
    IF (PRESENT(delete)) THEN
      IF (delete/=ioNoDelete)                                                  &
          Local_delete=ioDelete
    END IF

    IF (PRESENT(env_len)) THEN
      local_envlen=env_len
    ELSE
      local_envlen=LEN_TRIM(env)
    END IF



    Local_Error = 0
    CALL checkUnit(unit)

    IF (Local_nameinenv == ioNameInEnv) THEN
      CALL fort_get_env(env,local_envlen,filename,maxFileNameLength,info)
      fileNameLength=LEN_TRIM(filename)
      IF (info == -1 ) THEN
        WRITE(io_message,'(A,A,A)')'Environment variable ',TRIM(env),          &
            ' not set.'
        CALL io_ereport('io:file_close',abortCode,io_message,                  &
            unit=unit)
      END IF
    ELSE
      IF (local_envlen>maxFileNameLength) THEN
        WRITE(io_message,'(A,I10,A,I3)')'File name too long (',                &
            local_envlen,') must be less than ',maxFileNameLength
        CALL io_ereport('io:file_close',abortCode,                             &
            io_message,unit=unit)
      END IF
      fileNameLength=local_envlen
      filename      =env(1:fileNameLength)
    END IF

    !Blank the end of the filename string so we can use trim()
    ! with impunity
    filename(fileNameLength+1:)=' '


    IF (is_unit_open(unit)) THEN
      IF     (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
        CALL close_single(unit,TRIM(filename),LEN_TRIM(filename),1,            &
            Local_delete,Local_Error)
      ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)               &
          .AND. mype==0) THEN
        CALL close_single(unit,TRIM(filename),LEN_TRIM(filename),1,            &
            Local_delete,Local_Error)
      ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)              &
          .AND. mype==0) THEN
        CALL IOS_Close(unit,TRIM(filename),Local_delete)
      END IF

      IF (printstatus>=prstatus_oper) THEN
        WRITE(6,'(A,A,A,I3)')                                                  &
            'IO: Close: ',TRIM(fileAttributes(unit)%fileName),                 &
            ' on unit ',unit
        IF (fileAttributes(unit)%from_environment) THEN
          WRITE(6,'(A,A)')'IO: Originally from environment variable ',         &
              TRIM(fileAttributes(unit)%env)
        END IF
      END IF

      IF (TRIM(fileName) /= TRIM(fileAttributes(unit)%fileName)) THEN
        IF (printstatus>=prstatus_oper .AND. Local_delete /= 0) THEN
          CALL io_ereport('io:file_close ',warncode,                           &
              'Filename to delete does not match with open',                   &
              unit=unit)
        END IF
      END IF
    ELSE
      IF (printstatus>=prstatus_diag) THEN
        CALL io_ereport('io:file_close',warncode,                              &
            'Tried to close a unit which was not open',unit=unit)
      END IF
    END IF ! IS Unit Open?



    CALL resetFileAttributes(unit)

    IF ( PRESENT(error) ) error = Local_Error

    IF (lhook) CALL dr_hook('IO:FILE_CLOSE',zhook_out,zhook_handle)

  END SUBROUTINE file_close


!+ Parallel UM version of buffout
  SUBROUTINE buffout64_r                                                       &
      (unit,array,numWords,wordsWritten,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL        :: wordsWritten    ! Words written
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL        :: numWords        ! No. of words to write out
    REAL, INTENT(OUT),                                                         &
        OPTIONAL        :: errorCode
    REAL(KIND=real64),                                                         &
        INTENT(IN)      :: array(:)        ! Array to write out

! Local copies of opt. args
    INTEGER             :: Local_WordsWritten    ! Words written
    INTEGER             :: Local_NumWords
    REAL                :: Local_Error      ! Exit Code

    INTEGER             :: info
    INTEGER             :: abortCode=13
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsWritten=Local_NumWords

    IF (performOperation(unit,writing=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffout',abortCode,                                  &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffout64_single                                                    &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)                 &
        .AND.mype==0) THEN
      CALL buffout64_single                                                    &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)                &
        .AND.mype==0) THEN
      CALL ios_write64     (unit,fileAttributes(unit)%nextLocation,            &
          array(1:Local_NumWords),Local_NumWords)
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_out,zhook_handle)

  END SUBROUTINE buffout64_r


!+ Output 64 bit data to the given unit
  SUBROUTINE buffout64_i                                                       &
      (unit,array,numWords,wordsWritten,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL        :: wordsWritten    ! Words written
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL        :: numWords        ! no. of words to write out
    REAL, INTENT(OUT),                                                         &
        OPTIONAL        :: errorCode
    INTEGER(KIND=integer64),                                                   &
        INTENT(IN)      :: array(:)

! Local copies of opt. args
    INTEGER             :: Local_NumWords
    REAL                :: Local_Error       ! Exit Code       ! Exit Code
    INTEGER             :: Local_WordsWritten    ! Words written

    INTEGER             :: abortCode=14
    INTEGER             :: info
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsWritten=Local_NumWords

    IF (performOperation(unit,writing=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffout64',abortCode,                                &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffout64_single                                                    &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal).AND.            &
        mype==0) THEN
      CALL buffout64_single                                                    &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF(testAttribute(fileAttributes(unit)%state,ioRemote) .AND.           &
        mype==0) THEN
      CALL ios_write64     (unit,fileAttributes(unit)%nextLocation,            &
          array(1:Local_NumWords),Local_NumWords)
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_out,zhook_handle)

  END SUBROUTINE buffout64_i


!+ Output 64 bit data to the given unit
  SUBROUTINE buffout64_l                                                       &
      (unit,array,numWords,wordsWritten,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL        :: wordsWritten ! Words written
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL        :: numWords     ! no. of words to write out
    REAL, INTENT(OUT),                                                         &
        OPTIONAL        :: errorCode
    LOGICAL(KIND=logical64),                                                   &
        INTENT(IN)      :: array(:)     ! Array to write out

! Local copies of opt. args
    INTEGER             :: Local_WordsWritten
    INTEGER             :: Local_NumWords
    REAL                :: Local_Error

    INTEGER             :: info
    INTEGER             :: abortCode=15
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsWritten=Local_NumWords

    IF (performOperation(unit,writing=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffout64',abortCode,                                &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffout64_single                                                    &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal).AND.            &
        mype==0) THEN
      CALL buffout64_single                                                    &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF(testAttribute(fileAttributes(unit)%state,ioRemote) .AND.           &
        mype==0) THEN
      CALL ios_write64     (unit,fileAttributes(unit)%nextLocation,            &
          array(1:Local_NumWords),Local_NumWords)
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_out,zhook_handle)

  END SUBROUTINE buffout64_l


!+ Output 64 bit data to the given unit
  SUBROUTINE buffout64_r2D                                                     &
      (unit,array,numWords,wordsWritten,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL        :: numWords     ! no. of words to write out
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL        :: wordsWritten ! Words written
    REAL, INTENT(OUT),                                                         &
        OPTIONAL        :: errorCode
    REAL(KIND=real64),                                                         &
        INTENT(IN)      :: array(:,:)   ! Array to write data from

! Local copies of opt. args
    INTEGER             :: Local_NumWords
    REAL                :: Local_Error
    INTEGER             :: Local_WordsWritten

    INTEGER             :: info
    INTEGER             :: abortCode=16
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsWritten=Local_NumWords

    IF (performOperation(unit,writing=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffout64',abortCode,                                &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffout64_single                                                    &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal).AND.            &
        mype==0) THEN
      CALL buffout64_single                                                    &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote) .AND.          &
        mype==0) THEN
      CALL ios_write64     (unit,fileAttributes(unit)%nextLocation,            &
          array(:,:),Local_NumWords)
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_out,zhook_handle)

  END SUBROUTINE buffout64_r2D


!+ Output 64 bit data to the given unit
  SUBROUTINE buffout64_i2D                                                     &
      (unit,array,numWords,wordsWritten,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER(KIND=integer64),                                                   &
        INTENT(IN)      :: array(:,:)   ! Array to write data from
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL        :: numWords     ! no. of words to write out
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL        :: wordsWritten ! Words written
    REAL, INTENT(OUT),                                                         &
        OPTIONAL        :: errorCode

! Local copies of opt. args
    INTEGER             :: Local_NumWords
    REAL                :: Local_Error
    INTEGER             :: Local_WordsWritten

    INTEGER             :: info
    INTEGER             :: abortCode=17
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsWritten=Local_NumWords

    IF (performOperation(unit,writing=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffout64',abortCode,                                &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffout64_single                                                    &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal).AND.            &
        mype==0) THEN
      CALL buffout64_single                                                    &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF(testAttribute(fileAttributes(unit)%state,ioRemote) .AND.           &
        mype==0) THEN
      CALL ios_write64     (unit,fileAttributes(unit)%nextLocation,            &
          array(:,:),Local_NumWords)
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_out,zhook_handle)

  END SUBROUTINE buffout64_i2D


!+ Output 64 bit data to the given unit
  SUBROUTINE buffout64_l2D                                                     &
      (unit,array,numWords,wordsWritten,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    LOGICAL(KIND=logical64),                                                   &
        INTENT(IN)      :: array(:,:)   ! Array to write data from
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL        :: numWords     ! no. of words to write out
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL        :: wordsWritten ! Words written
    REAL, INTENT(OUT),                                                         &
        OPTIONAL        :: errorCode

    ! Local copies of opt. args
    INTEGER             :: Local_NumWords
    REAL                :: Local_Error
    INTEGER             :: Local_WordsWritten

    INTEGER             :: info
    INTEGER             :: abortCode=18
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsWritten=Local_NumWords

    IF (performOperation(unit,writing=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffout64',abortCode,                                &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffout64_single                                                    &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal).AND.            &
        mype==0) THEN
      CALL buffout64_single                                                    &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote) .AND.          &
        mype==0) THEN
      CALL ios_write64     (unit,fileAttributes(unit)%nextLocation,            &
          array(:,:),Local_NumWords)
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_out,zhook_handle)

  END SUBROUTINE buffout64_l2D


!+ Output 32 bit data to the given unit
  SUBROUTINE buffout32_r                                                       &
      (unit,array,numWords,wordsWritten,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL         :: numWords     ! no. of words to write out
    REAL(KIND=real32),                                                         &
        INTENT(IN)       :: array(:)     ! Array to write out
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL        :: wordsWritten ! Words written

! Local copies of opt. args
    REAL                 :: Local_Error
    INTEGER              :: Local_WordsWritten
    INTEGER              :: Local_NumWords

    INTEGER              :: info         ! Broadcast information code
    INTEGER              :: abortCode=19
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsWritten=Local_NumWords

    IF (performOperation(unit,writing=.TRUE.).AND.                             &
        Local_NumWords>2*SIZE(array)) THEN
      WRITE(6,'(A,I10,A)')'Request for ',Local_NumWords,' 32 bit words'
      WRITE(6,'(A,A,I10)')'is not supported by a 64 bit ',                     &
          'real buffer of size ',                                              &
          SIZE(array)
      CALL io_ereport('io:buffo32',abortCode,                                  &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffo32_single                                                      &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal).AND.            &
        mype==0) THEN
      CALL buffo32_single                                                      &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote).AND.           &
        mype==0) THEN
      CALL ios_write32   (unit,fileAttributes(unit)%nextLocation,              &
          array(1:Local_NumWords),Local_NumWords)
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_out,zhook_handle)

  END SUBROUTINE buffout32_r


!+ Output 32 bit data to the given unit
  SUBROUTINE buffout32_i                                                       &
      (unit,array,numWords,wordsWritten,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL         :: numWords     ! no. of words to write out
    INTEGER(KIND=integer32),                                                   &
        INTENT(IN)       :: array(:)     ! Array to write out
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL        :: wordsWritten ! Words written

! Local copies of opt. args
    INTEGER              :: Local_WordsWritten
    INTEGER              :: Local_NumWords
    REAL                 :: Local_Error

    INTEGER              :: info         ! Broadcast information code
    INTEGER              :: abortCode=20
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsWritten=Local_NumWords

    IF (performOperation(unit,writing=.TRUE.).AND.                             &
        Local_NumWords>2*SIZE(array)) THEN
      WRITE(6,'(A,I10,A)')'Request for ',Local_NumWords,' 32 bit words'
      WRITE(6,'(A,A,I10)')'is not supported by a 64 bit int ',                 &
          'buffer of size ',                                                   &
          SIZE(array)
      CALL io_ereport('io:buffo32',abortCode,                                  &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffo32_single                                                      &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal).AND.            &
        mype==0) THEN
      CALL buffo32_single                                                      &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote) .AND.          &
        mype==0) THEN
      CALL ios_write32   (unit,fileAttributes(unit)%nextLocation,              &
          array(1:Local_NumWords),Local_NumWords)
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_out,zhook_handle)

  END SUBROUTINE buffout32_i


!+ Output 32 bit data to the given unit
  SUBROUTINE buffout32_i2D                                                     &
      (unit,array,numWords,wordsWritten,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL         :: numWords     ! no. of words to write out
    INTEGER(KIND=integer32),                                                   &
        INTENT(IN)       :: array(:,:)   ! Array to write out
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL        :: wordsWritten ! Words written

! Local copies of opt. args
    INTEGER              :: Local_WordsWritten
    INTEGER              :: Local_NumWords
    REAL                 :: Local_Error

    INTEGER              :: info         ! Broadcast information code
    INTEGER              :: abortCode=21
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsWritten=Local_NumWords

    IF (performOperation(unit,writing=.TRUE.).AND.                             &
        Local_NumWords>2*SIZE(array)) THEN
      WRITE(6,'(A,I10,A)')'Request for ',Local_NumWords,' 32 bit words'
      WRITE(6,'(A,I10)')'is not supported by a 64 bit 2D int buffer of size ', &
          SIZE(array)
      CALL io_ereport('io:buffo32',abortCode,                                  &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffo32_single                                                      &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal).AND.            &
        mype==0) THEN
      CALL buffo32_single                                                      &
          (unit,array,Local_NumWords,Local_WordsWritten,Local_Error)
    ELSE IF(testAttribute(fileAttributes(unit)%state,ioRemote) .AND.           &
        mype==0) THEN
      CALL ios_write32   (unit,fileAttributes(unit)%nextLocation,              &
          array(:,:),Local_NumWords)
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF ( PRESENT(wordsWritten) ) wordsWritten = Local_WordsWritten
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFOUT',zhook_out,zhook_handle)

  END SUBROUTINE buffout32_i2D


!+ Read 64 bit data from the given unit
  SUBROUTINE buffin64_r(unit,array,numWords,wordsRead,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL         :: numWords     ! no. of words to read in
    REAL(KIND=real64),                                                         &
        INTENT(OUT)      :: array(:)     ! Array to read in
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
    INTEGER              :: Local_NumWords
    INTEGER              :: Local_WordsRead
    REAL                 :: Local_Error

    INTEGER              :: abortCode=22
    INTEGER              :: info         ! Broadcast information code
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsRead=Local_NumWords

    IF (performOperation(unit,reading=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffin',abortCode,                                   &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffin64_single                                                     &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)) THEN
      IF (mype==0)                                                             &
          CALL buffin64_single                                                 &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
      IF(broadcast_read(unit))                                                 &
          CALL gc_rbcast(1,Local_NumWords,0,nproc,info,array)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)) THEN
      IF(broadcast_read(unit)) THEN
        CALL ios_read64     (unit,fileAttributes(unit)%nextLocation,           &
            array(1:Local_NumWords),Local_NumWords,IOS_Read_Broadcast)
      ELSE
        CALL ios_read64     (unit,fileAttributes(unit)%nextLocation,           &
            array(1:Local_NumWords),Local_NumWords,IOS_Read_NoBroadcast)
      END IF
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
      WRITE(io_message,'(A,F4.2,A,I10,A,I10)')                                 &
          'Error in buffin errorCode=',Local_Error,                            &
          ' len=',Local_WordsRead,'/',Local_NumWords
      CALL io_ereport('io:buffin',abortCode,io_message,                        &
          unit=unit)
    END IF

    IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_out,zhook_handle)

  END SUBROUTINE buffin64_r


!+ Read 64 bit data from the given unit
  SUBROUTINE buffin64_i(unit,array,numWords,wordsRead,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL         :: numWords     ! no. of words to read in
    INTEGER(KIND=integer64),                                                   &
        INTENT(OUT)      :: array(:)     ! Array to read in
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL         :: wordsRead    ! No. of words read
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode

! Local copies of opt. args
    INTEGER              :: Local_NumWords
    INTEGER              :: Local_WordsRead
    REAL                 :: Local_Error

    INTEGER              :: info         ! Broadcast information code
    INTEGER              :: abortCode=24
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsRead=Local_NumWords

    IF (performOperation(unit,reading=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffin',abortCode,                                   &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffin64_single                                                     &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)) THEN
      IF (mype==0)                                                             &
          CALL buffin64_single                                                 &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
      IF(broadcast_read(unit))                                                 &
          CALL gc_ibcast(1,Local_NumWords,0,nproc,info,array)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)) THEN
      IF(broadcast_read(unit)) THEN
        CALL ios_read64     (unit,fileAttributes(unit)%nextLocation,           &
            array(1:Local_NumWords),Local_NumWords,IOS_Read_Broadcast)
      ELSE
        CALL ios_read64     (unit,fileAttributes(unit)%nextLocation,           &
            array(1:Local_NumWords),Local_NumWords,IOS_Read_NoBroadcast)
      END IF
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
      WRITE(io_message,'(A,F4.2,A,I10,A,I10)')                                 &
          'Error in buffin errorCode=',Local_Error,                            &
          ' len=',Local_WordsRead,'/',Local_NumWords
      CALL io_ereport('io:buffin',abortCode,io_message,                        &
          unit=unit)
    END IF

    IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_out,zhook_handle)

  END SUBROUTINE buffin64_i


!+ Read 64 bit data from the given unit
  SUBROUTINE buffin64_l(unit,array,numWords,wordsRead,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL         :: numWords     ! no. of words to read in
    LOGICAL(KIND=logical64),                                                   &
        INTENT(OUT)      :: array(:)     ! Array to read into
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
    INTEGER              :: Local_WordsRead
    INTEGER              :: Local_NumWords
    REAL                 :: Local_Error

    INTEGER              :: abortCode=25
    INTEGER              :: info         ! Broadcast information code
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsRead=Local_NumWords

    IF (performOperation(unit,reading=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffin',abortCode,                                   &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffin64_single                                                     &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)) THEN
      IF (mype==0)                                                             &
          CALL buffin64_single                                                 &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
      IF(broadcast_read(unit))                                                 &
          CALL gc_ibcast(1,Local_NumWords,0,nproc,info,array)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)) THEN
      IF(broadcast_read(unit)) THEN
        CALL ios_read64     (unit,fileAttributes(unit)%nextLocation,           &
            array(1:Local_NumWords),Local_NumWords,IOS_Read_Broadcast)
      ELSE
        CALL ios_read64     (unit,fileAttributes(unit)%nextLocation,           &
            array(1:Local_NumWords),Local_NumWords,IOS_Read_NoBroadcast)
      END IF
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
      WRITE(io_message,'(A,F4.2,A,I10,A,I10)')                                 &
          'Error in buffin errorCode=',Local_Error,                            &
          ' len=',Local_WordsRead,'/',Local_NumWords
      CALL io_ereport('io:buffin',abortCode,io_message,                        &
          unit=unit)
    END IF

    IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_out,zhook_handle)

  END SUBROUTINE buffin64_l


!+ Read 64 bit data from the given unit
  SUBROUTINE buffin64_i2D                                                      &
      (unit,array,numWords,wordsRead,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL         :: numWords     ! no. of words to read in
    INTEGER(KIND=integer64),                                                   &
        INTENT(OUT)      :: array(:,:)   ! Array to read into
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
    INTEGER              :: Local_WordsRead
    INTEGER              :: Local_NumWords
    REAL                 :: Local_Error

    INTEGER              :: info         ! Broadcast information code
    INTEGER              :: abortCode=26
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsRead=Local_NumWords

    IF (performOperation(unit,reading=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffin',abortCode,                                   &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffin64_single                                                     &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)) THEN
      IF (mype==0)                                                             &
          CALL buffin64_single                                                 &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
      IF(broadcast_read(unit))                                                 &
          CALL gc_ibcast(1,Local_NumWords,0,nproc,info,array)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)) THEN
      IF(broadcast_read(unit)) THEN
        CALL ios_read64     (unit,fileAttributes(unit)%nextLocation,           &
            array,Local_NumWords,IOS_Read_Broadcast)
      ELSE
        CALL ios_read64     (unit,fileAttributes(unit)%nextLocation,           &
            array,Local_NumWords,IOS_Read_NoBroadcast)
      END IF
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
      WRITE(io_message,'(A,F4.2,A,I10,A,I10)')                                 &
          'Error in buffin errorCode=',Local_Error,                            &
          ' len=',Local_WordsRead,'/',Local_NumWords
      CALL io_ereport('io:buffin',abortCode,io_message,                        &
          unit=unit)
    END IF

    IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_out,zhook_handle)

  END SUBROUTINE buffin64_i2D


!+ Read 64 bit data from the given unit
  SUBROUTINE buffin64_l2D                                                      &
      (unit,array,numWords,wordsRead,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL         :: numWords     ! no. of words to read in
    LOGICAL(KIND=logical64),                                                   &
        INTENT(OUT)      :: array(:,:)   ! Array to read into
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
    INTEGER              :: Local_WordsRead
    INTEGER              :: Local_NumWords
    REAL                 :: Local_Error

    INTEGER              :: info         ! Broadcast information code
    INTEGER              :: abortCode=27
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsRead=Local_NumWords

    IF (performOperation(unit,reading=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffin',abortCode,                                   &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffin64_single                                                     &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)) THEN
      IF (mype==0)                                                             &
          CALL buffin64_single                                                 &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
      IF(broadcast_read(unit))                                                 &
          CALL gc_ibcast(1,Local_NumWords,0,nproc,info,array)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)) THEN
      IF(broadcast_read(unit)) THEN
        CALL ios_read64     (unit,fileAttributes(unit)%nextLocation,           &
            array,Local_NumWords,IOS_Read_Broadcast)
      ELSE
        CALL ios_read64     (unit,fileAttributes(unit)%nextLocation,           &
            array,Local_NumWords,IOS_Read_NoBroadcast)
      END IF
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
      WRITE(io_message,'(A,F4.2,A,I10,A,I10)')                                 &
          'Error in buffin errorCode=',Local_Error,                            &
          ' len=',Local_WordsRead,'/',Local_NumWords
      CALL io_ereport('io:buffin',abortCode,io_message,                        &
          unit=unit)
    END IF

    IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_out,zhook_handle)

  END SUBROUTINE buffin64_l2D


!+ Read 64 bit data from the given unit
  SUBROUTINE buffin64_r2D                                                      &
      (unit,array,numWords,wordsRead,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL         :: numWords     ! no. of words to read in
    REAL(KIND=real64),                                                         &
        INTENT(OUT)      :: array(:,:)   ! Array to read into
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
    INTEGER              :: Local_WordsRead
    INTEGER              :: Local_NumWords
    REAL                 :: Local_Error

    INTEGER              :: info         ! Broadcast information code
    INTEGER              :: abortCode=28
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsRead=Local_NumWords

    IF (performOperation(unit,reading=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffin',abortCode,                                   &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffin64_single                                                     &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)) THEN
      IF (mype==0)                                                             &
          CALL buffin64_single                                                 &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
      IF(broadcast_read(unit))                                                 &
          CALL gc_rbcast(1,Local_NumWords,0,nproc,info,array)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)) THEN
      IF(broadcast_read(unit)) THEN
        CALL ios_read64     (unit,fileAttributes(unit)%nextLocation,           &
            array,Local_NumWords,IOS_Read_Broadcast)
      ELSE
        CALL ios_read64     (unit,fileAttributes(unit)%nextLocation,           &
            array,Local_NumWords,IOS_Read_NoBroadcast)
      END IF
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
      WRITE(io_message,'(A,F4.2,A,I10,A,I10)')                                 &
          'Error in buffin errorCode=',Local_Error,                            &
          ' len=',Local_WordsRead,'/',Local_NumWords
      CALL io_ereport('io:buffin',abortCode,io_message,                        &
          unit=unit)
    END IF

    IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_out,zhook_handle)

  END SUBROUTINE buffin64_r2D


!+ Read 32 bit data from the given unit
!
  SUBROUTINE buffin32_r(unit,array,numWords,wordsRead,errorCode)

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL         :: numWords     ! no. of words to read in
    REAL(KIND=real32),                                                         &
        INTENT(OUT)      :: array(:)     ! Array to read into
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
    INTEGER              :: Local_WordsRead
    INTEGER              :: Local_NumWords
    REAL                 :: Local_Error

    INTEGER              :: info         ! Broadcast information code
    INTEGER              :: abortCode=29
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsRead=Local_NumWords

    IF (performOperation(unit,reading=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffin32',abortCode,                                 &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffin32_single                                                     &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)) THEN
      IF (mype==0)                                                             &
          CALL buffin32_single                                                 &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
      IF(broadcast_read(unit))                                                 &
          CALL gc_rbcast(1,Local_NumWords/2,0,nproc,info,array)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)) THEN
      IF(broadcast_read(unit)) THEN
        CALL ios_read32     (unit,fileAttributes(unit)%nextLocation,           &
            array(1:Local_NumWords),Local_NumWords,IOS_Read_Broadcast)
      ELSE
        CALL ios_read32     (unit,fileAttributes(unit)%nextLocation,           &
            array(1:Local_NumWords),Local_NumWords,IOS_Read_NoBroadcast)
      END IF
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
      WRITE(io_message,'(A,F4.2,A,I10,A,I10)')                                 &
          'Error in buffin errorCode=',Local_Error,                            &
          ' len=',Local_WordsRead,'/',Local_NumWords
      CALL io_ereport('io:buffin32',abortCode,io_message,                      &
          unit=unit)
    END IF

    IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_out,zhook_handle)

  END SUBROUTINE buffin32_r


!+ Read 32 bit data from the given unit
  SUBROUTINE buffin32_i(unit,array,numWords,wordsRead,errorCode)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL         :: numWords     ! no. of words to read in
    INTEGER(KIND=integer32),                                                   &
        INTENT(OUT)      :: array(:)     ! Array to read into
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
    INTEGER              :: Local_WordsRead
    INTEGER              :: Local_NumWords
    REAL                 :: Local_Error

    INTEGER              :: info         ! Broadcast information code
    INTEGER              :: abortCode=29
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsRead=Local_NumWords

    IF (performOperation(unit,reading=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffin32',abortCode,                                 &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffin32_single                                                     &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)) THEN
      IF (mype==0)                                                             &
          CALL buffin32_single                                                 &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
      IF(broadcast_read(unit))                                                 &
          CALL gc_ibcast(1,Local_NumWords/2,0,nproc,info,array)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)) THEN
      IF(broadcast_read(unit)) THEN
        CALL ios_read32     (unit,fileAttributes(unit)%nextLocation,           &
            array(1:Local_NumWords),Local_NumWords,IOS_Read_Broadcast)
      ELSE
        CALL ios_read32     (unit,fileAttributes(unit)%nextLocation,           &
            array(1:Local_NumWords),Local_NumWords,IOS_Read_NoBroadcast)
      END IF
      fileAttributes(unit)%nextLocation=ioNoLocation
    ELSE
      WRITE(io_message,'(A,I10)')                                              &
          'Error in buffin - unable to resolve filestate ',                    &
          fileAttributes(unit)%state
      CALL io_ereport('io:buffin32',abortCode,io_message,                      &
          unit=unit)
    END IF

    IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
      WRITE(io_message,'(A,F4.2,A,I10,A,I10)')                                 &
          'Error in buffin errorCode=',Local_Error,                            &
          ' len=',Local_WordsRead,'/',Local_NumWords
      CALL io_ereport('io:buffin32',abortCode,io_message,                      &
          unit=unit)
    END IF

    IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_out,zhook_handle)

  END SUBROUTINE buffin32_i


!+ Read 32 bit data from the given unit
  SUBROUTINE buffin32_i2D(unit,array,numWords,wordsRead,errorCode)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL         :: numWords     ! no. of words to read in
    INTEGER(KIND=integer32),                                                   &
        INTENT(OUT)      :: array(:,:)   ! Array to read into
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
    INTEGER              :: Local_WordsRead
    INTEGER              :: Local_NumWords
    REAL                 :: Local_Error

    INTEGER              :: info         ! Broadcast information code
    INTEGER              :: abortCode=29
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsRead=Local_NumWords

    IF (performOperation(unit,reading=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffin32',abortCode,                                 &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffin32_single                                                     &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)) THEN
      IF (mype==0)                                                             &
          CALL buffin32_single                                                 &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
      IF(broadcast_read(unit))                                                 &
          CALL gc_ibcast(1,Local_NumWords/2,0,nproc,info,array)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)) THEN
      IF(broadcast_read(unit)) THEN
        CALL ios_read32     (unit,fileAttributes(unit)%nextLocation,           &
            array,Local_NumWords,IOS_Read_Broadcast)
      ELSE
        CALL ios_read32     (unit,fileAttributes(unit)%nextLocation,           &
            array,Local_NumWords,IOS_Read_NoBroadcast)
      END IF
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
      WRITE(io_message,'(A,F4.2,A,I10,A,I10)')                                 &
          'Error in buffin errorCode=',Local_Error,                            &
          ' len=',Local_WordsRead,'/',Local_NumWords
      CALL io_ereport('io:buffin32',abortCode,io_message,                      &
          unit=unit)
    END IF

    IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_out,zhook_handle)

  END SUBROUTINE buffin32_i2D


!+ Read 32 bit data from the given unit
  SUBROUTINE buffin32_l(unit,array,numWords,wordsRead,errorCode)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: unit         ! Fortran unit number
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL         :: numWords     ! no. of words to read in
    LOGICAL(KIND=logical32),                                                   &
        INTENT(OUT)      :: array(:)     ! Array to read into
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL         :: wordsRead    ! No. of words read

! Local copies of opt. args
    INTEGER              :: Local_WordsRead
    INTEGER              :: Local_NumWords
    REAL                 :: Local_Error

    INTEGER              :: info         ! Broadcast information code
    INTEGER              :: abortCode=29
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_in,zhook_handle)

    IF (PRESENT(numWords)) THEN
      Local_NumWords=numWords
    ELSE
      Local_NumWords=SIZE(array)
    END IF
    Local_Error=-1.0
    Local_WordsRead=Local_NumWords

    IF (performOperation(unit,reading=.TRUE.).AND.                             &
        Local_NumWords>SIZE(array)) THEN
      WRITE(6,'(A,I10,A,I10)')'Request for ',Local_NumWords,                   &
          ' words, is not supported by a buffer of size ',SIZE(array)
      CALL io_ereport('io:buffin32',abortCode,                                 &
          'Supplied buffer too small',unit=unit)
    END IF

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffin32_single                                                     &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)) THEN
      IF (mype==0)                                                             &
          CALL buffin32_single                                                 &
          (unit,array,Local_NumWords,Local_WordsRead,Local_Error)
      IF(broadcast_read(unit))                                                 &
          CALL gc_ibcast(1,Local_NumWords/2,0,nproc,info,array)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)) THEN
      IF(broadcast_read(unit)) THEN
        CALL ios_read32     (unit,fileAttributes(unit)%nextLocation,           &
            array(1:Local_NumWords),Local_NumWords,IOS_Read_Broadcast)
      ELSE
        CALL ios_read32     (unit,fileAttributes(unit)%nextLocation,           &
            array(1:Local_NumWords),Local_NumWords,IOS_Read_NoBroadcast)
      END IF
      fileAttributes(unit)%nextLocation=ioNoLocation
    END IF

    IF (Local_Error /= -1 .OR. Local_WordsRead /= Local_NumWords ) THEN
      WRITE(io_message,'(A,F4.2,A,I10,A,I10)')                                 &
          'Error in buffin errorCode=',Local_Error,                            &
          ' len=',Local_WordsRead,'/',Local_NumWords
      CALL io_ereport('io:buffin32',abortCode,io_message,                      &
          unit=unit)
    END IF

    IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_out,zhook_handle)

  END SUBROUTINE buffin32_l


!+ Read character data from the given unit
  SUBROUTINE buffin_character                                                  &
      (unit,array,numWords,csize,wordsRead,errorCode)

    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: unit            ! Fortran unit number
    INTEGER,  INTENT(IN) :: numWords        ! no. of elements to read
    INTEGER,  INTENT(IN) :: csize           ! Size of an element
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL         :: wordsRead    ! No. of words read
    REAL, INTENT(OUT),                                                         &
        OPTIONAL         :: errorCode
    CHARACTER(LEN=csize),                                                      &
        INTENT(OUT)      :: array(numWords) ! Array to read in

! Local copies of opt. args
    INTEGER              :: Local_WordsRead ! local copy
    REAL                 :: Local_Error     ! local copy

    INTEGER              :: info            ! Broadcast information code
    INTEGER              :: abortCode=23
    REAL(KIND=jprb)      :: zhook_handle

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_in,zhook_handle)

    Local_Error=-1.0
    Local_WordsRead=numWords*csize

    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL buffin8_single(unit,array,numWords*                                 &
          csize,Local_WordsRead,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)) THEN
      IF (mype==0)                                                             &
          CALL buffin8_single(unit,array,numWords*                             &
          csize,Local_WordsRead,Local_Error)
      IF(broadcast_read(unit))                                                 &
          CALL gc_cbcast(1,numWords*csize,0,nproc,info,array)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote).AND.           &
        mype==0) THEN
      CALL io_ereport('io:buffin_character',abortCode,                         &
          'Remote access via buffin8 not supported',unit=unit)
    END IF

    IF ( PRESENT(wordsRead) ) wordsRead = Local_WordsRead
    IF ( PRESENT(errorCode) ) errorCode = Local_Error

    IF (lhook) CALL dr_hook('IO:BUFFIN',zhook_out,zhook_handle)

  END SUBROUTINE buffin_character


!+ Set the file pointer to a position in the file
!  (a 64 bit word address)
  SUBROUTINE setpos(unit,ipos,error)

    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: unit   ! Fortran unit number
    INTEGER, INTENT(IN)     :: ipos   ! Position in file (word aligned)
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL            :: error
    INTEGER                 :: Local_Error  ! Return code
    INTEGER                 :: abortCode=30
    REAL(KIND=jprb)         :: zhook_handle

    IF (lhook) CALL dr_hook('IO:SETPOS',zhook_in,zhook_handle)

    Local_Error = 0
    CALL checkUnit(unit)
    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL setpos_single(unit,ipos,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal).AND.            &
        mype==0) THEN
      CALL setpos_single(unit,ipos,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote).AND.           &
        mype==0) THEN
      fileAttributes(unit)%nextLocation=ipos
!      CALL IOS_Setpos   (unit,ipos,Local_Error)
    END IF

    IF ( PRESENT(error) ) error = Local_Error

    IF (lhook) CALL dr_hook('IO:SETPOS',zhook_out,zhook_handle)

  END SUBROUTINE setpos


!+ Set the file pointer to a position in the file
!  (an 8 bit word address)
  SUBROUTINE setpos8(unit,ipos,error)
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: unit   ! Fortran unit number
    INTEGER, INTENT(IN)     :: ipos   ! Position in file (word aligned)
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL            :: error
    INTEGER                 :: Local_Error  ! Return code
    INTEGER                 :: abortCode=30
    REAL(KIND=jprb)         :: zhook_handle

    IF (lhook) CALL dr_hook('IO:SETPOS8',zhook_in,zhook_handle)

    Local_Error = 0
    CALL checkUnit(unit)
    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL setpos8_single(unit,ipos,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal).AND.            &
        mype==0) THEN
      CALL setpos8_single(unit,ipos,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote).AND.           &
        mype==0) THEN
      CALL io_ereport('io:setpos8',abortcode,                                  &
          'setpos8 not implemented for IO Servers!')
!      fileAttributes(unit)%nextLocation=ipos
!      CALL IOS_Setpos   (unit,ipos,Local_Error)
    END IF

    IF ( PRESENT(error) ) error = Local_Error

    IF (lhook) CALL dr_hook('IO:SETPOS8',zhook_out,zhook_handle)

  END SUBROUTINE setpos8


!+ Set the file pointer to a position in the file
!  (a 32 bit word address)
  SUBROUTINE setpos32(unit,ipos,error)
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: unit   ! Fortran unit number
    INTEGER, INTENT(IN)     :: ipos   ! Position in file (word aligned)
    INTEGER, INTENT(OUT),                                                      &
        OPTIONAL            :: error
    INTEGER                 :: abortCode=30
    INTEGER                 :: Local_Error
    REAL(KIND=jprb)         :: zhook_handle

    IF (lhook) CALL dr_hook('IO:SETPOS32',zhook_in,zhook_handle)

    Local_Error = 0
    CALL checkUnit(unit)
    IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      CALL setpos32_single(unit,ipos,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal).AND.            &
        mype==0) THEN
      CALL setpos32_single(unit,ipos,Local_Error)
    ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote).AND.           &
        mype==0) THEN
      CALL io_ereport('io:setpos32',abortcode,                                 &
          'setpos32 not implemented for IO Servers!')
!      fileAttributes(unit)%nextLocation=ipos
!      CALL IOS_Setpos   (unit,ipos,Local_Error)
    END IF

    IF ( PRESENT(error) ) error = Local_Error

    IF (lhook) CALL dr_hook('IO:SETPOS32',zhook_out,zhook_handle)

  END SUBROUTINE setpos32


!+ Set a bit in a vaiable
  SUBROUTINE setAttribute(var,mask)

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: var  ! Variable
    INTEGER, INTENT(IN)    :: mask ! Bitmask to set
    var=IOR(var,mask)
  END SUBROUTINE setAttribute


!+ Clear a bit in a variable
  SUBROUTINE clearAttribute(var,mask)

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: var  ! Variable
    INTEGER, INTENT(IN)    :: mask ! Bitmask to clear
    INTEGER, PARAMETER     :: minus_one=-1
    INTEGER                :: tmp
    tmp=IEOR(minus_one,mask)
    var=IAND(var,tmp)
  END SUBROUTINE clearAttribute


!+ Test whether a bit is set in a variable
  LOGICAL FUNCTION testAttribute(var,mask)

    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: var  ! Variable
    INTEGER, INTENT(IN)    :: mask ! Bitmask to compare
    testAttribute=.FALSE.
    IF (IAND(var,mask)==mask) testAttribute=.TRUE.
    RETURN
  END FUNCTION testAttribute


!+ Cause a unit to NOT broadcast read data to all tasks
  SUBROUTINE set_unit_bcast_flag(unit)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    CALL checkUnit(unit)
    CALL clearAttribute(fileAttributes(unit)%state,ioReadReplicate)
  END SUBROUTINE set_unit_bcast_flag


!+ Cause a unit to broadcast read data to all tasks
  SUBROUTINE clear_unit_bcast_flag(unit)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    CALL checkUnit(unit)
    CALL setAttribute(fileAttributes(unit)%state,ioReadReplicate)
  END SUBROUTINE clear_unit_bcast_flag


!+ Enquire as to whether a unit will broadcast read data
  LOGICAL FUNCTION broadcast_read(unit)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    CALL checkUnit(unit)
    broadcast_read=                                                            &
        testAttribute(fileAttributes(unit)%state,ioReadReplicate)
    RETURN
  END FUNCTION broadcast_read


!+ Pass a timestep event to the subsystem
!   (important for coupled models)
  SUBROUTINE io_timestep()
    ! Timestep notification

    IMPLICIT NONE
    INTEGER, SAVE   :: ts=0
    REAL(KIND=jprb) :: zhook_handle

    IF (lhook) CALL dr_hook('IO:TS',zhook_in,zhook_handle)
    ts=ts+1
    IF (l_ios_active().AND.mype==0) CALL IOS_Process(ts)
    IF (lhook) CALL dr_hook('IO:TS',zhook_out,zhook_handle)
  END SUBROUTINE io_timestep


!+ Commit data to the disk subsystem with optional wait on completion
!   (important for coupled models)
  SUBROUTINE ioDiskSynchronise(unit,wait)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit      ! Unit to sync to disk
    LOGICAL, OPTIONAL   :: wait      ! If true, the caller should block until
                                     ! the opertaion comletes. If the unit
                                     ! is open locally then this will have
                                     ! little effect, but remote IO units
                                     ! setting this flag will see a significant

    INTEGER             :: icode     ! Local error code
    LOGICAL             :: L_wait    ! local copy of optional argument
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:DISKSYNC',zhook_in,zhook_handle)

    L_Wait=.FALSE.
    IF (PRESENT(wait)) THEN
      IF (wait)L_wait=wait
    END IF

! As this can be expensive, lets have some documentation that it happened
    IF (printstatus>=prstatus_oper) THEN
      WRITE(6,'(A,I3,A)')'IO: Synchronising unit: ',unit,' with disk.'
      IF (l_wait) THEN
        WRITE(6,'(A,I3)')'IO: Disk sync blocking until completion...'
      END IF
    END IF

    IF (testAttribute(fileAttributes(unit)%state,ioOpen)) THEN
      IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
        CALL sync_single(unit,icode)             ! Flush 1 layers
      ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)               &
          .AND.mype==0) THEN
        CALL sync_single(unit,icode)
      ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)              &
          .AND.mype==0) THEN
        CALL IOS_DiskSync(unit,l_wait)
      END IF
    END IF

    IF (lhook) CALL dr_hook('IO:DISKSYNC',zhook_out,zhook_handle)

  END SUBROUTINE ioDiskSynchronise


!+ Generate a fence operation. Define a model epoch such as to make guarantees
!+ about IO operation ordering
!+
!+ Calling this is equivilent to the statement "Nothing further will happen
!+ regarding unit 'unit', until all other units instructions issued before
!+ point this have been completed". The epoch is described to have passed
!+ when this condition is satisfied.
!+
!+ Implementation notes:
!+
!+  Operations issued for other units after the epoch may still complete
!+  before the epoch completes.
!+
!+  Atmosphere does not know that the epoch has completed.
!+
!+  If 'unit' is open locally then Atmosphere will stall waiting for the
!+  epoch to pass. Therefore it is not recommended to use this method
!+  on a local file if IO servers are active and processing other IO channels.
!+
!+  If 'unit' is remote then the IO server responsible will stall. This may
!+  cause performance problems (ie a stall in Atmosphere) if his queue becomes
!+  full, and cannot receive new work from the model.

  SUBROUTINE ioFence(unit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER             :: root
    INTEGER             :: errorcode=41
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:FENCE',zhook_in,zhook_handle)
    IF (testAttribute(fileAttributes(unit)%state,ioOpen)) THEN
      IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
        CALL io_ereport('io:Fence',errorcode,                                  &
            'Cannot perform fence operation on a unit with many writers')
      ELSE IF  (mype==0) THEN
        IF (testAttribute(fileAttributes(unit)%state,ioLocal)) THEN
          root=mype
        ELSE
          root=-1
        END IF
        CALL IOS_Fence(unit,root)
      END IF

    END IF
    IF (lhook) CALL dr_hook('IO:FENCE',zhook_out,zhook_handle)
  END SUBROUTINE ioFence

  INTEGER FUNCTION readableWords(unit,wordLength)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER, OPTIONAL   :: wordLength
    INTEGER             :: local_wordLength
    TYPE(fileState)     :: state
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:READABLEWORDS',                                &
        zhook_in,zhook_handle)

! If not specified we'll assume they want default type lengths
    local_wordLength=umFortranRealSize()

    IF ( PRESENT(wordLength) ) local_wordLength=wordLength

    CALL ioFileState(unit,state)
    readableWords=(state%fileExtent-state%filePosition)/                       &
        local_wordLength

    IF (lhook) CALL dr_hook('IO:READABLEWORDS',                                &
        zhook_out,zhook_handle)

  END FUNCTION readableWords

  SUBROUTINE ioFileState(unit,state)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER             :: info ! generic error code
    TYPE(fileState)     :: state
    TYPE(IOS_Status)    :: rstate
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IO:FILESTATE',zhook_in,zhook_handle)

    state%unit         = unit
    state%fileExtent   = ioNotKnown
    state%filePosition = ioNotKnown
    info               = 0

    IF (testAttribute(fileAttributes(unit)%state,ioOpen)) THEN
      IF      (testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
        CALL getpos8  (unit,state%filePosition)
        CALL getextent(unit,state%fileExtent,info)
      ELSE IF (testAttribute(fileAttributes(unit)%state,ioLocal)               &
          .AND.mype==0) THEN
        CALL getpos8  (unit,state%filePosition)
        CALL getextent(unit,state%fileExtent,info)
      ELSE IF (testAttribute(fileAttributes(unit)%state,ioRemote)              &
          .AND.mype==0) THEN
        CALL IOS_fileState(unit,fileAttributes(unit)%nextLocation,rstate)
        state%filePosition =rstate%position
        state%fileExtent   =rstate%extent
      END IF
    END IF

    IF (.NOT.testAttribute(fileAttributes(unit)%state,ioAllLocal)) THEN
      IF(broadcast_read(unit)) THEN
        CALL gc_ibcast(1,filestate_len,0,nproc,info,state)
      END IF
    END IF

    IF (lhook) CALL dr_hook('IO:FILESTATE',zhook_out,zhook_handle)

  END SUBROUTINE ioFileState

  CHARACTER(LEN=maxFileNameLength) FUNCTION ReturnFileName(unit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit

    ReturnFileName=fileAttributes(unit)%filename

  END FUNCTION ReturnFileName

  LOGICAL FUNCTION isReadOnly(unit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit

    isReadOnly=testAttribute(fileAttributes(unit)%state,ioReadOnly)

  END FUNCTION isReadOnly

END MODULE io
