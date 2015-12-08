! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! A simple data module containing variables related to STASH buffering
!
! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C84
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

MODULE model_file

  USE yomhook,    ONLY : lhook, dr_hook
  USE parkind1,   ONLY : jprb, jpim
  USE ios, ONLY :                                                              &
      IOS_Init_PP_Lookup,                                                      &
      IOS_Set_Header,                                                          &
      IOS_Init_Header
  USE IOS_Communicators, ONLY :                                                &
      model_rank, model_procs
  USE IOS_Stash_Common, ONLY :                                                 &
      L_IO_Server,                                                             &
      isUsingAsyncStash
  USE io
  USE ereport_mod, ONLY : ereport
  USE UM_ParVars
  USE printstatus_mod
  USE lookup_addresses
  USE filenamelength_mod, ONLY : fileNameLength                                

  USE MPPIO_file_utils, ONLY : file_touch, file_delete

  USE io_configuration_mod, ONLY : io_external_control

  IMPLICIT NONE
  PRIVATE
! Profiling
  INTEGER(KIND=jpim),PRIVATE, PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim),PRIVATE, PARAMETER :: zhook_out = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Locations in the pp file for interesting items, otherwise
! see f3.ps from the umdp
!
! File addresses
  INTEGER, PARAMETER, PUBLIC   :: MF_Lookup_Address=150
  INTEGER, PARAMETER, PUBLIC   :: MF_Data_Address=160
! Other Data
  INTEGER, PARAMETER, PUBLIC   :: MF_Num_Preallocated_Headers=152
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CHARACTER (LEN=80),PRIVATE           :: Cmessage ! General purpose string

! This value is part of the file format, and other utilities expect it
  INTEGER, PARAMETER, PUBLIC   :: MF_Data_Missing     =-99
  INTEGER, PARAMETER, PUBLIC   :: MF_Unit_Max         = 179  ! Max unit number
  INTEGER, PARAMETER, PUBLIC   :: MF_Unit_Min         = 10   ! Min unit number

! Parts of the file
  INTEGER, PARAMETER, PUBLIC   :: FixedHeader         = 1
  INTEGER, PARAMETER, PUBLIC   :: IntConsts           = 2
  INTEGER, PARAMETER, PUBLIC   :: RealConsts          = 3
  INTEGER, PARAMETER, PUBLIC   :: RowConsts           = 4
  INTEGER, PARAMETER, PUBLIC   :: LevConsts           = 5
  INTEGER, PARAMETER, PUBLIC   :: ColConsts           = 6

! States of parts of the file
  INTEGER, PARAMETER           :: stUnset             = 1
  INTEGER, PARAMETER           :: stUpdFromModel      = 2
  INTEGER, PARAMETER           :: stUpdFromIOS        = 3
  INTEGER, PARAMETER           :: stPendingFromIOS    = 4

  INTEGER, PARAMETER   :: fixedHeaderLen      = 256

  TYPE mf_data_t
     LOGICAL           :: managed = .FALSE.
     INTEGER           :: mf_lookup_address
     INTEGER, POINTER  :: fixed_header(:)     =>NULL()!Fixed Length Header
     INTEGER, POINTER  :: mf_lookup_table(:,:)=>NULL()
  END TYPE mf_data_t

! Static buffers for PP data - 1 per permissible unit
  TYPE(mf_data_t),SAVE :: mf_data( mf_unit_min : mf_unit_max )

  PUBLIC model_file_init
  PUBLIC model_file_open
  PUBLIC model_file_close
  PUBLIC model_file_managed
  PUBLIC get_file_address
  PUBLIC get_mf_information
  PUBLIC InitHeader
  PUBLIC attachHeader
  PUBLIC SetHeader
  PUBLIC loadHeader
  PUBLIC storeHeader
  PUBLIC InitLookups
  PUBLIC attachLookups
  PUBLIC SetLookups
  PUBLIC StoreLookups
  PUBLIC StoreAllLookups
  PUBLIC SynchroniseAll
  PUBLIC setRecordDiskLength
  PUBLIC setRecordDataLength
  PUBLIC setRecordDiskStart
! pass through:
  PUBLIC ioFileTypeUM

! Semantics:
!
! Init           Initialise memory and reset, relays to IOS
! attach         Provide a pointer to the memory
! Set            Inform that data has changed, possibly providing data
!                relays to IOS
! Get            Recover data to user buffer
! load           load from disk
! store          save to disk

CONTAINS


! Determine whether an operation has a relay to IO server action associated
  LOGICAL FUNCTION shouldPropogate(unit)
    IMPLICIT NONE
    INTEGER, INTENT(IN)         :: unit
    shouldPropogate=.FALSE.

    IF (model_rank==0.AND.isRemote(unit)) shouldPropogate=.TRUE.

  END FUNCTION shouldPropogate


  SUBROUTINE mf_printstate(unit)
    IMPLICIT NONE
    INTEGER, INTENT(IN)         :: unit
    WRITE(6,'(A,I4,A)')'Managed unit ',unit,' *******************'
    WRITE(6,'(A,L1)')' mf_lookup_address = ',                                  &
        mf_data(unit)%mf_lookup_address
    IF (ASSOCIATED(mf_data(unit)%fixed_header)) THEN
      WRITE(6,'(A,I8)')' fixed_header(size)',                                  &
          SIZE(mf_data(unit)%fixed_header)
    ELSE
      WRITE(6,'(A)')' fixed_header not associated'
    END IF
    IF (ASSOCIATED(mf_data(unit)%mf_lookup_table)) THEN
      WRITE(6,'(A,I8)')' mf_lookup_table(size)',                               &
          SIZE(mf_data(unit)%mf_lookup_table)
    ELSE
      WRITE(6,'(A)')' mf_lookup_table not associated'
    END IF

  END SUBROUTINE mf_printstate

  SUBROUTINE mf_ereport(r,code,m,unit)
    USE io, ONLY : io_ereport
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN) :: r
    CHARACTER(LEN=*),INTENT(IN) :: m
    INTEGER, INTENT(IN)         :: code ! code to pass through to ereport
    INTEGER, INTENT(IN),                                                       &
        OPTIONAL                :: unit ! unit responsible for the error
    INTEGER                     :: lcode
    INTEGER                     :: i
    REAL(KIND=jprb)             :: zhook_handle

    IF (lhook) CALL dr_hook('MODEL_FILE:EREPORT',zhook_in,zhook_handle)
    lcode=code ! fix ereport_mod/ereport API mismatch
    WRITE(6,'(A,A)')                                                           &
        '****************** Model_File Error Report ',                         &
        '***************************'
    IF (PRESENT(unit))THEN
      WRITE(6,'(A,I5)')'Unit Generating error=',unit
!      IF (printstatus>=prstatus_diag.OR.code>0) THEN
      CALL mf_printstate(unit)
      DO i=mf_unit_min,mf_unit_max
        IF (mf_data(i)%managed) THEN
          CALL mf_printstate(i)
        END IF
      END DO

!      ELSE
!        WRITE(6,'(A,A)')&
!            '*** File states can be reported by ',&
!            'setting diagnostic output levels **'
!      END IF
    END IF
    CALL io_ereport(r,lcode,m,unit)

    IF (lhook) CALL dr_hook('MODEL_FILE:EREPORT',zhook_out,zhook_handle)
  END SUBROUTINE mf_ereport


!----------------------------------------------------------------------
! Function: get_file_address
! Get the address of an object for a pp file
!----------------------------------------------------------------------
  INTEGER FUNCTION get_file_address(unit,item)

    IMPLICIT NONE
    INTEGER,INTENT(IN) :: unit
    INTEGER,INTENT(IN) :: item
    CHARACTER (LEN=*), PARAMETER :: SubroutineName='MODEL_FILE:GET_FILE_ADDRESS'

    CALL checkUnit(unit,SubroutineName)

    IF (.NOT.ASSOCIATED(mf_data(unit)%fixed_header)) THEN
      CALL mf_ereport(SubroutineName, 99,                                      &
          'Request for pp information before initalisation',unit)
    END IF

    SELECT CASE(item)
    CASE (MF_Lookup_Address,MF_Data_Address)
      get_file_address=mf_data(unit)%fixed_header(item)-1
    CASE DEFAULT
      CALL mf_ereport(SubroutineName, 99,                                      &
          'Request for pp file address of an unknown type',unit)
    END SELECT

  END FUNCTION get_file_address


!----------------------------------------------------------------------
! Function: get_mf_information
! Provide information about a pp file
!----------------------------------------------------------------------
  INTEGER FUNCTION get_mf_information(unit,item)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(IN) :: item
    CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
        'MODEL_FILE:GET_MF_INFORMATION'

    CALL checkUnit(unit,SubroutineName)

    IF (.NOT.ASSOCIATED(mf_data(unit)%fixed_header)) THEN

      CALL mf_ereport(SubroutineName, 99,                                      &
          'Request for pp information before initalisation',unit)
    END IF

    IF ( item == MF_Num_Preallocated_Headers ) THEN
      get_mf_information=mf_data(unit)%fixed_header(item)
    ELSE

      CALL mf_ereport(SubroutineName, 99,                                      &
          'Request for pp information of an unknown type',unit)
    END IF

  END FUNCTION get_mf_information


!----------------------------------------------------------------------
! Subroutine: CheckUnit
! Sanity: ensure the unit is in the valid range
!----------------------------------------------------------------------
  SUBROUTINE CheckUnit(Unit,src)

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: Unit
    CHARACTER(LEN=*), INTENT(IN) :: src
    REAL(KIND=jprb)              :: zhook_handle
    INTEGER                      :: errCode
    CHARACTER (LEN=*), PARAMETER :: SubroutineName='MODEL_FILE:CHECKUNIT'

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)

    IF ( Unit < mf_unit_min .OR. Unit > mf_unit_max ) THEN
      WRITE(Cmessage,'(A,i10,A,i3,A,i3,A)')'Unit ',Unit,                       &
          ' supplied to a PP op was outside the allowed range [',              &
          mf_unit_min,'-',mf_unit_max,']'
      CALL mf_ereport(src//' via '//SubroutineName, 99,Cmessage,unit)
    END IF

    IF (.NOT. mf_data(Unit) % managed  ) THEN
      errCode=99
      WRITE(Cmessage,'(A,i4,A)')'Unit ',Unit,                                  &
          ' supplied to model_file is not managed by model_file'
      CALL mf_ereport(src//' via '//SubroutineName, errCode,Cmessage,unit)
    END IF

    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END SUBROUTINE CheckUnit

!----------------------------------------------------------------------
! Subroutine: Model_File_Managed
! Sanity: Report whether the unit is owned by model_file
!----------------------------------------------------------------------
  FUNCTION Model_File_Managed(Unit) result (r)

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: Unit
    LOGICAL                      :: r
    REAL(KIND=jprb)              :: zhook_handle
    CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
        'MODEL_FILE:MODEL_FILE_MANAGED'

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)

    IF ( Unit < mf_unit_min .OR. Unit > mf_unit_max ) THEN
      WRITE(Cmessage,'(A,i10,A,i3,A,i3,A)')'Unit ',Unit,                       &
          ' supplied to a PP op was outside the allowed range [',              &
          mf_unit_min,'-',mf_unit_max,']'
      CALL mf_ereport(SubroutineName, 99,Cmessage,unit)
    END IF

    r=.TRUE.
    IF ( .NOT.mf_data(unit) % managed ) THEN
      r=.FALSE.
    END IF

    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END FUNCTION Model_File_Managed


!----------------------------------------------------------------------
! Subroutine InitHeader
! Initialises fixed header structures
!----------------------------------------------------------------------
  SUBROUTINE InitHeader(unit, header, headerLength)

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: unit
    INTEGER, INTENT(IN)          :: header
    INTEGER, INTENT(IN), OPTIONAL:: headerLength
    INTEGER                      :: local_len
    REAL(KIND=jprb)              :: zhook_handle
    CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
        'MODEL_FILE:INITHEADER'

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)
    CALL checkUnit(unit,SubroutineName)

    IF (PRESENT(headerLength)) THEN
      local_len=headerLength
      IF (header==FixedHeader) THEN
        WRITE(6,'(A)')'Cannot mandate the size of Fixed Header in UM files'
        local_len=fixedHeaderLen
      END IF
    ELSE
      local_len=fixedHeaderLen
      IF (header/=FixedHeader) THEN
        CALL mf_ereport(SubroutineName,99,                                     &
            'Must specify the header length via the headLength argument')
      END IF
    END IF

    IF (shouldPropogate(unit)) CALL IOS_Init_Header(unit,local_len)

! Rank 0 in the atmos model needs to cause initialisation of the header
! on the remote server

    ALLOCATE(mf_data(unit) % fixed_header(local_len))
    mf_data(unit) % fixed_header(:)=MF_Data_Missing

    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END SUBROUTINE InitHeader


!----------------------------------------------------------------------
! Subroutine SetHeader
! Sets the content of a fixed header from an input array
!----------------------------------------------------------------------
  SUBROUTINE SetHeader(unit,header,data_in)

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: unit
    INTEGER, INTENT(IN)          :: header
    INTEGER, OPTIONAL            :: data_in(:)
    INTEGER, POINTER             :: fixhd_internal(:)
    REAL(KIND=jprb)              :: zhook_handle
    CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
        'MODEL_FILE:SETHEADER'

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)
    NULLIFY(fixhd_internal)
    CALL checkUnit(unit,SubroutineName)
    fixhd_internal=>attachHeader(unit,FixedHeader)

    IF (PRESENT(data_in)) THEN

      IF ( SIZE(fixhd_internal) /= SIZE(data_in) ) THEN
        WRITE(6,'(A,I3,A,I12,A,I12)')'ERROR: UNIT=',unit,                      &
            ' LEN INTERNAL=',                                                  &
            SIZE(fixhd_internal),' LEN INPUT=',                                &
            SIZE(data_in)

        CALL mf_ereport(SubroutineName, 99, 'bad parameters in setHeader',unit)
      END IF

      fixhd_internal(:)=data_in(:)

    END IF

    CALL IOS_Set_Header(unit,fixhd_internal)

    NULLIFY(fixhd_internal)
    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END SUBROUTINE SetHeader

  SUBROUTINE loadHeader(unit,header,data_out,icode)

    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: unit
    INTEGER, INTENT(IN)            :: header
    INTEGER, OPTIONAL, INTENT(OUT) :: data_out(:)
    INTEGER, OPTIONAL              :: icode
    INTEGER, POINTER               :: fixhd_internal(:)
    INTEGER                        :: mver_maj
    INTEGER                        :: mver_min
    REAL(KIND=jprb)                :: zhook_handle
    CHARACTER (LEN=*), PARAMETER   :: SubroutineName=                          &
        'MODEL_FILE:LOADHEADER'

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)
    IF (PRESENT(icode)) icode=0

    IF (.NOT. ASSOCIATED (mf_data(unit) % fixed_header))                       &
        CALL initHeader(unit,header)

    CALL setpos(unit,0)

    IF (PRESENT(icode)) THEN
      icode=0
      IF (readableWords(unit)<SIZE(mf_data(unit) % fixed_header)) THEN
        icode=-1
        WRITE(6,'(A,I3,A)')'loadHeader: The file on unit ',unit,               &
            ' does not have enough content to load the header'
        IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)
        RETURN
      END IF
    END IF

    CALL buffin(unit,mf_data(unit) % fixed_header)

    mver_maj=mf_data(unit) % fixed_header(12) / 100
    mver_min=mf_data(unit) % fixed_header(12) - 100 * mver_maj
    WRITE(6,'(A,I1,A,I2)')'loadHeader: Model Version:',mver_maj,'.',mver_min

    IF (PRESENT(data_out)) THEN
      data_out(:)=mf_data(unit) % fixed_header(:)
    END IF

    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END SUBROUTINE loadHeader

!----------------------------------------------------------------------
! SUBROUTINE storeHeader
! Not yet implemented
!----------------------------------------------------------------------

  SUBROUTINE storeHeader(unit,header)

    IMPLICIT NONE
    INTEGER, INTENT(IN)            :: unit
    INTEGER, INTENT(IN)            :: header
    REAL(KIND=jprb)                :: zhook_handle
    CHARACTER (LEN=*), PARAMETER   :: SubroutineName=                          &
        'MODEL_FILE:STOREHEADER'

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)
    CALL mf_ereport(SubroutineName,99,'Not Implemented',unit)
    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END SUBROUTINE storeHeader

!----------------------------------------------------------------------
! FUNCTION attachHeader
! Attaches a header to a given unit
!----------------------------------------------------------------------
  FUNCTION attachHeader(unit,header) result(hdr)

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: unit
    INTEGER, INTENT(IN)          :: header
    INTEGER, POINTER             :: hdr(:)
    REAL(KIND=jprb)              :: zhook_handle
    CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
        'MODEL_FILE:ATTACHHEADER'

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)
    NULLIFY(hdr)
    CALL checkUnit(unit,SubroutineName)


    IF (.NOT.ASSOCIATED(mf_data(unit) % fixed_header))                         &
        CALL mf_ereport(SubroutineName,99,                                     &
        'Request for fixed header which is not allocated yet',unit)

    hdr => mf_data(unit) % fixed_header
    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END FUNCTION attachHeader


!----------------------------------------------------------------------
! SUBROUTINE InitLookups
! Initialises lookup table structures
!----------------------------------------------------------------------
  SUBROUTINE InitLookups(ipplook, unit, len1,len2, address, step)

    IMPLICIT NONE
    INTEGER, POINTER                :: ipplook(:,:)
    INTEGER, INTENT(IN)             :: unit
    INTEGER, INTENT(IN)             :: step
    INTEGER, INTENT(IN)             :: len1
    INTEGER, INTENT(IN)             :: len2
    INTEGER, INTENT(IN)             :: address
    INTEGER                         :: ii, jj
    INTEGER                         :: icode
    REAL(KIND=jprb)                 :: zhook_handle
    INTEGER                         :: input_args(5)
    INTEGER, PARAMETER              :: current_io_pe = 0
    CHARACTER (LEN=*), PARAMETER    :: SubroutineName=                         &
        'MODEL_FILE:INITLOOKUPS'


    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)

    NULLIFY(ipplook)
    CALL checkUnit(unit,SubroutineName)

! Rank 0 in the atmos model needs to cause initialisation of the lookup
! table on the remote server
    IF ( shouldPropogate(unit) ) THEN
      input_args( 1 )=unit
      input_args( 2 )=len1
      input_args( 3 )=len2
      input_args( 4 )=address
      input_args( 5 )=step
      CALL IOS_Init_PP_Lookup(input_args)
    END IF

    IF ( step == 1 ) THEN

      IF (ASSOCIATED(mf_data(unit) % mf_lookup_table))                         &
          DEALLOCATE(mf_data(unit) % mf_lookup_table)

      IF ( mype == current_io_pe .OR. l_io_server ) THEN
        ALLOCATE(mf_data(unit) % mf_lookup_table(len1,len2))
      ELSE
        ! Save a bit of memory, we don't need lookups on
        ! general ranks - but we do want the table associated
        ALLOCATE(mf_data(unit) % mf_lookup_table(1,1))
      END IF

      ipplook => mf_data(unit) % mf_lookup_table
      ipplook(:,:) = MF_Data_Missing

    ELSE IF ( step == 2 ) THEN

      IF ( mype == current_io_pe .OR. l_io_server )                            &
          mf_data(unit) % mf_lookup_address = address

    ELSE   ! step /= 1 or 2
      Cmessage = 'InitLookups step should be 1 or 2'
      Icode = 10

      CALL mf_ereport(SubroutineName, Icode, Cmessage,unit)
    END IF ! Test on step

    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END SUBROUTINE InitLookups


!----------------------------------------------------------------------
! FUNCTION attachLookups
!
! This routine attaches a buffered lookup table to a pointer
! If we are using an IO server with A I O we send a corresponding
! instruction to the IO server too
!----------------------------------------------------------------------
  FUNCTION attachLookups(unit) result(lookup)

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: unit
    INTEGER, POINTER             :: lookup(:,:)
    REAL(KIND=jprb)              :: zhook_handle
    CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
        'MODEL_FILE:ATTACHLOOKUPS'

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)
    CALL checkUnit(unit,SubroutineName)

    NULLIFY (lookup)

    IF (ASSOCIATED(mf_data(unit) % mf_lookup_table)) THEN
      lookup => mf_data(unit) % mf_lookup_table
    ELSE
      WRITE(6,'(A)')'DEBUG: attachLookups: doesnt exist for unit ',unit
    END IF

    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END FUNCTION attachLookups

!----------------------------------------------------------------------
! Subroutine SetRecordDiskLength
! Set the length of a field on the lookup table
! (meaning the disk length, the bytes written to disk including any
!  padding)
!----------------------------------------------------------------------
  SUBROUTINE setRecordDiskLength(unit,RECORD,SIZE)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Unit
    INTEGER, INTENT(IN) :: RECORD
    INTEGER, INTENT(IN) :: SIZE
    INTEGER, POINTER    :: lookup(:,:) => NULL()
    lookup=>attachLookups(unit)
    lookup(lbnrec,RECORD )=SIZE
    NULLIFY(lookup)
  END SUBROUTINE setRecordDiskLength

!----------------------------------------------------------------------
! Subroutine SetRecordDataLength
! Set the length of a field on the lookup table
! (meaning the data length, the bytes written to disk excluding any
!  padding)
!----------------------------------------------------------------------
  SUBROUTINE setRecordDataLength(unit,RECORD,SIZE)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Unit
    INTEGER, INTENT(IN) :: RECORD
    INTEGER, INTENT(IN) :: SIZE
    INTEGER, POINTER    :: lookup(:,:) => NULL()
    lookup=>attachLookups(unit)
    lookup(lblrec,RECORD )=SIZE
    NULLIFY(lookup)

  END SUBROUTINE setRecordDataLength

!----------------------------------------------------------------------
! Subroutine SetDiskRecordStart
! Set the start address of a field in the lookup table
! (meaning the data start address of a field)
!----------------------------------------------------------------------
  SUBROUTINE setRecordDiskStart(unit,RECORD,location)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Unit
    INTEGER, INTENT(IN) :: RECORD
    INTEGER, INTENT(IN) :: Location
    INTEGER, POINTER    :: lookup(:,:) => NULL()

    lookup=>attachLookups(unit)
    lookup(lbegin,RECORD )=location
    NULLIFY(lookup)

  END SUBROUTINE setRecordDiskStart

!----------------------------------------------------------------------
! Subroutine Merge Lookup
! This routine merges a supplied lookup table into a stored one
! We overwrite all records with the supplied data except where
! the value is IDMI - in which case the original value is retained
!
! This is intended to be used by an IO server task to cope with
! situations where either the client or the server has data for a given
! field - e.g. data length, record length and location in the file, which
! may vary depending upon e.g. who did the compression.
!----------------------------------------------------------------------
  SUBROUTINE SetLookups(unit,data)

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: unit
    INTEGER, POINTER             :: data(:)
    INTEGER, POINTER             :: ipplook_in(:,:)
    INTEGER, POINTER             :: ipplook   (:,:)
    INTEGER                      :: i1,i2,ii,icode
    REAL(KIND=jprb)              :: zhook_handle
    CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
        'MODEL_FILE:SETLOOKUPS'

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)
    NULLIFY(ipplook_in,ipplook)
! Fetch the internal lookup
    ipplook=>attachLookups(unit)

! Check that the data input matches the dimention of the one we have
    IF ( SIZE(data) /= SIZE(ipplook) ) THEN
      WRITE(cmessage,'(A,I12,A,I12)')'supplied data size(',SIZE(data),         &
          ') mismatch with internal size (',SIZE(ipplook),')'
      icode=99

      CALL mf_ereport(SubroutineName, Icode, Cmessage,unit)
    END IF

!Allocate a 2D buffer to expand the 1D input buffer
    ALLOCATE( ipplook_in (SIZE(ipplook,1),SIZE(ipplook,2) ) )

!Expand to 2D to make f90 happy
    ii=1
    DO i2=1,SIZE(ipplook,2)
      DO i1=1,SIZE(ipplook,1)
        ipplook_in(i1,i2)=data(ii)
        ii=ii+1
      END DO
    END DO

!Merge the records
    DO i2=1,SIZE(ipplook,2)
      DO i1=1,SIZE(ipplook,1)
        IF ( ipplook(i1,i2) == MF_Data_Missing ) THEN
          ipplook(i1,i2)=ipplook_in(i1,i2)
        END IF
      END DO
    END DO
    DEALLOCATE( ipplook_in )
    NULLIFY(ipplook_in,ipplook)
    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END SUBROUTINE SetLookups

!----------------------------------------------------------------------
! Subroutine: StoreLookups
! Write the lookup table from memory into the file at
! the correct address, having checked it is open, and a
! pp file
!----------------------------------------------------------------------
  SUBROUTINE StoreLookups(unit)

    USE io, ONLY :                                                             &
        setpos,                                                                &
        buffout,                                                               &
        is_unit_open
    USE ios, ONLY :                                                            &
        IOS_Write_Lookup
    USE stwork_aux
    IMPLICIT NONE
    INTEGER,INTENT(IN)           :: unit
    INTEGER                      :: lookupAddress
    INTEGER                      :: icode
    INTEGER                      :: len_io
    REAL                         :: iostat
    REAL(KIND=jprb)              :: zhook_handle
    CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
        'MODEL_FILE:STORELOOKUPS'

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)
    CALL checkUnit(unit,SubroutineName)

    IF (.NOT. L_IO_Server) THEN
      ! Only do the flushing if we actually are storing a lookup table.
      IF (isActivePPFile(unit)) THEN
        lookupAddress=mf_data(unit) % mf_lookup_address
        IF (isUsingAsyncStash()) THEN

          ! Make sure all client side stash data (ie fields) are flushed
          ! to the IO server
          CALL flushPendingStashForUnit(unit)

          ! And then then send the lookup metadata after it
          CALL IOS_Write_Lookup                                       &
              (unit, lookupAddress, mf_data(unit) % mf_lookup_table)

        ELSE
          CALL setpos(unit,lookupAddress,icode)
          CALL buffout(unit,mf_data(unit) % mf_lookup_table,                   &
              SIZE (mf_data(unit) % mf_lookup_table ),                         &
              len_io,iostat)
        END IF
      END IF
    END IF

    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END SUBROUTINE StoreLookups


!----------------------------------------------------------------------
! Subroutine: storeAllLookups
! cycle over all units that potentially have cached lookups
! and flush to disk
!----------------------------------------------------------------------
  SUBROUTINE storeAllLookups()

    IMPLICIT NONE
    INTEGER                      :: i
    REAL(KIND=jprb)              :: zhook_handle
    CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
        'MODEL_FILE:STOREALLLOOKUPS'

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)
    DO i = mf_unit_min,mf_unit_max
      IF (model_file_managed(i)) CALL storeLookups(i)
    END DO
    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END SUBROUTINE storeAllLookups


!----------------------------------------------------------------------
! Subroutine: SynchroniseAll
! cycle over all units that potentially have cached lookups
! and flush to disk
!----------------------------------------------------------------------
  SUBROUTINE SynchroniseAll()

    USE io, ONLY :                                                             &
        ioDiskSynchronise,                                                     &
        is_unit_open
    IMPLICIT NONE
    INTEGER                      :: i
    REAL(KIND=jprb)              :: zhook_handle
    CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
        'MODEL_FILE:SYNCHRONISEALL'

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)
    DO i = mf_unit_min,mf_unit_max
      IF (is_unit_open(i)) THEN
        IF (ASSOCIATED(mf_data(i) % mf_lookup_table)) THEN
          CALL ioDiskSynchronise(i)
        END IF
      END IF
    END DO
    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END SUBROUTINE SynchroniseAll


!----------------------------------------------------------------------
! Subroutine: model_file_init
! init (ie nullify) pp data buffers
!----------------------------------------------------------------------
  SUBROUTINE model_file_init()

    IMPLICIT NONE
    INTEGER                      :: i
    REAL(KIND=jprb)              :: zhook_handle
    CHARACTER (LEN=*), PARAMETER :: SubroutineName=                            &
        'MODEL_FILE:MODEL_FILE_INIT'

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)

    DO i = mf_unit_min,mf_unit_max
      mf_data(i) % managed = .FALSE.
      mf_data(i) % mf_lookup_address=mf_data_missing
      NULLIFY (mf_data(i) % fixed_header )
      NULLIFY (mf_data(i) % mf_lookup_table )
    END DO

    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END SUBROUTINE model_file_init


!----------------------------------------------------------------------
! Function: isActivePPFile
! convenience function return whether a unit has cached lookups
! associated with it
!----------------------------------------------------------------------
  LOGICAL FUNCTION isActivePPFile(unit)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit

    isActivePPFile=.FALSE.
    IF (ASSOCIATED(mf_data(unit) % mf_lookup_table))                           &
        isActivePPFile=.TRUE.

  END FUNCTION isActivePPFile

  SUBROUTINE model_file_open (                                                 &
      unit,env,env_len,read_write,name_in_environ,error,                       &
      ioLocality,ioPolicy,allowRemap,fileType)
    USE IOS_Common, ONLY : IOS_local_ro_files

    IMPLICIT NONE
    INTEGER, INTENT(IN)               :: unit
    INTEGER, INTENT(IN), OPTIONAL     :: env_len
    INTEGER, INTENT(IN), OPTIONAL     :: read_write
    INTEGER, INTENT(IN), OPTIONAL     :: name_in_environ
    INTEGER, OPTIONAL                 :: ioPolicy
    INTEGER, OPTIONAL                 :: ioLocality
    LOGICAL, OPTIONAL                 :: allowRemap
    INTEGER, OPTIONAL                 :: error
    INTEGER, OPTIONAL                 :: fileType
    INTEGER                           :: errCode
    INTEGER                           :: lFileType
    CHARACTER(LEN=*)                  :: env
    REAL(KIND=jprb)     :: zhook_handle
    CHARACTER (LEN=*), PARAMETER                                               &
        :: SubroutineName=                                                     &
        'MODEL_FILE:MODEL_FILE_OPEN'
    CHARACTER(LEN=fileNameLength) :: filename

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)

    errCode=99

    IF (PRESENT(fileType)) THEN
      lFileType=fileType
    ELSE
      lFileType=ioFileTypeUM
    END IF

    IF (is_unit_open(unit)) THEN
      CALL mf_ereport(SubroutineName, errcode , 'Unit already open',unit)
    END IF

    IF (mf_data(unit) % managed ) THEN
      CALL mf_ereport(SubroutineName, errCode ,                                &
          'Unit already managed by model_file',unit)
    END IF

    CALL file_open (                                                           &
        unit,env,env_len,read_write,name_in_environ,error,                     &
        ioLocality,ioPolicy,allowRemap,lFileType)

    IF ((io_external_control).AND.(.NOT.(isReadOnly(unit))) ) THEN
      filename=TRIM(returnFilename(unit))
      CALL completion_open(unit,filename)
    END IF

    mf_data(unit) % managed = .TRUE.

    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END SUBROUTINE model_file_open


!----------------------------------------------------------------------
! Subroutine model_file_close
!----------------------------------------------------------------------
  SUBROUTINE model_file_close(unit,env,env_len,name_in_environ,delete,error)

    IMPLICIT NONE
    INTEGER, INTENT(IN)              :: unit
    CHARACTER(LEN=*), INTENT(IN)     :: env
    INTEGER, INTENT(IN), OPTIONAL    :: env_len  ! see io.F90 for details
    INTEGER, INTENT(IN), OPTIONAL    :: delete
    INTEGER, INTENT(IN), OPTIONAL    :: name_in_environ
    INTEGER, INTENT(OUT), OPTIONAL   :: error
    INTEGER                          :: errCode
    REAL(KIND=jprb)                  :: zhook_handle
    CHARACTER (LEN=*), PARAMETER     :: SubroutineName=                        &
        'MODEL_FILE:MODEL_FILE_CLOSE'

    LOGICAL :: test_complete_close
    CHARACTER(LEN=fileNameLength) :: filename

    IF (lhook) CALL dr_hook(SubroutineName,zhook_in,zhook_handle)

    IF (.NOT.mf_data(unit) % managed ) THEN
      errCode=99
      CALL mf_ereport(SubroutineName, errCode ,                                &
          'Unit not managed by model_file',unit)
    END IF

    CALL storeLookups(unit)

    mf_data(unit) % managed = .FALSE.

    IF ((io_external_control).AND.  &
        (.NOT.(isReadOnly(unit))) &
       ) THEN
       test_complete_close=.TRUE.
       filename=TRIM(returnfilename(unit))
    ELSE
       test_complete_close=.FALSE.
    END IF

    CALL file_close(unit,env,env_len,name_in_environ,delete,error)

    IF (test_complete_close) THEN
      CALL completion_close(unit,filename)
    END IF

    IF (lhook) CALL dr_hook(SubroutineName,zhook_out,zhook_handle)

  END SUBROUTINE model_file_close

!----------------------------------------------------------------------
! SUBROUTINE completion_open
! Suite Control : Remove any postfixed control file when file is opened
!----------------------------------------------------------------------

  SUBROUTINE completion_open(unit,filename)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    CHARACTER(LEN=fileNameLength), INTENT(IN) :: filename
    CHARACTER(LEN=fileNameLength+5) :: postfixed_name

    REAL(KIND=jprb)              :: zhook_handle

    IF (lhook) CALL dr_hook('COMPLETION_OPEN',zhook_in,zhook_handle)

    postfixed_name=TRIM(filename)//'.done'
    CALL file_delete(TRIM(postfixed_name),silent=.TRUE.)

    IF (lhook) CALL dr_hook('COMPLETION_OPEN',zhook_out,zhook_handle)

  END SUBROUTINE completion_open

!----------------------------------------------------------------------
! SUBROUTINE completion_close
! Suite Control : Generate a postfixed control file when a file is 
! closed
!----------------------------------------------------------------------

  SUBROUTINE completion_close(unit,filename)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    CHARACTER(LEN=fileNameLength), INTENT(IN) :: filename

    CHARACTER(LEN=fileNameLength+5) :: postfixed_name

    REAL(KIND=jprb)              :: zhook_handle

    IF (lhook) CALL dr_hook('COMPLETION_CLOSE',zhook_in,zhook_handle)

    postfixed_name=TRIM(filename)//'.done'

! File Op
    CALL file_touch(TRIM(postfixed_name))

    IF (lhook) CALL dr_hook('COMPLETION_CLOSE',zhook_out,zhook_handle)

  END SUBROUTINE completion_close

END MODULE model_file


