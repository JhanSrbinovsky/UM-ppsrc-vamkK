
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

!
! Handle asynchronous stash requests on an io server
!

MODULE IOS_Stash_Server
  USE io, ONLY :                                                               &
      setpos,                                                                  &
      buffout,                                                                 &
      buffo32
  USE mask_compression, ONLY :                                                 &
      compress_to_mask
  USE unite_output_files_mod, ONLY :                                           &
      unite_coex_files
  USE IOS_Queue_Mod
  USE IOS_Stash_Wgdos
  USE IOS_Common
  USE IOS_Stash_Common
  USE IOS_Model_Geometry
  USE UM_types
  USE mpl
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

  INTERFACE IOS_dump_pack
    MODULE PROCEDURE                                                           &
        IOS_dump_pack_1D,                                                      &
        IOS_dump_pack_2D
  END INTERFACE

  CHARACTER (LEN=*),                                                           &
      PARAMETER, PRIVATE :: RoutineName = 'IOS_Stash_Server'
  INTEGER, PARAMETER     :: stashLogUnit=10
  INTEGER, POINTER       :: atm_to_global(:)
  CHARACTER (LEN=80),                                                          &
      PRIVATE            :: iosStashServerMessage

! Debugging type parameters
  LOGICAL                :: disable_packing
  LOGICAL                :: disable_writes
  LOGICAL                :: disable_subdomaining

! Statistics for MPI comms
  INTEGER, PARAMETER     :: async_recv_calls  =1
  INTEGER, PARAMETER     :: async_recv_bytes  =2
  INTEGER, PARAMETER     :: async_recv_bytes2 =3
  INTEGER, POINTER       :: async_bytes(:,:)

! params/vars  for dr_hook
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CONTAINS

  LOGICAL FUNCTION isIOS(x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: x
    INTEGER             :: i
    INTEGER             :: j

    isIOS=.FALSE.
    DO i=1,IOS_Server_Groups
      DO j=1,model_procs
        IF (io_servers(i,j)==x)isIOS=.TRUE.
      END DO
    END DO
  END FUNCTION isIOS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialise the stash async subsystem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ios_stash_server_init(                                            &
      stash_active,                                                            &
      dumps_active,                                                            &
      disable_packing_in,                                                      &
      disable_writes_in,                                                       &
      disable_subdomaining_in)

    IMPLICIT NONE

    LOGICAL, INTENT(IN) :: disable_writes_in
    LOGICAL, INTENT(IN) :: stash_active
    LOGICAL, INTENT(IN) :: dumps_active
    LOGICAL, INTENT(IN) :: disable_packing_in
    LOGICAL, INTENT(IN) :: disable_subdomaining_in
    INTEGER             :: f
    INTEGER             :: l
    INTEGER             :: ierr
    INTEGER             :: myRank
    INTEGER             :: i
    INTEGER             :: k
    REAL(KIND=jprb)     :: zhook_handle
    CHARACTER (LEN=120) :: stashLogName

    IF (lhook) CALL dr_hook('IOS_STASH_SERVER:INIT',zhook_in,zhook_handle)

    disable_packing=disable_packing_in
    disable_writes=disable_writes_in
    disable_subdomaining=disable_subdomaining_in
    CALL useAsyncStash(stash_active)
    CALL useAsyncDump (dumps_active)

    IF (stash_active.OR.dumps_active) THEN

      atm_numprocs=global_procs-SIZE(io_servers)

      ALLOCATE(atm_to_global(0:atm_numprocs-1))!never deallocated

      IF (IOS_asyncDoStats) THEN
        ALLOCATE(async_bytes(3,0:atm_numprocs-1))!never deallocated
        async_bytes(:,:)=0
      END IF

      k=0
      DO i=0,global_procs-1
        IF (.NOT.isIOS(i)) THEN
          atm_to_global(k)=i
          k=k+1
        END IF
      END DO

      WRITE(6,'(A,I5,A,I3,A)')                                                 &
          'IOS: Info: Stash Server: Initialised: There are ',                  &
          atm_numprocs,                                                        &
          ' atm procs and ',                                                   &
          SIZE(io_servers),' io procs'
      CALL mpl_comm_rank(mpl_comm_world,myRank,ierr)

      IF (model_rank == 0) THEN
        WRITE(stashLogName,'(A,I4.4)')'ioserver_stash_log.',myRank
        WRITE(6,'(A,A)')'IOS: Info: Stash Server: Logging to ',                &
            TRIM(stashLogName)
        OPEN(unit=stashLogUnit,file=TRIM(stashLogName))

        WRITE(stashLogUnit,'(A,A)')                                            &
            '   time hnd unt    position  datasize ',                          &
            ' disksize  blk fld full   S    N    W    E pack'
      END IF
    END IF

    IF (lhook) CALL dr_hook('IOS_STASH_SERVER:INIT',zhook_out,zhook_handle)

  END SUBROUTINE ios_stash_server_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Shut down the stash async subsystem
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE ios_stash_server_fini()

    USE IOS_Server_Coupler, ONLY :                                             &
        proc_start,                                                            &
        proc_end

    IMPLICIT NONE

    INTEGER             :: proc
    REAL                :: ave
    REAL                :: var
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_STASH_SERVER:FINI',zhook_in,zhook_handle)

    IF (isUsingAsyncStash().OR.isUsingAsyncDumps()) THEN

      IF (model_rank == 0)                                                     &
          CLOSE(stashLogUnit)

      IF (IOS_asyncDoStats) THEN

        WRITE(6,'(A)')''
        WRITE(6,'(A)')'MPI counts for asynchronous operations'
        WRITE(6,'(A)')''
        WRITE(6,'(A4,A10,A12,A12)')'----','----------','------------',         &
            '------------'
        WRITE(6,'(A4,A10,A12,A12)')' CPU','  Receives','  Ave. Bytes',         &
            '        S.D.'
        WRITE(6,'(A4,A10,A12,A12)')'----','----------','------------',         &
            '------------'
        DO proc=proc_start,proc_end
          IF (async_bytes(async_recv_calls,proc)>0) THEN
            ave=1.0*async_bytes(async_recv_bytes,proc)/                        &
                async_bytes(async_recv_calls,proc)
            var=1.0*async_bytes(async_recv_bytes2,proc)/                       &
                async_bytes(async_recv_calls,proc)-ave*ave
          ELSE
            ave=0.0
            var=0.0
          END IF
          WRITE(6,'(I4,I10,F12.2,F12.2)')                                      &
              proc,                                                            &
              async_bytes(async_recv_calls,proc),                              &
              ave,                                                             &
              SQRT(var)
        END DO
        WRITE(6,'(A4,A10,A12,A12)')'----','----------','------------',         &
            '------------'
        DEALLOCATE(async_bytes)
        NULLIFY(async_bytes)
      END IF

      DEALLOCATE(atm_to_global)
      NULLIFY(atm_to_global)

    END IF
    WRITE(6,'(A)')'IOS: Info: Async service terminated.'

    IF (lhook) CALL dr_hook('IOS_STASH_SERVER:FINI',zhook_out,zhook_handle)
  END SUBROUTINE ios_stash_server_fini

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Post receive calls for inbound components of field
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE IOS_stash_server_init_recvs(node)!called by thr 0
    USE mpl, ONLY :                                                            &
        mpl_real,                                                              &
        mpl_request_null
    USE IOS_Server_Coupler, ONLY :                                             &
        proc_start,                                                            &
        proc_end

    IMPLICIT NONE

    TYPE(IOS_node_type),                                                       &
        INTENT(INOUT)  :: node
    INTEGER            :: receive_buffer_size
    INTEGER            :: proc
    INTEGER            :: FieldType
    INTEGER            :: j,i
    INTEGER            :: control_block_len
    INTEGER            :: num_records
    INTEGER            :: num_fields
    INTEGER            :: i_error
    REAL(KIND=jprb)    :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_STASH_SERVER:INIT_RECVS',zhook_in,zhook_handle)

    receive_buffer_size=getMaxFieldDomain()*levels_in_pack(node)

! we need to allow for the queue to accommodate the receive
! buffers even if we dont allocate them here
! otherwise we may see a deadlock at the actual allocation time

! DEPENDS ON : get_wallclock_time
    IF (IOS_Verbosity>=IOS_PrStatus_Diag                                       &
        .AND. model_rank ==0)                                                  &
        WRITE(stashLogUnit,'(F8.2,2I4,A)')                                     &
        get_wallclock_time()-IOS_Start_Time,                                   &
        node%metadata%handle,node%metadata%unit,                               &
        ' Timing: Receive transaction'


    NULLIFY (node%receive_tags)
    NULLIFY (node%receive_data_len)
    NULLIFY (node%receive_requests)
    NULLIFY (node%distributed_data)

    IF (.NOT.use_blocking_recvs) THEN
! If we are using non-blocking then we need to track the
! request objects
      ALLOCATE( node%receive_requests(proc_start:proc_end))
      !deallocate in fini_recvs
      ALLOCATE( node%distributed_data                                          &
          (receive_buffer_size,proc_start:proc_end) )
      !deallocate in stash_server_deallocate
    END IF

    ALLOCATE (node%receive_tags(proc_start:proc_end))
    ALLOCATE (node%receive_data_len(proc_start:proc_end))

! Pre calculation of tags/errors etc
    DO proc=proc_start,proc_end
      node%receive_data_len(proc) = payload_size(node,proc)
      node%receive_tags(proc)     = receive_tag(node,proc)

      IF (node%receive_data_len(proc)<0) THEN
        WRITE(iosStashServerMessage,'(A,I20,A)')                               &
            'recv ',                                                           &
            node%receive_data_len(proc),                                       &
            ' is negative '
        CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,              &
            md=node%metadata,unit=node%metadata%unit)
      END IF

      IF (node%receive_data_len(proc)>receive_buffer_size) THEN
        WRITE(iosStashServerMessage,'(A,I10,A,I6,A,I5,A,I5)')                  &
            'recv ',                                                           &
            node%receive_data_len(proc),                                       &
            ' from ',                                                          &
            proc,'/',atm_numprocs,                                             &
            ' exceeds local buffer size of ',receive_buffer_size
        CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,              &
            md=node%metadata,unit=node%metadata%unit)
      END IF

      IF (IOS_asyncDoStats) THEN
        IF (node%receive_data_len(proc)>0) THEN

          async_bytes(async_recv_calls,proc)=                                  &
              async_bytes(async_recv_calls,proc)+1
          async_bytes(async_recv_bytes,proc)=                                  &
              async_bytes(async_recv_bytes,proc)+                              &
              (node%receive_data_len(proc)*8)
          async_bytes(async_recv_bytes2,proc)=                                 &
              async_bytes(async_recv_bytes2,proc)+                             &
              (node%receive_data_len(proc)*8)*                                 &
              (node%receive_data_len(proc)*8)

        END IF
      END IF

    END DO

    IF (.NOT.use_blocking_recvs) THEN
      CALL acquire_lock_mpi()

      DO proc=proc_start,proc_end

        IF (node%receive_data_len(proc)>0) THEN
          CALL mpl_irecv( node%distributed_data(1,proc),                       &
              node%receive_data_len(proc),                                     &
              mpl_real,                                                        &
              atm_to_global(proc),                                             &
              node%receive_tags(proc),                                         &
              IOS_Async_Comm, node%receive_requests(proc),                     &
              i_error)
        ELSE
          node%receive_requests(proc)=mpl_request_null
        END IF

      END DO

      CALL release_lock_mpi()

      DEALLOCATE (node%receive_tags)
      DEALLOCATE (node%receive_data_len)
      NULLIFY(node%receive_data_len)
      NULLIFY(node%receive_tags)

    END IF

    IF (lhook)CALL dr_hook('IOS_STASH_SERVER:INIT_RECVS',zhook_out,zhook_handle)

  END SUBROUTINE IOS_stash_server_init_recvs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Wait for data to arrive for the given node
!
! This function to be called by thread 1. The function returns when
! receives associated with the async op have completed. This function
! does not actually perform the receives.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE IOS_stash_server_wait_for_data(node)

    IMPLICIT NONE
    TYPE(IOS_node_type),INTENT(INOUT) :: node
    LOGICAL                           :: flag
    REAL(KIND=jprb)                   :: zhook_handle


    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:WAIT_FOR_DATA',zhook_in,zhook_handle)
    END IF
! DEPENDS ON : get_wallclock_time
    IF (IOS_Verbosity>=IOS_PrStatus_Diag                                       &
        .AND. model_rank ==0)                                                  &
        WRITE(stashLogUnit,'(F8.2,2I4,A)')                                     &
        get_wallclock_time()-IOS_Start_Time,                                   &
        node%metadata%handle,node%metadata%unit,                               &
        ' Timing: Transaction at q head'

    IF (IOS_thread_0_calls_mpi) THEN
      ! We need to tell thread 0 that we need their help
      IOS_AsyncCompletionRequested=1
!$OMP FLUSH
      WRITE(6,'(A)')                                                           &
          'Waiting for thread 0 to complete a request for stash'
      DO WHILE  (IOS_AsyncCompletionRequested==1)
        CALL um_sleep(IOS_backoff_interval)
!$OMP FLUSH
      END DO
    ELSE
      ! This is called from thread 1, so flag can be ignored
      ! completion occurs on return
      flag=IOS_Stash_Server_Finish_Recvs(node)
    END IF
! DEPENDS ON : get_wallclock_time
    IF (IOS_Verbosity>=IOS_PrStatus_Diag .AND. model_rank ==0) THEN
      WRITE(stashLogUnit,'(F8.2,2I4,A)')                                       &
          get_wallclock_time()-IOS_Start_Time,                                 &
          node%metadata%handle,node%metadata%unit,                             &
          ' Timing: Data transmission complete'
    END IF
    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:WAIT_FOR_DATA',zhook_out,zhook_handle)
    END IF

  END SUBROUTINE IOS_stash_server_wait_for_data


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Complete receiving inbound components of field
!
! If we are called from thread 0 we may return .FALSE.
! if the operation is not yet complete and use_blocking_receives
! is not set. This is so that the listener thread can continue
! listening for other work whilst the comms are ongoing. If we are
! called from thread N!=0 we will block on completion of recvs and
! the function will return always return .TRUE. or deadlock waiting.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  LOGICAL FUNCTION IOS_stash_server_finish_recvs(node)

    USE IOS_Server_Coupler, ONLY :                                             &
        proc_start,                                                            &
        proc_end

    IMPLICIT NONE

    TYPE(IOS_node_type),INTENT(INOUT):: node

    INTEGER,POINTER    :: status(:,:)
    INTEGER            :: proc
    REAL               :: t1
    REAL               :: t2
    REAL               :: data_len
    INTEGER            :: thread
    LOGICAL            :: flag
    INTEGER            :: mpiStat(mpl_status_size)
    INTEGER            :: j
    INTEGER            :: mymax
    INTEGER            :: receive_buffer_size
    INTEGER            :: i_error
    REAL(KIND=jprb)    :: zhook_handle

    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:FINI_RECVS',zhook_in,zhook_handle)
    END IF

    thread=omp_get_thread_num()
    IOS_Stash_Server_Finish_Recvs=.FALSE.

    IF (use_blocking_recvs) THEN
      receive_buffer_size=getMaxFieldDomain()*levels_in_pack(node)

      ALLOCATE(node%distributed_data                                           &
          (receive_buffer_size,proc_start:proc_end))
      !deallocate in stash_server_deallocate

      IF (.NOT.ASSOCIATED(node%receive_tags)) THEN
        WRITE(iosStashServerMessage,'(A)')                                     &
            'MPI tag buffer not associated, exiting'
        CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,              &
            md=node%metadata,unit=node%metadata%unit)
      END IF
      IF (.NOT.ASSOCIATED(node%receive_data_len)) THEN
        WRITE(IosStashServerMessage,'(A)')                                     &
            'MPI length buffer not associated, exiting'
        CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,              &
            md=node%metadata,unit=node%metadata%unit)

      END IF
      CALL acquire_lock_mpi()

! DEPENDS ON : get_wallclock_time
      t1=get_wallclock_time()
      data_len=0
      DO proc=proc_start,proc_end
        IF (node%receive_data_len(proc)>0) THEN
          CALL mpl_recv( node%distributed_data(1,proc),                        &
              node%receive_data_len(proc),                                     &
              mpl_real, atm_to_global(proc),                                   &
              node%receive_tags(proc), IOS_async_comm,                         &
              mpiStat, i_error)

        END IF
        data_len=data_len+node%receive_data_len(proc)
      END DO
! DEPENDS ON : get_wallclock_time
      t2=get_wallclock_time()

      data_len=data_len*8.0/1024.0/1024.0
      IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
        WRITE(6,'(A,F6.2,A,I8.2,A)')                                           &
            'IOS: Info: Recv: ',data_len,'MB @',data_len/(t2-t1),'MB/s'
      END IF
      CALL release_lock_mpi()

      DEALLOCATE ( node%receive_data_len)!allocated in init_recvs
      DEALLOCATE ( node%receive_tags )!allocated in init_recvs
      NULLIFY(node%receive_data_len)
      NULLIFY(node%receive_tags)
      IOS_Stash_Server_Finish_Recvs=.TRUE.
    ELSE
      ALLOCATE (status(mpl_status_size,proc_start:proc_end))
                                       ! deallocate at end of block

      IF (.NOT.ASSOCIATED(node%receive_requests)) THEN
        WRITE(iosStashServerMessage,'(A)')                                     &
            'MPI request buffer not associated, exiting'
        CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,              &
            md=node%metadata,unit=node%metadata%unit)
      END IF
      IF (.NOT.ASSOCIATED(node%distributed_data)) THEN
        WRITE(IosStashServerMessage,'(A)')                                     &
            'MPI receive buffer not associated, exiting'
        CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,              &
            md=node%metadata,unit=node%metadata%unit)
      END IF
      IF (thread > 0) THEN
        CALL acquire_lock_mpi()
        CALL MPL_WaitAll(proc_end-proc_start+1,node%receive_requests,          &
            status,i_error)
        CALL release_lock_mpi()
        IOS_Stash_Server_Finish_Recvs=.TRUE.
      ELSE! a CALL from thread 0 shouldnt block, nor do we need to lock
        CALL MPL_Testall(proc_end-proc_start+1,node%receive_requests,          &
            flag,status,i_error)
        IOS_Stash_Server_Finish_Recvs=flag
      END IF
      IF (IOS_Stash_Server_Finish_Recvs) THEN
        ! If the operation completed deallocate the
        ! receive request objects
        DEALLOCATE(node%receive_requests)!allocated in init_recvs
        NULLIFY (node%receive_requests)
      END IF

      DEALLOCATE(status)! allocated at start of block
      NULLIFY (status)
    END IF

    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:FINI_RECVS',zhook_out,zhook_handle)
    END IF

  END FUNCTION IOS_stash_server_finish_recvs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Process a stash object from the IO queue
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE IOS_stash_server_process(node)

    USE model_file, ONLY :                                                     &
        MF_Data_Address,                                                       &
        Get_File_Address,                                                      &
        attachLookups
    USE IOS_Server_Coupler, ONLY :                                             &
        proc_start,                                                            &
        proc_end,                                                              &
        grid_row_start,                                                        &
        grid_row_end,                                                          &
        grid_point_start,                                                      &
        grid_point_end
    USE IOS_Constants
    USE IOS_geometry_utils

    USE lookup_addresses
    IMPLICIT NONE

    TYPE(IOS_node_type)  :: node
    INTEGER              :: FieldType
    INTEGER,                                                                   &
        ALLOCATABLE      :: dataOffset(:)

    INTEGER, POINTER     :: lookupTable(:,:)
    INTEGER              :: lenOut32
    INTEGER              :: lenOut64
    INTEGER              :: ControlBlockLength
    INTEGER              :: x_range
    INTEGER              :: y_range
    INTEGER              :: recordNumber
    INTEGER              :: diskLocation
    INTEGER              :: totalDiskSize
    INTEGER              :: totalDataSize
    INTEGER              :: currentTotalDataSize
    INTEGER              :: numLandPoints
    INTEGER, ALLOCATABLE :: ControlRecord(:)
    INTEGER, ALLOCATABLE :: AutoControlRecord(:)
    INTEGER, ALLOCATABLE :: partialLengths(:)
    INTEGER, ALLOCATABLE :: displacements(:)
    REAL, ALLOCATABLE    :: finalData(:)
    REAL, ALLOCATABLE    :: gatheredData(:)
    INTEGER              :: repeatCount
    LOGICAL              :: doWrite
    LOGICAL              :: firstRecordInPack
    REAL,                                                                      &
        POINTER          :: field(:,:) !Stores the global field
    REAL(KIND=real32),                                                         &
        ALLOCATABLE      :: packedPartialData32(:)
    REAL, ALLOCATABLE    :: landCompressedField(:)
    TYPE(box)            :: domainBox

! Parameters from the input record
    INTEGER              :: preprocess
    INTEGER              :: packing
    INTEGER              :: subdomainType
    INTEGER              :: packingType
    INTEGER              :: compressionAccuracy
    INTEGER              :: landmaskCompression
    INTEGER              :: fullField
    INTEGER              :: n,s,e,w !inclusive bounds for subfields
    INTEGER              :: DiskBlockSize  !write data in lumps of this
    INTEGER              :: DiskBlockStart !write address
    REAL                 :: missingDataIndicator

! identification of stash fields
    INTEGER              :: section
    INTEGER              :: code
    INTEGER              :: level
    INTEGER              :: dummy

! Various local counters/loop vars etc
    INTEGER              :: j
    INTEGER              :: iy
    INTEGER              :: xlow
    INTEGER              :: xhi
    INTEGER              :: ylow
    INTEGER              :: yhi
    INTEGER              :: proc

! Error codes, placeholders, and fluff
    INTEGER              :: len_io
    INTEGER              :: ErrorStatus
    INTEGER,                                                                   &
        PARAMETER        :: rootPE=0
    REAL                 :: iostat
    REAL(KIND=jprb)      :: zhook_handle
    CHARACTER (LEN=*),                                                         &
        PARAMETER        :: RoutineName = 'IOS_STASH_SERVER:PROCESS'


    IF (lhook) CALL dr_hook(RoutineName, zhook_in, zhook_handle)

    NULLIFY(lookupTable)
    NULLIFY(field)

    ALLOCATE(dataOffset(0:atm_numprocs-1))
    ALLOCATE(AutoControlRecord(1:ios_stash_control_auto_len))
    ALLOCATE(ControlRecord(1:ios_async_control_sz_max))
    ALLOCATE(partialLengths(model_procs))
    ALLOCATE(displacements(model_procs))

    dataOffset(0:atm_numprocs-1)=1
    firstRecordInPack=.TRUE.
    repeatCount=0


    ! Just keep moving through the control block, we don't know how long it is
    j=1
    DO WHILE (j<node%metadata%data_size)
      AutoControlRecord(1:ios_stash_control_auto_len) =                        &
          node%integer_data(j:j+ios_stash_control_auto_len-1)
      ! IF this word is a start marker....
      IF (AutoControlRecord(loc_record_type)                                   &
          ==ios_stash_record_start) THEN
        ControlBlockLength=                                                    &
            -1*AutoControlRecord(loc_record_len_control)
        j=j+ios_stash_control_auto_len!1st element provided by user
        IF (node%integer_data(j)==IOS_Repeat_Record) THEN
          repeatCount=repeatCount+1
        ELSE IF (node%integer_data(j)==                                        &
            ios_stash_distributed_field) THEN

          repeatCount=0
          ControlRecord(1:ios_async_control_sz_max) =                          &
              node%integer_data(j:j+ios_async_control_sz_max-1)
        ELSE
          WRITE(6,'(A,I4,I4)')'Bad Marker: ',                                  &
              node%integer_data(j),ControlRecord(loc_record_type)
        END IF
        IF (ControlRecord(loc_record_type)==                                   &
            ios_stash_distributed_field) THEN
          ! This looks like a field and we so we process it
          ! Lets get the contents into more meaningful variables
          FieldType           = ControlRecord(loc_fld_type)
          preprocess          = ControlRecord(loc_preprocess_flag)
          fullField           = ControlRecord(loc_subdomain_flag)
          packing             = ControlRecord(loc_packing_flag)
          packingType         = ControlRecord(loc_pack_type)
          compressionAccuracy = ControlRecord(loc_comp_accry)
          DiskBlockSize       = ControlRecord(loc_disk_block)
          DiskBlockStart      = ControlRecord(loc_disk_block_start)
          recordNumber        = ControlRecord(loc_seek_target)
          landmaskCompression = ControlRecord(loc_landmaskcompress)
          CALL IOS_unpack4(controlRecord(loc_boundary),n,s,e,w)
! recover that real value from the integer array...
          CALL um_memcpy_f                                                     &
              (MissingDataIndicator,ControlRecord(loc_dmi),1)

          ! Do not use the control record below here as we will modify values

          ! Because we may have a repeat record, update recordNumber
          recordNumber=recordNumber+repeatCount

          IF (FieldType < 1 .OR. FieldType > nfld_types ) THEN
            WRITE(iosStashServerMessage,'(A,I3)')                              &
                'ERROR INVALID FLD TYPE ',FieldType
            CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,          &
                md=node%metadata,unit=node%metadata%unit)
          END IF

          ALLOCATE(field(atm_global_points(1,FieldType),                       &
              atm_global_points(2,FieldType)))

          ! If the disable subdomaining flag is set, then
          ! zero the field so that parts we don't have data for
          ! (ie aren't in the requested domain) aren't garbage.

          IF (disable_subdomaining)                                            &
              field(:,:)=missingDataIndicator


          ! ****************************************************
          ! Extract into field data from everyone who sent us
          ! some component of the field
          ! ****************************************************

          !Assemble the field from the multitude of receive buffers
          domainBox%n = n
          domainBox%s = s
          domainBox%e = e
          domainBox%w = w

          DO proc=proc_start,proc_end

            IF (fullField == IOS_partial_field) THEN
              subdomainType=classify_subdomain(proc,FieldType,domainBox)
            ELSE
              subdomainType=complete_intersection
            END IF

            IF (subdomainType /= no_intersection                               &
                .OR. IOS_AsyncSendNull) THEN

              xlow=offset_map(1,FieldType,proc)
              xhi=offset_map(1,FieldType,proc)+                                &
                  size_map(1,FieldType,proc)-1
              ylow=offset_map(2,FieldType,proc)
              yhi=offset_map(2,FieldType,proc)+                                &
                  size_map(2,FieldType,proc)-1


              DO iy=ylow,yhi
                field(xlow:xhi,iy)=                                            &
                    node%distributed_data(dataOffset(proc):                    &
                    dataOffset(proc)+xhi-xlow,proc)
                dataOffset(proc)=dataOffset(proc)+xhi-xlow+1
              END DO
            END IF
          END DO

          IF (preprocess==ios_stash_preprocess) THEN
            ! Some kind of data processing required

            IF (node%metadata%action==IOS_Action_StashWritePPData              &
                .AND. disable_packing) THEN
              IF (packing==ios_packing) THEN
                WRITE(6,'(A)')                                                 &
                    'IOS: Warning: DEBUG: Disabling packing option active!'
                packing=ios_no_packing
                packingType=0
              END IF
            END IF

            doWrite=.TRUE.
            IF (disable_writes) THEN
              WRITE(6,'(A,I3)')                                                &
                  'IOS: Warning: DEBUG: Disabling data writes for unit ',      &
                  node%metadata%unit
              doWrite=.FALSE.
            END IF

            IF (node%metadata%action==IOS_Action_StashWritePPData              &
                .AND. disable_subdomaining) THEN
              IF (fullField==ios_partial_field) THEN
                WRITE(6,'(A)')                                                 &
                    'IOS: Warning: DEBUG: Disabling subdomaining option active!'
                fullField=ios_full_field
              END IF
            END IF

            ! ****************************************************
            ! From our assembled field construct a correctly sized
            ! array and copy in the data
            ! subdomain = INTERSECTION(req. domain, this IOS task)
            ! ****************************************************

            ! Ensure N,S,E,W are sane for full fields
            IF (fullField==ios_full_field) THEN
              S=1
              N=atm_global_points(2,FieldType)
              W=1
              E=atm_global_points(1,FieldType)
            END IF

            CALL IOS_stash_server_subdomain(                                   &
                field,        & ! on input a pointer to the global field array
                fullField,    & ! whether or not this is a full field request
                FieldType,    & !
                MIN(n,grid_row_end  (fieldtype)),                              &
                MAX(s,grid_row_start(fieldtype)),                              &
                E,                                                             &
                W,                                                             &
                x_range,                                                       &
                y_range)


            ! If there was a subdomain 'field' now points to the
            ! smaller region.

            IF (node%metadata%action==IOS_Action_StashWritePPData) THEN

              ! ****************************************************
              ! Operations for stash data
              ! ****************************************************


              IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
                CALL IOS_unpack4(controlRecord(loc_id),                        &
                    section,                                                   &
                    code,                                                      &
                    level,                                                     &
                    dummy)
                WRITE(6,'(A,I3,A,I6)')                                         &
                    'IOS: Info: Stash Server: unit=',                          &
                    node % metadata % unit,                                    &
                    ' field number=',recordNumber
                WRITE(6,'(A,I4,A,I4,A,I4)')                                    &
                    'IOS: Info: Stash Server: Section= ',                      &
                    section,' Code= ',code,                                    &
                    ' Level (relative to base)= ',level
              END IF

              ! Pack as needed

              ALLOCATE (packedPartialData32                                    &
                  (2*(x_range+coexHeaderWords32)*y_range+coexHeaderWords32+1))
                
              CALL IOS_stash_pack(                                             &
                  field,                                                       &
                  packedPartialData32,                                         &
                  packing,                                                     &
                  packingType,                                                 &
                  compressionAccuracy,                                         &
                  MissingDataIndicator,                                        &
                  lenOut32,                                                    &
                  lenOut64,                                                    &
                  FieldType,                                                   &
                  x_range,                                                     &
                  y_range,                                                     &
                  DiskBlockSize)

              CALL acquire_lock_mpi()
              CALL MPL_Gather(                                                 &
                  lenOut64,                                                    &
                  1,                                                           &
                  mpl_integer,                                                 &
                  partialLengths,                                              &
                  1,                                                           &
                  mpl_integer,                                                 &
                  rootPE,                                                      &
                  model_comm,                                                  &
                  errorStatus)
              CALL release_lock_mpi()

              IF (model_rank==0) THEN
                displacements(1)=0
                DO proc=2,model_procs
                  displacements(proc)=displacements(proc-1)+                   &
                      partialLengths(proc-1)
                END DO

                totalDataSize=SUM(partialLengths(:))

                ALLOCATE                                                       &
                    (finalData((totalDataSize)+DiskBlockSize))
                ALLOCATE                                                       &
                    (gatheredData((totalDataSize)+DiskBlockSize))
              ELSE
                ALLOCATE(finalData(1))
                ALLOCATE(gatheredData(1))
              END IF

              CALL acquire_lock_mpi()
              CALL MPL_GatherV(                                                &
                  packedPartialData32,                                         &
                  lenOut64,                                                    &
                  mpl_integer8,                                                &
                  gatheredData,                                                &
                  partialLengths,                                              &
                  displacements,                                               &
                  mpl_integer8,                                                &
                  rootPE,                                                      &
                  model_comm,                                                  &
                  errorStatus)
              CALL release_lock_mpi()

              IF (model_rank==0) THEN

                ! Move the file pointer to the right place
                diskLocation=                                                  &
                    IOS_stash_server_seek                                      &
                    (node%metadata%unit, recordNumber)

                IF ( packing == IOS_packing ) THEN

                  DO proc=1,model_procs
                    CALL unite_coex_files(                                     &
                        gatheredData(displacements(proc)+1),                   &
                        finaldata,                                             &
                        totalDataSize,                                         &
                        proc-1)
                  END DO

                  ! convert totalDataSize to whole 64 bit words
                  totalDataSize=(totalDataSize+1)/2

                  ! Write it onto disk
                  CALL IOS_stash_server_write(                                 &
                      node,                                                    &
                      finalData,                                               &
                      totalDataSize,                                           &
                      DiskBlockSize,                                           &
                      totalDiskSize,                                           &
                      doWrite)
                ELSE
                  CALL IOS_stash_server_write(                                 &
                      node,                                                    &
                      gatheredData,                                            &
                      totalDataSize,                                           &
                      DiskBlockSize,                                           &
                      totalDiskSize,                                           &
                      doWrite)

                END IF
                ! ****************************************************
                ! For stash, update the lookup tables
                ! ****************************************************

                ! now: output_size    = data length (words)
                !      totalDiskSize  = amount written to disk (words)
                !      diskLocation   = location in file (words)
                lookupTable=>attachLookups(node % metadata % unit)

                lookupTable(lbegin , recordNumber)=diskLocation
                lookupTable(lblrec , recordNumber)=totalDataSize
                lookupTable(lbnrec , recordNumber)=totalDiskSize

                ! Note for fields files NADDR is the same as LBEGIN
                lookupTable(naddr  , recordNumber)=diskLocation
                lookupTable(lbpack , recordNumber)=0
                IF (packing==ios_packing) THEN
                  lookupTable(lbpack , recordNumber)=1
                END IF

                ! If we disabled subdomains then we should also
                ! correct entries in the lookup which will be wrong
                ! compared to what the atmos model will later send.
                IF (disable_subdomaining) THEN
                  lookupTable(lbrow , recordNumber)=                           &
                      atm_global_points(2,FieldType)
                  lookupTable(lbnpt , recordNumber)=                           &
                      atm_global_points(1,FieldType)
                  lookupTable(bzx , recordNumber)= 0
                  lookupTable(bzy , recordNumber)= 0
                  lookupTable(lbhem , recordNumber)= 0
                  WRITE(6,'(A)')                                               &
                      'IOS: Warning: Exact details of grid are unknown'
                  WRITE(6,'(A,A)')                                             &
                      'IOS:          Field will be written',                   &
                      ' as global 0,0 origined'

                  ! We should set LBHEM, but we don't have enough data to do so.
                END IF

                NULLIFY(lookupTable)

                ! Make a note of what happened in the log.
                CALL IOS_stash_server_log(node,FieldType,                      &
                    preprocess, fullField,packing,s,n,w,e,                     &
                    MissingDataIndicator,DiskBlockSize,                        &
                    diskLocation,totalDiskSize*8,(totalDataSize+1)/2)

              END IF

              DEALLOCATE(packedPartialData32)
              DEALLOCATE(finalData)
              DEALLOCATE(gatheredData)

            ELSE IF (node%metadata%action ==                                   &
                IOS_Action_StashWriteDumpData) THEN

              ! ****************************************************
              ! Operations for dump data
              ! ****************************************************

              IF (landmaskCompression==ios_packing_type_landmask) THEN
                ! DiskBlockSize is the size of the record on disk
                ! not the sector size as it is for stash
                ALLOCATE(landCompressedField(x_range*y_range))
                CALL compress_to_mask(                                         &
                    field,                                                     &
                    landCompressedField,                                       &
                    land_mask(                                                 &
                    grid_point_start(fieldType):                               &
                    grid_point_end(fieldType)),                                &
                    x_range*y_range,                                           &
                    numLandPoints)
                ALLOCATE (packedPartialData32(2*numLandPoints))
                CALL IOS_dump_pack(                                            &
                    landCompressedField,                                       &
                    packedPartialData32,      &! field in/out
                    packing,                                                   &
                    packingType,              &! packing args
                    lenOut32,                                                  &
                    lenOut64,                 &! data lengths out
                    FieldType,                                                 &
                    numLandPoints,                                             &
                    DiskBlockSize)
                DEALLOCATE(landCompressedField)
              ELSE
                ALLOCATE (packedPartialData32(2*x_range*y_range))
                CALL IOS_dump_pack(                                            &
                    field,                                                     &
                    packedPartialData32,      & ! field in/out
                    packing,                  & ! packing args
                    packingType,                                               &
                    lenOut32,                 & ! data lengths out
                    lenOut64,                                                  &
                    FieldType,                                                 &
                    x_range,                                                   &
                    y_range,                                                   &
                    DiskBlockSize)

              END IF

              IF (model_rank==0) THEN
                IF (firstRecordInPack) THEN
                  ! we need to setpos
                  IF (IOS_Verbosity>=IOS_PrStatus_Diag)THEN
                    WRITE(6,'(A,A,I10)')                                       &
                        'IOS: Info: Stash Server: ',                           &
                        'Seek position for dump output=',                      &
                        DiskBlockStart
                  END IF
                  CALL setpos(node%metadata%unit,DiskBlockStart,               &
                      errorStatus)
                  firstRecordInPack=.FALSE.
                END IF
              END IF

              CALL acquire_lock_mpi()
              CALL MPL_Gather(                                                 &
                  lenOut32,                                                    &
                  1,                                                           &
                  mpl_integer,                                                 &
                  partialLengths,                                              &
                  1,                                                           &
                  mpl_integer,                                                 &
                  rootPE,                                                      &
                  model_comm,                                                  &
                  errorStatus)
              CALL release_lock_mpi()

              IF (model_rank==0) THEN
                displacements(1)=0
                DO proc=2,model_procs
                  displacements(proc)=                                         &
                      displacements(proc-1)+partialLengths(proc-1)
                END DO
                totalDataSize=SUM(partialLengths(:))
                ALLOCATE(finalData(DiskBlockSize))
                finaldata(:)=100.0
              ELSE
                ALLOCATE(finalData(1))
              END IF

              CALL acquire_lock_mpi()
              CALL MPL_GatherV(                                                &
                  packedPartialData32,                                         &
                  lenOut32,                                                    &
                  mpl_integer4,                                                &
                  finalData,                                                   &
                  partialLengths,                                              &
                  displacements,                                               &
                  mpl_integer4,                                                &
                  rootPE,                                                      &
                  model_comm,                                                  &
                  errorStatus)
              CALL release_lock_mpi()

              IF (packingType==ios_packing_type_pack21) THEN

                IF (model_rank==0) THEN
!DEPENDS ON : buffout32_f77
                  CALL buffout32_f77(node%metadata%unit,                       &
                      finalData,                                               &
                      DiskBlockSize*2,                                         &
                      len_io,iostat)
                END IF

              ELSE
                IF (model_rank==0) THEN
                  CALL buffout(node%metadata%unit,                             &
                      finalData,                                               &
                      DiskBlockSize,                                           &
                      len_io,iostat)
                END IF

              END IF

              DEALLOCATE(packedPartialData32)
              DEALLOCATE(finalData)

            ELSE ! Wrong action somehow
              CALL IOS_Ereport(RoutineName,99,                                 &
                  'Wrong action in Stash_Server',                              &
                  md=node%metadata,unit=node%metadata%unit)

            END IF

          ELSE

            WRITE(iosStashServerMessage,'(A,I10)')                             &
                'unknown flag passed for preprocessing: ',preprocess
            CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,          &
                md=node%metadata,unit=node%metadata%unit)
          END IF
          DEALLOCATE(field)
          NULLIFY (field)
        ELSE
          WRITE(IosStashServerMessage,'(A,I10)')                               &
              'Unknown data object in control record at ',j
          CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,            &
              md=node%metadata,unit=node%metadata%unit)
        END IF ! is a distributed field
        ! Advance to the next record
        j=j+ControlBlockLength
      ELSE
        !We didnt get a start of record marker where expected :-(
        WRITE(iosStashServerMessage,'(A,I10)')                                 &
            'process: unexpected lack of record ',j
        CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,              &
            md=node%metadata,unit=node%metadata%unit)
      END IF
      IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
        WRITE(6,'(A,I3,A)')'IOS: Info: Stash Server: unit=',                   &
            node % metadata % unit,                                            &
            ' done record'
      END IF
    END DO

    IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
      WRITE(6,'(A,I3,A)')'IOS: Info: Stash Server: unit=',                     &
          node % metadata % unit,                                              &
          ' done pack'
    END IF

    DEALLOCATE(dataOffset)
    DEALLOCATE(AutoControlRecord)
    DEALLOCATE(ControlRecord)
    DEALLOCATE(partialLengths)
    DEALLOCATE(displacements)

    IF (lhook) CALL dr_hook(RoutineName, zhook_out, zhook_handle)

  END SUBROUTINE IOS_stash_server_process

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set the file position to the location for a record
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  FUNCTION IOS_stash_server_seek(unit, recordNumber)                           &
      RESULT(seekAddress)

    USE model_file, ONLY :                                                     &
        get_file_address,                                                      &
        MF_data_address,                                                       &
        attachLookups
    USE lookup_addresses

    IMPLICIT NONE

    INTEGER, INTENT(IN)         :: unit
    INTEGER, INTENT(IN)         :: recordNumber
    INTEGER, POINTER            :: fixed_header(:)
    INTEGER, POINTER            :: ipplook(:,:)
    INTEGER                     :: seekAddress
    INTEGER                     :: dataAddress
    INTEGER                     :: ierror
    REAL(KIND=jprb)             :: zhook_handle

    NULLIFY(fixed_header)
    NULLIFY(ipplook)

    IF (lhook) CALL dr_hook('IOS_STASH_SERVER:SEEK',zhook_in,zhook_handle)

    ipplook=>attachLookups(unit)
    dataAddress=get_file_address(unit,MF_data_address)
    IF (recordNumber < 1) THEN
      CALL IOS_Ereport('Stash_seek',99,'invalid record number',                &
          unit=unit)

    ELSE IF (recordNumber==1) THEN
      seekAddress = dataAddress
    ELSE
      ! previous record location + prev record size
      seekAddress =                                                            &
          ipplook(lbegin, recordNumber-1)+                                     &
          ipplook(lbnrec, recordNumber-1)
    END IF

    IF (seekAddress < MAX(0,dataAddress)) THEN
      WRITE(iosStashServerMessage,'(A,I16,A,I16,A)')                           &
          'Seek value (',seekAddress,                                          &
          ') lower than data start address (',dataAddress,                     &
          ') from fixed header'
      CALL IOS_Ereport('stash_server_seek',99,                                 &
          iosStashServerMessage,unit=unit)
    END IF

    IF (IOS_Verbosity>=IOS_PrStatus_Diag)THEN
      WRITE(6,'(A,I8,A,I10)')                                                  &
          'IOS: Info: Stash Server: decoded record ',recordNumber,             &
          ' to disk address ',seekAddress
    END IF
    CALL setpos(unit,seekAddress, ierror)
    NULLIFY(ipplook)

    IF (lhook) CALL dr_hook('IOS_STASH_SERVER:SEEK',zhook_out,zhook_handle)

  END FUNCTION IOS_stash_server_seek

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Perform compression on a stash object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE IOS_stash_pack(                                                   &
      input,                                                                   &
      output,                                                                  &
      packing,                                                                 &
      packingType,                                                             &
      compressionAccuracy,                                                     &
      MissingDataIndicator,                                                    &
      lenOut32,                                                                &
      lenOut64,                                                                &
      FieldType,                                                               &
      x_range,                                                                 &
      y_range,                                                                 &
      DiskBlockSize)

    USE IOS_stash_common
    USE IOS_stash_wgdos

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: DiskBlockSize
    INTEGER, INTENT(IN)          :: compressionAccuracy
    INTEGER, INTENT(IN)          :: packingType
    INTEGER, INTENT(IN)          :: FieldType
    INTEGER, INTENT(IN)          :: x_range,y_range
    INTEGER, INTENT(INOUT)       :: packing
    INTEGER, INTENT(OUT)         :: lenOut32
    INTEGER, INTENT(OUT)         :: lenOut64
    REAL, INTENT(IN)             :: MissingDataIndicator
    REAL, POINTER                :: input(:,:)
    REAL(KIND=real32),                                                         &
        INTENT(OUT)              :: output(:)
    REAL(KIND=jprb)              :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_STASH_SERVER:STASHPACK',zhook_in,zhook_handle)

    IF (packing==ios_packing) THEN
      IF (packingType==ios_packing_type_wgdos) THEN
        IF (compressionAccuracy>-99) THEN
          CALL ios_stash_pack_wgdos(input,output,                              &
              lenOut32,compressionAccuracy,                                    &
              MissingDataIndicator)
          IF (lenout32==0) THEN
            WRITE(6,'(A,A)')'IOS: Info: Stsh: Packing Failed...',              &
                ' writing field unpacked.'
          END IF
        ELSE
          IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
            WRITE(6,'(A,I4)')                                                  &
                'IOS: Info: Stash Server: No packing: accuracy=',              &
                compressionAccuracy
          END IF
          packing=ios_no_packing
        END IF
      ELSE IF (packingType==ios_packing_type_pack21) THEN
        WRITE(6,'(A,A)')'IOS: Info: Stsh: PACK2D: pac21 packing...',           &
            'not implemented'
        packing=ios_no_packing
      ELSE
        WRITE(6,'(A)')'IOS: Info: Stsh: Unknown packing code'
        packing=ios_no_packing
      END IF
      lenOut64=(lenOut32+1)/2
    END IF

! We may enter this because WGDOS packing failed.
    IF (packing==ios_no_packing) THEN
      ! Literal copy of real data into integer output storage.
      CALL um_memcpy64(output,input,x_range*y_range)
      lenOut64=x_range*y_range
      lenOut32=lenOut64*2
    END IF

    IF (packing/=ios_packing.AND.packing/=ios_no_packing) THEN
      WRITE(iosStashServerMessage,'(A,I20)')                                   &
          'Packing flag incorrectly set ',packing
      CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage )
    END IF

    IF (lhook)CALL dr_hook('IOS_STASH_SERVER:STASHPACK',zhook_out,zhook_handle)

  END SUBROUTINE IOS_stash_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Perform compression on a dump object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE IOS_dump_pack_2D(                                                 &
      input,                                                                   &
      output32,                                                                &
      packing,                                                                 &
      packingType,                                                             &
      lenOut32,                                                                &
      lenOut64,                                                                &
      FieldType,                                                               &
      x_range,                                                                 &
      y_range,                                                                 &
      DiskBlockSize)

    USE IOS_stash_common
    USE IOS_stash_wgdos

    IMPLICIT NONE
    REAL, INTENT(IN)             :: input(:,:)
    REAL(KIND=real32)            :: output32(:)
    INTEGER, INTENT(OUT)         :: lenOut32
    INTEGER, INTENT(OUT)         :: lenOut64
    INTEGER, INTENT(INOUT)       :: packing
    INTEGER, INTENT(IN)          :: packingType
    INTEGER, INTENT(IN)          :: FieldType
    INTEGER, INTENT(IN)          :: x_range
    INTEGER, INTENT(IN)          :: y_range
    INTEGER                      :: DiskBlockSize
    REAL(KIND=jprb)              :: zhook_handle

    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:DUMPPACK2D',zhook_in,zhook_handle)
    END IF

    IF (packing==ios_packing) THEN
      IF  (packingType==ios_packing_type_pack21) THEN
        CALL pack21(x_range*y_range,input,output32)
      ELSE
        WRITE(6,'(A)')'unknown packing code'
        packing=ios_no_packing
      END IF
      lenOut32=x_range*y_range
      lenOut64=(lenOut32+1)/2
    END IF

    ! We may enter this because of an unknown packing code
    IF (packing==ios_no_packing)THEN
      ! Literal copy of real data into integer output storage.
      CALL um_memcpy64(output32,input,x_range*y_range)
      lenOut64=x_range*y_range
      lenOut32=lenOut64*2
    END IF

    IF (packing/=ios_packing .AND. packing/=ios_no_packing) THEN
      WRITE(iosStashServerMessage,'(A,I20)')                                   &
          'Packing flag incorrectly set ',packing
      CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage )
    END IF

    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:DUMPPACK2D',zhook_out,zhook_handle)
    END IF

  END SUBROUTINE IOS_dump_pack_2D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Perform compression on an object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE IOS_dump_pack_1D(input,output32,                                  &
      packing,packingType,                                                     &
      lenOut32,lenOut64,                                                       &
      FieldType,num_points,DiskBlockSize)
    USE IOS_stash_common
    USE ios_stash_wgdos

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: packingType
    INTEGER, INTENT(IN)          :: FieldType
    INTEGER, INTENT(IN)          :: num_points
    INTEGER, INTENT(IN)          :: DiskBlockSize
    INTEGER, INTENT(OUT)         :: lenOut32
    INTEGER, INTENT(OUT)         :: lenOut64
    INTEGER, INTENT(INOUT)       :: packing
    REAL, INTENT(IN)             :: input(:)
    REAL(KIND=real32)            :: output32(:)
    REAL(KIND=jprb)              :: zhook_handle

    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:DUMPPACK1D',zhook_in,zhook_handle)
    END IF

    IF (packing==ios_packing) THEN
      IF (packingType==ios_packing_type_pack21) THEN
        CALL pack21(num_points,input,output32)
      ELSE
        WRITE(6,'(A)')'Unknown packing code (pack 1D)'
        packing=ios_no_packing
      END IF
      lenOut32 = num_points
      lenOut64=(lenOut32+1)/2
    END IF

! We may enter this because WGDOS packing failed.
    IF (packing==ios_no_packing) THEN
      ! Literal copy of real data into integer output storage.

      CALL um_memcpy64(output32,input,num_points)
      lenOut64=num_points
      lenOut32=lenOut64*2

    END IF

    IF (packing/=ios_packing .AND. packing/=ios_no_packing) THEN
      WRITE(iosStashServerMessage,'(A,I20)')                                   &
          'Packing flag incorrectly set ',packing
      CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage )
    END IF

    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:DUMPPACK1D',zhook_out,zhook_handle)
    END IF

  END SUBROUTINE IOS_dump_pack_1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Write a stash record to disk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE IOS_stash_server_write(node,outData,data_length,                  &
      block_size,written,do_write)

    IMPLICIT NONE

    TYPE(IOS_node_type),                                                       &
        INTENT(INOUT)        :: node
    REAL, INTENT(IN)         :: outData(:)
    INTEGER, INTENT(IN)      :: data_length
    INTEGER, INTENT(IN)      :: block_size
    LOGICAL, INTENT(IN)      :: do_write
    INTEGER, INTENT(OUT)     :: written
    INTEGER                  :: io_len
    REAL                     :: rstat
    REAL(KIND=jprb)          :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_STASH_SERVER:WRITE',zhook_in,zhook_handle)

    io_len =((data_length+block_size-1)/block_size)*block_size

    IF (do_write) THEN
      CALL buffout( node%metadata%unit, outData,                               &
          io_len, written, rstat )
      IF (written /= io_len .OR. rstat /= -1.0) THEN
        IF (written /= io_len) THEN
          WRITE(iosStashServerMessage,'(A,A,I20,I20)')                         &
              'IOS:server:stash: ',                                            &
              'Mismatch between bytes written/requested',                      &
              written,io_len
          CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage )
        END IF
        IF (rstat /= -1.0) THEN
          WRITE(iosStashServerMessage,'(A,I4,I4)')                             &
              'IOS:server:stash: Bad Status:',rstat,-1.0
          CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage )
        END IF
      END IF
    ELSE
      WRITE(6,'(A)')'DEBUG: Skipping write per debug setting '
    END IF

    IF (lhook) CALL dr_hook('IOS_STASH_SERVER:WRITE',zhook_out,zhook_handle)

  END SUBROUTINE IOS_stash_server_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Deallocate memory associated with a stash object
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE IOS_stash_server_deallocate(node)

    USE IOS_Server_Coupler, ONLY :                                             &
        procs
    IMPLICIT NONE

    TYPE(IOS_node_type),INTENT(INOUT) :: node
    INTEGER                           :: receive_buffer_size
    REAL(KIND=jprb)                   :: zhook_handle

    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:DEALLOCATE',zhook_in,zhook_handle)
    END IF

    IF (ASSOCIATED(node%distributed_data)) THEN

      receive_buffer_size=                                                     &
          getMaxFieldDomain()*                                                 &
          levels_in_pack(node)*                                                &
          procs

      DEALLOCATE(node%distributed_data)
      NULLIFY(node%distributed_data)

      ! Release space associated with receive buffers
      CALL IOS_Increment_Queue_Length                                          &
          (-1*receive_buffer_size*IOS_BytesPerReal,need_lock=.TRUE.)
    ELSE
      WRITE(iosStashServerMessage,'(A)')                                       &
          'Deallocation requested, but no allocation was present'
      CALL IOS_Ereport( RoutineName, -99, IosStashServerMessage )
    END IF

    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:DEALLOCATE',zhook_out,zhook_handle)
    END IF

  END SUBROUTINE IOS_stash_server_deallocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the receive tag to use for a given request
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER FUNCTION receive_tag(node,processor)

    IMPLICIT NONE
    TYPE(IOS_node_type), INTENT(IN) :: node
    INTEGER, INTENT(IN)             :: processor

    receive_tag=node%metadata%handle

  END FUNCTION receive_tag


  INTEGER FUNCTION levels_in_pack(node)
    IMPLICIT NONE
    TYPE(IOS_node_type),                                                       &
        INTENT(IN)      :: node
    INTEGER             :: num_records
    INTEGER             :: num_fields
    INTEGER             :: fieldtype
    INTEGER             :: i
    INTEGER             :: j
    INTEGER             :: ControlBlockLength
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:LEVELS_IN_PACK',zhook_in,zhook_handle)
    END IF

    j=1
    levels_in_pack=0
    num_records=0
    num_fields=0
    DO WHILE(j<node%metadata%data_size)
      IF (node%integer_data(j)==ios_stash_record_start) THEN!new record
        num_records=num_records+1
        ControlBlockLength=-1*node%integer_data(j+loc_record_len_control-1)

        !Advance to next record
        j=j+ios_stash_control_auto_len+ControlBlockLength
        num_fields=num_fields+1
      ELSE
        WRITE(iosStashServerMessage,'(A,I10)')                                 &
            'levels_in_pack: Control Block FAULT: j=',j
        CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,              &
            md=node%metadata,unit=node%metadata%unit)
      END IF
    END DO

    levels_in_pack=num_records

    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:LEVELS_IN_PACK',zhook_out,zhook_handle)
    END IF

  END FUNCTION levels_in_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Convenience function to return the total payload accross all
! processors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER FUNCTION total_payload_size(node)

    IMPLICIT NONE
    TYPE(IOS_node_type),                                                       &
        INTENT(IN)      :: node
    INTEGER             :: processor
    total_payload_size=0
    DO processor=0,atm_numprocs-1
      total_payload_size=total_payload_size+payload_size(node,processor)
    END DO
  END FUNCTION total_payload_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Figure out how many items there are in a package from a
! given cpu
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER FUNCTION payload_size(node,processor)

    USE IOS_geometry_utils

    IMPLICIT NONE

    TYPE(IOS_node_type),                                                       &
        INTENT(IN)      :: node
    INTEGER, INTENT(IN) :: processor
    INTEGER             :: j,num_records,num_fields
    INTEGER             :: FieldType
    INTEGER             :: i
    INTEGER             :: n,s,e,w
    INTEGER             :: FieldDomain
    INTEGER             :: subdomainType
    INTEGER             :: ControlBlockLength
    INTEGER             :: last
    TYPE(box)           :: domainBox
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:PAYLOAD_SIZE',zhook_in,zhook_handle)
    END IF

    j=1
    payload_size=0
    num_records=0
    num_fields=0
    last=0
    DO WHILE (j < node%metadata%data_size)
      IF (node%integer_data(j)==ios_stash_record_start) THEN!new record
        num_records=num_records+1
        ControlBlockLength=-1*node%integer_data(j+loc_record_len_control-1)

        !Advance to 1st element provided by user
        j=j+ios_stash_control_auto_len

        IF (node%integer_data(j)==ios_stash_distributed_field) THEN

          FieldType   = node%integer_data(j+loc_fld_type-1)
          FieldDomain = node%integer_data(j+loc_subdomain_flag-1)
          CALL IOS_unpack4(node%integer_data(j+loc_boundary-1),                &
              domainBox%n,                                                     &
              domainBox%s,                                                     &
              domainBox%e,                                                     &
              domainBox%w)

          IF (fieldDomain == IOS_partial_field) THEN
            subdomainType=                                                     &
                classify_subdomain(processor,FieldType,domainBox)
          ELSE
            subdomainType=complete_intersection
          END IF

          IF (FieldType < 1 .OR. FieldType > nfld_types) THEN
            WRITE(IosStashServerMessage,'(A,I20)')                             &
                'ERROR INVALID FIELD TYPE ',FieldType
            DO i=1,SIZE(node%integer_data)
              WRITE(6,'(I8,A,I8,A,I16)')i,'/',                                 &
                  SIZE(node%integer_data),' value: ',                          &
                  node%integer_data(i)
            END DO
            CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,          &
                md=node%metadata,unit=node%metadata%unit)
          ELSE
            num_fields=num_fields+1
            IF (subdomainType /= no_intersection                               &
                .OR. IOS_AsyncSendNull) THEN
              payload_size=payload_size+                                       &
                  size_map(1,FieldType,processor)*                             &
                  size_map(2,FieldType,processor)
              last=                                                            &
                  size_map(1,FieldType,processor)*                             &
                  size_map(2,FieldType,processor)
            ELSE
              last=0
            END IF
          END IF
        ELSE IF (node%integer_data(j)==ios_repeat_record) THEN
          payload_size=payload_size+last
        ELSE
          WRITE(iosStashServerMessage,'(A,I10)')                               &
              'Unknown data object in control record at ',j
          DO i=1,SIZE(node%integer_data)
            WRITE(6,'(I8,A,I8,A,I16)')                                         &
                i,'/',SIZE(node%integer_data),' value: ',                      &
                node%integer_data(i)
          END DO
          CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,            &
              md=node%metadata,unit=node%metadata%unit)
        END IF
        j=j+ControlBlockLength
      ELSE
        WRITE(iosStashServerMessage,'(A,I10)')                                 &
            'payload_size: Control Block FAULT: j=',j
        DO i=1,SIZE(node%integer_data)
          WRITE(6,'(I8,A,I8,A,I16)')i,'/',SIZE(node%integer_data),' value: ',  &
              node%integer_data(i)
        END DO
        CALL IOS_Ereport( RoutineName, 99, IosStashServerMessage,              &
            md=node%metadata,unit=node%metadata%unit)
      END IF
    END DO

    IF (lhook) THEN
      CALL dr_hook('IOS_STASH_SERVER:PAYLOAD_SIZE',zhook_out,zhook_handle)
    END IF

  END FUNCTION payload_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Log to disk what we just did
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE IOS_stash_server_log(node,FieldType,preprocess,                   &
      fullField,packing,s,n,w,e,                                               &
      MissingDataIndicator,DiskBlockSize,location,sz,                          &
      data_len)

    IMPLICIT NONE

    TYPE(IOS_node_type) :: node
    INTEGER,INTENT(IN)  :: DiskBlockSize!write data in lumps of this
    INTEGER,INTENT(IN)  :: n,s,e,w
    INTEGER,INTENT(IN)  :: fieldtype
    INTEGER,INTENT(IN)  :: packing
    INTEGER,INTENT(IN)  :: fullField
    INTEGER,INTENT(IN)  :: preprocess
    INTEGER,INTENT(IN)  :: location
    INTEGER,INTENT(IN)  :: sz
    INTEGER,INTENT(IN)  :: data_len
    REAL   ,INTENT(IN)  :: MissingDataIndicator
    REAL                :: t
    REAL(KIND=jprb)     :: zhook_handle

! DEPENDS ON : get_wallclock_time
    t=get_wallclock_time()-IOS_Start_Time
    IF (IOS_Verbosity>=IOS_PrStatus_Oper .AND. model_rank ==0) THEN
      WRITE(stashLogUnit,'(F8.2,2I4,I12,I10,I10,I5,2I4,4I5,I4)')               &
          t,node%metadata%handle,node%metadata%unit,                           &
          location,data_len*8,sz,DiskBlockSize,fieldtype,                      &
          fullField,s,n,w,e,packing
    END IF

  END SUBROUTINE IOS_stash_server_log

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Pull out the part of the field wanted for subdomaining
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE  IOS_stash_server_subdomain(global,fullField,                     &
      FieldType,n,s,e,w,x_range,y_range)

    IMPLICIT NONE

    REAL, POINTER         :: working(:,:)
    REAL, POINTER         :: global (:,:)
    INTEGER,INTENT(IN)    :: n,s,w
    INTEGER,INTENT(INOUT) :: e
    INTEGER,INTENT(IN)    :: fullField
    INTEGER,INTENT(IN)    :: FieldType
    INTEGER,INTENT(OUT)   :: x_range,y_range
    INTEGER               :: first_slice
    INTEGER               :: second_slice
    REAL(KIND=jprb)       :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_STASH_SERVER:SUBDOMAIN',zhook_in,zhook_handle)

    NULLIFY(working)
    x_range=0
    y_range=0

    IF (fullField==ios_full_field.AND.model_procs==1) THEN
      ! the subdomain is the same as the global field.
      x_range=atm_global_points(1,FieldType)
      y_range=atm_global_points(2,FieldType)

    ELSE ! For other cases we need to create a new object
         ! representing the subdomain

      ! Normalise our E/W ordering.
      IF (w > e) THEN
        e=e+atm_global_points(1,FieldType)
      END IF

      IF (s > n) THEN ! The subdomain is empty

        y_range=0
        x_range=0
        ALLOCATE(working(1:0,1:0)) ! a zero lengthed array

      ELSE IF (e > atm_global_points(1,FieldType)) THEN

        x_range=e-w+1
        y_range=n-s+1
        first_slice  = atm_global_points(1,FieldType)-w+1
        second_slice = e-atm_global_points(1,FieldType)
        ALLOCATE(working(1:x_range,1:y_range))
        working(1:first_slice,:)=                                              &
            global(w:atm_global_points(1,FieldType),s:n)
        working(first_slice+1:x_range,:)=                                      &
            global(1:e-atm_global_points(1,FieldType),s:n)

      ELSE

        x_range=e-w+1
        y_range=n-s+1
        ALLOCATE(working(1:x_range,1:y_range))
        working(:,:)=global(w:e,s:n)

      END IF

      ! As the subdomain is now allocated and populated, retarget
      ! the pointer at the correctly sized object.
      DEALLOCATE (global)
      global => working
      NULLIFY (working)

    END IF

    IF (lhook)CALL dr_hook('IOS_STASH_SERVER:SUBDOMAIN',zhook_out,zhook_handle)

  END SUBROUTINE IOS_stash_server_subdomain

END MODULE IOS_Stash_Server


