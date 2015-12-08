! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

! Discharge events from an IOS queue
MODULE IO_Server_Writer

  USE IO
  USE IOS_Common
  USE IOS_Queue_Mod
  USE IOS_Stash_Server, ONLY :        &
      IOS_Stash_Server_Init_Recvs,    &
      IOS_Stash_Server_Wait_For_Data, &
      IOS_Stash_server_Process,       &
      IOS_Stash_server_Deallocate
  USE IOS_Model_Geometry, ONLY : &
      IOS_DumpInitModel
  USE MPPIO_job_control_common, ONLY :&
      jc_write,    &
      jc_filename, &
      jc_unit
  USE MPPIO_file_utils, ONLY : &
      file_action, &
      file_op_pseudo_unit
  USE MPL, ONLY :  &
      MPL_INTEGER, &
      MPL_CHARACTER
  USE model_file,     ONLY : &
      model_file_open,       &
      model_file_close,      &
      model_file_managed,    &
      setRecordDiskLength,   &
      setRecordDiskStart,    &
      attachLookups,         &
      attachHeader,          &
      initLookups,           &
      setLookups,            &
      initHeader,            &
      FixedHeader,           &
      setHeader
  USE yomhook, ONLY: lhook, dr_hook  
  USE parkind1, ONLY: jprb, jpim
  USE PrintStatus_mod, ONLY : printStatus
  USE lookup_addresses

  IMPLICIT NONE

  ! Due to the asynchronous nature of the server the listener and 
  ! writer timestep independently. Whilst we do not make use of 
  ! timestepping in this version it is an important concept in many 
  ! coupling frameworks.
  INTEGER :: timestep

  ! params/vars  for dr_hook
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CONTAINS
  SUBROUTINE IOS_writer()

    IMPLICIT NONE

    TYPE(IOS_node_type), POINTER      :: node
    INTEGER, POINTER                  :: ipplook(:,:)
    INTEGER, POINTER                  :: fixed_header(:)
    INTEGER                           :: ierror
    INTEGER                           :: len_io
    INTEGER                           :: record_id
    LOGICAL                           :: done
    TYPE(fileState)                   :: fState
    REAL                              :: io_stat
    REAL                              :: t
    REAL                              :: t1
    REAL                              :: t2
    REAL, EXTERNAL                    :: get_wallclock_time
    REAL(KIND=jprb)                   :: zhook_handle
    CHARACTER (LEN=*), PARAMETER      :: RoutineNameWriter = 'IOS_Writer'
    CHARACTER (LEN=132)               :: Ios_Wrt_Message


    IF (lhook) CALL dr_hook('IO_SERVER_WRITER',zhook_in,zhook_handle)
    timestep=0
    NULLIFY(ipplook)
    NULLIFY(fixed_header)
    done=.FALSE.

    DO WHILE(.NOT.done)
      
      IF (IOS_getQueueItems() == 0) THEN
        CALL um_sleep(IOS_backoff_interval)
      ELSE 
        
        CALL IOS_Get_Node_From_Queue( node )
        IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
! DEPENDS ON : get_wallclock_time
          t1=get_wallclock_time()
          WRITE(6,'(A,I8,A,I3,A,A)')                    &
              'IOS: Info: Writer: transaction: ',       &
              node%transaction_number,' is for unit ',  &
              node%metadata%unit,' :',                  &
              TRIM(IOS_ActionName(node%metadata%action))
! DEPENDS ON: um_fort_flush
          CALL um_fort_flush(6,ierror)
        END IF

        IF (node%metadata%unit>0) THEN
          IF (node%metadata%action /= IOS_Action_StashWritePPData   &
              .AND.                                                 &
              node%metadata%action /= IOS_Action_StashWriteDumpData &
              ) THEN
            IF (global_rank/=io_server_for_unit(node%metadata%unit)) THEN
              WRITE(ios_wrt_message,'(A,I3,A,I6)')                 &
                  'Transaction for unit ',                         &
                  node%metadata%unit,                              &
                  ' came to me, but I think it should be for PE:', &
                  io_server_for_unit(node%metadata%unit)              
              CALL IOS_Ereport( RoutineNameWriter, 10,             &
                  ios_wrt_message, md=node%metadata )
            END IF
          END IF
        END IF

        SELECT CASE (node%metadata%action)

        CASE (IOS_Action_assign_unit)
          io_server_for_unit_lookup(node%metadata%intent)= &
              node%metadata%address
          WRITE(6,'(A,I3,A,I4)')                               &
              'IOS: Info: Writer: unit ',node%metadata%intent, &
              ' is assigned to IO server rank ',               &
              node%metadata%address
! DEPENDS ON: um_fort_flush
          CALL um_fort_flush(6,ierror)
!$OMP FLUSH

          !--------------------------------------------------------------
          ! Process
          !--------------------------------------------------------------
          
        CASE (IOS_Action_Process)
          timestep=timestep+1
! DEPENDS ON : get_wallclock_time
          t=get_wallclock_time()
          IF ( IOS_Verbosity >= IOS_PrStatus_Oper )             &
              WRITE(6,'(A,I4,A,F10.3)')                         &
              'IOS: Info: Writer: Entering Timestep ',timestep, &
              ' at time=',t-IOS_start_time

          !--------------------------------------------------------------
          ! Write out data (64 bit)
          !--------------------------------------------------------------

        CASE (IOS_Action_Write64, &           
            IOS_Action_Write_PP_Prepacked)
          
          IF (node%metadata%address >= 0 ) THEN 
            ! Metadata contains embedded setpos address
            IF (node%metadata%action == IOS_Action_Write_PP_Prepacked) THEN
              
              ! We are given a stash record in the address field. 
              ! So we should decode and record its location in the lookup table.

              record_id=node%metadata%address

              IF (record_id==1) THEN
                fixed_header=>attachHeader(node%metadata%unit,FixedHeader)
                node%metadata%address=Fixed_Header(160)-1
                NULLIFY(fixed_header)
              ELSE
                ipplook=>attachLookups(node%metadata%unit)
                node%metadata%address=&
                    ipplook(LBEGIN, record_id-1)+ &
                    ipplook(LBNREC, record_id-1)
                NULLIFY(ipplook)
              END IF

              IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN                   
                WRITE(6,'(A,I8)')                                   &
                    'IOS: Info: Writer: prepacked data record=',    &
                    record_id
                WRITE(6,'(A,I8)')                                   &
                    'IOS: Info: Writer:  decoded to disk address ', &
                    node%metadata%address
              END IF

              CALL setRecordDiskStart( &
                  node%metadata%unit,  &
                  record_id,           &
                  node%metadata % address)
            END IF ! Whether this was stash
            IF (IOS_Verbosity>=IOS_PrStatus_Diag)               &
                WRITE(6,'(A,I12)')                              &
                'IOS: Info: Writer: Process embedded setpos: ', &
                node%metadata%address
            
            CALL Setpos( node%metadata%unit, &
                node%metadata%address, ierror)
          END IF ! Whether the address field was positive

          IF (.NOT.ASSOCIATED(node%real_data)) THEN
            WRITE(ios_wrt_message,'(A)')&
                'FAULT: there is no data on this (Action_Write) node'
            CALL IOS_Ereport( RoutineNameWriter, 10, Ios_Wrt_Message, &
                md=node%metadata )
          END IF

          IO_Stat = -1.0
          Len_IO  = node%metadata%data_size
          
          CALL Buffout(node%metadata%unit, node%real_data, &
              node%metadata%data_size, Len_IO, IO_Stat )

          IF (Len_IO /= node%metadata%data_size .OR. IO_Stat/=-1.0) THEN
            WRITE(Ios_Wrt_Message,'(A,I4,I10,I10,F6.2)') &
                'Error in Writing to file on unit ',     &
                node%metadata%unit, Len_IO,              &
                node%metadata%data_size, IO_Stat
            CALL IOS_Ereport( RoutineNameWriter, 11, Ios_Wrt_Message, &
                md=node%metadata )
          END IF

          ! For stash ops. We were given a record number, so we need to 
          ! update the correct record with the length of the write           
          IF (node%metadata%action==IOS_Action_Write_PP_Prepacked) &
              CALL setRecordDiskLength( &
              node%metadata%unit,       &
              record_id,                &
              node%metadata % data_size)
          
          !--------------------------------------------------------------
          ! Write out data (32 bit)
          !--------------------------------------------------------------

        CASE (IOS_Action_Write32)
 
          IF (node%metadata%address >= 0 ) THEN 
            ! Metadata contains embedded setpos address

            IF (IOS_Verbosity>=IOS_PrStatus_Diag)               &
                WRITE(6,'(A,I12)')                              &
                'IOS: Info: Writer: Process embedded setpos: ', &
                node%metadata%address
            
            CALL Setpos( node%metadata%unit, &
                node%metadata%address, ierror)
          END IF
            
          IF (.NOT.ASSOCIATED(node%real32_data)) THEN
            WRITE(ios_wrt_message,'(A)')&
                'FAULT: there is no data on this (Action_Write32) node'
            CALL IOS_Ereport( RoutineNameWriter, 10, Ios_Wrt_Message, &
                md=node%metadata )
          END IF
          IO_Stat = -1.0
          Len_IO  = node%metadata%data_size
          
          CALL Buffo32(node%metadata%unit, node%real32_data, &
              node%metadata%data_size, Len_IO, IO_Stat )
          
          IF (Len_IO /= node%metadata%data_size .OR. IO_Stat/=-1.0) THEN
            WRITE(Ios_Wrt_Message,'(A,I4,I10,I10,F6.2)') &
                'Error in Writing to file on unit ',     &
                node%metadata%unit, Len_IO,              &
                node%metadata%data_size, IO_Stat
            CALL IOS_Ereport( RoutineNameWriter, 12, Ios_Wrt_Message, &
                md=node%metadata )
          END IF

          !--------------------------------------------------------------
          ! Read Data
          !--------------------------------------------------------------
        
        CASE (IOS_Action_Read64)

          IF (node%metadata%address >= 0 ) THEN 
            ! Metadata contains embedded setpos address
            IF (IOS_Verbosity>=IOS_PrStatus_Diag)               &
                WRITE(6,'(A,I12)')                              &
                'IOS: Info: Writer: Process embedded setpos: ', &
                node%metadata%address
            
            CALL Setpos( node%metadata%unit,                &
                node%metadata%address, ierror)
          END IF

          CALL IOS_Increment_Queue_Length                   &
              (node%metadata%data_size*IOS_BytesPerWord64,  &
              need_lock=.TRUE.)
          ALLOCATE(node%real_data( node%metadata%data_size) )
          CALL Buffin( node%metadata%unit,                  &
              node%real_data,                               &
              node%metadata%data_size,                      &
              len_IO, IO_stat )
          IF (.NOT. IOS_thread_0_calls_mpi) THEN
            CALL acquire_lock_mpi()
            IF (node%metadata%intent==IOS_Read_Broadcast) THEN
              CALL MPL_BCast(                              &
                  node%real_data,                          &
                  node%metadata%data_size*IOS_TUPerWord64, &
                  IOS_TUType,                              &
                  bcast_procs-1,                           &
                  IOS_BCast_Server_Comm,                   &
                  ierror)
            ELSE
              CALL MPL_Send(                                  &
                  node%real_data,                             &
                  node%metadata%data_size*IOS_TUPerWord64,    &
                  IOS_TUType, node%metadata%client,           &
                  Node%payloadTag, Global_Comm,               &
                  ierror)
            END IF
            CALL release_lock_mpi()
          ELSE
            IOS_ReadCompletionRequested=1
!$OMP FLUSH
            IF (IOS_Verbosity>=IOS_PrStatus_Oper) &
                WRITE(6,'(A)') &
                'Waiting for thread 0 to complete a request for read' 
            DO WHILE (IOS_ReadCompletionRequested==1)
              CALL um_sleep(IOS_backoff_interval)
!$OMP FLUSH 
            END DO
          END IF

        CASE  (IOS_Action_Read32 )

          IF (node%metadata%address >= 0 ) THEN 
            ! Metadata contains embedded setpos address
            IF (IOS_Verbosity>=IOS_PrStatus_Diag)               &
                WRITE(6,'(A,I12)')                              &
                'IOS: Info: Writer: Process embedded setpos: ', &
                node%metadata%address
            
            CALL Setpos( node%metadata%unit,                &
                node%metadata%address, ierror)
          END IF

          CALL IOS_Increment_Queue_Length                   &
              (node%metadata%data_size*IOS_BytesPerWord32,  &
              need_lock=.TRUE.)
          ALLOCATE(node%real32_data( node %metadata%data_size) )
          CALL Buffin32( node%metadata%unit,                &
              node%real32_data, node%metadata%data_size,    &
              len_IO, IO_stat)
          IF (.NOT. IOS_thread_0_calls_mpi) THEN
            CALL acquire_lock_mpi()
            IF (node%metadata%intent==IOS_Read_Broadcast) THEN
              CALL MPL_BCast(                              &
                  node%real32_data,                        &
                  node%metadata%data_size*IOS_TUPerWord32, &
                  IOS_TUType,                              &
                  bcast_procs-1,                           &
                  IOS_BCast_Server_Comm,                   &
                  ierror)
            ELSE
              CALL MPL_Send(node%real32_data,                 &
                  node%metadata%data_size*IOS_TUPerWord32,    &
                  IOS_TUType,                                 &
                  node%metadata%client,                       &
                  Node%payloadTag, Global_Comm, ierror)
            END IF
            CALL release_lock_mpi()
          ELSE
            IOS_ReadCompletionRequested=1
!$OMP FLUSH 
            IF (IOS_Verbosity>=IOS_PrStatus_Oper) &
                WRITE(6,'(A)')                    &
                'Waiting for thread 0 to complete a request for read'
            DO WHILE (IOS_ReadCompletionRequested==1)
              CALL um_sleep(IOS_backoff_interval)
!$OMP FLUSH 
            END DO
          END IF

        CASE (IOS_Action_Read64_Integer)

          IF (node%metadata%address >= 0 ) THEN 
            ! Metadata contains embedded setpos address
            IF (IOS_Verbosity>=IOS_PrStatus_Diag)               &
                WRITE(6,'(A,I12)')                              &
                'IOS: Info: Writer: Process embedded setpos: ', &
                node%metadata%address
            
            CALL Setpos( node%metadata%unit,                &
                node%metadata%address, ierror)
          END IF

          CALL IOS_Increment_Queue_Length                   &
              (node%metadata%data_size*IOS_BytesPerWord64,  &
              need_lock=.TRUE.)
          ALLOCATE(node%integer_data( node%metadata%data_size) )
          CALL Buffin( node%metadata%unit,                  &
              node%integer_data,                            &
              node%metadata%data_size,  len_IO, IO_stat )
          IF (.NOT. IOS_thread_0_calls_mpi) THEN
            CALL acquire_lock_mpi()
            IF (node%metadata%intent==IOS_Read_Broadcast) THEN
              CALL MPL_BCast(                              &
                  node%integer_data,                       &
                  node%metadata%data_size*IOS_TUPerWord64, &
                  IOS_TUType,                              &
                  bcast_procs-1,                           &
                  IOS_BCast_Server_Comm,                   &
                  ierror)
            ELSE
              CALL MPL_Send(node%integer_data,                &
                  node%metadata%data_size*IOS_TUPerWord64,    &
                  IOS_TUType,                                 &
                  node%metadata%client, Node%payloadTag,      &
                  Global_Comm, ierror)
            END IF
            CALL release_lock_mpi()
          ELSE
            IOS_ReadCompletionRequested=1
            IF (IOS_Verbosity>=IOS_PrStatus_Oper) &
                WRITE(6,'(A)') &
                'Waiting for thread 0 to complete a request for read'
!$OMP FLUSH
            DO WHILE (IOS_ReadCompletionRequested==1)
              CALL um_sleep(IOS_backoff_interval)
!$OMP FLUSH
            END DO
          END IF

          !--------------------------------------------------------------
          ! Enquiry
          !--------------------------------------------------------------

        CASE (IOS_Action_Enquire)

          IF (node%metadata%address >= 0 ) THEN 
            ! Metadata contains embedded setpos address
            IF (IOS_Verbosity>=IOS_PrStatus_Diag)               &
                WRITE(6,'(A,I12)')                              &
                'IOS: Info: Writer: Process embedded setpos: ', &
                node%metadata%address
            
            CALL Setpos( node%metadata%unit, &
                node%metadata%address, ierror)
          END IF

          CALL ioFileState(node%metadata%unit,fState)
          ios_unit_status%extent   = fstate%fileExtent
          ios_unit_status%position = fstate%filePosition
          IF (.NOT. IOS_thread_0_calls_mpi) THEN
            CALL acquire_lock_mpi()
            CALL MPL_Send(ios_unit_status,              &
                IOS_State_Size*IOS_TUPerWord64,         &
                IOS_TUType,                             &
                node%metadata%client, Node%payloadTag,  &
                Global_Comm, ierror)
            CALL release_lock_mpi()
          ELSE
            IOS_EnqCompletionRequested=1
!$OMP FLUSH
            DO WHILE (IOS_EnqCompletionRequested==1)
              CALL um_sleep(IOS_backoff_interval)
!$OMP FLUSH
            END DO
          END IF


          !--------------------------------------------------------------
          ! Open File
          !--------------------------------------------------------------

        CASE (IOS_Action_Open)

          IF (node%metadata%address == ioFileTypeUM) THEN

            CALL Model_File_Open(node%metadata%unit,     &
                node%metadata%string,                    &
                node%metadata%name_length,               &
                node%metadata%intent,                    &
                ioNameProvided,                          &
                ierror,                                  &
                ioLocality=ioAllLocal,                   &
                fileType=node%metadata%address)
            
          ELSE

            CALL File_Open(node%metadata%unit,     &
                node%metadata%string,              &
                node%metadata%name_length,         &
                node%metadata%intent,              &
                ioNameProvided,                    &
                ierror,                            &
                ioLocality=ioAllLocal,             &
                fileType=node%metadata%address)
            
          END IF

          IF (ierror /= 0) THEN
            WRITE(Ios_Wrt_Message,'(A,I3)')       &
                'Error in Opening file on unit ', &
                node%metadata%unit
            CALL IOS_Ereport( RoutineNameWriter, 10, Ios_Wrt_Message, &
                md=node%metadata )
          END IF

          !--------------------------------------------------------------
          ! Close File
          !--------------------------------------------------------------

        CASE (IOS_Action_Close)

          IF (is_unit_open(node%metadata%unit)) THEN! Unit is Open

            IF (model_file_managed(node%metadata%unit)) THEN
              CALL Model_File_Close(node%metadata%unit,   &
                  node%metadata%string,             &
                  node%metadata%name_length,        &
                  1,                                &
                  node%metadata%delete,             &
                  ierror )
            ELSE
              CALL File_Close(node%metadata%unit,   &
                  node%metadata%string,             &
                  node%metadata%name_length,        &
                  1,                                &
                  node%metadata%delete,             &
                  ierror )
            END IF

            IF (ierror /= 0) THEN
              WRITE(Ios_Wrt_Message,'(A,I3)')       &
                  'Error in Closing File on unit ', &
                  node%metadata%unit
              CALL IOS_Ereport( RoutineNameWriter, 10, Ios_Wrt_Message, &
                  md=node%metadata )
            END IF
          ELSE
            IF (IOS_Verbosity>=IOS_PrStatus_Normal)         &
                WRITE(6,'(A,I3,A)')                         &
                'Tried to close unit (',node%metadata%unit, &
                ') which is not open'
          END IF          

          !--------------------------------------------------------------
          ! Setpos
          !--------------------------------------------------------------

        CASE (IOS_Action_Setpos)
          IF (Is_unit_open(node%metadata%unit)) THEN! Unit is Open
            CALL Setpos( node%metadata%unit, &
                node%metadata%address, ierror)
            IF (ierror /= 0) THEN
              WRITE(Ios_Wrt_Message,'(A,I3)') &
                  'Error in setpos on unit ', &
                  node%metadata%unit
              CALL IOS_Ereport( RoutineNameWriter, 10, Ios_Wrt_Message, &
                  md=node%metadata )
            END IF
          ELSE
            WRITE(Ios_Wrt_Message,'(A,I3)')        &
                'Setpos called on non-open unit:', &
                node%metadata%unit
            CALL IOS_Ereport( RoutineNameWriter, -999, Ios_Wrt_Message, &
                md=node%metadata )
          END IF

          !-------------------------------------------------------------------
          ! sync Ops
          !-------------------------------------------------------------------
          
        CASE (IOS_Action_Sync)
          CALL Sync_single(node%metadata%unit,ierror)

        CASE (IOS_Action_Sync_Barrier)
          CALL Sync_single(node%metadata%unit,ierror)
          IF (.NOT. IOS_thread_0_calls_mpi) THEN
            ! Send 1 integer to the client to indicate completion
            CALL MPL_Send(node%metadata,                &
              1 , MPL_INTEGER, node%metadata%client,    &
              Node%payloadTag, Global_comm,             &
              ierror)
          END IF

          !--------------------------------------------------------------
          ! Pipe Ops, open, close and release
          !--------------------------------------------------------------

        CASE (IOS_Action_Open_Pipe)
          jc_filename(1:node%metadata%name_length) &
              =node%metadata%string(1:node%metadata%name_length)
          jc_filename(node%metadata%name_length+1:)=' '
          IF (IOS_Verbosity>=IOS_PrStatus_Oper) &
              WRITE(6,'(A,A)')'JOB CONTROL: OPEN PIPE:',TRIM(jc_filename)
          OPEN(unit=jc_unit,file=jc_filename,iostat=ierror, action="write")
          IF (ierror /= 0) THEN
            WRITE(Ios_Wrt_Message,'(A)') 'JOB CONTROL: ERROR OPENING PIPE'
            CALL IOS_Ereport( RoutineNameWriter, ierror, Ios_Wrt_Message, &
                md=node%metadata )
          END IF
          
        CASE (IOS_Action_Close_Pipe)
          IF (IOS_Verbosity>=IOS_PrStatus_Oper) & 
              WRITE(6,'(A,A)')'JOB CONTROL: CLOSE PIPE:',TRIM(jc_filename)
          CLOSE(unit=jc_unit,iostat=ierror)
          IF (ierror /= 0) THEN
            WRITE(Ios_Wrt_Message,'(A)')'JOB CONTROL: ERROR CLOSING PIPE'
            CALL IOS_Ereport( RoutineNameWriter, ierror, Ios_Wrt_Message, &
                md=node%metadata )            
          END IF
          

        CASE (IOS_Action_Release)
          t=get_wallclock_time()
          IF (IOS_Verbosity>=IOS_PrStatus_Oper)                    &  
              WRITE(6,'(A,F10.3)')                                 &
              'IOS: Info: Writer: Release op, syncing at t=',      &
              t-IOS_Start_time
          IF (.NOT. IOS_thread_0_calls_mpi) THEN
            CALL acquire_lock_mpi()
            CALL IOS_SoftSync(io_server_for_unit(jc_unit))
            CALL release_lock_mpi()
          END IF
          t=get_wallclock_time()
          IF (IOS_Verbosity>=IOS_PrStatus_Oper)                    &  
              WRITE(6,'(A,F10.3)')                                 &
              'IOS: Info: Writer: Release op, syncing done at t=', &
              t-IOS_Start_time
          IF (global_rank==io_server_for_unit(jc_unit))            &
              CALL jc_write(node%metadata%intent,                  &
              node%metadata%string)
          
          !--------------------------------------------------------------
          ! File operations
          !--------------------------------------------------------------

        CASE (IOS_Action_FileOp)
          t=get_wallclock_time()
          IF (IOS_Verbosity>=IOS_PrStatus_Oper)                    &  
              WRITE(6,'(A,F10.3)')                                 &
              'IOS: Info: Writer: FileOp, syncing at t=',          &
              t-IOS_Start_time
          IF (.NOT. IOS_thread_0_calls_mpi) THEN
            CALL acquire_lock_mpi()
            CALL IOS_SoftSync(io_server_for_unit(file_op_pseudo_unit))
            CALL release_lock_mpi()
          END IF
          t=get_wallclock_time()
          IF (IOS_Verbosity>=IOS_PrStatus_Oper)                    &  
              WRITE(6,'(A,F10.3)')                                 &
              'IOS: Info: Writer: FileOp, syncing done at t=',     &
              t-IOS_Start_time
          IF (global_rank==io_server_for_unit(file_op_pseudo_unit)) THEN
            CALL file_action(TRIM(node%metadata%string))
          END IF

          !--------------------------------------------------------------
          ! Fence Ops
          !--------------------------------------------------------------
          
        CASE (IOS_Action_Fence)

          IF (.NOT. IOS_thread_0_calls_mpi) THEN

            t=get_wallclock_time()
            IF (IOS_Verbosity>=IOS_PrStatus_Oper)         &
                WRITE(6,'(A,I4,A,F10.3)')                 &
                'IOS: Info: Writer: Fence: server=',      &
                node%metadata%address,' at ',             &
                t-IOS_Start_time
            
            CALL acquire_lock_mpi()
            CALL IOS_SoftSync(node%metadata%address)
            CALL release_lock_mpi()
            
            t=get_wallclock_time()
            IF (IOS_Verbosity>=IOS_PrStatus_Oper)         &
                WRITE(6,'(A,F10.3)')                      &
                'IOS: Info: Writer: Fence completed at ', &
                t-IOS_Start_time
          ELSE
            IF (IOS_Verbosity>=IOS_PrStatus_Normal)       &
                WRITE(6,'(A,A)')'IOS: Info: Writer:',     &
                ' Fence occured previously due to MPI semantics'
          END IF

          !--------------------------------------------------------------
          ! Finish - close down server
          !--------------------------------------------------------------

        CASE (IOS_Action_Finish)

! DEPENDS ON : get_wallclock_time
          t=get_wallclock_time()
          WRITE(6,'(A)') &
              '********************************************************'
          WRITE(6,'(A,F10.3,A)') &
              '* IO SERVER:WRITER   RECEIVED FINISH CMD AT ', &
              t-IOS_Start_Time,' *'
          WRITE(6,'(A)') &
              '********************************************************'
          done=.TRUE.

          !-------------------------------------------------------------------
          ! Process a setpos. For Stash setpos we are given a 
          ! unit/record from which to compute a seek address.
          !-------------------------------------------------------------------
          
        CASE (IOS_Action_StashSetPos)
          
          IF( Is_unit_open( node%metadata % unit ) ) THEN
            fixed_header=>attachHeader(node%metadata%unit,FixedHeader)
            ipplook=>attachLookups(node%metadata%unit)
            IF (node%metadata%address==1) THEN
              node%metadata%address=Fixed_Header(160)-1
            ELSE
              node%metadata%address=&
                  IPPLOOK(LBEGIN, node%metadata%address-1)+ &
                  IPPLOOK(LBNREC, node%metadata%address-1)
            END IF
            
            CALL Setpos( node%metadata % unit,       &
                node%metadata % address, ierror)
            IF (ierror /= 0) THEN
              WRITE(Ios_Wrt_Message,'(A,I8)') 'Error in setpos on unit ', &
                  node%metadata % unit
              CALL IOS_Ereport( RoutineNameWriter, 50, Ios_Wrt_Message, &
                  md=node%metadata )
            END IF
          ELSE
            WRITE(Ios_Wrt_Message,'(A,I8)') 'Setpos called on non-open unit:', &
                node%metadata % unit
            CALL IOS_Ereport( RoutineNameWriter, -999, Ios_Wrt_Message, &
                  md=node%metadata )
          END IF
          NULLIFY(ipplook)
          NULLIFY(fixed_header)
          
          !-------------------------------------------------------------------
          ! Process The merging of PP lookup tables
          !-------------------------------------------------------------------
          
        CASE (IOS_ACTION_MergePPLookup)
          
          IF (.NOT.ASSOCIATED(node%integer_data)) THEN
            WRITE(ios_wrt_message,'(A)') &
                'FAULT: there is no data on this (Action_WritePPLook) node'
            CALL IOS_Ereport( RoutineNameWriter, 10, Ios_Wrt_Message, &
                md=node%metadata )
          END IF
          
          CALL setLookups(node%metadata%unit, node%integer_data)
          IF (IOS_Verbosity>=IOS_PrStatus_Oper)    &
              WRITE(6,'(A,A,I3)')'IOS: Writer: ',  &
              'Merged lookup records for unit ',   &
              node%metadata % unit          
          
          !-------------------------------------------------------------------
          ! Process The writing of PP lookup tables
          !-------------------------------------------------------------------
          
        CASE (IOS_ACTION_StashWritePPLookup)
          IF (.NOT.ASSOCIATED(node%integer_data)) THEN
            WRITE(ios_wrt_message,'(A)') &
                'FAULT: there is no data on this (Action_WritePPLook) node'
            CALL IOS_Ereport( RoutineNameWriter, 10, Ios_Wrt_Message, &
                md=node%metadata )
          END IF
          
          IF (node%metadata%address >= 0 ) THEN 
            ! Metadata contains a setpos address
            IF (IOS_Verbosity>=IOS_PrStatus_Diag)               &
                WRITE(6,'(A,I12)')                              &
                'IOS: Info: Writer: Process embedded setpos: ', &
                node%metadata%address
            
            CALL Setpos( node%metadata%unit, &
                node%metadata%address, ierror)
          END IF
          
          IO_Stat = -1.0
          Len_IO  = node%metadata % data_size
               
          CALL setLookups(node%metadata % unit,node%integer_data)
          IF (IOS_Verbosity>=IOS_PrStatus_Oper)    &
              WRITE(6,'(A,A,I3)')'IOS: Writer: ',  &
              'Merged lookup records for unit ',   &
              node%metadata % unit
          ipplook=>attachLookups(node%metadata % unit)

          CALL Buffout(node%metadata % unit,ipplook,     &
              node%metadata % data_size, Len_IO, IO_Stat )
          
          IF (Len_IO /= node%metadata % data_size .OR. &
              IO_Stat /= -1.0) THEN
            WRITE(Ios_Wrt_Message,'(A,I4,I10,I10,F6.2)') &
                'Error in Writing to file on unit ',     &
                node%metadata % unit, Len_IO,            &
                node%metadata % data_size,               &
                IO_Stat
            CALL IOS_Ereport( RoutineNameWriter, 10, Ios_Wrt_Message, &
                md=node%metadata )
          END IF
          
          NULLIFY(ipplook)

          !-------------------------------------------------------------------
          ! Process The initialisation of PP lookup buffers
          !-------------------------------------------------------------------

        CASE(IOS_Action_StashInitPPLookup)
          IF (.NOT.ASSOCIATED(node%integer_data)) THEN
            WRITE(ios_wrt_message,'(A)')&
                'IOS: FAULT: no data in Action_StashInitPPLookup node'
            CALL IOS_Ereport( RoutineNameWriter, 10, Ios_Wrt_Message, &
                md=node%metadata )
          END IF
          
          CALL InitLookups(ipplook,  &
              node%integer_data(1),  &
              node%integer_data(2),  &
              node%integer_data(3),  &
              node%integer_data(4),  &
              node%integer_data(5))
          NULLIFY(ipplook)

          !-------------------------------------------------------------------
          ! Process The initialisation of a fixed header
          !-------------------------------------------------------------------
          
        CASE(IOS_Action_StashInitHeader)
          CALL initHeader(node%metadata % unit, FixedHeader, &
              node%metadata%address)
          
          !-------------------------------------------------------------------
          ! Process The setting of a fixed header
          !-------------------------------------------------------------------

        CASE(IOS_Action_StashSetHeader)
          IF (.NOT.ASSOCIATED(node%integer_data)) THEN
            WRITE(ios_wrt_message,'(A)') &
                          'FAULT: no data in Action_StashSetFXH node'
            CALL IOS_Ereport( RoutineNameWriter, 10, Ios_Wrt_Message, &
                md=node%metadata )
          END IF
          CALL setHeader(node%metadata % unit, FixedHeader, node%integer_data)
          
          !--------------------------------------------------------------
          ! Dump initialisation
          !--------------------------------------------------------------

       CASE (IOS_Action_DumpInitModel)
         CALL IOS_DumpInitModel( &
             node%real_data,     &
             node%metadata%data_size )
          
          !--------------------------------------------------------------
          ! Write out pp data
          !--------------------------------------------------------------
          
        CASE (IOS_Action_StashWritePPData)
          CALL IOS_Stash_Server_Wait_For_Data(node)
          CALL IOS_Stash_server_process(node)
          CALL IOS_Stash_server_deallocate(node)

          !--------------------------------------------------------------
          ! Write out dump data
          !--------------------------------------------------------------
                    
       CASE (IOS_Action_StashWriteDumpData)
         CALL IOS_Stash_Server_Wait_For_Data(node)
         CALL IOS_Stash_server_process(node)
         CALL IOS_Stash_server_deallocate(node)
         
          !--------------------------------------------------------------
          ! Modify IOS behaviour
          !--------------------------------------------------------------
          
        CASE (IOS_Action_Config)
          SELECT CASE (node%metadata%intent)
            
          CASE ( IOS_PrStatus_Min, &
              IOS_PrStatus_Normal, &
              IOS_PrStatus_Oper,   &
              IOS_PrStatus_Diag )

            IF (IOS_acquire_model_prsts) THEN
                IOS_Verbosity=node%metadata%intent
                printstatus  =node%metadata%intent
                IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
                  WRITE(6,'(A,I2)')                               &
                      'IOS: Info: Writer: Setting verbosity to ', &
                      node%metadata%intent
                END IF
              ELSE
                IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
                  WRITE(6,'(A)')                                  &
                      'IOS: Info: Writer: Setting verbosity ignored.'
                END IF
            END IF
            
          CASE ( IOS_Timer_Activate )

            CALL io_timing_init()




          CASE DEFAULT
            WRITE(Ios_Wrt_Message,'(A,I8)') &
                'Configuration option not recognised:',   &
                node%metadata%intent
            CALL IOS_Ereport( RoutineNameWriter, -10, Ios_Wrt_Message, &
                md=node%metadata )
            
          END SELECT
        
        !--------------------------------------------------------------
        ! Unrecognised action
        !--------------------------------------------------------------

        CASE DEFAULT
          WRITE(Ios_Wrt_Message,'(A,I8)') &
              'Action not recognised:',   &
              node%metadata%action
          CALL IOS_Ereport( RoutineNameWriter, 10, Ios_Wrt_Message, &
              md=node%metadata )
        END SELECT

        IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
! DEPENDS ON : get_wallclock_time
          t2=get_wallclock_time()
          WRITE(6,'(A,I8,A,F8.3)')                         &
              'IOS: Info: Writer: transaction: ',          &
              node%transaction_number,' is completed in ', &
              t2-t1
! DEPENDS ON: um_fort_flush
          CALL um_fort_flush(6,ierror)
        END IF
        CALL IOS_Remove_Data_From_Queue()
      END IF
    END DO

    WRITE(6,'(A)') 'IOS: Info: Writer: Process closing'
! DEPENDS ON: um_fort_flush
    CALL um_fort_flush(6,ierror)
    IF (lhook) CALL dr_hook('IO_SERVER_WRITER',zhook_out,zhook_handle)

  END SUBROUTINE IOS_Writer
END MODULE IO_Server_Writer
