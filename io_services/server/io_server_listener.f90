! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

MODULE IO_Server_Listener
  
  USE MPPIO_Job_Control,      ONLY : &
      jc_unit
  USE IOS_Common  
  USE IOS_Queue_Mod
  USE IOS_Stash_common
  USE IOS_Stash_Server,       ONLY : &
      IOS_Stash_Server_Init_Recvs,   &
      levels_in_pack
  USE IOS_Model_Geometry,     ONLY : &
      IOS_Server_Geometry_Init,      &
      getMaxFieldDomain,             &
      atm_numprocs
  USE IOS_Server_Coupler, ONLY :     &
      IOS_Server_Coupler_Init,       &
      procs
  USE MPL
  USE yomhook,                ONLY : &
      lhook,                         &
      dr_hook  
  USE parkind1,               ONLY : &
      jprb,                          &
      jpim

  IMPLICIT NONE

  CHARACTER (LEN=132), PRIVATE :: IOS_listener_message  
  INTEGER                      :: timestep
! params/vars  for dr_hook
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CONTAINS

  LOGICAL FUNCTION ASSERT_PERMITTED_CLIENT(metadata)
    IMPLICIT NONE
    TYPE(IOS_metadata_type),INTENT(IN)      :: metadata
    IF ( metadata%client /= 0 ) &
        CALL IOS_ereport( 'IOS_listener',999,&
        'client rank !=0 not permitted',md=metadata) 
    ASSERT_PERMITTED_CLIENT=.TRUE.
  END FUNCTION ASSERT_PERMITTED_CLIENT

  SUBROUTINE IOS_Listener()
    IMPLICIT NONE

    CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_Listener'
    INTEGER                      :: ierror
    INTEGER                      :: status(MPL_STATUS_SIZE)
    INTEGER                      :: command_request !an MPL request obj
    INTEGER                      :: express_command_request !an MPL request obj
    INTEGER                      :: reserveSpace
    TYPE(IOS_metadata_type)      :: metadata
    TYPE(IOS_metadata_type)      :: express_metadata
    TYPE(IOS_node_type), POINTER :: node
    INTEGER, ALLOCATABLE         :: idata(:)
    INTEGER                      :: mdTag
    INTEGER                      :: exTag
    INTEGER                      :: strTag
    INTEGER                      :: payloadTag
    INTEGER                      :: tagOffset
    LOGICAL                      :: have_message
    LOGICAL                      :: OK
    LOGICAL                      :: process_express_item
    LOGICAL                      :: process_normal_item
    INTEGER                      :: sequenceID
    REAL                         :: t
    REAL, EXTERNAL               :: get_wallclock_time
    REAL(KIND=jprb):: zhook_handle1
    REAL(KIND=jprb):: zhook_handle2


    IF ( lhook ) CALL dr_hook('IOS_LISTENER',zhook_in,zhook_handle1)

    WRITE(6,'(A)') 'IOS: Info: Listener: Process started'
    timestep=0
    sequenceID=1;



    metadata%action = IOS_Action_Unset

! Start these flags as true so we post their initial receives
    process_express_item=.TRUE.
    process_normal_item=.TRUE.

    ! Keep listening until we get told to stop
    DO WHILE (metadata%action /= IOS_Action_Finish)
      NULLIFY(node)

! Calculate tags for this loop
      tagOffset=MOD(sequenceID,IOS_Request_Tag_Gap)
      mdTag      =IOS_Request_Tag_Base    +tagOffset
      exTag      =IOS_Request_Tag_Express 
      strTag     =IOS_Request_Tag_Str     +tagOffset
      payloadTag =IOS_Request_Tag_Payload +tagOffset

      !----------------------------------------------------------------
      ! Receive metadata. This is complicated because this thread may 
      ! not block. If he does, he may not notice that thread 1 is 
      ! asking for assistance in the case that a single threaded MPI 
      ! is in use, which would result in a deadlock.
      !----------------------------------------------------------------
      CALL acquire_lock_mpi()

! If we have not got an outstanding receive for an express item post one
      IF (process_express_item) THEN
        CALL MPL_iRecv(express_metadata , IOS_MD_Len, MPL_INTEGER, &
            MPL_ANY_SOURCE, exTag, Global_Comm,                    &
            express_command_request, ierror)   
        process_express_item=.FALSE.
      END IF

! If we have not got an outstanding receive for an regular item post one
      IF (process_normal_item) THEN
        CALL MPL_iRecv(metadata , IOS_MD_Len, MPL_INTEGER, &
            MPL_ANY_SOURCE, mdTag, Global_Comm,            &
            command_request, ierror)   
        process_normal_item=.FALSE.
      END IF
      CALL release_lock_mpi()

! Now wait for something to arrive
      have_message=.FALSE.
      IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
        WRITE(6,'(A,I6)')'IOS: Info: Listener: Waiting for operation ', &
            sequenceID
      END IF
      DO WHILE (.NOT.have_message)
        CALL acquire_lock_mpi()
! First check the express tag (note we don't release MPI - its express!)
        CALL MPL_Test(express_command_request, have_message, status, ierror)

! If we have an express message then we need to set the flag to post another
        IF (have_message) THEN
          process_express_item=.TRUE.

! Otherwise check on normal messages (and release MPI)
        ELSE 
          CALL MPL_Test(command_request, have_message, status, ierror)
          CALL release_lock_mpi()

! If we have a normal message then we need to set the flag to post another
! and update the sequence ID
          IF (have_message) THEN
            process_normal_item=.TRUE.
            sequenceID=sequenceID+1
          ELSE 
            IF ( IOS_thread_0_calls_mpi ) THEN
              !Check to see if thread 1 needs attention    
              CALL IOS_Aux_Thread1_Assist()
            END IF
          END IF
        END IF
        IF ( .NOT.have_message ) CALL um_sleep(IOS_backoff_interval)
      END DO

      IF (process_express_item) THEN
        IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
          WRITE(6,'(A,A)')                    &
              'IOS: Info: Listener: Received High Priority Action: ', &
              TRIM(IOS_ActionName(express_metadata%action))
! DEPENDS ON: um_fort_flush
          CALL um_fort_flush(6,ierror)
        END IF
      ELSE
        IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
          WRITE(6,'(A,A)')                    &
              'IOS: Info: Listener: Received a transaction: ', &
              TRIM(IOS_ActionName(metadata%action))
! DEPENDS ON: um_fort_flush
          CALL um_fort_flush(6,ierror)
        END IF
      END IF

      OK=assert_permitted_client(metadata)
      ! The client filled in the detail with his model id. 
      ! We may need to send back data, so we need to switch 
      ! that to his global ID after we checked it.
      metadata%client = status(MPL_SOURCE)
      
      ! For debugging we may want to pause here until 
      ! the writer isn't busy this will allow everything 
      ! to be processed in natural sequence with no queueing
      IF ( serialize_all_ops ) THEN
        CALL IOS_WaitForQueueDrain()
      END IF

!
! The remainder of this routine is a switch box to deal with inbound
! message types
!

      IF (process_express_item) THEN
        SELECT CASE (express_metadata%action)
          
        CASE (IOS_Action_LoadStatus)
 
        !--------------------------------------------------------------
        ! Process the request for my loading
        !--------------------------------------------------------------

          ALLOCATE( idata ( loadBalanceDataSize ) )
          idata(loadBalanceQueueLen)  = IOS_getQueueItems()
          idata(loadBalanceQueueData) = IOS_getQueuePayload()
          CALL MPL_Send(idata,                       &
              loadBalanceDataSize*IOS_TUPerWord64,   &
              IOS_TUType, express_metadata%client,   &
              exTag, Global_Comm, ierror)
          DEALLOCATE( idata  )
          
        CASE DEFAULT
          WRITE(IOS_listener_message,'(A,I3,A)') 'Express Action ', &
              express_metadata%action,' not recognised!' 
          CALL IOS_ereport( RoutineName, 60, IOS_listener_message,  &
              md=express_metadata )
          
        END SELECT

        ! Release the MPI lock from when we checked for the express message
        CALL release_lock_mpi()
        
      ELSE
        SELECT CASE (metadata%action)
          
        CASE (IOS_Action_Assign_Unit)
          
          !--------------------------------------------------------------
          ! Process the assignment of a unit
          !--------------------------------------------------------------
          
          node => make_new_node()
          node%metadata = metadata
          CALL IOS_Put_Node_In_Queue(node)

          !--------------------------------------------------------------
          ! Process Writes (64 bit/stash)
          !--------------------------------------------------------------

        CASE (IOS_Action_Write64,  &
            IOS_Action_Write_PP_Prepacked)

          CALL IOS_WaitForFreeSpace(metadata%data_size* &
              IOS_BytesPerWord64)

          node => make_new_node()
          node%metadata = metadata

          CALL IOS_Increment_Queue_Length                             &
              (metadata%data_size*IOS_BytesPerWord64)
          ALLOCATE( node%real_data( metadata%data_size) )

          CALL acquire_lock_mpi()
          CALL MPL_Recv(node%real_data, metadata%data_size, MPL_REAL, &
              metadata%client, payloadTag , Global_Comm,              &
              status, ierror)
          CALL release_lock_mpi()
          ! Put data on Queue
          CALL IOS_Put_Node_In_Queue(node)

          !--------------------------------------------------------------
          ! Process Writes (32 bit)
          !--------------------------------------------------------------

        CASE (IOS_Action_Write32)

          CALL IOS_WaitForFreeSpace(metadata%data_size* &
              IOS_BytesPerWord32)

          node => make_new_node()
          node%metadata = metadata

          CALL IOS_Increment_Queue_Length &
              (metadata%data_size*IOS_BytesPerWord32)
          ALLOCATE( node%real32_data( metadata%data_size) )

          CALL acquire_lock_mpi()
          CALL MPL_Recv(node%real32_data, metadata%data_size, MPL_REAL, &
              metadata%client, PayloadTag, Global_Comm,                 &
              status, ierror)
          CALL release_lock_mpi()
          ! Put data on Queue
          CALL IOS_Put_Node_In_Queue(node)

          !--------------------------------------------------------------
          ! Process Dump initialisation
          !--------------------------------------------------------------

        CASE ( IOS_Action_DumpInitModel )

          CALL IOS_WaitForFreeSpace(metadata%data_size* &
              IOS_BytesPerWord64)

          node => make_new_node()
          node%metadata = metadata

          CALL IOS_Increment_Queue_Length(metadata%data_size* &
              IOS_BytesPerInteger)
          ALLOCATE( node%real_data( metadata%data_size) )

          CALL acquire_lock_mpi()
          CALL MPL_Recv(node%real_data, metadata%data_size, MPL_REAL, &
              metadata%client, PayloadTag, Global_Comm,               &
              status, ierror)
          CALL release_lock_mpi()
          ! Put data on Queue
          CALL IOS_Put_Node_In_Queue(node)

          !--------------------------------------------------------------
          ! Process Reads
          !--------------------------------------------------------------

        CASE (IOS_Action_Read64,       &
            IOS_Action_Read32,         &
            IOS_Action_Read32_Integer, &
            IOS_Action_Read64_Integer)
          node => make_new_node()
          node%metadata   = metadata
          node%payloadTag = payloadTag
          CALL IOS_Put_Node_In_Queue(node)

          !--------------------------------------------------------------
          ! Process Enquiry
          !--------------------------------------------------------------

        CASE (IOS_Action_Enquire)
          node => make_new_node()
          node%metadata   = metadata
          node%payloadTag = payloadTag
          CALL IOS_Put_Node_In_Queue(node)

          !--------------------------------------------------------------
          ! Process A unit sync OP
          !--------------------------------------------------------------

        CASE (IOS_Action_Sync)
          CALL IOS_Put_Metadata_In_Queue(metadata)

        CASE (IOS_Action_Sync_Barrier)
          CALL IOS_Put_Metadata_In_Queue(metadata)
          IF ( IOS_thread_0_calls_mpi ) THEN
            CALL IOS_WaitForQueueDrain()
            ! Send 1 integer to the client to indicate completion
            CALL MPL_Send(metadata,                     &
                1, MPL_INTEGER, metadata%client,        &
                node%payloadTag, Global_Comm,           &
                ierror)
          END IF

          !--------------------------------------------------------------
          ! Process Open/Close/Release/Open the pipe
          !--------------------------------------------------------------

        CASE (IOS_Action_Open,IOS_Action_Close, &
            IOS_Action_Open_Pipe)
          ! Get extra metadata and put action on queue
          CALL acquire_lock_mpi()
          CALL MPL_Recv(metadata%string, IOS_string_max, MPL_CHARACTER, &
              metadata%client, strTag, Global_Comm,                     &
              status, ierror)
          CALL release_lock_mpi()
          CALL IOS_Put_Metadata_In_Queue(metadata)

        CASE (IOS_Action_Release,IOS_Action_FileOp)
          CALL acquire_lock_mpi()
          CALL MPL_Recv(metadata%string, IOS_string_max, MPL_CHARACTER, &
              metadata%client, strTag, Global_Comm,                     &
              status, ierror)
          CALL release_lock_mpi()

          IF ( IOS_thread_0_calls_mpi ) THEN
            IF ( IOS_Verbosity >= IOS_PrStatus_Oper )     &
                WRITE(6,'(A,A)')'IOS: Info: Listener: ',  &
                'Received Release or FileOp, stalling listener task'
            CALL IOS_WaitForQueueDrain() 
            IF ( IOS_Verbosity >= IOS_PrStatus_Oper ) THEN
              WRITE(6,'(A)')'IOS: Info: Listener: Queue drained'
              WRITE(6,'(A,A,I4)')'IOS: Info: Listener: ', &
                  'Calling SoftSync on server ',          &
                  io_server_for_unit(jc_unit)
            END IF

            CALL acquire_lock_mpi()
            CALL IOS_SoftSync(io_server_for_unit(jc_unit))
            CALL release_lock_mpi()

            IF ( IOS_Verbosity >= IOS_PrStatus_Oper )     &
                WRITE(6,'(A)')'IOS: Info: Listener: Sync Done!'
          END IF
          CALL IOS_Put_metadata_in_queue(metadata)

          !--------------------------------------------------------------
          ! Process Fence requests
          !--------------------------------------------------------------

        CASE (IOS_Action_Fence)
          IF ( IOS_thread_0_calls_mpi ) THEN
            CALL IOS_WaitForQueueDrain()
            CALL acquire_lock_mpi()
            CALL IOS_SoftSync(metadata%address)
            CALL release_lock_mpi()
          END IF
          CALL IOS_Put_metadata_in_queue(metadata)

          !--------------------------------------------------------------
          ! Process Config, Close the Pipe file
          !--------------------------------------------------------------

        CASE (IOS_Action_Config, IOS_Action_Close_Pipe)
          CALL IOS_Put_metadata_in_queue(metadata)

          !--------------------------------------------------------------
          ! Process Flush/Setpos
          !--------------------------------------------------------------

        CASE (IOS_Action_Flush, &
            IOS_Action_Setpos,  &
            IOS_Action_StashSetPos )
          CALL IOS_Put_Metadata_In_Queue(metadata)

          !--------------------------------------------------------------
          ! Process Finish
          !--------------------------------------------------------------

        CASE (IOS_Action_Finish)
          CALL IOS_Put_Metadata_In_Queue(metadata)

          ! Cancel the preemption receive request, so that we don't have any
          ! pending receives at shutdown.
          CALL MPL_Cancel(express_command_request, ierror)

! DEPENDS ON : get_wallclock_time
          t=get_wallclock_time()
          WRITE(6,'(A)') &
              '********************************************************'
          WRITE(6,'(A,F10.3,A)') &
              '* IO SERVER:LISTENER RECEIVED FINISH CMD AT ', &
              t-IOS_Start_Time,' *'
          WRITE(6,'(A)') &
              '********************************************************'

          !--------------------------------------------------------------
          ! Process Start of New Timestep
          !--------------------------------------------------------------

        CASE (IOS_Action_Process)
          timestep=timestep+1
          ! DEPENDS ON : get_wallclock_time
          t=get_wallclock_time()
          IF ( IOS_Verbosity >= IOS_PrStatus_Oper ) &
              WRITE(6,'(A,I4,A,F10.3)')                           &
              'IOS: Info: Listener: Entering Timestep ',timestep, &
              ' at time=',t-IOS_start_time
          CALL IOS_Put_Metadata_In_Queue(metadata)

          !-------------------------------------------------------------------
          ! Process Initialising a fixed Header
          !-------------------------------------------------------------------

        CASE (IOS_Action_StashInitHeader)        
          node => make_new_node()
          node % metadata = metadata
          CALL IOS_Put_Node_In_Queue(node)        

          !-------------------------------------------------------------------
          ! Process Setting the content of a fixed Header
          !-------------------------------------------------------------------

        CASE (IOS_Action_StashSetHeader)

          CALL IOS_WaitForFreeSpace(metadata % data_size* &
              IOS_BytesPerWord64)
          node => make_new_node()
          node % metadata = metadata
          CALL IOS_Increment_Queue_Length(metadata % data_size * &
              IOS_BytesPerInteger)
          ALLOCATE( node % integer_data( metadata % data_size) )

          CALL acquire_lock_mpi()
          CALL MPL_Recv(node % integer_data, metadata % data_size, &
              MPL_INTEGER,metadata%client, PayloadTag, Global_Comm,&
              status, ierror)
          CALL release_lock_mpi()

          CALL IOS_Put_Node_In_Queue(node)

          !-------------------------------------------------------------------
          ! Process Writes of PP lookup data
          !-------------------------------------------------------------------

        CASE (IOS_ACTION_StashWritePPLookup, &
              IOS_ACTION_MergePPLookup)

          CALL IOS_WaitForFreeSpace(metadata % data_size* &
              IOS_BytesPerWord64)
          node => make_new_node()
          node % metadata = metadata
          CALL IOS_Increment_Queue_Length(metadata % data_size * &
              IOS_BytesPerInteger)
          ALLOCATE( node % integer_data( metadata % data_size) )

          CALL acquire_lock_mpi()
          CALL MPL_Recv(node % integer_data, metadata % data_size, &
              MPL_INTEGER, metadata%client, PayloadTag, &
              Global_Comm, status, ierror)
          CALL release_lock_mpi()

          ! Put data on Queue
          CALL IOS_Put_Node_In_Queue(node)

          !-------------------------------------------------------------------
          ! Process The initialisation of PP lookup buffers
          !-------------------------------------------------------------------

        CASE(IOS_Action_StashInitPPLookup)
          CALL IOS_WaitForFreeSpace(metadata % data_size* &
              IOS_BytesPerWord64)
          node => make_new_node()
          node % metadata = metadata
          CALL IOS_Increment_Queue_Length(metadata % data_size * &
              IOS_BytesPerInteger)
          ALLOCATE( node % integer_data( metadata % data_size) )

          CALL acquire_lock_mpi()
          CALL MPL_Recv(node % integer_data, metadata % data_size, & 
              MPL_INTEGER, metadata%client, payloadTag,            &
              Global_Comm,                                         &
              status, ierror)
          CALL release_lock_mpi()

          ! Put data on Queue
          CALL IOS_Put_Node_In_Queue(node)

          !-------------------------------------------------------------------
          ! Process The initialisation of model geometry
          !-------------------------------------------------------------------

        CASE(IOS_Action_StashInitModel)

          ALLOCATE( idata ( metadata % data_size) )
          CALL acquire_lock_mpi()
          CALL MPL_Recv(idata, metadata % data_size, MPL_INTEGER,  &
              metadata%client, payloadTag,                         &
              Global_Comm,                                         &
              status, ierror)
          CALL release_lock_mpi()

          ! Do not put data on Queue, just call the stash server with the details
          ! directly but let us drain the queue first out of sheer paranoia
          ! Just in case the other thread is in the middle of doing something
          ! that we havn't thought carefully about....
          CALL IOS_WaitForQueueDrain()

          CALL acquire_lock_mpi()
          CALL IOS_Server_Geometry_Init(idata)
          CALL release_lock_mpi()
          CALL IOS_Server_Coupler_Init()
          DEALLOCATE(idata)

          !-------------------------------------------------------------------
          ! Process Write PP file data
          !-------------------------------------------------------------------

        CASE (IOS_Action_StashWritePPData,&
            IOS_Action_StashWriteDumpData)
          CALL IOS_WaitForFreeSpace(metadata % data_size* &
              IOS_BytesPerWord64)

          ! we need to allow for the queue to accommodate the receive
          ! buffers even if we do not allocate them here
          ! otherwise we may see a deadlock at the actual allocation time

          node => make_new_node()
          node % metadata = metadata
          ALLOCATE( node % integer_data( metadata % data_size) )
          CALL acquire_lock_mpi()
          CALL MPL_Recv(node % integer_data, metadata % data_size,      &
              MPL_INTEGER, metadata%client, payloadTag, Global_Comm,    &
              status, ierror)
          CALL release_lock_mpi()

          ! Cheap sanity check that the buffer we received looks like a stash
          ! descriptor block.....
          IF (node%integer_data(1) /= IOS_stash_record_start) THEN
            WRITE(IOS_listener_message,*)                               &
                'IOS check failure, 1st word of control buffer is not ',&
                IOS_stash_record_start
            CALL IOS_ereport( RoutineName, 999, IOS_listener_message, &
                md=metadata )
          END IF

          reserveSpace=             &
              getMaxFieldDomain()*  &
              levels_in_pack(node)* &
              procs
          CALL IOS_WaitForFreeSpace(reserveSpace* &
              IOS_BytesPerWord64)
          CALL IOS_Increment_Queue_Length(reserveSpace * &
              IOS_BytesPerReal) 

          CALL IOS_Stash_server_init_recvs(node)           
          CALL IOS_Put_Node_In_Queue(node)
          CALL IOS_Increment_Queue_Length(metadata % data_size * &
              IOS_BytesPerInteger)


          !--------------------------------------------------------------
          ! Process any other unrecognised action 
          !--------------------------------------------------------------

        CASE DEFAULT
          WRITE(IOS_listener_message,'(A,I3,A)') 'Action ', &
              metadata%action,' not recognised!' 
          CALL IOS_ereport( RoutineName, 60, IOS_listener_message, &
              md=metadata )

        END SELECT
      END IF
      NULLIFY(node)
    END DO

    WRITE(6,'(A)') 'IOS: Info: Listener: Exited listen loop'

    !We cannot just return at this point. In MPI models where only 
    !thread zero can call MPI, the writer thread may need to proxy MPI 
    !commands via this thread. We need to hang around and help him 
    !clear the queue. This is the same logic as we used at the top.
    IF ( lhook ) CALL dr_hook &
        ('IOS_LISTENER:TAILSPIN',zhook_in,zhook_handle2)

    IF ( IOS_thread_0_calls_mpi ) THEN
      DO WHILE (IOS_getQueueItems() > 0)
        CALL IOS_Aux_Thread1_Assist()
        CALL um_sleep(IOS_backoff_interval)
      END DO
    END IF

    IF ( lhook ) CALL dr_hook('IOS_LISTENER:TAILSPIN', &
        zhook_out,zhook_handle2)

    ! Final output
    WRITE(6,'(A)') 'IOS: Info: Listener: Process closing'
    IF ( lhook ) CALL dr_hook('IOS_LISTENER',zhook_out,zhook_handle1)

    RETURN
  END SUBROUTINE IOS_Listener

  SUBROUTINE construct_status(unit,status)
    ! Note we use the 1=true convention for the status items
    ! Status only knows about open/closed thus far
    USE IO, ONLY : is_unit_open
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(OUT):: status(1)    
    INTEGER             :: opn

    ! Is the file open?
    status(1)=0
    IF (is_unit_open(unit))status(1)=1 

  END SUBROUTINE construct_status

  SUBROUTINE IOS_Aux_Thread1_Assist()
    USE IOS_Stash_Server,       ONLY : &
        IOS_Stash_Server_Finish_Recvs
    IMPLICIT NONE

    TYPE(IOS_node_type), POINTER :: q_head
    LOGICAL                      :: test_status

    NULLIFY(q_head)

!$OMP FLUSH
    IF (IOS_ReadCompletionRequested==1) THEN
      !Thread 1 wants me to CALL MPI for him
      CALL IOS_Get_Node_From_Queue(q_head)
      CALL IOS_Aux_Read_Assist(q_head)
      NULLIFY (q_head)
      IOS_ReadCompletionRequested=0
!$OMP FLUSH
    ELSE IF (IOS_AsyncCompletionRequested==1) THEN
      !Thread 1 wants me to CALL MPI for him
      CALL IOS_Get_Node_From_Queue(q_head)
      test_status=IOS_Stash_Server_Finish_Recvs(q_head)
      NULLIFY (q_head)
      IF (test_status)IOS_AsyncCompletionRequested=0
!$OMP FLUSH
    ELSE IF (IOS_EnqCompletionRequested==1) THEN
      !Thread 1 wants me to CALL MPI for him
      CALL IOS_Get_Node_From_Queue(q_head)
      CALL IOS_Aux_Enq_Assist(q_head)
      NULLIFY (q_head)
      IOS_EnqCompletionRequested=0
!$OMP FLUSH
    END IF
  END SUBROUTINE IOS_Aux_Thread1_Assist


  SUBROUTINE IOS_Aux_Read_Assist(node)
    IMPLICIT NONE
    TYPE(IOS_node_type),POINTER  :: node
    INTEGER                      :: ierror
    REAL(KIND=jprb)              :: zhook_handle

    IF ( lhook ) CALL dr_hook('IOS_LISTENER:READ_ASSIST', &
        zhook_in,zhook_handle)

    SELECT CASE (node%metadata%action)

    CASE (IOS_Action_Read64)
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
        CALL MPL_Send(node%real_data,                                       &
            node%metadata%data_size*IOS_TUPerWord64,                        &
            IOS_TUType, node%metadata%client,                               &
            node%payloadTag, Global_Comm, ierror)
      END IF
      CALL release_lock_mpi()
    CASE  (IOS_Action_Read32 )
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
        CALL MPL_Send(node%real32_data,                                     &
            node%metadata%data_size*IOS_TUPerWord32,                        &
            IOS_TUType, node%metadata%client,                               &
            node%payloadTag, Global_Comm, ierror)
      END IF
      CALL release_lock_mpi()
    CASE (IOS_Action_Read64_Integer)
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
        CALL MPL_Send(node%integer_data,                                    &
            node%metadata%data_size*IOS_TUPerWord64,                        &
            IOS_TUType, node%metadata%client,                               &
            node%payloadTag, Global_Comm, ierror)
      END IF
      CALL release_lock_mpi()
    CASE DEFAULT
      CALL IOS_Ereport('IOS_Aux_Read_Assist', 10, &
          'Received a non-read operation',        &
           md=node%metadata)

    END SELECT

    IF ( lhook ) CALL dr_hook('IOS_LISTENER:READ_ASSIST', &
        zhook_out,zhook_handle)

  END SUBROUTINE IOS_Aux_Read_Assist

  SUBROUTINE IOS_Aux_Enq_Assist(node)
    IMPLICIT NONE
    TYPE(IOS_node_type),POINTER  :: node
    INTEGER                      :: ierror
    REAL(KIND=jprb)              :: zhook_handle
    
    IF ( lhook ) CALL dr_hook('IOS_LISTENER:ENQ_ASSIST', &
        zhook_in,zhook_handle)

    CALL acquire_lock_mpi()
    CALL MPL_Send(IOS_Unit_status,            &
        IOS_State_Size, MPL_INTEGER,          &
        node%metadata%client, node%payloadTag,&
        Global_Comm, ierror)
    CALL release_lock_mpi()

    IF ( lhook ) CALL dr_hook('IOS_LISTENER:ENQ_ASSIST', &
        zhook_out,zhook_handle)

  END SUBROUTINE IOS_Aux_Enq_Assist

  SUBROUTINE IOS_WaitForQueueDrain()
    IMPLICIT NONE
    LOGICAL                     :: done
    REAL(KIND=jprb)             :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_WATFORQUEUEDRAIN',zhook_in,zhook_handle)

    done=.FALSE.
    DO WHILE (.NOT.done)
      IF (.NOT.waitforDrain()) THEN
        !it returned before space was available so offer help
        !Check to see if thread 1 needs attention 
        CALL IOS_Aux_Thread1_Assist()   
      ELSE
        done=.TRUE.
      END IF
    END DO

    IF (lhook) CALL dr_hook('IOS_WATFORQUEUEDRAIN',zhook_out,zhook_handle)

  END SUBROUTINE IOS_WaitForQueueDrain

  SUBROUTINE IOS_WaitForFreeSpace(size)
    IMPLICIT NONE
    INTEGER, INTENT(IN)         :: size
    LOGICAL                     :: done
    REAL(KIND=jprb)             :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_WAITFORFREESPACE',zhook_in,zhook_handle)

    done=.FALSE.
    DO WHILE (.NOT.done)
      IF (.NOT.HasFreeSpace(size)) THEN
        !it returned before space was available so offer help
        !Check to see if thread 1 needs attention  
        CALL IOS_Aux_Thread1_Assist()   
      ELSE
        done=.TRUE.
      END IF
    END DO

    IF (lhook) CALL dr_hook('IOS_WAITFORFREESPACE',zhook_out,zhook_handle)

  END SUBROUTINE IOS_WaitForFreeSpace
 
END MODULE IO_Server_Listener
