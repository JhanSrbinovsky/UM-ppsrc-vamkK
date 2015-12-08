! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

! Implements a ring buffer for client side IOS operations
! and manages dispatch and data receipt into elements of that queue

MODULE IOS_Client_Queue

  USE IOS_Common
  USE yomhook, ONLY: lhook, dr_hook  
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE
  
  CHARACTER (LEN=*), &
      PARAMETER, PRIVATE          :: RoutineName = 'IOS_Client_Q'
  INTEGER, POINTER                :: ioGlobalNo(:,:)=>NULL()
  INTEGER                         :: IOS_Client_Queue_Slot
  REAL                            :: IOS_cl_queue_stall=0.0



  TYPE(IOS_async_object), POINTER :: IOS_Dispatch_Queue(:)

  CHARACTER (LEN=IOS_err_strlen), &
      PRIVATE                     :: IOS_clq_message  

  ! params/vars  for dr_hook
  INTEGER(KIND=jpim), &
      PARAMETER, PRIVATE          :: zhook_in  = 0
  INTEGER(KIND=jpim), &
      PARAMETER, PRIVATE          :: zhook_out = 1
  REAL, EXTERNAL                  :: get_wallclock_time
  
CONTAINS

  
  SUBROUTINE IOS_Client_Init()
    USE MPL, ONLY : MPL_REQUEST_NULL
    IMPLICIT NONE
    INTEGER :: i
    REAL(KIND=jprb)     :: zhook_handle
    
    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:CLIENT_INIT', &
        zhook_in,zhook_handle)


    ALLOCATE(IOS_Dispatch_Queue(IOS_Concurrency)) 
    ! Never deallocated.

    ALLOCATE(IOS_Sequence_id(IOS_Server_Groups,IOS_Tasks_per_server))  
    ! Never deallocated.
    ALLOCATE(ioGlobalNo     (IOS_Server_Groups,IOS_Tasks_per_server))
    ! Never deallocated.
    IOS_Sequence_id(:,:)=0
    IOS_Client_Queue_Slot=0
    DO i=1,IOS_Concurrency
      IOS_Dispatch_Queue(i)%state       =IOS_QUEUE_SLOT_UNUSED
      IOS_Dispatch_Queue(i)%request(:)  =MPL_REQUEST_NULL
    END DO

    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:CLIENT_INIT', &
        zhook_out,zhook_handle)
    
  END SUBROUTINE IOS_Client_Init


  SUBROUTINE IOS_assign_server_for_unit(unit)
    USE MPL, ONLY : MPL_INTEGER
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER, PARAMETER  :: loadBalalanceBcastRoot  = 0
    INTEGER, PARAMETER  :: loadBalalanceBcastItems = 1
    INTEGER             :: pe
    INTEGER             :: subtask
    INTEGER             :: qHandle
    INTEGER,SAVE        :: nextserver=0
    INTEGER             :: targetObject
    INTEGER             :: loads(loadBalanceDataSize)
    INTEGER             :: minLoad
    INTEGER             :: minSize
    INTEGER             :: curLoad
    INTEGER             :: curSize
    INTEGER             :: minLoad_server
    INTEGER             :: minSize_server
    INTEGER             :: errorCode
    INTEGER(IOS_TUKind), &
        POINTER         :: recvBuffer(:)
    LOGICAL             :: reallocation

    TYPE(IOS_metadata_type), &
        POINTER         :: metadata

    ErrorCode=-10
    targetObject=-1
    reallocation=.FALSE.

    IF (L_IOS_Active()) THEN

      SELECT CASE(IOS_Unit_Alloc_Policy)
      CASE(IOS_Unit_Alloc_Static)
        
        ! assignment was made at init time, 
        ! so this is a null op or an error

        IF (io_server_for_unit(unit) == IOS_No_Server )               &
            CALL IOS_Ereport('IOS_Assign_server_for_unit', ErrorCode, &
            'Unallocated unit was presented in a static scheme' )    

      CASE(IOS_Unit_Alloc_AtFirstUse)

        ! assignment was not made at init time, 
        ! but we may not reassign existing mappings

        IF (io_server_for_unit(unit) == IOS_No_Server ) THEN

          nextserver=nextserver+1
          IF (nextserver > IOS_Server_Groups) nextserver=1
          targetObject=io_servers(nextServer,1)

        END IF

      CASE(IOS_Unit_Alloc_Dynamic_Rotate)

        nextserver=nextserver+1
        IF (nextserver > IOS_Server_Groups) nextserver=1
        targetObject=io_servers(nextServer,1)

        IF (io_server_for_unit(unit) /= IOS_No_Server ) &
            reallocation=.TRUE.
        
      CASE(IOS_Unit_Alloc_Dynamic_LB)
        
        IF ( io_server_for_unit(unit) == IOS_No_Server ) THEN

          nextserver=nextserver+1
          IF (nextserver > IOS_Server_Groups) nextserver=1
          targetObject=io_servers(nextServer,1)

        ELSE
          
          IF (model_rank == 0 ) THEN
            minLoad=HUGE(minload)
            minSize=HUGE(minSize)
            minLoad_server=-1
            minSize_server=-1
            
            DO pe=1,IOS_Server_Groups
              qHandle=IOS_Init_MD(-1*io_servers(pe,1),    &
                  IOS_No_Location, IOS_Action_LoadStatus, &
                  datasize=loadBalanceDataSize)        
              CALL IOS_Send(qHandle) 
              recvBuffer => IOS_attach_recvBuffer(qHandle)
              
              loads(:)=-99
              CALL UM_MemCpy64(loads,recvBuffer,loadBalanceDataSize)
              
              IF (io_server_for_unit(unit)==io_servers(pe,1)) THEN
                IF (IOS_Verbosity>=IOS_prstatus_Oper) THEN
                  WRITE(6,'(A,F10.3,A,I4,A,I4)')                 &
                      'IOS: Info: Current queue size=',          &
                      loads(loadBalanceQueueData)/1024.0/1024.0, &
                      'MB items=',                               &
                      loads(loadBalanceQueueLen),                &
                      ' on server ',                             &
                      io_servers(pe,1)
                END IF
              END IF
                  
              IF (loads(loadBalanceQueueLen) < minLoad) THEN
                minload        = loads(loadBalanceQueueLen) 
                minLoad_server = pe
              END IF
              IF (loads(loadBalanceQueueData) < minSize) THEN
                minSize        = loads(loadBalanceQueueData) 
                minSize_server = pe
              END IF
            END DO
            
            IF (IOS_Verbosity>=IOS_prstatus_Oper) THEN
              WRITE(6,'(A,F8.2,A,I4,A,I4,A,I4,A)') &
                  'IOS: Info: Min queue size=',    &
                  minSize/1024.0/1024.0,           &
                  ' MB (',                         &
                  io_servers(minSize_Server,1),    &
                  ') Min queue load=',             &
                  minLoad,                         &
                  ' items (',                      &
                  io_servers(minLoad_Server,1),')'
            END IF
            
          END IF

          IF (Model_Procs > 1) THEN
            CALL MPL_Bcast(minSize_Server,            &
                loadBalalanceBcastItems, MPL_INTEGER, &
                loadBalalanceBcastRoot, Model_Comm,   &
                errorCode)
          END IF
          targetObject=io_servers(minSize_server,1)
          
          IF (io_server_for_unit(unit) == targetObject) THEN
            
            ! Do not waste the servers energy, if we are remapping 
            ! to the same place

            IF (IOS_Verbosity>=IOS_prstatus_Oper) THEN
              WRITE(6,'(A,I3,A,I4)')                       &
                  'IOS: Info: Unit ',unit,                 &
                  ' is not reassigned, and remains on ',   &
                  io_server_for_unit(unit)
            END IF
            reallocation=.FALSE.
            targetObject=-1
          ELSE
            reallocation=.TRUE.
          END IF

        END IF
        
      CASE DEFAULT

        ErrorCode=-10
        CALL IOS_Ereport('IOS_Assign_server_for_unit', ErrorCode, &
            'Unknown server allocation policy' )      

      END SELECT
      
      
      IF(targetObject >= 0) THEN
        io_server_for_unit_lookup(unit)=targetObject
        
        IF (IOS_Verbosity>=IOS_prstatus_Oper) THEN
          IF (reallocation) THEN
            WRITE(6,'(A,I3,A,I4)')                       &
                'IOS: Info: Unit ',unit,                 &
                ' is reassigned to IO server ',          &
                targetObject
          ELSE
            WRITE(6,'(A,I3,A,I4)')                       &
                'IOS: Info: Unit ',unit,                 &
                ' is assigned to IO server ',            &
                targetObject
          END IF
        END IF
        
        IF (model_rank == 0 ) THEN
          DO pe=1,IOS_Server_Groups
             DO subtask=1,IOS_tasks_per_server
                qHandle=IOS_Init_MD(-1*io_servers(pe,subtask), &
                    IOS_No_Location, IOS_Action_Assign_Unit)        
                metadata => IOS_Attach_Metadata(qHandle)
                metadata % intent  = unit
                metadata % address = targetObject
                CALL IOS_Send(qHandle)
             END DO
          END DO
        END IF
      END IF

    END IF

  END SUBROUTINE IOS_assign_server_for_unit


  FUNCTION IOS_getDestPe(handle) RESULT(pe)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: handle
    INTEGER             :: pe
    pe=IOS_Dispatch_Queue(handle)%pe
  END FUNCTION IOS_getDestPe


  INTEGER FUNCTION IOS_Init_MD                            &
      (targetObject, location, action, dataSize, bcast,   &
      rank) RESULT(handle)
    IMPLICIT NONE
    INTEGER, INTENT(IN)                  :: targetObject
    INTEGER, INTENT(IN)                  :: action
    INTEGER, INTENT(IN)                  :: location
    INTEGER, OPTIONAL                    :: dataSize
    INTEGER, OPTIONAL                    :: bcast
    INTEGER, OPTIONAL                    :: rank
    INTEGER                              :: pe
    INTEGER                              :: nextSlot
    LOGICAL                              :: OK
    LOGICAL                              :: flag
    LOGICAL                              :: l_bcast
    REAL                                 :: t1,t2
    REAL, SAVE                           :: t_stall=0.0
    INTEGER                              :: errorFlag
    REAL(KIND=jprb)                      :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:INIT_MD', &
        zhook_in,zhook_handle)

    OK=assert_client()

    IF (IOS_Verbosity>=IOS_prstatus_debug)             &    
        WRITE(6,'(A,A)')                               &
        'IOS: Info: Getting protocol queue slot for ', &
        TRIM(IOS_ActionName(action))

    NextSlot=IOS_Client_Queue_Slot+1
    IF (NextSlot>IOS_Concurrency)NextSlot=1
    
    IF (IOS_Dispatch_Queue(NextSlot)%state== &
        IOS_QUEUE_SLOT_INITIALIZED) THEN

      WRITE(IOS_clq_message,*) &
          'The next queue slot is already initialzed'
      errorFlag=99
      CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )      

    ELSE IF(IOS_Dispatch_Queue(NextSlot)%state== &
        IOS_QUEUE_SLOT_PARTFILLED) THEN

      WRITE(IOS_clq_message,*)'The next queue slot already filled '
      errorFlag=99
      CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )      

    ELSE IF(IOS_Dispatch_Queue(NextSlot)%state== &
        IOS_QUEUE_SLOT_DISPATCHED) THEN

      CALL IOS_WaitBuffer(NextSlot,'IOS_init_md :wating for buffer ')

    ELSE IF (IOS_Dispatch_Queue(NextSlot)%state/=IOS_QUEUE_SLOT_UNUSED) THEN

      WRITE(IOS_clq_message,'(A,I6)')'IOS_Client_Queue: UNKNOWN STATE=',&
          IOS_Dispatch_Queue(NextSlot)%state
      errorFlag=99
      CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )      

    END IF

    ! We can now safely reuse our target slot
    IOS_Client_Queue_Slot=nextSlot

    ! The handle for this request is the slot id
    handle=IOS_Client_Queue_Slot

    ! Set defaults from arguments
    IOS_Dispatch_Queue(handle)%md%action  = action
    IOS_Dispatch_Queue(handle)%md%unit    = targetObject
    IOS_Dispatch_Queue(handle)%md%address = location
    
    ! Blank other fields !
    IOS_Dispatch_Queue(handle)%md%name_length = 0
    IOS_Dispatch_Queue(handle)%md%intent      = 0
    IOS_Dispatch_Queue(handle)%md%delete      = 0
    IOS_Dispatch_Queue(handle)%md%data_size   = 0
    IOS_Dispatch_Queue(handle)%md%client      = model_rank
    IOS_Dispatch_Queue(handle)%md%handle      = 0
    IOS_Dispatch_Queue(handle)%md%Originating_Slot = -99

    NULLIFY(IOS_Dispatch_Queue(handle)%payload)
    WRITE(IOS_Dispatch_Queue(handle)%md%string,*)'NotSet'

    ! Figure out the destination pe
    IF (targetObject <= 0) THEN
      IOS_Dispatch_Queue(handle)%pe=-1*targetObject
    ELSE IF (targetObject > 0 .AND. targetObject < IOS_unit_first) THEN
      WRITE(IOS_clq_message,'(A,I5,A)') 'A target object of ', &
          targetObject,' is not allowed' 
      errorFlag=61
      CALL IOS_ereport( RoutineName, errorFlag, IOS_clq_message )
    ELSE      
      IOS_Dispatch_Queue(handle)%pe=io_server_for_unit(targetObject)
      IF (IOS_Dispatch_Queue(handle)%pe == IOS_No_Server) THEN
        WRITE(IOS_clq_message,'(A,I3,A)') 'No server for unit ', &
            targetObject,' has been assigned' 
        errorFlag=61
        CALL IOS_Ereport( RoutineName, errorFlag, IOS_clq_message )
      END IF

      IF (PRESENT(rank)) THEN
        IOS_Dispatch_Queue(handle)%pe= &
            io_servers(ioServerNo(IOS_Dispatch_Queue(handle)%pe),rank)
      END IF
    END IF

    ! Do we need to set the broadcast flag?
    l_bcast=.FALSE.
    IF (PRESENT(bcast)) THEN
      IF (bcast==IOS_Read_Broadcast) THEN
        l_bcast=.TRUE.
        IOS_Dispatch_Queue(handle)%md%intent=bcast
      END IF
    END IF

    ! Do we need a data buffer?
    IF (PRESENT(dataSize)) THEN
      IOS_Dispatch_Queue(handle)%md%data_size=dataSize

      IF (.NOT.l_bcast) THEN

        SELECT CASE(action)
        CASE (IOS_Action_Write32,      &
            IOS_Action_Read32_Integer, &
            IOS_Action_Read32          &
            )
          ALLOCATE (IOS_Dispatch_Queue(handle)%payload &
              (dataSize*IOS_TUPerWord32))
        CASE (                            &
            IOS_Action_Enquire,           &
            IOS_Action_LoadStatus,        &
            IOS_Action_Write64,           &
            IOS_Action_Read64,            &
            IOS_Action_Read64_Integer,    &
            IOS_Action_Sync_Barrier,      &
            IOS_Action_StashInitModel,    &
            IOS_Action_StashInitPPLookup, &
            IOS_Action_Write_PP_Prepacked,&
            IOS_Action_StashWritePPData,  & 
            IOS_Action_StashWritePPLookup,&
            IOS_Action_StashInitHeader,   &
            IOS_Action_StashSetHeader,    &
            IOS_Action_StashWriteDumpData,&
            IOS_Action_DumpInitModel     ,&
            IOS_Action_MergePPLookup)
          ALLOCATE (IOS_Dispatch_Queue(handle)%payload &
              (dataSize*IOS_TUPerWord64))        
        CASE Default
          WRITE(IOS_clq_message,'(A,I3,A)') 'Action ', &
              action,' not recognised!' 
          errorFlag=60
          CALL IOS_Ereport( RoutineName, errorFlag, IOS_clq_message )        
        END SELECT
      END IF
    END IF

    IF (IOS_Verbosity>=IOS_prstatus_debug)                & 
        WRITE(6,'(A,I4)')'IOS: Info: New queue slot @ ',handle

    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:INIT_MD', &
        zhook_out,zhook_handle)

  END FUNCTION IOS_Init_MD


  FUNCTION IOS_Attach_Metadata(handle) RESULT(md)
    IMPLICIT NONE
    INTEGER, INTENT(IN)              :: handle
    TYPE(IOS_Metadata_Type), POINTER :: md
    REAL(KIND=jprb)                  :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:ATTACH_MD', &
        zhook_in,zhook_handle)

    md => IOS_Dispatch_Queue(handle)%md

    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:ATTACH_MD', &
        zhook_out,zhook_handle)

  END FUNCTION IOS_Attach_Metadata


  FUNCTION IOS_Attach_SendBuffer(handle) RESULT(sb)
    IMPLICIT NONE
    INTEGER, INTENT(IN)               :: handle
    INTEGER(KIND=IOS_TUKind), POINTER :: sb(:)
    REAL(KIND=jprb)                   :: zhook_handle
    INTEGER                           :: errorFlag

    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:ATTACH_SB', &
        zhook_in,zhook_handle)

    errorFlag=99

    IF (.NOT.ASSOCIATED(IOS_Dispatch_Queue(handle)%payload))THEN
      WRITE(IOS_clq_message,*) &
          'Attach_SendBuffer Failed, no buffer allocated'
      CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )      
    END IF

    sb => IOS_Dispatch_Queue(handle)%payload(:)

    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:ATTACH_SB', &
        zhook_out,zhook_handle)

  END FUNCTION IOS_Attach_SendBuffer


  FUNCTION IOS_Attach_RecvBuffer(handle) RESULT(rb)
    IMPLICIT NONE
    INTEGER, INTENT(IN)               :: handle
    INTEGER(KIND=IOS_TUKind), POINTER :: rb(:)
    REAL(KIND=jprb)                   :: zhook_handle
    INTEGER                           :: errorFlag

    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:ATTACH_RB', &
        zhook_in,zhook_handle)
    
    IF (.NOT.ASSOCIATED(IOS_Dispatch_Queue(handle)%payload))THEN
      WRITE(IOS_clq_message,*) &
          'Attach_RecvBuffer Failed, no buffer allocated'
      errorFlag=99
      CALL IOS_Ereport( RoutineName,errorFlag , IOS_clq_message )      
    END IF

    rb => IOS_Dispatch_Queue(handle)%payload

    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:ATTACH_RB', &
        zhook_out,zhook_handle)

  END FUNCTION IOS_Attach_RecvBuffer


  SUBROUTINE IOS_Send(handle,hasString)
    USE MPL, ONLY :      &
        MPL_STATUS_SIZE, &
        MPL_INTEGER,     &
        MPL_CHARACTER
    
    IMPLICIT NONE
    INTEGER, INTENT(IN)   :: handle
    LOGICAL, OPTIONAL     :: hasString
    INTEGER               :: errorFlag
    INTEGER               :: pe
    INTEGER               :: ioServer
    INTEGER               :: ioServerR
    INTEGER               :: tagOffset
    INTEGER               :: mdtag
    INTEGER               :: pltag
    INTEGER               :: strtag
    INTEGER               :: status(MPL_STATUS_SIZE)
    REAL(KIND=jprb)       :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:SEND', &
        zhook_in,zhook_handle)

    pe        = ABS(IOS_Dispatch_Queue(handle)%pe)
    ioServer  = ioServerNo(pe)
    ioServerR = ioServerRank(pe)+1

    IF ( IOS_Dispatch_Queue(handle)%md%action == &
        IOS_Action_LoadStatus) THEN
      mdtag  = IOS_Request_Tag_Express
      pltag  = IOS_Request_Tag_Express
      strtag = IOS_Request_Tag_Express
    ELSE

      IOS_Sequence_ID(ioServer,ioServerR)= &
          IOS_Sequence_ID(ioServer,ioServerR)+1
      tagOffset=MOD(IOS_Sequence_ID(ioServer,ioServerR), &
          IOS_Request_Tag_Gap)

      mdtag  = IOS_Request_Tag_Base+tagOffset
      pltag  = IOS_Request_Tag_Payload+tagOffset
      strtag = IOS_Request_Tag_Str+tagOffset
    END IF

! Leave in handy for debugging.
!    WRITE(6,*)'IOS SEND: -------------------------------------------'
!    WRITE(6,*)'IOS SEND: targ pe=',pe,'  seq=',IOS_Sequence_ID(ioServer,ioServerR)
!    WRITE(6,*)'IOS SEND: server=',ioServer
!    WRITE(6,*)'IOS SEND: rank=',ioServerR
!    WRITE(6,*)'IOS SEND: handle=',handle
!    WRITE(6,*)'IOS SEND: tagOffset=',tagOffset,mdtag,pltag,strtag
!    WRITE(6,*)'IOS SEND: action=',IOS_Dispatch_Queue(handle)%md%action
!    WRITE(6,*)'IOS SEND: action=',TRIM(IOS_ActionName(IOS_Dispatch_Queue(handle)%md%action))
!    WRITE(6,*)'IOS SEND: intent=',IOS_Dispatch_Queue(handle)%md%intent


    CALL MPL_iSend(                            &
        IOS_Dispatch_Queue(handle)%md,         &
        IOS_MD_Len,                            &
        MPL_INTEGER,                           &
        IOS_Dispatch_Queue(handle)%pe,         &
        mdtag,                                 &
        Global_Comm,                           &
        IOS_Dispatch_Queue(handle)%request(1), &
        errorFlag)

    IF (PRESENT(hasString)) THEN
      IF (hasString) THEN
        CALL MPL_iSend(                            &
            IOS_Dispatch_Queue(handle)%md%string,  &
            IOS_String_Max,                        &
            MPL_CHARACTER,                         &
            IOS_Dispatch_Queue(handle)%pe,         &
            strtag,                                &
            Global_Comm,                           &
            IOS_Dispatch_Queue(handle)%request(2), &
            errorFlag)
      END IF
    END IF

    IF (ASSOCIATED(IOS_Dispatch_Queue(handle)%payload)) THEN
      SELECT CASE(IOS_Dispatch_Queue(handle)%md%action)
      CASE (                         &
          IOS_Action_Enquire,        &
          IOS_Action_LoadStatus,     &
          IOS_Action_Read32_Integer, &
          IOS_Action_Read32,         &
          IOS_Action_Read64_Integer, &
          IOS_Action_Read64          &
          )
        IF (IOS_Dispatch_Queue(handle)%md%intent/=IOS_Read_Broadcast) THEN
          CALL MPL_Recv(                                &            
              IOS_Dispatch_Queue(handle)%payload,       &
              SIZE(IOS_Dispatch_Queue(handle)%payload), &
              IOS_TUType,                               &
              IOS_Dispatch_Queue(handle)%pe,            &
              pltag,                                    &
              Global_Comm,                              &
              status,                                   &
              errorFlag)
        END IF
      CASE DEFAULT      
        CALL MPL_iSend(                               &
            IOS_Dispatch_Queue(handle)%payload,       &
            SIZE(IOS_Dispatch_Queue(handle)%payload), &
            IOS_TUType,                               &
            IOS_Dispatch_Queue(handle)%pe,            &
            pltag,                                    &
            Global_Comm,                              &
            IOS_Dispatch_Queue(handle)%request(3),    &
            errorFlag)
      END SELECT
    END IF
    IOS_Dispatch_Queue(handle)%state=IOS_QUEUE_SLOT_DISPATCHED

    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:SEND', &
        zhook_out,zhook_handle)

  END SUBROUTINE IOS_Send

  !---------------------------------------------------------------------
  ! Finalisation - terminate remote IO Servers
  !---------------------------------------------------------------------
  SUBROUTINE IOS_Finalise()
    IMPLICIT NONE
    INTEGER                            :: server
    INTEGER                            :: subtask
    TYPE(IOS_metadata_type),POINTER    :: metadata
    INTEGER                            :: qHandle

    IF (model_rank == 0 ) THEN
      DO server=1,IOS_Server_Groups
         DO subtask=1,IOS_tasks_per_server
            qHandle=IOS_Init_MD(-1*io_servers(server,subtask), &
                -1,IOS_Action_Finish)        
            CALL IOS_Send(qHandle)
         END DO
      END DO
    END IF

    RETURN
  END SUBROUTINE IOS_Finalise

  !---------------------------------------------------------------------
  ! Shutdown, Close down servers, and tidy up pending messages
  !---------------------------------------------------------------------
  SUBROUTINE IOS_Shutdown()

    IMPLICIT NONE
    INTEGER                              :: LastUsedSlot
    INTEGER                              :: NextSlot
    INTEGER                              :: errorFlag
    LOGICAL                              :: flag
    REAL                                 :: t1
    REAL                                 :: t2
    REAL(KIND=jprb)                      :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:SHUTDOWN', &
        zhook_in, zhook_handle)

    CALL IOS_Finalise()

    ! Wait completion of queue
    lastUsedSlot=IOS_Client_Queue_Slot
    NextSlot    =-1

    IF (NextSlot>IOS_Concurrency)NextSlot=1

! DEPENDS ON : get_wallclock_time
      t1=get_wallclock_time()
    
    DO WHILE (NextSlot /= LastUsedSlot)
      IF (NextSlot == -1) THEN
        NextSlot=lastUsedSlot+1
      ELSE
        NextSlot=NextSlot+1
      END IF

      IF (NextSlot>IOS_Concurrency)NextSlot=1

      IF (IOS_Verbosity>=IOS_prstatus_diag) &
        WRITE(6,'(A,I3)')'IOS: Shutdown: Waiting completion of Q slot: ', &
            NextSlot

      SELECT CASE(IOS_Dispatch_Queue(NextSlot)%state)

      CASE(IOS_QUEUE_SLOT_INITIALIZED)
        IF (IOS_Verbosity>=IOS_prstatus_diag) &
            WRITE(6,'(A)') 'IOS: Shutdown:  --> IOS_QUEUE_SLOT_INITIALIZED'
          
      CASE(IOS_QUEUE_SLOT_PARTFILLED)
        IF (IOS_Verbosity>=IOS_prstatus_diag) &
            WRITE(6,'(A)') 'IOS: Shutdown:  --> IOS_QUEUE_SLOT_PARTFILLED'

      CASE(IOS_QUEUE_SLOT_DISPATCHED)
        IF (IOS_Verbosity>=IOS_prstatus_diag) &
            WRITE(6,'(A)') 'IOS: Shutdown:  --> IOS_QUEUE_SLOT_DISPATCHED'
        CALL IOS_WaitBuffer(NextSlot,'IOS: Shutdown:  --> Waiting for slot ')
        
      CASE(IOS_QUEUE_SLOT_UNUSED)
        IF (IOS_Verbosity>=IOS_prstatus_diag)&
            WRITE(6,'(A)')'IOS: Shutdown:  --> IOS_QUEUE_SLOT_UNUSED'
        
      CASE DEFAULT
        WRITE(IOS_clq_message,*)'IOS_Client_Queue: UNKNOWN STATE=',&
            IOS_Dispatch_Queue(NextSlot)%state
        errorFlag=99
        CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )      
        
      END SELECT
        
    END DO
! DEPENDS ON : get_wallclock_time
    t2=get_wallclock_time()
      
    IF (IOS_Verbosity>=IOS_PrStatus_Normal) THEN
      WRITE(6,'(A,F8.3)')                                     &
          'IOS: Info: total stall time in protocol queue = ', &
          IOS_cl_queue_stall
      WRITE(6,'(A,F8.3)')                                     &
          'IOS: Info: total time waiting for shutdown    = ', &
          (t2-t1)
    END IF

    IF (lhook) CALL dr_hook('IOS_CLIENT_QUEUE:SHUTDOWN', &
        zhook_out,zhook_handle)
    
  END SUBROUTINE IOS_Shutdown

  SUBROUTINE IOS_WaitBuffer(slot,messg)
    USE MPL, ONLY :       &
        MPL_REQUEST_NULL, &
        MPL_STATUS_SIZE
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: slot
    CHARACTER(LEN=*)       :: messg
    INTEGER             :: errorFlag
    INTEGER             :: stat(MPL_STATUS_SIZE,3)
    REAL                :: t1
    REAL                :: t2
    LOGICAL             :: flag
    
    flag=.FALSE.
    CALL MPL_TESTALL &
        (3,IOS_Dispatch_Queue(Slot)%request,flag,stat,errorFlag)
    IF (.NOT.flag) THEN
      IF (IOS_Verbosity>=IOS_prstatus_oper) &
          WRITE(6,'(A,I3)')TRIM(messg),slot
! DEPENDS ON : get_wallclock_time
      t1=get_wallclock_time()
      DO WHILE (.NOT.flag)
        CALL MPL_Testall(3,IOS_Dispatch_Queue(Slot)%request, &
            flag,stat,errorFlag)
! DEPENDS ON : get_wallclock_time
        t2=get_wallclock_time()
        IF (t2-t1>IOS_Timeout) THEN
          WRITE(IOS_clq_message,'(A)') &
              'Timed out waiting for an IOS protocol buffer'
          errorFlag=99              
          CALL IOS_Ereport( RoutineName, errorFlag , IOS_clq_message )      
        END IF
      END DO
      IOS_cl_queue_stall=IOS_cl_queue_stall+(t2-t1)
      IF (IOS_Verbosity>=IOS_prstatus_oper)                   &
          WRITE(6,'(A,F10.3,A,F8.3)')                         &
          'IOS: Info: Stall getting protocol queue slot at ', &
          t1-IOS_Start_time,' of ',t2-t1
    END IF
    IF (ASSOCIATED(IOS_Dispatch_Queue(Slot)%payload)) &
        DEALLOCATE(IOS_Dispatch_Queue(Slot)%payload)
    IOS_Dispatch_Queue(Slot)%request(1:3)=MPL_REQUEST_NULL
    
  END SUBROUTINE IOS_WaitBuffer

END MODULE IOS_Client_Queue
