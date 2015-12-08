! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module defining interfaces for clients of an IO server
! usually called from io.F90 with exceptions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE IOS

  USE IOS_Common   
  USE IOS_Client_Queue
  USE MPL,ONLY :     &
      MPL_Integer,   &
      MPL_Integer8,  &
      MPL_Logical4,  &
      MPL_Logical8,  &
      MPL_REAL4,     &
      MPL_REAL8,     &
      MPL_REAL,      &
      MPL_CHARACTER, &
      MPL_STATUS_SIZE
  USE yomhook, ONLY: lhook, dr_hook  
  USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
  CHARACTER(LEN=80)  :: ios_message

  INTERFACE IOS_Write32
    MODULE PROCEDURE        &
        IOS_Write32_r,      &
        IOS_Write32_i,      &
        IOS_Write32_i2D
  END INTERFACE

  INTERFACE IOS_Write64
    MODULE PROCEDURE        &
        IOS_Write64_r,      &
        IOS_Write64_i,      &
        IOS_Write64_l,      &
        IOS_Write64_r2D,    &
        IOS_Write64_i2D,    &
        IOS_Write64_l2D
  END INTERFACE

  INTERFACE IOS_Read32
    MODULE PROCEDURE        &
        IOS_Read32_r,       &
        IOS_Read32_i,       &
        IOS_Read32_i2D,     &
        IOS_Read32_l
  END INTERFACE

  INTERFACE IOS_Read64
    MODULE PROCEDURE        &
        IOS_Read64_r,       &
        IOS_Read64_i,       &
        IOS_Read64_l,       &
        IOS_Read64_r2D,     &
        IOS_Read64_i2D,     &
        IOS_Read64_l2D
  END INTERFACE

! params/vars  for dr_hook
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

  INTEGER, PRIVATE :: ierror ! generic error code used everywhere

CONTAINS

  SUBROUTINE IOS_fileState(unit,location,state)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit    
    INTEGER, INTENT(IN)        :: location    
    INTEGER                    :: qHandle
    TYPE(IOS_Status)           :: state
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: recvBuffer(:)
    REAL(KIND=jprb)            :: zhook_handle

    IF ( lhook ) CALL dr_hook('IOS:FILESTATE',zhook_in,zhook_handle)

    qHandle=IOS_Init_MD(unit,location,IOS_Action_Enquire, &
        datasize=IOS_State_Size)
    CALL IOS_Send(qHandle)
    recvBuffer => IOS_attach_recvBuffer(qHandle)
    CALL um_memcpy64(state,recvBuffer,IOS_State_Size)
    
    IF ( lhook ) CALL dr_hook('IOS:FILESTATE',zhook_out,zhook_handle)
    
    RETURN
  END SUBROUTINE IOS_fileState


  LOGICAL FUNCTION IOS_Is_Unit_Open(unit)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit    
    TYPE(IOS_Status)           :: state
    REAL(KIND=jprb)            :: zhook_handle

    IF ( lhook ) CALL dr_hook('IOS:IS_UNIT_OPEN',zhook_in,zhook_handle)

    CALL IOS_fileState(unit,IOS_No_Location,state)
    IOS_Is_unit_open=state%isOpen

    IF ( lhook ) CALL dr_hook('IOS:IS_UNIT_OPEN',zhook_out,zhook_handle)

    RETURN
  END FUNCTION IOS_Is_Unit_Open

 
  SUBROUTINE IOS_Close(unit,name,delete)

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: unit   ! unit number
    INTEGER, INTENT(IN)          :: delete ! =0 read only, otherwise r/w
    CHARACTER (LEN=*), INTENT(IN):: name    ! nameironment var/filename
    INTEGER                      :: errorCode=99
    INTEGER                      :: qHandle
    CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_Close'
    TYPE(IOS_metadata_type), &
        POINTER                  :: metadata
    REAL(KIND=jprb):: zhook_handle

    IF ( lhook ) CALL dr_hook('IOS:CLOSE',zhook_in,zhook_handle)

    ! Check name is small enough to tranmit.
    IF ( LEN_TRIM(name) > IOS_String_Max ) THEN
      Ios_message    = 'IOS: String is too big to transmit'
      CALL IOS_Ereport(RoutineName, ErrorCode, Ios_message)
    END IF

    qHandle=IOS_Init_MD(unit,-1,IOS_Action_Close)
    metadata => IOS_Attach_Metadata(qHandle)
    metadata % name_length              = LEN_TRIM(name)
    metadata % delete                   = delete
    metadata % string(1:metadata % name_length) &
                                        = name(1:metadata % name_length)

    ! Determine the IO server responsible

    CALL IOS_Send(qHandle,hasString=.TRUE.)

    IF ( lhook ) CALL dr_hook('IOS:CLOSE',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Close


  SUBROUTINE IOS_Fence(unit,root)

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: unit
    INTEGER, INTENT(INOUT)       :: root
    INTEGER                      :: qHandle
    INTEGER                      :: pe
    REAL(KIND=jprb)              :: zhook_handle
    CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_Fence'
    TYPE(IOS_metadata_type)      :: metadata

    IF ( lhook ) CALL dr_hook('IOS:FENCE',zhook_in,zhook_handle)

    IF ( L_IOS_active() .AND. model_rank == 0 ) THEN

!The root is the global task id owning the unit, if its owned by
!atmosphere then its set to zero. 
      IF ( root == 0 ) THEN
        root=global_rank ! because atm rank 0 is not COMM_WORLD 0
      ELSE
        root=io_server_for_unit(unit)
      END IF

      WRITE(6,'(A,I3,A,I5)')'Fencing unit ',unit,' on rank ',root

      DO pe=1,IOS_Server_Groups

        !recycle address to encode 'owner' in this case
        qHandle=IOS_Init_MD(-1*io_servers(pe,1),root,IOS_Action_Fence)
        CALL IOS_Send(qHandle)

      END DO

! If, I, atmos rank 0, own the file I have to participate in the 
! operation. Of course I have no queue, so I don't need to wait

      IF ( root == global_rank ) THEN
        WRITE(6,'(A)')'Atmosphere waiting for fence'
        CALL IOS_SoftSync(root)
      END IF

    END IF

    IF ( lhook ) CALL dr_hook('IOS:FENCE',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Fence


  SUBROUTINE IOS_DiskSync(unit,wait)

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: unit! unit number
    LOGICAL, OPTIONAL            :: wait! length of env
    INTEGER                      :: errorCode=99
    INTEGER                      :: qHandle
    LOGICAL                      :: myWaitFlag
    REAL(KIND=jprb)              :: zhook_handle
    CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_Sync'
    TYPE(IOS_metadata_type), &
        POINTER                  :: metadata

    IF ( lhook ) CALL dr_hook('IOS:DISKSYNC',zhook_in,zhook_handle)

    myWaitFlag=.TRUE. ! Default is to block (ie slow but safe)
    IF ( PRESENT(wait) ) THEN
      IF ( .NOT.wait ) THEN
        myWaitFlag=.FALSE.
      END IF
    END IF

    IF ( myWaitFlag ) THEN
      qHandle=IOS_Init_MD(unit,-1,IOS_Action_Sync_Barrier,datasize=1)      
    ELSE
      qHandle=IOS_Init_MD(unit,-1,IOS_Action_Sync)
    END IF

    CALL IOS_Send(qHandle)

    ! Return code - assume OK, errors handled by IO Server

    IF ( lhook ) CALL dr_hook('IOS:DISKSYNC',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_DiskSync


  SUBROUTINE IOS_Config(verbosity,timer)

    IMPLICIT NONE
    INTEGER,OPTIONAL           :: verbosity
    LOGICAL,OPTIONAL           :: timer
    INTEGER                    :: qHandle
    INTEGER                    :: ios
    INTEGER                    :: i
    REAL(KIND=jprb)            :: zhook_handle
    TYPE(IOS_metadata_type), &
        POINTER                :: metadata

    IF ( lhook ) CALL dr_hook('IOS:CONFIG',zhook_in,zhook_handle)

    IF (L_IOS_active()) THEN
      IF (model_rank==0) THEN

        IF (PRESENT(verbosity)) THEN        
          IF (IOS_acquire_model_prsts) THEN
            IOS_Verbosity=verbosity          
            DO i=1,IOS_Server_Groups        
              qHandle=IOS_Init_MD(-1*io_servers(i,1),-1,IOS_Action_Config)
              metadata => IOS_Attach_Metadata(qHandle)
              metadata % intent = verbosity
              CALL IOS_Send(qHandle)
            END DO
          END IF
        END IF
      
        IF (PRESENT(timer)) THEN
          IF(timer) THEN            
            DO i=1,IOS_Server_Groups        
              qHandle=IOS_Init_MD(-1*io_servers(i,1),-1,IOS_Action_Config)
              metadata => IOS_Attach_Metadata(qHandle)
              metadata % intent = IOS_Timer_Activate
              CALL IOS_Send(qHandle)
            END DO            
          END IF
        END IF

      END IF
    END IF
    
    IF ( lhook ) CALL dr_hook('IOS:CONFIG',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Config


  FUNCTION IOS_Open(unit,name,read_write,filetype,allowRemap) RESULT(r)
 
    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: unit      ! unit number
    INTEGER, INTENT(IN)          :: read_write! ==0 read only, 
    !                                         ! otherwise r/w
    INTEGER, INTENT(IN)          :: fileType
    LOGICAL, OPTIONAL            :: allowRemap
    INTEGER                      :: r
    INTEGER                      :: qHandle
    INTEGER                      :: errorCode=99
    REAL(KIND=jprb)              :: zhook_handle
    CHARACTER (LEN=*), INTENT(IN):: name ! environment var/filename
    TYPE(IOS_metadata_type), &
        POINTER                  :: metadata
    CHARACTER (LEN=*), PARAMETER :: RoutineName = 'IOS_Open'

    IF ( lhook ) CALL dr_hook('IOS:OPEN',zhook_in,zhook_handle)

    ! Check name is small enough to transmit.
    IF (LEN_TRIM(name) > IOS_String_Max) THEN
      WRITE(IOS_Message,'(A,I6,A,I6,A)')                       &
          'Filename name too big (',                           &
          LEN_TRIM(name),' > maximum of ', IOS_String_Max,')'      
      CALL IOS_Ereport(RoutineName, ErrorCode, Ios_message)
    END IF

    IF (PRESENT(allowRemap)) THEN
      IF (allowRemap) THEN
        CALL IOS_assign_server_for_unit(unit) 
      END IF
    END IF

! If we are using a dynamic allocation, then we need to 
! check and assign a server for the unit. All ranks must 
! call this
    IF (io_server_for_unit(unit) == IOS_No_Server) &
        CALL IOS_assign_Server_for_unit(unit)

    r=-1
    IF (model_rank==0) THEN
      qHandle  = IOS_Init_MD(unit,-1,IOS_Action_Open)
      metadata => IOS_Attach_Metadata(qHandle)
      metadata % name_length              = LEN_TRIM(name)
      metadata % intent                   = read_write
      metadata % address                  = fileType
      metadata % string(1:metadata % name_length) &
                                          = name(1:metadata % name_length)
      r = IOS_GetDestPe(qHandle)
      CALL IOS_Send(qHandle,hasString=.TRUE.)      
    END IF

    IF ( lhook ) CALL dr_hook('IOS:OPEN',zhook_out,zhook_handle)

    RETURN
  END FUNCTION IOS_Open


  SUBROUTINE IOS_Process(timestep)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: timestep! current timestep
    INTEGER                    :: server
    INTEGER                    :: subtask
    INTEGER                    :: qHandle
    REAL, EXTERNAL             :: get_wallclock_time
    REAL                       :: t
    REAL(KIND=jprb)            :: zhook_handle
    TYPE(IOS_metadata_type), &
        POINTER                :: metadata

    IF ( lhook ) CALL dr_hook('IOS:PROCESS',zhook_in,zhook_handle)

    IF (model_rank == 0 ) THEN
      DO server=1,IOS_Server_Groups
        DO subtask=1,IOS_tasks_per_server
          qHandle=IOS_Init_MD(-1*io_servers(server,subtask), &
              -1,IOS_Action_Process)
          CALL IOS_Send(qHandle)
        END DO
      END DO
! DEPENDS ON : get_wallclock_time
      t=get_wallclock_time()
      IF (IOS_Verbosity>=IOS_prstatus_oper) THEN
        WRITE(6,'(A,I8,A,F10.3)') &
            'IOS: Info: Entering timestep ',timestep, &
            ' at time=',t-IOS_Start_time
      END IF
    END IF

    IF ( lhook ) CALL dr_hook('IOS:PROCESS',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Process


  SUBROUTINE IOS_Read32_r(unit, location, array, numWords, RdBroadCast)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    INTEGER, INTENT(IN)        :: RdBroadCast
    REAL(KIND=real32), &
        INTENT(OUT)            :: array(:)! data read
    INTEGER                    :: qHandle
    INTEGER                    :: errorCode
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: recvBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:READ32',zhook_in,zhook_handle)

    IF (model_rank==0) THEN
      qHandle=IOS_Init_MD(unit,location,IOS_Action_Read32, &
          dataSize=numWords,bcast=RdBroadCast)
      CALL IOS_Send(qHandle)
    END IF

    IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
      recvBuffer => IOS_attach_recvBuffer(qHandle)
      CALL um_memcpy32(array,recvBuffer,numWords)
    ELSE
      CALL MPL_BCast(array,                                       &
          numWords*IOS_TUPerWord32,                               &
          IOS_TUType,                                             &
          bcast_procs-1,                                          &
          IOS_BCast_Comm(ioServerNo(io_server_for_unit(unit))),   &
          errorCode)
    END IF

    IF ( lhook ) CALL dr_hook('IOS:READ32',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Read32_r


  SUBROUTINE IOS_Read32_i(unit, location, array, numWords, RdBroadCast)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    INTEGER, INTENT(IN)        :: RdBroadCast
    INTEGER                    :: errorCode
    INTEGER(KIND=integer32), &
        INTENT(OUT)            :: array(:)! data read
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: recvBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:READ32',zhook_in,zhook_handle)

    IF (model_rank==0) THEN
      qHandle=IOS_Init_MD(unit,location,IOS_Action_Read32, &
          dataSize=numWords,bcast=RdBroadCast)
      CALL IOS_Send(qHandle)
    END IF

    IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
      recvBuffer => IOS_attach_recvBuffer(qHandle)
      CALL um_memcpy32(array,recvBuffer,numWords)
    ELSE
      CALL MPL_BCast(array,                                       &
          numWords*IOS_TUPerWord32,                               &
          IOS_TUType,                                             &
          bcast_procs-1,                                          &
          IOS_BCast_Comm(ioServerNo(io_server_for_unit(unit))),   &
          errorCode)
    END IF

    IF ( lhook ) CALL dr_hook('IOS:READ32',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Read32_i


  SUBROUTINE IOS_Read32_i2D(unit, location, array, numWords, RdBroadCast)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    INTEGER, INTENT(IN)        :: RdBroadCast
    INTEGER(KIND=integer32), &
        INTENT(OUT)            :: array(:,:)! data read
    INTEGER                    :: qHandle
    INTEGER                    :: errorCode
    INTEGER                    :: status(MPL_STATUS_SIZE)
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: recvBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:READ32',zhook_in,zhook_handle)

    IF (model_rank==0) THEN
      qHandle=IOS_Init_MD(unit,location,IOS_Action_Read32, &
          dataSize=numWords,bcast=RdBroadCast)
      CALL IOS_Send(qHandle)
    END IF

    IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
      recvBuffer => IOS_attach_recvBuffer(qHandle)
      CALL um_memcpy32(array,recvBuffer,numWords)
    ELSE
      CALL MPL_BCast(array,                                       &
          numWords*IOS_TUPerWord32,                               &
          IOS_TUType,                                             &
          bcast_procs-1,                                          &
          IOS_BCast_Comm(ioServerNo(io_server_for_unit(unit))),   &
          errorCode)
    END IF

    IF ( lhook ) CALL dr_hook('IOS:READ32',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Read32_i2D


  SUBROUTINE IOS_Read32_l(unit, location, array, numWords, RdBroadCast)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    INTEGER, INTENT(IN)        :: RdBroadCast
    LOGICAL(KIND=logical32), &
        INTENT(OUT)            :: array(:)! data read
    INTEGER                    :: qHandle
    INTEGER                    :: errorCode
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: recvBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:READ32',zhook_in,zhook_handle)

    IF (model_rank==0) THEN
      qHandle=IOS_Init_MD(unit,location,IOS_Action_Read32, &
          dataSize=numWords,bcast=RdBroadCast)
      CALL IOS_Send(qHandle)
    END IF

    IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
      recvBuffer => IOS_attach_recvBuffer(qHandle)
      CALL um_memcpy32(array,recvBuffer,numWords)
    ELSE
      CALL MPL_BCast(array,                                       &
          numWords*IOS_TUPerWord32,                               &
          IOS_TUType,                                             &
          bcast_procs-1,                                          &
          IOS_BCast_Comm(ioServerNo(io_server_for_unit(unit))),   &
          errorCode)
    END IF

    IF ( lhook ) CALL dr_hook('IOS:READ32',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Read32_l

  
  SUBROUTINE IOS_Read64_r(unit, location, array, numWords, RdBroadCast)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location    
    INTEGER, INTENT(IN)        :: numWords
    INTEGER, INTENT(IN)        :: RdBroadCast
    REAL(KIND=real64), &
        INTENT(OUT)            :: array(:)! data read
    INTEGER                    :: qHandle
    INTEGER                    :: errorCode
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: recvBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:READ',zhook_in,zhook_handle)

    IF (model_rank==0) THEN
      qHandle=IOS_Init_MD(unit,location,IOS_Action_Read64, &
          dataSize=numWords,bcast=RdBroadCast)
      CALL IOS_Send(qHandle)
    END IF

    IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
      recvBuffer => IOS_attach_recvBuffer(qHandle)
      CALL um_memcpy64(array,recvBuffer,numWords)
    ELSE
      CALL MPL_BCast(array,                                       &
          numWords*IOS_TUPerWord64,                               &
          IOS_TUType,                                             &
          bcast_procs-1,                                          &
          IOS_BCast_Comm(ioServerNo(io_server_for_unit(unit))),   &
          errorCode)
    END IF

    IF ( lhook ) CALL dr_hook('IOS:READ',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Read64_r


  SUBROUTINE IOS_Read64_l(unit, location, array, numWords, RdBroadCast)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    INTEGER, INTENT(IN)        :: RdBroadCast
    LOGICAL(KIND=logical64), &
        INTENT(OUT)            :: array(:)! data read
    INTEGER                    :: qHandle
    INTEGER                    :: errorCode
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: recvBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:READ',zhook_in,zhook_handle)

    IF (model_rank==0) THEN
      qHandle=IOS_Init_MD(unit,location,IOS_Action_Read64, &
          dataSize=numWords,bcast=RdBroadCast)
      CALL IOS_Send(qHandle)
    END IF

    IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
      recvBuffer => IOS_attach_recvBuffer(qHandle)
      CALL um_memcpy64(array,recvBuffer,numWords)
    ELSE
      CALL MPL_BCast(array,                                       &
          numWords*IOS_TUPerWord64,                               &
          IOS_TUType,                                             &
          bcast_procs-1,                                          &
          IOS_BCast_Comm(ioServerNo(io_server_for_unit(unit))),   &
          errorCode)
    END IF

    IF ( lhook ) CALL dr_hook('IOS:READ',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Read64_l


  SUBROUTINE IOS_Read64_i(unit, location, array, numWords, RdBroadCast)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    INTEGER, INTENT(IN)        :: RdBroadCast
    INTEGER(KIND=integer64), &
        INTENT(OUT)            :: array(:)! data read
    INTEGER                    :: qHandle
    INTEGER                    :: errorCode
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: recvBuffer(:)
    TYPE(IOS_metadata_type), &
        POINTER                :: metadata

    IF ( lhook ) CALL dr_hook('IOS:READ',zhook_in,zhook_handle)

    IF (model_rank==0) THEN
      qHandle=IOS_Init_MD(unit,location,IOS_Action_Read64, &
          dataSize=numWords,bcast=RdBroadCast)
      CALL IOS_Send(qHandle)
    END IF

    IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
      recvBuffer => IOS_attach_recvBuffer(qHandle)
      CALL um_memcpy64(array,recvBuffer,numWords)
    ELSE
      CALL MPL_BCast(array,                                       &
          numWords*IOS_TUPerWord64,                               &
          IOS_TUType,                                             &
          bcast_procs-1,                                          &
          IOS_BCast_Comm(ioServerNo(io_server_for_unit(unit))),   &
          errorCode)
    END IF

    IF ( lhook ) CALL dr_hook('IOS:READ',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Read64_i


  SUBROUTINE IOS_Read64_r2D(unit, location, array, numWords, RdBroadCast)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    INTEGER, INTENT(IN)        :: RdBroadCast
    REAL(KIND=real64), &
        INTENT(OUT)            :: array(:,:)! data read
    INTEGER                    :: qHandle
    INTEGER                    :: errorCode
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: recvBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:READ',zhook_in,zhook_handle)

    IF (model_rank==0) THEN
      qHandle=IOS_Init_MD(unit,location,IOS_Action_Read64, &
          dataSize=numWords,bcast=RdBroadCast)
      CALL IOS_Send(qHandle)
    END IF

    IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
      recvBuffer => IOS_attach_recvBuffer(qHandle)
      CALL um_memcpy64(array,recvBuffer,numWords)
    ELSE
      CALL MPL_BCast(array,                                       &
          numWords*IOS_TUPerWord64,                               &
          IOS_TUType,                                             &
          bcast_procs-1,                                          &
          IOS_BCast_Comm(ioServerNo(io_server_for_unit(unit))),   &
          errorCode)
    END IF

    IF ( lhook ) CALL dr_hook('IOS:READ',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Read64_r2D


  SUBROUTINE IOS_Read64_l2D(unit, location, array, numWords, RdBroadCast)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    INTEGER, INTENT(IN)        :: RdBroadCast
    LOGICAL(KIND=logical64), &
        INTENT(OUT)            :: array(:,:)! data read
    INTEGER                    :: qHandle
    INTEGER                    :: errorCode
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: recvBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:READ',zhook_in,zhook_handle)

    IF (model_rank==0) THEN
      qHandle=IOS_Init_MD(unit,location,IOS_Action_Read64, &
          dataSize=numWords,bcast=RdBroadCast)
      CALL IOS_Send(qHandle)
    END IF

    IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
      recvBuffer => IOS_attach_recvBuffer(qHandle)
      CALL um_memcpy64(array,recvBuffer,numWords)
    ELSE
      CALL MPL_BCast(array,                                       &
          numWords*IOS_TUPerWord64,                               &
          IOS_TUType,                                             &
          bcast_procs-1,                                          &
          IOS_BCast_Comm(ioServerNo(io_server_for_unit(unit))),   &
          errorCode)
    END IF

    IF ( lhook ) CALL dr_hook('IOS:READ',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Read64_l2D


  SUBROUTINE IOS_Read64_i2D(unit, location, array, numWords, RdBroadCast)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    INTEGER, INTENT(IN)        :: RdBroadCast
    INTEGER(KIND=integer64), &
        INTENT(OUT)            :: array(:,:)! data read
    INTEGER                    :: qHandle
    INTEGER                    :: errorCode
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: recvBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:READ',zhook_in,zhook_handle)

    IF (model_rank==0) THEN
      qHandle=IOS_Init_MD(unit,location,IOS_Action_Read64, &
          dataSize=numWords,bcast=RdBroadCast)
      CALL IOS_Send(qHandle)
    END IF

    IF (RdBroadCast == IOS_Read_NoBroadcast) THEN
      recvBuffer => IOS_attach_recvBuffer(qHandle)
      CALL um_memcpy64(array,recvBuffer,numWords)
    ELSE
      CALL MPL_BCast(array,                                       &
          numWords*IOS_TUPerWord64,                               &
          IOS_TUType,                                             &
          bcast_procs-1,                                          &
          IOS_BCast_Comm(ioServerNo(io_server_for_unit(unit))),   &
          errorCode)
    END IF

    IF ( lhook ) CALL dr_hook('IOS:READ',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Read64_i2D


  SUBROUTINE IOS_Setpos(unit, pos, icode)
! Deprecated? completely obsolete?
    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: pos! file seek position
    INTEGER, INTENT(OUT)       :: icode! error flag
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    TYPE(IOS_metadata_type)    :: metadata

    IF ( lhook ) CALL dr_hook('IOS:SETPOS',zhook_in,zhook_handle)

    ! Set metadata
    qHandle=IOS_Init_MD(unit,pos,IOS_Action_Setpos)    
    CALL IOS_Send(qHandle)
    icode = 0

    IF ( lhook ) CALL dr_hook('IOS:SETPOS',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Setpos


  SUBROUTINE IOS_Write32_r(unit, location, array, numWords)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    REAL(KIND=real32),&
        INTENT(IN)             :: array(:)! data to write
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: sendBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:WRITE32',zhook_in,zhook_handle)

    qHandle=IOS_Init_MD(unit,location, & 
        IOS_Action_Write32, dataSize=numWords)
    sendBuffer => IOS_attach_SendBuffer(qHandle)
    CALL um_memcpy32(sendBuffer,array,numWords)
    CALL IOS_Send(qHandle)
    
    IF ( lhook ) CALL dr_hook('IOS:WRITE32',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Write32_r


  SUBROUTINE IOS_Write32_i(unit, location, array, numWords)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    INTEGER(KIND=integer32), &
        INTENT(IN)             :: array(:)! data to write
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: sendBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:WRITE32',zhook_in,zhook_handle)

    qHandle=IOS_Init_MD(unit,location, & 
        IOS_Action_Write32, dataSize=numWords)
    sendBuffer => IOS_attach_SendBuffer(qHandle)
    CALL um_memcpy32(sendBuffer,array,numWords)
    CALL IOS_Send(qHandle)

    IF ( lhook ) CALL dr_hook('IOS:WRITE32',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Write32_i


  SUBROUTINE IOS_Write32_i2D(unit, location, array, numWords)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    INTEGER(KIND=integer32),&
          INTENT(IN)           :: array(:,:)! data to write
    INTEGER                    :: qHandle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: sendBuffer(:)
    REAL(KIND=jprb)            :: zhook_handle

    IF ( lhook ) CALL dr_hook('IOS:WRITE32',zhook_in,zhook_handle)

    qHandle=IOS_Init_MD(unit,location, & 
        IOS_Action_Write32, dataSize=numWords)
    sendBuffer => IOS_attach_SendBuffer(qHandle)
    CALL um_memcpy32(sendBuffer,array,numWords)
    CALL IOS_Send(qHandle)

    IF ( lhook ) CALL dr_hook('IOS:WRITE32',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Write32_i2D


  SUBROUTINE IOS_Write64_r(unit, location, array, numWords)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    REAL(KIND=real64), &
        INTENT(IN)             :: array(:)! data to write
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: sendBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:WRITE',zhook_in,zhook_handle)

    qHandle=IOS_Init_MD(unit,location,IOS_Action_Write64,dataSize=numWords)
    sendBuffer => IOS_attach_SendBuffer(qHandle)
    CALL um_memcpy64(sendBuffer,array,numWords)
    CALL IOS_Send(qHandle)

    IF ( lhook ) CALL dr_hook('IOS:WRITE',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Write64_r


  SUBROUTINE IOS_Write64_i(unit, location, array, numWords)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    INTEGER(KIND=integer64), &
        INTENT(IN)             :: array(:)! data to write
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: sendBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:WRITE',zhook_in,zhook_handle)

    qHandle=IOS_Init_MD(unit,location,IOS_Action_Write64,dataSize=numWords)
    sendBuffer => IOS_attach_SendBuffer(qHandle)
    CALL um_memcpy64(sendBuffer,array,numWords)
    CALL IOS_Send(qHandle)

    IF ( lhook ) CALL dr_hook('IOS:WRITE',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Write64_i


  SUBROUTINE IOS_Write64_l(unit, location, array, numWords)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    LOGICAL(KIND=logical64), &
        INTENT(IN)             :: array(:)! data to write
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: sendBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:WRITE',zhook_in,zhook_handle)

    qHandle=IOS_Init_MD(unit,location,IOS_Action_Write64,dataSize=numWords)
    sendBuffer => IOS_attach_SendBuffer(qHandle)
    CALL um_memcpy64(sendBuffer,array,numWords)
    CALL IOS_Send(qHandle)

    IF ( lhook ) CALL dr_hook('IOS:WRITE',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Write64_l


  !2D...............
  SUBROUTINE IOS_Write64_r2D(unit, location, array, numWords)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    REAL(KIND=real64), &
        INTENT(IN)             :: array(:,:)! data to write
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: sendBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:WRITE',zhook_in,zhook_handle)

    qHandle=IOS_Init_MD(unit,location,IOS_Action_Write64,dataSize=numWords)
    sendBuffer => IOS_attach_SendBuffer(qHandle)
    CALL um_memcpy64(sendBuffer,array,numWords)
    CALL IOS_Send(qHandle)

    IF ( lhook ) CALL dr_hook('IOS:WRITE',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Write64_r2D


  SUBROUTINE IOS_Write64_i2D(unit, location, array, numWords)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    INTEGER(KIND=integer64), &
        INTENT(IN)             :: array(:,:)! data to write
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: sendBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:WRITE',zhook_in,zhook_handle)

    qHandle=IOS_Init_MD(unit,location,IOS_Action_Write64,dataSize=numWords)
    sendBuffer => IOS_attach_SendBuffer(qHandle)
    CALL um_memcpy64(sendBuffer,array,numWords)
    CALL IOS_Send(qHandle)

    IF ( lhook ) CALL dr_hook('IOS:WRITE',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Write64_i2D


  SUBROUTINE IOS_Write64_l2D(unit, location, array, numWords)

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: numWords
    LOGICAL(KIND=logical64), &
        INTENT(IN)             :: array(:,:)! data to write
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: sendBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:WRITE',zhook_in,zhook_handle)

    qHandle=IOS_Init_MD(unit,location,IOS_Action_Write64,dataSize=numWords)
    sendBuffer => IOS_attach_SendBuffer(qHandle)
    CALL um_memcpy64(sendBuffer,array,numWords)
    CALL IOS_Send(qHandle)

    IF ( lhook ) CALL dr_hook('IOS:WRITE',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Write64_l2D


  SUBROUTINE IOS_Init_Header(unit,len)

    USE IOS_STASH_COMMON, ONLY : isUsingAsyncStash

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit
    INTEGER, INTENT(IN)        :: len
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle

    IF ( lhook ) CALL dr_hook('IOS:INIT_HEADER',zhook_in,zhook_handle)

    IF (model_rank==0.AND.L_IOS_Active().AND.isUsingAsyncStash()) THEN
      qHandle=IOS_Init_MD(unit,len,IOS_Action_StashInitHeader)
      CALL IOS_Send(qHandle)   
    END IF

    IF ( lhook ) CALL dr_hook('IOS:INIT_HEADER',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Init_Header


  SUBROUTINE IOS_Init_PP_Lookup(args)

    USE IOS_STASH_COMMON, ONLY : isUsingAsyncStash

    IMPLICIT NONE
    INTEGER, PARAMETER         :: num_args=5
    INTEGER, INTENT(IN)        :: args(num_args)
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: sendBuffer(:)

    IF ( lhook ) &
        CALL dr_hook('IOS:INIT_PP_LOOKUP',zhook_in,zhook_handle)

    IF ( model_rank == 0 .AND. &
         L_IOS_Active()  .AND. &
         isUsingAsyncStash()) THEN

      qHandle=IOS_Init_MD(args(1),-1,IOS_Action_StashInitPPLookup, &
          dataSize=num_args)
      sendBuffer => IOS_attach_SendBuffer(qHandle)
      CALL um_memcpy64(sendBuffer,args,num_args)
      CALL IOS_Send(qHandle)

    END IF

    IF ( lhook ) &
        CALL dr_hook('IOS:INIT_PP_LOOKUP',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Init_PP_Lookup


  SUBROUTINE IOS_Set_Header(unit,header_data)

    USE IOS_STASH_COMMON, ONLY : isUsingAsyncStash

    IMPLICIT NONE
    INTEGER, INTENT(IN)          :: unit
    INTEGER, POINTER             :: header_data(:)
    INTEGER                      :: qHandle
    REAL(KIND=jprb)              :: zhook_handle
    INTEGER(KIND=IOS_TUKind), & 
        POINTER                  :: sendBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:SET_HEADER',zhook_in,zhook_handle)

    IF ( model_rank == 0 .AND.   &
         L_IOS_Active()  .AND.   &
         (.NOT. L_IO_SERVER) .AND. &
         isUsingAsyncStash()) THEN
      
      qHandle=IOS_Init_MD(unit,-1,IOS_Action_StashSetHeader, &
          dataSize=SIZE(header_data))
      sendBuffer => IOS_attach_SendBuffer(qHandle)
      CALL um_memcpy64(sendBuffer,header_data,SIZE(header_data))
      CALL IOS_Send(qHandle)

    END IF

    IF ( lhook ) CALL dr_hook('IOS:SET_HEADER',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Set_Header


  SUBROUTINE IOS_Setpos_Stash(unit, record, icode)

    USE IOS_STASH_COMMON, ONLY : isUsingAsyncStash

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: record! file seek position
    INTEGER, INTENT(OUT)       :: icode! error flag
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle

    IF ( lhook ) CALL dr_hook('IOS:SETPOS_STASH',zhook_in,zhook_handle)

    IF ( model_rank == 0 .AND. &
         L_IOS_Active()  .AND. &
         isUsingAsyncStash()) THEN

      qHandle=IOS_Init_MD(unit,record,IOS_Action_StashSetPos)
      CALL IOS_Send(qHandle)

    END IF

    icode = 0

    IF ( lhook ) CALL dr_hook('IOS:SETPOS_STASH',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Setpos_Stash


  SUBROUTINE IOS_MergeLookup(unit, array,isize)

    USE IOS_STASH_COMMON, ONLY : isUsingAsyncStash

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: isize! data size
    INTEGER, INTENT(IN)        :: array(:,:)! data to write
    INTEGER                    :: io_slave_pe
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    TYPE(IOS_metadata_type)    :: metadata
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: sendBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:MERGELOOKUP', &
        zhook_in,zhook_handle)

    IF ( model_rank == 0 .AND. &
         L_IOS_Active()  .AND. &
         isUsingAsyncStash()) THEN

      qHandle=IOS_Init_MD(unit,IOS_No_Location,IOS_ACTION_MergePPLookup, &
          dataSize=isize)
      sendBuffer => IOS_attach_SendBuffer(qHandle)
      CALL um_memcpy64(sendBuffer,array,isize)
      CALL IOS_Send(qHandle)
      
    END IF

    IF ( lhook ) CALL dr_hook('IOS:MERGELOOKUP', &
        zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_MergeLookup


  SUBROUTINE IOS_Write_Lookup(unit, location, array)

    USE IOS_STASH_COMMON, ONLY : isUsingAsyncStash

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: location
    INTEGER, INTENT(IN)        :: array(:,:)! data to write
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: sendBuffer(:)

    IF ( lhook ) CALL dr_hook('IOS:WRITE_LOOKUP', &
        zhook_in,zhook_handle)

    IF ( model_rank == 0 .AND. &
         L_IOS_Active()  .AND. &
         isUsingAsyncStash()) THEN

      qHandle=IOS_Init_MD(unit,location,IOS_ACTION_StashWritePPLookup, &
          dataSize=SIZE(array))
      sendBuffer => IOS_attach_SendBuffer(qHandle)
      CALL um_memcpy64(sendBuffer,array,SIZE(array))
      CALL IOS_Send(qHandle)
      
    END IF
    
    IF ( lhook ) CALL dr_hook('IOS:WRITE_LOOKUP', &
        zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Write_Lookup


  SUBROUTINE IOS_Write_Stash_PrePacked_Data &
      (unit ,record, array)

    USE IOS_STASH_COMMON, ONLY : isUsingAsyncStash

    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: unit! unit
    INTEGER, INTENT(IN)        :: record! the stash record
    REAL, INTENT(IN)           :: array(:)! data to write
    INTEGER                    :: qHandle
    REAL(KIND=jprb)            :: zhook_handle
    INTEGER(KIND=IOS_TUKind), &
        POINTER                :: sendBuffer(:)
    TYPE(IOS_metadata_type), &
        POINTER                :: metadata

    IF ( lhook ) CALL dr_hook('IOS:WRITE_STASH_PREPACKED_DATA', &
        zhook_in,zhook_handle)

    IF ( model_rank == 0 .AND. &
         L_IOS_Active()  .AND. &
         isUsingAsyncStash()) THEN

      ! Set metadata (note we use location to encode record)
      qHandle=IOS_Init_MD(unit,record,          &
          IOS_Action_Write_PP_Prepacked,        &
          datasize = SIZE(array))
      sendBuffer => IOS_attach_SendBuffer(qHandle)
      CALL um_memcpy64(sendBuffer,array,SIZE(array))
      CALL IOS_Send(qHandle)

    END IF
    
    IF ( lhook ) CALL dr_hook('IOS:WRITE_STASH_PREPACKED_DATA', &
        zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE IOS_Write_Stash_PrePacked_Data

END MODULE ios

