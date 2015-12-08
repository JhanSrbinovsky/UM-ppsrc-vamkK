! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

! Module defining interfaces that might be called by client or server tasks

MODULE IOS_Common
  USE IOS_types
  USE IOS_constants
  USE IOS_communicators
  USE missing_data_mod, ONLY: IMDI
  USE yomhook, ONLY: lhook, dr_hook  
  USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

  ! params/vars  for dr_hook
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1


  INTEGER :: IOS_Server_Groups=0
  INTEGER :: IOS_Tasks_Per_Server=IMDI  ! namelist item

! Sequencing: These items track the order of operations to ensure
! communications are always processed in order. There is a different 
! sequence for each IO server, that the client must track.
 
  INTEGER, PARAMETER :: IOS_Sequence_Max = 100000
  INTEGER, POINTER   :: IOS_Sequence_ID(:,:)
  
! Internally defined I/O Processor for a unit
  INTEGER            :: io_server_for_unit_lookup &
      (IOS_unit_first:IOS_unit_last)
  
! Which global ranks in the application do IO
  INTEGER, POINTER   :: io_servers(:,:) => NULL()

! Requests for MPI help for an operation.
  INTEGER            :: IOS_ReadCompletionRequested
  INTEGER            :: IOS_EnqCompletionRequested

! If I am an io server
  LOGICAL            :: l_io_server = .FALSE.
  LOGICAL            :: l_io_leader = .FALSE.
! File status 
  TYPE(IOS_STATUS)   :: IOS_unit_status

! Start time
  REAL               :: IOS_start_time

! WordLengths
  INTEGER            :: IOS_BytesPerReal
  INTEGER            :: IOS_BytesPerInteger

!Control the behaviour of IOS
!  Non-namelist items
  INTEGER            :: threading_model
  LOGICAL            :: use_blocking_recvs       = .FALSE. ! for debugging
  LOGICAL            :: serialize_all_ops                  ! for debugging
!  Namelist items
  INTEGER            :: IOS_backoff_interval     = IMDI
  INTEGER            :: IOS_Unit_Alloc_Policy    = IMDI
  INTEGER            :: IOS_concurrency          = IMDI
  INTEGER            :: IOS_Verbosity            = IMDI
  INTEGER            :: IOS_Timeout              = IMDI
!            How long we sit and wait before giving up and aborting the model
  LOGICAL            :: IOS_local_ro_files       = .FALSE.
  LOGICAL            :: IOS_acquire_model_prsts  = .FALSE.
  LOGICAL            :: IOS_serialise_mpi_calls  = .FALSE.
  LOGICAL            :: IOS_thread_0_calls_mpi   = .FALSE.
  LOGICAL            :: IOS_Lock_Meter           = .FALSE.

CONTAINS

  LOGICAL FUNCTION assert_client()
    IMPLICIT NONE
    INTEGER                         :: ErrorCode = 99
    CHARACTER (LEN=*), PARAMETER    :: RoutineName = &
        'IOS_Common:ASSERT_CLIENT'

    assert_client=.TRUE.
    IF (l_io_server) THEN
      assert_client=.FALSE.
      CALL IOS_ereport( RoutineName, ErrorCode, &
          'ASSERT_CLIENT FAILURE: ROUTINE CALLED BY IO SERVER PROCESS' )      
    END IF
  END FUNCTION assert_client

  LOGICAL FUNCTION assert_server()
    IMPLICIT NONE
    INTEGER                         :: ErrorCode = 99
    CHARACTER (LEN=*), PARAMETER    :: &
        RoutineName = 'IOS_Common:ASSERT_SERVER'
    assert_server=.TRUE.
    IF (.NOT.l_io_server) THEN
      assert_server=.FALSE.
      CALL IOS_ereport( RoutineName, ErrorCode, &
          'ASSERT_SERVER FAILURE: ROUTINE CALLED BY IO CLIENT PROCESS')
    END IF
  END FUNCTION assert_server
  
  SUBROUTINE IOS_ereport(r,code,m,md,unit)
    USE ereport_mod, ONLY : ereport
    IMPLICIT NONE
    INTEGER, INTENT(IN)               :: code
    TYPE(IOS_metadata_type), OPTIONAL :: md
    INTEGER, OPTIONAL                 :: unit
    CHARACTER(LEN=*), INTENT(IN)        :: r,m
    INTEGER                           :: lcode

    lcode=code 
    WRITE(6,*)'-------- IOS ERROR REPORT ---------------'
    IF (PRESENT(unit)) THEN
      WRITE(6,'(A,I8)')'Problem with unit ',unit
    ELSE IF (PRESENT(md)) THEN
      IF (md%unit > 0) &
          WRITE(6,'(A,I8)')'Problem with unit ',md%unit
    END IF

    IF (PRESENT(md))CALL IOS_REPORT_MD(md)
    CALL ereport(r,lcode,m)

  END SUBROUTINE IOS_ereport

  FUNCTION IOS_ActionName(i) RESULT(r)
    IMPLICIT NONE
    INTEGER, INTENT(IN)              :: i
    CHARACTER(LEN=IOS_Action_Strlen) :: r    

    WRITE(r,'(A)')IOS_Action_strings((i-1)*&
        IOS_Action_Strlen+1:i*IOS_Action_Strlen)
  END FUNCTION IOS_ActionName

  SUBROUTINE IOS_Report_MD(md)
    IMPLICIT NONE
    TYPE(IOS_metadata_type), INTENT(IN)  :: md
    WRITE(6,*)'------------ IOS METADATA REPORT -----------'
    WRITE(6,*)'   md%action=',md%action
    WRITE(6,*)'   md%unit=',md%unit
    WRITE(6,*)'   md%name_length=',md%name_length
    WRITE(6,*)'   md%intent=',md%intent
    WRITE(6,*)'   md%delete=',md%delete
    WRITE(6,*)'   md%data_size=',md%data_size
    WRITE(6,*)'   md%address=',md%address
    WRITE(6,*)'   md%client=',md%client
    WRITE(6,*)'   md%handle=',md%handle
    WRITE(6,*)'   md%Originating_Slot=',md%Originating_Slot
    WRITE(6,*)'   md%string=',md%string
    WRITE(6,*)'------------ IOS METADATA REPORT -----------'
  END SUBROUTINE IOS_Report_MD

  LOGICAL FUNCTION L_IOS_active()
    IMPLICIT NONE
    L_IOS_active=ASSOCIATED(io_servers)
  END FUNCTION L_IOS_active

  INTEGER FUNCTION io_server_for_unit(unit)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: unit
    INTEGER             :: errorCode=99
    CHARACTER (LEN=80)  :: Iosc_message

    IF (unit < IOS_unit_first.OR.unit > IOS_unit_last) THEN
      WRITE(iosc_message,'(A,I8,A,I3,A,I3,A)')               &
          'Supplied unit (',unit,') outside allowed range [',&
          IOS_unit_first,' - ',IOS_unit_last,']'
      CALL IOS_ereport("IOS_Common:io_server_for_unit",&
          errorCode,iosc_message)
    END IF

    IF (.NOT.l_IOS_Active()) THEN
      CALL IOS_ereport("IOS_Common:io_server_for_unit",&
          errorCode,'IOS appears deactivated')
    END IF

    io_server_for_unit=io_server_for_unit_lookup(unit)

  END FUNCTION io_server_for_unit


! Returns the server rank ID (i.e. 1..numServers) from
! a global rank

  FUNCTION ioServerNo(globRank) RESULT(server)
    IMPLICIT NONE
    INTEGER, INTENT(IN)              :: globRank
    INTEGER                          :: server
    INTEGER                          :: i
    INTEGER                          :: j
    INTEGER                          :: errorFlag
    REAL(KIND=jprb)                  :: zhook_handle
    CHARACTER (LEN=80)               :: Iosc_message
    
    IF (lhook) CALL dr_hook('IOS_COMMON:IOSERVERNO', &
        zhook_in,zhook_handle)
    
    errorFlag=99
    server=-1
    
    DO i=1,IOS_Server_Groups
       DO j=1,IOS_tasks_per_server
          IF (globRank==io_servers(i,j)) &
              server=i
       END DO
    END DO
    
    IF (server==-1) THEN
      WRITE(IOSC_message,'(A,I5,A)') &
          'Rank ',globRank,' is not an IO Server'
      CALL IOS_Ereport('IOS_COMMON:IOSERVERNO', &
          errorFlag , IOSC_message )      
    END IF
    
    IF (lhook) CALL dr_hook('IOS_COMMON:IOSERVERNO', &
        zhook_out,zhook_handle)
    
  END FUNCTION ioServerNo
  
  FUNCTION ioServerRank(globRank) RESULT(rank)
    IMPLICIT NONE
    INTEGER, INTENT(IN)              :: globRank
    INTEGER                          :: rank
    INTEGER                          :: i
    INTEGER                          :: j
    INTEGER                          :: errorFlag
    REAL(KIND=jprb)                  :: zhook_handle
    CHARACTER (LEN=80)               :: Iosc_message
    
    IF (lhook) CALL dr_hook('IOS_COMMON:IORANKNO', &
        zhook_in,zhook_handle)
    
    errorFlag=99
    rank=-1
    
    DO i=1,IOS_Server_Groups
       DO j=1,IOS_tasks_per_server
          IF (globRank==io_servers(i,j)) &
              rank=j-1
       END DO
    END DO
    
    IF (rank==-1) THEN
      WRITE(IOSC_message,'(A,I5,A)') &
          'Rank ',globRank,' is not an IO Server'
      CALL IOS_Ereport('IOS_COMMON:IORANKNO', &
          errorFlag , IOSC_message )      
    END IF
    
    IF (lhook) CALL dr_hook('IOS_COMMON:IORANKNO', &
        zhook_out,zhook_handle)
    
  END FUNCTION ioServerRank


  SUBROUTINE IOS_SoftSync(root)
    USE MPL
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: root
    INTEGER             :: server_pe
    INTEGER             :: i
    INTEGER             :: ierror
    INTEGER             :: req_num
    INTEGER             :: flag(1)
    INTEGER,ALLOCATABLE :: requests(:)
    INTEGER,ALLOCATABLE :: status(:,:)
    REAL(KIND=jprb)     :: zhook_handle

    IF (lhook) CALL dr_hook('IOS_SOFTSYNC',zhook_in,zhook_handle)

    IF (root<0 .OR. root>= global_procs) THEN
      CALL IOS_ereport('IOS_SoftSync', 99, &
          'Invalid root process provided' )      
    END IF

    flag(1)=1
    IF (global_rank==root) THEN
      ALLOCATE(requests(IOS_Server_Groups))
      ALLOCATE(status(MPL_STATUS_SIZE,IOS_Server_Groups))
      req_num=0
      DO i=1,IOS_Server_Groups
        server_pe=io_servers(i,1)        
        IF (global_rank/=server_pe) THEN
          req_num=req_num+1
          CALL MPL_iRecv(flag,1,MPL_INTEGER,&
              server_pe,IOS_ControlSync_Tag,&
              global_comm,requests(req_num),ierror)          
        END IF
      END DO
      IF (IOS_Verbosity>=IOS_PrStatus_Oper) &
          WRITE(6,'(A,I3,A)')'IOS: Info: SoftSync: Waiting for ', &
          req_num,' other IOS to pass epoch...'
      CALL MPL_WaitAll(req_num,requests,status,ierror) 
      IF (IOS_Verbosity>=IOS_PrStatus_Oper) & 
          WRITE(6,'(A)')'IOS: Info: epoch complete.'
      DEALLOCATE(requests)
      DEALLOCATE(status)
    ELSE 
      CALL MPL_SEND(flag,1,MPL_INTEGER,root,&
          IOS_ControlSync_Tag,global_comm,ierror)
    END IF

    IF (lhook) CALL dr_hook('IOS_SOFTSYNC',zhook_out,zhook_handle)

  END SUBROUTINE IOS_SoftSync

END MODULE IOS_Common
