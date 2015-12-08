! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: C96

MODULE IOS_Queue_Mod
  USE um_types
  USE IOS_Common
  USE mpl, ONLY : mpl_comm_world
!$ USE omp_lib
  USE yomhook,                ONLY :                                           &
      lhook,                                                                   &
      dr_hook
  USE parkind1,               ONLY :                                           &
      jprb,                                                                    &
      jpim

  IMPLICIT NONE

!-----------------------------------------------------------------------
! TYPE definitions
!-----------------------------------------------------------------------
  TYPE IOS_node_type
     ! Queue related
     INTEGER                           :: transaction_number
     INTEGER                           :: payloadTag
     TYPE(IOS_metadata_type)           :: metadata
     TYPE(IOS_node_type), POINTER      :: next
     TYPE(IOS_node_type), POINTER      :: prev
     
     ! Payload related:
     INTEGER, POINTER                  :: integer_data(:)
     INTEGER(KIND=integer32), POINTER  :: integer32_data(:)
     INTEGER, POINTER                  :: receive_requests(:)
     INTEGER, POINTER                  :: receive_tags(:)
     INTEGER, POINTER                  :: receive_data_len(:)
     REAL, POINTER                     :: distributed_data(:,:)
     REAL(KIND=real32), POINTER        :: real32_data(:)
     REAL, POINTER                     :: real_data(:)
  END TYPE IOS_node_type

!-----------------------------------------------------------------------
! Common variables
!-----------------------------------------------------------------------

  INTEGER, PRIVATE, PARAMETER   :: logunit                    = 9
  INTEGER, PRIVATE              :: IOS_transaction_number     = 0
  INTEGER, PRIVATE              :: IOS_transactions_completed = 0
  INTEGER, PRIVATE              :: IOS_buffer_max




  INTEGER, PRIVATE              :: Q_size                     = 0
  INTEGER, PRIVATE              :: Q_items                    = 0

  INTEGER                       :: IOS_AsyncCompletionRequested
  INTEGER (KIND=omp_lock_kind)  :: Q_lock_var
  INTEGER (KIND=omp_lock_kind)  :: IOS_comms_access_lock
  TYPE(IOS_node_type), POINTER  :: Q_first
  TYPE(IOS_node_type), POINTER  :: Q_last
  REAL                          :: get_wallclock_time
  CHARACTER (LEN=132)           :: qmessage
  CHARACTER (LEN=*),                                                           &
      PARAMETER, PRIVATE        :: RoutineName = 'IOS_queue'

! Lock metering
  INTEGER, PARAMETER, PRIVATE   :: words_per_cacheline        = 64
  INTEGER, PARAMETER, PRIVATE   :: num_locks                  = 2
  INTEGER, PARAMETER, PRIVATE   :: mpi_lock                   = 1
  INTEGER, PARAMETER, PRIVATE   :: Q_lock                     = 2
  INTEGER, POINTER, PRIVATE     :: IOS_Q_LockTries(:,:,:)
  REAL, POINTER, PRIVATE        :: IOS_Q_LockTime (:,:,:)

  EXTERNAL get_wallclock_time

  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in          = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out         = 1

CONTAINS

!-----------------------------------------------------------------------
! SUBROUTINE to put node into Queue (node allocated outside)
!-----------------------------------------------------------------------
  SUBROUTINE IOS_Put_Node_In_Queue(node)

    IMPLICIT NONE
    TYPE(IOS_node_type), POINTER :: node
    INTEGER                      :: errCode
    LOGICAL                      :: dummy

    REAL(KIND=jprb):: zhook_handle
    IF ( lhook ) CALL dr_hook('IOS_QUEUE:PUTNODEINQUEUE',                      &
        zhook_in,zhook_handle)

    IF ( .NOT.ASSOCIATED(node) ) THEN
      WRITE(qmessage,'(A)')'Tried to put an unassociated POINTER on the Q'
      CALL IOS_Ereport( RoutineName,99, Qmessage )
    END IF

    node%transaction_number=IOS_transaction_number
    IOS_transaction_number=IOS_transaction_number+1

    IF ( serialize_all_ops ) THEN
      IF ( IOS_thread_0_calls_mpi ) THEN
        WRITE(6,'(A,A)')                                                       &
            'IOS_queue: Warning: Cant serialize operations',                   &
            ' in this MPI mode'
      ELSE
        dummy=WaitForDrain()
      END IF
    END IF

    CALL acquire_lock()
    CALL IOS_Increment_Queue_Length(ios_md_len*umFortranIntegerSize(),         &
        need_lock=.FALSE.)
!$OMP FLUSH
    IF ( .NOT. ASSOCIATED( Q_first) ) THEN
      Q_first => node
!$OMP FLUSH
    END IF

    IF ( ASSOCIATED( Q_last ) ) THEN
      Q_last%next => node
      node%prev => Q_last
      NULLIFY (node%next)
    END IF

    Q_last => node
!$OMP FLUSH
    Q_items=Q_items+1
!$OMP FLUSH
    CALL release_lock()

    IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
      WRITE(qmessage,'(A,A,A,I10,A,I5,A)')                                     &
          'IOS: Info: Queue: Added action ',                                   &
          TRIM(IOS_ActionName(node%metadata%action)),                          &
          ' trns_no: ',node%transaction_number ,                               &
          ' to queue, now ',Q_items,' items'
      WRITE(6,'(A)')TRIM(qmessage)
! DEPENDS ON: um_fort_flush
      CALL um_fort_flush(6,errCode)
    END IF

    IF ( lhook ) CALL dr_hook('IOS_QUEUE:PUTNODEINQUEUE',                      &
        zhook_out,zhook_handle)

  END SUBROUTINE IOS_Put_Node_In_Queue

!-----------------------------------------------------------------------
! SUBROUTINE to put metadata into Queue
!-----------------------------------------------------------------------
  SUBROUTINE IOS_Put_Metadata_In_Queue(metadata)
    IMPLICIT NONE
    TYPE(IOS_metadata_type), INTENT(IN) :: metadata
    TYPE(IOS_node_type), POINTER        :: node
    CHARACTER(LEN=80) mes
    node => make_new_node()
    node%metadata = metadata
    CALL IOS_Put_Node_In_Queue(node)
  END SUBROUTINE IOS_Put_Metadata_In_Queue

!-----------------------------------------------------------------------
! SUBROUTINE to remove last request and data from Queue
!-----------------------------------------------------------------------
  SUBROUTINE IOS_Remove_Data_From_Queue()
    IMPLICIT NONE
    TYPE(IOS_node_type), POINTER         :: node
    INTEGER                              :: errCode
    REAL(KIND=jprb)                      :: zhook_handle

    IF ( lhook ) CALL dr_hook('IOS_QUEUE:REMOVEDATAFROMQUEUE',                 &
        zhook_in,zhook_handle)

    CALL acquire_lock()
    IF ( ASSOCIATED (Q_first) ) THEN
      node => Q_first
      IF ( ASSOCIATED (Q_first%next ) ) THEN
        Q_first => Q_first%next
      ELSE! There is only 1 thing in the queue and we are
          ! about to delete it
        NULLIFY(Q_first)
        NULLIFY(Q_last)
      END IF
!$OMP FLUSH
      Q_items=Q_items-1
!$OMP FLUSH


      IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
        IF ( Q_items == 0 ) THEN
          WRITE(qmessage,'(A,A,A,I10,A)')                                      &
              'IOS: Info: Queue: Removing ',                                   &
              TRIM(IOS_ActionName(node%metadata%action)),                      &
              ' Trns no ',node%transaction_number,                             &
              ' from queue, queue is empty'
        ELSE
          WRITE(qmessage,'(A,A,A,I10,A,I5,A)')                                 &
              'IOS: Info: Queue: Removing ',                                   &
              TRIM(IOS_ActionName(node%metadata%action)),                      &
              ' Trns no. ',node%transaction_number,                            &
              ' from queue, now ',Q_items,' items'
        END IF
        WRITE(6,'(A)')TRIM(qmessage)
! DEPENDS ON: um_fort_flush
        CALL um_fort_flush(6,errCode)
      END IF
      CALL destroy_node( node )
    ELSE
      CALL IOS_Ereport( RoutineName,99,                                        &
          'Remove data from queue called on empty queue')
    END IF
    CALL IOS_Increment_Queue_Length(-1*ios_md_len*umFortranIntegerSize(),      &
        need_lock=.FALSE.)
    CALL release_lock()

    IF ( lhook ) CALL dr_hook('IOS_QUEUE:REMOVEDATAFROMQUEUE',                 &
        zhook_out,zhook_handle)

  END SUBROUTINE IOS_Remove_Data_From_Queue

!-----------------------------------------------------------------------
! SUBROUTINE to Get Node from Queue (from the front ;-) )
!-----------------------------------------------------------------------
  SUBROUTINE IOS_Get_Node_From_Queue(node)
    IMPLICIT NONE
    TYPE(IOS_node_type), POINTER  :: node

!$OMP FLUSH
    IF ( ASSOCIATED (Q_first) ) THEN
      node => Q_first
    ELSE
      CALL IOS_Ereport( RoutineName,99,                                        &
          'Tried to get a node from an empty queue')
    END IF
  END SUBROUTINE IOS_Get_Node_From_Queue

!-----------------------------------------------------------------------
! FUNCTION to size of the queue (this is the payload size in words)
!-----------------------------------------------------------------------
  INTEGER FUNCTION IOS_getQueuePayload()
    IMPLICIT NONE
    CALL acquire_lock()
!$OMP FLUSH
    IOS_getQueuePayload=Q_size
    CALL release_lock()
    RETURN
  END FUNCTION IOS_getQueuePayload

!-----------------------------------------------------------------------
! Subroutine to modify the size of the queue in bytes
!-----------------------------------------------------------------------
  SUBROUTINE IOS_Increment_Queue_Length(amount,need_lock)
    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: amount
    LOGICAL, INTENT(IN), OPTIONAL :: need_lock
    LOGICAL                       :: l_need_lock
    INTEGER                       :: errCode
    REAL(KIND=real64)             :: timeStamp
    REAL(KIND=jprb):: zhook_handle
    IF ( lhook ) CALL dr_hook('IOS_QUEUE:INCREMENT_QUEUE_LENGTH',              &
        zhook_in,zhook_handle)

    L_need_lock=.TRUE.
    IF ( PRESENT(need_lock) ) L_need_lock=need_lock

    IF ( l_need_lock ) CALL acquire_lock()
!$OMP FLUSH
    Q_size=Q_size+amount
!$OMP FLUSH
    IF ( Q_size > IOS_buffer_max ) THEN
      WRITE(qmessage,'(A,I12,A,I12,A)')'Adding ',amount,                       &
          ' to the queue exceeded IOS_buffer_max (',                           &
          IOS_buffer_max,')'
      CALL IOS_Ereport( RoutineName,-99, Qmessage )
    ELSE IF ( Q_size < 0 ) THEN
      WRITE(qmessage,'(A,I10,A,I10,A)')                                        &
          ' Decrementing the Q length by ',-1*amount,                          &
          'resulted in a negative length (',Q_Size,')'
      CALL IOS_Ereport( RoutineName,-99, Qmessage )
    END IF
! DEPENDS ON: get_wallclock_time
    timeStamp=get_wallclock_time()
!$OMP FLUSH
    WRITE(logunit,'(F8.3,A,I10)')                                              &
        timeStamp-IOS_start_time,' ',Q_size
!$OMP FLUSH
    IF ((1.0*Q_Size)/(1.0*IOS_buffer_max) > 0.9) THEN
      IF (IOS_Verbosity >= IOS_PrStatus_Oper) THEN
        WRITE(6,'(A,F8.3,A)')'IOS: WARNING: Queue:',                           &
            (100.0*Q_Size)/(1.0*IOS_buffer_max),'% capacity'
! DEPENDS ON: um_fort_flush
        CALL um_fort_flush(6,errCode)
      END IF
    END IF


    IF ( l_need_lock )CALL release_lock()
    IF ( lhook ) CALL dr_hook('IOS_QUEUE:INCREMENT_QUEUE_LENGTH',              &
        zhook_out,zhook_handle)

  END SUBROUTINE IOS_Increment_Queue_Length

!-----------------------------------------------------------------------
! Subroutines to guard against concurrent updates to the queue
!-----------------------------------------------------------------------
  SUBROUTINE acquire_lock()
    IMPLICIT NONE
    REAL :: t1,t2
! DEPENDS ON: get_wallclock_time
    IF (IOS_Lock_meter) t1=get_wallclock_time()
    CALL   omp_set_lock(Q_lock_var)
    IF (IOS_Lock_meter) THEN
      t2=get_wallclock_time()
      IOS_Q_LockTries(1,q_lock,omp_get_thread_num())=                          &
          IOS_Q_LockTries(1,q_lock,omp_get_thread_num())+1
      IOS_Q_LockTime(1,q_lock,omp_get_thread_num())=                           &
          IOS_Q_LockTime(1,q_lock,omp_get_thread_num())+(t2-t1)
    END IF
  END SUBROUTINE acquire_lock

  SUBROUTINE release_lock()
    IMPLICIT NONE
    CALL  omp_unset_lock(Q_lock_var)
  END SUBROUTINE release_lock

!----------------------------------------------------------------------
! Subroutines to guard against concurrent calls to MPI
! (which some MPI doesn't like)
!----------------------------------------------------------------------
  SUBROUTINE acquire_lock_MPI()
    IMPLICIT NONE
    REAL :: t1,t2

    IF ( IOS_thread_0_calls_mpi .AND. omp_get_thread_num() /= 0 ) THEN
      WRITE(qmessage,'(A,I2,A)')                                               &
          'Thread ',omp_get_thread_num(),                                      &
          ' tried to call mpi, but only thread 0 is allowed'
      CALL IOS_Ereport( RoutineName,99, Qmessage )
    END IF

    IF ( IOS_serialise_mpi_calls ) THEN
! DEPENDS ON: get_wallclock_time
      IF (IOS_Lock_meter) t1=get_wallclock_time()
      CALL omp_set_lock(IOS_comms_access_lock)
      IF (IOS_Lock_meter) THEN
        t2=get_wallclock_time()
        IOS_Q_LockTries (1,mpi_lock,omp_get_thread_num())=                     &
            IOS_Q_LockTries (1,mpi_lock,omp_get_thread_num())+1
        IOS_Q_LockTime  (1,mpi_lock,omp_get_thread_num())=                     &
            IOS_Q_LockTime  (1,mpi_lock,omp_get_thread_num())+(t2-t1)
      END IF
    END IF
  END SUBROUTINE acquire_lock_MPI

  SUBROUTINE release_lock_MPI()
    IMPLICIT NONE

    IF ( IOS_serialise_mpi_calls ) CALL  omp_unset_lock(IOS_comms_access_lock)

  END SUBROUTINE release_lock_MPI

!-----------------------------------------------------------------------
! Subroutine to initialize the queue
!-----------------------------------------------------------------------
  SUBROUTINE IOS_Queue_initialise(buffer)
    IMPLICIT NONE
    INTEGER, INTENT(IN)           :: buffer
    INTEGER                       :: errCode
    CHARACTER (LEN=20)            :: logname
    CHARACTER (LEN=80)            :: buffer_env
    REAL(KIND=jprb):: zhook_handle
    IF ( lhook ) CALL dr_hook('IOS_QUEUE:QUEUE_INITIALISE',                    &
        zhook_in,zhook_handle)

! Initialise variables
    logname = "dummy filename"
    NULLIFY(Q_First)
    NULLIFY(Q_Last )
    Q_size = 0
    Q_items= 0
    IOS_AsyncCompletionRequested=0
    CALL omp_init_lock( Q_lock_var)
    CALL omp_init_lock( IOS_comms_access_lock)

! Set the maxumum payload size of the queue
! Unit as given is in MB - convert to bytes
    IOS_buffer_max = buffer  * 1024 * 1024

! Open a file to log usage to
    WRITE(logname,'(A,I4.4)')'ioserver_log.',global_rank
    WRITE(6,'(A,A)')'IOS_queue: Logging to: ',TRIM(logname)
    OPEN(unit=logunit,file=TRIM(logname))

    IF (IOS_Lock_meter) THEN
      ALLOCATE(IOS_Q_LockTime                                                  &
          (words_per_cacheline,num_locks,0:omp_get_max_threads()-1))
      ALLOCATE(IOS_Q_LockTries                                                 &
          (words_per_cacheline,num_locks,0:omp_get_max_threads()-1))
      IOS_Q_LockTime (1,1:num_locks,0:omp_get_max_threads()-1)=0.0
      IOS_Q_LockTries(1,1:num_locks,0:omp_get_max_threads()-1)=0
    END IF

! Record the start time
    IF ( lhook ) CALL dr_hook('IOS_QUEUE:QUEUE_INITIALISE',                    &
        zhook_out,zhook_handle)

  END SUBROUTINE IOS_Queue_initialise

!-----------------------------------------------------------------------
! Subroutine to finalise the queue
!-----------------------------------------------------------------------
  SUBROUTINE IOS_Queue_Finalise()
    IMPLICIT NONE
    INTEGER :: thread

    CALL omp_destroy_lock( Q_lock_var)
    CALL omp_destroy_lock( IOS_comms_access_lock)


    IF (IOS_Lock_meter) THEN
      WRITE(6,'(A)')''
      WRITE(6,'(A)')'IOS Queue lock metering'
      WRITE(6,'(A)')'-----------------------'
      WRITE(6,'(A)')''
      WRITE(6,'(A)')'MPI Access:'
      WRITE(6,'(A6,A12,A12)')'Thread ','Locks','Time'
      WRITE(6,'(A30)')'------------------------------'
      DO thread=0,omp_get_max_threads()-1
        WRITE(6,'(I6,I12,F12.2)')thread,                                       &
            IOS_Q_LockTries(1,mpi_lock,thread),                                &
            IOS_Q_LockTime(1,mpi_lock,thread)
      END DO
      WRITE(6,'(A)')''
      WRITE(6,'(A)')'IOS Queue Access:'
      WRITE(6,'(A6,A12,A12)')'Thread ','Locks','Time'
      WRITE(6,'(A30)')'------------------------------'
      DO thread = 0,omp_get_max_threads()-1
        WRITE(6,'(I6,I12,F12.2)')thread,                                       &
            IOS_Q_LockTries(1,q_lock,thread),                                  &
            IOS_Q_LockTime(1,q_lock,thread)
      END DO
      WRITE(6,'(A)')''

    END IF

    WRITE(6,'(A)')'IOS: Info: Queue service terminated.'

  END SUBROUTINE IOS_Queue_Finalise

!-----------------------------------------------------------------------
! FUNCTION that reports whether the queue has space to accommodate another
! len words. It may wait until the condition is met, and it may not,
! hence the return code.
!-----------------------------------------------------------------------
  LOGICAL FUNCTION HasFreeSpace(amountNeeded)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: amountNeeded
    INTEGER             :: start_size
    INTEGER             :: start_len
    REAL                :: t1,t2
    REAL(KIND=jprb):: zhook_handle
    IF ( lhook ) CALL dr_hook('IOS_QUEUE:HASFREESPACE',                        &
        zhook_in,zhook_handle)

    HasFreeSpace=.TRUE.

!$OMP FLUSH
    IF ( amountNeeded  > IOS_buffer_max ) THEN
      WRITE(Qmessage,'(A,I10,A)')                                              &
          'The request for ',amountNeeded,                                     &
          ' words of space on the Q cannot be satisfied'
      CALL IOS_Ereport( RoutineName,99, Qmessage )
    END IF

    IF ( Q_size + amountNeeded  > IOS_buffer_max ) THEN
      start_size= Q_size
      start_len = Q_items
! DEPENDS ON: get_wallclock_time
      t1=get_wallclock_time()
!$OMP FLUSH
      IF (IOS_Verbosity>=IOS_PrStatus_Oper)                                    &
          WRITE(6,'(A,I10,A)')'IOS: Info: Queue: Waiting for ',                &
          amountNeeded,' units of free space'
      DO WHILE(Q_size + amountNeeded  > IOS_buffer_max .AND.                   &
          IOS_AsyncCompletionRequested==0)
        CALL um_sleep(IOS_backoff_interval)
!$OMP FLUSH
      END DO
! DEPENDS ON: get_wallclock_time
      t2=get_wallclock_time()
      IF (IOS_Verbosity>=IOS_PrStatus_Oper)                                    &
          WRITE(6,'(A,F8.2,A)')'IOS: Info: Queue: Wait done: ',                &
          t2-t1,'s'
      IF ( IOS_AsyncCompletionRequested /= 0 ) HasFreeSpace=.FALSE.
    END IF
    IF ( lhook ) CALL dr_hook('IOS_QUEUE:HASFREESPACE',                        &
        zhook_out,zhook_handle)
  END FUNCTION HasFreeSpace

!-----------------------------------------------------------------------
! FUNCTION that reports whether the queue is fully drained
! It may wait until the condition is met, and it may not.
!-----------------------------------------------------------------------
  LOGICAL FUNCTION WaitForDrain()
    IMPLICIT NONE
    INTEGER            :: counter
    INTEGER, PARAMETER :: counter_targ = 300
    REAL               :: t1
    REAL               :: t2
    REAL(KIND=jprb)    :: zhook_handle

    IF ( lhook ) CALL dr_hook('IOS_QUEUE:WAITFORDRAIN',                        &
        zhook_in,zhook_handle)
    WaitForDrain=.TRUE.
! DEPENDS ON: get_wallclock_time
    t1=get_wallclock_time()
!$OMP FLUSH
    IF ( Q_items > 0 ) THEN
      IF (IOS_Verbosity>=IOS_PrStatus_Oper)                                    &
          WRITE(6,'(A)')'IOS: Info: Queue: Waiting for drain'
      DO WHILE(Q_items > 0 .AND. IOS_AsyncCompletionRequested==0)
        CALL um_sleep(IOS_backoff_interval)
        counter=counter+1
        IF (counter==counter_targ) THEN
! DEPENDS ON: get_wallclock_time
          counter = 0
          t2=get_wallclock_time()
          IF (IOS_Verbosity>=IOS_PrStatus_Diag)                                &
              WRITE(6,'(A,F8.2,A)')'IOS: Info: Queue: Waiting: ',              &
              t2-t1,'s'
        END IF
!$OMP FLUSH
      END DO
      IF ( IOS_AsyncCompletionRequested /= 0 ) WaitForDrain=.FALSE.
      IF (IOS_Verbosity>=IOS_PrStatus_Oper)                                    &
          WRITE(6,'(A)')'IOS: Info: Queue: done drain'
    END IF
    IF ( lhook ) CALL dr_hook('IOS_QUEUE:WAITFORDRAIN',                        &
        zhook_out,zhook_handle)
  END FUNCTION WaitForDrain

!-----------------------------------------------------------------------
! FUNCTION that reports the number of distinct items in the queue
!-----------------------------------------------------------------------
  INTEGER FUNCTION IOS_getQueueItems()
    IMPLICIT NONE
!$OMP FLUSH
    IOS_getQueueItems=Q_items
    RETURN
  END FUNCTION IOS_getQueueItems

!-----------------------------------------------------------------------
! FUNCTION: A factory that makes nodes to put in the queue
!-----------------------------------------------------------------------
  FUNCTION make_new_node() result(n)

    IMPLICIT NONE
    TYPE(IOS_node_type),POINTER :: n

    ALLOCATE(n)

    NULLIFY(n%next)
    NULLIFY(n%prev)
    NULLIFY(n%real_data)
    NULLIFY(n%real32_data)
    NULLIFY(n%distributed_data)
    NULLIFY(n%integer_data)
    NULLIFY(n%integer32_data)
    NULLIFY(n%receive_requests)
    NULLIFY(n%receive_tags)
    NULLIFY(n%receive_data_len)

    n%metadata%action=IOS_Action_unset
    n%transaction_number=-1

  END FUNCTION make_new_node

!-----------------------------------------------------------------------
! SUBROUTINE: A scrapyard for recyling old nodes safely
!-----------------------------------------------------------------------
  SUBROUTINE destroy_node(node)
    IMPLICIT NONE
    TYPE(IOS_node_type), POINTER :: node
    INTEGER                      :: wordsRemoved

    REAL(KIND=jprb):: zhook_handle
    IF ( lhook ) CALL dr_hook('IOS_QUEUE:DESTROY_NODE',                        &
        zhook_in,zhook_handle)
    wordsRemoved=0

    IF ( ASSOCIATED(node%real_data) ) THEN
      DEALLOCATE  (node%real_data)
      NULLIFY     (node%real_data)
      wordsRemoved=wordsRemoved-node%metadata%data_size*                       &
          IOS_BytesPerWord64
    END IF

    IF ( ASSOCIATED(node%integer_data) ) THEN
      DEALLOCATE  (node%integer_data)
      NULLIFY     (node%integer_data)
      wordsRemoved=wordsRemoved-node%metadata%data_size*                       &
          IOS_BytesPerInteger
    END IF

    IF ( ASSOCIATED(node%real32_data) ) THEN
      DEALLOCATE  (node%real32_data)
      NULLIFY     (node%real32_data)
      wordsRemoved=wordsRemoved-node%metadata%data_size*                       &
          IOS_BytesPerWord32
    END IF

    IF ( ASSOCIATED(node%integer32_data) ) THEN
      DEALLOCATE  (node%integer32_data)
      NULLIFY     (node%integer32_data)
      wordsRemoved=wordsRemoved-node%metadata%data_size*                       &
          IOS_BytesPerWord32
    END IF

    IF ( ASSOCIATED(node%receive_requests) ) THEN
      DEALLOCATE  (node%receive_requests)
      NULLIFY     (node%receive_requests)
    END IF

    IF ( ASSOCIATED(node%receive_tags) ) THEN
      DEALLOCATE (node%receive_tags)
      NULLIFY    (node%receive_tags)
    END IF

    IF ( ASSOCIATED(node%receive_data_len) ) THEN
      DEALLOCATE  (node%receive_data_len)
      NULLIFY     (node%receive_data_len)
    END IF

    IF ( ASSOCIATED(node%distributed_data) ) THEN
      DEALLOCATE  (node%distributed_data)
      NULLIFY     (node%distributed_data)
    END IF

!$OMP FLUSH
    CALL IOS_Increment_Queue_Length(wordsRemoved,need_lock=.FALSE.)
!$OMP FLUSH

    DEALLOCATE(node)
    NULLIFY(node)
    IF ( lhook ) CALL dr_hook('IOS_QUEUE:DESTROY_NODE',                        &
        zhook_out,zhook_handle)
  END SUBROUTINE destroy_node

!-----------------------------------------------------------------------
! SUBROUTINE: A means of dumping the state of the queue for debug
!-----------------------------------------------------------------------
  SUBROUTINE IOS_Queue_Report
    IMPLICIT NONE

!$OMP FLUSH
    CALL acquire_lock()
    WRITE(6,'(A,I6,A,I6,A,F5.2,A)')                                            &
        'IOS Q: ',                                                             &
        Q_size*8/1024.0/1024.0,' MB in ',                                      &
        Q_items,                                                               &
        ' items (',                                                            &
        (Q_size*100.0)/IOS_buffer_max,'%)'
    IF ( Q_items == 0 .AND. Q_size /= 0 ) THEN
      WRITE(qmessage,'(A)')'Q_items and Q_size mismatch'
      CALL IOS_Ereport( RoutineName,-99, Qmessage )
    END IF
    IF      ( Q_items == 0 ) THEN
      WRITE(6,'(A)')'Q contains no transactions '
    ELSE IF ( Q_items == 1 ) THEN
      WRITE(6,'(A,I8)')'Q contains transaction ',                              &
          Q_first%transaction_number
    ELSE IF ( Q_items >= 1 ) THEN
      WRITE(6,'(A,I8,A,I8)')'Q contains transactions ',                        &
          Q_first%transaction_number,                                          &
          ' to ',Q_last%transaction_number
    END IF
    CALL release_lock()
  END SUBROUTINE IOS_Queue_Report

END MODULE IOS_Queue_Mod
