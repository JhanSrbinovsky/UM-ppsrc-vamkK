! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: C96

!---------------------------------------------------------------------
! Initialisation of IOS, called whether or not IOS to be used
!---------------------------------------------------------------------

MODULE IOS_Init
  USE UM_Types
  USE IOS_Common
  ! C98_1A means OpenMP Compiles Only

  USE IOS_MPI_error_handlers
  USE missing_data_mod, ONLY : imdi
  USE ereport_mod, ONLY : ereport

  IMPLICIT NONE
  INTEGER, PRIVATE                :: IOS_spacing
  INTEGER, PRIVATE                :: IOS_buffer_size
  INTEGER, PRIVATE                :: IOS_Offset
  INTEGER, PRIVATE                :: IOS_as_concurrency
  INTEGER, PRIVATE                :: IOS_async_levs_per_pack
  INTEGER, PRIVATE                :: IOS_force_threading_mode
! where:
!  0 = single
!  1 = funneled
!  2 = serialized
!  3 = multiple
  LOGICAL, PRIVATE                :: IOS_Interleave
  LOGICAL, PRIVATE                :: IOS_use_async_stash
  LOGICAL, PRIVATE                :: IOS_use_async_dump
  LOGICAL, PRIVATE                :: IOS_debug_no_write
  LOGICAL, PRIVATE                :: IOS_debug_no_packing
  LOGICAL, PRIVATE                :: IOS_debug_no_subdomaining
  LOGICAL, PRIVATE                :: IOS_async_send_null
  LOGICAL, PRIVATE                :: IOS_async_stats
  LOGICAL, PRIVATE                :: IOS_use_helpers

  PRIVATE ioscntl
  NAMELIST / ioscntl /                                                         &
      IOS_tasks_per_server,                                                    &
      IOS_Spacing,                                                             &
      IOS_Interleave,                                                          &
      IOS_Offset,                                                              &
      IOS_buffer_size,                                                         &
      IOS_use_async_stash,                                                     &
      IOS_use_async_dump,                                                      &
      IOS_serialise_mpi_calls,                                                 &
      IOS_thread_0_calls_mpi,                                                  &
      IOS_force_threading_mode,                                                &
      IOS_debug_no_packing,                                                    &
      IOS_debug_no_write,                                                      &
      IOS_debug_no_subdomaining,                                               &
      IOS_Verbosity,                                                           &
      IOS_backoff_interval,                                                    &
      IOS_timeout,                                                             &
      IOS_acquire_model_prsts,                                                 &
      IOS_local_ro_files,                                                      &
      IOS_concurrency,                                                         &
      IOS_as_concurrency,                                                      &
      IOS_Unit_Alloc_Policy,                                                   &
      IOS_async_levs_per_pack,                                                 &
      IOS_async_send_null,                                                     &
      IOS_async_stats,                                                         &
      IOS_use_helpers,                                                         &
      IOS_Lock_Meter

CONTAINS

  LOGICAL FUNCTION IOS_Setup( numIOServers )
    USE mpl, ONLY :                                                            &
        mpl_comm_null,                                                         &
        mpl_thread_multiple,                                                   &
        mpl_thread_serialized,                                                 &
        mpl_thread_funneled,                                                   &
        mpl_thread_single
    USE FilenameLength_mod, ONLY :                                             &
        filenamelength
    USE model_file, ONLY :                                                     &
        Mf_unit_min,                                                           &
        Mf_unit_max

    IMPLICIT NONE

    ! Argument
    INTEGER, INTENT(IN)    :: numIOServers

    INTEGER, PARAMETER     :: atm_coloured=0
    INTEGER, PARAMETER     :: io_coloured=1
    INTEGER, PARAMETER     :: IOS_Namelist_Unit=127

! Note, formatting of output copes with 1,000,000 mpi ranks.
    INTEGER, PARAMETER     :: IOS_Max_Reported=8

    INTEGER                :: io_server_counter
    INTEGER                :: io_server_rank

    INTEGER                :: handler
    INTEGER                :: ierror
    INTEGER                :: num_threads
    INTEGER                :: listener_thread=-1
    INTEGER                :: writer_thread=-1
    INTEGER                :: reader_thread=-1
    INTEGER                :: this_thread 
    LOGICAL                :: isIOSCapable

    ! Various local vars for arithmetic
    INTEGER                :: atm_next_key
    INTEGER                :: num_units
    INTEGER                :: remaining_ios
    INTEGER                :: allocated_ios
    INTEGER                :: colour
    INTEGER                :: key
    INTEGER                :: subkey
    INTEGER                :: subcolour
    INTEGER                :: i,j
    LOGICAL                :: Flag

    ! Strings for IO
    CHARACTER (LEN=80)             :: IOS_ini_message
    CHARACTER (LEN=80)             :: IOS_capable_string
    CHARACTER (LEN=FileNameLength) :: NamelistFile
    CHARACTER (LEN=*), PARAMETER   :: IOS_Setup_Name = 'IOS_Setup'

    REAL, EXTERNAL                 :: get_wallclock_time

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

! Get current communicator from GCOM
    CALL GC_Get_Communicator(global_comm, ierror)
    CALL MPL_Comm_Rank(global_comm,global_rank  ,ierror)
    CALL MPL_Comm_Size(global_comm,global_procs ,ierror)

! Set some sensible defaults and initialise variables
! All namelist items not listed here are initialised in IOS_common
    IOS_Setup                   = .FALSE.
    IOS_force_threading_mode    = IMDI
    IOS_ReadCompletionRequested = 0
    IOS_EnqCompletionRequested  = 0
    IOS_spacing                 = IMDI
    IOS_Offset                  = IMDI
    IOS_Interleave              = .FALSE.
    IOS_use_async_stash         = .FALSE.
    IOS_use_async_dump          = .FALSE.
    l_io_server                 = .FALSE.
    l_io_leader                 = .FALSE.
    IOS_debug_no_packing        = .FALSE.
    IOS_debug_no_write          = .FALSE.
    IOS_debug_no_subdomaining   = .FALSE.
    IOS_async_send_null         = .FALSE.
    IOS_async_stats             = .FALSE.
    IOS_buffer_size             = IMDI
    IOS_concurrency             = IMDI
    IOS_as_concurrency          = IMDI
    IOS_async_levs_per_pack     = IMDI


! Read Namelist for control
    CALL fort_get_env(ft_environ(IOS_Namelist_Unit),                           &
        len_ft_envir(IOS_Namelist_Unit),                                       &
        NameListFile,FileNameLength,ierror)
    IF ( ierror /= 0 ) THEN
      WRITE(6,'(A,I4,A,A)')                                                    &
          'IOS: Info: Failed to get setup file unit=',                         &
          IOS_Namelist_Unit,' Name=',ft_environ(IOS_Namelist_Unit)

      NameListFile='IOSCNTL.namelist'
    END IF
    WRITE(6,'(A,A)')'IOS: Info: Control file:',TRIM(NameListFile)
    OPEN(IOS_Namelist_Unit,file=NameListFile,form='FORMATTED',                 &
        status='OLD',iostat=ierror)
    IF ( ierror == 0 ) THEN
      WRITE(6,'(A)')'IOS: Info: Reading IOS control file'
      READ(IOS_Namelist_Unit,ioscntl,iostat=ierror)
      IF ( ierror /= 0 ) THEN
        WRITE(6,'(A,A)')                                                       &
            'IOS: Warning: Error whilst reading IOS control file ',            &
            '(please check setup).'
      END IF

      CLOSE(IOS_Namelist_Unit)
    ELSE
      WRITE(6,'(A,A,A)')'IOS: Info: Failed to open IOS control file ',         &
          TRIM(NameListFile),' (please check setup).'
    END IF

! Some checks and resets:
    IF (numIOServers <= 0) THEN

! No I/O servers
      IOS_spacing = global_procs

    ELSE

      IF (IOS_Spacing == 0) THEN
! User has requested evenly distributed IO servers
        IOS_Spacing = global_procs/numIOServers
        IOS_Offset  = 0
      ELSE IF (IOS_Spacing < 0) THEN
! Spacing not set or invalid input
        WRITE(IOS_ini_message,'(A)')                                   &
          'Spacing of IO servers is negative or undefined'
        CALL IOS_Ereport(IOS_Setup_Name,10,IOS_ini_message)
      END IF

    END IF


! Set namelist flags in other modules
!$  IOS_AsyncNumSlots        = IOS_as_concurrency
!$  IOS_AsyncMaxFieldsInPack = IOS_async_levs_per_pack
!$  IOS_AsyncSendNull        = IOS_async_send_null
!$  IOS_AsyncDoStats         = IOS_async_stats

! Set default word lengths for types
    IOS_BytesPerInteger      = umFortranIntegerSize()
    IOS_BytesPerReal         = umFortranRealSize()

! Determine capability for IOS of runtime environment
    isIOSCapable=.FALSE.
    WRITE(IOS_capable_string,'(A)')'Not OpenMP compiled'
!$  isIOSCapable=.TRUE.  ! If openMP compiled
!$  WRITE(IOS_capable_string,'(A)')'IOS is possible'

! Override settings according to MPI's capabilities
    CALL mpl_query_thread (threading_model,ierror)

    IF (ios_force_threading_mode /= IMDI) THEN
      CALL Ereport("IOS_INIT", -10,                                    &
        "MPI's advertised threading capability is being overriden.") 
    END IF

! Check for a namelist override
    IF      (ios_force_threading_mode==0) THEN
      threading_model=MPL_THREAD_SINGLE
    ELSE IF (ios_force_threading_mode==1) THEN
      threading_model=MPL_THREAD_FUNNELED
    ELSE IF (ios_force_threading_mode==2) THEN
      threading_model=MPL_THREAD_SERIALIZED
    ELSE IF (ios_force_threading_mode==3) THEN
      threading_model=MPL_THREAD_MULTIPLE
    ELSE IF (ios_force_threading_mode/=IMDI) THEN
      CALL Ereport("IOS_INIT", 10,                                    &
        "Incorrect value set for ios_force_threading_mode.")
    END IF ! Any other value retains the result of MPL_QUERY_THREAD()

    IF      ( threading_model == mpl_thread_multiple ) THEN
      WRITE(6,'(A)')'IOS: Info: Full Multithreading available.'
    ELSE IF ( threading_model == mpl_thread_serialized ) THEN
      IOS_serialise_mpi_calls=.TRUE.
      WRITE(6,'(A)')'IOS: Info: Serialized threading available.'
    ELSE IF ( threading_model == mpl_thread_funneled ) THEN
      IOS_thread_0_calls_mpi=.TRUE.
      WRITE(6,'(A)')'IOS: Info: Funneled threading available.'
    ELSE IF ( threading_model == mpl_thread_single ) THEN
      isIOSCapable=.FALSE.
      WRITE(IOS_capable_string,'(A,A)')                                        &
          'MPI implementation not capable ',                                   &
          'of running threaded jobs.'
      WRITE(6,'(A,A)')  &
      'If you are using OpenMPI or think your MPI is thread capable, ', &
      'try setting ios_force_threading_mode appropriately in the UI' 
    ELSE
      isIOSCapable=.FALSE.
      WRITE(IOS_capable_string,'(A)')'Broken MPI cannot describe itself'
      WRITE(IOS_ini_message,'(A,I8)')                                          &
          'MPI reported an unknown threading model:',threading_model
      CALL IOS_Ereport(IOS_Setup_Name,-99,IOS_ini_message)
    END IF

! Saneness for single threaded MPI
    IF (IOS_thread_0_calls_mpi  .AND.    &
        (IOS_use_async_stash.OR.IOS_use_async_dump)) THEN
      IOS_use_async_stash=.FALSE.
      IOS_use_async_dump=.FALSE.
      WRITE(6,'(A,A)')'IOS: Info: IOS Accelerated dump/STASH disabled:',       &
          ' Require at least MPL_THREAD_SERIALIZED'
    END IF

! Print out parameters being used
    IF ( global_rank == 0 .AND.                                                &
        IOS_Verbosity >= IOS_PrStatus_Oper) THEN
      WRITE(6,'(A,I5,A)')'IOS: Info: Task spacing               = ',           &
          IOS_Spacing, ' tasks'
      WRITE(6,'(A,I5,A)')'IOS: Info: Task offset                = ',           &
          IOS_Offset, ' tasks'
      WRITE(6,'(A,L1,A)')'IOS: Info: Interleaved servers        = ',           &
          IOS_Interleave
      WRITE(6,'(A,I5,A)')'IOS: Info: Server size                = ',           &
          IOS_tasks_per_server, ' tasks'
      WRITE(6,'(A,I5,A)')'IOS: Info: Buffer size                = ',&
          IOS_buffer_size,' MB'
      WRITE(6,'(A,I1)')'IOS: Info: Verbosity                  = ',&
          IOS_Verbosity
      WRITE(6,'(A,L1)')'IOS: Info: Asynchronous stash         = ',             &
          IOS_Use_Async_Stash
      WRITE(6,'(A,L1)')'IOS: Info: Asynchronous dumps         = ',             &
          IOS_Use_Async_Dump
      WRITE(6,'(A,L1)')'IOS: Info: Serialise all mpi calls    = ',&
          IOS_serialise_mpi_calls
      WRITE(6,'(A,L1)')'IOS: Info: Only thread zero calls MPI = ',&
          IOS_thread_0_calls_mpi 
      WRITE(6,'(A,L1)')'IOS: Info: Read only files stay local = ',&
          IOS_local_ro_files
      WRITE(6,'(A,I1)')'IOS: Info: Unit allocation policy     = ',             &
          IOS_Unit_Alloc_Policy
      WRITE(6,'(A,I5)')'IOS: Info: Polling Interval           = ',             &
          IOS_backoff_interval
      WRITE(6,'(A,I3)')'IOS: Info: Timeout Interval           = ',             &
          IOS_timeout
      WRITE(6,'(A,L1)')'IOS: Info: Lock Metering              = ',             &
          IOS_lock_meter
      WRITE(6,'(A,L1)')'IOS: Info: Acquire model output level = ',             &
          IOS_acquire_model_prsts
      WRITE(6,'(A,L1)')'IOS: Info: Use helper threads         = ',             &
          IOS_use_helpers
      WRITE(6,'(A,I3)')'IOS: Info: Async Dispatch slots       = ',             &
          IOS_Concurrency
!$    WRITE(6,'(A,I3)')'IOS: Info: Async Stash Dispatch slots = ',             &
!$        IOS_AsyncNumSlots
!$    WRITE(6,'(A,I3)')'IOS: Info: Async fields/levs per pack = ',             &
!$        IOS_AsyncMaxFieldsInPack
!$    WRITE(6,'(A,L1)')'IOS: Info: Async send empty tiles     = ',             &
!$        IOS_AsyncSendNull
!$    WRITE(6,'(A,L1)')'IOS: Info: Async stats profiling      = ',             &
!$        IOS_AsyncDoStats
      IF (IOS_debug_no_packing)                                                &       
          WRITE(6,'(A)')'IOS: Info: Debug Option: No stash packing.'      
      IF (IOS_debug_no_subdomaining)                                           &
          WRITE(6,'(A)')'IOS: Info: Debug Option: No stash subdomaining.'
      IF (IOS_debug_no_write)                                                  &
          WRITE(6,'(A)')'IOS: Info: Debug Option: No data writes.'
          
    END IF

    IF(IOS_debug_no_packing .OR. IOS_debug_no_subdomaining .OR. &
        IOS_debug_no_write) THEN
      CALL IOS_Ereport(IOS_Setup_Name,-99, &
          "One or more debug options are set, output may be invalid")
    END IF

    IOS_server_groups=numIOServers/IOS_tasks_per_server

    IF (IOS_server_groups*IOS_tasks_per_server                                 &
        /= numIOServers) THEN
      WRITE(ios_ini_message,'(A,I3,A,I3,A)')                                   &
          'The number of IO tasks (',numIOServers,                             &
          ')must be a multiple of server size (',IOS_tasks_per_server,')'
      CALL IOS_Ereport(IOS_Setup_Name,99,IOS_ini_message)
    END IF

! Allocate Storage
    ALLOCATE(io_servers                                                        &
        (IOS_server_groups,IOS_tasks_per_Server))

! Allocate global processor numbers to IOS or ATM models
    atm_next_key      = 0
    io_server_counter = 0
    io_server_rank    = 0
    allocated_ios     = 0

! Loop over cpus and allocate them to a role
    DO i=0,global_procs-1

      remaining_ios=numIOServers-allocated_ios
      IF (( i-IOS_Offset >= 0                            .AND.                 &
          MOD(i-IOS_Offset,IOS_Spacing) == IOS_Spacing-1 .AND.                 &
          allocated_ios    < numIOServers )              .OR.                  &
          global_procs-i-1 < remaining_ios ) THEN

        IF ( global_rank == i ) THEN
          colour = io_coloured+io_server_counter
          key    = io_server_rank
        END IF

        io_servers(io_server_counter+1,io_server_rank+1)=i
        allocated_ios = allocated_ios+1

        IF (IOS_Interleave) THEN
          io_server_counter=io_server_counter+1
          IF (io_server_counter == IOS_Server_groups ) THEN
            io_server_counter = 0
            io_server_rank    = io_server_rank+1
          END IF
        ELSE
          io_server_rank = io_server_rank+1
          IF (io_server_rank == IOS_tasks_per_server) THEN
            io_server_rank    = 0
            io_server_counter = io_server_counter+1
          END IF
        END IF

      ELSE ! not an io server
        IF ( global_rank == i ) THEN
          colour     = atm_coloured
          key        = atm_next_key
        END IF
        atm_next_key = atm_next_key+1
      END IF
    END DO

    IF ( numIOServers > 0 ) THEN
      WRITE(6,'(A)')' '
      WRITE(ios_ini_message,*)                                                 &
          '(A,I5,A',(',I6',j=1,MIN(IOS_Max_Reported,IOS_tasks_per_server)),')'
      DO i=1,IOS_Server_groups
        ! This string is 30 + 6* ios_tasks_per_server
        ! so it can only support (30 is the number of literal chars)
        ! (80-30)/6 = 8 = IOS_Max_Reported servers
        WRITE(6,ios_ini_message)'IOS: Info: IO Server ',i,' is ',              &
            (io_servers(i,j),j=1,MIN(IOS_Max_Reported,IOS_tasks_per_server))
        IF (IOS_tasks_per_server > IOS_Max_Reported) THEN
          WRITE(6,'(A,I5,A)')'         :    and an additional',                &
              IOS_tasks_per_server-IOS_Max_Reported,' other processes '
        END IF
      END DO
    ELSE
      WRITE(6,'(A)')'IOS: Info: IO servers are not configured'
    END IF

!
! Set up all the communicators we want....
!
!
! Split the communicator into IO and ATM tasks
    subcolour=atm_coloured
    subkey=key
    IF (colour>atm_coloured) THEN
      subcolour = io_coloured
      subkey    = key+colour*IOS_tasks_per_server
    END IF
    CALL MPL_Comm_Split(global_comm, subcolour, subkey, io_comm, ierror)
    CALL MPL_Comm_Rank (io_comm, io_rank ,ierror)
    CALL MPL_Comm_Size (io_comm, io_procs ,ierror)

! IO communicator is useless on atmos ranks so zap it.
    IF (subcolour == atm_coloured) THEN
      io_comm  = mpl_comm_null
      io_procs = numIOServers
      io_rank  = -1
    END IF

! Split the communicator into leaders and non-leaders
    IF (key==0) THEN
      subcolour=0
      subkey=colour
    ELSE
      subcolour=1
      subkey=0
    END IF
    CALL MPL_Comm_Split(global_comm, subcolour, subkey, leader_comm,           &
        ierror)
    CALL MPL_Comm_Rank (leader_comm, leader_rank ,ierror)
    CALL MPL_Comm_Size (leader_comm, leader_procs ,ierror)
! leader communicator is useless on non-leader ranks so zap it.
    IF (subcolour == 1) THEN
      leader_comm  = mpl_comm_null
      leader_procs = IOS_Server_groups+1
      leader_rank  = -1
    END IF

! Split the communicator into the main IO Groups and ATM tasks
    CALL MPL_Comm_Split(global_comm, colour, key, model_comm, ierror)
    CALL MPL_Comm_Rank (model_comm , model_rank  ,ierror)
    CALL MPL_Comm_Size (model_comm , model_procs ,ierror)

! Lets keep funky communications in a funky communicator....
!$  CALL MPL_Comm_dup(global_comm,IOS_async_comm,ierror)

! Set up communicators for broadcasts....
    IF ( colour > atm_coloured ) THEN
      l_io_server = .TRUE.
      IF (leader_comm /= mpl_comm_null) THEN
        l_io_leader = .TRUE.
      END IF
    END IF

    ALLOCATE(IOS_BCast_Comm(numIOServers))
! Generate an MPI group suitable for broadcasting from an IO server

! set the server comm to something sane in case of inadvertant use.
    IOS_BCast_Server_Comm=mpl_comm_null

! There is a broadcast comm for each io server group.
    DO i=1,IOS_Server_groups
      IF (l_io_leader .AND. colour==i) THEN
        subcolour=1
        subkey=global_procs-io_procs
      ELSE IF (l_io_server) THEN
        subcolour=2
        subkey=model_rank
      ELSE
        subcolour=1
        subkey=model_rank
      END IF

      CALL MPL_Comm_Split(global_comm, subcolour,                              &
          subkey, IOS_bcast_comm(i), ierror)
      IF (subcolour /= 2) THEN
        CALL MPL_Comm_Size( IOS_bcast_comm(i), bcast_procs ,ierror)
        CALL MPL_Comm_Rank( IOS_bcast_comm(i), bcast_rank  ,ierror)
      END IF

      IF (l_io_leader .AND. colour == i) THEN
        IOS_BCast_Server_Comm=IOS_bcast_comm(i)
      END IF

    END DO






! Everyone synchronise their watches...
! DEPENDS ON: get_wallclock_time
    IOS_Start_Time=get_wallclock_time()

! Set the default communicator for GCOM
    CALL GC_Set_Communicator                                                   &
        (model_comm, model_rank, model_procs, ierror)

    IF (IOS_Verbosity >= IOS_PrStatus_Oper) THEN

      WRITE(6,'(A,I7,A,I7,A,I7)')                                              &
          'IOS: Info: Original Rank=',global_rank,' Global_Rank=',             &
          global_rank,' Model_Rank=',model_rank
      WRITE(6,'(A,I7,A,I7)')                                                   &
          'IOS: Info: Total Size=',global_procs,' Submodel size=',             &
          model_procs

      WRITE(6,'(A)')' '
      WRITE(6,'(A)')'Communicators....'
      WRITE(6,'(A)')' '
      WRITE(6,'(A12,A12,A12,A12)')'Function',                                  &
          'Comm ID','Processors','My Rank'
      WRITE(6,'(A12,A12,A12,A12)')'--------',                                  &
          '-------','----------','-------'
      WRITE(6,'(A12,I12,I12,I12)')'Global',                                    &
          global_comm,global_procs,global_rank
      WRITE(6,'(A12,I12,I12,I12)')'Model' ,                                    &
          model_comm,model_procs,model_rank
      WRITE(6,'(A12,I12,I12,I12)')'IO'    ,                                    &
          io_comm,io_procs,io_rank
      WRITE(6,'(A12,I12,I12,I12)')'Leader',                                    &
          leader_comm,leader_procs,leader_rank
      WRITE(6,'(A)')' '

    END IF

! The primary return of the function
    IOS_Setup = l_io_server

! register error handlers:

    CALL MPL_Comm_Create_ErrHandler(model_mpi_error_handler,                   &
        handler,ierror)
    CALL MPL_Comm_Set_ErrHandler(model_comm,handler,ierror)

!$  CALL MPL_Comm_Create_ErrHandler&
!$      (async_mpi_error_handler,handler,ierror)
!$  CALL MPL_Comm_Set_ErrHandler(IOS_async_comm,handler,ierror)

    CALL MPL_Comm_Create_ErrHandler(global_mpi_error_handler,                  &
        handler,ierror)
    CALL MPL_Comm_Set_ErrHandler(global_comm,handler,ierror)

    IF ( colour > atm_coloured ) THEN
      CALL MPL_Comm_Create_ErrHandler(io_mpi_error_handler,                    &
          handler,ierror)
      CALL MPL_Comm_Set_ErrHandler(io_comm,handler,ierror)
    END IF

    IF ( leader_rank >= 0 ) THEN
      CALL MPL_Comm_Create_ErrHandler(leader_mpi_error_handler,                &
          handler,ierror)
      CALL MPL_Comm_Set_ErrHandler(leader_comm,handler,ierror)
    END IF

! Check for no servers defined
    IF ( numIOServers <= 0 ) THEN
      isIOSCapable=.FALSE.
      WRITE(IOS_capable_string,'(A)')'No IO server processes assigned'
    END IF
! Deactivate if insufficient threads
    num_threads=2 ! allow for the 'not openmp message' to propogate
!$OMP PARALLEL DEFAULT(NONE) SHARED(num_threads)
!$OMP MASTER
!$    num_threads=omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
    IF ( num_threads < 2 .AND. isIOSCapable) THEN
      isIOSCapable=.FALSE.
      WRITE(IOS_capable_string,'(A,I3)')'Not enough OpenMP Threads:',          &
          num_threads
    END IF

! At present only 3 threads do anything, so we'll space them out
! over available threads...
    Listener_thread =0
    Writer_thread   =num_threads/2

    IF (num_threads>2 .AND. IOS_use_helpers)                                   &
        reader_thread = num_threads-1


    IF ( isIOSCapable ) THEN
! The following code depends on modules only available in OpenMP
! compilation so it should not compile on non-threaded builds (!$)
!
! --- Define the unit mapping - which IOS pes deal with which units ---
!
!$    IF (IOS_Unit_Alloc_Policy /= IOS_Unit_Alloc_Static) THEN
!$      io_server_for_unit_lookup(:)=IOS_No_Server
!$    ELSE
!$      DO i = IOS_unit_first, IOS_unit_last
!$        IF      ( i < Mf_unit_min ) THEN
!$          io_server_for_unit_lookup(i) = io_servers(1,1)
!$        ELSE IF ( i > Mf_unit_max ) THEN
!$          io_server_for_unit_lookup(i) = io_servers(IOS_Server_Groups,1)
!$        ELSE
!$          io_server_for_unit_lookup(i) =                                     &
!$              io_servers(1+MOD(i - Mf_unit_min, IOS_Server_groups),1)
!$        END IF
!$      END DO
!$    END IF
!$
!$    IF ( l_io_server ) THEN
!$      WRITE(6,'(A,I7,A)')                                                    &
!$          'IOS: Info: PE ', global_rank, ' is an I/O server'
!$      num_units=0
!$      IF ( l_io_leader ) THEN
!$        DO i = IOS_unit_first, IOS_unit_last
!$          IF ( io_server_for_unit(i) == global_rank ) THEN
!$              num_units=num_units+1
!$          END IF
!$        END DO
!$
!$        WRITE(6,'(A,I3,A)')                                                  &
!$            'IOS: Info: I am responsible for ',num_units,' units'
!$      END IF
!
! --- Initialise IOS queue structure ---
!
!$      CALL IOS_Queue_Initialise(IOS_buffer_size)
!
! --- Initialise IOS Stash support (server side) ---
! ---  this is incomplete initialisation, geometry will
! ---  need initialising later after the model is loaded
!
!$      CALL IOS_stash_server_init(                                            &
!$          IOS_use_async_stash,                                               &
!$          IOS_use_async_dump,                                                &
!$          IOS_debug_no_packing,                                              &
!$          IOS_debug_no_write,                                                &
!$          IOS_debug_no_subdomaining)
!


! --- Allocate subtasks to server threads ---
!
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(this_thread) SHARED(writer_thread,listener_thread,reader_thread)
!$      this_thread=omp_get_thread_num()       
!$      IF      ( this_thread == listener_thread ) THEN
!$        CALL IOS_listener()
!$      ELSE IF ( this_thread == writer_thread ) THEN
!$        CALL IOS_writer()
!$        WRITE(6,'(A)')'IOS: Info: Writer completed, detatching any helpers'
!$        CALL portioDettachAllHelpers()
!$        WRITE(6,'(A)')'IOS: Info: Writer completed, detatch request completed'
!$      ELSE IF ( this_thread == reader_thread ) THEN
!$        WRITE(6,'(A)')'IOS: Info: Reader thread helper attaching'
!$        CALL portioAddHelper(1)
!$        WRITE(6,'(A)')'IOS: Info: Reader thread has detatched'
!$      ELSE
!$        WRITE(6,'(A,I2,A)')'IOS: Info: Thread ',this_thread,                 &
!$            ' is inactive in IOS'
!$      END IF
!$OMP END PARALLEL
!
! --- Tidy up IOS Stash support (server side) ---
!
!$      CALL IOS_stash_server_fini()
!
!
! --- Close IOS queue structure ---
!
!$      CALL IOS_Queue_Finalise()
!
!$    ELSE
!
! --- Initialise IOS Stash support (client side) ---
!
!$      CALL IOS_client_init()
!$      CALL IOS_stash_client_init(IOS_use_async_stash, &
!$          IOS_use_async_dump,IOS_debug_no_subdomaining)
!$    END IF

    ELSE
      WRITE(6,'(A,A)')'IOS: Info: Configuration not IOS capable,',             &
          ' deactivating IOS'
      WRITE(6,'(A,A)')'IOS: Info: Reason:',TRIM(IOS_capable_string)

      DEALLOCATE(io_servers)
      NULLIFY(io_servers)
      DEALLOCATE(IOS_BCast_Comm)
      NULLIFY(IOS_BCast_Comm)
    END IF

    RETURN

  END FUNCTION IOS_Setup

END MODULE IOS_Init

