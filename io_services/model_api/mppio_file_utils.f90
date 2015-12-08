! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

MODULE MPPIO_file_utils

  ! Implement basic file operations. File ops originated through this
  ! module will execute immediately if no IO servers are present or will be
  ! queue on an IO server with fence semantics if IO servers are active. Fence
  ! semantics mandate that the operation will take place in-order with respect
  ! to any other IO server queued operations and the operation will be delayed
  ! until pending IO ops are completed.
  !
  ! Commands should normally be issued via high level interfaces;
  !    file_copy(<src>,<dest>)
  !    file_delete(<file>)
  !    create_directory(<name>)
  !
  ! which will be translated into a cannonical form and passed to
  ! the generic function;
  !    file_action()
  ! which in turn calls an implementation in c (portio).
  !    <error code>=file_op(<cannonical string>)
  !
  ! The optional force_local argument specifies that (if .TRUE.) the
  ! operation will be perfomed locally and immediately if the calling rank is
  ! rank 0 and the op can be assumed to have been completed on rank 0 upon
  ! return. Other ranks can not make this assumption.
  !
  ! An explicit initialisation on all ranks is required before use.

  USE io
  USE IOS_Client_Queue

  USE yomhook, ONLY :                                                          &
      lhook,                                                                   &
      dr_hook
  USE parkind1, ONLY :                                                         &
      jprb,                                                                    &
      jpim

  USE um_parvars
  USE ereport_mod
  USE PrintStatus_mod

  IMPLICIT NONE

  INTEGER, PARAMETER                    :: file_op_pseudo_unit=8

  ! Profiling
  INTEGER(KIND=jpim),PRIVATE, PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim),PRIVATE, PARAMETER :: zhook_out = 1

  ! Seperator character in cannonical command strings
  CHARACTER(LEN=*), PARAMETER             :: sep=":"

CONTAINS

  SUBROUTINE mppio_file_utils_init()
    IMPLICIT NONE
    REAL(KIND=jprb):: zhook_handle

    IF (lhook) CALL dr_hook('FILE_UTILS:INIT',zhook_in,zhook_handle)

    ! Ensure a mapping exists for file ops.
    ! Although these do not use a unit, we still need to define
    ! which IOS will handle them.
    IF (L_IOS_active() .AND. .NOT. L_IO_Server) THEN
      IF (io_server_for_unit(file_op_pseudo_unit)==IOS_No_Server)              &
          CALL IOS_assign_server_for_unit(file_op_pseudo_unit)
    END IF

    IF (printstatus>=prstatus_normal) THEN
      WRITE(6,'(A,I3)')'MPPIO_File_Utils: Initialised file utils using unit ', &
          file_op_pseudo_unit
    END IF

    IF (lhook) CALL dr_hook('FILE_UTILS:INIT',zhook_out,zhook_handle)

  END SUBROUTINE mppio_file_utils_init

  SUBROUTINE file_touch(theFile,force_local)
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN)  :: theFile  ! to read
    LOGICAL, OPTIONAL              :: force_local
    CHARACTER (LEN=*), PARAMETER   :: cmd="touch"
    CHARACTER (LEN=IOS_String_Max) :: lineBuffer
    REAL(KIND=jprb)                :: zhook_handle

    IF (lhook) CALL dr_hook('FILE_UTILS:FILE_TOUCH',zhook_in,zhook_handle)

    ! Check string is small enough to transmit.
    IF (LEN_TRIM(theFile) +LEN_TRIM(cmd)+ 3*LEN_TRIM(sep)                      &
        > IOS_String_Max) THEN
      WRITE(6,'(A,A)')'String too long theFile =',theFile
      CALL ereport('file_utils:file_delete',99,                                &
          'Arguments too long')
    END IF

    WRITE(lineBuffer,'(A,A,A,A,A)')                                            &
        TRIM(cmd),                                                             &
        TRIM(sep),                                                             &
        TRIM(sep),                                                             &
        TRIM(theFile),                                                         &
        TRIM(sep)

    CALL file_action(TRIM(lineBuffer),force_local)

    IF (lhook) CALL dr_hook('FILE_UTILS:FILE_TOUCH',zhook_out,zhook_handle)

  END SUBROUTINE file_touch

  SUBROUTINE file_copy(src,dest,force_local)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)    :: src  ! to read
    CHARACTER(LEN=*),INTENT(IN)    :: dest ! to write
    LOGICAL, OPTIONAL            :: force_local
    CHARACTER(LEN=*), PARAMETER    :: cmd="cp"
    CHARACTER (LEN=IOS_String_Max)  :: lineBuffer
    REAL(KIND=jprb)              :: zhook_handle

    IF (lhook) CALL dr_hook('FILE_UTILS:FILE_COPY',zhook_in,zhook_handle)

    ! Check string is small enough to transmit.
    IF (LEN_TRIM(src)+LEN_TRIM(dest)+LEN_TRIM(cmd)+3*LEN_TRIM(sep)             &
        > IOS_String_Max) THEN
      WRITE(6,'(A,A)')'Str too long src =',src
      WRITE(6,'(A,A)')'             dest=',dest
      CALL ereport('file_utils:file_copy',99,                                  &
          'Arguments too long')
    END IF

    WRITE(lineBuffer,'(A,A,A,A,A,A)')                                          &
        TRIM(cmd),                                                             &
        TRIM(sep),                                                             &
        TRIM(sep),                                                             &
        TRIM(src),                                                             &
        TRIM(sep),                                                             &
        TRIM(dest)

    CALL file_action(TRIM(lineBuffer),force_local)

    IF (lhook) CALL dr_hook('FILE_UTILS:FILE_COPY',zhook_out,zhook_handle)

  END SUBROUTINE file_copy

  SUBROUTINE dir_create(theDirectory,force_local)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)    :: theDirectory
    LOGICAL, OPTIONAL            :: force_local
    CHARACTER(LEN=*), PARAMETER    :: cmd="mkdir"
    CHARACTER (LEN=IOS_String_Max)  :: lineBuffer
    REAL(KIND=jprb)              :: zhook_handle

    IF (lhook) CALL dr_hook('FILE_UTILS:DIR_CREATE',zhook_in,zhook_handle)

    ! Check string is small enough to transmit.
    IF (LEN_TRIM(theDirectory) +LEN_TRIM(cmd)+ 3*LEN_TRIM(sep)                 &
        > IOS_String_Max) THEN
      WRITE(6,'(A,A)')'Str too long theDirectory =',theDirectory
      CALL ereport('file_utils:file_delete',99,                                &
          'Arguments too long')
    END IF

    WRITE(lineBuffer,'(A,A,A,A,A)')                                            &
        TRIM(cmd),                                                             &
        TRIM(sep),                                                             &
        TRIM(sep),                                                             &
        TRIM(theDirectory),                                                    &
        TRIM(sep)

    CALL file_action(TRIM(lineBuffer),force_local)

    IF (lhook) CALL dr_hook('FILE_UTILS:DIR_CREATE',zhook_out,zhook_handle)

  END SUBROUTINE dir_create

  SUBROUTINE file_delete(theFile,force_local,silent)
    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)    :: theFile  ! to delete
    LOGICAL, OPTIONAL              :: force_local
    LOGICAL, OPTIONAL              :: silent
    LOGICAL                        :: lsilent
    CHARACTER(LEN=*), PARAMETER    :: cmd="rm"
    CHARACTER(LEN=3)               :: option
    CHARACTER (LEN=IOS_String_Max) :: lineBuffer
    REAL(KIND=jprb)                :: zhook_handle

    IF (lhook) CALL dr_hook('FILE_UTILS:FILE_DEL',zhook_in,zhook_handle)

    ! Check string is small enough to transmit.
    IF (LEN_TRIM(theFile) + LEN_TRIM(cmd)+ 3*LEN_TRIM(sep)                     &
        > IOS_String_Max) THEN
      WRITE(6,'(A,A)')'Str too long theFile =',theFile
      CALL ereport('file_utils:file_delete',99,                                &
          'Arguments too long')
    END IF

    lsilent=.FALSE.
    IF (PRESENT(silent)) lsilent=silent

    WRITE(option,'(A)')''
    IF (lsilent) WRITE(option,'(A)')'-f'

    WRITE(lineBuffer,'(A,A,A,A,A,A)')                                          &
        TRIM(cmd),                                                             &
        TRIM(sep),                                                             &
        TRIM(option),                                                          &
        TRIM(sep),                                                             &
        TRIM(theFile),                                                         &
        TRIM(sep)

    CALL file_action(TRIM(lineBuffer),force_local)

    IF (lhook) CALL dr_hook('FILE_UTILS:FILE_DEL',zhook_out,zhook_handle)

  END SUBROUTINE file_delete

  SUBROUTINE file_action(command,force_local)
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)   :: command
    LOGICAL, OPTIONAL            :: force_local
    LOGICAL                      :: l_force_local ! localised copy of above
    INTEGER                      :: pe,i          ! loop counters
    INTEGER                      :: error         ! return code from c
    INTEGER                      :: qHandle       ! IOS operation handle
    INTEGER, EXTERNAL            :: file_op       ! c method to do the work
    TYPE(IOS_metadata_type),                                                   &
        POINTER                  :: metadata
    REAL(KIND=jprb)              :: zhook_handle

    IF (lhook) CALL dr_hook('FILE_UTILS:FILE_ACTION',zhook_in,zhook_handle)

    l_force_local=.FALSE.
    IF (PRESENT(force_local)) l_force_local=force_local

    IF (L_IOS_active()    .AND.                                                &
        .NOT. L_IO_Server .AND.                                                &
        .NOT. l_force_local) THEN

      IF (io_server_for_unit(file_op_pseudo_unit) ==                           &
          IOS_No_Server) THEN
        CALL ereport('file_utils:file_action',99,                              &
            'Module used before initialisation')
      END IF

      IF (model_rank == 0) THEN

        DO pe=1,IOS_Server_Groups
          qHandle=ios_init_md(-1*io_servers(pe,1),                             &
              IOS_No_location,                                                 &
              IOS_Action_FileOp)
          metadata => IOS_Attach_Metadata(qHandle)

          ! Command
          DO i=1,IOS_String_MAX
            metadata % string(i:i)=' '
          END DO
          WRITE(metadata % string,'(A)')command
          metadata % name_length=LEN_TRIM(command)
          metadata % intent = -1
          CALL IOS_Send(qHandle,hasString=.TRUE.)
        END DO
      END IF

    ELSE
      IF (model_rank == 0) THEN
        WRITE(6,'(A,A)')'MPPIO: file op: ',TRIM(command)
        error=file_op(command,LEN_TRIM(command))
        IF (error /= 0) THEN
          CALL ereport('file_utils:file_action',error,                         &
              'A file action failed in the C level implementation,'//          &
              ' see prior messages')
        ELSE
          WRITE(6,'(A)')'MPPIO: file op completed'
        END IF
      END IF
    END IF

    IF (lhook) CALL dr_hook('FILE_UTILS:FILE_ACTION',zhook_out,zhook_handle)

  END SUBROUTINE file_action

END MODULE MPPIO_file_utils
