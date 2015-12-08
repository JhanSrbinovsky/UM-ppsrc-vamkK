! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

MODULE MPPIO_Job_Control

! Implement operations for job control commuinication. We may or may not
! have an IO server present so we need to switch between local and remote ops
! to write to a pipe file.
!
! Method: We either write locally on pe 0 if there is no io server, or send the
! instruction to the IO servers. As the message must not arrive at the pipe
! file before associated data is commited, and because we don't know which IO
! server (or more accurately which unit) the message is really associated with
! (and generally to preserve ordering semantics in the pipe file) we send
! the instruction to all IO server tasks, such that they can synchronise with
! each other in an appropriate fashion (this module doesn't care how they do
! that, we simply discharge responsibility to IOS)

  USE io
  USE IOS_Client_Queue
  USE mpl
  USE mppio_job_control_common, ONLY :                                         &
      JC_Write,                                                                &
      JC_Unit,                                                                 &
      JC_Filename
  USE yomhook, ONLY :                                                          &
      lhook,                                                                   &
      dr_hook
  USE parkind1, ONLY :                                                         &
      jprb,                                                                    &
      jpim


  USE ereport_mod, ONLY : ereport
  USE UM_ParVars
  IMPLICIT NONE

! API for model code. Send a command/string to the pipe file
  INTERFACE jobCntrl
    MODULE PROCEDURE                                                           &
        jobCntrl_String,                                                       &
        jobCntrl_DescriptorString,                                             &
        jobCntrl_CommandString
  END INTERFACE ! CommandString is preffered

! Force clients to use the interface
  PRIVATE jobCntrl_String
  PRIVATE jobCntrl_DescriptorString
  PRIVATE jobCntrl_CommandString

! Profiling
  INTEGER(KIND=jpim),PRIVATE, PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim),PRIVATE, PARAMETER :: zhook_out = 1

CONTAINS

! Open the pipe file
  SUBROUTINE jobCntrlInit()
    IMPLICIT NONE
    INTEGER                      :: ierror
    INTEGER                      :: qHandle
    TYPE(IOS_metadata_type),                                                   &
        POINTER                  :: metadata
    REAL(KIND=jprb)              :: zhook_handle
    IF (lhook) CALL dr_hook('JOBCNTRL:INIT',zhook_in,zhook_handle)

    CALL get_file(jc_unit,jc_filename,256,ierror)

    IF (L_IOS_active()) THEN

      ! Check env is small enough to transmit.
      IF (LEN_TRIM(jc_filename) > IOS_String_Max) THEN
        CALL ereport('job_control:jobCntrlInit',99,                            &
            'filename too big to open')
      END IF

      IF (io_server_for_unit(jc_unit)==IOS_No_Server)                          &
          CALL IOS_assign_server_for_unit(jc_unit)

      IF (mype == 0) THEN

        qHandle=IOS_Init_MD(jc_unit,IOS_No_Location,IOS_Action_Open_Pipe)
        metadata => IOS_Attach_Metadata(qHandle)
        metadata % name_length           = LEN_TRIM(jc_filename)
        metadata % string(1:LEN_TRIM(jc_filename))                             &
            = jc_filename(1:LEN_TRIM(jc_filename))
        metadata % string(LEN_TRIM(jc_filename)+1:)=' '

        CALL IOS_Send(qHandle,hasString=.TRUE.)

        WRITE(6,'(A,A,A,I2,A)')                                                &
            'Opened pipe file ',TRIM(jc_filename),                             &
            ' using IOS for job control'
      END IF
    ELSE
      IF (mype == 0) THEN
        OPEN(jc_unit,FILE=jc_filename)
        WRITE(6,'(A,A,A)')'Opened pipe file ',                                 &
            TRIM(jc_filename),' locally for job control'
      END IF
    END IF
    IF (lhook) CALL dr_hook('JOBCNTRL:INIT',zhook_out,zhook_handle)

  END SUBROUTINE jobCntrlInit

  ! Close the pipe file
  SUBROUTINE jobCntrlFini()
    IMPLICIT NONE
    INTEGER                      :: ierror
    INTEGER                      :: qHandle
    REAL(KIND=jprb)              :: zhook_handle

    IF (lhook) CALL dr_hook('JOBCNTRL:FINI',zhook_in,zhook_handle)

    IF (mype == 0) THEN
      IF (L_IOS_active()) THEN

        qHandle=IOS_Init_MD(jc_unit,IOS_No_location,IOS_Action_Close_Pipe)

        CALL IOS_Send(qHandle)

      ELSE
        CLOSE(jc_unit)
      END IF
    END IF
    IF (lhook) CALL dr_hook('JOBCNTRL:FINI',zhook_out,zhook_handle)

  END SUBROUTINE jobCntrlFini

! Send a string to the pipe file
  SUBROUTINE jobCntrl_String(string)
    IMPLICIT NONE
    CHARACTER(LEN=*)               :: string ! to write
    INTEGER :: ierror
    INTEGER :: pe
    INTEGER :: qHandle
    TYPE(IOS_metadata_type),                                                   &
        POINTER                  :: metadata
    REAL(KIND=jprb)              :: zhook_handle

    IF (lhook) CALL dr_hook('JOBCNTRL',zhook_in,zhook_handle)

    IF (mype == 0) THEN
      IF (L_IOS_active()) THEN

        ! Check string is small enough to transmit.
        IF (LEN_TRIM(string) > IOS_String_Max) THEN
          WRITE(6,'(A,I8)')'Str too long len=',LEN_TRIM(string)
          WRITE(6,'(A,A)') '             str=',string
          CALL ereport('job_control:jobCntrlInit',99,                          &
              'Command String too long')
        END IF

        DO pe=1,IOS_Server_Groups
          qHandle=IOS_Init_MD(-1*io_servers(pe,1),                             &
              IOS_No_location,                                                 &
              IOS_Action_Release)
          metadata => IOS_Attach_Metadata(qHandle)
          metadata % string(1:LEN_TRIM(string))                                &
              = string(1:LEN_TRIM(string))
          metadata % string(LEN_TRIM(string)+1:)=' '
          metadata % name_length=LEN_TRIM(string)
          metadata % intent = -1
          CALL IOS_Send(qHandle,hasString=.TRUE.)
        END DO
      ELSE
        CALL jc_write(-1, TRIM(string))
      END IF
    END IF

    IF (lhook) CALL dr_hook('JOBCNTRL',zhook_out,zhook_handle)

  END SUBROUTINE jobCntrl_String

! Send a string to the pipe file as a formating descriptor and
! and argument to the descriptor - e.g. '('Please ',A,' me' )','feed'
  SUBROUTINE jobCntrl_DescriptorString(descriptor,string)
    IMPLICIT NONE
    CHARACTER(LEN=*)               :: descriptor ! fortran format descriptor
    CHARACTER(LEN=*)               :: string     ! to use with the descriptor
                                               ! probably a filename
    CHARACTER (LEN=IOS_String_Max)  :: lineBuffer
    REAL(KIND=jprb)              :: zhook_handle

    IF (lhook) CALL dr_hook('JOBCNTRL',zhook_in,zhook_handle)
    WRITE(lineBuffer,descriptor)string
    CALL jobCntrl(lineBuffer)
    IF (lhook) CALL dr_hook('JOBCNTRL',zhook_out,zhook_handle)

  END SUBROUTINE jobCntrl_DescriptorString

! Send a command to the pipe file as a numeric parameter
! and filename argument to that  - e.g. jc_delete,'an_unwanted_file'
  SUBROUTINE jobCntrl_CommandString(command,string)
    IMPLICIT NONE
    CHARACTER(LEN=*)               :: string  ! usually a filename
    INTEGER,INTENT(IN)           :: command ! parameterised values in
                                            ! mppio_job_control_common
    INTEGER                      :: ierror
    INTEGER                      :: pe
    INTEGER                      :: qHandle
    TYPE(IOS_metadata_type),                                                   &
        POINTER                  :: metadata
    REAL(KIND=jprb)              :: zhook_handle

    IF (lhook) CALL dr_hook('JOBCNTRL',zhook_in,zhook_handle)

    IF (mype == 0) THEN
      IF (L_IOS_active()) THEN

        ! Check string is small enough to transmit.
        IF (LEN_TRIM(string) > IOS_String_Max) THEN
          WRITE(6,'(A,I8)')'Str too long len=',LEN_TRIM(string)
          WRITE(6,'(A,A)') '             str=',string
          CALL ereport('job_control:jobCntrlInit',99,                          &
              'Command String too long')
        END IF

        DO pe=1,IOS_Server_Groups
          qHandle=IOS_Init_MD(-1*io_servers(pe,1),                             &
              IOS_No_location,                                                 &
              IOS_Action_Release)
          metadata => IOS_Attach_Metadata(qHandle)
          metadata % intent = command
          metadata % string(1:LEN_TRIM(string))                                &
              = string(1:LEN_TRIM(string))
          metadata % string(LEN_TRIM(string)+1:)=' '
          metadata % name_length=LEN_TRIM(string)
          CALL IOS_Send(qHandle,hasString=.TRUE.)
        END DO

      ELSE
        CALL jc_write(command,TRIM(string))
      END IF
    END IF
    IF (lhook) CALL dr_hook('JOBCNTRL',zhook_out,zhook_handle)
  END SUBROUTINE jobCntrl_CommandString

END MODULE MPPIO_Job_Control
