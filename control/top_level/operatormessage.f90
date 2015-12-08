! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Write a message to the operator for operational runs
!
! Subroutine Interface:
      SUBROUTINE OperatorMessage(nproc)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
! Description: OperatorMessage writes a message to the operator for
! operational runs to indicate that jobs have started correctly.
!
!
! Method: Use SHELL routine to fork a process without copying a memory
! image. This allows 1 command msgi to write a messsage to the
! operator. Only activated if environment variable UM_OPER_MESS="true".
! Note: fort_get_env fails if contents of environment variable is longer
! than the available storage, so length of UM_OPER_MESS (or RUNID)
! should not exceed env_char_len (=8 here).
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
! Declarations:

! Global variables (*CALLed COMDECKs etc...):

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER  nproc                  ! Number of processors

! Local parameters:

! Local scalars:
      CHARACTER(LEN=80) shell_arg     ! text for argument of shell call
      CHARACTER(LEN=60) message
      CHARACTER(LEN=8)  env_char      ! text of env. variable

      INTEGER shell_arg_len                                             &
                                 ! length of shell argument
     &       ,command_len                                               &
                                 ! length of shell command
     &       ,icode                                                     &
                                 ! error return code
     &       ,RUNID_char_len                                            &
                                 ! length of RUNID contents
     &       ,env_char_len       ! length of env_char array

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      Parameter(                                                        &
     &           RUNID_char_len=5)

!- End of header
      IF (lhook) CALL dr_hook('OPERATORMESSAGE',zhook_in,zhook_handle)
      message='"xxxx UM successfully running on xxxx PEs "' ! template
      shell_arg='msgi'            ! shell command
      command_len=4               ! length of shell command
      env_char_len=len(env_char)

! Check that env.var. UM_OPER_MESS="true" in top level 1 script
      env_char='false'
      CALL fort_get_env('UM_OPER_MESS',12,env_char,env_char_len,icode)

      IF(icode == 0 .AND. env_char(1:4) == 'true') THEN

! Construct message: place RUNID and no. of PEs in message
        CALL fort_get_env('RUNID',5,env_char,env_char_len,icode)
        IF(icode /= 0) THEN
          write(6,*) 'Routine OperatorMessage: Warning: problem in ',   &
     &               'reading environment variable RUNID, icode=',icode
        ENDIF

        message(2:RUNID_char_len+1)=env_char  ! substitute RUNID
        write(message(34:37),'(I4)') nproc    ! substitute no. of PEs
        shell_arg(command_len+2:)=message

! Send message to operator
        shell_arg_len=len(shell_arg)
        CALL shell(shell_arg,shell_arg_len)
        write(6,*) shell_arg                ! and standard output

      ENDIF                               ! Test on UM_OPER_MESS

      IF (lhook) CALL dr_hook('OPERATORMESSAGE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE OperatorMessage
