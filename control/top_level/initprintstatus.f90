! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialise print status for standard output
!
! Subroutine Interface:
      Subroutine InitPrintStatus

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE IOS, ONLY: IOS_Config

      USE PrintStatus_mod
      
      IMPLICIT NONE
!
! Description:
!   Initialises value of PrintStatus variable to control subsequent
!   printing of standard output to unit 6.
! Method:
! Use environment variable available from UMUI supplied to top level
! SCRIPT file to set switch PrintStatus, accessed locally by both
! fortran and C code.
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Declarations:
!   These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
! Global variables (*CALLed COMDECKs etc...):
! Print status information in CPRINTST:
! Subroutine arguments
! Local parameters:
      CHARACTER(LEN=*) RoutineName
      PARAMETER (   RoutineName='InitPrintStatus')
! Local scalars:
      INTEGER ErrorStatus       ! Error return code
      CHARACTER(LEN=80) PRINT_STATUS   ! Text content of env variable

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!- End of header

      IF (lhook) CALL dr_hook('INITPRINTSTATUS',zhook_in,zhook_handle)
      CALL Fort_Get_Env('PRINT_STATUS',12,PRINT_STATUS,80,ErrorStatus)

      IF(ErrorStatus /= 0) THEN
        WRITE(6,*) RoutineName,': Warning: problem in reading ',        &
     &  'environment variable PRINT_STATUS, ErrorStatus=',ErrorStatus
      ENDIF

      IF    (PRINT_STATUS(1:13) == 'PrStatus_Diag') THEN
        PrintStatus = PrStatus_Diag
      ELSEIF(PRINT_STATUS(1:15) == 'PrStatus_Normal') THEN
        PrintStatus = PrStatus_Normal
      ELSEIF(PRINT_STATUS(1:13) == 'PrStatus_Oper') THEN
        PrintStatus = PrStatus_Oper
      ELSE
        PrintStatus = PrStatus_Min
      ENDIF

! Set PrintStatus in C code for control of I/O messages
      CALL Set_Printstatus(PrintStatus)

      IF(PrintStatus >= PrStatus_Oper) THEN
        write(6,*) RoutineName,': PrintStatus initialised=',PrintStatus
      ENDIF

      CALL IOS_Config(verbosity=printstatus)

      IF (lhook) CALL dr_hook('INITPRINTSTATUS',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE InitPrintStatus
