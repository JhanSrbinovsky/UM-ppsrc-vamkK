! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

! This module provides constants and an output method for writing to
! the pipe unit. It is intended to be called only from mppio_job_control
! or an IO server component. Model code should only call mppio_job_control
! but obtain parameters from here.

MODULE mppio_job_control_common
  USE fileNameLength_Mod, ONLY : fileNameLength
  USE yomhook,  ONLY           : lhook, dr_hook
  USE parkind1, ONLY           : jprb, jpim

  IMPLICIT NONE
! Profiling
  INTEGER(KIND=jpim),PRIVATE, PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim),PRIVATE, PARAMETER :: zhook_out = 1

! The fortran unit that gets used for the pipe
  INTEGER, PARAMETER  :: JC_unit=8

! Parameterisations of commands for the server
  INTEGER, PARAMETER  :: JC_wakeup=1
  INTEGER, PARAMETER  :: JC_delete=2
  INTEGER, PARAMETER  :: JC_archiveDump=3
  INTEGER, PARAMETER  :: JC_archiveBndy=4
  INTEGER, PARAMETER  :: JC_archivePpnochart=5
  INTEGER, PARAMETER  :: JC_archive_dump=6
  INTEGER, PARAMETER  :: JC_archivePPChart=7
  INTEGER, PARAMETER  :: JC_plotPPChart=8
  INTEGER, PARAMETER  :: JC_WorkStationSend=9
  INTEGER, PARAMETER  :: JC_Quit=10
  INTEGER, PARAMETER  :: JC_strlen=20
  INTEGER, PARAMETER  :: JC_strnum=10

! String values for each command
  CHARACTER (LEN=*),PARAMETER::                                                &
      JC_strings=                                                              &
      "**WAKEUP**          "//                                                 &
      " DELETE             "//                                                 &
      " ARCHIVE DUMP       "//                                                 &
      " ARCHIVE BNDY       "//                                                 &
      " ARCHIVE PPNOCHART  "//                                                 &
      " ARCHIVE_DUMP       "//                                                 &
      " ARCHIVE PPCHART    "//                                                 &
      " PLOTONLY PPCHART   "//                                                 &
      " HPSEND             "//                                                 &
      "**QUIT**            "
  !    12345678901234567890 = JC_strlen = 20
  CHARACTER (LEN=fileNameLength)    :: JC_filename


CONTAINS

! Write an operation (command) to the pipe file
  SUBROUTINE jc_write(command,string)
    USE IOS_Common, ONLY : IOS_Start_Time
    IMPLICIT NONE
    CHARACTER(LEN=*)    :: string   ! The filename associated with the command,
                                    ! if needed.
    INTEGER, INTENT(IN) :: command  ! The command identifier, see parameterised
                                    ! list above
    INTEGER             :: ierror
    REAL                :: t        ! the time
    REAL(KIND=jprb)     :: zhook_handle
    REAL, EXTERNAL      :: get_wallclock_time

    IF (lhook) CALL dr_hook('JOBCNTRLCMMN:WRITE',zhook_in,zhook_handle)

! DEPENDS ON : get_wallclock_time
    t=get_wallclock_time()

    SELECT CASE(command)
    CASE (-1)
      WRITE(jc_unit,'(A)')TRIM(string)
      WRITE(6,'(A,F10.3,A,A)')'JOB CONTROL:',                                  &
          t-IOS_Start_Time,':',TRIM(string)
    CASE (0)
      WRITE(6,'(A,A)')'JOB CONTROL ERROR: Command field was 0'
    CASE (JC_wakeup,JC_Quit)
      WRITE(jc_unit,'(A)')                                                     &
          JC_strings((command-1)*JC_strlen+1:command*JC_strlen)
      WRITE(6,'(A,F10.3,A,A)')'JOB CONTROL:',                                  &
          t-IOS_Start_Time,':',                                                &
          JC_strings((command-1)*JC_strlen+1:command*JC_strlen)
    CASE Default
      IF (command > jc_strnum .OR. command < -1) THEN
        WRITE(6,'(A,i4)')'JOB Control: command not recognised:',command
      ELSE
        WRITE(jc_unit,"('%%%',1X,A,A)") TRIM(string),                          &
            JC_strings((command-1)*JC_strlen+1:command*JC_strlen)
        WRITE(6      ,"('JOB CONTROL:',F10.3,':%%%',1X,A,A)")                  &
            t-IOS_Start_Time,TRIM(string),                                     &
            JC_strings((command-1)*JC_strlen+1:command*JC_strlen)
      END IF
    END SELECT
! DEPENDS ON: um_fort_flush
    CALL um_fort_flush(8,ierror)
    CALL um_fort_flush(6,ierror)
    IF (lhook) CALL dr_hook('JOBCNTRLCMMN:WRITE',zhook_out,zhook_handle)

  END SUBROUTINE jc_write

END MODULE mppio_job_control_common
