! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routines to initialise and write to PEER output files which can then
! be analysed with the PEER utility.
!
      SUBROUTINE PEER_INITIALISE(                                       &
     &  ICODE,CMESSAGE)

! Description:
!   Creates output file to hold arrays and writes simple header.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Subroutine Arguments:

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE IO
      USE UM_ParVars
      IMPLICIT NONE

      INTEGER                                                           &
     &  ICODE                  ! INOUT: Error return code

      CHARACTER(LEN=80) CMESSAGE     ! INOUT: Error message

! COMDECKS and common blocks

!
! This Comdeck declares and stores the logical and integer
! variables used in the time-step control for writing general
! data.
!
!
! Switch which activates output of arrays for Peer utility
      LOGICAL L_PEER

! Switches which activate dump writing
      LOGICAL                                                           &
     &  L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                  &
     &  ,L_WRIT_INIT

! Timesteps for dump writing
      INTEGER                                                           &
     &  T_WRITD1_START                                                  &
                              ! First timestep
     &  ,T_WRITD1_END                                                   &
                              ! Last timestep
     &  ,T_WRITD1_INT         ! Timestep interval between dumps

      INTEGER                                                           &
     &  PEER_VN                  ! Version of PEER utility

      NAMELIST/NLSTWRITEDATA/                                           &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT

      COMMON/WRITEDATA/                                                 &
     &  L_PEER,PEER_VN                                                  &
     &  ,T_WRITD1_START  ,T_WRITD1_END  ,T_WRITD1_INT                   &
     &  ,L_WRIT_OCNSTEP ,L_WRIT_WAVSTEP                                 &
     &  ,L_WRIT_INIT
!

      INTEGER                                                           &
     &  PEER_ADDRESS            ! Address pointer of peer file

      COMMON /PEER_ADDRESS / PEER_ADDRESS

! Local variables:

      INTEGER                                                           &
     &  HEADER(4)               ! Header for output file

      REAL                                                              &
     &  IOSTAT                  ! Return code from I/O routines

      INTEGER                                                           &
     &  PUNIT                                                           &
                                ! Unit number for output data file
     &  ,SUNIT                                                          &
                                ! Unit number for output summary file
     &  ,IOREQ                                                          &
                                ! Stores number of words to write out
     &  ,LENIO                  ! Length actually written out

      CHARACTER(LEN=30)                                                      &
     &  PEER_FILENAME                                                   &
                                ! Filename for data
     &  ,SUMM_FILENAME          ! Filename for summary

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('PEER_INITIALISE',zhook_in,zhook_handle)
      IF (L_PEER) THEN
! Initialise pointer
        PEER_ADDRESS=1

! Get file names and unit numbers, and open files for writing
! DEPENDS ON: peer_file_details
        CALL PEER_FILE_DETAILS(PUNIT,PEER_FILENAME,SUNIT,SUMM_FILENAME)
! File that holds the data
        CALL FILE_OPEN(PUNIT,PEER_FILENAME,20,1,1,ICODE,&
            ioLocality=ioAllLocal)
        IF (ICODE /= 0)THEN
          WRITE(6,*)                                                    &
     &      'PEER: Error opening file ',PEER_FILENAME,' on unit ',PUNIT
          WRITE(6,*)'Error code ',ICODE
          CMESSAGE='PEER: Error opening file on one of the PEs'
          GOTO 9999
        ENDIF
! File that holds details on each field
        OPEN(UNIT=SUNIT,FILE=SUMM_FILENAME,FORM='FORMATTED',            &
     &    IOSTAT=ICODE)
        IF (ICODE /= 0)THEN
          WRITE(6,*)                                                    &
     &      'PEER: Error opening file ',SUMM_FILENAME,' on unit ',SUNIT
          WRITE(6,*)'Error code ',ICODE
          CMESSAGE='PEER: Error opening file on one of the PEs'
          GOTO 9999
        ENDIF
! Set the write position to the start of the file.

        WRITE(SUNIT,'(4I5)')PEER_VN,nproc_x,nproc_y,mype
      ENDIF ! IF (L_PEER)

 9999 CONTINUE
      IF (lhook) CALL dr_hook('PEER_INITIALISE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE PEER_INITIALISE
