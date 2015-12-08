! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Small execs
      SUBROUTINE READ_REC_ffread1a(FIELD,NUM_CRAY_WORDS,IWA,PPUNIT,     &
     &                    ICODE,CMESSAGE)
      USE IO
        USE ereport_mod, ONLY : ereport
        IMPLICIT NONE
      CHARACTER(LEN=*) cmessage
      INTEGER                                                           &
     &     ICODE                                                        &
                                  !OUT return code
     &    ,NUM_CRAY_WORDS                                               &
                                  !IN  No of CRAY words holding the data
     &    ,PPUNIT                                                       &
                                  !IN  FT no of the PP FILE
     &    ,IWA                    !IN  WORD address of field to be read
      REAL                                                              &
     &     FIELD(NUM_CRAY_WORDS)  !OUT array holding data
!    LOCAL VARIABLES
      INTEGER                                                           &
     &     I                                                            &
                                  ! local counter
     &    ,J                                                            &
                                  ! local counter
     &    ,IX                                                           &
                                  ! used in the UNIT command
     &    ,LEN_IO                 ! used for call to BUFFIN
      REAL                                                              &
     &     A_IO                   ! used for call to BUFFIN
      CALL SETPOS(PPUNIT,IWA,ICODE) ! C coded routine
      CALL BUFFIN(PPUNIT,FIELD,NUM_CRAY_WORDS,LEN_IO,A_IO)

      IF(ICODE  /=  0) CALL EREPORT("READ_REC", ICODE, CMESSAGE)
      RETURN
      END SUBROUTINE READ_REC_ffread1a
