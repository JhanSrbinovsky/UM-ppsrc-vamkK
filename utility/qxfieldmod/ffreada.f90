! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine: FFREAD and others (see below) ----------------------------
!
! Purpose: To read a   direct access PP file  and convert it to a
! sequential file read to be passed across to the IBM
!
! Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!
! Logical components covered:
!
! Project task:
!
! External documentation:
!
! -------------------------------------------------------------------
! Interface and arguments: ------------------------------------------
!
!    IEXTRA(1) == 0  ! unpacking is required
!    IEXTRA(2) == 0  ! lookup table entry deleted after access.
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Small execs
      SUBROUTINE FFREADA(IPROJ,FCT,ITYPE,INT_LEVEL,PPUNIT,FIELD,IDIM,   &
     &ILABEL,RLABEL,IEXTRA,TEST,PP_LEN2_LOOKUP,LEN1_LOOKUP,PP_FIXHD,    &
     &IWA,LEN1_RECORD,MAXFF,RECORD,PFNO,DATA_ADD,ICODE,CMESSAGE)
      USE IO
        USE ereport_mod, ONLY : ereport
        IMPLICIT NONE

      CHARACTER(LEN=*) cmessage
      INTEGER                                                           &
     &     MAXFF                                                        &
                                  !IN  Max number of opened files
     &    ,LEN1_LOOKUP                                                  &
                                  !IN  first dimension of the lookup
     &    ,LEN1_RECORD                                                  &
                                  !INOUT First dimension of record
     &    ,PP_LEN2_LOOKUP                                               &
                                  !IN  secnd dimension of the lookup
     &    ,IPROJ                                                        &
                                  !IN  map projection of data to read
     &    ,FCT                                                          &
                                  !IN  forecast period in hours
     &    ,ITYPE                                                        &
                                  !IN  M08 FIELD field type
     &    ,INT_LEVEL                                                    &
                                  !IN  LEVEL code (could be real)
     &    ,PPUNIT                                                       &
                                  !IN  unit no of required fieldsfile
     &    ,IDIM                                                         &
                                  !IN  dimension of FIELD
     &    ,ILABEL(45)                                                   &
                                  !OUT holds integer part of LOOKUP
     &    ,IEXTRA(10)                                                   &
                                  !IN  spare for future use
     &    ,ICODE                                                        &
                                  !OUT return code
     &    ,MAXPP                                                        &
                                  !    maximum number of unit number
     &    ,DATA_ADD                                                     &
                                  !IN  The word address of the data.
     &    ,PFNO                   !INOUT No of fields files opened
      INTEGER                                                           &
     &     PP_FIXHD(*),                                                 &
                                              !IN PPfile fixed header
     &     LOOKUP(LEN1_LOOKUP,PP_LEN2_LOOKUP) !OUTinteger lookup
      REAL                                                              &
     &     FIELD(IDIM)                                                  &
                                  !OUT array holding final output data.
     &    ,RLABEL(19)                                                   &
                                  !OUT holds real part of LOOKUP
     &    ,REAL_LEVEL             !IN  LEVEL code (could be real)
      LOGICAL                                                           &
     &     RECORD(LEN1_RECORD,MAXFF) !INOUT Record of the field no read
!     LOCAL VARIABLES
      INTEGER                                                           &
     &     I                                                            &
                                  ! local counter
     &    ,J                      ! local counter
      INTEGER                                                           &
     &     IX                                                           &
                                  ! used as a dummy variable in UNIT
     &    ,IWA                                                          &
                                  ! Word address in call SETPOS
     &    ,IK                                                           &
                                  ! Word address in call SETPOS
     &    ,ICOUNT                                                       &
                                  ! Counter
     &    ,LEN_IO                                                       &
                                  ! Length of data transferred from BUF
     &    ,LEN_IO_EXPECTED                                              &
                                  ! Length od data expected in transfer
     &    ,LENGTH_OF_DATA                                               &
                                  ! Length of a particular field
     &    ,ADDR                                                         &
                                  ! Address of a field in the data store
     &    ,IN_LBVC                ! Local copy of LBVC required to searc
      real                                                              &
     &     A_IO                   ! OUTPUT from UNIT command
      LOGICAL                                                           &
     &     TEST
!
!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS LOOKUP GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!     READ IN LOOKUP TABLE  IF FIRST TIME THRO
!
!     Read in the LOOKUP table.
!
        CALL SETPOS(PPUNIT,IWA,ICODE) ! C coded routine
        LEN_IO_EXPECTED=PP_LEN2_LOOKUP*LEN1_LOOKUP
        CALL BUFFIN(PPUNIT,LOOKUP,LEN_IO_EXPECTED,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_IO_EXPECTED) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('Buffer in lookup table   ',A_IO,LEN_IO,         &
     &    LEN_IO_EXPECTED)
          CMESSAGE='FFREADA: I/O error reading LOOKUP TABLE  '
          ICODE=3
          RETURN
        ENDIF
      IF(IEXTRA(2) == 0) THEN  ! Allows duplicate entries to be read
        IF(PP_LEN2_LOOKUP >  LEN1_RECORD) THEN
          CMESSAGE='FFREADA: LEN1_RECORD NOT LARGE ENOUGH    '
          ICODE=4
          RETURN
        ENDIF
         DO I=1,LEN1_RECORD
         IF(RECORD(I,PFNO)) THEN
            LOOKUP(14,I)=-99
          ENDIF
        ENDDO
      ENDIF
! DEPENDS ON: ffreadb
      CALL FFREADB      (IPROJ,FCT,ITYPE,INT_LEVEL,PPUNIT,FIELD,IDIM,   &
     &ILABEL,RLABEL,IEXTRA,PP_LEN2_LOOKUP,LEN1_LOOKUP,                  &
     &IWA,LEN1_RECORD,MAXFF,RECORD,PFNO,PP_FIXHD,LOOKUP,LOOKUP,DATA_ADD,&
     &ICODE,CMESSAGE)


      IF(ICODE  /=  0) CALL EREPORT("FFREADA", ICODE, CMESSAGE)

      RETURN
      END SUBROUTINE FFREADA
