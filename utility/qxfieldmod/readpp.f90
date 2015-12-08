! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Small execs
!Purpose: To read a direct access PP file and convert it to a sequential
!         file read to be passed across to the IBM
!
      SUBROUTINE READPP(LEN_INTHD,LEN_REALHD,LEN1_LEVDPC,LEN2_LEVDPC,   &
     &   LEN1_ROWDPC,LEN2_ROWDPC,LEN1_COLDPC,LEN2_COLDPC,               &
     &   LEN1_LOOKUP,LEN2_LOOKUP,LEN_FIXHD,PP_FIXHD,LOOKUP,ROOKUP,      &
     &   PPUNIT1,PPUNIT2,ICODE,CMESSAGE)
      USE filenamelength_mod, ONLY :                                    & 
          filenamelength
      USE IO
      IMPLICIT NONE
      INTEGER                                                           &
     &     LEN_FIXHD                                                    &
     &    ,LEN_INTHD                                                    &
     &    ,LEN_REALHD                                                   &
     &    ,LEN_LEVDPC                                                   &
     &    ,LEN_ROWDPC                                                   &
     &    ,LEN_COLDPC                                                   &
     &    ,LEN1_LEVDPC                                                  &
     &    ,LEN2_LEVDPC                                                  &
     &    ,LEN1_ROWDPC                                                  &
     &    ,LEN2_ROWDPC                                                  &
     &    ,LEN1_COLDPC                                                  &
     &    ,LEN2_COLDPC                                                  &
     &    ,LEN_LOOKUP                                                   &
     &    ,LEN1_LOOKUP                                                  &
     &    ,LEN2_LOOKUP                                                  &
     &    ,LEN1_LOOKNEW                                                 &
     &    ,LEN2_LOOKNEW                                                 &
     &    ,LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                              &
     &    ,PP_INTHD(LEN_INTHD)                                          &
     &    ,PP_FIXHD(LEN_FIXHD)                                          &
     &    ,LEN_IO                                                       &
     &    ,ICODE                                                        &
     &    ,PPUNIT1                                                      &
     &    ,PPUNIT2
      REAL                                                              &
     &     ROOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                              &
     &    ,PP_REALHD(LEN_REALHD)                                        &
     &    ,PP_LEVDPC(LEN1_LEVDPC*LEN2_LEVDPC+1)                         &
     &    ,PP_ROWDPC(LEN1_ROWDPC*LEN2_ROWDPC)                          &
     &    ,PP_COLDPC(LEN1_COLDPC*LEN2_COLDPC)                          &
     &    ,A_IO
      CHARACTER(LEN=*) cmessage
      CHARACTER(LEN=filenamelength) :: outfile
! Local variables
      INTEGER                                                           &
     &     START_BLOCK                                                  &
     &    ,NENT                                                         &
     &    ,K                                                            &
     &    ,Kk                                                           &
     &    ,iwa                                                          &
     &    ,RECL                                                         &
     &    ,IERR                                                         &
     &    ,PP_LEN_INTHD                                                 &
     &    ,PP_LEN_REALHD                                                 

!---------------------------------------------------------------------
      PP_LEN_REALHD=PP_FIXHD(106)
      PP_LEN_INTHD=PP_FIXHD(101)

!---------------------------------------------------------------------
      LEN_LEVDPC=LEN1_LEVDPC*LEN2_LEVDPC
      LEN_ROWDPC=LEN1_ROWDPC*LEN2_ROWDPC
      LEN_COLDPC=LEN1_COLDPC*LEN2_COLDPC
      LEN_LOOKUP=LEN1_LOOKUP*LEN2_LOOKUP
! The calculation of LEN_LEVDPC has PLUS 1 which is only true
! for PP headers and not model headers, hopefully the PLUS one will
! be removed as it is inconsistent)
      START_BLOCK=LEN_FIXHD+1
!L---------------------------------------------------------------
!L  Read in the integer constants
!L---------------------------------------------------------------
      IF(LEN_INTHD >  0) THEN  ! Integer constants to be read in
        IF(PP_FIXHD(100) /= START_BLOCK) THEN   ! Address incorrect
! DEPENDS ON: poserror
          CALL POSERROR('integer constants',START_BLOCK,100,            &
     &    PP_FIXHD(100))
          CMESSAGE=' READPP :  Adressing Conflict'
          ICODE=2
          RETURN
        ENDIF
        CALL BUFFIN(PPUNIT1,PP_INTHD,LEN_INTHD,LEN_IO,A_IO)
        WRITE(6,*)pp_inthd
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_INTHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR(' Buffer in of Integer constants',A_IO,LEN_IO    &
     &  ,  LEN_INTHD)
          CMESSAGE='READPP : I/O error'
          ICODE=3
          RETURN
        ENDIF
        START_BLOCK=START_BLOCK+LEN_INTHD
      ENDIF
!L---------------------------------------------------------------
!L  Read in the real constants
!L---------------------------------------------------------------
      IF(LEN_REALHD >  0) THEN  ! Real constants to be read in
        IF(PP_FIXHD(105) /= START_BLOCK) THEN   ! Address incorrect
! DEPENDS ON: poserror
          CALL POSERROR('Real constants',START_BLOCK,100,               &
     &    PP_FIXHD(105))
          CMESSAGE=' READPP :  Adressing Conflict'
          ICODE=4
          RETURN
        ENDIF
        CALL BUFFIN(PPUNIT1,PP_REALHD,LEN_REALHD,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_REALHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR(' Buffer in of Real constants',A_IO,LEN_IO       &
     &    ,LEN_REALHD)
          CMESSAGE='READPP : I/O error'
          ICODE=5
          RETURN
        ENDIF
        START_BLOCK=START_BLOCK+LEN_REALHD
      ENDIF
!L---------------------------------------------------------------
!L  Read in the level dependant constants
!L---------------------------------------------------------------
      IF(LEN_LEVDPC >  0) THEN  ! Level dep constants to be read in
        IF(PP_FIXHD(110) /= START_BLOCK) THEN   ! Address incorrect
! DEPENDS ON: poserror
          CALL POSERROR('Level depndt constants',START_BLOCK,100,       &
     &    PP_FIXHD(110))
          CMESSAGE=' READPP :  Adressing Conflict'
          ICODE=6
          RETURN
        ENDIF
        CALL BUFFIN(PPUNIT1,PP_LEVDPC,LEN_LEVDPC,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_LEVDPC) THEN
! DEPENDS ON: ioerror
          CALL IOERROR(' Buffer in of Level constants',A_IO,LEN_IO      &
     &    ,LEN_LEVDPC)
          CMESSAGE='READPP : I/O error'
          ICODE=7
          RETURN
        ENDIF
        START_BLOCK=START_BLOCK+LEN_LEVDPC
      ENDIF
!L---------------------------------------------------------------
!L  Read in the Row dependant constants
!L---------------------------------------------------------------
      IF(LEN_ROWDPC >  0) THEN  ! row dep constants to be read in
        IF(PP_FIXHD(115) /= START_BLOCK) THEN   ! Address incorrect
! DEPENDS ON: poserror
          CALL POSERROR('Row depndt constants',START_BLOCK,100,         &
     &                   PP_FIXHD(115))
          CMESSAGE=' READPP :  Adressing Conflict'
          ICODE=10
          RETURN
        ENDIF
        CALL BUFFIN(PPUNIT1,PP_ROWDPC,LEN_ROWDPC,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_ROWDPC) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('Buffer in of Row constants',A_IO,LEN_IO,        &
     &                  LEN_ROWDPC)
          CMESSAGE='READPP : I/O error'
          ICODE=11
          RETURN
        ENDIF
        START_BLOCK=START_BLOCK+LEN_ROWDPC
      ENDIF
!L---------------------------------------------------------------
!L  Read in the Col dependant constants
!L---------------------------------------------------------------
      IF(LEN_COLDPC >  0) THEN  ! col dep constants to be read in
        IF(PP_FIXHD(120) /= START_BLOCK) THEN   ! Address incorrect
! DEPENDS ON: poserror
          CALL POSERROR('Col depndt constants',START_BLOCK,100,         &
     &                   PP_FIXHD(120))
          CMESSAGE=' READPP :  Adressing Conflict'
          ICODE=20
          RETURN
        ENDIF
        CALL BUFFIN(PPUNIT1,PP_COLDPC,LEN_COLDPC,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_ROWDPC) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('Buffer in of Col constants',A_IO,LEN_IO,        &
     &                  LEN_COLDPC)
          CMESSAGE='READPP : I/O error'
          ICODE=21
          RETURN
        ENDIF
        START_BLOCK=START_BLOCK+LEN_COLDPC
      ENDIF      
!L---------------------------------------------------------------
!L  Read in the LOOKUP TABLE
!L---------------------------------------------------------------
      IF(LEN_LOOKUP >  0) THEN  ! Lookup Table to be read in
        IF(PP_FIXHD(150) /= START_BLOCK) THEN   ! Address incorrect
          WRITE(6,*) 'READPP : WARNING'
          WRITE(6,*) 'Conflict between start position of Lookup table'
          WRITE(6,*) 'block and pointer in fixed length header: ',      &
     &               'FIXHD(150) = ',PP_FIXHD(150)
          WRITE(6,*) 'Current position in file = ',START_BLOCK,         &
     &               ' words in'
          WRITE(6,*) 'Pointer moved to ',PP_FIXHD(150),' words in'
          CALL SETPOS(PPUNIT1,PP_FIXHD(150)-1,IERR)
        END IF
        CALL BUFFIN(PPUNIT1,LOOKUP,LEN_LOOKUP,LEN_IO,A_IO)
        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_LOOKUP) THEN
! DEPENDS ON: ioerror
          CALL IOERROR(' Buffer in of Lookup table   ',A_IO,LEN_IO      &
     &    ,LEN_LOOKUP)
          CMESSAGE='READPP : I/O error'
          ICODE=9
          RETURN
        ENDIF
        START_BLOCK=START_BLOCK+LEN_LOOKUP
      ENDIF
      WRITE(6,*)' ARRIVED HERE  ',START_BLOCK
      NENT=0
      DO K=1,LEN2_LOOKUP
        IF(LOOKUP(1,K) >  0) THEN
          NENT=NENT+1
        ELSE
          EXIT
        ENDIF
      END DO ! K

      WRITE(6,*)' VALUE OF NENT   ',NENT
      do k=nent-2,nent+1
        WRITE(6,*)'k=',k
        WRITE(6,*) (lookup(kk,k),kk=1,44)
      enddo
!-----------------------------------------------------------------
!    OPEN NEW TARGET FIELDSFILE INITIALISING BY CALLING INITPP
!-----------------------------------------------------------------
!L
!L        Open named file on unit 60
!L
        WRITE(6,*)"*** Opening new file on unit ",pPUNIT2
        CALL GET_FILE(PPUNIT2,OUTFILE,filenamelength,ICODE)
        CALL FILE_OPEN(PPUNIT2,OUTFILE,filenamelength,1,1,ICODE)
!
! DEPENDS ON: init_pp
      CALL INIT_PP(PPUNIT2,'p',LEN1_LOOKUP,LEN2_LOOKUP,PP_FIXHD,        &
     &             PP_INTHD,PP_REALHD,PP_LEVDPC,PP_ROWDPC,PP_COLDPC,    &
     &                       LEN_INTHD,LEN_REALHD,LEN1_LEVDPC,          &
     &             LEN2_LEVDPC,LEN1_ROWDPC,LEN2_ROWDPC,                 &
     &             LEN1_COLDPC,LEN2_COLDPC,PP_LEN_INTHD,                &
     &             PP_LEN_REALHD,ICODE,CMESSAGE)

      IF(ICODE /= 0) THEN
        WRITE(7,'(A,I2)')' ICODE EQUAL TO ',ICODE
        WRITE(7,'(A80)') CMESSAGE
        !FixMe : ereport here?
        RETURN
      ENDIF
      LEN1_LOOKNEW=LEN1_LOOKUP
      LEN2_LOOKNEW=LEN2_LOOKUP
! DEPENDS ON: control
      CALL CONTROL(PPUNIT1,PPUNIT2,LEN1_LOOKNEW,LEN2_LOOKNEW,           &
     &             LOOKUP,PP_INTHD,LEN_INTHD,                           &
     &             PP_FIXHD,LEN_FIXHD,ICODE,CMESSAGE,NENT)
      IF(ICODE /= 0) THEN
        WRITE(7,'(A,I2)') ICODE
        WRITE(7,'(A80)') CMESSAGE
        !FixMe : ereport here?
        RETURN
      ENDIF

      RETURN
      END SUBROUTINE READPP
