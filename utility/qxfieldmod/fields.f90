! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine: FIELDS ----------------------------------------------
!
! Purpose: To calculate fields from the Fields File such as those
! normaly derived in the Derived Printfile Program
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Small execs
!
! -------------------------------------------------------------------
! Interface and arguments: ------------------------------------------

      SUBROUTINE FIELDS(PP_FIXHD,LEN_FIXHD,LENBUF,LEN_FIELD,            &
     &                  LOOKUP,ROOKUP,LEN1_LOOKUP,LEN2_LOOKUP,NENT,     &
     &                  STIME_MOD,ETIME_MOD,NFIELDS_MOD,                &
     &                                        MTYP_MOD,MLEVS_MOD,AMULT, &
     &                  STIME_SEL,ETIME_SEL,NFIELDS_SEL,                &
     &                                        MTYP_SEL,MLEVS_SEL,       &
     &                  STIME_REJ,ETIME_REJ,NFIELDS_REJ,                &
     &                                        MTYP_REJ,MLEVS_REJ,       &
     &                  STIME_THI,ETIME_THI,NFIELDS_THI,                &
     &                      MTYP_THI,MLEVS_THI,IXXSTEP_THI,IYYSTEP_THI, &
     &                  MODIFY,SELECT,REJECT,THIN,OUTPUT_PACK_TYPE,     &
     &                  WIND_10M,WIND_10M_OROG,WIND_10M_SCALE,          &
     &                                                   PPUNIT_OROG,   &
     &                  PPUNIT1,PPUNIT2,ICODE,CMESSAGE)
      USE filenamelength_mod, ONLY :                                    & 
          filenamelength
      USE IO
      USE lookup_addresses
      IMPLICIT NONE
!LL  Stash variables
      INTEGER                                                           &
     &      LEN1_LOOKUP                                                 &
     &,     LEN2_LOOKUP                                                 &
     &,     LENBUF                                                      &
     &,     LEN_FIELD                                                   &
     &,     STIME_MOD                                                   &
     &,     ETIME_MOD                                                   &
     &,     NFIELDS_MOD                                                 &
     &,     MTYP_MOD(NFIELDS_mod)                                       &
     &,     MLEVS_MOD(NFIELDS_mod)                                      &
     &,     STIME_SEL                                                   &
     &,     ETIME_SEL                                                   &
     &,     NFIELDS_SEL                                                 &
     &,     MTYP_SEL(NFIELDS_sel)                                       &
     &,     MLEVS_SEL(NFIELDS_sel)                                      &
     &,     STIME_REJ                                                   &
     &,     ETIME_REJ                                                   &
     &,     NFIELDS_REJ                                                 &
     &,     MTYP_REJ(NFIELDS_rej)                                       &
     &,     MLEVS_REJ(NFIELDS_rej)                                      &
     &,     PPUNIT_OROG                                                 &
     &,     STIME_THI                                                   &
     &,     ETIME_THI                                                   &
     &,     NFIELDS_THI                                                 &
     &,     MTYP_THI(NFIELDS_THI)                                       &
     &,     MLEVS_THI(NFIELDS_THI)                                      &
     &,     IXXSTEP_THI(NFIELDS_THI)                                    &
     &,     IYYSTEP_THI(NFIELDS_THI)

      REAL                                                              &
     &      AMULT(NFIELDS_mod)                                          &
     &,     WIND_10M_OROG                                               &
     &,     WIND_10M_SCALE

      INTEGER                                                           &
     &      LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP),                            &
     &      LOOKNEW(LEN1_LOOKUP,LEN2_LOOKUP),                           &
     &      BUFOUT(LENBUF)

      CHARACTER(LEN=*) cmessage
      CHARACTER(LEN=filenamelength) ::  orogfile
      CHARACTER OUTPUT_PACK_TYPE*6
      LOGICAL LAST                !IN   indicates last record process
      LOGICAL                                                           &
     &      MODIFY                                                      &
     &,     REJECT                                                      &
     &,     SELECT                                                      &
     &,     WIND_10M                                                    &
     &,     THIN                                                        &
     &,     THIN_ALL

      INTEGER                                                           &
     &     LEN_FIXHD                                                    &
     &    ,PP_FIXHD(LEN_FIXHD)                                          &
     &    ,ICODE                                                        &
     &    ,PPUNIT1                                                      &
     &    ,PPUNIT2                                                      &
     &    ,DATA_ADDR                                                    &
                                  !    start address of data
     &    ,IEXTRA(10)                                                   &
                                  !IN  Used within FFREAD
     &    ,IER                                                          &
                                  !IN  error RETURN CODE from conversion
     &    ,ILABEL(45)                                                   &
                                  !IOUT  holds integet part of lookup
     &    ,NENT                                                         &
                                  !IN  NO. ENTRIES IN OLD LOOKUP
     &    ,ILABEL_OROG(45)
      REAL                                                              &
     &     ROOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                              &
     &    ,ROOKNEW(LEN1_LOOKUP,LEN2_LOOKUP)                             &
     &    ,RLABEL(19)                                                   &
                                  !OUT holds real part of LOOKUP
     &    ,FIELD(LEN_FIELD)                                             &
     &    ,RLABEL_OROG(19)                                              &
     &    ,MODEL_OROG(LENBUF)
!

      LOGICAL                                                           &
     &    PACKING                                                       &
     &,   READ                                                          &
     &,   CONVERT

      INTEGER                                                           &
     &     I                                                            &
                                  ! local counter
     &    ,J                                                            &
                                  ! local counter
     &    ,K                                                            &
                                  ! local counter
     &    ,IX                                                           &
                                  !
     &    ,IL                                                           &
                                  !
     &    ,BL                                                           &
                                  !
     &    ,TL                                                           &
                                  !
     &    ,IWL                                                          &
                                  !
     &    ,NLEV                                                         &
                                  !
     &    ,IWA                                                          &
                                  !
     &    ,IWB                                                          &
                                  !
     &    ,IENT                                                         &
                                  !
     &    ,IPROJ                                                        &
                                  !
     &    ,FCT                                                          &
                                  !
     &    ,ITYPE                                                        &
                                  !
     &    ,LEVEL                                                        &
                                  !
     &    ,IDIM                                                         &
                                  !
     &    ,LEN_LOOKUP                                                   &
                                  !
     &    ,LEN_IO                                                       &
                                  !
     &    ,LEN_BUF_WORDS                                                &
                                  !
     &    ,NUM_WORDS                                                    &
                                  !
     &    ,PACK_CODE                                                    &
     &    ,IXX                                                          &
                                  ! X dimension for THIN_FIELD
     &    ,IYY                                                          &
                                  ! Y dimension for THIN_FIELD
     &    ,IERR                   ! Error return from SETPOS
      REAL                                                              &
     &     A_IO                   !

      CHARACTER                                                         &
     &     PACK_TYPE(6)*6                                               &
     &    ,INPUT_PACK_TYPE*6


      PACK_TYPE(1)='NONE  '
      PACK_TYPE(2)='WGDOS '
      PACK_TYPE(3)='CRAY32'
      PACK_TYPE(4)='GRIB  '
      PACK_TYPE(5)='RUNLEN'
      PACK_TYPE(6)='      '

!
!    REMEMBER THAT BUFFER OUT STARTS AT ADDRESS 0 THUS LOOKUP GOES
!    FROM 0 to 262143 ie THE NEXT ADDRESS SHOULD BE IWA=262144 to
!    IWA=325119 then IWA=325120 to 388095 then 388096 etc
!

      LEN_LOOKUP=LEN1_LOOKUP*LEN2_LOOKUP

!----------------------- Section 4 ----------------------------------
!      Write to the PP file . First read in the  LOOKUP table.
!--------------------------------------------------------------------

      IWA=0
      CALL SETPOS(PPUNIT2,IWA,IERR)
      CALL BUFFIN(PPUNIT2,PP_FIXHD,LEN_FIXHD,LEN_IO,A_IO)
      IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
        CALL IOERROR('Buffer in fixed length header',A_IO,LEN_IO,       &
     &                                               LEN_FIXHD)
        ICODE=1
        CMESSAGE='REPLACE: I/O error'
        RETURN
      ENDIF
      IWL=PP_FIXHD(150)-1
      IWA=IWL
      DATA_ADDR = PP_FIXHD(160)

      CALL SETPOS(PPUNIT2,IWA,IERR)
      CALL BUFFIN(PPUNIT2,LOOKNEW,LEN_LOOKUP,LEN_IO,A_IO)
      IF(A_IO /= -1.0.OR.LEN_IO /= (PP_FIXHD(152)*PP_FIXHD(151)))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('Buffer in Lookup table   ',A_IO,LEN_IO,           &
     &                                    PP_FIXHD(152)*PP_FIXHD(151))
        ICODE=1
        CMESSAGE='Derived: I/O error in reading LOOKUP'
        RETURN
      ENDIF

      DO I=1,10
        IEXTRA(I)=0
      ENDDO

!     IF 10M WINDS TO BE FIXED GET MODEL OROGRAPHY FIELD FROM PP0
      IF(WIND_10M) THEN
        write(6,*) 'open unit',ppunit_orog
        CALL GET_FILE(PPUNIT_OROG,OROGFILE,filenamelength,ICODE)
        CALL FILE_OPEN(PPUNIT_OROG,OROGFILE,filenamelength,0,1,ICODE)
        IEXTRA(1) = 0
        IDIM  = LENBUF
        FCT   = 0
        IPROJ = LOOKUP(31,1)
        ITYPE = 73
        LEVEL = 9999
        write(6,*) ' read orography'
! DEPENDS ON: ffread
        CALL FFREAD(IPROJ,FCT,ITYPE,LEVEL,PPUNIT_OROG,MODEL_OROG,IDIM,  &
     &               ILABEL_OROG,RLABEL_OROG,IEXTRA,ICODE,CMESSAGE)
        write(6,*) 'close unit',ppunit_orog
        CLOSE(PPUNIT_OROG)
      ENDIF
!
!     loop through lookup read/write all fields
      IEXTRA(1) = 1      !DO NOT UNPACK
      DO IENT=1,NENT
        READ=.TRUE.
        CONVERT=.FALSE.
        IDIM=LENBUF
        FCT=LOOKUP(14,IENT)
        IPROJ=LOOKUP(31,IENT)
        ITYPE=LOOKUP(32,IENT)
        LEVEL=LOOKUP(33,IENT)
        PACK_CODE=MOD(LOOKUP(LBPACK,IENT),10)
        INPUT_PACK_TYPE=PACK_TYPE(PACK_CODE+1)
        IF(INPUT_PACK_TYPE /= OUTPUT_PACK_TYPE.AND.                     &
     &          PACK_CODE >  0) THEN ! leave unpacked data unpacked
          CONVERT=.TRUE.
        ENDIF
       WRITE(6,*)' pack code=',pack_code
       WRITE(6,*)input_pack_type,output_pack_type,convert
        IF(SELECT) THEN
          READ=.FALSE.
          IF(FCT >= STIME_SEL.AND.FCT <= ETIME_SEL) THEN
            DO J=1,NFIELDS_SEL
              IF(ITYPE == MTYP_SEL(J).AND.LEVEL == MLEVS_SEL(J)) THEN
                READ=.TRUE.
              ENDIF
            ENDDO
          ENDIF
        ENDIF
        IF(REJECT) THEN
          READ=.TRUE.
          IF(FCT >= STIME_REJ.AND.FCT <= ETIME_REJ) THEN
            DO J=1,NFIELDS_REJ
              IF(ITYPE == MTYP_REJ(J).AND.LEVEL == MLEVS_REJ(J)) THEN
                READ=.FALSE.
              ENDIF
            ENDDO
          ENDIF
        ENDIF

        WRITE(7,*) ' READ=',READ,IPROJ,FCT,ITYPE,LEVEL,PPUNIT1
        IF(READ) THEN
! DEPENDS ON: ffread
          CALL FFREAD(IPROJ,FCT,ITYPE,LEVEL,PPUNIT1,BUFOUT,IDIM,        &
     &                ILABEL,RLABEL,IEXTRA,ICODE,CMESSAGE)
          NUM_WORDS = ILABEL(15)
          LEN_BUF_WORDS = ILABEL(30)

          IF(STIME_THI == -9999) THEN
            THIN_ALL=.TRUE.
          ELSE
            THIN_ALL=.FALSE.
          ENDIF
          IF(THIN.OR.THIN_ALL) THEN
            IF((FCT >= STIME_THI.AND.FCT <= ETIME_THI)                  &
     &          .OR.THIN_ALL) THEN
              DO J=1,NFIELDS_THI
                IF((ITYPE == MTYP_THI(J).AND.LEVEL == MLEVS_THI(J))     &
     &              .OR.THIN_ALL) THEN
                  IYY = ILABEL(18)
                  IXX = ILABEL(19)
                  WRITE(7,*) ' THINNING FIELD,',ITYPE,LEVEL,FCT,        &
     &                                     IXXSTEP_THI(J),IYYSTEP_THI(J)
                  WRITE(6,*) ' THINNING FIELD,',ITYPE,LEVEL,FCT,        &
     &                                     IXXSTEP_THI(J),IYYSTEP_THI(J)
! DEPENDS ON: thin_field
                  CALL THIN_FIELD(BUFOUT,BUFOUT,NUM_WORDS,IXX,IYY,      &
     &                             IXXSTEP_THI(J),IYYSTEP_THI(J),       &
     &                               IDIM,PACK_CODE,RLABEL(18),         &
     &                               ILABEL(15),ICODE,CMESSAGE)
                  LEN_BUF_WORDS =((NUM_WORDS+511)/512)*512
                  ILABEL(15) = NUM_WORDS
                  ILABEL(30) = LEN_BUF_WORDS
                  ILABEL(18) = IYY
                  ILABEL(19) = IXX
                  rlabel(15) = rlabel(15) * IYYSTEP_THI(J)
                  rlabel(17) = rlabel(17) * IXXSTEP_THI(J)
                ENDIF
              ENDDO
            ENDIF
          ENDIF

          IF(MODIFY) THEN
            IF(FCT >= STIME_MOD.AND.FCT <= ETIME_MOD) THEN
              DO J=1,NFIELDS_MOD
                IF(ITYPE == MTYP_MOD(J).AND.LEVEL == MLEVS_MOD(J)) THEN
                  WRITE(7,*) ' SCALING FIELD,',ITYPE,LEVEL,FCT,AMULT(J)
                  WRITE(6,*) ' SCALING FIELD,',ITYPE,LEVEL,FCT,AMULT(J)
! DEPENDS ON: scale_field
                  CALL SCALE_FIELD(BUFOUT,BUFOUT,LEN_FIELD,AMULT(J),    &
     &                            ILABEL(15),IDIM,PACK_CODE,RLABEL(18), &
     &                            NUM_WORDS, ICODE, CMESSAGE)
                  LEN_BUF_WORDS =((NUM_WORDS+511)/512)*512
                  ILABEL(15) = NUM_WORDS
                  ILABEL(30) = LEN_BUF_WORDS
                ENDIF
              ENDDO
            ENDIF
          ENDIF
          IF(WIND_10M) THEN
            IF(ITYPE == 75.OR.ITYPE == 76) THEN
              write(6,*) 'call wind fix'
! DEPENDS ON: wind_10m_fix
              CALL WIND_10M_FIX(BUFOUT,BUFOUT,NUM_WORDS,                &
     &                          FCT,ITYPE,LEVEL,IPROJ,PPUNIT1,          &
     &                          WIND_10M_SCALE,WIND_10M_OROG,           &
     &                          MODEL_OROG,ILABEL_OROG,RLABEL_OROG,     &
     &                          IDIM,PACK_CODE,RLABEL(18))
              LEN_BUF_WORDS =((NUM_WORDS+511)/512)*512
              ILABEL(15) = NUM_WORDS
              ILABEL(30) = LEN_BUF_WORDS
            ENDIF
          ENDIF

          IF(CONVERT) THEN
! DEPENDS ON: conv_pack
            CALL CONV_PACK(ILABEL,RLABEL,PACK_CODE,                     &
     &                     INPUT_PACK_TYPE,OUTPUT_PACK_TYPE,            &
     &                     BUFOUT,IDIM,NUM_WORDS,                       &
     &                     PP_FIXHD,ICODE,CMESSAGE)
            LEN_BUF_WORDS =((NUM_WORDS+511)/512)*512
            ILABEL(15) = NUM_WORDS
            ILABEL(30) = LEN_BUF_WORDS
          ENDIF


! DEPENDS ON: fld_out
          CALL FLD_OUT(ICODE,CMESSAGE,BUFOUT,LENBUF,                    &
     &      LEN_BUF_WORDS,NUM_WORDS,                                    &
     &      PPUNIT2,LEN1_LOOKUP,LEN2_LOOKUP,LOOKNEW,LOOKNEW,            &
     &      ILABEL,RLABEL,IWL,DATA_ADDR)

!----------------------- Section 5 ----------------------------------
!          Output lookup table
!--------------------------------------------------------------------

          IWA=IWL
          CALL SETPOS(PPUNIT2,IWA,IERR)
          CALL BUFFOUT(PPUNIT2,LOOKNEW,LEN_LOOKUP,LEN_IO,A_IO)
!
          IF(A_IO /= -1.0.OR.LEN_IO /=                                  &
     &                              (PP_FIXHD(152)*PP_FIXHD(151)))THEN
! DEPENDS ON: ioerror
            CALL IOERROR('Buffer in fixed length header',A_IO,LEN_IO,   &
     &                    PP_FIXHD(151)*PP_FIXHD(152))
            ICODE=1
            CMESSAGE='Derived: I/O error in writing LOOKUP'
            RETURN
          ENDIF
        ENDIF
      ENDDO


 9999 CONTINUE
      RETURN
      END SUBROUTINE FIELDS
