! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL  Routine: READFF
!LL
!LL  Purpose: To read a   direct access PP file.
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered:
!LL
!LL  Project task:
!LL
!LL  External documentation:
!LL
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Small execs

      SUBROUTINE READFF_fieldcos(PPUNIT,FIELD,IDIM,ENTRY_NO,            &
     &ILABEL,RLABEL,IEXTRA,PP_LEN2_LOOKUP,LEN1_LOOKUP,                  &
     &PP_FIXHD,LOOKUP,ROOKUP,DATA_ADD,                                  &
     &MODEL_FLAG,MAX_LEN_ILABEL,MAX_LEN_RLABEL,                         &
     &LEN_ILABEL,LEN_RLABEL,                                            &
     &ICODE,CMESSAGE)

      USE lookup_addresses

      IMPLICIT NONE
!     arguments
      CHARACTER                                                         &
     &     CMESSAGE*(*)           !OUT error message
      LOGICAL                                                           &
     &     MODEL_FLAG             !IN  True => Dump False =>Fieldsfile
      INTEGER                                                           &
     &     LEN1_LOOKUP                                                  &
                                  !IN  first dimension of the lookup
     &    ,PP_LEN2_LOOKUP                                               &
                                  !IN  secnd dimension of the lookup
     &    ,PPUNIT                                                       &
                                  !IN  unit no of required fieldsfile
     &    ,IDIM                                                         &
                                  !IN  dimension of FIELD
     &    ,MAX_LEN_RLABEL                                               &
                                  !IN  max sixe of RLABEL
     &    ,MAX_LEN_ILABEL                                               &
                                  !IN  max sixe of ILABEL
     &    ,IEXTRA(10)                                                   &
                                  !IN  spare for future use
     &    ,DATA_ADD                                                     &
                                  !IN  The word address of the data.
     &    ,ENTRY_NO                                                     &
                                  !IN  Lookup entry no of the Field.
     &    ,PP_FIXHD(*)                                                  &
                                  !IN  PPfile fixed header
     &    ,LOOKUP(LEN1_LOOKUP,PP_LEN2_LOOKUP)                           &
                                              !IN integer lookup
     &    ,LEN_RLABEL                                                   &
                                  !OUT actual size of RLABEL
     &    ,LEN_ILABEL                                                   &
                                  !OUT actual size of ILABEL
     &    ,ILABEL(MAX_LEN_ILABEL)                                       &
                                  !OUT integer part of LOOKUP
     &    ,ICODE                  !OUT error code
      REAL                                                              &
     &     FIELD(IDIM)                                                  &
                                  !OUT array holding final output data.
     &    ,ROOKUP(LEN1_LOOKUP,PP_LEN2_LOOKUP)                           &
                                              !IN real lookup
     &    ,RLABEL(MAX_LEN_RLABEL) !OUT real part of LOOKUP

!     arguments for called routines
      INTEGER                                                           &
     &     PACK_TYPE                                                    &
                                  ! packing type N1 of LBPACK
     &    ,NUM_CRAY_WORDS                                               &
                                  ! number of words for field
     &    ,NVALS                                                        &
                                  ! number of points in a data field
     &    ,IWA                    ! Word address in call SETPOS
!*---------------------------------------------------------------------
!     LOCAL VARIABLES
      INTEGER                                                           &
     &     I                                                            &
                                  ! Local counter
     &    ,J                                                            &
                                  ! Local counter
     &    ,LENGTH_OF_DATA                                               &
                                  ! Length of a particular field
     &    ,ADDR                                                         &
                                  ! Address of a field in the data store
     &    ,IN_LBVC                                                      &
                                  ! Local copy of LBVC required to searc
     &    ,NUM_IBM_WORDS                                                &
                                  ! No of IBM words used to hold the dat
     &    ,POS_RLABEL                                                   &
                                  ! position of first REAL in PPhdr
     &    ,PACK_TYPE_I                                                  &
                                  ! packing type N1 of LBPACK
     &    ,DATA_COMP                                                    &
                                  ! data compression code
     &    ,DATA_COMP_DEF                                                &
                                  ! data compression definition
     &    ,NUMBER_FORMAT          ! number format
      REAL                                                              &
     &     AMDI                   ! Missing data indicator for lookup

! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

      AMDI=ROOKUP(BMDI,ENTRY_NO)
      IF (AMDI /= RMDI) WRITE(6,*)' NONE STANDARD MISSING DATA USED'

!
!     CALL PR_LOOK(LOOKUP(1,1),ROOKUP(1,1),ENTRY_NO)
!
!     decode LBPACK
      PACK_TYPE = MOD(LOOKUP(LBPACK,ENTRY_NO),10)
      DATA_COMP = MOD(LOOKUP(LBPACK,ENTRY_NO),100) - PACK_TYPE
      DATA_COMP_DEF = MOD(LOOKUP(LBPACK,ENTRY_NO),1000)                 &
     &                                      -DATA_COMP -PACK_TYPE
      NUMBER_FORMAT = LOOKUP(LBPACK,ENTRY_NO)/1000
!----------------------------------------------------------------------
!=== Reading a model type dump =======================================
!    A model dump has no direct addressing only relative.
!
      IF(MODEL_FLAG) THEN
! Old Format dumpfiles
        if((lookup(lbnrec,entry_no) == 0) .or.                          &
! Prog lookups in dump before vn3.2:
     &    ((lookup(lbnrec,entry_no) == imdi) .and.                      &
     &                             (pp_fixhd(12) <= 301))) then

        IF(PACK_TYPE == 2) THEN            ! 32 bit packing.
          NUM_CRAY_WORDS=(LOOKUP(LBLREC,ENTRY_NO)+1)/2
        ELSEIF(PACK_TYPE >  0) THEN
          NUM_CRAY_WORDS=LOOKUP(LBLREC,ENTRY_NO)/2
        ELSE
          NUM_CRAY_WORDS=LOOKUP(LBLREC,ENTRY_NO)
        ENDIF
        NVALS=LOOKUP(LBLREC,ENTRY_NO) ! No of data points
        ADDR=DATA_ADD
        IF(ENTRY_NO >  1) THEN
          DO I=1,ENTRY_NO-1
            PACK_TYPE_I = MOD(LOOKUP(LBPACK,I),10)
            IF(PACK_TYPE_I == 2) THEN ! 32 Bit packed
              LENGTH_OF_DATA=(LOOKUP(LBLREC,I)+1)/2
            ELSE
              LENGTH_OF_DATA=LOOKUP(LBLREC,I)
            ENDIF
            ADDR=ADDR+LENGTH_OF_DATA
          ENDDO
        ELSE       !  If the first entry.
          ADDR=DATA_ADD  !
          IF(PACK_TYPE == 2) THEN ! 32 Bit packed
            LENGTH_OF_DATA=(LOOKUP(LBLREC,1)+1)/2
          ELSE
            LENGTH_OF_DATA=LOOKUP(LBLREC,1)
          ENDIF
          WRITE(6,*)'  LENGTH_OF_DATA  ',LENGTH_OF_DATA
        ENDIF
        IWA=ADDR  ! Not -1 as this is already done in dump
      Else
! New format Dumpfiles (vn4.4 onwards)

        If(pack_type == 2) then            ! 32 bit packing.
          num_cray_words=(lookup(lblrec,entry_no)+1)/2
        Elseif(pack_type >  0) then
          num_cray_words=lookup(lblrec,entry_no)/2
        Else
          num_cray_words=lookup(lblrec,entry_no)
        Endif
        iwa = lookup(lbegin,entry_no)
        nvals = lookup(lbrow,entry_no) * lookup(lbnpt,entry_no)
      Endif
      ELSE
!=== Reading a PP type file.==========================================
        NUM_CRAY_WORDS=LOOKUP(LBLREC,ENTRY_NO) ! PP type file
        IWA=LOOKUP(LBEGIN,ENTRY_NO)
        NVALS=LOOKUP(LBROW,ENTRY_NO)*LOOKUP(LBNPT,ENTRY_NO)             &
     &         +LOOKUP(LBEXT,ENTRY_NO)
      ENDIF
!==============================================================
!       WRITE(6,107) ENTRY_NO,NUM_CRAY_WORDS,NVALS
  107 FORMAT(' ENTRY NO=',I5,'NUM_CRAY_WORDS= ',I6,'NVALS=',I6)
        IF(IDIM <  NUM_CRAY_WORDS) THEN
          ICODE=NUM_CRAY_WORDS
          CMESSAGE='READFF  Idim to small ICODE holds correct value'
          GOTO 9999
        ENDIF
      ICODE=0
!     RETURN
! DEPENDS ON: read_rec_fieldcos
      CALL READ_REC_FIELDCOS(FIELD,NUM_CRAY_WORDS,IWA,PPUNIT,ICODE,CMESSAGE)
 2212 FORMAT('  FIELDS FILE NUMBER ',I2,'  ON UNIT',I2,2X,'BEING READ')
      NUM_IBM_WORDS=NUM_CRAY_WORDS*2

      WRITE(7,106) ENTRY_NO,                                            &
                                                ! Field No
     &             LOOKUP(LBTYP,ENTRY_NO),                              &
                                                ! M08 Type
     &             LOOKUP(LBFC,ENTRY_NO),                               &
                                                ! PP Field Code
     &             LOOKUP(ITEM_CODE,ENTRY_NO),                          &
                                                ! Stash Code
     &             LOOKUP(LBLEV,ENTRY_NO),                              &
                                                ! M08 Level
     &             LOOKUP(LBFT,ENTRY_NO),                               &
                                                ! Forecast period
     &             LOOKUP(LBPROJ,ENTRY_NO),                             &
                                                ! M08 Projection no
     &             NUM_IBM_WORDS,                                       &
     &             NVALS,                                               &
     &             PACK_TYPE                    ! Packing Code

  106 FORMAT(' Field No ',I4,' M08/PP/Stash Code ',I3,I5,I6,            &
     &       ' Level ',I5,' Fcst ',I5,' Proj ',I3,                      &
     &       ' NWords=',I6,' NVals=',I5,' Pack Type=',I2)

        IF(ICODE == 0) THEN
          POS_RLABEL=MOD(LOOKUP(LBREL,ENTRY_NO),100)

          ! Treat lookup(45) as an integer to preserve submodel
          ! identifier in PP fields transferred between Cray and IBM.
          POS_RLABEL=46


          LEN_RLABEL=1+LEN1_LOOKUP-POS_RLABEL
          LEN_ILABEL=LEN1_LOOKUP-LEN_RLABEL
          DO I=1,LEN_ILABEL
            ILABEL(I)=LOOKUP(I,ENTRY_NO)
          ENDDO

!         check for valid release number
          if(ilabel(lbrel) <  1) then
            WRITE(6,*)' resetting LBREL from',ilabel(lbrel),' to 3'
            ilabel(lbrel)=3
          endif

!  test of header with position of reals

!         ilabel(lbrel)= 3*1000 + pos_rlabel
!         ilabel(lbrel)= 3
!         ilabel(lbsrce)=pos_rlabel

!  end of test

          DO I=1,LEN_RLABEL
            RLABEL(I)=ROOKUP(I+POS_RLABEL-1,ENTRY_NO)
          ENDDO
        ENDIF
!=======================================================================
! At this point FIELD holds the data either PACKED or UN-PACKED
! Is the packing indicator set and is un-packing required? If so then
! the data is temp un-packed into a work ARRAY of length IDIM
        IF(PACK_TYPE >  0) THEN                ! Is the field packed.
          IF(IEXTRA(1) == 0) THEN  ! unpacking is required
! DEPENDS ON: un_pack_fieldcos
            CALL UN_PACK_FIELDCOS(PACK_TYPE,IDIM,FIELD,NUM_CRAY_WORDS,  &
     &                   ILABEL,LEN_ILABEL,aMDI,PP_FIXHD,ICODE,CMESSAGE)
!           WRITE(7,*) ' NOW UNPACKED INTO ',ILABEL(LBLREC),' WORDS'
          ENDIF
        ELSEIF(LOOKUP(DATA_TYPE,ENTRY_NO) == 3) THEN !Fld is logical
! DEPENDS ON: logical_to_real_fieldcos
          CALL LOGICAL_TO_REAL_FIELDCOS(IDIM,FIELD,FIELD,NVALS,         &
     &                         ILABEL,ICODE,CMESSAGE)
        ELSEIF(LOOKUP(DATA_TYPE,ENTRY_NO) == 2) THEN !Fld is integer
! DEPENDS ON: integer_to_real_fieldcos
          CALL INTEGER_TO_REAL_FIELDCOS(IDIM,FIELD,FIELD,NVALS,         &
     &                         ILABEL,ICODE,CMESSAGE)
        ENDIF
!=======================================================================
 9999 CONTINUE
  100 FORMAT(//,32X,'   ARRAY        ',//,32(16F5.0/))
  101 FORMAT(//,32X,'   LOOKUP       ',//,32(16I5/))
  103 FORMAT('   LENIN  ',I12)
      RETURN
      END SUBROUTINE READFF_fieldcos
