! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL
!LL  Subroutine: READ_WRITE 
!LL
!LL  Purpose: To read a   direct access PP file  and convert it to a
!LL  sequential file read to be passed across to the IBM
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: C41
!LL
!LL  Project task: C4
!LL
!LL  External documentation: UM Documentation paper C4
!LL
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Small execs
      SUBROUTINE READ_WRITE_fieldcos(PP_LEN2_LOOKUP,LEN1_LOOKUP,        &
     &                      DATA_ADD,PP_FIXHD,                          &
     &                      IWA,UNPACK,FORMAT_OUT,PPUNIT,COS_PPUNIT,    &
     &                      IEXTRA,OPER,MASS,ICODE,CMESSAGE,LCAL360)
      USE filenamelength_mod, ONLY :                                    & 
          filenamelength
      USE IO
      USE lookup_addresses

      IMPLICIT NONE
!     Arguments
      CHARACTER                                                         &
     &     FORMAT_OUT*6                                                 &
                                ! IN    IBM/IEEE format for output
     &    ,CMESSAGE*80
                                ! OUT   error messages
      CHARACTER(LEN=filenamelength) :: outfile         
                                ! OUT   pathname of output file
      LOGICAL                                                           &
     &     UNPACK                                                       &
                                ! IN    indicates whether to unpack
     &    ,OPER                                                         &
                                ! IN    indicates whether operational
     &    ,MASS                                                         &
                                ! IN    indicates whether MASS
     &    ,LCAL360
      INTEGER                                                           &
     &     LEN1_LOOKUP                                                  &
                                ! IN    1st dimension of LOOKUP
     &    ,PP_LEN2_LOOKUP                                               &
                                ! IN    2nd dimension of LOOKUP
     &    ,PPUNIT                                                       &
                                ! IN    unit no of required fieldsfile
     &    ,COS_PPUNIT                                                   &
                                ! IN    unit no of COS output file
     &    ,DATA_ADD                                                     &
                                ! IN    word address of the data.
     &    ,IEXTRA(10)                                                   &
                                ! IN    Controls READFF
     &    ,IWA                                                          &
                                ! IN    Word address in call SETPOS
     &    ,PP_FIXHD(*)                                                  &
                                ! IN    PPfile fixed header
     &    ,ICODE                ! OUT   error code


!     arguments for called routines
      LOGICAL                                                           &
     &     MODEL_FLAG                                                   &
                                ! flag - set to true if model dump
     &    ,LAST                 ! indicates last record process
      INTEGER                                                           &
     &     LOOKUP(LEN1_LOOKUP,PP_LEN2_LOOKUP)                           &
                                               ! integer lookup
     &    ,NUM_VALUES                                                   &
                                ! No of data points in a field
     &    ,IDIM                                                         &
                                ! NUM_VALUES rounded to an even no
!                               !  used to dimension the output array
     &    ,IEXTRAW                                                      &
                                ! The number of words of "extra" data.
     &    ,ENTRY_NO                                                     &
                                ! lookkup entry no of the Field.
     &    ,LEN_IO                                                       &
                                ! actual no of words transferred by IO.
     &    ,LEN_IO_EXPECTED      ! expected no of words transferred by IO
      REAL                                                              &
     &     A_IO                 ! status returned by BUFFIN
!*---------------------------------------------------------------------
!    LOCAL VARIABLES
      INTEGER                                                           &
     &     I                                                            &
                                ! local counter
     &    ,J                                                            &
                                ! local counter
     &    ,IX                                                           &
                                ! used as a dummy variable in UNIT
     &    ,ICOUNT                                                       &
                                ! Counter
     &    ,NENT                                                         &
                                ! No of entries in the printfile
     &    ,TOTAL_WORDS          ! Total number of words output GRIB
                                ! option only

!L -----------Read in the LOOKUP table if first time thro------------
      CALL SETPOS(PPUNIT,IWA,ICODE)
      LEN_IO_EXPECTED=PP_LEN2_LOOKUP*LEN1_LOOKUP
      CALL BUFFIN(PPUNIT,LOOKUP,LEN_IO_EXPECTED,LEN_IO,A_IO)
      IF(A_IO /= -1.0.OR.LEN_IO /= LEN_IO_EXPECTED) THEN
! DEPENDS ON: ioerror
        CALL IOERROR('Buffer in lookup table   ',A_IO,LEN_IO,           &
     &            LEN_IO_EXPECTED )
        CMESSAGE='READ_W : I/O error reading LOOKUP TABLE  '
        ICODE=3
          WRITE(6,*)' I/O error reading LOOKUP TABLE'
          RETURN
      ENDIF

!L -----Having read the LOOKUP table Open the output COS File--------
      CALL GET_FILE(COS_PPUNIT,OUTFILE,filenamelength,ICODE)
      OPEN(UNIT=COS_PPUNIT,FILE=OUTFILE,FORM='UNFORMATTED')
      NENT=0

!L -----Calculate the number of fields in the File-------------------
      DO I=1,PP_LEN2_LOOKUP
        IF(LOOKUP(LBROW,I) /= -99) THEN
          NENT=NENT+1
        ELSE
          GOTO 2
        ENDIF
      ENDDO
    2 CONTINUE
      WRITE(6,*)' THE NUMBER OF FIELDS IN THE INPUT FILE IS ', NENT
      LAST=.FALSE.
!--------------------------------------------------------------------
! Note LBROW=18,LBNPT=19
! For a DUMP LBLREC will hold original no of data points.
! LBNREC will be set to zero.
!
! For a PP_file LBLREC will hold the no of CRAY words needed to hold
! the data. The original field size will be rows*columns.
! If the data is not packed then LBLREC=LBROW*LBNPT+LBEXT, where
! LBEXT will be greater than 0 for timeseries (which are never packed).
!  !! WARNING LBEXT - may be -32768 MISSING VALUE !!
!---------------------------------------------------------------------
!
!L -----Set MODEL_FLAG and reset UNPACK if DUMP ---------------------
      IF(PP_FIXHD(5) /= 3) THEN
        MODEL_FLAG=.TRUE.       ! Model dump
        UNPACK= .TRUE.          ! cray 32 bit packed data unpacked
        WRITE(6,*)'Model dump - UNPACK set TRUE '
      ELSE
        MODEL_FLAG=.FALSE.      ! Fieldsfile
      ENDIF
      IF(.NOT.UNPACK) IEXTRA(1)=1  ! DATA LEFT PACKED

!L -----Loop thro all the entries within the field ------------------
      DO I=1,NENT
        IF(I == NENT) LAST=.TRUE.
        IF(MODEL_FLAG) THEN
          NUM_VALUES=LOOKUP(LBLREC,I)    ! NCOLS*NROWS
        ELSE
          NUM_VALUES=LOOKUP(LBROW,I)*LOOKUP(LBNPT,I)+LOOKUP(LBEXT,I)
        ENDIF
        IEXTRAW=0
        IF(LOOKUP(LBEXT,I) >  0) THEN ! got some extra data
          IEXTRAW=LOOKUP(LBEXT,I)
! check to see that we don't have packing if we have extra data....
          IF(LOOKUP(LBROW,I)*LOOKUP(LBNPT,I)+LOOKUP(LBEXT,I)  /=        &
     &      LOOKUP(LBLREC,I)) THEN
            CMESSAGE='READ_WRT : Packing of extra data not supported'
            ICODE=1
            RETURN
          ENDIF
        ENDIF
        IDIM=((NUM_VALUES+1)/2)*2 ! Round to ensur an integer for IBM

!L---------------------------------------------------------------------
!L If packed simply read in the field ie LBLREC words for PP_type
!L files  & for Dump type read LBLREC/2 if packed and LBLREC if not.
!L All packed data is assumed real. If the data is to be un-packed
!L then it is un-packed into an array size IDIM. IDIM is NROWS*NCOLS+ext
!L rounded up to ensure it is even. If the data is not packed then it
!L could be REAL,LOGICAL or INTEGER .
!L--------------------------------------------------------------------

        ICODE=0
        ENTRY_NO=I
! Passing variable MASS only to CRAY_IEEE since it should be FALSE for all the
! other FORMAT_OUT types.  The check is done just after reading in namelists.
        IF(FORMAT_OUT == 'IBM') THEN
! DEPENDS ON: cray_ibm
          CALL CRAY_IBM(IDIM,NUM_VALUES,PPUNIT,                         &
     &                  LEN1_LOOKUP,PP_LEN2_LOOKUP,PP_FIXHD,LOOKUP,     &
     &                  LOOKUP,ENTRY_NO,DATA_ADD,MODEL_FLAG,            &
     &                  COS_PPUNIT,IEXTRA,IEXTRAW,LAST,OPER,            &
     &                  ICODE,CMESSAGE,LCAL360)
        ELSEIF(FORMAT_OUT == 'IEEE') THEN
! DEPENDS ON: cray_ieee
          CALL CRAY_IEEE(IDIM,NUM_VALUES,PPUNIT,                        &
     &                   LEN1_LOOKUP,PP_LEN2_LOOKUP,PP_FIXHD,LOOKUP,    &
     &                   LOOKUP,ENTRY_NO,DATA_ADD,MODEL_FLAG,           &
     &                   COS_PPUNIT,IEXTRA,IEXTRAW,LAST,OPER,MASS,      &
     &                   ICODE,CMESSAGE,LCAL360)
        ELSE
          ICODE=1
          CMESSAGE= ' OUTPUT FORMAT NOT YET AVAILABLE '
        ENDIF
        IF(ICODE /= 0) THEN
          RETURN
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE READ_WRITE_fieldcos
