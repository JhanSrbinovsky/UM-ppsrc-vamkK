! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!LL  Subroutine: CRAY_IEEE
!LL
!LL  Purpose: To read a   direct access PP file  and convert it to a
!LL  sequential file read to be passed in IEEE format
!LL
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: C41
!LL
!LL  Project task: C4
!LL
!LL  External documentation:
!LL
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Small execs

      SUBROUTINE CRAY_IEEE(IDIM,NUM_VALUES,PPUNIT,                      &
                   LEN1_LOOKUP,PP_LEN2_LOOKUP,PP_FIXHD,LOOKUP,          &
                   ROOKUP,ENTRY_NO,DATA_ADD,MODEL_FLAG,                 &
                   COS_PPUNIT,IEXTRA,IEXTRAW,LAST,OPER,MASS,            &
                   ICODE,CMESSAGE,LCAL360)
      USE lookup_addresses

      IMPLICIT NONE
!     Arguments
      CHARACTER                                                         &
           CMESSAGE*(*)           !OUT error messages
      LOGICAL                                                           &
           LAST                                                         &
                                  !IN indicates last record process
          ,OPER                                                         &
                                  !IN indicates whether operational
          ,MASS                                                         &
                                  !IN indicates whether using for MASS
          ,MODEL_FLAG                                                   &
                                  !IN True => dumps, False => fieldsfile
          ,LCAL360
      INTEGER                                                           &
           PPUNIT                                                       &
                                  !IN unit no of required fieldsfile
          ,COS_PPUNIT                                                   &
                                  !IN unit no of COS output file
          ,NUM_VALUES                                                   &
                                  !IN No of data points NROWS*NCOLS
          ,IDIM                                                         &
                                  !IN NUM_VALUES rounded to an even no
!                                 !  used to dimension The output array
          ,DATA_ADD                                                     &
                                  !IN The word address of the data.
     &    ,LEN1_LOOKUP                                                  &
                                  !IN First dimension of the lookup
     &    ,PP_LEN2_LOOKUP                                               &
                                  !IN Size of the LOOKUP on the file
     &    ,IEXTRA(10)                                                   &
                                  !IN Used within READFF
     &    ,IEXTRAW                                                      &
                                  !IN no of words of extra data.
     &    ,ENTRY_NO                                                     &
                                  !IN Lookup entry no of the Field.
     &    ,PP_FIXHD(*)                                                  &
                                  !IN PPfile fixed header
     &    ,LOOKUP(LEN1_LOOKUP,PP_LEN2_LOOKUP)                           &
                                               !IN integer lookup
     &    ,ICODE                  !OUT error code
      REAL                                                              &
     &     ROOKUP(LEN1_LOOKUP,PP_LEN2_LOOKUP)  !IN Real lookup
!*---------------------------------------------------------------------
!     Called routines
      EXTERNAL INT_FROM_REAL
      INTEGER IEEE2IEG
      INTEGER INT_FROM_REAL
!*---------------------------------------------------------------------
!     arguments for called routines
      INTEGER                                                           &
     &     MAX_LEN_ILABEL                                               &
     &    ,MAX_LEN_RLABEL
      PARAMETER (MAX_LEN_ILABEL=45,MAX_LEN_RLABEL=32)
      INTEGER                                                           &
     &     END_YEAR                                                     &
                           ! )
     &    ,END_MONTH                                                    &
                           ! )
     &    ,END_DAY                                                      &
                           ! )  arguments
     &    ,END_HOUR                                                     &
                           ! )
     &    ,END_MINUTE                                                   &
                           ! )     for
     &    ,END_SECOND                                                   &
                           ! )
     &    ,END_DAY_NUMBER                                               &
                           ! )
     &    ,END_TIME_DAYS                                                &
                           ! )
     &    ,END_TIME_SECS                                                &
                           ! )  date/time
     &    ,START_TIME_SECS                                              &
                           ! )
     &    ,START_TIME_DAYS                                              &
                           ! )
     &    ,DATA_YEAR                                                    &
                           ! )  conversion
     &    ,DATA_MONTH                                                   &
                           ! )
     &    ,DATA_DAY                                                     &
                           ! )     when
     &    ,DATA_HOUR                                                    &
                           ! )
     &    ,DATA_MINUTE                                                  &
                           ! )  OPER is TRUE
     &    ,DATA_SECOND                                                  &
                           ! )
     &    ,DATA_DAY_NUMBER                                              &
                           ! )
     &    ,ADDR                                                         &
                           ! address in fld, used to process extra data
     &    ,IEEE_ADDR                                                    &
                           ! address in ibm fld where extra data going.
     &    ,BIT_OFF                                                      &
                           ! what bit offset are we using
!                            (32 for odd, 0 for even values of addr)
     &    ,IER                                                          &
                           ! error RETURN CODE from conversion
     &    ,IV                                                           &
                           ! value of integer code for vectors
     &    ,LEN_ILABEL                                                   &
                           ! number of values in ILABEL
     &    ,LEN_RLABEL                                                   &
                           ! number of values in RLABEL
     &    ,DATA_VALUES                                                  &
                           ! number of values in real extra data
     &    ,ILABEL(MAX_LEN_ILABEL)                                       &
                                        ! holds integer part of LOOKUP
     &    ,IEEE_LABEL((LEN1_LOOKUP+1)/2)! holds IEEEconversion of LABEL

      REAL                                                              &
     &     FIELD(IDIM)                                                  &
                                  ! array holding data
     &    ,IEEE_FIELD(IDIM/2)                                           &
                                  ! array holding IEEE data
     &    ,RLABEL(MAX_LEN_RLABEL) ! holds real part of LOOKUP

!*---------------------------------------------------------------------
!    LOCAL VARIABLES
      INTEGER                                                           &
     &     I                                                            &
                          ! local counter
     &    ,PACK_TYPE                                                    &
                          ! packing type N1 of LBPACK
     &    ,DATA_COMP                                                    &
                          ! data compression code
     &    ,DATA_COMP_DEF                                                &
                          ! data compression definition
     &    ,NUMBER_FORMAT                                                &
                          ! number format
     &    ,FCST_PRD

      LOGICAL PACKED      ! indicates whether the data is packed



      DO I=1,IDIM        ! make sure FIELD is initialised. An odd
      FIELD(I)=0.0       ! number of points might upset conversion
      END DO
      PACKED=.FALSE.

!L access the Fields File.
! DEPENDS ON: readff_fieldcos
      CALL READFF_fieldcos(PPUNIT,FIELD,IDIM,ENTRY_NO,                  &
     &ILABEL,RLABEL,IEXTRA,PP_LEN2_LOOKUP,LEN1_LOOKUP,                  &
     &PP_FIXHD,LOOKUP,ROOKUP,DATA_ADD,                                  &
     &MODEL_FLAG,MAX_LEN_ILABEL,MAX_LEN_RLABEL,                         &
     &LEN_ILABEL,LEN_RLABEL,                                            &
     &ICODE,CMESSAGE)
!
      IF(ICODE /= 0) RETURN

!-----------------------------------------------------------------

! The data has now been read in and has 1) Been read in packed
! and left packed or 2) read in as packed and then un-packed or
! 3) The data was never packed at all.  If packed FIELD will have
! LBLREC/2 values if a DUMP and LBLREC values if a PP_FILE. If
! the data is not packed FIELD will have the no of data points
! length LBROW*LBNPT+LBEXT if a pp_file and LBLREC if a dump file.
!
! For a dump LBLREC will hold original no of data points.  For a
! pp_file LBLREC will hold the no of CRAY words needed to hold
! the data (if un-packed also no of data points)
!
! The value returned in ILABEL(LBLREC) may have to change because
! IEEE only has a 32 bit word length compared to the CRAY's 64
! bit word length. In IEEE ILABEL(LBLREC) will be no of 32 bit
! words needed to hold the data . If the data is not packed (or
! it has been un-packed) then this will be the no of data points.
! If the data is left packed the value of ILABEL(LBLREC)
! will have to be doubled as the no of 32bit words needed to
! hold the data will twice that on the CRAY.

! On output the data will either have been converted to IEEE
! numbers and stored in IEEE_FIELD or left packed in FIELD,but with
! 32 bit ieee numbers substituted into minimum values.
! If packed then LBLREC/2  words of FIELD are written as LBLREC is
! now the no of 32 bit words. If un-packed IEEE_FIELD which has size
! IDIM/2 (or NUM_VALUES/2) is written as it is.

!-----------------------------------------------------------------
!     decode LBPACK
      PACK_TYPE = MOD(ILABEL(LBPACK),10)
      DATA_COMP = MOD(ILABEL(LBPACK),100) - PACK_TYPE
      DATA_COMP_DEF = MOD(ILABEL(LBPACK),1000) -DATA_COMP -PACK_TYPE
      NUMBER_FORMAT = ILABEL(LBPACK)/1000

      IF(PACK_TYPE >  0) PACKED=.TRUE.
      IF(PACKED) THEN              ! Data left in packed form. Number of
        ILABEL(LBLREC)=ILABEL(LBLREC)*2 ! IEEE words needed is 2*CRAY
      ENDIF
! verify that don't have extra data and packing at once
      IF (IEXTRAW >  0.AND.PACKED) THEN
        CMESSAGE='FIELDCOS: Extra data with packing not supported'
        ICODE=1
        RETURN
      ENDIF

!L Convert ILABEL to IBM(Hitachi) integers.
! For either an accumulation or time mean (ie LBTIM /= 0) the start &
! end time are in a different order to the data and veri time for a
! snap shot type field. This anomaly has to be catered for operational
! use. Thus the PP package will not work properly on accum/time mn field
! for operational Fields files.
      IF(ILABEL(lbtim) /= 11.AND.OPER) THEN
        IF(ILABEL(LBYR) >  0)THEN
!       re -calculate the data time from the end time  and fcst period
!     First calculate the no of seconds from day 0
! if initial data time is missing, cannot recalculate data time
! for accum/time mean fields, hence, the PP Package would fail.
! Therefore, if the anomaly arises and temporal data is missing,
! this subroutine will be aborted and data from this field will
! not be written to COS file.
        END_YEAR=ILABEL(LBYRD)
        END_MONTH=ILABEL(LBMOND)
        END_DAY=ILABEL(LBDATD)
        END_HOUR=ILABEL(LBHRD)
        END_MINUTE=ILABEL(LBMIND)
        END_SECOND=ILABEL(LBSECD)
        FCST_PRD=ILABEL(LBFT)
!       WRITE(6,*)' START YR/MO/DA/HR/MIN/SEC BEFORE ',ILABEL(1), &
!       ILABEL(2),ILABEL(3),ILABEL(4),ILABEL(5),ILABEL(6)
!       WRITE(6,*)' END   YR/MO/DA/HR/MIN/SEC BEFORE ',ILABEL(7), &
!       ILABEL(8),ILABEL(9),ILABEL(10),ILABEL(11),ILABEL(12)
!       WRITE(6,*)' FCST_PRD BEFORE  ',FCST_PRD
! DEPENDS ON: time2sec
        CALL TIME2SEC (END_YEAR,END_MONTH,END_DAY,END_HOUR,             &
     &               END_MINUTE,END_SECOND,0,0,                         &
     &               END_TIME_DAYS,END_TIME_SECS,LCAL360)

!   Subtract forecast hours from end time in (days,seconds)

! DEPENDS ON: time_df
       CALL TIME_DF(END_TIME_DAYS,END_TIME_SECS,0,-FCST_PRD*3600,       &
     &              START_TIME_DAYS,START_TIME_SECS)

!     Go back and re-calculate Year/Month/Day/Hour/Sec.
! DEPENDS ON: sec2time
       CALL SEC2TIME(0,0,START_TIME_DAYS,START_TIME_SECS,               &
     &             DATA_YEAR,DATA_MONTH,DATA_DAY,                       &
     &             DATA_HOUR,DATA_MINUTE,DATA_SECOND,DATA_DAY_NUMBER,   &
     &             LCAL360)
        ILABEL(LBYRD)=DATA_YEAR
        ILABEL(LBMOND)=DATA_MONTH
        ILABEL(LBDATD)=DATA_DAY
        ILABEL(LBHRD)=DATA_HOUR
        ILABEL(LBMIND)=DATA_MINUTE
        ILABEL(LBSECD)=DATA_SECOND
        ILABEL(LBYR)=END_YEAR
        ILABEL(LBMON)=END_MONTH
        ILABEL(LBDAT)=END_DAY
        ILABEL(LBHR)=END_HOUR
        ILABEL(LBMIN)=END_MINUTE
        ILABEL(LBSEC)=END_SECOND
!       WRITE(6,*)' -----------------------------------------------'
!       WRITE(6,*)' Veri  YR/MO/DA/HR/MIN/SEC AFTER ',ILABEL(1), &
!       ILABEL(2),ILABEL(3),ILABEL(4),ILABEL(5),ILABEL(6)
!       WRITE(6,*)' Data  YR/MO/DA/HR/MIN/SEC AFTER ',ILABEL(7), &
!       ILABEL(8),ILABEL(9),ILABEL(10),ILABEL(11),ILABEL(12)
!       WRITE(6,*)' FCST_PRD AFTER   ',FCST_PRD
!       WRITE(6,*)' -----------------------------------------------'
!       WRITE(6,*)' -----------------------------------------------'
        ELSE
          WRITE(7,*) 'Missing temporal data.  Cannot recalculate '
          WRITE(7,*) 'data time.'
          WRITE(7,*) 'CONVERSION TO SEQUENTIAL PP FILE FAILED IN '
          WRITE(7,*) 'FIELD NO.',ENTRY_NO
          RETURN
        ENDIF


      ENDIF

!     data now in IEEE format
      ILABEL(LBPACK) = ILABEL(LBPACK) - NUMBER_FORMAT*1000 + 3000

!L Convert ILABEL to IEEE  Integer
      BIT_OFF=0
      IEEE_ADDR = 1
      IER=IEEE2IEG(2,LEN_ILABEL,IEEE_LABEL(IEEE_ADDR),                  &
     &            BIT_OFF,ILABEL,1,64,32)
        IF(IER /= 0) THEN
          ICODE=1
          CMESSAGE=' CRAY_IEEE error converting INT for IEEE_LABEL'
          RETURN
        ENDIF
!L Convert RLABEL to IEEE  Real.
      IEEE_ADDR=LEN_ILABEL/2
      IF(IEEE_ADDR*2 /= LEN_ILABEL) BIT_OFF=32
      IEEE_ADDR=IEEE_ADDR+1
      IER=IEEE2IEG(3,LEN_RLABEL,IEEE_LABEL(IEEE_ADDR),                  &
     &            BIT_OFF,RLABEL,1,64,32)
        IF(IER /= 0) THEN
          ICODE=1
          CMESSAGE=' CRAY_IEEE error converting REAL for IEEE_LABEL'
          RETURN
        ENDIF
      BIT_OFF=0
      IF(.NOT.PACKED) THEN
!L Convert Real DATA to IEEE Real if not packed.
        IF(ILABEL(DATA_TYPE) == 1) THEN          !Data Type Real
        IF(ilabel(32) == 74) THEN
          WRITE(6,*)'convert type 74 as logical and reset datatype'
          IER = IEEE2IEG(5,NUM_VALUES-IEXTRAW,IEEE_FIELD,               &
     &            BIT_OFF,FIELD,1,64,32)
          ILABEL(DATA_TYPE) = 3
        ELSE
          IER = IEEE2IEG(3,NUM_VALUES-IEXTRAW,IEEE_FIELD                &
     &            ,BIT_OFF,FIELD,1,64,32)
          IF(IER /= 0) THEN
            ICODE=1
            CMESSAGE='CRAY_IEEE error converting REAL for IEEE_FIELD'
            RETURN
          ENDIF
        ENDIF
!L Convert Integer data to IEEE Integer.
        ELSEIF(ILABEL(DATA_TYPE) == 2) THEN      !Data Type Integer
          IER = IEEE2IEG(2,NUM_VALUES-IEXTRAW,IEEE_FIELD,               &
     &            BIT_OFF,FIELD,1,64,32)
          IF(IER /= 0) THEN
            ICODE=1
            CMESSAGE='CRAY_IEEE error calling USICTI for IEEE_FIELD'
            RETURN
          ENDIF
        ELSEIF(ILABEL(DATA_TYPE) == 3) THEN      !Data Type Logical
          IER = IEEE2IEG(5,NUM_VALUES-IEXTRAW,IEEE_FIELD,               &
     &            BIT_OFF,FIELD,1,64,32)
          IF(IER /= 0) THEN
            ICODE=1
            CMESSAGE='CRAY_IEEE error converting logical forIEEE_FIELD'
            RETURN
          ENDIF
        ENDIF
      ELSE
        WRITE(6,*)'WARNING ! WGDOS packed data - contains IBM reals'

!       code to be added here to convert ibm reals in packed data
!       to ieee reals

      ENDIF

!L process extra data
      IF (IEXTRAW >  0) THEN ! process extra data as got some
!L init values for while loop
        ADDR=NUM_VALUES-IEXTRAW+1 ! start address in field for extra dat
        IEEE_ADDR=(ADDR+1)/2
        IF (IEEE_ADDR*2 == ADDR) THEN
          BIT_OFF=32
        ELSE
          BIT_OFF=0
        ENDIF

        DO WHILE (ADDR <  NUM_VALUES)
!L main while loop that works out code and then checks that code is ok.
!L if code is ok then data_values will contain the number of REAL entries
!L in the vector.
! DEPENDS ON: int_from_real
          IV=INT_FROM_REAL(FIELD(ADDR))
! DEPENDS ON: check_extra
          CALL CHECK_EXTRA(IV,DATA_VALUES,ICODE,CMESSAGE)
          IF (ICODE /= 0) THEN
            RETURN
          ENDIF
          IER = IEEE2IEG(2,1,IEEE_FIELD(IEEE_ADDR),BIT_OFF,             &
     &            FIELD(ADDR),1,64,32)
!         convert the integer from cray format to IEEE format
          IF (IER /= 0) THEN
            ICODE=1
            CMESSAGE='CRAY_IEEE: FAILED IN INTEGER CONV OF EXTRA DATA'
            RETURN
          ENDIF

!L         update bit_off, addr and IEEE_addr
          IF (BIT_OFF == 0) THEN
            BIT_OFF=32
          ELSE
            BIT_OFF=0
            IEEE_ADDR=IEEE_ADDR+1 ! gone on another word..
          ENDIF
          ADDR=ADDR+1             ! increment address
!L now to convert REAL vector to IEEE format.
          IER=IEEE2IEG(3,DATA_VALUES,IEEE_FIELD(IEEE_ADDR),             &
     &      BIT_OFF,FIELD(ADDR),1,64,32)
!         convert the real data values
          IF (IER /= 0) THEN
            ICODE=1
            CMESSAGE='CRAY_IEEE: FAILED IN REAL CONV OF EXTRA DATA'
            RETURN
          ENDIF
!L update loop variables.
          ADDR=ADDR+DATA_VALUES
          IEEE_ADDR=IEEE_ADDR+DATA_VALUES/2
          IF ((DATA_VALUES/2)*2 /= DATA_VALUES) THEN ! odd no. of values
            IF (BIT_OFF == 0) THEN
              BIT_OFF=32
            ELSE
              BIT_OFF=0
              IEEE_ADDR=IEEE_ADDR+1 ! gone on another word..
            ENDIF
          ENDIF
        ENDDO                   ! continue until run out of data....
!L Verify addr and IEEE_addr have correct values at end of whileloop
!L first check that addr is ok
        IF (ADDR /= NUM_VALUES+1) THEN
          WRITE(CMESSAGE,'(A,I5,A,I5)') 'CRAY_IEEE: ADDR',ADDR, &
          ' <> NUM_VALUES+1',NUM_VALUES+1
          ICODE=1
          RETURN
        ENDIF
!L and so is IEEE_addr
        IF (BIT_OFF == 0) IEEE_ADDR=IEEE_ADDR-1
        IF (IEEE_ADDR /= (NUM_VALUES+1)/2) THEN
          WRITE(CMESSAGE,'(A,I5,A,I5)') 'CRAY_IEEE: IEEE_ADDR ', &
          IEEE_ADDR,' <> (NUM_VALUES+1)/2',(NUM_VALUES+1)/2
          ICODE=1
          RETURN
        ENDIF
      ENDIF ! end processing of extra data

      IF(PACKED) THEN
        WRITE(COS_PPUNIT) IEEE_LABEL
        WRITE(COS_PPUNIT) (FIELD(I),I=1,ILABEL(LBLREC)/2)
      ELSE
        WRITE(COS_PPUNIT) IEEE_LABEL
        WRITE(COS_PPUNIT) IEEE_FIELD
      ENDIF
!
! The last field has been processed. An extra field is now written
! to act as a delimiter for the M08 software. This extra fields is
! a duplicate,but with a PP field code of -99 .  
! If MASS is true then last record (-99) is omitted.
      IF(LAST .AND. .NOT. MASS) THEN
        WRITE(6,'(A)') '  WRITING LAST RECORD IN THE COS FILE '
        ILABEL(23)=-99
!L Convert ILABEL to IEEE  Integer
        BIT_OFF=0
        IEEE_ADDR = 1
        IER=IEEE2IEG(2,LEN_ILABEL,IEEE_LABEL(IEEE_ADDR),                &
                     BIT_OFF,ILABEL,1,64,32)
        IF(IER /= 0) THEN
          ICODE=1
          CMESSAGE=' CRAY_IEEE error converting INT for IEEE_LABEL'
          RETURN
        ENDIF
!L Convert RLABEL to IEEE  Real.
        IEEE_ADDR=LEN_ILABEL/2
        IF(IEEE_ADDR*2 /= LEN_ILABEL) BIT_OFF=32
        IEEE_ADDR=IEEE_ADDR+1
        IER=IEEE2IEG(3,LEN_RLABEL,IEEE_LABEL(IEEE_ADDR),                &
                     BIT_OFF,RLABEL,1,64,32)
        IF(IER /= 0) THEN
          ICODE=1
          CMESSAGE=' CRAY_IEEE error converting REAL for IEEE_LABEL'
          RETURN
        ENDIF
        IF (PACKED) THEN
          WRITE(COS_PPUNIT) IEEE_LABEL
          WRITE(COS_PPUNIT) (FIELD(I),I=1,ILABEL(LBLREC)/2)
        ELSE
          WRITE(COS_PPUNIT) IEEE_LABEL
          WRITE(COS_PPUNIT) IEEE_FIELD
        ENDIF
      ENDIF
 9999 CONTINUE
      RETURN
      END SUBROUTINE CRAY_IEEE
