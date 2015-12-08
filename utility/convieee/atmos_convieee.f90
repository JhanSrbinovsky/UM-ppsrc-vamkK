! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL   PROGRAM MAIN_CONVIEEE and SUBROUTINE CONVIEEE ------------------
!LL
!LL
!    Purpose: Converts a dump, ancillary or fieldsfile
!             (atmosphere or ocean) between different
!             into 32-bit or 64-bit IEEE format or vice-versa.
!             The following conversions are supported:-
!               On a IEEE machine
!                 64-bit IEEE to 32-bit IEEE
!                 64-bit IEEE to 64-bit IEEE
!             In either case, WGDOS data will be unpacked if
!             requested.
!
!             MAIN_CONVIEEE reads in fixed length and integer
!             headers of UM file to be converted, extracts dimensions
!             of file and then passes these values to
!             subroutine CONVIEEE.
!
!            CONVIEEE reads in headers and data fields from unit NFTIN
!            converts them to IEEE format and writes them to NFTOUT.
!
!    Documentation: UM Doc Paper F5
!
!LLEND----------------------------------------------------------------
!*L  Arguments:-------------------------------------------------------


!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs

      SUBROUTINE ATMOS_CONVIEEE(NFTIN,NFTOUT,IEEE_TYPE,MAX_FIELD_SIZE   &
     &,                         LEN_FIXHD,LEN1_LOOKUP,LEN2_LOOKUP       &
     &,                         OLD_FIXHD_160,NEW_FIXHD_160             &
     &,                         LOOKUP_LBLREC                           &
     &,                         LOOKUP_LBEGIN,LOOKUP_LBNREC             &
     &,                         LOOKUP_LBEGIN_NEW,LOOKUP_LBNREC_NEW,    &
     &                          FIXHD,LOOKUP,WGDOS_EXPAND)
      USE IO
      USE ereport_mod, ONLY : ereport
      USE PrintStatus_mod
      USE lookup_addresses
! version_mod items required by cstash.h
      USE version_mod, ONLY :                                           &
          nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
          nlevp, npslevp, npslistp, outfile_s, outfile_e
      USE Submodel_Mod

      IMPLICIT NONE

      INTEGER IEEE_TYPE
      INTEGER LEN_FIXHD
      INTEGER LEN1_LOOKUP
      INTEGER LEN2_LOOKUP
      INTEGER MAX_FIELD_SIZE
      INTEGER NEW_FIXHD_160
      INTEGER NFTIN
      INTEGER NFTOUT
      INTEGER OLD_FIXHD_160
      INTEGER WGDOS_EXPAND

      INTEGER FIXHD(LEN_FIXHD)
      INTEGER LOOKUP_LBLREC(LEN2_LOOKUP)
      INTEGER LOOKUP_LBEGIN(LEN2_LOOKUP)
      INTEGER LOOKUP_LBNREC(LEN2_LOOKUP)
      INTEGER LOOKUP_LBEGIN_NEW(LEN2_LOOKUP)
      INTEGER LOOKUP_LBNREC_NEW(LEN2_LOOKUP)
      INTEGER LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)
! Local arrays:--------------------------------------------------------
      INTEGER D1(MAX_FIELD_SIZE) ! Data array used to read in each field
      REAL IEEE_64(MAX_FIELD_SIZE) !Array containing 64 bit IEEE data
!----------------------------------------------------------------------
! Local variables:-----------------------------------------------------
      INTEGER      I,J      ! Loop variables
      INTEGER      K        ! Return code from CRAY intrinsic functions
      INTEGER      ICODE    ! Error return code from READFLDS
      INTEGER      ITYPE    ! Conversion type
      INTEGER      LEN_IO   ! I/O length
      INTEGER      TOT_VALUES   ! Normal data + extra data (if applied)
      INTEGER      IEXTRAW      ! Extra data length
      INTEGER      NET_DATA_LEN ! Normal data length
      INTEGER      ADDR         ! Start address for extra data
      INTEGER      IEEE_ADDR    ! Start address for extra data in IEEE
      INTEGER      BIT_OFF      ! Bit offset in IEEE array
      INTEGER      CODE         ! Encoded info for extra data set
      INTEGER      DATA_VALUES  ! Decoded real extra data length
      INTEGER      PACK_CODE    ! Packing code
      INTEGER      IEEE2IEG
      INTEGER      INT_FROM_REAL ! Function to convert real to int
      REAL         A        ! Return code from BUFFIN; -1.0 = O.K.
      CHARACTER(LEN=80) CMESSAGE ! Character string returned if ICODE  /=  0
!----------------------------------------------------------------------
! Input arguments for decompose_smexe
      INTEGER                                                           &

     &  global_row_len,                                                 &
                          ! IN  :number of E-W points of entire model
     &  global_n_rows    ! IN  :number of P rows of entire mode

! COMDECK CSTASH
! Description:
!  Contains variables and arrays involved in STASH request specification
!  and STASH processing in the UM. Includes namelist STASH.
!
!
! Declarations:
! Imported global variables:
!    None, but see note above.

! Global parameters:

! Global scalars:
      INTEGER      NDIAG   ! No. of diagnostics
      INTEGER      NTPROF  ! No. of time profiles
      INTEGER      NSERIES ! No. of stash time series
      INTEGER      NDPROF  ! No. of domain profiles
      INTEGER      NUPROF  ! No. of useage profiles

! Global dynamic arrays:

!   STASH specification table (JSTASH file):
!   NDIAGPM set in VERSION_MOD module
      INTEGER      MODL_B(NDIAGPM)  ! Internal model no.
      INTEGER      ISEC_B(NDIAGPM)  ! Section
      INTEGER      ITEM_B(NDIAGPM)  ! Item
      INTEGER      ITIM_B(NDIAGPM)  ! Time profile number
      INTEGER      IDOM_B(NDIAGPM)  ! Domain profile number
      INTEGER      IUSE_B(NDIAGPM)  ! Useage profile number

!   Time profile information:

      CHARACTER(LEN=8)  TIMPRO(NPROFTP)         ! Name of profile
      INTEGER      ITYP_T(NPROFTP)         ! Type of profile
      INTEGER      INTV_T(NPROFTP)         ! Time Interval
      CHARACTER(LEN=2)  UNT1_T(NPROFTP)         ! Units for time interval
      INTEGER      ISAM_T(NPROFTP)         ! Sampling period
      CHARACTER(LEN=2)  UNT2_T(NPROFTP)         ! Units for sampling period
      INTEGER      IOPT_T(NPROFTP)         ! Output option
      INTEGER      ISTR_T(NPROFTP)         ! Output Start time
      INTEGER      IEND_T(NPROFTP)         ! Output End time
      INTEGER      ISDT_T(6, NPROFTP)      ! Output Start date
      INTEGER      IEDT_T(6, NPROFTP)      ! Output End date
      INTEGER      IFRE_T(NPROFTP)         ! Output frequency
      INTEGER      IOFF_T(NPROFTP)         ! Offset for sampling
      CHARACTER(LEN=2)  UNT3_T(NPROFTP)         ! Units for output times
      INTEGER      ITIM_T(NPROFTP)         ! No. of times in times table
      INTEGER      ISER_T(NTIMEP ,NPROFTP) ! Times table (with units)
      INTEGER      MODL_T(NPROFTP)         ! Indicates internal model
                                           !  for each times table

!   Domain profile information:

      CHARACTER(LEN=8) DOMPRO  (NPROFDP)           ! Name of domain profile
      INTEGER     IOPL_D  (NPROFDP)           ! Levels option
      INTEGER     LEVB_D  (NPROFDP)           ! Bottom level
      INTEGER     LEVT_D  (NPROFDP)           ! Top level
      INTEGER     IOPA_D  (NPROFDP)           ! Area option
      INTEGER     INTH_D  (NPROFDP)           ! North boundary
      INTEGER     ISTH_D  (NPROFDP)           ! South boundary
      INTEGER     IEST_D  (NPROFDP)           ! East boundary
      INTEGER     IWST_D  (NPROFDP)           ! West boundary
      INTEGER     IMSK_D  (NPROFDP)           ! Mask type
      INTEGER     IMN_D   (NPROFDP)           ! Meaning option
      INTEGER     IWT_D   (NPROFDP)           ! Weighting option
      CHARACTER(LEN=1) TS_D    (NPROFDP)           ! Time series profile
      INTEGER     IG_TS
      INTEGER     I1_TS
      INTEGER     I51_TS
      INTEGER     BLIM_TS (NTimSerP)
      INTEGER     TLIM_TS (NTimSerP)
      REAL        BLIMR_TS(NTimSerP)
      REAL        TLIMR_TS(NTimSerP)
      INTEGER     NLIM_TS (NTimSerP)
      INTEGER     SLIM_TS (NTimSerP)
      INTEGER     ELIM_TS (NTimSerP)
      INTEGER     WLIM_TS (NTimSerP)
      INTEGER     ILEV_D  (NPROFDP)           ! Output levels code
      INTEGER     LEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      REAL       RLEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      INTEGER     PLT_D   (NPROFDP)
      INTEGER     PLLEN_D (NPROFDP)
      INTEGER     PLPOS_D (NPROFDP)
      INTEGER     PSLIST_D(NPSLEVP ,NPSLISTP)
      INTEGER     NPSLISTS
      EQUIVALENCE        (RLEVLST_D,LEVLST_D)

! Useage information:

      CHARACTER(LEN=8) USEPRO(NPROFUP)   ! Name of useage profile
      INTEGER     LOCN_U(NPROFUP)   ! Storage location of profile
      INTEGER     IUNT_U(NPROFUP)   ! Unit no.

! Information from ppxref file:

      INTEGER      MODEL_ST       ! Internal model number
      INTEGER      ISPACE         ! Space code
      INTEGER      ITIMA          ! Time availability code
      INTEGER      IGP            ! Grid of data code
      INTEGER      ILEV           ! Level type code
      INTEGER      IBOT           ! First level code
      INTEGER      ITOP           ! Last level code
      INTEGER      IFLAG          ! Level compression flag
      INTEGER      IOPN(6)        ! Sectional option code
      INTEGER      VMSK           ! Integer equiv of bin vers mask
      INTEGER      IPSEUDO        ! Pseudo dimension type
      INTEGER      IPFIRST        ! First pseudo dim code
      INTEGER      IPLAST         ! Last pseudo dim code
      INTEGER      PTR_PROG       ! Section zero point back
      INTEGER      HALO_TYPE      ! Type of halo the field has

! PP output file units
      INTEGER      PPlen2LkUp(OUTFILE_S:OUTFILE_E)
      CHARACTER(LEN=1)  FTOutUnit (OUTFILE_S:OUTFILE_E)

! COMMON blocks:
      COMMON/STCHA/ TIMPRO,UNT1_T,UNT2_T,UNT3_T,DOMPRO,TS_D,            &
     &  USEPRO,FTOutUnit

      COMMON/STSH/                                                      &
     &  NDIAG   ,MODL_B  ,ISEC_B ,ITEM_B  ,ITIM_B  ,IDOM_B  ,IUSE_B,    &
     &  NTPROF  ,ITYP_T  ,INTV_T ,ISAM_T  ,ITIM_T  ,                    &
     &  IOPT_T  ,ISTR_T  ,IEND_T ,IFRE_T  ,IOFF_T, ISER_T  ,MODL_T  ,   &
     &  NDPROF  ,IOPL_D  ,LEVB_D ,ISDT_T  ,IEDT_T  ,                    &
     &  IOPA_D  ,INTH_D  ,ISTH_D ,IEST_D  ,IWST_D  ,                    &
     &  IMSK_D  ,IMN_D   ,IWT_D  ,                                      &
     &  LEVT_D  ,LEVLST_D,                                              &
     &  PLT_D   ,PLLEN_D ,PLPOS_D,PSLIST_D,NPSLISTS,                    &
     &  BLIM_TS ,TLIM_TS ,BLIMR_TS,TLIMR_TS,IG_TS   ,I1_TS   ,          &
     &  NLIM_TS ,SLIM_TS ,ELIM_TS ,WLIM_TS ,I51_TS  ,NSERIES ,          &
     &  NUPROF  ,LOCN_U  ,IUNT_U ,                                      &
     &  MODEL_ST,ISPACE  ,ITIMA  ,IGP     ,                             &
     &  ILEV    ,IBOT    ,ITOP   ,IFLAG   ,IOPN    ,VMSK    ,           &
     &  IPSEUDO ,IPFIRST ,IPLAST ,PTR_PROG, HALO_TYPE,                  &
     & PPlen2LkUp

! CSTASH end
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

! Loop over all fields
      DO I=1,LEN2_LOOKUP

! Check for the end of a PP format lookup table
        IF(LOOKUP(1,I) == -99) GOTO 2000

! Reset the headers in case WDGOS packing has altered them
        IF(WGDOS_EXPAND == 1) THEN
        LOOKUP(LBLREC,I)=LOOKUP_LBLREC(I)
        LOOKUP(LBEGIN,I)=LOOKUP_LBEGIN(I)
        LOOKUP(LBNREC,I)=LOOKUP_LBNREC(I)
        FIXHD(160)=OLD_FIXHD_160
        ENDIF

! Check if this field has already been converted - WGDOS only
        IF(MOD(LOOKUP(21,I),10) == 1 .AND. WGDOS_EXPAND /= 1) THEN
! Read in field
          IF(IEEE_TYPE == 32)THEN
! DEPENDS ON: readflds
            CALL READFLDS(NFTIN,1,I,LOOKUP,LEN1_LOOKUP,                 &
     &                    IEEE_64,MAX_FIELD_SIZE,FIXHD,                 &
     &                    WGDOS_EXPAND,ICODE,CMESSAGE)
          ELSEIF(IEEE_TYPE == 64)THEN
! Read in field
! DEPENDS ON: readflds
            CALL READFLDS(NFTIN,1,I,LOOKUP,LEN1_LOOKUP,                 &
     &                    IEEE_64,MAX_FIELD_SIZE,FIXHD,                 &
     &                    WGDOS_EXPAND,ICODE,CMESSAGE)
          END IF
! DEPENDS ON: abort_io
          IF(ICODE /= 0)CALL ABORT_IO('CONVIEEE',CMESSAGE,ICODE,NFTIN)

        ELSE IF(MOD(LOOKUP(21,I),10) == 1 .AND. WGDOS_EXPAND == 1 .AND. &
     &          IEEE_TYPE == 64 ) THEN
! Read in field
! DEPENDS ON: readflds
          CALL READFLDS(NFTIN,1,I,LOOKUP,LEN1_LOOKUP,                   &
     &                  IEEE_64,MAX_FIELD_SIZE,FIXHD,                   &
     &                  WGDOS_EXPAND,ICODE,CMESSAGE)
! DEPENDS ON: abort_io
          IF(ICODE /= 0)CALL ABORT_IO('CONVIEEE',CMESSAGE,ICODE,NFTIN)

        ELSE
! Read in field
! DEPENDS ON: readflds
          CALL READFLDS(NFTIN,1,I,LOOKUP,LEN1_LOOKUP,                   &
     &                  D1,MAX_FIELD_SIZE,FIXHD,                        &
     &                  WGDOS_EXPAND,ICODE,CMESSAGE)
! DEPENDS ON: abort_io
          IF(ICODE /= 0)CALL ABORT_IO('CONVIEEE',CMESSAGE,ICODE,NFTIN)

! Set data type
          IF(ABS(LOOKUP(DATA_TYPE,I)) == 1) THEN
! Type real
            IF(IEEE_TYPE == 32)THEN
              ITYPE=3
            ELSEIF(IEEE_TYPE == 64)THEN
              ITYPE=2
            ENDIF
          ELSE IF(ABS(LOOKUP(DATA_TYPE,I)) == 2) THEN
! Type integer
            IF(IEEE_TYPE == 32)THEN
              ITYPE=2
            ELSEIF(IEEE_TYPE == 64)THEN
              ITYPE=1
            ENDIF
          ELSE IF(ABS(LOOKUP(DATA_TYPE,I)) == 3) THEN
! Type logical
            ITYPE=5
          ELSE
! DEPENDS ON: pr_look
            CALL PR_LOOK(                                               &
     &                   LOOKUP,LOOKUP,LEN1_LOOKUP,I)
            ICODE=3
            CMESSAGE='CONVIEEE: Invalid code in LOOKUP(39,K)'
            RETURN
          ENDIF

! Convert to IEEE format and write to disk
          IF(ITYPE >= 0)THEN
            IF(IEEE_TYPE == 32)THEN

              TOT_VALUES=LOOKUP(LBLREC,I)
              IEXTRAW=0
              IF(LOOKUP(LBEXT,I) >  0) THEN ! got some extra data
                IEXTRAW=LOOKUP(LBEXT,I)

                ! Check that no pack_code 1 file left packed if we have
                ! extra data. (files with pack_code 2 & 3 are forced
                ! unpacked).  This version supports pack_code 4 data

                PACK_CODE=MOD((LOOKUP(LBPACK,I)),10)
                IF((LOOKUP(LBROW,I)*LOOKUP(LBNPT,I)+LOOKUP(LBEXT,I)  /= &
     &              LOOKUP(LBLREC,I)).AND.                              &
     &              (PACK_CODE /= 4)) THEN
                  CMESSAGE=                                             &
     &            'CONVIEE1 : Packing of extra data not supported'
                  WRITE(6,*)'Please use expand option'
                  ICODE=1
! DEPENDS ON: abort_io
                  CALL ABORT_IO('CONVIEEE',CMESSAGE,ICODE,NFTIN)
                ENDIF
              ENDIF

              ! Converting data into IEEE format consists of 3 stages
              ! a) convert normal grid data into IEEE format
              ! b) convert integer header of extra data vector into IEEE
              ! c) convert rest of extra data vector into IEEE format

              ! Process normal data

              NET_DATA_LEN=LOOKUP(LBLREC,I)-IEXTRAW
              IF (NET_DATA_LEN  ==  0) THEN
! In some cases (e.g., CovStats files), fields contain no data.  Instead
! of giving an error from the data conversion  routines, output a
! message and skip conversion.
                WRITE(6,*) 'Data length = 0 for field ',i,              &
     &                     ' - nothing to convert'
              ELSE
              K=IEEE2IEG(ITYPE,NET_DATA_LEN,IEEE_64,0,                  &
     &                 D1,1,64,IEEE_TYPE)
              IF(K /= 0)THEN
                WRITE(6,*)'Conversion Error - Return Code is ',K
                DO J=1,NET_DATA_LEN
                    WRITE(6,'(''Error converting field '',i5,           &
     &                        '' : Stash Code '',i5,                    &
     &                        '' : Point No. '',i5)')                   &
     &                        I, LOOKUP(ITEM_CODE,I),J
                    WRITE(6,*) 'Number unconvertable reset to RMDI'
                    IEEE_64(J)=RMDI
                END DO
              END IF
              END IF
              ! Process extra data

              ! About BIT OFFSET
              !
              ! 1         2         3         4         5 (addr)
              ! |---------|---------|---------|---------|  FIELD
              ! .        .         .
              ! .      .       .
              ! .    .    .
              ! |----|----|----|----|----|----|----|----|  IEEE_FIELD
              ! 1         2         3         4         5 (ieee_addr)
              !                     |    |
              ! <--------->         |    |
              !  a "word"           | bit_off=32
              !                 bit_off=0
              ! Example:
              !  if ADDR=2, IEEE_ADDR=3/2=1
              !  IEEE_ADDR*2 eq 1;  so BIT_OFF=32
              !


              IF (IEXTRAW >  0) THEN ! process extra data as got some
                ! init values for while loop
                ! start address in field for extra data
                ADDR=TOT_VALUES-IEXTRAW+1
                IEEE_ADDR=(ADDR+1)/2
                IF (IEEE_ADDR*2 == ADDR) THEN
                  BIT_OFF=32
                ELSE
                  BIT_OFF=0
                ENDIF

                DO WHILE (ADDR <  TOT_VALUES)
                ! CODE is integer header for extra data vector
                ! which contains encoded info for extra data -
                ! vector length & data type
                ! Decode CODE: data_values will be vector length
                ! Details about extra data, see Paper F3
                ! NB. integer header for extra data vector is converted
                !     to its real EQUIVALENCE during model run.  Hence,
                !     INT_FROM_REAL serves to convert it back to INTEGER

! DEPENDS ON: int_from_real
                  CODE=INT_FROM_REAL(D1(ADDR))
! DEPENDS ON: check_extra
                  CALL CHECK_EXTRA(CODE,DATA_VALUES,ICODE,CMESSAGE)
                  IF (ICODE /= 0) THEN
                    write(6,*)'Fail in CHECK_EXTRA'
                    RETURN
                  ENDIF

                  ! Convert integer extra_data head into IEEE
                  K=IEEE2IEG(2,1,IEEE_64(IEEE_ADDR),BIT_OFF,            &
     &                      D1(ADDR),1,64,IEEE_TYPE)

                  IF (K /= 0) THEN
                    ICODE=1
                    CMESSAGE=                                           &
     &              'CONVIEE1: failed in integer conv of extra data'
                    RETURN
                  ENDIF

                  ! update bit_off, addr and ibm_addr
                  IF (BIT_OFF == 0) THEN
                    BIT_OFF=32
                  ELSE
                    BIT_OFF=0
                    IEEE_ADDR=IEEE_ADDR+1 ! GONE ON ANOTHER WORD..
                  ENDIF
                  ADDR=ADDR+1           ! INCREMENT ADDRESS

                  ! Convert REAL vector to IEEE format.
                  K=IEEE2IEG(3,DATA_VALUES,IEEE_64(IEEE_ADDR),          &
     &                      BIT_OFF,D1(ADDR),1,64,IEEE_TYPE)

                  IF (K /= 0) THEN
                    ICODE=1
                    CMESSAGE=                                           &
     &              'CONVIEE1: FAILED IN REAL CONV OF EXTRA DATA'
                    RETURN
                  ENDIF

                  ! Update loop variables.
                  ADDR=ADDR+DATA_VALUES
                  IEEE_ADDR=IEEE_ADDR+DATA_VALUES/2
                  ! Odd no. of values.
                  IF ((DATA_VALUES/2)*2 /= DATA_VALUES) THEN
                    IF (BIT_OFF == 0) THEN
                      BIT_OFF=32
                    ELSE
                      BIT_OFF=0
                      IEEE_ADDR=IEEE_ADDR+1 ! GONE ON ANOTHER WORD..
                    ENDIF
                  ENDIF
                ENDDO                 ! continue until run out of data

                ! Verify addr and ieee_addr have correct values at end
                ! of whileloop. First check that addr is ok
                IF (ADDR /= TOT_VALUES+1) THEN
                  WRITE(CMESSAGE,109)ADDR,TOT_VALUES+1
 109              FORMAT('CONVIEE1: addr',i5,1x,'<> tot_values+1',i5)
                  ICODE=1
                  RETURN
                ENDIF
                ! and so is ieee_addr
                IF (BIT_OFF == 0) IEEE_ADDR=IEEE_ADDR-1
                IF (IEEE_ADDR /= (TOT_VALUES+1)/2) THEN
                  WRITE(CMESSAGE,110)IEEE_ADDR,(TOT_VALUES+1)/2
 110              FORMAT('CONVIEE1: ieee_addr ',i5,1x,                  &
     &                   ' <> (tot_values+1)/2',i5)
                  ICODE=1
                  RETURN
                ENDIF
              ENDIF ! end processing of extra data
            END IF

! 64 bit case - just copy into ieee_64 with appropriate type change
            IF (IEEE_TYPE == 64) THEN
              DO K = 1, LOOKUP(LBLREC,I)
                IEEE_64(K) = TRANSFER(D1(K),IEEE_64(1))
              END DO
            END IF

          ELSE
            DO K=1,LOOKUP(LBLREC,I)
              IEEE_64(k)=IAND(D1(K),1)
            END DO
          ENDIF
        ENDIF

! Write out field
        FIXHD(160)=NEW_FIXHD_160
        IF(IEEE_TYPE == 32)THEN
          CALL SETPOS32(NFTOUT, LOOKUP_LBEGIN_NEW(I), K)
          CALL BUFFOUT(NFTOUT, IEEE_64, (LOOKUP_LBNREC_NEW(I)+1)/2, LEN_IO, A)
          LEN_IO=LEN_IO*2
        ELSEIF(IEEE_TYPE == 64)THEN
          CALL SETPOS(NFTOUT, LOOKUP_LBEGIN_NEW(I), K)
          CALL BUFFOUT(NFTOUT, IEEE_64, LOOKUP_LBNREC_NEW(I), LEN_IO, A)
        ENDIF

! Check for I/O errors
        if(A /= -1.0.OR.LEN_IO /= LOOKUP_LBNREC_NEW(I)) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of data field',                      &
     &                 A,LEN_IO,LOOKUP(15,I))

          CALL EREPORT('ATMOS_CONVIEEE', 1002,                          &
     &     'Buffer out of data field wrong size')

        ENDIF

        If ( PrintStatus >= PrStatus_Diag ) Then
          WRITE(6,'(''Field '',i5,'' : Stash Code '',i5,                &
     &     '' : Level '',i4,                                            &
     &     '' : T + '',i3,                                              &
     &     '' : has been converted'')') I, LOOKUP(42,I)                 &
     &          , LOOKUP(33,I), LOOKUP(14,I)
        End If

! Reset the headers in case WDGOS packing has altered them
        IF(WGDOS_EXPAND == 1) THEN
        LOOKUP(LBLREC,I)=LOOKUP_LBLREC(I)
        LOOKUP(LBEGIN,I)=LOOKUP_LBEGIN(I)
        LOOKUP(LBNREC,I)=LOOKUP_LBNREC(I)
        FIXHD(160)=OLD_FIXHD_160
        ENDIF

      END DO
2000  CONTINUE

      RETURN
      END SUBROUTINE ATMOS_CONVIEEE
