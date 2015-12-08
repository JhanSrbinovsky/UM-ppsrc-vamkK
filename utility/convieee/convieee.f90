! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    SUBROUTINE CONVIEEE ------------------
!
!
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
!            CONVIEEE reads in headers and data fields from unit NFTIN
!            converts them to IEEE format and writes them to NFTOUT.
!
!    Documentation: UM Doc Paper F5
!
! Arguments:-------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      SUBROUTINE CONVIEEE                                               &
     &  (LEN_FIXHD,LEN_INTHD,LEN_REALHD,                                &
     &  LEN1_LEVDEPC,LEN2_LEVDEPC,LEN1_ROWDEPC,                         &
     &  LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                         &
     &  LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,                          &
     &  LEN_DUMPHIST,LEN_CFI1,LEN_CFI2,LEN_CFI3,                        &
     &  LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,P_FIELD,                       &
     &  P_ROWS,ROW_LENGTH,                                              &
     &  NFTIN,NFTOUT,IEEE_TYPE,                                         &
     &  MAX_FIELD_SIZE, WGDOS_EXPAND)

      USE io_configuration_mod, ONLY : &
          io_data_alignment,       &
          io_field_padding
      USE ereport_mod, ONLY : ereport
      USE Decomp_DB
      USE UM_ParVars
      USE writhead_mod
      USE lookup_addresses
      USE ppxlook_mod, ONLY : ppxrecs
      USE cppxref_mod, ONLY : ppx_grid_type
! version_mod items required by cstash.h
      USE version_mod, ONLY :                                           &
          nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
          nlevp, npslevp, npslistp, outfile_s, outfile_e
      USE Submodel_Mod

      IMPLICIT NONE

      INTEGER                                                           &

     & LEN_FIXHD                                                        &
                    !IN Length of fixed length header on input file
     &,LEN_INTHD                                                        &
                    !IN Length of integer header on input file
     &,LEN_REALHD                                                       &
                    !IN Length of real header on input file
     &,LEN1_LEVDEPC                                                     &
                    !IN 1st dim of lev dependent consts on input file
     &,LEN2_LEVDEPC                                                     &
                    !IN 2nd dim of lev dependent consts on input file
     &,LEN1_ROWDEPC                                                     &
                    !IN 1st dim of row dependent consts on input file
     &,LEN2_ROWDEPC                                                     &
                    !IN 2nd dim of row dependent consts on input file
     &,LEN1_COLDEPC                                                     &
                    !IN 1st dim of col dependent consts on input file
     &,LEN2_COLDEPC                                                     &
                    !IN 2nd dim of col dependent consts on input file
     &,LEN1_FLDDEPC                                                     &
                    !IN 1st dim of field dependent consts on input fi
     &,LEN2_FLDDEPC                                                     &
                    !IN 2nd dim of field dependent consts on input fi
     &,LEN_EXTCNST                                                      &
                    !IN Length of extra consts on input file
     &,LEN_DUMPHIST                                                     &
                    !IN Length of history header on input file
     &,LEN_CFI1                                                         &
                    !IN Length of index1 on input file
     &,LEN_CFI2                                                         &
                    !IN Length of index2 on input file
     &,LEN_CFI3                                                         &
                    !IN Length of index3 on input file
     &,LEN1_LOOKUP                                                      &
                    !IN 1st dim of LOOKUP on input file
     &,LEN2_LOOKUP                                                      &
                    !IN 2nd dim of LOOKUP on input file
     &,LEN_DATA                                                         &
                    !IN Length of data on input file
     &,P_FIELD                                                          &
                    !IN No of p-points per level on input file
     &,P_ROWS                                                           &
     &,ROW_LENGTH                                                       &
     &,MAX_FIELD_SIZE !Maximum field size on file
      integer wgdos_expand  ! set to 1 for expansion of WGDOS Fields

      INTEGER                                                           &
     & NFTIN                                                            &
                !IN Unit number of input UM dump
     &,NFTOUT                                                           &
                !IN Unit number of output IEEE dump
     &,IEEE_TYPE ! Output file precision

! Local arrays:---------------------------------------------------------
      INTEGER                                                           &
     & FIXHD(LEN_FIXHD),                                                &
                                                 !
     & INTHD(LEN_INTHD),                                                &
                                                 !\  integer
     & CFI1(LEN_CFI1+1),CFI2(LEN_CFI2+1),                               &
                                                 ! > file headers
     & CFI3(LEN_CFI3+1),                                                &
                                                 !/
     & LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP),                                 &
                                                 !
     & LOOKUP_21(LEN2_LOOKUP)                                           &
                                ! Holds values of input LOOKUP(21,K)
     &,LOOKUP_LBNREC(LEN2_LOOKUP)                                       &
     &,lookup_lblrec(len2_lookup)                                       &
     &,lookup_lbegin(len2_lookup)                                       &
                                       ! Old value of lbegin in lookup
     &,lookup_lblrec_new(len2_lookup)                                   &
                                       ! New value of lblrec in lookup
     &,lookup_lbnrec_new(len2_lookup)                                   &
                                       ! New value of lbnrec in lookup
     &,lookup_lbegin_new(len2_lookup)                                   &
                                       ! New value of lbegin in lookup
     &,disk_address                                                     &
                                       ! Current rounded disk address
     &,number_of_data_words_on_disk                                     &
                                       ! Number of data words on disk
     &,number_of_data_words_in_memory                                   &
                                       ! Number of Data Words in memory
     &,old_fixhd_160                                                    &
                                       ! Input value of FIXHD(160)
     &,new_fixhd_160                   ! Output value of FIXHD(160)

      REAL                                                              &
     & REALHD(LEN_REALHD),                                              &
     & LEVDEPC(1+LEN1_LEVDEPC*LEN2_LEVDEPC),                            &
                                                 !
     & ROWDEPC(1+LEN1_ROWDEPC*LEN2_ROWDEPC),                            &
                                                 !
     & COLDEPC(1+LEN1_COLDEPC*LEN2_COLDEPC),                            &
                                                 !\  real
     & FLDDEPC(1+LEN1_FLDDEPC*LEN2_FLDDEPC),                            &
                                                 ! > file headers
     & EXTCNST(LEN_EXTCNST+1),                                          &
                                                 !/
     & DUMPHIST(LEN_DUMPHIST+1)
      INTEGER                                                           &
     & D1(MAX_FIELD_SIZE)  ! Data array used to read in each field

!*----------------------------------------------------------------------
!*L  Local variables:---------------------------------------------------

      INTEGER                                                           &
     & ICODE                                                            &
                    ! Error return code from subroutines
     &,START_BLOCK                                                      &
                    ! READHEAD argument (not used)
     &,I,J,K,L                                                          &
                    ! Loop indices
     &,LEN_IO                                                           &
                    ! I/O length
     &,ITYPE                                                            &
                    ! Conversion type
     &,MODEL                                                            &
                        ! Internal model number
     &,SECTION                                                          &
                        ! Section number
     &,ITEM                                                             &
                        ! Item code
     &,JOC_NO_SEAPTS                                                    &
                        ! Number of points in compressed ocean field
     &,LEN_OCFLD                                                        &
                        ! Length of uncompressed ocean field
     &,INIT_FIXHD_161                                                   &
                        ! Initialised value of FIXHD(161)
     &,PPXREF_GRID_TYPE                                                 &
                        ! Grid type form ppxref
     &,LEN_BUF                                                          &
                        ! Record length of boundary dataset
     &,MAX_LEN_BUF                                                      &
                        ! Maximum record length of boundary dataset
     &,POS              ! Position of field in file

! Input arguments for decompose_smexe
      INTEGER                                                           &

     &  global_row_len,                                                 &
                          ! IN  :number of E-W points of entire model
     &  global_n_rows     ! IN  :number of P rows of entire mode

      INTEGER EXPPXI
      EXTERNAL EXPPXI
      INTEGER RowNumber


      REAL A        !Return code from BUFFIN; -1.0 = O.K.

      CHARACTER                                                         &
     & CMESSAGE*100 ! Character string returned if ICODE  /=  0
      INTEGER NFT1,NFT2
      PARAMETER (NFT1=22, NFT2=2)
!*----------------------------------------------------------------------
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

!L 0. Read in PPXREF
      cmessage = ' '
      ppxRecs=1
      RowNumber=0
      ICODE=0

! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_A',ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_A'
        WRITE(6,*) CMESSAGE

        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   'Error reading STASHmaster_A')
      END IF

! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_O',ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_O'
        WRITE(6,*) CMESSAGE

        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   'Error reading STASHmaster_O')
      END IF

! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_A',RowNumber,                  &
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE

        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   CMESSAGE)
      END IF

! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_O',RowNumber,                  &
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE

        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   CMESSAGE)
      END IF

!User STASHmaster
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(0,' ',ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE

        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   CMESSAGE)
      END IF

! DEPENDS ON: getppx
      CALL GETPPX(0,NFT2,' ',RowNumber,                                 &
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE

        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   CMESSAGE)
      END IF

!L 1. Read in file header

! DEPENDS ON: readhead
      CALL READHEAD(NFTIN,FIXHD,LEN_FIXHD,                              &
     &              INTHD,LEN_INTHD,REALHD,LEN_REALHD,                  &
     &              LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,                  &
     &              ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,                  &
     &              COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                  &
     &              FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC,                  &
     &              EXTCNST,LEN_EXTCNST,DUMPHIST,LEN_DUMPHIST,          &
     &              CFI1,LEN_CFI1,CFI2,LEN_CFI2,CFI3,LEN_CFI3,          &
     &              LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,            &
     &              START_BLOCK,ICODE,CMESSAGE)

      IF(ICODE /= 0)THEN
        WRITE(6,*)CMESSAGE,ICODE

        CALL EREPORT('CONVIEEE', ICODE,                                 &
     &   CMESSAGE)

      ENDIF

! 2: Check for PP format dataset if field to be expanded

! Preserve the original length values for re-use
      DO K=1,LEN2_LOOKUP
        LOOKUP_LBLREC(K)=LOOKUP(LBLREC,K)
        LOOKUP_LBEGIN(K)=LOOKUP(LBEGIN,K)
        LOOKUP_LBNREC(K)=LOOKUP(LBNREC,K)
      END DO

! Get decompostion information
! ----------------------------

        CALL decompose(ROW_LENGTH, P_ROWS,0,0,-99)
        CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
      IF(LOOKUP(LBNREC,1) >  0.AND.FIXHD(12) >  0)THEN
! Check for WGDOS expansion
        IF (WGDOS_EXPAND == 1)THEN
! Issue a message on why we are doing this
          write(6,'(//''***** Initial Scan for PP Format Dataset'',     &
     &     '' *****''/)')
          DO I=1,LEN2_LOOKUP
            IF(LOOKUP(1,I) == -99) GOTO 195
! DEPENDS ON: readflds
            CALL READFLDS(NFTIN,1,I,LOOKUP,LEN1_LOOKUP,                 &
     &                    D1,MAX_FIELD_SIZE,FIXHD,                      &
     &                    WGDOS_EXPAND,ICODE,CMESSAGE)
! DEPENDS ON: abort_io
            IF(ICODE /= 0)CALL ABORT_IO('CONVIEEE',CMESSAGE,ICODE,NFTIN)
          END DO
195       CONTINUE
        END IF
      END IF

! 3: Reset LOOKUP and FIXHD

      INIT_FIXHD_161=0
      OLD_FIXHD_160=FIXHD(160)
      DO K=1,LEN2_LOOKUP

        IF(LOOKUP(1,K) == -99)GOTO 200
! Set LOOKUP(LBNREC) = 0  in old dumps where UM version number
!   not in fixed length header
        IF(FIXHD(12) <  0.AND.FIXHD(5) /= 3)THEN
          LOOKUP(LBNREC,K)=0
        END IF

! Packing code = -2 now obselete, reset packing code to 2
        IF(LOOKUP(LBPACK,K) == -2)LOOKUP(LBPACK,K)=2

! Store values of packing indicator and set least significant
!   number in LOOKUP(LBPACK,K) to 0 to indicate no packing
        LOOKUP_21(K)=LOOKUP(LBPACK,K)
        LOOKUP(LBPACK,K)=MOD(LOOKUP(LBPACK,K),1000)
          IF((MOD(LOOKUP(LBPACK,K),10) /= 1) .AND.                      &
     &                     (MOD(LOOKUP(LBPACK,K),10) /= 4))THEN
            LOOKUP(LBPACK,K)=(LOOKUP(LBPACK,K)/10)*10
          ELSE IF(WGDOS_EXPAND == 1) THEN
            LOOKUP(LBPACK,K)=(LOOKUP(LBPACK,K)/10)*10
          END IF
! For option IEEE_TYPE=32 but the field is WGDOS packed but NOT
! being expanded then need to double the size of LOOKUP(LBLREC,K)
! before the call to set_dumpfile_address so that the output field
! size is calculated correctly
          IF (IEEE_TYPE  ==  32 .AND. WGDOS_EXPAND  ==  0               &
     &                  .AND. (MOD(LOOKUP(LBPACK,K),10)  ==  1)) THEN
            LOOKUP(LBLREC,K) = 2*LOOKUP(LBLREC,K)
          ENDIF
! Process compressed fields
          IF(MOD(LOOKUP(LBPACK, K),1000) == 110)THEN
            IF(K <= (INTHD(14)+2)*INTHD(8))THEN
! Calculate expanded field lengths for ocean compressed fields
              MODEL=LOOKUP(MODEL_CODE, K)
              ITEM=MOD(LOOKUP(ITEM_CODE, K),1000)
              SECTION=(LOOKUP(ITEM_CODE, K)-ITEM)/1000
! DEPENDS ON: exppxi
              PPXREF_GRID_TYPE=EXPPXI(MODEL,SECTION,ITEM,PPX_GRID_TYPE, &
     &                                ICODE,CMESSAGE)
              IF(PPXREF_GRID_TYPE == 36)THEN
! Ocean mass points.
                LOOKUP(LBNPT,K)   = INTHD(6)
                LOOKUP(LBROW,K)   = INTHD(7)
                LOOKUP(LBLREC, K) = INTHD(6)*INTHD(7)
              ELSEIF(PPXREF_GRID_TYPE == 37)THEN
! Ocean velocity points. One less row.
                LOOKUP(LBNPT,K)   = INTHD(6)
                LOOKUP(LBROW,K)   = INTHD(7)-1
                LOOKUP(LBLREC, K) = INTHD(6)*(INTHD(7)-1)
              END IF
              LOOKUP(LBPACK, K) = 0

            ELSE
! Field not compressed onto sea points. Correct packing code
              LOOKUP(LBPACK, K)=MOD(LOOKUP(LBPACK, K),10)
            END IF

          END IF
! Add to length of data
          INIT_FIXHD_161=INIT_FIXHD_161+LOOKUP(LBLREC,K)

        END DO

200     CONTINUE
        FIXHD(160)=FIXHD(150)+FIXHD(151)*FIXHD(152)
        FIXHD(161)=INIT_FIXHD_161
        LEN_DATA=INIT_FIXHD_161

        DO K=1,LEN2_LOOKUP
!  indicate output format.
          LOOKUP(LBPACK,K)=LOOKUP(LBPACK,K)+3000
        END DO

      IF(FIXHD(12) <  208)FIXHD(12)=208

! Boundary datasets are structured differently.
! Skip call to set_dumpfile_address for boundary datasets and
! Calculate addressing for well formed  boundary datasets explicitly.
      IF (FIXHD(5) /= 5)THEN

! Not a boundary dataset. Call set_dumpfile_address
!
!--reset the 32/64 bit lookup headers after packing, etc
!  has been removed
! DEPENDS ON: set_dumpfile_address
      call set_dumpfile_address(fixhd, len_fixhd,                       &
     &                          lookup, len1_lookup,                    &
     &                          len2_lookup,                            &
     &                          number_of_data_words_in_memory,         &
     &                          number_of_data_words_on_disk,           &
     &                          disk_address)
      ELSE

! Boundary  dataset. Calcuate start address from header and round it up
! to ensure we start on a sector boundary
        DISK_ADDRESS=FIXHD(160)-1
        DISK_ADDRESS=((DISK_ADDRESS+IO_DATA_ALIGNMENT-1)/                  &
     &                IO_DATA_ALIGNMENT)*IO_DATA_ALIGNMENT

! Loop over number of times for which data is present in dataset
        DO K=1,INTHD(3)

! Loop over number of different field types present
          LEN_BUF=0
          MAX_LEN_BUF=0
          DO I=1,INTHD(15)
            POS=(K-1)*INTHD(15)+I
            LOOKUP(LBEGIN,POS)=DISK_ADDRESS+LEN_BUF
            LOOKUP(LBNREC,POS)=LOOKUP(LBLREC,POS)
            LEN_BUF=LEN_BUF+LOOKUP(LBLREC,POS)
          END DO
          MAX_LEN_BUF=MAX0(LEN_BUF,MAX_LEN_BUF)

! Update disk address and ensure that next time starts
! on a sector boundary
          DISK_ADDRESS=DISK_ADDRESS+LEN_BUF
          DISK_ADDRESS=((DISK_ADDRESS+IO_FIELD_PADDING-1)/                &
     &                  IO_FIELD_PADDING) * IO_FIELD_PADDING

        END DO

      END IF
!--preserve the new length values for re-use
      new_fixhd_160=fixhd(160)
      do k=1,len2_lookup
        lookup_lblrec_new(k)=lookup(lblrec, k)
        lookup_lbnrec_new(k)=lookup(lbnrec, k)
        lookup_lbegin_new(k)=lookup(lbegin, k)
      end do
!L 1. Write out file header

      CALL WRITHEAD(NFTOUT,FIXHD,LEN_FIXHD,                             &
     &              INTHD,LEN_INTHD,REALHD,LEN_REALHD,                  &
     &              LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,                  &
     &              ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,                  &
     &              COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                  &
     &              FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC,                  &
     &              EXTCNST,LEN_EXTCNST,DUMPHIST,LEN_DUMPHIST,          &
     &              CFI1,LEN_CFI1,CFI2,LEN_CFI2,CFI3,LEN_CFI3,          &
     &              LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,            &
     &              IEEE_TYPE,                                          &
     &              START_BLOCK,ICODE,CMESSAGE)

      IF(ICODE /= 0)THEN
        WRITE(6,*)CMESSAGE,ICODE

        CALL EREPORT('CONVIEEE', ICODE,CMESSAGE)

      ENDIF

! Reset PP file indicator
      IF(LOOKUP_LBNREC(1) >  0)THEN
        DO K=1,LEN2_LOOKUP
          LOOKUP(LBNREC,K)=LOOKUP_LBNREC(K)
          lookup(lblrec,k)=lookup_lblrec(k)
          lookup(lbegin,k)=lookup_lbegin(k)
        ENDDO
      ENDIF
! Restore value of packing indicator

      DO K=1,LEN2_LOOKUP
        LOOKUP(21,K)=LOOKUP_21(K)
      ENDDO

!L 3. Read in each field, convert to IEEE format and write out
!L    results to new file

      IF (FIXHD(2) == 1)THEN

! Atmosphere file
! DEPENDS ON: atmos_convieee
        CALL ATMOS_CONVIEEE(NFTIN,NFTOUT,IEEE_TYPE,MAX_FIELD_SIZE       &
     &,                         LEN_FIXHD,LEN1_LOOKUP,LEN2_LOOKUP       &
     &,                         OLD_FIXHD_160,NEW_FIXHD_160             &
     &,                         LOOKUP_LBLREC                           &
     &,                         LOOKUP_LBEGIN,LOOKUP_LBNREC             &
     &,                         LOOKUP_LBEGIN_NEW,LOOKUP_LBNREC_NEW,    &
     &                          FIXHD,LOOKUP,WGDOS_EXPAND)

      ELSEIF (FIXHD(2) == 2)THEN

! Ocean file

! Calculate sizes of compressed and uncompressed ocean fields:
! First decide whether there are any compressed fields (use LBPACK
! to determine whether the first field contains sea points only)
        IF( MOD(LOOKUP(LBPACK, 1)/10,10)  ==  0) THEN
          JOC_NO_SEAPTS=1
          LEN_OCFLD    =1
        ELSE
        JOC_NO_SEAPTS=INTHD(11)
        LEN_OCFLD    =INTHD(6)*INTHD(7)*INTHD(8)
        END IF

! DEPENDS ON: ocean_convieee
        CALL OCEAN_CONVIEEE(NFTIN,NFTOUT,IEEE_TYPE,MAX_FIELD_SIZE       &
     &,                         LEN_FIXHD,LEN_INTHD                     &
     &,                         LEN_CFI1,LEN_CFI2,LEN_CFI3              &
     &,                         LEN1_LOOKUP,LEN2_LOOKUP                 &
     &,                         JOC_NO_SEAPTS,LEN_OCFLD                 &
     &,                         OLD_FIXHD_160,NEW_FIXHD_160             &
     &,                         LOOKUP_LBLREC,LOOKUP_LBLREC_NEW         &
     &,                         LOOKUP_LBEGIN,LOOKUP_LBEGIN_NEW         &
     &,                         LOOKUP_LBNREC,LOOKUP_LBNREC_NEW         &
     &,                         FIXHD,INTHD,LOOKUP,CFI1,CFI2,CFI3,      &
     &                          WGDOS_EXPAND)

      END IF

      WRITE(6,'(I4,'' fields have been converted'')') LEN2_LOOKUP

      RETURN
      END SUBROUTINE CONVIEEE
