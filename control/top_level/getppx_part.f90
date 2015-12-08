! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Reads required portion of PPXREF file into "look-up" arrays
!
!  Subroutine Interface:

      SUBROUTINE GETPPX_PART(NFT,NFTU,StmsrNam,Im_ident,RowNumber,      &
     &                       ErrorStatus,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE filenamelength_mod, ONLY :                                    & 
          filenamelength
      USE ppxlook_mod
      USE version_mod, ONLY:                                            &
          nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
          nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
          nlevp, npslevp, npslistp, outfile_s, outfile_e

      USE Submodel_Mod
      USE stextend_mod, ONLY: IN_S

      IMPLICIT NONE
!
!  Description:
!    Reads records from PPXREF file into arrays PPXI (for integer data)
!    and PPXC (for character data, i.e. name of diagnostic/prognostic).
!    Only those ppxref records corresponding to entries in the STASH
!    addresses array IN_S are read in. Also set up pointer array PPXPTR.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards UMDP 3, version 7.
!
!  System component covered:
!  System task:               Sub-Models Project
!
!  Global Variables:
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

!  Subroutine arguments
!    Scalar arguments with intent(in):
      INTEGER      NFT,NFTU      ! Unit nos. for STASHmaster files
      CHARACTER(LEN=13)  Stmsrnam      ! Names of stash master files

!    Array arguments with intent(out):
      CHARACTER(LEN=80)   CMESSAGE    ! Error return message

!    Error status:
      INTEGER        ErrorStatus ! Error return code

!  Local scalars:
      INTEGER      I,J,K,Model   ! Loop counters
      CHARACTER(LEN=filenamelength) :: stash_mstr    
                                 ! File name for STASH master
      INTEGER      Im_index      ! Internal model index (run dependent)
      INTEGER      Im_ident      ! Internal model identifier (absolute)
      INTEGER      Section,Sec   ! section no.
      INTEGER      Item,Itm      ! item no.
      INTEGER      RowNumber     ! Row no. counter for PPXI, PPXC arrays
      INTEGER      RowNum_U      ! Do. for PPXI_U, PPXC_U (user diags.)
      CHARACTER(LEN=36) NAME
      CHARACTER(LEN=1)  CHAR1
      INTEGER      FirstBlank
      INTEGER      IOStatus

!  Local arrays:
!  WARNING: must have PPXREF_CHARLEN=4*PPX_CHARWORD
!           to avoid overwriting
      CHARACTER DNAM (PPXREF_CHARLEN) ! For char part of ppx record
      INTEGER   CODES(PPXREF_CODELEN) ! For integer part of ppx record
      INTEGER   IMASK(20)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!- End of header -------------------------------------------------------
!
      IF (lhook) CALL dr_hook('GETPPX_PART',zhook_in,zhook_handle)
      ErrorStatus = 0
      IOStatus=0
!----------------------------------------------------------------------
! Check that the no. of requested diagnostics does not exceed max
! defined in comdecks VERSION and PPXLOOK.
!
      IF ( (ppxRecs  >   NDIAGP) .OR. (ppxRecs  >   NUM_DIAG_MAX) )     &
     &THEN
        WRITE(6,*) 'ERROR: no. of diags. reqested exceeds max'
        WRITE(6,*) 'ppxRecs=',ppxRecs,' NDIAGP=',NDIAGP,                &
     &             ' NUM_DIAG_MAX=',NUM_DIAG_MAX
        Errorstatus=104
        CMESSAGE= 'GTPPXPT1: ppxRecs GT (NDIAGP or NUM_DIAG_MAX)'
        GO TO 9999
      END IF
!----------------------------------------------------------------------
!Open STASHmaster file for current internal model
!  Get directory name for STASHmaster & append rest of filename
      CALL GET_FILE(NFT,STASH_MSTR,filenamelength,ErrorStatus)
      FirstBlank = 0
      DO I = 1,filenamelength
        IF (STASH_MSTR(I:I) == ' '.AND.FirstBlank == 0)                 &
     &                                   FirstBlank=I
      END DO
      STASH_MSTR(FirstBlank:FirstBlank)='/'
      STASH_MSTR(FirstBlank+1:FirstBlank+13)=StmsrNam
      OPEN(UNIT=NFT,FILE=STASH_MSTR,IOSTAT=IOStatus)

      IF(IOStatus /= 0) THEN
        WRITE(6,*) 'ERROR in routine GETPPX_PART'
        WRITE(6,*)                                                      &
     & 'CANNOT OPEN STASHmaster FILE, IOSTATUS=',IOStatus
        WRITE(6,*) 'UNIT=',NFT,' FILE=',STASH_MSTR
        ErrorStatus=100
        CMESSAGE=' GETPPX_PART: ERROR OPENING STASHmaster'
        GOTO 9999
      END IF

! Read the required ppxref records into PPXI, PPXC
      Im_index    = INTERNAL_MODEL_INDEX(Im_ident)
      DO Section  = 0,PPXREF_SECTIONS
        DO Item   = 1,PPXREF_ITEMS

! Check whether there is a stash entry
          IF (IN_S(1,Im_ident,Section,Item)  /=  0) THEN
! Assign pointer value
            PPXPTR(Im_index,Section,Item) = RowNumber

!  OriginFlag was compressed down at end of STASH_PROC,
!  to contain only those items requested.

            IF (OriginFlag(RowNumber) == 'U') THEN
!  Record is from user STASHmaster

!  GETPPX saved all userSTASHmaster records, not just
!  those requested, so search for correct record.
              DO  I = 1,NUM_USR_DIAG_MAX
                IF (PPXI_U(I,1) == Im_ident .and.                       &
     &              PPXI_U(I,2) == Section .and.                        &
     &              PPXI_U(I,3) == Item) THEN
!  Correct record found
                  RowNum_U = I
                END IF
              END DO
! Read user ppxref record from transfer arrays
              DO I=1,PPXREF_CHARLEN
                PPXC(RowNumber,I)=PPXC_U(RowNum_U,I)
              END DO
              DO I=1,PPXREF_CODELEN
                PPXI(RowNumber,I)=PPXI_U(RowNum_U,I)
              END DO
              IF ((PPXI(RowNumber,1) /= Im_ident).OR.                   &
     &            (PPXI(RowNumber,2) /= Section ).OR.                   &
     &            (PPXI(RowNumber,3) /= Item    )) THEN
                WRITE(6,*) 'ERROR, GETPPX_PART: '
                WRITE(6,*) 'Inconsistency in user ppxref transfer'
                WRITE(6,*) 'Model,Section,Item: ',                      &
     &                      Im_ident,Section,Item
                ErrorStatus=115
                GO TO 9999
              END IF

            ELSE IF (OriginFlag(RowNumber) == 'P') THEN

! Find appropriate record in STASHmaster file and read it in
 100          READ(NFT,'(A1)') CHAR1

              IF (CHAR1 == '1') THEN
                BACKSPACE NFT
                READ(NFT,'(2X,3(I5,2X))') Model,Sec,Itm

                IF (Model == -1) THEN
                  WRITE(6,*)                                            &
     &           'GETPPX_PART: End of STASHmaster file ',               &
     &            StmsrNam,' reached'
                  GO TO 1100
                END IF
                IF (Sec == Section .AND. Itm == Item) THEN
!   Correct record found
                  BACKSPACE NFT
! DEPENDS ON: readstm
                  CALL READSTM                                          &
     &           (IMASK,DNAM,CODES,NFT,ErrorStatus,CMESSAGE)
!   Transfer STASHmaster record to look-up arrays
                  DO I=1,PPXREF_CHARLEN
                    PPXC(RowNumber,I)=DNAM(I)
                  END DO
                  DO I=1,PPXREF_CODELEN
                    PPXI(RowNumber,I)=CODES(I)
                  END DO
                ELSE
                  GO TO 100
                END IF
              ELSE
                GO TO 100
              END IF
            ELSE IF (OriginFlag(RowNumber) /= ' ') THEN
              WRITE(6,*) 'ERROR, GETPPX_PART: INVALID OriginFlag'
              WRITE(6,*) 'Row number, Flag'
              WRITE(6,*) RowNumber, OriginFlag(RowNumber)
                ErrorStatus=135
                GO TO 9999
            END IF

            RowNumber = RowNumber + 1
 1100       CONTINUE

            IF ((RowNumber-1)  >   ppxRecs) THEN
              WRITE(6,*) 'Error in GETPPX_PART:'
              WRITE(6,*)                                                &
     &       ' PPXI row number exceeds total no. of ppx records'
              GO TO 9999
            END IF

          END IF   ! Stash entries
        END DO     ! Items
      END DO     ! Sections

 9999 CONTINUE

      CLOSE(UNIT=NFT)
      IF (lhook) CALL dr_hook('GETPPX_PART',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GETPPX_PART
