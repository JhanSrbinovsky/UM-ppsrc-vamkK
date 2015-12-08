! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Reads PPXREF file into "look-up" arrays
!
!  Subroutine Interface:

      SUBROUTINE GETPPX(NFTPPXREF,NFTSTMSTU,StmsrNam,RowNumber,         &
     &                       ErrorStatus,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars, ONLY: mype 
      USE filenamelength_mod, ONLY :                                    & 
          filenamelength
      USE nl_ustsnum_mod, ONLY : n_uSTASH, nrecs_uSTASH, ustsfils
      USE ppxlook_mod 
      USE cppxref_mod, ONLY:                                            &
          ppxref_codelen, ppxref_charlen,                               &
          ppx_model_number, ppx_section_number, ppx_item_number
! version_mod items required by cstash.h
      USE version_mod, ONLY:                                            &
          nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
          nlevp, npslevp, npslistp, outfile_s, outfile_e
      USE Submodel_Mod

      IMPLICIT NONE
!
!  Description:
!    Reads records from PPXREF file into arrays PPXI (for integer data)
!    and PPXC (for character data, i.e. name of diagnostic/prognostic).
!    The entire PPXREF file is read in (non-null records only).
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
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
      INTEGER        NFTPPXREF   ! Unit no. for PPXREF file
      INTEGER        NFTSTMSTU   ! Unit no. for user ppxref files
      CHARACTER(LEN=13)    StmsrNam    ! Names of stash master files

!    Array arguments with intent(out):
      CHARACTER(LEN=80) CMESSAGE    ! Error return message

!    Error status:
      INTEGER      ErrorStatus ! Error return code

!  Local scalars:
      INTEGER      I,J,IE,ID,II  ! Loop counters
      INTEGER      hashcount
      INTEGER      IFIL,IREC     ! Do.
      INTEGER      IOSTATUS
      CHARACTER(LEN=filenamelength) ::  UpsmFile      
                                 ! Full pathname for user psm files
      CHARACTER(LEN=filenamelength) ::  stash_mstr    
                                 ! Do. STASH master files
      CHARACTER(LEN=*)Routine       ! Subroutine name
      PARAMETER   (Routine = 'GETPPX')
      CHARACTER(LEN=1)  CHAR1
      INTEGER      Im_index      !
      INTEGER      Im_ident      !
      INTEGER      Section       !
      INTEGER      Item          !
      INTEGER      LModel  ,DM
      INTEGER      LSection,DS
      INTEGER      LItem   ,DI
      INTEGER      USTrow
      INTEGER      RowNumber     ! Row no. counter for PPXI, PPXC arrays
      INTEGER      FirstBlank    ! Used to append Upsm file name to dir
      INTEGER      RI            ! Row index
      INTEGER      NU_recs       ! No. of records in a user psm file
      LOGICAL      OVERWRITE ! Set T if a system stash master record
                             !  is being overwritten by a user rec.
!  Local arrays:
!  WARNING: must have PPXREF_CHARLEN=4*PPX_CHARWORD
!           to avoid overwriting
      CHARACTER DNAM (PPXREF_CHARLEN) ! For character part of ppx rec
      INTEGER   CODES(PPXREF_CODELEN) ! For integer part of ppx record
      INTEGER   IMASK(20)             ! For ver mask in user psm

! FLDCALC uses a different maginc number, parameterise 
! to keep the later if clean
      INTEGER, PARAMETER :: stashunit=22

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!- End of header -------------------------------------------------------
!
      IF (lhook) CALL dr_hook('GETPPX',zhook_in,zhook_handle)
      ErrorStatus = 0
      NU_recs     = 0
      IOStatus   =0
!----------------------------------------------------------------------
! Check that the no. of requested diagnostics does not exceed max
! defined in comdecks VERSION and PPXLOOK.
!
      IF ( (ppxRecs  >   NDIAGP) .OR. (ppxRecs  >   NUM_DIAG_MAX) )     &
     &THEN
        WRITE(6,*) 'ERROR: no. of diags. requested exceeds max'
        WRITE(6,*) 'ppxRecs=',ppxRecs,' NDIAGP=',NDIAGP,                &
     &             ' NUM_DIAG_MAX=',NUM_DIAG_MAX
        Errorstatus=104
        CMESSAGE= 'GETPPX: ppxRecs GT (NDIAGP or NUM_DIAG_MAX)'
        GO TO 9999
      END IF
!----------------------------------------------------------------------

      IF (NFTPPXREF == stashunit ) THEN

!----------------------------------------------------------------------
!Read in records from STASHmaster for current internal model
!----------------------------------------------------------------------
!Open STASHmaster file for current internal model
!  Get directory name for STASHmaster & append rest of filename
        CALL GET_FILE(NFTPPXREF,STASH_MSTR,filenamelength,ErrorStatus)
        FirstBlank = 0
        DO I = 1,filenamelength
          IF (STASH_MSTR(I:I) == ' '.AND.FirstBlank == 0)               &
     &                                   FirstBlank=I
        END DO
        STASH_MSTR(FirstBlank:FirstBlank)='/'
        STASH_MSTR(FirstBlank+1:FirstBlank+13)=StmsrNam
        OPEN(UNIT=NFTPPXREF,FILE=STASH_MSTR,IOSTAT=IOStatus)
        IF(IOStatus /= 0) THEN
          WRITE(6,*) 'ERROR in routine GETPPX'
          WRITE(6,*)                                                    &
     &   'CANNOT OPEN STASHmaster FILE, IOSTATUS=',IOStatus
          WRITE(6,*) 'UNIT=',NFTPPXREF,' FILE=',STASH_MSTR
          ErrorStatus=100
          CMESSAGE=' GETPPX: ERROR OPENING STASHmaster'
          GOTO 9999
        END IF

 100    READ(NFTPPXREF,'(A1)') CHAR1

        IF (CHAR1 == '1') THEN
!Read block of records
          BACKSPACE NFTPPXREF
! DEPENDS ON: readstm
          CALL READSTM(IMASK,DNAM,CODES,NFTPPXREF,ErrorStatus,CMESSAGE)
          Im_ident = CODES(ppx_model_number)
          Section  = CODES(ppx_section_number)
          Item     = CODES(ppx_item_number)
          IF (Im_ident == -1) THEN
!End of file reached
            CLOSE(UNIT=NFTPPXREF)
            GO TO 9999
          END IF
          Im_index= INTERNAL_MODEL_INDEX(Im_ident)
!   Increment row number
          RowNumber = RowNumber + 1
! Assign value to PPXPTR element corresponding to this record
          PPXPTR(Im_index,Section,Item) = RowNumber
!   Transfer data from ppx record to look-up arrays
          DO I=1,PPXREF_CHARLEN
            PPXC(RowNumber,I)=DNAM(I)
          END DO
          DO I=1,PPXREF_CODELEN
            PPXI(RowNumber,I)=CODES(I)
          END DO
!   Set row index - indicates values of model,sec,item for this row
          RowIndex  (RowNumber)=  Im_ident*100000                       &
     &                          + Section *1000                         &
     &                          + Item
!   Set flag to indicate record originated from ppxref file
          OriginFlag(RowNumber)='P'
          IF (RowNumber  >   ppxRecs) THEN
            WRITE(6,*) 'Error in GETPPX:'
            WRITE(6,*)                                                  &
     &    ' PPXI row number exceeds total no. of ppx records ',         &
     &      RowNumber
            GO TO 9999
          END IF
          GO TO 100  ! Back to READ
        ELSE
! Skip to next line
          GO TO 100
        END IF
      ELSE         ! NFTPPXREF /= 1
! ----------------------------------------------------------
! Insert user-defined diagnostics into ppxref look-up arrays
! ----------------------------------------------------------

        IF &
      (NRECS_USTASH >  0) &
        THEN
! There are user diagnostic records
      ErrorStatus=0
      IOStatus   =0
! Get directory name for Upsm files
      CALL GET_FILE(NFTSTMSTU,UpsmFile,filenamelength,ErrorStatus)
      FirstBlank = 0
      DO I = 1,filenamelength
        IF (UpsmFile(I:I) == ' '.AND.FirstBlank == 0) FirstBlank=I
      END DO

! Loop over user pre-stash master files
      DO IFIL = 1,N_USTASH
        UpsmFile(FirstBlank  :FirstBlank  )='.'
        UpsmFile(FirstBlank+1:FirstBlank+8)=USTSFILS(IFIL)

!   Open user stash master file
        OPEN(NFTSTMSTU,FILE=UpsmFile,IOSTAT=IOStatus)
        IF(IOStatus /= 0) THEN
          WRITE(6,*) 'CANNOT OPEN USER PPXREF FILE.IOSTATUS=',          &
     &                                             IOStatus
          WRITE(6,*) 'UNIT=',NFTSTMSTU,' FILE=',UpsmFile
          ErrorStatus=100
          CMESSAGE=' GETPPX: ERROR OPENING USER PPXREF'
          GOTO 9999
        END IF

!   Read number of records in this file
        READ(NFTSTMSTU,'(I3)') NU_recs

!   Read in records from user pre-stash master file
        DO IREC = 1,NU_recs
!   Initialise OVERWRITE switch
        OVERWRITE=.FALSE.
        hashcount=0
 200    READ(NFTSTMSTU,'(A1)') CHAR1
        IF (CHAR1 /= '1') THEN
          hashcount=hashcount+1
          IF (hashcount >  20) THEN
            Errorstatus=100
            CMESSAGE='INCORRECT FORMAT IN USER STASHmaster FILE'
            WRITE(6,*) 'INCORRECT FORMAT IN USER STASHmaster FILE'
            WRITE(6,*) 'GAP BETWEEN RECORDS TOO LARGE?'
            GO TO 9999
          ELSE
            GO TO 200
          END IF
        ELSE
!Read block of records
          BACKSPACE NFTSTMSTU
! DEPENDS ON: readstm
          CALL READSTM                                                  &
     &   (IMASK,DNAM,CODES,NFTSTMSTU,ErrorStatus,CMESSAGE)
          Im_ident = CODES(ppx_model_number)
          Section  = CODES(ppx_section_number)
          Item     = CODES(ppx_item_number)

!   Transfer data from ppx record to look-up arrays
!   No. of records extracted from STASHmaster file(s)= RowNumber.
          USTrow    =   0
          DO I=1,RowNumber
            RI      =   RowIndex(I)
            IF(RI <= 0) THEN   ! Check valid Rowindex
               Errorstatus=-1
               Cmessage= Routine//                                      &
     & ':Warning, invalid Rowindex for user STASHmaster record'
               write(6,*) Cmessage
               write(6,*) 'im_ident section item Rowindex=',            &
     &                     im_ident,section,item,RI
            ENDIF   ! Check valid Rowindex

!     Determine values of model,section,item for this row
            IF (USTrow >  0) THEN  ! Position of record found so
                exit               ! exit from I loop over RowNumber
            ENDIF

              LModel  =     RI/100000
              LSection=(RI-(RI/100000)*100000)/1000
              LItem   =(RI-(RI/1000  )*1000  )
!     Check whether previous item is being overwritten
              IF (Im_ident == LModel  .AND.                             &
     &            Section  == LSection.AND.                             &
     &            Item     == LItem        ) THEN
                IF      (OriginFlag(I) == 'P') THEN
                  OVERWRITE=.TRUE.
! Send userSTASH messages both to stdout and to the .stash file
                  IF (mype == 0) THEN
                    WRITE(200,'(A,A)')                                    &
             'A user-STASHmaster file has ', &
             'overwritten the following record:'
                    WRITE(200,'(A,I2,A,I4,A,I4)') '  Internal Model ',  &
                      Im_ident, ' Section ',Section,' Item ',Item
                  END IF
                  WRITE(6,'(A,A)')                                        &
             'A user-STASHmaster file has ',&
             'overwritten the following record:'
                  WRITE(6,'(A,I2,A,I4,A,I4)') '  Internal Model '  ,    &
                  Im_ident, ' Section ',Section,' Item ',Item
                ELSE IF (OriginFlag(I) == 'U') THEN
                  WRITE(6,'(A)') 'ERROR, GETPPX: '
                  WRITE(6,'(A)') 'User diagnostic duplicated'
                  WRITE(6,'(A,I2,I4,I4)') 'Model,Section,Item ',        &
     &                        Im_ident,Section,Item
                  ErrorStatus=100
                  CMESSAGE='ERROR,GETPPX:user diag duplicated'
                  GO TO 9999
                END IF
              END IF
!     Determine appropriate row number
              IF (LModel   == Im_ident.AND.                             &
     &            LSection == Section .AND.                             &
     &            LItem    == Item    .AND.USTrow == 0) THEN
                USTrow=I    ! Row number found
!     This record will overwrite a pre-existing record
!     Insert new record
                DO IE=1,PPXREF_CHARLEN
                  PPXC(USTrow,IE)=DNAM(IE)
                END DO
                DO IE=1,PPXREF_CODELEN
                  PPXI(USTRow,IE)=CODES(IE)
                END DO
!     Set flag to indicate record originated from user psm file
                OriginFlag(USTrow)='U'
              ELSE IF((LModel   >  Im_ident.AND.USTrow == 0) .OR.       &
     &                (LModel   == Im_ident.AND.                        &
     &                 LSection >  Section .AND.USTrow == 0) .OR.       &
     &                (LModel   == Im_ident.AND.                        &
     &                 LSection == Section .AND.                        &
     &                 LItem    >  Item    .AND.USTrow == 0)) THEN
                USTrow=I    ! Row number found
!     This record will be inserted between two pre-existing records
!     Create spare row - move all subsequent records up by one row
                DO ID = RowNumber+1,USTrow+1,-1
                  DO IE=1,PPXREF_CHARLEN
                    PPXC(ID,IE)=PPXC(ID-1,IE)
                  END DO
                  DO IE=1,PPXREF_CODELEN
                    PPXI(ID,IE)=PPXI(ID-1,IE)
                  END DO
                  RI            =RowIndex  (ID-1)
                  RowIndex  (ID)=RowIndex  (ID-1)
                  OriginFlag(ID)=OriginFlag(ID-1)
!     Determine values of model,section,item for this row
                  DM=     RI/100000
                  DS=(RI-(RI/100000)*100000)/1000
                  DI=(RI-(RI/1000  )*1000  )
!     Increment PPXPTR for record moved up
                  Im_index=INTERNAL_MODEL_INDEX(DM)
                  PPXPTR(Im_index,DS,DI)=PPXPTR(Im_index,DS,DI)+1
                END DO
!     Insert new record
                DO IE=1,PPXREF_CHARLEN
                  PPXC(USTrow,IE)=DNAM(IE)
                END DO
                DO IE=1,PPXREF_CODELEN
                  PPXI(USTRow,IE)=CODES(IE)
                END DO
!     Set row index - indicates model,sec,item for this row
                RowIndex  (USTrow)=  Im_ident*100000                    &
     &                             + Section *1000                      &
     &                             + Item
!     Set flag to indicate record originated from user psm file
                OriginFlag(USTrow)='U'
!     Set PPXPTR for the new record
                Im_index=INTERNAL_MODEL_INDEX(Im_ident)
                PPXPTR(Im_index,Section,Item)=USTrow


            ELSE IF (I == RowNumber) THEN ! last record and no match
                                          ! so append
!     This record will be added after all pre-existing records
              USTrow = I + 1      ! Set extra entry beyond current last
!     Add new record
              DO IE=1,PPXREF_CHARLEN
                PPXC(USTrow,IE)=DNAM(IE)
              END DO
              DO IE=1,PPXREF_CODELEN
                PPXI(USTrow,IE)=CODES(IE)
              END DO
!     Set row index - indicates model,sec,item for this row
              RowIndex  (USTrow)=  Im_ident*100000                      &
     &                           + Section *1000                        &
     &                           + Item
!     Set flag to indicate record originated from user psm file
              OriginFlag(USTrow)='U'
!     Set PPXPTR for the new record
              Im_index=INTERNAL_MODEL_INDEX(Im_ident)
              PPXPTR(Im_index,Section,Item)=USTrow
            END IF

          END DO  ! I=1,RowNumber : loop over current set of records

!     Increment RowNumber as UserSTASH record has been added.
!     don't increment it if a standard record has been overwritten.
        IF (.NOT.OVERWRITE) THEN
        RowNumber = RowNumber + 1
        END IF
        END IF        ! hashcount
        END DO        ! Loop over IREC recs in upsm file
      END DO          ! Loop over user psm files
      END IF          ! NRECS_USTASH >  0
! Copy user pre-stash master records to storage arrays -
!   for passing into to U_MODEL
! Note: OriginFlag will be compressed to requested items only
!  at the end of routine STASH_PROC (before used in GETPPX_PART)
      IF (NRECS_USTASH >  0) THEN
      RowNumber = 1
      DO I = 1,ppxRecs
        IF (OriginFlag(I) == 'U') THEN
          DO IE=1,PPXREF_CHARLEN
            PPXC_U(RowNumber,IE)=PPXC(I,IE)
          END DO
          DO IE=1,PPXREF_CODELEN
            PPXI_U(RowNumber,IE)=PPXI(I,IE)
          END DO
          RowNumber=RowNumber+1
        END IF
      END DO
      END IF

      END IF  !NFT (Standard STASHmstr or user STASHmstr)

 9999 CONTINUE
      IF (lhook) CALL dr_hook('GETPPX',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE GETPPX
