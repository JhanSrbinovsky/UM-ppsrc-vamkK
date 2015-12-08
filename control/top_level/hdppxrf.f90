! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL   SUBROUTINE HDPPXRF ---------------------------------------------
!LL
!LL   PROGRAM TO READ THE HEADER RECORD OF THE PPXREF FILE
!LL   CHECK THE VALUES AND RETURN THE FILE DIMENSIONS
!LL
!LL
!LL   PROGRAMMING STANDARD  UMDP 3
!LL
!LL   LOGICAL COMPONENT  R911
!LL
!LL   PROJECT TASK: C4
!LL
!LL   EXTERNAL DOCUMENT C4
!LL
!LLEND---------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Top Level

      SUBROUTINE HDPPXRF(NFT,StmsrNam,ICODE,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars, ONLY: mype
      USE application_description, ONLY : isSmallExec
      USE filenamelength_mod, ONLY :                                    & 
          filenamelength
      USE version_mod, ONLY: nprofdp, nlevlstsp, npslevp
      USE nl_ustsnum_mod, ONLY : &
             USTSNUM, n_uSTASH, nrecs_uSTASH, ustsfils
      USE Submodel_Mod

      USE ppxlook_mod, ONLY: ppxrecs
      USE cppxref_mod, ONLY: ppx_recordlen
! version_mod items required by cstash.h
      USE version_mod, ONLY:                                            &
          nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
          nlevp, npslevp, npslistp, outfile_s, outfile_e

      IMPLICIT NONE
      INTEGER NFT                    !IN:  UNIT NUMBER FOR FILE
      CHARACTER(LEN=80) CMESSAGE        !OUT: ERROR RETURN MESSAGE
      INTEGER ICODE                  !OUT: ERROR RETURN CODE

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

! Local parameters:
      INTEGER, PARAMETER ::                                             &
        no_version = -9999 ! Value returned by function get_umversion if
                           ! environment variable VN not set

! Local Scalars
      INTEGER :: get_um_version
      INTEGER LEN_IO
      INTEGER IU          ! Local - unit no. for stash control file
      INTEGER I
      INTEGER Int_Model_No
      INTEGER FirstBlank
      CHARACTER(LEN=13)  StmsrNam
      CHARACTER(LEN=filenamelength) ::  stash_mstr
      CHARACTER(LEN=1)  CHAR1
      INTEGER      IOStatus

      CHARACTER(LEN=8) c_um_version  !UM version as string
      CHARACTER(LEN=8) c_stm_version !STASHmaster version string
      INTEGER :: um_version     !Version of UM
      INTEGER  :: um_revision   !Revision of UM
      INTEGER :: stm_version    !Version of STASHmaster file
      INTEGER :: stm_revision   !Revision of STASHmaster file
      INTEGER   ocode           !Copy of the input value of ICODE
      LOGICAL found_version     !Indicates presence of STM version
      REAL STATUS
      INTEGER RECORD(PPX_RECORDLEN)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('HDPPXRF',zhook_in,zhook_handle)
      IOStatus=0
!    Check if an error has already been encountered, and get out
!    if it has.
      ocode = 0
      IF (icode  >   0) THEN
         GO TO 9999
      ELSE IF (icode  <   0) THEN
         ocode = icode
         icode = 0
      END IF

      IF &
      (NFT == 22) &
      THEN

!Open STASHmaster file for current internal model
!  Get directory name for STASHmaster & append rest of filename
        CALL GET_FILE(NFT,STASH_MSTR,filenamelength,ICODE)
        FirstBlank = 0
        DO I = 1,filenamelength
          IF (STASH_MSTR(I:I) == ' '.AND.FirstBlank == 0)               &
                                         FirstBlank=I
        END DO
        STASH_MSTR(FirstBlank:FirstBlank)='/'
        STASH_MSTR(FirstBlank+1:FirstBlank+13)=StmsrNam
        OPEN(UNIT=NFT,FILE=STASH_MSTR,IOSTAT=IOStatus)
        WRITE(6,'(A,A)') 'STASH_MSTR: ',STASH_MSTR
        IF(mype == 0 .AND. .NOT.isSmallExec() ) THEN
! Duplicate STASHmaster path to the .stash file for pe0
          WRITE(200,'(A,A)') 'STASH_MSTR: ',STASH_MSTR
        END IF
        
        IF (IOStatus /= 0) THEN
          CMESSAGE= 'Error opening STASHmaster file, routine HDPPXRF'
          WRITE(6,'(A,I4,A,A)')                                         &
         'HDPPXRF: Fortran Error Response = ',IOStatus,                 &
         ' Opening STASHmaster file ',StmsrNam
          ICODE=100
          GO TO 9999
        ENDIF

!    Get the UM version from the environment variable $VN.
! DEPENDS ON: get_um_version
      um_version = get_um_version()

      IF ( um_version  /=  no_version ) THEN


!     Now check through the header section of the STASHmaster
!     file looking for H3
      found_version = .FALSE.
      READ (nft, '(A1)') char1
      DO WHILE (char1  ==  'H' .or. char1  ==  '#')
         IF (char1  ==  'H') THEN
            BACKSPACE nft
            READ (nft, '(1X, A1)') char1
            IF (char1  ==  '3') THEN
!     This line starts with H3 and should
!     indicate the STASHmaster version. The line should look like
!     H3| UM_VERSION=X.Y
!     where X.Y is the UM version the STASHmaster is valid for.
               found_version = .TRUE.
               BACKSPACE nft
               READ (nft, '(15x,a8)') c_stm_version
               READ (c_stm_version, '(i1,1x,i1)')                       &
                    stm_version, stm_revision
               stm_version = stm_version*100 + stm_revision
!     Now perform the check against the UM version
               IF (stm_version  /=  um_version) THEN
                  WRITE (cmessage,'(A)')                                &
       'HDPPXRF : UM version and STASHmaster version differ'
                  WRITE (6,'(A,I4,A,I4,A,A)')                           &
                       'Version of STASHmaster file ('                  &
                       ,stm_version,                                    &
                       ') does not match UM version ('                  &
                       ,um_version,') in file ',StmsrNam
                  icode = 1
                  GO TO 9999
               END IF  ! version check
            END IF  ! char1 == '3'
         END IF                 ! char1 == 'H'
         READ (nft, '(A1)') char1
      END DO

      IF (.NOT. found_version) THEN
         WRITE (6,'(A)')                                                &
              'HDPPXRF : No STASHmaster version available; Unable to'
         WRITE (6,'(A,A)')                                              &
              'check against UM version for file ',StmsrNam
         cmessage = 'HDPPXRF : No STASHmaster version available'
         icode = -1
      END IF
!     For safety, rewind to the start of the STASH file.
      REWIND (nft)

      END IF
100   CONTINUE
!Count records - ppxRecs is counter
        READ(NFT,'(A1)') CHAR1
        IF (CHAR1 == '1') THEN
          BACKSPACE NFT
          READ(NFT,'(2X,I5)') Int_Model_No
          IF (Int_Model_No == -1) THEN
!End of file reached
!ppxRecs initialised to 1 before HDPPXRF - so subtract 1 now
            IF (StmsrNam(13:) == 'A') THEN
              IF (INTERNAL_MODEL_INDEX(A_IM) == 1) THEN
                ppxRecs=ppxRecs-1
              END IF
            END IF
            CLOSE(UNIT=NFT)
            GO TO 9999
          END IF
          ppxRecs = ppxRecs + 1
          GO TO 100
        ELSE
          GO TO 100
        END IF
      ELSE

! Read USTSNUM namelist from unit 5:
! Number of user stash files and total no. of user stash records
        IU = 5

!Initialisation
        nrecs_uSTASH = 0
        ustsfils = '        '
        n_uSTASH = 0
! Read namelist
        READ(IU,USTSNUM)
! Add no. of user stash records to ppxRecs
        ppxRecs = ppxRecs + nrecs_uSTASH
      END IF

 9999 CONTINUE
!     If we have found an error, leave it in icode.  If no error
!     occurred then check if the original input value of icode was
!     non-zero (a previous untrapped error/warning), and copy this
!     back into ICODE before eaving the routine.
      IF (icode  ==  0 .and. ocode  /=  0) THEN
         icode = ocode
      END IF
      IF (lhook) CALL dr_hook('HDPPXRF',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE HDPPXRF
