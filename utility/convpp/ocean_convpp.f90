! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE OCEAN_CONVPP-----------------------------------------
!LL
!LL Purpose: Converts UM ocean file to PP format.
!LL
!LL  Documentation: UM Doc Paper F5
!LL
!LL  System Tasks: F3,F4,F6
!LL
!LL  -----------------------------------------------------------------
!*L  Arguments:-------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      SUBROUTINE OCEAN_CONVPP                                           &
     &  (LEN_FIXHD,LEN_INTHD,LEN_REALHD,                                &
     &  LEN1_LEVDEPC,LEN2_LEVDEPC,LEN1_ROWDEPC,                         &
     &  LEN2_ROWDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                         &
     &  LEN1_FLDDEPC,LEN2_FLDDEPC,LEN_EXTCNST,                          &
     &  LEN_DUMPHIST,LEN_CFI1,LEN_CFI2,LEN_CFI3,                        &
     &  LEN1_LOOKUP,LEN2_LOOKUP,LEN_DATA,P_FIELD,                       &
     &  NFTIN,MAX_FIELD_SIZE,JOC_NO_SEAPTS,LEN_OCFLD)
!L
!L
      USE PP_HEADER_MANIPS, ONLY: HEADER_MANIP

      USE ereport_mod, ONLY : ereport
      USE Decomp_DB
      USE UM_ParVars
      USE ppxlook_mod, ONLY : ppxrecs
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
     &,MAX_FIELD_SIZE                                                   &
                      !IN Maximum field size on file
     &,JOC_NO_SEAPTS                                                    &
                      !IN Number of points in compressed array
     &,LEN_OCFLD      !IN Length of uncompressed ocean field

      INTEGER                                                           &
     & NFTIN


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
     & LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP)                                  &
                                                 !
     &,LOOKUP_OUT(LEN1_LOOKUP) ! Output lookup table

! Input arguments for decompose_smexe
      INTEGER                                                           &

     &  global_row_len,                                                 &
                          ! IN  :number of E-W points of entire model
     &  global_n_rows,                                                  &
                          ! IN  :number of P rows of entire mode
     &  tot_levels        ! IN  :total number of levels

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
     & DUMPHIST(LEN_DUMPHIST+1),                                        &
                                                 !
     & D1(MAX_FIELD_SIZE),                                              &
                           ! Array used to read in non-compressed fields
     & E1(MAX_FIELD_SIZE),                                              &
                           ! Array used to read in non-compressed fields
                           ! without wrap points
     & C1(JOC_NO_SEAPTS),                                               &
                           ! Array used to read in compressed fields
     & U1(LEN_OCFLD)       ! Array used to hold  uncompressed fields


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
!*----------------------------------------------------------------------
! Local Constants
      INTEGER, PARAMETER :: lenPpHead = 64, lenIntHead = 45

!*L  Local variables:---------------------------------------------------

      INTEGER                                                           &
     & ICODE                                                            &
                     ! Error return code from subroutines
     &,START_BLOCK                                                      &
                     ! READHEAD argument (not used)
     &,I,J,K,L                                                          &
                     ! Loop indices
     &,IJ_IN,IJ_OUT                                                     &
                     ! More loop indices
     &,NROWS                                                            &
                     ! Number of points north-south
     &,IROWS                                                            &
     &,NROWS_FIELD                                                      &
                     ! Number of rows in a field
     &,NCOLS_IN                                                         &
                     ! Number of points east-west
     &,ICOLS_IN                                                         &
     &,NCOLS_OUT                                                        &
                     ! Number of points east-west for pp fields
     &,NLEVS                                                            &
                     ! Number of points in vertical
     &,NT                                                               &
                     ! Number of tracers
     &,NCOMP                                                            &
                     ! Number of compressed fields
     &,RECNUM                                                           &
                     ! Record number of field in lookup table
     &,POSIN                                                            &
                     ! Start position of field within C1
     &,POSU1                                                            &
                     ! Start position of field within U1
     &,FIELD_CODE                                                       &
                     ! field code for this field
     &,LBPACK        ! packing indicator from lookup table

      CHARACTER                                                         &
     & CMESSAGE*100                                                     &
                    ! Character string returned if ICODE  /=  0
     &,STRING*20    ! Format control for header printout

      REAL                                                              &
     & RMDI         ! Real missing data indicator

      LOGICAL                                                           &
     & LL_CYCLIC_IN    ! T => cyclic ; F => not cyclic

      INTEGER RowNumber

      INTEGER NFT1,NFT2
      PARAMETER (NFT1=22, NFT2=2)

!*----------------------------------------------------------------------

!L 0. Read in PPXREF

      ppxRecs=1
      RowNumber=0
      CMESSAGE=''
      ICODE=0
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_A',ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_A'
        WRITE(6,*) CMESSAGE

        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   'Error reading STASHmaster_A')
      END IF

! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_O',ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_O'
        WRITE(6,*) CMESSAGE

        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   'Error reading STASHmaster_O')
      END IF

! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_A',RowNumber,                  &
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE

        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   CMESSAGE)
      ENDIF

! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_O',RowNumber,                  &
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE

        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   CMESSAGE)
      ENDIF

!User STASHmaster
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(0,'             ',ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE

        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   CMESSAGE)
      ENDIF

! DEPENDS ON: getppx
      CALL GETPPX(0,NFT2,'             ',RowNumber,                     &
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE

        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   CMESSAGE)
      ENDIF

!L 1. Read in file header

! DEPENDS ON: readhead
      CALL READHEAD(NFTIN,FIXHD,LEN_FIXHD,                              &
     &                INTHD,LEN_INTHD,                                  &
     &                REALHD,LEN_REALHD,                                &
     &                LEVDEPC,LEN1_LEVDEPC,LEN2_LEVDEPC,                &
     &                ROWDEPC,LEN1_ROWDEPC,LEN2_ROWDEPC,                &
     &                COLDEPC,LEN1_COLDEPC,LEN2_COLDEPC,                &
     &                FLDDEPC,LEN1_FLDDEPC,LEN2_FLDDEPC,                &
     &                EXTCNST,LEN_EXTCNST,                              &
     &                DUMPHIST,LEN_DUMPHIST,                            &
     &                CFI1,LEN_CFI1,                                    &
     &                CFI2,LEN_CFI2,                                    &
     &                CFI3,LEN_CFI3,                                    &
     &                LOOKUP,LEN1_LOOKUP,LEN2_LOOKUP,                   &
     &                LEN_DATA,                                         &
     &                START_BLOCK,ICODE,CMESSAGE)

      IF(ICODE /= 0)THEN
        WRITE(6,*)CMESSAGE,ICODE

        CALL EREPORT('OCEAN_CONVPP', ICODE,                             &
     &   CMESSAGE)
      ENDIF


      NROWS        = INTHD(7)
      NCOLS_IN     = INTHD(6)
      IROWS        = INTHD(7)
      ICOLS_IN     = INTHD(6)
      LBPACK       = 21
      NLEVS        = INTHD(8)
      NT           = INTHD(14)

      IF(LEN_REALHD >= 29)THEN
        RMDI         = REALHD(29)
      ENDIF

! Get decomposition information
! -----------------------------
      TOT_LEVELS=-99

      CALL decompose(ICOLS_IN,IROWS,0,0,TOT_LEVELS)
      CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)

! Determine whether input data is cyclic and number of columns to output

      IF ( MOD ( FIXHD(4), 100 )  /=  3 ) THEN
        LL_CYCLIC_IN = .TRUE.
      ELSE
        LL_CYCLIC_IN = .FALSE.
      ENDIF

!L 2. Read in compressed data

      RECNUM=1

! Decide whether there are any compressed fields and on number of
! compressed fields. Use LBPACK to work out whether the first field
! contains sea points only.

      IF ( MOD(LOOKUP(LBPACK,1)/10,10)  ==  0) THEN

       NCOMP = 0

      ELSE

       NCOMP = NT + 2

      DO L=1,NCOMP

! Loop over levels storing all levels in one 1-D array
        POSIN=1
        DO K=1,NLEVS

! DEPENDS ON: readflds
          CALL READFLDS(NFTIN,1,RECNUM,LOOKUP,LEN1_LOOKUP,C1(POSIN),    &
     &                  MAX_FIELD_SIZE,FIXHD,                           &
     &                  ICODE,CMESSAGE)
! DEPENDS ON: abort_io
          IF(ICODE /= 0)CALL ABORT_IO('CONVPP',CMESSAGE,ICODE,NFTIN)
          POSIN=POSIN+LOOKUP(15,K+(L-1)*NLEVS)
          RECNUM=RECNUM+1

        ENDDO

!L 3. Uncompress 3-D field
! DEPENDS ON: unpack
        CALL UNPACK(1,NROWS,1,NLEVS,NROWS,NLEVS,NCOLS_IN,NROWS,NLEVS,   &
     &          CFI1,CFI2,LEN_CFI1,CFI3,JOC_NO_SEAPTS,                  &
     &          C1,U1,RMDI,LL_CYCLIC_IN)

!L 4. Output data level by level
        DO K=1,NLEVS

! Fill output lookup table
          DO I=1,LEN1_LOOKUP
            LOOKUP_OUT(I)=LOOKUP(I,K+(L-1)*NLEVS)
          ENDDO

          FIELD_CODE = LOOKUP_OUT(23)

          IF (  FIELD_CODE  >   600 .AND. FIELD_CODE  <   700) THEN
            NROWS_FIELD = NROWS
          ELSE IF ( FIELD_CODE  >   699 .AND. FIELD_CODE  <   800 ) THEN
            NROWS_FIELD = NROWS - 1
          ELSE
            write(6,*) ' unknown field code : exiting '
            go to 9999
          END IF

! Determine number of columns to output

          IF ( LL_CYCLIC_IN ) THEN
            NCOLS_OUT = NCOLS_IN - 2
          ELSE
            NCOLS_OUT = NCOLS_IN
          ENDIF

          LOOKUP_OUT(15)=NROWS_FIELD*NCOLS_OUT
          LOOKUP_OUT(18)=NROWS_FIELD
          LOOKUP_OUT(19)=NCOLS_OUT
          LOOKUP_OUT(21)=0
! Note: call to header_manip included for consistency  - don't
!       think it should be used.
          CALL HEADER_MANIP(LOOKUP_OUT(1:lenIntHead))
          WRITE(10)(LOOKUP_OUT(I),I=1,lenPpHead)

          POSU1=(K-1)*NROWS*NCOLS_IN
          DO J=1,NROWS_FIELD
            DO I=1,NCOLS_OUT
              IJ_IN = I + (J-1) * NCOLS_IN
              IJ_OUT   = I + (J-1) * NCOLS_OUT
              E1(IJ_OUT) = U1(IJ_IN+POSU1)
            ENDDO
          ENDDO


          WRITE(10) (E1(I),I=1,NROWS_FIELD*NCOLS_OUT)

        ENDDO

      ENDDO

      END IF  ! LBPACK


!L 5.  Now processing non compressed fields
!  Print out individual fields
      DO L=RECNUM,LEN2_LOOKUP

        IF(LOOKUP(1,L) == -99)GOTO 100

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN,1,L,LOOKUP,LEN1_LOOKUP,                     &
     &                D1,MAX_FIELD_SIZE,FIXHD,                          &
     &                ICODE,CMESSAGE)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('CONVPP',CMESSAGE,ICODE,NFTIN)
        IF(FIXHD(5) /= 3)THEN

! Take off the extra columns if the dump is cyclic using the E1 array

          IF ( LL_CYCLIC_IN ) THEN
            NCOLS_OUT = NCOLS_IN - 2
          ELSE
            NCOLS_OUT = NCOLS_IN
          ENDIF

          DO J=1,NROWS
            DO I=1,NCOLS_OUT
              IJ_IN = I + (J-1) * NCOLS_IN
              IJ_OUT   = I + (J-1) * NCOLS_OUT
              E1(IJ_OUT) = D1(IJ_IN)
            ENDDO
          ENDDO

        ELSE

! Fieldsfile. NO cyclic columns

          DO I=1,LOOKUP(15,L)
            E1(I) = D1(I)
          ENDDO

        ENDIF
        DO K=1,LEN1_LOOKUP
          LOOKUP_OUT(K)=LOOKUP(K,L)
        ENDDO

        IF(FIXHD(5) /= 3)THEN
          LOOKUP_OUT(15)=NROWS*NCOLS_OUT
          LOOKUP_OUT(19)=NCOLS_OUT
        ENDIF
        LOOKUP_OUT(21)=MOD(LOOKUP_OUT(21),1000)

        CALL HEADER_MANIP(LOOKUP_OUT(1:lenIntHead)) 
        WRITE(10)(LOOKUP_OUT(K),K=1,lenPpHead)
        WRITE(10) (E1(K),K=1,LOOKUP_OUT(15))

      ENDDO

 100  CONTINUE
      WRITE(6,*)L-1,' pp fields written out'

9999  continue

      RETURN
      END SUBROUTINE OCEAN_CONVPP

