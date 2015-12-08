! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!    Program  MAIN_COMPARE and Subroutine COMPARE
!  
!    Purpose: Compares two UM atmosphere, ocean, or ancillary files.
!             MAIN_COMPARE reads in fixed length and integer
!             headers of UM files to be compared, extracts dimensions
!             of each file and then passes these values to
!             subroutine COMPARE.
!  
!              COMPARE subroutine:
!            Compares two UM atmosphere, ocean, or ancillary files.
!            COMPARE reads in headers and data fields from files on
!            NFTIN1 and NFTIN2, comparing values.
!            UNIT 6: If an exact compare is found the message 'OK'
!            is written out, otherwise
!            i)  if header, all differring values are printed
!            ii) if field, 1st 10 differring values are printed plus
!                the maximum difference between the fields.
!            iii) if field only present in one file, a warning message
!                 is displayed
!            UNIT 7: Number of differences displayed for each header.
!                    Number of fields with differences is also
!                    displayed along with the number of differences
!                    for each field which has differences
!  
!  
!    Programming standard:
!  
!    Logical components covered:
!  
!    System Tasks: F3,F4,F6
!  
!    Documentation: UM Doc Paper F5
!  
!    -----------------------------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Small execs
      PROGRAM MAIN_COMPARE
        USE IO
        USE ereport_mod, ONLY : ereport, ereport_finalise
        USE UM_Config, ONLY : &
            appInit,          &
            exe_cumf
! version_mod items required by cstash.h
        USE version_mod, ONLY :                                        &
            nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,      &
            nlevp, npslevp, npslistp, outfile_s, outfile_e

        USE Submodel_Mod
        
        USE filenamelength_mod, ONLY:                                     &
            filenamelength

        IMPLICIT NONE

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
      INTEGER                                                           &
     & FIXHD1(256)                                                      &
                         !Space for fixed length header file 1
     &,INTHD1(100)       !Space for integer header file 1

      INTEGER                                                           &
     & FIXHD2(256)                                                      &
                         !Space for fixed length header file 2
     &,INTHD2(100)       !Space for integer header file 2

      INTEGER                                                           &
     & LEN_FIXHD1                                                       &
                      !Length of fixed length header on file 1
     &,LEN_INTHD1                                                       &
                      !Length of integer header on file 1
     &,LEN_REALHD1                                                      &
                      !Length of real header on file 1
     &,LEN1_LEVDEPC1                                                    &
                      !1st dim of lev dependent consts on file 1
     &,LEN2_LEVDEPC1                                                    &
                      !2nd dim of lev dependent consts on file 1
     &,LEN1_ROWDEPC1                                                    &
                      !1st dim of row dependent consts on file 1
     &,LEN2_ROWDEPC1                                                    &
                      !2nd dim of row dependent consts on file 1
     &,LEN1_COLDEPC1                                                    &
                      !1st dim of col dependent consts on file 1
     &,LEN2_COLDEPC1                                                    &
                      !2nd dim of col dependent consts on file 1
     &,LEN1_FLDDEPC1                                                    &
                      !1st dim of field dependent consts on file 1
     &,LEN2_FLDDEPC1                                                    &
                      !2nd dim of field dependent consts on file 1
     &,LEN_EXTCNST1                                                     &
                      !Length of extra consts on file 1
     &,LEN_DUMPHIST1                                                    &
                      !Length of history header on file 1
     &,LEN_CFI11                                                        &
                      !Length of index1 on file 1
     &,LEN_CFI21                                                        &
                      !Length of index2 on file 1
     &,LEN_CFI31                                                        &
                      !Length of index3 on file 1
     &,LEN1_LOOKUP1                                                     &
                      !1st dim of LOOKUP on file 1
     &,LEN2_LOOKUP1                                                     &
                      !2nd dim of LOOKUP on file 1
     &,LEN_DATA1                                                        &
                      !Length of data on file 1
     &,ROW_LENGTH1                                                      &
                      !No of points E-W on file 1
     &,P_ROWS1                                                          &
                      !No of p-rows on file 1
     &,P_FIELD1                                                         &
                      !No of p-points per level on file 1
     &,MAX_FIELD_SIZE1 !Maximum field size on file 1

      INTEGER                                                           &
     & LEN_FIXHD2                                                       &
                      !Length of fixed length header on file 2
     &,LEN_INTHD2                                                       &
                      !Length of integer header on file 2
     &,LEN_REALHD2                                                      &
                      !Length of real header on file 2
     &,LEN1_LEVDEPC2                                                    &
                      !1st dim of lev dependent consts on file 2
     &,LEN2_LEVDEPC2                                                    &
                      !2nd dim of lev dependent consts on file 2
     &,LEN1_ROWDEPC2                                                    &
                      !1st dim of row dependent consts on file 2
     &,LEN2_ROWDEPC2                                                    &
                      !2nd dim of row dependent consts on file 2
     &,LEN1_COLDEPC2                                                    &
                      !1st dim of col dependent consts on file 2
     &,LEN2_COLDEPC2                                                    &
                      !2nd dim of col dependent consts on file 2
     &,LEN1_FLDDEPC2                                                    &
                      !1st dim of field dependent consts on file 2
     &,LEN2_FLDDEPC2                                                    &
                      !2nd dim of field dependent consts on file 2
     &,LEN_EXTCNST2                                                     &
                      !Length of extra consts on file 2
     &,LEN_DUMPHIST2                                                    &
                      !Length of history header on file 2
     &,LEN_CFI12                                                        &
                      !Length of index1 on file 2
     &,LEN_CFI22                                                        &
                      !Length of index2 on file 2
     &,LEN_CFI32                                                        &
                      !Length of index3 on file 2
     &,LEN1_LOOKUP2                                                     &
                      !1st dim of LOOKUP on file 2
     &,LEN2_LOOKUP2                                                     &
                      !2nd dim of LOOKUP on file 2
     &,LEN_DATA2                                                        &
                      !Length of data on file 2
     &,ROW_LENGTH2                                                      &
                      !No of points E-W on file 2
     &,P_ROWS2                                                          &
                      !No of p-rows on file 2
     &,P_FIELD2                                                         &
                      !No of p-points per level on file 2
     &,MAX_FIELD_SIZE2                                                  &
                      !Maximum field size on file 2
     &,num_lookup_ignore 
                    ! Number of items in lookup to ignore when calculating
                    ! number of differences

      INTEGER, ALLOCATABLE :: lookup_ignore(:)

      INTEGER                                                           &
     & LEN_IO                                                           &
                !Length of I/O returned by BUFFER IN
     &,I                                                                &
                !Loop index
     &,NFTIN1                                                           &
                !Unit number of input UM file 1
     &,NFTIN2                                                           &
                !Unit number of input UM file 2

     &,ERR                                                              &
                !Return code from OPEN
     &,ICODE    !Return code from setpos
      REAL A    !BUFFER IN UNIT function
!
      integer expand
      CHARACTER(LEN=100) DUMMY_ENV  ! replaces deprecated string which
                                    ! contained executable path
      INTEGER ME_GC,NPROC_GC

      CHARACTER (LEN=filenamelength) :: filename1
      CHARACTER (LEN=filenamelength) :: filename2

      CHARACTER(LEN=8) :: envname
                ! Environment variable name for ignored lookup items
      CHARACTER(LEN=2) :: envnumber
                ! Environment variable number for ignored lookup items
      CHARACTER(LEN=8) :: tmp_lookup_ignore
                ! Temporary variable to store lookup items from from environment

      INTEGER :: ig
                ! Loop counter
      LOGICAL :: ignore_missing_fields
      CHARACTER(LEN=8) :: c_num_lookup_ignore
      CHARACTER(LEN=1) :: c_ignore_missing

      DUMMY_ENV='dummy path'

      CALL GC_INIT(DUMMY_ENV,ME_GC,NPROC_GC)
! Initialise print status for standard output
! DEPENDS ON: initprintstatus
      CALL InitPrintStatus()
      CALL appInit(exe_cumf)
      CALL ioInit()

      expand=1

! Ascertain lookup items to ignore

! The number of items to ignore (required to allocate the array) is kept in
! the environment variable NUMIGNORE
      CALL fort_get_env('NUMIGNORE',9, c_num_lookup_ignore, 8, Err)

! If the environment variable is not set, use the default ignore list
      IF ( Err /=  0 ) THEN

!       Assign default items in lookup to ignore (none)
        num_lookup_ignore = 1
        ALLOCATE(lookup_ignore(num_lookup_ignore))
        lookup_ignore = (/ 0 /)
        
      ELSE

        READ( c_num_lookup_ignore, '(I8)') num_lookup_ignore

        ALLOCATE(lookup_ignore(num_lookup_ignore))

! Each item number to ignore (up to 99) are kept in consecutive environment
! variables IGNORE01, IGNORE02 etc. 
        DO ig = 1, num_lookup_ignore

          IF (ig < 10) THEN
            WRITE(envnumber, '("0",I1)') ig
          ELSE
            WRITE(envnumber, '(I2)') ig
          END IF

          WRITE(envname, '("IGNORE",A2)') envnumber

          CALL fort_get_env(envname,8, tmp_lookup_ignore, 8, Err)
          
          IF ( Err /= 0 ) THEN
            CALL EREPORT('COMPARE', Err,                        &
     &   'Error reading ignore environment variable')
          END IF

          READ( tmp_lookup_ignore, '(I8)') lookup_ignore(ig)
          
        END DO

      END IF

! Check to see if missing fields should be ignored or not
      CALL fort_get_env('IGNOREMISSING',13, c_ignore_missing, 1, Err)

! If the environment variable is not set, fail if missing fields present
      IF ( Err /=  0 ) THEN
        ignore_missing_fields = .FALSE.
      ELSE
        IF (c_ignore_missing == '1') THEN
          ignore_missing_fields = .TRUE.
        ELSE
          ignore_missing_fields = .FALSE.
        END IF
      END IF


!L 1. Assign unit numbers

      NFTIN1=20
      NFTIN2=21

      WRITE(6,*)' COMPARE - FULL MODE'
      WRITE(6,*)' -------------------'
      WRITE(6,*)' '

      CALL fort_get_env ('FILE1',5,filename1,filenamelength,err)
      CALL fort_get_env ('FILE2',5,filename2,filenamelength,err)

      write (6,*) ' Files being compared '
      write (6,*) ' -------------------- '
      write (6,*) ' File 1 : ',FileName1( 1:len_trim(FileName1) )
      write (6,*) ' File 2 : ',FileName2( 1:len_trim(FileName2) )
      write (6,*) ' '

      WRITE(6,'(20x,''FILE STATUS'')')
      WRITE(6,'(20x,''==========='')')
!     CALL OPEN(1,'PPXREF',6,0,0,ERR)
      CALL FILE_OPEN(NFTIN1,'FILE1',5,0,0,ERR)
      CALL FILE_OPEN(NFTIN2,'FILE2',5,0,0,ERR)

!L 2. Buffer in fixed length header record from file 1

      CALL BUFFIN(NFTIN1,FIXHD1,256,LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= 256)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of fixed length header of input file',  &
     &  A,LEN_IO,256)

      CALL EREPORT('MAIN_COMPARE', 1000,                                &
     & 'Buffer in of fixed length header wrong size')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,256
        IF(FIXHD1(I) <  0)FIXHD1(I)=0
      ENDDO

! Input file dimensions (ensure sizes are >= 0 for allocating)
      LEN_FIXHD1    = 256
      LEN_INTHD1    = MAX(FIXHD1(101),0)
      LEN_REALHD1   = MAX(FIXHD1(106),0)
      LEN1_LEVDEPC1 = MAX(FIXHD1(111),0)
      LEN2_LEVDEPC1 = MAX(FIXHD1(112),0)
      LEN1_ROWDEPC1 = MAX(FIXHD1(116),0)
      LEN2_ROWDEPC1 = MAX(FIXHD1(117),0)
      LEN1_COLDEPC1 = MAX(FIXHD1(121),0)
      LEN2_COLDEPC1 = MAX(FIXHD1(122),0)
      LEN1_FLDDEPC1 = MAX(FIXHD1(126),0)
      LEN2_FLDDEPC1 = MAX(FIXHD1(127),0)
      LEN_EXTCNST1  = MAX(FIXHD1(131),0)
      LEN_DUMPHIST1 = MAX(FIXHD1(136),0)
      LEN_CFI11     = MAX(FIXHD1(141),0)
      LEN_CFI21     = MAX(FIXHD1(143),0)
      LEN_CFI31     = MAX(FIXHD1(145),0)
      LEN1_LOOKUP1  = MAX(FIXHD1(151),0)
      LEN2_LOOKUP1  = MAX(FIXHD1(152),0)
      LEN_DATA1     = MAX(FIXHD1(161),0)

!L 3. Buffer in fixed length header record from file 2

      CALL BUFFIN(NFTIN2,FIXHD2,256,LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= 256)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of fixed length header of input file',  &
     &  A,LEN_IO,256)

      CALL EREPORT('MAIN_COMPARE', 1001,                                &
     & 'Buffer in of fixed length header wrong size')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,256
        IF(FIXHD2(I) <  0)FIXHD2(I)=0
      ENDDO

! Input file dimensions (ensure sizes are >= 0 for allocating)
      LEN_FIXHD2=256
      LEN_INTHD2    = MAX(FIXHD2(101),0)
      LEN_REALHD2   = MAX(FIXHD2(106),0)
      LEN1_LEVDEPC2 = MAX(FIXHD2(111),0)
      LEN2_LEVDEPC2 = MAX(FIXHD2(112),0)
      LEN1_ROWDEPC2 = MAX(FIXHD2(116),0)
      LEN2_ROWDEPC2 = MAX(FIXHD2(117),0)
      LEN1_COLDEPC2 = MAX(FIXHD2(121),0)
      LEN2_COLDEPC2 = MAX(FIXHD2(122),0)
      LEN1_FLDDEPC2 = MAX(FIXHD2(126),0)
      LEN2_FLDDEPC2 = MAX(FIXHD2(127),0)
      LEN_EXTCNST2  = MAX(FIXHD2(131),0)
      LEN_DUMPHIST2 = MAX(FIXHD2(136),0)
      LEN_CFI12     = MAX(FIXHD2(141),0)
      LEN_CFI22     = MAX(FIXHD2(143),0)
      LEN_CFI32     = MAX(FIXHD2(145),0)
      LEN1_LOOKUP2  = MAX(FIXHD2(151),0)
      LEN2_LOOKUP2  = MAX(FIXHD2(152),0)
      LEN_DATA2     = MAX(FIXHD2(161),0)


!L 4. Buffer in integer constants from file 1

       CALL BUFFIN(NFTIN1,INTHD1,LEN_INTHD1,LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= LEN_INTHD1)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of integer constants in input file 1',  &
     &  A,LEN_IO,LEN_INTHD1)

      CALL EREPORT('MAIN_COMPARE', 1002,                                &
     & 'Buffer in of integer constants wrong length')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,LEN_INTHD1
        IF(INTHD1(I) <  0)INTHD1(I)=0
      ENDDO

       ROW_LENGTH1=INTHD1(6)
       P_ROWS1=INTHD1(7)
       P_FIELD1=ROW_LENGTH1*P_ROWS1

!L Extract maximum field size from LOOKUP header
! DEPENDS ON: find_max_field_size
      CALL FIND_MAX_FIELD_SIZE(NFTIN1,LEN1_LOOKUP1,LEN2_LOOKUP1,FIXHD1    &
     &    ,max_field_size1, expand)

!L 5. Buffer in integer constants from file 2

       CALL BUFFIN(NFTIN2,INTHD2,LEN_INTHD2,LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= LEN_INTHD2)THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of integer constants in input file 2',  &
     &  A,LEN_IO,LEN_INTHD2)

      CALL EREPORT('MAIN_COMPARE', 1003,                                &
     & 'Buffer integer constants wrong length')

      ENDIF

! Set missing data indicator to zero
      DO  I=1,LEN_INTHD2
        IF(INTHD2(I) <  0)INTHD2(I)=0
      ENDDO

!L 6. Cause abort if files obviously different

      ROW_LENGTH2=INTHD2(6)
      P_ROWS2=INTHD2(7)
      P_FIELD2=ROW_LENGTH2*P_ROWS2

!L Extract maximum field size from LOOKUP header
! DEPENDS ON: find_max_field_size
      CALL FIND_MAX_FIELD_SIZE(NFTIN2,LEN1_LOOKUP2,LEN2_LOOKUP2,FIXHD2    &
     &    ,max_field_size2, expand)

      IF(P_FIELD1 /= P_FIELD2)THEN
       WRITE(6,*)'COMPARE: ERROR Dumps are at different resolutions'

       CALL EREPORT('MAIN_COMPARE', 1004,                               &
     &  'Dumps are at difference resolutions')

      ENDIF
      IF(LEN2_LOOKUP1 /= LEN2_LOOKUP2)THEN
       WRITE(6,*)                                                       &
     & 'COMPARE: WARNING Dumps have different number of fields'
      ENDIF

! Rewind files
      CALL SETPOS(NFTIN1,0,ICODE)
      CALL SETPOS(NFTIN2,0,ICODE)

!L 7. Call COMPARE

! DEPENDS ON: compare
      CALL COMPARE(LEN_FIXHD1,LEN_INTHD1,LEN_REALHD1,                   &
     &  LEN1_LEVDEPC1,LEN2_LEVDEPC1,LEN1_ROWDEPC1,                      &
     &  LEN2_ROWDEPC1,LEN1_COLDEPC1,LEN2_COLDEPC1,                      &
     &  LEN1_FLDDEPC1,LEN2_FLDDEPC1,LEN_EXTCNST1,                       &
     &  LEN_DUMPHIST1,LEN_CFI11,LEN_CFI21,LEN_CFI31,                    &
     &  LEN1_LOOKUP1,LEN2_LOOKUP1,LEN_DATA1,P_FIELD1,                   &
     &  P_ROWS1,P_ROWS2,ROW_LENGTH1,ROW_LENGTH2,                        &
     &  LEN_FIXHD2,LEN_INTHD2,LEN_REALHD2,                              &
     &  LEN1_LEVDEPC2,LEN2_LEVDEPC2,LEN1_ROWDEPC2,                      &
     &  LEN2_ROWDEPC2,LEN1_COLDEPC2,LEN2_COLDEPC2,                      &
     &  LEN1_FLDDEPC2,LEN2_FLDDEPC2,LEN_EXTCNST2,                       &
     &  LEN_DUMPHIST2,LEN_CFI12,LEN_CFI22,LEN_CFI32,                    &
     &  LEN1_LOOKUP2,LEN2_LOOKUP2,LEN_DATA2,P_FIELD2                    &
     & ,NFTIN1,NFTIN2,MAX_FIELD_SIZE1,MAX_FIELD_SIZE2,                  &
     & expand,num_lookup_ignore,lookup_ignore, ignore_missing_fields)

      DEALLOCATE(lookup_ignore)

      CALL ioShutdown()

      CALL ereport_finalise()

      CALL gc_exit()

      END PROGRAM MAIN_COMPARE
!*L  Arguments:-------------------------------------------------------

