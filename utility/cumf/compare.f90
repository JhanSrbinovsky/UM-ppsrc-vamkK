! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine COMPARE
!
!  Purpose: Compares two UM atmosphere, ocean, or ancillary files.
!           MAIN_COMPARE reads in fixed length and integer
!           headers of UM files to be compared, extracts dimensions
!           of each file and then passes these values to
!           subroutine COMPARE.
!
!          COMPARE subroutine:
!          Compares two UM atmosphere, ocean, or ancillary files.
!          COMPARE reads in headers and data fields from files on
!          NFTIN1 and NFTIN2, comparing values.
!          UNIT 6: If an exact compare is found the message 'OK'
!          is written out, otherwise
!          i)  if header, all differring values are printed
!          ii) if field, 1st 10 differring values are printed plus
!              the maximum difference between the fields.
!          iii) if field only present in one file, a warning message
!               is displayed
!          UNIT 7: Number of differences displayed for each header.
!                  Number of fields with differences is also
!                  displayed along with the number of differences
!                  for each field which has differences
!
!
!  Programming standard: UMDP 3
!
!  Logical components covered:
!
!  System Tasks: F3,F4,F6
!
!  Documentation: UM Doc Paper F5
!
!  -----------------------------------------------------------------
!    Arguments:-------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      SUBROUTINE COMPARE(LEN_FIXHD1,LEN_INTHD1,LEN_REALHD1,             &
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

      USE filenamelength_mod, ONLY :                                    & 
          filenamelength

      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Decomp_DB
      USE LBC_mod
      USE rimtypes
      USE lookup_addresses
      USE cppxref_mod, ONLY :                                          &
          ppxref_charlen, ppxref_codelen, ppx_grid_type,               &
          ppx_atm_lbc_theta, ppx_atm_lbc_orog, ppx_ocn_lbc_theta,      &
          ppx_atm_lbc_u, ppx_atm_lbc_v, ppx_ocn_lbc_u
      USE ppxlook_mod, ONLY : ppxrecs
! version_mod items required by cstash.h
      USE version_mod, ONLY :                                          &
          nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,        &
          nlevp, npslevp, npslistp, outfile_s, outfile_e

      USE Submodel_Mod
      IMPLICIT NONE

      INTEGER                                                           &
     & LEN_FIXHD1                                                       &
                    !IN Length of fixed length header on file 1
     &,LEN_INTHD1                                                       &
                    !IN Length of integer header on file 1
     &,LEN_REALHD1                                                      &
                    !IN Length of real header on file 1
     &,LEN1_LEVDEPC1                                                    &
                    !IN 1st dim of lev dependent consts on file 1
     &,LEN2_LEVDEPC1                                                    &
                    !IN 2nd dim of lev dependent consts on file 1
     &,LEN1_ROWDEPC1                                                    &
                    !IN 1st dim of row dependent consts on file 1
     &,LEN2_ROWDEPC1                                                    &
                    !IN 2nd dim of row dependent consts on file 1
     &,LEN1_COLDEPC1                                                    &
                    !IN 1st dim of col dependent consts on file 1
     &,LEN2_COLDEPC1                                                    &
                    !IN 2nd dim of col dependent consts on file 1
     &,LEN1_FLDDEPC1                                                    &
                    !IN 1st dim of field dependent consts on file 1
     &,LEN2_FLDDEPC1                                                    &
                    !IN 2nd dim of field dependent consts on file 1
     &,LEN_EXTCNST1                                                     &
                    !IN Length of extra consts on file 1
     &,LEN_DUMPHIST1                                                    &
                    !IN Length of history header on file 1
     &,LEN_CFI11                                                        &
                    !IN Length of index1 on file 1
     &,LEN_CFI21                                                        &
                    !IN Length of index2 on file 1
     &,LEN_CFI31                                                        &
                    !IN Length of index3 on file 1
     &,LEN1_LOOKUP1                                                     &
                    !IN 1st dim of LOOKUP on file 1
     &,LEN2_LOOKUP1                                                     &
                    !IN 2nd dim of LOOKUP on file 1
     &,LEN_DATA1                                                        &
                    !IN Length of data on file 1
     &,P_FIELD1                                                         &
                    !IN No of p-points per level on file 1
     &,P_ROWS1                                                          &
     &,ROW_LENGTH1                                                      &
     &,MAX_FIELD_SIZE1                                                  &
                       !IN Maximum field size on file 1
     &,expand ! IN set to 1 to expand WGDOS and run length encoded
              ! Fields for comparison.
!
      integer lblrec_1, lblrec_2, length_changed

      INTEGER                                                           &
     & LEN_FIXHD2                                                       &
                    !IN Length of fixed length header on file 2
     &,LEN_INTHD2                                                       &
                    !IN Length of integer header on file 2
     &,LEN_REALHD2                                                      &
                    !IN Length of real header on file 2
     &,LEN1_LEVDEPC2                                                    &
                    !IN 1st dim of lev dependent consts on file 2
     &,LEN2_LEVDEPC2                                                    &
                    !IN 2nd dim of lev dependent consts on file 2
     &,LEN1_ROWDEPC2                                                    &
                    !IN 1st dim of row dependent consts on file 2
     &,LEN2_ROWDEPC2                                                    &
                    !IN 2nd dim of row dependent consts on file 2
     &,LEN1_COLDEPC2                                                    &
                    !IN 1st dim of col dependent consts on file 2
     &,LEN2_COLDEPC2                                                    &
                    !IN 2nd dim of col dependent consts on file 2
     &,LEN1_FLDDEPC2                                                    &
                    !IN 1st dim of field dependent consts on file 2
     &,LEN2_FLDDEPC2                                                    &
                    !IN 2nd dim of field dependent consts on file 2
     &,LEN_EXTCNST2                                                     &
                    !IN Length of extra consts on file 2
     &,LEN_DUMPHIST2                                                    &
                    !IN Length of history header on file 2
     &,LEN_CFI12                                                        &
                    !IN Length of index1 on file 2
     &,LEN_CFI22                                                        &
                    !IN Length of index2 on file 2
     &,LEN_CFI32                                                        &
                    !IN Length of index3 on file 2
     &,LEN1_LOOKUP2                                                     &
                    !IN 1st dim of LOOKUP on file 2
     &,LEN2_LOOKUP2                                                     &
                    !IN 2nd dim of LOOKUP on file 2
     &,LEN_DATA2                                                        &
                    !IN Length of data on file 2
     &,P_FIELD2                                                         &
                    !IN No of p-points per level on file 2
     &,P_ROWS2                                                          &
     &,ROW_LENGTH2                                                      &
     &,MAX_FIELD_SIZE2 !IN Maximum field size on file 2

      INTEGER                                                           &
     & NFTIN1                                                           &
                    !IN Unit number for file 1
     &,NFTIN2       !IN Unit number for file 2

      INTEGER :: num_lookup_ignore
                    ! Number of items in lookup to ignore when calculating
                    ! number of differences
      INTEGER :: lookup_ignore(num_lookup_ignore)
                    ! Item numbers to ignore when calculating num. differences
      LOGICAL :: ignore_missing_fields
      

! Comdecks: ------------------------------------------------------------
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
!========================== COMDECK SX_SIZE ====================
!
!   Description:
!
!   This COMDECK contains COMMON blocks for landpt only and LBC
!   field for SX
!
! -------------------- COMMON BLOCKS --------------------
!
! =======================================================
! Common block for small executables
! =======================================================


! Variables and sizes for land-points only field
! ----------------------------------------------
      INTEGER:: LAND_FIELD           ! IN: No of land points in field
      INTEGER:: TR_VARS
      INTEGER:: TOT_LEVELS
      INTEGER:: global_ROWS          ! IN: No of global (theta) rows
      INTEGER:: global_ROW_LENGTH    ! IN: Points per global row
      INTEGER:: global_land_field    !  Global no of land points
      INTEGER:: local_land_field     !  Local no of land points
      COMMON /mpp_landpts/ global_land_field,local_land_field


! Variables and sizes for LBC
! ---------------------------

      ! Data structure sizes for OCEAN BOUNDARY file control routines
      INTEGER :: NRIM_TIMESO   ! IN: Max no of timelevels in rim flds

      ! Data structure sizes for WAVE BOUNDARY file control routines
      INTEGER :: RIMWIDTHW    ! IN: No of points width in rim fields
      INTEGER :: NRIM_TIMESW  ! IN: Max no of timelevels in rim flds

      ! Variables describing the Ocean Lateral Boundary Conditions
      INTEGER:: LENRIMO                ! Size of ocean LBC (theta)
      INTEGER:: LENRIMO_U              ! Size of ocean LBC (velocity)

      ! Variables that may be needed for vn5.2 but have not yet been
      ! dealt with at vn5.1
      INTEGER:: RIMFLDSO
      INTEGER:: global_LENRIMDATA_W
      INTEGER:: LENRIMDATA_O
      INTEGER:: LENRIMDATA_W
      INTEGER:: RIM_LOOKUPSO
      INTEGER:: RIM_LOOKUPSW
      INTEGER:: BOUND_LOOKUPSO
      INTEGER:: BOUND_LOOKUPSW

      COMMON/SX_NLSIZES/                                                &
     &  NRIM_TIMESO,                                                    &
     &  RIMWIDTHW, NRIM_TIMESW,                                         &
     &  LAND_FIELD, TR_VARS, TOT_LEVELS,                                &
     &  GLOBAL_ROW_LENGTH, GLOBAL_ROWS

      COMMON/DRSIZ_BO/                                                  &
      ! Atmosphere variables
      ! Wave model variables
     &  RIM_LOOKUPSW, LENRIMDATA_W, global_LENRIMDATA_W,                &
      ! Ocean variables
     &  LENRIMO, LENRIMO_U,                                             &
      ! Variables still to be dealt with
     &  RIMFLDSO,RIM_LOOKUPSO,          &
     &  BOUND_LOOKUPSO,BOUND_LOOKUPSW,                   &
     &  LENRIMDATA_O
! SX_SIZE end
!====================== COMDECK ATM_LSM ========================
! Description:
!   This comdeck contains a COMMON block which contains the
!   atmosphere land sea mask - both the full field, and the
!   local subdomain on this processor.
!   This data is required for various compression/decompression
!   algorithms.
!
!   Requires AMAXSIZE comdeck to be called first for Max2DFieldSize
!

      LOGICAL                                                           &
!  Full-grid land-sea mask:
     &  atmos_landmask(Max2DFieldSize)                                  &
! Local subdomain area land-sea mask:
     &, atmos_landmask_local(Max2DFieldSize)

      INTEGER atmos_number_of_landpts ! total number of land points

      COMMON /Atmos_LSM_Common/                                         &
     &  atmos_landmask                                                  &
     &, atmos_landmask_local                                            &
     &, atmos_number_of_landpts

! End of comdeck ATM_LSM

! Local arrays:---------------------------------------------------------
      INTEGER                                                           &
     & FIXHD1(LEN_FIXHD1),                                              &
                                                 !
     & INTHD1(LEN_INTHD1),                                              &
                                                 !\                    .
     & CFI11(LEN_CFI11+1),CFI21(LEN_CFI21+1),                           &
                                                 ! > file 1 headers
     & CFI31(LEN_CFI31+1),                                              &
                                                 !/
     & LOOKUP1(LEN1_LOOKUP1,LEN2_LOOKUP1)        !

      INTEGER                                                           &
     & FIXHD2(LEN_FIXHD2),                                              &
                                                 !
     & INTHD2(LEN_INTHD2),                                              &
                                                 !\                    .
     & CFI12(LEN_CFI12+1),CFI22(LEN_CFI22+1),                           &
                                                 ! > file 2 headers
     & CFI32(LEN_CFI32+1),                                              &
                                                 !/
     & LOOKUP2(LEN1_LOOKUP2,LEN2_LOOKUP2)        !

      REAL                                                              &
     & REALHD1(LEN_REALHD1),                                            &
                                                 !
     & LEVDEPC1(1+LEN1_LEVDEPC1*LEN2_LEVDEPC1),                         &
                                                 !
     & ROWDEPC1(1+LEN1_ROWDEPC1*LEN2_ROWDEPC1),                         &
                                                 !\                    .
     & COLDEPC1(1+LEN1_COLDEPC1*LEN2_COLDEPC1),                         &
                                                 ! > file 1 headers
     & FLDDEPC1(1+LEN1_FLDDEPC1*LEN2_FLDDEPC1),                         &
                                                 !/
     & EXTCNST1(LEN_EXTCNST1+1),                                        &
                                                 !
     & DUMPHIST1(LEN_DUMPHIST1+1),                                      &
                                                 !
     & R_D1(MAX_FIELD_SIZE1) ! REAL Array for field on file 1

      INTEGER                                                           &
     & I_D1(MAX_FIELD_SIZE1) ! INTEGER Array for field on file 1

      LOGICAL                                                           &
     & L_D1(MAX_FIELD_SIZE1) ! LOGICAL Array for field on file 1

      REAL                                                              &
     & REALHD2(LEN_REALHD2),                                            &
                                                 !
     & LEVDEPC2(1+LEN1_LEVDEPC2*LEN2_LEVDEPC2),                         &
                                                 !
     & ROWDEPC2(1+LEN1_ROWDEPC2*LEN2_ROWDEPC2),                         &
                                                 !\                    .
     & COLDEPC2(1+LEN1_COLDEPC2*LEN2_COLDEPC2),                         &
                                                 ! > file 2 headers
     & FLDDEPC2(1+LEN1_FLDDEPC2*LEN2_FLDDEPC2),                         &
                                                 !/
     & EXTCNST2(LEN_EXTCNST2+1),                                        &
                                                 !
     & DUMPHIST2(LEN_DUMPHIST2+1),                                      &
                                                 !
     & R_D2(MAX_FIELD_SIZE2) ! REAL Array for field on file 2

      INTEGER                                                           &
     & I_D2(MAX_FIELD_SIZE2) ! INTEGER Array for field on file 2

      LOGICAL                                                           &
     & L_D2(MAX_FIELD_SIZE2) ! LOGICAL Array for field on file 2

      INTEGER                                                           &
     & PP_XREF(PPXREF_CODELEN)  !PPXREF codes for a given section/item

      LOGICAL                                                           &
     &  LAND_MASK_FOUND  ! Is there a land mask in the dump

!*----------------------------------------------------------------------
!*L  Local variables:---------------------------------------------------
      REAL                                                              &
     & MAX_DIFF                                                         &
                 ! Maximum difference between two real fields
     &,RD1,RD2   ! Real variables to equivalent with LD1/ID1 & LD2/ID2
      REAL DIFF_PER,RMS_F1,RMS_F2,RMS_DIFF,A,RCODE

      INTEGER                                                           &
     & ICODE                                                            &
                    ! Error return code from subroutines
     &,START_BLOCK                                                      &
                    ! READHEAD argument (not used)
     &,I,J,K,N,ig                                                       &
                    ! Loop indices
     &,L,M,n_ignored                                                    &
                    ! Difference counters
     &,JMIN                                                             &
                    ! Minimum length of two headers
     &,S_ITEM_CODE                                                      &
                    ! STASH item code
     &,SECTION                                                          &
                    ! STASH section number
     &,ID1,ID2                                                          &
                    ! Integer variables to equivalent with RD1 and RD2
     &,N_DIFF                                                           &
                    ! No of differences to be listed
     &,IMAX_DIFF                                                        &
                    ! Maximum difference between two integer fields
     &,PACK_CODE1                                                       &
                    ! Packing code for LOOKUP table 1
     &,PACK_CODE2                                                       &
                    ! Packing code for LOOKUP table 2
     &,MAX_J        ! Location of max.diff


      INTEGER RowNumber
      INTEGER MODEL              !Internal model number
      INTEGER LEN_FIELD          !Number of points in field to be
                                 !compared
      INTEGER N1,N2
      INTEGER OFFSET1
      INTEGER OFFSET2
      INTEGER NUMREC1
      INTEGER NUMREC2
      INTEGER NMISSING1
      INTEGER NMISSING2
      INTEGER IROWDEPC1
      INTEGER IROWDEPC2
      INTEGER IMASK,IMASK_DUMMY,LEN_IO
      INTEGER HALOX1,HALOY1,HALOX2,HALOY2
      INTEGER LAND_POINTS
      INTEGER OLEN_FIELD
      INTEGER FIELD_ITEM1,FIELD_SECT1,FIELD_MODEL1
      INTEGER FIELD_ITEM2,FIELD_SECT2,FIELD_MODEL2
      INTEGER GRID_TYPE1,GRID_TYPE2
      INTEGER INDEX(LEN2_LOOKUP1)
      INTEGER NDIFFER(LEN2_LOOKUP1)
      LOGICAL LMISSING1(LEN2_LOOKUP1)
      LOGICAL LMISSING2(LEN2_LOOKUP2)
       INTEGER      EXPPXI
       CHARACTER(LEN=36) EXPPXC

      LOGICAL                                                           &
     & LD1,LD2      ! Logical variables to equivalent with RD1 and RD2

      CHARACTER                                                         &
     & CMESSAGE*100                                                     &
                    ! Character string returned if ICODE  /=  0
     &,PHRASE*(PPXREF_CHARLEN) ! Name of field
      CHARACTER(LEN=1) DIFF(MAX_FIELD_SIZE1)
      CHARACTER(LEN=200) KEY
      CHARACTER(LEN=filenamelength) :: filename ! Name of output file

      EQUIVALENCE (RD1,ID1,LD1) , (RD2,ID2,LD2)

      PARAMETER (N_DIFF=10)
      INTEGER NFT1,NFT2
      PARAMETER (NFT1=22, NFT2=2)

!*----------------------------------------------------------------------


!  0. Open PPXREF file

      ppxRecs=1
      RowNumber=0
      cmessage = ' '
      ICODE=0
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_A',ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_A'
        WRITE(6,*) CMESSAGE

        CALL EREPORT('COMPARE', ICODE,                                  &
     &   'Error reading STASHmaster_A')
      END IF

! DEPENDS ON: hdppxrf
      CALL HDPPXRF(NFT1,'STASHmaster_O',ICODE,CMESSAGE)
      IF(ICODE >  0)THEN
        WRITE(6,*) 'Error reading STASHmaster_O'
        WRITE(6,*) CMESSAGE

        CALL EREPORT('COMPARE', ICODE,                                  &
     &   'Error reading STASHmaster_O')
      END IF

! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_A',RowNumber,                  &
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE

        CALL EREPORT('COMPARE', ICODE,                                  &
     &   CMESSAGE)
      ENDIF

! DEPENDS ON: getppx
      CALL GETPPX(NFT1,NFT2,'STASHmaster_O',RowNumber,                  &
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE

        CALL EREPORT('COMPARE', ICODE,                                  &
     &   CMESSAGE)
      ENDIF

!User STASHmaster
! DEPENDS ON: hdppxrf
      CALL HDPPXRF(0,' ',ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE

        CALL EREPORT('COMPARE', ICODE,                                  &
     &   CMESSAGE)
      ENDIF

! DEPENDS ON: getppx
      CALL GETPPX(0,NFT2,' ',RowNumber,                                 &
     &            ICODE,CMESSAGE)
      IF(ICODE /= 0)THEN
        WRITE(6,*) CMESSAGE

        CALL EREPORT('COMPARE', ICODE,                                  &
     &   CMESSAGE)
      ENDIF

! 1: Open output files
!
! Open up unit 7: Summary part one
      CALL GET_FILE(7,FILENAME,filenamelength,ICODE)
      OPEN(7,FILE=FILENAME,STATUS='NEW',IOSTAT=ICODE)
      IF (ICODE /= 0) THEN
        WRITE(6,*) 'Can not write to ',FILENAME
        WRITE(CMESSAGE,'(A)') 'Error opening file on unit 7'
        CALL ereport('COMPARE', icode, cmessage)
      ELSE
        WRITE(6,*) 'OPEN: 7:',TRIM(FILENAME),'has been created'
      ENDIF
      WRITE(7,*)' COMPARE - SUMMARY MODE'
      WRITE(7,*)'-----------------------'
      WRITE(7,*)' '

! Open up unit 8: Summary part two
      CALL GET_FILE(8,FILENAME,filenamelength,ICODE)
      OPEN(8,FILE=FILENAME,STATUS='NEW',IOSTAT=ICODE)
      IF (ICODE /= 0) THEN
        WRITE(6,*) 'Can not write to ',FILENAME
        WRITE(CMESSAGE,'(A)') 'Error opening file on unit 8'
        CALL ereport('COMPARE', icode, cmessage)
      ELSE
        WRITE(6,*) 'OPEN: 8:',TRIM(FILENAME),'has been created'
      ENDIF

 ! Open up unit 10
      CALL GET_FILE(10,FILENAME,filenamelength,ICODE)
      OPEN(10,FILE=FILENAME,STATUS='NEW',IOSTAT=ICODE)
      IF (ICODE /= 0) THEN
        WRITE(6,*) 'Can not write to ',FILENAME
        WRITE(CMESSAGE,'(A)') 'Error opening file on unit 10'
        CALL ereport('COMPARE', icode, cmessage)
      ELSE
        WRITE(6,*) 'OPEN: 10:',TRIM(FILENAME),'has been created'
      ENDIF
      WRITE(10,*)' COMPARE - DIFFERENCE CHARTS'
      WRITE(10,*)'----------------------------'
      WRITE(10,*)' '

      WRITE(6,*)' '
      WRITE(6,*)'          FILE 1'
      WRITE(6,*)'          ------'
! DEPENDS ON: readhead
      CALL READHEAD(NFTIN1,FIXHD1,LEN_FIXHD1,                           &
     &                INTHD1,LEN_INTHD1,                                &
     &                REALHD1,LEN_REALHD1,                              &
     &                LEVDEPC1,LEN1_LEVDEPC1,LEN2_LEVDEPC1,             &
     &                ROWDEPC1,LEN1_ROWDEPC1,LEN2_ROWDEPC1,             &
     &                COLDEPC1,LEN1_COLDEPC1,LEN2_COLDEPC1,             &
     &                FLDDEPC1,LEN1_FLDDEPC1,LEN2_FLDDEPC1,             &
     &                EXTCNST1,LEN_EXTCNST1,                            &
     &                DUMPHIST1,LEN_DUMPHIST1,                          &
     &                CFI11,LEN_CFI11,                                  &
     &                CFI21,LEN_CFI21,                                  &
     &                CFI31,LEN_CFI31,                                  &
     &                LOOKUP1,LEN1_LOOKUP1,LEN2_LOOKUP1,                &
     &                LEN_DATA1,                                        &
     &                START_BLOCK,ICODE,CMESSAGE)

      IF(ICODE /= 0)THEN
        WRITE(6,*)CMESSAGE,ICODE

        CALL EREPORT('COMPARE', ICODE,                                  &
     &   CMESSAGE)

      ENDIF

!L 2. Read in file 2 header

      WRITE(6,*)' '
      WRITE(6,*)'          FILE 2'
      WRITE(6,*)'          ------'
! DEPENDS ON: readhead
      CALL READHEAD(NFTIN2,FIXHD2,LEN_FIXHD2,                           &
     &                INTHD2,LEN_INTHD2,                                &
     &                REALHD2,LEN_REALHD2,                              &
     &                LEVDEPC2,LEN1_LEVDEPC2,LEN2_LEVDEPC2,             &
     &                ROWDEPC2,LEN1_ROWDEPC2,LEN2_ROWDEPC2,             &
     &                COLDEPC2,LEN1_COLDEPC2,LEN2_COLDEPC2,             &
     &                FLDDEPC2,LEN1_FLDDEPC2,LEN2_FLDDEPC2,             &
     &                EXTCNST2,LEN_EXTCNST2,                            &
     &                DUMPHIST2,LEN_DUMPHIST2,                          &
     &                CFI12,LEN_CFI12,                                  &
     &                CFI22,LEN_CFI22,                                  &
     &                CFI32,LEN_CFI32,                                  &
     &                LOOKUP2,LEN1_LOOKUP2,LEN2_LOOKUP2,                &
     &                LEN_DATA2,                                        &
     &                START_BLOCK,ICODE,CMESSAGE)


      IF(ICODE /= 0)THEN
        WRITE(6,*)CMESSAGE,ICODE

        CALL EREPORT('COMPARE', ICODE,                                  &
     &   CMESSAGE)

      ENDIF

!L 3. Compare fixed length headers

      IF(FIXHD1(5) /= FIXHD2(5))THEN
        WRITE(6,'(''WARNING: FIXHD1(5)  = '',I3,'' FIXHD2(5)  = '',I3)')&
     &    FIXHD1(5),FIXHD2(5)
        WRITE(7,'(''WARNING: FIXHD1(5)  = '',I3,'' FIXHD2(5)  = '',I3)')&
     &    FIXHD1(5),FIXHD2(5)
        WRITE(6,'(''         File types are different'')')
        WRITE(7,'(''         File types are different'')')
      END IF
      WRITE(6,*)' '
      WRITE(6,*)'FIXED LENGTH HEADER:'

! Check length of fixed length headers
      JMIN=MIN0(LEN_FIXHD1,LEN_FIXHD2)
      IF(LEN_FIXHD1 /= LEN_FIXHD2)THEN
        WRITE(6,'(''WARNING: LEN_FIXHD1 = '',I3,'' LEN_FIXHD2 = '',I3)')&
     &    LEN_FIXHD1,LEN_FIXHD2
        WRITE(7,'(''WARNING: LEN_FIXHD1 = '',I3,'' LEN_FIXHD2 = '',I3)')&
     &    LEN_FIXHD1,LEN_FIXHD2
        WRITE(6,'(A)')'         Fixed length headers have different '//&
            'lengths'
        WRITE(7,'(A)')'         Fixed length headers have different '//    &
            'lengths'
        WRITE(6,'(''         Comparing first '',I3,''elements only'')') &
     &    JMIN
        WRITE(7,'(''         Comparing first '',I3,''elements only'')') &
     &    JMIN
      END IF

! Check fixed length header
      IF(FIXHD1(152) == FIXHD2(152))THEN
        IF(FIXHD1(160) /= FIXHD2(160))THEN
          WRITE(6,'(''WARNING: LEN1 = '',i9,'' and LEN2 = '',i9)')      &
     &      FIXHD1(160),FIXHD2(160)
          WRITE(7,'(''WARNING: LEN1 = '',i9,'' and LEN2 = '',i9)')      &
     &      FIXHD1(160),FIXHD2(160)
          WRITE(6,'(''         Data start address differs'')')
          WRITE(7,'(''         Data start address differs'')')
          WRITE(6,'(A)')'         Possibly due to comparing old and new '//&
              'format UM dumps or fieldsfiles'
          WRITE(7,'(A)')'         Possibly due to comparing old and new '//&
              'format UM dumps or fieldsfiles'
        ELSE IF(FIXHD1(161) /= FIXHD2(161))THEN
          WRITE(6,'(''WARNING: LEN1 = '',i9,'' and LEN2 = '',i9)')      &
     &      FIXHD1(161),FIXHD2(161)
          WRITE(7,'(''WARNING: LEN1 = '',i9,'' and LEN2 = '',i9)')      &
     &      FIXHD1(161),FIXHD2(161)
          WRITE(6,'(''         Length of data differs'')')
          WRITE(7,'(''         Length of data differs'')')
          WRITE(6,'(''         Possibly due to comparing old and new '',&
     &              ''format UM dumps or fieldsfiles'')')
          WRITE(7,'(''         Possibly due to comparing old and new '',&
     &              ''format UM dumps or fieldsfiles'')')
        END IF
      END IF

      K = 0
      length_changed=0
      DO I=1,JMIN
        IF(FIXHD1(I) /= FIXHD2(I))THEN
        WRITE(6,'(''ITEM = '',i4,''  Values = '',i11,'' and '',i11)')   &
     &    i, fixhd1(i), fixhd2(i)
        K = K + 1
        ENDIF
      ENDDO

      IF(K == 0) WRITE(6,*) 'OK'
      WRITE(8,'(a,a,i7)') 'FIXED LENGTH HEADER:        ',               &
     &           'Number of differences = ',K

!L 4. Compare integer headers

      IF(LEN_INTHD1 >  0.OR.LEN_INTHD2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'INTEGER HEADER:'
        IF(LEN_INTHD1 /= LEN_INTHD2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_INTHD1,' LEN2=',LEN_INTHD2
        ENDIF
        JMIN=MIN0(LEN_INTHD1,LEN_INTHD2)
        K=0
        DO I=1,JMIN
          IF(INTHD1(I) /= INTHD2(I))THEN
            K=K+1
            WRITE(6,*)'ITEM=',I,INTHD1(I),INTHD2(I)
          ENDIF
        ENDDO
      ENDIF

      IF(K == 0) WRITE(6,*) 'OK'
      WRITE(8,'(a,a,i7)') 'INTEGER HEADER:             ',               &
     &           'Number of differences = ',K
      L=K

!L 5. Compare real headers

      IF(LEN_REALHD1 >  0.OR.LEN_REALHD2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'REAL HEADER:'
        IF(LEN_REALHD1 /= LEN_REALHD2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_REALHD1,' LEN2=',LEN_REALHD2
        ENDIF
        JMIN=MIN0(LEN_REALHD1,LEN_REALHD2)
        K=0
        DO I=1,JMIN
          IF(REALHD1(I) /= REALHD2(I))THEN
            K=K+1
            WRITE(6,*)'ITEM=',I,REALHD1(I),REALHD2(I)
          ENDIF
        ENDDO
      ENDIF

      IF(K == 0) WRITE(6,*) 'OK'
      WRITE(8,'(a,a,i7)') 'REAL HEADER:                ',               &
     &           'Number of differences = ',K
      L=L+K

!L 6. Compare level dependent constants

      WRITE(6,*)' '
      WRITE(6,*)'LEVEL DEPENDENT CONSTS:'
      IF(FIXHD1(110) >  0 .AND. FIXHD2(110) >  0) THEN
      IF(LEN1_LEVDEPC1 /= LEN1_LEVDEPC2)THEN
        WRITE(6,*)'ERROR : different number of levels'
        WRITE(6,*)'LEV1=',LEN1_LEVDEPC1,' LEV2=',LEN1_LEVDEPC2

        CALL EREPORT('COMPARE', 1009,                                   &
     &   'Different number of levels')

      ELSEIF(LEN2_LEVDEPC1 >  0.OR.LEN2_LEVDEPC2 >  0)THEN
        IF(LEN2_LEVDEPC1 /= LEN2_LEVDEPC2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN2_LEVDEPC1,' LEN2=',LEN2_LEVDEPC2
        ENDIF
        JMIN=MIN0(LEN2_LEVDEPC1,LEN2_LEVDEPC2)
        K=0
        DO I=1,JMIN
          DO J=1,LEN1_LEVDEPC1
            IF(LEVDEPC1((I-1)*LEN1_LEVDEPC1+J) /=                       &
     &        LEVDEPC2((I-1)*LEN1_LEVDEPC1+J))THEN
              K=K+1
              WRITE(6,*)'LEVEL=',J,'ITEM=',I,                           &
     &        LEVDEPC1((I-1)*LEN1_LEVDEPC1+J),                          &
     &        LEVDEPC2((I-1)*LEN1_LEVDEPC1+J)
           ENDIF
          ENDDO
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,'(a,a,i7)') 'LEVEL DEPENDENT CONSTANTS:  ',             &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(110) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(110) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 7. Compare row dependent constants

      WRITE(6,*)' '
      WRITE(6,*)'ROW DEPENDENT CONSTS:'
      IF(FIXHD1(115) >  0 .AND. FIXHD2(115) >  0) THEN
      IF(LEN1_ROWDEPC1 /= LEN1_ROWDEPC2)THEN
        WRITE(6,*)'ERROR : different number of rows'
        WRITE(6,*)'ROW1=',LEN1_ROWDEPC1,' ROW2=',LEN1_ROWDEPC2

        CALL EREPORT('COMPARE', 1010,                                   &
     &   'Different number of rows')

      ELSEIF(LEN2_ROWDEPC1 >  0.OR.LEN2_ROWDEPC2 >  0)THEN
        IF(LEN2_ROWDEPC1 /= LEN2_ROWDEPC2)THEN
          WRITE(6,*)'WARNING different second dimension'
          WRITE(6,*)'LEN1=',LEN2_ROWDEPC1,' LEN2=',LEN2_ROWDEPC2
        ENDIF
        JMIN=MIN0(LEN2_ROWDEPC1,LEN2_ROWDEPC2)
        K=0
        DO I=1,JMIN
          DO J=1,LEN1_ROWDEPC1
            IF(ROWDEPC1((I-1)*LEN1_ROWDEPC1+J) /=                       &
     &         ROWDEPC2((I-1)*LEN1_ROWDEPC1+J))THEN
              K=K+1
              WRITE(6,*)'ROW=',I,'ITEM=',J,                             &
     &        ROWDEPC1((I-1)*LEN1_ROWDEPC1+J),                          &
     &        ROWDEPC2((I-1)*LEN1_ROWDEPC1+J)
            ENDIF
          ENDDO
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,'(a,a,i7)') 'ROW DEPENDENT CONSTANTS:    ',             &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(115) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(115) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 8. Compare column dependent constants

      WRITE(6,*)' '
      WRITE(6,*)'COLUMN DEPENDENT CONSTS:'
      IF(FIXHD1(120) >  0 .AND. FIXHD2(120) >  0) THEN
      IF(LEN1_COLDEPC1 /= LEN1_COLDEPC2)THEN
        WRITE(6,*)'ERROR : different number of columns.'
        WRITE(6,*)'COL1=',LEN1_COLDEPC1,' COL2=',LEN1_COLDEPC2

        CALL EREPORT('COMPARE', 1011,                                   &
     &   'Different number of columns')

      ELSEIF(LEN2_COLDEPC1 >  0.OR.LEN2_COLDEPC2 >  0)THEN
        IF(LEN2_COLDEPC1 /= LEN2_COLDEPC2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN2_COLDEPC1,' LEN2=',LEN2_COLDEPC2
        ENDIF
        JMIN=MIN0(LEN2_COLDEPC1,LEN2_COLDEPC2)
        K=0
        DO I=1,JMIN
          DO J=1,LEN1_COLDEPC1
            IF(COLDEPC1((I-1)*LEN1_COLDEPC1+J) /=                       &
     &        COLDEPC2((I-1)*LEN1_COLDEPC1+J))THEN
              K=K+1
              WRITE(6,*)'COL=',I,'ITEM=',J,                             &
     &        COLDEPC1((I-1)*LEN1_COLDEPC1+J),                          &
     &        COLDEPC2((I-1)*LEN1_COLDEPC1+J)
            ENDIF
          ENDDO
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,'(a,a,i7)') 'COLUMN DEPENDENT CONSTANTS: ',             &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(120) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(120) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 9. Compare field dependent constants

      WRITE(6,*)' '
      WRITE(6,*)'FIELD DEPENDENT CONSTS:'
      IF(FIXHD1(125) >  0 .AND. FIXHD2(125) >  0) THEN
      IF(LEN1_FLDDEPC1 /= LEN1_FLDDEPC2)THEN
        WRITE(6,*)'ERROR : different number of fields.'
        WRITE(6,*)'FLD1=',LEN1_FLDDEPC1,' FLD2=',LEN1_FLDDEPC2

        CALL EREPORT('COMPARE', 1012,                                   &
     &   'Different number of fields')

      ELSEIF(LEN2_FLDDEPC1 >  0.OR.LEN2_FLDDEPC2 >  0)THEN
        IF(LEN2_FLDDEPC1 /= LEN2_FLDDEPC2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN2_FLDDEPC1,' LEN2=',LEN2_FLDDEPC2
        ENDIF
        JMIN=MIN0(LEN2_FLDDEPC1,LEN2_FLDDEPC2)
        K=0
        DO I=1,JMIN
          DO J=1,LEN1_FLDDEPC1
            IF(FLDDEPC1((I-1)*LEN1_FLDDEPC1+J) /=                       &
     &        FLDDEPC2((I-1)*LEN1_FLDDEPC1+J))THEN
             K=K+1
             WRITE(6,*)'FIELD=',J,'ITEM=',I,                            &
     &       FLDDEPC1((I-1)*LEN1_FLDDEPC1+J),                           &
     &       FLDDEPC2((I-1)*LEN1_FLDDEPC1+J)
           ENDIF
          ENDDO
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,'(a,a,i7)') 'FIELD DEPENDENT CONSTANTS:  ',             &
     &             'Number of differences = ',K

        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(125) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(125) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 10. Compare extra constants

      WRITE(6,*)' '
      WRITE(6,*)'EXTRA CONSTANTS:'
      IF(FIXHD1(130) >  0 .AND. FIXHD2(130) >  0) THEN
      IF(LEN_EXTCNST1 >  0.OR.LEN_EXTCNST2 >  0)THEN
        IF(LEN_EXTCNST1 /= LEN_EXTCNST2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_EXTCNST1,' LEN2=',LEN_EXTCNST2
        ENDIF
        JMIN=MIN0(LEN_EXTCNST1,LEN_EXTCNST2)
        K=0
        DO I=1,JMIN
          IF(EXTCNST1(I) /= EXTCNST2(I))THEN
            K=K+1
            WRITE(6,*)'ITEM=',I,EXTCNST1(I),EXTCNST2(I)
          ENDIF
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,'(a,a,i7)') 'EXTRA CONSTANTS:            ',             &
     &             'Number of differences = ',K

        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(130) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(130) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 11. Compare dump history

      WRITE(6,*)' '
      WRITE(6,*)'HISTORY BLOCK:'
      IF(FIXHD1(135) >  0 .AND. FIXHD2(135) >  0) THEN
      IF(LEN_DUMPHIST1 >  0.OR.LEN_DUMPHIST2 >  0)THEN
        IF(LEN_DUMPHIST1 /= LEN_DUMPHIST2)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_DUMPHIST1,' LEN2=',LEN_DUMPHIST2
        ENDIF
        JMIN=MIN0(LEN_DUMPHIST1,LEN_DUMPHIST2)
        K=0
        DO I=1,JMIN
          IF(DUMPHIST1(I) /= DUMPHIST2(I))THEN
            K=K+1
            WRITE(6,*)'ITEM=',I,DUMPHIST1(I),DUMPHIST2(I)
          ENDIF
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,'(a,a,i7)') 'HISTORY BLOCK:              ',             &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(135) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(135) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 12. Compare compressed index 1

      WRITE(6,*)' '
      WRITE(6,*)'COMPRESSED INDEX 1:'
      IF(FIXHD1(140) >  0 .AND. FIXHD2(140) >  0) THEN
      IF(LEN_CFI11 >  0.OR.LEN_CFI12 >  0)THEN
        IF(LEN_CFI11 /= LEN_CFI12)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_CFI11,' LEN2=',LEN_CFI12
        ENDIF
        JMIN=MIN0(LEN_CFI11,LEN_CFI12)
        K=0
        DO I=1,JMIN
          IF(CFI11(I) /= CFI12(I))THEN
            K=K+1
            WRITE(6,*)'ITEM=',I,CFI11(I),CFI12(I)
          ENDIF
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,'(a,a,i7)') 'COMPRESSED INDEX 1:         ',             &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(140) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(140) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 13. Compare compressed index 2

      WRITE(6,*)' '
      WRITE(6,*)'COMPRESSED INDEX 2:'
      IF(FIXHD1(142) >  0 .AND. FIXHD2(142) >  0) THEN
      IF(LEN_CFI21 >  0.OR.LEN_CFI22 >  0)THEN
        IF(LEN_CFI21 /= LEN_CFI22)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_CFI21,' LEN2=',LEN_CFI22
        ENDIF
        JMIN=MIN0(LEN_CFI21,LEN_CFI22)
        K=0
        DO I=1,JMIN
          IF(CFI21(I) /= CFI22(I))THEN
            K=K+1
            WRITE(6,*)'ITEM=',I,CFI21(I),CFI22(I)
          ENDIF
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,'(a,a,i7)') 'COMPRESSED INDEX 2:         ',             &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(142) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(142) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 14. Compare compressed index 3

      WRITE(6,*)' '
      WRITE(6,*)'COMPRESSED INDEX 3:'
      IF(FIXHD1(144) >  0 .AND. FIXHD2(144) >  0) THEN
      IF(LEN_CFI31 >  0.OR.LEN_CFI32 >  0)THEN
        IF(LEN_CFI31 /= LEN_CFI32)THEN
          WRITE(6,*)'WARNING LEN1=',LEN_CFI31,' LEN2=',LEN_CFI32
        ENDIF
        JMIN=MIN0(LEN_CFI31,LEN_CFI32)
        K=0
        DO I=1,JMIN
          IF(CFI31(I) /= CFI32(I))THEN
            K=K+1
            WRITE(6,*)'ITEM=',I,CFI31(I),CFI32(I)
          ENDIF
        ENDDO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(8,'(a,a,i7)') 'COMPRESSED INDEX 3:         ',             &
     &             'Number of differences = ',K
        L=L+K
      ENDIF
      ELSE
        WRITE(6,*)'No comparison done'
        IF (FIXHD1(144) <= 0) WRITE(6,*)'No array in FILE 1'
        IF (FIXHD2(144) <= 0) WRITE(6,*)'No array in FILE 2'
      ENDIF

!L 15. Compare lookup tables

      IF(LEN1_LOOKUP1 /= LEN1_LOOKUP2)THEN
        WRITE(6,*)'ERROR first dimensions of lookup tables different'
        WRITE(6,*)'LEN1=',LEN1_LOOKUP1,' LEN2=',LEN1_LOOKUP2

        CALL EREPORT('COMPARE', 1013,                                   &
     &   'First dimensions of lookup tables different')

      ENDIF
      IF(LEN2_LOOKUP1 >  0.OR.LEN2_LOOKUP2 >  0)THEN
        WRITE(6,*)' '
        WRITE(6,*)'LOOKUP:'

! Check length of lookup tables
        IF(FIXHD1(5) == 3)THEN
          DO I=1,LEN2_LOOKUP1
            IF(LOOKUP1(1,I) /= -99)NUMREC1=I
          END DO
          DO I=1,LEN2_LOOKUP2
            IF(LOOKUP2(1,I) /= -99)NUMREC2=I
          END DO
        ELSE
          NUMREC1=LEN2_LOOKUP1
          NUMREC2=LEN2_LOOKUP2
        END IF

        IF(LEN2_LOOKUP1 /= LEN2_LOOKUP2)THEN
           WRITE(6,'(''WARNING LEN1 = '',i9,'' and LEN2 = '',i9)')      &
     &       len2_lookup1, len2_lookup2
          IF(FIXHD1(5) == 3)THEN
            WRITE(6,'(''Fieldsfile file1 contains '',i5,'' fields '',   &
     &                ''and '',I5,'' empty records'')')                 &
     &               NUMREC1,LEN2_LOOKUP1-NUMREC1
            WRITE(7,'(''Fieldsfile file1 contains '',i5,'' fields '',   &
     &                ''and '',I5,'' empty records'')')                 &
     &               NUMREC1,LEN2_LOOKUP1-NUMREC1
            WRITE(6,'(''Fieldsfile file2 contains '',i5,'' fields '',   &
     &                ''and '',I5,'' empty records'')')                 &
     &               NUMREC2,LEN2_LOOKUP2-NUMREC2
            WRITE(7,'(''Fieldsfile file2 contains '',i5,'' fields '',   &
     &                ''and '',I5,'' empty records'')')                 &
     &               NUMREC2,LEN2_LOOKUP2-NUMREC2
            IF( NUMREC1 == NUMREC2)THEN
              WRITE(6,*) 'Files contain same number of fields'
              WRITE(7,*) 'Files contain same number of fields'
              IF(LEN2_LOOKUP1 == NUMREC1)THEN
                WRITE(6,*) 'Empty records at the end of file1 ',        &
     &                     'have probably been removed by convieee'
                WRITE(7,*) 'Empty records at the end of file1 ',        &
     &                     'have probably been removed by convieee'
              ELSE IF(LEN2_LOOKUP2 == NUMREC2)THEN
                WRITE(6,*) 'Empty records at the end of file2 ',        &
     &                     'have probably been removed by convieee'
                WRITE(7,*) 'Empty records at the end of file2 ',        &
     &                     'have probably been removed by convieee'
              END IF
            END IF
          END IF
        END IF

! Build cross reference index

! Initialise arrays
        DO i = 1, len2_lookup1
          index(i) = IMDI
          lmissing1(i) = .TRUE.
        END DO
        
        DO i = 1, len2_lookup2
          lmissing2(i) = .TRUE.
        END DO
        

! For every value in the first look-up table...
        DO i = 1, len2_lookup1

!...check whether it exists in the second look-up table. If so, 
!   is the value in the second table already being matched by a prior value in
!   lookup1  (i.e. lmissing2 = FALSE). If not, store the 
!   index and set the lmissing1 and lmissing2 logicals to false
search1:  DO j = 1, len2_lookup2
            IF (lookup1(item_code,i) == lookup2(item_code,j)) THEN
              IF(lmissing1(i).AND.lmissing2(j)) THEN
                index(i)     = j
                lmissing1(i) = .FALSE.
                lmissing2(j) = .FALSE.
                EXIT search1
              END IF
            END IF
          END DO search1
  
        END DO


        NMISSING1=0
        DO I=1,LEN2_LOOKUP1
          IF(LMISSING1(I).AND.LOOKUP1(1,I) /= -99)THEN
            NMISSING1=NMISSING1+1
            WRITE(6,'(''WARNING: Field '',I5,'' of file1 '',            &
     &                            ''has no match in file2'')') I
          END IF
        END DO
        NMISSING2=0
        DO I=1,LEN2_LOOKUP2
          IF(LMISSING2(I).AND.LOOKUP2(1,I) /= -99)THEN
            NMISSING2=NMISSING2+1
            WRITE(6,'(''WARNING: Field '',I5,'' of file2 '',            &
     &                            ''has no match in file1'')') I
          END IF
        END DO

        n_ignored = 0
        K=0
        DO I=1,NUMREC1
          IF(.NOT.LMISSING1(I).AND.LOOKUP1(1,I) /= -99)THEN
            DO J=1,LEN1_LOOKUP1
              IF(LOOKUP1(J,I) /= LOOKUP2(J,INDEX(I)))THEN
                K=K+1

! Ignore certain items in the lookup
                DO ig = 1, num_lookup_ignore
                  IF (j == lookup_ignore(ig) ) THEN
                    n_ignored = n_ignored + 1
                  END IF
                END DO
                
                ID1=LOOKUP1(J,I)
                ID2=LOOKUP2(J,INDEX(I))
                IF (J >= 46 .AND. J <= 64) THEN
                  WRITE(6,'(''Header1: '',I5,'' Header2: '',I5,         &
     &                      '' Item: '',I3,'' Values: '',F12.5,F12.5)') &
     &                   I,INDEX(I),J,RD1,RD2
                ELSE
                  WRITE(6,'(''Header1: '',I5,'' Header2: '',I5,         &
     &                    '' Item: '',I3,'' Values: '',I10,I10)')         &
     &                   I,INDEX(I),J,ID1,ID2
                END IF
              END IF
            END DO
          END IF
        END DO

        IF(K == 0) WRITE(6,*) 'OK'
        WRITE(7,'(''Number of fields in file 1 = '',I5)') NUMREC1
        WRITE(7,'(''Number of fields in file 2 = '',I5)') NUMREC2
        WRITE(7,'(''Number of fields compared  = '',I5)')               &
     &        NUMREC1-NMISSING1
        IF(nmissing1 /= 0 .AND. ignore_missing_fields)THEN
          WRITE(7,'(''Number of fields from file 1 omitted from '',     &
     &              ''comparison = '',I5)') nmissing1
        ELSE IF (nmissing1 /= 0) THEN
          WRITE(7,'(''Number of fields from file 1 missing from '',     &
     &              ''file2 = '',I5)') nmissing1
        END IF
        IF(nmissing2 /= 0 .AND. ignore_missing_fields)THEN
          WRITE(7,'(''Number of fields from file 2 omitted from '',     &
     &              ''comparison = '',I5)') nmissing2
        ELSE IF (nmissing2 /= 0) THEN
          WRITE(7,'(''Number of fields from file 2 missing from '',     &
     &              ''file1 = '',I5)') nmissing2
        END IF
        IF (n_ignored /= 0) THEN
          WRITE (7,'(''Number of differences in lookup which are '',    &
     &              ''ignored = '',I5)') n_ignored
        END IF
        WRITE(8,'(a,a,i7)') 'LOOKUP:                     ',             &
     &             'Number of differences = ',K

! Number of lookup differences is K, total number of fatal differences is L
        L = L + K - n_ignored
! If we're not ignoring missing fields, add them in to the number of fatal
! differences
        IF (.NOT. ignore_missing_fields) THEN
          l = l + nmissing1 + nmissing2
        END IF
      END IF

!L 16. Compare data fields

!L Get decompostion information
!L ----------------------------
!L The two files should have same resolution
!L TOT_LEVELS not used in SX
        CALL decompose(ROW_LENGTH1, P_ROWS1,                      &
     &                       0,0,TOT_LEVELS)
        CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)

! If it is an instantenous atmospheric dump
      IF(((FIXHD1(5) == 1).AND.(FIXHD2(5) == 1)).AND.                   &
     &   ((FIXHD1(2) == 1).AND.(FIXHD2(2) == 1)))THEN

!L Pulling out land-sea mask
!  -------------------------
! This assumes land-sea mask of 2 files are identical
        LAND_MASK_FOUND=.FALSE.
!       Check whether one exists before trying to read it
        DO I=1,LEN2_LOOKUP1
          IF (LOOKUP1(ITEM_CODE,i)  ==  30) THEN
            LAND_MASK_FOUND=.TRUE.
          END IF
        END DO
        IF (LAND_MASK_FOUND) THEN
! DEPENDS ON: read_land_sea
          CALL READ_LAND_SEA(NFTIN1,RCODE,LOOKUP1,                      &
     &      LEN1_LOOKUP1,LEN2_LOOKUP1,FIXHD1,LEN_FIXHD1)


!         Check error code from read_land_sea
          IF (RCODE /= -1.0) THEN
            WRITE (6,*) 'Error in READ_LAND_SEA.'
            WRITE (6,*) 'Return code from READ_LAND_SEA ',rcode
            ICODE = 200
            WRITE (CMESSAGE,*) 'DRLANDF1 : Error in READ_LAND_SEA.'
            GO TO 9999
          ENDIF
        END IF
      ENDIF
! Extract information of grid type for each file from STASHmaster

      WRITE(6,*)' '
      WRITE(6,*)'DATA FIELDS:'

      M=0
      N=0
      NDIFFER(:) = 0
      DO I=1,NUMREC1   ! Begin loop over number of fields in file1

        S_ITEM_CODE=MOD(LOOKUP1(42,I),1000)
        SECTION=(LOOKUP1(42,I)-S_ITEM_CODE)/1000
        IF(FIXHD1(12) >= 305)THEN
          MODEL=LOOKUP1(45,I)
        ELSEIF(S_ITEM_CODE <= 100.OR.                                   &
     &        (S_ITEM_CODE >= 200.AND.S_ITEM_CODE <= 205))THEN
          MODEL=1
        ELSEIF((S_ITEM_CODE >  100.AND.S_ITEM_CODE <= 176).OR.          &
     &         (S_ITEM_CODE >= 180.AND.S_ITEM_CODE <= 200))THEN
          MODEL=2
        ELSEIF((S_ITEM_CODE >= 177.AND.S_ITEM_CODE <= 179).OR.          &
     &         (S_ITEM_CODE >= 210.AND.S_ITEM_CODE <= 212))THEN
          MODEL=3
        END IF

!       All diagnostics under model code of 10 are in section 20
!       of Atmos StashMaster file.
        If (MODEL == 10) Then
          MODEL = 1
        End If

! DEPENDS ON: exppxc
        PHRASE=EXPPXC(MODEL,SECTION,S_ITEM_CODE,                        &
     &                ICODE,CMESSAGE)
        IF(ICODE /= 0)THEN
           PHRASE='NON-STANDARD FIELD'
           ICODE = 0
        END IF

        IF(.NOT.LMISSING1(I))THEN
           M=INDEX(I)

        IF((LOOKUP1(42,I) == 28.OR.LOOKUP1(42,I) == 29).AND.            &
     &     (FIXHD1(12) /= FIXHD2(12).AND.                               &
     &     (FIXHD1(12) >= 400.OR.FIXHD2(12) >= 400)))THEN
         LEN_FIELD=MIN0(LOOKUP1(15,I),LOOKUP2(15,M))
        ELSE
          LEN_FIELD=LOOKUP1(15,I)
        END IF
        IF(FIXHD1(12) <  0.AND.FIXHD1(5) /= 3)LOOKUP1(30,I)=0
        IF(FIXHD2(12) <  0.AND.FIXHD2(5) /= 3)LOOKUP2(30,M)=0
        IF (LOOKUP1(1,I) /= -99 .AND. LOOKUP2(1,M) /= -99) THEN

        FIELD_ITEM1=MOD(LOOKUP1(42,I),1000)
        FIELD_SECT1=(LOOKUP1(42,I)-FIELD_ITEM1)/1000
        FIELD_MODEL1=LOOKUP1(45,I)

        If (FIELD_MODEL1 == 10) Then
          FIELD_MODEL1 = 1
        End If

! DEPENDS ON: exppxi
        GRID_TYPE1=EXPPXI(FIELD_MODEL1,FIELD_SECT1,FIELD_ITEM1,         &
     &                          ppx_grid_type,                          &
     &                          ICODE,CMESSAGE)
          IF ((GRID_TYPE1 <  0).OR.(GRID_TYPE1 >  100)) THEN
            GRID_TYPE1=1
!           WRITE(6,*)'COMPARE: CANNOT GET GRID_TYPE1 INFO'
!           WRITE(6,*)'         FROM STASHMASTER FOR '
!           WRITE(6,*)'         SECTION ',FIELD_SECT1,
!    &                ' ITEM ',FIELD_ITEM1
!           WRITE(6,*)'         GRID_TYPE1 SET TO 1 (NORMAL GRID)'
          ENDIF
        FIELD_ITEM2=MOD(LOOKUP2(42,M),1000)
        FIELD_SECT2=(LOOKUP2(42,M)-FIELD_ITEM2)/1000
        FIELD_MODEL2=LOOKUP2(45,M)

        If (FIELD_MODEL2 == 10) Then
          FIELD_MODEL2 = 1
        End If

! DEPENDS ON: exppxi
        GRID_TYPE2=EXPPXI(FIELD_MODEL2,FIELD_SECT2,FIELD_ITEM2,         &
     &                          ppx_grid_type,                          &
     &                          ICODE,CMESSAGE)
          IF ((GRID_TYPE2 <  0).OR.(GRID_TYPE2 >  100)) THEN
            GRID_TYPE2=1
!           WRITE(6,*)'COMPARE: CANNOT GET GRID_TYPE2 INFO'
!           WRITE(6,*)'         FROM STASHMASTER FOR '
!           WRITE(6,*)'         SECTION ',FIELD_SECT2,
!    &                ' ITEM ',FIELD_ITEM2
!           WRITE(6,*)'         GRID_TYPE2 SET TO 1 (NORMAL GRID)'
          ENDIF

        PACK_CODE1 = MOD(LOOKUP1(21,I),10)
        PACK_CODE2 = MOD(LOOKUP2(21,M),10)

       lblrec_1=lookup1(15, i)
       lblrec_2=lookup2(15, m)

        if (((pack_code1 == 1 .or. pack_code2 == 1) .or.                &
     &       (pack_code1 == 4 .or. pack_code2 == 4)) .and.              &
     &  expand /= 1) then

        ELSEIF (PACK_CODE1 == 3 .OR. PACK_CODE2 == 3) THEN

          WRITE(6,*)                                                    &
     &    'Field No ',I,' not compared. GRIB data not supported.'

        ELSE

!         Since, for LBC data, header definition for pre Vn5.0
!         and Vn5.2- is different, a check is to be done.
          IF ((GRID_TYPE1 == ppx_atm_lbc_theta).OR.                     &
     &        (GRID_TYPE1 == ppx_atm_lbc_u).or.                         &
     &        (GRID_TYPE1 == ppx_atm_lbc_v).or.                         &
     &        (GRID_TYPE1 == ppx_atm_lbc_orog).or.                      &
     &        (GRID_TYPE1 == ppx_ocn_lbc_theta).or.                     &
     &        (GRID_TYPE1 == ppx_ocn_lbc_u)) then

            IF ((FIXHD1(12) >= 502).AND.                                &
     &          (FIXHD2(12) >= 502)) THEN
              IF ((LOOKUP1(41,i)/10000) /=                              &
     &            (LOOKUP2(41,i)/10000)) THEN
                write(6,*)'COMPARE: Fields with different rims',        &
     &                    '         not yet supported'

                CALL EREPORT('COMPARE', 1014,                           &
     &           'Fields with different rims not yet supported')

              ENDIF

              RIMWIDTHA(Rima_type_norm)=LOOKUP1(41,i)/10000
              haloY1=MOD(LOOKUP1(41,i),10000)
              haloY1=haloY1/100
              haloX1=MOD(haloY1,100)

              haloY2=MOD(LOOKUP2(41,i),10000)
              haloY2=haloY2/100
              haloX2=MOD(haloY2,100)

            ENDIF
          ENDIF
      IF((LOOKUP1(39,I) ==  1 .AND. LOOKUP2(39,M) ==  1).OR.            &
     &   (LOOKUP1(39,I) == -1 .AND. LOOKUP2(39,M) == -1))THEN
! This is a REAL field


!     field is atmos/ocean LB field
!     Decomposition is compulsory as each field may have
!     different halo size
      IF ((GRID_TYPE1 == ppx_atm_lbc_theta).OR.                         &
     &    (GRID_TYPE1 == ppx_atm_lbc_u).OR.                             &
     &    (GRID_TYPE1 == ppx_atm_lbc_v).OR.                             &
     &    (GRID_TYPE1 == ppx_atm_lbc_orog).OR.                          &
     &    (GRID_TYPE1 == ppx_ocn_lbc_theta).OR.                         &
     &    (GRID_TYPE1 == ppx_ocn_lbc_u)) THEN

        CURRENT_DECOMP_TYPE=-1

        CALL decompose(ROW_LENGTH1,P_ROWS1,                       &
     &                       haloX1, haloY1,                            &
     &                       TOT_LEVELS)
        CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
! DEPENDS ON: dervsize
        CALL DERVSIZE(ICODE,CMESSAGE)

      ENDIF


! DEPENDS ON: readflds
        CALL READFLDS(NFTIN1,1,I,LOOKUP1,LEN1_LOOKUP1,                  &
     &               R_D1,MAX_FIELD_SIZE1,FIXHD1,                       &
     &               expand,icode,cmessage)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('COMPARE',CMESSAGE,ICODE,NFTIN1)


!       field is atmos/ocean LB field
        IF ((GRID_TYPE2 == ppx_atm_lbc_theta).OR.                       &
     &      (GRID_TYPE2 == ppx_atm_lbc_u).or.                           &
     &      (GRID_TYPE2 == ppx_atm_lbc_v).or.                           &
     &      (GRID_TYPE2 == ppx_atm_lbc_orog).or.                        &
     &      (GRID_TYPE2 == ppx_ocn_lbc_theta).or.                       &
     &      (GRID_TYPE2 == ppx_ocn_lbc_u)) then

          CURRENT_DECOMP_TYPE=-1

          CALL decompose(ROW_LENGTH2,P_ROWS2,                     &
     &                         haloX2, haloY2,                          &
     &                         TOT_LEVELS)
          CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
! DEPENDS ON: dervsize
          CALL DERVSIZE(ICODE,CMESSAGE)
        ENDIF

! DEPENDS ON: readflds
        CALL READFLDS(NFTIN2,1,M,LOOKUP2,LEN1_LOOKUP2,                  &
     &                R_D2,MAX_FIELD_SIZE2,FIXHD2,                      &
     &               expand,icode,cmessage)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('COMPARE',CMESSAGE,ICODE,NFTIN1)

      ELSE IF((LOOKUP1(39,I) ==  2 .AND. LOOKUP2(39,M) ==  2).OR.       &
     &        (LOOKUP1(39,I) == -2 .AND. LOOKUP2(39,M) == -2))THEN
! This is an INTEGER field

!     field is atmos/ocean LB field
      IF ((GRID_TYPE1 == ppx_atm_lbc_theta).OR.                         &
     &    (GRID_TYPE1 == ppx_atm_lbc_u).OR.                             &
     &    (GRID_TYPE1 == ppx_atm_lbc_v).OR.                             &
     &    (GRID_TYPE1 == ppx_atm_lbc_orog).OR.                          &
     &    (GRID_TYPE1 == ppx_ocn_lbc_theta).OR.                         &
     &    (GRID_TYPE1 == ppx_ocn_lbc_u)) THEN

        CURRENT_DECOMP_TYPE=-1

        CALL decompose(ROW_LENGTH1,P_ROWS1,                       &
     &                       haloX1, haloY1,                            &
     &                       TOT_LEVELS)
        CALL CHANGE_DECOMPOSITION(4,ICODE)
        CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
! DEPENDS ON: dervsize
        CALL DERVSIZE(ICODE,CMESSAGE)

      ENDIF


! DEPENDS ON: readflds
        CALL READFLDS(NFTIN1,1,I,LOOKUP1,LEN1_LOOKUP1,                  &
     &                I_D1,MAX_FIELD_SIZE1,FIXHD1,                      &
     &               expand,icode,cmessage)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('COMPARE',CMESSAGE,ICODE,NFTIN1)

!       field is atmos/ocean LB field
        IF ((GRID_TYPE2 == ppx_atm_lbc_theta).OR.                       &
     &      (GRID_TYPE2 == ppx_atm_lbc_u).or.                           &
     &      (GRID_TYPE2 == ppx_atm_lbc_v).or.                           &
     &      (GRID_TYPE2 == ppx_atm_lbc_orog).or.                        &
     &      (GRID_TYPE2 == ppx_ocn_lbc_theta).or.                       &
     &      (GRID_TYPE2 == ppx_ocn_lbc_u)) then

          CURRENT_DECOMP_TYPE=-1

          CALL decompose(ROW_LENGTH2,P_ROWS2,                     &
     &                         haloX2, haloY2,                          &
     &                         TOT_LEVELS)
          CALL CHANGE_DECOMPOSITION(4,ICODE)
! DEPENDS ON: dervsize
          CALL DERVSIZE(ICODE,CMESSAGE)
        ENDIF

! DEPENDS ON: readflds
       CALL READFLDS(NFTIN2,1,M,LOOKUP2,LEN1_LOOKUP2,                   &
     &                I_D2,MAX_FIELD_SIZE2,FIXHD2,                      &
     &               expand,icode,cmessage)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('COMPARE',CMESSAGE,ICODE,NFTIN1)

      ELSE IF((LOOKUP1(39,I) ==  3 .AND. LOOKUP2(39,M) ==  3).OR.       &
     &        (LOOKUP1(39,I) == -3 .AND. LOOKUP2(39,M) == -3))THEN
! This is an LOGICAL field

!     field is atmos/ocean LB field
      IF ((GRID_TYPE1 == ppx_atm_lbc_theta).OR.                         &
     &    (GRID_TYPE1 == ppx_atm_lbc_u).OR.                             &
     &    (GRID_TYPE1 == ppx_atm_lbc_v).OR.                             &
     &    (GRID_TYPE1 == ppx_atm_lbc_orog).OR.                          &
     &    (GRID_TYPE1 == ppx_ocn_lbc_theta).OR.                         &
     &    (GRID_TYPE1 == ppx_ocn_lbc_u)) THEN

        CURRENT_DECOMP_TYPE=-1

        CALL decompose(ROW_LENGTH1,P_ROWS1,                       &
     &                       haloX1, haloY1,                            &
     &                       TOT_LEVELS)
        CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
! DEPENDS ON: dervsize
        CALL DERVSIZE(ICODE,CMESSAGE)
      ENDIF


! DEPENDS ON: readflds
        CALL READFLDS(NFTIN1,1,I,LOOKUP1,LEN1_LOOKUP1,                  &
     &                L_D1,MAX_FIELD_SIZE1,FIXHD1,                      &
     &               expand,icode,cmessage)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('COMPARE',CMESSAGE,ICODE,NFTIN1)

!       field is atmos/ocean LB field
        IF ((GRID_TYPE2 == ppx_atm_lbc_theta).OR.                       &
     &      (GRID_TYPE2 == ppx_atm_lbc_u).or.                           &
     &      (GRID_TYPE2 == ppx_atm_lbc_v).or.                           &
     &      (GRID_TYPE2 == ppx_atm_lbc_orog).or.                        &
     &      (GRID_TYPE2 == ppx_ocn_lbc_theta).or.                       &
     &      (GRID_TYPE2 == ppx_ocn_lbc_u)) then

          CURRENT_DECOMP_TYPE=-1

          CALL decompose(ROW_LENGTH2,P_ROWS2,                     &
     &                         haloX2, haloY2,                          &
     &                         TOT_LEVELS)
          CALL CHANGE_DECOMPOSITION(decomp_smexe,ICODE)
! DEPENDS ON: dervsize
          CALL DERVSIZE(ICODE,CMESSAGE)
        ENDIF

! DEPENDS ON: readflds
       CALL READFLDS(NFTIN2,1,M,LOOKUP2,LEN1_LOOKUP2,                   &
     &                L_D2,MAX_FIELD_SIZE2,FIXHD2,                      &
     &               expand,icode,cmessage)
! DEPENDS ON: abort_io
        IF(ICODE /= 0)CALL ABORT_IO('COMPARE',CMESSAGE,ICODE,NFTIN1)

      ELSE
! This is an unrecognized field
          WRITE(6,*)                                                    &
     &     'Field No ',I,' not compared. Unrecognized type.'

      ENDIF
       if (((pack_code1 == 1 .or. pack_code2 == 1).or.                  &
     &      (pack_code1 == 4 .or. pack_code2 == 4)).and.                &
     &  expand == 1) then

         LEN_FIELD=LOOKUP1(15,I)
       endif
       lookup1(15, i)=lblrec_1
       lookup2(15, M)=lblrec_2


        WRITE(6,*)' '
        write(6,'(/''Field '',i5,'' : Stash Code '',i5,                 &
     &   '' : '',a,'' : Level '',i4/)')                                 &
     &   i, lookup1(42,i), phrase, lookup1(33,i)
        write(10,'(/''Field '',i5,'' : Stash Code '',i5,                &
     &   '' : '',a)') i, lookup1(42,i), phrase
        IF(GRID_TYPE1 >  100)THEN
          WRITE(10,*)'Grid type = information not available'
        ELSE
          WRITE(10,'(a,i4)')'Grid type = ',GRID_TYPE1
        ENDIF

        RMS_F1=0.0
        RMS_F2=0.0
        RMS_DIFF=0.0
        K=0
!       Real
        IF (LOOKUP1(39,I) == 1 .AND. LOOKUP2(39,M) == 1) THEN
          MAX_DIFF=0.
          DO J=1,LEN_FIELD
            DIFF(J)='.'
            IF(R_D1(J) /= R_D2(J))THEN
              k=k+1
              if(k <= 10) then
                write(6,'(a,i6,2(e25.15,'' ('',z16,'')''))')            &
     &            'ITEM=',j,r_d1(j),r_d1(j),r_d2(j),r_d2(j)
              endif
                RD1=R_D1(J)
                RD2=R_D2(J)
                MAX_DIFF=MAX(MAX_DIFF,ABS(RD1-RD2))
                IF(MAX_DIFF == ABS(RD1-RD2)) MAX_J = J
                if(rd1 == 0.) then
                  if(rd2 == 0.) then
                    diff_per=0.
                  else
                    diff_per=(abs(rd1-rd2)/abs(rd2))*100
                  endif
                else
                  diff_per=(abs(rd1-rd2)/abs(rd1))*100
                endif
              IF (DIFF_PER  >   10.0) DIFF(J)="#"
              IF (DIFF_PER  <   10.0) DIFF(J)="X"
              IF (DIFF_PER  <   1.0) DIFF(J)="O"
              IF (DIFF_PER  <   0.1) DIFF(J)="o"
              IF (DIFF_PER  <   0.01) DIFF(J)=":"
                RMS_F1=RMS_F1+(R_D1(J)*R_D1(J))
                RMS_F2=RMS_F2+(R_D2(J)*R_D2(J))
                RMS_DIFF=RMS_DIFF+(R_D1(J)-R_D2(J))*(R_D1(J)-R_D2(J))
            else
                RMS_F1=RMS_F1+(R_D1(J)*R_D1(J))
                RMS_F2=RMS_F2+(R_D2(J)*R_D2(J))
                RMS_DIFF=RMS_DIFF+(R_D1(J)-R_D2(J))*(R_D1(J)-R_D2(J))
            endif
          ENDDO
          IF (K /= 0) THEN
            WRITE(6,*)'NUMBER OF DIFFERENT VALUES = ',K
            WRITE(6,*)'MAXIMUM DIFFERENCE= ',MAX_DIFF,' AT PT. ',MAX_J
            RMS_F1=SQRT(RMS_F1/LEN_FIELD)
            RMS_F2=SQRT(RMS_F2/LEN_FIELD)
            RMS_DIFF=SQRT(RMS_DIFF/LEN_FIELD)
            WRITE(6,*) 'RMS FIELD1 : ',RMS_F1
            WRITE(6,*) 'RMS FIELD2 : ',RMS_F2
            DIFF_PER=ABS(RMS_F1-RMS_F2)
            WRITE(6,*) 'Difference: ',DIFF_PER                          &
     &                ,' RMS_difference: ',RMS_DIFF
            rd1=diff_per
            if(rms_f1 /= 0.) then
              diff_per=(diff_per/rms_f1)*100
              write(6,'(''Field '',i5,'' has a Difference between'',    &
     &         '' the RMS Values of '',e10.5,'' which is '',f10.3,      &
     &         '' Percent of Field 1, whose RMS Value is '',e10.5)')    &
     &         i, rd1, diff_per, rms_f1
              write(6,*) 'Difference as % of RMS FIELD1= ',DIFF_PER
            else if(rms_f2 /= 0.) then
              diff_per=(diff_per/rms_f2)*100
              write(6,'(''Field '',i5,'' has a Difference between'',    &
     &         '' the RMS Values of '',e10.5,'' which is '',f10.3,      &
     &         '' Percent of Field 2, whose RMS Value is '',e10.5)')    &
     &         i, rd1, diff_per, rms_f2
              write(6,*) 'Difference as % of RMS FIELD2= ',DIFF_PER
            endif
!
            if (diff_per  >   5) THEN
              WRITE(6,*)
              WRITE(6,*)
              WRITE(6,*) '************** WARNING ********************'
              WRITE(6,*) '***** LARGE DIFFERENCE ENCOUNTERED ********'
              WRITE(6,*) '*******************************************'
              WRITE(6,*)
            ENDIF

            KEY='# d>10%  ;  X 10%>d>1%  ;  O 1%>d>0.1%  ; '//          &
     &          'o 0.1%>d>0.01%  ;  : d<0.01%  ;  . d=0%  ; '//         &
     &          '~ no data (sea-points)'

! In some cases, such as non-standard fields, grid type is not
! specified. Random number result.  In these cases, diff map will
! be printed if field size is (73*96).
! Grid-type 21 will have its new converted diff map of size
! (nrow*ncol) which is the same as land-mask it's been applied to
! Grid-type other then 21 have its map dimension lookup(18)
! *lookup(19)
            IF(((GRID_TYPE1 <  60).AND.                                 &
     &         (GRID_TYPE1 /= 31).AND.                                  &
     &         (GRID_TYPE1 /= 32).AND.                                  &
     &         (GRID_TYPE1 /= 47).AND.                                  &
     &         (GRID_TYPE1 /= 51))                                      &
     &      .OR.                                                        &
     &        ((GRID_TYPE1 >  100).AND.                                 &
     &         (LOOKUP1(18,I) >  0).AND.(LOOKUP1(18,I) <  300).AND.     &
     &         (LOOKUP1(19,I) >  0).AND.(LOOKUP1(19,I) <  300)))THEN
! Only certain grid types are suitable for difference maps
! For real data
              IF((GRID_TYPE1 == 21).AND.(GRID_TYPE2 == 21)) THEN
! DEPENDS ON: print_dif_map
                CALL PRINT_DIF_MAP(DIFF,P_ROWS1,ROW_LENGTH1,KEY,        &
     &            GRID_TYPE1,GRID_TYPE2)
              ELSE IF(GRID_TYPE1 == GRID_TYPE2.and.                     &
     &            ((GRID_TYPE2 == ppx_atm_lbc_theta).OR.                &
     &            (GRID_TYPE2 == ppx_atm_lbc_u).or.                     &
     &            (GRID_TYPE2 == ppx_atm_lbc_v).or.                     &
     &            (GRID_TYPE2 == ppx_atm_lbc_orog).or.                  &
     &            (GRID_TYPE2 == ppx_ocn_lbc_theta).or.                 &
     &            (GRID_TYPE2 == ppx_ocn_lbc_u))) then

                write(10,*)'LBC field: use arbitrary grid dimension'
                write(10,*)'80 columns wide to display differences'
! DEPENDS ON: print_dif_map
               CALL PRINT_DIF_MAP(DIFF,LEN_FIELD/80,80,KEY,             &
     &           GRID_TYPE1,GRID_TYPE2)
              ELSE
! DEPENDS ON: print_dif_map
               CALL PRINT_DIF_MAP(DIFF,LOOKUP1(18,I),LOOKUP1(19,I),KEY, &
     &           GRID_TYPE1,GRID_TYPE2)
              ENDIF
            ELSE
              write(6,*) 'Difference map not printed'
              write(6,*) 'Grid Type not suitable for difference maps'
              write(6,*) 'Grid Type = ',LOOKUP1(16,I)
            ENDIF

          ELSE
            WRITE(6,*)'OK'
            WRITE(10,*)'OK'
          ENDIF
!       Integer
        ELSE IF (LOOKUP1(39,I) == 2 .AND. LOOKUP2(39,M) == 2) THEN
          IMAX_DIFF=0
          DO J=1,LEN_FIELD
            DIFF(J)='.'
            IF(I_D1(J) /= I_D2(J))THEN
              K=K+1
              if (k <= n_diff) write(6,*)'item=',j,i_d1(j),i_d2(j)
                ID1=I_D1(J)
                ID2=I_D2(J)
                IMAX_DIFF=MAX(IMAX_DIFF,ABS(ID1-ID2))

                IF (ID1  ==  0) THEN
                  IF (ID1  ==  ID2) THEN
                    DIFF_PER=0.0
                  ELSE
                    DIFF_PER=100.0
                  ENDIF
                ELSE
                  DIFF_PER=(REAL(ABS(ID1-ID2))/REAL(ABS(ID1)))*100.0
                ENDIF
              IF (DIFF_PER  >   10.0) DIFF(J)="#"
              IF (DIFF_PER  <   10.0) DIFF(J)="X"
              IF (DIFF_PER  <   1.0) DIFF(J)="O"
              IF (DIFF_PER  <   0.1) DIFF(J)="o"
              IF (DIFF_PER  <   0.01) DIFF(J)=":"
            ENDIF
          ENDDO
          IF (K /= 0) THEN
            write(6,'(''Field '',i5,'' has '',i5,'//&
                ' '' INTEGER Differences'', '//&
                ' '' with a Maximum Difference of '',i20)') &
                i, k, imax_diff
            WRITE(6,*)'NUMBER OF DIFFERENT VALUES = ',K
            WRITE(6,*)'MAXIMUM DIFFERENCE= ',IMAX_DIFF
            KEY='# d>10%  ;  X 10%>d>1%  ;  O 1%>d>0.1%  ; '//          &
     &          'o 0.1%>d>0.01%  ;  : d<0.01%  ;  . d=0%  ; '//         &
     &          '~ no data (sea-points)'

! In some cases, such as non-standard fields, grid type is not
! specified. Random number result.  In these cases, diff map will
! be printed if field size is (73*96).
! Grid-type 21 will have its new converted diff map of size
! (nrow*ncol) which is the same as land-mask it's been applied to
! Grid-type other then 21 have its map dimension lookup(18)
! *lookup(19)
            IF(((GRID_TYPE1 <  60).AND.                                 &
     &         (GRID_TYPE1 /= 31).AND.                                  &
     &         (GRID_TYPE1 /= 32).AND.                                  &
     &         (GRID_TYPE1 /= 47).AND.                                  &
     &         (GRID_TYPE1 /= 51))                                      &
     &      .OR.                                                        &
     &        ((GRID_TYPE1 >  100).AND.                                 &
     &         (LOOKUP1(18,I) >  0).AND.(LOOKUP1(18,I) <  300).AND.     &
     &         (LOOKUP1(19,I) >  0).AND.(LOOKUP1(19,I) <  300)))THEN
! Only certain grid types are suitable for difference maps
! For integer data
              IF((GRID_TYPE1 == 21).AND.(GRID_TYPE2 == 21)) THEN
! DEPENDS ON: print_dif_map
                CALL PRINT_DIF_MAP(DIFF,P_ROWS1,ROW_LENGTH1,KEY,        &
     &            GRID_TYPE1,GRID_TYPE2)
              ELSE
! DEPENDS ON: print_dif_map
                CALL PRINT_DIF_MAP(DIFF,LOOKUP1(18,I),LOOKUP1(19,I),KEY,&
     &            GRID_TYPE1,GRID_TYPE2)
              ENDIF
            ENDIF

          ELSE
           write(6,'(''Field '',i5,'//&
               ' '' has '',i5,'' INTEGER Differences'')') i, k
            WRITE(6,*)'OK'
            WRITE(10,*)'OK'
          ENDIF
!       Logical
        ELSE IF (LOOKUP1(39,I) == 3 .AND. LOOKUP2(39,M) == 3) THEN
          DO J=1,LEN_FIELD
            DIFF(J)='.'
            IF (L_D1(J).NEQV.L_D2(J)) THEN
              K=K+1
              LD1=L_D1(J)
              LD2=L_D2(J)
              IF (K <= N_DIFF) WRITE(6,*)'ITEM=',J,LD1,LD2
              DIFF(J)="#"
            ENDIF
          ENDDO
          IF (K /= 0) THEN
            write(6,'(''Field '',i5,'' has '',i5,'//&
                ' '' LOGICAL Differences'')') i, k
            WRITE(6,*)'NUMBER OF DIFFERENT VALUES = ',K
            KEY='# Different values  ;  . identical  ; '//              &
     &          '~ no data (sea-points)'

! In some cases, such as non-standard fields, grid type is not
! specified. Random number result.  In these cases, diff map will
! be printed if field size is (73*96).
! Grid-type 21 will have its new converted diff map of size
! (nrow*ncol) which is the same as land-mask it's been applied to
! Grid-type other then 21 have its map dimension lookup(18)
! *lookup(19)
            IF(((GRID_TYPE1 <  60).AND.                                 &
     &         (GRID_TYPE1 /= 31).AND.                                  &
     &         (GRID_TYPE1 /= 32).AND.                                  &
     &         (GRID_TYPE1 /= 47).AND.                                  &
     &         (GRID_TYPE1 /= 51))                                      &
     &      .OR.                                                        &
     &        ((GRID_TYPE1 >  100).AND.                                 &
     &         (LOOKUP1(18,I) >  0).AND.(LOOKUP1(18,I) <  300).AND.     &
     &         (LOOKUP1(19,I) >  0).AND.(LOOKUP1(19,I) <  300)))THEN
! Only certain grid types are suitable for difference maps
! For logical data
              IF((GRID_TYPE1 == 21).AND.(GRID_TYPE2 == 21)) THEN
! DEPENDS ON: print_dif_map
                CALL PRINT_DIF_MAP(DIFF,P_ROWS1,ROW_LENGTH1,KEY,        &
     &            GRID_TYPE1,GRID_TYPE2)
              ELSE
! DEPENDS ON: print_dif_map
                CALL PRINT_DIF_MAP(DIFF,LOOKUP1(18,I),LOOKUP1(19,I),KEY,&
     &            GRID_TYPE1,GRID_TYPE2)
              ENDIF
            ENDIF
          ELSE
           write(6,'(''Field '',i5,'' has '',i5,'//&
               ' '' LOGICAL Differences'')') i, k
            WRITE(6,*)'OK'
            WRITE(10,*)'OK'
          ENDIF
!       Real Timeseries
        ELSE IF (LOOKUP1(39,I) == -1 .AND. LOOKUP2(39,M) == -1) THEN
          MAX_DIFF=0.
          DO J=1,LEN_FIELD
            IF(R_D1(J) /= R_D2(J))THEN
              MAX_DIFF=AMAX1(MAX_DIFF,ABS(R_D1(J)-R_D2(J)))
              K=K+1
              IF(K <= 10)WRITE(6,*)'ITEM=',J,R_D1(J),R_D2(J)
            ENDIF
          ENDDO
          IF (K /= 0) THEN
            WRITE(6,*)'NUMBER OF DIFFERENT VALUES = ',K
            WRITE(6,*)'MAXIMUM DIFFERENCE= ',MAX_DIFF
          ELSE
            WRITE(6,*)'OK'
          ENDIF
!       Integer Timeseries
        ELSE IF (LOOKUP1(39,I) == -2 .AND. LOOKUP2(39,M) == -2) THEN
          IMAX_DIFF=0
          DO J=1,LEN_FIELD
            IF (I_D1(J) /= I_D2(J)) THEN
              K=K+1
              ID1=I_D1(J)
              ID2=I_D2(J)
              IMAX_DIFF=MAX(IMAX_DIFF,IABS(ID1-ID2))
              IF (K <= N_DIFF) WRITE(6,*)'ITEM=',J,ID1,ID2
            ENDIF
          ENDDO
          IF (K /= 0) THEN
            WRITE(6,*)'NUMBER OF DIFFERENT VALUES = ',K
            WRITE(6,*)'MAXIMUM DIFFERENCE= ',IMAX_DIFF
          ELSE
            WRITE(6,*)'OK'
          ENDIF
!       Logical Timeseries
        ELSE IF (LOOKUP1(39,I) == -3 .AND. LOOKUP2(39,M) == -3) THEN
          DO J=1,LEN_FIELD
            IF (L_D1(J).NEQV.L_D2(J)) THEN
              K=K+1
              LD1=L_D1(J)
              LD2=L_D2(J)
              IF (K <= N_DIFF) WRITE(6,*)'ITEM=',J,LD1,LD2
            ENDIF
          ENDDO
          IF (K /= 0) THEN
            WRITE(6,*)'NUMBER OF DIFFERENT VALUES = ',K
          ELSE
            WRITE(6,*)'OK'
          ENDIF
        ELSE
          WRITE(6,*)                                                    &
     &    'Field No ',I,' not compared. Different Data Type Numbers ?'
        ENDIF
        WRITE(6,*)' '
        IF(K /= 0)THEN
          NDIFFER(I)=K
          N=N+1
        END IF
        L=L+K
        END IF

      ENDIF
      ENDIF
      ENDDO     !End loop over number of fields

! Output remainder of summary information
      WRITE(8,'(a,a,i7)') 'DATA FIELDS:                ',               &
     &             'Number of fields with differences = ',N
      DO I = 1,NUMREC1   ! Begin loop over number of fields in file1
        IF(LOOKUP1(1,I) /= -99)THEN
          S_ITEM_CODE=MOD(LOOKUP1(42,I),1000)
          SECTION=(LOOKUP1(42,I)-S_ITEM_CODE)/1000
          IF(FIXHD2(12) >= 305)THEN
            MODEL=LOOKUP1(45,I)
          ELSEIF(S_ITEM_CODE <= 100.OR.                                 &
     &          (S_ITEM_CODE >= 200.AND.S_ITEM_CODE <= 205))THEN
            MODEL=1
          ELSEIF((S_ITEM_CODE >  100.AND.S_ITEM_CODE <= 176).OR.        &
     &           (S_ITEM_CODE >= 180.AND.S_ITEM_CODE <= 200))THEN
            MODEL=2
          ELSEIF((S_ITEM_CODE >= 177.AND.S_ITEM_CODE <= 179).OR.        &
     &           (S_ITEM_CODE >= 210.AND.S_ITEM_CODE <= 212))THEN
            MODEL=3
          END IF

!         All diagnostics under model code of 10 are in section 20
!         of Atmos StashMaster file.
          If (MODEL == 10) Then
            MODEL = 1
          End If

! DEPENDS ON: exppxc
          PHRASE=EXPPXC(MODEL,SECTION,S_ITEM_CODE,                      &
     &                  ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
             PHRASE='NON-STANDARD FIELD'
             ICODE = 0
          END IF
          IF(LMISSING1(I))THEN
            WRITE(8,'(/''Field '',i5,'' : Stash Code '',i5,'' : '',a,'//   &
                ' '' : No equivalent in file2'')')                 &
                I,LOOKUP1(42,I),PHRASE
          ELSE IF(NDIFFER(I) /= 0)THEN
            WRITE(8,'(/''Field '',I5,'' : Stash Code '',I5, '//&
                ' '' : '',A,'' : Number of differences = '',I8)')  &
                I, LOOKUP1(42,I), PHRASE, NDIFFER(I)
          END IF
        END IF
      END DO
      DO I = 1,NUMREC2   ! Begin loop over number of fields in file2
        IF(LMISSING2(I).AND.LOOKUP2(1,I) /= -99)THEN
          S_ITEM_CODE=MOD(LOOKUP2(42,I),1000)
          SECTION=(LOOKUP2(42,I)-S_ITEM_CODE)/1000
          IF(FIXHD2(12) >= 305)THEN
            MODEL=LOOKUP2(45,I)
          ELSEIF(S_ITEM_CODE <= 100.OR.                                 &
     &          (S_ITEM_CODE >= 200.AND.S_ITEM_CODE <= 205))THEN
            MODEL=1
          ELSEIF((S_ITEM_CODE >  100.AND.S_ITEM_CODE <= 176).OR.        &
     &           (S_ITEM_CODE >= 180.AND.S_ITEM_CODE <= 200))THEN
            MODEL=2
          ELSEIF((S_ITEM_CODE >= 177.AND.S_ITEM_CODE <= 179).OR.        &
     &           (S_ITEM_CODE >= 210.AND.S_ITEM_CODE <= 212))THEN
            MODEL=3
          END IF

          If (MODEL == 10) Then
            MODEL = 1
          End If

! DEPENDS ON: exppxc
          PHRASE=EXPPXC(MODEL,SECTION,S_ITEM_CODE,                      &
     &                  ICODE,CMESSAGE)
          IF(ICODE /= 0)THEN
             PHRASE='NON-STANDARD FIELD'
             ICODE = 0
          END IF
          WRITE(8,'(/''Field '',i5,'' : Stash Code '',i5,'' : '',a,'//     &
              ' '' : No equivalent in file1'')')                   &
              I,LOOKUP2(42,I),PHRASE
        ENDIF
      END DO
      CLOSE(10)
      IF(L == 0)THEN
        IF (num_lookup_ignore > 0) THEN
          ! Use format statement which has maximum length of lookup.
          WRITE(8,'(a,64i5)') ' files compare, ignoring Fixed Length Header' &
                              //' and lookup items ',lookup_ignore
        ELSE
         WRITE(8,'(A)') ' files compare, ignoring Fixed Length Header'
        END IF
      ELSE
        WRITE(8,'(a)')' files DO NOT compare'
      ENDIF
      WRITE(7,*)' '
      CLOSE(7)
      CLOSE(8)

 9999 CONTINUE
      RETURN
      END SUBROUTINE COMPARE

