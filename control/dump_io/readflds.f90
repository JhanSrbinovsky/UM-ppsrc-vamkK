! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Reads in a number of fields from UM format file
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Dump I/O

! Subroutine Interface:
SUBROUTINE readflds(nftin, number_of_fields,                      &
                    first_field,lookup,len1_lookup_arg,           &
                    d1,disused,fixhd,                             &




                    icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IO
USE PrintStatus_mod
USE UM_ParVars
USE rimtypes
USE lbc_mod
USE lookup_addresses
USE Submodel_Mod

USE cppxref_mod, ONLY:                                             & 
    ppx_grid_type,ppx_halo_type,                                   &
    ppx_atm_rim,ppx_wam_rim,                                       &
    ppx_atm_lbc_orog,ppx_atm_compressed,                           &
    ppx_atm_tall,ppx_atm_cuall,ppx_atm_cvall

IMPLICIT NONE

! Description:
!  Reads in NUMBER_OF_FIELDS fields from file on unit NFTIN,
!  starting at field number FIRST_FIELD. The data is returned
!  in the D1 array.

! Subroutine Arguments:

INTEGER                                                           &
  nftin                                                           &
                         ! IN: unit number to read data from
, number_of_fields                                                &
                         ! IN: number of fields to read in
, first_field                                                     &
                         ! IN: first field to read in
, len1_lookup_arg                                                 &
                         ! IN: first dimension of LOOKUP table
, lookup(len1_lookup_arg,*)                                       &
                         ! IN: lookup table starting at field 1
, disused                                                         &
                         ! IN: Not used since move to MPP
, fixhd(*)                                                        &
                         ! IN: fixed length header
, expand                                                          &
                         ! IN: (=1 if WGDOS or RLE packed data
                         !      is to be expanded)
                         ! Only used for small execs etc
, icode                  ! OUT: return code

REAL                                                              &
  d1(*)                  ! OUT: array to return the data in

CHARACTER(LEN=80)                                                      &
  cmessage               ! OUT: Error message if ICODE <> 0

! COMMON blocks and PARAMETERs
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
! D1_ADDR start
      ! Information for accessing D1 addressing array
      ! Number of items of info needed for each object and maximum
      ! number of objects in D1 -

      ! Number of items of information in D1 addressing array
      INTEGER,PARAMETER:: D1_LIST_LEN=17

! Names of items in D1 addressing array. Update D1_LIST_LEN above if
! items added

      ! Prognostic, Diagnostic, Secondary or other
      INTEGER,PARAMETER:: d1_object_type    = 1 ! Internal model id
      INTEGER,PARAMETER:: d1_imodl          = 2  ! Internal model id
      INTEGER,PARAMETER:: d1_section        = 3  ! Section
      INTEGER,PARAMETER:: d1_item           = 4  ! Item
      INTEGER,PARAMETER:: d1_address        = 5  ! Address in D1
      INTEGER,PARAMETER:: d1_length         = 6  ! Record length
      INTEGER,PARAMETER:: d1_grid_type      = 7  ! Grid type
      INTEGER,PARAMETER:: d1_no_levels      = 8  ! Number of levels

      ! Stash list number for diags. -1 for progs
      INTEGER,PARAMETER:: d1_stlist_no      = 9

      ! Pointer to dump header lookup table
      INTEGER,PARAMETER:: d1_lookup_ptr     = 10

      INTEGER,PARAMETER:: d1_north_code     = 11 ! Northern row
      INTEGER,PARAMETER:: d1_south_code     = 12 ! Southern row
      INTEGER,PARAMETER:: d1_east_code      = 13 ! Eastern row
      INTEGER,PARAMETER:: d1_west_code      = 14 ! Western row
      INTEGER,PARAMETER:: d1_gridpoint_code = 15 ! gridpoint info
      INTEGER,PARAMETER:: d1_proc_no_code   = 16 ! Processing Code
      INTEGER,PARAMETER:: d1_halo_type      = 17 ! Halo width type

      ! Types of items for d1_type

      INTEGER,PARAMETER:: prognostic = 0
      INTEGER,PARAMETER:: diagnostic = 1
      INTEGER,PARAMETER:: secondary  = 2
      INTEGER,PARAMETER:: other      = 3

! D1_ADDR end
! TYPSIZE start
!   Description:
!     This file contains sizes needed for dynamic allocation of
!   main data arrays within the model. Sizes read in from the user
!   interface via NAMELISTs are passed by /COMMON/. Other control
!   sizes that are fundamental in the definition of data structures
!   are assigned by PARAMETER statements.
!
!   Declarations for the NLSIZES namelist are also held in the module
!   nlsizes_namelist_mod. That module is currently only used by the
!   reconfiguration, while the UM uses this include file.
!
! All sizes
! Not dependent on sub-model
! DATA IN NAMLST#x MEMBER OF THE JOB LIBRARY
! ATMOS START
! Main sizes of fields for each submodel
! Grid-related sizes for ATMOSPHERE submodel.
INTEGER :: ROW_LENGTH           ! No of points per local row
INTEGER :: global_ROW_LENGTH    ! Points per global row
INTEGER :: ROWS                 ! No of local (theta) rows
INTEGER :: global_ROWS          ! No of global (theta) rows
INTEGER :: MODEL_LEVELS         ! No of model levels
INTEGER :: LAND_FIELD           ! No of land points in field
INTEGER :: NTILES               ! No of land surface tiles
INTEGER :: NICE                 ! No. of sea ice thickness categories
INTEGER :: NICE_USE             ! No. of sea ice categories used fully
                                !  in surface exchange and radiation
                                !  (If nice>1 & nice_use=1, categories only 
                                !  partially used in surface exchange)

! Physics-related sizes for ATMOSPHERE submodel
INTEGER :: WET_LEVELS          ! No of moist-levels
INTEGER :: CLOUD_LEVELS        ! No of cloud-levels
INTEGER :: ST_LEVELS           ! No of soil temperature levels
INTEGER :: SM_LEVELS           ! No of soil moisture levels
INTEGER :: BL_LEVELS           ! No of boundary-layer-levels
INTEGER :: OZONE_LEVELS        ! No of ozone-levels
INTEGER :: TPPS_OZONE_LEVELS   ! No of tropopause-ozone-levels
INTEGER :: RIVER_ROWS          ! No of rows for river routing
INTEGER :: RIVER_ROW_LENGTH    ! Row length for river routing
! Dynamics-related sizes for ATMOSPHERE submodel

INTEGER :: TR_LEVELS            ! No of tracer-levels
INTEGER :: TR_VARS              ! No of passive tracers
INTEGER :: TR_LBC_VARS          ! No of tracers in lbcs 
INTEGER :: TR_UKCA              ! No of UKCA tracers
INTEGER :: TR_LBC_UKCA          ! No of UKCA tracer lbcs 

! For Small executables

! Grid related sizes for data structure
! Data structure sizes for ATMOSPHERE submodel
INTEGER :: A_PROG_LOOKUP     ! No of prognostic fields
INTEGER :: A_PROG_LEN        ! Total length of prog fields
INTEGER :: A_LEN_INTHD       ! Length of INTEGER header
INTEGER :: A_LEN_REALHD      ! Length of REAL header
INTEGER :: A_LEN2_LEVDEPC    ! No of LEVEL-dependent arrays
INTEGER :: A_LEN2_ROWDEPC    ! No of ROW-dependent arrays
INTEGER :: A_LEN2_COLDEPC    ! No of COLUMN-dependent arrays
INTEGER :: A_LEN2_FLDDEPC    ! No of FIELD arrays
INTEGER :: A_LEN_EXTCNST     ! No of EXTRA scalar constants
INTEGER :: A_LEN_CFI1        ! Length of compressed fld index 1
INTEGER :: A_LEN_CFI2        ! Length of compressed fld index 2
INTEGER :: A_LEN_CFI3        ! Length of compressed fld index 3
! atmos end

! Data structure sizes for ATMOSPHERE ANCILLARY file control
! routines
INTEGER :: NANCIL_LOOKUPSA  ! Max no of fields to be read

! Data structure sizes for ATMOSPHERE INTERFACE file control
! routines
INTEGER :: N_INTF_A          ! No of atmosphere interface areas
INTEGER :: MAX_INTF_MODEL_LEVELS ! Max no of model levs in all areas
INTEGER :: MAX_LBCROW_LENGTH ! Max no of lbc row length in all areas
INTEGER :: MAX_LBCROWS ! Max no of lbc rows in all areas

!  Data structure sizes for ATMOSPHERE BOUNDARY file control
! routines

! Sizes applicable to all configurations (DUMPS/FIELDSFILE)

INTEGER :: PP_LEN_INTHD   ! Length of PP file integer header
INTEGER :: PP_LEN_REALHD  ! Length of PP file real    header


      ! Grid related sizes for COUPLING between ATMOS and OCEAN
      ! submodels [For MPP, sizes are 'global' values over all
      ! PEs.Also needed for river routing]
      INTEGER:: AOCPL_IMT                ! Ocean rowlength
      INTEGER:: AOCPL_JMT                ! Ocean no. of rows
      INTEGER:: AOCPL_ROW_LENGTH         ! Atmos rowlength
      INTEGER:: AOCPL_P_ROWS             ! Atmos no. of p rows

      COMMON/SIZE_AOCPL/                                                &
        AOCPL_IMT, AOCPL_JMT, AOCPL_ROW_LENGTH, AOCPL_P_ROWS

! Other sizes passed from namelist into common blocks
! Any additions to this common block must be mirrored in nlsizes_namelist_mod.
COMMON/NLSIZES/                                                     &
    ROW_LENGTH,global_ROW_LENGTH,ROWS,global_ROWS,                  &
    LAND_FIELD,MODEL_LEVELS,WET_LEVELS,                             &
    NTILES, NICE, NICE_USE,                                         &
    CLOUD_LEVELS,TR_LEVELS,ST_LEVELS,SM_LEVELS,BL_LEVELS,           &
    OZONE_LEVELS,TPPS_OZONE_LEVELS,TR_VARS,TR_LBC_VARS,             &
    TR_UKCA,TR_LBC_UKCA,RIVER_ROWS,RIVER_ROW_LENGTH,                &
    A_PROG_LOOKUP,A_PROG_LEN,                                       &
    A_LEN_INTHD,A_LEN_REALHD,                                       &
    A_LEN2_LEVDEPC,A_LEN2_ROWDEPC,A_LEN2_COLDEPC,                   &
    A_LEN2_FLDDEPC,A_LEN_EXTCNST,                                   &
    A_LEN_CFI1,A_LEN_CFI2,A_LEN_CFI3,                               &    
    NANCIL_LOOKUPSA,                                                &    
    N_INTF_A, MAX_INTF_MODEL_LEVELS, MAX_LBCROW_LENGTH,             &
    MAX_LBCROWS, PP_LEN_INTHD,PP_LEN_REALHD

!-----------------------------------------------------------------
! data in STASHC#x member of the job library

! Data structure sizes for ATMOSPHERE submodel (config dependent)
INTEGER :: A_LEN2_LOOKUP   ! Total no of fields (incl diags)
INTEGER :: A_LEN_DATA      ! Total no of words of data
INTEGER :: A_LEN_D1        ! Total no of words in atmos D1

! Size of main data array for this configuration

INTEGER :: LEN_TOT             ! Length of D1 array
INTEGER :: N_OBJ_D1_MAX         ! No of objects in D1 array

COMMON/STSIZES/                                                     &
    A_LEN2_LOOKUP,A_LEN_DATA,A_LEN_D1,                              &
    LEN_TOT,N_OBJ_D1_MAX
! global (ie. dump version) of *_LEN_DATA
INTEGER :: global_A_LEN_DATA

COMMON /MPP_STSIZES_extra/ global_A_LEN_DATA
! Sizes of Stash Auxillary Arrays and associated index arrays
! Initialised in UMINDEX and UMINDEX_A/O/W
INTEGER :: LEN_A_IXSTS
INTEGER :: LEN_A_SPSTS

COMMON /DSIZE_STS/                                                  &
    LEN_A_IXSTS, LEN_A_SPSTS
!     The number of land points is computed for each PE
!     before the addressing section. All prognostics on land
!     points in the D1 space are now dimensioned by the local
!     no of land points rather than the global no of land points.

      INTEGER:: global_land_field    !  Global no of land points
      INTEGER:: local_land_field     !  Local no of land points
      COMMON /mpp_landpts/ global_land_field,local_land_field
      ! ----------------------------------------------------------------
      ! extra variables not passed through user interface

      ! fundamental data sizes :
      ! Fundamental parameter  sizes of data structure
      ! Sizes applicable to all configurations (HISTORY FILE)

      ! Length of history file in dump
      INTEGER, PARAMETER :: LEN_DUMPHIST = 0

      ! Sizes applicable to all configurations (DUMPS/FIELDSFILE)
      ! Length of dump fixed header
      INTEGER, PARAMETER :: LEN_FIXHD = 256

      ! Size of a single LOOKUP header
      INTEGER, PARAMETER :: LEN1_LOOKUP  = 64
      INTEGER, PARAMETER :: MPP_LEN1_LOOKUP= 2

      ! Size of compressed LBC LOOKUP (only used internally and
      ! contains just the items which change between each set of LBCs
      INTEGER, PARAMETER :: LEN1_LBC_COMP_LOOKUP = 8

      ! Sizes applicable to all configurations (STASH)
      ! Moved to typstsz.h

      INTEGER:: INTF_LEN2_LEVDEPC !1st dim of interface out lev dep cons
      INTEGER:: INTF_LEN2_ROWDEPC !2nd dim of interface out Row dep cons
      INTEGER:: INTF_LEN2_COLDEPC !2nd dim of interface out Col dep cons
      
      COMMON /DSIZE/                                                    &
        INTF_LEN2_LEVDEPC,INTF_LEN2_ROWDEPC,INTF_LEN2_COLDEPC
      ! sub-model atmosphere   :
      ! Data structure sizes derived from grid size
      INTEGER:: A_LEN1_LEVDEPC ! IN: 1st dim of level  dep const
      INTEGER:: A_LEN1_ROWDEPC ! IN: 1st dim of row    dep const
      INTEGER:: A_LEN1_COLDEPC ! IN: 1st dim of column dep const
      INTEGER:: A_LEN1_FLDDEPC ! IN: 1st dim of field  dep const

      ! Data structure sizes for ATMOSPHERE INTERFACE file control
      ! routines
      INTEGER :: INTF_LOOKUPSA        ! No of interface lookups.
      COMMON /DSIZE_A/                                                  &
        A_LEN1_LEVDEPC,A_LEN1_FLDDEPC,A_LEN1_ROWDEPC,A_LEN1_COLDEPC,    &
        INTF_LOOKUPSA

      ! sub-model atmosphere   : derived sizes
      ! derived from model grid/levels. Arakawa B-grid

                                  ! Size of fields on THETA grid:
      INTEGER :: THETA_FIELD_SIZE     ! IN: with no halos
      INTEGER :: THETA_OFF_SIZE       ! IN: with simple halos
      INTEGER :: THETA_HALO_SIZE      ! IN: with extended halos

                                  ! Size of fields on U grid:
      INTEGER :: U_FIELD_SIZE         ! IN: with no halos
      INTEGER :: U_OFF_SIZE           ! IN: with simple halos
      INTEGER :: U_HALO_SIZE          ! IN: with extended halos

                                  ! Size of fields on V grid
      INTEGER :: V_FIELD_SIZE         ! IN: with no halos
      INTEGER :: V_OFF_SIZE           ! IN: with simple halos
      INTEGER :: V_HALO_SIZE          ! IN: with extended halos

      INTEGER :: N_ROWS               ! IN: No of V-rows
      INTEGER :: N_CCA_LEV            ! IN: No of CCA Levels
      COMMON/DRSIZE_A/                                                  &
        N_ROWS,N_CCA_LEV,THETA_FIELD_SIZE,U_FIELD_SIZE,V_FIELD_SIZE,    &
        THETA_OFF_SIZE,U_OFF_SIZE,V_OFF_SIZE,                           &
        THETA_HALO_SIZE,U_HALO_SIZE,V_HALO_SIZE
      ! boundary updating      : derived values
      ! Variables describing the Atmosphere Lateral Boundary Conditions
      ! Local (per processor) information


! TYPSIZE end

INTEGER, PARAMETER :: fh_lookupsize2 = 152

! Local variables

INTEGER                                                           &
  k                                                               &
                         ! loop over fields to read in
, d1_off                                                          &
                         ! local offset into D1 for this field
, pack_code                                                       &
                         ! packing code for field
, field_start                                                     &
                         ! location of field on disk
, data_size                                                       &
                         ! number of words of data on disk
                         ! (including padding for WFIO)
, data_read_size                                                  &
                         ! number of words to read from disk
, data_full_size                                                  &
                         ! number of words after any unpacking
, len_io                                                          &
                         ! number of words read from disk
, num_cray_words                                                  &
, num_unpack_values                                               &
, len_full_word                                                   &
, field_item                                                      &
                         ! Item number of field
, field_sect                                                      &
                         ! Section number of field
, field_model                                                     &
                         ! Model ID of field
, grid_type                                                       &
                         ! grid type code
, fld_type                                                        &
                         ! field type (P,U or V)
, halo_type                                                       &
                         ! halo type code
, i                                                               &
                         ! loop index
, local_len              ! size of field section on this PE


INTEGER                                                           &
  fake_d1_addr(d1_list_len)    ! Fake D1_ADDR record to
                               ! pass into read_multi
INTEGER unset                  ! unset values
PARAMETER (unset=-1)

! Functions

INTEGER get_fld_type,exppxi

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------

IF (lhook) CALL dr_hook('READFLDS',zhook_in,zhook_handle)

IF (fixhd(12) <  403) THEN
  WRITE(6,*) 'READFLDS: file created by UM version ', fixhd(12)
  icode=1
  cmessage='READFLDS: Cannot read files from before vn4.3'
  GO TO 9999
END IF

d1_off=0

DO k=first_field,first_field+number_of_fields-1

  pack_code=MOD((lookup(lbpack,k)),10)

! Set up the location of the field on disk
  field_start=lookup(lbegin,k) ! position of field in file
  IF (field_start <= 0) THEN
    WRITE(6,*) 'READFLDS: start address =',field_start
    icode=1
    cmessage='READFLDS: start address of field not given'
    GO TO 9999
  END IF

! data_size contains the number of words of data used to store the
! field on disk (needs to be halved if 32 bit packing has been used)
  IF (pack_code  ==  2) THEN
    data_size=(lookup(lblrec,k)+1)/2
  ELSE
    data_size=lookup(lblrec,k)
  END IF

! data_full_size is the number of words required to store the field
! in memory after any unpacking is done.
  ! This is to give buf the correct size in RDUNPCK, as
  ! buf will be the final expanded size of whole field
  ! including extra data
  ! warning lbext - may be -32768 missing value !
  IF ((pack_code == 4).AND.(lookup(lbext, k) >  0)) THEN
    data_full_size=MAX(lookup(lbrow, k)*lookup(lbnpt, k)          &
      +lookup(lbext, k) ,lookup(lblrec,k))
  ELSE
    data_full_size=MAX(lookup(lbrow, k)*lookup(lbnpt, k)          &
      ,lookup(lblrec,k))
  END IF

  IF ((lookup(lbrow,k) <  0).OR.                                  &
      (lookup(lbnpt,k) <  0)) THEN
    data_full_size=lookup(lblrec,k)
  END IF
! data_read_size contains the number of words to data that need to
! be read in for a field. Each field has extra words of dummy data
! added at the end to ensure each field starts on a disk sector
! boundary.   
  data_read_size = lookup(lbnrec,K)

! The last field on a dump does not have these extra words
! added.  So check against number of lookups and if next lookup is -99 assume
! we are at the end.  Otherwise assume we are the last field.
  IF (k < fixhd(fh_lookupsize2)) THEN
    IF (lookup(lbyr,K+1) == -99) THEN
      data_Read_Size = data_size
    END iF
  ELSE
    data_read_size = data_size
  END IF

  IF (data_read_size <  0) THEN
    WRITE(6,*) 'READFLDS: number of words to read =',             &
     data_read_size
    icode=1
    cmessage='READFLDS: number of words to read not given'
    GO TO 9999
  END IF

! data_full_size needs to be at least as big as data_read_size since
! it is used to dimension the BUF array in READ_MULTI.

  data_full_size = MAX(data_full_size, data_read_size)


! Move file pointer to the start of the field

  CALL setpos(nftin,field_start,icode)
  IF (icode  /=  0) THEN
    WRITE(6,*)                                                    &
     'READFLDS - SETPOS failed to move file pointer to ',         &
     field_start,' on unit ',nftin
    WRITE(6,*) 'SETPOS returned error code ',icode
    icode=100
    cmessage='SETPOS failed in READFLDS. See output.'
    GO TO 9999
  END IF

! Get some information about this field

  field_item=MOD(lookup(42,k),1000)
  field_sect=(lookup(42,k)-field_item)/1000
  field_model=lookup(45,k)
  IF (fixhd(5)  >=  6 .AND. fixhd(5)  <=  9) THEN
! Set grid_type and halo_type for ACobs and VARobs, Cx and CovStats
! (FIXHD(5)=6, 7, 8 and 9, respectively).
    grid_type = 1
    halo_type = halo_type_no_halo
  ELSE
! DEPENDS ON: exppxi
    grid_type=exppxi(field_model,field_sect,field_item,           &
                          ppx_grid_type,                          &
                          icode,cmessage)
! DEPENDS ON: exppxi
    halo_type=exppxi(field_model,field_sect,field_item,           &
                          ppx_halo_type,                          &
                          icode,cmessage)
    IF (icode  /=  0) THEN
      WRITE(6,*)                                                  &
        'READFLDS - EXPPXI failed to get PPXREF information ',    &
        'for field: '
      WRITE(6,*) 'Model ID : ',field_model
      WRITE(6,*) 'Section  : ',field_sect
      WRITE(6,*) 'Item     : ',field_item
      WRITE(6,*) 'Error code was ',icode
      WRITE(6,*) 'Error message was ',cmessage
      icode=1
      cmessage='READFLDS failed to get PPXREF information'
      GO TO 9999
    END IF
  END IF
! For Utitlies (except Makebc) we dont want to know about haloes.
! DEPENDS ON: get_fld_type
  fld_type=get_fld_type(grid_type) ! field type P,U or V

! Set up fake D1_ADDR record to describe data to be read in
! Only set those items actually required by read_multi
! Assume that no diagnostic type fields will be read.

  DO i=1,d1_list_len
    fake_d1_addr(i)=unset
  END DO

  fake_d1_addr(d1_object_type)=prognostic
  fake_d1_addr(d1_imodl)=field_model
  fake_d1_addr(d1_section)=field_sect
  fake_d1_addr(d1_item)=field_item
  fake_d1_addr(d1_halo_type)=halo_type
  fake_d1_addr(d1_length)=lasize(1,fld_type,halo_type)*           &
                          lasize(2,fld_type,halo_type)
  fake_d1_addr(d1_no_levels)=1

! Grid type - for LBCs we need some special logic...
  IF ( (lookup(lbhem,k)  >=  99) .AND.                            &
                                         ! This is a LBC field
       (lookup(lbhem,k)  <   1000)) THEN

    IF (field_model  ==  atmos_im) THEN
      IF (lookup(lbhem,k)  ==  99) THEN ! Old style LBCs
        fake_d1_addr(d1_grid_type)=ppx_atm_rim
      ELSE ! New style LBCs with different field types

        fake_d1_addr(d1_grid_type)=grid_type

        IF (grid_type  ==  ppx_atm_lbc_orog) THEN
          fake_d1_addr(d1_length)=                                &
            lenrima(fld_type,halo_type,rima_type_orog)*           &
            (lookup(lbhem,k)-100)
        ELSE
          fake_d1_addr(d1_length)=                                &
            lenrima(fld_type,halo_type,rima_type_norm)*           &
            (lookup(lbhem,k)-100)
        END IF ! IF (grid_type  ==  ppx_atm_lbc_orog)
        fake_d1_addr(d1_no_levels)=lookup(lbhem,k)-100
      END IF ! IF (LOOKUP(LBHEM,K)  ==  99)
    ELSE
      icode=2
      WRITE(6,*) 'READFLDS: Cannot process LBC for model type ',  &
                 field_model
      cmessage='READFLDS : Cannot read LBCS for this model type'
      GO TO 9999
    END IF

  ELSE IF (fixhd(5) == 4  .AND.                                   &
                                           ! Ancillary File
          (MOD( lookup(lbpack,k)/10, 10 ) == 0 )  .AND.           &
           grid_type == ppx_atm_compressed ) THEN

    ! Compressed in stashmaster but uncompressed in header
    SELECT CASE( fld_type )
      CASE( fld_type_p )
        fake_d1_addr(d1_grid_type) = ppx_atm_tall

      CASE( fld_type_u )
        fake_d1_addr(d1_grid_type) = ppx_atm_cuall

      CASE( fld_type_v )
        fake_d1_addr(d1_grid_type) = ppx_atm_cvall

    END SELECT

  ELSE ! not an LBC
    fake_d1_addr(d1_grid_type)=grid_type
  END IF

! Read the field from disk and distribute it over the processors

! DEPENDS ON: read_multi
  CALL read_multi(nftin,d1(d1_off+1),data_read_size,              &
                  data_full_size,len_io,local_len,                &
                  lookup(1,k),fixhd(12),fake_d1_addr,             &
                  fake_d1_addr(d1_no_levels),                     &
                  icode,cmessage)

! Check for I/O errors

  IF ((icode  /=  0) .OR. (len_io  /=  data_read_size)) THEN
    WRITE(6,*) 'READFLDS : Error reading in field ',k
    WRITE(6,*) 'UNIT ',nftin
    WRITE(6,*) 'MODEL ID ',field_model
    WRITE(6,*) 'SECTION ',field_sect
    WRITE(6,*) 'ITEM ',field_item
    IF (fixhd(5) <  6 .OR. fixhd(5) >  10) THEN

! DEPENDS ON: pr_look
      CALL pr_look(                                               &
              lookup,lookup,len1_lookup_arg,k)

    END IF
    cmessage='READFLDS: I/O error'
    GO TO 9999
  END IF

  d1_off=d1_off+local_len
  IF(local_len == 0)THEN
    d1_off=d1_off+lookup(lblrec,k)
  END IF
END DO ! K : loop over fields

9999 CONTINUE

IF (lhook) CALL dr_hook('READFLDS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE readflds

