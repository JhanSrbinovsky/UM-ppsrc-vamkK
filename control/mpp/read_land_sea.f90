! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Parallel UM : Reads in the local section of Land-Sea Mask.

! Subroutine Interface:
SUBROUTINE read_land_sea(nft,icode,lookup,loc_len1_lookup,       &
                         loc_len2_lookup,fixhd,loc_len_fixhd)


USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE IO
USE ereport_mod, ONLY : ereport
USE UM_ParVars
USE domain_params
USE lookup_addresses

IMPLICIT NONE

! Description:
!  This routine reads the land-sea mask (LSM) from the dump and puts
!  it in a COMMON block defined in IOVARS. It is required for
!  unpacking and packing fields which are stored compressed to
!  land points.

! Method:
!  The position of the LSM within the dump is found from examining
!  the LOOKUP headers, it is then read in, and the relevant part
!  of the field sent to each processor. The local number of land
!  points is counted, and the LAND_FIELD variable is reset to this
!  new value.
!  Note : Halos can contain land points - but only those halos
!         which are updated by SWAPBNDS.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Subroutine Arguments:

INTEGER, INTENT(IN) :: nft              ! IN : FORTRAN unit number
INTEGER, INTENT(IN) :: loc_len1_lookup  ! IN : Dimension of the LOOKUP array
INTEGER, INTENT(IN) :: loc_len2_lookup  ! IN : Dimension of the LOOKUP array
INTEGER, INTENT(IN) :: loc_len_fixhd    ! IN : Dimension of the FIXHD array

INTEGER, INTENT(IN) :: lookup(loc_len1_lookup,loc_len2_lookup)
                                            ! IN : LOOKUP array from dump header
INTEGER, INTENT(IN) :: fixhd(loc_len_fixhd) ! IN : FIXHD array from dump header

REAL, INTENT(OUT) ::   icode               ! OUT : Return code

! Parameters and Common blocks

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

! Local variables

INTEGER :: i,j,k,word_address,ipts,iproc,info,len_io
INTEGER :: landpts_local,local_off,global_off

LOGICAL::read_mask

INTEGER                      ::  errorstatus = 0
CHARACTER (LEN=80)           ::  cmessage = ' '
INTEGER                      ::  i_icode
REAL io_ret(3)

! Local Parameters
CHARACTER (LEN=*), PARAMETER ::  routinename = 'READ_LAND_SEA'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! --------------------------------------------------------------------

IF (lhook) CALL dr_hook('READ_LAND_SEA',zhook_in,zhook_handle)
icode=-1.0
ipts=0
len_io=0

! Find location of LSM in the dump

! disable read-broadcast for buffin by pe 0
CALL set_unit_bcast_flag(nft)
IF (mype  ==  0) THEN

  read_mask=.FALSE.

  DO i=1,loc_len2_lookup
    IF (lookup(item_code,i)  ==  30) THEN
      read_mask=.TRUE.
      EXIT
    END IF
  END DO

  IF (.NOT. read_mask) THEN
    WRITE(6,'(a,i4)')                                               &
    'Error in READ_LAND_SEA_MASK: Missing Field of Type ',item_code
    errorstatus = 1
    cmessage = 'Missing Field or Wrong Dumpfile'

    CALL ereport(routinename,errorstatus,cmessage)
  END IF

  k=i
  word_address=1
! Old Format dumpfiles
  IF((lookup(lbnrec,k) == 0) .OR.                                 &
! Prog lookups in dump before vn3.2:
    ((lookup(lbnrec,k) == imdi) .AND. (fixhd(12) <= 301))) THEN
! Dump and ancillary files
  word_address=1
  IF (i  >   1) THEN
    DO k=2,i
      IF(MOD((lookup(lbpack,k-1)),10) == 2) THEN
        ipts=(lookup(lblrec,k-1)+1)/2
      ELSE
        ipts=(lookup(lblrec,k-1))
      END IF
      word_address=word_address+ipts
    END DO
  END IF
  word_address=fixhd(160)+word_address-2
    ipts=lookup(lblrec, i)
  ELSE
! PP type files and new format Dumpfiles (vn4.4 onwards)
    word_address=lookup(lbegin,i)
! Use the stored round-up value
    ipts=lookup(lbnrec,i)
  END IF

  CALL setpos(nft,word_address,i_icode)
  IF (i_icode  /=  0) THEN
    WRITE(6,*) 'READ_LAND_SEA: Error Return from SETPOS',  &
               ' Status is ',i_icode

    errorstatus = 3
    cmessage = 'Error from SETPOS - see output for Status'

    CALL ereport(routinename,errorstatus,cmessage)
 
  END IF

! Read the LSM in to PE 0

!--check that there is space to read the data
  IF(ipts >  max2dfieldsize) THEN
    WRITE(6,'(a,i10,a,i10)')                                            &
    'READ_LAND_SEA_MASK: The number of Words to be Read ',ipts,       &
    ' is larger than the Buffer Size ',max2dfieldsize
    WRITE(6,'(a,i10)') 'Record length is ',lookup(lblrec, i)

    errorstatus = 4
    cmessage = 'Insufficient Space for Land Sea Mask'

    CALL ereport(routinename,errorstatus,cmessage)

  END IF

  CALL buffin(nft,atmos_landmask,ipts,                     &
                     len_io,icode)

END IF   ! (mype == 0)
CALL clear_unit_bcast_flag(nft)! Restore broadcast flag

! Broadcast the I/O Status to the other PE's

io_ret(1)=len_io
io_ret(2)=icode
io_ret(3)=ipts

CALL gc_rbcast(99, 3, 0, nproc, info, io_ret)

len_io=NINT(io_ret(1))
icode=io_ret(2)
ipts=io_ret(3)

! Check the I/O Status on all PE'e

IF ((icode /= -1.0).OR.(len_io /= ipts)) THEN
  WRITE(6,*)'ERROR READING DUMP ON UNIT ',nft
! DEPENDS ON: ioerror
  CALL ioerror('BUFFER IN FROM READ_LAND_SEA_MASK',               &
   icode,len_io,ipts)
  errorstatus = 1
  WRITE(cmessage,'(A,I5)') 'Error reading dump on using ',nft
  CALL ereport(routinename,errorstatus,cmessage)
END IF

! Broadcast the global LSM to all processors

CALL gc_ibcast(100,glsize(1,fld_type_p)*glsize(2,fld_type_p),     &
               0,nproc,info,atmos_landmask)


DO i=1,lasize(1,fld_type_p,halo_type_no_halo)*                    &
       lasize(2,fld_type_p,halo_type_no_halo)
    atmos_landmask_local(i)=.FALSE.
END DO

! Copy my local part of the full LSM into atmos_landmask_local

DO j=1,blsize(2,fld_type_p)

  local_off=(j-1)*lasize(1,fld_type_p,halo_type_no_halo)
  global_off=(j-1+datastart(2)-1)*glsize(1,fld_type_p)+           &
               datastart(1)-1

  DO i=1,blsize(1,fld_type_p)

    atmos_landmask_local(local_off+i)=                            &
      atmos_landmask(global_off+i)

  END DO ! i
END DO ! j


! Count the number of global land points

atmos_number_of_landpts=0
DO i=1,glsize(1,fld_type_p)*glsize(2,fld_type_p)
  IF (atmos_landmask(i))                                          &
      atmos_number_of_landpts=atmos_number_of_landpts+1
END DO

! Do a swap to get land points in halo areas

landpts_local=0
DO i=1,lasize(1,fld_type_p,halo_type_no_halo)*                    &
       lasize(2,fld_type_p,halo_type_no_halo)
  IF (atmos_landmask_local(i))                                    &
    landpts_local=landpts_local+1
END DO


IF (landpts_local  /=  land_field) THEN
  WRITE(6,*) 'PE ',mype,' : LAND_FIELD is being reset from ',     &
             land_field,' to ',landpts_local
  land_field=landpts_local
END IF

IF (lhook) CALL dr_hook('READ_LAND_SEA',zhook_out,zhook_handle)
RETURN
END SUBROUTINE read_land_sea

