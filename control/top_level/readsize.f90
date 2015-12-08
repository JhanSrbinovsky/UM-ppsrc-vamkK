! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine READSIZE --------------------------------------------
!LL
!LL    Purpose:
!LL  This routine reads the length required for the dimensioning
!LL  of all main data arrays that have dimensions set via questions
!LL  in the USER INTERFACE.
!LL  It is called by the shell program UM_SHELL and passes the
!LL  data via common block to that routine. Data is then passed by
!LL  argument lists to allow for proper dynamic allocation in the
!LL  portable model.
!LL
!LL Programming standard:
!LL
!LL ------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Top Level

SUBROUTINE READSIZE()


  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
! JULES
  USE ancil_info, ONLY: nsmax
  USE Check_IOStat_mod
  USE ereport_mod, ONLY : ereport
  USE Atmos_Max_Sizes
  USE UM_ParParams
  USE lbc_mod
  USE IO
  USE model_file
  USE UM_ParVars
  USE calc_ntiles_mod, ONLY : calc_ntiles
  USE read_jules_namelists_mod, ONLY: read_jules_nstypes
  USE switches, ONLY: l_aggregate, l_sice_multilayers, l_ctile
  USE bl_option_mod, ONLY: l_us_blsol
  USE nstypes, ONLY: npft, nnvg
  USE mym_option_mod, ONLY: tke_levels, shcu_levels
  USE rad_input_mod, ONLY: cusack_aero_hgt, aero_bl_levels
! module for aerosol emissions options
  USE run_aerosol_mod, ONLY: run_aerosol_check 

  IMPLICIT NONE

!  Local variables

      INTEGER                       :: ERRORSTATUS
      CHARACTER (LEN=*), PARAMETER  :: ROUTINENAME='READSIZE'
      CHARACTER (LEN=80)            :: CMESSAGE

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!LL  COMDECK: CNLSIZES
!LL
!LL  This calls the :                        COMDECK
!LL                 type declarations        TYPSIZE
!LL                 COMMON blocks            TYPSIZE
!LL                 NAMELISTs                NAMSIZE
!LL  for the array size variables defined by the User Interface.
!LL
!LL  Additions need to be added in all 3 positions in TYPSIZE and
!LL  NAMSIZE COMDECKS for sizes obtained from the UI. Other sizes are
!LL  derived from the UI set and only need entries in type declarations
!LL  and common blocks in TYPSIZE.
!LL
!LL  DATA is split into sections denoting which elements belong to
!LL  the set provided in the two UI provided files: STASHC# and NAMLST#.
!LL
!LL  A Common Block is used to allow for data to be read from
!LL  a NAMELIST.
!LL
!LL  This COMDECK should only be called by UM_SHELL and READSIZE.
!LL  Portability pre-processing does not allow arrays to be dimensioned
!LL  via variables on a common block within the same routine, but sizes
!LL  can be passed down as arguments for dynamic allocation at a lower
!LL  level of subroutine.
!LL
!LL  THIS COMDECK CANNOT CALL THE TYPSIZE COMDECKS AS *DEF CALLS WOULD
!LL  CHANGE THE CONTENTS OF THE NAMELIST. ALL NAMELIST ELEMENTS ARE
!LL  PROVIDED BY THE UI.
!LL
!LL  Model            Modification history:
!LL version  Date
!LL  3.2   05/05/93   M.CARTER: DECK CNLSIZES on MOD SET MC050593
!LL                   creation: added for dynamic allocation
!LL
!----------------------------------------------------------------------
!
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
!
! Description:
! NLSIZES namelist. Set up by UMUI and put in SIZES member of Job.
!
! This namelist is now duplicated in nlsizes_namelist_mod, along with the
! variable declarations required for it. Access to the NLSIZES namelist from 
! the reconfiguration should be made via that module, not here.
! The UM continues to use this file.
! ANY CHANGES TO THIS NAMELIST MUST BE MIRRORED IN NLSIZES_NAMELIST_MOD.
!
      NAMELIST/NLSIZES/                                                 &
     &  global_ROW_LENGTH,global_ROWS,LAND_FIELD,                       &
     &  MODEL_LEVELS,WET_LEVELS,                                        &
     &  CLOUD_LEVELS,TR_LEVELS,ST_LEVELS,BL_LEVELS,                     &
     &  OZONE_LEVELS,                                                   &
     &  PP_LEN_INTHD,PP_LEN_REALHD

      NAMELIST/ATM_SIZES/                                               &
     &  NICE, NSMAX, NANCIL_LOOKUPSA, NRIM_TIMESA

! NAMSIZE end
!*L --------------------- Comdeck: CENVIR   ----------------------------
!
!   Purpose: COMDECK defining Character enviroment variables used
!            by portable IO to open and close files
!
!
! Type declarations
!
      CHARACTER(LEN=8) FT_ENVIRON(199)  ! Array holding enviroment variables
!                                  for filenames
      INTEGER     LEN_FT_ENVIR(199) ! character length of each variable
!
!
!Common Blocks for character and integer arrays
!
      COMMON/CENVIR/FT_ENVIRON
      COMMON/CLENVIR/LEN_FT_ENVIR
!

      INTEGER, POINTER ::   fixhdr(:)
      INTEGER   len_io
      REAL      a
      INTEGER   icode
      INTEGER, ALLOCATABLE  :: inthdr(:) 
      INTEGER, PARAMETER    :: dumpunit=21 ! staple unit number for ASTART

!L 1.1 Read size information from namelists

      IF (lhook) CALL dr_hook('READSIZE',zhook_in,zhook_handle)

!     Initialise TPPS_OZONE_LEVELS from OZONE_LEVELS
!     Prevents failures when using 'preset=nan'
!     Currently not set through UMUI
      TPPS_OZONE_LEVELS = 0

! New sea ice scheme is still in development so set scheme to off here
! (This is not set in the UMUI at present)
      NICE_USE = 1

!     FILE ON UNIT 5 DOES NOT NEED TO BE OPENED.
!     READ THE UI PROVIDED NAMELIST OF ARRAY SIZES.
      READ (UNIT=5, NML=ATM_SIZES, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist ATM_SIZES")

! Read restart file header to get the value of the size of a number of arrays
      CALL model_file_open(dumpunit,ft_environ(dumpunit),                      &
          len_ft_envir(dumpunit),ioOpenReadOnly,ioNameInEnv)

      CALL loadHeader(dumpunit,FixedHeader)
      fixhdr=>attachHeader(dumpunit,FixedHeader)

      ALLOCATE (inthdr(fixhdr(101)))
      IF (fixhdr(100) >  0) THEN
        CALL buffin(dumpunit,inthdr,fixhdr(101))
      END IF

      global_ROW_LENGTH=inthdr(6)
      global_ROWS=inthdr(7)
      MODEL_LEVELS=inthdr(8)
! inthdr(21) contains the value used for "mdi"
      IF(inthdr(25).EQ.inthdr(21))THEN
       LAND_FIELD=0
      ELSE
       LAND_FIELD=inthdr(25)
      END IF
      IF(inthdr(9).EQ.inthdr(21))THEN
       WET_LEVELS=0
      ELSE
       WET_LEVELS=inthdr(9)
      END IF
      IF(inthdr(11).EQ.inthdr(21))THEN
       CLOUD_LEVELS=WET_LEVELS
      ELSE
       CLOUD_LEVELS=inthdr(11)
      END IF
      IF(inthdr(10).EQ.inthdr(21))THEN
       ST_LEVELS=0
      ELSE
       ST_LEVELS=inthdr(10)
      END IF
      IF(inthdr(28).EQ.inthdr(21))THEN
       SM_LEVELS=0
      ELSE
       SM_LEVELS=inthdr(28)
      END IF
      IF(inthdr(13).EQ.inthdr(21))THEN
       BL_LEVELS=0
      ELSE
       BL_LEVELS=inthdr(13)
      END IF
      IF(inthdr(26).EQ.inthdr(21))THEN
       OZONE_LEVELS=0
      ELSE
       OZONE_LEVELS=inthdr(26)
      END IF
      IF(inthdr(12).EQ.inthdr(21))THEN
       TR_LEVELS=0
      ELSE
       TR_LEVELS=inthdr(12)
      END IF
      IF(inthdr(14).EQ.inthdr(21))THEN 
       TR_VARS=0
      ELSE
       TR_VARS=inthdr(14)
      END IF
      IF(inthdr(19).EQ.inthdr(21))THEN 
       RIVER_ROW_LENGTH=0
      ELSE
       RIVER_ROW_LENGTH=inthdr(19)
      END IF
      IF(inthdr(20).EQ.inthdr(21))THEN 
       RIVER_ROWS=0
      ELSE
       RIVER_ROWS=inthdr(20)
      END IF

      A_LEN_INTHD=fixhdr(101)
      A_LEN_REALHD=fixhdr(106)
      PP_LEN_INTHD=fixhdr(101)
      PP_LEN_REALHD=fixhdr(106)
      IF(fixhdr(112).EQ.inthdr(21))THEN 
       A_LEN2_LEVDEPC=0
      ELSE
       A_LEN2_LEVDEPC=fixhdr(112)
      END IF
      IF(fixhdr(117).EQ.inthdr(21))THEN 
       A_LEN2_ROWDEPC=0
      ELSE
       A_LEN2_ROWDEPC=fixhdr(117)
      END IF
      IF(fixhdr(122).EQ.inthdr(21))THEN 
       A_LEN2_COLDEPC=0
      ELSE
       A_LEN2_COLDEPC=fixhdr(122)
      END IF
      IF(fixhdr(127).EQ.inthdr(21))THEN 
       A_LEN2_FLDDEPC=0
      ELSE
       A_LEN2_FLDDEPC=fixhdr(127)
      END IF
      IF(fixhdr(131).EQ.inthdr(21))THEN 
       A_LEN_EXTCNST=0
      ELSE
       A_LEN_EXTCNST=fixhdr(131)
      END IF
      IF(fixhdr(141).EQ.inthdr(21))THEN 
       A_LEN_CFI1=0
      ELSE
       A_LEN_CFI1=fixhdr(141)
      END IF
      IF(fixhdr(143).EQ.inthdr(21))THEN 
       A_LEN_CFI2=0
      ELSE
       A_LEN_CFI2=fixhdr(143)
      END IF
      IF(fixhdr(145).EQ.inthdr(21))THEN 
       A_LEN_CFI3=0
      ELSE
       A_LEN_CFI3=fixhdr(145)
      END IF

      DEALLOCATE(inthdr)

      CALL model_file_close(dumpunit,ft_environ(dumpunit),                     &
          len_ft_envir(dumpunit),ioNameInEnv)




! Check that dynamically provided sizes are not larger than their
! statically defined maximums.
      IF (GLOBAL_ROW_LENGTH > ROW_LENGTH_MAX) THEN
        ERRORSTATUS = 10
        CMESSAGE = 'Row length is larger than maximum defined in'       &
                 //' Atmos_Max_Sizes'

        CALL EREPORT( ROUTINENAME, ERRORSTATUS, CMESSAGE )
      END IF

      IF (GLOBAL_ROWS > ROWS_MAX) THEN
        ERRORSTATUS = 20
        CMESSAGE = 'Number of rows is larger than maximum defined in'   &
                 //' Atmos_Max_Sizes'

        CALL EREPORT( ROUTINENAME, ERRORSTATUS, CMESSAGE )
      END IF

      IF (MODEL_LEVELS > MODEL_LEVELS_MAX) THEN
        ERRORSTATUS = 30
        CMESSAGE = 'Number of levels is larger than maximum defined in' &
                 //' Atmos_Max_Sizes'

        CALL EREPORT( ROUTINENAME, ERRORSTATUS, CMESSAGE )
      END IF

      ntiles=9
      CALL read_jules_nstypes( 5 )
      CALL calc_ntiles(l_aggregate,npft,nnvg,ntiles)

! In preparation for NICE_USE being set in UMUI
      IF (NICE_USE /= 1 .AND. NICE_USE /= NICE ) THEN
        ERRORSTATUS = 40
        CMESSAGE = 'Problem with sea ice category settings' &
                 //'NICE_USE must either equal 1 or NICE'

        CALL EREPORT( ROUTINENAME, ERRORSTATUS, CMESSAGE )
      END IF

      ! Checks for sea ice settings - in preparation for these logicals being
      ! set in namelist. Perform this check in readsize() after nice has been
      ! read in by namelist and l_ctile and l_sice_multilayers have been read
      ! in from JULES_SWITCHES in readlsta()
      
      IF (.NOT. l_ctile .AND. nice_use > 1) THEN
        errorstatus = 10
        cmessage = 'Problem with sea ice settings:'                  &
             // 'Either set l_ctile=T or nice_use=1'
        CALL ereport( routinename, errorstatus, cmessage )
      END IF
      IF (l_sice_multilayers .AND. nice_use /= nice) THEN
        errorstatus = 20
        cmessage = 'Problem with sea ice settings:'                   &
             //'Either set l_sice_multilayers=F or nice_use=nice'
        CALL ereport( routinename, errorstatus, cmessage )
      END IF
      ! Check sea ice settings are compatible with boundary layer scheme
      IF ( (l_sice_multilayers .OR. nice_use > 1) .AND.              &
           .NOT. l_us_blsol ) THEN
        errorstatus = 103 
        cmessage    = 'Problem with sea ice settings:'//             &
             'l_us_blsol setting not compatible'
        CALL ereport(routinename, errorstatus, cmessage)
      END IF

! A negative value for "tke_levels" means it should default to bl_levels.
      IF (tke_levels < 0) THEN
        tke_levels = bl_levels
        cmessage = 'WARNING: the value of tke_levels has been '       &
                 //'reset to bl_levels'
        errorstatus = -1
        CALL ereport(RoutineName, errorstatus, cmessage)
      END IF
! A negative value for "shcu_levels" means it should default to tke_levels.
      IF (shcu_levels < 0) THEN
        shcu_levels = tke_levels
        cmessage = 'WARNING: the value of shcu_levels has been '       &
                 //'reset to tke_levels'
        errorstatus = -1
        CALL ereport(RoutineName, errorstatus, cmessage)
      END IF

! Set radiation aerosol switch that depends on bl_levels
      IF (cusack_aero_hgt==2) aero_bl_levels = bl_levels

! Check run_aerosol
      CALL run_aerosol_check

      IF (lhook) CALL dr_hook('READSIZE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE READSIZE
