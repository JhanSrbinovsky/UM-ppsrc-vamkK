! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Program MAKEBC : Top-level program to create boundary dataset
!                   from model analyses/dumps.
!
! Method : For each dump, boundary conditions are generated through
!          GEN_INTF for the area specified in the INFTCNST namelist.
!          This routine initialises various variables in TYPSIZE
!          before it can be used in the lower routines.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MakeBC

      PROGRAM main_makebc

      USE makebc_constants_mod
      USE check_iostat_mod

      USE ereport_mod, ONLY : ereport, ereport_finalise
      USE PrintStatus_mod
      USE UM_ParVars
      USE UM_Config, ONLY : &
          appInit, &
          exe_makebc, &
          exe_frames

      USE domain_params
      USE dust_parameters_mod, ONLY: l_dust
      USE um_input_control_mod,  ONLY:                                  &
           l_so2,             l_dms,             l_so4_aitken,          &
           l_so4_accu,        l_so4_diss,        l_nh3,                 &
           l_soot,            l_biomass,         l_ocff,                &
           l_nitrate,         lcal360

      USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
      USE cloud_inputs_mod, ONLY: l_pc2
      USE murk_inputs_mod,  ONLY: l_murk
      USE ppxlook_mod, ONLY: ppxrecs
      USE Submodel_Mod

      USE chsunits_mod, ONLY : nunits

 
      USE control_max_sizes, ONLY: max_n_intf_a
      USE missing_data_mod, ONLY: imdi
      USE filenamelength_mod, ONLY :                                    & 
          filenamelength
      USE version_mod, ONLY:                                            &
          nitemp, nsectp, nproftp, ntimep

      IMPLICIT NONE

      INTEGER :: internal_model   !  Internal Model Identifier
      INTEGER :: icode            !  Error code
      INTEGER :: errorstatus
      INTEGER :: me_gc
      INTEGER :: nproc_gc

      CHARACTER(LEN=80)           :: cmessage    !  Error Message
      CHARACTER(LEN=80),PARAMETER :: routinename='main_makebc'

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
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS


! TYPSTSZ end

      CHARACTER(LEN=*), PARAMETER :: ProgName = "MakeBC"


!
! Free and UKCA tracer related variables and namelist
!
      INTEGER            :: i
      INTEGER            :: no_tracers, no_tr_ukca
      INTEGER            :: tracers_active(max_tracers)
      INTEGER            :: tr_ukca_active(max_tr_ukca)

! new bits shifting dump2bound read to this routine
!     DUMP2BOUND namelist for MAKEBC Program
      INTEGER  :: n_dumps      !  No of model dumps
      INTEGER  :: nhours       !  No of hours between dumps
      INTEGER  :: um_versn     !  UM Version Boundary Dataset for

      INTEGER :: sub_hr_int   ! No of lbcs desired per hour
      INTEGER :: no_lams      ! No of lams to generate lbcs for
                              ! in this run

      INTEGER :: ndustbin_in  = 6   ! Number of dust bins in input data
      INTEGER :: ndustbin_out = 6   ! Number of dist bins in output data

      NAMELIST /DUMP2BOUND/ n_dumps,nhours,um_versn,lcal360,            &
              no_lams,sub_hr_int,                                       &
              l_pc2,l_murk,l_mcr_qcf2,l_mcr_qrain,l_mcr_qgraup,         &
              l_dust,l_so2,l_so4_aitken,l_so4_accu,l_so4_diss,l_dms,    &
              l_nh3,l_soot,l_biomass,l_ocff,l_nitrate,                  &
              tracers_active,tr_ukca_active,                            &
              ndustbin_in, ndustbin_out

! Number of optional lbcs (ie pc2, murk) required
      INTEGER :: num_optional_lbcs_out

      ! Input filename
      CHARACTER(LEN=filenamelength) :: nml_file

!-  End of header

! DEPENDS ON: timer
      CALL Timer ( ProgName, 1 )
! DEPENDS ON: initprintstatus
      CALL InitPrintStatus

      CALL gc_init(' ', me_gc, nproc_gc)
      CALL appInit(exe_makebc)
      WRITE (6,'(A)') ' ##########################################'
      WRITE (6,'(A)') ' Running MAKEBC Utility to create a'
      WRITE (6,'(A)') ' Boundary Dataset from Model Output'
      WRITE (6,'(A)') ' ##########################################'

      WRITE (6,'(A)') ' '
      icode = 0
      
      ! Close UNIT 5 and reopen as file.
      CLOSE(5)
      CALL get_file(5,nml_file,filenamelength,icode)
      OPEN(UNIT=5,FILE=nml_file)

!     Only Atmosphere Model catered for in um_submodel_init for makebc
! DEPENDS ON: um_submodel_init
      CALL um_submodel_init(icode)
      IF (icode >  0) THEN
        WRITE (cmessage,'(A)') 'Error in initialising sub-model information.'
        CALL ereport(routinename,icode,cmessage)
      END IF

!     Determine no of Atmos records in STASHmaster file
      ppxRecs=1
! DEPENDS ON: hdppxrf
      call hdppxrf(22,'STASHmaster_A',icode,cmessage)
      IF (icode >  0) THEN
        WRITE (6,*) 'Error in HDPPXRF for STASHmaster_A.'

        CALL ereport(routinename,icode,cmessage)
      END IF

!     Initialise variables in TYPSIZE as in wstlst
      nsects                = nsectp
      nitems                = nitemp
      n_req_items           = 1
      n_ppxrecs             = 1
      totitems              = 1
      nsttims               = ntimep
      nsttabl               = nproftp
      num_stash_pseudo      = 1
      num_pseudo_lists      = 1
      nstash_series_records = 1
      nstash_series_block   = 1
      ! n_obj_d1_max used by d1_addr later on.
      n_obj_d1_max          = 1

!     Dimensions of Headers in Boundary dataset
!     Integer/Real Constants
      pp_len_inthd   = 46
      pp_len_realhd  = 38
      a_len_inthd    = pp_len_inthd
      a_len_realhd   = pp_len_realhd

!     This sizes are used in dimensioning arrays in many include files but not
!     used (for example in typduma.h). These are eventually set in
!     loop_over_dumps.
      a_len1_levdepc = 1
      a_len2_levdepc = 1
      a_len1_rowdepc = 1
      a_len2_rowdepc = 1 
      a_len1_coldepc = 1 
      a_len2_coldepc = 1 
      a_len1_flddepc = 1 
      a_len2_flddepc = 1 
      a_len_extcnst  = 1 
      a_len_cfi1     = 1 
      a_len_cfi2     = 1 
      a_len_cfi3     = 1 
      a_len2_lookup  = 1 
      a_len_data     = 1 


!     Level Dependent Constants array (Second dimension)
      intf_len2_levdepc = 4
!     Row/Col Dependent Constants array (Second dimension)
      intf_len2_rowdepc = 2
      intf_len2_coldepc = 2
      
! Read in the DUMP2BOUND namelist
!     Defaults for DUMP2BOUND namelist
      n_dumps      = imdi
      nhours       = imdi
      um_versn     = 601
      lcal360      = .FALSE.
      sub_hr_int   = imdi
      no_lams      = 0
      l_pc2        = .FALSE.
      l_murk       = .FALSE.
      l_mcr_qcf2   = .FALSE.
      l_mcr_qrain  = .FALSE.
      l_mcr_qgraup = .FALSE.
! By default all aerosol LBCs are off       
      l_dust       = .FALSE.
      l_so2        = .FALSE.
      l_so4_aitken = .FALSE.
      l_so4_accu   = .FALSE.
      l_so4_diss   = .FALSE.
      l_dms        = .FALSE.
      l_nh3        = .FALSE.
      l_soot       = .FALSE.
      l_biomass    = .FALSE.
      l_ocff       = .FALSE.
      l_nitrate    = .FALSE.

! Set defaults for free/UKCA tracer variables 
      no_tracers        = 0
      tracers_active(:) = 0
      no_tr_ukca        = 0
      tr_ukca_active(:) = 0

!     Read in DUMP2BOUND namelist and print
      READ (UNIT=5, NML=DUMP2BOUND, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist DUMP2BOUND")
      WRITE (6,*) ' '
      WRITE (6,*) 'Namelist DUMPBOUND read in '
      WRITE (6,DUMP2BOUND)

      ! Print warning that these variables are now retired.
      IF (n_dumps    /= imdi .OR. &
          nhours     /= imdi .OR. &
          sub_hr_int /= imdi) THEN
       
        WRITE(cmessage,'(A)') &
          'Please remove n_dumps, nhours, and sub_hr_int from DUMP2BOUND.'
        icode = -1
        CALL ereport(routinename,icode,cmessage)
      
      END IF

      IF (no_lams > max_n_intf_a) THEN
        WRITE(cmessage,'(A,I2,A)') 'Too many LAMs (>',max_n_intf_a,') requested'
        errorstatus = 1
        CALL ereport(routinename, errorstatus, cmessage)
      ELSE IF (no_lams <= 0) THEN
        WRITE(cmessage,'(A,I2,A)') 'Zero (or negative) LAMs requested'
        errorstatus = 1
        CALL ereport(routinename, errorstatus, cmessage)
      ELSE
        ! No of areas requiring boundary conditions.
        ! Allow for multiple LAMs with number from dump2bound
        n_intf_a = no_lams
      END IF


!
! calculate number of tracer (free/UKCA) and print
      DO i=1,max_tracers
        IF(tracers_active(i)==1)THEN
          no_tracers=no_tracers+1
        END IF
      END DO
      DO i=1,max_tr_ukca
        IF(tr_ukca_active(i)==1)THEN
          no_tr_ukca=no_tr_ukca+1
        END IF
      END DO
      IF (printstatus >= prstatus_diag) THEN
         WRITE (6,*) ' '
         WRITE (6,*) 'tracers_active=',tracers_active
         WRITE (6,*) 'no_tracers=',no_tracers
         WRITE (6,*) 'tr_ukca_active=',tr_ukca_active
         WRITE (6,*) 'no_tr_ukca=',no_tr_ukca
      END IF



!     Derive data lengths.
! DEPENDS ON: derv_intf_a
      CALL derv_intf_a (max_intf_model_levels,max_lbcrow_length,max_lbcrows,  &
                        n_intf_a)

!     Length of super arrays.
      len_a_spsts =1
      len_a_ixsts =1

!     No of data types for which boundary conditions required
!     Assume no tracer variables and set default number of
!     optional lbcs to zero
      tr_vars = no_tracers
      tr_ukca = no_tr_ukca
      num_optional_lbcs_out=0

! Add the number of optional lbc in this run depending on l_pc2
! and l_murk
      IF(l_mcr_qcf2)THEN
        num_optional_lbcs_out=num_optional_lbcs_out+1
      END IF
      IF(l_mcr_qrain)THEN
        num_optional_lbcs_out=num_optional_lbcs_out+1
      END IF
      IF(l_mcr_qgraup)THEN
        num_optional_lbcs_out=num_optional_lbcs_out+1
      END IF
      IF(l_pc2)THEN
        num_optional_lbcs_out=num_optional_lbcs_out+3
      END IF
      IF(l_murk)THEN
        num_optional_lbcs_out=num_optional_lbcs_out+1
      END IF
      
! Call a subroutine which sets up all aerosol relevant logicals
! and increments num_optional_lbcs_out as appropriate
! DEPENDS ON: makebc_lbc_logic
      CALL makebc_lbc_logic(num_optional_lbcs_out,ndustbin_in, ndustbin_out)
!
! calculate the number of LBcs
! = mandatory LBCs + optional LBCs + free tracerLBCs + UKCA LBCs
      intf_lookupsa = num_req_lbcs + num_optional_lbcs_out +            &
        tr_vars + tr_ukca

! DEPENDS ON: makebc
      CALL makebc(um_versn,                                             &
                  tracers_active,no_tracers,                            &
                  tr_ukca_active,no_tr_ukca, ndustbin_in, ndustbin_out)

 9999 CONTINUE

      CALL ereport_finalise( )

      CLOSE(5)

! DEPENDS ON: timer
      CALL Timer ( ProgName, 2 )


      WRITE (6,'(A)') ' '
      WRITE (6,'(A)') ' ##################################'
      WRITE (6,'(A)') ' MAKEBC program completed normally.'
      WRITE (6,'(A)') ' ##################################'
      END PROGRAM main_makebc
