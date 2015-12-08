! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Routine: SETTSCTL -------------------------------------------------
!
!    Purpose: Sets timestep loop control switches and STASHflags.
!             Note STEP on entry is the values at the N+1 (ie. updated)
!             timelevel; control switches are therefore set using
!             A_STEP/O_STEP, whereas physical switches are set using
!             A_STEP-1/O_STEP-1 to ensure correct synchronisation.
!             Note also that step on entry is N (i.e. not updated) when
!             called from INITIAL.
!
!
!    Programming standard: UM Doc Paper 3, version 8.2
!
!    Project task: C0
!
!    External documentation: On-line UM document C0 - The top-level
!                            control system
!
!     ------------------------------------------------------------------
!    Interface and arguments: ------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Top Level

SUBROUTINE settsctl (                                             &
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
! ARGINFA Headers for atmosphere interface data sets
        fixhd_intfa, inthd_intfa, lookup_intfa,                         &
        realhd_intfa,levdepc_intfa,                                     &
      ! Row/Col DEPC for variable resolution LBCs
        rowdepc_intfa, coldepc_intfa,                                   &  
      ! Eta values for LBC levels
        lbc_eta_theta, lbc_eta_rho,                                     &
! ARGINFA end
       internal_model,linitial,meanlev,icode,cmessage)

USE rad_input_mod, ONLY:                                          &
    A_SW_RadStep,                                                 &
    A_LW_RadStep,                                                 &
    A_SW_RadStep_diag,                                            &
    A_SW_Radstep_prog,                                            &
    A_LW_RadStep_diag,                                            &
    A_LW_Radstep_prog,                                            &
    l_sw_radiate,      l_lw_radiate,                              &
    l_sw_radiate_diag, l_lw_radiate_diag,                         &
    l_sw_radiate_prog, l_lw_radiate_prog

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IAU_mod, ONLY : L_IAU
USE UM_ParParams
USE Control_Max_Sizes
USE um_input_control_mod, ONLY: lcal360
USE eng_corr_inputs_mod, ONLY: a_energysteps, lenergy, lflux_reset

USE lookup_addresses
USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim,     &
     dumpfreqim, dumptimesim, dumptimes_len1,                     &
     meanfreqim, meanfreq_len1, jobrel_stepim, jobrel_len1,       &
     jobrel_offsetim   

USE Submodel_Mod     

USE nlstcall_mod, ONLY : model_basis_time, &
                         ancil_reftime, &
                         ft_steps, &
                         ft_firststep, &
                         ft_laststep, &
                         LPP, &
                         lpp_select, &
                         ldump, &
                         lmean, &
                         lprint, &
                         linterface, &
                         lexit, &
                         ljobrelease, &
                         lancillary, &
                         lboundary, &
                         lassimilation, &
                         model_assim_mode, &
                         type_letter_1, &
                         type_letter_2, &
                         ft_output, &
                         run_assim_mode

USE chsunits_mod, ONLY : nunits

IMPLICIT NONE

INTEGER internal_model   ! IN  - internal model identifier
LOGICAL linitial         ! IN  - true in called from INITIAL
INTEGER meanlev          ! OUT - Mean level indicator
! Needed for N_INTERNAL_MODEL
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
! TYPDUMA needs TYPSIZE included first
!L --------------- Dump headers (atmosphere)-------------
      INTEGER :: A_FIXHD(LEN_FIXHD)    ! fixed length header
      INTEGER :: A_INTHD(A_LEN_INTHD)  ! integer header
      INTEGER :: A_CFI1(A_LEN_CFI1+1)  ! compress field index
      INTEGER :: A_CFI2(A_LEN_CFI2+1)  ! compress field index
      INTEGER :: A_CFI3(A_LEN_CFI3+1)  ! compress field index

      REAL::A_REALHD(A_LEN_REALHD)                    ! real header
      REAL::A_LEVDEPC(A_LEN1_LEVDEPC*A_LEN2_LEVDEPC+1)! level  dep const
      REAL::A_ROWDEPC(A_LEN1_ROWDEPC*A_LEN2_ROWDEPC+1)! row    dep const
      REAL::A_COLDEPC(A_LEN1_COLDEPC*A_LEN2_COLDEPC+1)! column dep const
      REAL::A_FLDDEPC(A_LEN1_FLDDEPC*A_LEN2_FLDDEPC+1)! field  dep const
      REAL::A_EXTCNST(A_LEN_EXTCNST+1)                ! extra constants
      REAL::A_DUMPHIST(LEN_DUMPHIST+1)                ! temp hist file

      ! Meaningful parameter names for integer constants header:
! ----------------------- include file: IHEADAPM -----------------------
! Description: Meaningful parameter names to index A_INTHD array in
!              atmosphere dump, ie INTEGER CONSTANTS, and reduce magic
!              numbers in code.
!
      INTEGER,PARAMETER:: ih_a_step          = 1  ! Timestep no.
      INTEGER,PARAMETER:: ih_rowlength       = 6  ! No. of points E-W
      INTEGER,PARAMETER:: ih_rows            = 7  ! No. of points N-S

      ! No. of model levels (0=surface)
      INTEGER,PARAMETER:: ih_model_levels    = 8

      ! No. of model levels with moisture
      INTEGER,PARAMETER:: ih_wet_levels      = 9

      ! No. of deep soil temperature levels
      INTEGER,PARAMETER:: ih_soilT_levels    = 10

      INTEGER,PARAMETER:: ih_cloud_levels    = 11 ! No. of cloud levels
      INTEGER,PARAMETER:: ih_tracer_levels   = 12 ! No. of tracer levels

      ! No. of boundary layer levels
      INTEGER,PARAMETER:: ih_boundary_levels = 13
      INTEGER,PARAMETER:: ih_N_types         = 15 ! No. of field types

       ! Height generation method
      INTEGER,PARAMETER:: ih_height_gen      = 17

      ! First rho level at which height is constant
      INTEGER,PARAMETER:: ih_1_c_rho_level   = 24

      INTEGER,PARAMETER:: ih_land_points     = 25 ! No. of land points
      INTEGER,PARAMETER:: ih_ozone_levels    = 26 ! No. of ozone levels

      ! No. of deep soil moisture levels
      INTEGER,PARAMETER:: ih_soilQ_levels    = 28

      ! Number of convective cloud levels
      INTEGER,PARAMETER:: ih_convect_levels  = 34
      INTEGER,PARAMETER:: ih_rad_step        = 35 ! Radiation timestep
      INTEGER,PARAMETER:: ih_AMIP_flag       = 36 ! Flag for AMIP run
      INTEGER,PARAMETER:: ih_AMIP_year       = 37 ! First AMIP year
      INTEGER,PARAMETER:: ih_AMIP_month      = 38 ! First AMIP month
      INTEGER,PARAMETER:: ih_AMIP_day        = 49 ! First AMIP day
      INTEGER,PARAMETER:: ih_ozone_month     = 40 ! Current ozone month
      INTEGER,PARAMETER:: ih_SH_zonal_quad   = 41 ! L_SH_zonal_quadratics
      INTEGER,PARAMETER:: ih_SH_zonal_begin  = 42 ! SH_zonal_begin
      INTEGER,PARAMETER:: ih_SH_zonal_period = 43 ! SH_zonal_period
      INTEGER,PARAMETER:: ih_SH_level_weight = 44 ! SuHe_level_weight
      INTEGER,PARAMETER:: ih_SH_sigma_cutoff = 45 ! SuHe_sigma_cutoff
      INTEGER,PARAMETER:: ih_friction_time   = 46 ! frictional_timescale

! IHEADAPM end
      ! Meaningful parameter names for real constants header:
! ----------------------- include file: RHEADAPM -----------------------
! Description: Meaningful parameter names to index A_REALHD array in
!              atmosphere dump, ie REAL CONSTANTS, and reduce magic
!              numbers in code.

      ! East-West   grid spacing in degrees
      INTEGER,PARAMETER:: rh_deltaEW         = 1

      ! North-South grid spacing in degrees
      INTEGER,PARAMETER:: rh_deltaNS         = 2

      ! Latitude  of first p point in degrees
      INTEGER,PARAMETER:: rh_baselat         = 3

      ! Longitude of first p point in degrees
      INTEGER,PARAMETER:: rh_baselong        = 4

      ! Latitude  of rotated N pole in degrees
      INTEGER,PARAMETER:: rh_rotlat          = 5

      ! Longitude of rotated N pole in degrees
      INTEGER,PARAMETER:: rh_rotlong         = 6

      ! Height of top theta level (m)
      INTEGER,PARAMETER:: rh_z_top_theta     =16

      ! total moisture of the atmosphere
      INTEGER,PARAMETER:: rh_tot_m_init      =18

      ! total mass of atmosphere
      INTEGER,PARAMETER:: rh_tot_mass_init   =19

      ! total energy of atmosphere
      INTEGER,PARAMETER:: rh_tot_energy_init =20

      ! energy correction = energy drift
      INTEGER,PARAMETER:: rh_energy_corr     =21

! RHEADAPM end
      ! Meaningful parameter names for fixed header:
! ----------------------- include file: FHEADAPM -----------------------
! Description: Meaningful parameter names to index A_FIXHD array in
!              atmosphere dump, ie REAL CONSTANTS, and reduce magic
!              numbers in code.
 
      ! Start of Row Dependent Constant
      INTEGER,PARAMETER:: fh_RowDepCStart   = 115

      ! Start of Col Dependent Constant
      INTEGER,PARAMETER:: fh_ColDepCStart   = 120

! FHEADAPM end
      ! PP headers

      INTEGER :: A_LOOKUP(LEN1_LOOKUP,A_LEN2_LOOKUP) ! lookup heads
      INTEGER :: A_MPP_LOOKUP(MPP_LEN1_LOOKUP,A_LEN2_LOOKUP)
      INTEGER :: a_ixsts(len_a_ixsts)     ! stash index array

      REAL    :: a_spsts(len_a_spsts)     ! atmos stash array
! TYPDUMA end
! Contains *CALL CPPXREF
! TYPSTS starts
! submodel_mod must be included before this file
!Applicable to all configurations
!STASH related variables for describing output requests and space
!management.
! Include sizes for dimensioning arrays in this deck
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

! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! Whether a calculation is needed for SF above
      LOGICAL :: SF_CALC(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
! TYPSTS end
! TYPINFA start
!
! This file needs file TYPSIZE to be included

      ! Headers for atmosphere interface data sets
      INTEGER::FIXHD_INTFA(LEN_FIXHD,N_INTF_A)        ! Fixed header
      INTEGER::INTHD_INTFA(PP_LEN_INTHD,N_INTF_A)     ! Integer header
      INTEGER::LOOKUP_INTFA(LEN1_LOOKUP,INTF_LOOKUPSA,N_INTF_A)! Lookups

      REAL :: REALHD_INTFA(PP_LEN_REALHD,N_INTF_A)   ! Real header
      REAL :: LEVDEPC_INTFA(MAX_INTF_MODEL_LEVELS,INTF_LEN2_LEVDEPC,    &
     &  N_INTF_A)
      REAL::ROWDEPC_INTFA(MAX_LBCROWS,INTF_LEN2_ROWDEPC,N_INTF_A)
      REAL::COLDEPC_INTFA(MAX_LBCROW_LENGTH,INTF_LEN2_COLDEPC,N_INTF_A)
 
      ! Eta Values for each area
      REAL :: LBC_ETA_THETA (MAX_INTF_MODEL_LEVELS+1,N_INTF_A)
      REAL :: LBC_ETA_RHO   (MAX_INTF_MODEL_LEVELS  ,N_INTF_A)
! TYPINFA end
INTEGER icode            ! Out - Return code
CHARACTER(LEN=80) cmessage  ! Out - Return error message

!*----------------------------------------------------------------------
!  Common blocks

! --------------------- Comdeck: CHISTORY ----------------------------
!
!  Purpose: COMMON block for history data needed by top level (C0)
!           routines, and passed from run to run.  Mostly set by
!           the User Interface.
!
!           Note that CHISTORY *CALLs ALL individual history comdecks
!
! --------------------------------------------------------------------
!
! ----------------------- Comdeck: IHISTO   ----------------------------
! Description: COMDECK defining Integer History variables for the
!              model overall.
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

!   Type declarations


      ! Array containing model data time (Same as MODEL_BASIS_TIME/MODEL
      ! ANALYSIS_HRS depending whether before/after assimilation)
      INTEGER :: model_data_time(6)

      ! Indicator of operational run type
      INTEGER :: run_indic_op

      ! Final target date for the run
      INTEGER :: run_resubmit_target(6)

      ! Last field written/read per FT unit
      INTEGER :: ft_lastfield(20:nunits)

      ! Number of automatically-resubmitted job chunks
      ! Used to name output file
      INTEGER :: run_job_counter

! History Common Block for overall model integers variables.

      COMMON /ihisto/                                                 &
         model_data_time,                                             &
         run_indic_op, run_job_counter,                               &
         run_resubmit_target, ft_lastfield

      NAMELIST /nlihisto/                                             &
         model_data_time,                                             &
         run_indic_op, run_job_counter,                               &
         run_resubmit_target, ft_lastfield

! IHISTO end
! ----------------------- Comdeck: CHISTO   ----------------------------
! Description: COMDECK defining Character History variables for the
!              model overall.
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

  CHARACTER(LEN=8) ::  run_type             ! Type of run
  CHARACTER(LEN=1) ::  ft_active(20:nunits) ! "Y" if file partly written

  LOGICAL :: newrun ! Set to true in NRUN to stop auto-resubmission

  ! History Common Block for overall model character variables.

  COMMON /chisto/                                     &
     run_type,                                        &
     newrun, ft_active

  NAMELIST /nlchisto/                                 &
     run_type,                                        &
     ft_active

! CHISTO end
! ----------------------- Comdeck: IHISTG   ----------------------------
! Description: COMDECK defining Integer History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! This file belongs in section: Top Level
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

!   Type declarations
      ! History block copy of A_STEP held in file CTIME
      INTEGER :: h_stepim(n_internal_model_max)

      ! No of means activated
      INTEGER :: mean_offsetim(n_internal_model_max)

      ! Offset between MEAN_REFTIME and model basis time(in model dumps)
      INTEGER :: offset_dumpsim(n_internal_model_max)

      ! No of mean periods chosen
      INTEGER :: mean_numberim(n_internal_model_max)

      ! Indicators used to correct logical units are used for
      ! atmos partial sum dump I/O
      INTEGER :: run_meanctl_indicim(4,n_internal_model_max)

      ! History Common Block for generic model integer variables.

      COMMON /ihistg/                                         &
         h_stepim, mean_offsetim, offset_dumpsim,             &
         mean_numberim, run_meanctl_indicim

      NAMELIST /nlihistg/                                     &
         h_stepim, mean_offsetim, offset_dumpsim,             &
         run_meanctl_indicim

! IHISTG end
! ----------------------- Comdeck: CHISTG   ----------------------------
! Description: COMDECK defining Character variables for
!              managing dump names
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 
!
!   Type declarations
!
! For keeping old restart dump name between creation of new restart dump
! and successful completion of climate means and history file.
CHARACTER(LEN=256) :: save_dumpname_im(n_internal_model_max)
! Name of current restart dump
CHARACTER(LEN=256) :: checkpoint_dump_im(n_internal_model_max)
! Blank name
CHARACTER(LEN=256) :: blank_file_name
!
! History Common Block for generic model characters variables.
!
COMMON /chistg/save_dumpname_im, checkpoint_dump_im, blank_file_name

NAMELIST /nlchistg/checkpoint_dump_im

! CHISTG end
!
!  Purpose: Defines unit numbers relevant to history file
!           and variables used to hold the logical to physical
!           file associations made within the model
!
!  Logical Filenames used in the model
!
      CHARACTER(LEN=256) hkfile,ppxref,config,stashctl,namelist,output,      &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,ftxx,    &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar
     

!
      CHARACTER(LEN=256) MODEL_FT_UNIT ! Array holding FORTRAN unit file
!                                 ! associations details for each unit
!
      INTEGER                                                           &
     &        MCTL_UNIT,                                                &
                                 ! Master control namelist file unit
     &        ICTL_UNIT,                                                &
                                 ! Interim control namelist file unit
     &        RSUB_UNIT,                                                &
                                 ! File indicating whether resub required
     &        XHIST_UNIT,                                               &
                                 ! Main history file unit
     &        THIST_UNIT,                                               &
                                 ! Backup history file unit
     &        HKFILE_UNIT,                                              &
                                 ! Operational houskeeping file unit    
     &        EG_UNIT            ! ENDGame diagnostics/info unit
!
! Parameters specifying unit numbers relevant to control/history tasks
!
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(RSUB_UNIT =10)
      PARAMETER(XHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)

!
! Parameters specifying unit numbers relevant to ENDGame diagnostics
!
      PARAMETER(EG_UNIT  = 55)

! UKCA unit numbers

      INTEGER, PARAMETER :: ukcafjxx_unit=170 ! Fast-J(X) cross section data
      INTEGER, PARAMETER :: ukcafjsc_unit=171 ! Fast-JX scattering data
      INTEGER, PARAMETER :: ukca2do3_unit=172 ! 2D top boundary O3 data 
      INTEGER, PARAMETER :: ukca2ch4_unit=173 ! 2D top boundary CH4 data
      INTEGER, PARAMETER :: ukca2noy_unit=174 ! 2D top boundary NOY data
      INTEGER, PARAMETER :: ukca2pho_unit=175 ! 2D photolysis input data
      INTEGER, PARAMETER :: ukcastrd_unit=176 ! Stratospheric model radiation field. 
      INTEGER, PARAMETER :: ukcasto3_unit=177 ! Strat standard atmosphere T and O3.
      INTEGER, PARAMETER :: ukcastar_unit=178 ! Stratospheric sulfate aerosol climatology 
      INTEGER, PARAMETER :: ukcafjar_unit=179 ! Sulfate aerosol cliamtology for Fast-JX
! Text output file for STASH-related information is assigned to UNIT 200

!
! Namelist of all permissible logical files.
!
      NAMELIST / nlcfiles /                                             &
                   hkfile,ppxref,config,stashctl,namelist,output,       &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,         &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar

!
!Common block definition
!
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
!
! Equivalence logical filenames within array MODEL_FT_UNIT
!
      EQUIVALENCE                                                       &
     &(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ), &
     &(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ), &
     &(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ), &
     &(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ), &
     &(ICTL      ,MODEL_FT_UNIT(9)  ),(RSUB       ,MODEL_FT_UNIT(10) ), &
     &(XHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ), &
     &(ICECALVE  ,MODEL_FT_UNIT(13) ),                                  &
     &(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ), &
     &                                (ASWAP      ,MODEL_FT_UNIT(18) ), &
     &(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ), &
     &(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ), &
     &(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ), &
     &(AOPSUM3   ,MODEL_FT_UNIT(25) )
!
      EQUIVALENCE                                                       &
     &(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ), &
     &(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ), &
     &(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ), &
     &(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ), &
     &(GENLAND   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ), &
     &(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ), &
     &(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ), &
     &(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ), &
     &(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ), &
     &(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ), &
     &(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ), &
     &(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ), &
     &(WFIN      ,MODEL_FT_UNIT(50) )
!
      EQUIVALENCE                                                       &
     &(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ), &
     &(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ), &
     &                                (FLUXCORR   ,MODEL_FT_UNIT(56) ), &
     &(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), &
     &(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ), &
     &(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ), &
     &(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ), &
     &(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ), &
     &(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ), &
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ), &
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ), &
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) ), &
     &(DUSTSOIL  ,MODEL_FT_UNIT(75) ),(BIOMASS    ,MODEL_FT_UNIT(76) ), &
     &(RIVSTOR   ,MODEL_FT_UNIT(77) ),(RIVCHAN    ,MODEL_FT_UNIT(78) ), &
     &(RIVER2A   ,MODEL_FT_UNIT(79) )
!
      EQUIVALENCE                                                       &
                                      (lwspectd   ,model_ft_unit(80) ), &
      (surgeou1  ,model_ft_unit(81) ),(surgeout   ,model_ft_unit(82) ), &
      (ppscreen  ,model_ft_unit(83) ),(ppsmc      ,model_ft_unit(84) ), &
      (wfout     ,model_ft_unit(85) ),(uarsout1   ,model_ft_unit(86) ), &
      (uarsout2  ,model_ft_unit(87) ),(icefout    ,model_ft_unit(88) ), &
      (mosout    ,model_ft_unit(89) ),(vert_lev   ,model_ft_unit(90) ), &
      (sstout    ,model_ft_unit(91) ),(siceout    ,model_ft_unit(92) ), &
      (curntout  ,model_ft_unit(93) ),(flxcrout   ,model_ft_unit(94) ), &
      (dmsconc   ,model_ft_unit(95) ),(orog       ,model_ft_unit(96) ), &
      (transp    ,model_ft_unit(97) ),(olabcin    ,model_ft_unit(98) ), &
      (ocndepth  ,model_ft_unit(99) ),                                  &
      (foamout1  ,model_ft_unit(100)),(foamout2   ,model_ft_unit(101)), &
      (cxbkgerr  ,model_ft_unit(102)),(rfmout     ,model_ft_unit(103)), &
      (idealise  ,model_ft_unit(106)),(tdf_dump   ,model_ft_unit(107)), &
      (iau_inc   ,model_ft_unit(108)),(murkfile   ,model_ft_unit(109)), &
      (sulpemis  ,model_ft_unit(110)),(usrancil   ,model_ft_unit(111)), &
      (usrmulti  ,model_ft_unit(112)),(ousrancl   ,model_ft_unit(113)), &
      (ousrmult  ,model_ft_unit(114)),(so2natem   ,model_ft_unit(115)), &
      (chemoxid  ,model_ft_unit(116)),(aerofcg    ,model_ft_unit(117)), &
      (co2emits  ,model_ft_unit(118)),(tppsozon   ,model_ft_unit(119)), &
      (landfrac  ,model_ft_unit(120)),(wlabcou1   ,model_ft_unit(121)), &
      (wlabcou2  ,model_ft_unit(122)),(wlabcou3   ,model_ft_unit(123)), &
      (wlabcou4  ,model_ft_unit(124)),(alabcin1   ,model_ft_unit(125)), &
      (alabcin2  ,model_ft_unit(126)),                                  &
      (ocffemis  ,model_ft_unit(128)),(horzgrid   ,model_ft_unit(129)), &
      (surfemis  ,model_ft_unit(130)),(aircrems   ,model_ft_unit(131)), &
      (stratems  ,model_ft_unit(132)),(extraems   ,model_ft_unit(133)), &
      (radonems  ,model_ft_unit(134)),(fracinit   ,model_ft_unit(135)), &
      (veginit   ,model_ft_unit(136)),(disturb    ,model_ft_unit(137)), &
      (cached    ,model_ft_unit(138)),(sootemis   ,model_ft_unit(139)), &
      (alabcou1  ,model_ft_unit(140)),(alabcou2   ,model_ft_unit(141)), &
      (alabcou3  ,model_ft_unit(142)),(alabcou4   ,model_ft_unit(143)), &
      (alabcou5  ,model_ft_unit(144)),(alabcou6   ,model_ft_unit(145)), &
      (alabcou7  ,model_ft_unit(146)),(alabcou8   ,model_ft_unit(147)), &
      (cariolo3  ,model_ft_unit(148)),(rpseed     ,model_ft_unit(149)), &
      (ppvar     ,model_ft_unit(150)),(pp10       ,model_ft_unit(151)), &
      (icfile    ,model_ft_unit(152)),(var_grid   ,model_ft_unit(153)), &
      (arclbiog  ,model_ft_unit(154)),(arclbiom   ,model_ft_unit(155)), &
      (arclblck  ,model_ft_unit(156)),(arclsslt   ,model_ft_unit(157)), &
      (arclsulp  ,model_ft_unit(158)),(arcldust   ,model_ft_unit(159)), &
      (arclocff  ,model_ft_unit(160)),(arcldlta   ,model_ft_unit(161)), &
      (topmean   ,model_ft_unit(162)),(topstdev   ,model_ft_unit(163)), &
      (ppmbc     ,model_ft_unit(164)),(ukcaprec   ,model_ft_unit(165)), &
      (ukcaacsw  ,model_ft_unit(166)),(ukcaaclw   ,model_ft_unit(167)), &
      (ukcacrsw  ,model_ft_unit(168)),(ukcacrlw   ,model_ft_unit(169)), &
      (ukcafjxx  ,model_ft_unit(170)),(ukcafjsc   ,model_ft_unit(171)), &
      (ukca2do3  ,model_ft_unit(172)),(ukca2ch4   ,model_ft_unit(173)), &
      (ukca2noy  ,model_ft_unit(174)),(ukca2pho   ,model_ft_unit(175)), &
      (ukcastrd  ,model_ft_unit(176)),(ukcasto3   ,model_ft_unit(177)), &
      (ukcastar  ,model_ft_unit(178)),(ukcafjar   ,model_ft_unit(179))
! Text output file for STASH-related information is assigned to UNIT 200

!
!LL  Comdeck: CCONTROL -------------------------------------------------
!LL
!LL  Purpose: COMMON block for top level switches and 2nd level switches
!LL           needed by the top level (C0) and 2nd level routines, but
!LL           not held in the history COMMON block.
!LL
!LLEND ---------------------------------------------------------------

!#include "cntlall.h"
! cntlgen.h was replaced by control/top_level/nlstgen_mod.F90
! #include "cntlgen.h"

! CTIME ----------------------------------------------------
!
!  Purpose: Derived model time/step information including start/end
!           step numbers and frequencies (in steps) of interface field
!           generation, boundary field updating, ancillary field
!           updating; and assimilation start/end times.
!           NB: Last three are set by IN_BOUND, INANCCTL, IN_ACCTL.
!           Also contains current time/date information, current
!           step number (echoed in history file) and steps-per-group.
!
!END -----------------------------------------------------------------

      INTEGER :: I_YEAR               ! Current model time (years)
      INTEGER :: I_MONTH              ! Current model time (months)
      INTEGER :: I_DAY                ! Current model time (days)
      INTEGER :: I_HOUR               ! Current model time (hours)
      INTEGER :: I_MINUTE             ! Current model time (minutes)
      INTEGER :: I_SECOND             ! Current model time (seconds)
      INTEGER :: I_DAY_NUMBER         ! Current model time (day no)
      INTEGER :: PREVIOUS_TIME(7)     ! Model time at previous step
      INTEGER :: IAU_DTResetStep      ! Data time reset step for IAU run

      INTEGER :: BASIS_TIME_DAYS  ! Integral no of days to basis time
      INTEGER :: BASIS_TIME_SECS  ! No of seconds-in-day at basis time

      LOGICAL :: L_C360DY

! UM6.5MODEL_ANALYSIS_HRS changed to REAL - 
!   requires FORECAST_HRS and DATA_MINUS_BASIS_HRS to REAL also 
      REAL    :: FORECAST_HRS     ! Hours since Data Time (ie T+nn)
      REAL    :: DATA_MINUS_BASIS_HRS ! Data time - basis time (hours)

      COMMON /CTIMED/ I_YEAR,I_MONTH,I_DAY,I_HOUR,I_MINUTE,I_SECOND,    &
        I_DAY_NUMBER,PREVIOUS_TIME,                                     &
        BASIS_TIME_DAYS,BASIS_TIME_SECS,                                &
        FORECAST_HRS,DATA_MINUS_BASIS_HRS,                              &
        IAU_DTResetStep, L_C360DY

      INTEGER :: STEPim(INTERNAL_ID_MAX)  ! Step no since basis time
      INTEGER :: GROUPim(INTERNAL_ID_MAX) ! Number of steps per group

      ! Finish step number this run
      INTEGER :: TARGET_END_STEPim(INTERNAL_ID_MAX)

      REAL :: SECS_PER_STEPim(INTERNAL_ID_MAX) ! Timestep length in secs

      ! Frequency of interface field generation in steps
      INTEGER :: INTERFACE_STEPSim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Start steps for interface field generation
      INTEGER :: INTERFACE_FSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! End steps for interface field generation
      INTEGER :: INTERFACE_LSTEPim(MAX_N_INTF_A,INTERNAL_ID_MAX)

      ! Frequency of  updating boundary fields in steps
      INTEGER :: BOUNDARY_STEPSim(INTERNAL_ID_MAX)

      ! No of steps from boundary data prior to basis time to model
      ! basis time
      INTEGER :: BNDARY_OFFSETim(INTERNAL_ID_MAX)

      ! Lowest frequency for updating of ancillary fields in steps
      INTEGER :: ANCILLARY_STEPSim(INTERNAL_ID_MAX)

      ! Start steps for assimilation
      INTEGER :: ASSIM_FIRSTSTEPim(INTERNAL_ID_MAX)

      ! Number of assimilation steps to analysis
      INTEGER :: ASSIM_STEPSim(INTERNAL_ID_MAX)

      ! Number of assimilation steps after analysis
      INTEGER :: ASSIM_EXTRASTEPSim(INTERNAL_ID_MAX)

      COMMON/CTIMEE/                                                    &
     &  STEPim,GROUPim,TARGET_END_STEPim,INTERFACE_STEPSim,             &
     &  INTERFACE_FSTEPim,INTERFACE_LSTEPim,BOUNDARY_STEPSim,           &
     &  BNDARY_OFFSETim,ANCILLARY_STEPSim,ASSIM_FIRSTSTEPim,            &
     &  ASSIM_STEPSim,ASSIM_EXTRASTEPSim,SECS_PER_STEPim

! CTIME end
! STPARAM
!
!  Purpose: Meaningful PARAMETER names for STASH processing routines.
!           Both a long name and short name have been declared, to
!           reduce the existence of "magic" numbers in STASH.
!           Format is that first the address of the item is declare in
!           both long and short form. example is;
!             integer st_item_code,s_item  !Item number (declaration)
!             parameter(st_item_code=3,s_item=3)
!
!  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!  Logical components covered: D70
!
!  Project task: D7
!
!  External documentation:
!    Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                 system (STASH)
!--------------------------------------------------------------

      ! Internal model number address
      INTEGER,PARAMETER:: st_model_code = 28
      INTEGER,PARAMETER:: s_modl        = 28

      ! Section Number address
      INTEGER,PARAMETER:: st_sect_no_code = 2
      INTEGER,PARAMETER:: s_sect          = 2
      INTEGER,PARAMETER:: st_sect_code    = 2

      INTEGER,PARAMETER:: st_item_code=1,s_item=1 ! Item number address

      ! Processing Code address
      INTEGER,PARAMETER:: st_proc_no_code=3,s_proc=3

      ! subsidiary codes for st_proc_no_code now
      INTEGER,PARAMETER:: st_replace_code=1
      INTEGER,PARAMETER:: st_accum_code=2
      INTEGER,PARAMETER:: st_time_mean_code=3
      INTEGER,PARAMETER:: st_time_series_code=4
      INTEGER,PARAMETER:: st_max_code=5
      INTEGER,PARAMETER:: st_min_code=6
      INTEGER,PARAMETER:: st_append_traj_code=7
      INTEGER,PARAMETER:: st_time_series_mean=8
      INTEGER,PARAMETER:: st_variance_code=9

      ! Frequency (Input & output) addres
      INTEGER,PARAMETER:: st_freq_code=4,s_freq=4

      ! Offset for sampling
      INTEGER,PARAMETER:: st_offset_code=30,s_offs=30

      ! start timestep address
      INTEGER,PARAMETER:: st_start_time_code=5,s_times=5

      ! end timestep address
      INTEGER,PARAMETER:: st_end_time_code=6,s_timee=6

      ! period in timesteps address
      INTEGER,PARAMETER:: st_period_code=7,s_period=7

      ! infinite end/period value
      INTEGER,PARAMETER:: st_infinite_time=-1

      INTEGER,PARAMETER:: st_end_of_list=-1 !end-of-list marker in times

      ! grid point stuff
      ! gridpoint info address
      INTEGER,PARAMETER:: st_gridpoint_code=8,s_grid=8

      ! now subsid grid point stuff
      ! no masking done
      INTEGER,PARAMETER:: stash_null_mask_code=1,s_nomask=1

      ! land mask conds
      INTEGER,PARAMETER:: stash_land_mask_code=2,s_lndms=2

      ! sea mask code
      INTEGER,PARAMETER:: stash_sea_mask_code=3,s_seams =3

      ! processing options

      ! size of block for gridpoint code
      INTEGER,PARAMETER:: block_size=10

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: extract_base=block_size*0

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: extract_top=block_size*1

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_base=block_size*1

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_top=block_size*2

      ! max code for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_base=block_size*2

      ! base codes for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_top=block_size*3

      ! max code for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_base=block_size*3

      ! base codes for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_top=block_size*4

      ! max code for field mean subroutine
      INTEGER,PARAMETER:: field_mean_base=block_size*4

      ! base codes for field mean subroutine
      INTEGER,PARAMETER:: field_mean_top=block_size*5

      ! max code for global mean subroutine
      INTEGER,PARAMETER:: global_mean_base=block_size*5

      ! base codes for global mean subroutine
      INTEGER,PARAMETER:: global_mean_top=block_size*6

      ! Weighting

      ! weighting info address
      INTEGER,PARAMETER:: st_weight_code=9,s_weight=9

      INTEGER,PARAMETER:: stash_weight_null_code  =0,s_noweight  =0
      INTEGER,PARAMETER:: stash_weight_area_code  =1,s_areaweight=1
      INTEGER,PARAMETER:: stash_weight_volume_code=2,s_volweight =2
      INTEGER,PARAMETER:: stash_weight_mass_code  =3,s_massweight=3

      ! Domain definition

      ! row addresses
      INTEGER,PARAMETER:: st_north_code=12,s_north=12
      INTEGER,PARAMETER:: st_south_code=13,s_south=13
      INTEGER,PARAMETER:: st_west_code =14,s_west =14
      INTEGER,PARAMETER:: st_east_code =15,s_east =15

      ! Levels

      ! input bottom level address
      INTEGER,PARAMETER:: st_input_bottom=10,s_bottom =10

      ! special code
      INTEGER,PARAMETER:: st_special_code=100,s_special=100

      ! input top level address
      INTEGER,PARAMETER:: st_input_top=11,s_top=11

      ! output bottom level address
      INTEGER,PARAMETER:: st_output_bottom=21,s_outbot=21

      ! output top level address
      INTEGER,PARAMETER:: st_output_top=22,s_outtop=22

      INTEGER,PARAMETER:: st_model_level_code=1,s_model=1

      ! code for pressure leve
      INTEGER,PARAMETER:: st_pressure_level_code=2,s_press=2

      ! code for height levels
      INTEGER,PARAMETER:: st_height_level_code=3,s_height=3

      ! input code addres
      INTEGER,PARAMETER:: st_input_code=16,s_input=16

      ! input length of diagnostic address
      INTEGER,PARAMETER:: st_input_length=17,s_length=17

      ! output code address
      INTEGER,PARAMETER:: st_output_code=18,s_output=18

      ! Pointer to D1 addressing information
      ! Pos of item in D1 for relevant submodel
      INTEGER,PARAMETER:: st_position_in_d1=29,st_d1pos=29

      ! Output destination options

      INTEGER,PARAMETER:: st_dump=1
      INTEGER,PARAMETER:: st_secondary=2

      ! output length of diagnostic address
      INTEGER,PARAMETER:: st_output_length=19,s_outlen=19
         integer st_dump_output_length,s_doutlen ! output length on
         parameter(st_dump_output_length=32,s_doutlen=32)  ! dump
         integer st_dump_level_output_length,s_dlevoutlen
         parameter(st_dump_level_output_length=33,s_dlevoutlen=33)
! output length of a single level on dump

         integer st_output_addr,s_outadd ! start locn of diag after stas
         parameter(st_output_addr=20,s_outadd=20)       ! output address
         integer st_dump_output_addr,s_doutadd ! output address on
         parameter(st_dump_output_addr=31,s_doutadd=31)  ! dump

      ! ptr to dump lookup header address
      INTEGER,PARAMETER:: st_lookup_ptr=23

      ! ptr into stash_series where control data address
      INTEGER,PARAMETER:: st_series_ptr=24

      ! subsid stuff for time series
      INTEGER,PARAMETER:: series_grid_type=1
      INTEGER,PARAMETER:: series_grid_code=0
      INTEGER,PARAMETER:: series_long_code=1
      INTEGER,PARAMETER:: series_size=2
      INTEGER,PARAMETER:: series_proc_code=3
      INTEGER,PARAMETER:: series_north=4
      INTEGER,PARAMETER:: series_south=5
      INTEGER,PARAMETER:: series_west=6
      INTEGER,PARAMETER:: series_east=7
      INTEGER,PARAMETER:: series_list_start=8
      INTEGER,PARAMETER:: series_list_end=9
      INTEGER,PARAMETER:: record_size=9

      ! Miscellaneous parameters

      ! system/user tag field in stlist address
      INTEGER,PARAMETER:: st_macrotag=25

      ! Pseudo-level list pointers

      ! pseudo-levels input list address
      INTEGER,PARAMETER:: st_pseudo_in=26

      ! pseudo-levels output list address
      INTEGER,PARAMETER:: st_pseudo_out=27

      ! Internal horizontal gridtype codes common to all diagnostics

      INTEGER,PARAMETER:: st_tp_grid =1 ! T-p grid
      INTEGER,PARAMETER:: st_uv_grid =2 ! u-v grid
      INTEGER,PARAMETER:: st_cu_grid =3 ! C-grid u point
      INTEGER,PARAMETER:: st_cv_grid =4 ! C-grid v point
      INTEGER,PARAMETER:: st_zt_grid =5 ! Zonal T-grid
      INTEGER,PARAMETER:: st_zu_grid =6 ! Zonal u-grid
      INTEGER,PARAMETER:: st_mt_grid =7 ! Meridional T-grid
      INTEGER,PARAMETER:: st_mu_grid =8 ! Meridional u-grid
      INTEGER,PARAMETER:: st_riv_grid= 23    ! river_routing grid
      INTEGER,PARAMETER:: st_scalar  =9 ! Scalar (ie. single value)
      INTEGER,PARAMETER:: st_wam_all= 60    ! Wam Field on Full Grid
      INTEGER,PARAMETER:: st_wam_sea= 62    ! Wam Field on Sea Points

! STPARAM end

!  Local variables

INTEGER i,nftunit        ! Loop counters
INTEGER itime            ! Loop index over selected dump times
INTEGER step             ! A_STEP or O_STEP for atmos/ocean
INTEGER is,ii,il,im,ntab,it,ie  ! STASH variables
INTEGER modl             ! Int model no, read from STASH list arra
INTEGER modjobrel        ! modulus of JOBREL_STEP()
INTEGER step_assim       ! model step rel. to assim start

INTEGER jintf            ! Interface area index
INTEGER ancil_ref_days,ancil_ref_secs
INTEGER ancil_offset_steps ! offset of ref. from basis time
INTEGER secs_per_step    ! seconds per timestep
INTEGER months_in        ! Number of months into forecast
INTEGER secs_per_period
INTEGER steps_per_period
INTEGER dumpfreq
INTEGER offset_dumps
INTEGER exitfreq
INTEGER target_end_step
INTEGER ancillary_steps
INTEGER boundary_steps
INTEGER bndary_offset
INTEGER dumptimes(dumptimes_len1)
INTEGER meanfreq (meanfreq_len1)
INTEGER jobrel_step(jobrel_len1)

LOGICAL iau_resetdt ! If .TRUE., reset data time.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('SETTSCTL',zhook_in,zhook_handle)

icode=0

! Initialise internal_model
internal_model = 1 ! Atmosphere

!  1. Set timestep loop top-level control switches
!
!  1.0  Initialise control switches which are shared in coupled runs
!

step=                        stepim(internal_model)
secs_per_period=  secs_per_periodim(internal_model)
steps_per_period=steps_per_periodim(internal_model)
secs_per_step=secs_per_period/steps_per_period
dumpfreq=                dumpfreqim(internal_model)
offset_dumps=        offset_dumpsim(internal_model)
target_end_step=  target_end_stepim(internal_model)
ancillary_steps=  ancillary_stepsim(internal_model)
boundary_steps =   boundary_stepsim(internal_model)
bndary_offset  =   bndary_offsetim(internal_model)

exitfreq=dumpfreq

DO i=1,dumptimes_len1
  dumptimes(i)=       dumptimesim(i,internal_model)
END DO ! i
DO i=1,meanfreq_len1
  meanfreq(i)=         meanfreqim(i,internal_model)
END DO ! i
DO i=1,jobrel_len1
  jobrel_step(i)=   jobrel_stepim(i,internal_model)
END DO ! i
meanlev=0

lassimilation=.FALSE.
ldump=        .FALSE.
lexit=        .FALSE.
lmean=        .FALSE.
lprint=       .FALSE.
lancillary=   .FALSE.
lboundary=    .FALSE.
linterface=   .FALSE.
ljobrelease = .FALSE.
!
!  1.1  Set up PPfile switches for the timestep
!
lpp =.FALSE.
  DO nftunit=20,nunits
    IF (ft_output(nftunit) == 'Y') THEN

      ! Initialise:
      lpp_select(nftunit) = .FALSE.

      ! Allow mix of real-month & regular reinit. periods on
      ! different units:
      IF (ft_steps(nftunit) <  0) THEN ! Real-month reinit.

        ! Select files to be reinitialised on this timestep:
        ! First initialisation may be on any timestep, not just 0
        IF ((type_letter_1(nftunit) == 'c' .AND.                   &
             step-1 == ft_firststep(nftunit))                      &
        .OR.(type_letter_1(nftunit) /= 'c' .AND.                   &
             step == ft_firststep(nftunit))) THEN

          lpp_select(nftunit) =.TRUE.

        ELSE IF (.NOT.lcal360                                                 &
                 .AND.(nftunit >= 60 .AND. nftunit <  70 .OR. nftunit == 151) &
                 .AND.(i_day     ==  1)                                       &
                 .AND.(i_hour    ==  secs_per_step/3600           )           &
                 .AND.(i_minute  ==  (MOD(secs_per_step,3600))/60 )           &
                 .AND.(i_second  ==   MOD(secs_per_step,60)       )           &
                 .AND.(step /= 1)) THEN

          ! Gregorian calendar: 
          ! Reinitialisation can only be selected when day=1 and
          ! hours, minutes and seconds equals exactly one time-step
          ! This code relies on the rule that the time-step length
          ! divides into a day.

          ! Additionally initialisation only takes place on 1-month,
          ! 3-month (season) or 12-month boundaries. So calculate months
          ! since start of job, and request reinitialisation at the 
          ! appropriate time according to which option selected.
          months_in = i_month - model_basis_time(2) +             &
             12 * (i_year - model_basis_time(1))

          IF (months_in >= ft_firststep(nftunit) .AND.             &
             (ft_laststep(nftunit) <  0 .OR.                       &
              months_in <= ft_laststep(nftunit))) THEN

            IF (ft_steps(nftunit) == -1) THEN
              lpp_select(nftunit) = .TRUE. ! Months
            ELSE IF (ft_steps(nftunit) == -3 .AND.                 &
                     MOD(months_in -                               &
                      ft_firststep(nftunit),3)  == 0) THEN
              lpp_select(nftunit) = .TRUE. ! Seasons
            ELSE IF (ft_steps(nftunit) == -12 .AND.                &
                     MOD(months_in -                               &
                      ft_firststep(nftunit),12) == 0) THEN
              lpp_select(nftunit) = .TRUE. ! Years
            END IF ! Of ft_steps(NFTUNIT) = reinit period

          END IF ! Of MONTHS_IN within reinitialisation limits

        END IF  ! of lcal360 and nftunit etc.

      ELSE IF (ft_steps(nftunit) >  0) THEN ! Regular reinit.

        ! Select files to be reinitialised on this timestep:
        ! First initialisation may be on any timestep, not just 0
        IF ((type_letter_1(nftunit) == 'c' .AND.                   &
             step-1 == ft_firststep(nftunit))                      &
        .OR.(type_letter_1(nftunit) /= 'c' .AND.                   &
             step == ft_firststep(nftunit))) THEN

          lpp_select(nftunit) = .TRUE.

        ! Deal with subsequent reinitialisation after end of period
        ELSE IF ((step-1) > ft_firststep(nftunit) .AND.            &
                (ft_laststep(nftunit) <  0 .OR.                    &
                (step-1) <= ft_laststep(nftunit))) THEN

          IF (MOD((step-1) -                                       &
              ft_firststep(nftunit),ft_steps(nftunit)) == 0)       &
            lpp_select(nftunit) = .TRUE.

          ! Do not reinitialise files on step 1 if they were
          ! initialised at step 0:
          IF (step == 1 .AND. ft_firststep(nftunit) == 0)          &
            lpp_select(nftunit) = .FALSE.

          ! Sub model id must correspond with TYPE_LETTER_2:
          lpp_select(nftunit) =  lpp_select(nftunit)               &
              .AND. internal_model == atmos_im                     &                 
              .AND. type_letter_2(nftunit) == 'a'

        END IF
!
!  1.1.1  Set up PPfile switches for boundary output files
!
        IF (type_letter_1(nftunit) == 'b') THEN ! Boundary File

!               Get interface area number
! DEPENDS ON: intf_area
          CALL intf_area( internal_model, nftunit, jintf)

          lpp_select(nftunit) = ( lpp_select(nftunit)              &
          .AND. .NOT.                                              &
!               Do not reinitialise first file on timestep ft_firststep
          (step-1-ft_firststep(nftunit) == 0) )                    &
          .OR.                                                     &
!               Initialise first file if start of sequence is offset
!               from beginning of run
          (step-1+interface_stepsim(jintf,internal_model) -        &
           ft_firststep(nftunit)  ==  0)                           &
          .OR.                                                     &
!               Select boundary file if incomplete on continuation run
          (linitial                                                &
          .AND.                                                    &
           step+interface_stepsim(jintf,internal_model) >          &
                                         ft_firststep(nftunit)     &
          .AND.                                                    &
           step <= interface_lstepim(jintf,internal_model)         &
          .AND.                                                    &
          (step-ft_firststep(nftunit) == 0 .OR.                    &
           MOD(step-ft_firststep(nftunit),ft_steps(nftunit)) /= 0) &
            )

        END IF

      ELSE  ! for files not reinitialised, ie. ft_steps(nftunit)=0

!           Initialise at step 0
        lpp_select(nftunit) = step == 0 .AND.                      &
        (ft_steps(nftunit) == 0.OR.ft_firststep(nftunit) == 0)

!           Select boundary file if incomplete on continuation run
        IF (linitial .AND.                                         &
          type_letter_1(nftunit) == 'b') THEN ! Boundary File

!             Get interface area number
! DEPENDS ON: intf_area
          CALL intf_area( internal_model, nftunit, jintf)

          lpp_select(nftunit) = lpp_select(nftunit) .OR.           &
          (step >  0.AND.step <=                                   &
                       interface_lstepim(jintf,internal_model))

        END IF

      END IF  ! of ft_steps(NFTUNIT) lt, gt or =0, ie. reinit. type
    ELSE ! of ft_output(NFTUNIT) =Y
      lpp_select(nftunit)=.FALSE.
    END IF
    lpp = lpp .OR. lpp_select(nftunit)
  END DO  ! of loop over nftunit from 20 to nunits

!
!  1.2   Set switches for general internal models.
!        For coupled models dump related switches can only be set when
!        the last internal model in a submodel has completed its group
!        of timesteps. For coupled models the only safe restart point
!        is at the completion of all groups within a model timestep.

IF(n_internal_model == 1.OR.(                                     &
                                     ! if not coupled model, or
   last_im_in_sm(internal_model).AND.                             &
                                      ! last model in submodel
   MOD(step,steps_per_period) == 0))                              &
                                      ! and last step in group
   THEN

!  ldump   : Write-up dump on this timestep

  IF (dumpfreq >  0) THEN
    ldump=       (MOD(step,dumpfreq)    == 0)
  ELSE
    ldump=.FALSE.
    DO itime=1,dumptimes_len1
      ldump=ldump.OR.(step == dumptimes(itime))
    END DO
  END IF

!  LMEAN   : Perform climate-meaning from dumps on this timestep
  lmean = .FALSE.
  IF (dumpfreq >  0.AND.meanfreq(1) >  0) THEN
    lmean=     (MOD(step,dumpfreq)       == 0)
  END IF

!  LEXIT   : Check for exit condition on this timestep
  IF (exitfreq >  0) THEN ! Implies climate meaning

    lexit=  ( (MOD(step,exitfreq)    == 0)  .OR.                  &
              (step  >=  target_end_step) )

  ELSE                    ! No climate meaning

    lexit=    (step  >=  target_end_step)

  END IF

!  lancillary: Update ancillary fields on this timestep

!   Convert ancillary reference time to days & secs
! DEPENDS ON: time2sec
  CALL time2sec(ancil_reftime(1),ancil_reftime(2),                &
                ancil_reftime(3),ancil_reftime(4),                &
                ancil_reftime(5),ancil_reftime(6),                &
                0,0,ancil_ref_days,ancil_ref_secs,lcal360)

!   Compute offset in timesteps of basis time from ancillary ref.time
! DEPENDS ON: tim2step
  CALL tim2step(basis_time_days-ancil_ref_days,                   &
                basis_time_secs-ancil_ref_secs,                   &
                steps_per_period,secs_per_period,ancil_offset_steps)

  IF (ancillary_steps >  0.AND.step >  0)                         &
    lancillary=(MOD(step+ancil_offset_steps,ancillary_steps) == 0)

END IF     ! Test for non-coupled or coupled + last step in group

!  lboundary    : Update boundary fields on this timestep

IF (boundary_steps >  0)                                          &
    lboundary= (MOD(step+bndary_offset,boundary_steps)  == 0)

!  ljobrelease  : Release user jobs on this timestep

DO i=1,jobrel_len1
  modjobrel=ABS(jobrel_step(i))
  IF (modjobrel /= 0) THEN
    ljobrelease = ljobrelease .OR.                                &
    (step == jobrel_step(i) .AND. jobrel_step(i) >  0).OR.        &
    (MOD(step,modjobrel) == jobrel_offsetim .AND. jobrel_step(i) <  0)
  END IF
END DO
!
!  1.2.1  Set switches for atmosphere timestep
!
IF (internal_model == atmos_im) THEN
!  1.2.2 Set switches for all cases

! Energy correction switches
! Set on the assumption that the energy correction is evaluated at the
! end of a timestep
  IF (a_energysteps >  0) THEN
! true if this is the step on which to do energy calculation
    lenergy =  (MOD(step,a_energysteps) == 0)
! True if this is the step after last energy calculation
    lflux_reset = (MOD(step-1,a_energysteps) == 0)
  END IF


!  lassimilation: Perform data assimilation on this timestep

  lassimilation= (step >  assim_firststepim(a_im) .AND.           &
   (step-assim_firststepim(a_im)  <=                              &
         assim_stepsim(a_im)+assim_extrastepsim(a_im)))
  lassimilation=lassimilation.AND.model_assim_mode == 'Atmosphere'
!
!       Reset Data Time fields in dump header at new analysis time
!       ( also in history block )
!       ( This now includes resetting prognostic field LOOKUP headers )
!       NB: done even if assimilation suppressed by run_assim_mode
!       Set MEANLEV to -1  and ldump to TRUE to force dump of analysis
!       - otherwise set MEANLEV to 0 (ie. instantaneous)
!
  iau_resetdt = .FALSE.
  IF (l_iau .OR. run_assim_mode  ==  "NoIAU     ") THEN
    IF (step  ==  iau_dtresetstep) iau_resetdt = .TRUE.
  END IF
  IF ( (step == assim_firststepim(a_im)+assim_stepsim(a_im)       &
        .AND.lassimilation) .OR. iau_resetdt) THEN

    DO i=21,27
      a_fixhd(i)=a_fixhd(i+7)
    END DO
    DO i=1,6
      model_data_time(i)=a_fixhd(20+i)
    END DO
    DO i=1,a_prog_lookup
      a_lookup(lbyrd ,i)=a_fixhd(21)
      a_lookup(lbmond,i)=a_fixhd(22)
      a_lookup(lbdatd,i)=a_fixhd(23)
      a_lookup(lbhrd ,i)=a_fixhd(24)
      a_lookup(lbmind,i)=a_fixhd(25)
      a_lookup(lbsecd,i)=a_fixhd(26)
    END DO
    meanlev=-1
  ELSE
    meanlev=0
  END IF
!
!       Suppress assimilation using run_assim_mode if necessary
  lassimilation=lassimilation .AND. run_assim_mode == 'Atmosphere'

!  L_SW_RADIATE : Perform SW-radiation on this timestep
!  L_LW_RADIATE : Perform LW-radiation on this timestep

  IF (a_sw_radstep >  0)                                          &
    l_sw_radiate=(MOD(step-1,a_sw_radstep)  == 0)
  IF (a_lw_radstep >  0)                                          &
    l_lw_radiate=(MOD(step-1,a_lw_radstep)  == 0)
  IF (a_sw_radstep_diag >  0)                                     &
    l_sw_radiate_diag=(MOD(step-1,a_sw_radstep_diag)  == 0)
  IF (a_lw_radstep_diag >  0)                                     &
    l_lw_radiate_diag=(MOD(step-1,a_lw_radstep_diag)  == 0)

  IF (a_sw_radstep_prog >  0)                                     &
    l_sw_radiate_prog=(MOD(step-1,a_sw_radstep_prog)  == 0)
  IF (a_lw_radstep_prog >  0)                                     &
    l_lw_radiate_prog=(MOD(step-1,a_lw_radstep_prog)  == 0)

!  linterface: Output interface fields on this timestep (>1 file poss.)

  DO jintf=1,n_intf_a
    IF (interface_stepsim(jintf,atmos_im) >  0) THEN
      linterface = linterface .OR.                                &
      (MOD(step-interface_fstepim(jintf,atmos_im),                &
           interface_stepsim(jintf,atmos_im)) == 0                &
      .AND. step >= interface_fstepim(jintf,atmos_im)             &
      .AND. step <= interface_lstepim(jintf,atmos_im) )
! To allow output of interface files from model
      IF (linitial) THEN
        lpp_select(jintf+139) = .TRUE.
      END IF
    END IF
  END DO
END IF ! Test for atmos_im
!
! ----------------------------------------------------------------------
!  2. Set STASHflags to activate diagnostics this timestep
!
!     IS is section number
!     IE is item number within section
!     IL is item number within STASHlist
!     II is counter within given section/item sublist for repeated items
!     IM is cumulative active item number within section

!   Clear all STASHflags

DO is=0,nsects
  DO im=0,nitems
    sf     (im,is) = .FALSE.
    sf_calc(im,is) = .FALSE.
  END DO
END DO

!   Loop over all items in STASHlist, enabling flags for diagnostics
!   which are active this step --
!     note that atmosphere/ocean diagnostics must be discriminated

DO il =1,totitems
  modl=stlist(st_model_code  ,il)
  is  =stlist(st_sect_no_code,il)
  im  =stlist(st_item_code   ,il)
! Skip diagnostics which are not active
  IF(stlist(st_proc_no_code,il) == 0) GO TO 200
! Skip diagnostics which don't correspond to the submodel this timestep
  IF(internal_model /= modl) GO TO 200

!     STASHflag is off by default.
!     But reset ...

!       ... if required every timestep between start and end

  IF (stlist(st_freq_code,il) == 1) THEN
    IF (step >= stlist(st_start_time_code,il).AND.                &
       (step <= stlist(st_end_time_code,il).OR.                   &
        stlist(st_end_time_code,il) == st_infinite_time))         &
    THEN
      sf(im,is)=.TRUE.
      sf(0 ,is)=.TRUE.
    END IF

!       ... if required at specified times and this is one of them

  ELSE IF(stlist(st_freq_code,il) <  0) THEN
    ntab=-stlist(st_freq_code,il)
    DO it=1,nsttims
      IF(sttabl(it,ntab) == st_end_of_list) EXIT
      IF(step == sttabl(it,ntab)) THEN
        sf(im,is)=.TRUE.
        sf(0 ,is)=.TRUE.
      END IF
    END DO

!       ... if required every N timesteps and this is one of them

  ELSE IF (stlist(st_freq_code,il) >  0) THEN
    IF  (MOD((step-stlist(st_start_time_code,il)),                &
              stlist(st_freq_code,il)) == 0.AND.                  &
         step >= stlist(st_start_time_code,il).AND.               &
        (step <= stlist(st_end_time_code,il).OR.                  &
         stlist(st_end_time_code,il) == st_infinite_time))        &
    THEN
      sf(im,is)=.TRUE.
      sf(0 ,is)=.TRUE.
    END IF
  END IF

! We now want just a subset of the above where we only want to calculate
  IF (sf(im,is)) THEN
    ! See whether we are performing a replace. 
    IF ( stlist(st_proc_no_code,il) == st_replace_code ) THEN
      ! We only want to know about replaces which are not dependent on other
      ! records.  If positive this STASH request requires a calculation.
      IF ( stlist(st_input_code,il) > 0 ) THEN
        sf_calc(im,is)  = .TRUE.
        sf_calc(0 ,is)  = .TRUE.
      END IF
    ! Everything not a replace we just let through.
    ELSE
      sf_calc(im,is)  = .TRUE.
      sf_calc(0 ,is)  = .TRUE.
    END IF

  END IF

!     Next item

  200 CONTINUE
END DO
! ----------------------------------------------------------------------
!  3. If errors occurred in setting STASHflags, set error return code
IF (icode == 1) cmessage='SETTSCTL: STASH code error'
 9999 CONTINUE
IF (lhook) CALL dr_hook('SETTSCTL',zhook_out,zhook_handle)
RETURN
! ----------------------------------------------------------------------
END SUBROUTINE settsctl
