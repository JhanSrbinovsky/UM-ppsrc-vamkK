! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine LOOP_OVER_DUMPS : Loop over dumps to get boundary data
!
! Subroutine Interface :
       SUBROUTINE loop_over_dumps (unit_no_bc,um_versn,                 &
! ARGINFA Headers for atmosphere interface data sets
        fixhd_intfa, inthd_intfa, lookup_intfa,                         &
        realhd_intfa,levdepc_intfa,                                     &
      ! Row/Col DEPC for variable resolution LBCs
        rowdepc_intfa, coldepc_intfa,                                   &  
      ! Eta values for LBC levels
        lbc_eta_theta, lbc_eta_rho,                                     &
! ARGINFA end
                  tracers_active,                                       &
                  no_tracers,tr_ukca_active,no_tr_ukca, ndustbin_in,    &
                  ndustbin_out )

! Introduce type for times
      USE makebc_time_mod, ONLY: &
        time
      USE makebc_constants_mod
 
!  needed by typ_atm_fields.h:
      USE atm_fields_bounds_mod

      USE IO, ONLY :  &
          file_open,  &
          file_close, &
          buffin,     &
          setpos
      USE ereport_mod, ONLY : ereport
      USE PrintStatus_mod
      USE UM_ParVars
      USE Control_Max_Sizes
      USE Decomp_DB
      USE domain_params
      USE um_input_control_mod,  ONLY: lcal360
      USE dynamics_input_mod, ONLY: l_endgame
      USE Submodel_Mod
      USE nlstcall_mod, ONLY : model_basis_time

      USE chsunits_mod, ONLY : nunits

      USE filenamelength_mod, ONLY: filenamelength
      USE missing_data_mod, ONLY: imdi
      IMPLICIT NONE
!
! Description : Loop over the dumps and get the boundary conditions
!
! Method : For each dump, GET_BC is called to read in the data from
!          the dump and generate the boundary conditions.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MakeBC
!
! Code Description :
! Language : FORTRAN 77 + common extensions
! This code is written to UMDP3 v6 programming standards.
!
! Declarations :
!
! Global Variables :
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
!*L------------------ COMDECK CINTFA ----------------------------------
!L CMAXSIZE should be called first.
!
!   Contains Variables, Headers and Index blocks for control of
!   generation of boundary information for the limited area model.
!
!   Interfaces to all other models are handled by STASH, and there is
!   no explicit coding written for them in the model.
!
! Interface variables initialised through INTFCNSTA
! Namelist read in the interface control routine INTF_CTL.

      INTEGER                                                           &
        intf_row_length                                                 &
                         ! Interface field row length
       ,intf_p_rows                                                     &
                         ! Interface field no of rows
       ,intf_p_levels                                                   &
                         ! Interface field no of levels
       ,intf_q_levels                                                   &
                         ! Interface field no of wet levels
       ,intf_tr_levels                                                  &
                         ! Interface field no of tracer levels
       ,intfwidtha                                                      &
                         ! Width of interface zone (atmosphere)
       ,intf_exthalo_ns                                                 &
                         ! Extended Halo in NS direction
       ,intf_exthalo_ew                                                 &
                         ! Extended Halo in EW direction
       ,a_intf_start_hr                                                 &
                         ! ) Start and End time in
       ,a_intf_freq_hr                                                  &
                         ! ) hours, Frequency in h,m,s for which
       ,a_intf_freq_mn                                                  &
                         ! ) atmosphere interface data
       ,a_intf_freq_sc                                                  &
                         ! ) is to be generated.
       ,a_intf_end_hr                                                   &
                         ! )
       ,intf_pack                                                       &
                         ! Packing Indicator for boundary data
       ,lbc_stream_a                                                    &
                         ! Output streams in UMUI
       ,lbc_unit_no_a                                                   &
                         ! Unit Nos for Atmos Boundary Dataset
       ,lbc_first_r_rho                                                 &
                         ! First rho level at which height is constant
       ,intf_v_int_order(max_n_intf_a)

      REAL                                                              &
        intf_ewspace                                                    &
                         ! E-W grid spacing (degrees)
       ,intf_nsspace                                                    &
                         ! N-S grid spacing (degrees)
       ,intf_firstlat                                                   &
                         ! Latitude of first row (degrees)
       ,intf_firstlong                                                  &
                         ! Longitude of first row (degrees)
       ,intf_polelat                                                    &
                         ! Real latitude of coordinate pole (degrees)
       ,intf_polelong                                                   &
                         ! Real longitude of coordinate pole (degrees)
       ,lbc_z_top_model                                                 &
                         ! Height of top of model
       ,lbc_q_min                                                       &
                         ! Minimum value for q
!
! VarRes grid spacing
      , lambda_intf_p(max_intf_lbcrow_length, max_n_intf_a)             &
      , lambda_intf_u(max_intf_lbcrow_length, max_n_intf_a)             &    
      , phi_intf_p(max_intf_lbcrows, max_n_intf_a)                      &
      , phi_intf_v(max_intf_lbcrows, max_n_intf_a)

      LOGICAL                                                           &
        intf_vert_interp                                                &
                         ! Switch to request vertical interpolation
       ,lnewbnd          ! True for initialising new boundary data file

! Switch for variable resolution LBC output
      LOGICAL  intf_l_var_lbc(max_n_intf_a)

! Switch to not rotate if input and output grids have same poles.
      LOGICAL intf_avoid_rot(MAX_N_INTF_A)

! Switch to output LBC for Endgame
      LOGICAL intf_l_eg_out(MAX_N_INTF_A)

! Files for VERTLEVS namelist     
      CHARACTER(LEN=256) :: intf_vertlevs

! Files for HorzGrid namelist  
      CHARACTER(LEN=256) :: intf_HorzGrid(max_n_intf_a)
!*----------------------------------------------------------------------
      COMMON /INTFCTL_ATMOS/                                            &
        intf_ewspace(max_n_intf_a)    ,intf_nsspace(max_n_intf_a),      &
        intf_firstlat(max_n_intf_a)   ,intf_firstlong(max_n_intf_a),    &
        intf_polelat(max_n_intf_a)    ,intf_polelong(max_n_intf_a),     &
        intf_row_length(max_n_intf_a) ,intf_p_rows(max_n_intf_a),       &
        intf_p_levels(max_n_intf_a)   ,intf_q_levels(max_n_intf_a),     &
        intf_tr_levels(max_n_intf_a)  ,intfwidtha(max_n_intf_a),        &
        intf_exthalo_ns(max_n_intf_a) ,intf_exthalo_ew(max_n_intf_a),   &
        a_intf_start_hr(max_n_intf_a) ,a_intf_freq_hr(max_n_intf_a),    &
        a_intf_freq_mn(max_n_intf_a)  ,a_intf_freq_sc(max_n_intf_a),    &
        a_intf_end_hr(max_n_intf_a)   ,                                 & 
        lnewbnd(max_n_intf_a)         ,intf_vert_interp(max_n_intf_a),  &
        intf_pack(max_n_intf_a)       ,lbc_stream_a(max_n_intf_a),      &
        lbc_unit_no_a(max_n_intf_a)   ,lbc_first_r_rho(max_n_intf_a),   &
        lbc_z_top_model(max_n_intf_a) ,                                 &
        intf_vertlevs(max_n_intf_a)   ,lbc_q_min,                       &
        intf_l_var_lbc                ,intf_horzgrid,                   &
        lambda_intf_p                 ,lambda_intf_u,                   &
        phi_intf_p                    ,phi_intf_v,                      &
        intf_avoid_rot                ,intf_v_int_order,                &
        intf_l_eg_out
!---------------------------------------------------------------------

! Subroutine arguments
!   Scalar arguments with intent(in) :

      INTEGER :: unit_no_bc  ! Unit No for Boundary dataset
      INTEGER :: um_versn    ! UM Version Boundary dataset for
      INTEGER, INTENT(IN) :: ndustbin_in   ! Number of dust bins in input file
      INTEGER, INTENT(IN) :: ndustbin_out ! Number of dust bins in output LBCs

!   Array arguments with intent(in) :



      INTEGER           :: errorstatus    ! Error code
      CHARACTER(LEN=80) :: cmessage       ! Error Message

!   Local scalars

      INTEGER                :: unit_no      !  Unit no for input dump
      INTEGER                :: read_only=0  !  Input dumps - read only
      INTEGER                :: read_write=1 !  Output Boundary File - r/w
      INTEGER,PARAMETER      :: len_env=7    !  Length of env. variable
      INTEGER                :: env_var=0    !  Filename is in env var
      CHARACTER(LEN=len_env) :: env          !  Input filename env. var.

      INTEGER :: j,jdump                   !  Loop indices
      INTEGER :: inthd(46)
      INTEGER :: elapsed_days              !  No of days elapsed
      INTEGER :: elapsed_secs              !  No of secs elapsed
      INTEGER :: len_actual                !  Length of data read in BUFFIN
      INTEGER :: p_rows

!   Local dynamic arrays

      INTEGER ::  fixhd(len_fixhd)  !  Fixed header from dump

! Declarations for vars that have been removed from TYPSIZE.
      INTEGER :: p_field
      INTEGER :: u_field

! Declarations for fieldsfile reading
      INTEGER                     :: test(64)
      INTEGER                     :: len_io
      INTEGER                     :: open_loop
      REAL                        :: a_io
      CHARACTER(LEN=80),PARAMETER :: routinename='loop_over_dumps'
      CHARACTER(LEN=len_env)      :: bc_filename
      ! Contents of environment variable.
      CHARACTER(LEN=filenamelength) :: env_contents

! Max_progs gives the size required for indexing
! arrays in uniform manner
      INTEGER :: max_progs

! Basis_time is the starting time.
      TYPE(time) :: basis_time
! The current time being processed.
      TYPE(time) :: current_time
! Flag to indicate all required data is present
      LOGICAL              :: alldata
!
! Free and UKCA tracer arguments
      INTEGER :: no_tracers, no_tr_ukca
      INTEGER :: tracers_active(max_tracers)
      INTEGER :: tr_ukca_active(max_tr_ukca)

!-  End of Header
! Set max_progs 
! Now use parameters from makebc_constants_mod instead of a magic number
      max_progs=num_req_lbcs+num_optional_lbcs_max+no_tracers+no_tr_ukca


! Initialise all items in basis_time to 0
      basis_time%year   = 0
      basis_time%month  = 0
      basis_time%day    = 0
      basis_time%hour   = 0
      basis_time%min    = 0
      basis_time%sec    = 0
      basis_time%day_no = 0
 
! Initialise all items in current_time to imdi
! This will identify if current_time is never set.
      current_time%year   = imdi 
      current_time%month  = imdi 
      current_time%day    = imdi 
      current_time%hour   = imdi 
      current_time%min    = imdi 
      current_time%sec    = imdi 
      current_time%day_no = imdi 

!     Open Boundary Dataset
      WRITE (6,*) ' '

!     Loop over number of LAMs, opening boundary file for each
      DO open_loop=1,n_intf_a
        WRITE(bc_filename,'(A,I2.2)')'BCFIL',open_loop
        IF(printstatus >= prstatus_diag)THEN
          WRITE(6,*)' bc_filename= ',bc_filename
        END IF
        CALL file_open                                                  &
           (unit_no_bc+(open_loop-1),bc_filename,len_env,read_write,    &
            env_var,errorstatus)
        IF (errorstatus /= 0) THEN
          WRITE (cmessage,'(A,I6)')                                     &
            'Error in opening Boundary Dataset on unit no ',            &
            unit_no_bc

          CALL ereport(routinename,errorstatus,cmessage)
        END IF
      END DO

!     Loop over model dumps
      jdump = 0
      DO
!       Increment file number to attempt to open.
        jdump=jdump + 1
!       Unit number for this dump
        unit_no = 30

!       Open the dump
        WRITE (env,'(A,I3.3)') 'FILE', unit_no + jdump

        CALL fort_get_env(env,len_env,env_contents,filenamelength, &
                          errorstatus)

        IF (errorstatus /= 0) THEN
          ! Reached end of files to input exit out and check we have all times
          ! outside of loop.
          EXIT
        END IF
 
        WRITE (6,'(A)') ' '
        WRITE (6,'(A,I3)') 'Processing dump no ',jdump


        WRITE (6,*) ' '
        CALL file_open (unit_no,env,len_env,read_only,env_var,errorstatus)
        IF (errorstatus /= 0) THEN
          WRITE (cmessage,'(A,I6)') 'Error in opening dump on unit no ',unit_no

          CALL ereport(routinename,errorstatus,cmessage)
        END IF

!       Read in fixed header from this dump
        CALL setpos (unit_no,0,errorstatus)
        IF (errorstatus >  0) THEN
          cmessage = 'Error in SETPOS for Fixed Header.'

          CALL ereport(routinename,errorstatus,cmessage)
        END IF

! DEPENDS ON: read_flh
        CALL read_flh (unit_no,fixhd,len_fixhd,errorstatus,cmessage)
        IF (errorstatus >  0) THEN
          ! Cmessage should be set from read_flh.
          WRITE (6,'(A,I6)') 'Error in READ_FLH for dump ',jdump

          CALL ereport(routinename,errorstatus,cmessage)
        END IF

        IF (fixhd(5) == 1)THEN
          WRITE(6,*)'this is a Dump'
        ELSEIF (fixhd(5) == 3)THEN
          WRITE(6,*)'this is a FieldsFile'
        ELSE
          cmessage= 'Invalid input file type'
          errorstatus=1

          CALL ereport( RoutineName, errorstatus, Cmessage )
        END IF
!       For first dump only
!       Set model basis time to be date/time in first dump and use
!       to initialise variables required for time processing
        IF (jdump == 1) THEN
          IF(fixhd(150) > 0)THEN
            CALL setpos(unit_no,fixhd(150)-1,errorstatus)
            IF (errorstatus >  0) THEN
              cmessage = 'Error in SETPOS for lookups.'

              CALL ereport(routinename,errorstatus,cmessage)
            END IF
            CALL buffin(unit_no, test, fixhd( 151 ),                      &
                                Len_IO, a_io )
            IF (a_io /= -1.0) THEN
              errorstatus=NINT(a_io)
              WRITE (cmessage,'(A,I9,A)') 'Problem reading ',             &
                len_actual, ' from first lookup.'
              CALL ereport(routinename,errorstatus,cmessage)
            END IF

            model_basis_time(1)=test(1)
            model_basis_time(2)=test(2)
            model_basis_time(3)=test(3)
            model_basis_time(4)=test(4)
            model_basis_time(5)=0
            model_basis_time(6)=0

            ! Now set the reference time variables which are used later to
            ! calculate times for headers in in_intf
            basis_time_days = 0
            basis_time_secs = 0
            CALL time2sec(model_basis_time(1), &
                          model_basis_time(2), &
                          model_basis_time(3), &
                          model_basis_time(4), &
                          model_basis_time(5), &
                          model_basis_time(6), &
                          basis_time_days, basis_time_secs, &
                          elapsed_days, elapsed_secs, lcal360)
            basis_time_days = elapsed_days
            basis_time_secs = elapsed_secs
          ELSE
            ! Exit since we cannot go any further.  
            errorstatus=2
            WRITE (cmessage,'(A)') ' No Lookups to get basis time from'
            CALL ereport(routinename,errorstatus,cmessage)
          END IF

! Set the basis_time from the model basis time
          basis_time%year = model_basis_time(1)
          basis_time%month= model_basis_time(2)
          basis_time%day  = model_basis_time(3)
          basis_time%hour = model_basis_time(4)
          basis_time%min  = model_basis_time(5)
          basis_time%sec  = model_basis_time(6)

        END IF

!       Remove negative dimensions, if any.
        DO j=100,256
          fixhd(j) = MAX(fixhd(j),0)
        END DO

! Get header dimensions from Fixed Header (ensure sizes are >= 0 for allocating)
        a_len_inthd    = MAX(fixhd(101),0)
        a_len_realhd   = MAX(fixhd(106),0)
        a_len1_levdepc = MAX(fixhd(111),0)
        a_len2_levdepc = MAX(fixhd(112),0)
        a_len1_rowdepc = MAX(fixhd(116),0)
        a_len2_rowdepc = MAX(fixhd(117),0)
        a_len1_coldepc = MAX(fixhd(121),0)
        a_len2_coldepc = MAX(fixhd(122),0)
        a_len1_flddepc = MAX(fixhd(126),0)
        a_len2_flddepc = MAX(fixhd(127),0)
        a_len_extcnst  = MAX(fixhd(131),0)
        a_len_cfi1     = MAX(fixhd(141),0)
        a_len_cfi2     = MAX(fixhd(143),0)
        a_len_cfi3     = MAX(fixhd(145),0)
        a_len2_lookup  = MAX(fixhd(152),0)
        a_len_data     = MAX(fixhd(161),0)

!       Get length of data in this dump
        len_tot  = a_len_data

!       Read in Integer Constants for this dump
        CALL setpos(unit_no,fixhd(100)-1,errorstatus)
        IF (errorstatus >  0) THEN
          cmessage = 'Error in SETPOS for Integer Constants array.'

          CALL ereport(routinename,errorstatus,cmessage)
        END IF

        CALL buffin (unit_no,inthd(1:),a_len_inthd,len_actual,a_io)
        IF (a_io /= -1.0) THEN
          errorstatus=NINT(a_io)
          WRITE (cmessage,'(A,I9,A)') 'Problem reading length ',          &
            len_actual, ' from integer constants array.'

          call ereport(routinename,errorstatus,cmessage)
        END IF

!       Get model grid for this dump
        row_length = test(19)   ! from lookup
        rows       = test(18)   ! from lookup
        IF (fixhd(9) == 3) THEN
          ! Arakawa C
          n_rows = rows-1
        ELSE IF (fixhd(9) == 6) THEN
          ! Alternative Arakawa C (ENDGAME)
          n_rows = rows+1
        END IF
        p_field    = row_length * rows
        u_field    = row_length * n_rows

!     Get model levels for this dump
        model_levels = inthd(8)
        wet_levels   = MAX(inthd(9),1)
        tr_levels    = MAX(inthd(12),1)
        bl_levels    = MAX(inthd(13),1)
        ozone_levels = MAX(inthd(26),1)
!     Set whether this dump is a endgame dump or not.
!     This is really checking the staggering but will have to do for now.
        l_endgame    = fixhd(9) == 6


        IF(printstatus >= prstatus_normal)THEN
          WRITE (6 , '(A)') ' '
          WRITE (6 , '(A)') ' Model Grid/Levels in this dump. '
          WRITE (6 , '(A,I9)') ' row_length   = ', row_length
          WRITE (6 , '(A,I9)') ' rows         = ', rows
          WRITE (6 , '(A,I9)') ' n_rows       = ', n_rows
          WRITE (6 , '(A,I9)') ' model_levels = ', model_levels
          WRITE (6 , '(A,I9)') ' wet_levels   = ', wet_levels
          WRITE (6 , '(A,I9)') ' tr_levels    = ', tr_levels
          WRITE (6 , '(A,I9)') ' bl_levels    = ', bl_levels
          WRITE (6 , '(A,I9)') ' ozone_levels = ', ozone_levels
          WRITE (6 , '(A,I9)') ' p_field      = ', p_field
          WRITE (6 , '(A,I9)') ' u_field      = ', u_field
        END IF

! Get decompostion information
! ----------------------------
! LEVELS not used in SX
        CALL decompose( row_length, rows,                           &
                           0, 0, model_levels)
! Formerly set as a side-effect of readflds, now specified here explicitly
        CALL change_decomposition(decomp_smexe, errorstatus)

! Set dimensions for ENDGame and NewDynamics
        CALL atm_fields_bounds_init(offx,offy,halo_i,halo_j,              &
                    row_length,rows,n_rows,model_levels,wet_levels,       &
                    tr_levels,bl_levels,ozone_levels)


!     Proceed to read model data from dump and get boundary conditions
! DEPENDS ON: get_bc
      call get_bc (jdump,unit_no,unit_no_bc,um_versn,                   &
     &  a_len2_lookup,                                                  &
     &  Fixhd_intfa,Inthd_intfa,Lookup_intfa,                           &
     &  Realhd_intfa,Levdepc_intfa,                                     &
     &  Rowdepc_intfa,Coldepc_intfa,                                    &
! Pass basis_time
          basis_time,                                                   &
          current_time,                                                 &
! Pass indicator alldata
          alldata,                                                        &
          lbc_eta_theta,lbc_eta_rho,                                      &
          max_progs,                                                      &
! Pass details of tracers in
          tracers_active, no_tracers, tr_ukca_active, no_tr_ukca,         &
          ndustbin_in, ndustbin_out )
!
!       Close the dump
        CALL file_close (unit_no,env,len_env,env_var,0,errorstatus)
        IF (errorstatus /= 0) THEN
          WRITE (cmessage,'(A,I6)') 'Error in closing dump on unit no ', unit_no

          CALL ereport(routinename,errorstatus,cmessage)
        END IF

        IF(alldata) THEN
!         If we have all the data then we can quit.
          EXIT
        END IF

      END DO   !    End of loop over dumps

      IF(.NOT. alldata)THEN
          WRITE(6,'(A)') "Missing time for:"
          WRITE(6,'(I6)') current_time
          cmessage='ERROR data missing for one or more LBC time'
          errorstatus=1
          CALL ereport(routinename,errorstatus,cmessage)
      END IF

      RETURN
      END SUBROUTINE loop_over_dumps
