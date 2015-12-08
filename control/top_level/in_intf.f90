! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL Subroutine IN_INTF
!LL
!LL Purpose : Takes as input, codes set by the user interface defining
!LL           the start time, end time, and interval for creating
!LL           interface data for a limited area model, and data
!LL           defining the limited area grid. The source model may also
!LL           be limited area. Sets up fixed length, integer & real
!LL           headers and level dependent constants for the interface
!LL           data set. All prognostic variables for which horizontal
!LL           differencing is performed are included, ie all tracers
!LL           but no surface fields.
!LL
!LL Programming standard; Unified Model Documentation Paper No. 3
!LL version no. 1, dated 15/01/90
!LL
!LL Logical components covered : D810
!LL
!LL System task : D81
!LL
!LL Documentation : Unified Model Documentation Paper No D8
!LL
!LLEND -------------------------------------------------------------
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Top Level

      SUBROUTINE IN_INTF (                                              &

!*L   Arguments:
! ARGDUMA Dump headers
        A_FIXHD, A_INTHD, A_CFI1, A_CFI2, A_CFI3, A_REALHD, A_LEVDEPC,  &
        A_ROWDEPC, A_COLDEPC, A_FLDDEPC, A_EXTCNST, A_DUMPHIST,         &
      ! PP lookup headers and Atmos stash array + index with lengths
        A_LOOKUP,A_MPP_LOOKUP,a_ixsts, a_spsts,                         &
! ARGDUMA end
! ARGINFA Headers for atmosphere interface data sets
        fixhd_intfa, inthd_intfa, lookup_intfa,                         &
        realhd_intfa,levdepc_intfa,                                     &
      ! Row/Col DEPC for variable resolution LBCs
        rowdepc_intfa, coldepc_intfa,                                   &  
      ! Eta values for LBC levels
        lbc_eta_theta, lbc_eta_rho,                                     &
! ARGINFA end
                 NFTOUT,ICODE,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE filenamelength_mod, ONLY:    filenamelength
      USE IO
      USE io_configuration_mod, ONLY : io_data_alignment
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Control_Max_Sizes
      USE um_input_control_mod, ONLY: lcal360
      USE nlstgen_mod, ONLY: steps_per_periodim, secs_per_periodim,     &
          llboutim
      USE Submodel_Mod


      USE nlstcall_mod, ONLY : ft_steps, &
                               ft_firststep

      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE

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

      INTEGER  NFTOUT  ! FTN Number to write interface data
      INTEGER  ICODE   ! Return code

      CHARACTER(LEN=80)  CMESSAGE  ! Error message

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

!    Local variables
      INTEGER    I,J,IJ    ! DO loop indices
      INTEGER    ROW, IRIM 
      INTEGER    DAYS,SECS ! Time difference relative to reference point
      INTEGER    LEVEL                                                    
      INTEGER    LEN_PP    ! Total length of PP headers
      INTEGER    LEN_IO    ! Total length of data read in
      INTEGER    NTIMES    ! Number of times when interface data is required
      INTEGER    ZERO
      INTEGER    JINTF     ! Interface area index
!L----------------------------------------------------------------------

      REAL  A_IO
      REAL, ALLOCATABLE :: Dumm1(:)   ! dummy array to fill file

      INTEGER YY,MM,DD,HH,MN,SEC,DAY_NO
      INTEGER A_SECS_PER_STEP ! secs per step for atmos. sub_model
      INTEGER A_INTF_FREQ_SECS  ! Interface frequency in seconds
      INTEGER :: ntime
      INTEGER :: lookup_start

      INTEGER lbc_len_inthd, lbc_len_realhd, lbc_lookupsa
      CHARACTER(LEN=filenamelength) :: string      ! MODEL_FT_UNIT value
      CHARACTER(LEN=14) PPNAME           ! Boundary data name
      Character (Len=*), Parameter :: RoutineName = 'IN_INTF'

      INTEGER :: num_rows

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('IN_INTF',zhook_in,zhook_handle)

! MakeBC only

!L    Internal Structure

      IF (LLBOUTim(A_IM)) THEN

! DEPENDS ON: intf_area
      CALL INTF_AREA (A_IM,NFTOUT,JINTF)

      A_INTF_FREQ_SECS=3600*A_INTF_FREQ_HR(JINTF) +                     &
        60*A_INTF_FREQ_MN(JINTF)+A_INTF_FREQ_SC(JINTF)
      IF (A_INTF_FREQ_SECS >  0) THEN

        IF (LNEWBND(JINTF)) THEN

!L      LNEWBND = true (New dataset to be set up)

! 1.1   Compute secs per step for atmosphere sub_model
        A_SECS_PER_STEP = secs_per_periodim(a_im)/                      &
                          steps_per_periodim(a_im)

!L  1.0 Set up headers

        IF (ft_steps(NFTOUT) == 0) THEN
          NTIMES = ((A_INTF_END_HR(JINTF)-A_INTF_START_HR(JINTF))*3600)/&
                    A_INTF_FREQ_SECS+1
        ELSE ! reinitialisation
          NTIMES = FT_STEPS(NFTOUT)/                                    &
                     (A_INTF_FREQ_SECS / A_SECS_PER_STEP)
          IF (STEPim(a_im)-1 <= FT_FIRSTSTEP(NFTOUT)) THEN ! 1st file
            NTIMES = NTIMES + 1
          ENDIF
        ENDIF

        lbc_len_inthd  = pp_len_inthd
        lbc_len_realhd = pp_len_realhd
        ! Since each 3D field is one lookup this is the number of LBC fields.
        lbc_lookupsa   = intf_lookupsa

        IF (intf_l_eg_out(jintf)) THEN
          num_rows = INTF_P_ROWS(JINTF)+1
        ELSE
          num_rows = INTF_P_ROWS(JINTF)
        END IF


!L  1.1 Fixed length header
!       Copy main dump header as first step

      DO I=1,LEN_FIXHD
        FIXHD_INTFA(I,JINTF)=A_FIXHD(I)
      ENDDO

!        Set boundary data identifiers

      FIXHD_INTFA(1,JINTF)  = IMDI
      FIXHD_INTFA(5,JINTF)  = 5
! Set staggering type to be different

      IF (intf_l_eg_out(jintf)) THEN
        ! Arakawa C alternative - VATPOLES
        FIXHD_INTFA(9,JINTF) = 6
      ELSE
        ! Arakawa C - New Dynamics
        FIXHD_INTFA(9,JINTF) = 3
      END IF
      FIXHD_INTFA(10,JINTF) = 1

!        Modify individual items

!L      Calculate start validity time

      IF (ft_steps(NFTOUT) == 0) THEN
        SECS = A_INTF_START_HR(JINTF)*3600
        DAYS = 0
      ELSE
       IF (STEPim(a_im)-1 <= ft_firststep(NFTOUT)) THEN   ! First file
! DEPENDS ON: stp2time
        CALL STP2TIME(FT_FIRSTSTEP(NFTOUT),steps_per_periodim(a_im),    &
                      secs_per_periodim(a_im),DAYS,SECS)
       ELSE ! not first file
! DEPENDS ON: stp2time
        CALL STP2TIME(STEPim(a_im)-1,steps_per_periodim(a_im),          &
                      secs_per_periodim(a_im),DAYS,SECS)
        SECS=SECS+ A_INTF_FREQ_SECS
       ENDIF
      ENDIF

! DEPENDS ON: sec2time
      CALL SEC2TIME(DAYS,SECS,BASIS_TIME_DAYS,BASIS_TIME_SECS,          &
                    YY,MM,DD,HH,MN,SEC,DAY_NO,LCAL360)

      FIXHD_INTFA(21,JINTF) = YY
      FIXHD_INTFA(22,JINTF) = MM
      FIXHD_INTFA(23,JINTF) = DD
      FIXHD_INTFA(24,JINTF) = HH
      FIXHD_INTFA(25,JINTF) = MN
      FIXHD_INTFA(26,JINTF) = SEC
      FIXHD_INTFA(27,JINTF) = DAY_NO

!     Data interval
      FIXHD_INTFA(35,JINTF) = 0
      FIXHD_INTFA(36,JINTF) = 0
      FIXHD_INTFA(37,JINTF) = A_INTF_FREQ_HR(JINTF)/24
      FIXHD_INTFA(38,JINTF) = MOD(A_INTF_FREQ_HR(JINTF),24)
      FIXHD_INTFA(39,JINTF) = A_INTF_FREQ_MN(JINTF)
      FIXHD_INTFA(40,JINTF) = A_INTF_FREQ_SC(JINTF)
      FIXHD_INTFA(41,JINTF) = A_INTF_FREQ_HR(JINTF)/24

!L      Calculate last validity time

      IF (ft_steps(NFTOUT) == 0) THEN
        SECS=A_INTF_END_HR(JINTF)*3600
        DAYS = 0
      ELSE
       IF (STEPim(a_im)-1 <= ft_firststep(NFTOUT)) THEN   ! First file
! DEPENDS ON: stp2time
        CALL STP2TIME(FT_FIRSTSTEP(NFTOUT),steps_per_periodim(a_im),    &
                      secs_per_periodim(a_im),DAYS,SECS)
        SECS=SECS + FT_STEPS(NFTOUT)*A_SECS_PER_STEP
       ELSE ! not first file
! DEPENDS ON: stp2time
        CALL STP2TIME(STEPim(a_im)-1+FT_STEPS(NFTOUT),                  &
            steps_per_periodim(a_im),secs_per_periodim(a_im),DAYS,SECS)
       ENDIF
      ENDIF

! DEPENDS ON: sec2time
      CALL SEC2TIME(DAYS,SECS,BASIS_TIME_DAYS,BASIS_TIME_SECS,          &
                    YY,MM,DD,HH,MN,SEC,DAY_NO,LCAL360)

      FIXHD_INTFA(28,JINTF) = YY
      FIXHD_INTFA(29,JINTF) = MM
      FIXHD_INTFA(30,JINTF) = DD
      FIXHD_INTFA(31,JINTF) = HH
      FIXHD_INTFA(32,JINTF) = MN
      FIXHD_INTFA(33,JINTF) = SEC
      FIXHD_INTFA(34,JINTF) = DAY_NO

!L      Modify header lengths

      FIXHD_INTFA(101,JINTF)=lbc_len_inthd
      FIXHD_INTFA(105,JINTF) =                                          &
      FIXHD_INTFA(100,JINTF)+FIXHD_INTFA(101,JINTF)
      FIXHD_INTFA(106,JINTF)=lbc_len_realhd
      FIXHD_INTFA(110,JINTF)=                                           &
      FIXHD_INTFA(105,JINTF)+FIXHD_INTFA(106,JINTF)

!L      Set length of level dependent constants

      FIXHD_INTFA(112,JINTF)=4
      IF (INTF_VERT_INTERP(JINTF)) THEN
        FIXHD_INTFA(111,JINTF)=INTF_P_LEVELS(JINTF)+1
      ELSE
        FIXHD_INTFA(111,JINTF)=MODEL_LEVELS+1
      END IF

!  NO row and column dependent constants
      FIXHD_INTFA(115,JINTF)=IMDI
      FIXHD_INTFA(116,JINTF)=IMDI
      FIXHD_INTFA(117,JINTF)=IMDI
      FIXHD_INTFA(120,JINTF)=IMDI
      FIXHD_INTFA(121,JINTF)=IMDI
      FIXHD_INTFA(122,JINTF)=IMDI
      
      IF (INTF_L_VAR_LBC(JINTF)) THEN
        FIXHD_INTFA(115,JINTF)=                                         &
                FIXHD_INTFA(110,JINTF)+                                 &
                FIXHD_INTFA(111,JINTF)*FIXHD_INTFA(112,JINTF)
        FIXHD_INTFA(116,JINTF)= num_rows
        FIXHD_INTFA(117,JINTF)=2
        FIXHD_INTFA(120,JINTF)=                                         &
                FIXHD_INTFA(115,JINTF)+                                 &
                FIXHD_INTFA(116,JINTF)*FIXHD_INTFA(117,JINTF)
        FIXHD_INTFA(121,JINTF)=INTF_ROW_LENGTH(JINTF)
        FIXHD_INTFA(122,JINTF)=2
      END IF
      
!  NO field_constants,extra constants,temp_historyfile or compressed
!  indexes
      DO I=125,145
        FIXHD_INTFA(I,JINTF)=IMDI
      ENDDO

!     Start address and length of look up table
      FIXHD_INTFA(150,JINTF) = FIXHD_INTFA(110,JINTF) +                 &
                        FIXHD_INTFA(111,JINTF) * FIXHD_INTFA(112,JINTF)

      IF (INTF_L_VAR_LBC(JINTF)) THEN
        FIXHD_INTFA(150,JINTF) = FIXHD_INTFA(120,JINTF) +               &
                        FIXHD_INTFA(121,JINTF) * FIXHD_INTFA(122,JINTF)
      END IF
       
      FIXHD_INTFA(152,JINTF) = NTIMES * LBC_LOOKUPSA
      FIXHD_INTFA(153,JINTF) = IMDI

!     Start address and length of data section
!--make sure the data starts on a sector bndry
      fixhd_intfa(160, jintf)=((fixhd_intfa(150, jintf)+                &
       fixhd_intfa(151, jintf)*fixhd_intfa(152, jintf)+                 &
       io_data_alignment-1)/io_data_alignment)*io_data_alignment+1
      FIXHD_INTFA(161,JINTF) = 0

!L  1.2 Integer header

      DO I=1,PP_LEN_INTHD
        INTHD_INTFA(I,JINTF) = A_INTHD(I)
      ENDDO

      INTHD_INTFA(1,JINTF) = INTERFACE_FSTEPim(JINTF,a_im)
      INTHD_INTFA(2,JINTF) = A_INTF_FREQ_HR(JINTF)
      INTHD_INTFA(3,JINTF) = NTIMES
      INTHD_INTFA(6,JINTF) = INTF_ROW_LENGTH(JINTF)
      INTHD_INTFA(7,JINTF) = INTF_P_ROWS(JINTF)
      INTHD_INTFA(8,JINTF) = INTF_P_LEVELS(JINTF)
      INTHD_INTFA(9,JINTF) = INTF_Q_LEVELS(JINTF)
      INTHD_INTFA(12,JINTF)= INTF_TR_LEVELS(JINTF)
      INTHD_INTFA(15,JINTF)= LBC_LOOKUPSA

      DO I=16,46
        INTHD_INTFA(I,JINTF) = IMDI
      END DO

!     Algorithm used for generating heights.
      INTHD_INTFA(17,JINTF) = 2

!     First rho level at which height is constant
      IF (Intf_Vert_Interp(jintf)) THEN
        INTHD_INTFA(24,JINTF) = LBC_FIRST_R_RHO(JINTF)
      ELSE
        INTHD_INTFA(24,JINTF) = A_INTHD(24)
      END IF

!L  1.3 Real header

      DO I=1,PP_LEN_REALHD
       REALHD_INTFA(I,JINTF) = A_REALHD(I)
      ENDDO

      REALHD_INTFA(1,JINTF) = INTF_EWSPACE(JINTF)
      REALHD_INTFA(2,JINTF) = INTF_NSSPACE(JINTF)
      REALHD_INTFA(3,JINTF) = INTF_FIRSTLAT(JINTF)
      REALHD_INTFA(4,JINTF) = INTF_FIRSTLONG(JINTF)
      REALHD_INTFA(5,JINTF) = INTF_POLELAT(JINTF)
      REALHD_INTFA(6,JINTF) = INTF_POLELONG(JINTF)

      DO I=7,38
        REALHD_INTFA(I,JINTF) = RMDI
      END DO

!     Height at top of model
      IF (Intf_Vert_Interp(jintf)) THEN
        REALHD_INTFA(16,JINTF)= LBC_Z_TOP_MODEL(JINTF)
      ELSE
        REALHD_INTFA(16,JINTF)= A_REALHD(16)
      END IF

!L  1.4 Level dependent constants

        IF (INTF_VERT_INTERP(JINTF)) THEN
!L      Set level dependent constants from namelist INTFCNST  if
!L      vertical interpolation required

          levdepc_intfa(:,:,jintf) = rmdi

          DO LEVEL=1,INTF_P_LEVELS(JINTF)+1
            LEVDEPC_INTFA(LEVEL,1,JINTF) = LBC_ETA_THETA(LEVEL,JINTF)
          ENDDO
          DO LEVEL=1,INTF_P_LEVELS(JINTF)
            LEVDEPC_INTFA(LEVEL,2,JINTF) = LBC_ETA_RHO(LEVEL,JINTF)
          ENDDO

        ELSE

!       Copy level dependent constants from source model

          LevDepC_Intfa(:,:,jintf) = rmdi

!         Eta - theta levels
          DO Level = 1, Model_Levels+1
            LevDepC_Intfa(Level,1,jintf) = A_LevDepc(Level)
          END DO

!         Eta - rho levels
          DO Level = 1, Model_Levels
            LevDepC_Intfa(Level,2,jintf) =                              &
            A_LevDepc(Model_Levels+1+Level)
          END DO

!         RH_Crit and Soil Moisture levels not copied.


        END IF

!L  1.5 Row/Col dependent constants
        
        rowdepc_intfa(:,:,jintf) = rmdi
        coldepc_intfa(:,:,jintf) = rmdi
        
        IF (INTF_L_VAR_LBC(JINTF)) THEN
 
          DO i = 1, num_rows
            ROWDEPC_INTFA(I,1,JINTF) =  PHI_INTF_P(I, JINTF)
            ROWDEPC_INTFA(I,2,JINTF) =  PHI_INTF_V(I, JINTF)
          ENDDO
          DO J=1, INTF_ROW_LENGTH(JINTF)
            COLDEPC_INTFA(J,1,JINTF) =  LAMBDA_INTF_P(J, JINTF)
            COLDEPC_INTFA(J,2,JINTF) =  LAMBDA_INTF_U(J, JINTF)
          ENDDO 
        END IF
        
!L  2.0 Write out headers

!L  2.1 Fixed length header


        CALL BUFFOUT(NFTOUT,FIXHD_INTFA(:,JINTF),LEN_FIXHD,LEN_IO,A_IO)

! Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of fixed length header',A_IO,LEN_IO, &
                        LEN_FIXHD)
          CMESSAGE='IN_INTF:I/O ERROR'
          ICODE=1
          IF (lhook) CALL dr_hook('IN_INTF',zhook_out,zhook_handle)
          RETURN
        END IF

!L  2.2 Integer header


        CALL BUFFOUT(NFTOUT,INTHD_INTFA(:,JINTF),                       &
                     LBC_LEN_INTHD,LEN_IO,A_IO)

!       Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= LBC_LEN_INTHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of integer header',A_IO,LEN_IO,      &
                       LBC_LEN_INTHD)
          CMESSAGE='IN_INTF:I/O ERROR'
          ICODE=2
          IF (lhook) CALL dr_hook('IN_INTF',zhook_out,zhook_handle)
          RETURN
        END IF

!L  2.3 Real header


        CALL BUFFOUT(NFTOUT,REALHD_INTFA(:,JINTF),                      &
                     LBC_LEN_REALHD,LEN_IO,A_IO)

!       Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= LBC_LEN_REALHD) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of real header',A_IO,LEN_IO,         &
                        LBC_LEN_REALHD)
          CMESSAGE='IN_INTF:I/O ERROR'
          ICODE=3
          IF (lhook) CALL dr_hook('IN_INTF',zhook_out,zhook_handle)
          RETURN
        END IF

!L  2.4 Level dependent constants

!       Write out each variable separately as second dimension
!       of LEVDEPC_INTFA is now a maximum dimension for all
!       interface areas being generated

        DO I=1,4


        CALL BUFFOUT(NFTOUT,LEVDEPC_INTFA(:,I,JINTF),                   &
                     FIXHD_INTFA(111,JINTF),LEN_IO,A_IO)

!       Check for I/O Errors

        IF (A_IO /= -1.0.OR.LEN_IO /= FIXHD_INTFA(111,JINTF)) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of level dependent constants',A_IO,  &
                        LEN_IO,FIXHD_INTFA(111,JINTF))
          CMESSAGE='IN_INTF : I/O ERROR'
          ICODE=4
          IF (lhook) CALL dr_hook('IN_INTF',zhook_out,zhook_handle)
          RETURN
        END IF

        ENDDO

!L  2.5 Row/Col dependent constants

        IF (INTF_L_VAR_LBC(JINTF)) THEN
! ROW first        
          DO I=1,2
        

            CALL BUFFOUT(NFTOUT,ROWDEPC_INTFA(:,I,JINTF),               &
                         FIXHD_INTFA(116,JINTF),LEN_IO,A_IO)

!           Check for I/O Errors

            IF (A_IO /= -1.0.OR.LEN_IO /= FIXHD_INTFA(116,JINTF)) THEN
! DEPENDS ON: ioerror
              CALL IOERROR('buffer out of row dependent constants',     &
                            A_IO,LEN_IO,FIXHD_INTFA(116,JINTF))
              CMESSAGE='IN_INTF : I/O ERROR'
              ICODE=5
              IF (lhook) CALL dr_hook('IN_INTF',zhook_out,zhook_handle)
              RETURN
            END IF

          ENDDO

! ColDepC        

          DO I=1,2
        

            CALL BUFFOUT(NFTOUT,COLDEPC_INTFA(:,I,JINTF),               &
                         FIXHD_INTFA(121,JINTF),LEN_IO,A_IO)

!           Check for I/O Errors

            IF (A_IO /= -1.0.OR.LEN_IO /= FIXHD_INTFA(121,JINTF)) THEN
! DEPENDS ON: ioerror
              CALL IOERROR('buffer out of col dependent constants',     &
                            A_IO,LEN_IO,FIXHD_INTFA(121,JINTF))
              CMESSAGE='IN_INTF : I/O ERROR'
              ICODE=6
              IF (lhook) CALL dr_hook('IN_INTF',zhook_out,zhook_handle)
              RETURN
            END IF

          ENDDO 
        END IF     ! INTF_L_VAR_LBC(JINTF) 

!L  2.6 Write dummy record to reserve space for PP headers

        LEN_PP =  (INTHD_INTFA(3,JINTF)-1)  * ( LBC_LOOKUPSA-1 )      &
                + LBC_LOOKUPSA                                        &
                + 1

        IF(LEN_PP >  LEN_TOT) THEN
          CMESSAGE='IN_INTF:Insufficient space for PP headers'
          ICODE=5
          IF (lhook) CALL dr_hook('IN_INTF',zhook_out,zhook_handle)
          RETURN
        END IF

        allocate (Dumm1(LEN_PP))

        CALL BUFFOUT(NFTOUT,Dumm1(:),LEN_PP,LEN_IO,A_IO)
        deallocate (Dumm1)

! Check for I/O Errors

        IF(A_IO /= -1.0.OR.LEN_IO /= LEN_PP) THEN
! DEPENDS ON: ioerror
          CALL IOERROR('buffer out of dummy PP headers',A_IO,LEN_IO,    &
                        LEN_PP)
          CMESSAGE='IN_INTF=I/O ERROR'
          ICODE=6
          IF (lhook) CALL dr_hook('IN_INTF',zhook_out,zhook_handle)
          RETURN
        END IF

!    Remainder of headers not used

       ELSE
!L      LNEWBND = False (dataset exists)

!L  3.0 Read in headers
!L      If reinitialised boundary output file to be processed
!L      then open it.
!L
        IF (ft_steps(NFTOUT) >  0) THEN
        STRING=MODEL_FT_UNIT(NFTOUT)
        PPNAME=STRING(18:31)

        CALL FILE_OPEN(NFTOUT,PPNAME,14,1,1,ICODE)
        IF (ICODE /= 0) THEN
          CMESSAGE="IN_INTF: Error opening preassigned boundary file"
          GO TO 9999   !  Return
        ENDIF
!
        ENDIF

!L  3.1  Fixed length header


        CALL BUFFIN(NFTOUT,FIXHD_INTFA(1:,JINTF),LEN_FIXHD,LEN_IO,A_IO)
! Check for I/O Errors
         IF(A_IO /= -1.0.OR.LEN_IO /= LEN_FIXHD) THEN
! DEPENDS ON: ioerror
           CALL IOERROR('buffer in of fixed length header',A_IO,LEN_IO, &
                         LEN_FIXHD)
           CMESSAGE='IN_INTF:I/O ERROR'
           ICODE=7
           IF (lhook) CALL dr_hook('IN_INTF',zhook_out,zhook_handle)
           RETURN
         ENDIF

!L  3.2  Integer header


         CALL BUFFIN(NFTOUT,INTHD_INTFA(1:,JINTF),PP_LEN_INTHD,          &
     &               LEN_IO,A_IO)

! Check for I/O Errors

         IF(A_IO /= -1.0.OR.LEN_IO /= PP_LEN_INTHD) THEN
! DEPENDS ON: ioerror
           CALL IOERROR('buffer in of integer header',A_IO,LEN_IO,      &
                         PP_LEN_INTHD)
           CMESSAGE='IN_INTF:I/O ERROR'
           ICODE=8
           IF (lhook) CALL dr_hook('IN_INTF',zhook_out,zhook_handle)
           RETURN
         END IF

!L  3.3  Real header


         CALL BUFFIN(NFTOUT,REALHD_INTFA(1:,JINTF),PP_LEN_REALHD,        &
                     LEN_IO,A_IO)

! Check for I/O Errors

         IF(A_IO /= -1.0.OR.LEN_IO /= PP_LEN_REALHD) THEN
! DEPENDS ON: ioerror
           CALL IOERROR('buffer in of real header',A_IO,LEN_IO,         &
                         PP_LEN_REALHD)
           CMESSAGE='IN_INTF:I/O ERROR'
           ICODE=9
           IF (lhook) CALL dr_hook('IN_INTF',zhook_out,zhook_handle)
           RETURN
         END IF

!   3.4  Look up Headers

!        Read in just one header ; this will be the last header that
!        corresponds to the last LBC variable before the start of a
!        CRUN. Required to continue disk and start addresses in LOOKUP.

         ntime = ft_lastfield(nftout)

         IF (ntime > 0) THEN

           lookup_start = fixhd_intfa(150,jintf) +                      &
                          fixhd_intfa(151,jintf) *                      &
                          (intf_lookupsa + (ntime-1)*(intf_lookupsa-1))

!          Point to last header

           CALL Setpos (nftout,lookup_start-len1_lookup-1,icode)

           IF (ICode /= 0) THEN
             CMESSAGE = 'Error with SETPOS for LBC Lookup Table'
             ICode    = 10

             CALL EReport (RoutineName, ICode, CMessage)
           END IF

!          Read in one lookup header

           CALL Buffin (nftout,lookup_intfa(1:,intf_lookupsa,jintf),     &
                        len1_lookup,len_io,a_io)

           IF (A_IO /= -1.0 .OR. LEN_IO /= len1_lookup) THEN
! DEPENDS ON: ioerror
             CALL IOERROR('buffer in of lookup header',A_IO,LEN_IO,     &
                         len1_lookup)
             CMESSAGE = 'I/O ERROR with buffin of LBC Lookup Header'
             ICode    = 11

             CALL EReport (RoutineName, ICode, CMessage)
           END IF

         END IF   !  If ntime > 0

!L Reset LNEWBND to true after reading in header information to allow
!L writing of new headers for subsequent reinitialised files.
!
         LNEWBND(JINTF) = .TRUE.
       END IF

      END IF

      END IF         !  LLBOUTim(A_IM)

!L  5.0  End of routine

9999  CONTINUE
      IF (lhook) CALL dr_hook('IN_INTF',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IN_INTF
