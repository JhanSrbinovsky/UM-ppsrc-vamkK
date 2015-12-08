! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Frames

SUBROUTINE frames_writflds(                                               &
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
lbc_size_interp,                                                          &
halo_x,                                                                   &
halo_y,                                                                   &
src_row_len,                                                              &
src_rows,                                                                 &
src_levels,                                                               &
lbc_index_bl,                                                             &
lbc_index_br,                                                             &
nftout,                                                                   &
jintf,                                                                    &
stash_num,                                                                &
field_type,                                                               &
d1_array                                                                  &
)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IO
USE um_parvars
USE frames_mod
USE control_max_sizes
USE io_configuration_mod
USE ereport_mod, ONLY: ereport
USE field_types, ONLY: fld_type_p, fld_type_u, fld_type_v
! This is a logical to determine whether the source data is endgame or new
! dynamics.
USE dynamics_input_mod, ONLY: l_endgame
USE Submodel_Mod

IMPLICIT NONE

! Arguments
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

INTEGER, INTENT(IN) :: lbc_size_interp           ! Size of interp. LBC array
INTEGER, INTENT(IN) :: halo_x                    ! Size of X halo
INTEGER, INTENT(IN) :: halo_y                    ! Size of y halo
INTEGER, INTENT(IN) :: src_row_len               ! Number of columns
INTEGER, INTENT(IN) :: src_rows                  ! Number of rows
INTEGER, INTENT(IN) :: src_levels                ! Number of levels

INTEGER, INTENT(IN) :: lbc_index_bl(lbc_size_interp,2)
INTEGER, INTENT(IN) :: lbc_index_br(lbc_size_interp,2)
                              ! Data points on the input grid required for 
                              ! interpolation to the LBC points (bottom left
                              ! and bottom right; top left and top right are
                              ! these + row_length)
INTEGER, INTENT(IN) :: nftout                      ! Output file unit number
INTEGER, INTENT(IN) :: jintf                     ! Index to interface area
INTEGER, INTENT(IN) :: stash_num                 ! From item_prog in gen_intf_a
INTEGER, INTENT(IN) :: field_type                ! Field Type (p, u, v)
REAL                :: d1_array(1-halo_x:src_row_len+halo_x,              &
                                1-halo_y:src_rows+halo_y,                 &
                                src_levels)

! Common
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



! Local variables

INTEGER :: ij,i,j,k, wrapped_i
INTEGER :: icode 
INTEGER :: ipt
INTEGER :: packed_size
INTEGER :: len_io
INTEGER :: jvar
INTEGER :: start_address
INTEGER :: isc
INTEGER :: frames_row_len
INTEGER :: frames_rows
INTEGER :: ie, iw, js, jn
INTEGER, SAVE :: orog_ie, orog_iw, orog_js, orog_jn
INTEGER, ALLOCATABLE :: frames_mask (:,:)
INTEGER, ALLOCATABLE :: Packed_Frames (:,:)
INTEGER, POINTER :: frames_lookup(:,:)
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
INTEGER :: v_offset

REAL :: bacc
REAL :: a_io
REAL :: mask_value
REAL :: src_zero_lat, src_zero_lon
REAL :: src_delta_lat, src_delta_lon
REAL :: frames_zero_lat, frames_zero_lon
REAL, ALLOCATABLE :: Data_Frames   (:,:,:)
REAL(KIND=jprb)               :: zhook_handle

LOGICAL, ALLOCATABLE :: col_mask(:)

CHARACTER (LEN=1) :: C_Keep_Pack
CHARACTER(LEN=80) cmessage    !  Error Message
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'frames_writflds'

IF (lhook) CALL dr_hook('FRAMES_WRITFLDS',zhook_in,zhook_handle)

frames_lookup => frames_lookups(jintf) % frames_lookup(:,:)
bacc = -99.
isc = NINT(bacc)
packed_size = 0
icode = 0
len_io = 0
jvar = frames_count(jintf)

! Find field in lookup
ipt = ppindex(stash_num,1)
mask_value = RMDI
IF (jvar == 0) THEN
  start_address = frames_fixhds(160,jintf) - 1
ELSE
  start_address = frames_lookup(29,jvar) + frames_lookup(30,jvar)
END IF

IF (l_endgame) THEN
  v_offset = 1
ELSE
  v_offset = -1
END IF

! Create mask - we -1 and +1 in the indices for any error in rotating grids.
! Also we +2 in the rows since we want +1 for any error and +1 for the top
! right and top left grid points.
ALLOCATE ( frames_mask (-2:src_row_len+2, -2:src_rows+2) )
frames_mask(:,:) = 0
ie = -1
iw = src_row_len+2
js = src_rows+2
jn = -1
DO ij = 1, lbc_size_interp

  j = lbc_index_bl(ij,2)
  i = lbc_index_bl(ij,1)
  ! This bit of code wraps a bigger box around the box which is required for 
  ! interpolation.
  ! * - * - * - *
  ! |   |   |   |
  ! * - + - + - *
  ! |   |   |   |
  ! * - @ - @ - *
  ! |   |   |   |
  ! * - * - * - *
  ! 
  ! In the above grid the @ points are the bottom left (bl) and right (br)
  ! points. While the + points are the top left and right points.
  ! The * points are the extra points to try and take into account any rounding
  ! issues when rotating the grid for the winds.
  frames_mask (i,j)     = 1
  frames_mask (i,j+1)   = 1
  frames_mask (i  ,j-1) = 1
  frames_mask (i-1,j-1) = 1
  frames_mask (i-1,j  ) = 1
  frames_mask (i-1,j+1) = 1
  frames_mask (i-1,j+2) = 1
  frames_mask (i  ,j+2) = 1
  ! Set western most point
  iw = MIN(i-1,iw)
  ! Set northern/southern most point
  js = MIN(j-1,js)
  jn = MAX(j+2,jn)
 
  j = lbc_index_br(ij,2)
  i = lbc_index_br(ij,1)
  ! Look above for description of what this is doing.
  frames_mask (i,j)     = 1
  frames_mask (i,j+1)   = 1
  frames_mask (i  ,j-1) = 1
  frames_mask (i+1,j-1) = 1
  frames_mask (i+1,j  ) = 1
  frames_mask (i+1,j+1) = 1
  frames_mask (i+1,j+2) = 1
  frames_mask (i  ,j+2) = 1
  ! Set eastern most point
  ie = MAX(i+1,ie)
  ! Set northern/southern most point
  js = MIN(j-1,js)
  jn = MAX(j+2,jn)
END DO

ALLOCATE(col_mask(src_row_len))
DO i = 1, src_row_len
! Check if column contains mask values
  col_mask(i) = ANY(frames_mask(i,1:src_rows) == 1)
END DO


! We probably have the data wrapped around the GM.
! We want the indices to go from the western boundary to the eastern boundary.
! We can use modular math later to wrap the index back again.
IF (col_mask(1) .AND. col_mask(src_row_len)) THEN
  DO i = 1, src_row_len
    IF (col_mask(i)) THEN
      ie = i + src_row_len
    ELSE
      EXIT
    END IF
  END DO
  DO i = src_row_len, 1, -1
    IF (col_mask(i)) THEN
      iw = i
    ELSE
      EXIT
    END IF
  END DO
END IF

IF (field_type == fld_type_p) THEN
! This is the largest grid.  We need this to be larger than everything else so
! lets make it that little bit larger.
  ie = ie + 1
  iw = iw - 1
  js = MAX(js - 1,1)
  jn = MIN(jn + 1,src_row_len)
END IF
  
! Orography is always peformed first.  We can use this to store the P grid to
! check u and v grid.
IF (stash_num == 33) THEN
  orog_ie = ie
  orog_iw = iw
  orog_js = js
  orog_jn = jn
END IF

IF (field_type == fld_type_p .OR. field_type == fld_type_u) THEN
  IF (orog_ie < ie .OR. orog_iw > iw .OR. &
      orog_js > js .OR. orog_jn < jn) THEN
    ! Call ereport
    CMessage = 'P or U grid is not within orog. grid when cutting out frames.'
    icode     = 1
    CALL Ereport ( RoutineName, icode, CMessage )
  ELSE
    ! We force the grid to be the same.
    ie = orog_ie
    iw = orog_iw
    js = orog_js
    jn = orog_jn
  END IF
ELSE IF (field_type == fld_type_v) THEN
  IF (orog_ie < ie .OR. orog_iw > iw .OR. &
      orog_js > js .OR. orog_jn + v_offset < jn) THEN
    ! Call ereport
    CMessage = 'V grid is not compatible with P grid when cutting out frames.'
    icode    = 1
    CALL Ereport ( RoutineName, icode, CMessage )
  ELSE
    jn = orog_jn + v_offset
  END IF 
  ! We force the grid to be compatible to the P grid.
  ie = orog_ie
  iw = orog_iw
  js = orog_js

END IF

! Size new grid.  Need to use module since we might be wrapped.
frames_row_len = MODULO(ie - iw,src_row_len) + 1
frames_rows    = jn - js + 1

!   -----------------
!   Get whether to keep packing
!   -----------------
IF (keep_pack == -1) THEN
  CALL Fort_Get_Env ('KEEP_PACK', 9, C_Keep_Pack, 1, ICode)
  IF (icode /= 0) THEN
    CMessage = 'Error in finding KEEP_PACK - setting to 0'
    Keep_pack = 0
    icode     = -1
    CALL Ereport ( RoutineName, icode, CMessage )
  ELSE
    READ(C_Keep_Pack, '(I1)') Keep_Pack 
  END IF
END IF

! Cut out minimal

ALLOCATE ( Data_Frames   ( 1:frames_row_len, js:jn, &
                           src_levels ) )
ALLOCATE ( Packed_Frames ( frames_row_len*frames_rows, src_levels ) )

DO k = 1, src_levels
  DO j = js, jn 
    DO i = iw, ie
      wrapped_i = MODULO(i-1,src_row_len)+1
      IF (frames_mask(wrapped_i,j) == 1) THEN
        Data_Frames (i-iw+1,j,k) = d1_array(wrapped_i,j,k)
      ELSE
        Data_Frames (i-iw+1,j,k) = mask_value
      END IF
    END DO
  END DO
END DO

DO k = 1, src_levels
  frames_count(jintf) = frames_count(jintf) + 1
  jvar = frames_count(jintf)

  DO i = 1, 64
    Frames_Lookup (i,jvar) = a_lookup (i,ipt+k-1)
  END DO

  ! Set the grid type correctly to LAM
  Frames_Lookup (17,jvar) = 3

  ! Set address pointers (seem to be the same)  
  Frames_Lookup (29,jvar) = Start_Address
  Frames_Lookup (40,jvar) = Start_Address

  ! Get rid of any extra data.
  Frames_Lookup (20,jvar) = 0
 
  ! Set new dimensions of grid.
  Frames_Lookup(15,jvar) = frames_row_len*frames_rows
  Frames_Lookup(18,jvar) = frames_rows
  Frames_Lookup(19,jvar) = frames_row_len
  
  ! Set new grid coordinates
  src_zero_lat  = TRANSFER(a_lookup(59,ipt+k-1),src_zero_lat)
  src_zero_lon  = TRANSFER(a_lookup(61,ipt+k-1),src_zero_lon)
  src_delta_lat = TRANSFER(a_lookup(60,ipt+k-1),src_delta_lat)
  src_delta_lon = TRANSFER(a_lookup(62,ipt+k-1),src_delta_lon)
  frames_zero_lat = src_zero_lat + (js-1)*src_delta_lat
  frames_zero_lon = src_zero_lon + (iw-1)*src_delta_lon
  Frames_Lookup(59,jvar)  = TRANSFER(frames_zero_lat,Frames_Lookup(59,jvar))
  Frames_Lookup(61,jvar)  = TRANSFER(frames_zero_lon,Frames_Lookup(61,jvar))

!     -----------------------------------
!     Write extracted data to Frames FFs
!     ----------------------------------
  IF (keep_pack == 1 .AND. MOD(Frames_Lookup (21,jvar),10) == 1 ) THEN

    bacc=TRANSFER(Frames_Lookup(51,jvar),bacc)
    isc=NINT(bacc)
    ! Set this to max for now.

! DEPENDS ON: coex
    CALL coex(Data_Frames(:,:,k), Frames_Lookup(15,jvar),                 &
      Packed_Frames(:,k), Frames_Lookup(15,jvar),                         &
      Frames_Lookup(19,jvar), Frames_Lookup(18,jvar),                     &
      Packed_Size, isc, .TRUE., RMDI, 64, icode, cmessage)
      
! Update packed data sizes
    Frames_Lookup(15,jvar) = Packed_Size
    
! The packed size is in 32 bit words.  Therefore need to divde by 2.
    Frames_Lookup (30,jvar) =                                             &
    ((((Packed_Size+1)/2)+io_field_padding-1)/io_field_padding)*io_field_padding
    
    Start_Address = Start_Address + Frames_Lookup (30,jvar)

    CALL SetPos ( nftout, Frames_Lookup(29,jvar), icode )

    IF (ICode /= 0) THEN
      CMessage = 'Error in SetPos for data in Frame FieldsFiles'
      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

    CALL buffout (nftout,                                                 &
                  Packed_Frames(:,k),                                     &
                  Frames_Lookup(15,jvar),                                 &
                  len_io, a_io)

  ELSE IF (keep_pack == 1 .AND. MOD(Frames_Lookup (21,jvar),10) == 2 ) THEN
    ! 32 bit packing
! DEPENDS ON: pack21
    CALL pack21(Frames_Lookup(15, jvar), Data_Frames(:,:,k), Packed_Frames(:,k))
    ! The number of 32 bit words.
    Packed_Size = Frames_Lookup(15, jvar)
! Update packed data sizes
! We want packing....
    
    ! The packed size is in 32 bit words.  Therefore need to divde by 2.
    Frames_Lookup (30,jvar) =                                             &
    ((((Packed_Size+1)/2)+io_field_padding-1)/io_field_padding)*io_field_padding
    
    Start_Address = Start_Address + Frames_Lookup (30,jvar)

    CALL SetPos ( nftout, Frames_Lookup(29,jvar), icode )
    CALL buffout (nftout,                                                 &
                  Packed_Frames(:,k),                                     &
                  Frames_Lookup(15,jvar),                                 &
                  len_io, a_io)
   
  ELSE
    Frames_Lookup (21,jvar) = 0
    
    Frames_Lookup (30,jvar) =                                             &
    ((Frames_Lookup (15,jvar) + io_field_padding - 1)/io_field_padding) *     &
      io_field_padding
    
    Start_Address = Start_Address + Frames_Lookup (30,jvar)
    CALL SetPos ( nftout, Frames_Lookup(29,jvar), icode )

    IF (ICode /= 0) THEN
      CMessage = 'Error in SetPos for data in Frame FieldsFiles'
      CALL Ereport ( RoutineName, icode, CMessage )
    END IF

    CALL buffout (nftout,                                                 &
                  Data_Frames(:,:,k),                                     &
                  Frames_Lookup(15,jvar),                                 &
                  len_io, a_io)

  END IF
  IF (A_IO /= -1.0 .OR. LEN_IO /= Frames_Lookup(15,jvar) ) THEN
! DEPENDS ON: IOERROR
    CALL IOERROR ('buffout of Frames data',                               & 
                   A_IO, LEN_IO, Frames_Lookup(15,jvar) )
    CMESSAGE = 'I/O ERROR with buffout of Frames data'
    ICode    = 11
    CALL EReport (RoutineName, ICode, CMessage)
  END IF

  ! Increment counter for the total number of fields.
END DO

DEALLOCATE ( frames_mask   )
DEALLOCATE ( Data_Frames   )
DEALLOCATE ( Packed_Frames )
DEALLOCATE ( col_mask      )

IF (lhook) CALL dr_hook('FRAMES_WRITFLDS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE frames_writflds
