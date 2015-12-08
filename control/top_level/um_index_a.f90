! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine: UM_INDEX_A---------------------------------------------
!LL
!LL  Purpose: Calculate addresses and lengths of atmosphere super arrays
!LL
!LL  Programming standard: UM Doc Paper 3, version 8.3
!LL
!LL  External documentation: On-line UM document C1 - The top-level
!LL                          dynamic allocation
!LL
!LL  -------------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Top Level

      SUBROUTINE UM_INDEX_A(                                            &
! ARGSZSPA super arrays lengths (atmosphere)
        A_SPDUM_LEN,A_SPPTR_LEN,A_SPCON_LEN,A_SPINF_LEN,A_SPANC_LEN,    &
        A_SPBND_LEN,A_SPSTS_LEN,                                        &
! ARGSZSPA end
     &              ICODE,CMESSAGE)
!
!*----------------------------------------------------------------------

      USE atm_fields_bounds_mod
      USE clmchfcg_scenario_mod, ONLY: nsulpat
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE Control_Max_Sizes
      USE lbc_mod
      IMPLICIT NONE
!
!  Subroutines called
!
!  Local variables
!
      INTEGER ICODE             ! Work - Internal return code
      CHARACTER(LEN=80)  CMESSAGE    ! Work - Internal error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
!  Configuration-dependent sizes for dynamic arrays
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
! Parameters for constants arrays
! CCONSTS start
! Description:
!    This file contains declarations for derived constants within
!   the atmospheric model. Where necessary PARAMETERS are defined to
!   dimension these constants. All constants are placed in the common
!   block CDERIVED, except hardwired constants, e.g. ETA_SPLIT and LENs.
!   file CMAXSIZE must be included first
!
!   The derived constants are calculated in the routine SETCONA.
!
      ! No of cloud types ie low/med/high
      INTEGER, PARAMETER :: NUM_CLOUD_TYPES = 3

      ! derived constants:
      INTEGER :: LOW_BOT_LEVEL      ! Bottom level of lowest cloud type
      INTEGER :: LOW_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: MED_BOT_LEVEL      ! Bottom   "    "  med      "    "
      INTEGER :: MED_TOP_LEVEL      ! Top      "    "   "       "    "
      INTEGER :: HIGH_BOT_LEVEL     ! Bottom   "    "  top      "    "
      INTEGER :: HIGH_TOP_LEVEL     ! Top      "    "   "       "    "

      ! height values to split model levels into l/m/h cloud
      REAL ::    h_split(NUM_CLOUD_TYPES+1)

      LOGICAL :: ELF                ! T if atmosphere model on LAM grid

      ! Constants for dynamics output independent of resolution but
      ! dependent on choice of levels for output.
      REAL :: REQ_THETA_PV_LEVS(MAX_REQ_THPV_LEVS)

      COMMON /CDERIVED/                                                 &
        h_split,LOW_BOT_LEVEL,LOW_TOP_LEVEL,MED_BOT_LEVEL,MED_TOP_LEVEL,&
        HIGH_BOT_LEVEL, HIGH_TOP_LEVEL,ELF,REQ_THETA_PV_LEVS
! CCONSTS end
!  Ancillary file parameters for ancillary length calculations
! CONANC start
!
! To be used in conjunction with file TYPANC, defining dimensions
! of headers from ancillary files. A separate file is needed to
! allow calculation of super array sizes in UM_INDEX.
      INTEGER, PARAMETER :: NANCIL_DATASETSA=99
      INTEGER, PARAMETER :: NANCIL_DATASETS=NANCIL_DATASETSA
! CONANC end
!
!  Super array sizes for dynamic allocation in U_MODEL
!
! TYPSZSPA super arrays lengths (atmosphere)
      INTEGER :: A_SPDUM_LEN
      INTEGER :: A_SPPTR_LEN
      INTEGER :: A_SPCON_LEN
      INTEGER :: A_SPINF_LEN
      INTEGER :: A_SPANC_LEN
      INTEGER :: A_SPBND_LEN
      INTEGER :: A_SPSTS_LEN
! TYPSZSPA end
!
!
!  Addresses of arrays in super arrays.
!
!---------------------Start of SPINDEX--------------------------------
!L
!L --------------- D1    array  -----------------------------------
!L (now includes extra copies of atmos/ocean D1 for MPP)
!L
      INTEGER     IXD1_LEN       ! No. of arrays
      PARAMETER(  IXD1_LEN =     4)
      INTEGER     IXD1           ! Addresses of arrays
      COMMON/    CIXD1/      IXD1(IXD1_LEN)
!L
!L --------------- Dump headers--------------------------
!L
!L ATMOS
      INTEGER   A_IXDUM_LEN       ! No. of arrays
      PARAMETER(A_IXDUM_LEN = 14)
      INTEGER   A_IXDUM           ! Addresses of arrays
      COMMON/  CA_IXDUM/ A_IXDUM(A_IXDUM_LEN)
      INTEGER   A_IXSTS_LEN
      PARAMETER(A_IXSTS_LEN = 12)
      INTEGER   A_IXSTS
      COMMON/  CA_IXSTS/ A_IXSTS(A_IXSTS_LEN)
!L
!L --------------- STASH arrays -----------------------------------
!L
      INTEGER     IXSTS_LEN       ! No. of arrays
      PARAMETER(  IXSTS_LEN = 14)
      INTEGER     IXSTS           ! Addresses of arrays
      COMMON/    CIXSTS/      IXSTS(IXSTS_LEN)
!L
!L --------------- Pointers in D1 array and row/level dependent ---
!L --------------- constants
!L ATMOS
      INTEGER   A_IXPTR_LEN       ! No. of arrays
      PARAMETER(A_IXPTR_LEN = 239)
      INTEGER   A_IXPTR           ! Addresses of arrays
      COMMON/  CA_IXPTR/ A_IXPTR(A_IXPTR_LEN)
!L
!L --------------- Pre-calculated arrays of constants -------------
!L
!L ATMOS
      INTEGER   A_IXCON_LEN       ! No. of arrays
! A_IXCON_LEN is so low b/c the bulk of constants that previously were
! stored in um_index.a are now allocated directly in SETCONA 
      PARAMETER(A_IXCON_LEN = 3)
      INTEGER   A_IXCON           ! Addresses of arrays
      COMMON/  CA_IXCON/ A_IXCON(A_IXCON_LEN)
!L
!L --------------- Headers for output interface datasets (boundary
!L                 conditions out)
!L ATMOS
      INTEGER   A_IXINF_LEN       ! No. of arrays
      PARAMETER(A_IXINF_LEN = 9 )
      INTEGER   A_IXINF           ! Addresses of arrays
      COMMON/  CA_IXINF/ A_IXINF(A_IXINF_LEN)
!L
!L --------------- Headers for ancillary files -------------------
!L
!L ATMOS
      INTEGER   A_IXANC_LEN       ! No. of arrays
      PARAMETER(A_IXANC_LEN = 4)
      INTEGER   A_IXANC           ! Addresses of arrays
      COMMON/  CA_IXANC/ A_IXANC(A_IXANC_LEN)
!L
!L --------------- Headers from input boundary files -------------
!L
!L NOT SUB-MODEL DEPENDENT
      INTEGER     IXBND_LEN       ! No. of arrays
      PARAMETER(  IXBND_LEN = 1)
      INTEGER     IXBND           ! Addresses of arrays
      COMMON/    CIXBND/ IXBND(IXBND_LEN)
!L
!L ATMOS
      INTEGER   A_IXBND_LEN       ! No. of arrays
      PARAMETER(A_IXBND_LEN = 5)
      INTEGER   A_IXBND           ! Addresses of arrays
      COMMON/  CA_IXBND/ A_IXBND(A_IXBND_LEN)
!L
!L --------------- Constant arrays needed for atmosphere-ocean----
!L --------------- coupling
!L
      INTEGER   AO_IXCPL_LEN      ! No. of arrays
      PARAMETER(AO_IXCPL_LEN = 10)
      INTEGER   AO_IXCPL          ! Addresses of arrays
      COMMON/ CAO_IXCPL/ AO_IXCPL(AO_IXCPL_LEN)
!L
!-------------------End of SPINDEX------------------------------------
!
!L----------------------------------------------------------------------
!L
      IF (lhook) CALL dr_hook('UM_INDEX_A',zhook_in,zhook_handle)
      ICODE=0

!L
!L 1.0 Calculate size of atmosphere stash array
!L

! [Note that size of a_ixsts must be the same as o_ixsts and any other
!  submodels, since lower level STASH routines assume the existence of
!  _ixsts(i) even if not used.]

      a_ixsts(1) = 1                              ! zseak_rho
      a_ixsts(2) = a_ixsts(1) + model_levels      ! Ck_rho
      a_ixsts(3) = a_ixsts(2) + model_levels      ! zseak_theta
      a_ixsts(4) = a_ixsts(3) + model_levels  +1  ! Ck_theta
      a_ixsts(5) = a_ixsts(4) + model_levels  +1  ! dummy
      a_ixsts(6) = a_ixsts(5) + 1                 ! dummy
      a_ixsts(7) = a_ixsts(6) + 1                 ! p     pointer in D1
      a_ixsts(8) = a_ixsts(7) + 1                 ! pstar pointer in D1
      a_ixsts(9) = a_ixsts(8) + 1                 ! cos_v_latitude
                                                  !         (no halos)
      a_ixsts(10)= a_ixsts(9) + row_length*n_rows ! cos_theta_latitude
                                                  !         (no halos)
      a_ixsts(11)= a_ixsts(10)+ row_length*rows   ! landsea mask
      a_ixsts(12)= a_ixsts(11)+ row_length*rows   ! land fraction
      a_spsts_len= a_ixsts(12)+ row_length*rows
!     Store A_SPSTS_LEN in TYPSIZE
      LEN_A_SPSTS = A_SPSTS_LEN

!L
!L 1.1   DUMP     super array
!L
!L          super array addresses
      A_IXDUM(1) =1
      A_IXDUM(2) =A_IXDUM(1) + LEN_FIXHD
      A_IXDUM(3) =A_IXDUM(2) + A_LEN_INTHD
      A_IXDUM(4) =A_IXDUM(3) + A_LEN_CFI1+1
      A_IXDUM(5) =A_IXDUM(4) + A_LEN_CFI2+1
      A_IXDUM(6) =A_IXDUM(5) + A_LEN_CFI3+1
      A_IXDUM(7) =A_IXDUM(6) + A_LEN_REALHD
      A_IXDUM(8) =A_IXDUM(7) + A_LEN1_LEVDEPC*A_LEN2_LEVDEPC+1
      A_IXDUM(9) =A_IXDUM(8) + A_LEN1_ROWDEPC*A_LEN2_ROWDEPC+1
      A_IXDUM(10)=A_IXDUM(9) + A_LEN1_COLDEPC*A_LEN2_COLDEPC+1
      A_IXDUM(11)=A_IXDUM(10)+ A_LEN1_FLDDEPC*A_LEN2_FLDDEPC+1
      A_IXDUM(12)=A_IXDUM(11)+ A_LEN_EXTCNST+1
      A_IXDUM(13)=A_IXDUM(12)+ LEN_DUMPHIST+1
      A_IXDUM(14)=A_IXDUM(13)+ LEN1_LOOKUP*A_LEN2_LOOKUP

!L
!L          super array length
      A_SPDUM_LEN  =A_IXDUM(14)+ MPP_LEN1_LOOKUP*A_LEN2_LOOKUP
      A_SPDUM_LEN  =A_SPDUM_LEN -1
!L
!L
!L 1.2   Pointers super array
!L
!L          super array addresses

! Comment at end of each line corresponds to the matching
! pointer in ARGPTRA. eg A_IXPTR(3) in ARTPTRA = jtheta in ARGPTRA

! For each line : A_IXPTR(n+1) = A_IXPTR(n) + n_levs
! where n_levs in the no of levels for pointer n.

      A_IXPTR(1) =1                                          ! ju
      A_IXPTR(2) =A_IXPTR(1) + udims%k_end                   ! jv
      A_IXPTR(3) =A_IXPTR(2) + vdims%k_end                   ! jw
      A_IXPTR(4) =A_IXPTR(3) + wdims%k_end +1                ! jrho
      A_IXPTR(5) =A_IXPTR(4) + pdims%k_end                   ! jtheta
      A_IXPTR(6) =A_IXPTR(5) + tdims%k_end -tdims%k_start +1 ! jq
      A_IXPTR(7) =A_IXPTR(6) + qdims%k_end -qdims%k_start +1 ! jqcl
      A_IXPTR(8) =A_IXPTR(7) + qdims%k_end -qdims%k_start +1 ! jqcf
      A_IXPTR(9) =A_IXPTR(8) + qdims%k_end -qdims%k_start +1 ! jexner_rho_lvls
      A_IXPTR(10)=A_IXPTR(9) + pdims%k_end +1                ! ju_adv
      A_IXPTR(11)=A_IXPTR(10)+ udims%k_end                   ! jv_adv
      A_IXPTR(12)=A_IXPTR(11)+ vdims%k_end                   ! jw_adv
      A_IXPTR(13)=A_IXPTR(12)+ wdims%k_end +1                ! jp
      A_IXPTR(14)=A_IXPTR(13)+ pdims%k_end +1                ! jp_theta_levels
      A_IXPTR(15)=A_IXPTR(14)+ tdims%k_end -tdims%k_start +1 ! jexner_theta_lvls
      A_IXPTR(16)=A_IXPTR(15)+ tdims%k_end -tdims%k_start +1 ! jcca
      A_IXPTR(17)=A_IXPTR(16)+ N_CCA_LEV                      ! jcf_area
      A_IXPTR(18)=A_IXPTR(17)+ qdims%k_end                    ! jcf_bulk
      A_IXPTR(19)=A_IXPTR(18)+ qdims%k_end -qdims%k_start +1 ! jcf_liquid
      A_IXPTR(20)=A_IXPTR(19)+ qdims%k_end -qdims%k_start +1 ! jcf_frozen
      A_IXPTR(21)=A_IXPTR(20)+ qdims%k_end -qdims%k_start +1 ! j_deep_soil_temp
      A_IXPTR(22)=A_IXPTR(21)+ ST_LEVELS                     ! jsmcl
      A_IXPTR(23)=A_IXPTR(22)+ SM_LEVELS                     ! jsthu
      A_IXPTR(24)=A_IXPTR(23)+ SM_LEVELS                     ! jsthf
      A_IXPTR(25)=A_IXPTR(24)+ SM_LEVELS                     ! jsw_incs
      A_IXPTR(26)=A_IXPTR(25)+ MODEL_LEVELS+2                ! jlw_incs
      A_IXPTR(27)=A_IXPTR(26)+ MODEL_LEVELS+1                ! jozone
      A_IXPTR(28)=A_IXPTR(27)+ o3dims2%k_end -o3dims2%k_start +1 ! jtracer
      A_IXPTR(29)=A_IXPTR(28)+ (trdims_xstl%k_end -trdims_xstl%k_start +1)  &
                                                 *(tr_vars+1)! jmurk
      A_IXPTR(30)=A_IXPTR(29)+ tdims%k_end -tdims%k_start +1 ! jmurk_source
      A_IXPTR(31)=A_IXPTR(30)+ tdims%k_end -tdims%k_start +1 ! jso2
      A_IXPTR(32)=A_IXPTR(31)+ tdims%k_end -tdims%k_start +1 ! jdms
      A_IXPTR(33)=A_IXPTR(32)+ tdims%k_end -tdims%k_start +1 ! jso4_aitken
      A_IXPTR(34)=A_IXPTR(33)+ tdims%k_end -tdims%k_start +1 ! jso4_accu
      A_IXPTR(35)=A_IXPTR(34)+ tdims%k_end -tdims%k_start +1 ! jso4_diss
      A_IXPTR(36)=A_IXPTR(35)+ tdims%k_end -tdims%k_start +1 ! jh2o2
      A_IXPTR(37)=A_IXPTR(36)+ tdims%k_end -tdims%k_start +1 ! jnh3
      A_IXPTR(38)=A_IXPTR(37)+ tdims%k_end -tdims%k_start +1 ! jsoot_new
      A_IXPTR(39)=A_IXPTR(38)+ tdims%k_end -tdims%k_start +1 ! jsoot_agd
      A_IXPTR(40)=A_IXPTR(39)+ tdims%k_end -tdims%k_start +1 ! jsoot_cld
      A_IXPTR(41)=A_IXPTR(40)+ tdims%k_end -tdims%k_start +1 ! jso2_natem
      A_IXPTR(42)=A_IXPTR(41)+ tdims%k_end -tdims%k_start +1 ! joh
      A_IXPTR(43)=A_IXPTR(42)+ tdims%k_end -tdims%k_start +1 ! jho2
      A_IXPTR(44)=A_IXPTR(43)+ tdims%k_end -tdims%k_start +1 ! jh2o2_limit
      A_IXPTR(45)=A_IXPTR(44)+ tdims%k_end -tdims%k_start +1 ! jo3_chem
      A_IXPTR(46)=A_IXPTR(45)+ tdims%k_end -tdims%k_start +1 ! jhadcm2_so4
      A_IXPTR(47)=A_IXPTR(46)+ NSULPAT               ! jco2
      A_IXPTR(48)=A_IXPTR(47)+ MODEL_LEVELS          ! juser_mult1
      A_IXPTR(49)=A_IXPTR(48)+ MODEL_LEVELS          ! juser_mult2
      A_IXPTR(50)=A_IXPTR(49)+ MODEL_LEVELS          ! juser_mult3
      A_IXPTR(51)=A_IXPTR(50)+ MODEL_LEVELS          ! juser_mult4
      A_IXPTR(52)=A_IXPTR(51)+ MODEL_LEVELS          ! juser_mult5
      A_IXPTR(53)=A_IXPTR(52)+ MODEL_LEVELS          ! juser_mult6
      A_IXPTR(54)=A_IXPTR(53)+ MODEL_LEVELS          ! juser_mult7
      A_IXPTR(55)=A_IXPTR(54)+ MODEL_LEVELS          ! juser_mult8
      A_IXPTR(56)=A_IXPTR(55)+ MODEL_LEVELS          ! juser_mult9
      A_IXPTR(57)=A_IXPTR(56)+ MODEL_LEVELS          ! juser_mult10
      A_IXPTR(58)=A_IXPTR(57)+ MODEL_LEVELS          ! juser_mult11
      A_IXPTR(59)=A_IXPTR(58)+ MODEL_LEVELS          ! juser_mult12
      A_IXPTR(60)=A_IXPTR(59)+ MODEL_LEVELS          ! juser_mult13
      A_IXPTR(61)=A_IXPTR(60)+ MODEL_LEVELS          ! juser_mult14
      A_IXPTR(62)=A_IXPTR(61)+ MODEL_LEVELS          ! juser_mult15
      A_IXPTR(63)=A_IXPTR(62)+ MODEL_LEVELS          ! juser_mult16
      A_IXPTR(64)=A_IXPTR(63)+ MODEL_LEVELS          ! juser_mult17
      A_IXPTR(65)=A_IXPTR(64)+ MODEL_LEVELS          ! juser_mult18
      A_IXPTR(66)=A_IXPTR(65)+ MODEL_LEVELS          ! juser_mult19
      A_IXPTR(67)=A_IXPTR(66)+ MODEL_LEVELS          ! juser_mult20
      A_IXPTR(68)=A_IXPTR(67)+ MODEL_LEVELS          ! j_orog_lbc
      A_IXPTR(69)=A_IXPTR(68)+ 1                     ! ju_lbc
      A_IXPTR(70)=A_IXPTR(69)+ 1                     ! ju_lbc_tend
      A_IXPTR(71)=A_IXPTR(70)+ 1                     ! jv_lbc
      A_IXPTR(72)=A_IXPTR(71)+ 1                     ! jv_lbc_tend
      A_IXPTR(73)=A_IXPTR(72)+ 1                     ! jw_lbc
      A_IXPTR(74)=A_IXPTR(73)+ 1                     ! jw_lbc_tend
      A_IXPTR(75)=A_IXPTR(74)+ 1                     ! jrho_lbc
      A_IXPTR(76)=A_IXPTR(75)+ 1                     ! jrho_lbc_tend
      A_IXPTR(77)=A_IXPTR(76)+ 1                     ! jtheta_lbc
      A_IXPTR(78)=A_IXPTR(77)+ 1                     ! jtheta_lbc_tend
      A_IXPTR(79)=A_IXPTR(78)+ 1                     ! jq_lbc
      A_IXPTR(80)=A_IXPTR(79)+ 1                     ! jq_lbc_tend
      A_IXPTR(81)=A_IXPTR(80)+ 1                     ! jqcl_lbc
      A_IXPTR(82)=A_IXPTR(81)+ 1                     ! jqcl_lbc_tend
      A_IXPTR(83)=A_IXPTR(82)+ 1                     ! jqcf_lbc
      A_IXPTR(84)=A_IXPTR(83)+ 1                     ! jqcf_lbc_tend
      A_IXPTR(85)=A_IXPTR(84)+ 1                     ! jexner_lbc
      A_IXPTR(86)=A_IXPTR(85)+ 1                     ! jexner_lbc_tend
      A_IXPTR(87)=A_IXPTR(86)+ 1                     ! ju_adv_lbc
      A_IXPTR(88)=A_IXPTR(87)+ 1                     ! ju_adv_lbc_tend
      A_IXPTR(89)=A_IXPTR(88)+ 1                     ! jv_adv_lbc
      A_IXPTR(90)=A_IXPTR(89)+ 1                     ! jv_adv_lbc_tend
      A_IXPTR(91)=A_IXPTR(90)+ 1                     ! jw_adv_lbc
      A_IXPTR(92)=A_IXPTR(91)+ 1                     ! jw_adv_lbc_tend
      A_IXPTR(93)=A_IXPTR(92)+ 1                     ! jtracer_lbc
      A_IXPTR(94)=A_IXPTR(93)+ MAX(TR_LBC_VARS,1)    ! jtracer_lbc_tend
!     tracer_lbc (93) and tracer_lbc_tend (94) have lengths of TR_LBC_VARS
!     which may be zero so use MAX(TR_LBC_VARS,1) for lengths.
      A_IXPTR(95)=A_IXPTR(94)+ MAX(TR_LBC_VARS,1)    ! jtr_ukca_lbc
      A_IXPTR(96)=A_IXPTR(95)+ MAX(TR_LBC_UKCA,1)    ! jtr_ukca_lbc_tend
!     tr_lbc_UKCA (95) and tr_lbc_UKCA_tend (96) have lengths 
!     TR_LBC_UKCA which may be zero so use MAX(TR_LBC_UKCA,1) 
!     for lengths.
      A_IXPTR(97)=A_IXPTR(96) + MAX(TR_LBC_UKCA,1)   ! jso2_lbc
      A_IXPTR(98)=A_IXPTR(97) + 1                    ! jso2_lbc_tend
      A_IXPTR(99)=A_IXPTR(98) + 1                    ! jdms_lbc
      A_IXPTR(100)=A_IXPTR(99) + 1                   ! jdms_lbc_tend
      A_IXPTR(101)=A_IXPTR(100) + 1                  ! jso4_aitken_lbc
      A_IXPTR(102)=A_IXPTR(101) + 1                  ! jso4_aitken_lbc_tend
      A_IXPTR(103)=A_IXPTR(102) + 1                  ! jso4_accu_lbc
      A_IXPTR(104)=A_IXPTR(103) + 1                  ! jso4_accu_lbc_tend
      A_IXPTR(105)=A_IXPTR(104) + 1                  ! jso4_diss_lbc
      A_IXPTR(106)=A_IXPTR(105) + 1                  ! jso4_diss_lbc_tend
      A_IXPTR(107)=A_IXPTR(106) + 1                  ! jnh3_lbc
      A_IXPTR(108)=A_IXPTR(107) + 1                  ! jnh3_lbc_tend
      A_IXPTR(109)=A_IXPTR(108) + 1                  ! jsoot_new_lbc
      A_IXPTR(110)=A_IXPTR(109) + 1                  ! jsoot_new_lbc_tend
      A_IXPTR(111)=A_IXPTR(110) + 1                  ! jsoot_agd_lbc
      A_IXPTR(112)=A_IXPTR(111) + 1                  ! jsoot_agd_lbc_tend
      A_IXPTR(113)=A_IXPTR(112) + 1                  ! jsoot_cld_lbc
      A_IXPTR(114)=A_IXPTR(113) + 1                  ! jsoot_cld_lbc_tend
      A_IXPTR(115)=A_IXPTR(114) + 1                  ! jtppsozone
!     tropopause ozone (115) has TPPS_OZONE_LEVELS which may be zero
!     so use MAX(TPPS_OZONE_LEVELS,1) for length.
      A_IXPTR(116)=A_IXPTR(115) + MAX(TPPS_OZONE_LEVELS,1)      ! jbmass_new
      A_IXPTR(117)=A_IXPTR(116) + tdims%k_end -tdims%k_start +1 ! jbmass_agd
      A_IXPTR(118)=A_IXPTR(117) + tdims%k_end -tdims%k_start +1 ! jbmass_cld      
      A_IXPTR(119)=A_IXPTR(118) + tdims%k_end -tdims%k_start +1 ! jbmass_new_lbc
      A_IXPTR(120)=A_IXPTR(119) + 1                  ! jbmass_new_lbc_tend
      A_IXPTR(121)=A_IXPTR(120) + 1                  ! jbmass_agd_lbc
      A_IXPTR(122)=A_IXPTR(121) + 1                  ! jbmass_agd_lbc_tend
      A_IXPTR(123)=A_IXPTR(122) + 1                  ! jbmass_cld_lbc
      A_IXPTR(124)=A_IXPTR(123) + 1                  ! jbmass_cld_lbc_tend
      A_IXPTR(125)=A_IXPTR(124) + 1                  ! jqcf2
      A_IXPTR(126)=A_IXPTR(125) + qdims%k_end -qdims%k_start +1 ! jqrain
      A_IXPTR(127)=A_IXPTR(126) + qdims%k_end -qdims%k_start +1 ! jqgraup
      A_IXPTR(128)=A_IXPTR(127) + qdims%k_end -qdims%k_start +1  ! jqcf2_lbc
      A_IXPTR(129)=A_IXPTR(128) + 1                  ! jqcf2_lbc_tend
      A_IXPTR(130)=A_IXPTR(129) + 1                  ! jqrain_lbc
      A_IXPTR(131)=A_IXPTR(130) + 1                  ! jqrain_lbc_tend
      A_IXPTR(132)=A_IXPTR(131) + 1                  ! jqgraup_lbc
      A_IXPTR(133)=A_IXPTR(132) + 1                  ! jqgraup_lbc_tend
      A_IXPTR(134)=A_IXPTR(133) + 1                  ! jdust_div1
      A_IXPTR(135)=A_IXPTR(134) + tdims%k_end -tdims%k_start +1 ! jdust_div2
      A_IXPTR(136)=A_IXPTR(135) + tdims%k_end -tdims%k_start +1 ! jdust_div3
      A_IXPTR(137)=A_IXPTR(136) + tdims%k_end -tdims%k_start +1 ! jdust_div4
      A_IXPTR(138)=A_IXPTR(137) + tdims%k_end -tdims%k_start +1 ! jdust_div5
      A_IXPTR(139)=A_IXPTR(138) + tdims%k_end -tdims%k_start +1 ! jdust_div6      
      A_IXPTR(140)=A_IXPTR(139) + tdims%k_end -tdims%k_start +1 ! jdust_div1_lbc
      A_IXPTR(141)=A_IXPTR(140) + 1                  ! jdust_div1_lbc_tend
      A_IXPTR(142)=A_IXPTR(141) + 1                  ! jdust_div2_lbc
      A_IXPTR(143)=A_IXPTR(142) + 1                  ! jdust_div2_lbc_tend
      A_IXPTR(144)=A_IXPTR(143) + 1                  ! jdust_div3_lbc
      A_IXPTR(145)=A_IXPTR(144) + 1                  ! jdust_div3_lbc_tend
      A_IXPTR(146)=A_IXPTR(145) + 1                  ! jdust_div4_lbc
      A_IXPTR(147)=A_IXPTR(146) + 1                  ! jdust_div4_lbc_tend
      A_IXPTR(148)=A_IXPTR(147) + 1                  ! jdust_div5_lbc
      A_IXPTR(149)=A_IXPTR(148) + 1                  ! jdust_div5_lbc_tend
      A_IXPTR(150)=A_IXPTR(149) + 1                  ! jdust_div6_lbc
      A_IXPTR(151)=A_IXPTR(150) + 1                  ! jdust_div6_lbc_tend
      A_IXPTR(152)=A_IXPTR(151) + 1                  ! jcf_bulk_lbc
      A_IXPTR(153)=A_IXPTR(152) + 1                  ! jcf_bulk_lbc_tend
      A_IXPTR(154)=A_IXPTR(153) + 1                  ! jcf_liquid_lbc
      A_IXPTR(155)=A_IXPTR(154) + 1                  ! jcf_liquid_lbc_tend
      A_IXPTR(156)=A_IXPTR(155) + 1                  ! jcf_frozen_lbc
      A_IXPTR(157)=A_IXPTR(156) + 1                  ! jcf_frozen_lbc_tend
      A_IXPTR(158)=A_IXPTR(157) + 1                  ! jdirpar
      A_IXPTR(159)=A_IXPTR(158) + 1                  ! jtr_ukca
      A_IXPTR(160)=A_IXPTR(159) + (tr_levels -tdims%k_start +1)  &
                                                   *(tr_ukca+1) ! jmurk_lbc
      A_IXPTR(161)=A_IXPTR(160) + 1                  ! jmurk_lbc_tend
      A_IXPTR(162)=A_IXPTR(161) + 1                  ! jOH_UKCA
      A_IXPTR(163)=A_IXPTR(162) + tdims%k_end -tdims%k_start +1 ! jHO2_UKCA
      A_IXPTR(164)=A_IXPTR(163) + tdims%k_end -tdims%k_start +1 ! jH2O2_UKCA
      A_IXPTR(165)=A_IXPTR(164) + tdims%k_end -tdims%k_start +1 ! jO3_UKCA
      A_IXPTR(166)=A_IXPTR(165) + tdims%k_end -tdims%k_start +1 ! jlcbase
      A_IXPTR(167)=A_IXPTR(166) + 1                  ! jccw_rad 
      A_IXPTR(168)=A_IXPTR(167) + qdims%k_end        ! jozone_tracer
      A_IXPTR(169)=A_IXPTR(168) + tdims%k_end -tdims%k_start +1 ! jO3_prod_loss
      A_IXPTR(170)=A_IXPTR(169) + tdims%k_end -tdims%k_start +1 ! jO3_P_L_VMR
      A_IXPTR(171)=A_IXPTR(170) + tdims%k_end -tdims%k_start +1 ! jO3_VMR
      A_IXPTR(172)=A_IXPTR(171) + tdims%k_end -tdims%k_start +1 ! jO3_P_L_temp
      A_IXPTR(173)=A_IXPTR(172) + tdims%k_end -tdims%k_start +1 ! jO3_temp
      A_IXPTR(174)=A_IXPTR(173) + tdims%k_end -tdims%k_start +1 ! jO3_P_L_colO3
      A_IXPTR(175)=A_IXPTR(174) + tdims%k_end -tdims%k_start +1 ! jO3_colO3
      A_IXPTR(176)=A_IXPTR(175) + tdims%k_end -tdims%k_start +1 ! jarclbiog_bg
      A_IXPTR(177)=A_IXPTR(176) + tdims%k_end -tdims%k_start +1 ! jarclbiom_fr
      A_IXPTR(178)=A_IXPTR(177) + tdims%k_end -tdims%k_start +1 ! jarclbiom_ag
      A_IXPTR(179)=A_IXPTR(178) + tdims%k_end -tdims%k_start +1 ! jarclbiom_ic
      A_IXPTR(180)=A_IXPTR(179) + tdims%k_end -tdims%k_start +1 ! jarclblck_fr
      A_IXPTR(181)=A_IXPTR(180) + tdims%k_end -tdims%k_start +1 ! jarclblck_ag
      A_IXPTR(182)=A_IXPTR(181) + tdims%k_end -tdims%k_start +1 ! jarclsslt_fi
      A_IXPTR(183)=A_IXPTR(182) + tdims%k_end -tdims%k_start +1 ! jarclsslt_je
      A_IXPTR(184)=A_IXPTR(183) + tdims%k_end -tdims%k_start +1 ! jarclsulp_ac
      A_IXPTR(185)=A_IXPTR(184) + tdims%k_end -tdims%k_start +1 ! jarclsulp_ak
      A_IXPTR(186)=A_IXPTR(185) + tdims%k_end -tdims%k_start +1 ! jarclsulp_di
      A_IXPTR(187)=A_IXPTR(186) + tdims%k_end -tdims%k_start +1 ! jarcldust_b1
      A_IXPTR(188)=A_IXPTR(187) + tdims%k_end -tdims%k_start +1 ! jarcldust_b2
      A_IXPTR(189)=A_IXPTR(188) + tdims%k_end -tdims%k_start +1 ! jarcldust_b3
      A_IXPTR(190)=A_IXPTR(189) + tdims%k_end -tdims%k_start +1 ! jarcldust_b4
      A_IXPTR(191)=A_IXPTR(190) + tdims%k_end -tdims%k_start +1 ! jarcldust_b5
      A_IXPTR(192)=A_IXPTR(191) + tdims%k_end -tdims%k_start +1 ! jarcldust_b6
      A_IXPTR(193)=A_IXPTR(192) + tdims%k_end -tdims%k_start +1 ! jarclocff_fr
      A_IXPTR(194)=A_IXPTR(193) + tdims%k_end -tdims%k_start +1 ! jarclocff_ag
      A_IXPTR(195)=A_IXPTR(194) + tdims%k_end -tdims%k_start +1 ! jarclocff_ic
      A_IXPTR(196)=A_IXPTR(195) + tdims%k_end -tdims%k_start +1 ! jarcldlta_dl
      A_IXPTR(197)=A_IXPTR(196) + tdims%k_end -tdims%k_start +1 ! jocff_new
      A_IXPTR(198)=A_IXPTR(197) + tdims%k_end -tdims%k_start +1 ! jocff_agd
      A_IXPTR(199)=A_IXPTR(198) + tdims%k_end -tdims%k_start +1 ! jocff_cld
      A_IXPTR(200)=A_IXPTR(199) + tdims%k_end -tdims%k_start +1 ! jocff_new_lbc
      A_IXPTR(201)=A_IXPTR(200) + 1                  ! jocff_new_lbc_tend
      A_IXPTR(202)=A_IXPTR(201) + 1                  ! jocff_agd_lbc
      A_IXPTR(203)=A_IXPTR(202) + 1                  ! jocff_agd_lbc_tend
      A_IXPTR(204)=A_IXPTR(203) + 1                  ! jocff_cld_lbc
      A_IXPTR(205)=A_IXPTR(204) + 1                  ! jocff_cld_lbc_tend
      A_IXPTR(206)=A_IXPTR(205) + 1                  ! jHNO3_UKCA
      A_IXPTR(207)=A_IXPTR(206) + tdims%k_end -tdims%k_start +1 ! jnitr_acc
      A_IXPTR(208)=A_IXPTR(207) + tdims%k_end -tdims%k_start +1 ! jnitr_diss
      A_IXPTR(209)=A_IXPTR(208) + tdims%k_end -tdims%k_start +1 ! jnitr_acc_lbc
      A_IXPTR(210)=A_IXPTR(209) + tdims_s%k_end-tdims_s%k_start+1 ! jnitr_acc_lbc_tend
      A_IXPTR(211)=A_IXPTR(210) + tdims_s%k_end-tdims_s%k_start+1 ! jnitr_diss_lbc
      A_IXPTR(212)=A_IXPTR(211) + tdims_s%k_end-tdims_s%k_start+1 ! jnitr_diss_lbc_tend
      A_IXPTR(213)=A_IXPTR(212) + tdims_s%k_end-tdims_s%k_start+1 ! je_trb      
      A_IXPTR(214)=A_IXPTR(213) + tdims%k_end -tdims%k_start +1 ! jtsq_trb      
      A_IXPTR(215)=A_IXPTR(214) + tdims%k_end -tdims%k_start +1 ! jqsq_trb      
      A_IXPTR(216)=A_IXPTR(215) + tdims%k_end -tdims%k_start +1 ! jcov_trb      
      A_IXPTR(217)=A_IXPTR(216) + tdims%k_end -tdims%k_start +1 ! jzhpar_shcu   
      A_IXPTR(218)=A_IXPTR(217) + 1 ! JDRYRHO
      A_IXPTR(219)=A_IXPTR(218) + pdims_s%k_end-pdims_s%k_start+1 ! JETADOT 
      A_IXPTR(220)=A_IXPTR(219) + wdims_s%k_end-wdims_s%k_start+1 ! JTHETAV 
      A_IXPTR(221)=A_IXPTR(220) + tdims_s%k_end-tdims_s%k_start+1 ! JPSIWS
      A_IXPTR(222)=A_IXPTR(221) + 1                               ! JPSIWL
      A_IXPTR(223)=A_IXPTR(222) + 1                               ! JMV
      A_IXPTR(224)=A_IXPTR(223) + qdims_s%k_end-qdims_s%k_start+1 ! JMCL
      A_IXPTR(225)=A_IXPTR(224) + qdims_s%k_end-qdims_s%k_start+1 ! JMCF
      A_IXPTR(226)=A_IXPTR(225) + qdims_s%k_end-qdims_s%k_start+1 ! JMCF2
      A_IXPTR(227)=A_IXPTR(226) + qdims_s%k_end-qdims_s%k_start+1 ! JMRAIN
      A_IXPTR(228)=A_IXPTR(227) + qdims_s%k_end-qdims_s%k_start+1 ! JMGRAUP
      A_IXPTR(229)=A_IXPTR(228) + qdims_s%k_end-qdims_s%k_start+1 ! 
      !CABLE
      A_IXPTR(230)=A_IXPTR(229) + SM_LEVELS          ! jtsoil_tile 
      A_IXPTR(231)=A_IXPTR(230) + SM_LEVELS          ! jsmcl_tile 
      A_IXPTR(232)=A_IXPTR(231) + SM_LEVELS          ! jsthf_tile 
      A_IXPTR(233)=A_IXPTR(232) + 3                  ! jsnow_depth3l 
      A_IXPTR(234)=A_IXPTR(233) + 3                  ! jsnow_mass3l 
      A_IXPTR(235)=A_IXPTR(234) + 3                  ! jsnow_tmp3l 
      A_IXPTR(236)=A_IXPTR(235) + 3                  ! jsnow_rho3l 
      A_IXPTR(237)=A_IXPTR(236) + 1                  ! jsnow_rho1l  
      A_IXPTR(238)=A_IXPTR(237) + 1                  ! jsnow_age 
      A_IXPTR(239)=A_IXPTR(238) + 1                  ! jsnow_flg3l 
!     Super array length. If length of last variable could be zero, use
!     max(nlevs,1) to avoid out of bounds warning messages

      A_SPPTR_LEN  =A_IXPTR(239) +1   
      A_SPPTR_LEN  =A_SPPTR_LEN -1
! 
!
!  1.3   Derived constants super array
!
!           super array addresses
!   The size of the increment on each line is actually the size of
!   the array referenced by the previous line,
!   e.g. (model_levels+1) is the size of eta_theta_levels.
      A_IXCON(1) = 1                         ! land_index
      A_IXCON(2) =A_IXCON(1) + land_field  ! land_ice_index
      A_IXCON(3) =A_IXCON(2) + land_field  ! soil_index
! other length info is now dealt with by ALLOCATE statements in SETCONA

      A_SPCON_LEN = A_IXCON(3) + land_field  ! = 3*land_field + 1
! Subtract 1 unless A_SPCON_LEN = 1
      IF (A_SPCON_LEN > 1) A_SPCON_LEN = A_SPCON_LEN -1
! 
!
!  1.4   Interface output (boundary conditions) super array
!
!           super array addresses
      ! FIXHD_INTFA
      a_ixinf(1) = 1
      ! INTHD_INTFA
      a_ixinf(2) = a_ixinf(1) + len_fixhd*n_intf_a
      ! LOOKUP_INTFA
      a_ixinf(3) = a_ixinf(2) + pp_len_inthd*n_intf_a
      ! REALHD_INTFA
      a_ixinf(4) = a_ixinf(3) + len1_lookup*intf_lookupsa*n_intf_a
      ! LEVDEPC_INTFA
      a_ixinf(5) = a_ixinf(4) + pp_len_realhd*n_intf_a
      ! ROWDEPC_INTFA
      a_ixinf(6) = a_ixinf(5) +                                        &
                  max_intf_model_levels * intf_len2_levdepc * n_intf_a
      ! coldepc_intfa
      a_ixinf(7) = a_ixinf(6) +                                        &
                   max(max_lbcrows * intf_len2_rowdepc * n_intf_a,1)
      ! lbc_eta_theta
      a_ixinf(8) = a_ixinf(7) +                                        &
                   max(max_lbcrow_length * intf_len2_coldepc * n_intf_a,1)
      ! lbc_eta_rho
      a_ixinf(9) = a_ixinf(8) + (max_intf_model_levels+1)*n_intf_a

!           super array length
      a_spinf_len = a_ixinf(9) + (max_intf_model_levels)*n_intf_a

      a_spinf_len = a_spinf_len -1

! 
!  1.5   Ancillary file super array
!
!           super array addresses
      A_IXANC(1) =1
      A_IXANC(2) =A_IXANC(1) + LEN_FIXHD*NANCIL_DATASETSA
      A_IXANC(3) =A_IXANC(2) + A_LEN_INTHD*NANCIL_DATASETSA
      A_IXANC(4) =A_IXANC(3) + LEN1_LOOKUP*NANCIL_LOOKUPSA
! 
!           super array length
      A_SPANC_LEN  =A_IXANC(4)+ A_LEN_REALHD*NANCIL_DATASETSA
      A_SPANC_LEN  =A_SPANC_LEN -1
! 
!
!  1.6   Input boundary constants super array
!
!           super array addresses
      A_IXBND(1) =1
      A_IXBND(2) =A_IXBND(1) + LEN_FIXHD
      A_IXBND(3) =A_IXBND(2) + A_LEN_INTHD
      A_IXBND(4) =A_IXBND(3) + LEN1_LOOKUP*RIM_LOOKUPSA
      A_IXBND(5) =A_IXBND(4) + LEN1_LBC_COMP_LOOKUP*BOUND_LOOKUPSA
! 
!           super array length
      A_SPBND_LEN  =A_IXBND(5)+ A_LEN_REALHD
      A_SPBND_LEN  =A_SPBND_LEN -1
! 
! ----------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook('UM_INDEX_A',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UM_INDEX_A
