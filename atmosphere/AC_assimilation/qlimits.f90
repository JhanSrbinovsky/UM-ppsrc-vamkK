! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Apply bounds to tropospheric and non-tropospheric humidities

MODULE qlimits_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE QLimits ( L_Diags,                & ! in
                     L_RmNonTropQIncs,       & ! in
                     trop_min_RH,            & ! in
                     nonTrop_max_q,          & ! in
                     nonTrop_min_q,          & ! in
                     nonTrop_max_RH,         & ! in
                     trop_min_p,             & ! in
                     trop_max_PV,            & ! in
                     nonTrop_max_p,          & ! in
                     q_orig,                 & ! in
                     pv_at_theta,            & ! in
                     theta,                  & ! in
                     p,                      & ! in
                     p_theta_levels,         & ! in
                     exner_theta_levels,     & ! in
                     q,                      & ! inout
                     qCL,                    & ! inout
                     qCF,                    & ! inout
                     area_cloud_fraction,    & ! inout
                     bulk_cloud_fraction,    & ! inout
                     cloud_fraction_liquid,  & ! inout
                     cloud_fraction_frozen )   ! inout

! Description:
!
!   A point is considered tropospheric if one of the following holds:
!
!   1. Its pressure value is greater than nonTrop_max_p.
!   2. Its potential vorticity (PV) value is no greater than trop_max_PV, and
!      its pressure value is at least trop_min_p.
!
!   If neither of these conditions holds, the point is considered
!   non-tropospheric.
!
!   For tropospheric points, we impose a lower limit trop_min_RH on the
!   relative humidities:
!
!     RH >= trop_min_RH
!
!   For non-tropospheric points, the specific humidity is limited to the range
!   [nonTrop_min_q, nonTrop_max_q], and the relative humidity is capped at
!   nonTrop_max_RH:
!
!     nonTrop_min_q <= q <= nonTrop_max_q
!     RH <= nonTrop_max_RH
!
!   An alternative for non-tropospheric points - activated via the switch
!   L_RmNonTropQIncs - is to restore the specific humidities to their original
!   values, as passed in via q_orig.
!
!   In all circumstances cloud amounts are set to zero outside the troposphere.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: AC Assimilation
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!
! Declarations:

USE atm_fields_bounds_mod

USE level_heights_mod, ONLY : &
    r_theta_levels

USE trignometric_mod, ONLY :  &
    cos_theta_latitude

USE earth_constants_mod, ONLY: g

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE UM_ParVars
USE Control_Max_Sizes
IMPLICIT NONE

! DEPENDS ON: qsat

! Common blocks:

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
!-----------------Start of TYPCONA------------------------------------

! Constants for routines independent of resolution.
! Requires Use Control_Max_Sizes for MAX_REQ_THPV_LEVS in cconsts.h
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

! typcona.h originally contained constants for the atmosphere.
! Almost all of these constants have moved to a set of modules:
! LEVEL_HEIGHTS_MOD, TRIGNOMETRIC_MOD, DYN_CORIOLIS_MOD, DYN_VAR_RES_MOD,
! DIFF_COEFF_MOD, RAD_MASK_TROP_MOD, ROT_COEFF_MOD

! The following common block did not correspond to constants specified in 
! argcona.h, so it remains here, even though argcona.h has been deleted.
! cderv_trig and CDERIVED in cconsts.h should be moved to modules too.

      ! Trigonometric co-ordinates in radians
      REAL:: Delta_lambda       ! EW (x) grid spacing in radians
      REAL:: Delta_phi          ! NS (y) grid spacing in radians
      REAL:: Base_phi           ! Lat of first theta point in radians
      REAL:: Base_lambda        ! Long of first theta point in radians
      REAL:: lat_rot_NP         ! Real lat of 'pseudo' N pole in radians
      REAL:: long_rot_NP        ! Real long of 'pseudo' N pole in radians

      COMMON/cderv_trig/                                                &
     &  Delta_lambda,Delta_phi,Base_phi,Base_lambda,                    &
     &  lat_rot_NP,long_rot_NP
!-------------End of TYPCONA---------------------------------------

! Subroutine arguments:

LOGICAL, INTENT(IN)    :: L_Diags          ! Write out diagnostics?
LOGICAL, INTENT(IN)    :: L_RmNonTropQIncs ! Remove non-trop q increments?

REAL,    INTENT(IN)    :: trop_min_RH    ! Lower limit to apply to trop RH
REAL,    INTENT(IN)    :: nonTrop_max_q  ! Upper limit to apply to non-trop q
REAL,    INTENT(IN)    :: nonTrop_min_q  ! Lower limit to apply to non-trop q
REAL,    INTENT(IN)    :: nonTrop_max_RH ! Upper limit to apply to non-trop RH

REAL,    INTENT(IN)    :: trop_min_p     ! Minimum tropospheric pressure
REAL,    INTENT(IN)    :: trop_max_PV    ! Maximum tropospheric ABS(PV)
REAL,    INTENT(IN)    :: nonTrop_max_p  ! Maximum non-trop     pressure

REAL,    INTENT(IN)    :: q_orig      ( qdims%i_start : qdims%i_end,  &
                                        qdims%j_start : qdims%j_end,  &
                                                    1 : qdims%k_end )

REAL,    INTENT(IN)    :: pv_at_theta ( tdims%i_start : tdims%i_end,  &
                                        tdims%j_start : tdims%j_end,  &
                                                    1 : tdims%k_end )

REAL,    INTENT(IN)    ::                                   &
  theta                 ( tdims_s%i_start : tdims_s%i_end,  &
                          tdims_s%j_start : tdims_s%j_end,  &
                          tdims_s%k_start : tdims_s%k_end )
      
REAL,    INTENT(IN)    ::                                   &
  p                     ( pdims_s%i_start : pdims_s%i_end,  &
                          pdims_s%j_start : pdims_s%j_end,  &
                          pdims_s%k_start : pdims_s%k_end + 1 )
      
REAL,    INTENT(IN)    ::                                   &
  p_theta_levels        ( tdims_s%i_start : tdims_s%i_end,  &
                          tdims_s%j_start : tdims_s%j_end,  &
                          tdims_s%k_start : tdims_s%k_end )
      
REAL,    INTENT(IN)    ::                                   &
  exner_theta_levels    ( tdims_s%i_start : tdims_s%i_end,  &
                          tdims_s%j_start : tdims_s%j_end,  &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT) ::                                   &
  q                     ( qdims_l%i_start : qdims_l%i_end,  &
                          qdims_l%j_start : qdims_l%j_end,  &
                          qdims_l%k_start : qdims_l%k_end )
      
REAL,    INTENT(INOUT) ::                                   &
  qCL                   ( qdims_l%i_start : qdims_l%i_end,  &
                          qdims_l%j_start : qdims_l%j_end,  &
                          qdims_l%k_start : qdims_l%k_end )
      
REAL,    INTENT(INOUT) ::                                   &
  qCF                   ( qdims_l%i_start : qdims_l%i_end,  &
                          qdims_l%j_start : qdims_l%j_end,  &
                          qdims_l%k_start : qdims_l%k_end )
      
REAL,    INTENT(INOUT) :: area_cloud_fraction (qdims%i_start:qdims%i_end,  &
                                               qdims%j_start:qdims%j_end,  &
                                                           1:qdims%k_end )
      
REAL,    INTENT(INOUT) ::                                   &
  bulk_cloud_fraction   ( qdims_l%i_start : qdims_l%i_end,  &
                          qdims_l%j_start : qdims_l%j_end,  &
                          qdims_l%k_start : qdims_l%k_end )
      
REAL,    INTENT(INOUT) ::                                   &
  cloud_fraction_liquid ( qdims_l%i_start : qdims_l%i_end,  &
                          qdims_l%j_start : qdims_l%j_end,  &
                          qdims_l%k_start : qdims_l%k_end )
      
REAL,    INTENT(INOUT) ::                                   &
  cloud_fraction_frozen ( qdims_l%i_start : qdims_l%i_end,  &
                          qdims_l%j_start : qdims_l%j_end,  &
                          qdims_l%k_start : qdims_l%k_end )

! Local variables:

INTEGER :: i, j, k, ICode

LOGICAL :: InTrop(tdims%i_start : tdims%i_end,  &
                  tdims%j_start : tdims%j_end)

LOGICAL :: Strat

REAL :: t    (tdims%i_start : tdims%i_end, &
              tdims%j_start : tdims%j_end)
REAL :: p_tmp(pdims%i_start : pdims%i_end, &
              pdims%j_start : pdims%j_end)
REAL :: q_sat(qdims%i_start : qdims%i_end, &
              qdims%j_start : qdims%j_end)
REAL :: q_new(qdims%i_start : qdims%i_end, &
              qdims%j_start : qdims%j_end)

REAL :: gb_mass                  ! Total air mass in gridbox
REAL :: total_mass               ! Total air mass
REAL :: trop_mass                ! Tropospheric air mass
REAL :: trop_frac                ! Fraction of mass in troposphere
REAL :: trop_vmass_before        ! Tropospheric     vapour mass on entry
REAL :: trop_vmass_after         ! Tropospheric     vapour mass on exit
REAL :: nonTrop_vmass_before     ! Non-tropospheric vapour mass on entry
REAL :: nonTrop_vmass_after      ! Non-tropospheric vapour mass on exit
REAL :: trop_vmass_percChange    ! Percentage change to     trop vapour mass
REAL :: nonTrop_vmass_percChange ! Percentage change to non-trop vapour mass

REAL :: stats(6, 1 : qdims%k_end)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!- End of header ---------------------------------------------------------------

!-------------------------------------------------------------------------------
! [1]: Apply tropospheric and non-tropospheric humidity limits, and if required
!      calculate local statistics for points on this PE.
!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('QLIMITS',zhook_in,zhook_handle)
DO k = 1, qdims%k_end

  ! Get temperature and pressure.
  t    (:,:) = exner_theta_levels(tdims%i_start : tdims%i_end, &
                                  tdims%j_start : tdims%j_end,k) &
             * theta             (tdims%i_start : tdims%i_end, &
                                  tdims%j_start : tdims%j_end,k)
  p_tmp(:,:) = p_theta_levels    (tdims%i_start : tdims%i_end, &
                                  tdims%j_start : tdims%j_end,k)

  ! Calculate saturated specific humidity.
  CALL QSAT (q_sat, t, p_tmp, row_length*rows)

  ! Apply limits:
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

      ! Is this point in the troposphere?
      InTrop(i,j) = (   p_theta_levels (i,j,k)  >  nonTrop_max_p .OR.  &
                      ( p_theta_levels (i,j,k)  >= trop_min_p    .AND. &
                        ABS(pv_at_theta(i,j,k)) <= trop_max_PV         &
                      )                                                &
                    )

      IF (InTrop(i,j)) THEN
        q_new(i,j) = MAX(q(i,j,k), q_sat(i,j)*trop_min_RH)
      ELSE
        ! Remove cloud:
        qCL(i,j,k) = 0.0
        qCF(i,j,k) = 0.0
        cloud_fraction_liquid(i,j,k) = 0.0
        cloud_fraction_frozen(i,j,k) = 0.0
        area_cloud_fraction  (i,j,k) = 0.0
        bulk_cloud_fraction  (i,j,k) = 0.0
        ! Adjust q:
        IF (L_RmNonTropQIncs) THEN
          q_new(i,j) = q_orig(i,j,k)
        ELSE
          q_new(i,j) = MAX(q    (i,j,k), nonTrop_min_q)
          q_new(i,j) = MIN(q_new(i,j),   nonTrop_max_q)
          q_new(i,j) = MIN(q_new(i,j),   q_sat(i,j)*nonTrop_max_RH)
        END IF
      END IF

    END DO ! i
  END DO ! j

  ! Calculate local statistics:
  IF (L_Diags) THEN

    ! Initialise stats.
    total_mass           = 0.0
    trop_mass            = 0.0
    trop_vmass_before    = 0.0
    trop_vmass_after     = 0.0
    nonTrop_vmass_before = 0.0
    nonTrop_vmass_after  = 0.0

  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

        ! Approximate mass in gridbox:
        gb_mass = ( p(i,j,k) - p(i,j,k+1) )  &
                  * cos_theta_latitude(i,j)  &
                  * delta_lambda * delta_phi &
                  * r_theta_levels(i,j,k)**2 &
                  / g

        total_mass = total_mass + gb_mass

        IF (InTrop(i,j)) THEN
          trop_mass            = trop_mass            + gb_mass
          trop_vmass_before    = trop_vmass_before    + gb_mass * q    (i,j,k)
          trop_vmass_after     = trop_vmass_after     + gb_mass * q_new(i,j)
        ELSE
          nonTrop_vmass_before = nonTrop_vmass_before + gb_mass * q    (i,j,k)
          nonTrop_vmass_after  = nonTrop_vmass_after  + gb_mass * q_new(i,j)
        END IF

      END DO
    END DO

    stats(1,k) = total_mass
    stats(2,k) = trop_mass
    stats(3,k) = trop_vmass_before
    stats(4,k) = trop_vmass_after
    stats(5,k) = nonTrop_vmass_before
    stats(6,k) = nonTrop_vmass_after

  END IF ! (L_Diags)

  ! Update q:
  q(qdims%i_start : qdims%i_end, &
    qdims%j_start : qdims%j_end,k) = q_new(:,:)

END DO ! k

!-------------------------------------------------------------------------------
! [2]: Obtain and write global diagnostics.
!-------------------------------------------------------------------------------

IF (L_Diags) THEN

  ! Convert local stats into global stats:
  CALL GCG_RSUMR (wet_levels*6, gc_all_proc_group, ICode, stats)

  WRITE (6,*) ''
  WRITE (6,*) '<><><><><><><><><><><><><><>'
  WRITE (6,*) 'Start of QLimits diagnostics'
  WRITE (6,*) ''
  WRITE (6,'(A,ES10.4)') &
    ' Lower limit to apply to trop RH:     ', trop_min_RH
  WRITE (6,'(A,ES10.4)') &
    ' Upper limit to apply to non-trop q:  ', nonTrop_max_q
  WRITE (6,'(A,ES10.4)') &
    ' Lower limit to apply to non-trop q:  ', nonTrop_min_q
  WRITE (6,'(A,ES10.4)') &
    ' Upper limit to apply to non-trop RH: ', nonTrop_max_RH
  WRITE (6,'(A,ES10.4)') &
    ' Minimum tropospheric pressure (hPa): ', trop_min_p    * 100.0
  WRITE (6,'(A,ES10.4)') &
    ' Maximum tropospheric ABS(PV) (PVU):  ', trop_max_PV   * 1.0E06
  WRITE (6,'(A,ES10.4)') &
    ' Maximum non-trop     pressure (hPa): ', nonTrop_max_p * 100.0
  WRITE (6,*) ''
  WRITE (6,*) 'Layer-by-layer diags:'
  WRITE (6,*) ''
  WRITE (6,*) '        Fraction of  Tropospheric   % change to '// &
                                 '  Non-trop       % change to'
  WRITE (6,*) '        air mass in  vapour mass    tropospheric'// &
                                 '  vapour mass    non-trop'
  WRITE (6,*) ' layer  troposphere  on entry (kg)  vapour mass '// &
                                 '  on entry (kg)  vapour mass'
  WRITE (6,*) ' -----  -----------  -------------  ------------'// &
                                 '  -------------  ------------'

  DO k = qdims%k_end, 1, -1
    trop_frac            = stats(2,k) / stats(1,k)
    trop_vmass_before    = stats(3,k)
    nonTrop_vmass_before = stats(5,k)
    IF (stats(3,k) /= 0.0) THEN
      trop_vmass_percChange    = 100.0 * (stats(4,k) - stats(3,k)) / stats(3,k)
    ELSE
      trop_vmass_percChange    = 0.0
    END IF
    IF (stats(5,k) /= 0.0) THEN
      nonTrop_vmass_percChange = 100.0 * (stats(6,k) - stats(5,k)) / stats(5,k)
    ELSE
      nonTrop_vmass_percChange = 0.0
    END IF
    WRITE (6,'(I6,ES13.3,2(ES15.5,ES14.4))')                       &
      k, trop_frac, trop_vmass_before,    trop_vmass_percChange,   &
                    nonTrop_vmass_before, nonTrop_vmass_percChange
  END DO

  WRITE (6,*) ''
  WRITE (6,*) 'End of QLimits diagnostics'
  WRITE (6,*) '<><><><><><><><><><><><><>'
  WRITE (6,*) ''

END IF ! (L_Diags)
IF (lhook) CALL dr_hook('QLIMITS',zhook_out,zhook_handle)
RETURN


END SUBROUTINE QLimits
END MODULE qlimits_mod
