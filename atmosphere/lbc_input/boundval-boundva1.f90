! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine BOUNDVAL
!
! Purpose : Checks whether a boundary incrementing step and increments
!           boundary values if required.
!
!
! ---------------------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: LBC Input
      SUBROUTINE boundval(                                              &
        lenrim                                                          &
      , l_mcr_qcf2_lbc, l_mcr_qrain_lbc, l_mcr_qgraup_lbc, l_pc2_lbc    &
      , l_murk_lbc, l_int_uvw_lbc                                       &   
      , l_dust_div1_lbc,l_dust_div2_lbc                                 &
      , l_dust_div3_lbc,l_dust_div4_lbc                                 &
      , l_dust_div5_lbc,l_dust_div6_lbc                                 &
      , l_SO2_lbc,l_dms_lbc,l_SO4_aitken_lbc                            &
      , l_SO4_accu_lbc,l_SO4_diss_lbc                                   &
      , l_nh3_lbc,l_soot_new_lbc,l_soot_agd_lbc                         &
      , l_soot_cld_lbc,l_bmass_new_lbc                                  &
      , l_bmass_agd_lbc,l_bmass_cld_lbc                                 &
      , l_ocff_new_lbc,l_ocff_agd_lbc,l_ocff_cld_lbc                    &
      , l_nitr_acc_lbc, l_nitr_diss_lbc                                 &
      , u_lbc, u_lbc_tend                                               &
      , v_lbc, v_lbc_tend                                               &
      , w_lbc, w_lbc_tend                                               &
      , rho_lbc, rho_lbc_tend                                           &
      , theta_lbc, theta_lbc_tend                                       &
      , q_lbc, q_lbc_tend                                               &
      , qcl_lbc, qcl_lbc_tend                                           &
      , qcf_lbc, qcf_lbc_tend                                           &
      , qcf2_lbc, qcf2_lbc_tend                                         &
      , qrain_lbc, qrain_lbc_tend                                       &
      , qgraup_lbc, qgraup_lbc_tend                                     &
      , cf_bulk_lbc, cf_bulk_lbc_tend                                   &
      , cf_liquid_lbc, cf_liquid_lbc_tend                               &
      , cf_frozen_lbc, cf_frozen_lbc_tend                               &
      , exner_lbc, exner_lbc_tend                                       &
      , u_adv_lbc, u_adv_lbc_tend                                       &
      , v_adv_lbc, v_adv_lbc_tend                                       &
      , w_adv_lbc, w_adv_lbc_tend                                       &
      , murk_lbc, murk_lbc_tend                                         &
      , dust_div1_lbc, dust_div1_lbc_tend                               &
      , dust_div2_lbc, dust_div2_lbc_tend                               &
      , dust_div3_lbc, dust_div3_lbc_tend                               &
      , dust_div4_lbc, dust_div4_lbc_tend                               &
      , dust_div5_lbc, dust_div5_lbc_tend                               &
      , dust_div6_lbc, dust_div6_lbc_tend                               &
      , SO2_lbc, SO2_lbc_tend                                           &
      , dms_lbc, dms_lbc_tend                                           &
      , SO4_aitken_lbc, SO4_aitken_lbc_tend                             &
      , SO4_accu_lbc, SO4_accu_lbc_tend                                 &
      , SO4_diss_lbc, SO4_diss_lbc_tend                                 &
      , NH3_lbc, NH3_lbc_tend                                           &
      , soot_new_lbc, soot_new_lbc_tend                                 &
      , soot_agd_lbc, soot_agd_lbc_tend                                 &
      , soot_cld_lbc, soot_cld_lbc_tend                                 &
      , bmass_new_lbc, bmass_new_lbc_tend                               &
      , bmass_agd_lbc, bmass_agd_lbc_tend                               &
      , bmass_cld_lbc, bmass_cld_lbc_tend                               &
      , ocff_new_lbc, ocff_new_lbc_tend                                 &
      , ocff_agd_lbc, ocff_agd_lbc_tend                                 &
      , ocff_cld_lbc, ocff_cld_lbc_tend                                 &
      , nitr_acc_lbc, nitr_acc_lbc_tend                                 &
      , nitr_diss_lbc, nitr_diss_lbc_tend                               &
      , tracer_lbc, tracer_lbc_tend                                     &
      , tracer_ukca_lbc, tracer_ukca_lbc_tend                           &     
      , IO1, IO2                                                        &
      , ICODE, CMESSAGE)

      USE Submodel_Mod
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE lbc_mod
      USE ereport_mod, ONLY : ereport
      USE UM_ParParams
      USE Control_Max_Sizes
      USE atm_fields_bounds_mod, ONLY: udims_s, vdims_s, wdims_s,       &
                                       pdims_s, tdims_s, qdims_l,       &
                                       trdims_s
      USE lbc_read_data_mod, ONLY: current_lbc_step
      USE lookup_addresses

      IMPLICIT NONE

! Parameters required for argument declarations
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
! TYPBND - needs TYPSIZE included first
      ! Headers from atmosphere boundary data sets
      ! Second index of header arrays = 1 Lateral boundary data
      ! = 2 Lower boundary data
      INTEGER :: FIXHD_BOUNDA(LEN_FIXHD,1)      ! Fixed header
      INTEGER :: INTHD_BOUNDA(A_LEN_INTHD,1)    ! Integer header

        ! Lookups for the first set of LBCs in the LBC file
      INTEGER :: LOOKUP_BOUNDA(LEN1_LOOKUP,RIM_LOOKUPSA)

      ! Varying items from the LOOKUP table for the entire LBC file
      INTEGER :: LOOKUP_COMP_BOUNDA(LEN1_LBC_COMP_LOOKUP,BOUND_LOOKUPSA)

      REAL ::  REALHD_BOUNDA(A_LEN_REALHD,1)   ! Real header

      !  Control data calculated from namelist
      INTEGER :: RIM_STEPSA      ! Set by IN_BOUND
      INTEGER :: NBOUND_LOOKUP(2)
      COMMON/CBND/                                                      &
       RIM_STEPSA
! TYPBND end
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

! Arguments:

      INTEGER                                                           &
        lenrim(Nfld_max,NHalo_Max)                                      &
                                  ! IN : Size of a level of LBC
      , IO1, IO2                  ! Offsets to allow for old or new lbcs

      LOGICAL, INTENT (IN) ::                                           &
        l_mcr_qcf2_lbc                                                  &
                         ! true if using second cloud ice lbcs
      , l_mcr_qrain_lbc                                                 &
                         ! true if using rain lbcs
      , l_mcr_qgraup_lbc                                                &
                         ! true if using graupel lbcs
      , l_pc2_lbc                                                       &
                         ! true if cloud fractions in lbcs
      , l_int_uvw_lbc                                                   &
                          ! true if using interpolated advecting winds 
                          !  in lateral boundaries
      , l_murk_lbc                                                      &
                          ! true if murk aerosol in lbcs
      , l_dust_div1_lbc                                                 &
      , l_dust_div2_lbc                                                 &
      , l_dust_div3_lbc                                                 &
      , l_dust_div4_lbc                                                 &
      , l_dust_div5_lbc                                                 &
      , l_dust_div6_lbc                                                 &
                         ! true if DUST in lbcs
      , l_SO2_lbc                                                       &
                         ! true if SO2 in lbcs
      , l_dms_lbc                                                       &
                         ! true if DMS in lbcs
      , l_SO4_aitken_lbc                                                &
                         ! true if SO4_AITKEN in lbcs
      , l_SO4_accu_lbc                                                  &
                         ! true if SO4_accU in lbcs           
      , l_SO4_diss_lbc                                                  &
                         ! true if SO4_DISS in lbcs            
      , l_nh3_lbc                                                       &
                         ! true if NH3 in lbcs       
      , l_soot_new_lbc                                                  &
      , l_soot_agd_lbc                                                  &
      , l_soot_cld_lbc                                                  &
                         ! true if soot in lbcs       
      , l_bmass_new_lbc                                                 &
      , l_bmass_agd_lbc                                                 &
      , l_bmass_cld_lbc                                                 &
                         ! true if biomass in lbcs
      , l_ocff_new_lbc                                                  &
      , l_ocff_agd_lbc                                                  &
      , l_ocff_cld_lbc                                                  &
                         ! true if fossil fuel aerosol in lbcs
      , l_nitr_acc_lbc                                                  &
      , l_nitr_diss_lbc   ! true if fossil fuel aerosol in lbcs

! Note: LBCs are the current value of the LBC (and will be updated)
!       LBC_tends are the value of the LBC at the end of the LBC
!       period, towards which the LBC will tend.

      REAL, INTENT (INOUT) ::                                           &
        u_lbc(lenrim(fld_type_u,halo_type_extended),                    &
              udims_s%k_start:udims_s%k_end)                            &
                                    ! IN/OUT : U LBC
      , u_lbc_tend(lenrim(fld_type_u,halo_type_extended),               &
                          udims_s%k_start:udims_s%k_end)                &
      , v_lbc(lenrim(fld_type_v,halo_type_extended),                    &
              vdims_s%k_start:vdims_s%k_end)                            &
                                    ! IN/OUT : V LBC
      , v_lbc_tend(lenrim(fld_type_v,halo_type_extended),               &
                   vdims_s%k_start:vdims_s%k_end)                       &
                                    ! IN : V LBC tendency
      , w_lbc(lenrim(fld_type_p,halo_type_extended),                    &
              wdims_s%k_start:wdims_s%k_end)                            &
                                    ! IN/OUT : V LBC
      , w_lbc_tend(lenrim(fld_type_p,halo_type_extended),               &
                   wdims_s%k_start:wdims_s%k_end)                       &
                                    ! IN : V LBC tendency
      , rho_lbc(lenrim(fld_type_p,halo_type_extended),                  &
                pdims_s%k_start:pdims_s%k_end)                          &
                                    ! IN/OUT : rho LBC
      , rho_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
                     pdims_s%k_start:pdims_s%k_end)                     &
                                    ! IN : rho LBC tendency
      , theta_lbc(lenrim(fld_type_p,halo_type_extended),                &
                  tdims_s%k_start:tdims_s%k_end)                        &
                                    ! IN/OUT : theta LBC
      , theta_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                       tdims_s%k_start:tdims_s%k_end)                   &
                                    ! IN : theta LBC tendency
      , q_lbc(lenrim(fld_type_p,halo_type_extended),                    &
              qdims_l%k_start:qdims_l%k_end)                            &
                                    ! IN/OUT : Q LBC
      , q_lbc_tend(lenrim(fld_type_p,halo_type_extended),               &
                   qdims_l%k_start:qdims_l%k_end)                       &
                                    ! IN : Q LBC tendency
      , qcl_lbc(lenrim(fld_type_p,halo_type_extended),                  &
                qdims_l%k_start:qdims_l%k_end)                          &
                                    ! IN/OUT : QCL LBC
      , qcl_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
                     qdims_l%k_start:qdims_l%k_end)                     &
                                    ! IN : QCL LBC tendency
      , qcf_lbc(lenrim(fld_type_p,halo_type_extended),                  &
                qdims_l%k_start:qdims_l%k_end)                          &
                                    ! IN/OUT : QCL LBC
      , qcf_lbc_tend(lenrim(fld_type_p,halo_type_extended),             &
                     qdims_l%k_start:qdims_l%k_end)                     &
                                    ! IN : QCL LBC tendency
      , qcf2_lbc(lenrim(fld_type_p,halo_type_extended),                 &
                 qdims_l%k_start:qdims_l%k_end)                         &
                                    ! IN/OUT : QCF2 LBC
      , qcf2_lbc_tend(lenrim(fld_type_p,halo_type_extended),            &
                      qdims_l%k_start:qdims_l%k_end)                    &
                                    ! IN : QCF2 LBC tendency
      , qrain_lbc(lenrim(fld_type_p,halo_type_extended),                &
                  qdims_l%k_start:qdims_l%k_end)                        &
                                    ! IN/OUT : QRAIN LBC
      , qrain_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                       qdims_l%k_start:qdims_l%k_end)                   &
                                    ! IN : QRAIN LBC tendency
      , qgraup_lbc(lenrim(fld_type_p,halo_type_extended),               &
                   qdims_l%k_start:qdims_l%k_end)                       &
                                    ! IN/OUT : QGRAUP LBC
      , qgraup_lbc_tend(lenrim(fld_type_p,halo_type_extended),          &
                        qdims_l%k_start:qdims_l%k_end)                  &
                                    ! IN : QGRAUP LBC tendency
      , cf_bulk_lbc(lenrim(fld_type_p,halo_type_extended),              &
                    qdims_l%k_start:qdims_l%k_end)                      &
                                    ! IN/OUT : CF_BULK LBC
      , cf_bulk_lbc_tend(lenrim(fld_type_p,halo_type_extended),         &
                         qdims_l%k_start:qdims_l%k_end)                 &
                                    ! IN : CF_BULK LBC tendency
      , cf_liquid_lbc(lenrim(fld_type_p,halo_type_extended),            &
                      qdims_l%k_start:qdims_l%k_end)                    &
                                    ! IN/OUT : CF_LIQUID LBC
      , cf_liquid_lbc_tend(lenrim(fld_type_p,halo_type_extended),       &
                                  qdims_l%k_start:qdims_l%k_end)        &
                                    ! IN : CF_LIQUID LBC tendency
      , cf_frozen_lbc(lenrim(fld_type_p,halo_type_extended),            &
                      qdims_l%k_start:qdims_l%k_end)                    &
                                    ! IN/OUT : CF_FROZEN LBC
      , cf_frozen_lbc_tend(lenrim(fld_type_p,halo_type_extended),       &
                           qdims_l%k_start:qdims_l%k_end)               &
                                    ! IN : CF_FROZEN LBC tendency
      , exner_lbc(lenrim(fld_type_p,halo_type_extended),                &
                  pdims_s%k_start:pdims_s%k_end+1)                      &
                                    ! IN/OUT : Exner LBC
      , exner_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                       pdims_s%k_start:pdims_s%k_end+1)                 &
                                         ! IN : Exner LBC tendency
      , u_adv_lbc(lenrim(fld_type_u,halo_type_extended),                &
                  udims_s%k_start:udims_s%k_end)                        &
                                    ! IN/OUT : u_adv LBC
      , u_adv_lbc_tend(lenrim(fld_type_u,halo_type_extended),           &
                       udims_s%k_start:udims_s%k_end)                   &
                                    ! IN : u_adv LBC tendency
      , v_adv_lbc(lenrim(fld_type_v,halo_type_extended),                &
                  vdims_s%k_start:vdims_s%k_end)                        &
                                    ! IN/OUT : v_adv LBC
      , v_adv_lbc_tend(lenrim(fld_type_v,halo_type_extended),           &
                       vdims_s%k_start:vdims_s%k_end)                   &
                                    ! IN : v_adv LBC tendency
      , w_adv_lbc(lenrim(fld_type_p,halo_type_extended),                &
                         wdims_s%k_start:wdims_s%k_end)                 &
                                    ! IN/OUT : W LBC
      , w_adv_lbc_tend(lenrim(fld_type_p,halo_type_extended),           &
                       wdims_s%k_start:wdims_s%k_end)                   &
                                       ! IN : W LBC tendency
      , murk_lbc(lenrim(fld_type_p,halo_type_single),                   &
                 tdims_s%k_start:tdims_s%k_end)                         &
                                    ! IN/OUT : MURK LBC
      , murk_lbc_tend(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end) ! IN : MURK LBC tendency

      REAL, INTENT (INOUT) ::                                           &    
        dust_div1_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : dust_div1 LBC
      , dust_div1_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : dust_div1 LBC tendency
      , dust_div2_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : dust_div2 LBC
      , dust_div2_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : dust_div2 LBC tendency
      , dust_div3_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : dust_div3 LBC
      , dust_div3_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : dust_div3 LBC tendency
      , dust_div4_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : dust_div4 LBC
      , dust_div4_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : dust_div4 LBC tendency
      , dust_div5_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : dust_div5 LBC
      , dust_div5_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : dust_div5 LBC tendency
      , dust_div6_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : dust_div6 LBC
      , dust_div6_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : dust_div6 LBC tendency
      , SO2_lbc(lenrim(fld_type_p,halo_type_single),                    &
                tdims_s%k_start:tdims_s%k_end)                          &
                                    ! IN/OUT : SO2 LBC
      , SO2_lbc_tend(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN : SO2 LBC tendency
      , dms_lbc(lenrim(fld_type_p,halo_type_single),                    &
                tdims_s%k_start:tdims_s%k_end)                          &
                                    ! IN/OUT : DMS LBC
      , dms_lbc_tend(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN : DMS LBC tendency
      , SO4_aitken_lbc(lenrim(fld_type_p,halo_type_single),             &
                       tdims_s%k_start:tdims_s%k_end)                   &
                                    ! IN/OUT : SO4_AITKEN LBC
      , SO4_aitken_lbc_tend(lenrim(fld_type_p,halo_type_single),        &
                            tdims_s%k_start:tdims_s%k_end)              &
                                    ! IN : SO4_AITKEN LBC tendency
      , SO4_accu_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : SO4_accU LBC
      , SO4_accu_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : SO4_accU LBC tendency
      , SO4_diss_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : SO4_DISS LBC
      , SO4_diss_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : SO4_DISS LBC tendency
      , NH3_lbc(lenrim(fld_type_p,halo_type_single),                    &
                tdims_s%k_start:tdims_s%k_end)                          &
                                    ! IN/OUT : NH3 LBC
      , NH3_lbc_tend(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN : NH3 LBC tendency
      , soot_new_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : soot_NEW LBC
      , soot_new_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : soot_NEW LBC tendency
      , soot_agd_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : soot_agd LBC
      , soot_agd_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : soot_agd LBC tendency
      , soot_cld_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : soot_CLD LBC
      , soot_cld_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : soot_CLD LBC tendency
      , bmass_new_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : bmass_NEW LBC
      , bmass_new_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : bmass_NEW LBC tendency
      , bmass_agd_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : bmass_agd LBC
      , bmass_agd_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : bmass_agd LBC tendency
      , bmass_cld_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : bmass_CLD LBC
      , bmass_cld_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : bmass_CLD LBC tendency
      , ocff_new_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : ocff_NEW LBC
      , ocff_new_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : ocff_NEW LBC tendency
      , ocff_agd_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : ocff_agd LBC
      , ocff_agd_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : ocff_agd LBC tendency
      , ocff_cld_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : ocff_CLD LBC
      , ocff_cld_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : ocff_CLD LBC tendency
      , nitr_acc_lbc(lenrim(fld_type_p,halo_type_single),               &
                     tdims_s%k_start:tdims_s%k_end)                     &
                                    ! IN/OUT : nitr_acc LBC
      , nitr_acc_lbc_tend(lenrim(fld_type_p,halo_type_single),          &
                          tdims_s%k_start:tdims_s%k_end)                &
                                    ! IN : nitr_acc LBC tendency
      , nitr_diss_lbc(lenrim(fld_type_p,halo_type_single),              &
                      tdims_s%k_start:tdims_s%k_end)                    &
                                    ! IN/OUT : nitr_DISS LBC
      , nitr_diss_lbc_tend(lenrim(fld_type_p,halo_type_single),         &
                           tdims_s%k_start:tdims_s%k_end)               &
                                    ! IN : nitr_DISS LBC tendency
      , tracer_lbc(lenrim(fld_type_p,halo_type_extended),               &
                          trdims_s%k_start:trdims_s%k_end,tr_lbc_vars)  &
                                           ! IN/OUT : Tracer LBCs
      , tracer_lbc_tend(lenrim(fld_type_p,halo_type_extended),          &
                        trdims_s%k_start:trdims_s%k_end,tr_lbc_vars)    &
                                           ! IN : Tracer LBCs tendency  
      , tracer_ukca_lbc(lenrim(fld_type_p,halo_type_extended),          &
                        trdims_s%k_start:trdims_s%k_end,tr_lbc_ukca)    &
                                           ! IN/OUT : UKCA Tracer LBCs
      , tracer_ukca_lbc_tend(lenrim(fld_type_p,halo_type_extended),     &
                             trdims_s%k_start:trdims_s%k_end,tr_lbc_ukca)
                                           ! IN : UKCA Tracer LBCs tend


      INTEGER ::                                                        &
        icode                           ! Error code

      CHARACTER(LEN=80) ::                                              &
        cmessage                        ! Error message

! Local variables

      INTEGER       :: steps_to_next_update

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! --------------------------------------------------------------------

      IF (lhook) CALL dr_hook('BOUNDVAL',zhook_in,zhook_handle)

      CMESSAGE=' '

! Check that the lbcs haven't lost synchronisation and then call
! increment_lbcs to increment them

      IF (current_lbc_step  /=  stepim(atmos_im) + IO1) THEN

        IF ((current_lbc_STEP /= STEPim(atmos_im) - IO2) .AND.          &
               (MOD(BNDARY_OFFSETim(atmos_im)+STEPim(atmos_im) + IO1,   &
                    RIM_STEPSA) /= 2) .AND.                             &
               (MOD(BNDARY_OFFSETim(atmos_im)+STEPim(atmos_im) + IO1,   &
                    RIM_STEPSA) /= 0)) THEN

          ICODE=10

          Call Ereport("BOUNDVA1 ", ICODE,                              &
                       "LBC values have lost synchronisation" )
        END IF

        steps_to_next_update = RIM_STEPSA -                             &
                     MOD(BNDARY_OFFSETim(atmos_im)+STEPim(atmos_im)-1,  &
                         RIM_STEPSA) + IO2

! DEPENDS ON: increment_atmos_lbcs
        CALL increment_atmos_lbcs(                                      &
              steps_to_next_update,                                     &
              lenrim,                                                   &
              model_levels,wet_levels,tr_lbc_vars,tr_levels,            &
              tr_lbc_ukca,                                              &
              l_mcr_qcf2_lbc, l_mcr_qrain_lbc, l_mcr_qgraup_lbc,        &
              l_pc2_lbc, l_murk_lbc, l_int_uvw_lbc,                     &
              l_dust_div1_lbc,l_dust_div2_lbc,                          &
              l_dust_div3_lbc,l_dust_div4_lbc,                          &
              l_dust_div5_lbc,l_dust_div6_lbc,                          &
              l_SO2_lbc,l_dms_lbc,l_SO4_aitken_lbc,                     &
              l_SO4_accu_lbc,l_SO4_diss_lbc,                            &
              l_nh3_lbc,l_soot_new_lbc,l_soot_agd_lbc,                  &
              l_soot_cld_lbc,l_bmass_new_lbc,                           &
              l_bmass_agd_lbc,l_bmass_cld_lbc,                          &
              l_ocff_new_lbc,l_ocff_agd_lbc,l_ocff_cld_lbc,             &
              l_nitr_acc_lbc,l_nitr_diss_lbc,                           &
              u_lbc,u_lbc_tend,                                         &
              v_lbc,v_lbc_tend,                                         &
              w_lbc,w_lbc_tend,                                         &
              rho_lbc,rho_lbc_tend,                                     &
              theta_lbc,theta_lbc_tend,                                 &
              q_lbc,q_lbc_tend,                                         &
              qcl_lbc,qcl_lbc_tend,                                     &
              qcf_lbc,qcf_lbc_tend,                                     &
              qcf2_lbc,qcf2_lbc_tend,                                   &
              qrain_lbc,qrain_lbc_tend,                                 &
              qgraup_lbc,qgraup_lbc_tend,                               &
              cf_bulk_lbc,cf_bulk_lbc_tend,                             &
              cf_liquid_lbc,cf_liquid_lbc_tend,                         &
              cf_frozen_lbc,cf_frozen_lbc_tend,                         &
              exner_lbc,exner_lbc_tend,                                 &
              u_adv_lbc,u_adv_lbc_tend,                                 &
              v_adv_lbc,v_adv_lbc_tend,                                 &
              w_adv_lbc,w_adv_lbc_tend,                                 &
              murk_lbc,murk_lbc_tend,                                   &
              dust_div1_lbc,dust_div1_lbc_tend,                         &
              dust_div2_lbc,dust_div2_lbc_tend,                         &
              dust_div3_lbc,dust_div3_lbc_tend,                         &
              dust_div4_lbc,dust_div4_lbc_tend,                         &
              dust_div5_lbc,dust_div5_lbc_tend,                         &
              dust_div6_lbc,dust_div6_lbc_tend,                         &
              SO2_lbc,SO2_lbc_tend,                                     &
              dms_lbc,dms_lbc_tend,                                     &
              SO4_aitken_lbc,SO4_aitken_lbc_tend,                       &
              SO4_accu_lbc,SO4_accu_lbc_tend,                           &
              SO4_diss_lbc,SO4_diss_lbc_tend,                           &
              NH3_lbc,NH3_lbc_tend,                                     &
              soot_new_lbc,soot_new_lbc_tend,                           &
              soot_agd_lbc,soot_agd_lbc_tend,                           &
              soot_cld_lbc,soot_cld_lbc_tend,                           &
              bmass_new_lbc,bmass_new_lbc_tend,                         &
              bmass_agd_lbc,bmass_agd_lbc_tend,                         &
              bmass_cld_lbc,bmass_cld_lbc_tend,                         &
              ocff_new_lbc,ocff_new_lbc_tend,                           &
              ocff_agd_lbc,ocff_agd_lbc_tend,                           &
              ocff_cld_lbc,ocff_cld_lbc_tend,                           &
              nitr_acc_lbc, nitr_acc_lbc_tend,                          &
              nitr_diss_lbc, nitr_diss_lbc_tend,                        &
              tracer_lbc,tracer_lbc_tend,                               &
              tracer_ukca_lbc,tracer_ukca_lbc_tend,                     &
              RIM_STEPSA,ICODE)

!       IF (ICODE  >   0) THEN
!             WRITE(6,*) 'Failure in INCREMENT_ATMOS_lbcS while ',      &
!    &                   'attempting to update LBCS for timestep ',     &
!    &                   STEPim(a_im)
!         GOTO 9999
!       ENDIF

        current_lbc_STEP = STEPim(atmos_im) + IO1

      END IF ! IF (current_lbc_STEP  /=  STEPim(atmos_im))

 9999 CONTINUE
      IF (lhook) CALL dr_hook('BOUNDVAL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE boundval
