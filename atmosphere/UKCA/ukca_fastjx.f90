! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Main routine for calculating online photolysis rates
!   using fastj-x. Borrows a lot from the fast-j main routine. 
!   Further developments that are still required include:
!    (i) use of MODE aerosol optical depths
!    (ii) use of more complete blocking for potential speedup
!    (iii) update to fast-jx 6.6 which will have some small
!    improvements (eg will remove the requirement that jaceto is last)
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      SUBROUTINE UKCA_FASTJX(                                          &
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
        dj,                                                            &
        p_layer_boundaries, pstar,                                     &
        p_theta_levels,                                                &
        t_theta_levels,                                                &
        rho_r2,                                                        &
        z_top_of_model,                                                &
        so4_aitken, so4_accum,                                         &
        qcl, qcf, rel_humid_frac, area_cloud_fraction,                 &    
        conv_cloud_lwp, conv_cloud_top, conv_cloud_base,               &
        conv_cloud_amount,                                             &
        surf_albedo,                                                   &
        tropopause_level,                                              &
        ozone,                                                         &
        land_fraction)

      USE FASTJX_DATA,    ONLY: fastjx_set_limits,                     &
                                fastjx_allocate_memory,                &
                                fastjx_deallocate_memory, cmessage,    &
                                jjpnl,nsl,                             &
                                Blocking,kpcx,jpcl,                    &
                                aer_swband,                            &
                                tau, daynumber,                        &
                                sza, szafac, sza_2d, szafac_2d, u0,    &
                                sa_block,                              &
                                rz_3d, rz_all, zzht,                   &
                                pz_3d, pz_all, tz_3d, sa_2d,           &
                                dm_3d, o3_3d,                          &
                                ods_3d, odw_3d, odi_3d
      USE ukca_option_mod,   ONLY: jppj, l_ukca
      USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels,       &
                                   eta_theta_levels
      USE trignometric_mod,  ONLY: true_longitude, sin_theta_latitude
      USE dyn_coriolis_mod,  ONLY: f3_at_u
      USE spec_sw_lw,        ONLY: sw_spectrum, spectrum
      USE rad_pcf,           ONLY: ip_aerosol_param_moist,             &
                                   ip_accum_sulphate,                  &
                                   ip_aitken_sulphate
      USE Control_Max_Sizes
      USE UM_ParVars
      USE UKCA_CONSTANTS,       ONLY: c_o3, m_s, m_nh42so4
      USE earth_constants_mod,  ONLY: g, two_omega
      USE water_constants_mod,  ONLY: tm
      USE conversions_mod,      ONLY: pi_over_180
      USE PrintStatus_mod
      USE ereport_mod,          ONLY: ereport
      USE yomhook,              ONLY: lhook, dr_hook
      USE parkind1,             ONLY: jprb, jpim
      USE um_input_control_mod, ONLY: lcal360
      USE cloud_inputs_mod,     ONLY: l_pc2
      USE cv_run_mod,           ONLY: l_3d_cca
      USE Submodel_Mod

      IMPLICIT NONE

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
! TYPLNDM
! Formerly integral part of TYPCONA, the variables below have been
! separated from the rest of TYPCONA as they are required by some
! of the Ocean routines in the Ocean-Atmosphere configuration of
! the UM whilest TYPCONA is not.

      ! Primary Arrays
      INTEGER::land_points     ! No. of land points  (can be 0)
      INTEGER::land_ice_points ! Number of land ice points
      INTEGER::soil_points     ! Number of soil points

      ! Do not allow these arrays to have zero size
      INTEGER::land_index    (max(1,land_field)) ! set from land_sea_mask
      INTEGER::land_ice_index(max(1,land_field)) ! Array of land ice points.
      INTEGER::soil_index    (max(1,land_field)) ! Array of soil points.

      ! Gets some sizes transported around the model :
      COMMON /land_soil_dimensions/                                     &
     &  land_points , land_ice_points , soil_points

! TYPLNDM end
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
 
! Output photolysis rates
      REAL, INTENT(INOUT) :: dj(row_length,rows,model_levels,jppj) 
! Pressure on rho levels
      REAL, INTENT(IN)    :: p_layer_boundaries(row_length,rows,        &
                                                model_levels+1) 
! Pressure on theta levels
      REAL, INTENT(IN)    :: p_theta_levels(row_length,rows,model_levels)
! Temperature
      REAL, INTENT(IN)    :: t_theta_levels(row_length,rows,model_levels)
! Sulphate aerosol (aitken mode) 
      REAL, INTENT(IN)    :: so4_aitken(row_length,rows,model_levels)
! Sulphate aerosol (accumulation mode)
      REAL, INTENT(IN)    :: so4_accum(row_length,rows,model_levels)
! Density * R * R
      REAL, INTENT(IN)    :: rho_r2(row_length,rows,model_levels)
! Top of model
      REAL, INTENT(IN)    :: z_top_of_model
! liquid water cloud
      REAL, INTENT(IN)    :: qcl(row_length,rows,wet_levels)
! ice water cloud
      REAL, INTENT(IN)    :: qcf(row_length,rows,wet_levels)
! Relative humidity fraction
      REAL, INTENT(IN)    :: rel_humid_frac(row_length,rows,wet_levels)
! cloud area fraction
      REAL, INTENT(IN)    :: area_cloud_fraction(row_length,rows,wet_levels)
! convective cloud amount
      REAL, INTENT(IN)    :: conv_cloud_amount(row_length,rows,wet_levels)
! ozone mmr
      REAL, INTENT(IN)    :: ozone(row_length,rows,ozone_levels)
! surface pressure
      REAL, INTENT(IN)    :: pstar(row_length,rows)
! Convective cloud LWP
      REAL, INTENT(IN)    :: conv_cloud_lwp(row_length,rows)
! Surface albedo
      REAL, INTENT(IN)    :: surf_albedo(row_length,rows)
! Tropopause level
      INTEGER, INTENT(IN) :: tropopause_level(row_length,rows)
! Land fraction
      REAL, INTENT(IN)    :: land_fraction(row_length,rows)
! Convective cloud top
      INTEGER, INTENT(IN) :: conv_cloud_top(row_length,rows)
! Convective cloud bottom
      INTEGER, INTENT(IN) :: conv_cloud_base(row_length,rows)


! Local variables

! Sulphate climatology surface area density        
      REAL :: so4_sa(row_length,rows,model_levels)
! Sulphate climatology effective radius
      REAL :: so4_reff(row_length,rows,model_levels)
! Sulphate climatology mass mixing ratio
      REAL :: so4_sa_mmr(row_length,rows,model_levels)
   
! Parameter to add stratospheric aerosol climatology to OD calculation
! Still not sure of the exact implementation so don't add it to cruntime yet
      LOGICAL, PARAMETER  :: L_use_stratclim = .FALSE.

      INTEGER :: icode

! Loop variables
      INTEGER                       :: i,j,k,l,ia, ix
! Effective diameter
      REAL                          :: d_eff
! Sine of true latitude
      REAL                          :: sin_true_latitude(row_length, rows)
! Timestep in hours
      REAL                          :: timej
! Latitude and Longitude of current point
      INTEGER                       :: nslat 
      INTEGER                       :: nslon 
! First time?
      LOGICAL                       :: first=.true.
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!----------------------------------------------
! Interface

      INTERFACE 

        SUBROUTINE FASTJX_INPHOT
          USE     FASTJX_DATA,         ONLY: w_
          USE     FASTJX_SPECS
          IMPLICIT NONE
        END SUBROUTINE FASTJX_INPHOT

        SUBROUTINE FASTJX_SOLAR2(timej,sinlat,longitude,lcal360)
          USE   FASTJX_DATA,    ONLY: daynumber,tau,szamax, &
                                      sza_2d,SZAFAC_2d
          IMPLICIT NONE
          LOGICAL,INTENT(IN)             :: lcal360      
          REAL,INTENT(IN)                :: timej
          REAL,DIMENSION(:,:),INTENT(IN) :: sinlat,longitude 
        END SUBROUTINE FASTJX_SOLAR2

        SUBROUTINE FASTJX_PHOTOJ(zpj)
          USE  FASTJX_DATA
          IMPLICIT NONE
          REAL,DIMENSION(:,:,:,:),INTENT(INOUT) :: zpj
        END SUBROUTINE FASTJX_PHOTOJ
      END INTERFACE

!**************************************
! End of Header

      IF (lhook) CALL dr_hook('UKCA_FASTJX',zhook_in,zhook_handle)

      IF (PrintStatus >= PrStatus_Diag) THEN
        write(6,*) 'UKCA_FASTJX inputs:'
        write(6,*) '1: ',minval(p_layer_boundaries),                    &
                         maxval(p_layer_boundaries),                    &
                         sum(p_layer_boundaries)/                       &
                         size(p_layer_boundaries)
        write(6,*) '2: ',minval(p_theta_levels),maxval(p_theta_levels), &
                         sum(p_theta_levels)/size(p_theta_levels)
        write(6,*) '3: ',minval(t_theta_levels),maxval(t_theta_levels), &
                         sum(t_theta_levels)/size(t_theta_levels)
        write(6,*) '4: ',minval(so4_aitken),maxval(so4_aitken),         &
                         sum(so4_aitken)/size(so4_aitken)
        write(6,*) '5: ',minval(so4_accum),maxval(so4_accum),           &
                         sum(so4_accum)/size(so4_accum)
        write(6,*) '6: ',minval(rho_r2),maxval(rho_r2),                 &
                         sum(rho_r2)/size(rho_r2)
        write(6,*) '7: ',z_top_of_model
        write(6,*) '8: ',minval(qcl),maxval(qcl),sum(qcl)/size(qcl)
        write(6,*) '9: ',minval(qcf),maxval(qcf),sum(qcf)/size(qcf)
        write(6,*) '10:',minval(rel_humid_frac),maxval(rel_humid_frac), &
                         sum(rel_humid_frac)/size(rel_humid_frac)
        write(6,*) '11:',minval(area_cloud_fraction),                   &
                         maxval(area_cloud_fraction),                   &
                         sum(area_cloud_fraction)/                      &
                         size(area_cloud_fraction)
        write(6,*) '12:',minval(conv_cloud_amount),                     &
                         maxval(conv_cloud_amount),                     &
                         sum(conv_cloud_amount)/size(conv_cloud_amount)
        write(6,*) '13:',minval(ozone),maxval(ozone),                   &
                         sum(ozone)/size(ozone)
        write(6,*) '14:',minval(pstar),maxval(pstar),                   &
                         sum(pstar)/size(pstar)
        write(6,*) '15:',minval(conv_cloud_lwp),maxval(conv_cloud_lwp), &
                         sum(conv_cloud_lwp)/size(conv_cloud_lwp)
        write(6,*) '16:',minval(surf_albedo),maxval(surf_albedo),       &
                         sum(surf_albedo)/size(surf_albedo)
        write(6,*) '17:',minval(tropopause_level),                      &
                         maxval(tropopause_level),                      &
                         sum(tropopause_level)/size(tropopause_level)
        write(6,*) '18:',minval(land_fraction),maxval(land_fraction),   &
                         sum(land_fraction)/size(land_fraction)
        write(6,*) '19:',minval(conv_cloud_top),maxval(conv_cloud_top), &
                         sum(conv_cloud_top)/size(conv_cloud_top)
        write(6,*) '20:',minval(conv_cloud_base),                       &
                         maxval(conv_cloud_base),                       &
                         sum(conv_cloud_base)/size(conv_cloud_base)
        icode=0
        call um_fort_flush(6,icode)
      END IF

! Initialise 
      dj(:,:,:,:) = 0. 

! Set Blocking mode:          0) Column-by-column
        !                     1) blocking 1 row 
        !                     2) blocking domain
        !                     3) compressed  (not implemented)
        !                     4) load balancing (not implemented)
        Blocking%Mode = 2

! Allocate arrays etc.
      CALL FASTJX_SET_LIMITS(row_length,rows,model_levels)
      CALL FASTJX_ALLOCATE_MEMORY 

! Read in data from files
      IF (first) THEN
! DEPENDS ON: FASTJX_INPHOT
        CALL FASTJX_INPHOT
        first=.false.
      END IF

! Initialise arrays for levels/units appropriate for fast-j
! Need to update to include aerosols
      CALL FASTJX_SET_ARRAYS

! Set variables concerning model time
! Convert timestep into hours
      timej              = secs_per_stepim(a_im)/3600.

        ! Day of the year
      daynumber          = i_day_number

        ! Time in hours
      tau                = i_hour*1.+i_minute/60.+i_second/3600.      &
                           - timej*0.5

! Calculate solar zenith angles

! DEPENDS ON: fastjx_solar2
      CALL FASTJX_SOLAR2 (timej, sin_true_latitude, true_longitude,   &
                            lcal360)

! Block the data appropriately and call photolysis routines
! Still need to add modes 3 (compressed) and 4 (load balancing)

! Initialise longitude counter to 1
      nslon = 1

      SELECT CASE (Blocking%Mode)

        ! if blocking point by point 
        CASE(0)

          ! Loop over rows
          DO j=1,rows

            ! Loop over row, setting longitude counter
            DO i=1,row_length

              ! Set latitude counter to row number
              nslat  = j
              nslon  = i

              ! Loop over row, setting longitude counter
              nsl(1,1)=nslon
              nsl(2,1)=nslat

              ! block using consistent approach
              DO ix=1, kpcx
                sza(ix) = sza_2d(nsl(1,ix),nsl(2,ix))
                szafac(ix) = szafac_2d(nsl(1,ix),nsl(2,ix))

                u0(ix) = sza(ix)*pi_over_180
                u0(ix) = COS(u0(ix))
              END DO

! DEPENDS ON: fastjx_photoj
              CALL FASTJX_PHOTOJ (dj)

            ENDDO
          ENDDO 

        ! if blocking row by row 
        CASE(1)

          ! Loop over rows
          DO j=1,rows

            ! Set latitude counter to row number
            nslat  = j

            ! Loop over row, setting longitude counter
            DO i=1,row_length
              nsl(1,i)=nslon+(i-1)
              nsl(2,i)=nslat
            ENDDO

            ! block using consistent approach
            DO ix=1, kpcx
              sza(ix) = sza_2d(nsl(1,ix),nsl(2,ix))
              szafac(ix) = szafac_2d(nsl(1,ix),nsl(2,ix))

              u0(ix) = sza(ix)*pi_over_180
              u0(ix) = COS(u0(ix))

              sa_block(ix) = sa_2d(nsl(1,ix),nsl(2,ix))
              sa_block(ix) = MIN(1.0, sa_block(ix))
              sa_block(ix) = MAX(0.0, sa_block(ix))

            END DO

! DEPENDS ON: fastjx_photoj
            CALL FASTJX_PHOTOJ (dj)

          ENDDO  

        !**********************************
        ! If blocking whole domain
        CASE(2)

          ! initialise latitude counter to 1
          nslat = 1

          ! Loop over rows
          DO j=1,rows

            ! Loop over columns
            DO i=1,row_length
 
              ! Calculate positions of longitude and latitude in blocked arrays      
              l=i+(j-1)*row_length

              nsl(1,l)=nslon+(i-1)
              nsl(2,l)=nslat+(j-1)
            END DO ! columns
          END DO ! rows

          ! block using consistent approach
          DO ix=1, kpcx
            sza(ix) = sza_2d(nsl(1,ix),nsl(2,ix))
            szafac(ix) = szafac_2d(nsl(1,ix),nsl(2,ix))

            u0(ix) = sza(ix)*pi_over_180
            u0(ix) = COS(u0(ix))

            sa_block(ix) = sa_2d(nsl(1,ix),nsl(2,ix))
            sa_block(ix) = MIN(1.0, sa_block(ix))
            sa_block(ix) = MAX(0.0, sa_block(ix))

          END DO

! DEPENDS ON: fastjx_photoj
          CALL FASTJX_PHOTOJ (dj)


        !**********************************
        ! Else exit with an error (need to implement blocking types 3 and 4)
        CASE DEFAULT

          cmessage='Blocking Mode does not Exist'
          CALL EREPORT('UKCA_FASTJX', Blocking%Mode,cmessage)          

      END SELECT

        ! Tidy up at the end
      CALL FASTJX_DEALLOCATE_MEMORY

      IF (lhook) CALL dr_hook('UKCA_FASTJX',zhook_out,zhook_handle)
      RETURN

      CONTAINS

! ######################################################################
      SUBROUTINE FASTJX_SET_ARRAYS

      USE UKCA_CONSTANTS,    ONLY: avogadro, m_air
      USE ukca_option_mod,   ONLY: L_ukca_use_background_aerosol

      USE chsunits_mod, ONLY : nunits

      IMPLICIT NONE

        ! Loop variables
      INTEGER                               :: i,j,k,l, ia

        ! Total mass in column
      REAL,DIMENSION(1:row_length,1:rows)                :: total_mass 

        ! Mass column per layer
      REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: d_mass

        ! Relative humidity
      REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: rh
      REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: qsat

        ! Sulphate in accumulation and aitken modes
      REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: sulph_accu
      REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: sulph_aitk
      REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: sulphur

        ! Cloud optical depths
      REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: odi
      REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: odw
      REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: ods

        ! Conversion factor from kg to molecules
      REAL                                               :: masfac

        ! Temp variables used to calc aerosol od
      INTEGER                                            :: i_humidity
      REAL                                               :: delta_humidity
      REAL                                               :: f

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


        !**************************************************
        ! EOH
      IF (lhook) CALL dr_hook('UKCA_FASTJX.FASTJX_SET_ARRAYS',zhook_in, &
                              zhook_handle)

        ! Calculate latitude
      sin_true_latitude = f3_at_u(1:row_length,1:rows) / two_Omega

        ! rz in cm
      rz_3d(:,:,1)              =                                     &
          r_Theta_levels(1:row_length,1:rows,0)*100.
      rz_3d(:,:,2:model_levels) =                                     &
          r_rho_levels  (1:row_length,1:rows,2:model_levels)*100.
      rz_3d(:,:,model_levels+1) =                                     &
          r_Theta_levels(1:row_length,1:rows,model_levels)*100.

        ! calculate pressure at box edges
      pz_3d(:,:,1)              = pstar
      pz_3d(:,:,2:model_levels) = p_layer_boundaries(:,:,2:model_levels)
      pz_3d(:,:,model_levels+1) = 0.

        ! Calculate mass in box from pressure differences 
        ! using hydrostatic approximation
      d_mass = (pz_3d(:,:,1:model_levels)                             &
             - pz_3d(:,:,2:model_levels+1))/g

        ! Calculate total mass within convective clouds
      total_mass = 0.0
      DO k = 1,model_levels
         WHERE ( k <= conv_cloud_top .AND. k >= conv_cloud_base) 
           total_mass = total_mass + d_mass(:,:,k)
         ENDWHERE
      END DO


        ! Initialise water vapour to large scale precipitation
      odw(:,:,1:wet_levels)                = qcl(:,:,1:wet_levels)
      odi(:,:,1:wet_levels)                = qcf(:,:,1:wet_levels)
      odw(:,:,(wet_levels+1):model_levels) = 0.0
      odi(:,:,(wet_levels+1):model_levels) = 0.0

!-----------------------------------------------------------------
! Only add the convective cloud if we are not using the PC2 cloud scheme
      IF (.NOT. L_PC2) THEN

        ! I L-3d_cca (cloud levels?)
      IF (l_3d_cca) THEN      
        DO k = 1,n_cca_lev

          WHERE (conv_cloud_top > 0 .AND. t_theta_levels(:,:,k) > tm)
! If above freezing point evaluate liquid water 
            odw(:,:,k) = odw(:,:,k) + (conv_cloud_lwp(:,:)*           &
                         conv_cloud_amount(:,:,k))/total_mass(:,:)
          ELSEWHERE (conv_cloud_top > 0)
! Else evaluate ice water
            odi(:,:,k) = odi(:,:,k) + (conv_cloud_lwp(:,:)*           &
                         conv_cloud_amount(:,:,k))/total_mass(:,:)
          ENDWHERE
 
        END DO !  k levels

      ELSE

          ! Else loop over model levels
        DO k=1,model_levels

          WHERE (conv_cloud_top >= k .AND. conv_cloud_base <= k .AND.   &
                 t_theta_levels(:,:,k) > tm)
! If above freezing point evaluate liquid water 
            odw(:,:,k) = odw(:,:,k)+(conv_cloud_lwp(:,:)*               &
                         conv_cloud_amount(:,:,1))/total_mass(:,:)
          ELSEWHERE (conv_cloud_top >= k .AND. conv_cloud_base <= k)
! Else evaluate ice water
            odi(:,:,k) = odi(:,:,k)+(conv_cloud_lwp(:,:)*               &
                         conv_cloud_amount(:,:,1))/total_mass(i,j)
          ENDWHERE

        END DO ! k levels
      END IF ! l_3d_cca

      END IF ! L_PC2

! Convert mass mixing ratios to column densities.
      odw = odw*d_mass
      odi = odi*d_mass

!--------------------------------------------------------------
! set effective radii for water drops, different for land&sea
! set effective diameter for ice crystals
! Formulae from John Edwards to convert from column densities to
! optical depths

      d_eff=100.0 ! in microns

      DO k=1,model_levels

        WHERE (land_fraction > 0.5) 
          odw(:,:,k) = odw(:,:,k)*(-8.86964+1.67373E3/6.)
        ELSEWHERE
          odw(:,:,k) = odw(:,:,k)*(-8.86964+1.67373E3/12.)
        ENDWHERE

      END DO

      odi = odi*(-2.189E-3+3.311E3/d_eff+3.611/d_eff**2)

! Use approach of Briegleb(1992), JGR 97, 7603.
! to account for random overlap of cloud layers
      odw(:,:,1:wet_levels)    =  odw(:,:,1:wet_levels) &
                   * (area_cloud_fraction(:,:,1:wet_levels))**(1.5)
      odi(:,:,1:wet_levels)    =  odi(:,:,1:wet_levels) &
                   * (area_cloud_fraction(:,:,1:wet_levels))**(1.5)

! Set optical depth in top layer to be the same as top+1 
! could just set to 0. For top level shouldn't make any difference
      odw_3d(:,:,1:model_levels)      = odw
      odw_3d(:,:,(model_levels+1))    = odw_3d(:,:,(model_levels))
      odi_3d(:,:,1:model_levels)      = odi
      odi_3d(:,:,(model_levels+1))    = odi_3d(:,:,(model_levels))

!***********************************************************
! Calculate aerosol od columns @600nm
! Note this uses the classic aerosol scheme. Work needs doing to
! incorporate the mode aerosol scheme


! Calculate optical depth from surface aerosol density

! Consistent approach with stratospheric heterogeneous chemistry
      IF( L_use_stratclim) THEN

! DEPENDS ON: ukca_read_aerosol
        CALL UKCA_READ_AEROSOL(i_year, i_month, row_length, rows,      &
            model_levels,                                              &
            sin_theta_latitude(1:row_length, 1:rows),                  &
            eta_theta_levels(1:model_levels) * z_top_of_model,         &
            L_ukca_use_background_aerosol,so4_sa)
          
! DEPENDS ON: ukca_read_reff
        CALL UKCA_READ_REFF(i_year, i_month, row_length, rows,         &
            model_levels,                                              &
            sin_theta_latitude(1:row_length, 1:rows),                  &
            eta_theta_levels(1:model_levels) * z_top_of_model,         &
            so4_reff)
          
! Convert surface aerosol density to mmr
! Note this routine uses assumptions that cause OD to be
! underestimated during Pinatubo period (reff too small)
! DEPENDS ON: ukca_strat_aero_clim_mmr
        CALL UKCA_STRAT_AERO_CLIM_MMR(row_length, rows, model_levels,  &
            r_rho_levels(1:row_length,1:rows,1:model_levels),          &
            rho_r2(1:row_length,1:rows,1:model_levels),                & 
            so4_sa(1:row_length,1:rows,1:model_levels),                &
            so4_reff(1:row_length,1:rows,1:model_levels),              &
            so4_sa_mmr(1:row_length,1:rows,1:model_levels))

          ! Add mmr from this process to accumulation mode
        DO i=1, rows
          DO j=1, row_length
            sulph_accu(j,i,(tropopause_level(j,i)):model_levels) =     &
            so4_sa_mmr(j,i,(tropopause_level(j,i)):model_levels)
          END DO
        END DO
      ELSE
! Use input aerosol concentrations
        sulph_accu(:,:,:) = so4_accum(:,:,:)
        sulph_aitk(:,:,:) = so4_aitken(:,:,:)
      END IF

!-----------------------------------------------------------
! Calculate relative humidity

      rh = rel_humid_frac !             *1.0E2
      IF (model_levels > wet_levels)                                    &
              rh(:,:,wet_levels+1:model_levels) = 0.01

! Adopt ukca_fastj approach (see ukca_fast.F90)
! Multiply by molecular weight ratio to convert from mass mixing ratio of 
! sulphur atoms to mass mixing ratio of ammonium sulphate.
! Use d_mass to convert to mass of ammonium sulphate.   
      sulph_accu = sulph_accu*d_mass*m_nh42so4/m_s
      sulph_aitk = sulph_aitk*d_mass*m_nh42so4/m_s

! Loop over aerosol types
      DO ia=1,sw_spectrum(1)%n_aerosol

          ! Select sulfate in aitken and accumulation modes
        IF (sw_spectrum(1)%type_aerosol(ia) == ip_accum_sulphate .OR.   &
            sw_spectrum(1)%type_aerosol(ia) == ip_aitken_sulphate) THEN

          IF (sw_spectrum(1)%i_aerosol_parametrization(ia) ==           &
                 ip_aerosol_param_moist) THEN ! moist aerosol

              ! evaluate size of humidity bins
            delta_humidity = sw_spectrum(1)%humidities(2,ia)-           &
                              sw_spectrum(1)%humidities(1,ia)

              ! Loop over grid cells calculating optical depths
            DO k=1,model_levels
              DO j = 1, rows
                DO i = 1, row_length

                   ! Determine humidity bin of parameterisation 
                  i_humidity = INT(rh(i,j,k)/delta_humidity)+1
 
                  ! fraction of bin
                  f = (rh(i,j,k)-                                       &
                       sw_spectrum(1)%humidities(i_humidity,ia))/       &
                       delta_humidity

                   ! Determine optical depths from parameterisation
                  ods(i,j,k) =                                          &
                  (sw_spectrum(1)%aerosol_absorption(i_humidity,ia,     &
                                                     aer_swband)+       &
                   sw_spectrum(1)%aerosol_scattering(i_humidity,ia,     &
                                                     aer_swband))*f+    &
                  (sw_spectrum(1)%aerosol_absorption(i_humidity+1,ia,   &
                                                     aer_swband)+       &
                   sw_spectrum(1)%aerosol_scattering(i_humidity+1,ia,   &
                                                     aer_swband))*(1-f)
                END DO
              END DO
            END DO

          ELSE ! dry aerosol
            DO k=1,model_levels
              DO j = 1, rows
                DO i = 1, row_length
                  ods(i,j,k) =                                          &
                    sw_spectrum(1)%aerosol_absorption(1,ia,aer_swband)+ &
                    sw_spectrum(1)%aerosol_scattering(1,ia,aer_swband)
                END DO
              END DO
            END DO
          END IF

          IF (sw_spectrum(1)%type_aerosol(ia) == ip_accum_sulphate)  &
              sulph_accu = sulph_accu*ods
          IF (sw_spectrum(1)%type_aerosol(ia) == ip_aitken_sulphate) &
              sulph_aitk = sulph_aitk*ods
        END IF   ! type_aerosol_sw
      END DO  ! ia

! Sum the aitkin and accumulation type optical depths
      sulphur = sulph_aitk + sulph_accu

! Copy sulphate od to global array
      ods_3d(:,:,1:model_levels)  = sulphur(:,:,1:model_levels) 
      ods_3d(:,:,(model_levels+1))= sulphur(:,:,model_levels)

      IF (printstatus >= prstatus_diag) THEN
        DO i=1, model_levels
          WRITE(6,*) 'fastjX converter max ', i,  maxval(so4_sa_mmr(:,:,i))
        END DO
      END IF

!*********************************************************
! Set other variables here

! Set surface albedo
      sa_2d                           = surf_albedo

! Calculate pressure at box edges (include box to TOA)
! NB convert Pa to hPa for fastj routines
      pz_all(:,:,1:(model_levels+1))  = p_layer_boundaries(:,:,:)/100.0
      pz_all(:,:,(model_levels+2))    = 0.0E0

! Calculate heights of box edges
      rz_all(:,:,1:model_levels+1)    = rz_3d(:,:,1:model_levels+1) 
      rz_all(:,:,model_levels+2)      = rz_all(:,:,model_levels+1)      &
                                        + zzht

! Set temperature to be temperature
! Use top temperature for top+1 level
      tz_3d(:,:,1:model_levels)       = t_theta_levels(:,:,:)
      tz_3d(:,:,(model_levels+1))     = t_theta_levels(:,:,model_levels)

! Copy air mass from local array
      dm_3d(:,:,1:model_levels)       = d_mass(:,:,1:model_levels) 
! Evaluate top box from difference with TOA p i.e. 0!
      dm_3d(:,:,(model_levels+1))     = pz_3d(:,:,model_levels+1)/g 

! Convert air mass from kg/m2 to molecules/cm2
      masfac                          = avogadro/(m_air*1.0E4)    
      dm_3d(:,:,:)                    = dm_3d(:,:,:)*masfac

! Convert ozone mmr to vmr and write to fj_ozone array
      o3_3d(:,:,1:model_levels)       = ozone/c_o3
      o3_3d(:,:,(model_levels+1))     = o3_3d(:,:,model_levels)

! Multiply by molecules per m2 to get it in units fastj wants
      o3_3d                           = o3_3d*dm_3d

      IF (lhook) CALL dr_hook('UKCA_FASTJX.FASTJX_SET_ARRAYS',zhook_out, &
                               zhook_handle)
      RETURN
      END SUBROUTINE FASTJX_SET_ARRAYS

!#######################################################################

      END SUBROUTINE UKCA_FASTJX
