! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Main routine for calculating online photolysis rates
!   using fast-j
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
        SUBROUTINE UKCA_FASTJ(                                         &
! ARGLNDM Constants for physics routines
        land_index, land_ice_index, soil_index,                         &
! ARGLNDM end
        dj,                                                            &
        p_rho_levels, p0,                                              &
        p_theta_levels,                                                &
        t, tstar,                                                      &
        so4_aitken, so4_accum,                                         &
        q, qcl, qcf, area_cloud_fraction,                              &
        conv_cloud_lwp, conv_cloud_top, conv_cloud_base,               &
        conv_cloud_amount,                                             &
        surf_alb,                                                      &
        nd_o3, um_ozone3d,                                             &
        land_frac_ctile,                                               &
        z_top_of_model,                                                &
        l_moses_ii)


      USE FASTJ_DATA,        ONLY: fastj_set_limits,                   &
                                fastj_allocate_memory,                 &
                                fastj_deallocate_memory,               &
                                Blocking,kpcx,jpcl,                    &
                                jjpnl,jppj,nsl,sa,p,rz_3d,             &
                                tau,month,iday,nslon,nslat,            &
                                sza,sza_2d,SZAFAC,SZAFAC_2d,u0,        &
                                od,odw,odi,sulphur,fj_ozone,           &
                                t_fastj
      USE FASTJ_MIE,         ONLY: nwb3,nwb2,fastj_mie_alloc_mem,      &
                                   fastj_mie_dealloc_mem
      USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels
      USE trignometric_mod,  ONLY: true_longitude
      USE dyn_coriolis_mod,  ONLY: f3_at_u
      USE spec_sw_lw,        ONLY: sw_spectrum, spectrum
      USE rad_pcf,           ONLY: ip_aerosol_param_moist,             &
                                   ip_accum_sulphate,                  &
                                   ip_aitken_sulphate
      USE UKCA_CONSTANTS,    ONLY: c_o3

      USE earth_constants_mod, ONLY: g, two_omega

      USE water_constants_mod, ONLY: tm

      USE conversions_mod,   ONLY: pi_over_180

      USE parkind1,          ONLY: jprb, jpim
      USE yomhook,           ONLY: lhook, dr_hook
      USE PrintStatus_mod
      USE UM_ParVars
      USE Control_Max_Sizes
      USE um_input_control_mod, ONLY: lcal360
      USE cloud_inputs_mod,     ONLY: l_pc2
      USE cv_run_mod,           ONLY: l_3d_cca
      USE submodel_mod
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


      REAL, DIMENSION(row_length,rows,model_levels,jppj), INTENT(INOUT) &
                                           :: dj
      REAL, DIMENSION(row_length,rows,model_levels+1), INTENT(IN)       &
                                           :: p_rho_levels
      REAL, DIMENSION(row_length,rows,model_levels), INTENT(IN)        &
                                           :: p_theta_levels
      REAL, DIMENSION(row_length,rows,model_levels)     :: t
      REAL, DIMENSION(row_length,rows,model_levels), INTENT(IN)        &
                                           :: so4_aitken
      REAL, DIMENSION(row_length,rows,model_levels), INTENT(IN)        &
                                           :: so4_accum
      REAL, DIMENSION(row_length,rows,wet_levels), INTENT(IN)          &
                                           :: q, qcl, qcf
      REAL, DIMENSION(row_length,rows,wet_levels), INTENT(IN)          &
                                           :: area_cloud_fraction
      REAL, DIMENSION(row_length,rows,wet_levels), INTENT(IN)          &
                                           :: conv_cloud_amount
      REAL, DIMENSION(row_length,rows,ozone_levels), INTENT(IN)        &
                                           :: um_ozone3d ! MMR
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: p0
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: tstar
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: conv_cloud_lwp
      REAL, DIMENSION(row_length,rows), INTENT(IN) :: surf_alb
      REAL, DIMENSION(land_points), INTENT(IN) :: land_frac_ctile
      INTEGER, DIMENSION(row_length,rows), INTENT(IN)                  &
                                           :: conv_cloud_top
      INTEGER, DIMENSION(row_length,rows), INTENT(IN)                  &
                                           :: conv_cloud_base


      REAL, INTENT(IN) :: z_top_of_model
      INTEGER, INTENT(IN) :: nd_o3
      LOGICAL, INTENT(IN) :: L_moses_ii

! Local variables

      INTEGER                           :: i_humidity
      INTEGER                           :: i,j,k,l,ia
      INTEGER                           :: ni,nj,ml
      INTEGER                           :: nr_light_points

      REAL, PARAMETER                   :: pi180=pi_over_180
      REAL, PARAMETER                   :: ADIFC = 0.06
      REAL                              :: d_eff
      REAL                              :: delta_humidity
      REAL                              :: f
      REAL                              :: timej
      REAL, DIMENSION(row_length,rows)  :: land_frac
      REAL, DIMENSION(row_length, rows) :: sin_true_latitude
      REAL, DIMENSION(row_length,rows,                                 &
                      model_levels)     :: qsat
      REAL, DIMENSION(row_length,rows,                                 &
                      model_levels)     :: rh


      LOGICAL, SAVE                         :: first=.true.
      LOGICAL, DIMENSION(row_length, rows)  :: is_day

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!     interface fastj routines

      INTERFACE 
        SUBROUTINE FASTJ_INPHOT
         USE     FASTJ_DATA,         ONLY: kpcx,dtaumax,dtausub,dsubdiv,&
                                           szamax,rad,zzht
         USE     FASTJ_SPECS,        ONLY: FASTJ_RD_JS, FASTJ_RD_TJPL
         IMPLICIT NONE
        END SUBROUTINE FASTJ_INPHOT

        SUBROUTINE FASTJ_SOLAR2(timej,sinlat,longitude,lcal360)
          USE   FASTJ_DATA,    ONLY: iday,tau,szamax,sza_2d,szafac_2d
          IMPLICIT NONE
          REAL,INTENT(IN)                :: timej
          REAL,DIMENSION(:,:),INTENT(IN) :: sinlat,longitude
          LOGICAL,INTENT(IN)             :: lcal360
        END SUBROUTINE FASTJ_SOLAR2

        SUBROUTINE FASTJ_PHOTOJ(zpj,timej)
          USE FASTJ_DATA
          IMPLICIT NONE
          REAL,DIMENSION(:,:,:,:),INTENT(OUT)         :: zpj
          REAL,INTENT(IN)                             :: timej
        END SUBROUTINE FASTJ_PHOTOJ
      END INTERFACE
!

      IF (lhook) CALL dr_hook('UKCA_FASTJ',zhook_in,zhook_handle)
      ni = row_length
      nj = rows
      ml = model_levels

      CALL FASTJ_SET_LIMITS(row_length,rows,model_levels)
      CALL FASTJ_ALLOCATE_MEMORY 
      CALL FASTJ_MIE_ALLOC_MEM (kpcx)

      IF (first) THEN
! DEPENDS ON: fastj_inphot
        CALL FASTJ_INPHOT
        first=.false.
      ENDIF

      CALL FASTJ_COMPUTE_PHOTJ_VALUES

!     set variables concerning model time

      timej              = secs_per_stepim(a_im)/3600.
      iday               = i_day_number
      month              = int(real(iday)*12.0/365.0)+1    !  Approximately
      tau                = i_hour*1.+i_minute/60.+i_second/3600.

      sa = surf_alb
      IF(Blocking%Mode <= 3)   THEN          !not in load balancing mode
        t_fastj(:,:,1:ml)  = t
        p                  = p_rho_levels/100.0  ! Pa to mbar
      END IF
      
!     Convert ozone mmr to vmr and write to fj_ozone array

      fj_ozone(:,:,1:ml) = um_ozone3d/c_o3
      fj_ozone(:,:,ml+1) = fj_ozone(:,:,ml)

! DEPENDS ON: fastj_solar2
      CALL FASTJ_SOLAR2 (timej,sin_true_latitude,true_longitude,       &
                         lcal360)

      IF (printstatus >= prstatus_diag) THEN 
        WRITE(6,*) 'so4_aitken: ',shape(so4_aitken),                     &
                 maxval(so4_aitken),minval(so4_aitken),minloc(so4_aitken)
        WRITE(6,*) 'so4_accum:  ',shape(so4_accum),                      &
                    maxval(so4_accum),minval(so4_accum)
        WRITE(6,*) 'sulphur:  ',shape(sulphur),                          &
                    maxval(sulphur),minval(sulphur)
        WRITE(6,*) 'sa:       ',shape(sa),                               &
                    maxval(sa),minval(sa)
        WRITE(6,*) 'p_rho_levels: ',shape(p_rho_levels),                 &
                    maxval(p_rho_levels),minval(p_rho_levels),           &
                    minloc(p_rho_levels)
        WRITE(6,*) 'fj_ozone: ',shape(fj_ozone),                         &
                    maxval(fj_ozone),minval(fj_ozone),minloc(fj_ozone)
        WRITE(6,*) 'odw:      ',shape(odw),                              &
                    maxval(odi),maxloc(odw),minval(odw)
        WRITE(6,*) 'odi:      ',shape(odi),                              &
                    maxval(odi),maxloc(odi),minval(odi)
        WRITE(6,*) 'od:       ',shape(od),                               &
                    maxval(od),maxloc(od),minval(od)
      ENDIF
      
      nslon = 1
      SELECT CASE (Blocking%Mode)

       CASE(1)
!       Blocking one row

        do j=1,rows
          SZA    = SZA_2d(:,j)
          SZAFAC = SZAFAC_2d(:,j)
          U0     = COS(sza*pi180)

          nslat  = j
          do i=1,row_length
            nsl(1,i)=nslon+(i-1)
            nsl(2,i)=nslat
          enddo

! DEPENDS ON: fastj_photoj
          CALL FASTJ_PHOTOJ (dj,timej)

        end do

       CASE(2)
!       Blocking whole domain

        nslat = 1
        do j=1,rows
          k = (j-1)*row_length+1
          SZA(k:k+rows-1)    = SZA_2d(:,j)
          SZAFAC(k:k+rows-1) = SZAFAC_2d(:,j)
          do i=1,row_length
            l=i+(j-1)*row_length
            nsl(1,l)=nslon+(i-1)
            nsl(2,l)=nslat+(j-1)
          end do
        enddo
        U0 = COS(sza*pi180)

! DEPENDS ON: fastj_photoj
        CALL FASTJ_PHOTOJ (dj,timej)

       CASE(3)
        l     = 0
        nslat = 1
        nr_light_points = 0.0
        do j=1,rows
          do i=1,row_length
            if( SZAFAC_2d(i,j)  > 0.001e0) then
              l = l+1
              SZA(l)    = SZA_2d(i,j)
              SZAFAC(l) = SZAFAC_2d(i,j)
              nsl(1,l)  = nslon+(i-1)
              nsl(2,l)  = nslat+(j-1)
            end if

            if(l == Blocking%BL_size  .OR.                             &
              (j == rows .and. i == row_length .and. l > 0))   then

              nr_light_points = nr_light_points+l
              kpcx = l
              jjpnl = jpcl*kpcx
              NWB3 = NWB2*kpcx
              u0(1:l) = COS(sza(1:l)*pi180)
! DEPENDS ON: fastj_photoj              
              CALL FASTJ_PHOTOJ (dj,timej)
              l = 0
            end if
          end do
        end do
        l = count(SZAFAC_2d > 0.001e0)

       CASE(4)
        nr_light_points = count(SZAFAC_2d > 0.001e0)
        CALL FASTJ_LOADBALANCE_PHOTOJ
      END SELECT

      IF (printstatus >= prstatus_diag) THEN 
        WRITE(6,*) 'tau: ',tau
        WRITE(6,*) 'dj:        ',sum(dj),maxval(dj),maxloc(dj)
      ENDIF
      CALL FASTJ_MIE_DEALLOC_MEM 
      CALL FASTJ_DEALLOCATE_MEMORY

      IF (lhook) CALL dr_hook('UKCA_FASTJ',zhook_out,zhook_handle)
      RETURN

!     internal subroutines

      CONTAINS

! ######################################################################
      SUBROUTINE FASTJ_COMPUTE_PHOTJ_VALUES
        USE PrintStatus_mod
        USE UM_ParVars
        USE Control_Max_Sizes
        IMPLICIT NONE
        INTEGER, PARAMETER                    :: sw_band_aer = 4
        INTEGER                               :: i,j,k,l

        REAL,DIMENSION(1:row_length,1:rows)                :: total_mass 
        REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: d_mass
        REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: sulph_accu
        REAL,DIMENSION(1:row_length,1:rows,1:model_levels) :: sulph_aitk
        REAL,DIMENSION(1:row_length,1:rows,1:model_levels+1) :: pj

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0 
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

        IF (lhook) CALL dr_hook('FASTJ_COMPUTE_PHOTJ_VALUES',zhook_in,zhook_handle)
        DO k=1,model_levels
          CALL qsat_wat(qsat(1:row_length,1:rows,k),                   &
                         t(1:row_length,1:rows,k),                     &
            p_theta_levels(1:row_length,1:rows,k),theta_field_size)
        ENDDO
        rh=q(1:row_length,1:rows,:)*(1.0-qsat)/                        &
         (qsat*(1.0-q(1:row_length,1:rows,:)))
        rh=MIN(rh,0.99) ! limit humidity to <= 99%
        rh=MAX(rh,0.01) ! limit humidity to >=  1%


! Expand land_fraction
        land_frac=0.
        DO l=1,land_points
          j=(land_index(l)-1)/row_length+1
          i=land_index(l)-(j-1)*row_length
          land_frac(i,j)=land_frac_ctile(l)
        END DO

        sin_true_latitude = f3_at_u(1:ni,1:nj) / two_omega

! rz in cm
        rz_3d(:,:,1)    = r_Theta_levels(1:ni,1:nj,0)*100.
        rz_3d(:,:,2:ml) = r_rho_levels(1:ni,1:nj,2:ml)*100.
        rz_3d(:,:,ml+1) = r_Theta_levels(1:ni,1:nj,ml)*100.

        pj(:,:,1)=p0
        pj(:,:,2:model_levels)=p_rho_levels(:,:,2:model_levels)
        pj(:,:,model_levels+1)=0.

!      We multiply by 4.125 to convert from mass mixing ratio of 
!      sulphur atoms to mass mixing ratio of ammonium sulphate.
!      Use d_mass to convert to mass of ammonium sulphate.

        d_mass=(pj(:,:,1:ml)-pj(:,:,2:ml+1))/g
        sulph_accu=so4_accum*4.125*d_mass
        sulph_aitk=so4_aitken*4.125*d_mass

        odw(:,:,1:wet_levels)    = qcl(:,:,1:wet_levels)
        odi(:,:,1:wet_levels)    = qcf(:,:,1:wet_levels)
        odw(:,:,wet_levels+1:ml) = 0.0
        odi(:,:,wet_levels+1:ml) = 0.0

 
!       Account for random cloud overlap using factor suggested by Briegleb et al (1992) 
        odw(:,:,1:wet_levels)    =  odw(:,:,1:wet_levels) & 
          * (area_cloud_fraction(:,:,1:wet_levels))**(1.5) 
        odi(:,:,1:wet_levels)    =  odi(:,:,1:wet_levels) & 
          * (area_cloud_fraction(:,:,1:wet_levels))**(1.5) 
 
!************************************************************************ 

!----------------------------------------------------------------- 
! Only add the convective cloud if we are not using the PC2 cloud scheme 
        IF(.NOT.L_PC2) THEN 

          total_mass = 0.0
          DO k = 1,model_levels
            WHERE ( k <= conv_cloud_top .and. k >= conv_cloud_base) 
              total_mass = total_mass + d_mass(:,:,k)
            ENDWHERE
          END DO

          IF (l_3d_cca) THEN
            DO k = 1,n_cca_lev
              DO j = 1, rows
                DO i = 1, row_length
                  IF (conv_cloud_top(i,j) > 0) THEN
                    IF (t(i,j,k) > tm)  THEN
                      odw(i,j,k) = odw(i,j,k) +                         &
                      (conv_cloud_lwp(i,j)*conv_cloud_amount(i,j,k))/   &
                       total_mass(i,j)
                    ELSE
                      odi(i,j,k) = odi(i,j,k) +                            &
                       (conv_cloud_lwp(i,j)*conv_cloud_amount(i,j,k))/     &
                       total_mass(i,j)
                    END IF
                  END IF
                END DO
              END DO
            END DO
          ELSE
            DO k=1,model_levels
              DO j = 1, rows
                DO i = 1, row_length
                  IF ( k <= conv_cloud_top(i,j)                           &
                          .and. k >= conv_cloud_base(i,j))   THEN
                    IF(t(i,j,k)>tm)  THEN
                      odw(i,j,k) = odw(i,j,k)+(conv_cloud_lwp(i,j)*       &
                          conv_cloud_amount(i,j,1))/total_mass(i,j)
                    ELSE
                      odi(i,j,k)= odi(i,j,k)+(conv_cloud_lwp(i,j)*         & 
                        conv_cloud_amount(i,j,1))/total_mass(i,j)
                    END IF
                  END IF
                END DO
              END DO
            END DO
          END IF

        ENDIF ! Not PC2 cloud scheme

!       Convert mass mixing ratios to column densities.

        odw = odw*d_mass
        odi = odi*d_mass

!       set effective radii for water drops, different for land&sea
!       set effective diameter for ice crystals
!       Formulae from John Edwards to convert from column densities to
!       optical depths

        d_eff=100.0 ! in microns
        DO k=1,model_levels
          WHERE (land_frac>0.5) 
           odw(:,:,k) = odw(:,:,k)*(-8.86964+1.67373E3/6.)
          ELSEWHERE
           odw(:,:,k) = odw(:,:,k)*(-8.86964+1.67373E3/12.)
          ENDWHERE
        END DO
        odi = odi*(-2.189E-3+3.311E3/d_eff+3.611/d_eff**2)

        DO ia=1,sw_spectrum(1)%n_aerosol
          IF (sw_spectrum(1)%type_aerosol(ia) == ip_accum_sulphate .OR. &
               sw_spectrum(1)%type_aerosol(ia)==ip_aitken_sulphate) THEN
            IF (sw_spectrum(1)%i_aerosol_parametrization(ia) ==         &
                   ip_aerosol_param_moist) THEN ! moist aerosol
              delta_humidity = sw_spectrum(1)%humidities(2,ia)-         &
                                sw_spectrum(1)%humidities(1,ia)
              DO k=1,model_levels
               DO j = 1, rows
                DO i = 1, row_length
                 i_humidity = INT(rh(i,j,k)/delta_humidity)+1
                 f=(rh(i,j,k)-sw_spectrum(1)%humidities(i_humidity,ia))/&
                    delta_humidity
                 od(i,j,k)=                                             &
                 (sw_spectrum(1)%aerosol_absorption(i_humidity,ia,      &
                                                      sw_band_aer)+     &
                  sw_spectrum(1)%aerosol_scattering(i_humidity,ia,      &
                                                      sw_band_aer))*f+  &
                 (sw_spectrum(1)%aerosol_absorption(i_humidity+1,ia,    &
                                                        sw_band_aer)+   &
                  sw_spectrum(1)%aerosol_scattering(i_humidity+1,ia,    &
                                                  sw_band_aer))*(1-f)
                END DO
               END DO
              END DO
            ELSE ! dry aerosol
              DO k=1,model_levels
               DO j = 1, rows
                DO i = 1, row_length
                 od(i,j,k) =                                            &
                   sw_spectrum(1)%aerosol_absorption(1,ia,sw_band_aer)+ &
                   sw_spectrum(1)%aerosol_scattering(1,ia,sw_band_aer)
                END DO
               END DO
              END DO
            END IF
            IF (sw_spectrum(1)%type_aerosol(ia) == ip_accum_sulphate)   &
              sulph_accu=sulph_accu*od
            IF (sw_spectrum(1)%type_aerosol(ia) == ip_aitken_sulphate)  &
              sulph_aitk=sulph_aitk*od
          END IF
        END DO

        sulphur = sulph_aitk + sulph_accu

      IF (lhook) CALL dr_hook('FASTJ_COMPUTE_PHOTJ_VALUES',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FASTJ_COMPUTE_PHOTJ_VALUES

! ######################################################################
      SUBROUTINE FASTJ_LOADBALANCE_PHOTOJ

        USE mpl, ONLY :     &
          MPL_ADDRESS_KIND, &
          MPL_REAL,         &
          MPL_INFO_NULL

        USE PrintStatus_mod
        USE UM_ParVars
        USE Control_Max_Sizes
        IMPLICIT NONE

!--  local variables

        INTEGER                                :: i,j,l,k,m,ii,istat
        INTEGER                                :: is,ie,ian,kr
        INTEGER                                :: np_new,np_old,me
        INTEGER                                :: vals_per_gp
        INTEGER                                :: vals_per_gp_f
        INTEGER                                :: vals_per_gp_b
        INTEGER                                :: np_trans
        INTEGER                                :: npes
        INTEGER                                :: ll
        INTEGER                                :: my_comm
        INTEGER,SAVE                           :: win_b
        
        INTEGER(kind=MPL_ADDRESS_KIND)         :: winsize
        INTEGER(kind=MPL_ADDRESS_KIND)         :: disp

        INTEGER,DIMENSION(0:nproc-1)           :: npoints
        INTEGER,DIMENSION(0:nproc-1,0:nproc-1) :: sr_matrix
        
        LOGICAL                                :: sender_PE,receiver_PE
        
        REAL,DIMENSION(ni*nj)                  :: sa_lb
        REAL,DIMENSION(ni*nj)                  :: SZA_lb,SZAFAC_lb
        REAL,DIMENSION(ni*nj,ml)               :: t_lb,sulphur_lb
        REAL,DIMENSION(ni*nj,ml)               :: odw_lb,odi_lb
        REAL,DIMENSION(ni*nj,ml+1)             :: p_lb,ozone_lb,rz_lb
        REAL,DIMENSION(ni*nj,ml,jppj)          :: dj_lb

!       Global array for MPI-2 one sided communication
        REAL,DIMENSION(:,:),allocatable        :: rbg

!CDIR GM_ARRAY(rbg)
  
!--     Interfaces 

        INTERFACE FASTJ_CALC_NEW_NUMBER_OF_POINTS
          SUBROUTINE FASTJ_CALC_NEW_NUMBER_OF_POINTS (                 &
                np,npes,npoints,sr_matrix,me,np_new,sender_PE)


      USE chsunits_mod, ONLY : nunits

            IMPLICIT NONE
            INTEGER,INTENT(IN)                       :: np
            INTEGER,INTENT(IN)                       :: npes
            INTEGER,INTENT(OUT),DIMENSION(0:npes-1)  :: npoints
            INTEGER,INTENT(OUT),                                       &
                DIMENSION(0:npes-1,0:npes-1)         :: sr_matrix
            INTEGER,INTENT(OUT)                      :: me
            INTEGER,INTENT(OUT)                      :: np_new
            LOGICAL,INTENT(OUT)                      :: sender_PE
          END SUBROUTINE FASTJ_CALC_NEW_NUMBER_OF_POINTS
        END INTERFACE FASTJ_CALC_NEW_NUMBER_OF_POINTS

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0 
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle

        IF (lhook) CALL dr_hook('FASTJ_LOADBALANCE_PHOTOJ',zhook_in,zhook_handle)
        npes = nproc

!       calculate balanced number of gridpoints

! DEPENDS ON: fastj_calc_new_number_of_points
        CALL FASTJ_CALC_NEW_NUMBER_OF_POINTS (nr_light_points,         &
            nproc, npoints,sr_matrix,me,np_new,sender_PE)
        receiver_PE = .not. sender_PE


        vals_per_gp_f = 4*ml+3*(ml+1)+3    !values per pridpoint forward
        vals_per_gp_b = jppj*ml            !backward
        vals_per_gp   = max(vals_per_gp_f,vals_per_gp_b)

!       Number of gridpoints to transport
        np_trans = maxval(abs(sr_matrix))
        np_trans = np_trans*vals_per_gp


        if(np_trans == 0)   then   !No data transfer in this case
           sender_PE   = .false.
           receiver_PE = .false.
           np_trans    = 1     !dummy windows length
        end if

!       Create RealGlobalBuffer for data to bet transported between PEs

        allocate(rbg(np_trans,0:npes-1))

!       Get communicator to use from GCOM (i.e. same communicator
!       as the rest of the model)
        call GC_Get_Communicator(MY_COMM, istat)

!       Create MPI2 windows for one-sided communication
!       to allow remote PEs to access data in window
        winsize = size(rbg)*8  ! in bytes
        call MPL_Win_create (rbg,winsize,8,MPL_INFO_NULL,              &
                                     MY_COMM,win_b, istat)
!       Initial setup
        np_old  = nr_light_points


!       First, copy and compress local arrays in large one-dimensional
!       spatial chunks and eliminate night pts

!       2-D arrays

        l = 0
        do j=1,nj
!CDIR NODEP
          do i=1,ni
            if( SZAFAC_2d(i,j)  > 0.001e0) then
              l = l+1
              sa_lb(l)     = sa(i,j)
              SZA_lb(l)    = SZA_2d(i,j)
              SZAFAC_lb(l) = SZAFAC_2d(i,j)
            end if
          end do
        end do
!       3-D arrays
        do k=1,ml
          l = 0
          do j=1,nj
!CDIR NODEP
            do i=1,ni
              if( SZAFAC_2d(i,j)  > 0.001e0) then
                l = l+1
                t_lb(l,k)       = t(i,j,k)
                odw_lb(l,k)     = odw(i,j,k)
                odi_lb(l,k)     = odi(i,j,k)
                sulphur_lb(l,k) = sulphur(i,j,k)
              end if
            end do
          end do
        end do
        do k=1,ml+1
          l = 0
          do j=1,nj
!CDIR NODEP
            do i=1,ni
              if( SZAFAC_2d(i,j)  > 0.001e0) then
                l = l+1
                p_lb(l,k)       = p_rho_levels(i,j,k)/100.0 ! Pa to mbar
                ozone_lb(l,k)   = fj_ozone(i,j,k)
                rz_lb(l,k)      = rz_3d(i,j,k)
              end if
            end do
          end do
        end do

!       synchronise data transfer into window (like barrier)
        call MPL_Win_Fence (0,win_b, istat)

!       Load balancing
!       move data from sender PEs to receiver PEs to establish equal 
!       number of columns on every PE.

        if(sender_PE)   then
!         Move data in one-dimensional buffer
!         All data is move between PEs in ONE MPI call/remote PE

          do i=0,npes-1

            if(sr_matrix(me,i) /= 0) then   !something to do
!             is and ie are start and end points of 1D array
              is  = np_old+sr_matrix(me,i)+1!+ because value is negative
              ie  = np_old
              ian = ie-is+1
   
              kr  = 1

!             Fill buffer
!             2-D arrays
              rbg(kr:kr+ian-1,i) = sa_lb(is:ie) ;          kr = kr+ian
              rbg(kr:kr+ian-1,i) = SZA_lb(is:ie) ;         kr = kr+ian
              rbg(kr:kr+ian-1,i) = SZAFAC_lb(is:ie) ;      kr = kr+ian
!             3-D arrays
              do k=1,ml
                rbg(kr:kr+ian-1,i) = t_lb(is:ie,k) ;       kr = kr+ian
                rbg(kr:kr+ian-1,i) = odw_lb(is:ie,k) ;     kr = kr+ian
                rbg(kr:kr+ian-1,i) = odi_lb(is:ie,k) ;     kr = kr+ian
                rbg(kr:kr+ian-1,i) = sulphur_lb(is:ie,k) ; kr = kr+ian
              end do
              do k=1,ml+1
                rbg(kr:kr+ian-1,i) = p_lb(is:ie,k) ;       kr = kr+ian
                rbg(kr:kr+ian-1,i) = ozone_lb(is:ie,k) ;   kr = kr+ian
                rbg(kr:kr+ian-1,i) = rz_lb(is:ie,k) ;      kr = kr+ian
              end do
              np_old = np_old-ian
            end if
          end do
        end if

        call MPL_Win_Fence (0,win_b, istat)

        if(receiver_PE)   then
          do i=0,npes-1  ! looping over remote/sender PEs
            if(sr_matrix(me,i) /= 0) then   !something to do
              is     = np_old+1
              ie     = np_old+sr_matrix(me,i)
              ian    = ie-is+1
              kr     = 1

              ll   = sr_matrix(me,i)*vals_per_gp_f
              ! me is receiver PE
              disp = me*np_trans ! displacement unit (not in bytes)
                                 ! offset in window 
              CALL MPL_Get(rbg(1,i),ll,MPL_REAL, i,disp,ll,            &
                                       MPL_REAL, win_b,istat)


!             get data from  buffer
!             2-D arrays
              sa_lb(is:ie)     = rbg(kr:kr+ian-1,i) ;     kr = kr+ian
              SZA_lb(is:ie)    = rbg(kr:kr+ian-1,i) ;     kr = kr+ian
              SZAFAC_lb(is:ie) = rbg(kr:kr+ian-1,i) ;     kr = kr+ian
!             3-D arrays
              do k=1,ml
                t_lb(is:ie,k)       = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
                odw_lb(is:ie,k)     = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
                odi_lb(is:ie,k)     = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
                sulphur_lb(is:ie,k) = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
              end do
              do k=1,ml+1
                p_lb(is:ie,k)       = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
                ozone_lb(is:ie,k)   = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
                rz_lb(is:ie,k)      = rbg(kr:kr+ian-1,i) ;  kr = kr+ian
              end do
              np_old = np_old+ian
            end if
          end do
        end if

        call MPL_Win_Fence (0,win_b, istat)

!       Call Photolysis in Block mode
!       data for an equal number of grid points is located in the
!       .._lb arrays. The data is coped in chunks of Blocking%BL_size
!       into aray locted in mo_fastj_data. 
!       photoj is called for thes chunks and the output is copied into
!       dj_lb.

        do m=1,np_new,Blocking%BL_size
           do ii=m,min(np_new,m+Blocking%BL_size-1)
             i = mod(ii-1,row_length)+1
             j = (ii-1)/row_length+1
             l = ii-m+1
             sa(i,j)   = sa_lb(ii)
             SZA(l)    = SZA_lb(ii)
             SZAFAC(l) = SZAFAC_lb(ii)
             nsl(1,l)  = i
             nsl(2,l)  = j
           end do
           do k=1,ml
             do ii=m,min(np_new,m+Blocking%BL_size-1)
               i = mod(ii-1,row_length)+1
               j = (ii-1)/row_length+1
               t_fastj(i,j,k)  = t_lb(ii,k)
               odw(i,j,k)      = odw_lb(ii,k)
               odi(i,j,k)      = odi_lb(ii,k)
               sulphur(i,j,k)  = sulphur_lb(ii,k)
             end do
           end do
           do k=1,ml+1
             do ii=m,min(np_new,m+Blocking%BL_size-1)
               i = mod(ii-1,row_length)+1
               j = (ii-1)/row_length+1
               p(i,j,k)        = p_lb(ii,k)
               fj_ozone(i,j,k) = ozone_lb(ii,k)
               rz_3d(i,j,k)    = rz_lb(ii,k)
             end do
           end do
           u0 = COS(sza*pi180)

           kpcx = l
           jjpnl = jpcl*kpcx
           NWB3 = NWB2*kpcx
! DEPENDS ON: fastj_photoj
           CALL FASTJ_PHOTOJ (dj,timej)

           DO l=1,jppj
             DO k=1,ml
               DO ii=m,min(np_new,m+Blocking%BL_size-1)
                 i = mod(ii-1,row_length)+1
                 j = (ii-1)/row_length+1
                 dj_lb(ii,k,l) = dj(i,j,k,l)
               END DO
             END DO
           END DO
        END DO

        CALL MPL_Win_Fence (0,win_b, istat)

        np_old  = nr_light_points

!       Move results back to origional PEs.

        IF (receiver_PE)   THEN

!         Move data in one-dimensional buffer

          DO i=0,npes-1

            IF(sr_matrix(me,i) /= 0) THEN   !something to do
              is     = np_old+1
              ie     = np_old+sr_matrix(me,i)
              ian    = ie-is+1

              kr  = 1

              DO m=1,jppj
                DO k=1,ml
                  rbg(kr:kr+ian-1,i) = dj_lb(is:ie,k,m) ;   kr = kr+ian
                END DO
              END DO

              np_old = np_old+ian
            END IF
          END DO
        END IF

        CALL MPL_Win_Fence (0,win_b, istat)

        IF (sender_PE)   THEN
          DO i=0,npes-1
            IF(sr_matrix(me,i) /= 0) THEN   !something to do
              is  = np_old+sr_matrix(me,i)+1!+ because value is negative
              ie  = np_old
              ian = ie-is+1

              kr     = 1

              ll   = abs(sr_matrix(me,i))*vals_per_gp_b
              disp = me*np_trans
              CALL MPL_Get(rbg(1,i),ll,MPL_REAL, i,disp,ll,            &
                                       MPL_REAL,  win_b,istat)

              DO m=1,jppj
                DO k=1,ml
                  dj_lb(is:ie,k,m) = rbg(kr:kr+ian-1,i) ;   kr = kr+ian
                END DO
              END DO
              np_old = np_old-ian
            END IF
          END DO
        END IF

!       Move data back from the 1-dimensional (spatial) buffer array
!       into the original output 2-D array.

        dj = 0.0
        DO m=1,jppj
          DO k=1,ml
            l = 0
            DO j=1,nj
              DO i=1,ni
                IF( SZAFAC_2d(i,j)  > 0.001e0) THEN
                  l = l+1
                  dj(i,j,k,m) = dj_lb(l,k,m)
                END IF
              END DO
            END DO
          END DO
        END DO
              
        CALL MPL_Win_Fence (0,win_b, istat)

!       Free MPI-2 recources

        CALL MPL_Win_free(win_b,istat)
        DEALLOCATE (rbg)

      IF (lhook) CALL dr_hook('FASTJ_LOADBALANCE_PHOTOJ',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FASTJ_LOADBALANCE_PHOTOJ

      END SUBROUTINE UKCA_FASTJ
