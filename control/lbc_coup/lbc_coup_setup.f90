! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Setup for LBC coupling
SUBROUTINE lbc_coup_setup(ltimer_lbc)

  USE um_input_control_mod, ONLY : model_domain
  USE domain_params, ONLY: mt_global

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE filenamelength_mod, ONLY :                                        & 
      filenamelength
  USE ereport_mod, ONLY : ereport
  USE PrintStatus_mod
  USE UM_ParVars
  USE Control_Max_Sizes
  ! This makes sure um_sleep gets compiled from portio
  USE io_dependencies
  IMPLICIT NONE


! Description
!  Does the setup in preparation for LBC coupling.

! Method
!  The main function of this coupling code is to allow the driving model
! (generating the LBCs) and the LAM model (receiving the LBCs) to run in
! parallel. In a nutshell, this coupling ensures that the LBCs are
! available every time it reaches the stage when it needs to read in more
! LBCs. Communication between the global and LAM was done through a
! communication file - which the global wrote to and the LAM read from it.

! This functionality was used for a number of years in the operational
! suite to enable the global and the old UK mesoscale run alongside rather
! than the UK Mes wait for the global to complete its run. It stopped being
! used in the oper suite during 2006 when the NAE took over the old
! UK Mes, and the global/NAE main runs are not run in parallel. The
! global/NAE update runs are run in parallel but the NAE uses a LBC file
! from an earlier global run so no need for coupling.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC_coup

! Code Description:
!   Language: FORTRAN 90 + common extensions.
!   This code is written to UMDP3 v6 programming standards.


  LOGICAL,INTENT(IN) :: ltimer_lbc  ! Timer logical

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
  !L       Control for boundary updating.
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

  INTEGER           :: um_lbc_coup ! LBC Coupling Switch : 1/0 is on/off
  CHARACTER(LEN=filenamelength) ::                                      &
                       filename    ! Filename of communication file.
  CHARACTER(LEN=8)  :: c_lbc_coup  ! Character variable to read env var
  CHARACTER(LEN=8)  :: ch_date2    ! Date returned from date_and_time
  CHARACTER(LEN=10) :: ch_time2    ! Time returned from date_and_time
  CHARACTER(LEN=5)  :: runid_char  ! RUNID for job
  CHARACTER(LEN=4)  :: runtype_char! Run Type (ie. NRUN, CRUN)
  INTEGER           :: lbc_ntimes  ! No of BCs in communication file.
  INTEGER           :: len_wait_tot! Total wait for availability of BCs
  INTEGER           :: um_lbc_wait_usec ! Total number of microseconds to sleep.
  INTEGER           :: j           ! Loop counter

  LOGICAL           :: l_exist     !  T : Communication File exists
  LOGICAL           :: l_active    !  T : Output stream active for LBCs.

  ! Error reporting
  INTEGER           :: icode       ! =0 normal exit; >0 error exit
  CHARACTER(LEN=256):: cmessage    ! Error message
  CHARACTER(LEN=*)        routinename

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle
  PARAMETER (routinename='LBC_SETUP_FILES')

! --------------------- COMDECK LBC_COUP --------------------------
!    Description:
!       This COMDECK stores the variables connected with the
!       Parallel Running between the Global and Mesoscale models.
!       The Mesoscale has to be held back until there are sufficient
!       Boundary Conditions (BCs) to proceed.
!
      logical l_lbc_coup        ! T : Global/Meso Coupling on
      integer um_lbc_stream     ! Output Stream generating BCs.
      
      ! UM 6.5 -  MODEL_ANALYSIS_HRS changed to REAL - 
      !                requires LBC_FC_HRS to REAL also
      real    lbc_fc_hrs        ! Forecast time w.r.t analysis time
      integer um_lbc_wait       ! Wait time between re-tries if BCs
                                ! not available.
      integer um_lbc_wait_max   ! Maximum wait time.
      CHARACTER(LEN=256) lbc_filename ! Name of file that communicates between
                                 ! Global and Meso.

      COMMON /LBC_COUP/ l_lbc_coup, um_lbc_stream, lbc_fc_hrs,          &
     &                  um_lbc_wait, um_lbc_wait_max, lbc_filename
! ----------------------------------------------------------------------

  ! Find out if LBC Coupling has been switched on this run
  ! from the env. variable UM_LBC_COUP.

  IF (lhook) CALL dr_hook('LBC_COUP_SETUP',zhook_in,zhook_handle)
  CALL fort_get_env('UM_LBC_COUP',11,c_lbc_coup,8,icode)
  IF (icode /= 0) THEN
     um_lbc_coup = 0 ! No coupling
     IF (printstatus  >=  prstatus_normal) THEN
        WRITE (6,*) ' Env Var UM_LBC_COUP not set.'
        WRITE (6,*) ' Setting UM_LBC_COUP to ',um_lbc_coup
     END IF
  ELSE
     READ(c_lbc_coup,'(i8)') um_lbc_coup
     IF (printstatus  >=  prstatus_normal) THEN
        WRITE (6,*) ' UM_LBC_COUP is set to ',um_lbc_coup
     END IF
  END IF
  IF (um_lbc_coup == 0 .OR. um_lbc_coup == 1) THEN
     l_lbc_coup = um_lbc_coup == 1
  ELSE
     WRITE (6,*) ' Invalid value given to UM_LBC_COUP ',                       &
          um_lbc_coup
     WRITE (6,*) ' Valid values are 0 or 1'
     WRITE (6,*) ' L_LBC_COUP set to F. No LBC Coupling ',                     &
          'in this run.'
     cmessage = 'U_MODEL : Invalid value given to UM_LBC_COUP'
     icode = 100

     CALL ereport(routinename,icode,cmessage)
  END IF

  IF (printstatus  >=  prstatus_normal) THEN
     IF (l_lbc_coup) THEN
        WRITE (6,*) ' LBC COUPLING switched on in this run.'
     ELSE
        WRITE (6,*) ' LBC COUPLING switched off in this run.'
     END IF
  END IF

  IF (l_lbc_coup) THEN

     IF (model_domain == mt_global) THEN

     ! Find out which LBC output stream is providing the data
     ! from the env. variable UM_LBC_STREAM.

     CALL fort_get_env('UM_LBC_STREAM',13,c_lbc_coup,8,icode)
     IF (icode /= 0) THEN
        um_lbc_stream = 0 ! No coupling
        IF (printstatus  >=  prstatus_min) THEN
           WRITE (6,*) ' gl : Env Var UM_LBC_STREAM not set.'
           WRITE (6,*) ' gl : Setting UM_LBC_STREAM to ',um_lbc_stream
        END IF
     ELSE
        READ(c_lbc_coup,'(i8)') um_lbc_stream
        IF (printstatus  >=  prstatus_normal) THEN
           WRITE (6,*) ' gl : UM_LBC_STREAM is set to ',um_lbc_stream
        END IF
     END IF

     ! Check validity of UM_LBC_STREAM

     IF (um_lbc_stream <  1.OR.um_lbc_stream >  max_n_intf_a) THEN
        WRITE (6,*) ' gl : UM_LBC_STREAM = ',um_lbc_stream,                    &
             ' is an invalid value.'
        WRITE (6,*) ' gl : Valid values are 1-',max_n_intf_a
        cmessage = 'U_MODEL : Invalid value given to UM_LBC_STREAM'
        icode = 101

        CALL ereport(routinename,icode,cmessage)
     END IF

     ! Check if this output stream is active.
     l_active = .FALSE.
     DO j=1,n_intf_a
        l_active = l_active .OR. um_lbc_stream == lbc_stream_a(j)
     END DO
     IF (.NOT.l_active) THEN
        WRITE (6,*) ' gl : Output LBC stream ',um_lbc_stream,                  &
             ' is inactive. Check UM_LBC_STREAM.'
        WRITE (6,*) ' gl : Active LBC streams are ',                           &
             (lbc_stream_a(j),j=1,n_intf_a)
        cmessage = 'U_MODEL : Output LBC stream is inactive.'
        icode = 101

        CALL ereport(routinename,icode,cmessage)
     END IF

     ELSE

     ! Find out how long the mesoscale is to wait if there
     ! are insufficient boundary conditions to proceed.

     CALL fort_get_env('UM_LBC_WAIT',11,c_lbc_coup,8,icode)
     IF (icode /= 0) THEN
        um_lbc_wait = 0 ! No waiting
        IF (printstatus  >=  prstatus_min) THEN
           WRITE (6,*) ' ms : Env Var UM_LBC_WAIT not set.'
           WRITE (6,*) ' ms : Setting UM_LBC_WAIT to ',um_lbc_wait
        END IF
     ELSE
        READ(c_lbc_coup,'(i8)') um_lbc_wait
        IF (printstatus  >=  prstatus_normal) THEN
           WRITE (6,*) ' ms : UM_LBC_WAIT is set to ',um_lbc_wait
        END IF
     END IF

     ! Find out maximum wait if there are insufficient
     ! boundary conditions to proceed.

     CALL fort_get_env('UM_LBC_WAIT_MAX',15,c_lbc_coup,8,icode)
     IF (icode /= 0) THEN
        um_lbc_wait_max = 0 ! No waiting
        IF (printstatus  >=  prstatus_min) THEN
           WRITE (6,*) ' ms : Env Var UM_LBC_WAIT_MAX not set.'
           WRITE (6,*) ' ms : Setting UM_LBC_WAIT_MAX to ',                    &
                um_lbc_wait_max
        END IF
     ELSE
        READ(c_lbc_coup,'(i8)') um_lbc_wait_max
        IF (printstatus  >=  prstatus_normal) THEN
           WRITE (6,*) ' ms : UM_LBC_WAIT_MAX is set to ',um_lbc_wait_max
        END IF
     END IF
     END IF  ! if GLOBAL

  END IF  !  if l_lbc_coup

  icode=0  ! (incase FORT_GET_ENV has left with non-zero value)
  IF (model_domain == mt_global) THEN

  IF (l_lbc_coup) THEN

     IF (mype == 0) THEN

        ! Get filename attached to Unit 190
        CALL get_file(190,filename,filenamelength,icode)

        IF (icode /= 0) THEN
           WRITE (6,*) ' gl : Problem with GET_FILE',                          &
                ' for Unit No 190.'
           WRITE (6,*) ' gl : Return code from GET_FILE ',icode
           WRITE (cmessage,*)                                                  &
                'U_MODEL : Error in GET_FILE for Unit No 190.'
           icode = 102

           CALL ereport(routinename,icode,cmessage)
        END IF

        IF (printstatus  >=  prstatus_normal) THEN
           WRITE (6,*) ' gl : Filename for unit no 190 ',filename
        END IF

        ! Open the file with WRITE permission.
        OPEN(UNIT=190,FILE=filename,ACTION="write",IOSTAT=icode)

        IF (icode /= 0) THEN
           WRITE (6,*) ' gl : Problem with OPEN for Unit 190.'
           WRITE (6,*) ' gl : Return code from OPEN ',icode
           WRITE (cmessage,*)                                                  &
                'U_MODEL : Problem with OPEN for Unit No 190.'
           icode = 103

           CALL ereport(routinename,icode,cmessage)
        END IF

        IF (printstatus  >=  prstatus_normal) THEN
           WRITE (6,*) ' gl : File(unit 190) has been opened.'
        END IF

        ! Send info to meso that L_LBC_COUP=T in Global.

        lbc_ntimes = 1000 + um_lbc_stream
        WRITE (190,*) lbc_ntimes

! DEPENDS ON: um_fort_flush
        CALL um_fort_flush (190,icode)

        IF (icode /= 0) THEN
           WRITE (6,*) 'Return Code from FLUSH ',icode
           icode = 104
           WRITE (cmessage,*) 'U_MODEL : Error flushing out '//                &
                'contents for Unit 190.'

           CALL ereport(routinename,icode,cmessage)
        END IF

        IF (printstatus  >=  prstatus_normal) THEN
           WRITE (6,*) ' gl : lbc_ntimes ', lbc_ntimes,                        &
                ' sent to LBC_FILE.'
        END IF

        ! Unit 191 : File with information for operators to
        ! monitor progress with Boundary Conditions generated.

        ! Get filename attached to Unit 191
        CALL get_file(191,filename,filenamelength,icode)

        IF (icode /= 0) THEN
           WRITE (6,*) ' gl : Problem with GET_FILE for Unit 191.'
           WRITE (6,*) ' gl : Return code from GET_FILE ',icode
           WRITE (cmessage,*)                                                  &
                'U_MODEL : Error in GET_FILE for Unit No 191.'
           icode = 105

           CALL ereport(routinename,icode,cmessage)
        END IF

        IF (printstatus  >=  prstatus_normal) THEN
           WRITE (6,*) ' gl : Filename for unit no 191 ',filename
        END IF

        ! Open the file with WRITE permission only.
        OPEN(UNIT=191,FILE=filename,ACTION="write",IOSTAT=icode)

        IF (icode /= 0) THEN
           WRITE (6,*) ' gl : Problem with OPEN for Unit 191.'
           WRITE (6,*) ' gl : Return code from OPEN ',icode
           WRITE (cmessage,*)                                                  &
                'U_MODEL : Problem with OPEN for Unit No 191.'
           icode = 106

           CALL ereport(routinename,icode,cmessage)
        END IF

        WRITE (6,*) ' gl : File opened on unit 191.'

        ! Send RUNID and date to file on unit 191
        CALL fort_get_env('RUNID',5,runid_char,5,icode)
        IF (icode /= 0) THEN
           WRITE (6,*) ' Problem with FORT_GET_ENV for RUNID.'
           WRITE (cmessage,*)                                                  &
                'U_MODEL : Problem with FORT_GET_ENV for RUNID.'
           icode = 107

           CALL ereport(routinename,icode,cmessage)
        END IF

        CALL fort_get_env('TYPE',4,runtype_char,5,icode)
        IF (icode /= 0) THEN
           WRITE (6,*) ' Problem with FORT_GET_ENV for TYPE.'
           WRITE (cmessage,*)                                                  &
                'U_MODEL : Problem with FORT_GET_ENV for TYPE.'
           icode = 108

           CALL ereport(routinename,icode,cmessage)
        END IF

        CALL DATE_AND_TIME(ch_date2, ch_time2)
        WRITE (191,*) ' RUNID : ',TRIM(runid_char),                            &
             '  RUN TYPE : ',TRIM(runtype_char),                               &
             '  on ',ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4)

! DEPENDS ON: um_fort_flush
        CALL um_fort_flush (191,icode)

        IF (icode /= 0) THEN
           WRITE (6,*) 'Return Code from FLUSH ',icode
           icode = 109
           WRITE (cmessage,*) 'U_MODEL : Error flushing out '//                &
                'contents for Unit 191.'

           CALL ereport(routinename,icode,cmessage)
        END IF

     END IF  !  if mype == 0

  END IF  !   if l_lbc_coup

  ELSE

  IF (l_lbc_coup) THEN

! DEPENDS ON: timer
     IF (ltimer_lbc) CALL timer ('LBC_WAIT',3)

     IF (mype == 0) THEN

        ! Get filename attached to Unit 190
        CALL get_file(190,lbc_filename,filenamelength,icode)
        IF (printstatus  >=  prstatus_normal) THEN
           WRITE (6,*) ' ms : Filename from GET_FILE ',lbc_filename
        END IF

        IF (icode /= 0) THEN
           WRITE (6,*) ' Return code from GET_FILE ',icode
           icode = 600
           WRITE (cmessage,*) 'U_MODEL : Problem with GET_FILE '//             &
                'for Unit No 190.'

           CALL ereport(routinename,icode,cmessage)
        END IF

        CALL DATE_AND_TIME(ch_date2, ch_time2)

        IF (printstatus  >=  prstatus_normal) THEN
           WRITE(6,*) 'LBC_COUP: ',                                            &
                ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),' on ',      &
                ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),             &
                ' Wait to see if file exists.'
        END IF

        len_wait_tot = 0
149     CONTINUE

        ! Check that the file exists.
        INQUIRE (FILE=lbc_filename,EXIST=l_exist,IOSTAT=icode)

        IF (l_exist) THEN  !  file exists

           CALL DATE_AND_TIME(ch_date2, ch_time2)
           IF (printstatus  >=  prstatus_normal) THEN
              WRITE(6,*) 'LBC_COUP: ',                                         &
                   ch_time2(1:2),':',ch_time2(3:4),':',ch_time2(5:6),' on ',   &
                   ch_date2(7:8),'/',ch_date2(5:6),'/',ch_date2(1:4),          &
                   ' File exists - Proceed to open file.'
           END IF

           ! Open the file with READ ONLY permission.
           OPEN (UNIT=190,FILE=lbc_filename,ACTION="read",                     &
                IOSTAT=icode)

           ! Check return code from OPEN.
           IF (icode /= 0) THEN
              WRITE (6,*) ' Return code from OPEN ',icode
              icode = 601
              WRITE (cmessage,*) 'U_MODEL : Problem with OPEN '//              &
                   'for Unit No 190.'

              CALL ereport(routinename,icode,cmessage)
           END IF

        ELSE  !  file does not exist

           IF (len_wait_tot >= um_lbc_wait_max) THEN

              ! Maximum wait time has been reached/exceeded.

              WRITE (6,*) ' ms : lbc_file does not exist.'
              WRITE (6,*) ' ms : Maximum wait time reached',                   &
                   ' after ',um_lbc_wait_max,' seconds.'
              icode = 602
              cmessage = 'U_MODEL : LBC_FILE does not exist.'

              CALL ereport(routinename,icode,cmessage)

           ELSE

              ! Wait for um_lbc_wait seconds before another attempt
              ! to see if file exists.

              IF (printstatus  >=  prstatus_normal) THEN
                 WRITE (6,*) ' ms : lbc_file does not exist yet.'
                 WRITE (6,*) ' ms : wait for ',um_lbc_wait,                    &
                      ' seconds and retry.'
              END IF
              um_lbc_wait_usec = 1000000*um_lbc_wait
              CALL um_sleep(um_lbc_wait_usec)
              len_wait_tot = len_wait_tot+um_lbc_wait
              IF (printstatus  >=  prstatus_normal) THEN
                 WRITE (6,*) ' ms : Total Wait so far ',len_wait_tot,          &
                      ' seconds.'
              END IF
              GO TO 149   !  Retry to see if LBC_FILE exists

           END IF

        END IF  ! if l_exist

     END IF  !  if mype=0


     ! Check that LBC_COUPLING has been switched on in Global.
     IF (mype == 0) THEN

        len_wait_tot = 0
150     CONTINUE

        ! Close the communication file and re-open
        CLOSE(190)
        OPEN (190,FILE=lbc_filename,ACTION="read",IOSTAT=icode)

        ! Check retrun code from OPEN
        IF (icode /= 0) THEN
           WRITE (6,*) ' Return code from OPEN ',icode
           icode = 603
           WRITE (cmessage,*) 'U_MODEL : Problem with OPEN '//                 &
                'for Unit No 190.'

           CALL ereport(routinename,icode,cmessage)
        END IF

        ! Read in the first value
        READ (190,*,IOSTAT=icode) lbc_ntimes

        ! Check return code from READ
        IF (icode /= 0) THEN

           IF (printstatus  >=  prstatus_normal) THEN
              WRITE (6,*) ' ms : Return code from READ ',icode
           END IF

           IF (len_wait_tot >= um_lbc_wait_max) THEN

              ! Maximum wait time has been reached or exceeded.
              ! Give up waiting and abort.

              WRITE (6,*) ' ms : Required LBC_NTIMES not read in',             &
                   ' after ',um_lbc_wait_max,' seconds.'
              icode = 604
              cmessage = 'U_MODEL : Required LBC_NTIMES '//                    &
                      'not found in LBC_FILE.'

              CALL ereport(routinename,icode,cmessage)

           ELSE

              ! Wait for um_lbc_wait seconds abefore another attempt
              ! to read a value.

              IF (printstatus  >=  prstatus_normal) THEN
                 WRITE (6,*) ' ms : wait for ',um_lbc_wait,                    &
                      ' seconds and retry.'
              END IF
              um_lbc_wait_usec = 1000000*um_lbc_wait
              CALL um_sleep(um_lbc_wait_usec)

              len_wait_tot = len_wait_tot+um_lbc_wait
              IF (printstatus  >=  prstatus_normal) THEN
                 WRITE (6,*) ' ms : Total Wait so far ',len_wait_tot,          &
                      ' seconds.'
              END IF
              GO TO 150   !  Retry to see if required LBC_NTIMES exists

           END IF

        END IF  !  if icode /= 0

        ! The first value in the file is >1000.
        IF (lbc_ntimes >  1000) THEN
           um_lbc_stream = lbc_ntimes - 1000
           IF (printstatus  >=  prstatus_normal) THEN
              WRITE (6,*) ' ms : l_lbc_coup = T in Global'
              WRITE (6,*) ' ms : global output stream is ',um_lbc_stream
           END IF
        END IF  !  if l_lbc_ntimes

     END IF  !  if mype == 0

! DEPENDS ON: timer
     IF (ltimer_lbc) CALL timer ('LBC_WAIT',4)

  END IF  !  if l_lbc_coup

  END IF  !  if GLOBAL

  IF (lhook) CALL dr_hook('LBC_COUP_SETUP',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE lbc_coup_setup
