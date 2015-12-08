! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Computes the grid sizes for one level of LBC data for each grid
!
! Method:
!   <Say how it does it: include references to external documentation>
!   <If this routine is very complex, then include a "pseudo code"
!    description of it to make its structure and method clear>
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: LBC Output

      SUBROUTINE LBC_Grid_Sizes ( jintf, l_lbc_eg )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE Control_Max_Sizes
      USE lbc_mod
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
! Include file : parlbcs.h
!
! Must be called after parvars.h
!
! Description:
!   Contains variables in connection with generating LBCs.
!
! -----------------------------------------------------------
! Stash Codes for LBCs in Section 32 (and section 31 except tracers) 
!
      Integer, Parameter :: lbc_stashcode_orog    = 1 
      Integer, Parameter :: lbc_stashcode_u       = 2 
      Integer, Parameter :: lbc_stashcode_v       = 3 
      Integer, Parameter :: lbc_stashcode_w       = 4 
      Integer, Parameter :: lbc_stashcode_density = 5 
      Integer, Parameter :: lbc_stashcode_theta   = 6 
      Integer, Parameter :: lbc_stashcode_q       = 7 
      Integer, Parameter :: lbc_stashcode_qcl     = 8 
      Integer, Parameter :: lbc_stashcode_qcf     = 9 
      Integer, Parameter :: lbc_stashcode_exner   = 10 
      Integer, Parameter :: lbc_stashcode_u_adv   = 11 
      Integer, Parameter :: lbc_stashcode_v_adv   = 12 
      Integer, Parameter :: lbc_stashcode_w_adv   = 13 
      Integer, Parameter :: lbc_stashcode_qcf2    = 14 
      Integer, Parameter :: lbc_stashcode_qrain   = 15 
      Integer, Parameter :: lbc_stashcode_qgraup  = 16 
      Integer, Parameter :: lbc_stashcode_cf_bulk = 17 
      Integer, Parameter :: lbc_stashcode_cf_liquid = 18 
      Integer, Parameter :: lbc_stashcode_cf_frozen = 19 
      Integer, Parameter :: lbc_stashcode_murk      = 20 
      Integer, Parameter :: lbc_stashcode_free_tracer = 21 
      Integer, Parameter :: lbc_stashcode_ukca_tracer = 22 
      Integer, Parameter :: lbc_stashcode_dust_div1 = 23
      Integer, Parameter :: lbc_stashcode_dust_div2 = 24
      Integer, Parameter :: lbc_stashcode_dust_div3 = 25
      Integer, Parameter :: lbc_stashcode_dust_div4 = 26
      Integer, Parameter :: lbc_stashcode_dust_div5 = 27
      Integer, Parameter :: lbc_stashcode_dust_div6 = 28
      Integer, Parameter :: lbc_stashcode_so2      = 29
      Integer, Parameter :: lbc_stashcode_dms      = 30
      Integer, Parameter :: lbc_stashcode_so4_aitken = 31
      Integer, Parameter :: lbc_stashcode_so4_accu = 32
      Integer, Parameter :: lbc_stashcode_so4_diss = 33
      Integer, Parameter :: lbc_stashcode_nh3      = 35
      Integer, Parameter :: lbc_stashcode_soot_new = 36
      Integer, Parameter :: lbc_stashcode_soot_agd = 37
      Integer, Parameter :: lbc_stashcode_soot_cld = 38
      Integer, Parameter :: lbc_stashcode_bmass_new = 39
      Integer, Parameter :: lbc_stashcode_bmass_agd = 40
      Integer, Parameter :: lbc_stashcode_bmass_cld = 41
      Integer, Parameter :: lbc_stashcode_ocff_new = 42
      Integer, Parameter :: lbc_stashcode_ocff_agd = 43
      Integer, Parameter :: lbc_stashcode_ocff_cld = 44
      Integer, Parameter :: lbc_stashcode_nitr_acc = 45
      Integer, Parameter :: lbc_stashcode_nitr_diss = 46

! -----------------------------------------------------------
!     Data Time for LBC data
      Integer :: LBC_DT_Year
      Integer :: LBC_DT_Month
      Integer :: LBC_DT_Day
      Integer :: LBC_DT_Hour
      Integer :: LBC_DT_Min
      Integer :: LBC_DT_Sec
      Integer :: LBC_DT_DayNo

      COMMON /LBC_DT/ LBC_DT_Year, LBC_DT_Month, LBC_DT_Day,            &
         LBC_DT_Hour, LBC_DT_Min,  LBC_DT_Sec,   LBC_DT_DayNo

! -----------------------------------------------------------

!     Validity Time for LBC data
      Integer :: LBC_VT_Year
      Integer :: LBC_VT_Month
      Integer :: LBC_VT_Day
      Integer :: LBC_VT_Hour
      Integer :: LBC_VT_Min
      Integer :: LBC_VT_Sec
      Integer :: LBC_VT_DayNo

      COMMON /LBC_VT/ LBC_VT_Year, LBC_VT_Month, LBC_VT_Day,            &
         LBC_VT_Hour, LBC_VT_Min,  LBC_VT_Sec,   LBC_VT_DayNo

! -----------------------------------------------------------

      Integer, Parameter :: P_Src_Grid = 2
      Integer, Parameter :: P_LBC_Grid = 4

!     1 : Start Latitude
!     2 : Start Longitude
!     3 : Row Length
!     4 : Rows

      Real :: Src_Grid (Nfld_max, P_Src_Grid)
      Real :: LBC_Grid (Nfld_max, P_LBC_Grid)

      COMMON /LBC_Grids/ Src_Grid, LBC_Grid

! -------------------------------------------------------------

      Integer :: LBC_Global_LenRimA (Nfld_max, Nhalo_max)
      Integer :: LBC_Interp_LenRimA (Nfld_max, Nhalo_max)

      COMMON /LBC_Sizes/ LBC_Global_LenRimA, LBC_Interp_LenRimA

! -------------------------------------------------------------

      INTEGER  ::  jintf            !  Index to interface area
      LOGICAL  ::  l_lbc_eg         !  Endgame LBC?

! Local variables

      INTEGER :: lbc_row_len        ! Row length including haloes
      INTEGER :: lbc_rows           ! No of rows including haloes
      INTEGER :: Ifld, Ihalo, Iside !  Loop indices

!*---------------------------------------------------------------------
      INTEGER intf_halosize(2,Nhalo_max)

      INTEGER :: u_adjust
      INTEGER :: v_adjust
      INTEGER :: v_extra_ew &
              &, v_extra_ns &
              &, u_extra_ew &
              &, u_extra_ns

      INTEGER   ::  ErrorStatus
      CHARACTER (Len=80)           :: CMessage
      CHARACTER (Len=*), PARAMETER :: RoutineName= 'LBC_Grid_Sizes'

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!      -------------------------

      IF (lhook) CALL dr_hook('LBC_GRID_SIZES',zhook_in,zhook_handle)
      ErrorStatus = 0
      CMESSAGE=' '

!     Set up intf_halosize for this area
      intf_halosize(1,1)=1
      intf_halosize(2,1)=1

      intf_halosize(1,2)=intf_exthalo_ew(jintf)
      intf_halosize(2,2)=intf_exthalo_ns(jintf)

      intf_halosize(1,3)=0
      intf_halosize(2,3)=0

!     Setup offset for fields
      IF (l_lbc_eg) THEN
        ! For ND we have one less column for U, but for EG
        ! all fields have the same row length
        u_adjust      = 0
        ! This is the number of v rows compared to p.
        v_adjust       = 1
        ! This is the number of extra points needed to allow interpolation onto
        ! p points safely.
        v_extra_ew    = 2
        v_extra_ns    = 1
        u_extra_ew    = 1
        u_extra_ns    = 1
      ELSE
        u_adjust      = -1
        v_adjust      = -1
        v_extra_ew    = 0
        v_extra_ns    = 1
        u_extra_ew    = 1
        u_extra_ns    = 0
      END IF

! --------------------------------------------------
! Determine global length of lbc field for one level
! --------------------------------------------------

      DO ifld=1,nfld_max
        DO ihalo=1,nhalo_max

          lbc_global_lenrima(ifld,ihalo) = 0

          DO iside=1,4

            IF ( (iside == PNorth) .OR. (iside == PSouth) ) THEN

! ND:
!     Full row length for p and v grid
!     One point less for u grid
! EG:
!     P, U and V will all have same row length
!     (equal to the intf_row_length value from the
!     namelist).

              IF (ifld == fld_type_u) THEN
                lbc_row_len = intf_row_length(jintf) + u_adjust
              ELSE
                lbc_row_len = intf_row_length(jintf)
              END IF
              lbc_row_len = lbc_row_len + 2 * intf_halosize(1,ihalo)

!     No of rows is rimwidth plus halo size

              lbc_rows   = intfwidtha(jintf) + intf_halosize(2,ihalo)

            ELSE ! East or West Boundaries

!     Row length is rimwidth plus halo size

              lbc_row_len = intfwidtha(jintf) + intf_halosize(1,ihalo)

! ND:
!     All rows used for p and u grid
!     One row less for v grid
! EG:
!     All rows used for p and u grid
!     One row MORE for v grid

              IF (ifld == fld_type_v) THEN
                lbc_rows = intf_p_rows(jintf) + v_adjust
              ELSE
                lbc_rows = intf_p_rows (jintf)
              END IF
              lbc_rows   = lbc_rows - 2 * intfwidtha(jintf)

            END IF

!     Accumulate grid size

            IF (intfwidtha(jintf) > 0 ) THEN
              lbc_global_lenrima(ifld,ihalo) =                          &
     &        lbc_global_lenrima(ifld,ihalo) + lbc_row_len * lbc_rows
            END IF

          ENDDO  ! iside
        ENDDO    ! ihalo
      ENDDO      ! ifld

!  ------------------------------------------------------------
!  Calculate the sizes of the fields for interpolation
!  u and v fields differ from the actual sizes.
!  An extra P col or row is required to interpolate u/v to
!  the p grid before it is interpolated to the u or v grid.
!  ------------------------------------------------------------

      DO ifld=1,nfld_max
        DO ihalo=1,nhalo_max

          lbc_interp_lenrima(ifld,ihalo) = 0

          DO iside=1,4

            IF ( (iside == PNorth) .OR. (iside == PSouth) ) THEN

! ND:
!     P grid is suitable for interpolation to U points in North/South
!     section without modification. Extra row needed at bottom of
!     North section and top of South section for interpolation from P to V.
! EG:
!     Need extra row at top of North section and bottom of South section
!     for interpolation from P to V grid. Need extra column of P points at
!     left of North/South sections for interpolation to U grid.
             
             IF (ifld == fld_type_u .OR. ifld == fld_type_v) THEN
              lbc_row_len = intf_row_length(jintf) + u_extra_ns
             ELSE
              lbc_row_len = intf_row_length(jintf)
             END IF
             
              lbc_row_len = lbc_row_len + 2 * intf_halosize(1,ihalo)

              lbc_rows   = intfwidtha(jintf) + intf_halosize(2,ihalo)
              
              IF (ifld == fld_type_u .OR. ifld == fld_type_v) THEN
                lbc_rows = lbc_rows + v_extra_ns
              END IF

            ELSE ! East or West Boundaries

! ND:
!     P grid is suitable for interpolation to V points in East/West
!     section without modification. Extra column needed at left/right
!     for interpolation from P to U grid.
! EG:
!     Need extra rows at top and bottom of East/West section to interpolate
!     from P to V grid. Need extra column at left of East/West section
!     to interpolate from P to U grid.

              lbc_row_len = intfwidtha(jintf) + intf_halosize(1,ihalo)
              
              IF (ifld == fld_type_u .OR. ifld == fld_type_v) THEN
                lbc_row_len = lbc_row_len + u_extra_ew
              END IF

              IF (ifld == fld_type_u .OR. ifld == fld_type_v) THEN
                lbc_rows = intf_p_rows(jintf) + v_extra_ew
              ELSE
                lbc_rows = intf_p_rows(jintf)
              END IF
              lbc_rows = lbc_rows - 2 * intfwidtha(jintf)

            END IF

            IF (intfwidtha(jintf) > 0 ) THEN
              lbc_interp_lenrima(ifld,ihalo) =                          &
              lbc_interp_lenrima(ifld,ihalo) + lbc_row_len * lbc_rows
            END IF

          ENDDO  ! iside
        ENDDO    ! ihalo
      ENDDO      ! ifld

      IF (lhook) CALL dr_hook('LBC_GRID_SIZES',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE LBC_Grid_Sizes
