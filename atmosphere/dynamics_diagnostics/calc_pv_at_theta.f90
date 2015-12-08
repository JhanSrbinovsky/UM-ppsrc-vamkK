! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Calculate PV at all theta points.

      SUBROUTINE Calc_PV_at_theta ( u, v, theta, rho,                   &
                                                                  ! in
     &                              r_theta_levels, r_rho_levels,       &
                                                                  ! in
     &                              r_at_u, r_at_v,                     &
                                                                  ! in
     &                              sec_v_latitude,                     &
                                                                  ! in
     &                              tan_v_latitude,                     &
                                                                  ! in
     &                              sec_theta_latitude,                 &
                                                                  ! in
     &                              f3_at_v,                            &
                                                                  ! in
     &                              delta_lambda, delta_phi,            &
                                                                  ! in
     &                              Model_domain,                       &
                                                                  ! in
     &                              pv_at_theta )                 ! out

! Description:
!
!   Calculate PV at all theta points.
!
! Method:
!
!   1. Call Calc_PV to obtain PV midway E-W between v points on rho
!      levels.
!   2. Add haloes to resulting field.
!   3. Interpolate horizontally and vertically to theta points.
!      (PV at top theta level is set to PV at top rho level.)
!   4. Reset polar rows to their mean values.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Diagnostics
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!
! Declarations:

      USE dynamics_grid_mod, ONLY: l_vatpoles

      USE atm_fields_bounds_mod, ONLY:                                  &
          udims, vdims, udims_s, vdims_s, tdims_s, pdims_s,             &
          udims_l, vdims_l

      USE global_2d_sums_mod, ONLY : global_2d_sums
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParVars
      USE domain_params
      IMPLICIT NONE

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

! Subroutine arguments:

      REAL, INTENT(IN) ::                                               &
      
        u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,  &
          udims_s%k_start:udims_s%k_end),                               &
      
        v(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,  &
          vdims_s%k_start:vdims_s%k_end),                               &
      
        theta(tdims_s%i_start:tdims_s%i_end,                            &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end),                          &
      
        rho(pdims_s%i_start:pdims_s%i_end,                              &
             pdims_s%j_start:pdims_s%j_end,                             &
             pdims_s%k_start:pdims_s%k_end),                            &
      
        r_theta_levels     ( 1 - halo_i : row_length + halo_i,          &
                             1 - halo_j : rows       + halo_j,          &
                             0          : model_levels ),               &
      
        r_rho_levels       ( 1 - halo_i : row_length + halo_i,          &
                             1 - halo_j : rows       + halo_j,          &
                             1          : model_levels ),               &
      
        r_at_u (udims_l%i_start:udims_l%i_end,                          &
                udims_l%j_start:udims_l%j_end,                          &
                udims_l%k_start:udims_l%k_end),                         &
      
        r_at_v (vdims_l%i_start:vdims_l%i_end,                          &
                vdims_l%j_start:vdims_l%j_end,                          &
                vdims_l%k_start:vdims_l%k_end),                         &
      
        sec_v_latitude (vdims_s%i_start:vdims_s%i_end,                  &
                        vdims_s%j_start:vdims_s%j_end),                 &
      
        tan_v_latitude(vdims%i_start:vdims%i_end,                       &
                       vdims%j_start:vdims%j_end),                      &
      
        sec_theta_latitude ( 1 - offx   : row_length + offx,            &
                             1 - offy   : rows       + offy ),          &
      
        f3_at_v (vdims_s%i_start:vdims_s%i_end,                         &
                 vdims_s%j_start:vdims_s%j_end),                        &
      
        delta_lambda,                                                   &
        delta_phi

      INTEGER, INTENT(IN) ::                                            &
      
     &  Model_domain

      REAL, INTENT(OUT) ::                                              &
      
     &  pv_at_theta (row_length, rows, model_levels)

! Local variables:

      INTEGER :: i, j, k
      INTEGER :: ICode

      REAL :: pv (udims%i_start:udims%i_end,                            &
                  vdims%j_start:vdims%j_end, model_levels)

      REAL :: pv_plus_haloes ( udims_s%i_start:udims_s%i_end,           &
                               vdims_s%j_start:vdims_s%j_end,           &
     &                         1          : model_levels )

      REAL :: polar_sums (model_levels)
      REAL :: polar_means(model_levels)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! External subroutines called:

      EXTERNAL Calc_PV, Swap_Bounds

!- End of header ------------------------------------------------------

!----------------------------------------------------------------------
! [1]: Calculate PV midway E-W between v points on rho levels.
!----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('CALC_PV_AT_THETA',zhook_in,zhook_handle)
! DEPENDS ON: calc_pv
      CALL Calc_PV  ( u, v, theta, rho,                                 &
                                                      ! in
     &                r_theta_levels, r_rho_levels,                     &
                                                      ! in
     &                r_at_u, r_at_v,                                   &
                                                      ! in
     &                sec_v_latitude, tan_v_latitude,                   &
                                                      ! in
     &                sec_theta_latitude, f3_at_v,                      &
                                                      ! in
     &                delta_lambda, delta_phi,                          &
                                                      ! in
     &                row_length, rows, n_rows,                         &
                                                      ! in
     &                model_levels,                                     &
                                                      ! in
     &                offx, offy, halo_i, halo_j,                       &
                                                      ! in
     &                at_extremity,                                     &
                                                      ! in
     &                pv )                            ! out

!----------------------------------------------------------------------
! [2]: Add haloes.
!----------------------------------------------------------------------

      pv_plus_haloes(:,:,:) = 0.0

      DO k = 1, model_levels
        DO j = vdims%j_start,vdims%j_end
          DO i = udims%i_start,udims%i_end
            pv_plus_haloes(i,j,k) = pv(i,j,k)
          END DO
        END DO
      END DO

! DEPENDS ON: swap_bounds
      CALL Swap_Bounds ( pv_plus_haloes,                                &
                                                             ! inout
     &                   row_length, n_rows, model_levels,              &
                                                             ! in
     &                   offx, offy, fld_type_v, .FALSE. )   ! in

!----------------------------------------------------------------------
! [3]: Interpolate to theta points.
!----------------------------------------------------------------------

      DO k = 1, model_levels - 1
        DO j = 1, rows
          DO i = 1, row_length
            pv_at_theta(i,j,k) = 0.25                                   &
     &                           * ( ( pv_plus_haloes(i,  j,  k)        &
     &                               + pv_plus_haloes(i,  j-1,k)        &
     &                               + pv_plus_haloes(i-1,j,  k)        &
     &                               + pv_plus_haloes(i-1,j-1,k)        &
     &                               )                                  &
     &                             * ( r_rho_levels  (i,j,k+1)          &
     &                               - r_theta_levels(i,j,k)            &
     &                               )                                  &
     &                             + ( pv_plus_haloes(i,  j,  k+1)      &
     &                               + pv_plus_haloes(i,  j-1,k+1)      &
     &                               + pv_plus_haloes(i-1,j,  k+1)      &
     &                               + pv_plus_haloes(i-1,j-1,k+1)      &
     &                               )                                  &
     &                             * ( r_theta_levels(i,j,k)            &
     &                               - r_rho_levels  (i,j,k)            &
     &                             ) )                                  &
     &                           / ( r_rho_levels(i,j,k+1)              &
     &                             - r_rho_levels(i,j,k)                &
     &                             )
          END DO
        END DO
      END DO

      ! Set PV at top theta level equal to PV at top rho level.
      k = model_levels
      DO j = 1, rows
        DO i = 1, row_length
          pv_at_theta(i,j,k) = 0.25                                     &
     &                       * ( pv_plus_haloes(i,  j,  k)              &
     &                         + pv_plus_haloes(i,  j-1,k)              &
     &                         + pv_plus_haloes(i-1,j,  k)              &
     &                         + pv_plus_haloes(i-1,j-1,k)              &
     &                         )
        END DO
      END DO

!----------------------------------------------------------------------
! [4]: Reset polar rows to their mean values.
!----------------------------------------------------------------------
      IF (.NOT. l_vatpoles) THEN
      IF (Model_domain == mt_global) THEN

        IF (at_extremity(PNorth)) THEN

          CALL global_2d_sums(pv_at_theta(:,rows:rows,:), row_length,   &
                              1, 0, 0, model_levels, polar_sums,        &
                              gc_proc_row_group)

           polar_means(:) = polar_sums(:) / REAL(global_row_length)

           DO k = 1, model_levels
             pv_at_theta(:,rows,k) = polar_means(k)
           END DO

        END IF   ! PNorth

        IF (at_extremity(PSouth)) THEN

          CALL global_2d_sums(pv_at_theta(:,1:1,:), row_length,         &
                              1, 0, 0, model_levels, polar_sums,        &
                              gc_proc_row_group)

           polar_means(:) = polar_sums(:) / REAL(global_row_length)

           DO k = 1, model_levels
             pv_at_theta(:,1,k) = polar_means(k)
           END DO

        END IF

      END IF ! (Model_domain == mt_global)
      END IF  ! vatpoles

      IF (lhook) CALL dr_hook('CALC_PV_AT_THETA',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Calc_PV_at_theta
