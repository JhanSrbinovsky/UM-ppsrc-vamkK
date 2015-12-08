! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
!  Description:
!    To evaluate concentrations from rate coefficients, J values,
!    and dry and wet deposition rates.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      SUBROUTINE UKCA_DERIV_AERO(nr_therm, nr_phot, n_be_calls,         &
                                 rc, wdep, ddep,                        &
                                 dj, k_dms, h2o, m, o2,                 &
                                 dts, y,                                &
                                 delta_SO2_wetox_H2O2,                  &
                                 delta_SO2_wetox_O3,                    &
                                 delta_SO2_dryox_OH)

      USE UKCA_CONSTANTS,   ONLY: avogadro
      USE ukca_option_mod, ONLY: jpspec, jppj, nit
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE UM_ParVars
      USE Control_Max_Sizes
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

      INTEGER, INTENT(IN) :: nr_therm      ! No of thermal reactions
      INTEGER, INTENT(IN) :: nr_phot       ! No of photolytic reactions
      INTEGER, INTENT(IN) :: n_be_calls    ! No diagnostics

      REAL, INTENT(IN) :: rc(theta_field_size,nr_therm)   ! Rate coeffs
      REAL, INTENT(IN) :: wdep(theta_field_size,jpspec)   ! Wet dep rates
      REAL, INTENT(IN) :: ddep(theta_field_size,jpspec)   ! Dry dep rates
      REAL, INTENT(IN) :: dj(theta_field_size,jppj)       ! Phot rates
      REAL, INTENT(IN) :: k_dms(theta_field_size,5)       ! DMS rate coeffs
      REAL, INTENT(IN) :: h2o(theta_field_size)           ! H2O
      REAL, INTENT(IN) :: m(theta_field_size)             ! Air density
      REAL, INTENT(IN) :: o2(theta_field_size)            ! Oxygen
      REAL, INTENT(IN) :: dts                   ! timestep

! Species concentrations
      REAL, INTENT(INOUT) :: y(theta_field_size,jpspec)

! Fluxes of SO2 by oxidation pathway
      REAL, INTENT(OUT) :: delta_SO2_wetox_H2O2(theta_field_size)
      REAL, INTENT(OUT) :: delta_SO2_wetox_O3(theta_field_size)
      REAL, INTENT(OUT) :: delta_SO2_dryox_OH(theta_field_size)

! Local variables
      REAL :: p1(theta_field_size)        ! Production and loss terms
      REAL :: p2(theta_field_size)
      REAL :: p(theta_field_size)
      REAL :: l(theta_field_size)
      REAL :: l1(theta_field_size)
      REAL :: l2(theta_field_size)
      REAL :: l3(theta_field_size)
      REAL :: r1(theta_field_size)
      REAL :: r2(theta_field_size)
      REAL :: f_so2(theta_field_size)     ! fraction of dms to SO2
      REAL :: f_so4(theta_field_size)     ! fraction of dms to SO4
      REAL :: f_msa(theta_field_size)     ! fraction of dms to MSA
      REAL :: yp(theta_field_size,jpspec) ! initial concentrations

      INTEGER :: i,j,n              ! loop counters

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('UKCA_DERIV_AERO',zhook_in,zhook_handle)

! Initialise sulphur oxidation fluxes
       delta_SO2_wetox_H2O2(:) = 0.0
       delta_SO2_wetox_O3(:) = 0.0
       delta_SO2_dryox_OH(:) = 0.0

      DO n = 1, n_be_calls

        DO j = 1,jpspec
          DO i = 1, theta_field_size
            yp(i,j) = y(i,j)
          END DO
        END DO

!      Set constant fields

       YP(:,12) = 5.0e-7*M(:)     ! H2
       YP(:,16) = 350.0e-6*M(:)   ! CO2
       YP(:,19) = H2O(:)          ! H2O
       YP(:,22) = O2(:)           ! O2
       YP(:,23) = M(:) - O2(:)    ! N2

       Y (:,12) = YP(:,12)        ! H2
       Y (:,16) = YP(:,16)        ! CO2
       Y (:,19) = YP(:,19)        ! H2O
       Y (:,22) = YP(:,22)        ! O2
       Y (:,23) = YP(:,23)        ! N2


! Iteration start
       DO I=1,NIT

!
!   This section written automatically by MECH9GEN from the file ukca_tropch3
!   with 129 equations and 46 species.
!
!          O(3P)        Y( 1)
      P(:) =                                                            &
      +(DJ(:,15) *Y(:,3 ))                                              &
      +(DJ(:,12) *Y(:,5 ))        +(DJ(:,13) *Y(:,22)*2.00)             &
      +(RC(:,76) *Y(:,10)*Y(:,10))+(DJ(:,10) *Y(:,6 ))                  &
      +(RC(:,49) *Y(:,2 )*Y(:,23))+(RC(:,50) *Y(:,2 )*Y(:,22))
      L(:) = DDEP(:,1) + WDEP(:,1)                                      &
      +(RC(:,44) *Y(:,3 ))+(RC(:,83) *Y(:,6 ))+(RC(:,96) *Y(:,22))
      Y(:, 1) = P(:)/L(:)
!
!          O(1D)        Y( 2)
      P(:) =                                                            &
      +(DJ(:,14) *Y(:,3 ))
      L(:) = DDEP(:,2) + WDEP(:,2)                                      &
      +(RC(:,48) *Y(:,19))+(RC(:,49) *Y(:,23))+(RC(:,50) *Y(:,22))      &
      +(RC(:,45) *Y(:,14))+(RC(:,46) *Y(:,14))+(RC(:,47) *Y(:,14))
      Y(:, 2) = P(:)/L(:)
!
!          OH           Y(10)
      P(:) =                                                            &
      +(DJ(:,1)  *Y(:,34))        +(DJ(:,1)  *Y(:,39))                  &
      +(DJ(:,17) *Y(:,21))        +(DJ(:,1)  *Y(:,33))                  &
      +(DJ(:,6)  *Y(:,9 ))        +(DJ(:,1)  *Y(:,20))                  &
      +(DJ(:,1)  *Y(:,26))        +(DJ(:,2)  *Y(:,13)*2.00)             &
      +(RC(:,82) *Y(:,10)*Y(:,34))+(RC(:,85) *Y(:,45)*Y(:,44))          &
      +(RC(:,72) *Y(:,10)*Y(:,39))+(RC(:,80) *Y(:,10)*Y(:,33))          &
      +(RC(:,57) *Y(:,10)*Y(:,26))+(RC(:,66) *Y(:,10)*Y(:,20))          &
      +(RC(:,45) *Y(:,2 )*Y(:,14))+(RC(:,48) *Y(:,2 )*Y(:,19)*2.00)     &
      +(RC(:,3)  *Y(:,11)*Y(:,3 ))+(RC(:,10) *Y(:,11)*Y(:,28))          &
      +(RC(:,1)  *Y(:,11)*Y(:,4 ))+(RC(:,2)  *Y(:,11)*Y(:,5 ))
      L(:) = DDEP(:,10) + WDEP(:,10)                                    &
      +(RC(:,99) *Y(:,10))                                              &
      +(RC(:,97) *Y(:,4 ))+(RC(:,98) *Y(:,6 ))+(RC(:,99) *Y(:,10))      &
      +(RC(:,80) *Y(:,33))+(RC(:,81) *Y(:,34))+(RC(:,82) *Y(:,34))      &
      +(RC(:,77) *Y(:,29))+(RC(:,78) *Y(:,40))+(RC(:,79) *Y(:,33))      &
      +(RC(:,75) *Y(:,3 ))+(RC(:,76) *Y(:,10))+(RC(:,76) *Y(:,10))      &
      +(RC(:,72) *Y(:,39))+(RC(:,73) *Y(:,27))+(RC(:,74) *Y(:,5 ))      &
      +(RC(:,69) *Y(:,37))+(RC(:,70) *Y(:,37))+(RC(:,71) *Y(:,39))      &
      +(RC(:,66) *Y(:,20))+(RC(:,67) *Y(:,20))+(RC(:,68) *Y(:,41))      &
      +(RC(:,63) *Y(:,8 ))+(RC(:,64) *Y(:,9 ))+(RC(:,65) *Y(:,21))      &
      +(RC(:,60) *Y(:,13))+(RC(:,61) *Y(:,17))+(RC(:,62) *Y(:,11))      &
      +(RC(:,57) *Y(:,26))+(RC(:,58) *Y(:,26))+(RC(:,59) *Y(:,12))      &
      +(RC(:,54) *Y(:,30))+(RC(:,55) *Y(:,15))+(RC(:,56) *Y(:,35))      &
      +(RC(:,51) *Y(:,14))+(RC(:,52) *Y(:,24))+(RC(:,53) *Y(:,30))      &
      +(RC(:,103)*Y(:,48))+(RC(:,111)*Y(:,51))+(RC(:,112)*Y(:,52))      &
      +(RC(:,105)*Y(:,47))+(RC(:,106)*Y(:,47))
      Y(:,10) = P(:)/L(:)
!
!          O3           Y( 3)
      P(:) =                                                            &
      +(RC(:,96) *Y(:,1 )*Y(:,22))                                      &
      +(RC(:,9)  *Y(:,11)*Y(:,28))+(RC(:,14) *Y(:,11)*Y(:,36))
      L(:) = DDEP(:,3) + WDEP(:,3)                                      &
      +(DJ(:,15) )                                                      &
      +(RC(:,44) *Y(:,1 ))+(RC(:,75) *Y(:,10))+(DJ(:,14) )              &
      +(RC(:,3)  *Y(:,11))+(RC(:,37) *Y(:,4 ))+(RC(:,38) *Y(:,6 ))      &
      +(RC(:,115)*Y(:,48))+(RC(:,113)*Y(:,52))
      Y(:, 3) = (YP(:, 3)+DTS*P(:))/(1.0+DTS*L(:))
!
!          NO           Y( 4)
      P(:) =                                                            &
      +(DJ(:,11) *Y(:,5 ))        +(DJ(:,17) *Y(:,21))                  &
      +(RC(:,83) *Y(:,1 )*Y(:,6 ))+(DJ(:,10) *Y(:,6 ))                  &
      -(RC(:,111)*Y(:,51)*Y(:,10))     ! from NH3+OH
      L(:) = DDEP(:,4) + WDEP(:,4)                                      &
      +(RC(:,36) *Y(:,5 ))+(RC(:,37) *Y(:,3 ))+(RC(:,97) *Y(:,10))      &
      +(RC(:,30) *Y(:,32))+(RC(:,32) *Y(:,36))+(RC(:,34) *Y(:,38))      &
      +(RC(:,23) *Y(:,25))+(RC(:,26) *Y(:,28))+(RC(:,28) *Y(:,31))      &
      +(RC(:,1)  *Y(:,11))+(RC(:,16) *Y(:,18))+(RC(:,17) *Y(:,18))
      Y(:, 4) = (YP(:, 4)+DTS*P(:))/(1.0+DTS*L(:))
!
!          NO2          Y( 6)
      P(:) =                                                            &
      +(DJ(:,20) *Y(:,41))                                              &
      +(DJ(:,16) *Y(:,29))        +(DJ(:,16) *Y(:,40))                  &
      +(DJ(:,9)  *Y(:,7 ))        +(DJ(:,12) *Y(:,5 ))                  &
      +(DJ(:,5)  *Y(:,8 ))        +(DJ(:,6)  *Y(:,9 ))                  &
      +(RC(:,94) *Y(:,7 ))        +(RC(:,101)*Y(:,40))                  &
      +(RC(:,91) *Y(:,8 ))        +(RC(:,93) *Y(:,29))                  &
      +(RC(:,77) *Y(:,10)*Y(:,29))+(RC(:,78) *Y(:,10)*Y(:,40))          &
      +(RC(:,68) *Y(:,10)*Y(:,41))+(RC(:,74) *Y(:,10)*Y(:,5 ))          &
      +(RC(:,63) *Y(:,10)*Y(:,8 ))+(RC(:,65) *Y(:,10)*Y(:,21))          &
      +(RC(:,36) *Y(:,4 )*Y(:,5 )*2.00)+(RC(:,37) *Y(:,4 )*Y(:,3 ))     &
      +(RC(:,34) *Y(:,38)*Y(:,4 ))+(RC(:,35) *Y(:,38)*Y(:,5 ))          &
      +(RC(:,32) *Y(:,36)*Y(:,4 ))+(RC(:,33) *Y(:,36)*Y(:,5 ))          &
      +(RC(:,30) *Y(:,32)*Y(:,4 ))+(RC(:,31) *Y(:,32)*Y(:,5 ))          &
      +(RC(:,28) *Y(:,31)*Y(:,4 ))+(RC(:,29) *Y(:,31)*Y(:,5 ))          &
      +(RC(:,26) *Y(:,28)*Y(:,4 ))+(RC(:,27) *Y(:,28)*Y(:,5 ))          &
      +(RC(:,23) *Y(:,25)*Y(:,4 ))+(RC(:,24) *Y(:,25)*Y(:,5 ))          &
      +(RC(:,16) *Y(:,18)*Y(:,4 ))+(RC(:,18) *Y(:,18)*Y(:,5 ))          &
      +(RC(:,1)  *Y(:,11)*Y(:,4 ))+(RC(:,2)  *Y(:,11)*Y(:,5 ))
      L(:) = DDEP(:,6) + WDEP(:,6)                                      &
      +(RC(:,100)*Y(:,36))+(DJ(:,10) )                                  &
      +(RC(:,92) *Y(:,28))+(RC(:,95) *Y(:,5 ))+(RC(:,98) *Y(:,10))      &
      +(RC(:,38) *Y(:,3 ))+(RC(:,83) *Y(:,1 ))+(RC(:,90) *Y(:,11))
      Y(:, 6) = (YP(:, 6)+DTS*P(:))/(1.0+DTS*L(:))
!
!          NO3          Y( 5)
!      P(:) =
!     &+(RC(:,94) *Y(:,7 ))        +(DJ(:,9)  *Y(:,7 ))
!     &+(RC(:,38) *Y(:,6 )*Y(:,3 ))+(RC(:,64) *Y(:,10)*Y(:,9 ))
!      L(:) = DDEP(:,5) + WDEP(:,5)
!     &+(DJ(:,11) )        +(DJ(:,12) )
!     &+(RC(:,42) *Y(:,37))+(RC(:,74) *Y(:,10))+(RC(:,95) *Y(:,6 ))
!     &+(RC(:,39) *Y(:,17))+(RC(:,40) *Y(:,27))+(RC(:,41) *Y(:,35))
!     &+(RC(:,33) *Y(:,36))+(RC(:,35) *Y(:,38))+(RC(:,36) *Y(:,4 ))
!     &+(RC(:,27) *Y(:,28))+(RC(:,29) *Y(:,31))+(RC(:,31) *Y(:,32))
!     &+(RC(:,2)  *Y(:,11))+(RC(:,18) *Y(:,18))+(RC(:,24) *Y(:,25))
!     &+(RC(:,107) *Y(:,47))
!      Y(:, 5) = (YP(:, 5)+DTS*P(:))/(1.0+DTS*L(:))
!
!          N2O5         Y( 7)
!      P(:) =
!     &+(RC(:,95) *Y(:,6 )*Y(:,5 ))
!      L(:) = DDEP(:,7) + WDEP(:,7)
!     &+(RC(:,43) *Y(:,19))+(RC(:,94) )        +(DJ(:,9)  )
!      Y(:, 7) = (YP(:, 7)+DTS*P(:))/(1.0+DTS*L(:))
!
!        NO3/N2O5  Y(:, 5)/Y(:, 7)
       P1(:) =                                                          &
      +(RC(:,38) *Y(:,6 )*Y(:,3 ))+(RC(:,64) *Y(:,10)*Y(:,9 ))
       L(:) = DDEP(:,5) + WDEP(:,5)                                     &
      +(DJ(:,11) )        +(DJ(:,12) )                                  &
      +(RC(:,42) *Y(:,37))+(RC(:,74) *Y(:,10))+(RC(:,95) *Y(:,6 ))      &
      +(RC(:,39) *Y(:,17))+(RC(:,40) *Y(:,27))+(RC(:,41) *Y(:,35))      &
      +(RC(:,33) *Y(:,36))+(RC(:,35) *Y(:,38))+(RC(:,36) *Y(:,4 ))      &
      +(RC(:,27) *Y(:,28))+(RC(:,29) *Y(:,31))+(RC(:,31) *Y(:,32))      &
      +(RC(:,2)  *Y(:,11))+(RC(:,18) *Y(:,18))+(RC(:,24) *Y(:,25))      &
      +(RC(:,107)*Y(:,47))+(RC(:,114)*Y(:,52))
      P2(:) = 0.0
      R1(:) = RC(:,94) + DJ(:,9)
      R2(:) = RC(:,95) *Y(:,6 )
      L1(:) = DDEP(:,7) + WDEP(:,7)                                     &
      +(RC(:,43) *Y(:,19))+(RC(:,94) )        +(DJ(:,9)  )
      L2(:) = 1.0+L(:)*DTS
      L3(:) = 1.0+L1(:)*DTS
      Y(:,5) = (L3(:)*(YP(:,5)+P1(:)*DTS)+R1(:)*DTS*YP(:,7))/           &
             ((L3(:)*L2(:))-R1(:)*R2(:)*DTS**2)
      Y(:,7) = (YP(:,7) + P2(:)*DTS + R2(:)*DTS*Y(:,5))/L3(:)
!
!          HO2NO2       Y( 8)
      P(:) =                                                            &
      +(RC(:,90) *Y(:,11)*Y(:,6 ))
      L(:) = DDEP(:,8) + WDEP(:,8)                                      &
      +(RC(:,63) *Y(:,10))+(RC(:,91) )        +(DJ(:,5)  )
      Y(:, 8) = (YP(:, 8)+DTS*P(:))/(1.0+DTS*L(:))
!
!          HONO2        Y( 9)
      P(:) =                                                            &
      +(RC(:,43) *Y(:,7 )*Y(:,19)*2.00)+(RC(:,98) *Y(:,10)*Y(:,6 ))     &
      +(RC(:,41) *Y(:,5 )*Y(:,35))+(RC(:,42) *Y(:,5 )*Y(:,37))          &
      +(RC(:,39) *Y(:,5 )*Y(:,17))+(RC(:,40) *Y(:,5 )*Y(:,27))          &
      +(RC(:,107)*Y(:,5 )*Y(:,47))
      L(:) = DDEP(:,9) + WDEP(:,9)                                      &
      +(RC(:,64) *Y(:,10))+(DJ(:,6)  )
      Y(:, 9) = (YP(:, 9)+DTS*P(:))/(1.0+DTS*L(:))
!
!          HO2          Y(11)
      P(:) =                                                            &
      +(DJ(:,1)  *Y(:,34))        +(DJ(:,20) *Y(:,41))                  &
      +(DJ(:,18) *Y(:,35))        +(DJ(:,1)  *Y(:,33))                  &
      +(DJ(:,7)  *Y(:,27))        +(DJ(:,1)  *Y(:,20))                  &
      +(DJ(:,3)  *Y(:,17)*2.00)        +(DJ(:,5)  *Y(:,8 ))             &
      +(RC(:,91) *Y(:,8 ))        +(DJ(:,1)  *Y(:,26))                  &
      +(RC(:,75) *Y(:,10)*Y(:,3 ))+(RC(:,84) *Y(:,46)*Y(:,44))          &
      +(RC(:,61) *Y(:,10)*Y(:,17))+(RC(:,74) *Y(:,10)*Y(:,5 ))          &
      +(RC(:,59) *Y(:,10)*Y(:,12))+(RC(:,60) *Y(:,10)*Y(:,13))          &
      +(RC(:,47) *Y(:,2 )*Y(:,14)*2.00)+(RC(:,55) *Y(:,10)*Y(:,15))     &
      +(RC(:,31) *Y(:,32)*Y(:,5 ))+(RC(:,39) *Y(:,5 )*Y(:,17))          &
      +(RC(:,29) *Y(:,31)*Y(:,5 ))+(RC(:,30) *Y(:,32)*Y(:,4 ))          &
      +(RC(:,25) *Y(:,25)*Y(:,28))+(RC(:,28) *Y(:,31)*Y(:,4 ))          &
      +(RC(:,23) *Y(:,25)*Y(:,4 ))+(RC(:,24) *Y(:,25)*Y(:,5 ))          &
      +(RC(:,20) *Y(:,18)*Y(:,18)*2.00)+(RC(:,21) *Y(:,18)*Y(:,28))     &
      +(RC(:,16) *Y(:,18)*Y(:,4 ))+(RC(:,18) *Y(:,18)*Y(:,5 ))
      L(:) = DDEP(:,11) + WDEP(:,11)                                    &
      +(RC(:,89) *Y(:,11))+(RC(:,90) *Y(:,6 ))                          &
      +(RC(:,15) *Y(:,38))+(RC(:,62) *Y(:,10))+(RC(:,89) *Y(:,11))      &
      +(RC(:,12) *Y(:,32))+(RC(:,13) *Y(:,36))+(RC(:,14) *Y(:,36))      &
      +(RC(:,9)  *Y(:,28))+(RC(:,10) *Y(:,28))+(RC(:,11) *Y(:,31))      &
      +(RC(:,6)  *Y(:,18))+(RC(:,7)  *Y(:,25))+(RC(:,8)  *Y(:,28))      &
      +(RC(:,4)  *Y(:,11))+(RC(:,4)  *Y(:,11))+(RC(:,5)  *Y(:,18))      &
      +(RC(:,1)  *Y(:,4 ))+(RC(:,2)  *Y(:,5 ))+(RC(:,3)  *Y(:,3 ))
      Y(:,11) = (YP(:,11)+DTS*P(:))/(1.0+DTS*L(:))
!
!          H2           Y(12)
!      P(:) =                                                            &
!     &+(RC(:,46) *Y(:,2 )*Y(:,14))+(DJ(:,4)  *Y(:,17))
!      L(:) = DDEP(:,12) + WDEP(:,12)
!     &+(RC(:,59) *Y(:,10))
!      Y(:,12) = (YP(:,12)+DTS*P(:))/(1.0+DTS*L(:))
       Y(:,12) = YP(:,12)
!
!          H2O2         Y(13)
      P(:) =                                                            &
      +(RC(:,99) *Y(:,10)*Y(:,10))                                      &
      +(RC(:,4)  *Y(:,11)*Y(:,11))+(RC(:,89) *Y(:,11)*Y(:,11))
      L(:) = DDEP(:,13) + WDEP(:,13)                                    &
      +(RC(:,60) *Y(:,10))+(DJ(:,2)  )                                  &
      +(RC(:,104)*Y(:,48))
      Y(:,13) = (YP(:,13)+DTS*P(:))/(1.0+DTS*L(:))
!
!          CH4          Y(14)
      P(:) =                                                            &
      +(DJ(:,8)  *Y(:,27))
      L(:) = DDEP(:,14) + WDEP(:,14)                                    &
      +(RC(:,51) *Y(:,10))                                              &
      +(RC(:,45) *Y(:,2 ))+(RC(:,46) *Y(:,2 ))+(RC(:,47) *Y(:,2 ))
      Y(:,14) = (YP(:,14)+DTS*P(:))/(1.0+DTS*L(:))
!
!          CO           Y(15)
      P(:) =                                                            &
      +(DJ(:,18) *Y(:,35))                                              &
      +(DJ(:,7)  *Y(:,27))        +(DJ(:,8)  *Y(:,27))                  &
      +(DJ(:,3)  *Y(:,17))        +(DJ(:,4)  *Y(:,17))                  &
      +(RC(:,39) *Y(:,5 )*Y(:,17))+(RC(:,61) *Y(:,10)*Y(:,17))
      L(:) = DDEP(:,15) + WDEP(:,15)                                    &
      +(RC(:,55) *Y(:,10))
      Y(:,15) = (YP(:,15)+DTS*P(:))/(1.0+DTS*L(:))
!
!          CO2          Y(16)
!      P(:) =                                                            &
!     &+(RC(:,32) *Y(:,36)*Y(:,4 ))+(RC(:,33) *Y(:,36)*Y(:,5 ))
!     &+(RC(:,26) *Y(:,28)*Y(:,4 ))+(RC(:,27) *Y(:,28)*Y(:,5 ))
!      L(:) = DDEP(:,16) + WDEP(:,16)
      Y(:,16) = YP(:,16)
!
!          HCHO         Y(17)
      P(:) =                                                            &
      +(DJ(:,20) *Y(:,41))                                              &
      +(DJ(:,1)  *Y(:,20))        +(DJ(:,1)  *Y(:,39))                  &
      +(RC(:,68) *Y(:,10)*Y(:,41))+(RC(:,77) *Y(:,10)*Y(:,29))          &
      +(RC(:,47) *Y(:,2 )*Y(:,14))+(RC(:,66) *Y(:,10)*Y(:,20))          &
      +(RC(:,35) *Y(:,38)*Y(:,5 ))+(RC(:,46) *Y(:,2 )*Y(:,14))          &
      +(RC(:,22) *Y(:,18)*Y(:,28))+(RC(:,34) *Y(:,38)*Y(:,4 ))          &
      +(RC(:,20) *Y(:,18)*Y(:,18)*2.00)+(RC(:,21) *Y(:,18)*Y(:,28))     &
      +(RC(:,18) *Y(:,18)*Y(:,5 ))+(RC(:,19) *Y(:,18)*Y(:,18))          &
      +(RC(:,6)  *Y(:,11)*Y(:,18))+(RC(:,16) *Y(:,18)*Y(:,4 ))
      L(:) = DDEP(:,17) + WDEP(:,17)                                    &
      +(DJ(:,4)  )                                                      &
      +(RC(:,39) *Y(:,5 ))+(RC(:,61) *Y(:,10))+(DJ(:,3)  )
      Y(:,17) = (YP(:,17)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MeOO         Y(18)
      P(:) =                                                            &
      +(DJ(:,7)  *Y(:,27))        +(DJ(:,19) *Y(:,37))                  &
      +(RC(:,51) *Y(:,10)*Y(:,14))+(RC(:,67) *Y(:,10)*Y(:,20))          &
      +(RC(:,27) *Y(:,28)*Y(:,5 ))+(RC(:,45) *Y(:,2 )*Y(:,14))          &
      +(RC(:,25) *Y(:,25)*Y(:,28))+(RC(:,26) *Y(:,28)*Y(:,4 ))          &
      +(RC(:,10) *Y(:,11)*Y(:,28))+(RC(:,21) *Y(:,18)*Y(:,28))
      L(:) = DDEP(:,18) + WDEP(:,18)                                    &
      +(RC(:,21) *Y(:,28))+(RC(:,22) *Y(:,28))                          &
      +(RC(:,19) *Y(:,18))+(RC(:,20) *Y(:,18))+(RC(:,20) *Y(:,18))      &
      +(RC(:,17) *Y(:,4 ))+(RC(:,18) *Y(:,5 ))+(RC(:,19) *Y(:,18))      &
      +(RC(:,5)  *Y(:,11))+(RC(:,6)  *Y(:,11))+(RC(:,16) *Y(:,4 ))
      Y(:,18) = (YP(:,18)+DTS*P(:))/(1.0+DTS*L(:))
!
!          H2O          Y(19)
!      P(:) =                                                            &
!     &+(RC(:,81) *Y(:,10)*Y(:,34))+(RC(:,86) *Y(:,43)*Y(:,19))          &
!     &+(RC(:,79) *Y(:,10)*Y(:,33))+(RC(:,80) *Y(:,10)*Y(:,33))          &
!     &+(RC(:,77) *Y(:,10)*Y(:,29))+(RC(:,78) *Y(:,10)*Y(:,40))          &
!     &+(RC(:,73) *Y(:,10)*Y(:,27))+(RC(:,76) *Y(:,10)*Y(:,10))          &
!     &+(RC(:,70) *Y(:,10)*Y(:,37))+(RC(:,71) *Y(:,10)*Y(:,39))          &
!     &+(RC(:,68) *Y(:,10)*Y(:,41))+(RC(:,69) *Y(:,10)*Y(:,37))          &
!     &+(RC(:,66) *Y(:,10)*Y(:,20))+(RC(:,67) *Y(:,10)*Y(:,20))          &
!     &+(RC(:,64) *Y(:,10)*Y(:,9 ))+(RC(:,65) *Y(:,10)*Y(:,21))          &
!     &+(RC(:,62) *Y(:,10)*Y(:,11))+(RC(:,63) *Y(:,10)*Y(:,8 ))          &
!     &+(RC(:,60) *Y(:,10)*Y(:,13))+(RC(:,61) *Y(:,10)*Y(:,17))          &
!     &+(RC(:,58) *Y(:,10)*Y(:,26))+(RC(:,59) *Y(:,10)*Y(:,12))          &
!     &+(RC(:,56) *Y(:,10)*Y(:,35))+(RC(:,57) *Y(:,10)*Y(:,26))          &
!     &+(RC(:,53) *Y(:,10)*Y(:,30))+(RC(:,54) *Y(:,10)*Y(:,30))          &
!     &+(RC(:,51) *Y(:,10)*Y(:,14))+(RC(:,52) *Y(:,10)*Y(:,24))
!      L(:) = DDEP(:,19) + WDEP(:,19)                                    &
!     &+(RC(:,43) *Y(:,7 ))+(RC(:,48) *Y(:,2 ))+(RC(:,86) *Y(:,43))
      Y(:,19) = YP(:,19)
!
!          MeOOH        Y(20)
      P(:) =                                                            &
      +(RC(:,5)  *Y(:,11)*Y(:,18))
      L(:) = DDEP(:,20) + WDEP(:,20)                                    &
      +(RC(:,66) *Y(:,10))+(RC(:,67) *Y(:,10))+(DJ(:,1)  )
      Y(:,20) = (YP(:,20)+DTS*P(:))/(1.0+DTS*L(:))
!
!          HONO         Y(21)
      P(:) =                                                            &
      +(RC(:,97) *Y(:,10)*Y(:,4 ))
      L(:) = DDEP(:,21) + WDEP(:,21)                                    &
      +(RC(:,65) *Y(:,10))+(DJ(:,17) )
      Y(:,21) = (YP(:,21)+DTS*P(:))/(1.0+DTS*L(:))
!
!          O2           Y(22)
!      P(:) =                                                            &
!     &+(DJ(:,15) *Y(:,44))
!     &+(DJ(:,15) *Y(:,3 ))        +(DJ(:,14) *Y(:,44))
!     &+(DJ(:,11) *Y(:,5 ))        +(DJ(:,14) *Y(:,3 ))
!     &+(RC(:,88) *Y(:,43)*Y(:,22))+(RC(:,89) *Y(:,11)*Y(:,11))
!     &+(RC(:,84) *Y(:,46)*Y(:,44))+(RC(:,85) *Y(:,45)*Y(:,44))
!     &+(RC(:,75) *Y(:,10)*Y(:,3 ))+(RC(:,83) *Y(:,1 )*Y(:,6 ))
!     &+(RC(:,44) *Y(:,1 )*Y(:,3 )*2.00)+(RC(:,50) *Y(:,2 )*Y(:,22))
!     &+(RC(:,3)  *Y(:,11)*Y(:,3 ))+(RC(:,13) *Y(:,11)*Y(:,36))
!      L(:) = DDEP(:,22) + WDEP(:,22)
!     &+(RC(:,102)*Y(:,42))+(DJ(:,13) )
!     &+(RC(:,50) *Y(:,2 ))+(RC(:,88) *Y(:,43))+(RC(:,96) *Y(:,1 ))
      Y(:,22) = YP(:,22)
!
!          N2           Y(23)
!      P(:) =                                                            &
!     &+(RC(:,49) *Y(:,2 )*Y(:,23))+(RC(:,87) *Y(:,43)*Y(:,23))
!      L(:) = DDEP(:,23) + WDEP(:,23)
!     &+(RC(:,49) *Y(:,2 ))+(RC(:,87) *Y(:,43))
      Y(:,23) = YP(:,23)
!
!          C2H6         Y(24)
      P(:) = 0.0
      L(:) = DDEP(:,24) + WDEP(:,24)                                    &
      +(RC(:,52) *Y(:,10))
      Y(:,24) = (YP(:,24)+DTS*P(:))/(1.0+DTS*L(:))
!
!          EtOO         Y(25)
      P(:) =                                                            &
      +(DJ(:,18) *Y(:,35))                                              &
      +(RC(:,52) *Y(:,10)*Y(:,24))+(RC(:,58) *Y(:,10)*Y(:,26))          &
      +(RC(:,32) *Y(:,36)*Y(:,4 ))+(RC(:,33) *Y(:,36)*Y(:,5 ))
      L(:) = DDEP(:,25) + WDEP(:,25)                                    &
      +(RC(:,25) *Y(:,28))                                              &
      +(RC(:,7)  *Y(:,11))+(RC(:,23) *Y(:,4 ))+(RC(:,24) *Y(:,5 ))
      Y(:,25) = (YP(:,25)+DTS*P(:))/(1.0+DTS*L(:))
!
!          EtOOH        Y(26)
      P(:) =                                                            &
      +(RC(:,7)  *Y(:,11)*Y(:,25))
      L(:) = DDEP(:,26) + WDEP(:,26)                                    &
      +(RC(:,57) *Y(:,10))+(RC(:,58) *Y(:,10))+(DJ(:,1)  )
      Y(:,26) = (YP(:,26)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MeCHO        Y(27)
      P(:) =                                                            &
      +(RC(:,78) *Y(:,10)*Y(:,40))+(DJ(:,1)  *Y(:,26))                  &
      +(RC(:,25) *Y(:,25)*Y(:,28))+(RC(:,57) *Y(:,10)*Y(:,26))          &
      +(RC(:,23) *Y(:,25)*Y(:,4 ))+(RC(:,24) *Y(:,25)*Y(:,5 ))
      L(:) = DDEP(:,27) + WDEP(:,27)                                    &
      +(DJ(:,8)  )                                                      &
      +(RC(:,40) *Y(:,5 ))+(RC(:,73) *Y(:,10))+(DJ(:,7)  )
      Y(:,27) = (YP(:,27)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MeCO3        Y(28)
      P(:) =                                                            &
      +(DJ(:,19) *Y(:,37))        +(DJ(:,1)  *Y(:,39))                  &
      +(RC(:,93) *Y(:,29))        +(DJ(:,16) *Y(:,29))                  &
      +(RC(:,40) *Y(:,5 )*Y(:,27))+(RC(:,73) *Y(:,10)*Y(:,27))          &
      +(RC(:,34) *Y(:,38)*Y(:,4 ))+(RC(:,35) *Y(:,38)*Y(:,5 ))
      L(:) = DDEP(:,28) + WDEP(:,28)                                    &
      +(RC(:,26) *Y(:,4 ))+(RC(:,27) *Y(:,5 ))+(RC(:,92) *Y(:,6 ))      &
      +(RC(:,21) *Y(:,18))+(RC(:,22) *Y(:,18))+(RC(:,25) *Y(:,25))      &
      +(RC(:,8)  *Y(:,11))+(RC(:,9)  *Y(:,11))+(RC(:,10) *Y(:,11))
      Y(:,28) = (YP(:,28)+DTS*P(:))/(1.0+DTS*L(:))
!
!          PAN          Y(29)
      P(:) =                                                            &
      +(RC(:,92) *Y(:,28)*Y(:,6 ))
      L(:) = DDEP(:,29) + WDEP(:,29)                                    &
      +(RC(:,77) *Y(:,10))+(RC(:,93) )        +(DJ(:,16) )
      Y(:,29) = (YP(:,29)+DTS*P(:))/(1.0+DTS*L(:))
!
!          C3H8         Y(30)
      P(:) = 0.0
      L(:) = DDEP(:,30) + WDEP(:,30)                                    &
      +(RC(:,53) *Y(:,10))+(RC(:,54) *Y(:,10))
      Y(:,30) = (YP(:,30)+DTS*P(:))/(1.0+DTS*L(:))
!
!          n-PrOO       Y(31)
      P(:) =                                                            &
      +(RC(:,53) *Y(:,10)*Y(:,30))+(RC(:,79) *Y(:,10)*Y(:,33))
      L(:) = DDEP(:,31) + WDEP(:,31)                                    &
      +(RC(:,11) *Y(:,11))+(RC(:,28) *Y(:,4 ))+(RC(:,29) *Y(:,5 ))
      Y(:,31) = (YP(:,31)+DTS*P(:))/(1.0+DTS*L(:))
!
!          i-PrOO       Y(32)
      P(:) =                                                            &
      +(RC(:,54) *Y(:,10)*Y(:,30))+(RC(:,81) *Y(:,10)*Y(:,34))
      L(:) = DDEP(:,32) + WDEP(:,32)                                    &
      +(RC(:,12) *Y(:,11))+(RC(:,30) *Y(:,4 ))+(RC(:,31) *Y(:,5 ))
      Y(:,32) = (YP(:,32)+DTS*P(:))/(1.0+DTS*L(:))
!
!          n-PrOOH      Y(33)
      P(:) =                                                            &
      +(RC(:,11) *Y(:,11)*Y(:,31))
      L(:) = DDEP(:,33) + WDEP(:,33)                                    &
      +(RC(:,79) *Y(:,10))+(RC(:,80) *Y(:,10))+(DJ(:,1)  )
      Y(:,33) = (YP(:,33)+DTS*P(:))/(1.0+DTS*L(:))
!
!          i-PrOOH      Y(34)
      P(:) =                                                            &
      +(RC(:,12) *Y(:,11)*Y(:,32))
      L(:) = DDEP(:,34) + WDEP(:,34)                                    &
      +(RC(:,81) *Y(:,10))+(RC(:,82) *Y(:,10))+(DJ(:,1)  )
      Y(:,34) = (YP(:,34)+DTS*P(:))/(1.0+DTS*L(:))
!
!          EtCHO        Y(35)
      P(:) =                                                            &
      +(RC(:,80) *Y(:,10)*Y(:,33))+(DJ(:,1)  *Y(:,33))                  &
      +(RC(:,28) *Y(:,31)*Y(:,4 ))+(RC(:,29) *Y(:,31)*Y(:,5 ))
      L(:) = DDEP(:,35) + WDEP(:,35)                                    &
      +(RC(:,41) *Y(:,5 ))+(RC(:,56) *Y(:,10))+(DJ(:,18) )
      Y(:,35) = (YP(:,35)+DTS*P(:))/(1.0+DTS*L(:))
!
!          EtCO3        Y(36)
      P(:) =                                                            &
      +(RC(:,101)*Y(:,40))        +(DJ(:,16) *Y(:,40))                  &
      +(RC(:,41) *Y(:,5 )*Y(:,35))+(RC(:,56) *Y(:,10)*Y(:,35))
      L(:) = DDEP(:,36) + WDEP(:,36)                                    &
      +(RC(:,33) *Y(:,5 ))+(RC(:,100)*Y(:,6 ))                          &
      +(RC(:,13) *Y(:,11))+(RC(:,14) *Y(:,11))+(RC(:,32) *Y(:,4 ))
      Y(:,36) = (YP(:,36)+DTS*P(:))/(1.0+DTS*L(:))
!
!          Me2CO        Y(37)
      P(:) =                                                            &
      +(RC(:,82) *Y(:,10)*Y(:,34))+(DJ(:,1)  *Y(:,34))                  &
      +(RC(:,30) *Y(:,32)*Y(:,4 ))+(RC(:,31) *Y(:,32)*Y(:,5 ))
      L(:) = DDEP(:,37) + WDEP(:,37)                                    &
      +(DJ(:,19) )                                                      &
      +(RC(:,42) *Y(:,5 ))+(RC(:,69) *Y(:,10))+(RC(:,70) *Y(:,10))
      Y(:,37) = (YP(:,37)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MeCOCH2O     Y(38)
      P(:) =                                                            &
      +(RC(:,70) *Y(:,10)*Y(:,37))+(RC(:,71) *Y(:,10)*Y(:,39))          &
      +(RC(:,42) *Y(:,5 )*Y(:,37))+(RC(:,69) *Y(:,10)*Y(:,37))
      L(:) = DDEP(:,38) + WDEP(:,38)                                    &
      +(RC(:,15) *Y(:,11))+(RC(:,34) *Y(:,4 ))+(RC(:,35) *Y(:,5 ))
      Y(:,38) = (YP(:,38)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MeCOCH2O     Y(39)
      P(:) =                                                            &
      +(RC(:,15) *Y(:,11)*Y(:,38))
      L(:) = DDEP(:,39) + WDEP(:,39)                                    &
      +(RC(:,71) *Y(:,10))+(RC(:,72) *Y(:,10))+(DJ(:,1)  )
      Y(:,39) = (YP(:,39)+DTS*P(:))/(1.0+DTS*L(:))
!
!          PPAN         Y(40)
      P(:) =                                                            &
      +(RC(:,100)*Y(:,36)*Y(:,6 ))
      L(:) = DDEP(:,40) + WDEP(:,40)                                    &
      +(RC(:,78) *Y(:,10))+(RC(:,101))        +(DJ(:,16) )
      Y(:,40) = (YP(:,40)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MeONO2       Y(41)
      P(:) =                                                            &
      +(RC(:,17) *Y(:,18)*Y(:,4 ))
      L(:) = DDEP(:,41) + WDEP(:,41)                                    &
      +(RC(:,68) *Y(:,10))+(DJ(:,20) )
      Y(:,41) = (YP(:,41)+DTS*P(:))/(1.0+DTS*L(:))
!
!          O(3P)S       Y(42)
      P(:) =                                                            &
      +(DJ(:,15) *Y(:,44))                                              &
      +(RC(:,87) *Y(:,43)*Y(:,23))+(RC(:,88) *Y(:,43)*Y(:,22))
      L(:) = DDEP(:,42) + WDEP(:,42)                                    &
      +(RC(:,102)*Y(:,22))
      Y(:,42) = (YP(:,42)+DTS*P(:))/(1.0+DTS*L(:))
!
!          O(1D)S       Y(43)
      P(:) =                                                            &
      +(DJ(:,14) *Y(:,44))
      L(:) = DDEP(:,43) + WDEP(:,43)                                    &
      +(RC(:,86) *Y(:,19))+(RC(:,87) *Y(:,23))+(RC(:,88) *Y(:,22))
      Y(:,43) = (YP(:,43)+DTS*P(:))/(1.0+DTS*L(:))
!
!          O3S          Y(44)
      P(:) =                                                            &
      +(RC(:,102)*Y(:,42)*Y(:,22))
      L(:) = DDEP(:,44) + WDEP(:,44)                                    &
      +(DJ(:,15) )                                                      &
      +(RC(:,84) *Y(:,46))+(RC(:,85) *Y(:,45))+(DJ(:,14) )
      Y(:,44) = (YP(:,44)+DTS*P(:))/(1.0+DTS*L(:))
!
!          OHS          Y(45)
      P(:) = 0.0
      L(:) = DDEP(:,45) + WDEP(:,45)                                    &
      +(RC(:,85) *Y(:,44))
      Y(:,45) = (YP(:,45)+DTS*P(:))/(1.0+DTS*L(:))
!
!          HO2S         Y(46)
      P(:) = 0.0
      L(:) = DDEP(:,46) + WDEP(:,46)                                    &
      +(RC(:,84) *Y(:,44))
      Y(:,46) = (YP(:,46)+DTS*P(:))/(1.0+DTS*L(:))
!
! Fraction of DMS oxidation to products
      f_so2(:) = k_dms(:,1)/                                            &
                 (k_dms(:,1)+k_dms(:,2)*y(:,3)+k_dms(:,3)*y(:,6))
      f_so4(:) = (1.0-f_so2)*k_dms(:,5)/(k_dms(:,5)+k_dms(:,4)*y(:,11))
      f_msa(:) = (1.0-f_so2)*k_dms(:,4)*y(:,11)/                        &
                 (k_dms(:,5)+k_dms(:,4)*y(:,11))
!
!          DMS          Y(47)
      P(:) = 0.0
      L(:) = DDEP(:,47) + WDEP(:,47)                                    &
      +(RC(:,105) *Y(:,10))+(RC(:,106) *Y(:,10))                        &
      +(RC(:,107) *Y(:,5))
      Y(:,47) = (YP(:,47)+DTS*P(:))/(1.0+DTS*L(:))
!
!          SO2          Y(48)
      P(:) = 0.0                                                        &
      +(RC(:,105)*Y(:,10)+RC(:,106)*Y(:,10)                             &
      + RC(:,107)*Y(:,5))*Y(:,47)*f_so2(:)
      L(:) = 0.0 + DDEP(:,48) + WDEP(:,48)                              &
      +(RC(:,103) * Y(:,10))+(RC(:,104) * Y(:,13))                      &
      +(RC(:,115) * Y(:,3))
      Y(:,48) = (YP(:,48)+DTS*P(:))/(1.0+DTS*L(:))
!
!          H2SO4        Y(49)
      P(:) = 0.0                                                        &
      +(RC(:,103) * Y(:,48) * Y(:,10))                                  &
      +(RC(:,105)*Y(:,10)+RC(:,106)*Y(:,10)                             &
      + RC(:,107)*Y(:,5))*Y(:,47)*f_so4(:)
      L(:) = 0.0
      Y(:,49) = (YP(:,49)+DTS*P(:))/(1.0+DTS*L(:))
!
!          MSA          Y(50)
      P(:) = 0.0                                                        &
      +(RC(:,105)*Y(:,10)+RC(:,106)*Y(:,10)                             &
      + RC(:,107)*Y(:,5))*Y(:,47)*f_msa(:)
      L(:) = 0.0
      Y(:,50) = (YP(:,50)+DTS*P(:))/(1.0+DTS*L(:))

!          NH3          Y(51)
      P(:) = 0.0
      L(:) = DDEP(:,51) + WDEP(:,51)                                    &
      +(RC(:,111) * Y(:,10))
      Y(:,51) = (YP(:,51)+DTS*P(:))/(1.0+DTS*L(:))

!          MONOTERP     Y(52)
      P(:) = 0.0
      L(:) = DDEP(:,52)                                                 &
      +(RC(:,112) * Y(:,10)) + (RC(:,113) * Y(:,3))                     &
      +(RC(:,114) * Y(:,5))
      Y(:,52) = (YP(:,52)+DTS*P(:))/(1.0+DTS*L(:))

!          SEC_ORG      Y(53)
      P(:) = 0.0                                                        &
      +(RC(:,112)*Y(:,52)*Y(:,10)*0.13 + RC(:,113)*Y(:,52)*Y(:,3)*0.13) &
      +(RC(:,114)*Y(:,52)*Y(:,5)*0.13)
      L(:) = DDEP(:,53) + WDEP(:,53)
      Y(:,53) = (YP(:,53)+DTS*P(:))/(1.0+DTS*L(:))

!      iteration loop stop
      ENDDO

! Calculate flux terms at end of iteration.

! Fluxes to aqueous sulphate (always required for MODE):
      delta_SO2_wetox_H2O2(:)= delta_SO2_wetox_H2O2(:) +                &
                               rc(:,104)*y(:,13)*y(:,48)*dts
      delta_SO2_wetox_O3(:)  = delta_SO2_wetox_O3(:) +                  &
                               rc(:,115)*y(:,3)*y(:,48)*dts
      delta_SO2_dryox_OH(:)  = delta_SO2_dryox_OH(:) +                  &
                               rc(:,103)*y(:,10)*y(:,48)*dts

      END DO  ! n_be_calls

      IF (lhook) CALL dr_hook('UKCA_DERIV_AERO',zhook_out,zhook_handle)
      RETURN


      END SUBROUTINE UKCA_DERIV_AERO
