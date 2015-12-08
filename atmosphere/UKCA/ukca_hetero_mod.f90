! *****************************COPYRIGHT*******************************
! 
! (c) [University of Cambridge] [2008]. All rights reserved. 
! This routine has been licensed to the Met Office for use and 
! distribution under the UKCA collaboration agreement, subject  
! to the terms and conditions set out therein. 
! [Met Office Ref SC138]  
! 
! *****************************COPYRIGHT*******************************
!
!  Description:
!    To give calculate heterogeneous rates for UKCA.
!    Contains the following routines:
!     UKCA_HETERO
!     UKCA_SOLIDPHASE
!     UKCA_CALCKPSC
!     UKCA_EQCOMP
!     UKCA_POSITION
!     UKCA_PSCPRES
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!          Originally used in SLIMCAT CTM.
!
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      MODULE ukca_hetero_mod


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Control_Max_Sizes
      IMPLICIT NONE

      CONTAINS

      SUBROUTINE ukca_hetero(n_points, stratflag)
! Description:
!
! Changed from version by Peter Breasicke to allow for dynamical limitation of
! reaction rates as functions of the abundances of the reaction partners.
! Also remove a bug with the ordering of reactions, and allow for changing
! of the order of reactions without recoding. Finally, introduce limits
! on pressure, temperature and latitude where heterogeneous chemistry is
! performed.
!
!
! Method: Heterogeneous reaction rates are specified as pseudo-bimolecular
!         reactions.
!
!
! Code Description:
! Language: FORTRAN 90 + common extensions.
!
! Declarations:
! These are of the form:-

      USE ASAD_MOD
      USE UKCA_D1_DEFS
      USE ukca_option_mod, ONLY: jpctr, jphk 
      USE ereport_mod, ONLY : ereport
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

! Subroutine interface
      INTEGER, INTENT(IN) :: n_points
! logical to indicate whether point is in stratosphere
      LOGICAL, INTENT(IN) :: stratflag(theta_field_size)

! Local variables
      REAL, PARAMETER :: limit = 1.0E-12

      LOGICAL, SAVE :: gpsa
      LOGICAL, SAVE :: gphocl
      LOGICAL, SAVE :: gppsc
      LOGICAL, SAVE :: gpsimp
      LOGICAL :: L_ukca_sulphur
      LOGICAL :: L_ukca_presaer

      INTEGER :: js
      INTEGER :: jh
! Tracer names
      INTEGER, SAVE :: ih2o=0
      INTEGER, SAVE :: ihno3=0
      INTEGER, SAVE :: ihcl=0
      INTEGER, SAVE :: iclono2=0
      INTEGER, SAVE :: ihocl=0
      INTEGER, SAVE :: in2o5=0
! reaction names
      INTEGER, SAVE :: n_clono2_hcl=0
      INTEGER, SAVE :: n_clono2_h2o=0
      INTEGER, SAVE :: n_n2o5_h2o=0
      INTEGER, SAVE :: n_n2o5_hcl=0
      INTEGER, SAVE :: n_hocl_hcl=0

      LOGICAL, SAVE :: first = .TRUE.

      REAL :: zp(theta_field_size)
      REAL :: zt(theta_field_size)
      REAL :: zhno3(theta_field_size)
      REAL :: zh2o(theta_field_size)
      REAL :: zhcl(theta_field_size)
      REAL :: zclono2(theta_field_size)
      REAL :: zn2o5(theta_field_size)
      REAL :: zhocl(theta_field_size)
      REAL :: psc1(theta_field_size)
      REAL :: psc2(theta_field_size)
      REAL :: psc3(theta_field_size)
      REAL :: psc4(theta_field_size)
      REAL :: psc5(theta_field_size)
      REAL :: hk(theta_field_size,5)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!

      IF (lhook) CALL dr_hook('UKCA_HETERO_MOD:UKCA_HETERO',zhook_in,zhook_handle)
      IF (first) THEN
        DO js = 1, jpctr
          SELECT CASE (advt(js))
            CASE ('H2O       ','H2OS      ')
              ih2o = js
            CASE ('HONO2     ')
              ihno3 = js
            CASE ('HCl       ')
              ihcl = js
            CASE ('ClONO2    ')
              iclono2 = js
            CASE ('N2O5      ')
              in2o5 = js
            CASE ('HOCl      ')
              ihocl = js
          END SELECT
        END DO

        DO jh = 1, jphk
          SELECT CASE (sph(jh,1))
            CASE ('H2O       ','H2OS      ')
              IF (sph(jh,2) == 'ClONO2    ') n_clono2_h2o = nhrkx(jh)
              IF (sph(jh,2) == 'N2O5      ') n_n2o5_h2o   = nhrkx(jh)
            CASE ('ClONO2    ')
              IF (sph(jh,2) == 'HCl       ') n_clono2_hcl = nhrkx(jh)
              IF (sph(jh,2) == 'H2O       ') n_clono2_h2o = nhrkx(jh)
              IF (sph(jh,2) == 'H2OS      ') n_clono2_h2o = nhrkx(jh)
            CASE ('HCl       ')
              IF (sph(jh,2) == 'ClONO2    ') n_clono2_hcl = nhrkx(jh)
              IF (sph(jh,2) == 'N2O5      ') n_n2o5_hcl   = nhrkx(jh)
              IF (sph(jh,2) == 'HOCl      ') n_hocl_hcl   = nhrkx(jh)
            CASE ('N2O5      ')
              IF (sph(jh,2) == 'H2O       ') n_n2o5_h2o   = nhrkx(jh)
              IF (sph(jh,2) == 'H2OS      ') n_n2o5_h2o   = nhrkx(jh)
              IF (sph(jh,2) == 'HCl       ') n_n2o5_hcl   = nhrkx(jh)
            CASE ('HOCl      ')
              IF (sph(jh,2) == 'HCl       ') n_hocl_hcl   = nhrkx(jh)
          END SELECT
        END DO

        first = .FALSE.
        l_ukca_sulphur=.FALSE.
        L_ukca_presaer=.TRUE.

! activate sulphur chemistry according to UKCA d1_defs
        gpsa = L_ukca_sulphur .OR. L_ukca_presaer
! DO include HCl + HOCl reaction on SA aerosols
        gphocl = .TRUE.
! Use heterogeneous chemistry on NAT and ice PSCs
        gppsc  = .TRUE.
! Use full not simplified scheme for PSCs.
        gpsimp = .FALSE.
      END IF
!
! pressure in hPa here.
      zp(1:n_points)       = p(1:n_points) / 100.0
      zt(1:n_points)       = t(1:n_points)

! copy tracers. Make sure they have been found correctly.
      IF (ihno3 > 0) THEN
        zhno3 = f(:,ihno3)
      ELSE
        zhno3 = 0.
      END IF
! if water vapour tracer is not present, use special water
! vapour field.
      IF (ih2o > 0) THEN
        zh2o = f(:,ih2o)
      ELSE
        zh2o = wp
      END IF
      IF (ihcl > 0) THEN
        zhcl = f(:,ihcl)
      ELSE
        zhcl = 0.
      END IF
      IF (iclono2 > 0) THEN
        zclono2 = f(:,iclono2)
      ELSE
        zclono2 = 0.
      END IF
      IF (in2o5 > 0) THEN
        zn2o5 = f(:,in2o5)
      ELSE
        zn2o5 = 0.
      END IF
      IF (ihocl > 0) THEN
        zhocl = f(:,ihocl)
      ELSE
        zhocl = 0.
      END IF

! Remove tropospheric ice clouds. They would cause model instability!
      WHERE (.NOT.(stratflag)) sph2o = 0.
!
! calculate the amount of hno3 and h2o in the solid phase and return
! the residual gas phase concentration
!
      CALL ukca_pscpres(zt(1:n_points),zp(1:n_points),tnd(1:n_points),  &
                    zh2o(1:n_points), zhno3(1:n_points), 1, n_points,   &
                    n_points, sph2o(1:n_points))

      IF (ihno3 > 0) f(:,ihno3) = zhno3
      IF (ih2o  > 0) f(:,ih2o)  = zh2o
!
      CALL ukca_calckpsc( za(1:n_points), zt(1:n_points),               &
                     zh2o(1:n_points), zhcl(1:n_points),                &
                     zclono2(1:n_points), zn2o5(1:n_points),            &
                     zhocl(1:n_points),                                 &
                     psc1(1:n_points), psc2(1:n_points),                &
                     psc3(1:n_points), psc4(1:n_points),                &
                     psc5(1:n_points), gpsa, gphocl,                    &
                     gppsc, gpsimp, n_points, 1, n_points, cdt )
!
! divide rates by h2o or hcl as asad treats psc reactions as bimolecular
!
      WHERE ( zhcl > peps )
        hk(:,2) = psc1 / zhcl
        hk(:,3) = psc5 / zhcl
        hk(:,5) = psc4 / zhcl
      ELSEWHERE
        hk(:,2) = 0.0
        hk(:,3) = 0.0
        hk(:,5) = 0.0
      ENDWHERE
      WHERE ( zh2o > peps )
        hk(:,1) = psc2 / zh2o
        hk(:,4) = psc3 / zh2o
      ELSEWHERE
        hk(:,1) = 0.0
        hk(:,4) = 0.0
      ENDWHERE
!
! copy the relevant hk's to rk's
! Introduce dynamical upper limit. Consider A + B -> C. Throughput through
! reaction rk*[A]*[B]*dt should be less than 0.5*min([A],[B])
! Also introduce flexible numbering (allow for reordering of reactions
! in rath.d
! Olaf Morgenstern  18/10/2004
! Do not do limiting in the case of non-families chemistry
!
! 1. ClONO2 + H2O --> HOCl + HNO3

      IF (n_clono2_h2o > 0) THEN
        rk(:,n_clono2_h2o) = hk(:,1)
      END IF

      IF (n_clono2_hcl > 0) THEN
! 2. ClONO2 + HCl --> Cl2 + HNO3
        rk(:,n_clono2_hcl) = hk(:,2)
      END IF

      IF (n_hocl_hcl > 0) THEN
! 3. HOCl + HCl --> Cl2 + H2O
        rk(:,n_hocl_hcl) = hk(:,3)
      END IF

      IF (n_n2o5_h2o > 0) THEN
! 4. N2O5 + H2O -> 2 HNO3
        rk(:,n_n2o5_h2o) = hk(:,4)
      END IF

      IF (n_n2o5_hcl > 0) THEN
! 5. N2O5 + HCl -> ClNO2 + HNO3
        rk(:,n_n2o5_hcl) = hk(:,5)
      END IF

! save the solid phase hno3 to add back after end of the
! chemistry timestep

      sphno3 = shno3

      IF (lhook) CALL dr_hook('UKCA_HETERO_MOD:UKCA_HETERO',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ukca_hetero

! ######################################################################
      SUBROUTINE ukca_solidphase(n_points)

! Description:
!
!**** *solidphase* - adds HONO2 and H2O in solid state back to main
!                    arrays.
!
!     Purpose.
!     --------
!         The amount of HONO2 and H2O in the solid phase when PSCs are
!     on is calculated and stored during a chemical step. Only the
!     amount still in the gas phase is used to calc. species tendencies.
!
!     Interface
!     ---------
!         This routine *MUST* be called at the end of each chemical
!     substep to add the solid phase HONO2 and H2O back to the main
!     ASAD arrays.
!
!     Method
!     ------
!          See comments in routine 'hetero'.
!
! The return of ice to the gasphase is disactivated if UM_ICE is set
! because the UM has an explicit ice tracer which we don't want to
! affect. Also, returning HNO3S is disactivated in case NAT PSC
! sedimentation is selected because that is done outside of ASAD.
!
!
! Code Description:
! Language: FORTRAN 90 + common extensions.
!
! Declarations:
! These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable
!
!---------------------------------------------------------------------
!

      USE ASAD_MOD,        ONLY: advt, f, sphno3
      USE ukca_option_mod, ONLY: jpctr
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Control_Max_Sizes
      IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: n_points

! local variables
      CHARACTER(LEN=72) :: cmessage
      INTEGER :: errcode                ! Variable passed to ereport
      INTEGER :: js
      INTEGER, SAVE :: ihno3=0
      LOGICAL, SAVE :: firstcall = .TRUE.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!     ----------------------------------------------------------------
!           1.  Add amount of HONO2 back to main ASAD arrays.
!               --- ------ -- ----- ---- -- ---- ---- -------
!
      IF (lhook) CALL dr_hook('UKCA_HETERO_MOD:UKCA_SOLIDPHASE',zhook_in,zhook_handle)
      IF (firstcall) THEN
        firstcall = .FALSE.
        DO js = 1, jpctr
          IF (advt(js) == 'HONO2     ')  ihno3 = js
        END DO
        IF (ihno3 == 0) THEN
          errcode=1
          cmessage='Select HONO2 as advected tracer.'
          CALL ereport('SOLIDPHASE',errcode,cmessage)
        END IF
      END IF

      f(1:n_points,ihno3) = f(1:n_points,ihno3) +                       &
                       sphno3(1:n_points)

      IF (lhook) CALL dr_hook('UKCA_HETERO_MOD:UKCA_SOLIDPHASE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE ukca_solidphase

! ######################################################################
      SUBROUTINE UKCA_CALCKPSC(SASA,T,TH2O,THCL,TCNIT,TN2O5,THOCL,      &
                          AKPSC1,AKPSC2,AKPSC3,AKPSC4,AKPSC5,           &
                          LPSA,LPHOCL,LPPSC,LPSIMP,                     &
                          KCHMLEV,KSTART,KEND,DT)
!
!     CALCKPSC - CALCULATION OF HETEROGENEOUS REACTION RATES
!
!         PETER GOOD, UNIVERSITY OF CAMBRIDGE, 9/11/93
!         BASED ON CODE BY LEE GRENFELL AND MARTYN CHIPPERFIELD
!
!     PURPOSE
!     -------
!         THIS ROUTINE CALCULATES FIRST ORDER RATES FOR REACTIONS
!     OCCURING ON TYPE 1 AND 2 PSC'S, AND AEROSOL.
!
!     INTERFACE
!     ---------
!         ARGUMENTS IN :  SASA   - (SULPHATE) AEROSOL SURFACE AREA
!                                  PER UNIT VOLUME (cm2 cm-3)
!                         T      - TEMPERATURE
!                         TH2O   - H2O NUMBER DENSITY (cm-3)
!                         THCL   - HCl NUMBER DENSITY (cm-3)
!                         TCNIT  - ClONO2 NUMBER DENSITY (cm-3)
!                         TN2O5  - N2O5 NUMBER DENSITY (cm-3)
!                         THOCL  - HOCl NUMBER DENSITY (cm-3)
!                         LPSA   - IF .TRUE., AEROSOL REACTIONS ARE ON
!                         LPHOCL - IF .TRUE., HOCL+HCL OCCURS ON AEROSOL
!                         LPPSC  - IF .TRUE., PSC REACTIONS ARE ON
!                         LPSIMP - IF .TRUE., USE SIMPLE PSC SCHEM
!                         KCHMLEV- LEVEL DIMENSION OF ARRAYS
!                         KSTART - TOP LEVEL OF CHEMISTRY
!                         KEND   - BOTTOM LEVEL OF CHEMISTRY
!                         DT     - MODEL TIMESTEP (s)
!
!     RESULTS
!     -------
!         FIRST ORDER RATE COEFFICIENTS AKPSC1 ... AKPSC5:
!
!           AKPSC1: HCl + ClONO2 -> Cl2 + HNO3 -> 2ClOX + HNO3
!           AKPSC2: ClONO2 + H2O -> HOCl + HNO3
!           AKPSC3: N2O5 + H2O   -> 2HNO3
!           AKPSC4: N2O5 + HCl   -> ClNO2 + HNO3
!           AKPSC5: HOCL + HCL   -> H2O + CL2    -> H2O + 2CLOX
!
!-----------------------------------------------------------------------
!
!     METHOD NOTES
!     ------ -----
!
!         AEROSOL REACTIONS 2 AND 5:
!
!              AKPSC2: ClONO2 + H2O -> HOCl + HNO3
!              The rate of this reaction is a strong function of the
!              sulphate acid wt.percent, and hence of temperature.
!
!              AKPSC5: HOCL + HCL -> CL2 + H2O
!              The rate of this reaction is proportional to the Henry's
!              law coefficients for HCl(aq) + Cl-(aq) and HOCl(aq)
!              That for HCl(aq) + Cl-(aq) is estimated from the sulphuri
!              acid wt.percent and temperature; the other is only known
!              for 60% w/w H2SO4.
!
!        See Cox, MacKenzie, Muller, Peter, Crutzen 1993
!
!----------------------------------------------------------------------
!
!     The collision frequency, v, of gas molecules with the
!     reacting surface is calculated using:
!
!     v=(A/4)*SQRT(8kT/pi*M)
!
!     A, Reacting surface concentration per unit volume.
!     k, Boltzmann's constant.
!     T, Temperature in kelvin.
!     M, Molecular mass of air molecules (=RMM*U).
!
!-----------------------------------------------------------------------
!
      USE ASAD_MOD,        ONLY: fpsc1, fpsc2, shno3, sh2o
      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Control_Max_Sizes
      IMPLICIT NONE

! Subroutine interface
      LOGICAL, INTENT(IN) :: lpsa
      LOGICAL, INTENT(IN) :: lphocl
      LOGICAL, INTENT(IN) :: lppsc
      LOGICAL, INTENT(IN) :: lpsimp

      INTEGER, INTENT(IN) :: kchmlev
      INTEGER, INTENT(IN) :: kstart
      INTEGER, INTENT(IN) :: kend

      REAL, INTENT(IN)  :: dt
      REAL, INTENT(IN)  :: t(kchmlev)
      REAL, INTENT(IN)  :: sasa(kchmlev)
      REAL, INTENT(IN)  :: thcl(kchmlev)
      REAL, INTENT(IN)  :: th2o(kchmlev)
      REAL, INTENT(IN)  :: tcnit(kchmlev)
      REAL, INTENT(IN)  :: tn2o5(kchmlev)
      REAL, INTENT(IN)  :: thocl(kchmlev)

      REAL, INTENT(OUT) :: AKPSC1(KCHMLEV)
      REAL, INTENT(OUT) :: AKPSC2(KCHMLEV)
      REAL, INTENT(OUT) :: AKPSC3(KCHMLEV)
      REAL, INTENT(OUT) :: AKPSC4(KCHMLEV)
      REAL, INTENT(OUT) :: AKPSC5(KCHMLEV)

! Local variables
      REAL, PARAMETER :: BOLTZ=1.38066E-23
      REAL, PARAMETER :: AVOGAD=6.02E+23
      REAL, PARAMETER :: PI=3.14159265
      REAL, PARAMETER :: U=1.66056E-27
!
      REAL, PARAMETER  :: RHO1=1.35
      REAL, PARAMETER  :: RHO2=0.928
      REAL, PARAMETER  :: RAD1=1.0E-4
      REAL, PARAMETER  :: RAD2=10.0E-4
      REAL, PARAMETER  :: RADSA=.1E-4
!
!     Bimolecular rate coefficient (HOCl+HCl*)(aq) (mol-1 m3 s-1)
      REAL, PARAMETER :: CPK1=1.0E+2
!
      REAL, PARAMETER :: GAM1A=0.3
      REAL, PARAMETER :: GAM1B=0.006
      REAL, PARAMETER :: GAM1C=0.0006
      REAL, PARAMETER :: GAM1D=0.003
      REAL, PARAMETER :: GAM1E=0.3
      REAL, PARAMETER :: GAM2A=0.3
      REAL, PARAMETER :: GAM2B=0.3
      REAL, PARAMETER :: GAM2C=0.03
      REAL, PARAMETER :: GAM2D=0.03
      REAL, PARAMETER :: GAM2E=0.3
      REAL, PARAMETER :: GAM3C=0.1
!
      INTEGER :: jl

      REAL :: zfcnit
      REAL :: zfn2o5
      REAL :: zfhocl
      REAL :: vol
      REAL :: order2
      REAL :: zrate
      REAL :: zdhcl
      REAL :: zfact

      REAL :: GAM3B(KCHMLEV)
      REAL :: HSHCL(KCHMLEV)
      REAL :: HHOCL(KCHMLEV)
      REAL :: CCNIT(KCHMLEV)
      REAL :: CN2O5(KCHMLEV)
      REAL :: CHOCL(KCHMLEV)
      REAL :: PSC1SA(KCHMLEV)
      REAL :: PSC2SA(KCHMLEV)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!
!
!-----------------------------------------------------------------------
!          1.  INITIALISE RATES TO ZERO
!              ---------- ----- -- ----
!
      IF (lhook) CALL dr_hook('UKCA_HETERO_MOD:UKCA_CALCKPSC',zhook_in,zhook_handle)
      AKPSC1(kstart:kend) = 0.0
      AKPSC2(kstart:kend) = 0.0
      AKPSC3(kstart:kend) = 0.0
      AKPSC4(kstart:kend) = 0.0
      AKPSC5(kstart:kend) = 0.0
!
!     If heterogeneous processes are on.
      IF ( LPPSC .OR. LPSA ) THEN
!
!-----------------------------------------------------------------------
!          2.  GENERAL TERMS IN COLLISION FREQUENCY EXPRESSION
!              ------- ----- -- --------- --------- ----------
!
         ZFCNIT=(8.0*BOLTZ)/(PI*97.5*U)
         ZFN2O5=(8.0*BOLTZ)/(PI*108.*U)
         ZFHOCL=(8.0*BOLTZ)/(PI*52.5*U)
!
         CCNIT(kstart:kend) = 0.25*SQRT(ZFCNIT*T(kstart:kend))
         CN2O5(kstart:kend) = 0.25*SQRT(ZFN2O5*T(kstart:kend))
         CHOCL(kstart:kend) = 0.25*SQRT(ZFHOCL*T(kstart:kend))
!
!-----------------------------------------------------------------------
!          3.  PSC RATES
!              --- -----
!
         IF (LPPSC) THEN
!
!              3.1  SIMPLE PSC SCHEME
!
         IF (LPSIMP) THEN
!
!         Zero order PSC rates.
           AKPSC1(kstart:kend)=4.6E-5*FPSC1(kstart:kend)
           AKPSC2(kstart:kend)=4.6E-5*FPSC1(kstart:kend)
           AKPSC3(kstart:kend)=4.6E-5*FPSC1(kstart:kend)
!
         ELSE
!
!              3.2  CALCULATE SURFACE AREA OF PSC'S
!
!           TYPE 1
           PSC1SA(kstart:kend)=                                         &
              1.85*63.0*U*3.0E5*SHNO3(kstart:kend)/(RHO1*RAD1)
!           TYPE 2
           PSC2SA(kstart:kend)=                                         &
                   18.0*U*3.0E5*SH2O (kstart:kend)/(RHO2*RAD2)
!
           AKPSC1(kstart:kend) =                                        &
             CCNIT(kstart:kend)*(PSC1SA(kstart:kend)*GAM1A +            &
                                    PSC2SA(kstart:kend)*GAM2A)
           AKPSC2(kstart:kend) =                                        &
             CCNIT(kstart:kend)*(PSC1SA(kstart:kend)*GAM1B +            &
                                    PSC2SA(kstart:kend)*GAM2B)
           AKPSC3(kstart:kend) =                                        &
             CN2O5(kstart:kend)*(PSC1SA(kstart:kend)*GAM1C +            &
                                    PSC2SA(kstart:kend)*GAM2C)
           AKPSC4(kstart:kend) =                                        &
             CN2O5(kstart:kend)*(PSC1SA(kstart:kend)*GAM1D +            &
                                    PSC2SA(kstart:kend)*GAM2D)
           AKPSC5(kstart:kend) =                                        &
             CHOCL(kstart:kend)*(PSC1SA(kstart:kend)*GAM1E +            &
                                    PSC2SA(kstart:kend)*GAM2E)
!
         END IF
!
         END IF
!
!----------------------------------------------------------------------
!          4.  AEROSOL REACTIONS
!              ------- ---------
!
!        If required, include the reactions on the sulphate aerosols.
         IF ( LPSA ) THEN
!
            CALL UKCA_EQCOMP(T,TH2O,KSTART,KEND,KCHMLEV,LPHOCL,         &
                          GAM3B,HSHCL,HHOCL)
!
            AKPSC2(kstart:kend) = AKPSC2(kstart:kend) +                 &
           CCNIT(kstart:kend)*100.0*SASA(kstart:kend)*GAM3B(kstart:kend)
            AKPSC3(kstart:kend) = AKPSC3(kstart:kend) +                 &
           CN2O5(kstart:kend)*100.0*SASA(kstart:kend)*GAM3C
!
!           If required, include HOCl + HCl -> Cl2 + H2O
            IF (LPHOCL)THEN
             DO JL = KSTART , KEND
!
!             Estimate specific volume from surface area
              VOL = SASA(JL)*RADSA/3.0
!
!             Add aerosol rate, converted to pseudo first order
              ORDER2=1.E6*CPK1*VOL*AVOGAD*HSHCL(JL)*HHOCL(JL)*          &
                           (BOLTZ*T(JL))**2
              AKPSC5(JL) = AKPSC5(JL) + ORDER2*THCL(JL)
             END DO
            END IF
!
!----------------------------------------------------------------------
!
         END IF
!
!        5.  CHECK FOR LOW HCL
!            ----- --- --- ---
!
         DO JL=KSTART, KEND
            ZRATE=MAX(1.,                                               &
                      DT*(AKPSC1(JL)*TCNIT(JL)+                         &
                      AKPSC4(JL)*TN2O5(JL)+                             &
                      AKPSC5(JL)*THOCL(JL)))
            ZDHCL=MIN(ZRATE,THCL(JL))
            ZFACT=ZDHCL/ZRATE
            AKPSC1(JL)=ZFACT*AKPSC1(JL)
            AKPSC4(JL)=ZFACT*AKPSC4(JL)
            AKPSC5(JL)=ZFACT*AKPSC5(JL)
         ENDDO
!
!----------------------------------------------------------------------
!
      END IF

      IF (lhook) CALL dr_hook('UKCA_HETERO_MOD:UKCA_CALCKPSC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CALCKPSC

! ######################################################################
      SUBROUTINE UKCA_EQCOMP(T,TH2O,KSTART,KEND,KCHMLEV,LPHOCL,         &
          RPCNCL,HSHCL,HHOCL)
!
!-----------------------------------------------------------------------
!
!     EQCOMP    - CALCULATIONS INVOLVING EQUILIBRIUM AEROSOL COMPOSITION
!
!            PETER GOOD, UNIVERSITY OF CAMBRIDGE, 10/11/93
!            BASED ON CODE BY SLIMANE BEKKI
!
!     PURPOSE
!     -------
!         CALCULATE THOSE TERMS IN THE AEROSOL REACTION RATE EXPRESSIONS
!     WHICH DEPEND ON THE AEROSOL COMPOSITION.
!
!     INTERFACE
!     ---------
!         ARGUMENTS IN :
!              T   -    TEMPERATURE (K)
!              TH2O   - H2O NUMBER DENSITY (cm-3)
!              KSTART - INDEX OF FIRST CHEMSITRY LEVEL
!              KEND   - LAST CHEMISTRY LEVEL
!              KCHMLEV- LEVEL DIMENSION OF ARRAYS
!              LPHOCL - IF .TRUE. THEN HOCL+HCL IN AEROSOL IS TURNED ON
!
!     RESULTS
!     -------
!              RPCNCL - REACTION PROBABILITY FOR ClONO2 + H2O
!              HSHCL  - MODIFIED HENRY'S LAW COEFFICIENT (mol/Nm)
!              HHOCL  - HENRY'S LAW COEFFICIENT
!
!     REFERENCES
!     ----------
!         HAMILL & STEELE  (TABLE FOR RPCNCL)
!         ZHANG ET AL. 1993   (DATA FOR EVALUATING HSHCL AND HHOCL)
!
!-----------------------------------------------------------------------
!
        USE ereport_mod, ONLY : ereport
        USE UM_ParVars
        USE Control_Max_Sizes
        IMPLICIT NONE
        
! Subroutine interface

      INTEGER, INTENT(IN)  :: kchmlev
      INTEGER, INTENT(IN)  :: kstart
      INTEGER, INTENT(IN)  :: kend
      LOGICAL, INTENT(IN)  :: LPHOCL
      REAL   , INTENT(IN)  :: T(KCHMLEV)
      REAL   , INTENT(IN)  :: TH2O(KCHMLEV)
      REAL   , INTENT(OUT) :: HSHCL(KCHMLEV)
      REAL   , INTENT(OUT) :: RPCNCL(KCHMLEV)
      REAL   , INTENT(OUT) :: HHOCL(KCHMLEV)

! Local variables
      INTEGER, PARAMETER :: NOPH2O=16
      INTEGER, PARAMETER :: NOTEMP =28
      INTEGER, PARAMETER :: NOCOMP= 4
      REAL   , PARAMETER :: BOLTZ=1.38066E-23
      INTEGER :: jxp1
      INTEGER :: jyp1
      INTEGER :: jx
      INTEGER :: jy
      INTEGER :: ierx
      INTEGER :: iery
      INTEGER :: i
      INTEGER :: jl
      REAL :: sxy
      REAL :: sx1y
      REAL :: sx1y1
      REAL :: sxy1
      REAL :: ta
      REAL :: tb
      REAL :: tt
      REAL :: ua
      REAL :: ub
      REAL :: u
      REAL :: a
      REAL :: b
      REAL :: ajx
      REAL :: ajxp1
      REAL :: bjx
      REAL :: bjxp1
      REAL :: ph2o
!
      REAL :: COMP(KCHMLEV)
      REAL :: COMPINDX(NOCOMP)
      REAL :: AINDEX(NOCOMP)
      REAL :: BINDEX(NOCOMP)
      REAL :: PINDX(NOPH2O)
      REAL :: TINDX(NOTEMP)
      REAL :: C(NOPH2O,NOTEMP)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!     Composition index for HSHCL
      DATA COMPINDX/35.0, 40.0, 50.0, 60.0/
!
!     Intercept and gradient look-up table for HSHCL
      DATA AINDEX/-11.35, -11.95, -13.80, -14.68/
      DATA BINDEX/6.10E+03, 5.88E+03, 5.53E+03, 4.92E+03/
!
!     P:Ambient partial pressure of water (mb)
      DATA PINDX/0.1000E-05 , 0.1000E-04 , 0.5000E-04 , 0.1000E-03 ,    &
           0.1500E-03 , 0.2000E-03 , 0.3000E-03 , 0.5000E-03 ,          &
           0.6000E-03 , 0.8000E-03 , 0.1000E-02 , 0.1200E-02 ,          &
           0.1500E-02 , 0.2000E-02 , 0.3000E-02 , 0.1000E-01/
!
!     T:Ambient temperature (K)
      DATA TINDX/0.1750E+03 , 0.1800E+03 , 0.1850E+03 , 0.1900E+03 ,    &
           0.1950E+03 , 0.2000E+03 , 0.2050E+03 , 0.2100E+03 ,          &
           0.2150E+03 , 0.2200E+03 , 0.2250E+03 , 0.2300E+03 ,          &
           0.2350E+03 , 0.2400E+03 , 0.2450E+03 , 0.2500E+03 ,          &
           0.2550E+03 , 0.2600E+03 , 0.2650E+03 , 0.2700E+03 ,          &
           0.2750E+03 , 0.2800E+03 , 0.2850E+03 , 0.2900E+03 ,          &
           0.2950E+03 , 0.3000E+03 , 0.3050E+03 , 0.3100E+03/
!
!     C:Composition (wt)
      DATA (C(I, 1),I=1,16)/                                            &
             0.9400E+02, 0.7800E+02, 0.5957E+02, 0.4363E+02,            &
             0.3430E+02, 0.2768E+02, 0.1836E+02, 0.4000E+02,            &
             0.4000E+02, 0.4000E+02, 0.4000E+02, 0.4000E+02,            &
             0.4000E+02, 0.4000E+02, 0.4000E+02, 0.9000E+01/
      DATA (C(I, 2),I=1,16)/                                            &
             0.9400E+02, 0.7900E+02, 0.6152E+02, 0.4723E+02,            &
             0.3887E+02, 0.3294E+02, 0.2458E+02, 0.1272E+02,            &
             0.4000E+02, 0.4000E+02, 0.1142E+02, 0.4000E+02,            &
             0.1341E+02, 0.4000E+02, 0.4000E+02, 0.9000E+01/
      DATA (C(I, 3),I=1,16)/                                            &
             0.9400E+02, 0.8000E+02, 0.6348E+02, 0.5084E+02,            &
             0.4344E+02, 0.3820E+02, 0.3080E+02, 0.1929E+02,            &
             0.1196E+02, 0.1681E+02, 0.1728E+02, 0.1300E+02,            &
             0.1928E+02, 0.1442E+02, 0.4000E+02, 0.9000E+01/
      DATA (C(I, 4),I=1,16)/                                            &
             0.9400E+02, 0.8100E+02, 0.6543E+02, 0.5444E+02,            &
             0.4801E+02, 0.4345E+02, 0.3702E+02, 0.2585E+02,            &
             0.1538E+02, 0.2541E+02, 0.2315E+02, 0.1806E+02,            &
             0.2515E+02, 0.1988E+02, 0.1615E+02, 0.9000E+01/
      DATA (C(I, 5),I=1,16)/                                            &
             0.9400E+02, 0.8200E+02, 0.6935E+02, 0.6165E+02,            &
             0.5715E+02, 0.5396E+02, 0.4946E+02, 0.4226E+02,            &
             0.3935E+02, 0.3402E+02, 0.2902E+02, 0.2313E+02,            &
             0.3102E+02, 0.2535E+02, 0.2334E+02, 0.9000E+01/
      DATA (C(I, 6),I=1,16)/                                            &
             0.9400E+02, 0.8300E+02, 0.7125E+02, 0.6594E+02,            &
             0.6283E+02, 0.6062E+02, 0.5751E+02, 0.5278E+02,            &
             0.5073E+02, 0.4693E+02, 0.4369E+02, 0.4086E+02,            &
             0.3689E+02, 0.3082E+02, 0.3052E+02, 0.9000E+01/
      DATA (C(I, 7),I=1,16)/                                            &
             0.9400E+02, 0.8367E+02, 0.7395E+02, 0.6976E+02,            &
             0.6731E+02, 0.6557E+02, 0.6312E+02, 0.5955E+02,            &
             0.5811E+02, 0.5561E+02, 0.5344E+02, 0.5144E+02,            &
             0.4863E+02, 0.4449E+02, 0.3771E+02, 0.9000E+01/
      DATA (C(I, 8),I=1,16)/                                            &
             0.9400E+02, 0.8420E+02, 0.7626E+02, 0.7284E+02,            &
             0.7084E+02, 0.6942E+02, 0.6742E+02, 0.6455E+02,            &
             0.6341E+02, 0.6147E+02, 0.5983E+02, 0.5838E+02,            &
             0.5646E+02, 0.5369E+02, 0.4849E+02, 0.1209E+02/
      DATA (C(I, 9),I=1,16)/                                            &
             0.9400E+02, 0.8519E+02, 0.7841E+02, 0.7548E+02,            &
             0.7377E+02, 0.7256E+02, 0.7085E+02, 0.6845E+02,            &
             0.6752E+02, 0.6594E+02, 0.6462E+02, 0.6347E+02,            &
             0.6196E+02, 0.5983E+02, 0.5640E+02, 0.3239E+02/
      DATA (C(I,10),I=1,16)/                                            &
             0.9400E+02, 0.8603E+02, 0.8020E+02, 0.7768E+02,            &
             0.7621E+02, 0.7517E+02, 0.7370E+02, 0.7163E+02,            &
             0.7083E+02, 0.6949E+02, 0.6839E+02, 0.6743E+02,            &
             0.6619E+02, 0.6447E+02, 0.6175E+02, 0.4271E+02/
      DATA (C(I,11),I=1,16)/                                            &
             0.9424E+02, 0.8691E+02, 0.8179E+02, 0.7959E+02,            &
             0.7830E+02, 0.7738E+02, 0.7609E+02, 0.7429E+02,            &
             0.7360E+02, 0.7244E+02, 0.7148E+02, 0.7066E+02,            &
             0.6960E+02, 0.6815E+02, 0.6589E+02, 0.5007E+02/
      DATA (C(I,12),I=1,16)/                                            &
             0.9433E+02, 0.8780E+02, 0.8323E+02, 0.8127E+02,            &
             0.8012E+02, 0.7930E+02, 0.7815E+02, 0.7656E+02,            &
             0.7595E+02, 0.7493E+02, 0.7410E+02, 0.7338E+02,            &
             0.7245E+02, 0.7119E+02, 0.6925E+02, 0.5567E+02/
      DATA (C(I,13),I=1,16)/                                            &
             0.9445E+02, 0.8860E+02, 0.8451E+02, 0.8275E+02,            &
             0.8172E+02, 0.8099E+02, 0.7996E+02, 0.7853E+02,            &
             0.7798E+02, 0.7708E+02, 0.7633E+02, 0.7570E+02,            &
             0.7489E+02, 0.7377E+02, 0.7207E+02, 0.6017E+02/
      DATA (C(I,14),I=1,16)/                                            &
             0.9478E+02, 0.8945E+02, 0.8571E+02, 0.8411E+02,            &
             0.8317E+02, 0.8250E+02, 0.8156E+02, 0.8027E+02,            &
             0.7977E+02, 0.7896E+02, 0.7829E+02, 0.7772E+02,            &
             0.7699E+02, 0.7600E+02, 0.7449E+02, 0.6392E+02/
      DATA (C(I,15),I=1,16)/                                            &
             0.9568E+02, 0.9057E+02, 0.8700E+02, 0.8546E+02,            &
             0.8456E+02, 0.8392E+02, 0.8302E+02, 0.8183E+02,            &
             0.8138E+02, 0.8063E+02, 0.8002E+02, 0.7951E+02,            &
             0.7885E+02, 0.7795E+02, 0.7659E+02, 0.6707E+02/
      DATA (C(I,16),I=1,16)/                                            &
             0.9695E+02, 0.9190E+02, 0.8836E+02, 0.8684E+02,            &
             0.8595E+02, 0.8532E+02, 0.8443E+02, 0.8327E+02,            &
             0.8284E+02, 0.8215E+02, 0.8158E+02, 0.8111E+02,            &
             0.8050E+02, 0.7969E+02, 0.7845E+02, 0.6977E+02/
      DATA (C(I,17),I=1,16)/                                            &
             0.9900E+02, 0.9374E+02, 0.9000E+02, 0.8840E+02,            &
             0.8746E+02, 0.8679E+02, 0.8585E+02, 0.8467E+02,            &
             0.8425E+02, 0.8357E+02, 0.8303E+02, 0.8258E+02,            &
             0.8202E+02, 0.8126E+02, 0.8012E+02, 0.7214E+02/
      DATA (C(I,18),I=1,16)/                                            &
             0.9900E+02, 0.9563E+02, 0.9170E+02, 0.9001E+02,            &
             0.8902E+02, 0.8832E+02, 0.8733E+02, 0.8610E+02,            &
             0.8566E+02, 0.8497E+02, 0.8444E+02, 0.8399E+02,            &
             0.8344E+02, 0.8272E+02, 0.8164E+02, 0.7408E+02/
      DATA (C(I,19),I=1,16)/                                            &
             0.9900E+02, 0.9753E+02, 0.9341E+02, 0.9163E+02,            &
             0.9059E+02, 0.8985E+02, 0.8881E+02, 0.8753E+02,            &
             0.8707E+02, 0.8637E+02, 0.8585E+02, 0.8540E+02,            &
             0.8486E+02, 0.8418E+02, 0.8316E+02, 0.7602E+02/
      DATA (C(I,20),I=1,16)/                                            &
             0.9900E+02, 0.9900E+02, 0.9511E+02, 0.9324E+02,            &
             0.9215E+02, 0.9138E+02, 0.9029E+02, 0.8896E+02,            &
             0.8848E+02, 0.8777E+02, 0.8726E+02, 0.8681E+02,            &
             0.8628E+02, 0.8564E+02, 0.8468E+02, 0.7796E+02/
      DATA (C(I,21),I=1,16)/                                            &
             0.9900E+02, 0.9900E+02, 0.9681E+02, 0.9486E+02,            &
             0.9372E+02, 0.9291E+02, 0.9177E+02, 0.9039E+02,            &
             0.8989E+02, 0.8917E+02, 0.8867E+02, 0.8822E+02,            &
             0.8770E+02, 0.8710E+02, 0.8620E+02, 0.7990E+02/
      DATA (C(I,22),I=1,16)/                                            &
             0.9900E+02, 0.9900E+02, 0.9851E+02, 0.9647E+02,            &
             0.9528E+02, 0.9444E+02, 0.9325E+02, 0.9182E+02,            &
             0.9130E+02, 0.9057E+02, 0.9008E+02, 0.8963E+02,            &
             0.8912E+02, 0.8856E+02, 0.8772E+02, 0.8184E+02/
      DATA (C(I,23),I=1,16)/                                            &
             0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9809E+02,            &
             0.9685E+02, 0.9597E+02, 0.9473E+02, 0.9325E+02,            &
             0.9271E+02, 0.9197E+02, 0.9149E+02, 0.9104E+02,            &
             0.9054E+02, 0.9002E+02, 0.8924E+02, 0.8378E+02/
      DATA (C(I,24),I=1,16)/                                            &
             0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02,            &
             0.9842E+02, 0.9750E+02, 0.9621E+02, 0.9468E+02,            &
             0.9412E+02, 0.9337E+02, 0.9290E+02, 0.9245E+02,            &
             0.9196E+02, 0.9148E+02, 0.9076E+02, 0.8572E+02/
      DATA (C(I,25),I=1,16)/                                            &
             0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02,            &
             0.9900E+02, 0.9900E+02, 0.9769E+02, 0.9611E+02,            &
             0.9553E+02, 0.9477E+02, 0.9431E+02, 0.9386E+02,            &
             0.9338E+02, 0.9294E+02, 0.9228E+02, 0.8766E+02/
      DATA (C(I,26),I=1,16)/                                            &
             0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02,            &
             0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9754E+02,            &
             0.9694E+02, 0.9617E+02, 0.9572E+02, 0.9527E+02,            &
             0.9480E+02, 0.9440E+02, 0.9380E+02, 0.8960E+02/
      DATA (C(I,27),I=1,16)/                                            &
             0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02,            &
             0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9897E+02,            &
             0.9835E+02, 0.9757E+02, 0.9713E+02, 0.9668E+02,            &
             0.9622E+02, 0.9586E+02, 0.9532E+02, 0.9154E+02/
      DATA (C(I,28),I=1,16)/                                            &
             0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02,            &
             0.9900E+02, 0.9900E+02, 0.9900E+02, 0.9900E+02,            &
             0.9900E+02, 0.9897E+02, 0.9854E+02, 0.9809E+02,            &
             0.9764E+02, 0.9732E+02, 0.9684E+02, 0.9348E+02/
!
!-----------------------------------------------------------------------
!          1.  REACTION PROBABILITY FOR CLONO2 + H2O
!              -------- ----------- --- ------ - ---
!
      IF (lhook) CALL dr_hook('UKCA_HETERO_MOD:UKCA_EQCOMP',zhook_in,zhook_handle)

      DO JL=KSTART, KEND

        PH2O=1.0E4*TH2O(JL)*BOLTZ*T(JL)
        CALL UKCA_POSITION(PINDX,NOPH2O,PH2O, JX,IERX,0)
        CALL UKCA_POSITION(TINDX,NOTEMP,T(JL),JY,IERY,0)
!
!     temperature first criteria
        IF ( JY == 0 ) THEN
          COMP(JL) = 9.0
        ELSE IF ( IERY == 1 ) THEN
          COMP(JL) = 95.0
        ELSE IF ( JX == 0 ) THEN
          COMP(JL) = 95.0
        ELSE IF ( IERX == 1 ) THEN
          COMP(JL) = 9.0
        ELSE
          JXP1 = JX + 1
          JYP1 = JY + 1
          SXY = C(JX,JY)
          SX1Y = C(JXP1,JY)
          SX1Y1 = C(JXP1,JYP1)
          SXY1 = C(JX,JYP1)
          TA = PH2O - PINDX(JX)
          TB = PINDX(JXP1) - PINDX(JX)
          TT = TA/TB
          UA = T(JL) - TINDX(JY)
          UB = TINDX(JYP1) - TINDX(JY)
          U = UA/UB
          COMP(JL) = (1.0-TT)*(1.0-U)*SXY + TT*(1.0-U)*SX1Y + TT*U*SX1Y1+&
                (1.0-TT)*U*SXY1

          IF ( COMP(JL)<9.0 ) COMP(JL) = 9.0
          IF ( COMP(JL)>95.0 ) COMP(JL) = 95.0
!
        END IF
!
!       WMO 1991
        RPCNCL(JL)=10.0**(1.87-(0.074*COMP(JL)))
        IF(RPCNCL(JL)>1.0) RPCNCL(JL)=1.0
      END DO
!
!-----------------------------------------------------------------------
!          2.  HENRY'S LAW COEFFICIENTS HSHCL AND HHOCL
!              ------- --- ------------ ----- --- -----
!
      IF (LPHOCL) THEN
        DO JL=KSTART,KEND
!
          CALL UKCA_POSITION(COMPINDX,NOCOMP,COMP(JL),JX,IERX,0)
!
          IF ( JX == 0 ) THEN
            A = -11.35
            B = 6.10E+03
          ELSE IF ( IERX == 1 ) THEN
            A = -14.68
            B = 4.92E+03
          ELSE
            JXP1 = JX + 1
            AJX = AINDEX(JX)
            AJXP1 = AINDEX(JXP1)
            BJX = BINDEX(JX)
            BJXP1 = BINDEX(JXP1)
            TA = COMP(JL) - COMPINDX(JX)
            TB = COMPINDX(JXP1) - COMPINDX(JX)
            TT = TA/TB
            A = (1.0-TT)*AJX + TT*AJXP1
            B = (1.0-TT)*BJX + TT*BJXP1
          END IF
!
          HSHCL(JL) = EXP(A + B/T(JL))/101.3
          HHOCL(JL) = EXP(0.71 + 1633.0/T(JL))/101.3
        END DO
      END IF
!
      IF (lhook) CALL dr_hook('UKCA_HETERO_MOD:UKCA_EQCOMP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_EQCOMP

! ######################################################################
      SUBROUTINE UKCA_POSITION(XC,N,X,JX,IER,IORDER)
!
!     Auxiliary subroutine for Eqcomp, Slimane Bekki, Nov. 1991
!
!-------------------------------------------------------------------------
!
        USE ereport_mod, ONLY : ereport
        USE UM_ParVars
        USE Control_Max_Sizes
        IMPLICIT NONE

! Subroutine interface
      INTEGER, INTENT(IN) :: n
      INTEGER, INTENT(IN) :: iorder
      REAL   , INTENT(IN) :: x
      REAL   , INTENT(IN) :: xc(n)

      INTEGER, INTENT(OUT) :: ier
      INTEGER, INTENT(OUT) :: jx

! Local variables
      INTEGER :: i

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
!     --------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook('UKCA_HETERO_MOD:UKCA_POSITION',zhook_in,zhook_handle)
      IER = 1
      IF ( IORDER == 0 ) THEN
        IF ( X < XC(1) ) THEN
          JX = 0
          ier = 0
        ELSE
          DO I = 1 , N
            IF ( X < XC(I) ) THEN
              ier = 0
              EXIT
            END IF
          END DO
          JX = I - 1
        END IF
      ELSE IF ( X > XC(1) ) THEN
        JX = 0
        ier = 0
      ELSE
        DO I = 1 , N
          IF ( X > XC(I) ) THEN
            ier = 0
            EXIT
          END IF
        ENDDO
        JX = I - 1
      END IF
!
!     --------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook('UKCA_HETERO_MOD:UKCA_POSITION',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_POSITION

! ######################################################################
      SUBROUTINE UKCA_PSCPRES(T,P,TND,TH2O,THNO3,                       &
                         KSTART,KEND,KCHMLEV, SPH2O)

      USE ASAD_MOD,           ONLY: shno3, sh2o, fpsc1, fpsc2
      USE ukca_d1_defs

      USE ereport_mod, ONLY : ereport
      USE UM_ParVars
      USE Control_Max_Sizes
      IMPLICIT NONE

! Subroutine interface

      INTEGER, INTENT(IN) :: kstart
      INTEGER, INTENT(IN) :: kend
      INTEGER, INTENT(IN) :: kchmlev
      REAL,    INTENT(IN) :: T(KCHMLEV)
      REAL,    INTENT(IN) :: P(KCHMLEV)
      REAL,    INTENT(IN) :: TND(KCHMLEV)
      REAL,    INTENT(INOUT) :: TH2O(KCHMLEV)
      REAL,    INTENT(INOUT) :: THNO3(KCHMLEV)
      REAL,    INTENT(INOUT) :: SPH2O(KCHMLEV)

! Local variables
      INTEGER :: j
      INTEGER :: jl
      REAL :: zh2ot
      REAL :: ztpsc
      REAL :: zmt
      REAL :: zbt
      REAL :: zhno3eq
      REAL :: zh2oeq

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
      IF (lhook) CALL dr_hook('UKCA_HETERO_MOD:UKCA_PSCPRES',zhook_in,zhook_handle)
      SHNO3(kstart:kend) = 0.0
      SH2O (kstart:kend) = 0.0
      FPSC1(kstart:kend) = 0.0
      FPSC2(kstart:kend) = 0.0
!
      DO JL=KSTART,KEND
!
!     Type 1 PSCs
!
        ZH2OT = 0.75*P(JL)*TH2O(JL)/TND(JL)
        ZH2OT = MAX (ZH2OT,1.0E-12)
        ZTPSC = T(JL)
        ZMT     = -2.7836 - 0.00088*ZTPSC
        ZBT     = 38.9855 - 11397.0/ZTPSC + 0.009179*ZTPSC
        ZHNO3EQ = 10.0**(ZMT*ALOG10(ZH2OT)+ ZBT)
!
!     Convert to number density.
!
        ZHNO3EQ = ZHNO3EQ*133.3*TND(JL)/(100.0*P(JL))
!
        IF (THNO3(JL) > ZHNO3EQ) THEN
          FPSC1(JL) = 1.0
          SHNO3(JL) = THNO3(JL)-ZHNO3EQ
          THNO3(JL) = ZHNO3EQ
        ELSE
          FPSC1(JL) = 0.0
        END IF
!
!     Type 2 PSCs
!
        IF (.TRUE.) THEN
! just sh2o from volume mixing ratio to number density and set
! FPSC2 flag
          IF (sph2o(jl) > 0.) THEN
            sh2o(jl) = sph2o(jl) * tnd(jl)
            fpsc2(jl) = 1.
          ELSE
            fpsc2(jl) = 0.
            sh2o(jl) = 0.
          END IF
        ELSE
! calculate water ice number density locally
          ZH2OEQ=610.78*EXP(21.875*(T(JL)-273.16)/(T(JL)-7.66))
!
!     Convert to number density.
!
          ZH2OEQ = ZH2OEQ*TND(JL)/(100.0*P(JL))
!
          IF (TH2O(JL) > ZH2OEQ) THEN
            FPSC2(JL) = 1.0
            SH2O(JL) = TH2O(JL)-ZH2OEQ
            TH2O(JL)  = ZH2OEQ
          ELSE
            FPSC2(JL) = 0.0
          END IF
        END IF
      END DO
!
      IF (lhook) CALL dr_hook('UKCA_HETERO_MOD:UKCA_PSCPRES',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_PSCPRES
! =======================================================================
      END MODULE
