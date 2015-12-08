! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!  To calculate the reaction rate co-efficients for use in the
!  Backward Euler solver
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method: Calculates rate coefficients
!                 Units : cm, molecule, s
!
!   Inputs  : tc,m,h2o,o2
!   Outputs : rc
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------
!
      SUBROUTINE UKCA_CHEMCO(nr,n_pnts,tc,m,h2o,o2,                     &
                             clw,fcloud,fdiss,k_dms,rc)

      USE UKCA_CONSTANTS,   ONLY: avc, H_plus, rhow, m_air
      USE ukca_option_mod,  ONLY: L_ukca_aerchem, jpspec
      USE ASAD_MOD,         ONLY: peps, jpeq
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

      INTEGER, INTENT(IN) :: nr                ! No of reactions
      INTEGER, INTENT(IN) :: n_pnts            ! No of points

      REAL, INTENT(IN)    :: tc(n_pnts)        ! Temperature
      REAL, INTENT(IN)    :: m(n_pnts)         ! Air density
      REAL, INTENT(IN)    :: h2o(n_pnts)       ! Water vapour
      REAL, INTENT(IN)    :: o2(n_pnts)        ! Oxygen
      REAL, INTENT(IN)    :: clw(n_pnts)       ! Cloud Liquid water
      REAL, INTENT(IN)    :: fcloud(n_pnts)    ! Cloud fraction
! Dissolved fraction - allows for jpeq+1 ions:
      REAL, INTENT(IN)    :: fdiss(n_pnts,jpspec,jpeq+1)

! Rate coeffs to allow parameterisation of DMS products:
      REAL, INTENT(OUT)   :: k_dms(n_pnts,5)

!     Thermal rate coefficients:
      REAL, INTENT(OUT)   :: rc(n_pnts,nr)

!     Local variables

      INTEGER :: i, j

      REAL :: vr(n_pnts)                ! Volume ratio of cloud water
      REAL :: z1(n_pnts)
      REAL :: z2(n_pnts) 
      REAL :: z3(n_pnts) 
      REAL :: z4(n_pnts) 
      REAL :: ratiob2total(n_pnts) 
      REAL :: ratioa2b(n_pnts) 
      REAL :: at(7)
      REAL :: t300(n_pnts)
      REAL :: zo(n_pnts), zi(n_pnts), arg2(n_pnts)
      REAL :: zfc(n_pnts), zr(n_pnts)
      REAL :: faq(n_pnts,jpspec)           ! total dissolved fraction

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



      IF (lhook) CALL dr_hook('UKCA_CHEMCO',zhook_in,zhook_handle)
      IF (L_ukca_aerchem) THEN
        faq(:,:) = fdiss(:,:,1) + fdiss(:,:,2)
! convert clw in kg/kg to volume ratio
        vr(:) = clw(:)*m(:)*m_air*(1e6/avc)/rhow
      ENDIF

!     Initialise reaction rate coefficients 

      DO j = 1, nr
        DO i = 1, n_pnts
          rc(i,j) = 0.0
        END DO
      END DO

!     Calculate bimolecular reaction rate coefficients

!     R1: HO2 + NO = OH + NO2              3.60E-12  0.00   -270.0
!     IUPAC no change

      RC(:,1) = 3.60E-12*exp(270.0/TC(:))

!     R2: HO2 + NO3 = OH + NO2             4.00E-12  0.00      0.0
!     IUPAC no change

      RC(:,2) = 4.00E-12

!     R3: HO2 + O3 = OH + O2               2.03E-16  4.57   -693.0
!     IUPAC no change

      RC(:,3) = 2.03E-16*(tc(:)/300.0)**4.57*exp(693.0/tc(:))

!     R4: HO2 + HO2 = H2O2                 2.20E-13  0.00   -600.0
!     IUPAC no change - Note 4

      RC(:,4) = 2.20E-13*exp(600.0/tc(:))
      RC(:,4) = RC(:,4)*(1.0 +                                       &
     &          1.4E-21*h2o(:)*exp(2200.0/tc(:)) )

!     R5: HO2 + MeOO = MeOOH               3.80E-13  0.00   -780.0
!     IUPAC 2005 update - Note 5 - Brnch A

      ratiob2total(:) = 1.0 / (1.0 + 498.0*EXP(-1160.0/tc(:)))
      RC(:,5)         = 3.80E-13*exp(780.0/tc(:))
      RC(:,5)         = RC(:,5)*(1.0-ratiob2total(:))

!     R6: HO2 + MeOO = HCHO                3.80E-13  0.00   -780.0
!     IUPAC 2005 update - Note 6 - Brnch B

      RC(:,6) = 3.80E-13*exp(780.0/tc(:))
      RC(:,6) = RC(:,6)*ratiob2total(:)

!     R7: HO2 + EtOO = EtOOH               3.80E-13  0.00   -900.0
!     IUPAC no change

      RC(:,7) = 3.80E-13*exp(900.0/tc(:))

!     R8: HO2 + MeCO3 = MeCO3H             2.08E-13  0.00   -980.0
!     IUPAC 2005 update - Note 8

      RC(:,8) = 2.08E-13*exp(980.0/tc(:))

!     R9: HO2 + MeCO3 = MeCO2H + O3        1.04E-13  0.00   -980.0
!     IUPAC 2005 update - note 9

      RC(:,9) = 1.04E-13*exp(980.0/tc(:))

!     R10 HO2 + MeCO3 = OH + MeOO          2.08E-13  0.00   -980.0
!     IUPAC 2005 update - Note 10

      RC(:,10) = 2.08E-13*exp(980.0/tc(:))

!     R11 HO2 + n-PrOO = n-PrOOH           1.51E-13  0.00  -1300.0
!     MCM no change

      RC(:,11) = 1.51E-13*exp(1300.0/tc(:))

!     R12 HO2 + i-PrOO = i-PrOOH           1.51E-13  0.00  -1300.0
!     MCM no change

      RC(:,12) = 1.51E-13*exp(1300.0/tc(:))

!     R13 HO2 + EtCO3 = O2 + EtCO3H        3.05E-13  0.00  -1040.0
!     MCM no change

      RC(:,13) = 3.05E-13*exp(1040.0/tc(:))

!     R14 HO2 + EtCO3 = O3 + EtCO2H        1.25E-13  0.00  -1040.0
!     MCM no change

      RC(:,14) = 1.25E-13*exp(1040.0/tc(:))

!     R15 HO2 + MeCOCH2OO = MeCOCH2OOH     1.36E-13  0.00  -1250.0
!     MCM no change - Note 15

      RC(:,15) = 1.36E-13*exp(1250.0/tc(:))

!     R16 MeOO + NO = HO2 + HCHO + NO2     2.95E-12  0.00   -285.0
!     IUPAC no change

      RC(:,16) = 2.95E-12*exp(285.0/tc(:))

!     R17 MeOO + NO = MeONO2               2.95E-15  0.00   -285.0
!     IUPAC no change - note 17

      RC(:,17) = 2.95E-15*exp(285.0/tc(:))

!     R18 MeOO + NO3 = HO2 + HCHO + NO2    1.30E-12  0.00      0.0
!     IUPAC no change

      RC(:,18) = 1.30E-12

!     R19 MeOO + MeOO = MeOH + HCHO        1.03E-13  0.00   -365.0
!     IUPAC no change - note 19 - Brnch A

      ratiob2total(:) = 1.0/(1.0+(EXP(1300.0/tc(:)))/33.0)
      RC(:,19)        = 1.03E-13*EXP(365.0/tc(:))
      RC(:,19)        = RC(:,19)*(1.0-ratiob2total(:))

!     R20 MeOO + MeOO = 2HO2 + 2HCHO       1.03E-13  0.00   -365.0
!     IUPAC no change - note 20 - Brnch B

      RC(:,20) = 1.03E-13*EXP(365.0/tc(:))
      RC(:,20) = RC(:,20)*ratiob2total(:)

!     R21 MeOO + MeCO3 = HO2 + HCHO + MeOO 1.80E-12  0.00   -500.0
!     IUPAC no change - note 21

      RC(:,21) = 1.80E-12*EXP(500.0/tc(:))

!     R22 MeOO + MeCO3 = MeCO2H + HCHO     2.00E-13  0.00   -500.0
!     IUPAC no change - note 22

      RC(:,22) = 2.00E-13*EXP(500.0/tc(:))

!     R23 EtOO + NO = MeCHO + HO2 + NO2    2.60E-12  0.00   -380.0
!     IUPAC no change - note 23

      RC(:,23) = 2.60E-12*EXP(380.0/tc(:))

!     R24 EtOO + NO3 = MeCHO + HO2 + NO2   2.30E-12  0.00      0.0
!     IUPAC no change

      RC(:,24) = 2.30E-12

!     R25 EtOO + MeCO3 = MeCHO + HO2 + MeOO 4.40E-13  0.00  -1070.0
!     IUPAC no change - note 25

      RC(:,25) = 4.40E-13*EXP(1070.0/tc(:))

!     R26 MeCO3 + NO = MeOO + CO2 + NO2    7.50E-12  0.00   -290.0
!     IUPAC no change

      RC(:,26) = 7.50E-12*EXP(290.0/tc(:))

!     R27 MeCO3 + NO3 = MeOO + CO2 + NO2   4.00E-12  0.00      0.0
!     MCM no change

      RC(:,27) = 4.00E-12

!     R28 n-PrOO + NO = EtCHO + HO2 + NO2  2.90E-12  0.00   -350.0
!     IUPAC no change - note 28

      RC(:,28) = 2.90E-12*EXP(350.0/tc(:))

!     R29 n-PrOO + NO3 = EtCHO + HO2 + NO2 2.50E-12  0.00      0.0
!     MCM no change

      RC(:,29) = 2.50E-12

!     R30 i-PrOO + NO = Me2CO + HO2 + NO2  2.70E-12  0.00   -360.0
!     IUPAC no change - note 30

      RC(:,30) = 2.70E-12*EXP(360.0/tc(:))

!     R31 i-PrOO + NO3 = Me2CO + HO2 + NO2 2.50E-12  0.00      0.0
!     MCM no change

      RC(:,31) = 2.50E-12

!     R32 EtCO3 + NO = EtOO + CO2 + NO2    6.70E-12  0.00   -340.0
!     IUPAC no change

      RC(:,32) = 6.70E-12*EXP(340.0/tc(:))

!     R33 EtCO3 + NO3 = EtOO + CO2 + NO2   4.00E-12  0.00      0.0
!     MCM no change

      RC(:,33) = 4.00E-12

!     R34 MeCOCH2OO + NO = MeCO3+HCHO+NO2  2.80E-12  0.00   -300.0
!     Tyndall et al. - note 34

      RC(:,34) = 2.80E-12*EXP(300.0/tc(:))

!     R35 MeCOCH2OO + NO3 = MeCO3+HCHO+NO2 2.50E-12  0.00      0.0
!     MCM no change

      RC(:,35) = 2.50E-12

!     R36 NO + NO3 = NO2 + NO2             1.80E-11  0.00   -110.0
!     IUPAC no change

      RC(:,36) = 1.80E-11*EXP(110.0/tc(:))

!     R37 NO + O3 = NO2                    1.40E-12  0.00   1310.0
!     IUPAC no change

      RC(:,37) = 1.40E-12*EXP(-1310.0/tc(:))

!     R38 NO2 + O3 = NO3                   1.40E-13  0.00   2470.0
!     IUPAC no change

      RC(:,38) = 1.40E-13*EXP(-2470.0/tc(:))

!     R39 NO3 + HCHO = HONO2 + HO2 + CO    2.00E-12  0.00   2440.0
!     IUPAC no change - note 39

      RC(:,39) = 2.00E-12*EXP(-2440.0/tc(:))

!     R40 NO3 + MeCHO = HONO2 + MeCO3      1.40E-12  0.00   1860.0
!     IUPAC no change

      RC(:,40) = 1.40E-12*EXP(-1860.0/tc(:))

!     R41 NO3 + EtCHO = HONO2 + EtCO3      3.46E-12  0.00   1862.0
!     MCM no change - note 41

      RC(:,41) = 3.46E-12*EXP(-1862.0/tc(:))

!     R42 NO3 + Me2CO = HONO2 + MeCOCH2OO  3.00E-17  0.00      0.0
!     IUPAC no change - note 42

      RC(:,42) = 3.00E-17

!     R43 N2O5 + H2O = HONO2 + HONO2       2.50E-22  0.00      0.0
!     IUPAC no change - note 43

      RC(:,43) = 2.50E-22

!     R44 O(3P) + O3 = O2 + O2             8.00E-12  0.00   2060.0
!     IUPAC no change

      RC(:,44) = 8.00E-12*EXP(-2060.0/tc(:))

!     R45 O(1D) + CH4 = OH + MeOO          1.05E-10  0.00      0.0
!     IUPAC no change - note 45

      RC(:,45) = 1.05E-10

!     R46 O(1D) + CH4 = HCHO + H2          7.50E-12  0.00      0.0
!     IUPAC no change - note 46

      RC(:,46) = 7.50E-12

!     R47 O(1D) + CH4 = HCHO + HO2 + HO2   3.45E-11  0.00      0.0
!     IUPAC no change - note 47

      RC(:,47) = 3.45E-11

!     R48 O(1D) + H2O = OH + OH            2.20E-10  0.00      0.0
!     IUPAC no change

      RC(:,48) = 2.20E-10

!     R49 O(1D) + N2 = O(3P) + N2          2.10E-11  0.00   -115.0
!     Ravishankara et al. - note 49

      RC(:,49) = 2.10E-11*EXP(115.0/tc(:))

!     R50 O(1D) + O2 = O(3P) + O2          3.20E-11  0.00    -67.0
!     IUPAC no change

      RC(:,50) = 3.20E-11*EXP(67.0/tc(:))

!     R51 OH + CH4 = H2O + MeOO            1.85E-12  0.00   1690.0
!     IUPAC no change

      RC(:,51) = 1.85E-12*EXP(-1690.0/tc(:))

!     R52 OH + C2H6 = H2O + EtOO           6.90E-12  0.00   1000.0
!     IUPAC no change

      RC(:,52) = 6.90E-12*EXP(-1000.0/tc(:))

!     R53 OH + C3H8 = n-PrOO + H2O         7.60E-12  0.00    585.0
!     IUPAC no change - note 53 - Brnch A

      ratioa2b(:) = 226.0 * tc(:)**(-0.64)*EXP(-816.0/tc(:))
      RC(:,53) = 7.60E-12*EXP(-585.0/tc(:))
      RC(:,53) = RC(:,53)*(ratioa2b(:)/(ratioa2b(:)+1.0))

!     R54 OH + C3H8 = i-PrOO + H2O         7.60E-12  0.00    585.0
!     IUPAC no change - note 54 - Brnch B

      RC(:,54) = 7.60E-12*EXP(-585.0/tc(:))
      RC(:,54) = RC(:,54)/(ratioa2b(:)+1.0)

!     R55 OH + CO = HO2                    1.44E-13  0.00      0.0
!     IUPAC 2005 update - note 55

      RC(:,55) = 1.44E-13*(1.0+m(:)/4.2E19)

!     R56 OH + EtCHO = H2O + EtCO3         5.10E-12  0.00   -405.0
!     IUPAC no change

      RC(:,56) = 5.10E-12*EXP(405.0/tc(:))

!     R57 OH + EtOOH = H2O + MeCHO + OH    8.01E-12  0.00      0.0
!     MCM no change

      RC(:,57) = 8.01E-12

!     R58 OH + EtOOH = H2O + EtOO          1.90E-12  0.00   -190.0
!     MCM no change

      RC(:,58) = 1.90E-12*EXP(190.0/tc(:))

!     R59 OH + H2 = H2O + HO2              7.70E-12  0.00   2100.0
!     IUPAC no change

      RC(:,59) = 7.70E-12*EXP(-2100.0/tc(:))

!     R60 OH + H2O2 = H2O + HO2            2.90E-12  0.00    160.0
!     IUPAC no change

      RC(:,60) = 2.90E-12*EXP(-160.0/tc(:))

      IF (L_ukca_aerchem) THEN

!       Reduce rate to account for dissolved H2O2
        RC(:,60) = RC(:,60)*(1.0-(faq(:,13)*fcloud(:)))

      ENDIF

!     R61 OH + HCHO = H2O + HO2 + CO       5.40E-12  0.00   -135.0
!     IUPAC 2004 update - note 61

      RC(:,61) = 5.40E-12*EXP(135.0/tc(:))

!     R62 OH + HO2 = H2O                   4.80E-11  0.00   -250.0
!     IUPAC no change

      RC(:,62) = 4.80E-11*EXP(250.0/tc(:))

!     R63 OH + HO2NO2 = H2O + NO2          1.90E-12  0.00   -270.0
!     IUPAC no change

      RC(:,63) = 1.90E-12*EXP(270.0/tc(:))

!     R64 OH + HONO2 = H2O + NO3           1.50E-13  0.00      0.0
!     IUPAC no change - note 64

      z1(:)    = 2.4e-14 * exp( 460.0/tc(:))
      z3(:)    = 6.5e-34 * exp(1335.0/tc(:))
      z4(:)    = 2.7e-17 * exp(2199.0/tc(:))
      z2(:)    = z3(:)*m(:) / (1.0+z3(:)*m(:)/z4(:))
      RC(:,64) = z1(:) + z2(:)

!     R65 OH + HONO = H2O + NO2            2.50E-12  0.00   -260.0
!     IUPAC no change

      RC(:,65) = 2.50E-12*exp(260.0/tc(:))

!     R66 OH + MeOOH = H2O + HCHO + OH     1.02E-12  0.00   -190.0
!     IUPAC no change - note 66

      RC(:,66) = 1.02E-12*EXP(190.0/tc(:))

!     R67 OH + MeOOH = H2O + MeOO          1.89E-12  0.00   -190.0
!     IUPAC no change - note 67

      RC(:,67) = 1.89E-12*EXP(190.0/tc(:))

!     R68 OH + MeONO2 = HCHO + NO2 + H2O   4.00E-13  0.00    845.0
!     IUPAC no change

      RC(:,68) = 4.00E-13*EXP(-845.0/tc(:))

!     R69 OH + Me2CO = H2O + MeCOCH2OO     8.80E-12  0.00   1320.0
!     IUPAC no change - note 69

      RC(:,69) = 8.80E-12*EXP(-1320.0/tc(:))
      
!     R70 OH + Me2CO = H2O + MeCOCH2OO     1.70E-14  0.00   -420.0
!     IUPAC no change - note 70

      RC(:,70) = 1.70E-14*EXP(420.0/tc(:))

!     R71 OH + MeCOCH2OOH = H2O+MeCOCH2OO  1.90E-12  0.00   -190.0
!     MCM no change - note 71

      RC(:,71) = 1.90E-12*EXP(190.0/tc(:))

!     R72 OH + MeCOCH2OOH = OH + MGLY      8.39E-12  0.00      0.0
!     MCM no change - note 72

      RC(:,72) = 8.39E-12
      
!     R73 OH + MeCHO = H2O + MeCO3         4.40E-12  0.00   -365.0
!     IUPAC no change

      RC(:,73) = 4.40E-12*EXP(365.0/tc(:))

!     R74 OH + NO3 = HO2 + NO2             2.00E-11  0.00      0.0
!     IUPAC no change

      RC(:,74) = 2.00E-11

!     R75 OH + O3 = HO2 + O2               1.70E-12  0.00    940.0
!     IUPAC no change

      RC(:,75) = 1.70E-12*EXP(-940.0/tc(:))

!     R76 OH + OH = H2O + O(3P)            6.31E-14  2.60   -945.0
!     IUPAC no change - note 76

      RC(:,76) = 6.31E-14*((tc(:)/300.0)**2.6)*EXP(945.0/tc(:))

!     R77 OH + PAN = HCHO + NO2 + H2O      3.00E-14  0.00      0.0
!     IUPAC no change - note 77

      RC(:,77) = 3.00E-14

!     R78 OH + PPAN = MeCHO + NO2 + H2O    1.27E-12  0.00      0.0
!     MCM no change

      RC(:,78) = 1.27E-12

!     R79 OH + n-PrOOH = n-PrOO + H2O      1.90E-12  0.00   -190.0
!     MCM no change

      RC(:,79) = 1.90E-12*EXP(190.0/tc(:))

!     R80 OH + n-PrOOH = EtCHO + H2O + OH  1.10E-11  0.00      0.0  
!     MCM no change

      RC(:,80) = 1.10E-11

!     R81 OH + i-PrOOH = i-PrOO + H2O      1.90E-12  0.00   -190.0
!     MCM no change

      RC(:,81) = 1.90E-12*EXP(190.0/tc(:))

!     R82 OH + i-PrOOH = Me2CO + OH        1.66E-11  0.00      0.0
!     MCM no change

      RC(:,82) = 1.66E-11

!     R83 O(3P) + NO2 = NO + O2            5.50E-12  0.00   -188.0
!     IUPAC no change

      RC(:,83) = 5.50E-12*EXP(188.0/tc(:))

!     R84 HO2S + O3S = HO2S + O2           2.03E-16  4.57   -693.0
!     IUPAC no change - note 84

      RC(:,84) = 2.03E-16*((tc(:)/300.0)**4.57)*EXP(693.0/tc(:))

!     R85 OHS + O3S = OHS + O2             1.70E-12  0.00    940.0
!     IUPAC no change

      RC(:,85) = 1.70E-12*EXP(-940.0/tc(:))

!     R86 O(1D)S + H2O = H2O               2.20E-10  0.00      0.0
!     IUPAC no change

      RC(:,86) = 2.20E-10

!     R87 O(1D)S + N2 = O(3P)S + N2        2.10E-11  0.00   -115.0
!     Ravishankara et al. - note 49

      RC(:,87) = 2.10E-11*EXP(115.0/tc(:))

!     R88 O(1D)S + O2 = O(3P)S + O2        3.20E-11  0.00    -67.0
!     IUPAC no change

      RC(:,88) = 3.20E-11*EXP(67.0/tc(:))
      
!     R89 HO2 + HO2 = H2O2 + O2     0.00 1.90E-33  0.00   -980.0 0.00E+00  0.00      0.0
!     IUPAC no change - note 1

      at(1) = 0.00
      at(2) = 1.90E-33
      at(3) = 0.00
      at(4) = -980.0
      at(5) = 0.00E+00
      at(6) = 0.00
      at(7) = 0.0

      t300(:) = tc(:)/300.0
      zo(:)   = at(2)*EXP(at(3)*LOG(t300(:)))*EXP(-at(4)/tc(:))*m(:)
      zi(:)   = at(5)*EXP(at(6)*LOG(t300(:)))*EXP(-at(7)/tc(:))

      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,89) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,89) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((log10(zr(i)))**2))
          RC(i,89) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO
      RC(:,89) = RC(:,89)*(1.0 + 1.4E-21*h2o(:)*exp(2200.0/tc(:)) )

      
!     R90 HO2 + NO2 = HO2NO2 + m    0.60 1.80E-31 -3.20  0.0 4.70E-12  0.00      0.0
!     IUPAC no change

      at(1) = 0.60
      at(2) = 1.80E-31
      at(3) = -3.20
      at(4) = 0.0
      at(5) = 4.70E-12
      at(6) = 0.00
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))
      
      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,90) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,90) = zo(i)
        ELSE
          IF( at(1) <= 1.0 )THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = EXP( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,90) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R91 HO2NO2 + m = HO2 + NO2  0.60 4.10E-05  0.00  10650.0 4.80E+15  0.00  11170.0
!     IUPAC no change

      at(1) = 0.60
      at(2) = 4.10E-5
      at(3) = 0.0
      at(4) = 10650.0
      at(5) = 4.80E+15
      at(6) = 0.00
      at(7) = 11170.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))
      
      DO i=1,n_pnts
        IF( zo(i) < peps )THEN
          RC(i,91) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,91) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,91) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R92 MeCO3 + NO2 = PAN + m  0.30 2.70E-28 -7.10   0.0 1.20E-11 -0.90      0.0
!     IUPAC no change

      at(1) = 0.30
      at(2) = 2.70E-28
      at(3) = -7.10
      at(4) = 0.0
      at(5) = 1.21E-11
      at(6) = -0.90
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))
      
      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,92) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,92) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,92) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R93 PAN + m = MeCO3 + NO2    0.30 4.90E-03  0.00  12100.0 5.40E+16  0.00  13830.0
!     IUPAC no change

      at(1) = 0.30
      at(2) = 4.90E-03
      at(3) = 0.00
      at(4) = 12100.0
      at(5) = 5.40E+16
      at(6) = 0.00
      at(7) = 13830.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))
      
      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,93) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,93) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,93) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R94 N2O5 + m = NO2 + NO3  0.35 1.30E-03 -3.50  11000.0 9.70E+14  0.10  11080.0
!     IUPAC no change - note 6

      at(1) = 0.35
      at(2) = 1.30E-03
      at(3) = -3.50
      at(4) = 11000.0
      at(5) = 9.70E+14
      at(6) = 0.10
      at(7) = 11080.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,94) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,94) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,94) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R95 NO2 + NO3 = N2O5 + m   0.35 3.60E-30 -4.10   0.0 1.90E-12  0.20      0.0
!     IUPAC no change - note 7

      at(1) = 0.35
      at(2) = 3.60E-30
      at(3) = -4.10
      at(4) = 0.0
      at(5) = 1.90E-12
      at(6) = 0.20
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,95) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,95) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,95) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R96 O(3P) + O2 = O3 + m      0.00 5.70E-34 -2.60  0.0 0.00E+00  0.00      0.0
!     IUPAC no change - note 8

      at(1) = 0.00
      at(2) = 5.70E-34
      at(3) = -2.60
      at(4) = 0.0
      at(5) = 0.00E+00
      at(6) = 0.00
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF ( zo(i) < peps ) THEN
          RC(i,96) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,96) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,96) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

! R97 OH + NO = HONO + m       1420.00 7.40E-31 -2.40      0.0 3.30E-11 -0.30      0.0
! IUPAC no change - note 9

      at(1) = 1420.0
      at(2) = 7.40E-31
      at(3) = -2.40
      at(4) = 0.0
      at(5) = 3.30E-11
      at(6) = -0.30
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,97) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,97) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,97) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R98 OH + NO2 = HONO2 + m   0.40 3.30E-30 -3.00   0.0 4.10E-11  0.00      0.0
!     IUPAC no change - note 10

      at(1) = 0.40
      at(2) = 3.30E-30
      at(3) = -3.00
      at(4) = 0.0
      at(5) = 4.10E-11
      at(6) = 0.00
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,98) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,98) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,98) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R99 OH + OH = H2O2 + m      0.50 6.90E-31 -0.80   0.0 2.60E-11  0.00      0.0
!     IUPAC no change

      at(1) = 0.50
      at(2) = 6.90E-31
      at(3) = -0.80
      at(4) = 0.0
      at(5) = 2.60E-11
      at(6) = 0.00
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))
      
      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,99) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,99) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)    = zo(i) / zi(i)
          arg2(i)  = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,99) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R100 EtCO3 + NO2 = PPAN + m    0.30 2.70E-28 -7.10   0.0 1.20E-11 -0.90      0.0
!     MCM v3.1 no change - note 12

      at(1) = 0.30
      at(2) = 2.70E-28
      at(3) = -7.10
      at(4) = 0.0
      at(5) = 1.20E-11
      at(6) = -0.90
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,100) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,100) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)     = zo(i) / zi(i)
          arg2(i)   = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,100) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R101 PPAN + m = EtCO3 + NO2  0.36 1.70E-03  0.00  11280.0 8.30E+16  0.00  13940.0
!     IUPAC no change

      at(1) = 0.36
      at(2) = 1.70E-03
      at(3) = 0.00
      at(4) = 11280.0
      at(5) = 8.30E+16
      at(6) = 0.00
      at(7) = 13940.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF( zo(i) < peps ) THEN
          RC(i,101) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,101) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)     = zo(i) / zi(i)
          arg2(i)   = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,101) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

!     R102 O(3P)S + O2 = O3S + m   0.00 5.70E-34 -2.60   0.0 0.00E+00  0.00      0.0
!     IUPAC no change - note 8

      at(1) = 0.00
      at(2) = 5.70E-34
      at(3) = -2.60
      at(4) = 0.0
      at(5) = 0.00E+00
      at(6) = 0.00
      at(7) = 0.0

      zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
      zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

      DO i=1,n_pnts
        IF ( zo(i) < peps ) THEN
          RC(i,102) = zi(i)
        ELSE IF( zi(i) < peps ) THEN
          RC(i,102) = zo(i)
        ELSE
          IF( at(1) <= 1.0 ) THEN
            zfc(i) = at(1)
          ELSE
            zfc(i) = exp( -tc(i)/at(1) )
          END IF
          zr(i)     = zo(i) / zi(i)
          arg2(i)   = 1.0/(1.0+((alog10(zr(i)))**2))
          RC(i,102) = (zo(i)/(1.0+zr(i)))*exp(arg2(i)*alog(zfc(i)))
        END IF
      ENDDO

      IF(L_ukca_aerchem) THEN

!       R103  SO2      OH     m    H2SO4   Rate data from NASA/JPL.

        at(1) = 0.60
        at(2) = 3.00E-31
        at(3) = -3.30
        at(4) = 0.0
        at(5) = 1.50E-12
        at(6) = 0.00
        at(7) = 0.0

        zo(:) = at(2)*exp(at(3)*alog(t300(:)))*exp(-at(4)/tc(:))*m(:)
        zi(:) = at(5)*exp(at(6)*alog(t300(:)))*exp(-at(7)/tc(:))

        DO i=1,n_pnts
          IF(zo(i) < peps ) THEN
             RC(i,103) = zi(i)
          ELSE IF( zi(i) < peps ) THEN
            RC(i,103) = zo(i)
          ELSE
            IF( at(1) <= 1.0 ) THEN
               zfc(i) = at(1)
            ELSE
              zfc(i) = EXP( -tc(i)/at(1) )
            END IF
            zr(i) = zo(i) / zi(i)
            arg2(i) = 1.0/(1.0+((ALOG10(zr(i)))**2))
            rc(i,103) = (zo(i)/(1.0+zr(i)))*EXP(arg2(i)*ALOG(zfc(i)))
            rc(i,103) = RC(i,103) * (1.0 - (faq(i,48) * fcloud(i)))
          END IF
        ENDDO

!       R104:  HSO3(-)(aq) + H2O2(aq) =>  H2SO4
!       Rate data from Eqn. 7 of Bower et al (1991, Atm.Env. 25A, 2401-2418)
!       Adjustments to allow for use of aqueous chemistry in model with
!       gas-phase reactions from Bergel et al (2004, JGR 109).

        WHERE(clw(:) < 1E-10)
          rc(:,104) = 0.0
        ELSEWHERE
          rc(:,104) = 1.10E+12 * (EXP(-3655.3 / tc(:))) *               &
                      (H_plus/(0.1 + H_plus))*fcloud(:)*fdiss(:,48,2)*  &
                      fdiss(:,13,1) * 1000.0 / (avc * vr(:))
        ENDWHERE

!       DMS Oxidation from IUPAC March 2007.
!       R105 : DMS + OH = DMSO2

        rc(:,105) = 1.12E-11*EXP(-250.0/TC(:))

!       R106 : DMS + OH = DMSO2

        rc(:,106) = 9.3E-39*o2(:)*EXP(5270.0/tc(:))/                  &
                  (1.0 + 7.4E-29*o2(:)*EXP(5610.0/tc(:)))

!       R107 : DMS + NO3 = DMSO2 + HNO3

        rc(:,107) = 1.9E-13*EXP(520.0/tc(:))

!       R108, R109, R110 : parameterised DMS

!       R111 : NH3 + OH = NH2 + H2O  (NH2 + NO = N2 + H2O)
        rc(:,111) = 3.5E-12*EXP(-925.0/tc(:))

!       R112 : MONOTERP + OH = SEC_ORG
        rc(:,112) = 1.2E-11*EXP(444.0/tc(:))            ! MM (2008)

!       R113 : MONOTERP + O3 = SEC_ORG
        rc(:,113) = 1.01E-15*EXP(-732.0/tc(:))          ! MM (2008)

!       R114 : MONOTERP + NO3 = SEC_ORG
        rc(:,114) = 1.19E-12*EXP(490.0/tc(:))           ! MM (2008)

!       R115 : SO3 + O3 (Aq) =  SO4                     !
        WHERE(clw(:) < 1E-10)
          rc(:,115) = 0.0
        ELSEWHERE
!          rc(:,115) = ((4.28E+13*(EXP(-5532.8/tc(:)))*fdiss(:,48,2))+   &
! Only so3 + o3 for now.
           rc(:,115) = (7.43E+16*(EXP(-5280.0/tc(:)))*fdiss(:,48,3))*   &
                      fcloud(:)*fdiss(:,3,1)*1000.0/(avc*vr(:))
        ENDWHERE

!       These are not used directly, but to calculate product ratios

!       CH3SO2 => CH3 + SO2
        k_dms(:,1) = 100.0                       ! Karl et al 2007

!       CH3SO2 + O3 => CH3SO3
        k_dms(:,2) = 6.3E-13                     ! Karl et al 2007

!       CH3SO2 + NO2 => CH3SO3
        k_dms(:,3) = 2.2E-11                     ! Karl et al 2007

!       CH3SO3 + HO2 => MSA
        k_dms(:,4) = 5.0E-11                     ! Karl et al 2007

!       CH3SO3 => CH3 + SO3
        k_dms(:,5) = 1.2e-3                      ! Karl et al 2007

      END IF    ! L_ukca_aerchem


      IF (lhook) CALL dr_hook('UKCA_CHEMCO',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CHEMCO
