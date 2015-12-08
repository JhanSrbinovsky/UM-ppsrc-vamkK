! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Module defining ASAD arrays, variables, and parameters
!
!  Part of the UKCA model, a community model supported by the
!  Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!   Language:  FORTRAN 90 
!   This code is written to UMDP3 v6 programming standards.
!
! ----------------------------------------------------------------------

MODULE ASAD_MOD

USE ukca_option_mod, ONLY: jpctr, jpspec, jpdw, jpdd, jpnr, jpbk, jptk,  &
                           jphk, jppj, L_ukca_strattrop, L_ukca_strat,   &
                           L_ukca_stratcfc
USE ukca_chem_schemes_mod, ONLY: int_method_nr
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE Control_Max_Sizes
IMPLICIT NONE
SAVE
PUBLIC

! Interface section
INTERFACE asad_mod_init
  MODULE PROCEDURE asad_mod_init
END INTERFACE asad_mod_init

INTERFACE asad_mod_final
  MODULE PROCEDURE asad_mod_final
END INTERFACE asad_mod_final

REAL, ALLOCATABLE :: wp(:)           ! water vapour field
REAL, ALLOCATABLE :: dpd(:,:)
REAL, ALLOCATABLE :: dpw(:,:)
REAL, ALLOCATABLE :: emr(:,:)
REAL, ALLOCATABLE :: fj(:,:,:)       ! Full jacobian
REAL, ALLOCATABLE :: qa(:,:)
REAL, ALLOCATABLE :: ratio(:,:)
REAL, ALLOCATABLE :: p(:)
REAL, ALLOCATABLE :: t(:)
REAL, ALLOCATABLE :: t300(:)
REAL, ALLOCATABLE :: tnd(:)          ! total number density
REAL, ALLOCATABLE :: pmintnd(:)
REAL, ALLOCATABLE :: f(:,:)          ! tracer concentrations
REAL, ALLOCATABLE :: fdot(:,:)
REAL, TARGET, ALLOCATABLE :: pd(:,:)
REAL, POINTER     :: prod(:,:)
REAL, POINTER     :: slos(:,:)
REAL, ALLOCATABLE :: y(:,:)
REAL, ALLOCATABLE :: ydot(:,:)
REAL, ALLOCATABLE :: ftilde(:,:)     ! lower order solution 
REAL, ALLOCATABLE :: f_no(:,:)       ! partioning factors
REAL, ALLOCATABLE :: f_no2(:,:)      !    between
REAL, ALLOCATABLE :: f_no3(:,:)      ! NO, NO2, and NO3
REAL, ALLOCATABLE :: ej(:,:)
REAL, ALLOCATABLE :: rk(:,:)
REAL, ALLOCATABLE :: prk(:,:)
REAL, ALLOCATABLE :: deriv(:,:,:)
REAL, ALLOCATABLE :: za(:)           ! Aerosol surface area
REAL, ALLOCATABLE :: co3(:)          ! Column ozone
REAL, ALLOCATABLE :: lati(:)         ! Latitude
REAL, ALLOCATABLE :: sphno3(:)       ! Amount of HNO3 in solid phase
REAL, ALLOCATABLE :: sph2o(:)        ! Amount of H2O in solid phase
REAL, ALLOCATABLE :: depvel(:,:,:)
REAL, ALLOCATABLE :: k298(:)
REAL, ALLOCATABLE :: dhr(:)
REAL, ALLOCATABLE :: kd298(:,:)
REAL, ALLOCATABLE :: ddhr(:,:)
REAL, ALLOCATABLE :: ab(:,:)
REAL, ALLOCATABLE :: at(:,:)
REAL, ALLOCATABLE :: aj(:,:)
REAL, ALLOCATABLE :: ah(:,:)
REAL, ALLOCATABLE :: ztabpd(:,:)
REAL, ALLOCATABLE :: shno3(:)        ! No. density type 1 psc solid phase hno3
REAL, ALLOCATABLE :: sh2o(:)         ! No. density type 2 psc solid phase h2o
REAL, ALLOCATABLE :: fpsc1(:)        ! 1.0 if type 1 psc's are present, else 0
REAL, ALLOCATABLE :: fpsc2(:)        ! 1.0 if type 2 psc's are present, else 0
REAL, ALLOCATABLE :: spfj(:,:)       ! Sparse full Jacobian

INTEGER, ALLOCATABLE :: madvtr(:)
INTEGER, ALLOCATABLE :: majors(:)
INTEGER, ALLOCATABLE :: moffam(:)
INTEGER, ALLOCATABLE :: nodd(:)
INTEGER, ALLOCATABLE :: nltrf(:)
INTEGER, ALLOCATABLE :: nltr3(:)
INTEGER, ALLOCATABLE :: ipa(:,:)     ! Pivot information for solving jacobian
INTEGER, ALLOCATABLE :: ipa2(:)
INTEGER, ALLOCATABLE :: nltrim(:,:)
INTEGER, ALLOCATABLE :: nlpdv(:,:)
INTEGER, ALLOCATABLE :: nfrpx(:)     ! index to fractional product array nfrpx:
! one entry for each reaction. If zero, there are no fractional products. If nonzero,
! contains the array element in frpx for the first coefficient for that reaction.
INTEGER, ALLOCATABLE :: ntabfp(:,:)  ! Table used for indexing the fractional
!    products:  ntabfp(i,1) contains the species no.,
!               ntabfp(i,2) contains the reaction no., and
!               ntabfp(i,3) contains the array location in the frpx array.
INTEGER, ALLOCATABLE :: ntabpd(:,:)
INTEGER, ALLOCATABLE :: npdfr(:,:)
INTEGER, ALLOCATABLE :: ngrp(:,:)
INTEGER, ALLOCATABLE :: njcgrp(:,:)
INTEGER, ALLOCATABLE :: nprdx3(:,:,:)
INTEGER, ALLOCATABLE :: nprdx2(:,:)
INTEGER, ALLOCATABLE :: nprdx1(:)
INTEGER, ALLOCATABLE :: njacx3(:,:,:)
INTEGER, ALLOCATABLE :: njacx2(:,:)
INTEGER, ALLOCATABLE :: njacx1(:)
INTEGER, ALLOCATABLE :: nmpjac(:)
INTEGER, ALLOCATABLE :: npjac1(:,:)
INTEGER, ALLOCATABLE :: nbrkx(:)
INTEGER, ALLOCATABLE :: ntrkx(:)
INTEGER, ALLOCATABLE :: nprkx(:)
INTEGER, ALLOCATABLE :: nhrkx(:)
INTEGER, ALLOCATABLE :: nlall(:)
INTEGER, ALLOCATABLE :: nlstst(:)
INTEGER, ALLOCATABLE :: nlf(:)
INTEGER, ALLOCATABLE :: nlmajmin(:)
INTEGER, ALLOCATABLE :: nldepd(:)
INTEGER, ALLOCATABLE :: nldepw(:)
INTEGER, ALLOCATABLE :: nlemit(:)
INTEGER, ALLOCATABLE :: nldepx(:)
INTEGER, ALLOCATABLE :: njcoth(:,:)
INTEGER, ALLOCATABLE :: nmzjac(:)
INTEGER, ALLOCATABLE :: nzjac1(:,:)
INTEGER, ALLOCATABLE :: njcoss(:,:)
INTEGER, ALLOCATABLE :: nmsjac(:)
INTEGER, ALLOCATABLE :: nsjac1(:,:)
INTEGER, ALLOCATABLE :: nsspt(:)
INTEGER, ALLOCATABLE :: nspi(:,:)
INTEGER, ALLOCATABLE :: nsspi(:,:)
INTEGER, ALLOCATABLE :: nssi(:)
INTEGER, ALLOCATABLE :: nssrt(:)
INTEGER, ALLOCATABLE :: nssri(:,:)
INTEGER, ALLOCATABLE :: nssrx(:,:)
REAL, ALLOCATABLE :: frpb(:)         ! fractional product array (bimol)
REAL, ALLOCATABLE :: frpt(:)         ! fractional product array (trimol)
REAL, ALLOCATABLE :: frpj(:)         ! fractional product array (phot)
REAL, ALLOCATABLE :: frph(:)         ! fractional product array (het)
REAL, ALLOCATABLE :: frpx(:)         ! fractional product array (total)
! sparse algebra
INTEGER, ALLOCATABLE :: pointer(:,:)  ! Map of nonzero entries
INTEGER, ALLOCATABLE :: pointer1(:,:) ! modified map (after decomposition)
INTEGER, ALLOCATABLE :: pointer2(:,:) ! Map of nonzero entries, before reordering
INTEGER, ALLOCATABLE :: ro(:)         ! reordering of tracers to minimize fill-in

INTEGER, ALLOCATABLE :: ilcf(:)
INTEGER, ALLOCATABLE :: ilss(:)
INTEGER, ALLOCATABLE :: ilct(:)
INTEGER, ALLOCATABLE :: ilftr(:)
INTEGER, ALLOCATABLE :: ilft(:)
INTEGER, ALLOCATABLE :: ilstmin(:)

LOGICAL, ALLOCATABLE :: linfam(:,:)
LOGICAL, ALLOCATABLE :: ldepd(:)     ! T for dry deposition
LOGICAL, ALLOCATABLE :: ldepw(:)     ! T for wet deposition
LOGICAL, ALLOCATABLE :: lemit(:)     ! T for emission

CHARACTER(len=10), ALLOCATABLE :: advt(:)       ! advected tracers
CHARACTER(len=10), ALLOCATABLE :: nadvt(:)      ! non-advected species 
CHARACTER(len=10), ALLOCATABLE :: family(:)     ! family
CHARACTER(len=10), ALLOCATABLE :: speci(:)      ! species names
CHARACTER(len=2),  ALLOCATABLE :: ctype(:)      ! species type
CHARACTER(len=10), ALLOCATABLE :: spb(:,:)      ! species from bimolecular rates
CHARACTER(len=10), ALLOCATABLE :: spt(:,:)      ! species from termolecular rates
CHARACTER(len=10), ALLOCATABLE :: spj(:,:)      ! species from photolysis rates
CHARACTER(len=10), ALLOCATABLE :: sph(:,:)      ! species from heterogenous rates

REAL, PARAMETER    :: pmin = 1.0e-20
REAL, PARAMETER    :: ptol = 1.0e-5            ! tolerance for time integration
REAL, PARAMETER    :: ftol = 1.0e-3            ! tolerance in family member iteration

INTEGER, PARAMETER :: kfphot=0
INTEGER, PARAMETER :: jpss = 16
INTEGER, PARAMETER :: jpssr = 51
INTEGER, PARAMETER :: nvar = 17
INTEGER, PARAMETER :: nllv = 17
INTEGER, PARAMETER :: ninv = 23
INTEGER, PARAMETER :: nout=71
!     nout:        Fortran channel for output in subroutine OUTVMR
INTEGER, PARAMETER :: jpem = 9      ! IS THIS A CONSTANT/USED?
INTEGER, PARAMETER :: jpeq = 2      ! dimension for dissociation arrays 
INTEGER, PARAMETER :: jddept = 6
!     jddept:      Number of time periods used in dry deposition i.e.
!                  summer(day,night,24h ave), winter(day,night,24h ave)
INTEGER, PARAMETER :: jddepc = 5
!     jddepc:      Number of land use categories used in dry dep.
INTEGER, PARAMETER :: jpdwio = 56
!     jpdwio       Fortran i/o unit to read/write anything to do with
!                  wet/dry deposition
INTEGER, PARAMETER :: jpemio = 57
!     jpemio       Fortran i/o unit to read in emissions
INTEGER, PARAMETER :: jpfrpd=100
INTEGER, PARAMETER :: jpab    = 3
INTEGER, PARAMETER :: jpat    = 7
INTEGER, PARAMETER :: jpaj    = 3
INTEGER, PARAMETER :: jpah    = 3
INTEGER, PARAMETER :: jpspb   = 6
INTEGER, PARAMETER :: jpspt   = 4
INTEGER, PARAMETER :: jpspj   = 6
INTEGER, PARAMETER :: jpsph   = 6
INTEGER, PARAMETER :: jpmsp   = jpspb
INTEGER, PARAMETER :: jppjac  = 10
INTEGER, PARAMETER :: jpkargs = 10
INTEGER, PARAMETER :: jprargs = 10
INTEGER, PARAMETER :: jpcargs = 1


! Fractional product parameters - these are initialsed to jp values in asad_mod_init
INTEGER :: jpfrpb
INTEGER :: jpfrpt
INTEGER :: jpfrpj
INTEGER :: jpfrph
INTEGER :: jpfrpx

INTEGER, PARAMETER :: jpcio   = 55
INTEGER, PARAMETER :: spfjsize_max =1000 ! maximum number of
                                         ! nonzero matrix elements


CHARACTER(LEN=2),  PARAMETER :: jpfm = 'FM'     ! Family member
CHARACTER(LEN=2),  PARAMETER :: jpif = 'FT'     ! Family member depending in timestep
CHARACTER(LEN=2),  PARAMETER :: jpsp = 'TR'     ! Independent tracer
CHARACTER(LEN=2),  PARAMETER :: jpna = 'SS'     ! Steady-state species
CHARACTER(LEN=2),  PARAMETER :: jpco = 'CT'     ! Constant
CHARACTER(LEN=2),  PARAMETER :: jpcf = 'CF'     ! Constant with spatial field

LOGICAL, PARAMETER :: lvmr=.true.    ! T for volume mixing ratio
LOGICAL      :: o1d_in_ss      ! T for steady state,
LOGICAL      :: o3p_in_ss      ! these are set in routine:
LOGICAL      :: n_in_ss        ! asad_mod_init
LOGICAL      :: h_in_ss        ! 
INTEGER, PARAMETER :: nss_o1d=1      ! indicies of deriv array
INTEGER, PARAMETER :: nss_o3p=2      !    "         "
INTEGER, PARAMETER :: nss_n=3        !    "         "
INTEGER, PARAMETER :: nss_h=4        !    "         "

REAL :: fch4,fco2,fh2,fn2,fo2
REAL    :: cdt                       ! chemistry timestep
REAL    :: peps                      ! 
REAL    :: dtime                     ! model timestep
REAL, PARAMETER :: tslimit = 1200.0  ! timestep limit for some solvers

INTEGER, PARAMETER :: kcdt = 3600    ! timestep for N-R solver !! MAY NOT WANT HERE
INTEGER :: nrsteps
INTEGER :: nitnr          ! Iterations in ftoy for IMPACT solver
INTEGER :: nitfg          ! Max no of iterations in ftoy
INTEGER :: ntrf           ! Counter for tracers
INTEGER :: ntr3
INTEGER :: nnaf           ! Counter for non-advected tracers
INTEGER :: nuni
INTEGER :: nsst           ! No of steady-state species
INTEGER :: ncsteps        ! No of chemical steps
INTEGER :: jlp            ! level to be integrated
INTEGER :: nit0=20        ! ftoy iterations with method=0
INTEGER :: nfphot
INTEGER :: jsubs
INTEGER :: jlst
INTEGER :: method         ! chemistry integration method
INTEGER :: interval       ! interval in timesteps between calls to chemistry
INTEGER :: nnfrp          ! Total number of fractional products
INTEGER :: nfrpd
INTEGER :: nstst          ! No of steady state species
INTEGER :: nf
INTEGER :: ndepd          ! No of dry deposited species
INTEGER :: ndepw          ! No of wet deposited species
INTEGER :: nemit          ! No of emitted species
INTEGER :: ntro3, ntroh, ntrho2, ntrno 
INTEGER :: nspo1d, nspo3p, nspo3, nspoh
INTEGER :: nspho2, nspno, nspn, nsph
INTEGER :: ih_o3, ih_h2o2, ih_so2, ih_hno3  ! index for soluble species
INTEGER :: ih_ho2, ih_n2o5                  ! index for species involved in het. chem
INTEGER :: ihso3_h2o2                       ! Index for HSO3- + H2O2(aq) reaction
INTEGER :: ihso3_o3                         ! Index for HSO3- + O3(aq) reaction
INTEGER :: iso3_o3                          ! Index for SO3-- + O3(aq) reaction
INTEGER :: ih2so4_hv                        ! Index for H2SO4 + hv reaction
INTEGER :: iso2_oh                          ! Index for SO2 + OH reaction
INTEGER :: ih2o2_oh                         ! Index for H2O2 + OH reaction
INTEGER :: ihno3_oh                         ! Index for HNO3 + OH reaction
INTEGER :: in2o5_h                          ! Index for N2O5 => HONO2 heterog. reaction
INTEGER :: iho2_h                           ! Index for HO2 + HO2 => H2O2 heterog. "

LOGICAL :: lsvjac         ! Flag for saving jacobian if not recalculated
LOGICAL :: ljacx

! ltrig set to debug slow convergence systems
! shared between asad_spimpmjp and asad_spmjpdriv
LOGICAL :: ltrig

CONTAINS

! ######################################################################
SUBROUTINE asad_mod_init

! To allocate and initialise ASAD arrays and variables

USE UKCA_CHEM_DEFS_MOD,       ONLY:  chch_t, chch_defs
USE UM_ParVars
USE Control_Max_Sizes
USE Ereport_mod,              ONLY: ereport
USE ukca_option_mod,          ONLY:l_ukca, ukca_int_method
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

INTEGER :: i
INTEGER :: errcode
CHARACTER(LEN=72) :: cmessage

LOGICAL, SAVE :: firstcall=.true.
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('ASAD_MOD:ASAD_MOD_INIT',zhook_in,zhook_handle)
! nullify prod and slos on firstcall to give DISSASSOCIATED attribute
IF (firstcall) THEN
   nullify(prod)
   nullify(slos)
   firstcall = .false.
end IF



! Fractional product parameters - set total number allowed to be
! total number of reactions * 2 as each fractional product has 4 potentials!!
jpfrpb  = (jpspb-2)*jpbk
jpfrpt  = (jpspt-2)*jptk   
jpfrpj  = (jpspj-2)*jppj   
jpfrph  = (jpsph-2)*jphk   
jpfrpx  = jpfrpb + jpfrpt + jpfrpj + jpfrph  ! Total FPs



! pd is a TARGET
IF (.NOT. ALLOCATED(pd)) ALLOCATE(pd(theta_field_size,2*jpspec))
! prod and slos are pointers
IF (.NOT. ASSOCIATED(prod)) ALLOCATE(prod(theta_field_size,jpspec))
IF (.NOT. ASSOCIATED(slos)) ALLOCATE(slos(theta_field_size,jpspec))


IF (.NOT. ALLOCATED(wp)) ALLOCATE(wp(theta_field_size))
IF (.NOT. ALLOCATED(linfam)) ALLOCATE(linfam(theta_field_size,0:jpctr))
IF (.NOT. ALLOCATED(madvtr)) ALLOCATE(madvtr(jpspec))
IF (.NOT. ALLOCATED(majors)) ALLOCATE(majors(jpctr))
IF (.NOT. ALLOCATED(moffam)) ALLOCATE(moffam(jpspec))
IF (.NOT. ALLOCATED(nodd)) ALLOCATE(nodd(jpspec))
IF (.NOT. ALLOCATED(nltrf)) ALLOCATE(nltrf(jpctr))
IF (.NOT. ALLOCATED(nltr3)) ALLOCATE(nltr3(jpctr))
IF (.NOT. ALLOCATED(advt)) ALLOCATE(advt(jpctr))
IF (.NOT. ALLOCATED(nadvt)) ALLOCATE(nadvt(jpspec-jpctr))
IF (.NOT. ALLOCATED(family)) ALLOCATE(family(jpspec))
IF (.NOT. ALLOCATED(speci)) ALLOCATE(speci(jpspec))
IF (.NOT. ALLOCATED(ctype)) ALLOCATE(ctype(jpspec))
IF (.NOT. ALLOCATED(dpd)) ALLOCATE(dpd(theta_field_size,jpspec))
IF (.NOT. ALLOCATED(dpw)) ALLOCATE(dpw(theta_field_size,jpspec))
IF (.NOT. ALLOCATED(emr)) ALLOCATE(emr(theta_field_size,jpspec))
IF (.NOT. ALLOCATED(ipa)) ALLOCATE(ipa(theta_field_size,jpctr))
IF (.NOT. ALLOCATED(ipa2)) ALLOCATE(ipa2(jpctr))
IF (.NOT. ALLOCATED(fj)) ALLOCATE(fj(theta_field_size,jpctr,jpctr))
IF (.NOT. ALLOCATED(qa)) ALLOCATE(qa(theta_field_size,jpspec))
IF (.NOT. ALLOCATED(ratio)) ALLOCATE(ratio(theta_field_size,jpspec))
IF (.NOT. ALLOCATED(p)) ALLOCATE(p(theta_field_size))
IF (.NOT. ALLOCATED(t)) ALLOCATE(t(theta_field_size))
IF (.NOT. ALLOCATED(t300)) ALLOCATE(t300(theta_field_size))
IF (.NOT. ALLOCATED(tnd)) ALLOCATE(tnd(theta_field_size))
IF (.NOT. ALLOCATED(pmintnd)) ALLOCATE(pmintnd(theta_field_size))
IF (.NOT. ALLOCATED(f)) ALLOCATE(f(theta_field_size,jpctr))
IF (.NOT. ALLOCATED(fdot)) ALLOCATE(fdot(theta_field_size,jpctr))
IF (.NOT. ALLOCATED(y)) ALLOCATE(y(theta_field_size,jpspec))
IF (.NOT. ALLOCATED(ydot)) ALLOCATE(ydot(theta_field_size,jpspec))
IF (.NOT. ALLOCATED(ftilde)) ALLOCATE(ftilde(theta_field_size, jpctr))
IF (.NOT. ALLOCATED(f_no)) ALLOCATE(f_no(theta_field_size,model_levels))
IF (.NOT. ALLOCATED(f_no2)) ALLOCATE(f_no2(theta_field_size,model_levels))
IF (.NOT. ALLOCATED(f_no3)) ALLOCATE(f_no3(theta_field_size,model_levels))          
IF (.NOT. ALLOCATED(ej)) ALLOCATE(ej(theta_field_size,jpctr))
IF (.NOT. ALLOCATED(rk)) ALLOCATE(rk(theta_field_size,jpnr))
IF (.NOT. ALLOCATED(prk)) ALLOCATE(prk(theta_field_size,jpnr))
IF (.NOT. ALLOCATED(deriv)) ALLOCATE(deriv(theta_field_size,4,4))
IF (.NOT. ALLOCATED(nspi)) ALLOCATE(nspi(jpnr,jpmsp))
IF (.NOT. ALLOCATED(nsspt)) ALLOCATE(nsspt(jpss))
IF (.NOT. ALLOCATED(nsspi)) ALLOCATE(nsspi(jpss,jpssr))
IF (.NOT. ALLOCATED(nssi)) ALLOCATE(nssi(jpss))
IF (.NOT. ALLOCATED(nssrt)) ALLOCATE(nssrt(jpss))
IF (.NOT. ALLOCATED(nssri)) ALLOCATE(nssri(jpss,jpssr))
IF (.NOT. ALLOCATED(nssrx)) ALLOCATE(nssrx(jpss,jpssr)) 
IF (.NOT. ALLOCATED(za)) ALLOCATE(za(theta_field_size))
IF (.NOT. ALLOCATED(co3)) ALLOCATE(co3(theta_field_size))
IF (.NOT. ALLOCATED(lati)) ALLOCATE(lati(theta_field_size))
IF (.NOT. ALLOCATED(sphno3)) ALLOCATE(sphno3(theta_field_size))
IF (.NOT. ALLOCATED(sph2o)) ALLOCATE(sph2o(theta_field_size))
IF (.NOT. ALLOCATED(depvel)) ALLOCATE(depvel(jddept,jddepc,jpdd))
IF (.NOT. ALLOCATED(k298)) ALLOCATE(k298(jpdw))
IF (.NOT. ALLOCATED(dhr)) ALLOCATE(dhr(jpdw))
IF (.NOT. ALLOCATED(kd298)) ALLOCATE(kd298(jpdw,jpeq))
IF (.NOT. ALLOCATED(ddhr)) ALLOCATE(ddhr(jpdw,jpeq))
IF (.NOT. ALLOCATED(ldepd)) ALLOCATE(ldepd(jpspec))
IF (.NOT. ALLOCATED(ldepw)) ALLOCATE(ldepw(jpspec))
IF (.NOT. ALLOCATED(lemit)) ALLOCATE(lemit(jpspec))
IF (.NOT. ALLOCATED(nltrim)) ALLOCATE(nltrim(0:jpctr,3))
IF (.NOT. ALLOCATED(nlpdv)) ALLOCATE(nlpdv((jpspj-2)*jppj,2))
IF (.NOT. ALLOCATED(ab)) ALLOCATE(ab(jpbk+1,jpab))
IF (.NOT. ALLOCATED(at)) ALLOCATE(at(jptk+1,jpat))
IF (.NOT. ALLOCATED(aj)) ALLOCATE(aj(jppj+1,jpaj))
IF (.NOT. ALLOCATED(ah)) ALLOCATE(ah(jphk+1,jpah))
IF (.NOT. ALLOCATED(spb)) ALLOCATE(spb(jpbk+1,jpspb))
IF (.NOT. ALLOCATED(spt)) ALLOCATE(spt(jptk+1,jpspt))
IF (.NOT. ALLOCATED(spj)) ALLOCATE(spj(jppj+1,jpspj))
IF (.NOT. ALLOCATED(sph)) ALLOCATE(sph(jphk+1,jpsph))
IF (.NOT. ALLOCATED(frpb)) ALLOCATE(frpb(jpfrpb))
IF (.NOT. ALLOCATED(frpt)) ALLOCATE(frpt(jpfrpt))
IF (.NOT. ALLOCATED(frpj)) ALLOCATE(frpj(jpfrpj))
IF (.NOT. ALLOCATED(frph)) ALLOCATE(frph(jpfrph))
IF (.NOT. ALLOCATED(frpx)) ALLOCATE(frpx(jpfrpx))
IF (.NOT. ALLOCATED(ztabpd)) ALLOCATE(ztabpd(jpfrpd,2))
IF (.NOT. ALLOCATED(nfrpx)) ALLOCATE(nfrpx(jpnr))
IF (.NOT. ALLOCATED(ntabfp)) ALLOCATE(ntabfp(jpfrpx,3))
IF (.NOT. ALLOCATED(ntabpd)) ALLOCATE(ntabpd(jpfrpd,3))
IF (.NOT. ALLOCATED(npdfr)) ALLOCATE(npdfr(jpnr,2))
IF (.NOT. ALLOCATED(ngrp)) ALLOCATE(ngrp(2*jpspec,3))
IF (.NOT. ALLOCATED(njcgrp)) ALLOCATE(njcgrp(jpctr,3))
IF (.NOT. ALLOCATED(nprdx3)) ALLOCATE(nprdx3(3,(jpnr/(3*3))+3*3,2*jpspec))
IF (.NOT. ALLOCATED(nprdx2)) ALLOCATE(nprdx2(2,2*jpspec))
IF (.NOT. ALLOCATED(nprdx1)) ALLOCATE(nprdx1(2*jpspec))
IF (.NOT. ALLOCATED(njacx3)) ALLOCATE(njacx3(3,(jpnr/(3*3))+3*3,jpctr))
IF (.NOT. ALLOCATED(njacx2)) ALLOCATE(njacx2(2,jpctr))
IF (.NOT. ALLOCATED(njacx1)) ALLOCATE(njacx1(jpctr))
IF (.NOT. ALLOCATED(nmpjac)) ALLOCATE(nmpjac(jpctr))
IF (.NOT. ALLOCATED(npjac1)) ALLOCATE(npjac1(jppjac,jpctr))
IF (.NOT. ALLOCATED(nbrkx)) ALLOCATE(nbrkx(jpbk+1))
IF (.NOT. ALLOCATED(ntrkx)) ALLOCATE(ntrkx(jptk+1))
IF (.NOT. ALLOCATED(nprkx)) ALLOCATE(nprkx(jppj+1))
IF (.NOT. ALLOCATED(nhrkx)) ALLOCATE(nhrkx(jphk+1))
IF (.NOT. ALLOCATED(nlall)) ALLOCATE(nlall(jpspec))
IF (.NOT. ALLOCATED(nlstst)) ALLOCATE(nlstst(jpspec))
IF (.NOT. ALLOCATED(nlf)) ALLOCATE(nlf(jpspec))
IF (.NOT. ALLOCATED(nlmajmin)) ALLOCATE(nlmajmin(jpspec))
IF (.NOT. ALLOCATED(nldepd)) ALLOCATE(nldepd(jpspec))
IF (.NOT. ALLOCATED(nldepw)) ALLOCATE(nldepw(jpspec))
IF (.NOT. ALLOCATED(nlemit)) ALLOCATE(nlemit(jpspec))
IF (.NOT. ALLOCATED(nldepx)) ALLOCATE(nldepx(jpspec))
IF (.NOT. ALLOCATED(njcoth)) ALLOCATE(njcoth(jpnr,jpmsp))
IF (.NOT. ALLOCATED(nmzjac)) ALLOCATE(nmzjac(jpctr))
IF (.NOT. ALLOCATED(nzjac1)) ALLOCATE(nzjac1(jpnr,jpctr))
IF (.NOT. ALLOCATED(njcoss)) ALLOCATE(njcoss(jpnr,jpmsp))
IF (.NOT. ALLOCATED(nmsjac)) ALLOCATE(nmsjac(jpctr))
IF (.NOT. ALLOCATED(nsjac1)) ALLOCATE(nsjac1(jpnr,jpctr))
IF (.NOT. ALLOCATED(shno3)) ALLOCATE(shno3(theta_field_size))
IF (.NOT. ALLOCATED(sh2o)) ALLOCATE(sh2o(theta_field_size))
IF (.NOT. ALLOCATED(fpsc1)) ALLOCATE(fpsc1(theta_field_size))
IF (.NOT. ALLOCATED(fpsc2)) ALLOCATE(fpsc2(theta_field_size))

! the following had save attribs.
IF (.NOT. ALLOCATED(ilcf)) ALLOCATE(ilcf(jpspec))
IF (.NOT. ALLOCATED(ilss)) ALLOCATE(ilss(jpspec))
IF (.NOT. ALLOCATED(ilct)) ALLOCATE(ilct(jpspec))
IF (.NOT. ALLOCATED(ilftr)) ALLOCATE(ilftr(jpspec))
IF (.NOT. ALLOCATED(ilft)) ALLOCATE(ilft(jpspec))
IF (.NOT. ALLOCATED(ilstmin)) ALLOCATE(ilstmin(jpspec))

! EQUIVALENCE ( pd(1,1), prod(1,1) )
! EQUIVALENCE ( pd(1,jpspec+1), slos(1,1) )
prod => pd(:,1:jpspec)
slos => pd(:,jpspec+1:2*jpspec)

! Set integration method (1 = IMPACT; 3 = N-R solver; 5 = Backward-Euler)
method = ukca_int_method

! Initialize variables that may be changed in cinit
nrsteps = 45            ! No of N-R steps, set > 50 to debug convergence failures
nitnr   = 10            ! Iterations in ftoy
nitfg   = 10            ! Max number of iterations in ftoy

! Initialise arrays
njcoth(:,:) = 0
dtime   = 1200.0         ! model timestep, reset in chemistry_ctl

deriv(:,:,:) = 1.0     ! Temp fix for deriv being uninitialised in first
                       ! solver iteration 

IF (method == int_method_NR) THEN
  IF (.NOT. ALLOCATED(pointer))  ALLOCATE(pointer(jpctr, jpctr))
  IF (.NOT. ALLOCATED(pointer1))  ALLOCATE(pointer1(jpctr, jpctr))
  IF (.NOT. ALLOCATED(pointer2))  ALLOCATE(pointer2(jpctr, jpctr))
  IF (.NOT. ALLOCATED(ro))  ALLOCATE(ro(jpctr))
  IF (.NOT. ALLOCATED(spfj))  ALLOCATE(spfj(theta_field_size,spfjsize_max))
END IF

! Find out which species are in steady state (for N-R solver)
o1d_in_ss = .FALSE.
o3p_in_ss = .FALSE.
n_in_ss   = .FALSE.
h_in_ss   = .FALSE.
DO i=1,jpspec
 IF (chch_defs(i)%speci=='O(1D)     ' .AND.                             &
     chch_defs(i)%ctype(1:2)==jpna) o1d_in_ss=.TRUE.
 IF (chch_defs(i)%speci=='O(3P)     ' .AND.                             &
     chch_defs(i)%ctype(1:2)==jpna) o3p_in_ss=.TRUE.
 IF (chch_defs(i)%speci=='N         ' .AND.                             &
     chch_defs(i)%ctype(1:2)==jpna) n_in_ss=.TRUE.
 IF (chch_defs(i)%speci=='H         ' .AND.                             &
     chch_defs(i)%ctype(1:2)==jpna) h_in_ss=.TRUE.
ENDDO

! Currently N-R solver assumes that O(1D) is in steady-state
IF (.NOT. o1d_in_ss .AND. (L_ukca_strat .OR. &   ! L_ukca_wachem .OR.  &
    L_ukca_stratcfc .OR. L_ukca_strattrop)) THEN
  cmessage=' O(1D) is not a Steady-State species'
  errcode = 1
  CALL EREPORT('ASAD_MOD_INIT',errcode,cmessage)
ENDIF

IF (lhook) CALL dr_hook('ASAD_MOD:ASAD_MOD_INIT',zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_mod_init

! ######################################################################
SUBROUTINE asad_mod_final

! To deallocate ASAD arrays

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('ASAD_MOD:ASAD_MOD_FINAL',zhook_in,zhook_handle)
! pd is a TARGET
IF (ALLOCATED(pd)) DEALLOCATE(pd)
! prod and slos are POINTERs
IF (ASSOCIATED(prod)) DEALLOCATE(prod)
IF (ASSOCIATED(slos)) DEALLOCATE(slos)

IF (ALLOCATED(nltrf)) DEALLOCATE(nltrf)
IF (ALLOCATED(nltr3)) DEALLOCATE(nltr3)
IF (ALLOCATED(advt)) DEALLOCATE(advt)
IF (ALLOCATED(nadvt)) DEALLOCATE(nadvt)
IF (ALLOCATED(family)) DEALLOCATE(family)
IF (ALLOCATED(speci)) DEALLOCATE(speci)
IF (ALLOCATED(ctype)) DEALLOCATE(ctype)
IF (ALLOCATED(dpd)) DEALLOCATE(dpd)
IF (ALLOCATED(dpw)) DEALLOCATE(dpw)
IF (ALLOCATED(emr)) DEALLOCATE(emr)
IF (ALLOCATED(ipa)) DEALLOCATE(ipa)
IF (ALLOCATED(ipa2)) DEALLOCATE(ipa2)
IF (ALLOCATED(fj)) DEALLOCATE(fj)
IF (ALLOCATED(qa)) DEALLOCATE(qa)
IF (ALLOCATED(ratio)) DEALLOCATE(ratio)
IF (ALLOCATED(p)) DEALLOCATE(p)
IF (ALLOCATED(t)) DEALLOCATE(t)
IF (ALLOCATED(t300)) DEALLOCATE(t300)
IF (ALLOCATED(tnd)) DEALLOCATE(tnd)
IF (ALLOCATED(pmintnd)) DEALLOCATE(pmintnd)
IF (ALLOCATED(f)) DEALLOCATE(f)
IF (ALLOCATED(fdot)) DEALLOCATE(fdot)
IF (ALLOCATED(y)) DEALLOCATE(y)
IF (ALLOCATED(ydot)) DEALLOCATE(ydot)
IF (ALLOCATED(ftilde)) DEALLOCATE(ftilde)
IF (ALLOCATED(f_no)) DEALLOCATE(f_no)
IF (ALLOCATED(f_no2)) DEALLOCATE(f_no2)
IF (ALLOCATED(f_no3)) DEALLOCATE(f_no3)
IF (ALLOCATED(ej)) DEALLOCATE(ej)
IF (ALLOCATED(rk)) DEALLOCATE(rk)
IF (ALLOCATED(prk)) DEALLOCATE(prk)
IF (ALLOCATED(deriv)) DEALLOCATE(deriv)
IF (ALLOCATED(nspi)) DEALLOCATE(nspi)
IF (ALLOCATED(nsspt)) DEALLOCATE(nsspt)
IF (ALLOCATED(nsspi)) DEALLOCATE(nsspi)
IF (ALLOCATED(nssi)) DEALLOCATE(nssi)
IF (ALLOCATED(nssrt)) DEALLOCATE(nssrt)
IF (ALLOCATED(nssri)) DEALLOCATE(nssri)
IF (ALLOCATED(nssrx)) DEALLOCATE(nssrx)
IF (ALLOCATED(za)) DEALLOCATE(za)
IF (ALLOCATED(co3)) DEALLOCATE(co3)
IF (ALLOCATED(lati)) DEALLOCATE(lati)
IF (ALLOCATED(sphno3)) DEALLOCATE(sphno3)
IF (ALLOCATED(sph2o)) DEALLOCATE(sph2o)
IF (ALLOCATED(depvel)) DEALLOCATE(depvel)
IF (ALLOCATED(k298)) DEALLOCATE(k298)
IF (ALLOCATED(dhr)) DEALLOCATE(dhr)
IF (ALLOCATED(kd298)) DEALLOCATE(kd298)
IF (ALLOCATED(ddhr)) DEALLOCATE(ddhr)
IF (ALLOCATED(ldepd)) DEALLOCATE(ldepd)
IF (ALLOCATED(ldepw)) DEALLOCATE(ldepw)
IF (ALLOCATED(lemit)) DEALLOCATE(lemit)
IF (ALLOCATED(nltrim)) DEALLOCATE(nltrim)
IF (ALLOCATED(nlpdv)) DEALLOCATE(nlpdv)
IF (ALLOCATED(ab)) DEALLOCATE(ab)
IF (ALLOCATED(at)) DEALLOCATE(at)
IF (ALLOCATED(aj)) DEALLOCATE(aj)
IF (ALLOCATED(ah)) DEALLOCATE(ah)
IF (ALLOCATED(spb)) DEALLOCATE(spb)
IF (ALLOCATED(spt)) DEALLOCATE(spt)
IF (ALLOCATED(spj)) DEALLOCATE(spj)
IF (ALLOCATED(sph)) DEALLOCATE(sph)
IF (ALLOCATED(frpb)) DEALLOCATE(frpb)
IF (ALLOCATED(frpt)) DEALLOCATE(frpt)
IF (ALLOCATED(frpj)) DEALLOCATE(frpj)
IF (ALLOCATED(frph)) DEALLOCATE(frph)
IF (ALLOCATED(frpx)) DEALLOCATE(frpx)
IF (ALLOCATED(ztabpd)) DEALLOCATE(ztabpd)
IF (ALLOCATED(nfrpx)) DEALLOCATE(nfrpx)
IF (ALLOCATED(ntabfp)) DEALLOCATE(ntabfp)
IF (ALLOCATED(ntabpd)) DEALLOCATE(ntabpd)
IF (ALLOCATED(npdfr)) DEALLOCATE(npdfr)
IF (ALLOCATED(ngrp)) DEALLOCATE(ngrp)
IF (ALLOCATED(njcgrp)) DEALLOCATE(njcgrp)
IF (ALLOCATED(nprdx3)) DEALLOCATE(nprdx3)
IF (ALLOCATED(nprdx2)) DEALLOCATE(nprdx2)
IF (ALLOCATED(nprdx1)) DEALLOCATE(nprdx1)
IF (ALLOCATED(njacx3)) DEALLOCATE(njacx3)
IF (ALLOCATED(njacx2)) DEALLOCATE(njacx2)
IF (ALLOCATED(njacx1)) DEALLOCATE(njacx1)
IF (ALLOCATED(nmpjac)) DEALLOCATE(nmpjac)
IF (ALLOCATED(npjac1)) DEALLOCATE(npjac1)
IF (ALLOCATED(nbrkx)) DEALLOCATE(nbrkx)
IF (ALLOCATED(ntrkx)) DEALLOCATE(ntrkx)
IF (ALLOCATED(nprkx)) DEALLOCATE(nprkx)
IF (ALLOCATED(nhrkx)) DEALLOCATE(nhrkx)
IF (ALLOCATED(nlall)) DEALLOCATE(nlall)
IF (ALLOCATED(nlstst)) DEALLOCATE(nlstst)
IF (ALLOCATED(nlf)) DEALLOCATE(nlf)
IF (ALLOCATED(nlmajmin)) DEALLOCATE(nlmajmin)
IF (ALLOCATED(nldepd)) DEALLOCATE(nldepd)
IF (ALLOCATED(nldepw)) DEALLOCATE(nldepw)
IF (ALLOCATED(nlemit)) DEALLOCATE(nlemit)
IF (ALLOCATED(nldepx)) DEALLOCATE(nldepx)
IF (ALLOCATED(njcoth)) DEALLOCATE(njcoth)
IF (ALLOCATED(nmzjac)) DEALLOCATE(nmzjac)
IF (ALLOCATED(nzjac1)) DEALLOCATE(nzjac1)
IF (ALLOCATED(njcoss)) DEALLOCATE(njcoss)
IF (ALLOCATED(nmsjac)) DEALLOCATE(nmsjac)
IF (ALLOCATED(nsjac1)) DEALLOCATE(nsjac1)
IF (ALLOCATED(shno3)) DEALLOCATE(shno3)
IF (ALLOCATED(sh2o)) DEALLOCATE(sh2o)
IF (ALLOCATED(fpsc1)) DEALLOCATE(fpsc1)
IF (ALLOCATED(fpsc2)) DEALLOCATE(fpsc2)



IF (method == int_method_NR) THEN       ! sparse vars
   IF (ALLOCATED(pointer))   DEALLOCATE(pointer)
   IF (ALLOCATED(pointer1))   DEALLOCATE(pointer1)
   IF (ALLOCATED(pointer2))   DEALLOCATE(pointer2)
   IF (ALLOCATED(ro))   DEALLOCATE(ro)
   IF (ALLOCATED(spfj))   DEALLOCATE(spfj)
ENDIF


IF (lhook) CALL dr_hook('ASAD_MOD:ASAD_MOD_FINAL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE asad_mod_final

END MODULE
