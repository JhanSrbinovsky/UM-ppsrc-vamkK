! *****************************COPYRIGHT*******************************

! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]

! *****************************COPYRIGHT*******************************

!  Description:
!   Module to contain routines concerned with stratospheric chemistry
!   Contains subroutines: relax_ozone, conserve, ukca_strat_photol,
!   and ukca_calc_ozonecol.

!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA

!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################

MODULE ukca_strat_update

USE ukca_option_mod,       ONLY: i_ukca_photol, fastjx_mode
USE ukca_photo_scheme_mod, ONLY: i_ukca_fastjx

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE um_parvars
USE control_max_sizes
USE printstatus_mod
IMPLICIT NONE

CONTAINS

!  UKCA stratospheric photolysis main routine.
!  Test version
! ######################################################################

! Subroutine Interface:

!---------------------------------------------------------------------------
! Subroutine STRAT_PHOTOL
!------------------------------------------------------------------------

! This routine computes stratospheric photolysis rates and merges the
! rates, where necessary, with the tropospheric rates. This is done for
! one level at a time. The stratospheric photolysis routines are taken
! from SLIMCAT.


SUBROUTINE ukca_strat_photol(pressure,                                         &
  temp,                                                                        &
  ozonecol,                                                                    &
  cos_zenith_angle,                                                            &
  photrates)

! SLIMCAT stratospheric photolysis routines
      USE ukca_option_mod, ONLY: jppj
USE calcjs_mod, ONLY: calcjs
USE inijtab_mod, ONLY: inijtab
USE ukca_dissoc    ! Holds photolysis rates
USE asad_mod,      ONLY: jpab, jpat, jpaj, jpah, jpspb, jpspt,                 &
  jpspj, jpsph, spj
USE um_parvars
USE control_max_sizes
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


REAL, INTENT(IN) :: pressure(row_length, rows)
REAL, INTENT(IN) :: temp(row_length, rows)
REAL, INTENT(IN) :: ozonecol(row_length, rows)
REAL, INTENT(IN) :: cos_zenith_angle(row_length, rows)

REAL, INTENT(INOUT) :: photrates(row_length, rows, jppj)

! local variables
INTEGER :: i
INTEGER :: j
INTEGER :: l
INTEGER :: k

REAL :: frac
LOGICAL, SAVE :: firstcall = .TRUE.

INTEGER, SAVE :: ih2o2=0
INTEGER, SAVE :: ihchoa=0
INTEGER, SAVE :: ihchob=0
INTEGER, SAVE :: iho2no2a=0
INTEGER, SAVE :: ihno3=0
INTEGER, SAVE :: imeooh=0
INTEGER, SAVE :: in2o5=0
INTEGER, SAVE :: ino2=0
INTEGER, SAVE :: ino3a=0
INTEGER, SAVE :: ino3b=0
INTEGER, SAVE :: io2a=0
INTEGER, SAVE :: io2b=0
INTEGER, SAVE :: io3a=0
INTEGER, SAVE :: io3b=0
INTEGER, SAVE :: io3sa=0
INTEGER, SAVE :: io3sb=0
INTEGER, SAVE :: ich4=0
INTEGER, SAVE :: ih2o=0
INTEGER, SAVE :: ino=0
INTEGER, SAVE :: in2o=0
INTEGER, SAVE :: if11=0
INTEGER, SAVE :: if12=0
INTEGER, SAVE :: iclono2a=0
INTEGER, SAVE :: iclono2b=0
INTEGER, SAVE :: ihcl=0
INTEGER, SAVE :: ihocl=0
INTEGER, SAVE :: ioclo=0
INTEGER, SAVE :: icl2o2=0
INTEGER, SAVE :: ibro=0
INTEGER, SAVE :: ihobr=0
INTEGER, SAVE :: ibrono2a=0
INTEGER, SAVE :: ibrcl=0
INTEGER, SAVE :: imebr=0
INTEGER, SAVE :: iccl4=0
INTEGER, SAVE :: if113=0
INTEGER, SAVE :: imecl=0
INTEGER, SAVE :: imcfm=0
INTEGER, SAVE :: if22=0
INTEGER, SAVE :: ih1211=0
INTEGER, SAVE :: ih1301=0
INTEGER, SAVE :: icof2=0
INTEGER, SAVE :: icofcl=0
INTEGER, SAVE :: ico2=0
INTEGER, SAVE :: ibrono2b=0
INTEGER, SAVE :: iho2no2b=0
INTEGER, SAVE :: icos=0
INTEGER, SAVE :: idbrm=0
INTEGER, SAVE :: ihono=0
INTEGER, SAVE :: imeono2=0
INTEGER, SAVE :: ichbr3=0
INTEGER, SAVE :: ics2=0
INTEGER, SAVE :: ih2so4=0
INTEGER, SAVE :: iso3=0


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


! upon first entry initialize positions of photolysis reactions
IF (lhook) CALL dr_hook('UKCA_STRAT_UPDATE:UKCA_STRAT_PHOTOL',                 &
                          zhook_in,zhook_handle)
IF (firstcall) THEN
  DO k=1,jppj
    ! Merge tropospheric and stratospheric photolysis rates
    SELECT CASE (spj(k,1))

      ! H2O2   + h nu -> 2 OH
    CASE ('H2O2      ')
      ih2o2 = k

      ! consider branching in case of HCHO
    CASE ('HCHO      ')
      ! HCHO   + h nu -> H2 + CO
      IF ((spj(k,3) == 'H2        ') .OR.                                      &
        (spj(k,4) == 'H2        ')) THEN
        ihchob = k
      ELSE
        ! HCHO   + h nu -> H + CHO -> CO + 2 HO2
        ihchoa = k
      END IF

      ! HO2NO2 + h nu -> NO2 + HO2
    CASE ('HO2NO2    ')
      IF ((spj(k,3) == 'HO2       ')                                           &
        .OR.(spj(k,4) == 'HO2       ')) THEN
        ! HO2NO2 + h nu -> NO2 + HO2
        iho2no2a = k
      ELSE
        ! HO2NO2 + h nu -> NO3 + OH
        iho2no2b = k
      END IF

      ! HNO3   + h nu -> NO2 + OH
    CASE ('HONO2     ')
      ihno3 = k

      ! MeOOH  + h nu -> MeO + OH -> HCHO + HO2 + OH
    CASE ('MeOOH     ')
      imeooh = k

      ! N2O5   + h nu -> NO2 + NO3
    CASE ('N2O5      ')
      in2o5 = k

      ! NO2    + h nu -> NO  + O
    CASE ('NO2       ')
      ino2 = k

      ! consider branching for NO3
    CASE ('NO3       ')
      IF ((spj(k,3) == 'O(3P)     ')                                           &
        .OR.(spj(k,4) == 'O(3P)     ')) THEN
        ! NO3    + h nu -> NO2 + O
        ino3b = k
      ELSE
        ! NO3    + h nu -> NO  + O2
        ino3a = k
      END IF

    CASE ('O2        ')
      ! O2     + h nu -> 2 O / O + O(1D)
      IF ((spj(k,3) == 'O(1D)     ')                                           &
        .OR.(spj(k,4) == 'O(1D)     ')) THEN
        io2b = k
      ELSE
        io2a = k
      END IF

      ! consider branching for O3
    CASE ('O3        ')
      IF ((spj(k,3) == 'O(1D)     ')                                           &
        .OR.(spj(k,4) == 'O(1D)     ')) THEN
        ! O3     + h nu -> O2  + O(1D)
        io3a = k
      ELSE
        ! O3     + h nu -> O2  + O
        io3b = k
      END IF

      ! consider branching for O3S
    CASE ('O3S       ')
      IF ((spj(k,3) == 'O(1D)S    ')                                           &
        .OR.(spj(k,4) == 'O(1D)S    ')) THEN
        ! O3S    + h nu -> O2  + O(1D)S
        io3sa = k
      ELSE
        ! O3S    + h nu -> O2  + O(3P)S
        io3sb = k
      END IF

      ! CH4    + h nu -> CH3   + H -> MeOO + HO2
    CASE ('CH4       ')
      ich4 = k

      ! H2O    + h nu -> OH    + H -> OH   + HO2
    CASE ('H2O       ','H2OS      ')
      ih2o = k

      ! NO     + h nu -> O     + N -> 2O    + NO
    CASE ('NO        ')
      ino = k

      ! N2O    + h nu -> O(1D) + N2
    CASE ('N2O       ')
      in2o = k

      ! F11    + h nu -> 3Cl
    CASE ('CFCl3     ')
      if11 = k

      ! F12    + h nu -> 2Cl
    CASE ('CF2Cl2    ')
      if12 = k

      ! ClONO2 + h nu -> Cl   + NO3 / ClO + NO2
    CASE ('ClONO2    ')
      IF ((spj(k,3) == 'Cl        ')                                           &
        .OR.(spj(k,4) == 'Cl        ')) THEN
        ! ClONO2 + h nu -> Cl + NO3
        iclono2a = k
      ELSE
        ! ClONO2 + h nu -> ClO + NO2
        iclono2b = k
      END IF

      ! HCl    + h nu -> H    + Cl
    CASE ('HCl       ')
      ihcl = k

      ! HOCl   + h nu -> OH   + Cl
    CASE ('HOCl      ')
      ihocl = k

      ! OClO   + h nu -> O    + ClO
    CASE ('OClO      ')
      ioclo = k

      ! Cl2O2  + h nu -> 2Cl  + O2
    CASE ('Cl2O2     ')
      icl2o2 = k

      ! BrO    + h nu -> Br   + O
    CASE ('BrO       ')
      ibro = k

      ! HOBr   + h nu -> Br   + OH
    CASE ('HOBr      ')
      ihobr = k

      ! BrONO2 + h nu -> Br   + NO3 / BrO + NO2
    CASE ('BrONO2    ')
      IF ((spj(k,3) == 'Br        ')                                           &
        .OR.(spj(k,4) == 'Br        ')) THEN
        ! BrONO2 + h nu -> Br + NO3
        ibrono2a = k
      ELSE
        ! BrONO2 + h nu -> BrO + NO2
        ibrono2b = k
      END IF

      ! BrCl   + h nu -> Br   + Cl
    CASE ('BrCl      ')
      ibrcl = k

      ! MeBr   + h nu -> Br
    CASE ('MeBr      ')
      imebr = k

      ! CCl4   + h nu -> 4 Cl
    CASE ('CCl4      ')
      iccl4 = k

      ! CF2ClCFCl2+h nu->3 Cl
    CASE ('CF2ClCFCl2')
      if113 = k

      ! CH3Cl +h nu   -> Cl
    CASE ('MeCl      ')
      imecl = k

      ! CH3CCl3 +h nu -> 3 Cl
    CASE ('MeCCl3    ')
      imcfm = k

      ! CHF2Cl +h nu -> Cl
    CASE ('CHF2Cl    ')
      if22 = k

      ! CBrClF2 +h nu -> Cl + Br
    CASE ('CF2ClBr   ')
      ih1211 = k

      ! CBrF3 +h nu   -> Br
    CASE ('CF3Br     ')
      ih1301 = k

      ! COF2 +h nu    -> 2F + CO
    CASE ('COF2      ')
      icof2 = k

      ! COFCl +h nu   -> F + Cl
    CASE ('COFCl     ')
      icofcl = k

      ! CO2 +h nu     -> CO + O(3P)
    CASE ('CO2       ')
      ico2 = k

      ! COS + h nu    -> CO + S
    CASE ('COS       ')
      icos = k

      ! HONO + h nu   -> OH + NO
    CASE ('HONO      ')
      ihono = k

      ! MeONO2 + h nu -> HO2 + HCHO + NO2
    CASE ('MeONO2    ')
      imeono2 = k
      ! CHBr3 + h nu  -> HBr + 2Br
    CASE ('CHBr3     ')
      ichbr3 = k

      ! CH2Br2 + h nu -> H2O + 2Br
    CASE ('CH2Br2    ')
      idbrm = k

      ! CS2 + h nu -> COS + SO2
    CASE ('CS2     ')
      ics2 = k
      ! H2SO4 + h nu -> SO3 + OH
    CASE ('H2SO4     ')
      ih2so4 = k
      ! SO3 + h nu -> SO2 + O(3P)
    CASE ('SO3     ')
      iso3 = k

    END SELECT
  END DO

  CALL inijtab(mype, ( (i_ukca_photol == i_ukca_fastjx) .AND.                  &
                       (fastjx_mode /= 1)) )
  firstcall = .FALSE.

END IF

! CALCJS fills the stratospheric photolysis arrays AJxyz with sensible
! values.

CALL calcjs(1, theta_field_size,                                               &
  RESHAPE(cos_zenith_angle,(/theta_field_size/)),                              &
  RESHAPE(pressure,(/theta_field_size/)),                                      &
  RESHAPE(temp,(/theta_field_size/)),                                          &
  RESHAPE(ozonecol,(/theta_field_size/)),                                      &
  theta_field_size)

! here: use only existing photolysis reactions where pressure is less than
! 300 hPa, with a linear transition into stratospheric rates

DO i=1,rows
  DO j=1,row_length
    l=(i-1) * row_length + j

    IF (pressure(j,i) < 30000.) THEN
      IF (pressure(j,i) < 20000.) THEN
        frac = 1.
      ELSE
        frac = (30000. - pressure(j,i))/10000.
      END IF

      ! Merge tropospheric and stratospheric photolysis rates
      ! H2O2   + h nu -> 2 OH
      photrates(j,i,ih2o2) = frac  * ajh2o2(l)                                 &
        + (1.-frac) * photrates(j,i,ih2o2)

      ! consider branching in case of HCHO
      ! HCHO   + h nu -> H2 + CO
      photrates(j,i,ihchob) = frac  * ajc2ob(l)                                &
        + (1.-frac) * photrates(j,i,ihchob)

      ! HCHO   + h nu -> H + CHO -> CO + 2 HO2
      photrates(j,i,ihchoa) = frac  * ajc2oa(l)                                &
        + (1.-frac) * photrates(j,i,ihchoa)

      ! HO2NO2 + h nu. Branching ratio from JPL (2002)
      IF (iho2no2a > 0) THEN
        IF (iho2no2b > 0) THEN
          ! HO2NO2 + h nu -> NO2 + HO2
          photrates(j,i,iho2no2a) = 0.667 * frac  * ajpna(l)                   &
            + (1.-frac) * photrates(j,i,iho2no2a)
          ! HO2NO2 + h nu -> NO3 + OH
          photrates(j,i,iho2no2b) = 0.333 * frac  * ajpna(l)                   &
            + (1.-frac) * photrates(j,i,iho2no2b)
        ELSE
          photrates(j,i,iho2no2a) = frac  * ajpna(l)                           &
            + (1.-frac) * photrates(j,i,iho2no2a)
        END IF
      END IF

      ! HNO3   + h nu -> NO2 + OH
      photrates(j,i,ihno3) = frac  * ajhno3(l)                                 &
        + (1.-frac) * photrates(j,i,ihno3)

      ! MeOOH  + h nu -> MeO + OH -> HCHO + HO2 + OH
      photrates(j,i,imeooh) = frac  * ajmhp(l)                                 &
        + (1.-frac) * photrates(j,i,imeooh)

      ! N2O5   + h nu -> NO2 + NO3
      photrates(j,i,in2o5) = frac  * ajn2o5(l)                                 &
        + (1.-frac) * photrates(j,i,in2o5)

      ! NO2    + h nu -> NO  + O
      photrates(j,i,ino2) = frac  * ajno2(l)                                   &
        + (1.-frac) * photrates(j,i,ino2)

      ! consider branching for NO3
      ! NO3    + h nu -> NO2 + O
      photrates(j,i,ino3b) = frac  * ajno32(l)                                 &
        + (1.-frac) * photrates(j,i,ino3b)

      ! NO3    + h nu -> NO  + O2
      photrates(j,i,ino3a) = frac  * ajno31(l)                                 &
        + (1.-frac) * photrates(j,i,ino3a)

      ! O2     + h nu -> 2 O
      photrates(j,i,io2a) = frac  * aj2a(l)                                    &
        + (1.-frac) * photrates(j,i,io2a)

      ! consider branching for O3
      ! O3     + h nu -> O2  + O(1D)
      photrates(j,i,io3a) = frac  * aj3a(l)                                    &
        + (1.-frac) * photrates(j,i,io3a)
      ! O3     + h nu -> O2  + O
      photrates(j,i,io3b) = frac  * aj3(l)                                     &
        + (1.-frac) * photrates(j,i,io3b)

      ! consider branching for O3S
      ! O3S    + h nu -> O2  + O(1D)S
      IF (io3sa > 0)                                                           &
        photrates(j,i,io3sa) = frac  * aj3a(l)                                 &
        + (1.-frac) * photrates(j,i,io3sa)
      ! O3S    + h nu -> O2  + O(3P)S
      IF (io3sb > 0)                                                           &
        photrates(j,i,io3sb) = frac  * aj3(l)                                  &
        + (1.-frac) * photrates(j,i,io3sb)

    END IF    ! pressure < 3000

    ! purely stratospheric photolysis rates
    ! SLIMCAT specific photolysis reactions

    ! O2     + h nu -> O + O(1D)
    IF (io2b > 0) THEN
      photrates(j,i,io2b) = aj2b(l)
    ELSE
      photrates(j,i,io2a) = photrates(j,i,io2a) + aj2b(l)
    END IF
    ! CH4    + h nu -> CH3   + H -> MeOO + HO2
    IF (ich4 > 0)    photrates(j,i,ich4) = ajch4(l)

    ! H2O    + h nu -> OH    + H -> OH   + HO2
    IF (ih2o > 0)    photrates(j,i,ih2o) = ajh2o(l)

    ! NO     + h nu -> O     + N -> 2O    + NO
    IF (ino > 0)     photrates(j,i,ino)  = ajno(l)

    ! N2O    + h nu -> O(1D) + N2
    IF (in2o > 0)    photrates(j,i,in2o) = ajn2o(l)

    ! F11    + h nu -> 3Cl
    IF (if11 > 0)    photrates(j,i,if11) = ajf11(l)

    ! F12    + h nu -> 2Cl
    IF (if12 > 0)    photrates(j,i,if12) = ajf12(l)

    ! ClONO2 + h nu -> Cl   + NO3
    IF (iclono2a > 0) THEN
      photrates(j,i,iclono2a) = ajcnita(l)
      IF (iclono2b == 0)                                                       &
        photrates(j,i,iclono2a) = photrates(j,i,iclono2a)                      &
        + ajcnitb(l)
    END IF

    ! ClONO2 + h nu -> ClO  + NO2
    IF (iclono2b > 0) THEN
      photrates(j,i,iclono2b) = ajcnitb(l)
      IF (iclono2a == 0)                                                       &
        photrates(j,i,iclono2b) = photrates(j,i,iclono2b)                      &
        + ajcnita(l)
    END IF

    ! HCl    + h nu -> H    + Cl
    IF (ihcl > 0)    photrates(j,i,ihcl) = ajhcl(l)

    ! HOCl   + h nu -> OH   + Cl
    IF (ihocl > 0)   photrates(j,i,ihocl) = ajhocl(l)

    ! OClO   + h nu -> O    + ClO
    IF (ioclo > 0)   photrates(j,i,ioclo) = ajoclo(l)

    ! Cl2O2  + h nu -> 2Cl  + O2
    IF (icl2o2 > 0)  photrates(j,i,icl2o2) = ajcl2o2(l)

    ! BrO    + h nu -> Br   + O
    IF (ibro > 0)    photrates(j,i,ibro) = ajbro(l)

    ! HOBr   + h nu -> Br   + OH
    IF (ihobr > 0)   photrates(j,i,ihobr) = ajhobr(l)

    ! BrONO2 + h nu -> Br   + NO3 / BrO + NO2
    ! Consider branching ratio (JPL, 2002)
    IF (ibrono2a*ibrono2b > 0) THEN
      photrates(j,i,ibrono2a) = 0.29*ajbrno3(l) ! Br + NO3 channel
      photrates(j,i,ibrono2b) = 0.71*ajbrno3(l) ! BrO + NO2 channel
    ELSE IF (ibrono2a > 0) THEN ! no branching
      photrates(j,i,ibrono2a) = ajbrno3(l)
    ELSE IF (ibrono2b > 0) THEN
      photrates(j,i,ibrono2b) = ajbrno3(l)
    END IF

    ! BrCl   + h nu -> Br   + Cl
    IF (ibrcl > 0)   photrates(j,i,ibrcl) = ajbrcl(l)

    ! MeBr   + h nu -> Br
    IF (imebr > 0)   photrates(j,i,imebr) = ajch3br(l)

    ! CCl4   + h nu -> 4 Cl
    IF (iccl4 > 0)   photrates(j,i,iccl4) = ajccl4(l)

    ! CF2ClCFCl2+h nu->3 Cl
    IF (if113 > 0)   photrates(j,i,if113) = ajf113(l)

    ! CH3Cl +h nu   -> Cl
    IF (imecl > 0)   photrates(j,i,imecl) = ajch3cl(l)

    ! CH3CCl3 +h nu -> 3 Cl
    IF (imcfm > 0)   photrates(j,i,imcfm) = ajmcfm(l)

    ! CHF2Cl +h nu -> Cl
    IF (if22 > 0)    photrates(j,i,if22) = ajf22(l)

    ! CBrClF2 +h nu -> Cl + Br
    IF (ih1211 > 0)  photrates(j,i,ih1211) = ajf12b1(l)

    ! CBrF3 +h nu   -> Br
    IF (ih1301 > 0)  photrates(j,i,ih1301) = ajf13b1(l)

    ! COF2 +h nu    -> 2F + CO
    IF (icof2 > 0)   photrates(j,i,icof2) = ajcof2(l)

    ! COFCl +h nu   -> F + Cl
    IF (icofcl > 0)  photrates(j,i,icofcl) = ajcofcl(l)

    ! CO2 +h nu     -> CO + O(3P)
    IF (ico2 > 0)    photrates(j,i,ico2) = ajco2(l)

    ! COS +h nu     -> CO + S
    IF (icos > 0)    photrates(j,i,icos) = ajcos(l)

    ! CHBr3 +h nu   -> HBr + 2Br
    IF (ichbr3 > 0)  photrates(j,i,ichbr3) = ajchbr3(l)

    ! CH2Br2 +h nu  -> H2O + 2Br
    IF (idbrm > 0)   photrates(j,i,idbrm) = ajdbrm(l)

    ! CS2 +h nu     -> COS + SO2
    IF (ics2 > 0)    photrates(j,i,ics2) = ajcs2(l)
    ! H2SO4 +h nu     -> SO3 + OH
    IF (ih2so4 > 0)  photrates(j,i,ih2so4) = ajh2so4(l)
    ! SO3 +h nu     -> SO2 + O(3P)
    IF (iso3 > 0)    photrates(j,i,iso3) = ajso3(l)

  END DO
END DO


IF (lhook) CALL dr_hook('UKCA_STRAT_UPDATE:UKCA_STRAT_PHOTOL',                 &
                          zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_strat_photol

! ######################################################################
!------------------------------------------------------------------
! Subroutine CALC_OZONECOL
!------------------------------------------------------------------

! This routine calculates overhead ozone columns in molecules/cm^2 needed
! for photolysis calculations using the Lary scheme of SLIMCAT. The
! contribution from the area above the model top is correct for a model top
! of 65 km but is not changed here for simplicity.

SUBROUTINE ukca_calc_ozonecol(model_levels, rows, row_length,                  &
  z_top_of_model, p_layer_boundaries,                                          &
  p_layer_centres, ozone_vmr, ozonecol)

USE um_parvars
USE control_max_sizes
USE ereport_mod,          ONLY: ereport
IMPLICIT NONE

! Subroutine interface

! Model dimensions
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: row_length

REAL, INTENT(IN) :: z_top_of_model       ! model top (m)
REAL, INTENT(IN) :: p_layer_boundaries(row_length, rows,                       &
  0:model_levels)
REAL, INTENT(IN) :: p_layer_centres(row_length, rows, model_levels)
REAL, INTENT(IN) :: ozone_vmr(row_length, rows, model_levels)

REAL, INTENT(OUT) :: ozonecol(row_length, rows, model_levels)
! Ozone column above level in molecules/cm^2

! local variables

! Ozone column at model top. At 38 km, calculated from 60-level ozone clima
! tology. At 85 km, extrapolated (hence a wild guess but should not matter
! too much...).  In molecules/cm^2 units

REAL, PARAMETER :: ozcol_39km = 5.0e17
REAL, PARAMETER :: ozcol_85km = 6.7e13

!REAL, PARAMETER :: colfac = 2.132e20 !! OLD VALUE
REAL, PARAMETER :: colfac = 2.117e20      ! Na/(g*M_air*1e4)

INTEGER :: l
INTEGER           :: errcode      ! Error code
CHARACTER(LEN=72) :: cmessage     !   "   message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('UKCA_STRAT_UPDATE:UKCA_CALC_OZONECOL',                &
                          zhook_in,zhook_handle)
ozonecol = 0.
DO l=model_levels-1,1,-1

  ! compute the contributions from layers above l
  ozonecol(:,:,l) = ozonecol(:,:,l+1) +                                        &
    ozone_vmr(:,:,l+1) *                                                       &
    (p_layer_boundaries(:,:,l) - p_layer_boundaries(:,:,l+1))                  &
    * colfac
END DO

! add contribution within top of level L

DO l=1,model_levels-1
  ozonecol(:,:,l) = ozonecol(:,:,l) +                                          &
    ozone_vmr(:,:,model_levels) *                                              &
    (p_layer_centres(:,:,model_levels) -                                       &
    p_layer_boundaries(:,:,model_levels)) * colfac
END DO


! add contribution above model top. For 38 levels assume model top at
! 39 km with a climatological ozone column there of 5E17 molecules/cm^2.
! For L60 & L85, the top is assumed to be at 85 km with an ozone column of
! 6.7E13 molecules/cm^2.

IF (z_top_of_model > 38500. .AND. z_top_of_model < 39500.) THEN
  ozonecol = ozonecol + ozcol_39km
ELSE IF (z_top_of_model > 84000. .AND. z_top_of_model < 85500.)                &
  THEN
  !        Both L60 (84132 km) and L85 (85000 km) versions
  ozonecol = ozonecol + ozcol_85km
ELSE
  cmessage = ' Ozone column undefined for model top'
  errcode = 1
  WRITE(6,'(A)') cmessage
  WRITE(6,'(A,E12.5)') 'z_top_of_model: ',z_top_of_model
  CALL ereport('UKCA_STRAT_UPDATE:UKCA_CALC_OZONECOL',errcode,                 &
    cmessage)
END IF

IF (lhook) CALL dr_hook('UKCA_STRAT_UPDATE:UKCA_CALC_OZONECOL',zhook_out,      &
                          zhook_handle)
RETURN
END SUBROUTINE ukca_calc_ozonecol

! ######################################################################

SUBROUTINE ukca_conserve(row_length, rows, model_levels,                       &
  ntracers, tracers,                                                           &
  pres, drain, crain, um_ozone_3d, direction)
! Description:

! This routine calculates and conserves total chlorine, bromine, and
! hydrogen. For these elements closed chemistry should be prescribed.
! Called before chemistry, with direction = .TRUE., it calculates
! total bromine, chlorine, and hydrogen as 3-D fields. Called afer
! chemistry, with direction = .FALSE., it rescales the chlorine, bromine
! and hydrogen containing compounds so that total chlorine, bromine
! and hydrogen are conserved under chemistry. Where a compound contains
! more than one of the 3 elements. e.g, BrCl, it is only scaled for the
! less abundant of the two constituents, Br. It is then subtracted from
! the total chlorine.

! Method: Rescaling of tracer variables.


! Code Description:
! Language: FORTRAN 90 + common extensions.

! Declarations:
! These are of the form:-
!     INTEGER      ExampleVariable      !Description of variable

! Global variables (*CALLed COMDECKs etc...):



USE asad_mod,           ONLY: advt
USE ukca_constants,     ONLY: c_cf2cl2, c_cfcl3, c_clo, c_cl,                  &
  c_cl2o2, c_hcl, c_hocl, c_oclo, c_brcl,                                      &
  c_clono2, c_chf2cl, c_meccl3, c_mecl,                                        &
  c_cf2clbr, c_ccl4, c_ch2br2, c_cf2clcfcl2,                                   &
  c_bro, c_br, c_hobr, c_brono2, c_mebr, c_hbr,                                &
  c_cf3br, c_h2o, c_ch4, c_h2,                                                 &
  c_ho2no2, c_hono2, c_h2o2, c_hcho, c_meooh,                                  &
  c_hono, c_c2h6, c_etooh, c_mecho, c_pan,                                     &
  c_c3h8, c_prooh, c_etcho, c_mecoch2ooh,                                      &
  c_ppan, c_meono2, c_h, c_oh, c_ho2, c_meoo,                                  &
  c_etoo, c_proo, c_etco3, c_mecoch2oo
      USE ukca_option_mod, ONLY: jpctr, L_ukca_h2o_feedback       
USE um_parvars
USE control_max_sizes
USE ereport_mod,        ONLY: ereport
IMPLICIT NONE

! ----------------------- Header file CRUNTIMC  -----------------------
! Description: Run-time constants for the Atmosphere model (read only).
!              Contains variables that define parametrization values
!              chosen for atmosphere physics and dynamics schemes.
!              [Note that CNTLATM holds accompanying control switches
!              needed for addressing.]
!
! This file belongs in section: Top Level

!
!------------------   Physics:   --------------------------------------
! Generalised physics switches:

!------------------   End of Physics   ---------------------------------

      INTEGER :: rpemax ! array size needed for diagnostic printing
      INTEGER :: rpemin ! array size needed for diagnostic printing
      INTEGER :: rpesum ! array size needed for diagnostic printing
      INTEGER :: ipesum ! array size needed for diagnostic printing
      INTEGER :: time_theta1_min      ! Timestep of min level 1 theta
      INTEGER :: time_w_max(model_levels_max) ! Timestep of max w
      INTEGER :: time_div_max(model_levels_max) ! Timestep of max div
      INTEGER :: time_div_min(model_levels_max) ! Timestep of min div
      INTEGER :: time_lapse_min(model_levels_max) ! Timestep of min
      INTEGER :: time_max_shear(model_levels_max) !Timestep max shear
      INTEGER :: time_max_wind(model_levels_max) ! Timestep of max wind
      INTEGER :: time_KE_max(model_levels_max) ! Timestep of max KE
      INTEGER :: time_KE_min(model_levels_max) ! Timestep of min KE
      INTEGER :: time_noise_max(model_levels_max) ! Timestep of max

      REAL:: frictional_timescale(model_levels_max) ! For idealised case
      REAL :: tropics_deg  ! define latitude for tropics
      REAL :: min_theta1_run                  ! Min theta level 1
      REAL :: dtheta1_run   ! Largest -ve delta theta at min theta1
      REAL :: max_w_run(0:model_levels_max) ! Max w at a level
      REAL :: max_div_run(model_levels_max) ! Max divergence at a level
      REAL :: min_div_run(model_levels_max) ! Min divergence at a level
      REAL :: min_lapse_run(model_levels_max) ! Min dtheta/dz at a level
      REAL :: max_shear_run(model_levels_max) ! Max shear at a level
      REAL :: max_wind_run(model_levels_max) ! Max wind at a level
      REAL :: max_KE_run(model_levels_max)   ! Max KE at a level
      REAL :: min_KE_run(model_levels_max)   ! Min KE at a level
      REAL :: max_noise_run(model_levels_max) ! Max noise at a level

!     Problem_number not set here  ! Now controlled by namelist input
!     Instability_diagnostics      ! Now controlled by namelist input
!     frictional_timescale         ! Now intitialised in SETCONA
!------------------   Diagnostics:   --------------------------------

      COMMON  /RUN_Diagnostics/                                         &
        rpemax, rpemin, ipesum, rpesum,                                 &
        max_w_run, min_theta1_run, dtheta1_run,                         &
        max_div_run, min_div_run, min_lapse_run,                        &
        max_shear_run, max_wind_run,                                    &
        max_noise_run, max_KE_run, min_KE_run,                          &
        time_KE_max, time_KE_min,                                       &
        time_w_max, time_div_max, time_div_min, time_lapse_min,         &
        time_max_shear, time_max_wind,                                  &
        time_theta1_min, time_noise_max

!------------------   Dynamics:   --------------------------------------
! Suarez-Held variables
      REAL :: SuHe_newtonian_timescale_ka
      REAL :: SuHe_newtonian_timescale_ks
      REAL :: SuHe_pole_equ_deltaT
      REAL :: SuHe_static_stab
      REAL :: base_frictional_timescale
      REAL :: SuHe_sigma_cutoff
      REAL :: SuHe_level_weight(model_levels_max)
      REAL :: friction_level(model_levels_max)

      INTEGER :: SuHe_relax
      INTEGER :: SuHe_fric

      LOGICAL :: L_SH_Williamson

      COMMON/Run_Dyncore/                                              &
       SuHe_newtonian_timescale_ka, SuHe_newtonian_timescale_ks,       &
       SuHe_pole_equ_deltaT, SuHe_static_stab,                         &
       base_frictional_timescale, SuHe_sigma_cutoff,                   &
       L_SH_Williamson, SuHe_relax, SuHe_fric,                         &
       SuHe_level_weight, frictional_timescale, friction_level

!------------------  Idealised model   ----------------------------

      INTEGER,PARAMETER:: max_num_profile_data = 100
      INTEGER,PARAMETER:: max_num_force_times = 100
      INTEGER,PARAMETER:: idl_max_num_bubbles = 3

! Idealised  variables
      REAL :: h_o
      REAL :: h_o_actual  ! height of growing mountain
      REAL :: h_o_per_step ! height change per step of growing mountain
      REAL :: lambda_fraction
      REAL :: phi_fraction
      REAL :: half_width_x
      REAL :: half_width_y
      REAL :: Witch_power
      REAL :: plat_size_x
      REAL :: plat_size_y
      REAL :: height_domain
      REAL :: delta_x
      REAL :: delta_y
      REAL :: big_factor
      REAL :: mag
      REAL :: vert_grid_ratio
      REAL :: first_theta_height
      REAL :: thin_theta_height
      REAL :: p_surface
      REAL :: theta_surface
      REAL :: dtheta_dz1(3)
      REAL :: height_dz1(3)
      REAL :: Brunt_Vaisala
      REAL :: u_in(4)
      REAL :: v_in(4)
      REAL :: height_u_in(3)
      REAL :: u_ramp_start
      REAL :: u_ramp_end
      REAL :: ujet_lat
      REAL :: ujet_width
      REAL :: t_horizfn_data(10)
      REAL :: q1
      REAL :: theta_ref(model_levels_max)
      REAL :: rho_ref(model_levels_max)
      REAL :: exner_ref(model_levels_max + 1)
      REAL :: q_ref(model_levels_max)
      REAL :: u_ref(model_levels_max)
      REAL :: v_ref(model_levels_max)
      REAL :: z_orog_print(0:model_levels_max)
      REAL :: f_plane
      REAL :: ff_plane
      REAL :: r_plane
      REAL :: zprofile_data(max_num_profile_data)
      REAL :: tprofile_data(max_num_profile_data)
      REAL :: qprofile_data(max_num_profile_data)
      REAL :: z_uvprofile_data(max_num_profile_data)
      REAL :: uprofile_data(max_num_profile_data)
      REAL :: vprofile_data(max_num_profile_data)
      REAL :: tforce_time_interval
      REAL :: qforce_time_interval
      REAL :: uvforce_time_interval
      REAL :: newtonian_timescale
      REAL :: z_tforce_data(max_num_profile_data)
      REAL :: z_qforce_data(max_num_profile_data)
      REAL :: z_uvforce_data(max_num_profile_data)
      REAL :: tforce_data(max_num_profile_data, max_num_force_times)
      REAL :: qforce_data(max_num_profile_data, max_num_force_times)
      REAL :: uforce_data(max_num_profile_data, max_num_force_times)
      REAL :: vforce_data(max_num_profile_data, max_num_force_times)
      REAL :: tforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: qforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: uforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: vforce_data_modlev(model_levels_max, max_num_force_times)
      REAL :: pforce_time_interval
      REAL :: p_surface_data(max_num_force_times)
      REAL :: perturb_factor
      REAL :: perturb_magnitude_t
      REAL :: perturb_magnitude_q
      REAL :: perturb_height(2)
      REAL :: orog_hgt_lbc
      REAL :: zprofile_orog
      REAL :: hf
      REAL :: cool_rate
      REAL :: IdlSurfFluxSeaParams(10) ! Idealised surface flux params
      REAL :: roughlen_z0m   
      REAL :: roughlen_z0h
      ! Idealised bubbles
      REAL :: idl_bubble_max(idl_max_num_bubbles) ! Bubble max amplitude
      REAL :: idl_bubble_height(idl_max_num_bubbles)  ! Bubble height
      REAL :: idl_bubble_width(idl_max_num_bubbles)   ! Bubble width
      REAL :: idl_bubble_depth(idl_max_num_bubbles)   ! Bubble depth
      ! Bubble x-offset, y-offset in normalised units (0:1)
      ! (0.5=domain centre)
      REAL :: idl_bubble_xoffset(idl_max_num_bubbles)
      REAL :: idl_bubble_yoffset(idl_max_num_bubbles)
      REAL :: DMPTIM, HDMP, ZDMP   ! Damping layer values
      REAL :: u_geo, v_geo         ! Geostrophic wind

! ENDGAME
      REAL :: T_surface
      REAL :: Eccentricity
      ! Following two variables used only if L_rotate_grid=.true.
      REAL :: grid_NP_lon ! Longitude (degrees) of grid's north pole
      REAL :: grid_NP_lat ! Latitude (degrees) of grid's north pole
      REAL :: AA_jet_u0   ! See QJRMS 133,1605--1614
      REAL :: AA_jet_A    !
      REAL :: theta_pert
      REAL :: ring_height
      REAL :: angular_velocity ! Planet's angular velocity
      REAL :: T0_P, T0_E ! deep atmosphere baroclinic wave surface temperatures
      INTEGER :: Trefer_number
      INTEGER :: tstep_plot_frequency
      INTEGER :: tstep_plot_start
      INTEGER :: AA_jet_m  ! See QJRMS 133,1605--1614
      INTEGER :: AA_jet_n  !
      INTEGER :: chain_number ! Run continuation number
      LOGICAL :: L_rotate_grid    ! .true. for rotating North pole of grid
      LOGICAL :: L_baro_Perturbed ! Used for baroclinic test to specify
                                  ! pert or steady
      LOGICAL :: L_shallow, L_const_grav, L_HeldSuarez,L_HeldSuarez1_drag,  &
                 L_HeldSuarez2_drag,                                        &
                 L_baro_inst, L_isothermal, L_exact_profile, L_balanced,    &
                 L_solid_body
      LOGICAL :: L_deep_baro_inst ! deep atmosphere baroclinic wave switch          


      INTEGER :: surface_type
      INTEGER :: grow_steps
      INTEGER :: grid_number
      INTEGER :: grid_flat
      INTEGER :: tprofile_number
      INTEGER :: qprofile_number
      INTEGER :: uvprofile_number
      INTEGER :: num_profile_data
      INTEGER :: num_uvprofile_data
      INTEGER :: t_horizfn_number
      INTEGER :: uv_horizfn_number
      INTEGER :: pforce_option
      INTEGER :: num_pforce_times
      INTEGER :: tforce_option
      INTEGER :: qforce_option
      INTEGER :: uvforce_option
      INTEGER :: num_tforce_levels
      INTEGER :: num_tforce_times
      INTEGER :: num_qforce_levels
      INTEGER :: num_qforce_times
      INTEGER :: num_uvforce_levels
      INTEGER :: num_uvforce_times
      INTEGER :: IdlSurfFluxSeaOption  ! Idealised surface flux option
      INTEGER :: first_constant_r_rho_level_new
      INTEGER :: big_layers
      INTEGER :: transit_layers
      INTEGER :: mod_layers
      INTEGER :: idl_bubble_option(idl_max_num_bubbles) ! Bubble option
      INTEGER :: idl_interp_option  ! Profile interpolation option
      INTEGER :: perturb_type
      INTEGER :: b_const, k_const ! deep atmosphere baroclinic wave parameters

      LOGICAL :: L_initialise_data
      LOGICAL :: L_constant_dz
      LOGICAL :: L_trivial_trigs !.false. for Cartesian coords (lat=0.0)
      LOGICAL :: L_idl_bubble_saturate(idl_max_num_bubbles)
      LOGICAL :: L_fixed_lbcs
      LOGICAL :: L_fix_orog_hgt_lbc
      LOGICAL :: L_pressure_balance
      LOGICAL :: L_wind_balance
      LOGICAL :: L_rotate_winds
      LOGICAL :: L_polar_wind_zero
      LOGICAL :: L_vert_Coriolis
      LOGICAL :: L_rotating     ! .true. for Earth's rotation
      LOGICAL :: L_perturb      ! add random perturb. to surface theta
      LOGICAL :: L_code_test    ! User switch for testing code
      LOGICAL :: L_pforce
      LOGICAL :: L_baroclinic
      LOGICAL :: L_cyclone
      LOGICAL :: L_force
      LOGICAL :: L_force_lbc
      LOGICAL :: L_perturb_t
      LOGICAL :: L_perturb_q
      LOGICAL :: L_perturb_correlate_tq
      LOGICAL :: L_perturb_correlate_vert
      LOGICAL :: L_perturb_correlate_time
      LOGICAL :: L_damp      ! Logical for damping layer
      LOGICAL :: L_geo_for ! Logical for geostrophic wind forcing
      LOGICAL :: L_bomex     ! Logical for BOMEX set up
      LOGICAL :: L_spec_z0   ! specification of roughness length    

      COMMON  /RUN_Ideal/                                              &
       h_o, h_o_actual, h_o_per_step,                                  &
       lambda_fraction, phi_fraction, half_width_x, half_width_y,      &
       Witch_power, plat_size_x, plat_size_y,                          &
       height_domain, delta_x, delta_y, big_factor, mag, vert_grid_ratio, &
       first_theta_height, thin_theta_height, p_surface,               &
       theta_surface, dtheta_dz1, height_dz1, Brunt_Vaisala,           &
       u_in, v_in, height_u_in, u_ramp_start, u_ramp_end, q1,          &
       ujet_lat, ujet_width,                                           &
       t_horizfn_number, t_horizfn_data, uv_horizfn_number,            &
       u_ref, v_ref, theta_ref, exner_ref, rho_ref, q_ref,             &
       z_orog_print, grow_steps,                                       &
       surface_type, grid_number, grid_flat,                           &
       tprofile_number, qprofile_number, uvprofile_number,             &
       num_profile_data, num_uvprofile_data,                           &
       tforce_option, qforce_option, uvforce_option,                   &
       num_tforce_levels, num_tforce_times,                            &
       num_qforce_levels, num_qforce_times,                            &
       num_uvforce_levels, num_uvforce_times,                          &
       L_pforce, pforce_option, num_pforce_times,                      &
       first_constant_r_rho_level_new,                                 &
       big_layers, transit_layers, mod_layers,                         &
       zprofile_data, tprofile_data, qprofile_data,                    &
       z_uvprofile_data, uprofile_data, vprofile_data,                 &
       tforce_time_interval, qforce_time_interval,                     &
       uvforce_time_interval, newtonian_timescale,                     &
       z_tforce_data, z_qforce_data, z_uvforce_data,                   &
       tforce_data, qforce_data, uforce_data, vforce_data,             &
       tforce_data_modlev, qforce_data_modlev,                         &
       uforce_data_modlev, vforce_data_modlev,                         &
       pforce_time_interval, p_surface_data,                           &
       L_initialise_data,                                              &
       L_perturb_t, perturb_magnitude_t,                               &
       L_perturb_q, perturb_magnitude_q,                               &
       L_perturb_correlate_tq,                                         &
       L_perturb_correlate_vert,                                       &
       L_perturb_correlate_time,                                       &
       perturb_type, perturb_height,                                   &
       L_constant_dz, L_polar_wind_zero,                               &
       L_wind_balance, L_rotate_winds,                                 &
       IdlSurfFluxSeaOption, IdlSurfFluxSeaParams,                     &
       L_spec_z0, roughlen_z0m, roughlen_z0h,                          &
       L_pressure_balance, L_vert_Coriolis,                            &
       cool_rate, L_force, L_force_lbc,                                &
       zprofile_orog, idl_interp_option, hf,                           &
       L_fix_orog_hgt_lbc, orog_hgt_lbc,                               &
       L_trivial_trigs, f_plane, ff_plane, r_plane,                    &
       idl_bubble_option, idl_bubble_max                               &
      , idl_bubble_height, idl_bubble_width, idl_bubble_depth          &
      , idl_bubble_xoffset,idl_bubble_yoffset                          &
      , L_idl_bubble_saturate,                                         &
       L_rotating, L_fixed_lbcs, L_code_test,                          &
       L_baroclinic, L_cyclone,                                        &
       L_damp, L_geo_for, L_bomex,                                     &
       DMPTIM, HDMP, ZDMP,                                             &
       u_geo, v_geo,                                                   &
!ENDGAME
       T_surface, chain_number,                                        &
       Trefer_number,                                                  &
       tstep_plot_frequency, tstep_plot_start, Eccentricity,           &
       L_rotate_grid, grid_NP_lon, grid_NP_lat,                        &
       AA_jet_m, AA_jet_n, AA_jet_u0, AA_jet_A, L_baro_Perturbed,      &
       L_shallow, L_const_grav, L_HeldSuarez,L_HeldSuarez1_drag,       &
       L_HeldSuarez2_drag,                                             &
       L_baro_inst, L_deep_baro_inst, T0_P, T0_E, b_const, k_const,    &      
       ring_height, theta_pert, L_isothermal,                          &
       L_exact_profile, L_balanced, L_solid_body, angular_velocity
! CRUNTIMC end

! Subroutine interface
INTEGER, INTENT(IN) :: row_length        ! no of points E-W
INTEGER, INTENT(IN) :: rows              ! no of points N-S
INTEGER, INTENT(IN) :: model_levels      ! no of levels
INTEGER, INTENT(IN) :: ntracers          ! no of tracers

LOGICAL, INTENT(IN) :: direction         ! T to calculate total Cl etc
! F to rescale total Cl into compouds

REAL, INTENT(IN) :: pres(row_length, rows, model_levels)  ! pressure (Pa)

REAL, INTENT(IN) :: um_ozone_3d(row_length, rows, model_levels) ! O3 MMR
REAL, INTENT(IN) :: drain(row_length, rows, model_levels)
REAL, INTENT(IN) :: crain(row_length, rows, model_levels)

REAL, INTENT(INOUT) :: tracers(row_length, rows, model_levels,                 &
  ntracers)    ! tracer array

! Local variables

! Maximum number of Br, Cl, and H compounds permitted
INTEGER, PARAMETER :: max_comp = 50

REAL, PARAMETER :: adjust_level = 500. ! pressure below which allow
! for changes of total hydrogen due to dehydration, and above which only
! advective changes are allowed.

REAL, PARAMETER :: washout_limit = 10000. ! pressure limit below which
! hydrogen conservation is not enforced.

INTEGER :: m,i,j,k

! Reservoir tracer for chlorine
INTEGER, SAVE :: n_hcl=0
! Reservoir tracer for bromine
INTEGER, SAVE :: n_brx=0
INTEGER, SAVE :: n_toth=0      ! position of total hydrogen tracer
INTEGER, SAVE :: n_h2o=0         !             H2O
INTEGER, SAVE :: n_h2os=0        !             H2O(S)
INTEGER, SAVE :: n_ox=0          !             Ox or O3
INTEGER, SAVE :: n_n2o=0         !             N2O

LOGICAL, SAVE :: firstcall = .TRUE. ! flag for first call of subr.

INTEGER, SAVE ::       ncl_tracers     ! number of Cl tracers
INTEGER, SAVE ::       nbr_tracers     ! number of Br tracers
INTEGER, SAVE ::       nh_tracers      ! number of H tracers

REAL, ALLOCATABLE, SAVE :: total_cl(:,:,:) ! total chlorine VMR
REAL, ALLOCATABLE, SAVE :: total_br(:,:,:) ! total bromine VMR
REAL, ALLOCATABLE, SAVE :: total_h (:,:,:) ! total hydrogen VMR


INTEGER, SAVE ::  cl_tracers(max_comp) ! positions of Cl tracers
INTEGER, SAVE ::  br_tracers(max_comp) ! positions of Br tracers
INTEGER, SAVE ::   h_tracers(max_comp) ! positions of hydrogen tracers

REAL   , SAVE ::c_cl_tracers(max_comp) ! conversion factors VMR/MMR
REAL   , SAVE ::c_br_tracers(max_comp) ! conversion factors
REAL   , SAVE ::c_h_tracers(max_comp)  ! conversion factors

REAL   , SAVE :: cl_validity(max_comp) ! number of Cl atoms per mol.
REAL   , SAVE :: br_validity(max_comp) ! number of Br atoms per mol.

REAL   , SAVE :: h_validity(max_comp)  ! number of H atoms per mol.

LOGICAL, SAVE :: contains_bromine(max_comp) ! flag for Cl compounds which
! also contain bromine.
LOGICAL, SAVE :: do_not_change(max_comp)    ! leave unchanged.

! Correction factor to achieve Cl or Br conservation
REAL, ALLOCATABLE :: corrfac(:,:,:)

INTEGER           :: errcode        ! error code
CHARACTER(LEN=72) :: cmessage       !   "   message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('UKCA_STRAT_UPDATE:CONSERVE',zhook_in,zhook_handle)
IF (firstcall) THEN

  IF (ntracers /= jpctr) THEN
    errcode = 1
    cmessage=' Number of tracers does not agree with jpctr'
    WRITE(6,'(1X,A,1X,I4,1X,I4)') cmessage,ntracers,jpctr
    CALL ereport('UKCA_STRAT_UPDATE:UKCA_CONSERVE',errcode,                    &
      cmessage)
  END IF

  ! Calculate the number of Cl and Br tracers, their positions, mmr/vmr
  ! conversion ratios, and validities (numbers of Cl/Br atoms per molecule).

  ncl_tracers = 0
  nbr_tracers = 0
  nh_tracers = 0
  contains_bromine = .FALSE.
  do_not_change = .FALSE.
  cl_tracers = 0
  br_tracers = 0
  h_tracers = 0
  cl_validity = 1.
  br_validity = 1.
  h_validity = 1.
  c_br_tracers = 0.
  c_cl_tracers = 0.
  c_h_tracers = 0.
  DO m=1,jpctr
    SELECT CASE (advt(m))

      ! chlorine tracers
    CASE ('CF2Cl2    ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_cf2cl2
      cl_validity(ncl_tracers) = 2.
    CASE ('CFCl3     ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_cfcl3
      cl_validity(ncl_tracers) = 3.
    CASE ('Clx       ','ClO       ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_clo
    CASE ('Cl        ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_cl
    CASE ('Cl2O2     ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_cl2o2
      cl_validity(ncl_tracers) = 2.
    CASE ('HCl       ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_hcl
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hcl
      do_not_change(nh_tracers) = .TRUE.
      n_hcl = m
    CASE ('HOCl      ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_hocl
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hocl
      do_not_change(nh_tracers) = .TRUE.
    CASE ('OClO      ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_oclo
    CASE ('BrCl      ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_brcl
      contains_bromine(ncl_tracers) = .TRUE.
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_brcl
    CASE ('ClONO2    ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_clono2
    CASE ('CF2ClCFCl2')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_cf2clcfcl2
      cl_validity(ncl_tracers) = 3.
    CASE ('CHF2Cl    ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_chf2cl
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_chf2cl
      do_not_change(nh_tracers) = .TRUE.
    CASE ('MeCCl3    ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_meccl3
      cl_validity(ncl_tracers) = 3.
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_meccl3
      do_not_change(nh_tracers) = .TRUE.
    CASE ('CCl4      ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_ccl4
      cl_validity(ncl_tracers) = 4.
    CASE ('MeCl      ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_mecl
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_mecl
      do_not_change(nh_tracers) = .TRUE.
    CASE ('CF2ClBr   ')
      ncl_tracers = ncl_tracers + 1
      cl_tracers(ncl_tracers) = m
      c_cl_tracers(ncl_tracers) = c_cf2clbr
      contains_bromine(ncl_tracers) = .TRUE.
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_cf2clbr
      ! Bromine tracers
    CASE ('Brx       ','BrO       ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_bro
      n_brx = m
    CASE ('Br        ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_br
    CASE ('HOBr      ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_hobr
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hobr
      do_not_change(nh_tracers) = .TRUE.
    CASE ('BrONO2    ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_brono2
    CASE ('MeBr      ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_mebr
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_mebr
      do_not_change(nh_tracers) = .TRUE.
    CASE ('HBr       ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_hbr
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hbr
      do_not_change(nh_tracers) = .TRUE.
    CASE ('CF3Br     ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_cf3br
    CASE ('CH2Br2    ')
      nbr_tracers = nbr_tracers + 1
      br_tracers(nbr_tracers) = m
      c_br_tracers(nbr_tracers) = c_ch2br2
      br_validity(nbr_tracers) = 2.
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_ch2br2
      h_validity(nh_tracers) = 2.
      do_not_change(nh_tracers) = .TRUE.
      ! hydrogen tracers
    CASE ('TOTH      ')
      n_toth = m
    CASE ('H2O       ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_h2o
      h_validity(nh_tracers) = 2.
      n_h2o = m
    CASE ('CH4       ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_ch4
      h_validity(nh_tracers) = 4.
      do_not_change(nh_tracers) = .TRUE.
    CASE ('H2        ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_h2
      h_validity(nh_tracers) = 2.
    CASE ('HO2NO2    ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_ho2no2
      do_not_change(nh_tracers) = .TRUE.
    CASE ('HONO2     ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hono2
      do_not_change(nh_tracers) = .TRUE.
    CASE ('H2O2      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_h2o2
      h_validity(nh_tracers) = 2.
    CASE ('HCHO      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hcho
      h_validity(nh_tracers) = 2.
    CASE ('MeOOH     ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_meooh
      h_validity(nh_tracers) = 4.
    CASE ('HONO      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_hono
      do_not_change(nh_tracers) = .TRUE.
    CASE ('C2H6      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_c2h6
      h_validity(nh_tracers) = 6.
    CASE ('EtOOH     ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_etooh
      h_validity(nh_tracers) = 6.
    CASE ('MeCHO     ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_mecho
      h_validity(nh_tracers) = 4.
    CASE ('PAN       ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_pan
      h_validity(nh_tracers) = 3.
      do_not_change(nh_tracers) = .TRUE.
    CASE ('C3H8      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_c3h8
      h_validity(nh_tracers) = 8.
    CASE ('n-PrOOH   ','i-PrOOH   ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_prooh
      h_validity(nh_tracers) = 8.
    CASE ('EtCHO     ','Me2CO     ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_etcho
      h_validity(nh_tracers) = 6.
    CASE ('MeCOCH2OOH')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_mecoch2ooh
      h_validity(nh_tracers) = 6.
    CASE ('PPAN      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_ppan
      h_validity(nh_tracers) = 5.
      do_not_change(nh_tracers) = .TRUE.
    CASE ('MeONO2    ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_meono2
      h_validity(nh_tracers) = 3.
      do_not_change(nh_tracers) = .TRUE.
    CASE ('H         ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_h
    CASE ('OH        ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_oh
    CASE ('HO2       ','HOx       ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_ho2
    CASE ('MeOO      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_meoo
      h_validity(nh_tracers) = 3.
    CASE ('EtOO      ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_etoo
      h_validity(nh_tracers) = 5.
    CASE ('i-PrOO    ','n-PrOO    ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_proo
      h_validity(nh_tracers) = 7.
    CASE ('EtCO3     ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_etco3
      h_validity(nh_tracers) = 5.
    CASE ('MeCOCH2OO ')
      nh_tracers = nh_tracers + 1
      h_tracers(nh_tracers) = m
      c_h_tracers(nh_tracers) = c_mecoch2oo
      h_validity(nh_tracers) = 5.
      ! special tracers
    CASE('Ox        ','O3        ')
      n_ox = m
    CASE('N2O       ')
      n_n2o = m
    CASE default
      cmessage = ' CONSERVE: SPECIES '//advt(m)//                              &
        ' not treated in CASE'
      WRITE(6,'(1X,I4,1X,A)') m,advt(m)
      errcode = -1*m
      CALL ereport('UKCA_STRAT_UPDATE:CONSERVE',errcode,cmessage)
    END SELECT
  END DO

  firstcall = .FALSE.

END IF  ! if firstcall

IF (direction) THEN
  ! Calculate total chlorine and bromine vmrs.
  IF (.NOT. ALLOCATED(total_br))                                               &
    ALLOCATE(total_br(row_length, rows, model_levels))
  IF (.NOT. ALLOCATED(total_cl))                                               &
    ALLOCATE(total_cl(row_length, rows, model_levels))

  total_br = 0.
  total_cl = 0.

  DO m=1,nbr_tracers
    total_br = total_br + tracers(:,:,:,br_tracers(m))*                        &
      br_validity(m) / c_br_tracers(m)
  END DO
  DO m=1,ncl_tracers
    total_cl = total_cl + tracers(:,:,:,cl_tracers(m))*                        &
      cl_validity(m) / c_cl_tracers(m)
  END DO
  ! Do not do hydrogen conservation unless at least two hydrogen
  ! reservoirs (H2O, CH4) are defined.
  IF (l_ukca_h2o_feedback) THEN
    ALLOCATE(total_h(row_length, rows, model_levels))
    total_h = 0.
    DO m=1,nh_tracers
      total_h = total_h + tracers(:,:,:,h_tracers(m))*                         &
        h_validity(m) / c_h_tracers(m)
    END DO
    IF (n_toth > 0) THEN
      WHERE (pres > adjust_level)                                              &
        tracers(:,:,:,n_toth) = total_h
      WHERE (pres <= adjust_level)                                             &
        total_h = tracers(:,:,:,n_toth)
    END IF
  END IF

ELSE    ! direction = F
  IF (model_levels > 38) THEN
    tracers(:,:,model_levels-1:model_levels,n_ox) =                            &
      um_ozone_3d(:,:,model_levels-1:model_levels)
  END IF
  ! Set N2O to 0 at the top, to avoid initialization problem
  IF (n_n2o > 0)                                                               &
    tracers(:,:,model_levels,n_n2o) = 0.

  ! Adjust tracers to match
  ALLOCATE(corrfac(row_length, rows, model_levels))
  IF (nbr_tracers > 0) THEN
    corrfac = 0.

    ! Calculate new total bromine
    DO m=1,nbr_tracers
      corrfac = corrfac + tracers(:,:,:,br_tracers(m))*                        &
        br_validity(m) / c_br_tracers(m)
    END DO

    ! Adjust bromine tracers to match total bromine computed before
    ! chemistry

    corrfac = total_br / corrfac
    DO m=1,nbr_tracers
      DO i=1,row_length
        DO j=1,rows
          DO k=2,model_levels
            IF (drain(i,j,k) + crain(i,j,k) == 0.) THEN
              tracers(i,j,k,br_tracers(m)) =                                   &
                tracers(i,j,k,br_tracers(m)) * corrfac(i,j,k)
            ELSE
              corrfac(i,j,k) = 1.
            END IF
          END DO
        END DO
      END DO
    END DO
    IF (((MINVAL(corrfac) < 0.9) .OR. (MAXVAL(corrfac) > 1.1))                 &
      .AND. (printstatus >= prstatus_normal)) THEN
      WRITE (6,'(A,2F8.4)') 'Correct bromine ',MINVAL(corrfac),MAXVAL(corrfac)
    END IF
  END IF     ! nbr_tracers > 0

  ! Calculate new total chlorine, excluding BrCl and CF2ClBr

  IF (ncl_tracers > 0) THEN
    corrfac = 0.

    DO m=1,ncl_tracers
      IF (.NOT.(contains_bromine(m))) THEN
        corrfac = corrfac + tracers(:,:,:,cl_tracers(m))*                      &
          cl_validity(m) / c_cl_tracers(m)
      ELSE
        total_cl = total_cl - tracers(:,:,:,cl_tracers(m))*                    &
          cl_validity(m) / c_cl_tracers(m)
      END IF
    END DO

    ! Adjust chlorine species to match total chlorine computed before.
    ! Leave BrCl and CF2ClBr alone.

    corrfac = total_cl / corrfac
    DO m=1,ncl_tracers
      IF (.NOT.(contains_bromine(m))) THEN
        DO i=1,row_length
          DO j=1,rows
            DO k=2,model_levels
              IF (crain(i,j,k) + drain(i,j,k) == 0.) THEN
                tracers(i,j,k,cl_tracers(m)) =                                 &
                  tracers(i,j,k,cl_tracers(m)) * corrfac(i,j,k)
              ELSE
                corrfac(i,j,k) = 1.
              END IF
            END DO
          END DO
        END DO
      END IF
    END DO
    IF (((MINVAL(corrfac) < 0.9) .OR. (MAXVAL(corrfac) > 1.1))                 &
      .AND. (printstatus >= prstatus_normal)) THEN
      WRITE(6,'(A,2F8.4)')'Correct chlorine ',MINVAL(corrfac),MAXVAL(corrfac)
    END IF
  END IF   ! ncl_tracers > 0

  IF (l_ukca_h2o_feedback) THEN
    corrfac = 0.

    ! Calculate new total hydrogen
    DO m=1,nh_tracers
      IF (.NOT.(do_not_change(m))) THEN
        corrfac = corrfac + tracers(:,:,:,h_tracers(m))*                       &
          h_validity(m) / c_h_tracers(m)
      ELSE
        total_h = total_h - tracers(:,:,:,h_tracers(m))*                       &
          h_validity(m) / c_h_tracers(m)
      END IF
    END DO

    ! adjust upper boundary for hydrogen tracers only. This is needed to
    ! prevent model instability if too much water is present at the model
    ! top.
    ! Adjust hydrogen tracers to match total hydrogen computed before
    ! chemistry

    corrfac = total_h / corrfac

    ! Do not enforce hydrogen conservation in the troposphere, due to
    ! washout and dry deposition.

    WHERE (pres > washout_limit) corrfac = 1.

    DO m=1,nh_tracers
      IF (.NOT.(do_not_change(m)))                                             &
        tracers(:,:,:,h_tracers(m)) = tracers(:,:,:,h_tracers(m))              &
        * corrfac
    END DO
    IF (ALLOCATED(total_h)) DEALLOCATE(total_h)
  END IF ! L_ukca_h2o_feedback

  ! Deallocate fields
  IF (ALLOCATED(total_br)) DEALLOCATE(total_br)
  IF (ALLOCATED(total_cl)) DEALLOCATE(total_cl)
  IF (ALLOCATED(corrfac)) DEALLOCATE(corrfac)
END IF   ! direction

IF (lhook) CALL dr_hook('UKCA_STRAT_UPDATE:CONSERVE',zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_conserve

! ######################################################################

SUBROUTINE ukca_relax_ozone(rows, row_length, model_levels,                    &
  pressure, ozone_mmr, um_ozone_3d)


!  This subroutine replaces chemical ozone with climatological ozone
!  above 0.2 hPa, with a transition region between 0.2 and 0.3 hPa.
!  This is meant to improve ozone in this region and allow for better
!  transmission of UV radiation through the mesosphere.

USE um_parvars
USE control_max_sizes
IMPLICIT NONE

! Subroutine interface

INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: model_levels
REAL, INTENT(IN) :: pressure(row_length, rows, model_levels)
REAL, INTENT(IN) :: um_ozone_3d(row_length, rows, model_levels)

REAL, INTENT(INOUT) :: ozone_mmr(row_length, rows, model_levels)

! Local vaiables:

REAL, PARAMETER :: maxpres = 30.
REAL, PARAMETER :: minpres = 20.

INTEGER :: i
INTEGER :: j
INTEGER :: k

REAL :: frac

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook('UKCA_STRAT_UPDATE:UKCA_RELAX_OZONE',                  &
  zhook_in,zhook_handle)
DO k=1,model_levels
  IF (MINVAL(pressure(:,:,k)) < maxpres) THEN
    DO i=1,rows
      DO j=1,row_length
        IF (pressure(j,i,k) < maxpres) THEN
          IF (pressure(j,i,k) < minpres) THEN
            ozone_mmr(j,i,k) = um_ozone_3d(j,i,k)
          ELSE
            frac = (maxpres-pressure(j,i,k))/(maxpres - minpres)
            ozone_mmr(j,i,k) = frac  *    um_ozone_3d(j,i,k) +                 &
              (1. - frac) * ozone_mmr(j,i,k)
          END IF
        END IF
      END DO
    END DO
  END IF
END DO

IF (lhook) CALL dr_hook('UKCA_STRAT_UPDATE:UKCA_RELAX_OZONE',                  &
  zhook_out,zhook_handle)
RETURN
END SUBROUTINE ukca_relax_ozone

END MODULE ukca_strat_update
