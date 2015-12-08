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
! Description:
!  Main driver routine for chemistry
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!   Called from UKCA_MAIN1.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!------------------------------------------------------------------
!
      SUBROUTINE UKCA_CHEMISTRY_CTL(i_month, i_day_number, i_hour,     &
                      i_minute, secs_per_step,                         &
                      ntracers,                                        &
                      ndiags,                                          &
                      sinlat,                                          &
                      coslat,                                          &
                      true_longitude,                                  &
                      pres, temp, q,                                   &
                      qcf, qcl, rh,                                    &
                      p_layer_boundaries,                              &
                      r_theta_levels,                                  &
                      z_top_of_model,                                  &
                      cos_zenith_angle,                                &
                      tracer,                                          &
                      user_diags,                                      &
                      t_surf, dzl, z0m, u_s,                           &
                      drain, crain,                                    &
                      cloud_frac,                                      &
                      fastj_dj,                                        &
                      volume, mass,                                    &
                      land_points, land_index,                         &
                      tile_pts, tile_index, tile_frac,                 &
                      zbl, surf_hf, seaice_frac, stcon,                &
                      soilmc_lp, fland, laift_lp, canhtft_lp,          &
                      z0tile_lp, t0tile_lp, canwctile_lp,              &
                      pv_at_theta,                                     &
                      theta,                                           &
                      um_ozone3d,                                      &
                      uph2so4inaer,                                    &
                      delso2_wet_h2o2,                                 &
                      delso2_wet_o3,                                   &
                      delh2so4_chem,                                   &
                      delso2_drydep,                                   &
                      delso2_wetdep,                                   &
                      so4_sa                                           &
                      )

      USE conversions_mod, ONLY: pi_over_180, pi
      USE nstypes, ONLY: ntype, npft
      USE ASAD_MOD
      USE ASAD_CHEM_FLUX_DIAGS
      USE UKCA_D1_DEFS
      USE UKCA_CSPECIES
      USE UKCA_tropopause
      USE UKCA_strat_update
      USE UKCA_CONSTANTS,       ONLY: c_o3, c_h2o, c_hono2, c_o1d
      USE UKCA_dissoc
      USE UKCA_phot2d,          ONLY: UKCA_PHOTIN, UKCA_CURVE,         &
                                      UKCA_INPR2D, nolev, nlphot,      &
                                      ntphot, pjin, pr2d
      USE ukca_option_mod,      ONLY: L_ukca_raq, L_ukca_trop,         &
                                      L_ukca_trophet, L_ukca_use_2dtop,&
                                      L_ukca_h2o_feedback, L_ukca_chem,&
                                      L_ukca_het_psc, L_ukca_achem,    &
                                      L_ukca_advh2o, L_ukca_aerchem,   &
                                      L_ukca_intdd,      &
                                      i_ukca_photol, fastjx_prescutoff,&
                                      dts0
      USE ukca_photo_scheme_mod, ONLY: i_ukca_phot2d, i_ukca_fastj,    &
                                       i_ukca_fastjx
      USE earth_constants_mod,  ONLY: g
                                           
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE PrintStatus_mod
      USE UM_ParVars
      USE Control_Max_Sizes
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
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------
!-------------------COMDECK C_SULCHM--------------------------------
! Parameters for Sulphur Cycle Chemistry
      REAL                                                              &
     &     EVAPTAU,                                                     &
                          ! timescale for dissolved SO4 to evaporate
     &     NUCTAU,                                                      &
                          ! timescale for accumulation mode particles
!                           to nucleate once they enter a cloud.
     &     DIFFUSE_AIT,                                                 &
                          ! diffusion coefficient of Aitken particles
     &     K_SO2OH_HI,                                                  &
                                  ! high pressure reaction rate limit
     &     K_DMS_OH,                                                    &
                                  ! reaction rate for DMS+OH  cc/mcl/s
     &      K4_CH3SO2_O3,                                               &
                             ! Rate coeff for CH3SO2+O3 -> CH3SO3+O2
     &      K5_CH3SO3_HO2,                                              &
                             ! Rate coeff for CH3SO3+HO2 -> MSA+O2
     &      RMM_O3,                                                     &
                             ! relative molecular mass O3
     &     BRAT_SO2,                                                    &
                                  ! branching ratio for SO2 in DMS oxidn
     &     BRAT_MSA,                                                    &
                                  ! branching ratio for MSA in DMS oxidn
     &     AVOGADRO,                                                    &
                                 ! no. of molecules in 1 mole
     &     RMM_S,                                                       &
                                 ! relative molecular mass S kg/mole
     &     RMM_H2O2,                                                    &
                                 ! relative molecular mass H2O2 kg/mole
     &     RMM_HO2,                                                     &
                                 ! relative molecular mass HO2 kg/mole
     &     RMM_AIR,                                                     &
                                 ! relative molecular mass dry air
     &     RMM_W,                                                       &
                                 ! relative molecular mass water
     &     RELM_S_H2O2,                                                 &
                                 ! rel atomic mass sulphur/RMM_H2O2
     &     RELM_S_2N,                                                   &
                              ! rel atomic mass Sulphur/2*Nitrogen
     &     PARH,                                                        &
                                ! power of temp dependence of K_SO2OH_LO
     &     K1,                                                          &
                                ! parameters for calcn of K_SO2OH_LO
     &     T1,                                                          &
                                !
     &     FC,                                                          &
                                ! parameters for interpolation between
     &     FAC1,                                                        &
                                !   LO and HI reaction rate limits
     &     K2,K3,K4,                                                    &
                                ! parameters for calcn of K_HO2_HO2
     &     T2,T3,T4,                                                    &
                                !
     &     CLOUDTAU,                                                    &
                                  ! air parcel lifetime in cloud
     &     CHEMTAU,                                                     &
                                  ! chem lifetime in cloud before oxidn
     &     O3_MIN,                                                      &
                              ! min mmr of O3 required for oxidn
     &     THOLD                  ! threshold for cloud liquid water
!
!
      PARAMETER (                                                       &
     &           EVAPTAU = 300.0,                                       &
                                              ! secs  (=5 mins)
     &             NUCTAU = 30.0,                                       &
                                          ! secs
     &       DIFFUSE_AIT = 1.7134E-9,                                   &
                                             ! sq m/s
     &        K_SO2OH_HI = 2.0E-12,                                     &
                                       ! cc/mcl/s from STOCHEM model
     &           K_DMS_OH = 9.1E-12,                                    &
                                          ! cc/mcl/s
     &       K4_CH3SO2_O3 = 1.0E-14,                                    &
                                        ! cc/mcl/s
     &      K5_CH3SO3_HO2 = 4.0E-11,                                    &
     &             RMM_O3 = 4.8E-2,                                     &
                                        ! kg/mole
     &          BRAT_SO2 = 0.9,                                         &
     &           BRAT_MSA = 1.0-BRAT_SO2,                               &
     &           AVOGADRO = 6.022E23,                                   &
                                          ! per mole
     &           RMM_S    = 3.20E-2,                                    &
                                          ! kg/mole
     &           RMM_H2O2 = 3.40E-2,                                    &
                                          ! kg/mole
     &           RMM_HO2  = 3.30E-2,                                    &
                                          ! kg/mole
     &            RMM_AIR = 2.896E-2,                                   &
                                          ! kg/mole
     &              RMM_W = 1.8E-2,                                     &
                                          ! kg/mole
     &        RELM_S_H2O2 = 3.206/3.40,                                 &
     &           RELM_S_2N = 3.206/2.80,                                &
     &               PARH = 3.3,                                        &
     &                K1 = 4.0E-31,                                     &
                                       ! (cc/mcl)2/s from STOCHEM
     &                 T1 = 300.0,                                      &
                                          ! K
     &                FC = 0.45,                                        &
                                        ! from STOCHEM model
     &              FAC1 = 1.1904,                                      &
                                    ! 0.75-1.27*LOG10(FC) from STOCHEM
     &                 K2 = 2.2E-13,                                    &
                                          ! cc/mcl/s
     &                 K3 = 1.9E-33,                                    &
                                          ! (cc/mcl)2/s
     &                 K4 = 1.4E-21,                                    &
                                          ! cc/mcl
     &                 T2 = 600.0,                                      &
                                          ! K
     &                 T3 = 890.0,                                      &
                                          ! K
     &                 T4 = 2200.0,                                     &
                                          ! K
     &           CLOUDTAU = 1.08E4,                                     &
                                          ! secs (=3 hours)
     &            CHEMTAU = 9.0E2,                                      &
                                          ! secs (=15 mins)
     &              O3_MIN = 1.6E-8,                                    &
                                        !(kg/kg, equiv. 10ppbv)
     &              THOLD = 1.0E-8                                      &
                                          ! kg/kg
     &          )
!
      REAL RAD_AIT,                                                     &
                            ! median radius of Aitken mode particles
     &     DIAM_AIT,                                                    &
                            !   "    diameter    "
     &     RAD_ACC,                                                     &
                            ! median radius of acccumulation mode
     &     DIAM_ACC,                                                    &
                            !   "    diameter    "
     &     CHI,                                                         &
                            ! mole fraction of S in particle
     &     RHO_SO4,                                                     &
                            ! density of  SO4 particle
     &     SIGMA,                                                       &
                            ! standard devn of particle size distn
!                                 for accumulation mode
     &     E_PARM,                                                      &
                            ! param relating size distns of Ait & Acc
     &     NUM_STAR         ! threshold concn of accu mode particles
                            !  below which PSI=1
!
      PARAMETER (                                                       &
     &           RAD_AIT = 6.5E-9,                                      &
                                             ! m
     &          DIAM_AIT = 2.0*RAD_AIT,                                 &
     &           RAD_ACC = 95.0E-9,                                     &
                                             ! m
     &          DIAM_ACC = 2.0*RAD_ACC,                                 &
     &               CHI = 32.0/132.0,                                  &
     &           RHO_SO4 = 1769.0,                                      &
                                              ! kg/m3
     &             SIGMA = 1.4,                                         &
     &            E_PARM = 0.9398,                                      &
     &          NUM_STAR = 1.0E6                                        &
                                             ! m-3
     &          )
!
      REAL BOLTZMANN       !Boltzmanns constant.
      REAL MFP_REF         !Reference value of mean free path.
      REAL TREF_MFP        !Reference temperature for mean free path.
      REAL PREF_MFP        !Reference pressure for mean free path.
      REAL SIGMA_AIT       !Geometric standard deviation of the Aitken
!                             mode distribution.
!
      PARAMETER (BOLTZMANN = 1.3804E-23)  ! J K-1
      PARAMETER (MFP_REF = 6.6E-8                                       &
                                          ! m
     &        ,  TREF_MFP = 293.15                                      &
                                          ! K
     &        ,  PREF_MFP = 1.01325E5)    ! Pa
      PARAMETER (SIGMA_AIT = 1.30)
!
!*---------------------------------------------------------------------
! TYPSTS starts
! submodel_mod must be included before this file
!Applicable to all configurations
!STASH related variables for describing output requests and space
!management.
! Include sizes for dimensioning arrays in this deck
! TYPSTSZ start
!  Sizes derived from STASHC file of UMUI job, and includes those
!  sizes needed to dimension arrays in TYPSTS .h deck.

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: LEN_STLIST   = 33

      ! No of items per timeseries recd
      INTEGER, PARAMETER :: TIME_SERIES_REC_LEN = 9

      INTEGER :: NSECTS               ! Max no of diagnostic sections
      INTEGER :: N_REQ_ITEMS          ! Max item number in any section
      INTEGER :: NITEMS               ! No of distinct items requested
      INTEGER :: N_PPXRECS            ! No of PP_XREF records this run
      INTEGER :: TOTITEMS             ! Total no of processing requests
      INTEGER :: NSTTIMS              ! Max no of STASHtimes in a table
      INTEGER :: NSTTABL              ! No of STASHtimes tables
      INTEGER :: NUM_STASH_LEVELS     ! Max no of levels in a levelslist
      INTEGER :: NUM_LEVEL_LISTS      ! No of levels lists
      INTEGER :: NUM_STASH_PSEUDO     ! Max no of pseudo-levs in a list
      INTEGER :: NUM_PSEUDO_LISTS     ! No of pseudo-level lists
      INTEGER :: NSTASH_SERIES_BLOCK  ! No of blocks of timeseries recds
      INTEGER :: NSTASH_SERIES_RECORDS! Total no of timeseries records

      COMMON/STSIZES_TYPSTS/                                            &
     &  NSECTS,N_REQ_ITEMS,NITEMS,N_PPXRECS,TOTITEMS,NSTTABL,           &
     &  NUM_STASH_LEVELS,NUM_LEVEL_LISTS,NUM_STASH_PSEUDO,              &
     &  NUM_PSEUDO_LISTS,NSTTIMS,NSTASH_SERIES_BLOCK,                   &
     &        NSTASH_SERIES_RECORDS


! TYPSTSZ end

! This file is needed to get ppxref_codelen to dimension PP_XREF
      ! sizes in STASH used for defining local array dimensions at a
      ! lower level.
      INTEGER :: MAX_STASH_LEVS  ! Max no of output levels for any diag
      INTEGER :: PP_LEN2_LOOKUP  ! Max no of LOOKUPs needed in STWORK
      COMMON/CARGST/MAX_STASH_LEVS,PP_LEN2_LOOKUP

      ! STASHflag (.TRUE. for processing this timestep). SF(0,IS) .FALSE.
      ! if no flags on for section IS.
      LOGICAL :: SF(0:NITEMS,0:NSECTS)

      ! Whether a calculation is needed for SF above
      LOGICAL :: SF_CALC(0:NITEMS,0:NSECTS)

      ! STASH list index
      INTEGER :: STINDEX(2,NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! List of STASH output requests
      INTEGER :: STLIST (LEN_STLIST,TOTITEMS)

      ! Address of item from generating plug compatible routine (often
      ! workspace)
      INTEGER :: SI     (  NITEMS,0:NSECTS,N_INTERNAL_MODEL)

      ! STASH times tables
      INTEGER :: STTABL (NSTTIMS,NSTTABL)

      ! Length of STASH workspace required in each section
      INTEGER:: STASH_MAXLEN       (0:NSECTS,N_INTERNAL_MODEL          )
      INTEGER:: PPINDEX            (  NITEMS,N_INTERNAL_MODEL          )
      INTEGER:: STASH_LEVELS       (NUM_STASH_LEVELS+1,NUM_LEVEL_LISTS )
      INTEGER:: STASH_PSEUDO_LEVELS(NUM_STASH_PSEUDO+1,NUM_PSEUDO_LISTS)
      INTEGER:: STASH_SERIES(TIME_SERIES_REC_LEN,NSTASH_SERIES_RECORDS)
      INTEGER:: STASH_SERIES_INDEX(2,NSTASH_SERIES_BLOCK)
! TYPSTS end

        INTEGER, INTENT(IN) :: ntracers          ! no. of tracers
        INTEGER, INTENT(IN) :: ndiags            ! no. of diagnostics
        INTEGER, INTENT(IN) :: i_month           ! month
        INTEGER, INTENT(IN) :: i_day_number      ! day
        INTEGER, INTENT(IN) :: i_hour            ! hour
        INTEGER, INTENT(IN) :: i_minute          ! minute
        INTEGER, INTENT(IN) :: uph2so4inaer      ! flag for H2SO4 updating

!       Variables for interactive dry deposition scheme
        INTEGER, INTENT(IN) :: land_points
        INTEGER, INTENT(IN) :: land_index(land_points)
        INTEGER, INTENT(IN) :: tile_pts(ntype)
        INTEGER, INTENT(IN) :: tile_index(land_points,ntype)

        REAL, INTENT(IN) :: secs_per_step                      ! time step
        REAL, INTENT(IN) :: z_top_of_model                     ! top of model (m)
        REAL, INTENT(INOUT) :: sinlat(row_length, rows)        ! sin(latitude)
        REAL, INTENT(IN) :: coslat(row_length, rows)           ! cos(latitude)
        REAL, INTENT(IN) :: true_longitude(row_length,rows)    ! longitude
        REAL, INTENT(IN) :: pres(row_length,rows,model_levels) ! pressure
        REAL, INTENT(IN) :: p_layer_boundaries(row_length,rows,    &
                                               0:model_levels) ! pressure
        REAL, INTENT(IN) :: r_theta_levels(row_length,rows,             &
                                           0:model_levels)
        REAL, INTENT(IN) :: temp(row_length,rows,model_levels) ! actual temp
        REAL, INTENT(IN) :: dzl(row_length, rows, bl_levels)   ! thickness
        REAL, INTENT(IN) :: u_s(row_length, rows)              ! ustar
        REAL, INTENT(IN) :: z0m(row_length, rows)              ! roughness
        REAL, INTENT(IN) :: t_surf(row_length, rows)           ! surface temp
        REAL, INTENT(IN) :: drain(row_length,rows,model_levels) ! 3-D LS rain
        REAL, INTENT(IN) :: crain(row_length,rows,model_levels) ! 3-D convec
        REAL, INTENT(IN) :: cos_zenith_angle(row_length, rows) ! cosine of ZA
        REAL, INTENT(IN) :: volume(row_length,rows,model_levels) ! cell vol.
        REAL, INTENT(IN) :: mass(row_length, rows, model_levels) ! cell mass
        REAL, INTENT(IN) :: pv_at_theta(row_length, rows, model_levels) ! PV
        REAL, INTENT(IN) :: theta(row_length, rows, model_levels)! theta
        REAL, INTENT(IN) :: um_ozone3d(row_length, rows, model_levels) ! O3
        REAL, INTENT(IN) :: qcf(row_length, rows, wet_levels)  ! qcf
        REAL, INTENT(IN) :: qcl(row_length, rows, wet_levels)  ! qcl
        REAL, INTENT(IN) :: rh(row_length, rows, wet_levels)    ! RH frac
        REAL, INTENT(IN) :: cloud_frac(row_length, rows, wet_levels)
        REAL, INTENT(IN) :: so4_sa(row_length,rows,model_levels)  ! aerosol
!                                                        surface area

!       Variables for interactive dry deposition scheme

        REAL, INTENT(IN) :: tile_frac(land_points,ntype)
        REAL, INTENT(IN) :: zbl(row_length,rows)
        REAL, INTENT(IN) :: surf_hf(row_length,rows)
        REAL, INTENT(IN) :: seaice_frac(row_length,rows)
        REAL, INTENT(IN) :: stcon(row_length,rows,npft)
        REAL, INTENT(IN) :: soilmc_lp(land_points)
        REAL, INTENT(IN) :: fland(land_points)
        REAL, INTENT(IN) :: laift_lp(land_points,npft)
        REAL, INTENT(IN) :: canhtft_lp(land_points,npft)
        REAL, INTENT(IN) :: z0tile_lp(land_points,ntype)
        REAL, INTENT(IN) :: t0tile_lp(land_points,ntype)
        REAL, INTENT(IN) :: canwctile_lp(land_points,ntype)

        REAL, INTENT(INOUT) :: fastj_dj(row_length,rows,model_levels,   &
                                        jppj)
        REAL, INTENT(INOUT) :: q(row_length,rows,model_levels)   ! water vapo
        REAL, INTENT(INOUT) :: tracer(row_length,rows,                  &
                                      model_levels,ntracers)     ! tracer MMR

        REAL, INTENT(INOUT) :: user_diags(row_length,rows,              &
                                        model_levels,ndiags)     ! chem diags

! SO2 increments
      REAL, INTENT(INOUT) :: delSO2_wet_H2O2(row_length,rows,           &
                                             model_levels)
      REAL, INTENT(INOUT) :: delSO2_wet_O3(row_length,rows,model_levels)
      REAL, INTENT(INOUT) :: delh2so4_chem(row_length,rows,model_levels)
      REAL, INTENT(INOUT) :: delSO2_drydep(row_length,rows,model_levels)
      REAL, INTENT(INOUT) :: delSO2_wetdep(row_length,rows,model_levels)


! Local variables
      INTEGER :: nlev_in_bl(row_length, rows)     ! No levs in bl
      INTEGER :: nlev_in_bl2(theta_field_size)    ! No levs in bl

      INTEGER, SAVE :: nr            ! no of rxns for BE
      INTEGER, SAVE :: n_be_calls    ! no of call to BE solver
      INTEGER, SAVE :: first_row
      INTEGER, SAVE :: first_column

      INTEGER :: i             ! Loop variable
      INTEGER :: j             ! loop variable
      INTEGER :: js            ! loop variable
      INTEGER :: k             ! loop variable
      INTEGER :: l             ! loop variable
      INTEGER :: ll            ! loop variable
      INTEGER :: m             ! loop variable
      INTEGER :: n             ! loop variable
      INTEGER :: sp            ! loop variable
      INTEGER :: n_pnts        ! no. of pts in 2D passed to CDRIVE
      INTEGER, EXTERNAL :: asad_findreaction ! integer function

      INTEGER           :: ierr                     ! Error code: asad diags routines
      INTEGER           :: errcode                  ! Error code: ereport
      CHARACTER(LEN=72) :: cmessage                 ! Error message
      CHARACTER(LEN=10)      :: prods(2)                 ! Products
      LOGICAL           :: blmask(theta_field_size) ! mask

      REAL, PARAMETER :: fxb = 23.45 * Pi_Over_180 ! tropic of capricorn
      REAL, PARAMETER :: fxc = 24.0/ pi
      REAL, SAVE      :: DTS                       ! B. Euler timestep

!       Pressure level above which top boundary conditions are applied
        REAL, PARAMETER :: p_above = 7000.           ! Pa

        REAL :: tgmt               ! GMT time (decimal represent
        REAL :: declin             ! declination
        REAL :: total_water        ! Total water
        REAL :: const              ! constant
        REAL, PARAMETER :: limit = 1200.   ! seconds. For timesteps
                                           ! greater than this we halve
                                           ! the chemical timestep
                                           ! (for IMPACT).

! SO2 increments in molecules/cm^3
      REAL :: SO2_wetox_H2O2(theta_field_size)
      REAL :: SO2_wetox_O3(theta_field_size)
      REAL :: SO2_dryox_OH(theta_field_size)

      REAL, ALLOCATABLE :: BE_rc(:,:)       ! 1-D Rate coeff array
      REAL, ALLOCATABLE :: zfnatr(:,:)      ! 1-D array of non-transported tracers
      REAL, ALLOCATABLE :: ystore(:)        ! array to store H2SO4 when updated in MODE
      REAL :: zftr(theta_field_size,jpctr)  ! 1-D array of tracers
      REAL :: zp  (theta_field_size)        ! 1-D pressure
      REAL :: zt  (theta_field_size)        ! 1-D temperature
      REAL :: zclw(theta_field_size)        ! 1-D cloud liquid water
      REAL :: zfcloud(theta_field_size)     ! 1-D cloud fraction
      REAL :: cdot(theta_field_size,jpctr)  ! 1-D chem. tendency
      REAl :: zq(theta_field_size)          ! 1-D water vapour mmr
      REAL :: pjinda(rows, ntphot, jppj)    ! PJIN at one level
      REAL :: zprt(row_length, rows, jppj)  ! 2-D photolysis rates
      REAL :: zprt1d(theta_field_size,jppj) ! 1-D photolysis rates
      REAL :: tloc  (row_length, rows)      ! local time
      REAL :: daylen(row_length, rows)      ! local daylength
      REAL :: cs_hour_ang(row_length, rows) ! cosine hour angle
      REAL :: tanlat(row_length, rows)      ! tangens of latitude
      REAL :: zdryrt(row_length, rows, jpdd)                ! dry dep rate
      REAL :: zdryrt2(theta_field_size, jpdd)               ! dry dep rate
      REAL :: zwetrt(row_length, rows, model_levels, jpdw)  ! wet dep rate
      REAL :: zwetrt2(theta_field_size, jpdw)               ! wet dep rat
      REAL :: zwetrt3(theta_field_size, model_levels, jpdw) ! wet dep rat
      REAL :: zfrdiss2(theta_field_size,jpdw,jpeq+1)        ! dissolved fraction
      REAL :: zfrdiss(row_length, rows, model_levels, jpdw, jpeq+1)
      REAL :: rc_het(theta_field_size,2)                ! heterog rates for trop chem
      REAL :: kp_nh(row_length, rows, model_levels)     ! Dissociation const
      REAL :: kp_nh2(theta_field_size)                  ! Dissociation const
      REAL :: ozonecol(row_length, rows, model_levels)  ! for strat chem
      REAL :: BE_tnd(theta_field_size)                  ! total no density
      REAL :: BE_h2o(theta_field_size)                  ! water vapour concn
      REAL :: BE_o2(theta_field_size)                   ! oxygen concn
      REAL :: BE_vol(theta_field_size)                  ! gridbox volume
      REAL :: BE_wetrt(theta_field_size,jpspec)         ! wet dep rates (s-1)
      REAL :: BE_dryrt(theta_field_size,jpspec)         ! dry dep rates (s-1)
      REAL :: BE_deprt(theta_field_size,jpspec)         ! dep rates (s-1)
      REAL :: BE_frdiss(theta_field_size,jpspec,jpeq+1) ! dissolved fraction
      REAL :: BE_y (theta_field_size,jpspec)            !
      REAL :: strat_ch4loss(theta_field_size,model_levels) ! for strat ch4 loss
      REAL :: k_dms(theta_field_size,5)                 ! dms rate coeffs
      REAL :: pr2dj(nlphot)                             ! 2D photolysis level pres
      REAL, SAVE :: first_lat
      REAL, SAVE :: dellat
      REAL, SAVE :: first_lon
      REAL, SAVE :: dellon

        LOGICAL, SAVE :: firstcall = .true.

! Variables for heterogeneous chemistry
        REAL, ALLOCATABLE :: shno3_3d(:,:,:)
        LOGICAL :: stratflag(theta_field_size)

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


        IF (lhook) CALL dr_hook('UKCA_CHEMISTRY_CTL',zhook_in,zhook_handle)
        n_pnts = rows * row_length
        nr = nr_therm + nr_phot
        IF (.NOT. ALLOCATED(zfnatr)) ALLOCATE(zfnatr(n_pnts,nnaf)) 
        IF (.NOT. ALLOCATED(ystore) .AND. uph2so4inaer == 1)            &
                                     ALLOCATE(ystore(theta_field_size)) 

        IF (firstcall) THEN

!         Check that theta_field_size = n_pnts

          IF (theta_field_size /= n_pnts) THEN
            cmessage='theta_field_size not equal to n_pnts'
            CALL EREPORT('UKCA_CHEMISTRY_CTL',n_pnts,cmessage)
          END IF

! Determine where the domain is in latitude
          first_lat = ASIN(sinlat(1,1))
          dellat = ASIN(sinlat(1,2)) - first_lat
          first_row = INT((first_lat + 0.5*pi) / dellat + 0.00001) + 1

! Determine where the domain is in longitude
          first_lon = true_longitude(1,2)
          dellon = true_longitude(2,2) - first_lon
          first_column = INT(first_lon / dellon + 0.00001) + 1

!         Update ASAD timestep variables

          dtime = secs_per_step

!         Backward Euler timestep variables

          n_be_calls = INT(dtime/dts0)
! Ensure we call BE at least once
          IF (n_be_calls == 0) n_be_calls = 1
! Calculate the BE Timestep
          dts=dtime/n_be_calls 

          IF ( printstatus >= prstatus_oper ) THEN
            WRITE (6,*) 'n_be_calls, dts= ',n_be_calls, dts
          ENDIF

!         Check whether water vapour is advective tracer. Then,
!         check whether UM and ASAD advective tracers correspond
!         to each other.

          IF ((L_ukca_advh2o) .AND. (n_h2o == 0)) THEN
            cmessage='No tracer for advected water vapour'
            errcode = 4
            CALL EREPORT('UKCA_CHEMISTRY_CTL',errcode,cmessage)
          END IF

! Allocate stratospheric photolysis arrays
          IF (L_ukca_strat .OR. L_ukca_stratcfc .OR. L_ukca_strattrop)  &
            CALL UKCA_STRAT_PHOTOL_INIT()

! Identify the SO2+OH rate coeff, the products are alternatives depending
!  on whether the H2SO4 tracer updating is to be done in ASAD or in MODE.
        iso2_oh = 0
        ih2so4_hv = 0
        IF (L_ukca_achem) THEN
          prods = (/'H2SO4     ','HO2       '/)
          iso2_oh = asad_findreaction( 'SO2       ', 'OH        ',      &
                                 prods, 2, spt, ntrkx, jptk+1, jpspt )

          IF (iso2_oh == 0) THEN   ! check for stratospheric sulphur chemistry
            prods = (/'SO3       ','H2O       '/)
            iso2_oh = asad_findreaction( 'SO2       ', 'OH        ',    &
                                 prods, 2, spt, ntrkx, jptk+1, jpspt )
            prods = (/'SO3       ','OH        '/)
            ih2so4_hv = asad_findreaction( 'H2SO4     ', 'PHOTON    ',  &
                                 prods, 2, spj, nprkx, jppj+1, jpspj )

          END IF

          IF (iso2_oh == 0 .AND. ih2so4_hv == 0) THEN
            cmessage=' Sulphur chemistry reactions not found'
            write(6,*) cmessage
            write(6,*) 'iso2_oh: ',iso2_oh,' ih2so4_hv: ',ih2so4_hv
            errcode = 1
            CALL EREPORT('UKCA_CHEMISTRY_CTL',errcode,cmessage)
          END IF
        END IF   ! L_ukca_achem

      END IF  ! of initialization of chemistry subroutine (firstcall)


!       Calculate local time as function of longitude

        tgmt = real(i_hour) + real(i_minute)/60.0                      &
                            + secs_per_step * 0.5 / 3600.0
!        IF (tgmt < 0.) tgmt = tgmt + 24.

        tloc = tgmt + 24.0 * true_longitude/pi/2.0
        WHERE (tloc > 24.0) tloc = tloc - 24.0

!       Calculate Declination Angle and Daylength for each row for
!       current day of the year.
!       Ensure COS of HOUR ANGLE does not exceed + or - 1, and set DAY
!       LENGTH of 1st & last rows to be the same as that of adjacent rows
!       to avoid possible problems at the poles (tan(90)=infinity).

        daylen = 0.0
        DO i=1,rows
          DO j=1,row_length
            IF (ABS(coslat(j,i)) < 1e-10) THEN
              IF (sinlat(j,i) >= 0.0) THEN
                tanlat(j,i) = 1.0e20
              ELSE
                tanlat(j,i) = -1.0e20
              ENDIF
            ELSE
              tanlat(j,i)=sinlat(j,i)/coslat(j,i)
            END IF
          END DO
        END DO

        declin = fxb * SIN(pi_over_180*(266.0+i_day_number))
        cs_hour_ang = -tanlat * TAN(declin)
        WHERE (cs_hour_ang < -1.0) cs_hour_ang = -1.0
        WHERE (cs_hour_ang >  1.0) cs_hour_ang =  1.0
        daylen = fxc * acos(cs_hour_ang)


! Eliminate spurious negative values from sinlat
        IF (MAXVAL(sinlat) > 0.0 .AND. MINVAL(sinlat) < 0.0) THEN
          WHERE (sinlat < 1e-10) sinlat = 0.0
        ENDIF


!       Call routine to calculate dry deposition rates.


        zdryrt  = 0.0
        zdryrt2 = 0.0
        IF (ndepd /= 0) THEN

          IF (L_ukca_intdd) THEN           ! Call interactive dry dep

! DEPENDS ON: ukca_ddepctl
            CALL UKCA_DDEPCTL(row_length, rows, bl_levels,             &
              land_points, land_index, tile_pts, tile_index,           &
              secs_per_step, sinlat, tile_frac, t_surf,                &
              p_layer_boundaries(:,:,0), dzl, zbl, surf_hf, u_s,       &
              rh, stcon, soilmc_lp, fland, seaice_frac, laift_lp,      &
              canhtft_lp, z0tile_lp, t0tile_lp, canwctile_lp,          &
              nlev_in_bl, zdryrt)

          ELSE                             ! Call prescribed dry dep

! DEPENDS ON: ukca_ddeprt
            CALL UKCA_DDEPRT(daylen, tloc, n_pnts, dzl, bl_levels,     &
                                     z0m, u_s, t_surf,                 &
                                     sinlat, i_month,                  &
                                     1, n_pnts,                        &
                                     zdryrt)

          ENDIF
        ENDIF

!       Call routine to calculate wet deposition rates.

        zwetrt  = 0.0
        zwetrt2 = 0.0
        IF (ndepw /= 0) THEN

! DEPENDS ON: ukca_wdeprt
          CALL UKCA_WDEPRT(drain, crain, n_pnts, model_levels, temp,    &
                           sinlat, secs_per_step, 1,n_pnts, zwetrt)
        ENDIF

! Calculate dissolved fraction
      IF (L_ukca_aerchem) THEN
! DEPENDS ON: ukca_fracdiss
        CALL UKCA_FRACDISS(row_length, rows, model_levels, wet_levels,  &
                         temp, pres, rh, qcl, zfrdiss, kp_nh)
      ENDIF

!       Call routine to read in 2D photolysis rates once per day and
!       interpolate to model latitude and levels.


        IF ( (i_ukca_photol == i_ukca_phot2d) .AND.                     &
             ((i_hour == 0 .AND. i_minute == 0) .OR. firstcall) ) THEN
          CALL UKCA_PHOTIN(i_day_number, row_length, rows,              &
                           model_levels, n_pnts,                        &
                           first_row, global_row_length, global_rows,   &
                           RESHAPE(sinlat,(/theta_field_size/)),        &
                           pres, jppj)
        ELSE IF ((i_ukca_photol == i_ukca_fastj) .OR.                   &
                 (i_ukca_photol == i_ukca_fastjx) ) THEN
          CALL UKCA_INPR2D(pr2d,pr2dj)  ! for stratosphere only
        ENDIF

!       Calculate ozone column for stratospheric photolysis

        IF (L_ukca_strat .OR. L_ukca_stratcfc .OR. L_ukca_strattrop)    &
          THEN
          IF (n_o3 > 0) THEN
            CALL UKCA_CALC_OZONECOL(model_levels, rows, row_length,     &
                           z_top_of_model, p_layer_boundaries, pres,    &
                           tracer(:,:,:,n_o3)/c_o3,                     &
                           ozonecol)
          ENDIF

! Calculate total chlorine and total bromine before chemistry

          IF (nn_cl > 0)  THEN
            CALL UKCA_CONSERVE(row_length, rows, model_levels, ntracers,&
                 tracer, pres, drain, crain, um_ozone3d, .TRUE.)
          ENDIF
        ENDIF    ! L_ukca_strat etc

! if heterogeneous chemistry is selected, allocate solid HNO3 array
      IF (L_ukca_het_psc) THEN
        IF (.NOT. ALLOCATED(shno3_3d)) &
             ALLOCATE(shno3_3d(row_length, rows, model_levels))
        shno3_3d = 0.
      ENDIF

!       Initialize budget variables
      strat_ch4loss          = 0.0

! Wet deposition rates for asad
      DO k=1,model_levels
        DO l=1,jpdw
          zwetrt3(:,k,l) = RESHAPE(zwetrt(:,:,k,l),(/theta_field_size/))
        END DO
      END DO

! Need this line here for ASAD interactive DD.
      nlev_in_bl2(:) = RESHAPE(nlev_in_bl(:,:),(/theta_field_size/))

! Model levels loop
      DO k=1,model_levels

! Copy water vapour and ice field into 1-D arrays
        IF (L_ukca_het_psc) THEN 
          IF (k <= wet_levels) THEN
            sph2o(:) = RESHAPE(qcf(:,:,k),(/theta_field_size/))/c_h2o 
          ELSE
            sph2o(:) = 0.0
          END IF
        END IF 

        zdryrt2(:,:) = 0.0e0
        IF (L_ukca_intdd) THEN
           ! Interactive scheme extracts from levels in boundary layer
           blmask(:) = (k <= nlev_in_bl2(:))
           DO l=1,jpdd
              WHERE (blmask(:))
                 zdryrt2(:,l) = &
                      RESHAPE(zdryrt(:,:,l),(/theta_field_size/))
              ENDWHERE
           END DO
        ELSE    ! non-interactive
           IF (k == 1) THEN
              DO l=1,jpdd
                 zdryrt2(:,l) = &
                      RESHAPE(zdryrt(:,:,l),(/theta_field_size/))
              END DO
           END IF
        END IF

!       Put pressure, temperature and tracer mmr into 1-D arrays
!       for use in ASAD chemical solver

        zp(:) = RESHAPE(pres(:,:,k),(/theta_field_size/))
        zt(:) = RESHAPE(temp(:,:,k),(/theta_field_size/))
        zq(:) = RESHAPE(q(:,:,k),(/theta_field_size/))

        IF (L_UKCA_aerchem) THEN 
          IF (k <= wet_levels) THEN
            zclw(:) = RESHAPE(qcl(:,:,k),(/theta_field_size/))
            zfcloud(:) = RESHAPE(cloud_frac(:,:,k),(/theta_field_size/))
          ELSE
            zclw(:) = 0.0
            zfcloud(:) = 0.0
          END IF 
        END IF 

! Convert mmr into vmr for tracers
        DO js=1,jpctr
          zftr(:,js) = RESHAPE(tracer(:,:,k,js),                        &
                              (/theta_field_size/))/c_species(js)
        END DO

! Update non-advected species
        DO js=1,nnaf 
          zfnatr(:,js) = RESHAPE(user_diags(:,:,k,js),(/n_pnts/))/      &
                         c_na_species(js) 
        ENDDO 

! Calculate photolysis. 
! Calculate TOMCAT-heritage photolysis if PHOT2D is selected.
        IF (i_ukca_photol == i_ukca_phot2d) THEN

!         Interpolate 2D photolysis rates as function of longitude
!         per model step.

          pjinda = pjin(1:rows,k,:,:)
          CALL UKCA_CURVE(pjinda, RESHAPE(tloc,(/n_pnts/)),            &
                          RESHAPE(daylen,(/n_pnts/)), n_pnts,          &
                          model_levels, rows, row_length, zprt1d)
          zprt = RESHAPE(zprt1d,(/row_length, rows, jppj/))

        ELSE
          zprt = 0.
        END IF

! Calculate SLIMCAT-heritage photolysis if PHOT2D is selected or for
! middle-atmosphere chemistry in the mesosphere. STRAT_PHOTOL merges
! the SLIMCAT rates with TOMCAT rates as calculated above.
! If FASTJX is selected, only photolysis for wavelengths less than 175 nm
! is calculated, and the result is added to the FAST-JX rates. This 
! only matters at heights > 60 km as below the short wavelengths can be ignored.

        IF (L_ukca_strat .OR. L_ukca_strattrop .OR. L_ukca_stratcfc)    &
          THEN

          IF ((i_ukca_photol == i_ukca_fastj) .OR.                      &
              (i_ukca_photol == i_ukca_fastjx)) THEN

            IF(minval(pres(:,:,k)) < fastjx_prescutoff .AND.            &
               fastjx_mode /= 3) THEN 
              CALL UKCA_STRAT_PHOTOL(pres(:,:,k), temp(:,:,k),          &
                   ozonecol(:,:,k), cos_zenith_angle,                   &
                   zprt)

! Add up j-rates for lambda < 177 nm with those for lambda >= 177 nm, 
! calculated within FAST-JX.
              IF(fastjx_mode /= 1) THEN
                zprt = zprt + fastj_dj(:,:,k,:)
              ENDIF

              fastj_dj(:,:,k,:) = zprt


            ELSE
! Below cutoff simply use FAST-JX rates. Ignore short wavelengths.
              zprt = fastj_dj(:,:,k,:)
            END IF
          ELSE 
! No FAST-JX here.
            CALL ukca_strat_photol(pres(:,:,k), temp(:,:,k),            &
                 ozonecol(:,:,k), cos_zenith_angle,                     &
                 zprt)
            fastj_dj(:,:,k,:) = zprt 
          END IF

        ELSE 

! tropospheric chemistry selected here. Use FAST-JX rates or previously
! calculated 2-D rates.

          IF ((i_ukca_photol == i_ukca_fastj) .OR.                      &
              (i_ukca_photol == i_ukca_fastjx)) THEN
            zprt = fastj_dj(:,:,k,:)
          END IF
        END IF
      
! Map photolysis rates onto 1-D array.
        DO l=1,jppj
          zprt1d(:,l) = RESHAPE(zprt(:,:,l),(/theta_field_size/))
        END DO

!       Call ASAD routines to do chemistry integration
!       In lowest levels choose half the dynamical timestep for
!       chemistry. If dynamical timestep > 20 min, use half and
!       quarter of dynamical timestep for chemistry.

        IF (.NOT.(L_ukca_trop .OR. L_ukca_aerchem .OR.                  &
                  L_ukca_raq)) THEN                      ! Not B-E

! fill stratospheric flag indicator and SO4 surface area
! retrieve tropospheric heterogenous rates from previous time step
          IF(L_ukca_trophet) THEN
            rc_het(:,1) = RESHAPE(user_diags(:,:,k,10),                 &
                                  (/theta_field_size/))       ! N2O5
            rc_het(:,2) = RESHAPE(user_diags(:,:,k,11),                 &
                                  (/theta_field_size/))       ! HO2+HO2
          ENDIF
          stratflag(:) = (.NOT. RESHAPE(L_troposphere(:,:,k),           &
                                    (/theta_field_size/)) )
          za(:) = RESHAPE(so4_sa(:,:,k),(/theta_field_size/))

          IF (method == 3) THEN
! timestep for Newton-Raphson solver: 1 hour. Note that in this case
! UKCA_CHEMISTRY_CTL is only called every 2nd/3rd dynamical timestep
! (for a 30 / 20 minutes stepsize), and dtime should equal 1 h.
! In case of non-convergence Newton-Raphson contains automatic stepsize
! refinement.
! secs_per_step is the *CHEMISTRY* timestep and not the model timestep, 
! and equals INTERVAL*TIMESTEP, where INTERVAL=3600/INT(TIMESTEP) and 
! determines when chemistry is called 
            dtime = secs_per_step !*2  ! changed 2010/02/08 20:40 
            ncsteps = 1
            cdt = 3600.0 !dtime
          ELSEIF (method == 1) THEN
! use about 15 or 10 minutes, depending on dynamical timestep
            IF (secs_per_step < limit) THEN
              cdt = dtime
              ncsteps = 1
            ELSE
              cdt = 0.5*dtime
              ncsteps = 2
            END IF
          END IF


          IF (uph2so4inaer == 1) THEN 
! H2SO4 will be updated in MODE, so store old value here 
            ystore(:) = y(:,nn_h2so4) 
          END IF 

! DEPENDS ON: asad_cdrive
          CALL ASAD_CDRIVE(cdot, zftr, zp, zt, zq,                      &
                           RESHAPE(cloud_frac(:,:,k),(/n_pnts/)),       &
                           RESHAPE(qcl(:,:,k),(/n_pnts/)),              &
                           k,zdryrt2, zwetrt3, rc_het,                  &
                           zprt1d, n_pnts, stratflag)

          IF (L_ukca_het_psc) THEN
! Save MMR of NAT PSC particles into 3-D array for PSC sedimentation.
! Note that sphno3 is NAT in number density of HNO3.
            IF (ANY(sphno3(1:theta_field_size) > 0.)) THEN
              shno3_3d(:,:,k) = RESHAPE(sphno3(:)/tnd(:),               &
                                   (/row_length,rows/))*c_hono2
            ELSE
              shno3_3d(:,:,k) = 0.
            END IF
          END IF

          ! copy over O(1D) into chem_diags
          IF (LVMR) THEN
             user_diags(1:row_length,1:rows,k,1) = &
                  RESHAPE(y(:,nn_o1d)/tnd(:),(/row_length,rows/))*c_o1d
          ELSE
             user_diags(1:row_length,1:rows,k,1) = &
                  RESHAPE(y(:,nn_o1d),(/row_length,rows/))*c_o1d
          END IF


          IF (L_ukca_chem .AND. L_ukca_achem) THEN
! Calculate chemical fluxes for MODE
            IF (ihso3_h2o2 > 0) delSO2_wet_H2O2(:,:,k) =                &
              delSO2_wet_H2O2(:,:,k) + RESHAPE(rk(:,ihso3_h2o2)*        &
              y(:,nn_so2)*y(:,nn_h2o2),(/row_length,rows/))*cdt
            IF (ihso3_o3 > 0) delSO2_wet_O3(:,:,k) =                    &
              delSO2_wet_O3(:,:,k) + RESHAPE(rk(:,ihso3_o3)*            &
              y(:,nn_so2)*y(:,nn_o3),(/row_length,rows/))*cdt
            IF (iso3_o3 > 0) delSO2_wet_O3(:,:,k) =                     &
              delSO2_wet_O3(:,:,k) + RESHAPE(rk(:,iso3_o3)*             &
              y(:,nn_so2)*y(:,nn_o3),(/row_length,rows/))*cdt
            IF (iso2_oh > 0 .AND. ih2so4_hv > 0) THEN  ! net H2SO4 production
              delh2so4_chem(:,:,k) = delh2so4_chem(:,:,k) +             &
               (RESHAPE(rk(:,iso2_oh)*y(:,nn_so2)*y(:,nn_oh),           &
               (/row_length,rows/)) - RESHAPE(rk(:,ih2so4_hv)*          &
                y(:,nn_h2so4),(/row_length,rows/)))*cdt
            ELSE IF (iso2_oh > 0) THEN
              delh2so4_chem(:,:,k) = delh2so4_chem(:,:,k) +             &
               RESHAPE(rk(:,iso2_oh)*y(:,nn_so2)*y(:,nn_oh),            &
               (/row_length,rows/))*cdt
            END IF
! Restore H2SO4 tracer as it will be updated in MODE using delh2so4_chem
            IF (uph2so4inaer == 1) y(:,nn_h2so4) = ystore(:) 
          END IF

! 3D flux diagnostics
          IF (L_asad_use_chem_diags .AND.                               &
               ((L_asad_use_flux_rxns .OR. L_asad_use_rxn_rates) .OR.   &
               (L_asad_use_wetdep .OR. L_asad_use_drydep)))             &
               CALL ASAD_CHEMICAL_DIAGNOSTICS(row_length,rows,          &
               model_levels,k,secs_per_step,volume,ierr) 

! PSC diagnostics
          IF (L_asad_use_chem_diags .AND. L_asad_use_psc_diagnostic)    &
               CALL ASAD_PSC_DIAGNOSTIC(row_length,rows,k,ierr) 

! Bring results back from vmr to mmr. Treat water
! vapour separately if necessary.

          IF (L_ukca_h2o_feedback) THEN
            DO l=1,jpctr
              tracer(:,:,k,l) = RESHAPE(zftr(:,l),(/row_length,rows/))  &
                                        *c_species(l)
            END DO
            q(:,:,k) = RESHAPE(zftr(:,n_h2o),(/row_length,rows/))*c_h2o
          ELSE
            IF (n_h2o < jpctr) THEN
              DO l=n_h2o+1,jpctr
                tracer(:,:,k,l) = RESHAPE(zftr(:,l),                    &
                                (/row_length,rows/))*c_species(l)
              END DO
            END IF

            IF (n_h2o > 1) THEN
              DO l=1,n_h2o-1
                tracer(:,:,k,l) = RESHAPE(zftr(:,l),                    &
                                (/row_length,rows/))*c_species(l)
              END DO
            END IF
          END IF

! Set SS species concentrations for output (stratospheric configurations)
          IF (O1D_in_ss) user_diags(:,:,k,1) =                          &
            RESHAPE(y(:,nn_o1d)/tnd(:),(/row_length,rows/))*            &
            c_species(nn_o1d)                              ! O1D mmr
          IF (O3P_in_ss) user_diags(:,:,k,2) =                          &
            RESHAPE(y(:,nn_o3p)/tnd(:),(/row_length,rows/))*            &
            c_species(nn_o3p)                              ! O3P mmr

        ELSE  ! Backward Euler with non-families

!         Calculate total number  density, o2, h2o, and tracer
!         concentrations for Backward Euler solver

          const = Boltzmann*1.0e6
          nlev_in_bl2(:) = RESHAPE(nlev_in_bl(:,:),(/theta_field_size/))
          BE_tnd(:) = zp(:) / (const * zt(:) )
          BE_o2(:)  = 0.2095 * BE_tnd(:)
          BE_h2o(:) = zq(:)*BE_tnd(:)
          BE_vol(:) = RESHAPE(volume(:,:,k),(/theta_field_size/))*1.0e6 !  m3->cm3
          DO l=1,jpdw
            zwetrt2(:,l) = RESHAPE(zwetrt(:,:,k,l),(/theta_field_size/))
          END DO
          DO l=1,jpdd
            zdryrt2(:,l) = RESHAPE(zdryrt(:,:,l),(/theta_field_size/))
          END DO
          DO l=1,jpctr
            zftr(:,l) = zftr(:,l) * BE_tnd(:)
          END DO

! Update non-advected species (only for B-E solvers)
          IF (.NOT. ALLOCATED(zfnatr)) ALLOCATE(zfnatr(n_pnts,nnaf)) 
          DO js=1,nnaf 
            zfnatr(:,js) = RESHAPE(user_diags(:,:,k,js),(/n_pnts/))/    &
                           c_na_species(js) 
          ENDDO 
          DO l=1,nnaf
            zfnatr(:,l) = zfnatr(:,l) * BE_tnd(:)
          END DO

          IF (L_ukca_aerchem) THEN 
            DO l=1,jpdw
              DO j=1,jpeq+1
                zfrdiss2(:,l,j) = RESHAPE(zfrdiss(:,:,k,l,j),           &
                                  (/theta_field_size/)) 
              END DO
            END DO
            kp_nh2(:) = RESHAPE(kp_nh(:,:,k),(/theta_field_size/)) 
          END IF

!         Assign wet and dry deposition rates to species

! DEPENDS ON: ukca_be_wetdep
          CALL UKCA_BE_WETDEP(n_pnts, zwetrt2, be_wetrt)

! DEPENDS ON: ukca_be_drydep
          CALL UKCA_BE_DRYDEP(k, n_pnts, nlev_in_bl2, zdryrt2, be_dryrt) 
 
! Assign fractional dissociation
          IF(L_ukca_aerchem) THEN
            DO i=1,jpeq+1
! DEPENDS ON: ukca_be_wetdep
              CALL UKCA_BE_WETDEP(n_pnts, zfrdiss2(:,:,i),              &
                                  be_frdiss(:,:,i))
            END DO
          ENDIF

!         Calculate reaction rate coefficients

          IF (.NOT. ALLOCATED(BE_rc)) &
               ALLOCATE(BE_rc(theta_field_size,nr_therm))
          
          IF (L_ukca_raq) THEN
! DEPENDS ON: ukca_chemco_raq
            CALL UKCA_CHEMCO_RAQ(nr_therm, n_pnts, zt(1:n_pnts),        &
                            BE_tnd(1:n_pnts), BE_h2o(1:n_pnts),         &
                            BE_o2(1:n_pnts), BE_rc(1:n_pnts,:))
          ELSE
  
! DEPENDS ON: ukca_chemco
          CALL UKCA_CHEMCO(nr_therm, n_pnts, zt(1:n_pnts),              &
                           BE_tnd(1:n_pnts), BE_h2o(1:n_pnts),          &
                           BE_o2(1:n_pnts), zclw(1:n_pnts),             &
                           zfcloud(1:n_pnts), BE_frdiss(1:n_pnts,:,:),  &
                           k_dms(1:n_pnts,:), BE_rc(1:n_pnts,:))
          END IF

!         Assign tracer concentrations to species concentrations

          BE_y(1:n_pnts,1:jpspec) = 0.0
          DO i = 1,jpspec
            DO j = 1,jpctr
              IF (speci(i) == advt(j)) THEN
                BE_y(1:n_pnts,i) = zftr(1:n_pnts,j)
                EXIT
              ENDIF
            ENDDO
          ENDDO

! Assign non-advected concentrations to species concentrations
          DO i = 1,jpspec 
            DO j = 1,nnaf 
              IF (speci(i) == nadvt(j)) THEN 
                BE_y(1:n_pnts,i) = zfnatr(1:n_pnts,j)  
                EXIT 
              ENDIF 
            ENDDO 
          ENDDO 

          IF (L_ukca_aerchem) THEN
            SO2_wetox_H2O2(:) = 0.0
            SO2_dryox_OH(:)   = 0.0
            SO2_wetox_O3(:)   = 0.0
          ELSE
            delSO2_wet_H2O2(:,:,:) = 0.0
            delSO2_wet_O3(:,:,:)   = 0.0
            delh2so4_chem(:,:,:)   = 0.0
            delSO2_drydep(:,:,:)   = 0.0
            delSO2_wetdep(:,:,:)   = 0.0
          ENDIF

!         Call Backward Euler solver
!         N.B. Emissions already added, via call to TR_MIX from
!         UKCA_EMISSION_CTL

          IF (L_ukca_aerchem) THEN

! DEPENDS ON: ukca_deriv_aero
            CALL UKCA_DERIV_AERO(nr_therm, nr_phot, n_be_calls,         &
                          BE_rc(1:n_pnts,:),                            &
                          BE_wetrt(1:n_pnts,:), BE_dryrt(1:n_pnts,:),   &
                          zprt1d(1:n_pnts,:), k_dms(1:n_pnts,:),        &
                          BE_h2o(1:n_pnts), BE_tnd(1:n_pnts),           &
                          BE_o2(1:n_pnts),                              &
                          dts, BE_y(1:n_pnts,:),                        &
                          SO2_wetox_H2O2(1:n_pnts),                     &
                          SO2_wetox_O3(1:n_pnts),                       &
                          SO2_dryox_OH(1:n_pnts) )
          ELSE IF (L_ukca_raq) THEN

! DEPENDS ON: ukca_deriv_raq
            CALL UKCA_DERIV_RAQ(nr_therm, n_be_calls,                   &
                       n_pnts,                                          &
                       BE_rc(1:n_pnts,:), BE_wetrt(1:n_pnts,:),         &
                       BE_dryrt(1:n_pnts,:), zprt1d(1:n_pnts,:),        &
                       BE_h2o(1:n_pnts), BE_tnd(1:n_pnts),              &
                       BE_o2(1:n_pnts),                                 &
                       dts, BE_y(1:n_pnts,:) )
          ELSE

! DEPENDS ON: ukca_deriv
            CALL UKCA_DERIV(nr, n_be_calls, n_pnts,                    &
                       BE_rc(1:n_pnts,:), BE_wetrt(1:n_pnts,:),        &
                       BE_dryrt(1:n_pnts,:), zprt1d(1:n_pnts,:),       &
                       BE_h2o(1:n_pnts), BE_tnd(1:n_pnts),             &
                       BE_o2(1:n_pnts),                                &
                       dts, BE_y(1:n_pnts,:) )
          ENDIF
          IF (ALLOCATED(BE_rc)) DEALLOCATE(BE_rc)

          IF (k == model_levels .OR. k == model_levels-1 .OR.          &
              k == model_levels-2 .AND. (.NOT. L_ukca_use_2dtop)) THEN
! DEPENDS ON: ukca_ch4_stratloss
            CALL UKCA_CH4_STRATLOSS(n_be_calls, n_pnts,                &
                       BE_vol(1:n_pnts), dts,                          &
                       BE_y(1:n_pnts,nn_ch4), strat_ch4loss(1:n_pnts,k))
          ENDIF

          DO j = 1,jpctr
            DO i = 1,jpspec
              IF (advt(j) == speci(i)) THEN
                zftr(1:n_pnts,j) = BE_y(1:n_pnts,i)
                EXIT
              ENDIF
            ENDDO
          ENDDO

          DO js = 1, jpctr
            tracer(:,:,k,js) = RESHAPE(zftr(:,js)/BE_tnd(:),           &
                                   (/row_length,rows/))*c_species(js)
          END DO

! Convert non-advected tracers back to mmr
          DO j = 1,nnaf
            DO i = 1,jpspec
              IF (nadvt(j) == speci(i)) THEN
                zfnatr(:,j) = BE_y(1:n_pnts,i)/BE_tnd(1:n_pnts)
                EXIT
              END IF
            ENDDO
          ENDDO

          IF (n_chem_diags >= nnaf) THEN
            DO js = 1, nnaf
              user_diags(:,:,k,js)=RESHAPE(zfnatr(:,js),                &
                                  (/row_length,rows/))*c_na_species(js)
            END DO
          ELSE
            cmessage=' Not enough space in user_diags for '//           &
                     'non-advected tracers'
            errcode = 1
            CALL EREPORT('UKCA_CHEMISTRY_CTL',errcode,cmessage)
          END IF

          IF (n_chem_diags >= 22) THEN
            IF (.NOT.(L_ukca_raq)) THEN    
! trop CH4 burden in mol
              user_diags(:,:,k,16) = RESHAPE(zftr(:,n_ch4),             &
                                  (/row_length,rows/))*volume(:,:,k)*   &
                                   1.0e6/avogadro
! trop O3 burden in mol
              user_diags(:,:,k,17) = RESHAPE(zftr(:,n_o3),              &
                                  (/row_length,rows/))*volume(:,:,k)*   &
                                   1.0e6/avogadro
! trop OH burden in mol
              user_diags(:,:,k,18) = RESHAPE(BE_y(:,nn_oh),             &
                                  (/row_length,rows/))*volume(:,:,k)*   &
                                   1.0e6/avogadro
! strat CH4 burden in mol
              user_diags(:,:,k,19) = user_diags(:,:,k,16)
! pv-theta surface
              user_diags(:,:,k,20) = p_tropopause(:,:)      ! pv-theta surface
! strat ch4-oh rxn flux
              user_diags(:,:,k,21) = 0.0
! strat ch4 loss
              user_diags(:,:,k,22) = RESHAPE(strat_ch4loss(:,k),        &
                                  (/row_length,rows/))

              WHERE (L_troposphere(:,:,:))
                user_diags(:,:,:,19) = 0.0   ! strat ch4 burden in mol
                user_diags(:,:,:,21) = 0.0   ! strat ch4-oh rxn flux
              ELSEWHERE
                user_diags(:,:,:,16) = 0.0   ! trop ch4 burden in mol
                user_diags(:,:,:,17) = 0.0   ! trop o3  burden in mol
                user_diags(:,:,:,18) = 0.0   ! trop oh  burden in mol
              ENDWHERE
          
            ELSE   ! There are 18 non-advected species in RAQ chemistry
! trop CH4 burden in mol
              user_diags(:,:,k,19) = RESHAPE(zftr(:,n_ch4),             &
                                  (/row_length,rows/))*volume(:,:,k)*   &
                                   1.0e6/avogadro
! trop O3 burden in mol
              user_diags(:,:,k,20) = RESHAPE(zftr(:,n_o3),              &
                                  (/row_length,rows/))*volume(:,:,k)*   &
                                   1.0e6/avogadro
! trop OH burden in mol
              user_diags(:,:,k,21) = RESHAPE(BE_y(:,nn_oh),             &
                                  (/row_length,rows/))*volume(:,:,k)*   &
                                   1.0e6/avogadro
! pv-theta surface
              user_diags(:,:,k,22) = p_tropopause(:,:)   ! pv-theta surface
            
              WHERE (.NOT. L_troposphere(:,:,:))
                user_diags(:,:,:,19) = 0.0   ! trop ch4 burden in mol
                user_diags(:,:,:,20) = 0.0   ! trop o3  burden in mol
                user_diags(:,:,:,21) = 0.0   ! trop oh  burden in mol
              ENDWHERE
            ENDIF     ! L_ukca_RAQ
          ELSE
            cmessage=' Not enough space for tropospheric diags. '//     &
                     'in user_diags array.'
            errcode = 1
            CALL EREPORT('UKCA_CHEMISTRY_CTL',errcode,cmessage)
          END IF        ! n_chem_diags >= 22
    
          IF (L_ukca_aerchem) THEN
            delSO2_wet_H2O2(:,:,k)=RESHAPE(SO2_wetox_H2O2(:),           &
                                   (/row_length,rows/))
            delSO2_wet_O3(:,:,k)=RESHAPE(SO2_wetox_O3(:),               &
                                   (/row_length,rows/))
            delh2so4_chem(:,:,k)=RESHAPE(SO2_dryox_OH(:),               &
                                   (/row_length,rows/))
            delSO2_drydep(:,:,k) = RESHAPE(BE_dryrt(:,nn_so2)*          &
                        BE_y(:,nn_so2),(/row_length,rows/)) * DTS
            delSO2_wetdep(:,:,k) = RESHAPE(BE_wetrt(:,nn_so2)*          &
                        BE_y(:,nn_so2),(/row_length,rows/)) * DTS
          ENDIF

        END IF
      END DO            ! level loop (k)

! Rescale bromine and chlorine tracers to guarantee conservation of total
! chlorine, bromine, and hydrogen over timestep. Only makes sense if at least
! chlorine chemistry is present.

          IF (nn_cl > 0) THEN
            CALL ukca_conserve(row_length, rows, model_levels, ntracers,&
                 tracer, pres, drain, crain, um_ozone3d,.false.)
          END IF

      IF (L_ukca_strat .OR. L_ukca_stratcfc .OR. L_ukca_strattrop)      &
         THEN
        IF (L_ukca_het_psc) THEN
! Do NAT PSC sedimentation

! take NAT out of gasphase again
          tracer(:,:,:,n_hono2) = tracer(:,:,:,n_hono2) - shno3_3d

! DEPENDS ON: ukca_sediment
          CALL ukca_sediment(rows, row_length, model_levels, shno3_3d,  &
                   qcf, r_theta_levels, mass, secs_per_step)

! add solid-phase HNO3 back to gasphase HNO3
          tracer(:,:,:,n_hono2) = tracer(:,:,:,n_hono2) + shno3_3d
        END IF

! Again calculate ozone column, this time from actual post-chemistry
! ozone, for diagnostic purposes
        IF (L_ukca_strat .OR. L_ukca_strattrop .OR. L_ukca_stratcfc)    &
          THEN
          IF (n_o3 > 0) THEN
            CALL UKCA_CALC_OZONECOL(model_levels, rows, row_length,     &
                               z_top_of_model,p_layer_boundaries,pres,  &
                               tracer(:,:,:,n_o3)/c_o3,                 &
                               user_diags(:,:,:,icd_O3col))
            ! convert to Dobson Units from molecules/cm2
            user_diags(:,:,:,icd_O3col) = user_diags(:,:,:,icd_O3col)/  &
                                          2.685E16
          END IF   ! n_o3

! Copy NAT MMR into user_diagostics
          IF (L_ukca_het_psc) user_diags(:,:,:,icd_NAT)=shno3_3d(:,:,:)
        END IF

! Tracer overwrites required to stop accumulation of tracer mass in the uppermost layers
        DO i=1,rows
          DO j=1,row_length
            tracer(j,i,model_levels  ,:) = tracer(j,i,model_levels-2,:)
            tracer(j,i,model_levels-1,:) = tracer(j,i,model_levels-2,:)
          END DO
        END DO

      ELSE     ! tropospheric chemistry
! Call routine to overwrite O3, CH4 and NOy species once per day
! above level defined by p_above. Only for tropospheric chemistry

! DEPENDS ON: ukca_stratf
        CALL UKCA_STRATF(i_day_number, row_length,rows, model_levels,  &
                      n_pnts, first_row, global_row_length,            &
                      global_rows, jpctr, sinlat, pres,                &
                      um_ozone3d, p_above,                             &
                      tracer(1:row_length,1:rows,1:model_levels,       &
                             1:jpctr))

      END IF     ! L_ukca_strat etc


      IF (L_ukca_het_psc .AND. ALLOCATED(shno3_3d)) DEALLOCATE(shno3_3d)
      IF (ALLOCATED(zfnatr)) DEALLOCATE(zfnatr)
      IF (ALLOCATED(ystore)) DEALLOCATE(ystore)

      IF (firstcall) firstcall = .false.


      IF (lhook) CALL dr_hook('UKCA_CHEMISTRY_CTL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE UKCA_CHEMISTRY_CTL
