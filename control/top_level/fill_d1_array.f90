! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  SUBROUTINE FILL_D1_ARRAY------------------------------------------
!LL
!LL  PURPOSE: Fill D1 addressing array with useful information.
!
!             Code Owner: See Unified Model Code Owners HTML page
!             This file belongs in section: Top Level

      SUBROUTINE FILL_D1_ARRAY(                                         &
! ARGSTS Applicable to all configurations. STASH related variables for
! describing output requests and space management.
     & SF,STINDEX,STLIST,SI,STTABL,STASH_MAXLEN,PPINDEX,STASH_LEVELS,        &
     & STASH_PSEUDO_LEVELS,STASH_SERIES,STASH_SERIES_INDEX,SF_CALC,          &
! ARGSTS end
! ARGD1 start
      ! IN/OUT:Addressing of D1 & D1 array
     &  D1_ADDR,D1,LD1,ID1,                                             &
! ARGD1 end
     &                  ICODE,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE PrintStatus_mod
      USE UM_ParParams 
      USE cppxref_mod, ONLY: ppx_grid_type, ppx_halo_type
      USE version_mod, ONLY:                                            &
          nitemp, nelemp, nrecdp, nsectp, nlevp_s, nlevlstsp,           &
          nproftp, nprofdp, nprofup, ndiagpm, ntimep, NTimSerP,         &
          nlevp, npslevp, npslistp, outfile_s, outfile_e

      USE Submodel_Mod  
      USE chsunits_mod, ONLY : nunits
      USE stextend_mod, ONLY: IN_S, D1_PADDR, N_OBJ_D1,             &
                  d1_type, d1_im, d1_extra_info, d1_levs, d1_sect,  &
                  prog, diag, seco, extra_d1, levlst_s, lenplst

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
! --------------------- Comdeck: CHISTORY ----------------------------
!
!  Purpose: COMMON block for history data needed by top level (C0)
!           routines, and passed from run to run.  Mostly set by
!           the User Interface.
!
!           Note that CHISTORY *CALLs ALL individual history comdecks
!
! --------------------------------------------------------------------
!
! ----------------------- Comdeck: IHISTO   ----------------------------
! Description: COMDECK defining Integer History variables for the
!              model overall.
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

!   Type declarations


      ! Array containing model data time (Same as MODEL_BASIS_TIME/MODEL
      ! ANALYSIS_HRS depending whether before/after assimilation)
      INTEGER :: model_data_time(6)

      ! Indicator of operational run type
      INTEGER :: run_indic_op

      ! Final target date for the run
      INTEGER :: run_resubmit_target(6)

      ! Last field written/read per FT unit
      INTEGER :: ft_lastfield(20:nunits)

      ! Number of automatically-resubmitted job chunks
      ! Used to name output file
      INTEGER :: run_job_counter

! History Common Block for overall model integers variables.

      COMMON /ihisto/                                                 &
         model_data_time,                                             &
         run_indic_op, run_job_counter,                               &
         run_resubmit_target, ft_lastfield

      NAMELIST /nlihisto/                                             &
         model_data_time,                                             &
         run_indic_op, run_job_counter,                               &
         run_resubmit_target, ft_lastfield

! IHISTO end
! ----------------------- Comdeck: CHISTO   ----------------------------
! Description: COMDECK defining Character History variables for the
!              model overall.
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

  CHARACTER(LEN=8) ::  run_type             ! Type of run
  CHARACTER(LEN=1) ::  ft_active(20:nunits) ! "Y" if file partly written

  LOGICAL :: newrun ! Set to true in NRUN to stop auto-resubmission

  ! History Common Block for overall model character variables.

  COMMON /chisto/                                     &
     run_type,                                        &
     newrun, ft_active

  NAMELIST /nlchisto/                                 &
     run_type,                                        &
     ft_active

! CHISTO end
! ----------------------- Comdeck: IHISTG   ----------------------------
! Description: COMDECK defining Integer History variables for
!              generic aspects of internal models
!              Generic means values likely to be common to the control
!              of any sub-model/internal model.
!
! This file belongs in section: Top Level
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 

!   Type declarations
      ! History block copy of A_STEP held in file CTIME
      INTEGER :: h_stepim(n_internal_model_max)

      ! No of means activated
      INTEGER :: mean_offsetim(n_internal_model_max)

      ! Offset between MEAN_REFTIME and model basis time(in model dumps)
      INTEGER :: offset_dumpsim(n_internal_model_max)

      ! No of mean periods chosen
      INTEGER :: mean_numberim(n_internal_model_max)

      ! Indicators used to correct logical units are used for
      ! atmos partial sum dump I/O
      INTEGER :: run_meanctl_indicim(4,n_internal_model_max)

      ! History Common Block for generic model integer variables.

      COMMON /ihistg/                                         &
         h_stepim, mean_offsetim, offset_dumpsim,             &
         mean_numberim, run_meanctl_indicim

      NAMELIST /nlihistg/                                     &
         h_stepim, mean_offsetim, offset_dumpsim,             &
         run_meanctl_indicim

! IHISTG end
! ----------------------- Comdeck: CHISTG   ----------------------------
! Description: COMDECK defining Character variables for
!              managing dump names
!
! This file belongs in section: common
!
! Code description: 
!  Language: Fortran 95. 
!  This code is written to UMDP3 standards. 
!
!   Type declarations
!
! For keeping old restart dump name between creation of new restart dump
! and successful completion of climate means and history file.
CHARACTER(LEN=256) :: save_dumpname_im(n_internal_model_max)
! Name of current restart dump
CHARACTER(LEN=256) :: checkpoint_dump_im(n_internal_model_max)
! Blank name
CHARACTER(LEN=256) :: blank_file_name
!
! History Common Block for generic model characters variables.
!
COMMON /chistg/save_dumpname_im, checkpoint_dump_im, blank_file_name

NAMELIST /nlchistg/checkpoint_dump_im

! CHISTG end
!
!  Purpose: Defines unit numbers relevant to history file
!           and variables used to hold the logical to physical
!           file associations made within the model
!
!  Logical Filenames used in the model
!
      CHARACTER(LEN=256) hkfile,ppxref,config,stashctl,namelist,output,      &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,ftxx,    &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar
     

!
      CHARACTER(LEN=256) MODEL_FT_UNIT ! Array holding FORTRAN unit file
!                                 ! associations details for each unit
!
      INTEGER                                                           &
     &        MCTL_UNIT,                                                &
                                 ! Master control namelist file unit
     &        ICTL_UNIT,                                                &
                                 ! Interim control namelist file unit
     &        RSUB_UNIT,                                                &
                                 ! File indicating whether resub required
     &        XHIST_UNIT,                                               &
                                 ! Main history file unit
     &        THIST_UNIT,                                               &
                                 ! Backup history file unit
     &        HKFILE_UNIT,                                              &
                                 ! Operational houskeeping file unit    
     &        EG_UNIT            ! ENDGame diagnostics/info unit
!
! Parameters specifying unit numbers relevant to control/history tasks
!
      PARAMETER(HKFILE_UNIT= 1)
      PARAMETER(MCTL_UNIT  = 8)
      PARAMETER(ICTL_UNIT  = 9)
      PARAMETER(RSUB_UNIT =10)
      PARAMETER(XHIST_UNIT =11)
      PARAMETER(THIST_UNIT =12)

!
! Parameters specifying unit numbers relevant to ENDGame diagnostics
!
      PARAMETER(EG_UNIT  = 55)

! UKCA unit numbers

      INTEGER, PARAMETER :: ukcafjxx_unit=170 ! Fast-J(X) cross section data
      INTEGER, PARAMETER :: ukcafjsc_unit=171 ! Fast-JX scattering data
      INTEGER, PARAMETER :: ukca2do3_unit=172 ! 2D top boundary O3 data 
      INTEGER, PARAMETER :: ukca2ch4_unit=173 ! 2D top boundary CH4 data
      INTEGER, PARAMETER :: ukca2noy_unit=174 ! 2D top boundary NOY data
      INTEGER, PARAMETER :: ukca2pho_unit=175 ! 2D photolysis input data
      INTEGER, PARAMETER :: ukcastrd_unit=176 ! Stratospheric model radiation field. 
      INTEGER, PARAMETER :: ukcasto3_unit=177 ! Strat standard atmosphere T and O3.
      INTEGER, PARAMETER :: ukcastar_unit=178 ! Stratospheric sulfate aerosol climatology 
      INTEGER, PARAMETER :: ukcafjar_unit=179 ! Sulfate aerosol cliamtology for Fast-JX
! Text output file for STASH-related information is assigned to UNIT 200

!
! Namelist of all permissible logical files.
!
      NAMELIST / nlcfiles /                                             &
                   hkfile,ppxref,config,stashctl,namelist,output,       &
                   output2,mctl,ictl,rsub,xhist,thist,icecalve,         &
                   cache1,cache2,aswap,oswap,                           &
                   ainitial,astart,arestart,aopsum1,aopsum2,aopsum3,    &
                   aopsum4,aomean,ssu,                                  &
                   ozone,smcsnowd,dsoiltmp,soiltype,genland,sstin,      &
                   sicein,perturb,mask,                                 &
                   oinitial,ostart,orestart,aopstmp1,aopstmp2,aopstmp3, &
                   aopstmp4,                                            &
                   wfin,hfluxin,pmein,icefin,airtmp,                    &
                   swspectd,                                            &
                   pp0,pp1,pp2,pp3,pp4,pp5,pp6,pp7,pp8,pp9,             &
                   ppvar,pp10,                                          &
                   obs01,obs02,obs03,obs04,obs05,                       &
                   dustsoil,biomass,rivstor,rivchan,river2a,            &
                   surfemis, aircrems, stratems, extraems, radonems,    &
                   lwspectd,surgeou1,surgeout,ppscreen,ppsmc,wfout,     &
                   uarsout1,uarsout2,icefout,mosout,sstout,siceout,     &
                   curntout,flxcrout,dmsconc,orog,olabcin,ocndepth,     &
                   curntin,fluxcorr,slabhcon,atmanl,ocnanl,bas_ind,     &
                   transp,atracer,otracer,sulpemis,usrancil,usrmulti,   &
                   ousrancl,ousrmult,murkfile,                          &
                   alabcin1,alabcin2,                                   &
                   alabcou1,alabcou2,alabcou3,alabcou4,                 &
                   alabcou5,alabcou6,alabcou7,alabcou8,cariolo3,        &
                   foamout1,foamout2,cxbkgerr,rfmout,                   &
                   wlabcou1,wlabcou2,wlabcou3,wlabcou4,horzgrid,        &
                   tdf_dump,iau_inc,                                    &
                   landfrac,                                            &
                   so2natem,chemoxid,aerofcg,fracinit,veginit,disturb,  &
                   cached,sootemis,                                     &
                   co2emits,tppsozon,                                   &
                   vert_lev,var_grid,                                   &
                   idealise,icfile,                                     &
                   arclbiog,arclbiom,arclblck,arclsslt,arclsulp,        &
                   arcldust,arclocff,arcldlta,rpseed,ocffemis,          &
                   topmean,topstdev,ppmbc,                              &
                   ukcaprec,ukcaacsw,ukcaaclw,ukcacrsw,ukcacrlw,        &
                   ukcafjxx,ukcafjsc,ukca2do3,ukca2ch4,ukca2noy,        &
                   ukca2pho,ukcastrd,ukcasto3,ukcastar,ukcafjar

!
!Common block definition
!
      COMMON/CLFHIST/MODEL_FT_UNIT(NUNITS)
!
! Equivalence logical filenames within array MODEL_FT_UNIT
!
      EQUIVALENCE                                                       &
     &(HKFILE    ,MODEL_FT_UNIT(1)  ),(PPXREF     ,MODEL_FT_UNIT(2)  ), &
     &(CONFIG    ,MODEL_FT_UNIT(3)  ),(STASHCTL   ,MODEL_FT_UNIT(4)  ), &
     &(NAMELIST  ,MODEL_FT_UNIT(5)  ),(OUTPUT     ,MODEL_FT_UNIT(6)  ), &
     &(OUTPUT2   ,MODEL_FT_UNIT(7)  ),(MCTL       ,MODEL_FT_UNIT(8)  ), &
     &(ICTL      ,MODEL_FT_UNIT(9)  ),(RSUB       ,MODEL_FT_UNIT(10) ), &
     &(XHIST     ,MODEL_FT_UNIT(11) ),(THIST      ,MODEL_FT_UNIT(12) ), &
     &(ICECALVE  ,MODEL_FT_UNIT(13) ),                                  &
     &(CACHE1    ,MODEL_FT_UNIT(15) ),(CACHE2     ,MODEL_FT_UNIT(16) ), &
     &                                (ASWAP      ,MODEL_FT_UNIT(18) ), &
     &(OSWAP     ,MODEL_FT_UNIT(19) ),(AINITIAL   ,MODEL_FT_UNIT(20) ), &
     &(ASTART    ,MODEL_FT_UNIT(21) ),(ARESTART   ,MODEL_FT_UNIT(22) ), &
     &(AOPSUM1   ,MODEL_FT_UNIT(23) ),(AOPSUM2    ,MODEL_FT_UNIT(24) ), &
     &(AOPSUM3   ,MODEL_FT_UNIT(25) )
!
      EQUIVALENCE                                                       &
     &(AOPSUM4   ,MODEL_FT_UNIT(26) ),(AOMEAN     ,MODEL_FT_UNIT(27) ), &
     &(ATMANL    ,MODEL_FT_UNIT(28) ),(SSU        ,MODEL_FT_UNIT(29) ), &
     &(OZONE     ,MODEL_FT_UNIT(30) ),(SMCSNOWD   ,MODEL_FT_UNIT(31) ), &
     &(DSOILTMP  ,MODEL_FT_UNIT(32) ),(SOILTYPE   ,MODEL_FT_UNIT(33) ), &
     &(GENLAND   ,MODEL_FT_UNIT(34) ),(SSTIN      ,MODEL_FT_UNIT(35) ), &
     &(SICEIN    ,MODEL_FT_UNIT(36) ),(PERTURB    ,MODEL_FT_UNIT(37) ), &
     &(CURNTIN   ,MODEL_FT_UNIT(38) ),(MASK       ,MODEL_FT_UNIT(39) ), &
     &(OINITIAL  ,MODEL_FT_UNIT(40) ),(OSTART     ,MODEL_FT_UNIT(41) ), &
     &(ORESTART  ,MODEL_FT_UNIT(42) ),(AOPSTMP1   ,MODEL_FT_UNIT(43) ), &
     &(AOPSTMP2  ,MODEL_FT_UNIT(44) ),(AOPSTMP3   ,MODEL_FT_UNIT(45) ), &
     &(AOPSTMP4  ,MODEL_FT_UNIT(46) ),(OCNANL     ,MODEL_FT_UNIT(47) ), &
     &(ATRACER   ,MODEL_FT_UNIT(48) ),(OTRACER    ,MODEL_FT_UNIT(49) ), &
     &(WFIN      ,MODEL_FT_UNIT(50) )
!
      EQUIVALENCE                                                       &
     &(HFLUXIN   ,MODEL_FT_UNIT(51) ),(PMEIN      ,MODEL_FT_UNIT(52) ), &
     &(ICEFIN    ,MODEL_FT_UNIT(53) ),(AIRTMP     ,MODEL_FT_UNIT(54) ), &
     &                                (FLUXCORR   ,MODEL_FT_UNIT(56) ), &
     &(SWSPECTD  ,MODEL_FT_UNIT(57) ),(BAS_IND    ,MODEL_FT_UNIT(58) ), &
     &(SLABHCON  ,MODEL_FT_UNIT(59) ),(PP0        ,MODEL_FT_UNIT(60) ), &
     &(PP1       ,MODEL_FT_UNIT(61) ),(PP2        ,MODEL_FT_UNIT(62) ), &
     &(PP3       ,MODEL_FT_UNIT(63) ),(PP4        ,MODEL_FT_UNIT(64) ), &
     &(PP5       ,MODEL_FT_UNIT(65) ),(PP6        ,MODEL_FT_UNIT(66) ), &
     &(PP7       ,MODEL_FT_UNIT(67) ),(PP8        ,MODEL_FT_UNIT(68) ), &
     &(PP9       ,MODEL_FT_UNIT(69) ),(OBS01      ,MODEL_FT_UNIT(70) ), &
     &(OBS02     ,MODEL_FT_UNIT(71) ),(OBS03      ,MODEL_FT_UNIT(72) ), &
     &(OBS04     ,MODEL_FT_UNIT(73) ),(OBS05      ,MODEL_FT_UNIT(74) ), &
     &(DUSTSOIL  ,MODEL_FT_UNIT(75) ),(BIOMASS    ,MODEL_FT_UNIT(76) ), &
     &(RIVSTOR   ,MODEL_FT_UNIT(77) ),(RIVCHAN    ,MODEL_FT_UNIT(78) ), &
     &(RIVER2A   ,MODEL_FT_UNIT(79) )
!
      EQUIVALENCE                                                       &
                                      (lwspectd   ,model_ft_unit(80) ), &
      (surgeou1  ,model_ft_unit(81) ),(surgeout   ,model_ft_unit(82) ), &
      (ppscreen  ,model_ft_unit(83) ),(ppsmc      ,model_ft_unit(84) ), &
      (wfout     ,model_ft_unit(85) ),(uarsout1   ,model_ft_unit(86) ), &
      (uarsout2  ,model_ft_unit(87) ),(icefout    ,model_ft_unit(88) ), &
      (mosout    ,model_ft_unit(89) ),(vert_lev   ,model_ft_unit(90) ), &
      (sstout    ,model_ft_unit(91) ),(siceout    ,model_ft_unit(92) ), &
      (curntout  ,model_ft_unit(93) ),(flxcrout   ,model_ft_unit(94) ), &
      (dmsconc   ,model_ft_unit(95) ),(orog       ,model_ft_unit(96) ), &
      (transp    ,model_ft_unit(97) ),(olabcin    ,model_ft_unit(98) ), &
      (ocndepth  ,model_ft_unit(99) ),                                  &
      (foamout1  ,model_ft_unit(100)),(foamout2   ,model_ft_unit(101)), &
      (cxbkgerr  ,model_ft_unit(102)),(rfmout     ,model_ft_unit(103)), &
      (idealise  ,model_ft_unit(106)),(tdf_dump   ,model_ft_unit(107)), &
      (iau_inc   ,model_ft_unit(108)),(murkfile   ,model_ft_unit(109)), &
      (sulpemis  ,model_ft_unit(110)),(usrancil   ,model_ft_unit(111)), &
      (usrmulti  ,model_ft_unit(112)),(ousrancl   ,model_ft_unit(113)), &
      (ousrmult  ,model_ft_unit(114)),(so2natem   ,model_ft_unit(115)), &
      (chemoxid  ,model_ft_unit(116)),(aerofcg    ,model_ft_unit(117)), &
      (co2emits  ,model_ft_unit(118)),(tppsozon   ,model_ft_unit(119)), &
      (landfrac  ,model_ft_unit(120)),(wlabcou1   ,model_ft_unit(121)), &
      (wlabcou2  ,model_ft_unit(122)),(wlabcou3   ,model_ft_unit(123)), &
      (wlabcou4  ,model_ft_unit(124)),(alabcin1   ,model_ft_unit(125)), &
      (alabcin2  ,model_ft_unit(126)),                                  &
      (ocffemis  ,model_ft_unit(128)),(horzgrid   ,model_ft_unit(129)), &
      (surfemis  ,model_ft_unit(130)),(aircrems   ,model_ft_unit(131)), &
      (stratems  ,model_ft_unit(132)),(extraems   ,model_ft_unit(133)), &
      (radonems  ,model_ft_unit(134)),(fracinit   ,model_ft_unit(135)), &
      (veginit   ,model_ft_unit(136)),(disturb    ,model_ft_unit(137)), &
      (cached    ,model_ft_unit(138)),(sootemis   ,model_ft_unit(139)), &
      (alabcou1  ,model_ft_unit(140)),(alabcou2   ,model_ft_unit(141)), &
      (alabcou3  ,model_ft_unit(142)),(alabcou4   ,model_ft_unit(143)), &
      (alabcou5  ,model_ft_unit(144)),(alabcou6   ,model_ft_unit(145)), &
      (alabcou7  ,model_ft_unit(146)),(alabcou8   ,model_ft_unit(147)), &
      (cariolo3  ,model_ft_unit(148)),(rpseed     ,model_ft_unit(149)), &
      (ppvar     ,model_ft_unit(150)),(pp10       ,model_ft_unit(151)), &
      (icfile    ,model_ft_unit(152)),(var_grid   ,model_ft_unit(153)), &
      (arclbiog  ,model_ft_unit(154)),(arclbiom   ,model_ft_unit(155)), &
      (arclblck  ,model_ft_unit(156)),(arclsslt   ,model_ft_unit(157)), &
      (arclsulp  ,model_ft_unit(158)),(arcldust   ,model_ft_unit(159)), &
      (arclocff  ,model_ft_unit(160)),(arcldlta   ,model_ft_unit(161)), &
      (topmean   ,model_ft_unit(162)),(topstdev   ,model_ft_unit(163)), &
      (ppmbc     ,model_ft_unit(164)),(ukcaprec   ,model_ft_unit(165)), &
      (ukcaacsw  ,model_ft_unit(166)),(ukcaaclw   ,model_ft_unit(167)), &
      (ukcacrsw  ,model_ft_unit(168)),(ukcacrlw   ,model_ft_unit(169)), &
      (ukcafjxx  ,model_ft_unit(170)),(ukcafjsc   ,model_ft_unit(171)), &
      (ukca2do3  ,model_ft_unit(172)),(ukca2ch4   ,model_ft_unit(173)), &
      (ukca2noy  ,model_ft_unit(174)),(ukca2pho   ,model_ft_unit(175)), &
      (ukcastrd  ,model_ft_unit(176)),(ukcasto3   ,model_ft_unit(177)), &
      (ukcastar  ,model_ft_unit(178)),(ukcafjar   ,model_ft_unit(179))
! Text output file for STASH-related information is assigned to UNIT 200

!
! STPARAM
!
!  Purpose: Meaningful PARAMETER names for STASH processing routines.
!           Both a long name and short name have been declared, to
!           reduce the existence of "magic" numbers in STASH.
!           Format is that first the address of the item is declare in
!           both long and short form. example is;
!             integer st_item_code,s_item  !Item number (declaration)
!             parameter(st_item_code=3,s_item=3)
!
!  Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!  Logical components covered: D70
!
!  Project task: D7
!
!  External documentation:
!    Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                 system (STASH)
!--------------------------------------------------------------

      ! Internal model number address
      INTEGER,PARAMETER:: st_model_code = 28
      INTEGER,PARAMETER:: s_modl        = 28

      ! Section Number address
      INTEGER,PARAMETER:: st_sect_no_code = 2
      INTEGER,PARAMETER:: s_sect          = 2
      INTEGER,PARAMETER:: st_sect_code    = 2

      INTEGER,PARAMETER:: st_item_code=1,s_item=1 ! Item number address

      ! Processing Code address
      INTEGER,PARAMETER:: st_proc_no_code=3,s_proc=3

      ! subsidiary codes for st_proc_no_code now
      INTEGER,PARAMETER:: st_replace_code=1
      INTEGER,PARAMETER:: st_accum_code=2
      INTEGER,PARAMETER:: st_time_mean_code=3
      INTEGER,PARAMETER:: st_time_series_code=4
      INTEGER,PARAMETER:: st_max_code=5
      INTEGER,PARAMETER:: st_min_code=6
      INTEGER,PARAMETER:: st_append_traj_code=7
      INTEGER,PARAMETER:: st_time_series_mean=8
      INTEGER,PARAMETER:: st_variance_code=9

      ! Frequency (Input & output) addres
      INTEGER,PARAMETER:: st_freq_code=4,s_freq=4

      ! Offset for sampling
      INTEGER,PARAMETER:: st_offset_code=30,s_offs=30

      ! start timestep address
      INTEGER,PARAMETER:: st_start_time_code=5,s_times=5

      ! end timestep address
      INTEGER,PARAMETER:: st_end_time_code=6,s_timee=6

      ! period in timesteps address
      INTEGER,PARAMETER:: st_period_code=7,s_period=7

      ! infinite end/period value
      INTEGER,PARAMETER:: st_infinite_time=-1

      INTEGER,PARAMETER:: st_end_of_list=-1 !end-of-list marker in times

      ! grid point stuff
      ! gridpoint info address
      INTEGER,PARAMETER:: st_gridpoint_code=8,s_grid=8

      ! now subsid grid point stuff
      ! no masking done
      INTEGER,PARAMETER:: stash_null_mask_code=1,s_nomask=1

      ! land mask conds
      INTEGER,PARAMETER:: stash_land_mask_code=2,s_lndms=2

      ! sea mask code
      INTEGER,PARAMETER:: stash_sea_mask_code=3,s_seams =3

      ! processing options

      ! size of block for gridpoint code
      INTEGER,PARAMETER:: block_size=10

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: extract_base=block_size*0

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: extract_top=block_size*1

      ! max code for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_base=block_size*1

      ! base codes for vertical mean subroutine
      INTEGER,PARAMETER:: vert_mean_top=block_size*2

      ! max code for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_base=block_size*2

      ! base codes for zonal mean subroutine
      INTEGER,PARAMETER:: zonal_mean_top=block_size*3

      ! max code for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_base=block_size*3

      ! base codes for meridional mean subroutine
      INTEGER,PARAMETER:: merid_mean_top=block_size*4

      ! max code for field mean subroutine
      INTEGER,PARAMETER:: field_mean_base=block_size*4

      ! base codes for field mean subroutine
      INTEGER,PARAMETER:: field_mean_top=block_size*5

      ! max code for global mean subroutine
      INTEGER,PARAMETER:: global_mean_base=block_size*5

      ! base codes for global mean subroutine
      INTEGER,PARAMETER:: global_mean_top=block_size*6

      ! Weighting

      ! weighting info address
      INTEGER,PARAMETER:: st_weight_code=9,s_weight=9

      INTEGER,PARAMETER:: stash_weight_null_code  =0,s_noweight  =0
      INTEGER,PARAMETER:: stash_weight_area_code  =1,s_areaweight=1
      INTEGER,PARAMETER:: stash_weight_volume_code=2,s_volweight =2
      INTEGER,PARAMETER:: stash_weight_mass_code  =3,s_massweight=3

      ! Domain definition

      ! row addresses
      INTEGER,PARAMETER:: st_north_code=12,s_north=12
      INTEGER,PARAMETER:: st_south_code=13,s_south=13
      INTEGER,PARAMETER:: st_west_code =14,s_west =14
      INTEGER,PARAMETER:: st_east_code =15,s_east =15

      ! Levels

      ! input bottom level address
      INTEGER,PARAMETER:: st_input_bottom=10,s_bottom =10

      ! special code
      INTEGER,PARAMETER:: st_special_code=100,s_special=100

      ! input top level address
      INTEGER,PARAMETER:: st_input_top=11,s_top=11

      ! output bottom level address
      INTEGER,PARAMETER:: st_output_bottom=21,s_outbot=21

      ! output top level address
      INTEGER,PARAMETER:: st_output_top=22,s_outtop=22

      INTEGER,PARAMETER:: st_model_level_code=1,s_model=1

      ! code for pressure leve
      INTEGER,PARAMETER:: st_pressure_level_code=2,s_press=2

      ! code for height levels
      INTEGER,PARAMETER:: st_height_level_code=3,s_height=3

      ! input code addres
      INTEGER,PARAMETER:: st_input_code=16,s_input=16

      ! input length of diagnostic address
      INTEGER,PARAMETER:: st_input_length=17,s_length=17

      ! output code address
      INTEGER,PARAMETER:: st_output_code=18,s_output=18

      ! Pointer to D1 addressing information
      ! Pos of item in D1 for relevant submodel
      INTEGER,PARAMETER:: st_position_in_d1=29,st_d1pos=29

      ! Output destination options

      INTEGER,PARAMETER:: st_dump=1
      INTEGER,PARAMETER:: st_secondary=2

      ! output length of diagnostic address
      INTEGER,PARAMETER:: st_output_length=19,s_outlen=19
         integer st_dump_output_length,s_doutlen ! output length on
         parameter(st_dump_output_length=32,s_doutlen=32)  ! dump
         integer st_dump_level_output_length,s_dlevoutlen
         parameter(st_dump_level_output_length=33,s_dlevoutlen=33)
! output length of a single level on dump

         integer st_output_addr,s_outadd ! start locn of diag after stas
         parameter(st_output_addr=20,s_outadd=20)       ! output address
         integer st_dump_output_addr,s_doutadd ! output address on
         parameter(st_dump_output_addr=31,s_doutadd=31)  ! dump

      ! ptr to dump lookup header address
      INTEGER,PARAMETER:: st_lookup_ptr=23

      ! ptr into stash_series where control data address
      INTEGER,PARAMETER:: st_series_ptr=24

      ! subsid stuff for time series
      INTEGER,PARAMETER:: series_grid_type=1
      INTEGER,PARAMETER:: series_grid_code=0
      INTEGER,PARAMETER:: series_long_code=1
      INTEGER,PARAMETER:: series_size=2
      INTEGER,PARAMETER:: series_proc_code=3
      INTEGER,PARAMETER:: series_north=4
      INTEGER,PARAMETER:: series_south=5
      INTEGER,PARAMETER:: series_west=6
      INTEGER,PARAMETER:: series_east=7
      INTEGER,PARAMETER:: series_list_start=8
      INTEGER,PARAMETER:: series_list_end=9
      INTEGER,PARAMETER:: record_size=9

      ! Miscellaneous parameters

      ! system/user tag field in stlist address
      INTEGER,PARAMETER:: st_macrotag=25

      ! Pseudo-level list pointers

      ! pseudo-levels input list address
      INTEGER,PARAMETER:: st_pseudo_in=26

      ! pseudo-levels output list address
      INTEGER,PARAMETER:: st_pseudo_out=27

      ! Internal horizontal gridtype codes common to all diagnostics

      INTEGER,PARAMETER:: st_tp_grid =1 ! T-p grid
      INTEGER,PARAMETER:: st_uv_grid =2 ! u-v grid
      INTEGER,PARAMETER:: st_cu_grid =3 ! C-grid u point
      INTEGER,PARAMETER:: st_cv_grid =4 ! C-grid v point
      INTEGER,PARAMETER:: st_zt_grid =5 ! Zonal T-grid
      INTEGER,PARAMETER:: st_zu_grid =6 ! Zonal u-grid
      INTEGER,PARAMETER:: st_mt_grid =7 ! Meridional T-grid
      INTEGER,PARAMETER:: st_mu_grid =8 ! Meridional u-grid
      INTEGER,PARAMETER:: st_riv_grid= 23    ! river_routing grid
      INTEGER,PARAMETER:: st_scalar  =9 ! Scalar (ie. single value)
      INTEGER,PARAMETER:: st_wam_all= 60    ! Wam Field on Full Grid
      INTEGER,PARAMETER:: st_wam_sea= 62    ! Wam Field on Sea Points

! STPARAM end
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
! COMDECK CSTASH
! Description:
!  Contains variables and arrays involved in STASH request specification
!  and STASH processing in the UM. Includes namelist STASH.
!
!
! Declarations:
! Imported global variables:
!    None, but see note above.

! Global parameters:

! Global scalars:
      INTEGER      NDIAG   ! No. of diagnostics
      INTEGER      NTPROF  ! No. of time profiles
      INTEGER      NSERIES ! No. of stash time series
      INTEGER      NDPROF  ! No. of domain profiles
      INTEGER      NUPROF  ! No. of useage profiles

! Global dynamic arrays:

!   STASH specification table (JSTASH file):
!   NDIAGPM set in VERSION_MOD module
      INTEGER      MODL_B(NDIAGPM)  ! Internal model no.
      INTEGER      ISEC_B(NDIAGPM)  ! Section
      INTEGER      ITEM_B(NDIAGPM)  ! Item
      INTEGER      ITIM_B(NDIAGPM)  ! Time profile number
      INTEGER      IDOM_B(NDIAGPM)  ! Domain profile number
      INTEGER      IUSE_B(NDIAGPM)  ! Useage profile number

!   Time profile information:

      CHARACTER(LEN=8)  TIMPRO(NPROFTP)         ! Name of profile
      INTEGER      ITYP_T(NPROFTP)         ! Type of profile
      INTEGER      INTV_T(NPROFTP)         ! Time Interval
      CHARACTER(LEN=2)  UNT1_T(NPROFTP)         ! Units for time interval
      INTEGER      ISAM_T(NPROFTP)         ! Sampling period
      CHARACTER(LEN=2)  UNT2_T(NPROFTP)         ! Units for sampling period
      INTEGER      IOPT_T(NPROFTP)         ! Output option
      INTEGER      ISTR_T(NPROFTP)         ! Output Start time
      INTEGER      IEND_T(NPROFTP)         ! Output End time
      INTEGER      ISDT_T(6, NPROFTP)      ! Output Start date
      INTEGER      IEDT_T(6, NPROFTP)      ! Output End date
      INTEGER      IFRE_T(NPROFTP)         ! Output frequency
      INTEGER      IOFF_T(NPROFTP)         ! Offset for sampling
      CHARACTER(LEN=2)  UNT3_T(NPROFTP)         ! Units for output times
      INTEGER      ITIM_T(NPROFTP)         ! No. of times in times table
      INTEGER      ISER_T(NTIMEP ,NPROFTP) ! Times table (with units)
      INTEGER      MODL_T(NPROFTP)         ! Indicates internal model
                                           !  for each times table

!   Domain profile information:

      CHARACTER(LEN=8) DOMPRO  (NPROFDP)           ! Name of domain profile
      INTEGER     IOPL_D  (NPROFDP)           ! Levels option
      INTEGER     LEVB_D  (NPROFDP)           ! Bottom level
      INTEGER     LEVT_D  (NPROFDP)           ! Top level
      INTEGER     IOPA_D  (NPROFDP)           ! Area option
      INTEGER     INTH_D  (NPROFDP)           ! North boundary
      INTEGER     ISTH_D  (NPROFDP)           ! South boundary
      INTEGER     IEST_D  (NPROFDP)           ! East boundary
      INTEGER     IWST_D  (NPROFDP)           ! West boundary
      INTEGER     IMSK_D  (NPROFDP)           ! Mask type
      INTEGER     IMN_D   (NPROFDP)           ! Meaning option
      INTEGER     IWT_D   (NPROFDP)           ! Weighting option
      CHARACTER(LEN=1) TS_D    (NPROFDP)           ! Time series profile
      INTEGER     IG_TS
      INTEGER     I1_TS
      INTEGER     I51_TS
      INTEGER     BLIM_TS (NTimSerP)
      INTEGER     TLIM_TS (NTimSerP)
      REAL        BLIMR_TS(NTimSerP)
      REAL        TLIMR_TS(NTimSerP)
      INTEGER     NLIM_TS (NTimSerP)
      INTEGER     SLIM_TS (NTimSerP)
      INTEGER     ELIM_TS (NTimSerP)
      INTEGER     WLIM_TS (NTimSerP)
      INTEGER     ILEV_D  (NPROFDP)           ! Output levels code
      INTEGER     LEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      REAL       RLEVLST_D(NLEVP   ,NPROFDP ) ! Levels list
      INTEGER     PLT_D   (NPROFDP)
      INTEGER     PLLEN_D (NPROFDP)
      INTEGER     PLPOS_D (NPROFDP)
      INTEGER     PSLIST_D(NPSLEVP ,NPSLISTP)
      INTEGER     NPSLISTS
      EQUIVALENCE        (RLEVLST_D,LEVLST_D)

! Useage information:

      CHARACTER(LEN=8) USEPRO(NPROFUP)   ! Name of useage profile
      INTEGER     LOCN_U(NPROFUP)   ! Storage location of profile
      INTEGER     IUNT_U(NPROFUP)   ! Unit no.

! Information from ppxref file:

      INTEGER      MODEL_ST       ! Internal model number
      INTEGER      ISPACE         ! Space code
      INTEGER      ITIMA          ! Time availability code
      INTEGER      IGP            ! Grid of data code
      INTEGER      ILEV           ! Level type code
      INTEGER      IBOT           ! First level code
      INTEGER      ITOP           ! Last level code
      INTEGER      IFLAG          ! Level compression flag
      INTEGER      IOPN(6)        ! Sectional option code
      INTEGER      VMSK           ! Integer equiv of bin vers mask
      INTEGER      IPSEUDO        ! Pseudo dimension type
      INTEGER      IPFIRST        ! First pseudo dim code
      INTEGER      IPLAST         ! Last pseudo dim code
      INTEGER      PTR_PROG       ! Section zero point back
      INTEGER      HALO_TYPE      ! Type of halo the field has

! PP output file units
      INTEGER      PPlen2LkUp(OUTFILE_S:OUTFILE_E)
      CHARACTER(LEN=1)  FTOutUnit (OUTFILE_S:OUTFILE_E)

! COMMON blocks:
      COMMON/STCHA/ TIMPRO,UNT1_T,UNT2_T,UNT3_T,DOMPRO,TS_D,            &
     &  USEPRO,FTOutUnit

      COMMON/STSH/                                                      &
     &  NDIAG   ,MODL_B  ,ISEC_B ,ITEM_B  ,ITIM_B  ,IDOM_B  ,IUSE_B,    &
     &  NTPROF  ,ITYP_T  ,INTV_T ,ISAM_T  ,ITIM_T  ,                    &
     &  IOPT_T  ,ISTR_T  ,IEND_T ,IFRE_T  ,IOFF_T, ISER_T  ,MODL_T  ,   &
     &  NDPROF  ,IOPL_D  ,LEVB_D ,ISDT_T  ,IEDT_T  ,                    &
     &  IOPA_D  ,INTH_D  ,ISTH_D ,IEST_D  ,IWST_D  ,                    &
     &  IMSK_D  ,IMN_D   ,IWT_D  ,                                      &
     &  LEVT_D  ,LEVLST_D,                                              &
     &  PLT_D   ,PLLEN_D ,PLPOS_D,PSLIST_D,NPSLISTS,                    &
     &  BLIM_TS ,TLIM_TS ,BLIMR_TS,TLIMR_TS,IG_TS   ,I1_TS   ,          &
     &  NLIM_TS ,SLIM_TS ,ELIM_TS ,WLIM_TS ,I51_TS  ,NSERIES ,          &
     &  NUPROF  ,LOCN_U  ,IUNT_U ,                                      &
     &  MODEL_ST,ISPACE  ,ITIMA  ,IGP     ,                             &
     &  ILEV    ,IBOT    ,ITOP   ,IFLAG   ,IOPN    ,VMSK    ,           &
     &  IPSEUDO ,IPFIRST ,IPLAST ,PTR_PROG, HALO_TYPE,                  &
     & PPlen2LkUp

! CSTASH end
! MODEL Defines model-dependent quantities used by data addressing and
! STASH
!
! submodel_mod must be used before this one
! VERSION_MOD module is required for nsectp, outfile_s and outfile_e
!
      INTEGER, PARAMETER :: AASSETS    = 9
      INTEGER, PARAMETER :: MEAD_TYPES = 4
      INTEGER, PARAMETER :: A_MAX_TRVARS=150 !Max.no.of tracers allowed
      INTEGER, PARAMETER :: A_MAX_UKCAVARS=150 ! Max.no.of UKCA allowed
      INTEGER, PARAMETER :: MAX_AOBS=100

      REAL :: H_A_EWSPACE
      REAL :: H_A_NSSPACE
      REAL :: H_A_FIRSTLAT
      REAL :: H_A_FIRSTLONG
      REAL :: H_A_POLELAT
      REAL :: H_A_POLELONG

      INTEGER :: H_A_GROUP
      INTEGER :: H_OROG_ROUGH
      INTEGER :: A_ASSMGRPS
      INTEGER :: NUM_PVPR

      LOGICAL :: A_RECON
      LOGICAL :: H_OROG_GRAD
      LOGICAL :: ATMODS
      LOGICAL :: CMODS
      LOGICAL :: LMESO

      LOGICAL :: TRACER_A (0:A_MAX_TRVARS)
      LOGICAL :: TR_UKCA_A (0:A_MAX_UKCAVARS)
      LOGICAL :: AASSET   (AASSETS)
      INTEGER :: AASPF    (AASSETS)
      INTEGER :: AASPL    (AASSETS)
      INTEGER :: AOBINC   (MAX_AOBS)
      INTEGER :: AOBGRP   (MAX_AOBS)
      INTEGER :: RUN_TARGET_END( 6)

      COMMON/MODELA/ H_A_EWSPACE,H_A_NSSPACE,H_A_FIRSTLAT,H_A_FIRSTLONG,&
     &  H_A_POLELAT,H_A_POLELONG,A_ASSMGRPS,NUM_PVPR ,A_RECON,H_A_GROUP,&
     &  H_OROG_GRAD,ATMODS,CMODS,LMESO,TRACER_A,TR_UKCA_A,              &
     &  AASSET,AASPF,AASPL

!Total data length for primary fields for each submodel data partition
      INTEGER      LPRIM(N_SUBMODEL_PARTITION_MAX)
! Global (ie. dump on disk) version of LPRIM
      INTEGER      global_LPRIM(N_SUBMODEL_PARTITION_MAX)
!Total data length for primary fields for each internal model
      INTEGER      LPrimIM(N_INTERNAL_MODEL_MAX)
!Total data length for diagnostic flds for each submodel data partition
! Global (ie. dump on disk) version of LPrimIM
      INTEGER      global_LPrimIM(N_INTERNAL_MODEL_MAX)
      INTEGER      LDUMP(N_SUBMODEL_PARTITION_MAX)
! Global (ie. dump on disk) version of LDUMP
      INTEGER      global_LDUMP(N_SUBMODEL_PARTITION_MAX)
!Total data length for diagnostic flds for each internal model
      INTEGER      LDumpIM(N_INTERNAL_MODEL_MAX)
! Global (ie. dump on disk) version of LDumpIM
      INTEGER      global_LDumpIM(N_INTERNAL_MODEL_MAX)
!Total data length for secondary flds for each submodel data partition
      INTEGER      LSECD(N_SUBMODEL_PARTITION_MAX)
!Total data length for secondary flds for each internal model
      INTEGER      LSecdIM(N_INTERNAL_MODEL_MAX)
!Total workspace length for each submodel data partition
      INTEGER      LWORK(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each submodel data partition
      INTEGER      NHeadSub(N_SUBMODEL_PARTITION_MAX)
!Total number of headers (i.e. levels) for each internal model
      INTEGER      NHEAD(N_INTERNAL_MODEL_MAX)
!Total length of extra space for each submod. data part.
      INTEGER      LEXTRA(N_SUBMODEL_PARTITION_MAX)
!Data length for dual-time level ocean fields
      INTEGER      LPRIM_O2
      INTEGER      ITEM_MAX_REQ
      INTEGER      ITEM_MAX_ALL

      INTEGER      NRECS_S
      INTEGER      NTIMES_S
      INTEGER      NSERBLK_S
      INTEGER      NSERREC_S
      INTEGER      NLEVL_S
      INTEGER      NMAXLEV_S
      INTEGER      NPSLISTS_S
      INTEGER      NMAXPSL_S
      INTEGER      NHEAD_FILE(OUTFILE_S:OUTFILE_E)
      LOGICAL      LSTUSER

      COMMON/STRET/                                                     &
     &  LPRIM,LDUMP,LSECD,LWORK,NHEAD,LEXTRA,LPRIM_O2,LPrimIM,LDumpIM,  &
     &  LSecdIM,NHeadSub,ITEM_MAX_REQ,ITEM_MAX_ALL,NSERBLK_S,NSERREC_S, &
     &  NLEVL_S,NMAXLEV_S,NPSLISTS_S,NMAXPSL_S,LSTUSER,NRECS_S,NTIMES_S,&
     &  NHEAD_FILE,                                                     &
     &  global_LPRIM,global_LPrimIM,global_LDUMP,global_LDumpIM
      CHARACTER(LEN=1)  H_ATMOS
      CHARACTER(LEN=1)  H_FLOOR
      CHARACTER(LEN=1)  H_STRAT
      CHARACTER(LEN=1)  H_GLOBAL(N_INTERNAL_MODEL_MAX         )
      INTEGER      H_VERS  (N_INTERNAL_MODEL_MAX,0:NSECTP)

      COMMON/CHOICE/ H_ATMOS,H_GLOBAL,H_FLOOR,H_STRAT

      COMMON/HVERS/ H_VERS

! These are set in SETMODL:
      INTEGER MEAN_NUMBER(N_INTERNAL_MODEL_MAX)
      COMMON/MODLMEAN/ MEAN_NUMBER


! Variables read in by namelist and used in SETMODL
      INTEGER      OCAAA 
      REAL         EWSPACEA,NSSPACEA
      REAL         FRSTLATA,FRSTLONA

      LOGICAL      ZonAvOzone
      LOGICAL      ZonAvTppsOzone
      REAL         LATS
      REAL         LONS
      INTEGER      LWBND
      INTEGER      OCALB
      REAL         POLELATA
      REAL         POLELONA
      INTEGER      SWBND
      INTEGER      TCA(A_MAX_TRVARS)
      INTEGER      TCA_LBC(A_MAX_TRVARS)  ! =1 if tracer in lbc file 
      INTEGER      TC_UKCA(A_MAX_UKCAVARS)
      INTEGER      TC_LBC_UKCA(A_MAX_UKCAVARS) ! =1 if tr in lbc file 
      INTEGER      StLevGWdrag
      INTEGER      BotVDiffLev
      INTEGER      TopVDiffLev


      COMMON/STSHCOMM/                                                  &
     &  RUN_TARGET_END,                                                 &
     &  OCAAA,EWSPACEA,POLELATA,FRSTLATA,LATS,                          &
     &  NSSPACEA,POLELONA,FRSTLONA,LONS,                                &
     &  SWBND,LWBND,                                                    &
     &  ZonAvOzone ,ZonAvTppsOzone, AOBINC,  AOBGRP,                    &
     &  StLevGWdrag, BotVDiffLev,TopVDiffLev,                           &
     &  OCALB,TCA,TCA_LBC,TC_UKCA,TC_LBC_UKCA


      CHARACTER(LEN=1) :: LFLOOR
      CHARACTER(LEN=1) :: OROGR
      CHARACTER(LEN=1) :: SWMCR
      CHARACTER(LEN=1) :: MESO

      COMMON/STSHCHAR/                                                  &
     &     LFLOOR,                                                      &
     &  OROGR,   SWMCR, MESO

      NAMELIST/STSHCOMP/                                                &
        RUN_TARGET_END,                                                 &
        OCAAA       ,EWSPACEA    ,POLELATA ,FRSTLATA  ,LATS   ,         &
                     NSSPACEA    ,POLELONA ,FRSTLONA  ,LONS   ,         &
        SWBND       ,LWBND                            ,OROGR  ,         &
        ZonAvOzone  ,SWMCR       ,MESO     ,                            &
        OCALB       ,LFLOOR      ,AOBINC   ,TCA,                        &
        TCA_LBC     ,TC_UKCA     ,TC_LBC_UKCA   ,AOBGRP          
  
! MODEL end
! For accessing D1 addressing array
! TYPD1 Common block containing the ALT_N_SUBMODEL_PARTITION variables
! CALTSUBM
! TYPD1 needs access to N_SUBMODEL_PARTITION/_MAX in CSUBMODL. However,
! they are not always called in the same decks and in the right order.
! Therefore, copy the values to another file and include it from TYPD1

      INTEGER ALT_N_SUBMODEL_PARTITION

      INTEGER, PARAMETER :: ALT_N_SUBMODEL_PARTITION_MAX=1

      COMMON/CALTSUBM/ALT_N_SUBMODEL_PARTITION
! CALTSUBM end
! This file needs TYPSIZE included first

      REAL    ::  D1(LEN_TOT)       ! IN/OUT: Main data array
      LOGICAL :: LD1(LEN_TOT)       ! IN/OUT: Main data array (logical)
      INTEGER :: ID1(LEN_TOT)       ! I/OUT: Main data array (integer)

! D1_ADDR start
      ! Information for accessing D1 addressing array
      ! Number of items of info needed for each object and maximum
      ! number of objects in D1 -

      ! Number of items of information in D1 addressing array
      INTEGER,PARAMETER:: D1_LIST_LEN=17

! Names of items in D1 addressing array. Update D1_LIST_LEN above if
! items added

      ! Prognostic, Diagnostic, Secondary or other
      INTEGER,PARAMETER:: d1_object_type    = 1 ! Internal model id
      INTEGER,PARAMETER:: d1_imodl          = 2  ! Internal model id
      INTEGER,PARAMETER:: d1_section        = 3  ! Section
      INTEGER,PARAMETER:: d1_item           = 4  ! Item
      INTEGER,PARAMETER:: d1_address        = 5  ! Address in D1
      INTEGER,PARAMETER:: d1_length         = 6  ! Record length
      INTEGER,PARAMETER:: d1_grid_type      = 7  ! Grid type
      INTEGER,PARAMETER:: d1_no_levels      = 8  ! Number of levels

      ! Stash list number for diags. -1 for progs
      INTEGER,PARAMETER:: d1_stlist_no      = 9

      ! Pointer to dump header lookup table
      INTEGER,PARAMETER:: d1_lookup_ptr     = 10

      INTEGER,PARAMETER:: d1_north_code     = 11 ! Northern row
      INTEGER,PARAMETER:: d1_south_code     = 12 ! Southern row
      INTEGER,PARAMETER:: d1_east_code      = 13 ! Eastern row
      INTEGER,PARAMETER:: d1_west_code      = 14 ! Western row
      INTEGER,PARAMETER:: d1_gridpoint_code = 15 ! gridpoint info
      INTEGER,PARAMETER:: d1_proc_no_code   = 16 ! Processing Code
      INTEGER,PARAMETER:: d1_halo_type      = 17 ! Halo width type

      ! Types of items for d1_type

      INTEGER,PARAMETER:: prognostic = 0
      INTEGER,PARAMETER:: diagnostic = 1
      INTEGER,PARAMETER:: secondary  = 2
      INTEGER,PARAMETER:: other      = 3

! D1_ADDR end
      ! D1 addressing array and number of objects in each submodel
      INTEGER :: D1_ADDR(D1_LIST_LEN,N_OBJ_D1_MAX,                      &
     &  ALT_N_SUBMODEL_PARTITION)

      INTEGER :: NO_OBJ_D1(ALT_N_SUBMODEL_PARTITION_MAX)

      COMMON/common_D1_ADDRESS/ NO_OBJ_D1
! TYPD1 end

      INTEGER                                                           &
     &  II,                                                             &
                 ! Addresses preliminary array
     &  SM,                                                             &
                 ! Addresses final array=1 for 1st submod =2 for 2nd
                 ! submodel etc
     &  TYPE,                                                           &
                 ! Code for prognostic, diagnostic, secondary or other
     &  IOBJ,                                                           &
                 ! Addresses final array
     &  ISEC,                                                           &
                 ! Section number
     &  ITM,                                                            &
                 ! Item number
     &  LEVS,                                                           &
                 ! No of levels
     &  INF,                                                            &
              ! Diagnostic STASHlist number or prognosic item number
     &  Im_ident,                                                       &
     &  Sm_ident,                                                       &
     &  LOOKUP_PTR,                                                     &
                    ! Pointer to lookup table
     &  EXT_ADDR,                                                       &
                  ! Temporary pointer
     &  ICODE                   ! OUT: Error return code
!
      CHARACTER(LEN=256)                                                     &
     &    CMESSAGE               ! OUT: Error return message

      INTEGER EXPPXI

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
      EXTERNAL EXPPXI

! Initialise array
      IF (lhook) CALL dr_hook('FILL_D1_ARRAY',zhook_in,zhook_handle)
! Remove superfluous DO Sm_ident=1,N_SUBMODEL_PARTITION
      Sm_ident=1   ! N_SUBMODEL_PARTITION set in UM_Submodel_Init
        DO II=1,N_OBJ_D1_MAX
          DO INF=1,D1_LIST_LEN
            D1_ADDR(INF,II,Sm_ident)=-1
            NO_OBJ_D1(Sm_ident)=0
          END DO
        END DO

      IF (PrintStatus >= PrStatus_Oper) THEN
! Set up addressing of D1
      WRITE(6,'(/,A)') '********************************************'// &
                       '***********************************'
      WRITE(6,'(A)')'Addressing of D1 array'
      WRITE(6,'(A)')'Key to Type:'
      WRITE(6,'(A)')'Type=0: Prognostic'
      WRITE(6,'(A)')'Type=1: Diagnostics in dump'
      WRITE(6,'(A)')'Type=2: Secondary diagnostics'
      WRITE(6,'(A)')'Type=3: Others (eg P_EXNER)'
      END IF  ! PrintStatus test
      SM=0
! Remove superfluous DO Sm_ident=1,N_SUBMODEL_PARTITION_MAX
      Sm_ident=1   ! N_SUBMODEL_PARTITION_MAX is PARAMETER in csubmmax.h
        IOBJ=0
        SM=SUBMODEL_FOR_SM(Sm_ident)
        IF (SM /= 0) THEN
         IF (NO_OBJ_D1(SM) == 0) THEN
          NO_OBJ_D1(SM)=N_OBJ_D1(Sm_ident)
          WRITE(6,'(A,I4)')'Submodel id ',Sm_ident
          WRITE(6,'(A,I4)')'Submodel Number ',SM
          WRITE(6,'(A,I6)')'No of objects in this submodel: ',NO_OBJ_D1(SM)
! Address if submodel not empty and not already addressed
          DO II=1,NO_OBJ_D1(SM)
!           Preliminary array held in D1_PADDR - full array in D1_ADDR
!           Index II in D1_PADDR goes into index IOBJ of D1_ADDR
!           First add prognostics followed by diagnostics...
            Im_ident=D1_PADDR(d1_im,II,Sm_ident)
            INF=D1_PADDR(d1_extra_info,II,Sm_ident)
            ISEC=D1_PADDR(d1_sect,II,Sm_ident)
            TYPE=D1_PADDR(d1_type,II,Sm_ident)
            IF (TYPE == prog) THEN
              IOBJ=IOBJ+1
              D1_ADDR(d1_stlist_no,IOBJ,SM)=INF
              D1_ADDR(d1_section,IOBJ,SM)=ISEC
              D1_ADDR(d1_no_levels,IOBJ,SM)=D1_PADDR(d1_levs,II,Sm_ident)
              D1_ADDR(d1_object_type,IOBJ,SM)=prognostic
              D1_ADDR(d1_imodl,IOBJ,SM)  = Im_ident
              D1_ADDR(d1_address,IOBJ,SM)= IN_S(1,Im_ident,ISEC,INF)
            ELSE IF (TYPE == diag) THEN
              IOBJ=IOBJ+1
              D1_ADDR(d1_stlist_no,IOBJ,SM)=INF
              D1_ADDR(d1_object_type,IOBJ,SM)=diagnostic
              D1_ADDR(d1_imodl,IOBJ,SM)  = Im_ident
              D1_ADDR(d1_address,IOBJ,SM)= STLIST(st_output_addr,INF)
            END IF
          END DO

!         Extra data between primary and secondary diagnostics
          DO II=1,NO_OBJ_D1(SM)
            ISEC=D1_PADDR(d1_sect,II,Sm_ident)
            TYPE=D1_PADDR(d1_type,II,Sm_ident)
            IF (TYPE == extra_d1) THEN
              Im_ident=D1_PADDR(d1_im,II,Sm_ident)
              INF=D1_PADDR(d1_extra_info,II,Sm_ident)
              IOBJ=IOBJ+1
              D1_ADDR(d1_stlist_no,IOBJ,SM)=INF
              D1_ADDR(d1_section,IOBJ,SM)=ISEC
              D1_ADDR(d1_no_levels,IOBJ,SM)=D1_PADDR(d1_levs,II,Sm_ident)
              D1_ADDR(d1_object_type,IOBJ,SM)=other
              D1_ADDR(d1_imodl,IOBJ,SM)  = Im_ident
!               NOT OCEAN: Address was calculated in ADDRES
              D1_ADDR(d1_address,IOBJ,SM)=IN_S(1,Im_ident,ISEC,INF)
            END IF
          END DO
!         Finally add secondary diagnostics
          DO II=1,NO_OBJ_D1(SM)
!           Preliminary array held in D1_PADDR - full array in D1_ADDR
            Im_ident=D1_PADDR(d1_im,II,Sm_ident)
            INF=D1_PADDR(d1_extra_info,II,Sm_ident)
            TYPE=D1_PADDR(d1_type,II,Sm_ident)
            IF (TYPE == seco) THEN
              IOBJ=IOBJ+1
              D1_ADDR(d1_stlist_no,IOBJ,SM)=INF
              D1_ADDR(d1_object_type,IOBJ,SM)=secondary
              D1_ADDR(d1_imodl,IOBJ,SM)  = Im_ident
              D1_ADDR(d1_address,IOBJ,SM)= STLIST(st_output_addr,INF)
            END IF
          END DO

          LOOKUP_PTR=0
          DO II=1,NO_OBJ_D1(SM)
            TYPE= D1_ADDR(d1_object_type,II,SM)
            ISEC= D1_ADDR(d1_section,II,SM)
            INF = D1_ADDR(d1_stlist_no,II,SM)
            Im_ident = D1_ADDR(d1_imodl,II,SM)
            IF ((TYPE == prognostic).OR.(TYPE == other))THEN
! Prognostics don't have STASHlist numbers
              D1_ADDR(d1_stlist_no,II,SM)= -1
              D1_ADDR(d1_item,II,SM)   = INF
              D1_ADDR(d1_length,II,SM) = IN_S(2,Im_ident,ISEC,INF)
              ISEC = D1_ADDR(d1_section,II,SM)
              ITM  = INF
!-------------------------------------------------------------------
! Prognostic items:
! Additional items can be added to the array here. Its code (eg
! d1_item, d1_levels) should be added to the TYPD1 comdeck and
! set as a parameter. The D1_LIST_LEN parameter should be changed
! as required
!-------------------------------------------------------------------
            ELSE
              D1_ADDR(d1_section,II,SM)= STLIST(st_sect_code,INF)
              D1_ADDR(d1_item,II,SM)   = STLIST(st_item_code,INF)
              D1_ADDR(d1_length,II,SM) = STLIST(st_output_length,INF)
              ISEC=D1_ADDR(d1_section,II,SM)
              ITM=D1_ADDR(d1_item,II,SM)
! STASH list pointer to D1 address information
              STLIST(st_position_in_d1,INF) = II
!-------------------------------------------------------------------
! Diagnostic items
! Add items as per prognostics
!-------------------------------------------------------------------
              D1_ADDR(d1_north_code,II,SM)    =STLIST(st_north_code,INF)
              D1_ADDR(d1_south_code,II,SM)    =STLIST(st_south_code,INF)
              D1_ADDR(d1_east_code,II,SM)     =STLIST(st_east_code,INF)
              D1_ADDR(d1_west_code,II,SM)     =STLIST(st_west_code,INF)
              D1_ADDR(d1_gridpoint_code,II,SM)=STLIST(s_grid,INF)
              D1_ADDR(d1_proc_no_code,II,SM)  =STLIST(s_proc,INF)
! 1. Number of levels
              IF (STLIST(st_output_bottom,INF) == 100) THEN
! Special levels
                LEVS=1
              ELSE IF(STLIST(st_series_ptr,INF) /= 0) THEN
! Time series domain
                LEVS=1
              ELSE IF(STLIST(st_gridpoint_code,INF) >= 10               &
     &            .AND.STLIST(st_gridpoint_code,INF) <  20) THEN
! Vertical ave.
                LEVS=1
              ELSE  IF(STLIST(st_output_bottom,INF) <  0) THEN
! Levels list
                LEVS=LEVLST_S(1,-STLIST(st_output_bottom,INF))
              ELSE
! Range of model levels
                LEVS=STLIST(st_output_top   ,INF)                       &
     &            -STLIST(st_output_bottom,INF)+1
              END IF

              IF (STLIST(st_pseudo_out,INF) >  0) THEN
! Pseudo levels
                LEVS=LEVS*LENPLST(STLIST(st_pseudo_out,INF))
              END IF
              D1_ADDR(d1_no_levels,II,SM) = LEVS
            END IF
!-------------------------------------------------------------------
! Items whose settings are common to progs and diags (eg from PPXREF)
! Add items as per prognostics
! ISEC and ITM set above
!-------------------------------------------------------------------
            D1_ADDR(d1_grid_type,II,SM) =                               &
! DEPENDS ON: exppxi
     &        EXPPXI(Im_ident,ISEC,ITM,ppx_grid_type,                   &
                    ICODE, CMESSAGE)
            D1_ADDR(d1_halo_type,II,SM) =                               &
! DEPENDS ON: exppxi
     &        EXPPXI(Im_ident,ISEC,ITM,ppx_halo_type,                   &
                    ICODE, CMESSAGE)
            LOOKUP_PTR=LOOKUP_PTR+D1_ADDR(d1_no_levels,II,SM)
            D1_ADDR(d1_lookup_ptr,II,SM)=LOOKUP_PTR
          END DO
          IF (PrintStatus >= PrStatus_Normal) THEN
          WRITE(6,'(A)')                                                &
     &'      Type Modl Sect Item   Address   Length Levels Gridtype',   &
     &' Halotype'
          DO II=1,NO_OBJ_D1(SM)
            WRITE(6,'(5I5,I11,I9,I6,I7,I8)')                            &
     &     II,D1_ADDR(d1_object_type,II,SM),D1_ADDR(d1_imodl,IOBJ,SM),  &
     &        D1_ADDR(d1_section,II,SM),D1_ADDR(d1_item,II,SM),         &
     &        D1_ADDR(d1_address,II,SM),D1_ADDR(d1_length,II,SM),       &
     &        D1_ADDR(d1_no_levels,II,SM),D1_ADDR(d1_grid_type,II,SM),  &
     &        D1_ADDR(d1_halo_type,II,SM)

          END DO
          END IF  ! PrintStatus test
         END IF ! IF (NO_OBJ_D1(SM) == 0) THEN
        END IF


      WRITE(6,'(A,/)') '********************************************'// &
                       '***********************************'
      IF (lhook) CALL dr_hook('FILL_D1_ARRAY',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FILL_D1_ARRAY
