! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
      SUBROUTINE READCNTL(UNIT,ICODE,CMESSAGE)

USE readwritd_mod, ONLY: readwritd

USE um_input_control_mod, ONLY:                                              &
  check_nlstcatm,                                                            &
  l_so2                ,l_dms              ,                                 &
  l_so4_aitken         ,l_so4_accu         ,l_so4_diss      ,l_nh3          ,&
  l_soot_new           ,l_soot_agd         ,l_soot_cld      ,l_bmass_new    ,&
  l_bmass_agd          ,l_bmass_cld        ,l_ocff_new      ,l_ocff_agd     ,&
  l_ocff_cld           ,l_nitr_acc         ,l_nitr_diss     ,l_dust_div1_lbc,&
  l_dust_div2_lbc      ,l_dust_div3_lbc    ,l_dust_div4_lbc ,l_dust_div5_lbc,&
  l_dust_div6_lbc      ,l_so2_lbc          ,l_dms_lbc       ,                &
  l_so4_aitken_lbc     ,l_so4_accu_lbc     ,l_so4_diss_lbc  ,l_nh3_lbc      ,&
  l_soot_new_lbc       ,l_soot_agd_lbc     ,l_soot_cld_lbc  ,l_bmass_new_lbc,&
  l_bmass_agd_lbc      ,l_bmass_cld_lbc    ,l_ocff_new_lbc  ,l_ocff_agd_lbc ,&
  l_ocff_cld_lbc       ,l_nitr_acc_lbc     ,l_nitr_diss_lbc ,                &
  l_dust_div6_lbc_out  ,l_so2_lbc_out      ,l_dms_lbc_out   ,                &
  l_so4_accu_lbc_out   ,l_so4_diss_lbc_out ,l_nh3_lbc_out   ,                &
  l_soot_agd_lbc_out   ,l_soot_cld_lbc_out ,                                 &
  l_bmass_cld_lbc_out  ,l_ocff_new_lbc_out ,l_use_seasalt_pm,                &
  l_nitr_acc_lbc_out   ,l_nitr_diss_lbc_out,nlstcatm                         &
                       ,lcal360         ,                                    &
  npmsl_height         ,                                                     &
  l_pmsl_sor           ,                                                     &
  l_dust_div1_lbc_out  ,                                                     &
  l_dust_div5_lbc_out  ,                                                     &
  l_so4_aitken_lbc_out ,                                                     &
  l_soot_new_lbc_out   ,l_ocff_agd_lbc_out ,                                 &
  l_bmass_agd_lbc_out  ,l_dust_div2_lbc_out,                                 &
  l_ocff_cld_lbc_out   ,l_dust_div3_lbc_out,                                 &
  l_bmass_new_lbc_out  ,l_dust_div4_lbc_out



USE check_iostat_mod
USE ereport_mod, ONLY : ereport

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

! JULES
USE switches, ONLY: l_360

USE filenamelength_mod, ONLY :                                    & 
          filenamelength

! MPP
USE mpp_conf_mod, ONLY: nlst_mpp

USE domain_params
USE UM_ParVars
USE IO
USE model_file
USE lbc_mod
USE ancilcta_namelist_mod
USE umsections_mod
USE model_id_mod, ONLY: itab, configid
USE nlstgen_mod
USE Submodel_Mod

USE nlstcall_mod

USE chsunits_mod

USE nlstcall_pop_mod, ONLY : pop_nlstcall

IMPLICIT NONE
!
! Description:
!  Reads overall, generic and model-specific control variables.
!
! Method:
!  Reads namelist containing all overall, generic and model-specific
!  control variables not in History file.  Control variables are then
!  available from COMMON block storage. Declarations, namelist and
!  common are all in included comdecks.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.3 programming standards.
!
! Declarations:
!
! Global variables (*CALLed COMDECKs etc...):
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

! Subroutine arguments
!   Scalar arguments with intent(in):
!   Array  arguments with intent(in):
!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):
      INTEGER UNIT  ! namelist input unit no.
      INTEGER ICODE ! return code
      INTEGER ERRORSTATUS
      CHARACTER (LEN=*), PARAMETER  :: ROUTINENAME='READCNTL'
      CHARACTER(LEN=80) CMESSAGE ! error message
!   Array  arguments with intent(out):

! Local parameters:

! Local scalars:
      CHARACTER(LEN=filenamelength) :: filename ! physical filename

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

! Local dynamic arrays:

!*L --------------------- Comdeck: CENVIR   ----------------------------
!
!   Purpose: COMDECK defining Character enviroment variables used
!            by portable IO to open and close files
!
!
! Type declarations
!
      CHARACTER(LEN=8) FT_ENVIRON(199)  ! Array holding enviroment variables
!                                  for filenames
      INTEGER     LEN_FT_ENVIR(199) ! character length of each variable
!
!
!Common Blocks for character and integer arrays
!
      COMMON/CENVIR/FT_ENVIRON
      COMMON/CLENVIR/LEN_FT_ENVIR
!

      INTEGER   :: NFTIN     ! FTN number for read
      INTEGER   :: ERROR     ! Error code returned by OPEN
      INTEGER   :: fixhdr(len_fixhd)
      INTEGER   :: len_io
      REAL      :: a
      INTEGER, ALLOCATABLE  :: lookup(:,:) 

!- End of header

      IF (lhook) CALL dr_hook('READCNTL',zhook_in,zhook_handle)

!  Read UM section choices
      READ (UNIT=UNIT, NML=umsections, IOSTAT=ErrorStatus)
      CALL check_iostat(ErrorStatus, "namelist UMSECTIONS")
      
!  Assign UM sections in ROSE format to atmos_sr and indep_sr
      CALL assign_umsections
      
!  Read overall control data into COMMON/CNTLCALL/
      pp_len2_look(:) = 0
      READ (UNIT=UNIT, NML=NLSTCALL, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist NLSTCALL")

      CALL pop_nlstcall(type_letter_1,type_letter_3)

!  Set control information for WRITD1 diagnostics
      CALL READWRITD()

!  Read configuration id, is itab set explictly as an input?
      READ (UNIT=UNIT, NML=configid, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist configid")

!  Read generic control data from NLSTCGEN/
      READ (UNIT=UNIT, NML=NLSTCGEN, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist NLSTCGEN")

! Read MPP configuration data 
      READ (UNIT=UNIT, NML=NLST_MPP, IOSTAT=ErrorStatus)
      CALL check_iostat(errorstatus, "namelist NLST_MPP")
      
      REWIND(UNIT)

!  Read atmospheric control data into COMMON/CNTLCATM/
!  Continuation runs only read in NLSTCATM once, from UNIT05
      IF (UNIT == 5) THEN

!     Initialise Section Zero Tracer LBC logicals
        L_SO2                = .FALSE.
        L_DMS                = .FALSE.
        L_SO4_AITKEN         = .FALSE.
        L_SO4_ACCU           = .FALSE.
        L_SO4_DISS           = .FALSE.
        L_NH3                = .FALSE.
        L_SOOT_NEW           = .FALSE.
        L_SOOT_AGD           = .FALSE.
        L_SOOT_CLD           = .FALSE.
        L_BMASS_NEW          = .FALSE.
        L_BMASS_AGD          = .FALSE.
        L_BMASS_CLD          = .FALSE.
        L_OCFF_NEW           = .FALSE.
        L_OCFF_AGD           = .FALSE.
        L_OCFF_CLD           = .FALSE.
        L_NITR_ACC           = .FALSE.
        L_NITR_DISS          = .FALSE.
        L_DUST_DIV1_LBC      = .FALSE.
        L_DUST_DIV2_LBC      = .FALSE.
        L_DUST_DIV3_LBC      = .FALSE.
        L_DUST_DIV4_LBC      = .FALSE.
        L_DUST_DIV5_LBC      = .FALSE.
        L_DUST_DIV6_LBC      = .FALSE.
        L_SO2_LBC            = .FALSE.
        L_DMS_LBC            = .FALSE.     
        L_SO4_AITKEN_LBC     = .FALSE.
        L_SO4_ACCU_LBC       = .FALSE.
        L_SO4_DISS_LBC       = .FALSE.
        L_NH3_LBC            = .FALSE.
        L_SOOT_NEW_LBC       = .FALSE.
        L_SOOT_AGD_LBC       = .FALSE.
        L_SOOT_CLD_LBC       = .FALSE.
        L_BMASS_NEW_LBC      = .FALSE.
        L_BMASS_AGD_LBC      = .FALSE.
        L_BMASS_CLD_LBC      = .FALSE.
        L_OCFF_NEW_LBC       = .FALSE.
        L_OCFF_AGD_LBC       = .FALSE.
        L_OCFF_CLD_LBC       = .FALSE.
        L_NITR_ACC_LBC       = .FALSE.
        L_NITR_DISS_LBC      = .FALSE.
        L_DUST_DIV1_LBC_OUT  = .FALSE.
        L_DUST_DIV2_LBC_OUT  = .FALSE.
        L_DUST_DIV3_LBC_OUT  = .FALSE.
        L_DUST_DIV4_LBC_OUT  = .FALSE.
        L_DUST_DIV5_LBC_OUT  = .FALSE.
        L_DUST_DIV6_LBC_OUT  = .FALSE.
        L_SO2_LBC_OUT        = .FALSE.
        L_DMS_LBC_OUT        = .FALSE.
        L_SO4_AITKEN_LBC_OUT = .FALSE.
        L_SO4_ACCU_LBC_OUT   = .FALSE.
        L_SO4_DISS_LBC_OUT   = .FALSE.
        L_NH3_LBC_OUT        = .FALSE.
        L_SOOT_NEW_LBC_OUT   = .FALSE.
        L_SOOT_AGD_LBC_OUT   = .FALSE.
        L_SOOT_CLD_LBC_OUT   = .FALSE.
        L_BMASS_NEW_LBC_OUT  = .FALSE.
        L_BMASS_AGD_LBC_OUT  = .FALSE.
        L_BMASS_CLD_LBC_OUT  = .FALSE.
        L_OCFF_NEW_LBC_OUT   = .FALSE.
        L_OCFF_AGD_LBC_OUT   = .FALSE.
        L_OCFF_CLD_LBC_OUT   = .FALSE.
        L_NITR_ACC_LBC_OUT   = .FALSE.
        L_NITR_DISS_LBC_OUT  = .FALSE.
!     **********    Aerosol scheme defaults   *********
        L_USE_SEASALT_PM  = .FALSE.

        READ (UNIT=UNIT, NML=NLSTCATM, IOSTAT=ErrorStatus)
        CALL check_iostat(errorstatus, "namelist NLSTCATM")
        CALL check_nlstcatm()

        L_SSTANOM=.FALSE.
        LAMIPII=.FALSE.

        READ (UNIT=UNIT, NML=ANCILCTA, IOSTAT=ErrorStatus)
        CALL check_iostat(errorstatus, "namelist ANCILCTA")
        REWIND(UNIT)
      
! overwrite the JULES module values
        l_360         = lcal360

        IF(Num_ALBCs /= 0 )THEN

          NFTIN = 125

          CALL MODEL_FILE_OPEN(NFTIN,FT_ENVIRON(NFTIN),                     &
            LEN_FT_ENVIR(NFTIN),0,0)

          CALL buffin(nftin,fixhdr,len_fixhd,len_io,a)

! Check for I/O errors
          IF (a /= -1.0 .OR. len_io /= len_fixhd) THEN
! DEPENDS ON: ioerror
            CALL ioerror('buffer in of fixed length header',a,len_io  &
               ,len_fixhd)
            cmessage='READCNTL: I/O error'
            icode=1
            IF (lhook) CALL dr_hook('READCNTL',zhook_out,zhook_handle)
            RETURN
          END IF

! Broadcast the fixed length header from proc 0 to all PEs
          CALL gc_ibcast(77,len_fixhd,0,nproc_max,icode,fixhdr)

          IF (icode /= 0) THEN
            Cmessage = 'Error in broadcast'      
            CALL Ereport(RoutineName, icode, cmessage)
          END IF

! Move to start of Look Up Table
          CALL SETPOS(NFTIN,fixhdr(150)-1,ICODE)

          ALLOCATE (lookup(fixhdr(151),fixhdr(152)))

! Read in fields from LOOKUP table
          CALL BUFFIN(NFTIN,lookup(:,:),fixhdr(151)*fixhdr(152),LEN_IO,A)

          IF (a /= -1.0 .OR. len_io /= fixhdr(151)*fixhdr(152)) THEN
! DEPENDS ON: ioerror
            CALL ioerror('buffer in of lookup table',a,len_io,   &
                   fixhdr(151)*fixhdr(152))
            cmessage='READCNTL: I/O error'
            icode=3
            IF (lhook) CALL dr_hook('READCNTL',zhook_out,zhook_handle)
            RETURN
          END IF

          rimwidtha(1)=lookup(41,1)/10000

! Broadcast the rimwidth from proc 0 to all PEs
          CALL gc_ibcast(77,nrima_max,0,nproc_max,icode,rimwidtha)
          IF (icode /= 0) THEN
            Cmessage = 'Error in broadcast'      
            CALL Ereport(RoutineName, icode, cmessage)
          END IF

          CALL MODEL_FILE_CLOSE(NFTIN,FT_ENVIRON(NFTIN),                    &
                      LEN_FT_ENVIR(NFTIN),0,0,ICODE)

          DEALLOCATE (lookup)
        ELSE
          rimwidtha(1)=0
        END IF

      END IF ! IF UNIT == 5

      IF (lhook) CALL dr_hook('READCNTL',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE READCNTL
