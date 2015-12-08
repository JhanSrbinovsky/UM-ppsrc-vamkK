! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Passes out atmosphere LBCs to processors at boundaries
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: MPP

! Subroutine Interface
SUBROUTINE scatter_atmos_lbcs(                                    &
  full_lbc,decomp_lbc,                                            &
  full_lbc_size,decomp_lbc_size,                                  &
  full_lbc_levels,decomp_lbc_levels,                              &
  fld_type,halo_type,rim_type,                                    &
  pe_for_level,                                                   &
  icode,cmessage)

USE dynamics_grid_mod, ONLY: l_vatpoles

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE gcom_mod
USE UM_ParVars
USE lbc_mod
USE mpl, ONLY: mpl_real
IMPLICIT NONE

! Description:
! Scatters atmosphere LBCs to relevant processors at the grid boundaries

! Method:
! I - If master task (0)
! I - Master Task (0) Allocates Buffer
! I  L--Loop over all processors : iproc
! I L L--Loop over sides (North,East,South,West) : iside
! I L L I--IF processor iproc has an extremity at iside THEN
! I L L I  Copy the data to the buffer (data_lbc)
! I L L I--ENDIF
! I L L -- End Loop iside
! I L -- End Loop iproc
! I - ENDIF - master task
! Scatter data ( MPL_Scatterv collective operation )

! Subroutine Arguments:

INTEGER, INTENT(IN) ::  full_lbc_size 
                               ! IN single level size of the FULL_LBC array
INTEGER, INTENT(IN) ::  decomp_lbc_size 
                               ! IN single level size of the DECOMP_LBC
                               !    array
INTEGER, INTENT(IN) :: full_lbc_levels      
                               ! IN number of levels of FULL_LBC on this
                               !    processor
INTEGER, INTENT(IN) ::  decomp_lbc_levels
                               ! IN number of levels of DECOMP_LBC
INTEGER, INTENT(IN) ::  fld_type  
                               ! IN Which fld_type is the LBC?
INTEGER, INTENT(IN) ::  halo_type 
                               ! IN Which halo_type is the LBC?
INTEGER, INTENT(IN) ::  rim_type
                               ! IN Which rim_type is the LBC?
INTEGER, INTENT(IN) ::  pe_for_level(decomp_lbc_levels)
                               ! IN which level of FULL_LBC is on
                               !    which processor

REAL, INTENT(IN) ::  full_lbc(full_lbc_size,full_lbc_levels)
                               ! IN Some levels of the full LBC
REAL, INTENT(OUT) :: decomp_lbc(decomp_lbc_size,decomp_lbc_levels)
                               ! OUT All levels of the decomposed LBC on
                               !     this processor

INTEGER, INTENT(INOUT) :: icode  ! Return code not used at present.

CHARACTER(LEN=80) ::  cmessage ! OUT Error message


! Parameters and COMMON

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

! Local variables

INTEGER ::                                                        &
  iproc                                                           &
                     ! loop counter for loop over processors
, iside                                                           &
                     ! loop counter for loop over sides
, k                                                               &
                     ! loop counter for levels
, full_lbc_row_len                                                &
                     ! length of a row of data on the full LBC
, full_lbc_nrows                                                  &
                     ! number of rows of data on the full LBC
, decomp_lbc_row_len                                              &
                     ! length of a row of data on the
                     ! decomposed LBC
, decomp_lbc_nrows                                                &
                     ! number of rows of data on the
                     ! decomposed LBC
, first_lbc_pt                                                    &
                     ! first point in full LBC to start
                     ! copying from
, first_lbc_row                                                   &
                     ! first row in fill LBC to start
                     ! copying from
, level_index_pe(0:nproc-1)                                       &
                     ! How many levels on each PE
, level_index(decomp_lbc_levels)                                  &
                     ! Which level full_LBC corresponds
                     ! to the real level
, full_lbc_start_pt                                               &
                     ! First point index on a level of the
                     ! full LBC to start sending
, decomp_lbc_start_pt                                             &
                     ! First point index on a level of the
                     ! decomposed LBC to start receiving
, info                                                            
                     ! GCOM return code

INTEGER :: g_decomp_lbc_size(0:nproc-1)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! New data buffer to fill the decomposition to scatter

INTEGER :: size_data_lbc, data_lbc_start
REAL, ALLOCATABLE :: data_lbc(:)

! Control data needed for the scattering
INTEGER :: recvcount
INTEGER :: displs(0:nproc-1), sendcounts(0:nproc-1)
INTEGER :: my_comm ! communicator

! Other variables

INTEGER :: i, j, index_s, index_d

!--------------------------------------------------------------------
!--------------------------------------------------------------------
! 1.0 Set up indexing describing where each level of full LBC data is
!     held

IF (lhook) CALL dr_hook('SCATTER_ATMOS_LBCS',zhook_in,zhook_handle)
icode = 0
g_decomp_lbc_size=0
g_decomp_lbc_size(mype)=decomp_lbc_size
CALL gc_imax(nproc,nproc,info,g_decomp_lbc_size(0:nproc-1))

level_index_pe(:)=0

! Calculate indexing of full LBC data on the processors.
! A given level "k" of full LBC data will be held on
! processor "PE_FOR_LEVEL(k)" and the index of the level
! on this processor will be "level_index(k)"

DO k=1,decomp_lbc_levels
  level_index_pe(pe_for_level(k))=                                &
    level_index_pe(pe_for_level(k))+1
  level_index(k)=level_index_pe(pe_for_level(k))
END DO

IF (mype == 0) THEN

! Allocate the buffer

  size_data_lbc = SUM(g_decomp_lbc_size) * full_lbc_levels
  allocate(data_lbc(size_data_lbc))
  displs(0)=0
  DO iproc=1, nproc-1
    displs(iproc)=displs(iproc-1) +            &
      g_decomp_lbc_size(iproc-1) * full_lbc_levels
  END DO

  DO iproc=0, nproc-1
    DO iside=1,4
      IF (g_at_extremity(iside,iproc)) THEN

! This processor is at edge type iside and so needs LBC data
! In order to copy data from full_lbc array to data_lbc ( that
! contains data to be scattered ) the following variables need to
! be calculated

! full_lbc_row_len : East-West dimension of the full lbc side
! full_lbc_nrows : North-South dimension of the full lbc side
! decomp_lbc_row_len : East-West dimension of decomp lbc side
! decomp_lbc_nrows : North-South dimension of decomp lbc side
! full_lbc_start_pt : First point of the decomposed lbc side inside
!                     the 1d full_lbc array
! decomp_lbc_start_pt : First point of the decomposed lbc side inside
!                       the 1d decomposed lbc array
! data_lbc_start : First point inside the data lbc array ( corresponds
!                  to the first point of the distributed lbc_comp )
        IF ((iside == pnorth).OR.(iside==psouth)) THEN
                ! Calculate size of FULL_LBC

          IF (.NOT. l_vatpoles) THEN
          IF (fld_type  ==  fld_type_u) THEN
            full_lbc_row_len=glsize(1,fld_type)-1
          ELSE
            full_lbc_row_len=glsize(1,fld_type)
          END IF
          ELSE
          full_lbc_row_len=glsize(1,fld_type) 
          END IF  ! vatpoles

          full_lbc_row_len=full_lbc_row_len+2*halosize(1,halo_type)

          full_lbc_nrows=halosize(2,halo_type)+rimwidtha(rim_type)

! Calculate size of DECOMP_LBC

          IF (.NOT. l_vatpoles) THEN
          IF ((fld_type  ==  fld_type_u) .AND.                      &
            (g_at_extremity(peast,iproc))) THEN
            decomp_lbc_row_len=                                     &
            g_lasize(1,fld_type,halo_type,iproc)-1
          ELSE
            decomp_lbc_row_len=                                     &
            g_lasize(1,fld_type,halo_type,iproc)
          END IF
          ELSE
          decomp_lbc_row_len=                                       & 
          g_lasize(1,fld_type,halo_type,iproc) 
          END IF  ! vatpoles

          decomp_lbc_nrows=halosize(2,halo_type)+                   &
                         rimwidtha(rim_type)

! Calculate first point of DECOMP_LBC in FULL_LBC

          first_lbc_pt=g_datastart(1,iproc)
          first_lbc_row=1

          full_lbc_start_pt=                                          &
            global_lbc_starta(iside,fld_type,halo_type,rim_type)+     &
            (first_lbc_row-1)*full_lbc_row_len +                      &
            first_lbc_pt-1

          decomp_lbc_start_pt=                                        &
          g_lbc_starta(iside,fld_type,halo_type,rim_type,iproc)


        ELSE IF ((iside==peast).OR.(iside==pwest)) THEN

! Calculate size of FULL_LBC

          full_lbc_row_len=halosize(1,halo_type)+                   &
                         rimwidtha(rim_type)
          full_lbc_nrows=glsize(2,fld_type)-2*rimwidtha(rim_type)

        ! Calculate size of DECOMP_LBC

          decomp_lbc_row_len=halosize(1,halo_type)+                 &
                             rimwidtha(rim_type)
          decomp_lbc_nrows=g_lasize(2,fld_type,halo_type,iproc)
          IF (g_at_extremity(pnorth,iproc))                         &
            decomp_lbc_nrows=decomp_lbc_nrows-                      &
                             halosize(2,halo_type)-               &
                             rimwidtha(rim_type)
          IF (g_at_extremity(psouth,iproc))                         &
            decomp_lbc_nrows=decomp_lbc_nrows-                      &
                             halosize(2,halo_type)-               &
                             rimwidtha(rim_type)

        ! Calculate first point of DECOMP_LBC in FULL_LBC

          first_lbc_pt=1
          first_lbc_row=g_datastart(2,iproc)
          IF (.NOT. g_at_extremity(psouth,iproc))                   &
            first_lbc_row=first_lbc_row-                            &
                          rimwidtha(rim_type)-                    &
                          halosize(2,halo_type)

          full_lbc_start_pt=                                          &
            global_lbc_starta(iside,fld_type,halo_type,rim_type)+     &
            (first_lbc_row-1)*full_lbc_row_len +                      &
            first_lbc_pt-1

          decomp_lbc_start_pt=                                        &
          g_lbc_starta(iside,fld_type,halo_type,rim_type,iproc)

        END IF ! N/S or E/W test       

        data_lbc_start = displs(iproc)

! Now we can do the copy

        DO k=1, decomp_lbc_levels
          DO j=1, decomp_lbc_nrows
            DO i=1, decomp_lbc_row_len
               index_s=full_lbc_start_pt + (i-1) + (j-1) * full_lbc_row_len
               index_d=data_lbc_start +                              &
                 decomp_lbc_start_pt +                       &
                 (i-1) + (j-1) * decomp_lbc_row_len +                &
                 (k-1) * g_decomp_lbc_size (iproc)
               data_lbc(index_d)=full_lbc(index_s,k)
            END DO
          END DO
        END DO

      END IF ! At extremity
    END DO ! iside loop       
  END DO ! iproc loop

ELSE
  ALLOCATE(data_lbc(1))
END IF ! Master PE Section

! And Now we can do the scatter

CALL gc_get_communicator(my_comm, info)

sendcounts(:)=g_decomp_lbc_size(:) * decomp_lbc_levels
recvcount=decomp_lbc_size*decomp_lbc_levels

CALL mpl_scatterv(data_lbc, sendcounts, displs, mpl_real, &
                decomp_lbc, recvcount, mpl_real, 0, my_comm, info)

deallocate(data_lbc)

IF (lhook) CALL dr_hook('SCATTER_ATMOS_LBCS',zhook_out,zhook_handle)
RETURN
END SUBROUTINE scatter_atmos_lbcs
