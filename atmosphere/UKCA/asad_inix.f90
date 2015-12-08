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
! Purpose: To set up the arrays holding the indicies of reactions,
!     species etc used in the main parts of the code e.g. prls.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CINIT
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!     Interface
!     ---------
!     On exit from this routine, the following arrays will have
!     been set: ngrp, nprdx3, nprdx2, nprdx1, nlall, nlstst, nlf,
!     nldepd, nldepw and nlemit.
!
!     IMPORTANT! This routine MUST be called after inrats since it
!     expects the nspi array to be initialised.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_INIX

        USE ASAD_MOD
        USE ukca_option_mod, ONLY: jpspec, jpnr
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE ereport_mod, ONLY : ereport
        USE UM_ParParams
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

!       Local variables

        INTEGER :: j                    ! Loop variable
        INTEGER :: jr                   ! Loop variable
        INTEGER :: js                   ! Loop variable
        INTEGER :: jp                   ! Loop variable
        INTEGER :: iprd
        INTEGER :: igrp3
        INTEGER :: iloss
        INTEGER :: ngmax
        INTEGER :: is
        INTEGER :: ispec
        INTEGER :: jdry
        INTEGER :: jwet
        INTEGER :: jpnpx3

        INTEGER :: igrp(3)
        INTEGER :: idepd(jpspec)
        INTEGER :: idepw(jpspec)

        CHARACTER (LEN=2)  :: itype
        CHARACTER (LEN=72) :: cmessage  ! Error message

        LOGICAL :: gdry
        LOGICAL :: gwet

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


!       1.  Initialise the types of specie arrays.
!           ---------- --- ----- -- ------ -------

        IF (lhook) CALL dr_hook('ASAD_INIX',zhook_in,zhook_handle)
        jpnpx3 = (jpnr/(3*3))+3*3
        nstst = 0
        nf = 0
        DO js = 1, jpspec
          itype = ctype(js)
          nlall(js) = js
          IF (itype==jpfm .OR. itype==jpif .OR. itype==jpna) THEN
            nstst = nstst + 1
            nlstst(nstst) = js
          END IF
          IF ( itype==jpfm .OR. itype==jpif .OR. itype==jpsp ) THEN
            nf = nf + 1
            nlf(nf) = js
          END IF
        END DO

!       2.  Fill in the 3, 2 & 1 term product index arrays.
!           ---- -- --- -- - - - ---- ------- ----- -------

        DO js = 1, jpspec

!         2.1  For each species, scan the reactions for groups
!              of 3, 2 & 1 product terms and loss terms.

          iprd  = 0
          igrp3 = 0
          DO jr = 1, jpnr
            IF ( nfrpx(jr) == 0 ) THEN

!             If reaction does NOT have fractional products.

              DO jp = 3, jpmsp
                IF ( nspi(jr,jp) == js ) THEN
                  iprd = iprd + 1
                  igrp(iprd) = jr
                  IF ( iprd == 3 ) THEN
                    iprd = 0
                    igrp3 = igrp3 + 1
                    DO j = 1, 3
                      nprdx3(j,igrp3,js) = igrp(j)
                    END DO
                  END IF
                END IF
              END DO
            END IF
          END DO
          ngrp(js,3) = igrp3
          IF ( iprd == 2 ) THEN
            ngrp(js,2) = 1
            nprdx2(1,js) = igrp(1)
            nprdx2(2,js) = igrp(2)
          ELSE IF ( iprd == 1 ) THEN
            ngrp(js,1) = 1
            nprdx1(js) = igrp(1)
          END IF

          iloss = 0
          igrp3 = 0
          DO jr = 1, jpnr
            DO jp = 1, 2
              IF ( nspi(jr,jp) == js ) THEN
                iloss = iloss + 1
                igrp(iloss) = jr
                IF ( iloss == 3 ) THEN
                  iloss = 0
                  igrp3 = igrp3 + 1
                  DO j = 1, 3
                    nprdx3(j,igrp3,jpspec+js) = igrp(j)
                  END DO
                END IF
              END IF
            END DO
          END DO

          ngrp(jpspec+js,3) = igrp3
          IF ( iloss == 2 ) THEN
            ngrp(js+jpspec,2) = 1
            nprdx2(1,jpspec+js) = igrp(1)
            nprdx2(2,jpspec+js) = igrp(2)
          ELSE IF ( iloss == 1 ) THEN
            ngrp(js+jpspec,1) = 1
            nprdx1(jpspec+js) = igrp(1)
          END IF

        END DO

!       2.2  Check array size is ok.

        ngmax = ngrp(1,3)
        DO js = 2, jpspec
          ngmax = max( ngmax, ngrp(js,3) )
        END DO
        IF ( ngmax > jpnpx3 ) THEN
           WRITE (6,*) '** ASAD ERROR: The parameter jpnpx3 is too',   &
         ' small. Please change it not less than ',ngmax
          cmessage = ' Parameter jpnpx3 is too small'

          CALL EREPORT('ASAD_INIX',ngmax,cmessage)
        END IF

!       2.3  Build table for fractional products

        nnfrp = 0
        DO jr = 1, jpnr
          IF ( nfrpx(jr) /= 0 ) THEN
            DO jp = 3, jpspb
              is = nspi(jr,jp)
              IF ( is /= 0 ) THEN
                nnfrp = nnfrp + 1
                IF ( nnfrp > jpfrpx ) THEN
                  WRITE (6,*) '*ASAD ERROR: The parameter jpfrpx is',  &
                 ' too small. Please increase it.'
                  cmessage=' Parameter jpfrpx is too small for no. '// &
                            'of fractional products'

                  CALL EREPORT('ASAD_INIX',nnfrp,cmessage)
                END IF
                ntabfp(nnfrp,1) = is
                ntabfp(nnfrp,2) = jr
                ntabfp(nnfrp,3) = nfrpx(jr) + (jp-3)
              END IF
            END DO
          END IF
        END DO

!       3.   Set up lists for deposition and emissions.
!            --- -- ----- --- ---------- --- ----------

        ndepd = 0
        ndepw = 0
        nemit = 0
        DO  js = 1, jpspec
          IF ( ldepd(js) ) THEN
            ndepd = ndepd + 1
            nldepd(ndepd) = js
          END IF
          IF ( ldepw(js) ) THEN
            ndepw = ndepw + 1
            nldepw(ndepw) = js
          END IF
          IF ( lemit(js) ) THEN
            nemit = nemit + 1
            nlemit(nemit) = js
          END IF
        END DO

!       3.1.  Now form the list used in prls.f

        DO js = 1, jpspec
          idepd(js) = nldepd(js)
          idepw(js) = nldepw(js)
        END DO

!       First look for species with both dry and wet on.

        nldepx(1) = 7
        nldepx(2) = nldepx(1) - 1

        DO ispec = 1,jpspec
          gdry = .false.
          gwet = .false.

          DO j = 1, jpspec
            IF ( idepd(j) == ispec ) THEN
              gdry = .true.
              jdry = j
            END IF
            IF ( idepw(j) == ispec ) THEN
              gwet = .true.
              jwet = j
            END IF
          END DO

          IF ( gdry .AND. gwet ) THEN
            nldepx(2)         = nldepx(2) + 1
            nldepx(nldepx(2)) = ispec
            idepd(jdry)       = 0
            idepw(jwet)       = 0
          END IF
        END DO   ! ispec loop

!       Now get the remaining species with dry dep. only

        nldepx(3) = nldepx(2) + 1
        nldepx(4) = nldepx(3) - 1
        DO j = 1, jpspec
          IF ( idepd(j) /= 0 ) THEN
            nldepx(4) = nldepx(4) + 1
            nldepx(nldepx(4)) = idepd(j)
          END IF
        END DO

!       Remaining species with wet deposition on.

        nldepx(5) = nldepx(4) + 1
        nldepx(6) = nldepx(5) - 1
        DO j = 1, jpspec
          IF ( idepw(j) /= 0 ) THEN
            nldepx(6) = nldepx(6) + 1
            nldepx(nldepx(6)) = idepw(j)
          END IF
        END DO

        IF (lhook) CALL dr_hook('ASAD_INIX',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_INIX
