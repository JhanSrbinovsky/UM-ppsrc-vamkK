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
! Purpose: Reads species and reaction data. Combines reactions into one
!          array and reorders them to put single reactant reactions first
!          to improve code in the prls routine.
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
!     Reads the information:
!                  - Chosen chemistry , contains information
!                    on species involved in the chemistry and
!                    the families/tracers to which they belong.
!                  - Data for bimolecular reactions.
!                  - Data for trimolecular reactions.
!                  - Data for photolysis reactions.
!                  - Data for heterogeneous reactions.
!     from ukca_chem1 module
!
!     Method
!     ------
!     The file specifies the species types using 2 letter
!     codes for easier reading.
!
!             ctype         Meaning
!             'FM'          Family member
!             'FT'          Tracer but will be put into a family
!                           if lifetime becomes short.
!             'TR'          Tracer, advected by calling model.
!             'SS'          Steady state species.
!             'CT'          Constant species.
!
!     Local variables
!     ---------------
!     iadv         Counter of number of model tracers/families.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_INRATS

        USE ASAD_MOD
        USE UKCA_CHEM_DEFS_MOD
        USE UKCA_D1_DEFS,         ONLY: n_chem_diags
        USE ukca_option_mod,      ONLY: jpctr, jpspec, jpbk, jptk, jphk, &
                                        jppj, jpnr 
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE ereport_mod, ONLY : ereport
        USE PrintStatus_mod
        USE UM_ParVars
        USE Control_Max_Sizes
        USE domain_params
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

        INTEGER :: errcode                ! Variable passed to ereport

        INTEGER :: ispb(jpbk+1,jpspb)
        INTEGER :: ispt(jptk+1,jpspt)
        INTEGER :: ispj(jppj+1,jpspj)
        INTEGER :: isph(jphk+1,jpsph)
        INTEGER :: ifrpbx(jpbk+1)
        INTEGER :: ifrpjx(jppj+1)
        INTEGER :: ifrptx(jptk+1)
        INTEGER :: ifrphx(jphk+1)
        INTEGER :: ilmin(jpspec)
        INTEGER :: ilmaj(jpspec)
        INTEGER :: ixdumt(jptk+1)
        INTEGER :: ixdumj(jppj+1)
        INTEGER :: ixdumh(jphk+1)
        INTEGER :: ierror
        INTEGER :: iadv                  ! Counter for advected species
        INTEGER :: inadv                 ! Counter for non-advected species
        INTEGER :: imajor                ! Counter
        INTEGER :: iminor                ! Counter
        INTEGER :: ix                    ! Counter
        INTEGER :: icount                ! Counter
        INTEGER :: ifam                  ! Index
        INTEGER :: imaj                  ! Index
        INTEGER :: iflag                 ! Used to test family order
        INTEGER :: idummy                ! Dummy variable
        INTEGER :: istat                 ! Tag for communication
        INTEGER :: j                     ! Loop variable
        INTEGER :: jf                    ! Loop variable
        INTEGER :: jadv                  ! Loop variable
        INTEGER :: jb                    ! Loop variable
        INTEGER :: jctr                  ! Loop variable
        INTEGER :: jh                    ! Loop variable
        INTEGER :: jj                    ! Loop variable
        INTEGER :: jp                    ! Loop variable
        INTEGER :: jr                    ! Loop variable
        INTEGER :: js                    ! Loop variable
        INTEGER :: jspb                  ! Loop variable
        INTEGER :: jsph                  ! Loop variable
        INTEGER :: jspj                  ! Loop variable
        INTEGER :: jspt                  ! Loop variable
        INTEGER :: jt                    ! Loop variable
        INTEGER :: k                     ! Loop variable
        INTEGER :: ind                   ! Loop index

        REAL :: zdumt(1)                 ! Dummy array for trimol file
        REAL :: zdumj(1)                 ! Dummy array for photol file
        REAL :: zdumh(1)                 ! Dummy array for heter file

        CHARACTER (LEN=10) :: cmntb(jpbk+1)
        CHARACTER (LEN=10) :: cmntt(jptk+1)
        CHARACTER (LEN=10) :: cmntj(jppj+1)
        CHARACTER (LEN=10) :: cmnth(jphk+1)
        CHARACTER (LEN=10), PARAMETER :: nullx='          '
        CHARACTER (LEN=72) :: cmessage        ! Error message

        LOGICAL :: gtype3
        LOGICAL :: gdebug
        LOGICAL :: L_exist
        LOGICAL :: L_fa

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


!       1.  Read chosen chemistry file and determine chemistry
!           ---- ------ --------- ---- --- --------- ---------

      IF (lhook) CALL dr_hook('ASAD_INRATS',zhook_in,zhook_handle)

! Initialise local counters, see asad_cinit for (e.g.) spb and frpb 
        ispb(:,:) = 0
        ifrpbx(:) = 0
        ispt(:,:) = 0
        ifrptx(:) = 0
        ispj(:,:) = 0
        ifrpjx(:) = 0
        isph(:,:) = 0
        ifrphx(:) = 0

!       1.1.1  Set types and find which species are in families.

! Replace read statement by module:

        IF (size(chch_defs) /= jpspec) THEN
          errcode=1
          cmessage=' jpspec and chch_defs are inconsistent'
          WRITE(6,*) cmessage

          CALL EREPORT('ASAD_INRATS',errcode,cmessage)
        ENDIF
        DO k=1,jpspec
          speci(k)  = chch_defs(k)%speci
          nodd(k)   = chch_defs(k)%nodd
          ctype(k)  = chch_defs(k)%ctype
          family(k) = chch_defs(k)%family
        ENDDO

        iadv = 0
        inadv = 0
        ntrf = 0
        nnaf = 0
        ntr3 = 0
        DO js = 1, jpspec
          gtype3 = .false.
          IF ( ctype(js) /= jpfm .AND. ctype(js) /= jpif )              &
            family(js)='          '
          IF ( ctype(js) == jpsp ) THEN               ! tracers
            iadv = iadv + 1
            ntrf = ntrf + 1
            IF ( iadv > jpctr ) THEN
              WRITE (6,*) '** ASAD ERROR in subroutine inrats'
              WRITE (6,*) '** Parameter jpctr is too low; found',iadv
              WRITE (6,*) '** tracers so far with ',jpspec-js
              WRITE (6,*) '** species to check.'
              cmessage = 'ASAD ERROR: jpctr is too low'

              CALL EREPORT('ASAD_INRATS',iadv,cmessage)
            ENDIF
            advt(iadv)  = speci(js)
            nltrf(ntrf) = iadv
          ELSE IF( ctype(js) == jpna )THEN        ! non-advected species 
            inadv = inadv + 1 
            nnaf = nnaf + 1 
            IF ( inadv > n_chem_diags ) THEN 
              cmessage = 'ASAD ERROR: Too many non-advected tracers' 
 
              CALL EREPORT('ASAD_INRATS',iadv,cmessage) 
            ENDIF 
            nadvt(inadv)  = speci(js) 
          ELSE IF( ctype(js) == jpfm .OR. ctype(js) == jpif )THEN
            cmessage = ' Family chemistry not available in this version'
            CALL EREPORT('ASAD_INRATS',js,cmessage)
            IF ( ctype(js) == jpif ) THEN
              iadv = iadv + 1
              ntr3 = ntr3 + 1
              IF ( iadv  >   jpctr ) THEN
                WRITE (6,*) '** ASAD ERROR in subroutine inrats'
                WRITE (6,*) '** Parameter jpctr is too low; found',iadv
                WRITE (6,*) '** tracers so far with ',jpspec-js
                WRITE (6,*) '** species to check.'
                cmessage = 'ERROR in jpctr'

                CALL EREPORT('ASAD_INRATS',js,cmessage)
              ENDIF
              advt(iadv)  = speci(js)
              nltr3(ntr3) = iadv
            ENDIF
            l_fa=.true.
            DO jadv = 1, iadv
              IF ( family(js) == advt(jadv) ) THEN
                l_fa=.false.
                EXIT
              ENDIF
            ENDDO
            IF (l_fa) THEN
              iadv = iadv + 1
              ntrf = ntrf + 1
              IF ( iadv  >   jpctr ) then
                WRITE (6,*) '***** ASAD ERROR in subroutine inrats'
                WRITE (6,*) '***** Param jpctr is too low; found',iadv
                WRITE (6,*) '***** tracers so far with ',jpspec-js
                WRITE (6,*) '***** species to check.'
                cmessage = 'INRATS ERROR : jpctr is too low'

                CALL EREPORT('ASAD_INRATS',iadv,cmessage)
              ENDIF
              advt(iadv) = family(js)
              nltrf(ntrf) = iadv
            ENDIF      ! l_fa
          ENDIF
        ENDDO

!       1.2 Find major species of families

        DO jadv = 1, iadv
          DO js = 1, jpspec
            IF( family(js) == advt(jadv) .and. ctype(js) /= jpif )     &
              majors(jadv) = js
            IF( speci(js) == advt(jadv) ) majors(jadv) = js
          ENDDO
        ENDDO

!       1.3 Allocate families to species

        DO js = 1, jpspec
          moffam(js) = 0
          madvtr(js) = 0
          DO jadv = 1, iadv
            IF (family(js) == advt(jadv) ) moffam(js) = jadv
            IF (speci(js)  == advt(jadv) ) madvtr(js) = jadv
            IF (family(js) == advt(jadv) .AND. js > majors(jadv)) THEN
              WRITE (6,*) '** ASAD ERROR: '
              WRITE (6,*) 'RE-ORDER SPECIES FILE SO THAT THE MAJOR '
              WRITE (6,*) 'SPECIES OF A FAMILY OCCURS AFTER THE OTHERS'
              cmessage = 'INRATS ERROR : Order of species is incorrect'

              CALL EREPORT('ASAD_INRATS',jadv,cmessage)
            ENDIF
          ENDDO
        ENDDO

!       1.4  Build the list of major and minor species

        nlmajmin(1) = 5
        imajor      = 0
        iminor      = 0
        DO js = 1, jpspec
          ifam = moffam(js)
          imaj = 0
          IF ( ifam /= 0 ) THEN
           imaj = majors(ifam)
           IF ( imaj /= js ) THEN
              iminor = iminor + 1
              ilmin(iminor) = js
            ELSE
              imajor = imajor + 1
              ilmaj(imajor) = js
            ENDIF
          ENDIF
        ENDDO
        nlmajmin(2) = nlmajmin(1) + imajor - 1
        nlmajmin(3) = nlmajmin(1) + imajor
        nlmajmin(4) = nlmajmin(3) + iminor - 1
        DO j = nlmajmin(1), nlmajmin(2)
          nlmajmin(j) = ilmaj(j-nlmajmin(1)+1)
        ENDDO
        DO j = nlmajmin(3), nlmajmin(4)
          nlmajmin(j) = ilmin(j-nlmajmin(3)+1)
        ENDDO

        IF ( iadv /= jpctr ) then
          WRITE (6,*) '** ASAD ERROR: Number of advected tracers',     &
           ' specified in chch_defs does not match jpctr'
          WRITE (6,*) 'Found ',iadv,' but expected ',jpctr
          cmessage = 'INRATS ERROR : iadv and jpctr do not match'

          CALL EREPORT('ASAD_INRATS',iadv,cmessage)
        ENDIF

!       2.  Write details of chemistry selection to log file
!           ----- ------- -- --------- --------- -- --- ----

        IF (mype == 0 .AND. printstatus >= prstatus_oper) THEN
          WRITE(6,*)
          WRITE(6,*)'  ***  CHEMISTRY INFORMATION  ***'
          WRITE(6,*)
          WRITE(6,*)'ASAD IS TREATING ADVECTED TRACERS IN THE ORDER:'
          WRITE(6,*)
          WRITE(6,'(5(2x,i2,1x,a10))')(jctr,advt(jctr), jctr= 1,jpctr)
          WRITE(6,*)
          WRITE(6,*)'ASAD IS TREATING NON-ADVECTED TRACERS IN THE ORDER:' 
          WRITE(6,*) 
          WRITE(6,'(5(2x,i2,1x,a10))')(jctr,nadvt(jctr), jctr= 1,nnaf) 
          WRITE(6,*) 
          WRITE(6,*)'IF THE TRACERS WERE NOT INITIALISED IN THIS '
          WRITE(6,*)'ORDER THEN THE MODEL RESULTS ARE WORTHLESS '
          WRITE(6,*)

          iflag = 0
          DO jctr = 1, jpctr
            IF (advt(jctr) /= speci(majors(jctr)) ) THEN
              IF (iflag == 0 )then
                WRITE(6,*)'THE MAJOR MEMBER OF EACH OF THE FAMILIES'
                WRITE(6,*)'IS GIVEN BELOW. IF THIS IS NOT ACCEPTABLE,'
                WRITE(6,*)'THEN YOU MUST REORDER THE SPECIES IN'//      &
                          ' CHCH_DEFS'
                WRITE(6,*)'SO THE MAJOR SPECIES FOLLOWS THE OTHERS.'
                WRITE(6,*)
                iflag = 1
              ENDIF
              WRITE(6,'(a10,1x,a10)') advt(jctr), speci(majors(jctr))
            ENDIF
          ENDDO
        ENDIF     ! End of IF mype statement

!       3.  Bimolecular ratefile
!           ----------- --------

!       Get bimolecular rates from module

        IF (size(ratb_defs) /= jpbk) THEN 
          errcode=1
          cmessage='size of ratb_defs is inconsistent with jpbk'

          CALL EREPORT('ASAD_INRATS',errcode,cmessage)
        END IF 
        icount=1
        DO k=1,jpbk
          spb(k,1) = ratb_defs(k)%react1
          spb(k,2) = ratb_defs(k)%react2
          spb(k,3) = ratb_defs(k)%prod1
          spb(k,4) = ratb_defs(k)%prod2
          spb(k,5) = ratb_defs(k)%prod3
          spb(k,6) = ratb_defs(k)%prod4
          ab(k,1)  = ratb_defs(k)%K0
          ab(k,2)  = ratb_defs(k)%alpha
          ab(k,3)  = ratb_defs(k)%beta
          IF (ratb_defs(k)%pyield1 > 1e-18) THEN
            ifrpbx(k)     = icount
            frpb(icount)  = ratb_defs(k)%pyield1
            IF (spb(k,4) /= nullx) frpb(icount+1) = ratb_defs(k)%pyield2
            IF (spb(k,5) /= nullx) frpb(icount+2) = ratb_defs(k)%pyield3
            IF (spb(k,6) /= nullx) frpb(icount+3) = ratb_defs(k)%pyield4
            icount = icount + 4
          ENDIF
        ENDDO

        DO jb = 1, jpbk
          DO js = 1, jpspec
            DO jspb = 1, jpspb
              IF ( speci(js) == spb(jb,jspb) ) ispb(jb,jspb) = js
            ENDDO
          ENDDO
        ENDDO

! Load in bimol fractional prod coefs to frpx array
        DO jf = 1, jpfrpb
          frpx(jf) = frpb(jf)
        END DO

!       4.  Trimolecular ratefile
!           ------------ --------

!       Get trimolecular rates from module for UM version

        IF (size(ratt_defs) /= jptk) THEN
          errcode=1
          cmessage='size of ratt_defs is inconsistent with jptk'

          CALL EREPORT('ASAD_INRATS',errcode,cmessage)
        ENDIF
        icount=1
        DO k=1,jptk
          spt(k,1) = ratt_defs(k)%react1
          spt(k,2) = ratt_defs(k)%react2
          spt(k,3) = ratt_defs(k)%prod1
          spt(k,4) = ratt_defs(k)%prod2
          at(k,1)  = ratt_defs(k)%F
          at(k,2)  = ratt_defs(k)%K1
          at(k,3)  = ratt_defs(k)%alpha1
          at(k,4)  = ratt_defs(k)%beta1
          at(k,5)  = ratt_defs(k)%K2
          at(k,6)  = ratt_defs(k)%alpha2
          at(k,7)  = ratt_defs(k)%beta2
          IF (ratt_defs(k)%pyield1 > 1e-18) THEN
            ifrptx(k)      = icount
            frpt(icount)   = ratt_defs(k)%pyield1
            IF (spt(k,4) /= nullx) frpt(icount+1) = ratt_defs(k)%pyield2
            icount = icount + 2
          ENDIF
        ENDDO

        DO jt = 1, jptk
          DO js = 1, jpspec
            DO jspt = 1, jpspt
              IF (speci(js) == spt(jt,jspt) ) ispt(jt,jspt) = js
            ENDDO
          ENDDO
        ENDDO

! Load in trimol fractional prod coefs to frpx array
        DO jf = 1, jpfrpt
          ind = jf + jpfrpb
          frpx(ind) = frpt(jf)
        END DO


!       5.  Photolysis ratefile
!           ---------- --------

!       use module to get spj

        IF (size(ratj_defs) /= jppj) THEN
          errcode=1
          cmessage='size of ratj_defs is not equal to jppj'

          CALL EREPORT('ASAD_INRATS',errcode,cmessage)
        ENDIF
        icount=1
        DO k=1,jppj
          spj(k,1) = ratj_defs(k)%react1
          spj(k,2) = ratj_defs(k)%react2
          spj(k,3) = ratj_defs(k)%prod1
          spj(k,4) = ratj_defs(k)%prod2
          spj(k,5) = ratj_defs(k)%prod3
          spj(k,6) = ratj_defs(k)%prod4
          IF (ratj_defs(k)%pyield1 > 1e-18) THEN
            ifrpjx(k)     = icount
            frpj(icount)  = ratj_defs(k)%pyield1
            IF (spj(k,4) /= nullx) frpj(icount+1) = ratj_defs(k)%pyield2
            IF (spj(k,5) /= nullx) frpj(icount+2) = ratj_defs(k)%pyield3
            IF (spj(k,6) /= nullx) frpj(icount+3) = ratj_defs(k)%pyield4
            icount = icount + 4
          ENDIF
        ENDDO

        DO jj = 1, jppj
          DO js = 1, jpspec
            DO jspj = 1, jpspj
              IF (speci(js) == spj(jj,jspj) ) ispj(jj,jspj) = js
            ENDDO
          ENDDO
        ENDDO

! Load in photol fractional prod coefs to frpx array
        DO jf = 1, jpfrpj
          ind = jf + jpfrpb + jpfrpt
          frpx(ind) = frpj(jf)
        END DO

!       6.  Heterogeneous ratefile
!           ------------- --------

        IF (jphk > 0) THEN

!         use module to get sph

          IF (size(rath_defs) /= jphk) THEN
            errcode=1
            cmessage='size of rath_defs is not equal to jphk'

            CALL EREPORT('ASAD_INRATS',errcode,cmessage)
          ENDIF
          icount=1
          DO k=1,jphk
            sph(k,1) = rath_defs(k)%react1
            sph(k,2) = rath_defs(k)%react2
            sph(k,3) = rath_defs(k)%prod1
            sph(k,4) = rath_defs(k)%prod2
            sph(k,5) = rath_defs(k)%prod3
            sph(k,6) = rath_defs(k)%prod4
            IF (rath_defs(k)%pyield1 > 1e-18) THEN
              ifrphx(k)     = icount
              frph(icount)  = rath_defs(k)%pyield1
              IF (sph(k,4) /= nullx) frph(icount+1)=rath_defs(k)%pyield2
              IF (sph(k,5) /= nullx) frph(icount+2)=rath_defs(k)%pyield3
              IF (sph(k,6) /= nullx) frph(icount+3)=rath_defs(k)%pyield4
              icount = icount + 4
            ENDIF
          ENDDO

          DO jh = 1, jphk
            DO js = 1, jpspec
              DO jsph = 1, jpsph
                IF (speci(js) == sph(jh,jsph) ) isph(jh,jsph) = js
              ENDDO
            ENDDO
          ENDDO

! Load in het fractional prod coefs to frpx array
          DO jf = 1, jpfrph
            ind = jf + jpfrpb + jpfrpt + jpfrpj
            frpx(ind) = frph(jf)
          END DO

        END IF       ! jphk > 0

!       7.  Reorder reactions, putting single reactants first.
!           ------- ---------- ------- ------ --------- ------

        nuni = 0

!       7.1  Single reactants; scan ratefiles in turn.

!
        DO jr = 1, jpbk
          IF (ispb(jr,2) == 0 ) THEN
            nuni = nuni + 1
            nbrkx(jr) = nuni
            DO jp = 1, jpspb
              nspi(nuni,jp) = ispb(jr,jp)
            ENDDO
            IF ( ifrpbx(jr) /= 0 ) nfrpx(nuni) = ifrpbx(jr)
          ENDIF
        ENDDO

        DO jr = 1, jptk
          IF ( ispt(jr,2) == 0 ) THEN
            nuni = nuni + 1
            ntrkx(jr) = nuni
            DO jp = 1, jpspt
              nspi(nuni,jp) = ispt(jr,jp)
            ENDDO
            IF (ifrptx(jr) /= 0) nfrpx(nuni) = ifrptx(jr)+jpfrpb
          ENDIF
        ENDDO

        DO jr = 1, jppj
          IF ( ispj(jr,2) == 0 ) THEN
            nuni = nuni + 1
            nprkx(jr) = nuni
            DO jp = 1, jpspj
              nspi(nuni,jp) = ispj(jr,jp)
            ENDDO
            IF (ifrpjx(jr) /= 0) nfrpx(nuni) = ifrpjx(jr)+jpfrpb+jpfrpt
          ENDIF
        ENDDO

        IF ( jphk > 0 ) THEN
          DO jr = 1, jphk
            IF ( isph(jr,2) == 0 ) THEN
              nuni = nuni + 1
              nhrkx(jr) = nuni
              DO jp = 1, jpsph
                nspi(nuni,jp) = isph(jr,jp)
              ENDDO
              IF (ifrphx(jr) /= 0) nfrpx(nuni) = ifrphx(jr)+jpfrpb+         &
                                                 jpfrpt+jpfrpj
            ENDIF
          ENDDO
        ENDIF

!       7.2  Two reactants; copy remaining reactions

        ix = nuni
        DO jr = 1, jpbk
          IF ( ispb(jr,2) /= 0 ) THEN
            ix = ix + 1
            nbrkx(jr) = ix
            DO jp = 1, jpspb
              nspi(ix,jp) = ispb(jr,jp)
            ENDDO

            IF ( ifrpbx(jr) /= 0 ) nfrpx(ix) = ifrpbx(jr)
          ENDIF
        ENDDO

        DO jr = 1, jptk
          IF ( ispt(jr,2) /= 0 ) THEN
            ix = ix + 1
            ntrkx(jr) = ix
            DO jp = 1, jpspt
              nspi(ix,jp) = ispt(jr,jp)
            ENDDO
            IF (ifrptx(jr) /= 0) nfrpx(ix) = ifrptx(jr)+jpfrpb
          ENDIF
        ENDDO

        DO jr = 1, jppj
          IF (ispj(jr,2) /= 0 ) THEN
            ix = ix + 1
            nprkx(jr) = ix
            DO jp = 1, jpspj
              nspi(ix,jp) = ispj(jr,jp)
            ENDDO
            IF (ifrpjx(jr) /= 0) nfrpx(ix) = ifrpjx(jr)+jpfrpb+jpfrpt
          ENDIF
        ENDDO

        IF ( jphk > 0 ) THEN
          DO jr = 1, jphk
            IF ( isph(jr,2) /= 0 ) THEN
              ix = ix + 1
              nhrkx(jr) = ix
              DO jp = 1, jpsph
                nspi(ix,jp) = isph(jr,jp)
              ENDDO
              IF (ifrphx(jr) /= 0) nfrpx(ix) = ifrphx(jr)+jpfrpb+       &
                                               jpfrpt+jpfrpj
            ENDIF
          ENDDO
        ENDIF

        IF ( ix /= jpnr ) THEN
          WRITE (6,*) '*** INTERNAL ASAD ERROR: Number of reactions',   &
                      ' placed in nspi array does not equal jpnr. '
          WRITE (6,*) '                         Check that reaction',   &
                      ' files and value of jpnr in UKCA namelist are' 
          WRITE (6,*) '                         consistent. Found: ',   &
                      ix,' jpnr: ',jpnr
          cmessage = 'No of rxns in nspi array is not equal to jpnr'

          CALL EREPORT('ASAD_INRATS',ix,cmessage)
        ENDIF

        IF (lhook) CALL dr_hook('ASAD_INRATS',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_INRATS
