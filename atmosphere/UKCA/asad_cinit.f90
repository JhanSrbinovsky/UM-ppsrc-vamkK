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
! Purpose: To initialize variables used in the chemistry
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_CHEMISTR_CTL before call to 
!          ASAD_CDRIVE
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!     Method
!     ------
!     Input arguments are checked and copied to common. Species
!     and reaction data are read in. Other variables used in the
!     chemistry are initialised.
!
!     -- Calculation of peps
!
!     At several places in the ASAD code, we needed a min. value
!     to use to guard against zero divides. We have to compute
!     this to allow for all possible computer hardwares and precisions
!     that ASAD might be run at.
!
!     Externals
!     ---------
!     inrats      - Reads species and reaction data.
!     inphot      - Initialises photolysis scheme.
!     ukca_inwdep - Initialises wet deposition data
!                   (user supplied routine)
!     ukca_inddep - Initialises dry deposition data
!                   (user supplied routine)
!     inemit      - Initialises emission data
!                   (user supplied routine).
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_CINIT(p_field)

        USE ASAD_MOD
        USE ukca_option_mod, ONLY: jpnr, jpctr, jpspec
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE ereport_mod, ONLY : ereport
        USE UM_ParParams
        USE Control_Max_Sizes
        USE domain_params
        USE PrintStatus_mod
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

        INTEGER, INTENT(IN) :: p_field

!       Local variables

        INTEGER :: j                      ! Loop variable
        INTEGER :: jc                     ! Loop variable
        INTEGER :: jf                     ! Loop variable
        INTEGER :: jg                     ! Loop variable
        INTEGER :: jl                     ! Loop variable
        INTEGER :: jp                     ! Loop variable
        INTEGER :: jr                     ! Loop variable
        INTEGER :: js                     ! Loop variable
        INTEGER :: jtr                    ! Loop variable
        INTEGER :: jx                     ! Loop variable
        INTEGER :: jpnpx3                 ! Loop variable
        INTEGER :: errcode                ! Variable passed to ereport
        INTEGER, PARAMETER :: nrsteps_max = 200  ! max steps
        INTEGER, PARAMETER :: nit0_max    = 50   ! max
        INTEGER, PARAMETER :: nitfg_max   = 50   ! max
        INTEGER, PARAMETER :: nitnr_max   = 50   ! max

        REAL :: sfmin

        CHARACTER (LEN=10), PARAMETER :: nullx='          '
        CHARACTER (LEN=72) :: cmessage    ! Error message

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


!       Initialisations for COMMON variables.

        IF (lhook) CALL dr_hook('ASAD_CINIT',zhook_in,zhook_handle)
        ljacx = .true.
        
        jpnpx3 = (jpnr/(3*3))+3*3


!       1.  Copy input arguments to COMMON and check them.
!           ---- ----- --------- -- ------ --- ----- -----

!       Logical arguments
        DO j = 1, jpspec
          lemit(j) = .false.
        ENDDO

        IF (nrsteps < 0 .OR. nrsteps > 200) THEN
          cmessage = ' NRSTEPS IS OUT OF RANGE, RESETTING'
          write(6,*) 'NRSTEPS = ',nrsteps,' Reset to: ',nrsteps_max
          nrsteps = nrsteps_max
          errcode=-1

          CALL EREPORT('ASAD_CINIT',errcode,cmessage)
        ENDIF

        IF (nit0 < 0 .OR. nit0 > 50) THEN
          cmessage = ' NIT0 IS OUT OF RANGE, RESETTING'
          write(6,*) 'NIT0 = ',nit0,' Reset to: ',nit0_max
          nit0 = nit0_max
          errcode=-1

          CALL EREPORT('ASAD_CINIT',errcode,cmessage)
        ENDIF

        IF (nitfg < 0 .OR. nitfg > 50) THEN
          cmessage = ' NITFG IS OUT OF RANGE, RESETTING'
          write(6,*) 'NITFG = ',nitfg,' Reset to: ',nitfg_max
          nitfg = nitfg_max
          errcode=-1

          CALL EREPORT('ASAD_CINIT',errcode,cmessage)
       ENDIF

       IF (nitnr < 0 .OR. nitnr > 50 )then
         cmessage = ' NITNR IS OUT OF RANGE, RESETTING'
         write(6,*) 'NITNR = ',nitnr,' Reset to: ',nitnr_max
         nitnr = nitnr_max
         errcode=-1

         CALL EREPORT('ASAD_CINIT',errcode,cmessage)
       ENDIF

!      2.  Assign chemical timestep values for ASAD
!          ----------------------------------------

! Set up timestep counting. Interval depends on solver: IMPACT and Rosenbrock
! are run every dynamical timestep. Newton-Raphson is run every 2 / 3
! timesteps for a 30/20 minutes dynamical timestep.

        IF (method == 3) THEN          ! N-R solver
          ncsteps = 1
          cdt = REAL(kcdt)
          interval = NINT(cdt/dtime)
        ELSE IF (method == 1) THEN     ! IMPACT
! use about 15 or 10 minutes, depending on dynamical timestep
          IF (dtime < tslimit) THEN
            cdt = dtime
            ncsteps = 1
            interval = 1
          ELSE
            cdt = 0.5*dtime
            ncsteps = 2
            interval = 1
          END IF
        ELSE                           ! Other solvers
          cdt = dtime
          ncsteps = 1
          interval = 1
        END IF

        IF (printstatus >= prstatus_oper) THEN
          WRITE(6,*) 'Interval for chemical solver set to: ', interval
          WRITE(6,*) 'Timestep for chemical solver set to: ', cdt
          WRITE(6,*) 'No steps for chemical solver set to: ', ncsteps
        END IF

        IF (ABS(cdt*ncsteps - dtime*interval) > 1e-4) THEN
          cmessage=' chemical timestep does not fit dynamical timestep'
          errcode = kcdt
          CALL EREPORT('ASAD_CINIT',errcode,cmessage)
        END IF

!       2.1  Set photolysis frequency.

        IF (kfphot < 0 .AND. abs(kfphot) > dtime) THEN
          write (6,*) '**CINIT WARNING: VALUE OF KFPHOT ',kfphot,      &
          'EXCEEDS THE MODEL TIMESTEP. ROUTINE PHOTOL WILL ONLY',      &
          ' BE CALLED ONCE.'
          nfphot = 0
        ELSEIF ( kfphot > 0 .AND. kfphot > ncsteps ) THEN
          write (6,*) '**CINIT WARNING: FREQUENCY KFPHOT ',kfphot,     &
           ' EXCEEDS THE TOTAL NUMBER OF CHEMICAL SUBSTEPS. ROUTINE ', &
           ' PHOTOL WILL BE CALLED ONCE ONLY.'
          nfphot = 0
        ELSEIF (kfphot < 0) THEN
          nfphot = int( abs(kfphot)/cdt )
        ELSE
          nfphot = kfphot
        ENDIF

!       2.2  Compute minimum safe value (see Method above)

        sfmin = tiny(1.0d0)
        sfmin = 10.0**(int(log10(sfmin))+1)
        peps  = 1.0e19 * sfmin

!       3.  Set fixed vmrs (Now done in UKCA_MAIN1)

!       4.  Clear the species arrays

        f      = 0.0
        fdot   = 0.0
        ej     = 0.0
        linfam = .false.
        linfam = .false.

        y    = 0.0
        ydot = 0.0
        prod = 0.0
        slos = 0.0
        dpd  = 0.0
        dpw  = 0.0
        emr  = 0.0

!       5.   Clear the rates and index arrays.
!            ----- --- ----- --- ----- -------


        rk   = 0.0
        prk  = 0.0
        nspi = 0

        DO js = 1, jpspec
          ngrp(js,1)          = 0
          ngrp(js,2)          = 0
          ngrp(js,3)          = 0
          nprdx2(1,js)        = 0
          nprdx2(2,js)        = 0
          nprdx1(js)          = 0
          ngrp(js+jpspec,1)   = 0
          ngrp(js+jpspec,2)   = 0
          ngrp(js+jpspec,3)   = 0
          nprdx2(1,js+jpspec) = 0
          nprdx2(2,js+jpspec) = 0
          nprdx1(js+jpspec)   = 0
          nlall(js)           = 0
          nlstst(js)          = 0
          nlf(js)             = 0
          nlmajmin(js)        = 0
          nldepd(js)          = 0
          nldepw(js)          = 0
          nlemit(js)          = 0
          nldepx(js)          = 0
        ENDDO

        DO js = 1, 2*jpspec
          DO jx = 1, jpnpx3
            nprdx3(1,jx,js) = 0
            nprdx3(2,jx,js) = 0
            nprdx3(3,jx,js) = 0
          ENDDO
        ENDDO

        nbrkx = 0
        ntrkx = 0
        nprkx = 0
        nhrkx = 0

        njacx3(1,:,:) = 0
        njacx3(2,:,:) = 0
        njacx3(3,:,:) = 0

        njcgrp(:,1) = 0
        njcgrp(:,2) = 0
        njcgrp(:,3) = 0
        njacx2(1,:) = 0
        njacx2(2,:) = 0
        njacx1(:)   = 0
        nltrf(:)    = 0
        nltr3(:)    = 0

        DO jc = 1, jpctr
          nmpjac(jc) = 0
          DO jp = 1, jppjac
            npjac1(jp,jc) = 0
          ENDDO
        ENDDO

! Initialise the character arrays
        spb(:,:) = nullx
        spt(:,:) = nullx
        spj(:,:) = nullx
        sph(:,:) = nullx

! Initialise the fractional product arrays
        frpb(:)  = 0.0
        frpt(:)  = 0.0
        frpj(:)  = 0.0
        frph(:)  = 0.0
        frpx(:)  = 0.0
        nfrpx    = 0

        ntabfp(1:jpfrpx,1) = 0
        ntabfp(1:jpfrpx,2) = 0
        ntabfp(1:jpfrpx,3) = 0
        nmzjac = 0
        nmsjac = 0
        nzjac1 = 0
        nsjac1 = 0
        ntabpd = 0
        ztabpd = 0.0
        npdfr  = 0

!       6.  Read chemistry data
!           ---- --------- ----

! DEPENDS ON: asad_inrats
        CALL asad_inrats

! Check that deposition and emission is not on for constant species
        DO js = 1, jpspec
          IF ( ldepd(js) .and. ctype(js)(1:1)  ==  'C' ) THEN
            cmessage='Dry deposition turned on for constant species'

            CALL EREPORT('ASAD_CINIT',js,cmessage)
          ENDIF
          IF ( ldepw(js) .and. ctype(js)(1:1)  ==  'C' ) THEN
            cmessage='Wet deposition turned on for constant species'

            CALL EREPORT('ASAD_CINIT',js,cmessage)
          ENDIF
          IF ( lemit(js) .and. ctype(js)(1:1)  ==  'C' ) THEN
            cmessage='Emission turned on for constant species'

            CALL EREPORT('ASAD_CINIT',js,cmessage)
          ENDIF
        ENDDO


        IF ( method == 3) THEN   ! For Newton-Raphson solver only
! DEPENDS ON: asad_setsteady
          CALL asad_setsteady    ! Initialize steady-state species
        ENDIF

        IF ( method >= 10 ) THEN
          DO j = 1, jpspec
            IF ( ctype(j)  ==  jpfm .or. ctype(j)  ==  jpif ) THEN
              WRITE(6,*) '*** ASAD ERROR: You cannot use families ',   &
            ' with one of the stiff integrators. If method  >=  10 ',  &
            ' you cannot have species specified as ',jpfm,' or ',jpif
              cmessage = 'ASAD ABORTED'

              CALL EREPORT('ASAD_CINIT',j,cmessage)
            ENDIF
          ENDDO
        ENDIF

!       7.  Set up the index arrays.
!           --- -- --- ----- -------

! DEPENDS ON: asad_inix
        CALL asad_inix
! DEPENDS ON: asad_inijac
        CALL asad_inijac

!       8.  Initialise photolysis and heterogeneous chemistry
!           ---------- ---------- --- ------------- ---------

! These are dummy routines at vn7.0 of the UM
!! DEPENDS ON: asad_inphot
!        CALL asad_inphot
!! DEPENDS ON: asad_inhet
!        CALL asad_inhet

!       9.  Read deposition and emission data
!           ---- ---------- --- -------- ----

! DEPENDS ON: ukca_inwdep
        IF ( ndepw /= 0 ) CALL UKCA_INWDEP
! DEPENDS ON: ukca_inddep
        IF ( ndepd /= 0 ) CALL UKCA_INDDEP
! DEPENDS ON: asad_inemit
        IF ( nemit /= 0 ) CALL asad_inemit

        IF (lhook) CALL dr_hook('ASAD_CINIT',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_CINIT
