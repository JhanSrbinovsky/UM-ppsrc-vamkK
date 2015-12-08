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
! Purpose: Main chemistry driver routine
!
!     Chemistry driver routine. If necessary, model concentrations are
!     converted from vmr to number density. If the chemistry is to be
!     "process-split" from the transport, the chemistry tendencies are
!     integrated by the chosen method and an average chemistry tendency
!     over the model timestep is returned. If the chemistry is not to
!     be integrated, instantaneous chemistry tendencies are returned.
!     The tracer tendencies returned to the calling routine are the
!     final values held in the chemistry.
!
!     Note: It is important for conservation that the concentrations
!     passed to this routine are non-negative.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from UKCA_CHEMISTRY_CTL
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!     Arguments:
!        cdot        - Tracer tendencies due to chemistry.
!        ftr         - Tracer concentrations.
!        pp          - Pressure (Nm-2).
!        pt          - Temperature (K).
!        pq          - Water vapor field (vmr).
!        nlev        - Model level
!        dryrt       - Dry deposition rates (s-1)
!        wetrt       - Wet deposition rates (s-1)
!        n_points    - No. of points calculations be done.
!
!     Method
!     ------
!     Since the heterogeneous rates may depend on species
!     concentrations, the call to hetero is made inside the
!     loop over chemistry sub-steps.
!
!     Photolysis rates may be computed at frequency set by the
!     variable nfphot and so the call to photol is also inside
!     the chemical sub-step loop.
!
!     Externals
!     ---------
!     asad_bimol  - Calculates bimolecular rate coefficients.
!     asad_trimol - Calculates trimolecular rate coefficients.
!     ukca_photol - Calculates photolysis rate coefficients.
!     ukca_drydep - Calculates dry deposition rates.
!     ukca_wetdep - Calculates wet deposition rates.
!     asad_emissn - Calculates emission rates.
!     asad_totnud - Calculates total number densities.
!     asad_impact - IMPACT time integration scheme.
!     asad_hetero - Calculates heterogeneous rate coefficients.
!     asad_posthet- Performs housekeeping after heterogeneous chemistry.
!     asad_ftoy   - Partitions families.
!     asad_diffun - Calculates chemistry tendencies.
!
!     Local variables
!     ---------------
!     ifam       Family index of in/out species.
!     itr        Tracer index of in/out species.
!     iodd       Number of odd atoms in in/out species.
!     gfirst     .true. when the species need to be
!                initialised on the first chemical step.
!     gphot      .true. if photol needs to be called.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_CDRIVE(cdot, ftr, pp, pt, pq, cld_f, cld_l,     &
                               nlev, dryrt, wetrt, rc_het, prt,         &
                               n_points,stratflag)
        USE ASAD_MOD
        USE UKCA_HETERO_MOD
        USE ukca_option_mod, ONLY: L_ukca_achem, L_ukca_trophet,        &
                                   L_ukca_het_psc, jpctr, jpspec,       &
                                   jpdd, jppj 
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
        USE ereport_mod, ONLY : ereport
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

! Subroutine interface
        INTEGER, INTENT(IN) :: n_points              ! No of points
        INTEGER, INTENT(IN) :: nlev                  ! Model level

        REAL, INTENT(IN) :: prt(n_points,jppj)       ! Photolysis rates
        REAL, INTENT(IN) :: dryrt(n_points,jpdd)     ! Dry dep rates
        REAL, INTENT(IN) :: wetrt(n_points,model_levels,jpdw) ! Wet dep rates
        REAL, INTENT(IN) :: rc_het(n_points,2)       ! Hetero. Chemistry rates
        REAL, INTENT(IN) :: pp(n_points)             ! Pressure
        REAL, INTENT(IN) :: pt(n_points)             ! Temperature
        REAL, INTENT(IN) :: pq(n_points)             ! Water vapour
        REAL, INTENT(IN) :: cld_f(n_points)          ! Cloud fraction
        REAL, INTENT(IN) :: cld_l(n_points)          ! Cloud liquid water (kg/kg)
        LOGICAL, INTENT(IN) :: stratflag(n_points)   ! Strat indicator

        REAL, INTENT(INOUT) :: ftr(n_points,jpctr)   ! Tracer concs

        REAL, INTENT(OUT)   :: cdot(n_points,jpctr)  ! Tracer tendencies

!       Local variables

        INTEGER :: errcode               ! Variable passed to ereport

        INTEGER :: jtr                                ! Loop variable
        INTEGER :: jl                                 ! Loop variable
        INTEGER :: js                                 ! Loop variable
        INTEGER :: nl
        INTEGER :: ifam
        INTEGER :: itr
        INTEGER :: iodd

        LOGICAL :: gfirst
        LOGICAL :: gphot

        CHARACTER(len=72) :: cmessage          ! Error message

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


!       1.  Initialise variables and arrays

!       1.1   Clear tendencies to avoid contributions from levels
!             on which no chemistry is performed

        IF (lhook) CALL dr_hook('ASAD_CDRIVE',zhook_in,zhook_handle)
        DO jtr = 1, jpctr
          DO jl = 1, n_points
            cdot(jl,jtr) = 0.0
          ENDDO
        ENDDO

!       1.2  Copy pressure and temperature to common

        DO jl = 1, n_points
          p(jl) = pp(jl)
          t(jl) = pt(jl)
        ENDDO

!       1.2.1 Copy water vapor to common

        DO jl = 1, n_points
          wp(jl) = pq(jl)
        ENDDO

!       2.  Calculate total number densities

! DEPENDS ON: asad_totnud
        CALL asad_totnud(n_points)

!       3.  Read model tracer concentrations into working array,
!           and if necessary, convert vmr to number densities

        IF ( lvmr ) THEN
          do jtr = 1, jpctr
            do jl  = 1, n_points
              ftr(jl,jtr) = ftr(jl,jtr) * tnd(jl)
              f(jl,jtr)   = ftr(jl,jtr)
            ENDDO
          ENDDO
        ELSE
          do jtr = 1, jpctr
            do jl  = 1, n_points
              f(jl,jtr)   = ftr(jl,jtr)
            ENDDO
          ENDDO
        ENDIF

!       4.  Calculate reaction rate coefficients
!           --------- -------- ---- ------------

! DEPENDS ON: asad_bimol
        CALL asad_bimol (n_points)
! DEPENDS ON: asad_trimol
        CALL asad_trimol(n_points)

! Calculate aqueous-phase SO2 oxdn. and tropospheric heterogenous rates
        IF (L_ukca_achem .OR. L_ukca_trophet)                           &
          THEN
! DEPENDS ON: asad_hetero
          CALL ASAD_HETERO(n_points, cld_f, cld_l, rc_het)
        END IF

!       5.  Calculate deposition and emission rates
!           --------- ---------- --- -------- -----

! DEPENDS ON: ukca_wetdep
        IF ( ndepw /= 0 ) CALL UKCA_WETDEP(nlev, wetrt, n_points)
! DEPENDS ON: ukca_drydep
        IF ( ndepd /= 0 ) CALL UKCA_DRYDEP(nlev, dryrt, n_points)
! DEPENDS ON: asad_emissn
        IF ( nemit /= 0 ) CALL asad_emissn

!       6.  Integrate chemistry by chosen method. Otherwise,
!           simply calculate tendencies due to chemistry
!           ------ --------- ---------- --- -- ---------

        gphot = .true.
        IF ( method /= 0 ) THEN

          DO jsubs = 1, ncsteps

            gfirst = jsubs  ==  1
            gphot  = gfirst
            IF ( nfphot /= 0 .AND. .NOT.gfirst )                       &
                          gphot = mod(jsubs-1,nfphot) == 0

!           ---------------------------------------------------
!           NON-STIFF integrators take values in the range 1-9.
!           ---------------------------------------------------

!           6.1  IMPACT integration: first compute heterogeneous
!                and photolysis rates, species and tendencies.
!                ===============================================

            nl = n_points
            IF ( method == 1 ) THEN
! DEPENDS ON: ukca_photol
              IF ( gphot ) CALL UKCA_PHOTOL(prt,n_points)
              IF (L_ukca_het_psc)  CALL ukca_hetero(n_points,stratflag)
! DEPENDS ON: asad_ftoy
              CALL asad_ftoy( gfirst, nitfg, n_points )
! DEPENDS ON: asad_diffun
              CALL asad_diffun( nl )
! DEPENDS ON: asad_jac
              CALL asad_jac( n_points )
! DEPENDS ON: asad_impact
              CALL asad_impact( n_points )

!             6.2.  Quasi-steady state scheme.
!             ================================

            ELSEIF ( method == 2 ) THEN
              cmessage='QSSA not in UM6.5 build'
              errcode=1
              CALL EREPORT('ASAD_CDRIVE',errcode,cmessage)

!              6.3   Sparse Newton-Raphson solver
!              ==================================

            ELSEIF ( method == 3 ) THEN
! DEPENDS ON: ukca_photol
              IF ( gphot ) CALL UKCA_PHOTOL(prt,n_points)
              IF (L_ukca_het_psc)  CALL UKCA_HETERO(n_points,stratflag)
! DEPENDS ON: asad_ftoy
              CALL ASAD_FTOY( gfirst, nitfg, n_points )
               
! DEPENDS ON: asad_spmjpdriv
              CALL ASAD_SPMJPDRIV(nlev,n_points)

!              6.5   Backward Euler solver
!              ===========================

            ELSEIF ( method == 5 ) then
! DEPENDS ON: ukca_photol
              IF ( gphot ) CALL UKCA_PHOTOL(prt,n_points)
              IF (L_ukca_het_psc)  CALL UKCA_HETERO(n_points,stratflag)
! DEPENDS ON: asad_ftoy
              CALL ASAD_FTOY( gfirst, nitfg, n_points )
! DEPENDS ON: asad_bedriv
              CALL ASAD_BEDRIV(nlev,mype,n_points)

!           -------------------------------------------------
!           STIFF integrators take values in the range 10-19.
!           -------------------------------------------------

!           6.10  NAG BDF stiff integrator.

            ELSEIF ( method == 10 ) THEN
              cmessage='NAG BDF not in UM6.5 build'
              errcode=1
              CALL EREPORT('ASAD_CDRIVE',errcode,cmessage)

!             6.11  SVODE ODE stiff integrator from NETLIB.

            ELSEIF ( method == 11 ) THEN
              cmessage='SVODE not in UM6.5 build'
              errcode=1
              CALL EREPORT('ASAD_CDRIVE',errcode,cmessage)

            ENDIF

!           6.12  Do any final work for the heterogeneous
!                 chemistry before the end of the time loop.
!                 =========================================

            IF ( L_ukca_het_psc ) CALL ukca_solidphase(n_points)
! DEPENDS ON: asad_posthet
!            IF ( L_ukca_het_psc ) CALL asad_posthet

          ENDDO        ! End of looping over chemical timesteps

        ELSE           ! Method is equal to zero

!       6.99  Not integrating: just compute tendencies.

          method = 0
! DEPENDS ON: ukca_photol
          IF ( gphot ) CALL ukca_photol(prt,n_points)
          IF ( L_ukca_het_psc )  CALL ukca_hetero(n_points,stratflag)
! DEPENDS ON: asad_ftoy
          CALL asad_ftoy( .true., nit0, n_points )
! DEPENDS ON: asad_diffun
          CALL asad_diffun( nl )
          IF ( L_ukca_het_psc ) CALL ukca_solidphase( n_points )
! DEPENDS ON: asad_posthet
          IF ( L_ukca_het_psc ) CALL asad_posthet

        ENDIF        ! End of IF statement for method


!       7.  Determine concentrations and tendencies to be returned to
!           the model depending on whether or not the chemistry has
!           been integrated  -- -- ------- -- --- --- --------- ---
!           ---- ----------

!       7.1  Obtain model family concentrations by subtracting
!            concentrations of in/out species from ASAD families

        DO js = 1, jpspec
          IF ( ctype(js) == jpif ) THEN
            ifam = moffam(js)
            itr  = madvtr(js)
            iodd = nodd(js)
            DO jl = 1, n_points
              IF ( linfam(jl,itr) ) f(jl,ifam) =                       &
                                    f(jl,ifam) - iodd*f(jl,itr)
            ENDDO
          ENDIF
        ENDDO

!       7.2  Returned values of concentration and chemical tendency

        DO jtr = 1, jpctr
          IF ( method /= 0 ) THEN
            DO jl = 1, n_points
              cdot(jl,jtr) = ( f(jl,jtr)-ftr(jl,jtr)) / (cdt*ncsteps)
              ftr(jl,jtr)  = f(jl,jtr)
            ENDDO
          ELSE
            DO jl = 1, n_points
              cdot(jl,jtr) = fdot(jl,jtr)
              ftr(jl,jtr)  = f(jl,jtr)
            ENDDO
          ENDIF
        ENDDO


!       8.  If necessary, convert from number densities back to vmr
!           -- ---------- ------- ---- ------ --------- ---- -- ---

        IF ( lvmr ) THEN
          DO jtr = 1, jpctr
            DO jl = 1, n_points
              ftr(jl,jtr)  = ftr(jl,jtr)  / tnd(jl)
              cdot(jl,jtr) = cdot(jl,jtr) / tnd(jl)
            ENDDO
          ENDDO
        ENDIF

        IF (lhook) CALL dr_hook('ASAD_CDRIVE',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_CDRIVE
