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
! Purpose: Sets the species concentrations prior to calculating the
!          chemistry tendencies. This includes partitioning families.
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE and ASAD_IMPACT
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!     Interface
!     ---------
!     Must always be called before diffun.
!
!     Arguments:
!        ofirst     Should be .true. if this is first call to ftoy
!                   during the current model/dynamical timestep.
!        iter       Maximum number of iterations.
!
!     Method
!     ------
!    After initialisation, the family member concentrations and
!    steady state species are iterated until they converge to a
!    solution. Firstly, the ratios for family members and steady
!    state concentrations for the non-family steady state species
!    are computed. Then, the family member concentrations are
!    evaluated using the ratios. A convergence test is performed
!    for all species at all spatial points. If convergence is not
!    achieved for any species at any point, the iterations continue
!    until either convergence is achieved or the maximum number of
!    iterations is reached.
!
!    In/out species are determined to be in or out of the family
!    on the first call to ftoy during a model timestep according
!    to their lifetime compared with some threshold value. If the
!    species is found to be 'in' the family, its concentration
!    is added to the family f concentration before the sum
!    is partitioned amongst the members. The in/out test is
!    performed at all spatial points. This species then stays in
!    the family for the rest of the call to cdrive. The species is
!    subtracted from the family at the end of cdrive.
!
!    If the species self reacts, the steady state concentration
!    is found by solving a quadratic. Otherwise, it is simply
!    found by dividing the production rate by the loss rate.
!    The self reacting term, qa, has previously been calculated
!    in fyself. Note qa is +ve.
!
!    The ratio of a minor family member to the major member is
!    determined by dividing the steady state concentration of the
!    minor member from fystst by the family member concentration:
!                 R1m = Y1/Ym, R2m = Y2/Ym ....
!    The ratio of the major member to the family as a whole is
!    determined from the minor member ratios weighted by the number
!    of odd atoms in each species:
!                 Rmf = 1 / ( Nm + N1*R1m +N2*R2m + .....)
!    The concentrations are then found from:
!                 Ym = Rmf*Z, Y1 = R1m*Ym, Y2 = R2m*Ym, .....
!    where Z is the remaining family concentration after the
!    appropriate in/out species concentrations have been subtracted.
!
!    Although ASAD supports the user of deposition and emissions,
!    we do not include these terms in the calculation of the ratios
!    even though they represent a loss or production from or in to
!    the family. This is because, family members are assumed to be
!    in steady state, generally taken to be a 'chemical' steady
!    state, one that is controlled by the concentrations of other
!    species. If for instance, a strong emission source was present,
!    for NO, this would give a high bogus value for the NO ratio,
!    which would be wrong because by assuming steady state we assume
!    that the extra NO has already reacted to reach a new equilibrium
!    with its family members.
!
!    Externals
!    ---------
!    fyinit         Sets initial species concentrations.
!    fyself         Calculates self-reacting terms.
!    fyfixr         Calculates family member concentrations
!                   using ratios computed on a previous call.
!    fystst         Calculates concentration of a species
!                   in steady state.
!    prls           Calculates production & loss terms.
!
!    Local variables
!    ---------------
!    ifam           Index of family to which species belongs.
!    imaj           Index of major member of family to which
!                   species belongs.
!    itr            Index of model tracer to which species
!                   corresponds.
!    gconv          .true. if convergence has been achieved.
!    gonce          .true. for only the first call to this routine.
!    zthresh        Threshold value of loss rate to determine
!                   whether in/out species are in or out of family.
!    zy             Value of y on previous iteration.
!    ilstmin        List of species which are either steady state,
!                   'SS', family members ('FM' and 'FT') but EXCLUDING
!                   major family species.
!    istmin         No. of entries in ilstmin.
!    ilft           List of species of type 'FT'.
!    ift            No. of entries in ilft.
!    zb, zc, zd     Variables used in the quadratic. Note that
!                   zb is really '-b' in the quadratic.
!    zb             The loss rate less self reacting terms.
!    zc             The production rate.
!    zd             "b^2 - 4*a*c"
!
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_FTOY(ofirst,iter,n_points)

        USE ASAD_MOD,            ONLY: f, y, prod, slos, ratio, qa,     &
                                       peps, cdt, pmintnd, nstst,       &
                                       jpfm, jpif, jpna, moffam,        &
                                       majors, ilstmin, ilft, nodd,     &
                                       nlmajmin, madvtr, linfam,        &
                                       nlstst, ctype, ftol
        USE ukca_option_mod,     ONLY: jpspec
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

        INTEGER, INTENT(IN)    :: n_points  ! No of spatial points
        INTEGER, INTENT(INOUT) :: iter      ! Max no of iterations

        LOGICAL, INTENT(IN) :: ofirst

!       Local variables

        INTEGER       :: j           ! Loop variable
        INTEGER       :: jit         ! Loop variable
        INTEGER       :: jl          ! Loop variable
        INTEGER       :: js          ! Loop variable
        INTEGER       :: ifam
        INTEGER       :: imaj
        INTEGER       :: nl
        INTEGER       :: istart
        INTEGER       :: iend
        INTEGER       :: iodd
        INTEGER       :: itr
        INTEGER       :: icode       ! Error code

        INTEGER, SAVE :: istmin
        INTEGER, SAVE :: ift
!        INTEGER       :: ilstmin(jpspec)      ! Now in ASAD_MOD
!        INTEGER       :: ilft(jpspec)         ! Now in ASAD_MOD

        CHARACTER (LEN=72) :: cmessage     ! Error message

        REAL          :: zthresh
        REAL          :: sl
        REAL          :: zy(theta_field_size,jpspec)
        REAL          :: zb(theta_field_size)
        REAL          :: zc(theta_field_size)
        REAL          :: zd(theta_field_size)

        LOGICAL       :: gconv
        LOGICAL, SAVE :: gonce  = .true.
        LOGICAL, SAVE :: gdepem = .false.

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


!       1.  Initialise.
!           -----------

        IF (lhook) CALL dr_hook('ASAD_FTOY',zhook_in,zhook_handle)
        IF ( gonce ) THEN
          gonce  = .false.
          istmin = 0
          ift    = 0
          DO j = 1, jpspec
            ilstmin(j) = 0
            ilft(j) = 0
          ENDDO

          DO j = 1, nstst
            js = nlstst(j)
            IF ( ctype(js) == jpfm ) THEN
              ifam = moffam(js)
              imaj = majors(ifam)
              IF ( imaj /= js ) THEN
                istmin          = istmin + 1
                ilstmin(istmin) = js
              ENDIF
            ELSEIF ( ctype(js) == jpif ) THEN
              istmin          = istmin + 1
              ilstmin(istmin) = js
              ift             = ift + 1
              ilft(ift)       = js
            ELSEIF ( ctype(js) == jpna ) THEN
              istmin          = istmin + 1
              ilstmin(istmin) = js
            ELSE
              WRITE (6,'(A)') '**** ASAD ERROR in ASAD_FTOY!! '
              WRITE (6,'(A)') 'ASAD_FTOY found an unexpected species type'
              WRITE (6,'(A)') 'in the species list nlstst ', nlstst
              cmessage='Found unexpected species type'
              CALL EREPORT('ASAD_FTOY',js,cmessage)
            ENDIF
          ENDDO
        ENDIF   ! End of IF (gonce) statement

!       1.1 Set concentrations/initialise species

! DEPENDS ON: asad_fyinit
        CALL ASAD_FYINIT(ofirst,n_points)

!       1.2 If there are no family or steady species then exit

        IF ( nstst == 0) THEN 
          IF (lhook) CALL dr_hook('ASAD_FTOY',zhook_out,zhook_handle)
          RETURN
        END IF

!       1.3 Initialise local variables and do sanity checks.

        nl = n_points

        DO j = 1, nstst
          js = nlstst(j)
          DO jl = 1, n_points
            zy(jl,js) = 0.0
          ENDDO
        ENDDO

        IF (ofirst .AND. iter < 5) THEN
          iter = 5
          icode = -1
          cmessage = 'iter too low on first call, resetting to 5'
          CALL EREPORT('ASAD_FTOY',icode,cmessage)
        END IF

        zthresh = 2.0 / cdt

!       2.  Calculate self-reacting terms
!           --------- ------------- -----

! DEPENDS ON: asad_fyself
        IF ( ofirst ) CALL ASAD_FYSELF(n_points)

!       3.  Calculate family members using previous ratios
!           --------- ------ ------- ----- -------- ------

        IF (iter == 0) THEN
! DEPENDS ON: asad_fyfixr
          CALL ASAD_FYFIXR(n_points)
          IF (lhook) CALL dr_hook('ASAD_FTOY',zhook_out,zhook_handle)
          RETURN
        ENDIF

!       4.  Iterate family members (including in/out species)
!           ------- ------ ------- --------------------------
!           and non-model steady-state species
!           --- --------- ------------ -------

        DO jit = 1, iter

!         check to see if a species has gone zero with a nonzero
!         production.

          IF ( jit > 1 ) THEN
            DO j = 1, nstst
              js = nlstst(j)
              DO jl = 1, n_points
                IF (abs(y(jl,js)) < peps .AND.                         &
                    prod(jl,js) > peps ) THEN
                  y(jl,js) = peps
                ENDIF
              ENDDO
            ENDDO
          ENDIF

!         4.1 Calculate production and loss terms.
!             Note that the effect of deposition (wet & dry) and
!             emissions are not included in the ratios calculation.
!             Regardless of whether the user has them turned on or not.
!             See method above.

! DEPENDS ON: asad_prls
          CALL ASAD_PRLS( nl, istmin, ilstmin, gdepem )

!         4.2 Initialise ratios

          DO j = 1, nstst
            js = nlstst(j)
            DO jl = 1, n_points
             ratio(jl,js) = 0.0
            ENDDO
          ENDDO

!         4.3  Compute species values for species in steady state
!              (minor members of family, 'SS' and 'FT' species)
!              N.B. because of the way asad computes slos, the y
!              value can go slightly negative during the quadratic
!              due to loss of precision. We forcibly fix it here.
!              Also note that we cannot permit y=0.0 during the ftoy
!              iteration because of the need to get 'sl'.

          DO j = 1, istmin
            js = ilstmin(j)
            DO jl = 1, n_points
              IF ( y(jl,js) < peps ) THEN
                sl = 0.0
              ELSE
                sl = slos(jl,js) / y(jl,js)
              ENDIF
              IF ( qa(jl,js) > peps ) THEN
                zb(jl) = sl - qa(jl,js) * y(jl,js)
                zc(jl) = prod(jl,js)
                zd(jl) = zb(jl)*zb(jl) + 4.0*qa(jl,js)*zc(jl)
                IF ( zd(jl) > 0.0 ) THEN
                  y(jl,js) = (zb(jl) - sqrt(zd(jl)))/(-2.0*qa(jl,js))
                  if ( y(jl,js)  <   0.0 ) y(jl,js) = 10.0*peps
                ELSE
                  y(jl,js) = zb(jl) / ( -2.0 * qa(jl,js) )
                ENDIF
              ELSE IF (qa(jl,js) <= peps .AND. sl > peps ) THEN
                y(jl,js) = prod(jl,js) / sl
              ELSE
                y(jl,js) = 0.0
              END IF
            ENDDO
          ENDDO

!         4.4  Now compute ratios for minor species members.
!              ** could possibly take out 1/y(imaj) and convert to '*'

          istart = nlmajmin(3)
          iend   = nlmajmin(4)
          DO j = nlmajmin(3), nlmajmin(4)
            js = nlmajmin(j)
            iodd = nodd(js)
            ifam = moffam(js)
            imaj = 0
            IF ( ifam /= 0 ) imaj = majors(ifam)
            DO jl = 1, n_points
              IF ( y(jl,imaj) > peps ) THEN
                ratio(jl,js)   = y(jl,js) / y(jl,imaj)
                ratio(jl,imaj) = ratio(jl,imaj) +                      &
                                 iodd * ratio(jl,js)
              ELSE
                ratio(jl,js) = 0.0
              END IF
            ENDDO
          ENDDO

!         4.5  Now compute ratios from species of type 'FT' if
!              applicable. If this is the first call, we also set
!              whether the species is in or out of the family.

          DO j = 1, ift
            js = ilft(j)
            iodd = nodd(js)
            ifam = moffam(js)
            itr  = madvtr(js)
            imaj = majors(ifam)
            IF ( ofirst .and. jit == 1 ) THEN
              DO jl = 1, n_points
                IF ( y(jl,js) > peps ) then
                  sl = slos(jl,js) / y(jl,js)
                ELSE
                  sl = 0.0
                ENDIF
                linfam(jl,itr) = sl  >   zthresh
                IF ( linfam(jl,itr) )                                  &
                     f(jl,ifam) = f(jl,ifam) + iodd*f(jl,itr)
              ENDDO
            END IF

            DO jl = 1, n_points
              if ( linfam(jl,itr) ) then
                if ( y(jl,imaj)  >   peps ) then
                  ratio(jl,js)   = y(jl,js) / y(jl,imaj)
                  ratio(jl,imaj) = ratio(jl,imaj) +                    &
                                   iodd * ratio(jl,js)
                else
                  ratio(jl,js) = 0.0
                endif
              else
                y(jl,js) = f(jl,itr)
              endif
            ENDDO
          ENDDO

!         4.6  Finally compute ratio of major species to family and
!              hence set the minor family members.

          istart = nlmajmin(1)
          iend   = nlmajmin(2)
          DO j = istart, iend
            js = nlmajmin(j)
            iodd = nodd(js)
            ifam = moffam(js)
            DO jl = 1, n_points
              ratio(jl,js) = 1.0 / (iodd + ratio(jl,js))
              y(jl,js)     = f(jl,ifam) * ratio(jl,js)
            ENDDO
          ENDDO

          istart = nlmajmin(3)
          iend   = nlmajmin(4)
          DO j = istart, iend
            js = nlmajmin(j)
            ifam = moffam(js)
            imaj = majors(ifam)
            IF ( ctype(js) /= jpif ) THEN
              DO jl = 1, n_points
                y(jl,js) = y(jl,imaj) * ratio(jl,js)
              ENDDO
            else
              itr = madvtr(js)
              DO jl = 1, n_points
                if ( linfam(jl,itr) ) y(jl,js) =                       &
                     y(jl,imaj) * ratio(jl,js)
              ENDDO
            ENDIF
          ENDDO

!         4.7 Test for convergence.

          gconv = .true.
          DO j = 1,nstst
            js = nlstst(j)
              DO jl = 1, n_points
                IF ( ABS(y(jl,js)-zy(jl,js)) >  ftol*y(jl,js)          &
                .AND. y(jl,js) >  pmintnd(jl) ) gconv=.false.
              ENDDO
            IF ( .NOT. gconv ) EXIT
          ENDDO
          IF (gconv) THEN 
            IF (lhook) CALL dr_hook('ASAD_FTOY',zhook_out,zhook_handle)
            RETURN
          END IF

          DO j = 1, nstst
            js = nlstst(j)
            DO jl = 1, n_points
              zy(jl,js) = y(jl,js)
            ENDDO
          ENDDO

        ENDDO    ! End of iterations


!       9. Convergence achieved or max. iterations reached.

        IF (lhook) CALL dr_hook('ASAD_FTOY',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_FTOY
