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
! Purpose: Controlling routine for the IMPACT integration scheme
!
!  Part of the UKCA model, a community model supported by
!  The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
!          Called from ASAD_CDRIVE
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!     Interface
!     ---------
!     Called from chemistry driver routine cdrive.
!
!     Arguments:
!          lphot - Logical array. One entry per gridpoint. If set .T.
!                  then photolysis is occurring and the scheme will
!                  adjust the accuracy of the Jacobian. See method.
!
!     All other variables are obtained via COMMON.
!
!     Method
!     ------
!     The IMPACT scheme used the Euler backward implicit timescheme
!     as its starting point. This is a first order conservative
!     scheme. A predictor-corrector approach is used to solve this
!     by employing a Newton-Raphson type iteration. However, rather
!     than compute and invert a near full Jacobian matrix, this scheme
!     only considers the diagonal terms only. For mildly stiff chemistry
!     this should work fine. However, if the scheme is not converging
!     rapidly enough during daylight, this is a symptom that the
!     off-diagonal production terms, which are usually ignored, are
!     important. In this case, set the switch ljacx .T. and make sure
!     that the lphot array is set .true. where photolysis is on.
!     Since, only approximate forms of the Jacobian are considered, if
!     switched on, the method will add in the terms in the Jacobian
!     due to photolysis. It does this in an approximate way to avoid
!     an matrix inversion but this greatly improves the performance of
!     the scheme under stiffer conditions.
!
!     The IMPACT scheme is fully described in the forthcoming paper;
!     Carver and Stott, 1997, 'IMPACT: an implicit time integration
!     scheme for chemical species and families', to be submitted to
!     Annales Geophysicae. Preprints available.
!
!     Externals
!     ---------
!     ftoy         Converts f concentrations to y (partitions
!                  families).
!     diffun       Computes tendencies due to chemistry and
!                  main diagonal Jacobian elements.
!
!     Local variables
!     ---------------
!     gconv         .true. if convergence in Newton-Rhapson has been
!                   achieved.
!     ifi           Number of ftoy iterations.
!     zf            Values of f at start of chemistry step.
!     zprf          Value of f on previous Newton-Rhapson iteration.
!     binv          Inverse of the approximate (main diagonal only)
!                   of the matrix ( I - dt.J) where J is Jacobian.
!     dely          Represents the change in variable, y, produced by
!                   multiplying binv to the function we are trying to
!                   find the root of.
!
! Code description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
! ---------------------------------------------------------------------
!
        SUBROUTINE ASAD_IMPACT(n_points)

        USE ASAD_MOD
        USE ukca_option_mod, ONLY: jpctr
        USE parkind1, ONLY: jprb, jpim
        USE yomhook, ONLY: lhook, dr_hook
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

        INTEGER, INTENT(IN ) :: n_points   ! No of spatial points

!       Local variables

        INTEGER :: lphot(theta_field_size)
        INTEGER :: inl
        INTEGER :: nl
        INTEGER :: j                        ! Loop variable
        INTEGER :: jit                      ! Loop variable
        INTEGER :: jl                       ! Loop variable
        INTEGER :: jr                       ! Loop variable
        INTEGER :: jtr                      ! Loop variable
        INTEGER :: jt1                      ! Loop variable
        INTEGER :: isp                      ! Index
        INTEGER :: ifi
        INTEGER :: ntr1
        INTEGER :: ntr2
        INTEGER :: ipos0
        INTEGER :: ipos
        INTEGER :: ir
        INTEGER :: irk
        INTEGER :: njr
        INTEGER :: ireac
        INTEGER :: iprod

        REAL :: zf(theta_field_size,jpctr)
        REAL :: zprf(theta_field_size,jpctr)
        REAL :: binv(theta_field_size,jpctr)
        REAL :: dely(theta_field_size,jpctr)
        REAL :: dd(theta_field_size)
        REAL :: corrn(theta_field_size)

        LOGICAL :: gconv
        LOGICAL :: gfam
        LOGICAL, SAVE :: gfirst = .true.

        INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
        INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
        REAL(KIND=jprb)               :: zhook_handle


!       1.  Predictor step (first guess)
!           --------- ---- ------ ------

!       Crude way of determining if photolysis is on or not.
!       Needs improving.

        IF (lhook) CALL dr_hook('ASAD_IMPACT',zhook_in,zhook_handle)
        inl = 0
        DO jl = 1, n_points
          IF ( rk(jl,nprkx(1)) > peps ) THEN
            inl = inl + 1
            lphot(inl) = jl
          endif
        ENDDO

        IF ( gfirst ) THEN
          gfirst = .false.
! DEPENDS ON: asad_inimpct
          CALL ASAD_INIMPCT
        ENDIF

!       1.1  Do the linearised first guess to give first approx.
!            solution at y(n+1).

        nl = n_points
        DO jtr = 1, jpctr
          isp = majors(jtr)
          DO jl = 1, n_points
            zf(jl,jtr) = f(jl,jtr)

!           1.2  Test Jacobian has not gone positive. Possible since
!                we only use approximate form. If it has, use
!                an explicit step.

            if ( ej(jl,jtr) > 0.0 ) then
              f(jl,jtr) = zf(jl,jtr) + cdt*fdot(jl,jtr)
            else
              f(jl,jtr) = zf(jl,jtr) + ( cdt*fdot(jl,jtr) )            &
                               / ( 1.0 - cdt*ej(jl,jtr) )
            endif
            if( linfam(jl,jtr) ) f(jl,jtr) = y(jl,isp)

          ENDDO
        ENDDO

!       2.  Newton-Rhapson iteration (corrector step).
!           -------------- --------- ---------- ------

        DO jit = 1, nrsteps

!         2.1  Decide if we are recomputing ratios or not.

          IF ( jit < 4 ) THEN
            ifi = 0
          ELSE
            ifi = nitnr
          ENDIF

!         2.2  Work out rates of change and main diagonal of J.

! DEPENDS ON: asad_ftoy
          CALL ASAD_FTOY(.false., ifi, n_points)
! DEPENDS ON: asad_diffun
          CALL ASAD_DIFFUN( nl )
! DEPENDS ON: asad_jac
          CALL ASAD_JAC(n_points)

!         2.3  Do the normal Newton-Raphson iteration with
!              Jacobian approximated to main diagonal.

          DO jtr = 1, jpctr
            isp = majors(jtr)
            DO jl = 1, n_points
              zprf(jl,jtr) = f(jl,jtr)
              binv(jl,jtr) = 1.0 / ( 1.0 - cdt*ej(jl,jtr) )
              dely(jl,jtr) = binv(jl,jtr) *                            &
                     ( (zf(jl,jtr) - f(jl,jtr)) + cdt*fdot(jl,jtr) )
              f(jl,jtr) = f(jl,jtr) + dely(jl,jtr)

              if( linfam(jl,jtr) ) f(jl,jtr) = y(jl,isp)
            ENDDO
          ENDDO

!         2.4  If requested for all points where photolysis is
!              occurring, include the terms in the Jacobian
!              arising from just the photolysis terms. Note
!              that we only do this for points where photolysis
!              is actually turned on.

          IF ( ljacx .AND. inl > 0 ) THEN
            ntr1 = nltrim(0,1)

!           2.4.1  For each tracer that needs correcting (ie.
!                  loop over the rows in the Jacobian matrix).

            DO jt1 = 1, ntr1

!             2.4.2  Compute the contribution from non-zero
!                    elements in the Jacobian due to photolysis.
!                    ie. work our way along the columns. We first
!                    compute the contribution from ordinary tracers
!                    'TR' and then from families.

              ntr2  = nltrim(jt1,2)
              ipos0 = nltrim(jt1,3)
              ipos  = ipos0
              DO jl = 1, n_points
                corrn(jl) = 0.0
              ENDDO

              DO

              ir   = nlpdv(ipos,1)
              irk  = nprkx(ir)
              isp  = nspi(irk,1)
              njr  = nlpdv(ipos,2)
              ireac = madvtr(isp)
              gfam = ireac  ==  0
              DO j = 1, inl
                dd(j) = 0.0
              ENDDO

              IF ( gfam ) THEN
                ireac = moffam(isp)
                DO jr = 1, njr
                  ir  = nlpdv(ipos,1)
                  irk = nprkx(ir)
                  isp = nspi(irk,1)
                  DO j = 1, inl
                    jl    = lphot(j)
                    dd(j) = dd(j) + rk(jl,irk)*y(jl,isp)
                  ENDDO
                  ipos = ipos + 1
                ENDDO
                DO j = 1, inl
                  jl    = lphot(j)
                  dd(j) = dd(j) / zprf(jl,ireac)
                ENDDO

              ELSE

!             2.4.3.  Species of type 'FT' will come here.
!                     Check, at each point, whether it's gone
!                     into the family or not.

                DO jr = 1, njr
                  ir  = nlpdv(ipos,1)
                  irk = nprkx(ir)
                  isp = nspi(irk,1)
                  DO j = 1, inl
                    jl = lphot(j)
                    if ( linfam(jl,ireac) ) then
                      dd(j) = dd(j)+rk(jl,irk)*y(jl,isp)/zprf(jl,ireac)
                    else
                      dd(j) = dd(j)+rk(jl,irk)
                    endif
                  ENDDO
                  ipos = ipos + 1
                ENDDO
              ENDIF

              DO j = 1, inl
                jl = lphot(j)
                corrn(j) = corrn(j) + dd(j)*dely(jl,ireac)
              ENDDO

              IF ( ipos - ipos0 >= ntr2 ) EXIT
              ENDDO

!             2.4.3  Now add correction to the tracer (rows). If the
!                    tracer is of type 'FT', then for each pt, we must
!                    check whether the tracer has been put into the f
!                    or not. If it has, we add the correction to the
!                    family and not to the tracer.

              iprod = nltrim(jt1,1)
              isp   = majors(iprod)
              IF ( ctype(isp) == jpif ) THEN
                DO j = 1, inl
                  jl    = lphot(j)
                  iprod = nltrim(jt1,1)
                  IF ( linfam(jl,iprod) ) iprod = moffam(isp)
                  f(jl,iprod) = f(jl,iprod)+cdt*binv(jl,iprod)*corrn(j)
                ENDDO
              ELSE
                DO j = 1, inl
                  jl = lphot(j)
                  f(jl,iprod) = f(jl,iprod)+cdt*binv(jl,iprod)*corrn(j)
                ENDDO
              ENDIF

            ENDDO
          ENDIF      ! End of IF ( ljacx .AND. inl > 0 ) statement

!         9. Check for convergence
!            ----- --- -----------

          gconv = .true.
          DO jtr = 1,jpctr
            DO jl = 1, n_points
              IF ( ABS(f(jl,jtr)-zprf(jl,jtr)) >  ptol*f(jl,jtr)       &
              .AND. f(jl,jtr) >  pmintnd(jl) ) gconv=.false.
            ENDDO
            IF ( .NOT. gconv ) EXIT
          ENDDO
          IF (gconv) THEN
            IF (lhook) CALL dr_hook('ASAD_IMPACT',zhook_out,zhook_handle)
            RETURN
          END IF
        ENDDO           ! End of NR (jit) loop

!       9.1  Convergence achieved or max. no. of iterations reached

        IF (lhook) CALL dr_hook('ASAD_IMPACT',zhook_out,zhook_handle)
        RETURN
        END SUBROUTINE ASAD_IMPACT
