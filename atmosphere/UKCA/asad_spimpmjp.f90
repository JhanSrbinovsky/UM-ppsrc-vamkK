! *****************************COPYRIGHT*******************************
!
! Copyright (c) 2008, Regents of the University of California
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are
! met:
!
!     * Redistributions of source code must retain the above copyright
!       notice, this list of conditions and the following disclaimer. 
!     * Redistributions in binary form must reproduce the above
!       copyright notice, this list of conditions and the following
!       disclaimer in the documentation and/or other materials provided
!       with the distribution. 
!     * Neither the name of the University of California, Irvine nor the
!       names of its contributors may be used to endorse or promote
!       products derived from this software without specific prior
!       written permission.
!
!       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
!       IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!       PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!       OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!       EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!       PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!       NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!    Driver for fully implicit ODE integrator.
!    Part of the ASAD chemical solver
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!   Called from asad_spmjpdriv
!
!     Purpose.
!     --------
!     To organise the integration of the chemical rate equations using
!     the MJP implicit integrator.
!
!     Interface
!     ---------
!     Called from chemistry driver routine via the mjpdriv driver.
!
!     Method.
!     -------
!     Solves the non-linear system of equations via a Newton-Raphson
!     technique. The full Jacobian 'fj' is constructed in 'fuljac', and
!     the net change in family concentration is calculated in 'linslv'.
!     The first few iterations (currently 7) are controlled to prevent
!     very rapid initial changes leading to divergence. Convergence is
!     determined by checking that the total concentration error 'errxo'
!     is less than tolerance 'raferr' or that the total rate error 'errpl'
!     is less than tolerance 'rafpml'. If convergence is not achieved
!     after 'nrsteps', or if divergence is encountered, the routine resets
!     the family concentrations to their initial values, and exits with a
!     non-zero value for ndxraf.
!
!     Global variables
!     ----------------
!     rafpml - tolerance - set to  1.0E-10 in input file
!     rafmin - limit for first few iterations, 0.1 in input file
!     rafmax - limit for first few iterations, 1.0E+04 in input file
!     raferr - tolerance (again?) - set to 1.0E-06 in input file
!     rafbig - maximum concentration above which divergence is assumed
!     rafeps - small non-zero concentration
!
!     Local variables
!     ---------------
!     ifi           Number of ftoy iterations.
!     zf            Values of f at start of chemistry step.
!     zprf          Value of f on previous Newton-Rhapson iteration.
!     ndxraf        Convergence exit code
!     damp1         Damping factor to apply to the first iteration
!     deltt         Reciprocal time step length  (1./cdt)
!
!  Changes for whole-atmosphere chemistry:
!  a. Increase rafeps to sqrt(peps)
!  b. Disactivate crash if too many negatives occur. They are now
!     fine during the iteration.
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      SUBROUTINE asad_spimpmjp(ndxraf, nlev, n_points)

      USE ASAD_MOD,              ONLY: ptol, peps, cdt,                 &
                                       f, fdot, nrsteps, nitnr, nstst,  &
                                       y, prod, slos, fj, lsvjac, ltrig
      USE ASAD_SPARSE_VARS,      ONLY: spfj, pointer2, spfuljac,        &
                                       splinslv2, spresolv2
      USE ukca_option_mod,       ONLY: jpctr
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE PrintStatus_mod
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
      INTEGER, INTENT(IN) :: n_points
      INTEGER, INTENT(IN) :: nlev
      INTEGER, INTENT(OUT):: ndxraf

! Local variables
      INTEGER, PARAMETER :: maxneg=2000     ! Max No. negatives allowed
      INTEGER :: jtr
      INTEGER :: jit
      INTEGER :: ifi
      INTEGER :: i
      INTEGER :: itr
      INTEGER :: nl
      INTEGER :: jl
      INTEGER :: kr
      INTEGER :: j
      INTEGER :: ip

      REAL :: ztmp
      REAL :: rafpml
      REAL :: rafmin
      REAL :: rafmax
      REAL :: raferr
      REAL :: rafbig
      REAL :: rafeps
      REAL :: deltt
      REAL :: damp1
      REAL :: errxo
      REAL :: errpl

      LOGICAL :: errfl80

      REAL :: zf(theta_field_size,jpctr)
      REAL :: xoo(theta_field_size,jpctr)
      REAL :: fxo(theta_field_size,jpctr)
      REAL :: zsum(jpctr)
      REAL :: tmprc(theta_field_size,jpctr)
      REAL :: bx(jpctr)
!
      INTEGER, PARAMETER :: ltrig_jit=51    ! Set to nrsteps if want LTRIG

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!

      IF (lhook) CALL dr_hook('ASAD_SPIMPMJP',zhook_in,zhook_handle)
      rafpml=1.0e-10
      rafmin=1.0e-01
      rafmax=1.0e+04
      raferr=ptol*10.   !   care...

      rafbig=1.0/SQRT(peps)
      rafeps=SQRT(peps)

      ndxraf = 0
      errxo = 1.
      lsvjac = .FALSE.
      deltt = 1./cdt
      damp1 = 0.5

      nl = n_points
!
!  Save values of f at start of step and make linearised first guess
      zf = f
      WHERE (f<rafeps) f = rafeps

! Call ASAD_STEADY at start of step to initialise deriv properly
! DEPENDS ON: asad_steady
      CALL asad_steady( nl )
      
      CALL spfuljac(n_points)
      DO jtr=1,jpctr
        ip = pointer2(jtr,jtr)
        DO jl=1,n_points
          IF  (spfj(jl,ip) > 0.) THEN
            f(jl,jtr) = zf(jl,jtr) + cdt*fdot(jl,jtr)
          ELSE
            f(jl,jtr) = zf(jl,jtr) + (cdt*fdot(jl,jtr))                 &
                                     /(1.0-cdt*spfj(jl,ip))
          END IF
          IF (f(jl,jtr) < rafeps) f(jl,jtr) = rafeps
        END DO
      END DO

!  Start Loop - ensure mixing ratios positive and non-zero on entry
      DO jit=1,nrsteps
        ifi = 0
        IF(jit == 1) ifi=nitnr
!
! DEPENDS ON: asad_ftoy
        CALL asad_ftoy(.FALSE.,ifi, n_points)
        IF (nstst /= 0) THEN
! DEPENDS ON: asad_steady
          if(ifi == 0) CALL asad_steady( nl )
        END IF
!
        IF (ltrig .AND. printstatus >= prstatus_oper) THEN
          DO jl=1,n_points
            WRITE(6,*) 'Point: ',jl
            IF (jit == 1) THEN
             WRITE(6,"(1x,i2,20(1x,1pG12.4))")                          &
                         jit-1,(zf(jl,jtr),jtr=1,jpctr)
             WRITE(6,"(1x,i2,20(1x,1pG12.4))") jit-1,                   &
                         ((zf(jl,jtr)+cdt*fdot(jl,jtr)),jtr=1,jpctr)
            END IF
            WRITE(6,"(1x,i2,20(1x,1pG12.4))")                           &
                            jit-1,(f(jl,jtr),jtr=1,jpctr),              &
                            (y(jl,i),i=1,2)
          END DO
        END IF
!
! DEPENDS ON: asad_diffun
        CALL asad_diffun( nl )
!
!  Temporary prod+loss array
        tmprc(1:n_points,:) = prod(1:n_points,1:jpctr) +                &
                              slos(1:n_points,1:jpctr)
!
        IF (errxo < raferr) THEN
          IF(jit >= ltrig_jit) THEN               ! 51, Set to 50 if want LTRIG  !!
            ndxraf = 4
            IF (printstatus >= prstatus_oper)                           &
            WRITE (6,                                                   &
      "('Convergence problems (',i3,1x,'iter) at lev=',i3,' pe=',i3)")  &
              jit, nlev, mype
          END IF
          GOTO 9999
        END IF

        IF (jit == nrsteps) THEN
          IF (printstatus >= prstatus_oper) WRITE(6,                    &
          "('Convergence not achieved in spimpmjp (iter',i3,') lev=',"//&
          "i3,'  pe=',i3,'; halving step')") jit, nlev, mype
          ndxraf = 2
          f = zf
          IF(jit >= ltrig_jit) THEN               ! 51, Set to 50 if want LTRIG  !!
            ndxraf = 4
            IF (printstatus >= prstatus_oper) WRITE(6,                  &
      "('Convergence problems (',i3,1x,'iter) at lev=',i3,'  pe=',i3)") &
              jit, nlev, mype
          END IF
          GOTO 9999
        END IF
!
!  Calculate fxo
        fxo = (f - zf)*deltt - fdot
!
!  Test for convergence
        errfl80 = .FALSE.
        errpl = 0.
        DO jl=1,n_points
          DO jtr=1,jpctr
            IF(ABS(tmprc(jl,jtr)) > rafeps)                             &
                   errpl=MAX(errpl,ABS(fxo(jl,jtr)/tmprc(jl,jtr)))
            IF (f(jl,jtr) > rafbig)  THEN
              errfl80 = .TRUE.
            END IF
          END DO
        END DO
        IF (errfl80) THEN
          ndxraf = 3
          f = zf
          IF(jit >= ltrig_jit) THEN               ! Set to 50 if want LTRIG  !!
            ndxraf = 4
            IF (printstatus >= prstatus_oper) WRITE(6,                  &
      "('Convergence problems (',i3,1x,'iter) at lev=',i3,'  pe=',i4)") &
              jit, nlev, mype
          END IF
          GOTO 9999
        END IF
        IF (errpl < rafpml) THEN
          IF(jit >= ltrig_jit) THEN               ! Set to 50 if want LTRIG  !!
            ndxraf = 4
            IF (printstatus >= prstatus_oper) WRITE(6,                  &
      "('Convergence problems (',i3,1x,'iter) at lev=',i3,'  pe=',i4)")&
              jit, nlev, mype
          END IF
          GOTO 9999
        END IF
!
!  Fill in the Jacobian, or just solve if lsvjac = .true.
        IF (lsvjac) THEN
! Sparse resolve routine
          CALL spresolv2(fxo,xoo,n_points,rafeps)
        ELSE
          CALL spfuljac(n_points)
!
          IF (ltrig .AND. printstatus == PrStatus_Diag) THEN
            WRITE(6,*) 'Iteration ',jit
            DO jl=1,n_points
              WRITE(6,*) 'Point: ',jl
              fj(jl,:,:) = 0.
              DO jtr = 1,jpctr
                DO itr = 1,jpctr
                  IF (pointer2(jtr,itr) > 0)                            &
                    fj(jl,jtr,itr) = spfj(jl,pointer2(jtr,itr))
                END DO
              END DO
              DO jtr=1,jpctr
                WRITE(6,"(1x,i2,20(1x,1pG12.4))") jtr,                 &
                   (fj(jl,jtr,itr),itr=1,jpctr)
              END DO
            END DO
          END IF

          xoo(:,:)=0.0     ! initialisation required
          CALL splinslv2(fxo,xoo,n_points,rafeps,rafbig)

          IF (ltrig .AND. printstatus == PrStatus_Diag) THEN
            DO jl=1,n_points
              WRITE(6,*) 'Point: ',jl
              WRITE(6,"(1x,a3,20(1x,1pG12.4))")                        &
                  'fxo',(fxo(jl,jtr),jtr=1,jpctr)
              WRITE(6,"(1x,a3,20(1x,1pG12.4))")                        &
                  'fdt',(fdot(jl,jtr),jtr=1,jpctr)
              WRITE(6,"(1x,a3,20(1x,1pG12.4))")                        &
                  'del',((f(jl,jtr)-zf(jl,jtr))*deltt,jtr=1,jpctr)
              WRITE(6,"(1x,a3,20(1x,1pG12.4))")                        &
                  'f  ',(f(jl,jtr),jtr=1,jpctr)
              WRITE(6,"(1x,a3,20(1x,1pG12.4))")                        &
                  'xoo',(xoo(jl,jtr),jtr=1,jpctr)
            END DO
! DEPENDS ON: asad_fuljac
            CALL asad_fuljac(n_points)
            WRITE(6,*) 'Itertion: ',jit
            DO jl=1,n_points
              DO jtr=1,jpctr
                zsum(jtr) = 0.
                DO itr=1,jpctr
                  zsum(jtr)=zsum(jtr)+fj(jl,jtr,itr)*xoo(jl,itr)
                END DO
                WRITE(6,"(1x,i2,20(1x,1pG12.4))") jtr,                 &
                        (fj(jl,jtr,itr)*xoo(jl,itr),itr=1,jpctr)
              END DO
              WRITE(6,"(1x,a3,20(1x,1pG12.4))") 'sum',                 &
                (zsum(jtr),jtr=1,jpctr)
            END DO
          END IF
!
        END IF
!
        errxo=0.
        ndxraf=0
        DO jtr=1,jpctr
          DO jl=1,n_points
!  Filter for negatives
!  Special kick for troublesome convergence
            IF (jit == 1) xoo(jl,jtr) = damp1*xoo(jl,jtr)
!  Calculate error
            xoo(jl,jtr) = MIN(MAX(xoo(jl,jtr),-rafbig),rafbig)
            IF (ABS(xoo(jl,jtr)) > 1.e-16)                              &
                errxo = MAX(errxo,ABS(xoo(jl,jtr)/                      &
                   MAX(f(jl,jtr),rafeps)))
!  New mixing ratios
            ztmp = f(jl,jtr) + xoo(jl,jtr)
!  Put limit on MAXimum correction for first few (6) iterations
            IF (jit < 7)                                                &
              ztmp = MAX(rafmin*f(jl,jtr),MIN(rafmax*f(jl,jtr),ztmp))
!  Filter negatives and zeros
            IF (ztmp == 0.) ztmp = rafeps
            IF (ztmp < 0.) THEN
              ztmp = rafeps
              ndxraf = ndxraf+1
            END IF
!  Final mixing ratios
            f(jl,jtr) = ztmp
          END DO
        END DO
! if 5 or more negatives, drop out and halve step
        IF(ndxraf > maxneg) THEN
          write(6,*) ndxraf, ' Negatives - exceeds maxneg'
          WRITE(6,                                                      &
            "(1x,'Too many negatives (>',i4,') in spimpmjp (iter',i3,"//&
            "')  lon=',i3,'  lat=',i3,'; halving step')")               &
            maxneg,jit, nlev, mype
          ndxraf=2
          f = zf
          IF(jit >= ltrig_jit) THEN               ! Set to 50 if want LTRIG  !!
            ndxraf = 4
            WRITE (6,"(1x,'Convergence problems (',i3,1x,"//            &
             "'iter) at lon=',i3,'  lat=',i3)")                         &
              jit, nlev, mype
          END IF
          GOTO 9999
        END IF

        ndxraf=0
      END DO
!
!  Failure to Converge - reset f's and exit
      WRITE(6,                                                          &
        "('Convergence not achieved in spimpmjp (iter',i3,')  lev=',"// &
        "i3,'  pe=',i3,'; halving step')") jit, nlev, mype
      ndxraf = 2
      f = zf
      IF(jit >= ltrig_jit) THEN               ! Set to 50 if want LTRIG  !!
         ndxraf = 4
         WRITE (6,                                                    &
      "('Convergence problems (',i3,1x,'iter) at lev=',i3,'  pe=',i3)")&
              jit, nlev, mype
      END IF

 9999 CONTINUE
      IF (lhook) CALL dr_hook('ASAD_SPIMPMJP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE asad_spimpmjp
