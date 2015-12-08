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
!  Description:
!    Computes steady-state species for Newton-Raphson integrator.
!    Part of the ASAD chemical solver.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!     ASAD: ycn                      Version: steady.f 4.1 20/7/07
!
!     Purpose
!     -------
!
!  Routine to explicitly define steady state expressions - prevents
!  generality, but removes need for sluggish iteration round family
!  stuff.
!
!  Important Notes:
!     1) Ordering of calculations is important - need to avoid feedbacks!
!     2) Needs to be rewritten whenever reaction numbering is changed.
!
!  Additions:
!     To improve generality, reactions involving steady state species
!  are selected in 'setsteady' and loaded into nss* integer arrays,
!  which are then used in this routine. This causes a very slight
!  increase in CPU time, but removes the need to rewrite the routine
!  whenever a new species is added or a reaction changed.
!
!                                            Oliver   (3 Feb 1998)
! We add more general terms for the steady-state species. It is assumed
! that O(1D), O(3P), H and N can be put in steady state.
!
!
!     Method
!     ------
!
!   We sum up production and loss terms for SS species, and divide.
!   Moreover, corresponding terms in the Jacobian are calculated that
!   account for the dependence of steady-state variables on tracer
!   variables.
!
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      SUBROUTINE asad_steady( kl )

      USE ASAD_MOD,               ONLY:  deriv, y, rk, peps,            &
                                         nspi, nssi, nssrt, nssrx,      &
                                         nssri, nsspt, nsspi, nsst,     &
                                         jlst, nspo1d, nspo3, nspoh,    &
                                         nspo3p, nsph, nuni,            &
                                         nspho2, nspno, nspn, nss_o3p,  &
                                         nss_o1d, nss_n, nss_h,         &
                                         o3p_in_ss, n_in_ss, h_in_ss
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
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
      INTEGER, INTENT(IN) :: kl            ! No. of points

! Local variables
      INTEGER, PARAMETER :: n_o3=1           ! indicies for dssden etc
      INTEGER, PARAMETER :: n_oh=2
      INTEGER, PARAMETER :: n_ho2=3
      INTEGER, PARAMETER :: n_no=4

      INTEGER :: jl
      INTEGER :: jl1
      INTEGER :: jl2
      INTEGER :: jr
      INTEGER :: ix
      INTEGER :: i
      INTEGER :: j

      REAL :: ssnum(theta_field_size)
      REAL :: ssden(theta_field_size)

! Add here derivatives w.r.t. O3, OH, and HO2, and NO of numerator and
! denominator of steady state species
      REAL :: dssnum(theta_field_size,4)
      REAL :: dssden(theta_field_size,4)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!
! Set up loops correctly
      IF (lhook) CALL dr_hook('ASAD_STEADY',zhook_in,zhook_handle)
      jl1 = 1
      jl2 = kl
      IF (kl == 1) THEN
        jl1 = jl1+jlst
        jl2 = jl1
      END IF

! Initialise DERIV - now done in ASAD_INIT, to 1.0
!     deriv(:,:,:) = 0.0

! Loop through steady state species

      DO ix = 1,nsst
        ssnum = 0.
        ssden = 0.
        dssden = 0.
        dssnum = 0.

! Production terms
        DO jr = 1,nsspt(ix)
          i = nsspi(ix,jr)
          IF (i <= nuni) THEN
            ssnum(jl1:jl2) = ssnum(jl1:jl2) +                           &
                rk(jl1:jl2,i)*y(jl1:jl2,nspi(i,1))

            IF ((ix < 5) .AND. (nspi(i,1) == nspo3 ))                   &
! add terms to derivative for d(j[O3])/d[O3] = j_o3
              dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +             &
                                       rk(jl1:jl2,i)
            IF ((ix < 5) .AND. (nspi(i,1) == nspno ))                   &
! add terms to derivative for d(j[NO])/d[NO] = j_no
              dssnum(jl1:jl2,n_no) = dssnum(jl1:jl2,n_no) +             &
                                       rk(jl1:jl2,i)
          ELSE
            ssnum(jl1:jl2) = ssnum(jl1:jl2) +                           &
                rk(jl1:jl2,i)*y(jl1:jl2,nspi(i,1))*y(jl1:jl2,nspi(i,2))
            IF (ix < 5) THEN

! add terms for derivative w.r.t. ozone.
              IF (nspi(i,1) == nspo1d)                                  &
                dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +           &
                     rk(jl1:jl2,i)                                      &
                    *y(jl1:jl2,nspi(i,2))*deriv(jl1:jl2,nss_o1d,n_o3)

              IF (nspi(i,2) == nspo1d)                                  &
                  dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))*deriv(jl1:jl2,nss_o1d,n_o3)

              IF (nspi(i,1) == nspo3p)                                  &
                  dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +         &
                     rk(jl1:jl2,i)                                      &
                    *y(jl1:jl2,nspi(i,2))*deriv(jl1:jl2,nss_o3p,n_o3)

              IF (nspi(i,2) == nspo3p)                                  &
                  dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))*deriv(jl1:jl2,nss_o3p,n_o3)

              IF (nspi(i,1) == nspo3)                                   &
                dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))

              IF (nspi(i,2) == nspo3)                                   &
                dssnum(jl1:jl2,n_o3) = dssnum(jl1:jl2,n_o3) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))

! add terms for derivative w.r.t OH
              IF (nspi(i,1) == nspo3p)                                  &
! add terms to derivates for d(a[A][B])
                dssnum(jl1:jl2,n_oh) = dssnum(jl1:jl2,n_oh) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))*deriv(jl1:jl2,nss_o3p,n_oh)

              IF (nspi(i,2) == nspo3p)                                  &
! add terms to derivates for d(a[O2][O1D])/d[O3] and b[N2][O1D]
                dssnum(jl1:jl2,n_oh) = dssnum(jl1:jl2,n_oh) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))*deriv(jl1:jl2,nss_o3p,n_oh)

              IF (nspi(i,1) == nspoh)                                   &
! add terms to derivates for d(a[O2][O1D])/d[O3] and b[N2][O1D]
                dssnum(jl1:jl2,n_oh) = dssnum(jl1:jl2,n_oh) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))

              IF (nspi(i,2) == nspoh)                                   &
! add terms to derivates for d(a[O2][O1D])/d[O3] and b[N2][O1D]
                dssnum(jl1:jl2,n_oh) = dssnum(jl1:jl2,n_oh) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))

! add terms for derivative w.r.t HO2
              IF (nspi(i,1) == nspo3p)                                  &
                dssnum(jl1:jl2,n_ho2) = dssnum(jl1:jl2,n_ho2) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))*deriv(jl1:jl2,nss_o3p,n_ho2)

              IF (nspi(i,2) == nspo3p)                                  &
                dssnum(jl1:jl2,n_ho2) = dssnum(jl1:jl2,n_ho2) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))*deriv(jl1:jl2,nss_o3p,n_ho2)

              IF (nspi(i,1) == nspho2)                                  &
                dssnum(jl1:jl2,n_ho2) = dssnum(jl1:jl2,n_ho2) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))

              IF (nspi(i,2) == nspho2)                                  &
                dssnum(jl1:jl2,n_ho2) = dssnum(jl1:jl2,n_ho2) +         &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))

! add terms for derivative w.r.t NO
              IF (nspi(i,1) == nspno)                                   &
                dssnum(jl1:jl2,n_no) = dssnum(jl1:jl2,n_no) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,2))

              IF (nspi(i,2) == nspno)                                   &
                dssnum(jl1:jl2,n_no) = dssnum(jl1:jl2,n_no) +           &
                    rk(jl1:jl2,i)                                       &
                    *y(jl1:jl2,nspi(i,1))

            END IF
          END IF
        END DO       ! jr
!
! Destruction terms
        DO jr = 1,nssrt(ix)
          i = nssri(ix,jr)
          j = nssrx(ix,jr)
          IF (i <= nuni) THEN
            ssden(jl1:jl2) = ssden(jl1:jl2) + rk(jl1:jl2,i)
          ELSE
            ssden(jl1:jl2) = ssden(jl1:jl2) +                           &
              rk(jl1:jl2,i) * y(jl1:jl2,nspi(i,j))
            IF (ix < 5) THEN
              IF (nspi(i,j) == nspo3 )                                  &
                dssden(jl1:jl2,n_o3 ) = dssden(jl1:jl2,n_o3 ) +         &
                  rk(jl1:jl2,i)
              IF (nspi(i,j) == nspoh )                                  &
                dssden(jl1:jl2,n_oh ) = dssden(jl1:jl2,n_oh ) +         &
                  rk(jl1:jl2,i)
              IF (nspi(i,j) == nspho2)                                  &
                dssden(jl1:jl2,n_ho2) = dssden(jl1:jl2,n_ho2) +         &
                  rk(jl1:jl2,i)
              IF (nspi(i,j) == nspno )                                  &
                dssden(jl1:jl2,n_no ) = dssden(jl1:jl2,n_no ) +         &
                  rk(jl1:jl2,i)
            END IF
          END IF
        END DO    ! jr
!
! Steady state and derivatives of steady state
        y(jl1:jl2,nssi(ix)) = ssnum(jl1:jl2)/ssden(jl1:jl2)
        IF (ix < 5) THEN
          DO jr =1,4
            deriv(jl1:jl2,ix,jr) =                                      &
                (ssden(jl1:jl2)*dssnum(jl1:jl2,jr) -                    &
                 ssnum(jl1:jl2)*dssden(jl1:jl2,jr))/                    &
                 (ssden(jl1:jl2) * ssden(jl1:jl2))
          END DO    ! jr
        END IF
      END DO  ! ix

! rescale deriv to mean [O3]/[O] * d[O]/d[O3], where [O] = [O(1D)] or [O(3P)]
! for O(1D), and O(3P), N and H when these are SS species

      WHERE (y(jl1:jl2,nspo1d) > peps)
        deriv(jl1:jl2,nss_o1d,n_o3) = deriv(jl1:jl2,nss_o1d,n_o3 )*     &
                                  y(jl1:jl2,nspo3 )/y(jl1:jl2,nspo1d)
        deriv(jl1:jl2,nss_o1d,n_oh) = deriv(jl1:jl2,nss_o1d,n_oh )*     &
                                  y(jl1:jl2,nspoh )/y(jl1:jl2,nspo1d)
        deriv(jl1:jl2,nss_o1d,n_ho2)= deriv(jl1:jl2,nss_o1d,n_ho2)*     &
                                  y(jl1:jl2,nspho2)/y(jl1:jl2,nspo1d)
        deriv(jl1:jl2,nss_o1d,n_no )= deriv(jl1:jl2,nss_o1d,n_no )*     &
                                  y(jl1:jl2,nspno )/y(jl1:jl2,nspo1d)
      ELSEWHERE
        deriv(jl1:jl2,nss_o1d,n_o3)  = 1.
        deriv(jl1:jl2,nss_o1d,n_oh)  = 1.
        deriv(jl1:jl2,nss_o1d,n_ho2) = 1.
        deriv(jl1:jl2,nss_o1d,n_no)  = 1.
      ENDWHERE

      IF (o3p_in_ss) THEN
        WHERE (y(jl1:jl2,nspo3p) > peps)
          deriv(jl1:jl2,nss_o3p,n_o3 )= deriv(jl1:jl2,nss_o3p,n_o3 )*   &
                                  y(jl1:jl2,nspo3 )/y(jl1:jl2,nspo3p)
          deriv(jl1:jl2,nss_o3p,n_oh )= deriv(jl1:jl2,nss_o3p,n_oh )*   &
                                  y(jl1:jl2,nspoh )/y(jl1:jl2,nspo3p)
          deriv(jl1:jl2,nss_o3p,n_ho2)= deriv(jl1:jl2,nss_o3p,n_ho2)*   &
                                  y(jl1:jl2,nspho2)/y(jl1:jl2,nspo3p)
          deriv(jl1:jl2,nss_o3p,n_no )= deriv(jl1:jl2,nss_o3p,n_no )*   &
                                  y(jl1:jl2,nspno )/y(jl1:jl2,nspo3p)
        ELSEWHERE
          deriv(jl1:jl2,nss_o3p,n_o3)  = 1.
          deriv(jl1:jl2,nss_o3p,n_oh)  = 1.
          deriv(jl1:jl2,nss_o3p,n_ho2) = 1.
          deriv(jl1:jl2,nss_o3p,n_no)  = 1.
        ENDWHERE
      END IF

      IF (n_in_ss) THEN
        WHERE (y(jl1:jl2,nspn  ) > peps)
          deriv(jl1:jl2,nss_n,n_o3) = deriv(jl1:jl2,nss_n,n_o3 )*       &
                                y(jl1:jl2,nspo3 )/y(jl1:jl2,nspn)
          deriv(jl1:jl2,nss_n,n_oh) = deriv(jl1:jl2,nss_n,n_oh )*       &
                                y(jl1:jl2,nspoh )/y(jl1:jl2,nspn)
          deriv(jl1:jl2,nss_n,n_ho2)= deriv(jl1:jl2,nss_n,n_ho2)*       &
                                y(jl1:jl2,nspho2)/y(jl1:jl2,nspn)
          deriv(jl1:jl2,nss_n,n_no )= deriv(jl1:jl2,nss_n,n_no )*       &
                                y(jl1:jl2,nspno )/y(jl1:jl2,nspn)
        ELSEWHERE
          deriv(jl1:jl2,nss_n,n_o3)  = 1.
          deriv(jl1:jl2,nss_n,n_oh)  = 1.
          deriv(jl1:jl2,nss_n,n_ho2) = 1.
          deriv(jl1:jl2,nss_n,n_no)  = 1.
        ENDWHERE
      END IF

      IF (h_in_ss) THEN
        WHERE (y(jl1:jl2,nsph  ) > peps)
          deriv(jl1:jl2,nss_h,n_o3) = deriv(jl1:jl2,nss_h,n_o3 )*       &
                              y(jl1:jl2,nspo3 )/y(jl1:jl2,nsph  )
          deriv(jl1:jl2,nss_h,n_oh) = deriv(jl1:jl2,nss_h,n_oh )*       &
                              y(jl1:jl2,nspoh )/y(jl1:jl2,nsph  )
          deriv(jl1:jl2,nss_h,n_ho2)= deriv(jl1:jl2,nss_h,n_ho2)*       &
                              y(jl1:jl2,nspho2)/y(jl1:jl2,nsph  )
          deriv(jl1:jl2,nss_h,n_no )= deriv(jl1:jl2,nss_h,n_no )*       &
                              y(jl1:jl2,nspno )/y(jl1:jl2,nsph  )
        ELSEWHERE
          deriv(jl1:jl2,nss_h,n_o3)  = 1.
          deriv(jl1:jl2,nss_h,n_oh)  = 1.
          deriv(jl1:jl2,nss_h,n_ho2) = 1.
          deriv(jl1:jl2,nss_h,n_no)  = 1.
        ENDWHERE
      END IF

      IF (lhook) CALL dr_hook('ASAD_STEADY',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE asad_steady
