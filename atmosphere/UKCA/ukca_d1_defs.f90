! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved. 
! For further details please refer to the file COPYRIGHT.txt 
! which you should have received as part of this distribution. 
! *****************************COPYRIGHT*******************************
!
! Description:
!  Module to define fields from D1 needed for UKCA
!  these are set in UKCA_SETD1DEFS
!
!  Part of the UKCA model. UKCA is a community model supported
!  by The Met Office and NCAS, with components provided initially
!  by The University of Cambridge, University of Leeds and
!  The Met. Office.  See www.ukca.ac.uk
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
      MODULE UKCA_D1_DEFS
      IMPLICIT NONE

!     No of Prognostics and Diagnostics required

      INTEGER, SAVE :: Nukca_D1items      ! Size of UkcaD1codes array
      INTEGER, SAVE :: n_all_tracers=150  ! max number of tracers
      INTEGER, SAVE :: n_use_tracers      ! no of tracers used
      INTEGER, SAVE :: n_use_emissions    ! no of emissions used
      INTEGER, SAVE :: n_3d_emissions     ! no of 3-D emissions used
      INTEGER, SAVE :: n_in_progs         ! No of prognostics required
      INTEGER, SAVE :: n_in_diags0        ! No of diagnostics (sect 0)
      INTEGER, SAVE :: n_in_diags1        ! No of diagnostics (sect 1)
      INTEGER, SAVE :: n_in_diags2        ! No of diagnostics (sect 2)
      INTEGER, SAVE :: n_in_diags3        ! No of diagnostics (sect 3)
      INTEGER, SAVE :: n_in_diags4        ! No of diagnostics (sect 4)
      INTEGER, SAVE :: n_in_diags5        ! No of diagnostics (sect 5)
      INTEGER, SAVE :: n_in_diags8        ! No of diagnostics (sect 8)
      INTEGER, SAVE :: n_in_diags15       ! No of diagnostics (sect 15)
      INTEGER, SAVE :: n_in_diags30       ! No of diagnostics (sect 30)
      INTEGER, SAVE :: n_in_diags33       ! No of diagnostics (sect 33)
      INTEGER, SAVE :: n_in_diags34       ! No of diagnostics (sect 34)
      INTEGER, SAVE :: n_in_diags38       ! No of diagnostics (sect 38)
      INTEGER, SAVE :: n_in_diags50       ! No of diagnostics (sect 50)
      INTEGER, SAVE :: n_out_diags        ! Not used
      INTEGER, SAVE :: n_emiss_first      ! Position of emissions in section
      INTEGER, SAVE :: n_emiss_last       ! Position of emissions in section
      INTEGER, SAVE :: idiag_first        ! Position of first chem diag
      INTEGER, SAVE :: idiag_last         ! Position of last chem diag
      INTEGER, SAVE :: iemiss_first       ! item for 1st emission diag
      INTEGER, SAVE :: iemiss_last        ! item for last emission diag
      INTEGER, SAVE :: istrat_first       ! item for 1st strat flux diag
      INTEGER, SAVE :: istrat_last        ! item for last strat flux diag
      INTEGER, SAVE :: imode_first        ! item for 1st MODE diag
      INTEGER, SAVE :: imode_last         ! item for last MODE diag
      INTEGER, PARAMETER :: UKCA_sect=34       ! stash section for UKCA
      INTEGER, PARAMETER :: UKCA_ems_sect=0    ! stash section for UKCA ems
      INTEGER, PARAMETER :: MODE_diag_sect=38  ! stash section for MODE diags
      INTEGER, PARAMETER :: ASAD_diag_sect=50  ! stash section for ASAD diags
      INTEGER, PARAMETER :: item1_chem_diags = 151  ! 1st item No. for Chem diags
      INTEGER, PARAMETER :: item2_chem_diags = 172  ! last item No. for Chem diags  
      INTEGER, PARAMETER :: item1_asad_diags = 1    ! 1st item No. for ASAD diags (s50)
      INTEGER, PARAMETER :: item1_stratflux  = 1    ! 1st item No. Strat fluxes (s38)
      INTEGER, PARAMETER :: item1_mode_diags = 201  ! 1st item No. for MODE diags (s38) 
      INTEGER, PARAMETER :: n_boundary_vals  = 7    ! No. lower boundary vals for strat chem

! Indices for chem_diags array  (use for N-R only)
      INTEGER, PARAMETER :: icd_o3p      = 1      ! index for O(3P)
      INTEGER, PARAMETER :: icd_o1d      = 2      ! index for O(1D)
      INTEGER, PARAMETER :: icd_no2      = 3      ! index for NO2
      INTEGER, PARAMETER :: icd_bro      = 4      ! index for BrO
      INTEGER, PARAMETER :: icd_hcl      = 5      ! index for HCl
      INTEGER, PARAMETER :: icd_cly      = 6      ! index for Cly
      INTEGER, PARAMETER :: icd_surfarea = 7      ! index for aerosol surface area
      INTEGER, PARAMETER :: icd_NAT      = 8      ! index for NAT
      INTEGER, PARAMETER :: icd_O3col    = 9      ! index for Ozone column
      INTEGER, PARAMETER :: icd_rc_het1  = 10     ! index for rate for heterogeneous N2O5 loss
      INTEGER, PARAMETER :: icd_rc_het2  = 11     ! index for rate for heterogeneous HO2 loss
      INTEGER, PARAMETER :: icd_cdnc     = 12     ! index for CDNC
      INTEGER, PARAMETER :: icd_cdnc3    = 13     ! index for (CDNC)^1/3

      TYPE CODE
       INTEGER :: section        ! section code
       INTEGER :: item           ! item code
       INTEGER :: n_levels       ! number of levels
       INTEGER :: address        ! address in D1
       INTEGER :: length         ! length of field
       INTEGER :: halo_type      ! field halo type
       INTEGER :: grid_type      ! grid type
       INTEGER :: field_type     ! field grid type
       INTEGER :: len_dim1       ! length of array dim1
       INTEGER :: len_dim2       ! length of array dim2
       INTEGER :: len_dim3       ! length of array dim3
       LOGICAL :: prognostic     ! prognostic t/f
       LOGICAL :: required       ! t/f
       CHARACTER(len=10) :: name ! species name
      ENDTYPE CODE

!     Number of tracers, emissions and diagnostics, set according
!     to UMUI choices in UKCA_SETD1DEFS

      INTEGER, SAVE   :: n_chem_tracers     ! No. tracers for chemistry
      INTEGER, SAVE   :: n_aero_tracers     ! No. tracers for chemistry
      INTEGER, SAVE   :: n_nonchem_tracers  ! No. tracers non-chemistry
      INTEGER, SAVE   :: n_chem_emissions   ! "   emissions        "
      INTEGER, SAVE   :: n_chem_diags       ! "   diagnostics      "
      INTEGER, SAVE   :: nmax_strat_fluxdiags  ! Max no strat flux diags
      INTEGER, SAVE   :: nmax_mode_diags    ! Max no MODE diags
      INTEGER, SAVE   :: n_strat_fluxdiags  ! No. strat flux diags
      INTEGER, SAVE   :: n_MODE_tracers     ! No. tracers for MODE
      INTEGER, SAVE   :: n_MODE_emissions   ! "   emissions   "
      INTEGER, SAVE   :: n_MODE_diags       ! "   diagnostics "
      INTEGER, SAVE   :: n_dust_tracers     ! No. tracers for dust
      INTEGER, SAVE   :: n_dust_emissions   ! "   emissions   "
      INTEGER, SAVE   :: n_dust_diags       ! "   diagnostics "
      INTEGER, SAVE   :: nr_therm           ! "   thermal reactions
      INTEGER, SAVE   :: nr_phot            ! "   photolytic reactions

!     list of prognostic/diagnostics from/to D1

      TYPE(CODE), DIMENSION(:), ALLOCATABLE, SAVE :: UkcaD1Codes

!     Names for tracers which have surface emissions

      CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE, SAVE :: em_chem_spec
      CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE, SAVE :: em_MODE_spec
      CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE, SAVE :: em_dust_spec
! Lower BCs for stratospheric species
      CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE, SAVE :: lbc_spec   ! species names
      REAL,              DIMENSION(:), ALLOCATABLE, SAVE :: lbc_mmr    ! mixing ratios

!     Names for defined tracers

      CHARACTER(LEN=10), DIMENSION(:), ALLOCATABLE, SAVE :: nm_spec

!     Index arrays for tracers and emissions

      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE           :: tr_index
      INTEGER, DIMENSION(:), ALLOCATABLE, SAVE           :: em_index

!     UKCA Logicals - will eventually be replaced with umui/namelist
!      The flux logicals are now turned on via STASH panel requests.
      LOGICAL, SAVE :: L_ukca_stratflux = .false.  ! T for stratospheric fluxes
      LOGICAL, SAVE :: L_ukca_mode_diags= .false.  ! T for MODE diags.
      LOGICAL, SAVE :: L_ukca_emiss_diags= .false. ! T for emission diags.

! STASH Codes of species used for  coupling of UKCA with radiation scheme 
! and sulphur cycle of UM
 
!! Note: TO BE UPDATED FOR ANY SEC 34 ITEM CHANGES     !!

! species used for feedback to radiation scheme
      INTEGER, PARAMETER :: i_ukca_grg_o3  =  1, i_ukca_grg_ch4 = 9 ,           &
                            i_ukca_grg_n2o = 49, i_ukca_grg_f11 = 55,           &
                            i_ukca_grg_f12 = 56, i_ukca_grg_f113= 64,           &
                            i_ukca_grg_h22 = 65, i_ukca_grg_cf3chf2 = -1,       &
                            i_ukca_grg_chf2chf2 = -1

! Item numbers for UKCA oxidants used in the CLASSIC sulphur cycle. 
! Set in ukca_setup_chem - address of OH and HO2 varies with solver

!       001 : O3 MASS MIXING RATIO AFTER TSTEP 
!       007 : HONO2 MASS MIXING RATIO AFTER TSTEP 
!       008 : H2O2 MASS MIXING RATIO AFTER TSTEP 
!           : OH MASS MIXING RATIO  AFTER TSTEP 
!           : HO2 MASS MIXING RATIO AFTER TSTEP 
      INTEGER, SAVE :: ukca_item_sulpc(5)

      END MODULE UKCA_D1_DEFS
