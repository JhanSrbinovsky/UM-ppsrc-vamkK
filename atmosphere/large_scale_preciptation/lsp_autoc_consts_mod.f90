! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lsp_autoc_consts_mod

! Description:
! Holds autoconversion constants required by the large-scale
! precipitation scheme
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

  IMPLICIT NONE

! ----------------------------------------------------------------------
!      AUTOCONVERSION TERMS
! ----------------------------------------------------------------------

!     LOGICAL, PARAMETER :: L_AUTOCONV_MURK is set in UMUI
! Set to .TRUE. to calculate droplet concentration from MURK aerosol,
! which will override L_USE_SULPHATE_AUTOCONV (second indirect effect
! of sulphate aerosol). If both are .FALSE., droplet concentrations
! from comdeck C_MICRO are used to be consistent with the values
      ! used in the radiation scheme.

      ! This next set of parameters is to allow the 3B scheme to
      ! be replicated at 3C/3D
        ! Inhomogeneity factor for autoconversion rate
  REAL,PARAMETER:: inhomog_rate=1.0

        ! Inhomogeneity factor for autoconversion limit
  REAL,PARAMETER:: inhomog_lim=1.0

        ! Threshold droplet radius for autoconversion
  REAL,PARAMETER:: r_thresh=7.0e-6
      ! End of 3B repeated code

      !Do not alter R_AUTO and N_AUTO since these values are effectively
      ! hard wired into a numerical approximation in the autoconversion
      ! code. EC_AUTO will be multiplied by CONSTS_AUTO

      ! Threshold radius for autoconversion
  REAL, PARAMETER :: r_auto=20.0e-6

      ! Critical droplet number for autoconversion
  REAL, PARAMETER :: n_auto=1000.0

      ! Collision coalesence efficiency for autoconversion
!      REAL, PARAMETER :: EC_AUTO is set in UMUI

      ! The autoconversion powers define the variation of the rate with
      ! liquid water content and droplet concentration. The following are
      ! from Tripoli and Cotton

      !  Dependency of autoconversion rate on droplet concentration
  REAL, PARAMETER :: power_droplet_auto=-0.33333

      ! Dependency of autoconversion rate on water content
  REAL, PARAMETER :: power_qcl_auto=2.33333

      ! Dependency of autoconversion rate on air density
  REAL, PARAMETER :: power_rho_auto=1.33333

      ! CONSTS_AUTO = (4 pi)/( 18 (4 pi/3)^(4/3)) g /  mu (rho_w)^(1/3)
      ! See UM documentation paper 26, equation P26.132

      ! Combination of physical constants
  REAL, PARAMETER :: consts_auto=5907.24

      ! Quantites for calculation of drop number by aerosols.
      ! Need only set if L_AUTOCONV_MURK=.TRUE.  See file C_VISBTY

      ! Scaling concentration (m-3) in aerosol number concentration
      ! Value from Haywood et al (2008) QJRMS paper
  REAL, PARAMETER :: n0_haywood = 200.0e7
      ! Value from Clark et al (2008) QJRMS paper
  REAL, PARAMETER :: n0_clark   = 500.0e6

      ! Scaling mass (kg/kg) in aerosol number calculation from aersol mass
      ! M0_MURK = 4/3 * pi * N0_MURK * (r0)^3 * (rho_aerosol/rho_air)
      ! r0 = 0.11E-6 in Haywood et al (2008) Scheme
      ! rho_aerosol = 1700.0 in Haywood et al (2008) scheme
      ! rho_air = 1.0 in Haywood et al (2008) scheme
     ! Value from Haywood et al (2008) QJRMS paper
  REAL, PARAMETER :: m0_haywood = 1.8956e-8
      ! Value from Clark et al (2008) QJRMS paper
  REAL, PARAMETER :: m0_clark   = 1.4584e-8 

      ! Power in droplet number calculation from aerosols
  REAL, PARAMETER :: power_murk=0.5

      ! Ice water content threshold for graupel autoconversion (kg/m^3)
  REAL, PARAMETER :: auto_graup_qcf_thresh = 3.e-4

      ! Temperature threshold for graupel autoconversion (degC)
  REAL, PARAMETER :: auto_graup_t_thresh = -4.0

      ! Temperature threshold for graupel autoconversion
  REAL, PARAMETER :: auto_graup_coeff = 0.5

!-----------------------------------------------------------------------
! New parameters added for the droplet tapering process
! Put in this module for now, eventually may plumb from the UMUI if
! there is a user requirement for these.
!-----------------------------------------------------------------------

      ! Altitude (m) of where the droplet number goes to its minimum
      ! value
  REAL, PARAMETER :: z_low_nd     = 2000.0

      ! Droplet number at altitudes of z_low_nd and above
  REAL, PARAMETER :: min_drop_alt = 100.0e6

!-----------------------------------------------------------------------
! Parameters for autoconversion scheme if l_warm_new=.true.
! Currently these are for Khairoutdinov & Kogan (2000, MWR)
!-----------------------------------------------------------------------
      ! Pre-factor
  REAL, PARAMETER :: aut_pref = 1350.0
      ! qcl power
  REAL, PARAMETER :: aut_qc   = 2.47
      ! n_drop power
  REAL, PARAMETER :: aut_nc   = -1.79

END MODULE lsp_autoc_consts_mod
