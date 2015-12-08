! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module file ukca_radaer_struct_mod
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA

! UKCA_RADAER:
!
! Defines maximum dimensions
! Defines type ukca_radaer_struct, the structure used by UKCA_RADAER
!
MODULE ukca_radaer_struct_mod

      IMPLICIT NONE

      !
      ! Maximum dimensions for automatic arrays:
      !  Number of modes
      !  Number of component per modes
      !  Total number of components
      ! Numbers should match nmodes and ncp in UKCA_MODE_SETUP (checked
      ! by routine ukca_radaer_init).
      !
      INTEGER, PARAMETER :: npd_ukca_mode         = 7
      INTEGER, PARAMETER :: npd_ukca_cpnt_in_mode = 6
      INTEGER, PARAMETER :: npd_ukca_cpnt =                             &
                              npd_ukca_mode * npd_ukca_cpnt_in_mode
      
      !
      ! Internal IDs for the types of UKCA aerosol modes.
      !
      INTEGER, PARAMETER :: IP_UKCA_MODE_NUCLEATION = 0
      INTEGER, PARAMETER :: IP_UKCA_MODE_AITKEN     = 1
      INTEGER, PARAMETER :: IP_UKCA_MODE_ACCUM      = 2
      INTEGER, PARAMETER :: IP_UKCA_MODE_COARSE     = 3
      
      !
      ! Internal IDs for the types of UKCA aerosol components.
      !
      ! Note: When adding an aerosol component, also increase 
      !       npd_ukca_maxcomptype in module ukca_radaer_precalc_mod and 
      !       update the pre-computed file read by ukca_radaer_read_precalc. 
      !       Any new aerosol component should be added before IP_UKCA_WATER.
      ! 
      ! Water is not an aerosol species, but is included here
      ! as it behaves like one in ukca_radaer_band_average(). 
      ! However, no aerosol component should be of type IP_UKCA_WATER.
      !
      INTEGER, PARAMETER :: IP_UKCA_SULPHATE      = 1
      INTEGER, PARAMETER :: IP_UKCA_BLACKCARBON   = 2
      INTEGER, PARAMETER :: IP_UKCA_ORGANICCARBON = 3
      INTEGER, PARAMETER :: IP_UKCA_SEASALT       = 4
      INTEGER, PARAMETER :: IP_UKCA_DUST          = 5
      INTEGER, PARAMETER :: IP_UKCA_SECONDORGANIC = 6
      
      INTEGER, PARAMETER :: IP_UKCA_WATER         = 7
      
      !
      ! Main structure holding all the variables needed for 
      ! interacting UKCA aerosols with radiation.
      
      TYPE ukca_radaer_struct
        
        !
        ! Information about UKCA aerosol modes
        !
        
        ! Actual number of modes
        INTEGER n_mode
        
        ! Type of mode (i.e. nucleation, Aitken, accum, or coarse)
        INTEGER, DIMENSION(npd_ukca_mode) :: i_mode_type
        
        ! Solubility of mode (soluble if true)
        LOGICAL, DIMENSION(npd_ukca_mode) :: l_soluble
        
        ! Lower and upper limits on the geometric mean diameter (m)
        ! in each mode
        REAL, DIMENSION(npd_ukca_mode) :: d0low
        REAL, DIMENSION(npd_ukca_mode) :: d0up
        
        ! Geometric standard deviation in this mode
        REAL, DIMENSION(npd_ukca_mode) :: sigma

        ! STASH code, D1 address, number of levels and total length 
        ! of the D1 fields corresponding to modal dry diameter.
        INTEGER, DIMENSION(npd_ukca_mode) :: stashcode_dry,  &
                                             d1_address_dry, &
                                             d1_nlevs_dry,   &
                                             d1_length_dry
        
        ! STASH code, D1 address, number of levels and total length 
        ! of the D1 fields corresponding to modal wet diameter.
        INTEGER, DIMENSION(npd_ukca_mode) :: stashcode_wet,  &
                                             d1_address_wet, &
                                             d1_nlevs_wet,   &
                                             d1_length_wet
        
        ! STASH code, D1 address, number of levels and total length
        ! of the D1 fields corresponding to modal density.
        INTEGER, DIMENSION(npd_ukca_mode) :: stashcode_rho,  &
                                             d1_address_rho, &
                                             d1_nlevs_rho,   &
                                             d1_length_rho

        ! STASH code, D1 address, number of levels and total length
        ! of the D1 fields corresponding to water volume in each mode.
        INTEGER, DIMENSION(npd_ukca_mode) :: stashcode_wtv,  &
                                             d1_address_wtv, &
                                             d1_nlevs_wtv,   &
                                             d1_length_wtv
        
        ! STASH code, D1 address, number of levels, total length 
        ! and halo type of the D1 fields corresponding to 
        ! modal number concentrations.
        INTEGER, DIMENSION(npd_ukca_mode) :: stashcode_nbr,  &
                                             d1_address_nbr, &
                                             d1_nlevs_nbr,   &
                                             d1_length_nbr,  &
                                             d1_halo_type_nbr
        
        ! Number of components in each mode and index of each
        ! component in array ukca_cpnt_info
        INTEGER, DIMENSION(npd_ukca_mode) :: n_cpnt_in_mode
        INTEGER, DIMENSION(npd_ukca_cpnt_in_mode, npd_ukca_mode) ::     &
          i_cpnt_index
        
        ! Modal diameter of the dry aerosol (m)
        REAL, POINTER, DIMENSION(:, :, :, :) :: dry_diam
        
        ! Modal diameter of the wet aerosol (m)
        REAL, POINTER, DIMENSION(:, :, :, :) :: wet_diam
        
        ! Modal densities (kg/m3)
        REAL, POINTER, DIMENSION(:, :, :, :) :: modal_rho
        
        ! Modal volumes (including water for soluble modes)
        REAL, POINTER, DIMENSION(:, :, :, :) :: modal_vol

        ! Fractional volume of water in each mode
        REAL, POINTER, DIMENSION(:, :, :, :) :: modal_wtv
        
        ! Modal number concentrations
        REAL, POINTER, DIMENSION(:, :, :, :) :: modal_nbr
        
        !
        ! Information about UKCA aerosol components
        !
        
        ! Actual number of components
        INTEGER n_cpnt
        
        ! Type of component (e.g. sulphate)
        INTEGER, DIMENSION(npd_ukca_cpnt) :: i_cpnt_type
        
        ! Mass density of each component (kg/m3)
        REAL, DIMENSION(npd_ukca_cpnt) :: density
        
        ! Array index of the mode this component belongs to
        INTEGER, DIMENSION(npd_ukca_cpnt) :: i_mode
        
        ! STASH code, D1 address, number of levels, total length 
        ! and halo type of the D1 fields corresponding to the 
        ! mass-mixing ratio of each component.
        INTEGER, DIMENSION(npd_ukca_cpnt) :: stashcode_mmr,   &
                                             d1_address_mmr,  &
                                             d1_nlevs_mmr,    &
                                             d1_length_mmr,   &
                                             d1_halo_type_mmr
        
        ! STASH code, D1 address, number of levels and total length
        ! of the D1 fields corresponding to the fractional volume of
        ! each component.
        INTEGER, DIMENSION(npd_ukca_cpnt) :: stashcode_cvl,  &
                                             d1_address_cvl, &
                                             d1_nlevs_cvl,   &
                                             d1_length_cvl
        
        ! Component mass-mixing ratio (kg/kg)
        REAL, POINTER, DIMENSION(:, :, :, :) :: mix_ratio
        
        ! Component volumes
        REAL, POINTER, DIMENSION(:, :, :, :) :: comp_vol
      
      END TYPE ukca_radaer_struct

END MODULE ukca_radaer_struct_mod
