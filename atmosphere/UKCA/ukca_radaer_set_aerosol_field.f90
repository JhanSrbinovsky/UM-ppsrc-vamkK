! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!+ Copy UKCA-MODE aerosol fields to arrays that follow the radiation
!+ convention (top-to-bottom in the vertical, gathered horizontal
!+ points)
!
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
SUBROUTINE ukca_radaer_set_aerosol_field(                               &
   i_gather, L_EXTRA_TOP,                                               &
   n_layer, n_profile,                                                  &
   ukca_dim1, ukca_dim2,                                                &
   ukca_mmr, ukca_cvl,                                                  &
   ukca_dry, ukca_wet, ukca_rho, ukca_vol, ukca_wtv, ukca_nbr,          &
   ukcaaer_mix_ratio, ukcaaer_comp_vol,                                 &
   ukcaaer_dry_diam, ukcaaer_wet_diam, ukcaaer_modal_rho,               &
   ukcaaer_modal_vol, ukcaaer_modal_wtv, ukcaaer_modal_nbr,             &
   NPD_FIELD, NPD_PROFILE, NPD_LAYER,                                   &
   n_ukca_cpnt, n_ukca_mode                                             &
   )

USE parkind1, ONLY: jpim, jprb 
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE

!
! Arguments
!
!
! Fixed array dimensions
!
  ! Field size
  INTEGER NPD_FIELD
  
  ! Size of array of profiles
  INTEGER NPD_PROFILE
  
  ! Maximum number of layers
  INTEGER NPD_LAYER
  
  ! Total number of UKCA aerosol components
  INTEGER n_ukca_cpnt
  
  ! Total number of UKCA aerosol modes
  INTEGER n_ukca_mode  
!
! Actual array dimensions
!
  ! Number of profiles
  INTEGER n_profile
  
  ! Number of layers seen by radiation
  INTEGER n_layer
  
  ! Dimensions for input UKCA arrays
  INTEGER ukca_dim1
  INTEGER ukca_dim2

!
! With intent in
!
  ! model switch to include an extra top layer in the radiation scheme
  LOGICAL L_EXTRA_TOP
  
  ! gathering array
  INTEGER, DIMENSION(NPD_FIELD) :: i_gather
  
  ! UKCA component mass mixing ratios
  REAL, DIMENSION(ukca_dim1, ukca_dim2, n_ukca_cpnt) :: ukca_mmr
  
  ! UKCA component volumes
  REAL, DIMENSION(ukca_dim1, ukca_dim2, n_ukca_cpnt) :: ukca_cvl
  
  ! UKCA modal dry and wet diameters
  REAL, DIMENSION(ukca_dim1, ukca_dim2, n_ukca_mode) :: ukca_dry
  REAL, DIMENSION(ukca_dim1, ukca_dim2, n_ukca_mode) :: ukca_wet
  
  ! UKCA modal densities
  REAL, DIMENSION(ukca_dim1, ukca_dim2, n_ukca_mode) :: ukca_rho
  
  ! UKCA modal volumes
  REAL, DIMENSION(ukca_dim1, ukca_dim2, n_ukca_mode) :: ukca_vol
  
  ! UKCA modal volume of water
  REAL, DIMENSION(ukca_dim1, ukca_dim2, n_ukca_mode) :: ukca_wtv
  
  ! UKCA modal number concentrations
  REAL, DIMENSION(ukca_dim1, ukca_dim2, n_ukca_mode) :: ukca_nbr
  
!
! With intent out
!
  ! UKCA component mass mixing ratios on radiation code domain
  REAL, DIMENSION(NPD_PROFILE, NPD_LAYER, n_ukca_cpnt) :: &
    ukcaaer_mix_ratio
    
  ! UKCA component volumes on radiation code domain
  REAL, DIMENSION(NPD_PROFILE, NPD_LAYER, n_ukca_cpnt) :: &
    ukcaaer_comp_vol
    
  ! UKCA modal dry and wet diameters on radiation code domain
  REAL, DIMENSION(NPD_PROFILE, NPD_LAYER, n_ukca_mode) :: &
    ukcaaer_dry_diam
  REAL, DIMENSION(NPD_PROFILE, NPD_LAYER, n_ukca_mode) :: &
    ukcaaer_wet_diam
    
  ! UKCA modal densities on radiation code domain
  REAL, DIMENSION(NPD_PROFILE, NPD_LAYER, n_ukca_mode) :: &
    ukcaaer_modal_rho
  
  ! UKCA modal volumes
  REAL, DIMENSION(NPD_PROFILE, NPD_LAYER, n_ukca_mode) :: &
    ukcaaer_modal_vol
    
  ! UKCA modal volumes of water
  REAL, DIMENSION(NPD_PROFILE, NPD_LAYER, n_ukca_mode) :: &
    ukcaaer_modal_wtv
    
  ! UKCA modal number concentrations
  REAL, DIMENSION(NPD_PROFILE, NPD_LAYER, n_ukca_mode) :: &
    ukcaaer_modal_nbr
    
!
! Local variables
!
  INTEGER i, &
          j, &
          l
  INTEGER i_top_copy
  INTEGER lg

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!
!
  IF (lhook) THEN
    CALL dr_hook('UKCA_RADAER_SET_AEROSOL_FIELD', zhook_in, zhook_handle)
  END IF
  
  IF (L_EXTRA_TOP) THEN
    i_top_copy = 2
  ELSE
    i_top_copy = 1
  END IF

  DO j = 1, n_ukca_mode

    ! Here would a good place to switch off the direct effect
    ! of a given mode depending on a model switch.
  
    DO i = i_top_copy, n_layer
      DO l = 1, n_profile
        lg = i_gather(l)
        ukcaaer_dry_diam(l, i, j)  = ukca_dry(lg, n_layer+1-i, j)
        ukcaaer_wet_diam(l, i, j)  = ukca_wet(lg, n_layer+1-i, j)
        ukcaaer_modal_rho(l, i, j) = ukca_rho(lg, n_layer+1-i, j)
        ukcaaer_modal_vol(l, i, j) = ukca_vol(lg, n_layer+1-i, j)
        ukcaaer_modal_wtv(l, i, j) = ukca_wtv(lg, n_layer+1-i, j)
        ukcaaer_modal_nbr(l, i, j) = ukca_nbr(lg, n_layer+1-i, j)
      END DO ! l
    END DO ! i

    ! If using an extra top layer extrapolate the modal properties
    ! from the adjacent layer.  
    IF (L_EXTRA_TOP) THEN
      DO l = 1, n_profile
        ukcaaer_dry_diam(l, 1, j)  = ukcaaer_dry_diam(l, 2, j)
        ukcaaer_wet_diam(l, 1, j)  = ukcaaer_wet_diam(l, 2, j)
        ukcaaer_modal_rho(l, 1, j) = ukcaaer_modal_rho(l, 2, j)
        ukcaaer_modal_vol(l, 1, j) = ukcaaer_modal_vol(l, 2, j)
        ukcaaer_modal_wtv(l, 1, j) = ukcaaer_modal_wtv(l, 2, j)
        ukcaaer_modal_nbr(l, 1, j) = ukcaaer_modal_nbr(l, 2, j)
      END DO ! l
    END IF
  
  END DO ! j

  DO j = 1, n_ukca_cpnt

    ! Here would be a good place to switch off the direct effect
    ! of a given component depending on a model switch.
  
    DO i = i_top_copy, n_layer
      DO l = 1, n_profile
        lg = i_gather(l)
        ukcaaer_mix_ratio(l, i, j) = ukca_mmr(lg, n_layer+1-i, j)
        ukcaaer_comp_vol(l, i, j)  = ukca_cvl(lg, n_layer+1-i, j)
      END DO ! l
    END DO ! i
  
    ! If using an extra top layer extrapolate the mixing ratio
    ! from the adjacent layer.
    IF (L_EXTRA_TOP) THEN
      DO l = 1, n_profile
        ukcaaer_mix_ratio(l, 1, j) = ukcaaer_mix_ratio(l, 2, j)
        ukcaaer_comp_vol(l, 1, j)  = ukcaaer_comp_vol(l, 2, j)
      END DO ! l
    END IF
  
  END DO ! j

  IF (lhook) THEN
    CALL dr_hook('UKCA_RADAER_SET_AEROSOL_FIELD', zhook_out, zhook_handle)
  END IF
  
RETURN
END SUBROUTINE ukca_radaer_set_aerosol_field
