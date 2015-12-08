! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!+ Re-arrange UKCA-MODE input to match the expectations of
!+ routine ukca_radaer_band_average().
!
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
SUBROUTINE ukca_radaer_prepare(                                         &
      ! Actual array dimensions
      n_profile, n_layer, n_ukca_mode, n_ukca_cpnt,                     &
      ! UKCA_RADAER structure
      ukca_radaer,                                                      &
      ! Component mass-mixing ratios
      ukca_mix_ratio,                                                   &
      ! Modal mass-mixing ratios 
      ukca_modal_mixr,                                                  &
      ! Input modal number concentrations
      ukca_modal_nbr,                                                   &
      ! Output modal number concentrations
      ukca_modal_number,                                                &
      ! Pressure and temperature
      pressure, temperature,                                            &
      ! Fixed array dimensions
      npd_profile, npd_layer                                            &
  )

USE parkind1, ONLY: jpim, jprb 
USE yomhook,  ONLY: lhook, dr_hook
USE ukca_radaer_struct_mod

IMPLICIT NONE

!
! Arguments with intent(in)
!
!
! Fixed array dimensions
!
INTEGER :: npd_profile, &
           npd_layer
!
! Actual array dimensions
!
INTEGER :: n_profile,   &
           n_layer,     &
           n_ukca_mode, &
           n_ukca_cpnt
!
! Structure for UKCA/radiation interaction
!
TYPE (ukca_radaer_struct) :: ukca_radaer

!
! Component mass-mixing ratios
!
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_cpnt) :: ukca_mix_ratio

!
! Modal number concentrations divided by molecular concentration of air
!
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: ukca_modal_nbr

!
! Pressure and temperature fields.
REAL, DIMENSION(npd_profile, npd_layer) :: pressure, &
                                           temperature
!
!
! Arguments with intent(out)
!
!
!
! Modal mass-mixing ratios
!
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: ukca_modal_mixr

!
! Modal number concentrations (in m-3)
!
REAL, DIMENSION(npd_profile, npd_layer, n_ukca_mode) :: ukca_modal_number

!
! Local variables
!
INTEGER i, &
        j, &
        k, &
        l
INTEGER this_cpnt

!
! Boltzmann constant
!
REAL, PARAMETER :: k_boltzmann = 1.3807E-23 ! J/K

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('UKCA_RADAER_PREPARE', zhook_in, zhook_handle)

!
! Modal mass-mixing ratios.
!
! Simply sum the mixing ratios of all components included in a 
! given mode to get the mixing ratio for that mode.
!
DO j = 1, n_ukca_mode

  DO k = 1, n_layer
  
    DO l = 1, n_profile
    
      ukca_modal_mixr(l, k, j) = 0.0
    
    END DO ! l
  
  END DO ! k
  
  DO i = 1, ukca_radaer%n_cpnt_in_mode(j)
  
    this_cpnt = ukca_radaer%i_cpnt_index(i, j)
    
    DO k = 1, n_layer
    
      DO l = 1, n_profile
      
        ukca_modal_mixr(l, k, j) = &
          ukca_modal_mixr(l, k, j) + ukca_mix_ratio(l, k, this_cpnt)
      
      END DO ! l
    
    END DO ! k
  
  END DO ! i
  
END DO ! j

!
! Modal number concentrations
!
! Multiply by the molecular concentration of air (p/kT) to obtain
! the acutal aerosol number concentrations.
!
DO j = 1, n_ukca_mode

  DO k = 1, n_layer
  
    DO l = 1, n_profile
      
      ukca_modal_number(l, k, j) = ukca_modal_nbr(l, k, j) * &
                                   pressure(l, k) / &
                                   (k_boltzmann * temperature(l, k))
      
    END DO ! l
  
  END DO ! k

END DO ! j

IF (lhook) CALL dr_hook('UKCA_RADAER_PREPARE', zhook_out, zhook_handle)

RETURN
END SUBROUTINE ukca_radaer_prepare
