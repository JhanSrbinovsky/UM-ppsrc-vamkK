! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!+ Returns basic information on the spectral files loaded in memory.
!
!
! Subroutine Interface:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: UKCA
!
      SUBROUTINE ukca_get_specinfo(                                     &
             nbr_band_sw,                                               &
             nbr_band_lw,                                               &
             nbr_aod_wavel                                              &
        )

      USE parkind1, ONLY: jpim, jprb 
      USE yomhook,  ONLY: lhook, dr_hook
      USE spec_sw_lw

      IMPLICIT NONE
     
!     Arguments with intent(out)
      INTEGER nbr_band_sw
      INTEGER nbr_band_lw
      INTEGER nbr_aod_wavel

!     Local variables
      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) THEN
        CALL dr_hook('UKCA_RADAER_GET_SPECINFO', zhook_in, zhook_handle)
      END IF
      
!     Copy the required info

      nbr_band_sw   = sw_spectrum(1)%n_band
      nbr_band_lw   = lw_spectrum(1)%n_band
      nbr_aod_wavel = lw_spectrum(1)%n_aod_wavel

      IF (lhook) THEN
        CALL dr_hook('UKCA_RADAER_GET_SPECINFO', zhook_out, zhook_handle)
      END IF
      
      RETURN
      END SUBROUTINE ukca_get_specinfo
