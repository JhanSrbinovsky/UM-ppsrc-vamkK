! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine to set the default values of the control structure.

! Description:
!   This module defines the controlling structure for LW calculations
!   and performs suitable default initializations.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90

SUBROUTINE lw_control_default(lw_control)

      USE control_struc
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

      TYPE (control_option) :: lw_control
!       The block of controlling options for the code

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('LW_CONTROL_DEFAULT',zhook_in,zhook_handle)

      lw_control%spectral_file=''

!     Default setting of range of bands
      lw_control%first_band=1
      lw_control%last_band=1

!     Miscellaneous options
      lw_control%isolir=2
      lw_control%i_solar_src=3
!       Index of solar source function used (the tail may be used in LW).
!       i_solar_src = 1  Original Kurucz function
!       i_solar_src = 2  AER version of Kurucz
!       i_solar_src = 3  Kurucz function used in UK model
!       i_solar_src = 4  UK reduced version
!       i_solar_src = 5  Kurucz modtran
!       i_solar_src = 6  Labs Neckel
!                             Z. Sun
      lw_control%i_scatter_method=1
      lw_control%l_ir_source_quad=.TRUE.
      lw_control%l_extra_top=.FALSE.
      lw_control%l_rad_deg=.FALSE.
      lw_control%l_subsample=.FALSE.

!     Properties of clouds
      lw_control%i_cloud=6
      lw_control%i_cloud_representation=4
      lw_control%i_st_water=0
      lw_control%i_cnv_water=0
      lw_control%i_st_ice=0
      lw_control%i_cnv_ice=0
      lw_control%l_local_cnv_partition=.TRUE.
      lw_control%l_global_cloud_top=.TRUE.
      lw_control%i_fsd=0
      lw_control%i_inhom=0
      lw_control%i_overlap=0

!     Physical processes
      lw_control%l_microphysics=.TRUE.
      lw_control%l_gas=.TRUE.
      lw_control%l_rayleigh=.FALSE.
      lw_control%l_continuum=.TRUE.
      lw_control%l_cloud=.TRUE.
      lw_control%l_drop=.TRUE.
      lw_control%l_ice=.TRUE.
      lw_control%l_aerosol=.TRUE.
      lw_control%l_aerosol_ccn=.TRUE.
      lw_control%l_solar_tail_flux=.FALSE.

!     Gaseous absorption
      lw_control%i_gas_overlap=5
      lw_control%l_o2=.FALSE.
      lw_control%l_n2o=.FALSE.
      lw_control%l_ch4=.FALSE.
      lw_control%l_cfc11=.FALSE.
      lw_control%l_cfc12=.FALSE.
      lw_control%l_cfc113=.FALSE.
      lw_control%l_cfc114=.FALSE.
      lw_control%l_hcfc22=.FALSE.
      lw_control%l_hfc125=.FALSE.
      lw_control%l_hfc134a=.FALSE.

!     Angular integration (including algorithmic options)
      lw_control%i_angular_integration=1
      lw_control%i_2stream=12
      lw_control%i_solver=14
      lw_control%n_order_gauss=0
      lw_control%i_truncation=3
      lw_control%i_sph_algorithm=1
      lw_control%n_order_phase_solar=1
      lw_control%ls_global_trunc=9
      lw_control%ms_min=0
      lw_control%ms_max=0
      lw_control%ls_brdf_trunc=0
      lw_control%accuracy_adaptive=1.0e-04
      lw_control%l_rescale=.TRUE.
      lw_control%l_henyey_greenstein_pf=.TRUE.
      lw_control%i_sph_mode=1
      lw_control%l_euler_trnf=.FALSE.

!     Satellite data
      lw_control%l_geostationary=.FALSE.
      lw_control%sat_desc="                    "  // &
                          "                    "  // &
                          "                    "  // &
                          "                    "
      lw_control%sat_hgt=0.0
      lw_control%sat_lon=0.0
      lw_control%sat_lat=0.0
      lw_control%max_view_lon=0.0
      lw_control%min_view_lon=0.0
      lw_control%max_view_lat=0.0
      lw_control%min_view_lat=0.0

      lw_control%l_layer=.TRUE.
      lw_control%l_cloud_layer=.TRUE.
      lw_control%l_2_stream_correct=.FALSE.

      IF (lhook) CALL dr_hook('LW_CONTROL_DEFAULT',zhook_out,zhook_handle)

END SUBROUTINE lw_control_default
