! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine to set the default values of the control structure.

! Description:
!   This module declares the controlling structure for SW radiation
!   and performs suitable default initializations.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90

SUBROUTINE sw_control_default(sw_control)

      USE control_struc
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

      TYPE (control_option) :: sw_control
!       The block of controlling options for the code

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('SW_CONTROL_DEFAULT',zhook_in,zhook_handle)

      sw_control%spectral_file=''

!     Range of bands
      sw_control%first_band=1
      sw_control%last_band=1

!     Miscellaneous options
      sw_control%isolir=1
      sw_control%i_solar_src=3
!       Index of solar source function used.
!       i_solar_src = 1  Original Kurucz function
!       i_solar_src = 2  AER version of Kurucz
!       i_solar_src = 3  Kurucz function used in UK model
!       i_solar_src = 4  UK reduced version
!       i_solar_src = 5  Kurucz modtran
!       i_solar_src = 6  Labs Neckel
!                             Z. Sun
      sw_control%i_scatter_method=1
      sw_control%l_ir_source_quad=.FALSE.
      sw_control%l_extra_top=.FALSE.
      sw_control%l_rad_deg=.FALSE.
      sw_control%l_subsample=.FALSE.

!     Properties of clouds
      sw_control%i_cloud=6
      sw_control%i_cloud_representation=4
      sw_control%i_st_water=0
      sw_control%i_cnv_water=0
      sw_control%i_st_ice=0
      sw_control%i_cnv_ice=0
      sw_control%l_local_cnv_partition=.TRUE.
      sw_control%l_global_cloud_top=.TRUE.
      sw_control%i_fsd=0
      sw_control%i_inhom=0
      sw_control%i_overlap=0

!     Physical processes
      sw_control%l_microphysics=.TRUE.
      sw_control%l_gas=.TRUE.
      sw_control%l_rayleigh=.TRUE.
      sw_control%l_continuum=.TRUE.
      sw_control%l_cloud=.TRUE.
      sw_control%l_drop=.TRUE.
      sw_control%l_ice=.TRUE.
      sw_control%l_aerosol=.TRUE.
      sw_control%l_aerosol_ccn=.TRUE.

!     Gaseous absorption
      sw_control%i_gas_overlap=5
      sw_control%l_o2=.FALSE.
      sw_control%l_n2o=.FALSE.
      sw_control%l_ch4=.FALSE.
      sw_control%l_cfc11=.FALSE.
      sw_control%l_cfc12=.FALSE.
      sw_control%l_cfc113=.FALSE.
      sw_control%l_cfc114=.FALSE.
      sw_control%l_hcfc22=.FALSE.
      sw_control%l_hfc125=.FALSE.
      sw_control%l_hfc134a=.FALSE.

!     Angular integration (including algorithmic options)
      sw_control%i_angular_integration=1
      sw_control%i_2stream=16
      sw_control%i_solver=14
      sw_control%n_order_gauss=0
      sw_control%i_truncation=1
      sw_control%i_sph_algorithm=1
      sw_control%n_order_phase_solar=1
      sw_control%ls_global_trunc=9
      sw_control%ms_min=0
      sw_control%ms_max=0
      sw_control%ls_brdf_trunc=0
      sw_control%accuracy_adaptive=1.0e-04
      sw_control%l_rescale=.TRUE.
      sw_control%l_henyey_greenstein_pf=.TRUE.
      sw_control%i_sph_mode=1
      sw_control%l_euler_trnf=.FALSE.

!     Satellite data
      sw_control%l_geostationary=.FALSE.
      sw_control%sat_desc="                    "  // &
                          "                    "  // &
                          "                    "  // &
                          "                    "
      sw_control%sat_hgt=0.0
      sw_control%sat_lon=0.0
      sw_control%sat_lat=0.0
      sw_control%max_view_lon=0.0
      sw_control%min_view_lon=0.0
      sw_control%max_view_lat=0.0
      sw_control%min_view_lat=0.0

      sw_control%l_layer=.TRUE.
      sw_control%l_cloud_layer=.TRUE.
      sw_control%l_2_stream_correct=.FALSE.

      IF (lhook) CALL dr_hook('SW_CONTROL_DEFAULT',zhook_out,zhook_handle)

END SUBROUTINE sw_control_default
