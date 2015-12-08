! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to find single scattering properties of all regions.
!
! Method:
!   Straightforward.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiance Core
!
!- ---------------------------------------------------------------------
SUBROUTINE single_scattering_all(i_scatter_method_band                  &
!                 Atmospheric Properties
    , n_profile, n_layer, d_mass                                        &
!                 Cloudy Properties
    , l_cloud, n_cloud_top, n_cloud_type                                &
!                 Optical Properties
    , ss_prop, k_gas_abs                                                &
!                 Dimensions of Arrays
    , nd_profile, nd_layer, nd_layer_clr, id_ct, nd_cloud_type          &
    )


  USE realtype_rd, ONLY: RealK
  USE def_ss_prop
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Sizes of arrays.
  INTEGER, INTENT(IN) ::                                                &
      nd_profile                                                        &
!       Size allocated for profiles
    , nd_layer                                                          &
!       Size allocated for layers
    , nd_layer_clr                                                      &
!       Size allocated for completely clear layers
    , id_ct                                                             &
!       Topmost declared cloudy layer
    , nd_cloud_type
!       Size allocated for types of clouds


! Dummy variables.
  INTEGER, INTENT(IN) ::                                                &
      i_scatter_method_band
!       Treatment of scattering in the band

!                 Atmospheric properties
  INTEGER, INTENT(IN) ::                                                &
      n_profile                                                         &
!       Number of profiles
    , n_layer
!       Number of layers
  REAL (RealK), INTENT(IN) ::                                           &
      d_mass(nd_profile, nd_layer)
!       Mass thickness of each layer

!                 Cloudy properties
  LOGICAL, INTENT(IN) ::                                                &
      l_cloud
!       Flag for clouds
  INTEGER, INTENT(IN) ::                                                &
      n_cloud_top                                                       &
!       Topmost cloudy layer
    , n_cloud_type
!       Number of types of clouds

!                 Optical properties
  TYPE(STR_ss_prop), INTENT(INOUT) :: ss_prop
!       Single scattering properties of the atmosphere
  REAL (RealK), INTENT(IN) ::                                           &
      k_gas_abs(nd_profile, nd_layer)
!       Gaseous extinction


! Local variables.
  INTEGER                                                               &
      k
!       Loop variable

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('SINGLE_SCATTERING_ALL',zhook_in,zhook_handle)

! Clear-sky properties:

! In the following call K_GAS_ABS can be used as if it had the
! smaller dimension ND_LAYER_CLR as long as the last dimension
! is over atmospheric layers.
! DEPENDS ON: single_scattering
  CALL single_scattering(i_scatter_method_band                          &
    , n_profile, 1, n_cloud_top-1                                       &
    , d_mass                                                            &
    , ss_prop%k_grey_tot_clr, ss_prop%k_ext_scat_clr, k_gas_abs         &
    , ss_prop%tau_clr, ss_prop%omega_clr                                &
    , nd_profile, nd_layer, 1, nd_layer_clr                             &
    )
  CALL single_scattering(i_scatter_method_band                          &
    , n_profile, n_cloud_top, n_layer                                   &
    , d_mass                                                            &
    , ss_prop%k_grey_tot(1, id_ct, 0)                                   &
    , ss_prop%k_ext_scat(1, id_ct, 0)                                   &
    , k_gas_abs                                                         &
    , ss_prop%tau(1, id_ct, 0), ss_prop%omega(1, id_ct, 0)              &
    , nd_profile, nd_layer, id_ct, nd_layer                             &
    )

  IF (l_cloud) THEN
    DO k=1, n_cloud_type
      CALL single_scattering(i_scatter_method_band                      &
        , n_profile, n_cloud_top, n_layer                               &
        , d_mass                                                        &
        , ss_prop%k_grey_tot(1, id_ct, k)                               &
        , ss_prop%k_ext_scat(1, id_ct, k)                               &
        , k_gas_abs                                                     &
        , ss_prop%tau(1, id_ct, k), ss_prop%omega(1, id_ct, k)          &
        , nd_profile, nd_layer, id_ct, nd_layer                         &
        )
    END DO
  END IF


  IF (lhook) CALL dr_hook('SINGLE_SCATTERING_ALL',zhook_out,zhook_handle)

END SUBROUTINE single_scattering_all
