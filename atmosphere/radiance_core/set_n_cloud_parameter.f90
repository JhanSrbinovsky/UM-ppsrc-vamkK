! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Function to set number of cloudy parameters.
!
! Method:
!       Straightforward
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------
FUNCTION set_n_cloud_parameter(i_scheme, i_component, n_phase_term)


  USE realtype_rd, ONLY: RealK
  USE rad_pcf
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE


! Dummy arguments.
  INTEGER, INTENT(IN) ::                                                &
      i_scheme                                                          &
!       Parametrization scheme
    , i_component                                                       &
!       Component in cloud
    , n_phase_term
!       Number of terms in the phase function

  INTEGER ::                                                            &
      set_n_cloud_parameter
!       Returned number of coefficients in parametrization

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle


  IF (lhook) CALL dr_hook('SET_N_CLOUD_PARAMETER',zhook_in,zhook_handle)

  IF ( (i_component == ip_clcmp_st_water).OR.                           &
       (i_component == ip_clcmp_cnv_water) ) THEN

    IF (i_scheme == ip_slingo_schrecker) THEN
      set_n_cloud_parameter=6
    ELSE IF (i_scheme == ip_ackerman_stephens) THEN
      set_n_cloud_parameter=9
    ELSE IF (i_scheme == ip_drop_pade_2) THEN
      set_n_cloud_parameter=16
    ELSE IF (i_scheme == ip_slingo_schr_phf) THEN
      set_n_cloud_parameter=4+2*n_phase_term
    ELSE IF (i_scheme == IP_drop_pade_2_PHF) THEN
      set_n_cloud_parameter=11+5*n_phase_term
    ELSE IF (i_scheme == IP_ps_size_PHF) THEN
      set_n_cloud_parameter=11+5*n_phase_term
    END IF

  ELSE IF ( (i_component == ip_clcmp_st_ice).OR.                        &
            (i_component == ip_clcmp_cnv_ice) ) THEN

    IF (i_scheme == ip_slingo_schrecker_ice) THEN
      set_n_cloud_parameter=6
    ELSE IF (i_scheme == ip_ice_adt) THEN
      set_n_cloud_parameter=30
    ELSE IF (i_scheme == ip_ice_adt_10) THEN
      set_n_cloud_parameter=36
    ELSE IF (i_scheme == ip_sun_shine_vn2_vis) THEN
      set_n_cloud_parameter=6
    ELSE IF (i_scheme == ip_sun_shine_vn2_ir) THEN
      set_n_cloud_parameter=0
    ELSE IF (i_scheme == ip_ice_fu_solar) THEN
      set_n_cloud_parameter=14
    ELSE IF (i_scheme == ip_ice_fu_ir) THEN
      set_n_cloud_parameter=10
    ELSE IF (i_scheme == ip_slingo_schr_ice_phf) THEN
      set_n_cloud_parameter=4+2*n_phase_term
    ELSE IF (i_scheme == ip_ice_fu_phf) THEN
      set_n_cloud_parameter=9+5*n_phase_term
    ELSE IF (i_scheme == ip_ice_t_iwc) THEN
      set_n_cloud_parameter=9
    END IF

  END IF


  IF (lhook) CALL dr_hook('SET_N_CLOUD_PARAMETER',zhook_out,zhook_handle)
END FUNCTION set_n_cloud_parameter
