! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for lw radiation.

! Description:
!   Module containing control as used by the lw radiation code.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP

MODULE lw_rad_input_mod

USE sat_opt_mod
USE lw_control_struct, ONLY: lw_control, n_lwcall
USE missing_data_mod, ONLY: imdi 

IMPLICIT NONE

CHARACTER  (LEN=132) :: spectral_file_lw    !       Spectral file

!     Options for trace gases:
LOGICAL :: l_n2o_lw            ! flag for nitrous oxide 
LOGICAL :: l_ch4_lw            ! flag for methane
LOGICAL :: l_cfc11_lw          ! flag for cfc11 
LOGICAL :: l_cfc12_lw          ! flag for cfc12
LOGICAL :: l_cfc113_lw         ! flag for cfc113 
LOGICAL :: l_cfc114_lw         ! flag for cfc114
LOGICAL :: l_hcfc22_lw         ! flag for hcfc22 
LOGICAL :: l_hfc125_lw         ! flag for hfc125 
LOGICAL :: l_hfc134a_lw        ! flag for hfc134a 

LOGICAL :: l_solar_tail_flux   ! Flag for adding solar tail flux to LW 
LOGICAL :: l_microphysics_lw   ! Flag for microphysics in LW 

INTEGER :: i_scatter_method_lw ! treatment of scattering
INTEGER :: i_gas_overlap_lw    ! treatment of gaseous overlaps 

! types of droplets or ice crystals used for parametrisations
INTEGER :: i_st_water_lw             ! type for stratiform water 
INTEGER :: i_cnv_water_lw            ! type for convective water 
INTEGER :: i_st_ice_lw               ! type for stratiform ice  
INTEGER :: i_cnv_ice_lw              ! type for convective ice 

! Second set of control options for the second radiation call: 
CHARACTER  (LEN=132) :: spectral_file_lw2 
LOGICAL :: l_n2o_lw2
LOGICAL :: l_ch4_lw2 
LOGICAL :: l_cfc11_lw2
LOGICAL :: l_cfc12_lw2 
LOGICAL :: l_cfc113_lw2
LOGICAL :: l_cfc114_lw2 
LOGICAL :: l_hcfc22_lw2
LOGICAL :: l_hfc125_lw2
LOGICAL :: l_hfc134a_lw2
LOGICAL :: l_solar_tail_flux_2
LOGICAL :: l_microphysics_lw2 
INTEGER :: i_scatter_method_lw2
INTEGER :: i_gas_overlap_lw2 
INTEGER :: i_st_water_lw2 
INTEGER :: i_cnv_water_lw2 
INTEGER :: i_st_ice_lw2
INTEGER :: i_cnv_ice_lw2

! --------------------

NAMELIST/r2lwncal/n_lwcall

NAMELIST/r2lwclnl/                                                      &
     spectral_file_lw, spectral_file_lw2,                               &
     l_solar_tail_flux, l_solar_tail_flux_2,                            &
     i_gas_overlap_lw, i_gas_overlap_lw2,                               &
     l_n2o_lw, l_n2o_lw2, l_ch4_lw, l_ch4_lw2,                          &
     l_cfc11_lw, l_cfc11_lw2, l_cfc12_lw, l_cfc12_lw2,                  &
     l_cfc113_lw, l_cfc113_lw2, l_cfc114_lw, l_cfc114_lw2,              &
     l_hcfc22_lw, l_hcfc22_lw2, l_hfc125_lw, l_hfc125_lw2,              &
     l_hfc134a_lw, l_hfc134a_lw2,                                       &
     i_st_water_lw, i_st_water_lw2, i_cnv_water_lw, i_cnv_water_lw2,    &
     i_st_ice_lw, i_st_ice_lw2, i_cnv_ice_lw, i_cnv_ice_lw2,            &
     i_scatter_method_lw, i_scatter_method_lw2,                         &
     l_microphysics_lw, l_microphysics_lw2


! ----------------------------------------------------------------
CONTAINS

! Subroutine to set the input values of the lw control structure.

SUBROUTINE lw_input
USE check_iostat_mod
USE rad_pcf, ONLY: ip_cloud_mix_max, ip_cloud_mix_random,               &
  ip_cloud_triple, ip_cloud_part_corr, ip_cloud_part_corr_cnv,          &
  ip_cloud_mcica,                                                       &
  ip_solver_homogen_direct, ip_solver_mix_direct_hogan,                 &
  ip_solver_triple_hogan, ip_solver_triple_app_scat,                    &
  ip_solver_mix_app_scat,                                               &
  ip_cloud_homogen, ip_cloud_ice_water, ip_cloud_conv_strat,            &
  ip_cloud_csiw,                                                        &
  ip_scaling, ip_mcica,                                                 &
  ip_max_rand, ip_exponential_rand, ip_rand,                            &
  ip_scatter_approx
USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE ereport_mod, ONLY : ereport
USE rad_input_mod

IMPLICIT NONE

INTEGER :: j,                        &
     errorstatus      ! Return code : 0 Normal Exit : >0 Error

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*) RoutineName
PARAMETER (   RoutineName='sw_rad_input_mod')
CHARACTER(LEN=256) :: CMessage         ! Error message if Errorstatus >0

IF (lhook) CALL dr_hook('lw_input',zhook_in,zhook_handle)

! Initialisation: set namelist items to MDI / FALSE.
  spectral_file_lw           = ''
  l_solar_tail_flux          = .FALSE.
  i_gas_overlap_lw           = imdi
  l_n2o_lw                   = .FALSE.
  l_ch4_lw                   = .FALSE.
  l_cfc11_lw                 = .FALSE.
  l_cfc12_lw                 = .FALSE.
  l_cfc113_lw                = .FALSE.
  l_cfc114_lw                = .FALSE.
  l_hcfc22_lw                = .FALSE.
  l_hfc125_lw                = .FALSE.
  l_hfc134a_lw               = .FALSE.
  i_st_water_lw              = imdi
  i_cnv_water_lw             = imdi
  i_st_ice_lw                = imdi
  i_cnv_ice_lw               = imdi
  i_scatter_method_lw        = imdi
  l_microphysics_lw          = .FALSE.

  spectral_file_lw2          = ''
  l_solar_tail_flux_2        = .FALSE.
  i_gas_overlap_lw2          = imdi
  l_n2o_lw2                  = .FALSE.
  l_ch4_lw2                  = .FALSE.
  l_cfc11_lw2                = .FALSE.
  l_cfc12_lw2                = .FALSE.
  l_cfc113_lw2               = .FALSE.
  l_cfc114_lw2               = .FALSE.
  l_hcfc22_lw2               = .FALSE.
  l_hfc125_lw2               = .FALSE.
  l_hfc134a_lw2              = .FALSE.
  i_st_water_lw2             = imdi
  i_cnv_water_lw2            = imdi
  i_st_ice_lw2               = imdi
  i_cnv_ice_lw2              = imdi
  i_scatter_method_lw2       = imdi
  l_microphysics_lw2         = .FALSE.

! Options for the longwave



  READ (UNIT=5, NML= r2lwclnl, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist R2LWCLNL")


DO j=1, n_lwcall
! Set default values of control variables.
! Note: to allow a simple interface for ROSE only the commonly used
! items are now set via the namelist. In order to set further variables
! (eg. for radiance diagnostics) a branch must now be used.

! DEPENDS ON: lw_control_default
  CALL lw_control_default(lw_control(j))

! Transfer namelist items to the data structure.
  IF (j==2) THEN
    lw_control(j)%spectral_file          = spectral_file_lw2
    lw_control(j)%l_solar_tail_flux      = l_solar_tail_flux_2
    lw_control(j)%i_gas_overlap          = i_gas_overlap_lw2
    lw_control(j)%i_cloud_representation = i_cloud_representation_2
    lw_control(j)%i_fsd                  = i_fsd_2
    lw_control(j)%i_overlap              = i_overlap_2
    lw_control(j)%i_inhom                = i_inhom_2
    lw_control(j)%l_n2o                  = l_n2o_lw2
    lw_control(j)%l_ch4                  = l_ch4_lw2
    lw_control(j)%l_cfc11                = l_cfc11_lw2
    lw_control(j)%l_cfc12                = l_cfc12_lw2
    lw_control(j)%l_cfc113               = l_cfc113_lw2
    lw_control(j)%l_cfc114               = l_cfc114_lw2
    lw_control(j)%l_hcfc22               = l_hcfc22_lw2
    lw_control(j)%l_hfc125               = l_hfc125_lw2
    lw_control(j)%l_hfc134a              = l_hfc134a_lw2
    lw_control(j)%i_st_water             = i_st_water_lw2
    lw_control(j)%i_cnv_water            = i_cnv_water_lw2
    lw_control(j)%i_st_ice               = i_st_ice_lw2
    lw_control(j)%i_cnv_ice              = i_cnv_ice_lw2
    lw_control(j)%i_scatter_method       = i_scatter_method_lw2
    lw_control(j)%l_microphysics         = l_microphysics_lw2
  ELSE
    lw_control(j)%spectral_file          = spectral_file_lw
    lw_control(j)%l_solar_tail_flux      = l_solar_tail_flux
    lw_control(j)%i_gas_overlap          = i_gas_overlap_lw
    lw_control(j)%i_cloud_representation = i_cloud_representation
    lw_control(j)%i_fsd                  = i_fsd
    lw_control(j)%i_overlap              = i_overlap
    lw_control(j)%i_inhom                = i_inhom
    lw_control(j)%l_n2o                  = l_n2o_lw
    lw_control(j)%l_ch4                  = l_ch4_lw
    lw_control(j)%l_cfc11                = l_cfc11_lw
    lw_control(j)%l_cfc12                = l_cfc12_lw
    lw_control(j)%l_cfc113               = l_cfc113_lw
    lw_control(j)%l_cfc114               = l_cfc114_lw
    lw_control(j)%l_hcfc22               = l_hcfc22_lw
    lw_control(j)%l_hfc125               = l_hfc125_lw
    lw_control(j)%l_hfc134a              = l_hfc134a_lw
    lw_control(j)%i_st_water             = i_st_water_lw
    lw_control(j)%i_cnv_water            = i_cnv_water_lw
    lw_control(j)%i_st_ice               = i_st_ice_lw
    lw_control(j)%i_cnv_ice              = i_cnv_ice_lw
    lw_control(j)%i_scatter_method       = i_scatter_method_lw
    lw_control(j)%l_microphysics         = l_microphysics_lw
  END IF
END DO

! Set i_cloud from combination of i_overlap, i_fsd and 
! i_representation
DO j=1,n_lwcall
  IF (lw_control(j)%i_cloud_representation                          &
        == ip_cloud_homogen) THEN
    lw_control(j)%i_solver=ip_solver_homogen_direct
    IF (lw_control(j)%i_inhom == ip_mcica) THEN
      lw_control(j)%i_cloud=ip_cloud_mcica
    ELSE
      IF (lw_control(j)%i_overlap == ip_max_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_mix_max
      ELSE IF (lw_control(j)%i_overlap == ip_exponential_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_part_corr
      ELSE IF (lw_control(j)%i_overlap == ip_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage) 
      END IF
    END IF  

  ELSE IF (lw_control(j)%i_cloud_representation                      &
       == ip_cloud_ice_water) THEN
    IF (lw_control(j)%i_inhom == ip_mcica) THEN
      lw_control(j)%i_cloud=ip_cloud_mcica
      lw_control(j)%i_solver=ip_solver_homogen_direct
    ELSE
      IF (lw_control(j)%i_scatter_method == ip_scatter_approx) THEN
        lw_control(j)%i_solver=ip_solver_mix_app_scat
      ELSE
        lw_control(j)%i_solver=ip_solver_mix_direct_hogan
      END IF
      IF (lw_control(j)%i_overlap == ip_max_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_mix_max
      ELSE IF (lw_control(j)%i_overlap == ip_exponential_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_part_corr
      ELSE IF (lw_control(j)%i_overlap == ip_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage) 
      END IF
    END IF  

  ELSE IF ((lw_control(j)%i_cloud_representation                     &
           == ip_cloud_conv_strat) .OR.                              &
           (lw_control(j)%i_cloud_representation                     &
           == ip_cloud_csiw)) THEN
    IF (lw_control(j)%i_inhom == ip_mcica) THEN
      ErrorStatus = 100
      CMessage = 'McICA is not compatible with the selected'//       &
                 ' cloud representation'
      CALL ereport(RoutineName, ErrorStatus, CMessage) 
    ELSE
      IF (lw_control(j)%i_scatter_method == ip_scatter_approx) THEN
        lw_control(j)%i_solver=ip_solver_triple_app_scat
      ELSE
        lw_control(j)%i_solver=ip_solver_triple_hogan
      END IF
       IF (lw_control(j)%i_overlap == ip_max_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_triple
      ELSE IF (lw_control(j)%i_overlap == ip_exponential_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_part_corr_cnv
      ELSE IF (lw_control(j)%i_overlap == ip_rand) THEN
        lw_control(j)%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage) 
      END IF
    END IF
  END IF

  IF (lw_control(j)%i_inhom == ip_scaling) THEN
    l_inhom_cloud=.TRUE.
  END IF

END DO

IF (lhook) CALL dr_hook('lw_input',zhook_out,zhook_handle)

END SUBROUTINE lw_input

END MODULE lw_rad_input_mod
