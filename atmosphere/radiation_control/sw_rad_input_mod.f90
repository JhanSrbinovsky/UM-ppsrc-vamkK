! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for sw radiation.

! Description:
!   Module containing control as used by the sw radiation code.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP

MODULE sw_rad_input_mod

USE sat_opt_mod
USE sw_control_struct, ONLY: sw_control, n_swcall
USE missing_data_mod, ONLY: imdi 

IMPLICIT NONE

CHARACTER (LEN=132) :: spectral_file_sw     !       Spectral file

LOGICAL :: l_o2_sw                   ! flag for oxygen
LOGICAL :: l_ch4_sw                  ! flag for methen
LOGICAL :: l_n2o_sw                  ! flag for N2O

INTEGER :: i_gas_overlap_sw          ! treatment of gaseous overlaps

! types of droplets or ice crystals used for parametrizations
INTEGER :: i_st_water_sw             ! type for stratiform water
INTEGER :: i_cnv_water_sw            ! type for convective water
INTEGER :: i_st_ice_sw               ! type for stratiform ice
INTEGER :: i_cnv_ice_sw              ! type for convective ice

! Second set of control options for the second radiation call: 
CHARACTER (LEN=132) :: spectral_file_sw2 
LOGICAL :: l_o2_sw2 
LOGICAL :: l_ch4_sw2 
LOGICAL :: l_n2o_sw2 
INTEGER :: i_gas_overlap_sw2 
INTEGER :: i_st_water_sw2 
INTEGER :: i_cnv_water_sw2 
INTEGER :: i_st_ice_sw2 
INTEGER :: i_cnv_ice_sw2 

! -----------------------
!     Elements of namelists for options in the radiation code.
NAMELIST/r2swncal/n_swcall

NAMELIST/r2swclnl/                                                      &
     spectral_file_sw, spectral_file_sw2,                               &
     i_gas_overlap_sw, i_gas_overlap_sw2,                               &
     l_o2_sw, l_o2_sw2, l_ch4_sw, l_ch4_sw2, l_n2o_sw, l_n2o_sw2,       &
     i_st_water_sw, i_st_water_sw2, i_cnv_water_sw, i_cnv_water_sw2,    &
     i_st_ice_sw, i_st_ice_sw2, i_cnv_ice_sw, i_cnv_ice_sw2



! ----------------------------------------------------------------
CONTAINS

! Subroutine to set the input values of the sw control structure

SUBROUTINE sw_input

USE check_iostat_mod

USE rad_pcf, ONLY: ip_cloud_mix_max, ip_cloud_mix_random,               &
  ip_cloud_triple, ip_cloud_part_corr, ip_cloud_part_corr_cnv,          &
  ip_cloud_mcica,                                                       &
  ip_solver_homogen_direct, ip_solver_mix_direct_hogan,                 &
  ip_solver_triple_hogan,                                               &
  ip_cloud_homogen, ip_cloud_ice_water, ip_cloud_conv_strat,            &
  ip_cloud_csiw,                                                        &
  ip_scaling, ip_mcica,                                                 &
  ip_max_rand, ip_exponential_rand, ip_rand
USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook
USE ereport_mod, ONLY : ereport
USE rad_input_mod

IMPLICIT NONE

INTEGER :: j,                                                     &
           errorstatus      ! Return code : 0 Normal Exit : >0 Error

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*) RoutineName
PARAMETER (   RoutineName='sw_rad_input_mod')
CHARACTER(LEN=256) :: CMessage         ! Error message if Errorstatus >0

IF (lhook) CALL dr_hook('sw_input',zhook_in,zhook_handle)

! Initialisation: set namelist items to MDI / FALSE.
  spectral_file_sw           = ''
  i_gas_overlap_sw           = imdi
  l_o2_sw                    = .FALSE.
  l_n2o_sw                   = .FALSE.
  l_ch4_sw                   = .FALSE.
  i_st_water_sw              = imdi
  i_cnv_water_sw             = imdi
  i_st_ice_sw                = imdi
  i_cnv_ice_sw               = imdi

  spectral_file_sw2          = ''
  i_gas_overlap_sw2          = imdi
  l_o2_sw2                   = .FALSE.
  l_n2o_sw2                  = .FALSE.
  l_ch4_sw2                  = .FALSE.
  i_st_water_sw2             = imdi
  i_cnv_water_sw2            = imdi
  i_st_ice_sw2               = imdi
  i_cnv_ice_sw2              = imdi

! Options for the shortwave



  READ (UNIT=5, NML=r2swclnl, IOSTAT=errorstatus)
  CALL check_iostat(errorstatus, "namelist R2SWCLNL")


DO j=1,n_swcall
! Set default values of control variables.
! Note: to allow a simple interface for ROSE only the commonly used
! items are now set via the namelist. In order to set further variables
! (eg. for radiance diagnostics) a branch must now be used.

! DEPENDS ON: sw_control_default
  CALL sw_control_default(sw_control(j))

! Transfer namelist items to the data structure.
  IF (j==2) THEN
    sw_control(j)%spectral_file          = spectral_file_sw2
    sw_control(j)%i_gas_overlap          = i_gas_overlap_sw2
    sw_control(j)%l_o2                   = l_o2_sw2
    sw_control(j)%l_n2o                  = l_n2o_sw2
    sw_control(j)%l_ch4                  = l_ch4_sw2
    sw_control(j)%i_st_water             = i_st_water_sw2
    sw_control(j)%i_cnv_water            = i_cnv_water_sw2
    sw_control(j)%i_st_ice               = i_st_ice_sw2
    sw_control(j)%i_cnv_ice              = i_cnv_ice_sw2
    sw_control(j)%i_cloud_representation = i_cloud_representation_2
    sw_control(j)%i_overlap              = i_overlap_2
    sw_control(j)%i_inhom                = i_inhom_2
    sw_control(j)%i_fsd                  = i_fsd_2
  ELSE
    sw_control(j)%spectral_file          = spectral_file_sw
    sw_control(j)%i_gas_overlap          = i_gas_overlap_sw
    sw_control(j)%l_o2                   = l_o2_sw
    sw_control(j)%l_n2o                  = l_n2o_sw
    sw_control(j)%l_ch4                  = l_ch4_sw
    sw_control(j)%i_st_water             = i_st_water_sw
    sw_control(j)%i_cnv_water            = i_cnv_water_sw
    sw_control(j)%i_st_ice               = i_st_ice_sw
    sw_control(j)%i_cnv_ice              = i_cnv_ice_sw
    sw_control(j)%i_cloud_representation = i_cloud_representation
    sw_control(j)%i_overlap              = i_overlap
    sw_control(j)%i_inhom                = i_inhom
    sw_control(j)%i_fsd                  = i_fsd
  END IF
END DO

! Set i_cloud from combination of i_overlap, i_fsd and 
! i_representation
DO j=1,n_swcall
  IF (sw_control(j)%i_cloud_representation                          &
        == ip_cloud_homogen) THEN

    sw_control(j)%i_solver=ip_solver_homogen_direct
    IF (sw_control(j)%i_inhom == ip_mcica) THEN
      sw_control(j)%i_cloud=ip_cloud_mcica
    ELSE
      IF (sw_control(j)%i_overlap == ip_max_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_mix_max
      ELSE IF (sw_control(j)%i_overlap == ip_exponential_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_part_corr
      ELSE IF (sw_control(j)%i_overlap == ip_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage) 
      END IF
    END IF  

  ELSE IF (sw_control(j)%i_cloud_representation                      &
       == ip_cloud_ice_water) THEN
    IF (sw_control(j)%i_inhom == ip_mcica) THEN
      sw_control(j)%i_cloud=ip_cloud_mcica
      sw_control(j)%i_solver=ip_solver_homogen_direct
    ELSE
      sw_control(j)%i_solver=ip_solver_mix_direct_hogan
      IF (sw_control(j)%i_overlap == ip_max_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_mix_max
      ELSE IF (sw_control(j)%i_overlap == ip_exponential_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_part_corr
      ELSE IF (sw_control(j)%i_overlap == ip_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage) 
      END IF
    END IF  

  ELSE IF ((sw_control(j)%i_cloud_representation                     &
           == ip_cloud_conv_strat) .OR.                              &
           (sw_control(j)%i_cloud_representation                     &
           == ip_cloud_csiw)) THEN
    IF (sw_control(j)%i_inhom == ip_mcica) THEN
      ErrorStatus = 100
      CMessage = 'McICA is not compatible with the selected'//       &
                 'cloud representation'
      CALL ereport(RoutineName, ErrorStatus, CMessage) 
    ELSE
      sw_control(j)%i_solver=ip_solver_triple_hogan
      IF (sw_control(j)%i_overlap == ip_max_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_triple
      ELSE IF (sw_control(j)%i_overlap == ip_exponential_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_part_corr_cnv
      ELSE IF (sw_control(j)%i_overlap == ip_rand) THEN
        sw_control(j)%i_cloud=ip_cloud_mix_random
      ELSE
        ErrorStatus = 100
        CMessage = 'The selected cloud overlap is not available'
        CALL ereport(RoutineName, ErrorStatus, CMessage) 
      END IF
    END IF
  END IF
  IF (sw_control(j)%i_inhom == ip_scaling) THEN
    l_inhom_cloud=.TRUE.
  END IF

END DO

IF (lhook) CALL dr_hook('sw_input',zhook_out,zhook_handle)

END SUBROUTINE sw_input

END MODULE sw_rad_input_mod
