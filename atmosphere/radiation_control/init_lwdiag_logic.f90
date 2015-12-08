! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Initialise all the LW diagnostics to false
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
! Description of Code:
!    FORTRAN 77/90  with extensions listed in documentation.
!
!-----------------------------------------------------------------------
!
Subroutine init_lwdiag_logic(j_lw)

  Use lw_diag_mod
  USE parkind1, ONLY: jprb, jpim
  USE yomhook, ONLY: lhook, dr_hook
  Implicit None

! arguments with intent in/out

  Integer, intent(in) :: j_lw      ! call to LW radiation

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('INIT_LWDIAG_LOGIC',zhook_in,zhook_handle)

!
! Switches enabled as STASHflags: LW
!
! Fluxes and Heating Rate
!

       LW_diag(j_lw)%l_flux_up                       = .FALSE.
       LW_diag(j_lw)%l_flux_down                     = .FALSE.
       LW_diag(j_lw)%l_flux_up_clear                 = .FALSE.
       LW_diag(j_lw)%l_flux_down_clear               = .FALSE.
       LW_diag(j_lw)%l_total_cloud_cover             = .false.
       LW_diag(j_lw)%l_clear_olr                     = .false.
       LW_diag(j_lw)%l_surf_down_clr                 = .false.
       LW_diag(j_lw)%l_clear_hr                      = .false.
       LW_diag(j_lw)%L_net_flux_trop                 = .false.
       LW_diag(j_lw)%L_down_flux_trop                = .false.

! Extinction and absorptivity diagnostics

       LW_diag(j_lw)%L_total_cloud_on_levels         = .false.
       LW_diag(j_lw)%L_cloud_absorptivity            = .false.
       LW_diag(j_lw)%L_cloud_weight_absorptivity     = .false.
       LW_diag(j_lw)%L_ls_cloud_absorptivity         = .false.
       LW_diag(j_lw)%L_ls_cloud_weight_absorptivity  = .false.
       LW_diag(j_lw)%L_cnv_cloud_absorptivity        = .false.
       LW_diag(j_lw)%L_cnv_cloud_weight_absorptivity = .false.

!  Isccp diagnostics

       LW_diag(j_lw)%L_isccp_weights                 = .false.
       LW_diag(j_lw)%L_isccp_cf                      = .false.
       LW_diag(j_lw)%L_isccp_cf_tau_0_to_p3          = .false.
       LW_diag(j_lw)%L_isccp_cf_tau_p3_to_1p3        = .false.
       LW_diag(j_lw)%L_isccp_cf_tau_1p3_to_3p6       = .false.
       LW_diag(j_lw)%L_isccp_cf_tau_3p6_to_9p4       = .false.
       LW_diag(j_lw)%L_isccp_cf_tau_9p4_to_23        = .false.
       LW_diag(j_lw)%L_isccp_cf_tau_23_to_60         = .false.
       LW_diag(j_lw)%L_isccp_cf_tau_ge_60            = .false.
       LW_diag(j_lw)%L_meanalbedocld                 = .false.
       LW_diag(j_lw)%L_meantaucld                    = .false.
       LW_diag(j_lw)%L_meanptop                      = .false.
       LW_diag(j_lw)%L_totalcldarea                  = .false.

! Grid-box mean cloud diagnostics
    
       LW_diag(j_lw)%L_ls_qcl_rad                    = .false.
       LW_diag(j_lw)%L_ls_qcf_rad                    = .false.
       LW_diag(j_lw)%L_cc_qcl_rad                    = .false.
       LW_diag(j_lw)%L_cc_qcf_rad                    = .false.
       LW_diag(j_lw)%L_ls_cl_rad                     = .false.
       LW_diag(j_lw)%L_ls_cf_rad                     = .false.
       LW_diag(j_lw)%L_cc_cl_rad                     = .false.
       LW_diag(j_lw)%L_cc_cf_rad                     = .false.

! Aerosol optical depth diagnostics

       LW_diag(j_lw)%L_aod_sulphate                  = .false.
       LW_diag(j_lw)%L_aod_dust                      = .false.
       LW_diag(j_lw)%L_aod_seasalt                   = .false.
       LW_diag(j_lw)%L_aod_soot                      = .false.
       LW_diag(j_lw)%L_aod_biomass                   = .false.
       LW_diag(j_lw)%L_aod_biogenic                  = .false.
       LW_diag(j_lw)%L_aod_ocff                      = .false.
       LW_diag(j_lw)%L_aod_delta                     = .false.
       LW_diag(j_lw)%L_aod_nitrate                   = .false.
       LW_diag(j_lw)%L_aod_total_radn                = .false.
       LW_diag(j_lw)%L_angst_total_radn              = .false.
       LW_diag(j_lw)%L_aod_prog_sulphate             = .false.
       LW_diag(j_lw)%L_aod_prog_dust                 = .false.
       LW_diag(j_lw)%L_aod_prog_seasalt              = .false.
       LW_diag(j_lw)%L_aod_prog_soot                 = .false.
       LW_diag(j_lw)%L_aod_prog_biomass              = .false.
       LW_diag(j_lw)%L_aod_prog_ocff                 = .false.
       LW_diag(j_lw)%L_aod_prog_nitrate              = .false.

! UKCA aerosol optical depth diagnostics 
 
       LW_diag(j_lw)%L_aod_ukca_ait_sol              = .false. 
       LW_diag(j_lw)%L_aod_ukca_acc_sol              = .false. 
       LW_diag(j_lw)%L_aod_ukca_cor_sol              = .false. 
       LW_diag(j_lw)%L_aod_ukca_ait_ins              = .false. 
       LW_diag(j_lw)%L_aod_ukca_acc_ins              = .false. 
       LW_diag(j_lw)%L_aod_ukca_cor_ins              = .false.

! Radiance

       LW_diag(j_lw)%L_toa_radiance                  = .false.

  IF (lhook) CALL dr_hook('INIT_LWDIAG_LOGIC',zhook_out,zhook_handle)
  RETURN

End Subroutine init_lwdiag_logic
