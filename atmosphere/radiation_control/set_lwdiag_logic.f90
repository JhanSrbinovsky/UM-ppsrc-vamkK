! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Set the LW radiation diagnostic flags to true if necessary
!
! Method:
!
! If a radiation diagnostic has been chosen in STASH then the
! flag of the corresponding radiation diagnostic in the structure
! LW_diag is set to true.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
!-----------------------------------------------------------------------
!
SUBROUTINE set_lwdiag_logic(sf,nitems,nsects, l_radiance, j_lw, i_off)

      Use lw_diag_mod
      USE rad_input_mod, ONLY: l_rad_perturb
      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      Implicit None

! arguments with intent in

      integer, intent(in) :: nitems       ! STASH item number
      integer, intent(in) :: nsects       ! STASH section number
      integer, intent(in) :: j_lw         ! call to LW radiation
      integer, intent(in) :: i_off        ! offset for diagnostics

      logical, intent(in) :: l_radiance   ! Flag for radiances

      logical, intent(in) :: sf(0:nitems,0:nsects) ! STASH Flags

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('SET_LWDIAG_LOGIC',zhook_in,zhook_handle)



! Switches enabled as STASHflags: LW

      if (l_radiance.and.(j_lw > 1)) then

        LW_diag(j_lw)%L_toa_radiance      = sf(297+i_off,2)

      else

! Fluxes

        LW_diag(j_lw)%l_flux_up           = sf(217+i_off,2)
        LW_diag(j_lw)%l_flux_down         = sf(218+i_off,2)
        LW_diag(j_lw)%L_net_flux_trop     = sf(237+i_off,2)
        LW_diag(j_lw)%L_down_flux_trop    = sf(238+i_off,2)

! Cloud amount

        LW_diag(j_lw)%L_total_cloud_cover     = sf(204+i_off,2)
        LW_diag(j_lw)%L_total_cloud_on_levels = sf(261+i_off,2)

!  Isccp diagnostics

        if (j_lw == 1) then
          LW_diag(1)%L_isccp_weights                 = sf(269,2)
          LW_diag(1)%L_isccp_cf                      = sf(270,2)
          LW_diag(1)%L_isccp_cf_tau_0_to_p3          = sf(271,2)
          LW_diag(1)%L_isccp_cf_tau_p3_to_1p3        = sf(272,2)
          LW_diag(1)%L_isccp_cf_tau_1p3_to_3p6       = sf(273,2)
          LW_diag(1)%L_isccp_cf_tau_3p6_to_9p4       = sf(274,2)
          LW_diag(1)%L_isccp_cf_tau_9p4_to_23        = sf(275,2)
          LW_diag(1)%L_isccp_cf_tau_23_to_60         = sf(276,2)
          LW_diag(1)%L_isccp_cf_tau_ge_60            = sf(277,2)
          LW_diag(1)%L_meanalbedocld                 = sf(290,2)
          LW_diag(1)%L_meantaucld                    = sf(291,2)
          LW_diag(1)%L_meanptop                      = sf(292,2)
          LW_diag(1)%L_totalcldarea                  = sf(293,2)
        end if

! Grid-box mean cloud diagnostics as seen by radiation:

        LW_diag(j_lw)%L_ls_qcl_rad                    = sf(308+i_off,2)
        LW_diag(j_lw)%L_ls_qcf_rad                    = sf(309+i_off,2)
        LW_diag(j_lw)%L_cc_qcl_rad                    = sf(310+i_off,2)
        LW_diag(j_lw)%L_cc_qcf_rad                    = sf(311+i_off,2)

        if (312+i_off <= nitems) then
          LW_diag(j_lw)%L_ls_cl_rad                   = sf(312+i_off,2)
        end if
        if (313+i_off <= nitems) then
          LW_diag(j_lw)%L_ls_cf_rad                   = sf(313+i_off,2)
        end if
        if (314+i_off <= nitems) then
          LW_diag(j_lw)%L_cc_cl_rad                   = sf(314+i_off,2)
        end if
        if (315+i_off <= nitems) then
          LW_diag(j_lw)%L_cc_cf_rad                   = sf(315+i_off,2)
        end if


      If (.not.(l_rad_perturb.and.(j_lw == 2))) Then
!       For the incremental time-stepping scheme (l_rad_perturb) many
!       of the diagnostics are not calculated on the "cloud only" 
!       radiation calls (j_lw==2).

! Clear-sky Fluxes and Heating Rates

        LW_diag(j_lw)%l_flux_up_clear     = sf(219+i_off,2)
        LW_diag(j_lw)%l_flux_down_clear   = sf(220+i_off,2)
        LW_diag(j_lw)%L_clear_olr         = sf(206+i_off,2)
        LW_diag(j_lw)%L_surf_down_clr     = sf(208+i_off,2)
        LW_diag(j_lw)%L_clear_hr          = sf(233+i_off,2)

! Extinction and absorptivity diagnostics

        LW_diag(j_lw)%L_cloud_absorptivity            = sf(262+i_off,2)
        LW_diag(j_lw)%L_cloud_weight_absorptivity     = sf(263+i_off,2)
        LW_diag(j_lw)%L_ls_cloud_absorptivity         =(sf(264+i_off,2) &
          .OR.(LW_diag(j_lw)%L_isccp_weights))
        LW_diag(j_lw)%L_ls_cloud_weight_absorptivity  =(sf(265+i_off,2) &
          .OR.(LW_diag(j_lw)%L_isccp_weights))
        LW_diag(j_lw)%L_cnv_cloud_absorptivity        =(sf(266+i_off,2) &
          .OR.(LW_diag(j_lw)%L_isccp_weights))
        LW_diag(j_lw)%L_cnv_cloud_weight_absorptivity =(sf(267+i_off,2) &
          .OR.(LW_diag(j_lw)%L_isccp_weights))

! Aerosol optical depth diagnostics
        LW_diag(j_lw)%L_aod_sulphate                  = sf(284+i_off,2)
        LW_diag(j_lw)%L_aod_dust                      = sf(285+i_off,2)
        LW_diag(j_lw)%L_aod_seasalt                   = sf(286+i_off,2)
        LW_diag(j_lw)%L_aod_soot                      = sf(287+i_off,2)
        LW_diag(j_lw)%L_aod_biomass                   = sf(288+i_off,2)
        LW_diag(j_lw)%L_aod_biogenic                  = sf(289+i_off,2)
        LW_diag(j_lw)%L_aod_ocff                      = sf(295+i_off,2)
        LW_diag(j_lw)%L_aod_delta                     = sf(296+i_off,2)
        LW_diag(j_lw)%L_aod_nitrate                   = sf(297+i_off,2)
        LW_diag(j_lw)%L_aod_total_radn                = sf(298+i_off,2)
        LW_diag(j_lw)%L_angst_total_radn              = sf(299+i_off,2)
! additional logic required to swtich these off radtv forcing diags are on
        if (421+i_off <= nitems) then
          LW_diag(j_lw)%L_aod_prog_sulphate           = sf(421+i_off,2)
        end if
        if (422+i_off <= nitems) then
          LW_diag(j_lw)%L_aod_prog_dust               = sf(422+i_off,2)
        end if
        if (423+i_off <= nitems) then
          LW_diag(j_lw)%L_aod_prog_seasalt            = sf(423+i_off,2)
        end if
        if (424+i_off <= nitems) then
          LW_diag(j_lw)%L_aod_prog_soot               = sf(424+i_off,2)
        end if
        if (425+i_off <= nitems) then
          LW_diag(j_lw)%L_aod_prog_biomass            = sf(425+i_off,2)
        end if
        if (426+i_off <= nitems) then
          LW_diag(j_lw)%L_aod_prog_ocff               = sf(426+i_off,2)
        end if
        if (427+i_off <= nitems) then
          LW_diag(j_lw)%L_aod_prog_nitrate            = sf(427+i_off,2)
        end if

! UKCA aerosol optical depth diagnostics 
        LW_diag(j_lw)%L_aod_ukca_ait_sol              = sf(300+i_off,2) 
        LW_diag(j_lw)%L_aod_ukca_acc_sol              = sf(301+i_off,2) 
        LW_diag(j_lw)%L_aod_ukca_cor_sol              = sf(302+i_off,2) 
        LW_diag(j_lw)%L_aod_ukca_ait_ins              = sf(303+i_off,2)
        LW_diag(j_lw)%L_aod_ukca_acc_ins              = sf(304+i_off,2) 
        LW_diag(j_lw)%L_aod_ukca_cor_ins              = sf(305+i_off,2)

      End If ! .not.(l_rad_perturb.and.(j_lw == 2))

      end if

      IF (lhook) CALL dr_hook('SET_LWDIAG_LOGIC',zhook_out,zhook_handle)
      RETURN

end subroutine set_lwdiag_logic
