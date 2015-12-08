! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
Module lwrdiag_mod
!+ ---------------------------------------------------------------------
!  Module to define diagnostic structures for LW diagnostics.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------
!
!     Note for coders: Adding new diagnostics to the radiation
!                      code from 5.3
!
!     As the number of diagnostics required from the model has grown
!     the old practice of returning diagnostics from the radiation
!     code has resulted in our reaching the limit on the number of
!     continuation lines allowed. To obviate this problem diagnostics
!     are now bundled into a structure. For convenience, those, such
!     as the SW flux below 690 nm, that may be used to advance the model
!     are not included; those not calculated within the radiation code
!     itself are also omitted. Neither of these points is an absolute
!     requirement.
!
!     New diagnostics should be added within this structure by creating
!     a flag and a pointer array for the diagnostic. Under control of
!     STASH, allocate space for the diagnostic and initialize to 0
!     in the SW (to fill unlit points) in RAD_CTL and remember to
!     deallocate the space later. Below RAD_CTL, the element of the
!     structure can be treated as a normal array. 2-D indexing is
!     used to match the structure of the main fields in the model.

      IMPLICIT NONE

      TYPE StrLWDiag

        LOGICAL :: l_flux_up                          = .FALSE.
        LOGICAL :: l_flux_down                        = .FALSE.
        LOGICAL :: l_flux_up_clear                    = .FALSE.
        LOGICAL :: l_flux_down_clear                  = .FALSE.
        Logical :: l_total_cloud_cover                = .false.
        Logical :: l_clear_olr                        = .false.
        Logical :: l_surf_down_clr                    = .false.
        Logical :: l_clear_hr                         = .false.
        Logical :: L_net_flux_trop                    = .false.
        Logical :: L_down_flux_trop                   = .false.
        Logical :: L_total_cloud_on_levels            = .false.
        Logical :: L_cloud_absorptivity               = .false.
        Logical :: L_cloud_weight_absorptivity        = .false.
        Logical :: L_ls_cloud_absorptivity            = .false.
        Logical :: L_ls_cloud_weight_absorptivity     = .false.
        Logical :: L_cnv_cloud_absorptivity           = .false.
        Logical :: L_cnv_cloud_weight_absorptivity    = .false.
        Logical :: L_isccp_weights                    = .false.
        Logical :: L_isccp_cf                         = .false.
        Logical :: L_isccp_cf_tau_0_to_p3             = .false.
        Logical :: L_isccp_cf_tau_p3_to_1p3           = .false.
        Logical :: L_isccp_cf_tau_1p3_to_3p6          = .false.
        Logical :: L_isccp_cf_tau_3p6_to_9p4          = .false.
        Logical :: L_isccp_cf_tau_9p4_to_23           = .false.
        Logical :: L_isccp_cf_tau_23_to_60            = .false.
        Logical :: L_isccp_cf_tau_ge_60               = .false.
        Logical :: L_meanalbedocld                    = .false.
        Logical :: L_meantaucld                       = .false.
        Logical :: L_meanptop                         = .false.
        Logical :: L_totalcldarea                     = .false.
        Logical :: L_ls_qcl_rad                       = .false.
        Logical :: L_ls_qcf_rad                       = .false.
        Logical :: L_cc_qcl_rad                       = .false.
        Logical :: L_cc_qcf_rad                       = .false.
        Logical :: L_ls_cl_rad                        = .false.
        Logical :: L_ls_cf_rad                        = .false.
        Logical :: L_cc_cl_rad                        = .false.
        Logical :: L_cc_cf_rad                        = .false.
        Logical :: L_aod_sulphate                     = .false.
        Logical :: L_aod_dust                         = .false.
        Logical :: L_aod_seasalt                      = .false.
        Logical :: L_aod_soot                         = .false.
        Logical :: L_aod_biomass                      = .false.
        Logical :: L_aod_biogenic                     = .false.
        Logical :: L_aod_ocff                         = .false.
        Logical :: L_aod_delta                        = .false.
        Logical :: L_aod_nitrate                      = .false.
        Logical :: L_aod_total_radn                   = .false.
        Logical :: L_angst_total_radn                 = .false.
        Logical :: L_aod_prog_sulphate                = .false.
        Logical :: L_aod_prog_dust                    = .false.
        Logical :: L_aod_prog_seasalt                 = .false.
        Logical :: L_aod_prog_soot                    = .false.
        Logical :: L_aod_prog_biomass                 = .false.
        Logical :: L_aod_prog_ocff                    = .false.
        Logical :: L_aod_prog_nitrate                 = .false.
        LOGICAL :: L_aod_ukca_ait_sol                 = .false.  
        LOGICAL :: L_aod_ukca_acc_sol                 = .false.  
        LOGICAL :: L_aod_ukca_cor_sol                 = .false.  
        LOGICAL :: L_aod_ukca_ait_ins                 = .false.  
        LOGICAL :: L_aod_ukca_acc_ins                 = .false.  
        LOGICAL :: L_aod_ukca_cor_ins                 = .false.
! Logical for Radiances
        Logical :: L_toa_radiance                     = .false.


        REAL, POINTER :: flux_up(:, :, :)                     => NULL()
!                          Upward fluxes on model levels
        REAL, POINTER :: flux_down(:, :, :)                   => NULL()
!                          Downward fluxes on model levels
        REAL, POINTER :: flux_up_clear(:, :, :)               => NULL()
!                          Clear-sky upward fluxes on model levels
        REAL, POINTER :: flux_down_clear(:, :, :)             => NULL()
!                          Clear-sky downward fluxes on model levels
        Real, pointer :: total_cloud_cover(:, :)              => NULL()
!                          Total cloud cover at all grid-points
        Real, pointer :: clear_olr(:, :)                      => NULL()
!                          Clear-sky outgoing LW radiation calculated
!                          at all grid-points omitting cloud
!                          (This is known as Method II)
        Real, pointer :: surf_down_clr(:, :)                  => NULL()
!                          Clear-sky downward LW flux at the surface
!                          (Method II)
        Real, pointer :: clear_hr(:, :, :)                    => NULL()
!                          Clear-sky heating rates (N. B. These are
!                          the heating rates calculated by removing
!                          all clouds, which are not necessarily
!                          the same as the heating rates in the clear
!                          portions of grid-boxes)
        Real, pointer :: net_flux_trop(:, :)                  => NULL()
!                          The net downward flux at the tropopause
        Real, pointer :: down_flux_trop(:, :)                 => NULL()
!                          The actual downward flux at the tropopause
        Real, pointer :: total_cloud_on_levels(:, :, :)       => NULL()
!                           Total cloud on model layers
        Real, pointer :: cloud_absorptivity(:, :, :)          => NULL()
!                           Mean absorption coefficient in clouds,
!                           weighted by the cloud amount and the
!                           clear sky flux
        Real, pointer :: cloud_weight_absorptivity(:, :, :)   => NULL()
!                           Weighting factor for absorption in clouds :
!                           the product of the cloud amount and the
!                           clear-sky direct flux.
        Real, pointer :: ls_cloud_absorptivity(:, :, :)       => NULL()
!                           Mean absorption coefficient in layer clouds,
!                           weighted by the cloud amount and the
!                           clear sky flux
        Real, pointer :: ls_cloud_weight_absorptivity(:, :, :)=> NULL()
!                           Weighting factor for absorption in layer
!                           clouds : the product of the cloud amount
!                           and the clear-sky direct flux.
        Real, pointer :: cnv_cloud_absorptivity(:, :, :)      => NULL()
!                           Mean absorption coefficient in conv. clouds,
!                           weighted by the cloud amount and the
!                           clear sky flux
        Real, pointer :: cnv_cloud_weight_absorptivity(:, :, :)=> NULL()
!                           Weighting factor for absorption in conv.
!                           clouds : the product of the cloud amount
!                           and the clear-sky direct flux.
        Real, pointer :: isccp_weights(:, :)                  => NULL()
!                           Amount of isccp cloud
        Real, pointer :: isccp_cf(:, :, :)                    => NULL()
!                           Amount of isccp cloud
!                           tau >= 0.3
        Real, pointer :: isccp_cf_tau_0_to_p3(:, :, :)        => NULL()
!                           Amount of isccp cloud
!                           tau < 0.3
        Real, pointer :: isccp_cf_tau_p3_to_1p3(:, :, :)      => NULL()
!                           Amount of isccp cloud
!                           0.3 <= tau < 1.3
        Real, pointer :: isccp_cf_tau_1p3_to_3p6(:, :, :)     => NULL()
!                           Amount of isccp cloud
!                           1.3 <= tau < 3.6
        Real, pointer :: isccp_cf_tau_3p6_to_9p4(:, :, :)     => NULL()
!                           Amount of isccp cloud
!                           3.6 <= tau < 9.4
        Real, pointer :: isccp_cf_tau_9p4_to_23(:, :, :)      => NULL()
!                           Amount of isccp cloud
!                           9.4 <= tau < 23
        Real, pointer :: isccp_cf_tau_23_to_60(:, :, :)       => NULL()
!                           Amount of isccp cloud
!                           23 <= tau < 60
        Real, pointer :: isccp_cf_tau_ge_60(:, :, :)          => NULL()
!                           Amount of isccp cloud
!                           60 <= tau
        Real, pointer :: meanalbedocld(:, :)                  => NULL()
!                           Mean cloud albedo
        Real, pointer :: meantaucld(:, :)                     => NULL()
!                           Mean cloud optical depth
        Real, pointer :: meanptop(:, :)                       => NULL()
!                           Mean cloud top pressure
        Real, pointer :: totalcldarea(:, :)                   => NULL()
!                           Total cloud area
        Real, pointer :: ls_qcl_rad(:, :, :)                  => NULL()
!                           Grid-box mean LS liquid water mixing ratio
        Real, pointer :: ls_qcf_rad(:, :, :)                  => NULL()
!                           Grid-box mean LS ice water mixing ratio
        Real, pointer :: cc_qcl_rad(:, :, :)                  => NULL()
!                           Grid-box mean CONV liquid water mixing ratio
        Real, pointer :: cc_qcf_rad(:, :, :)                  => NULL()
!                           Grid-box mean CONV ice water mixing ratio
        Real, pointer :: ls_cl_rad(:, :, :)                   => NULL()
!                           Grid-box mean LS liquid cloud fraction
        Real, pointer :: ls_cf_rad(:, :, :)                   => NULL()
!                           Grid-box mean LS ice cloud fraction
        Real, pointer :: cc_cl_rad(:, :, :)                   => NULL()
!                           Grid-box mean CONV liquid cloud fraction
        Real, pointer :: cc_cf_rad(:, :, :)                   => NULL()
!                           Grid-box mean CONV ice cloud fraction
        Real, pointer :: aod_sulphate(:, :, :)                => NULL()
!                           Sulphate aerosol optical depth
        Real, pointer :: aod_dust(:, :, :)                    => NULL()
!                           Mineral dust aerosol optical depth
        Real, pointer :: aod_seasalt(:, :, :)                 => NULL()
!                           Sea salt aerosol optical depth
        Real, pointer :: aod_soot(:, :, :)                    => NULL()
!                           Soot (black-carbon) aerosol optical depth
        Real, pointer :: aod_biomass(:, :, :)                 => NULL()
!                           Biomass-burning aerosol optical depth
        Real, pointer :: aod_biogenic(:, :, :)                => NULL()
!                           Biogenic aerosol optical depth
        Real, pointer :: aod_ocff(:, :, :)                    => NULL()
!                           Fossil-fuel org carb aerosol optical depth
        Real, pointer :: aod_delta(:, :, :)                   => NULL()
!                           Delta aerosol optical depth
        Real, pointer :: aod_nitrate(:, :, :)                 => NULL()
!                           Nitrate aerosol optical depth
        Real, pointer :: aod_total_radn(:, :, :)              => NULL()
!                           Total aerosol optical depth in radiation
        Real, pointer :: angst_total_radn(:, :, :)             => NULL()
!                           Angstom Exp from Total AOD in radiation
        Real, pointer :: aod_prog_sulphate(:, :, :)           => NULL()
!                Prognostic Sulphate aerosol optical depth (radn or not)
        Real, pointer :: aod_prog_dust(:, :, :)               => NULL()
!                Prognostic Mineral dust aerosol optical depth
        Real, pointer :: aod_prog_seasalt(:, :, :)            => NULL()
!                Seasalt AOD diagnosed from prognostics
        Real, pointer :: aod_prog_soot(:, :, :)               => NULL()
!                Prognostic Soot (black-carbon) aerosol optical depth
        Real, pointer :: aod_prog_biomass(:, :, :)            => NULL()
!                Prognostic Biomass-burning aerosol optical depth
        Real, pointer :: aod_prog_ocff(:, :, :)               => NULL()
!                Prognostic Fossil-fuel org carb aerosol optical depth
        Real, pointer :: aod_prog_nitrate(:, :, :)            => NULL()
!                Prognostic Nitrate aerosol optical depth
        REAL, POINTER :: aod_ukca_ait_sol(:, :, :)            => NULL()  
!                           UKCA Aitken soluble mode optical depth  
        REAL, POINTER :: aod_ukca_acc_sol(:, :, :)            => NULL()  
!                           UKCA accum. soluble mode optical depth  
        REAL, POINTER :: aod_ukca_cor_sol(:, :, :)            => NULL()  
!                           UKCA coarse soluble mode optical depth  
        REAL, POINTER :: aod_ukca_ait_ins(:, :, :)            => NULL()  
!                           UKCA Aitken insoluble mode optical depth  
        REAL, POINTER :: aod_ukca_acc_ins(:, :, :)            => NULL()  
!                           UKCA accum. insoluble mode optical depth  
        REAL, POINTER :: aod_ukca_cor_ins(:, :, :)            => NULL()  
!                           UKCA coarse insoluble mode optical depth
        Real, Pointer :: toa_radiance(:, :, :)                => NULL()
!                          The radiance observed at the top of
!                          the atmosphere
!
!
      END TYPE StrLWDiag
!     ------------------------------------------------------------------
End Module lwrdiag_mod
