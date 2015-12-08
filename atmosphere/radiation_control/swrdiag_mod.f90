! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE swrdiag_mod
!+ ---------------------------------------------------------------------
!  Module to define diagnostic structures for SW diagnostics.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!- ----------------------------------------------------------------------

!     Note for coders: Adding new diagnostics to the radiation
!                      code from 5.3

!     As the number of diagnostics required from the model has grown
!     the old practice of returning diagnostics from the radiation
!     code has resulted in our reaching the limit on the number of
!     continuation lines allowed. To obviate this problem diagnostics
!     are now bundled into a structure. For convenience, those, such
!     as the SW flux below 690 nm, that may be used to advance the model
!     are not included; those not calculated within the radiation code
!     itself are also omitted. Neither of these points is an absolute
!     requirement.

!     New diagnostics should be added within this structure by creating
!     a flag and a pointer array for the diagnostic. Under control of
!     STASH, allocate space for the diagnostic and initialize to 0
!     in the SW (to fill unlit points) in RAD_CTL and remember to
!     deallocate the space later. Below RAD_CTL, the element of the
!     structure can be treated as a normal array. 2-D indexing is
!     used to match the structure of the main fields in the model.

      IMPLICIT NONE

      TYPE strswdiag

        LOGICAL ::  l_flux_up                      = .FALSE.
        LOGICAL ::  l_flux_down                    = .FALSE.
        LOGICAL ::  l_flux_up_clear                = .FALSE.
        LOGICAL ::  l_flux_down_clear              = .FALSE.
        LOGICAL ::  l_solar_out_toa                = .FALSE.
        LOGICAL ::  l_solar_out_clear              = .FALSE.
        LOGICAL ::  l_surface_down_flux            = .FALSE.
        LOGICAL ::  l_surf_down_clr                = .FALSE.
        LOGICAL ::  l_surf_up_clr                  = .FALSE.
        LOGICAL ::  l_clear_hr                     = .FALSE.
        LOGICAL ::  l_net_flux_trop                = .FALSE.
        LOGICAL ::  l_up_flux_trop                 = .FALSE.
        LOGICAL ::  l_flux_direct                  = .FALSE.
        LOGICAL ::  l_flux_diffuse                 = .FALSE.
        LOGICAL ::  re_conv_flag                   = .FALSE.
        LOGICAL ::  re_strat_flag                  = .FALSE.
        LOGICAL ::  wgt_conv_flag                  = .FALSE.
        LOGICAL ::  wgt_strat_flag                 = .FALSE.
        LOGICAL ::  lwp_strat_flag                 = .FALSE.
        LOGICAL ::  weighted_re_flag               = .FALSE.
        LOGICAL ::  sum_weight_re_flag             = .FALSE.
        LOGICAL ::  wgtd_warm_re_flag              = .FALSE.
        LOGICAL ::  sum_wgt_warm_re_flag           = .FALSE.
        LOGICAL ::  ntot_diag_flag                 = .FALSE.
        LOGICAL ::  strat_lwc_diag_flag            = .FALSE.
        LOGICAL ::  so4_ccn_diag_flag              = .FALSE.
        LOGICAL ::  cond_samp_wgt_flag             = .FALSE.
        LOGICAL ::  seasalt_film_flag              = .FALSE.
        LOGICAL ::  seasalt_jet_flag               = .FALSE.
        LOGICAL ::  nc_diag_flag                   = .FALSE.
        LOGICAL ::  nc_weight_flag                 = .FALSE.
        LOGICAL ::  l_FlxSolBelow690nmSurf         = .FALSE.
        LOGICAL ::  l_FlxSeaBelow690nmSurf         = .FALSE.
        LOGICAL ::  l_cloud_extinction             = .FALSE.
        LOGICAL ::  l_cloud_weight_extinction      = .FALSE.
        LOGICAL ::  l_ls_cloud_extinction          = .FALSE.
        LOGICAL ::  l_ls_cloud_weight_extinction   = .FALSE.
        LOGICAL ::  l_cnv_cloud_extinction         = .FALSE.
        LOGICAL ::  l_cnv_cloud_weight_extinction  = .FALSE.
! Logicals for Orography
        LOGICAL ::  l_orog_corr                    = .FALSE.
        LOGICAL ::  l_sol_bearing                  = .FALSE.
! Logical for Radiances
        LOGICAL ::  l_toa_radiance                 = .FALSE.
! Logicals for UV-Fluxes
        LOGICAL ::  l_uvflux_direct                = .FALSE.
        LOGICAL ::  l_uvflux_up                    = .FALSE.
        LOGICAL ::  l_uvflux_net                   = .FALSE.
        LOGICAL ::  l_surf_uv                      = .FALSE.
        LOGICAL ::  l_surf_uv_clr                  = .FALSE.
! Logicals for Albedos 
        LOGICAL :: l_direct_albedo                 = .FALSE. 
        LOGICAL :: l_diffuse_albedo                = .FALSE. 
        LOGICAL :: l_vis_albedo_sc                 = .FALSE. 
        LOGICAL :: l_nir_albedo_sc                 = .FALSE. 


        REAL, POINTER :: flux_up(:, :, :)                     => NULL()
!                          Upward fluxes on model levels
        REAL, POINTER :: flux_down(:, :, :)                   => NULL()
!                          Downward fluxes on model levels
        REAL, POINTER :: flux_up_clear(:, :, :)               => NULL()
!                          Clear-sky upward fluxes on model levels
        REAL, POINTER :: flux_down_clear(:, :, :)             => NULL()
!                          Clear-sky downward fluxes on model levels
        REAL, POINTER :: solar_out_toa(:, :)                  => NULL()
!                          Reflected SW flux at the top
!                          of the atmosphere
        REAL, POINTER :: solar_out_clear(:, :)                => NULL()
!                          Reflected clear-sky SW flux
!                          at the top of the atmosphere:
!                          this is calculated at all points
!                          omitting cloud (Method II)
        REAL, POINTER :: surface_down_flux(:, :)              => NULL()
!                          Downward SW flux at the surface
!                          (not net)
        REAL, POINTER :: surf_down_clr(:, :)                  => NULL()
!                          Clear-sky downward SW flux at the
!                          surface (Method II)
        REAL, POINTER :: surf_up_clr(:, :)                    => NULL()
!                          Clear-sky upward SW flux at the
!                          surface (Method II)
        REAL, POINTER :: net_flux_trop(:, :)                  => NULL()
!                          Net downward flux at the tropopause
        REAL, POINTER :: up_flux_trop(:, :)                   => NULL()
!                          Actual upward flux at the tropopause
        REAL, POINTER :: flux_direct(:, :, :)                 => NULL()
!                          Direct downward SW flux
        REAL, POINTER :: flux_diffuse(:, :, :)                => NULL()
!                          Diffuse downward SW flux
        REAL, POINTER :: clear_hr(:, :, :)                    => NULL()
!                          Clear-sky heating rates, calculated
!                          by ignoring all clouds (Method II):
!                          these are not necessarily the same
!                          as the heating rates in the cloud-
!                          free parts of a grid-box
        REAL, POINTER :: re_strat(:, :, :)                    => NULL()
!                          The weighted effective radius in
!                          stratiform clouds multiplied by 10^6
!                          to avoid packing problems
        REAL, POINTER :: wgt_strat(:, :, :)                   => NULL()
!                          The weighting factor for re_strat
!                          and lwp_strat
        REAL, POINTER :: lwp_strat(:, :, :)                   => NULL()
!                          The liquid water path in stratiform
!                          cloud weighted by wgt_strat
        REAL, POINTER :: re_conv(:, :, :)                     => NULL()
!                          The weighted effective radius in
!                          Convective clouds multiplied by 10^6
!                          to avoid packing problems
        REAL, POINTER :: wgt_conv(:, :, :)                    => NULL()
!                          The weighting factor for re_conv
        REAL, POINTER :: ntot_diag(:, :, :)                   => NULL()
!                          The number concentration of droplets
!                          multiplied by stratiform weighting factor
        REAL, POINTER :: strat_lwc_diag(:, :, :)              => NULL()
!                          The liquid water content of stratiform
!                          clouds multiplied by stratiform weighting
!                          factor
        REAL, POINTER :: so4_ccn_diag(:, :, :)                => NULL()
!                          The mass concentration of SO4 CCN multiplied
!                          by the conditional sampling weight
        REAL, POINTER :: cond_samp_wgt(:, :, :)               => NULL()
!                          The conditional sampling weight for
!                          so4_ccn_diag
        REAL, POINTER :: weighted_re(:, :)                    => NULL()
!                          The effective radius as seen from space
!                          multiplied by an appropriate weight
        REAL, POINTER :: sum_weight_re(:, :)                  => NULL()
!                          The weighting factor for the effective
!                          radius as viewed from space
        REAL, POINTER :: weighted_warm_re(:, :)               => NULL()
!                          The effective radius as seen from space
!                          for warm clouds only (T>273K) multiplied
!                          by an appropriate weight
        REAL, POINTER :: sum_weight_warm_re(:, :)             => NULL()
!                          The weighting factor for the warm-cloud-
!                          only effective radius as viewed from space
        REAL, POINTER :: nc_diag(:, :)                        => NULL()
!                          The column-integrated cloud droplet number
!                          multiplied by an appropriate weight
        REAL, POINTER :: nc_weight(:, :)                      => NULL()
!                          The weighting factor for the column-
!                          integrated cloud droplet number
        REAL, POINTER :: FlxSolBelow690nmSurf(:, :)           => NULL()
!                          The grid-box mean flux below 690 nm
!                          into the solid surface
        REAL, POINTER :: FlxSeaBelow690nmSurf(:, :)           => NULL()
!                          The grid-box mean flux below 690 nm
!                          into the sea surface
        REAL, POINTER :: cloud_extinction(:, :, :)            => NULL()
!                           Mean extinction coefficient in clouds,
!                           weighted by the cloud amount and the
!                           clear sky flux
        REAL, POINTER :: cloud_weight_extinction(:, :, :)     => NULL()
!                           Weighting factor for extinction in clouds :
!                           the product of the cloud amount and the
!                           clear-sky direct flux.
        REAL, POINTER :: ls_cloud_extinction(:, :, :)         => NULL()
!                           Mean extinction coefficient in layer clouds,
!                           weighted by the cloud amount and the
!                           clear sky flux
        REAL, POINTER :: ls_cloud_weight_extinction(:, :, :)  => NULL()
!                           Weighting factor for extinction in layer
!                           clouds : the product of the cloud amount
!                           and the clear-sky direct flux.
        REAL, POINTER :: cnv_cloud_extinction(:, :, :)        => NULL()
!                           Mean extinction coefficient in conv. clouds,
!                           weighted by the cloud amount and the
!                           clear sky flux
        REAL, POINTER :: cnv_cloud_weight_extinction(:, :, :) => NULL()
!                           Weighting factor for extinction in conv.
!                           clouds : the product of the cloud amount
!                           and the clear-sky direct flux.
        REAL, POINTER :: orog_corr(:, :)                      => NULL()
!                          Correction factor for the direct solar flux
!                          reaching the surface for sloping terrain.
        REAL, POINTER :: sol_bearing(:, :)                    => NULL()
!                          Mean local bearing of the sun over the time
!                          step, in radians clockwise from grid north.
        REAL, POINTER :: toa_radiance(:, :, :)                => NULL()
!                          The radiance observed at the top of
!                          the atmosphere
        REAL, POINTER :: uvflux_direct(:, :, :)               => NULL()
!                          direct UV-flux
        REAL, POINTER :: uvflux_net(:, :, :)                  => NULL()
!                          net UV-flux
        REAL, POINTER :: uvflux_up(:, :, :)                   => NULL()
!                          upward UV-flux
        REAL, POINTER :: surf_uv(:, :)                        => NULL()
!                          Surface down UV flux
        REAL, POINTER :: surf_uv_clr(:, :)                    => NULL()
!                          Clear-sky surface down UV flux
        REAL, POINTER :: direct_albedo(:, :, :)               => NULL() 
!                          direct albedo on SW bands 
        REAL, POINTER :: diffuse_albedo(:, :, :)              => NULL() 
!                          diffuse albedo on SW bands 
        REAL, POINTER :: vis_albedo_sc(:, :, :)               => NULL() 
!                          sclaing to VIS albedo obs on land tiles 
        REAL, POINTER :: nir_albedo_sc(:, :, :)               => NULL() 
!                          sclaing to NIR albedo obs on land tiles 

      END TYPE strswdiag
! ----------------------------------------------------------------------
END MODULE swrdiag_mod
