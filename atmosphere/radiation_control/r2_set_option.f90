! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Subroutine to set algorithmic options.
!
! Purpose:
!   Algorithmic options and array sizes to be set interactively
!   are determined.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------
      SUBROUTINE r2_set_option(ierr                                     &
        , n_layer, spec                                                 &
        , isolir                                                        &
        , i_gas_overlap_in, i_gas_overlap                               &
        , l_aerosol_enabled, l_climat_aerosol                           &
        , l_use_sulpc_direct, l_use_soot_direct, l_use_biogenic         &
        , l_use_dust, l_use_bmass_direct, l_use_ocff_direct             &
        , l_use_nitrate_direct                                          &
        , l_use_seasalt_direct, l_murk_rad, l_aerosol                   &
        , l_use_sulpc_indirect, l_aerosol_ccn, n_arcl_species           &
        , l_global_cloud_top, global_cloud_top, n_cloud_top_global      &
        , l_clear, i_angular_integration, i_solver, i_solver_clear      &
        , l_rescale, n_order_forward                                    &
        , i_truncation, ls_global_trunc, l_euler_trnf, euler_factor     &
        , i_scatter_method, i_scatter_method_band                       &
        , l_rad_tile                                                    &
        , weight_band                                                   &
        , nd_overlap_coeff, nd_2sg_profile, nd_layer_clr                &
        , nd_source_coeff, nd_max_order, nd_sph_coeff                   &
        , nd_region                                                     &
        , nd_profile                                                    &
        )



!     Modules included
      USE rad_pcf
      USE dec_spec
      USE tileid3z

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY: ereport
      IMPLICIT NONE


!     Dummy variables:

!     Dimensions of arrays:
      INTEGER                                                           &
                !, intent(in)
          nd_profile
!           Maximum number of atmospheric profiles

      INTEGER                                                           &
                !, intent(out)
          ierr
!           Error flag

!     Actual sizes used:
      INTEGER                                                           &
                !, intent(in)
          n_layer
!           Number of atmospheric layers used in radiation

!     Spectral data:
      type (spectrum) spec

!     Aerosol Flags:
      LOGICAL                                                           &
                !, intent(in)
          l_aerosol_enabled                                             &
!           Generic radiative switch for aerosol effects
        , l_climat_aerosol                                              &
!           Flag for climatological aerosols
        , l_use_sulpc_direct                                            &
!           Flag to include the direct effect of sulphate aerosols
        , l_use_soot_direct                                             &
!           Flag to include the direct effect of soot aerosols
        , l_use_sulpc_indirect                                          &
!           Flag to include the indirect effect of sulphate aerosols
        , l_use_seasalt_direct                                          &
!           Flag to include the direct effect of seasalt aerosols
        , l_use_biogenic                                                &
!           Flag to include the direct effect of biogenic aerosols
        , l_use_dust                                                    &
!           Flag to use direct effect of mineral dust
        , l_use_bmass_direct                                            &
!           Flag to use direct radiative effect of biomass smoke
        , l_use_ocff_direct                                             &
!           Flag to use direct radiative effect of fossil-fuel oc
        , l_use_nitrate_direct                                          &
!           Flag to use direct radiative effect of nitrate aerosols
        , l_murk_rad
!           Flag to include urban aerosols

      INTEGER n_arcl_species
!           Number of species from the NWP aerosol climatology.
!           Zero is the NWP climatology is not used.

      LOGICAL                                                           &
                !, intent(out)
          l_aerosol                                                     &
!           Flag for direct aerosol effects passed to the radiation
!           code
        , l_aerosol_ccn
!           Flag for indirect aerosol effects passed to the radiation
!           code

      INTEGER                                                           &
                 !, intent(in)
          i_angular_integration                                         &
!           Method of angular integration
        , i_solver                                                      &
!           Solver selected
        , global_cloud_top                                              &
!           Cloud top in ascending order
        , i_gas_overlap_in                                              &
!           Overall option for gaseous overlaps
        , isolir                                                        &
!           Spectral region
        , i_scatter_method
!           Treatment of scattering supplied in the controlling list
      LOGICAL                                                           &
                 !, intent(in)
          l_rescale                                                     &
!           Rescaling flag
        , l_euler_trnf
!           Flag to apply Euler's transformation to alternating series
      LOGICAL                                                           &
                !, intent(in)
          l_global_cloud_top
!           Flag to use a global cloud-top


!     Dimensions defined by the routine
      INTEGER                                                           &
                !, intent(out)
          nd_source_coeff                                               &
!           Size allocated for two-stream source coefficients
        , nd_max_order                                                  &
!           Size allocated for orders of spherical harmonics
        , nd_sph_coeff                                                  &
!           Size allocated for coefficients of spherical harmonics
!           (includes polar and azimuthal orders)
        , nd_overlap_coeff                                              &
!           Size allocated for cloud overlap coefficients
        , nd_region                                                     &
!           Size allocated for aggregated cloudy regions
        , nd_2sg_profile                                                &
!           Size allocated for profiles in two-stream flux arrays
        , nd_layer_clr
!           Size allocated for totally clear layers
      INTEGER                                                           &
                 !, intent(out)
          i_solver_clear                                                &
!           Clear-sky solver
        , n_cloud_top_global                                            &
!           Inverted global value for the cloud top
        , i_gas_overlap(spec%npd_band)                                  &
!           Gaseous overlaps in each band
        , ls_global_trunc                                               &
!           Order of global truncation
        , i_truncation                                                  &
!           Type of spherical truncation
        , n_order_forward
!           Order of term in phase function used to define forward
!           scattering fraction
      INTEGER                                                           &
                !, intent(out)
          i_scatter_method_band(spec%npd_band)
!           Treatment of scattering in each band
      REAL                                                              &
                !, intent(out)
          euler_factor                                                  &
!           Weighting applied to the last term of alternating series
        , weight_band(spec%npd_band)
!           Weightings applied to each band
      LOGICAL                                                           &
                !, intent(out)
          l_clear                                                       &
!           Flag for clear-sky calculations
        , l_rad_tile
!           Flag to enable surface tiling



!     Local variables.
      INTEGER                                                           &
          i
!           Loop variable

      CHARACTER (LEN=*), PARAMETER :: RoutineName = 'r2_set_option'
      CHARACTER (LEN=240) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('R2_SET_OPTION',zhook_in,zhook_handle)

!     Decide on the final options for aerosols:
      l_aerosol=l_aerosol_enabled.AND.                                  &
        (l_use_sulpc_direct.OR.                                         &
         l_use_dust.OR.                                                 &
         l_use_soot_direct.OR.                                          &
         l_use_bmass_direct.OR.                                         &
         l_use_ocff_direct.OR.                                          &
         l_use_nitrate_direct.OR.                                       &
         l_use_seasalt_direct.OR.                                       &
         l_use_biogenic.OR.                                             &
         l_murk_rad.OR.                                                 &
         l_climat_aerosol.OR.                                           &
         (n_arcl_species > 0))

!     Whilst l_aerosol_ccn is a generic flag for determining CCN
!     from aerosol, the view is currently taken that sulphate aerosols
!     must be included with all indirect effects, other aerosols
!     being additional, so l_aerosol_ccn is assigned solely from
!     l_use_sulpc_indirect.
      l_aerosol_ccn=l_use_sulpc_indirect

      IF (i_angular_integration == ip_two_stream) THEN

        IF ( (i_solver /= ip_solver_pentadiagonal).AND.                 &
             (i_solver /= ip_solver_mix_app_scat).AND.                  &
             (i_solver /= ip_solver_mix_direct).AND.                    &
             (i_solver /= ip_solver_mix_direct_hogan).AND.              &
             (i_solver /= ip_solver_homogen_direct).AND.                &
             (i_solver /= ip_solver_triple).AND.                        &
             (i_solver /= ip_solver_triple_hogan).AND.                  &
             (i_solver /= ip_solver_triple_app_scat)                    &
          ) THEN
          cmessage =                                                    &
            '*** error: an invalid solver has been selected '
          ierr=i_err_fatal
          GO TO 9999
        END IF

        nd_2sg_profile=nd_profile
        nd_source_coeff=2
        nd_max_order=1
        nd_sph_coeff=0
        IF ( (i_solver == ip_solver_triple).OR.                         &
             (i_solver == ip_solver_triple_hogan).OR.                   &
             (i_solver == ip_solver_triple_app_scat) ) THEN
          nd_overlap_coeff=18
        ELSE IF ( (i_solver == ip_solver_mix_direct).OR.                &
             (i_solver == ip_solver_mix_direct_hogan).OR.               &
             (i_solver == ip_solver_mix_app_scat)) THEN
          nd_overlap_coeff=8
        ELSE
          nd_overlap_coeff=0
        END IF

        IF ( (i_solver == ip_solver_triple).OR.                         &
             (i_solver == ip_solver_triple_hogan).OR.                   &
             (i_solver == ip_solver_triple_app_scat) ) THEN
          nd_region=3
        ELSE
          nd_region=2
        END IF

        IF (l_rescale) n_order_forward=2

        IF (l_clear) THEN

!         Select a clear-sky solver to match the main solver.
          IF (i_solver == ip_solver_pentadiagonal) THEN
            i_solver_clear=ip_solver_pentadiagonal
          ELSE IF (i_solver == ip_solver_mix_app_scat) THEN
            i_solver_clear=ip_solver_homogen_direct
          ELSE IF (i_solver == ip_solver_mix_direct) THEN
            i_solver_clear=ip_solver_homogen_direct
          ELSE IF (i_solver == ip_solver_mix_direct_hogan) THEN
            i_solver_clear=ip_solver_homogen_direct
          ELSE IF (i_solver == ip_solver_homogen_direct) THEN
            i_solver_clear=ip_solver_homogen_direct
          ELSE IF (i_solver == ip_solver_triple) THEN
            i_solver_clear=ip_solver_homogen_direct
          ELSE IF (i_solver == ip_solver_triple_hogan) THEN
            i_solver_clear=ip_solver_homogen_direct
          ELSE IF (i_solver == ip_solver_triple_app_scat) THEN
            i_solver_clear=ip_solver_homogen_direct
          ELSE
            cmessage =                                                  &
              '*** error: no clear-sky counterpart has been ' //        &
              'specified for this solver.'
            ierr=i_err_fatal
            GO TO 9999
          END IF

        END IF

!       We permit tiling of sea-ice points only with the two-stream
!       option at present. Tiling is only of use in separating
!       different components of the fluxes at the surface and in
!       particular is not relevent to the calculation of TOA radiances.
        l_rad_tile=.TRUE.

      ELSE IF (i_angular_integration == ip_spherical_harmonic) THEN

        nd_2sg_profile=1
        nd_source_coeff=0
        nd_overlap_coeff=0
        nd_region=0
        nd_max_order=ls_global_trunc+2
        IF (i_truncation == ip_trunc_triangular) THEN
          nd_sph_coeff=(ls_global_trunc+3)*(ls_global_trunc+4)/2
        ELSE IF (i_truncation == ip_trunc_azim_sym) THEN
          nd_sph_coeff=ls_global_trunc+2
        ELSE
          cmessage = 'error: illegal truncation'
          ierr=i_err_fatal
          GO TO 9999
        END IF

        IF (l_rescale) n_order_forward=ls_global_trunc+1

!       As currently implemented, Euler's transformation is applied
!       only in its most basic form, adding just half of the last
!       term in an alternating series.
        IF (l_euler_trnf) THEN
          euler_factor=0.5
        ELSE
          euler_factor=1.0
        END IF

!       Clear-sky fluxes are not available from the spherical harmonic
!       code in the same call as cloudy fluxes yet. If required, they
!       should be diagnosed by using a separate call to the code with
!       clouds switched off.
        IF ( l_clear ) THEN
          cmessage =                                                    &
            'Clear-sky fluxes not directly available in harmonics'
          ierr=i_err_fatal
          GO TO 9999
        END IF

!       We permit tiling of sea-ice points only with the two-stream
!       option at present.
        l_rad_tile=.FALSE.

      END IF



!     Set properties for individual bands.
      DO i=1, spec%n_band
        weight_band(i)=1.0
        i_gas_overlap(i)=i_gas_overlap_in
        IF (ANY(spec%i_scale_fnc(i,:) == ip_scale_ses2)) THEN
!         If SES2 scaling is used in this band then the
!         overlap must also use SES2:
          i_gas_overlap(i)=ip_overlap_mix_ses2
        END IF
!       Extend the treatment of scattering from the control structure
!       to each band.
        i_scatter_method_band(i)=i_scatter_method
      END DO


!     Invert the topmost cloudy layer if using a global value.
      IF (l_global_cloud_top) THEN
        n_cloud_top_global=n_layer+1-global_cloud_top
        nd_layer_clr=n_cloud_top_global-1
      END IF


 9999 CONTINUE
! Check error condition
      IF (ierr /= i_normal) THEN
        CALL ereport(RoutineName, ierr, cmessage)
      END IF

      IF (lhook) CALL dr_hook('R2_SET_OPTION',zhook_out,zhook_handle)
      END SUBROUTINE r2_set_option
