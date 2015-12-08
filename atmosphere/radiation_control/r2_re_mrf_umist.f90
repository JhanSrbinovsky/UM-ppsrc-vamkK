! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to implement the MRF UMIST parametrization.

! Purpose:
!   Effective Radii are calculated in accordance with this
!   parametrization.

! Method:
!   The number density of CCN is found from the concentration
!   of aerosols, if available. This yields the number density of
!   droplets: if aerosols are not present, the number of droplets
!   is fixed. Effective radii are calculated from the number of
!   droplets and the LWC. Limits are applied to these values. In
!   deep convective clouds fixed values are assumed.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!- ---------------------------------------------------------------------
      SUBROUTINE r2_re_mrf_umist(n_profile, n_layer, nclds              &
         , i_gather                                                     &
         , l_pc2, l_aerosol_ccn                                         &
         , l_biomass_ccn, l_ocff_ccn, l_nitrate_ccn                     &
         , sea_salt_film, sea_salt_jet                                  &
         , l_seasalt_ccn, salt_dim_a, salt_dim_b                        &
         , l_use_biogenic, biogenic, biogenic_dim1, biogenic_dim2       &
         , accum_sulphate, diss_sulphate, aitken_sulphate               &
         , aged_bmass, cloud_bmass                                      &
         , aged_ocff, cloud_ocff                                        &
         , accum_nitrate, diss_nitrate                                  &
         , lying_snow_g                                                 &
         , i_cloud_representation                                       &
         , land_g, flandg_g                                             &
         , density_air, condensed_mix_ratio, cc_depth                   &
         , n_drop_pot, condensed_re                                     &
         , d_mass                                                       &
         , strat_liq_cloud_fraction                                     &
         , conv_liq_cloud_fraction                                      &
         , nc_diag_flag                                                 &
         , nc_diag_g                                                    &
         , nc_weight_g                                                  &
         , ntot_diag_g, ntot_land, ntot_sea                             &
         , strat_lwc_diag_g                                             &
         , so4_ccn_diag_g                                               &
         , sulp_dim1, sulp_dim2                                         &
         , bmass_dim1, bmass_dim2                                       &
         , ocff_dim1, ocff_dim2                                         &
         , nitrate_dim1, nitrate_dim2                                   &
         , npd_field, npd_profile, npd_layer, npd_aerosol_species       &
         , id_ct0                                                       &
         )


      USE rad_pcf
      USE rad_input_mod, ONLY: l_use_ndrop, lrad_ccrad

      USE water_constants_mod, ONLY: rho_water

      USE conversions_mod, ONLY: pi

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE c_micro_mod, ONLY: kparam_land, kparam_sea, dconre_land,      &
                             dconre_sea, deep_convection_limit_land,    &
                             deep_convection_limit_sea
      USE dimfix3a_mod, ONLY: npd_cloud_component, npd_cloud_type,      &
                          npd_cloud_representation, npd_overlap_coeff,  &
                          npd_source_coeff, npd_region
      USE number_droplet_mod, ONLY : number_droplet
      IMPLICIT NONE


!     COMDECKS INCLUDED:


!     DUMMY ARGUMENTS:

!     SIZES OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
           npd_field                                                    &
!             SIZE OF INPUT FIELDS TO THE RADIATION
         , npd_profile                                                  &
!             MAXIMUM NUMBER OF PROFILES
         , npd_layer                                                    &
!             MAXIMUM NUMBER OF LAYERS
         , npd_aerosol_species                                          &
!             MAXIMUM NUMBER OF AEROSOL SPECIES
         , id_ct0                                                       &
!             Topmost declared cloudy layer (only used with 3Z)
         , sulp_dim1                                                    &
!             1ST DIMENSION OF ARRAYS OF SULPHATE
         , sulp_dim2                                                    &
!             2ND DIMENSION OF ARRAYS OF SULPHATE
         , bmass_dim1                                                   &
!             1ST DIMENSION OF ARRAYS OF BIOMASS SMOKE
         , bmass_dim2                                                   &
!             2ND DIMENSION OF ARRAYS OF BIOMASS SMOKE
         , ocff_dim1                                                    &
!             1ST DIMENSION OF ARRAYS OF FOSSIL-FUEL ORGANIC CARBON
         , ocff_dim2                                                    &
!             2ND DIMENSION OF ARRAYS OF FOSSIL-FUEL ORGANIC CARBON
         , salt_dim_a                                                   &
!             1ST DIMENSION OF ARRAYS OF SEA-SALT
         , salt_dim_b                                                   &
!             2ND DIMENSION OF ARRAYS OF SEA-SALT
         , biogenic_dim1                                                &
!             1ST DIMENSION OF BIOGENIC AEROSOL ARRAY
         , biogenic_dim2                                                &
!             2ND DIMENSION OF BIOGENIC AEROSOL ARRAY
         , nitrate_dim1                                                 &
!             1ST DIMENSION OF NITRATE AEROSOL ARRAY
         , nitrate_dim2
!             2ND DIMENSION OF NITRATE AEROSOL ARRAY

      INTEGER                                                           &
                !, INTENT(IN)
           n_profile                                                    &
!             NUMBER OF ATMOSPHERIC PROFILES
         , n_layer                                                      &
!             Number of layers seen in radiation
         , nclds
!             NUMBER OF CLOUDY LEVELS

      INTEGER                                                           &
                !, INTENT(IN)
           i_gather(npd_field)
!             LIST OF POINTS TO BE GATHERED
      LOGICAL                                                           &
                !, INTENT(IN)
           land_g(npd_profile)
!             GATHERED MASK FOR LAND POINTS
      INTEGER                                                           &
                !, INTENT(IN)
           i_cloud_representation
!             REPRESENTATION OF CLOUDS

!     VARIABLES FOR PC2
      LOGICAL                                                           &
                !, INTENT(IN)
           l_pc2
!             PC2 cloud scheme is in use

!     VARIABLES FOR AEROSOLS
      LOGICAL                                                           &
                !, INTENT(IN)
           l_aerosol_ccn                                                &
!             FLAG TO USE AEROSOLS TO FIND CCN.
         , l_seasalt_ccn                                                &
!             FLAG TO USE SEA-SALT AEROSOL FOR CCN
         , l_use_biogenic                                               &
!             FLAG TO USE BIOGENIC AEROSOL FOR CCN
         , l_biomass_ccn                                                &
!             FLAG TO USE BIOMASS SMOKE AEROSOL FOR CCN
         , l_ocff_ccn                                                   &
!             FLAG TO USE FOSSIL-FUEL ORGANIC CARBON AEROSOL FOR CCN
         , l_nitrate_ccn                                                &
!             FLAG TO USE NITRATE AEROSOL FOR CCN
         , nc_diag_flag
!             FLAG TO DIAGNOSE COLUMN-INTEGRATED DROPLET NUMBER

      REAL                                                              &
                !, INTENT(IN)
           accum_sulphate(sulp_dim1, sulp_dim2)                         &
!             MIXING RATIOS OF ACCUMULATION MODE SULPHATE
         , aitken_sulphate(sulp_dim1, sulp_dim2)                        &
!             Mixing ratios of Aitken-mode sulphate
         , diss_sulphate(sulp_dim1, sulp_dim2)                          &
!             MIXING RATIOS OF DISSOLVED SULPHATE
         , aged_bmass(bmass_dim1, bmass_dim2)                           &
!             MIXING RATIOS OF AGED BIOMASS SMOKE
         , cloud_bmass(bmass_dim1, bmass_dim2)                          &
!             MIXING RATIOS OF IN-CLOUD BIOMASS SMOKE
         , aged_ocff(ocff_dim1, ocff_dim2)                              &
!             MIXING RATIOS OF AGED FOSSIL-FUEL ORGANIC CARBON
         , cloud_ocff(ocff_dim1, ocff_dim2)                             &
!             MIXING RATIOS OF IN-CLOUD FOSSIL-FUEL ORGANIC CARBON
         , sea_salt_film(salt_dim_a, salt_dim_b)                        &
!             NUMBER CONCENTRATION OF FILM-MODE SEA-SALT AEROSOL
         , sea_salt_jet(salt_dim_a, salt_dim_b)                         &
!             NUMBER CONCENTRATION OF JET-MODE SEA-SALT AEROSOL
         , biogenic(biogenic_dim1, biogenic_dim2)                       &
!             M.M.R. OF BIOGENIC AEROSOL
         , accum_nitrate(nitrate_dim1, nitrate_dim2)                    &
!             M.M.R. OF ACCUMULATION NITRATE AEROSOL
         , diss_nitrate(nitrate_dim1, nitrate_dim2)
!             M.M.R. OF DISSOLVED NITRATE AEROSOL

      REAL                                                              &
                !, INTENT(IN)
           density_air(npd_profile, npd_layer)
!             DENSITY OF AIR

      REAL                                                              &
                !, INTENT(IN)
           condensed_mix_ratio(npd_profile, id_ct0:npd_layer            &
              , npd_cloud_component)                                    &
!             MIXING RATIOS OF CONDENSED SPECIES
         , cc_depth(npd_profile)                                        &
!             DEPTH OF CONVECTIVE CLOUD
         , d_mass(npd_profile, npd_layer)                               &
!             MASS THICKNESS OF LAYER
         , strat_liq_cloud_fraction(npd_profile, npd_layer)             &
!             STRATIFORM LIQUID CLOUD COVER IN LAYERS (T>273K)
         , conv_liq_cloud_fraction(npd_profile, npd_layer)
!             CONVECTIVE LIQUID CLOUD COVER IN LAYERS (T>273K)



      REAL                                                              &
                !, Intent(IN)
           ntot_land                                                    &
                               ! Number of droplets over land / m-3
         , ntot_sea            ! Number of droplets over sea / m-3


      REAL, INTENT(IN) :: n_drop_pot(npd_field, nclds)

      REAL                                                              &
                !, INTENT(OUT)
           condensed_re(npd_profile,id_ct0:npd_layer,npd_cloud_component)
!             EFFECTIVE RADII OF CONDENSED COMPONENTS OF CLOUDS

      REAL                                                              &
                !, INTENT(OUT)
           ntot_diag_g(npd_profile, npd_layer)                          &
!             DIAGNOSTIC ARRAY FOR NTOT (GATHERED)
         , strat_lwc_diag_g(npd_profile, npd_layer)                     &
!             DIAGNOSTIC ARRAY FOR STRATIFORM LWC (GATHERED)
         , so4_ccn_diag_g(npd_profile, npd_layer)                       &
!             DIAGNOSTIC ARRAY FOR SO4 CCN MASS CONC (GATHERED)
         , nc_diag_g(npd_profile)                                       &
!             DIAGNOSTIC ARRAY FOR INTEGRATED DROPLET NUMBER (GATHERED)
         , nc_weight_g(npd_profile)
!             DIAGNOSTIC ARRAY FOR INT DROP NO. SAMPLING WGT (GATHERED)


      REAL                                                              &
           lying_snow_g(npd_profile)                                    &
!             GATHERED SNOW DEPTH (>5000m = LAND ICE SHEET)
         , flandg_g(npd_profile)
!             GATHERED GLOBAL LAND FRACTION

!     LOCAL VARIABLES:
      INTEGER                                                           &
           i                                                            &
!             LOOP VARIABLE
         , j                                                            &
!             LOOP VARIABLE
         , l                                                            &
!             LOOP VARIABLE
         , lg                                                           &
!             Index for gathering
         , idum(npd_field)                                              &
!             Temporary integer used to locate cloudy levels
         , sulphate_ptr_a                                               &
         , sulphate_ptr_b                                               &
!             POINTERS FOR SULPHATE ARRAYS
         , seasalt_ptr_a                                                &
         , seasalt_ptr_b                                                &
!             POINTERS FOR SEA-SALT ARRAYS
         , biomass_ptr_a                                                &
         , biomass_ptr_b                                                &
!             POINTERS FOR BIOMASS SMOKE ARRAYS
         , ocff_ptr_a                                                   &
         , ocff_ptr_b                                                   &
!             POINTERS FOR FOSSIL-FUEL ORGANIC CARBON ARRAYS
         , biogenic_ptr_a                                               &
         , biogenic_ptr_b                                               &
!             POINTERS FOR BIOGENIC ARRAY
         , nitrate_ptr_a                                                &
         , nitrate_ptr_b
!             POINTERS FOR NITRATE ARRAY


! Temporary arrays for gathered sets of data in calls to 
! number_droplet
      REAL :: accum_sulphate_g(n_profile, 1, n_layer)
      REAL :: diss_sulphate_g(n_profile, 1, n_layer)
      REAL :: sea_salt_film_g(n_profile, 1, n_layer)
      REAL :: sea_salt_jet_g(n_profile, 1, n_layer)
      REAL :: biogenic_g(n_profile, 1, n_layer)
      REAL :: aged_bmass_g(n_profile, 1, n_layer)
      REAL :: cloud_bmass_g(n_profile, 1, n_layer)
      REAL :: aged_ocff_g(n_profile, 1, n_layer)
      REAL :: cloud_ocff_g(n_profile, 1, n_layer)
      REAL :: accum_nitrate_g(n_profile, 1, n_layer)
      REAL :: diss_nitrate_g(n_profile, 1, n_layer)
      REAL :: density_air_g(n_profile, 1, n_layer)
      REAL :: lying_snow_g_loc(n_profile,1)
      REAL :: flandg_g_loc(n_profile,1)
      REAL :: n_drop_g(n_profile, 1, n_layer)

      REAL &
           total_mix_ratio_st(npd_profile)                              &
!             TOTAL MIXING RATIO OF WATER SUBSTANCE IN STRATIFORM CLOUD
         , total_mix_ratio_cnv(npd_profile)                             &
!             TOTAL MIXING RATIO OF WATER SUBSTANCE IN STRATIFORM CLOUD
         , total_strat_liq_cloud_fraction(npd_profile)                  &
!             TOTAL STRATIFORM LIQUID CLOUD COVER (T>273K)
         , total_conv_liq_cloud_fraction(npd_profile)
!             TOTAL CONVECTIVE LIQUID CLOUD COVER (T>273K)

      REAL                                                              &
           n_drop(npd_profile, npd_layer)                               &
!             NUMBER DENSITY OF DROPLETS
         , kparam, kparam_arr(npd_profile)                              &
!             RATIO OF CUBES OF VOLUME RADIUS TO EFFECTIVE RADIUS
         , temp
!             Temporary in calculation of effective radius

!     FIXED CONSTANTS OF THE PARAMETRIZATION:
      REAL                                                              &
           deep_convective_cloud
!             THRESHOLD VALUE FOR DEEP CONVECTIVE CLOUD
      PARAMETER(                                                        &
           deep_convective_cloud=5.0e+02                                &
         )


      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('R2_RE_MRF_UMIST',zhook_in,zhook_handle)


!     CALCULATE THE NUMBER DENSITY OF DROPLETS
      IF (l_use_ndrop) THEN

         DO i=n_layer+1-nclds, n_layer
            DO l=1, n_profile
              n_drop(l, i) = n_drop_pot(i_gather(l), n_layer+1-i)
            END DO
         END DO

      ELSE IF (l_aerosol_ccn) THEN

         DO l = 1, n_profile
           flandg_g_loc(l,1) = flandg_g(l)
           lying_snow_g_loc(l,1) = lying_snow_g(l)
         END DO

         DO i=n_layer+1-nclds, n_layer
            DO l=1, n_profile

               density_air_g(l,1,i) = density_air(l,i)
               sulphate_ptr_a=i_gather(l)
               sulphate_ptr_b=n_layer+1-i

               accum_sulphate_g(l,1,i) =                                  &
                       accum_sulphate(sulphate_ptr_a, sulphate_ptr_b)
               diss_sulphate_g(l,1,i) =                                   &
                       diss_sulphate(sulphate_ptr_a, sulphate_ptr_b)
!  Diagnose SO4 aerosol concentrations. Mass mixing ratio of ammonium
!  sulphate is converted to microgrammes of the sulphate ion per m3
!  for diagnostic purposes.

               so4_ccn_diag_g(l, i)=                                    &
                  (aitken_sulphate(sulphate_ptr_a, sulphate_ptr_b)      &
                   +accum_sulphate(sulphate_ptr_a, sulphate_ptr_b)      &
                    +diss_sulphate(sulphate_ptr_a, sulphate_ptr_b))     &
                    * density_air(l, i) * (96./132.) * 1.0e+09

              IF (l_seasalt_ccn) THEN
                 seasalt_ptr_a=i_gather(l)
                 seasalt_ptr_b=n_layer+1-i
                 sea_salt_film_g(l,1,i) =                                   &
                         sea_salt_film(seasalt_ptr_a, seasalt_ptr_b)
                 sea_salt_jet_g(l,1,i) =                                    &
                         sea_salt_jet(seasalt_ptr_a, seasalt_ptr_b)
              END IF

              IF (l_biomass_ccn) THEN
                 biomass_ptr_a=i_gather(l)
                 biomass_ptr_b=n_layer+1-i
                 aged_bmass_g(l,1,i) =                                      &
                          aged_bmass(biomass_ptr_a, biomass_ptr_b)
                 cloud_bmass_g(l,1,i) =                                     &
                          cloud_bmass(biomass_ptr_a, biomass_ptr_b)
              END IF

              IF (l_ocff_ccn) THEN
                 ocff_ptr_a=i_gather(l)
                 ocff_ptr_b=n_layer+1-i
                 aged_ocff_g(l,1,i) =                                      &
                          aged_ocff(ocff_ptr_a, ocff_ptr_b)
                 cloud_ocff_g(l,1,i) =                                     &
                          cloud_ocff(ocff_ptr_a, ocff_ptr_b)
              END IF

              IF (l_use_biogenic) THEN
                 biogenic_ptr_a=i_gather(l)
                 biogenic_ptr_b=n_layer+1-i
                 biogenic_g(l,1,i) =                                      &
                          biogenic(biogenic_ptr_a, biogenic_ptr_b)
              END IF

              IF (l_nitrate_ccn) THEN
                 nitrate_ptr_a=i_gather(l)
                 nitrate_ptr_b=n_layer+1-i
                 accum_nitrate_g(l,1,i) =                                 &
                          accum_nitrate(nitrate_ptr_a, nitrate_ptr_b)
                 diss_nitrate_g(l,1,i) =                                  &
                          diss_nitrate(nitrate_ptr_a, nitrate_ptr_b)
              END IF

            END DO
         END DO

        CALL number_droplet(                                        &
             1, n_profile,                                          &
             1, 1,                                                  &
             1, n_layer,                                            &
             n_layer+1-nclds, n_layer,                              &
             l_aerosol_ccn, .TRUE.,                                 &
             accum_sulphate_g,                                      &
             diss_sulphate_g,                                       &
             l_seasalt_ccn,                                         &
             sea_salt_film_g,                                       &
             sea_salt_jet_g,                                        &
             l_use_biogenic,                                        &
             biogenic_g,                                            &
             l_biomass_ccn,                                         &
             aged_bmass_g,                                          &
             cloud_bmass_g,                                         &
             l_ocff_ccn,                                            &
             aged_ocff_g,                                           &
             cloud_ocff_g,                                          &
             l_nitrate_ccn,                                         &
             accum_nitrate_g,                                       &
             diss_nitrate_g,                                        &
             density_air_g,                                         &
             lying_snow_g_loc,                                      &
             flandg_g_loc,                                          &
             ntot_land, ntot_sea,                                   &
             n_drop_g )

        DO i=n_layer+1-nclds, n_layer
          DO l=1, n_profile
            n_drop(l,i) = n_drop_g(l,1,i)
          END DO
        END DO


      ELSE

!        Without aerosols, the number of droplets is fixed.

         DO l=1, n_profile
            IF (flandg_g(l) >= 0.5) THEN
               n_drop(l, n_layer+1-nclds:n_layer) = ntot_land
            ELSE
               n_drop(l, n_layer+1-nclds:n_layer) = ntot_sea
            END IF
         END DO

      END IF

!  Diagnose column-integrated cloud droplet number if required.

      IF (nc_diag_flag) THEN

! DEPENDS ON: r2_calc_total_cloud_cover
         CALL r2_calc_total_cloud_cover(n_profile, n_layer, nclds       &
            , i_cloud_representation, strat_liq_cloud_fraction          &
            , total_strat_liq_cloud_fraction                            &
            , npd_profile, npd_layer)

! DEPENDS ON: r2_calc_total_cloud_cover
         CALL r2_calc_total_cloud_cover(n_profile, n_layer, nclds       &
            , i_cloud_representation, conv_liq_cloud_fraction           &
            , total_conv_liq_cloud_fraction                             &
            , npd_profile, npd_layer)

! DEPENDS ON: r2_column_droplet_conc
         CALL r2_column_droplet_conc(npd_profile, npd_layer             &
            , n_profile, n_layer, nclds                                 &
            , strat_liq_cloud_fraction                                  &
            , total_strat_liq_cloud_fraction                            &
            , conv_liq_cloud_fraction                                   &
            , total_conv_liq_cloud_fraction                             &
            , n_drop, d_mass, density_air                               &
            , nc_diag_g, nc_weight_g)

      END IF


      IF (i_cloud_representation == ip_cloud_ice_water) THEN

        ! We only need to find land across n_profile
        DO i = 1, n_profile
          IF (land_g(i)) THEN
            kparam_arr(i)=kparam_land
          ELSE
            kparam_arr(i)=kparam_sea
          END IF
        END DO
        DO i=n_layer+1-nclds, n_layer

!         Find the total mixing ratio of water substance in the cloud.
          DO l=1, n_profile
            total_mix_ratio_st(l)                                       &
              = condensed_mix_ratio(l, i, ip_clcmp_st_water)
            total_mix_ratio_cnv(l) = 0.0e+00

            condensed_re(l, i, ip_clcmp_cnv_water) = 0.0e+00
            condensed_re(l, i, ip_clcmp_st_water)  = MAX(0.0e+00,       &
              3.0e+00*total_mix_ratio_st(l)*density_air(l, i)           &
              /(4.0e+00*pi*rho_water*kparam_arr(l)*n_drop(l, i)))       &
              **(1.0e+00/3.0e+00)

            ntot_diag_g(l, i)=n_drop(l, i)*1.0e-06
            strat_lwc_diag_g(l, i)                                      &
              = total_mix_ratio_st(l)*density_air(l, i)*1.0e03
          END DO

        END DO

      ELSE

      DO i=n_layer+1-nclds, n_layer

!        FIND THE TOTAL MIXING RATIO OF WATER SUBSTANCE IN THE CLOUD
!        AS IMPLIED BY THE REPRESENTATION.
         IF (i_cloud_representation == ip_cloud_conv_strat) THEN
            DO l=1, n_profile
               total_mix_ratio_st(l)                                    &
                  =condensed_mix_ratio(l, i, ip_clcmp_st_water)         &
                  +condensed_mix_ratio(l, i, ip_clcmp_st_ice)
               total_mix_ratio_cnv(l)                                   &
                  =condensed_mix_ratio(l, i, ip_clcmp_cnv_water)        &
                  +condensed_mix_ratio(l, i, ip_clcmp_cnv_ice)
            END DO
         ELSE IF (i_cloud_representation == ip_cloud_csiw) THEN
            DO l=1, n_profile
               total_mix_ratio_st(l)                                    &
                  =condensed_mix_ratio(l, i, ip_clcmp_st_water)
               total_mix_ratio_cnv(l)                                   &
                  =condensed_mix_ratio(l, i, ip_clcmp_cnv_water)
            END DO
         END IF

         IF (l_pc2) THEN
           DO l=1, n_profile
              IF (land_g(l)) THEN
                 kparam=kparam_land
              ELSE
                 kparam=kparam_sea
              END IF
              temp                                                      &
               =(3.0e+00*total_mix_ratio_cnv(l)*density_air(l, i)       &
               /(4.0e+00*pi*rho_water*kparam*n_drop(l, i)))
              IF (temp  >=  0.0) THEN
                condensed_re(l, i, ip_clcmp_cnv_water)                  &
                = temp**(1.0e+00/3.0e+00)
              ELSE
                condensed_re(l, i, ip_clcmp_cnv_water) = 0.0
              END IF
              temp                                                      &
               =(3.0e+00*total_mix_ratio_st(l)*density_air(l, i)        &
               /(4.0e+00*pi*rho_water*kparam*n_drop(l, i)))
              IF (temp  >=  0.0) THEN
                condensed_re(l, i, ip_clcmp_st_water)                   &
                = temp**(1.0e+00/3.0e+00)
              ELSE
                condensed_re(l, i, ip_clcmp_st_water) = 0.0
              END IF
           END DO
         ELSE
           DO l=1, n_profile
              IF (land_g(l)) THEN
                 kparam=kparam_land
              ELSE
                 kparam=kparam_sea
              END IF
              condensed_re(l, i, ip_clcmp_cnv_water)                    &
               =(3.0e+00*total_mix_ratio_cnv(l)*density_air(l, i)       &
               /(4.0e+00*pi*rho_water*kparam*n_drop(l, i)))             &
               **(1.0e+00/3.0e+00)
              condensed_re(l, i, ip_clcmp_st_water)                     &
               =(3.0e+00*total_mix_ratio_st(l)*density_air(l, i)        &
               /(4.0e+00*pi*rho_water*kparam*n_drop(l, i)))             &
               **(1.0e+00/3.0e+00)
           END DO
         END IF
         DO l=1, n_profile
            ntot_diag_g(l, i)=n_drop(l, i)*1.0e-06
            strat_lwc_diag_g(l, i)                                      &
               =total_mix_ratio_st(l)*density_air(l, i)*1.0e03
         END DO
      END DO

      !-----------------------------------------------------------------------
      ! Reset the effective radii for deep convective clouds.
      !-----------------------------------------------------------------------

      IF (lrad_ccrad) THEN

        ! Loop over all gridpoints which are affected by radiation
        DO l=1, n_profile

          cc_depth(l) = 0.0e+00
          lg          = i_gather(l)
          idum        = 0

          ! Loop from top layer down i.e. radiation indexes from the top of
          ! atmosphere down
          DO i=n_layer+1-nclds, n_layer


            IF (condensed_mix_ratio(l,i,ip_clcmp_cnv_water) > 0.0e+00) THEN

              IF (idum(lg) == 0) THEN
                ! Mark 1st cloud top (i.e highest cloud top in profile)
                idum(lg) = i
              END IF

            ELSE

              IF (idum(lg) /= 0) THEN
                ! Cloud top already found so this is cloud base and have found
                ! a contigious cloud bank, begin to calculate depth

                DO j=idum(lg), i-1
                  cc_depth(l) = cc_depth(l)                                   &
                              + ( d_mass(l,j)/density_air(l,j) )
                END DO


                ! Reset effective radii if clouds are physically deep
                DO j=idum(lg), i-1
                  IF (land_g(l)) THEN

                    ! LAND POINT
                    IF ((cc_depth(l) > deep_convection_limit_land) .OR.      &
                        (condensed_re(l,j,ip_clcmp_cnv_water) > dconre_land))&
                        THEN

                      condensed_re(l,j,ip_clcmp_cnv_water) = dconre_land
                    END IF

                  ELSE

                    ! SEA POINT
                    IF ((cc_depth(l) > deep_convection_limit_sea) .OR.       &
                        (condensed_re(l,j,ip_clcmp_cnv_water) > dconre_sea)) &
                        THEN

                      condensed_re(l,j,ip_clcmp_cnv_water) = dconre_sea

                    END IF

                  END IF      ! (LAND_G)
                END DO      ! J

                ! Reset cloud top and carry on looking for cloud
                idum(lg) = 0

              END IF      ! (idum)
            END IF      ! (CONDENSED_MIX_RATIO)
          END DO      ! I (N_LAYER)
        END DO      ! L (N_PROFILE)

      ELSE ! original code

        DO i=n_layer+1-nclds, n_layer
          DO l=1, n_profile

            IF (land_g(l)) THEN

              IF (cc_depth(l) > deep_convection_limit_land) THEN
                condensed_re(l,i,ip_clcmp_cnv_water) = dconre_land
              ELSE
                IF (condensed_re(l,i,ip_clcmp_cnv_water) > dconre_land) THEN
                  condensed_re(l,i,ip_clcmp_cnv_water) = dconre_land
                END IF
              END IF

            ELSE

              IF (cc_depth(l) > deep_convection_limit_sea) THEN
                condensed_re(l,i,ip_clcmp_cnv_water) = dconre_sea
              ELSE
                IF (condensed_re(l,i,ip_clcmp_cnv_water) > dconre_sea) THEN
                  condensed_re(l,i,ip_clcmp_cnv_water) = dconre_sea
                END IF
              END IF

            END IF      ! (LAND_G)

          END DO      ! L (N_PROFILE)
        END DO      ! I (N_LAYER)
      END IF      ! lrad_ccrad

      END IF ! I_CLOUD_REPRESENTATION == IP_CLOUD_ICE_WATER


      IF (lhook) CALL dr_hook('R2_RE_MRF_UMIST',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE r2_re_mrf_umist
