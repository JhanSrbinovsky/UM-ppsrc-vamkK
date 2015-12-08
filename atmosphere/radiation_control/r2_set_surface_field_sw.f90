! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to set surface fields.

! Purpose:
!   The albedos and emissivity of the surface are set.

! Method:
!   Straightforward. Though the arrays passed to the code may depend
!   on the spectral band, the input arrays have no spectral dependence.
!   Note that BRDFs are treated as Lambertian here.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!- ---------------------------------------------------------------------
      SUBROUTINE r2_set_surface_field_sw(ierr                           &
        , n_band, ls_brdf_trunc                                         &
        , nlit, list                                                    &
        , l_ctile, l_use_spec_sea                                       &
        , land, land0p5, open_sea_albedo                                &
        , sea_ice_albedo                                                &
        , flandg, ice_fraction                                          &
        , land_albedo, weight_690nm                                     &
        , land0p5_g, flandg_g                                           &
        , n_brdf_basis_fnc, f_brdf, rho_alb                             &
        , l_rad_tile, n_point_tile, n_tile, list_tile, rho_alb_tile     &
        , list_tile_outer, index_tile                                   &
        , nd_field, nd_profile, nd_band, nd_brdf_basis_fnc              &
        , nd_brdf_trunc, nd_point_tile, nd_tile                         &
        )


      USE rad_pcf
      USE tileid3z
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE


!     Dummy variables:

!     Dimensions of arrays:
      INTEGER                                                           &
                !, intent(in)
          nd_field                                                      &
!           Size of input fields
        , nd_profile                                                    &
!           Maximum number of atmospheric profiles
        , nd_band                                                       &
!           Maximum number of spectral bands
        , nd_brdf_basis_fnc                                             &
!           Maximum number of BRDF basis functions
        , nd_brdf_trunc                                                 &
!           Maximum order of BRDF terms
        , nd_point_tile                                                 &
!           Size allocated for points to be tiled
        , nd_tile
!           Size allocated for number of tiles

      INTEGER                                                           &
                !, intent(out)
          ierr
!           Error flag

!     Actual sizes used:
      INTEGER                                                           &
                !, intent(in)
          n_band
!           Number of spectral bands

!     Lit points:
      INTEGER                                                           &
                !, intent(in)
          nlit                                                          &
!           Number of lit points
        , list(nd_field)
!           List of sunlit points

!     Surface options
      LOGICAL, INTENT(IN) :: l_ctile
!       Coastal tiling is used
      LOGICAL, INTENT(IN) :: l_use_spec_sea
!       Flag for spectrally dependent sea albedos

!     Physical properties of surfaces:
      INTEGER                                                           &
                !, intent(in)
          ls_brdf_trunc
!           Order of truncation applied to BRDFs
      LOGICAL, INTENT(IN) :: land(nd_field)
!           Land mask
      LOGICAL, INTENT(IN) :: land0p5(nd_field)
!           Land mask, true where the land fraction exceeds 0.5
      REAL                                                              &
                !, intent(in)
          open_sea_albedo(nd_field, 2)                                  &
!           Diffuse albedo field
        , flandg(nd_field)                                              &
!           Fraction of land in a grid-box
        , sea_ice_albedo(nd_field, 4)                                   &
!           Sea-ice albedos
        , land_albedo(nd_field, 4)                                      &
!           Land surface albedo fields
        , weight_690nm(nd_band)                                         &
!           Weights for each band for region below 690 nm
        , ice_fraction(nd_field)
!           Fraction of sea ice


!     Surface properties set.
      INTEGER                                                           &
                !, intent(out)
          n_brdf_basis_fnc
!           Number of basis functions for BRDFs
      REAL                                                              &
                !, intent(out)
          f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                  &
            , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)                     &
!           Basis functions for the surface
        , rho_alb(nd_profile, nd_brdf_basis_fnc, nd_band)
!           Weights of the surface BRDF basis functions

!     Gathered surface fields
      LOGICAL, INTENT(OUT), DIMENSION(nd_profile) :: land0p5_g
!       Gathered land mask: .TRUE. if land fraction > 0.5
      REAL, INTENT(OUT), DIMENSION(nd_profile) :: flandg_g
!       Gathered land fraction

!     Arrays related to tiling of the surface
      LOGICAL                                                           &
                !, intent(in)
          l_rad_tile
!           Local to allow tiling options
      INTEGER                                                           &
                !, intent(out)
          n_point_tile                                                  &
!           Number of points to tile
        , n_tile                                                        &
!           Number of tiles used
        , list_tile(nd_point_tile)                                      &
!           List of points with surface tiling
        , list_tile_outer(nd_point_tile)                                &
!           List of points with surface tiling in the full list
        , index_tile(npd_tile_type)
!           The indexing number of tiles of the given type
      REAL                                                              &
                !, intent(out)
          rho_alb_tile(nd_point_tile, nd_brdf_basis_fnc                 &
            , nd_tile, nd_band)
!           Weights for the basis functions of the BRDFs
!           at the tiled points


!     Local variables.
      INTEGER                                                           &
          i                                                             &
!           Loop variable
        , l                                                             &
!           Loop variable
        , ll
!           Loop variable

      REAL :: SpectralSea(n_band)

      CHARACTER (LEN=*), PARAMETER :: RoutineName = 'r2_set_surface_field_sw'
      CHARACTER (LEN=240) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('R2_SET_SURFACE_FIELD_SW',zhook_in,zhook_handle)

!     Spectrally dependent sea albedos
      SpectralSea = 1.0

      IF (l_use_spec_sea .AND. n_band == 6) THEN
         SpectralSea = (/ 0.000000, 1.444205, 1.799420,                 &
                         0.5470291, 0.000000, 0.000000 /)
      END IF


      DO l=1, nlit
        land0p5_g(l)=land0p5(list(l))
        flandg_g(l)=flandg(list(l))
      END DO


!     Define weightings for the basis functions of the surface
!     BRDFs: in effect these are surafce albedos. Each grid-box
!     may contain land, sea or sea-ice, which are treated by tiles
!     within the radiation scheme (which are used more generally
!     then simply when the coastal tiling scheme is on).

!     Note: Without coastal tiling the land albedo contains the
!     albedo of the solid surface. If coastal tiling is enabled,
!     the land albedo refers only to the land surface and there
!     is a separate sea-ice albedo.

!     The land_albedo array allows for a spectral dependence,
!     split between the VIS and NIR, as well as direct/diffuse.

      DO i=1, n_band
        DO l=1, nlit


!         Oceanic surface.
          IF (flandg(list(l)) < 1.0) THEN
            rho_alb(l, ip_surf_alb_diff, i)                             &
              =(weight_690nm(i)*sea_ice_albedo(list(l), 2)              &
              +(1.0-weight_690nm(i))*sea_ice_albedo(list(l), 4))        &
              *ice_fraction(list(l))                                    &
              + MAX(0.0, MIN(1.0, open_sea_albedo(list(l), 2)*          &
              SpectralSea(i) ) )*(1.0-ice_fraction(list(l)))
            rho_alb(l, ip_surf_alb_dir, i)                              &
              =(weight_690nm(i)*sea_ice_albedo(list(l), 1)              &
              +(1.0-weight_690nm(i))*sea_ice_albedo(list(l), 3))        &
              *ice_fraction(list(l))                                    &
              + MAX(0.0, MIN(1.0, open_sea_albedo(list(l), 1)*          &
              SpectralSea(i) ) )*(1.0-ice_fraction(list(l)))
          ELSE
            rho_alb(l, ip_surf_alb_diff, i) = 0.0
            rho_alb(l, ip_surf_alb_dir, i) = 0.0
          END IF

!         Add contributions from the land.
          IF (flandg(list(l)) > 0.0) THEN
              IF (l_ctile) THEN

                rho_alb(l, ip_surf_alb_diff, i)                         &
                  =(1.0-flandg_g(l))*rho_alb(l, ip_surf_alb_diff, i)    &
                  +flandg_g(l) * (weight_690nm(i)                       &
                  *land_albedo(list(l), 2)                              &
                  +(1.0-weight_690nm(i))*land_albedo(list(l), 4))
                rho_alb(l, ip_surf_alb_dir, i)                          &
                  =(1.0-flandg_g(l))*rho_alb(l, ip_surf_alb_dir, i)     &
                  +flandg_g(l) * (weight_690nm(i)                       &
                  *land_albedo(list(l), 1)                              &
                  +(1.0-weight_690nm(i))*land_albedo(list(l), 3))
              ELSE
                rho_alb(l, ip_surf_alb_diff, i)                         &
                  =weight_690nm(i)*land_albedo(list(l), 2)              &
                  +(1.0-weight_690nm(i))*land_albedo(list(l), 4)
                rho_alb(l, ip_surf_alb_dir, i)                          &
                  =weight_690nm(i)*land_albedo(list(l), 1)              &
                  +(1.0-weight_690nm(i))*land_albedo(list(l), 3)
              END IF
          END IF

        END DO


      END DO

! Check albedo is in physical limits:
      DO i=1, n_band 
        DO l=1, nlit 
          rho_alb(l, ip_surf_alb_dir, i) = MAX(0.0, MIN(1.0,            & 
                                 rho_alb(l, ip_surf_alb_dir, i) ) ) 
          rho_alb(l, ip_surf_alb_diff, i) = MAX(0.0, MIN(1.0,           & 
                                 rho_alb(l, ip_surf_alb_diff, i) ) ) 
        END DO 
      END DO 

!     Set the surface basis functions for a Lambertian surface.
      n_brdf_basis_fnc=1
!     By setting F_{1,0,0,0} equal to 4 we can set rho_alb equal to
!     the diffuse albedo.
      f_brdf(1, 0, 0, 0)=4.0
      IF (ls_brdf_trunc /= 0) THEN
        cmessage = 'error: the order of surface truncation is too high.'
        ierr=i_err_fatal
        GO TO 9999
      END IF


      IF (l_rad_tile) THEN

!       Set up the surface tiling variables. There are multiple levels
!       of indexing. Over all points in the domain only those in the
!       array list require radiative calculations and of these points
!       only those in the array list_file require tiling, this array
!       being indexed over points where radiative calculations are to
!       be done. list_tile_outer is indexed over the tiled points and
!       gives the index in the whole domain.

        IF (l_ctile) THEN

!         With coastal tiling we can have land, open sea or sea ice
!         in the grid-box.

          n_tile=3
          index_tile(ip_ocean_tile)=1
          index_tile(ip_seaice_tile)=2
          index_tile(ip_land_tile)=3
          n_point_tile=0
          DO ll=1, nlit
            l=list(ll)
            IF ( (flandg(l) < 1.0) .AND.                                &
                 ( (flandg(l) > 0.0) .OR.                               &
                   ( (ice_fraction(l) > 0.0) .AND.                      &
                     (ice_fraction(l) < 1.0) ) ) ) THEN
              n_point_tile=n_point_tile+1
              list_tile(n_point_tile)=ll
              list_tile_outer(n_point_tile)=l
            END IF
          END DO

        ELSE

!         Without coastal tiling we have only open sea or sea ice
!         forming the coastal tiling.

          n_tile=2
          index_tile(ip_ocean_tile)=1
          index_tile(ip_seaice_tile)=2
          n_point_tile=0
          DO ll=1, nlit
            l=list(ll)
            IF ( (.NOT.land(l)).AND.                                    &
                 (ice_fraction(l) >  0.0).AND.                          &
                 (ice_fraction(l) <  1.0) ) THEN
              n_point_tile=n_point_tile+1
              list_tile(n_point_tile)=ll
              list_tile_outer(n_point_tile)=l
            END IF
          END DO

        END IF

!       Now assign the tiled surface properties at points where tiling
!       is active.
        DO i=1, n_band

!         The oceanic surface.
          rho_alb_tile(1:n_point_tile                                   &
            , ip_surf_alb_dir, ip_ocean_tile, i)                        &
            =open_sea_albedo(list_tile_outer(1:n_point_tile), 1)        &
              *SpectralSea(i)
          rho_alb_tile(1:n_point_tile                                   &
            , ip_surf_alb_diff, ip_ocean_tile, i)                       &
            =open_sea_albedo(list_tile_outer(1:n_point_tile), 2)        &
              *SpectralSea(i)

            IF (l_ctile) THEN

!             With coastal tiling there is a real distinction between
!             land and seaice.

!             Seaice
              rho_alb_tile(1:n_point_tile                               &
                , ip_surf_alb_dir, ip_seaice_tile, i)                   &
                =weight_690nm(i)                                        &
                *sea_ice_albedo(list_tile_outer(1:n_point_tile), 1)     &
                +(1.0-weight_690nm(i))                                  &
                *sea_ice_albedo(list_tile_outer(1:n_point_tile), 3)
              rho_alb_tile(1:n_point_tile                               &
                , ip_surf_alb_diff, ip_seaice_tile, i)                  &
                =weight_690nm(i)                                        &
                *sea_ice_albedo(list_tile_outer(1:n_point_tile), 2)     &
                +(1.0-weight_690nm(i))                                  &
                *sea_ice_albedo(list_tile_outer(1:n_point_tile), 4)

!             Land
              rho_alb_tile(1:n_point_tile                               &
                , ip_surf_alb_dir, ip_land_tile, i)                     &
                =weight_690nm(i)                                        &
                *land_albedo(list_tile_outer(1:n_point_tile), 1)        &
                +(1.0-weight_690nm(i))                                  &
                *land_albedo(list_tile_outer(1:n_point_tile), 3)
              rho_alb_tile(1:n_point_tile                               &
                , ip_surf_alb_diff, ip_land_tile, i)                    &
                =weight_690nm(i)                                        &
                *land_albedo(list_tile_outer(1:n_point_tile), 2)        &
                +(1.0-weight_690nm(i))                                  &
                *land_albedo(list_tile_outer(1:n_point_tile), 4)
            ELSE

!             The land albedo fields contain the values for sea-ice.
              rho_alb_tile(1:n_point_tile                               &
                , ip_surf_alb_dir, ip_seaice_tile, i)                   &
                =weight_690nm(i)                                        &
                *land_albedo(list_tile_outer(1:n_point_tile), 1)        &
                +(1.0-weight_690nm(i))                                  &
                *land_albedo(list_tile_outer(1:n_point_tile), 3)
              rho_alb_tile(1:n_point_tile                               &
                , ip_surf_alb_diff, ip_seaice_tile, i)                  &
                =weight_690nm(i)                                        &
                *land_albedo(list_tile_outer(1:n_point_tile), 2)        &
                +(1.0-weight_690nm(i))                                  &
                *land_albedo(list_tile_outer(1:n_point_tile), 4)

            END IF

        END DO

! Check albedo is in physical limits:
        DO i=1, n_band 
          DO l=1, n_point_tile 
            DO ll =1, n_tile 
              rho_alb_tile(l, ip_surf_alb_dir, ll, i) = MAX(0.0,MIN(1.0,& 
                            rho_alb_tile(l, ip_surf_alb_dir, ll, i) ) ) 
              rho_alb_tile(l, ip_surf_alb_diff, ll, i) =MAX(0.0,MIN(1.0,& 
                            rho_alb_tile(l, ip_surf_alb_diff, ll, i) ) ) 
            END DO 
          END DO 
        END DO 

      END IF

 9999 CONTINUE
! Check error condition
      IF (ierr /= i_normal) THEN

        CALL ereport(RoutineName, ierr, cmessage)
      END IF

      IF (lhook) CALL dr_hook('R2_SET_SURFACE_FIELD_SW',zhook_out,zhook_handle)
      END SUBROUTINE r2_set_surface_field_sw
