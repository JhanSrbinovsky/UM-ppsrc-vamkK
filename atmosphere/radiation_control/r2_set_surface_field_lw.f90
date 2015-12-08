! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to set surface fields.

! Purpose:
!   The albedos and emissivity of the surface are set.

! Method:
!   Straightforward.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!- ---------------------------------------------------------------------
      SUBROUTINE r2_set_surface_field_lw(ierr                           &
        , n_points, list, n_band, ls_brdf_trunc                         &
        , l_ctile, flandg, emis_land                                    &
        , n_brdf_basis_fnc, f_brdf, rho_alb                             &
        , land, ice_fraction, t_rad_land, t_rad_sice, t_rad_sea         &
        , l_rad_tile, n_point_tile, n_tile, list_tile, rho_alb_tile     &
        , frac_tile, t_tile                                             &
        , list_tile_outer, index_tile                                   &
        , nd_field, nd_profile, nd_band, nd_brdf_basis_fnc              &
        , nd_brdf_trunc, nd_point_tile, nd_tile                         &
        )



      USE rad_pcf
      USE tileid3z
      USE surf_param, ONLY: emis_sea, emis_sice
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE


!     Dummy variables:

!     Dimensions of arrays:
      INTEGER, INTENT(IN) ::                                            &
          nd_field                                                      &
!           Allocated size of fields of data
        , nd_profile                                                    &
!           Maximum number of atmospheric profiles
        , nd_band                                                       &
!           Maximum number of spectral bands
        , nd_brdf_basis_fnc                                             &
!           Maximum number of BRDF basis functions
        , nd_brdf_trunc                                                 &
!           Maximum order of truncation of BRDF
        , nd_point_tile                                                 &
!           Size allocated for points where the surface is tiled
        , nd_tile
!           Size allocated for surface tiles

      INTEGER, INTENT(INOUT) ::                                         &
          ierr
!           Error flag

!     Actual sizes used:
      INTEGER, INTENT(IN) ::                                            &
          n_points                                                      &
!           Number of atmospheric points
        , n_band
!           Number of spectral bands

!     Points to be treated:
      INTEGER, INTENT(IN) ::                                            &
          list(nd_field)
!           List of sunlit points

      LOGICAL, INTENT(IN) ::                                            &
          land(nd_field)
!           Land flag
      LOGICAL, INTENT(IN) :: l_ctile
!           Flag for coastal tiling
      REAL, INTENT(IN) ::                                               &
          flandg(nd_field)                                              &    
!           land fraction in grid box    
        , emis_land(nd_field)    
!           Mean land emissivity in a gridbox
      REAL, INTENT(IN) ::                                               &
          ice_fraction(nd_field)                                        &
!           Ice fractions in oceanic grid-boxes
        , t_rad_land(nd_field)                                          &
!           Effective radiative temperature of land
        , t_rad_sice(nd_field)                                          &
!           Effective radiative temperature of sea-ice
        , t_rad_sea(nd_field)
!           Effective radiative temperature of sea

!     Physical properties of surfaces:
      INTEGER, INTENT(IN) ::                                            &
          ls_brdf_trunc
!           Order of truncation applied to BRDFs

!     Surface properties set.
      INTEGER, INTENT(OUT) ::                                           &
          n_brdf_basis_fnc
!           Number of basis functions for BRDFs
      REAL, INTENT(OUT) ::                                              &
          f_brdf(nd_brdf_basis_fnc, 0: nd_brdf_trunc/2                  &
            , 0: nd_brdf_trunc/2, 0: nd_brdf_trunc)                     &
!           Basis functions for the surface
        , rho_alb(nd_profile, nd_brdf_basis_fnc, nd_band)
!           Weights of the surface BRDF basis functions

!     Arrays related to tiling of the surface
      LOGICAL, INTENT(IN) ::                                            &
          l_rad_tile
!           Local to allow tiling options
      INTEGER, INTENT(OUT) ::                                           &
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
      REAL, INTENT(OUT) ::                                              &
          rho_alb_tile(nd_point_tile, nd_brdf_basis_fnc                 &
            , nd_tile, nd_band)                                         &
!           Weights for the basis functions of the BRDFs
!           at the tiled points
        , t_tile(nd_point_tile, nd_tile)                                &
!           Local surface temperatures on individual tiles
        , frac_tile(nd_point_tile, nd_tile)
!           Fractions of each tiled grid-point occupied by tiles
!           of the appropriate type


!     Local variables.

      INTEGER ::                                                        &
          i                                                             &
!           Loop variable
        , ib                                                            &
!           Loop variable
        , l                                                             &
!           Loop variable
        , ll
!           Loop variable

      CHARACTER (LEN=*), PARAMETER :: RoutineName = 'r2_set_surface_field_lw'
      CHARACTER (LEN=240) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


      IF (lhook) CALL dr_hook('R2_SET_SURFACE_FIELD_LW',zhook_in,zhook_handle)

!     Set the radiative characteristics of the surface.
    
!     Zero the irrelevant direct albedo.  
      rho_alb(1:n_points, ip_surf_alb_dir, 1:n_band) = 0.0    
  
      DO ib=1,n_band    
        DO i=1,n_points    
          rho_alb(i, ip_surf_alb_diff, ib)                            &    
            = flandg(list(i)) * (1.0 - emis_land(list(i))) +          &  
              (1.0 - flandg(list(i))) * (                             &  
              (1.0 - ice_fraction(list(i))) * (1.0 - emis_sea) +      &  
              ice_fraction(list(i)) * (1.0 - emis_sice)               &  
              )  
        END DO    
      END DO    

!     Set the surface basis functions for a Lambertian surface.
      n_brdf_basis_fnc=1
!     By defining F_{1,0,0,0} to be 4, RHO_ALB beomes equal to the
!     diffuse albedo.
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
          DO ll=1, n_points
            l=list(ll)
            IF ( (flandg(l) < 1.0) .AND.                                &
                 ( (flandg(l) > 0.0) .OR.                               &
                   ( (ice_fraction(l) > 0.0) .AND.                      &
                     (ice_fraction(l) < 1.0) ) ) ) THEN
              n_point_tile=n_point_tile+1
              list_tile(n_point_tile)=ll
              list_tile_outer(n_point_tile)=l
!             Assign tiled fractions consistent with the indices above.
              frac_tile(n_point_tile, 1)                                &
                =(1.0-flandg(l))*(1.0-ice_fraction(l))
              frac_tile(n_point_tile, 2)                                &
                =(1.0-flandg(l))*ice_fraction(l)
              frac_tile(n_point_tile, 3)                                &
                =flandg(l)
            END IF
          END DO

        ELSE

!         Without coastal tiling we have only open sea or sea ice
!         forming the coastal tiling.

          n_tile=2
          index_tile(ip_ocean_tile)=1
          index_tile(ip_seaice_tile)=2
          n_point_tile=0
          DO ll=1, n_points
            l=list(ll)
            IF ( (.NOT.land(l)).AND.                                    &
                 (ice_fraction(l) >  0.0).AND.                          &
                 (ice_fraction(l) <  1.0) ) THEN
              n_point_tile=n_point_tile+1
              list_tile(n_point_tile)=ll
              list_tile_outer(n_point_tile)=l
              frac_tile(n_point_tile, 1)                                &
                =(1.0-ice_fraction(l))
              frac_tile(n_point_tile, 2)                                &
                =ice_fraction(l)
            END IF
          END DO

        END IF

!  
!       Now assign the tiled surface properties at points where tiling  
!       is active. Open sea and sea ice always need to be set with  
!       radiative tiling, but specific land fields are required only  
!       if coastal tiling is on.  
  
!       Zero the irrelevant direct albedos
        rho_alb_tile(1:n_point_tile, ip_surf_alb_dir                    &
          , 1:n_tile, 1:n_band) = 0.0
 
!       The oceanic surface.
        rho_alb_tile(1:n_point_tile, ip_surf_alb_diff                   &
          , ip_ocean_tile, 1:n_band) = (1.0 - emis_sea)  
  
!       Sea-ice.  
        rho_alb_tile(1:n_point_tile, ip_surf_alb_diff                   &  
          , ip_seaice_tile, 1:n_band) = (1.0 - emis_sice)  
  
!       Land points. The test on l_ctile is required only
!       because without coastal tiling land points will be wholly land
!       and do not require radiative tiling.
        IF ( l_ctile ) THEN    
          DO i = 1, n_point_tile    
            rho_alb_tile(i, ip_surf_alb_diff, ip_land_tile,             &    
              1:n_band) = 1.0 - emis_land(list_tile_outer(i))    
          END DO    
        END IF

!       Tiled temperatures are required to handle emission from the
!       surface. Ensure that the indexing of the second subscript is
!       consistent with the assignment of index_tile above.
        t_tile(1:n_point_tile, 1)                                       &
          = t_rad_sea(list_tile_outer(1:n_point_tile))
        IF (l_ctile) THEN
          t_tile(1:n_point_tile, 2)                                     &
            = t_rad_sice(list_tile_outer(1:n_point_tile))
          t_tile(1:n_point_tile, 3)                                     &
            = t_rad_land(list_tile_outer(1:n_point_tile))
        ELSE
!         Without coastal tiling, the solid part of the grid-box 
!         can only be sea ice.
          t_tile(1:n_point_tile, 2)                                     &
            = t_rad_sice(list_tile_outer(1:n_point_tile))
        END IF

      END IF

 9999 CONTINUE
! Check error condition
      IF (ierr /= i_normal) THEN

        CALL ereport(RoutineName, ierr, cmessage)
      END IF

      IF (lhook) CALL dr_hook('R2_SET_SURFACE_FIELD_LW',zhook_out,zhook_handle)
      END SUBROUTINE r2_set_surface_field_lw
