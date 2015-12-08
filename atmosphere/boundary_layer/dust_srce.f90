! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Calculates flux of mineral dust entrained into atmosphere
!
! Method:
!   Calculates mineral dust flux as a function of friction velocity (U*),
!   soil moisture and soil porperties developed from description in:
!     "Modelling the atmospheric lifecycle..."
!     Woodward, JGR106, D16, pp18155-18166, 2001.
!   Treatment of soil moisture as in:
!     "Parametrization of the increase of..."
!     Fecan et al, Ann. Geophysicae 17, 149-157, 1999.!
!   NB: Currently only calculates flux from bare soil tiles
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!
! Code description:
!  Language: Fortran 90
!  Programming standard : unified model documentation paper No 3
!
! Subroutine Interface:
SUBROUTINE dust_srce(                                                   &
! IN arguments
     land_pts, ntiles, sm_levels, tile_pts, tile_index,                 &
     fland, tstar_tile, rhostar_land, soil_layer_moisture, snow_tile,   &
     u_s_std_tile, mrel_land, clay_land, sand_land, ho2r2_orog,         &
! OUT arguments
     dust_flux_tile, u_s_t_tile, u_s_t_dry_tile )

  USE nstypes, ONLY: ntype, soil
  USE earth_constants_mod, ONLY: g
  USE conversions_mod, ONLY: zerodegc
  USE dust_parameters_mod, ONLY:                                        &
       ! number of divisions that can be lifted from the surface and the
       ! number that can be blown horizontally along the surface and 
       ! contribute to the lifting of the 1st NDIV divisions and the 
       ! number of these that contain data in the dust soil ancillaries
       ndiv, ndivh, ndivl,                                              &
       ! impact U*t derived from Bagnold (1941)
       ustd_bas,                                                        &
       ! additional and multiplicative tunings to U* and U*t
       us_aa, us_am, ust_aa, ust_am,                                    &
       ! multiplicative tuning factor to level 1 soil moisture 
       sm_corr,                                                         &
       ! maximum clay fraction used in calculation of vertical flux
       clay_max,                                                        & 
       ! various limits used in diagnosis of horizontal flux
       u_s_min, snowmin, h_orog_limit, fland_lim,                       &
       ! constants used in horiz. and vert. fluxes (see module for details)
       horiz_c, horiz_d, vert_a, vert_b, vert_c,                        &
       ! switch to diagnose vertical flux using a fixed size distribution
       ! if off the vertical flux is proportional to the horizontal flux
       l_fix_size_dist,                                                 &
       ! proportion of flux from each bin if using the fixed distribution
       size_dist,                                                       &
       ! Two-bin dust flag
       l_twobin_dust
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

  ! Subroutine arguments
  INTEGER, INTENT(IN) :: land_pts       ! No. land points in whole grid
  INTEGER, INTENT(IN) :: ntiles         ! No. land-surface tile types
  INTEGER, INTENT(IN) :: sm_levels      ! No. soil levels
  INTEGER, INTENT(IN) :: tile_pts(ntype)! Total number of tiles
  INTEGER, INTENT(IN) :: tile_index(land_pts,ntype)
                                      ! Index of tiles on landpts
  REAL, INTENT(IN) :: fland(land_pts)   ! Land fraction on land points
  REAL, INTENT(IN) :: tstar_tile(land_pts,ntiles)
                                      ! Surface temp on tiles
  REAL, INTENT(IN) :: rhostar_land(land_pts)
                                      ! Surface air density on land pts
  REAL, INTENT(IN) :: soil_layer_moisture(land_pts,sm_levels)
                                      ! Soil moisture (kg m-2)
  REAL, INTENT(IN) :: snow_tile(land_pts,ntiles)
                                      ! Lying snow on tiles (kg m-2)
  REAL, INTENT(IN) :: u_s_std_tile(land_pts,ntiles)
                                      ! Friction velocity on tiles
  REAL, INTENT(IN) :: mrel_land(land_pts,ndivl)
                                      ! Relative soil mass per size div
  REAL, INTENT(IN) :: clay_land(land_pts)
                                      ! Soil clay fraction
  REAL, INTENT(IN) :: sand_land(land_pts)
                                      ! Soil sand fraction
  REAL, INTENT(IN) :: ho2r2_orog(land_pts)
                                      ! Half peak to trough ht of orog
! Outputs
  REAL, INTENT(OUT):: dust_flux_tile(land_pts,ntiles,ndiv)
                                      ! Dust flux (kg m-2 s-1)
  REAL, INTENT(OUT):: u_s_t_tile(land_pts,ntiles,ndivh)
                                      ! Thresh. friction vel. per tile
                                      ! (all 9 divisions)
  REAL, INTENT(OUT):: u_s_t_dry_tile(land_pts,ntiles,ndivh)
                                      ! Thresh. frict. vel. per tile
                                      ! excluding soil moisture effects
! Local variables
  INTEGER :: l                          !index of land pt
  INTEGER :: m                          !loop counter, tile types
  INTEGER :: n                          !loop counter, tile points
  INTEGER :: idiv                       !loop counter, dust divisions

  REAL :: ratio                         ! u_s_t/u_s_std
  REAL :: mrel7(land_pts)               ! mrel for 31.6 - 100 um radius
  REAL :: mrel8(land_pts)               ! mrel for 100 - 316 um
  REAL :: mrel9(land_pts)               ! mrel for 316 - 1000 um
  REAL :: mrel_land_all(land_pts,ndivh) !all the mrels


  REAL :: horiz_flux_tile(land_pts,ntiles,ndivh)
                                      !horizontal flux per tile
  REAL :: tot_horiz_flux_tile(land_pts,ntiles)
                                      !total over all divs
  REAL :: horiz_flux_789(land_pts,ntiles)
                                      !total over div 7,8,9
  REAL :: horiz_flux_1to6(land_pts,ntiles)
                                      !total over divs 1 to 6
  REAL :: us(land_pts)                  ! "tuned" U* on soil tiles
  REAL :: smt(land_pts)                 ! Soil moisture term in U*t calc
  REAL :: smt_work                      ! Temp variable in smt calc

  INTEGER :: soil_tile                     ! Tile on which to calculate
                                           ! bare soil dust flux

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  IF (lhook) CALL dr_hook('DUST_SRCE',zhook_in,zhook_handle)

! Initialisation
  DO idiv = 1, ndiv
    DO m = 1, ntiles
      DO n = 1, tile_pts(m)
        l = tile_index(n,m)
        dust_flux_tile(l,m,idiv) = 0.
      END DO !tile_pts
    END DO !ntiles
  END DO !ndiv

  DO m = 1, ntiles
    DO n = 1, tile_pts(m)
      l = tile_index(n,m)
      tot_horiz_flux_tile(l,m) = 0.
      horiz_flux_789(l,m) = 0.
      horiz_flux_1to6(l,m) = 0.
    END DO
  END DO

  DO idiv = 1, ndivh
    DO m = 1, ntiles
      DO n = 1, tile_pts(m)
        l = tile_index(n,m)
        u_s_t_tile(l,m,idiv) = 0.
        u_s_t_dry_tile(l,m,idiv) = 0.
        horiz_flux_tile(l,m,idiv) = 0.
      END DO
    END DO
  END DO

! In 1 bin scheme, bare soil dust flux is calculated on single aggregated tile
  IF (ntiles > 1) THEN
    soil_tile = soil
  ELSE
    soil_tile = 1
  END IF

! Calculate relative mass for all divs and put into single array
  DO l = 1, land_pts
    mrel7(l)=sand_land(l)*0.312
    mrel8(l)=sand_land(l)*0.312
    mrel9(l)=sand_land(l)*0.312
    mrel_land_all(l,7)=mrel7(l)
    mrel_land_all(l,8)=mrel8(l)
    mrel_land_all(l,9)=mrel9(l)
  END DO !land_pts
  
! This is explicitly specified as ndivl=6, even for the two-bin scheme  
  DO idiv = 1, ndivl
    DO l = 1, land_pts
      mrel_land_all(l,idiv)=mrel_land(l,idiv)
    END DO
  END DO

! Increase U*t for "tuned" soil moisture (w) > a threshold w' (fecan'99)
! smt_work = w-w'
! w is % volmetric sm, but approximate as same as s_l_m in top moses level
! ( w = s_l_m * 100(for %age) * (100cm/10cm) / rho=998 ~ s_l_m)
  DO l = 1, land_pts
    smt_work = soil_layer_moisture(l,1) * sm_corr -                     &
         14.*clay_land(l)*clay_land(l) - 17.*clay_land(l)
    IF (smt_work > 0.) THEN
      smt(l) =  ((1.0 + 1.21* smt_work **.68 )**.5)
    ELSE
      smt(l) = 1.0
    END IF
  END DO

! Loop over points in each div calculating dust flux, if any
  DO idiv = 1, ndivh
    DO m = 1, ntiles
      DO n = 1, tile_pts(m)
        l = tile_index(n,m)
        !
        ! Calculate "tuned" U* and/or U*t for this tile
        !
        u_s_t_dry_tile(l,m,idiv)=ustd_bas(idiv)*ust_am+ust_aa
        u_s_t_tile(l,m,idiv)=ustd_bas(idiv)*smt(l)*ust_am+ust_aa

        ! use the U* on bare soil, regardless of whether this is a bare
        ! soil tile or not - we're looking at the U* of any bare soil in
        ! the tile, the fraction of which is determined by 
        ! dust_calc_emiss_frac
        ! For the 1 tile case, u_s_std_tile is calculated as bare soil only
        us(l)=u_s_std_tile(l,soil_tile)*us_am+us_aa

        !
        ! Horizontal flux
        !
        IF ( (u_s_t_tile(l,m,idiv) < us(l)) .AND.                         &
             (fland(l) >  fland_lim) .AND.                                & 
             (tstar_tile(l,m) > zerodegc) .AND.                           &
             (snow_tile(l,m)  <   snowmin) .AND.                          &
             (mrel_land_all(l,idiv)  >   0.) .AND.                        &
             (ho2r2_orog(l)  <=  h_orog_limit) .AND.                      &
             ((soil_layer_moisture(l,1)*sm_corr <                         &
             (clay_land(l)+.12)/.03)) .AND.                               &
             (us(l) > u_s_min) ) THEN

          ratio = u_s_t_tile(l,m,idiv) / us(l)
          horiz_flux_tile(l,m,idiv) = horiz_c * rhostar_land(l) *         &
               horiz_d * us(l)**3. * (1.0 + ratio) *                      &
               (1.0 - ratio*ratio) * mrel_land_all(l,idiv) / g
          tot_horiz_flux_tile(l,m)=                                       &
               tot_horiz_flux_tile(l,m)+horiz_flux_tile(l,m,idiv)

        END IF
      END DO !tile_pts
    END DO ! tiles
  END DO !ndiv

  DO m = 1, ntiles
    DO n = 1, tile_pts(m)
      l = tile_index(n,m)
      horiz_flux_789(l,m)=horiz_flux_tile(l,m,7)+                         &
           horiz_flux_tile(l,m,8)+horiz_flux_tile(l,m,9)
      horiz_flux_1to6(l,m)=horiz_flux_tile(l,m,1)+horiz_flux_tile(l,m,2)+ &
           horiz_flux_tile(l,m,3)+horiz_flux_tile(l,m,4)+                 &
           horiz_flux_tile(l,m,5)+horiz_flux_tile(l,m,6)
    END DO
  END DO

! Vertical flux
! Again, this calculation is for bare soil only, but could be
! extended to other tile types
  IF (l_fix_size_dist.OR.l_twobin_dust) THEN
!   Lift dust using fixed size distribution
    DO idiv = 1, ndiv
      DO m = 1, ntiles
        DO n = 1, tile_pts(m)
          l = tile_index(n,m)
          IF (snow_tile(l,m) < snowmin) THEN
!           Calculate proportion of vertical flux to come from each bin
            dust_flux_tile(l,m,idiv) = vert_c * size_dist(idiv) *         &
                 (horiz_flux_1to6(l,m)+horiz_flux_789(l,m))          
!           Calculate vertical flux
            dust_flux_tile(l,m,idiv) =  dust_flux_tile(l,m,idiv) *        &
               10**(vert_a * MIN(clay_land(l),clay_max) + vert_b)
          END IF
        END DO
      END DO
    END DO
  ELSE
!   Lift dust using diagnosed size distribution
    DO idiv = 1, ndiv
      DO m = 1, ntiles
        DO n = 1, tile_pts(m)
          l = tile_index(n,m)
!         Calculate proportion of vertical flux to come from each bin
          IF ((snow_tile(l,m) < snowmin) .AND.                             &
              (mrel_land(l,idiv) > 0.) .AND.                               &
              (horiz_flux_1to6(l,m) > 0.)) THEN
            dust_flux_tile(l,m,idiv) = vert_c *                            &
                 (1.0 + horiz_flux_789(l,m)/horiz_flux_1to6(l,m)) *        &
                 horiz_flux_tile(l,m,idiv)
          END IF
!         Calculate vertical flux
          dust_flux_tile(l,m,idiv) =  dust_flux_tile(l,m,idiv) *           &
             10**(vert_a * MIN(clay_land(l),clay_max) + vert_b)
        END DO
      END DO
    END DO
  END IF

! Correct fluxes for land fraction - for coastal points
  DO idiv = 1, ndiv
    DO m = 1, ntiles
      DO n = 1, tile_pts(m)
        l = tile_index(n,m)
      ! Horizontal flux not yet output as a diagnositic
      ! If added, do the following calc (but for 9 divisions)
!     horiz_flux_tile(l,m,idiv)=horiz_flux_tile(l,m,idiv)*fland(l)
        dust_flux_tile(l,m,idiv) = dust_flux_tile(l,m,idiv)*fland(l)
      END DO !tile_pts
    END DO !ntiles
  END DO !ndiv

  IF (lhook) CALL dr_hook('DUST_SRCE',zhook_out,zhook_handle)
  RETURN

END SUBROUTINE dust_srce
