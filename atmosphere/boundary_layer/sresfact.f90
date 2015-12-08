! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!----------------------------------------------------------------------
!
!  Purpose: Calculate gridbox-mean resistance factor used by
!           BL_TRMIX_DD to calculate dry deposition of tracers.
!           For 8A (MOSES II tiled land surface) bdy lyr versions,
!           adapted from RESTILE8A.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Boundary Layer
!
! Code description:
!   Language: FORTRAN 77  + common extensions
!   This code is written to UMDP3 v6 programming standards.
!
! System Components covered:
!
! System task:
!
! Documentation:  UMDP 20, HCTN 30
!
!----------------------------------------------------------------------
SUBROUTINE sresfact (land_pts, land_index,                              &
                     ntiles, tile_index, tile_pts, soluble,             &
                     canopy, catch, gs_tile, tile_frac, snow_tile,      &
                     aresist, aresist_tile, resb, ress,                 &
                     resist_b_tile, res_factor_land)

  USE atm_fields_bounds_mod, ONLY: tdims
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  IMPLICIT NONE

! Arguments with intent IN:

  INTEGER, INTENT(IN) ::                                                &
   land_pts,                                                            &
                                    !Total number of land points
   land_index(land_pts),                                                &
                                    !Index of land points.
! For MOSES II
   ntiles,                                                              &
                                    !Number of land tiles.
   tile_index(land_pts,ntiles),                                         &
                                    !Index of tile points.
   tile_pts(ntiles)             !Number of tile points.

  LOGICAL, INTENT(IN) ::                                                &
   soluble                      !.TRUE. for soluble aerosols

  REAL, INTENT(IN) ::                                                   &
   aresist(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                                    !Aerodynamic resistance Ra (s/m)
   resb,                                                                &
                                    !Rb(aerosol) / Rb(H2O).
   ress,                                                                &
                                    !Rs(aerosol) / Rs(H2O).
! For MOSES II
   aresist_tile(land_pts,ntiles),                                       &
                                    ! 1/(CD_STD*VSHR) on land tiles.
   canopy(land_pts,ntiles),                                             &
                                    !Surface water on land tiles (kg/m2)
   catch(land_pts,ntiles),                                              &
                                    !Surface capacity (max. surface
                                    !    water) of land tiles (kg/m2)
   gs_tile(land_pts,ntiles),                                            &
                                    !Surface conductance for land tiles
   snow_tile(land_pts,ntiles),                                          &
                                    !Snow mass on land tiles (kg/m2).
   tile_frac(land_pts,ntiles)       !Tile fractions.

! Arguments with intent INOUT
  REAL, INTENT(INOUT) ::                                                &
   resist_b_tile(land_pts,ntiles)
                                    !(1/CH-1/CD_STD)/VSHR on land tiles.

! Arguments with intent OUT:
  REAL, INTENT(OUT) ::                                                  &
   res_factor_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                    !Ra/(Ra+Rb+Rs) mean over land
!                                         part of grid box

! Local variables:

  INTEGER ::                                                            &
   i, j,                                                                &
                   !Loop counters (horizontal field index).
   k,                                                                   &
                   !Loop counter (tile field index).
   l,                                                                   &
                   !Loop counter (land point field index).
   n           !Loop counter (tile index).

  REAL ::                                                               &
   damp_factor(land_pts,ntiles),                                        &
                                     !Canopy moistening factor
   rs_tile(land_pts,ntiles),                                            &
                                     !Surface resistance for land tiles

   str_resist_b,                                                        &
                                     !Rb for aerosol.
   str_resist_s,                                                        &
                                     !Rs for aerosol.
   snow_f                        !Snow cover fraction.

! Parameters:
  REAL ::                                                               &
   asnow,                                                               &
                                     !Parameter for snow fraction calcn
   cond_lim,                                                            &
                                     !Low limit for canopy conductance.
   r_snow                        !Resistance to dry dep. over snow

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

  PARAMETER (asnow=0.2,                                                 &
             cond_lim=1.0e-3,                                           &
             r_snow=1.0e3)

!   Initialise DAMP_FACTOR to 1.0
!   Note that for RESIST_B_TILE, negative values have already been
!   set to in SFEXCH8A, but we repeat it here for safety.
!   For STR_RESIST_S values depend on surface type (land, sea,
!   snow, ice) as well as tracer identity. First calculate stomatal
!   resistance (=1/conductance, avoiding dividing by zero).

  IF (lhook) CALL dr_hook('SRESFACT',zhook_in,zhook_handle)
  DO n = 1, ntiles
    DO k = 1, tile_pts(n)
      l = tile_index(k,n)
      damp_factor(l,n) = 1.0
      IF (resist_b_tile(l,n)  <   0.0)  THEN
        resist_b_tile(l,n) = 0.0
      END IF
      IF (gs_tile(l,n)  >   cond_lim) THEN
        rs_tile(l,n) = 1. / gs_tile(l,n)
      ELSE
        rs_tile(l,n) = 1. / cond_lim
      END IF
    END DO
  END DO

!  For SOLUBLE species (SO2 and NH3) reduce the surface resistance by up
!  to two-thirds if the canopy is damp (the value of 2/3 is empirical).
!  Two special cases need to be trapped here. The canopy capacity
!  (CATCH) is zero at land ice points, so exclude these from the
!  calculation. Also, there is a possibility that canopy water may
!  exceed canopy capacity due to leaves having fallen, so take care
!  of this too.

  IF (soluble) THEN
! Loop over all land tiles:
    DO n = 1, ntiles
      DO k = 1, tile_pts(n)
        l = tile_index(k,n)

        IF( (catch(l,n)  >   0.01) .AND.                                &
            (canopy(l,n)  >   0.0) ) THEN
          IF( canopy(l,n)  <=  catch(l,n) ) THEN
            damp_factor(l,n) = 1. - 0.66667*canopy(l,n)/catch(l,n)
          ELSE
            damp_factor(l,n) = 0.33333
          END IF
        END IF

      END DO
    END DO

  END IF

! Initialise RES_FACTOR_LAND to zero:

  DO j=tdims%j_start,tdims%j_end
    DO i=tdims%i_start,tdims%i_end
      res_factor_land(i,j) = 0.
    END DO
  END DO

!  Need to set STR_RESIST_S to suitable values over snow and ice.
!  Where there is snow cover, calculate an approximate snow fraction
!  for the tile using the formula 1-exp(-ASNOW*SNODEP)
!  Note that for atmospheric model run there should not be any sea
!  points with SNODEP >  0, and land_ice points (including Antarctica)
!  should all have large values of SNOW_TILE.

!  Note that the value of ARESIST used here in the calculation of
!  RES_FACTOR_LAND must be the same as that incorporated in RHO_ARESIST
!  (i.e. a grid-box mean including land-sea averaging) used in TR_MIX.

  DO n = 1, ntiles
    DO k = 1, tile_pts(n)
      l = tile_index(k,n)
      j = (land_index(l)-1)/tdims%i_end + 1
      i = land_index(l) - (j-1)*tdims%i_end
! Loop over all land tiles:
!  (Note: Routine TILEPTS sets array elements TILE_PTS(N) to the no. of
!   gridboxes including surface type N, and TILE_INDEX(K,N) to the land
!   array index of the k-th gridbox containing surface type N.
!   See HCTN 30 p 25.)

      str_resist_b = resb*resist_b_tile(l,n)
      str_resist_s = ress*rs_tile(l,n)*damp_factor(l,n)

      IF (snow_tile(l,n) >  0.0 .AND.                                   &
                             str_resist_s >  0.0) THEN
        snow_f = 1. - EXP(-asnow*snow_tile(l,n))
        str_resist_s = 1./                                              &
            (snow_f/r_snow + (1.-snow_f)/str_resist_s)
      END IF

      res_factor_land(i,j) = res_factor_land(i,j) +                     &
                             aresist(i,j)*tile_frac(l,n) /              &
                     (aresist_tile(l,n)+str_resist_b+str_resist_s)

    END DO
  END DO

  IF (lhook) CALL dr_hook('SRESFACT',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE sresfact
