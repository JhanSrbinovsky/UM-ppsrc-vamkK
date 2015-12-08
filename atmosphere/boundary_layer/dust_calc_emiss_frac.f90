! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description:
!   Calculates the dust emission fraction for all land surface tiles, or 
!   provides an aggregate value in the case of a single tile scheme.
!
! Method:
!   For a given input tile it calculates the bare soil emission fraction
!   for that tile, using the required inputs. The different methods for 
!   doing this calculation are dependent on the dust_veg_emiss switch:
!   0: Do not allow dust emission from tiles other than the bare soil tile
!   1: Allow emission from all vegetated tiles, according to their bare soil
!      radiative fraction (used to calculate the surface albedo)
!   2: As 1, but does not allow dust emission from trees
!   It is expected that other options will be implemented as and when the
!   dust emission from vegetation is better constrained by obsevations and
!   more sophisticated models become available.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: 3: Boundary Layer
!
! Code Description:
!   Language: Fortran 90.
!   Programming standard : unified model documentation paper No 3
SUBROUTINE dust_calc_emiss_frac(land_pts, ntiles, tile_pts, tile_index, &
    frac, lai_ft, dust_emiss_frac)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE nstypes, ONLY: npft, ntype, soil
USE pftparm, ONLY: dust_veg_scj
USE dust_parameters_mod, ONLY: dust_veg_emiss

IMPLICIT NONE

! Subroutine arguments  
! IN arguments
  INTEGER, INTENT(IN) :: land_pts       ! No. land points in whole grid
  INTEGER, INTENT(IN) :: ntiles         ! No. land-surface tile types
  INTEGER, INTENT(IN) :: tile_pts(ntype)! Total number of tiles
  INTEGER, INTENT(IN) :: tile_index(land_pts,ntype)
                                      ! Index of tiles on landpts
  REAL,    INTENT(IN) :: frac(land_pts,ntype)
                                      ! IN Fractions of surface types.
  REAL,    INTENT(IN) :: lai_ft(land_pts,npft)
                                   ! IN Leaf area index
! OUT arguments
  REAL, INTENT(OUT):: dust_emiss_frac(land_pts,ntiles)

! Local variables
  INTEGER :: l                          !index of land pt
  INTEGER :: m                          !loop counter, tiles or tile types
  INTEGER :: n                          !loop counter, tile points
  REAL    :: rfracbs(land_pts,npft)     !raditive fraction of bare soil
  REAL    :: dust_veg_sc(npft)          !scaling factor for each PFT

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header

IF (lhook) CALL dr_hook('dust_calc_emiss_frac',zhook_in,zhook_handle)

! explicitly initialise the dust_emiss_frac array as all zeroes as it will be 
! calculated as a sum over tiles
dust_emiss_frac(:,:) = 0.0
! and test to see if this does anything...
rfracbs(:,:) = 0.0

! start with the NULL option, where the emission is only from the bare soil
! tile, and none from vegetation
IF (dust_veg_emiss == 0) THEN

   IF (ntiles == 1) THEN
      m = soil
      DO n = 1, tile_pts(m)
         l = tile_index(n,m)
         dust_emiss_frac(l,1) = frac(l,soil)
      END DO

   ELSE
      m = soil
      DO n = 1, tile_pts(m)
         l = tile_index(n,m)
         dust_emiss_frac(l,m) = frac(l,soil)
      END DO
     
   END IF

ELSE IF (dust_veg_emiss == 1) THEN
! Option 1 is to include the bare soil radiative fraciton of bare soil
! of all vegetated tiles, using the leaf area index:

! Get the radiative fraction of bare soil within each of the vegetated tiles:
! assuming that the vegetated ones are always first (which they should be)
   DO m=1,npft
      DO n=1,tile_pts(m)
         l = tile_index(n,m)
         rfracbs(l,m) = EXP(-lai_ft(l,m)/2.0)
      END DO
   END DO

! get the sclaing parameters for each PFT, for JULES or non-JULES:
   dust_veg_sc(:)=dust_veg_scj(:)

! For 1 tile model, the output is an aggregated emission fraction of all tiles
   IF (ntiles == 1) THEN
! start with the bare soil fraction, where there are tile_pts for that fraction:
      m = soil
      DO n = 1, tile_pts(m)
         l = tile_index(n,m)
         dust_emiss_frac(l,1) = frac(l,m)
      END DO

! Add on the radiative fractions of bare soil over the vegetated tiles
      DO m = 1, npft
         DO n = 1, tile_pts(m)
            l = tile_index(n,m)
            ! output on 1 tile, so (l,1)
            dust_emiss_frac(l,1) = (dust_veg_sc(m)*frac(l,m)*rfracbs(l,m)) +&
                                   dust_emiss_frac(l,1)
         END DO
      END DO
! ends the 1 tile version
   ELSE
! For the multiple tile code, basically the same but each tile has an array:

! start with the vegetated tiles (assumes they are first in the tile arrays)
      DO m = 1, npft
         DO n = 1, tile_pts(m)
            l = tile_index(n,m)
            dust_emiss_frac(l,m) = dust_veg_sc(m)*frac(l,m)*rfracbs(l,m)
         END DO
      END DO

! then do the bare soil tile
      m = soil
      DO n = 1, tile_pts(m)
         l = tile_index(n,m)
         dust_emiss_frac(l,m) = frac(l,m)
      END DO

! other tiles (urban, lake, glacier) are left as zero, so do not emit dust

   END IF ! ends 1 tile or multiple tile test

END IF ! Ends test on dust_veg_emiss switch (currently 0 or 1).
! Additional values of dust_veg_emiss could be added later, using more 
! sophisticated models of emission from partial vegetation.

IF (lhook) CALL dr_hook('dust_calc_emiss_frac',zhook_out,zhook_handle)

RETURN
END SUBROUTINE dust_calc_emiss_frac
