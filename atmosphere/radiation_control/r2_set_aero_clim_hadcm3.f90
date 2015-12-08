! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to set fields of climatological aerosols in HADCM3.

! Purpose:
!   This routine sets the mixing ratios of climatological aerosols.
!   A separate subroutine is used to ensure that the mixing ratios
!   of these aerosols are bit-comparable with earlier versions of
!   the model where the choice of aerosols was more restricted:
!   keeping the code in its original form reduces the opportunity
!   for optimizations which compromise bit-reproducibilty.
!   The climatoogy used here is the one devised for HADCM3.

! Method:
!   Straightforward.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!- ---------------------------------------------------------------------
      SUBROUTINE r2_set_aero_clim_hadcm3(n_profile, nlevs, n_layer      &
         , i_gather, l_extra_top                                        &
         , l_clim_aero_hgt, bl_depth, t, n_levels_bl                    &
         , land, lying_snow, pstar, p_layer_boundaries, trindx          &
         , alat, previous_time                                          &
         , aerosol_mix_ratio_clim                                       &
         , npd_field, npd_profile, npd_layer, first_layer               &
         )

      USE rad_input_mod, ONLY:                                          &
                      l_rad_use_clim_volc,                              &
                      clim_rad_volc_eruption_year,                      &
                      clim_rad_volc_eruption_month,                     &
                      clim_rad_volc_eruption_weight

      USE earth_constants_mod, ONLY: g

      USE atmos_constants_mod, ONLY: r

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE


!     COMDECKS INCLUDED.

!     DUMMY ARGUMENTS.

!     SIZES OF ARRAYS:
      INTEGER                                                           &
                !, INTENT(IN)
           npd_field                                                    &
!             FIELD SIZE IN CALLING PROGRAM
         , npd_profile                                                  &
!             SIZE OF ARRAY OF PROFILES
         , npd_layer                                                    &
!             MAXIMUM NUMBER OF LAYERS
         , first_layer
!             First layer for some variables
!             0 for flux_calc(3A/C), 1 for radiance_calc(3Z)

!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
           n_profile                                                    &
!             NUMBER OF PROFILES
         , nlevs                                                        &
!             Number of atmospheric layers used outside the radiation
!             scheme
         , n_layer
!             Number of atmospheric layers seen in radiation

!     Variables related to options for setting the field
      LOGICAL, INTENT(IN) ::                                            &
           l_clim_aero_hgt                                              &
!             Flag to use the depth of the boundary layer to set
!             the layers in which a boundary-layer aerosol is used
         , l_extra_top
!             Flag to use an extra top layer in radiative calculations
      INTEGER, INTENT(IN) ::                                            &
           n_levels_bl
!             Common number of layers taken to be occupied by the
!             boundary-layer aerosol if this the boundary layer
!             depth is not used to determine the number separately
!             at each grid-point
      REAL, INTENT(IN) ::                                               &
           bl_depth(npd_field)                                          &
!             Depth of the boundary layer
         , t(npd_profile, first_layer: npd_layer)
!             Temperatures of atmospheric layers

!     GATHERING ARRAY:
      INTEGER                                                           &
                !, INTENT(IN)
           i_gather(npd_field)
!             LIST OF POINTS TO GATHER

!     GENERAL ATMOSPHERIC PROPERTIES:
      INTEGER                                                           &
                !, INTENT(IN)
           trindx(npd_field)                                            &
!             LAYER BOUNDARY OF TROPOPAUSE
      ,    previous_time(7)
!             Time

      REAL                                                              &
                !, INTENT(IN)
           pstar(npd_field)                                             &
!             SURFACE PRESSURES
      ,    p_layer_boundaries(npd_field,0:nlevs)                        &
!             PRESSURE AT BOUNDARIES OF LAYERS
      ,    alat(npd_profile)

!     SURFACE FIELDS
      LOGICAL                                                           &
                !, INTENT(IN)
           land(npd_field)
!             LAND-SEA MASK
      REAL                                                              &
                !, INTENT(IN)
           lying_snow(npd_field)
!             DEPTH OF LYING SNOW

      REAL                                                              &
                !, INTENT(OUT)
           aerosol_mix_ratio_clim(npd_profile, first_layer:npd_layer, 5)
!             MIXING RATIOS OF CLIMATOLOGICAL AEROSOLS



!     LOCAL VARIABLES:
      INTEGER                                                           &
           i                                                            &
!             LOOP VARIABLE
         , i_um                                                         &
!             Index of a level in the UM's upward convention
         , j                                                            &
!             LOOP VARIABLE
         , l                                                            &
!             LOOP VARIABLE
         , lg                                                           &
!             INDEX FOR GATHERING
         , bl_top_lyr(npd_profile)                                      &
!             Topmost layer occupied by boundary layer aerosol,
!             indexed using the top-down convention of radiation.
         , n_points_set_z                                               &
!             Number of points where the top of the boundary layer
!             has not yet been reached
         , imon                                                         &
!             Number of months into the volcanic aerosol period
         , j_up_volc                                                    &
!             Top level number of volcanic aerosol layer
         , volcts(192)
!             Array of 0.55 micron optical depths for notional
!             eruption*10000. 4 quarters of globe for 12 months
!             for 4 years.

      LOGICAL                                                           &
           l_set_z(npd_profile)                                         &
!             Array to flag points where the top of the boundary layer
!             has not yet been reached
         , l_volc_per
!             Flag for whether in post-volcanic 4 years
      REAL                                                              &
           pressure_wt(npd_field)                                       &
!             ARRAY FOR SCALING AEROSOL AMOUNTS FOR DIFFERENT SURFACE
!             PRESSURES
         , p_toa(npd_profile)                                           &
!             Pressure at the top of the atmosphere seen by radiation
         , z_bottom(npd_profile)                                        &
!             Height of the bottom of the current layer above
!             the surface
         , dz                                                           &
!             Depth of the current layer
         , mascon
!             Coversion factor from optical depth to mass loading

      REAL, PARAMETER :: volc_up_press = 3000.0
!             Pressure level for top of volcanic aerosol layer
      REAL, PARAMETER :: delta_lat = 1.0e-6
!             Tolerance for checking latitude


!     TOTAL COLUMN MASS (KG M-2) OF EACH AEROSOL SPECIES IN
!     THE BOUNDARY LAYER, THE FREE TROPOSPHERE AND THE STRATOSPHERE
!     RESPECTIVELY. THIS MODEL ASSUMES THAT THERE ARE FIVE AEROSOLS.
      REAL                                                              &
           bl_oceanmass(5)                                              &
         , bl_landmass(5)                                               &
         , freetrop_mass(5)                                             &
         , strat_mass(5)                                                &
         , volc_mass

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!     INITIALIZATION FOR THE CLIMATOLOGICAL AEROSOL MODEL
      DATA bl_landmass/2.77579e-5, 6.70018e-5, 0.0, 9.57169e-7, 0.0/
      DATA bl_oceanmass/1.07535e-5, 0.0, 2.043167e-4, 0.0, 0.0/
      DATA freetrop_mass/3.46974e-6, 8.37523e-6, 0.0, 1.19646e-7, 0.0/
      DATA strat_mass/0.0, 0.0, 0.0, 0.0, 1.86604e-6/
      DATA volcts/90, 188, 177, 66, 108, 471,                           &
       337, 105, 175, 827, 633, 197, 236, 1035, 845, 329, 378, 1211,    &
       1000, 489, 525, 1316, 1133, 645, 663, 1257, 1092, 809, 768, 1317,&
       1235, 938, 869, 1289, 1308, 1049, 989, 1204, 1223, 1035, 1022,   &
       1164, 1172, 1013, 1054, 1126, 1124, 1004, 1074, 1017, 1015, 969, &
       1023, 982, 985, 943, 982, 916, 926, 921, 894, 840, 898, 872, 841,&
       802, 906, 854, 801, 738, 858, 819, 751, 639, 739, 753, 722, 585, &
       672, 712, 679, 526, 611, 650, 649, 473, 547, 600, 621, 430, 502, &
       560, 582, 399, 466, 515, 545, 363, 426, 477, 500, 340, 396, 438, &
       458, 332, 376, 412, 414, 312, 351, 387, 380, 296, 327, 373, 352, &
       279, 305, 346, 324, 252, 272, 317, 301, 233, 253, 291, 288, 215, &
       235, 272, 270, 201, 212, 265, 257, 178, 197, 251, 245, 168, 184, &
       235, 225, 154, 169, 216, 205, 140, 154, 196, 184, 126, 138, 177, &
       164, 112, 123, 157, 143, 98, 108, 137, 123, 84, 92, 118, 102, 70,&
       77, 98, 82, 56, 61, 79, 61, 42, 46, 59, 41, 28, 31, 39, 20, 14,  &
       15, 20, 1, 1, 1, 1 /

      IF (lhook) CALL dr_hook('R2_SET_AERO_CLIM_HADCM3',zhook_in,zhook_handle)

!     TROPOSPHERIC AEROSOL LOADING IS A SIMPLE FUNCTION OF SURFACE
!     PRESSURE: HALVING PSTAR HALVES THE TROPOSPHERIC AEROSOL BURDEN.
!     THE STRATOSPHERIC BURDEN IS INDEPENDENT OF PSTAR.  NOTE THE
!     FACTOR MULTIPLING AEROSOL AMOUNTS USES A REFERENCE PRESSURE
!     OF 1013 mbars.
      DO l=1, n_profile
        pressure_wt(l)=pstar(i_gather(l))*(1.0/1.013e5)
      END DO

!     For each of the 5 aerosol species, the column amount in the
!     boundary layer, free troposphere and stratosphere are known for
!     a standard atmosphere over ocean and land. These can be used
!     to find mixing ratios for the UM by dividing total aerosol by
!     total air mass (and using pressure weighting in the
!     troposphere).

!     Firstly, mixing ratios are set for the 5 aerosol species in the
!     stratosphere.
      IF (l_extra_top) THEN
!       With an extra layer the radiative atmosphere is notionally
!       extended to zero pressure.
        DO l=1, n_profile
          p_toa(l)=0.0e+00
        END DO
      ELSE
!       Otherwise the top of the atmosphere seen elsewhere in the
!       model is used.
        DO l=1, n_profile
          p_toa(l)=p_layer_boundaries(i_gather(l), nlevs)
        END DO
      END IF

!     Loop over aersosol species for stratosphere

      DO i=1,5

!     Set climatology

        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio_clim(l, n_layer+1-trindx(lg), i)            &
            = strat_mass(i)*g/                                          &
              (p_layer_boundaries(lg,trindx(lg)-1)                      &
               - p_toa(l))
        END DO

        DO l=1, n_profile
          lg=i_gather(l)
          DO j=(trindx(lg)+1),n_layer
            aerosol_mix_ratio_clim(l, n_layer+1-j, i)=                  &
              aerosol_mix_ratio_clim(l, n_layer+1-trindx(lg), i)
          END DO
        END DO

        IF (l_rad_use_clim_volc      .AND.                              &
            i  ==  5                     ) THEN

!     Use time-varying aerosol amounts representing an idealised
!     volcanic eruption. The eruption effect lasts 4 years after
!     the date specified by clim_rad_volc_eruption_year and
!     clim_rad_volc_eruption_month. Volcanic amounts are only
!     applied within this period, which is checked by setting
!     l_volc_per:

!     Find out how far the current time is through the period:

          imon = (previous_time(1) - clim_rad_volc_eruption_year)*12 +  &
                 (previous_time(2) - clim_rad_volc_eruption_month)
          IF (imon  >=  0 .AND. imon  <   48) THEN
            l_volc_per = .TRUE.
          ELSE
            l_volc_per = .FALSE.
          END IF

          mascon = 0.76*clim_rad_volc_eruption_weight*0.0001*1.16898/&
                     (6.09715e+3 + 4.60178e-4)
!     0.76 scales the numbers in the VOLCTS so the maximum optical
!     depth is equal to 0.1 (moderate sized eruption).
!     clim_rad_volc_eruption_weight is a multiplier specified by the
!     user based on the size of the eruption. 0.0001 converts VOLCTS to
!     optical depth. 1.16898 converts .55 micron optical depth to the
!     Slingo (1989) band-1 optical depth (from the HadCM2 mod). Division
!     is by the total k (scattering+absorption) for bands 2 & 3,
!     originally from spec3a_sw_3_tknlaero but still the same in
!     spec3a_sw_hadgem1_2.


          IF (l_volc_per) THEN

            DO l=1, n_profile
              lg=i_gather(l)
              IF      (alat(l)  >    30.0) THEN
                volc_mass = REAL(volcts(imon*4 + 1))
              ELSE IF (alat(l)  >     0.0 .AND. alat(l)  <   30.0) THEN
                volc_mass = REAL(volcts(imon*4 + 2))
              ELSE IF (alat(l)  >   -30.0 .AND. alat(l)  <    0.0) THEN
                volc_mass = REAL(volcts(imon*4 + 3))
              ELSE IF (alat(l)  <   -30.0) THEN
                volc_mass = REAL(volcts(imon*4 + 4))
              END IF

              IF      (ABS(alat(l)-30.0)  <   delta_lat) THEN
                volc_mass = 0.5*REAL(volcts(imon*4 + 1) +               &
                                     volcts(imon*4 + 2))
              ELSE IF (ABS(alat(l))       <   delta_lat) THEN
                volc_mass = 0.5*REAL(volcts(imon*4 + 2) +               &
                                     volcts(imon*4 + 3))
              ELSE IF (ABS(alat(l)+30.0)  <   delta_lat) THEN
                volc_mass = 0.5*REAL(volcts(imon*4 + 3) +               &
                                     volcts(imon*4 + 4))
              END IF
              volc_mass = volc_mass*mascon
              IF (volc_mass  <   strat_mass(i)) THEN
                volc_mass = strat_mass(i)
              END IF

!     Locate level of upper limit of volcanic aerosol layer

              DO j=(trindx(lg)+1),n_layer
                IF (p_layer_boundaries(lg,j-1)  >=  volc_up_press .AND. &
                    p_layer_boundaries(lg,j)    <   volc_up_press) THEN
                  j_up_volc = j
                END IF
              END DO

              aerosol_mix_ratio_clim(l, n_layer+1-trindx(lg), i)        &
                  = volc_mass*g/                                        &
                    (p_layer_boundaries(lg,trindx(lg)-1)                &
                     - p_layer_boundaries(lg,j_up_volc))

              DO j=(trindx(lg)+1),j_up_volc
                aerosol_mix_ratio_clim(l, n_layer+1-j, i)=              &
                 aerosol_mix_ratio_clim(l, n_layer+1-trindx(lg), i)
              END DO
            END DO

          END IF
        END IF

      END DO


!      At each point the number of layers considered to contain
!      boundary layer aerosol must be determined. In the original
!      form of the scheme this number is fixed, but in the newer
!      form it is determined from the boundary layer depth.
       IF (l_clim_aero_hgt) THEN

!        Initialize:
         DO l=1, n_profile
           bl_top_lyr(l)=n_layer
           l_set_z(l)=.TRUE.
           z_bottom(l)=0.0e+00
         END DO
         n_points_set_z=n_profile
           i=n_layer

         DO WHILE (n_points_set_z >  0)

!          Assign the UM's indexing over layers: the UM now indexes the
!          boundaries of layers from 0 at the bottom to NLEVS at the
!          top, while the radiation code uses the same range but starts
!          at the top (possibly with an extra layer at the top of
!          the model). I and I_UM refer to layers, not to the edges of
!          layers. BL_TOP_LYR now holds the topmost layer containing
!          boundary layer aerosol indexed as in the radiation code.
           i_um=n_layer+1-i
           DO l=1, n_profile
             IF (l_set_z(l)) THEN
               lg=i_gather(l)
               dz=r*t(l, i)*LOG(p_layer_boundaries(lg, i_um-1)          &
                 /p_layer_boundaries(lg, i_um))/g
               IF ( (MAX(bl_depth(lg), 5.0e+02) >=                      &
                     z_bottom(l)+0.5e+00*dz).AND.                       &
                    (i_um <= trindx(lg)-2) ) THEN
!                The top of the boundary layer is above the middle
!                of this layer, so we take the layer to contain
!                boundary layer aerosol, provisionally setting
!                the index of the top of the boundary layer to this
!                level: it will be overwritten if higher layers are
!                also in the boundary layer. A upper limit is applied
!                to the number of the layers that can be filled with
!                the boundary layer aerosol so that there is at least
!                one layer in the free troposphere.
                 bl_top_lyr(l)=i
!                Increment the height of the bottom of the layer
!                for the next pass.
                 z_bottom(l)=z_bottom(l)+dz
               ELSE
!                The top of the boundary layer has been passed at this
!                point: it does not need to be considered further.
                 l_set_z(l)=.FALSE.
                 n_points_set_z=n_points_set_z-1
               END IF
             END IF
           END DO
           i=i-1
         END DO
       ELSE
         DO l=1, n_profile
!          We do not allow the stratosphere to meet the boundary
!          layer, so we ensure that, counting up from the surface,
!          the highest layer allowed to be filled with boundary-
!          layer aerosol is the TRINDX-2nd: as a result, there must
!          be at least one layer in the free troposphere.
           bl_top_lyr(l)=n_layer+1                                      &
             -MIN(n_levels_bl, trindx(i_gather(l))-2)
         END DO
       END IF

!      Now, the mixing ratios are set for the 5 aerosol species
!      in the free troposphere. Initially set the mixing ratio
!      in the lowest layer of the free troposphere.
       DO i=1,5
         DO l=1, n_profile
           lg=i_gather(l)
           aerosol_mix_ratio_clim(l, bl_top_lyr(l)-1, i)                &
             =freetrop_mass(i)*g*                                       &
             pressure_wt(l)/                                            &
             (p_layer_boundaries(lg, n_layer+1-bl_top_lyr(l))-          &
              p_layer_boundaries(lg,trindx(lg)-1) )
         END DO
       END DO
!      Fill in the remaining levels where necessary.
       DO l=1, n_profile
         lg=i_gather(l)
         DO i=1,5
           DO j=n_layer+2-trindx(lg), bl_top_lyr(l)-2
             aerosol_mix_ratio_clim(l, j, i)=                           &
               aerosol_mix_ratio_clim(l, bl_top_lyr(l)-1, i)
           END DO
         END DO
       END DO

!      Now, the boundary layer mixing ratios are set for the
!      5 aerosol species. A continental aerosol is used over most land
!      areas, but not over ice sheets, which are identified by the
!      criterion used in the boundary layer scheme that the mass of
!      lying snow exceeds 5000 kgm-2. Over ice sheets a maritime
!      aerosol is used.
       DO i=1,5
         DO l=1, n_profile
           lg=i_gather(l)
           IF ( land(lg).AND.(lying_snow(lg) <  5.0e+03) ) THEN
             aerosol_mix_ratio_clim(l, bl_top_lyr(l), i)                &
               =bl_landmass(i)*g*pressure_wt(l)                         &
                 /(pstar(lg)                                            &
                 -p_layer_boundaries(lg, n_layer+1-bl_top_lyr(l)))
           ELSE
             aerosol_mix_ratio_clim(l, bl_top_lyr(l), i)                &
               =bl_oceanmass(i)*g*pressure_wt(l)                        &
                 /(pstar(lg)                                            &
                 -p_layer_boundaries(lg, n_layer+1-bl_top_lyr(l)))
           END IF
         END DO
       END DO
       DO i=1,5
         DO l=1, n_profile
           DO j=bl_top_lyr(l)+1, n_layer
            aerosol_mix_ratio_clim(l,j,i)=                              &
               aerosol_mix_ratio_clim(l, bl_top_lyr(l), i)
           END DO
         END DO
       END DO



      IF (lhook) CALL dr_hook('R2_SET_AERO_CLIM_HADCM3',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE r2_set_aero_clim_hadcm3
