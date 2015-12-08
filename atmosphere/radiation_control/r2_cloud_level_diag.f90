! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Subroutine to calculate an observed effective radius.

! Purpose:
!   An effective radius as observed from above the cloud-top is
!   calculated.

! Method:
!   For each type of cloud containing water in any layer the effective
!   radius is weighted with the product of the area of the cloud and the
!   probability that light emitted from the cloud reaches the observing
!   instrument.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.

!- ---------------------------------------------------------------------
      SUBROUTINE r2_cloud_level_diag(ierr, n_profile, n_layer, nclds    &
         , i_gather                                                     &
         , i_cloud, i_cloud_representation                              &
         , t, w_cloud, frac_cloud, l_all_temps                          &
         , condensed_mix_ratio, condensed_re                            &
         , l_observed_re, weighted_re, sum_weight_re                    &
         , col_list, row_list, row_length, rows                         &
         , npd_field, npd_profile, npd_layer                            &
         , first_layer, id_ct0, id_ct1                                  &
         )


      USE rad_pcf
      USE water_constants_mod, ONLY: tm
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE ereport_mod, ONLY : ereport
      USE dimfix3a_mod, ONLY: npd_cloud_component, npd_cloud_type,      &
                          npd_cloud_representation, npd_overlap_coeff,  &
                          npd_source_coeff, npd_region
      IMPLICIT NONE


!     COMDECKS INCLUDED.


      INTEGER                                                           &
                !, INTENT(OUT)
           ierr
!             ERROR FLAG

!     DIMENSIONS OF ARRAYS:
      INTEGER, INTENT(IN) :: row_length
!                              Number of grid-points in EW-direction
!                              in the local domain
      INTEGER, INTENT(IN) :: rows
!                              Number of grid-points in NS-direction
!                              in the local domain
      INTEGER                                                           &
                !, INTENT(IN)
           npd_field                                                    &
!             SIZE OF ARRAY OF ARRAYS PASSED FROM MAIN CODE
         , npd_profile                                                  &
!             SIZE OF ARRAY OF PROFILES
         , npd_layer                                                    &
!             MAXIMUM NUMBER OF LAYERS
         , first_layer                                                  &
!             First layer for some variables
!             0 for flux_calc (3A/C), 1 for radiance_calc (3Z)
         , id_ct0, id_ct1
!             Topmost declared cloudy layer (only used with 3Z)

!     ACTUAL SIZES USED:
      INTEGER                                                           &
                !, INTENT(IN)
           n_profile                                                    &
!             NUMBER OF PROFILES
         , n_layer                                                      &
!             Number of atmospheric layers in radiation
         , nclds                                                        &
!             NUMBER OF CLOUDY LEVELS
         , i_gather(npd_field)
!             LIST OF GATHERED POINTS
      INTEGER, INTENT(IN) :: col_list(npd_field)
!                              EW indices of gathered points in the
!                              2-D domain
      INTEGER, INTENT(IN) :: row_list(npd_field)
!                              NS indices of gathered points in the
!                              2-D domain

!     LOGICAL FLAGS FOR DIAGNOSTICS
      LOGICAL                                                           &
                !, INTENT(IN)
           l_observed_re                                                &
!             FLAG TO ENABLE DIAGNOSIS OF EFFECTIVE RADIUS SEEN FROM
!             SPACE (N.B. THE ROUTINE IS AT PRESENT CALLED ONLY IF
!             THIS IS TRUE, BUT ITS PRESENCE HERE ALLOWS FOR POSSIBLE
!             FUTURE EXTENSION OF THE ROUTINE).
         , l_all_temps
!             IF TRUE, THE ROUTINE HAS BEEN CALLED TO OBTAIN DIAGNOSTICS
!             FOR CLOUDS CONSISTING OF LIQUID WATER AT ANY TEMPERATURE
!             (AS DONE IN MODIS RETRIEVALS). IF FALSE, ONLY CLOUDS WITH
!             TEMPERATURES ABOVE FREEZING ARE TO BE DIAGNOSED (AS DONE
!             IN AVHRR RETRIEVALS).

!     REPRESENTATION OF CLOUDS
      INTEGER                                                           &
                !, INTENT(IN)
           i_cloud_representation                                       &
!             REPRESENTATION OF CLOUDS
         , i_cloud
!             TREATMENT OF OVERLAPS

      REAL                                                              &
                !, INTENT(IN)
           w_cloud(npd_profile, id_ct1:npd_layer)                       &
!             TOTAL AMOUNTS OF CLOUD
         , frac_cloud(npd_profile, id_ct1:npd_layer, npd_cloud_type)    &
!             FRACTION OF TYPES OF CLOUD
         , condensed_re(npd_profile, id_ct0:npd_layer                   &
            , npd_cloud_component)                                      &
!             EFFECTIVE RADII OF CLOUDY COMPONENTS
         , condensed_mix_ratio(npd_profile, id_ct0: npd_layer           &
            , npd_cloud_component)                                      &
!             MASS MIXING RATIOS OF CONDENSED COMPONENTS
         , t(npd_profile, first_layer:npd_layer)
!             TEMPERATURE

      REAL                                                              &
                !, INTENT(OUT)
           weighted_re(row_length, rows)                                &
!             WEIGHTED SUM OF EFFECTIVE RADIUS AND WEIGHTING FUNCTION
         , sum_weight_re(row_length, rows)
!             SUM OF WEIGHTS FOR EFFECTIVE RADIUS



!     LOCAL VARIABLES:
      INTEGER                                                           &
           i                                                            &
!             LOOP VARIABLE
         , l                                                            &
!             LOOP VARIABLE
         , i_inv
!             INVERTED LOOP INDEX
      REAL                                                              &
           trans_overlying_space(npd_profile)                           &
!             PROBABILITY OF A PHOTON IN CLEAR AIR IN THE LEVEL ABOVE
!             THE CURRENT ONE REACHING SPACE
         , area_exposed(npd_profile)                                    &
!             TOTAL AREA OF CLOUD IN THE CURRENT LAYER EXPOSED TO
!             CLEAR AIR IN THE LAYER ABOVE
         , area_exposed_st(npd_profile)                                 &
!             TOTAL AREA OF STRATIFORM CLOUD IN THE CURRENT LAYER
!             EXPOSED TO CLEAR AIR IN THE LAYER ABOVE
         , area_exposed_cnv(npd_profile)                                &
!             TOTAL AREA OF CONVECTIVE CLOUD IN THE CURRENT LAYER
!             EXPOSED TO CLEAR AIR IN THE LAYER ABOVE
         , area_clear_above(npd_profile)                                &
!             AREA OF THE CLEAR SKY REGION IN THE LAYER ABOVE
         , area_strat(npd_profile)                                      &
!             AREA OF STRATIFORM CLOUD IN THE CURRENT LAYER
         , area_strat_above(npd_profile)                                &
!             AREA OF STRATIFORM CLOUD IN THE LAYER ABOVE
         , area_conv(npd_profile)                                       &
!             AREA OF CONVECTIVE CLOUD IN THE CURRENT LAYER
         , area_conv_above(npd_profile)                                 &
!             AREA OF CONVECTIVE CLOUD IN THE LAYER ABOVE
         , area_clear_clear(npd_profile)                                &
!             AREA OF BOUNDARY WHERE CLEAR SKY OVERLIES CLEAR SKY
         , area_clear(npd_profile)                                      &
!             AREA OF CLEAR SKY IN THE CURRENT LAYER
!             DOWN TO A LEVEL
         , area_uncorrelated(npd_profile)                               &
!             UNCORRELATED REGION ON THE INTERFACE
         , weighted_re_g(npd_profile)                                   &
!             WEIGHTED SUM OF EFFECTIVE RADIUS AND WEIGHTING FUNCTION
         , sum_weight_re_g(npd_profile)
!             SUM OF WEIGHTS FOR EFFECTIVE RADIUS

!     VARIABLES FOR GATHERING
      INTEGER                                                           &
           n_list                                                       &
!             NUMBER OF POINTS IN LIST
         , l_list(npd_profile)
!             INDICES OF POINTS IN LIST

!     INDICATOR FUNCTION
      REAL                                                              &
           chi_cnv(npd_profile)                                         &
!             CONVECTIVE INDICATOR FUNCTION
         , chi_st(npd_profile)
!             STRATIFORM INDICATOR FUNCTION

      CHARACTER (LEN=*), PARAMETER :: RoutineName = 'r2_cloud_level_diag'
      CHARACTER (LEN=240) :: cmessage

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle




!     INITIALIZATION OF LOCAL FIELDS.

      IF (lhook) CALL dr_hook('R2_CLOUD_LEVEL_DIAG',zhook_in,zhook_handle)

      IF (l_observed_re) THEN
         DO l=1, n_profile
            weighted_re_g(l)=0.0e+00
            sum_weight_re_g(l)=0.0e+00
         END DO
      END IF

!     INITIALIZE THE TRANSMISION ABOVE CLOUDS.
      DO l=1, n_profile
         trans_overlying_space(l)=1.0e+00
         area_clear_above(l)=1.0e+00
      END DO
      IF (l_observed_re.AND.(i_cloud == ip_cloud_triple)) THEN
         DO l=1, n_profile
            area_strat_above(l)=0.0e+00
            area_conv_above(l)=0.0e+00
         END DO
      END IF

!     STEP DOWN THROUGH THE ATMOSPHERE CALCULATING CONTRIBUTIONS TO
!     THE DIAGNOSTICS AND SUBSEQUENTLY ALLOWING FOR TRANSMISSION
!     THROUGH THE CURRENT LAYER.

      DO i=nclds, 1, -1
         i_inv=n_layer+1-i

         DO l=1, n_profile
            area_clear(l)=1.0e+00-w_cloud(l, i_inv)
         END DO

!        CALCULATE THE LOCAL AREA OF CLOUD RADIATING INTO CLEAR AIR.
         IF (i_cloud == ip_cloud_mix_random) THEN
            DO l=1, n_profile
               area_exposed(l)=w_cloud(l, i_inv)                        &
                  *area_clear_above(l)
            END DO
         ELSE IF ( (i_cloud == ip_cloud_mix_max).OR.                    &
                   (i_cloud == ip_cloud_triple) ) THEN
            DO l=1, n_profile
               area_exposed(l)=MAX(0.0e+00, (w_cloud(l, i_inv)          &
                  +area_clear_above(l)-1.0e+00))
            END DO
         END IF




         IF (l_observed_re) THEN



            IF (i_cloud_representation == ip_cloud_conv_strat) THEN


               IF ( (i_cloud == ip_cloud_mix_max).OR.                   &
                    (i_cloud == ip_cloud_mix_random) ) THEN

!                 IF THE OVERLAP OF CONVECTIVE CLOUD IS NOT ASSUMED
!                 TO BE COHERENT THE OVERALL EXPOSED AREA MAY BE
!                 PARTITIONED ACCORDING TO THE FRACTIONAL
!                 CONTRIBUTIONS OF CLOUD IN THE CURRENT LAYER.

                  DO l=1, n_profile
                     area_exposed_st(l)=area_exposed(l)                 &
                        *frac_cloud(l, i_inv, ip_cloud_type_strat)
                     area_exposed_cnv(l)=area_exposed(l)                &
                        *frac_cloud(l, i_inv, ip_cloud_type_conv)
                  END DO

               ELSE IF (i_cloud == ip_cloud_triple) THEN

!                 HERE, THE DIFFERENT TYPES OF CLOUDS OVERLAP
!                 COHERENTLY SO STRATIFORM CLOUD WILL BE EXPOSED
!                 ONLY IF THERE IS LESS STRATIFORM CLOUD IN THE
!                 LAYER ABOVE AND MORE CLEAR AIR IN THE LAYER ABOVE:
!                 UNDER THESE CONDITIONS THE NON-CORRELATED AREAS
!                 OVERLAP RANDOMLY.

                  DO l=1, n_profile
                     area_strat(l)=w_cloud(l, i_inv)                    &
                        *frac_cloud(l, i_inv, ip_cloud_type_strat)
                     area_conv(l)=w_cloud(l, i_inv)                     &
                       *frac_cloud(l, i_inv, ip_cloud_type_conv)
                     area_uncorrelated(l)=1.0e+00                       &
                       -MIN(area_clear(l), area_clear_above(l))         &
                       -MIN(area_strat(l), area_strat_above(l))         &
                       -MIN(area_conv(l), area_conv_above(l))
!                    First find the area of uncorrelated
!                    stratiform cloud.
                     area_exposed_st(l)=MAX(0.0e+00                     &
                        , (area_strat(l)-area_strat_above(l)))
                     area_exposed_st(l)=MAX(0.0e+00, area_exposed_st(l) &
                        *(area_clear_above(l)-area_clear(l)))
!                    Now normalize within the uncorrelated region.
!                    If the uncorrelated area is 0 the exposed area
!                    must be 0, so no second branch of the IF-test
!                    is required.
                     IF (area_uncorrelated(l) >  0.0e+00)               &
                       area_exposed_st(l)                               &
                         =area_exposed_st(l)/area_uncorrelated(l)
                     area_exposed_cnv(l)                                &
                        =area_exposed(l)-area_exposed_st(l)
                  END DO
               ELSE
                  cmessage =                                            &
                     '*** ERROR: THE DIAGNOSTIC OF OBSERVED RE HAS NOT '&
                     //'BEEN IMPLEMENTED WITH THIS OVERLAP OPTION.'
                  ierr=i_err_fatal
                  GO TO 9999
               END IF

!              THE INDICATOR FUNCTIONS FOR LIQUID WATER IN
!              CONVECTIVE OR STRAIFORM CLOUDS ARE SET TO 1
!              IF THERE IS ANY LIQUID WATER AND TO 0 OTHERWISE.
               DO l=1, n_profile
                  IF (condensed_mix_ratio(l, i_inv, ip_clcmp_cnv_water) &
                      >  0.0e+00) THEN
                     chi_cnv(l)=1.0e+00
                  ELSE
                     chi_cnv(l)=0.0e+00
                  END IF
                  IF (condensed_mix_ratio(l, i_inv, ip_clcmp_st_water)  &
                      >  0.0e+00) THEN
                     chi_st(l)=1.0e+00
                  ELSE
                     chi_st(l)=0.0e+00
                  END IF
               END DO

!              INCLUDE CONTRIBUTIONS FROM CONVECTIVE AND STRATIFORM
!              WATER CLOUDS.
               DO l=1, n_profile
                  weighted_re_g(l)=weighted_re_g(l)                     &
                     +trans_overlying_space(l)                          &
                     *(area_exposed_cnv(l)*chi_cnv(l)                   &
                     *condensed_re(l, i_inv, ip_clcmp_cnv_water)        &
                     +area_exposed_st(l)*chi_st(l)                      &
                     *condensed_re(l, i_inv, ip_clcmp_st_water))
                  sum_weight_re_g(l)=sum_weight_re_g(l)                 &
                     +trans_overlying_space(l)                          &
                     *(area_exposed_cnv(l)*chi_cnv(l)                   &
                     +area_exposed_st(l)*chi_st(l))
               END DO

            ELSE IF ((i_cloud_representation == ip_cloud_csiw).OR.      &
                    (i_cloud_representation == ip_cloud_ice_water)) THEN

               IF ( (i_cloud == ip_cloud_mix_max).OR.                   &
                    (i_cloud == ip_cloud_mix_random) ) THEN

!                 IF THE OVERLAP OF CONVECTIVE CLOUD IS NOT ASSUMED
!                 TO BE COHERENT THE OVERALL EXPOSED AREA MAY BE
!                 PARTITIONED ACCORDING TO THE FRACTIONAL
!                 CONTRIBUTIONS OF CLOUD IN THE CURRENT LAYER.
!                 THE EXPOSED AREAS INCLUDE ONLY THE PARTS OF THE
!                 CLOUDS CONTAINING WATER DROPLETS.

                  DO l=1, n_profile
                     area_exposed_st(l)=area_exposed(l)                 &
                        *frac_cloud(l, i_inv, ip_cloud_type_sw)
                     area_exposed_cnv(l)=area_exposed(l)                &
                        *frac_cloud(l, i_inv, ip_cloud_type_cw)
                  END DO

               ELSE IF (i_cloud == ip_cloud_triple) THEN

!                 HERE, THE DIFFERENT TYPES OF CLOUDS OVERLAP
!                 COHERENTLY SO STRATIFORM CLOUD WILL BE EXPOSED
!                 ONLY IF THERE IS LESS STRATIFORM CLOUD IN THE
!                 LAYER ABOVE AND MORE CLEAR AIR IN THE LAYER ABOVE:
!                 UNDER THESE CONDITIONS THE NON-CORRELATED AREAS
!                 OVERLAP RANDOMLY.
!                 THE ACTUAL EXPOSED AREAS OF CONVECTIVE OR
!                 STRATIFORM CLOUD MUST THEN BE WEIGHTED BY FACTORS
!                 REPRESENTING THE LIQUID PORTION OF EACH CLOUD, SINCE
!                 NOTHING IS RETRIEVED OVER ICE. (THE HORIZONTAL
!                 ARRANGEMENT OF ICE AND WATER WITHIN EITHER TYPE OF
!                 CLOUD IS RANDOM).

                  DO l=1, n_profile

                     area_strat(l)=w_cloud(l, i_inv)                    &
                        *(frac_cloud(l, i_inv, ip_cloud_type_sw)        &
                        +frac_cloud(l, i_inv, ip_cloud_type_si))
                     area_conv(l)=w_cloud(l, i_inv)                     &
                        *(frac_cloud(l, i_inv, ip_cloud_type_cw)        &
                        +frac_cloud(l, i_inv, ip_cloud_type_ci))
                     area_uncorrelated(l)=1.0e+00                       &
                        -MIN(area_clear(l), area_clear_above(l))        &
                        -MIN(area_strat(l), area_strat_above(l))        &
                        -MIN(area_conv(l), area_conv_above(l))
                     area_exposed_st(l)=MAX(0.0e+00                     &
                        , (area_strat(l)-area_strat_above(l)))
                     IF (area_uncorrelated(l) >  0.0e+00) THEN
                        area_exposed_st(l)                              &
                           =MAX(0.0e+00, area_exposed_st(l)             &
                           *(area_clear_above(l)-area_clear(l)))        &
                           /area_uncorrelated(l)
                     ELSE
                        area_exposed_st(l)=0.0e+00
                     END IF
                     area_exposed_cnv(l)                                &
                        =area_exposed(l)-area_exposed_st(l)

                     IF (frac_cloud(l, i_inv, ip_cloud_type_cw)         &
                         >  0.0e+00) THEN
                        area_exposed_cnv(l)=area_exposed_cnv(l)         &
                           /(1.0e+00                                    &
                           +frac_cloud(l, i_inv, ip_cloud_type_ci)      &
                           /frac_cloud(l, i_inv, ip_cloud_type_cw))
                     ELSE
                        area_exposed_cnv(l)=0.0e+00
                     END IF

                     IF (frac_cloud(l, i_inv, ip_cloud_type_sw)         &
                         >  0.0e+00) THEN
                        area_exposed_st(l)=area_exposed_st(l)           &
                           /(1.0e+00                                    &
                           +frac_cloud(l, i_inv, ip_cloud_type_si)      &
                           /frac_cloud(l, i_inv, ip_cloud_type_sw))
                     ELSE
                        area_exposed_st(l)=0.0e+00
                     END IF

                  END DO
               ELSE
                  cmessage =                                            &
                     '*** ERROR: THE DIAGNOSTIC OF OBSERVED RE HAS NOT '&
                     //'BEEN IMPLEMENTED WITH THIS OVERLAP OPTION.'
                  ierr=i_err_fatal
                  GO TO 9999
               END IF


               DO l=1, n_profile

                IF ((t(l, i_inv)  >   tm) .OR. l_all_temps) THEN
                  weighted_re_g(l)=weighted_re_g(l)                     &
                     +trans_overlying_space(l)                          &
                     *(area_exposed_cnv(l)                              &
                     *condensed_re(l, i_inv, ip_clcmp_cnv_water)        &
                     +area_exposed_st(l)                                &
                     *condensed_re(l, i_inv, ip_clcmp_st_water))
                  sum_weight_re_g(l)=sum_weight_re_g(l)                 &
                     +trans_overlying_space(l)                          &
                     *(area_exposed_cnv(l)+area_exposed_st(l))
                END IF
               END DO

            END IF


         END IF



!        ADVANCE THE STORED QUANTITIES REFFERRING TO OVERLYING LAYERS.


!        THE TRANSMISSION TO SPACE CURRENTLY HOLDS THE PROBABILITY THAT
!        A PHOTON TRAVELLING UPWARDS IN THE CLEAR AIR IN THE LAYER ABOVE
!        WILL ESCAPE TO SPACE WITHOUT ENCOUNTERING A CLOUD. TO ADVANCE
!        THIS TO THE CURRENT LAYER IT MUST BE MULTIPLIED BY A FACTOR
!        REPRESENTING THE OVERLAP ASSUMPTION AT THE TOP OF THE PRESENT
!        LAYER.

         IF (i_cloud == ip_cloud_mix_random) THEN

            DO l=1, n_profile
               trans_overlying_space(l)=trans_overlying_space(l)        &
                  *area_clear_above(l)
            END DO

         ELSE IF ( (i_cloud == ip_cloud_mix_max).OR.                    &
                   (i_cloud == ip_cloud_triple) ) THEN

            DO l=1, n_profile
               area_clear_clear(l)=MIN(area_clear(l)                    &
                  , area_clear_above(l))
               IF (area_clear(l) >  0.0e+00) THEN
                  trans_overlying_space(l)=trans_overlying_space(l)     &
                     *area_clear_clear(l)/area_clear(l)
               ELSE
                  trans_overlying_space(l)=0.0e+00
               END IF
            END DO

         END IF

!        ADVANCE THE AREAS OF CLOUD.
         DO l=1, n_profile
            area_clear_above(l)=area_clear(l)
         END DO
         IF (i_cloud_representation == ip_cloud_conv_strat) THEN
            DO l=1, n_profile
               area_strat_above(l)=w_cloud(l, i_inv)                    &
                  *frac_cloud(l, i_inv, ip_cloud_type_strat)
            END DO
         ELSE IF ((i_cloud_representation == ip_cloud_csiw).OR.          &
                  (i_cloud_representation == ip_cloud_ice_water)) THEN
            DO l=1, n_profile
               area_strat_above(l)=area_strat(l)
               area_conv_above(l)=area_conv(l)
            END DO
         END IF

      END DO



      IF (l_observed_re) THEN
!        SCATTER THE DIAGNOSTICS BACK TO THE OUTPUT ARRAYS AND CONVERT
!        TO MICRONS (TO AVOID FIELDS BEING CORRUPTED BY PACKING).
         DO l=1, n_profile
            weighted_re(col_list(l), row_list(l))                       &
              =1.0e+06*weighted_re_g(l)
            sum_weight_re(col_list(l), row_list(l))                     &
              =sum_weight_re_g(l)
         END DO
      END IF


 9999 CONTINUE
! Check error condition
      IF (ierr /= i_normal) THEN

        CALL ereport(RoutineName, ierr, cmessage)
      END IF

      IF (lhook) CALL dr_hook('R2_CLOUD_LEVEL_DIAG',zhook_out,zhook_handle)
      END SUBROUTINE r2_cloud_level_diag
