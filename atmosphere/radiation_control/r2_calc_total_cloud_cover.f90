! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to calculate the total cloud cover.

! Purpose:
!   The total cloud cover at all grid-points is determined.

! Method:
!   A separate calculation is made for each different assumption about
!   the overlap.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

!- ---------------------------------------------------------------------
      SUBROUTINE r2_calc_total_cloud_cover(n_profile, n_layer, nclds    &
         , i_cloud, w_cloud_in, total_cloud_cover                       &
         , npd_profile, npd_layer                                       &
         )



      USE rad_pcf
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE


!     DECLARATION OF ARRAY SIZES.
      INTEGER                                                           &
                !, INTENT(IN)
           npd_profile                                                  &
!             MAXIMUM NUMBER OF PROFILES
         , npd_layer
!             MAXIMUM NUMBER OF LAYERS


!     DUMMY ARGUMENTS.
      INTEGER                                                           &
                !, INTENT(IN)
           n_profile                                                    &
!             NUMBER OF PROFILES
         , n_layer                                                      &
!             Number of layers seen in radiation
         , nclds                                                        &
!             NUMBER OF CLOUDY LAYERS
         , i_cloud
!             CLOUD SCHEME EMPLOYED
      REAL                                                              &
                !, INTENT(IN)
           w_cloud_in(npd_profile, npd_layer)
!             CLOUD AMOUNTS

      REAL                                                              &
                !, INTENT(OUT)
           total_cloud_cover(npd_profile)
!             TOTAL CLOUD COVER


!     LOCAL VARIABLES.
      INTEGER                                                           &
           l                                                            &
!             LOOP VARIABLE
         , i
!             LOOP VARIABLE
      REAL                                                              &
           w_cloud(npd_profile, npd_layer)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



!     COPY W_CLOUD_IN TO W_CLOUD AND THEN CHECK THAT VALUES ARE
!     BETWEEN 0 AND 1.

      IF (lhook) CALL dr_hook('R2_CALC_TOTAL_CLOUD_COVER',zhook_in,zhook_handle)
         DO i=n_layer+1-nclds, n_layer
            DO l=1, n_profile
              w_cloud(l,i) = MIN( MAX(w_cloud_in(l,i),0.0) , 1.0 )
            END DO
         END DO

!     DIFFERENT OVERLAP ASSUMPTIONS ARE CODED INTO EACH SOLVER.

      IF (i_cloud == ip_cloud_mix_max) THEN

!        USE THE TOTAL CLOUD COVER TEMPORARILY TO HOLD THE CLEAR-SKY
!        FRACTION AND CONVERT BACK TO CLOUD COVER LATER.
!        WE CALCULATE THIS QUANTITY BY IMAGINING A TOTALLY TRANSPARENT
!        ATMOSPHERE CONTAINING TOTALLY OPAQUE CLOUDS AND FINDING THE
!        TRANSMISSION.
         DO l=1, n_profile
            total_cloud_cover(l)=1.0e+00-w_cloud(l, n_layer+1-nclds)
         END DO
         DO i=n_layer+1-nclds, n_layer-1
            DO l=1, n_profile
               IF (w_cloud(l, i+1) >  w_cloud(l, i)) THEN
                  total_cloud_cover(l)=total_cloud_cover(l)             &
                     *(1.0e+00-w_cloud(l, i+1))/(1.0e+00-w_cloud(l, i))
               END IF
            END DO
         END DO
         DO l=1, n_profile
            total_cloud_cover(l)=1.0e+00-total_cloud_cover(l)
         END DO

      ELSE IF (i_cloud == ip_cloud_mix_random) THEN

!        USE THE TOTAL CLOUD COVER TEMPORARILY TO HOLD THE CLEAR-SKY
!        FRACTION AND CONVERT BACK TO CLOUD COVER LATER.
         DO l=1, n_profile
            total_cloud_cover(l)=1.0e+00
         END DO
         DO i=n_layer+1-nclds, n_layer
            DO l=1, n_profile
               total_cloud_cover(l)=total_cloud_cover(l)                &
                  *(1.0e+00-w_cloud(l, i))
            END DO
         END DO
         DO l=1, n_profile
            total_cloud_cover(l)=1.0e+00-total_cloud_cover(l)
         END DO

      ELSE IF (i_cloud == ip_cloud_column_max) THEN

         DO l=1, n_profile
            total_cloud_cover(l)=0.0e+00
         END DO
         DO i=n_layer+1-nclds, n_layer
            DO l=1, n_profile
               total_cloud_cover(l)=MAX(total_cloud_cover(l)            &
                  , w_cloud(l, i))
            END DO
         END DO

      ELSE IF (i_cloud == ip_cloud_triple) THEN

!        USE THE TOTAL CLOUD COVER TEMPORARILY TO HOLD THE CLEAR-SKY
!        FRACTION AND CONVERT BACK TO CLOUD COVER LATER.
!        WE CALCULATE THIS QUANTITY BY IMAGINING A TOTALLY TRANSPARENT
!        ATMOSPHERE CONTAINING TOTALLY OPAQUE CLOUDS AND FINDING THE
!        TRANSMISSION.
         DO l=1, n_profile
            total_cloud_cover(l)=1.0e+00-w_cloud(l, n_layer+1-nclds)
         END DO
         DO i=n_layer+1-nclds, n_layer-1
            DO l=1, n_profile
               IF (w_cloud(l, i+1) >  w_cloud(l, i)) THEN
                  total_cloud_cover(l)=total_cloud_cover(l)             &
                     *(1.0e+00-w_cloud(l, i+1))/(1.0e+00-w_cloud(l, i))
               END IF
            END DO
         END DO
         DO l=1, n_profile
            total_cloud_cover(l)=1.0e+00-total_cloud_cover(l)
         END DO

      ELSE IF (i_cloud == ip_cloud_clear) THEN

         DO l=1, n_profile
            total_cloud_cover(l)=0.0e+00
         END DO

      END IF



      IF (lhook) CALL dr_hook('R2_CALC_TOTAL_CLOUD_COVER',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE r2_calc_total_cloud_cover
