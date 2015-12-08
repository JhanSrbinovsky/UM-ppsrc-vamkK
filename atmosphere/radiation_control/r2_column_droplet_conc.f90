! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Subroutine to calculate column-integrated cloud droplet number.

! Purpose:
!   To calculate a diagnostic of column-integrated cloud droplet
!   number which may be validated aginst satellite data.

! Method:
!   Column cloud droplet concentration (i.e. number of droplets per
!   unit area) is calculated as the vertically integrated droplet
!   number concentration averaged over the portion of the gridbox
!   covered by stratiform and convective liquid cloud with T>273K.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Radiation Control

! Description of Code:
!   FORTRAN 77  with extensions listed in documentation.

!- ---------------------------------------------------------------------
      SUBROUTINE r2_column_droplet_conc(npd_profile, npd_layer          &
         , n_profile, n_layer, nclds                                    &
         , strat_liq_cloud_fraction                                     &
         , total_strat_liq_cloud_fraction                               &
         , conv_liq_cloud_fraction                                      &
         , total_conv_liq_cloud_fraction                                &
         , n_drop, d_mass, density_air                                  &
         , nc_diag, nc_weight)




      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE



!  Input variables:

      INTEGER                                                           &
           npd_profile                                                  &
!             Maximum number of profiles
         , npd_layer                                                    &
!             Maximum number of layers
         , n_profile                                                    &
!             Number of atmospheric profiles
         , n_layer                                                      &
!             Number of layers seen in radiation
         , nclds
!             Number of cloudy layers

      REAL                                                              &
           strat_liq_cloud_fraction(npd_profile, npd_layer)             &
!             Stratiform liquid (T>273K) cloud cover in layers
         , total_strat_liq_cloud_fraction(npd_profile)                  &
!             Total liquid (T>273K) stratiform cloud cover
         , conv_liq_cloud_fraction(npd_profile, npd_layer)              &
!             Convective liquid (T>273K) cloud cover in layers
         , total_conv_liq_cloud_fraction(npd_profile)                   &
!             Total liquid (T>273K) convective cloud cover
         , n_drop(npd_profile, npd_layer)                               &
!             Number concentration of cloud droplets (m-3)
         , d_mass(npd_profile, npd_layer)                               &
!             Mass thickness of layer (kg m-2)
         , density_air(npd_profile, npd_layer)
!             Air density (kg m-3)


!  Output variables:

      REAL                                                              &
           nc_diag(npd_profile)                                         &
!             Column-integrated droplet number diagnostic (m-2)
         , nc_weight(npd_profile)
!             Weighting factor for column droplet number


!  Local variables:

      INTEGER                                                           &
           i, l
!             Loop counters

      REAL                                                              &
           n_sum_s                                                      &
         , n_sum_c                                                      &
!             Temporary sums
         , ncol_s                                                       &
         , wgt_s                                                        &
         , ncol_c                                                       &
         , wgt_c
!             N-column values and weights for stratiform
!             and convective clouds separately.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



      IF (lhook) CALL dr_hook('R2_COLUMN_DROPLET_CONC',zhook_in,zhook_handle)
      DO l=1, n_profile

         n_sum_s=0.0
         DO i=n_layer+1-nclds, n_layer
            n_sum_s=n_sum_s+(strat_liq_cloud_fraction(l, i)             &
                      *n_drop(l, i)*d_mass(l, i)/density_air(l, i))
         END DO

         n_sum_c=0.0
         DO i=n_layer+1-nclds, n_layer
            n_sum_c=n_sum_c+(conv_liq_cloud_fraction(l, i)              &
                      *n_drop(l, i)*d_mass(l, i)/density_air(l, i))
         END DO

         IF (total_strat_liq_cloud_fraction(l)  >   0.0) THEN
            ncol_s=n_sum_s/total_strat_liq_cloud_fraction(l)
            wgt_s=total_strat_liq_cloud_fraction(l)
         ELSE
            ncol_s=0.0
            wgt_s=0.0
         END IF

         IF (total_conv_liq_cloud_fraction(l)  >   0.0) THEN
            ncol_c=n_sum_c/total_conv_liq_cloud_fraction(l)
            wgt_c=total_conv_liq_cloud_fraction(l)
         ELSE
            ncol_c=0.0
            wgt_c=0.0
         END IF

         IF ((wgt_s+wgt_c)  >   0.0) THEN
            nc_diag(l)=((ncol_s*wgt_s)+(ncol_c*wgt_c))/(wgt_s+wgt_c)
            nc_weight(l)=wgt_s+wgt_c
         ELSE
            nc_diag(l)=0.0
            nc_weight(l)=0.0
         END IF

      END DO

!     Note: weighting is done later in R2_SET_CLOUD_FIELD



      IF (lhook) CALL dr_hook('R2_COLUMN_DROPLET_CONC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE r2_column_droplet_conc
