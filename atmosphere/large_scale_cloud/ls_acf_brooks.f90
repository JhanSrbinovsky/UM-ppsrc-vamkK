! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!     Subroutine Interface:

SUBROUTINE ls_acf_brooks (                                              &
! in coordinate information
  delta_lambda, delta_phi,                                              &
! trig arrays
  fv_cos_theta_latitude,                                                &
! in data fields
  bulk_cloud_fraction, cloud_fraction_liquid,                           &
  cloud_fraction_frozen,                                                &
! in logical control
  cumulus,                                                              &
! out data fields
  area_cloud_fraction)

!     Description:
!       Calculates area_cloud_fraction from bulk_cloud_fraction

!     Method:
!       The calculation is  based on the parametrisation in
!       Brooks 2005 equations 2-3.
!       (Brooks et al, July 2005, JAS vol 62 pp 2248-2260)
!       The initial parametrisation uses the values in equations
!       4,5,7 and 8 of the paper for ice and liquid cloud without
!       wind shear.  For mixed phase clouds the maximum of the two
!       area_cloud_fractions resulting will be used.
!       Only area_cloud_fraction will be updated.
!       Grid box size is needed to be known.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
!     Code Description:
!       FORTRAN 77 with extensions recommended in the Met. Office
!       F77 Standard.

  USE atm_fields_bounds_mod, ONLY: tdims_l, tdims_s, qdims
  USE level_heights_mod, ONLY:  r_theta_levels
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

! Array Arguments with INTENT(IN)
! Co-ordinate arrays:
  REAL ::                                                               &
   fv_cos_theta_latitude (tdims_s%i_start:tdims_s%i_end,                & 
                          tdims_s%j_start:tdims_s%j_end),               &
!       Finite volume cos(lat)
   delta_lambda,                                                        &
!       EW (x) grid spacing in radians
   delta_phi        
!       NS (y) grid spacing in radians

! Data arrays:
  REAL ::                                                               &
   bulk_cloud_fraction (  qdims%i_start:qdims%i_end,                    & 
                          qdims%j_start:qdims%j_end,                    &
                                      1:qdims%k_end),                   &
!       Cloud fraction at processed levels (decimal fraction).
   cloud_fraction_liquid (qdims%i_start:qdims%i_end,                    & 
                          qdims%j_start:qdims%j_end,                    &
                                      1:qdims%k_end),                   &
!       Liquid cloud fraction at processed levels (decimal fraction).
   cloud_fraction_frozen (qdims%i_start:qdims%i_end,                    & 
                          qdims%j_start:qdims%j_end,                    &
                                      1:qdims%k_end)
!       Frozen cloud fraction at processed levels (decimal fraction).

  LOGICAL, INTENT(IN)::                                                 &
   cumulus(               qdims%i_start:qdims%i_end,                    & 
                          qdims%j_start:qdims%j_end)

! Arguments with INTENT(OUT):
! Data arrays:
  REAL ::                                                               &
   area_cloud_fraction (  qdims%i_start:qdims%i_end,                    & 
                          qdims%j_start:qdims%j_end,                    &
                                      1:qdims%k_end)
!       Cloud fraction at processed levels (decimal fraction).

! Local Parameters:

! Parameters for liquid clouds, from Brooks 2005 equations 7 and 8
  REAL, PARAMETER :: power_law_gradient_liquid =  0.1635 ! A
  REAL, PARAMETER :: vert_fit_liquid           =  0.6694 ! alpha
  REAL, PARAMETER :: horiz_fit_liquid          = -0.1882 ! beta

! Parameters for frozen clouds, from Brooks 2005 equations 4 and 5
  REAL, PARAMETER :: power_law_gradient_frozen =  0.0880 ! A
  REAL, PARAMETER :: vert_fit_frozen           =  0.7679 ! alpha
  REAL, PARAMETER :: horiz_fit_frozen          = -0.2254 ! beta

! Local Scalars:
! Loop counters
  INTEGER ::                                                            &
   i, j, k

  REAL ::                                                               &
   symmetric_adjustment_liquid,                                         &
!    function f in eqn 7 in Brooks 2005
   symmetric_adjustment_frozen,                                         &
!    function f in eqn 4 in Brooks 2005
   horiz_scale,                                                         &
!    horizontal scale size of the grid box (m)
   vert_scale
!    vertical scale size of the grid box (m)

!  Local Arrays:
  REAL ::                                                               &
   acf_liquid (           qdims%i_start:qdims%i_end,                    & 
                          qdims%j_start:qdims%j_end,                    &
                                      1:qdims%k_end),                   &
!    area cloud fraction based on liquid parameters
   acf_frozen (           qdims%i_start:qdims%i_end,                    & 
                          qdims%j_start:qdims%j_end,                    &
                                      1:qdims%k_end)
!    area cloud fraction based on frozen parameters

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!-    End of header
! ----------------------------------------------------------------------

  IF (lhook) CALL dr_hook('LS_ACF_BROOKS',zhook_in,zhook_handle)

! ==Main Block==--------------------------------------------------------

! Initialise arrays and local variables to zero
  DO k = 1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        area_cloud_fraction(i,j,k) = 0.0
        acf_liquid(i,j,k) = 0.0
        acf_frozen(i,j,k) = 0.0
      END DO
    END DO
  END DO
  horiz_scale = 0.0
  vert_scale = 0.0
  symmetric_adjustment_liquid = 0.0
  symmetric_adjustment_frozen = 0.0

  DO k = 1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end

! Test if bulk_cloud_fraction is within bounds of possibility

        IF ( bulk_cloud_fraction(i,j,k) <= 0.0 ) THEN

          area_cloud_fraction(i,j,k) = 0.0

        ELSE IF ( bulk_cloud_fraction(i,j,k) >= 1.0 ) THEN

          area_cloud_fraction(i,j,k) = 1.0

        ELSE IF (cumulus(i,j)) THEN
              ! This is a convective point so do not apply the
              ! area cloud representation
          area_cloud_fraction(i,j,k) =                                  &
            bulk_cloud_fraction(i,j,k)

        ELSE

! Only calculate area_cloud_fraction if the bulk_cloud_fraction
! is between (not equal to) 0.0 and 1.0

! Calculate horizontal and vertical scales.
! The horizontal scale is taken as the square root of the
! area of the grid box.
! The vertical scale is taken as the difference in radius
! from the centre of the Earth between the upper and lower
! boundaries of the grid box.

          horiz_scale = SQRT (                                          &
                            r_theta_levels(i,j,k)                       &
                            * r_theta_levels(i,j,k)                     &
                            * delta_lambda * delta_phi                  &
                            * fv_cos_theta_latitude(i,j) )
          IF (k  ==  qdims%k_end) THEN
                ! Assume top layer thickness is the same as the
                ! thickness of the layer below
            vert_scale =  r_theta_levels(i,j,k)                         &
                         - r_theta_levels(i,j,k-1)
          ELSE
            vert_scale =  r_theta_levels(i,j,k+1)                       &
                         - r_theta_levels(i,j,k)
          END IF  ! k eq wet_model_levels

! Calculate the symmetric_adjustment (f).
! This parameter controls the extent to which the area cloud fraction
! is greater than the bulk cloud fraction.  If f = 0, they are equal.

          symmetric_adjustment_liquid =                                 &
                     power_law_gradient_liquid                          &
                     * ( vert_scale ** vert_fit_liquid )                &
                     * ( horiz_scale ** horiz_fit_liquid )
          symmetric_adjustment_frozen =                                 &
                     power_law_gradient_frozen                          &
                     * ( vert_scale ** vert_fit_frozen )                &
                     * ( horiz_scale ** horiz_fit_frozen )

! Calculate the area cloud fractions for liquid and frozen cloud
! Calculate the liquid and frozen fractions separately to
! allow for greatest flexibility in future choice of decisions
! regarding mixed phase cloud.

          acf_liquid(i,j,k) = 1./                                       &
              ( 1. + ( EXP(-1.*symmetric_adjustment_liquid)             &
                       * ( 1./bulk_cloud_fraction(i,j,k) - 1.) ) )
          acf_frozen(i,j,k) = 1./                                       &
              ( 1. + ( EXP(-1.*symmetric_adjustment_frozen)             &
                       * ( 1./bulk_cloud_fraction(i,j,k) - 1.) ) )

! Calculate the final area cloud fraction for each grid box
! Currently this is based on which there is more of, ice or liquid.

          IF ( cloud_fraction_frozen(i,j,k) == 0.0 ) THEN
            IF ( cloud_fraction_liquid(i,j,k) == 0.0 ) THEN

! If there is no liquid or frozen cloud, there should be no area cloud
              area_cloud_fraction(i,j,k) = 0.0

            ELSE

! If there is no frozen cloud but there is liquid cloud,
! then the area cloud fraction is given by the liquid parametrisation
! 0 no cloud, 1 either, 2 liq, 3 ice'
              area_cloud_fraction(i,j,k) = acf_liquid(i,j,k)
            END IF

          ELSE ! cloud_fraction_frozen

            IF ( cloud_fraction_liquid(i,j,k) == 0.0 ) THEN

! If there is frozen cloud but there is no liquid cloud,
! then the area cloud fraction is given by the frozen parametrisation
              area_cloud_fraction(i,j,k) = acf_frozen(i,j,k)

            ELSE

! If there is frozen cloud and there is liquid cloud,
! then the area cloud fraction is given by the maximum of the two
! parametrisations
              area_cloud_fraction(i,j,k) =                              &
                 MAX( acf_liquid(i,j,k),acf_frozen(i,j,k) )

            END IF

          END IF ! cloud_fraction_frozen

        END IF ! bulk_cloud_fraction between 0.0 and 1.0

      END DO
    END DO
  END DO

  IF (lhook) CALL dr_hook('LS_ACF_BROOKS',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE ls_acf_brooks
! ======================================================================
