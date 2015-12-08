! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ PC2 Cloud Scheme: Reset clouds for extreme relative total humidity.
! Subroutine Interface:
SUBROUTINE pc2_checks2(                                                 &
!      Pressure related fields
 p_theta_levels, rhcrit,                                                &
!      Array dimensions
 rhc_row_length,rhc_rows,                                               &
!      Prognostic Fields
 t, cf, cfl, cff, q, qcl,                                               &
!      Logical control
 l_mixing_ratio)

  USE water_constants_mod,   ONLY: lc
  USE atmos_constants_mod,   ONLY: cp, r, repsilon
  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim
  USE atm_fields_bounds_mod, ONLY: pdims, qdims, tdims
  USE pc2_constants_mod,     ONLY: cloud_pc2_tol, cloud_pc2_tol_2

  IMPLICIT NONE

! Purpose:
!   Check that cloud fraction is either zero or one when relative
!   total humidity is small or large.

! Method:
!   Calculate relative total humidity, compare to RHcrit and adjust
!   the cloud variables appropriately.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation

!  Subroutine Arguments:------------------------------------------------
  INTEGER ::                                                            &
                        !, INTENT(IN)
   rhc_row_length,rhc_rows
!       Dimensions of the RHCRIT variable.

  REAL ::                                                               &
                        !, INTENT(IN)
   p_theta_levels(pdims%i_start:pdims%i_end,                            & 
                  pdims%j_start:pdims%j_end,                            &  
                  pdims%k_start:pdims%k_end),                           &
!       pressure at all points (Pa)
   rhcrit(        rhc_row_length,                                       &
                  rhc_rows,                                             &
                              1:qdims%k_end),                           &
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for the give
!       set of levels.
   cff(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                              1:qdims%k_end)
!       Ice cloud fraction (no units)

  LOGICAL ::                                                            &
                        !, INTENT(IN)
   l_mixing_ratio   ! Use mixing ratio formulation

  REAL ::                                                               &
                        !, INTENT(INOUT)
   t(             tdims%i_start:tdims%i_end,                            & 
                  tdims%j_start:tdims%j_end,                            &  
                              1:tdims%k_end),                           &
!       Temperature (K)
   cf(            qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                              1:qdims%k_end),                           &
!       Total cloud fraction (no units)
   cfl(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                              1:qdims%k_end),                           &
!       Liquid cloud fraction (no units)
   q(             qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                              1:qdims%k_end),                           &
!       Vapour content (kg water per kg air)
   qcl(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                              1:qdims%k_end)
!       Liquid content (kg water per kg air)

!  External functions:

!  Local parameters and other physical constants------------------------
  REAL, PARAMETER :: lcrcp=lc/cp
!       Latent heat of condensation divided by heat capacity of air.
  REAL, PARAMETER :: c_thresh_low = cloud_pc2_tol
  REAL, PARAMETER :: c_thresh_low_2 = cloud_pc2_tol_2
!       Low cloud fraction thresholds
  REAL, PARAMETER :: c_thresh_high  = 1.0 - cloud_pc2_tol
  REAL, PARAMETER :: c_thresh_high_2 = 1.0 - cloud_pc2_tol_2
!       High cloud fraction thresholds

!  Local scalars--------------------------------------------------------

!  (a)  Scalars effectively expanded to workspace by the Cray (using
!       vector registers).
  REAL ::                                                               &
   alpha,                                                               &
!       Rate of change of saturation specific humidity with
!       temperature calculated at dry-bulb temperature (kg kg-1 K-1)
   al,                                                                  &
!       1 / (1 + alpha L/cp)  (no units)
   rht,                                                                 &
!       Relative total humidity
   sd         
!       Saturation deficit

!  (b)  Others.
  INTEGER :: k,i,j,                                                     &
!       Loop counters: K - vertical level index
!       I,J - horizontal position index
   irhi,irhj,                                                           &
!       Indices for RHcrit array
   multrhc,                                                             &
!       Zero if (rhc_row_length*rhc_rows) le 1, else 1
   npt

!  Local dynamic arrays-------------------------------------------------
!    3 blocks of real workspace are required.
  REAL ::                                                               &
   qsl_t(     (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
!       Saturated specific humidity for temperature T
   qsl_tl(    (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
!       Saturated specific humidity for liquid temperature TL
   tl_c(      (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
!    Liquid temperature (= T - L/cp QCL)  (K)
   t_c(       (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
   p_c(       (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) )
!    Temperature and pressure on compressed points for qsat call
  INTEGER ::                                                            &
   index_npt(qdims%i_start:qdims%i_end,                                 & 
             qdims%j_start:qdims%j_end)

!- End of Header

! Set up a flag to state whether RHcrit is a single parameter or defined
! on all points.

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL   (KIND=jprb)            :: zhook_handle

  IF (lhook) CALL dr_hook('PC2_CHECKS2',zhook_in,zhook_handle)

  IF (rhc_row_length * rhc_rows > 1) THEN
    multrhc=1
  ELSE
    multrhc=0
  END IF

! ==Main Block==--------------------------------------------------------

! Loop round levels to be processed
! Levels_do1:

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k,    & 
!$OMP& npt, irhi, irhj, rht, alpha, al, sd, index_npt, qsl_t, qsl_tl,   &
!$OMP& t_c, p_c, tl_c)
  DO k = 1, qdims%k_end

! Copy points into compressed arrays
    npt = 0
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        IF ( cfl(i,j,k)  >=  c_thresh_high .OR.                         &
             ( cfl(i,j,k)  <=  c_thresh_low .AND.                       &
               cfl(i,j,k) > 0. ) ) THEN
          npt = npt + 1
          index_npt(i,j) = npt
          t_c(npt) = t(i,j,k)
          tl_c(npt) = t(i,j,k)-lcrcp*qcl(i,j,k)
          p_c(npt) = p_theta_levels(i,j,k)
        END IF
      END DO
    END DO

! ----------------------------------------------------------------------
! 2. Calculate Saturated Specific Humidity with respect to liquid water
!    for dry bulb temperature and liquid temperature.
! ----------------------------------------------------------------------
    IF (npt > 0) THEN
! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(qsl_t, t_c, p_c, npt, l_mixing_ratio)
! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(qsl_tl, tl_c, p_c, npt, l_mixing_ratio)
    END IF

    DO j = qdims%j_start, qdims%j_end

      DO i = qdims%i_start, qdims%i_end

        IF ( cfl(i,j,k)  >=  c_thresh_high .OR.                         &
             ( cfl(i,j,k)  <=  c_thresh_low .AND.                       &
               cfl(i,j,k) > 0. ) ) THEN

! Set up index pointers to critical relative humidity value

          irhi = (multrhc * (i - 1)) + 1
          irhj = (multrhc * (j - 1)) + 1

! Calculate relative total humidity with respect to the liquid
! temperature and threshold relative humidity
          rht=(q(i,j,k)+qcl(i,j,k))/qsl_tl(index_npt(i,j))

! ----------------------------------------------------------------------
! 3. Determine whether resetting is required and calculate appropriate
!    increments if this is the case.
! ----------------------------------------------------------------------

! Is the relative total humidity greater than the threshold amount?
! If so, evaporate some liquid to set the saturation deficit to zero.

          IF ( (rht  >   (2.0-rhcrit(irhi,irhj,k)) .AND.                &
            (cfl(i,j,k) >= c_thresh_high)) .OR.                         &
            (cfl(i,j,k)  >=  c_thresh_high_2) ) THEN
            cfl(i,j,k) = 1.0
            cf(i,j,k)  = 1.0

! Calculate the saturation deficit

            alpha = repsilon * lc * qsl_t(index_npt(i,j)) /             &
                  (r * t_c(index_npt(i,j)) ** 2)
            al    = 1.0 / (1.0 + lcrcp * alpha)
            sd    = al * (qsl_t(index_npt(i,j)) - q(i,j,k))

! Update the water contents

            qcl(i,j,k) = qcl(i,j,k) - sd
            q(i,j,k)   = q(i,j,k)   + sd
            t(i,j,k)   = t(i,j,k)   - sd * lcrcp

          END IF

! Is the relative total humidity less than the threshold amount?
! If so, evaporate all the liquid.

          IF ( rht  <   (rhcrit(irhi,irhj,k)) .AND.                     &
            (cfl(i,j,k)  <=  c_thresh_low) )  THEN
            cfl(i,j,k) = 0.0
            cf(i,j,k)  = cff(i,j,k)
            q(i,j,k)   = q(i,j,k) + qcl(i,j,k)
            t(i,j,k)   = t(i,j,k) - qcl(i,j,k) * lcrcp
            qcl(i,j,k) = 0.0
          END IF

        END IF ! cloud fraction

        IF ((cfl(i,j,k)  <=  c_thresh_low_2)  ) THEN
          cfl(i,j,k) = 0.0
          cf(i,j,k)  = cff(i,j,k)
          q(i,j,k)   = q(i,j,k) + qcl(i,j,k)
          t(i,j,k)   = t(i,j,k) - qcl(i,j,k) * lcrcp
          qcl(i,j,k) = 0.0
        END IF

      END DO

    END DO

  END DO
!$OMP END PARALLEL DO 

! End of the subroutine

  IF (lhook) CALL dr_hook('PC2_CHECKS2',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE pc2_checks2
