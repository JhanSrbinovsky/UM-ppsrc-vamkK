! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Area cloud parameterisation for use with PC2 Cloud Scheme.

SUBROUTINE pc2_hom_arcld(                                               &
!      Pressure related fields
 p_layer_centres, p_layer_boundaries,                                   &
!      Array dimensions
 large_levels, levels_per_level,                                        &
!      Prognostic Fields
 cf_area, t, cf, cfl, cff, q, qcl, qcf,                                 &
!      Logical control
 l_mixing_ratio)

  USE water_constants_mod,  ONLY: lc
  USE atmos_constants_mod,  ONLY: cp, r
  USE yomhook,              ONLY: lhook, dr_hook
  USE parkind1,             ONLY: jprb, jpim
  USE atm_fields_bounds_mod,ONLY: pdims,qdims,tdims

  IMPLICIT NONE

! Description:
!   Area cloud parameterisation for use with PC2:
!   Cusack-like vertical interpolation onto three sub-levels, to obtain
!   increments to P,T,Q,QCL for use with homogeneous forcing routine.
!   Area cloud fraction is then maximum of the 3 bulk values,
!   nothing else is changed.
!
! Method:
!   See the PC2 documentation.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
! Code Description:
!   Language: fortran 77 + common extensions
!   This code is written to UMDP3 v6 programming standards.

!  Global Variables:----------------------------------------------------

!  Subroutine Arguments:------------------------------------------------
  INTEGER ::                                                            &
                        !, INTENT(IN)
   large_levels,                                                        &
!       Total no. of sub-levels being processed by cloud scheme.
!       Currently ((levels - 2)*levels_per_level) + 2
   levels_per_level
!       No. of sub-levels being processed by area cloud scheme.
!       Want an odd number of sublevels per level.
!       NB: levels_per_level = 3 is currently hardwired in the do loops

  LOGICAL ::                                                            &
                           !, INTENT(IN)
   l_mixing_ratio           ! Use mixing ratio formulation

  REAL ::                                                               &
                        !, INTENT(IN)
   p_layer_centres(   pdims%i_start:pdims%i_end,                        & 
                      pdims%j_start:pdims%j_end,                        &  
                                  0:pdims%k_end),                       &
!       pressure at all points, on theta levels (Pa).
   p_layer_boundaries(pdims%i_start:pdims%i_end,                        & 
                      pdims%j_start:pdims%j_end,                        &  
                                  0:pdims%k_end),                       &
!       pressure at all points, on u,v levels (Pa).
   cff(               qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end),                       &
!       Ice cloud fraction (no units)
   qcf(               qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end),                       &
!       Cloud ice content at processed levels (kg water per kg air).
   t(                 tdims%i_start:tdims%i_end,                        & 
                      tdims%j_start:tdims%j_end,                        &  
                                  1:tdims%k_end),                       &
!       Temperature (K)
   cf(                qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end),                       &
!       Total cloud fraction (no units)
   cfl(               qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end),                       &
!       Liquid cloud fraction (no units)
   q(                 qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end),                        &
!       Vapour content (kg water per kg air)
   qcl(               qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end)
!       Liquid content (kg water per kg air)

  REAL ::                                                               &
                        !, INTENT(INOUT)
   cf_area(           qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end)
!       Area cloud fraction

! --------------------------------------------------------------------
! Local variables
! ---------------------------------------------------------------------
  INTEGER :: i,j,k      ! Loop counters:  k - vertical level index.
!                                       i,j - horizontal field index.
  INTEGER :: k_index    ! Extra loop counter for large arrays.

  REAL ::                                                               &
    inverse_level,                                                      &
!       Set to (1. / levels_per_level)
    qt_norm_next,                                                       &
!       Temporary space for qT_norm
    stretcher,                                                          &
    delta_p        
!       Layer pressure thickness * inverse_level

  REAL ::                                                               &
    qsl(              qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end),                       &
!       Saturated specific humidity for temp TL or T.
    tl(               qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end),                       &
!       Liquid temperature (TL) (K).
    qt_norm(          qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end),                       &
!       Total water content normalized to qSAT_WAT.
    p_large(          pdims%i_start:pdims%i_end,                        & 
                      pdims%j_start:pdims%j_end,                        &
                      large_levels),                                    &
!       Values of quantities on large_levels
    t_large(          tdims%i_start:tdims%i_end,                        & 
                      tdims%j_start:tdims%j_end,                        &
                      large_levels),                                    &

    q_large(          qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &

    qcl_large(        qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &

    cf_large(         qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &

    cfl_large(        qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &

    cff_large(        qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &

    dldt_large(       qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &

    dtdt_large(       tdims%i_start:tdims%i_end,                        & 
                      tdims%j_start:tdims%j_end,                        &
                      large_levels),                                    &

    dqdt_large(       qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &

    dpdt_large(       pdims%i_start:pdims%i_end,                        & 
                      pdims%j_start:pdims%j_end,                        &
                      large_levels)

  LOGICAL ::                                                            &
   linked(            qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end)
!       True for sub-layers that have similar supersaturation properties

!  Local parameters and other physical constants------------------------
  REAL, PARAMETER :: lcrcp        = lc/cp 
!       Latent heat of condensation divided by heat capacity of air.
  REAL, PARAMETER  :: drat_thresh = 3.0e-1
!       Test for continuity of sub-levels
  REAL, PARAMETER  :: tol_test    = 1.0e-11
!       Tolerance for non-zero humidities

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL   (KIND=jprb)            :: zhook_handle

! ---------------------------------------------------------------------
! Code starts
! ---------------------------------------------------------------------
  IF (lhook) CALL dr_hook('PC2_HOM_ARCLD',zhook_in,zhook_handle)
  inverse_level = 1. / levels_per_level

! create new arrays for TL and current qcl

  DO k = 1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        tl(i,j,k) = t(i,j,k) - lcrcp*qcl(i,j,k)
      END DO
    END DO
  END DO

! Test for continuity between adjacent layers based on supersaturation
! (qt - qsl) / qsl : as we take differences the - qsl term drops out.
! DEPENDS ON: qsat_wat_mix
  CALL qsat_wat_mix( qsl, tl(qdims%i_start,qdims%j_start,1),            &
           p_layer_centres(qdims%i_start,qdims%j_start,1),              &
           (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start), &
           l_mixing_ratio)

  DO j = qdims%j_start, qdims%j_end
    DO i = qdims%i_start, qdims%i_end
      qt_norm(i,j) =(q(i,j,1)+qcl(i,j,1)                                &
                        +qcf(i,j,1))/qsl(i,j)
    END DO
  END DO

! DEPENDS ON: qsat_wat_mix
  CALL qsat_wat_mix( qsl, tl(qdims%i_start,qdims%j_start,2),            &
           p_layer_centres(qdims%i_start,qdims%j_start,2),              &
           (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start), &
           l_mixing_ratio)

! Do nothing to top and bottom layers

  DO j = qdims%j_start, qdims%j_end
    DO i = qdims%i_start, qdims%i_end

      p_large   (i,j,1) = p_layer_centres(i,j,1)
      t_large   (i,j,1) = t              (i,j,1)
      q_large   (i,j,1) = q              (i,j,1)
      qcl_large (i,j,1) = qcl            (i,j,1)
      cf_large  (i,j,1) = cf             (i,j,1)
      cfl_large (i,j,1) = cfl            (i,j,1)
      cff_large (i,j,1) = cff            (i,j,1)
      dtdt_large(i,j,1) = 0.0
      dqdt_large(i,j,1) = 0.0
      dpdt_large(i,j,1) = 0.0
      dldt_large(i,j,1) = 0.0

      p_large   (i,j,large_levels) = p_layer_centres(i,j,pdims%k_end)
      t_large   (i,j,large_levels) = t              (i,j,tdims%k_end)
      q_large   (i,j,large_levels) = q              (i,j,qdims%k_end)
      qcl_large (i,j,large_levels) = qcl            (i,j,qdims%k_end)
      cf_large  (i,j,large_levels) = cf             (i,j,qdims%k_end)
      cfl_large (i,j,large_levels) = cfl            (i,j,qdims%k_end)
      cff_large (i,j,large_levels) = cff            (i,j,qdims%k_end)
      dtdt_large(i,j,large_levels) = 0.0
      dqdt_large(i,j,large_levels) = 0.0
      dpdt_large(i,j,large_levels) = 0.0
      dldt_large(i,j,large_levels) = 0.0

! Test for continuity (assumed if linked is .true.)
      qt_norm_next  = ( q(i,j,2)+qcl(i,j,2)+qcf(i,j,2) ) / qsl(i,j)
      linked(i,j,1) = (drat_thresh >= ABS(qt_norm(i,j) - qt_norm_next))
      qt_norm(i,j)  = qt_norm_next

    END DO !i
  END DO !j

  DO k = 2, (qdims%k_end - 1)
    k_index = 3 + (levels_per_level * (k-2))

! DEPENDS ON: qsat_wat_mix
    CALL qsat_wat_mix( qsl, tl(qdims%i_start,qdims%j_start,(k+1)),      &
          p_layer_centres(qdims%i_start,qdims%j_start,(k+1)),           &
          (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start),  &
          l_mixing_ratio)

    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
! Test for continuity (assumed if linked = .true.)
        qt_norm_next = (q(i,j,(k+1)) + qcl(i,j,(k+1)) + qcf(i,j,(k+1)) )&
                       / qsl(i,j)
        linked(i,j,k)= (drat_thresh >= ABS(qt_norm(i,j) - qt_norm_next))
        qt_norm(i,j) = qt_norm_next

! Select interpolated pressure levels
        delta_p = (p_layer_boundaries(i,j,(k-1)) -                      &
                   p_layer_boundaries(i,j,k))    * inverse_level
        IF (p_layer_centres(i,j,k) >=                                   &
           (p_layer_boundaries(i,j,k) + delta_p)) THEN
          p_large(i,j,k_index) = p_layer_centres(i,j,k)
        ELSE
          p_large(i,j,k_index) = 0.5*(p_layer_boundaries(i,j,k) +       &
                                 p_layer_boundaries(i,j,(k-1)))
        END IF
        p_large(i,j,(k_index-1)) = p_large(i,j,k_index)
        p_large(i,j,(k_index+1)) = p_large(i,j,k_index)

! Select variable values at layer centres
        t_large   (i,j,k_index) = t(i,j,k)
        q_large   (i,j,k_index) = q(i,j,k)
        qcl_large (i,j,k_index) = qcl(i,j,k)
        cf_large  (i,j,k_index) = cf(i,j,k)
        cfl_large (i,j,k_index) = cfl(i,j,k)
        cff_large (i,j,k_index) = cff(i,j,k)
        dtdt_large(i,j,k_index) = 0.0
        dqdt_large(i,j,k_index) = 0.0
        dpdt_large(i,j,k_index) = 0.0
        dldt_large(i,j,k_index) = 0.0

! Calculate increment in variable values, pressure interpolation
! NB: Using X_large(i,j,(k_index+1)) as store for X increments
! Lsarc_if2:
        IF ( linked(i,j,(k-1)) ) THEN
          IF ( linked(i,j,k) ) THEN
! Interpolate from level k-1 to k+1
            stretcher = delta_p /                                       &
           (p_layer_centres(i,j,k-1)-p_layer_centres(i,j,k+1))

            t_large(i,j,(k_index+1)) = stretcher *                      &
           (t(i,j,(k+1)) - t(i,j,(k-1)))

            q_large(i,j,(k_index+1)) = stretcher *                      &
           (q(i,j,(k+1)) - q(i,j,(k-1)))

            qcl_large(i,j,(k_index+1)) = stretcher *                    &
           (qcl(i,j,(k+1)) - qcl(i,j,(k-1)))

          ELSE
! Interpolate from level k-1 to k
            stretcher = delta_p /                                       &
           (p_layer_centres(i,j,k-1) - p_large(i,j,k_index))

            t_large(i,j,(k_index+1)) = stretcher *                      &
           (t_large(i,j,k_index) - t(i,j,(k-1)))
 
            q_large(i,j,(k_index+1)) = stretcher *                      &
           (q_large(i,j,k_index) - q(i,j,(k-1)))
 
           qcl_large(i,j,(k_index+1)) = stretcher *                    &
           (qcl_large(i,j,k_index) - qcl(i,j,(k-1)))

          END IF

        ELSE
          IF ( linked(i,j,k) ) THEN
! Interpolate from level k to k+1
            stretcher = delta_p /                                       &
            (p_large(i,j,k_index) - p_layer_centres(i,j,k+1))

            t_large(i,j,(k_index+1)) = stretcher *                      &
           (t(i,j,(k+1)) - t_large(i,j,k_index))

            q_large(i,j,(k_index+1)) = stretcher *                      &
           (q(i,j,(k+1)) - q_large(i,j,k_index))

            qcl_large(i,j,(k_index+1)) = stretcher *                    &
           (qcl(i,j,(k+1)) - qcl_large(i,j,k_index))
          ELSE
! No interpolation, freeze at level k
            t_large(i,j,(k_index+1)) = 0.
            q_large(i,j,(k_index+1)) = 0.
            qcl_large(i,j,(k_index+1)) = 0.
          END IF

        END IF

! Protect against q or qcl going negative (T would imply blow-up anyway)
        IF (q_large(i,j,k_index)  <                                     &
                   (ABS(q_large(i,j,(k_index+1)))+tol_test))            &
                        q_large(i,j,(k_index+1)) = 0.
 
        IF (qcl_large(i,j,k_index)  <                                   &
                     (ABS(qcl_large(i,j,(k_index+1)))+tol_test))        &
                          qcl_large(i,j,(k_index+1)) = 0.

! Select variable values at level below layer centre
        t_large   (i,j,(k_index-1)) = t_large(i,j,k_index)
        q_large   (i,j,(k_index-1)) = q_large(i,j,k_index)
        qcl_large (i,j,(k_index-1)) = qcl_large(i,j,k_index)
        cf_large  (i,j,(k_index-1)) = cf(i,j,k)
        cfl_large (i,j,(k_index-1)) = cfl(i,j,k)
        cff_large (i,j,(k_index-1)) = cff(i,j,k)
        dtdt_large(i,j,(k_index-1)) = -t_large(i,j,(k_index+1))
        dqdt_large(i,j,(k_index-1)) = -q_large(i,j,(k_index+1))
        dpdt_large(i,j,(k_index-1)) = delta_p
        dldt_large(i,j,(k_index-1)) = -qcl_large(i,j,(k_index+1))

! Select variable values at level above layer centre
! NB: CEASE using X_large(i,j,(k_index+1)) as store for X increments.
        dtdt_large(i,j,(k_index+1)) = t_large(i,j,(k_index+1))
        dqdt_large(i,j,(k_index+1)) = q_large(i,j,(k_index+1))
        dpdt_large(i,j,(k_index+1)) = -delta_p
        dldt_large(i,j,(k_index+1)) = qcl_large(i,j,(k_index+1))
        t_large   (i,j,(k_index+1)) = t_large(i,j,k_index)
        q_large   (i,j,(k_index+1)) = q_large(i,j,k_index)
        qcl_large (i,j,(k_index+1)) = qcl_large(i,j,k_index)
        cf_large  (i,j,(k_index+1)) = cf(i,j,k)
        cfl_large (i,j,(k_index+1)) = cfl(i,j,k)
        cff_large (i,j,(k_index+1)) = cff(i,j,k)

      END DO !i
    END DO !j
  END DO !k

! DEPENDS ON: pc2_homog_plus_turb
  CALL pc2_homog_plus_turb(p_large,large_levels, 0.0,                   &
    t_large,cf_large,cfl_large,cff_large,q_large,qcl_large,             &
    dtdt_large,dqdt_large,dldt_large,dpdt_large,                        &
    0.0,0.0,l_mixing_ratio)

  DO j = qdims%j_start, qdims%j_end
    DO i = qdims%i_start, qdims%i_end

      cf_area     (i,j,1) = cf_large(i,j,1)
      cf_area(i,j,qdims%k_end) = cf_large(i,j,large_levels)

! Check CF_area isn't greater than 1 or less than 0
      cf_area          (i,j,1) = MAX( MIN( cf_area(i,j,1), 1.0), 0.0)
      cf_area(i,j,qdims%k_end) = MAX(                                   &
         MIN( cf_area(i,j,qdims%k_end), 1.0 ) ,0.0 )

    END DO
  END DO

! Output variables for remaining layers
  DO k = 2, (qdims%k_end - 1)
    k_index = 3 + (levels_per_level * (k-2))
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
! Area cloud fraction is maximum of sub-layer cloud fractions
        cf_area(i,j,k) =                                                &
        MAX( cf_large(i,j,k_index),                                     &
             (MAX(cf_large(i,j,(k_index+1)),                            &
                  cf_large(i,j,(k_index-1)))) )

! Check CF_area isn't greater than 1 or less than 0
        cf_area(i,j,k) = MAX(MIN(cf_area(i,j,k),1.0),0.0)
      END DO !i
     END DO !j
  END DO !k

  IF (lhook) CALL dr_hook('PC2_HOM_ARCLD',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE pc2_hom_arcld
! ======================================================================
