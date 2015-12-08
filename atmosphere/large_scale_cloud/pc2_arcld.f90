! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Area cloud parameterisation for use with PC2 Cloud Scheme.

SUBROUTINE pc2_arcld(                                                   &
!      Pressure related fields
 p_layer_centres, p_layer_boundaries, ccb, cumulus, rhcrit,             &
!      Array dimensions
 rhc_row_length, rhc_rows,                                              &
 large_levels, levels_per_level,                                        &
!      Prognostic Fields
 cf_area, t, cf, cfl, cff, q, qcl, qcf, rhts,                           &
 tlts, qtts, ptts,                                                      &
!      Logical control
 l_mixing_ratio)

  USE water_constants_mod,   ONLY: lc
  USE atmos_constants_mod,   ONLY: cp
  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim
  USE atm_fields_bounds_mod, ONLY: qdims, pdims, tdims

  IMPLICIT NONE

! Description:
!   Cusack-like vertical interpolation onto 3 sub-levels to calculate
!   cloud fraction and condensate in PC2 initiation. Also area cloud
!   parameterisation for use with radiation scheme.
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
   rhc_row_length, rhc_rows,                                            &
!       Dimensions of the rhcrit variable.
   ccb(               qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end),                       &
!       convective cloud base
   large_levels,                                                        &
!       Total no. of sub-levels being processed by cloud scheme.
!       Currently ((levels - 2)*levels_per_level) + 2
   levels_per_level
!       No. of sub-levels being processed by area cloud scheme.
!       Want an odd number of sublevels per level.
!       NB: levels_per_level = 3 is currently hardwired in the do loops

  LOGICAL ::                                                            &
                        !, INTENT(IN)
   cumulus(           qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end),                       &
!       Is this a boundary layer cumulus point
   l_mixing_ratio  
!       Use mixing ratio formulation

  REAL ::                                                               &
                        !, INTENT(IN)
   p_layer_centres(   pdims%i_start:pdims%i_end,                        & 
                      pdims%j_start:pdims%j_end,                        &  
                      0:pdims%k_end),                                   &
!       Pressure at all points, on theta levels (Pa).
   p_layer_boundaries(pdims%i_start:pdims%i_end,                        & 
                      pdims%j_start:pdims%j_end,                        &  
                      0:pdims%k_end),                                   &
!       Pressure at all points, on u,v levels (Pa).
   rhcrit(            rhc_row_length,                                   &
                      rhc_rows,                                         &
                                  1:qdims%k_end),                       &
!       Critical relative humidity.  See the the paragraph incorporating
!       eqs P292.11 to P292.14; the values need to be tuned for the give
!       set of levels.
   cff(               qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end),                       &
!       Ice cloud fraction (no units)
   qcf(               qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end)
!       Cloud ice content at processed levels (kg water per kg air).

  REAL ::                                                               &
                        !, INTENT(INOUT)
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
                                  1:qdims%k_end),                       &
!       Vapour content (kg water per kg air)
   qcl(               qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end),                       &
!       Liquid content (kg water per kg air)
   rhts(              qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end),                       &
!       Variable carrying initial RHT wrt TL from start of timestep
   tlts(              qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end),                       &
!       TL at start of timestep
   qtts(              qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end),                       &
!       qT at start of timestep
   ptts(              pdims%i_start:pdims%i_end,                        & 
                      pdims%j_start:pdims%j_end,                        &  
                      pdims%k_start:pdims%k_end),                       &
!       Pressure at theta levels at start of timestep
   cf_area(           qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end) 
!       Area cloud fraction

! --------------------------------------------------------------------
! Local variables
! ---------------------------------------------------------------------
  INTEGER :: i,j,k   ! Loop counters: k   - vertical level index.
!                                     i,j - horizontal field index.
  INTEGER :: k_index ! Extra loop counter for large arrays.

  REAL ::                                                               &
    inverse_level,                                                      &
                     ! Set to (1. / levels_per_level)
    qt_norm_next,                                                       &
                     ! Temporary space for qT_norm
    stretcher,                                                          &
    delta_p          ! Layer pressure thickness * inverse_level

  REAL ::                                                               &
    qsl(              qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end),                       &
!       Saturated specific humidity for temp TL or T.
    qcl_latest(       qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end),                       &
!       Cloud liquid content at processed levels (kg water per kg air).
    tl(               qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end),                       &
!       Liquid temperature (TL) (K).
    qt_norm(          qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end),                       &
!       Total water content normalized to qSAT_WAT.
    rhcrit_large(rhc_row_length,rhc_rows,large_levels),                 &
!
!       Values of quantities on large_levels
    p_large(          pdims%i_start:pdims%i_end,                        & 
                      pdims%j_start:pdims%j_end,                        &
                      large_levels),                                    &
!
    t_large(          tdims%i_start:tdims%i_end,                        & 
                      tdims%j_start:tdims%j_end,                        &
                      large_levels),                                    &
!
    q_large(          qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &
!
    qcl_large(        qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &
!
    cf_large(         qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &
!
    cfl_large(        qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &
!
    cff_large(        qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &
!
    rhts_large(       qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &
!
    tlts_large(       qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &
!
    qtts_large(       qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &
                      large_levels),                                    &
!
    ptts_large(       pdims%i_start:pdims%i_end,                        & 
                      pdims%j_start:pdims%j_end,                        &
                      large_levels)

  LOGICAL ::                                                            &
   linked(            qdims%i_start:qdims%i_end,                        & 
                      qdims%j_start:qdims%j_end,                        &  
                                  1:qdims%k_end)
!       True for sub-layers that have similar supersaturation properties

!  Local parameters and other physical constants------------------------
  REAL, PARAMETER :: lcrcp=lc/cp
!       Latent heat of condensation divided by heat capacity of air.
  REAL, PARAMETER :: drat_thresh =3.0e-1
!       Test for continuity of sub-levels
  REAL, PARAMETER :: tol_test  =1.0e-11
!       Tolerance for non-zero humidities

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL   (KIND=jprb)            :: zhook_handle

! ---------------------------------------------------------------------
! Code starts
! ---------------------------------------------------------------------
  IF (lhook) CALL dr_hook('PC2_ARCLD',zhook_in,zhook_handle)
  inverse_level = 1. / levels_per_level

! Create new arrays for TL and current qcl

  DO k = 1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        tl(i,j,k)         = t(i,j,k) - lcrcp*qcl(i,j,k)
        qcl_latest(i,j,k) = qcl(i,j,k)
      END DO !i
    END DO !j
  END DO !k

! Test for continuity between adjacent layers based on supersaturation
! (qt - qsl) / qsl : as we take differences the - qsl term drops out.

! DEPENDS ON: qsat_wat_mix
  CALL qsat_wat_mix( qsl, tl(qdims%i_start,qdims%j_start,1),            &
           p_layer_centres(qdims%i_start,qdims%j_start,1),              &
           (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start), &
           l_mixing_ratio)

  DO j = qdims%j_start, qdims%j_end
    DO i = qdims%i_start, qdims%i_end
      qt_norm(i,j) =(q(i,j,1)+qcl_latest(i,j,1)+qcf(i,j,1))/qsl(i,j)
    END DO
  END DO

! DEPENDS ON: qsat_wat_mix
  CALL qsat_wat_mix( qsl, tl(qdims%i_start,qdims%j_start,2),            &
           p_layer_centres(qdims%i_start,qdims%j_start,2),              &
           (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start), &
           l_mixing_ratio)

! Do nothing to top and bottom layers
  DO j = 1, rhc_rows
    DO i = 1, rhc_row_length
      rhcrit_large(i,j,1)            = rhcrit(i,j,1)
      rhcrit_large(i,j,large_levels) = rhcrit(i,j,qdims%k_end)
    END DO
  END DO

  DO j = qdims%j_start, qdims%j_end
    DO i = qdims%i_start, qdims%i_end

      p_large   (i,j,1) = p_layer_centres(i,j,1)
      t_large   (i,j,1) = t(i,j,1)
      q_large   (i,j,1) = q(i,j,1)
      qcl_large (i,j,1) = qcl_latest(i,j,1)
      cf_large  (i,j,1) = cf(i,j,1)
      cfl_large (i,j,1) = cfl(i,j,1)
      cff_large (i,j,1) = cff(i,j,1)
      tlts_large(i,j,1) = tlts(i,j,1)
      qtts_large(i,j,1) = qtts(i,j,1)
      ptts_large(i,j,1) = ptts(i,j,1)
      rhts_large(i,j,1) = rhts(i,j,1)

      p_large   (i,j,large_levels) = p_layer_centres(i,j,pdims%k_end)
      t_large   (i,j,large_levels) = t(i,j,tdims%k_end)
      q_large   (i,j,large_levels) = q(i,j,qdims%k_end)
      qcl_large (i,j,large_levels) = qcl_latest(i,j,qdims%k_end)
      cf_large  (i,j,large_levels) = cf(i,j,qdims%k_end)
      cfl_large (i,j,large_levels) = cfl(i,j,qdims%k_end)
      cff_large (i,j,large_levels) = cff(i,j,qdims%k_end)
      tlts_large(i,j,large_levels) = tlts(i,j,large_levels)
      qtts_large(i,j,large_levels) = qtts(i,j,large_levels)
      ptts_large(i,j,large_levels) = ptts(i,j,large_levels)
      rhts_large(i,j,large_levels) = rhts(i,j,large_levels)

! Test for continuity (assumed if linked is .true.)
      qt_norm_next=(q(i,j,2)+qcl_latest(i,j,2)+qcf(i,j,2))/qsl(i,j)
      linked(i,j,1) =                                                   &
            (drat_thresh >= ABS(qt_norm(i,j) - qt_norm_next))
      qt_norm(i,j) = qt_norm_next

    END DO !i
  END DO ! j

  DO k = 2, (qdims%k_end - 1)
    k_index = 3 + (levels_per_level * (k-2))

! DEPENDS ON: qsat_wat_mix
    CALL qsat_wat_mix( qsl, tl(qdims%i_start,qdims%j_start,(k+1)),      &
           p_layer_centres(qdims%i_start,qdims%j_start,(k+1)),          &
           (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start), &
           l_mixing_ratio)

! Select associated rhcrit values
    DO j = 1, rhc_rows
      DO i = 1, rhc_row_length
        rhcrit_large(i,j,k_index-1) = rhcrit(i,j,k)
        rhcrit_large(i,j,k_index)   = rhcrit(i,j,k)
        rhcrit_large(i,j,k_index+1) = rhcrit(i,j,k)
      END DO
    END DO

    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
! Test for continuity (assumed if linked = .true.)
        qt_norm_next  = (q(i,j,(k+1)) + qcl_latest(i,j,(k+1))           &
                      + qcf(i,j,(k+1))) / qsl(i,j)
        linked(i,j,k) =                                                 &
            (drat_thresh >= ABS(qt_norm(i,j) - qt_norm_next))
        qt_norm(i,j)  = qt_norm_next
!
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
        p_large(i,j,(k_index-1)) = p_large(i,j,k_index)+delta_p
        p_large(i,j,(k_index+1)) = p_large(i,j,k_index)-delta_p
!
! Select variable values at layer centres
        t_large   (i,j,k_index) = t(i,j,k)
        q_large   (i,j,k_index) = q(i,j,k)
        qcl_large (i,j,k_index) = qcl_latest(i,j,k)
        cf_large  (i,j,k_index) = cf(i,j,k)
        cfl_large (i,j,k_index) = cfl(i,j,k)
        cff_large (i,j,k_index) = cff(i,j,k)
        tlts_large(i,j,k_index) = tlts(i,j,k)
        qtts_large(i,j,k_index) = qtts(i,j,k)
        ptts_large(i,j,k_index) = ptts(i,j,k)
        rhts_large(i,j,k_index) = rhts(i,j,k)

! Calculate increment in variable values, pressure interpolation
! NB: Using X_large(i,j,(k_index+1)) as store for X increments
        IF ( linked(i,j,(k-1)) ) THEN
          IF ( linked(i,j,k) ) THEN
!               Interpolate from level k-1 to k+1
            stretcher = delta_p /                                       &
           (p_layer_centres(i,j,k-1)-p_layer_centres(i,j,k+1))
            t_large(i,j,(k_index+1)) = stretcher *                      &
           (t(i,j,(k+1)) - t(i,j,(k-1)))
            q_large(i,j,(k_index+1)) = stretcher *                      &
           (q(i,j,(k+1)) - q(i,j,(k-1)))
            qcl_large(i,j,(k_index+1)) = stretcher *                    &
           (qcl_latest(i,j,(k+1)) - qcl_latest(i,j,(k-1)))
            tlts_large(i,j,(k_index+1)) = stretcher *                   &
           (tlts(i,j,(k+1)) - tlts(i,j,(k-1)))
            qtts_large(i,j,(k_index+1)) = stretcher *                   &
           (qtts(i,j,(k+1)) - qtts(i,j,(k-1)))
            ptts_large(i,j,(k_index+1)) = stretcher *                   &
           (ptts(i,j,(k+1)) - ptts(i,j,(k-1)))
          ELSE
!               Interpolate from level k-1 to k
            stretcher = delta_p /                                       &
           (p_layer_centres(i,j,k-1) - p_large(i,j,k_index))
            t_large(i,j,(k_index+1)) = stretcher *                      &
           (t_large(i,j,k_index) - t(i,j,(k-1)))
            q_large(i,j,(k_index+1)) = stretcher *                      &
           (q_large(i,j,k_index) - q(i,j,(k-1)))
            qcl_large(i,j,(k_index+1)) = stretcher *                    &
           (qcl_large(i,j,k_index) - qcl_latest(i,j,(k-1)))
            tlts_large(i,j,(k_index+1)) = stretcher *                   &
           (tlts_large(i,j,k_index) - tlts(i,j,(k-1)))
            qtts_large(i,j,(k_index+1)) = stretcher *                   &
           (qtts_large(i,j,k_index) - qtts(i,j,(k-1)))
            ptts_large(i,j,(k_index+1)) = stretcher *                   &
           (ptts_large(i,j,k_index) - ptts(i,j,(k-1)))
          END IF

        ELSE
          IF ( linked(i,j,k) ) THEN
!               Interpolate from level k to k+1
            stretcher = delta_p /                                       &
            (p_large(i,j,k_index) - p_layer_centres(i,j,k+1))
            t_large(i,j,(k_index+1)) = stretcher *                      &
           (t(i,j,(k+1)) - t_large(i,j,k_index))
            q_large(i,j,(k_index+1)) = stretcher *                      &
           (q(i,j,(k+1)) - q_large(i,j,k_index))
            qcl_large(i,j,(k_index+1)) = stretcher *                    &
           (qcl_latest(i,j,(k+1)) - qcl_large(i,j,k_index))
            tlts_large(i,j,(k_index+1)) = stretcher *                   &
           (tlts(i,j,(k+1)) - tlts_large(i,j,k_index))
            qtts_large(i,j,(k_index+1)) = stretcher *                   &
           (qtts(i,j,(k+1)) - qtts_large(i,j,k_index))
            ptts_large(i,j,(k_index+1)) = stretcher *                   &
           (ptts(i,j,(k+1)) - ptts_large(i,j,k_index))
          ELSE
!               No interpolation, freeze at level k
            t_large   (i,j,(k_index+1)) = 0.
            q_large   (i,j,(k_index+1)) = 0.
            qcl_large (i,j,(k_index+1)) = 0.
            tlts_large(i,j,(k_index+1)) = 0.
            qtts_large(i,j,(k_index+1)) = 0.
            ptts_large(i,j,(k_index+1)) = 0.
          END IF

        END IF

! Protect against q or qcl going negative (T would imply blow-up anyway)
        IF (q_large(i,j,k_index)  <                                     &
                   (ABS(q_large(i,j,(k_index+1)))+tol_test))            &
                        q_large(i,j,(k_index+1)) = 0.
        IF (qcl_large(i,j,k_index)  <                                   &
                     (ABS(qcl_large(i,j,(k_index+1)))+tol_test))        &
                          qcl_large(i,j,(k_index+1)) = 0.
        IF (qtts_large(i,j,k_index) <                                   &
                    (ABS(qtts_large(i,j,(k_index+1)))+tol_test))        &
                         qtts_large(i,j,(k_index+1)) = 0.

! Select variable values at level below layer centre
        t_large   (i,j,(k_index-1)) = t_large(i,j,k_index) -            &
                                        t_large(i,j,(k_index+1))
        q_large   (i,j,(k_index-1)) = q_large(i,j,k_index) -            &
                                        q_large(i,j,(k_index+1))
        qcl_large (i,j,(k_index-1)) = qcl_large(i,j,k_index)-           &
                                        qcl_large(i,j,(k_index+1))
        cf_large  (i,j,(k_index-1)) = cf(i,j,k)
        cfl_large (i,j,(k_index-1)) = cfl(i,j,k)
        cff_large (i,j,(k_index-1)) = cff(i,j,k)
        tlts_large(i,j,(k_index-1)) = tlts_large(i,j,k_index) -         &
                                        tlts_large(i,j,(k_index+1))
        qtts_large(i,j,(k_index-1)) = qtts_large(i,j,k_index) -         &
                                        qtts_large(i,j,(k_index+1))
        ptts_large(i,j,(k_index-1)) = ptts_large(i,j,k_index) -         &
                                        ptts_large(i,j,(k_index+1))

! Select variable values at level above layer centre
! NB: CEASE using X_large(i,j,(k_index+1)) as store for X increments.
        t_large   (i,j,(k_index+1)) = t_large(i,j,(k_index+1)) +        &
                                        t_large(i,j,k_index)
        q_large   (i,j,(k_index+1)) = q_large(i,j,(k_index+1)) +        &
                                        q_large(i,j,k_index)
        qcl_large (i,j,(k_index+1)) = qcl_large(i,j,(k_index+1)) +      &
                                        qcl_large(i,j,k_index)
        cf_large  (i,j,(k_index+1)) = cf(i,j,k)
        cfl_large (i,j,(k_index+1)) = cfl(i,j,k)
        cff_large (i,j,(k_index+1)) = cff(i,j,k)
        tlts_large(i,j,(k_index+1)) = tlts_large(i,j,(k_index+1)) +     &
                                        tlts_large(i,j,k_index)
        qtts_large(i,j,(k_index+1)) = qtts_large(i,j,(k_index+1)) +     &
                                        qtts_large(i,j,k_index)
        ptts_large(i,j,(k_index+1)) = ptts_large(i,j,(k_index+1)) +     &
                                        ptts_large(i,j,k_index)

      END DO
    END DO

! Calculate RH above and below layer centres
! DEPENDS ON: qsat_wat_mix
    CALL qsat_wat_mix(qsl,                                              &
           tlts_large(qdims%i_start,qdims%j_start,k_index-1),           &
           ptts_large(qdims%i_start,qdims%j_start,k_index-1),           &
           (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start), &
           l_mixing_ratio)

    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        rhts_large(i,j,k_index-1)=qtts_large(i,j,k_index-1)             &
                                    /qsl(i,j)
      END DO
    END DO

! DEPENDS ON: qsat_wat_mix
    CALL qsat_wat_mix(qsl,                                              &
           tlts_large(qdims%i_start,qdims%j_start,k_index+1),           &
           ptts_large(qdims%i_start,qdims%j_start,k_index+1),           &
           (1+qdims%i_end-qdims%i_start)*(1+qdims%j_end-qdims%j_start), &
           l_mixing_ratio)

    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        rhts_large(i,j,k_index+1)=qtts_large(i,j,k_index+1)             &
                                    /qsl(i,j)
      END DO
    END DO

  END DO

! DEPENDS ON: pc2_initiate
  CALL pc2_initiate(p_large,ccb,cumulus,rhcrit_large,                   &
    large_levels, rhc_row_length,rhc_rows,                              &
    t_large,cf_large,cfl_large,cff_large,q_large,qcl_large,             &
    rhts_large,l_mixing_ratio)

  DO j = qdims%j_start, qdims%j_end
    DO i = qdims%i_start, qdims%i_end
      t         (i,j,1)           = t_large(i,j,1)
      q         (i,j,1)           = q_large(i,j,1)

      cf_area   (i,j,1)           = cf_large(i,j,1)
      cf        (i,j,1)           = cf_large(i,j,1)

      qcl_latest(i,j,1)           = qcl_large(i,j,1)
      cfl       (i,j,1)           = cfl_large(i,j,1)

      t         (i,j,qdims%k_end) = t_large(i,j,large_levels)
      q         (i,j,qdims%k_end) = q_large(i,j,large_levels)

      cf_area   (i,j,qdims%k_end) = cf_large(i,j,large_levels)
      cf        (i,j,qdims%k_end) = cf_large(i,j,large_levels)

      qcl_latest(i,j,qdims%k_end) = qcl_large(i,j,large_levels)
      cfl       (i,j,qdims%k_end) = cfl_large(i,j,large_levels)
    END DO
  END DO

! Output variables for remaining layers
  DO k = 2, (qdims%k_end - 1)

    k_index = 3 + (levels_per_level * (k-2))

    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
! Area cloud fraction is maximum of sub-layer cloud fractions
        cf_area(i,j,k) =                                                &
              MAX(cf_large(i,j,k_index),                                &
             (MAX(cf_large(i,j,(k_index+1)),                            &
                  cf_large(i,j,(k_index-1)))) )

! Bulk cloud fraction is mean of sub-layer cloud fractions : strictly
! this is a pressure weighted mean being used to approximate a volume
! mean. Over the depth of a layer the difference should not be large.
        cf(i,j,k) = inverse_level *                                     &
           ( cf_large(i,j,(k_index-1)) +                                &
             cf_large(i,j, k_index)    +                                &
             cf_large(i,j,(k_index+1)) )

! The pressure weighted mean of qcf is the input qcf: do not update.

! Qcl is the pressure weighted mean of qcl from each sub-layer.
        qcl_latest(i,j,k) = inverse_level *                             &
      ( qcl_large(i,j,(k_index-1)) +                                    &
        qcl_large(i,j,k_index) + qcl_large(i,j,(k_index+1)) )

! Liq. cloud fraction is mean of sub-layer cloud fractions : strictly
! this is a pressure weighted mean being used to approximate a volume
! mean. Over the depth of a layer the difference should not be large.
        cfl(i,j,k) = inverse_level *                                    &
      ( cfl_large(i,j,(k_index-1)) +                                    &
        cfl_large(i,j,k_index)     +                                    &
        cfl_large(i,j,(k_index+1)) )

! Update Q
! Update T
! Move qcl_latest into qcl.
        q(i,j,k)   = q(i,j,k) + qcl(i,j,k) - qcl_latest(i,j,k)
        t(i,j,k)   = t(i,j,k) - (qcl(i,j,k)*lcrcp) +                    &
                                (qcl_latest(i,j,k) * lcrcp)
        qcl(i,j,k) = qcl_latest(i,j,k)

      END DO
    END DO

  END DO

  IF (lhook) CALL dr_hook('PC2_ARCLD',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE pc2_arcld
! ======================================================================
