! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ PC2 Cloud Scheme: Homogeneous forcing and Turbulence
! Subroutine Interface:
SUBROUTINE pc2_homog_plus_turb(                                         &
!   Pressure related fields
 p_theta_levels,                                                        &
!   Array dimensions
 nlevels,                                                               &
!   Timestep
 timestep,                                                              &
!   Prognostic Fields
 t, cf, cfl, cff, q, qcl,                                               &
!   Forcing quantities for driving the homogeneous forcing
 dtdt, dqdt, dldt, dpdt,                                                &
!   Other quantities for the turbulence
 dbsdtbs0, dbsdtbs1,                                                    &
!   Model switches
 l_mixing_ratio)

  USE water_constants_mod,   ONLY: lc
  USE pc2_constants_mod,     ONLY: cloud_rounding_tol
  USE atmos_constants_mod,   ONLY: cp, r, repsilon
  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim
  USE atm_fields_bounds_mod, ONLY: pdims, qdims, tdims
  USE pc2_constants_mod,     ONLY: dbsdtbs_exp, pdf_power,              &
                                   pdf_merge_power
  USE cloud_inputs_mod,      ONLY: l_fixbug_pc2_qcl_incr,l_fixbug_pc2_mixph

  IMPLICIT NONE

! Purpose:
!   This subroutine calculates the change in liquid content, liquid
!   cloud fraction and total cloud fraction as a result of homogeneous
!   forcing of the gridbox with temperature, pressure, vapour and liquid
!   increments and turbulent mixing effects.

! Method:
!   Uses the method in Gregory et al (2002, QJRMS 128 1485-1504) and
!   Wilson and Gregory (2003, QJRMS 129 967-986)
!   which considers a probability density distribution whose
!   properties are only influenced by a change of width due to
!   turbulence.
!   There is a check to ensure we cannot remove more condensate than was 
!   there to start with.
!   There is a check to ensure we remove all liquid if all 
!   fraction is removed.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
! Documentation: PC2 Cloud Scheme Documentation

!  Subroutine Arguments:------------------------------------------------
  INTEGER ::                                                            &
                        !, INTENT(IN)
   nlevels
!    No. of levels being processed.

  REAL ::                                                               &
                        !, INTENT(IN)
   timestep,                                                            &
!    Model timestep (s)
   p_theta_levels(pdims%i_start:pdims%i_end,                            & 
                  pdims%j_start:pdims%j_end,                            &  
                  nlevels),                                             &
!    pressure at all points (Pa)
   cff(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!    Ice cloud fraction (no units)
   dtdt(          tdims%i_start:tdims%i_end,                            & 
                  tdims%j_start:tdims%j_end,                            &  
                  nlevels),                                             &
!    Increment of temperature from forcing mechanism (K)
   dqdt(          qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!    Increment of vapour from forcing mechanism (kg kg-1)
   dldt(          qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!    Increment of liquid from forcing mechanism (kg kg-1)
   dpdt(          pdims%i_start:pdims%i_end,                            & 
                  pdims%j_start:pdims%j_end,                            &  
                  nlevels),                                             &
!    Increment in pressure from forcing mechanism (Pa)
   dbsdtbs0,                                                            &
!    Value of dbs/dt / bs which is independent of forcing (no units)
   dbsdtbs1
!    Value of dbs/dt / bs which is proportional to the forcing
!    ( (kg kg-1 s-1)-1 )

  LOGICAL ::                                                            &
                         !, INTENT(IN)
   l_mixing_ratio        ! Use mixing ratio formulation

  REAL ::                                                               &
                        !, INTENT(INOUT)
   t(  tdims%i_start:tdims%i_end,                                       & 
       tdims%j_start:tdims%j_end,                                       &  
       nlevels),                                                        &
!    Temperature (K)
   cf( qdims%i_start:qdims%i_end,                                       & 
       qdims%j_start:qdims%j_end,                                       &  
       nlevels),                                                        &
!    Total cloud fraction (no units)
   cfl(qdims%i_start:qdims%i_end,                                       & 
       qdims%j_start:qdims%j_end,                                       &  
       nlevels),                                                        &
!    Liquid cloud fraction (no units)
   q(  qdims%i_start:qdims%i_end,                                       & 
       qdims%j_start:qdims%j_end,                                       &  
       nlevels),                                                        &
!    Vapour content (kg water per kg air)
   qcl(qdims%i_start:qdims%i_end,                                       & 
       qdims%j_start:qdims%j_end,                                       &  
       nlevels)
!    Liquid content (kg water per kg air)

!    External functions:

!    Local parameters and other physical constants-----------------------
  REAL, PARAMETER :: lcrcp=lc/cp
!    Latent heat of condensation divided by heat capacity of air.
  REAL, PARAMETER :: b_factor=(pdf_power+1.0)/(pdf_power+2.0)
!    Premultiplier to calculate the amplitude of the probability
!    density function at the saturation boundary (G).
  REAL, PARAMETER :: smallp=1.0e-10
!    Small positive value for use in if tests

!  Local scalars--------------------------------------------------------
  REAL ::                                                               &
   alpha,                                                               &
                  ! Rate of change of saturation specific humidity with
                  ! temperature calculated at dry-bulb temperature
                  ! (kg kg-1 K-1)
   alpha_p,                                                             &
                  ! Rate of change of saturation specific humidity with
                  ! pressure calculated at dry-bulb temperature (Pa K-1)
   al,                                                                  &
                  ! 1 / (1 + alpha L/cp)  (no units)
   c_1,                                                                 &
                  ! Mid-timestep liquid cloud fraction (no units)
   dbsdtbs,                                                             &
                  ! Relative rate of change of distribution width (s-1)
   dqcdt,                                                               &
                  ! Forcing of QC (kg kg-1 s-1)
   deltal,                                                              &
                  ! Change in liquid content (kg kg-1)
   g,                                                                   &
                  ! Amplitude of the probability density function at
                  ! the saturation boundary (kg kg-1)-1
   qc,                                                                  &
                  ! aL (q + l - qsat(TL) )  (kg kg-1)
   sd             ! Saturation deficit (= aL (q - qsat(T)) )  (kg kg-1)

!  (b) Others.
  INTEGER :: k,i,j,                                                     &
                  ! Loop counters: K - vertical level index
                  ! I,J - horizontal position index
             npt  ! Number of point on which to perform calculations

!  Local arrays---------------------------------------------------------
  REAL ::                                                               &
   qsl_t(     (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
!    Saturated specific humidity for dry bulb temperature T
   qsl_tl(    (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
!    Saturated specific humidity for liquid temperature TL
   tl_c(      (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
!    Liquid temperature (= T - L/cp QCL)  (K)
   t_c(       (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
   p_c(       (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
!    Temperature and pressure on compressed points for qsat call
   cf_c(      (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
!    Total cloud fraction on compressed points
   cfl_c(     (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
!    Liquid cloud fraction on compressed points
   cff_c(     (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
!    Ice cloud fraction on compressed points
   deltac_c(  (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
!    Change in total cloud fraction (no units)
   deltacl_c( (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
!    Change in liquid cloud fraction (no units)
   deltacf_c( (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) )
!    Change in ice cloud fraction (no units)
  INTEGER ::                                                            &
   index_npt(qdims%i_start:qdims%i_end,                                 & 
             qdims%j_start:qdims%j_end)


  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL   (KIND=jprb)            :: zhook_handle

!- End of Header

  IF (lhook) CALL dr_hook('PC2_HOMOG_PLUS_TURB',zhook_in,zhook_handle)

! ==Main Block==--------------------------------------------------------



! Loop round levels to be processed

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(cfl_c,      &
!$OMP& index_npt, npt, cf_c, cff_c, t_c, tl_c, p_c, deltacf_c, qsl_t,   &
!$OMP& qsl_tl, alpha, al, alpha_p, sd, g, dqcdt, dbsdtbs, qc, deltal, i,&
!$OMP& j, k, c_1, deltacl_c)
  DO k = 1, nlevels

! Initialisation of change in frozen cloud fraction, which is always
! zero from this routine.
    deltacf_c(:) = 0.0

! Copy points into compressed arrays
    npt = 0
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        IF (cfl(i,j,k) > cloud_rounding_tol ) THEN
          npt = npt + 1
          index_npt(i,j) = npt
          cfl_c(npt) = cfl(i,j,k)
          cf_c(npt)  = cf(i,j,k)
          cff_c(npt) = cff(i,j,k)
          t_c(npt) = t(i,j,k)
          tl_c(npt) = t(i,j,k)-lcrcp*qcl(i,j,k)
          p_c(npt) = p_theta_levels(i,j,k)
        END IF
      END DO
    END DO

! ----------------------------------------------------------------------
! 2. Calculate Saturated Specific Humidity with respect to liquid water
!    for both dry bulb and wet bulb temperatures.
! ----------------------------------------------------------------------
    IF (npt > 0) THEN
! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(qsl_t, t_c, p_c, npt, l_mixing_ratio)
! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(qsl_tl, tl_c, p_c, npt, l_mixing_ratio)
    END IF

    DO j = qdims%j_start, qdims%j_end

      DO i = qdims%i_start, qdims%i_end

! There is no need to perform the total cloud fraction calculation in
! this subroutine if there is no, or full, liquid cloud cover.

        IF (cfl(i,j,k) >        cloud_rounding_tol .AND.                &
            cfl(i,j,k) < (1.0 - cloud_rounding_tol)) THEN

! ----------------------------------------------------------------------
! 3. Calculate the parameters relating to the probability density func.
! ----------------------------------------------------------------------

! Need to estimate the rate of change of saturated specific humidity
! with respect to temperature (alpha) first, then use this to calculate
! factor aL. Also estimate the rate of change of qsat with pressure.
          alpha=repsilon*lc*qsl_t(index_npt(i,j))/(r*t_c(index_npt(i,j))**2)
          al=1.0/(1.0+lcrcp*alpha)
          alpha_p = -qsl_t(index_npt(i,j))/p_c(index_npt(i,j))

! Calculate the saturation deficit SD

          sd=al*(qsl_t(index_npt(i,j))-q(i,j,k))

! Calculate the amplitude of the probability density function at the
! saturation boundary.

          IF (qcl(i,j,k) > smallp .AND. sd > smallp) THEN

            g=b_factor*(   cfl(i,j,k)**pdf_merge_power                  &
             *(1.0-cfl(i,j,k))**2/sd                                    &
             +         (1.0-cfl(i,j,k))**pdf_merge_power                &
             *cfl(i,j,k)**2/qcl(i,j,k)  )                               &
             /(   cfl(i,j,k)**pdf_merge_power+(1.0-cfl(i,j,k))          &
             **pdf_merge_power   )

          ELSE
            g=0.0
          END IF

! Calculate the rate of change of Qc due to the forcing

          dqcdt=al * ( dqdt(i,j,k)-alpha*dtdt(i,j,k)                    &
                      -alpha_p*dpdt(i,j,k) ) + dldt(i,j,k)

! Calculate Qc

          qc = al * (q(i,j,k) + qcl(i,j,k) - qsl_tl(index_npt(i,j)))

! Calculate the relative rate of change of width of the distribution
! dbsdtbs from the forcing rate
          IF (dbsdtbs0 /= 0.0 .or. dbsdtbs1 /= 0.0) THEN
            dbsdtbs = (dbsdtbs0 * timestep + dqcdt * dbsdtbs1) *        &
               EXP(-dbsdtbs_exp * qc / (al * qsl_tl(index_npt(i,j))))
          ELSE
            dbsdtbs = 0.0
          END IF
! ----------------------------------------------------------------------
! 4. Calculate the change of liquid cloud fraction. This uses the
! arrival value of QC for better behaved numerics.
! ----------------------------------------------------------------------

! DQCDT is the homogeneous forcing part, (QC+DQCDT)*DBSDTBS is the
! width narrowing part

          deltacl_c(index_npt(i,j)) = g * ( dqcdt - (qc + dqcdt) * dbsdtbs)

! Calculate the condensation amount DELTAL. This uses a mid value
! of cloud fraction for better numerical behaviour.

          c_1 = cfl(i,j,k) + deltacl_c(index_npt(i,j))
          IF (c_1  >   1.0) THEN
            deltacl_c(index_npt(i,j)) = 1.0 - cfl(i,j,k)
            c_1=1.0
          ELSE IF (c_1  <   0.0) THEN
            c_1=0.0
            deltacl_c(index_npt(i,j)) = (- cfl(i,j,k) )
          END IF
          c_1 = 0.5 * (c_1 + cfl(i,j,k))
          deltal = c_1 * dqcdt + (qcl(i,j,k) - qc * c_1) * dbsdtbs
!
! If we have removed all fraction, remove all liquid
          IF (l_fixbug_pc2_qcl_incr) THEN
             IF (cfl(i,j,k)+deltacl_c(index_npt(i,j)) == 0.0) THEN
                deltal = -qcl(i,j,k)
             END IF
          END IF
!
        ELSE IF ( ( cfl(i,j,k)  >=  (1.0 - cloud_rounding_tol) .AND.     &
                    cfl(i,j,k)  <=  (1.0 + cloud_rounding_tol) ) .OR.    &
! this if test is wrong, it should be cfl >= 1
! add fix on a switch to preserve bit-comparison
                  ( cfl(i,j,k)  >=  (1.0 - cloud_rounding_tol) .AND.     &
                  l_fixbug_pc2_mixph ) ) THEN

! Cloud fraction is 1

          alpha=repsilon*lc*qsl_t(index_npt(i,j))/(r*t_c(index_npt(i,j))**2)
          al=1.0/(1.0+lcrcp*alpha)
          alpha_p = -qsl_t(index_npt(i,j))/p_c(index_npt(i,j))
          deltal=al * (dqdt(i,j,k)-alpha*dtdt(i,j,k)                    &
                       -alpha_p*dpdt(i,j,k)) + dldt(i,j,k)

          deltacl_c(index_npt(i,j)) = 0.0

        ELSE

! Cloud fraction is 0 (actually CFL < cloud_rounding_tol)

          deltal = 0.0

        END IF

! Update water contents
        IF (l_fixbug_pc2_qcl_incr) THEN
! Don't allow more water to be removed than was there to start with.
          IF (qcl(i,j,k) + deltal < 0.0) THEN  
            deltal = -qcl(i,j,k)
! Set qcl to be exactly zero in this case.  
            qcl(i,j,k) = 0.0
          ELSE
            qcl(i,j,k) = qcl(i,j,k) + deltal
          ENDIF

        ELSE
          qcl(i,j,k) = qcl(i,j,k) + deltal
        ENDIF

! Update vapour content
! Q = input Q + Forcing - Condensation
        q(i,j,k)   = q(i,j,k)   + dqdt(i,j,k)                           &
                    - (deltal - dldt(i,j,k))

! Update temperature due to latent heating
        t(i,j,k)   = t(i,j,k)   + dtdt(i,j,k)                           &
                      + lcrcp * (deltal - dldt(i,j,k))

      END DO !i

    END DO !j

! ----------------------------------------------------------------------
! 5. Now update cloud fractions.
! ----------------------------------------------------------------------

! Calculate change in total cloud fraction.

    IF (npt > 0) THEN
! DEPENDS ON: pc2_total_cf
      CALL pc2_total_cf(                                                &
            npt,cfl_c,cff_c,deltacl_c,deltacf_c,cf_c)
    END IF

    DO j = qdims%j_start, qdims%j_end

      DO i = qdims%i_start, qdims%i_end

        IF (cfl(i,j,k) >        cloud_rounding_tol .AND.                &
            cfl(i,j,k) < (1.0 - cloud_rounding_tol)) THEN
! Update cloud fractions
          cf(i,j,k)  = cf_c(index_npt(i,j))
          cfl(i,j,k) = cfl(i,j,k) + deltacl_c(index_npt(i,j))
        END IF ! End if for CFL gt 0 and CFL lt 1

      END DO !i

    END DO !j

  END DO
!$OMP END PARALLEL DO

! End of the subroutine

  IF (lhook) CALL dr_hook('PC2_HOMOG_PLUS_TURB',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE pc2_homog_plus_turb
