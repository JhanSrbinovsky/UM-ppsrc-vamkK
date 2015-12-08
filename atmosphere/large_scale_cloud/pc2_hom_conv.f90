! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ Cloud Scheme: Homogenous forcing and Turbulence (non-updating)
! Subroutine Interface:
SUBROUTINE pc2_hom_conv(                                                &
!      Pressure related fields
 p_theta_levels,                                                        &
!      Array dimensions
 nlevels,                                                               &
!      Timestep
 timestep,                                                              &
!      Prognostic Fields
 t, q, qcl, cf, cfl, cff,                                               &
!      Forcing quantities for driving the homogenous forcing
 dtin, dqin, dqclin, dpdt, dcflin,                                      &
!      Cloud fraction in layers above and below.  
 cf_above, cf_below,                                                    &
!      Output increments to the prognostic fields
 dtpc2, dqpc2, dqclpc2, dcfpc2, dcflpc2,                                &
!      Other quantities for the turbulence
 pc2mixingrate, dbsdtbs1, l_mixing_ratio)

  USE water_constants_mod,   ONLY: lc
  USE atmos_constants_mod,   ONLY: cp, r, repsilon
  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim
  USE atm_fields_bounds_mod, ONLY: pdims, qdims, tdims
  USE cloud_inputs_mod,      ONLY: i_pc2_erosion_method,                &
       l_fixbug_pc2_qcl_incr,l_fixbug_pc2_mixph
  USE pc2_constants_mod,     ONLY: pc2eros_exp_rh,                      &
       pc2eros_hybrid_allfaces, pc2eros_hybrid_sidesonly,               &
       pdf_power, pdf_merge_power, dbsdtbs_exp, cloud_rounding_tol
  USE ereport_mod, ONLY : ereport

  IMPLICIT NONE

! Description:
!   This subroutine calculates the change in liquid content, liquid
!   cloud fraction and total cloud fraction as a result of homogenous
!   forcing of the gridbox with temperature, pressue, vapour and liquid
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
! Documentation: PC2 cloud scheme documentation

! Subroutine Arguments:------------------------------------------------

! Arguments with intent in. ie: input variables.

  INTEGER ::                                                            &
                        !, INTENT(IN)
   nlevels
!       No. of levels being processed.

  REAL ::                                                               &
                        !, INTENT(IN)
   timestep,                                                            &
!       Model timestep (s)
   pc2mixingrate,                                                            &
!       If i_pc2_erosion_method==pc2eros_exp_rh
!       Value of dbs/dt / bs which is independent of forcing (no units)
!       If i_pc2_erosion_method==pc2eros_hybrid_allfaces OR
!       i_pc2_erosion_method==pc2eros_hybrid_sidesonly:
!       Erosion rate to use in hybrid erosion method.
   dbsdtbs1
!       Used if i_pc2_erosion_method==pc2eros_exp_rh
!       Value of dbs/dt / bs which is proportional to the forcing
!       ( (kg kg-1 s-1)-1 )

  REAL ::                                                               &
                        !, INTENT(IN)
   p_theta_levels(pdims%i_start:pdims%i_end,                            & 
                  pdims%j_start:pdims%j_end,                            &  
                  nlevels),                                             &
!       Pressure at all points (Pa)
   t(             tdims%i_start:tdims%i_end,                            & 
                  tdims%j_start:tdims%j_end,                            &  
                  nlevels),                                             &
!       Temperature (K)
   q(             qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!       Vapour content (kg water per kg air)
   qcl(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!       Liquid content (kg water per kg air)
   cf(            qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!       Total cloud fraction (no units)
   cfl(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!       Liquid cloud fraction (no units)
   cff(           qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!       Ice cloud fraction (no units)
   dtin(          tdims%i_start:tdims%i_end,                            & 
                  tdims%j_start:tdims%j_end,                            &  
                  nlevels),                                             &
!       Increment of temperature from forcing mechanism (K)
   dqin(          qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!       Increment of vapour from forcing mechanism (kg kg-1)
   dqclin(        qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!       Increment of liquid from forcing mechanism (kg kg-1)
   dpdt(          pdims%i_start:pdims%i_end,                            & 
                  pdims%j_start:pdims%j_end,                            &  
                  nlevels),                                             &
!       Increment in pressure from forcing mechanism (Pa)
   dcflin(        qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!       Increment in liquid cloud fraction (no units)
   cf_above(      qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!       Total cloud fraction in layer above (no units)
   cf_below(      qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels)
!       Total cloud fraction in layer below (no units)

  LOGICAL ::                                                            &
                         !, INTENT(IN)
   l_mixing_ratio        ! Use mixing ratio formulation

! Arguments with intent out. ie: output variables.

  REAL ::                                                               &
                        !, INTENT(OUT)
   dtpc2(         tdims%i_start:tdims%i_end,                            & 
                  tdims%j_start:tdims%j_end,                            &  
                  nlevels),                                             &
!       PC2 Increment to Temperature (K)
   dqpc2(         qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!       PC2 Increment to Vapour content (kg water per kg air)
   dqclpc2(       qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!       PC2 Increment to Liquid content (kg water per kg air)
   dcfpc2(        qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels),                                             &
!       PC2 Increment to Total cloud fraction (no units)
   dcflpc2(       qdims%i_start:qdims%i_end,                            & 
                  qdims%j_start:qdims%j_end,                            &  
                  nlevels)
!       PC2 Increment to Liquid cloud fraction (no units)

!  External subroutine calls: ------------------------------------------

!  Local parameters and other physical constants------------------------
  REAL, PARAMETER :: lcrcp = lc / cp
!       Latent heat of condensation divided by heat capacity of air.
  REAL, PARAMETER :: b_factor = (pdf_power+1.0) / (pdf_power+2.0)
!       Premultiplier to calculate the amplitude of the probability
!       density function at the saturation boundary (G_MQC).
  REAL, PARAMETER :: smallp = 1.0e-10
!       Small positive value for use in if tests

!  Local scalars--------------------------------------------------------

  REAL ::                                                               &
   alpha,                                                               &
!       Rate of change of saturation specific humidity with
!       temperature calculated at dry-bulb temperature (kg kg-1 K-1)
   alpha_p,                                                             &
!       Rate of change of saturation specific humidity with
!       pressure calculated at dry-bulb temperature (Pa K-1)
   al,                                                                  &
!       1 / (1 + alpha L/cp)  (no units)
   c_1,                                                                 &
!       Mid-timestep liquid cloud fraction (no units)
   dbsdtbs,                                                             &
!       Relative rate of change of distribution width (s-1)
   dqcdt,                                                               &
!       Forcing of QC (kg kg-1 s-1)
   deltal,                                                              &
!       Change in liquid content (kg kg-1)
   cfl_to_m,                                                            &
!       CFL(i,j,k)**PDF_MERGE_POWER
   sky_to_m,                                                            &
!       (1-CFL(i,j,k))**PDF_MERGE_POWER
   g_mqc,                                                               &
!       Amplitude of the probability density function at
!       the saturation boundary (kg kg-1)-1
   qc,                                                                  &
!       aL (q + l - qsat(TL) )  (kg kg-1)
   sd,                                                                  &
!       Saturation deficit (= aL (q - qsat(T)) )  (kg kg-1)
   dcs        
!       Injected cloud fraction

!  (b)  Others.
  INTEGER :: k,i,j ! Loop counters:   K - vertical level index
!                                   I,J - horizontal position index
  INTEGER :: npt

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL   (KIND=jprb)            :: zhook_handle

!  Local arrays---------------------------------------------------------

  REAL ::                                                               &
   qsl_t(     (1+qdims%i_end-qdims%i_start)*                            &
              (1+qdims%j_end-qdims%j_start) ),                          &
!       Saturated specific humidity for dry bulb temperature T
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

  REAL :: tmp
  REAL :: dqcl, dcl, midpoint_qcl
  REAL :: exposed_area, side_face, satdiff, top_face, bot_face

  CHARACTER(LEN=256)       :: message
  INTEGER                  :: errorstatus

!- End of Header

! ==Main Block==--------------------------------------------------------

  IF (lhook) CALL dr_hook('PC2_HOM_CONV',zhook_in,zhook_handle)

! Loop round levels to be processed
! Levels_do1:
  DO k = 1, nlevels

! copy points into compressed arrays
    npt = 0
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        IF (cfl(i,j,k) > cloud_rounding_tol ) THEN
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
          alpha   = repsilon*lc*qsl_t(index_npt(i,j)) /                 &
                    (r*t_c(index_npt(i,j))**2)
          al      = 1.0 / ( 1.0 + lcrcp * alpha )
          alpha_p = -qsl_t(index_npt(i,j)) / p_c(index_npt(i,j))

! Calculate the saturation deficit SD

          sd      = al * ( qsl_t(index_npt(i,j)) - q(i,j,k) )

! Calculate the amplitude of the probability density function at the
! saturation boundary.

          IF (qcl(i,j,k) > smallp .AND. sd > smallp) THEN

            cfl_to_m = cfl(i,j,k)**pdf_merge_power
            sky_to_m = (1.0 - cfl(i,j,k))**pdf_merge_power

            g_mqc = b_factor * ( (1.0-cfl(i,j,k))**2 *                  &
             cfl_to_m / (sd * (cfl_to_m + sky_to_m))                    &
                  +                   cfl(i,j,k)**2 *                   &
             sky_to_m / (qcl(i,j,k) * (cfl_to_m + sky_to_m)) )

          ELSE
            g_mqc = 0.0
          END IF

! Calculate the rate of change of Qc due to the forcing

          dqcdt = al * ( dqin(i,j,k) - alpha*dtin(i,j,k)                &
                 -alpha_p*dpdt(i,j,k) ) + dqclin(i,j,k)

! For the background homogeneous forcing from the convection there is
! an additional term because the detrained plume must be saturated.
! This can also be written as a forcing.

!              IF (CFL(i,j,k)  >   0.0 .AND. CFL(i,j,k)  <   1.0) THEN
! This if test is already guaranteed
          IF (dcflin(i,j,k) > 0.0) THEN
            dcs = dcflin(i,j,k) / (1.0 - cfl(i,j,k))
          ELSE IF (dcflin(i,j,k) < 0.0) THEN
            dcs = dcflin(i,j,k) / ( - cfl(i,j,k))
          ELSE
            dcs = 0.0
          END IF
! Limit DCS to 0 and 1
          dcs   = MAX( MIN(dcs,1.0) ,0.0)

          dqcdt = dqcdt - al * dcs * (qsl_t(index_npt(i,j))-q(i,j,k))

! Calculate Qc

          qc    = al * ( q(i,j,k) + qcl(i,j,k) - qsl_tl(index_npt(i,j)) )

! Calculate the relative rate of change of width of the distribution
! dbsdtbs from the forcing rate

          IF (i_pc2_erosion_method == pc2eros_exp_rh) THEN
            ! Original Wilson et al (2008) formulation: rate of  
            ! narrowing related to RH via an ad-hoc exponetial. 
            dbsdtbs = (pc2mixingrate * timestep + dqcdt * dbsdtbs1) *   &
                   EXP(-dbsdtbs_exp * qc / (al * qsl_tl(index_npt(i,j))))
          ELSE IF (i_pc2_erosion_method == pc2eros_hybrid_allfaces .OR. &
                   i_pc2_erosion_method == pc2eros_hybrid_sidesonly) THEN
            ! Hybrid method 
            dbsdtbs = 0.0
            ! By setting this to zero, the next bit of code will not 
            ! do any width-narrowing to represent erosion. 
            ! So need to represent erosion in some other way. 
          ELSE
            errorstatus=10
            message='Attempting to use undefined i_pc2_erosion_method'
            CALL ereport('PC2_hom_conv',errorstatus,message)
          END IF 

! ----------------------------------------------------------------------
! 4. Calculate the change of liquid cloud fraction. This uses the
! arrival value of QC for better behaved numerics.
! ----------------------------------------------------------------------

! DQCDT is the homogeneous forcing part, (QC+DQCDT)*DBSDTBS is the
! width narrowing part

          dcflpc2(i,j,k) = g_mqc * ( dqcdt - (qc + dqcdt)*dbsdtbs)

! Calculate the condensation amount DELTAL. This uses a mid value
! of cloud fraction for better numerical behaviour.

          c_1 = MAX( 0., MIN( (cfl(i,j,k) + dcflpc2(i,j,k)), 1.) )

          dcflpc2(i,j,k) = c_1 - cfl(i,j,k)

          c_1 = 0.5 * (c_1 + cfl(i,j,k))

          IF (l_fixbug_pc2_qcl_incr) THEN
            ! Calculate increment here, without checking it for
            ! potential removal of more QCL than there is.
            ! A check is carried out at end of routine.
            deltal = (c_1 * dqcdt) +                                    &
                   ! The homogeneous forcing part 
                     ( ( qcl(i,j,k) - (qc * c_1) ) * dbsdtbs)
                   ! The PDF width-narrowing part.

            ! If we have removed all fraction, remove all liquid
            IF(cfl(i,j,k)+dcflpc2(i,j,k) == 0.0) THEN
               deltal = -qcl(i,j,k)
            END IF
          ELSE
            ! Original code. Note that only the width-narrowing 
            ! (erosion) bit is limited to not remove more QCL than there
            ! is. The homog forcing part can remove too much QCL.
            deltal = c_1 * dqcdt +                                      &
                MAX( (qcl(i,j,k) - qc * c_1) * dbsdtbs , (-qcl(i,j,k)) )
          ENDIF

          IF (i_pc2_erosion_method == pc2eros_exp_rh) THEN
! Only calculate change in total cloud fraction here if not using 
! hybrid method as need to add increments from hybrid erosion. Change in
! total cloud fraction will be done later.

! ----------------------------------------------------------------------
! 5. Calculate change in total cloud fraction.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! The following If test is a copy of the PC2_TOTAL_CF subroutine.
! ----------------------------------------------------------------------
            IF (dcflpc2(i,j,k) > 0.0) THEN
! ...  .AND. CFL(i,j,k)  <   1.0 already assured.
               IF (l_fixbug_pc2_mixph) THEN
! minimum overlap, this is consistent with pc2_totalcf
                  dcfpc2(i,j,k) = MIN(dcflpc2(i,j,k),(1.0-cf(i,j,k)))
               ELSE
! random overlap, this is inconsistent with pc2_totalcf
                  dcfpc2(i,j,k) = dcflpc2(i,j,k) * (1.0 - cf(i,j,k)) /  &
                                                   (1.0 - cfl(i,j,k))
               END IF
            ELSE IF (dcflpc2(i,j,k) < 0.0) THEN
! ...  .AND. CFL(i,j,k)  >   0.0 already assured.
               IF (l_fixbug_pc2_mixph) THEN
! minimum overlap, this is consistent with pc2_totalcf
                  dcfpc2(i,j,k) = MAX(dcflpc2(i,j,k),(cff(i,j,k)-cf(i,j,k)))
               ELSE
! random overlap, this is inconsistent with pc2_totalcf
                  dcfpc2(i,j,k) = dcflpc2(i,j,k) *                      &
                      (cf(i,j,k)- cff(i,j,k)) / cfl(i,j,k)
               END IF
            ELSE
              dcfpc2(i,j,k) = 0.0
            END IF
! ----------------------------------------------------------------------
          END IF

        ELSE IF ( ( cfl(i,j,k)  >=  (1.0 - cloud_rounding_tol) .AND.     &
                    cfl(i,j,k)  <=  (1.0 + cloud_rounding_tol) ) .OR.    &
! this if test is wrong, it should be cfl >= 1
! add fix on a switch to preserve bit-comparison
                  ( cfl(i,j,k)  >=  (1.0 - cloud_rounding_tol) .AND.     &
                  l_fixbug_pc2_mixph ) ) THEN

! Cloud fraction is 1

          dcfpc2(i,j,k)  = 0.0
          dcflpc2(i,j,k) = 0.0

          alpha   = repsilon * lc * qsl_t(index_npt(i,j)) /             &
                   (r * t_c(index_npt(i,j))**2)
          al      = 1.0 / (1.0 + lcrcp*alpha)
          alpha_p = -qsl_t(index_npt(i,j)) / p_c(index_npt(i,j))
          deltal  = al * (dqin(i,j,k) - alpha*dtin(i,j,k)               &
                -alpha_p*dpdt(i,j,k)) + dqclin(i,j,k)

        ELSE

! Cloud fraction is 0

          dcfpc2(i,j,k)  = 0.0
          dcflpc2(i,j,k) = 0.0

          deltal         = 0.0

        END IF

!===========================================================
! Hybrid PC2 erosion
!===========================================================
        IF (i_pc2_erosion_method == pc2eros_hybrid_allfaces .OR.        &
            i_pc2_erosion_method == pc2eros_hybrid_sidesonly) THEN
          ! Although this alternative method of doing erosion
          ! is being done in a separate, subsequent bit of code to the
          ! homogeneous forcing, it is effectively making its 
          ! calculations in parallel since it uses the input values of
          ! cloud fields and thermodynamics.

          IF (i_pc2_erosion_method == pc2eros_hybrid_allfaces) THEN
            ! Assuming maximum overlap between layers where cloud is  
            ! present in both layers
            top_face=max(0.0, cf(i,j,k)-cf_above(i,j,k))
            bot_face=max(0.0, cf(i,j,k)-cf_below(i,j,k))
          END IF

          IF (cfl(i,j,k) > cloud_rounding_tol .AND.                     &
              cfl(i,j,k) < (1.0-cloud_rounding_tol)) THEN
            ! Erosion of cloud is assumed to only happen from the  
            ! cloud surface area exposed to clear sky. 
 
            ! Calculate the difference from saturation
            satdiff = qsl_t(index_npt(i,j)) - q(i,j,k)

            ! Calculate exposed lateral surface area. Define a function 
            ! which is an upside-down U shape, going to zero at CFL=0  
            ! and CFL=1 and with a peak value of 0.5 at CFL=0.5. 
            side_face = (2.0*cfl(i,j,k)) - (2.0*cfl(i,j,k)*cfl(i,j,k))

            IF (i_pc2_erosion_method == pc2eros_hybrid_allfaces) THEN
              ! Exposed surface area is made up of contribution
              ! from lateral sides and from top and bottom faces.
              exposed_area = side_face + top_face + bot_face
            ELSE !i_pc2_erosion_method == pc2eros_hybrid_sidesonly
              exposed_area = side_face
            END IF

            dqcl = - exposed_area * pc2mixingrate * satdiff * timestep

            ! Find the value of qcl half way through the erosion  
            ! process for better numerical behaviour. 
            midpoint_qcl = qcl(i,j,k) + ( 0.5 * dqcl ) 

            ! Assume that the change in QCL is due to 
            ! width-narrowing. Use that width-narrowing rate to find 
            ! the consistent change in cloud fraction. 
            tmp = (midpoint_qcl - (qc * cfl(i,j,k)) ) 
 
            dcl = - g_mqc * qc * dqcl / tmp 

          ELSE IF ( ( cfl(i,j,k)  >=  (1.0 - cloud_rounding_tol) .AND.  &
                    cfl(i,j,k)  <=  (1.0 + cloud_rounding_tol) ) .OR.   &
! this if test is wrong, it should be cfl >= 1
! add fix on a switch to preserve bit-comparison
                  ( cfl(i,j,k)  >=  (1.0 - cloud_rounding_tol) .AND.    &
                  l_fixbug_pc2_mixph ) ) THEN
            ! No contribution from lateral exposed area, 
            ! just from top and bottom.

            IF (i_pc2_erosion_method == pc2eros_hybrid_allfaces) THEN

              exposed_area = top_face + bot_face
     
              ! Calculate the difference from saturation
              satdiff = qsl_t(index_npt(i,j)) - q(i,j,k)

              dqcl = - exposed_area * pc2mixingrate * satdiff * timestep

              ! If the sky is overcast, and there is no lateral edges 
              ! where mixing can occur, then mixing can only happen at
              ! top and bottom, this is only likely to reduce the QCL 
              ! of the cloud and not reduce its lateral extent, so in 
              ! the case of CFL=1 assume 
              ! no change in cloud fraction due to erosion.
              dcl = 0.0

            ELSE !i_pc2_erosion_method == pc2eros_hybrid_sidesonly
              ! If we assume mixing is only happening from the sides
              ! and not the top and bottom, then when cfl=1 there
              ! will be no sides exposed and so no mixing at all.
              dqcl = 0.0
              dcl  = 0.0
            END IF

          ELSE
            ! Cloud fraction = 0, so no cloud there to remove. 
            ! Set sink terms to zero. 
            dcl=0.0
            dqcl=0.0
          END IF

          ! Add the the erosion increments to the homogeneous forcing ones
          dcflpc2(i,j,k) = dcflpc2(i,j,k) + dcl

          deltal = deltal + dqcl

          ! Test for silly things using original value and total increments

          IF (qcl(i,j,k)+deltal <= 0.0) THEN
            ! If we are about to removed all QCL remove all CFL.
            dcflpc2(i,j,k) = -cfl(i,j,k)
            deltal = -qcl(i,j,k)
          END IF

          IF ( cfl(i,j,k)+dcflpc2(i,j,k) <= 0.0 ) THEN
            ! If we are about to removed all CFL remove all QCL.
            dcflpc2(i,j,k) = -cfl(i,j,k)
            deltal = -qcl(i,j,k)
          END IF

          IF ( cfl(i,j,k)+dcflpc2(i,j,k) > 1.0 ) THEN
            ! If we are about to make CFL>1 then only increase it to 1.
            dcflpc2(i,j,k) = 1.0-cfl(i,j,k)
          END IF

          ! Calculate change in total cloud fraction.

! The following only needs to be done if 0<CFL<1, but is it worth enforcing that?
! as dcflpc2 will only be non-zero if thing have been done to it in 0<cfl<1 bits

! ----------------------------------------------------------------------
! The following If test is a copy of the PC2_TOTAL_CF subroutine.
! ----------------------------------------------------------------------
            IF (dcflpc2(i,j,k) > 0.0) THEN
! ...  .AND. CFL(i,j,k)  <   1.0 already assured.
               IF (l_fixbug_pc2_mixph) THEN
! minimum overlap, this is consistent with pc2_totalcf
                  dcfpc2(i,j,k) = MIN(dcflpc2(i,j,k),(1.0-cf(i,j,k)))
               ELSE
! random overlap, this is inconsistent with pc2_totalcf
                  dcfpc2(i,j,k) = dcflpc2(i,j,k) * (1.0 - cf(i,j,k)) /  &
                                                   (1.0 - cfl(i,j,k))
               END IF
            ELSE IF (dcflpc2(i,j,k) < 0.0) THEN
! ...  .AND. CFL(i,j,k)  >   0.0 already assured.
               IF (l_fixbug_pc2_mixph) THEN
! minimum overlap, this is consistent with pc2_totalcf
                  dcfpc2(i,j,k) = MAX(dcflpc2(i,j,k),(cff(i,j,k)-cf(i,j,k)))
               ELSE
! random overlap, this is inconsistent with pc2_totalcf
                  dcfpc2(i,j,k) = dcflpc2(i,j,k) *                      &
                    (cf(i,j,k)- cff(i,j,k)) / cfl(i,j,k)
               END IF
            ELSE
              dcfpc2(i,j,k) = 0.0
            END IF
! ----------------------------------------------------------------------

        END IF ! (i_pc2_erosion_method == hybrid_method)

!===========================================================

! Increment water contents and temperature due to latent heating.
! This subroutine will output only the condensation increments
! hence we comment out updates to qcl, q and t
!           QCL(i,j,k) = QCL(i,j,k) + DELTAL
! Q = input Q + Forcing - Condensation
!           Q(i,j,k)   = Q(i,j,k) + DQIN(i,j,k) - (DELTAL - DLIN(i,j,k))
!           T(i,j,k)   = T(i,j,k) + DTIN(i,j,k)
!    &                            + LCRCP * (DELTAL - DLIN(i,j,k))

! These are the condensation increments
        dqclpc2(i,j,k) = deltal - dqclin(i,j,k)

        IF (l_fixbug_pc2_qcl_incr) THEN
          ! Ensure QCL increment cannot make QCL go negative.
          IF ( qcl(i,j,k) + dqclpc2(i,j,k) < 0.0) THEN 
            dqclpc2(i,j,k) = -qcl(i,j,k) 
          ENDIF 
        ENDIF

        dqpc2(i,j,k)   = - dqclpc2(i,j,k)
        dtpc2(i,j,k)   = lcrcp * dqclpc2(i,j,k)

      END DO  ! i
    END DO  ! j
  END DO  ! k

! End of the subroutine

  IF (lhook) CALL dr_hook('PC2_HOM_CONV',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE pc2_hom_conv
