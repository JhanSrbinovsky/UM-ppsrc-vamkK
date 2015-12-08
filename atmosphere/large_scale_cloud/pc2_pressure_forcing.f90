! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!+ PC2 Cloud Scheme: Forcing due to advection (adiabatic cooling)
! Subroutine Interface:
SUBROUTINE pc2_pressure_forcing(                                        &
!      Pressure related fields
 p, pstar, p_theta_levels,                                              &
!      Array dimensions
 rhc_row_length, rhc_rows,                                              &
!      timestep and rhcrit values
 timestep, rhcpt,                                                       &
!      Prognostic Fields
 theta, cf, cfl, cff, q, qcl, qcf,                                      &
!      Forcing quantities for driving the homogeneous forcing
 exner_departure, exner_current,                                        &
!      Convective cloud base and convection flag
 ccb, cumulus, rhts, tlts, qtts, ptts, cf_area,                         &
!      Diagnostics
  t_inc, q_inc, qcl_inc, qcf_inc, cf_inc, cfl_inc, cff_inc,             &
  t_dini,q_dini,qcl_dini,qcf_dini,cf_dini,cfl_dini,cff_dini,            &
!      Switches
  l_mixing_ratio, l_acf_cusack, l_cld_area                              &
 )

  USE atmos_constants_mod,   ONLY: recip_kappa, pref
  USE yomhook,               ONLY: lhook, dr_hook
  USE parkind1,              ONLY: jprb, jpim
  USE atm_fields_bounds_mod, ONLY: qdims, pdims,tdims, pdims_s, tdims_s,&
                                   qdims_l

  IMPLICIT NONE

! Purpose:
!   This subroutine interfaces to the homogeneous forcing routine and
!   calculates the temperature forcing associated with adiabatic
!   changes in pressure across a timestep. This is mainly due to
!   large-scale vertical advection. It also checks values of water
!   contents and cloud fractions to check that they are sensible.
!
! Method:
!   Uses the departure and end values of exner to calculate the
!   pressure and, from theta, the temperature and temperature
!   forcing. Calls the homogeneous forcing routine to update values
!   of temperature, moisture and cloud, and converts the temperature
!   change to a theta increment.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Large Scale Cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation

!  Global Variables:----------------------------------------------------

!  Subroutine Arguments:------------------------------------------------
  INTEGER ::                                                            &
                        !, INTENT(IN)
   rhc_row_length, rhc_rows,                                            &
!       Row length and number of rows in the rhcpt variable
   ccb(            qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end)
!       Convective cloud base

  LOGICAL ::                                                            &
                        !, INTENT(IN)
  cumulus(         qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end),                          &
  l_cld_area,                                                           &
!       Switch for area cloud fraction param
  l_acf_cusack,                                                         &
!       ... to select Cusack
  l_mixing_ratio
!       Use mixing ratio formulation

  REAL ::                                                               &
                        !, INTENT(IN)
   p_theta_levels( pdims_s%i_start:pdims_s%i_end,                       & 
                   pdims_s%j_start:pdims_s%j_end,                       &  
                   pdims_s%k_start:pdims_s%k_end),                      &
!       Pressure at all points (Pa)
   p(              pdims_s%i_start:pdims_s%i_end,                       & 
                   pdims_s%j_start:pdims_s%j_end,                       &  
                   pdims_s%k_start:pdims_s%k_end+1),                    &
   pstar(          pdims%i_start:pdims%i_end,                           & 
                   pdims%j_start:pdims%j_end),                          &
   timestep,                                                            &
!       Model timestep (s)
   exner_departure(pdims_s%i_start:pdims_s%i_end,                       & 
                   pdims_s%j_start:pdims_s%j_end,                       &  
                   pdims_s%k_start:pdims_s%k_end),                      &
!       Departure value of exner (before the advection)
   exner_current(  pdims_s%i_start:pdims_s%i_end,                       & 
                   pdims_s%j_start:pdims_s%j_end,                       &  
                   pdims_s%k_start:pdims_s%k_end),                      &
!       Current value of exner
   rhcpt(rhc_row_length,rhc_rows,qdims%k_end),                          &
!       Values of critical relative humidity
   rhts(           qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end),                          &
!       Relative total humidity at start of timestep (time level n)
   tlts(           qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end),                          &
!       TL at start of timestep
   qtts(           qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end),                          &
!       qT at start of timestep
   ptts(           pdims%i_start:pdims%i_end,                           & 
                   pdims%j_start:pdims%j_end,                           &  
                   pdims%k_start:pdims%k_end)
!       Pressure at theta levels at start of timestep

  REAL ::                                                               &
                        !, INTENT(INOUT)
   theta(          tdims_s%i_start:tdims_s%i_end,                       & 
                   tdims_s%j_start:tdims_s%j_end,                       &  
                   tdims_s%k_start:tdims_s%k_end),                      &
!       Potential temperature (K)
   cf(             qdims_l%i_start:qdims_l%i_end,                       & 
                   qdims_l%j_start:qdims_l%j_end,                       &  
                   qdims_l%k_start:qdims_l%k_end),                      &
!       Total cloud fraction (no units)
   cfl(            qdims_l%i_start:qdims_l%i_end,                       & 
                   qdims_l%j_start:qdims_l%j_end,                       &  
                   qdims_l%k_start:qdims_l%k_end),                      &
!       Liquid cloud fraction (no units)
   cff(            qdims_l%i_start:qdims_l%i_end,                       & 
                   qdims_l%j_start:qdims_l%j_end,                       &  
                   qdims_l%k_start:qdims_l%k_end),                      &
!       Ice cloud fraction (no units)
   q(              qdims_l%i_start:qdims_l%i_end,                       & 
                   qdims_l%j_start:qdims_l%j_end,                       &  
                   qdims_l%k_start:qdims_l%k_end),                      &
!       Vapour content (kg water per kg air)
   qcl(            qdims_l%i_start:qdims_l%i_end,                       & 
                   qdims_l%j_start:qdims_l%j_end,                       &  
                   qdims_l%k_start:qdims_l%k_end),                      &
!       Liquid content (kg water per kg air)
   qcf(            qdims_l%i_start:qdims_l%i_end,                       & 
                   qdims_l%j_start:qdims_l%j_end,                       &  
                   qdims_l%k_start:qdims_l%k_end),                      &
!       Ice content (kg water per kg air)
   cf_area(          qdims%i_start:qdims%i_end,                         & 
                     qdims%j_start:qdims%j_end,                         &  
                                 1:qdims%k_end)
!       Area cloud fraction

  REAL ::                                                               &
                        !, INTENT(OUT)
   t_inc(          tdims%i_start:tdims%i_end,                           & 
                   tdims%j_start:tdims%j_end,                           &  
                               1:tdims%k_end),                          &
!       Change in temperature due to CONDENSATION (K)
   q_inc(          qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end),                          &
!       Change in vapour due to CONDENSATION (kg kg-1)
   qcl_inc(        qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end),                          &
!       Change in liquid due to CONDENSATION (kg kg-1)
   qcf_inc(        qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end),                          &
!       Change in ice due to CONDENSATION (kg kg-1)
   cf_inc(         qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end),                          &
!       Change in total cloud fraction due to CONDENSATION (no units)
   cfl_inc(        qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end),                          &
!       Change in liquid cloud fraction due to CONDENSATION (no units)
   cff_inc(        qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end)
!       Change in ice cloud fraction due to CONDENSATION (no units)

  REAL ::                                                               &
                        !, INTENT(OUT)
   t_dini(         tdims%i_start:tdims%i_end,                           & 
                   tdims%j_start:tdims%j_end,                           &
                               1:tdims%k_end),                          &
!       Change in temperature due to initiation CONDENSATION (K)
   q_dini(         qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &
                               1:qdims%k_end),                          &
!       Change in vapour due to initiation CONDENSATION (kg kg-1)
   qcl_dini(       qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &
                               1:qdims%k_end),                          &
!       Change in liquid due to initiation CONDENSATION (kg kg-1)
   qcf_dini(       qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &
                               1:qdims%k_end),                          &
!       Change in ice due to initiation CONDENSATION (kg kg-1)
   cf_dini(        qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &
                               1:qdims%k_end),                          &
!       Change in total cloud fraction due to initiation CONDENSATION
   cfl_dini(       qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &
                               1:qdims%k_end),                          &
!       Change in liquid cloud fraction due to initiation CONDENSATION
   cff_dini(       qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &
                               1:qdims%k_end)
!       Change in ice cloud fraction due to initiation CONDENSATION

!  External functions:


!  Local scalars--------------------------------------------------------

  REAL ::                                                               &
   t_departure,    & ! Temperature at the departure point (K)
   p_departure       ! Pressure at the departure point (Pa)

  INTEGER :: k,i,j   ! Loop counters: K   - vertical level index
!                                     I,J - horizontal position index

!  Local dynamic arrays-------------------------------------------------
  REAL ::                                                               &
   t(              tdims%i_start:tdims%i_end,                           &
                   tdims%j_start:tdims%j_end,                           &  
                               1:tdims%k_end),                          &
!       Temperature (K)
   deltat(         tdims%i_start:tdims%i_end,                           & 
                   tdims%j_start:tdims%j_end,                           &  
                               1:tdims%k_end),                          &
!       Change in temperature between departure and arrival points (K)
   deltap(         pdims%i_start:pdims%i_end,                           & 
                   pdims%j_start:pdims%j_end,                           &  
                   pdims%k_start:pdims%k_end),                          &
!       Change in pressure between departure and arrival points (Pa)
   zeros(          tdims%i_start:tdims%i_end,                           & 
                   tdims%j_start:tdims%j_end,                           &  
                               1:tdims%k_end),                          &
!       Array of zero values

!  Local arrays which do not contain halos
   p_theta_levels_no_halos(                                             &
                   pdims%i_start:pdims%i_end,                           & 
                   pdims%j_start:pdims%j_end,                           &  
                   pdims%k_start:pdims%k_end),                          &
   cf_no_halos(    qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end),                          &
   cfl_no_halos(   qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end),                          &
   cff_no_halos(   qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end),                          &
   q_no_halos(     qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end),                          &
   qcl_no_halos(   qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end),                          &
   qcf_no_halos(   qdims%i_start:qdims%i_end,                           & 
                   qdims%j_start:qdims%j_end,                           &  
                               1:qdims%k_end)

!  Local dummy variables for SCM Diagnostics to be passed to
!  PC2_initiation_ctl
  INTEGER, PARAMETER :: nSCMDpkgs = 12

  LOGICAL ::                                                            &
   L_SCMDiags(nSCMDpkgs)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL   (KIND=jprb)            :: zhook_handle

!- End of Header

! ==Main Block==--------------------------------------------------------

  IF (lhook) CALL dr_hook('PC2_PRESSURE_FORCING',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! 1. Calculate temperature at departure and current (arrival) points
! ----------------------------------------------------------------------

!$OMP  PARALLEL DO DEFAULT(SHARED) SCHEDULE(STATIC) PRIVATE(i, j, k,     &
!$OMP& t_departure, p_departure) 

  DO k = 1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end       
      DO i = qdims%i_start, qdims%i_end

        t_departure   = theta(i,j,k) * exner_departure(i,j,k)
        t(i,j,k)      = theta(i,j,k) * exner_current(i,j,k)
        deltat(i,j,k) = t(i,j,k)     - t_departure
        p_departure   = exner_departure(i,j,k) ** recip_kappa * pref
        deltap(i,j,k) = p_theta_levels(i,j,k) - p_departure
        zeros(i,j,k)  = 0.0

! Copy primary variables into non-halo variables
        p_theta_levels_no_halos(i,j,k) = p_theta_levels(i,j,k)
        cf_no_halos            (i,j,k) = cf(i,j,k)
        cfl_no_halos           (i,j,k) = cfl(i,j,k)
        cff_no_halos           (i,j,k) = cff(i,j,k)
        q_no_halos             (i,j,k) = q(i,j,k)
        qcl_no_halos           (i,j,k) = qcl(i,j,k)
        qcf_no_halos           (i,j,k) = qcf(i,j,k)

! Copy variables into increment variables for diagnostics
        t_inc  (i,j,k) = t(i,j,k)
        q_inc  (i,j,k) = q(i,j,k)
        qcl_inc(i,j,k) = qcl(i,j,k)
        qcf_inc(i,j,k) = qcf(i,j,k)
        cf_inc (i,j,k) = cf(i,j,k)
        cfl_inc(i,j,k) = cfl(i,j,k)
        cff_inc(i,j,k) = cff(i,j,k)

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! 2. Call homogeneous forcing routine with just the temperature forcing
! ----------------------------------------------------------------------

! DEPENDS ON: pc2_homog_plus_turb
  CALL pc2_homog_plus_turb(p_theta_levels_no_halos,                     &
                           qdims%k_end,                                 &
                           timestep, t,                                 &
                           cf_no_halos, cfl_no_halos, cff_no_halos,     &
                           q_no_halos, qcl_no_halos,                    &
                           deltat, zeros, zeros, deltap, 0.0, 0.0,      &
                           l_mixing_ratio)

! ----------------------------------------------------------------------
! 3. Remember that the homogeneous forcing routine will include the
!    forcing increments in its changes. Deltat has already been included
!    in the model, so we need to subtract it again.
! ----------------------------------------------------------------------

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(i, j , k)     &
!$OMP& SHARED(qdims, t, deltat)
  DO k = 1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end
        t(i,j,k) = t(i,j,k) - deltat(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! 4. Call initiation routine
! ----------------------------------------------------------------------

! Set dummy variables for SCM Diagnostics
  L_SCMDiags(1:nSCMDpkgs) = .FALSE.

! NB if you are changing the argument list to PC2_initiation_ctl
! please do an equivalent change in routine scm_main to keep the
! Single Column Model consistent.

! DEPENDS ON: pc2_initiation_ctl
  CALL pc2_initiation_ctl(                                              &
    rhc_row_length, rhc_rows,                                           &
    .FALSE., l_mixing_ratio,                                            &
    l_acf_cusack, l_cld_area,                                           &
    timestep,                                                           &
    nSCMDpkgs, L_SCMDiags,                                              &
    t, q_no_halos, qcl_no_halos, qcf_no_halos,                          &
    cf_no_halos,cfl_no_halos,cff_no_halos,rhts,                         &
    tlts, qtts, ptts, cf_area, p, pstar,                                &
    p_theta_levels_no_halos, ccb, cumulus,                              &
    rhcpt,t_dini,q_dini,qcl_dini,qcf_dini,cf_dini,cfl_dini,             &
    cff_dini)

! ----------------------------------------------------------------------
! 5. Update variables
! ----------------------------------------------------------------------

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k) SCHEDULE(STATIC)
  DO k = 1, qdims%k_end
    DO j = qdims%j_start, qdims%j_end
      DO i = qdims%i_start, qdims%i_end

! Convert change of temperature to change of theta
        theta(i,j,k) = t(i,j,k) / exner_current(i,j,k)

! Copy changed non-halo variables into primary variables
        cf (i,j,k) = cf_no_halos(i,j,k)
        cfl(i,j,k) = cfl_no_halos(i,j,k)
        cff(i,j,k) = cff_no_halos(i,j,k)
        q  (i,j,k) = q_no_halos(i,j,k)
        qcl(i,j,k) = qcl_no_halos(i,j,k)
        qcf(i,j,k) = qcf_no_halos(i,j,k)

! Increment in temperature due to the CONDENSATION associated with the
! change in pressure. This is not the same as the temperature change
! due to vertical advection, which is in delta_t
        t_inc  (i,j,k) = t(i,j,k)   - t_inc(i,j,k)  -  t_dini(i,j,k)

! Other increment variables
        q_inc  (i,j,k) = q(i,j,k)   - q_inc(i,j,k)  -  q_dini(i,j,k)
        qcl_inc(i,j,k) = qcl(i,j,k) - qcl_inc(i,j,k)-qcl_dini(i,j,k)
        qcf_inc(i,j,k) = qcf(i,j,k) - qcf_inc(i,j,k)-qcf_dini(i,j,k)
        cf_inc (i,j,k) = cf(i,j,k)  - cf_inc(i,j,k) - cf_dini(i,j,k)
        cfl_inc(i,j,k) = cfl(i,j,k) - cfl_inc(i,j,k)-cfl_dini(i,j,k)
        cff_inc(i,j,k) = cff(i,j,k) - cff_inc(i,j,k)-cff_dini(i,j,k)

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO


! End of the subroutine

  IF (lhook) CALL dr_hook('PC2_PRESSURE_FORCING',zhook_out,zhook_handle)
  RETURN
END SUBROUTINE pc2_pressure_forcing
