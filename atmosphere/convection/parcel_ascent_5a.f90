! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates parcel ascent
!
! Subroutine Interface:
SUBROUTINE parcel_ascent (npnts, wet_levels, nSCMDpkgs,                  &
                          nlcl, k_plume,                                 &
                          l_mixing_ratio, l_dilute, L_SCMDiags,          &
                          sl_plume, qw_plume,                            &
                          T, Tl, q, qcl, qcf, qw, t_dens_env,            &
                          p_theta_lev, exner_theta_levels,               &
                          z_theta, entrain_fraction,                     &
                          T_parc, t_parc_dil,                            &
                          buoyancy, buoyancy_dil,                        &
                          env_svl, par_svl, par_svl_dil,                 &
                          denv_bydz, dpar_bydz, dqsatdz )
                          
USE earth_constants_mod, ONLY: g

USE cv_run_mod, ONLY:                                                    &
  plume_water_load, dil_plume_water_load, qlmin, fac_qsat

USE cv_param_mod, ONLY:                                                  &
  qlcrit

USE cv_derived_constants_mod, ONLY:                                      &
  ls, lsrcp, lcrcp, gamma_dry

  USE atmos_constants_mod, ONLY:                                        &
      cp, r, repsilon, c_virtual

USE water_constants_mod, ONLY: lc, lf, tm
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   This routine calculates a parcel ascent from k_plume.
!   An undilute ascent is always calculated. If the option l_dilute is set 
!   to .true. a dilute ascent is also calculated. The dilute ascent uses
!   the entrainment rate held in entrain_fraction to mix in environmental air.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------

! Subroutine arguments
 
INTEGER, INTENT(IN) :: &
  npnts                & ! Number of points
 ,wet_levels           & ! Number of wet model levels 
 ,nSCMDpkgs              ! SCM  - No of diagnostics packages

INTEGER, INTENT(IN) :: &
  nlcl(npnts)          & ! Lifting condensation level
 ,k_plume(npnts)         ! Starting model level for plume ascent

LOGICAL, INTENT(IN) :: & 
  l_mixing_ratio       & ! .TRUE. if input q a mixing ratio otherwise 
                         ! assumes input q a specific humidity.
 ,l_dilute             & ! .TRUE. if a dilute parcel ascent also required.
 ,L_SCMDiags(nSCMDpkgs)  ! SCM - Logicals for diagnostics packages

REAL, INTENT(IN) ::           &
  T(npnts,wet_levels)   & ! Temperature on model levels (K)
 ,Tl(npnts,wet_levels)  & ! Liquid water temperature on model lev (K)
 ,q(npnts,wet_levels)   & ! water vapour on model levels (kg/kg)
 ,qcl(npnts,wet_levels) & ! cloud liquid water on model levels (kg/kg)
 ,qcf(npnts,wet_levels) & ! cloud ice water on model levels (kg/kg)
 ,qw(npnts,wet_levels)  & ! total water (kg/kg)
 ,t_dens_env(npnts,wet_levels)  
                             ! Density potential temperature of environment (K)

REAL, INTENT(IN) ::                          &
  p_theta_lev(npnts,wet_levels)        & ! Pressure on theta levels (Pa)
 ,exner_theta_levels(npnts,wet_levels) & ! Exner Pressure on theta levels 
 ,z_theta(npnts,wet_levels)            & ! Height of theta levels  (m)
 ,entrain_fraction(npnts,wet_levels)     ! fraction of environmental air 
                                               ! to mix with parcel
 
REAL, INTENT(INOUT) :: &
  sl_plume(npnts)      & ! SL at start of plume
 ,qw_plume(npnts)        ! total water at plume start (kg/kg)

REAL, INTENT(OUT) ::                   &
  t_parc(npnts,wet_levels)       & ! Parcel temperature  (K)
 ,t_parc_dil(npnts,wet_levels)   & ! Dilute Parcel temperature (K)
 ,buoyancy(npnts,wet_levels)     & ! Parcel buoyancy  (K)
 ,buoyancy_dil(npnts,wet_levels) & ! Dilute parcel buoyancy (K)
 ,env_svl(npnts,wet_levels)      & ! Density (virtual) static energy
                                         ! over CP for layer.
 ,par_svl(npnts,wet_levels)      & ! Density (virtual) static energy
                                         ! over CP of parcel for level.
 ,par_svl_dil(npnts,wet_levels)  & ! Density (virtual) static energy
                                         ! over CP of parcel for level.
 ,denv_bydz(npnts, wet_levels)   & ! Gradient of density potential
                                         ! temperature in the environment.
 ,dpar_bydz(npnts, wet_levels)   & ! Gradient of density potential
 ,dqsatdz(npnts, wet_levels)       ! dqsat/dz along an adiabat undilute 
                                         ! parcel

! Local variables

INTEGER ::        & 
 ii,k                ! loop counters

REAL ::           &
  q_liq_env       &  ! Condensed water content of environment.
 ,dq_sat_env      &  ! DQSAT/DT for environment
 ,lrcp_const      &  ! lc or lc+lf over cp
 ,lrcp_const_env  &  ! lc or lc+lf over cp
 ,lrcp_const_parc &  ! lc or lc+lf over cp
 ,l_const         &  ! lc or lc+lf
 ,l_const_env     &  ! lc or lc+lf
 ,dz              &  ! layer depth
 ,dtdz            &  ! temperature gradient along undilute parcel ascent
 ,z_pr            &  ! used in estimating th_ref at next level
 ,th_par          &  ! theta value for parcel
 ,dq_sat_par      &  ! dqsat/dT for parcel
 ,dq_sat_par_dil  &  ! dqsat/dT for parcel
 ,temp_parc       &  ! average temperature of parcel after entrainment
 ,q_vap_parc      &  ! Vapour content of undilute parcel
 ,q_liq_parc      &  ! Liquid water content of undilute parcel
 ,qcl_parc        &  ! parcel qcl - dilute parcel cal
 ,qcf_parc        &  ! parcel qcf - dilute parcel cal 
 ,ql_remove          ! ql removed from plume


! parcel calculation

REAL ::                                &
  t_ref(npnts)                         & ! reference temperature
 ,th_ref(npnts)                        & ! reference potential temperature
 ,th_par_km1(npnts)                    & ! theta refernce for level below
 ,qsat_lev(npnts)                      & ! qsat for reference temperature
 ,qsat_env(npnts)                      & ! qsat for environment temperature
 ,t_dens_parc(npnts, wet_levels)         ! Density potential temperature
                                         ! of parcel.
REAL ::                                &
  max_qw_plume                           ! maximum plume water content

! Arrays added for dilute parcel calculation
      
REAL ::                           &
  t_ref_dil(npnts)                & ! dilute parcel reference temperature 
 ,th_ref_dil(npnts)               & ! reference potential temperature
 ,th_par_km_dil(npnts)            & ! reference potential temperature 2nd
 ,qsat_lev_dil(npnts)             & ! qsat for dilute parcel
 ,qw_parc_dil(npnts,wet_levels)   & ! parcel total water dilute plume
 ,qw_parc(npnts,wet_levels)       & ! parcel total water undilute plume
 ,qv_parc(npnts,wet_levels)       & ! parcel water undilute plume (array not
                                    ! essential but may require in future)
 ,sl_parc(npnts,wet_levels)       & ! parcel SL undilute
 ,sl_parc_dil(npnts,wet_levels)   & ! parcel SL dilute
 ,ql_parc(npnts,wet_levels)       & ! parcel water
 ,t_dens_parc_dil(npnts,wet_levels) & ! dilute parcel t_dens
 ,ql_parc_dil(npnts,wet_levels)       ! dilute parcel liquid water

LOGICAL ::     &
  l_keep_water    ! if true keeps water loading in plume
                  ! false removed if water exceeds 1g/kg

! Model constants



INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('PARCEL_ASCENT',zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 1.0 Initialisation
!-----------------------------------------------------------------------

l_keep_water  = .false. ! water loading restricted to 1g/kg

! Initialise parcel reference theta and value for level below

DO ii=1, npnts
  th_ref(ii) = tl(ii,k_plume(ii))/exner_theta_levels(ii,k_plume(ii))
  th_par_km1(ii) = th_ref(ii)
END DO

IF (l_dilute) THEN     ! dilute parcel ascent 
  DO ii=1, npnts
    th_ref_dil(ii)    = th_ref(ii)
    th_par_km_dil(ii) = th_ref_dil(ii)
  END DO
END IF

!-----------------------------------------------------------------------
! 2.0 Parcel ascent 
!-----------------------------------------------------------------------
! Parcel starts from level k_plume and is lifted upwards to top of model.
! Dilute parcel ascent - mix in environmental air above lifting 
! condensation level
! 
!-----------------------------------------------------------------------
! Calculate parcel water by linearising qsat about the Parcel's
! temperature extrapolated up to the next grid_level.
!-----------------------------------------------------------------------

! Level loop calculating parcel ascent

DO  K = 1,wet_levels

  ! Require t_ref on all point for qsat call

!$OMP  PARALLEL DO DEFAULT(NONE) SHARED(t_ref, th_ref,                  &
!$OMP& exner_theta_levels, k, npnts) PRIVATE(ii)
  DO ii=1, npnts
    t_ref(ii)     = th_ref(ii)*exner_theta_levels(ii,k)
  END DO
!$OMP END PARALLEL DO

! DEPENDS ON: qsat_mix
  call qsat_mix(qsat_lev,t_ref,p_theta_lev(1,k),npnts,l_mixing_ratio)

! DEPENDS ON: qsat_mix
  call qsat_mix(qsat_env,t(1,k),p_theta_lev(1,k),npnts,l_mixing_ratio)

  IF (l_dilute) THEN     ! dilute parcel ascent 

!$OMP  PARALLEL DO DEFAULT(NONE) SHARED(t_ref_dil, th_ref_dil,            &
!$OMP& exner_theta_levels, k, npnts) PRIVATE(ii)
    DO ii=1, npnts
      t_ref_dil(ii) = th_ref_dil(ii)*exner_theta_levels(ii,k)
    END DO
!$OMP END PARALLEL DO

! DEPENDS ON: qsat_mix
    call qsat_mix(qsat_lev_dil,t_ref_dil,p_theta_lev(1,k)                 &
                                                ,npnts,l_mixing_ratio)
  END IF                 ! dilute parcel

! Undilute parcel calculation always required

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(ii, lrcp_const, lrcp_const_env,  &
!$OMP& l_const, dq_sat_env, dq_sat_par, q_liq_parc, q_liq_env, z_pr,     &
!$OMP& lrcp_const_parc, ql_remove, max_qw_plume, q_vap_parc, th_par,     &
!$OMP& temp_parc, qcl_parc, qcf_parc, l_const_env, dq_sat_par_dil)

!$OMP DO SCHEDULE(STATIC)
  DO ii=1, npnts

    If(T_ref(ii) >  TM) THEN
      lrcp_const = lcrcp
      l_const    = lc
    ELSE
      lrcp_const = lsrcp
      l_const    = ls
    END IF
    If(T(ii,k) >  TM) THEN
      lrcp_const_env = lcrcp
      l_const_env    = lc
    ELSE
      lrcp_const_env = lsrcp
      l_const_env    = ls
    END IF

    dq_sat_env = repsilon*l_const_env*qsat_env(ii)/(R*T(ii,k)**2)
    dq_sat_par = repsilon*l_const    *qsat_lev(ii)/(R*T_ref(ii)**2)

    q_liq_parc = MAX( 0.0, ( qw_plume(ii) - qsat_lev(ii)                     &
              -dq_sat_par*( sl_plume(ii)-gamma_dry*z_theta(ii,K)-T_ref(ii) ) &
                                 ) / (1.0+lrcp_const*dq_sat_par) )

    q_liq_env  = MAX( 0.0, ( qw(ii,K) - qsat_env(ii)                      &
                -dq_sat_env*( TL(ii,K)               - T(ii,k) )          &
                                 ) / (1.0+Lrcp_const_env*dq_sat_env) )
!
! add on the difference in the environment's ql as calculated by the
! UM cloud scheme (using some RH_CRIT value) and what it
! would be If RH_CRIT=1. This THEN imitates partial condensation
! in the parcel.
!
    ql_parc(ii,k) = q_liq_parc + qcl(ii,k) + qcf(ii,k) - q_liq_env 

    T_PARC(ii,k)=sl_plume(ii)-gamma_dry*z_theta(ii,K) +lrcp_const*ql_parc(ii,k)  


! May need to recalculate if T_parc is > Tm and T_ref < Tm

    IF (T_ref(ii) <= TM .and. T_parc(ii,k) >  TM) THEN

! recalculate using corrected latent heats
      lrcp_const_parc = lcrcp

      q_liq_parc = MAX( 0.0, ( qw_plume(ii) - qsat_lev(ii)                    &
           -dq_sat_par*( sl_plume(ii)-gamma_dry*z_theta(ii,K)-T_ref(ii) )     &
                          ) / (1.0+lrcp_const_parc*dq_sat_par) )

      ql_parc(ii,k) = q_liq_parc + qcl(ii,k) + qcf(ii,k)- q_liq_env 

! revised at parcel calculation

      T_PARC(ii,k)=sl_plume(ii)-gamma_dry*z_theta(ii,K)                       &
                                       +lrcp_const_parc*ql_parc(ii,k)

    END IF       ! test on T_ref and t_parc

    ! Add water removal from undilute plume 
    IF (plume_water_load == 1) THEN

     ! water removed from parcel after condensation if > 0.001 kg/kg
      IF (ql_parc(ii,k) > 0.001) THEN
        ql_remove = ql_parc(ii,k) -  0.001
        ql_parc(ii,k) = 0.001
        qw_plume(ii) = qw_plume(ii) - ql_remove
        ! Also adjust sl_plume as well as altered energy
        sl_plume(ii) = sl_plume(ii) + lrcp_const*ql_remove
      END IF
 
    ELSE IF (plume_water_load == 2) THEN

      ! water remaining depends on qsat environment 
      max_qw_plume = fac_qsat*qsat_env(ii)
      max_qw_plume = MIN (qlcrit,max_qw_plume)  ! max
      max_qw_plume = MAX (qlmin,max_qw_plume)   ! min value 
      IF (ql_parc(ii,k) > max_qw_plume) THEN
        ql_remove = ql_parc(ii,k) - max_qw_plume 
        ql_parc(ii,k) = max_qw_plume
        qw_plume(ii) = qw_plume(ii) - ql_remove
        ! Also adjust sl_plume as well as altered energy
        sl_plume(ii) = sl_plume(ii) + lrcp_const*ql_remove

      END IF           

    END IF    ! end test on plume_water_load

    q_vap_parc=qw_plume(ii)-ql_parc(ii,k)
    qv_parc(ii,k) = q_vap_parc 

    t_dens_parc(ii,k)=T_PARC(ii,k)*(1.0+c_virtual*q_vap_parc-ql_parc(ii,k))


! calculate t_ref for next level
    IF (k >  1 .and. k <   wet_levels-1) THEN
      z_pr = (z_theta(ii,k+1)-z_theta(ii,k))/(z_theta(ii,k)-z_theta(ii,K-1))
      th_par = t_parc(ii,k)/exner_theta_levels(ii,k)
      th_ref(ii) = th_par*(1.+z_pr) - th_par_km1(ii)*z_pr

! Check sensible value otherwise set to previous reference value
! Problems can occur near top of model where calculation are nolonger 
! important.
      IF (th_ref(ii) < 0.0) THEN
        th_ref(ii) = th_par_km1(ii)
      END IF
      IF (th_par > 0.0) THEN   
        th_par_km1(ii) = th_par
      END IF
    END IF

!-----------------------------------------------------------------------
! Dilute parcel ascent
!-----------------------------------------------------------------------
    IF (l_dilute) THEN 

      IF (k <= nlcl(ii)) THEN  ! Dilute parcel same as undilute
                               ! parcel ascent (no entrainment)

        sl_parc(ii,k) = sl_plume(ii)
        qw_parc(ii,k) = qw_plume(ii)
        t_parc_dil(ii,k)     = t_parc(ii,k)
        ql_parc_dil(ii,k)    = ql_parc(ii,k)         
        t_dens_parc_dil(ii,k)= t_dens_parc(ii,k)
        th_ref_dil(ii)       = th_ref(ii)
        th_par_km_dil(ii)    = th_par_km1(ii)

      ELSE                      ! Dilute parcel ascent now required

        If(t_ref_dil(ii) >  TM) THEN
          lrcp_const = lcrcp
          l_const    = lc
        ELSE
          lrcp_const = lsrcp
          l_const    = ls
        END IF

!-----------------------------------------------------------------------
! Dilute parcel
!-----------------------------------------------------------------------
! Mix in entrain_fraction from environmental air from level below and 
! raise this to current level.
! Assume mix in fraction of mass from environment.
! Estimate parcel properties after mixing air from environment with 
! parcel. Temperature given approximately by average

         temp_parc = (t_parc_dil(ii,k-1)                                   &
                                + entrain_fraction(ii,k)*t(ii,k-1))        &
                         /(1.+entrain_fraction(ii,k))

     
         qw_parc(ii,k) = (qw_parc(ii,k-1) +                                &
                               entrain_fraction(ii,k)*qw(ii,k-1))          & 
                          /(1.+entrain_fraction(ii,k)) 

         qcl_parc = (ql_parc_dil(ii,k-1)   +                               &    
                               entrain_fraction(ii,k)*qcl(ii,k-1))         & 
                          /(1.+entrain_fraction(ii,k)) 

         qcf_parc = (0.0     +                                             &    
                               entrain_fraction(ii,k)*qcf(ii,k-1))         & 
                          /(1.+entrain_fraction(ii,k)) 

! All condensed water either ice or liquid based on t_ref 
         sl_parc(ii,k) = temp_parc - lrcp_const*(qcl_parc+qcf_parc)        &
                                 +gamma_dry*z_theta(ii,k-1) 

         dq_sat_par_dil = repsilon*l_const*qsat_lev_dil(ii)                &
                                     /(R*t_ref_dil(ii)**2)

         q_liq_parc = MAX( 0.0, ( qw_parc(ii,k) - qsat_lev_dil(ii)            &
       -dq_sat_par_dil*( sl_parc(ii,k)-gamma_dry*z_theta(ii,K)-t_ref_dil(ii)) &
                                    ) / (1.0+lrcp_const*dq_sat_par_dil) )

! add on the dIfference in the environment's ql as calculated by the
! UM cloud scheme (using some RH_CRIT value) and what it
! would be If RH_CRIT=1. This THEN imitates partial condensation
! in the parcel.
!
         ql_parc_dil(ii,k) = q_liq_parc + qcl(ii,k) + qcf(ii,k) - q_liq_env
 
         t_parc_dil(ii,k) = sl_parc(ii,k)-gamma_dry*z_theta(ii,K)               &
                                       +lrcp_const*ql_parc_dil(ii,k)

! May need to recalculate if T_parc is > Tm and T_ref < Tm

         IF (t_ref_dil(ii) <= TM.and.t_parc_dil(ii,k) >  TM) THEN

! recalculate using corrected latent heats
           lrcp_const_parc = lcrcp

           q_liq_parc = MAX( 0.0, ( qw_parc(ii,k) - qsat_lev_dil(ii)          &
        -dq_sat_par_dil*(sl_parc(ii,k)-gamma_dry*z_theta(ii,K)-t_ref_dil(ii)) &
                          ) / (1.0+lrcp_const_parc*dq_sat_par_dil) )

           ql_parc_dil(ii,k) = q_liq_parc + qcl(ii,k)                      &
                                               + qcf(ii,k)- q_liq_env

! revised at parcel calculation

           t_parc_dil(ii,k)=sl_parc(ii,k)-gamma_dry*z_theta(ii,K)               &
                                    +lrcp_const_parc*ql_parc_dil(ii,k)

         END IF   ! test on t_ref

         q_vap_parc=qw_parc(ii,k)-ql_parc_dil(ii,k)

! Water loading in plume  ( option 0 water remains in plume)

         IF (dil_plume_water_load == 1) THEN

           ! water removed from parcel after condesation if > 0.001 kg/kg
           IF (ql_parc_dil(ii,k) > 0.001) THEN
             ql_parc_dil(ii,k) = 0.001
             qw_parc(ii,k) = q_vap_parc + ql_parc_dil(ii,k)
           END IF

         ELSE IF (dil_plume_water_load == 2) THEN

           ! water remaining depends on qsat environment 
           max_qw_plume = 0.5*qsat_env(ii)
           max_qw_plume = MIN (qlcrit,max_qw_plume)  ! max
           max_qw_plume = MAX (qlmin,max_qw_plume)   ! min value 
           IF (ql_parc_dil(ii,k) > max_qw_plume) THEN
             ql_parc_dil(ii,k) = max_qw_plume
             qw_parc(ii,k) = q_vap_parc + ql_parc_dil(ii,k)
           END IF           

         END IF    ! end test on dil_plume_water_load

         t_dens_parc_dil(ii,k)=t_parc_dil(ii,k)*                           &
                         (1.0+c_virtual*q_vap_parc-ql_parc_dil(ii,k))

! calculate dilute t_ref for next level
         IF (k >  1 .and. k <   wet_levels-1) THEN
           z_pr = (z_theta(ii,k+1)-z_theta(ii,k))                          &
                                       /(z_theta(ii,k)-z_theta(ii,K-1))

           th_par = t_parc_dil(ii,k)/exner_theta_levels(ii,k)
           th_ref_dil(ii) = th_par*(1.+z_pr) - th_par_km_dil(ii)*z_pr
! Check new reference sensible
           IF (th_ref_dil(ii) < 0.0) THEN
             th_ref_dil(ii) = th_par_km_dil(ii)
           END IF
           IF (th_par > 0.0) THEN   
             th_par_km_dil(ii) = th_par
           END IF
         END IF      ! k level test

       END IF   ! test on LCL
     END IF    ! test on L_dilute
 

     buoyancy(ii,k) = t_dens_parc(ii,k) - t_dens_env(ii,k)

     env_svl(ii,k)  = t_dens_env(ii,k)  + gamma_dry*z_theta(ii,K)

     par_svl(ii,k)  = t_dens_parc(ii,k) + gamma_dry*z_theta(ii,K)

     IF (L_dilute) THEN
       buoyancy_dil(ii,k) = t_dens_parc_dil(ii,k) - t_dens_env(ii,k)
       par_svl_dil(ii,k)  = t_dens_parc_dil(ii,k) + gamma_dry*z_theta(ii,K)
     END IF

   END DO     ! Loop over ii
!$OMP END DO

!$OMP END PARALLEL 
! Gradient calculations

   IF (k >= 2) THEN



! DEPENDS ON: qsat_mix
     Call qsat_mix(qsat_lev,T_parc(1,k),p_theta_lev(1,k)                  &
                                                 ,npnts,l_mixing_ratio)

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(dz, l_const, dtdz, ii)
     DO ii=1,npnts
     !-------------------------------------------------------------
     ! Find vertical gradients in parcel and environment SVL
     ! (using values from level below (i.e. K-1)).
     !-------------------------------------------------------------

       dz = z_theta(ii,K) - z_theta(ii,K-1)

       dpar_bydz(ii,k) = (par_svl(ii,k) - par_svl(ii,k-1))/dz

       denv_bydz(ii,k) = (env_svl(ii,k) - env_svl(ii,k-1))/dz

     !-----------------------------------------------------------------
     ! Temperature gradient and qsat/dz along an undilute parcel ascent
     !-----------------------------------------------------------------
 
       dtdz = (T_parc(ii,k) - T_parc(ii,k-1))/dz

       IF (T_parc(ii,k).gt.TM) THEN
         l_const=lc
       ELSE
         l_const=ls
       END IF

       dqsatdz(ii,k) = (qsat_lev(ii)/(R*T_parc(ii,k)))                    &
                             * (dtdz*repsilon*l_const/T_parc(ii,k) + g) 

     END DO    ! ii loop
!$OMP END PARALLEL DO

   ELSE

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) PRIVATE(ii)             &
!$OMP& SHARED(dpar_bydz, denv_bydz, dqsatdz, k, npnts)
     DO ii=1,npnts
       dpar_bydz(ii,k) = 0.0
       denv_bydz(ii,k) = 0.0
       dqsatdz(ii,k) = 0.0
     END DO    ! ii loop
!$OMP END PARALLEL DO

   END IF   ! test on k

END DO      ! level loop

!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook('PARCEL_ASCENT',zhook_out,zhook_handle)
END SUBROUTINE parcel_ascent
