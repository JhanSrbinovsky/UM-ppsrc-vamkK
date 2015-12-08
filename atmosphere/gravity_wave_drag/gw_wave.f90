! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine gw_wave to calculate the vertical profile of gravity wave
!            stress and associated wind increments

SUBROUTINE GW_WAVE(levels,points,u,v,rho,nsq,ulow,vlow,rholow,psi1,&
                   psilow,nlow,modu,ktop,rho_levels,theta_levels,  &
                   delta_lambda,delta_phi,latitude,mt_high,sd_orog,&
                   slope,zb,banis,canis,dudt,dvdt,dtdt,dt,         &
                   l_dynbeta,l_nonhydro,l_smooth,fsat,Gsharp,      &
                   l_drag,l_gw_heating,                            & 
!diagnostics
                   du_dt_satn,points_du_dt_satn,du_dt_satn_on,     &
                   du_dt_satn_p_on,                                &
                   dv_dt_satn,points_dv_dt_satn,dv_dt_satn_on,     &
                   stress_ud,points_stress_ud ,stress_ud_on,       &
                   stress_ud_p_on,                                 &
                   stress_vd,points_stress_vd ,stress_vd_on,       &
                   stress_ud_satn,points_stress_ud_satn,           &
                   stress_ud_satn_on,                              &
                   stress_vd_satn,points_stress_vd_satn,           &
                   stress_vd_satn_on,                              &
                   tausx_d, tausx_d_on   , points_tausx_d  ,       &
                   tausy_d, tausy_d_on   , points_tausy_d)

  USE earth_constants_mod, ONLY: earth_radius
  USE conversions_mod, ONLY: pi

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE c_gwave_mod, ONLY: nsigma, amplitude_saturation, stress_saturation,  &
                         beta_fix, frac_wl, lambdaz_min, lambdaz_max,      &
                         nsq_neutral, zav_converge, zav_iterate
  USE atmos_constants_mod, ONLY: cp
  IMPLICIT NONE    

! Description:
!     calculates the gwd stress profiles and wind increments.
!     1. calculate stress profile and wind increments for
!        linear hydrostatic waves. 
!     2. stress may be deposited either over a single model 
!        level or over a vertical gravity wave wavelength.    
!     3. non-hydrostatic waves may be allowed to propagate outside 
!        of the column. 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Gravity Wave Drag
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
! Local constants

!----------------------
! Intent(in) variables
!----------------------      
  INTEGER, INTENT(IN) ::        &!block for intent(in)
    levels,                     &!num of vertical levels 
    points                       !num of land points

  INTEGER, INTENT(IN) ::        &!block for intent(in)
    ktop(points)                 !first model level above mountain top

  INTEGER, INTENT(IN) ::        &! intent(in) for diags
    points_stress_ud,           &
    points_stress_vd,           &
    points_stress_ud_satn,      &
    points_stress_vd_satn,      &
    points_du_dt_satn,          &
    points_dv_dt_satn,          &
    points_tausx_d,             &
    points_tausy_d

  REAL, INTENT(IN)    ::        &!block for intent(in)
    u(points,levels),           &!zonal wind
    v(points,levels),           &!meridional wind
    rho(points,levels),         &!density
    nsq(points,levels),         &!brunt-vaisala freq squared
    rho_levels(points,levels),  &!height(theta_levels)
    theta_levels(points,levels),&!height (rho_levels)
    slope(points),              &!sso params - slope
    mt_high(points),            &!sso params - height
    sd_orog(points),            &!standard deviation of orography (m)
    banis(points),              &!const (function of anisotropy)
    canis(points),              &!const (function of anisotropy)
    ulow(points),               &!u averaged from z=0.5mt_high to z=mt_high  
    vlow(points),               &!v averaged from z=0.5mt_high to z=mt_high
    modu(points),               &!modulus of horizontal wind ( sqrt(u^2+v^2))
    nlow(points),               &!N bulk averaged from 0.5mt_high to mt_high
    psilow(points),             &!psi averaged from z=0.5mt_high to z=mt_high
    rholow(points),             &!rho averaged from z=0.5mt_high to z=mt_high
    psi1(points),               &!atan(vlow/ulow)
    zb(points),                 &!depth of flow blocking layer
    latitude(points),           &!latitude
    delta_lambda,               &!spacing between points in the i direction.
    delta_phi,                  &!spacing between points in the j direction.
    fsat,                       &!saturation Froude number 
    Gsharp,                     &!function of mtn sharpness
    dt                           !time-step

   LOGICAL, INTENT(IN) ::      &
     l_drag(points),           &!whether point has a non-zero stress or not
     l_dynbeta,                &!dynamically adjusting beta(angle of group vel)
     l_nonhydro,               &!nonhydro scheme
     l_smooth,                 &!lambda_z smoothing of acc
     l_gw_heating               !calculate heating tendency if true

  LOGICAL, INTENT(IN)  ::      &!intent(in) for diags  
    stress_ud_on,              &!u stress
    stress_ud_p_on,            &!u stress on pressure points
    stress_vd_on,              &!v stress
    stress_ud_satn_on,         &!u satn stress
    stress_vd_satn_on,         &!v satn stress
    du_dt_satn_on,             &!u accel (saturation)
    du_dt_satn_p_on,           &!u accel (saturation) on pressure points
    dv_dt_satn_on,             &!v accel (saturation)
    tausx_d_on,                &!tausx_d switch
    tausy_d_on                  !tausy_d switch
!----------------------
! Intent(inout) variables
!----------------------      
  REAL, INTENT(INOUT) ::       &!block for intent(inout)
    dudt(points,levels),       &
    dvdt(points,levels),       &!profiles of acceleration
    dtdt(points,levels)         !profiles of heating


  REAL,INTENT(INOUT)  ::                             &!intent(inout) diags  
    stress_ud (points_stress_ud,0:levels),           &
    stress_vd (points_stress_vd,0:levels)
!----------------------
! Intent(out) variables
!----------------------      
  REAL,INTENT(OUT)    ::                             &!intent(out) diags  
    stress_ud_satn (points_stress_ud_satn, 0:levels),&
    stress_vd_satn (points_stress_vd_satn, 0:levels),&                      
    du_dt_satn (points_du_dt_satn,levels),           &               
    dv_dt_satn (points_dv_dt_satn,levels),           &                      
    tausx_d(points_tausx_d),                         &               
    tausy_d(points_tausy_d)

!----------------------
! Local variables
!----------------------      
   REAL ::                  &!block for local variables
    wave_amp(points,levels),&!Profile of wave amplitude
    lambdaz(points,levels), &!Vertical wavelength
    tau(points,levels),     &!stress profile 
    beta(points,levels),    &!ratio of vert group vel to hrz group vel
    nonhyd(points,levels),  &!used to calc frac of stress leaving grid-column
    k_wave(points),         &!horizontal wavenumber
    heff(points),           &!cut-off mtn height (h-zb)
    pd1(points),            &!Directional term 
    pd2(points),            &!Directional term 
    pdmod(points),          &!Directional term 
    tau_sfc(points),        &!used in calc of sfc stress  
    spd(points,levels),     &!wind resolved in direction of low-level stress
    n(points,levels),       &!buoyancy freq sqrt(nsq)
    deltaz(points,levels),  &!z vertical grid-length 
    u_n(points,levels),     &!wind speed div by n
    maxgwdlev,              &!max height for gwd
    deltax(points,levels),  &!x horizontal grid-length 
    deltay(points,levels),  &!y horizontal grid-length    
    uacc,                   &!dudt at level k (prior to spreading over lambdaz)
    vacc,                   &!dvdt at level k (prior to spreading over lambdaz)
    rhoav,                  &!Average density over lambdaz
    rhob,                   &!rho on theta grid at level k 
    rhob_l,                 &!rho on theta grid at level k-1 
    ub,                     &!u on theta grid at level k 
    vb,                     &!v on theta grid at level k 
    spsi,                   &!(sin(psi))**2  
    scpsi,                  &!(sin(psi))*(cos(psi))  
    cpsi1_dzdy,             &!cos(psi1)*deltaz(k)*deltay
    spsi1_dzdx,             &!sin(psi1)*deltaz(k)*deltax
    ztau,                   &!sfc stress(before directional terms)  
    zvt1,zvt2,              &
    ztemp,                  &
    amplow,                 &!Wave amplitude at level k-1
    acrit,                  &!Critical wave amplitude
    dzb,dzu,dzl,            &!layer thickness used to calc average u etc.
    uhat,                   &!u after gravity wave drag increment
    vhat,                   &!v after gravity wave drag increment
    ududt,                  &!Rate of change of kinetic energy after wind incs
    vdvdt                    !Rate of change of kinetic energy after wind incs



   REAL ::                  &!block for local variables used to calc diags
    dudt_gw(points,levels),      &!u acc due to gravity waves
    dvdt_gw(points,levels),      &!v acc due to gravity waves 
    local_xstress(points,0:levels),&!for calculation of stress diagnostics
    local_ystress(points,0:levels),&!for calculation of stress diagnostics
    dz                            !delta z on theta grid for stress diags

  INTEGER :: i, k, t, kl, kk,kmtn
  INTEGER :: kbot,khigh(points,levels),klow(points,levels)

  LOGICAL :: l_cont(points), l_cont2(points)   
 
  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!---------------------------------------------------------------------
! 1.0 initialisation
!---------------------------------------------------------------------
    IF (lhook) CALL dr_hook('GW_WAVE',zhook_in,zhook_handle)
    IF (stress_ud_satn_on) THEN
      DO k=0,levels 
        DO i=1,points
          stress_ud_satn(i,k) = 0.0
        END DO
      END DO
    END IF
    IF (stress_vd_satn_on) THEN
      DO k=0,levels 
        DO i=1,points
          stress_vd_satn(i,k) = 0.0 
        END DO
      END DO
    END IF
    IF (du_dt_satn_on .OR. du_dt_satn_p_on) THEN
      DO k=1,levels 
        DO i=1,points
          du_dt_satn(i,k) = 0.0
        END DO
      END DO
    END IF
    IF (dv_dt_satn_on) THEN
      DO k=1,levels 
        DO i=1,points
          dv_dt_satn(i,k) = 0.0
       END DO
      END DO
    END IF
    IF( stress_ud_on .OR. stress_ud_p_on .OR. stress_ud_satn_on      &
                     .OR. tausx_d_on ) THEN
      DO k=0,levels 
        DO i=1,points
          local_xstress(i,k) =0.0
        END DO
      END DO
    END IF
    IF( stress_vd_on .OR. stress_vd_satn_on .OR. tausy_d_on ) THEN
      DO k=0,levels 
        DO i=1,points
          local_ystress(i,k) =0.0
        END DO
      END DO
    END IF
    DO k=1,levels
      DO i=1,points
        tau(i,k)     = 0.0
        nonhyd(i,k)  = 0.0
        dudt_gw(i,k) = 0.0
        dvdt_gw(i,k) = 0.0
      END DO 
    END DO
    DO i = 1, points
      tau_sfc(i) = 0.0
    END DO
!
! Set maxgwdlev
!
    maxgwdlev = 1.e10

    IF ( L_nonhydro ) THEN
      DO k = 1, levels
        DO i = 1, points
          deltay(i,k) = (earth_radius+theta_levels(i,k))
          deltax(i,k) = COS(latitude(i))*delta_lambda*deltay(i,k)
          deltay(i,k) = deltay(i,k)*delta_phi
        END DO
      END DO
    END IF   

    DO i = 1, points 
      IF (l_drag(i)) THEN 
        Heff(i)   = mt_high(i)-zb(i)
        scpsi     = SIN(psilow(i))*COS(psilow(i))
        spsi      = (SIN(psilow(i)))**2  
        pd1(i)    = banis(i)-(banis(i)-canis(i))*spsi
        pd2(i)    = (banis(i)-canis(i))*scpsi
        pdmod(i)  = SQRT(pd1(i)**2+pd2(i)**2)
      END IF
    END DO

    DO k=1, levels-1 
      DO i = 1, points 
        IF (l_drag(i)) THEN
!Resolve wind into direction of sfc stress
!on theta (theta_levels) grid
          dzl          = theta_levels(i,k)    -  rho_levels(i,k)
          dzu          = rho_levels(i,k+1)   -  theta_levels(i,k)
          dzb          = rho_levels(i,k+1)   -   rho_levels(i,k)
          ub           = (dzu*u(i,k)       + dzl*u(i,k+1))/dzb
          vb           = (dzu*v(i,k)       + dzl*v(i,k+1))/dzb
          zvt1         =  ulow(i)*ub + vlow(i)*vb
          zvt2         = -vlow(i)*ub + ulow(i)*vb
          spd(i,k)     = (zvt1*pd1(i) + zvt2*pd2(i))/(modu(i)*pdmod(i))
          IF (nsq(i,k) > nsq_neutral) THEN
! first consider stable cases 
            n(i,k)   = SQRT(nsq(i,k)) 
            u_n(i,k) = spd(i,k)/n(i,k)
          ELSE
! now consider neutral (and near neutral) cases 
            u_n(i,k) = spd(i,k)/SQRT(nsq_neutral)
          END IF
! limit u_n to sensible values
          u_n(i,k)     = max(lambdaz_min,u_n(i,k))
          u_n(i,k)     = min(lambdaz_max,u_n(i,k))
          lambdaz(i,k) = frac_wl*(2.0*pi*u_n(i,k))
          deltaz(i,k)  = rho_levels(i,k+1)-rho_levels(i,k)
! set beta as beta_fix
          beta(i,k)    = beta_fix 
        END IF
      END DO
    END DO

    IF (l_dynbeta) THEN
      DO i = 1, points 
        IF (l_drag(i)) THEN 
           k_wave(i)  = slope(i)/mt_high(i) !1/L, could do 2pi/L here
        END IF
      END DO
      DO k=1, levels-1 
        DO i = 1, points 
          IF (l_drag(i)) THEN
!Calculate dynamically adjusting group velocity angle at each height
            IF (spd(i,k) /= 0.) THEN
              IF ((k_wave(i)**2*u_n(i,k)**2) < 1.) THEN
                beta(i,k) = SQRT(1.-k_wave(i)**2*u_n(i,k)**2)&
                                 /(k_wave(i)*u_n(i,k))
              ELSE
                beta(i,k) = 1.
              END IF
            ELSE 
              beta(i,k) = 1.e16
            END IF  
          END IF   
        END DO
      END DO
    END IF

    DO i = 1, points 
      IF (l_drag(i)) THEN
!-----------------------------------
! Calculate surface gravity wave stress 
! following LM97 (but with cut-off mtn)
!-----------------------------------
! here assume that maximum depth of zb is mt_high
        ztau        = rholow(i)*modu(i)*nlow(i)*(0.25*Heff(i)**2)&
                    *(slope(i)/sd_orog(i))*Gsharp
        tau_sfc(i)  = ztau*pdmod(i)
      END IF
    END DO

! set wave_amp equal to cut-off mtn height
    DO k = 1, levels-1 
      DO i = 1, points 
        IF (l_drag(i)) THEN
          IF (k < ktop(i)) THEN 
!set stress and wave amp to zero if encounter a neutral layer 
!or a critical layer below the mountain top so that wave 
!propagation loop doesn't do calc for cases where 
!there is a neutral layer or critical layer at k=ktop-1
            IF ((nsq(i,k) <= nsq_neutral) .OR.     &
               (spd(i,k+1)*spd(i,k) <= 0.0)) THEN
               tau(i,k)      = 0.0
               wave_amp(i,k) = 0.0 
            ELSE
!otherwise set wave amp and stress as follows
            wave_amp(i,k) = max(0.,Heff(i))
            tau(i,k) = tau_sfc(i)
            END IF
          ELSE IF (k == ktop(i)) THEN 
            wave_amp(i,k) = max(0.,Heff(i))
          ELSE IF (k > ktop(i)) THEN
            wave_amp(i,k) = 0.0
          END IF
        END IF
      END DO
    END DO

!-----------------------------------------------------------------
! Calculate stress lost out of grid-box due to non-hydrostatic prop.
! and deposit over height that it is lost
!-----------------------------------------------------------------
   IF (l_nonhydro) THEN 
     DO k = 1, levels-1 
       DO i = 1, points 
         IF (l_drag(i)) THEN
           cpsi1_dzdy  = ABS(COS(psi1(i)))*deltaz(i,k)*deltay(i,k)
           spsi1_dzdx  = ABS(SIN(psi1(i)))*deltaz(i,k)*deltax(i,k)
           nonhyd(i,k) = (cpsi1_dzdy + spsi1_dzdx)         &
                         /(beta(i,k)*deltax(i,k)*deltay(i,k) &
                          + cpsi1_dzdy + spsi1_dzdx)
         END IF
       END DO
     END DO
   END IF !l_nonhydro

!-----------------------------------
! Wave propagation following McFarlane
!-----------------------------------

! Start loop at mountain top 
! Just use sfc stress as launch stress from mtn top
 
    DO k=2, levels-1
      DO i = 1, points 
        IF (l_drag(i)) THEN
          IF (k >= ktop(i) .AND. theta_levels(i,k) <=maxgwdlev) THEN
            acrit       = 0.0      
            kl          = k-1
            IF ((wave_amp(i,kl)  /=  0.0) .AND. (tau(i,kl) /= 0.0)) THEN
!   All wave stress deposited if dth/dz < 0 or wind more than
!   pi/2 to direction of surface stress
              IF ((nsq(i,k) <= nsq_neutral) .OR.                &
                  (spd(i,k) <= 0.0)) THEN
                tau(i,k) = 0.0
              ELSE
! Test whether wave amplitude exceeds critical amplitude, 
! acrit, for wave breaking 
                dzl           = theta_levels(i,k)   -  rho_levels(i,k)
                dzu           = rho_levels(i,k+1)   -  theta_levels(i,k)
                dzb           = rho_levels(i,k+1)   -   rho_levels(i,k) 
                rhob          = (dzu*rho(i,k)       + dzl*rho(i,k+1))/dzb

                dzl           = theta_levels(i,k-1) -  rho_levels(i,k-1)
                dzu           = rho_levels(i,k)     -  theta_levels(i,k-1)
                dzb           = rho_levels(i,k)     -   rho_levels(i,k-1) 
                rhob_l        = (dzu*rho(i,k-1)     + dzl*rho(i,k))  /dzb
                amplow        = rhob_l*n(i,kl)*spd(i,kl)
                wave_amp(i,k) = wave_amp(i,kl)*        &
                                SQRT(amplow/(rhob*spd(i,k)*n(i,k)))
                acrit         = fsat*(spd(i,k)/n(i,k))
                IF (wave_amp(i,k) > acrit) THEN
                  wave_amp(i,k) = acrit
                END IF
                tau(i,k) = tau(i,kl)*                                   &
                  (wave_amp(i,k) / wave_amp(i,kl))**2*                  &
                  (rhob*n(i,k)*spd(i,k)) / (rhob_l*n(i,kl)*spd(i,kl)) * &
                  (1.-nonhyd(i,k))
              END IF  !  n(i,k)<= 0 and spd(i,k)<= 0
            END IF !((wave_amp(i,kl)  /=  0.0) .AND. (tau(i,kl) /= 0.0))
          END IF !  (k >= ktop(i) .AND. theta_levels(i,k) <=maxgwdlev)
        END IF !(l_drag(i))
      END DO !  i = 1, points 
    END DO !k=1, levels-1


! For this if test should i have some epsilon so we only do stress 
! calcs where neccessary and never diagnose accs due to rounding error?
    IF (l_smooth) THEN 
      DO k=1, levels-1
        DO i = 1, points 
          IF (l_drag(i)) THEN
            IF (k >= ktop(i) .AND. theta_levels(i,k) <=maxgwdlev) THEN
              kl = k-1
              IF ((wave_amp(i,kl)  /=  0.0) .AND. (tau(i,kl) /= 0.0)) THEN
                IF (tau(i,k) /= tau(i,kl)) THEN
! Find region to apply stress over (vertical wavelength)
                  l_cont(i)  = .TRUE.
                  l_cont2(i) = .TRUE.
                  khigh(i,k) = levels-1
                  klow(i,k)  = 2
                  DO kk = k, 1, -1
                    IF ((theta_levels(i,kk) <= &
                        (theta_levels(i,k)-0.5*lambdaz(i,k))) .AND. &
                        (l_cont(i))) THEN
                      klow(i,k) =  kk
                      l_cont(i) = .FALSE.
                    END IF
                  END DO !kk, level k down to surface 
                  DO kk = k, levels
                    IF ((theta_levels(i,kk) >= &
                        (theta_levels(i,k)+0.5*lambdaz(i,k))) .AND. &
                        (l_cont2(i))) THEN
                      khigh(i,k) =  kk
                      l_cont2(i) = .FALSE.
                    END IF
                  END DO !kk, level k up to model top 
                  IF(klow(i,k)  <=  ktop(i)) THEN
                    klow(i,k)    = 1
! lambdaz is hydrostatic vertical wavelength on theta levels  
                    lambdaz(i,k) = theta_levels(i,khigh(i,k)) 
                  ELSE
                    lambdaz(i,k) = theta_levels(i,khigh(i,k)) - &
                                   theta_levels(i,klow(i,k))
                  END IF !(klow(i,k) <=ktop(i))
                END IF ! (tau(i,k) /= tau(i,kl))
              END IF !((wave_amp(i,kl)  /=  0.0) .AND. (tau(i,kl) /= 0.0))
            END IF !  (k >= ktop(i) .AND. theta_levels(i,k) <=maxgwdlev)
          END IF !(l_drag(i))
        END DO !  i = 1, points 
      END DO !k=1, levels-1
    ELSE !if l_smooth is false
      DO k=1, levels-1
        DO i = 1, points 
          IF (l_drag(i)) THEN
            IF (k >= ktop(i) .AND. theta_levels(i,k) <=maxgwdlev) THEN
              kl          = k-1
              IF ((wave_amp(i,kl)  /=  0.0) .AND. (tau(i,kl) /= 0.0)) THEN
                IF (tau(i,k) /= tau(i,kl)) THEN
                  khigh(i,k)  = k
                  klow(i,k)   = k
                  lambdaz(i,k)= deltaz(i,k)
                  IF(k  ==  ktop(i)) THEN
                    klow(i,k)    = 1
                    lambdaz(i,k) = theta_levels(i,k)  
                  END IF !(k  ==  ktop(i))
                END IF !(tau(i,k) /= tau(i,kl))
              END IF !((wave_amp(i,kl)  /=  0.0) .AND. (tau(i,kl) /= 0.0))
            END IF !  (k >= ktop(i) .AND. theta_levels(i,k) <=maxgwdlev)
          END IF !(l_drag(i))
        END DO !  i = 1, points 
      END DO !k=1, levels-1
    END IF !(l_smooth)
!------------------------------------------------------------------
! Calculate drag from vertical stress divergence
! Stress in theta (theta_levels) levels so stress divergence 
! and hence acceleration are on rho (z) levels
!------------------------------------------------------------------

! need to resolve tau into x and y direction
    DO k=1, levels-1
      DO i = 1, points 
        IF (l_drag(i)) THEN
          IF (k >= ktop(i) .AND. theta_levels(i,k) <=maxgwdlev) THEN
            kl          = k-1
            IF ((wave_amp(i,kl)  /=  0.0) .AND. (tau(i,kl) /= 0.0)) THEN 
              IF (tau(i,k) /= tau(i,kl)) THEN
                ztemp=(tau(i,k) - tau(i,kl))  /  lambdaz(i,k)
                uacc =((ulow(i)*pd1(i)-vlow(i)*pd2(i))* ztemp/pdmod(i))/modu(i)
                vacc =((vlow(i)*pd1(i)+ulow(i)*pd2(i))* ztemp/pdmod(i))/modu(i)
                rhoav = 0.0
                IF (klow(i,k) == 1) THEN
                  rhoav = rho(i,1) * theta_levels(i,1)
                  DO kk = 2, khigh(i,k)
                    rhoav = rhoav + rho(i,kk)*               &
                           (theta_levels(i,kk)-theta_levels(i,kk-1))
                  END DO
                ELSE
                  DO kk = klow(i,k), khigh(i,k)
                    rhoav = rhoav + rho(i,kk)*               &
                           (theta_levels(i,kk)-theta_levels(i,kk-1))
                  END DO
                END IF !klow(i,k) == 1
                rhoav = rhoav / lambdaz(i,k)
                DO kk = klow(i,k), khigh(i,k)
                  dudt_gw(i,kk) = uacc/rhoav  + dudt_gw(i,kk)
                  dvdt_gw(i,kk) = vacc/rhoav  + dvdt_gw(i,kk)
                END DO
              END IF !tau(i,k) /= tau(i,kl)
            END IF    ! wave_amp(kl)  /=  0.0
          END IF !k >= ktop(i) and (theta_levels(i,k) <=maxgwdlev)
        END IF !(l_drag(i)
      END DO ! Loop over points
    END DO ! Loop over levels

!update arrays with total acceleration (gravity wave + flow blocking)
    DO k=1, levels-1
      DO i = 1, points 
        dudt(i,k) = dudt(i,k)  + dudt_gw(i,k)
        dvdt(i,k) = dvdt(i,k)  + dvdt_gw(i,k)
      END DO ! Loop over points
    END DO ! Loop over levels

!-----------------------------------------------------------------
! Calculate heating due to gravity wave dissipation
!-----------------------------------------------------------------
    IF ( l_gw_heating ) THEN

     DO k = 1, levels-1
       DO i = 1, points

           dzb    = theta_levels(i,k)   -  rho_levels(i,k)
           dzu    = rho_levels(i,k+1)   -  theta_levels(i,k)
           dzl    = rho_levels(i,k+1)   -   rho_levels(i,k)

!          u and v on theta_level(k)
           uhat   =  dzu * u(i,k) + dzb * u(i,k+1)
           vhat   =  dzu * v(i,k) + dzb * v(i,k+1)

!          u*du/dt abd v*dv/dt on theta_level(k)
           ududt  = uhat *( dzu * dudt_gw(i,k) + dzb * dudt_gw(i,k+1) )
           vdvdt  = vhat *( dzu * dvdt_gw(i,k) + dzb * dvdt_gw(i,k+1) )

!          dT/dt on theta_level(k)
           dtdt(i,k)  = dtdt(i,k) - (ududt + vdvdt) / ( cp * dzl * dzl )

       END DO !Loop over points
     END DO !Loop over levels

    END IF !l_gw_heating





!------------------------------------------------------------------
! diagnostics
!------------------------------------------------------------------

    IF( du_dt_satn_on .OR. du_dt_satn_p_on ) THEN
      DO k=1, levels-1
        DO i = 1, points 
          IF (l_drag(i)) THEN
            du_dt_satn(i,k) = dudt_gw(i,k)
          END IF !(l_drag(i)
        END DO ! Loop over points
      END DO ! Loop over levels
    END IF

    IF( dv_dt_satn_on ) THEN
      DO k=1, levels-1
        DO i = 1, points 
          IF (l_drag(i)) THEN
            dv_dt_satn(i,k) = dvdt_gw(i,k)
          END IF !(l_drag(i)
        END DO ! Loop over points
      END DO ! Loop over levels
    END IF

    IF( stress_ud_on .OR. stress_ud_p_on .OR. stress_ud_satn_on      &
                     .OR. tausx_d_on ) THEN
      DO k =levels-1,0,-1  
        DO i=1,points
         IF (l_drag(i)) THEN
           IF ( k  ==  0 ) THEN
             dz = theta_levels(i,k+1)
           ELSE
             dz = theta_levels(i,k+1) - theta_levels(i,k)
           END IF
           local_xstress(i,k) = local_xstress(i,k+1) -               &
                               (dudt_gw(i,k+1)*rho(i,k+1)*dz)
         END IF
        END DO
      END DO
      IF( stress_ud_on .OR. stress_ud_p_on) THEN
        DO k =0, levels-1  
          DO i=1,points
            IF (l_drag(i)) THEN
              stress_ud(i,k) = stress_ud(i,k) + local_xstress(i,k) 
            END IF
          END DO 
        END DO
        DO i=1,points
          stress_ud(i,levels) = stress_ud(i,levels-1)
        END DO
      END IF
      IF( stress_ud_satn_on ) THEN
        DO k =0, levels-1  
          DO i=1,points
            IF (l_drag(i)) THEN
              stress_ud_satn(i,k) = local_xstress(i,k) 
            END IF
          END DO 
        END DO
        DO i=1,points
          stress_ud_satn(i,levels) = stress_ud_satn(i,levels-1)
        END DO
      END IF
      IF ( tausx_d_on ) THEN
        DO i=1,points
          IF (l_drag(i)) THEN
            tausx_d(i) = tausx_d(i) + local_xstress(i,0) 
          END IF
        END DO
      END IF
    END IF
 
    IF( stress_vd_on .OR. stress_vd_satn_on .OR. tausy_d_on ) THEN
      DO k =levels-1,0,-1  
        DO i=1,points
         IF (l_drag(i)) THEN
           IF ( k  ==  0 ) THEN
             dz = theta_levels(i,k+1)
           ELSE
             dz = theta_levels(i,k+1) - theta_levels(i,k)
           END IF
           local_ystress(i,k) = local_ystress(i,k+1) -                &
                                (dvdt_gw(i,k+1)*rho(i,k+1)*dz)
         END IF
        END DO
      END DO
      IF( stress_vd_on ) THEN
        DO k =0, levels-1  
          DO i=1,points
            IF (l_drag(i)) THEN
              stress_vd(i,k) =  stress_vd(i,k) + local_ystress(i,k)
            END IF
          END DO
        END DO
        DO i=1,points
          stress_vd(i,levels) = stress_vd(i,levels-1)
        END DO
      END IF
      IF( stress_vd_satn_on ) THEN
        DO k =0, levels-1  
          DO i=1,points
            IF (l_drag(i)) THEN
              stress_vd_satn(i,k) = local_ystress(i,k)
            END IF
          END DO 
        END DO
        DO i=1,points
          stress_vd_satn(i,levels) = stress_vd_satn(i,levels-1)
        END DO
      END IF 
      IF ( tausy_d_on ) THEN
        DO i=1,points
          IF (l_drag(i)) THEN
            tausy_d(i) = tausy_d(i) + local_ystress(i,0) 
          END IF
        END DO
      END IF
    END IF

    IF (lhook) CALL dr_hook('GW_WAVE',zhook_out,zhook_handle)
END SUBROUTINE GW_WAVE
