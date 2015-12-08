! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine gw_block to calculate the vertical profile of stress and 
!            associated wind increments due to flow blocking
!
SUBROUTINE GW_BLOCK(levels,points,dt,u,v,rho,nsq,ulow,vlow,rholow,psilow,  &
                    modulow,rho_levels,theta_levels,mt_high,sd_orog,       &
                    slope,anis,                                            &
                    mtdir,zb,banis,canis,dudt,dvdt,dtdt,                   &
                    fbcd,fcrit,l_drag,l_fb_heating,                        &
!diagnostics
                    du_dt_wake,points_du_dt_wake,du_dt_wake_on,            &
                    dv_dt_wake,points_dv_dt_wake,dv_dt_wake_on,            &
                    stress_ud,stress_ud_on,points_stress_ud,               &
                    stress_ud_p_on,                                        &
                    stress_vd,stress_vd_on,points_stress_vd,               &
                    stress_ud_wake,points_stress_ud_wake,stress_ud_wake_on,&
                    stress_vd_wake,points_stress_vd_wake,stress_vd_wake_on,&
                    fr_d, fr_d_on,points_fr_d,                             &
                    bld_d, bld_d_on,points_bld_d,                          &
                    bldt_d, bldt_d_on, points_bldt_d,                      & 
                    tausx_d,tausx_d_on,points_tausx_d,                     &
                    tausy_d,tausy_d_on,points_tausy_d)

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE c_gwave_mod, ONLY: nsigma, amplitude_saturation, stress_saturation,  &
                         beta_fix, frac_wl, lambdaz_min, lambdaz_max,      &
                         nsq_neutral, zav_converge, zav_iterate
  USE atmos_constants_mod, ONLY: cp
  IMPLICIT NONE      

! Description:
!     calculates the flow blocking stress profiles and wind increments.
!     1. calculate stress profile and wind increments for
!        the blocked flow. Calculations based on Lott and Miller (1997)
!        with modification to blocking layer calculation 
!        from Vosper et al (2009).
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
  INTEGER, INTENT(IN) :: &!block for intent(in)
    levels,              &!num of vertical levels 
    points                !num of land points

  INTEGER, INTENT(IN)  :: &! intent(in) for diags
    points_stress_ud_wake,&
    points_stress_vd_wake,&
    points_stress_ud,     &
    points_stress_vd,     &
    points_du_dt_wake,    &
    points_dv_dt_wake,    &
    points_fr_d,          &
    points_bld_d,         &
    points_bldt_d,        &
    points_tausx_d,       &
    points_tausy_d        

  REAL, INTENT(IN) ::           &!block for intent(in)
    u(points,levels),           &!zonal wind
    v(points,levels),           &!meridional wind
    rho(points,levels),         &!density
    nsq(points,levels),         &!brunt-vaisala freq squared 
    rho_levels(points,levels),  &!height (rho_levels)
    theta_levels(points,levels),&!height (theta_levels)
    slope(points),              &!sso params - slope
    anis(points),               &!sso params - anisotropy
    mtdir(points),              &!sso params - angle of major axis 
    mt_high(points),            &!sso params - height
    sd_orog(points),            &!standard deviation of orography (m)
    ulow(points),               &!u averaged from z=0.5mt_high to z=mt_high  
    vlow(points),               &!v averaged from z=0.5mt_high to z=mt_high
    modulow(points),            &!modulus of hrz wind (sqrt(ulow^2+vlow^2)
    psilow(points),             &!psi averaged from z=0.5mt_high to z=mt_high
    rholow(points),             &!rho averaged from z=0.5mt_high to z=mt_high
    banis(points),              &!const (function of anisotropy)
    canis(points),              &!const (function of anisotropy)
    fcrit,                      &!critical Froude number  
    fbcd,                       &!flow blocking drag coefficient
    dt                           !time-step

  LOGICAL, INTENT(IN) ::   &!
    l_drag(points),        &!whether point has a non-zero stress or not
    l_fb_heating            !calculate heating tendency if true

! Below are the stash flags for calculating diagnostics:
  LOGICAL, INTENT(IN) ::  &!intent(in) for diags  
    stress_ud_wake_on,    &!u wake stress
    stress_vd_wake_on,    &!v wake stress
    stress_ud_on,         &!total u stress    
    stress_ud_p_on,       &!total u stress on p-points
    stress_vd_on,         &!total v stress
    du_dt_wake_on,        &!u accel blocked flow
    dv_dt_wake_on,        &!v accel blocked flow
    fr_d_on,              &!fr_d switch
    bld_d_on,             &!bld_d switch
    bldt_d_on,            &!bldt_d switch
    tausx_d_on,           &!tausx_d switch
    tausy_d_on             !tausy_d switch
!----------------------
! Intent(out) variables
!----------------------      
  REAL, INTENT(OUT) ::     &!block for intent(out)
    dudt(points,levels),   &
    dvdt(points,levels),   &!profiles of acceleration
    dtdt(points,levels),   &!profiles of heating
    zb(points)              !depth of flow blocking layer
 
  REAL, INTENT(OUT)  ::                                 &
    stress_ud      ( points_stress_ud      , 0:levels ),&
    stress_vd      ( points_stress_vd      , 0:levels ),&   
    stress_ud_wake ( points_stress_ud_wake , 0:levels ),&
    stress_vd_wake ( points_stress_vd_wake , 0:levels ),&   
    du_dt_wake ( points_du_dt_wake , levels ),          &                
    dv_dt_wake ( points_dv_dt_wake , levels ),          &
    fr_d(points_fr_d),                                  &                    
    bld_d(points_bld_d),                                &                     
    bldt_d(points_bldt_d),                              &                      
    tausx_d(points_tausx_d),                            &               
    tausy_d(points_tausy_d)

!----------------------
! Local variables
!----------------------      
 
   REAL ::                         &!block for local variables
    psi(points,levels),            &!wind direction rel. to major axis of SSO
    local_xstress(points,0:levels),&!used in stress calculation for diags
    local_ystress(points,0:levels)

   REAL ::          &!block for local variables
    zneu(points),   &!depth of near surface neutral layer
    fav(points),    &!Froude number for calc zb (zb=max(0,mt_high(fcrit-fav))
    zav(points),    &!depth used to calc fav  
    zav1(points),   &!used to calculate zav1 
    zav_new(points),&!used in calculation of zav
    u_n(points),    &!wind speed div by buoyancy freq
    nav(points),    &!bulk averaged n^2 from z=0 to z=zav
    uav(points),    &!u averaged from z=0 to z=zav
    vav(points)      !v averaged from z=0 to z=zav

   REAL ::          &!block for local variables
    cpsi,        &!(cos(psi))**2  
    spsi,        &!(sin(psi))**2  
    ratio,       &!aspect ratio of SSO as seen by incident flow
    modu,        &!modulus of horizontal wind (i.e. sqrt(u^2+v^2)
    width,       &!mountain width seen by incident flow (i.e.sqrt((zb-z)/(z+h))
    zzd1,        &!direction of blocking drag (i.e.bcos(psi)^2+csin(psi)^2)
    dblk,        &!flow blocking drag = dblk*u
    wind,        &!wind speed resolved in direction of low-level flow. 
    dzt,         &!layer thickness used to calc average u etc.
    dzb,         &!layer thickness used to calc average u etc.
    rdt,         &!reciprocal of timestep (ie. 1/dt)
    delta_z,     &!for stress diag calc 
    dzz,         &!layer thickness used to calculate heating increment
    ududt,       &!Rate of change of kinetic energy after wind increments
    vdvdt,       &!Rate of change of kinetic energy after wind increments
    uhat,        &!u on theta level
    vhat          !v on theta level
  INTEGER :: i, ii, k, t 
  INTEGER :: ktop_fb(points)

  LOGICAL :: l_cont(points),  &
             l_cont2(points), &
             l_cont3(points)

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  REAL(KIND=jprb)               :: zhook_handle

!-------------------------------------------------------------------
!   1.0 start  preliminaries
! initialise increment and increment diagnostics
!------------------------------------------------------------
    IF (lhook) CALL dr_hook('GW_BLOCK',zhook_in,zhook_handle)
    IF (tausx_d_on) THEN
      DO i=1,points
        tausx_d(i) = 0.0
      END DO
    END IF
    IF (tausy_d_on) THEN
      DO i=1,points
        tausy_d(i) = 0.0
      END DO
    END IF
    IF (fr_d_on) THEN
      DO i = 1, points
        fr_d(i)   = 0.0
      END DO
    END IF
    IF (bld_d_on) THEN
      DO i = 1, points
        bld_d(i)  = 0.0
      END DO  
    END IF
    IF (bldt_d_on) THEN
      DO i = 1, points
        bldt_d(i) = 0.0
      END DO
    END IF
    IF (du_dt_wake_on) THEN
      DO k = 1, levels
        DO i = 1, points
          du_dt_wake(i,k)=0.0
        END DO
      END DO
    END IF
    IF (dv_dt_wake_on) THEN
      DO k = 1, levels
        DO i = 1, points
          dv_dt_wake(i,k)=0.0
        END DO
      END DO
    END IF
    IF (stress_ud_wake_on) THEN
      DO k = 0, levels
        DO i = 1, points
          stress_ud_wake(i,k)=0.0
        END DO
      END DO
    END IF
    IF (stress_vd_wake_on) THEN
      DO k = 0, levels
        DO i = 1, points
          stress_vd_wake(i,k)=0.0
        END DO
      END DO
    END IF
    IF (stress_ud_on .OR. stress_ud_p_on) THEN
      DO k = 0, levels
        DO i = 1, points
          stress_ud(i,k)=0.0
        END DO
      END DO
    END IF
    IF (stress_vd_on) THEN
      DO k = 0, levels
        DO i = 1, points
          stress_vd(i,k)=0.0
        END DO
      END DO
    END IF
     IF (stress_ud_on .OR. stress_ud_p_on .OR. stress_ud_wake_on .OR. tausx_d_on) THEN
      DO k = 0, levels
        DO i = 1, points
          local_xstress(i,k)=0.0
        END DO
      END DO
    END IF       
    IF (stress_vd_on .OR. stress_vd_wake_on .OR. tausy_d_on) THEN
      DO k = 0, levels
        DO i = 1, points
          local_ystress(i,k)=0.0
        END DO
      END DO
    END IF       
    DO k = 1, levels
      DO i = 1, points
        dudt(i,k) = 0.0
        dvdt(i,k) = 0.0
        dtdt(i,k) = 0.0
      END DO
    END DO

    DO i = 1, points
      zav_new(i) = 0.0
      zneu(i)    = 0.0
      uav(i)     = 0.0
      vav(i)     = 0.0
      nav(i)     = 0.0
      fav(i)     = -1.0
      zb(i)      = 0.0
      ktop_fb(i) = 2
      l_cont(i)  = .TRUE.
      l_cont2(i) = .TRUE.
      l_cont3(i) = .TRUE.
    END DO
    rdt = 1.0 / dt

! ------------------------------------------------
!  calculate zav following Vosper et al (2009) 
! ------------------------------------------------

! calculate depth of neutral layer - zneu
    DO k = 1, levels 
      DO i = 1, points
        IF (l_drag(i)) THEN
          IF (l_cont(i)) THEN 
            IF (nsq(i,k) < nsq_neutral) THEN
              zneu(i)    = rho_levels(i,k)
            ELSE
              l_cont(i)  = .FALSE.
            END IF !nsq(i,k) < nsq_neutral
          END IF !l_cont 
        END IF !l_drag(i)
      END DO!i = 1, points
    END DO !k = 1, levels 
 
! find max(mt_high,zneu) for zav 
    DO i = 1, points
      IF (l_drag(i)) THEN
        IF (zneu(i) > mt_high(i)) THEN 
          zav1(i) = zneu(i)
        ELSE
          zav1(i) = mt_high(i)
        END IF   
        zav(i)     = zav1(i)
      END IF!l_drag(i)
    END DO !i = 1, points

    DO ii = 1, zav_iterate
!Reset l_cont3 (logical to test if rho_levels(k) lt zav)
      DO i = 1, points
        l_cont3(i) = .TRUE.
      END DO !i = 1, points
      DO k = 2, levels-1
        DO i = 1, points
          IF (l_drag(i)) THEN
            IF (l_cont2(i)) THEN !l_cont2(i) tests if zav is converged
!-----------------------------------------------------
! need an if test for case where zav doesn't 
! converge in zav_iterate iterations?
!-----------------------------------------------------
              IF (l_cont3(i)) THEN !l_cont3(i) tests if rho_levels(k) lt zav
                IF (theta_levels(i,k) < zav(i)) THEN 
                  dzt = theta_levels(i,k)
                  IF (k == 2) THEN
                    dzb    = dzt
                  ELSE
                    dzb    = theta_levels(i,k) -  theta_levels(i,k-1)
                  END IF !(k==2)
                ELSE 
                  dzt    = zav(i)
                  IF (k == 2) THEN
                    dzb    = dzt
                  ELSE                  
                  dzb     = zav(i) - theta_levels(i,k-1)
                  ktop_fb(i) = k 
                  l_cont3(i)  = .FALSE.
                  END IF !(k==2) 
                END IF !theta_levels(i,k)<zav(i)
!-----------------------------------------------------
! average u,v and n from z = 0 to current level 
! which when k = ktop_fb become z = 0 tO z = zav(i) averages
! nav(i) calc is equivalent to bulk average n
! i.e. n^2 = sqrt(g/theta0*thetaav-theta0/zav(i))
!-----------------------------------------------------
                uav(i) = (uav(i)*theta_levels(i,k-1) + u(i,k)*dzb)  / dzt
                vav(i) = (vav(i)*theta_levels(i,k-1) + v(i,k)*dzb)  / dzt
                nav(i) = (nav(i)*theta_levels(i,k-1) + nsq(i,k)*dzb)/ dzt
              END IF    ! l_cont3
            END IF    ! l_cont2
          END IF !l_drag
        END DO ! i=1, points
      END DO   ! loop over k levels  
      DO i = 1, points
        IF (l_drag(i)) THEN
          IF (l_cont2(i)) THEN !l_cont2(i) tests if zav is converged
! resolve wind in the direction of the low-level flow
            wind    = (ulow(i)*uav(i) + vlow(i)*vav(i)) /modulow(i)  
            wind    = ABS(wind)
            IF (nav(i) > nsq_neutral) THEN
! first consider stable cases 
              u_n(i) = wind/SQRT(nav(i))
            ELSE
! now consider neutral (and near neutral) cases 
              u_n(i) = wind/SQRT(nsq_neutral)
            END IF
! limit u_n to sensible values
            u_n(i) = max(lambdaz_min,u_n(i))
            u_n(i) = min(lambdaz_max,u_n(i))
            zav_new(i) = zav1(i) + u_n(i)
! currently set zav_converge to 0.05
            IF ((zav(i) < zav_new(i)*(1.+zav_converge)) .AND.  &
               ( zav(i) > zav_new(i)*(1.-zav_converge))) THEN 
              l_cont2(i)  =  .FALSE.
            END IF ! test if zav(i) is convereged 
            zav(i)    = zav_new(i)
          END IF  ! l_cont2(i)
        END IF!l_drag(i)
      END DO!i = 1, points
    END DO  !ii= 1, zav_iterate

! -----------------------------------------
!  calculate fav and blocked layer depth 
! -----------------------------------------
    DO i = 1, points
      l_cont(i)  =  .TRUE.
    END DO
    DO i = 1, points
      IF (l_drag(i)) THEN
! calculate froude number
! prevent div by zero if nav(i) =0.
        IF (nav(i) > 0.) THEN
!limit fav with wavelength too
          fav(i) = u_n(i)/mt_high(i)
        ELSE
          fav(i) = -1.
        END IF 
! find zb and variable drag coefficient
        IF ((fav(i) > 0.) .AND. ((fcrit - fav(i)) > 0.)) THEN 
          zb(i) = mt_high(i)*(1 - fav(i)/fcrit)
        ELSE
          zb(i) = 0.
        END IF  
! find level at top of blocked layer
        ktop_fb(i) = 0
      END IF!l_drag(i)
    END DO !i = 1, points 
    DO k = 1, levels
      DO i = 1, points
        IF ((l_drag(i)) .AND. (zb(i) > 0.)) THEN
          IF (l_cont(i)) THEN 
            IF (rho_levels(i,k) >= zb(i)) THEN 
              ktop_fb(i)   = k 
              l_cont(i) = .FALSE.
            END IF  !rho_levels(i,k)>=zb(i)
! calculate psi at each level
            psi(i,k) = ATAN2 (v(i,k),u(i,k))
            psi(i,k) = mtdir(i) - psi(i,k)
          END IF  ! l_cont 
        END IF!l_drag(i)
      END DO!i = 1, points
    END DO !k=1, levels 
! --------------------------------
! calculate flow blocking drag 
! --------------------------------
    DO k = 1, levels
      DO i = 1, points
        IF ((l_drag(i)) .AND. (zb(i) > 0.)) THEN
          IF (k <  ktop_fb(i)) THEN
            cpsi  = (COS(psi(i,k)))**2
            spsi  = 1.0-cpsi
            ratio = 0.0 
!Note ratio expression revised from LM97 due to mistake in derivation in LM97
            IF ((anis(i)**2*cpsi + spsi) /=0) THEN
              ratio = SQRT((cpsi + anis(i)**2*spsi) / (anis(i)**2*cpsi + spsi))
            END IF
            IF (ratio /=0) THEN
              ratio =  2. - 1./ratio
            END IF
            IF (ratio < 0) THEN 
              ratio = 0. 
            END IF
! calculate modu using uav and vav at previous time-step         
            modu  = SQRT(uav(i)**2 + vav(i)**2)
            width = SQRT((zb(i) - rho_levels(i,k)) / &
                        (rho_levels(i,k) + sd_orog(i)))
            zzd1  = banis(i)*cpsi + canis(i)*spsi
!------------------------------------------------------------------
! calculate tendencies using partially implicit formulation.
!------------------------------------------------------------------
            dblk         = -fbcd*ratio*slope(i)/(2.*sd_orog(i)) &
                           *width*zzd1*0.5
            dblk         = ((1./(1.-dblk*modu*dt))-1.)*rdt
            dudt(i,k)    =  u(i,k) * dblk
            dvdt(i,k)    =  v(i,k) * dblk
          END IF !(k <  ktop_fb(i))
        END IF !(l_drag(i))
      END DO !Loop over points
    END DO !Loop over levels

!-----------------------------------------------------------------
! Calculate heating due to dissipation
!-----------------------------------------------------------------
    IF ( l_fb_heating ) THEN

     DO k = 1, levels-1
       DO i = 1, points
         IF ((l_drag(i)) .AND. (zb(i) > 0.)) THEN

           dzb          = theta_levels(i,k)   -  rho_levels(i,k)
           dzt          = rho_levels(i,k+1)   -  theta_levels(i,k)
           dzz          = rho_levels(i,k+1)   -   rho_levels(i,k)

!          u and v on theta_level(k)
           uhat       =  dzt * u(i,k) + dzb * u(i,k+1)
           vhat       =  dzt * v(i,k) + dzb * v(i,k+1)

!          u*du/dt abd v*dv/dt on theta_level(k)
           ududt      = uhat *( dzt * dudt(i,k) + dzb * dudt(i,k+1) )
           vdvdt      = vhat *( dzt * dvdt(i,k) + dzb * dvdt(i,k+1) )

!          dT/dt on theta_level(k)
           dtdt(i,k)  = - (ududt + vdvdt) / ( cp * dzz * dzz )

         END IF !(l_drag(i))
       END DO !Loop over points
     END DO !Loop over levels
 
    END IF !l_fb_heating


!-----------------------------------------------------------------
! 4 diagnostics
!-----------------------------------------------------------------

    IF ( bld_d_on ) THEN
      DO i=1,points
         IF (l_drag(i)) THEN
           bld_d(i) = zb(i)
         END IF 
      END DO
    END IF
    IF ( bldt_d_on ) THEN
      DO i=1,points
        IF ((dudt(i,1) /=  0.0) .OR. (dvdt(i,1) /= 0.0 )) THEN
          bldt_d(i) = 100.
        ELSE
          bldt_d(i) = 0.0
        END IF
      END DO
    END IF   
    IF ( fr_d_on ) THEN
      DO i=1,points
        IF (l_drag(i)) THEN 
          fr_d(i) = fav(i)
        END IF
      END DO
    END IF
    IF( du_dt_wake_on ) THEN
      DO k=1, levels
        DO i=1,points
          IF( l_drag(i) .and. k <= ktop_fb(i) ) THEN
            du_dt_wake(i,k) = dudt(i,k) 
          END IF
        END DO
      END DO
    END IF
    IF( dv_dt_wake_on ) THEN
      DO k=1, levels
        DO i=1,points
          IF( l_drag(i) .AND. k <= ktop_fb(i) ) THEN
            dv_dt_wake(i,k) = dvdt(i,k)
          END IF
        END DO
      END DO
    END IF

    IF( stress_ud_wake_on .OR. stress_ud_on .OR. &
        tausx_d_on .OR. stress_ud_p_on ) THEN
      DO k=levels-1, 0, -1
        DO i=1,points
          IF( l_drag(i) .AND. k <= ktop_fb(i) ) THEN
            IF ( k  ==  0 ) THEN
              delta_z = theta_levels(i,k+1)
            ELSE
              delta_z = theta_levels(i,k+1) - theta_levels(i,k)
            END IF
            local_xstress(i,k) = local_xstress(i,k+1) -       &
                                  (dudt(i,k+1)*rho(i,k+1)*delta_z) 
          ELSE 
            local_xstress(i,k) = 0.0
          END IF
        END DO
      END DO
      IF( stress_ud_wake_on) THEN
        DO k=0, levels-1
          DO i=1,points
              stress_ud_wake(i,k) = local_xstress(i,k)
          END DO
        END DO
        DO i=1,points
          stress_ud_wake(i,levels) = stress_ud_wake(i,levels-1)
        END DO
      END IF
      IF( stress_ud_on .OR. stress_ud_p_on) THEN
        DO k=0, levels-1
          DO i=1,points
             stress_ud(i,k) = local_xstress(i,k)
          END DO  
        END DO
        DO i=1,points
          stress_ud(i,levels) = stress_ud(i,levels-1)
        END DO
      END IF
      IF (tausx_d_on) THEN
        DO i=1,points
           IF( l_drag(i)) THEN
            tausx_d(i) = local_xstress(i,0) 
           END IF
        END DO
      END IF
    END IF

    IF( stress_vd_wake_on .OR. stress_vd_on .OR. tausy_d_on ) THEN
      DO k=levels-1, 0, -1
        DO i=1,points
          IF( l_drag(i) .AND. k <= ktop_fb(i) ) THEN
            IF ( k  ==  0 ) THEN
              delta_z = theta_levels(i,k+1)
            ELSE
              delta_z = theta_levels(i,k+1) - theta_levels(i,k)
            END IF
            local_ystress(i,k) = local_ystress(i,k+1) -       &
                                  (dvdt(i,k+1)*rho(i,k+1)*delta_z) 
          ELSE 
            local_ystress(i,k) = 0.0
          END IF
        END DO
      END DO
      IF( stress_vd_wake_on) THEN
        DO k=0, levels-1
          DO i=1,points
              stress_vd_wake(i,k) = local_ystress(i,k)
          END DO
        END DO
        DO i=1,points
          stress_vd_wake(i,levels) = stress_vd_wake(i,levels-1)
        END DO
      END IF
      IF( stress_vd_on ) THEN
        DO k=0, levels-1
          DO i=1,points
             stress_vd(i,k) = local_ystress(i,k)
          END DO  
        END DO
        DO i=1,points
          stress_vd(i,levels) = stress_vd(i,levels-1)
        END DO
      END IF
      IF (tausy_d_on) THEN
        DO i=1,points
           IF( l_drag(i)) THEN
            tausy_d(i) = local_ystress(i,0) 
           END IF
        END DO
      END IF
    END IF


    IF (lhook) CALL dr_hook('GW_BLOCK',zhook_out,zhook_handle)
END SUBROUTINE GW_BLOCK
