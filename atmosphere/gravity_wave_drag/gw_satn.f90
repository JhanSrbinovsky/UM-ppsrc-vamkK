! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE GW_SATN: Calculates stress profile for linear hydrostatic
!                     Waves
!
      SUBROUTINE gw_satn(                                               &
     &   rho,r_rho_levels,r_theta_levels                                &
     &  ,theta,u,v,s_x_stress,s_y_stress,levels                         &
     &  ,points,kay,sd_orog                                             &
     &  ,s_x_orog,s_y_orog                                              &
     &  ,du_dt,dv_dt                                                    &
     &  ,k_top,k_top_max,fr,l_drag,rho_s,l_fix_gwsatn,l_gwd_40km        &
     &  ,sat_scheme,fsat,lift                                           &
! diagnostics
     &  ,stress_ud,points_stress_ud,stress_ud_on,stress_ud_p_on         &
     &  ,stress_vd,points_stress_vd,stress_vd_on                        &
     &  ,stress_ud_satn,points_stress_ud_satn,stress_ud_satn_on         &
     &  ,stress_vd_satn,points_stress_vd_satn,stress_vd_satn_on         &
     &  ,du_dt_satn,points_du_dt_satn,du_dt_satn_on,du_dt_satn_p_on     &
     &  ,dv_dt_satn,points_dv_dt_satn,dv_dt_satn_on )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      
      USE earth_constants_mod, ONLY: g
      
      USE c_gwave_mod, ONLY: nsigma, amplitude_saturation, stress_saturation,&
                         beta_fix, frac_wl, lambdaz_min, lambdaz_max,   &
                         nsq_neutral, zav_converge, zav_iterate
      IMPLICIT NONE

! Description:
!
!   Calculate stress profiles and hence drags due to linear
!   hydrostatic gravity waves. in contrast to the 3a/3b schemes,
!   linear hydrostatic waves may occur for any Froude number now.
!   When the Froude number drops below Frc (as specified in c_gwave.h),
!   the wave amplitude is reduced to be that of the air that is not
!   blocked, rather than be the full depth of the sub-grid mountains.
!   Thus, as the Froude number reduces to zero, so the linear
!   hydrostatic wave amplitude reduces to zero.
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Gravity Wave Drag
!
! Language: fortran 90
! This code is written to umdp3 v6 programming standards.
! Suitable for single column use,rotated grids


! Local constants

! subroutine arguements;

      INTEGER                                                           &
                            !,intent(in):
     & levels               ! number of model levels

      INTEGER                                                           &
                            !,intent(in):
     & points               ! number of points

      INTEGER                                                           &
                            !,intent(in):
     & k_top(points)                                                    &
                            ! model level at mountain tops
     &,k_top_max            ! max(k_top)

!
! The integers below are set in GWD_CTL2 to points if the diagnostics
! are called and to 1 if not.
!
      INTEGER                                                           &
                            !,intent(in):
     & points_stress_ud                                                 &
     &,points_stress_vd                                                 &
     &,points_stress_ud_satn                                            &
     &,points_stress_vd_satn                                            &
     &,points_du_dt_satn                                                &
     &,points_dv_dt_satn

      REAL                                                              &
                            !,intent(in):
     & r_rho_levels(points,levels)
!                           ! heights above z=0 on rho levels

      REAL                                                              &
                            !,intent(in):
     & r_theta_levels(points,levels)
!                           ! heights above z=0 on theta levels

      REAL                                                              &
                            !,intent(in):
     & rho(points,levels)                                               &
                            ! density
     &,rho_s(points)        ! surface density (calculated in gwsurf)

      REAL                                                              &
                            !,intent(in):
     & theta(points,levels) ! theta field

      REAL                                                              &
                            !,intent(in):
     & u(points,levels)                                                 &
                            ! u field
     &,v(points,levels)     ! v field

      REAL                                                              &
                            !,intent(in):
     & s_x_stress(points)                                               &
                            ! linear surface x_stress
     &,s_y_stress(points)   ! linear surface y_stress

      REAL                                                              &
                            !,intent(in):
     & s_x_orog(points)                                                 & 
                            ! 'surface' x_orog
     &,s_y_orog(points)     ! 'surface' y_orog

      REAL                                                              &
                            !,intent(in):
     & sd_orog(points)      ! standard deviation of orography

      REAL                                                              &
                            !,intent(in):
     & kay                  ! GWD Constant (m-1)

      REAL                                                              &
                            !,intent(in):
     & fr(points)           !in low level Froude number

      REAL                                                              &
                            !,intent(in):
       fsat                 ! Froude number used to scale critical
!                           ! wave amplitude for breaking

      REAL                                                              &
                              !,intent(in):
     & lift(points)           ! blocked layer depth

      REAL                                                              &
                            !,intent(out):
     & du_dt(points,levels)                                             &
                            ! total GWD du/dt
     &,dv_dt(points,levels) ! total GWD dv/dt

! Diagnostics

      REAL                                                              &
                            !,intent(out):
     & du_dt_satn(points_du_dt_satn,levels)                             &
!                           ! Saturation du/dt
     &,dv_dt_satn(points_dv_dt_satn,levels)
!                           ! Saturation dv/dt

      REAL                                                              &
                            !,intent(out):
     & stress_ud(points_stress_ud,0:levels)                             &
!                           ! u-stress diagnostic
     &,stress_vd(points_stress_vd,0:levels)
!                           ! v-stress diagnostic

      REAL                                                              &
                            !,intent(out):
     & stress_ud_satn(points_stress_ud_satn,0:levels)                   &
!                           ! x saturation stress diag
     &,stress_vd_satn(points_stress_vd_satn,0:levels)
!                           ! y saturation stress diag

      LOGICAL                                                           &
                            !,intent(in):
     & l_drag(points)       ! whether a point has a
!                           !non-zero surface stress

      LOGICAL                                                           &
                            !,intent(in):
     & l_fix_gwsatn         ! Switch at vn6.2 - if true then
!                           ! minor big fixes invoked

      LOGICAL                                                           &
                            !,intent(in):
     & l_gwd_40km           ! Switch at vn6.6 - if true then
!                           ! GWD not applied above z=40km

      INTEGER                                                           &
                            !,intent(in):
       sat_scheme           ! Switch at vn7.1
!                           ! If =1 use amplitude based saturation test
!                           ! If =0    old stress based saturation test


!
!  These logical switches determine whether a diagnostic is to be
!  calculated.
!
      LOGICAL                                                           &
                            !,intent(in):
     & stress_ud_on                                                     &
     &,stress_ud_p_on                                                   &
     &,stress_vd_on                                                     &
     &,stress_ud_satn_on                                                &
     &,stress_vd_satn_on                                                &
     &,du_dt_satn_on                                                    &
     &,du_dt_satn_p_on                                                  &
     &,dv_dt_satn_on


!-----------------------------------------------------------------
! LOCAL ARRAYS AND SCALARS
!-----------------------------------------------------------------

      INTEGER                                                           &
     & i,k                  ! loop counters in routine

      INTEGER                                                           &
     & kk,kl,ku             ! level counters in routine

      REAL                                                              &
     & x_stress(points,2)                                               &
                            ! x_stresses (layer boundaries)
     &,y_stress(points,2)   ! y_stresses (layer boundaries)

      REAL                                                              &
     & x_s_const(points)                                                &
                            ! level independent constants for
     &,y_s_const(points)    ! calculation of critical stresses.

      REAL                                                              &
       wave_amp(points,2)   ! Wave amplitude (layer boundaries)

      REAL                                                              &
     & dzkrho(points)       ! depth of level k *rho
!                           ! (for diagnostic calculations)

      REAL                                                              &
     & dzb                                                              &
                            ! height difference across rho   layer
     &,delta_z                                                          &
                            ! height difference across theta layer
     &,dzt                  ! height difference across 2 theta layers

      REAL                                                              &
     & dzu,dzl              ! height differences in half layers

      REAL                                                              &
     & ub                                                               &
                            ! u-wind at u layer boundary
     &,vb                                                               &
                            ! v-wind at u layer boundary
     &,rhob                                                             &
     &,rhob_l               ! Density at layer boundaries


      REAL                                                              &
     & n                                                                &
                            ! Brunt Vaisala frequency
     &,n_sq                                                             &
                            ! square of Brunt Vaisala frequency
     &,n_l                  ! Brunt Vaisala frequency at lower layer

      REAL                                                              &
     & c_x_stress                                                       &
                            ! critical x_stress
     &,c_y_stress                                                       &
                            ! critical y_stress
     &,c_stress                                                         &
                            ! critical stress amplitude
     &,stress                                                           &
                            ! current  stress amplitude
     &,r_stress                                                         &
                            ! inverse of current  stress amplitude
     &,s_stress_sq                                                      &
                            ! current  stress amplitude squared
     &,stress_factor        ! 1. or 1.000001 when using more
!                           ! numerically robust option

      REAL                                                              &
     & speedcalc                                                        &
     &,spd,spd_l                                                        &
                            ! Dot product calculations for speed/stress
      ,acrit                ! Critical wave amplitude for wave breaking


      REAL                                                              &
     & maxgwdlev            ! set to 40km if l_gwd_40km is true
!                           ! set to 1e10 if l_gwd_40km is false 
!                           ! (so that if test is never true)  

      INTEGER ii, i_block   ! omp block variable

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! function and subroutine calls
! none

!-------------------------------------------------------------------
!   1. preliminaries
!-------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GW_SATN',zhook_in,zhook_handle)
      i_block = 4
!
! Set critical stress_factor according to switch l_fix_gwsatn
!

      If ( L_fix_gwsatn ) Then
        stress_factor = 1.000001
      Else
        stress_factor = 1.0
      End If

!
! Set maxgwdlev according to l_gwd_40km
!
      If ( L_gwd_40km ) Then
        maxgwdlev = 40000.
      Else
        maxgwdlev = 1.e10
      End If


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(kl,ku,s_stress_sq,r_stress, dzl, &
!$OMP& dzu,ub,vb,speedcalc,dzt,n_sq,n,dzb,delta_z,kk,stress,spd,rhob,   &
!$OMP& spd_l,rhob_l,n_l,acrit,ii,c_y_stress,c_x_stress,c_stress, i, k)

      kl=1
      ku=2

!$OMP DO SCHEDULE(STATIC)
      Do i=1,points
        x_stress(i,kl) = s_x_stress(i)
        y_stress(i,kl) = s_y_stress(i)

!       Set wave amplitude to be the cut-off mountain height
        if(l_drag(i)) then
         wave_amp(i,kl)=max(0.,nsigma*sd_orog(i)-lift(i))
        else
        wave_amp(i,kl)=0.
        endif

      End do
!$OMP END DO

!
!  Calculate x_s_const.
!  This includes all terms in the calculation of the critical
!  stress that are independent of model level
!

!$OMP DO SCHEDULE(DYNAMIC)
      Do i=1,points
        If ( l_drag(i) .and. kay  /=  0.0 ) Then
          s_stress_sq = sqrt( s_x_stress(i) * s_x_stress(i)             &
     &                 + s_y_stress(i) * s_y_stress(i) )
!
!  2 is from rho average in following loop over levels
!  fr is passed in from gw_satn and scales the critical stress
!  so that it is equal to the surface stress if rho, U and N
!  were the same as in the surface stress calculation
!
          r_stress = 1. / ( 2. *  s_stress_sq**3 * nsigma * sd_orog(i)  &
     &                      * nsigma * sd_orog(i) * fr(i) * fr(i)     )
          x_s_const(i) = s_x_orog(i) * r_stress
          y_s_const(i) = s_y_orog(i) * r_stress
        Else
          x_s_const(i) = 0.0
          y_s_const(i) = 0.0
        End if
      End do
!$OMP END DO


      If( stress_ud_on .or. stress_ud_p_on ) Then
!$OMP DO SCHEDULE(STATIC)
        Do i=1,points
          stress_ud(i,0) = s_x_stress(i)
        End do
!$OMP END DO
      End if

      If( stress_vd_on ) Then
!$OMP DO SCHEDULE(STATIC)
        Do i=1,points
          stress_vd(i,0) = s_y_stress(i)
        End do
!$OMP END DO
      End if

      If( stress_ud_satn_on ) Then
!$OMP DO SCHEDULE(STATIC)
        Do i=1,points
          stress_ud_satn(i,0) = s_x_stress(i)
        End do
!$OMP END DO
      End if

      If( stress_vd_satn_on ) Then
!$OMP DO SCHEDULE(STATIC)
        Do i=1,points
          stress_vd_satn(i,0) = s_y_stress(i)
        End do
!$OMP END DO
      End if

      If ( sat_scheme == Stress_saturation) Then
!------------------------------------------------------------------
!   2(a) Use critical stress for wave saturation test
!------------------------------------------------------------------

!------------------------------------------------------------------
!        Loop over levels calculating critical stresses at each theta
!        level and reducing the actual stress to that level if it
!        exceeds the critical stress. The change in stress
!        across the model level is then used to determine the
!        deceleration of the wind.
!
!        note: k is looping over rho levels
!------------------------------------------------------------------

!$OMP DO SCHEDULE(DYNAMIC)
        Do ii=1, points, i_block
          kl=1
          ku=2
          Do k=2,levels-1 ! min value of k_top is 2 so first check must
!                     ! be at the top of that level

            Do i=ii,MIN(ii+i_block-1,points)
            
              x_stress(i,ku) = x_stress(i,kl)
              y_stress(i,ku) = y_stress(i,kl)
              
!
! First level for wave breaking test is the first
! one above the mountain tops.
! Top level for wavebreaking now either 40km or level below model lid. 
!
              If ( k>=k_top(i) .and. r_theta_levels(i,k)<=maxgwdlev ) Then
              
                If  (  ( x_stress(i,kl)  /=  0.0)                       &
                     &         .or.( y_stress(i,kl)  /=  0.0) )   Then
                  
                  dzl    = r_theta_levels(i,k  ) -   r_rho_levels(i,k  )
                  dzu    = r_rho_levels  (i,k+1) - r_theta_levels(i,k  )
! ub here actually is ub * dzb. Only divide speedcalc by dzb if
! needed in if test below
                  ub         = dzu*u(i,k)       + dzl*u(i,k+1)
                  vb         = dzu*v(i,k)       + dzl*v(i,k+1)
                  speedcalc  = ub*s_x_stress(i) + vb*s_y_stress(i)
!
!   All wave stress deposited if dth/dz < 0 or wind more than
!   pi/2 to direction of surface stress
!
                  If ( theta(i,k+1)  <=  theta(i,k-1) .or.              &
                    &         speedcalc  <=  0.0                 ) Then
                    
                    x_stress(i,ku) = 0.0
                    y_stress(i,ku) = 0.0
                    
                  Else
                    
                    dzt        = r_theta_levels(i,k+1) -                &
                         &                       r_theta_levels(i,k-1)
                    n_sq       = g * ( theta(i,k+1) - theta(i,k-1) )    &
                         &                         / ( theta(i,k) * dzt )
                    n          = SQRT( n_sq )
                    
                    dzb        = r_rho_levels(i,k+1) - r_rho_levels(i,k)
                    speedcalc  = speedcalc / dzb
                    
!
! rho here is just the average of the full level values. Should
! interpolate in z if this error wasn't so relatively small!
!

                    c_y_stress =  kay*(rho(i,k)+rho(i,k+1))*            &
                         &       stress_factor * speedcalc**3/n
                    c_x_stress = x_s_const(i)   * c_y_stress
                    c_y_stress = y_s_const(i)   * c_y_stress
                    
                    stress     = x_stress(i,kl) * x_stress(i,kl) +      &
                         &             y_stress(i,kl) * y_stress(i,kl)
                    c_stress   = c_x_stress     * c_x_stress     +      &
                         &             c_y_stress     * c_y_stress
                    
                    If ( stress  >   c_stress ) Then
                      x_stress(i,ku) = x_stress(i,kl)                   &
                           &    *sqrt(c_stress/stress)
                      y_stress(i,ku) = y_stress(i,kl)                   &
                           &    *sqrt(c_stress/stress)
                    Else
                      x_stress(i,ku) = x_stress(i,kl)
                      y_stress(i,ku) = y_stress(i,kl)
                    End if!      stress  <   c_stress
                    
                  End if    !  n_sq > 0 and speedcalc > 0

!------------------------------------------------------------------
! Calculate drag from vertical stress convergence
! Note that any stress drop at the first level above the mountains
! is applied uniformly from that level down to the ground.
!------------------------------------------------------------------
                  If( k  ==  k_top(i) ) Then
                    delta_z = r_theta_levels(i,k)
                    If ( L_fix_gwsatn ) Then
                      delta_z = delta_z * rho_s(i)/ rho(i,k)
                    End If
                  Else
                    delta_z = r_theta_levels(i,k) - r_theta_levels(i,k-1)
                  End if
                  
                  du_dt(i,k) = (x_stress(i,ku) - x_stress(i,kl))        &
                       &                    /(delta_z*rho(i,k))
                  dv_dt(i,k) = (y_stress(i,ku) - y_stress(i,kl))        &
                       &                    /(delta_z*rho(i,k))
                End if    ! stress x or y ne 0
                
              End if    ! k>=k_top(i) and r_theta_levels(i,k)<=maxgwdlev
              
            End do ! Loop over points
        

! diagnostics
            If( stress_ud_on .or. stress_ud_p_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                If ( k  >=  k_top(i)   ) Then
                  stress_ud(i,k) = x_stress(i,ku)
                End if
              End do
            End if
            
            If( stress_vd_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                If ( k  >=  k_top(i)   ) Then
                  stress_vd(i,k) = y_stress(i,ku)
                End if
              End do
            End if
            
            If( stress_ud_satn_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                If ( k  >=  k_top(i)   ) Then
                  stress_ud_satn(i,k) = x_stress(i,ku)
                End if
              End do
            End if
            
            If( stress_vd_satn_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                If ( k  >=  k_top(i)   ) Then
                  stress_vd_satn(i,k) = y_stress(i,ku)
                End if
              End do
            End if
            
            If( du_dt_satn_on .or. du_dt_satn_p_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                du_dt_satn(i,k) = du_dt(i,k)
              End do
            End if
            
            If( dv_dt_satn_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                dv_dt_satn(i,k) = dv_dt(i,k)
              End do
            End if
            
! swap storage for lower and upper layers
            kk=kl
            kl=ku
            ku=kk
            
          End do !  Loop over levels
        End do
!$OMP END DO

        
      Else If (sat_scheme == Amplitude_saturation) Then
!------------------------------------------------------------------
!   2(b) Use wave amplitude for saturation test
!------------------------------------------------------------------

!------------------------------------------------------------------
!        Loop over levels calculating wave amplitude at each theta
!        level and reducing the amplitude (and stress) if it exceeds
!        exceeds the critical amplitude. 
!-------------------------------------------------------------------

!$OMP DO SCHEDULE(DYNAMIC)
        Do ii=1, points, i_block
          Do k=2,levels-1 ! min value of k_top is 2 so first check must
!                     ! be at the top of that level

            Do i=ii,MIN(ii+i_block-1,points)
              
              x_stress(i,ku) = x_stress(i,kl)
              y_stress(i,ku) = y_stress(i,kl)
              wave_amp(i,ku) = wave_amp(i,kl)
            
!
! First level for wave breaking test is the first
! one above the mountain tops.
! Top level for wavebreaking now either 40km or level below model lid. 
!
              If ( k>=k_top(i) .and. r_theta_levels(i,k)<=maxgwdlev ) Then
                
                If  (  ( x_stress(i,kl)  /=  0.0)                       &
                     .or.( y_stress(i,kl)  /=  0.0) )   Then
                  
                  dzl        = r_theta_levels(i,k ) - r_rho_levels(i,k)
                  dzu        = r_rho_levels(i,k+1) - r_theta_levels(i,k)
! ub here actually is ub * dzb. Only divide speedcalc by dzb if
! needed in if test below
                  ub         = dzu*u(i,k)       + dzl*u(i,k+1)
                  vb         = dzu*v(i,k)       + dzl*v(i,k+1)
                  speedcalc  = ub*s_x_stress(i) + vb*s_y_stress(i)

!
!   All wave stress deposited if dth/dz < 0 or wind more than
!   pi/2 to direction of surface stress
!
                  If ( theta(i,k+1)  <=  theta(i,k-1) .or.              &
                       speedcalc  <=  0.0                 ) Then
                    
                    x_stress(i,ku) = 0.0
                    y_stress(i,ku) = 0.0
                    
                  Else
                    
                    dzt        = r_theta_levels(i,k+1) -                &
                         r_theta_levels(i,k-1)
                    n_sq       = g * ( theta(i,k+1) - theta(i,k-1) )    &
                         / ( theta(i,k) * dzt )
                    n          = SQRT( n_sq )
                    dzb        = r_rho_levels(i,k+1) - r_rho_levels(i,k)
                    speedcalc  = speedcalc / dzb
                    rhob       = (dzu*rho(i,k)   + dzl*rho(i,k+1))/dzb
                    stress     =sqrt(s_x_stress(i)*s_x_stress(i)+       &
                         s_y_stress(i)*s_y_stress(i))
                    spd        =speedcalc/stress
                    
!               Calculate density and wind resolved in surface
!               stress direction at level below
                    dzl        = r_theta_levels(i,k-1)- r_rho_levels(i,k-1)
                    dzu        = r_rho_levels(i,k)-r_theta_levels(i,k-1)
                    ub         = dzu*u(i,k-1) + dzl*u(i,k)
                    vb         = dzu*v(i,k-1) + dzl*v(i,k)
                    dzb        = r_rho_levels(i,k) - r_rho_levels(i,k-1)
                    spd_l      = (ub*s_x_stress(i)+vb*s_y_stress(i))/   &
                         (dzb*stress)
                    rhob_l     = (dzu*rho(i,k-1) + dzl*rho(i,k))/dzb
                    
                    
                    If ( theta(i,k)  <=  theta(i,k-2) .or.              &
                         spd*spd_l <= 0.    ) Then
                      
                      x_stress(i,ku) = 0.0
                      y_stress(i,ku) = 0.0
                      
                    Else

!                  Calculate Brunt Vaisala frequency at level below
                      If ( k > 2) Then
                        dzt        = r_theta_levels(i,k) -              &
                             r_theta_levels(i,k-2)
                        n_sq       = g * ( theta(i,k) - theta(i,k-2) )  &
                             / ( theta(i,k-1) * dzt )
                        n_l        = SQRT( n_sq )
                      Else
                        n_l=n
                      Endif
                      
!                  Test whether wave amplitude exceeds critical amplitude, 
!                  acrit, for wave breaking
                      wave_amp(i,ku)=wave_amp(i,kl)*                    &
                           sqrt((rhob_l*n_l*spd_l)/(rhob*n*spd))
                      acrit=(spd/n)*fsat
                      If(wave_amp(i,ku) > acrit) wave_amp(i,ku)=acrit
                      x_stress(i,ku)=x_stress(i,kl)*                    &
                           (wave_amp(i,ku)/wave_amp(i,kl))**2*          &
                           (rhob*n*spd)/(rhob_l*n_l*spd_l)
                      y_stress(i,ku)=y_stress(i,kl)*                    &
                           (wave_amp(i,ku)/wave_amp(i,kl))**2*          &
                           (rhob*n*spd)/(rhob_l*n_l*spd_l)
                      
                    Endif !theta(i,k) <= theta(i,k-2) .or. spd*spd_l <= 0. 
                    
                  End if    !  n_sq > 0 and speedcalc > 0
                
!------------------------------------------------------------------
! Calculate drag from vertical stress convergence
!------------------------------------------------------------------
                  If( k  ==  k_top(i) ) Then
                    delta_z = r_theta_levels(i,k) * rho_s(i) / rho(i,k)
                  Else
                    delta_z = r_theta_levels(i,k) - r_theta_levels(i,k-1)
                  End if
                  
                  du_dt(i,k) = (x_stress(i,ku) - x_stress(i,kl))        &
                       /(delta_z*rho(i,k))
                  dv_dt(i,k) = (y_stress(i,ku) - y_stress(i,kl))        &
                       /(delta_z*rho(i,k))
                End if    ! stress x or y ne 0
                
              End if    ! k>=k_top(i) and  r_theta_levels(i,k)<=maxgwdlev
              
            End do ! Loop over points
            
! diagnostics
            If( stress_ud_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                If ( k  >=  k_top(i)   ) Then
                  stress_ud(i,k) = x_stress(i,ku)
                End if
              End do
            End if
            
            If( stress_vd_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)          
                If ( k  >=  k_top(i)   ) Then
                  stress_vd(i,k) = y_stress(i,ku)
                End if
              End do
            End if
            
            If( stress_ud_satn_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                If ( k  >=  k_top(i)   ) Then
                  stress_ud_satn(i,k) = x_stress(i,ku)
                End if
              End do
            End if
            
            If( stress_vd_satn_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                If ( k  >=  k_top(i)   ) Then
                  stress_vd_satn(i,k) = y_stress(i,ku)
                End if
              End do
            End if
            
            If( du_dt_satn_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                du_dt_satn(i,k) = du_dt(i,k)
              End do
            End if
            
            If( dv_dt_satn_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                dv_dt_satn(i,k) = dv_dt(i,k)
              End do
            End if
            
! swap storage for lower and upper layers
            kk=kl
            kl=ku
            ku=kk
            
          End do !  Loop over levels
        End do
!$OMP END DO

      Endif ! Test on sat_scheme 

!
! Top level diagnostics - assume now that no drag is applied at
! the top level
!
      If( stress_ud_on .or. stress_ud_p_on ) Then
!$OMP DO SCHEDULE(STATIC)
        Do i=1,points
          stress_ud(i,levels) = x_stress(i,kl)
        End do
!$OMP END DO
      End if

      If( stress_vd_on ) Then
!$OMP DO SCHEDULE(STATIC)
        Do i=1,points
          stress_vd(i,levels) = y_stress(i,kl)
        End do
!$OMP END DO
      End if

      If( stress_ud_satn_on ) Then
!$OMP DO SCHEDULE(STATIC)
        Do i=1,points
          stress_ud_satn(i,levels) = x_stress(i,kl)
        End do
!$OMP END DO
      End if

      If( stress_vd_satn_on ) Then
!$OMP DO SCHEDULE(STATIC)
        Do i=1,points
          stress_vd_satn(i,levels) = y_stress(i,kl)
        End do
!$OMP END DO
      End if

!--------------------------------------------------------------------
! 3.1 set drags and stresses below k_top
!     can only be Done now because need drag at k_top first
!--------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
      Do ii=1, points, i_block
        Do k=1,k_top_max-1
          Do i=ii,MIN(ii+i_block-1,points)
            If ( k  <   k_top(i) ) Then
              du_dt(i,k) = du_dt(i,k_top(i))
              dv_dt(i,k) = dv_dt(i,k_top(i))
            End if
          End do
        End do
      End do
!$OMP END DO

!
! saturation drag diagnostic
!
      If( du_dt_satn_on .or. du_dt_satn_p_on ) Then

!$OMP DO SCHEDULE(STATIC)
        Do ii=1, points, i_block
          Do k=1,k_top_max-1
            Do i=ii,MIN(ii+i_block-1,points)
              If ( k  <   k_top(i) ) Then
                du_dt_satn(i,k) = du_dt_satn(i,k_top(i))
              End if
            End do
          End do
        End do
!$OMP END DO
      End if

      If( dv_dt_satn_on ) Then

!$OMP DO SCHEDULE(STATIC)
        Do ii=1, points, i_block
          Do k=1,k_top_max-1
            Do i=ii,MIN(ii+i_block-1,points)
              If ( k  <   k_top(i) ) Then
                dv_dt_satn(i,k) = dv_dt_satn(i,k_top(i))
              End if
            End do
          End do
        End do
!$OMP END DO
      End if

!
!  stress diagnostics
!
      If( stress_ud_on .or. stress_ud_p_on .or. stress_ud_satn_on .or.  &
     &    stress_vd_on .or. stress_vd_satn_on) Then

!$OMP DO SCHEDULE(STATIC)
        Do ii=1, points, i_block
          Do k=1,k_top_max-1
            Do i=ii,MIN(ii+i_block-1,points)
              If ( k  <   k_top(i) ) Then
                If ( k  ==  1 ) Then
                  dzkrho(i) =  r_theta_levels(i,1)* rho(i,k)
                Else
                  dzkrho(i) = (r_theta_levels(i,k)                      &
              &                        -r_theta_levels(i,k-1))* rho(i,k)
                End if
              End if
            End do
            
            If( stress_ud_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                If ( k  <   k_top(i) ) Then
                  stress_ud(i,k) = stress_ud(i,k-1)+du_dt(i,k)*dzkrho(i)
                End if
              End do
            End if
            
            If( stress_ud_satn_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                If ( k  <   k_top(i) ) Then
                  stress_ud_satn(i,k)= stress_ud_satn(i,k-1)              &
                       &               + du_dt(i,k) * dzkrho(i)
                End if
              End do
            End if
            
            If( stress_vd_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                If ( k  <   k_top(i) ) Then
                  stress_vd(i,k)= stress_vd(i,k-1) + dv_dt(i,k) *dzkrho(i)
                End if
              End do
            End if
            
            If( stress_vd_satn_on ) Then
              Do i=ii,MIN(ii+i_block-1,points)
                If ( k  <   k_top(i) ) Then
                  stress_vd_satn(i,k)= stress_vd_satn(i,k-1)              &
                       &                    + dv_dt(i,k) * dzkrho(i)
                End if
              End do
            End if
            
          End do ! k=1,k_top_max-1
        End do
!$OMP END DO

      End if   ! stress_ud_on .or. stress_ud_p_on .or. stess_ud_satn_on .or.
!              ! stress_vd_on .or. stess_vd_satn_on

!$OMP END PARALLEL


      IF (lhook) CALL dr_hook('GW_SATN',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE gw_satn
