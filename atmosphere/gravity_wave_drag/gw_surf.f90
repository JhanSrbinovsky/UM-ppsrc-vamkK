! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  subroutine gw_surf to calculate surface stress vector for gwd.
!             also, calculate blocked layer surface stress.
!
      SUBROUTINE gw_surf(                                               &
     &   r_theta_levels,rho,                                            &
     &   theta,u,v,timestep,sd_orog,                                    &
     &   sigma_xx,sigma_xy,sigma_yy,                                    &
     &   s_x_lin_stress,s_y_lin_stress,                                 &
     &   s_x_wake_stress,s_y_wake_stress,                               &
     &   s_x_orog,s_y_orog,levels,                                      &
     &   points,kay,rho_s,l_taus_scale,                                 &
     &   k_top,k_top_max,lift,l_drag,fr,frc,                            &
     &   u_s_d   ,u_s_d_on   ,points_u_s_d   ,                          &
     &   v_s_d   ,v_s_d_on   ,points_v_s_d   ,                          &
     &   nsq_s_d ,nsq_s_d_on ,points_nsq_s_d ,                          &
     &   fr_d    , fr_d_on   ,points_fr_d    ,                          &
     &   bld_d   ,bld_d_on   ,points_bld_d   ,                          &
     &   bldt_d  ,bldt_d_on  ,points_bldt_d  ,                          &
     &   num_lim_d  ,num_lim_d_on  ,points_num_lim_d ,                  &
     &   num_fac_d  ,num_fac_d_on  ,points_num_fac_d ,                  &
     &   tausx_d    ,tausx_d_on    ,points_tausx_d   ,                  &
     &   tausy_d    ,tausy_d_on    ,points_tausy_d   ,                  &
     &   taus_scale_d ,    taus_scale_d_on           ,                  &
     &                                     points_taus_scale_d )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      
      USE earth_constants_mod, ONLY: g
      
      USE c_gwave_mod, ONLY: nsigma, amplitude_saturation,              &
                             stress_saturation, beta_fix, frac_wl,      &
                             lambdaz_min, lambdaz_max, nsq_neutral,     &
                             zav_converge, zav_iterate
      IMPLICIT NONE

! Description:
!              calculates the total surface stress using the
!              linear hydrostatic GWD surface stress equation
!              for a 2-d hill in the absence of rotation and
!              friction. The blocked layer depth is then used
!              to partition this stress into a linear hydrostatic
!              GWD component and a blocked flow component.
!              nsigma*sd_orog is taken to be the top of the sub-grid
!              mountains and the low level U and N are calculated
!              as averages over this layer. This U and N are used
!              in the calculation of the surface stress and in the
!              determination of the blocked layer and low level
!              Froude number.
!              From vn6.2 it is possible to allow the surface drag to
!              be dependent on the low level Froude Number. This is
!              motivated by the results of
!              Wells et al.2005 (QJRMS, 131, 1321-1338).
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Gravity Wave Drag
!
! language: fortran 90
! this code is written to umdp3 programming standards.
! suitable for single column use,rotated grids

! global variables

! local constants

!
! SUBROUTINE ARGUMENTS
!
      INTEGER                                                           &
                              !,intent(in):
     & levels                 ! number of model levels

      INTEGER                                                           &
                              !,intent(in):
     & points                 ! number of gwd points on pe

      INTEGER                                                           &
                              !,intent(out):
     & k_top(points)          ! topmost model level including mountains.
!                             ! defined to be the kth rho level in
!                             ! which the nsigma*sd_orog height lies,
!                             ! i.e.nsigma*sd_orog lies between
!                             ! r_theta_levels(k-1) and
!                             ! r_theta_levels(k).

      INTEGER                                                           &
                              !,intent(out):
     & k_top_max              ! max(k_top)

!
! Dimensions of diagnostic arrays.
! dimension = points if diag called,
!           = 1      if diag not called.
! These are set in GWD_CTL2
!
      INTEGER                                                           &
                              !,intent(in):
     & points_nsq_s_d                                                   &
     &,points_u_s_d                                                     &
     &,points_v_s_d                                                     &
     &,points_fr_d                                                      &
     &,points_bld_d                                                     &
     &,points_bldt_d                                                    &
     &,points_num_lim_d                                                 &
     &,points_num_fac_d                                                 &
     &,points_tausx_d                                                   &
     &,points_tausy_d                                                   &
     &,points_taus_scale_d

      REAL                                                              &
                              !,intent(in):
     & r_theta_levels(points,levels)
!                             ! heights on theta levels

      REAL                                                              &
                              !,intent(in):
     & rho(points,levels)     ! density on rho levels

      REAL                                                              &
                              !,intent(in):
     & theta(points,levels)   ! theta on theta levels

      REAL                                                              &
                              !,intent(in):
     & u(points,levels)       ! u on rho levels

      REAL                                                              &
                              !,intent(in):
     & v(points,levels)       ! v on rho levels

      REAL                                                              &
                              !,intent(in):
     & timestep               ! timestep (s)

      REAL                                                              &
                              !,intent(in):
     & kay                    !  gwd surface stress constant (m-1)

      REAL                                                              &
                              !,intent(in):
     & frc                    !  critical froude number below which
!                             !  hydraulic jumps are triggered

      REAL                                                              &
                              !,intent(in):
     & sd_orog(points)        ! standard deviation of the sub-grid
!                             ! orography

      REAL                                                              &
                              !,intent(in):
     & sigma_xx(points)       ! (dh/dx)^2 grid box average of the
!                             ! the hi-res source dataset. here,
!                             ! dh is the height change of the hi-res
!                             ! data and dx is the distance between
!                             ! adjacent data points.

      REAL                                                              &
                              !,intent(in):
     & sigma_xy(points)       ! (dh/dx)*(dh/dy)

      REAL                                                              &
                              !,intent(in):
     & sigma_yy(points)       ! (dh/dy)^2


      REAL                                                              &
                              !,intent(out):
     & s_x_lin_stress(points)                                           & 
                              ! 'surface' linear stress in x-dirn
     &,s_y_lin_stress(points) ! 'surface' linear stress in y-dirn

      REAL                                                              &
                              !,intent(out):
     & s_x_wake_stress(points)                                          &
                              ! wake stress in x-dirn
     &,s_y_wake_stress(points)! wake stress in y-dirn

      REAL                                                              &
                              !,intent(out):
     & s_x_orog(points)       ! 'surface' stress/orog in x-dirn

      REAL                                                              &
                              !,intent(out):
     & s_y_orog(points)       ! 'surface' stress/orog in y-dirn

      REAL                                                              &
                              !,intent(out):
     & rho_s(points)          ! density - av from z=0 to nsigma*sd_orog

      REAL                                                              &
                              !,intent(out):
     & lift(points)           ! blocked layer depth

      REAL                                                              &
                              !,intent(out):
     & fr(points)             ! low level froude number for linear
!                             ! hydrostatic waves - so nsigma*sd_orog
!                             ! is replaced by nsigma*sd_orog-lift

      REAL                                                              &
                              !,intent(out):
     & u_s_d(points_u_s_d)                                              &
                              ! 0-nsigma*sd_orog u_s diag
     &,v_s_d(points_v_s_d)    ! 0-nsigma_sd_orog v_s diag

      REAL                                                              &
                              !,intent(out):
     & nsq_s_d(points_nsq_s_d)! 0-nsigma*sd_orog nsq_s diagnostic

      REAL                                                              &
                              !,intent(out):
     & fr_d(points_fr_d)      ! Froude no. diagnostic

      REAL                                                              &
                              !,intent(out):
     & bld_d(points_bld_d)    ! blocked layer depth diagnostic

      REAL                                                              &
                              !,intent(out):
     & bldt_d(points_bldt_d)  ! %  of time blocked flow param. invoked

      REAL                                                              &
                              !,intent(out):
     & num_lim_d(points_num_lim_d)
                              ! % of time numerical limiter invoked

      REAL                                                              &
                              !,intent(out):
     & num_fac_d(points_num_fac_d)
                              ! % reduction of flow blocking stress
                              ! after numerical limiter invoked

      REAL                                                              &
                              !,intent(out):
     &tausx_d(points_tausx_d) ! x-component of total surface stress

      REAL                                                              &
                              !,intent(out):
     &tausy_d(points_tausy_d) ! y-component of total surface stress

      REAL                                                              &
                              !,intent(out):
     & taus_scale_d(points_taus_scale_d)
!                             ! scaling factor for surface stress when
!                             ! Fr dependence of surface stress invoked.

      LOGICAL                                                           &
                              !,intent(in):
     & l_taus_scale           ! true allows variation of surface stress
!                             ! on the low level Froude number

      LOGICAL                                                           &
                              !,intent(out):
     & l_drag(points)         ! whether a point has a non-zero surface
!                             ! stress or not

!
! Diagnostic switches
!
      LOGICAL                                                           &
                              !,intent(in):
     & u_s_d_on                                                         &
     &,v_s_d_on                                                         &
     &,nsq_s_d_on                                                       &
     &,fr_d_on                                                          &
     &,bld_d_on                                                         &
     &,bldt_d_on                                                        &
     &,num_lim_d_on                                                     &
     &,num_fac_d_on                                                     &
     &,tausx_d_on                                                       &
     &,tausy_d_on                                                       &
     &,taus_scale_d_on

!------------------------------------------------------------------
! LOCAL ARRAYS AND SCALARS
!------------------------------------------------------------------

      REAL                                                              &
                              !
     & u_s(points)                                                      &
                              ! u-winds - av from z=0 to nsigma*sd_orog
     &,v_s(points)            ! v-winds - av from z=0 to nsigma*sd_orog

      REAL                                                              &
     & nsq_s(points)          ! n squared av from z=0 to nsigma*sd_orog

      REAL                                                              &
                              !
     & fr_for_diag(points)    ! low level froude number

      REAL                                                              &
     & num_lim(points)        ! 0 if limiter not invoked, 100 if it is.

      REAL                                                              &
     & num_fac(points)        ! percentage reduction of the
!                             ! flow-blocking stress after the limiter
!                             ! was invoked

      REAL                                                              &
     & taus_scale(points)     ! factor by which surf stress is scaled
!                             ! when Froude no. dependency is invoked.

      INTEGER                                                           &
     & i,k                    ! loop counter in routine

      REAL                                                              &
     & speed                                                            &
                              ! wind speed / wind speed in dirn stress
     &,speedcalc                                                        &
                              ! numerator of calcuation for speed
     &,s_stress_sq                                                      &
                              ! denominater of calculation for speed
     &,s_wake_stress                                                    &
                              ! amplitude of surface wake stress
     &,s_wake_limit           ! numerical limit for s_wake_stress

      REAL                                                              &
     & n                      ! brunt_vaisala frequency

      REAL                                                              &
     & r_frc                  ! reciprocal of the critical froude number

      REAL                                                              &
     & calc                                                             &
                              ! calculation for surface stress magnitude
     &,calc1                  ! as per calc

      REAL                                                              &
     & n_squared              ! n^2 on rho levels

      REAL                                                              &
     & dzb                                                              &
                              ! relevant depth to vertical average
!                             ! of current theta level
     &,dzt                    ! height from surface to kth theta level
!                             ! or nsigma*sd_orog

      LOGICAL                                                           &
     & l_cont(points)         ! level continue

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! function and subroutine calls: none


!---------------------------------------------------------------------
! 1.0 initialisation
!---------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GW_SURF',zhook_in,zhook_handle)
      k_top_max = 1
      r_frc     = 1./frc

      Do i=1,points
        nsq_s(i)       = 0.0
      End do

      Do i=1,points
        rho_s(i)       = 0.0
      End do

      Do i=1,points
        u_s(i)         = 0.0
      End do

      Do i=1,points
        v_s(i)         = 0.0
      End do

      Do i=1,points
        l_cont(i)     =.true.
      End do

      Do i=1,points
        s_x_lin_stress(i)  = 0.0
      End do

      Do i=1,points
        s_y_lin_stress(i)  = 0.0
      End do

      Do i=1,points
        s_x_wake_stress(i) = 0.0
      End do

      Do i=1,points
        s_y_wake_stress(i) = 0.0
      End do

      Do i=1,points
        l_drag(i) = .true.
      End do

      Do i=1,points
        num_lim(i) = 0.0
      End do

      Do i=1,points
        num_fac(i) = 0.0
      End do

      Do i=1,points
        taus_scale(i) = 1.0
      End do

!------------------------------------------------------------------
! 2.
!    Calculation of the average surface quantities and
!    the nsigma*sd_orog point values. The surface
!    is the average from 0 up to nsigma*sd_orog.
!    All interpolations now done in terms of height.
!    k_top is defined to be the rho level containing nsigma*sd_orog
!    dzb=dzt when k=2 because include level 1 depth in z, but
!    don't calculate u and n at rho level 1.
!
!---------------------------------------------------------------------
      Do i=1,points
        If (sd_orog(i)  <=  0.0 ) Then
          l_drag(i) = .false.
          l_cont(i) = .false.
          k_top(i)  = 2
        End if
      End do

      Do k=2,levels-1
        Do i=1,points
          If ( l_cont(i) ) Then

            If ( r_theta_levels(i,k)  <   nsigma*sd_orog(i) ) Then
              dzt = r_theta_levels(i,k)
              If ( k  ==  2 ) Then
                dzb = dzt
              Else
                dzb = r_theta_levels(i,k) -  r_theta_levels(i,k-1)
              End if
            Else

              dzt = nsigma * sd_orog(i)
              If ( k  ==  2 ) Then
                dzb = dzt
              Else
                dzb = nsigma*sd_orog(i) -  r_theta_levels(i,k-1)
              End if

              k_top(i)   =  k
              If( k_top_max  <   k ) k_top_max=k
              l_cont(i) = .false.

            End if   ! r_theta_levels(i,k) >  nsigma*sd_orog(i)

            n_squared = 2.*g*(theta(i,k)-theta(i,k-1) )                 &
     &                 /( (theta(i,k)+theta(i,k-1)) *                   &
     &                    (r_theta_levels(i,k)-r_theta_levels(i,k-1)) )
!
! Form z=0 to current level averages which
! when k=k_top become 0-nsigma*sd_orog averages
!
            u_s(i)   = ( u_s(i)   * r_theta_levels(i,k-1) +             &
     &                   u(i,k)   * dzb    ) / dzt
            v_s(i)   = ( v_s(i)   * r_theta_levels(i,k-1) +             &
     &                   v(i,k)   * dzb    ) / dzt
            nsq_s(i) = ( nsq_s(i) * r_theta_levels(i,k-1) +             &
     &                   n_squared* dzb    ) / dzt
            rho_s(i) = ( rho_s(i) * r_theta_levels(i,k-1) +             &
     &                   rho(i,k) * dzb    ) / dzt

          End if  !            l_cont

        End do

      End do

!------------------------------------------------------------------
! 3.1 calculation of total surface stress and blocked
!     layer depth, and the split of the total stress into a linear
!     hydrostatic gravity wave stress and a blocked flow stress
!-----------------------------------------------------------------
      Do i=1,points

!
! First calculate the speed of the wind for the linear
! hydrostatic code
!
        speed = u_s(i)*u_s(i) + v_s(i)*v_s(i)

        If ( speed  <=  0.0 ) Then
          s_x_orog(i) = 0.0
          s_y_orog(i) = 0.0
        Else
          speed = SQRT(speed)
          s_x_orog(i)= (u_s(i)*sigma_xx(i)+v_s(i)*sigma_xy(i)) /speed
          s_y_orog(i)= (u_s(i)*sigma_xy(i)+v_s(i)*sigma_yy(i)) /speed
          s_stress_sq= s_x_orog(i)*s_x_orog(i)+s_y_orog(i)*s_y_orog(i)
          speedcalc  = u_s(i)*s_x_orog(i) + v_s(i)*s_y_orog(i)
          If ( s_stress_sq  <=  0.0 ) Then
            speed    = 0.0
          Else
! speed is the component of the wind perpendicular to the major
! axis of the orography.
            speed    = speedcalc / SQRT( s_stress_sq )
          End if
        End if

        If ( nsq_s(i) >  0.0 .and. speed >  0.0 .and. l_drag(i) ) Then

          n              = SQRT( nsq_s(i) )
          fr_for_diag(i) = speed / (n*nsigma*sd_orog(i))
          lift(i)        = nsigma * sd_orog(i) - r_frc * speed/n
          If ( lift(i)  <   0.0 ) Then
            lift(i)      = 0.0
          End if

!
!  fr here is for the linear hydrostatic surface stress and
!  then also for the critical stress calculation in gw_satn
!
          fr(i) = speed / ( n*(nsigma*sd_orog(i)-lift(i)))
!
! Calculate factor by which to scale surface stresses when dependence
! of surface stress on Froude number is invoked. Empirical expression
! is a simple fit to the experimental results plotted in Fig. 7 of
! Wells et al. (2005).
!
! This function was modified at vn6.6 to keep ratio of drag:linear drag 
! equal to 1 for Froude numbers>2 
!
          If ( l_taus_scale ) Then

            If ( fr_for_diag(i) <= 1. ) Then
              taus_scale(i) = 1.6*fr_for_diag(i)
            Else If ( fr_for_diag(i) <= 2. ) Then
              taus_scale(i) = 2.2 - 0.6*fr_for_diag(i)
            End if

          End if

          calc  = kay * rho_s(i) * speed**3 * taus_scale(i) /           &
     &            (n*nsigma*nsigma*sd_orog(i)*sd_orog(i)*fr(i)*fr(i))
          s_x_lin_stress(i) = s_x_orog(i) * calc
          s_y_lin_stress(i) = s_y_orog(i) * calc

!
! s_x_orog(i)*calc1 = linear hydrostatic prediction of the stress
!
          calc1              = kay*rho_s(i)* speed * n * taus_scale(i)
          s_x_wake_stress(i) = s_x_orog(i) * calc1 - s_x_lin_stress(i)
          s_y_wake_stress(i) = s_y_orog(i) * calc1 - s_y_lin_stress(i)
!
! Limit wake_stress so that numerical instability is not permitted
!
          s_wake_stress = sqrt( s_x_wake_stress(i)*s_x_wake_stress(i)   &
     &                         +s_y_wake_stress(i)*s_y_wake_stress(i))
          s_wake_limit  = 2*nsigma*sd_orog(i)*speed*rho_s(i)/timestep

          If (s_wake_stress  >   s_wake_limit) Then
            s_x_wake_stress(i) = s_x_wake_stress(i) * s_wake_limit      &
     &                                              /s_wake_stress
            s_y_wake_stress(i) = s_y_wake_stress(i) * s_wake_limit      &
     &                                              /s_wake_stress
            num_lim(i) = 100.0
            num_fac(i) = 100. * (s_wake_stress - s_wake_limit )         &
     &                          / s_wake_stress
          End If

        Else

          l_drag(i)      = .false.
          lift(i)        = 0.0
          fr_for_diag(i) = 0.0
          fr(i)          = 0.0

        End if     ! speed or n or orog  >   0.0

      End do      ! i=points


!-----------------------------------------------------------------
! 4 diagnostics
!-----------------------------------------------------------------
      If ( u_s_d_on ) Then
        Do i=1,points
          u_s_d(i) = u_s(i)
        End do
      End if

      If ( v_s_d_on ) Then
        Do i=1,points
          v_s_d(i) = v_s(i)
        End do
      End if

      If ( nsq_s_d_on ) Then
        Do i=1,points
          If ( nsq_s(i)  <   0.0 ) nsq_s(i) = 0.0
          nsq_s_d(i) = SQRT( nsq_s(i) )
        End do
      End if

      If ( bld_d_on ) Then
        Do i=1,points
          bld_d(i) = lift(i)
        End do
      End if

      If ( bldt_d_on ) Then
        Do i=1,points
          If( s_x_wake_stress(i)  /=  0.0 .or.                          &
     &        s_y_wake_stress(i) /= 0.0 )  Then
            bldt_d(i) = 100.
          Else
            bldt_d(i) = 0.0
          End if
        End do
      End if

      If ( fr_d_on ) Then
        Do i=1,points
          fr_d(i) = fr_for_diag(i)
!  limit Fr to sensible-ish values as previously had problems
!  when creating time means which include very large values.
          If (fr_d(i)  >   1000.) fr_d(i) = 1000.
        End do
      End if

      If ( num_lim_d_on ) Then
        Do i=1,points
          num_lim_d(i) = num_lim(i)
        End do
      End if

      If ( num_fac_d_on ) Then
        Do i=1,points
          num_fac_d(i) = num_fac(i)
        End do
      End if

      If ( tausx_d_on ) Then
        Do i=1,points
          tausx_d(i) = s_x_lin_stress(i) + s_x_wake_stress(i)
        End do
      End if

      If ( tausy_d_on ) Then
        Do i=1,points
          tausy_d(i) = s_y_lin_stress(i) + s_y_wake_stress(i)
        End do
      End if

      If ( taus_scale_d_on ) Then
        Do i=1,points
          taus_scale_d(i) = taus_scale(i)
        End do
      End if


      IF (lhook) CALL dr_hook('GW_SURF',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE gw_surf
