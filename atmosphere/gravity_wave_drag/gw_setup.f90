! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine gw_setup to set up variables for gravity wave and 
!            flow blocking schemes
!
SUBROUTINE GW_SETUP(levels,points,u,v,rho,theta,                 &
                    nsq,ulow,vlow,rholow,nlow,psilow,psi1,modu,  &
                    rho_levels,theta_levels,                     &
                    sd_orog,grad_xx,grad_xy,grad_yy,mt_high,     &
                    slope,anis,banis,canis,mtdir,ktop,l_drag,    & 
!diagnostics
                    u_s_d,u_s_d_on,points_u_s_d,                 &
                    v_s_d,v_s_d_on,points_v_s_d,                 &
                    nsq_s_d,nsq_s_d_on,points_nsq_s_d)

!Code to set up variables for gravity wave and flow blocking schemes
!Calculates
!(i)   Buoyancy frequency squared
!(ii)  Low-level average variables (for FB and GW schemes)
!(iii) SSO variables
!(iv)  ktop
  USE earth_constants_mod, ONLY: g

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE c_gwave_mod, ONLY: nsigma, amplitude_saturation, stress_saturation,  &
                         beta_fix, frac_wl, lambdaz_min, lambdaz_max,      &
                         nsq_neutral, zav_converge, zav_iterate
  USE um_types, ONLY: real32
  IMPLICIT NONE      
! Description:
!     code to set up variables for gravity wave and flow blocking schemes
!     calculates:
!     1. Buoyancy frequency squared
!     2. Low-level average variables (for FB and GW schemes)
!     3. Sub-gridscale orography (SSO) variables      
!     4. First model level above mountain top 
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

  INTEGER, INTENT(IN) ::  &!intent(in)
    points_nsq_s_d,       &                                             
    points_u_s_d,         &                                                
    points_v_s_d         

  REAL, INTENT(IN) ::           &!block for intent(in)
    u(points,levels),           &!zonal wind
    v(points,levels),           &!meridional wind
    rho(points,levels),         &!density
    theta(points,levels),       &!potential temperature
    theta_levels(points,levels),&!height(theta_levels)
    rho_levels(points,levels),  &!height (rho_levels)
    sd_orog(points),            &!standard deviation of orography (m)
    grad_xx(points),            &!dh/dx squared gradient orography
    grad_xy(points),            &!(dh/dx)(dh/dy) gradient orography
    grad_yy(points)              !dh/dy squared gradient orography

  LOGICAL, INTENT(IN) ::  &!intent(in)
    u_s_d_on,             &                                                  
    v_s_d_on,             &                                                   
    nsq_s_d_on      
                                       
!----------------------
! Intent(out) variables
!----------------------      
   INTEGER, INTENT(OUT) ::  &!block for intent(out)
     ktop(points)            !First model level above mountain top

   REAL, INTENT(OUT) ::  &!block for intent(out)
     mt_high(points),    &!sso height (n_sigma*sd_orog)
     slope(points),      &!sso params - slope
     anis(points),       &!sso params - anisotropy
     mtdir(points),      &!sso params - angle of major axis 
     banis(points),      &!const (function of anisotropy)
     canis(points),      &!const (function of anisotropy)
     nsq(points,levels), &!brunt-vaisala freq squared
     ulow(points),       &!u averaged from z=0.5mt_high to z=mt_high  
     vlow(points),       &!v averaged from z=0.5mt_high to z=mt_high
     modu(points),       &!modulus of horizontal wind (i.e. sqrt(u^2+v^2)
     nlow(points),       &!N bulk averaged from z=0.5mt_high to z=mt_high
     psilow(points),     &!psi averaged from z=0.5mt_high to z=mt_high
     rholow(points),     &!rho averaged from z=0.5mt_high to z=mt_high
     psi1(points)         !atan(vlow/ulow)

  REAL, INTENT(OUT)  ::         &!intent(out) diags 
    u_s_d(points_u_s_d),        &!0.5mt_high-mt_high av u (ulow) diag
    v_s_d(points_v_s_d),        &!0.5mt_high-mt_high av v (vlow) diag
    nsq_s_d(points_nsq_s_d)      !0.5mt_high-mt_high av n (nlow) diag

  LOGICAL, INTENT(OUT)  ::  &
     l_drag(points)          !whether point has a non-zero stress or not

!------------------------
! Diagnostic variables
!- ----------------------
              
!----------------------
! Local variables
!----------------------      
   INTEGER :: i,k
   INTEGER :: kbot(points)
   REAL    :: smallp !Small positive number
   REAL(KIND=real32)  :: real4  !dummy variable real*4
   REAL    :: lmk,lml,lmm,&!Used to calculate the sso params
              plm,mlm      !Used to calculate the sso params
   REAL    ::       &
        dzt,&
        dzb,&
        dzu

   INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
   INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
   REAL(KIND=jprb)               :: zhook_handle
! ----------------------------------------
! initialise arrays
! ----------------------------------------
   IF (lhook) CALL dr_hook('GW_SETUP',zhook_in,zhook_handle)

   smallp = EPSILON(real4) / 100.0
   DO k = 1, levels
     DO i =1, points
        nsq(i,k)=0.0
     END DO
   END DO
   IF (u_s_d_on) THEN 
     DO i=1,points
       u_s_d(i)   = 0.0
     END DO
   END IF
   IF (v_s_d_on) THEN 
     DO i=1,points
       v_s_d(i)   = 0.0
     END DO
   END IF
   IF (nsq_s_d_on) THEN
     DO i=1,points
       nsq_s_d(i) = 0.0
     END DO
   END IF
   DO i=1,points
     ulow(i)    = 0.0
     vlow(i)    = 0.0
     rholow(i)  = 0.0 
     psilow(i)  = 0.0
     psi1(i)    = 0.0
     nlow(i)    = 0.0
     anis(i)    = 1.0
     banis(i)   = 0.0
     canis(i)   = 0.0
     slope(i)   = 0.0
     mtdir(i)   = 0.0
     ktop(i)    = 2
     kbot(i)    = 1
     mt_high(i) = nsigma*sd_orog(i) 
     l_drag(i)  = .TRUE.
!u(i,1)=0 implies at north/south pole in global model
     IF ((mt_high(i) <= 0.) .OR. (u(i,1) == 0.)) THEN
       l_drag(i) = .FALSE.  
     END IF
   END DO

! ----------------------------------------
!  calculate sso parameters and constants
! ----------------------------------------
   DO i=1,points
     IF (l_drag(i)) THEN 
       lmk          = 0.5*(grad_xx(i)+grad_yy(i))
       lml          = 0.5*(grad_xx(i)-grad_yy(i))
       lmm          = grad_xy(i)
!Denominator of equation (1) in technical documentation
       plm          = lmk + SQRT(lml**2+lmm**2)
!Numerator of equation (1) in technical documentation
       mlm          = ABS(lmk - SQRT(lml**2+lmm**2))
!Note do not take ATAN(slope) as slope is always used as TAN(slope) 
!in scheme (i.e. keep slope as non-dim number rather than convert 
!to angle
!Equation (3) in technical documentation
       slope(i)     = SQRT(plm)
       IF ((slope(i) == 0.) ) THEN
         l_drag(i) = .FALSE.  
       END IF
       IF ((ABS(lmm) <= smallp) .AND. (ABS(lml) <= smallp)) THEN
         mtdir(i) = 0.0
       ELSE
!Use ATAN2 to give correct range, as need mtdir to be in range 
![-pi/2,pi/2] (using ATAN would give range [-pi/4,pi/4]
!Equation (2) in technical documentation
         mtdir(i) = 0.5*ATAN2(lmm,lml)
       END IF
       
       IF (plm <= smallp) THEN
         anis(i) = 1.0
       ELSE 
!Equation (1) in technical documentation
         anis(i) = SQRT(mlm/plm)
       END IF

       banis(i)    = 1-0.18*anis(i) - 0.04*anis(i)**2
       canis(i)    =   0.48*anis(i) + 0.3 *anis(i)**2

     END IF !(l_drag(i)).
   END DO ! i=1 ,points  

     DO k = 2, levels-1
       DO i = 1, points 
         IF (l_drag(i)) THEN 
           nsq(i,k) = 2.*g*  ( theta(i,k) - theta(i,k-1) )                 &
                          /( ( theta(i,k) + theta(i,k-1) )                 &
                            *( theta_levels(i,k) - theta_levels(i,k-1) ) )
         END IF!(l_drag(i)).
       END DO !i = 1, points 
     END DO !k= 2, levels-1

     DO i = 1, points
!Set nsq(1)=nsq(2) as nsq is undefined on level 1
       nsq(i,1)       = nsq(i,2)
!Set nsq(levels)=nsq(levels-1) as nsq is undefined on top level
       nsq(i,levels)  = nsq(i,levels-1)
     END DO !i = 1, points 

     DO k = 1,levels 
       DO i = 1, points
         IF (l_drag(i)) THEN 
           IF (theta_levels(i,k) <= 0.5*mt_high(i)) THEN
             kbot(i) = k
           END IF  
           IF (theta_levels(i,k) < mt_high(i)) THEN
             ktop(i) = k+1
           END IF   
         END IF !(l_drag(i)).
       END DO !i = 1, points
     END DO !k = 1,levels 
     
     DO k = 2, levels-1
       DO i = 1, points
         IF (l_drag(i)) THEN 
           IF ((k > kbot(i)) .AND. (k <= ktop(i))) THEN 
             dzt = theta_levels(i,k)   -  0.5 * mt_high(i)
             dzb = theta_levels(i,k)   -  theta_levels(i,k-1)
             IF (k  ==  kbot(i)+1) THEN
               dzu = 0
               dzb = dzt
             ELSE
               dzu = theta_levels(i,k-1) - 0.5 * mt_high(i)
             END IF ! (k == bot(i)+1)
             IF (k == ktop(i)) THEN
               IF (k  ==  kbot(i)+1) THEN
                 dzb = dzt
               ELSE
                 dzt = mt_high(i) -  0.5 * mt_high(i)
                 dzb = mt_high(i) -  theta_levels(i,k-1)
               END IF !(k  ==  kbot(i)+1)
             END IF ! (k == ktop(i))
!------------------------------------------------------
! average u,v and rho from z = 0.5h to current level 
!-----------------------------------------------------
             ulow(i)  = (ulow(i)   * dzu  + u(i,k)  *dzb)/dzt
             vlow(i)  = (vlow(i)   * dzu  + v(i,k)  *dzb)/dzt
             rholow(i)= (rholow(i) * dzu  + rho(i,k)*dzb)/dzt
! bulk average N from z = 0.5mt_high tO z = mt_high
             nlow(i)  = (nlow(i)   * dzu  + nsq(i,k)*dzb)/dzt
           END IF ! (k > kbot(i) .AND k <= ktop(i))
         END IF !(l_drag(i))
       END DO !i = 1, points
     END DO  ! k = 2, levels-1

     DO i = 1, points
       IF (l_drag(i)) THEN 
         modu(i)    = SQRT(ulow(i)**2 + vlow(i)**2)
         psi1(i)    = ATAN2(vlow(i),ulow(i))
         psilow(i)  = mtdir(i) - psi1(i)
         IF (nlow(i) > 0.) THEN
           nlow(i) = SQRT(nlow(i))
         ELSE
           nlow(i)   = 0.0
           l_drag(i) = .FALSE.  
         END IF !(nlow(i) > 0.)
       END IF !(l_drag(i)).
     END DO !i=1,points
!-----------------------------------------------------------------
! 4 diagnostics
!-----------------------------------------------------------------
   IF ( u_s_d_on ) THEN
     DO i=1,points
       u_s_d(i) = ulow(i)
     END DO
   END IF

   IF ( v_s_d_on ) THEN
     DO i=1,points
       v_s_d(i) = vlow(i)
     END DO
   END IF

   IF ( nsq_s_d_on ) THEN
     DO i=1,points
       nsq_s_d(i) = nlow(i)
     END DO
   END IF  

   IF (lhook) CALL dr_hook('GW_SETUP',zhook_out,zhook_handle)
END SUBROUTINE GW_SETUP
