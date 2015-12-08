! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculate wet bulb temperatute and potential temperature
!
! Subroutine Interface:
MODULE thetaw_mod

IMPLICIT NONE

CONTAINS

      SUBROUTINE Thetaw(                                                &
! In:
     & n,                                                               &
     & T,q,                                                             &
     & pressure, l_potential,                                           &
! Out:
     & TW)
     
      USE conversions_mod, ONLY: zerodegc
      USE water_constants_mod, ONLY: lc
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
! Description:
!   Thetaw generates the wet bulb temperature or wet bulb potential
!   temperature for a field of points which on a surface defined
!   by array of pressures p, with input temperature and humidity
!   already interpolated on the pressure surface.
!
! Method:
! ---------------------------------------------------------------------
! 1. Calculate wet bulb temperature at requested pressure level
! 1.1 Determine g function for first guess of wet bulb T
! 1.2 Solve for wet bulb T using Newton iteration
! [Note that solution is to within a convergence criterion for each
! point individually, using flags so that results don't change if there
! is a different set of accompanying points.]
! ---------------------------------------------------------------------
! 2. Descend wet bulb line to 1000mb to determine wet bulb theta
! 2.1 Loop over delta(pressure) intervals
! 2.2 Integrate dT/dp from requested pressure level to 1000mb using
!     Runge-Kutta method:
!     a1= Dp dT/dp(p0     ,T0     ) a2= Dp dT/dp(p0+Dp/2,T0+a1/2)
!     a3= Dp dT/dp(p0+Dp/2,T0+a2/2) a4= Dp DT/dp(p0+Dp  ,T0+a3  )
!     T1 = T0 + a1/6 + a2/3 + a3/3 + a4/6
! ---------------------------------------------------------------------

!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Physics Diagnostics
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!   Documentation: UMDP 80

! Subroutine arguments
!   Scalar arguments with intent(in):
      INTEGER, INTENT(IN) ::                                            &
     & n                   ! number of points
!   Array  arguments with intent(in):
      REAL, INTENT(IN) ::                                               &
     & T(n)                                                             &
                           ! temperature on defined surface
     &,q(n)                                                             &
                           ! humidity on defined surface
     &,pressure(n)         ! pressure on defined surface
      LOGICAL ::                                                        &
     & l_potential         ! T=Output wet bulb potential temperature
                           ! F=Output wet bulb temperature
!   Scalar arguments with intent(InOut):
!   Array  arguments with intent(InOut):
!   Scalar arguments with intent(out):
!   Array  arguments with intent(out):
      REAL, INTENT(OUT) ::                                              &
     & TW(n)               ! wet bulb potential temperature

! Local parameters:
      REAL, PARAMETER ::                                                &
     & L_coeff = 2.34e3                                                 &
                           ! latent heat temperature dependence
     &,Cp_v    = 1850.                                                  &
                           ! specific heat of water vapour
     &,Cp_d    = 1005.                                                  &
                           ! specific heat of dry air
     &,m_v     = 0.01801                                                &
                           ! Mol wt of water vapour KG/MOL
     &,m_d     = 0.02896                                                &
                           ! Mol wt of dry air      KG/MOL
     &,Rstar   = 8.314                                                  &
                           ! Universal gas constant R
     &,R_d     = Rstar/m_d                                              &
                           ! R for dry air
     &,R_v     = Rstar/m_v                                              &
                           ! R for water vapour
     &,delta_TW_tol = 0.005                                             &
                           ! tolerance for TW convergence (degrees)
     &,delta_p_max_RK = 300.E2                                          &
                               ! deltap threshold for extra R-K loops
     &,p_1000mb = 1.0E5    ! pressure at 1000mb
      INTEGER, PARAMETER ::                                             &
     & loop_max = 10                                                    &
                           ! max loops for TW iterations
     &,loop_RK  = 5        ! 5 Iterations of Runge-Kutta are
                           ! sufficient for convergence upto 250Mb.
                           ! Significantly more are required for
                           ! higher levels, but are expensive.

! Local scalars:
      INTEGER                                                           &
     & i,loop
      LOGICAL                                                           &
     & all_converged       ! test for TW convergence of Newton method
      REAL                                                              &
     & dgbydt                                                           &
                           ! derivative of g with respect to T
     &,delta_TW                                                         &
                           ! increment wet bulb T
     &,g                                                                &
                           ! function of wet bulb T
     &,dTbydP                                                           &
                           ! derivative of T with respect to p
     &,LH                                                               &
                           ! latent heat scalar
     &,a1,a2,a3,a4         ! Runge-Kutta coefficients

! Local dynamic arrays:
      LOGICAL                                                           &
     & converged(n)        ! test each TW convergence of Newton method
      REAL                                                              &
     & p(n)                                                             &
                           ! pressure
     &,delta_p(n)                                                       &
                                ! Pressure Increment
     &,qs(n)                                                            &
                           ! saturated humidity
     &,g_TW(n)                                                          &
                           ! g(TW)
     &,L(n)                                                             &
                           ! latent heat
     &,Cp_moist(n)                                                      &
                           ! specific heat of moist air
     &,TW_0(n),TW_1(n)     ! TW on pressure surface temporary arrays

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! Function & Subroutine calls:
      External Qsat
!- End of header

      IF (lhook) CALL dr_hook('THETAW',zhook_in,zhook_handle)

! ---------------------------------------------------------------------
! 1. Calculate wet bulb temperature at requested pressure level
! ---------------------------------------------------------------------

! Populate pressure array for this (and subsequent) qsat calls
      DO i=1,n
       p(i) = pressure(i)
      ENDDO

! Get saturated humidity with respect to water/ice
! DEPENDS ON: qsat
      CALL Qsat(qs,T,p,n)

! ---------------------------------------------------------------------
! 1.1 Determine g function for first guess of wet bulb T
!     g(T) = qs(T) L + C T
! ---------------------------------------------------------------------

      DO i=1,n
       TW(i) = T(i)                              ! first guess
       L(i) = LC - L_coeff*(T(I) - zerodegC)     ! latent heat
       Cp_moist(i) = (1.0-q(i))*Cp_d + q(i)*Cp_v ! specific heat
       g_TW(i) = L(i)*q(i)  + Cp_moist(i) * T(i)
       converged(i) = .false.
      ENDDO

! ---------------------------------------------------------------------
! 1.2 Solve for wet bulb T using Newton iteration,
!     where g(TW) = Cp T0 + L q = constant
! ---------------------------------------------------------------------

      DO loop=1,loop_max
       all_converged=.true.

       DO i=1,n
        IF(.NOT.converged(i)) THEN
          g      = L(i)*qs(i) + Cp_moist(i) * TW(i)
          dgbydt = L(i)*L(i)*qs(i) / (R_v*TW(i)*TW(i)) + Cp_moist(i)
          delta_TW = (g_TW(i) - g) / dgbydt
          TW(i) = TW(i) + delta_TW
          IF(abs(delta_TW) >  delta_TW_tol) THEN
             all_converged=.false.
          ELSE
             converged(i)=.true.
          ENDIF
        ENDIF ! not converged

       ENDDO ! i

       IF(all_converged) exit

! DEPENDS ON: qsat
       CALL Qsat(qs,TW,p,n)     ! revise qs values

      ENDDO ! loop

      IF (l_potential) THEN
! ---------------------------------------------------------------------
! 2. Descend wet bulb line to 1000mb to determine wet bulb theta
! ---------------------------------------------------------------------

      DO i=1,n
        delta_p(i)= (p_1000mb - pressure(i))/loop_RK
      END DO

! ---------------------------------------------------------------------
! 2.1 Loop over delta(pressure) intervals
! ---------------------------------------------------------------------

      DO loop=1,loop_RK

! Determine first guess for dT/dp at wet bulb temperature
      DO i=1,n
       TW_0(i)=TW(i)
      ENDDO

! ---------------------------------------------------------------------
! 2.2 Integrate dT/dp from requested pressure level to 1000mb using
!     Runge-Kutta method:
!     a1= Dp dT/dp(p0     ,T0     ) a2= Dp dT/dp(p0+Dp/2,T0+a1/2)
!     a3= Dp dT/dp(p0+Dp/2,T0+a2/2) a4= Dp DT/dp(p0+Dp  ,T0+a3  )
!     T1 = T0 + a1/6 + a2/3 + a3/3 + a4/6
! ---------------------------------------------------------------------

! A1 = dp * DTbyDP (p0      ,T0  )
! DEPENDS ON: qsat
      CALL Qsat(qs,TW,p,n)
      DO i=1,n
        LH = LC - L_coeff*(TW(I) - zerodegC)
        dTbydp = (LH*qs(i) + R_d*TW(i)) /                               &
     &  (p(i) * (Cp_moist(i) + LH*LH*qs(i)/(R_v*TW(i)*TW(i)) ))
        p(i) = p(i) + delta_p(i)*0.5
        a1 = delta_p(i)*dTbydp
        TW(i) = TW_0(i) + a1*0.5
        TW_1(i) = TW_0(i) + a1/6.0 ! build up Runge-Kutta estimate
      ENDDO

! A2 = dp * DTbyDP (p0 + dp/2,T0 + a1/2  )
! DEPENDS ON: qsat
      CALL Qsat(qs,TW,p,n)
      DO i=1,n
        LH = LC - L_coeff*(TW(I) - zerodegC)
        dTbydp = (LH*qs(i) + R_d*TW(i)) /                               &
     &  (p(i) * (Cp_moist(i) + LH*LH*qs(i)/(R_v*TW(i)*TW(i)) ))
        a2 = delta_p(i)*dTbydp
        TW(i) = TW_0(i) + a2*0.5
        TW_1(i) = TW_1(i) + a2/3.0 ! build up Runge-Kutta estimate
      ENDDO

! A3 = dp * DTbyDP (p0 + dp/2,T0 + a2/2  )
! DEPENDS ON: qsat
      CALL Qsat(qs,TW,p,n)
      DO i=1,n
        LH = LC - L_coeff*(TW(I) - zerodegC)
        dTbydp = (LH*qs(i) + R_d*TW(i)) /                               &
     &  (p(i) * (Cp_moist(i) + LH*LH*qs(i)/(R_v*TW(i)*TW(i)) ))
        p(i) = p(i) + delta_p(i)*0.5
        a3 = delta_p(i)*dTbydp
        TW(i) = TW_0(i) + a3
        TW_1(i) = TW_1(i) + a3/3.0 ! build up Runge-Kutta estimate
      ENDDO

! A4 = dp * DTbyDP (p0 + dp  ,T0 + a3  )
! DEPENDS ON: qsat
      CALL Qsat(qs,TW,p,n)
      DO i=1,n
        LH = LC - L_coeff*(TW(I) - zerodegC)
        dTbydp = (LH*qs(i) + R_d*TW(i)) /                               &
     &  (p(i) * (Cp_moist(i) + LH*LH*qs(i)/(R_v*TW(i)*TW(i)) ))
        a4 = delta_p(i)*dTbydp
        TW(i) = TW_1(i) + a4/6.0 ! final Runge-Kutta estimate
      ENDDO

      ENDDO ! Loop over Runge-Kutta

      END IF
      IF (lhook) CALL dr_hook('THETAW',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE Thetaw
END MODULE thetaw_mod
