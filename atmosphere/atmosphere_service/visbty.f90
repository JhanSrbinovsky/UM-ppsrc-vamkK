! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate visibility (in metres) from temperature, specific humidity, 
! cloud liquid water or ThetaL and qt and Murk aerosol (if present).
!
! Subroutine Interface:
      SUBROUTINE VISBTY(                                                &
     &             p_layer,T,Q,QCL,QCF                                  &
                                                      !INPUT
     &           ,AEROSOL, PROB, RHCRIT, L_MURK                         &
                                                      !INPUT
     &           ,P_FIELD                                               &
                                                      !INPUT
     &           ,VISIBILITY)  
                                                      !OUTPUT
! Modules 
      USE water_constants_mod, ONLY: rho_water
      USE visbty_constants_mod      
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE
!
! Description:   
!   Process fields of temperature, specific humidity, cloud liquid
!   water or ThetaL and qt and Murk aerosol (if present) to give 
!   visibility in metres. 
!   Calculated at a single model level or level within surface layer 
!   e.g. screen height (1.5m)
!
! Documentation:
!    Wright, B. J., 1997: Improvements to the Nimrod Visibility
!       Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!    Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!       for Nimrod. Met. Office FR Tech Rep., No. 222.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Atmosphere Service
!
! Code description:
!   Programming standard:  Unified Model Documentation Paper No 3,
!                          Version 5, dated 08/12/92.

      INTEGER                                                           &
     &        P_FIELD                   ! IN NO. points in field.
      REAL                                                              &
     &  p_layer(p_field)                                                &
     &       ,T(P_FIELD)                                                &
                                        ! IN Temperature
     &       ,Q(P_FIELD)                                                &
                                        ! IN Qt
     &       ,QCL(P_FIELD)                                              &
                                        ! IN cloud water array.
     &       ,QCF(P_FIELD)                                              &
                                        ! IN cloud ice array.
     &       ,AEROSOL(P_FIELD)                                          &
                                        ! IN Aerosol mixing ratio(ug/kg)
     &       ,PROB                                                      &
                                      ! IN Probability level ( e.g 0.5
                                      !    corresponds to median).
     &       ,RHCRIT                  ! IN Critical RH (determines
                                      !    width of distribiution)
      LOGICAL                                                           &
     &   L_MURK                        ! IN : Aerosol present
!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
      REAL                                                              &
     &      VISIBILITY(P_FIELD)         ! OUT visibility array.
!*--------------------------------------------------------------------
!*L-------------------------------------------------------------------
! Local varables:-----------------------------------------------------
!---------------------------------------------------------------------
      REAL                                                              &
     &       QT(P_FIELD)                                                &
                                      ! total of cloud water and vapour
     &      ,P(P_FIELD)                                                 &
                                      ! pressure of level
     &      ,Qs(P_FIELD)              ! saturation vapour pressure
!*L  External subroutine called ----------------------------------------
      EXTERNAL QSAT_WAT
!---------------------------------------------------------------------
!
! Local parameter variables
!
      REAL                                                              &
     &  RHmax                                                           &
                        ! Maximum value of relative humidity
                        !  which is allowed to feed into the
                        !  calculation of the 'fog' droplet radius
     &, RHmin                                                           &
                        ! Minimum value of relative humidity which
                        !  is allowed to feed into the calculation
                        !   of the 'fog' droplet radius
     &, Weight                                                          &
                        ! Weighting on new value for iterative
                        !  solution of droplet radius
     &, Delta_radius_star                                               &
                         ! Convergence required for iterative
                        !  solution of droplet radius
     &, N                                                               &
                        ! Local number density
     &, qt_limit                                                        &
                        ! Smallest Qt value allowed
     &, radius_star_min                                                 &
                        !
     &, radius_star_max                                                 &
                        !
     &, radius_star_factor!

      INTEGER                                                           &
     &  Niterations     !  Maximum number of iteration used to
                        !   estimate the water droplet radius
      PARAMETER (    RHmin = 0.001                                      &
     &,              RHmax = 0.99                                       &
     &,             Weight = 0.75                                       &
     &,  Delta_radius_star = 0.001                                      &
     &,        Niterations = 20                                         &
     &,           qt_limit = 0.0001                                     &
     &,    radius_star_min = 1.0                                        &
     &,    radius_star_max = 1000.0                                     &
     &, radius_star_factor = 4.0 )
!
! Local workspace variables
!
       INTEGER                                                          &
     &  Point                                                           &
                         !  Loop variable for points
     &, Iteration        !  Loop variable iterations used to estimate
                         !   the water droplet radius

      REAL                                                              &
     &  m_over_m0                                                       &
                         !  Ratio of  aerosol mass mixing ratio and
                         !   the standard aerosol mass mixing ratio
     &, RecipVis                                                        &
                         !  Recipirical of the visibility
     &, radius_dry                                                      &
                         !  Radius of dry aerosol particle (m)
     &, radius                                                          &
                         !  Radius of fog droplets (m)
     &, radius_star1                                                    &
                         !  Previous estimate of water droplet radius
                         !   divided by the dry radius
     &, radius_star2                                                    &
                         !  Current best estimate of water droplet
                         !   radius divided by the dry radius
     &, radius_act                                                      &
                         !  Activation droplet radius
     &, radius_star_act                                                 &
                         !  Activation droplet rad divided by dry rad
     &, A                                                               &
                         !  A0 divided by the dry radius
     &, RH_lim                                                          &
                         !  Limited RH value (fractional)
     &, Fn                                                              &
                         !  Value of droplet radius function
     &, Deriv                                                           &
                         !  Derivative of droplet radius function
     &, radius_star_diff                                                &
                         !  Absolute value of radius_star1 minus
                         !    radius_star2
     &, RHterm                                                          &
                         !  Relative humidity term in function to be
                         !   minimised to find the droplet radius
     &, qLterm                                                          &
                         !  Liquid water term in function to be
                         !   minimised to find the droplet radius
     &, RHderiv                                                         &
                         !  Derivative of relative humidity term
     &, qLderiv                                                         &
                         !  Derivative of liquid water term
     &, bs                                                              &
                         !  Width of distribution in total water
                         !   mixing ratio space (kg/kg)
     &, qt_mod                                                          &
                         !  Modified total water value based on the
                         !   probability of the value occurring
                         !   assuming a triangular distriubtion
                         !   of width bs.
     &, qt_mod_factor    !  Factor to multiply bs to modify qt

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!     Check Prob is legal
      IF (lhook) CALL dr_hook('VISBTY',zhook_in,zhook_handle)
      IF ( Prob  <   0.0 .OR.  Prob  >   1.0 ) THEN
        Write(6,*)"INVALID PROBABILITY VALUE in VISBTY",Prob
        Prob=MIN(MAX(Prob,0.0),1.0)
      ENDIF
!     Create factor to multiply bs by to modify qt
      IF ( Prob  ==  0.5 ) THEN
        qt_mod_factor = 0.0
      ELSE IF ( Prob  >=  0.0 .AND. Prob  <   0.5 ) THEN
        qt_mod_factor = ( 1.0 - SQRT( 2.0 * Prob ) )
      ELSE IF ( Prob  >=  0.5 .AND. Prob  <=  1.0 ) THEN
        qt_mod_factor = - ( 1.0 - SQRT( 2.0 * (1.0-Prob) ) )
      END IF
! ----------------------------------------------------------------------
! For the new cloud and precipitation scheme only use the liquid content
! 1. Calculate total of water vapour and liquid water contents, P and
! limit aerosol
! ----------------------------------------------------------------------
      DO Point=1,P_FIELD
        QT(Point) = Q(Point)+QCL(Point)
      END DO

!     Calculate pressure
      DO Point=1,P_FIELD
        P(Point)=P_layer(Point)
      ENDDO

! DEPENDS ON: qsat_wat
      CALL QSAT_WAT (QS,T,P,P_FIELD)

      DO Point = 1 , P_FIELD

!-------------------------------------------------------------------
!* 2. Calculate the ratio of the aerosol mass mixing ratio to the
!*    standard mass mixing ratio, m_over_m0, and the aerosol number
!*    density, N, the dry radius, radius_dry:
!*                      p
!*                  (m )
!*           r = r0 (--)
!*            d     (m0)
!*
!*
!*    And the activation radius:
!*
!*                             1/2
!*                  (       3 )
!*                  ( 3 B0 r  )
!*           r    = ( ------d-)
!*            act   (   A0    )
!*
!*    and A (A0 divided by the dry radius).
! N.B. AEROSOL is in ug/kg, m in kg/kg
! If not available, use 10 ug/kg
!-------------------------------------------------------------------

        if (L_MURK) then
! Ensure that assumed aerosol conc. is at least Aero0:
          m_over_m0 = max(Aerosol(Point)/m0*1.0E-9,                     &
     &                             Aero0/m0*1.0E-9, 0.0001)
        else
          m_over_m0 = max(10.0/m0*1.0E-9, 0.0001)
        endif

        N = N0 * m_over_m0**(1.0-3*power)

        radius_dry = radius0 * (m_over_m0)**power
        A = A0 / radius_dry

        radius_act = SQRT( (3 * B0 * radius_dry**3) / A0 )
        radius_star_act =  radius_act/radius_dry

!-------------------------------------------------------------
!* 3. Calculate the width of the total water
!*    distribution and a modified value of total water, based
!*    on a probability.
!-------------------------------------------------------------

        bs = (1.0-RHcrit) * qs(Point)

        qt_mod = MAX( qt_limit, qt(Point)+ qt_mod_factor* bs)


!====================================================================
!* 4.  Use Newton-Raphson to iteratively improve on a first-guess
!*     droplet radius, using the droplet growth equation and the
!*     geometric relation between liquid water and droplet radius.
!====================================================================
!* 4.1 Calculate a first guess relative humidity, qt/qs, but limit it
!*     to be in the range 0.001 -> 0.999.
!*     From this calculate a first-guess normalised radius using a
!*     simplified version of the droplet growth equation:
!*
!*                              1/3
!*                (       B0   )
!*           r  = ( 1 - ------ )
!*            *   (     ln(RH) )
!*
!----------------------------------------------------------------------

        RH_lim = MIN( MAX( qt_mod/qs(Point), RHmin ) , RHmax )
        radius_star2 = (1.0-B0/LOG(RH_lim))**OneThird

!----------------------------------------------------------------------
!* 4.2 Initialise the iteration counter, the normalised radius
!*     difference, and the updated normalised radius value.
!----------------------------------------------------------------------

        Iteration = 0
        radius_star_diff = 1.0
        radius_star1 = radius_star2

        Do While ( Iteration  <   Niterations .AND.                     &
     &               radius_star_diff  >   Delta_radius_star )

!----------------------------------------------------------------------
!* 4.3 Update the iteration counter and the normalised radius value.
!----------------------------------------------------------------------

          Iteration = Iteration + 1
          radius_star1 = Weight * radius_star2                          &
     &                 + ( 1.0 - Weight ) * radius_star1

!----------------------------------------------------------------------
!* 4.4 Calculate the relative humidity term:
!*
!*                      ( A        B0   )
!*          RHterm = exp( --  -  ------ )
!*                      ( r       3     )
!*                      (  *     r  - 1 )
!*                      (         *     )
!*
!*      and its derivative with respect to the normalised radius:
!*
!*                    (                 2    )
!*                    (   A       3 B0 r     )
!*          RHderiv = ( - --  +  -------*- 2 ) * RHterm
!*                    (    2     (  3     )  )
!*                    (   r      ( r  - 1 )  )
!*                    (    *     (  *     )  )
!*
!----------------------------------------------------------------------

          If ( radius_star1  <   radius_star_act ) then
            RHterm  = EXP( A/radius_star1                               &
     &                     - B0/(radius_star1**3-1.0) )* qs(Point)
            RHderiv = - RHterm * ( -A/(radius_star1**2)                 &
     &                + (3.0*B0*radius_star1**2)                        &
     &                /(radius_star1**3-1.0)**2 )
          Else
            RHterm  = EXP( A/radius_star_act                            &
     &                     - B0/(radius_star_act**3-1.0) ) * qs(Point)
            RHderiv = 0.0
          Endif


!----------------------------------------------------------------------
!* 4.5 Calculate the liquid water mixing ratio term:
!*
!*
!*                   4             3 (  3     )
!*          qLterm = - Pi rho_w N r  ( r  - 1 )
!*                   3             d (  *     )
!*
!*      and its derivative with respect to the normalised radius:
!*
!*                                  3  2
!*          qLderiv = 4 Pi rho_w N r  r
!*                                  d  *
!*
!----------------------------------------------------------------------

          qLterm  = N * FourThirds * Pi * RHO_WATER * radius_dry**3     &
     &                * ( radius_star1**3 - 1.0 )
          qLderiv  = - N * 4.0 * Pi * RHO_WATER                         &
     &                * radius_dry**3 * radius_star1**2

!----------------------------------------------------------------------
!* 4.6 Calculate the function, Fn, and its derivative, Deriv, and
!*     an improved estimate of the normalised radius,
!*     using Newton Raphson:
!*
!*          Fn = qt - RHterm - qLterm
!*
!*          Deriv = RHderiv + qLderiv
!*
!*                          Fn
!*          r      = r  -  -----
!*           * new    *    Deriv
!*
!*     The new estimate of the normalised radius is limited lie between
!*     prescribed maximum and minimum values and within a factor of the
!*     previous value to ensure that the soultion does not diverge.
!----------------------------------------------------------------------

          Fn    = qt_mod - RHterm - qLterm
          Deriv = RHderiv + qLderiv

          radius_star2 = radius_star1 - Fn/Deriv

          IF ( radius_star2  <   radius_star_min )                      &
     &        radius_star2 = radius_star_min
          IF ( radius_star2  >   radius_star_max )                      &
     &        radius_star2 = radius_star_max
          IF ( radius_star2  >   radius_star_factor * radius_star1 )    &
     &        radius_star2 = radius_star_factor * radius_star1
          IF ( radius_star2  <   radius_star1 / radius_star_factor )    &
     &        radius_star2 = radius_star1 / radius_star_factor

!---------------------------------------------------------------------
!* 4.7 Calculate difference between the old and the new values of the
!*     normalised radius.
!---------------------------------------------------------------------

          radius_star_diff = ABS( radius_star1 - radius_star2 )

        END DO

!---------------------------------------------------------------------
!* 5.  Calculate the radius from the final normalised radius.
!---------------------------------------------------------------------

        radius = radius_star2 * radius_dry

!---------------------------------------------------------------------
!* 6. Calculate the visibility, Vis, using the equation:
!*
!*                 ln(liminal contrast)
!*           Vis = -------------2------
!*                     Beta0 N r
!*
!*    (An extra term RecipVisAir is included in the recipical of
!*     visibility to limit visibilities to 100km in clean air).
!---------------------------------------------------------------------

        RecipVis = (N * radius**2) / VisFactor + RecipVisAir
        Visibility(Point) = 1/RecipVis

      END DO

      IF (lhook) CALL dr_hook('VISBTY',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE VISBTY
