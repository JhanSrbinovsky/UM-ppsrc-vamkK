! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Invert relationship between aerosol, visibility and water content
!
! Subroutine Interface:
      SUBROUTINE VISTOQT                                                &
     &           (VISIBILITY                                            &
                                                      !INPUT
     &           ,Qs                                                    &
                                                      !INPUT
     &           ,AEROSOL                                               &
                                                      !INPUT
     &           ,L_MURK                                                &
                                                      !INPUT
     &           ,Npoints                                               &
                                                      !INPUT
     &           ,qt )      
       
! Modules      
      USE visbty_constants_mod
      USE water_constants_mod, ONLY: rho_water
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE 
!
! Description:  
!   Invert relationship between aerosol, visibility and water
!   content. This is needed for fog probability calculation.
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

!*L-------------------------------------------------------------------
! input variables-----------------------------------------------------
!---------------------------------------------------------------------
      INTEGER                                                           &
     &        Npoints             ! IN NO. points in field.
      REAL                                                              &
     &        VISIBILITY(Npoints)                                       &
                                 ! IN visibility
     &       ,Qs(Npoints)                                               &
                                  !  Saturated humidity mixing ratio
     &       ,AEROSOL(Npoints)    ! IN Aerosol mixing ratio(ug/kg)
      LOGICAL                                                           &
     &        L_MURK                    ! IN : Aerosol present
!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
      REAL                                                              &
     &        qt(Npoints)         ! OUT Total water mixing ratio (kg/kg)
!---------------------------------------------------------------------
! constants for visibility calculation used to be set here but now
! set in COMDECK.
!---------------------------------------------------------------------
!
! Local parameter variables
!
      REAL                                                              &
     &  qt_limit       ! Smallest total water mixing ratio value allowed
      PARAMETER (   qt_limit = 0.0001 )
!
! Local workspace variables
!
       INTEGER                                                          &
     &  Point            !  Loop variable for points

      REAL                                                              &
     &  qL                                                              &
                         !  Liquid water mixing ratio (Kg/Kg).
     &, radius_dry                                                      &
                         !  Dry particle radius for aerosol (m)
     &, radius                                                          &
                         !  Radius of fog droplets (m)
     &, radius_star                                                     &
                         !  Water droplet radius divided by dry radius
     &, radius_act                                                      &
                         !  Activation droplet radius
     &, radius_star_act                                                 &
                         !  Activation droplet radius divided by the dry
     &, radius_star_used                                                &
                         !  Water droplet radius divided by the dry
                         !   radius actually used for the relative
                         !   humidity calculation
     &, RH                                                              &
                         !  Relative humidity derived from visibility
     &, A                                                               &
                         !  A0 divided by the dry radius
     &, m_over_m0                                                       &
                         !  Ratio of the aerosol mass mixing ratio and
                         !   the standard aerosol mass mixing ratio
     &, N                !  Number density of aerosol particles (/m3)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



!**

      IF (lhook) CALL dr_hook('VISTOQT',zhook_in,zhook_handle)
      Do Point = 1 , Npoints

!---------------------------------------------------------------------
!* 1.  Calculate the ratio of the aerosol mass mixing ratio to the
!*     standard mass mixing ratio, m_over_m0, and the aerosol number
!*     density, N:
!*
!*                      (1-3p)
!*                  (m )
!*           N = N0 (--)
!*                  (m0)
!*
!*     And the dry radius, radius_dry:
!*                      p
!*                  (m )
!*           r = r0 (--)
!*            d     (m0)
!*
!*     And A (A0 divided by the dry radius).
!*
!*     And the activation radius:
!*
!*                             1/2
!*                  (       3 )
!*                  ( 3 B0 r  )
!*           r    = ( ------d-)
!*            act   (   A0    )
!*
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

!----------------------------------------------------------------------
!* 2.  Calculate a water droplet radius, from the visibility:
!*
!*                                1/2
!*               ( ln( epsilon ) )
!*           r = (---------------)
!*               (  Vis N Beta0  )
!*
!*    (An extra term RecipVisAir is included in the recipical of
!*     visibility to limit visibilities to 100km in clean air).
!----------------------------------------------------------------------

              radius =(VisFactor/N) *                                   &
     &           (1.0/Visibility(Point) - RecipVisAir)
              IF (radius  >   0.0) THEN
                radius = SQRT( radius )
              ELSE
                radius = radius_dry
              ENDIF

!----------------------------------------------------------------------
!* 3.  Provided the diagnosed radius is greater than the dry radius,
!*     calculate the normalised droplet radius, and the saturated
!*     humidity mixing ratio.
!----------------------------------------------------------------------

        If ( radius  >   radius_dry ) then

          radius_star = radius / radius_dry

!----------------------------------------------------------------------
!* 5.  Calculate the corresponding liquid water mixing ratio:
!*
!*
!*               4            (  3     3 )
!*          qL = - Pi rho_w N ( r  - r   )
!*               3            (       d  )
!*
!----------------------------------------------------------------------

          qL = FourThirds * Pi * rho_water * N  *                       &
     &         ( radius**3 - radius_dry**3 )

!----------------------------------------------------------------------
!* 6.  Calculate the relative humidity:
!*
!*                  ( A        B0   )
!*          RH = exp( --  -  ------ )
!*                  ( r       3     )
!*                  (  *     r  - 1 )
!*                  (         *     )
!*
!----------------------------------------------------------------------

          If ( radius_star  <   radius_star_act ) then
            RH = EXP( A/radius_star                                     &
     &                - B0 /( radius_star **3 - 1.0 ) )
          Else
            RH = EXP( A/radius_star_act                                 &
     &                  - B0 /( radius_star_act **3 - 1.0 ) )
          Endif

!----------------------------------------------------------------------
!* 7.  Calculate the total water mixing ratio:  qt = RH * qs(T) + qL
!----------------------------------------------------------------------

          qt(Point) = MAX( RH * qs(Point) + qL, qt_limit )

!----------------------------------------------------------------------
!* 8. If the droplet radius is less than the dry radius, then set the
!*    total water mixing ratio to the minimum value.
!----------------------------------------------------------------------

        Else

          qt(Point) = qt_limit

        End if


      End Do

      IF (lhook) CALL dr_hook('VISTOQT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE VISTOQT
