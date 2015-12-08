! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates fraction of grid-box covered with fog.
!
! Subroutine Interface:
      SUBROUTINE FOG_FR(                                                &
     & p_layer,RHCRIT,LEVELS,PFIELD,                                    &
     & T,AEROSOL,L_MURK,Q,QCL,QCF,VIS,FF,NVIS)

! Modules 
      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE
!
! Description:   
!   Calculates the fraction of a gridsquare with visibility
!   less than threshold, Vis_thresh, given the total water
!   mixing ratio (qt), temperature (T), pressure (p), and the
!   (Murk) aerosol mass mixing ratio (m), assuming a triangular
!   distribution of states about the median, characterised by
!   a critical relative humdity value, RHcrit.
!   NB:  Throughout, levels are counted from the bottom up,
!   i.e. the lowest level under consideration is level 1, the
!   next lowest level 2, and so on.
!
!   Suitable for single-column use.
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
     & LEVELS                                                           &
                           ! IN No. of levels being processed.
     &,PFIELD                                                           &
                           ! IN No. of points in field (at one level).
     &,NVIS                ! IN No. of visibility thresholds
      REAL                                                              &
     & p_layer(PFIELD,LEVELS)                                           &
                              ! IN pressure (Pa) at processed levels.
     &,RHCRIT(LEVELS)                                                   &
                           ! IN Critical relative humidity.  See the
!                          !    the paragraph incorporating eqs P292.11
!                          !    to P292.14; the values need to be tuned
!                          !    for the given set of levels.
     &,Q(PFIELD,LEVELS)                                                 &
                           ! IN Specific Humidity
!                          !    (kg per kg air).
     &,QCL(PFIELD,LEVELS)                                               &
                           ! Cloud liquid water content at
!                          !     processed levels (kg per kg air).
     &,QCF(PFIELD,LEVELS)                                               &
                           ! Cloud ice content at processed levels
!                          !    (kg per kg air).
     &,T(PFIELD,LEVELS)                                                 &
                           ! IN Temperature (K).
     &,AEROSOL(PFIELD,LEVELS)                                           &
                              ! IN Aerosol mixing ratio(ug/kg)
     &,VIS(PFIELD,LEVELS,NVIS)  ! Visibility thresholds
      LOGICAL                                                           &
     &   L_MURK               ! IN : Aerosol present

      REAL                                                              &
     & FF(PFIELD,LEVELS,NVIS)   ! OUT Vis prob at processed levels
!                          !     (decimal fraction).
!
!*--------------------------------------------------------------------
!*L  Workspace usage----------------------------------------------------
      REAL                                                              &
                           ! "Automatic" arrays on Cray.
     & P(PFIELD)                                                        &
     &,QT(PFIELD)                                                       &
                           ! total of cloud water and vapour
     &,QS(PFIELD)                                                       &
                           ! Saturated spec humidity for temp T
     &,qt_thresh(PFIELD)                                                &
                           ! modified qt
     &,bs
!*L  External subroutine called ----------------------------------------
      EXTERNAL QSAT_WAT,VISTOQT
!* Local, including SAVE'd, storage------------------------------------
!
      INTEGER K,I,J     ! Loop counters: K - vertical level index.
!                       !                I - horizontal field index.
                        !                J - Vis threshold index.

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
!
!-----------------------------------------------------------------------
!L Subroutine structure :
!L Loop round levels to be processed.
!-----------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook('FOG_FR',zhook_in,zhook_handle)
      DO K=1,LEVELS
!
!-----------------------------------------------------------------------
!L 1. Calculate Pressure and initialise temporary arrays
!-----------------------------------------------------------------------
!
        DO I=1,PFIELD
          P(I)=p_layer(I,K)
          QT(I)=Q(I,K)+QCL(I,K)
        ENDDO ! Loop over points

!-----------------------------------------------------------------------
!* 2.  Calculate total water threshold corresponding to visibility
!      Since Qs is needed more than once, pre-calculate and pass it
!-----------------------------------------------------------------------

! DEPENDS ON: qsat_wat
        CALL QSAT_WAT (Qs,T(1,K),P,PFIELD)

        DO J=1,NVIS

! DEPENDS ON: vistoqt
          Call VISTOQT( VIS(1,K,J), Qs, AEROSOL(1,K), L_MURK,           &
     &                PFIELD, qt_thresh )


!-----------------------------------------------------------------------
!* 3.  Calculate the width of the distribution in total water space, bs:
!*
!*           bs = ( 1 - RHcrit ) * qs(T)
!*
!-----------------------------------------------------------------------

          Do I = 1 , PFIELD

            bs = (1.0-RHcrit(K)) * qs(I)

!=======================================================================
!* 4.  Calculate the fraction of states in a triangular
!*     distribution which exceed the total water threshold.
!=======================================================================

!-----------------------------------------------------------------------
!* 4.1 If total water threshold value is less than the total water value
!*     minus the width of the distribution, then all of the states have
!*     a total water value exceeding the threshold, so set the
!*     visibility fraction to 1.0
!-----------------------------------------------------------------------

            if ( qt_thresh(I)  <=  qt(I)-bs ) then

              FF(I,K,J) = 1.0

!-----------------------------------------------------------------------
!* 4.2 If total water threshold value is greater than the total water
!*     value minus the width of the distribution, but less than the
!*     total water value then the visibility fraction, VF, is given by:
!*
!*                                                    2
!*                             ( qt       - qt + bs  )
!*            VF = 1.0 - 0.5 * (    thresh           )
!*                             ( ------------------- )
!*                             (          bs         )
!*
!-----------------------------------------------------------------------

             Else if ( qt_thresh(I)  >   qt(I)-bs .AND.                 &
     &                 qt_thresh(I)  <=  qt(I) ) then

               FF(I,K,J) = 1.0 - 0.5 *                                  &
     &              (( qt_thresh(I) - qt(I) + bs )/ bs)**2

!-----------------------------------------------------------------------
!* 4.3 If total water threshold value is greater than the total water
!*     value, but less than the total water value plus the width of the
!*     distribution, then the visibility fraction, VF, is given by:
!*
!*                                              2
!*                       ( qt + bs - qt        )
!*            VF = 0.5 * (             thresh  )
!*                       ( ------------------- )
!*                       (          bs         )
!*
!-----------------------------------------------------------------------

             Else if ( qt_thresh(I)  >   qt(I) .AND.                    &
     &                 qt_thresh(I)  <=  qt(I)+bs    ) then

                FF(I,K,J)= 0.5 * (( qt(I) + bs - qt_thresh(I))/bs)**2

!-----------------------------------------------------------------------
!* 4.4 If total water threshold value is greater than the total water
!*     value plus the width of the distribution, then non of the states
!*     have a total water value exceeding the threshold, so set the
!*     visibility fraction to 0.0
!-----------------------------------------------------------------------

             Else

               FF(I,K,J) = 0.0

            End if

          End Do ! Loop over PFIELD I

        End Do ! Loop over VIS J

      ENDDO ! Loop over levels
!
      IF (lhook) CALL dr_hook('FOG_FR',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE FOG_FR
