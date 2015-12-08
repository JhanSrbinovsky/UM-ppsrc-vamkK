! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculate contribution of precipitation to extinction for visibility.
!
      SUBROUTINE VIS_PRECIP                                             &
     &           (Vis_No_Precip                                         &
                                                      !INPUT
     &           ,LCA,CCA,PCT                                           &
                                                      !INPUT
     &           ,Beta_LS_Rain, Beta_LS_Snow                            &
                                                      !INPUT
     &           ,Beta_C_Rain, Beta_C_Snow                              &
                                                      !INPUT
     &           ,P_FIELD,POINTS,K1STPT                                 &
                                                      !INPUT
     &           ,Vis_overall,Vis_LSP,Vis_CP                            &
                                                      !OUTPUT
     &           ,ERROR)                              !OUTPUT

! Modules 
      USE visbty_constants_mod, ONLY : LnLiminalContrast 
      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE
!
! Description: 
!   Process fields of precipitation intensity to give scattering coefft
!   in 1/metres. This is added to an input visibility to give a total.
!   Calculated at a single model level or level within surface layer 
!   e.g. screen height (1.5m)
!
! Documentation:
!   Forecasting Research Scientific Paper NO.4
!   Diagnosis of visibility in the UK Met Office Mesoscale Model
!   and the use of a visibility analysis to constrain initial
!   conditions.  SP Ballard, BJ Wright, BW Golding    1992
!     NIMROD diagnostic:
!   Wright, B. J., 1997: Improvements to the Nimrod Visibility
!     Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!   Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!     for Nimrod. Met. Office FR Tech Rep., No. 222.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Atmosphere Service
!
! Code description:
!   Programming standard:  Unified Model Documentation Paper No 3,
!                          Version 5, dated 08/12/92.
!

!---------------------------------------------------------------------
! input variables-----------------------------------------------------
!---------------------------------------------------------------------
      INTEGER                                                           &
     & P_FIELD                                                          &
                                        ! IN NO. points in field.
     &,POINTS                                                           &
                ! IN Number of gridpoints being processed.
     &,K1STPT                                                           &
                ! IN First gridpoint processed within complete field.
     &,ERROR    ! OUT Error code
      REAL                                                              &
     & Vis_No_Precip(P_FIELD)                                           &
                                        ! IN Vis outside precip.
     &,LCA(P_FIELD)                                                     &
                                        ! IN Total Layer Cloud.
     &,CCA(P_FIELD)                                                     &
                                        ! IN Convective Cloud.
     &,Beta_LS_Rain(P_FIELD)                                            &
                                        ! IN Scattering in LS Rain.
     &,Beta_LS_Snow(P_FIELD)                                            &
                                        ! IN Scattering in LS Snow.
     &,Beta_C_Rain(P_FIELD)                                             &
                                        ! IN Scattering in Conv Rain
     &,Beta_C_Snow(P_FIELD)             ! IN Scattering in Conv Snow

      LOGICAL                                                           &
     & PCT                              ! IN T:Cloud amounts are in %
!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
      REAL                                                              &
     & Vis_overall(P_FIELD)                                             &
                                       ! OUT Visibility overall
     &,Vis_LSP(P_FIELD)                                                 &
                                       ! OUT Visibility in LS Precip.
     &,Vis_CP(P_FIELD)                 ! OUT Visibility in Conv Precip.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
! Local varables:------------------------------------------------------
!  Define local variables ---------------------------------------------
      INTEGER I       ! Loop counters: I - horizontal field index;

      REAL                                                              &
     & Beta_No_Precip                                                   &
     &,P_LSP(P_FIELD)                                                   &
     &,P_CP(P_FIELD)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!---------------------------------------------------------------------
!  External subroutine called ----------------------------------------
!---------------------------------------------------------------------
!
      IF (lhook) CALL dr_hook('VIS_PRECIP',zhook_in,zhook_handle)
      ERROR=0
      IF((K1STPT+POINTS-1) >  P_FIELD)THEN
        ERROR=1
        GOTO 9999
      ENDIF
      IF(PCT) THEN
        DO I=K1STPT,K1STPT+POINTS-1
          P_CP(I)=CCA(I)/100.0
          P_LSP(I)=(1.0-P_CP(I))*LCA(I)/100.0
        ENDDO
      ELSE
        DO I=K1STPT,K1STPT+POINTS-1
          P_CP(I)=CCA(I)
          P_LSP(I)=(1.0-P_CP(I))*LCA(I)
        ENDDO
      ENDIF


      DO I=K1STPT,K1STPT+POINTS-1

        Beta_No_Precip=-LnLiminalContrast/Vis_No_Precip(I)

        IF(P_LSP(I)  >   0.0) THEN
          Vis_LSP(I) = -LnLiminalContrast /                             &
     &      (Beta_No_Precip +                                           &
     &       Beta_LS_Rain(I) + Beta_LS_Snow(I))
        ELSE
          Vis_LSP(I)=Vis_No_Precip(I)
        ENDIF

        IF(P_CP(I)  >   0.0) THEN
          Vis_CP(I) = -LnLiminalContrast /                              &
     &      (Beta_No_Precip +                                           &
     &       Beta_C_Rain(I) + Beta_C_Snow(I))
        ELSE
          Vis_CP(I)=Vis_No_Precip(I)
        ENDIF

! Ensure no rounding problems lead to vis > vis in clear air
        Vis_overall(I) = MIN((1.0-P_CP(I)-P_LSP(I))*Vis_No_Precip(I) +  &
     &                   P_LSP(I)*Vis_LSP(I) +P_CP(I)*Vis_CP(I),        &
     &                   Vis_No_Precip(I))

      ENDDO

 9999 Continue

      IF (lhook) CALL dr_hook('VIS_PRECIP',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE VIS_PRECIP
