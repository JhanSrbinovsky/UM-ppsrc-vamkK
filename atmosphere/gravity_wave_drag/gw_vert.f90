! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine gw_vert to calculate the vertical profile of gravity wave
!            stress and associated wind increments
!
      SUBROUTINE gw_vert(                                               &
     &  rho,r_rho_levels,r_theta_levels                                 &
     & ,theta,u,v,levels,points,kay,sd_orog                             &
     & ,s_x_lin_stress,s_y_lin_stress,s_x_wake_stress,s_y_wake_stress   &
     & ,s_x_orog,s_y_orog,du_dt,dv_dt                                   &
     & ,k_top,k_top_max,lift,l_drag,fr,rho_s,l_fix_gwsatn,l_gwd_40km    &
     & ,sat_scheme,fsat                                                 &
! diagnostics
     &  ,stress_ud,points_stress_ud,stress_ud_on,stress_ud_p_on         &
     &  ,stress_vd,points_stress_vd,stress_vd_on                        &
     &  ,stress_ud_satn,points_stress_ud_satn,stress_ud_satn_on         &
     &  ,stress_vd_satn,points_stress_vd_satn,stress_vd_satn_on         &
     &  ,stress_ud_wake,points_stress_ud_wake,stress_ud_wake_on         &
     &  ,stress_vd_wake,points_stress_vd_wake,stress_vd_wake_on         &
     &  ,du_dt_satn,points_du_dt_satn,du_dt_satn_on,du_dt_satn_p_on     &
     &  ,dv_dt_satn,points_dv_dt_satn,dv_dt_satn_on                     &
     &  ,du_dt_wake,points_du_dt_wake,du_dt_wake_on                     &
     &  ,dv_dt_wake,points_dv_dt_wake,dv_dt_wake_on )


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

! Description:
!     calculates the gwd stress profiles and wind increments.
!     1. calculate stress profile and wind increments for
!        linear hydrostatic waves.
!     2. calculate stress profile and wind increments for
!        the blocked flow.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Gravity Wave Drag
!
! code description:
! language: fortran 90
! this code is written to umdp3 v6 programming standards.
! suitable for single column use,rotated grids


! local constants
! none

!
! SUBROUTINE ARGUMENTS:
!

      INTEGER                                                           &
                             !,intent(in):
     & levels                ! number of model levels

      INTEGER                                                           &
                             !,intent(in):
     & points                ! number of points

      INTEGER                                                           &
                             !,intent(in):
     & k_top(points)                                                    &
                             ! model level at mountain tops
!                            ! full definition in gwsurf
     &,k_top_max             ! max(k_top)

!
! Integers below determine size of diagnostic arrays. They are
! set to points if called and to 1 if not. This is done in GWD_CTL2
!
      INTEGER                                                           &
                             !,intent(in):
     & points_stress_ud                                                 &
     &,points_stress_vd                                                 &
     &,points_stress_ud_satn                                            &
     &,points_stress_vd_satn                                            &
     &,points_stress_ud_wake                                            &
     &,points_stress_vd_wake                                            &
     &,points_du_dt_satn                                                &
     &,points_dv_dt_satn                                                &
     &,points_du_dt_wake                                                &
     &,points_dv_dt_wake

      REAL                                                              &
                             !,intent(in):
     & r_rho_levels(points,levels)
!                            ! heights on rho levels

      REAL                                                              &
                             !,intent(in):
     & r_theta_levels(points,levels)
!                            ! heights on theta levels

      REAL                                                              &
                             !,intent(in):
     & rho(points,levels)
!                            ! density on rho levels

      REAL                                                              &
                             !,intent(in):
     & theta(points,levels)  ! theta field

      REAL                                                              &
                             !,intent(in):
     & u(points,levels)                                                 &
                             ! u field
     &,v(points,levels)      ! v field

      REAL                                                              &
                             !,intent(in):
     & s_x_lin_stress(points)                                           & 
                             ! 'surface' lin  x_stress
     &,s_y_lin_stress(points)! 'surface' lin  y_stress

      REAL                                                              &
                             !,intent(in):
     & s_x_wake_stress(points)                                          &
!                            ! 'surface' x_wake_stress
     &,s_y_wake_stress(points)
!                            ! 'surface' y_wake_stress

      REAL                                                              &
                             !,intent(in):
     & s_x_orog(points)                                                 & 
                             ! 'surface' x_stress term
     &,s_y_orog(points)      ! 'surface' y_stress term

      REAL                                                              &
                             !,intent(in):
     & lift(points)          ! blocked layer depth


      REAL                                                              &
                              !,intent(in):
     & sd_orog(points)        ! standard deviation of the sub-grid
!                             ! orography

      REAL                                                              &
                             !,intent(in):
     & fr(points)            ! low level froude number

      REAL                                                              &
                             !,intent(in):
     & kay                   ! GWD hydrostatic constant

      REAL                                                              &
                             !,intent(in):
     & rho_s(points)         ! Low level density calculated in gwsurf

      REAL                                                              &
                             !,intent(out):
     & du_dt(points,levels)                                             &
                             ! total GWD du/dt
     &,dv_dt(points,levels)  ! total GWD dv/dt

!
! Diagnostics
!
      REAL                                                              &
                             !,intent(out):
     & stress_ud     (points_stress_ud,0:levels)                        &
!                            ! x total stress diag
     &,stress_vd     (points_stress_vd,0:levels)
!                            ! y total stress diag

      REAL                                                              &
                             !,intent(out):
     & stress_ud_satn(points_stress_ud_satn,0:levels)                   &
!                            ! x saturation stress diag
     &,stress_vd_satn(points_stress_vd_satn,0:levels)
!                            ! y saturation stress diag

      REAL                                                              &
                             !,intent(out):
     & stress_ud_wake(points_stress_ud_wake,0:levels)                   &
!                            ! x blocked flow stress diag
     &,stress_vd_wake(points_stress_vd_wake,0:levels)
!                            ! y blocked flow stress diag

      REAL                                                              &
                             !,intent(out):
     & du_dt_satn(points_du_dt_satn,levels)                             &
!                            ! du/dt diagnostic (saturation)
     &,dv_dt_satn(points_dv_dt_satn,levels)
!                            ! dv/dt diagnostic (saturation)

      REAL                                                              &
                             !,intent(out):
     & du_dt_wake(points_du_dt_wake,levels)                             &
!                            ! du/dt diagnostic (blocked flow)
     &,dv_dt_wake(points_dv_dt_wake,levels)
!                            ! dv/dt diagnostic (blocked flow)

      REAL                                                              &
                            !,intent(in):
     & fsat                 ! Froude number used to scale critical
!                           ! wave amplitude for breaking

      LOGICAL                                                           &
                             !,intent(in):
     & l_drag(points)        ! true if a non-zero surface stress

      LOGICAL                                                           &
                             !,intent(in):
     & l_fix_gwsatn                                                     &
                             ! switch to include minor bug fixes  
     &,l_gwd_40km            ! switch to turn off GWD above 40km

      INTEGER                                                           &
                            !,intent(in):
     & sat_scheme           ! Switch at vn7.1
!                           ! If =1 use amplitude based saturation test
!                           ! If =0    old stress based saturation test


!
!  Diagnostic switches
!
      LOGICAL                                                           &
                             !,intent(in):
     & stress_ud_on                                                     &
     &,stress_ud_p_on                                                   &
     &,stress_vd_on                                                     &
     &,stress_ud_satn_on                                                &
     &,stress_vd_satn_on                                                &
     &,stress_ud_wake_on                                                &
     &,stress_vd_wake_on                                                &
     &,du_dt_satn_on                                                    &
     &,du_dt_satn_p_on                                                  &
     &,dv_dt_satn_on                                                    &
     &,du_dt_wake_on                                                    &
     &,dv_dt_wake_on

!--------------------------------------------------------------------
! LOCAL ARRAYS AND SCALARS
!--------------------------------------------------------------------

      INTEGER                                                           &
     & i,k                   ! loop counter in routine

      REAL                                                              &
     & unit_x(points)                                                   &
                             ! x_compnt of unit stress vector
     &,unit_y(points)        ! y_compnt of unit stress vector

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! function and subroutine calls:
      EXTERNAL gw_satn,gw_wake

!-------------------------------------------------------------------
!   1.0 start  preliminaries
! initialise increment and increment diagnostics
!------------------------------------------------------------
      IF (lhook) CALL dr_hook('GW_VERT',zhook_in,zhook_handle)
      Do k=1,levels

        Do i=1,points
          du_dt(i,k)=0.0
          dv_dt(i,k)=0.0
        End do

        If( du_dt_satn_on .or. du_dt_satn_p_on ) Then
          Do i=1,points
            du_dt_satn(i,k)=0.0
          End do
        End if

        If( dv_dt_satn_on ) Then
          Do i=1,points
            dv_dt_satn(i,k)=0.0
          End do
        End if

        If( du_dt_wake_on ) Then
          Do i=1,points
            du_dt_wake(i,k)=0.0
          End do
        End if

        If( dv_dt_wake_on ) Then
          Do i=1,points
            dv_dt_wake(i,k)=0.0
          End do
        End if

      End do ! levels

!
!  and now for stress diagnostics
!
      Do k=0,levels

        If (stress_ud_on .or. stress_ud_p_on ) Then
          Do i=1,points
            stress_ud(i,k) = 0.0
          End do
        End if

        If (stress_vd_on ) Then
          Do i=1,points
            stress_vd(i,k) = 0.0
          End do
        End if

        If (stress_ud_satn_on ) Then
          Do i=1,points
            stress_ud_satn(i,k) = 0.0
          End do
        End if

        If (stress_vd_satn_on ) Then
          Do i=1,points
            stress_vd_satn(i,k) = 0.0
          End do
        End if

        If (stress_ud_wake_on ) Then
          Do i=1,points
            stress_ud_wake(i,k) = 0.0
          End do
        End if

        If (stress_vd_wake_on ) Then
          Do i=1,points
            stress_vd_wake(i,k) = 0.0
          End do
        End if

      End do

!---------------------------------------------------------------------
!  2. launch linear hydrostatic waves from level k_top and
!     calculate stress profile due to wave saturation effects.
!---------------------------------------------------------------------
! DEPENDS ON: gw_satn
      CALL gw_satn(                                                     &
     &   rho,r_rho_levels,r_theta_levels                                &
     &  ,theta,u,v,s_x_lin_stress,s_y_lin_stress,levels                 &
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

!
!------------------------------------------------------------------
! 3 apply uniform stress reduction between surface and k_top
!     for points where blocked flow was diagnosed.
!------------------------------------------------------------------

! DEPENDS ON: gw_wake
       CALL gw_wake(                                                    &
     &   s_x_wake_stress,s_y_wake_stress,levels                         &
     &   ,rho,r_rho_levels,r_theta_levels                               &
     &   ,points,k_top,l_drag,du_dt,dv_dt                               &
! diagnostics
     &  ,stress_ud,points_stress_ud,stress_ud_on,stress_ud_p_on         &
     &  ,stress_vd,points_stress_vd,stress_vd_on                        &
     &  ,stress_ud_wake,points_stress_ud_wake,stress_ud_wake_on         &
     &  ,stress_vd_wake,points_stress_vd_wake,stress_vd_wake_on         &
     &  ,du_dt_wake,points_du_dt_wake,du_dt_wake_on                     &
     &  ,dv_dt_wake,points_dv_dt_wake,dv_dt_wake_on )


      IF (lhook) CALL dr_hook('GW_VERT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE gw_vert

