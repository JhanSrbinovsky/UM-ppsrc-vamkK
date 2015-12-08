! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! subroutine gw_wake to deposit blocked flow stress uniformly
! from ground up to sub-grid mountain tops
!
      SUBROUTINE gw_wake(                                               &
     &   s_x_stress,s_y_stress,levels                                   &
     &   ,rho,r_rho_levels,r_theta_levels                               &
     &   ,points,k_top,l_drag,du_dt,dv_dt                               &
! diagnostics
     &  ,stress_ud,points_stress_ud,stress_ud_on, stress_ud_p_on        &
     &  ,stress_vd,points_stress_vd,stress_vd_on                        &
     &  ,stress_ud_wake,points_stress_ud_wake,stress_ud_wake_on         &
     &  ,stress_vd_wake,points_stress_vd_wake,stress_vd_wake_on         &
     &  ,du_dt_wake,points_du_dt_wake,du_dt_wake_on                     &
     &  ,dv_dt_wake,points_dv_dt_wake,dv_dt_wake_on  )



      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

! Description: deposit blocked flow stress uniformly between
!              the surface and the blocked layer top (k_top)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Gravity Wave Drag
!
! code description:
! language: fortran 77 + common extensions
! this code is written to umdp3 v6 programming standards.
! suitable for single column use,rotated grids

! local constants
! none

!
! SUBROUTINE ARGUMENTS
!
      INTEGER                                                           &
                           !,intent(in):
     & levels              ! number of model levels

      INTEGER                                                           &
                           !,intent(in):
     & points              ! number of GWD points

!
! integers below set diagnostic array sizes to points if diagnostic
! is called or to 1 if not. These are set in GWD_CTL2
!
      INTEGER                                                           &
                           !,intent(in):
     & points_stress_ud                                                 &
     &,points_stress_vd                                                 &
     &,points_stress_ud_wake                                            &
     &,points_stress_vd_wake                                            &
     &,points_du_dt_wake                                                &
     &,points_dv_dt_wake

      INTEGER                                                           &
                           !,intent(in):
     & k_top(points)       ! level of blocked layer top

      REAL                                                              &
                           !,intent(in):
     & s_x_stress(points)                                               &
                           ! surface x_stress
     &,s_y_stress(points)  ! surface y_stress

      REAL                                                              &
                           !,intent(in):
     & r_rho_levels(points,levels)
!                          ! heights above z=0 on rho levels

      REAL                                                              &
                           !,intent(in):
     & r_theta_levels(points,levels)
!                          ! heights above z=0 on theta levels

      REAL                                                              &
                           !,intent(in):
     & rho(points,levels)  ! density on rho levels

      REAL                                                              &
                           !,intent(inout):
     & du_dt(points,levels)                                             &
                           ! total GWD du/dt
     &,dv_dt(points,levels)! total GWD dv/dt

!
! Diagnostic arrays
!
      REAL                                                              &
     & du_dt_wake(points_du_dt_wake,levels)                             &
                                            ! u-acceln  diagnostic
     &,dv_dt_wake(points_dv_dt_wake,levels)                             &
                                            ! v-acceln  diagnostic
     &,stress_ud(points_stress_ud,0:levels)                             &
                                            ! u stress  diagnostic
     &,stress_vd(points_stress_vd,0:levels)                             &
                                            ! v stress  diagnostic
     &,stress_ud_wake(points_stress_ud_wake,0:levels)                   &
!                                           ! u blocked flow stress diag
     &,stress_vd_wake(points_stress_vd_wake,0:levels)
!                                           ! v blocked flow stress diag


      LOGICAL                                                           &
                           !,intent(in):
     & l_drag(points)      ! true if a non-zero surface stress

!
! switches for diagnostics.
!
      LOGICAL                                                           &
                           !,intent(in):
     & stress_ud_on                                                     &
     &,stress_ud_p_on                                                   &
     &,stress_vd_on                                                     &
     &,stress_ud_wake_on                                                &
     &,stress_vd_wake_on                                                &
     &,du_dt_wake_on                                                    &
     &,dv_dt_wake_on


!----------------------------------------------------------------
! LOCAL ARRAYS AND SCALARS
!----------------------------------------------------------------
      REAL                                                              &
     & x_stress(points,2)                                               &
                            ! x_stresses (layer boundaries)
     &,y_stress(points,2)   ! y_stresses (layer boundaries)

      REAL                                                              &
     & dz_x_stress(points)                                              &
                            ! x component of stress gradient
     &,dz_y_stress(points)  ! y component of stress gradient

      REAL                                                              &
     & delta_z              ! difference in height across layer(s)

      INTEGER                                                           &
     & i,k                  ! loop counters

      INTEGER                                                           &
     & kk,kl,ku             ! level counters in routine

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



! function and subroutine calls
! none

!-------------------------------------------------------------------
!   1. start level  preliminaries
!-------------------------------------------------------------------

      IF (lhook) CALL dr_hook('GW_WAKE',zhook_in,zhook_handle)
      kl=1
      ku=2

      Do i=1,points
        If( l_drag(i) ) Then

          delta_z = r_theta_levels(i,k_top(i))

          dz_x_stress(i) = s_x_stress(i) / delta_z
          dz_y_stress(i) = s_y_stress(i) / delta_z

        End if ! if l_drag

      End do  ! points


      If( stress_ud_on .or. stress_ud_p_on .or. stress_ud_wake_on) Then

        Do i=1,points
          If( l_drag(i) ) Then
            x_stress(i,kl) = s_x_stress(i)
          End if
        End do

        If( stress_ud_on .or. stress_ud_p_on) Then
          Do i=1,points
            If( l_drag(i) ) Then
              stress_ud(i,0) = stress_ud(i,0) + x_stress(i,kl)
            End if
          End do
        End if

        If( stress_ud_wake_on ) Then
          Do i=1,points
            If( l_drag(i) ) Then
              stress_ud_wake(i,0) = x_stress(i,kl)
            End if
          End do
        End if

      End if  !  stress_ud_on .or. stress_ud_p_on .or. stress_ud_wake_on


      If( stress_vd_on .or. stress_vd_wake_on) Then

        Do i=1,points
          If( l_drag(i) ) Then
            y_stress(i,kl) = s_y_stress(i)
          End if
        End do

        If( stress_vd_on ) Then
          Do i=1,points
            If( l_drag(i) ) Then
              stress_vd(i,0) = stress_vd(i,0) + y_stress(i,kl)
            End if
          End do
        End if

        If( stress_vd_wake_on ) Then
          Do i=1,points
            If( l_drag(i) ) Then
              stress_vd_wake(i,0) = y_stress(i,kl)
            End if
          End do
        End if

      End if  !  stress_vd_on .or. stress_vd_wake_on

!------------------------------------------------------------------
!    2 loop levels
!      k is the rho level counter
!------------------------------------------------------------------

      Do k=1,levels


        Do i=1,points
          If( l_drag(i) .and. k <= k_top(i) ) Then
!
!  note that a constant stress drop across a level in height coords =>
!  a non-uniform drag because of rho variation.
!
            du_dt(i,k) =  du_dt(i,k) - dz_x_stress(i)/rho(i,k)
            dv_dt(i,k) =  dv_dt(i,k) - dz_y_stress(i)/rho(i,k)

          End if   ! if l_drag(i) .and. k<=k_top(i)

        End do ! points

! diagnostics
        If( du_dt_wake_on ) Then
          Do i=1,points
            If( l_drag(i) .and. k <= k_top(i) ) Then
              du_dt_wake(i,k) = - dz_x_stress(i)/rho(i,k)
            End if
          End do
        End if

        If( dv_dt_wake_on ) Then
          Do i=1,points
            If( l_drag(i) .and. k <= k_top(i) ) Then
              dv_dt_wake(i,k) = - dz_y_stress(i)/rho(i,k)
            End if
          End do
        End if

        If( stress_ud_on .or. stress_ud_p_on .or. stress_ud_wake_on) Then
          Do i=1,points
            If( l_drag(i) .and. k <= k_top(i) ) Then
              If ( k  ==  1 ) Then
                delta_z = r_theta_levels(i,k)
              Else
                delta_z= r_theta_levels(i,k) - r_theta_levels(i,k-1)
              End if
              x_stress(i,ku) = x_stress(i,kl)
              x_stress(i,ku) = x_stress(i,kl)-dz_x_stress(i)*delta_z
            End if   ! l_drag & k < k_top
          End do

          If( stress_ud_on .or. stress_ud_p_on) Then
            Do i=1,points
              If( l_drag(i) .and. k <= k_top(i) ) Then
                stress_ud(i,k) = stress_ud(i,k) + x_stress(i,ku)
              End if   ! l_drag & k < k_top
            End do
          End if       ! stress_ud on

          If( stress_ud_wake_on ) Then
            Do i=1,points
              If( l_drag(i) .and. k <= k_top(i) ) Then
                stress_ud_wake(i,k) = x_stress(i,ku)
              End if   ! l_drag & k < k_top
            End do
          End if       ! stress_ud_wake on

        End if      ! stress_ud_on .or. stress_ud_p_on .or. stress_wake_ud_on


        If( stress_vd_on .or. stress_vd_wake_on) Then
          Do i=1,points
            If( l_drag(i) .and. k <= k_top(i) ) Then
              If ( k  ==  1 ) Then
                delta_z = r_theta_levels(i,k)
              Else
                delta_z= r_theta_levels(i,k) - r_theta_levels(i,k-1)
              End if
              y_stress(i,ku) = y_stress(i,kl)
              y_stress(i,ku) = y_stress(i,kl)-dz_y_stress(i)*delta_z
            End if   ! l_drag & k < k_top
          End do

          If( stress_vd_on ) Then
            Do i=1,points
              If( l_drag(i) .and. k <= k_top(i) ) Then
                stress_vd(i,k) = stress_vd(i,k) + y_stress(i,ku)
              End if   ! l_drag & k < k_top
            End do
          End if       ! stress_vd on

          If( stress_vd_wake_on ) Then
            Do i=1,points
              If( l_drag(i) .and. k <= k_top(i) ) Then
                stress_vd_wake(i,k) = y_stress(i,ku)
              End if   ! l_drag & k < k_top
            End do
          End if       ! stress_vd_wake on

        End if         ! stress_vd_on .or. stress_wake_vd_on


! swap storage for lower and upper layers
        kk=kl
        kl=ku
        ku=kk

      End do
!   end loop levels

      IF (lhook) CALL dr_hook('GW_WAKE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE gw_wake
