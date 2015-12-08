! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_FORCE_SUBS
!

      Subroutine IDL_Force_subs(                                        &
     &             row_length, rows, model_levels                       &
     &,            off_x, off_y                                         &
     &,            halo_i, halo_j                                       &
     &,            timestep                                             &
     &,            r_theta_levels, r_rho_levels                         &
     &,            L_bomex                                              &
     &,            q_star, theta_star                                   &
     &,            q, theta                                             &
     &,            th_inc_subs, q_inc_subs                              &
     &,            th_inc_ls, q_inc_ls )


! Purpose: To apply a large scale vertical velocity forcing w_ls
!          and large scale pot. temperature and moisture forcing
!
! Method:  Calculate subsidence and theta and q forcing.
!          Update large-scale forcing.             
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE PrintStatus_mod
      IMPLICIT NONE

! Variables with Intent (In)

      Integer, Intent(In) ::                                            &
     &  row_length                                                      &
                            ! Number of points on a row
     &, rows                                                            &
                            ! Number of rows
     &, model_levels                                                    &
                            ! Number of model levels
     &, off_x                                                           &
                            ! Size of (small) halo in i
     &, off_y                                                           &
                            ! Size of halo in j.
     &, halo_i                                                          &
     &, halo_j

      Real, Intent(In) ::                                               &
     &  timestep            ! Length of timestep in seconds

      Real, Intent(In) ::                                               &
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, model_levels)

      Logical, Intent(In) :: L_bomex
     

! Variables with Intent (InOut)
      Real, Intent (InOut) ::                                           &
     & q_inc_subs(row_length, rows, model_levels)                       &
     &,th_inc_subs(row_length, rows, model_levels)                      &
                                  ! subsidence increments
     &,q_inc_ls(row_length, rows, model_levels)                         &
     &,th_inc_ls(row_length, rows, model_levels)                      
                                  ! large scale increments

      Real, Intent(InOut) ::                                            &
     & q_star(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &         model_levels)                                            &
                                 ! latest q
     &,q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,              &
     &      model_levels)

      Real, Intent(InOut) ::                                            &
     & theta_star(1-off_x:row_length+off_x, 1-off_y:rows+off_y,         &
     &         model_levels)                                            &
                                 ! latest theta
     &,theta(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
     &         model_levels)     ! initial time theta


      Real                                                              &
     &  w_ls(model_levels)                                              &
                            ! large scale vertical velocity
     &, sthls(model_levels)                                             &
                             ! theta forcing 
     &, sqls(model_levels) ! q forcing

! Local variables
      Real                                                              &
     &  z_uv, z_tq                                                      &
     &, rsecsinday                                                      
                             ! reciprocal seconds in a day      
                             
      Real, Parameter :: secsinday = 86400.  ! second in a day
 
      Integer                                                           &
     &  i,j,k      ! loop counters

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



!---------------------------------------------------------------------
      IF (lhook) CALL dr_hook('IDL_FORCE_SUBS',zhook_in,zhook_handle)
      rsecsinday = 1.0/secsinday
 
      If (L_bomex) Then 
        If ( PrintStatus >= PrStatus_Normal ) Then
          Write (6,*)' BOMEX forcing on '
        End If

        ! Calculate subsidence and theta and q forcing
        Do k = 1, model_levels
          ! Height of uv levels
          z_uv= r_rho_levels(1,1,k)   -   r_theta_levels(1,1,0)
          ! Height or tq levels
          z_tq= r_theta_levels(1,1,k) -   r_theta_levels(1,1,0)

     
          ! BOMEX subsidence on uv levels
          If (z_uv <= 1500.0) Then
            w_ls(k)=-0.0065*z_uv/1500.0
          Else If ( (z_uv>1500.0) .AND. (z_uv<=2100.0) )Then
            w_ls(k)=-0.0065+0.0065*(z_uv-1500.0)/(2100.-1500.)
          Else
            w_ls(k)=0.0
          End If


          ! BOMEX temperature forcing
          If ( z_tq <= 1500.0) Then
            STHLS(K)=-2.0*RSECSINDAY
          Else If ( (z_tq > 1500.0) .AND. (z_tq <= 2500.0) ) Then
            STHLS(K)=(-2.0+(z_tq-1500.0)*2.0/1000.0)                     &
     &                *RSECSINDAY
          Else
            STHLS(K)=0.0
          End If


          ! BOMEX drying forcing
          If(z_tq <= 300.0)Then
            SQLS(K)=-1.2E-08
          Else If ((z_tq > 300.0) .AND. (z_tq <= 500.0)) Then
            SQLS(K)=-1.2E-08+1.2E-08*(z_tq-300.)/(500.-300.)
          Else
            SQLS(K)=0.0
          End If    
         
        End Do ! loop over k

        ! BOMEX subsidence forcing

        Do k = 2, model_levels
          Do j = 1,rows
            Do i = 1,row_length

            ! subsidence forcing
            theta_star(i,j,k) = theta_star(i,j,k) -                     &
     &        timestep*w_ls(k)*                                         &
     &        ( theta(i,j,k+1)  - theta(i,j,k) )/                       &
     &        ( r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k) )

            q_star(i,j,k) = q_star(i,j,k) -                             &
     &        timestep*w_ls(k)*                                         &
     &        ( q(i,j,k+1)  - q(i,j,k) )/                               &
     &        ( r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k) )


            End Do
          End Do
        End Do  ! loop over levels

        Do k = 1, model_levels
          Do j = 1,rows
            Do i = 1,row_length
              th_inc_subs(i,j,k)=theta_star(i,j,k) - th_inc_subs(i,j,k)
              q_inc_subs(i,j,k)= q_star(i,j,k) - q_inc_subs(i,j,k)
            End Do
          End Do
        End Do  ! loop over levels

        Do k = 1, model_levels
          Do j = 1,rows
            Do i = 1,row_length

              ! Large scale forcings
              theta_star(i,j,k) = theta_star(i,j,k) + timestep*STHLS(K)
              q_star(i,j,k)     = q_star(i,j,k) + timestep*SQLS(K)
              th_inc_ls(i,j,k)= timestep*STHLS(K)
              q_inc_ls(i,j,k)=  timestep*SQLS(K)
            End Do
          End Do
        End Do  ! loop over levels

      End If ! L_bomex

! End of routine.
      IF (lhook) CALL dr_hook('IDL_FORCE_SUBS',zhook_out,zhook_handle)
      RETURN
      End Subroutine IDL_Force_subs

