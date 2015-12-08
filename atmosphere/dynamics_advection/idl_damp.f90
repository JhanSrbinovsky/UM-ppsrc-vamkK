! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_Damp
!

      Subroutine IDL_Damp(                                              &
     &             row_length, rows, n_rows, model_levels               &
     &,            off_x, off_y                                         &
     &,            halo_i, halo_j                                       &
     &,            timestep, height_domain                              &
     &,            eta_theta_levels,  eta_rho_levels                    &
     &,            DMPTIM, HDMP, ZDMP                                   &
     &,            u, v, theta, q                                       &
     &,            u_ref, v_ref, theta_ref, q_ref                       &
     &,            R_u, R_v, theta_star, q_star)



! Purpose: To apply Newtonian relaxation back to initial
!          (reference) values for
!          u, v, q and theta above a certain height,
!          in a similar way to LEM
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None


! Variables with Intent (In)

      Integer                                                           &
     &  row_length                                                      &
                            ! Number of points on a row
     &, rows                                                            &
                            ! Number of rows (p points)
     &, n_rows                                                          &
                            ! Number of rows (v points)
     &, model_levels                                                    &
                            ! Number of model levels
     &, off_x                                                           &
                            ! Size of halo in i
     &, off_y                                                           &
                            ! Size of halo in j.
     &, halo_i                                                          &
                         ! Size of large halo in i
     &, halo_j           ! Size of large halo in j.
      Real                                                              &
     &  timestep                                                        &
                               ! Timestep interval (secs)
     &, height_domain          ! height of the domain

      Real                                                              &
     &  eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

      Real                                                              &
     &  theta (1-off_x:row_length+off_x, 1-off_y:rows+off_y,            &
     &         model_levels)                                            &
     &, q(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,             &
     &      model_levels)                                               &
     &, u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                 &
     &      model_levels)                                               &
     &, v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,               &
     &      model_levels)

! Reference initial values of u, v, theta, q
      Real                                                              &
     &  u_ref(model_levels)                                             &
     &, v_ref(model_levels)                                             &
     &, theta_ref(model_levels)                                         &
     &, q_ref(model_levels)

      Real                                                              &
     & DMPTIM, HDMP, ZDMP   ! Damping layer values

! Latest values of variables
      Real, Intent (InOut) ::                                           &
     & theta_star (1-off_x:row_length+off_x, 1-off_y:rows+off_y,        &
     &         model_levels)                                            &
     &,q_star(1-off_x:row_length+off_x, 1-off_y:rows+off_y,             &
     &      model_levels)                                               &
     &,R_u(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
     &        model_levels)                                             &
     &,R_v(1-off_x:row_length+off_x, 1-off_y:n_rows+off_y,              &
     &        model_levels)

! Local variables

      Integer                                                           &
     &  i,j,k                                                           &
                    ! Loop counters
     &, KDMPMIN     ! starting level of damping layer.

      Real                                                              &
     & DMPCO_TQ(model_levels)                                           &
     &,DMPCO_UV(model_levels)                                           &
                            ! Damping coeffs (on theta and uv levs)
     &,z_tq(model_levels)                                               &
     &,z_uv(model_levels) 
                            ! Height of theta and uv levels

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle



!---------------------------------------------------------------------
!      Damping layer coefficients
!      DMPTIM  Reciprocal of damping time scale (s^-1)
!      HDMP    Height scale of damping layer (m)
!      ZDMP    Starting height of damping layer (m)

      ! Set up damping coefficients and height variables
      ! If L_damp
     
     IF (lhook) CALL dr_hook('IDL_DAMP',zhook_in,zhook_handle)
      Do k = 1, model_levels
        z_tq(k) = eta_theta_levels(k)*height_domain
        z_uv(k) = eta_rho_levels(k)*height_domain
        If (z_uv(k) <= ZDMP) Then
          KDMPMIN=k
          DMPCO_TQ(k)=0.0
          DMPCO_UV(k)=0.0
        Else
          DMPCO_TQ(K)=DMPTIM*( EXP((Z_TQ(K)-ZDMP)/HDMP)-1.0 )
          DMPCO_UV(K)=DMPTIM*( EXP((Z_UV(K)-ZDMP)/HDMP)-1.0 )
        End If

      End Do

        !-------------------------------------------------------------
        ! Add relaxation term to latest values
        !-------------------------------------------------------------


      Do k= KDMPMIN, model_levels
        Do j = 1,rows
          Do i = 1,row_length
            R_u(i,j,k)= R_u(i,j,k)                                      &
     &        - DMPCO_UV(k)*(u(i,j,k)-u_ref(k)) * timestep
            theta_star(i,j,k)= theta_star(i,j,k)                        &
     &        - DMPCO_TQ(k)*(theta(i,j,k)-theta_ref(k)) * timestep
            q_star(i,j,k)= q_star(i,j,k)                                &
     &        - DMPCO_TQ(k)*(q(i,j,k)-q_ref(k)) * timestep
          End Do
        End Do
        Do j = 1,n_rows
          Do i = 1,row_length
            R_v(i,j,k)= R_v(i,j,k)                                      &
     &       - DMPCO_UV(k)*(v(i,j,k)-v_ref(k)) * timestep
          End Do
        End Do

      End Do



! End of routine.
      IF (lhook) CALL dr_hook('IDL_DAMP',zhook_out,zhook_handle)
      RETURN
      End Subroutine IDL_Damp

