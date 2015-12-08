! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine IDL_setup_def_front

      Subroutine IDL_setup_def_front(                                   &
     &                      row_length, rows, n_rows, model_levels      &
     &,                     global_row_length, global_rows              &
     &,                     l_datastart                                 &
     &,                     halo_i, halo_j, off_x, off_y                &
     &,                     delta_lambda, delta_phi                     &
     &,                     lambda_p, phi_p, lambda_u, phi_v            &
     &,                     L_regular                                   &
     &,                     base_phi                                    &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     q, qcl, qcf, u_adv, v_adv, w_adv            &
     &,                     exner_theta_levels                          &
     &,                     f3_at_u                                     &
     &,                     u_ref, t_horizfn_data                       &
     &,                     height_domain                               &
     &,                     L_trivial_trigs                             &
     &                      )

! Purpose:
!          Sets up "balanced" deformation wind field and zonal jet
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + extensions
!   This code is written to UMDP3 programming standards.
!

      USE earth_constants_mod, ONLY: g, earth_radius
      USE atmos_constants_mod, ONLY:                                    &
          cp, repsilon  

      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

      Integer                                                           &
     &  row_length                                                      &
                             ! number of points on a row
     &, rows                                                            &
                             ! number of rows in a theta-field
     &, n_rows                                                          &
                             ! number of rows in a v-field
     &, model_levels                                                    &
                             ! number of model levels
     &, halo_i                                                          &
                             ! Size of halo in i direction.
     &, halo_j                                                          &
                             ! Size of halo in j direction.
     &, off_x                                                           &
     &, off_y                                                           &
     &, global_row_length                                               &
     &, global_rows                                                     &
     &, l_datastart(2)       ! First gridpoints held by this processor

      Real                                                              &
     &  height_domain                                                   &
     &, delta_lambda                                                    &
     &, delta_phi                                                       &
     &, base_phi


      Real                                                              &
           ! VarRes horizontal co-ordinate information
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, phi_p(1-halo_j:rows+halo_j)                                     &
     &, lambda_u(1-halo_i:row_length+halo_i)                            &
     &, phi_v(1-halo_j:n_rows+halo_j)

      Logical, Intent(In) :: L_regular

      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, model_levels)                &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)                                    &
     &, f3_at_u (1-off_x:row_length+off_x,                              &
     &                      1-off_y:rows+off_y)

!Output Arrays from this routine
      Real                                                              &
     &  q(1-halo_i:row_length+halo_i                                    &
     &,       1-halo_j:rows+halo_j, model_levels)                       &
     &, qcl(1-halo_i:row_length+halo_i                                  &
     &,       1-halo_j:rows+halo_j, model_levels)                       &
     &, qcf(1-halo_i:row_length+halo_i                                  &
     &,       1-halo_j:rows+halo_j, model_levels)                       &
     &, u_adv(1-halo_i:row_length+halo_i                                &
     &,       1-halo_j:rows+halo_j, model_levels)                       &
     &, v_adv(1-halo_i:row_length+halo_i                                &
     &,       1-halo_j:n_rows+halo_j, model_levels)                     &
     &, w_adv(1-halo_i:row_length+halo_i                                &
     &,       1-halo_j:rows+halo_j, 0:model_levels)                     &
     &, exner_theta_levels(1-halo_i:row_length+halo_i                   &
     &,                    1-halo_j:rows+halo_j, model_levels)

      Real                                                              &
     &  u_ref(model_levels)

      Real                                                              &
     &  t_horizfn_data(10)  ! Data values describing horizontal t fn

      Logical                                                           &
     &  L_trivial_trigs

! local variables
      Integer                                                           &
     &  i, j, k, gj

      Real                                                              &
     &  weight                                                          &
     &, temp, temp1, temp2  

      Real                                                              &
     &  work1(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, work2(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, work3(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, work4(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)         &
     &, work5(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j)

      Real                                                              &
           ! components of coriolis force.
     &  f3_at_p(row_length, rows)

      Real geop,theta2,theta3,alpha
      Real ujet_j(1-halo_j:global_rows+halo_j)
      Real thetabal(1-halo_j:global_rows+halo_j)

      Integer jetmaxlevel

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!---------------------------------------------------------------------
! Set horizontal scaling functions
!---------------------------------------------------------------------
! Deformation wind field
! Set up deformation flow u=-ax, v=ay
! Define x=0, y=0 centre point of domain and setup a 2D function
! with distance from the centre point

      IF (lhook) CALL dr_hook('IDL_SETUP_DEF_FRONT',zhook_in,zhook_handle)

        ! Set up linear ramp fn for u-wind in E-W direction
        temp1 = l_datastart(1)-1-global_row_length/2.0
        temp2 = l_datastart(1)-1-(global_row_length+1)/2.0

        If (L_trivial_trigs) Then

          Do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i

              ! Calculate E-W distance from grid centre for p,theta cols
              work2(i,j) = (temp1 + i + 0.5) *                          &
     &                        delta_lambda * Earth_radius

              ! Calculate N-S distance from grid centre for p,theta rows
              work3(i,j) = (temp2 + j + 0.5) *                          &
     &                        delta_phi * Earth_radius

              ! Calculate E-W distance from grid centre for u-cols
              work4(i,j) = (temp1 + i) *                                &
     &                        delta_lambda * Earth_radius

            End Do
          End Do

        Else ! .NOT. trivial_trigs

          ! cos_theta_latitude(i,j) only extends to off_x,off_y
          ! so need to recalculate if not trivial trigs !
          Do j = 1-halo_j, rows+halo_j
            gj = l_datastart(2) + j - 1
            Do i = 1-halo_i, row_length+halo_i

              ! latitude of u,p,theta points
              work1(i,j) = (base_phi + (gj-1) * delta_phi)

              ! Calculate E-W distance from grid centre for p,theta cols
              work2(i,j) = (temp1 + i + 0.5) * cos(work1(i,j)) *        &
     &                        delta_lambda * Earth_radius

              ! Calculate N-S distance from grid centre for p,theta rows
              work3(i,j) = (temp2 + j + 0.5) *                          &
     &                        delta_phi * Earth_radius

              ! Calculate E-W distance from grid centre for u-cols
              work4(i,j) = (temp1 + i) * cos(work1(i,j)) *              &
     &                        delta_lambda * Earth_radius
            End Do
          End Do

        End If ! on trivial_trigs

        ! Calculate N-S distance from grid centre for v-rows
        Do j = 1-halo_j, n_rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            work5(i,j) = (temp2 + j) * delta_phi * Earth_radius
          End Do
        End Do

!----------------------------------------------------------------------
! Define deformation wind field
!----------------------------------------------------------------------
! Distance from centre point function is zero at centre point
! and is scaled by 1.e6 so that u_in refers to wind speed change
! in m/s per 1000km

        ! Define 3D varying u-wind field
        Do k = 1,model_levels
          Do j = 1-halo_j, rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              u_adv(i,j,k) = u_adv(i,j,k) * work4(i,j) / 1.e6
            End Do
          End Do
        End Do

        ! Define 3D varying v-wind field
        Do k = 1,model_levels
          Do j = 1-halo_j, n_rows+halo_j
            Do i = 1-halo_i, row_length+halo_i
              v_adv(i,j,k) = -1.0 * v_adv(i,j,k) * work5(i,j) / 1.e6
            End Do
          End Do
        End Do

! ---------------------------------------------------------------------
! Calculate geostrophic pressure due to deformation wind field
! ---------------------------------------------------------------------
! NOTES:
! qcl(1:model_levels) is theta(1:model_levels)
! w_adv(0:model_levels) is exner_rho_levels(1:model_levels+1)
!
! This section of code assumes an f-plane !
!
!     Interpolate f3_at_u to p-grid
      Do j = 1, rows
        Do i = 1,row_length
          f3_at_p(i,j)= 0.5 * ( f3_at_u(i,j) + f3_at_u(i-1,j) )
        End Do
      End Do

      Do k = 2, model_levels

        ! u_ref is in m/s/1000km
        ! Divide by 1.e6 to get alpha in m/s/m
        alpha = -1.0*u_ref(k) /1.e6

        Do j = 1-halo_j, rows+halo_j
          Do i= 1-halo_i, row_length+halo_i

            ! Calculate geopotential change on pressure surface

            ! Use f3_at_p(1,1) at present because f3_at_p not defined
            ! in halos. Could change this in the future.
            geop = (f3_at_p(1,1)*alpha*alpha*work2(i,j)*work3(i,j)      &
     &             - alpha**2 * (work2(i,j)**2 + work3(i,j)**2)/2.)

            ! Calculate theta on rho/u/p level
            work1(i,j) = ( (r_rho_levels(i,j,k) -                       &
     &                                   r_theta_levels(i,j,k-1)) *     &
     &                                   qcl(i,j,k) +                   &
     &                                   (r_theta_levels(i,j,k) -       &
     &                                   r_rho_levels(i,j,k)) *         &
     &                                   qcl(i,j,k-1) ) /               &
     &                                   (r_theta_levels(i,j,k) -       &
     &                                    r_theta_levels(i,j,k-1))
            ! Calculate average theta
            theta2 = work1(i,j) +                                       &
     &             (qcl(i,j,k)-qcl(i,j,k-1))*geop/                      &
     &             (g*(r_theta_levels(i,j,k)-r_theta_levels(i,j,k-1)))

            theta3 = work1(i,j) *                                       &
     &               (1. + (1./repsilon -1.) * q(i,j,k))


            ! Calculate pressure change on height surface
            ! Store temporarily in qcf
            qcf(i,j,k-1) = geop/(Cp * theta3)

            ! If top level set top level of qcf
            If (k  ==  model_levels) qcf(i,j,model_levels)=0.0

          End Do
        End Do
      End Do

! =====================================================================
!
!                     Set up Upper Level Zonal Jet
!
! =====================================================================
!

        ! t_horizfn_data(1) = jet core maximum  (e.g. 60m/s)
        ! t_horizfn_data(2) = width of jet (e.g. 800000m)
        ! t_horizfn_data(3) = height of maximum jet (e.g. 10000m)
        ! t_horizfn_data(4) = height jet goes to zero (e.g. 15000m)

        ! Force jet to go to zero at the top rho level (model_levels+1)
        k = model_levels
        temp = (eta_theta_levels(k) - eta_rho_levels(k) +               &
     &          eta_theta_levels(k)) * height_domain
        t_horizfn_data(4) = temp

        ! Set up vertical varying scaling factor
        Do k = 1, model_levels
          temp = eta_rho_levels(k) * height_domain
          If (temp  <   t_horizfn_data(3)) Then
            u_ref(k) = temp*t_horizfn_data(1)/t_horizfn_data(3)
            jetmaxlevel = k
          Else If (temp  <   t_horizfn_data(4)) Then
            u_ref(k) = (t_horizfn_data(4)-temp)*                        &
     &                  t_horizfn_data(1)/                              &
     &                 (t_horizfn_data(4)-t_horizfn_data(3))
          Else
            u_ref(k) = 0.0
          End If
        End Do ! k=1, model_levels

        ! Set up ramp function for theta across domain
        ! centred on domain centre with magnitude t_horizfn_data(1)
        ! and width t_horizfn_data(2)

        Do j = 1-halo_j, global_rows+halo_j
          ! Calculate distance from the domain centre point
          temp = (-1.0*(global_rows+1)/2.0 + j-1)                       &
     &            *delta_phi*earth_radius
          If (temp  <   -0.5*t_horizfn_data(2)) Then
            ujet_j(j) = 0.0
          Else If (temp  <   0.5*t_horizfn_data(2)) Then
            ujet_j(j) = COS(Pi*temp/t_horizfn_data(2))
          Else
            ujet_j(j) = 0.0
          End If

        End Do

! Add jet to u-wind field
        Do k = 1, model_levels
          Do j = 1-halo_j, rows+halo_j
            Do i= 1-halo_i, row_length+halo_i
              u_adv(i,j,k) = u_adv(i,j,k) +                             &
     &                       u_ref(k)*ujet_j(l_datastart(2)+j-1)
            End Do
          End Do
        End Do

! ----------------------------------------------------------------------
! Set up a temperature in thermal wind balance with the jet
! ----------------------------------------------------------------------
!NOTES:
! Starts from the southern (warm) boundary,
! so the temperature decreases towards the North
! (Only correct if there is no orography)

        temp = delta_phi*earth_radius
        Do k = 1, model_levels-1
          i=1
          thetabal(1-halo_j) = qcl(1,1,k)
          Do j = 1-halo_j+1, global_rows+halo_j
            thetabal(j) = thetabal(j-1) -                               &
     &         ((u_ref(k+1)*ujet_j(j))-(u_ref(k)*ujet_j(j)))*           &
     &         (f3_at_u(1,1)*qcl(i,j,k)*temp)                           &
     &         /(g*(r_rho_levels(1,1,k+1) - r_rho_levels(1,1,k)))
          End Do

          Do j = 1-halo_j, rows+halo_j
            Do i= 1-halo_i, row_length+halo_i
              qcl(i,j,k) = thetabal(l_datastart(2)+j-1)
            End Do
          End Do

        End Do

        ! Top level, assumes jet goes to zero at model_levels+1
        k = model_levels
        i=1
        ! temp1 is height of top rho level (model_levels+1)
        temp1 = (eta_theta_levels(k) - eta_rho_levels(k) +              &
     &           eta_theta_levels(k)) * height_domain
        thetabal(1-halo_j) = qcl(1,1,k)
        Do j = 1-halo_j+1, global_rows+halo_j
          thetabal(j) = thetabal(j-1) +                                 &
     &         u_ref(k)*ujet_j(j)*                                      &
     &         f3_at_u(1,1)*qcl(i,j,k)*temp                             &
     &         /(g*(temp1 - eta_rho_levels(k)*height_domain))
        End Do

        Do j = 1-halo_j, rows+halo_j
          Do i= 1-halo_i, row_length+halo_i
            qcl(i,j,k) = thetabal(l_datastart(2)+j-1)
          End Do
        End Do


! ---------------------------------------------------------------------
! Set up balanced pressure profile (from hydrostatic relation)
! ---------------------------------------------------------------------

! w_adv(0:model_levels) is exner_rho_levels(1:model_levels+1)
! qcl(1:model_levels) is theta(1:model_levels)

      ! Set first level pressure
      do j = 1-halo_j, rows+halo_j
        do i = 1-halo_i, row_length+halo_i

          w_adv(i,j,0) =  1.0 -                                         &
     &           g*(r_rho_levels(i,j,1) - r_theta_levels(i,j,0))        &
     &            /(Cp * qcl(i,j,1))

        enddo
      enddo

      ! Set all intermediate levels
      do k = 2, model_levels
        do j = 1-halo_j, rows+halo_j
          do i = 1-halo_i, row_length+halo_i

            w_adv(i,j,k-1) =  w_adv(i,j,k-2) -                          &
     &           g*(r_rho_levels(i,j,k) - r_rho_levels(i,j,k-1))        &
     &            /(Cp * qcl(i,j,k-1))

          enddo
        enddo
        Write(6,*) 'p:',k,w_adv(1,1,k-1),qcl(1,1,k-1)
      enddo

      ! Set top level
      k = model_levels + 1
        do j = 1-halo_j, rows+halo_j
          do i = 1-halo_i, row_length+halo_i
          ! top exner level is same distance above top theta level as
          ! the rho level is below.

          w_adv(i,j,k-1)= w_adv(i,j,k-2)-                               &
     &         g*2.0*(r_theta_levels(i,j,k-1)-r_rho_levels(i,j,k-1))    &
     &         /(Cp*qcl(i,j,k-1))

        enddo
      enddo

! ---------------------------------------------------------------------
! Add deformation pressure balance increments
! ---------------------------------------------------------------------

      Do k = 1, model_levels
        Do j = 1-halo_j, rows+halo_j
          Do i = 1-halo_i, row_length+halo_i
            w_adv(i,j,k) = w_adv(i,j,k) + qcf(i,j,k)
            ! Set lowest level
            If (k  ==  1) w_adv(i,j,0)=w_adv(i,j,0)+qcf(i,j,1)
          End Do
        End Do
      End Do

! ---------------------------------------------------------------------
! Set up exner_theta_levels from exner_rho_levels
! ---------------------------------------------------------------------

       Do k = 1, model_levels - 1
         Do j = 1-halo_j, rows+halo_j
           Do i = 1-halo_i, row_length+halo_i
           weight = (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))/    &
     &                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
           exner_theta_levels(i,j,k) =   weight *                       &
     &                                   w_adv(i,j,k-1) +               &
     &                  (1.0 - weight) * w_adv(i,j,k)
           End Do
         End Do
        End Do

         k =  model_levels
         Do j = 1-halo_j, rows+halo_j
           Do i = 1-halo_i, row_length+halo_i
           exner_theta_levels(i,j,k) =  0.5 * w_adv(i,j,k-1) +          &
     &                                  0.5 * w_adv(i,j,k)
          End Do
         End Do

      IF (lhook) CALL dr_hook('IDL_SETUP_DEF_FRONT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_setup_def_front
