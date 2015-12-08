! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine IDL_BAROCLINIC

      Subroutine IDL_BAROCLINIC(                                        &
     &                      row_length, rows, n_rows, model_levels      &
     &,                     global_row_length, global_rows              &
     &,                     l_datastart                                 &
     &,                     halo_i, halo_j, off_x, off_y                &
     &,                     delta_lambda, delta_phi                     &
     &,                     base_phi                                    &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     q, theta, u_adv, v_adv, exner               &
     &,                     exner_theta_levels                          &
     &,                     f3_at_u                                     &
     &,                     t_horizfn_data                              &
     &,                     theta_ref                                   &
     &,                     L_trivial_trigs                             &
     &                      )

      USE earth_constants_mod, ONLY: g, earth_radius

      USE atmos_constants_mod, ONLY: cp
      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      Implicit None

!  Description:
!           A routine to initialise a limited area domain
!           with a baroclinic jet and tropopause perturbation
!           in order to simulate baroclinic instability.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
! Declarations:

! input scalars
      Integer, intent(in) ::                                            &
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

      Real, intent(in) ::                                               &
     &  delta_lambda                                                    &
     &, delta_phi                                                       &
     &, base_phi


      ! vertical co-ordinate information
      Real, intent(in) ::                                               &
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, model_levels)                &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)                                    &
     &, f3_at_u (1-off_x:row_length+off_x,                              &
     &                      1-off_y:rows+off_y)

      Real, intent(in) ::                                               &
     & theta_ref(model_levels)

      Real, intent(in) ::                                               &
     &  t_horizfn_data(10)  ! Data values describing horizontal t fn

      Logical, intent(in) ::                                            &
     &  L_trivial_trigs

!inout arrays
      Real, intent(inout) ::                                            &
     &  q(1-halo_i:row_length+halo_i                                    &
     &,       1-halo_j:rows+halo_j, model_levels)                       &
     &, theta(1-halo_i:row_length+halo_i                                &
     &,       1-halo_j:rows+halo_j, model_levels)                       &
     &, u_adv(1-halo_i:row_length+halo_i                                &
     &,       1-halo_j:rows+halo_j, model_levels)                       &
     &, v_adv(1-halo_i:row_length+halo_i                                &
     &,       1-halo_j:n_rows+halo_j, model_levels)                     &
     &, exner(1-halo_i:row_length+halo_i                                &
     &,       1-halo_j:rows+halo_j, 0:model_levels)                     &
     &, exner_theta_levels(1-halo_i:row_length+halo_i                   &
     &,                    1-halo_j:rows+halo_j, model_levels)


! Local constants
      Real ::                                                           &
     & th_00, l_y, h_t, x_0, corr

! indicies of bottom left corner on this PE
      Integer ::                                                        &
     &  gi, gj

      Real ::                                                           &
     &  weight                                                          &
     &, dx,dy, temp1, temp2

      INTEGER :: index,indexi


      REAL :: mm,ll,eta,aa,wi

! Local variables
      Real ::                                                           &
     &  z                                                               &
            ! vertical height
     &, r   ! radial distance

      Integer ::                                                        &
     &  i, j, k

! Fields for the entire domain
      REAL ::                                                           &
     & theta_big(1-halo_i:global_row_length+halo_i                      &
     &,       1-halo_j:global_rows+halo_j, model_levels)                &
     &, u_big(1-halo_i:global_row_length+halo_i                         &
     &,       1-halo_j:global_rows+halo_j, model_levels)                &
     &, v_big(1-halo_i:global_row_length+halo_i                         &
     &,       1-halo_j:global_rows-1+halo_j, model_levels)              &
     &, pi_big(1-halo_i:global_row_length+halo_i                        &
     &,       1-halo_j:global_rows+halo_j, 0:model_levels)


      Real ::                                                           &
     &  work1(1-halo_i:global_row_length+halo_i,                        &
     &        1-halo_j:global_rows+halo_j)                              &
     &, x_th(1-halo_i:global_row_length+halo_i,                         &
     &       1-halo_j:global_rows+halo_j)                               &
     &, y_u(1-halo_i:global_row_length+halo_i,                          &
     &      1-halo_j:global_rows+halo_j)                                &
     &, x_u(1-halo_i:global_row_length+halo_i,                          &
     &      1-halo_j:global_rows+halo_j)                                &
     &, x_v(1-halo_i:global_row_length+halo_i,                          &
     &      1-halo_j:global_rows-1+halo_j)                              &
     &, y_v(1-halo_i:global_row_length+halo_i,                          &
     &      1-halo_j:global_rows-1+halo_j)

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!-----------------------------------------------------------------------

! Centre of domain is (indexi,index)
      IF (lhook) CALL dr_hook('IDL_BAROCLINIC',zhook_in,zhook_handle)
      index=INT(REAL(global_rows)/2.)
      indexi=INT(REAL(global_row_length)/2.)


! (temp1,temp2) is the coordinate of the bottom left hand corner
! for this processor relative  to the centre of the domain
      temp1 = 1 - INT(global_row_length/2.0)
      temp2 = 1 - INT(global_rows/2.0)

      IF (L_trivial_trigs) THEN

        DO j = 1-halo_j, global_rows+halo_j
          DO i = 1-halo_i, global_row_length+halo_i

            x_th(i,j) = (temp1 + i - 1 ) *                              &
     &                        delta_lambda * Earth_radius


         ! Calculate N-S distance from grid centre for p,theta rows
            y_u(i,j) = (temp2 + j - 1 ) *                               &
     &                        delta_phi * Earth_radius


         ! Calculate E-W distance from grid centre for u-cols
            x_u(i,j) = (temp1 + i - 0.5 ) *                             &
     &                        delta_lambda * Earth_radius

          END DO
        END DO

      ELSE ! .NOT. trivial_trigs

        ! cos_theta_latitude(i,j) only extends to off_x,off_y
        ! so need to recalculate if not trivial trigs !
        DO j = 1-halo_j, rows+halo_j
          gj = l_datastart(2) + j - 1
          DO i = 1-halo_i, row_length+halo_i

            ! latitude of u,p,theta points
            work1(i,j) = (base_phi + gj * delta_phi)

            ! Calculate E-W distance from grid centre for p,theta cols
            x_th(i,j) = (temp1 + i + 0.5) * cos(work1(i,j)) *           &
     &                        delta_lambda * Earth_radius

            ! Calculate N-S distance from grid centre for p,theta rows
            y_u(i,j) = (temp2 + j + 0.5) *                              &
     &                        delta_phi * Earth_radius

            ! Calculate E-W distance from grid centre for u-cols
            x_u(i,j) = (temp1 + i) * cos(work1(i,j)) *                  &
     &                        delta_lambda * Earth_radius
          END DO
        END DO

      END IF ! on trivial_trigs

      ! Calculate N-S distance from grid centre for v-rows
      DO j = 1-halo_j, global_rows-1+halo_j
        DO i = 1-halo_i, global_row_length+halo_i
          y_v(i,j) = (temp2 + j -0.5) * delta_phi * Earth_radius
          x_v(i,j) = (temp1 + i - 1) *                                  &
     &                 delta_lambda * Earth_radius
        END DO
      END DO

!----------------------------------------------------------------------
! Define Baroclinic Jet
! Based on Beare et al 2003 QJRMS
!----------------------------------------------------------------------

! Width of jet
      l_y=t_horizfn_data(2)
! Tropopause height
      h_t=9.0e3
      mm=0.7
      ll=2*pi/l_y
      eta=t_horizfn_data(6)*1.0e6
      aa=1.863
! width of perturbation
      wi=t_horizfn_data(3)
! EW  location of tropopause perturbation
      x_0=t_horizfn_data(4)

      DO k = 1,model_levels
        z=(r_rho_levels(1,1,k) - r_rho_levels(1,1,1))
        DO j = 1-halo_j, global_rows+halo_j
          DO i = 1-halo_i, global_row_length+halo_i
            r=ABS(y_u(i,j))

            IF ((z <= h_t) .AND. (r <= l_y/2.0) ) THEN

              u_big(i,j,k) =  (t_horizfn_data(1)/2.)                    &
     &           * (sinh(mm*z/h_t)/sinh(mm))                            &
     &           * (1 - cos(ll*(y_u(i,j)+l_y/2.0)))

            ELSE IF ((z > h_t).AND.(r <= l_y/2.0)) THEN

              u_big(i,j,k) =  (t_horizfn_data(1)/2.)                    &
     &           * exp(mm*2.0*(1.0-z/h_t))                              &
     &           * (1 - cos(ll*(y_u(i,j)+l_y/2.0)))

            ELSE
              u_big(i,j,k) = 0.0

            END IF
! Vortex perturbation
            r=sqrt((x_u(i,j)-x_0)**2+y_u(i,j)**2)
            IF ((z <= h_t).AND.(r <= wi/2.)) THEN
              u_big(i,j,k) = u_big(i,j,k) -t_horizfn_data(5)*           &
     &           exp(200.0*pi*(z-h_t)/wi) * (y_u(i,j)/MAX(r,0.1))       &
     &           * sin(2*pi*r/wi)
            ELSE IF ((z > h_t).AND.(r <= wi/2.)) THEN
              u_big(i,j,k) = u_big(i,j,k) -t_horizfn_data(5)*           &
     &           exp(400.0*pi*(h_t-z)/wi) * (y_u(i,j)/MAX(r,0.1))       &
     &           * sin(2*pi*r/wi)
            END IF

          END DO
        END DO
      END DO



      DO k = 1,model_levels
        z=(r_rho_levels(1,1,k) - r_rho_levels(1,1,1))
        DO j = 1-halo_j, global_rows-1+halo_j
          DO i = 1-halo_i, global_row_length+halo_i

            v_big(i,j,k) = 0.0

! Vortex perturbation
            r=sqrt((x_v(i,j)-x_0)**2+y_v(i,j)**2)
            IF ((z <= h_t).AND.(r <= wi/2.)) THEN
              v_big(i,j,k) = v_big(i,j,k) +t_horizfn_data(5)*           &
     &          exp(200.0*pi*(z-h_t)/wi) * ((x_v(i,j)-x_0)/MAX(r,0.1))  &
     &          * sin(2*pi*r/wi)
            ELSE IF  ((z > h_t).AND.(r <= wi/2.)) THEN
              v_big(i,j,k) = v_big(i,j,k) +t_horizfn_data(5)*           &
     &          exp(400.0*pi*(h_t-z)/wi) * ((x_v(i,j)-x_0)/MAX(r,0.1))  &
     &          * sin(2*pi*r/wi)
            END IF

          END DO
        END DO
      END DO

      dx = delta_lambda*earth_radius
      dy = delta_phi*earth_radius
      th_00=theta_ref(1)

! x integration of thermal wind balance

      DO k = 1, model_levels-1
        DO j = 1-halo_j, global_rows-1+halo_j
          theta_big(global_row_length+halo_i,j,k)=theta_ref(k)
          DO i= global_row_length+halo_i-1,  1-halo_i,-1
              theta_big(i,j,k) =  theta_big(i+1,j,k)  -                 &
     &         (v_big(i,j,k+1)-v_big(i,j,k))*                           &
     &         (f3_at_u(1,1)*th_00*dx)                                  &
     &         /(g*(r_rho_levels(1,1,k+1) - r_rho_levels(1,1,k)))
          END DO
        END DO
      END DO

! y integration of thermal wind balance
      DO k = 1, model_levels-1
        DO j = 1-halo_j+1, global_rows+halo_j
          DO i = 1-halo_i, global_row_length+halo_i
             theta_big(i,j,k) =  theta_big(i,j-1,k)  -                  &
     &           (u_big(i,j,k+1)-u_big(i,j,k))*                         &
     &           (f3_at_u(1,1)*th_00*dy)                                &
     &           /(g*(r_rho_levels(1,1,k+1) - r_rho_levels(1,1,k)))
          END DO
        END DO
      END DO

! top level theta
      DO j = 1-halo_j, global_rows+halo_j
        DO i= 1-halo_i, global_row_length+halo_i
          theta_big(i,j,model_levels) =  theta_ref(model_levels)
        END DO
      END DO
! Correction to make surface potential
! temperature in centre of domain=theta_ref(1)
      corr = theta_ref(1)-theta_big(indexi,index,1)
      DO k = 1, model_levels
        DO j = 1-halo_j, global_rows+halo_j
          DO i = 1-halo_i, global_row_length+halo_i
           theta_big(i,j,k)=theta_big(i,j,k) +corr
          END DO
        END DO
      END DO

! ---------------------------------------------------------------------
! Set up balanced pressure profile (from hydrostatic relation)
! ---------------------------------------------------------------------

! Set final  level pressure
      DO j = 1-halo_j, global_rows+halo_j
        DO i = 1-halo_i, global_row_length+halo_i

          pi_big(i,j,model_levels) = 0.0

        END DO
      END DO

      ! Set all intermediate levels
      DO k = model_levels-1,1,-1
        DO j = 1-halo_j, global_rows+halo_j
          DO i = 1-halo_i, global_row_length+halo_i

            pi_big(i,j,k) =  pi_big(i,j,k+1) +                          &
     &           g*(r_rho_levels(1,1,k+1) - r_rho_levels(1,1,k))        &
     &            /(Cp * theta_big(i,j,k))

          END DO
        END DO
      END DO

      DO j = 1-halo_j, global_rows+halo_j
        DO i = 1-halo_i, global_row_length+halo_i

            pi_big(i,j,0) =  pi_big(i,j,1) +                            &
     &           g*(r_rho_levels(1,1,1) - r_theta_levels(1,1,0))        &
     &            /(Cp * theta_big(i,j,1))

        END DO
      END DO

!make pi_big 1 at surface in middle  of domain
      corr=1.0-pi_big(indexi,index,0)
      DO k = 0,model_levels
        DO j = 1-halo_j, global_rows+halo_j
          DO i = 1-halo_i, global_row_length+halo_i

            pi_big(i,j,k) =  pi_big(i,j,k) +  corr

          END DO
        END DO
      END DO

! Distribute section of the _big fields to this processor
! with starting index (l_datastart(1),l_datastart(2))

      DO k = 1, model_levels
        DO j = 1-halo_j, rows+halo_j
          gj=l_datastart(2)+j-1
          DO i = 1-halo_i, row_length+halo_i
             gi=l_datastart(1)+i-1
             u_adv(i,j,k)=u_big(gi,gj,k)
             theta(i,j,k)=theta_big(gi,gj,k)
          END DO
        END DO
      END DO

      DO k = 1, model_levels
        DO j = 1-halo_j, n_rows+halo_j
          gj=l_datastart(2)+j-1
          DO i = 1-halo_i, row_length+halo_i
             gi=l_datastart(1)+i-1
             v_adv(i,j,k)=v_big(gi,gj,k)
          END DO
        END DO
      END DO

      DO k = 0, model_levels
        DO j = 1-halo_j, rows+halo_j
          gj=l_datastart(2)+j-1
          DO i = 1-halo_i, row_length+halo_i
             gi=l_datastart(1)+i-1
             exner(i,j,k)=pi_big(gi,gj,k)
          END DO
        END DO
      END DO


! ---------------------------------------------------------------------
! Set up exner_theta_levels from exner_rho_levels
! ---------------------------------------------------------------------

      DO k = 1, model_levels - 1
        DO j = 1-halo_j, rows+halo_j
          DO i = 1-halo_i, row_length+halo_i
          weight = (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))/     &
     &                (r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k))
          exner_theta_levels(i,j,k) =   weight *                        &
     &                                   exner(i,j,k-1) +               &
     &                  (1.0 - weight) * exner(i,j,k)
          END DO
        END DO
      END DO

      k =  model_levels
      DO j = 1-halo_j, rows+halo_j
        DO i = 1-halo_i, row_length+halo_i
          exner_theta_levels(i,j,k) =  0.5 * exner(i,j,k-1) +           &
     &                                  0.5 * exner(i,j,k)
        END DO
      END DO

      IF (lhook) CALL dr_hook('IDL_BAROCLINIC',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_BAROCLINIC
