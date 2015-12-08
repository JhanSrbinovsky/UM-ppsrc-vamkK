! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine IDL_pr_balance

      Subroutine IDL_pr_balance(                                        &
     &                      model_domain, row_length, rows, n_rows      &
     &,                     model_levels                                &
     &,                     delta_x, delta_y, halo_i, halo_j            &
     &,                     me, l_datastart                             &
     &,                     delta_lambda, delta_phi                     &
     &,                     lambda_p, phi_p                             &
     &,                     L_regular                                   &
     &,                     Base_phi, base_lambda                       &
     &,                     eta_theta_levels, eta_rho_levels            &
     &,                     r_theta_levels, r_rho_levels                &
     &,                     theta, exner_rho_levels, u_adv, v_adv       &
     &,                     theta_ref, exner_ref                        &
! Profile settings
     &,                     height_domain, theta_surface, Brunt_Vaisala &
     &,                     u_in, v_in, u_ramp_start, u_ramp_end        &
     &,                     ujet_lat, ujet_width                        &
     &,                     f_plane, r_plane                            &
!  Options
     &,                     L_trivial_trigs, L_code_test)

! Purpose:
!         For a given wind field, derives balanced (Exner) pressures and
!         potential temperatures by inverting geostrophic relation.
!        To work on Cartesian grid  (L_trivial_trigs .true.)
!          set f_plane in idealise NAMELIST (default is -89.0 - S Pole)
!          Both u and v may be set
!         NB Model will behave as on Cartesian grid but ancillaries
!           with real orography will not be adjusted and METVIEW will
!           treat output as if on spherical grid.
!         To work in spherical coordinates v = 0 and  either
!               (in if test order)
!             1.   u and f = constant
!                     set f_plane  >   -89.0 in IDEALISE NAMELIST
!
!         or  2.   u = const.* (cos(lat) - cos(ujet_lat + ujet_width))/
!                         (cos(ujet_lat) - cos(ujet_lat + ujet_width))
!      for  ujet_lat < latitude <  ujet_lat + ujet_width
!                  u = const.* (cos(lat) - cos(ujet_lat - ujet_width))/
!                         (cos(ujet_lat) - cos(ujet_lat - ujet_width))
!      for  ujet_lat - ujet_width < latitude <  ujet_lat
!      and    u = 0 for  latitude > ujet_lat + ujet_width
!             u = 0 for  latitude < ujet_lat - ujet_width
!             u = 0 for  latitude < ujet_lat - ujet_width
!                    set  ujet_lat > -90.0 and ujet_width =90.0
!                     in IDEALISE NAMELIST
!                     (hence only use for GLOBAL configurations)
!
!         or  3.   u = const.*
!                      (cos^2(lat) - cos^2(ujet_lat + ujet_width))/
!                      (cos^2(ujet_lat) - cos^2(ujet_lat + ujet_width))
!      for  ujet_lat < latitude <  ujet_lat + ujet_width
!                  u = const.*
!                      (cos^2(lat) - cos^2(ujet_lat - ujet_width))/
!                      (cos^2(ujet_lat) - cos^2(ujet_lat - ujet_width))
!      for  ujet_lat - ujet_width < latitude <  ujet_lat
!      and    u = 0 for  latitude > ujet_lat + ujet_width
!             u = 0 for  latitude < ujet_lat - ujet_width
!             u = 0 for  latitude < ujet_lat - ujet_width
!                    set  ujet_lat > -90.0 and ujet_width < 89.0
!                         in IDEALISE NAMELIST
!
!        or   4.   u = constant and variable  f
!                  do NOT set f_plane or ujet_lat in IDEALISE NAMELIST
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      USE earth_constants_mod, ONLY: g, earth_radius, two_omega

      USE atmos_constants_mod, ONLY:                                    &
          r, cp, p_zero 

      USE conversions_mod, ONLY: pi_over_180,pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE


      Integer                                                           &
     &  model_domain                                                    &
                         ! holds integer code for model domain
     &, row_length                                                      &
                         ! number of points on a row
     &, rows                                                            &
                         ! number of rows in a theta field
     &, n_rows                                                          &
                             ! Local number of rows in a v field
     &, model_levels                                                    &
                         ! number of model levels
     &, halo_i                                                          &
                             ! Size of halo in i direction.
     &, halo_j               ! Size of halo in j direction.

      Real                                                              &
     &  theta_surface                                                   &
     &, height_domain                                                   &
     &, u_in(4)                                                         &
                   ! Input values of zonal u
     &, v_in(4)                                                         &
                   ! Input values of southerly wind v
     &, ujet_lat                                                        &
                    ! To specify centre latitude (degrees) of jet core
     &, ujet_width                                                      &
                      ! To specify width (degrees) of jet
     &, Brunt_Vaisala                                                   &
     &, delta_x, delta_y                                                &
                              ! Resolution at equator
     &, u_ramp_start                                                    &
                        ! ramping starting latitude for u
     &, u_ramp_end      ! ramping ending latitude for u

      Logical                                                           &
     &  L_code_test                                                     &
                       ! user switch
     &, L_trivial_trigs    !  makes grid Cartesian if .true.

      Integer                                                           &
     &  me                                                              &
                   ! My processor number
     &, l_datastart(2)       ! First gridpoints held by this processor


      Real                                                              &
           ! horizontal co-ordinate information
        delta_lambda                                                    &
      , delta_phi                                                       &
      , base_phi                                                        &
      , base_lambda                                                     &
      , f_plane                                                         &
                      ! reference latitude f_plane
      , r_plane       ! reference latitude for row 1 (bottom row)

      Real                                                              &
           !VarRes horizontal co-ordinate information
     &  lambda_p(1-halo_i:row_length+halo_i)                            &
     &, phi_p    ( 1-halo_i : row_length + halo_i                       &
     &,            1-halo_j : rows + halo_j )

      Logical, Intent(In) :: L_regular

      Real                                                              &
           ! vertical co-ordinate information
     &  r_theta_levels(1-halo_i:row_length+halo_i,                      &
     &                 1-halo_j:rows+halo_j,0:model_levels)             &
     &, r_rho_levels(1-halo_i:row_length+halo_i,                        &
     &               1-halo_j:rows+halo_j, model_levels)                &
     &, eta_theta_levels(0:model_levels)                                &
     &, eta_rho_levels(model_levels)

! Primary Arrays with large haloes used as work arrays
! Needed so that external halo values in LAMs can be set correctly
! Use w_adv for exner_rho_levels, qcl for theta
! u_adv may be modified if jet-like option chosen in NAMELIST
      Real                                                              &
     &  exner_rho_levels(1-halo_i:row_length+halo_i                     &
     &,                  1-halo_j:rows+halo_j, model_levels +1)         &
     &, theta(1-halo_i:row_length+halo_i                                &
     &,       1-halo_j:rows+halo_j, model_levels)                       &
     &, u_adv(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,         &
     &        model_levels)                                             &
     &, v_adv(1-halo_i:row_length+halo_i, 1-halo_j:n_rows+halo_j,       &
     &        model_levels)                                             &
     &, theta_ref(model_levels)                                         &
                                 !theta profile for use in sponge & lbcs
     &, exner_ref(model_levels + 1)  ! Exner profile for use in lbcs

! local variables
      Integer                                                           &
     &  i, j, k                                                         &
     &, gj, gi     ! loop counters relative to global origin

                   ! Cartesian coordinates
      Real                                                              &
     &  x(1-halo_i : row_length+halo_i)                                 &
     &, y(1-halo_j : rows+halo_j)                                       &
     &, lat(1-halo_j : rows+halo_j)                                     &
     &, coslat(1-halo_j : rows+halo_j)

      Real                                                              &
        recip_delta_lambda                                              &
      , recip_delta_phi                                                 &
      , BV_squared                                                      &
      , BV_squared_over_g                                               &
      , tolerance                                                       &
      , foverg                                                          &
      , usquared                                                        &
      , biggterm                                                        &
      , corterm                                                         &
      , temp                                                            &
      , weight                                                          &
      , cos_ramp_end                                                    &
      , ujet_rad                                                        &
                    ! centre latitude of jet core in radians
      , ujet_rwidth                                                     &
                       !  width  of jet in radians
      , constant                                                        &
      , conN                                                            &
      , conS                                                            &
      , clat                                                            &
      , coslat2                                                         &
      , coslat4                                                         &
      , cosM                                                            &
      , cosN                                                            &
      , cosS                                                            &
      , cosX                                                            &
      , denom                                                           &
      , M2, M4, N2, N4, S2, S4, X2, X4, DN, DS

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


! ----------------------------------------------------------------------
! Section 0.  Initialise Data fields.
! ----------------------------------------------------------------------

      IF (lhook) CALL dr_hook('IDL_PR_BALANCE',zhook_in,zhook_handle)
      tolerance = 0.0001

      if(me  ==  0)then
        print*,' Balanced pressure and temperature fields '
        print*,' Zonal wind must be a function of latitude only '
        print*,' Southerly wind (v-component) must be zero '
        print*,' For global configurations restricted options '
        print*,' Do not set f_plane and set values for ujet_lat '
        print*,' and  ujet_width '
      Endif   !(me  ==  0)

      BV_squared =  Brunt_Vaisala * Brunt_Vaisala
      BV_squared_over_g =  BV_squared / g
      usquared = u_in(1)* u_in(1)
      biggterm = g *  g  / (BV_squared * Cp * theta_surface)
      recip_delta_lambda = 1./delta_lambda
      recip_delta_phi = 1./delta_phi
      ujet_rad = ujet_lat * Pi / 180.0
      ujet_rwidth = ujet_width * Pi / 180.0
!   cos_ramp_end defines latitude polewards from which  u = 0
      cos_ramp_end   = 0.0
      if( u_ramp_start  >   0.1)then
        cos_ramp_end   = cos(Pi_over_180 * u_ramp_end)
      endif ! u_ramp_start  >   0.1

!-----------------------------------------------------------------
!  Section 2 Derive pressure field
! Either balance pressure field  (Geostrophic - see comments below)
!-------------------------------------------------------------------

      if(L_trivial_trigs ) then

!       foverg = f3_at_u(1,1) / g
        foverg = two_omega * sin(f_plane * Pi_over_180) / g

!    Origin (reference point) is bottom left corner
! Haloed arrays are needed to set external haloes for lbcs
        If (L_regular) Then
          do j = 1-halo_j, rows+halo_j
            gj = l_datastart(2) + j - 1
            y(j) = gj * delta_y
          enddo
          do i = 1-halo_i, row_length+halo_i
            gi = l_datastart(1) + i - 1
            x(i) = gi * delta_x
          enddo
        else !  variable resolution
          do j = 1-halo_j, rows+halo_j
            y(j) = (phi_p(1,j) - base_phi) * Earth_radius
          enddo
          do i = 1-halo_i, row_length+halo_i
            x(i) = (lambda_p(i) - base_lambda) * Earth_radius
          enddo
        end If !  L_regular
        do k = 1, model_levels
          do j = 1-halo_j, rows+halo_j
            do i = 1-halo_i, row_length+halo_i
              theta(i,j,k) = theta_surface * exp ( BV_squared_over_g *  &
     &                      (r_theta_levels(i,j,k) - Earth_radius +     &
     &                    foverg * (u_in(1) * y(j) - v_in(1) * x(i))) )
              exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /     &
     &                       exp ( BV_squared_over_g *                  &
     &                      (r_rho_levels(i,j,k) - Earth_radius +       &
     &                    foverg * (u_in(1) * y(j) - v_in(1) * x(i))) )
            enddo
          enddo
        enddo  !  k = 1, model_levels

        k = model_levels + 1
          do j = 1-halo_j, rows+halo_j
            do i = 1-halo_i, row_length+halo_i
              exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /     &
     &                                    exp ( BV_squared_over_g *     &
     &                            ( 2.0 * r_theta_levels(i,j,k-1) -     &
     &                       r_rho_levels(i,j,k-1) - Earth_radius +     &
     &                   foverg * (u_in(1) * y(j) - v_in(1) * x(i))) )
            enddo
          enddo

        do k = 1, model_levels
          do j = 1-halo_j, rows+halo_j
            do i = 1-halo_i, row_length+halo_i
              u_adv(i,j,k) = u_in(1)
            enddo
          enddo
          do j = 1-halo_j, n_rows+halo_j
            do i = 1-halo_i, row_length+halo_i
              v_adv(i,j,k) = v_in(1)
            enddo
          enddo
        enddo    !  k = 1, model_levels

      else    !   spherical coordinate options
!   In spherical geometry, require v=0 everywhere
        do k = 1, model_levels
          do j = 1-halo_j, n_rows+halo_j
            do i = 1-halo_i, row_length+halo_i
              v_adv(i,j,k) = 0.0
            enddo
          enddo
        enddo

      if( f_plane  >   -89.0 ) then   ! Option 1
!  u and  f must be constant
!    therefore will only work for unrotated LAM on f-plane
!    and no need to average these quantities to p-points etc.
        if(me  ==  0)then
          print*,' Zonal wind u and Coriolis f must be constant '
          print*,' Must be a LAM domain; f_plane = ',f_plane
        Endif   !(me  ==  0)

!   v = 0 and u = const. , f = const. everywhere
!  For f = const, set f_plane latitude in NAMELIST
!       foverg = f3_at_u(1,1) / g
        foverg = two_omega * sin(f_plane * Pi_over_180) / g

!    Origin (reference point) is bottom left corner
! Haloed arrays are needed to set external haloes for lbcs
        If (L_regular) Then
          do j = 1-halo_j, rows+halo_j
            gj = l_datastart(2) + j - 1
            coslat(j) = cos(Base_phi + gj * delta_phi)
            y(j) = gj * delta_y
          end do !  j = 1-halo_j, rows+halo_j
        Else  ! variable resolution
          do j = 1-halo_j, rows+halo_j
            coslat(j) = cos( phi_p(1,j) )
            y(j) = (phi_p(1,j) - base_phi) * Earth_radius
          end do !  j = 1-halo_j, rows+halo_j
        Endif  !  L_regular
        do k = 1, model_levels
          do j = 1-halo_j, rows+halo_j
            If (coslat(j)  >   cos_ramp_end) then
              temp =  u_in(1) * ( foverg * y(j)  -                     &
     &                         u_in(1) * log(coslat(j)) / g )
            else
              temp =  0.0
            endIf  ! cos_lat(j)  <=  cos_ramp_end
            do i = 1-halo_i, row_length+halo_i
              u_adv(i,j,k) = u_in(1)
              theta(i,j,k) = theta_surface * exp ( BV_squared_over_g *  &
     &                      (r_theta_levels(i,j,k) - Earth_radius       &
     &                                    + temp ) )
              exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /     &
     &                       exp ( BV_squared_over_g *                  &
     &                      (r_rho_levels(i,j,k) - Earth_radius         &
     &                                    + temp ) )
            enddo
          enddo
        enddo !  k = 1, model_levels

        k = model_levels + 1
        do j = 1-halo_j, rows+halo_j
          If (coslat(j)  >   cos_ramp_end) then
            temp =  u_in(1) * ( foverg * y(j)  -                       &
     &                         u_in(1) * log(coslat(j)) / g )
          else
            temp =  0.0
          endIf  !cos_lat  <=  cos_ramp_end
          do i = 1-halo_i, row_length+halo_i
            exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /       &
     &                      exp ( BV_squared_over_g *                   &
     &                ( 2.0 * r_theta_levels(i,j,k-1) -                 &
     &                 r_rho_levels(i,j,k-1) - Earth_radius             &
     &                                         + temp ) )
          enddo
        enddo

      elseif( ujet_lat  >   -89.0 ) then   ! Option 2 or 3
        if( r_plane   <   -89.0 ) then   ! Option 2
! variable Coriolis,
!  u has cosine jet-like structure centred at ujet_lat
          if(me  ==  0)then
            print*,' u has cosine jet-like profile  '
            print*,' ujet_lat = ', ujet_lat,' ujet_width = ', ujet_width
            print*,' Full vertical Coriolis term used '
          Endif   !(me  ==  0)

        corterm = u_in(1) * two_omega * Earth_radius
        cosM = cos(ujet_rad)
        M2 = cosM * cosM
        cosN = cos(ujet_rad + ujet_rwidth )
        N2 = cosN * cosN
        DN = 1.0 / (cosM - cosN)
        conN = u_in(1)* u_in(1) * DN * DN * N2 * (log(cosN) - 1.5)      &
     &                                      - 0.5 * corterm * DN * N2
        cosS = cos(ujet_rad - ujet_rwidth )
        S2 = cosS * cosS
        M2 = cosM * cosM
        DS = 1.0 / (cosM - cosS)
        conS = conN + corterm * DN * (cosN * cosM - 0.5 * M2)           &
     &        -  u_in(1)* u_in(1) * DN * DN *                           &
     &           (N2 * log(cosM) - 2.0 * cosN * cosM + 0.5 * M2)        &
     &              - corterm * DS * (cosS * cosM - 0.5 * M2)           &
     &        +  u_in(1)* u_in(1) * DS * DS *                           &
     &           (S2 * log(cosM) - 2.0 * cosS * cosM + 0.5 * M2)

!    Origin (reference point) is bottom left corner
!    Haloed arrays are needed to set external haloes for lbcs
        If (L_regular) Then
          do j = 1-halo_j, rows+halo_j
            gj = l_datastart(2) + j - 1
            lat(j) = Base_phi + gj * delta_phi
          enddo
        Else  ! variable resolution
          do j = 1-halo_j, rows+halo_j
            lat(j) =  phi_p(1,j)
          enddo
        Endif  !  L_regular
        do k = 1, model_levels
          do j = 1-halo_j, rows+halo_j
            clat = cos(lat(j))
            coslat2 = clat * clat
            if( lat(j)  <   ujet_rad - ujet_rwidth .or.                 &
     &          lat(j)  >   ujet_rad + ujet_rwidth ) then
              do i = 1-halo_i, row_length+halo_i
                u_adv(i,j,k) = 0.0
                theta(i,j,k) = theta_surface * exp ( BV_squared_over_g *&
     &               (r_theta_levels(i,j,k) - Earth_radius))
                exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /   &
     &                                  exp ( BV_squared_over_g *       &
     &                      (r_rho_levels(i,j,k) - Earth_radius))
              end do

            else   !  lat is within jet bounds
! the following code should not be used at equator due to logcoslat term
              if( lat(j)  >   ujet_rad) then
                cosX = cosN
                constant = conN
              else   ! lat(j)  <   ujet_rad
                cosX = cosS
                constant = conS
              endif  ! lat(j)  <   ujet_rad
              denom = 1.0 / (cosM - cosX)
              temp = corterm * denom * clat * (cosX - 0.5 * clat)      &
     &                - u_in(1)* u_in(1) * denom * denom *              &
     &                      ( 0.5 * coslat2 - 2.0 * clat * cosX +       &
     &                                   cosX * cosX  * log(clat) )     &
     &                   + constant
              do i = 1-halo_i, row_length+halo_i
                u_adv(i,j,k) = u_in(1) * (clat - cosX) * denom
                theta(i,j,k) = theta_surface * exp ( BV_squared_over_g *&
     &                     ( r_theta_levels(i,j,k) - Earth_radius +     &
     &                           temp/g ) )
          exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /         &
     &                                   exp ( BV_squared_over_g *      &
     &                 ( r_rho_levels(i,j,k) - Earth_radius +           &
     &                           temp/g ) )
              end do
            endif    !  lat outside jet bounds
          end do
        end do   ! k = 1, model_levels

        k = model_levels + 1
        do j = 1-halo_j, rows+halo_j
          clat = cos(lat(j))
          coslat2 = clat * clat
          if( lat(j)  <   ujet_rad - ujet_rwidth .or.                   &
     &        lat(j)  >   ujet_rad + ujet_rwidth ) then
            do i = 1-halo_i, row_length+halo_i
              exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /     &
     &                                  exp ( BV_squared_over_g *       &
     &              ( 2.0 * r_theta_levels(i,j,k-1) -                   &
     &             r_rho_levels(i,j,k-1) - Earth_radius))
            end do

          else   !  lat(j) is within jet bounds
! the following code should not be used at equator due to logcoslat term
            if( lat(j)  >   ujet_rad) then
              cosX = cosN
              constant = conN
            else   ! lat(j)  <   ujet_rad
              cosX = cosS
              constant = conS
            endif  ! lat(j)  <   ujet_rad
            denom = 1.0 / (cosM - cosX)
              temp = corterm * denom * clat * (cosX - 0.5 * clat)      &
     &                - u_in(1)* u_in(1) * denom * denom *              &
     &                      ( 0.5 * coslat2 - 2.0 * clat * cosX +       &
     &                                   cosX * cosX  * log(clat) )     &
     &                   + constant
            do i = 1-halo_i, row_length+halo_i
          exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /         &
     &                                   exp ( BV_squared_over_g *      &
     &              (( 2.0 * r_theta_levels(i,j,k-1) -                  &
     &             r_rho_levels(i,j,k-1)) - Earth_radius +              &
     &                           temp/g ) )
            end do
          endif    !  lat outside jet bounds
        end do

       else    ! Option 3  r_plane > -89.0 used as switch  Option 3
! variable Coriolis,
!  u has cos squared jet-like structure centred at ujet_lat
        if(me  ==  0)then
         print*,' u has cosine-squared jet-like profile  '
         print*,' ujet_lat = ', ujet_lat,' ujet_width = ', ujet_width
         print*,' Full vertical Coriolis term used '
        Endif   !(me  ==  0)

      corterm = u_in(1) * two_omega * Earth_radius / g
       cosM = cos(ujet_rad)
       M2 = cosM * cosM
       M4 = M2 * M2
       clat = cos(ujet_rad + ujet_rwidth )
       N2 = clat * clat
       N4 = N2 * N2
       DN = 1.0 / (M2 - N2)
       conN = u_in(1)* u_in(1) * DN * DN * N4 *                         &
     &                         ( log(clat) - 0.75 ) / g -               &
     &           2.0 * corterm * DN * clat * N2 /3.0
       clat = cos(ujet_rad - ujet_rwidth )
       S2 = clat * clat
       S4 = S2 * S2
       DS = 1.0 / (M2 - S2)
       conS = u_in(1)* u_in(1) * DS * DS *                              &
     &                ( S4 * log(cosM) - S2 * M2 + 0.25 * M4 ) / g -    &
     &                corterm * DS * cosM * (S2 - M2/3.0) -             &
     &        u_in(1)* u_in(1) * DN * DN *                              &
     &                ( N4 * log(cosM) - N2 * M2 + 0.25 * M4 ) / g +    &
     &                corterm * DN * cosM * (N2 - M2/3.0) + conN
!    Origin (reference point) is bottom left corner
!    Haloed arrays are needed to set external haloes for lbcs
        If (L_regular) Then
          do j = 1-halo_j, rows+halo_j
            gj = l_datastart(2) + j - 1
            lat(j) = Base_phi + gj * delta_phi
          end do
        Else  ! variable resolution
          do j = 1-halo_j, rows+halo_j
            lat(j) =  phi_p(1,j)
          end do
        Endif !  L_regular)
      do k = 1, model_levels
        do j = 1-halo_j, rows+halo_j
          clat = cos(lat(j))
          coslat2 = clat * clat
          coslat4 = coslat2 * coslat2
          if( lat(j)  <   ujet_rad - ujet_rwidth .or.                   &
     &        lat(j)  >   ujet_rad + ujet_rwidth ) then
            do i = 1-halo_i, row_length+halo_i
              u_adv(i,j,k) = 0.0
              theta(i,j,k) = theta_surface * exp ( BV_squared_over_g *  &
     &                       (r_theta_levels(i,j,k) - Earth_radius))
              exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /     &
     &                                  exp ( BV_squared_over_g *       &
     &                       (r_rho_levels(i,j,k) - Earth_radius))
            end do ! i = 1-halo_i, row_length+halo_i

          else   !  lat(j) is within jet bounds
! the following code should not be used at equator due to logcoslat term
            if( lat(j)  >   ujet_rad) then
              X2 = N2
              X4 = N4
              constant = conN
            else   ! lat(j)  <   ujet_rad
              X2 = S2
              X4 = S4
              constant = conS
            endif  ! lat(j)  <   ujet_rad
            denom = 1.0 / (M2 - X2)
            temp =  corterm * denom * clat * (X2 - coslat2/3.0) -      &
     &                     u_in(1)* u_in(1) * denom * denom *           &
     &                      ( 0.25 * coslat4 - coslat2 * X2 +           &
     &                         X4 * log(clat) ) / g + constant
            do i = 1-halo_i, row_length+halo_i
              u_adv(i,j,k) = u_in(1) * (coslat2 - X2) * denom
              theta(i,j,k) = theta_surface * exp ( BV_squared_over_g *  &
     &                  ( r_theta_levels(i,j,k) - Earth_radius + temp))
              exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /     &
     &                                    exp ( BV_squared_over_g *     &
     &                    ( r_rho_levels(i,j,k) - Earth_radius + temp))
            end do  !  i = 1-halo_i, row_length+halo_i
          endif    !  lat outside jet bounds
        end do  ! j = 1-halo_j, rows+halo_j
      end do   ! k = 1, model_levels

        k = model_levels + 1
        do j = 1-halo_j, rows+halo_j
          clat = cos(lat(j))
          coslat2 = clat * clat
          coslat4 = coslat2 * coslat2
          if( lat(j)  <   ujet_rad - ujet_rwidth .or.                   &
     &        lat(j)  >   ujet_rad + ujet_rwidth ) then
            do i = 1-halo_i, row_length+halo_i
              exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /     &
     &                                    exp ( BV_squared_over_g *     &
     &                            ( 2.0 * r_theta_levels(i,j,k-1) -     &
     &                        r_rho_levels(i,j,k-1) - Earth_radius))
            end do !  i = 1-halo_i, row_length+halo_i

          else   !  lat(j) is within jet bounds
! the following code should not be used at equator due to logcoslat term
            if( lat(j)  >   ujet_rad) then
              X2 = N2
              X4 = N4
              constant = conN
            else   ! lat(j)  <   ujet_rad
              X2 = S2
              X4 = S4
              constant = conS
            endif  ! lat(j)  <   ujet_rad
            denom = 1.0 / (M2 - X2)
            temp =  corterm * denom * clat * (X2 - coslat2/3.0) -      &
     &                            u_in(1)* u_in(1) * denom * denom *    &
     &                             ( 0.25 * coslat4 - coslat2 * X2 +    &
     &                             X4 * log(clat) ) / g + constant
            do i = 1-halo_i, row_length+halo_i
              exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /     &
     &                                    exp ( BV_squared_over_g *     &
     &                           (( 2.0 * r_theta_levels(i,j,k-1) -     &
     &                      r_rho_levels(i,j,k-1)) - Earth_radius +     &
     &                                                       temp))
            end do  !  i = 1-halo_i, row_length+halo_i
          endif    !  lat outside jet bounds
        end do  ! j = 1-halo_j, rows+halo_j

       endif  ! r_plane  <   -89.0    ! Option 2 or 3


      else     ! variable Coriolis, u = constant  Option 4

!   v = 0 and u = const. , f =2omega sinlat
        if(me  ==  0)then
          print*,' variable Coriolis, u = constant '
          print*,' Full vertical Coriolis term used '
        Endif   !(me  ==  0)

!      corterm = u*2omega*Earth_radius g
        corterm = u_in(1) * two_omega * Earth_radius / g
        If (L_regular) Then
          do j = 1-halo_j, rows+halo_j
            gj = l_datastart(2) + j - 1
            coslat(j) = cos(Base_phi + gj * delta_phi)
            y(j) = gj * delta_y
          end do
        Else  ! variable resolution
          do j = 1-halo_j, rows+halo_j
            coslat(j) = cos( phi_p(1,j) )
            y(j) = ( phi_p(1,j)- base_phi ) * Earth_radius
          end do
        Endif  !  L_regular

        do k = 1, model_levels
          do j = 1-halo_j, rows+halo_j
            If (coslat(j)  >   tolerance) then
              temp = u_in(1)* u_in(1) * log(coslat(j)) / g
            else   ! coslat(j)  <   tolerance
              temp = 0.0
            endif  ! coslat(j)  >   tolerance
            do i = 1-halo_i, row_length+halo_i
              If (coslat(j)  >   cos_ramp_end) then
                theta(i,j,k) = theta_surface * exp ( BV_squared_over_g *&
     &             (r_theta_levels(i,j,k) - Earth_radius +              &
     &              corterm * coslat(j) - temp ) )
                exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /   &
     &                            exp ( BV_squared_over_g *             &
     &                 (r_rho_levels(i,j,k) - Earth_radius +            &
     &                            corterm * coslat(j) - temp ) )
              else
                theta(i,j,k) = theta_surface * exp ( BV_squared_over_g *&
     &            (r_theta_levels(i,j,k) - Earth_radius))
                exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /   &
     &                                    exp ( BV_squared_over_g *     &
     &            (r_rho_levels(i,j,k) - Earth_radius))
              endIf !cos_lat  <=  cos_ramp_end
            end do
          end do
        end do
        k = model_levels + 1
        do j = 1-halo_j, rows+halo_j
          If (coslat(j)  >   tolerance) then
            temp = u_in(1)* u_in(1) * log(coslat(j)) / g
          else   ! coslat(j)  <   tolerance
            temp = 0.0
          endif  ! coslat(j)  >   tolerance
          do i = 1-halo_i, row_length+halo_i
            If (coslat(j) >   cos_ramp_end) then
              exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /     &
     &                                    exp ( BV_squared_over_g *     &
     &                            ( 2.0 * r_theta_levels(i,j,k-1) -     &
     &                       r_rho_levels(i,j,k-1) - Earth_radius +     &
     &                              corterm * coslat(j) - temp ) )
            else
              exner_rho_levels(i,j,k) = 1.0 - biggterm + biggterm /     &
     &                                    exp ( BV_squared_over_g *     &
     &                            ( 2.0 * r_theta_levels(i,j,k-1) -     &
     &                    r_rho_levels(i,j,k-1) - Earth_radius))
            endif !  cos_lat(j)  <=  cos_ramp_end
          end do
        end do

      endif     !  f_plane  >   -89.0

      endif !  L_trivial_trigs


      do k = 1, model_levels
        theta_ref(k) = theta_surface *                                  &
     &                     exp ( BV_squared_over_g *                    &
     &                     eta_theta_levels(k) * height_domain )
        exner_ref(k) = 1.0 - biggterm + biggterm /                      &
     &                     exp ( BV_squared_over_g *                    &
     &                     eta_rho_levels(k) * height_domain )
      end do

      k = model_levels + 1
      exner_ref(k) = 1.0 - biggterm + biggterm /                        &
     &            exp ( BV_squared_over_g * height_domain *             &
     &             (2.0 * eta_theta_levels(k-1) - eta_rho_levels(k-1)))

      IF (lhook) CALL dr_hook('IDL_PR_BALANCE',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE IDL_pr_balance
!  End subroutine IDL_pr_balance
