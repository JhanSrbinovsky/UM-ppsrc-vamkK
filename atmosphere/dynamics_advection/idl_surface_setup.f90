! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine IDL_Surface_setup

      SUBROUTINE IDL_Surface_setup(                                     &
                            row_length, rows, model_levels,             &
                            global_row_length, global_rows,             &
                            halo_i, halo_j,                             &
                            me, n_proc, at_extremity, model_domain,     &
                            l_datastart, all_proc_group,                &
                            delta_lambda, delta_phi, Base_phi,          &
                            n_rows, base_lambda,                        &
!  VarRes Grid Spacing
                            lambda_p, phi_p, lambda_u, phi_v,           &
                            lambda_p_end, phi_p_end,                    &
                            L_regular,                                  &
                            delta_x, delta_y,                           &
                            orography, orog_haloes,                     &
                            L_fix_orog_hgt_lbc, orog_hgt_lbc, rimwidth, &
                            surface_type,                               &
                            h_o, lambda_fraction, phi_fraction,         &
                            half_width_x, half_width_y,                 &
                            xp, yp, Witch_power,                        &
                            L_code_test)

! Purpose:
!          Sets up surface
!          surface_number defines options
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Advection
!
! Code Description:
!   Language: FORTRAN 77 + CRAY extensions
!   This code is written to UMDP3 programming standards.

      USE dynamics_grid_mod, ONLY: l_vatpoles

      USE earth_constants_mod, ONLY: earth_radius

      USE conversions_mod, ONLY: pi
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE UM_ParParams
      USE PrintStatus_mod
      IMPLICIT NONE

      INTEGER                                                           &
        row_length                                                      &
                         ! number of points on a processor row
      , rows                                                            &
                         ! number of rows in a processor theta field
      , model_levels                                                    &
                         ! number of model levels
      , halo_i                                                          &
                             ! Size of halo in i direction.
      , halo_j                                                          &
                             ! Size of halo in j direction.
      , n_rows                                                          &
                         ! number of rows in a processor v field
      , rimwidth         ! IN : Size of boundary region

      REAL                                                              &
        orography(row_length,rows)                                      &
      , orog_haloes(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j)   &
!    orog_haloes is r_theta_levels(0)
      , h_o                                                             &
                              ! mountain height
      , orog_hgt_lbc                                                    &
                              ! Height of orog in lbc zone
      , lambda_fraction                                                 &
                              ! hill position as fraction of domain
      , phi_fraction                                                    &
                              ! hill position as fraction of domain
      , half_width_x                                                    &
                             ! hill half-width East-West
      , half_width_y                                                    &
                             ! hill half-width North-South
      , xp                                                              &
                         ! Plateau EW half-width  (surface_type 4)
      , yp                                                              &
                         ! Plateau NS half-width  (surface_type 4)
      , Witch_power                                                     &
                         ! Exponent power in Witch definition
      , delta_x, delta_y

       INTEGER                                                          &
        model_domain                                                    &
                         ! holds integer code for model domain
      , global_row_length                                               &
                                ! number of points on a row
      , global_rows                                                     &
                                ! number of rows in a theta field
      , surface_type


      INTEGER                                                           &
        me                                                              &
                   ! My processor number
      , n_proc                                                          &
                   ! Total number of processors
      , all_proc_group                                                  &
                       ! Group identifier for all processors.
      , l_datastart(2)       ! First gridpoints held by this processor

      LOGICAL                                                           &
        at_extremity(4)                                                 &
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
      , L_code_test                                                     &
                         ! User switch
      , L_fix_orog_hgt_lbc ! Switch to fix orog hgt in lbc zone


      REAL                                                              &
           ! horizontal co-ordinate information
        delta_lambda                                                    &
      , delta_phi                                                       &
      , Base_phi                                                        &
      , Base_lambda                                                     &
      , lambda_p_end                                                    &
      , phi_p_end

      REAL                                                              &
           ! VarRes grid information
        lambda_p(1-halo_i:row_length+halo_i)                            &
      , lambda_u(1-halo_i:row_length+halo_i)                            &
      , phi_p    ( 1-halo_i : row_length + halo_i                       &
      ,            1-halo_j : rows + halo_j )                           &
      , phi_v    ( 1-halo_i : row_length + halo_i                       &
      ,            1-halo_j : n_rows+halo_j )

      LOGICAL, INTENT(IN) :: L_regular


! local variables

      INTEGER                                                           &
        i, j, gi ,gj                                                    &
      , haloi, haloj                                                    &
                      ! Local haloes - different for LAM/global
      , info, power

      REAL                                                              &
        x, y                                                            &
      , xo, yo                                                          &
      , width_x, width_y                                                &
      , twopi                                                           &
      , rad_to_deg                                                      &
      , lambda_o, phi_o                                                 &
      , lambda_deg, phi_deg                                             &
      , h_agnesi                                                        &
      , half_circumf                                                    &
      , dist_x                                                          &
      , dist_y                                                          &
      , h_max

      REAL :: offset

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle      


! Description: COMDECK containing surface types
!  for use in idealised problems
!

      INTEGER, PARAMETER :: surface_zero=0
      INTEGER, PARAMETER :: surface_ellipse=1
      INTEGER, PARAMETER :: surface_ridge=2
      INTEGER, PARAMETER :: surface_plateau=3
      INTEGER, PARAMETER :: surface_massif=4
      INTEGER, PARAMETER :: surface_mask=5
      INTEGER, PARAMETER :: surface_gauss=6
      INTEGER, PARAMETER :: surface_ridge_series=7
! ENDGAME-only parameters
      INTEGER, PARAMETER :: surface_schar_ridge=8
      INTEGER, PARAMETER :: surface_baroclinic=9
! End of ENDGAME-only parameters
      INTEGER, PARAMETER :: surface_dump=10

! No External routines

! ----------------------------------------------------------------------
! Section 1. Set domain information and check settings
! ----------------------------------------------------------------------
      IF (lhook) CALL dr_hook('IDL_SURFACE_SETUP',zhook_in,zhook_handle)

      rad_to_deg = 180.0 / Pi
      half_circumf = Earth_radius * Pi  !no need to vary with latitude
      IF (Witch_power  <   1.1) THEN
        power = NINT(Witch_power)
      ELSE
        power = 2
      END IF !   Witch_power  <   1.1


      delta_x = Earth_radius * delta_lambda
      delta_y = Earth_radius * delta_phi

      IF (l_vatpoles) THEN
         offset = 0.5
      ELSE
         offset = 1.0
      END IF ! vatpoles

!  xo is equator value -
!  latitude accounted for in orog_haloes calculation
      IF (model_domain  ==  mt_global) THEN
!  For global domains, do not do haloes (apply swap bounds at end)
        haloi = 0
        haloj = 0
        lambda_o =  lambda_fraction * 2.0 * Pi
        phi_o = phi_fraction * Pi
        xo = Earth_radius * lambda_o
        yo = Earth_radius * phi_o
        phi_o = phi_o - 0.5 * Pi
        lambda_deg = lambda_o * rad_to_deg
        phi_deg = phi_o * rad_to_deg
      ELSE IF (model_domain  ==  mt_cyclic_LAM .OR.                     &
             model_domain  ==  mt_bi_cyclic_LAM) THEN

!  For cyclic domains, do not do haloes (apply swap bounds at end)

        haloi = 0
        IF(model_domain  ==  mt_cyclic_LAM) THEN
           haloj=halo_j
        ELSE
           haloj=0
        END IF

        IF (L_regular) THEN
          xo = lambda_fraction *  global_row_length * delta_x
          IF (l_vatpoles) THEN
            yo = phi_fraction * global_rows * delta_y
            phi_o =  (phi_fraction - 0.5) * global_rows * delta_phi
          ELSE
            yo = phi_fraction * (global_rows - 1.0) * delta_y
            phi_o =  (phi_fraction - 0.5) * (global_rows - 1.0) * delta_phi
          END IF ! vatpoles
          lambda_o = lambda_fraction * global_row_length * delta_lambda
          lambda_deg = lambda_o * rad_to_deg
          phi_deg = phi_o * rad_to_deg
        ELSE
          lambda_o = lambda_fraction * (lambda_p_end-base_lambda)
          lambda_deg = lambda_o * rad_to_deg
          phi_o =  (phi_fraction - 0.5) * (phi_p_end-base_phi)
          phi_deg = phi_o * rad_to_deg
          xo = lambda_o * Earth_radius
          yo = phi_fraction *(phi_p_end-base_phi) * Earth_radius
        END IF     ! L_regular
      ELSE IF (model_domain  ==  mt_lam .AND. L_regular .AND.           &
                    lambda_fraction  <   0.49) THEN
!  For LAMs, fill haloes so extended halo regions have correct values
        haloi = halo_i
        haloj = halo_j
        xo = lambda_fraction *  global_row_length * delta_x
        IF (l_vatpoles) THEN
           yo = phi_fraction * (global_rows - 1.0) * delta_y
           phi_o =  (phi_fraction - 0.5) * global_rows * delta_phi
        ELSE
           yo = phi_fraction * (global_rows - 1.0) * delta_y
           phi_o =  (phi_fraction - 0.5) * (global_rows - 1.0) * delta_phi
        END IF ! vatpoles
        lambda_o = lambda_fraction * global_row_length * delta_lambda
        lambda_deg = lambda_o * rad_to_deg
        phi_deg = phi_o * rad_to_deg
      ELSE IF (model_domain  ==  mt_lam .AND. L_regular .AND.           &
                    lambda_fraction  >   0.49) THEN
!  For LAMs, fill haloes so extended halo regions have correct values
        haloi = halo_i
        haloj = halo_j
        IF (l_vatpoles) THEN
          xo = lambda_fraction * global_row_length * delta_x
          yo = phi_fraction * global_rows * delta_y
          lambda_o = lambda_fraction * global_row_length * delta_lambda
          phi_o =  (phi_fraction - 0.5) * global_rows * delta_phi
        ELSE
          xo = lambda_fraction *  (global_row_length - 1.0) * delta_x
          yo = phi_fraction * (global_rows - 1.0) * delta_y
          lambda_o = lambda_fraction * (global_row_length - 1.0) *        &
                                              delta_lambda
          phi_o =  (phi_fraction - 0.5) * (global_rows - 1.0) * delta_phi
        END IF ! vatpoles
        lambda_deg = lambda_o * rad_to_deg
        phi_deg = phi_o * rad_to_deg
      ELSE IF (model_domain  ==  mt_lam .AND. .NOT. L_regular) THEN
        haloi = halo_i
        haloj = halo_j
        lambda_o = lambda_fraction * (lambda_p_end-base_lambda)
        lambda_deg = lambda_o * rad_to_deg
        phi_o =  (phi_fraction - 0.5) * (phi_p_end-base_phi)
        phi_deg = phi_o * rad_to_deg
        xo = lambda_o * Earth_radius
        yo = phi_fraction *(phi_p_end-base_phi) * Earth_radius
      END IF

      IF (me  ==  0) THEN
        PRINT*,' Model_domain = ',model_domain
        IF ( model_domain  ==  mt_global ) THEN
          PRINT*,' Global model_domain  =  ',model_domain
        ELSE IF ( model_domain  ==  mt_LAM ) THEN
          PRINT*,' Limited-area model_domain  =  ',model_domain
        ELSE IF ( model_domain  ==  mt_cyclic_LAM ) THEN
          PRINT*,' Cyclic (East-West) Limited-area model_domain  =  '   &
                ,model_domain
        ELSE IF ( model_domain  ==  mt_bi_cyclic_LAM ) THEN
          PRINT*,' Bi-Cyclic Limited-area model_domain  =  '            &
                ,model_domain
        ELSE
          PRINT*,' DANGER model domain ',model_domain,' NOT SUPPORTED'
          PRINT*,' DANGER  MODEL WILL FAIL '
        END IF ! model_domain  ==  mt_global
      END IF ! (me == 0)

! ----------------------------------------------------------------------
! Section 2.  Insert mountain and set up orography arrays
! ----------------------------------------------------------------------

      IF (surface_type  ==  surface_zero ) THEN

        IF(me  ==  0) THEN
          PRINT*,' Flat surface_type = ',surface_type
          PRINT*,' Orography set to 0.0 everywhere '
        END IF ! (me == 0)

!  Make sure h_o = 0 for diagnostic printing
        h_o = 0.0

        DO j = 1-haloj, rows+haloj
          DO i = 1-haloi, row_length+haloi
            orog_haloes(i,j) = 0.0
          END DO
        END DO

      ELSE IF(surface_type  ==  surface_ellipse) THEN
        IF(me  ==  0) THEN
          PRINT*,'surface_type = ',surface_type,' ** Witch of Agnesi **'
          PRINT*,'Hill height = ',h_o,' metres. '
          PRINT*,'Centre at ',lambda_deg,' degrees East or '            &
          ,lambda_o,' radians, ',xo,' metres from Greenwich meridian'
          IF (phi_deg  >=  0.0) THEN
            PRINT*,'Centre at ',phi_deg,' degrees North or '            &
                  ,phi_o,' radians, ',yo,' metres from South Pole '
          ELSE
            PRINT*,'Centre at ',phi_deg,' degrees South or '            &
                  ,phi_o,' radians, ',yo,' metres from South Pole '
          END IF    !(phi_deg  >=  0.0)
          IF(half_width_x  ==  half_width_y)THEN
            PRINT*,' Circular hill half-width = ',half_width_x,' metres'
          ELSE
            PRINT*,'Elliptical hill half_width_x = ',half_width_x,      &
                   ' metres',' half_width_y = ',half_width_y,' metres'
          END IF     !(half_width_x  ==  half_width_y)
        END IF ! (me == 0)

        IF (L_regular) THEN
          DO j = 1-haloj, rows+haloj
            DO i = 1-haloi, row_length+haloi
              gi = l_datastart(1) + i - 1
              gj = l_datastart(2) + j - 1
              x = (gi-offset) * delta_x
              y = (gj-offset) * delta_y
              dist_x = ABS(x-xo)
!  Spherical term cos_latitude not needed in x since it cancels with
!   same term required for half_width_x
            IF(dist_x  >   half_circumf) THEN
              dist_x = 2.0 * half_circumf - x + xo
            END IF     ! dist_x  >   half_circumf

              IF(power  ==  1)THEN
                orog_haloes(i,j) = h_o /(1. +                             &
                     (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2)
              ELSE
                orog_haloes(i,j) = h_o /(1. +                             &
                    (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2)  &
                     ** Witch_power
              END IF   !(power  ==  1)
            END DO
          END DO
        ELSE
          DO j = 1-haloj, rows+haloj
            DO i = 1-haloi, row_length+haloi
              x = (lambda_p(i)-base_lambda) * Earth_radius
              y = (phi_p(i,j) - base_phi) * Earth_radius
              dist_x = ABS(x-xo)

!  Spherical term cos_latitude not needed in x since it cancels with
!   same term required for half_width_x

              IF(dist_x  >   half_circumf) THEN
                dist_x = 2.0 * half_circumf - x + xo
              END IF     ! dist_x  >   half_circumf

              IF(power  ==  1)THEN
                orog_haloes(i,j) = h_o /(1. +                             &
                   (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2)
              ELSE
                orog_haloes(i,j) = h_o /(1. +                             &
                  (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2)    &
                   ** Witch_power
              END IF   !(power  ==  1)
            END DO
          END DO
        END IF   ! L_regular

      ELSE IF(surface_type  ==  surface_ridge) THEN
!
! code to create a cos**2 ridge with amplitude 4*half_width_x
! Sets to 0 outside region
!
        width_x = 2.0 * half_width_x
        width_y = 2.0 * half_width_y

        IF(me  ==  0)THEN
          PRINT*,' surface_type = ',surface_type
          PRINT*,' North-South ridge'
          PRINT*,' Ridge height = ',h_o,' metres. Centre at '           &
                ,xo,' metres'
          PRINT*,' Ridge half-width = ',half_width_x,' metres'
        END IF ! (me == 0)

        IF (L_regular) THEN
          DO j = 1-haloj, rows+haloj
            DO i = 1-haloi, row_length+haloi
               gi = l_datastart(1) + i - 1
               x =  (gi-offset) * delta_x
              IF( ABS(x-xo)  <=  width_x ) THEN
                orog_haloes(i,j) = h_o * COS(Pi*(x-xo)/(2.*width_x))**2
              ELSE
                orog_haloes(i,j) = 0.0
              END IF
            END DO
          END DO
        ELSE
          DO j = 1-haloj, rows+haloj
            DO i = 1-haloi, row_length+haloi
              x = (lambda_p(i)-base_lambda) * Earth_radius
              IF( ABS(x-xo)  <=  width_x ) THEN
                orog_haloes(i,j) = h_o * COS(Pi*(x-xo)/(2.*width_x))**2
              ELSE
                orog_haloes(i,j) = 0.0
              END IF
            END DO
          END DO
        END IF   ! L_regular

      ELSE IF(surface_type  ==  surface_ridge_series) THEN

        width_x = 2.0 * half_width_x
        width_y = 2.0 * half_width_y
        twopi = 2.0 * Pi
        IF(me  ==  0)THEN
          PRINT*,' surface_type = ',surface_type
          PRINT*,' North-South ridge series'
          PRINT*,' Ridge height = ',h_o,' metres. Centre at '           &
                ,xo,' metres'
          PRINT*,' Ridge half-width = ',half_width_x,' metres'
          gi = l_datastart(1)

          IF (L_regular) THEN
            x = (gi-offset) * delta_x
          ELSE
            x = (lambda_p(1)-base_lambda) * Earth_radius
          END IF
          orog_haloes(1,1) = h_o * COS(twopi * (x-xo) / width_x)
          PRINT*,' First point gi value = ',gi,' First point x = ',x
          PRINT*,' orog_haloes(1,1) = ',orog_haloes(1,1)
        END IF ! (me == 0)

        IF (L_regular) THEN
          DO j = 1-haloj, rows+haloj
            DO i = 1-haloi, row_length+haloi
              gi = l_datastart(1) + i - 1
              x =  (gi-offset)*delta_x
              IF( ABS(x-xo)  <=  half_width_x ) THEN
                orog_haloes(i,j) = h_o * COS(twopi*(x-xo)/width_x)
              END IF
            END DO
          END DO
        ELSE
          DO j = 1-haloj, rows+haloj
            DO i = 1-haloi, row_length+haloi
              x = (lambda_p(i)-base_lambda) * Earth_radius
              IF( ABS(x-xo)  <=  half_width_x ) THEN
                orog_haloes(i,j) = h_o * COS(twopi*(x-xo)/width_x)
              END IF
            END DO
          END DO
        END IF   ! L_regular

      ELSE IF(surface_type  ==  surface_plateau )THEN
        IF(me  ==  0)THEN
          PRINT*,'surface_type = ',surface_type,                        &
           ' ** FLATTENED  Witch of Agnesi **'
          PRINT*,'** Constant height from centre to half-width '
          PRINT*,'Plateau height = ',h_o,' metres. '
          PRINT*,'Centre at ',lambda_deg,' degrees East or '            &
            ,lambda_o,' radians, ',xo,' metres from Greenwich meridian'
          IF (phi_deg  >=  0.0) THEN
            PRINT*,'Centre at ',phi_deg,' degrees North or '            &
            ,phi_o,' radians, ',yo,' metres from South Pole '
          ELSE
            PRINT*,'Centre at ',phi_deg,' degrees South or '            &
            ,phi_o,' radians, ',yo,' metres from South Pole '
          END IF    !(phi_deg  >=  0.0)

          IF(half_width_x  ==  half_width_y)THEN
            PRINT*,' Circular hill half-width = ',half_width_x,' metres'
          ELSE
         PRINT*,'Elliptical hill half_width_x = ',half_width_x,' metres' &
                  , ' half_width_y = ',half_width_y,' metres'
          END IF     !(half_width_x  ==  half_width_y)
        END IF ! (me == 0)

        IF (L_regular) THEN
          h_agnesi = 2.0 * h_o
          DO j = 1-haloj, rows+haloj
            DO i = 1-haloi, row_length+haloi
              gi = l_datastart(1) + i - 1
              gj = l_datastart(2) + j - 1
              x = (gi-offset) * delta_x
              y = (gj-offset) * delta_y
              dist_x = ABS(x-xo)

!  Spherical term cos_latitude not needed in x since it cancels with
!   same term required for half_width_x
              IF(dist_x  >   half_circumf) THEN
                dist_x = 2.0 * half_circumf - x + xo
              END IF 

              IF(power  ==  1)THEN
                orog_haloes(i,j) = h_agnesi /(1. +                        &
                     (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2)
              ELSE
                orog_haloes(i,j) = h_agnesi /(1. +                        &
                     (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2) &
                     ** Witch_power
              END IF 

              IF(orog_haloes(i,j)  >   h_o )THEN
                orog_haloes(i,j) = h_o
              END IF
            END DO
          END DO
        ELSE
          h_agnesi = 2.0 * h_o
          DO j = 1-haloj, rows+haloj
            DO i = 1-haloi, row_length+haloi
              x = (lambda_p(i)-base_lambda) * Earth_radius
              y = (phi_p(i,j) - base_phi) * Earth_radius
              dist_x = ABS(x-xo)
!  Spherical term cos_latitude not needed in x since it cancels with
!   same term required for half_width_x

              IF(dist_x  >   half_circumf) THEN
                dist_x = 2.0 * half_circumf - x + xo
              END IF

              IF(power  ==  1)THEN
                orog_haloes(i,j) = h_agnesi /(1. +                        &
                     (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2)
              ELSE
                orog_haloes(i,j) = h_agnesi /(1. +                        &
                     (dist_x/half_width_x)**2 + ((y-yo)/half_width_y)**2) &
                     ** Witch_power
              END IF

              IF(orog_haloes(i,j)  >   h_o )THEN
                orog_haloes(i,j) = h_o
              END IF
            END DO
          END DO
        END IF   ! L_regular

       ELSE IF(surface_type  ==  surface_massif)THEN
       IF(me  ==  0)THEN
          PRINT*,'surface_type = ',surface_type
          PRINT*,' ** Massif plateau based upon 4 Witches of Agnesi **'
          PRINT*,'Hill height = ',h_o,' metres. '
          PRINT*,'Centre at ',lambda_deg,' degrees East or '            &
            ,lambda_o,' radians, ',xo,' metres from Greenwich meridian'

          IF (phi_deg  >=  0.0) THEN
            PRINT*,'Centre at ',phi_deg,' degrees North or '            &
                  ,phi_o,' radians, ',yo,' metres from South Pole '
          ELSE
            PRINT*,'Centre at ',phi_deg,' degrees South or '            &
                  ,phi_o,' radians, ',yo,' metres from South Pole '
          END IF

          IF(half_width_x  ==  half_width_y)THEN
            PRINT*,' Circular hill half-width = ',half_width_x,' metres'
          ELSE
            PRINT*,'Elliptical hill half_width_x = ',half_width_x,      &
                   ' metres', ' half_width_y = ',half_width_y,' metres'
          END IF     !(half_width_x  ==  half_width_y)
       END IF

! There are nine sub-regions to set count 1-9 from bottom (S-Pole)
! and Eastwards
!               *            *
!           7   *     8      *    9
!               *            *
!     *************************************
!               * 5          * |
!           4   *   (xo,yo)  * 2yp   6
!               *            * |
!     *************************************
!               *<-- 2xp  -->*
!           1   *            *  3
!               *     2      *

        DO j = 1-haloj, rows+haloj
          DO i = 1-haloi, row_length+haloi
            IF (L_regular) THEN
              gi = l_datastart(1) + i - 1
              gj = l_datastart(2) + j - 1
              x = (gi-offset) * delta_x
              y = (gj-offset) * delta_y
            ELSE
              x = (lambda_p(i) - base_lambda) * Earth_radius
              y = (phi_p(i,j) - base_phi) * Earth_radius
            END IF
            IF(x  <   xo - xp) THEN

              IF( y  <   yo - yp)THEN       ! region 1
                dist_x =  ABS(x-xo+xp)
                dist_y =  ABS(y-yo+yp)
              ELSE IF( y  >   yo + yp)THEN   ! region 7
                dist_x =  ABS(x-xo+xp)
                dist_y =  ABS(y-yo-yp)
              ELSE                         ! region 4
                dist_x =  ABS(x-xo+xp)
                dist_y =  0.0
              END IF

            ELSE IF(x  >   xo + xp) THEN
              IF( y  <   yo - yp)THEN       ! region 3
                dist_x =  ABS(x-xo-xp)
                dist_y =  ABS(y-yo+yp)
              ELSE IF( y  >   yo + yp)THEN   ! region 9
                dist_x =  ABS(x-xo-xp)
                dist_y =  ABS(y-yo-yp)
              ELSE                         ! region 6
                dist_x =  ABS(x-xo-xp)
                dist_y =  0.0
              END IF

            ELSE
              IF( y  <   yo - yp)THEN       ! region 2
                dist_x =  0.0
                dist_y =  ABS(y-yo+yp)
              ELSE IF( y  >   yo + yp)THEN   ! region 8
                dist_x =  0.0
                dist_y =  ABS(y-yo-yp)
              ELSE                         ! region 5
                dist_x =  0.0
                dist_y =  0.0
              END IF
            END IF 

!  Spherical term cos_latitude not needed in x since it cancels with
!   same term required for half_width_x

            IF(dist_x  >   half_circumf ) THEN
              IF(x  <   xo - xp) THEN
                dist_x = 2.0 * half_circumf - x + xo - xp
              ELSE IF(x  >   xo + xp) THEN
                dist_x = 2.0 * half_circumf - x + xo + xp
              END IF
            END IF

            IF(power  ==  1)THEN
              orog_haloes(i,j) = h_o /(1. +                             &
                  (dist_x/half_width_x)**2 + (dist_y/half_width_y)**2)
            ELSE
              orog_haloes(i,j) = h_o /(1. +                             &
                  (dist_x/half_width_x)**2 + (dist_y/half_width_y)**2)  &
                   ** Witch_power
            END IF
          END DO
        END DO

      ELSE IF(surface_type  ==  surface_mask)THEN

       IF(me  ==  0)THEN
        PRINT*,'surface_type = ',surface_type,' ** Real orography **'
        PRINT*,'within box East-West half-width ',half_width_x,' metres'&
      , ' North-South half-width ',half_width_y,' metres'
        PRINT*,'Centre at ',lambda_deg,' degrees East or '              &
       ,lambda_o,' radians, '
        PRINT*,xo,' metres from Greenwich meridian'

         IF (phi_deg  >=  0.0) THEN
           PRINT*,'Centre at ',phi_deg,' degrees North or '             &
                              ,phi_o,' radians '
           PRINT*,yo,' metres from South Pole '
         ELSE
           PRINT*,'Centre at ',phi_deg,' degrees South or '             &
                              ,phi_o,' radians, '
           PRINT*,yo,' metres from South Pole '
         END IF
       END IF

!  Zero r_theta_levels(i,j,0)
      DO j = 1-halo_j, rows+halo_j
        DO i = 1-halo_i, row_length+halo_i
          orog_haloes(i,j)  = 0.0
        END DO
      END DO
!  Copy orography into r_theta_levels(i,j,0)
      DO j = 1, rows
        DO i = 1, row_length
          orog_haloes(i,j)  = orography(i,j)
        END DO
      END DO

! DEPENDS ON: swap_bounds
      CALL SWAP_BOUNDS(orog_haloes, row_length, rows, 1,                &
                       halo_i, halo_j, fld_type_p,.FALSE.)


!   Reset h_o and use to store max orography
        h_max = 0.0

        DO j = 1-haloj, rows+haloj
          DO i = 1-haloi, row_length+haloi
            IF (L_regular) THEN
              gi = l_datastart(1) + i - 1
              gj = l_datastart(2) + j - 1
              x = (gi-offset) * delta_x
              y = (gj-offset) * delta_y
            ELSE
              x = (lambda_p(i) - base_lambda) * Earth_radius
              y = (phi_p(i,j) - base_phi) * Earth_radius
            END IF

            dist_x = ABS(x-xo)
            dist_y = ABS(y-yo)
            IF(dist_x  >   half_width_x .OR.                            &
               dist_y  >   half_width_y) THEN
              orog_haloes(i,j) = 0.0
            END IF

!  Set h_o to max orography found
            IF (orog_haloes(i,j)  >   h_max)THEN
              h_max = orog_haloes(i,j)
            END IF
          END DO
        END DO

! find max  over all processors
        CALL gc_rmax(1, n_proc, info, h_max)
! Put max in h_o so that z_orog_print can calculate grid over max h_o
        h_o = h_max

      ELSE IF(surface_type  ==  surface_gauss)THEN

        IF(me  ==  0)THEN
          PRINT*,' surface_type = ',surface_type,' ** Gaussian **'
          PRINT*,' Hill height = ',h_o,' metres. East-West centre at '  &
                ,xo,' metres, North-South centre at ',yo,' metres'
          IF(half_width_x  ==  half_width_y)THEN
            PRINT*,' Circular hill half-width = ',half_width_x,' metres'
          ELSE
           PRINT*,'Elliptical hill half_width_x = ',                    &
                  half_width_x,' metres',                               &
                 ' half_width_y = ',half_width_y,' metres'
          END IF
        END IF

        DO j = 1-haloj, rows+haloj
          DO i = 1-haloi, row_length+haloi
            IF (L_regular) THEN
              gi = l_datastart(1) + i - 1
              gj = l_datastart(2) + j - 1
              x = (gi-offset) * delta_x
              y = (gj-offset) * delta_y
            ELSE
              x = (lambda_p(i) - base_lambda) * Earth_radius
              y = (phi_p(i,j) - base_phi) * Earth_radius
            END IF

            dist_x = ABS(x-xo)
            IF(dist_x  >   half_circumf) THEN
              dist_x = 2.0 * half_circumf - x + xo
            END IF     ! dist_x  >   half_circumf
            orog_haloes(i,j) = h_o * exp( -1.0 * (                      &
                         ( dist_x/half_width_x)**2. +                   &
                               ((y-yo)/half_width_y)**2.))
          END DO
        END DO
        
!-----------------------------------------------------------------------
! Schar ridge surface
!-----------------------------------------------------------------------

      ELSE IF(surface_type  ==  surface_schar_ridge) THEN

        IF(me  ==  0 .AND. PrintStatus >= PrStatus_Normal)THEN
          PRINT*, ' surface_type = ',surface_type
          PRINT*, ' North-South Schar ridge'
          PRINT*, ' ridge height = ',h_o, ' metres. Centre at ',        &
                    360.0*lambda_fraction,' degrees'
          PRINT*, ' ridge half-width = ',half_width_x, ' radians'
          PRINT*,' ' 
        END IF ! (me == 0)

! theta-point orog
        DO j=1-halo_j, rows+halo_j
          DO i=1-halo_i, row_length+halo_i
            IF (L_regular) THEN
              gi = l_datastart(1) + i - 1
              gj = l_datastart(2) + j - 1
              x = (gi-offset) * delta_x
              y = (gj-offset) * delta_y
            ELSE
              x = (lambda_p(i) - base_lambda) * Earth_radius
              y = (phi_p(i,j) - base_phi) * Earth_radius
            END IF        
            x = x - 2.0*pi*lambda_fraction
            
            orog_haloes(i,j) = h_o * COS(pi*x/(half_width_x))**2          &
                           *EXP(-(x/half_width_y)**2) 
          END DO
        END DO
        
      ELSE IF(surface_type  ==  surface_dump)THEN

!  Zero r_theta_levels(i,j,0)
        DO j = 1-halo_j, rows+halo_j
          DO i = 1-halo_i, row_length+halo_i
            orog_haloes(i,j)  = 0.0
          END DO
        END DO
!  Copy orography into r_theta_levels(i,j,0)
        DO j = 1, rows
          DO i = 1, row_length
            orog_haloes(i,j)  = orography(i,j)
          END DO
        END DO

! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(orog_haloes, row_length, rows, 1,              &
     &                 halo_i, halo_j, fld_type_p,.FALSE.)

        !--------------------------------------------------------
        !
        !          Set orography in the lbc zone
        !
        !--------------------------------------------------------

        ! If (L_fix_orog_hgt_lbc) is true and this is a limited
        ! area model domain, then set the height of
        ! the orography in the lateral boundary zone to a constant
        ! (external halo and rimwidth) and linearly interpolate
        ! in a zone in the interior to the input orography
        ! for a smooth transition (currently hardwired to 4 points).

        IF (L_fix_orog_hgt_lbc .AND.                                    &
            model_domain == mt_cyclic_LAM) THEN

          IF (me == 0) THEN
            PRINT*, ' '
            PRINT*, ' L_fix_orog_hgt_lbc = .True.'
            PRINT*, ' Setting constant height orography'
            PRINT*, ' in the lbc zone.',orog_hgt_lbc,' (m)'
          END IF

          IF (at_extremity(PSouth)) THEN
            DO j = 1-halo_j,  rimwidth
              DO i = 1-halo_i, row_length+halo_i
                orog_haloes(i,j) = orog_hgt_lbc
              END DO
            END DO
          END IF

          IF (at_extremity(PNorth)) THEN
            DO j = 1-halo_j,  rimwidth
              DO i = 1-halo_i, row_length+halo_i
                orog_haloes(i,rows-j+1) = orog_hgt_lbc
              END DO
            END DO
          END IF

          IF (at_extremity(PWest)) THEN
            DO j = 1, rows
              DO i = 1-halo_i, rimwidth
                orog_haloes(i,j) = orog_hgt_lbc
              END DO
            END DO
          END IF

          IF (at_extremity(PEast)) THEN
            DO j = 1, rows
              DO i = 1-halo_i, rimwidth
                orog_haloes(rows-i+1, j) = orog_hgt_lbc
              END DO
            END DO
          END IF

          IF (at_extremity(PSouth)) THEN
            DO j = rimwidth + 1,  rimwidth + 3
              DO i =  rimwidth + 1, row_length - rimwidth
                dist_x = j - rimwidth
                h_max = orography(i, rimwidth + 4) - orog_hgt_lbc
                orog_haloes(i,j) = orog_hgt_lbc + 0.25 * dist_x * h_max
              END DO
            END DO
          END IF

          IF (at_extremity(PNorth)) THEN
            DO j = rimwidth + 1,  rimwidth + 3
              DO i =  rimwidth + 1, row_length - rimwidth
                dist_x = j - rimwidth
                gj = rows-j+1
                y = orography(i, rows - rimwidth - 3) - orog_hgt_lbc
                orog_haloes(i,gj) = orog_hgt_lbc + 0.25 * dist_x * y
              END DO
            END DO
          END IF

          IF (at_extremity(PWest)) THEN
            DO j = rimwidth + 4, rows - rimwidth - 3
              DO i =  rimwidth + 1,  rimwidth + 3
                dist_x = i - rimwidth
                h_max = orography(rimwidth + 4 , j) - orog_hgt_lbc
                orog_haloes(i,j) = orog_hgt_lbc + 0.25 * dist_x * h_max
              END DO
            END DO
          END IF

          IF (at_extremity(PEast)) THEN
            DO j = rimwidth + 4, rows - rimwidth - 3
              DO i =  rimwidth + 1,  rimwidth + 3
                dist_x = i - rimwidth
                gi = row_length-i+1
                y = orography(row_length-rimwidth-3, j) - orog_hgt_lbc
                orog_haloes(gi, j)= orog_hgt_lbc + 0.25 * dist_x * y
              END DO
            END DO
          END IF

          IF (at_extremity(PSouth) .AND. at_extremity(PEast)) THEN
            h_max = orography(rimwidth+4, rimwidth+4) - orog_hgt_lbc
            DO i = 1, 3
              gi = rimwidth+i
              dist_x = i
              DO j = gi, rimwidth+4, 1
                orog_haloes(j, gi) = orog_hgt_lbc + 0.25 * dist_x*h_max
                orog_haloes(gi, j) = orog_hgt_lbc + 0.25 * dist_x*h_max
              END DO
            END DO
          END IF

          IF (at_extremity(PSouth) .AND. at_extremity(PWest)) THEN
            h_max = orography(row_length-rimwidth-3, rimwidth + 4)      &
                            - orog_hgt_lbc
            DO i = 1, 3
              gi = row_length-rimwidth+1-i
              gj = rimwidth+i
              dist_x = i
              DO j = gi, row_length-rimwidth-1,-1
                orog_haloes(j, gj) = orog_hgt_lbc + 0.25 * dist_x*h_max
              END DO
              DO j = gj, rimwidth+4 ,1
                orog_haloes(gi, j) = orog_hgt_lbc + 0.25 * dist_x*h_max
              END DO
            END DO
          END IF

          IF (at_extremity(PNorth) .AND. at_extremity(PWest)) THEN
            h_max = orography(row_length-rimwidth-3, rows-rimwidth-3)   &
                            - orog_hgt_lbc
            DO i = 1, 3
              gi = row_length-rimwidth+1-i
              gj = rows-rimwidth+1-i
              dist_x = i
              DO j = gi, row_length-rimwidth-1,-1
                orog_haloes(j, gj) = orog_hgt_lbc + 0.25 * dist_x*h_max
              END DO
              DO j = gj, rows-rimwidth-1,-1
                orog_haloes(gi, j) = orog_hgt_lbc + 0.25 * dist_x*h_max
              END DO
            END DO
          END IF

          IF (at_extremity(PNorth) .AND. at_extremity(PEast)) THEN
            h_max = orography(rimwidth+4,rows-rimwidth-3)-orog_hgt_lbc
            DO i = 1, 3
              gi = rimwidth +i
              gj = rows-rimwidth+1-i
              dist_x = i
              DO j = gi, rimwidth+4, 1
                orog_haloes(j, gj) = orog_hgt_lbc + 0.25 * dist_x*h_max
              END DO
              DO j = gj, rows-rimwidth-1,-1
                orog_haloes(gi, j) = orog_hgt_lbc + 0.25 * dist_x*h_max
              END DO
            END DO
          END IF

        END IF

!   Reset h_o and use to store max orography
        h_max = 0.0

        DO j=1,rows
          DO i=1,row_length
            IF (L_regular) THEN
              gi = l_datastart(1) + i - 1
              gj = l_datastart(2) + j - 1
              x = (gi-offset) * delta_x
              y = (gj-offset) * delta_y
            ELSE
              x = (lambda_p(i) - base_lambda) * Earth_radius
              y = (phi_p(i,j) - base_phi) * Earth_radius
            END IF

            dist_x = ABS(x-xo)
            dist_y = ABS(y-yo)
!  Set h_o to max orography found
            IF (orography(i,j)  >   h_max)THEN
              h_max = orography(i,j)
            END IF
          END DO
        END DO

! find max  over all processors
        CALL gc_rmax(1, n_proc, info, h_max)
! Put max in h_o so that z_orog_print can calculate grid over max h_o

        h_o = h_max
        IF(me  ==  0)THEN
          PRINT*,'surface_type = ',surface_type
          PRINT*,' ** Orography in input dump used everywhere **'
          PRINT*,'Maximum orographic height = ',h_o,' metres'
        END IF ! (me == 0)

      ELSE

        PRINT*,' surface_type = ',surface_type,' option unavailable '

      END IF        ! on surface_type

! DEPENDS ON: swap_bounds
        CALL SWAP_BOUNDS(orog_haloes, row_length, rows, 1,              &
                       halo_i, halo_j, fld_type_p,.FALSE.)

! Next series of loops adjust the mountain profile so as to ensure that
! it is zero along Northern and Southern LAM boundaries.
!  Not needed for ridge (type 2) nor input orography (type 10)

      IF(    surface_type /=  surface_ridge                             &
       .AND. surface_type /=  surface_dump                              &
       .AND. surface_type /=  surface_schar_ridge )THEN

        PRINT*,' Northern and Southern boundaries set to zero'

        IF (at_extremity(PSouth)) THEN
          DO j = 1, 2
            DO i = 1-haloi, row_length+haloi
              orog_haloes(i,j) = 0.0
            END DO
          END DO
        END IF

        IF (at_extremity(PNorth)) THEN
          DO j = rows-1, rows
            DO i = 1-haloi, row_length+haloi
              orog_haloes(i,j) = 0.0
            END DO
          END DO
        END IF

      END IF 

!  Copy orog_haloes into D1 array

      DO j = 1, rows
        DO i = 1, row_length
          orography(i,j) = orog_haloes(i,j)
        END DO
      END DO

      IF (lhook) CALL dr_hook('IDL_SURFACE_SETUP',zhook_out,zhook_handle)

      END SUBROUTINE IDL_Surface_setup

