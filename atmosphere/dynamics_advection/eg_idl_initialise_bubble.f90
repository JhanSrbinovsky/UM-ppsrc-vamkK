! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_idl_initialise_bubble_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_idl_initialise_bubble(                                  &
                      row_length, rows, halo_i, halo_j,               &
                      offx, offy, model_levels,                       &
                      delta_xi1, delta_xi2, base_xi1, base_xi2,       &
                      xi1_p, xi2_p, xi3_at_rho, xi3_at_theta,         &
                      intw_w2rho, intw_rho2w,                         &
                      me,                                             &
                      l_datastart, global_row_length, global_rows,    &
                      idl_max_num_bubbles, idl_bubble_width,          &
                      idl_bubble_height, idl_bubble_depth,            &
                      idl_bubble_xoffset, idl_bubble_yoffset,         &
                      dt_bubble, theta, rho, exner,                   &
                      exner_star,                                     &
                      t_surface, p_surface, l_cartesian,              &
                      bubble_option )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atmos_constants_mod, ONLY : r,kappa,p_zero,cp
USE earth_constants_mod, ONLY : earth_radius, g
USE atm_fields_bounds_mod
USE eg_idl_1d_profs_mod
USE horiz_grid_mod, ONLY : glob_xi1_p,glob_xi2_p
USE conversions_mod, ONLY: pi
USE timestep_mod, ONLY : timestep_number
IMPLICIT NONE
!
! Description:
!   Adds a 3D temperature anomaly (bubble) to the potential temperature
!   field (warm or cold). There are a number of options for different
!   shaped bubbles. There is also the option to saturate the bubble.
!  
!
! Method:
!   Sets up x,y,z functions containing the distance from the centre
!   of the bubble, then applies a function to generate the bubble
!   amplitude at all points, and adds this to the model field.
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle



! Array Dimensions
INTEGER, INTENT(IN) :: row_length   ! no. points on a row
INTEGER, INTENT(IN) :: rows         ! no. rows
INTEGER, INTENT(IN) :: halo_i       ! Size of halo x-direction
INTEGER, INTENT(IN) :: halo_j       ! Size of halo y-direction
INTEGER, INTENT(IN) :: offx         ! Small halo x-direction
INTEGER, INTENT(IN) :: offy         ! Small halo y-direction
INTEGER, INTENT(IN) :: model_levels ! number of model levels

! Horizontal co-ordinate information


REAL, INTENT(IN) ::                                                   &
  xi3_at_rho(1-halo_i:row_length+halo_i, 1-halo_j:rows+halo_j,        &
               model_levels),                                         &
  xi3_at_theta(1-halo_i:row_length+halo_i,                            &
                            1-halo_j:rows+halo_j,0:model_levels)


! Multi-processor array variables
INTEGER, INTENT(IN) :: me                ! My processor number
INTEGER, INTENT(IN) :: global_row_length ! global number of pts
INTEGER, INTENT(IN) :: global_rows       ! global number of rows
INTEGER, INTENT(IN) :: l_datastart(3)    ! Start row for pe


! Idealised options
INTEGER, INTENT(IN)  :: idl_max_num_bubbles
                    ! Number of bubbles: currently 1
REAL, INTENT(IN)  :: dt_bubble
                    ! Bubble maximum amplitude (K)
REAL, INTENT(IN)  :: idl_bubble_width
                    ! Essentially a scaling factor (m)
REAL, INTENT(IN)  :: idl_bubble_height
                    ! Height of bubble centre (m)
REAL, INTENT(IN)  :: idl_bubble_depth
                    ! Bubble's flat radius (m)
REAL, INTENT(IN)  :: idl_bubble_xoffset
         ! Bubble x-offset (normalised units: 0.5 = domain centre)
REAL, INTENT(IN)  :: idl_bubble_yoffset
         ! Bubble y-offset (normalised units: 0.5 = domain centre)         
INTEGER, INTENT(IN) :: bubble_option
         ! bubble choice (1) Guassian, (2) Sin in the vertical

REAL, INTENT(IN) ::     delta_xi1,                                    &
                        delta_xi2,                                    &
                        base_xi1,                                     &
                        base_xi2
                        
REAL, INTENT(IN) :: t_surface, p_surface                      

REAL, INTENT(IN) ::                                                   &
  xi1_p(1-halo_i:row_length+halo_i),                                  &
  xi2_p(1-halo_j:rows+halo_j) ,                                       &
  intw_w2rho(model_levels,2), intw_rho2w(model_levels,2)

! Data arrays
REAL, INTENT (INOUT) ::                                                  &
  exner(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1),         &
  exner_star(1-offx:row_length+offx, 1-offy:rows+offy)

REAL, INTENT (INOUT) ::                                               &
  theta(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
                      ! Potential Temperature (K)
  rho(1-offx:row_length+offx,1-offy:rows+offy,model_levels)  
  
LOGICAL, INTENT(IN)   :: l_cartesian

! Local Variables

INTEGER ::             i,                                             &
                       j,                                             &
                       k   ! Loop indices

INTEGER :: b         ! bubble number

REAL ::                 x0,                                           &
                        y0,                                           &
                        z0,                                           &
                        dist,                                         &
                        tmp,                                          &
                        t,                                            &
                        bubble_flat_radius,                           &
                        a0

REAL :: bubblefn(row_length,rows,0:model_levels)

REAL :: grav_1d(0:model_levels)

REAL :: exner_1d(0:model_levels+1)

REAL :: dtheta_dz1(3) 

REAL :: t_surf, p_surf 

INTEGER :: gaussian_bubble = 1
INTEGER :: sin_bubble = 2
!- End of header

! 1.0 Start of subroutine code: perform the calculation.

IF ( timestep_number /= 1 ) RETURN

IF (lhook) CALL dr_hook('EG_IDL_INITIALISE_BUBBLE',zhook_in,zhook_handle)


! ---------------------------------------------------------------------
! Set up background atmosphere 
dtheta_dz1 = 0.0
DO k = 0, model_levels
   grav_1d(k)    = g
END DO
t_surf = t_surface
p_surf = p_surface
a0 = earth_radius

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
! Adjust p_surf for non-zero orgraphy
!     t_surf = t_surface -6.0/1000.0*(xi3_at_theta(i,j,0)-a0)
      p_surf = p_surface*EXP(-(xi3_at_theta(i,j,0)-a0)                &
              *g/(r*t_surf))

      exner_star(i,j) = (p_surf/p_zero)**kappa

      CALL eg_idl_1d_profs(                                           &
                       exner_1d,                                      &
                       rho(i,j,:),                                    &
                       theta(i,j,:),                                  &
                       2,                                             &
                       t_surf, p_surf,                                &
                       xi1_p(i), xi2_p(j),                            &
                       intw_rho2w, intw_w2rho,                        &
                       xi3_at_theta(i,j,:), xi3_at_rho(i,j,:),        &
                       1.0, earth_radius, dtheta_dz1, grav_1d,        &
                       model_levels)
      DO k = 1,model_levels+1
        exner(i,j,k) = exner_1d(k)
      END DO
    END DO
  END DO
 
! WRITE(6,*) 'Initial atmosphere'
!   DO j = pdims%j_start, pdims%j_end
!     DO i = pdims%i_start, pdims%i_end
!       WRITE(6,fmt='(3I4,2E16.8)') i,j,0,0.0,theta(i,j,0)
!       DO k = 1,model_levels
!         WRITE(6,fmt='(3I4,3E16.8)') i,j,k,exner(i,j,k),theta(i,j,k),rho(i,j,k)
!       END DO
!     END DO
!   END DO
! WRITE(6,*)

!IF (me == 0) THEN
! WRITE (6,*) ' '
! WRITE (6,*) ' BUBBLE PERTURBATIONS '
!END IF

! Initialise bubble function to 0
bubblefn(:,:,:) = 0.0


IF ( bubble_option == gaussian_bubble ) THEN
!-------------------------------------------------

! Calculate a "Gaussian" bubble perturbation
!  EXP(-(dist/scale)**2)

!-------------------------------------------------

!  IF (me == 0) THEN
!   WRITE (6,*) ' Option 1: Gaussian bubble'
!   WRITE (6,FMT='(A,F8.2,A2)')                                       &
!          '   Initiating bubble with centre at height ',             &
!          idl_bubble_height,' m'
!   WRITE (6,FMT='(A,F5.2,A2)')                                       &
!          '   Max potential temperature perturbation = ',            &
!                              dt_bubble,' K'
!  END IF
! Calculate bubble 3D Gaussian function
  x0 = glob_xi1_p(1) + idl_bubble_xoffset                             &
      *(glob_xi1_p(global_row_length)-glob_xi1_p(1))
  y0 = glob_xi2_p(1) + idl_bubble_yoffset                             &
      *(glob_xi2_p(global_rows) - glob_xi2_p(1))
  z0 = idl_bubble_height 

!  IF (me == 0) THEN
!   WRITE(6,fmt='(A,3E)')  'bubble centre:',x0, y0, z0
!   WRITE(6,fmt='(A,6E16.8)')  'domain dimensions:',                   &
!           xi1_p(pdims%i_start),xi1_p(pdims%i_end),                   &
!           xi2_p(pdims%j_start),xi2_p(pdims%j_end),                   &
!           MINVAL(xi3_at_theta(:,:,0)),                               &
!           MAXVAL(xi3_at_theta(:,:,model_levels)) 
!   WRITE (6,*) ' '          
!  END IF  

  bubble_flat_radius = idl_bubble_depth

  DO k = 0, model_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        dist = SQRT( (xi1_p(i)-x0)**2 + (xi2_p(j)-y0)**2 +            &
                     (xi3_at_theta(pdims%i_start,j,k)-a0-z0)**2 )

        IF ( dist > bubble_flat_radius ) THEN
          tmp = ((dist-bubble_flat_radius)/idl_bubble_width)**2
          bubblefn(i,j,k) = EXP(-tmp)
        ELSE
          bubblefn(i,j,k) = 1.0
        END IF
! Add bubble temperature perturbation to the initial state
        theta(i,j,k) = theta(i,j,k) + dt_bubble*bubblefn(i,j,k)
      END DO
    END DO
  END DO

ELSE IF ( bubble_option == sin_bubble) THEN
!-------------------------------------------------

! Calculate a "Gaussian x-y, sin z" bubble perturbation

!-------------------------------------------------

!  IF (me == 0) THEN
!    WRITE (6,*) ' Option 2: Gaussian x-y, sin z bubble'
!    WRITE (6,FMT='(A,2F8.2,A2)')                                       &
!           '   Initiating bubble with width (x,y) ',             &
!           idl_bubble_width,idl_bubble_depth,' m'
!    WRITE (6,FMT='(A,F5.2,A2)')                                       &
!           '   Max potential temperature perturbation = ',            &
!                               dt_bubble,' K'
!  END IF
! Calculate bubble 3D Gaussian function
  x0 = glob_xi1_p(1) + idl_bubble_xoffset                             &
      *(glob_xi1_p(global_row_length)-glob_xi1_p(1))
  y0 = glob_xi2_p(1) + idl_bubble_yoffset                             &
      *(glob_xi2_p(global_rows) - glob_xi2_p(1))
  z0 = idl_bubble_height

!  IF (me == 0) THEN
!   WRITE(6,fmt='(A,3E)')  'bubble centre:',x0, y0, z0
!   WRITE(6,fmt='(A,6E16.8)')  'domain dimensions:',                   &
!           xi1_p(pdims%i_start),xi1_p(pdims%i_end),                   &
!           xi2_p(pdims%j_start),xi2_p(pdims%j_end),                   &
!           MINVAL(xi3_at_theta(:,:,0)),                               &
!           MAXVAL(xi3_at_theta(:,:,model_levels)) 
!   WRITE (6,*) ' '          
!  END IF

  bubble_flat_radius = idl_bubble_depth

  DO k = 0, model_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        dist = 1.0 + (xi1_p(i)-x0)**2/idl_bubble_width**2             & 
                   + (xi2_p(j)-y0)**2/idl_bubble_depth**2
        z0= (xi3_at_theta(i,j,k)-a0)                                  &
           /(xi3_at_theta(i,j,model_levels) - a0)

        bubblefn(i,j,k) = SIN(pi*z0)/dist

! Add bubble temperature perturbation to the initial state
        theta(i,j,k) = theta(i,j,k) + dt_bubble*bubblefn(i,j,k)
      END DO
    END DO
  END DO
!ELSE
!  WRITE(6,*) 'No bubble option specified'
END IF


! Correct rho
DO k = 1, model_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      t     = (intw_w2rho(k,1)*theta(i,j,k) +                         &
               intw_w2rho(k,2)*theta(i,j,k-1))*exner(i,j,k)
      rho(i,j,k) = p_zero*exner(i,j,k)**(1.0/kappa)/(r*t)
    END DO
  END DO
END DO

IF (lhook) CALL dr_hook('EG_IDL_INITIALISE_BUBBLE',zhook_out,zhook_handle)

END SUBROUTINE eg_idl_initialise_bubble
END MODULE eg_idl_initialise_bubble_mod
