! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine Calc_PMSL
MODULE calc_pmsl_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE calc_pmsl(                                             &
                   theta, exner_theta_levels, p,                  &
                   row_length, rows, model_levels,                &
                   boundary_layer_levels,                         &
                   off_x, off_y, halo_i, halo_j,                  &
                   pmsl,p_star,                                   &

!-- New inputs to be added to call from diagnostics_end--------

                       me, n_proc, n_procx, n_procy,              &
                       neighbour, at_extremity,                   &
                       all_proc_group,model_domain,               &
                       delta_lambda, delta_phi,                   &
                       npmsl_height, l_pmsl_sor,                  &
                       sec_theta_latitude,                        &
                       global_row_length, global_rows)

!-------------------------------------------------------------


! Purpose:
!          Calculates Pressure at mean sea level.

! Method:
!          Modified version of that described in equations 3.9 and 3.11
!          in Unified Model Documentation
!          Paper Number S1. A. Dickinson

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids


! Code Description:
!   Language: FORTRAN 90 
!   This code is written to UMDP3 programming standards. v8.2

USE atm_fields_bounds_mod, ONLY: pdims_s, tdims_s, wdims_l, pdims_l

USE earth_constants_mod, ONLY: g, earth_radius
USE atmos_constants_mod, ONLY: r, lapse

USE level_heights_mod, ONLY:                                      &
                r_theta_levels, r_rho_levels,                     &
                eta_theta_levels

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE conversions_mod, ONLY: upperheight
USE calc_npmsl_mod, ONLY: calc_npmsl
USE calc_npmsl_redbl_mod, ONLY: calc_npmsl_redbl
IMPLICIT NONE

! Arguments with Intent IN. i.e.: Input variables.

INTEGER, INTENT(IN)  ::  row_length    ! number of points on a row
INTEGER, INTENT(IN)  ::  rows          ! number of rows of data
INTEGER, INTENT(IN)  ::  model_levels  ! number of levels of data
INTEGER, INTENT(IN)  ::  boundary_layer_levels 
                                       ! number of boundary layer levels
                                       
INTEGER, INTENT(IN)  ::  off_x 
INTEGER, INTENT(IN)  ::  off_y
INTEGER, INTENT(IN)  ::  halo_i
INTEGER, INTENT(IN)  ::  halo_j
INTEGER, INTENT(IN)  ::  me            !IN. Processor number
INTEGER, INTENT(IN)  ::  n_proc 
INTEGER, INTENT(IN)  ::  n_procx
INTEGER, INTENT(IN)  ::  n_procy
INTEGER, INTENT(IN)  ::  neighbour(4)
                         ! Array with the Ids of the four neighbours
                         ! in the horizontal plane
INTEGER, INTENT(IN)  ::  all_proc_group ! Group id for all processors
INTEGER, INTENT(IN)  ::  model_domain   
                         ! indicator as to model type, ie global, lam
INTEGER, INTENT(IN)  ::  global_row_length
                         ! NUMBER OF points on a global row
INTEGER, INTENT(IN)  ::  global_rows    !NUMBER OF global rows

LOGICAL, INTENT(IN)  ::  at_extremity(4)  
                         ! Indicates if this processor is at north,
                         ! south, east or west of the processor grid
                         
LOGICAL, INTENT(IN)  :: l_pmsl_sor
                         ! if false use more scalable smoothing
                        
REAL, INTENT(IN)  ::   delta_lambda
REAL, INTENT(IN)  ::   delta_phi 
REAL, INTENT(IN)  ::   npmsl_height 
                          ! Orographic height above which relaxation occurs
                          
                          

REAL, INTENT(IN) :: p(pdims_s%i_start:pdims_s%i_end,                   &
                     pdims_s%j_start:pdims_s%j_end,                    &
                     pdims_s%k_start:pdims_s%k_end)
REAL, INTENT(IN) :: theta(tdims_s%i_start:tdims_s%i_end,               &
                     tdims_s%j_start:tdims_s%j_end,                    &
                     tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN) :: exner_theta_levels                                 &
                     (tdims_s%i_start:tdims_s%i_end,                   &
                      tdims_s%j_start:tdims_s%j_end,                   &
                      tdims_s%k_start:tdims_s%k_end)
REAL, INTENT(IN) :: sec_theta_latitude (1-off_x:row_length+off_x,      &
                                                  1-off_y:rows+off_y)



REAL, INTENT(OUT) :: pmsl (row_length, rows)
REAL, INTENT(OUT) :: p_star (row_length, rows)


! Local variables

INTEGER  :: i, j  
INTEGER  :: k, levelupper

REAL     :: t_ref_level_1 
REAL     :: t_at_mean_sea_level (row_length, rows)
REAL     :: power 
REAL     :: phi_star (row_length, rows)
REAL     :: cos_p_latitude (row_length, rows) 
REAL     :: height_domain

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ----------------------------------------------------------------------
! Section 1. Calculate reference temperature at lowest pressure level
!            from value at level just above top of boundary layer and
!            assuming a fixed lapse rate with height.
!            Uses equation 3.9 from UMDP S1.
!            Calculate PMSL using equation 3.11 from UMDP S1 but with
!            level 1 values replacing * ones. the expression on the top
!            of the quotient in equation 3.11 is in fact simply the
!            temperature at mean sea level.
! ----------------------------------------------------------------------

IF (lhook) CALL dr_hook('CALC_PMSL',zhook_in,zhook_handle)

! find upper level (replacing level just above top of boundary layer!)
! height_domain is the top model height and is the same everywhere
height_domain = r_theta_levels(1,1,model_levels) - earth_radius

DO k = 1, model_levels
  levelupper = k
  IF ( height_domain * eta_theta_levels(k) > upperheight) EXIT
END DO

power = g / (r * lapse)


DO j = 1, rows
  DO i = 1, row_length

    t_ref_level_1 = theta(i,j,levelupper) *                       &
                 exner_theta_levels(i,j,levelupper)               &
                  + lapse *                                       &
                  (r_theta_levels(i,j,levelupper) -               &
                   r_rho_levels(i,j,1) )

    t_at_mean_sea_level(i,j) = t_ref_level_1                      &
                          + lapse *                               &
                           (r_rho_levels(i,j,1) -                 &
                            earth_radius )

    pmsl(i,j) = p(i,j,1) *                                        &
               (t_at_mean_sea_level(i,j) / t_ref_level_1 )        &
               ** power

  END DO
END DO

!   New method for pmsl from Unified Model


DO j = 1, rows
  DO i = 1, row_length
    phi_star(i,j)=g * (r_theta_levels(i,j,0)-earth_radius)
    cos_p_latitude(i,j)= 1./sec_theta_latitude(i,j)
  END DO
END DO

IF (l_pmsl_sor) THEN

  CALL  calc_npmsl(pmsl,p_star,                                   &
                 phi_star,theta,t_at_mean_sea_level,              &
                 cos_p_latitude,delta_lambda,delta_phi,           &
                 row_length,rows,                                 &
                 global_row_length,global_rows,                   &
                 me, n_proc, n_procx, n_procy,                    &
                 off_x, off_y,                                    &
                 neighbour, at_extremity,                         &
                 all_proc_group,model_domain,                     &
                 npmsl_height)

ELSE

  CALL  calc_npmsl_redbl(pmsl,p_star,                             &
                 phi_star,theta,t_at_mean_sea_level,              &
                 cos_p_latitude,delta_lambda,delta_phi,           &
                 row_length,rows,                                 &
                 global_row_length,global_rows,                   &
                 me, n_proc, n_procx, n_procy,                    &
                 off_x, off_y, halo_i, halo_j,                    &
                 neighbour, at_extremity,                         &
                 all_proc_group,model_domain,                     &
                 npmsl_height)

END IF

IF (lhook) CALL dr_hook('CALC_PMSL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_pmsl


END MODULE calc_pmsl_mod
