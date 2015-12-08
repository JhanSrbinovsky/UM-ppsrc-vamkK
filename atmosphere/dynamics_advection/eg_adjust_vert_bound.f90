! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:

MODULE eg_adjust_vert_bound_mod_old
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_adjust_vert_bound(                                      &
            row_length, rows, n_rows, model_levels, halo_i,halo_j,    &
            offx, offy, nproc_x, global_row_length, datastart,        &
            at_extremity, mype, g_i_pe, gc_proc_row_group,            &
            l_sl_halo_reprod, model_domain, pnt_type, dep_row_len,    &
            dep_rows, off_i, off_j, off_k, offz,                      &
            interp_vertical_search_tol, check_bottom_levels,          &
            l_regular, timestep,eta_theta_levels, eta_rho_levels,     &
            etadot, etadot_np1, depart_xi1, depart_xi2,               &
            depart_eta )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod, ONLY : ereport
USE Field_Types
USE horiz_grid_mod
USE atm_fields_bounds_mod
USE proc_info_mod,            ONLY: n_proc

IMPLICIT NONE
!
! Description:
!  
!   Adjust departure points crossing the top and bottom full level using
!   the method described by eqs (10.51) and (10.52) of ENDGame formulation.
!
! Method: ENDGame formulation version 1.03, section 10.
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



INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels,        &
                       model_domain, pnt_type

! row_length, rows and indexing offset for u or v or w type point

INTEGER, INTENT(IN) :: dep_row_len, dep_rows, off_i, off_j,           &
                       off_k, offz

! MPP options
INTEGER, INTENT(IN) ::                                                &
  mype,                                                               &
                     ! Processor number
  nproc_x,                                                            &
                     ! Number of processors in longitude
  halo_i,                                                             &
                     ! Size of halo in i.
  halo_j,                                                             &
                     ! Size of halo in j.
  offx,                                                               &
                     ! Size of small halo in i
  offy,                                                               &
                     ! Size of small halo in j.
  datastart(3),                                                       &
                     ! First gridpoints held by this processor.
  gc_proc_row_group,                                                  &
                     ! Group id for processors on the same row
  global_row_length,                                                  &
                     ! global number of points on a row
  g_i_pe(1-halo_i:global_row_length+halo_i)


LOGICAL, INTENT(IN) ::                                                &
  at_extremity(4)  ! Indicates if this processor is at north,
                   ! south, east or west of the processor grid

! Loop index bounds for arrays defined on p, u, v points respectively


! Integer parameters for advection

INTEGER, INTENT(IN) ::                                                &
  interp_vertical_search_tol,                                         &
                      ! used in interpolation code.
  check_bottom_levels ! used in interpolation code, and is
                      ! the number of levels to check to see
                      ! if the departure point lies inside the
                      ! orography.

! Logical switches for advection

LOGICAL, INTENT(IN) ::                                                &
  l_regular,                                                          &
                   ! True, then regular grid to be interpolated
  l_sl_halo_reprod
                   ! True, SL code bit repoducible with
                   !       any sensible halo size

REAL, INTENT(IN) :: timestep

REAL, INTENT(IN) ::                                                   &
  eta_theta_levels(0:model_levels),                                   &
  eta_rho_levels(1-offz:model_levels+offz)

REAL, INTENT(IN) ::                                                   &
  etadot(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
  etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,                 &
             0:model_levels)

! Departure point coordinates

REAL, INTENT(OUT) ::                 &
    depart_xi1(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
    1-off_k:model_levels)
REAL, INTENT(OUT) ::                 &
    depart_xi2(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
    1-off_k:model_levels)
REAL, INTENT(OUT) ::                 &
    depart_eta(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
    1-off_k:model_levels)
! Local variables

INTEGER :: i, j, k
REAL :: deta_1, eta_min, eta_max,                                     &
        nlmax, nlmin, tau_lk, numax, numin, tau_uk

REAL ::                                                               &
  etadot_d(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j),         &
  etadot_1(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j),          &
  etadot_nm1(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j)

INTEGER :: k_max(2),ierr

! End of header

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_ADJUST_VERT_BOUND',zhook_in,zhook_handle)

k_max(1) = 0
k_max(2) = 0

! choose method of boundary adjustment
IF ( interp_vertical_search_tol == 0 ) THEN

  SELECT CASE ( pnt_type )

    CASE ( fld_type_u, fld_type_v, fld_type_p )
      eta_min = eta_rho_levels(1)
      eta_max = eta_rho_levels(model_levels)
    CASE ( fld_type_w )
      eta_min = eta_theta_levels(0)
      eta_max = eta_theta_levels(model_levels)
    CASE DEFAULT

     Call ereport("eg_adjust_vert_bound", 1,                          &
                  "Invalid grid point type" )
  END SELECT

! 0-level in w-points ignored as it is clamped (see end of sub)

!$omp parallel do private(i,j)
  DO k=1, model_levels
    DO j=1-off_j, dep_rows-off_j
      DO i=1-off_i, dep_row_len-off_i
        IF ( depart_eta(i,j,k) < eta_min ) THEN
          depart_eta(i,j,k) = eta_min
        END IF
        IF ( depart_eta(i,j,k) > eta_max ) THEN
          depart_eta(i,j,k) = eta_max
        END IF
      END DO
    END DO
  END DO
!$omp end parallel do

ELSE

! Copy etadot at lev 1 and ml-1 at single level extended arrays for use
! by eg_bi_linear_h()

  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      etadot_1(i,j)   = etadot(i,j,1) 
      etadot_nm1(i,j) = etadot(i,j,model_levels-1) 
    END DO
  END DO

!DEPENDS ON: swap_bounds
  CALL swap_bounds( etadot_1,row_length,rows,1,halo_i,                &
                                         halo_j,fld_type_p,.FALSE. )
!DEPENDS ON: swap_bounds
  CALL swap_bounds( etadot_nm1,row_length,rows,1,halo_i,              &
                                         halo_j,fld_type_p,.FALSE. )

! Adjust lower boundary.

  deta_1  = eta_theta_levels(1) - eta_theta_levels(0)

  DO k=1, check_bottom_levels

! DEPENDS ON: eg_bi_linear_h
    CALL eg_bi_linear_h( etadot_1, depart_xi1(:,:,k:k),               &
            depart_xi2(:,:,k:k), fld_type_p, row_length, rows, 1,     &
            dep_row_len, dep_rows, 1, delta_xi1, delta_xi2,           &
            base_xi1, base_xi2, xi1_p, xi2_p, model_domain,           &
            mype, nproc_x, halo_i, halo_j, datastart,                 &
            global_row_length, g_i_pe, at_extremity, 2,               &
            gc_proc_row_group, l_sl_halo_reprod, l_regular,           &
            etadot_d )

    SELECT CASE ( pnt_type )

      CASE ( fld_type_w )
         nlmax = eta_theta_levels(k)
         nlmin = eta_theta_levels(1)
      CASE ( fld_type_u, fld_type_v, fld_type_p )
         nlmax = max ( eta_rho_levels(k)  ,  eta_theta_levels(1) )
         nlmin = min ( eta_rho_levels(k)  ,  eta_theta_levels(1) )
      CASE DEFAULT
      Call ereport("eg_adjust_vert_bound", 1,                         &
                   "Invalid grid point type" )

    END SELECT

    DO j=1-off_j, dep_rows-off_j
      DO i=1-off_i, dep_row_len-off_i

        IF ( depart_eta(i,j,k) < nlmin ) THEN

          k_max(1) = max(k_max(1),k)

          tau_lk = (nlmax -  eta_theta_levels(1)) /                   &
                   (nlmax-max(eta_theta_levels(0),depart_eta(i,j,k) ))

          depart_eta(i,j,k) = eta_theta_levels(0) +                   &
                              (nlmin-eta_theta_levels(0)) *           &
                              EXP(-timestep*(1. - tau_lk)/( deta_1 )  &
                              * ((1. - tau_lk)/2.*etadot_np1(i,j,1) + &
                              (1. + tau_lk)/2.*etadot_d(i,j)))

          depart_eta(i,j,k) = MIN( depart_eta(i,j,k) ,                &
                                   eta_theta_levels(1) )

        END IF

      END DO
    END DO

  END DO

! Adjust upper boundary.

  deta_1 = eta_theta_levels(model_levels) -                           &
           eta_theta_levels(model_levels-1)

  DO k= model_levels-interp_vertical_search_tol, model_levels-1

! DEPENDS ON: eg_bi_linear_h
    CALL eg_bi_linear_h(etadot_nm1, depart_xi1(:,:,k:k),              &
            depart_xi2(:,:,k:k), fld_type_p, row_length, rows, 1,     &
            dep_row_len, dep_rows, 1, delta_xi1, delta_xi2,           &
            base_xi1, base_xi2, xi1_p, xi2_p, model_domain,           &
            mype, nproc_x, halo_i, halo_j, datastart,                 &
            global_row_length, g_i_pe, at_extremity, 2,               &
            gc_proc_row_group, l_sl_halo_reprod, l_regular,           &
            etadot_d )

    SELECT CASE ( pnt_type )

      CASE ( fld_type_w )
          numax = eta_theta_levels(model_levels-1)
          numin = eta_theta_levels(k)
      CASE ( fld_type_u, fld_type_v, fld_type_p )
          numax = max ( eta_rho_levels(k)  ,                          &
                        eta_theta_levels(model_levels-1) )
          numin = min ( eta_rho_levels(k)  ,                          &
                        eta_theta_levels(model_levels-1) )
      CASE DEFAULT
      Call ereport("eg_adjust_vert_bound", 1,                         &
                   "Invalid grid point type" )

    END SELECT

    DO j=1-off_j, dep_rows-off_j
      DO i=1-off_i, dep_row_len-off_i

        IF ( depart_eta(i,j,k) >  numax ) THEN

          k_max(2) = max(k_max(2),ABS(k-model_levels+1))

          tau_uk = (eta_theta_levels(model_levels-1) - numin) /       &
                   (min(eta_theta_levels(model_levels),               &
                    depart_eta(i,j,k) )-numin)

          depart_eta(i,j,k) = eta_theta_levels(model_levels) -        &
              ( eta_theta_levels(model_levels) - numax ) *            &
                exp( timestep * (1.-tau_uk)/( deta_1 ) *              &
                ((1.-tau_uk)/2.*etadot_np1(i,j,model_levels-1) +      &
                               (1.+tau_uk)/2.*etadot_d(i,j)))

          depart_eta(i,j,k) = MAX( depart_eta(i,j,k) ,                &
                                   eta_theta_levels(model_levels-1) )

        END IF

      END DO
    END DO
 
  END DO          

  CALL gc_imax(2,n_proc,ierr,k_max)

  IF((k_max(1).ge.check_bottom_levels-1)) THEN

      Call ereport("eg_adjust_vert_bound", -k_max(1),                 &
                   "search tolerances at limit (bottom)!" )
  END IF

  IF((k_max(2).ge.interp_vertical_search_tol-1)) THEN

      Call ereport("eg_adjust_vert_bound", -k_max(2),                 &
                   "search tolerances at limit (top)!" )
  END IF

END IF ! IF interp_vertical_search_tol == 0          



IF (lhook) CALL dr_hook('EG_ADJUST_VERT_BOUND',zhook_out,zhook_handle)

RETURN
END SUBROUTINE eg_adjust_vert_bound
END MODULE
