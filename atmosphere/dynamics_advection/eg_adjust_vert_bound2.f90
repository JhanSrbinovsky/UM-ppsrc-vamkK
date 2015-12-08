! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:

MODULE eg_adjust_vert_bound_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_adjust_vert_bound(                                      &
            row_length, rows, model_levels,g_i_pe,                    &
            pnt_type, dep_row_len,                                    &
            dep_rows, off_i, off_j, off_k, offz,                      &
            etadot, etadot_np1, depart_xi1, depart_xi2,               &
            depart_eta )

USE dynamics_input_mod,    ONLY : l_fast_vert_adjust,                 &
                                  l_sl_bc_correction
USE level_heights_mod,     ONLY : eta_theta_levels, eta_rho_levels
USE timestep_mod,          ONLY : timestep
USE um_parvars,            ONLY : halo_i, halo_j,                     &
                                  datastart,at_extremity

USE parkind1,              ONLY : jpim, jprb       !DrHook
USE yomhook,               ONLY : lhook, dr_hook   !DrHook

USE ereport_mod,           ONLY : ereport
USE Field_Types
USE horiz_grid_mod
USE atm_fields_bounds_mod
USE proc_info_mod,         ONLY : n_proc
USE PrintStatus_mod
USE proc_info_mod,         ONLY : global_row_length,                  &
                                  gc_proc_row_group,                  &
                                  model_domain,                       &
                                  mype=>me, nproc_x=>n_procx

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



INTEGER, INTENT(IN) :: row_length, rows, model_levels, pnt_type

! row_length, rows and indexing offset for u or v or w type point

INTEGER, INTENT(IN) :: dep_row_len, dep_rows, off_i, off_j,           &
                       off_k, offz

! MPP options
INTEGER, INTENT(IN) :: g_i_pe(1-halo_i:global_row_length+halo_i)

! Logical switches for advection

REAL, INTENT(IN) ::          etadot(wdims_s%i_start:wdims_s%i_end,    &
                                    wdims_s%j_start:wdims_s%j_end,    &
                                    wdims_s%k_start:wdims_s%k_end)
REAL, INTENT(IN) ::      etadot_np1(wdims_s%i_start:wdims_s%i_end,    &
                                    wdims_s%j_start:wdims_s%j_end,    &
                                    wdims_s%k_start:wdims_s%k_end)
! Departure point coordinates

REAL, INTENT(OUT) ::                                                  &
     depart_xi1(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
                1-off_k:model_levels)

REAL, INTENT(OUT) ::                                                  &
     depart_xi2(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
                1-off_k:model_levels)

REAL, INTENT(OUT) ::                                                  &
     depart_eta(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j,     &
                1-off_k:model_levels)

! Local variables

INTEGER :: i, j, k
REAL :: deta_1,                                                        &
        nlmax(model_levels), nlmin(model_levels), tau_lk,              &
        numax(model_levels), numin(model_levels), tau_uk

REAL ::                                                               &
  etadot_d(1-off_i:dep_row_len-off_i,1-off_j:dep_rows-off_j),         &
  etadot_1(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j),          &
  etadot_nm1(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j)

REAL    :: eta_low, eta_upp

INTEGER :: ierr, Nk

INTEGER check(model_levels)

! End of header

  
! 1.0 Start of subroutine code: perform the calculation.

  IF (lhook) CALL dr_hook('EG_ADJUST_VERT_BOUND2',zhook_in,zhook_handle)

  IF( pnt_type == fld_type_w ) THEN
     Nk = model_levels-1
  ELSE
     Nk = model_levels
  END IF

  IF( l_fast_vert_adjust .OR.                                         &
      (l_sl_bc_correction .AND. pnt_type /= fld_type_w) ) THEN
     SELECT CASE ( pnt_type )
       CASE ( fld_type_w )
          eta_low  = eta_theta_levels(0)
          eta_upp  = eta_theta_levels(model_levels)
        CASE ( fld_type_u, fld_type_v, fld_type_p )
          eta_low  = eta_rho_levels(1)
          eta_upp  = eta_rho_levels(model_levels)
        CASE DEFAULT
           CALL ereport("eg_adjust_vert_bound", 1, "Invalid grid point type" )
     END SELECT

     DO k=1, Nk
       DO j=1-off_j, dep_rows-off_j
          DO i=1-off_i, dep_row_len-off_i
             depart_eta(i,j,k) = MAX(depart_eta(i,j,k), eta_low)
             depart_eta(i,j,k) = MIN(depart_eta(i,j,k), eta_upp)
          END DO
       END DO
     END DO
 
  GO TO 9999 
  END IF

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

  check = 0

  DO k = 1, Nk
! bottom 
    SELECT CASE ( pnt_type )

      CASE ( fld_type_w )
         nlmax(k) = eta_theta_levels(k)
         nlmin(k) = eta_theta_levels(1)
      CASE ( fld_type_u, fld_type_v, fld_type_p )
         nlmax(k) = MAX ( eta_rho_levels(k)  ,  eta_theta_levels(1) )
         nlmin(k) = MIN ( eta_rho_levels(k)  ,  eta_theta_levels(1) )
      CASE DEFAULT
      CALL ereport("eg_adjust_vert_bound", 1,                         &
                   "Invalid grid point type" )
    END SELECT

! top
    SELECT CASE ( pnt_type )

      CASE ( fld_type_w )
          numax(k) = eta_theta_levels(model_levels-1)
          numin(k) = eta_theta_levels(k)
      CASE ( fld_type_u, fld_type_v, fld_type_p )
          numax(k) = MAX ( eta_rho_levels(k)  ,                          &
                        eta_theta_levels(model_levels-1) )
          numin(k) = MIN ( eta_rho_levels(k)  ,                          &
                        eta_theta_levels(model_levels-1) )
      CASE DEFAULT
          CALL ereport("eg_adjust_vert_bound", 1,                        &
                   "Invalid grid point type" )

    END SELECT

    IF ( ANY(depart_eta(1-off_i: dep_row_len-off_i,                      &
                        1-off_j: dep_rows-off_j,k) < nlmin(k)) )         &
         check(k) = 1

    IF ( ANY(depart_eta(1-off_i: dep_row_len-off_i,                      &
                        1-off_j: dep_rows-off_j,k) > numax(k)) )         &
         check(k) = 2

  END DO

  CALL gc_imax(Nk,n_proc,ierr,check)

  IF( mype == 0 .AND. PrintStatus > PrStatus_Normal) THEN
     WRITE(6,*)  "***********************************************"
     WRITE(6,*)  "* vertical adjustment performed at levels     *"
     WRITE(6,*)  "* =========================================== *"
     WRITE(6,*)  "* level:  where: 1=adjust bottom 2=adjust top *"
     DO k=1, model_levels-1
       IF(check(k)/=0) WRITE(6,fmt='(A,2I4,A)') ' *  ',k,&
          check(k), "                                   *"
     END DO
     WRITE(6,*)  "* if no output above then no adjustments  *"
     WRITE(6,*)  "*******************************************"
  END IF

  DO k=1, Nk

ne0:IF(check(k)/=0) THEN

chkk1:IF (check(k)==1) THEN

         deta_1  = eta_theta_levels(1) - eta_theta_levels(0)

! DEPENDS ON: eg_bi_linear_h
      CALL eg_bi_linear_h( etadot_1, depart_xi1(:,:,k:k),             &
            depart_xi2(:,:,k:k), fld_type_p, row_length, rows, 1,     &
            dep_row_len, dep_rows, 1, model_domain,                   &
            mype, nproc_x, halo_i, halo_j, datastart,                 &
            global_row_length, g_i_pe, at_extremity,                  &
            gc_proc_row_group, etadot_d )

        DO j=1-off_j, dep_rows-off_j
          DO i=1-off_i, dep_row_len-off_i

            IF ( depart_eta(i,j,k) < nlmin(k) ) THEN

              tau_lk = (nlmax(k) -  eta_theta_levels(1)) /            &
                (nlmax(k)-max(eta_theta_levels(0),depart_eta(i,j,k) ))

              depart_eta(i,j,k) = eta_theta_levels(0) +               &
                              (nlmin(k)-eta_theta_levels(0)) *        &
                              EXP(-timestep*(1. - tau_lk)/( deta_1 )  &
                              * ((1. - tau_lk)/2.*etadot_np1(i,j,1) + &
                              (1. + tau_lk)/2.*etadot_d(i,j)))

              depart_eta(i,j,k) = MIN( depart_eta(i,j,k) ,            &
                                   eta_theta_levels(1) )

            END IF
          END DO
        END DO
      END IF chkk1


chkk2:IF ( check(k)==2 ) THEN

    ! Adjust upper boundary.

        deta_1 = eta_theta_levels(model_levels) -                     &
                 eta_theta_levels(model_levels-1)

! DEPENDS ON: eg_bi_linear_h
      CALL eg_bi_linear_h( etadot_nm1, depart_xi1(:,:,k:k),           &
            depart_xi2(:,:,k:k), fld_type_p, row_length, rows, 1,     &
            dep_row_len, dep_rows, 1, model_domain,                   &
            mype, nproc_x, halo_i, halo_j, datastart,                 &
            global_row_length, g_i_pe, at_extremity,                  &
            gc_proc_row_group, etadot_d )

        DO j=1-off_j, dep_rows-off_j
          DO i=1-off_i, dep_row_len-off_i

            IF ( depart_eta(i,j,k) >  numax(k) ) THEN

              tau_uk = (eta_theta_levels(model_levels-1) - numin(k)) /&
                   (min(eta_theta_levels(model_levels),               &
                    depart_eta(i,j,k) )-numin(k))

              depart_eta(i,j,k) = eta_theta_levels(model_levels) -    &
                ( eta_theta_levels(model_levels) - numax(k) ) *       &
                exp( timestep * (1.-tau_uk)/( deta_1 ) *              &
                ((1.-tau_uk)/2.*etadot_np1(i,j,model_levels-1) +      &
                               (1.+tau_uk)/2.*etadot_d(i,j)))

              depart_eta(i,j,k) = MAX( depart_eta(i,j,k) ,            &
                                   eta_theta_levels(model_levels-1) )

            END IF
          END DO
        END DO

      END IF chkk2
    END IF ne0

  END DO          

9999 CONTINUE

IF (lhook) CALL dr_hook('EG_ADJUST_VERT_BOUND2',zhook_out,zhook_handle)

END SUBROUTINE eg_adjust_vert_bound
END MODULE
