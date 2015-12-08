! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
      MODULE eg_BiCGStab_mod
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE eg_BiCGStab(x,r,tol,pre_type, l_rel_tol,             &
                 row_length, rows, n_rows, model_levels,              &
                 sc_err_min, init_err, fin_err, no_its, l_inc, ICODE)

      USE um_parvars,          ONLY : offx, offy, halo_i, halo_j,     &
                                      datastart,at_extremity
      USE proc_info_mod,       ONLY : global_row_length,              &
                                      model_domain
      USE global_2d_sums_mod,  ONLY : global_2d_sums

      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE atm_fields_bounds_mod
      USE eg_Inner_Prod_mod
      USE gmres1_coef_mod
      USE eg_calc_ax_mod
      USE eg_precon_mod

      IMPLICIT NONE

!
! Description: Postconditioned BiCGstab algorithm from Wesseling (2002)
!              to Solve Ax = b
!
! Method:
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

      INTEGER, PARAMETER   :: ItMax = 999
      REAL,    PARAMETER   :: small = 1.0e-16

! Array dimensions 

      INTEGER,                             INTENT(IN)    :: row_length
      INTEGER,                             INTENT(IN)    :: rows
      INTEGER,                             INTENT(IN)    :: model_levels
      INTEGER,                             INTENT(IN)    :: n_rows
      INTEGER,                             INTENT(IN)    :: pre_type
      INTEGER,                             INTENT(INOUT) :: ICODE

! Data for parallel code and domain decomposition


      LOGICAL, INTENT(IN)  :: l_rel_tol, l_inc

      REAL,    INTENT(IN)  :: sc_err_min
      REAL,    INTENT(OUT) :: init_err, fin_err
      INTEGER, INTENT(OUT) :: no_its


      REAL, INTENT(INOUT) ::                                                   &
          r(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      REAL, INTENT(OUT)   ::                                                   &
          x(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Residual error and "Stabilizer" used in the iterations

      REAL                ::                                                   &
          cr(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Temp work space arrays

      REAL                ::                                                   &
          Ax(1-offx:row_length+offx,1-offy:rows+offy,model_levels),            &
          p(1-offx:row_length+offx,1-offy:rows+offy,model_levels),             &
          t(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      REAL                ::                                                   &
          v(1-offx:row_length+offx,1-offy:rows+offy,model_levels),             &
          s(1-offx:row_length+offx,1-offy:rows+offy,model_levels),             &
          cs(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Absolute tolerance 
      REAL                :: tol

! Temp variables
      REAL                :: rho, alf
      REAL                :: omg, bet
      REAL                :: nrm

! Variables used to define the normalized error
      REAL                :: err, sc_err, fin_tol
! Iteration count
      INTEGER             :: it, i, j, k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('EG_BICGSTAB',zhook_in,zhook_handle)

! Define error with respect to the magnitude of the forcing

      sc_err = SQRT(eg_Inner_Prod(r,r))

! don't let sc_err get too small
      sc_err = MAX(sc_err, sc_err_min)

! Calculate initial residual

      IF( .NOT. l_inc ) THEN

         CALL eg_Calc_Ax(Ax,x)

!$OMP PARALLEL DO PRIVATE(i,j) SHARED(model_levels,rows,row_length,r,Ax,cr)
         DO k = 1, model_levels
            DO j = 1, rows
               DO i = 1, row_length
                  r(i,j,k)  = r(i,j,k) - Ax(i,j,k)
                  cr(i,j,k) = r(i,j,k)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO

      ELSE
!$OMP PARALLEL DO PRIVATE(i,j) SHARED(model_levels,rows,row_length,r,Ax,cr)
         DO k = 1, model_levels
            DO j = 1, rows
               DO i = 1, row_length
                  cr(i,j,k) = r(i,j,k)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO
      END IF

! Calculate initial error
      rho = eg_Inner_Prod(r,r)
      err = SQRT(rho)
       
      init_err = err/sc_err
      fin_tol  = tol
      IF( l_rel_tol ) fin_tol = init_err*tol

      alf = 1.0
      omg = 1.0
      nrm = 1.0            

      DO it = 1, itmax

         IF( it > 1 ) rho = eg_Inner_Prod(r,cr)

         bet = ( rho/nrm )*( alf/omg )

         IF( it == 1 ) THEN

            CALL eg_precon(p,r,pre_type)
         ELSE

!$OMP PARALLEL DO PRIVATE(i,j) SHARED(model_levels,rows,row_length,t,r,v,bet,omg)
            DO k = 1, model_levels
               DO j = 1, rows
                  DO i = 1, row_length
                     t(i,j,k)   = r(i,j,k) - bet*omg*v(i,j,k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO

            CALL eg_precon(s,t,pre_type)

!$OMP PARALLEL DO PRIVATE(i,j) SHARED(model_levels,rows,row_length,p,s,bet)
            DO k = 1, model_levels
               DO j = 1, rows
                  DO i = 1, row_length
                     p(i,j,k) = s(i,j,k) + bet*p(i,j,k)
                  END DO
               END DO
            END DO
!$OMP END PARALLEL DO
         END IF

         CALL eg_Calc_Ax(v,p)

         nrm = eg_Inner_Prod(cr,v)

         alf = rho/nrm

!$OMP PARALLEL DO PRIVATE(i,j) SHARED(model_levels,rows,row_length,s,r,alf,v)
         DO k = 1, model_levels
            DO j = 1, rows
               DO i = 1, row_length
                  s(i,j,k)   = r(i,j,k) - alf*v(i,j,k)
               END DO
            END DO
         END DO
!$OMP END PARALLEL DO

         CALL eg_precon(cs,s,pre_type)

         CALL eg_calc_Ax(t,cs)

         omg = gmres1_coef(s,t,err)

!$OMP PARALLEL DO PRIVATE(i,j) SHARED(model_levels,rows,row_length,x,alf,       &
!$OMP&             p,omg,cs,r,s,t)
         DO k = 1, model_levels
           x(:,:,k)   = x(:,:,k) + alf*p(:,:,k) + omg*cs(:,:,k)
           DO j = 1, rows
               DO i = 1, row_length
                  r(i,j,k)   = s(i,j,k) - omg*t(i,j,k)
               END DO
           END DO
         END DO
!$OMP END PARALLEL DO

         nrm = rho

         IF( abs(omg) < small ) THEN
            ICODE = 11
            IF (lhook) CALL dr_hook('EG_BICGSTAB',zhook_out,zhook_handle)
            RETURN
         END IF

! Check residual for convergence

         err = SQRT(err)/sc_err
         IF( err < fin_tol ) EXIT

      END DO

      IF( err > fin_tol ) ICODE = 1

      fin_err = err
      no_its  = it

      IF (lhook) CALL dr_hook('EG_BICGSTAB',zhook_out,zhook_handle)

      END SUBROUTINE eg_BiCGStab
END MODULE eg_BiCGStab_mod
