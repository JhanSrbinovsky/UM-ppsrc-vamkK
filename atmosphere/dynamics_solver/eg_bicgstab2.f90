! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      SUBROUTINE eg_BiCGStab2(x,b,tol,pre_type, l_rel_tol,                     &
                 row_length, rows, n_rows, model_levels, halo_i, halo_j,       &
                 offx, offy,                                                   &
                 mype, nproc, nproc_x,                                         &
                 nproc_y, global_row_length, global_rows,                      &
                 datastart, at_extremity, g_i_pe,                              &
                 gc_proc_row_group, gc_proc_col_group, model_domain,           &
                 sc_err_min, init_err, fin_err, no_its, ICODE)


      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE proc_info_mod, ONLY : me
      USE integrity_mod
      USE atm_fields_bounds_mod
      USE eg_Inner_Prod_mod
      USE eg_calc_ax_mod
      USE eg_precon_mod

      IMPLICIT NONE

!
! Description: Preconditioned BiCGstab(2) algorithm
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

     INTEGER, PARAMETER                                 :: ItMax = 999

! Array dimensions

     INTEGER,                             INTENT(IN)    :: offx, offy
     INTEGER,                             INTENT(IN)    :: row_length
     INTEGER,                             INTENT(IN)    :: rows
     INTEGER,                             INTENT(IN)    :: model_levels
     INTEGER,                             INTENT(IN)    :: halo_i,halo_j
     INTEGER,                             INTENT(IN)    :: n_rows
     INTEGER,                             INTENT(IN)    :: pre_type
     INTEGER,                             INTENT(INOUT) :: ICODE

     REAL,                                INTENT(IN)    :: tol

! Data for parallel code and domain decomposition

     INTEGER, INTENT(IN) ::  mype, nproc, nproc_x, nproc_y,                   &
                              global_row_length, global_rows,                  &
                              datastart(3),                                    &
                              g_i_pe(1-halo_i:global_row_length+halo_i),       &
                              gc_proc_row_group, gc_proc_col_group,            &
                              model_domain

      LOGICAL, INTENT(IN) :: at_extremity(4), l_rel_tol

      REAL,    INTENT(IN)  :: sc_err_min
      REAL,    INTENT(OUT) :: init_err, fin_err
      INTEGER, INTENT(OUT) :: no_its

      REAL, INTENT(IN)    ::                                                   &
          b(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      REAL, INTENT(OUT)   ::                                                   &
          x(1-offx:row_length+offx,1-offy:rows+offy,model_levels)


! Local workspace arrays

      REAL                ::                                                   &
          r(1-offx:row_length+offx,1-offy:rows+offy,model_levels),                &
         cr(1-offx:row_length+offx,1-offy:rows+offy,model_levels),                &
         Ax(1-offx:row_length+offx,1-offy:rows+offy,model_levels),                &
          u(1-offx:row_length+offx,1-offy:rows+offy,model_levels),                &
          t(1-offx:row_length+offx,1-offy:rows+offy,model_levels),                &
          v(1-offx:row_length+offx,1-offy:rows+offy,model_levels),                &
          s(1-offx:row_length+offx,1-offy:rows+offy,model_levels),                &
          w(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      REAL    :: rho0, rho1, alf, omg1, omg2, bet
      REAL    :: gma,  ss, st, tt
      REAL    :: err, sc_err, fin_tol
      INTEGER :: it, i, j, k

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('EG_BICGSTAB2',zhook_in,zhook_handle)

      IF (integrity_test) CALL check_hash_m(b,SIZE(b),'rhs__')

! Define error with respect to the magnitude of the forcing

      sc_err = SQRT(eg_inner_prod(b,b))


! don't let sc_err get too small
      sc_err = MAX(sc_err, sc_err_min)

! Calculate initial residual

      CALL eg_calc_ax(Ax,x)

      DO k = 1, model_levels
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            r(i,j,k)  = b(i,j,k) - Ax(i,j,k)
          END DO
        END DO
      END DO

! Calculate initial error
      err = SQRT(eg_Inner_Prod(r,r))

      init_err = err/sc_err
      fin_tol  = tol
      IF( l_rel_tol ) fin_tol = init_err*tol

      CALL eg_precon(cr,r,pre_type)

      DO k = 1, model_levels
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            r(i,j,k) = cr(i,j,k)
            u(i,j,k) = 0.0
          END DO
        END DO
      END DO


      rho0 = 1.0
      alf  = 0.0
      omg2 = 1.0

      DO it = 1, ItMax
        rho0 = -omg2*rho0

! Odd step

        rho1 = eg_Inner_Prod(cr,r)

        bet  = alf*rho1/rho0
        rho0 = rho1

        DO k = 1, model_levels
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              u(i,j,k) = r(i,j,k) - bet*u(i,j,k)
            END DO
          END DO
        END DO

        CALL eg_Calc_Ax(Ax,u)

        CALL eg_precon(v,Ax,pre_type)

        gma  = eg_Inner_Prod(v,cr)


        alf  = rho0/gma
        DO k = 1, model_levels
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              r(i,j,k) = r(i,j,k) - alf*v(i,j,k)
            END DO
          END DO
        END DO

        CALL eg_calc_ax(Ax,r)

        CALL eg_precon(s,Ax,pre_type)

        DO k = 1, model_levels
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              x(i,j,k) = x(i,j,k) + alf*u(i,j,k)
            END DO
          END DO
        END DO

! Even step

        rho1 = eg_Inner_Prod(cr,s)

        bet  = alf*rho1/rho0
        rho0 = rho1

        DO k = 1, model_levels
          DO j = pdims%j_start, pdims%j_end
            DO i = pdims%i_start, pdims%i_end
              v(i,j,k) = s(i,j,k) - bet*v(i,j,k)
            END DO
          END DO
        END DO

        CALL eg_calc_ax(Ax,v)

        CALL eg_precon(w,Ax,pre_type)

         gma = eg_Inner_Prod(w,cr)
         alf = rho0/gma
         DO k = 1, model_levels
           DO j = pdims%j_start, pdims%j_end
             DO i = pdims%i_start, pdims%i_end
               u(i,j,k) = r(i,j,k) - bet*u(i,j,k)
               r(i,j,k) = r(i,j,k) - alf*v(i,j,k)
               s(i,j,k) = s(i,j,k) - alf*w(i,j,k)
             END DO
           END DO
         END DO

         CALL eg_calc_ax(Ax,s)

         CALL eg_precon(t,Ax,pre_type)

! GCR(2) step

         omg1 = eg_Inner_Prod(r,s)

         ss   = eg_Inner_Prod(s,s)

         st   = eg_Inner_Prod(s,t)

         tt   = eg_Inner_Prod(t,t)

         omg2 = eg_Inner_Prod(r,t)


         tt   = tt - st*st/ss
         omg2 = ( omg2 - st*omg1/ss )/tt
         omg1 = ( omg1 - st*omg2 )/ss

         DO k = 1, model_levels
           DO j = pdims%j_start, pdims%j_end
             DO i = pdims%i_start, pdims%i_end
               x(i,j,k) = x(i,j,k) + omg1*r(i,j,k) + omg2*s(i,j,k)             &
                                   + alf*u(i,j,k)
               r(i,j,k) = r(i,j,k) - omg1*s(i,j,k) - omg2*t(i,j,k)
               u(i,j,k) = u(i,j,k) - omg1*v(i,j,k) - omg2*w(i,j,k)
             END DO
           END DO
         END DO

! Check residual for convergence

         CALL eg_Calc_Ax(Ax,x)

         DO k = 1, model_levels
           DO j = pdims%j_start, pdims%j_end
             DO i = pdims%i_start, pdims%i_end
               Ax(i,j,k) = b(i,j,k) - Ax(i,j,k)
             END DO
           END DO
         END DO

         err = sqrt(eg_Inner_Prod(Ax,Ax))/sc_err


         IF( err < fin_tol ) EXIT

      END DO

      IF( err > fin_tol ) ICODE = 2

      fin_err = err
      no_its  = it

      IF (lhook) CALL dr_hook('EG_BICGSTAB2',zhook_out,zhook_handle)

      END SUBROUTINE eg_BiCGStab2
