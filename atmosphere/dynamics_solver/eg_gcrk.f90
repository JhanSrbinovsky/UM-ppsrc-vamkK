! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      SUBROUTINE eg_GCRk(x,b,k,tol,pre_type, l_rel_tol,                        &
             row_length, rows, n_rows, model_levels, halo_i, halo_j,           &
             offx, offy,                                                       &
                 mype, nproc, nproc_x,                                         &
                 nproc_y, global_row_length, global_rows,                      &
                 datastart, at_extremity, g_i_pe,                              &
                 gc_proc_row_group, gc_proc_col_group, model_domain,           &
                 sc_err_min, init_err, fin_err, no_its, ICODE )

      USE yomhook,  ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE atm_fields_bounds_mod
      USE eg_inner_prod_mod
      USE eg_calc_ax_mod

      IMPLICIT NONE

!
! Description: Preconditioned GCR(k) algorithm from Wesseling (2002)
!              to Solve Ax = b
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

      INTEGER, PARAMETER                              :: itmax = 999

! Array dimensions

      INTEGER,                             INTENT(IN)    :: offx, offy
      INTEGER,                             INTENT(IN)    :: row_length
      INTEGER,                             INTENT(IN)    :: rows
      INTEGER,                             INTENT(IN)    :: model_levels
      INTEGER,                             INTENT(IN)    :: halo_i,halo_j
      INTEGER,                             INTENT(IN)    :: n_rows
      INTEGER,                             INTENT(IN)    :: pre_type
      INTEGER,                             INTENT(IN)    :: k
      INTEGER,                             INTENT(INOUT) :: ICODE

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


      REAL,   INTENT(IN)    :: tol
      REAL,   INTENT(IN)    ::                                                 &
          b(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
      REAL,    INTENT(OUT)   ::                                                &
          x(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

! Temp work space arrays

      REAL                  ::                                                 &
          Ax(1-offx:row_length+offx,1-offy:rows+offy,model_levels),            &
           r(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

      REAL                  ::                                                 &
           v(1-offx:row_length+offx,1-offy:rows+offy,model_levels,k),      &
           s(1-offx:row_length+offx,1-offy:rows+offy,model_levels,k)

      REAL                  :: alf
      REAL                  :: err, sc_err, fin_tol
      INTEGER               :: it, m, n

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

      IF (lhook) CALL dr_hook('EG_GCRK',zhook_in,zhook_handle)

      sc_err = SQRT(eg_Inner_Prod(b,b))

! don't let sc_err get too small
      sc_err = MAX(sc_err, sc_err_min)

! Calculate initial residual

      CALL eg_Calc_Ax(Ax,x,                                                    &
             row_length, rows, n_rows, model_levels, halo_i, halo_j,           &
             offx, offy, model_domain)

      r = b - Ax

! Calculate initial error
      err = SQRT(eg_Inner_Prod(r,r))

      init_err = err/sc_err
      fin_tol  = tol
      IF( l_rel_tol ) fin_tol = init_err*tol

      DO it = 1, itmax
        DO m = 1, k

! DEPENDS ON: eg_Precon
          CALL eg_Precon(s(:,:,:,m),r,pre_type,                                &
             row_length, rows, n_rows, model_levels, halo_i, halo_j,           &
             offx, offy,model_domain,datastart)

          CALL eg_Calc_Ax(v(:,:,:,m),s(:,:,:,m),                               &
             row_length, rows, n_rows, model_levels, halo_i, halo_j,           &
             offx, offy, model_domain)

          DO n = 1, m-1

            alf = eg_Inner_Prod(v(:,:,:,m),v(:,:,:,n))

            v(:,:,:,m) = v(:,:,:,m) - alf*v(:,:,:,n)
            s(:,:,:,m) = s(:,:,:,m) - alf*s(:,:,:,n)
          END DO

          alf = SQRT( eg_Inner_Prod(v(:,:,:,m),v(:,:,:,m)))
          v(:,:,:,m) = v(:,:,:,m)/alf
          s(:,:,:,m) = s(:,:,:,m)/alf

          alf = eg_Inner_Prod( r,v(:,:,:,m))

          x = x + alf*s(:,:,:,m)
          r = r - alf*v(:,:,:,m)
        END DO

        err = SQRT(eg_Inner_Prod(r,r))/sc_err

        IF( err < fin_tol ) EXIT

      END DO

      IF( it >= itmax .and. err > fin_tol ) ICODE = 3

      fin_err = err
      no_its  = it

     IF (lhook) CALL dr_hook('EG_GCRK',zhook_out,zhook_handle)

     END SUBROUTINE eg_GCRk
