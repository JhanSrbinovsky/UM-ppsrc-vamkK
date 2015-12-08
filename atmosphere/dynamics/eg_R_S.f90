! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_r_s_mod
  IMPLICIT NONE
CONTAINS
SUBROUTINE eg_r_s(                                                    &
           theta_star,q_star,qcl_star,qcf_star,qcf2_star,             &
           qrain_star, qgraup_star , theta_star_n ,m_star_n,          &
           mcl_star_n, mcf_star_n, mcf2_star_n, mrain_star_n,         &
           mgraup_star_n, l_physics, l_skeb2)

USE eg_implicit_horz_drag_mod, ONLY : eg_implicit_horz_drag
USE eg_horz_drag_mod,          ONLY : l_impl_horz_drag
USE parkind1,                  ONLY : jpim, jprb       !DrHook
USE yomhook,                   ONLY : lhook, dr_hook   !DrHook
USE atmos_constants_mod,       ONLY : epsln => repsilon
USE proc_info_mod,       ONLY : me,n_proc, model_domain
USE eg_q_to_mix_mod
USE atm_fields_bounds_mod
USE integrity_mod
USE field_types
USE um_input_control_mod, ONLY: l_mr_physics2
USE mphys_inputs_mod,     ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain

USE PrintStatus_mod

USE domain_params
USE eg_swap_bounds_mod
USE um_parvars,           ONLY: offx,offy
USE fields_rhs_mod,       ONLY: r_u_p2, r_v_p2, r_w_p2 , r_u_p2_n,    &
                                r_v_p2_n , r_w_p2_n, s_u,s_v,s_w,     &
                                s_thetav,s_m_v,s_m_cl,s_m_cf,         &
                                s_m_cf2,s_m_r,s_m_gr,                 &
                                r_u_skeb,r_v_skeb

IMPLICIT NONE
!
! Description:
!  
!    Computes the fast increment S from the increment R and the
!     * state as returned from atmos_physics2

!    theta* is converted into theta_vd. The increment is
!    then computed by subtracting the value before the
!    call to atmos_physics2 saved in R_X_SP
!
!    At the end a simple implementation of boundary layer drag
!    is optionally added.
!
! Method:
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments


LOGICAL l_physics

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


!increments from atmos_physics2:
REAL ::      theta_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::          m_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::        mcl_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::        mcf_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::      mcf2_star (tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::      mrain_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::     mgraup_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::          q_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::        qcl_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::        qcf_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::       qcf2_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::      qrain_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::     qgraup_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)
                               
! previously stored increments
! (before atmos physics 2 aka F1SP predictor:
REAL ::        m_star_n(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::      mcl_star_n(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::     mcf_star_n (tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)                 

REAL ::     mcf2_star_n(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::    mrain_star_n(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::   mgraup_star_n(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::    theta_star_n(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)


REAL    :: inv_epsilon
INTEGER :: i,j,k
INTEGER :: ierr

! local diagnostics
REAL :: max_s_theta,  min_s_theta,   av_s_theta,   has_nan_s_theta
REAL :: max_s_m_v,    min_s_m_v,     av_s_m_v,     has_nan_s_m_v
REAL :: max_s_m_cl,   min_s_m_cl,    av_s_m_cl,    has_nan_s_m_cl
REAL :: max_s_m_cf,   min_s_m_cf,    av_s_m_cf,    has_nan_s_m_cf
REAL :: max_s_m_cf2,  min_s_m_cf2,   av_s_m_cf2,   has_nan_s_m_cf2
REAL :: max_s_m_rain, min_s_m_rain,  av_s_m_rain,  has_nan_s_m_rain
REAL :: max_s_m_graup,min_s_m_graup, av_s_m_graup, has_nan_s_m_graup
REAL :: min_s_u,      max_s_u,       av_s_u,       has_nan_s_u
REAL :: min_s_v,      max_s_v,       av_s_v,       has_nan_s_v
REAL :: min_s_w,      max_s_w,       av_s_w,       has_nan_s_w

REAL :: fieldsize
LOGICAL :: l_skeb2

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_R_S',zhook_in,zhook_handle)

inv_epsilon = 1./epsln

IF (.NOT.l_physics) THEN
  IF (lhook) CALL dr_hook('EG_R_S',zhook_out,zhook_handle)
  RETURN
END IF

IF (l_mr_physics2) THEN
  DO k=tdims_s%k_start,tdims_s%k_end
    m_star  (:,:,k) = q_star  (:,:,k)
    mcl_star(:,:,k) = qcl_star(:,:,k)
    mcf_star(:,:,k) = qcf_star(:,:,k)
  END DO

  IF (l_mcr_qcf2  ) THEN
    DO k=tdims_s%k_start,tdims_s%k_end
      mcf2_star  (:,:,k) = qcf2_star  (:,:,k)
    END DO
  END IF

  IF (l_mcr_qrain ) THEN
    DO k=tdims_s%k_start,tdims_s%k_end
      mrain_star (:,:,k) = qrain_star (:,:,k)
    END DO
  END IF

  IF (l_mcr_qgraup) THEN
    DO k=tdims_s%k_start,tdims_s%k_end
      mgraup_star(:,:,k) = qgraup_star(:,:,k)
    END DO
  END IF
ELSE
  CALL eg_q_to_mix                                                    &
                  (tdims_s,tdims_s,                                   &
                   q_star, qcl_star, qcf_star,                        &
                   qcf2_star, qrain_star, qgraup_star,                &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,             &
                   m_star, mcl_star, mcf_star,                        &
                   mcf2_star, mrain_star, mgraup_star,.FALSE.)
END IF

!$OMP PARALLEL  PRIVATE(i,j,k) SHARED(udims, s_u,r_u_p2,              &
!$OMP&            r_u_p2_n,vdims,s_v,r_v_p2,r_v_p2_n,pdims,s_w,       &
!$OMP&            r_w_p2,r_w_p2_n,r_u_skeb,r_v_skeb,tdims_s,          &
!$OMP&            theta_star,m_star,inv_epsilon,s_m_v,s_m_cl,mcl_star,&
!$OMP&            mcl_star_n,s_m_cf,mcf_star,mcf_star_n,s_thetav,     &
!$OMP&            tdims,mrain_star,mrain_star_n,mgraup_star,          &
!$OMP&            l_mcr_qcf2,l_mcr_qrain, l_mcr_qgraup,s_m_r,s_m_gr,  &
!$OMP&            s_m_cf2,l_skeb2,theta_star_n,m_star_n,              &
!$OMP&            mgraup_star_n,mcf2_star,mcf2_star_n) DEFAULT(NONE)

IF (l_skeb2) THEN
!$OMP DO SCHEDULE (STATIC)
  DO k=udims%k_start,udims%k_end
    DO j=udims%j_start,udims%j_end
      DO i=udims%i_start,udims%i_end
        s_u(i,j,k) = r_u_p2(i,j,k) - r_u_p2_n(i,j,k) + r_u_skeb(i,j,k)
      END DO
    END DO

    DO j=vdims%j_start,vdims%j_end
      DO i=vdims%i_start,vdims%i_end
        s_v(i,j,k) = r_v_p2(i,j,k) - r_v_p2_n(i,j,k) + r_v_skeb(i,j,k)
      END DO
    END DO

!   zeroth level for s_w is done below (outside of OMP section)
    DO j=pdims%j_start,pdims%j_end
      DO i=pdims%i_start,pdims%i_end
        s_w(i,j,k) = r_w_p2(i,j,k) - r_w_p2_n(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
ELSE
!$OMP DO SCHEDULE (STATIC)
  DO k=udims%k_start,udims%k_end
    DO j=udims%j_start,udims%j_end
      DO i=udims%i_start,udims%i_end
        s_u(i,j,k) = r_u_p2(i,j,k) - r_u_p2_n(i,j,k)
      END DO
    END DO

    DO j=vdims%j_start,vdims%j_end
      DO i=vdims%i_start,vdims%i_end
        s_v(i,j,k) = r_v_p2(i,j,k) - r_v_p2_n(i,j,k)
      END DO
    END DO

!   zeroth level for s_w is done below (outside of OMP section)
    DO j=pdims%j_start,pdims%j_end
      DO i=pdims%i_start,pdims%i_end
        s_w(i,j,k) = r_w_p2(i,j,k) - r_w_p2_n(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!$OMP DO SCHEDULE (STATIC)
DO k=tdims_s%k_start,tdims_s%k_end

      s_thetav(tdims%i_start:tdims%i_end,                             &
               tdims%j_start:tdims%j_end,k) =                         &
    theta_star(tdims%i_start:tdims%i_end,                             &
               tdims%j_start:tdims%j_end,k)                           &
      * ( 1.0 + m_star(tdims%i_start:tdims%i_end,                     &
                       tdims%j_start:tdims%j_end,k)    * inv_epsilon)-&
          theta_star_n(tdims%i_start:tdims%i_end,                     &
                       tdims%j_start:tdims%j_end,k)                   &
    * ( 1.0 + m_star_n(tdims%i_start:tdims%i_end,                     &
                       tdims%j_start:tdims%j_end,k) * inv_epsilon)


      s_m_v (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,k)=  &
     m_star (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,k)-  &
   m_star_n (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,k)

      s_m_cl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,k)=  &
    mcl_star(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,k)-  &
  mcl_star_n(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,k)

      s_m_cf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,k)=  &
    mcf_star(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,k)-  &
 mcf_star_n (tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,k)
ENDDO
!$OMP END DO NOWAIT

IF (l_mcr_qcf2 ) THEN
!$OMP DO SCHEDULE(STATIC)
   DO k=tdims_s%k_start,tdims_s%k_end
     DO j=pdims%j_start,pdims%j_end
       DO i=pdims%i_start,pdims%i_end
          s_m_cf2(i,j,k) = mcf2_star(i,j,k) - mcf2_star_n(i,j,k)
       END DO
     END DO
   END DO
!$OMP END DO NOWAIT
ENDIF

IF (l_mcr_qrain ) THEN
!$OMP DO SCHEDULE(STATIC)
   DO k=tdims_s%k_start,tdims_s%k_end
     DO j=pdims%j_start,pdims%j_end
       DO i=pdims%i_start,pdims%i_end
         s_m_r(i,j,k) = mrain_star(i,j,k) - mrain_star_n(i,j,k)
       END DO
     END DO
   END DO
!$OMP END DO NOWAIT
ENDIF

IF (l_mcr_qgraup) THEN
!$OMP DO SCHEDULE(STATIC)
   DO k=tdims_s%k_start,tdims_s%k_end
      DO j=pdims%j_start,pdims%j_end
         DO i=pdims%i_start,pdims%i_end
            s_m_gr(i,j,k) = mgraup_star(i,j,k) - mgraup_star_n(i,j,k)
         END DO
      END DO
   END DO
!$OMP END DO NOWAIT
ENDIF

!$OMP END PARALLEL

     s_w(pdims%i_start:pdims%i_end,                                   &
         pdims%j_start:pdims%j_end,0)                                 &
       =                                                              &
         r_w_p2(pdims%i_start:pdims%i_end,                            &
                pdims%j_start:pdims%j_end,0) -                        &
       r_w_p2_n(pdims%i_start:pdims%i_end,                            &
                pdims%j_start:pdims%j_end,0)

!
!
! apply simple boundary layer drag if switched on
! this should use the fast physics predictors as u and v!
IF (l_impl_horz_drag) CALL eg_implicit_horz_drag(s_u,s_v,             &
                                                 r_u_p2_n,r_v_p2_n)



IF (printstatus >=  prstatus_diag) THEN

  fieldsize       = REAL((udims%i_end-udims%i_start+1) *              &
                         (udims%j_end-udims%j_start+1) *              &
                         (udims%k_end-udims%k_start+1) * n_proc)
  max_s_u         = MAXVAL(s_u(udims%i_start:udims%i_end,             &
                            udims%j_start:udims%j_end,                &
                            udims%k_start:udims%k_end))
  min_s_u         = MINVAL(s_u(udims%i_start:udims%i_end,             &
                            udims%j_start:udims%j_end,                &
                            udims%k_start:udims%k_end))
   av_s_u         = SUM(s_u(udims%i_start:udims%i_end,                &
                            udims%j_start:udims%j_end,                &
                            udims%k_start:udims%k_end))/fieldsize
  has_nan_s_u = 0.
!  IF(ANY(s_u.ne.s_u)) has_nan_s_u = 1.

  fieldsize       = REAL((vdims%i_end-vdims%i_start+1) *              &
                         (vdims%j_end-vdims%j_start+1) *              &
                         (vdims%k_end-vdims%k_start+1) * n_proc)
  max_s_v         = MAXVAL(s_v(vdims%i_start:vdims%i_end,             &
                            vdims%j_start:vdims%j_end,                &
                            vdims%k_start:vdims%k_end))
  min_s_v         = MINVAL(s_v(vdims%i_start:vdims%i_end,             &
                            vdims%j_start:vdims%j_end,                &
                            vdims%k_start:vdims%k_end))
   av_s_v         = SUM(s_v(vdims%i_start:vdims%i_end,                &
                            vdims%j_start:vdims%j_end,                &
                            vdims%k_start:vdims%k_end))/fieldsize
  has_nan_s_v = 0.
!  IF(ANY(s_v.ne.s_v)) has_nan_s_v = 1.

  fieldsize       = REAL((wdims%i_end-wdims%i_start+1) *              &
                         (wdims%j_end-wdims%j_start+1) *              &
                         (wdims%k_end-wdims%k_start+1) * n_proc)
  max_s_w         = MAXVAL(s_w(wdims%i_start:wdims%i_end,             &
                            wdims%j_start:wdims%j_end,                &
                            wdims%k_start:wdims%k_end))
  min_s_w         = MINVAL(s_w(wdims%i_start:wdims%i_end,             &
                            wdims%j_start:wdims%j_end,                &
                            wdims%k_start:wdims%k_end))
   av_s_w         = SUM(s_w(wdims%i_start:wdims%i_end,                &
                            wdims%j_start:wdims%j_end,                &
                            wdims%k_start:wdims%k_end))/fieldsize

  has_nan_s_w = 0.
!  IF(ANY(s_w.ne.s_w)) has_nan_s_w = 1.

  fieldsize       = REAL((tdims%i_end-tdims%i_start+1) *              &
                         (tdims%j_end-tdims%j_start+1) *              &
                         (tdims%k_end-tdims%k_start+1) * n_proc)
  max_s_theta     = MAXVAL(s_thetav(tdims%i_start:tdims%i_end,        &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
  min_s_theta     = MINVAL(s_thetav(tdims%i_start:tdims%i_end,        &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
   av_s_theta     = SUM(s_thetav(tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))/fieldsize
  has_nan_s_theta = 0.
!  IF(ANY(s_thetav.ne.s_thetav)) has_nan_s_theta = 1.


  max_s_m_v       = MAXVAL(s_m_v(tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
  min_s_m_v       = MINVAL(s_m_v(tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
   av_s_m_v       = SUM(s_m_v  (tdims%i_start:tdims%i_end,            &
                                tdims%j_start:tdims%j_end,            &
                                tdims%k_start:tdims%k_end))/fieldsize
  has_nan_s_m_v = 0.
!  IF(ANY(s_m_v.ne.s_m_v)) has_nan_s_m_v = 1.

  max_s_m_cl      = MAXVAL(s_m_cl(tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
  min_s_m_cl      = MINVAL(s_m_cl(tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
   av_s_m_cl      = SUM(s_m_cl (tdims%i_start:tdims%i_end,            &
                                tdims%j_start:tdims%j_end,            &
                                tdims%k_start:tdims%k_end))/fieldsize
  has_nan_s_m_cl = 0.
!  IF(ANY(s_m_cl.ne.s_m_cl)) has_nan_s_m_cl = 1.


  max_s_m_cf      = MAXVAL(s_m_cf(tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
  min_s_m_cf      = MINVAL(s_m_cf(tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
   av_s_m_cf      = SUM(s_m_cf (tdims%i_start:tdims%i_end,            &
                                tdims%j_start:tdims%j_end,            &
                                tdims%k_start:tdims%k_end))/fieldsize
  has_nan_s_m_cf = 0.
!  IF(ANY(s_m_cf.ne.s_m_cf)) has_nan_s_m_cf = 1.


  IF(l_mcr_qcf2) THEN
    max_s_m_cf2   = MAXVAL(s_m_cf2(tdims%i_start:tdims%i_end,         &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
    min_s_m_cf2   = MINVAL(s_m_cf2(tdims%i_start:tdims%i_end,         &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
     av_s_m_cf2   = SUM(s_m_cf2(tdims%i_start:tdims%i_end,            &
                                tdims%j_start:tdims%j_end,            &
                                tdims%k_start:tdims%k_end))/fieldsize
    has_nan_s_m_cf2 = 0.
!    IF(ANY(s_m_cf2.ne.s_m_cf2)) has_nan_s_m_cf2 = 1.

  END IF

  IF(l_mcr_qrain) THEN
    max_s_m_rain  = MAXVAL(s_m_r(tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
    min_s_m_rain  = MINVAL(s_m_r(tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
     av_s_m_rain  = SUM(s_m_r  (tdims%i_start:tdims%i_end,            &
                                tdims%j_start:tdims%j_end,            &
                                tdims%k_start:tdims%k_end))/fieldsize
    has_nan_s_m_rain = 0.
!    IF(ANY(s_m_r.ne.s_m_r)) has_nan_s_m_rain = 1.

  END IF

  IF(l_mcr_qgraup) THEN
    max_s_m_graup = MAXVAL(s_m_gr(tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
    min_s_m_graup = MINVAL(s_m_gr(tdims%i_start:tdims%i_end,          &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))
     av_s_m_graup = SUM(s_m_gr (tdims%i_start:tdims%i_end,            &
                                tdims%j_start:tdims%j_end,            &
                                tdims%k_start:tdims%k_end))/fieldsize
    has_nan_s_m_graup = 0.
!    IF(ANY(s_m_gr.ne.s_m_gr)) has_nan_s_m_graup = 1.

  END IF

  CALL gc_rmax(1,n_proc,ierr,max_s_u)
  CALL gc_rmin(1,n_proc,ierr,min_s_u)
  CALL gc_rmax(1,n_proc,ierr,has_nan_s_u)
  CALL gc_rsum(1,n_proc,ierr, av_s_u)
  CALL gc_rmax(1,n_proc,ierr,max_s_v)
  CALL gc_rmin(1,n_proc,ierr,min_s_v)
  CALL gc_rmax(1,n_proc,ierr,has_nan_s_v)
  CALL gc_rsum(1,n_proc,ierr, av_s_v)
  CALL gc_rmax(1,n_proc,ierr,max_s_w)
  CALL gc_rmin(1,n_proc,ierr,min_s_w)
  CALL gc_rmax(1,n_proc,ierr,has_nan_s_w)
  CALL gc_rsum(1,n_proc,ierr, av_s_w)
  CALL gc_rmax(1,n_proc,ierr,max_s_theta)
  CALL gc_rmin(1,n_proc,ierr,min_s_theta)
  CALL gc_rmax(1,n_proc,ierr,has_nan_s_theta)
  CALL gc_rsum(1,n_proc,ierr, av_s_theta)
  CALL gc_rmax(1,n_proc,ierr,max_s_m_v)
  CALL gc_rmin(1,n_proc,ierr,min_s_m_v)
  CALL gc_rmax(1,n_proc,ierr,has_nan_s_m_v)
  CALL gc_rsum(1,n_proc,ierr, av_s_m_v)
  CALL gc_rmax(1,n_proc,ierr,max_s_m_cl)
  CALL gc_rmin(1,n_proc,ierr,min_s_m_cl)
  CALL gc_rmax(1,n_proc,ierr,has_nan_s_m_cl)
  CALL gc_rsum(1,n_proc,ierr, av_s_m_cl)
  CALL gc_rmax(1,n_proc,ierr,max_s_m_cf)
  CALL gc_rmin(1,n_proc,ierr,min_s_m_cf)
  CALL gc_rmax(1,n_proc,ierr,has_nan_s_m_cf)
  CALL gc_rsum(1,n_proc,ierr, av_s_m_cf)
  IF(l_mcr_qcf2) THEN
    CALL gc_rmax(1,n_proc,ierr,max_s_m_cf2)
    CALL gc_rmin(1,n_proc,ierr,min_s_m_cf2)
    CALL gc_rmax(1,n_proc,ierr,has_nan_s_m_cf2)
    CALL gc_rsum(1,n_proc,ierr, av_s_m_cf2)
  END IF
  IF(l_mcr_qrain) THEN
    CALL gc_rmax(1,n_proc,ierr,max_s_m_rain)
    CALL gc_rmin(1,n_proc,ierr,min_s_m_rain)
    CALL gc_rmax(1,n_proc,ierr,has_nan_s_m_rain)
    CALL gc_rsum(1,n_proc,ierr, av_s_m_rain)
  END IF
  IF(l_mcr_qgraup) THEN
    CALL gc_rmax(1,n_proc,ierr,max_s_m_graup)
    CALL gc_rmin(1,n_proc,ierr,min_s_m_graup)
    CALL gc_rmax(1,n_proc,ierr,has_nan_s_m_graup)
    CALL gc_rsum(1,n_proc,ierr, av_s_m_graup)
  END IF

  IF(me.eq.0) THEN
    WRITE(6,fmt='(A)')   '********************************************'//&
                         '********************************************'//&
                         '****'
    WRITE(6,fmt='(A)')   'Fast physics sources for ENDGame from atmos_'//&
                         'physics2:'
    WRITE(6,fmt='(A)')   '             min                      max   '//&
                         '                   average (non-bit reproduc'//&
                         'ing  (1=has NaN 0= no NaN)'
    WRITE(6,fmt='(A,4E32.16)') 's_u      :',min_s_u    ,max_s_u            &
                                      , av_s_u    ,has_nan_s_u 
    WRITE(6,fmt='(A,4E32.16)') 's_v      :',min_s_v    ,max_s_v            &
                                      , av_s_v    ,has_nan_s_v 
    WRITE(6,fmt='(A,4E32.16)') 's_w      :',min_s_w    ,max_s_w            &
                                      , av_s_w    ,has_nan_s_w 
    WRITE(6,fmt='(A,4E32.16)') 's_thetav :',min_s_theta,max_s_theta        &
                                      , av_s_theta    ,has_nan_s_theta
    WRITE(6,fmt='(A,4E32.16)') 's_m_v    :',min_s_m_v  ,max_s_m_v          &
                                      , av_s_m_v    ,has_nan_s_m_v 
    WRITE(6,fmt='(A,4E32.16)') 's_m_cl   :',min_s_m_cl ,max_s_m_cl         &
                                      , av_s_m_cl    ,has_nan_s_m_cl
    WRITE(6,fmt='(A,4E32.16)') 's_m_cf   :',min_s_m_cf ,max_s_m_cf         &
                                      , av_s_m_cf    ,has_nan_s_m_cf
    IF(l_mcr_qcf2)  WRITE(6,fmt='(A,4E32.16)') 's_m_cf2  :',min_s_m_cf2    &
                      ,max_s_m_cf2    , av_s_m_cf2    ,has_nan_s_m_cf2
    IF(l_mcr_qgraup)WRITE(6,fmt='(A,4E32.16)') 's_m_graup:',min_s_m_graup  &
                  ,max_s_m_graup  , av_s_m_graup    ,has_nan_s_m_graup
    IF(l_mcr_qrain) WRITE(6,fmt='(A,4E32.16)') 's_m_rain :',min_s_m_rain   &
                    ,max_s_m_rain   , av_s_m_rain    ,has_nan_s_m_rain
    WRITE(6,fmt='(A)')   '********************************************'//&
                         '********************************************'//&
                         '****'
  END IF
END IF

IF (integrity_test)                                                   &
  CALL update_hash_m(S_u,           SIZE(S_u),             's_u__',   &
                     S_v,           SIZE(S_v),             's_v__',   &
                     S_w,           SIZE(S_w),             's_w__',   &
                     S_thetav,      SIZE(S_thetav),        's_the',   &
                     S_m_v,         SIZE(S_m_v),           's_m_v',   &
                     S_m_cl,        SIZE(S_m_cl),          's_mcl',   &
                     S_m_cf,        SIZE(S_m_cf),          's_mcf',   &
                     S_m_r,         SIZE(S_m_r),           's_m_r',   &
                     S_m_gr,        SIZE(S_m_gr),          's_mgr',   &
                     S_m_cf2,       SIZE(S_m_cf2),         'smcf2')

IF (lhook) CALL dr_hook('EG_R_S',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_r_s
END MODULE eg_r_s_mod
