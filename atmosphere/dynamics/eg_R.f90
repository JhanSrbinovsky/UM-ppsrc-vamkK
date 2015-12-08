! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_r_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE eg_r(theta_star, m_star, mcl_star,mcf_star,mgraup_star,    &
                mrain_star, mcf2_star ,m_v, theta, l_physics )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atmos_constants_mod, ONLY : epsln=>repsilon
USE proc_info_mod,       ONLY : me,n_proc
USE atm_fields_bounds_mod
USE integrity_mod
USE field_types

USE PrintStatus_mod

USE domain_params
USE eg_swap_bounds_mod

USE fields_rhs_mod,  ONLY:  r_m_v ,r_m_cl  ,r_m_cf  ,r_m_gr,       &
                            r_m_r  ,r_m_cf2,r_theta

USE um_input_control_mod,  ONLY: model_domain
USE mphys_inputs_mod,      ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
USE science_fixes_mod,     ONLY: l_exclude_fix_theta_source


IMPLICIT NONE
!
! Description: compute the slow physics source terms for the dynamics
!  
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


REAL ::      theta_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::      theta     (tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::          m_star(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL ::             m_v(tdims_s%i_start:tdims_s%i_end,                &
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

REAL inv_epsln

INTEGER :: rows, row_length, model_levels, offx, offy

INTEGER ierr

! local diagnostics
REAL max_r_theta  ,min_r_theta,   av_r_theta,   has_nan_r_theta
REAL max_r_m_v    ,min_r_m_v,     av_r_m_v,     has_nan_r_m_v
REAL max_r_m_cl   ,min_r_m_cl,    av_r_m_cl,    has_nan_r_m_cl
REAL max_r_m_cf   ,min_r_m_cf,    av_r_m_cf,    has_nan_r_m_cf
REAL max_r_m_cf2  ,min_r_m_cf2,   av_r_m_cf2,   has_nan_r_m_cf2
REAL max_r_m_rain ,min_r_m_rain,  av_r_m_rain,  has_nan_r_m_rain
REAL max_r_m_graup,min_r_m_graup, av_r_m_graup, has_nan_r_m_graup

REAL fieldsize

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_R',zhook_in, zhook_handle)

inv_epsln = 1./epsln

r_theta   = 0.
r_m_v     = 0.
r_m_cl    = 0.
r_m_cf    = 0.
r_m_cf2   = 0.
r_m_r     = 0.
r_m_gr    = 0.

IF (.NOT.l_physics) THEN
  IF (lhook) CALL dr_hook('EG_R',zhook_out,zhook_handle)
  RETURN
END IF

!
! atmos_physics1 returns increments.
! So no conversion is required (apart from dry -> virtual dry)
! 

IF (.NOT.l_exclude_fix_theta_source) THEN
      r_theta(tdims%i_start:tdims%i_end,                              &
              tdims%j_start:tdims%j_end,                              &
              tdims%k_start:tdims%k_end) =                            &
   theta_star(tdims%i_start:tdims%i_end,                              &
              tdims%j_start:tdims%j_end,                              &
              tdims%k_start:tdims%k_end) *                            &
(1. + (m_star(tdims%i_start:tdims%i_end,                              &
              tdims%j_start:tdims%j_end,                              &
              tdims%k_start:tdims%k_end) +                            &
           m_v(tdims%i_start:tdims%i_end,                             &
               tdims%j_start:tdims%j_end,                             &
               tdims%k_start:tdims%k_end))*inv_epsln) +               &
        theta(tdims%i_start:tdims%i_end,                              &
              tdims%j_start:tdims%j_end,                              &
              tdims%k_start:tdims%k_end) *                            &
       m_star(tdims%i_start:tdims%i_end,                              &
              tdims%j_start:tdims%j_end,                              &
              tdims%k_start:tdims%k_end)*inv_epsln
ELSE
      r_theta(tdims%i_start:tdims%i_end,                              &
              tdims%j_start:tdims%j_end,                              &
              tdims%k_start:tdims%k_end) =                            &
   theta_star(tdims%i_start:tdims%i_end,                              &
              tdims%j_start:tdims%j_end,                              &
              tdims%k_start:tdims%k_end) *                            &
(1. + (m_star(tdims%i_start:tdims%i_end,                              &
              tdims%j_start:tdims%j_end,                              &
              tdims%k_start:tdims%k_end) +                            &
           m_v(tdims%i_start:tdims%i_end,                             &
               tdims%j_start:tdims%j_end,                             &
               tdims%k_start:tdims%k_end))*inv_epsln)
END IF

         r_m_v(tdims%i_start:tdims%i_end,                             &
               tdims%j_start:tdims%j_end,                             &
               tdims%k_start:tdims%k_end) =                           &
        m_star(tdims%i_start:tdims%i_end,                             &
               tdims%j_start:tdims%j_end,                             &
               tdims%k_start:tdims%k_end)

        r_m_cl(tdims%i_start:tdims%i_end,                             &
               tdims%j_start:tdims%j_end,                             &
               tdims%k_start:tdims%k_end) =                           &
      mcl_star(tdims%i_start:tdims%i_end,                             &
               tdims%j_start:tdims%j_end,                             &
               tdims%k_start:tdims%k_end)

        r_m_cf(tdims%i_start:tdims%i_end,                             &
               tdims%j_start:tdims%j_end,                             &
               tdims%k_start:tdims%k_end) =                           &
      mcf_star(tdims%i_start:tdims%i_end,                             &
               tdims%j_start:tdims%j_end,                             &
               tdims%k_start:tdims%k_end)

IF (l_mcr_qcf2)                                                       &
    r_m_cf2(tdims%i_start:tdims%i_end,                                &
            tdims%j_start:tdims%j_end,                                &
            tdims%k_start:tdims%k_end) =                              &
  mcf2_star(tdims%i_start:tdims%i_end,                                &
            tdims%j_start:tdims%j_end,                                &
            tdims%k_start:tdims%k_end)

IF (l_mcr_qrain)                                                      &
   r_m_r   (tdims%i_start:tdims%i_end,                                &
            tdims%j_start:tdims%j_end,                                &
            tdims%k_start:tdims%k_end) =                              &
 mrain_star(tdims%i_start:tdims%i_end,                                &
            tdims%j_start:tdims%j_end,                                &
            tdims%k_start:tdims%k_end)

IF (l_mcr_qgraup)                                                     &
  r_m_gr   (tdims%i_start:tdims%i_end,                                &
            tdims%j_start:tdims%j_end,                                &
            tdims%k_start:tdims%k_end) =                              &
mgraup_star(tdims%i_start:tdims%i_end,                                &
            tdims%j_start:tdims%j_end,                                &
            tdims%k_start:tdims%k_end)


CALL eg_swap_bounds(r_theta ,tdims_s,fld_type_p,.FALSE.) 
CALL eg_swap_bounds(r_m_v   ,tdims_s,fld_type_p,.FALSE.) 
CALL eg_swap_bounds(r_m_cl  ,tdims_s,fld_type_p,.FALSE.) 
CALL eg_swap_bounds(r_m_cf  ,tdims_s,fld_type_p,.FALSE.) 
IF (l_mcr_qcf2   ) CALL eg_swap_bounds(r_m_cf2,tdims_s,fld_type_p,.FALSE.) 
IF (l_mcr_qrain  ) CALL eg_swap_bounds(r_m_r     ,tdims_s,fld_type_p,.FALSE.) 
IF (l_mcr_qgraup ) CALL eg_swap_bounds(r_m_gr    ,tdims_s,fld_type_p,.FALSE.) 


IF (printstatus >=  prstatus_diag) THEN

  fieldsize       = REAL((tdims%i_end-tdims%i_start+1) *              &
                         (tdims%j_end-tdims%j_start+1) *              &
                         (tdims%k_end-tdims%k_start+1) * n_proc)

  max_r_theta     = MAXVAL(r_theta)
  min_r_theta     = MINVAL(r_theta)
   av_r_theta     = SUM(r_theta (tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))          &
                                 /fieldsize
  has_nan_r_theta = 0.
  IF(ANY(r_theta.ne.r_theta)) has_nan_r_theta = 1.

  max_r_m_v       = MAXVAL(r_m_v)
  min_r_m_v       = MINVAL(r_m_v)
   av_r_m_v       = SUM(r_m_v  (tdims%i_start:tdims%i_end,            &
                                tdims%j_start:tdims%j_end,            &
                                tdims%k_start:tdims%k_end))           &
                                /fieldsize
  has_nan_r_m_v = 0.
  IF(ANY(r_m_v.ne.r_m_v)) has_nan_r_m_v = 1.

  max_r_m_cl      = MAXVAL(r_m_cl)
  min_r_m_cl      = MINVAL(r_m_cl)
   av_r_m_cl      = SUM(r_m_cl (tdims%i_start:tdims%i_end,            &
                                tdims%j_start:tdims%j_end,            &
                                tdims%k_start:tdims%k_end))           &
                                /fieldsize
  has_nan_r_m_cl = 0.
  IF(ANY(r_m_cl.ne.r_m_cl)) has_nan_r_m_cl = 1.

  max_r_m_cf      = MAXVAL(r_m_cf)
  min_r_m_cf      = MINVAL(r_m_cf)
   av_r_m_cf      = SUM(r_m_cf (tdims%i_start:tdims%i_end,            &
                                tdims%j_start:tdims%j_end,            &
                                tdims%k_start:tdims%k_end))           &
                                /fieldsize
  has_nan_r_m_cf = 0.
  IF(ANY(r_m_cf.ne.r_m_cf)) has_nan_r_m_cf = 1.

  IF(l_mcr_qcf2) THEN
    max_r_m_cf2   = MAXVAL(r_m_cf2)
    min_r_m_cf2   = MINVAL(r_m_cf2)
     av_r_m_cf2   = SUM(r_m_cf2(tdims%i_start:tdims%i_end,            &
                                tdims%j_start:tdims%j_end,            &
                                tdims%k_start:tdims%k_end))           &
                                 /fieldsize
  has_nan_r_m_cf2 = 0.
  IF(ANY(r_m_cf2.ne.r_m_cf2)) has_nan_r_m_cf2 = 1.

  END IF
  IF(l_mcr_qrain) THEN
    max_r_m_rain  = MAXVAL(r_m_r   )
    min_r_m_rain  = MINVAL(r_m_r   )
     av_r_m_rain  = SUM(r_m_r   (tdims%i_start:tdims%i_end,           &
                                 tdims%j_start:tdims%j_end,           &
                                 tdims%k_start:tdims%k_end))          &
                                 /fieldsize
  has_nan_r_m_rain = 0.
  IF(ANY(r_m_r.ne.r_m_r)) has_nan_r_m_rain = 1.

  END IF
  IF(l_mcr_qgraup) THEN
    max_r_m_graup = MAXVAL(r_m_gr   )
    min_r_m_graup = MINVAL(r_m_gr   )
     av_r_m_graup = SUM(r_m_gr    (tdims%i_start:tdims%i_end,         &
                                   tdims%j_start:tdims%j_end,         &
                                   tdims%k_start:tdims%k_end))        &
                                   /fieldsize
  has_nan_r_m_graup = 0.
  IF(ANY(r_m_gr.ne.r_m_gr)) has_nan_r_m_graup = 1.

  END IF

  CALL gc_rmax(1,n_proc,ierr,max_r_theta)
  CALL gc_rmin(1,n_proc,ierr,min_r_theta)
  CALL gc_rmax(1,n_proc,ierr,has_nan_r_theta)
  CALL gc_rsum(1,n_proc,ierr, av_r_theta)
  CALL gc_rmax(1,n_proc,ierr,max_r_m_v)
  CALL gc_rmin(1,n_proc,ierr,min_r_m_v)
  CALL gc_rmax(1,n_proc,ierr,has_nan_r_m_v)
  CALL gc_rsum(1,n_proc,ierr, av_r_m_v)
  CALL gc_rmax(1,n_proc,ierr,max_r_m_cl)
  CALL gc_rmin(1,n_proc,ierr,min_r_m_cl)
  CALL gc_rmax(1,n_proc,ierr,has_nan_r_m_cl)
  CALL gc_rsum(1,n_proc,ierr, av_r_m_cl)
  CALL gc_rmax(1,n_proc,ierr,max_r_m_cf)
  CALL gc_rmin(1,n_proc,ierr,min_r_m_cf)
  CALL gc_rmax(1,n_proc,ierr,has_nan_r_m_cf)
  CALL gc_rsum(1,n_proc,ierr, av_r_m_cf)
  IF(l_mcr_qcf2) THEN
    CALL gc_rmax(1,n_proc,ierr,max_r_m_cf2)
    CALL gc_rmin(1,n_proc,ierr,min_r_m_cf2)
    CALL gc_rmax(1,n_proc,ierr,has_nan_r_m_cf2)
    CALL gc_rsum(1,n_proc,ierr, av_r_m_cf2)
  END IF
  IF(l_mcr_qrain) THEN
    CALL gc_rmax(1,n_proc,ierr,max_r_m_rain)
    CALL gc_rmin(1,n_proc,ierr,min_r_m_rain)
    CALL gc_rmax(1,n_proc,ierr,has_nan_r_m_rain)
    CALL gc_rsum(1,n_proc,ierr, av_r_m_rain)
  END IF
  IF(l_mcr_qgraup) THEN
    CALL gc_rmax(1,n_proc,ierr,max_r_m_graup)
    CALL gc_rmin(1,n_proc,ierr,min_r_m_graup)
    CALL gc_rmax(1,n_proc,ierr,has_nan_r_m_graup)
    CALL gc_rsum(1,n_proc,ierr, av_r_m_graup)
  END IF

  IF(me.eq.0) THEN
    WRITE(6,fmt='(A,4E32.16)') 'r_thetav :',min_r_theta,max_r_theta,       &
                                            av_r_theta,has_nan_r_theta
    WRITE(6,fmt='(A,4E32.16)') 'r_m_v    :',min_r_m_v  ,max_r_m_v,         &
                                              av_r_m_v ,has_nan_r_m_v
    WRITE(6,fmt='(A,4E32.16)') 'r_m_cl   :',min_r_m_cl ,max_r_m_cl,        &
                                              av_r_m_cl,has_nan_r_m_cl
    WRITE(6,fmt='(A,4E32.16)') 'r_m_cf   :',min_r_m_cf ,max_r_m_cf,        &
                                               av_r_m_cf,has_nan_r_m_cf
    IF(l_mcr_qcf2)  WRITE(6,fmt='(A,4E32.16)') 'r_m_cf2  :',min_r_m_cf2,   &
                                                       max_r_m_cf2,   &
                                             av_r_m_cf2,has_nan_r_m_cf2
    IF(l_mcr_qgraup)WRITE(6,fmt='(A,4E32.16)') 'r_m_graup:',min_r_m_graup, &
                                                       max_r_m_graup ,&
                                         av_r_m_graup,has_nan_r_m_graup
    IF(l_mcr_qrain) WRITE(6,fmt='(A,4E32.16)') 'r_m_rain :',min_r_m_rain , &
                                                       max_r_m_rain , &
                                            av_r_m_rain,has_nan_r_m_rain

    WRITE(6,fmt='(A)') '============================================', &
                       '========================================'
   END IF
END IF

IF (integrity_test)                                                  &
  CALL update_hash_m(r_theta,        SIZE(r_theta),         'r_the', &
                     r_m_v,          SIZE(r_m_v),           'r_m_v', &
                     r_m_cl,         SIZE(r_m_cl),          'r_mcl', &
                     r_m_cf,         SIZE(r_m_cf),          'r_mcf', &
                     r_m_r,          SIZE(r_m_r),           'r_m_r', &
                     r_m_gr,         SIZE(r_m_gr),          'r_mgr', &
                     r_m_cf2,        SIZE(r_m_cf2),         'rmcf2')

IF (lhook) CALL dr_hook('EG_R',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_r
END MODULE eg_r_mod
