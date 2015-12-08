! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE eg_dry_static_adj_mod
!
! Description:
!          Enforces static stability on input theta field.
!          Uses very simple sort algorithm resulting in no mixing of
!          temperature just a simple re-arrangement.
!  
!
! Method:
!          Adopted from ND version
!  
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section:  Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

CONTAINS

SUBROUTINE eg_dry_static_adj(theta,rho,exner,intw_w2rho)

USE atm_fields_bounds_mod, ONLY : pdims,pdims_s,tdims_s
USE atmos_constants_mod,   ONLY : p_zero,r,kappa

USE proc_info_mod,         ONLY : me
USE ereport_mod,           ONLY: ereport

USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook
USE integrity_mod

USE PrintStatus_mod

IMPLICIT NONE

! Subroutine arguments


REAL                                                                  &
  theta (tdims_s%i_start:tdims_s%i_end,                               &
         tdims_s%j_start:tdims_s%j_end,                               &
         tdims_s%k_start:tdims_s%k_end)

REAL                                                                  &
  rho  (pdims_s%i_start:pdims_s%i_end,                                &
        pdims_s%j_start:pdims_s%j_end,                                &
        pdims_s%k_start:pdims_s%k_end)

REAL                                                                  &
  exner  (pdims_s%i_start:pdims_s%i_end,                              &
          pdims_s%j_start:pdims_s%j_end,                              &
          pdims_s%k_start:pdims_s%k_end)

REAL, INTENT (IN) :: intw_w2rho(pdims_s%k_start:pdims_s%k_end,2)


! Local variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER                                                               &
  i,                                                                  &
  j,                                                                  &
  k,                                                                  &
  kk,                                                                 &
  list,                                                               &
                    ! Loop indices
  cnt,                                                                &
  count_old,                                                          &
  iterations,                                                         &
  n_search,                                                           &
  i_search,                                                           &
  n_match

REAL                                                                  &
  temp
! Local arrays

INTEGER                                                               &
  unstable    (pdims%i_end*pdims%j_end*pdims%k_end),                  &
  unstable_old(pdims%i_end*pdims%j_end*pdims%k_end),                  &
  i_prev      (pdims%i_end*pdims%j_end),                              &
  j_prev      (pdims%i_end*pdims%j_end)

INTEGER maxit                     ! Maximum of iteration for trying 
                                  ! to remove instabillity

INTEGER errorstatus

! No External Routines


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_DRY_STATIC_ADJ',zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 1. Check for inertial instability and remove.
!            Iterative sort procedure is used.
! ----------------------------------------------------------------------

maxit = 2*pdims%k_end

iterations = 0
cnt = 0

DO k = pdims%k_start,pdims%k_end
  DO j = pdims%j_start,pdims%j_end
    DO i = pdims%i_start,pdims%i_end

      IF (theta(i,j,k)  <   theta(i,j,k-1)) THEN
        temp = theta(i,j,k)
        theta(i,j,k) = theta(i,j,k-1)
        theta(i,j,k-1) = temp

        DO kk=MAX(1,k-1),k+1
!         recompute density to satisfy equation of state
          rho (i,j,kk) = p_zero/(R*                                   &
                          (intw_w2rho(kk,1)*theta(i,j,kk)             &
                         + intw_w2rho(kk,2)*theta(i,j,kk-1) ))        &
                         *(exner(i,j,kk))**((1.0-kappa)/kappa)
        END DO
        cnt = cnt + 1
        unstable(cnt) = k*1e6+j*1000 + i
      END IF
    END DO
  END DO
END DO

  IF ( PrintStatus == PrStatus_Diag .AND. me == 0  .AND. cnt  /=  0) THEN
    write(6,fmt='(A)') '================================================'
    write(6,fmt='(A)') 'static instabillity at (i,j,k,count):'
    write(6,fmt='(A)') '------------------------------------------------'

    DO list = 1, min(10,cnt)

      k = unstable(list) * 1e-6
      j = (unstable(list) - k * 1e6) *1e-3
      i = unstable(list) - k * 1e6 - j*1e3

      write(6,fmt='(4I4)') i,j,k,list

    END DO

    IF (cnt.gt.10) THEN
      write(6,fmt='(A,I7)')' further printing suppressed. total count:',cnt
    END IF

    write(6,fmt='(A)') '================================================'

  END IF

DO WHILE (cnt  >   0 .AND. iterations  <   maxit)
  count_old = cnt
  DO list = 1, cnt
    unstable_old(list) = unstable(list)
  END DO
  cnt = 0
  iterations = iterations + 1

  n_search = 0
  n_match = 0
  DO list = 1, count_old
    k = unstable_old(list) * 1e-6
    j = (unstable_old(list) - k * 1e6) *1e-3
    i = unstable_old(list) - k * 1e6 - j*1e3
    DO i_search = 1, n_search
      IF (i  ==  i_prev(i_search) .AND.                               &
          j  ==  j_prev(i_search)) THEN
        n_match = 1
      END IF
    END DO
    IF (n_match  ==  0) THEN
      n_search = n_search + 1
      i_prev(n_search) = i
      j_prev(n_search) = j
    END IF
    n_match = 0
  END DO

  DO list = 1, n_search
    i = i_prev(list)
    j = j_prev(list)
    DO k = 1, pdims%k_end
      IF (theta(i,j,k)  <   theta(i,j,k-1)) THEN

        temp = theta(i,j,k)
        theta(i,j,k) = theta(i,j,k-1)
        theta(i,j,k-1) = temp

        DO kk=MAX(1,k-1),k+1
!   recompute density to satisfy equation of state
          rho (i,j,kk) = p_zero/(R*                                   &
                          (intw_w2rho(kk,1)*theta(i,j,kk)             &
                         + intw_w2rho(kk,2)*theta(i,j,kk-1) ))        &
                         *(exner(i,j,kk))**((1.0-kappa)/kappa)
        END DO

        cnt = cnt + 1
        unstable(cnt) = k*1e6+j*1000 + i
      END IF
    END DO
  END DO

  IF ( PrintStatus == PrStatus_Diag .AND. me == 0  .AND. cnt  /=  0) THEN
    write(6,fmt='(A)') '================================================'
    write(6,fmt='(A)') 'static instabillity at (i,j,k,count):'
    write(6,fmt='(A)') '------------------------------------------------'

    DO list = 1, min(10,cnt)

      k = unstable(list) * 1e-6
      j = (unstable(list) - k * 1e6) *1e-3
      i = unstable(list) - k * 1e6 - j*1e3

      write(6,fmt='(4I4)') i,j,k,list

    END DO

    IF (cnt.gt.10) THEN
      write(6,fmt='(A,I7)')' further printing suppressed. total count:',cnt
    END IF

    write(6,fmt='(A)') '================================================'

  END IF

END DO

IF (iterations  ==  maxit .AND. cnt  /=  0) THEN
     errorstatus = 1
     Call ereport("eg_dry_static_adj",errorstatus ,                   &
                  " unable to remove static instability" )
END IF

IF (iterations > 0) THEN

  errorstatus = -iterations
  Call ereport("eg_dry_static_adj",errorstatus ,                      &
                  " static instability removed" )

END IF

IF (integrity_test)                                                   &
  CALL update_hash_m(                                                 &
                  rho,             SIZE(rho),             'r_np1',    &
                  exner,           SIZE(exner),           'pinp1',    &
                  theta,           SIZE(theta),           'tvnp1')

! End of routine
IF (lhook) CALL dr_hook('DRY_STATIC_ADJ',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_dry_static_adj
END MODULE eg_dry_static_adj_mod
