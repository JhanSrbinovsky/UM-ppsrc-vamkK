! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_thetav_theta_mod

IMPLICIT NONE

CONTAINS
SUBROUTINE eg_thetav_theta                                            &
                  (thetav,theta,mix_v                                 &
                  ,p, pstar, p_theta_levels                           &
                  ,exner, exner_star, exner_theta_levels              &
                   )

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE atmos_constants_mod, ONLY: recip_epsilon, recip_kappa, p_zero
USE integrity_mod

USE atm_fields_bounds_mod
USE horiz_grid_mod, ONLY : intw_rho2w

IMPLICIT NONE
!
! Description:
!          Convert from potential virtual dry temperature to
!          potential temperature and calculate pressure from exner
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

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER                                                               &
  row_length                                                          &
                  ! number of points on a row
, rows                                                                &
                  ! number of rows of data
, model_levels                                                        &
, offx                                                                &
, offy

REAL, INTENT (IN) ::                                                  &
       mix_v(tdims_s%i_start:tdims_s%i_end,                           &
             tdims_s%j_start:tdims_s%j_end,                           &
             tdims_s%k_start:tdims_s%k_end),                          &
      thetav(tdims_s%i_start:tdims_s%i_end,                           &
             tdims_s%j_start:tdims_s%j_end,                           &
             tdims_s%k_start:tdims_s%k_end),                          &
  exner_star(pdims_s%i_start:pdims_s%i_end,                           &
             pdims_s%j_start:pdims_s%j_end)

REAL ::           exner(pdims_s%i_start:pdims_s%i_end,                &
                        pdims_s%j_start:pdims_s%j_end,                &
                        pdims_s%k_start:pdims_s%k_end+1),             &
                      p(pdims_s%i_start:pdims_s%i_end,                &
                        pdims_s%j_start:pdims_s%j_end,                &
                        pdims_s%k_start:pdims_s%k_end+1),             &
                  pstar(pdims%i_start:pdims%i_end,                    &
                        pdims%j_start:pdims%j_end),                   &
         p_theta_levels(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end),               &
     exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

REAL, INTENT (OUT) ::                                                 &
                theta  (tdims_s%i_start:tdims_s%i_end,                &
                        tdims_s%j_start:tdims_s%j_end,                &
                        tdims_s%k_start:tdims_s%k_end)

INTEGER :: i, j, k 

REAL :: inv_epsilon
REAL :: rkp


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook('EG_THETAV_THETA',zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(pdims_s,p,exner)
DO  k = pdims_s%k_start, pdims_s%k_end+1
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
       p(i,j,k) = p_zero*exner(i,j,k)**recip_kappa
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(tdims_s,theta,thetav,mix_v)
DO  k = tdims_s%k_start, tdims_s%k_end
  DO j = tdims_s%j_start, tdims_s%j_end
    DO i = tdims_s%i_start, tdims_s%i_end
       theta (i,j,k) = thetav(i,j,k) /                                &
                       (1.0+ mix_v (i,j,k)*recip_epsilon)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! Exner at surface

DO j = pdims_s%j_start, pdims_s%j_end
  DO i = pdims_s%i_start, pdims_s%i_end
      exner_theta_levels(i,j,tdims_s%k_start) = exner_star(i,j)
   END DO
END DO

DO  k = tdims_s%k_start+1,tdims_s%k_end
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
      exner_theta_levels(i,j,k) = (intw_rho2w(k,2)*exner(i,j,k)+      &
                                   intw_rho2w(k,1)*exner(i,j,k+1))
    END DO
  END DO
END DO

! Convert to pressure

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
     pstar(i,j) = p_zero*exner_star(i,j)**recip_kappa
     p_theta_levels(i,j,tdims_s%k_start) = pstar(i,j)
  END DO
END DO

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(tdims_s,pdims_s,p_theta_levels,exner_theta_levels)
DO  k = tdims_s%k_start+1, tdims_s%k_end
  DO j = pdims_s%j_start, pdims_s%j_end
    DO i = pdims_s%i_start, pdims_s%i_end
    p_theta_levels(i,j,k) = p_zero*exner_theta_levels(i,j,k)**recip_kappa
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (INTEGRITY_TEST)                                                   &
  CALL update_hash_m( exner,           SIZE(exner),           'pi___')

! end of routine

IF (lhook) CALL dr_hook('EG_THETAV_THETA',zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_thetav_theta
END MODULE eg_thetav_theta_mod
