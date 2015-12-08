! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ takes thv and mse as inputs and works out theta and qv
!

SUBROUTINE dts_deduce_thetaandqv(n_dp,nlev,k_ad,q,mse,newthv,newqhyd,newmse &
                                 ,theta,wmse,wthv,wthv_smth,                &
                                 newqvap,newtheta)

USE atmos_constants_mod, ONLY: cp, c_virtual
USE water_constants_mod, ONLY: lc
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

!-----------------------------------------------------------------------
! Description:
!
! solves equations: thv = theta +0.61 theta qv - theta qhyd =
!                                               theta(1+0.61 qvap-qhyd)
!                   mse = theta + lvap/cp qvap
! for theta, qv via a quadratic (-b +/- sqrt(b^2-4ac))/2a
!
! Inputs:
! -------
! theta_v (= newthv)
! q_hyd (= newqhyd)
! mse (= newmse)
!
! Outputs:
! --------
! theta (= newtheta)
! qvap (= newqvap)
!
! Other information:
! ------------------
! called from FLUX CALCULATION SECTION of deep_turb_conv
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Convection

! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP 3 programming standards version vn8.2.
!-----------------------------------------------------------------------

! Subroutine Arguments

INTEGER, INTENT(IN) ::      &
  n_dp                      & ! No. of deep convection points
 ,nlev                      & ! No. of model layers
 ,k_ad(n_dp)

REAL, INTENT(IN)    ::      &
  newthv(n_dp,nlev)         &
 ,newqhyd(n_dp,nlev)        &
 ,newmse(n_dp,nlev)         &
 ,q(n_dp,nlev)              &
 ,mse(n_dp,nlev)            &
 ,wmse(n_dp,nlev)           &
 ,wthv(n_dp,nlev)           &
 ,wthv_smth(n_dp,nlev)      &
 ,theta(n_dp,nlev)          
             
REAL, INTENT(OUT)    ::     &
  newqvap(n_dp,nlev)        &
 ,newtheta(n_dp,nlev)       


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
INTEGER ::         &
  i_dp,k              ! Loop counters

REAL ::            &
  a(n_dp,nlev)     &  ! a in solution of quadratic
 ,b(n_dp,nlev)     &  ! b in solution of quadratic 
 ,c(n_dp,nlev)     &  ! c in solution of quadratic
 ,alpha            &
 ,wq(n_dp,nlev)       ! Calulated but not actually used   

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

!This is a quadratic...:
IF (lhook) CALL dr_hook('DTS_DEDUCE_THETAANDQV',zhook_in,zhook_handle)
  newqvap(:,:) =  q(:,:) ! initialise in case can't solve

  a(:,:) = c_virtual*lc/cp 

!nb removing dependence on newqhyd 27/1/09
!        b(:,:) = (lc/cp)*(1.0-newqhyd(:,:)) - c_virtual*newmse(:,:)
!        c(:,:) = newthv(:,:)-newmse(:,:)*(1.0-newqhyd(:,:))

  b(:,:) = (lc/cp) - c_virtual*newmse(:,:)
  c(:,:) = newthv(:,:)-newmse(:,:)
        
  DO k=1,nlev
    DO i_dp=1,n_dp
!NB putting in a fudge that ensures discrepancies between wh and wthv
! don't lead to unphysical increments to the water vapour term 

      IF(b(i_dp,k)**2 >= 4.0*a(i_dp,k)*c(i_dp,k)) THEN 
        newqvap(i_dp,k) = (- b(i_dp,k) + sqrt(b(i_dp,k)**2 -          &
                             4.0*a(i_dp,k)*c(i_dp,k)))/               &
                              (2.0*a(i_dp,k)) 
        ! because 4ac is <0, neg option gives neg qvap

      END IF

      ! apportion the same fraction of the mse as was on the
      ! previous timestep to the new water vapour NB this is
      ! a dodgy fudge to get around tiny differences in the
      ! buoyancy flux leading to big q increments...
      ! may lead to supersaturation...
      IF(newqvap(i_dp,k) < 0.0) THEN 
        alpha = (lc/cp)*q(i_dp,k)/mse(i_dp,k)
              
        IF(alpha > 0.0) THEN 
          newqvap(i_dp,k) = alpha*newmse(i_dp,k)/(lc/cp)
        ELSE
          newqvap(i_dp,k) = 0.0
        END IF
      END IF ! newqvap < 0

      newtheta(i_dp,k) = newmse(i_dp,k) - (lc/cp)*newqvap(i_dp,k)

    END DO
  END DO

! for checking purposes only at this stage: calculate wq
  wq(:,:) = (wmse(:,:)*(1+0.6*q(:,:))-(wthv(:,:)+wthv_smth(:,:)) )&
             /(lc/cp*(1.0+c_virtual*q(:,:))-0.6*theta(:,:))

IF (lhook) CALL dr_hook('DTS_DEDUCE_THETAANDQV',zhook_out,zhook_handle)
RETURN

END SUBROUTINE dts_deduce_thetaandqv
