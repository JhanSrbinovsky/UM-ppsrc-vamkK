! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Stochastic Physics (sect35) Random Parameters Ver. 2
MODULE stph_rp_pert_mod

IMPLICIT NONE

CONTAINS


 SUBROUTINE stph_rp_pert(l_firstcall, rp_rand0, rp_rand, rp_0, rp_max,  &
                         rp_min, rp_pert)

 USE stochastic_physics_run_mod, ONLY:                                  &
     ran_max, ran_count
     
 USE yomhook, ONLY: lhook, dr_hook           
 USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
!
! Description:
!   Performs calculation of new perturbed value of parameter passed
!   from stph_rp2
!
! Method:
!   Increment perturbation using AR1 process
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Stochastic Physics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

!
! IN variables
!
 REAL, INTENT(IN)    :: rp_max            ! Max value of random parameter
 REAL, INTENT(IN)    :: rp_min            ! Min value of random parameter
 REAL, INTENT(IN)    :: rp_rand0(ran_max) ! Random number for firstcall
 REAL, INTENT(IN)    :: rp_rand(ran_max)  ! Random number
 LOGICAL, INTENT(IN) :: l_firstcall       ! 1st call to RP?
!
! IN/OUT variables
!
 REAL, INTENT(INOUT) :: rp_pert           ! Perturbed parameter value
 REAL, INTENT(INOUT) :: rp_0              ! Default (deterministic value)
!
! Local variables
!
 REAL                :: shock             ! shock term in perturbation
 REAL                :: maxran, minran    ! max and min values of paramter
!
! Variables associated with autoregression model
!
 REAL, PARAMETER     :: corr=0.95         ! Correlation value for the
                                          ! first-order autoregression

 INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
 INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
 REAL(KIND=jprb)               :: zhook_handle

 IF (lhook) CALL dr_hook('STPH_RP_PERT',zhook_in,zhook_handle)

!----------------------------------------------------------------------
! On first call set UMUI default field and randomise initial pert value
!----------------------------------------------------------------------
 IF (l_firstcall) THEN
   rp_0 = rp_pert
   rp_pert = rp_rand0(ran_count) * (rp_max - rp_min) + rp_min
 END IF

!----------------------------------------------------------------------
! Perturb parameter using AR1 process
!----------------------------------------------------------------------
 maxran = (rp_max - rp_min)/3.0
 minran = SIGN(maxran, -1.0)
 shock  = rp_rand(ran_count) * (maxran - minran) + minran

 rp_pert = rp_0 + ((rp_pert - rp_0)*corr) + shock
 rp_pert = MAX( rp_pert, rp_min )
 rp_pert = MIN( rp_pert, rp_max )

! Increment random array counter (limit by maximum size)
 ran_count = ran_count + 1
 ran_count = MIN(ran_count, ran_max)

 IF (lhook) CALL dr_hook('STPH_RP_PERT',zhook_out,zhook_handle)
 RETURN
END SUBROUTINE stph_rp_pert
END MODULE stph_rp_pert_mod
