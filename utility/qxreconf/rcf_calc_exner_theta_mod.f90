! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************

!+ Calculates exner at theta levels from values on rho levels.

Module Rcf_Calc_exner_theta_mod

!  Subroutine Rcf_Calc_exner_theta
!
! Description: Calculates exner at theta levels.
!
! Method: Identical to that previously used for pressure instead of
!         exner, see "A semi-Implicit scheme for the Unified Model".
!          F.R. Division working paper No 154.
!          M. J. P. Cullen and T. Davies.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   27/07/01   Original code.  S. Cusack
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_calc_exner_theta( &
      p_field,               &! Intent(IN) Number of points in field
      model_levels,          &! Intent(IN) Number of model levels
      r_theta_levels,        &! Intent(IN) Height at theta levels
      r_rho_levels,          &! Intent(IN) Height at rho levels
      exner_rho,             &! Intent(IN) Exner at rho levels
      exner_theta_levels     &! Intent(OUT) Exner at theta levels
                               )

Implicit None

! Arguments with Intent IN. ie: Input variables.
Integer, Intent(In)     :: p_field       ! Number of points in field
Integer, Intent(In)     :: model_levels  ! number of model levels.

Real, Intent(In)        :: exner_rho(p_field,model_levels+1)
Real, Intent(In)        :: r_theta_levels(p_field,0:model_levels)
Real, Intent(In)        :: r_rho_levels(p_field, model_levels)

! Arguments with Intent OUT. ie: Output variables.
Real, Intent(Out)       :: exner_theta_levels(p_field, model_levels)

! Local Variables.
Integer                 :: i, k      ! Loop indices

! --------------------------------------------------------------------
! Section 1.   Calculate pressure at desired theta levels.
! --------------------------------------------------------------------
Do k = 1, model_levels - 1
  Do i = 1, p_field
    exner_theta_levels(i,k) =                                         &
        (exner_rho(i,k) *(r_rho_levels(i,k+1) - r_theta_levels(i,k))  &
       + exner_rho(i,k+1) * (r_theta_levels(i,k) - r_rho_levels(i,k)))&
                         / (r_rho_levels(i,k+1) - r_rho_levels(i,k) )
  End Do
End Do


k = model_levels

!  extra pressure level above top theta level is same height above
!  as is the pressure below - hence weights are 0.5 and there is
!  no need to store the r_rho_level for the extra pressure

Do i = 1, p_field
  exner_theta_levels(i,k) = 0.5 * ( exner_rho(i,k) + exner_rho(i,k+1) )
End Do


! End of routine
Return
End Subroutine Rcf_Calc_exner_theta

End Module Rcf_Calc_exner_theta_mod
