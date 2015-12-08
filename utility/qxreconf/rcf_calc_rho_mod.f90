! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Calculates RHO for the output dump

Module Rcf_Calc_Rho_Mod

!  Subroutine Rcf_Calc_Rho - calculates tho
!
! Description:
!   Calculates RHO for the 5.0/5.1 dump (rho should not be interpolated
!   so as to maintain dynamical balance)
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Calc_Rho( theta, q, exner, p, theta_heights,        &
                         rho_heights, rho )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Field_Equals_Mod, Only : &
    Rcf_Field_Equals

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE UM_ParVars, Only : &
    mype

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

Use Rcf_Calc_Exner_Theta_Mod, Only : &
    Rcf_Calc_Exner_Theta

USE atmos_constants_mod, ONLY: r, repsilon

Implicit None

! Arguments
Type( field_type ), Intent( In )     :: theta
Type( field_type ), Intent( In )     :: q
Type( field_type ), Intent( In )     :: p           ! on rho levels
Type( field_type ), Intent( In )     :: exner       ! on rho levels
Type( field_type ), Intent( InOut )  :: rho

Real, Intent(In)                   :: rho_heights( :, 0: )
Real, Intent(In)                   :: theta_heights( :, 0: )

! Local variables
Type( field_type )                 :: exner_theta ! on theta levels
Integer                            :: i
Integer                            :: k
Integer                            :: k_off
Integer                            :: k_rho
Real                               :: weight1
Real                               :: weight2
Real                               :: weight3
Real                               :: temp
Real                               :: work_real ( theta % level_size )
Real                               :: work_real2( theta % level_size )

!-------------------------------------------------------------------
! Need to calculate exner on theta levels
!-------------------------------------------------------------------
Call Rcf_field_equals( exner_theta, exner )
exner_theta % levels = theta % levels
Call Rcf_Alloc_Field( exner_theta )

IF ( theta % bottom_level == 0 ) THEN
  k_off = 1
  ! Have to handle zeroth level theta differently
  Call Rcf_calc_exner_theta( exner % level_size, exner_theta % levels-1, &
                             theta_heights, rho_heights(:,1:),          &
                             exner % Data, exner_theta % Data(:,2:) )
  exner_theta % Data (:,1) = exner_theta % Data (:,2)
ELSE
  k_off = 0
  Call Rcf_calc_exner_theta( exner % level_size, exner_theta % levels, &
                             theta_heights, rho_heights(:,1:),        &
                             exner % Data, exner_theta % Data )
END IF

!--------------------------------------------------------------------
! Now do the calculation of rho
!--------------------------------------------------------------------

If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
  Write (6,*) 'Calculating Rho'
End If

k = 1+k_off
Do i = 1, theta % level_size
  k_rho = k-k_off
! calculate thetav
  work_real(i) = theta % Data(i,k) * (1. +                             &
                                   (1./repsilon -1.) * q % Data(i,k) )
! calculate rho
  rho % Data(i,k_rho) = rho_heights(i,k_rho) * rho_heights(i,k_rho) * &
                 p % Data(i,k_rho) / (R * work_real(i) * exner % Data(i,k_rho))
End Do

Do k = 2+k_off, theta % levels
  k_rho = k-k_off
  Do i = 1, theta % level_size
    work_real2(i) = work_real(i)
  End Do

  If (k .le. q % levels ) Then
    Do i = 1, theta % level_size
    work_real(i) = theta % Data(i,k) * (1. +                           &
                                   (1./repsilon -1.) * q % Data(i,k) )
    End Do

  Else
    Do i = 1, theta % level_size
      work_real(i) = theta % Data(i,k)
    End Do
  End If

  If (k .ne. theta % levels) Then
    Do i = 1, theta % level_size
      weight1 = rho_heights(i,k_rho) - theta_heights(i,k_rho-1)
      weight2 = theta_heights(i,k_rho) - rho_heights(i,k_rho)
      weight3 = theta_heights(i,k_rho) - theta_heights(i,k_rho-1)

!      temp = ( weight2 * work_real(i) +             &
!               weight1 * work_real2(i) ) /          &
!               weight3

      temp = ( weight1 * work_real(i) +             &
               weight2 * work_real2(i) ) /          &
               weight3

      rho % Data(i,k_rho) = rho_heights(i,k_rho) * rho_heights(i,k_rho)   &
                      * p % Data(i,k_rho) / (R * temp * exner % Data(i,k_rho))
    End Do

  Else

    Do i= 1, theta % level_size

      temp = work_real2(i) *                                   &
             exner_theta % Data(i, exner_theta % levels - 1) / &
             exner % Data(i,k_rho)

      rho % Data(i,k_rho) = rho_heights(i,k_rho) * rho_heights(i,k_rho) *      &
                        p % Data(i,k_rho) / (R * temp * exner % Data(i,k_rho) )

    End Do

  End If
End Do

!------------------------------------------------------------------
! Tidy up
!------------------------------------------------------------------
Call Rcf_DeAlloc_Field( exner_theta )

Return
End Subroutine Rcf_Calc_Rho
End Module Rcf_Calc_Rho_Mod
