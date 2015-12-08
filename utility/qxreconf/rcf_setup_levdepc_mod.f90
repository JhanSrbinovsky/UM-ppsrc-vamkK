! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ sets up the output dump level dependent constants

Module Rcf_Setup_LevDepC_Mod

!  Subroutine Rcf_Setup_LevDepC - sets up the output dump
!                                 level dependent contants.
!
! Description:
!   The level dependent constants for the output dump are constructed.
!
! Method:
!   The Level dependent constants are constructed from the output
!   grid details - see UMDP F3 for full  details.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Setup_LevDepC( Hdr_Out, Hdr_In, Grid )

Use Rcf_Grid_Type_Mod, Only : &
    Grid_Type

Use Rcf_UMhead_Mod, Only : &
    UM_Header_Type

Use Rcf_Generate_Heights_Mod, Only : &
    height_gen_smooth

Use Rcf_Headaddress_Mod, Only :            &
    LDC_EtaTheta,            LDC_EtaRho,   &
    LDC_RHCrit,              SoilDepths,   &
    LDC_ZseaTheta,           LDC_CkTheta,  &
    LDC_ZseaRho,             LDC_CkRho

Implicit None

! Arguments
Type( UM_Header_Type ), Intent(In)    :: Hdr_In
Type( UM_Header_Type ), Intent(InOut) :: Hdr_Out
Type( Grid_Type ), Intent(In)         :: Grid      ! Output Grid

! Local Variables
Integer                               :: i        ! Looper
Integer                               :: j        ! Looper
Integer                               :: Len2Min  ! Minimum length



! Eta Theta Values
Do i = 0, Grid % model_levels
  Hdr_Out % LevDepC(i+1,LDC_EtaTheta) = Grid % eta_theta_levels( i )
End Do

! Eta Rho Values
Do i = 1, Grid % model_levels
  Hdr_Out % LevDepC(i,LDC_EtaRho) = Grid % eta_rho_levels( i )
End Do

! RHCrit values
Do i = 1, Grid % model_levels
  Hdr_Out % LevDepC(i,LDC_RHCrit) = Grid % rhcrit(i)
End Do

! Soil Depths
Do i = 1, Grid % sm_levels
  Hdr_Out % LevDepC(i,SoilDepths) = Grid % soil_depths( i )
End Do

! Only need level dependent constants 5-8 if we are using the
! new height generation method
If ( Grid % Height_Gen_Method == height_gen_smooth ) Then

  ! Zsea values for Theta
  Do i = 0, Grid % model_levels
    Hdr_Out % LevDepC(i+1,LDC_ZseaTheta) = Grid % eta_theta_levels(i)&
                                         * Grid % z_top_of_model
  End Do

  ! Ck values for Theta
  Do i = 0, grid % first_constant_r_rho_level - 1
    Hdr_Out % LevDepC(i+1,LDC_CkTheta) =(1.0 -                       &
            grid % eta_theta_levels( i )/                            &
            grid % eta_rho_levels(grid % first_constant_r_rho_level))&
             ** 2
  End Do

  Do i = grid % first_constant_r_rho_level , grid % model_levels
    Hdr_Out % LevDepC(i+1,LDC_CkTheta) = 0.0
  End Do

  ! Zsea values for Rho
  Do i = 1, Grid % model_levels
    Hdr_Out % LevDepC(i,LDC_ZseaRho) = Grid % eta_rho_levels( i ) *  &
                                       Grid % z_top_of_model
  End Do

  ! Ck values for Rho
  Do i = 1, grid % first_constant_r_rho_level
    Hdr_Out % LevDepC(i,LDC_CkRho) = (1.0 -                          &
            grid % eta_rho_levels( i ) /                             &
            grid % eta_rho_levels(grid % first_constant_r_rho_level))&
            ** 2
  End Do

  Do i = grid % first_constant_r_rho_level + 1, grid % model_levels
    Hdr_Out % LevDepC(i,LDC_CkRho) = 0.0
  End Do

End If

Return
End Subroutine Rcf_Setup_LevDepC

End Module Rcf_Setup_LevDepC_Mod
