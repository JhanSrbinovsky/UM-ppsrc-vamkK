! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

Module Rcf_Calc_2D_CCA_Mod

!  Subroutine Rcf_Calc_2D_CCA - calculates 2D cca from 3D cca
!
! Description: Calculates a 2D convective cloud amount from the 2D
!              convective cloud amount array.
!
! Method: If an anvil is detected (ie there is more convective cloud
!         at the top than the base), the base value is divided by
!         the tower_factor. Otherwise, the value is the one at
!         the cloud base. Note that the tower_factor used to genenerate
!         the 3D field is not available, so the most common one
!         in practice is used.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.3   30/08/01   Original code.  P.Selwood
!   6.0   11/09/03   Change argument of Tiny(min_val) to Tiny(0.1)
!                    for IBM (provided by Zoe Chaplin)   P.Dando
!   6.2   29/07/05   Extended change of Tiny(min_val) to Portland
!                    GNU/Linux compiler. T.Edwards
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Reconfiguration

Contains

Subroutine Rcf_Calc_2D_CCA( cca_3d, ccb, cct, cca_2d )

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Implicit None

! Arguments
Type( field_type ), Intent(In)    :: cca_3d
Type( field_type ), Intent(In)    :: ccb
Type( field_type ), Intent(In)    :: cct
Type( field_type ), Intent(InOut) :: cca_2d

! Comdecks

! Local variables
Integer                           :: i          ! Looper

! The most commonly used value for tower_factor is set as a parameter
! as it is not supplied to the reconfiguration
Real, Parameter                   :: tower_factor = 0.25

! min_val is the amount of difference in convective cloud which
! constitutes an anvil.



Real, Parameter                   :: min_val = Tiny(min_val)


! max_2d_cca is the maximum amount of 2d cca permitted
Real, Parameter                   :: max_2d_cca = 0.5

! Loop over all points
Do i=1, cca_3d % level_size

  ! Do we have any cca to convert?
  If ( ccb % Data_Int(i,1) == 0 .AND. cct % Data_Int(i,1) == 0) Then

    cca_2d % Data(i,1) = 0.0

  ! Do we have an anvil?
  Else If ((cca_3d % Data( i, cct % Data_Int(i,1) - 1) -        &
       cca_3d % Data( i, ccb % Data_Int(i,1))) >  min_val) Then
    cca_2d % Data(i,1) = cca_3d % Data(i, ccb % Data_Int(i,1))  &
                       / tower_factor

  Else
    cca_2d % Data(i,1) = cca_3d % Data(i, ccb % Data_Int(i,1))
  End If

  ! Ensure that the maximum 2d cca isn't breached.
  If (cca_2d % Data(i,1) > max_2d_cca) Then
    cca_2d % Data(i,1) = max_2d_cca
  End If

End Do

Return
End Subroutine Rcf_Calc_2D_CCA
End Module Rcf_Calc_2D_CCA_Mod
