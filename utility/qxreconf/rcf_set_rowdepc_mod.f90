! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ sets up the output dump row dependent constants

Module Rcf_Setup_RowDepC_Mod

!  Subroutine Rcf_Setup_RowDepC - sets up the output dump
!                                 row dependent contants.
!
! Description:
!   The row dependent constants for the output dump are constructed.
!
! Method:
!   The row dependent constants are constructed from the output
!   grid details - see UMDP F3 for full  details.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains
Subroutine Rcf_Setup_RowDepC( Hdr_Out, Grid )

Use Rcf_Grid_Type_Mod, Only : &
    Grid_Type

Use Rcf_UMhead_Mod, Only : &
    UM_Header_Type

Use Rcf_Headaddress_Mod, Only :  &
    RDC_Phi_input_p, RDC_Phi_input_v

Use Rcf_readnl_horizont_mod, Only : &
    Phi_input_p, Phi_input_v
 
Implicit None

! Arguments
Type( UM_Header_Type ), Intent(InOut) :: Hdr_Out
Type( Grid_Type ), Intent(In)         :: Grid      ! Output Grid

! Local Variables
Integer                          :: i        ! Looper

  Do i = 1,  Hdr_Out % Len1RowDepC
     Hdr_Out % RowDepC( i, RDC_Phi_input_p) =  Phi_input_p(i)
     Hdr_Out % RowDepC( i, RDC_Phi_input_v) =  Phi_input_v(i)
  End Do

Return
End Subroutine Rcf_Setup_RowDepC
End Module Rcf_Setup_RowDepC_Mod
