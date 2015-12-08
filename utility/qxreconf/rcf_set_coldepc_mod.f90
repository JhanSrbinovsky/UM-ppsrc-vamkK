! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ sets up the output dump column dependent constants

Module Rcf_Setup_ColDepC_Mod

!  Subroutine Rcf_Setup_ColDepC - sets up the output dump
!                                 column dependent contants.
!
! Description:
!   The column dependent constants for the output dump are constructed.
!
! Method:
!   The column dependent constants are constructed from the output
!   grid details - see UMDP F3 for full  details.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains
Subroutine Rcf_Setup_ColDepC( Hdr_Out, Grid )

Use Rcf_Grid_Type_Mod, Only : &
    Grid_Type

Use Rcf_UMhead_Mod, Only : &
    UM_Header_Type

Use Rcf_Headaddress_Mod, Only :  &
    CDC_Lambda_input_p, CDC_Lambda_input_u

Use Rcf_readnl_horizont_mod, Only : &
    Lambda_input_p, Lambda_input_u

Implicit None

! Arguments
Type( UM_Header_Type ), Intent(InOut) :: Hdr_Out
Type( Grid_Type ), Intent(In)         :: Grid      ! Output Grid

! Local Variables
Integer                         :: i        ! Looper

  Do i = 1,  Hdr_Out % Len1ColDepC
     Hdr_Out % ColDepC( i, CDC_Lambda_input_p) =  Lambda_input_p(i)
     Hdr_Out % ColDepC( i, CDC_Lambda_input_u) =  Lambda_input_u(i)
  End Do

Return
End Subroutine Rcf_Setup_ColDepC
End Module Rcf_Setup_ColDepC_Mod
