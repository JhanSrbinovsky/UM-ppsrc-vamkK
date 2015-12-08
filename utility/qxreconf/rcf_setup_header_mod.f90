! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ sets up the output dump header

Module Rcf_Setup_Header_Mod

!  Subroutine Rcf_Setup_Header - sets up the output dump header.
!
! Description:
!   The header for the output dump is constructed.
!
! Method:
!   The header is setup section by section - see UMDP F3 for
!   details.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Setup_Header( Hdr_In, Hdr_Out )

Use Rcf_UMhead_Mod, Only : &
    Um_header_type,    &
    LenFixHd

Use Rcf_Setup_FixHd_Mod, Only : &
    Rcf_Setup_FixHd

Use Rcf_Setup_IntC_Mod, Only : &
    Rcf_Setup_IntC

Use Rcf_Setup_RealC_Mod, Only : &
    Rcf_Setup_RealC

Use Rcf_Setup_Lookup_Mod, Only : &
    Rcf_Setup_Lookup

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_AllocHdr_Mod, Only : &
    Rcf_AllocHdr

Use Rcf_Setup_LevDepC_Mod, Only : &
    Rcf_Setup_LevDepC

Use Rcf_Setup_RowDepC_Mod, Only : &
    Rcf_Setup_RowDepC

Use Rcf_Setup_ColDepC_Mod, Only : &
    Rcf_Setup_ColDepC

Implicit None

! Arguments
Type (Um_Header_Type), Intent(In)    :: Hdr_In
Type (Um_Header_Type), Intent(InOut) :: Hdr_Out

! Comdecks
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! Local arguments
Integer                              :: i          ! looper
Integer                              :: Len1_min   ! Min(Len1)
Integer                              :: Len2_min   ! Min(Len2) 

!--------------------------------------------------------------
! Fixed Header
!--------------------------------------------------------------
Allocate( Hdr_Out % FixHd( LenFixHd ) )
Call Rcf_Setup_FixHd( Hdr_In, Hdr_Out )
Call Rcf_AllocHdr( Hdr_Out )

!--------------------------------------------------------------
! Integer Constants
!--------------------------------------------------------------
Call Rcf_Setup_IntC( Hdr_Out )

!--------------------------------------------------------------
! Real Constants
!--------------------------------------------------------------
Call Rcf_Setup_RealC( Hdr_In, Hdr_Out )

!-------------------------------------------------------------
! Below here things could be put into their own routines -
! I've not done this for simplicity, but it should be done if
! there is significant expansion of any of the sections.
!-------------------------------------------------------------

!-------------------------------------------------------------
! Level dependent constants
!-------------------------------------------------------------
! Initialise to RMDI
Hdr_Out % LevDepC(:,:) = RMDI

Call Rcf_Setup_LevDepC( Hdr_Out, Hdr_In, Output_Grid )

!--------------------------------------------------------------
! row dependent constants
!--------------------------------------------------------------
Call Rcf_Setup_RowDepC( Hdr_Out, Output_Grid )

!----------------------------------------------------------------
! Copy column dependent constants
!----------------------------------------------------------------
Call Rcf_Setup_ColDepC( Hdr_Out, Output_Grid ) 

!-----------------------------------------------------------------
! Copy fields of constants
!-----------------------------------------------------------------
Hdr_Out % FldsOfC( : ) = Hdr_In % FldsOfC( : )

!-----------------------------------------------------------------
! Copy extra constants
!-----------------------------------------------------------------
Do i = 1, Hdr_Out % LenExtraC
  Hdr_Out % ExtraC( i ) = Hdr_In % ExtraC( i )
End Do

!-----------------------------------------------------------------
! Copy History block
!-----------------------------------------------------------------
Hdr_Out % HistFile( : ) = Hdr_In % HistFile( : )

!-----------------------------------------------------------------
! Copy Compressed field indexes - not needed since ocean removed.
!-----------------------------------------------------------------

!--------------------------------------------------------------
! Lookup Tables
!--------------------------------------------------------------
 Call Rcf_Setup_Lookup( Hdr_In, Hdr_Out )

Return
End Subroutine Rcf_Setup_Header
End Module Rcf_Setup_Header_Mod
