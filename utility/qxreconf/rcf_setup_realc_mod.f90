! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up the real constants in the output header

Module Rcf_setup_RealC_Mod

!  Subroutine Rcf_Setup_RealC - initialises real constants in header.
!
! Description:
!   Sets the real constants in the output dump header.
!
! Method:
!    Uses namelist information to set header.
!    UMDP F3 defines the real constants.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Setup_RealC( Hdr_In, Hdr_Out )

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Rcf_UMhead_Mod, Only : &
    Um_Header_type

Use rcf_headers_mod, Only : &
    RelHd

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

USE rcf_model_mod, ONLY : &
    delta_lon, delta_lat, &
    frstlona,  frstlata,  &
    polelona,  polelata

Use Rcf_HeadAddress_Mod, Only :                    &
    RC_PoleLong,                  RC_LongSpacing,  &
    RC_LatSpacing,                RC_PoleLat,      &
    RC_FirstLong,                 RC_FirstLat,     &
    RC_ModelTop,                  RC_SWLDEG,       &
    RC_WEdgeDeg

Implicit None

! Arguments
Type (Um_Header_Type), Intent( In ) :: Hdr_In
Type (Um_Header_Type), Target       :: Hdr_Out

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

! Local vars
Real, Pointer                  :: RealC(:)
Integer                        :: i

!---------------------------------------------------------------
! Clean start - RMDI for everything
!---------------------------------------------------------------
RealC => Hdr_Out % RealC
RealC(:) = RMDI

!--------------------------------------------------------------
! Set values we have numbers for
!--------------------------------------------------------------

RealC( RC_LongSpacing ) = delta_lon
RealC( RC_LatSpacing  ) = delta_lat
RealC( RC_FirstLat    ) = frstlata
RealC( RC_FirstLong   ) = frstlona
RealC( RC_PoleLat     ) = polelata
RealC( RC_PoleLong    ) = polelona

! Submodel specifics
RealC( RC_ModelTop    ) = Output_Grid % z_top_of_model

!--------------------------------------------------------------
! Overrides from namelists
!--------------------------------------------------------------
Do i = 1, Hdr_Out % LenRealC
  If( RelHd(i) /= RMDI ) Then
    If (PrintStatus >= PrStatus_Oper) Then
      Write (6,*) 'RealC(',i,') has been reset from ', RealC(i), &
                  ' to ', RelHd(i)
    End If

    RealC(i) = RelHd(i)
  End If
End Do

Return
End Subroutine Rcf_Setup_RealC
End Module Rcf_Setup_RealC_Mod
