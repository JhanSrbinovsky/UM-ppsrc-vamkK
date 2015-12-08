! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up the integer constants for output dump

Module Rcf_Setup_IntC_Mod

!  Subroutine Rcf_Setup_IntC - initialisation of Integer Constants
!
! Description:
!   Sets up the output dump integer constants in the header.
!
! Method:
!   Uses namelist variables to define.
!   UMDP F3 defines the integer constants.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

Contains

Subroutine Rcf_Setup_IntC( Hdr_Out )

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_UMhead_Mod, Only : &
    um_header_type

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Rcf_HeadAddress_Mod, Only : &
    IC_XLen,              IC_YLen,            IC_PLevels,       &
    IC_WetLevels,         IC_SoilTlevels,     IC_NoCloudLevels, &
    IC_NoSeaPts,          IC_TracerLevs,      IC_BLevels,       &
    IC_NumLandPoints,     IC_NumOzoneLevs,    IC_ConvectLevs,   &
    IC_MDI,               IC_SoilMoistLevs,   IC_1stConstRho,   &
    IC_TracerVars,        IC_HeightMethod,    IC_RiverRows,     &
    IC_RiverRowLength

Use rcf_headers_Mod, Only : &
    inthd

USE river_inputs_mod, ONLY : &
    l_rivers

USE nlsizes_namelist_mod, ONLY : &
    tr_vars

Implicit None

! Arguments
Type (Um_header_type),Target     :: Hdr_Out

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

! Local variables
Integer, Pointer         :: IntC(:)
Integer                  :: i


!------------------------------------------------------------------
! Start with clean sheet - IMDI all around
!------------------------------------------------------------------
IntC => Hdr_Out % IntC
IntC(:) = IMDI

!------------------------------------------------------------------
! Set values we have numbers for
!------------------------------------------------------------------
IntC( IC_XLen        )   = Output_Grid % glob_p_row_length
IntC( IC_YLen        )   = Output_Grid % glob_p_rows
IntC( IC_PLevels     )   = Output_Grid % model_levels

IntC( IC_WetLevels   )   = Output_Grid % wet_levels
IntC( IC_SoilTLevels )   = Output_Grid % st_levels
IntC( IC_NoCloudLevels ) = Output_Grid % cloud_levels
IntC( IC_TracerVars    ) = tr_vars
IntC( IC_BLevels       ) = Output_Grid % bl_levels
IntC( IC_NumOzoneLevs  ) = Output_Grid % ozone_levels
IntC( IC_SoilMoistLevs ) = Output_Grid % sm_levels
IntC( IC_1stConstRho   ) = Output_Grid % first_constant_r_rho_level
IntC( IC_ConvectLevs   ) = Output_Grid % conv_levels
IntC( IC_NumLandPoints ) = Output_Grid % glob_land_field
IntC( IC_HeightMethod  ) = Output_Grid % height_gen_method

If ( L_RIVERS ) Then
  IntC( IC_RiverRows     )  = Output_Grid % glob_r_rows
  IntC( IC_RiverRowLength)  = Output_Grid % glob_r_row_length
Else
  IntC( IC_RiverRows     )  = IMDI
  IntC( IC_RiverRowLength)  = IMDI
End If

IntC( IC_TracerLevs    ) = Output_Grid % tr_levels
IntC( IC_MDI           ) = IMDI

! Note AMIP/other IntC (35 +) isn't set yet.
! Note also that IntC(15) - number of land tile types - has not
! been set

!------------------------------------------------------------------
! Namelist overrides
!------------------------------------------------------------------
Do i = 1, Hdr_Out % LenIntC
  If ( IntHd(i) /= IMDI ) Then
    If ( PrintStatus >= PrStatus_Oper ) Then
      Write (6,*) 'IntC(',i,') has been reset from ', IntC(i), &
                  ' to ', IntHd(i)
    End If

    IntC( i ) = IntHd( i )
  End If
End Do

Return
End Subroutine Rcf_Setup_IntC
End Module Rcf_Setup_IntC_Mod
