! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Convert volumetric soil moisture to soil moisture amount

Module Rcf_Grib_Spcl_SoilM_Mod

! SUBROUTINE Rcf_Grib_Spcl_SoilM
!
! Description: This routine handles the conversion of volumetric
!              soil moisture into soil moisture amounts.
!              At present there is a double check to ensure the data
!              came from ECMWF as not all input data will be volumetric.
!
! Method: Multiply volumetric quantity by density and layer depth.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
contains
Subroutine Rcf_Grib_Spcl_SoilM(Current,FieldData)

Use Rcf_GRIB_Block_Params_Mod, Only : &
    Grib_Record,     &
    LenArrayMax,     &
    p_Orig_cntr,     &
    p_Lvl_Desc_1,    &
    p_Lvl_Desc_2

Use rcf_GRIB_lookups_Mod, Only : &
    GrbOrigECMWF

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag

USE UM_ParVars, Only : &
    mype

Use EReport_Mod, Only :     &
    EReport

USE water_constants_mod, ONLY: rho_water

IMPLICIT NONE

! Subroutine arguments

!< Array  arguments with intent(InOut):>
Type (Grib_Record),Pointer          :: Current
Real, Intent(InOut)                 :: FieldData(LenArrayMax)

! Local constants
Character (Len=*), Parameter     :: RoutineName='Rcf_GRIB_Spcl_SoilM'

! Local variables

Character (Len=80)               :: Cmessage(2)   ! used for EReport
Integer                          :: ErrorStatus   ! used for EReport
Integer                          :: ilayer        ! Index of soil layer
Real                             :: layer_depth   ! Depth of soil layer


!=======================================================================
!  Routine Code Start :
!=======================================================================

! double check it's an ECMWF field
If (Current % Block_1(p_Orig_cntr) == GrbOrigECMWF) Then

  If ( PrintStatus >= PrStatus_Diag  ) Then
    If ( mype == 0 ) Then
      Write (6,'(A)') "Converting vol. SoilM to SoilM Amount"
    End If
  End If

  ! Hard wired for time being - need to derive these some how
  ! Note: Although Block_1(p_Lvl_Desc_1/2) have been modified
  !       in rcf_grib_check, the second read of the data
  !       overwrites the existing header information. Hence
  !       the values are as they are in the original grib file.
  Select Case (Current % Block_1(p_Lvl_Desc_1))
    Case (0)
      layer_depth = 0.10
    Case (7)
      layer_depth = 0.25
    Case (28)
      layer_depth = 0.65
    Case (100)
      layer_depth = 2.00
  End Select

  ! Multiply volumetric values by density and layer depth  
  FieldData(:) = layer_depth * rho_water * FieldData(:)

Else
  Write (Cmessage(1),*) 'Data did not come from ECMWF'
  ErrorStatus = 10
  Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
End If

Return

End Subroutine Rcf_Grib_Spcl_SoilM
End Module Rcf_Grib_Spcl_SoilM_Mod
