! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Convert log(pstar) into pstar

Module Rcf_Grib_Spcl_LPstar_Mod

! SUBROUTINE Rcf_Grib_Spcl_LPstar
!
! Description: This routine handles the conversion of log(pstar) to
!              pstar.
!              At present there is a double check to ensure the data
!              came from ECMWF as log(pstar) is not a standard variable.
!
! Method: Take the exponential of every entry in the data field.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
contains
Subroutine Rcf_Grib_Spcl_LPstar(Current,FieldData)

Use Rcf_GRIB_Block_Params_Mod, Only : &
    Grib_Record,     &
    LenArrayMax,     &
    p_Orig_cntr,     &
    p_Param_ID,      &
    EID_Surf_Press

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

IMPLICIT NONE

! Subroutine arguments

!< Array  arguments with intent(InOut):>
Type (Grib_Record),Pointer          :: Current
Real, Intent(InOut)                 :: FieldData(LenArrayMax)

! Local constants
Character (Len=*), Parameter     :: RoutineName='Rcf_GRIB_Spcl_LPstar'

! Local variables

Character (Len=80)               :: Cmessage(2)   ! used for EReport
Integer                          :: ErrorStatus   ! used for EReport

!=======================================================================
!  Routine Code Start :
!=======================================================================

! double check it's an ECMWF field
If (Current % Block_1(p_Orig_cntr) == GrbOrigECMWF) Then

  If ( PrintStatus >= PrStatus_Diag  ) Then
    If ( mype == 0 ) Then
      Write (6,'(A)') "Converting log(pstar) to pstar"
    End If
  End If
  FieldData(:) = exp( FieldData(:) )
  Current % Block_1(p_Param_ID) = EID_Surf_Press
Else
  Write (Cmessage(1),*) 'Data did not come from ECMWF'
  ErrorStatus = 10
  Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
End If

Return

End Subroutine Rcf_Grib_Spcl_LPstar
End Module Rcf_Grib_Spcl_LPstar_Mod
