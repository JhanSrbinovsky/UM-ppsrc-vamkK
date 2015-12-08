! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Convert geopotential field to Orography

Module Rcf_Grib_Spcl_Orog_Mod

! SUBROUTINE Rcf_Grib_Spcl_Orog
!
! Description: This routine handles the conversion of Geopotential to
!              Orography.
!              At present there is a double check to ensure the data
!              came from ECMWF but if all centers return geopotential
!              instead of Orography then this check can be removed.
!
! Method: Divide every entry in the data field by G, the gravitational
!         constant.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
contains
Subroutine Rcf_Grib_Spcl_Orog(Current,FieldData)

Use Rcf_GRIB_Block_Params_Mod, Only : &
    Grib_Record,     &
    LenArrayMax,     &
    p_Orig_cntr

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

USE earth_constants_mod, ONLY: g

IMPLICIT NONE


! Subroutine arguments

!< Array  arguments with intent(InOut):>
Type (Grib_Record),Pointer          :: Current
Real, Intent(InOut)                 :: FieldData(LenArrayMax)

! Local constants
Character (Len=*), Parameter     :: RoutineName='Rcf_GRIB_Spcl_Orog'

! Local variables

Character (Len=80)               :: Cmessage(2)   ! used for EReport
Integer                          :: ErrorStatus   ! used for EReport

!=======================================================================
!  Routine Code Start :
!=======================================================================

! double check it's an ECMWF field (do other sources spec geopot ?)
If (Current % Block_1(p_Orig_cntr) == GrbOrigECMWF) Then

  If ( PrintStatus >= PrStatus_Diag  ) Then
    If ( mype == 0 ) Then
      Write (6,'(A)') "Converting Geopotential to Orography"
    End If
  End If
  FieldData(:) = FieldData(:) / G
Else
  Write (Cmessage(1),*) 'Data did not come from ECMWF'
  ErrorStatus = 10
  Call EReport( RoutineName, ErrorStatus, Cmessage(1) )
End If

Return

End Subroutine Rcf_Grib_Spcl_Orog
End Module Rcf_Grib_Spcl_Orog_Mod
