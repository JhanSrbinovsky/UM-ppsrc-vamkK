! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to calculate weather diagnostics

SUBROUTINE PrSym( DynRain,      &  ! in
                  ConRain,      &  ! in
                  DynSnow,      &  ! in
                  ConSnow,      &  ! in
                  SurfTemp,     &  ! in
                  SymCode,      &  ! inout
                  ErrorStatus )    ! inout

! Description:
!   Calculate precipation symbol.
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE conversions_mod, ONLY: zerodegc

USE IO_Mod, ONLY:         &
  PP_Field_type,          &
  PP_Header_type
USE Err_Mod, ONLY:        &
  StatusOK
USE FldCodes_Mod, ONLY:   &
  ST_PrSym,  MO8_PrSym,  PP_PrSym
IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: DynRain
TYPE(PP_Field_type), INTENT(IN) :: ConRain
TYPE(PP_Field_type), INTENT(IN) :: DynSnow
TYPE(PP_Field_type), INTENT(IN) :: ConSnow
TYPE(PP_Field_type), INTENT(IN) :: SurfTemp

TYPE(PP_Field_type), INTENT(INOUT) :: SymCode
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "PrSym"

! Local Variables:
REAL :: PrecipDyn (DynRain%Hdr%NumCols, DynRain%Hdr%NumRows) ! Total
REAL :: PrecipCon (DynRain%Hdr%NumCols, DynRain%Hdr%NumRows) !   precip.
REAL :: PrecipTot (DynRain%Hdr%NumCols, DynRain%Hdr%NumRows) !   rates
REAL :: PrecipSnow(DynRain%Hdr%NumCols, DynRain%Hdr%NumRows)

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( SymCode % RData ) ) THEN
  DEALLOCATE( SymCode % RData )
END IF
SymCode % Hdr = DynRain % Hdr
SymCode % Hdr % PPCode   =  PP_PrSym
SymCode % Hdr % MO8Type  = MO8_PrSym
SymCode % Hdr % MO8Level = 8888
SymCode % Hdr % STCode   =  ST_PrSym
ALLOCATE( SymCode % RData(DynRain % Hdr % NumCols, &
                          DynRain % Hdr % NumRows) )
SymCode % RData(:,:) = 0.0

! Set precip codes
PrecipDyn  = DynRain % RData + DynSnow % RData
PrecipCon  = ConRain % RData + ConSnow % RData
PrecipSnow = DynSnow % RData + ConSnow % RData
PrecipTot  = PrecipDyn + PrecipCon

WHERE( PrecipTot > 0.0 )
  SymCode % RData = PrecipTot + 1000.0
END WHERE
WHERE( PrecipSnow > 0.0 )
  SymCode % RData = PrecipTot + 2000.0
END WHERE
WHERE( (ConSnow % RData > 0.0) .AND. &
       (SurfTemp % RData > 2.0+ZeroDegC) )
  SymCode % RData = PrecipTot + 3000.0
END WHERE

WHERE( (PrecipTot > 0.0) .AND. (PrecipDyn > PrecipCon) )
  SymCode % RData = SymCode % RData + 10000.0
END WHERE
WHERE( (PrecipTot > 0.0) .AND. (PrecipDyn <= PrecipCon) )
  SymCode % RData = SymCode % RData + 20000.0
END WHERE

9999 CONTINUE

END SUBROUTINE PrSym

!=======================================================================


