! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to read requested fields from a UM fieldsfile.
!=======================================================================

SUBROUTINE GetB_Flds ( NumFlds,     &  ! in
                       MaxFlds,     &  ! in
                       STCode,      &  ! in
                       MO8Level,    &  ! in
                       FCTime,      &  ! in
                       MinsPastHr,  &  ! in
                       MinsDif,     &  ! in 
                       Store,       &  ! in
                       PPHdrMod,    &  ! in
                       B_Hdr,       &  ! in
                       UMHdr,       &  ! in
                       Fields,      &  ! inout
                       ErrorStatus )   ! inout

! Description:
!   This subroutine reads in a
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumFlds
INTEGER, INTENT(IN) :: MaxFlds
INTEGER, INTENT(IN) :: STCode  (NumFlds)
INTEGER, INTENT(IN) :: MO8Level(NumFlds)
INTEGER, INTENT(IN) :: FCTime  (NumFlds)
INTEGER, INTENT(IN) :: Store   (NumFlds)
INTEGER, INTENT(IN) :: MinsPastHr (NumFlds)
INTEGER, INTENT(IN) :: MinsDif    (NumFlds) 
LOGICAL, INTENT(IN) :: PPHdrMod
TYPE(PP_Header_type), INTENT(IN) :: B_Hdr   ! Template for B-grid
TYPE(UM_Header_type), INTENT(IN) :: UMHdr

TYPE(PP_Field_type), INTENT(INOUT) :: Fields(MaxFlds)
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "GetB_Flds"

! Local Variables:
INTEGER :: i                      ! local counter
INTEGER :: ifld
CHARACTER(LEN=80) :: ErrMessage
TYPE(PP_Field_type) :: TempField(1)

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

Nullify ( TempField(1) % RData )

! Attempt to store all the fields in the fields array.
! DEPENDS ON: getflds
CALL GetFlds( NumFlds, MaxFlds, -STCode, MO8Level, FCTime,         &
              (/ (-1,i=1,NumFlds) /),                              &
              MinsPastHr, MinsDif, Store, PPHdrMod, UMHdr, Fields, &
              ErrorStatus )
IF ( ErrorStatus /= StatusOK ) THEN
  GO TO 9999
END IF

!-------------------------------
! Loop through required fields
!-------------------------------
DO i = 1,NumFlds
  ifld = Store(i)     ! Array index where field should be stored

  IF ( Fields(ifld) % Hdr % STcode == -STCode(i) ) THEN
    ! Already found previously interpolated field otherwise it would be a
    ! positive number.
    CYCLE
  END IF

  ALLOCATE(TempField(1) % Rdata(Fields(ifld) % Hdr % NumCols, &
                                Fields(ifld) % Hdr % NumRows)) 

! Retrieve the appropriate C-grid field to TempField
  TempField(1) % Hdr       = Fields(ifld) % Hdr
  TempField(1) % RData     = Fields(ifld) % Rdata
  TempField(1) % ArrayPos  = Fields(ifld) % Arraypos
  TempField(1) % LookupPos = Fields(ifld) % Lookuppos

! DEPENDS ON: ctobgrid
  CALL CtoBgrid( TempField(1), B_Hdr, Fields(ifld), ErrorStatus )

  DEALLOCATE( TempField(1) % RData )
  NULLIFY( TempField(1) % RData )

  IF ( ErrorStatus /= StatusOK ) THEN
    GO TO 9999
  END IF

END DO

9999 CONTINUE
! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE GetB_Flds

