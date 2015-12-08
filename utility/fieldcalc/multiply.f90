! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE Multiply( NumFlds,       &  ! in
                     Fields1,       &  ! in
                     Fields2,       &  ! in
                     ProductFields, &  ! inout
                     ErrorStatus)      ! inout

! Description:
!   Multiplies two vectors of fields together.
!
! Method:
!   The output fields are initialised with the same headers as the first
!   vector of input fields. The stash codes are then adjusted and the
!   multiplication performed.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY : ereport

IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumFlds
TYPE(PP_Field_type), INTENT(IN)    :: Fields1(NumFlds)
TYPE(PP_Field_type), INTENT(IN)    :: Fields2(NumFlds)
TYPE(PP_Field_type), INTENT(INOUT) :: ProductFields(NumFlds)
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Multiply"

! Local Variables:
INTEGER :: i

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

DO i = 1, NumFlds

  IF ( (Fields1(i) % Hdr % NumCols /=                                 &
        Fields2(i) % Hdr % NumCols) .OR.                              &
       (Fields1(i) % Hdr % NumRows /=                                 &
        Fields2(i) % Hdr % NumRows) ) THEN
    ErrorStatus = StatusWarning

    CALL EReport( RoutineName, ErrorStatus,                           &
                  "Cannot multiply fields of different dimensions" )
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF

  IF ( ASSOCIATED( ProductFields(i) % RData ) .AND.                &
      (ProductFields(i) % ArrayPos /= Fields1(i) % ArrayPos) .AND. &
      (ProductFields(i) % ArrayPos /= Fields2(i) % ArrayPos) ) THEN
    DEALLOCATE( ProductFields(i) % RData )
    NULLIFY( ProductFields(i) % RData )
  END IF

  IF ( .NOT. ASSOCIATED( ProductFields(i) % RData ) ) THEN
    ALLOCATE( ProductFields(i) % RData(Fields1(i) % Hdr % NumCols,    &
                                       Fields1(i) % Hdr % NumRows) )
    ProductFields(i) % Hdr = Fields1(i) % Hdr
  END IF

END DO


WHERE ( (Fields1 % Hdr % STCode < 50000) )
  ProductFields % Hdr % STCode = Fields1 % Hdr % STCode + 50000
END WHERE

DO i = 1, NumFlds

  WHERE ( (Fields1(i) % RData /= Fields1(i) % Hdr % BMDI) .AND.       &
          (Fields2(i) % RData /= Fields2(i) % Hdr % BMDI) )
    ProductFields(i) % RData = Fields1(i) % RData * Fields2(i) % RData
  ELSEWHERE
    ProductFields(i) % RData = ProductFields(i) % Hdr % BMDI
  END WHERE

END DO

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE Multiply
