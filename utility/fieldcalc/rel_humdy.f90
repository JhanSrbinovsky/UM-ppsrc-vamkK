! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate relative humidity

SUBROUTINE Rel_humdy( NumFlds,       &  ! in
                      Fields1,       &  ! in - specific humidity
                      Fields2,       &  ! in - saturated mixing ratio
                      ProductFields, &  ! inout - relative humidity
                      ErrorStatus)      ! inout

! Description:
!   Calculates the relative humidity given the specific humidity and SVP
!   as input
!
! Method:
!   The headers for the output fields are copied from the first vector of
!   input fields.
!   Calculates RH as:
!   RH = x/x_{sat} = (SH/(1-SH))/x_{sat}
!   where RH = relative humidity, SH = specific humidity, and
!         x_{sat} = saturated mixing ratio
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
IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumFlds                               ! No. levels
TYPE(PP_Field_type), INTENT(IN) :: Fields1(NumFlds)          ! specific humidity
TYPE(PP_Field_type), INTENT(IN) :: Fields2(NumFlds)          ! saturated mixing
                                                             !   ratio

TYPE(PP_Field_type), INTENT(INOUT) :: ProductFields(NumFlds) ! relative
                                                             !   humidity
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Rel_humdy"
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

! Local Variables:
INTEGER :: i                                                 ! Loop counter

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )


IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF


! Allocate output fields
DO i = 1, NumFlds

  IF ( (Fields1(i) % Hdr % NumCols /= Fields2(i) % Hdr % NumCols) .OR.     &
       (Fields1(i) % Hdr % NumRows /= Fields2(i) % Hdr % NumRows) ) THEN

    CALL EReport( RoutineName, ErrorStatus,                        &
                "Cannot multiply fields of different dimensions")
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF

  IF ( ASSOCIATED( ProductFields(i) % RData ) ) THEN
    DEALLOCATE( ProductFields(i) % RData )
    NULLIFY( ProductFields(i) % RData )
  END IF

  IF ( .NOT. ASSOCIATED( ProductFields(i) % RData ) ) THEN
    ALLOCATE( ProductFields(i) % RData(Fields1(i) % Hdr % NumCols,         &
                                       Fields1(i) % Hdr % NumRows) )
  END IF

END DO


! Set header
ProductFields % Hdr = Fields1 % Hdr


! Calculate relative humidity
DO i = 1, NumFlds

  WHERE ( (Fields1(i) % RData /= Fields1(i) % Hdr % BMDI) .AND.            &
          (Fields2(i) % RData /= Fields2(i) % Hdr % BMDI) )
    ProductFields(i) % RData = (Fields1(i) % RData / &
                               (1.0 - Fields1(i) % RData)) &
                               / Fields2(i) % RData
  ELSEWHERE
    ProductFields(i) % RData = ProductFields(i) % Hdr % BMDI
  END WHERE

END DO


9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE Rel_humdy
