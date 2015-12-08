! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+Routine to calculate saturated mixing ratio


SUBROUTINE Sat_MRatio( NumFlds,       &  ! in
                       Fields1,       &  ! in
                       Fields2,       &  ! in
                       SatMR,         &  ! inout
                       ErrorStatus)      ! inout

! Description:
!   Calculates the saturated mixing ratio (for use in relative humidity
!   calculation) given temperature and pressure as input.
!
! Method:
!   The headers for the output fields are copied from the first vector of input
!   fields (i.e. pressure). This should be of no consequence unless x_sat is
!   being calculated in its own right, as in the relative humidity calculation
!   the headers are taken from the specific humidity field.
!
!   x_{sat} is calculated as:
!   x_{sat} = ( 0.622 * e_{s} ) / ( p - e_{s} )
!   where e_{s} is the saturated vapour pressure (SVP):
!         e_{s} = 6.11 * 10^( 7.5*T / (237.7 + T) )
!   using the Murray formulation, with T in deg C and p in hPa
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6


USE atmos_constants_mod, ONLY: repsilon

USE conversions_mod, ONLY: zerodegc

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumFlds
TYPE(PP_Field_type), INTENT(IN) :: Fields1(NumFlds)        ! pressure
TYPE(PP_Field_type), INTENT(IN) :: Fields2(NumFlds)        ! temperature

TYPE(PP_Field_type), INTENT(INOUT) :: SatMR(NumFlds)       ! saturated
                                                           !   mixing ratio
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Sat_MRatio"
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
INTEGER :: i                                               ! Loop counter

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

    CALL EReport( RoutineName, ErrorStatus,                                &
                  "Cannot multiply fields of different dimensions")
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF

  IF ( ASSOCIATED( SatMR(i) % RData ) ) THEN
    DEALLOCATE( SatMR(i) % RData )
    NULLIFY( SatMR(i) % RData )
  END IF

  IF ( .NOT. ASSOCIATED( SatMR(i) % RData ) ) THEN
    ALLOCATE( SatMR(i) % RData(Fields1(i) % Hdr % NumCols,                 &
                               Fields1(i) % Hdr % NumRows) )
  END IF

END DO


! Set header
SatMR % Hdr = Fields1 % Hdr


! Calculate saturated mixing ratio using the two equations described
! in the header combined into one. Temperature is converted into
! degrees C and pressure into hPa, hence in the equation
!   T = Fields2(i)%RData - 273.15
!   p = Fields1(i)%RData / 100.0
! Note also the part of the equation (237.7+T) has been reduced to
!   Fields2(i)%RData - 35.45
! which is equivalent to 237.7 + Fields2(i)%RData - 273.15

DO i = 1, NumFlds
  WHERE ( (Fields1(i) % RData /= Fields1(i) % Hdr % BMDI) .AND.  &
          (Fields2(i) % RData /= Fields2(i) % Hdr % BMDI) )

    SatMR(i) % RData = ( repsilon * &
                       ( 6.11 * &
                       ( 10.0**( (7.5 * (Fields2(i) % RData - ZeroDegC)) &
                                      / (Fields2(i) % RData - 35.45) ) )  ) &
!
                        /  &
!
                        ( (Fields1(i) % RData / 100.0) - &
                          ( 6.11 * &
                          ( 10.0**( (7.5 * (Fields2(i) % RData - ZeroDegC)) &
                                      / (Fields2(i) % RData - 35.45) ) ) ) ) )


  ELSEWHERE
    SatMR(i) % RData = SatMR(i) % Hdr % BMDI
  END WHERE

END DO


9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE Sat_MRatio
