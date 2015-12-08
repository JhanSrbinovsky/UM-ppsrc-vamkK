! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE ICOnPLyr( NumLevs ,      &  ! in
                     NumLyrs,       &  ! in
                     LyrMO8L,       &  ! in
                     LyrLwrB,       &  ! in
                     LyrUprB,       &  ! in
                     LevFields ,    &  ! in
                     PFields ,      &  ! in
                     LyrFields,     &  ! inout
                     ErrorStatus)      ! inout

! Description:
!   Calculates a vector of pressure layer fields of icing potential
!   from a vector of model level fields of icing potential.
!
! Method:
!   Two vectors of pressures are supplied.  These contain the lower and
!   upper boundaries of a set of layers of the atmosphere.  For each of
!   these layers the corresponding data on the model level fields is
!   found using pressures on model levels supplied.  Two fields for
!   every layer are output by the subroutine.  These are contained in
!   the field array LyrFields. To each gridpoint in the fields
!   comprising the first half of the output the mean of the
!   corresponding data is assigned (this is contained in
!   LyrFields(1:NumLyrs)).  To each gridpoint in the fields comprising
!   the second half of the output the maximum of the corresponding
!   data is assigned (this is contained in
!   LyrFields(NumLyrs+1:2*NumLyrs)).
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
USE FldCodes_mod, ONLY:   &
  ST_MnIceP, ST_MxIceP,   &
  MO8_MnIceP, MO8_MxIceP, &
  PP_MnIceP, PP_MxIceP,   &
  VC_MnIceP, VC_MxIceP
USE ereport_mod, ONLY : ereport
IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs
INTEGER, INTENT(IN) :: NumLyrs
INTEGER, INTENT(IN) :: LyrMO8L(NumLyrs)
REAL,    INTENT(IN) :: LyrLwrB(NumLyrs)
REAL,    INTENT(IN) :: LyrUprB(NumLyrs)
TYPE(PP_Field_type), INTENT(IN)    :: LevFields(NumLevs)
TYPE(PP_Field_type), INTENT(IN)    :: PFields(NumLevs)
TYPE(PP_Field_type), INTENT(INOUT) :: LyrFields(2*NumLyrs)
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "ICOnPLyr"

! Local Variables:
INTEGER :: i, j, k, n, l
REAL    :: RData(NumLevs), Press(NumLevs)
LOGICAL :: NotMDI(NumLevs), Mask(NumLevs)

! new local variables for vectorization
integer :: n_count(LyrFields(1) % Hdr % NumCols, &
                   LyrFields(1) % Hdr % NumRows,NumLyrs)
real :: max_value(LyrFields(1) % Hdr % NumCols, &
                  LyrFields(1) % Hdr % NumRows,NumLyrs)
real :: sum_value(LyrFields(1) % Hdr % NumCols, &
                  LyrFields(1) % Hdr % NumRows,NumLyrs)
real :: Rdata_3(LyrFields(1) % Hdr % NumCols, &
                LyrFields(1) % Hdr % NumRows,NumLevs)
! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF (NumLevs < NumLyrs) THEN
  ErrorStatus = StatusWarning

  CALL EReport( RoutineName, ErrorStatus,                             &
                "Too few model level fields supplied" )
  ErrorStatus = StatusWarning
  GO TO 9999
END IF

DO k = 1, NumLevs

  IF ( (LevFields(k) % Hdr % NumCols   /=                             &
        PFields(k) % Hdr % NumCols) .OR.                              &
       (LevFields(k) % Hdr % NumRows   /=                             &
        PFields(k) % Hdr % NumRows) .OR.                              &
       (LevFields(k) % Hdr % ZerothLon /=                             &
        PFields(k) % Hdr % ZerothLon) .OR.                            &
       (LevFields(k) % Hdr % ZerothLat /=                             &
        PFields(k) % Hdr % ZerothLat) ) THEN
    ErrorStatus = StatusWarning

    CALL EReport( RoutineName, ErrorStatus,                           &
                  "Model level fields supplied are on different grids" )
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF

END DO

DO k = 1, 2*NumLyrs

  IF ( ASSOCIATED( LyrFields(k) % RData ) ) THEN
    DEALLOCATE( LyrFields(k) % RData )
    NULLIFY( LyrFields(k) % RData )
  END IF

  IF ( .NOT. ASSOCIATED( LyrFields(k) % RData ) ) THEN
    ALLOCATE( LyrFields(k) % RData(LevFields(1) % Hdr % NumCols,     &
                                   LevFields(1) % Hdr % NumRows) )
  END IF

END DO


!Set the headers for the output fields held in LyrFields. Do this by
! using the values from the headers of the input fields (in LevFields)
! but re-set the fieldcodes as approporiate. Also set the parameters
! in the headers to specify the bounderies of the layers of atmosphere.

LyrFields(1:NumLyrs) % Hdr           = LevFields(1:NumLyrs) % Hdr
LyrFields(NumLyrs+1:2*NumLyrs) % Hdr = LevFields(1:NumLyrs) % Hdr
LyrFields(1:NumLyrs) % Hdr % STCode  =  ST_MnIceP
LyrFields(NumLyrs+1:2*NumLyrs) % Hdr % STCode  =  ST_MxIceP

LyrFields(1:NumLyrs) % Hdr % MO8Type = MO8_MnIceP
LyrFields(NumLyrs+1:2*NumLyrs) % Hdr % MO8Type = MO8_MxIceP

LyrFields(1:NumLyrs) % Hdr % PPCode  =  PP_MnIceP
LyrFields(NumLyrs+1:2*NumLyrs) % Hdr % PPCode  =  PP_MxIceP

LyrFields(1:NumLyrs) % Hdr % LBVC  =  VC_MnIceP
LyrFields(NumLyrs+1:2*NumLyrs) % Hdr % LBVC  =  VC_MxIceP

LyrFields(1:NumLyrs) % Hdr % MO8Level           = LyrMO8L
LyrFields(NumLyrs+1:2*NumLyrs) % Hdr % MO8Level = LyrMO8L

LyrFields(1:NumLyrs) % Hdr % RLevel             = LyrUprB
LyrFields(NumLyrs+1:2*NumLyrs) % Hdr % RLevel   = LyrUprB

LyrFields(1:NumLyrs) % Hdr % RefLevel           = LyrLwrB
LyrFields(NumLyrs+1:2*NumLyrs) % Hdr % RefLevel = LyrLwrB

LyrFields % Hdr % LBVC   = 8
LyrFields % Hdr % LBRVC  = 8

LyrFields % Hdr % BULev  = 0.0
LyrFields % Hdr % BHULev = 0.0
LyrFields % Hdr % BHLev  = 0.0
LyrFields % Hdr % BHRLev = 0.0

n_count=0
max_value=0.0
sum_value=0.0

DO k=1,NumLevs
   DO j = 1, LyrFields(1) % Hdr % NumRows
      DO i = 1, LyrFields(1) % Hdr % NumCols
         Rdata_3(i,j,k)= LevFields(k) % RData(i,j)
      END DO
   END DO
END DO


!Within the loops below:
!1) the number of model theta levels falling within each desired atmosphere
!   layer (for each gridpoint) are counted;
!2) the values of the input field occuring in this layer are summed;
!3) and the maximum value of the input field in the layer is found.

!First loop over gridpoints
DO j = 1, LyrFields(1) % Hdr % NumRows
  DO i = 1, LyrFields(1) % Hdr % NumCols

  !Now loop over desired atmosphere layers
    DO l=1, NumLyrs

      !For this gridpoint and atmosphere layer, loop through the
      !model levels and see if the pressure on these levels
      !occur within the atmosphere layer
      DO k=1,NumLevs

        IF ( ( RData_3(i,j,k) /=  LevFields(1) % Hdr % BMDI ).AND. &
          (  PFields(k) % RData(i,j) /= PFields(1) % Hdr % BMDI )  &
          .AND.   &
          (  0.01*PFields(k) % RData(i,j) <= LyrLwrB(l)) .AND.     &
          (  0.01*PFields(k) % RData(i,j) >= LyrUprB(l) ))  THEN

          n_count(i,j,l)=n_count(i,j,l)+1
          IF (max_value(i,j,l).le.RData_3(i,j,k)) THEN
            max_value(i,j,l)=RData_3(i,j,k)
          END IF
          sum_value(i,j,l)=sum_value(i,j,l)+ RData_3(i,j,k)
        END IF
      END DO
    END DO
  END DO
END DO



!For each gridpoint and desired layer:
!  1) calculate the mean input field value and put this into the first
!     set of fields contained in LyrFields
!  2) copy the maximum value of the input field in this layer into the
!     second set of fields contained in LyrFields
!  3) if there are no input field values in the desired layer, set the
!     output mean and maximum field values to be equal to the missing data
!     indicator

!Loop over gridpoints
DO j = 1, LyrFields(1) % Hdr % NumRows
  DO i = 1, LyrFields(1) % Hdr % NumCols

    !Loop over desired atmosphere layers
    DO l=1, NumLyrs

      IF (n_count(i,j,l) > 0) then
          LyrFields(l) % RData(i,j) = sum_value(i,j,l)
          LyrFields(l) % RData(i,j) = LyrFields(l) % RData(i,j) &
                                        /real(n_count(i,j,l))

          LyrFields(l+NumLyrs) % RData(i,j) = max_value(i,j,l)
      ELSE
          LyrFields(l) % RData(i,j) = LyrFields(l) % Hdr % BMDI
          LyrFields(l+NumLyrs) % RData(i,j) =                   &
                               LyrFields(l+NumLyrs) % Hdr % BMDI
      END IF

    END DO

  END DO
END DO


9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE ICOnPLyr
