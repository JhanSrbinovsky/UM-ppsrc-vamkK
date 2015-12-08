! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
SUBROUTINE Fill_Pressure( Input_fld,    &  !inout
                          Input_mask,   &  !in
                          Output_fld,   &  !inout
                          ErrorStatus)     !inout

! Description:
!   Interpolate field values within Cb mask area.
!
! Method:
! Performs linear interpolation horizontally E-W and vertically N-S
! to fill in the field values (pressure or heights of Cb bases and tops
! or horizontal extent of Cb) at points within the Cb mask that have
! missing values.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6
!
!
!
USE IO_mod, ONLY:        &
    PP_Header_type,      &
    PP_Field_type
USE Err_Mod, ONLY:       &
    StatusOK

IMPLICIT None


TYPE(PP_Field_type), INTENT(INOUT) :: Input_fld  ! Input/output field index
TYPE(PP_Field_type), INTENT(IN)    :: Input_mask ! Mask field index
TYPE(PP_Field_type), INTENT(INOUT) :: Output_fld ! Dummy field index (can
                                                 ! also be used for output)
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Fill_Pressure"
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

INTEGER :: i, j, k                                  ! loop counters
REAL    :: a                                        ! holds adjusted field value

!  Variables p & q contain the input field values in a 3x3 box around the point
!  of interest (p5/q5)
!
!  p1,p2,p3,
!  p4,p5,p6,
!  p7,p8,p9,
!
!  q1,q2,q3,
!  q4,q5,q6,
!  q7,q8,q9
!
! Only p5, q3, q5, q6, q9 are used in this routine
!
REAL    :: p5, q3, q5, q6, q9


! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF



! Copy input field to dummy field ready for checking the field and making
! adjustments.
Output_fld % RData(:,:)=Input_fld % RData(:,:)


!---- Check through the field and adjust according to neighbouring points ----
!
! If input field has a missing value at point of interest but mask field
! indicates Cb, check if there is an input field value at the southerly,
! northery or easterly neighbouring gridpoints. If so, set value at point
! of interest to be this value. Repeat check through field a further 9 times.
!
! Optimisation comments:
! kk Exchanging I and J loop avoids bank conflicts AND allows vectorization.

DO k = 1, 10
  DO j = 1, Input_mask % Hdr % NumRows
    DO i = 1, Input_mask % Hdr % NumCols

      p5 = Input_mask % RData(i,j)
      q3 = Output_fld % RData(i+1,j+1)
      q5 = Output_fld % RData(i,j)
      q6 = Output_fld % RData(i+1,j)
      q9 = Output_fld % RData(i+1,j-1)

      IF ( (q5 == RMDI ) .AND. ( p5 /= RMDI)) THEN
        a=RMDI
        IF (q9 /= RMDI )  THEN
          a=q9
        END IF
        IF (q3 /= RMDI ) THEN
          a=q3
        END IF
        IF (q6 /= RMDI ) THEN
          a=q6
        END IF
        Output_fld % RData(i,j)=a
      END IF
    END DO
  END DO
END DO


! Copy adjusted field back to Pfields(Inf1) ready for output
Input_fld % Rdata(:,:) = Output_fld % RData(:,:)


9999 CONTINUE

END SUBROUTINE Fill_Pressure
