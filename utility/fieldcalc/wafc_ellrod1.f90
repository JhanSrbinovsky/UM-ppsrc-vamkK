! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE WAFC_ELLROD1( UStdLev,     &  ! in
                         VStdLev,     &  ! in
                         dUdX,        &  ! in
                         dUdY,        &  ! in
                         dVdX,        &  ! in
                         dVdY,        &  ! in
                         dWdZ,        &  ! in
                         CATPred,     &  ! inout
                         ErrorStatus )   ! inout

! Description:
!   Routine to calculate Ellrods 1st CAT index TI1.
!
! Method:
!   1) Calculate deformation and vertical wind shear.
!   2) Use Ellrod's empirical formula to calculate the index:
!        TI1 = vert. wind shear * Deformation
!        (see Ellrod and Knapp, 1991)
!   3) Convert TI1 values to percentage probability using  equation
!         prob_TI1 = (1383149.3*TI1) + 2.7258455
!        which was obtained by correlating TI1 with CAT probability.
!   4) Lastly, set all values < 3.5%  to zero.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK
IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: UStdLev       ! U & V wind comps on
TYPE(PP_Field_type), INTENT(IN) :: VStdLev       !   std level
TYPE(PP_Field_type), INTENT(IN) :: dUdX, dUdY    ! Horiz. derivs of u
TYPE(PP_Field_type), INTENT(IN) :: dVdX, dVdY    ! Horiz. derivs of v
TYPE(PP_Field_type), INTENT(IN) :: dWdZ          ! Vert. deriv of wspeed

TYPE(PP_Field_type), INTENT(INOUT) :: CATPred     ! CAT diagnostic
!TYPE(PP_Field_type), INTENT(INOUT) :: CATPred_dup ! CAT diagnostic duplica
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "WAFC_Ellrod1"

! Local Variables:
INTEGER :: i, j                  ! Loop counters
REAL    :: alat                  ! Latitude
REAL    :: ShearV                ! Vertical wind shear   (VWS)
REAL    :: St_Def                ! Stretching deformation
REAL    :: Sh_Def                ! Shearing deformation
REAL    :: Def                   ! Total deformation
REAL    :: TI1                   ! Ellrod turbulence indicator TI1
REAL    :: prob_TI1              ! TI1 converted to percentage probability
                                 !   of encountering turbulence

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF


!------------------------------------------------------------------------
! Calculate shearing & stretching deformation and TI1 for each point
! Then convert TI1 to percentage probability.

DO j = 1, CATPred % Hdr % NumRows

  ! calculate latitude
  alat = UStdLev % Hdr % ZerothLat + &
         UStdLev % Hdr % LatInt * FLOAT(j)

  DO i = 1,CATPred % Hdr % NumCols

    IF ( (dUdX % RData(i,j) == dUdX % Hdr % BMDI) .OR.  &
         (dUdY % RData(i,j) == dUdY % Hdr % BMDI) .OR.  &
         (dVdX % RData(i,j) == dVdX % Hdr % BMDI) .OR.  &
         (dVdY % RData(i,j) == dVdY % Hdr % BMDI) ) THEN
      CATPred % RData(i,j) = CATPred % Hdr % BMDI
    ELSE

      ShearV = dWdZ % RData(i,j)

      !Calculate shearing & stretching deformation, then TI1
      St_Def = dUdX % RData(i,j) - dVdY % RData(i,j)
      Sh_Def = dVdX % RData(i,j) + dUdY % RData(i,j)
      Def = SQRT( (St_Def**2) + (Sh_Def**2) )
      TI1 = ShearV * Def

      !Convert to percentage probability
      prob_TI1 = (1383149.3*TI1) + 2.7258455
      IF (prob_TI1 < 3.5) THEN
        prob_TI1=0.0
      END IF

      CATPred % RData(i,j) = prob_TI1

    END IF
  END DO
END DO

9999 CONTINUE

END SUBROUTINE WAFC_ELLROD1
