! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to calculate surface dust conc and ~2000-5000ft conc

SUBROUTINE DUST_CON (Numlevs,       & !in
                     Orog,          & !in
                     Zfields,       & !in
                     dust,          & !in
                     dustc,         & !inout
                     dusty,         & !inout
                     ErrorStatus )    !inout

! Description: Calculates dust concentrations.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6
!

USE atmos_constants_mod, ONLY: r

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK
USE FldCodes_Mod, ONLY:             &
  ST_Ptheta, ST_Ttheta,             &
  ST_Dust,                          &
  ST_DUSTC,  MO8_DUSTC,  PP_DUSTC,  &
  ST_DUSTY,  MO8_DUSTY,  PP_DUSTY,  &
  VC_Surface,                       &
  LV_Surface,  LV_Special

IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs
TYPE(PP_Field_type), INTENT(IN) :: Orog              ! Model Orography
TYPE(PP_Field_type), INTENT(IN) :: ZFields(NumLevs)  ! height on rho levels
TYPE(PP_Field_type), INTENT(IN) :: dust(NumLevs)     ! and dust conc on theta
                                                     ! (microg/m3)

TYPE(PP_Field_type), INTENT(INOUT) :: dustc  ! surface dust concentration (g/m3)
TYPE(PP_Field_type), INTENT(INOUT) :: dusty  ! 2000-5000ft dust concentration (g/m3)

INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "DUST_CON"

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

! Layer dimensions in metres.
REAL, PARAMETER :: top_of_layer = 1524.0
REAL, PARAMETER :: bottom_of_layer = 609.0
REAL, PARAMETER :: layer_thickness = 915.0

! Local Variables:
INTEGER :: i,j,k           ! Loop counters

!  thickness of dust layer
REAL :: thick (dust(1)%Hdr%NumCols, dust(1)%Hdr%NumRows)

! End of header -----------------------------------------------

CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( dustc % RData ) ) THEN
  DEALLOCATE( dustc % RData )
END IF
IF ( ASSOCIATED( dusty % RData ) ) THEN
  DEALLOCATE( dusty % RData )
END IF

dustc % Hdr = dust(1) % Hdr
dustc % Hdr % LBVC     = VC_Surface
dustc % Hdr % MO8Level = LV_Surface
dustc % Hdr % BULEV    = 0.0
dustc % Hdr % BHULEV   = 0.0
dustc % Hdr % RLevel   = 0.0
dustc % Hdr % RefLevel = 0.0
dustc % Hdr % BHLEV    = 0.0
dustc % Hdr % BHRLEV   = 0.0
dustc % Hdr % BMDI     = RMDI
dustc % Hdr % PPCode   =  PP_DUSTC
dustc % Hdr % MO8Type  = MO8_DUSTC
dustc % Hdr % STCode   =  ST_DUSTC

dusty % Hdr = dustc % Hdr
dusty % Hdr % LBPROC   = 2048  !  weighted mean between two levels
dusty % Hdr % LBVC     = VC_Surface
dusty % Hdr % MO8Level = LV_Special
dusty % Hdr % RLevel   = 5000.0 ! 5000ft (1524m)
dusty % Hdr % RefLevel = 2000.0 ! 2000ft (609m)
dusty % Hdr % PPCode   =  PP_DUSTY
dusty % Hdr % MO8Type  = MO8_DUSTY
dusty % Hdr % STCode   =  ST_DUSTY

ALLOCATE( dustc % RData(dustc % Hdr % NumCols, &
                        dustc % Hdr % NumRows) )
ALLOCATE( dusty % RData(dusty % Hdr % NumCols, &
                        dusty % Hdr % NumRows) )

! loop over model levels
! here we consider conc between 2000 and 5000ft agl
! so approx 609 and 1524m and mask out data if outside this layer

dusty%RData(:,:) = 0.0

DO k = 2,Numlevs-1

  thick=( ZFields(k+1)%RData - ZFields(k)%RData )

  WHERE ( (ZFields(k+1)%RData - Orog%RData) <  bottom_of_layer )
        thick(:,:)=0.0
  END WHERE

  WHERE ( ZFields(k)%RData - Orog%RData >  top_of_layer )
       thick(:,:)=0.0
  END WHERE

  DO j = 1,dustc % Hdr % NumRows
    DO i = 1,dustc % Hdr % NumCols

       dusty%RData(i,j) = dusty%RData(i,j) +    &
                         ( thick(i,j) * dust(k)%RData(i,j) )

     END DO
  END DO

END DO

DO j = 1,dustc % Hdr % NumRows
  DO i = 1,dustc % Hdr % NumCols
    ! convert from microg/m3 to g/m3

    ! for final conc scale by ~thickness layer
    dusty%RData(i,j) = 1.0E-6 * dusty%RData(i,j) / layer_thickness

    ! surface conc much simpler here we take level 1
    dustc%RData(i,j) = 1.0E-6 * dust(1)%RData(i,j)

  END DO
END DO


9999 CONTINUE

CALL Timer( RoutineName, 4 )

END SUBROUTINE DUST_CON
