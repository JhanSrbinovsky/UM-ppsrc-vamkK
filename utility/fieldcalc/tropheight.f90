! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate tropopause height

SUBROUTINE TropHeight( NumLevs,      &  ! in
                       PFields,      &  ! in
                       TFields,      &  ! in
                       ZFields,      &  ! in
                       TropSTCode,   &  ! in
                       TropMax,     &  ! in
                       TropP,        &  ! inout
                       TropT,        &  ! inout
                       TropZ,        &  ! inout
                       ErrorStatus )    ! inout

! Description:
! Find interval (k-1 to k) in which tropopause lies by checking the lapse
! rate of surrounding intervals, and then calculate the corresponding 
! height, pressure and temperature at the tropopause.
!
! Method:
! 1) Set up PP-headers and field codes for output fields.
! 2) Calculate lapse rate field for model level intervals 1~2 and 2~3 
!    (1 here is the lowest level passed to the subroutine, not 
!     necessarily the lowest in the model)
! 3) Loop over k from 3 to number of levels (NumLevs-1)
!   (NB. operational version loops to NumLevs; this version to NumLevs-1
!    to prevent k=NumLevs which would cause problems later)
!    i)  Calculate lapse rate field for interval k~k+1
!    ii) For each horizontal grid point, WHERE:
!          Pressure in interval k-1~k is between 50000 and 5000Pa, and
!          Temperature is less than -30 deg C (243.0 K), and
!          Lapse rate of current interval (k-1~k) and interval above (k~k+1) 
!           is less than tropopause lapse rate (lapse_trop), and
!          Lapse rate of interval below (k-2~k-1) is > 0, and
!          Tropopause not previously found (at a level below), 
!        then set tropopause interval level to k (i.e. k-1~k) 
! 4) Loop through each horizontal grid point
!    i) If tropopause level not found, set to level of previous point used 
!       (if first point, set to NumLevs-2)
!    ii)Find tropopause height by finding point at which the temperature 
!       profiles (assumed to be straight lines) below and above tropopause 
!       intersect.
!       Equation of first line: 
!         in form y=mx+c where m=-LapseLwr and c = TLwr + LapseLwr*ZLwr
!         T = - LapseLwr*Z + (TLwr + LapseLwr*ZLwr) 
!       Equation of second line:
!         in form y=mx+c where m=-LapseUpr and c = TUpr + LapseUpr*ZUpr
!         T = - LapseUpr*Z + (TUpr + LapseUpr*ZUpr)
!       Point of intersection is given by:
!         - LapseUpr*Z + (TUpr + LapseUpr*ZUpr) = 
!                             - LapseLwr*Z + (TLwr + LapseLwr*ZLwr)
!       => (LapseLwr - LapseUpr)*Z = 
!                             (TLwr + LapseLwr*ZLwr) - (TUpr + LapseUpr*ZUpr)
!       =>  Z = (TLwr + LapseLwr*ZLwr) - 
!                             (TUpr + LapseUpr*ZUpr)/(LapseLwr - LapseUpr)
!
!     iii) Calculate tropopause temperature by substuting trop height into
!          equation of first line:
!          TropT = TLwr - LapseLwr * (TropZ - ZLwr)
!            
!     iv) Calculate tropopause pressure:
!          TropP = PLwr * (TropT/TLwr)**(G_over_R/LapseLwr)
!
! 5) Set arbitrary max tropopause as follows:
!      if pressure < 10100Pa, set pressure to 10100Pa, 
!      temperature to 199K, height to 16180m 
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE atmos_constants_mod, ONLY: r, lapse_trop

USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type
USE Err_Mod, ONLY:        &
  StatusOK
USE FldCodes_Mod, ONLY:         &
  ST_TropP,  MO8_TropP,  PP_P,  &
  ST_TropT,  MO8_TropT,  PP_T,  &
  ST_TropZ,  MO8_TropZ,  PP_Z,  &
  ST_TropP_GRIB2,  MO8_TropP_GRIB2,  &
  ST_TropT_GRIB2,  MO8_TropT_GRIB2,  &
  ST_TropZ_GRIB2,  MO8_TropZ_GRIB2,  &

  VC_Trop,   LV_Special


USE earth_constants_mod, ONLY: g

IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs
TYPE(PP_Field_type), INTENT(IN) :: PFields(NumLevs) ! Pressure, temp
TYPE(PP_Field_type), INTENT(IN) :: TFields(NumLevs) !   and height
TYPE(PP_Field_type), INTENT(IN) :: ZFields(NumLevs) !   on theta levels

INTEGER, INTENT(IN) :: TropSTCode                   ! STASH for which trop.
REAL, INTENT(IN) :: TropMax(3)                      ! Height of trop
TYPE(PP_Field_type), INTENT(INOUT) :: TropP         ! Pressure, temp
TYPE(PP_Field_type), INTENT(INOUT) :: TropT         !   and height
TYPE(PP_Field_type), INTENT(INOUT) :: TropZ         !   at Tropopause
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "TropHeight"
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
REAL, PARAMETER :: G_over_R = G / R
REAL, PARAMETER :: VSmall = 1.e-6

! Local Variables:
INTEGER :: i,j,k           ! Loop counters
INTEGER :: TLev(PFields(1)%Hdr%NumCols,PFields(1)%Hdr%NumRows)
REAL :: PLwr, PUpr
REAL :: TLwr, TUpr
REAL :: ZLwr, ZUpr
REAL :: LapseUpr, LapseLwr
! Lapse rates below, at and above current layer
REAL :: Lapse_b(PFields(1)%Hdr%NumCols,PFields(1)%Hdr%NumRows)
REAL :: Lapse_ (PFields(1)%Hdr%NumCols,PFields(1)%Hdr%NumRows)
REAL :: Lapse_a(PFields(1)%Hdr%NumCols,PFields(1)%Hdr%NumRows)
REAL :: delta_lapse        ! Lapse_b-Lapse_a

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( TropP % RData ) ) THEN
  DEALLOCATE( TropP % RData )
END IF
IF ( ASSOCIATED( TropT % RData ) ) THEN
  DEALLOCATE( TropT % RData )
END IF
IF ( ASSOCIATED( TropZ % RData ) ) THEN
  DEALLOCATE( TropZ % RData )
END IF

! Setup TropP header (info doesnt depend on limit of Tropopause)
TropP % Hdr = PFields(1) % Hdr
TropP % Hdr % LBVC     = VC_Trop
TropP % Hdr % MO8Level = LV_Special
TropP % Hdr % BULEV    = 0.0
TropP % Hdr % BHULEV   = 0.0
TropP % Hdr % RLevel   = 0.0
TropP % Hdr % RefLevel = 0.0
TropP % Hdr % BHLEV    = 0.0
TropP % Hdr % BHRLEV   = 0.0
TropP % Hdr % BMDI     = RMDI

! Copy TropP header into TropT and TropZ headers
TropT % Hdr = TropP % Hdr
TropZ % Hdr = TropP % Hdr

! Not setup information dependent on limit.
IF(TropSTCode == ST_TropP) THEN
  TropP % Hdr % PPCode   =  PP_P
  TropP % Hdr % MO8Type  = MO8_TropP
  TropP % Hdr % STCode   =  ST_TropP
  TropT % Hdr % PPCode   =  PP_T
  TropT % Hdr % MO8Type  = MO8_TropT
  TropT % Hdr % STCode   =  ST_TropT
  TropZ % Hdr % PPCode   =  PP_Z
  TropZ % Hdr % MO8Type  = MO8_TropZ
  TropZ % Hdr % STCode   =  ST_TropZ
ELSE IF(TropSTCode == ST_TropP_GRIB2) THEN
  TropP % Hdr % PPCode   =  PP_P
  TropP % Hdr % MO8Type  = MO8_TropP_GRIB2
  TropP % Hdr % STCode   =  ST_TropP_GRIB2
  TropT % Hdr % PPCode   =  PP_T
  TropT % Hdr % MO8Type  = MO8_TropT_GRIB2
  TropT % Hdr % STCode   =  ST_TropT_GRIB2
  TropZ % Hdr % PPCode   =  PP_Z
  TropZ % Hdr % MO8Type  = MO8_TropZ_GRIB2
  TropZ % Hdr % STCode   =  ST_TropZ_GRIB2
ELSE
  TropP % Hdr % PPCode   =  IMDI
  TropP % Hdr % MO8Type  =  IMDI
  TropP % Hdr % STCode   =  IMDI
  TropT % Hdr % PPCode   =  IMDI
  TropT % Hdr % MO8Type  =  IMDI
  TropT % Hdr % STCode   =  IMDI
  TropZ % Hdr % PPCode   =  IMDI
  TropZ % Hdr % MO8Type  =  IMDI
  TropZ % Hdr % STCode   =  IMDI
END IF


ALLOCATE( TropP % RData(TropP % Hdr % NumCols, &
                        TropP % Hdr % NumRows) )
ALLOCATE( TropT % RData(TropT % Hdr % NumCols, &
                        TropT % Hdr % NumRows) )
ALLOCATE( TropZ % RData(TropZ % Hdr % NumCols, &
                        TropZ % Hdr % NumRows) )

TLev (:,:) = 0
Lapse_  = (TFields(1) % RData - TFields(2) % RData) / &
                 (ZFields(2) % RData - ZFields(1) % RData)
Lapse_a = (TFields(2) % RData - TFields(3) % RData) / &
                 (ZFields(3) % RData - ZFields(2) % RData)

DO k = 3,NumLevs
  Lapse_b = Lapse_
  Lapse_  = Lapse_a
  IF ( k < NumLevs ) THEN
    Lapse_a = (TFields(k) % RData - TFields(k+1) % RData) / &
                     (ZFields(k+1) % RData - ZFields(k) % RData)
  END IF

  DO j = 1, PFields(1) % Hdr % NumRows
    DO i = 1, PFields(1) % Hdr % NumCols
      IF (TLev(i,j)                == 0          .AND. &
          Lapse_(i,j)               < lapse_trop .AND. &
          Lapse_a(i,j)              < lapse_trop .AND. &
          Lapse_b(i,j)              > 0.0        .AND. &
          PFields(k)   % RData(i,j) < 50000.0    .AND. &
          PFields(k-1) % RData(i,j) > 5000.0     .AND. &
          TFields(k)   % RData(i,j) < 243.0 ) THEN
        TLev(i,j) = k
      END IF
    END DO ! i-loop
  END DO ! j-loop
END DO ! k-loop

! Optimisation introduced at vn6.2 results in differences in Global
! as well as NAE and is reversed pending further investigation
k = NumLevs - 2

! Optimisation added at vn6.2
! DO j = 1,TropP % Hdr % NumRows
!   DO i = 1,TropP % Hdr % NumCols

DO j = 1,TropP % Hdr % NumRows
  DO i = 1,TropP % Hdr % NumCols

! Optimisation added at vn6.2 - scalar variable must be initialised
! inside the loop to allow vectorisation
!    k = NumLevs - 2

    IF ( TLev(i,j) /= 0 ) THEN   ! This is a TEMPORARY fix to avoid
      k = TLev(i,j)              !  falling over when trop is not
    END IF                       !  found.  Will have to rewrite this.

    PLwr = PFields(k-1) % RData(i,j)
    PUpr = PFields(k  ) % RData(i,j)
    TLwr = TFields(k-1) % RData(i,j)
    TUpr = TFields(k  ) % RData(i,j)
    ZLwr = ZFields(k-1) % RData(i,j)
    ZUpr = ZFields(k  ) % RData(i,j)
    LapseUpr = (TUpr - TFields(k+1) % RData(i,j)) / &
                      (ZFields(k+1) % RData(i,j) - ZUpr)
    LapseLwr = (TFields(k-2) % RData(i,j) - TLwr) / &
                      (ZLwr - ZFields(k-2) % RData(i,j))

    ! Find Z at tropopause
    delta_lapse = LapseLwr - LapseUpr
    IF ( ABS(delta_lapse) < VSmall ) THEN
      IF ( delta_lapse >= 0 ) delta_lapse =  VSmall
      IF ( delta_lapse <  0 ) delta_lapse = -VSmall
    END IF
    TropZ % RData(i,j) = ( (TLwr+(LapseLwr*ZLwr)) -   &
                                (TUpr+(LapseUpr*ZUpr)) ) / delta_lapse

    IF ( TropZ % RData(i,j) < ZLwr ) THEN
      TropZ % RData(i,j) = ZLwr  ! ensure trop level doesn't undershoot
    END IF
    IF ( TropZ % RData(i,j) > ZUpr )  THEN
      TropZ % RData(i,j) = ZUpr  ! or overshoot
    END IF

    ! T at tropopause
    TropT % RData(i,j) = TLwr - LapseLwr * (TropZ % RData(i,j) - ZLwr)

    IF ( ABS(LapseLwr) < VSmall ) THEN
      IF ( LapseLwr >= 0 ) LapseLwr =  VSmall
      IF ( LapseLwr <  0 ) LapseLwr = -VSmall
    END IF

    ! P at tropopause
    TropP % RData(i,j) = PLwr *  &
         (TropT % RData(i,j)/TLwr)**(G_over_R/LapseLwr)

  END DO
END DO

WHERE (TropP % RData < TropMax(1))      ! set arbitrary max tropopause
  TropP % RData = TropMax(1)
  TropT % RData = TropMax(2)
  TropZ % RData = TropMax(3)
END WHERE

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE TropHeight

