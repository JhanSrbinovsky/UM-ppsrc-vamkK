! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE MtnStress_PRef( NumLevs,      &  ! in
                           Params,       &  ! in
                           PRef,         &  ! in
                           WindUPRef,    &  ! in
                           WindVPRef,    &  ! in
                           UStress,      &  ! inout
                           VStress,      &  ! inout
                           PFields,      &  ! inout
                           MWTPred,      &  ! inout
                           ErrorStatus )    ! inout

! Description/Method:
!
!This is a copy of the Mtn_Stress.f90 routine but altered so it calculates
!mountain wave predictor (held in MWTPred) for a certain pressure level (PR
!instead of determining the maximum mountain wave stress. It does this by:
!
!1) the pressure levels corresponding to the model levels which UStress and
!   VStress are on are passed into the subroutine (these are in hPa).
!
!2) the pressure level (in hPa) for which the mountain wave stress is to be
!   calculated is also passed in to the subroutine
!
!3) New fields UStress_PRef and VStress_PRef are set up to hold the
!   interpolated values of U and V gravity wave stress at the pressure
!   level of interest.
!
!4) For each gridpoint, the values of UStress and VStress at PRef are found
!   by linear interpolation.
!
!5) For each gridpoint, the magnitude of UStress and VStress is determined
!   and put into the field Stress.
!
!6) Then advect the stress/MWPred field the same way as the Max stress is
!   advected in the original version (this part hasn't been changed).
!
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_mod, ONLY:          &
  PP_Header_type,          &
  PP_Field_type
USE Err_Mod, ONLY:         &
  StatusOK
USE FldCodes_Mod, ONLY:    &
  ST_MnCATPt, MO8_MnCATPt, &
  PP_MnCATPt, VC_MnCATPt,  &
  LV_Special

IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs
REAL, INTENT(IN)    :: Params(4)     ! Params(1:3) change the best fit
                                     ! parameters
                                     ! Params(4) switches the contour code
                                     ! on (1) or off (0)
REAL, INTENT(IN)    :: PRef                   ! Pressure level (in Pa) where
                                              !  MWPred is to be calculated at
TYPE(PP_Field_type), INTENT(IN) :: WindUPRef  !  wind U-component at PRef
TYPE(PP_Field_type), INTENT(IN) :: WindVPRef  !  wind V-component at PRef
TYPE(PP_Field_type), INTENT(IN) :: UStress(NumLevs)
TYPE(PP_Field_type), INTENT(IN) :: VStress(NumLevs)
TYPE(PP_Field_type), INTENT(IN) :: PFields(NumLevs)  !pressure fields (Pa)

TYPE(PP_Field_type), INTENT(INOUT) :: MWTPred
INTEGER, INTENT(INOUT) :: ErrorStatus
!TYPE(PP_Field_type) :: MWTPred_prob

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "MtnStress_Ref"
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
INTEGER, PARAMETER :: NPts = 2                 ! No. interpolation coeffs
INTEGER, PARAMETER :: NumContours = 3
REAL,    PARAMETER :: Contours(NumContours) =  &
                          (/0.007,0.065,0.25/) ! The contour thresholds

! Local Variables:
INTEGER :: jlwr, jupr, jmid  ! Variables holding upper, lower & mid levels
INTEGER :: halfn             ! half no. of interpolation points
                             !   (use to set lowest and highest level)
INTEGER :: i, j, k           ! DO loop variables
INTEGER :: NumRows           ! The number of latitudes = header(18)
INTEGER :: NumCols           ! The number of longitudes = header(19)
INTEGER :: NumAbvCutOff      ! Number of points above Contours
REAL :: m_UStress            ! m in y=mx+c for UStress interpolation
REAL :: c_UStress            ! c in y=mx+c for UStress interpolation
REAL :: m_VStress            ! m in y=mx+c for VStress interpolation
REAL :: c_VStress            ! c in y=mx+c for VStress interpolation
REAL :: CutOff               ! One of the Contours values
REAL :: MinWSpeed            ! Wind threshold for advection
REAL :: prob                 ! Holds probability

INTEGER, ALLOCATABLE :: PtPosn(:,:)    ! Lat Lon array of
                                       ! stress > cutoff
INTEGER, ALLOCATABLE :: AdvPtPosn(:,:) ! advected LL array of
                                       ! stress > cutoff

TYPE(PP_Field_type) :: WindFPRef    ! wind speed at this press level
TYPE(PP_Field_type) :: WindDPRef    ! wind direction at this press level
TYPE(PP_Field_type) :: UStress_PRef ! UStress at PRef
TYPE(PP_Field_type) :: VStress_PRef ! VStress at PRef
TYPE(PP_Field_type) :: Stress       ! Magnitude of gravity wave stress at
                                    ! this pressure level

! New variables for vectorization

integer :: Lower(UStress(1) % Hdr % NumCols,   &
                   UStress(1) % Hdr % NumRows)
integer :: Upper(UStress(1) % Hdr % NumCols,   &
                   UStress(1) % Hdr % NumRows)

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( MWTPred % RData ) ) THEN
  DEALLOCATE( MWTPred % RData )
END IF
MWTPred % Hdr = WindUPRef % Hdr
MWTPred % Hdr % PPCode   =  PP_MnCATPt
MWTPred % Hdr % LBVC     =  VC_MnCATPt
MWTPred % Hdr % MO8Type  =  MO8_MnCATPt
MWTPred % Hdr % STCode   =  ST_MnCATPt
MWTPred % Hdr % RLevel   = 0.0
MWTPred % Hdr % RefLevel = 0.0
MWTPred % Hdr % BHLEV    = 0.0
MWTPred % Hdr % BHRLEV   = 0.0
MWTPred % Hdr % BULEV    = 0.0
MWTPred % Hdr % BHULEV   = 0.0
MWTPred % Hdr % BMDI     = RMDI
ALLOCATE( MWTPred % RData(MWTPred % Hdr % NumCols, &
                          MWTPred % Hdr % NumRows) )

! Initialise the local variables
NumCols = MWTPred % Hdr % NumCols
NumRows = MWTPred % Hdr % NumRows
NULLIFY( Stress % RData )
NULLIFY( UStress_PRef % RData )
NULLIFY( VStress_PRef % RData )
NULLIFY( WindFPRef % RData )
NULLIFY( WindDPRef % RData )



UStress_PRef % Hdr = UStress(1) % Hdr
ALLOCATE( UStress_PRef % RData(UStress_PRef % Hdr % NumCols, &
                               UStress_PRef % Hdr % NumRows) )
VStress_PRef % Hdr = VStress(1) % Hdr
ALLOCATE( VStress_PRef % RData(VStress_PRef % Hdr % NumCols, &
                               VStress_PRef % Hdr % NumRows) )


! DEPENDS ON: vecmag
CALL VecMag( WindUPRef, WindVPRef, WindFPRef, ErrorStatus )
! DEPENDS ON: vecdir
CALL VecDir( WindUPRef, WindVPRef, WindDPRef, ErrorStatus )


!--------------------------------------------------------------------------
! Determine model levels above and below this pressure level
! Then determine UStress & VStress at PRef by linear interpolation of
! UStress & VStress from the level above and below. Do this by first
! checking if Ustress is the same above and below PRef (in which case
! m=0) if so, set UStress equal to this value at PRef; similarly for
! VStress. Otherwise perform linear interpolation - first find the
! m and c coefficents for the equation y=mx+c for the Ustress and
! Vstress ready for interpolation.


halfn = NPts/2

DO j = 1,UStress(1) % Hdr % NumRows
   DO i = 1,UStress(1) % Hdr % NumCols
      IF (PRef >= PFields(halfn) % RData(i,j)) THEN !all smaller
         Lower(i,j) = halfn
         Upper(i,j) = halfn+1
      ELSE IF (PRef < PFields(NumLevs+1-halfn) % RData(i,j)) THEN
!all lower
         Lower(i,j) = NumLevs-halfn
         Upper(i,j) =  NumLevs+1-halfn
      ELSE
         Lower (i,j) = -1 ! One level to find
      END IF
   END DO
END DO

DO k=halfn , NumLevs+1-halfn
  DO j = 1,UStress(1) % Hdr % NumRows
     DO i = 1,UStress(1) % Hdr % NumCols

       IF ( ( PFields(k) % RData(i,j) > Pref  ) &
          .AND.( PFields(k+1) % RData(i,j) <= Pref  ) &
          .AND.(Lower(i,j).eq.-1) ) THEN
          Lower(i,j) = k
          Upper(i,j) = k+1
       END IF

     END DO
  END DO
END DO


DO j = 1,UStress(1) % Hdr % NumRows
  DO i = 1,UStress(1) % Hdr % NumCols


    m_UStress = 0.0
    m_VStress = 0.0
    c_UStress = 0.0
    c_UStress = 0.0

    jlwr=Lower(i,j)
    jupr=Upper(i,j)

    IF (UStress(jlwr)%RData(i,j) == UStress(jupr)%RData(i,j)) THEN
      UStress_PRef%RData(i,j) = UStress(jlwr)%RData(i,j)
    ELSE

      m_UStress =(UStress(jlwr)%RData(i,j)-UStress(jupr)%RData(i,j))  &
               / (PFields(jlwr)%RData(i,j)- PFields(jupr)%RData(i,j) )
      c_UStress = UStress(jlwr)%RData(i,j) - &
                  (m_UStress * PFields(jlwr)%RData(i,j))
      UStress_PRef%RData(i,j) = (m_UStress * PRef )+ c_UStress

    END IF


    IF (VStress(jlwr)%RData(i,j) == VStress(jupr)%RData(i,j)) THEN
      VStress_PRef%RData(i,j) = VStress(jlwr)%RData(i,j)
    ELSE
     m_VStress =(VStress(jlwr)%RData(i,j)-VStress(jupr) % RData(i,j)) &
                /(PFields(jlwr)%RData(i,j) - PFields(jupr)%RData(i,j) )
      c_VStress = VStress(jlwr)%RData(i,j) - &
                  (m_VStress * PFields(jlwr)%RData(i,j))
      VStress_PRef%RData(i,j) = (m_VStress * PRef )+ c_VStress
    END IF

  END DO
END DO



! Now have 2 fields that hold the interpolated values for Ustress and
! Vstress for this PRef for all grid points. Determine the magnitude at
! each gridpoint of the gravity wave stress.

! DEPENDS ON: vecmag
CALL VecMag( UStress_PRef, VStress_PRef, Stress, ErrorStatus )



! Set the MW predictor to the magnitude of the stress at each grid point.
! Then determine whether the stress or wind is strong enough to advect
! the stress and modify the MW predictor if needbe.


ALLOCATE( PtPosn(2,NumCols*NumRows) )

MWTPred % RData(:,:) = Stress % RData(:,:)

! ***************************************************************************
! * Until 70 level effects have been investigated this code can be switched *
! * off by setting params(4) to be 1                                        *
! ***************************************************************************
IF(INT(params(4))==0) THEN
  DO k = 1, NumContours     ! loop over the contour threshold values

    !  Determine the arrays of values above the thresholds.
    CutOff = Contours(k)
    NumAbvCutOff = 0       ! Initialise to 0
    PtPosn(:,:) = 0
    DO j = 1, NumRows
      DO i = 1, NumCols
        IF ( Stress % RData(i,j) > CutOff ) THEN
          ! Increase number of points by 1 and store in PtPosn
          NumAbvCutOff = NumAbvCutOff + 1
          PtPosn( 1:2, NumAbvCutOff ) = (/i,j/)
        END IF
      END DO
    END DO

    ! If NumAbvCutOff > 0, the 'new' mwt points need calculating.
    IF ( NumAbvCutOff > 0 ) THEN

      ALLOCATE( AdvPtPosn(2,NumAbvCutOff) )
      AdvPtPosn(:,:) = PtPosn(:,1:NumAbvCutOff)

      ! Determine the array of next points. Point 2 in Appendix.
      MinWSpeed = 0.0
  ! DEPENDS ON: advectmws_ext
      CALL AdvectMWS_Ext( NumAbvCutOff, CutOff,            &
                          MinWSpeed,                       &
                          PtPosn    (1:2,1:NumAbvCutOff),  &
                          Stress,                          &
                          WindFPRef,     WindDPRef,        &
                          AdvPtPosn (1:2,1:NumAbvCutOff),  &
                          MWTPred,                         &
                          ErrorStatus )

      !*** need to swap halos of AdvPtPosn if MPP***
      ! Determine the second array of next points. Point 4 in Appendix.
      MinWSpeed = 30.0  ! This is the wind speed required to advect a
                        ! distance of 2 grid squares.
  ! DEPENDS ON: advectmws_ext
      CALL AdvectMWS_Ext( NumAbvCutOff, CutOff,            &
                          MinWSpeed,                       &
                          PtPosn    (1:2,1:NumAbvCutOff),  &
                          Stress,                          &
                          WindFPRef,     WindDPRef,        &
                          AdvPtPosn (1:2,1:NumAbvCutOff),  &
                          MWTPred,                         &
                          ErrorStatus )

      ! If an extra point is required downwind for extra resolution
      ! or strong flow, calculate an extra point. Point 6 in Appendix.
      ! Possibly do this by calling AdvectMWS from within a loop.
      ! Calculate a MinWindSpd threshold from the grid resolution and pass
      ! it into the subroutine.

      DEALLOCATE(AdvPtPosn)
    END IF
  END DO
END IF

DO j = 1,MWTPred % Hdr % NumRows
  DO i = 1,MWTPred % Hdr % NumCols
    IF (MWTPred%RData(i,j) < Params(1))THEN
      prob = 0.0
    ELSE
      prob = Params(2)*MWTPred%RData(i,j) + Params(3)
    END IF
! By default these values are worked out from extrapolating the points
! 0.0645=>4.5% and 0.25=>15.5% thus when these are extrpoltaed
! we get an equation y=mx+c or y=42.667x+1.8333
! this extends up past 100% which needs to stop thus x canot be
! greater than 2.300764
! Currently need to scale by 0.1 since it was too strong otherwise.
    MWTPred%RData(i,j) = 0.1*MIN(prob,100.0)
  END DO
END DO


DEALLOCATE( WindFPRef % RData )
DEALLOCATE( WindDPRef % RData )
DEALLOCATE( PtPosn )
DEALLOCATE( UStress_PRef % RData )
DEALLOCATE( VStress_PRef % RData )
DEALLOCATE( Stress % RData )


9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE MtnStress_PRef
