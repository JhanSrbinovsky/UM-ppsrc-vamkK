! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to calculate Clear Air Turbulence diagnostics.


SUBROUTINE CATurb( NumLevs,      &  ! in
                   UStdLev,      &  ! in
                   VStdLev,      &  ! in
                   UFields,      &  ! in
                   VFields,      &  ! in
                   PFields,      &  ! in
                   ZFields,      &  ! in
                   CATPred,      &  ! inout
                   ErrorStatus )    ! inout

! Description:
!   Routine to calculate Dutton's CAT indicator and convert to a CAT
!   predictor with a value between 0 and 7.5 for display on Av SigWX
!   charts.
!
! Method:
!   Calculate horizontal and vertical derivatives of wind components and
!   combine these using Dutton's algorithm.
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
USE FldCodes_mod, ONLY:   &
  ST_CAT, MO8_CAT, PP_CAT
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumLevs                      ! No. levels
TYPE(PP_Field_type), INTENT(IN) :: UStdLev          ! U & V wind comps
TYPE(PP_Field_type), INTENT(IN) :: VStdLev          !   on std level
TYPE(PP_Field_type), INTENT(IN) :: UFields(NumLevs) ! U & V wind comp,
TYPE(PP_Field_type), INTENT(IN) :: VFields(NumLevs) !   pressure and
TYPE(PP_Field_type), INTENT(IN) :: PFields(NumLevs) !   height on rho
TYPE(PP_Field_type), INTENT(IN) :: ZFields(NumLevs) !   levels

TYPE(PP_Field_type), INTENT(INOUT) :: CATPred       ! CAT diagnostic
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "CATurb"
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
REAL :: PRef
TYPE(PP_Field_type) :: dUdX, dUdY, dVdX, dVdY      ! Horiz. Derivatives
TYPE(PP_Field_type) :: dZdP, dUdP, dVdP            ! Vert. Deriv (wrt p)
TYPE(PP_Field_type) :: dWdZ                        ! Vert. Deriv (wrt z)

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

Nullify ( dZdP % RData )
Nullify ( dUdP % RData )
Nullify ( dVdP % RData )
Nullify ( dUdX % RData )
Nullify ( dVdX % RData )
Nullify ( dUdY % RData )
Nullify ( dVdY % RData )

IF ( ASSOCIATED( CATPred % RData ) ) THEN
  DEALLOCATE( CATPred % RData )
END IF
! Set up header & allocate memory
CATPred % Hdr = UStdLev % Hdr
CATPred % Hdr % PPCode  =  PP_CAT
CATPred % Hdr % MO8Type = MO8_CAT
CATPred % Hdr % STCode  =  ST_CAT
CATPred % Hdr % BMDI    = RMDI
ALLOCATE( CATPred % RData(CATPred % Hdr % NumCols, &
                          CATPred % Hdr % NumRows) )

!---- VERTICAL DERIVATIVES ----
! Find vertical derivatives of z, u & v wrt pressure using cubic spline
! Note that spline values must be monotonically increasing.
PRef = UStdLev % Hdr % RLevel * 100.0
! DEPENDS ON: diffp
CALL DiffP( NumLevs, PRef, ZFields, PFields, dZdP, ErrorStatus )
! DEPENDS ON: diffp
CALL DiffP( NumLevs, PRef, UFields, PFields, dUdP, ErrorStatus )
! DEPENDS ON: diffp
CALL DiffP( NumLevs, PRef, VFields, PFields, dVdP, ErrorStatus )

dWdZ % Hdr = UFields(1) % Hdr
ALLOCATE( dWdZ % RData(dWdZ % Hdr % NumCols, &
                       dWdZ % Hdr % NumRows) )
dWdZ % RData(:,:) = SQRT(dUdP % RData(:,:)**2 + dVdP % RData(:,:)**2) &
                                                / ABS(dZdP % RData(:,:))

!---- HORIZONTAL DERIVATIVES ----
! latitudinal derivatives
! DEPENDS ON: diffx
CALL DiffX( UStdLev, dUdX, ErrorStatus )
! DEPENDS ON: diffx
CALL DiffX( VStdLev, dVdX, ErrorStatus )
! longitudinal derivatives
! DEPENDS ON: diffy
CALL DiffY( UStdLev, dUdY, ErrorStatus )
! DEPENDS ON: diffy
CALL DiffY( VStdLev, dVdY, ErrorStatus )

!---- Calculate CAT indicator ----
! DEPENDS ON: duttoncat
CALL DuttonCAT( UStdLev, VStdLev,             &
                dUdX, dUdY, dVdX, dVdY, dWdZ, &
                CATPred, ErrorStatus )

DEALLOCATE( dUdP % RData )
DEALLOCATE( dVdP % RData )
DEALLOCATE( dZdP % RData )
DEALLOCATE( dUdX % RData )
DEALLOCATE( dVdX % RData )
DEALLOCATE( dUdY % RData )
DEALLOCATE( dVdY % RData )
DEALLOCATE( dWdZ % RData )

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE CATurb

!=======================================================================


!=======================================================================


!=======================================================================


