! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE WAFC_CATurb( NumLevs,      &  ! in
                        UStdLev,      &  ! in
                        VStdLev,      &  ! in
                        UFields,      &  ! in
                        VFields,      &  ! in
                        PFields,      &  ! in
                        ZFields,      &  ! in
                        CATPred,      &  ! inout
                        ErrorStatus )    ! inout

! Description:
!   Routine to calculate Ellrod's first CAT indicator for display on Av
!   SigWXcharts.
!
! Method:
!   Calculate horizontal and vertical derivatives of wind components and
!   then call subroutine Ellrod1 to combine these using Ellrod's first
!   algorithm
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_mod, ONLY:           &
  PP_Header_type,           &
  PP_Field_type
USE Err_Mod, ONLY:          &
  StatusOK
USE FldCodes_mod, ONLY:     &
  ST_MnCATPt, MO8_MnCATPt,  &
  PP_MnCATPt, VC_MnCATPt,   &
  ST_MxCATPt, MO8_MxCATPt,  &
  PP_MxCATPt, VC_MxCATPt
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN)                :: NumLevs          ! No. levels
TYPE(PP_Field_type), INTENT(IN)    :: UStdLev          ! U & V wind comps
TYPE(PP_Field_type), INTENT(IN)    :: VStdLev          !   on std level
TYPE(PP_Field_type), INTENT(IN)    :: UFields(NumLevs) ! U & V wind comp,
TYPE(PP_Field_type), INTENT(IN)    :: VFields(NumLevs) !   pressure and
TYPE(PP_Field_type), INTENT(IN)    :: PFields(NumLevs) !   height on rho
TYPE(PP_Field_type), INTENT(IN)    :: ZFields(NumLevs) !   levels
TYPE(PP_Field_type), INTENT(INOUT) :: CATPred          ! CAT predictor
INTEGER, INTENT(INOUT)             :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "WAFC_CATurb"
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
REAL                :: PRef
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
Nullify ( dWdZ % RData )


!------------------------------------------------------------------------
!Set up/ allocate CAT predictor field
IF ( ASSOCIATED( CATPred % RData ) ) THEN
  DEALLOCATE( CATPred % RData )
END IF
! Set up header & allocate memory
CATPred % Hdr = UStdLev % Hdr
CATPred % Hdr % PPCode  =  PP_MnCATPt
CATPred % Hdr % MO8Type = MO8_MnCATPt
CATPred % Hdr % STCode  =  ST_MnCATPt
CATPred % Hdr % LBVC    =  VC_MnCATPt

CATPred % Hdr % BMDI    = RMDI
ALLOCATE( CATPred % RData(CATPred % Hdr % NumCols, &
                          CATPred % Hdr % NumRows) )


!-----------------------------------------------------------------------
! Calc. vertical derivatives of z, u & v wrt pressure using cubic spline
! Note that spline values must be monotonically increasing.
! Then calculate vertical wind shear dWdz
PRef = UStdLev % Hdr % RLevel * 100.0
! DEPENDS ON: diffp
CALL DiffP( NumLevs, PRef, ZFields, PFields, dZdP, ErrorStatus )
! DEPENDS ON: diffp
CALL DiffP( NumLevs, PRef, UFields, PFields, dUdP, ErrorStatus )
! DEPENDS ON: diffp
CALL DiffP( NumLevs, PRef, VFields, PFields, dVdP, ErrorStatus )

dWdZ % Hdr = UFields(1) % Hdr
ALLOCATE( dWdZ % RData(dWdZ % Hdr % NumCols, dWdZ % Hdr % NumRows) )

dWdZ % RData(:,:) = SQRT(dUdP % RData(:,:)**2 + dVdP % RData(:,:)**2) &
                        / ABS(dZdP % RData(:,:))


!-----------------------------------------------------------------------
! Calculate horizontal derivatives

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


!-----------------------------------------------------------------------
! Calculate CAT predictor (Ellrod's index TI1)
! DEPENDS ON: wafc_ellrod1
CALL WAFC_ELLROD1( UStdLev, VStdLev, &
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

END SUBROUTINE WAFC_CATurb
