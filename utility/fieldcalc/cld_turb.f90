! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
SUBROUTINE cld_turb( PRef,           &  ! in
                     NumLevs,        &  ! in
                     CldFields,      &  ! in
                     theta_e_Fields, &  ! in
                     PFields,        &  ! in
                     ZFields,        &  ! in
                     Cld_turb_Pred,  &  ! inout
                     ErrorStatus )      ! inout

! Description:
!    Calculates the in cloud turbulence predictor (Cld_turb_Pred).
!
! Method:
!   1) Set up (allocate) the in-cloud predictor field Cld_turb_Pred
!   2) Loop through the gridpoints and for each find the model level
!      immediately above and below the pressure level of interest.
!   3) Construct a cloud mask array indicating if there is cloud on the
!      model level above or below each grid point.
!   4) Find vertical derivatives of height and equivalent potential temp
!      w.r.t. pressure using cubic spline. Then calculate the verical
!      derviative of equivalent potential temperature with height.
!   5) Then, for each gridpoint, set the in cloud turbulence predictor
!      (cld_turb_Pred). If orography makes calculation impossible, set
!      cld_turb_Pred to missing data indicator. If calculation is possible,
!      check that the cloud mask=1. If so, check sign of dthetaedZ. If
!      negative mark predictor as turbulent (set to dthetadz); if positive
!      mark as non-turbulent (zero).
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
USE FldCodes_mod, ONLY:   &
  ST_MnICTP, MO8_MnICTP,  &
  PP_MnICTP, VC_MnICTP
USE Err_Mod, ONLY:        &
  StatusOK
IMPLICIT None

! Subroutine Arguments:
INTEGER, INTENT(IN) :: PRef                      ! input pressure in hPa
INTEGER, INTENT(IN) :: NumLevs                   ! number of model levels

TYPE(PP_Field_type), INTENT(IN)   :: CldFields(NumLevs)      !Cloud on
                                                             ! model levs
TYPE(PP_Field_type), INTENT(IN)   :: theta_e_Fields(NumLevs) !Equivalent
                                                             ! pot. temp
TYPE(PP_Field_type), INTENT(IN)   :: PFields(NumLevs)        !Pressure
TYPE(PP_Field_type), INTENT(IN)   :: ZFields(NumLevs)        !Height
TYPE(PP_Field_type), INTENT(INOUT):: Cld_turb_Pred           !Cloud turb.
                                                             ! predictor
INTEGER, INTENT(INOUT)            :: ErrorStatus


! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "cld_turb"
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
INTEGER              :: i, j, k           ! Loop counters
INTEGER              :: jlwr, jupr, jmid  ! Variables used to find model
                                          !  levels above and below level
                                          !  of interest
INTEGER, ALLOCATABLE :: Lower (:,:)       ! Lower model level
INTEGER, ALLOCATABLE :: Upper (:,:)       ! Upper model level
INTEGER, ALLOCATABLE :: MD_flag (:,:)     ! Missing data flag (used when
                                          !   level is beneath first model
                                          !   level)
INTEGER, ALLOCATABLE :: cld_mask (:,:)    ! Cloud mask
REAL                 :: PRef_Pa           ! Pressure reference level (Pa)
TYPE(PP_Field_type)  :: dZdP              ! Vert. Deriv of height wrt p
TYPE(PP_Field_type)  :: dthetaedP         ! Vert. Deriv of equivalent pot.
                                          !   temp wrt p
TYPE(PP_Field_type)  :: dthetaedZ         ! Vert. Deriv of equivalent pot.
                                          !   temp wrt p


! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

Nullify ( dZdP % RData )
Nullify ( dthetaedP % RData )
Nullify ( dthetaedZ % RData )


PRef_Pa = REAL(PRef) * 100.0

!--------------------------------------------------------------------------
!Set up output field(s)
IF ( ASSOCIATED( cld_turb_Pred % RData ) ) THEN
  DEALLOCATE( cld_turb_Pred % RData )
END IF

! Set up header & allocate memory
cld_turb_Pred % Hdr = CldFields(1) % Hdr
cld_turb_Pred % Hdr % STCode   =  ST_MnICTP
cld_turb_Pred % Hdr % MO8Type  = MO8_MnICTP
cld_turb_Pred % Hdr % PPCode   =  PP_MnICTP
cld_turb_Pred % Hdr % LBVC     =  VC_MnICTP
cld_turb_Pred % Hdr % MO8Level = PRef
cld_turb_Pred % Hdr % RLevel   = REAL(PRef)
cld_turb_Pred % Hdr % RefLevel = REAL(PRef)
cld_turb_Pred % Hdr % BULEV    = 0.0
cld_turb_Pred % Hdr % BHULEV   = 0.0
cld_turb_Pred % Hdr % BHLEV    = 0.0
cld_turb_Pred % Hdr % BHRLEV   = 0.0
cld_turb_Pred % Hdr % BMDI     = RMDI

ALLOCATE( cld_turb_Pred % RData(cld_turb_Pred % Hdr % NumCols, &
                               cld_turb_Pred % Hdr % NumRows) )


!------------------------------------------------------
!Find the model level immediately above and below each gridpoint


ALLOCATE ( Lower ( PFields(1) % Hdr % NumCols,     &
                   PFields(1) % Hdr % NumRows ) )

ALLOCATE ( Upper ( PFields(1) % Hdr % NumCols,     &
                   PFields(1) % Hdr % NumRows ) )

ALLOCATE ( MD_flag ( PFields(1) % Hdr % NumCols,   &
                     PFields(1) % Hdr % NumRows ) )

ALLOCATE ( cld_mask ( PFields(1) % Hdr % NumCols,  &
                      PFields(1) % Hdr % NumRows ) )


DO j = 1,PFields(1) % Hdr % NumRows
  DO i = 1,PFields(1) % Hdr % NumCols
    IF (PRef_Pa >= PFields(1) % RData(i,j)) THEN
      Lower(i,j)=1
      Upper(i,j)=2
    ELSE IF ( PRef_Pa < PFields(NumLevs) % RData (i,j)) THEN
      Lower(i,j)=NumLevs-1
      Upper(i,j)=NumLevs
    ELSE
      Lower(i,j)=-1 ! one level to find
    END IF
  END DO
END DO

DO k=2,NumLevs
   DO j = 1,PFields(1) % Hdr % NumRows
      DO i = 1,PFields(1) % Hdr % NumCols
         IF ((PFields(k) % RData(i,j) <= PRef_Pa).and. &
             (Lower(i,j).eq.-1)) then
            Lower(i,j)=k-1
            Upper(i,j)=k
         END IF
      END DO
   END DO
END DO

DO j = 1,PFields(1) % Hdr % NumRows
   DO i = 1,PFields(1) % Hdr % NumCols
      IF (  PFields(Lower(i,j)) % RData(i,j) < PRef_Pa &
           .OR. Lower(i,j) < 3 ) THEN
         MD_flag(i,j) = 1
      ELSE
         MD_flag(i,j) = 0
      END IF
   END DO
END DO

!------------------------------------------------------
!Construct cloud mask. Do this by checking that, for each gridpoint, there
! either the model level above or below this pressure level; if so mark cld
! if not mark as 0.

DO j = 1,CldFields(1) % Hdr % NumRows
  DO i = 1,CldFields(1) % Hdr % NumCols

    IF (CldFields(Lower(i,j))%RData(i,j) == CldFields(1)%Hdr%BMDI.OR. &
         CldFields(Upper(i,j))%RData(i,j) == CldFields(1)%Hdr%BMDI.OR.&
         MD_flag(i,j) == 1 ) THEN
      cld_mask(i,j) = IMDI
    ELSE IF ( CldFields( Lower (i,j)) % RData(i,j) > 0.0 .OR.         &
              CldFields( Upper (i,j)) % RData(i,j) > 0.0 ) THEN
      cld_mask(i,j) = 1
    ELSE
      cld_mask(i,j) = 0
    END IF

  END DO
END DO


!----------------------------------------------------------
! Find vertical derivatives of z and equivalent pot temp wrt pressure using
! cubic spline. Note that spline values must be monotonically increasing.
! Pressure level supplied to DiffP must be in Pa.

! DEPENDS ON: diffp
CALL DiffP( NumLevs, PRef_Pa, ZFields, PFields, dZdP, ErrorStatus )
! DEPENDS ON: diffp
CALL DiffP( NumLevs, PRef_Pa, theta_e_Fields, PFields, dthetaedP, &
            ErrorStatus )

dthetaedZ % Hdr = theta_e_Fields(1) % Hdr
ALLOCATE( dthetaedZ % RData(dthetaedZ % Hdr % NumCols, &
                            dthetaedZ % Hdr % NumRows) )

WHERE (  cld_mask /=  IMDI )
  dthetaedZ % RData = dthetaedP % RData / dZdP % RData
ELSEWHERE
  dthetaedZ % RData = RMDI
END WHERE


!--------------------------------------------------------------
!Loop over all gridpoints at this level & set cloud turbulence predictor.
! If orography makes calculation impossible, set cld_turb_Pred
! to missing data indicator. If calculation is possible,
! check cloud mask=1. If so, check sign of dthetaedZ. If negative,
! mark as turbulent (set to dthetadz); if positive mark as
! non-turbulent (zero).

DO j = 1, cld_turb_Pred % Hdr % NumRows
  DO i = 1, cld_turb_Pred % Hdr % NumCols

    IF (cld_mask(i,j) == IMDI ) THEN
      cld_turb_Pred % RData(i,j) = cld_turb_Pred % Hdr % BMDI
    ELSE IF (cld_mask (i,j) == 1 .AND. &
             dthetaedZ % RData(i,j) < 0.0 ) THEN
      cld_turb_Pred % RData(i,j) = ABS(dthetaedZ % RData(i,j))
    ELSE
      cld_turb_Pred % RData(i,j) = 0.0
    END IF

  END DO
END DO


DEALLOCATE( dZdP % RData )
DEALLOCATE( dthetaedP % RData )
DEALLOCATE( dthetaedZ % RData )
DEALLOCATE( Lower )
DEALLOCATE( Upper )
DEALLOCATE( MD_flag )
DEALLOCATE( cld_mask )

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )


END SUBROUTINE cld_turb
