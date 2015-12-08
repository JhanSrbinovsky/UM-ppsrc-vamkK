! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routines to calculate Clear Air Turbulence diagnostics.
!=======================================================================

SUBROUTINE MaxCATurb( CAT300,       &  ! in
                      CAT250,       &  ! in
                      CAT200,       &  ! in
                      MaxCAT,       &  ! inout
                      MaxCATLev,    &  ! inout
                      ErrorStatus )    ! inout

! Description:
!   Subroutine to find the maximum CAT value from the 3 standard levels
!
! Method:
!   The 3 standard CAT levels are passed in already calculated, and this
!   subroutine simply looks at each field in turn to find the maximum
!   value.
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
USE FldCodes_mod, ONLY:             &
  ST_MxCAT,   MO8_MxCAT,   PP_CAT,  &
  ST_MxCATP,  MO8_MxCATP,  PP_P,    &
  VC_Turb,    LV_Special
IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: CAT300          ! Cat predictor at
TYPE(PP_Field_type), INTENT(IN) :: CAT250          !   300, 250 and
TYPE(PP_Field_type), INTENT(IN) :: CAT200          !   200 hPa

TYPE(PP_Field_type), INTENT(INOUT) :: MaxCAT       ! Max CAT value
TYPE(PP_Field_type), INTENT(INOUT) :: MaxCATLev    ! Max CAT level
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "MaxCATurb"
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

! End of header --------------------------------------------------------

! DEPENDS ON: timer
CALL Timer( RoutineName, 3 )

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( MaxCAT % RData ) ) THEN
  DEALLOCATE( MaxCAT % RData )
END IF
IF ( ASSOCIATED( MaxCATLev % RData ) ) THEN
  DEALLOCATE( MaxCATLev % RData )
END IF
! Set up headers & allocate memory
MaxCAT % Hdr = CAT300 % Hdr
MaxCAT % Hdr % LBVC     =  VC_Turb
MaxCAT % Hdr % PPCode   =  PP_CAT
MaxCAT % Hdr % MO8Type  =  MO8_MxCAT
MaxCAT % Hdr % MO8Level =  LV_Special
MaxCAT % Hdr % STCode   =  ST_MxCAT
MaxCAT % Hdr % RLevel   =  0.0
MaxCAT % Hdr % BMDI     =  RMDI
MaxCATLev % Hdr = MaxCAT % Hdr
MaxCATLev % Hdr % PPCode   =  PP_P
MaxCATLev % Hdr % MO8Type  = MO8_MxCATP
MaxCATLev % Hdr % STCode   =  ST_MxCATP
ALLOCATE( MaxCAT % RData(MaxCAT % Hdr % NumCols, &
                         MaxCAT % Hdr % NumRows) )
ALLOCATE( MaxCATLev % RData(MaxCATLev % Hdr % NumCols, &
                            MaxCATLev % Hdr % NumRows) )

MaxCAT    % RData(:,:) = CAT300 % RData(:,:)
MaxCATLev % RData(:,:) = 100.0 * CAT300 % Hdr % RLevel    ! Pa (SI unit)
WHERE ( CAT250 % RData(:,:) > MaxCAT % RData(:,:) )
  MaxCAT    % RData(:,:) = CAT250 % RData(:,:)
  MaxCATLev % RData(:,:) = 100.0 * CAT250 % Hdr % RLevel  ! Pa
END WHERE
WHERE ( CAT200 % RData(:,:) > MaxCAT % RData(:,:) )
  MaxCAT    % RData(:,:) = CAT200 % RData(:,:)
  MaxCATLev % RData(:,:) = 100.0 * CAT200 % Hdr % RLevel  ! Pa
END WHERE

9999 CONTINUE

! DEPENDS ON: timer
CALL Timer( RoutineName, 4 )

END SUBROUTINE MaxCATurb

