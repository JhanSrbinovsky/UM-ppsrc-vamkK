! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to calculate snow probability diagnostic

SUBROUTINE SnowProb( ZStd1000,     &  ! in
                     ZStd0850,     &  ! in
                     SnowPr,       &  ! inout
                     ErrorStatus )    ! inout

! Description:
!
! Method:
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
  StatusOK
USE FldCodes_Mod, ONLY:   &
  ST_SnowPr, MO8_SnowPr,  &
  PP_Prob,   VC_Surface,  &
  LV_Special
IMPLICIT None

! Subroutine Arguments:
TYPE(PP_Field_type), INTENT(IN) :: ZStd1000   ! 1000hPa Height
TYPE(PP_Field_type), INTENT(IN) :: ZStd0850   !  850hPa Height

TYPE(PP_Field_type), INTENT(INOUT) :: SnowPr
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "SnowProb"
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

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( ASSOCIATED( SnowPr % RData ) ) THEN
  DEALLOCATE( SnowPr % RData )
END IF
SnowPr % Hdr = ZStd1000 % Hdr
SnowPr % Hdr % PPCode   =  PP_Prob
SnowPr % Hdr % LBVC     =  VC_Surface
SnowPr % Hdr % MO8Type  = MO8_SnowPr
SnowPr % Hdr % MO8Level =  LV_Special
SnowPr % Hdr % STCode   =  ST_SnowPr
SnowPr % Hdr % RLevel   = 0.0
SnowPr % Hdr % BMDI     = RMDI
ALLOCATE( SnowPr % RData(SnowPr % Hdr % NumCols, &
                         SnowPr % Hdr % NumRows) )

SnowPr % RData = 5220.0 + 3.86666 * ZStd1000 % RData &
                        -     4.0 * ZStd0850 % RData

WHERE ( ZStd1000 % RData >=  1.e8 )
  SnowPr % RData = 0.0
END WHERE

! Limit percentage probability
WHERE (  SnowPr  % RData <=    0.0 )
  Snowpr % RData = 0.0
END WHERE
WHERE (  SnowPr  % RData >=  100.0 )
  SnowPr % RData = 100.0
END WHERE

9999 CONTINUE

END SUBROUTINE SnowProb
