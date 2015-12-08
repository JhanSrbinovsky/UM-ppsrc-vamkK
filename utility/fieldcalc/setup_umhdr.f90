! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Routine to allocate space for a UM header.

SUBROUTINE Setup_UMHdr ( UMHdr )    ! inout

! Description:
!   This is a subroutine to allocate memory for the various arrays
!   required by a UM header.
!
! Method:
!   The UM header passed in to the subroutine must already contain the
!   fixed length header, either read from a file or defined by another
!   subroutine.  The values within the FLH are then used to allocate
!   arrays to store the rest of the UM header.  The lookup table and
!   data arrays are allocated elsewhere, to save space.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY: &
  LenFixHd,       &
  UM_Header_type

IMPLICIT None

! Subroutine Arguments:
TYPE(UM_Header_type), INTENT(INOUT) :: UMHdr

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Setup_UMHdr"

! End of header --------------------------------------------------------

IF ( UMHdr % FixHd(100) > 0 ) THEN     ! Integer Consts array present
  UMHdr % LenIntC = UMHdr % FixHd(101)
ELSE                                   ! No Integer Consts array
  UMHdr % LenIntC = 0
END IF
ALLOCATE ( UMHdr % IntC(UMHdr % LenIntC) )

IF ( UMHdr % FixHd(105) > 0 ) THEN     ! Real Consts array present
  UMHdr % LenRealC = UMHdr % FixHd(106)
ELSE                                   ! No Real Consts array
  UMHdr % LenRealC = 0
END IF
ALLOCATE ( UMHdr % RealC(UMHdr % LenRealC) )

IF ( UMHdr % FixHd(110) > 0 ) THEN     ! Lev Dep Consts array present
  UMHdr % Len1LevDepC = UMHdr % FixHd(111)
  UMHdr % Len2LevDepC = UMHdr % FixHd(112)
ELSE                                   ! No Lev Dep Consts array
  UMHdr % Len1LevDepC = 0
  UMHdr % Len2LevDepC = 0
END IF
ALLOCATE ( UMHdr % LevDepC(1 + UMHdr % Len1LevDepC * &
                               UMHdr % Len2LevDepC) )

IF ( UMHdr % FixHd(115) > 0 ) THEN     ! Row Dep Consts array present
  UMHdr % Len1RowDepC = UMHdr % FixHd(116)
  UMHdr % Len2RowDepC = UMHdr % FixHd(117)
ELSE                                   ! No Row Dep Consts array
  UMHdr % Len1RowDepC = 0
  UMHdr % Len2RowDepC = 0
END IF
ALLOCATE ( UMHdr % RowDepC(1 + UMHdr % Len1RowDepC * &
                               UMHdr % Len2RowDepC) )

IF ( UMHdr % FixHd(120) > 0 ) THEN     ! Col Dep Consts array present
  UMHdr % Len1ColDepC = UMHdr % FixHd(121)
  UMHdr % Len2ColDepC = UMHdr % FixHd(122)
ELSE                                   ! No Col Dep Consts array
  UMHdr % Len1ColDepC = 0
  UMHdr % Len2ColDepC = 0
END IF
ALLOCATE ( UMHdr % ColDepC(1 + UMHdr % Len1ColDepC * &
                               UMHdr % Len2ColDepC) )

IF ( UMHdr % FixHd(125) > 0 ) THEN     ! Flds Of Consts array present
  UMHdr % Len1FldsOfC = UMHdr % FixHd(126)
  UMHdr % Len2FldsOfC = UMHdr % FixHd(127)
ELSE                                   ! No Flds Of Consts array
  UMHdr % Len1FldsOfC = 0
  UMHdr % Len2FldsOfC = 0
END IF
ALLOCATE ( UMHdr % FldsOfC(1 + UMHdr % Len1FldsOfC * &
                               UMHdr % Len2FldsOfC) )

IF ( UMHdr % FixHd(130) > 0 ) THEN     ! Extra Consts array present
  UMHdr % LenExtraC = UMHdr % FixHd(131)
ELSE                                   ! No Extra Consts array
  UMHdr % LenExtraC = 0
END IF
ALLOCATE ( UMHdr % ExtraC(1 + UMHdr % LenExtraC) )

IF ( UMHdr % FixHd(135) > 0 ) THEN     ! Hist file array present
  UMHdr % LenHistFile = UMHdr % FixHd(136)
ELSE                                   ! No Hist file array
  UMHdr % LenHistFile = 0
END IF
ALLOCATE ( UMHdr % HistFile(1 + UMHdr % LenHistFile) )

IF ( UMHdr % FixHd(140) > 0 ) THEN     ! Comp Field Index1 array present
  UMHdr % LenCompFldI1 = UMHdr % FixHd(141)
ELSE                                   ! No Comp Field Index1 array
  UMHdr % LenCompFldI1 = 0
END IF
ALLOCATE ( UMHdr % CompFldI1(1 + UMHdr % LenCompFldI1) )

IF ( UMHdr % FixHd(142) > 0 ) THEN     ! Comp Field Index2 array present
  UMHdr % LenCompFldI2 = UMHdr % FixHd(143)
ELSE                                   ! No Comp Field Index2 array
  UMHdr % LenCompFldI2 = 0
END IF
ALLOCATE ( UMHdr % CompFldI2(1 + UMHdr % LenCompFldI2) )

IF ( UMHdr % FixHd(143) > 0 ) THEN     ! Comp Field Index3 array present
  UMHdr % LenCompFldI3 = UMHdr % FixHd(144)
ELSE                                   ! No Comp Field Index3 array
  UMHdr % LenCompFldI3 = 0
END IF
ALLOCATE ( UMHdr % CompFldI3(1 + UMHdr % LenCompFldI3) )

IF ( UMHdr % FixHd(150) > 0 ) THEN     ! Lookup Table present
  UMHdr % Len1Lookup = UMHdr % FixHd(151)
  UMHdr % Len2Lookup = UMHdr % FixHd(152)
  NULLIFY ( UMHdr % Lookup )
  ! Allocate the Lookup table later, to save space
ELSE                                   ! No Lookup Table
  UMHdr % Len1Lookup = 0
  UMHdr % Len2Lookup = 0
  NULLIFY ( UMHdr % Lookup )
  ! Possibly report error here?  No fields to read?
END IF

IF ( UMHdr % FixHd(160) > 0 ) THEN     ! Data present
  UMHdr % StartData = UMHdr % FixHd(160)
! LenData is rarely used and doesnt get updated (unless we ever write dumps).
  UMHdr % LenData   = UMHdr % FixHd(161)
ELSE                                   ! No Data present
  UMHdr % StartData = 0
  UMHdr % LenData   = 0
END IF

END SUBROUTINE Setup_UMHdr
