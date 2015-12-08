! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+ Module containing derived types used for IO in Fieldcalc

MODULE IO_Mod

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

USE Ereport_mod, ONLY : Ereport
IMPLICIT NONE

! Global Constants:
INTEGER, PARAMETER :: LenFixHd = 256   ! Length of Fixed_Length_Header
INTEGER, PARAMETER :: LenWord  = 64

! Global Type Definitions:
TYPE PP_Header_type
  SEQUENCE
  INTEGER :: ValidYear    = -99  ! 1  LBYR   : Year       \ .
  INTEGER :: ValidMonth   = -99  ! 2  LBMON  : Month       \ .
  INTEGER :: ValidDate    = -99  ! 3  LBDAT  : Day          \  Validity time
  INTEGER :: ValidHour    = -99  ! 4  LBHR   : Hour         / .
  INTEGER :: ValidMin     = -99  ! 5  LBMIN  : Minute      / .
  INTEGER :: ValidSec     = -99  ! 6  LBSEC  : Second     / .
  INTEGER :: DataYear     = -99  ! 7  LBYRD  : Year       \ .
  INTEGER :: DataMonth    = -99  ! 8  LBMOND : Month       \ .
  INTEGER :: DataDate     = -99  ! 9  LBDATD : Day          \ Data time
  INTEGER :: DataHour     = -99  ! 10 LBHRD  : Hour         / .
  INTEGER :: DataMin      = -99  ! 11 LBMIND : Minute      / .
  INTEGER :: DataSec      = -99  ! 12 LBSECD : Second     / .
  INTEGER :: LBTim        = -99  ! 13          Time indicator
  INTEGER :: FCRange      = -99  ! 14 LBFT   : Forecast range (hrs)
  INTEGER :: LBLRec       = -99  ! 15          Length of data record
  INTEGER :: LBCode       = -99  ! 16          Grid type code
  INTEGER :: LBHem        = -99  ! 17          Hemisphere indicator
  INTEGER :: NumRows      = -99  ! 18 LBROW  : Number of rows in grid
  INTEGER :: NumCols      = -99  ! 19 LBNPT  : Number of columns in grid
  INTEGER :: LBExt        = -99  ! 20          Length of extra data
  INTEGER :: LBPack       = -99  ! 21          Packing method indicator
  INTEGER :: LBRel        = -99  ! 22          Header release number
  INTEGER :: PPCode       = -99  ! 23 LBFC   : Field code
  INTEGER :: LBCFC        = -99  ! 24          Second field code
  INTEGER :: LBProc       = -99  ! 25          Processing code
  INTEGER :: LBVC         = -99  ! 26          Vertical co-ordinate type
  INTEGER :: LBRVC        = -99  ! 27          Vert co-ord type for reference
  INTEGER :: LBExp        = -99  ! 28          Experiment number
  INTEGER :: DataPos      = -99  ! 29 LBEGIN : Word address of data
  INTEGER :: LBNRec       = -99  ! 30          Disk length / No. records
  INTEGER :: MO8Proj      = -99  ! 31 LBPROJ : MetO8 Projection number
  INTEGER :: MO8Type      = -99  ! 32 LBTYP  : MetO8 Field Type
  INTEGER :: MO8Level     = -99  ! 33 LBLEV  : MetO8 Level code
  INTEGER :: LBRsvd(4)    = -99  ! 34-37       PP-Package use
  INTEGER :: LBSrce       = -99  ! 38          Version & Model
  INTEGER :: LBUser1      = -99  ! 39          Data type (real, etc.) (UM only)
  INTEGER :: LBUser2      = -99  ! 40          Start address in DATA (UM only)
  INTEGER :: LBUser3      = -99  ! 41          No. periods for timeseries
  INTEGER :: STCode       = -99  ! 42          STASH code (UM only)
  INTEGER :: LBUser5      = -99  ! 43          STASH pseudo dimension (UM only)
  INTEGER :: LBUser6      = -99  ! 44          Free (used by FIELDCOS)
  INTEGER :: LBUser7      = -99  ! 45          Sub-model number (UM only)
  REAL    :: BULev        = -99. ! 46          Upper layer boundary
  REAL    :: BHULev       = -99. ! 47          Upper layer boundary
  REAL    :: BRsvd(2)     = -99. ! 48-49       PP-Package use
  REAL    :: BDatum       = -99. ! 50          Datum value
  REAL    :: BAcc         = -99. ! 51          Packing accuracy
  REAL    :: RLevel       = -99. ! 52 BLEV   : Level (B or Zsea for hybrid)
  REAL    :: RefLevel     = -99. ! 53 BRLEV  : Reference level - lower layer bd
  REAL    :: BHLev        = -99. ! 54          Level (A or C for hybrid)
  REAL    :: BHRLev       = -99. ! 55          Lower layer boundary
  REAL    :: PseudoLat    = -99. ! 56 BPLAT  : Real lat of 'pseudo' N Pole
  REAL    :: PseudoLon    = -99. ! 57 BPLON  : Real lon of 'pseudo' N Pole
  REAL    :: BGOR         = -99. ! 58          Grid orientation
  REAL    :: ZerothLat    = -99. ! 59 BZY    : Zeroth latitude
  REAL    :: LatInt       = -99. ! 60 BDY    : Latitude interval
  REAL    :: ZerothLon    = -99. ! 61 BZX    : Zeroth longitude
  REAL    :: LonInt       = -99. ! 62 BDX    : Longitude interval
  REAL    :: BMDI         = -99. ! 63          Missing data indicator
  REAL    :: BMKS         = -99. ! 64          M.K.S. scaling factor

END TYPE PP_Header_type

TYPE PP_Field_type

  INTEGER              :: LookupPos   ! No of header in lookup
  INTEGER              :: ArrayPos    ! No of header in Fields array
  REAL, POINTER        :: RData(:,:) => NULL()
  TYPE(PP_Header_type) :: Hdr

END TYPE PP_Field_type

TYPE UM_Header_type

  INTEGER :: LenIntC       ! Length of Integer_Constants array
  INTEGER :: LenRealC      ! Length of Real_Constants array
  INTEGER :: Len1LevDepC   ! 1st dim \ Level_Dependent_Constants
  INTEGER :: Len2LevDepC   ! 2nd dim / array
  INTEGER :: Len1RowDepC   ! 1st dim \ Row_Dependent_Constants
  INTEGER :: Len2RowDepC   ! 2nd dim / array
  INTEGER :: Len1ColDepC   ! 1st dim \ Column_Dependent_Constants
  INTEGER :: Len2ColDepC   ! 2nd dim / array
  INTEGER :: Len1FldsOfC   ! 1st dim \ Fields_Of_Constants
  INTEGER :: Len2FldsOfC   ! 2nd dim / array
  INTEGER :: LenExtraC     ! Length of Extra_Constants array
  INTEGER :: LenHistFile   ! Length of Temp_History_File
  INTEGER :: LenCompFldI1  ! Length of Compressed_Field_Index1
  INTEGER :: LenCompFldI2  ! Length of Compressed_Field_Index2
  INTEGER :: LenCompFldI3  ! Length of Compressed_Field_Index3
  INTEGER :: Len1Lookup    ! 1st dim \ Lookup table
  INTEGER :: Len2Lookup    ! 2nd dim /
  INTEGER :: LenData       ! Length of Data array
  INTEGER :: StartData     ! Position of start of Data array
  INTEGER :: NumFlds       ! Number of Data fields
  INTEGER :: UnitNum       ! Unit number associated with UM dump

  INTEGER, POINTER :: FixHd(:)     => NULL()   ! Fixed_Length_Header
  INTEGER, POINTER :: IntC(:)      => NULL()   ! Integer_Constants array
  INTEGER, POINTER :: CompFldI1(:) => NULL()   ! Compressed_Field_Index1 array
  INTEGER, POINTER :: CompFldI2(:) => NULL()   ! Compressed_Field_Index2 array
  INTEGER, POINTER :: CompFldI3(:) => NULL()   ! Compressed_Field_Index3 array
  TYPE(PP_Header_type), POINTER :: Lookup(:) => NULL() ! Lookup table

  REAL, POINTER :: RealC(:)    => NULL()
  REAL, POINTER :: LevDepC(:)  => NULL()
  REAL, POINTER :: RowDepC(:)  => NULL()
  REAL, POINTER :: ColDepC(:)  => NULL()
  REAL, POINTER :: FldsOfC(:)  => NULL()
  REAL, POINTER :: ExtraC(:)   => NULL()
  REAL, POINTER :: HistFile(:) => NULL()

  CHARACTER(LEN=16) :: FileNameEnv

END TYPE UM_Header_type

! Needed to unpack land/sea-packed fields - requires SAVE attribute due to
! default initialisation within PP_Field_type and is inside a module.
TYPE(PP_Field_type), SAVE  :: LSMField       ! Land-sea mask PP-Field

CONTAINS
  SUBROUTINE buffin_um_lookup(header,length,error)
    USE IO, ONLY :              &
        buffin,                 &
        setpos
    USE Err_Mod, ONLY:          &
        StatusOK,               &
        StatusWarning
    USE ereport_mod, ONLY : ereport
    IMPLICIT NONE
    TYPE(UM_Header_type), INTENT(INOUT) :: header
    INTEGER, INTENT(OUT) :: length
    REAL, INTENT(OUT)    :: error
    INTEGER              :: ErrorStatus

    INTEGER, ALLOCATABLE :: buffer(:)

    CALL setpos ( header % UnitNum, header % FixHd(150)-1, ErrorStatus )
    IF ( ErrorStatus /= StatusOK ) THEN
      CALL EReport( "io_mod::buffin_um_lookup",   &
          ErrorStatus, "Failure in SETPOS" )
    END IF

    ALLOCATE(buffer(header % Len1Lookup * header % Len2Lookup))
    CALL buffin(header % UnitNum, buffer,         &
        header % Len1Lookup * header % Len2Lookup,length,error)
    CALL um_memcpy_f( header % lookup , buffer ,  &
        header % Len1Lookup * header % Len2Lookup)
    DEALLOCATE(buffer)

  END SUBROUTINE buffin_um_lookup


  SUBROUTINE buffout_um_lookup(header,length,error)
    USE IO, ONLY :              &
        buffout,                &
        setpos
    USE Err_Mod, ONLY:          &
        StatusOK,               &
        StatusWarning
    USE ereport_mod, ONLY : ereport
    IMPLICIT NONE
    TYPE(UM_Header_type), INTENT(IN) :: header
    INTEGER, INTENT(OUT) :: length
    REAL, INTENT(OUT)    :: error
    INTEGER              :: ErrorStatus

    INTEGER, ALLOCATABLE :: buffer(:)

    CALL setpos ( header % UnitNum, header % FixHd(150)-1, ErrorStatus )
    IF ( ErrorStatus /= StatusOK ) THEN
      CALL EReport( "io_mod::buffout_um_lookup",   &
          ErrorStatus, "Failure in SETPOS" )
    END IF

! Even though we only have "header % NumFlds" values we also want to write
! out the unused lookup to allow other programs to recognise when the lookups
! end.
    ALLOCATE(buffer(header % Len1Lookup * header % Len2Lookup))
    CALL um_memcpy_f(buffer, header % lookup ,     &
        header % Len1Lookup * header % Len2Lookup)
    CALL buffout(header % UnitNum, buffer,         &
        header % Len1Lookup * header % Len2Lookup,length,error)
    DEALLOCATE(buffer)

  END SUBROUTINE buffout_um_lookup

END MODULE IO_Mod

