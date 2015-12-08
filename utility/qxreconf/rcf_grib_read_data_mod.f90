! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Read, from file, and Decode 1 field's worth of GRIB data

MODULE Rcf_Grib_Read_Data_Mod

IMPLICIT NONE

! SUBROUTINE Rcf_Grib_Read_Data : Read raw info from file
!
! Description:
!   This routine is used to handle the reading of raw GRIB encoded data
!   in. This data is then passed to DECODE to be decoded.
!
! Method:
!   Recieve position markers and Storage arrays (or Types) from
!   calling routine.
!   Use SETPOS8 to set position in file to start of GRIB record.
!   Read a block of Raw data.
!   Call DECODE to decode raw Data.
!   Store the decoded data
!   Use the decoded length to calculate the next start point and move on
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
CONTAINS
SUBROUTINE Rcf_Grib_Read_Data(Unit_Num ,Current,Data,Len_Data,        &
                              Pos_in_File)

USE Rcf_GRIB_Block_Params_Mod, Only : &
    Grib_Record,     &       ! derived type for storing a GRIB header
    LenArrayMax,     &       ! Max size of buffer for data
    p_B4Undef

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Min,           &       ! =1 Minimum output
    PrStatus_Normal,        &       ! =2 Short informative output
    PrStatus_Oper,          &       ! =3 Full informative output
    PrStatus_Diag                   ! =4 Extra Diagnostic output

USE UM_ParVars, Only : &
    mype

USE EReport_Mod, Only :     &
    EReport                         ! Subroutine to generate error msgs

USE IO  

USE ukmo_grib_mod, Only : ukmo_decode

IMPLICIT NONE

! Subroutine arguments

!< Scalar arguments with intent(in):>
INTEGER, INTENT(IN)              :: Unit_Num  ! unit number file is on
INTEGER, INTENT(IN)              :: Len_Data  ! Length of data to read
INTEGER, INTENT(IN)              :: Pos_in_File ! position to start
                                                ! reading data from

!< Array  arguments with intent(InOut):>
TYPE (Grib_Record),POINTER       :: Current        !\ Pointer to current
                                                   !/ grib record
!< Array  arguments with intent(out):>
REAL, INTENT(OUT)                :: Data(Len_Data) ! Array deocded field
                                                   ! is returned in.

! Local constants
! This was usedfor debugging. It allowed all non debugging messages to
! be turned off.
INTEGER , PARAMETER              :: PrStatus_Debug = 0 ! debug

! Local variables

CHARACTER (LEN=*), PARAMETER     :: RoutineName='Rcf_Grib_Read_Data'

CHARACTER (LEN=80)               :: Cmessage      ! used for EReport
INTEGER                          :: ErrorStatus   ! used for EReport
INTEGER                          :: len_read_data ! used for buffin
INTEGER                          :: Error         ! Error Status
REAL                             :: RError         ! Error Status
INTEGER                          :: act_io_len

! Comdecks
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

!=======================================================================
!  Variables Usedfor calls to 'DECODE'
!=======================================================================

INTEGER , PARAMETER    :: LenVertCord = 3000    ! Len of Vert_Coords
INTEGER , PARAMETER    :: LenBitmap   = LenArrayMax ! Len of bitmap used
INTEGER , PARAMETER    :: LenQuasi    =    1    ! Len of array Quasi
INTEGER , PARAMETER    :: Word_Size   =   64    ! Word size for machine
                                                ! 64 bit on Cray else 32

! LenWork1 and LenWorkR should be >= to no. of columns in grid
! LenWork2 should be >= twice the no. of rows in grid.
INTEGER , PARAMETER    :: LenWork1    =  288   ! Len of Work array 1
INTEGER , PARAMETER    :: LenWork2    =  500   ! Len of Work array 2
INTEGER , PARAMETER    :: LenWorkR    =  288   ! Len of Work array R
INTEGER , PARAMETER    :: LenPosn     =    4   ! Len of array Posn

!Variables used when calling DECODE
INTEGER                :: Width                ! No. of bits used to
                                               ! encode field els
INTEGER                :: MsgLvl               ! Message level in
                                               ! DECODE

INTEGER                :: Words                ! no of words
INTEGER, PARAMETER     :: Error_Unit =     6   ! Error unit used by
                                               ! DECODE

INTEGER                :: Bitmap(LenBitmap)    ! Array for bitmap of
                                               ! non full field data
INTEGER                :: Quasi(LenQuasi)      ! desc of Quasi-reg grid
INTEGER                :: I_Record(Len_Data)   ! Int array holding
                                               ! Character array of
                                               !encoded GRIB data
INTEGER                :: Work_Int1(LenWork1)  ! Work Array
INTEGER                :: Work_Int2(LenWork2)  ! Work Array

INTEGER                :: Off                  ! Word Offset - not used

REAL                   :: FpWork(Len_Data)      ! Work Array
REAL                   :: VertCoords(LenVertCord) ! Vertical coord params
INTEGER                :: Posn(LenPosn)        ! Not Used- reqd 4 decode
REAL                   :: Work_Re1(LenWorkR)   ! Work Array
CHARACTER(LEN=1)       :: C_Record(Len_Data) ! Encoded data read into
                                               ! this using Buffin
! Lets size this to 2 8-byte INTEGERs to support GRIB1 and GRIB2.
INTEGER                :: section0(2)
INTEGER                :: grib_version

!=======================================================================
!  Read binary form of record into CRecord work array
!=======================================================================

C_Record(:) = " "                           ! Make sure buffer is clear
act_io_len  = 0

!Set MsgLvl used by Decode - almost the reverse of rcf system
SELECT CASE (PrintStatus)
CASE (PrStatus_Min,PrStatus_Normal)
  MsgLvl = 3          ! No messages
CASE (PrStatus_Oper)
  MsgLvl = 1          ! Errors and Warnings
CASE (PrStatus_Diag)
  MsgLvl = 0          ! Errors, Warnings and Notes
END SELECT
! there is also level 2, 'Errors Only' available.

! Set Position in file to be read
CALL SetPos8( Unit_Num, Pos_in_File, Error)

IF ( Error /= 0 ) THEN
  WRITE(Cmessage,'(A,I8)') 'Failed trying to SetPos to ',        &
                               Pos_in_File
  ErrorStatus = 30
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

! Read the section 0 block of data into a 8 byte integer.
CALL Buffin (Unit_Num, section0,                 &
             2, act_io_len, RError)

IF ( RError > 0.0 ) THEN
  Cmessage    = 'Failed trying to read section 0 from GRIB file'
  ErrorStatus = 40
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

len_read_data = 0
grib_version = IBITS(section0(1),0,8)
IF (grib_version == 1) THEN
! Within section 0 of a GRIB file 5:7 bytes contains the length.
  len_read_data = IBITS(section0(1), 8, 24)
ELSE IF (grib_version == 2) THEN
! Within section 0 of a GRIB2 file 9:16 bytes contains the length.
  len_read_data = section0(2)
ELSE
  WRITE(6,'(A,I0)') "GRIB version found is ", grib_version
  Cmessage    = 'Unknown GRIB version whilst extracting length of message.'
  ErrorStatus = 41
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

IF (len_read_data > len_data) THEN
  WRITE(6,'(A,I0)') "Length of message is ", len_read_data
  WRITE(6,'(A,I0)') "Maximum buffer size is ", len_data
  Cmessage    = 'GRIB message is bigger than expected maximum buffer size'
  ErrorStatus = 41
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

! Set Position in file to be read again to reread data.
CALL SetPos8( Unit_Num, Pos_in_File, Error)

! Read the block of data into CRecord
CALL Buffin (Unit_Num, C_Record(:),                 &
             len_read_data, 1, act_io_len, RError)

IF ( RError > 0.0 ) THEN
  Cmessage    = 'Failed trying to read record from GRIB file'
  ErrorStatus = 42
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

! Report on progress
IF ( PrintStatus >= PrStatus_Diag .AND. mype == 0 ) THEN
    WRITE (6,'(A,I8,A,I8)') "Whilst trying to read full record tried for ", &
                 Len_read_Data, " and got ", act_io_len + 8
END IF


!=======================================================================
!  Prepare for and Call DECODE routine
!=======================================================================
! Initialise some of the arrays used in the call to DECODE
Posn(:)               = 0
Words                 = 0
Off                   = 0
Error                 = 0
I_Record(:)           = 0
Current % Num_Fp      = 0
Current % Num_vert    = 0
Current % Num_Bitmap  = 0
Current % Num_Quasi   = 0

! Despite not being used DECODE expects this variable to be 0
Current % Block_4(p_B4Undef) = 0

! Transfer the character array that the encoded data is in into
! an Integer array which DECODE expects
! After some initial problems with Transfer it was found to be more
! reliable when 'size' was supplied. careful note must be made of the
! relative sizes of the 2 types. On the Cray Int >= Char(len=1) so
! this should not cause problems.

I_Record(:Len_read_Data) = Transfer(C_Record(:Len_read_Data),                 &
                                    I_Record(1), Len_read_Data)
CALL ukmo_decode(Data, FpWork, Len_Data,    &
            Current % Num_Fp,               &
            VertCoords, LenVertCord,        &  !
            Current % Num_vert,             &
            Bitmap, LenBitmap,              &
            Current % Num_Bitmap,           &
            Quasi, LenQuasi,                &
            Current % Num_Quasi,            &
            Width, Word_Size,               &
            Current % Block_0,              &  !
            Current % Block_1,              &  !
            Current % Block_2,              &  !
            Current % Block_3,              &  !
            Current % Block_4,              &
            Current % Block_R,              &  !
            I_Record, act_io_len, Posn, Words, Off,  &
            Error, Work_Int1, Work_Int2, Work_Re1,  &
            Error_Unit, MsgLvl )

IF ( Error /= 0 ) THEN
  WRITE (Cmessage,'(A,I2,A)') "Error ",Error," returned by DECODE"

  IF (Error == 3) THEN       ! DECODE is reporting garbage out
    ErrorStatus = 50
    Cmessage = TRIM(Cmessage) // " -Garbage Out-"
    CALL EReport( RoutineName, ErrorStatus, Cmessage )

  ELSE IF (Error == 2) THEN  ! Decode is complaining it's an extended
                             ! table 2 (e.g. ECMWF)
    IF ( PrintStatus >= PrStatus_Diag ) THEN
      ErrorStatus = -55        ! DECODE is reporting warnings
      Cmessage = TRIM(Cmessage) // " -Extended Table 2 (ECMWF ?)-"
      CALL EReport( RoutineName, ErrorStatus, Cmessage )
    END IF
  ELSE
    ErrorStatus = -50        ! DECODE is reporting warnings
    CALL EReport( RoutineName, ErrorStatus, Cmessage )
  END IF
END IF

! If vertical coordinates exist, store in grib record
IF (Current % Num_vert /= 0) THEN
  IF ( mype == 0 ) WRITE(6,'(A)') "Storing vertical coordinates"
  ALLOCATE(Current % VertCoords(Current % Num_vert))
  Current % VertCoords=VertCoords(1:Current % Num_vert)
END IF

! Set missing data values properly in bitmaps
IF (Current % Num_Bitmap /= 0) THEN
  IF ( mype == 0 ) WRITE(6,'(A)') "Setting MDI properly for bitmap fields"
  WHERE (Data == Current % Block_R(2)) Data=rmdi
END IF

RETURN

END SUBROUTINE Rcf_Grib_Read_Data
END MODULE Rcf_Grib_Read_Data_Mod
