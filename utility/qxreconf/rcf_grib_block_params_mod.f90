! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Define the 'Grib_record' type and the parameters for the 'blocks'

Module rcf_GRIB_Block_Params_Mod

! Description: Parameters defining the elements of the block data read
!              from GRIB files.
!              Plus the derived types in which the data is stored
!              within reconfiguration
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!

!Pull in some parameters from other modules

Use Rcf_GRIB_Lookups_Mod, Only      :    &
  LenDesc

IMPLICIT NONE

!=======================================================================
!  Block Sizes
!=======================================================================

Integer , Parameter    :: Len_Block_0 =    4   !\ The Blocks used
Integer , Parameter    :: Len_Block_1 = 1000   !/ to store the
Integer , Parameter    :: Len_Block_2 =   20   !\ GRIB header
Integer , Parameter    :: Len_Block_3 =    2   !/ Info
Integer , Parameter    :: Len_Block_4 =    2   !\                      .
Integer , Parameter    :: Len_Block_R =   20   !/                      .

!=======================================================================
!  Variables Usedfor calls to 'DECODE'
!=======================================================================

Integer , Parameter    :: LenArrayMax =5000000 ! Max Array size used
                                               ! for work arrays and
                                               ! such in DECODE

Real                   :: FpData(LenArrayMax)  ! Array decoded field is
                                               ! returned in.
Logical                :: LgData(LenArrayMax)  ! Logical Array for
                                               ! Logical fields (LSM)

!=======================================================================
!  The Derived Type used to hold the field header Information
!=======================================================================

Type Grib_Record

  Integer               :: Start_pos            !\ Start pos in file
  Integer               :: StashCode            !  Stash Code for field
  Integer               :: Data_Type            !  Type of data

  Integer               :: Num_Fp               !\ Number of elements
  Integer               :: Num_Vert             !/ returned by DECODE
  Integer               :: Num_Bitmap           !\ in the respective
  Integer               :: Num_Quasi            !/ arrays.

  Integer               :: Block_0(Len_Block_0) !\ The blocks holding
  Integer               :: Block_1(Len_Block_1) !/ the header info
  Integer               :: Block_2(Len_Block_2) !\ returned by DECODE
  Integer               :: Block_3(Len_Block_3) !/
  Integer               :: Block_4(Len_Block_4) !

  Real                  :: Block_R(Len_Block_R)

  Real,Pointer          :: VertCoords(:)        ! Points to dynamically
                                                ! allocated array

  Character(len=lenDesc) :: Desc                ! Description of param

  Type(Grib_Record), Pointer     :: Next        !\ Pointers to the next
  Type(Grib_Record), Pointer     :: Prev        !/ and prev in list

End Type Grib_Record

Integer , Parameter    :: Grb_Data_Real =  1   !/
Integer , Parameter    :: Grb_Data_Int  =  2   !/
Integer , Parameter    :: Grb_Data_Log  =  3   !/
!=======================================================================
!  The Derived Type to /mark/ the lists of fields
!=======================================================================
! note- if the list has only one member then begin and end will point
! to the same place.

Type List_Marker
  Type(Grib_Record),Pointer      :: Begin
  Type(Grib_Record),Pointer      :: End
  Integer                        :: LstCount
End Type List_Marker

!=======================================================================
!  Block 0
!=======================================================================

Integer, Parameter  :: p_Ed_no       =  1  ! GRIB edition number
Integer, Parameter  :: p_Tbl_Vers_No =  2  ! Tables version no.
Integer, Parameter  :: p_Mes_Len     =  3  ! Total message (record) len
Integer, Parameter  :: p_B0_undef    =  4  ! currently undefined

!=======================================================================
!  Block 1
!=======================================================================

Integer, Parameter  :: p_Orig_cntr   =  1  ! Originating center
Integer, Parameter  :: p_Gen_Proc_ID =  2  ! Generating process ID
Integer, Parameter  :: p_Grd_ID      =  3  ! Grid ID no.
Integer, Parameter  :: p_Blck_ID     =  4  ! Block ID flags (table 1)
Integer, Parameter  :: p_Param_ID    =  5  ! Parameter ID (Table 2)
Integer, Parameter  :: p_Lvl_Type    =  6  ! Type of Level ID (table 3)
Integer, Parameter  :: p_Lvl_Desc_1  =  7  ! 1st level desription param
Integer, Parameter  :: p_Lvl_Desc_2  =  8  ! 2nd level desription param
                                           ! either are 0 if not needed
                                           ! otherwise specify hieght,
                                           ! pressure etc
Integer, Parameter  :: p_Ref_Year    =  9  ! Year of century, ref time
Integer, Parameter  :: p_Ref_Month   = 10  ! Month of ref time
Integer, Parameter  :: p_Ref_Day     = 11  ! Day of Reference Time
Integer, Parameter  :: p_Ref_Hour    = 12  !
Integer, Parameter  :: p_Ref_Min     = 13  !
Integer, Parameter  :: p_Time_Unit   = 14  ! Indicator of time unit
                                           ! (Table 4)
Integer, Parameter  :: p_Time_Int_1  = 15  ! 1st time interval if reqd
Integer, Parameter  :: p_Time_Int_2  = 16  ! 2nd time interval (P2)
Integer, Parameter  :: p_Time_Range  = 17  ! Time range indicator
                                           ! Table 5
Integer, Parameter  :: p_No_in_Avg   = 18  ! No. included in Average
                                           ! if reqd for time range
Integer, Parameter  :: p_Ref_Cent    = 19  ! Century of ref time
Integer, Parameter  :: p_Dec_Scale   = 20  ! Decimal Scale factor
Integer, Parameter  :: p_Len_Blck_1  = 21  ! Length of block 1 in octets

! elements 22-33 are reserved for future use

Integer, Parameter  :: p_Sub_C_ID    = 34  ! Sub center ID no.

! octets 36 onwards contain originating center octets.
! 1 octet is returned in each array element. The number of originating
! center octets is Length of block 1 - 40
! If center = 74 and Sub center Id is 2 then 36-79 hold superstash pp
! data in 32 bit words.

!=======================================================================
!  Block 2
!=======================================================================

! The 1st 3 elements of block 2 are always the same.
! The rest vary depending on grid type.
! At present, only lat/long, Gaussian, mercator grid params
! are specified here. i.e polar stereographic are not mentioned

Integer, Parameter  :: p_V_Coord_Ps =  1  ! no. Vertical coord params
Integer, Parameter  :: p_PV_PL      =  2  ! PV or PL if level is hybrid
                                          ! or quasi regular
Integer, Parameter  :: p_Data_Rep   =  3  ! Data representation type

! Assuming Lat/long, Gaussian, Mercator or Rotated Lat/long

Integer, Parameter  :: p_Pnts_Prll  =  4  ! No. of points along parall
Integer, Parameter  :: p_Pnts_Merid =  5  ! No. of points along merid
Integer, Parameter  :: p_LatGrdPnt1 =  6  ! Lat of 1st Grid point
Integer, Parameter  :: p_LonGrdPnt1 =  7  ! Long of 1st Grid point
Integer, Parameter  :: p_Res_Flg    =  8  ! resolution and cmpnt flags
                                          ! as defined in Table 7
Integer, Parameter  :: p_LatExtrmPt =  9  ! Lat of extreme point
Integer, Parameter  :: p_LonExtrmPt = 10  ! Long of extreme point
Integer, Parameter  :: p_IncrPrll   = 11  ! Increment along Parrallel
Integer, Parameter  :: p_IncrMerid  = 12  ! Increment along Meridian
Integer, Parameter  :: p_ScanModeF  = 13  ! Scanning mode flags as
                                          ! defined in Table 8
Integer, Parameter  :: p_S_Pole_Lat = 14  ! Lat of S Pole -Iff rotated
Integer, Parameter  :: p_S_Pole_Lon = 15  ! Long of S Pole -Iff rotated

! Elements 16 to 20 are undefined

!=======================================================================
!  Block 3
!=======================================================================

! Block 3 is currently not held in the GRIB_Record derived data type
! and is not defined in/by DECODE

!=======================================================================
!  Block 4
!=======================================================================

Integer, Parameter  :: p_PackFlag  =  1  ! Packing Flag
Integer, Parameter  :: p_B4Undef   =  2  ! Input value, unchanged

!=======================================================================
!  Block R
!=======================================================================

Real,    Parameter  :: p_RotAngl   =  1  ! Angle of rotation for rotated
                                         ! Lat long grid
Real,    Parameter  :: p_GMDI      =  2  ! Missing data Indicator used
                                         ! for bitmap grids

! Elements 3 to 20 are currently undefined

!=======================================================================
!  Miscellaneous values set up as parameters
!=======================================================================

! Values to do with dump header addressing
Integer, Parameter  :: p_Len_IntC    = 46 ! normally from start dump
Integer, Parameter  :: p_Len_RealC   = 38 ! ditto -or- Recon namelist
Integer, Parameter  :: p_Len2_LevDepC = 8 ! ditto -or- Recon namelist
Integer, Parameter  :: p_Len1_RowDepC = 0 ! ditto -or- Recon namelist
Integer, Parameter  :: p_Len2_RowDepC = 0 ! ditto -or- Recon namelist
Integer, Parameter  :: p_Len1_ColDepC = 0 ! ditto -or- Recon namelist
Integer, Parameter  :: p_Len2_ColDepC = 0 ! ditto -or- Recon namelist
Integer, Parameter  :: p_Len1_FldsofC = 0 ! ditto -or- Recon namelist
Integer, Parameter  :: p_Len2_FldsofC = 0 ! ditto -or- Recon namelist
Integer, Parameter  :: p_Len_ExtraC   = 0 ! ditto -or- Recon namelist
Integer, Parameter  :: p_Len_HistFile = 0 ! ditto -or- Recon namelist
Integer, Parameter  :: p_Len_CompFldI1= 0 ! ditto -or- Recon namelist
Integer, Parameter  :: p_Len_CompFldI2= 0 ! ditto -or- Recon namelist
Integer, Parameter  :: p_Len_CompFldI3= 0 ! ditto -or- Recon namelist
Integer, Parameter  :: p_Len1_Lookup = 64 ! ditto -or- Recon namelist

!=======================================================================
!  Values Used in GRIB Table 2 (Parameter IDs)
!=======================================================================

!Those used by ECMWF (the non-standard ones)
Integer, Parameter  :: EID_Geopotential   = 129
Integer, Parameter  :: EID_Temperature    = 130
Integer, Parameter  :: EID_U_Wind         = 131
Integer, Parameter  :: EID_V_Wind         = 132
Integer, Parameter  :: EID_Spec_Humidity  = 133
Integer, Parameter  :: EID_Surf_Press     = 134
Integer, Parameter  :: EID_W_Wind         = 135
Integer, Parameter  :: EID_Surf_Temp      = 139
Integer, Parameter  :: EID_Soil_Moisture  = 140
Integer, Parameter  :: EID_Log_Surf_Press = 152

!=======================================================================
!  Values Used in GRIB Table 3 (Type and Value of level)
!=======================================================================

Integer, Parameter  :: Tb3_Surface  = 1   ! Ground or Water Surface
Integer, Parameter  :: Tb3_Pressure = 100 ! Isobaric Level type
Integer, Parameter  :: Tb3_Mid_IsoB = 101 ! layer between Isobaric
Integer, Parameter  :: Tb3_Mean_Sea = 102 ! Mean Sea Level
Integer, Parameter  :: Tb3_Hybrid   = 109 ! Hybrid Level
Integer, Parameter  :: Tb3_GrndLay  = 112 ! Below ground layer

!=======================================================================
!  Done
!=======================================================================

End Module rcf_GRIB_Block_Params_Mod
