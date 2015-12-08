! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Magic numbers for header components

MODULE Rcf_HeadAddress_Mod

! Description:
!   Magic numbers for Fixed, Integer, Real, etc header components.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

! Fixed length header constants

INTEGER,PARAMETER :: FH_Version               = 1
INTEGER,PARAMETER :: FH_SubModel              = 2
INTEGER,PARAMETER :: FH_VertCoord             = 3
INTEGER,PARAMETER :: FH_HorizGrid             = 4
INTEGER,PARAMETER :: FH_Dataset               = 5
INTEGER,PARAMETER :: FH_RunId                 = 6
INTEGER,PARAMETER :: FH_ExptNo                = 7
INTEGER,PARAMETER :: FH_CalendarType          = 8
INTEGER,PARAMETER :: FH_GridStagger           = 9
INTEGER,PARAMETER :: FH_AncilDataId           = 10
INTEGER,PARAMETER :: FH_ProjNo                = 11
INTEGER,PARAMETER :: FH_OcnBndyConds          = 11     ! ocean boundary condition
INTEGER,PARAMETER :: FH_ModelVersion          = 12
INTEGER,PARAMETER :: FH_OceanDynamics         = 12
INTEGER,PARAMETER :: FH_ObsFileType           = 14

INTEGER,PARAMETER :: FH_DTYear                = 21
INTEGER,PARAMETER :: FH_DTMonth               = 22
INTEGER,PARAMETER :: FH_DTDay                 = 23
INTEGER,PARAMETER :: FH_DTHour                = 24
INTEGER,PARAMETER :: FH_DTMinute              = 25
INTEGER,PARAMETER :: FH_DTSecond              = 26
INTEGER,PARAMETER :: FH_DTDayNo               = 27

INTEGER,PARAMETER :: FH_VTYear                = 28
INTEGER,PARAMETER :: FH_VTMonth               = 29
INTEGER,PARAMETER :: FH_VTDay                 = 30
INTEGER,PARAMETER :: FH_VTHour                = 31
INTEGER,PARAMETER :: FH_VTMinute              = 32
INTEGER,PARAMETER :: FH_VTSecond              = 33
INTEGER,PARAMETER :: FH_VTDayNo               = 34

INTEGER,PARAMETER :: FH_CTYear                = 35
INTEGER,PARAMETER :: FH_CTMonth               = 36
INTEGER,PARAMETER :: FH_CTDay                 = 37
INTEGER,PARAMETER :: FH_CTHour                = 38
INTEGER,PARAMETER :: FH_CTMinute              = 39
INTEGER,PARAMETER :: FH_CTSecond              = 40
INTEGER,PARAMETER :: FH_CTDayNo               = 41

INTEGER,PARAMETER :: FH_IntCStart             = 100
INTEGER,PARAMETER :: FH_IntCSize              = 101
INTEGER,PARAMETER :: FH_RealCStart            = 105
INTEGER,PARAMETER :: FH_RealCSize             = 106
INTEGER,PARAMETER :: FH_LevDepCStart          = 110
INTEGER,PARAMETER :: FH_LevDepCSize1          = 111
INTEGER,PARAMETER :: FH_LevDepCSize2          = 112
INTEGER,PARAMETER :: FH_RowDepCStart          = 115
INTEGER,PARAMETER :: FH_RowDepCSize1          = 116
INTEGER,PARAMETER :: FH_RowDepCSize2          = 117
INTEGER,PARAMETER :: FH_ColDepCStart          = 120
INTEGER,PARAMETER :: FH_ColDepCSize1          = 121
INTEGER,PARAMETER :: FH_ColDepCSize2          = 122
INTEGER,PARAMETER :: FH_FldsOfCStart          = 125
INTEGER,PARAMETER :: FH_FldsOfCSize1          = 126
INTEGER,PARAMETER :: FH_FldsOfCSize2          = 127
INTEGER,PARAMETER :: FH_ExtraCStart           = 130
INTEGER,PARAMETER :: FH_ExtraCSize            = 131
INTEGER,PARAMETER :: FH_HistStart             = 135
INTEGER,PARAMETER :: FH_HistSize              = 136
INTEGER,PARAMETER :: FH_CompFldI1Start        = 140
INTEGER,PARAMETER :: FH_CompFldI1Size         = 141
INTEGER,PARAMETER :: FH_CompFldI2Start        = 142
INTEGER,PARAMETER :: FH_CompFldI2Size         = 143
INTEGER,PARAMETER :: FH_CompFldI3Start        = 144
INTEGER,PARAMETER :: FH_CompFldI3Size         = 145
INTEGER,PARAMETER :: FH_LookupStart           = 150
INTEGER,PARAMETER :: FH_LookupSize1           = 151
INTEGER,PARAMETER :: FH_LookupSize2           = 152
INTEGER,PARAMETER :: FH_NumProgFields         = 153
INTEGER,PARAMETER :: FH_DataStart             = 160
INTEGER,PARAMETER :: FH_DataSize              = 161
INTEGER,PARAMETER :: FH_MaxDataSize           = 162

INTEGER,PARAMETER :: FH_Version_Value         = 20
INTEGER,PARAMETER :: FH_SubModel_Atmos        = 1
INTEGER,PARAMETER :: FH_SubModel_Ocean        = 2
INTEGER,PARAMETER :: FH_SubModel_Wave         = 4

INTEGER,PARAMETER :: FH_VertCoord_Hybrid      = 1
INTEGER,PARAMETER :: FH_VertCoord_Sigma       = 2
INTEGER,PARAMETER :: FH_VertCoord_Pressure    = 3
INTEGER,PARAMETER :: FH_VertCoord_Depth       = 4
INTEGER,PARAMETER :: FH_VertCoord_CP          = 5
INTEGER,PARAMETER :: FH_VertCoord_Wave        = 6

INTEGER,PARAMETER :: FH_HorizGrid_Global      = 0
INTEGER,PARAMETER :: FH_HorizGrid_NH          = 1
INTEGER,PARAMETER :: FH_HorizGrid_SH          = 2
INTEGER,PARAMETER :: FH_HorizGrid_LamNoWrap   = 3
INTEGER,PARAMETER :: FH_HorizGrid_LamWrap     = 4
INTEGER,PARAMETER :: FH_HorizGrid_Eq          = 100
INTEGER,PARAMETER :: FH_HorizGrid_LamNoWrapEq = 103
INTEGER,PARAMETER :: FH_HorizGrid_LamWrapEq   = 104

INTEGER,PARAMETER :: FH_GridStagger_A         = 1
INTEGER,PARAMETER :: FH_GridStagger_B         = 2
INTEGER,PARAMETER :: FH_GridStagger_C         = 3
INTEGER,PARAMETER :: FH_GridStagger_D         = 4
INTEGER,PARAMETER :: FH_GridStagger_E         = 5
INTEGER,PARAMETER :: FH_GridStagger_Endgame   = 6

INTEGER,PARAMETER :: FH_Dataset_InstDump      = 1
INTEGER,PARAMETER :: FH_Dataset_MeanDump      = 2
INTEGER,PARAMETER :: FH_Dataset_FF            = 3
INTEGER,PARAMETER :: FH_Dataset_Ancil         = 4
INTEGER,PARAMETER :: FH_Dataset_Boundary      = 5
INTEGER,PARAMETER :: FH_Dataset_ACOBS         = 6
INTEGER,PARAMETER :: FH_Dataset_VAROBS        = 7
INTEGER,PARAMETER :: FH_Dataset_CX            = 8
INTEGER,PARAMETER :: FH_Dataset_COV           = 9
INTEGER,PARAMETER :: FH_Dataset_OBSTORE       = 10

INTEGER,PARAMETER :: FH_ObsFileType_Atmos     = 1
INTEGER,PARAMETER :: FH_ObsFileType_Ocean     = 2
INTEGER,PARAMETER :: FH_ObsFileType_SST       = 3
INTEGER,PARAMETER :: FH_ObsFileType_Wave      = 4

INTEGER, PARAMETER :: IC_TorTheta             = 1 !location in header
INTEGER, PARAMETER :: IC_TorTheta_T           = 1 !value of above if T
INTEGER, PARAMETER :: IC_TorTheta_Theta       = 2 !value of above if Theta
INTEGER, PARAMETER :: IC_XLen                 = 6
INTEGER, PARAMETER :: IC_YLen                 = 7
INTEGER, PARAMETER :: IC_PLevels              = 8
INTEGER, PARAMETER :: IC_WetLevels            = 9
INTEGER, PARAMETER :: IC_SoilTLevels          = 10
INTEGER, PARAMETER :: IC_NoCloudLevels        = 11 ! ATMOS only
INTEGER, PARAMETER :: IC_NoSeaPts             = 11 ! OCEAN only
INTEGER, PARAMETER :: IC_TracerLevs           = 12
INTEGER, PARAMETER :: IC_BLevels              = 13
INTEGER, PARAMETER :: IC_TracerVars           = 14
INTEGER, PARAMETER :: IC_HeightMethod         = 17 !method for creating heights
INTEGER, PARAMETER :: IC_RiverRowLength       = 19 !river-routing row-length
INTEGER, PARAMETER :: IC_RiverRows            = 20 !river-routing rows
INTEGER, PARAMETER :: IC_MDI                  = 21
INTEGER, PARAMETER :: IC_1stConstRho          = 24
INTEGER, PARAMETER :: IC_NumLandPoints        = 25
INTEGER, PARAMETER :: IC_NumOzoneLevs         = 26
INTEGER, PARAMETER :: IC_SoilMoistLevs        = 28
INTEGER, PARAMETER :: IC_NumObsTotal          = 28
INTEGER, PARAMETER :: IC_LenObCol             = 29
INTEGER, PARAMETER :: IC_LenCxCol             = 30 ! Varobs, not acobs
INTEGER, PARAMETER :: IC_ObsGroup             = 31 ! "
INTEGER, PARAMETER :: IC_ObsRelease           = 32 ! "
INTEGER, PARAMETER :: IC_NumMetaMax           = 33 ! "
INTEGER, PARAMETER :: IC_ConvectLevs          = 34
INTEGER, PARAMETER :: IC_NumItemMax           = 34 ! "
INTEGER, PARAMETER :: IC_NumObVarMax          = 35
INTEGER, PARAMETER :: IC_NumObPItemMax        = 36
INTEGER, PARAMETER :: IC_NumCxPItemMax        = 37
INTEGER, PARAMETER :: IC_NumCxSFVarMax        = 38
INTEGER, PARAMETER :: IC_NumCxUaVarMax        = 39
INTEGER, PARAMETER :: IC_NumMeta              = 40
INTEGER, PARAMETER :: IC_NumItem              = 41
INTEGER, PARAMETER :: IC_NumObVar             = 42
INTEGER, PARAMETER :: IC_NumObPItem           = 43
INTEGER, PARAMETER :: IC_NumCxPItem           = 44
INTEGER, PARAMETER :: IC_NumCxSfVar           = 45
INTEGER, PARAMETER :: IC_NumCxUaVar           = 46
INTEGER, PARAMETER :: IC_NumObLev             = 47
INTEGER, PARAMETER :: IC_NumCxLev             = 48
INTEGER, PARAMETER :: IC_NumVarBatches        = 49

INTEGER, PARAMETER :: RC_LongSpacing          = 1
INTEGER, PARAMETER :: RC_LatSpacing           = 2
INTEGER, PARAMETER :: RC_FirstLat             = 3
INTEGER, PARAMETER :: RC_FirstLong            = 4
INTEGER, PARAMETER :: RC_PoleLat              = 5
INTEGER, PARAMETER :: RC_PoleLong             = 6
INTEGER, PARAMETER :: RC_SWLDEG               = 7 ! Ocean - lat of South wall
INTEGER, PARAMETER :: RC_WEDGEDEG             = 8 !   " = long of West bdy
INTEGER, PARAMETER :: RC_ModelTop             = 16
INTEGER, PARAMETER :: RC_PressureTop          = 17
INTEGER, PARAMETER :: RC_AtmMoist             = 18
INTEGER, PARAMETER :: RC_AtmMass              = 19
INTEGER, PARAMETER :: RC_AtmEnergy            = 20
INTEGER, PARAMETER :: RC_EnergyCorr           = 21

INTEGER, PARAMETER :: CC_Meta_Latitude        = 1 ! Used in varobs
INTEGER, PARAMETER :: CC_Meta_Longitude       = 2 !      "
INTEGER, PARAMETER :: CC_Meta_Time            = 3 !      "
INTEGER, PARAMETER :: CC_Meta_Type            = 4 !      "
INTEGER, PARAMETER :: CC_Meta_Call            = 5 !      "
INTEGER, PARAMETER :: CC_Meta_Level           = 6 !      "
INTEGER, PARAMETER :: CC_Meta_RepPGE          = 7 !      "

INTEGER, PARAMETER :: CC_Item_Value           = 1 ! Used in varobs
INTEGER, PARAMETER :: CC_Item_Error           = 2 !      "
INTEGER, PARAMETER :: CC_Item_PGE             = 3 !      "
INTEGER, PARAMETER :: LDC_EtaTheta  = 1
INTEGER, PARAMETER :: LDC_Pressure  = 1
INTEGER, PARAMETER :: LDC_MLIndex   = 1
INTEGER, PARAMETER :: LDC_EtaRho    = 2
INTEGER, PARAMETER :: LDC_RHCrit    = 3
INTEGER, PARAMETER :: SoilDepths    = 4
INTEGER, PARAMETER :: LDC_ZseaTheta = 5
INTEGER, PARAMETER :: LDC_ak_hybrid = 5
INTEGER, PARAMETER :: LDC_CkTheta   = 6
INTEGER, PARAMETER :: LDC_bk_hybrid = 6
INTEGER, PARAMETER :: LDC_ZseaRho   = 7
INTEGER, PARAMETER :: LDC_CkRho     = 8

INTEGER, PARAMETER :: RDC_Phi_input_p         = 1 ! New Row and Column
INTEGER, PARAMETER :: RDC_Phi_input_v         = 2 ! Dependent Constants

INTEGER, PARAMETER :: CDC_Lambda_input_p      = 1
INTEGER, PARAMETER :: CDC_Lambda_input_u      = 2
END MODULE Rcf_HeadAddress_Mod
