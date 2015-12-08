! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! A module to contain information about submodels

MODULE Submodel_Mod

! Description:
!   Data module to contain information about submodels. Largely
!   unnecessary but used for consistency with UM.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

!  1. Maximum internal model/submodel array sizes for this version.
! Max no. of internal models (1 more for expansion according to csubmmax.h)
  INTEGER, PARAMETER  ::  N_Internal_Model_Max = 2

! Max no. of subm. dump parts (1 more for expansion according to csubmmax.h)
  INTEGER, PARAMETER  ::  N_Submodel_Partition_Max = 1

! Max value of int. model id
  INTEGER, PARAMETER  ::  Internal_Id_Max = N_Internal_Model_Max

! Max value of subm. dump id
  INTEGER, PARAMETER  ::  Submodel_Id_Max = N_Submodel_Partition_Max


!  2. Internal Model identifiers in Long and Short form
  INTEGER, PARAMETER  ::  Atmos_IM = 1  
  INTEGER, PARAMETER  ::  A_IM     = 1

! Atmos SM partition in Long and Short form
  INTEGER, PARAMETER  ::   Atmos_SM = 1 
  INTEGER, PARAMETER  ::   A_SM     = 1

!  3. Lists of internal models and their submodel dump partitions -
!     initialised by the user interface - experiment specific.
  INTEGER, SAVE  :: N_Internal_Model      ! No. of internal models
  INTEGER, SAVE  :: N_Submodel_Partition  ! No. of submodel partitions
  INTEGER, SAVE  :: Internal_Model_List(N_Internal_Model_Max) 

! Submodel identifier for each internal model in list
  INTEGER, SAVE       :: Submodel_For_IM(N_Internal_Model_Max) 

! Submodel number for each submodel id
  INTEGER, SAVE       :: Submodel_For_SM(N_Internal_Model_Max) 

! Namelist for information in 3.
  NAMELIST/NSUBMODL/                                                   &
        N_Internal_Model, N_Submodel_Partition,                        &
        Internal_Model_List, Submodel_For_IM

!  4. Lists calculated in model from user interface supplied arrays -
!     - experiment specific.
! No of internal models in  each submodel partition indexed by sm identifier
  INTEGER, SAVE  :: N_Internal_For_SM(Submodel_ID_Max)

! List of submodel partition identifiers
  INTEGER, SAVE  :: Submodel_Partition_List(N_Submodel_Partition_Max)

! Submodel partition identifier indexed by internal model identifier
  INTEGER, SAVE  :: Submodel_Partition_Index(Internal_ID_Max)

! Sequence number of internal model indexed by internal model identifier
! required to map from id to STASH internal model sequence
  INTEGER, SAVE  :: Internal_Model_Index(Internal_ID_Max)

! Last internal model within a submodel partition if .TRUE.,
! indexed by internal model id.
  LOGICAL, SAVE  :: Last_IM_IN_SM(Internal_ID_Max)

  INTEGER, SAVE  ::  Submodel_Ident

END MODULE Submodel_Mod 
