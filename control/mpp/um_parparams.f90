! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP.

MODULE UM_ParParams
  USE Field_Types
  USE domain_params

  IMPLICIT NONE

  !   Description:
  !   Two sets of parameters are set up -
  !     i)  for the MPP-UM itself.
  !     ii) for the interface to the Message Passing Software.
  !

  !=================================================================
  ! Parameters needed for the MPP-UM
  !=================================================================

  ! maximum number of spatial dimensions
  INTEGER,PARAMETER :: Ndim_max           = 3 ! 3d data
  
  ! number of different halo types
  INTEGER,PARAMETER :: NHalo_max          = 3 ! for N.D. atmos. model
  
  INTEGER,PARAMETER :: halo_type_single   = 1
  INTEGER,PARAMETER :: halo_type_extended = 2
  INTEGER,PARAMETER :: halo_type_no_halo  = 3
  
  ! Used in addressing to indicate if calculation is for a local or
  ! global (ie. disk dump) size  
  INTEGER,PARAMETER :: local_data         = 1
  INTEGER,PARAMETER :: global_dump_data   = 2
    
  !=================================================================
  ! Parameters needed for the Message Passing Software
  !=================================================================
  
  ! Processor addresses in the neighbour array
  INTEGER,PARAMETER :: NoDomain    = -1
  INTEGER,PARAMETER :: PNorth      = 1
  INTEGER,PARAMETER :: PEast       = 2
  INTEGER,PARAMETER :: PSouth      = 3
  INTEGER,PARAMETER :: PWest       = 4
    
  INTEGER,PARAMETER :: BC_STATIC   = 1 ! Static boundary conditions
  INTEGER,PARAMETER :: BC_CYCLIC   = 2 ! Cyclic boundary conditions
  INTEGER,PARAMETER :: BC_OVERPOLE = 3 ! Transfer over pole

END MODULE UM_ParParams
