! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module: REGRID_TYPES -------------------------------------------------
!
!  Purpose: Defines data types used for parallel regridding routines 
!  all x,y grid point information are store in global coordinates 

!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: MPP

MODULE REGRID_TYPES

 
  IMPLICIT NONE   
  
  ! to store regridding information needed to regrid 
  ! using area averaging method
  ! see regrid_methods.F90
  TYPE CONCERN
    
    INTEGER :: proc_num = -1
    ! the process to receive or send resource associate to this instance 
    ! of concern 
    
    INTEGER, POINTER :: x(:) => NULL() , y(:) => NULL()
    ! grid points in i.e. field(i) -> (x(i), y(i)) 
    
    REAL, POINTER :: field(:) => NULL()
    ! value at grid points for concern
    
    INTEGER :: size = 0
    ! size of field and grid points
    
  END TYPE CONCERN
  
  ! to store regridding information needed to regrid 
  ! using 'map max weight src to target' method
  ! see regrid_methods.F90
  TYPE CONCERN_MAX
    
    INTEGER :: proc_num = -1
    
    ! grid points in global indexing
    INTEGER, POINTER :: x(:) => NULL() , y(:) => NULL()
    ! source field grid points
    
    INTEGER, POINTER :: xtarg(:) => NULL(), ytarg(:) => NULL()
    ! target field grid points
    
    REAL, POINTER :: weight(:) => NULL()
    ! weight of source field point
    
    LOGICAL, POINTER :: contribute(:) => NULL()
    ! whether given source grid point contributes to target grid point
    
    REAL, POINTER :: field(:) => NULL()
    ! source field value at grid point
    
    INTEGER :: size = 0
    ! size data above, for convenience
    
  END TYPE CONCERN_MAX
  
  ! stores regridding information for a single 
  ! target grid point. The information represents
  ! all source points needed to regrid it and 
  ! respective weights 
  TYPE CONTRIBUTION_INFO
    
    INTEGER size 
    ! number on contributing targets to source subdomain
    
    INTEGER, POINTER :: x(:) => NULL()
    INTEGER, POINTER :: y(:) => NULL()
    ! src grid point
    
    INTEGER, POINTER :: contrib_proc(:) => NULL()
    ! what processor is contributing
    
    REAL, POINTER :: weight(:) => NULL()
    ! src point weight 
    
  END TYPE CONTRIBUTION_INFO
  
END MODULE REGRID_TYPES
