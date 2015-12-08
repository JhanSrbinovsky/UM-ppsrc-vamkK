! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Land-Sea mask data

Module Rcf_Lsm_Mod

! Description:
!   Stores land sea masks and related sizes/data & coastal adjustment
!   gather indexes etc.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   03/10/01   Addition of Implicit None. R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Implicit None

! Storage for land-sea masks
Logical, Allocatable, Target, Save   :: glob_lsm_in(:)
Logical, Allocatable, Target, Save   :: local_lsm_in(:)
Logical, Allocatable, Target, Save   :: glob_lsm_out(:)
Logical, Allocatable, Target, Save   :: local_lsm_out(:)
Integer, Target, Save                :: local_land_in
Integer, Target, Save                :: local_land_out
Integer, Target, Save                :: glob_land_in
Integer, Target, Save                :: glob_land_out

! Pointers to active lsm fields
Logical, Pointer, Save       :: glob_atmos_landmask(:)
Logical, Pointer, Save       :: local_atmos_landmask(:)
Integer, Pointer, Save       :: local_land_field
Integer, Pointer, Save       :: glob_land_field

! Storage for gather indexes etc for coastal adjustment
Integer, Save                :: N_Coastal_Points
Integer, Save                :: N_land_points_unres
Integer, Save                :: N_sea_points_unres
Integer, Allocatable, Save   :: Coast_Index_In(:)
Integer, Allocatable, Save   :: Coast_Index_Out(:)
Integer, Allocatable, Save   :: Index_Targ_Land(:)
Integer, Allocatable, Save   :: Index_Targ_Sea(:)
Integer, Allocatable, Save   :: Land_Unres_Index(:)
Integer, Allocatable, Save   :: Sea_Unres_Index(:)

! Logical for coast adjustment
Logical, Save                :: Cyclic

End Module Rcf_Lsm_Mod
