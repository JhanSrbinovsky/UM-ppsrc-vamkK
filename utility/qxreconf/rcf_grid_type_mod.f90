! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Module to define a grid data type.

Module Rcf_Grid_Type_Mod

! Description:
!   Defines the grid_type data type
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Implicit None

Type grid_type

  ! General attributes
  Logical      :: global       ! global == .true., LAM ==.false.
  Logical      :: rotated      ! .true. if rotated LAM
  Integer      :: grid_stagger ! As read from fixed header - same values.
  

  ! Grid dimensions - horizontally
  ! 1st the Global variables (in MPP sense)
  Integer      :: glob_p_row_length   ! length of 'pressure' rows
  Integer      :: glob_u_row_length   ! length of U rows
  Integer      :: glob_v_row_length   ! length of V rows
  Integer      :: glob_r_row_length   ! length of River rows
  Integer      :: glob_p_rows         ! No. of 'pressure' rows
  Integer      :: glob_u_rows         ! No. of U rows
  Integer      :: glob_v_rows         ! No. of V rows
  Integer      :: glob_r_rows         ! No. of River rows
  Integer      :: glob_p_field        ! Total size of P field
  Integer      :: glob_u_field        ! Total size of U field
  Integer      :: glob_v_field        ! Total size of V field
  Integer      :: glob_r_field        ! Total size of River field
  Integer      :: glob_land_field     ! No. of land points in field

 ! Now the Local variables (in MPP sense)
  Integer      :: loc_p_row_length   ! length of 'pressure' rows
  Integer      :: loc_u_row_length   ! length of U rows
  Integer      :: loc_v_row_length   ! length of V rows
  Integer      :: loc_r_row_length   ! length of River rows
  Integer      :: loc_p_rows         ! No. of 'pressure' rows
  Integer      :: loc_u_rows         ! No. of U rows
  Integer      :: loc_v_rows         ! No. of V rows
  Integer      :: loc_r_rows         ! No. of River rows
  Integer      :: loc_p_field        ! Total size of P field
  Integer      :: loc_u_field        ! Total size of U field
  Integer      :: loc_v_field        ! Total size of V field
  Integer      :: loc_r_field        ! Total size of River field
  Integer      :: loc_land_field     ! No of land points in field

  ! Grid dimensions - vertically
  ! These are *all* Global...
  Integer      :: model_levels ! No. of pressure levels
  Integer      :: wet_levels   ! No. of moist levels
  Integer      :: cloud_levels ! No. of cloud levels
  Integer      :: st_levels    ! No. of soil temp. levels
  Integer      :: sm_levels    ! No. of soil moisture levels
  Integer      :: bl_levels    ! No. of boundary-layer levels
  Integer      :: ozone_levels ! No. of ozone levels
  Integer      :: tr_levels    ! No. of tracer levels
  Integer      :: conv_levels  ! No. of convective levels

  ! Stuff for calculation of heights
  Integer       :: height_gen_method   ! method used to generate heights
  Integer       :: first_constant_r_rho_level
  Real          :: z_top_of_model
  Real, Pointer :: eta_theta_levels( : )
  Real, Pointer :: eta_rho_levels( : )
  Real, Pointer :: rhcrit( : )
  Real, Pointer :: soil_depths( : )

  ! VarRes horizontal grid spacing
  Real, Pointer ::  Lambda_p( : )
  Real, Pointer ::  Lambda_u( : )
  Real, Pointer ::  Phi_p( : )
  Real, Pointer ::  Phi_v( : )

  ! Rotated pole information
  REAL :: lambda_npole
  REAL :: phi_npole

End Type grid_type

Type (grid_type), Save :: Input_grid
Type (grid_type), Save :: Output_grid
End Module Rcf_Grid_Type_Mod
