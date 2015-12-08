! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Performs horizontal interpolation and related tasks

Module Rcf_horizontal_Mod
IMPLICIT NONE

!  Subroutine Rcf_horizontal - horizontal interpolation
!
! Description:
!   This module contains a wrapper subroutine for horizontal
!   interpolation. The interpolation is done level-by-level
!   but on a vertically decomposed grid. Thus we need to gather
!   data by levels onto the compute PEs and then rescatter it onto
!   the output grids. This won't be the fastest routine in the world,
!   but should be reliable and be easily verifiable for
!   bit-reproducibility etc.
!
! Method:
!   Note that land compressed fields are uncompressed and recompressed,
!   copying is done if that is all that is required, coastal
!   adjustment is performed and polar row averaging is done.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_horizontal( field_in, field_out, grid_in, grid_out )

Use Ereport_mod, Only : &
    Ereport

USE mask_compression, ONLY: compress_to_mask, expand_from_mask

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_Field_Type_Mod, Only : &
    field_type

USE UM_ParVars, Only :  &
    nproc,              &
    mype,               &
    gc_all_proc_group,  &
    change_decomposition

USE decomp_params, ONLY : &
    decomp_rcf_input,     &
    decomp_rcf_output

Use Rcf_select_weights_mod, Only : &
    Rcf_select_weights

Use Rcf_average_polar_mod, Only : &
    Rcf_average_polar

Use Rcf_H_Int_Ctl_Mod, Only : &
    Rcf_H_Int_Ctl

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal,        &
    LTimer

Use Rcf_Lsm_Mod, Only : &
    n_coastal_points,                n_land_points_unres, &
    n_sea_points_unres,              coast_index_in,      &
    coast_index_out,                 index_targ_land,     &
    index_targ_sea,                  land_unres_index,    &
    sea_unres_index,                 local_lsm_in,        &
    local_lsm_out,                   glob_lsm_out,        &
    cyclic

Use Rcf_Recon_Mod, Only : &
    Lspiral_S

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_h_only,                   interp_v_only,       &
    interp_all,                      interp_copy,         &
    interp_no_op

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_prog_sec,    &
    stashcode_tstar,       &
    stashcode_vol_smc_wilt,&
    stashcode_vol_smc_cri, &
    stashcode_vol_smc_sat

Use Rcf_Scatter_Zonal_Field_Mod, Only : &
    Rcf_Scatter_Zonal_Field

Use Rcf_Scatter_Field_Mod, Only : &
    Rcf_Scatter_Field

Use Rcf_Gather_Field_Mod, Only : &
    Rcf_Gather_Field

Use Rcf_Gather_Zonal_Field_Mod, Only : &
    Rcf_Gather_Zonal_Field

Use Rcf_Interp_Weights_Mod  ! All of it

Use Rcf_HeadAddress_Mod, Only : &
    FH_GridStagger

Use cppxref_mod, ONLY :                     &
    ppx_type_real,      ppx_type_int,       &
    ppx_type_log,       ppx_atm_compressed, &
    ppx_atm_ozone,      ppx_atm_tall,       &
    ppx_atm_tzonal,     ppx_atm_tsea

IMPLICIT NONE

! Derived Type Arguments
Type (field_type), Target, Intent(InOut) :: field_in  ! Input data field
Type (field_type), Target, Intent(InOut) :: field_out !Output data field
Type (grid_type), Intent(In)     :: grid_in   ! Input grid sizes
Type (grid_type), Intent(In)     :: grid_out  ! Output grid sizes

! Local data
Integer             :: div         ! for calculation of partition
Integer             :: rem         ! for calculation of partition
Integer             :: pe          ! PE from which to gather/scatter
Integer             :: i,j,k       ! Looping
Integer             :: size        ! check value for lsm
Integer             :: stat        ! GCOM status
Integer             :: field_averaged ! 1/0 status of field
Integer             :: orig_h_int_method ! Stores original interp method
                                         ! since soil can use different one.
Integer             :: sea_points_unres_tmp  ! tmp stuff for coast aj
Integer             :: land_points_unres_tmp ! tmp stuff for coast aj
Integer             :: lsm_tmp( grid_out % glob_p_rows * &
                                grid_out % glob_p_row_length )
Integer             :: Index_targ_land_tmp( grid_out % glob_p_rows * &
                                          grid_out % glob_p_row_length )
Integer             :: Index_targ_sea_tmp( grid_out % glob_p_rows * &
                                          grid_out % glob_p_row_length )
Integer             :: maxdim      ! a size parameter for coast aj
Integer             :: ErrorStatus
Logical             :: averaged    ! return from polar_average
Character (Len=*), Parameter :: RoutineName = 'Rcf_horzontal'
Character (Len=80)           :: Cmessage

                                       ! Single level after scatter
Real, Allocatable   :: level_field_in( : )

                                        ! Single level after interp.
Real, Allocatable   :: level_field_out( : )

Type (field_type), Target     :: field_in_tmp
Type (field_type), Target     :: field_out_tmp
Type (field_type), Pointer    :: ptr_field_in
Type (field_type), Pointer    :: ptr_field_out

! Pointers for choice of data etc
Integer, Pointer              :: ptr_bl_index_b_l(:)
Integer, Pointer              :: ptr_bl_index_b_r(:)
Integer, Pointer              :: ptr_aw_index_targ_lhs(:)
Integer, Pointer              :: ptr_aw_index_targ_top(:)
Real, Pointer                 :: ptr_aw_colat_t(:)
Real, Pointer                 :: ptr_aw_long_l(:)
Real, Pointer                 :: ptr_weight_b_l(:)
Real, Pointer                 :: ptr_weight_b_r(:)
Real, Pointer                 :: ptr_weight_t_l(:)
Real, Pointer                 :: ptr_weight_t_r(:)

External Intf_Coast_AJ, Timer

averaged = .false.

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3)

! If input and output datatypes are different, we have a problem.
If (  field_in % stashmaster % data_type /=                     &
     field_out % stashmaster % data_type ) Then
  Cmessage = 'Input and Output field datatypes differ!'
  ErrorStatus = 10
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Make sure that field_averaged is initialised
field_averaged = 0

!-------------------------------------------------------------------
! Is interpolation activated? If not, copy data across and that's
! all we will do.
!-------------------------------------------------------------------
Select Case( field_in % interp )
  Case( interp_copy, interp_v_only )     ! copy data

  ! Sizes should be the same, but will check...
  If ( field_in % level_size /= field_out % level_size ) Then
    write (6,*) "Aborting due to mismatch in local datasizes "
    write (6,*) "Input dump data size =",field_in % level_size
    write (6,*) "Output dump data size =",field_out % level_size
    Cmessage = 'No interpolation required, but input and output '//&
        'data fields are different sizes'
    ErrorStatus = 20
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  Select Case( field_in % stashmaster % data_type )
    Case ( ppx_type_real )
      field_out % Data(:,:) = field_in % Data(:,:)

    Case ( ppx_type_int )
      If ( Associated( field_in %  Data ) ) Then
        field_out % Data(:,:) = field_in % Data(:,:)
      Else
        field_out % Data_Int(:,:) = field_in % Data_Int(:,:)
      End If

    Case ( ppx_type_log )
      If ( Associated( field_in % Data ) ) Then
        field_out % Data(:,:) = field_in % Data(:,:)
      Else
        field_out % Data_Log(:,:) = field_in % Data_Log(:,:)
      End If

    Case Default
      Cmessage = 'Unsupported Data-Type'
      ErrorStatus = -30
      Call Ereport( RoutineName, ErrorStatus, Cmessage )

  End Select

  Case( interp_all, interp_h_only )

!---------------------------------------------------------------
! If data is compressed onto land points, need to expand it
! To do this neatly requires some fiddling with pointers.
! ptr_field_in and ptr_field_out are the ones to work with. Set
! up here.
!--------------------------------------------------------------
  If ( field_in % stashmaster % grid_type == ppx_atm_compressed) Then
    Allocate( field_in_tmp % Data( grid_in % loc_p_rows *           &
                             grid_in % loc_p_row_length ,           &
                             field_in % levels ) )

    Allocate( field_out_tmp % Data( grid_out % loc_p_rows *         &
                              grid_out % loc_p_row_length,          &
                              field_out % levels ) )

    field_in_tmp % levels       = field_in % levels
    field_in_tmp % rows         = grid_in % loc_p_rows
    field_in_tmp % row_len      = grid_in % loc_p_row_length
    field_in_tmp % level_size   = field_in_tmp % rows *              &
                                  field_in_tmp % row_len
    field_in_tmp % glob_rows    = grid_in % glob_p_rows
    field_in_tmp % glob_row_len = grid_in % glob_p_row_length
    field_in_tmp % glob_level_size = field_in_tmp % glob_rows *      &
                                     field_in_tmp % glob_row_len
    field_in_tmp % stashmaster => field_in % stashmaster

    field_out_tmp % levels       = field_out % levels
    field_out_tmp % rows         = grid_out % loc_p_rows
    field_out_tmp % row_len      = grid_out % loc_p_row_length
    field_out_tmp % level_size   = field_out_tmp % rows *              &
                                   field_out_tmp % row_len
    field_out_tmp % glob_rows    = grid_out % glob_p_rows
    field_out_tmp % glob_row_len = grid_out % glob_p_row_length
    field_out_tmp % glob_level_size = field_out_tmp % glob_rows *      &
                                      field_out_tmp % glob_row_len
    field_out_tmp % stashmaster => field_out % stashmaster

    ! Expand level by level
    Do i = 1, field_in % levels
      Call expand_from_mask( field_in_tmp % Data(1:,i),                   &
                             field_in % Data(1:,i),                    &
                             local_lsm_in,  field_in_tmp % level_size,&
                             field_in % level_size )
    End Do

    ptr_field_in  => field_in_tmp
    ptr_field_out => field_out_tmp
  Else If ( field_in % stashmaster % grid_type == ppx_atm_ozone .AND. &
           ( field_in % glob_row_len  /= 1                      .OR. &
             field_out % glob_row_len /= 1) ) THEN
  ! We still want to make zonal to zonal work.  This is to fix anything to
  ! do with full field ozone.
    Allocate( field_in_tmp % Data( grid_in % loc_p_rows *           &
                             grid_in % loc_p_row_length ,           &
                             field_in % levels ) )

    Allocate( field_out_tmp % Data( grid_out % loc_p_rows *         &
                              grid_out % loc_p_row_length,          &
                              field_out % levels ) )

    field_in_tmp % levels       = field_in % levels
    field_in_tmp % rows         = grid_in % loc_p_rows
    field_in_tmp % row_len      = grid_in % loc_p_row_length
    field_in_tmp % level_size   = field_in_tmp % rows *              &
                                  field_in_tmp % row_len
    field_in_tmp % glob_rows    = grid_in % glob_p_rows
    field_in_tmp % glob_row_len = grid_in % glob_p_row_length
    field_in_tmp % glob_level_size = field_in_tmp % glob_rows *      &
                                     field_in_tmp % glob_row_len
    ! Create new stashmaster entry to change from ozone to full field.
    ALLOCATE(field_in_tmp % stashmaster)
    field_in_tmp % stashmaster   = field_in % stashmaster
    field_in_tmp % stashmaster % grid_type = ppx_atm_tall

    field_out_tmp % levels       = field_out % levels
    field_out_tmp % rows         = grid_out % loc_p_rows
    field_out_tmp % row_len      = grid_out % loc_p_row_length
    field_out_tmp % level_size   = field_out_tmp % rows *              &
                                   field_out_tmp % row_len
    field_out_tmp % glob_rows    = grid_out % glob_p_rows
    field_out_tmp % glob_row_len = grid_out % glob_p_row_length
    field_out_tmp % glob_level_size = field_out_tmp % glob_rows *      &
                                      field_out_tmp % glob_row_len
    ! Create new stashmaster entry to change from ozone to full field.
    ALLOCATE(field_out_tmp % stashmaster)
    field_out_tmp % stashmaster   = field_out % stashmaster
    field_out_tmp % stashmaster % grid_type = ppx_atm_tall

    ! Expand level by level
    If (field_in % glob_row_len == 1) THEN
      ! Zonal ozone
      Do k = 1, field_in_tmp % levels
        Do j = 1, field_in_tmp % rows
          Do i = 1, field_in_tmp % row_len
            field_in_tmp % data((j-1)*field_in_tmp%row_len+i,k) =          &
                                                        field_in % data(j,k)
          End Do
        End Do
      End Do
    Else
      ! Full field ozone
      field_in_tmp % data(:,:) = field_in % data(:,:)
    End If

    ptr_field_in  => field_in_tmp
    ptr_field_out => field_out_tmp

  Else
    ptr_field_in  => field_in
    ptr_field_out => field_out
  End If

!------------------------------------------------------------------
! Can now allocate space for local levels for interpolation
!------------------------------------------------------------------
  Allocate( level_field_in( ptr_field_in % glob_level_size ) )
  Allocate( level_field_out( ptr_field_out % glob_level_size ) )

!-------------------------------------------------------------------
! Need to work out which weights we wish to use. This depends on
! the type of field that is being interpolated.
!-------------------------------------------------------------------

! For soil moisture fields we can have nearest neighbour interpolation.
! This is assuming that bilinear is being used as rcf_select_weights
! enforces this condition.
  orig_h_int_method = h_int_method
  IF ( ptr_field_in % stashmaster % section  == stashcode_prog_sec      .AND.  &
       ( ptr_field_in % stashmaster % item   == stashcode_vol_smc_wilt  .OR.   &
         ptr_field_in % stashmaster % item   == stashcode_vol_smc_cri   .OR.   &
         ptr_field_in % stashmaster % item   == stashcode_vol_smc_sat ) .AND.  &
       smcp_int_nearest_neighbour ) THEN
    h_int_method = nearest_neighbour
  END IF

  Call Rcf_select_weights( ptr_bl_index_b_l, ptr_bl_index_b_r, &
                           ptr_weight_b_l, ptr_weight_b_r,     &
                           ptr_weight_t_l, ptr_weight_t_r,     &
                           ptr_aw_index_targ_lhs,              &
                           ptr_aw_index_targ_top,              &
                           ptr_aw_colat_t, ptr_aw_long_l,      &
                           ptr_field_in % stashmaster % grid_type, &
                           ptr_field_in % stashmaster % section,   &
                           ptr_field_in % stashmaster % item )

!--------------------------------------------------------------------
! We need to gather <nproc> levels onto pes.
!--------------------------------------------------------------------

  div = ptr_field_in % levels / nproc
  rem = Mod( ptr_field_in % levels, nproc )
  pe = 0


  Do i = 1, div

    Call Change_Decomposition( decomp_rcf_input )

    Do j = ((i-1) * nproc) + 1, i * nproc
      ! Will gather level j on pe

      If (ptr_field_in % glob_row_len == 1) Then ! Zonal Data
        Call Rcf_Gather_Zonal_Field( ptr_field_in % Data(:,j),         &
                                     level_field_in,                   &
                                     ptr_field_in % level_size,        &
                                     ptr_field_in % glob_level_size, 1,&
                                     ppx_atm_tzonal, pe )

      Else

        Call Rcf_Gather_Field( ptr_field_in % Data(:,j),              &
                               level_field_in,                        &
                               ptr_field_in % row_len,                &
                               ptr_field_in % rows,                   &
                               ptr_field_in % glob_row_len,           &
                               ptr_field_in % glob_rows, pe,          &
                               gc_all_proc_group )
      End If

      pe = pe + 1
      if (pe == nproc) pe = 0
    End Do

!-------------------------------------------------------------------
! All PEs are currently full of level data - so do the interpolation
!-------------------------------------------------------------------
    Call Rcf_H_Int_Ctl( ptr_field_out % glob_level_size,                 &
                        ptr_field_in % glob_row_len,                     &
                        ptr_field_out % glob_row_len,                    &
                        ptr_field_in % glob_rows,                        &
                        ptr_field_out % glob_rows,                       &
                        grid_out % global, ptr_aw_index_targ_lhs,        &
                        ptr_aw_index_targ_top, ptr_bl_index_b_l,         &
                        ptr_bl_index_b_r, ptr_aw_colat_t, ptr_aw_long_l, &
                        level_field_in, ptr_weight_t_r, ptr_weight_b_r,  &
                        ptr_weight_t_l, ptr_weight_b_l, level_field_out )
  
!-------------------------------------------------------------------
! Coastal Adjustment for land-only, sea-only and T*  fields
!-------------------------------------------------------------------
! DEPENDS ON: timer
    If (LTimer) Call Timer( 'Coastal adjustment', 103)
  
    If (ptr_field_out % stashmaster % grid_type == ppx_atm_compressed.OR.   &
        ptr_field_out % stashmaster % grid_type == ppx_atm_tsea .OR.        &
        ( ptr_field_out % stashmaster % section == stashcode_prog_sec .AND. &
          ptr_field_out % stashmaster % item    == stashcode_tstar )) Then
      Do j = 1, n_coastal_points
        level_field_out( coast_index_out( j ) ) = &
                         level_field_in( coast_index_in( j ) )
      End Do
  
      If (LSpiral_S) Then     ! Spiral adjustment
        maxdim = min(grid_out % glob_p_rows, grid_out % glob_p_row_length)
        Do j = 1, grid_out % glob_p_rows * grid_out % glob_p_row_length
          index_targ_sea_tmp( j )  = index_targ_sea( j )
          index_targ_land_tmp( j ) = index_targ_land( j )
          If ( glob_lsm_out( j ) ) Then
            lsm_tmp( j ) = 1
          Else
            lsm_tmp( j ) = 0
          End If
        End Do
  
        sea_points_unres_tmp  = n_sea_points_unres
        land_points_unres_tmp = n_land_points_unres
  
        If (ptr_field_out % stashmaster % grid_type /=                   &
                                          ppx_atm_compressed) Then
! Only do the sea coastal points if the field isn't land only.
! DEPENDS ON: intf_coast_aj
          Call Intf_Coast_AJ( lsm_tmp, index_targ_sea_tmp,               &
               sea_points_unres_tmp, grid_out % glob_p_rows,             &
               grid_out % glob_p_row_length, level_field_out, 0, cyclic, &
               maxdim )
        End If

! DEPENDS ON: intf_coast_aj
        Call Intf_Coast_AJ( lsm_tmp, index_targ_land_tmp,              &
             land_points_unres_tmp, grid_out % glob_p_rows,            &
             grid_out % glob_p_row_length, level_field_out, 1, cyclic, &
             maxdim )

      Else    ! Non-spiral adjustment

        Do j = 1, n_land_points_unres
          level_field_out( index_targ_land( j ) ) = &
                           level_field_out( land_unres_index( j ) )
        End Do

        Do j = 1, n_sea_points_unres
          level_field_out( index_targ_sea( j ) ) = &
                         level_field_out( sea_unres_index( j ) )
        End Do

      End If ! Spiral
    End If

! DEPENDS ON: timer
    If (LTimer) Call Timer( 'Coastal adjustment', 104)

!-----------------------------------------------------------------
! Average the polar rows
!-----------------------------------------------------------------
    If ( ptr_field_out % stashmaster % grid_type <= 3 .OR.               &
         ptr_field_out % stashmaster % grid_type == 21 ) Then
      If ( ptr_field_out % glob_rows > grid_out % glob_v_rows ) Then
        Call Rcf_Average_polar(level_field_out, ptr_field_out % glob_rows, &
                             ptr_field_out % glob_row_len,                 &
                             grid_out % global, averaged )
      End If

      If ( averaged .AND. PrintStatus >= PrStatus_Normal) Then
        field_averaged = 1
      End If
    End If
    If ( field_out % stashmaster % grid_type == ppx_atm_ozone ) Then
      If (field_out % glob_row_len == 1 .AND. field_in % glob_row_len /= 1) Then
        ! Need to use field_out since we have change stashmaster and row_len to
        ! full field.
        Do j = 1, field_out % rows
          ! We can reuse polar rows and average over 1 row
          Call Rcf_Average_polar(                                        &
                      level_field_out((j-1)*ptr_field_out % glob_row_len:&
                                      j*ptr_field_out % glob_row_len),   &
                                  1,                                     &
                                  ptr_field_out % glob_row_len,          &
                                  grid_out % global, averaged )

        End Do
      End If

    End If

!-------------------------------------------------------------------
! And re-scatter the data back to original PEs
!------------------------------------------------------------------
    Call Change_Decomposition( decomp_rcf_output )
    Do j = ((i-1) * nproc) + 1, i * nproc

      If (ptr_field_out % glob_row_len == 1) Then ! Zonal Data
        Call Rcf_Scatter_Zonal_Field( ptr_field_out % Data(:,j),        &
                                  level_field_out,                      &
                                  ptr_field_out % level_size,           &
                                  ptr_field_out % glob_level_size, 1,   &
                                  ppx_atm_tzonal, pe )
      Else
        Call Rcf_Scatter_Field( ptr_field_out % Data(:,j),              &
                                level_field_out,                        &
                                ptr_field_out % row_len,                &
                                ptr_field_out % rows,                   &
                                ptr_field_out % glob_row_len,           &
                                ptr_field_out % glob_rows, pe,          &
                                gc_all_proc_group )
      End If
      pe = pe + 1
      If (pe == nproc) pe = 0
    End Do
  End Do

!-------------------------------------------------------------------
! There are rem levels left to process. Will do these now.
!-------------------------------------------------------------------
  Call Change_Decomposition( decomp_rcf_input )
  pe = 0
  Do i = 1, rem
    j = nproc * div + i

    If (ptr_field_in % glob_row_len == 1) Then ! Zonal Data
      Call Rcf_Gather_Zonal_Field( ptr_field_in % Data(:,j),            &
                                   level_field_in,                      &
                                   ptr_field_in % level_size,           &
                                   ptr_field_in % glob_level_size, 1,   &
                                   ppx_atm_tzonal, pe )
    Else
      Call Rcf_Gather_Field( ptr_field_in % Data(:,j),                  &
                             level_field_in,                            &
                             ptr_field_in % row_len,                    &
                             ptr_field_in % rows,                       &
                             ptr_field_in % glob_row_len,               &
                             ptr_field_in % glob_rows, pe,              &
                             gc_all_proc_group )
    End If

    pe = pe + 1
  End Do

  If (mype < pe) Then
    Call Rcf_H_Int_Ctl( ptr_field_out % glob_level_size,                 &
                        ptr_field_in % glob_row_len,                     &
                        ptr_field_out % glob_row_len,                    &
                        ptr_field_in % glob_rows,                        &
                        ptr_field_out % glob_rows,                       &
                        grid_out % global, ptr_aw_index_targ_lhs,        &
                        ptr_aw_index_targ_top, ptr_bl_index_b_l,         &
                        ptr_bl_index_b_r, ptr_aw_colat_t, ptr_aw_long_l, &
                        level_field_in, ptr_weight_t_r, ptr_weight_b_r,  &
                        ptr_weight_t_l, ptr_weight_b_l, level_field_out )

!-------------------------------------------------------------------
! Coastal Adjustment for land-only, sea-only and T*  fields
!-------------------------------------------------------------------
! DEPENDS ON: timer
    If (LTimer) Call Timer( 'Coastal adjustment', 103)
    If (ptr_field_out % stashmaster % grid_type == ppx_atm_compressed.OR.   &
        ptr_field_out % stashmaster % grid_type == ppx_atm_tsea .OR.        &
        ( ptr_field_out % stashmaster % section == stashcode_prog_sec .AND. &
          ptr_field_out % stashmaster % item    == stashcode_tstar )) Then
      Do j = 1, n_coastal_points
        level_field_out( coast_index_out( j ) ) = &
                         level_field_in( coast_index_in( j ) )
      End Do

      If (LSpiral_S) Then     ! Spiral adjustment
        maxdim = min( grid_out % glob_p_rows,grid_out % glob_p_row_length)
        Do j = 1, grid_out % glob_p_rows * grid_out % glob_p_row_length
          index_targ_sea_tmp( j )  = index_targ_sea( j )
          index_targ_land_tmp( j ) = index_targ_land( j )
          If ( glob_lsm_out( j ) ) Then
            lsm_tmp( j ) = 1
          Else
            lsm_tmp( j ) = 0
          End If
        End Do

        sea_points_unres_tmp  = n_sea_points_unres
        land_points_unres_tmp = n_land_points_unres

        If (ptr_field_out % stashmaster % grid_type /=                   &
                                          ppx_atm_compressed) Then
! Only do the sea coastal points if the field isn't land only.
! DEPENDS ON: intf_coast_aj
          Call Intf_Coast_AJ( lsm_tmp, index_targ_sea_tmp,               &
               sea_points_unres_tmp, grid_out % glob_p_rows,             &
               grid_out % glob_p_row_length, level_field_out, 0, cyclic, &
               maxdim )
        End If

! DEPENDS ON: intf_coast_aj
        Call Intf_Coast_AJ( lsm_tmp, index_targ_land_tmp,              &
             land_points_unres_tmp, grid_out % glob_p_rows,            &
             grid_out % glob_p_row_length, level_field_out, 1, cyclic, &
             maxdim )

      Else    ! Non-spiral adjustment

        Do j = 1, n_land_points_unres
          level_field_out( index_targ_land( j ) ) = &
                           level_field_out( land_unres_index( j ) )
        End Do

        Do j = 1, n_sea_points_unres
          level_field_out( index_targ_sea( j ) ) = &
                           level_field_out( sea_unres_index( j ) )
        End Do

      End If ! Spiral
    End If

! DEPENDS ON: timer
    If (LTimer) Call Timer( 'Coastal adjustment', 104)

!-----------------------------------------------------------------
! Average the polar rows
!-----------------------------------------------------------------
    If ( ptr_field_out % stashmaster % grid_type <= 3 .OR.               &
         ptr_field_out % stashmaster % grid_type == 21 ) Then
! Might need to do something here for ENDGAME grid... theta points not at poles
! anymore.
      If ( ptr_field_out % glob_rows > grid_out % glob_v_rows ) Then
        Call Rcf_Average_polar( level_field_out, ptr_field_out % glob_rows, &
                                ptr_field_out % glob_row_len,               &
                                grid_out % global, averaged )
      End If

      If ( averaged .AND. PrintStatus >= PrStatus_Normal) Then
        field_averaged = 1
      End If
    End If
    If ( field_out % stashmaster % grid_type == ppx_atm_ozone ) Then
      If (field_out % glob_row_len == 1 .AND. field_in % glob_row_len /= 1) Then
        ! Need to use field_out since we have change stashmaster and row_len to
        ! full field.
        Do j = 1, field_out % rows
          ! We can reuse polar rows and average over 1 row
          Call Rcf_Average_polar(                                        &
                      level_field_out((j-1)*ptr_field_out % glob_row_len:&
                                      j*ptr_field_out % glob_row_len),   &
                                  1,                                     &
                                  ptr_field_out % glob_row_len,          &
                                  grid_out % global, averaged )

        End Do
      End If
    End If

  End If

!-------------------------------------------------------------------
! And re-scatter data
!-------------------------------------------------------------------

  Call Change_Decomposition( decomp_rcf_output )
  pe = 0
  Do i = 1, rem
    j = nproc * div + i
    If (ptr_field_out % glob_row_len == 1) Then ! Zonal Data
      Call Rcf_Scatter_Zonal_Field( ptr_field_out % Data(:,j),        &
                                level_field_out,                      &
                                ptr_field_out % level_size,           &
                                ptr_field_out % glob_level_size, 1,   &
                                ppx_atm_tzonal , pe )
    Else
      Call Rcf_Scatter_Field( ptr_field_out % Data(:,j),              &
                              level_field_out,                        &
                              ptr_field_out % row_len,                &
                              ptr_field_out % rows,                   &
                              ptr_field_out % glob_row_len,           &
                              ptr_field_out % glob_rows, pe,          &
                              gc_all_proc_group )
    End If
    pe = pe + 1
  End Do

!----------------------------------------------------------------
! Print out a message if required for polar averaged fields
!----------------------------------------------------------------

  Call GC_IMAX( 1, nproc, stat, field_averaged )

  If ( field_averaged == 1 .AND. mype == 0 .AND.                       &
       PrintStatus >= PrStatus_Normal ) Then
    Write (6,*) 'Interpolation has required averaging of polar rows ', &
                'for section: ', ptr_field_out % stashmaster % section,&
                ', item: ', ptr_field_out % stashmaster % item
  End If

!----------------------------------------------------------------
! Deallocate levels
!----------------------------------------------------------------
  Deallocate( level_field_in )
  Deallocate( level_field_out )

!-----------------------------------------------------------------
! Reverse process above, take fields back to land points
!-----------------------------------------------------------------

  If (field_out % stashmaster % grid_type == ppx_atm_compressed) Then
    Do i = 1, field_out % levels
      Call compress_to_mask( field_out_tmp % Data(1:,i),                 &
                           field_out % Data(1:,i),                     &
                           local_lsm_out, field_out_tmp % level_size, &
                           size )

      If (size /= field_out % level_size ) Then
        Cmessage = 'Recompression onto land points - sizes mismatch'
        ErrorStatus = 60
        Call Ereport( RoutineName, ErrorStatus, Cmessage )
      End If

    End Do

  ! Release temp. memory
    Deallocate( field_out_tmp % Data )
    Deallocate( field_in_tmp % Data )
  Else If ( field_in % stashmaster % grid_type == ppx_atm_ozone .AND. &
           (field_in % glob_row_len  /= 1                        .OR. &
            field_out % glob_row_len /= 1 )) THEN
    If (field_out % glob_row_len == 1) Then
      ! Zonal ozone
      Do k = 1, field_out % levels
        Do j = 1, field_out % rows
          ! This has been averaged earlier.
          field_out % data(j,k) =                                     &
                field_out_tmp % data(1+(j-1)*field_out_tmp%row_len,k)
        End Do
      End Do
    Else
      ! Full field
      field_out % data(:,:) = field_out_tmp % data(:,:)
    End If
    ! Release temporary stashmasters
    Deallocate( field_out_tmp % stashmaster)
    Deallocate( field_in_tmp % stashmaster)
    ! Release temp. memory
    Deallocate( field_out_tmp % Data )
    Deallocate( field_in_tmp % Data )
  End If

  Nullify( ptr_field_out )
  Nullify( ptr_field_in  )

 ! Earlier check if using soil could have changed this.  Set it back for
 ! future fields. 
  h_int_method = orig_h_int_method

  Case( interp_no_op)
    ! do nothing

End Select

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4)

Return
End Subroutine Rcf_horizontal

End Module Rcf_horizontal_Mod
