! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Choose weights for horizontal interpolation.

Module Rcf_select_weights_mod

!  Subroutine Rcf_Select_Weights - choose weights for horiz. interp.
!
! Description:
! This module calculates which part of the weights array is
! required by a given horizontal interpolation method, based
! on the type of field being used.
!
! Method:
! The correct weights are chosen on grid-code and interpolation
! method chosen. Pointers are set to relevant weight array slices.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_select_weights( ptr_bl_index_b_l, ptr_bl_index_b_r, &
                               ptr_weight_b_l, ptr_weight_b_r,     &
                               ptr_weight_t_l, ptr_weight_t_r,     &
                               ptr_aw_index_targ_lhs,              &
                               ptr_aw_index_targ_top,              &
                               ptr_aw_colat_t, ptr_aw_long_l,      &
                               gridcode, section, item )

Use Rcf_Interp_Weights_Mod, Only : &
    bilinear,                  &
    area_weighted,             &
    nearest_neighbour,         &
    h_int_method,              &
    weight_b_l, weight_b_r,    &
    weight_t_l, weight_t_r,    &
    bl_index_b_l, bl_index_b_r,&
    aw_index_targ_lhs,         &
    aw_index_targ_top,         &
    aw_colat_t, aw_long_l

Use Rcf_Stashcodes_Mod, Only :                        &
    stashcode_prog_sec,      stashcode_tracer_sec,    &
    stashcode_ukca_sec,                               &
    stashcode_u,             stashcode_surf_z_curr,   &
    stashcode_w,             stashcode_3d_cca,        &
    stashcode_v,             stashcode_surf_m_curr,   &
    stashcode_riv_sequence,  stashcode_3d_ccw,        &
    stashcode_u_compnt_pert, stashcode_v_compnt_pert, &
    stashcode_3d_nat_so2_em, stashcode_3d_oh_conc

Use Ereport_mod, Only : &
    Ereport

Use cppxref_mod, Only :                               &
    ppx_atm_tall,            ppx_atm_uall,            &
    ppx_atm_cuall,           ppx_atm_cvall,           &
    ppx_atm_tsea,            ppx_atm_usea,            &
    ppx_atm_ozone,           ppx_atm_compressed

IMPLICIT NONE

! Arguments
Integer, Intent(In)           :: gridcode
Integer, Intent(In)           :: section
Integer, Intent(In)           :: item
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

! Local Data
Integer                      :: index
Integer                      :: ErrorStatus
Character (Len=*), Parameter :: RoutineName='Select_Weights'
Character (Len=80)           :: Cmessage

!------------------------------------------------------------------
! Choose the index based on the grid code
!------------------------------------------------------------------
Select Case ( gridcode )
  Case( ppx_atm_tall, ppx_atm_tsea, ppx_atm_compressed )
    index = 1

  Case( ppx_atm_uall, ppx_atm_cuall, ppx_atm_usea )

    Select Case ( section )

    Case ( stashcode_prog_sec )
      Select Case ( item )
! Removed items 135, 139, 148 at vn6.6 since no STASHmaster record exists
      Case (stashcode_u,             & 
            stashcode_surf_z_curr,   &
            stashcode_3d_nat_so2_em, &
            stashcode_w,             &
            stashcode_u_compnt_pert, &
            stashcode_3d_cca )
        index = 2

! Removed items 136, 140 and 149 at vn6.6 since no STASHmaster record exists
      Case (stashcode_v,             & 
            stashcode_surf_m_curr,   &
            stashcode_3d_oh_conc,    &
            stashcode_riv_sequence,  &
            stashcode_v_compnt_pert, &
            stashcode_3d_ccw )
        index = 3

      Case Default
        index = 2
      End Select
    Case Default
      ! Set other sections to default.
      index = 2
    End Select

  Case( ppx_atm_cvall )
    index = 3

  Case( ppx_atm_ozone )
    ! This assumes zonal ozone - not sure how to distinguish otherwise
    index = 4

  Case Default

End Select

!------------------------------------------------------------------
! Set pointers for relevant weight slices
!------------------------------------------------------------------


Select Case ( h_int_method )
  Case ( bilinear, nearest_neighbour )
    ptr_bl_index_b_l => bl_index_b_l( :, index )
    ptr_bl_index_b_r => bl_index_b_r( :, index )
    ptr_weight_b_l   => weight_b_l  ( :, index )
    ptr_weight_b_r   => weight_b_r  ( :, index )
    ptr_weight_t_l   => weight_t_l  ( :, index )
    ptr_weight_t_r   => weight_t_r  ( :, index )
    Nullify( ptr_aw_index_targ_lhs )
    Nullify( ptr_aw_index_targ_top )
    Nullify( ptr_aw_colat_t )
    Nullify( ptr_aw_long_l )


  Case ( area_weighted )
    ptr_aw_index_targ_lhs => aw_index_targ_lhs( :, index )
    ptr_aw_index_targ_top => aw_index_targ_top( :, index )
    ptr_aw_colat_t        => aw_colat_t       ( :, index )
    ptr_aw_long_l         => aw_long_l        ( :, index )
    Nullify( ptr_bl_index_b_l )
    Nullify( ptr_bl_index_b_r )
    Nullify( ptr_weight_b_l)
    Nullify( ptr_weight_b_r)
    Nullify( ptr_weight_t_l)
    Nullify( ptr_weight_t_r)
  Case Default
    Cmessage = 'Interpolation method not yet supported'
    ErrorStatus = 20
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select

Return
End Subroutine Rcf_select_weights

End Module Rcf_select_weights_mod
