! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Rcf_Calc_P_Star  - calculates P* from rho and exner

Module Rcf_Calc_P_Star_Mod

!  Subroutine Rcf_Calc_P_Star - calculates P* from rho and exner
!
! Description:
!    This subroutine calculates P* (surface pressure)
!
! Method:
!     Code derived New Dynamics code. Uses rho and exner and
!     calculates heights internally based on grid parameters and
!     orography.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Calc_P_Star( grid, fields, field_count, hdr, decomp,  &
                            orog, p_star )

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_rho,             &
    stashcode_exner,           &
    stashcode_p,               &
    stashcode_prog_sec

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_Grid_Type_Mod, Only : &
    grid_type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Exner_P_Convs_Mod, Only : &
    Rcf_Conv_Exner_P

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Generate_Heights_Mod, Only : &
    height_gen_original,             &
    height_gen_smooth

Use Rcf_field_equals_mod, Only : &
    Rcf_Field_Equals

USE earth_constants_mod, ONLY: g, earth_radius

IMPLICIT NONE

! Arguments
Integer, Intent(In)                :: field_count   ! number of fields
Integer, Intent(In)                :: decomp        ! decomposition
Type( grid_type ), Intent(In)      :: grid
Type( field_type ), Pointer        :: fields(:)
Type( field_type ), Intent(In)     :: orog          ! orography field
Type( UM_header_type ), Intent(In) :: hdr           ! dump header
Type( field_type ), Intent(InOut)  :: p_star        ! surface pressure

! Local variables
Integer                             :: i             ! looper
Integer                             :: pos           ! field position
Real                                :: first_rho_level
Type( field_type )                  :: exner
Type( field_type ), Pointer         :: rho

Integer                             :: ErrorStatus
Character (Len=*), Parameter        :: RoutineName = 'Rcf_Calc_Pstar'
Character (Len=80)                  :: Cmessage

!-------------------------------------------------------------------
! Find rho in the fields array and read it in
!-------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_rho,                  &
                 fields, field_count, pos )
rho => fields(pos)
Call Rcf_Alloc_Field( rho )
Call Rcf_Read_Field( rho, hdr, decomp )

!------------------------------------------------------------------
! Find, read and convert exner to P - OR find and read P.
! Will be P in exner variable eventually however.
!------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_exner,                &
                 fields, field_count, pos, .TRUE. )

If (pos /= 0) Then
  Call Rcf_Field_Equals( exner, fields(pos) )
  Call Rcf_Alloc_Field( exner )
  Call Rcf_Read_Field( exner, hdr, decomp )
  Call Rcf_Conv_Exner_P( exner )
Else        ! have P in dump
  Call Rcf_Locate( stashcode_prog_sec, stashcode_p,                  &
                   fields, field_count, pos )
  Call Rcf_Field_Equals( exner, fields(pos) )
  Call Rcf_Alloc_Field( exner )
  Call Rcf_Read_Field( exner, hdr, decomp )
End If

!-------------------------------------------------------------------
! Calculate the first rho level height and thus P*
!-------------------------------------------------------------------
Do i = 1, orog % level_size
  If ( grid % height_gen_method == height_gen_original ) Then

  first_rho_level = orog % Data(i,1) + Earth_Radius +                &
                    grid % eta_rho_levels(1) * grid % z_top_of_model

  Else If ( grid % height_gen_method == height_gen_smooth ) Then

    first_rho_level = Earth_Radius + grid % eta_rho_levels(1) *        &
                      grid % z_top_of_model + orog % Data(i,1) *       &
                      (1.0 - grid % eta_rho_levels(1) /                &
            grid % eta_rho_levels(grid % first_constant_r_rho_level))**2

  Else
    ErrorStatus = 10
    Cmessage = 'Height generation method unknown'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If


  p_star % Data(i,1) = exner % Data(i,1) + G * rho % Data(i,1) *    &
                      (first_rho_level -                            &
                      (orog % Data(i,1) + Earth_Radius)) /          &
                      ( first_rho_level * first_rho_level )
End Do

!-------------------------------------------------------------------
! Remove the no longer needed exner and rho fields
!-------------------------------------------------------------------
Call Rcf_Dealloc_Field( exner )
Call Rcf_Dealloc_Field( rho )
Nullify( rho )

Return
End Subroutine Rcf_Calc_P_Star
End Module Rcf_Calc_P_Star_Mod

