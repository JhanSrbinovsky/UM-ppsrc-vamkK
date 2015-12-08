! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+  Rcf_Calc_Output_Exner_Mod - calculate output dump exner pressure

Module Rcf_Calc_Output_Exner_Mod

!  Subroutine Rcf_Calc_Output_Exner
!
! Description:
!    This module calculates exner pressure for the output dump
!
! Method:
!    This code is derived from that used in the New Dynamics
!    reconfiguration program. First P* is calculated on the input
!    grid and then interpolated horizontally onto the output grid.
!    This is then used to calculate the 1st level of exner which
!    in turn is used to calculate the higher levels.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Calc_Output_Exner( fields_in, field_count_in, orog_in, &
                                  hdr_in, orog_out, T, Q, Exner,      &
                                  orog_source, r_rho_levels,          &
                                  r_theta_levels )
USE UM_ParVars, Only : &
    mype

USE PrintStatus_mod, Only : &
    PrintStatus,                &
    PrStatus_Normal

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_UMhead_Mod, Only : &
    UM_header_type

Use Rcf_Grid_Type_Mod, Only: &
    Input_Grid,          &
    Output_Grid

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_pstar,           &
    stashcode_prog_sec

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_h_only,                   &
    interp_copy,                     &
    interp_no_op

Use Submodel_Mod, Only :             &
    Atmos_IM

Use Rcf_Calc_P_Star_Mod, Only : &
    Rcf_Calc_P_Star

USE decomp_params, ONLY : &
    Decomp_rcf_input

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Interpolate_Mod, Only : &
    Rcf_Interpolate

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Adjust_Pstar_Mod, Only : &
    Rcf_Adjust_Pstar

Use Rcf_Data_Source_Mod, Only : &
    Ancillary_File

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

USE rcf_field_equals_mod, ONLY : &
    rcf_field_equals

USE earth_constants_mod, ONLY: g    

USE atmos_constants_mod, ONLY: cp, kappa, pref, repsilon

Implicit None

! Arguments
Integer, Intent(In)                :: field_count_in
Type( field_type ), Pointer        :: fields_in(:)
Type( field_type ), Intent(InOut)  :: orog_in
Type( field_type ), Intent(In)     :: orog_out
Type( field_type ), Intent(In)     :: T
Type( field_type ), Intent(In)     :: Q
Type( field_type ), Intent(InOut)  :: Exner
Type( UM_header_type ), Intent(In) :: hdr_in
Real, Intent(In)                   :: r_rho_levels(               &
                                      output_grid % loc_p_field,  &
                                      0 : output_grid % model_levels+1)
Real, Intent(In)                   :: r_theta_levels(             &
                                      output_grid % loc_p_field,  &
                                      0 : output_grid % model_levels+1)
Integer, Intent(In)                :: orog_source

! Local variables
Integer                            :: i          ! Looper
Integer                            :: k          ! Looper
Integer                            :: pos        ! position in fields_in
Integer                            :: k_off      ! Offset due to 0th level in T
Real                               :: temp       ! temporary value
Real                               :: weight1    ! averaging weight
Real                               :: weight2    ! averaging weight
Real                               :: deltaz     ! vertical spacing
Type( field_type )                 :: pstar_in
Type( field_type )                 :: pstar_out
Type( field_type )                 :: dummy

! Comdecks
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

If (mype == 0 .AND. PrintStatus >= PrStatus_Normal) Then
  Write (6,*) 'Calculating Exner'
End If

!--------------------------------------------------------------------
! Setup p_star_in field - either read or calculate from input dump
!--------------------------------------------------------------------
! Check input dump for p_star
Call Rcf_Locate( stashcode_prog_sec, stashcode_Pstar,                &
                 fields_in, field_count_in, pos, .True. )

! If present, read in, else calculate.
If (pos /= 0) Then
  Pstar_in  = fields_in(pos)
  Call Rcf_Alloc_Field( Pstar_in )
  Call Rcf_Read_Field( Pstar_in, Hdr_In, decomp_rcf_input )

Else

  CALL rcf_field_equals(pstar_in, input_grid)
  pstar_in % stashmaster     => Rcf_Exppx(Atmos_IM, 0, stashcode_pstar, &
                                          stash_in_arg = .TRUE. )

  Call Rcf_Alloc_Field( pstar_in )
  Call Rcf_Calc_P_Star( Input_Grid, fields_in, field_count_in, hdr_in, &
                        decomp_rcf_input, orog_in, pstar_in )

End If

! Only need to do horizontal interpolation of P* if h_int_active
If (h_int_active) Then
  pstar_in % interp          = interp_h_only
Else
  pstar_in % interp          = interp_copy
End If

!--------------------------------------------------------------------
! Setup p_star_out field - interpolated from p_star_in
!--------------------------------------------------------------------
CALL rcf_field_equals(pstar_out, output_grid)
pstar_out % stashmaster     => Rcf_Exppx( Atmos_IM, 0, stashcode_pstar )

Call Rcf_Alloc_Field( pstar_out )
Call Rcf_Interpolate( pstar_in, pstar_out, Input_Grid, Output_Grid, &
                      dummy, dummy)

!--------------------------------------------------------------------
! If required, adjust the p* field for orographic effects
!--------------------------------------------------------------------
If ( orog_source == Ancillary_File ) Then
  Call Rcf_Adjust_Pstar( pstar_out, t, orog_in, orog_out, Input_Grid, &
                         Output_Grid, r_theta_levels )
End If

!--------------------------------------------------------------------
! calculate exner
! equation is
!
! exner(k) (Cp T_v(k-1) + g weight1 (r(k)-r(k-1)) =
! exner(k-1) (Cp T_v(k-1) - (1-weight1) (r(k)-r(k-1))
!
! where the weights represent the interpolation coefficients
! required to calculate exner at a T point.
! Note: In the new integration scheme Theta below level 1 is assumed
!       to be equal to theta at level 1. The equation is thus
!
! exner(1) (Cp T_v(1) + g (r(k)-r(k-1)) =
! exner* (Cp T_v(1) )
!--------------------------------------------------------------------

IF (T % bottom_level == 0 .AND. q % bottom_level == 0) THEN
  k_off = 1
ELSE
  k_off = 0
END IF

!--------------------------------------------------------------------
! Rho level 1
!--------------------------------------------------------------------

k = 1
Do i = 1, pstar_out % level_size

  temp = T % Data(i,k+k_off) * (1. + (1./repsilon -1.) * q % Data(i,k+k_off))
  exner % Data(i,k) = (CP * temp *                                &
                      (pstar_out % Data(i,1) / pref)**kappa ) /   &
                      (CP * temp + g * (r_rho_levels(i,k) -       &
                        r_theta_levels(i,k-1) ) )

End Do

!--------------------------------------------------------------------
! Rho levels 2 to model_levels + 1
!--------------------------------------------------------------------
Do k = 2, exner % levels
  If ( k <= Output_Grid % wet_levels + 1 ) Then  ! Wet Level
    Do i = 1, exner % level_size
      If (k > Output_Grid % model_levels) Then

        ! extra pressure levels is same height above top theta as
        ! pressure below is below top theta - hence weights are 0.5
        weight1 = 0.5
        deltaz  = 2.0 * (r_theta_levels(i,k-1) - r_rho_levels(i,k-1))

      Else

        deltaz  = r_rho_levels(i,k) - r_rho_levels(i,k-1)
        weight1 = (r_rho_levels(i,k) - r_theta_levels(i,k-1) ) / deltaz

      End If ! (K > model_levels)

      weight2 = 1. - weight1

      temp = T % Data(i,k-1+k_off) * (1. + (1./repsilon -1.)                 &
                             * q % Data(i,k-1+k_off))

      exner % Data(i,k) = (CP * temp - g * weight1 * deltaz) *         &
                          exner % Data(i,k-1) /                        &
                          (CP * temp + g * weight2 * deltaz)

    End Do
  Else    ! Dry Levels
    Do i = 1, exner % level_size
      If (k > Output_Grid % model_levels ) Then

        ! extra pressure levels is same height above top theta as
        ! pressure below is below top theta - hence weights are 0.5
        weight1 = 0.5
        deltaz = 2.0 * (r_theta_levels(i,k-1) - r_rho_levels(i,k-1))

      Else
        deltaz  = r_rho_levels(i,k) - r_rho_levels(i,k-1)
        weight1 = (r_rho_levels(i,k) - r_theta_levels(i,k-1) ) / deltaz

      End If ! (k > model_levels)

      weight2 = 1. - weight1

      temp = T % Data(i,k-1+k_off)
      exner % Data(i,k) = (CP * temp - g * weight1 * deltaz) *       &
                          exner % Data(i,k-1) /                      &
                          (CP * temp + g * weight2 * deltaz)

    End Do
  End If ! wet/dry layer
End Do ! k


!--------------------------------------------------------------------
! Clear up memory for the pstar fields
!--------------------------------------------------------------------
Call Rcf_Dealloc_Field( pstar_in )
Call Rcf_Dealloc_Field( pstar_out )

Return
End Subroutine Rcf_Calc_Output_Exner
End Module Rcf_Calc_Output_Exner_Mod
