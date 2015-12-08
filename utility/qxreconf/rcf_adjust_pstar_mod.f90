! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Performs orographic adjustment of P*

Module Rcf_Adjust_Pstar_Mod

!  Subroutine Rcf_Adjust_Pstar    Adjusts P* when an ancillary
!                                 orography has been used with an
!                                 interpolated P*
!
! Description:
!   A similar method to that used in 4.5 is employed, with a
!   change to a height based vertical co-ordinate frame.
!
! Method:
!   Assuming a constant lapse rate, a surface temperature (Ts) is
!   calculated based on a temperature at a reference height (Tr)
!              Ts = Tr + lapse * (Zr - Z0_at_old_orography)
!
!   P* is then adjusted thus:-
!     P*' = P* [ Ts - lapse ( Zr - Z0_new_orography) ]  ** g/(R lapse)
!              [ ----------------------------------- ]
!              [                 Ts                  ]
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Adjust_Pstar( pstar, t, orog_in, orog_out, Input_Grid,  &
                             Output_Grid, r_theta_levels)

Use Rcf_Grid_Type_Mod, Only : &
    Grid_Type

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Set_Interp_Flags_Mod, Only : &
    interp_h_only,                   &
    interp_copy,                     &
    interp_no_op

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Interpolate_Mod, Only : &
    Rcf_Interpolate

Use Rcf_Field_Equals_Mod, Only : &
    Rcf_Field_Equals

USE earth_constants_mod, ONLY: g, earth_radius

USE atmos_constants_mod, ONLY: r, lapse

USE conversions_mod, ONLY : &
    upperheight

Implicit None

! Arguments
Type( Grid_Type ), Intent(In)       :: Input_Grid
Type( Grid_Type ), Intent(In)       :: Output_Grid
Type( Field_Type ), Intent(InOut)   :: pstar
Type( Field_Type ), Intent(InOut)   :: orog_in
Type( Field_Type ), Intent(In)      :: orog_out
Type( Field_Type ), Intent(In)      :: t
Real, Intent(In)                    :: r_theta_levels(                 &
                                           output_grid % loc_p_field,  &
                                       0 : output_grid % model_levels )

! Local Variables/Parameters
Type( Field_Type )    :: orog_interp         ! interpolated input
                                             ! orography
Type( Field_Type )    :: dummy               ! dummy field
Real                  :: g_over_lapse_r
Real                  :: Tr                  ! T at ref. level.
Real                  :: Zr                  ! Height at ref. level.
Real                  :: Z0i                 ! Interpolated orog.
Real                  :: Z0o                 ! Output orog.
Integer               :: ref_lev             ! first theta level above
                                             ! upperheight
Integer               :: i                   ! Looper

! Comdecks

!--------------------------------------------------------------------
! Set up the orography fields to obtain interpolated orography
!--------------------------------------------------------------------
Call Rcf_Field_Equals( orog_interp, orog_out )

If (h_int_active) Then
  orog_in % interp   = interp_h_only
Else
  orog_in % interp   = interp_copy
End If

orog_interp % interp = interp_no_op

Call Rcf_Alloc_Field( orog_interp )
Call Rcf_Interpolate( orog_in, orog_interp, Input_Grid, Output_Grid, &
                      dummy, dummy )

! Similar calculation in calc_pmsl (could be merged into a common function?)
DO i = 1, output_grid % model_levels
! Upperheight is approximated to be where temperature is free from surface
! effects.  NB Model levels are terrain following near ground so height would be
! lower over land.
  ref_lev = i
  IF ( output_grid % z_top_of_model * output_grid % eta_theta_levels(i) &
       > upperheight ) EXIT
END DO

!--------------------------------------------------------------------
! Can now do the p* adjustment calculation as specified above
!--------------------------------------------------------------------
g_over_lapse_r = g/(lapse * R)
Do i = 1, pstar % level_size
  Tr  = t % Data( i, ref_lev )
  Zr  = r_theta_levels( i, ref_lev ) - Earth_Radius
  Z0i = orog_interp % Data( i, 1 )     ! interpolated orography
  Z0o = orog_out % Data( i, 1 )        ! output orography

  pstar % Data( i, 1 ) = pstar % Data( i, 1 ) *                       &
                      ( ( Tr + lapse * (Zr - Z0o) ) /                 &
                        ( Tr + lapse * (Zr - Z0i) ) ) ** g_over_lapse_r
End Do


!--------------------------------------------------------------------
! Tidy up
!--------------------------------------------------------------------
Call Rcf_Dealloc_Field( orog_interp )


Return
End Subroutine Rcf_Adjust_Pstar

End Module Rcf_Adjust_Pstar_Mod
