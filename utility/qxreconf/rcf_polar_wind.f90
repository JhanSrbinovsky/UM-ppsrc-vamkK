! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Averaging of Polar wind rows

Module Rcf_Polar_Wind_Mod

!  Subroutine Rcf_Polar_Wind - Averaging of Polar wind rows
!
! Description:
!   Performs a vector mean on Polar wind rows.
!
! Method:
!   Magnitude and direction of non-polar rows are calculated and used to
!   set adjecent polar rows.  Normally this is using v to estimate u at pole
!   but ENDGAME requires v at poles so uses u to calculate v.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

Subroutine Rcf_Polar_Wind( polar_wind, non_polar_wind, delta_lon)

Use Rcf_Field_Type_Mod, Only : &
field_type

USE UM_ParVars, Only : &
mype,               &
atNorth,            &
atSouth,            &
datastart,          &
gc_proc_row_group

USE PrintStatus_mod, Only : &
PrintStatus,            &
PrStatus_Normal

USE conversions_mod, ONLY : &
pi_over_180

Implicit None

! Arguments
Type( field_type), Intent(InOut)  :: polar_wind
Type( field_type), Intent(In)     :: non_polar_wind
Real,              Intent(In)     :: delta_lon

! Local Variables.

Integer      :: i
Integer      :: k
Integer      :: info
Integer      :: npw_start_pos, pw_start_pos, pole_dir
Integer      :: num_pole, pole

Real         :: a_p
Real         :: b_p
Real         :: longitude
Real         :: mag_vector
Real         :: dir_vector

Real         :: wrk( 2 * non_polar_wind % levels)
Real         :: rwrk(non_polar_wind % row_len, 2 * non_polar_wind % levels)
Real         :: delta_lon_rad

External gcg_rvecsumr


If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Finding vector mean for polar wind rows'
End If

delta_lon_rad = pi_over_180 * delta_lon
! ----------------------------------------------------------------------
! Section 1.   Calculate magnitude and direction of polar vector wind
!              from non-polar component of wind on row around pole.
! ----------------------------------------------------------------------

! Only need data at poles

If (atNorth .AND. atSouth) Then
  num_pole = 2
Else
  num_pole = 1
End If

If (atNorth .OR. atSouth) Then
  Do pole = 1, num_pole

! Calculate start index in data depending on which pole we are.
! Pole_dir gives the sign of the result (cos(t)+sin(t)) after 
! rotating a vector.
    If (atNorth .AND. pole == 1 ) Then
      npw_start_pos = (non_polar_wind % rows - 1) * non_polar_wind % row_len
      pw_start_pos = (polar_wind % rows - 1) * polar_wind % row_len
      pole_dir = 1.0
    Else If (atSouth) Then
      npw_start_pos = 0
      pw_start_pos = 0
      pole_dir = -1.0
    End If

    Do k = 1, 2 * non_polar_wind % levels
      wrk(k) = 0.
    End Do

    Do k = 1, non_polar_wind % levels
      Do i = 1, non_polar_wind % row_len
! This strange way of writing things is to fox the Cray
! optimiser and thus retain bit comparison (compiler bug)
        longitude = (datastart(1) + i - 2) * delta_lon_rad
        rwrk(i,2*(k-1)+1) = non_polar_wind % Data( npw_start_pos + i, k) * &
                            cos( longitude )
        rwrk(i,2*k)       = non_polar_wind % Data( npw_start_pos + i, k) * &
                            sin( longitude )
      End Do
    End Do

    Call gcg_rvecsumr(non_polar_wind % row_len, non_polar_wind % row_len, &
                      1, 2 * non_polar_wind % levels, rwrk,               &
                      gc_proc_row_group, info, wrk)

! If estimating v using u so rotation is needed in the opposite
! direction.
    If (polar_wind % stashmaster % grid_type == 19) Then
      pole_dir = -1.0 * pole_dir
    End If

    Do k = 1, non_polar_wind % levels

      a_p = 2. * wrk(2*(k-1)+1) / non_polar_wind % glob_row_len
      b_p = 2. * wrk(2*k) / non_polar_wind % glob_row_len

      mag_vector = sqrt(a_p*a_p + b_p*b_p)

      If (a_p .eq. 0. .and. b_p .eq. 0.) Then
        dir_vector = 0.
      Else
        dir_vector = atan2 (b_p, a_p)
      End If

      Do i = 1, polar_wind % row_len
        polar_wind % Data( pw_start_pos + i, k) = pole_dir *   &
          ( mag_vector * sin( (datastart(1) + i - 1 - .5) *    &
          delta_lon_rad - dir_vector ) )
      End Do
    End Do

  End Do
End If

Return
End Subroutine Rcf_Polar_Wind

End Module Rcf_Polar_Wind_Mod
