! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initial setting of interpolation logical switches.

Module Rcf_Set_Interp_Logicals_Mod
IMPLICIT NONE

!  Subroutine Rcf_Set_Interp_Logicals
!
! Description:
!    Sets the h_int_active, v_int_active and v_int_active_soil
!    logical switches that control horizontal and vertical
!    interpolation. This is just an initial setting and may be
!    changed at a later time.
!
! Method:
!    h_int_active - Checks for domain and row/row-length changes
!    v_int_active - Checks for # of model/wet level changes and
!                   changes to eta values or height field constants
!    v_int_active_soil - Checks for changes in number and values of
!                        soil levels.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

PRIVATE :: rcf_check_horiz_change

Contains

Subroutine Rcf_Set_Interp_Logicals( Input_Grid, Output_Grid, Hdr_In)

Use Rcf_Grid_Type_Mod, Only : &
    Grid_Type

Use Rcf_UMhead_Mod, Only : &
    UM_Header_Type

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active,                  &
    h_int_active_u,                &
    h_int_active_v

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_active,         &
    v_int_active_soil

USE PrintStatus_mod, Only : &
    PrintStatus,                &
    PrStatus_Normal,            &
    PrStatus_Oper

Use Rcf_HeadAddress_Mod, Only :                                &
    RC_LongSpacing,        RC_LatSpacing,         RC_FirstLat, &
    RC_FirstLong,          RC_PoleLat,            RC_PoleLong

USE rcf_model_mod, ONLY : &
    frstlona, frstlata,   &
    polelona, polelata

Use Rcf_Generate_Heights_Mod, Only : &
    height_gen_original

USE UM_ParVars, Only : &
    mype

Implicit None

! Arguments
Type( Grid_Type ), Intent(In)       :: Input_Grid
Type( Grid_Type ), Intent(In)       :: Output_Grid
Type( UM_Header_Type ), Intent(In)  :: Hdr_In

! Local variables
Integer                             :: i        ! Looper
Character (Len=40)                  :: reason   ! why is interp.
                                                ! switched on.

!--------------------------------------------------------------------
! Horizontal Interpolation
!--------------------------------------------------------------------
! Default for all gridtypes.
h_int_active = .FALSE.
h_int_active_u = .FALSE.
h_int_active_v = .FALSE.

! First check P grid.
CALL rcf_check_horiz_change(input_grid, output_grid, "p grid", h_int_active)

! If P grid has changed we assume all gridtypes must change.
IF (h_int_active) THEN
  WRITE (6,'(A)') 'Horizontal interpolation is switched ON for u and v grid.'
  IF (PrintStatus >= PrStatus_Oper) THEN
    WRITE (6,'(A)') 'Because of p grid being interpolated.'
  END IF
  h_int_active_u = .TRUE.
  h_int_active_v = .TRUE.
ELSE
! If P grid has not changed we can check the other grids to make sure
! staggering has not changed.  ND LAM -> EG LAM makes this possible.
  CALL rcf_check_horiz_change(input_grid, output_grid, "u grid", h_int_active_u)
  CALL rcf_check_horiz_change(input_grid, output_grid, "v grid", h_int_active_v)
END IF

!---------------------------------------------------------------------
! Vertical interpolation
!---------------------------------------------------------------------
v_int_active = .FALSE.

! Turn on if number of levels varies
If (Input_Grid % model_levels /= Output_Grid % model_levels .OR.    &
    Input_Grid % wet_levels /= Output_Grid % wet_levels) Then

  v_int_active = .TRUE.
  reason = 'changed number of model or wet levels'

Else   ! Turn on if eta levels vary
  Do i = 1, Hdr_In % Len1LevDepC
    If ( Abs( Hdr_In % LevDepC( i, 1 ) -                            &
         Output_Grid % eta_theta_levels( i-1)) > Epsilon( 1.0 )) Then

      v_int_active = .TRUE.
      reason = 'change in eta_theta levels'
     End If
  End Do

  Do i = 1, Hdr_In % Len1LevDepC - 1
    If ( Abs( Hdr_In % LevDepC( i, 2 ) -                            &
         Output_Grid % eta_rho_levels( i ) ) > Epsilon( 1.0 ) ) Then

      v_int_active = .TRUE.
      reason = 'change in eta_rho levels'
    End If
  End Do
End If

! Turn on if height specifying constants vary
If ( (Input_Grid  % first_constant_r_rho_level /=     &
      Output_Grid % first_constant_r_rho_level ) .OR. &
     (Input_Grid % z_top_of_model - Output_Grid % z_top_of_model ) > &
                                    Epsilon( 1.0 ) ) Then
  v_int_active = .TRUE.
  reason = 'change in height defining constants'
End If

! Turn on if number of boundary levels change
If ( Input_Grid % Height_Gen_Method == height_gen_original .AND. &
     Input_Grid % BL_Levels /= Output_Grid % BL_Levels ) Then
  v_int_active = .TRUE.
  reason = 'change in number of boundary levels'
End If

! Turn on if change in way heights are calculated
If ( Input_Grid  % Height_Gen_Method /=     &
     Output_Grid % Height_Gen_Method ) Then
  v_int_active = .TRUE.
  reason = 'change in method of height generation'
End If

If (v_int_active) Then
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
    Write (6,*) 'Vertical interpolation is switched ON'
    If (PrintStatus >= PrStatus_Oper) Then
      Write (6,*) 'Because of ', reason
    End If
  End If
Else
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
    Write (6,*) 'Vertical interpolation is switched OFF'
  End If
End If

!---------------------------------------------------------------------
! Vertical interpolation for soil levels
!---------------------------------------------------------------------
v_int_active_soil = .FALSE.
If ( Input_Grid % sm_levels /= Output_Grid % sm_levels) Then
  v_int_active_soil = .TRUE.
  reason = 'changed number of soil depths'
Else
  Do i = 1, Input_Grid % sm_levels
    If ( Abs( Input_Grid % soil_depths(i) -                          &
              Output_Grid % soil_depths(i) ) > Epsilon(1.0) ) Then

      v_int_active_soil = .TRUE.
      reason = 'changed soil depth'
    End If
  End Do
End If

If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  If (v_int_active_soil) Then
    Write (6,*) 'Vertical interpolation between soil depths is '//&
                 'switched ON'
    If (PrintStatus >= PrStatus_Oper) Then
      Write (6,*) 'Because of ', reason
    End If
  Else
    Write (6,*) 'Vertical interpolation between soil depths is '//&
                  'switched OFF'
  End If
    Write (6,*)
End If

Return
End Subroutine Rcf_set_interp_logicals

!  Subroutine Rcf_check_horiz_change
!
! Description:
!   Given a input and output grid and a gridtype to check it will set a flag
!   corresponding to whether the grid has changed.
!
SUBROUTINE rcf_check_horiz_change(input_grid, output_grid, gridtype, flag)

Use Rcf_Grid_Type_Mod, Only : &
    Grid_Type
USE UM_ParVars, Only : &
    mype
USE PrintStatus_mod, Only : &
    PrintStatus,                &
    PrStatus_Normal,            &
    PrStatus_Oper

IMPLICIT NONE
! Arguments
Type( Grid_Type ), Intent(In)       :: Input_Grid
Type( Grid_Type ), Intent(In)       :: Output_Grid
CHARACTER (LEN=*)                   :: gridtype
LOGICAL          , INTENT(OUT)      :: flag

REAL, POINTER :: lambda_in (:)
REAL, POINTER :: phi_in    (:)
REAL, POINTER :: lambda_out(:)
REAL, POINTER :: phi_out   (:)
INTEGER       :: row_length_in
INTEGER       :: rows_in
INTEGER       :: row_length_out
INTEGER       :: rows_out
INTEGER       :: i
CHARACTER (LEN=60)                  :: reason   ! why is interp.
                                                ! switched on.

flag = .FALSE.

IF (gridtype == "p grid") THEN
  row_length_in  =  input_grid  % glob_p_row_length
  rows_in        =  input_grid  % glob_p_rows
  phi_in         => input_grid  % phi_p
  lambda_in      => input_grid  % lambda_p
  row_length_out =  output_grid % glob_p_row_length
  rows_out       =  output_grid % glob_p_rows
  phi_out        => output_grid % phi_p
  lambda_out     => output_grid % lambda_p
ELSE IF(gridtype == "u grid") THEN
  row_length_in  =  input_grid  % glob_u_row_length
  rows_in        =  input_grid  % glob_u_rows
  phi_in         => input_grid  % phi_p
  lambda_in      => input_grid  % lambda_u
  row_length_out =  output_grid % glob_u_row_length
  rows_out       =  output_grid % glob_u_rows
  phi_out        => output_grid % phi_p
  lambda_out     => output_grid % lambda_u
ELSE IF(gridtype == "v grid") THEN
  row_length_in  =  input_grid  % glob_v_row_length
  rows_in        =  input_grid  % glob_v_rows
  phi_in         => input_grid  % phi_v
  lambda_in      => input_grid  % lambda_p
  row_length_out =  output_grid % glob_v_row_length
  rows_out       =  output_grid % glob_v_rows
  phi_out        => output_grid % phi_v
  lambda_out     => output_grid % lambda_p
END IF


IF ( .NOT. Output_Grid % Global ) THEN
  IF ( ABS( phi_in(1) - phi_out(1) ) > EPSILON( 1.0 ) ) THEN
    flag = .TRUE.
    reason = 'differing first latitude for ' // gridtype
  END IF
END IF

IF ( ABS( lambda_in(1) - lambda_out(1) ) > EPSILON( 1.0 ) .OR. &
     ABS( input_grid % lambda_npole - output_grid % lambda_npole ) >       &
                                            EPSILON( 1.0 ) .OR. &
     ABS( input_grid % phi_npole    - output_grid % phi_npole )    >       &
                                            EPSILON( 1.0 ) ) THEN
  flag = .TRUE.
  reason = 'changed domain for ' // gridtype
END IF

! Check horizontal grid resolution.
IF ( rows_in       /= rows_out       .OR. & 
     row_length_in /= row_length_out) THEN
  flag = .TRUE.
  reason = 'changed horizontal resolution for ' // gridtype
ELSE
! If the input and output grids have the same number of rows and points on rows
! we can compare the actual grid coordinates due to variable resolution is
! possible.
  DO i = 1, rows_in
    IF (ABS(phi_in(i) - phi_out(i)) > EPSILON(1.0)) THEN
      flag = .TRUE.
      reason = 'changed phi grid information for ' // gridtype
      EXIT
    END IF
  END DO
  
  DO i = 1, row_length_in
    IF (ABS(lambda_in(i) - lambda_out(i)) > EPSILON(1.0)) THEN
      flag = .TRUE.
      reason = 'changed lambda grid information for ' // gridtype
      EXIT
    END IF
  END DO
END IF

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE (6,'(A)') ''
  IF (flag) THEN
    WRITE (6,'(A)') 'Horizontal interpolation is switched ON for ' // gridtype
    IF (PrintStatus >= PrStatus_Oper) THEN
      WRITE (6,'(A)') 'Because of ', reason
    END IF
  ELSE
    WRITE (6,'(A)') 'Horizontal interpolation is switched OFF for ' // gridtype
  END IF
END IF


END SUBROUTINE rcf_check_horiz_change
End Module Rcf_Set_Interp_Logicals_Mod
