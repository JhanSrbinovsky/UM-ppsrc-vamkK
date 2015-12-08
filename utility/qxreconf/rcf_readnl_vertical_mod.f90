! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ reads the vertical levels namelists

Module rcf_readnl_vertical_Mod
IMPLICIT NONE

!  Subroutine Rcf_Readnl_vertical - reads the vertical levels namelists
!
! Description:
! Module to read in the vertical levels namelists
! VERTICAL, VERTLEVS and JULES_SOIL_PARAM.
! Note that it *must* be called after readnl_recona as
! some sizing is required!
!
! Method:
!  Namelists read in and sets Output_Grid variables appopriately
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


! Some data we need to carry around before it can be assigned
! properly to the output


CONTAINS

Subroutine Rcf_Readnl_Vertical( nft_vertical, nft_vertlevs, nft_shared )

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_order

USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Oper

Use Ereport_Mod, Only : &
    Ereport

USE UM_ParVars, Only : &
    mype

Use Rcf_Generate_Heights_Mod, Only : &
    height_gen_smooth,               &
    height_gen_original

USE vertnamelist_mod, ONLY:  first_constant_r_rho_level,            &
    z_top_of_model, eta_theta, eta_rho, vertlevs

USE Atmos_Max_Sizes, ONLY : model_levels_max

USE missing_data_mod, ONLY : rmdi

USE interpor_mod, ONLY :      &
  interp_order_linear,        &
  interp_order_linear_noex,   &
  interp_order_cubic,         &
  interp_order_quintic

USE read_jules_namelists_mod, ONLY: read_jules_soil_param

USE soil_param, ONLY: dzsoil_io

IMPLICIT NONE

! Subroutine arguments
Integer, Intent(In)          :: nft_vertical  ! File unit - VERTICAL
Integer, Intent(In)          :: nft_vertlevs  ! File unit - VERTLEVS
Integer, Intent(In)          :: nft_shared    ! File unit - JULES_SOIL_PARAM

! local variables/constants
Character (Len=*), Parameter :: RoutineName = 'Rcf_Readnl_Vertical'
Character (Len=80)           :: Cmessage
Integer                      :: ErrorStatus
Integer                      :: levels_theta  ! number of levels of
Integer                      :: levels_rho    ! values read in.
Integer                      :: i             ! looper
Real                         :: last_value    ! used in checking code


! default parameter which is now always used but can be changed if required.
INTEGER, PARAMETER           :: height_method = height_gen_smooth

NAMELIST /VERTICAL/ v_int_order

! Set VERTLEVS defaults
eta_theta(:)          = rmdi
eta_rho(:)            = rmdi

! Set JULES_SOIL_PARAM defaults
dzsoil_io(:)   = 0.0
dzsoil_io(1:4) = (/0.1, 0.25, 0.65, 2.0/)

! Quick error check to make sure parameter model_levels_max isn't too small
If ( model_levels_max < Output_Grid % model_levels ) Then
  ErrorStatus = 10
  Cmessage = 'Internal parameter model_levels_max is too small!'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Read VERTICAL Namelist
Read( Unit = nft_vertical, Nml=vertical )

! Read VERTLEVS Namelist
Read( Unit = nft_vertlevs, Nml=vertlevs )

! Read JULES_SOIL_PARAM Namelist
! We only require dzsoil_io (= soil_depth); the rest is discarded
CALL read_jules_soil_param (nft_shared)

! Write out namelist for diagnostic
If (PrintStatus >= PrStatus_Oper) Then
  If ( mype == 0 ) Then
    Write ( Unit = 6, Nml = vertical)

    Write ( Unit = 6, Fmt = * ) 'Values in VERTLEVS Namelist.'
    Write ( Unit = 6, Fmt = * ) 'z_top_of_model ',z_top_of_model
    Write ( Unit = 6, Fmt = * ) 'first_constant_r_rho_level ',    &
                                 first_constant_r_rho_level
    Write ( Unit = 6, Fmt = * ) 'Eta_Theta'
    Write ( Unit = 6, Fmt = '(5F10.7)' )                          &
          ( eta_theta(i),i=1,Output_Grid % Model_levels+1)
    Write ( Unit = 6, Fmt = * ) 'Eta_Rho'
    Write ( Unit = 6, Fmt = '(5F10.7)' )                          &
          ( eta_rho(i),i=1,Output_Grid % Model_levels)

  End If
End If

! Some sanity checks before we allocate space
If ( Output_Grid % model_levels < 1 ) Then
  Cmessage = 'Model levels for output grid is < 1!'
  ErrorStatus = 20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )

Else If ( Output_Grid % model_levels > 1000 ) Then
  Cmessage = 'Model Levels for output grid is > 1000 - Is this correct?'
  ErrorStatus = -30
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

IF ( v_int_order /= interp_order_linear      .AND. &
     v_int_order /= interp_order_linear_noex .AND. &
     v_int_order /= interp_order_cubic       .AND. &
     v_int_order /= interp_order_quintic ) THEN
  cmessage = 'Vertical Interpolation Order not recognised'
  errorstatus = 40
  CALL ereport( routinename, errorstatus, cmessage )
END IF

! Check that the number of eta values correspond to the number
! of model levels
levels_theta = 0
levels_rho   = 0
Do i = 1, model_levels_max
  If (eta_theta(i) /= RMDI) Then
    levels_theta = levels_theta + 1
  End If

  If (eta_rho(i) /= RMDI) Then
    levels_rho = levels_rho + 1
  End If
End Do

If ( levels_theta /= Output_Grid % model_levels + 1) Then
  Cmessage = 'Mismatch in number of model levels and number of '//&
        'eta_theta values'
  ErrorStatus = 50
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

If ( levels_rho /= Output_Grid % model_levels ) Then
  Cmessage = 'Mismatch in number of model levels and number of '//&
        'eta_rho values'
  ErrorStatus = 60
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Check that eta_rho levels and eta_theta levels are monotone
! ascending.
last_value = RMDI
Do i = 1, Output_Grid % model_levels + 1
  If ( eta_theta(i) <= last_value ) Then
    Cmessage = 'Eta_Theta values are not monotone ascending'
    ErrorStatus = 70
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  last_value = eta_theta(i)
End Do

last_value = RMDI
Do i = 1, Output_Grid % model_levels
  If ( eta_rho(i) <= last_value ) Then
    Cmessage = 'Eta_Rho values are not monotone ascending'
    ErrorStatus = 80
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  last_value = eta_rho(i)
End Do

  If ( height_method /= height_gen_original .AND.   &
       height_method /= height_gen_smooth ) Then
    ErrorStatus = 90
    Cmessage = 'Unrecognised height generation method'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

! Allocate space and fill up relevant parts of Output_Grid
Allocate (Output_Grid % eta_theta_levels                   &
                       ( 0 : Output_Grid % model_levels) )
Allocate ( Output_Grid % eta_rho_levels( Output_Grid % model_levels ) )
Allocate ( Output_Grid % soil_depths( Output_Grid % sm_levels ) )

Output_Grid % eta_theta_levels( 0 : Output_Grid % model_levels ) = &
                     eta_theta( 1 : Output_Grid % model_levels + 1)
Output_Grid % eta_rho_levels( 1 : Output_Grid % model_levels ) =   &
                     eta_rho( 1 : Output_Grid % model_levels )
Output_Grid % soil_depths( 1 : Output_Grid % sm_levels ) =         &
                dzsoil_io( 1 : Output_Grid % sm_levels )
Output_Grid % z_top_of_model = z_top_of_model
Output_Grid % first_constant_r_rho_level = first_constant_r_rho_level
Output_Grid % height_gen_method = height_method


Return
End Subroutine Rcf_Readnl_Vertical
End Module Rcf_Readnl_Vertical_Mod
