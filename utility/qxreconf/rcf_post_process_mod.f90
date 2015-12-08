! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Post `main loop' processing of data

Module Rcf_Post_Process_Mod

!  Subroutine Rcf_Post_Process - process fields after `main loop'
!
! Description:
! This subroutine performs post-processing of data using transforms
! etc that require a number of output-grid fields. This involves
! some transforms to original data types. Performed after the main
! field data creation loop.
!
! Method:
!  The relevant processing is just run through - this may need to be
!  be made more modular in future.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

Contains

Subroutine Rcf_Post_Process_Atmos( fields_in, field_count_in, orog_in,&
                                   fields_out, field_count_out, orog, &
                                   hdr_in, hdr_out, data_source)

Use Rcf_Field_Equals_Mod, Only : &
    Rcf_Field_Equals

Use Rcf_Exppx_Mod, Only : &
    Rcf_Exppx

Use Submodel_Mod, Only :  &
    Atmos_IM

Use Rcf_Calc_Rho_mod, Only : &
    Rcf_Calc_Rho

Use Rcf_Set_Interp_Flags_Mod, Only  : &
    interp_done,                      &
    interp_no_op

Use Rcf_Freeze_Soil_Mod, Only : &
    Rcf_Freeze_Soil

Use Rcf_Data_Source_Mod, Only : &
    Data_Source_Type,           &
    Ancillary_File,             &
    Input_Dump

Use Rcf_Conv_Cld_Chk_Mod, Only : &
    Rcf_Conv_Cld_Chk

Use Rcf_Cloud_Frac_Chk_Mod, Only : &
    Rcf_Cloud_Frac_Chk

Use Rcf_Snow_Amount_Chk_Mod, Only : &
    Rcf_Snow_Amount_Chk

Use Rcf_Sea_Ice_Frac_Chk_Mod, Only : &
    Rcf_Sea_Ice_Frac_Chk

Use Rcf_Soil_Moist_Chk_Mod, Only : &
    Rcf_Soil_Moist_Chk

Use Rcf_WriteUMhdr_Mod, Only : &
    Rcf_WriteUMhdr

Use Rcf_Field_Type_Mod, Only : &
    field_type

Use Rcf_UMhead_Mod, Only :&
    um_header_type

Use Rcf_Grid_Type_Mod, Only : &
    output_grid

USE decomp_params, ONLY : &
    decomp_rcf_output

Use Rcf_HeadAddress_Mod, Only : &
    RC_PressureTop,         &
    RC_AtmMoist,            &
    RC_AtmMass,             &
    RC_AtmEnergy,           &
    RC_EnergyCorr,          &
    LDC_RHCrit,             &
    FH_GridStagger,         &
    FH_GridStagger_Endgame

Use Rcf_V_Int_Ctl_Mod, Only : &
    v_int_active,             &
    v_int_active_soil

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_active

Use Rcf_Locate_mod, Only : &
    Rcf_Locate

Use Rcf_Alloc_Field_mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

Use Rcf_generate_heights_mod, Only : &
    Rcf_generate_heights

Use rcf_theta_t_convs_mod, Only : &
    Rcf_conv_theta_t,             &
    Rcf_conv_t_theta

Use Rcf_read_field_mod, Only : &
    Rcf_read_field

Use Rcf_Calc_Output_Exner_Mod, Only : &
    Rcf_Calc_Output_Exner

Use Rcf_Write_field_mod, Only : &
    Rcf_write_field

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_Exner_P_Convs_Mod, Only : &
    Rcf_Conv_Exner_P,             &
    Rcf_Conv_P_Exner

USE PrintStatus_mod, Only : &
    LTimer

Use Rcf_Stashcodes_Mod, Only :&
    stashcode_theta,          stashcode_q,       &
    stashcode_rho,            stashcode_exner,        &
    stashcode_p,              stashcode_soil_moist,   &
    stashcode_soil_temp,      stashcode_orog,         &
    stashcode_prog_sec

Use Rcf_Soilstress_To_Soilmoist_Mod,only:&
    Rcf_Soilstress_to_soilmoist

USE rcf_soilconc_to_soilmoist_mod,ONLY:&
    rcf_soilconc_to_soilmoist

Use Rcf_Recon_Mod, Only : &
    use_smc_stress, polar_check

USE Rcf_Polar_Rows_Chk_Mod, ONLY:  &
    Rcf_Polar_Rows_Chk

USE Rcf_Lsm_Mod, Only : &
    local_lsm_out

USE cppxref_mod, ONLY : &
    ppx_atm_tall,       &
    ppx_theta_level,    &
    ppx_rho_level

USE rcf_lsh_land_ice_chk_mod, ONLY : &
    rcf_lsh_land_ice_chk

IMPLICIT NONE

! Arguments
Type( field_type ), Pointer           :: fields_in(:)
Type( field_type ), Pointer           :: fields_out(:)
Type( field_type ), Intent(InOut)     :: orog_in
Type( field_type ), Intent(In)        :: orog
Type( data_source_type ), Pointer     :: data_source(:)
Type( um_header_type ), Intent(In)    :: hdr_in
Type( um_header_type ), Intent(InOut) :: hdr_out
Integer, Intent(In)                   :: field_count_in
Integer, Intent(In)                   :: field_count_out

! Local variables
Integer                       :: pos
Integer                       :: pos_st   !} positions in array of
Integer                       :: pos_sm   !} soil temp and moisture
Integer                       :: pos_orog !} and of orography
Integer                       :: ErrorStatus
Character (Len=*), Parameter  :: RoutineName = 'Rcf_Post_Process'
Character (Len=80)            :: Cmessage

Real                          :: theta_heights(                     &
                                    output_grid % loc_p_field,      &
                                    0 : output_grid % model_levels + 1)
                                    ! Height of theta levels
Real                          :: rho_heights(                       &
                                    output_grid % loc_p_field,      &
                                    0 : output_grid % model_levels + 1)
                                    ! Height of rho levels
Logical                       :: exner_out ! is exner in the output?
Logical                       :: l_soil_change   ! is soil changed

Type( field_type ), Pointer   :: exner
Type( field_type ), Pointer   :: theta
Type( field_type ), Pointer   :: q
Type( field_type ), Pointer   :: rho
Type( field_type ), Pointer   :: p
Type( field_type ), Target    :: pressure  ! exner OR p depending on
                                           ! circumstance

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3)

l_soil_change = .FALSE.     ! initialise soil moisture correction flag 

Call Rcf_Locate( stashcode_prog_sec, stashcode_soil_moist,           &
                 fields_out, field_count_out, pos_sm, .TRUE. )

IF (data_source( pos_sm ) % source == input_dump) THEN
  IF (use_smc_stress) THEN
     CALL rcf_soilstress_to_soilmoist( fields_out, field_count_out,  &
                                     output_grid,decomp_rcf_output,  &
                                     hdr_out)
  ELSE
    IF (v_int_active_soil) THEN
     CALL rcf_soilconc_to_soilmoist( fields_out, field_count_out,    &
                                     output_grid,decomp_rcf_output,  &
                                     hdr_out)
    END IF
  END IF
End If

!-------------------------------------------------------------------
! Find and setup Theta (will hold T if interpolated)
!-------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_theta,                &
                 fields_out, field_count_out, pos)
theta => fields_out( pos )
Call Rcf_Alloc_Field( theta )
Call Rcf_Read_Field( theta, hdr_out, decomp_rcf_output )

If ( h_int_active .OR. v_int_active ) Then
!-------------------------------------------------------------------
! Calculate required heights
!-------------------------------------------------------------------
  Call rcf_generate_heights( output_grid, orog,                &
                             ppx_atm_tall, ppx_theta_level,    &
                             theta_heights, theta % level_size )

  Call rcf_generate_heights( output_grid, orog,               &
                             ppx_atm_tall, ppx_rho_level,     &
                             rho_heights,  theta % level_size )
End If

!--------------------------------------------------------------------
! Find and read Q
!--------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_q,                    &
                 fields_out, field_count_out, pos )
q => fields_out( pos )
Call Rcf_Alloc_Field( q )
Call Rcf_Read_Field( q, hdr_out, decomp_rcf_output )

!-------------------------------------------------------------------
! Find exner and P (as required)
! Note assumption that will have exner *OR* P available
!-------------------------------------------------------------------
Call Rcf_Locate( stashcode_prog_sec, stashcode_exner,                &
                 fields_out, field_count_out, pos, .TRUE. )

If (pos /= 0) Then
  exner => fields_out( pos )
  Call Rcf_Alloc_Field( exner )
  exner_out = .TRUE.

  ! Also need space for P
  p => pressure
  Call Rcf_Field_Equals( p, exner )
  Call Rcf_Alloc_Field( p )
  p % stashmaster => Rcf_Exppx( Atmos_IM, 0, stashcode_p )

Else    ! exner not in output dump

  Call Rcf_Locate( stashcode_prog_sec, stashcode_p,                  &
                   fields_out, field_count_out, pos )
  p => fields_out( pos )
  Call Rcf_Alloc_Field( p )
  exner_out = .FALSE.

  ! Also need space for exner
  exner => pressure
  Call Rcf_Field_Equals( exner, p )
  Call Rcf_Alloc_Field( exner )
  exner % stashmaster => Rcf_Exppx( Atmos_IM, 0, stashcode_exner )

End If

!--------------------------------------------------------------------
! Either read or calculate exner (as appropriate)
!--------------------------------------------------------------------
If ( exner % interp /= interp_no_op ) Then
  If ( exner_out ) Then
    Call Rcf_Read_Field( exner, hdr_out, decomp_rcf_output )
  Else
    Call Rcf_Read_Field( p, hdr_out, decomp_rcf_output )
    exner % Data(:,:) = p % Data(:,:)
    exner % stashmaster => Rcf_Exppx( Atmos_IM, 0, stashcode_p )
    Call Rcf_Conv_P_Exner( exner )
  End If

Else
  If ( theta % interp /= interp_done ) Then
    ErrorStatus = 10
    Cmessage = 'Only have Theta available - need T to calculate exner'
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  Call Rcf_Locate( stashcode_prog_sec, stashcode_orog,               &
                   fields_out, field_count_out, pos_orog )
  Call Rcf_Calc_Output_Exner( fields_in, field_count_in, orog_in,     &
                              hdr_in, orog, theta, q, exner,          &
                              Data_Source( pos_orog ) % source,       &
                              rho_heights, theta_heights )

  If (exner_out) Then
    Call Rcf_Write_Field( exner, hdr_out, decomp_rcf_output )
  End If

  ! Also calculate P - needed internally. Note need to temporarily
  ! fool with stashcode
  p % Data(:,:) = exner % Data(:,:)
  p % stashmaster => Rcf_Exppx( Atmos_IM, 0, stashcode_exner )
  Call Rcf_Conv_Exner_P( p )

  ! Write this newly calcuated field to dump
  If (.NOT. exner_out) Then
    Call Rcf_Write_Field( p, hdr_out, decomp_rcf_output )
  End If

End If

!********************************************************************
! Convert T (in dump) back to theta
!********************************************************************
If ( theta % interp == interp_done ) Then
  Call Rcf_Conv_T_Theta( theta, fields_out, field_count_out, hdr_out, &
                         decomp_rcf_output, rho_heights, theta_heights)

  Call Rcf_Write_Field( theta, hdr_out, decomp_rcf_output )
End If


!*******************************************************************
! Calculate rho
!*******************************************************************
Call Rcf_Locate( stashcode_prog_sec, stashcode_rho,                  &
                 fields_out, field_count_out, pos)
rho => fields_out(pos)
If ( rho % interp == interp_no_op ) Then

  ! No interpolation for rho need to calculate
  Call Rcf_Alloc_Field( rho )

  Call Rcf_Calc_Rho( theta, q, exner, p, theta_heights, rho_heights, &
                     rho)

  Call Rcf_Write_Field( rho, hdr_out, decomp_rcf_output )
  Call Rcf_DeAlloc_Field( rho )

End if

!--------------------------------------------------------------------
! Tidy up
!--------------------------------------------------------------------
Call Rcf_DeAlloc_Field( p )
Call Rcf_DeAlloc_Field( q )
Call Rcf_DeAlloc_Field( exner )
Call Rcf_DeAlloc_Field( theta )

!********************************************************************
! Reset Soil Moisture to match vegetation if interpolated
!********************************************************************
If (h_int_active) Then
  Call Rcf_Soil_Moist_Chk( fields_out, field_count_out, output_grid, &
                           decomp_rcf_output, hdr_out, l_soil_change) 
End If

!********************************************************************
! Check some LSH fields for consistency with land-ice.
!********************************************************************
! We should always call this to make sure any pre-existing dumps are consistent
! even when not interpolating.
CALL rcf_lsh_land_ice_chk( fields_out, field_count_out,              &
                           decomp_rcf_output, hdr_out)


!********************************************************************
! Adjust frozen/unfrozen soil if updated from ancillary
!********************************************************************
Call Rcf_Locate( stashcode_prog_sec, stashcode_soil_temp,            &
                 fields_out, field_count_out, pos_st, .TRUE. )


! Only need to correct if soil moisture/temperature are both in dump
If (pos_st /= 0 .AND. pos_sm /= 0) Then
  If ( Data_Source( pos_st ) % Source == Ancillary_File  .OR. &
       Data_Source( pos_sm ) % Source == Ancillary_File  .OR. &
       l_soil_change ) Then

    Call Rcf_Freeze_Soil( fields_out, field_count_out, hdr_out, &
                          output_grid )
  End If
End If

!********************************************************************
! Check convective cloud base and top are sensible
!********************************************************************
Call Rcf_Conv_Cld_Chk( fields_out, field_count_out, output_grid,  &
                       decomp_rcf_output, hdr_out )

!********************************************************************
! Check that cloud fraction fields are consistent with q fields.
!********************************************************************
Call Rcf_Cloud_Frac_Chk( fields_out, field_count_out, output_grid,  &
                         decomp_rcf_output, hdr_out )

!********************************************************************
! Check for negative snow amounts
!********************************************************************
If (h_int_active) Then
  Call Rcf_Snow_Amount_Chk( fields_out, field_count_out,             &
                            output_grid,  decomp_rcf_output, hdr_out )
End If

!********************************************************************
! Check for small sea-ice fractions if interpolated
!********************************************************************
If (h_int_active) Then
  Call Rcf_Sea_Ice_Frac_Chk( fields_out, field_count_out,            &
                             decomp_rcf_output, hdr_out, data_source)
End If

!********************************************************************
! If appropriate, ensure scalar fields in the polar rows are uniform
!********************************************************************
IF ( output_grid % global .AND. polar_check ) THEN
  IF ( Hdr_Out % Fixhd(FH_GridStagger) == FH_GridStagger_Endgame ) THEN
    ErrorStatus = -20
    WRITE(cmessage,*) 'Polar row checking activated for an invalid grid ',  &
                      '- disabling test'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  ELSE
    CALL Rcf_Polar_Rows_Chk( fields_out, field_count_out, output_grid,      &
                             decomp_rcf_output, hdr_out, local_lsm_out,     &
                             data_source )
  END IF
END IF

!----------------------------------------
! Perform some modifications to headers.  
!----------------------------------------
! This would be better earlier in the code when the headers are setup 
! but there is logic after the headers are setup to turn on interpolation 
! due to change in orography ancillary.  This is the next best place to 
! perform this task.
IF ( .NOT. v_int_active ) THEN
  IF ( .NOT. h_int_active ) THEN
! Want to keep energy correction information if no interpolation 
    hdr_out % realc( rc_atmmoist   ) = hdr_in % realc( rc_atmmoist   )
    hdr_out % realc( rc_atmmass    ) = hdr_in % realc( rc_atmmass    )
    hdr_out % realc( rc_atmenergy  ) = hdr_in % realc( rc_atmenergy  )
    hdr_out % realc( rc_energycorr ) = hdr_in % realc( rc_energycorr )
  END IF
! We want to copy the rhcrit values across if not vertically interpolating.
  hdr_out % levdepc(:,ldc_rhcrit) = hdr_in % levdepc(:,ldc_rhcrit)
!---------------------------------------------------------------
! Write out the header here again due to changes above.
!---------------------------------------------------------------
  CALL rcf_writeumhdr( hdr_out )
END IF

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4)

Return
End Subroutine Rcf_Post_Process_Atmos
End Module Rcf_Post_Process_Mod
