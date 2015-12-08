! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Set up Fsat/Fwetl for large-scale hydrology

Module Rcf_Calc_Fsat_Mod

! Subroutine Rcf_Calc_Fsat
!
! Description:
!   Initialises the surface saturation and/or wetland fraction
!
! Method:
!   The fields are calculated from the initial soil parameters and
!   soil moisture content and water table depth.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Reconfiguration

IMPLICIT NONE

Contains

Subroutine Rcf_Calc_Fsat( fields_out, field_count_out, &
                          frac_pos, hdr_out )

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

USE decomp_params, ONLY : &
    decomp_rcf_output

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_HeadAddress_Mod, Only : &
    IC_SoilMoistLevs

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_clapp_hb,        &
    stashcode_ksat,            &
    stashcode_vol_smc_sat,     &
    stashcode_unfrozen_soil,   &
    stashcode_frozen_soil,     &
    stashcode_ti_mean,         &
    stashcode_ti_sig,          &
    stashcode_fexp,            &
    stashcode_gamtot,          &
    stashcode_zw,              &
    stashcode_fsat,            &
    stashcode_fwetl,           &
    stashcode_prog_sec

USE PrintStatus_mod, Only : &
    PrintStatus,                &
    PrStatus_Normal,            &
    LTimer

USE UM_ParVars, Only : &
    mype

Use Rcf_Read_Field_Mod, Only : &
    Rcf_Read_Field

Use Rcf_Alloc_Field_Mod, Only : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

Use Rcf_Grid_Type_Mod, Only : &
    Output_Grid

Use Ereport_mod, Only : &
    Ereport

USE switches, ONLY: l_top

IMPLICIT NONE

! Arguments
Type( field_type ), Pointer          :: fields_out(:)
Type( um_header_type ), Intent( In ) :: hdr_out

Integer, Intent( In )                :: field_count_out
INTEGER, INTENT( IN )                :: frac_pos

! Local variables
Type( field_type ), Pointer          :: frac_field
Type( field_type ), Pointer          :: clapp_hb
Type( field_type ), Pointer          :: ksat
Type( field_type ), Pointer          :: vol_smc_sat
Type( field_type ), Pointer          :: unfrozen_soil
Type( field_type ), Pointer          :: frozen_soil
Type( field_type ), Pointer          :: ti_mean
Type( field_type ), Pointer          :: ti_sig
Type( field_type ), Pointer          :: fexp
Type( field_type ), Pointer          :: gamtot
Type( field_type ), Pointer          :: zw

Integer                              :: soil_index  &
                                       (fields_out(frac_pos) % level_size)
Integer                              :: soil_pts
Integer                              :: size
Integer                              :: pos  ! field position
Integer                              :: i    ! looper
Integer                              :: j    ! looper
Integer                              :: n    ! looper

Real          :: wutot    (fields_out(frac_pos) % level_size)
Real          :: top_crit (fields_out(frac_pos) % level_size)
Real          :: qbase    (fields_out(frac_pos) % level_size)

Real          :: qbase_l (fields_out(frac_pos) % level_size,                    &
                               hdr_out % IntC( IC_SoilMoistLevs ) +1 )
Real          :: Ksz     (fields_out(frac_pos) % level_size,                    &
                               0: hdr_out % IntC( IC_SoilMoistLevs ) )
Real          :: zdepth  (0: hdr_out % IntC( IC_SoilMoistLevs) )

! Save the results if we call this routine again.
Real, Allocatable, Save    :: fsat     (:)
Real, Allocatable, Save    :: fwetl    (:)

Integer       :: ErrorStatus

Character (Len=80)           :: Cmessage
Character (Len=*), Parameter :: RoutineName='Rcf_Calc_Fsat'

Logical :: L_Gamtot=.false.

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Setting up surface saturation/wetland fraction'
End If

If (.NOT.L_TOP) Then
   ErrorStatus = 10
   Write(CMessage, '(A)') 'Field only required in TOPMODEL-based hydrology'
   call Ereport ( RoutineName, ErrorStatus, CMessage)
End If

! Setup pointer to output field.
frac_field => fields_out(frac_pos)

! Only calculate on first call when both fields are not allocated.
IF (.NOT. ALLOCATED(fsat) .AND. .NOT. ALLOCATED(fwetl)) THEN
  ! Allocate the SAVE arrays for results.
  ALLOCATE(fsat(frac_field % level_size))
  ALLOCATE(fwetl(frac_field % level_size))
!----------------------------------------------------------------------
! Find required fields in output dump and read them in:
!----------------------------------------------------------------------
  Call Rcf_Locate( stashcode_prog_sec, stashcode_clapp_hb,           &
                   fields_out, field_count_out, pos)
  clapp_hb => fields_out(pos)
  Call Rcf_Alloc_Field( clapp_hb )
  Call Rcf_Read_Field( clapp_hb, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_ksat,               &
                   fields_out, field_count_out, pos)
  ksat => fields_out(pos)
  Call Rcf_Alloc_Field( ksat )
  Call Rcf_Read_Field( ksat, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_vol_smc_sat,        &
                   fields_out, field_count_out, pos)
  vol_smc_sat => fields_out(pos)
  Call Rcf_Alloc_Field( vol_smc_sat )
  Call Rcf_Read_Field( vol_smc_sat, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_unfrozen_soil,      &
                   fields_out, field_count_out, pos)
  unfrozen_soil => fields_out(pos)
  Call Rcf_Alloc_Field( unfrozen_soil )
  Call Rcf_Read_Field( unfrozen_soil, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_frozen_soil,        &
                   fields_out, field_count_out, pos)
  frozen_soil => fields_out(pos)
  Call Rcf_Alloc_Field( frozen_soil )
  Call Rcf_Read_Field( frozen_soil, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_ti_mean,            &
                   fields_out, field_count_out, pos)
  ti_mean => fields_out(pos)
  Call Rcf_Alloc_Field( ti_mean )
  Call Rcf_Read_Field( ti_mean, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_ti_sig,             &
                   fields_out, field_count_out, pos)
  ti_sig => fields_out(pos)
  Call Rcf_Alloc_Field( ti_sig )
  Call Rcf_Read_Field( ti_sig, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_fexp,               &
                   fields_out, field_count_out, pos)
  fexp => fields_out(pos)
  Call Rcf_Alloc_Field( fexp )
  Call Rcf_Read_Field( fexp, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_gamtot,             &
                   fields_out, field_count_out, pos)
  gamtot => fields_out(pos)
  Call Rcf_Alloc_Field( gamtot )
  Call Rcf_Read_Field( gamtot, hdr_out, decomp_rcf_output )

  Call Rcf_Locate( stashcode_prog_sec, stashcode_zw,                 &
                   fields_out, field_count_out, pos)
  zw => fields_out(pos)
  Call Rcf_Alloc_Field( zw )
  Call Rcf_Read_Field( zw, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Set up soil index:
!----------------------------------------------------------------------
  soil_pts=0
  soil_index(:)=0
  Do i=1 , vol_smc_sat % level_size
    If (vol_smc_sat % Data(i,1) > 0.0) Then
      soil_pts = soil_pts + 1
      soil_index(soil_pts) = i
    End If
  End do

  zdepth(0) = 0.0
  Do n = 1 , unfrozen_soil % levels
    zdepth(n) = zdepth(n-1) + Output_Grid % soil_depths(n)
  End Do

!----------------------------------------------------------------------
! Set up saturated hydraulic conductivity for each soil layer:
!----------------------------------------------------------------------
  Do j = 1 , soil_pts
    i = soil_index(j)
    Do n = 0 , unfrozen_soil % levels
      Ksz(i,n) = ksat % Data(i,1)
    End do
  End do

!----------------------------------------------------------------------
! Initialise variables to zero:
!----------------------------------------------------------------------
  fsat(:) = 0
  fwetl(:) = 0
  qbase(:) = 0
  qbase_l(:,:) = 0
  frac_field % Data (:,:) = 0

!----------------------------------------------------------------------
! Calculate Baseflow:
!----------------------------------------------------------------------
! DEPENDS ON: calc_baseflow
      Call Calc_Baseflow(                    &
                      soil_pts,              &
                      soil_index,            &
                      ti_mean % level_size,  &
                      unfrozen_soil % levels,&
                      zdepth,                &
                      Ksz,                   &
                      clapp_hb % Data,       &
                      fexp % Data,           &
                      ti_mean % Data,        &
                      zw % Data,             &
                      frozen_soil % Data,    &
                      unfrozen_soil % Data,  &
                      wutot,                 &
                      top_crit,              &
                      qbase,                 &
                      qbase_l)

!----------------------------------------------------------------------
! Calculate surface saturation and wetland fraction:
!----------------------------------------------------------------------
! DEPENDS ON: calc_fsat
      Call Calc_Fsat(                        &
                      L_gamtot,              &
                      soil_pts,              &
                      soil_index,            &
                      ti_mean % level_size,  &
                      ti_mean % Data,        &
                      ti_sig % Data,         &
                      wutot,                 &
                      top_crit,              &
                      gamtot % Data,         &
                      fsat,                  &
                      fwetl)


!----------------------------------------------------------------------
! Tidy Up
!----------------------------------------------------------------------
  Call Rcf_Dealloc_Field( clapp_hb )
  Call Rcf_Dealloc_Field( ksat )
  Call Rcf_Dealloc_Field( vol_smc_sat )
  Call Rcf_Dealloc_Field( unfrozen_soil )
  Call Rcf_Dealloc_Field( frozen_soil )
  Call Rcf_Dealloc_Field( ti_mean )
  Call Rcf_Dealloc_Field( ti_sig )
  Call Rcf_Dealloc_Field( fexp )
  Call Rcf_Dealloc_Field( gamtot )
  Call Rcf_Dealloc_Field( zw )

END IF

!----------------------------------------------------------------------
! Set up output to surface saturation or wetland fraction:
!----------------------------------------------------------------------
SELECT CASE(frac_field % stashmaster % item)
CASE(stashcode_fsat)
  frac_field % Data(:,1)=fsat(:)
  DEALLOCATE(fsat)
CASE(stashcode_fwetl)
  frac_field % Data(:,1)=fwetl(:)
  DEALLOCATE(fwetl)
END SELECT

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4)
RETURN
End Subroutine Rcf_Calc_Fsat

End Module Rcf_Calc_Fsat_Mod
