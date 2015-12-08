! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Sets up Fitting parameters for large-scale hydrology

Module Rcf_Fit_Fsat_Mod

! Subroutine Rcf_Fit_Fsat
!
! Description:
!   Calls Calc_Fit_Fsat which calculates the fitting parameters for
!   LSH model and stores them.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Reconfiguration

IMPLICIT NONE

Contains

Subroutine Rcf_Fit_Fsat( fields_out, field_count_out,        &
                         field_pos, hdr_out )

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

USE decomp_params, ONLY : &
    decomp_rcf_output

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Stashcodes_Mod, Only :  &
    stashcode_sthzw,            &
    stashcode_vol_smc_sat,      &
    stashcode_unfrozen_soil,    &
    stashcode_ti_mean,          &
    stashcode_ti_sig,           &
    stashcode_gamtot,           &
    stashcode_fexp,             &
    stashcode_a_fsat,           &
    stashcode_c_fsat,           &
    stashcode_a_fwet,           &
    stashcode_c_fwet,           &
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

INTEGER, INTENT( IN )                :: field_pos
Integer, Intent( In )                :: field_count_out

! Local variables
Type( field_type ), Pointer          :: vol_smc_sat
Type( field_type ), Pointer          :: unfrozen_soil
Type( field_type ), Pointer          :: ti_mean
Type( field_type ), Pointer          :: ti_sig
Type( field_type ), Pointer          :: gamtot
Type( field_type ), Pointer          :: fexp

Integer                              :: soil_index  &
                                       (fields_out(field_pos) % level_size)
Integer                              :: soil_pts

Integer                              :: pos  ! field position
Integer                              :: i    ! looper
Integer                              :: n    ! looper
Integer                              :: j    ! looper

! These fields will be saved
Real, Allocatable, Save  :: a_fsat (:,:)
Real, Allocatable, Save  :: c_fsat (:,:)
Real, Allocatable, Save  :: a_fwet (:,:)
Real, Allocatable, Save  :: c_fwet (:,:)

Integer      :: ErrorStatus

Character (Len=80)           :: Cmessage
Character (Len=*), Parameter :: RoutineName='Rcf_Fit_Fsat'

Real                                 :: soil_depth

! Comdecks:
! C_TOPOG start
!
! Topographic index increment:
      REAL,PARAMETER :: DTI = 0.2
! Maximum topographic index considered:
      REAL,PARAMETER :: TI_MAX = 10.0
! Maximum allowed water table depth (m):
      REAL,PARAMETER :: ZW_MAX = 6.0
! Standard deviation of LOG(Ksat(0)):
      REAL,PARAMETER :: SIGMA_LOGK = 0.0
! Parameter to remove very high water tables
! from the calculated wetland fraction:
      REAL,PARAMETER :: TI_WETL = 1.5
!
! C_TOPOG end

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
If (.NOT.L_TOP) Then
  ErrorStatus = 30
  Write(CMessage, '(A)') 'Field only required in TOPMODEL-based hydrology'
  Call Ereport ( RoutineName, ErrorStatus, CMessage)
End If

! Only calculate when all fields are not allocated.
IF (.NOT. ALLOCATED(a_fsat) .AND. &
    .NOT. ALLOCATED(c_fsat) .AND. &
    .NOT. ALLOCATED(a_fwet) .AND. &
    .NOT. ALLOCATED(c_fwet) ) THEN
  
  If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
    Write (6,*) 'Setting up call for calc_fit_fsat in rcf_fit_fsat'
  End If

! Allocate all fields.
  ALLOCATE(a_fsat( fields_out(field_pos) % level_size, &
                   fields_out(field_pos) % levels ) )
  ALLOCATE(c_fsat( fields_out(field_pos) % level_size, &
                   fields_out(field_pos) % levels ) )
  ALLOCATE(a_fwet( fields_out(field_pos) % level_size, &
                   fields_out(field_pos) % levels ) )
  ALLOCATE(c_fwet( fields_out(field_pos) % level_size, &
                   fields_out(field_pos) % levels ) )
  
!----------------------------------------------------------------------
! Find required fields in output dump and read them in:
!----------------------------------------------------------------------
  Call Rcf_Locate( stashcode_prog_sec, stashcode_vol_smc_sat,            &
                   fields_out, field_count_out, pos)
  vol_smc_sat => fields_out(pos)
  Call Rcf_Alloc_Field( vol_smc_sat )
  Call Rcf_Read_Field( vol_smc_sat, hdr_out, decomp_rcf_output )
  
  Call Rcf_Locate( stashcode_prog_sec, stashcode_unfrozen_soil,          &
                   fields_out, field_count_out, pos)
  unfrozen_soil => fields_out(pos)
  Call Rcf_Alloc_Field( unfrozen_soil )
  Call Rcf_Read_Field( unfrozen_soil, hdr_out, decomp_rcf_output )
  
  Call Rcf_Locate( stashcode_prog_sec, stashcode_ti_mean,                &
                   fields_out, field_count_out, pos)
  ti_mean => fields_out(pos)
  Call Rcf_Alloc_Field( ti_mean )
  Call Rcf_Read_Field( ti_mean, hdr_out, decomp_rcf_output )
  
  Call Rcf_Locate( stashcode_prog_sec, stashcode_ti_sig,                 &
                   fields_out, field_count_out, pos)
  ti_sig => fields_out(pos)
  Call Rcf_Alloc_Field( ti_sig )
  Call Rcf_Read_Field( ti_sig, hdr_out, decomp_rcf_output )
  
  Call Rcf_Locate( stashcode_prog_sec, stashcode_gamtot,                 &
                   fields_out, field_count_out, pos)
  gamtot => fields_out(pos)
  Call Rcf_Alloc_Field( gamtot )
  Call Rcf_Read_Field( gamtot, hdr_out, decomp_rcf_output )
  
  Call Rcf_Locate( stashcode_prog_sec, stashcode_fexp,                   &
                   fields_out, field_count_out, pos)
  fexp => fields_out(pos)
  Call Rcf_Alloc_Field( fexp )
  Call Rcf_Read_Field( fexp, hdr_out, decomp_rcf_output )
  
!----------------------------------------------------------------------
! Initialise fitting variables to zero:
!----------------------------------------------------------------------
  a_fsat(:,:) = 0.0
  c_fsat(:,:) = 0.0
  a_fwet(:,:) = 0.0
  c_fwet(:,:) = 0.0
  
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
  End Do
  
  Soil_depth=0.0
  Do n = 1,unfrozen_soil % levels
    soil_depth=soil_depth+Output_Grid % soil_depths(n)
  End Do
  
!----------------------------------------------------------------------
! Set up variables needed for call to calc_fit_fsat
!----------------------------------------------------------------------
  
  Do j = 1, soil_pts
    i = soil_index(j)
    If (gamtot % Data(i,1) <= 0.0) Then
      ErrorStatus = 20
      Write (CMessage, '(A)')'Gamma integral should be positive and non-zero.'
      Call Ereport ( RoutineName, ErrorStatus, CMessage)
    End If
 
    If (fexp % Data(i,1) <= 0.0) Then
      ErrorStatus = 21
      Write (CMessage, '(A)')'Exp. decay in Vol_Smc_Sat wrongly/not set'
      Call Ereport ( RoutineName, ErrorStatus, CMessage)
    End If
  
  End Do
  
!----------------------------------------------------------------------
! Get fitting parameters:
!----------------------------------------------------------------------
  
! DEPENDS ON: calc_fit_fsat
  Call Calc_Fit_Fsat(                       &
                     soil_pts,              &
                     soil_index,            &
                     ti_mean % level_size,  &
                     fexp % Data(:,1),      &
                     ti_mean % Data(:,1),   &
                     ti_sig % Data(:,1),    &
                     gamtot % Data(:,1),    &
                     soil_depth,            &
                     a_fsat,                &
                     c_fsat,                &
                     a_fwet,                &
                     c_fwet)
  
!----------------------------------------------------------------------
! Tidy Up
!----------------------------------------------------------------------
  Call Rcf_Dealloc_Field( vol_smc_sat )
  Call Rcf_Dealloc_Field( unfrozen_soil )
  Call Rcf_Dealloc_Field( ti_mean )
  Call Rcf_Dealloc_Field( ti_sig )
  Call Rcf_Dealloc_Field( gamtot )
  Call Rcf_Dealloc_Field( fexp )
END IF

Select Case( fields_out( field_pos ) % stashmaster % item )

  Case( stashcode_a_fsat)
    fields_out(field_pos) % Data (:,:) = a_fsat (:,:)
    DEALLOCATE(a_fsat)

  Case( stashcode_c_fsat )
    fields_out(field_pos) % Data (:,:) = c_fsat (:,:)
    DEALLOCATE(c_fsat)

  Case( stashcode_a_fwet)
    fields_out(field_pos) % Data (:,:) = a_fwet (:,:)
    DEALLOCATE(a_fwet)

  Case( stashcode_c_fwet )
    fields_out(field_pos) % Data (:,:) = c_fwet (:,:)
    DEALLOCATE(c_fwet)

End Select

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4)

RETURN
End Subroutine Rcf_Fit_Fsat

End Module Rcf_Fit_Fsat_Mod
