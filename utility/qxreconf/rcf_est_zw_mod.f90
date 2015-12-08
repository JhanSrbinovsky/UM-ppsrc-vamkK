! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Set up Zw for large-scale hydrology

Module Rcf_Est_Zw_Mod

! Subroutine Rcf_Est_Zw
!
! Description:
!   Initialises gridbox mean water table depth.
!
! Method:
!   Estimates the gridbox mean water table depth by assuming
!   that the soil drainage and baseflow balance, i.e:
!   Ksat*THETA^C=Ksat/f*exp(-f*Z-Ti)
!   The estimate of the exact  water table depth if it is within the
!   soil model is not important as it is diagnosed at this depth anyway.
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Reconfiguration

IMPLICIT NONE

Contains

Subroutine Rcf_Est_Zw( fields_out, field_count_out, &
                             zw, hdr_out )

Use Rcf_Locate_Mod, Only : &
    Rcf_Locate

USE decomp_params, ONLY : &
    decomp_rcf_output

Use Rcf_Field_Type_Mod, Only : &
    Field_Type

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Stashcodes_Mod, Only : &
    stashcode_clapp_hb,        &
    stashcode_vol_smc_sat,     &
    stashcode_unfrozen_soil,   &
    stashcode_frozen_soil,     &
    stashcode_zw,              &
    stashcode_ti_mean,         &
    stashcode_fexp,            &
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
Type( field_type ), Intent( InOut )  :: zw
Type( um_header_type ), Intent( In ) :: hdr_out

Integer, Intent( In )                :: field_count_out

! Local variables
Type( field_type ), Pointer          :: clapp_hb
Type( field_type ), Pointer          :: vol_smc_sat
Type( field_type ), Pointer          :: unfrozen_soil
Type( field_type ), Pointer          :: frozen_soil
Type( field_type ), Pointer          :: ti_mean
Type( field_type ), Pointer          :: fexp

Integer                              :: soil_index  &
                                       (zw % level_size)
Integer                              :: soil_pts

Integer                              :: size
Integer                              :: pos  ! field position
Integer                              :: i    ! looper
Integer                              :: n    ! looper
Integer                              :: j    ! looper

Integer      :: ErrorStatus

Character (Len=80)           :: Cmessage
Character (Len=*), Parameter :: RoutineName='Rcf_Est_Zw'

Real                                 :: fnstuf
Real                                 :: temp_fn
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
If (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) Then
  Write (6,*) 'Estimating the initial water table depth'
End If

If (.NOT.L_TOP) Then
   ErrorStatus = 10
   Write(CMessage, '(A)') 'Field only required in TOPMODEL-based hydrology'
   call Ereport ( RoutineName, ErrorStatus, CMessage)
End If

!----------------------------------------------------------------------
! Find required fields in output dump and read them in:
!----------------------------------------------------------------------

  Call Rcf_Locate( stashcode_prog_sec, stashcode_clapp_hb,           &
                   fields_out, field_count_out, pos)
  clapp_hb => fields_out(pos)
  Call Rcf_Alloc_Field( clapp_hb )
  Call Rcf_Read_Field( clapp_hb, hdr_out, decomp_rcf_output )

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

  Call Rcf_Locate( stashcode_prog_sec, stashcode_fexp,               &
                   fields_out, field_count_out, pos)
  fexp => fields_out(pos)
  Call Rcf_Alloc_Field( fexp )
  Call Rcf_Read_Field( fexp, hdr_out, decomp_rcf_output )

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

  Soil_depth=0.0
  Do n = 1,unfrozen_soil % levels
    soil_depth=soil_depth+Output_Grid % soil_depths(n)
  End do

!----------------------------------------------------------------------
! Estimate zw by assuming equilibruim:
!----------------------------------------------------------------------
      zw % Data(:,:) = Zw_Max

      Do j = 1, soil_pts
        i = soil_index(j)

        If (fexp % Data(i,1) <= 0.0) Then
           ErrorStatus = 20
           write (CMessage, '(A)')                                   &
             'Exponent. decay in Vol_Smc_Sat wrongly/not set'
           call Ereport ( RoutineName, ErrorStatus, CMessage)
        End If

        If (frozen_soil % Data(i,unfrozen_soil % levels) < 1.0) Then

          fnstuf = unfrozen_soil % Data(i,unfrozen_soil % levels) /  &
            (1.0 - frozen_soil % Data(i,unfrozen_soil % levels))

          temp_fn = fnstuf**(2. * clapp_hb % Data(i,1) + 3.0) *      &
            fexp % Data(i,1) * exp(ti_mean % Data(i,1)) +            &
            exp(-fexp % Data(i,1) * (Zw_Max - Soil_Depth))

          If (temp_fn > 0.0) Then
            zw % Data(I,1) =                        &
            -1./Fexp % Data(i,1) * ALOG(Temp_Fn) + Soil_Depth
          End If

          If (zw % Data(I,1) > Zw_Max) Then
            Zw % Data(I,1) = Zw_Max
          End If
          If (zw % Data(I,1) < 0.0) Then
            Zw % Data(I,1) = 0.0
          End If

        End If

      End do

!----------------------------------------------------------------------
! Tidy Up
!----------------------------------------------------------------------
  Call Rcf_Dealloc_Field( clapp_hb )
  Call Rcf_Dealloc_Field( vol_smc_sat )
  Call Rcf_Dealloc_Field( unfrozen_soil )
  Call Rcf_Dealloc_Field( frozen_soil )
  Call Rcf_Dealloc_Field( ti_mean )
  Call Rcf_Dealloc_Field( fexp )

! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 4)

RETURN
End Subroutine Rcf_Est_Zw

End Module Rcf_Est_Zw_Mod
