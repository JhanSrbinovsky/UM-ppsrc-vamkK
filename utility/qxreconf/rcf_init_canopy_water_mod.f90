! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialisations for tiled canopy water field

MODULE Rcf_Init_Canopy_Water_Mod
IMPLICIT NONE

! Subroutine Rcf_Init_Canopy_Water
!
! Description:
!   Initialises canopy water on tiles.
!
! Method:
!   Initialises canopy water on tiles from mean canopy water.  This
!   latter field is only found in the input dump, so needs to be
!   interpolated.
!
!   Code Owner: See Unified Model Code Owners HTML page
!   This file belongs in section: Reconfiguration

CONTAINS

SUBROUTINE Rcf_Init_Canopy_Water( fields_in,  field_count_in,  hdr_in, &
                                  canopy_water)

USE Rcf_Field_Type_Mod, ONLY : &
    Field_Type

USE Rcf_UMhead_Mod, ONLY : &
    um_header_type

USE Rcf_Stashcodes_Mod, ONLY : &
    stashcode_mean_canopyw,    &
    stashcode_prog_sec

USE PrintStatus_mod, ONLY : &
    PrintStatus,                &
    PrStatus_Normal

USE UM_ParVars, ONLY : &
    mype

USE Rcf_Read_Field_Mod, ONLY : &
    Rcf_Read_Field

USE Rcf_Interpolate_Mod, ONLY : &
    Rcf_Interpolate

USE Rcf_Field_Equals_Mod, ONLY : &
    Rcf_Field_Equals

USE decomp_params, ONLY : &
    decomp_rcf_input

USE Rcf_Grid_Type_Mod, ONLY : &
    Input_Grid,               &
    Output_Grid

USE Rcf_Alloc_Field_Mod, ONLY : &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Interp_Weights_Mod, ONLY : &
    h_int_active

USE Rcf_Set_Interp_Flags_Mod, ONLY : &
    interp_h_only,                   &
    interp_copy

USE Rcf_Locate_Mod, ONLY : &
    Rcf_Locate

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER          :: fields_in(:)
TYPE( field_type ), INTENT( InOut )  :: canopy_water
TYPE( um_header_type ), INTENT( In ) :: hdr_in

INTEGER, INTENT( In )                :: field_count_in

! Local variables
TYPE (field_type), POINTER           :: mean_canopyw_in
TYPE (field_type)                    :: mean_canopyw_out
TYPE (field_type)                    :: dummy  ! pretend orography
INTEGER                              :: pos    ! field position
INTEGER                              :: i      ! looper
!jhan: deal with this hard-wiring
INTEGER, parameter                   :: iurban=15 ! CABLE
INTEGER, parameter                   :: npft=13 ! CABLE


!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE (6,'(A)') 'Initialising tile canopy '//&
      'water from input mean canopy water'
END IF

!----------------------------------------------------------------------
! Find mean canopy water in input fields and read it in
!----------------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_mean_canopyw,         &
                 fields_in, field_count_in, pos)
mean_canopyw_in => fields_in(pos)
CALL Rcf_Alloc_Field( mean_canopyw_in )
CALL Rcf_Read_Field( mean_canopyw_in, hdr_in, decomp_rcf_input )

!----------------------------------------------------------------------
! Initialise mean_canopyw_out - and allocate space
!----------------------------------------------------------------------
CALL Rcf_Field_Equals( mean_canopyw_out, mean_canopyw_in )
mean_canopyw_out % rows            = canopy_water % rows
mean_canopyw_out % row_len         = canopy_water % row_len
mean_canopyw_out % level_size      = canopy_water % level_size
mean_canopyw_out % glob_rows       = canopy_water % glob_rows
mean_canopyw_out % glob_row_len    = canopy_water % glob_row_len
mean_canopyw_out % glob_level_size = canopy_water % glob_level_size

IF (h_int_active) THEN
  mean_canopyw_in % interp = interp_h_only
ELSE
  mean_canopyw_in % interp = interp_copy
END IF

CALL Rcf_Alloc_Field( mean_canopyw_out )

!----------------------------------------------------------------------
! Interpolate mean_canopyw and copy result to canopy_water levels 1-6
! (levels 1-5 are vegetation, level 6 is urban).  On other levels
! (water, bare soil, ice) set canopy_water to 0.
!----------------------------------------------------------------------
CALL Rcf_Interpolate( mean_canopyw_in, mean_canopyw_out, input_grid, &
                      output_grid, dummy, dummy )

!CABLE: option for ntiles=17 (hard-wired reflecting input dataset)
IF(canopy_water % levels == 1) THEN
  canopy_water % DATA(:,1) = mean_canopyw_out % DATA(:,1)
ELSE IF(canopy_water % levels == 9) THEN
  DO i = 1, 6
    canopy_water % DATA(:,i) = mean_canopyw_out % DATA(:,1)
  END DO
  DO i = 7, 9
    canopy_water % DATA(:,i) = 0.0
  END DO
ELSE 
  DO i = 1, npft
    canopy_water % DATA(:,i) = mean_canopyw_out % DATA(:,1)
  END DO
  DO i = npft+1, canopy_water % levels
    canopy_water % DATA(:,i) = 0.0
  END DO
  i=iurban 
  canopy_water % DATA(:,i) = mean_canopyw_out % DATA(:,1)
ENDIF

!----------------------------------------------------------------------
! Tidy up
!----------------------------------------------------------------------
CALL Rcf_Dealloc_Field( mean_canopyw_out )
CALL Rcf_Dealloc_Field( mean_canopyw_in )

END SUBROUTINE Rcf_Init_Canopy_Water

END MODULE Rcf_Init_Canopy_Water_Mod
