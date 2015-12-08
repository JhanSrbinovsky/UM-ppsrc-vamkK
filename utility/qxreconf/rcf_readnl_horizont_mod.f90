! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Read the HORIZONT namelist

MODULE rcf_readnl_horizont_mod

!  Subroutine Rcf_Readnl_Horizont - Read the HORIZONT namelist
!
! Description:
!   Reads the HORIZONT namelist controlling horizontal interpolation
!   and domains.
!
! Method:
!   Data  read and Output_Grid set accordingly.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

USE vrhoriz_grid_mod, ONLY: lambda_input_p, lambda_input_u,       &
    phi_input_p, phi_input_v, horizgrid
    
USE missing_data_mod, ONLY : imdi, rmdi   

IMPLICIT NONE
 
 
! Data for namelist reads - mostly for LAM
INTEGER, PARAMETER     :: iproj = imdi    ! unset as no required.          
INTEGER, SAVE          :: orog_blend_width = imdi     !} For orography blending
REAL, POINTER, SAVE    :: blend_weights(:)            !} zone
 
CONTAINS

SUBROUTINE rcf_readnl_horizont( nft, nft_horizgrid )
    
USE Rcf_Interp_Weights_Mod, ONLY : &
    h_int_method,              &
    bilinear,                  &
    area_weighted,             &
    nearest_neighbour,         &
    smcp_int_nearest_neighbour,&
    l_limit_rotations

USE rcf_model_mod, ONLY : &
    ewspacea,  nsspacea,  &
    delta_lon, delta_lat, &
    frstlona,  frstlata,  &
    polelata,  polelona

USE PrintStatus_mod, ONLY : &
    PrintStatus,            &
    PrStatus_Oper

USE Ereport_Mod, ONLY : &
    Ereport

USE Rcf_Grid_Type_Mod, ONLY : &
    Output_Grid

USE UM_ParVars, ONLY : &
    mype

USE Rcf_cntlatm_mod, ONLY : &
    model_domain,           &
    l_endgame,              &
    l_regular

USE Rcf_FortranIO_Mod, ONLY : &
    Rcf_Free_Unit,            &
    Max_Filename_Len

USE Rcf_HeadAddress_Mod, ONLY : &
    FH_GridStagger_C,           &
    FH_GridStagger_Endgame

USE domain_params

USE missing_data_mod, ONLY : &
    imdi, rmdi

USE nlsizes_namelist_mod, ONLY : &
    a_len2_rowdepc,              &
    a_len2_coldepc

USE check_iostat_mod

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)              :: nft             ! File Unit
INTEGER, INTENT(IN)              :: nft_horizgrid   ! File Unit
LOGICAL                          :: l_exist

! Local vars/params
CHARACTER (LEN=*), PARAMETER     :: RoutineName = 'Rcf_readnl_horizont'
CHARACTER (LEN=Max_Filename_Len) :: FileName
CHARACTER (LEN=Max_Filename_Len) :: FileName_horizgrid
CHARACTER (LEN=80)               :: Cmessage
INTEGER                          :: ErrorStatus
INTEGER, PARAMETER               :: orog_blend_max = 40
INTEGER                          :: i

! Temp storage for namelist information
REAL                             :: orog_blend_weights(orog_blend_max)
INTEGER                          :: icount


NAMELIST /HORIZONT/ h_int_method,                                      &
                    orog_blend_weights,                                &
                    smcp_int_nearest_neighbour, l_limit_rotations

! Set defaults
orog_blend_weights(:)  = rmdi
lambda_input_p(:)      = rmdi
lambda_input_u(:)      = rmdi
phi_input_p(:)         = rmdi
phi_input_v(:)         = rmdi
 
 
! Read horizont Namelist
READ( Unit = nft, Nml=horizont, IOSTAT=ErrorStatus)
! check that inputs are as expected
! for example number of blending_weights <= orog_blend_max
! in such an event the READ would fail producing an errorcode /= 0
! the actual value dept upon the compiler used
CALL check_iostat(errorstatus, "namelist horizont")

! Write out namelist for diagnostic
IF (PrintStatus >= PrStatus_Oper) THEN
  IF ( mype == 0 ) THEN
    WRITE ( 6, horizont )
  END IF
END IF

! Some checks
IF ( h_int_method /= bilinear      .AND. &
     h_int_method /= area_weighted .AND. &
     h_int_method /= nearest_neighbour ) THEN
  Cmessage = 'Unsupported horizontal interpolation method'
  ErrorStatus = 10
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! find the apparent size of the blending width set in the inputs.
orog_blend_width=0
DO icount = 1, orog_blend_max
  IF (orog_blend_weights(icount) /= rmdi ) THEN
    orog_blend_width=orog_blend_width+1
  END IF  
END DO

! Can't have an orography blending for a global model.
! Set width to zero.
IF ( model_domain == mt_global ) orog_blend_width = 0

! Allocate and set the internal weights for orography blending
IF (orog_blend_width > 0) THEN
  DO icount = 1,orog_blend_width
    IF (orog_blend_weights(icount) == rmdi) THEN
      WRITE (Cmessage, '(A,A)')                                &
        'Not all blending weights are set appropriately',       &
        ' Please check orog_blend_weights in HORIZONT'
      ErrorStatus = 51
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )     
    END IF 
  END DO
  
  ALLOCATE( blend_weights( orog_blend_width ) )

  ! blending zone weights come from namelist
  blend_weights( 1 : orog_blend_width ) = &
                 orog_blend_weights( 1 : orog_blend_width )
END IF


! Fill in gaps about Output Grid
Output_Grid % global = ( model_domain == mt_global ) 

! Set gridstagger information
IF (l_endgame) THEN
  Output_Grid % grid_stagger      = FH_GridStagger_Endgame
ELSE
! Assume other than ENDGAME we output C grid.
  Output_Grid % grid_stagger      = FH_GridStagger_C
END IF


! Set up grid u,v,p row_length and rows
Output_Grid % glob_u_row_length = Output_Grid % glob_p_row_length
Output_Grid % glob_v_row_length = Output_Grid % glob_p_row_length
Output_Grid % glob_u_rows       = Output_Grid % glob_p_rows

IF (l_endgame) THEN ! New grid for endgame
  Output_Grid % glob_v_rows       = Output_Grid % glob_p_rows + 1
ELSE
  Output_Grid % glob_v_rows       = Output_Grid % glob_p_rows - 1
END IF

! set up whether rotated or not and compensate for cyclic LAM  

IF ( output_grid % global )  THEN    
    Output_Grid % rotated           = .FALSE.
ELSE    ! Atmos C grid LAM
  IF (polelata /= 90 .OR. polelona /= 0 ) THEN
    Output_Grid % Rotated         = .TRUE.
  ELSE
    Output_Grid % Rotated         = .FALSE.
  END IF
  IF (model_domain == mt_bi_cyclic_lam) THEN
    Output_Grid % glob_v_rows       = Output_Grid % glob_p_rows
  END IF
END IF

! set up field sizes, u,v,p

Output_Grid % glob_p_field = Output_Grid % glob_p_row_length *    &
                             Output_Grid % glob_p_rows
Output_Grid % glob_u_field = Output_Grid % glob_u_row_length *    &
                             Output_Grid % glob_u_rows
Output_Grid % glob_v_field = Output_Grid % glob_v_row_length *    &
                             Output_Grid % glob_v_rows
Output_Grid % glob_r_field = Output_Grid % glob_r_row_length *    &
                             Output_Grid % glob_r_rows

! Set some information which is not set via the namelists for global.
IF (output_grid % global) THEN
  If (l_endgame) Then
    ! Read from stshcomp but are RMDI for global jobs
    ewspacea = 360./Real( Output_Grid % glob_p_row_length )
    nsspacea = 180./Real( Output_Grid % glob_p_rows )
    delta_lon = ewspacea
    delta_lat = nsspacea
  Else
    ewspacea = 360./Real( Output_Grid % glob_p_row_length )
    nsspacea = 180./Real( Output_Grid % glob_p_rows - 1 )
    delta_lon = ewspacea
    delta_lat = nsspacea
  End If
END IF

ALLOCATE(output_grid % lambda_p(output_grid % glob_p_row_length))
ALLOCATE(output_grid % lambda_u(output_grid % glob_u_row_length))
ALLOCATE(output_grid % phi_p   (output_grid % glob_p_rows))
ALLOCATE(output_grid % phi_v   (output_grid % glob_v_rows))

! ------------------------
! Horizontal Grid Namelist
! ------------------------
! Lets default the row and column dependent constants to have 0 entries
! for fixed resolution.
a_len2_rowdepc = 0
a_len2_coldepc = 0

IF (.NOT. L_regular) THEN
! Lets set the row dependent and column dependent sizes - originally set
! with the SIZES namelist but due to rationalisation we can just calculate it.
! Row dependent constants storing p and u
  a_len2_rowdepc = 2
! Column dependent constants storing p and v
  a_len2_coldepc = 2

! Find the filename containing horizontal grid from Env Vars
  CALL Fort_Get_Env('VAR_GRID', 8, FileName_horizgrid,              &
                     Max_Filename_Len, ErrorStatus)

  IF ( ErrorStatus /= 0 ) THEN
    ErrorStatus = 30
    Cmessage =  &
    'Unable to Obtain horizontal grid Filename from Environment'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF

  FileName_horizgrid = TRIM( FileName_horizgrid )

! Check file exists - Do all this by hand (not using file_open)
! as wish to do this on *all* PEs together
  INQUIRE( file=FileName_horizgrid, exist=l_exist )

  IF ( .NOT. l_exist ) THEN
    WRITE (6,'(2A)') 'Horizontal grid file: ', TRIM(FileName_horizgrid)
    ErrorStatus = 40
    Cmessage = ' Horizontal grid Namelist file does not exist!'
    CALL Ereport( RoutineName, ErrorStatus, Cmessage )
  END IF
 
! Read horizgrid Namelist 
! Open the file containing horizontal grid
  OPEN( UNIT=nft_horizgrid, FILE=FileName_horizgrid, IOstat=Errorstatus ) 

! Write out namelist for diagnostic
  If ( PrintStatus >= PrStatus_Oper ) Then
    If ( mype == 0 ) Then  
      Write (6,'(2A)') 'horizontal grid file: ',TRIM(FileName_horizgrid)
    End If
  End If 
  Read( Unit = nft_horizgrid, Nml=horizgrid )
  ! Make sure phi_input_p and phi_input_v are consistent
  phi_input_v(Output_Grid % glob_v_rows + 1) = rmdi
  phi_input_p(Output_Grid % glob_p_rows + 1) = rmdi

  ! Check grid definition matches namelist sizes (similar check in makebc)
  IF (lambda_input_p(Output_Grid % glob_p_row_length+1) /= RMDI   .OR.   &
      lambda_input_p(Output_Grid % glob_p_row_length)   == RMDI   .OR.   &
      phi_input_p(Output_Grid % glob_p_rows+1)  /= RMDI           .OR.   &
      phi_input_p(Output_Grid % glob_p_rows)    == RMDI)  THEN
      WRITE(cmessage, '(A)') 'Dimensions of &HORIZGRID namelist'&
          // ' do not match those in &HORIZONT'
      errorstatus = 45
      CALL ereport(routinename, errorstatus, cmessage)
  END IF

! Close file but let calling routine free the unit since it "owns" the unit.
  CLOSE( UNIT=nft_horizgrid)

! Set the following to BMDI for var grids model.
  delta_lon = rmdi
  delta_lat = rmdi
  frstlona  = rmdi
  frstlata  = rmdi
ELSE
! Endgame and New Dynamics has slightly different staggering so we need to
! make sure its setup correctly.
  IF (l_endgame) THEN
    DO i = 1, output_grid % glob_p_row_length
      lambda_input_p(i) = frstlona + (i-0.5)*delta_lon
    END DO
    DO i = 1, output_grid % glob_u_row_length
      lambda_input_u(i) = frstlona + (i-1)*delta_lon
    END DO
    DO i = 1, output_grid % glob_p_rows
      phi_input_p(i) = frstlata + (i-0.5)*delta_lat
    END DO
    DO i = 1, output_grid % glob_v_rows
      phi_input_v(i) = frstlata + (i-1)*delta_lat
    END DO
  ELSE
    DO i = 1, output_grid % glob_p_row_length
      lambda_input_p(i) = frstlona + (i-1)*delta_lon
    END DO
    DO i = 1, output_grid % glob_u_row_length
      lambda_input_u(i) = frstlona + (i-0.5)*delta_lon
    END DO
    DO i = 1, output_grid % glob_p_rows
      phi_input_p(i) = frstlata + (i-1)*delta_lat
    END DO
    DO i = 1, output_grid % glob_v_rows
      phi_input_v(i) = frstlata + (i-0.5)*delta_lat
    END DO
  END IF
END IF

! Copy grid information into reconfiguration versions rather than using the
! input variables.
output_grid % lambda_p(:) = lambda_input_p(1:output_grid % glob_p_row_length)
output_grid % lambda_u(:) = lambda_input_u(1:output_grid % glob_u_row_length)
output_grid % phi_p   (:) = phi_input_p   (1:output_grid % glob_p_rows)
output_grid % phi_v   (:) = phi_input_v   (1:output_grid % glob_v_rows)

! Set pole coordinates
output_grid % phi_npole    = polelata
output_grid % lambda_npole = polelona

RETURN
END SUBROUTINE rcf_readnl_horizont
END MODULE rcf_readnl_horizont_mod
