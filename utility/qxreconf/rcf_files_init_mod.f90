! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Initialisation of input and output files

Module Rcf_Files_Init_Mod

!  Subroutine Rcf_Files_Init - initialisation of files
!
! Description:
!   Open files to be opened, setup input header and output header
!   default sizes. Set up Input_Grid datatype and copy input
!   dump if rotation required.
!
! Method:
!   Input and output dumps opened with File_Open. Input header
!   setup and read in. Output header default sizes setup from
!   input sizes.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

Contains

Subroutine Rcf_Files_Init( hdr_in, hdr_out )

Use Rcf_HeadAddress_Mod, Only : &
    IC_1stConstRho,         IC_BLevels,            IC_TracerLevs,    &
    IC_XLen,                IC_YLen,               IC_NumLandPoints, &
    IC_PLevels,             IC_WetLevels,          IC_NoCloudLevels, &
    IC_SoilTLevels,         IC_SoilMoistLevs,      IC_NumOzoneLevs,  &
    LDC_EtaRho,             LDC_EtaTheta,          RC_ModelTop,      &
    FH_HorizGrid_Global,    FH_HorizGrid,          IC_ConvectLevs,   &
    FH_HorizGrid_NH,                                                 &
    FH_HorizGrid_LamWrap,   FH_HorizGrid_LamWrapEQ,FH_HorizGrid_SH,  &
    IC_HeightMethod,        IC_RiverRows,          IC_RiverRowLength,&
    FH_GridStagger,         FH_GridStagger_Endgame,FH_GridStagger_A, &
    FH_HorizGrid_LamNoWrapEq, FH_HorizGrid_LamNoWrap,                &
    RC_PoleLat,             RC_PoleLong,                             &
    RC_FirstLat,            RC_FirstLong,                            &
    RC_LatSpacing,          RC_LongSpacing

Use Rcf_Recon_Mod, Only : &
    LEN_DUMPHIST,      GRIB, &
    l_interp_input_only

USE nlsizes_namelist_mod, ONLY:                              &
    a_len_inthd,        a_len_realhd,      a_len2_levdepc,   &
    a_len2_rowdepc,     a_len2_coldepc,    a_len2_flddepc,   &
    a_len_extcnst

Use Rcf_Generate_Heights_Mod, Only : &
    height_gen_original,             &
    height_gen_smooth,               &
    height_gen_ECMWF_Press,          &
    height_gen_ECMWF_Hybrd

Use Submodel_Mod, Only : &
    Internal_Model_List, &
    N_Internal_Model

Use Rcf_NRecon_Mod, Only : &
    DumpProgLevs,      &
    PrimDataLen

Use Rcf_UMhead_Mod, Only : &
    um_header_type

Use Rcf_Grid_Type_Mod, Only : &
    Input_Grid,        &
    Output_Grid

Use Ereport_Mod, Only : &
    Ereport

Use Rcf_FortranIO_Mod, Only : &
    Rcf_Get_Unit,             &
    Max_Filename_Len

USE UM_ParVars, Only : &
    mype

Use Rcf_cntlatm_mod, Only : &
    model_domain


USE PrintStatus_mod, Only : &
    PrintStatus,            &
    PrStatus_Normal

USE rcf_model_mod, ONLY : &
    polelona, polelata

Use Rcf_interp_weights_mod, Only : &
    l_limit_rotations, &
    l_same_rotation


USE IO
USE model_file, Only : model_file_open, model_file_close
USE domain_params
Use Rcf_Address_Mod, Only :       &
    Rcf_Address
USE rcf_readumhdr_mod, ONLY : &
    rcf_readumhdr

iMPLICIT NONE

! Arguments
Type (Um_Header_type), Intent(InOut)    :: hdr_in
Type (Um_Header_type), Intent(InOut)    :: hdr_out

! Comdecks
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! Local variables
Character (Len=*), Parameter            :: RoutineName='Rcf_Files_Init'
Character (Len=80)                      :: Cmessage
Character (Len=Max_Filename_Len)        :: DumpName
Integer                                 :: ErrorStatus
Integer                                 :: err
Integer                                 :: i

! Initialise errorstatus
errorstatus = 0
!------------------------------------------------------------------
! First open and setup Input Dump
!------------------------------------------------------------------
Call Rcf_Get_Unit( hdr_in % UnitNum )
If ( Grib ) Then
  Call File_Open( hdr_in % UnitNum, 'RECONTMP', 8, 0, 0, err)
Else
  Call File_Open( hdr_in % UnitNum, 'AINITIAL', 8, 0, 0, err)
End If

If ( err /= 0 ) Then
  Cmessage    = 'Failed to Open Start Dump'
  ErrorStatus = 10
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If
If (PrintStatus >= PrStatus_Normal .and. mype == 0) Then
  
  write (6,'(A)') ''
  If ( Grib ) Then
    Call Fort_Get_Env( 'RECONTMP', 8, DumpName, Max_Filename_Len, err)
    write (6,'(2A)') 'Input GRIB data : ',TRIM(DumpName)
  Else
    Call Fort_Get_Env( 'AINITIAL', 8, DumpName, Max_Filename_Len, err)
    write (6,'(2A)') 'Input dump : ',TRIM(DumpName)
  End If

End If
Call Rcf_ReadUMhdr( hdr_in )

!-----------------------------------------------------------------
! Set addressing information for output dump.
!-----------------------------------------------------------------
! Define submodel and section/version configuration
! DEPENDS ON: setmodl
Call SETMODL(ErrorStatus,CMESSAGE)
If (ErrorStatus  /=  0) Then
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Output dump addressing
IF (l_interp_input_only) THEN
  CALL rcf_address(hdr_in % lookup)
ELSE
  CALL rcf_address()
END IF


! Question: Do we need to zero IMDIs in the FLH?

!-------------------------------------------------------------------
! Setup Output Dump header sizes (from Input) and open file
!-------------------------------------------------------------------

Call Rcf_Get_Unit( hdr_out % UnitNum )


Call Model_File_Open ( hdr_out % UnitNum, 'ASTART', 6, 1, 0, err)

If ( err /= 0 ) Then
  Cmessage    = 'Failed to Open Output Dump'
  ErrorStatus = 20
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

If (PrintStatus >= PrStatus_Normal .and. mype == 0) Then

  Call Fort_Get_Env( 'ASTART', 6, DumpName, Max_Filename_Len, err)
 
  write (6,'(A)')  ''
  write (6,'(2A)') 'Output dump : ',TRIM(DumpName)


End If

! Initially set some values as in input dump (ensure sizes are >= 0 for 
! allocating)
hdr_out % LenIntC      = MAX(hdr_in % LenIntC,0)
hdr_out % LenRealC     = MAX(hdr_in % LenRealC,0)
hdr_out % Len2LevDepC  = MAX(hdr_in % Len2LevDepC,0)
hdr_out % Len2RowDepC  = MAX(hdr_in % Len2RowDepC,0)
hdr_out % Len1Lookup   = MAX(hdr_in % Len1Lookup,0)
hdr_out % LenCompFldI1 = MAX(hdr_in % LenCompFldI1,0)
hdr_out % LenCompFldI2 = MAX(hdr_in % LenCompFldI2,0)
hdr_out % LenCompFldI3 = MAX(hdr_in % LenCompFldI3,0)
hdr_out % Len2ColDepC  = MAX(hdr_in % Len2ColDepC,0)
hdr_out % Len1FldsOfC  = MAX(hdr_in % Len1FldsOfC,0)
hdr_out % Len2FldsOfC  = MAX(hdr_in % Len2FldsOfC,0)
hdr_out % LenExtraC    = MAX(hdr_in % LenExtraC,0)
hdr_out % LenHistFile  = MAX(hdr_in % LenHistFile,0)

! Other values may be overwritten from Namelist RECON
If ( a_len_inthd /= IMDI ) THEN
  hdr_out % LenIntC = a_len_inthd
End If
If ( a_len_realhd /= IMDI ) THEN
  hdr_out % LenRealC = a_len_realhd
End If
If ( a_len2_levdepc /= IMDI ) THEN
  hdr_out % Len2LevDepC = a_len2_levdepc
End If
If ( a_len2_rowdepc /= IMDI ) THEN
  hdr_out % Len2RowDepC = a_len2_rowdepc
End If
If ( a_len2_coldepc /= IMDI ) THEN
  hdr_out % Len2ColDepC = a_len2_coldepc
End If
If ( a_len2_flddepc /= IMDI ) THEN
  hdr_out % Len2FldsOfC = a_len2_flddepc
End If
If ( a_len_extcnst /= IMDI ) THEN
  hdr_out % LenExtraC = a_len_extcnst
End If
If ( LEN_DUMPHIST /= IMDI ) THEN
  hdr_out % LenHistFile = LEN_DUMPHIST
End If

! Check that sizes of integer constants are big enough for vn5.0
!  dumps
If (hdr_out % LenIntC < 46) Then
  ErrorStatus = 30
  Cmessage = 'Length of Integer Constants needs to be at least&
             & 46 for vn5.0 and higher dumps'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! Other values can be deduced from Output_Grid
! Row length is always same between P and U - even for LAMs where we need to
! keep similar structure to global.
hdr_out % Len1ColDepC  = Output_Grid % glob_p_row_length
! Rows are different between Endgame and ND due to V rows being greater than P.
hdr_out % Len1RowDepC  = MAX(Output_Grid % glob_p_rows, &
                             Output_Grid % glob_v_rows)

hdr_out % Len1LevDepC  = Output_Grid % model_levels + 1

! Others have to be deduced from STASH generated addressing.
hdr_out % Len2Lookup   = 0
hdr_out % LenData      = 0
Do i = 1, N_Internal_Model
  hdr_out % Len2Lookup = hdr_out % Len2Lookup +                    &
                         DumpProgLevs( Internal_Model_List( i ) )
  hdr_out % LenData    = hdr_out % LenData +                       &
                         PrimDataLen( Internal_Model_List( i ) )
End Do

!------------------------------------------------------------------
! Set the input grid data type with info from header
! Assume C grid throughout - Ocean however is B grid (?)
!------------------------------------------------------------------
If ( hdr_in % FixHd( FH_HorizGrid ) == FH_HorizGrid_Global ) Then
  Input_Grid % global   =  .TRUE.
  Input_Grid % Rotated  =  .FALSE.
Else
  Input_Grid % global   = .FALSE.
  If (  hdr_in % RealC( RC_PoleLat ) == polelata    .AND. &
        hdr_in % RealC( RC_PoleLong) == polelona    .AND. &
        l_limit_rotations ) Then
    If (PrintStatus >= PrStatus_Normal .and. mype == 0) Then
      Write (6,'(A)') 'LAM Poles are the same. Will limit rotations performed.'
    End If
    l_same_rotation = .TRUE.
  End If
  If ( hdr_in % RealC( RC_PoleLat ) /= 90 .OR. &
            hdr_in % RealC( RC_PoleLong) /= 0 ) Then
    Input_Grid % Rotated = .TRUE.
  Else
    Input_Grid % Rotated = .FALSE.
  End If
End If

Input_Grid % glob_p_row_length = hdr_in % IntC( IC_XLen )
Input_Grid % glob_p_rows       = hdr_in % IntC( IC_YLen )
! Set grid stagger information
Input_Grid % grid_stagger      = Hdr_in % fixhd(FH_GridStagger)

IF ( Hdr_in % fixhd(FH_GridStagger) == FH_GridStagger_A ) THEN   
! GRIB data on A grid ( v rows not 1 short)
  Input_Grid % glob_u_row_length = hdr_in % IntC( IC_XLen )
  Input_Grid % glob_u_rows       = hdr_in % IntC( IC_YLen )
  Input_Grid % glob_v_row_length = hdr_in % IntC( IC_XLen )
  Input_Grid % glob_v_rows       = hdr_in % IntC( IC_YLen )
  

! The Fixhd has no means of identifying whether a grid is a cyclic
! or bicyclic lam. A wrapping LAM could be E-W N-S or both and the code
! treats wrapping as seperate to cyclic.
! The RCF does not support cyclic LAMs and those users appear
! thankfully not to use this code.

! Switching logic for the LAM/global for use by ND and EG.  

ELSE     ! C grid (global or LAM) Atmos 
  IF (Hdr_In % FixHd(FH_GridStagger) == FH_GridStagger_Endgame) THEN
    Input_Grid % glob_u_row_length = hdr_in % IntC( IC_XLen )
    Input_Grid % glob_u_rows       = hdr_in % IntC( IC_YLen )
    Input_Grid % glob_v_row_length = hdr_in % IntC( IC_XLen )
    Input_Grid % glob_v_rows       = hdr_in % IntC( IC_YLen ) + 1
  ELSE   ! New Dynamics
    Input_Grid % glob_u_row_length = hdr_in % IntC( IC_XLen )
    Input_Grid % glob_u_rows       = hdr_in % IntC( IC_YLen )
    Input_Grid % glob_v_row_length = hdr_in % IntC( IC_XLen )
    Input_Grid % glob_v_rows       = hdr_in % IntC( IC_YLen ) - 1
  END IF
  IF (Hdr_In % FixHd(FH_HorizGrid) == FH_HorizGrid_LamWrap .OR. &
      Hdr_In % FixHd(FH_HorizGrid) == FH_HorizGrid_LamWrapEq) THEN
    Input_Grid % glob_v_rows       = hdr_in % IntC( IC_YLen )
  END IF
END IF


Input_Grid % glob_land_field   = hdr_in % IntC( IC_NumLandPoints )

Input_Grid % model_levels      = hdr_in % IntC( IC_PLevels )

Input_Grid % wet_levels        = hdr_in % IntC( IC_WetLevels )
Input_Grid % cloud_levels      = hdr_in % IntC( IC_NoCloudLevels )
Input_Grid % st_levels         = hdr_in % IntC( IC_SoilTLevels )
Input_Grid % sm_levels         = hdr_in % IntC( IC_SoilMoistLevs )
Input_Grid % bl_levels         = hdr_in % IntC( IC_BLevels )
Input_Grid % ozone_levels      = hdr_in % IntC( IC_NumOzoneLevs )
Input_Grid % tr_levels         = hdr_in % IntC( IC_TracerLevs )
Input_Grid % conv_levels       = hdr_in % IntC( IC_ConvectLevs )
Input_Grid % z_top_of_model    = hdr_in % RealC( RC_ModelTop )
Input_Grid % first_constant_r_rho_level = &
                                 hdr_in % IntC( IC_1stConstRho)

! Original height generation method is set if Integer Const is
! either 1 or MDI (ie unset)

If ( hdr_in % IntC( IC_HeightMethod ) == IMDI .OR.                 &
     hdr_in % IntC( IC_HeightMethod ) == height_gen_original ) Then

  Input_Grid % height_gen_method = height_gen_original

Else If ( hdr_in % IntC( IC_HeightMethod ) == height_gen_smooth ) Then
  Input_Grid % height_gen_method = height_gen_smooth

Else If (hdr_in%IntC(IC_HeightMethod) == height_gen_ECMWF_Press ) Then
  Input_Grid % height_gen_method = height_gen_ECMWF_Press

Else If (hdr_in%IntC(IC_HeightMethod) == height_gen_ECMWF_Hybrd ) Then
  Input_Grid % height_gen_method = height_gen_ECMWF_Hybrd

Else
  ErrorStatus = 40
  Cmessage = 'Input dump has unknown height generation method'
  Call Ereport( RoutineName, ErrorStatus, Cmessage )
End If

! If available set the river routing rows/row_length from header.
! Otherwise set to 0 to allow decomposition to work correctly.
If ( hdr_in % IntC( IC_RiverRows ) /= IMDI ) Then
  Input_Grid % glob_r_rows       = hdr_in % IntC( IC_RiverRows )
  Input_Grid % glob_r_row_length = hdr_in % IntC( IC_RiverRowLength )
Else
  Input_Grid % glob_r_rows       = 0
  Input_Grid % glob_r_row_length = 0
End If

Allocate( Input_grid % eta_theta_levels( 0:Input_Grid % model_levels))
Allocate( Input_grid % eta_rho_levels( Input_Grid % model_levels))
Allocate( Input_grid % rhcrit( Input_Grid % model_levels ) )
Allocate( Input_grid % soil_depths( Input_Grid % sm_levels ) )

Input_grid % eta_theta_levels( 0 : Input_Grid % model_levels ) =  &
             hdr_in % LevDepC( 1 : Input_Grid % model_levels + 1, 1)

Input_grid % eta_rho_levels( 1 : Input_Grid % model_levels ) =    &
           hdr_in % LevDepC( 1 : Input_Grid % model_levels, 2)

Input_Grid % rhcrit(  1 : Input_Grid % model_levels ) =           &
    hdr_in % LevDepC( 1 : Input_Grid % model_levels, 3 )

Input_Grid % soil_depths( 1 : Input_Grid % sm_levels ) =          &
        hdr_in % LevDepC( 1 : Input_Grid % sm_levels, 4 )

Input_Grid % glob_p_field = Input_Grid % glob_p_row_length *        &
                            Input_Grid % glob_p_rows
Input_Grid % glob_u_field = Input_Grid % glob_u_row_length *        &
                            Input_Grid % glob_u_rows
Input_Grid % glob_v_field = Input_Grid % glob_v_row_length *        &
                            Input_Grid % glob_v_rows

ALLOCATE(input_grid % lambda_p(Input_Grid % glob_p_row_length))
ALLOCATE(input_grid % lambda_u(Input_Grid % glob_u_row_length))
ALLOCATE(input_grid % phi_p   (Input_Grid % glob_p_rows))
ALLOCATE(input_grid % phi_v   (Input_Grid % glob_v_rows))

! If input variable grid we can copy the dimensions across.
IF (hdr_in % Len2RowDepc == 2 .AND. hdr_in % Len2ColDepc == 2) THEN
  input_grid % lambda_p(:) = &
    hdr_in % coldepc(1:input_grid % glob_p_row_length,1)
  input_grid % lambda_u(:) = &
    hdr_in % coldepc(1:input_grid % glob_u_row_length,2)
  input_grid % phi_p   (:) = hdr_in % rowdepc(1:input_grid % glob_p_rows,1)
  input_grid % phi_v   (:) = hdr_in % rowdepc(1:input_grid % glob_v_rows,2)
ELSE
! Assume we need to create the grid if we are not a variable grid.
  IF ( Hdr_in % fixhd(FH_GridStagger) == FH_GridStagger_A ) THEN
    DO i = 1, input_grid % glob_p_row_length
      input_grid % lambda_p(i) = Hdr_In % RealC( RC_FirstLong ) + &
        (i-1)*Hdr_In % RealC( RC_LongSpacing )
    END DO
    DO i = 1, input_grid % glob_u_row_length
      input_grid % lambda_u(i) = Hdr_In % RealC( RC_FirstLong ) + &
        (i-1)*Hdr_In % RealC( RC_LongSpacing )
    END DO
    DO i = 1, input_grid % glob_p_rows
      input_grid % phi_p(i) = Hdr_In % RealC( RC_FirstLat ) + &
        (i-1)*Hdr_In % RealC( RC_LatSpacing )
    END DO
    DO i = 1, input_grid % glob_v_rows
      input_grid % phi_v(i) = Hdr_In % RealC( RC_FirstLat ) + &
        (i-1)*Hdr_In % RealC( RC_LatSpacing )
    END DO
  ELSE IF (Hdr_In % FixHd(FH_GridStagger) == FH_GridStagger_Endgame) THEN
    DO i = 1, input_grid % glob_p_row_length
      input_grid % lambda_p(i) = Hdr_In % RealC( RC_FirstLong ) + &
        (i-0.5)*Hdr_In % RealC( RC_LongSpacing )
    END DO
    DO i = 1, input_grid % glob_u_row_length
      input_grid % lambda_u(i) = Hdr_In % RealC( RC_FirstLong ) + &
        (i-1)*Hdr_In % RealC( RC_LongSpacing )
    END DO
    DO i = 1, input_grid % glob_p_rows
      input_grid % phi_p(i) = Hdr_In % RealC( RC_FirstLat ) + &
        (i-0.5)*Hdr_In % RealC( RC_LatSpacing )
    END DO
    DO i = 1, input_grid % glob_v_rows
      input_grid % phi_v(i) = Hdr_In % RealC( RC_FirstLat ) + &
        (i-1)*Hdr_In % RealC( RC_LatSpacing )
    END DO
  ELSE
    DO i = 1, input_grid % glob_p_row_length
      input_grid % lambda_p(i) = Hdr_In % RealC( RC_FirstLong ) + &
        (i-1)*Hdr_In % RealC( RC_LongSpacing )
    END DO
    DO i = 1, input_grid % glob_u_row_length
      input_grid % lambda_u(i) = Hdr_In % RealC( RC_FirstLong ) + &
        (i-0.5)*Hdr_In % RealC( RC_LongSpacing )
    END DO
    DO i = 1, input_grid % glob_p_rows
      input_grid % phi_p(i) = Hdr_In % RealC( RC_FirstLat ) + &
        (i-1)*Hdr_In % RealC( RC_LatSpacing )
    END DO
    DO i = 1, input_grid % glob_v_rows
      input_grid % phi_v(i) = Hdr_In % RealC( RC_FirstLat ) + &
        (i-0.5)*Hdr_In % RealC( RC_LatSpacing )
    END DO
  END IF
END IF

! Set the pole coordinates.
input_grid % phi_npole    = hdr_in % RealC( RC_PoleLat  )
input_grid % lambda_npole = hdr_in % RealC( RC_PoleLong )

!--------------------------------------------------------------------
! If input grid is rotated, the winds will be unrotated and
! written back to file, so will copy the file here and
! make sure the input is not overwritten
!--------------------------------------------------------------------
If ( Input_Grid % Rotated .AND. &
     .NOT. GRIB           .AND. &
     .NOT. l_same_rotation) Then

  Call File_Close( hdr_in % UnitNum, 'AINITIAL', 8, 0, 0, err)
  
  If ( err /= 0 ) Then
    Cmessage    = 'Failed to Close Start Dump'
    ErrorStatus = 50
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

  If (mype == 0 ) Then
    Call Shell('cp $AINITIAL $RECONTMP; chmod +rw $RECONTMP',43)
  End If

  Call File_Open( hdr_in % UnitNum, 'RECONTMP', 8, 1, 0, err)

  If ( err /= 0 ) Then
    Cmessage    = 'Failed to Open Copied Start Dump'
    ErrorStatus = 50
    Call Ereport( RoutineName, ErrorStatus, Cmessage )
  End If

End If

Return
End Subroutine Rcf_Files_Init
End Module Rcf_Files_Init_Mod
