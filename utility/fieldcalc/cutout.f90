! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Given header grid information returns parameters of fixed portion of grid

SUBROUTINE cutout(input_hdr, output_hdr, factor, PPHdrMod)

USE IO, ONLY:             &
  file_close

USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal
  
USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type,         &
  LenFixHd

! Missing data
USE missing_data_mod

! Maximum number of fields to be written
USE maxfields_mod, ONLY: MaxFldsOut

USE stashmaster_mod, ONLY: stm_get_grid_type
USE ereport_mod, ONLY: ereport
USE grdtypes_mod, ONLY: &
  gt_thetamass,         &
  gt_velocity,          &
  gt_u_c,               &
  gt_v_c,               &
  gt_full


IMPLICIT NONE

!
! Description:
!   Cuts out a fixed-grid region of a variable resolution domain
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.


! Subroutine arguments  
TYPE(UM_Header_type), INTENT(IN)    :: input_hdr   !   input file header
REAL, INTENT(IN)                    :: factor(6)   !   domains of interest
LOGICAL, INTENT(IN)                 :: PPHdrMod


TYPE(UM_Header_type), INTENT(INOUT) :: output_hdr  !   output file header

INTEGER, PARAMETER                  :: lbcode_reg_lat_lon     = 1
INTEGER, PARAMETER                  :: lbcode_reg_lat_lon_rot = 101

! Internal variables for each field

INTEGER :: nx, ny             ! No. of points in x/y dir fixed grid
INTEGER :: xstart, ystart     ! Array posns of start of fixed grid
REAL    :: origin_x, origin_y ! Lat/long of first fixed res point
REAL    :: dx, dy             ! Fixed res grid spacing
INTEGER :: grid_type          ! The grid type for the field
INTEGER :: model_id
! Internal variables - grid (1 is the P points and 2 is the B grids winds"

INTEGER :: nx1, ny1             ! No. of points in x/y dir fixed p- grid
INTEGER :: xstart1, ystart1     ! Array posns of start of fixed p-grid
REAL    :: origin_x1, origin_y1 ! Lat/long of first fixed res p-point
REAL    :: dx1, dy1             ! Fixed res grid spacing

INTEGER :: nx2, ny2             ! No. of points in x/y dir fixed b-grid
INTEGER :: xstart2, ystart2     ! Array posns of start of fixed b-grid
REAL    :: origin_x2, origin_y2 ! Lat/long of first fixed res b-point
REAL    :: dx2, dy2             ! Fixed res grid spacing

! Optional specified versions of some of the above.
INTEGER :: nx1_nl, ny1_nl          ! No. of points in x/y dir fixed p- grid
INTEGER :: xstart1_nl, ystart1_nl  ! Array posns of start of fixed p-grid

! x-y of fixed grid of interest (converted from Factor)
INTEGER :: domainx, domainy 

! Arguments to gt_decode to derive grid information.
INTEGER :: gt_model_type
INTEGER :: gt_content
INTEGER :: gt_coverage
INTEGER :: gt_domain
INTEGER :: gt_cyclic

! Info to record whether grid is supported.
LOGICAL :: l_valid_grid
INTEGER :: num_skipped

INTEGER :: NumLookups         ! Number of items in lookup

INTEGER :: i, ii, j, jj, k    ! DO loop iterators
INTEGER :: xpos, ypos         ! Index positions
INTEGER :: kk, totaldata      ! Counter and max for data points
LOGICAL :: LDecode            ! Unpack on read indicator

REAL :: PackAcc               ! Original packing accuracy of field
INTEGER :: LBPack

! Fields for internal use
REAL, POINTER :: TempArray(:,:)
TYPE(PP_Field_type) :: InputField, OutputField
TYPE(PP_Field_type) :: XpndField, CompField
TYPE(UM_Header_type) :: temp_hdr  !   temporary file header
INTEGER :: numcols
INTEGER :: numrows
INTEGER :: max_numcols
INTEGER :: max_numrows

CHARACTER(LEN=*), PARAMETER :: CPackType = "WGDOS"
CHARACTER(LEN=*), PARAMETER :: PPFormat = "p"

! Error reporting variables

CHARACTER(LEN=*), PARAMETER :: RoutineName = "cutout"
INTEGER :: ErrorStatus                 ! Program status monitor
CHARACTER(LEN=80) :: ErrMessage

! End of header

! Set error status to OK
ErrorStatus = StatusOK

! Ascertain fixed grid portion from namelist
domainx    = NINT(Factor(1))
domainy    = NINT(Factor(2))
xstart1    = 0
xstart1_nl = NINT(Factor(3))
ystart1    = 0
ystart1_nl = NINT(Factor(4))
nx1        = 0
nx1_nl     = NINT(Factor(5))
ny1        = 0
ny1_nl     = NINT(Factor(6))

! Calculate GRID for New Dynamics to begin with.
! Check what type of file we are - variable or fixed resolution.
IF ( input_hdr % FixHd(115) > 1 .OR. input_hdr % FixHd(120) > 1 ) THEN
!DEPENDS ON: ascertain_fixed_grid
  CALL ascertain_fixed_grid(input_hdr % Len1ColDepC,                         &
                            input_hdr % Len1RowDepC,                         &
                            input_hdr % ColDepC(1:input_hdr % Len1ColDepC),  &
                            input_hdr % RowDepC(1:input_hdr % Len1RowDepC),  &
                            origin_x1, origin_y1,                            &
                            dx1, dy1, nx1, ny1, xstart1, ystart1,            &
                            domainx, domainy)
!DEPENDS ON: ascertain_fixed_grid
  CALL ascertain_fixed_grid(input_hdr % Len1ColDepC,                         &
                            input_hdr % Len1RowDepC,                         &
                            input_hdr % ColDepC(input_hdr % Len1ColDepC+1:   &
                                                2*input_hdr % Len1ColDepC),  &
                            input_hdr % RowDepC(input_hdr % Len1RowDepC+1:   &
                                                2*input_hdr % Len1RowDepC),  &
                            origin_x2, origin_y2,                            &
                            dx2, dy2, nx2, ny2, xstart2, ystart2,            &
                            domainx, domainy)

ELSE
  ! Fixed grid.
  ! We are not specifying cutout area.
  IF ( xstart1_nl == 0) THEN
    ErrorStatus = StatusFatal
    CALL EReport(RoutineName, ErrorStatus,                                     &
                 "Need to specify cutout grid for fixed resolution file.")
  END IF
END IF


! We have calculated fixed part of VR now see if we have specified a grid.
IF (xstart1_nl /= 0) THEN
  ! Check if we calculated inner fixed resolution section.
  IF (xstart1 /= 0) THEN
    IF (xstart1_nl >= xstart1 .AND. xstart1_nl + nx1_nl <= xstart1 + nx1 .AND. &
        ystart1_nl >= ystart1 .AND. ystart1_nl + ny1_nl <= ystart1 + ny1) THEN
      ! Reset calculated origin.
      origin_x1 = MODULO(origin_x1 + (xstart1_nl-xstart1)*dx,360.0)
      origin_y1 = origin_y1 + (ystart1_nl-ystart1)*dy
      xstart1 = xstart1_nl
      ystart1 = ystart1_nl
      nx1     = nx1_nl
      ny1     = ny1_nl
    ELSE
      ErrorStatus = StatusFatal
      CALL EReport(RoutineName, ErrorStatus,                                 &
                   "Specified grid is not within fixed part of variable grid.")
     
    END IF
  ELSE
    ! All this is required for fixed resolution information.
    xstart1   = xstart1_nl
    ystart1   = ystart1_nl
    nx1       = nx1_nl
    ny1       = ny1_nl
    dx1       = input_hdr % RealC(1) ! E-W grid spacing in degrees
    dy1       = input_hdr % RealC(2) ! N-S grid spacing in degrees
    origin_y1 = input_hdr % RealC(3) ! Latitude of first position
    origin_x1 = input_hdr % RealC(4) ! Longitude of first position

    ! Add the start position on for theta points
    origin_x1 = MODULO(origin_x1 + (xstart1-1)*dx1,360.0)
    origin_y1 = origin_y1 + (ystart1-1)*dy1
  END IF
  ! All the following is the same for fixed and the variable grid.
  ! Find origin.  Assume New Dynamics and correct if Endgame.
 
  ! Setup wind grid.
  xstart2   = xstart1
  ystart2   = ystart1
  nx2       = nx1
  ! One less row than P grid.
  ny2       = ny1 - 1
  dx2       = dx1
  dy2       = dy1
  ! Calculate offset
  origin_x2 = origin_x1 + 0.5*dx1
  origin_y2 = origin_y1 + 0.5*dy1

  ! Check if we are endgame
  IF (input_hdr % Fixhd(9) == 6) THEN
    ny2       = ny1 + 1
    origin_x2 = origin_x1 - 0.5*dx1
    origin_y2 = origin_y1 - 0.5*dy1
  END IF
END IF

WRITE(6,'(A25)') 'Initial fixed p-grid information'
WRITE(6,'("xstart1: ",I5,1x,"Origin1: ",ES14.6)') xstart1,origin_x1
WRITE(6,'("ystart1: ",I5,1x,"Origin1: ",ES14.6)') ystart1,origin_y1
WRITE(6,'("dx1: ",ES12.4," dy1: ",ES12.4)') dx1,dy1
WRITE(6,'("nx1: ",I6," ny1: ",I6)') nx1,ny1

WRITE(6,'(A25)') 'Initial fixed b-grid information'
WRITE(6,'("xstart2: ",I5,1x,"Origin2: ",ES14.6)') xstart2,origin_x2
WRITE(6,'("ystart2: ",I5,1x,"Origin2: ",ES14.6)') ystart2,origin_y2
WRITE(6,'("dx2: ",ES14.6," dy2: ",ES14.6)') dx2,dy2
WRITE(6,'("nx2: ",I6," ny2: ",I6)') nx2,ny2


! Fixes for New Dynamics.
IF (input_hdr % Fixhd(9) == 3) THEN
  ! Check domain is similar to what a fixed resolution output would be.
  IF (origin_x2 < origin_x1 .OR. origin_y2 < origin_y1) THEN
    ErrorStatus = StatusWarning
    CALL EReport(RoutineName, ErrorStatus,                                   &
                "Fixing origin to be P box corner.")
    origin_x2 = origin_x1 + 0.5*dx1
    origin_y2 = origin_y1 + 0.5*dy1
    xstart2   = xstart1
    ystart2   = ystart1
  END IF
  ! Check for number of rows.
  IF ( ny2 == ny1 .OR. ny2 == ny1 + 1) THEN
    ! For a C grid (New Dynamics) we have 1 less v row.
    ny2 = ny1 - 1
    ! Report warning.
    ErrorStatus = StatusWarning
    CALL EReport(RoutineName, ErrorStatus,                                   &
                 "Fixing grid to have consistent rows.")

  END IF

  IF (origin_x1 > origin_x2 .OR. origin_y1 > origin_y2) THEN
    Errorstatus = StatusFatal
    CALL EReport(RoutineName, ErrorStatus,                                   &
                 "Origin incorrect for staggering")
  END IF

  IF (ny2 /= ny1 - 1) THEN
    Errorstatus = StatusFatal
    CALL EReport(RoutineName, ErrorStatus,                                   &
                 "Wrong number of rows for staggering")
  END IF

! Checks for Endgame.
ELSE IF (input_hdr % Fixhd(9) == 6) THEN
  ! Check domain is similar to what a fixed resolution output would be.
  IF (origin_x2 > origin_x1 .OR. origin_y2 > origin_y1) THEN
    ErrorStatus = StatusWarning
    CALL EReport(RoutineName, ErrorStatus,                                   &
                 "Fixing origin to be P box center.")

    ! Lets make sure Endgame gives same P points as ND.  Winds may be taken
    ! from variable spacing but probably doesnt matter too much.
    origin_x2 = origin_x1 - 0.5*dx1
    origin_y2 = origin_y1 - 0.5*dy1
    xstart2   = xstart1
    ystart2   = ystart1
  END IF

  ! Check number of rows
  IF (ny2 == ny1 .OR. ny2 == ny1 - 1) THEN
    ! If we are endgame always make the winds the first grid point even if its
    ! technically still inside variable grid.
    ny2       = ny1 + 1
    Errorstatus = StatusWarning
    CALL EReport(RoutineName, ErrorStatus,                                   &
                 "Fixing grid to have consistent rows.")

  END IF

  IF (origin_x2 > origin_x1 .OR. origin_y2 > origin_y1) THEN
    Errorstatus = StatusFatal
    CALL EReport(RoutineName, ErrorStatus,                                   &
                 "Origin incorrect for staggering")
  END IF

  IF (ny2 /= ny1 + 1) THEN
    Errorstatus = StatusFatal
    CALL EReport(RoutineName, ErrorStatus,                                   &
                 "Wrong number of rows for staggering")
  END IF
ELSE
  ErrorStatus = StatusFatal
  CALL EReport(RoutineName, ErrorStatus,                                     &
              "Unsupported grid staggering.")
END IF

! Check length of rows
IF ( nx2 /= nx1 ) THEN
  ! The row length has to be equal for the P and u/v grid.
  IF (nx2 == nx1 - 1 .OR. nx2 == nx1 + 1) THEN
    ! We need equal columns.
    nx1 = nx2
  ELSE
  ! Grid is too different to calculate fixed resolution part.
    ErrorStatus = StatusFatal
    CALL EReport(RoutineName, ErrorStatus,                                   &
                 "Fixed grid appears to have inconsistent row lengths.")

  END IF
END IF

! Check we can cutout grid from source grid.
IF (input_hdr % IntC(6) < xstart1 + nx1-1 .AND. input_hdr % Fixhd(4) > 0) THEN
  ErrorStatus = StatusFatal
  CALL EReport(RoutineName, ErrorStatus,                                     &
               "Target grid E-W is outside the source grid")
ELSE IF (input_hdr % IntC(7) < ystart1 + ny1-1) THEN
  ErrorStatus = StatusFatal
  CALL EReport(RoutineName, ErrorStatus,                                     &
               "Target grid N-S is outside the source grid.")
END IF


WRITE(6,'(A25)') 'Final fixed p-grid information'
WRITE(6,'("xstart1: ",I5,1x,"Origin1: ",ES14.6)') xstart1,origin_x1
WRITE(6,'("ystart1: ",I5,1x,"Origin1: ",ES14.6)') ystart1,origin_y1
WRITE(6,'("dx1: ",ES12.6," dy1: ",ES12.6)') dx1,dy1
WRITE(6,'("nx1: ",I6," ny1: ",I6)') nx1,ny1

WRITE(6,'(A25)') 'Final fixed b-grid information'
WRITE(6,'("xstart2: ",I5,1x,"Origin2: ",ES14.6)') xstart2,origin_x2
WRITE(6,'("ystart2: ",I5,1x,"Origin2: ",ES14.6)') ystart2,origin_y2
WRITE(6,'("dx2: ",ES14.6," dy2: ",ES14.6)') dx2,dy2
WRITE(6,'("nx2: ",I6," ny2: ",I6)') nx2,ny2

! Setup output headers

! If input header is global we need to make sure some headers are set to LAM
IF (input_hdr % Fixhd(4) == 0) THEN
  output_hdr % Fixhd(4)    = 3
END IF

! Most items same as input header except for:
output_hdr % IntC(6)    = nx1       ! Number of points E-W
output_hdr % IntC(7)    = ny1       ! Number of points N-S
output_hdr % RealC(1)   = dx1       ! E-W grid spacing in degrees
output_hdr % RealC(2)   = dy1       ! N-S grid spacing in degrees
output_hdr % RealC(3)   = origin_y1 ! Latitude of first position
output_hdr % RealC(4)   = origin_x1 ! Longitude of first position
output_hdr % FixHd(115) = IMDI     ! No row-dependent constants
output_hdr % FixHd(116) = IMDI     ! 
output_hdr % FixHd(117) = IMDI     ! 
output_hdr % FixHd(120) = IMDI     ! No column-dependent constants
output_hdr % FixHd(121) = IMDI     ! 
output_hdr % FixHd(122) = IMDI     ! 
output_hdr % RowDepC    = 0.0
output_hdr % ColDepC    = 0.0      
output_hdr % Len1RowDepC = 0
output_hdr % Len2RowDepC = 0
output_hdr % Len1ColDepC = 0 
output_hdr % Len2ColDepC = 0


! Rewrite the header by reopening the file again with the new values

 CALL File_Close ( output_hdr % UnitNum,                &
                  output_hdr % FileNameEnv,            &
                  LEN_TRIM(output_hdr % FileNameEnv),  &
                  0, 0, ErrorStatus )

IF ( ErrorStatus /= StatusOK ) THEN
  ErrMessage = 'Problem closing output file'

  CALL EReport( RoutineName, ErrorStatus, ErrMessage )
  ErrorStatus = StatusWarning
END IF
 
 temp_hdr = output_hdr

!DEPENDS ON: new_umhdr
 CALL New_UMHdr( output_hdr,        &  ! in
                      MaxFldsOut,      &  ! in
                      temp_hdr,     &  ! inout
                      ErrorStatus )    ! inout
output_hdr = temp_hdr

IF ( ErrorStatus /= StatusOK ) THEN
  ErrMessage = 'Problem writing new fixed-grid header'

  CALL EReport( RoutineName, ErrorStatus, ErrMessage )
  ErrorStatus = StatusWarning
END IF



! Copy each field
WRITE(6,'(A22)') 'CUTOUT: Copying fields'

NumLookups = input_hdr % NumFlds
! Allocate temporary field this is rather wasteful but it is what is performed
! elsewhere.  We should be more clever with the reading of fields from files.
max_numcols = MAXVAL( input_hdr % Lookup(1:NumLookups) % NumCols)
max_numrows = MAXVAL( input_hdr % Lookup(1:NumLookups) % NumRows)

ALLOCATE(TempArray(max_numcols * max_numrows, 1 ))

num_skipped = 0
lookups:  DO i = 1, NumLookups

    ! We require model id - fieldcalc sets this to 10 which but
    ! is really atmosphere (which is 1).
    model_id = input_hdr % lookup(i) % lbuser7
    IF (model_id == 10) THEN
      model_id = 1
    END IF
    ! Find the grid type
    grid_type = stm_get_grid_type(model_id,                            &
                                  input_hdr % lookup(i) % STCode/1000, &
                                  MOD(input_hdr % lookup(i) % STCode,1000))

    ! Default to P grid.
    l_valid_grid = .TRUE.
    ny = ny1
    nx = nx1
    origin_y  = origin_y1
    dy = dy1
    origin_x  = origin_x1
    dx = dx1
    xstart = xstart1
    ystart = ystart1 
 
! Decode the grid_type to give more useful information.  
! DEPENDS ON: gt_decode
    CALL gt_decode(grid_type, gt_model_type, gt_content, gt_coverage,  &
                   gt_domain, gt_cyclic)
! If we are B grid or U on C grid then set col information.
    IF (gt_content /= gt_thetamass) THEN
      l_valid_grid = .FALSE.
      IF (gt_content == gt_velocity .OR. gt_content == gt_u_c) THEN
        nx = nx2
        origin_x = origin_x2
        dx = dx2
        xstart = xstart2
        l_valid_grid = .TRUE.
      END IF  
! If we are B grid or V on C grid then set row information
      IF (gt_content == gt_velocity .OR. gt_content == gt_v_c) THEN
        ny = ny2
        origin_y = origin_y2
        dy = dy2
        ystart = ystart2
        l_valid_grid = .TRUE.
      END IF
    END IF
! Do we really want to try and support more than just full grids.  We would
! need a different LSM for input and output.  Possible but not required.
    IF (gt_domain /= gt_full) THEN
      l_valid_grid = .FALSE.
    END IF

! Check lookup for a lat/lon grid.
    IF (input_hdr % lookup(i) % lbcode /= lbcode_reg_lat_lon .AND. &
        input_hdr % lookup(i) % lbcode /= lbcode_reg_lat_lon_rot) THEN
      l_valid_grid = .FALSE.
    END IF

! We just print out the fields we are skipping.
    IF (.NOT. l_valid_grid) THEN
      WRITE(6,'(A,I6,A)')                                               &
        "INFO: Field with STASH ", input_hdr % lookup(i) % STCode,      &
        " not on a supported grid for cutting out - skipping."
      num_skipped = num_skipped + 1
      CYCLE
    END IF


    ALLOCATE( OutputField % RData( nx, ny) )

! Set up field headers
    ErrorStatus = StatusOK
    InputField % LookupPos = i
    InputField % Hdr = input_hdr % Lookup( InputField % LookupPos )
    OutputField % LookupPos = i
    OutputField % Hdr = output_hdr % Lookup( OutputField % LookupPos )
    
! Save packing accuracy since readflds will unpack the data and reset it.  Also
! save lbpack value.
 
    PackAcc = InputField % Hdr % BAcc
    LBPack  = InputField % Hdr % Lbpack
    
! Read field 
    LDecode = .TRUE.
! Set the input field data to temporary array.
    InputField % RData => TempArray
! DEPENDS ON: readfld
    CALL ReadFld( input_hdr, LDecode, PPHdrMod, InputField,                  &
                    ErrorStatus )
    IF ( ErrorStatus /= StatusOK ) THEN
      WRITE( ErrMessage, '(A36,I4)' )                                        &
             "Could not read input field position ", i
      Errorstatus = StatusWarning

      CALL EReport( RoutineName, ErrorStatus, ErrMessage )
      CYCLE lookups
    END IF

! Reshape array from (nx*ny,1) to (nx,ny)
    NULLIFY(InputField % RData)
    numcols = input_hdr % lookup(i) % numcols
    numrows = input_hdr % lookup(i) % numrows
    ALLOCATE( InputField % RData( numcols, numrows ) )
    InputField % RData = RESHAPE(SOURCE = TempArray,            &
                                 SHAPE  = (/numcols,numrows/))
! TempArray is always pointing to largest buffer required.
    
! Get correct part of field
    DO jj = 1, ny
      DO ii = 1, nx
        xpos = MOD(xstart+ii-2,numcols)+1
        ypos = MOD(ystart+jj-2,numrows)+1
        OutputField % RData(ii,jj) = &
          InputField % RData (xpos, ypos)
      END DO
    END DO

! Setup output field header correctly
    OutputField % Hdr              = InputField % Hdr
    OutputField % Hdr % NumRows    = ny
    OutputField % Hdr % NumCols    = nx
    OutputField % Hdr % ZerothLat  = origin_y - dy    ! Zeroth, so deduct dy
    OutputField % Hdr % LatInt     = dy
    OutputField % Hdr % ZerothLon  = origin_x - dx
    OutputField % Hdr % LonInt     = dx
    OutputField % Hdr % LBLRec     = nx*ny
    OutputField % Hdr % LBHem      = MOD(Output_hdr % Fixhd(4),100)

! Pack field if the original was packed
    IF (PackAcc > -99.0 .AND. MOD(LBPack,10) == 1) THEN
      
      totaldata = nx * ny
      ALLOCATE( XpndField % RData( nx*ny, 1) )
      XpndField % Hdr = OutputField % Hdr
      kk = 0
      DO k = 1, OutputField % Hdr % NumRows
        DO j = 1, OutputField % Hdr % NumCols
          kk = kk+1
          XpndField % RData(kk,1) = OutputField % RData (j,k)
        END DO
      END DO
       
! DEPENDS ON: pack_single
      CALL Pack_Single( PackAcc, CPackType, XpndField,                       &
                          CompField, ErrorStatus )

      IF ( ErrorStatus /= StatusOK ) THEN
        WRITE( ErrMessage, '(A36,I5)' )                                      &
             "Pack Single reported problem for ", OutputField % Hdr %STCode
        Errorstatus = StatusFatal

        CALL EReport( RoutineName, ErrorStatus, ErrMessage )
        CYCLE lookups
      END IF

      DEALLOCATE(OutputField % RData)
      NULLIFY(OutputField % RData)
      ALLOCATE(OutputField % RData (SIZE(CompField % RData),1))

      OutputField % Hdr = CompField % Hdr
      OutputField % RData( 1:SIZE(CompField % RData), :) =                   &
                                                CompField % RData
      DEALLOCATE( CompField % RData )
      NULLIFY( CompField % RData )
      DEALLOCATE( XpndField % RData )
      NULLIFY( XpndField % RData )
    END IF ! Packed field?

! Write field to disk

! DEPENDS ON: FldOut
    CALL FldOut( OutputField, output_hdr, ErrorStatus )

! Deallocate arrays

    IF(ASSOCIATED(OutputField % Rdata)) THEN
      DEALLOCATE( OutputField % RData )
      NULLIFY   ( OutputField % Rdata )
    END IF
    IF(ASSOCIATED(InputField % Rdata)) THEN
      DEALLOCATE( InputField % RData )
      NULLIFY   ( InputField % Rdata )
    END IF

END DO lookups

IF (ASSOCIATED(TempArray)) THEN
  DEALLOCATE(TempArray)
  NULLIFY(TempArray)
END IF
   
WRITE(6,'(A,I6,A,I6,A)') 'CUTOUT: Attempted ', NumLookups, ' and skipped ', &
                         num_skipped, ' fields.'
WRITE(6,'(A)') 'CUTOUT: Complete.'

RETURN
END SUBROUTINE cutout
