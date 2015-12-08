! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel RCF: Transform from global to local co-ordinates:

MODULE Rcf_Global_To_Local_Mod
  
!  Subroutine Global_To_Local_Subdomain : global subdomain boundaries
!                                         to local ones
!  Subroutine global_to_local_RC: converts global row,column co-ords to
!                                 processor co-ordinates plus local
!                                 co-ordinates within the processor.
!   Function Rcf_Get_Fld_Type : Determines P, U or V type of field
!
! Description:
!   Takes a global definition of a subdomain region (in terms of
!   model gridpoints) and translates it into local numbers.
!   This effectively means local co-ordinates of the region of the
!   subdomain which intersects with this processor's area.
!
! Method:
!   Use the datastart variable in PARVARS to see if the requested
!   subdomain intersects with this processor's area, if it does
!   then use datastart to convert to local co-ordinate and do a bit
!   of logic using MAX and MIN to ensure the local co-ordinates
!   actually lie within the local area  Then make any corrections
!   necessary to account for a subdomain which crosses over the
!   0 longitude line. Finally, if L_include_halos is set to
!   .TRUE. - include any relevant halo regions.
!
! Derived from UM 4.5 code
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration

USE field_types
IMPLICIT NONE

CONTAINS

SUBROUTINE RCF_GLOBAL_TO_LOCAL_SUBDOMAIN(                           &
                                   L_include_halosEW,               &
                                   L_include_halosNS,               &
                                   grid_code,procid,                &
                                   global_north_in,global_east_in,  &
                                   global_south_in,global_west_in,  &
                                   local_north,local_east,          &
                                   local_south,local_west)
USE UM_ParVars

IMPLICIT NONE

! Subroutine arguments:
LOGICAL, INTENT(In)   :: L_include_halosEW  ! If true, include East-West
                                            ! halos in local region
LOGICAL, INTENT(In)   :: L_include_halosNS  ! if true incl. North-South
                                            ! halos in local region
INTEGER, INTENT(In)   :: grid_code          ! STASH grid type
INTEGER, INTENT(In)   :: procid             ! process result wanted for
INTEGER, INTENT(In)   :: global_north_in    ! global northern boundary
INTEGER, INTENT(In)   :: global_east_in     ! global eastern boundary
INTEGER, INTENT(In)   :: global_south_in    ! global southern boundary
INTEGER, INTENT(In)   :: global_west_in     ! global western boundary

INTEGER, INTENT(Out)  :: local_north        ! local northern boundary
INTEGER, INTENT(Out)  :: local_south        ! local sothern boundary
INTEGER, INTENT(Out)  :: local_east         ! local eastern boundary
INTEGER, INTENT(Out)  :: local_west         ! local westernboundary

! Parameters and Common blocks
INTEGER, PARAMETER    :: st_no_data = -3    ! Magic number

! Local variables
! Copies of the input arguments, that can be modified for
! wrap-around calculations
INTEGER   :: global_north, global_east, global_south, global_west
INTEGER   :: fld_type           ! is field on P or U or V grid?
INTEGER   :: row_len_nh         ! row length when halos are removed
INTEGER   :: nrows_nh           ! number of rows when halos are removed
INTEGER   :: first_global_pt_EW ! global point number of first and last
INTEGER   :: last_global_pt_EW  ! local points in local area
INTEGER   :: first_global_pt_NS ! in the East-West and
INTEGER   :: last_global_pt_NS  ! North-South directions

! Logicals indicating if this processor contains part of a
! subdomain
LOGICAL   :: NS_intersect, EW_intersect
LOGICAL   :: wrap      ! set to .TRUE. if the subdomain passes over the
                       ! the 0 degree longitude line
LOGICAL   :: fullfield ! if the field is NOT a subdomain

! ------------------------------------------------------------------

! Copy the global_in variables into local variables

global_north=global_north_in
global_east=global_east_in
global_south=global_south_in
global_west=global_west_in

! Find out if the data is on a mass or velocity grid

fld_type = Rcf_Get_Fld_Type(grid_code)

IF (fld_type .EQ. fld_type_unknown) THEN
  WRITE(6,*) 'GLOBAL_TO_LOCAL_SUBDOMAIN encountered ',     &
             'field with gridtype code ',grid_code
  WRITE(6,*) 'Unable to process this field.'
  local_north=st_no_data
  local_south=st_no_data
  local_east=st_no_data
  local_west=st_no_data
  GOTO 9999
ENDIF

! Set up logical indicating if this is a full field, or just
! a subdomain

fullfield= ((( global_west .EQ. 1 ) .AND.            &
             ( global_south .EQ. 1 )) .AND.          &
           (((fld_type .EQ. fld_type_p ) .AND.       &
             ( global_north .EQ. glsizep(2) ) .AND.   &
             ( global_west  .EQ. glsizep(1) )) .OR.   &
            ((fld_type .EQ. fld_type_u ) .AND.       &
             ( global_north .EQ. glsizeu(2) ) .AND.  &
             ( global_west  .EQ. glsizeu(1))) .OR.   &
            ((fld_type .EQ. fld_type_v ) .AND.       &
             ( global_north .EQ. glsizev(2) ) .AND.  &
             ( global_west  .EQ. glsizev(1) ))))

! If this is a fullfield (ie. not a subdomain) the local addressing
! is easy:

IF (fullfield) THEN

  IF (L_include_halosNS) THEN
    local_north = g_lasize(2,fld_type_p, &
        halo_type_single,procid)
    local_south = 1
  ELSE
    local_north = g_lasize(2,fld_type_p, &
        halo_type_single,procid)-Offy
    local_south = 1+Offy
  ENDIF
  IF (L_include_halosEW) THEN
    local_west=1
    local_east=g_lasize(1,fld_type_p, &
        halo_type_single,procid)
  ELSE
    local_west=1+Offx
    local_east=g_lasize(1,fld_type_p, &
        halo_type_single,procid)-Offx
  ENDIF

ELSE ! a subdomain requires some careful analysis:

  IF (fld_type .EQ. fld_type_p) THEN
    row_len_nh=g_blsizep(1,procid)
    nrows_nh=g_blsizep(2,procid)
  ELSE IF (fld_type .EQ. fld_type_u) THEN
    row_len_nh=g_blsizeu(1,procid)
    nrows_nh=g_blsizeu(2,procid)
  ELSE
    row_len_nh=g_blsizev(1,procid)
    nrows_nh=g_blsizev(2,procid)
  ENDIF

! Set up variables giving the global point numbers of the
! start and end of this processor's subdomain

  first_global_pt_EW=g_datastart(1,procid)
  last_global_pt_EW=first_global_pt_EW+row_len_nh-1

  first_global_pt_NS=g_datastart(2,procid)
  last_global_pt_NS=first_global_pt_NS+nrows_nh-1

! If global_east is greater than the global row length, this
! indicates a wrap around - but not in the format this code
! wants - where it expects a wrap around to be indicated by
! the east column being less than the west column.

  IF (global_east .LT. global_west) THEN
    wrap=.TRUE.
  ELSEIF (global_east .GT. glsizep(1)) THEN
    wrap=.TRUE.
    global_east=global_east-glsizep(1)
  ELSE
    wrap=.FALSE.
  ENDIF

  EW_intersect =                                         &
         (( .NOT. wrap) .AND.                            &
          ((global_east .GE. first_global_pt_EW) .AND.   &
           (global_west .LE. last_global_pt_EW)))        &
         .OR.                                            &
         ((wrap) .AND.                                   &
          ((global_west .LE. last_global_pt_EW) .OR.     &
           (global_east .GE. first_global_pt_EW)))

  NS_intersect =                                         &
         ((global_south .LE. last_global_pt_NS) .AND.    &
          (global_north .GE. first_global_pt_NS))

  IF (NS_intersect) THEN

    IF ((global_south .GE. first_global_pt_NS) .AND.     &
        (global_south .LE. last_global_pt_NS)) THEN
! This processor contains the NS start of the subarea
      local_south=global_south-first_global_pt_NS+Offy+1
    ELSE
! This processor is to the North of the start of the subarea
      local_south=1+Offy
    ENDIF

    IF ((global_north .GE. first_global_pt_NS) .AND.     &
        (global_north .LE. last_global_pt_NS)) THEN
! This processor contains the NS end of the subarea
      local_north=global_north-first_global_pt_NS+Offy+1
    ELSE
! This processor is to the South of the subarea
      local_north=nrows_nh+Offy
    ENDIF

  ELSE

    local_north=st_no_data
    local_south=st_no_data

  ENDIF

  IF (EW_intersect) THEN

    IF ((global_west .GE. first_global_pt_EW) .AND.     &
        (global_west .LE. last_global_pt_EW)) THEN
! This processor contains the EW start of the subarea
      local_west=global_west-first_global_pt_EW+Offx+1
    ELSE
! This processor is to the right of the start of the subarea
      local_west=1+Offx
    ENDIF

    IF ((global_east .GE. first_global_pt_EW) .AND.     &
        (global_east .LE. last_global_pt_EW)) THEN
! This processor contains the EW end of the subarea
      local_east=global_east-first_global_pt_EW+Offx+1
    ELSE
! This processor is to the left of the end of the subarea
      local_east=Offx+row_len_nh
    ENDIF

  ELSE

    local_east=st_no_data
    local_west=st_no_data

  ENDIF

ENDIF ! is this a fullfield?

 9999 CONTINUE

RETURN
END SUBROUTINE Rcf_Global_To_Local_Subdomain

!--------------------------------------------------------------------
! This routine is included for completeness
!--------------------------------------------------------------------

! Subroutine Interface:
SUBROUTINE Rcf_Global_To_Local_RC( grid_code,                      &
                                   global_column_in , global_row,  &
                                   processor_x , processor_y,      &
                                   local_column, local_row)

! Description:
! Takes a global co-ordinate, in model gridpoints, and returns
! the processor co-ordinate of the processor containing that
! point, and the local co-ordinates of the point on that processor.

USE UM_ParVars

IMPLICIT NONE

INTEGER, INTENT(In)  :: grid_code          ! STASH grid type code
INTEGER, INTENT(In)  :: global_column_in   ! global column number
INTEGER, INTENT(In)  :: global_row         ! global row number
INTEGER, INTENT(Out) :: processor_x        ! processor X (EW) co-ord.
                                           !       (0->nproc_x)
INTEGER, INTENT(Out) :: processor_y        ! processor Y (NS) co-ord.
                                           ! (0->nproc_y)
INTEGER, INTENT(Out) :: local_column       ! local column no. on proc.
INTEGER, INTENT(Out) :: local_row          ! local row number on proc.

! Parameters and COMMON blocks
INTEGER, PARAMETER   :: st_no_data = -3


! Local variables

INTEGER   :: global_column ! modified version of global_column_in which
!                          ! takes account of domains wrapping over
!                          ! 0 degree longitude
INTEGER   :: fld_type      ! field stored on P grid or U grid?
INTEGER   :: row_len_nh,nrows_nh ! row_len and n_rows when halos removed
INTEGER   :: proc  ! loop counter for loop over processors

! global column and row numbers delimiting a processors area
INTEGER   :: start_col, end_col, start_row, end_row

! ------------------------------------------------------------------

! Find out if the data is on a mass or velocity grid

fld_type=Rcf_Get_Fld_Type(grid_code)

IF (fld_type .EQ. fld_type_unknown) THEN
  WRITE(6,*) 'GLOBAL_TO_LOCAL_RC encountered ',       &
             'field with gridtype code ',grid_code
  WRITE(6,*) 'Unable to process this field.'
  processor_x=st_no_data
  processor_y=st_no_data
  local_column=st_no_data
  local_row=st_no_data
  GOTO 9999
ENDIF

! If global_column_in is more than the global row length, perform
! a wrap around to ensure it falls within the global bounds

IF (global_column_in .GT. glsizep(1)) THEN
  global_column=MOD(global_column_in+1,glsizep(1))-1
ELSE
  global_column=global_column_in
ENDIF

IF ((global_column .LT. 1) .OR.                        &
    (global_row .LT. 1) .OR.                           &
    (global_row .GT. glsizep(2))) THEN

  WRITE(6,*) 'GLOBAL_TO_LOCAL_RC encountered ',            &
             'impossible global row/column co-ordinates ', &
             'row: ',global_row,' column: ',global_column

  processor_x=st_no_data
  processor_y=st_no_data
  local_column=st_no_data
  local_row=st_no_data

ENDIF

! Make a first guess at the processor co-ordinates

processor_x=MIN(global_column/(glsizep(1)/gridsize(1)), nproc_x-1)
processor_y=MIN(global_row/(glsizep(2)/gridsize(2)), nproc_y-1)

proc=processor_x+processor_y*gridsize(1)

IF (fld_type .EQ. fld_type_p) THEN
  row_len_nh=g_blsizep(1,proc)
  nrows_nh=g_blsizep(2,proc)
ELSE IF (fld_type .EQ. fld_type_u) THEN
  row_len_nh=g_blsizeu(1,proc)
  nrows_nh=g_blsizeu(2,proc)
ELSE
  row_len_nh=g_blsizev(1,proc)
  nrows_nh=g_blsizev(2,proc)
ENDIF

start_col=g_datastart(1,proc)
end_col=start_col+row_len_nh-1
start_row=g_datastart(2,proc)
end_row=start_row+nrows_nh-1

! Now iterate around these processors until we hit the right one

DO WHILE  (((global_column .LT. start_col) .OR.     &
            (global_column .GT. end_col  ))         &
        .OR.                                        &
           ((global_row .LT. start_row) .OR.        &
            (global_row .GT. end_row)))


  IF (global_column .LT. start_col) THEN
    processor_x=processor_x-1
  ELSEIF (global_column .GT. end_col) THEN
    processor_x=processor_x+1
  ENDIF

  IF (global_row .LT. start_row) THEN
    processor_y=processor_y-1
  ELSEIF (global_row .GT. end_row) THEN
    processor_y=processor_y+1
  ENDIF

  proc=processor_x+processor_y*gridsize(1)

IF (fld_type .EQ. fld_type_p) THEN
  row_len_nh=g_blsizep(1,proc)
  nrows_nh=g_blsizep(2,proc)
ELSE IF (fld_type .EQ. fld_type_u) THEN
  row_len_nh=g_blsizeu(1,proc)
  nrows_nh=g_blsizeu(2,proc)
ELSE
  row_len_nh=g_blsizev(1,proc)
  nrows_nh=g_blsizev(2,proc)
ENDIF
  start_col=g_datastart(1,proc)
  end_col=start_col+row_len_nh-1
  start_row=g_datastart(2,proc)
  end_row=start_row+nrows_nh-1

ENDDO

! Now we have the processor co-ordinates, we can calculate the
! local co-ordinates.

local_column=Offx+global_column-start_col+1
local_row=Offy+global_row-start_row+1

 9999 CONTINUE

RETURN
END SUBROUTINE Rcf_Global_To_Local_RC

!-------------------------------------------------------------------
!
! Rcf only version
!
!-------------------------------------------------------------------

! Function Interface
INTEGER FUNCTION RCF_GET_FLD_TYPE (grid_type_code)

USE UM_ParVars
USE cppxref_mod, ONLY :                                            &
    ppx_atm_tall,       ppx_atm_tmerid,       ppx_atm_umerid,      &
    ppx_atm_uall,       ppx_atm_cuall,        ppx_atm_cvall,       &
    ppx_atm_river,      ppx_atm_ozone,        ppx_atm_compressed,  &
    ppx_ocn_tall,       ppx_ocn_uall,         ppx_ocn_cuall,       &
    ppx_ocn_cvall,      ppx_ocn_tzonal,       ppx_ocn_uzonal,      &
    ppx_ocn_tmerid,     ppx_ocn_umerid,       ppx_ocn_tfield,      &
    ppx_ocn_ufield 
IMPLICIT NONE

!
! Description:
! Takes a STASH grid type code, and returns which type of
! grid this is - mass or wind grid.

! Subroutine arguments:

INTEGER, INTENT(In)   :: grid_type_code     ! IN : STASH grid type code

! Parameters

IF (((grid_type_code .GE. ppx_atm_tall)    .AND.    &
     (grid_type_code .LE. ppx_atm_tmerid)) .OR.     &
     (grid_type_code .EQ. ppx_atm_ozone)   .OR.     &
     (grid_type_code .EQ. ppx_atm_compressed) .OR.  &
     (grid_type_code .EQ. ppx_ocn_tall)    .OR.     &
     (grid_type_code .EQ. ppx_ocn_cuall)  .OR.      &
     (grid_type_code .EQ. ppx_ocn_tfield)  .OR.     &
     (grid_type_code .EQ. ppx_ocn_tzonal)  .OR.     &
     (grid_type_code .EQ. ppx_ocn_tmerid)) THEN
  Rcf_Get_Fld_Type=fld_type_p
ELSEIF                                              &
   (((grid_type_code .GE. ppx_atm_uall)    .AND.    &
     (grid_type_code .LE. ppx_atm_umerid)) .OR.     &
     (grid_type_code .EQ. ppx_atm_cuall)   .OR.     &
     (grid_type_code .EQ. ppx_ocn_uall)    .OR.     &
     (grid_type_code .EQ. ppx_ocn_ufield)  .OR.     &
     (grid_type_code .EQ. ppx_ocn_uzonal)  .OR.     &
     (grid_type_code .EQ. ppx_ocn_umerid)) THEN
  Rcf_Get_Fld_Type=fld_type_u
ELSEIF                                              &
   (((grid_type_code .EQ. ppx_atm_cvall)   .OR.     &
     (grid_type_code .EQ. ppx_ocn_cvall))) THEN
  Rcf_Get_Fld_Type=fld_type_v
ELSEIF                                              &
   (grid_type_code .EQ. ppx_atm_river) THEN
  Rcf_Get_Fld_Type=fld_type_r
ELSE
  Rcf_Get_Fld_Type=fld_type_unknown
ENDIF

RETURN

END FUNCTION Rcf_Get_Fld_Type

END MODULE Rcf_Global_To_Local_Mod



