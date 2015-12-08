! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Parallel UM : Transform from global to local co-ordinates:
! GLOBAL_TO_LOCAL_SUBDOMAIN: converts global subdomain boundaries
!                            to local subdomain boundaries
!
! Subroutine Interface:
SUBROUTINE global_to_local_subdomain(                             &
                              l_include_halosew,                  &
                              l_include_halosns,                  &
                              grid_code,halo_type,procid,         &
                              global_start_row_in,                &
                              global_end_col_in,                  &
                              global_end_row_in,                  &
                              global_start_col_in,                &
                              local_start_row,local_end_col,      &
                              local_end_row,local_start_col)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE UM_ParVars
USE sterr_mod, ONLY: st_upper_less_lower, st_not_supported,       &
                     st_no_data,st_nd, st_bad_array_param,        &
                     st_bad_address, st_unknown,                  &
                     st_bad_wraparound, st_illegal_weight,        &
                     unknown_weight, unknown_mask,                &
                     unknown_processing, nonsense
IMPLICIT NONE

! Description:
! Takes a global definition of a subdomain region (in terms of
! model gridpoints) and translates it into local numbers.
! This effectively means local co-ordinates of the region of the
! subdomain which intersects with this processor's area.

! Method:
! Use the datastart variable in PARVARS to see if the requested
! subdomain intersects with this processor's area, if it does
! then use datastart to convert to local co-ordinate and do a bit
! of logic using MAX and MIN to ensure the local co-ordinates
! actually lie within the local area  Then make any corrections
! necessary to account for a subdomain which crosses over the
! 0 longitude line. Finally, if L_include_halos is set to
! .TRUE. - include any relevant halo regions.

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Subroutine arguments:

LOGICAL, INTENT(IN) :: l_include_halosew 
                            ! IN : include East-West halos in local
                            !      region if set to .TRUE.
LOGICAL, INTENT(IN) :: l_include_halosns  
                            ! IN : include North-South halos in local
                            !      region if set to .TRUE.
                            
INTEGER, INTENT(IN) :: grid_code     ! IN : STASH grid type of field
INTEGER, INTENT(IN) :: halo_type     ! IN : which type of halo has the field got
INTEGER, INTENT(IN) :: procid        ! IN : processor to produce result for
INTEGER, INTENT(IN) :: global_start_row_in 
                                     ! IN : first row of global subdomain
INTEGER, INTENT(IN) :: global_end_col_in 
                                     ! IN : last column of global subdomain
INTEGER, INTENT(IN) :: global_end_row_in   
                                     ! IN : last row of global subdomain
INTEGER, INTENT(IN) :: global_start_col_in  
                                     ! IN : first column of global subdomain

INTEGER, INTENT(OUT) :: local_start_row   
                                     ! OUT : first row of local subdomain
INTEGER, INTENT(OUT) :: local_end_col   
                                     ! OUT : last column of local subdomain
INTEGER, INTENT(OUT) :: local_end_row  
                                     ! OUT : last row of local subdomain
INTEGER, INTENT(OUT) :: local_start_col    
                                     ! OUT : first column of local subdomain

! Parameters and Common blocks


! Local variables
INTEGER ::                                                        &
! Copies of the input arguments, that can be modified for
! wrap-around calculations
  global_start_row,global_end_col                                 &
, global_end_row,global_start_col                                 &
, fld_type                                                        &
                     ! is field on P or U grid?
, row_len_nh                                                      &
                     ! row length when halos are removed
, nrows_nh                                                        &
                     ! number of rows when halos are removed
, first_global_pt_ew                                              &
                     ! global point number of first and last
, last_global_pt_ew                                               &
                     ! local points in local area
, first_global_pt_ns                                              &
                     ! in the East-West and
, last_global_pt_ns  ! North-South directions

LOGICAL  ::                                                       &
! Logicals indicating if this processor contains part of a
! subdomain
  ns_intersect,ew_intersect                                       &
, wrap                                                            &
                     ! set to .TRUE. if the subdomain passes over the
!                      the 0 degree longitude line
, fullfield ! if the field is NOT a subdomain

INTEGER :: get_fld_type  ! function

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ------------------------------------------------------------------

! Copy the global_in variables into local variables

IF (lhook) CALL dr_hook('GLOBAL_TO_LOCAL_SUBDOMAIN',zhook_in,zhook_handle)
global_start_row=global_start_row_in
global_end_col=global_end_col_in
global_end_row=global_end_row_in
global_start_col=global_start_col_in

! Find out if the data is on a mass or velocity grid

! DEPENDS ON: get_fld_type
fld_type=get_fld_type(grid_code)

IF (fld_type  ==  fld_type_unknown) THEN
  WRITE(6,*) 'GLOBAL_TO_LOCAL_SUBDOMAIN encountered ',            &
    'field with gridtype code ',grid_code
  WRITE(6,*) 'Unable to process this field.'
  local_start_row=st_no_data
  local_end_col=st_no_data
  local_end_row=st_no_data
  local_start_col=st_no_data
  GO TO 9999
END IF

! Set up logical indicating if this is a full field, or just
! a subdomain

fullfield= ((global_start_col  ==  1) .AND.                       &
            (global_end_col  ==  glsize(1,fld_type)) .AND.        &
            (global_start_row  ==  1) .AND.                       &
            (global_end_row  ==  glsize(2,fld_type)))

! If this is a fullfield (ie. not a subdomain) the local addressing
! is easy:

IF (fullfield) THEN

  IF (l_include_halosns) THEN
    local_start_row=1
    local_end_row=g_lasize(2,fld_type,halo_type,procid)
  ELSE
    local_start_row=1+halosize(2,halo_type)
    local_end_row=g_lasize(2,fld_type,halo_type,procid)-          &
                  halosize(2,halo_type)
  END IF
  IF (l_include_halosew) THEN
    local_start_col=1
    local_end_col=g_lasize(1,fld_type,halo_type,procid)
  ELSE
    local_start_col=1+halosize(1,halo_type)
    local_end_col=g_lasize(1,fld_type,halo_type,procid)-          &
                 halosize(1,halo_type)
  END IF

ELSE ! a subdomain requires some careful analysis:

  row_len_nh=g_blsize(1,fld_type,procid)
  nrows_nh=g_blsize(2,fld_type,procid)

! Set up variables giving the global point numbers of the
! start and end of this processor's subdomain

  first_global_pt_ew=g_datastart(1,procid)
  last_global_pt_ew=first_global_pt_ew+row_len_nh-1

  first_global_pt_ns=g_datastart(2,procid)
  last_global_pt_ns=first_global_pt_ns+nrows_nh-1

! If global_east is greater than the global row length, this
! indicates a wrap around - but not in the format this code
! wants - where it expects a wrap around to be indicated by
! the east column being less than the west column.

  IF (global_end_col  <   global_start_col) THEN
    wrap=.TRUE.
  ELSE IF (global_end_col  >   glsize(1,fld_type)) THEN
    wrap=.TRUE.
    global_end_col=global_end_col-glsize(1,fld_type)
  ELSE
    wrap=.FALSE.
  END IF

  ew_intersect =                                                  &
    (( .NOT. wrap) .AND.                                          &
     ((global_end_col  >=  first_global_pt_ew) .AND.              &
      (global_start_col  <=  last_global_pt_ew)))                 &
    .OR.                                                          &
    ((wrap) .AND.                                                 &
     ((global_start_col  <=  last_global_pt_ew) .OR.              &
      (global_end_col  >=  first_global_pt_ew)))

  ns_intersect =                                                  &
    ((global_start_row  <=  last_global_pt_ns) .AND.              &
     (global_end_row  >=  first_global_pt_ns))

  IF (ns_intersect) THEN

    IF ((global_start_row  >=  first_global_pt_ns) .AND.          &
        (global_start_row  <=  last_global_pt_ns)) THEN
! This processor contains the NS start of the subarea
      local_start_row=global_start_row-first_global_pt_ns+        &
                    halosize(2,halo_type)+1
    ELSE
! This processor is to the North of the start of the subarea
      local_start_row=1+halosize(2,halo_type)
    END IF

    IF ((global_end_row  >=  first_global_pt_ns) .AND.            &
        (global_end_row  <=  last_global_pt_ns)) THEN
! This processor contains the NS end of the subarea
      local_end_row=global_end_row-first_global_pt_ns+            &
                    halosize(2,halo_type)+1
    ELSE
! This processor is to the South of the subarea
      local_end_row=halosize(2,halo_type)+nrows_nh
    END IF

  ELSE

    local_start_row=st_no_data
    local_end_row=st_no_data

  END IF

  IF (ew_intersect) THEN

    IF ((global_start_col  >=  first_global_pt_ew) .AND.          &
        (global_start_col  <=  last_global_pt_ew)) THEN
! This processor contains the EW start of the subarea
      local_start_col=global_start_col-first_global_pt_ew+        &
                   halosize(1,halo_type)+1
    ELSE
! This processor is to the right of the start of the subarea
      local_start_col=1+halosize(1,halo_type)
    END IF

    IF ((global_end_col  >=  first_global_pt_ew) .AND.            &
        (global_end_col  <=  last_global_pt_ew)) THEN
! This processor contains the EW end of the subarea
      local_end_col=global_end_col-first_global_pt_ew+            &
                   halosize(1,halo_type)+1
    ELSE
! This processor is to the left of the end of the subarea
      local_end_col=halosize(1,halo_type)+row_len_nh
    END IF

  ELSE

    local_start_col=st_no_data
    local_end_col=st_no_data

  END IF

END IF ! is this a fullfield?

 9999 CONTINUE

IF (lhook) CALL dr_hook('GLOBAL_TO_LOCAL_SUBDOMAIN',zhook_out,zhook_handle)
RETURN
END SUBROUTINE global_to_local_subdomain

! Subroutine Interface:

! Function Interface

