! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Parallel UM : Transform from global to local co-ordinates:
! GLOBAL_TO_LOCAL_RC: converts global row,column co-ordinates to
!                     processor co-ordinates plus local
!                     co-ordinates within the processor.

! Subroutine Interface:
SUBROUTINE global_to_local_rc(grid_code,halo_type,                &
                              global_column_in , global_row,      &
                              processor_x , processor_y,          &
                              local_column, local_row)


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
! Takes a global co-ordinate, in model gridpoints, and returns
! the processor co-ordinate of the processor containing that
! point, and the local co-ordinates of the point on that processor.


! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Subroutine arguments:

INTEGER, INTENT(IN) ::  grid_code         ! IN : STASH grid type code
INTEGER, INTENT(IN) ::  halo_type         ! IN : which type of halo has the field got
INTEGER, INTENT(IN) ::  global_column_in  ! IN : global column number
INTEGER, INTENT(IN) ::  global_row        ! IN : global row number
INTEGER, INTENT(OUT) ::  processor_x      ! OUT : processor X (EW) co-ordinate
                                          !       (0->nproc_x)
INTEGER, INTENT(OUT) ::  processor_y      ! OUT : processor Y (NS) co-ordinate
                                          !       (0->nproc_y)
INTEGER, INTENT(OUT) ::  local_column     ! OUT : local column number on processor
INTEGER, INTENT(OUT) ::  local_row        ! OUT : local row number on processor

! Parameters and COMMON blocks

! Local variables

INTEGER ::  global_column ! modified version of global_column_in which
                          ! takes account of domains wrapping over
                          ! 0 degree longitude
INTEGER ::  fld_type      ! field stored on P grid or U grid?
INTEGER ::  row_len_nh,nrows_nh   ! row_len and n_rows when halos removed
INTEGER ::  proc          ! loop counter for loop over processors

INTEGER ::  get_fld_type  ! function

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! ------------------------------------------------------------------

IF (lhook) CALL dr_hook('GLOBAL_TO_LOCAL_RC',zhook_in,zhook_handle)

! Find out if the data is on a mass or velocity grid

! DEPENDS ON: get_fld_type
fld_type=get_fld_type(grid_code)

IF (fld_type  ==  fld_type_unknown) THEN
  WRITE(6,*) 'GLOBAL_TO_LOCAL_RC encountered ',                   &
    'field with gridtype code ',grid_code
  WRITE(6,*) 'Unable to process this field.'
  processor_x=st_no_data
  processor_y=st_no_data
  local_column=st_no_data
  local_row=st_no_data
  GO TO 9999
END IF

! If global_column_in is more than the global row length, perform
! a wrap around to ensure it falls within the global bounds

IF (global_column_in  >   glsize(1,fld_type)) THEN
  global_column=MOD(global_column_in+1,glsize(1,fld_type))-1
ELSE
  global_column=global_column_in
END IF

IF ((global_column  <   1) .OR.                                   &
    (global_row  <   1) .OR.                                      &
    (global_row  >   glsize(2,fld_type))) THEN

  WRITE(6,*) 'GLOBAL_TO_LOCAL_RC encountered ',                   &
  'impossible global row/column co-ordinates ',                   &
  'row: ',global_row,' column: ',global_column

  processor_x=st_no_data
  processor_y=st_no_data
  local_column=st_no_data
  local_row=st_no_data

END IF

processor_x=g_pe_index_ew(global_column)
processor_y=g_pe_index_ns(global_row)

proc=processor_x+processor_y*gridsize(1)

local_column=halosize(1,halo_type)+global_column-                 &
             g_datastart(1,proc)+1
local_row=halosize(2,halo_type)+global_row-                       &
          g_datastart(2,proc)+1
 9999 CONTINUE

IF (lhook) CALL dr_hook('GLOBAL_TO_LOCAL_RC',zhook_out,zhook_handle)
RETURN
END SUBROUTINE global_to_local_rc

! Function Interface

