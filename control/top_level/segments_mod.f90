! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Defines derived types and subroutines for segmentation control.
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Top Level

MODULE segments_mod

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE atm_fields_bounds_mod, ONLY: array_dims

IMPLICIT NONE
PRIVATE ! Everything is private, unless explicitly exposed later

!------------------------------------------------------------------------------
! Public derived types
!------------------------------------------------------------------------------

! Derived type to contain information about individual segments

!------------------------------------------------------------------------------
! NAME
!   segment_type
!
! COMPONENTS
!   fp          -- Index of the first point of the segment, where the index runs
!                  from 1 to the number of points in the collapsed array.
!   first_x     -- East-West index of the first point in the segment, where the
!                  index runs from 1 to the number of columns.
!   first_y     -- North-South index of the first point in the segment, where
!                  the index runs from 1 to the number of rows.
!   start_index -- The index in list_points on which this segment starts. Where
!                  the list_points array is not used, start_index and fp should
!                  be numerically equal.
!   end_index   -- The index in list_points on which this segment ends.
!   use_points  -- The number of points in this segment. The last segment size
!                  is calculated from the number of grid points to handle for
!                  this MPI task.
!   seg_points  -- The number of points in this segment. The last segment size
!                  is calculated from the total number of horizontal data points
!                  on the MPI domain on this MPI task.
!   start       -- Index of the first point in the segment of collapsed array,
!                  where the index runs from 1 to the number of points in the
!                  collapsed array.
!   finish      -- Index of the last point in the segment of collapsed array.
!
! NOMENCLATURE
!   collapsed array -- Refers to arrays of rank 2 (or rank 2 slices of
!                      higher-ranked arrays) representing horizontal data
!                      points as if they were collapsed into a single rank 1
!                      array with the standard Fortran data ordering pattern
!                      (column-major).
!
! NOTES
!   The components fp and start are numerically equivalent. They have been kept
!   as two separate names because they were named differently and calculated
!   separately in segmentation control code when it resided in individual model
!   code sections themselves.  The different names help the developer to see
!   that the original segmentation functionality and the functionality in this
!   module are equivalent. It may be wise to remove fp in the future as this
!   module becomes more established.
!
!   There may be a case for reviewing the need for both use_points and
!   seg_points.
!------------------------------------------------------------------------------
TYPE segment_type
  INTEGER             :: fp
  INTEGER             :: first_x
  INTEGER             :: first_y
  INTEGER             :: start_index
  INTEGER             :: end_index
  INTEGER             :: use_points
  INTEGER             :: seg_points
  INTEGER             :: start
  INTEGER             :: finish
END TYPE segment_type

! Derived type to contain information about segments as a whole, not individual
! segments.
TYPE meta_segment_type
  INTEGER :: num_parallel !The number of members in the parallel team.
  INTEGER :: npnts !The number of grid points to be handled by
                   !each member of the parallel team.
  INTEGER :: fp    !The first grid point to be handled by each
                   !member of the parallel team.
  INTEGER :: lp    !The last grid point to be handled by each
                   !member of the parallel team.
  INTEGER :: num_segments !The total number of segments to be 
                          !handled by each member in the parallel team.
  INTEGER :: step   !Step from a point in one
                             !segment to the equivalent point in the next.
END TYPE meta_segment_type

! routines, particularly where arrays with levels are passed through as
! subroutine arguments.

TYPE segment_dims_type
  TYPE(array_dims) :: udims
  TYPE(array_dims) :: vdims
  TYPE(array_dims) :: pdims
  TYPE(array_dims) :: qdims
  TYPE(array_dims) :: tdims
  TYPE(array_dims) :: wdims
  TYPE(array_dims) :: udims_s
  TYPE(array_dims) :: vdims_s
  TYPE(array_dims) :: pdims_s
  TYPE(array_dims) :: qdims_s
  TYPE(array_dims) :: tdims_s
  TYPE(array_dims) :: wdims_s
END TYPE segment_dims_type

! Make types publicly accessible
PUBLIC :: segment_type
PUBLIC :: meta_segment_type
PUBLIC :: segment_dims_type

PUBLIC :: array_dims !Defined in atm_fields_bounds_mod. 
                     !Ought to be there, really.

!------------------------------------------------------------------------------
! Public routines
!------------------------------------------------------------------------------

PUBLIC :: segments_mod_seg_meta        ! Initialises meta_segments variables
PUBLIC :: segments_mod_seg_dims        ! Fills in segment dimensions arrays
PUBLIC :: segments_mod_segments        ! Initialises segments variables
PUBLIC :: segments_mod_copy_3d         ! Copy data into array without halos
PUBLIC :: segments_mod_copy_2d         ! Copy data into array without halos

!------------------------------------------------------------------------------
! Private variables
!------------------------------------------------------------------------------

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER :: zhook_in     = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out    = 1

CONTAINS 

!-----------------------------------------------------------------------------
! Calculate segment properties
!-----------------------------------------------------------------------------

!------------------------------------------------------------------------------
! SYNOPSIS
!   segments_mod_seg_meta(meta_segments)
!
! DESCRIPTION
!   Sets up meta-information about the segmentation. The information is placed
!   into a variable of type meta_segment_type.
!
! ARGUMENTS
!   meta_segments -- Output high-level information about segmentation, such as
!                    the number of segments.
!   grid_points   -- The number of grid points to be handled by this MPI task.
!   segment_size  -- The target segment size, if greater than 0.
!   num_segments  -- The number of segments into which to divide the data. This
!                    argument has no effect if segment_size has been set to a
!                    value greater than zero.
!------------------------------------------------------------------------------

SUBROUTINE segments_mod_seg_meta(meta_segments, ipar, num_parallel,  &
                                grid_points, segment_size, num_segments)
  IMPLICIT NONE

  !Arguments
  INTEGER,                  INTENT(IN) :: ipar
  INTEGER,                  INTENT(IN) :: num_parallel
  INTEGER,                  INTENT(IN) :: grid_points
  INTEGER,                  INTENT(IN) :: segment_size
  INTEGER,                  INTENT(IN) :: num_segments
  TYPE(meta_segment_type), INTENT(OUT) :: meta_segments

  !Internal variables
  INTEGER :: num_points
  INTEGER :: remainder_points
  INTEGER :: i
  INTEGER :: running_fp
  INTEGER :: running_lp

  !DrHook
  REAL(KIND=jprb) :: zhook_handle

  IF(lhook) CALL dr_hook('SEGMENTS_MOD:SEGMENTS_MOD_SEG_META', &
                         zhook_in, zhook_handle)

  !Set the number of parallel team members in the derived type.
  meta_segments%num_parallel = num_parallel

  !Divide the total number of grid points between the number of members of the
  !parallel team.
  num_points = grid_points/num_parallel
  meta_segments%npnts = num_points

  !Find the remainder points
  remainder_points = MOD(grid_points, num_parallel)

  !Find the start point and the end point for this parallel team member. The
  !Set up the first one first.
  running_fp = 0
  running_lp = 0

  !Loop over the remaining team members, accumulating the first and last points
  !as we go. If there is only one team member, this loop will not execute at
  !all.
  DO i = 1,ipar
    running_fp = running_lp + 1
    running_lp = running_lp + num_points
    IF (i <= remainder_points) THEN
      running_lp = running_lp + 1
    END IF
  END DO

  !Set up first points and last points in meta_segments type.
  meta_segments%fp = running_fp
  meta_segments%lp = running_lp

  !If this parallel team member is numbered up-to and including the number of
  !remainder points, add another point to this member. Just an exercise in
  !spreading out the remainder points over different team members.
  meta_segments%npnts = num_points
  IF (ipar <= remainder_points) THEN
    meta_segments%npnts = num_points + 1
  END IF

  !Different behaviour depending on the segment size
  IF (segment_size > 0) THEN
    meta_segments%step = segment_size
    meta_segments%num_segments =                   &
      CEILING( REAL(meta_segments%npnts)           &
      / REAL(meta_segments%step ))
  ELSE
    meta_segments%num_segments =                   &
      MIN(meta_segments%npnts, num_segments)
    meta_segments%step = meta_segments%npnts       &
      / meta_segments%num_segments
  END IF

  IF(lhook) CALL dr_hook('SEGMENTS_MOD:SEGMENTS_MOD_SEG_META', &
                         zhook_out, zhook_handle)

  RETURN
END SUBROUTINE segments_mod_seg_meta

!------------------------------------------------------------------------------
! SYNOPSIS
!   segments_mod_seg_dims(segment_dims,row_length,rows,model_levels,wet_levels)
!
! DESCRIPTION
!   Fills the segment_dims variable with array dimension information needed by
!   some segmented routines to carry the total number of points given to a
!   particular MPI task, so that changing the level number results in the right
!   memory address being accessed.
!
! ARGUMENTS
!   segment_dims -- The segment dimension variable to be filled with
!                   appropriate values.
!   row_length   -- The row length for this MPI task.
!   rows         -- The number of rows on this MPI task.
!   model_levels -- The number of model levels.
!   wet_levels   -- The number of wet levels.
!------------------------------------------------------------------------------

SUBROUTINE segments_mod_seg_dims(segment_dims,row_length,rows,           &
                                            model_levels,wet_levels)
  IMPLICIT NONE

  !Arguments
  INTEGER,                 INTENT(IN)  :: row_length
  INTEGER,                 INTENT(IN)  :: rows
  INTEGER,                 INTENT(IN)  :: model_levels
  INTEGER,                 INTENT(IN)  :: wet_levels
  TYPE(segment_dims_type), INTENT(OUT) :: segment_dims

  !Internal variables
  INTEGER :: np_field  !The total number of horizontal points in this MPI task.

  !DrHook
  REAL(KIND=jprb) :: zhook_handle

  IF(lhook) CALL dr_hook('SEGMENTS_MOD:SEGMENTS_MOD_SEG_DIMS', &
                         zhook_in, zhook_handle)

  !Initialisations
  np_field = row_length*rows

  !Set up one set of dimensions. Others will be copies of this.
  segment_dims % pdims % i_start = 1
  segment_dims % pdims % i_end   = np_field
  segment_dims % pdims % j_start = 1
  segment_dims % pdims % j_end   = 1
  segment_dims % pdims % k_start = 1
  segment_dims % pdims % k_end   = model_levels

  !No halos
  segment_dims % qdims = segment_dims % pdims
  segment_dims % wdims = segment_dims % pdims
  segment_dims % tdims = segment_dims % pdims
  segment_dims % udims = segment_dims % pdims
  segment_dims % vdims = segment_dims % pdims

  !Small halos. Same as no halos for now.
  segment_dims % qdims_s = segment_dims % pdims
  segment_dims % wdims_s = segment_dims % pdims
  segment_dims % tdims_s = segment_dims % pdims
  segment_dims % udims_s = segment_dims % pdims
  segment_dims % vdims_s = segment_dims % pdims

  !Handle any deviations from the values set above
  segment_dims % qdims   % k_end = wet_levels
  segment_dims % qdims_s % k_end = wet_levels

  IF(lhook) CALL dr_hook('SEGMENTS_MOD:SEGMENTS_MOD_SEG_DIMS', &
                         zhook_out, zhook_handle)

  RETURN
END SUBROUTINE segments_mod_seg_dims

!------------------------------------------------------------------------------
! SYNOPSIS
!   segments_mod_segments(segments,meta_segments,row_length,rows,list_points)
!
! DESCRIPTION
!   Calculates the number of points and the start and end points for each
!   segment, and stores them in a derived type array. Uses some pre-calculated
!   meta-information. The information about each segment is stored in the
!   segments array.
!
! ARGUMENTS
!   segments      -- The output derived-type array containing information about
!                    each individual segment.
!   meta_segments -- Derived-type variable containing high-level information
!                    about segmentation, such as the number of segments.
!   row_length    -- The number row length in the whole domain for this
!                    MPI task.
!   rows          -- The number of rows in the whole domain for this MPI task.
!   list_points (optional) -- Array of points on which to work.
!
! NOTES
!   If a list_points array (optional) is supplied to the routine, then the
!   total number of grid points may not be the same as the sum of the segment
!   sizes.
!------------------------------------------------------------------------------

SUBROUTINE segments_mod_segments (segments, meta_segments, grid_points, &
     row_length, rows,  list_points)
  IMPLICIT NONE

  !Arguments
  TYPE(meta_segment_type), INTENT(IN)    :: meta_segments
  INTEGER                , INTENT(IN)    :: grid_points
  INTEGER                , INTENT(IN)    :: rows
  INTEGER                , INTENT(IN)    :: row_length
  TYPE(segment_type) , INTENT(INOUT) :: segments(1:meta_segments%num_segments)
  INTEGER, OPTIONAL                      :: list_points(:)

  !Local variables
  LOGICAL :: use_list_points
  INTEGER :: first_point
  INTEGER :: last_point
  INTEGER :: seg_start
  INTEGER :: seg_finish
  INTEGER :: start_index
  INTEGER :: end_index
  INTEGER :: seg_points
  INTEGER :: use_points
  INTEGER :: i
  INTEGER :: j

  !Load-balancing across threads. 
  INTEGER :: num_threads
  INTEGER :: rem_points

  !DrHook
  REAL(kind=jprb) :: zhook_handle

  IF(lhook) CALL dr_hook('SEGMENTS_MOD:SEGMENTS_MOD_SEGMENTS', &
                         zhook_in, zhook_handle)

  !Determine whether a list_points array was provided to the subroutine,
  !and set a logical flag accordingly.
  use_list_points = .FALSE.
  IF(PRESENT(list_points)) THEN
    use_list_points = .TRUE.
    first_point = list_points(meta_segments%fp)
    last_point  = list_points(meta_segments%lp)
  ELSE
    first_point = meta_segments%fp
    last_point  = meta_segments%lp
  END IF

  seg_start   = first_point
  seg_finish  = last_point

  !Handle segmentation for the calling parallel team member only.
  DO i = 1, meta_segments%num_segments 

    !Note that end_index may be different for the last segment, and is modified
    !later.
    start_index = meta_segments%fp + (i-1) * meta_segments%step
    end_index   = start_index + meta_segments%step -1

    !Find the number of points in this segment, and also the number of points to
    !use in this segment. This IF block may reduce further, but kept as it is
    !for the sake of simplicity.
    IF(use_list_points) THEN
      IF (i == meta_segments%num_segments) THEN
        use_points = meta_segments%npnts -         &
                  meta_segments%step*(meta_segments%num_segments-1)
        end_index  = start_index + use_points - 1
      ELSE
        use_points = meta_segments%step
      END IF

      seg_points = list_points(end_index) - seg_start + 1

    ELSE !not using list_points array
      IF (i == meta_segments%num_segments) THEN
        use_points = meta_segments%npnts -           &
          meta_segments%step * (meta_segments%num_segments-1)
        end_index  = start_index + use_points - 1
      ELSE
        use_points = meta_segments%step
      END IF

      seg_points = use_points

    END IF

    !Set first points
    segments(i)%fp         = first_point
    segments(i)%first_y    = (first_point-1)/row_length + 1
    segments(i)%first_x    = first_point-(segments(i)%first_y-1)*row_length

    !Set number of points in segment
    segments(i)%use_points = use_points
    segments(i)%seg_points = seg_points

    !Set start and end point in segments array
    segments(i)%start  = seg_start
    segments(i)%finish = seg_start + seg_points - 1

    !Set start and end indices
    segments(i)%start_index = start_index
    segments(i)%end_index   = end_index

    !Increments
    seg_start   = seg_start   + seg_points
    first_point = first_point + seg_points

  END DO


  IF(lhook) CALL dr_hook('SEGMENTS_MOD:SEGMENTS_MOD_SEGMENTS', &
                         zhook_out, zhook_handle)

  RETURN
END SUBROUTINE segments_mod_segments

!------------------------------------------------------------------------------
! SYNOPSIS
!   segments_mod_copy_3d(input_array,output_array,input_dims,segment,
!                        row_length, rows,model_levels)
!
! DESCRIPTION
!   Packs a portion of a 3D domain contiguously into a 1D array. Allows for
!   2D segmentation only; all vertical levels are copied.
!
! ARGUMENTS
!   input_array  -- The input 2D array.
!   output_array -- The output 1D array.
!   input_dims   -- Array dimensions of the input data array.
!   segment      -- Information about the particular segment to be copied.
!   row_length   -- The row length of the whole domain on this MPI task.
!   rows         -- The number of rows in the whole domain on this MPI task.
!   model_levels -- The number model levels.
!------------------------------------------------------------------------------

SUBROUTINE segments_mod_copy_3d (input_array,output_array,input_dims,   &
                                segment, row_length,rows,model_levels)
  IMPLICIT NONE

  !Arguments
  TYPE(array_dims),   INTENT(IN) :: input_dims
  TYPE(segment_type), INTENT(IN) :: segment
  INTEGER,            INTENT(IN) :: row_length
  INTEGER,            INTENT(IN) :: rows
  INTEGER,            INTENT(IN) :: model_levels

  REAL, INTENT(IN) :: input_array(input_dims%i_start:input_dims%i_end,  &
                                  input_dims%j_start:input_dims%j_end,  &
                                  input_dims%k_start:input_dims%k_end)

  REAL, INTENT(OUT) :: output_array(row_length*rows, 1, 1:model_levels)

  !Internal variables
  INTEGER :: this_x
  INTEGER :: this_y
  INTEGER :: k
  INTEGER :: s

  INTEGER :: s_first
  INTEGER :: s_last

  !DrHook
  REAL(KIND=jprb) :: zhook_handle

  IF(lhook) CALL dr_hook('SEGMENTS_MOD:SEGMENTS_MOD_COPY_3D', &
                         zhook_in, zhook_handle)

  !Get first and last point information from the segment derived type.
  s_first = segment%start
  s_last  = segment%finish

  !Loop over points in this segment
  DO k = 1, model_levels
    DO s = s_first, s_last

      !Find the row and column in which this point resides
      this_y = ((s-1)/row_length) + 1
      this_x = s-(this_y-1)*row_length

      !Copy the array element into the output array
      output_array(s,1,k) = input_array(this_x, this_y, k)

    END DO
  END DO

  IF(lhook) CALL dr_hook('SEGMENTS_MOD:SEGMENTS_MOD_COPY_3D', &
                         zhook_out, zhook_handle)

  RETURN
END SUBROUTINE segments_mod_copy_3d

!------------------------------------------------------------------------------
! SYNOPSIS
!   segments_mod_copy_2d(input_array, output_array, input_dims, segment,
!                        row_length, rows)
!
! DESCRIPTION
!   Packs a portion of a 2D domain contiguously into a 1D array.
!
! ARGUMENTS
!   input_array  -- The input 2D array.
!   output_array -- The output 1D array.
!   input_dims   -- Array dimensions of the input data array.
!   segment      -- Information about the particular segment to be copied.
!   row_length   -- The row length of the whole domain on this MPI task.
!   rows         -- The number of rows in the whole domain on this MPI task.
!------------------------------------------------------------------------------

SUBROUTINE segments_mod_copy_2d (input_array, output_array, input_dims,   &
                                 segment, row_length, rows)
  IMPLICIT NONE

  !Arguments
  TYPE(array_dims),   INTENT(IN) :: input_dims
  TYPE(segment_type), INTENT(IN) :: segment
  INTEGER,            INTENT(IN) :: row_length
  INTEGER,            INTENT(IN) :: rows

  REAL, INTENT(IN) :: input_array(input_dims%i_start:input_dims%i_end,   &
                                  input_dims%j_start:input_dims%j_end)

  REAL, INTENT(OUT) :: output_array(row_length*rows,1)

  !Internal variables
  INTEGER :: this_x
  INTEGER :: this_y
  INTEGER :: k
  INTEGER :: s

  INTEGER :: s_first
  INTEGER :: s_last

  !DrHook
  REAL(KIND=jprb) :: zhook_handle

  IF(lhook) CALL dr_hook('SEGMENTS_MOD:SEGMENTS_MOD_COPY_2D', &
                         zhook_in, zhook_handle)

  !Get first and last point information from the segment derived type.
  s_first = segment%start
  s_last  = segment%finish

  !Loop over points in this segment
  DO s = s_first, s_last

    !Find the row and column in which this point resides
    this_y = ((s-1)/row_length) + 1
    this_x = s-(this_y-1)*row_length     

    !Copy the array element into the output array
    output_array(s,1) = input_array(this_x,this_y)

  END DO

  IF(lhook) CALL dr_hook('SEGMENTS_MOD:SEGMENTS_MOD_COPY_2D', &
                         zhook_out, zhook_handle)

  RETURN
END SUBROUTINE segments_mod_copy_2d

END MODULE segments_mod
