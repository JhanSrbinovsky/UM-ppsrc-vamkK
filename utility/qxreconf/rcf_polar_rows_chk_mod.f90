! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE Rcf_Polar_Rows_Chk_Mod

!
! Subroutine Rcf_Polar_Rows_Chk
!
! Description:
!   This subroutine takes a global dump and checks the north and south
!   polar rows of every level of every field on the P grid.
!   If it finds that such a row contains non-uniform data it aborts the
!   job, otherwise nothing happens.
!
! Method:
!   Assume the rows are uniform and run through each field level by level,
!   comparing entries against a base value and noting any inconsistencies.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!   Roddy Sharp
!
! Code description:
!   Language: Fortran 90
!   Standard: UMDP3 v8.2

CONTAINS

SUBROUTINE Rcf_Polar_Rows_Chk( fields, field_count, grid, decomp, hdr, &
                               local_lsm, data_source )

USE Rcf_Field_Type_Mod, ONLY:       &
    field_type

USE Rcf_Global_To_Local_Mod, ONLY:  &
    Rcf_Get_Fld_Type

USE Rcf_data_source_Mod, ONLY:      &
    data_source_type,               &
    ancillary_file

USE UM_ParVars, ONLY:          &
    fld_type_p,                     &
    atNorth, atSouth,               &
    gc_proc_row_group,              &
    nproc_x, nproc

USE Ereport_Mod, ONLY:              &
    Ereport

USE Rcf_Read_Field_Mod, ONLY :      &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY :     &
    Rcf_Alloc_Field,                &
    Rcf_Dealloc_Field

USE Rcf_UMhead_Mod, ONLY :          &
    um_header_type

USE Rcf_Grid_Type_Mod, ONLY:        &
    grid_type

USE mask_compression, ONLY: expand_from_mask

USE cppxref_mod, ONLY :             &
    ppx_atm_compressed,             &
    ppx_type_real,                  &
    ppx_type_int,                   &
    ppx_type_log

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER           :: fields(:)
INTEGER, INTENT(IN)                   :: field_count
TYPE (grid_type), TARGET, INTENT(IN)  :: grid
TYPE( um_header_type ), INTENT(IN)    :: hdr
INTEGER, INTENT(IN)                   :: decomp
LOGICAL, TARGET                       :: local_lsm(:)
TYPE( data_source_type ), POINTER     :: data_source(:)

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
INTEGER                       :: i               ! looper (field)
INTEGER                       :: j               ! looper (level)
INTEGER                       :: k               ! looper (column)
INTEGER                       :: pole            ! looper (pole)
INTEGER                       :: fld_type        ! field on P, U or V grid?
INTEGER                       :: info            ! error return from GCOM calls
INTEGER                       :: ErrorStatus     ! error reporting
INTEGER                       :: faults          ! flag for non-uniform rows
INTEGER                       :: row_start(1:2)  ! start index of searchable row
INTEGER                       :: row_end(1:2)    ! end index of searchable row
INTEGER                       :: counter         ! index counter

REAL, ALLOCATABLE             :: row_data(:)  ! temp store for polar row data
REAL                          :: polar_max    ! maximum along polar row
REAL                          :: polar_min    ! minimum along polar row
REAL                          :: MDI          ! generic Missing Data Indicator

CHARACTER (LEN=20), PARAMETER :: RoutineName = 'Rcf_Check_Polar_Rows'
CHARACTER (LEN=80)            :: Cmessage

TYPE (field_type), POINTER    :: field        ! temp structure for each level
TYPE (field_type), TARGET     :: field_decmp  ! for uncompressed fields

! Initialise faults (assume everything is uniform)
faults = 0

DO i = 1, field_count               ! loop over all fields

! Only need to check ancillaries
! Everything else should be an invalid dump or already interpolated
  IF ( data_source( i ) % source /= ancillary_file ) CYCLE

  fld_type = Rcf_Get_Fld_Type( fields(i) % stashmaster % grid_type )

! Only want to check fields lying on the P grid
  IF ( fld_type == fld_type_p ) THEN

! Get field data
    field => fields(i)
    CALL Rcf_Alloc_Field( field )
    CALL Rcf_Read_Field( field, hdr, decomp )

! Compressed data needs unpacking first:
    IF (field % stashmaster % grid_type == ppx_atm_compressed) THEN

      field_decmp % levels          = field % levels
      field_decmp % rows            = grid % loc_p_rows
      field_decmp % row_len         = grid % loc_p_row_length
      field_decmp % level_size      = field_decmp % rows * field_decmp % row_len
      field_decmp % glob_rows       = grid % glob_p_rows
      field_decmp % glob_row_len    = grid % glob_p_row_length
      field_decmp % glob_level_size = field_decmp % glob_rows *                &
                                       field_decmp % glob_row_len
      field_decmp % stashmaster => field % stashmaster

      CALL Rcf_Alloc_Field( field_decmp )

      DO j = 1, field % levels
        CALL expand_from_mask(field_decmp % data(:,j),                         &
            field % data(:, j),                                                &
            local_lsm, field_decmp % level_size,                               &
            field % level_size )
      END DO

      field => field_decmp

    END IF

! Grid data required for searching:
    row_start(1) = 1
    row_end(1)   = field % row_len
    row_start(2) = field % level_size - field % row_len + 1
    row_end(2)   = field % level_size

    ALLOCATE( row_data( field % row_len ) )

! Check data is of a type suitable for testing
    IF (field % stashmaster % data_type == ppx_type_real .OR.                  &
        field % stashmaster % data_type == ppx_type_int  .OR.                  &
        field % stashmaster % data_type == ppx_type_log) THEN

      DO j = 1, field % levels   ! iterate across levels

        DO pole = 1, 2    ! loop over poles

! Set up an array containing the row data, whatever type it is:
          SELECT CASE ( field % stashmaster % data_type )

          CASE (ppx_type_real)                  ! real fields
            row_data(1:field % row_len) =                                      &
               field % data(row_start(pole):row_end(pole), j)
            MDI = RMDI

          CASE (ppx_type_int)                   ! integer fields
            row_data(1:field % row_len) =                                      &
               REAL( field % data_int(row_start(pole):row_end(pole), j) )
            MDI = REAL(IMDI)

! Logicals: set values individually, true = 1, false = 0
          CASE (ppx_type_log)
            counter = 1
            DO k = row_start(pole), row_end(pole)
              IF ( field % data_log(k, j) ) THEN
                row_data(counter) = 1.0
              ELSE
                row_data(counter) = 0.0
              END IF
              counter = counter + 1
            END DO
            MDI = -1.0           ! Dummy value

          END SELECT

! Begin checking polar rows
          IF ( (pole == 1 .AND. atSouth) .OR. (pole == 2 .AND. atNorth) ) THEN

! Initialise max/min (in case all values are MDIs)
            polar_max = 0.0
            polar_min = 0.0

! Search row data for maxima:
            IF (field % row_len > 0) THEN  ! avoid zero-length arrays
              polar_max = MAXVAL(row_data, MASK = row_data /= MDI)
              polar_min = MINVAL(row_data, MASK = row_data /= MDI)
            END IF

! Find maxima across all processors along row:
            CALL gcg_rmax(1, gc_proc_row_group, info, polar_max)
            CALL gcg_rmin(1, gc_proc_row_group, info, polar_min)

! If max > min then row is not uniform; flag up the error
            IF ( polar_max > polar_min )        &
               faults = 1

          END IF          ! end of polar row checking

! If found, ensure all PEs are aware of a fault 
          CALL gc_imax(1, nproc, info, faults)

        END DO     ! end loop over poles

      END DO     ! end of loop across levels

    ELSE

! Not real, int or log data
      ErrorStatus = -10
      Cmessage = 'Unsupported datatype'
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )

    END IF


    IF ( faults > 0 ) THEN
      ErrorStatus = 20

      WRITE(Cmessage, '(A, I3, A, I4)')                                     &
       'Non-uniform polar row in output dump: STASH section ',              &
        field % stashmaster % section, ' item ', field % stashmaster % item

      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF


    DEALLOCATE( row_data )

    IF (field % stashmaster % grid_type == ppx_atm_compressed)                &
      CALL Rcf_Dealloc_Field( field_decmp )

    CALL Rcf_Dealloc_Field( field )

  END IF     ! end test over P grid

END DO     ! end loop over fields

RETURN
END SUBROUTINE Rcf_Polar_Rows_Chk
END MODULE Rcf_Polar_Rows_Chk_Mod
