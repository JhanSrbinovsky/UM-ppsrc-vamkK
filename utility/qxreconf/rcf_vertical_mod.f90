! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Wrapper for vertical interpolation

MODULE Rcf_vertical_Mod

IMPLICIT NONE

!  Subroutine Rcf_Vertical - wrapper for vertical interpolation
!
! Description:
! This module contains a wrapper subroutine for vertical
! interpolation. Data is left on the current decomposition
! as there should be no spatial dependencies.
!
! Method:
!   If no interp. required, data is copied. Otherwise, loop
!   through levels calling relevant interpolation routine for each.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CONTAINS

SUBROUTINE Rcf_vertical( field_in, field_out, grid_in, grid_out, &
                         heights_in, heights_out )

USE Ereport_Mod, ONLY : &
    Ereport

USE Rcf_Field_Type_Mod, ONLY : &
    field_type

USE Rcf_Grid_Type_Mod, ONLY : &
    grid_type

USE Rcf_V_Int_Ctl_Mod, ONLY : &
    v_int_order

Use Submodel_Mod, Only :      &
    submodel_ident

USE Rcf_Set_Interp_Flags_Mod, ONLY : &
    interp_all,                      &
    interp_v_only,                   &
    interp_h_only,                   &
    interp_copy,                     &
    interp_no_op

USE PrintStatus_mod, ONLY : &
    LTimer

USE vert_interp_mod, ONLY: vert_interp

USE cppxref_mod, ONLY :             &
    ppx_type_real,                  &
    ppx_type_int,                   &
    ppx_type_log

IMPLICIT NONE

! Arguments
TYPE (field_type), INTENT(INOUT)  :: field_in
TYPE (field_type), INTENT(INOUT)  :: field_out
TYPE (grid_type), INTENT(IN)      :: grid_in
TYPE (grid_type), INTENT(IN)      :: grid_out
REAL                              :: heights_in(field_in % level_size,&
                                       0 : grid_in % model_levels+1)
REAL                              :: heights_out(field_out %level_size,&
                                       0 : grid_out % model_levels+1)

! Local Data
CHARACTER (LEN=*), PARAMETER      :: RoutineName='Rcf_vertical'
CHARACTER (LEN=80)                :: Cmessage
INTEGER                           :: ErrorStatus
INTEGER                           :: i
INTEGER                           :: j
INTEGER                           :: k
INTEGER                           :: start_level_in
INTEGER                           :: start_level_out
INTEGER                           :: end_level_in
INTEGER                           :: end_level_out


! DEPENDS ON: timer
If (LTimer) Call Timer( RoutineName, 3 )

! Find bottom and top level for output field.
IF (field_out % bottom_level >= 0) THEN
  start_level_out = field_out % bottom_level
ELSE
  start_level_out = 1
END IF

IF (field_out % top_level >= 0) THEN
  end_level_out = field_out % top_level
ELSE
  end_level_out = field_out % levels
END IF

! Find bottom and top level for input field.
IF (field_in % bottom_level >= 0) THEN
  start_level_in = field_in % bottom_level
ELSE
  start_level_in = 1
END IF

IF (field_in % top_level >= 0) THEN
  end_level_in = field_in % top_level
ELSE
  end_level_in = field_in % levels
END IF


IF (end_level_in - start_level_in + 1 /= field_in % levels) THEN
  Cmessage = 'Input field information is inconsistent!'
  WRITE(6,'(A)') cmessage
  WRITE(6,'(A,I3,A,I3)') "STASH section ",field_in % stashmaster % section, &
                         " Item no. ",field_in % stashmaster % item
  WRITE(6,'(A,I3,A,I3,A,I3,A)') "Referencing ",start_level_in, " to ",      &
              end_level_in, " does not equate to ", field_in % levels, "levels"
  ErrorStatus = 5
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

IF (end_level_out - start_level_out + 1 /= field_out % levels) THEN
  Cmessage = 'Output field information is inconsistent!'
  WRITE(6,'(A)') cmessage
  WRITE(6,'(A,I3,A,I3)') "STASH section ",field_in % stashmaster % section, &
                         " Item no. ",field_in % stashmaster % item
  WRITE(6,'(A,I3,A,I3,A,I3,A)') "Referencing ",start_level_out, " to ",     &
        end_level_out, " does not equate to ", field_out % levels, "levels"
  ErrorStatus = 6
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF


! sizes should be the same, but will check
IF ( field_in % level_size /= field_out % level_size ) THEN
  Cmessage = 'Field sizes have different horizontal sizes!'
  WRITE(6,'(A)') cmessage
  WRITE(6,'(A,I3,A,I3)') "STASH section ",field_in % stashmaster % section, &
                         " Item no. ",field_in % stashmaster % item
  WRITE(6,'(A,I0,A,I0)') "Referencing ",field_in % level_size,              &
                         " is not equal to ", field_out % level_size
  ErrorStatus = 7
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

! start levels should be either 0 or 1
IF ( ABS(start_level_in - start_level_out) > 1 ) THEN
  Cmessage = 'Fields have incompatible starting levels!'
  WRITE(6,'(A)') cmessage
  WRITE(6,'(A,I3,A,I3)') "STASH section ",field_in % stashmaster % section, &
                         " Item no. ",field_in % stashmaster % item
  WRITE(6,'(A,I0,A,I0)') "Referencing input level ", start_level_in,        &
                         " is not 1 level from output level ", start_level_out
  ErrorStatus = 8
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

!-----------------------------------------------------------------
! Is interpolation activated? If not, copy data across is all we
! will do.
!-----------------------------------------------------------------
SELECT CASE( field_in % interp )
  CASE( interp_copy, interp_h_only )         ! Do a copy for vertical

    ! levels should be a supported transition (e.g. ND->EG)
    IF ( ( start_level_out   == 0 .AND. start_level_in == 1 .AND.       &
           field_in % levels /= field_out % levels - 1           ) .OR. &
         ( start_level_out   == 1 .AND. start_level_in == 0 .AND.       &
           field_in % levels /= field_out % levels + 1           ) .OR. &
         ( start_level_out   == start_level_in              .AND.       &
           field_in % levels /= field_out % levels               ) ) THEN
      Cmessage = 'No interpolation, but data field sizes/levels are '// &
                 'different!'
      ErrorStatus = 11
      CALL Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF


    SELECT CASE( field_in % stashmaster % data_type )
      CASE( ppx_type_real )
        CALL copy_input_field_real(field_in, field_out, start_level_in,   &
                                   start_level_out, end_level_out  )
       
      CASE( ppx_type_int )
        IF ( Associated( field_in %  Data ) ) THEN
          CALL copy_input_field_real(field_in, field_out, start_level_in, &
                                     start_level_out, end_level_out  )
                 
        ELSE
          CALL copy_input_field_int(field_in, field_out, start_level_in,  &
                                    start_level_out, end_level_out  )
         
        END IF

      CASE ( ppx_type_log )
        IF ( Associated( field_in % Data ) ) Then
          CALL copy_input_field_real(field_in, field_out, start_level_in,  &
                                     start_level_out, end_level_out  )
          
        ELSE
          CALL copy_input_field_log(field_in, field_out, start_level_in,   &
                                    start_level_out, end_level_out  )
         
        END IF

      CASE DEFAULT
        Cmessage = 'Unsupported Data-Type'
        ErrorStatus = -20
        Call Ereport( RoutineName, ErrorStatus, Cmessage )

    END SELECT

  CASE( interp_all, interp_v_only )      ! Do vertical interpolation

! Rcf_generate_heights makes sure we are using the correct height
! information for the type of level we are on.

    ! Check input level is valid since it is used as index within
    ! heights_in.
    IF ( start_level_in < 0 ) THEN
      Cmessage = 'Input start level is negative.'
      ErrorStatus=45
      Call Ereport( RoutineName, ErrorStatus, Cmessage )
    END IF

    ! Loop through levels of output grid
    j = 1
    DO i = start_level_out, end_level_out

      Call vert_interp( field_in % Data, field_in % level_size,      &
                        field_in % levels, heights_out(1:,i),        &
                        heights_in(:,start_level_in:), v_int_order,  &
                        field_out % Data(1:,j) )

      j = j + 1

    END DO

    ! If we are using theta level 0 (i.e. STASH bottom level code is 38) then 
    ! lets copy theta level 1 to theta level 0.
    ! For vertical wind this is set to 0 
    IF ( field_out % bottom_level == 0 ) THEN
      field_out % Data(:,1) = field_out % data(:,2) 
    END IF


  CASE( interp_no_op )
  ! Do nothing

END SELECT

! DEPENDS ON: timer
IF (LTimer) CALL Timer( RoutineName, 4 )

RETURN

CONTAINS

! Small helper routines to copy the required information.  It uses a different 
! array inside the field type depending on datatype.
SUBROUTINE copy_input_field_real(field_in, field_out, start_level_in, &
                                 start_level_out, end_level_out  )

IMPLICIT NONE
TYPE (field_type), INTENT(INOUT)  :: field_in
TYPE (field_type), INTENT(INOUT)  :: field_out
INTEGER, INTENT(IN)               :: start_level_in
INTEGER, INTENT(IN)               :: start_level_out
INTEGER, INTENT(IN)               :: end_level_out

INTEGER :: i
INTEGER :: j
INTEGER :: k

j = 1
DO i = start_level_out, end_level_out
  k = start_level_out - start_level_in + j
  k = MAX(k,1)
  k = MIN(k,field_in % levels)
  field_out % Data(:,j) = field_in % Data(:,k)
  j = j + 1
END DO

RETURN
END SUBROUTINE copy_input_field_real

SUBROUTINE copy_input_field_int(field_in, field_out, start_level_in, &
                                 start_level_out, end_level_out  )

IMPLICIT NONE
TYPE (field_type), INTENT(INOUT)  :: field_in
TYPE (field_type), INTENT(INOUT)  :: field_out
INTEGER, INTENT(IN)               :: start_level_in
INTEGER, INTENT(IN)               :: start_level_out
INTEGER, INTENT(IN)               :: end_level_out

INTEGER :: i
INTEGER :: j
INTEGER :: k

j = 1
DO i = start_level_out, end_level_out
  k = start_level_out - start_level_in + j
  k = MAX(k,1)
  k = MIN(k,field_in % levels)
  field_out % Data_int(:,j) = field_in % Data_int(:,k)
  j = j + 1
END DO

RETURN
END SUBROUTINE copy_input_field_int

SUBROUTINE copy_input_field_log(field_in, field_out, start_level_in, &
                                 start_level_out, end_level_out  )

IMPLICIT NONE
TYPE (field_type), INTENT(INOUT)  :: field_in
TYPE (field_type), INTENT(INOUT)  :: field_out
INTEGER, INTENT(IN)               :: start_level_in
INTEGER, INTENT(IN)               :: start_level_out
INTEGER, INTENT(IN)               :: end_level_out

INTEGER :: i
INTEGER :: j
INTEGER :: k

j = 1
DO i = start_level_out, end_level_out
  k = start_level_out - start_level_in + j
  k = MAX(k,1)
  k = MIN(k,field_in % levels)
  field_out % Data_log(:,j) = field_in % Data_log(:,k)
  j = j + 1
END DO

RETURN
END SUBROUTINE copy_input_field_log

END SUBROUTINE Rcf_vertical

END MODULE Rcf_vertical_Mod
