! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96
!
! Description:
!   This provides minimal parvars like functionality for io servers.


MODULE IOS_Server_Coupler

  INTEGER              :: proc_start
  INTEGER              :: proc_end
  INTEGER              :: procs
  INTEGER              :: proc_row_start
  INTEGER              :: proc_row_end
  INTEGER              :: proc_rows
  INTEGER, POINTER     :: grid_row_start(:)
  INTEGER, POINTER     :: grid_row_end(:)
  INTEGER, POINTER     :: grid_rows(:)
  INTEGER, POINTER     :: grid_point_start(:)
  INTEGER, POINTER     :: grid_point_end(:)
  INTEGER, POINTER     :: grid_points(:)

CONTAINS

  SUBROUTINE IOS_Server_Coupler_Init()
    USE IOS_Model_Geometry
    USE IOS_Common
    USE IOS_Decompose, ONLY : distribute_range
    IMPLICIT NONE

    INTEGER :: gridpoints
    INTEGER :: fld_type
    INTEGER :: proc
    INTEGER :: errCode
    CHARACTER (LEN=132)               :: errorMessage

    IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
      WRITE(6,'(A)')''
      WRITE(6,'(A)')'Coupling relationship: This IOS rank owns'
      WRITE(6,'(A)')''
      WRITE(6,'(A15,A12,A10,A10,A10)')'Item','Field Type','From', &
          'To','Total'
      WRITE(6,'(A15,A12,A10,A10,A10)')'----','----------','----', &
          '--','-----'
    END IF

    CALL distribute_range(1,ios_model_atm_nprocy, &
        proc_row_start,proc_row_end,              &
        model_rank,model_procs)

    proc_rows=proc_row_end-proc_row_start+1

    proc_start=(proc_row_start-1)*ios_model_atm_nprocx
    proc_end  =(proc_row_end    )*ios_model_atm_nprocx-1
    procs=proc_end-proc_start+1

    IF (proc_rows <= 0) THEN
      WRITE(errorMessage,'(A,A,I3,A,I3)')                    &
          'Bad decomposition: This IO Server has no rows, ', &
          'reduce IOS parallelism from ',model_procs,' to ', &
          ios_model_atm_nprocy
      errCode=10
      CALL IOS_Ereport('IOS_Server_coupler',errCode,errorMessage)
    END IF

    IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
      WRITE(6,'(A15,A12,I10,I10,I10)')'processors','   ',     &
          proc_start,proc_end,proc_end-proc_start+1
      WRITE(6,'(A15,A12,I10,I10,I10)')'processor rows','   ', &
          proc_row_start,proc_row_end,                        &
          proc_row_end-proc_row_start+1
    END IF

    ALLOCATE(grid_row_start  (nfld_types))
    ALLOCATE(grid_row_end    (nfld_types))
    ALLOCATE(grid_rows       (nfld_types))
    ALLOCATE(grid_point_start(nfld_types))
    ALLOCATE(grid_point_end  (nfld_types))
    ALLOCATE(grid_points     (nfld_types))
    
    DO fld_type=1,nfld_types
      grid_row_start(fld_type)=offset_map(2,Fld_Type,proc_start)
      grid_row_end(fld_type)  =offset_map(2,Fld_Type,proc_end)+  &
          size_map(2,Fld_Type,proc_end)-1 
      grid_rows(fld_type)=                                       &
          grid_row_end(fld_type)-grid_row_start(fld_type)+1
      
      IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
        WRITE(6,'(A15,I12,I10,I10,I10)')'grid rows',fld_type,    &
            grid_row_start(fld_type),grid_row_end(fld_type),     &
            grid_row_end(fld_type)-grid_row_start(fld_type)+1
      END IF
      gridpoints=0
      IF (proc_start>0) THEN
        DO proc=0,proc_start-1
          gridpoints=gridpoints+          &
              size_map(1,fld_type,proc) * &
              size_map(2,fld_type,proc)
        END DO
      END IF
      grid_point_start(fld_type)=gridpoints+1
      DO proc=proc_start,proc_end
        gridpoints=gridpoints+            &
            size_map(1,fld_type,proc) *   &
            size_map(2,fld_type,proc)
      END DO
      grid_point_end(fld_type)=gridpoints
      grid_points(fld_type)= &
          grid_point_end(fld_type)-grid_point_start(fld_type)+1
      
      IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
        WRITE(6,'(A15,I12,I10,I10,I10)')'grid points',fld_type,  &
            grid_point_start(fld_type),grid_point_end(fld_type), &
            grid_point_end(fld_type)-grid_point_start(fld_type)+1
      END IF
    END DO

    IF (IOS_Verbosity>=IOS_PrStatus_Oper) THEN
      WRITE(6,'(A)')''
    END IF
  END SUBROUTINE IOS_Server_Coupler_Init

END MODULE IOS_Server_Coupler
