! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

! Description:
!   Provides error handlers for MPI for use by different communicators 
!   in the IOS universe
!
! Method:
!   We can't call mpl_abort as we are already in an MPI fail condition,
!   provide traceback if possible and abort (providing core as the OS 
!   specifies)
 
MODULE ios_mpi_error_handlers

  USE MPL
  IMPLICIT NONE
  CONTAINS

  SUBROUTINE ios_mpi_error_common(comm,code,string)
    IMPLICIT NONE
! These need to use Kind types of the underlying MPI since it will be MPI
! calling this routine not GCOM.
    INTEGER(KIND=mpl_int_kind), INTENT(IN) :: comm
    INTEGER(KIND=mpl_int_kind), INTENT(IN) :: code
    INTEGER                                :: my_code    
    INTEGER                                :: string_length
    INTEGER                                :: ierr
    CHARACTER (LEN=mpl_max_error_string)   :: error_string
    CHARACTER(LEN=*)                          :: string

    my_code = code              ! convert to model kind for MPL call
    CALL mpl_error_string(my_code, error_string, string_length, ierr)
    WRITE(6,'(A)')'An error occured inside the MPI library during an operation'
    WRITE(6,'(A,A)')'on the communicator: ',string
    WRITE(6,'(A,I16,A,I16,A)')'MPI_COMMUNICATOR=',comm,         &
                              ' MPI_ERROR_CODE=',code,          &
                              ' aborting...'
    WRITE(6,'(A)') TRIM(error_string)



    CALL abort()
  END SUBROUTINE ios_mpi_error_common

  SUBROUTINE global_mpi_error_handler(comm,code)
    IMPLICIT NONE
    INTEGER(KIND=mpl_int_kind), INTENT(IN) :: comm
    INTEGER(KIND=mpl_int_kind), INTENT(IN) :: code
    CALL ios_mpi_error_common(comm,code,'global')
  END SUBROUTINE global_mpi_error_handler

  SUBROUTINE model_mpi_error_handler(comm,code)
    IMPLICIT NONE
    INTEGER(KIND=mpl_int_kind), INTENT(IN) :: comm
    INTEGER(KIND=mpl_int_kind), INTENT(IN) :: code
    CALL ios_mpi_error_common(comm,code,'model')
  END SUBROUTINE model_mpi_error_handler

  SUBROUTINE io_mpi_error_handler(comm,code)
    IMPLICIT NONE
    INTEGER(KIND=mpl_int_kind), INTENT(IN) :: comm
    INTEGER(KIND=mpl_int_kind), INTENT(IN) :: code
    CALL ios_mpi_error_common(comm,code,'io')
  END SUBROUTINE io_mpi_error_handler

  SUBROUTINE leader_mpi_error_handler(comm,code)
    IMPLICIT NONE
    INTEGER(KIND=mpl_int_kind), INTENT(IN) :: comm
    INTEGER(KIND=mpl_int_kind), INTENT(IN) :: code
    CALL ios_mpi_error_common(comm,code,'leader')
  END SUBROUTINE leader_mpi_error_handler

  SUBROUTINE async_mpi_error_handler(comm,code)
    IMPLICIT NONE
    INTEGER(KIND=mpl_int_kind), INTENT(IN) :: comm
    INTEGER(KIND=mpl_int_kind), INTENT(IN) :: code
    CALL ios_mpi_error_common(comm,code,'async')
  END SUBROUTINE async_mpi_error_handler

END MODULE ios_mpi_error_handlers
