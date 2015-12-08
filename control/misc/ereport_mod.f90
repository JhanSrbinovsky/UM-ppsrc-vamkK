! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc.

MODULE Ereport_Mod

  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim
  USE UM_ParCore

  IMPLICIT NONE

  ! Default private
  PRIVATE

  INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
  
  ! Create an interface to allow 32 or 64 bit errorstatuses to be used.
  INTERFACE ereport
    MODULE PROCEDURE ereport64, ereport32
  END INTERFACE ereport
  
  ! Only make the interface available publicly together with finalise
  PUBLIC :: ereport, ereport_finalise

  ! Counter for number of warnings
  INTEGER :: numwarnings = 0

CONTAINS

  !  Subroutine Ereport - handles errors/warnings
  !
  ! Description:
  !   The standard reconfiguration error and warning reporting routine.
  !   Prints the relevant error and exits if required.
  !
  ! Method:
  !   ErrorStatus > 0 is a hard error
  !   ErrorStatus < 0 is a warning
  !   ErrorStatus = 0 is hunky-dory

  SUBROUTINE Ereport64 ( RoutineName, ErrorStatus, Cmessage )

    USE UM_Types, ONLY : &
      integer64
    IMPLICIT NONE

! Arguments
    CHARACTER (LEN=*), INTENT( IN )   :: RoutineName
    CHARACTER (LEN=*), INTENT( IN )   :: Cmessage

    INTEGER(KIND=Integer64)           :: ErrorStatus
    ! We can't intent this because the UM passes constants.
    
! Local scalars
! Unset is just a temporary variable which is not checked since if we cant
! flush stdout/stderr no way of reporting problems anyway.  However cannot be a
! parameter since um_fort_flush can change it.
    INTEGER                       :: unset = -99
    INTEGER                       :: unit,i
    INTEGER                       :: out_units(2)                     
    REAL(KIND=jprb)               :: zhook_handle
    CHARACTER (LEN=*), PARAMETER  :: astline = REPEAT("?",80)
    CHARACTER (LEN=*), PARAMETER  :: msg ='Job aborted from ereport.'

    ! Temporary parameters to avoid causing problems with triple ? in the
    ! AIX Fortran preprocessor
    CHARACTER (LEN=3), PARAMETER  :: qmark = REPEAT('?',3)
    CHARACTER (LEN=3), PARAMETER  :: exmark = REPEAT('!',3)


    DATA out_units /0,6/
    IF (lhook) CALL dr_hook('EREPORT64',zhook_in,zhook_handle)

    IF ( ErrorStatus > 0 ) THEN
      DO i=1,size(out_units)
        unit=out_units(i)
        IF (mype == 0 .OR. unit /= 0) THEN ! don't get mangled-up 
                                         ! messages in STDERR
          WRITE(unit,'(/,A)')    astline
          WRITE(unit,'(A)') REPEAT(qmark//exmark,6) // ' ERROR ' //     &
                            REPEAT(qmark//exmark,6) // '?'
          WRITE(unit,'(A,A)')  '? Error in routine: ', &
              RoutineName( 1 : Len_Trim(RoutineName) )
          WRITE(unit,'(A,I5)') '? Error Code: ',ErrorStatus
          WRITE(unit,'(A,A)')  '? Error Message: ', &
              Cmessage( 1 : Len_Trim(Cmessage) )
          WRITE(unit,'(A,I5)') '? Error generated from processor: ', mype
          WRITE(unit,'(A,I3,A)') '? This run generated ',numwarnings,' warnings'
          WRITE(unit,'(A,/)')    astline
! DEPENDS ON: um_fort_flush
          CALL Um_Fort_Flush(unit, unset)
        END IF
      END DO
      
      CALL GC_Abort( mype, nproc_max, msg )

    ELSE IF ( ErrorStatus < 0) THEN
      numwarnings = numwarnings + 1
      DO i=1,size(out_units)
        unit=out_units(i)
        IF (unit /= 0) THEN ! don't write warnings to stderr
          WRITE(unit,'(/,A)')    astline
          WRITE(unit,'(A)')  REPEAT("?",35) // ' WARNING ' // REPEAT("?",36)
          WRITE(unit,'(A,A)')  '? Warning in routine: ', &
              RoutineName( 1 : Len_Trim(RoutineName) )
          WRITE(unit,'(A,I5)') '? Warning Code: ',ErrorStatus
          WRITE(unit,'(A,A)')  '? Warning Message: ', &
              Cmessage( 1 : LEN_TRIM(Cmessage) )
          WRITE(unit,'(A,I5)') '? Warning generated from processor: ', mype
          WRITE(unit,'(A,/)')    astline
        ENDIF
      END DO
    END IF
    
    ! Reset ErrorStatus
    ErrorStatus = 0

    IF (lhook) CALL dr_hook('EREPORT64',zhook_out,zhook_handle)

    RETURN
  END SUBROUTINE Ereport64

  ! The 32 bit interface to the above routine
  SUBROUTINE Ereport32 ( RoutineName, ErrorStatus32, Cmessage )

    USE UM_Types, ONLY : &
      integer32,         &
      integer64
    IMPLICIT NONE

! Arguments
    CHARACTER (LEN=*), INTENT( IN )   :: RoutineName
    CHARACTER (LEN=*), INTENT( IN )   :: Cmessage

    INTEGER(KIND=Integer32)           :: ErrorStatus32

! Local variable
    INTEGER(KIND=Integer64)           :: ErrorStatus64

    REAL(KIND=jprb)                   :: zhook_handle
    
    IF (lhook) CALL dr_hook('EREPORT32',zhook_in,zhook_handle)
    
! Copy 32 bit errorstatus to 64 bit temporary
    ErrorStatus64 = ErrorStatus32

! Call the ereport64 routine which does all the work.
    CALL Ereport64(RoutineName, ErrorStatus64, Cmessage)

! Copy the 64 bit errorstatus (can be reset above) to the 32 bit version.
    ErrorStatus32 = ErrorStatus64

    IF (lhook) CALL dr_hook('EREPORT32',zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE Ereport32


! Subroutine to display the number of warnings generated
  SUBROUTINE Ereport_finalise( )
  
    USE UM_ParCore 
    IMPLICIT NONE

    INTEGER                       :: unset = -99
    INTEGER                       :: unit,i
    INTEGER                       :: out_units(2)                     

    CHARACTER (LEN=*), PARAMETER  :: astline = REPEAT("?",80)
    REAL(KIND=jprb)                   :: zhook_handle
    
    DATA out_units /0,6/
    IF (lhook) CALL dr_hook('EREPORT_FINALISE',zhook_in,zhook_handle)

    IF (numwarnings > 0) THEN    
      DO i=1,size(out_units)
        unit=out_units(i)
        IF (mype == 0 .OR. unit /= 0) THEN ! don't get mangled-up 
                                         ! messages in STDERR
          WRITE(unit,'(A)')    astline
          WRITE(unit,'(A)')    '!!'// REPEAT("?",32) // &
                   ' ATTENTION '//REPEAT("?",32) // '!!'
          WRITE(unit,'(A,I3,A)')'? This run generated ',numwarnings,' warnings'
          WRITE(unit,'(A,/)')    astline
! DEPENDS ON: um_fort_flush
          CALL Um_Fort_Flush(unit, unset)
        END IF
      END DO
    END IF
      
    IF (lhook) CALL dr_hook('EREPORT_FINALISE',zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE Ereport_finalise
  

END MODULE Ereport_Mod
