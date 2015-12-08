MODULE OASIS3_split_comm_mod
! *****************************COPYRIGHT****************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *****************************COPYRIGHT****************************************
IMPLICIT NONE

CONTAINS

SUBROUTINE OASIS3_split_comm(colour, atm_coloured)
!
! Description: This routinere defines the communicator used by OASIS
!              (specifically OASIS3-MCT) when excluding certain
!              processes from coupling (typically IO server procs.)  
!
! Code Owner: See Unified Model Code Owners HTML page 
! This file belongs in section: Coupling 
!==================================================================

USE oasis3_atmos_init_mod, ONLY: PRISM_create_couplcomm
 
USE mpl, ONLY: mpl_undefined
USE oasis_atm_data_mod, ONLY: comm_in, icpl
USE um_types
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

INTEGER :: colour        ! The "colour" of this process defined by IOS
                         ! (or potentailly any other sub-division of
                         ! processes not involved in coupling. 
INTEGER :: atm_coloured  ! The standard "colour" of atmosphere processes

! Local variables...
INTEGER (KIND=integer32) :: ierror ! Return code for OASIS3 calls

INTEGER (KIND=integer32) :: comm_in32 ! 32 bit version of comm_in
                                      ! - the communicator defined
                                      ! by prism/oasis_get_localcomm 
INTEGER (KIND=integer32) :: couplcomm ! OASIS communicator 
                                           
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('OASIS3_split_comm',zhook_in,zhook_handle)

! We need to create a special communicator for OASIS which excludes the 
! IOS processes. 
! Set default to "not involved in coupling". It's critical that this value
! matches what OASIS is expecting in terms of MPI_UNDEFINED. 
icpl = mpl_undefined

! Set all atmos-only procs to be involved in coupling
IF ( colour == atm_coloured ) THEN
  icpl = 1
END IF

! We need a 32 bit version of the original global communicator
comm_in32 = comm_in

! Tell OASIS to create its special communicator. In principle we could
! use an existing communicator and call PRISM_set_couplecomm
! instead, but nothing suitable is available in the UM
! at this point. The point being that we need a communicator
! in which some of the processors (the IOS ones) are undefined, but this
! and certain other oasis initialisation calls are collective.  
CALL PRISM_create_couplcomm(icpl, comm_in32, couplcomm, ierror)

IF (lhook) CALL dr_hook('OASIS3_split_comm',zhook_out,zhook_handle)

RETURN

END SUBROUTINE OASIS3_split_comm

END MODULE OASIS3_split_comm_mod
