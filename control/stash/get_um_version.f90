! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Function read the UM version number from environment variable.
!
! Function interface:
      INTEGER FUNCTION GET_UM_VERSION()
!
! Description:
!   Function reads environment variable $VN and returns the UM version
!   number as an integer
!
! Method:
!   Reads the two components of the UM version number from environment
!   variable $VN and combine them into a four digit integer.
!   e.g If $VN = 5.3, function returns 503.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH
!
!
! Code Description:
!   Language: FORTRAN 77 + some CRAY extensions
!   This code is written to UMDP3 v6 programming standards.
!

      USE parkind1, ONLY: jprb, jpim
      USE yomhook, ONLY: lhook, dr_hook
      USE ereport_mod, ONLY : ereport
      IMPLICIT NONE

! Declarations:

! Local parameters:
      INTEGER, PARAMETER :: env_name_length = 2
      INTEGER, PARAMETER :: env_length      = 8
      INTEGER, PARAMETER :: no_version      = -9999

      CHARACTER (LEN=*), PARAMETER :: routinename ='get_um_version'
      CHARACTER (LEN=*), PARAMETER :: env_name    = 'VN'

! Local scalars:
      INTEGER ::                                                        &
     &  icode                                                           &
                                ! Return code  =0 Normal exit  >1 Error
     &, um_version                                                      &
                                ! Version of UM
     &, um_revision             ! Revision of UM
      CHARACTER (LEN=8) ::                                              &
     & c_um_version             ! UM version as string
      CHARACTER (LEN=80) ::                                             &
     & cmessage                 ! Error message

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of header
 
!    Get the UM version from the environment variable $VN.
      IF (lhook) CALL dr_hook('GET_UM_VERSION',zhook_in,zhook_handle)
      CALL FORT_GET_ENV(env_name, env_name_length, c_um_version,        &
     &                  env_length, icode)

      IF (icode  /=  0) THEN     ! $VN was not set

         WRITE (6,*)                                                    &
     &        'GET_UM_VERSION : WARNING : Environment variable VN not ',&
     &        'set or not obtainable; skipping version checking.'

         cmessage = 'Environment variable VN not set, no version '//    &
     &        'checking performed'
         get_um_version = no_version

      ELSE
         READ (c_um_version, '(i1,1x,i1)') um_version, um_revision
         get_um_version = um_version*100 + um_revision
      END IF

      IF (icode /= 0) THEN
        CALL Ereport(RoutineName,icode,Cmessage)
      END IF

      IF (lhook) CALL dr_hook('GET_UM_VERSION',zhook_out,zhook_handle)
      RETURN
      END FUNCTION GET_UM_VERSION
