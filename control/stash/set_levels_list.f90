! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SET_LEVELS_LIST
!
! Purpose : To set up a list of levels at which a diagnostic is
!           required, using information in the STASH list.
! Service routine  version for Cray YMP
!
! Programming Standard : Unified Model Documentation paper number 4
!                      : Version no 2, dated 18/01/90
!
! System components covered : D3
!
! System task : P0
!
! Documentation: U.M. Documentation paper number P0,C4
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH

      SUBROUTINE SET_LEVELS_LIST(LEVELS,                                &
     &                    LEN_STLIST,STLIST,LIST,STASH_LEVELS,          &
     &      LEN_STASHLEVELS,ICODE,CMESSAGE)


      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE

      INTEGER                                                           &
     &       LEVELS,                                                    &
                          ! IN Number of levels in input data
     &       LEN_STLIST,                                                &
                          ! IN
     &       STLIST(LEN_STLIST),                                        &
                                 ! IN STASH list
     &  LEN_STASHLEVELS,                                                &
     &  STASH_LEVELS(LEN_STASHLEVELS,*),                                &
                                          ! IN - list of levels required
     &       ICODE        ! OUT Return code =0 Normal exit
!                                       >1 Error message

      LOGICAL                                                           &
     &       LIST(LEVELS) ! OUT List of levels required.

      CHARACTER(LEN=80) CMESSAGE ! Error message

!L Local variables

      INTEGER                                                           &
     &      K,                                                          &
     &      KOUT

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle


!L Initialise levels list to false

      IF (lhook) CALL dr_hook('SET_LEVELS_LIST',zhook_in,zhook_handle)
      DO K=1,LEVELS
        LIST(K)= .FALSE.
      END DO

!L Check for method of levels selection
!L Levels list must be present.

      IF(STLIST(10) <  0) THEN

! Set logical array list to identify levels required.

        DO KOUT=2,STASH_LEVELS(1,-STLIST(10))+1
          IF((STASH_LEVELS(KOUT,-STLIST(10)) >= 1).AND.                 &
     &    (STASH_LEVELS(KOUT,-STLIST(10)) <= LEVELS)) THEN
!         LEVEL IS IN THE RANGE OF LIST.
              LIST(STASH_LEVELS(KOUT,-STLIST(10))) =.TRUE.
          ELSE
!         LEVEL IS OUT OF THE RANGE OF LIST.
              CMESSAGE=  ' SET_LEVELS_LIST: level out of range'
              WRITE(6,*) ' SET_LEVELS_LIST: level out of range'
              WRITE(6,*) ' level=',STASH_LEVELS(KOUT,-STLIST(10))
              WRITE(6,*) ' Section, Item =',STLIST(2),STLIST(1)
              ICODE=2
          END IF
        END DO


      ELSE if (STLIST(10) /= 100) THEN

!L Set list of levels according to its definition in stlist :
!l If stlist(10) positive, input on range of model levels starting
!l at level stlist(10), finishing at level stlist(11)

         DO K=stlist(10),stlist(11)
            LIST(K)= .true.
         END DO


      ELSE

!L Illegal control data

        ICODE=1
        CMESSAGE='SET_LEVELS_LIST: Illegal control data'
      WRITE(6,*) 'Illegal control data SET_LEVELS_LIST,STLIST(10,11)=', &
     &         STLIST(10) ,STLIST(11)
      WRITE(6,*) 'Section and item numbers ',STLIST(2),STLIST(1)
      IF (lhook) CALL dr_hook('SET_LEVELS_LIST',zhook_out,zhook_handle)
      RETURN

      END IF

      IF (lhook) CALL dr_hook('SET_LEVELS_LIST',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE SET_LEVELS_LIST
