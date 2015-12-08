! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!+
! Subroutine Interface:
      SUBROUTINE LEVSRT(TYPE,NLEVS,IL,RL)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
! Description:
!
! Method:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc
!
!  Code description:
!    FORTRAN 77 + common Fortran 90 extensions.
!    Written to UM programming standards version 7.
!
!  System component covered:
!  System task:               Sub-Models Project

! Subroutine arguments:

!   Scalar arguments with intent(in):

      CHARACTER(LEN=1) TYPE
      INTEGER     NLEVS

!   Array arguments with intent(inout):

      REAL        RL(NLEVS)
      INTEGER     IL(NLEVS)

! Local variables:

      LOGICAL     LSWAP
      INTEGER     I
      INTEGER     J
      INTEGER     ILT
      REAL        RLT

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of Header ----------------------------------------------------

      IF (lhook) CALL dr_hook('LEVSRT',zhook_in,zhook_handle)
      DO I=1,NLEVS
        LSWAP=.FALSE.
        DO J=1,NLEVS-1
          IF(TYPE == 'I') THEN
            IF(IL(J) >  IL(J+1)) THEN
              LSWAP=.TRUE.
              ILT=IL(J)
              IL(J)=IL(J+1)
              IL(J+1)=ILT
            END IF
          ELSE
            IF(RL(J) <  RL(J+1)) THEN
              LSWAP=.TRUE.
              RLT=RL(J)
              RL(J)=RL(J+1)
              RL(J+1)=RLT
            END IF
          END IF
        END DO
        IF(.NOT.LSWAP) THEN 
          GOTO 9999
        END IF
      END DO

 9999 CONTINUE
      IF (lhook) CALL dr_hook('LEVSRT',zhook_out,zhook_handle)
      RETURN
      END SUBROUTINE LEVSRT
