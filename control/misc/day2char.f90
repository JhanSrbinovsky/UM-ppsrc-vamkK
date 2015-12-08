! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Subroutine : DAY2CHAR  -------------------------------------------
!LL
!LL  Purpose: Convert days to the character to represent the period.
!LL
!LL  Programming standard: UM Doc Paper 3, version 1 (15/1/90)
!LL
!LL  Logical components covered: S51
!LL
!LL  Project task: S51
!LL
!LL  External documentation: UM documentation paper 7 - Filenaming
!LL                          conventions for the Unified Model
!LL
!LLEND -----------------------------------------------------------------
!*L  Interface and arguments: ------------------------------------------
!
!
!Code Owner: See Unified Model Code Owners HTML page
!This file belongs in section: Misc
      SUBROUTINE DAY2CHAR(NDAYS,DAYCHAR)
!

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      IMPLICIT NONE
!
      INTEGER       NDAYS        ! IN  - number of days in period
      CHARACTER(LEN=1)   DAYCHAR      ! OUT - character for period

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!
!*----------------------------------------------------------------------
!  Common blocks
!
!
!  Local variables
!
!
!  IF period a multiple of years
!
      IF (lhook) CALL dr_hook('DAY2CHAR',zhook_in,zhook_handle)
        IF (MOD(NDAYS,360) == 0)THEN
         IF (NDAYS == 360) THEN        ! 1 year mean
           DAYCHAR='y'
         ELSE IF (NDAYS == 1800) THEN  ! 5 year means
           DAYCHAR='v'
         ELSE IF (NDAYS == 3600) THEN  ! 10 year means
           DAYCHAR='x'
         ELSE IF (NDAYS == 18000) THEN ! 50 year means
           DAYCHAR='l'
         ELSE IF (NDAYS == 36000) THEN ! 100 year means
           DAYCHAR='u'
         ELSE IF (NDAYS == 360000) THEN ! 1000 year means
           DAYCHAR='z'
         ELSE                           ! not a special period
           DAYCHAR='0'
         ENDIF
!
        ELSE
! periods less than one year
!
         IF (NDAYS == 5) THEN          ! 5 days means
           DAYCHAR='p'
         ELSE IF (NDAYS == 7) THEN     ! weekly means
           DAYCHAR='w'
         ELSE IF (NDAYS == 10) THEN    ! 10 day means
           DAYCHAR='t'
         ELSE IF (NDAYS == 14) THEN    ! fortnightly means
           DAYCHAR='r'
         ELSE IF (NDAYS == 30) THEN    ! monthly means
           DAYCHAR='m'
         ELSE IF (NDAYS == 90) THEN    ! seasonal means
           DAYCHAR='s'
         ELSE                          ! not a special period
           DAYCHAR='0'
         ENDIF
        ENDIF
!
        IF (lhook) CALL dr_hook('DAY2CHAR',zhook_out,zhook_handle)
        RETURN
!L----------------------------------------------------------------------
      END SUBROUTINE DAY2CHAR
