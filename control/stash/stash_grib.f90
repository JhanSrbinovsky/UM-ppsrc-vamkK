! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  SUBROUTINE STASH_GRIB---------------------------------------------
!
!  Purpose:
!   STASH_GRIB is a subroutine to code the stash section number and
!   parameter value in elements of the grib header.
!
!   octet 4 of the grib product definition section is the version
!   number of the table 2 (parameter code description)
!   values from 128 to 254 areavailable for local use, and we
!   use them to describe the stash section number of the field. for
!   each stash section number there are two octet 4 values. the first
!   is for stash parameter values from 0 to 255, the second for values
!   256 to 511.
!   octet 9 is the code value in table 2, ie stash parameter value, or
!   stash parameter value -256 if it is more than 255.
!
!  Programming standard: Unified Model Documentation Paper No 3
!
!  System component:
!
!  System task:
!
!  Documentation:
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH
      SUBROUTINE STASH_GRIB(STASH_SECTION_NUMBER,STASH_ITEM_NUMBER,     &
     &                      GRIB_BLOCK1_OCTET4,GRIB_BLOCK1_OCTET9,      &
     &                      ERROR)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim

      IMPLICIT NONE

      INTEGER                                                           &
     &   STASH_SECTION_NUMBER                                           &
                              ! STASH SECTION NUMBER         INPUT
     &  ,STASH_ITEM_NUMBER                                              &
                              ! STASH PARAMETER VALUE        INPUT
     &  ,GRIB_BLOCK1_OCTET4                                             &
                              ! OCTET 4 FROM GRIB PDB        OUTPUT
     &  ,GRIB_BLOCK1_OCTET9                                             &
                              ! OCTET 9 FROM GRIB PDB        OUTPUT
     &  ,ERROR                ! ERROR OUTPUT CODE            OUTPUT
!     LOCAL VARIABLES
      INTEGER                                                           &
     &   CARRY   ! CARRY VALUE FROM ODD VALUES OF GRIB_BLOCK1_OCTET4

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle
!****
      IF (lhook) CALL dr_hook('STASH_GRIB',zhook_in,zhook_handle)
      IF(STASH_ITEM_NUMBER >  511.OR.STASH_ITEM_NUMBER <  0) THEN
        ERROR = 999
        IF (lhook) CALL dr_hook('STASH_GRIB',zhook_out,zhook_handle)
        RETURN
      ELSE IF(STASH_ITEM_NUMBER >  255) THEN
        CARRY = 1
        GRIB_BLOCK1_OCTET9 = STASH_ITEM_NUMBER - 256
      ELSE
        CARRY = 0
        GRIB_BLOCK1_OCTET9 = STASH_ITEM_NUMBER
      ENDIF
      IF((STASH_SECTION_NUMBER >= 0).AND.                               &
     &   (STASH_SECTION_NUMBER <= 62)) THEN
        GRIB_BLOCK1_OCTET4 = STASH_SECTION_NUMBER*2 + 128 + CARRY
      ELSE
        ERROR = 999
      ENDIF
      IF (lhook) CALL dr_hook('STASH_GRIB',zhook_out,zhook_handle)
      RETURN
!****
      END SUBROUTINE STASH_GRIB
