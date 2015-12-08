! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Subroutine DERV_LAND_FIELD : Computes no of land points in MPP jobs

! Subroutine Interface :
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: MPP

SUBROUTINE calc_land_field (unit_no,fixhd,                        &
                            len1_lookup,len2_lookup,              &
                            icode,cmessage)

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE IO
IMPLICIT NONE

!     Arguments
INTEGER, INTENT(IN) ::  unit_no       ! IN Unit Number
INTEGER, INTENT(IN) ::  fixhd(256)    ! IN Fixed header
INTEGER, INTENT(IN) ::  len1_lookup   ! IN First dimension of lookup table
INTEGER, INTENT(IN) ::  len2_lookup   ! IN Seconf dimension of lookup table
INTEGER, INTENT(OUT) :: icode         ! OUT Return code

CHARACTER(LEN=*) cmessage ! OUT Error message

!     Local variables
INTEGER :: len_io        ! length of data returned from buffin
INTEGER :: lookup(len1_lookup,len2_lookup)   !  Lookup table
REAL    :: rcode            ! Real return code

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CALC_LAND_FIELD',zhook_in,zhook_handle)

!     Position atmos dump to read in lookup-table

CALL setpos (unit_no,fixhd(150)-1,icode)

!     Check error code from setpos
IF (icode >  0) THEN
  WRITE (6,*) 'Error in SETPOS called from CALC_LAND_FIELD.'
  WRITE (6,*) 'Trying to point to start of lookup table '//       &
              'in atmos dump.'
  WRITE (cmessage,*) 'DRLANDF1 : Error in SETPOS.'
  GO TO 9999   !  Return
END IF

!     Read in the look-up table

CALL buffin (unit_no,lookup,len1_lookup*len2_lookup,len_io,rcode)

!     Check error code from buffin
IF (rcode /= -1.0) THEN
  WRITE (6,*) 'Error in BUFFIN called from CALC_LAND_FIELD.'
  WRITE (6,*) 'Trying to read lookup table from atmos dump.'
  WRITE (6,*) 'Return code from BUFFIN ',rcode
  icode = 100
  WRITE (cmessage,*) 'DRLANDF1 : Error in BUFFIN.'
  GO TO 9999   !  Return
END IF

!     Read in land-sea mask and then
!     compute the number of land points for each PE
! DEPENDS ON: read_land_sea
CALL read_land_sea (unit_no,rcode,lookup,len1_lookup,len2_lookup, &
                    fixhd,256)

!     Check error code from read_land_sea
IF (rcode /= -1.0) THEN
  WRITE (6,*) 'Error in READ_LAND_SEA.'
  WRITE (6,*) 'Return code from READ_LAND_SEA ',rcode
  icode = 200
  WRITE (cmessage,*) 'DRLANDF1 : Error in READ_LAND_SEA.'
  GO TO 9999   !  Return
END IF

 9999 CONTINUE
IF (lhook) CALL dr_hook('CALC_LAND_FIELD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE calc_land_field
