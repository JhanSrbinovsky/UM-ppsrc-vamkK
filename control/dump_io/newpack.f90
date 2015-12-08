! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE NEWPACK:--------------------------------------
!
!    Purpose: Packing codes stored in LOOKUP(21,K) & LOOKUP(39,K)
!             are changed from pre vn2.8 values to
!             specification required at release 2.8
!
!
!    Documentation: UM Documentation Paper F3
!
!      -----------------------------------------------------------------
!
!      Code Owner: See Unified Model Code Owners HTML page
!      This file belongs in section: Dump I/O

SUBROUTINE newpack                                                &
(lookup,len1_lookup,len2_lookup)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

INTEGER                                                           &
 len1_lookup                                                      &
,len2_lookup                                                      &
,lookup(len1_lookup,len2_lookup)

INTEGER                                                           &
 n1                                                               &
,n2                                                               &
,n3                                                               &
,k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('NEWPACK',zhook_in,zhook_handle)
DO k=1,len2_lookup
  n1=0
  n2=0
  n3=0
  IF (lookup(21,k) == -2) n1=2
! Ocean field packed using index array
  IF (lookup(21,k) >  9 .AND. lookup(21,k) <  100) THEN
    n2=1
    n3=lookup(21,k)-10
  END IF
! Ocean field compressed using bit mask
  IF (lookup(21,k) > 99) THEN
    n2=2
    n3=lookup(21,k)-100
  END IF
! Real field stored at land pts
  IF (lookup(39,k) == 4) THEN
    lookup(39,k)=1
    n2=2
    n3=1
  END IF
! Integer field stored at land pts
  IF (lookup(39,k) == 5) THEN
    lookup(39,k)=2
    n2=2
    n3=1
  END IF

  lookup(21,k)=100*n3+10*n2+n1
END DO

IF (lhook) CALL dr_hook('NEWPACK',zhook_out,zhook_handle)
RETURN
END SUBROUTINE newpack
