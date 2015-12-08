! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! SUBROUTINE COEX
!
! PURPOSE:   PACK TO AND UNPACK FROM WGDOS FORMAT
!
! -------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: STASH

SUBROUTINE coex(field,m,icomp,n,ix,iy,num,isc,oco,rmdi,lword,     &
                icode,cmessage)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Subroutine arguments
INTEGER :: n
INTEGER :: m
INTEGER :: ix
INTEGER :: iy
INTEGER :: num
INTEGER :: isc
INTEGER :: lword
INTEGER :: icode
INTEGER :: icomp(n)
LOGICAL :: oco
REAL    :: field(m)
REAL    :: rmdi
CHARACTER :: cmessage*80

! Local variables
INTEGER :: ier
INTEGER :: ocode

! External functions used
INTEGER :: ieee2ibm
INTEGER :: ibm2ieee

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle



!                     OCO=.TRUE.                 OCO=.FALSE.

!      FIELD   =  FIELD FOR COMPRESSING     RETURNED EXPANDED DATA
!          M   =  SIZE OF FIELD             SIZE OF FIELD
!      ICOMP   =  RETURNED COMPRESSED DATA  FIELD FOR EXPANSION
!          N   =  SIZE OF COMP                 -------
!         IX   =  X DIMENSION OF FIELD      RETURNED X DIMENSION
!         IY   =  Y DIMENSION OF FIELD      RETURNED Y DIMENSION
!        NUM   =  TOTAL NO. OF COMPRESSED      -------
!                 (32 BIT) WORDS RETURNED
!        ISC   =  ACCURACY IN POWER OF 2
!        OCO   =  .TRUE. FOR COMPRESSION    .FALSE. FOR EXPANSION

!        USERS OF THIS ROUTINE SHOULD ENSURE THAT THE ARRAY 'COMP'
!      WILL BE BIG ENOUGH TO HOLD THE COMPRESSED DATA.
!        IF INPUT ARRAY IS ONE DIMENSIONAL PUT IY=1 AND IX=DIMENSION,
!      WHEN USING WITH PRINTFILE FIELDS USE IX=192,IY=121  NOT IX=23232,
!      IY=1 THIS WILL MAKE THE COMPRESSION MORE EFFICIENT.
!      FOR MOST PRINTFILE FIELDS USING AN ACCURACY OF AN EIGHTH (ISC=-3)
!      A COMPRESSION FACTOR OF 4 CAN BE ACHIEVED. FOR DDDFFF FIELDS
!      USERS ARE ADVISED TO SPLIT THE FIELD INTO U AND V
!      COMPONENTS.

!     CRAY ROUTINE     64 BIT WORDS - CRAY FULLWORD

! Check if an error has already been encountered, and get out
! if it has.
IF (lhook) CALL dr_hook('COEX',zhook_in,zhook_handle)
ocode = 0
IF (icode  >   0) THEN
   GO TO 9999
ELSE IF (icode  <   0)THEN
   ocode = icode
   icode = 0
END IF

IF (oco) THEN
! COMPRESSION OF DATA
! INSERT SCALE FACTOR, COLS AND ROWS INTO HEADER

    ier=ieee2ibm(2,1,icomp(1),32,isc,1,64,32)
    ier=ieee2ibm(2,1,icomp(2),0,ix,1,64,16)
    ier=ieee2ibm(2,1,icomp(2),16,iy,1,64,16)


! Well vectorized version compressing the whole field at once
! DEPENDS ON: cmps_all
    CALL cmps_all(field,icomp,n,ix,iy,num,isc,rmdi,icode,cmessage)

! END OF COMPRESSION SECTION

ELSE

! EXPANSION SECTION
! EXTRACT SCALE FACTOR, COLS AND ROWS FROM HEADER

    ier=ibm2ieee(2,1,icomp(1),32,isc,1,64,32)
    ier=ibm2ieee(2,1,icomp(2),0,ix,1,64,16)
    ier=ibm2ieee(2,1,icomp(2),16,iy,1,64,16)


! Vectorized version expanding the whole field at once
! This version also has OpenMP parallelisation
! DEPENDS ON: xpnd_all
    CALL xpnd_all(field,icomp,n,ix,iy,isc,rmdi,                   &
                  icode,cmessage)

END IF


IF (icode  ==  0 .AND. ocode  /=  0) THEN
   icode = ocode
END IF

 9999 CONTINUE

IF (lhook) CALL dr_hook('COEX',zhook_out,zhook_handle)
RETURN

END SUBROUTINE coex
