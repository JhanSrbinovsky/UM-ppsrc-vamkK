! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE PR_LOOK----------------------------------------
!
!    Purpose: Prints out Kth 64-word PP header
!
!    Programming standard:  Unified Model Documentation Paper No 3
!
!    Documentation:  Unified Model Documentation Paper No F3
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Dump I/O

SUBROUTINE pr_look(                                               &
                   lookup,rlookup,len1_lookup,k)





USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE lookup_addresses

IMPLICIT NONE

INTEGER len1_lookup     ! IN First dimension of Look Up Table
INTEGER k               ! IN Field number in Look Up Table
INTEGER                                                           &
 lookup(len1_lookup,*)  ! IN Integer equivalence of PP LOOKUP
REAL                                                              &
 rlookup(len1_lookup,*) ! IN Real equivalence of PP LOOKUP

CHARACTER(LEN=36) exppxc

! Local variables:---------------------------------------------
INTEGER icode             !Error code
INTEGER item              !STASH item number
INTEGER section           !STASH section number
INTEGER model             !Internal model number
INTEGER i                 !Index
INTEGER lsec,lsecd        !local seconds values (0 if lbrel<3)

CHARACTER(LEN=36) phrase       !Character part of PPXREF record
CHARACTER(LEN=80) cmessage     !Error message





INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------------

!  Write time and field type
IF (lhook) CALL dr_hook('PR_LOOK',zhook_in,zhook_handle)
item=MOD(lookup(42,k),1000)
section=(lookup(42,k)-item)/1000
model=lookup(45,k)

!       All diagnostics under model code of 10 are in section 20
!       of Atmos StashMaster file.
IF ( model == 10 ) THEN
  model = 1
END IF
icode = 0
! DEPENDS ON: exppxc
phrase=exppxc(model,section,item,                                 &
              icode,cmessage)
IF (icode /= 0) THEN
  phrase='NON-STANDARD FIELD'
END IF

WRITE (6,'('' FIELD NO.'', i5,4x,a)') k,phrase

! For lookup release < 3, lbsec(d) held day_number(s), so output zero

IF (lookup(lbrel,k) < 3) THEN
  lsec = 0
  lsecd = 0
ELSE
  lsec = lookup(lbsec,k)
  lsecd = lookup(lbsecd,k)
END IF

WRITE (6,'('' VALID AT: '',2(i2.2,'':''),i2.2,''Z  '',2(i2.2,''/''),'//&
    'i4.4,''    DATA TIME: '',2(i2.2,'':''),i2.2,''Z  '',2(i2.2,''/''),'//&
    'i4.4)') lookup(lbhr,k),lookup(lbmin,k),lsec,  &
    lookup(lbdat,k),lookup(lbmon,k),lookup(lbyr,k),  &
    lookup(lbhrd,k),lookup(lbmind,k),lsecd,          &
    lookup(lbdatd,k),lookup(lbmond,k),lookup(lbyrd,k)

!  Rest of header
WRITE(6,'(''   LBTIM   LBFT    LBLREC LBCODE  LBHEM  LBROW'',         '//&
    ' ''  LBNPT  LBEXT LBPACK'',/, 1x, 2i7, i10, 6i7,/,               '//&
    ' ''   LBREL   LBFC  LBCFC LBPROC   LBVC  LBRVC  LBEXP'',         '//&
    ' ''   LBBEGIN    LBNREC'',/, 1x, 7i7, 2i10,/,                    '//&
    ' ''  LBPROJ  LBTYP  LBLEV LBRSVD LBRSVD LBRSVD LBRSVD   LBSRCE'','//&
    ' /,1x,7i7,i9,/,                                                  '//&
    ' ''  DATA_TYPE     NADDR    LBUSER ITEM_CODE    LBPLEV'',        '//&
    ' ''    LBUSER MODEL_CODE'',/, 1x, 6i10, i11,/,                   '//&
    ' 9x,''BULEV'',7x,''BHULEV'',5x,''BRSVD(3)'',5x,''BRSVD(4)'',     '//&
    ' 7x,''BDATUM'',/, 1x, 1p, 5e13.4,/,                              '//&
    ' 10x,''BACC'',9x,''BLEV'',8x,''BRLEV'',8x,''BHLEV'',7x,''BHRLEV'''//&
    ' ,/, 1x, 1p, 5e13.4, /,                                          '//&
    ' 9x,''BPLAT'',8x,''BPLON'',9x,''BGOR'',10x,''BZY'',10x,''BDY'',/,'//&
    ' 1x, 1p, 5e13.4,/,11x,''BZX'',10x,''BDX'',9x,''BMDI'',9x,''BMKS'''//&
    ' ,/,1x,1p,4e13.4)') (lookup(i,k),i=13,45),(rlookup(i,k),i=46,64)

WRITE(6,'('' '')')

IF (lhook) CALL dr_hook('PR_LOOK',zhook_out,zhook_handle)
RETURN
END SUBROUTINE pr_look

