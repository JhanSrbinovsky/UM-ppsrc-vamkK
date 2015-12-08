! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!     SUBROUTINE CMPS_ALL -----------------
!
!     PURPOSE:   PACK TO AND UNPACK FROM WGDOS FORMAT
!
!
!     PROGRAMMING STANDARD: UNIFIED MODEL DOCUMENTATION PAPER NO. 3,
!
!    Logical component number: S72
!
!     SYSTEM TASK: P7
!     -------------------------------------------------------------
!
!     Code Owner: See Unified Model Code Owners HTML page
!     This file belongs in section: STASH
SUBROUTINE cmps_all(field,icomp64,n,ix,iy,num,isc,rmdi,           &
                    icode,cmessage)

!     Better vectorized version of CMPS compressing the whole field
!     at once.

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE ereport_mod, ONLY : ereport
IMPLICIT NONE

!     Subroutine arguments

INTEGER, INTENT(IN)  :: n, ix, iy, isc
INTEGER, INTENT(OUT) :: num
REAL,    INTENT(IN)  :: field(ix,iy)    ! Input uncompressed data
REAL,    INTENT(IN)  :: rmdi
INTEGER, INTENT(OUT) :: icomp64(n)      ! Compressed field
INTEGER              :: icode           ! Not returned
CHARACTER            :: cmessage*80     ! Not returned

!     Local variables

INTEGER :: i, j, npoint, nshft, ival
INTEGER :: isig, iexp, iman
INTEGER :: is1, is2, is3
INTEGER :: nbits_bmap, nwords_bmap, nwords_data
INTEGER :: nbits_pack, nvals_pack
INTEGER :: i1, i2, j1, j2

REAL    :: aprec, bprec

LOGICAL :: obtmis, obtzer

INTEGER :: ibase(iy)              ! Base values
INTEGER :: imax(iy)               ! Max values
INTEGER :: ibit(iy)
INTEGER :: nmiss(iy)              ! missing-data bitmaps
INTEGER :: nzero(iy)              ! zero bitmaps
INTEGER :: ibm(iy)                ! IBM representation of base
INTEGER :: itmp(2*ix+32)          ! temporary storage
INTEGER :: icomp(2*ix+8,MAX(iy,1))! temporary storage for compression
INTEGER :: iword(iy)              ! per row sizes
INTEGER :: istart(iy)             ! start position in a row for data

REAL    :: atmp(ix)
REAL    :: base(iy)
REAL    :: fmax(iy)

CHARACTER(LEN=*), PARAMETER :: routinename='cmps_all'
INTEGER, PARAMETER :: mask32 = z'FFFFFFFF'
INTEGER, PARAMETER :: mask_expt_ibm = z'7F000000'
INTEGER, PARAMETER :: mask_sign_ibm = z'80000000'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook('CMPS_ALL',zhook_in,zhook_handle)

! GENERAL REMARK:
! All lengths and alignments in WGDOS packing are for 32-bit words,
! so life gets much easier when we treat the packed data as 32-bit
! words.
! So we gather all compressed data in the low 32 bits of icomp
! and compress this array at the end to 64 bits

! Scale factor
aprec = 2.**isc
bprec = 1./aprec

! Parallelisation is over rows - these can be be compressed
! independently of each other and combined into a single buffer
! at the end
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& PRIVATE(j, i, i2, isig, iexp, iman, itmp, obtmis, obtzer,        &
!$OMP&         nbits_bmap, npoint, is1, is2, is3, nshft, nvals_pack,    &
!$OMP&         nbits_pack, nwords_data, ival, atmp, nwords_bmap,        &
!$OMP&         icode, cmessage)                                         &
!$OMP& SHARED(iy, base, fmax, nmiss, nzero, field, ix, rmdi, ibase, ibm,&
!$OMP&        imax, ibit, aprec, bprec, istart, iword, icomp)
DO j=1,iy

!     Find minimum and maximum value for every row,
!     count number of missing and zero numbers
  base(j) = HUGE(0.0)
  fmax(j) = -HUGE(0.0)
  nmiss(j) = 0
  nzero(j) = 0
  DO i=1,ix
    IF(field(i,j)/=rmdi) THEN
      base(j) = MIN(base(j),field(i,j))
      fmax(j) = MAX(fmax(j),field(i,j))
    ELSE
      nmiss(j) = nmiss(j)+1
    END IF
  END DO
! If we have a row of rmdi then lets set fmax and base to -1.0.  This is not
! really defined anywhere - we want something sensible though but not rmdi or
! 0.0 since they are special already.
  IF (nmiss(j) == ix) THEN
    fmax(j) = -1.0
    base(j) = -1.0
  END IF

!     ROUND BASE TO PRECISION REQUIRED
  ibase(j) = NINT(base(j)*bprec)
  base(j)  = ibase(j)*aprec

!     IBM floating point representation of base:
  IF(base(j)==0.) THEN
    ibm(j) = 0
  ELSE
    isig = SIGN(1.,base(j))
    base(j) = ABS(base(j))
    iexp = EXPONENT(base(j)) + 256
    iman = FRACTION(base(j)) * 16777216.
    i = MOD(iexp,4)
    IF(i==1) iman = ISHFT(iman,-3)
    IF(i==2) iman = ISHFT(iman,-2)
    IF(i==3) iman = ISHFT(iman,-1)
    iexp = (iexp+3)/4
    ibm(j) = IOR(IAND(ISHFT(iexp,24),mask_expt_ibm),iman)
    IF(isig<0) ibm(j)=IOR(ibm(j),mask_sign_ibm)
  END IF

!     Find maximum scaled value
  imax(j) = 0
  IF(nmiss(j)<ix) THEN
    imax(j) = MAX(NINT(fmax(j)*bprec)-ibase(j), 0)
  END IF


!     FIND NUMBER OF BITS REQUIRED TO STORE MAX DIFFERENCE
  ibit(j) = 0
  i2 = 1
  DO i = 1,32
    IF (imax(j) >= i2) ibit(j) = i
    i2 = i2*2
  END DO


  IF (imax(j) > 2147483647) THEN
    icode = 2
    cmessage = 'COEX: Unable to WGDOS pack to this accuracy'
    CALL ereport(routinename, icode, cmessage)
  END IF

! Now pack the non-MDI data.
  npoint = 0
  DO i = 1, ix
    IF (field(i,j) /= rmdi) THEN
      npoint = npoint + 1
      ! This may pack the data to 0.0 which would have been missed before.
      atmp(i) = NINT(field(i,j)*bprec)
      IF (atmp(i) == 0.0) nzero(j) = nzero(j) + 1
    ELSE
      atmp(i) = rmdi
    END IF
  END DO

!     Fill data into output array
  istart(j) = 1

  ! iword is the addressing (number of words) for this row
  iword(j) = istart(j)+2 ! the two header words are inserted at end

  ! Check which bitmaps we need
  obtmis = (nmiss(j)>0)

  ! Check if it is worthwile to use zero-bitmap:
  obtzer = (ibit(j)*nzero(j) > ix)

  ! Set itmp with the bitmap pattern
  nbits_bmap = 0
  IF (obtmis) THEN
    itmp(1:ix) = MERGE(1,0,atmp(:)==rmdi)
    nbits_bmap = ix
  END IF

  ! The bitmap is actually non-zero rather than zeroes.
  IF (obtzer) THEN
    itmp(nbits_bmap+1:nbits_bmap+ix) = MERGE(1,0,atmp(:)/=0.0)
    nbits_bmap = nbits_bmap+ix
  END IF

  ! Insert bitmap - this is done row-by-row into icomp.
  IF (nbits_bmap>0) THEN

    ! add 1's to the end since bitmaps should be padded with 1's
    itmp(nbits_bmap+1:nbits_bmap+31) = 1

    ! The number of words to be used for bitmap
    nwords_bmap = (nbits_bmap+31)/32

    ! Compress itmp

    ! Combine 4 contiguous 1-bit-items to one 4-bit-item
    DO i=1,nwords_bmap*8
      itmp(i) = IOR(IOR(ISHFT(itmp(4*i-3),3),                     &
                        ISHFT(itmp(4*i-2),2)),                    &
                    IOR(ISHFT(itmp(4*i-1),1),itmp(4*i)))
    END DO

    ! Combine 4 contiguous 4-bit-items to one 16-bit-item
    DO i=1,nwords_bmap*2
      itmp(i) = IOR(IOR(ISHFT(itmp(4*i-3),12),                    &
                        ISHFT(itmp(4*i-2), 8)),                   &
                    IOR(ISHFT(itmp(4*i-1), 4),itmp(4*i)))
    END DO

    ! Combine 2 contiguous 16-bit-items to the final destination
    DO i=1,nwords_bmap
      icomp(iword(j)+i-1,j) = IOR(ISHFT(itmp(2*i-1),16),itmp(2*i))
    END DO

    iword(j) = iword(j) + nwords_bmap

  END IF

  ! Insert data

  IF (ibit(j)>0) THEN

    ! Get rid of missing values
    IF (obtmis) THEN
      npoint = 0
      DO i=1,ix
        IF (atmp(i)/=rmdi) THEN
          npoint = npoint+1
          atmp(npoint) = atmp(i)
        END IF
      END DO
    ELSE
      npoint = ix
    END IF

    ! Now get rid of zero values
    IF (obtzer) THEN
      ival = npoint
      npoint = 0
!CDIR NODEP
      DO i=1,ival
        IF(atmp(i)/=0.) THEN
          npoint = npoint+1
          atmp(npoint) = atmp(i)
        END IF
      END DO
    END IF

    ! Number of words used for the compressed data
    nwords_data = (npoint*ibit(j)+31)/32

    ! Use scaled value and find difference from base
    DO i=1,npoint
      itmp(i) = MAX(NINT(atmp(i))-ibase(j), 0)
    END DO

    ! As long as ibit(j) is <=16 we can combine two contiguous
    ! items to one with the double number of bits, halving the
    ! number of words to be compressed
    nbits_pack = ibit(j)
    nvals_pack = npoint

    DO WHILE (nbits_pack <= 16)
      itmp(nvals_pack+1) = 0 ! for odd numbers
      DO i=1,(nvals_pack+1)/2
        itmp(i) = IOR(ISHFT(itmp(2*i-1),nbits_pack),itmp(2*i))
      END DO
      nbits_pack = 2*nbits_pack
      nvals_pack = (nvals_pack+1)/2
    END DO

    IF (nbits_pack == 32) THEN
      ! This is the case if ibit(j) is 1, 2, 4, 8 or 16
      ! We have not much to do, just copy itmp
      DO i=1,nwords_data
        icomp(iword(j)+i-1,j) = itmp(i)
      END DO

    ELSE

      ! Shift every value in itmp to the left and append
      ! the bits of the following 2 words
      is1 = 64-nbits_pack   ! amount to shift itmp(i)
      is2 = 64-2*nbits_pack ! amount to shift itmp(i+1)
      is3 = 64-3*nbits_pack ! amount to shift itmp(i+2)

      itmp(nvals_pack+1) = 0
      itmp(nvals_pack+2) = 0

      DO i=1,nvals_pack
        itmp(i) = IOR(IOR(ISHFT(itmp(i  ),is1),                   &
                          ISHFT(itmp(i+1),is2)),                  &
                          ISHFT(itmp(i+2),is3))
      END DO

      ! Now itmp contains enough data so that we can cut out
      ! the compressed data words
      DO i=1,nwords_data

        ! Word which contains compressed data word:
        ival = itmp(((i-1)*32)/nbits_pack + 1)

        ! Number of bits we have to shift to the left
        ! so that we have the compressed data word left packed:
        nshft = MOD((i-1)*32,nbits_pack)

        ival = ISHFT(ival,nshft)

        ! Shift to the right half and mask out upper bits
        ! (for the case that ISHFT does an arithmetic shift)
        icomp(iword(j)+i-1,j) = IAND(ISHFT(ival,-32),mask32)

      END DO
    END IF

    iword(j) = iword(j) + nwords_data

  END IF

  ! Now insert the header for this row:
  ! First word of compressed data: IBM representation of base
  ! Second word of compressed data:
  ! 16 bits: ibit(j) + flags
  ! 16 bits: number of words of data following
  icomp(istart(j),j) = ibm(j)
  IF (obtzer) ibit(j) = ibit(j) + 128
  IF (obtmis) ibit(j) = ibit(j) + 32
  icomp(istart(j)+1,j) = IOR(ISHFT(ibit(j),16),iword(j)-istart(j)-2)

END DO ! j
!$OMP END PARALLEL DO

! Each row of icomp now has compressed data, in low bits. Need
! to fill the return buffer, icomp64, with both high and low bits
! of the final data.

! Fill first 3 words of compressed data (overall header)
num = SUM(iword(1:iy))+3-iy
itmp(1) = num
itmp(2) = IAND(isc,mask32)
itmp(3) = IOR(ISHFT(ix,16),iy)
IF (iy>0) icomp(iword(iy),iy) = 0

!     Compress to 64 bit words
IF ((num+1)/2 > n) THEN
  icode = 2
  cmessage='COEX: Dimension of ICOMP too small'
  GO TO 9999
END IF


! First put the header information and the first word into icomp64
icomp64(1) = IOR(ISHFT(itmp(1),32),IAND(itmp(2),mask32))
icomp64(2) = IOR(ISHFT(itmp(3),32),IAND(icomp(1,1),mask32))

IF (iy>0) THEN

! Setup addressing
i1 = 2           ! column for first word to pack
j1 = 1           ! row for first word to pack
IF (iword(1) > 3) THEN
  i2 = 3         ! column for second word if on same row
  j2 = 1
ELSE
  i2 = 1         ! column/row for second word if on next row
  j2 = 2
END IF

! Put the rest of the data into icomp64
DO i = 3, (num+1)/2
  icomp64(i) = IOR(ISHFT(icomp(i1,j1),32),IAND(icomp(i2,j2),mask32))

  ! Increment addressing by 2 as we've put 2 words of icomp into icomp64
  i1 = i1 + 2
  i2 = i2 + 2

  ! If we've gone past the end of a row for 1st word, adjust addressing
  IF (i1 >= iword(j1)) THEN
    i1 = i1 - iword(j1) + 1
    j1 = j1 + 1
  END IF

  ! If we've gone past the end of a row for 2nd word, adjust addressing
  IF (i2 >= iword(j2) .AND. j2 /= iy) THEN
    i2 = i2 - iword(j2) + 1
    j2 = j2 + 1
  END IF
END DO

END IF

9999  CONTINUE
IF (lhook) CALL dr_hook('CMPS_ALL',zhook_out,zhook_handle)
RETURN

END SUBROUTINE cmps_all
