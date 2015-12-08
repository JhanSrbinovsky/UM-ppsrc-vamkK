! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: See Unified Model Code owners HTML page
! This file belongs in section: C96

!
! An F90 wrapper to coex packing
!

MODULE ios_stash_wgdos

  USE UM_Types
  USE yomhook, ONLY: lhook, dr_hook
  USE parkind1, ONLY: jprb, jpim

  IMPLICIT NONE

  INTEGER, PARAMETER :: coexHeaderWords32=3


! params for dr_hook
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
  INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CONTAINS

  SUBROUTINE ios_stash_pack_wgdos(                                             &
      field,                                                                   &
      compressedField,                                                         &
      num,                                                                     &
      accuracy,                                                                &
      missingDataIndicator)

    USE ios
    USE IOS_Common
    USE UM_Types

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: accuracy ! accuracy is accuracy as a power of two
    REAL,    INTENT(IN)  :: field(:,:)          ! is the field to compress
    REAL,    INTENT(IN)  :: MissingDataIndicator! parameter to the algorithm
    INTEGER, INTENT(OUT) :: num
    REAL(KIND=real32),                                                         &
        INTENT(OUT)      :: compressedField(:)  ! is the output buffer

    INTEGER              :: ix  ! ix,iy are the field dimensions
    INTEGER              :: iy
    INTEGER              :: icode               ! is the error code
    INTEGER              :: lword               ! a bug
    INTEGER              :: i
    CHARACTER(LEN=80)       :: message             ! is the error message
    icode=0

! Note that input buffer is 64 bit, but output buffer is 32 bit, leading to
! occasional odd factors of two in specifying buffer sizes.

    ix=SIZE(field,1)
    iy=SIZE(field,2)

    IF ( SIZE(compressedField) < 2*(ix+coexHeaderWords32)*iy+                  &
        coexHeaderWords32 )  THEN
      WRITE(message,'(A,I10,A,I10,A)')                                         &
          'size of compressedField array (',SIZE(compressedField),             &
          ') smaller than potential need  (',2*(ix+coexHeaderWords32)*iy,')'
      CALL IOS_Ereport("IO SERVER WGDOS MODULE",99,message)
    END IF

! DEPENDS ON: coex
    CALL coex(                                                                 &
        field,    &  ! data
        SIZE(field),    &  ! data size
        compressedField,&! data out
        SIZE(compressedField),                                                 &
                      ! but we will abort if coex claims to expand the
                      ! field
        ix,       &  ! x size
        iy,       &  ! y size
        num,      &  ! count of 32 bit words in the compressed buffer (output)
        accuracy, &  ! accrcy
        .TRUE.,   &  ! Packing (ie not unpacking)
        MissingDataIndicator,                                                  &
        lword,                                                                 &
        icode,                                                                 &
        message)

    ! Compression factor - note factor of 2 as num is in 32 bit words
    IF (IOS_Verbosity>=IOS_PrStatus_Diag) THEN
      IF (ix*iy>0) THEN !Avoid div0
      WRITE(6,'(A,F5.2,A,I4)')                                                 &
          'IOS: Info: Stash Server: WGDOS Compression is ',                    &
          (100.0*num)/(2*ix*iy),'%, for accuracy of ',accuracy
      ELSE
        WRITE(6,'(A)') 'IOS: Info: Stash Server: WGDOS Compression: No Data'
      END IF
    END IF

    IF (num >= SIZE(compressedField)) THEN
      WRITE(6,'(A)')                                                           &
          'A buffer overrun may have occurred'
      WRITE(6,'(A)')
      WRITE(6,'(A)')       'Arguments after call to coex:'
      WRITE(6,'(A,I12)')   'data size: ',ix*iy
      WRITE(6,'(A,I12)')   'compressed field size: ',SIZE(compressedField)
      WRITE(6,'(A,I12)')   'x:',ix
      WRITE(6,'(A,I12)')   'y:',iy
      WRITE(6,'(A,I12)')   'num (32 bit words after compression:',num
      WRITE(6,'(A,I12)')   'accuracy:',accuracy
      WRITE(6,'(A,F32.16)')'MDI:',MissingDataIndicator
      WRITE(6,'(A,I12)')   'lword:',lword
      WRITE(6,'(A,I12)')   'icode:',icode

      DO i=1,iy
        WRITE(6,*)field(:,i)
      END DO
      CALL IOS_Ereport("IO SERVER WGDOS MODULE",99,                            &
          'Potential buffer overrun due to -ve compression')
    END IF

    RETURN
  END SUBROUTINE  ios_stash_pack_wgdos

END MODULE ios_stash_wgdos
