! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Integer Function to extract data from lookup array PPXI
!
! Function Interface:

!---------------------------------------------------------------------
!+ Character Function to extract names from lookup array PPXC
!
! Function Interface:
      CHARACTER(LEN=36) FUNCTION EXPPXC(Im_ident,section,item,               &
                                   ErrorStatus ,CMESSAGE)
      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Submodel_Mod
      Use ppxlook_mod, ONLY: ppxc, ppxptr
      USE cppxref_mod, ONLY: ppxref_charlen

      IMPLICIT NONE
!
! Description:
!   Extracts a diagnostic name from ppxref lookup array PPXC.
!
! Method:
!   The required name is identified by the function arguments
!   Im_ident, section, item. The appropriate row in PPXC is found
!   from the 3-d pointer array PPXPTR as PPXPTR(m,s,i).
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Misc
!
! Code description:
!   FORTRAN 77 + common Fortran 90 extensions.
!   Written to UM programming standards version 7.
!
! System component covered:
! System task:               Sub-Models Project
!
! Global Variables:

! Function arguments:
!   Scalar arguments with intent(in):
      INTEGER Im_ident    ! Internal model identifier (absolute)
      INTEGER section     ! STASH section no.
      INTEGER item        ! STASH item no.

!   Scalar arguments with intent(out):
      CHARACTER(LEN=80) CMESSAGE

! Local scalars
      INTEGER row         ! Row no. in PPXC array
      INTEGER I           ! Loop counter
      INTEGER Im_index    ! Internal model index

! Error status:
      INTEGER ErrorStatus !+ve = fatal error

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of Header ---------------------------------------------------

      IF (lhook) CALL dr_hook('EXPPXC',zhook_in,zhook_handle)

      IF (Im_ident <= 0 .OR. section <  0 .OR. item <= 0) THEN
        IF (Im_ident <= 0) WRITE(6,*) 'EXPPXC: INVALID Im_ident'
        IF (section  <  0) WRITE(6,*) 'EXPPXC: INVALID SECTION NO.'
        IF (item     <= 0) WRITE(6,*) 'EXPPXC: INVALID ITEM NO.'
        WRITE(6,'(3(A,I5))')                                   &
     & 'Im_ident ',Im_ident,' section ',section,' item ',item
        ErrorStatus=1
        CMESSAGE='ERROR EXPPXC: INVALID STASH RECORD ID'
      ELSE

! Obtain row no. in PPXC array




        Im_index = INTERNAL_MODEL_INDEX(Im_ident)
        row = PPXPTR(Im_index,section,item)
        IF (row <= 0) THEN
          WRITE(6,'(A,I5)') 'ERROR in EXPPXC: INVALID row VALUE: ',row
          WRITE(6,'(A,3I5)') 'Model,Sec,Item: ',Im_ident,section,item
        END IF


! Obtain required name
        IF (row >  0) THEN
          DO I = 1,PPXREF_CHARLEN
            EXPPXC(I:I) = PPXC(row,I)
          END DO
        ELSE
! Invalid record: return an invalid value.
          exppxc = 'ERROR: invalid STASH record'
        END IF
      END IF
      IF (lhook) CALL dr_hook('EXPPXC',zhook_out,zhook_handle)
      RETURN
      END FUNCTION EXPPXC
