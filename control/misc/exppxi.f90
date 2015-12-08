! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Integer Function to extract data from lookup array PPXI
!
! Function Interface:
      INTEGER FUNCTION EXPPXI(Im_ident,section,item,element,            &
                              ErrorStatus ,CMESSAGE)

      USE yomhook, ONLY: lhook, dr_hook
      USE parkind1, ONLY: jprb, jpim
      USE Submodel_Mod
      USE ppxlook_mod, ONLY: ppxi,ppxptr

      IMPLICIT NONE
!
! Description:
!   Extracts an individual data value from ppxref lookup array PPXI.
!
! Method:
!   The required data element is identified by the function arguments
!   Im_ident, section, item, element. The appropriate row in PPXI is
!   found from the 3-d pointer array PPXPTR as PPXPTR(m,s,i). The
!   address of the required element in PPXI is then given by
!   (row, element).
!
!
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
! History:
! Version  Date  Comment
!  3.4   18/5/94 Add PP missing data indicator. J F Thomson
!  5.1    6/3/00 Convert to Free/Fixed format. P Selwood
!*L------------------COMDECK C_MDI-------------------------------------
      ! PP missing data indicator (-1.0E+30)
      Real, Parameter    :: RMDI_PP  = -1.0E+30

      ! Old real missing data indicator (-32768.0)
      Real, Parameter    :: RMDI_OLD = -32768.0

      ! New real missing data indicator (-2**30)
      Real, Parameter    :: RMDI     = -32768.0*32768.0

      ! Integer missing data indicator
      Integer, Parameter :: IMDI     = -32768
!*----------------------------------------------------------------------

! Function arguments:
!   Scalar arguments with intent(in):
      INTEGER Im_ident    ! Internal model identifier (absolute)
      INTEGER section     ! STASH section no.
      INTEGER item        ! STASH item no.
      INTEGER element     ! Position of required value in PPXI row

!   Scalar arguments with intent(out):
      CHARACTER(LEN=80) CMESSAGE

! Error status:
      INTEGER ErrorStatus !+ve = fatal error

! Local scalars
      INTEGER row         ! Row no. in PPXI array
      INTEGER Im_index    ! Internal model index

      INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
      INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
      REAL(KIND=jprb)               :: zhook_handle

!- End of Header ---------------------------------------------------

      IF (lhook) CALL dr_hook('EXPPXI',zhook_in,zhook_handle)
      ErrorStatus = 0
      IF (Im_ident <= 0 .OR. section <  0 .OR. item <= 0) THEN
        IF (Im_ident <= 0) WRITE(6,*) 'EXPPXI: INVALID Im_ident'
        IF (section  <  0) WRITE(6,*) 'EXPPXI: INVALID SECTION NO.'
        IF (item     <= 0) WRITE(6,*) 'EXPPXI: INVALID ITEM NO.'
        WRITE(6,'(3(A,I5))')                                            &
     & 'Im_ident ',Im_ident,' section ',section,' item ',item
        ErrorStatus=1
        CMESSAGE='ERROR EXPPXI: INVALID STASH RECORD ID'
      ELSE

! Obtain row no. in PPXI array
        Im_index = INTERNAL_MODEL_INDEX(Im_ident)
        row = PPXPTR(Im_index,section,item)
        IF (row <= 0) THEN
          WRITE(6,'(A,I5)') 'ERROR EXPPXI: INVALID row VALUE: ',row
          WRITE(6,'(A,3I5)') 'Im_ident,Sec,Item: ',Im_ident,section,item
        END IF

! Obtain required data value
        IF (row >  0) THEN
          EXPPXI = PPXI(row,element)
        ELSE
! Invalid record: return an invalid value.
          exppxi = IMDI
        END IF
      END IF
      IF (lhook) CALL dr_hook('EXPPXI',zhook_out,zhook_handle)
      RETURN
      END FUNCTION EXPPXI
