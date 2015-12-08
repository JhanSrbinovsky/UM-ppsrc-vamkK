! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    Subroutine FIND_MAX_FIELD_SIZE ---------------------------------
!  
!    Purpose:  Reads in and searches LOOKUP header for maximum field
!              size
!  
!  
!    Logical component number: E5
!  
!    External Documentation: None
!  
!     
!    Arguments:--------------------------------------------------------
!
!    Code Owner: See Unified Model Code Owners HTML page
!    This file belongs in section: Small execs

      SUBROUTINE FIND_MAX_FIELD_SIZE(                                   &




     &  NFTIN,LEN1_LOOKUP,LEN2_LOOKUP,FIXHD,MAX_FIELD_SIZE              &

     &  )

      USE IO
      USE io_configuration_mod, ONLY : io_field_padding
      USE ereport_mod, ONLY : ereport

      USE lookup_addresses
      IMPLICIT NONE

      INTEGER                                                           &
     & NFTIN                                                            &
                      !IN Unit number of file
     &,LEN1_LOOKUP                                                      &
                      !IN 1st dim of LOOKUP array
     &,LEN2_LOOKUP                                                      &
                      !IN 2nd dim of LOOKUP array
     &,FIXHD(*)                                                         &
                      !IN Fixed length header
     &,MAX_FIELD_SIZE !OUT Maximum size of field held on file



! -------------------------------------------------------------
! Workspace usage:---------------------------------------------

      INTEGER                                                           &
     & LOOKUP(LEN1_LOOKUP,LEN2_LOOKUP) ! Lookup header

! -------------------------------------------------------------
! Local variables:---------------------------------------------

      INTEGER                                                           &
     & LEN_IO                                                           &
                 ! No of words transferred by BUFFIN
     &,K                                                                &
                 ! Loop index

     &,ICODE     !Return code from setpos
      REAL                                                              &
     & A         ! BUFFIN error code





!*-------------------------------------------------------------

!  Internal structure: none

! Move to start of Look Up Table
      CALL SETPOS(NFTIN,FIXHD(150)-1,ICODE)

! Read in fields from LOOKUP table
      CALL BUFFIN(NFTIN,LOOKUP(:,:),FIXHD(151)*FIXHD(152),LEN_IO,A)

! Check for I/O errors
      IF(A /= -1.0.OR.LEN_IO /= FIXHD(151)*FIXHD(152))THEN
! DEPENDS ON: ioerror
        CALL IOERROR('buffer in of lookup table',A,LEN_IO,              &
     &               FIXHD(151)*FIXHD(152))

        CALL EREPORT('FIND_MAX_FIELD_SIZE', 1000,                       &
     &   'buffer in of lookup table wrong size')

      ENDIF

! Find maximum field size
      MAX_FIELD_SIZE=0
      DO K=1,LEN2_LOOKUP
        MAX_FIELD_SIZE=MAX0(MAX_FIELD_SIZE,LOOKUP(15,K))
        max_field_size=max(max_field_size,lookup(30,k))
      ENDDO

      RETURN
      END SUBROUTINE FIND_MAX_FIELD_SIZE
