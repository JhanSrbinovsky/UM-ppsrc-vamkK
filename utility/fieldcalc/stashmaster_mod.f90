! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module to interface with the STASHmaster information.

MODULE stashmaster_mod

! Description:
!   Defines the following subroutines:
!    * init_stm()
!    * stm_get_grid_type(model,section,item)
!
!   It can be easily expanded to cater for other information required from the
!   STASHmaster.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: FieldCalc
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.

USE ereport_mod, ONLY: ereport
USE cppxref_mod, ONLY: ppx_grid_type
USE ppxlook_mod, ONLY: ppxrecs

IMPLICIT NONE

PUBLIC :: init_stm, stm_get_grid_type

CONTAINS
SUBROUTINE init_stm()
! This has to be unique.  Lets give it a nice high number.
INTEGER, PARAMETER :: nft1 = 40
INTEGER, PARAMETER :: nft2 = 2
CHARACTER(len=*), PARAMETER :: progname = "init_stm"
INTEGER :: rownumber

CHARACTER(len=80) :: cmessage
INTEGER :: icode

ppxrecs=1
rownumber=0
cmessage = ' '
icode=0
! DEPENDS ON: hdppxrf
CALL hdppxrf(nft1,'STASHmaster_A',icode,CMESSAGE)
IF (icode >  0) THEN
  WRITE(6,*) 'Error reading STASHmaster_A'
  WRITE(6,*) cmessage

  CALL ereport(progname, icode, 'Error reading STASHmaster_A')
END IF

! DEPENDS ON: getppx
CALL getppx(nft1,nft2,'STASHmaster_A',rownumber,                               &
  icode,cmessage)
IF (icode /= 0) THEN

  CALL ereport(progname, icode, cmessage)
END IF

! Comment this out - this requires a large modification to how to input
! to fieldcalc.
!!User STASHmaster
!! DEPENDS ON: hdppxrf
!CALL HDPPXRF(0,' ',ICODE, CMESSAGE)
!IF(ICODE /= 0)THEN
!  WRITE(6,*) CMESSAGE
!
!  CALL EREPORT(progname, ICODE, CMESSAGE)
!ENDIF
!
!! DEPENDS ON: getppx
!CALL GETPPX(0,NFT2,' ',RowNumber,                                             &
!            ICODE,CMESSAGE)
!IF(ICODE /= 0)THEN
!
!  CALL EREPORT(progname, ICODE, CMESSAGE)
!ENDIF

RETURN
END SUBROUTINE init_stm

FUNCTION stm_get_grid_type(model,section,item)

INTEGER :: stm_get_grid_type
INTEGER :: model
INTEGER :: section
INTEGER :: item

CHARACTER(len=80) :: cmessage
INTEGER :: icode
INTEGER :: l_stm_get_grid_type
CHARACTER(len=*), PARAMETER :: progname = "stm_get_grid_type"
! Function to extract STASHmaster value
INTEGER :: exppxi
EXTERNAL exppxi

cmessage = " "
icode    = 0
l_stm_get_grid_type = 0
! DEPENDS ON: exppxi
l_stm_get_grid_type=exppxi(model,section,item,                                 &
                           ppx_grid_type,                                      &
                           icode,cmessage)
IF (icode /= 0) THEN

  CALL ereport(progname, icode, cmessage)
END IF

stm_get_grid_type = l_stm_get_grid_type
RETURN
END FUNCTION stm_get_grid_type

END MODULE stashmaster_mod
