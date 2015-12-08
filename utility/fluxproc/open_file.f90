! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: Open_file
!
! Purpose: Flux processing routine.
!          Opens an external file
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine Open_file(Unit,Form,Stat,ICode)

      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn, &
                           cerr, cstd, csub
      implicit none

!  declaration of argument list
      INTEGER      Unit      !IN      Unit Number on which to open file
      CHARACTER(LEN=11) Form      !IN      Format with which to open file
      CHARACTER(LEN=7)  Stat      !IN      Status with which to open file
      INTEGER      ICode     !IN/OUT  Return Code from file open

! declaration of glbbals

! declaration of local scalars
      CHARACTER(LEN=1)  CUn1D       ! Unit Number as Characters(1 digit)
      CHARACTER(LEN=2)  CUn2D       ! Unit Number as Characters(2 digits)
      CHARACTER(LEN=5)  FName       ! Filename for open statement
!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'Open_file'  ! subroutine name for error messages


! 2. Construct file name

      IF (Unit <= 9) THEN
        CUn1D='0'
        WRITE(CUn1D,'(I1)') Unit
        FName='ftn'//CUn1D
      ELSE
        CUn2D='00'
        WRITE(CUn2D,'(I2)') Unit
        FName='ftn'//CUn2D
      ENDIF

! 3. Open file

      OPEN(UNIT=Unit,FILE=FName,FORM=Form,STATUS=Stat,IOSTAT=ICode)

9999  continue
      return
      END SUBROUTINE Open_file
!----------------------------------------------------------------------
