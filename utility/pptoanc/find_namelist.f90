! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

      FUNCTION FIND_NAMELIST(iunit, namelist_name)

      implicit none
!
! Description:
!  This routine searches the input stream given by 'iunit'
!  to find the NAMELIST given by 'namelist_name'.  The
!  input file is then correctly positioned to let F90
!  library routines read the namelist.
!
!  The namelist is assumed to be contained within the
!  first 24 characters of the record
!
!  Return values:
!
!  -1  error - could not find the namelist - file is now
!      at end-of-file.
!   0  namelist ready to be processed.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
!
! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.
!
!
! Declarations:
! 1.0 Subroutine arguments
!   1.1 Scalar arguments with intent(in):
      integer iunit

!   1.2 Scalar arguments with intent(out):
      CHARACTER(LEN=*) namelist_name

      integer find_namelist

! 2.0 Local scalars:
      integer i, j

      CHARACTER(LEN=24) chvar

!- End of header

1000  continue
      read(iunit, '(a)', end=9000, err=9000) chvar

!  Check for leading '&' for namelist
      i=index(chvar, '&')
!  Found '&' - check for the name we want
      if(i /= 0) then
        j=index(chvar, namelist_name)
!  Not the name we want - print skipped message
        if(j == 0) then
          if(index('endEndeNdenDENdEnDeNDEND',                          &
     &     chvar(i+1:i+3)) == 0) then
            write(0,*)'- Skipped record named: ',                       &
     &       chvar(i+1:),' On Unit:',iunit,                             &
     &       ' - f90 version'
          endif
          goto 1000
        endif
        goto 1100
      endif
      goto 1000

!  Found the namelist we want - backspace to position correctly
1100  continue
      backspace iunit
      FIND_NAMELIST=0
      return

!  Cannot find the namelist we want
9000  continue
      FIND_NAMELIST=1
      return
      END FUNCTION FIND_NAMELIST
