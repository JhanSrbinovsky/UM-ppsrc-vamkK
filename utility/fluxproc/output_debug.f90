! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Programming standard: Unified Model Documentation Paper No 3
!                       Version No 1 15/1/90
!----------------------------------------------------------------------
! contains routines: output_debug
!
! Purpose: Flux processing routine.
!          Writes out field values at selected grid points
!----------------------------------------------------------------------
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Small execs
      subroutine output_debug(CMessage, nrows, ncols, Field)

      USE cmess_mod, ONLY: outunitdbg, unwarn, unerr, unstd, cwarn,     &
                           cerr, cstd, csub
      implicit none

! argument list: all intent IN
      CHARACTER(LEN=*) CMessage    ! message to output
      integer    nrows          ! input dimension; number of rows
      integer    ncols          ! input dimension; number of columns
      real       Field(ncols,nrows)   ! field to output

! globals
!----------------------------------------------------------------------
! comdeck: CDEBUG
! Purpose: declares and stores information needed to output debugging
!          diagnostics at user selected points as fields are processed.
!----------------------------------------------------------------------
! common block:
      COMMON / CDebug /                                                 &
     &  NoDbgPts,IColDbg,JRowDbg,l_winds_dbg,l_heat_dbg,l_moisture_dbg, &
     &  l_sea_ice_dbg, l_references_dbg, l_pressure_dbg,l_windspd_dbg

      common / CCDebug / CValues

! declarations:
      integer MaxNoDbgPts   ! parameter: max. number of points to output
      parameter ( MaxNoDbgPts = 20)

! points to output
      integer NoDbgPts      ! actual number of points to output
      integer IColDbg(MaxNoDbgPts)   ! column of each point
      integer JRowDbg(MaxNoDbgPts)   ! row    of each point

! character array for output
      CHARACTER(LEN=11) CValues(MaxNoDbgPts) ! values to write out

! debug logical for each output file
      LOGICAL :: l_winds_dbg
      LOGICAL :: l_heat_dbg
      LOGICAL :: l_moisture_dbg
      LOGICAL :: l_sea_ice_dbg
      LOGICAL :: l_references_dbg
      LOGICAL :: l_pressure_dbg
      LOGICAL :: l_windspd_dbg

! CDEBUG end

! local scalars
      integer IDbg   ! loop index over selected points to debug

!----------------------------------------------------------------------
! 0. Preliminaries
      CSub = 'output_debug'  ! subroutine name for error messages


! 1. write out field descriptor
      if ( CMessage  /=  '         ') then
        write(OutUnitDbg, *) CMessage
      end if

! 2. convert values to output
      do IDbg = 1, NoDbgPts

        if ( IColDbg(IDbg)  <=  ncols .and.                             &
     &       JRowDbg(IDbg)  <=  nrows       ) then
           write(CValues(IDbg), '(G11.4)' )                             &
     &                          Field( IColDbg(IDbg), JRowDbg(IDbg) )
         else
           CValues(IDbg) = ' OOB '
         end if

      end do ! IDbg

! 3. output values
        write(OutUnitDbg, '(11A11)' ) (CValues(IDbg), IDbg=1,NoDbgPts)

      return
      END SUBROUTINE output_debug
!----------------------------------------------------------------------
