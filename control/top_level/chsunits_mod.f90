! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ 

!  Description:  Defines the number of i/o units

! copied from include file chsunits.h

! Code Owner: See Unified Model Code Owner's HTML page
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90

MODULE chsunits_mod

  IMPLICIT NONE

   ! These values must be consistent with OUTFILE_S, OUTFILE_L
   ! and OUTFILE_E in file VERSION.
  INTEGER,PARAMETER::nunits=179   ! No. of I/O units
  ! length of most unit no arrays
  INTEGER,PARAMETER::nunits_LEN=nunits-19

END MODULE chsunits_mod
