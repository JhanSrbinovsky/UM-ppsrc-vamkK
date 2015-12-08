! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Data module defineing data held by trans namelist

Module Rcf_trans_Mod

! Description:
!   Defines arrays to hold data from the trans namelist
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! History:
! Version   Date     Comment
! -------   ----     -------
!   5.1   16/03/00   Original code.  P.Selwood
!   5.3   03/10/01   Addition of Implicit None. R.Sharp
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


Implicit None

Integer, Allocatable      :: ItemC_array(:)   ! item code array
Integer, Allocatable      :: SctnC_array(:)   ! section code array
Integer, Allocatable      :: Lev1_array(:)    ! bottom level
Integer, Allocatable      :: Lev2_array(:)    ! top level
Integer, Allocatable      :: Col1_array(:)    ! 1st column
Integer, Allocatable      :: Col2_array(:)    ! 2nd column
Integer, Allocatable      :: Row1_array(:)    ! 1st row
Integer, Allocatable      :: Row2_array(:)    ! 2nd row

Integer                   :: Num_Trans     ! number of transplant items

End Module Rcf_trans_Mod
