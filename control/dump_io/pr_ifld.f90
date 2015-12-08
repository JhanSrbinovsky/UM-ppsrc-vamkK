! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!    SUBROUTINE PR_IFLD----------------------------------------
!
!    Purpose: Prints out selected values from integer data
!             using information from associated PP header.
!
!
!    Programming standard:
!             Unified Model Documentation Paper No 3
!
!    Documentation:
!             Unified Model Documentation Paper No F3
!
!             Code Owner: See Unified Model Code Owners HTML page
!             This file belongs in section: Dump I/O

SUBROUTINE pr_ifld(lookup,rlookup,id1,k)


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE lookup_addresses

IMPLICIT NONE

INTEGER                                                           &
 k                                                                &
               !IN Field number ie position in 2nd dim of LOOKUP
,lookup(64,*)                                                     &
               !IN Integer equivalence of PP LOOKUP
,id1(*)        !IN Kth field in integer equiv of data array

REAL                                                              &
 rlookup(64,*) !IN Real equivalence of PP LOOKUP

! -------------------------------------------------------------
!   Local control constants:-----------------------------------
INTEGER                                                           &
 ns_pts                                                           &
               !PARAM No of points down to print
,ew_pts        !PARAM No of points across to print
PARAMETER(ns_pts=6,ew_pts=5)
! -------------------------------------------------------------
! Workspace usage:---------------------------------------------
REAL lon(ew_pts)     ! Longitudes printed out
INTEGER i(ew_pts)    ! Index of values printed out
CHARACTER(LEN=12) dash(ew_pts)  !Stores dashed lines
!*-------------------------------------------------------------
! Local variables:---------------------------------------------
INTEGER                                                           &
 n_rows                                                           &
             ! No of rows in field
,n_cols                                                           &
             ! No of colums in field
,points                                                           &
,row                                                              &
             ! Row number
,r_inc,f_inc                                                      &
             ! No of rows/points between printed lines
,j,l                                                              &
             ! Loop counts
,ew_print                                                         &
             ! No of E-W values printed out
,pos_min                                                          &
             ! Position of Minimum value of field
,pos_max                                                          &
             ! Position of Maximum value of field
,f_min                                                            &
             ! Minimum value of field
,f_max       ! Maximum value of field

REAL                                                              &
 lat         ! Latitude

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!--------------------------------------------------------------

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

! Initialise string used to create table boundaries
IF (lhook) CALL dr_hook('PR_IFLD',zhook_in,zhook_handle)
DO j=1,ew_pts
  dash(j)='------------'
END DO

IF (lookup(lbcode,k) == imdi) THEN
!       IF LBCODE IS MISSING DATA, ASSUME THAT THE FIELD IN DUMP
!       HAS NOT BEEN WRITTEN TO BY STASH.
!       THIS SHOULD ONLY OCCUR TO DIAGNOSTIC PARTS OF THE DUMP BEFORE
!       FIRST WRITE BY STASH TO THAT AREA/HEADER.
  WRITE(6,*) 'MESSAGE FROM PR_IFLD'
  WRITE(6,*) 'LBCODE NOT SET; ASSUME DATA NOT SET. NO PRINT'
  IF (lhook) CALL dr_hook('PR_IFLD',zhook_out,zhook_handle)
  RETURN
END IF

! No of rows and columns in field
n_rows=lookup(lbrow,k)
n_cols=lookup(lbnpt,k)


IF (n_cols /= 0.AND.n_cols /= imdi) THEN

! No of E-W values to be printed
  ew_print=MIN(n_cols,ew_pts)

! Calculate longitudes and addresses of values to be printed from 1st ro
  i(1)=1
  lon(1)=rlookup(bzx,k)+rlookup(bdx,k)
  DO j=1,ew_pts-2
    i(j+1)=i(j)+n_cols/(ew_pts-1)
    lon(j+1)=lon(j)+rlookup(bdx,k)*(n_cols/(ew_pts-1))
  END DO
  i(ew_pts)=n_cols
  lon(ew_pts)=rlookup(bzx,k)+rlookup(bdx,k)*n_cols

! Initialise row and field pointers
  row=1
  lat=rlookup(bzy,k)+rlookup(bdy,k)
  r_inc=n_rows/(ns_pts-1)
  f_inc=r_inc*n_cols

! Print 1st row
  WRITE(6,'(14X,9A12)')(dash(j),j=1,ew_print)
  WRITE(6,'('' FIELD NO'',i4,'':'',9(F10.3,2X))')                 &
        k,(lon(j),j=1,ew_print)
  WRITE(6,'(14X,9A12)')(dash(j),j=1,ew_print)

! Print remaining rows except last
  DO l=1,ns_pts-1
    WRITE(6,'(1x,i3,'':'',f8.3,'':'',3X,9(I9,3X))')row,lat,       &
          (id1(i(j)),j=1,ew_print)
    DO j=1,ew_pts
      i(j)=i(j)+f_inc
    END DO
    row=row+r_inc
    lat=lat+r_inc*rlookup(bdy,k)
  END DO

! Calculate addresses used to print values for last row
  i(1)=1+(n_rows-1)*n_cols
  DO j=1,ew_pts-2
    i(j+1)=i(j)+n_cols/(ew_pts-1)
  END DO
  i(ew_pts)=n_rows*n_cols

! Set row pointers to last row
  lat=rlookup(bzy,k)+rlookup(bdy,k)*n_rows
  row=n_rows

! Print last row
  WRITE(6,'(1x,i3,'':'',f8.3,'':'',3X,9(I9,3X))')row,lat,         &
        (id1(i(j)),j=1,ew_print)
  WRITE(6,'(14X,9A12)')(dash(j),j=1,ew_print)
ELSE

! Print out summary of non standard fields

  ew_print=MIN(ew_pts,lookup(lblrec,k))
  WRITE(6,'(14X,9A12)')(dash(j),j=1,ew_print)
  WRITE(6,'(A,I4,A)') ' FIELD NO',k,                              &
  ':  DATA NOT ON MODEL GRID SO FIRST FEW VALUES PRINTED'
  WRITE(6,'(14X,9A12)')(dash(j),j=1,ew_print)
  WRITE(6,'(1x,3x,'':'',8x,'':'',3X,9(I9,3X))')                   &
        (id1(j),j=1,ew_print)
  WRITE(6,'(14X,9A12)')(dash(j),j=1,ew_print)

END IF


! Print out max and min values of field
f_min=id1(1)
f_max=id1(1)
pos_max=1
pos_min=1
IF (lookup(lbext,k) >  0) THEN
  points=n_rows*n_cols
ELSE
  points=lookup(lblrec,k)
END IF
DO j=1,points
  IF (id1(j) >  f_max) THEN
    f_max=id1(j)
    pos_max=j
  END IF
  IF (id1(j) <  f_min) THEN
    f_min=id1(j)
    pos_min=j
  END IF
END DO

WRITE(6,'(A,I8,A,I8,A,I8,A,I8)')                                &
     ' MINIMUM=',f_min,' POSITION=',pos_min,                    &
     ' MAXIMUM=',f_max,' POSITION=',pos_max   

WRITE(6,'('' '')')

IF (lhook) CALL dr_hook('PR_IFLD',zhook_out,zhook_handle)
RETURN
END SUBROUTINE pr_ifld

