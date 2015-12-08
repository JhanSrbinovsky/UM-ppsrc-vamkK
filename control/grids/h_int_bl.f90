! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE h_int_bl_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE h_int_bl(rows_in,row_length_in,len_field_out           &
,                   index_b_l,index_b_r,data_in                   &
,                   weight_b_l,weight_b_r,weight_t_l,weight_t_r   &
,                   data_out)

!    System component: S121
!
!    System task: S1
!
!    Purpose:
!
!    Documentation:
!              The interpolation formulae are described in
!              unified model on-line documentation paper S1.
!

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! Description:
!   Carries out bi-linear horizontal interpolation using coefficients
!   and gather indices calculated in subroutine H_INT_CO

! Method:
!   See UMDP S1 for full desciption

! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Grids

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v8.2 programming standards.

! System component covered: S121
! System Task:              S1


! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN)  ::  rows_in              !No of P rows on source grid
INTEGER, INTENT(IN)  ::  row_length_in        !No of pts per row on source grid
INTEGER, INTENT(IN)  ::  len_field_out        !No of points on target grid

!   Array  arguments with intent(in):
INTEGER, INTENT(IN)  ::  index_b_l(len_field_out)
                               !Index of bottom lefthand corner
                               !  of source gridbox
INTEGER, INTENT(IN)  ::  index_b_r(len_field_out)
                               !Index of bottom righthand corner
                               !  of source gridbox
REAL, INTENT(IN)  ::     data_in(rows_in*row_length_in)
                                !Data before interpolation
REAL, INTENT(IN)  ::     weight_b_l(len_field_out)
                               !Weight applied to value at bottom
                               !lefthand corner of source gridbox
REAL, INTENT(IN)  ::     weight_b_r(len_field_out)
                               !Weight applied to value at bottom
                               !righthand corner of source gridbox
REAL, INTENT(IN)  ::     weight_t_l(len_field_out)
                               !Weight applied to value at top
                               !lefthand corner of source gridbox
REAL, INTENT(IN)  ::     weight_t_r(len_field_out)
                               !Weight applied to value at top
                               !righthand corner of source gridbox

!   Array  arguments with intent(out):
REAL, INTENT(OUT)  ::     data_out(len_field_out) !Data after interpolation

! Local scalars:
INTEGER   ::     i

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

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
! Function & Subroutine calls:
!     External None

! End of header

IF (lhook) CALL dr_hook('H_INT_BL',zhook_in,zhook_handle)

!     1. Carry out horizontal interpolation using equation (2.1)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                        &
!$OMP&         SHARED(len_field_out, data_out, weight_b_l, weight_b_r,  &
!$OMP&                weight_t_l, weight_t_r, index_b_l, index_b_r,     &
!$OMP&                data_in, row_length_in)                           &
!$OMP&            PRIVATE(i)
DO i=1,len_field_out

  data_out(i)=weight_b_l(i)*data_in(index_b_l(i))                   &
             +weight_b_r(i)*data_in(index_b_r(i))                   &
             +weight_t_l(i)*data_in(index_b_l(i)+row_length_in)     &
             +weight_t_r(i)*data_in(index_b_r(i)+row_length_in)

END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook('H_INT_BL',zhook_out,zhook_handle)
RETURN
END SUBROUTINE h_int_bl

END MODULE h_int_bl_mod
