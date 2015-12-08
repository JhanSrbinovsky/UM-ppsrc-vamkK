! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Performs Bi-linear horizitontal interpolation
! Subroutine Interface:
MODULE h_int_bl_halo_mod

IMPLICIT NONE

CONTAINS

SUBROUTINE h_int_bl_halo(rows_in,row_length_in,len_field_out           &
,                        index_b_l,index_b_r,data_in                   &
,                        weight_b_l,weight_b_r,weight_t_l,weight_t_r   &
,                        data_out, source_halo_x, source_halo_y)

!    Purpose:
!
!    Documentation:
!              The interpolation formulae are described in
!              unified model on-line documentation paper S1.
!
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

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER, INTENT(IN)  ::  rows_in              !No of P rows on source grid
INTEGER, INTENT(IN)  ::  row_length_in        !No of pts per row on source grid
INTEGER, INTENT(IN)  ::  len_field_out        !No of points on target grid
INTEGER, INTENT(IN)  ::  source_halo_x, source_halo_y ! Source haloes

!   Array  arguments with intent(in):
INTEGER, INTENT(IN)  ::  index_b_l(len_field_out,2)
                               !Index of bottom lefthand corner
                               !  of source gridbox
INTEGER, INTENT(IN)  ::  index_b_r(len_field_out,2)
                               !Index of bottom righthand corner
                               !  of source gridbox
REAL, INTENT(IN)  ::     data_in(1-source_halo_x:row_length_in+source_halo_x, &
                                 1-source_halo_y:rows_in+source_halo_y)
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
!INTEGER   ::     bl_x,bl_y,br_x,br_y


! End of header

!     1. Carry out horizontal interpolation using equation (2.1)

!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                        &
!$OMP&         SHARED(len_field_out, data_out, weight_b_l, weight_b_r,  &
!$OMP&                weight_t_l, weight_t_r, index_b_l, index_b_r,     &
!$OMP&                data_in)                                          &
!$OMP&            PRIVATE(i)
DO i=1,len_field_out
  data_out(i)=weight_b_l(i)*data_in(index_b_l(i,1),index_b_l(i,2))     &
             +weight_b_r(i)*data_in(index_b_r(i,1),index_b_r(i,2))     &
             +weight_t_l(i)*data_in(index_b_l(i,1),index_b_l(i,2)+1)   &
             +weight_t_r(i)*data_in(index_b_r(i,1),index_b_r(i,2)+1)


END DO
!$OMP END PARALLEL DO

RETURN
END SUBROUTINE h_int_bl_halo
END MODULE h_int_bl_halo_mod
