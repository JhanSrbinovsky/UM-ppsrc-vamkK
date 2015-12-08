! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Controls horizontal interpolation

Module Rcf_H_Int_Nearest_Mod

!  Subroutine Rcf_H_Int_Nearest_Mod - controls horizontal interpolation.
!
! Description:
!   Add neareat method when h_int_method = 3. Changgui Wang (24/9/07)
!
! Method:
!   Chooses method based on h_int_method variable.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

Contains

SUBROUTINE Rcf_H_Int_Nearest(ROWS_IN,ROW_LENGTH_IN,LEN_FIELD_OUT,          &
                             BL_INDEX_B_L,BL_INDEX_B_R,DATA_IN,            &
                             WEIGHT_B_L,WEIGHT_B_R, WEIGHT_T_L,WEIGHT_T_R, &
                             DATA_OUT)

IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER      ROWS_IN          !No of rows on source grid
INTEGER      ROW_LENGTH_IN    !No of pts per row on source grid
INTEGER      LEN_FIELD_OUT    !No of points on target grid

!   Array  arguments with intent(in):
INTEGER, POINTER ::      BL_INDEX_B_L(:)
                              !Gather index for bottom l.h.c of
                              !source grid box. 1=P-pts; 2=UV-pts
INTEGER, POINTER ::      BL_INDEX_B_R(:)
                              !Gather index for bottom r.h.c of
                              !source grid box. 1=P-pts; 2=UV-pts
REAL          ::         DATA_IN(:)
                              !Data before interpolation
REAL, POINTER ::         WEIGHT_B_L(:)  !\ Weights used in
REAL, POINTER ::         WEIGHT_B_R(:)  ! \bilinear horizontal
REAL, POINTER ::         WEIGHT_T_L(:)  ! /interpolation
REAL, POINTER ::         WEIGHT_T_R(:)  !/ 1=P-pts; 2=U-pts
                                        !  3=V-pts; 4=zonal
                                        !             means   


!   Array  arguments with intent(out):
REAL         DATA_OUT(*)
                              !Data after interpolation
Integer                      :: ErrorStatus
Character (Len=*), Parameter :: RoutineName = 'Rcf_H_Int_Nearest'
Character (Len=80)           :: Cmessage

! Local
real        ::nest_mask(4)
integer     ::i,nearest(1)

!- End of header

!      1. Carry out horizontal interpolation using nearest neighbour

DO I=1,LEN_FIELD_OUT
     nest_mask(1)          = WEIGHT_B_L(I)
     nest_mask(2)          = WEIGHT_B_R(I)
     nest_mask(3)          = WEIGHT_T_L(I)
     nest_mask(4)          = WEIGHT_T_R(I)
     nearest               = MAXLOC( nest_mask )
     nest_mask             = 0.0
     nest_mask(nearest(1)) = 1.0

     DATA_OUT(I) = nest_mask(1) * DATA_IN(BL_INDEX_B_L(I)) +  &
                   nest_mask(2) * DATA_IN(BL_INDEX_B_R(I)) +  &
                   nest_mask(3) * DATA_IN(BL_INDEX_B_L(I)  +  &
                   ROW_LENGTH_IN)                          +  &
                   nest_mask(4) * DATA_IN(BL_INDEX_B_R(I)  +  &
                   ROW_LENGTH_IN)
END DO

RETURN
END Subroutine Rcf_H_Int_Nearest
End Module Rcf_H_Int_Nearest_Mod
