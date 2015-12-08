! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!+ Controls horizontal interpolation

Module Rcf_H_Int_Ctl_Mod

!  Subroutine rcf_h_int_ctl - controls horizontal interpolation.
!
! Description:
!   Wrapper routine choosing correct horizontal interpolation method
!
! Method:
!   Chooses method based on h_int_method variable.
!
! Code Owner: See Unified Model Code Owners HTML page
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

Contains

SUBROUTINE Rcf_H_INT_CTL(LEN_FIELD_OUT, ROW_LENGTH_IN, ROW_LENGTH_OUT, &
                        ROWS_IN, ROWS_OUT,GLOBAL,AW_INDEX_TARG_LHS,    &
                        AW_INDEX_TARG_TOP, BL_INDEX_B_L,BL_INDEX_B_R,  &
                        AW_COLAT_T,AW_LONG_L,DATA_IN,WEIGHT_T_R,       &
                        WEIGHT_B_R,WEIGHT_T_L,WEIGHT_B_L, DATA_OUT)

Use Rcf_Interp_Weights_Mod, Only : &
    h_int_method,                  &
    bilinear,                      &
    area_weighted,                 &
    nearest_neighbour

Use Rcf_H_Int_Nearest_Mod, Only : &
    Rcf_H_Int_Nearest

Use Ereport_mod, Only : &
    Ereport

USE h_int_bl_mod, ONLY: h_int_bl
IMPLICIT NONE

! Subroutine arguments
!   Scalar arguments with intent(in):
INTEGER      LEN_FIELD_OUT    !No of points on target grid
INTEGER      ROWS_IN          !No of rows on source grid
INTEGER      ROWS_OUT         !No of rows on target grid
INTEGER      ROW_LENGTH_IN    !No of pts per row on source grid
INTEGER      ROW_LENGTH_OUT   !No of pts per row on target grid
LOGICAL      GLOBAL           !True if global area required

!   Array  arguments with intent(in):
INTEGER, POINTER ::      AW_INDEX_TARG_LHS(:)
                              !Index of source box overlapping
                              !lhs of target grid-box
INTEGER, POINTER ::      AW_INDEX_TARG_TOP(:)
                              !Index of source box overlapping
                              !top of target grid-box
INTEGER, POINTER ::      BL_INDEX_B_L(:)
                              !Gather index for bottom l.h.c of
                              !source grid box. 1=P-pts; 2=UV-pts
INTEGER, POINTER ::      BL_INDEX_B_R(:)
                              !Gather index for bottom r.h.c of
                              !source grid box. 1=P-pts; 2=UV-pts
REAL, POINTER ::         AW_COLAT_T(:)
                              !Colatitude of top of target grd-box
                              ! (in units of DELTA_LAT_SRCE)
REAL, POINTER ::         AW_LONG_L(:)
                              !Left longitude of target grid-box
                              ! (in units of DELTA_LONG_SRCE)
REAL          ::         DATA_IN(:)
                              !Data before interpolation
REAL, POINTER ::         WEIGHT_T_R(:) !\ Weights used in
REAL, POINTER ::         WEIGHT_B_R(:) ! \bilinear horizontal
REAL, POINTER ::         WEIGHT_T_L(:) ! /interpolation
REAL, POINTER ::         WEIGHT_B_L(:) !/ 1=P-pts; 2=U-pts
                                       !  3=V-pts; 4=zonal
                                       !             means
!   Array  arguments with intent(out):
REAL         DATA_OUT(*)
                              !Data after interpolation
Integer                      :: ErrorStatus
Character (Len=*), Parameter :: RoutineName = 'H_INT_CTL'
Character (Len=80)           :: Cmessage

!- End of header

Select Case ( h_int_method )
  Case ( bilinear )
  
  ! Bi-linear interpolation requested

      CALL H_INT_BL(ROWS_IN,ROW_LENGTH_IN,LEN_FIELD_OUT,          &
                    BL_INDEX_B_L,BL_INDEX_B_R,DATA_IN,            &
                    WEIGHT_B_L,WEIGHT_B_R, WEIGHT_T_L,WEIGHT_T_R, &
                    DATA_OUT)

  Case ( nearest_neighbour )

  ! Nearest neighbour 
      Call Rcf_H_Int_Nearest(ROWS_IN,ROW_LENGTH_IN,LEN_FIELD_OUT,          &
                             BL_INDEX_B_L,BL_INDEX_B_R,DATA_IN,            &
                             WEIGHT_B_L,WEIGHT_B_R, WEIGHT_T_L,WEIGHT_T_R, &
                             DATA_OUT)
 
  Case ( area_weighted )

  !  Area weighted interpolation

! DEPENDS ON: h_int_aw
    CALL H_INT_AW(ROWS_IN,ROWS_OUT,                    &
                  ROW_LENGTH_IN,ROW_LENGTH_OUT,GLOBAL, &
                  AW_INDEX_TARG_LHS,AW_INDEX_TARG_TOP, &
                  AW_COLAT_T,AW_LONG_L,DATA_IN,DATA_OUT)




  Case Default
    Cmessage = 'Unsupported interpolation method'
    ErrorStatus = 10
    Call Ereport( RoutineName, ErrorStatus, Cmessage )

End Select

RETURN
END Subroutine Rcf_H_Int_Ctl



End Module Rcf_H_Int_Ctl_Mod
