! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!LL  Routines: OA_PACK, OA_UNPACK and OA_LEV_CMP
!LL
!LL
!LL  Logical components covered:
!LL
!LL
!LL  Programming standard : FOAM Doc Paper 3/2/1
!LL
!
!  Code Owner: See Unified Model Code Owners HTML page
!  This file belongs in section: Small execs
      SUBROUTINE OA_PACK(ICODE, CMESSAGE, LL_AC_TIM,                    &
     & NO_ROWS_M, NO_COLS_M, NO_LEVS_M, NO_SEG, NO_CMP_KLEV,            &
     & INDX_CMP, INDX_EXP, INDX_TO_ROWS, NO_CMP, REAL_MDI,              &
     & KLEV, I_FLD_TYP, LL_CYC_M, FLD_EXP, FLD_CMP)
!*
!LL  Purpose:  This subroutine packs one full field of model data
!LL            at a given level  into compressed form.
!
      IMPLICIT NONE
!
!*L  ARGUMENT LIST
!
      INTEGER ICODE      ! return error code
      CHARACTER(LEN=256) CMESSAGE      ! return error message
      LOGICAL LL_AC_TIM  ! T => time output on input and exit
! Dimensions
      INTEGER NO_ROWS_M  ! IN number of rows on model grid
      INTEGER NO_COLS_M  ! IN number of columns on expanded model grid
      INTEGER NO_LEVS_M  ! IN number of levels on model grid
      INTEGER NO_SEG     ! IN total number of sea segments
      INTEGER NO_CMP_KLEV! IN number of sea points on this level
! Compression and expansion indices
      INTEGER INDX_CMP(NO_SEG) ! IN contains position in compressed arra
!                                   of start of each sea segment
      INTEGER INDX_EXP(NO_SEG) ! IN contains position in expanded array
!                                   of start of each sea segment
      INTEGER INDX_TO_ROWS(NO_ROWS_M*NO_LEVS_M) ! IN contains number of
!                          first/next sea segment for each row and level
      INTEGER NO_CMP           ! IN total no of points in a 3D
!                                   compressed array
      REAL REAL_MDI            ! IN missing data indicator
! Level and field type
      INTEGER KLEV             ! IN model level
      INTEGER I_FLD_TYP        ! IN =0 for tracers, =1 for currents
      LOGICAL LL_CYC_M         ! T => FLD_EXP is cyclic in W-E
! Fields
      REAL FLD_EXP(NO_COLS_M*NO_ROWS_M) ! IN  data on expanded grid
      REAL FLD_CMP(NO_CMP_KLEV)         ! OUT data on compressed grid
!*
!L   NO PARAMETERS
!
!L*  NO WORK ARRAYS
!
!    NO EXTERNAL SUBROUTINES CALLED
!
!    LIST OF OTHER VARIABLES
!
      INTEGER ICYC        ! number of columns of expanded field; sec. 1.
      INTEGER ICOUNT      ! index for loop over points in segment
      INTEGER INC_CYC     ! extra columns on cyclic grid; sec. 1
      INTEGER IPT_CMP     ! pointer to location in compressed field
      INTEGER IPT_EXP     ! pointer to location in 3D expanded field
      INTEGER IPT_EXP_CYC ! pointer to location in 2D expanded cyclic fi
      INTEGER IPT_SEG     ! segment index number
      INTEGER ISEG        ! index for loop over segments in row
      INTEGER ISEG_ST     ! first segment on this level
      INTEGER IST_CMP_M1  ! index of 1st point in compressed field at th
!                           level, minus one
      INTEGER IST_EXP_M1  ! index of 1st point in expanded field at this
!                           level (for non-cyclic grid), minus one
      INTEGER J           ! index for loop over rows on level
      INTEGER JPT         ! point index number
      INTEGER LEN_SEG     ! number of grid points in current segment
      INTEGER NO_SEG_ROW  ! number of segments in row
!*
!-----------------------------------------------------------------------
!
!L 1. Prelminaries
!
!L 1.1 Set the number of columns of distinct data
      IF (LL_CYC_M) THEN
        INC_CYC = 2
        ICYC = NO_COLS_M - 2
      ELSE
        INC_CYC = 0
        ICYC = NO_COLS_M
      END IF
!
!L 1.2 Set offsets for compressed and expanded fields at this level
      ISEG_ST = INDX_TO_ROWS(NO_ROWS_M*(KLEV-1) + 1)
      IST_CMP_M1 = INDX_CMP(ISEG_ST) - 1
      IST_EXP_M1 = (KLEV-1)*NO_ROWS_M*ICYC
!
!L 2. Loop over rows (index J)
!
!
!L 2.2 Start loop over rows and define the pointer to the row
!
      DO J = 1, NO_ROWS_M
        JPT = (KLEV - 1)*NO_ROWS_M + J
!
!L 2.3 Calculate the number of sea segments in the row
!
         IF (JPT  ==  NO_LEVS_M*NO_ROWS_M ) THEN
            NO_SEG_ROW = NO_SEG - INDX_TO_ROWS(JPT) + 1
         ELSE
            NO_SEG_ROW = INDX_TO_ROWS(JPT+1) - INDX_TO_ROWS(JPT)
         END IF
!
!L 2.4 Start loop over sea segments and define pointer to segment
         DO ISEG = 1, NO_SEG_ROW
            IPT_SEG = INDX_TO_ROWS(JPT) + ISEG - 1
!
!L 2.5 Calculate the length of the present sea segment
!
            IF (IPT_SEG  <   NO_SEG) THEN
               LEN_SEG = INDX_CMP(IPT_SEG+1) - INDX_CMP(IPT_SEG)
            ELSE
               LEN_SEG = NO_CMP - INDX_CMP(IPT_SEG) + 1
            END IF
!
!L 2.6 Calculate FLD_CMP for all points in the segment
!L     (the last point may be overwritten in 2.7 if I_FLD_TYP = 1)
!
            DO ICOUNT = 1, LEN_SEG
              IPT_EXP = INDX_EXP(IPT_SEG) + ICOUNT - 1
              IPT_EXP_CYC = IPT_EXP - IST_EXP_M1 + INC_CYC*(J-1)
              IPT_CMP = INDX_CMP(IPT_SEG) + ICOUNT - 1
              FLD_CMP(IPT_CMP - IST_CMP_M1) = FLD_EXP(IPT_EXP_CYC)
            END DO  ! index ICOUNT
!
!L 2.7 Case of current field:
!
            IF (I_FLD_TYP  ==  1) THEN
!
!L  Last value in segment is only retained if grid is cyclic and
!L  the first point on the row is ICYC-1 points before it.
              IF(LL_CYC_M)THEN
                IF(IPT_EXP-INDX_EXP(IPT_SEG+1-ISEG)  /=  ICYC-1)THEN
                  FLD_CMP(IPT_CMP - IST_CMP_M1) = REAL_MDI
                END IF
              ELSE
                FLD_CMP(IPT_CMP - IST_CMP_M1) = REAL_MDI
              END IF
            END IF
!
         END DO  ! index ISEG
!
      END DO  ! index J
!
!L  End loop over rows
!
      RETURN
      END SUBROUTINE OA_PACK
!
